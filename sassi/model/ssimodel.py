"""The in-memory database of one numbered model (ARCHITECTURE section 8, spec 08 section 1.7).

:class:`SSIModel` is what the command interpreter edits and what AFWRITE reads to build the
module decks.  It holds raw user data only; physical interpretation is deferred to
CHECK/AFWRITE (D-MDL-10, D-AFW-06).  Helper methods give the derived views every downstream
package needs (global node coordinates, element node tables by type, nodal masses in mass units,
the frequency step, a content hash).

Bookkeeping rules implemented here:

* nodes are stored in the coordinate system active at definition (D-MDL-04);
* "last defined" for nodes and loads is the definition order, kept in per-kind histories where a
  redefinition moves the id to the end and a deletion removes it (D-MDL-05; spec 08 section 1.5);
  for property tables it is the highest index;
* NDEL removes the loads, masses and mass-unit flags of the deleted nodes; element references are
  left dangling for CHECK Error 41 (D-MDL-11).
"""
from __future__ import annotations

import copy as _copy
import hashlib
import json
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np

from .entities import (BeamSection, CoordSys, ElementRequest, Element, Group, MatrixProp, Material,
                       NodalLoad, NodalRequest, Node, RelDispRequest, SoilLayer, SpringProp)
from .materials import G_DEFAULT, to_mass
from .options import MoptRecord, OptionStore

SCHEMA_VERSION = "1.0"   # major.minor of the SAVE/RESUME representation (D-MDL-02)

# Argument positions (1-based) of the shared variables in their owning records (spec 07 section 3)
HOUSE_GRAVITY, HOUSE_GELEV = 1, 2
SITE_FSTEP, SITE_DELT, SITE_NFT, SITE_FREQ = 3, 12, 13, 14


class History:
    """Definition order of ids ("last defined", D-MDL-05; spec 08 section 1.5).

    A (re)definition moves the id to the end, a deletion removes it; both are O(1) (an
    insertion-ordered dict), so building a model with many nodes stays linear.  Behaves like a
    read-only list for iteration, ``len``, indexing (``h[-1]`` = last defined) and ``==``.
    """

    __slots__ = ("_d",)

    def __init__(self, ids: Iterable[int] = ()):
        self._d: Dict[int, None] = dict.fromkeys(int(i) for i in ids)

    def touch(self, i: int) -> None:
        self._d.pop(i, None)
        self._d[i] = None

    def discard(self, ids: Iterable[int]) -> None:
        for i in ids:
            self._d.pop(i, None)

    def __len__(self) -> int:
        return len(self._d)

    def __iter__(self) -> Iterator[int]:
        return iter(self._d)

    def __contains__(self, i) -> bool:
        return i in self._d

    def __getitem__(self, idx):
        if isinstance(idx, int) and idx < 0 and -idx <= len(self._d):
            for k, v in enumerate(reversed(self._d), start=1):   # reversed(dict): Python >= 3.8
                if k == -idx:
                    return v
        return list(self._d)[idx]

    def __eq__(self, other) -> bool:
        if isinstance(other, History):
            return list(self._d) == list(other._d)
        if isinstance(other, (list, tuple)):
            return list(self._d) == list(other)
        return NotImplemented

    def __repr__(self) -> str:
        return f"History({list(self._d)})"


def touch_history(hist: History, i: int) -> None:
    """Move ``i`` to the end of a definition history (D-MDL-05: "last defined" = definition order)."""
    hist.touch(i)


def _hist_remove(hist: History, ids: Iterable[int]) -> None:
    hist.discard(ids)


@dataclass
class SSIModel:
    """One model of the interpreter (models are numbered 0, 1, ...; ACTM selects the active one)."""

    # identity
    name: str = ""                 # MDL / MDLNAME
    path: str = ""                 # MDL model directory
    title: str = ""                # TIT
    # geometry
    csys: Dict[int, CoordSys] = field(default_factory=dict)
    csys_active: int = 0
    nodes: Dict[int, Node] = field(default_factory=dict)
    node_history: History = field(default_factory=History)
    # elements
    groups: Dict[int, Group] = field(default_factory=dict)
    group_active: Optional[int] = None
    mact: int = 1                  # D-MDL-17: initial 1, global
    ract: int = 1
    # property tables
    materials: Dict[int, Material] = field(default_factory=dict)       # M
    layers: Dict[int, SoilLayer] = field(default_factory=dict)         # L
    sections: Dict[int, BeamSection] = field(default_factory=dict)     # R
    springs: Dict[int, SpringProp] = field(default_factory=dict)       # SC
    matrices: Dict[int, MatrixProp] = field(default_factory=dict)      # MXR/MXI/MXM
    # loads and masses
    forces: Dict[int, NodalLoad] = field(default_factory=dict)         # F
    moments: Dict[int, NodalLoad] = field(default_factory=dict)        # MM
    tmass: Dict[int, List[float]] = field(default_factory=dict)        # MT
    rmass: Dict[int, List[float]] = field(default_factory=dict)        # MR
    mass_units: Dict[int, int] = field(default_factory=dict)           # MUNITS (absent = 1, weight)
    force_history: History = field(default_factory=History)
    moment_history: History = field(default_factory=History)
    tmass_history: History = field(default_factory=History)
    rmass_history: History = field(default_factory=History)
    # frequency sets and lists (list append/delete class)
    freq_sets: Dict[int, List[int]] = field(default_factory=dict)      # FREQ
    damp: List[float] = field(default_factory=list)                    # DAMP
    topl: List[int] = field(default_factory=list)                      # TOPL
    amp: Dict[int, List[float]] = field(default_factory=dict)          # AMP per motion
    # output requests (request append class)
    nout: List[NodalRequest] = field(default_factory=list)             # NOUT
    eout: List[ElementRequest] = field(default_factory=list)           # EOUT
    rdnd: List[RelDispRequest] = field(default_factory=list)           # RDND
    # module option records (record / indexed / string setters)
    options: OptionStore = field(default_factory=OptionStore)
    # JSON-able data of other packages (saved by SAVE; written by WRITE through writer hooks)
    extensions: Dict[str, Any] = field(default_factory=dict)
    # GUI state saved by SAVE but not by WRITE (hide sets, D-UI-11)
    ui_state: Dict[str, Any] = field(default_factory=dict)

    MUNITS_DEFAULT = 1   # D-MDL-08

    # ================================================================== shared variables
    @property
    def gravity(self) -> float:
        """SSI acceleration of gravity: HOUSE ``<gravity>`` (also set by GRAVITY); default 32.2."""
        rec = self.options.record("HOUSE")
        v = rec.number(HOUSE_GRAVITY) if rec is not None else None
        return G_DEFAULT if v is None else v

    @property
    def ground_elevation(self) -> float:
        """HOUSE ``<gelev>`` (also set by GROUNDELEV); default 0."""
        rec = self.options.record("HOUSE")
        v = rec.number(HOUSE_GELEV) if rec is not None else None
        return 0.0 if v is None else v

    @property
    def mopt(self) -> MoptRecord:
        """MOPT record (defaults when MOPT was never given)."""
        rec = self.options.record("MOPT")
        if rec is None:
            return MoptRecord()
        return rec  # type: ignore[return-value]

    def frequency_step(self) -> float:
        """Resolved frequency step: SITE ``<fstep>`` if > 0, else ``1/(delt*nft)`` (spec 07 section 9.2.15)."""
        from ..conventions import frequency_step
        from ..io.decks import SCHEMAS
        s = SCHEMAS["SITE"]
        rec = self.options.record("SITE")
        fstep = rec.number(SITE_FSTEP) if rec is not None else None
        delt = rec.number(SITE_DELT) if rec is not None else None
        nft = rec.integer(SITE_NFT) if rec is not None else None
        fstep = s.param("fstep").default if fstep is None else fstep
        delt = s.param("delt").default if delt is None else delt
        nft = s.param("nft").default if nft is None else nft
        return frequency_step(delt, nft, fstep)

    def frequencies(self, set_no: int) -> Tuple[List[int], List[float]]:
        """Sorted frequency numbers of a FREQ set and their values in Hz (f = n * df)."""
        nums = sorted(self.freq_sets.get(set_no, []))
        df = self.frequency_step()
        return nums, [n * df for n in nums]

    # ================================================================== nodes and coordinates
    def define_node(self, nid: int, xyz: Sequence[float], csys: int) -> bool:
        """Create or redefine node ``nid`` with coordinates ``xyz`` in system ``csys``.

        A redefinition replaces the coordinates and the system but keeps the fixities and INT
        flags (spec 08 section 3.10).  Returns True when the node existed before.
        """
        n = self.nodes.get(nid)
        existed = n is not None
        if n is None:
            n = Node(nid)
            self.nodes[nid] = n
        n.x, n.y, n.z = (float(xyz[0]), float(xyz[1]), float(xyz[2]))
        n.csys = int(csys)
        touch_history(self.node_history, nid)
        return existed

    def delete_nodes(self, ids: Iterable[int]) -> int:
        """Delete nodes with the D-MDL-11 cascade; returns the number of dangling element references."""
        ids = [i for i in ids if i in self.nodes]
        s = set(ids)
        for i in ids:
            del self.nodes[i]
            self.forces.pop(i, None)
            self.moments.pop(i, None)
            self.tmass.pop(i, None)
            self.rmass.pop(i, None)
            self.mass_units.pop(i, None)
        for h in (self.node_history, self.force_history, self.moment_history, self.tmass_history,
                  self.rmass_history):
            _hist_remove(h, s)
        dangling = 0
        if s:
            for g in self.groups.values():
                for e in g.elements.values():
                    dangling += sum(1 for n in e.nodes if n in s)
        return dangling

    def node_ids(self) -> List[int]:
        return sorted(self.nodes)

    def system(self, s: int) -> Optional[CoordSys]:
        return self.csys.get(s)

    def to_global(self, xyz: Sequence[float], s: int) -> np.ndarray:
        """Global coordinates of a point given in system ``s`` (0 = global)."""
        p = np.asarray(xyz, float)
        if s == 0:
            return p.copy()
        cs = self.csys.get(s)
        if cs is None:
            raise KeyError(f"coordinate system {s} is not defined")
        return cs.to_global(p)

    def to_system(self, XYZ: Sequence[float], s: int) -> np.ndarray:
        """Coordinates in system ``s`` of a point given in global coordinates."""
        P = np.asarray(XYZ, float)
        if s == 0:
            return P.copy()
        cs = self.csys.get(s)
        if cs is None:
            raise KeyError(f"coordinate system {s} is not defined")
        return cs.to_local(P)

    def node_global(self, nid: int) -> np.ndarray:
        """Global coordinates of node ``nid`` (local systems resolved, spec 08 section 2.2)."""
        n = self.nodes[nid]
        return self.to_global(n.xyz, n.csys)

    def node_in_system(self, nid: int, s: int) -> np.ndarray:
        """Coordinates of node ``nid`` in system ``s``; exact (no round trip) when stored in ``s``."""
        n = self.nodes[nid]
        if n.csys == s:
            return np.array(n.xyz, float)
        return self.to_system(self.to_global(n.xyz, n.csys), s)

    def global_coordinates(self, ids: Optional[Sequence[int]] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Vectorised global coordinates: returns ``(ids (n,), xyz (n, 3))`` sorted by id."""
        if ids is None:
            ids = self.node_ids()
        ids_a = np.asarray(list(ids), dtype=np.int64)
        xyz = np.zeros((len(ids_a), 3))
        cs = np.zeros(len(ids_a), dtype=np.int64)
        for k, i in enumerate(ids_a):
            n = self.nodes[int(i)]
            xyz[k] = n.xyz
            cs[k] = n.csys
        for s in np.unique(cs):
            if s == 0:
                continue
            m = cs == s
            sysobj = self.csys.get(int(s))
            if sysobj is None:
                raise KeyError(f"coordinate system {int(s)} is not defined")
            xyz[m] = sysobj.to_global(xyz[m])
        return ids_a, xyz

    def nodes_in_system(self, s: int) -> List[int]:
        return sorted(i for i, n in self.nodes.items() if n.csys == s)

    # ================================================================== groups and elements
    @property
    def active_group(self) -> Optional[Group]:
        return self.groups.get(self.group_active) if self.group_active is not None else None

    def iter_elements(self) -> Iterator[Tuple[Group, Element]]:
        """All (group, element) pairs, groups and elements in ascending number."""
        for gid in sorted(self.groups):
            g = self.groups[gid]
            for e in g.sorted_elements():
                yield g, e

    def n_elements(self) -> int:
        return sum(len(g.elements) for g in self.groups.values())

    def elements_by_type(self) -> Dict[int, List[Tuple[int, Element]]]:
        """Group type code -> list of (group id, element)."""
        out: Dict[int, List[Tuple[int, Element]]] = {}
        for g, e in self.iter_elements():
            out.setdefault(g.type, []).append((g.id, e))
        return out

    def element_node_table(self, type_code: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Element node lists of one element type as arrays.

        Returns ``(group_ids (n,), element_ids (n,), nodes (n, m))`` with ``m`` the largest node
        count of that type in the model; missing slots are 0.
        """
        rows = self.elements_by_type().get(type_code, [])
        m = max((len(e.nodes) for _, e in rows), default=0)
        gids = np.array([g for g, _ in rows], dtype=np.int64)
        eids = np.array([e.id for _, e in rows], dtype=np.int64)
        nodes = np.zeros((len(rows), m), dtype=np.int64)
        for k, (_, e) in enumerate(rows):
            nodes[k, :len(e.nodes)] = e.nodes
        return gids, eids, nodes

    def used_nodes(self) -> List[int]:
        """Nodes referenced by any element (beam K nodes included)."""
        s = set()
        for g in self.groups.values():
            for e in g.elements.values():
                s.update(n for n in e.nodes if n)
        return sorted(s)

    # ================================================================== masses and loads
    def mass_unit(self, nid: int) -> int:
        """MUNITS flag of a node: 0 mass, 1 weight (default 1, D-MDL-08)."""
        return self.mass_units.get(nid, self.MUNITS_DEFAULT)

    def nodal_masses(self, gravity: Optional[float] = None) -> Dict[int, np.ndarray]:
        """Nodal masses in **mass units**: node -> [mx, my, mz, mxx, myy, mzz].

        Values entered in weight units (MUNITS 1) are divided by g (HOUSE gravity unless
        ``gravity`` is given); spec 08 section 7.19, UT-15.
        """
        g = self.gravity if gravity is None else gravity
        out: Dict[int, np.ndarray] = {}
        for nid in sorted(set(self.tmass) | set(self.rmass)):
            m = np.zeros(6)
            if nid in self.tmass:
                m[:3] = self.tmass[nid]
            if nid in self.rmass:
                m[3:] = self.rmass[nid]
            out[nid] = to_mass(m, self.mass_unit(nid), g)
        return out

    def general_matrices(self, prop: int, gravity: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
        """GENERAL element property ``prop``: complex stiffness and mass in mass units (MOPT ``<matrix>``)."""
        p = self.matrices[prop]
        g = self.gravity if gravity is None else gravity
        return p.stiffness(), p.mass(self.mopt.get("matrix"), g)

    # ================================================================== summaries
    def counts(self) -> Dict[str, int]:
        return dict(nodes=len(self.nodes), csys=len(self.csys), groups=len(self.groups),
                    elements=self.n_elements(), materials=len(self.materials), layers=len(self.layers),
                    sections=len(self.sections), springs=len(self.springs), matrices=len(self.matrices),
                    forces=len(self.forces), moments=len(self.moments), tmass=len(self.tmass),
                    rmass=len(self.rmass), interaction=sum(1 for n in self.nodes.values() if 0 in n.flags),
                    freq_sets=len(self.freq_sets))

    def is_empty(self) -> bool:
        return not (self.nodes or self.groups or self.materials or self.layers or self.options.names())

    # ================================================================== serialisation
    def to_json(self) -> Dict[str, Any]:
        """Complete JSON-able state (SAVE, D-MDL-02)."""
        return {
            "schema": SCHEMA_VERSION,
            "name": self.name, "path": self.path, "title": self.title,
            "csys": [self.csys[k].to_json() for k in sorted(self.csys)],
            "csys_active": self.csys_active,
            "nodes": [self.nodes[k].to_json() for k in sorted(self.nodes)],
            "node_history": list(self.node_history),
            "groups": [self.groups[k].to_json() for k in sorted(self.groups)],
            "group_active": self.group_active,
            "mact": self.mact, "ract": self.ract,
            "materials": [self.materials[k].to_json() for k in sorted(self.materials)],
            "layers": [self.layers[k].to_json() for k in sorted(self.layers)],
            "sections": [self.sections[k].to_json() for k in sorted(self.sections)],
            "springs": [self.springs[k].to_json() for k in sorted(self.springs)],
            "matrices": [self.matrices[k].to_json() for k in sorted(self.matrices)],
            "forces": [self.forces[k].to_json() for k in sorted(self.forces)],
            "moments": [self.moments[k].to_json() for k in sorted(self.moments)],
            "tmass": [[k, list(self.tmass[k])] for k in sorted(self.tmass)],
            "rmass": [[k, list(self.rmass[k])] for k in sorted(self.rmass)],
            "mass_units": [[k, self.mass_units[k]] for k in sorted(self.mass_units)],
            "force_history": list(self.force_history), "moment_history": list(self.moment_history),
            "tmass_history": list(self.tmass_history), "rmass_history": list(self.rmass_history),
            "freq_sets": [[k, list(self.freq_sets[k])] for k in sorted(self.freq_sets)],
            "damp": list(self.damp), "topl": list(self.topl),
            "amp": [[k, list(self.amp[k])] for k in sorted(self.amp)],
            "nout": [r.to_json() for r in self.nout],
            "eout": [r.to_json() for r in self.eout],
            "rdnd": [r.to_json() for r in self.rdnd],
            "options": self.options.to_json(),
            "extensions": _copy.deepcopy(self.extensions),
            "ui_state": _copy.deepcopy(self.ui_state),
        }

    @classmethod
    def from_json(cls, d: Dict[str, Any]) -> "SSIModel":
        major = str(d.get("schema", "1.0")).split(".")[0]
        if major != SCHEMA_VERSION.split(".")[0]:
            raise ValueError(f"model schema version {d.get('schema')} is not supported "
                             f"(this build reads {SCHEMA_VERSION.split('.')[0]}.x)")
        m = cls(name=d.get("name", ""), path=d.get("path", ""), title=d.get("title", ""))
        for c in d.get("csys", []):
            cs = CoordSys.from_json(c)
            m.csys[cs.id] = cs
        m.csys_active = int(d.get("csys_active", 0))
        for n in d.get("nodes", []):
            node = Node.from_json(n)
            m.nodes[node.id] = node
        m.node_history = History(d.get("node_history", []))
        for g in d.get("groups", []):
            grp = Group.from_json(g)
            m.groups[grp.id] = grp
        ga = d.get("group_active")
        m.group_active = None if ga is None else int(ga)
        m.mact, m.ract = int(d.get("mact", 1)), int(d.get("ract", 1))
        for key, klass, tab in (("materials", Material, m.materials), ("layers", SoilLayer, m.layers),
                                ("sections", BeamSection, m.sections), ("springs", SpringProp, m.springs)):
            for row in d.get(key, []):
                obj = klass.from_json(row)
                tab[obj.id] = obj
        for row in d.get("matrices", []):
            p = MatrixProp.from_json(row)
            m.matrices[p.id] = p
        for row in d.get("forces", []):
            ld = NodalLoad.from_json(row)
            m.forces[ld.node] = ld
        for row in d.get("moments", []):
            ld = NodalLoad.from_json(row)
            m.moments[ld.node] = ld
        m.tmass = {int(k): [float(x) for x in v] for k, v in d.get("tmass", [])}
        m.rmass = {int(k): [float(x) for x in v] for k, v in d.get("rmass", [])}
        m.mass_units = {int(k): int(v) for k, v in d.get("mass_units", [])}
        for key in ("force_history", "moment_history", "tmass_history", "rmass_history"):
            setattr(m, key, History(d.get(key, [])))
        m.freq_sets = {int(k): [int(x) for x in v] for k, v in d.get("freq_sets", [])}
        m.damp = [float(v) for v in d.get("damp", [])]
        m.topl = [int(v) for v in d.get("topl", [])]
        m.amp = {int(k): [float(x) for x in v] for k, v in d.get("amp", [])}
        m.nout = [NodalRequest.from_json(r) for r in d.get("nout", [])]
        m.eout = [ElementRequest.from_json(r) for r in d.get("eout", [])]
        m.rdnd = [RelDispRequest.from_json(r) for r in d.get("rdnd", [])]
        m.options = OptionStore.from_json(d.get("options", {}))
        m.extensions = _copy.deepcopy(d.get("extensions", {}))
        m.ui_state = _copy.deepcopy(d.get("ui_state", {}))
        return m

    def copy(self) -> "SSIModel":
        """Deep copy (CPMODEL)."""
        return SSIModel.from_json(json.loads(json.dumps(self.to_json())))

    #: keys of :meth:`to_json` that are not part of the *model content*
    IDENTITY_KEYS = ("name", "path")
    HISTORY_KEYS = ("node_history", "force_history", "moment_history", "tmass_history", "rmass_history")
    UI_KEYS = ("ui_state",)

    def canonical(self, identity: bool = False, history: bool = False, ui: bool = False,
                  title: bool = True) -> Dict[str, Any]:
        """Canonical state for comparisons (UT-03 deep equality) and hashing.

        By default the identity (name, path), the definition histories and the GUI state are left
        out: a WRITE -> INP round trip reproduces the model *content*; MDL is only written on
        request and the replay order of N commands is canonical, not historical.
        """
        d = self.to_json()
        d.pop("schema", None)
        drop = []
        if not identity:
            drop += list(self.IDENTITY_KEYS)
        if not history:
            drop += list(self.HISTORY_KEYS)
        if not ui:
            drop += list(self.UI_KEYS)
        if not title:
            drop.append("title")
        for k in drop:
            d.pop(k, None)
        return d

    def same_state(self, other: "SSIModel", **kw) -> bool:
        return json.dumps(self.canonical(**kw), sort_keys=True) == json.dumps(other.canonical(**kw), sort_keys=True)

    def model_hash(self) -> str:
        """SHA-256 of the model content (identity, title, histories and GUI state excluded).

        Modules record it in the meta data of the binary files they write so that a restart can
        detect a model that changed since (D-GEN-03).
        """
        text = json.dumps(self.canonical(title=False), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(text.encode("utf-8")).hexdigest()
