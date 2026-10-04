"""Model generation and combination algorithms (requirements section 3.4.I; spec 09 sections 2 and 4).

Pure functions on :class:`sassi.model.SSIModel`.  The command handlers of
:mod:`sassi.prep.commands.generation` (ETYPEGEN, INTGEN, RADIUS, GLB2LOC, GROUPMAT, RMVUNUSED, NCOM,
GCOM, WELD, MERGE, MERGESOIL, MERGEGROUP, ROTATE, TRANSLATE, EXCAV, SOILMESH) call them and print
the collected :class:`Notes`; Python scripts and the GUI can call them directly.

Background for the structural engineer
--------------------------------------
In the flexible-volume (FV) method of Lysmer et al. (1981) the SSI model is the structure plus the
*excavated soil* -- an FE model of the soil volume that the embedded structure displaces, with the
free-field layer properties.  The far-field soil acts on the model only through the *interaction
nodes* (R1 section 4.4).  The substructuring variants differ only in which excavation nodes are
interaction nodes (R1 section 4.5; spec 01 section 6 rules 6-8):

* **FV** (direct method): every node of the excavated soil;
* **FI-EVBN** / modified subtraction (MSM): the nodes on the whole boundary of the excavated volume,
  ground-surface face included;
* **FI-FSIN** / subtraction (SM): only the foundation-soil interface (lateral and bottom faces);
* **FFV** (fast flexible volume): EVBN plus internal horizontal node levels at a regular skip;
* **surface** foundations: the foundation nodes at grade.

INTGEN derives these sets from the excavation elements; EXCAV builds the excavation mesh from the
basement; MERGESOIL joins a structure model and an excavation model (shared nodes or stiff
springs at the interface); WELD / RMVUNUSED / NCOM clean up node numbering afterwards.

Geometric tolerance (D-GEN-06): ``tol = max(1e-6 L_ref, 1e-9)`` of :class:`sassi.prep.check.ModelView`
(``EDUOPT,GEOMTOL`` overrides).  Every function validates its input before it changes anything
and raises :class:`GenerationError` (the model is then unchanged).
"""
from __future__ import annotations

import bisect
import math
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

from ..model import SSIModel
from ..model.entities import (BeamSection, CoordSys, Element, ElementRequest, Group, MatrixProp, Material,
                              NodalLoad, NodalRequest, Node, RelDispRequest, SoilLayer, SpringProp)
from ..model.geometry import euler_from_matrix, loc_matrix, local_from_points, rot_x, rot_y, rot_z
from ..model.materials import elastic_constants
from ..model.options import make_record
from ..model.ssimodel import History
from ..model.values import fmt_num
from .check import (BEAMS, GENERAL, PLANE, SHELL, SOLID, SOLID_FACES, SPRING, TSHELL, ElemRef,
                    ModelView)
from .options import eduopt_float

#: element types that carry an ETYPE (spec 09 section 1.3)
ETYPE_TYPES = (SOLID, PLANE, SHELL, TSHELL)
#: element types that use the M table (SOLID/PLANE only when structural)
MATERIAL_TYPES = (SOLID, PLANE, SHELL, TSHELL, BEAMS)
#: maximum node number of the IKTR9-compatible solver build (spec 09 section 1.6)
NODE_LIMIT = 99_999

#: INTGEN method names (spec 09 section 2.19)
INTGEN_METHODS = {0: "clear", 1: "FV (flexible volume)", 2: "FI-EVBN (modified subtraction, MSM)",
                  3: "FI-FSIN (subtraction, SM)", 4: "surface foundation", 5: "FFV (fast flexible volume)"}


class GenerationError(ValueError):
    """Invalid input of a generation command; the model is left unchanged."""


class Notes:
    """Warnings and information collected by a generation function (the handler prints them)."""

    def __init__(self) -> None:
        self.warnings: List[str] = []
        self.infos: List[str] = []

    def warn(self, text: str) -> None:
        self.warnings.append(text)

    def info(self, text: str) -> None:
        self.infos.append(text)


def _notes(notes: Optional[Notes]) -> Notes:
    return notes if notes is not None else Notes()


def _short(ids: Sequence[int], k: int = 10) -> str:
    ids = list(ids)
    s = ", ".join(str(i) for i in ids[:k])
    return s + (f", ... ({len(ids)} in total)" if len(ids) > k else "")


def check_node_limit(m: SSIModel, notes: Notes) -> None:
    """Warn when node numbers exceed the 99,999 limit of the IKTR9 build (spec 09 section 1.6)."""
    if m.nodes and max(m.nodes) > NODE_LIMIT:
        notes.warn(f"node numbers up to {max(m.nodes)} exceed {NODE_LIMIT:,}: run RMVUNUSED and NCOM")


# ======================================================================================
# Small geometry helpers
# ======================================================================================
def cluster_values(values: Sequence[float], thr: float) -> Tuple[np.ndarray, np.ndarray]:
    """Cluster values into levels (spec 09 section 4.1 step 2, D-MDL-19).

    The values are sorted; a new level starts whenever ``value - first value of the current level
    > thr``.  Returns ``(level means ascending, label of every input value)``.
    """
    v = np.asarray(values, float)
    if v.size == 0:
        return np.zeros(0), np.zeros(0, dtype=np.int64)
    order = np.argsort(v, kind="stable")
    labels = np.empty(v.size, dtype=np.int64)
    k, start = 0, v[order[0]]
    for idx in order:
        if v[idx] - start > thr:
            k += 1
            start = v[idx]
        labels[idx] = k
    means = np.bincount(labels, weights=v) / np.bincount(labels)
    return means, labels


def signed_area(xy: np.ndarray) -> float:
    """Signed (shoelace) area of a polygon in the XY plane; > 0 counter-clockwise seen from above."""
    p = np.asarray(xy, float)
    if len(p) < 3:
        return 0.0
    x, y = p[:, 0], p[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def convex_hull(xy: np.ndarray) -> List[int]:
    """Indices of the convex-hull vertices of 2D points, counter-clockwise (Andrew's monotone chain;
    collinear points on the hull are excluded)."""
    p = np.asarray(xy, float)
    n = len(p)
    if n < 3:
        return list(range(n))
    order = sorted(range(n), key=lambda i: (p[i, 0], p[i, 1]))

    def cross(o, a, b):
        return (p[a, 0] - p[o, 0]) * (p[b, 1] - p[o, 1]) - (p[a, 1] - p[o, 1]) * (p[b, 0] - p[o, 0])

    lower: List[int] = []
    for i in order:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], i) <= 0:
            lower.pop()
        lower.append(i)
    upper: List[int] = []
    for i in reversed(order):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], i) <= 0:
            upper.pop()
        upper.append(i)
    return lower[:-1] + upper[:-1]


def hull_area(xy: np.ndarray) -> float:
    """Area of the convex hull of 2D points (0 for fewer than 3 non-collinear points)."""
    p = np.asarray(xy, float)
    if len(p) < 3:
        return 0.0
    h = convex_hull(p)
    return abs(signed_area(p[h])) if len(h) >= 3 else 0.0


def unique_points(xy: np.ndarray, tol: float) -> np.ndarray:
    """Label coincident points (distance <= tol, transitively); labels are 0..k-1 in order of
    first appearance."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    p = np.asarray(xy, float)
    n = len(p)
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    pairs = cKDTree(p).query_pairs(r=tol, output_type="ndarray")
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)) if len(pairs) else \
        coo_matrix((n, n))
    _, lab = connected_components(g, directed=False)
    first: Dict[int, int] = {}
    out = np.empty(n, dtype=np.int64)
    for i, l in enumerate(lab):
        out[i] = first.setdefault(int(l), len(first))
    return out


def set_geomtol(m: SSIModel, tol: float) -> None:
    """Store ``EDUOPT,GEOMTOL,<tol>`` (the D-GEN-06 override of the geometric tolerance) in a model;
    WRITE writes it, so it survives WRITE -> INP."""
    m.options.set_entry("EDUOPT", "GEOMTOL", make_record("EDUOPT", ["GEOMTOL", fmt_num(float(tol))]))


def _ordered_unique(seq: Iterable[int]) -> List[int]:
    """Polygon node list without repeated nodes (consecutive or not), order kept."""
    return list(dict.fromkeys(n for n in seq if n))


def below_grade(v: ModelView, nodes: Sequence[int]) -> bool:
    """D-HOU-01: an element is below grade when every node has ``z <= gelev + tol`` and the centroid
    of its distinct nodes is strictly below ``gelev`` (``gelev - tol``)."""
    ns = [n for n in nodes if n]
    if not ns or not v.defined(ns):
        return False
    z = np.array([v.P[n][2] for n in ns])
    uniq = np.array([v.P[n][2] for n in dict.fromkeys(ns)])
    return bool(np.all(z <= v.gelev + v.tol) and float(uniq.mean()) < v.gelev - v.tol)


def _layer_index(v: ModelView, depth: float) -> Tuple[Optional[int], bool]:
    """TOPL position (0-based) of the layer containing ``depth`` (None without TOPL) and whether the
    depth lies below the last TOPL layer."""
    depths = v.interfaces()
    if depths is None:
        return None, False
    k = int(np.searchsorted(depths, depth, side="right")) - 1
    nl = len(depths) - 1
    if k >= nl:
        return nl - 1, True
    return max(k, 0), False


def _off_interface_levels(v: ModelView, levels: Sequence[float]) -> List[float]:
    depths = v.interfaces()
    if depths is None:
        return []
    t = v.interface_tol(depths)
    return [z for z in levels if v.interface_of(z, depths)[1] > t]


# ======================================================================================
# ETYPEGEN (spec 09 section 2.6)
# ======================================================================================
def etypegen(m: SSIModel, etype: int, reset: bool = False, notes: Optional[Notes] = None) -> Dict[int, int]:
    """Set the ETYPE of every SOLID, PLANE, SHELL and TSHELL element (beams, springs and GENERAL
    elements are untouched; ETYPE means nothing for them).

    * ``etype`` 1 / 2: all elements structure / excavated soil (shells: buried shells, whose nodes
      INTGEN makes interaction nodes);
    * ``etype`` 0: the implicit classification is made explicit (spec 09 Q-G4): SOLID/PLANE below
      grade by D-HOU-01 -> 2, otherwise 1; SHELL/TSHELL -> 1 (implicit shells are structure);
    * ``reset``: every element back to the implicit ETYPE 0.

    Returns the number of elements set to each ETYPE value.
    """
    notes = _notes(notes)
    if etype not in (0, 1, 2):
        raise GenerationError("<type> must be 0 (by location), 1 (structure) or 2 (excavated soil)")
    v = ModelView(m) if (etype == 0 and not reset) else None
    counts = {0: 0, 1: 0, 2: 0}
    buried_reset = new_buried = 0
    for g, e in m.iter_elements():
        if g.type not in ETYPE_TYPES:
            continue
        if reset:
            new = 0
        elif etype in (1, 2):
            new = etype
        elif g.type in (SOLID, PLANE):
            new = 2 if below_grade(v, e.nodes) else 1
        else:
            new = 1
        if g.type in (SHELL, TSHELL):
            if e.etype == 2 and new != 2:
                buried_reset += 1
            if new == 2 and e.etype != 2:
                new_buried += 1
        e.etype = new
        counts[new] += 1
    if new_buried:
        notes.warn(f"{new_buried} SHELL/TSHELL elements are now buried shells (ETYPE 2): INTGEN makes all "
                   f"their nodes interaction nodes")
    if buried_reset:
        notes.warn(f"{buried_reset} buried shells (ETYPE 2) are now structural shells")
    return counts


# ======================================================================================
# INTGEN (spec 09 section 2.19; spec 01 section 6 rules 5-10; R1 section 4.5)
# ======================================================================================
@dataclass
class IntgenResult:
    """Interaction set of one INTGEN option (computed before it is applied)."""
    option: int
    nodes: Set[int] = field(default_factory=set)      # the method's interaction set
    added: List[int] = field(default_factory=list)    # nodes newly flagged by the command
    total: int = 0                                     # interaction nodes after the command
    excavation_elements: int = 0
    implicit: int = 0                                  # excavation elements with implicit ETYPE 0
    buried: int = 0                                    # nodes of buried shells
    levels: List[float] = field(default_factory=list)  # FFV: z-levels of the excavation (bottom up)
    selected_levels: List[int] = field(default_factory=list)   # FFV: selected internal level indices


def ffv_levels(v: ModelView, nodes: Iterable[int], skip: int) -> Tuple[np.ndarray, Dict[int, int], List[int]]:
    """Horizontal node levels of the excavation and the FFV selection (spec 09 Q-G9, decision in the
    module docstring of :mod:`sassi.prep.commands.generation`).

    Levels are the excavation-node elevations clustered with ``tol``, numbered j = 0 (bottom) ..
    L-1 (top).  Internal levels j = 1..L-2 are selected when ``j mod (skip + 1) == 0``: ``skip = 1``
    keeps every other internal level, ``skip = 0`` keeps all of them (identical to FV).
    """
    ns = sorted(nodes)
    z = np.array([v.P[n][2] for n in ns])
    means, lab = cluster_values(z, v.tol)
    level = {n: int(l) for n, l in zip(ns, lab)}
    L = len(means)
    sel = [j for j in range(1, L - 1) if j % (skip + 1) == 0]
    return means, level, sel


def interaction_set(m: SSIModel, option: int, skip: int = 1, v: Optional[ModelView] = None) -> IntgenResult:
    """The interaction-node set of INTGEN option 1..5 (nothing is changed).

    Excavation elements are the SOLID/PLANE elements with resolved ETYPE 2 (explicit 2, or implicit
    0 below grade by D-HOU-01).  Their boundary faces are the faces used by exactly one excavation
    element (:meth:`ModelView.excavation`); the ground-surface top faces are those with every node
    at ``gelev``.  Nodes of buried shells (SHELL/TSHELL ETYPE 2, spec 01 rule 10) are always added.
    """
    if option not in (1, 2, 3, 4, 5):
        raise GenerationError("<type> must be 0 (clear), 1 (FV), 2 (FI-EVBN), 3 (FI-FSIN), 4 (surface) or 5 (FFV)")
    if skip < 0:
        raise GenerationError("[level skip] must be >= 0")
    v = v or ModelView(m)
    res = IntgenResult(option)
    buried = v.buried_shell_nodes() & set(v.P)
    res.buried = len(buried)
    if option == 4:
        # surface foundation: foundation nodes at grade (no excavation needed; FV = FI = FFV).  Nodes
        # of SOLID/PLANE/SHELL/TSHELL elements and the I/J nodes of BEAMS (rigid beam grids); never
        # beam K nodes, springs or GENERAL elements (FIXROT ground nodes are fixed)
        s = {n for r in v.elems if r.type in ETYPE_TYPES + (BEAMS,) for n in r.dof_nodes
             if n in v.P and abs(v.P[n][2] - v.gelev) <= v.tol}
        if not s and not buried:
            raise GenerationError(f"no SOLID, PLANE, SHELL or BEAMS node lies at the ground elevation {v.gelev:g} "
                                  f"(set GROUNDELEV)")
        res.nodes = s | buried
        return res
    exc = v.excavation()
    X = exc["elements"]
    if not X:
        raise GenerationError("no excavated-soil elements (SOLID/PLANE with ETYPE 2): define the excavation "
                              "volume with ETYPE or ETYPEGEN")
    res.excavation_elements = len(X)
    res.implicit = sum(1 for r in X if r.elem.etype == 0)
    if option == 1:
        s = set(exc["nodes"])
    elif option == 2:
        s = set(exc["boundary"])
    elif option == 3:
        s = set(exc["fsin"])
    else:
        means, level, sel = ffv_levels(v, exc["nodes"], skip)
        res.levels = [float(z) for z in means]
        res.selected_levels = sel
        ss = set(sel)
        s = set(exc["boundary"]) | {n for n, j in level.items() if j in ss}
    res.nodes = s | buried
    return res


def set_interaction(m: SSIModel, nodes: Iterable[int]) -> List[int]:
    """Set the INT code 0 (interaction) of nodes, exactly as ``INT,n,n,1,1,0``; returns the nodes
    that were not interaction nodes before (D-MDL-07)."""
    added = []
    for n in sorted(nodes):
        nd = m.nodes.get(n)
        if nd is not None and 0 not in nd.flags:
            nd.flags.add(0)
            added.append(n)
    return added


def clear_interaction(m: SSIModel) -> int:
    """INTGEN,0: reset the interaction flag (code 0) of every node; returns the number reset."""
    k = 0
    for nd in m.nodes.values():
        if 0 in nd.flags:
            nd.flags.discard(0)
            k += 1
    return k


def interaction_diagnostics(m: SSIModel, v: ModelView, nodes: Iterable[int], notes: Notes) -> None:
    """Checks after generation (spec 09 section 2.19 "After generation", spec 01 rules 5 and 6)."""
    nodes = sorted(n for n in nodes if n in v.P)
    depths = v.interfaces()
    if depths is not None:
        t = v.interface_tol(depths)
        off = [n for n in nodes if v.interface_of(v.P[n][2], depths)[1] > t]
        if off:
            notes.warn(f"{len(off)} interaction nodes are not on a soil-layer interface (EDU-01): {_short(off)}")
    above = [n for n in nodes if v.P[n][2] > v.gelev + v.tol]
    if above:
        notes.warn(f"{len(above)} interaction nodes lie above the ground elevation {v.gelev:g}: {_short(above)}")
    fixed = [n for n in nodes if any(m.nodes[n].fix[:3])]
    if fixed:
        notes.warn(f"{len(fixed)} interaction nodes have fixed translations (FIXEDINT): {_short(fixed)}")
    allint = sorted(i for i, nd in m.nodes.items() if 0 in nd.flags and i in v.P)
    if len(allint) > 1:
        z = np.array([v.P[i][2] for i in allint])
        if np.any(np.diff(z) < -v.tol):
            notes.warn("interaction nodes are not numbered bottom-up (ascending numbers from the lowest level "
                       "to the ground surface); this matters for incoherent analysis")


def intgen(m: SSIModel, option: int, skip: int = 1, notes: Optional[Notes] = None) -> IntgenResult:
    """INTGEN,<type>,[level skip]: generate interaction nodes (spec 09 section 2.19).

    Option 0 clears every interaction flag; options 1-5 add the method's set to the existing
    interaction nodes (union, never removal).  Flags are set exactly as ``INT,n,n,1,1,0`` so
    AFWRITE writes them.
    """
    notes = _notes(notes)
    if option == 0:
        res = IntgenResult(0)
        clear_interaction(m)
        return res
    v = ModelView(m)
    res = interaction_set(m, option, skip, v)
    if res.implicit:
        notes.warn(f"{res.implicit} excavation elements have the implicit ETYPE 0 and were classified by the "
                   f"ground elevation {v.gelev:g} (D-HOU-01); the manual requires an explicit ETYPE (ETYPEGEN)")
    res.added = set_interaction(m, res.nodes)
    res.total = sum(1 for nd in m.nodes.values() if 0 in nd.flags)
    interaction_diagnostics(m, v, res.nodes, notes)
    return res


# ======================================================================================
# RADIUS (spec 09 section 2.25, D-PNT-05)
# ======================================================================================
def excavation_radii(m: SSIModel, scale: float = 0.9, v: Optional[ModelView] = None) -> List[Tuple[int, int, float]]:
    """POINT central-zone radius of every excavation element: ``r_e = Scale * sqrt(A_plan)``
    (D-PNT-05), ``A_plan`` = area of the convex hull of the element's plan (XY) projection; PLANE
    elements use their horizontal width ``h``, ``r_e = Scale * h``.

    For a uniform square mesh of spacing h, Scale = 0.9 gives the manual's rule R0 = 0.90 h
    (R1 section 3.4).  Returns ``[(group, element, r_e)]``.
    """
    v = v or ModelView(m)
    out = []
    for r in v.excavation()["elements"]:
        ns = [n for n in dict.fromkeys(r.elem.nodes) if n]
        pts = np.array([v.P[n] for n in ns])
        if r.type == PLANE:
            h = float(pts[:, 0].max() - pts[:, 0].min())
        else:
            h = math.sqrt(hull_area(pts[:, :2]))
        out.append((r.group, r.id, scale * h))
    return out


# ======================================================================================
# GLB2LOC (spec 09 section 2.15)
# ======================================================================================
def glb2loc(m: SSIModel, ids: Iterable[int], sysno: int) -> int:
    """Store nodes in Cartesian system ``sysno`` (0 = global): ``x_loc = R_s^T (x_g - O_s)``; nodes
    in another local system are converted to global first.  Returns the number of nodes converted."""
    if sysno != 0 and sysno not in m.csys:
        raise GenerationError(f"coordinate system {sysno} is not defined")
    ids = [i for i in ids if i in m.nodes]
    k = 0
    for i in ids:
        n = m.nodes[i]
        if n.csys == sysno:
            continue
        x = m.to_system(m.node_global(i), sysno)
        n.x, n.y, n.z = (float(x[0]), float(x[1]), float(x[2]))
        n.csys = sysno
        k += 1
    return k


# ======================================================================================
# GROUPMAT (spec 09 section 2.18; P2 Option NON preparation)
# ======================================================================================
def groupmat(m: SSIModel, notes: Optional[Notes] = None) -> Dict[int, int]:
    """One material per group: material k (k = 1, 2, ... in group order) copies the material of the
    group's first element and is assigned to all its elements; the original M table is replaced.

    Groups without materials (SPRING, GENERAL) and excavated soil elements (their MSET points to the
    L table) are skipped (spec 09 Q-G8).  Returns ``{group: new material}``.
    """
    notes = _notes(notes)
    v = ModelView(m)
    excavated = {(r.group, r.id) for r in v.elems if r.excavated}
    plan = []
    for gid in sorted(m.groups):
        g = m.groups[gid]
        if g.type not in MATERIAL_TYPES:
            continue
        els = [e for e in g.sorted_elements() if (gid, e.id) not in excavated]
        if not els:
            continue
        m0 = els[0].mat
        if m0 not in m.materials:
            raise GenerationError(f"material {m0} of group {gid}, element {els[0].id} is not defined")
        plan.append((gid, m0, els))
    mats: Dict[int, Material] = {}
    out: Dict[int, int] = {}
    for k, (gid, m0, els) in enumerate(plan, start=1):
        s = m.materials[m0]
        mats[k] = Material(k, s.val1, s.val2, s.weight, s.pdamp, s.sdamp, s.mtype)
        other = sum(1 for e in els if e.mat != m0)
        if other:
            notes.warn(f"group {gid}: {other} elements had another material than {m0}; they now use material {k}")
        for e in els:
            e.mat = k
        out[gid] = k
    m.materials = mats
    return out


# ======================================================================================
# Node references: RMVUNUSED, NCOM (spec 09 sections 2.22, 2.26)
# ======================================================================================
def _symm_nodes(m: SSIModel) -> Set[int]:
    s: Set[int] = set()
    for _, rec in m.options.entries("SYMM"):
        for k in (3, 4, 5):
            v = rec.integer(k)
            if v:
                s.add(v)
    return s


def rmvunused(m: SSIModel, notes: Optional[Notes] = None) -> List[int]:
    """Delete the nodes that no element references (beam / GENERAL K nodes count as used) and that
    are not interaction nodes; their masses, loads and fixities go with them (NDEL cascade,
    D-MDL-11), and they are removed from the NOUT / RDND output requests.  Node numbers and
    connectivity are not changed (use NCOM).  Nodes referenced by SYMM are kept."""
    notes = _notes(notes)
    keep = set(m.used_nodes()) | {i for i, n in m.nodes.items() if 0 in n.flags} | _symm_nodes(m)
    remove = sorted(set(m.nodes) - keep)
    if not remove:
        return []
    rs = set(remove)
    loaded = [n for n in remove if n in m.tmass or n in m.rmass or n in m.forces or n in m.moments]
    if loaded:
        notes.warn(f"{len(loaded)} removed nodes carried masses or loads (dropped): {_short(loaded)}")
    m.delete_nodes(remove)
    k = 0
    for r in m.nout:
        n0 = len(r.nodes)
        r.nodes = [n for n in r.nodes if n not in rs]
        k += n0 - len(r.nodes)
    m.nout = [r for r in m.nout if r.nodes]
    n0 = len(m.rdnd)
    m.rdnd = [r for r in m.rdnd if r.node not in rs]
    k += n0 - len(m.rdnd)
    if k:
        notes.warn(f"{k} NOUT/RDND output requests of removed nodes dropped")
    return remove


def renumber_nodes(m: SSIModel, mp: Dict[int, int]) -> None:
    """Apply an old -> new node-number map to every node reference of the model: the node table and
    its definition history, element connectivity (K nodes included), LOCAL defining nodes, fixities
    and INT flags (they belong to the node), F/MM loads, MT/MR masses, MUNITS flags and their
    histories, NOUT/RDND requests and the SYMM nodes.  Numbers not in ``mp`` are kept."""
    def f(n: int) -> int:
        return mp.get(n, n)

    nodes: Dict[int, Node] = {}
    for old in sorted(m.nodes):
        nd = m.nodes[old]
        nd.id = f(old)
        nodes[nd.id] = nd
    m.nodes = nodes
    m.node_history = History(f(i) for i in m.node_history)
    for g in m.groups.values():
        for e in g.elements.values():
            e.nodes = [f(n) if n else 0 for n in e.nodes]
    for cs in m.csys.values():
        if cs.nodes:
            cs.nodes = tuple(f(n) for n in cs.nodes)
    for name in ("forces", "moments"):
        tab = {}
        for k, ld in getattr(m, name).items():
            ld.node = f(k)
            tab[ld.node] = ld
        setattr(m, name, tab)
    m.tmass = {f(k): v for k, v in m.tmass.items()}
    m.rmass = {f(k): v for k, v in m.rmass.items()}
    m.mass_units = {f(k): v for k, v in m.mass_units.items()}
    for name in ("force_history", "moment_history", "tmass_history", "rmass_history"):
        setattr(m, name, History(f(i) for i in getattr(m, name)))
    for r in m.nout:
        r.nodes = [f(n) for n in r.nodes]
    for r in m.rdnd:
        r.node = f(r.node)
    for _, rec in m.options.entries("SYMM"):
        for k in (3, 4, 5):
            v = rec.integer(k)
            if v:
                rec.set_arg(k, f(v))


def ncom(m: SSIModel, notes: Optional[Notes] = None) -> Dict[int, int]:
    """Renumber the nodes 1..N without gaps, keeping their relative order, and remap every node
    reference (:func:`renumber_nodes`; multiple-excitation ME node ranges map to the new numbers of
    the first and last existing nodes of the range).  Returns the map old -> new."""
    notes = _notes(notes)
    dangling = sorted({n for g in m.groups.values() for e in g.elements.values() for n in e.nodes
                       if n and n not in m.nodes})
    if dangling:
        raise GenerationError(f"elements reference undefined nodes {_short(dangling)} (CHECK Error 41): "
                              f"correct them before NCOM")
    ids = sorted(m.nodes)
    mp = {old: new for new, old in enumerate(ids, start=1)}
    # ME node ranges (order preserving map: first/last existing node of the range)
    me_updates = []
    for key, rec in m.options.entries("ME"):
        a, b = rec.integer(2), rec.integer(3)
        if not a or not b:
            continue
        lo, hi = min(a, b), max(a, b)
        i = bisect.bisect_left(ids, lo)
        j = bisect.bisect_right(ids, hi) - 1
        if i <= j:
            me_updates.append((rec, mp[ids[i]], mp[ids[j]]))
        else:
            notes.warn(f"ME {key}: no node in the range {lo}..{hi}; range kept")
    renumber_nodes(m, mp)
    for rec, a, b in me_updates:
        rec.set_arg(2, a)
        rec.set_arg(3, b)
    return mp


# ======================================================================================
# GCOM, MERGEGROUP (spec 09 sections 2.13, 4.3)
# ======================================================================================
def gcom(m: SSIModel) -> Dict[int, int]:
    """Renumber the groups 1..G keeping their order (titles and elements move with them); EOUT
    requests follow.  Element numbers are not changed (ECOMPR).  Returns the map old -> new."""
    mp = {old: new for new, old in enumerate(sorted(m.groups), start=1)}
    groups = {}
    for old in sorted(m.groups):
        g = m.groups[old]
        g.id = mp[old]
        groups[g.id] = g
    m.groups = groups
    if m.group_active is not None:
        m.group_active = mp.get(m.group_active, m.group_active)
    for r in m.eout:
        r.group = mp.get(r.group, r.group)
    return mp


def mergegroup(m: SSIModel, dest: int, others: Sequence[int], notes: Optional[Notes] = None) \
        -> Dict[Tuple[int, int], Tuple[int, int]]:
    """Merge groups into ``dest`` (spec 09 section 4.3): destination elements keep their numbers,
    the elements of each added group (argument order) are numbered after the current maximum, and
    the added groups are removed (GCOM closes the gaps).  All groups must have the same type.
    Returns the element map ``(group, element) -> (dest, new element)``."""
    notes = _notes(notes)
    if dest not in m.groups:
        raise GenerationError(f"group {dest} is not defined")
    seen: List[int] = []
    for g in others:
        if g == dest:
            notes.warn(f"group {g} is the destination; ignored")
            continue
        if g in seen:
            notes.warn(f"group {g} given twice; merged once")
            continue
        if g not in m.groups:
            raise GenerationError(f"group {g} is not defined")
        if m.groups[g].type != m.groups[dest].type:
            raise GenerationError(f"group {g} is {m.groups[g].type_name}, group {dest} is "
                                  f"{m.groups[dest].type_name}: all groups must have the same element type")
        seen.append(g)
    D = m.groups[dest]
    nxt = max(D.elements, default=0) + 1
    emap: Dict[Tuple[int, int], Tuple[int, int]] = {}
    for g in seen:
        for e in m.groups[g].sorted_elements():
            D.elements[nxt] = e.copy(new_id=nxt)
            emap[(g, e.id)] = (dest, nxt)
            nxt += 1
        del m.groups[g]
        if m.group_active == g:
            m.group_active = dest
    removed = set(seen)
    for r in m.eout:
        if r.group in removed:
            old = r.group
            r.elements = [emap[(old, e)][1] for e in r.elements if (old, e) in emap]
            r.group = dest
    m.eout = [r for r in m.eout if r.elements]
    return emap


# ======================================================================================
# WELD (spec 09 section 4.9, D-MDL-15)
# ======================================================================================
def _constrained_clusters(n: int, pairs: np.ndarray, ends: Dict[int, Set[int]]) -> np.ndarray:
    """Union-find over the index ``pairs`` (processed in ascending (i, j) order) with cannot-link
    constraints: ``ends[i]`` holds the ids of the constraints (protected springs) that point i is an
    end of; a union is skipped when it would put both ends of one constraint in the same cluster.
    Returns the cluster root (lowest index) of every point."""
    parent = list(range(n))
    cons: Dict[int, Set[int]] = {i: set(c) for i, c in ends.items()}     # root -> constraint ids inside

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]:
        a, b = find(int(i)), find(int(j))
        if a == b:
            continue
        ca, cb = cons.get(a, set()), cons.get(b, set())
        if ca & cb:                       # both ends of a protected spring would coincide: keep apart
            continue
        lo, hi = (a, b) if a < b else (b, a)
        parent[hi] = lo
        if ca or cb:
            cons[lo] = ca | cb
            cons.pop(hi, None)
    return np.array([find(i) for i in range(n)], dtype=np.int64)


def coincident_groups(ids: Sequence[int], xyz: np.ndarray, tol: float,
                      cannot_link: Sequence[Tuple[int, int]] = ()) -> Dict[int, int]:
    """Substitution table of coincident points: every point within ``tol`` of another (transitively)
    maps to the lowest id of its cluster.

    ``cannot_link`` lists id pairs that must stay in different clusters (the two ends of a
    zero-length SPRING, D-MDL-15).  Clusters touching such a pair are built by a union-find over
    the coincident pairs in ascending id order that skips any union joining the two ends of a
    constrained pair; every other coincident point is still welded (to the lower-numbered
    cluster it meets first)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    order = np.argsort(np.asarray(ids))
    ids = np.asarray(ids)[order]
    X = np.asarray(xyz, float)[order]
    n = len(ids)
    if n < 2:
        return {}
    pairs = cKDTree(X).query_pairs(r=tol, output_type="ndarray")
    if len(pairs) == 0:
        return {}
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    # constrained clusters: components that contain both ends of a cannot-link pair
    pos = {int(i): k for k, i in enumerate(ids)}
    ends: Dict[int, Set[int]] = {}
    for c, (a, b) in enumerate(cannot_link):
        ka, kb = pos.get(int(a)), pos.get(int(b))
        if ka is None or kb is None or ka == kb or lab[ka] != lab[kb]:
            continue
        ends.setdefault(ka, set()).add(c)
        ends.setdefault(kb, set()).add(c)
    if ends:
        comps = {int(lab[k]) for k in ends}
        mask = np.isin(lab[pairs[:, 0]], list(comps))
        root = _constrained_clusters(n, pairs[mask], ends)
        # unconstrained components keep the plain connected-component labels (offset to stay disjoint)
        lab = np.where(np.isin(lab, list(comps)), root, n + lab)
    _, first = np.unique(lab, return_index=True)       # ids are sorted: first index = lowest id
    rep_of = dict(zip(np.unique(lab).tolist(), ids[first].tolist()))
    sub = {}
    for k in range(n):
        r = int(rep_of[int(lab[k])])
        if r != int(ids[k]):
            sub[int(ids[k])] = r
    return sub


def weld(m: SSIModel, force: bool = False, notes: Optional[Notes] = None) -> Dict[int, int]:
    """Merge coincident nodes in the element connectivity, keeping the lowest number (D-MDL-15).

    The substituted nodes stay in the node list, unattached (RMVUNUSED removes them); masses,
    fixities and interaction flags are not changed (manual).  A zero-length SPRING element (FIXROT
    drilling springs, MERGESOIL interface springs) is never collapsed unless ``force``: its two end
    nodes are kept in different clusters (cannot-link), but any *other* node at the same point is
    still welded to one of them -- e.g. the matching node of a second mesh joins the lower-numbered
    spring end, so no unintended crack is left (D-MDL-15 forbids collapsing the spring only).
    Returns the substitution table ``{old: kept}``.
    """
    notes = _notes(notes)
    v = ModelView(m)
    ids = [i for i in sorted(v.P) if np.all(np.isfinite(v.P[i]))]
    if len(ids) < 2:
        return {}
    keep_apart: List[Tuple[int, int]] = []
    if not force:
        for r in v.elems:
            if r.type == SPRING:
                ns = [n for n in r.elem.nodes if n]
                if len(ns) == 2 and v.defined(ns) and ns[0] != ns[1] \
                        and float(np.linalg.norm(v.P[ns[0]] - v.P[ns[1]])) <= v.tol:
                    keep_apart.append((ns[0], ns[1]))
    xyz = np.array([v.P[i] for i in ids])
    sub = coincident_groups(ids, xyz, v.tol, keep_apart)
    if keep_apart:
        notes.info(f"{len(keep_apart)} zero-length SPRING elements (FIXROT / MERGESOIL interface springs) kept: "
                   f"their two end nodes are not welded together (D-MDL-15); other coincident nodes are welded "
                   f"to one of the ends (WELD,FORCE would collapse the springs too)")
    if not sub:
        return {}
    degenerate = collapsed = 0
    for g, e in m.iter_elements():
        old = list(e.nodes)
        new = [sub.get(n, n) if n else 0 for n in old]
        if new != old:
            if len(set(n for n in new if n)) < len(set(n for n in old if n)):
                if g.type == SPRING:
                    collapsed += 1
                else:
                    degenerate += 1
            e.nodes = new
    if degenerate:
        notes.warn(f"{degenerate} elements now have repeated nodes that were distinct before (check the "
                   f"degenerate shapes)")
    if collapsed:
        notes.warn(f"{collapsed} SPRING elements now connect a node to itself")
    moved = [n for n in sub if m.nodes[n].flags or any(m.nodes[n].fix) or n in m.tmass or n in m.rmass
             or n in m.forces or n in m.moments]
    if moved:
        notes.warn(f"{len(moved)} substituted nodes carry masses, loads, fixities or INT flags, which are not "
                   f"transferred: {_short(sorted(moved))}")
    return sub


# ======================================================================================
# ROTATE, TRANSLATE (spec 09 sections 4.6, 4.8; D-MDL-14)
# ======================================================================================
def rotation_matrix(rxy: float, ryz: float, rzx: float) -> np.ndarray:
    """ROTATE matrix (D-MDL-14): Rz(rxy) first, then Rx(ryz), then Ry(rzx), right-hand positive about
    the fixed global axes: ``T = Ry(rzx) Rx(ryz) Rz(rxy)``."""
    return rot_y(rzx) @ rot_x(ryz) @ rot_z(rxy)


def _set_loc(cs: CoordSys, origin: np.ndarray, R: np.ndarray) -> None:
    """Store a LOC system from origin and axes so that WRITE -> INP reproduces it exactly: the
    angles are recomputed and ``R`` is rebuilt from them (as LOC does)."""
    a, b, c = euler_from_matrix(R)
    cs.params = [float(origin[0]), float(origin[1]), float(origin[2]), a, b, c]
    cs.origin = np.array(cs.params[:3])
    cs.R = loc_matrix(a, b, c)


def _set_local_points(cs: CoordSys, pts: np.ndarray) -> None:
    cs.points = [[float(x) for x in p] for p in pts]
    cs.origin = np.array(cs.points[0])
    cs.R = local_from_points(*cs.points)


def rotate_model(m: SSIModel, center: Sequence[float], rxy: float, ryz: float, rzx: float,
                 notes: Optional[Notes] = None) -> int:
    """Rotate the model about ``center`` (D-MDL-14): ``x' = c + T (x - c)``.  Local coordinate
    systems rotate too, so nodes stored in them move with them; beam K nodes rotate with all
    nodes, so beam orientations are kept.  Global SC spring constants, 2-node GENERAL matrices,
    F/MM load factors, D fixities and MT/MR masses are global-axis quantities and are NOT rotated:
    a warning names the data whose physical meaning changes (:func:`global_axis_data_changed`).
    Returns the number of nodes moved."""
    notes = _notes(notes)
    T = rotation_matrix(rxy, ryz, rzx)
    c = np.asarray(center, float)
    for cs in m.csys.values():
        if cs.kind == "LOCAL" and len(cs.points) == 3:
            _set_local_points(cs, (np.asarray(cs.points, float) - c) @ T.T + c)
        else:
            _set_loc(cs, T @ (np.asarray(cs.origin, float) - c) + c, T @ np.asarray(cs.R, float))
            if cs.kind == "LOCAL":
                cs.kind = "LOC"
    k = 0
    for n in m.nodes.values():
        k += 1
        if n.csys != 0:
            continue
        p = T @ (np.array(n.xyz) - c) + c
        n.x, n.y, n.z = (float(p[0]), float(p[1]), float(p[2]))
    used_sc = {e.prop for g, e in m.iter_elements() if g.type == SPRING}
    aniso = [p for p in sorted(used_sc) if p in m.springs and
             (len(set(m.springs[p].k[:3])) > 1 or len(set(m.springs[p].k[3:])) > 1)]
    if aniso and not np.allclose(T, np.eye(3)):
        notes.warn(f"SC spring properties {_short(aniso)} are global and anisotropic: they are not rotated")
    gm2 = sum(1 for g, e in m.iter_elements() if g.type == GENERAL and len([n for n in e.nodes if n]) < 3)
    if gm2:
        notes.warn(f"{gm2} GENERAL elements use global-axis matrices (2 nodes): they are not rotated")
    if (m.forces or m.moments) and not np.allclose(T, np.eye(3)):
        notes.warn("F/MM load factors are global and are not rotated")
    fixed, masses = global_axis_data_changed(m, T)
    if fixed:
        notes.warn(f"{len(fixed)} nodes have D fixities on global axes that restrain other physical directions "
                   f"after the rotation (e.g. a FIXROT drilling restraint ROTZ becomes a bending restraint): "
                   f"fixities are not rotated; re-run FIXROT or redefine D after ROTATE: {_short(fixed)}")
    if masses:
        notes.warn(f"{len(masses)} nodes have direction-dependent MT/MR masses (global axes): they are not "
                   f"rotated: {_short(masses)}")
    if ryz % 360 or rzx % 360:
        notes.info("the ground elevation is not changed by ROTATE")
    return k


def _axes_invariant(T: np.ndarray, weights: Sequence[float], atol: float = 1e-9) -> bool:
    """True when the global-axis diagonal tensor ``diag(weights)`` is unchanged by the rotation:
    ``T D T^T == D``.  For 0/1 fixity flags this says that the span of the restrained axes is
    mapped onto itself (all-or-none always is; ``UZ`` alone is for a rotation about Z); for masses
    that the anisotropic mass looks the same in the rotated model."""
    w = np.asarray(weights, float)
    D = np.diag(w)
    scale = max(1.0, float(np.max(np.abs(w))))
    return bool(np.allclose(T @ D @ T.T, D, rtol=0.0, atol=atol * scale))


def global_axis_data_changed(m: SSIModel, T: np.ndarray) -> Tuple[List[int], List[int]]:
    """Nodes whose global-axis data would change meaning under the rotation ``T`` (ROTATE keeps
    them on the global axes): D fixities (translation and rotation sets tested separately, rotations
    transform like vectors under a proper rotation) and MT/MR masses with unequal components.
    Returns ``(fixity nodes, mass nodes)``."""
    if np.allclose(T, np.eye(3), rtol=0.0, atol=1e-12):
        return [], []
    fixed = [i for i in sorted(m.nodes) if any(m.nodes[i].fix)
             and not (_axes_invariant(T, [1.0 if f else 0.0 for f in m.nodes[i].fix[:3]])
                      and _axes_invariant(T, [1.0 if f else 0.0 for f in m.nodes[i].fix[3:6]]))]
    masses = sorted({i for tab in (m.tmass, m.rmass) for i, val in tab.items()
                     if len(val) >= 3 and not _axes_invariant(T, val[:3])})
    return fixed, masses


def translate_model(m: SSIModel, d: Sequence[float], notes: Optional[Notes] = None) -> int:
    """Move every node by ``d`` (spec 09 section 4.8); nodes stored in a local system get
    ``x_loc += R^T d``, so the systems themselves do not move.  No other data changes (the ground
    elevation is not shifted).  Returns the number of nodes moved."""
    notes = _notes(notes)
    d = np.asarray(d, float)
    k = 0
    for n in m.nodes.values():
        if n.csys == 0:
            dl = d
        else:
            cs = m.csys.get(n.csys)
            if cs is None:
                notes.warn(f"node {n.id} is stored in the undefined system {n.csys}; not moved")
                continue
            dl = np.asarray(cs.R, float).T @ d
        n.x, n.y, n.z = (n.x + float(dl[0]), n.y + float(dl[1]), n.z + float(dl[2]))
        k += 1
    if d[2] != 0:
        notes.warn("the ground elevation is not shifted: check GROUNDELEV")
    return k


# ======================================================================================
# MERGE, MERGESOIL (spec 09 sections 4.2, 4.5; spec 04 section 7; D-MDL-12)
# ======================================================================================
@dataclass
class AppendMaps:
    """Number maps of a model appended by :func:`append_model`."""
    node: Dict[int, int] = field(default_factory=dict)
    group: Dict[int, int] = field(default_factory=dict)
    material: Dict[int, int] = field(default_factory=dict)
    section: Dict[int, int] = field(default_factory=dict)
    spring: Dict[int, int] = field(default_factory=dict)
    matrix: Dict[int, int] = field(default_factory=dict)
    csys: Dict[int, int] = field(default_factory=dict)
    node_offset: int = 0
    group_offset: int = 0


def _top(tab) -> int:
    return max(tab) if tab else 0


def _layer_same(a: SoilLayer, b: SoilLayer, thickness: bool = True) -> bool:
    """Same soil-layer data (relative 1e-9); ``thickness=False`` compares the material data only."""
    fields = SoilLayer.FIELDS if thickness else SoilLayer.FIELDS[1:]
    return all(math.isclose(getattr(a, f), getattr(b, f), rel_tol=1e-9, abs_tol=1e-12) for f in fields)


def append_model(dest: SSIModel, src: SSIModel, shift: Sequence[float] = (0.0, 0.0, 0.0),
                 mset: Optional[Callable[[Group, Element], int]] = None,
                 materials: Optional[Set[int]] = None, notes: Optional[Notes] = None) -> AppendMaps:
    """Append a copy of ``src`` to ``dest`` with numbering offsets (MERGE, D-MDL-12).

    Node, group, material (M), section (R), spring (SC), matrix (MX) and coordinate-system numbers
    of ``src`` are offset by the maxima of ``dest`` (new numbers start at max + 1); element numbers
    stay (the groups are new).  The L table is **not** offset (far-field layers are shared): missing
    layers are added, a same-number layer that differs keeps the ``dest`` value with a warning.
    ``src`` nodes are translated by ``shift`` (global nodes directly, local systems by their
    origin).  Element MSET is mapped by ``mset(group, element)`` (default: + material offset),
    RSET by group type (BEAMS -> R, SPRING -> SC, GENERAL -> MX).  Fixities and INT flags,
    masses, loads and the NOUT/EOUT/RDND requests are copied and remapped.  ``materials`` limits
    the copied M entries (MERGESOIL).  Analysis options, frequency sets and lists stay those of
    ``dest``.
    """
    notes = _notes(notes)
    d = np.asarray(shift, float)
    moved = bool(np.any(d != 0.0))
    mp = AppendMaps()
    noff, goff = _top(dest.nodes), _top(dest.groups)
    moff, roff, soff, xoff, coff = (_top(dest.materials), _top(dest.sections), _top(dest.springs),
                                    _top(dest.matrices), _top(dest.csys))
    mp.node_offset, mp.group_offset = noff, goff

    def nf(n: int) -> int:
        return n + noff if n else 0

    mp.node = {i: i + noff for i in src.nodes}
    mp.csys = {0: 0}
    mp.csys.update({s: s + coff for s in src.csys})
    # coordinate systems
    for sid in sorted(src.csys):
        cs = CoordSys.from_json(src.csys[sid].to_json())
        cs.id = mp.csys[sid]
        if cs.kind == "LOCAL":
            cs.nodes = tuple(nf(n) for n in cs.nodes)
        if moved:
            if cs.kind == "LOCAL" and len(cs.points) == 3:
                _set_local_points(cs, np.asarray(cs.points, float) + d)
            elif cs.kind == "LOC" and len(cs.params) >= 3:
                p = list(cs.params) + [0.0] * (6 - len(cs.params))
                p[0], p[1], p[2] = p[0] + d[0], p[1] + d[1], p[2] + d[2]
                cs.params = p
                cs.origin = np.array(p[:3])
            else:
                cs.origin = np.asarray(cs.origin, float) + d
        dest.csys[cs.id] = cs
    # nodes
    for i in sorted(src.nodes):
        n = src.nodes[i]
        xyz = np.array(n.xyz)
        if n.csys == 0 and moved:
            xyz = xyz + d
        dest.nodes[mp.node[i]] = Node(mp.node[i], float(xyz[0]), float(xyz[1]), float(xyz[2]),
                                      mp.csys.get(n.csys, n.csys), list(n.fix), set(n.flags))
    for i in src.node_history:
        if i in mp.node:
            dest.node_history.touch(mp.node[i])
    # property tables
    for k in sorted(src.materials):
        if materials is not None and k not in materials:
            continue
        s = src.materials[k]
        mp.material[k] = k + moff
        dest.materials[k + moff] = Material(k + moff, s.val1, s.val2, s.weight, s.pdamp, s.sdamp, s.mtype)
    for k in sorted(src.sections):
        s = src.sections[k]
        mp.section[k] = k + roff
        dest.sections[k + roff] = BeamSection(k + roff, s.axial, s.shear2, s.shear3, s.tors, s.flex2, s.flex3)
    for k in sorted(src.springs):
        s = src.springs[k]
        mp.spring[k] = k + soff
        dest.springs[k + soff] = SpringProp(k + soff, s.scx, s.scy, s.scz, s.scxx, s.scyy, s.sczz, s.damp)
    for k in sorted(src.matrices):
        p = MatrixProp.from_json(src.matrices[k].to_json())
        p.id = k + xoff
        mp.matrix[k] = p.id
        dest.matrices[p.id] = p
    for k in sorted(src.layers):
        s = src.layers[k]
        if k in dest.layers:
            if not _layer_same(dest.layers[k], s):
                notes.warn(f"soil layer {k} differs between the models; the first model's layer is kept "
                           f"(the L table is shared, D-MDL-12)")
            continue
        dest.layers[k] = SoilLayer(k, s.thick, s.weight, s.vp, s.vs, s.pdamp, s.sdamp)
    # groups and elements
    for gid in sorted(src.groups):
        g = src.groups[gid]
        ng = Group(gid + goff, g.type, g.title)
        mp.group[gid] = ng.id
        for e in g.sorted_elements():
            ne = e.copy(nodes=[nf(n) for n in e.nodes])
            ne.mat = mset(g, e) if mset is not None else e.mat + moff
            if g.type == BEAMS:
                ne.prop = e.prop + roff
            elif g.type == SPRING:
                ne.prop = e.prop + soff
            elif g.type == GENERAL:
                ne.prop = e.prop + xoff
            ng.elements[ne.id] = ne
        dest.groups[ng.id] = ng
    # loads and masses
    for name, hist in (("forces", "force_history"), ("moments", "moment_history")):
        tab = getattr(dest, name)
        for k, ld in getattr(src, name).items():
            tab[nf(k)] = NodalLoad(nf(k), list(ld.factor), list(ld.arrival))
        for k in getattr(src, hist):
            getattr(dest, hist).touch(nf(k))
    for name, hist in (("tmass", "tmass_history"), ("rmass", "rmass_history")):
        tab = getattr(dest, name)
        for k, val in getattr(src, name).items():
            tab[nf(k)] = list(val)
        for k in getattr(src, hist):
            getattr(dest, hist).touch(nf(k))
    for k, val in src.mass_units.items():
        dest.mass_units[nf(k)] = val
    # output requests
    for r in src.nout:
        dest.nout.append(NodalRequest(r.dir, list(r.codes), [nf(n) for n in r.nodes]))
    for r in src.eout:
        dest.eout.append(ElementRequest(list(r.codes), r.group + goff, list(r.elements)))
    for r in src.rdnd:
        dest.rdnd.append(RelDispRequest(nf(r.node), list(r.flags)))
    return mp


def merge_models(m1: SSIModel, m2: SSIModel, shift: Sequence[float] = (0.0, 0.0, 0.0),
                 notes: Optional[Notes] = None) -> Tuple[SSIModel, AppendMaps]:
    """MERGE,<Mdl1>,<Mdl2>,<X>,<Y>,<Z> (spec 09 section 4.2, D-MDL-12): a copy of Mdl1 with Mdl2
    appended by :func:`append_model` and translated by (X, Y, Z).  No coincident nodes are merged
    (WELD / MERGESOIL).  Excavated SOLID/PLANE elements of Mdl2 keep their MSET (L table, shared);
    other elements get the material offset."""
    notes = _notes(notes)
    out = m1.copy()
    v2 = ModelView(m2)
    exc2 = {(r.group, r.id) for r in v2.elems if r.excavated}
    moff = _top(out.materials)
    maps = append_model(out, m2, shift, mset=lambda g, e: e.mat if (g.id, e.id) in exc2 else e.mat + moff,
                        notes=notes)
    # implicit ETYPE 0 elements classified differently in the merged model (other gelev, shift)
    v = ModelView(out)
    newg = set(maps.group.values())
    changed = sum(1 for r in v.elems if r.group in newg and r.type in (SOLID, PLANE) and r.elem.etype == 0
                  and ((r.group - maps.group_offset, r.id) in exc2) != r.excavated)
    if changed:
        notes.warn(f"{changed} SOLID/PLANE elements of the second model with implicit ETYPE 0 are classified "
                   f"differently in the merged model (ground elevation {v.gelev:g}); set ETYPE explicitly")
    check_node_limit(out, notes)
    return out, maps


@dataclass
class MergeSoilResult:
    model: SSIModel
    maps: AppendMaps
    pairs: List[Tuple[int, int]]                 # (soil node, structure node), soil numbering of the result
    node_map: Dict[int, int]                     # soil-model node -> node in the result (Mapping file)
    springs: int = 0
    spring_group: Optional[int] = None


def _soil_layer_from_material(mat: Material, gravity: float, thick: float) -> SoilLayer:
    """L entry equivalent to an M entry: ``Vs = sqrt(G/rho)``, ``Vp = sqrt(M/rho)`` with
    ``rho = weight/g``; damping ratios as entered (spec 04 section 7.2)."""
    k = elastic_constants(mat.mtype, mat.val1, mat.val2, mat.weight, gravity)
    return SoilLayer(mat.id, thick, mat.weight, k.Vp, k.Vs, mat.pdamp, mat.sdamp)


def mergesoil(struct: SSIModel, soil: SSIModel, mode: int = 1, stiff: float = 1.0e7, stiff2: float = 10.0,
              seplevel: Optional[float] = None, notes: Optional[Notes] = None) -> MergeSoilResult:
    """MERGESOIL (spec 09 section 4.5, spec 04 section 7.2, D-MDL-12).

    1. The structure model is copied unchanged; the soil (excavation) model is appended as by MERGE
       (no translation).  Every soil SOLID/PLANE/SHELL/TSHELL element gets ETYPE 2 (excavated soil /
       buried shell) and the MSET of the soil SOLID/PLANE elements become soil layers, whatever
       their ETYPE in the soil model: MSET m is the soil model's L entry m when defined, otherwise
       its M entry m converted to L layer m (``Vs = sqrt(G/rho)``, ``Vp = sqrt(M/rho)``; the manual
       workflow INP, ETYPEGEN,2, MERGESOIL of spec 04 section 7.2); a warning when both exist and
       differ, an error when neither does.  A layer m already in the structure model is kept (with a
       warning when it differs, D-MDL-12).
    2. Interface pairs: each soil node s on the foundation-soil interface of the excavation (the
       lateral and bottom boundary faces of the soil SOLID/PLANE elements, at or below grade
       ``z <= gelev + tol``) and a structure node t at the same location (``|x_s - x_t| <= tol``,
       ``tol`` the larger of the two models' geometric tolerances -- an EXCAV model with delta > 0
       carries ``EDUOPT,GEOMTOL = delta``, which is then also set on the merged model; used structure
       nodes first, the lowest exactly coincident number, else the nearest).  Interior excavation
       nodes and interior ground-surface nodes are never joined, even when they coincide with
       structure nodes (floor slabs): excavated soil connects to the structure only at the
       foundation-soil interface (spec 01 rules 4 and 11; the EXCSTRCHK condition).  Interface nodes
       without a structure partner are reported (warning in Modes 1-3, spec 01 rules 8 and 11).
    3. Mode 0 unbonded (nothing), 1 merge (s replaced by the lower-numbered t in the soil elements;
       the interaction flag of s moves to t; s becomes unused), 2 stiff SPRING elements (s, t)
       ``scx = scy = scz = Stiff`` (rotations 0, damping 0) in a new SPRING group, 3 as 2 with
       ``Stiff2`` above ``SepLevel`` (``z <= SepLevel + tol`` counts as below).
    """
    notes = _notes(notes)
    if mode not in (0, 1, 2, 3):
        raise GenerationError("[Mode] must be 0 (unbonded), 1 (merge nodes), 2 (stiff springs) or 3 (stiff below "
                              "SepLevel, soft above)")
    if mode == 3 and seplevel is None:
        raise GenerationError("Mode 3 needs [SepLevel]")
    if mode in (2, 3) and (stiff <= 0 or (mode == 3 and stiff2 <= 0)):
        raise GenerationError("spring stiffnesses must be > 0")
    if not soil.groups:
        raise GenerationError("the soil model has no elements")
    vq = ModelView(soil)
    # materials of the soil SOLID/PLANE elements -> soil layers, layer number = MSET (D-MDL-12):
    # whatever the element ETYPE, MSET m is the soil model's L entry m when it exists (EXCAV models:
    # excavation layers = far-field layers), otherwise its M entry m converted to a layer (models
    # built with structural materials, e.g. converted from ANSYS, then ETYPEGEN,2: spec 04 section
    # 7.2 "their materials become soil layers").  Neither defined -> error.
    layer_src: Dict[int, Tuple[str, object]] = {}
    keep_mats: Set[int] = set()
    zext: Dict[int, List[float]] = {}
    for g, e in soil.iter_elements():
        if g.type in (SOLID, PLANE):
            if e.mat not in layer_src:
                lay, mat = soil.layers.get(e.mat), soil.materials.get(e.mat)
                if lay is None and mat is None:
                    raise GenerationError(f"soil model group {g.id} element {e.id}: MSET {e.mat} is neither a "
                                          f"soil layer (L) nor a material (M) of the soil model")
                if lay is not None:
                    layer_src[e.mat] = ("L", lay)
                    if mat is not None:
                        try:
                            same = _layer_same(_soil_layer_from_material(mat, soil.gravity, lay.thick), lay,
                                               thickness=False)
                        except ValueError:
                            same = False
                        if not same:
                            notes.warn(f"soil model: MSET {e.mat} is both soil layer {e.mat} and material {e.mat}, "
                                       f"which differ; the soil layer is used (D-MDL-12)")
                else:
                    layer_src[e.mat] = ("M", mat)
            zs = [vq.P[n][2] for n in e.nodes if n in vq.P]
            if zs:
                zext.setdefault(e.mat, []).extend([min(zs), max(zs)])
        elif g.type in MATERIAL_TYPES:
            keep_mats.add(e.mat)
    others = sorted({g.type_name for g in soil.groups.values() if g.type not in (SOLID, PLANE)})
    if others:
        notes.warn(f"the soil model also contains {', '.join(others)} elements (SHELL/TSHELL become buried "
                   f"shells)")
    implicit_s = sum(1 for g, e in struct.iter_elements() if g.type in (SOLID, PLANE) and e.etype == 0)
    if implicit_s:
        notes.warn(f"the structure model has {implicit_s} SOLID/PLANE elements with implicit ETYPE 0 (use "
                   f"ETYPEGEN,1 on the structure model)")
    out = struct.copy()
    gelev = out.ground_elevation
    if not math.isclose(soil.ground_elevation, gelev, rel_tol=0.0, abs_tol=1e-9 * max(1.0, abs(gelev))):
        notes.warn(f"ground elevations differ (structure {gelev:g}, soil {soil.ground_elevation:g}); the "
                   f"structure's value is used")
    moff = _top(out.materials)

    def mset(g: Group, e: Element) -> int:
        return e.mat if g.type in (SOLID, PLANE) else e.mat + moff

    maps = append_model(out, soil, mset=mset, materials=keep_mats, notes=notes)
    newg = set(maps.group.values())
    for gid in newg:
        g = out.groups[gid]
        if g.type in ETYPE_TYPES:
            for e in g.elements.values():
                e.etype = 2
    # soil layers from the soil materials
    for lid in sorted(layer_src):
        kind, src = layer_src[lid]
        if kind == "L":
            lay = src                                    # already merged by append_model
            if lid not in out.layers:
                out.layers[lid] = SoilLayer(lid, lay.thick, lay.weight, lay.vp, lay.vs, lay.pdamp, lay.sdamp)
            continue
        z = zext.get(lid, [0.0, 1.0])
        thick = (max(z) - min(z)) if max(z) > min(z) else 1.0
        try:
            lay = _soil_layer_from_material(src, soil.gravity, thick)     # type: ignore[arg-type]
        except ValueError as exc:
            raise GenerationError(f"soil material {lid}: {exc}") from None
        if lid in out.layers:
            if not _layer_same(lay, out.layers[lid], thickness=False):
                notes.warn(f"soil material {lid} differs from the existing soil layer {lid}, which is kept "
                           f"(D-MDL-12)")
        else:
            out.layers[lid] = lay
    # foundation-soil interface of the soil model (lateral + bottom faces, gelev of the structure)
    qq = soil.copy()
    for g in qq.groups.values():
        if g.type in (SOLID, PLANE):
            for e in g.elements.values():
                e.etype = 2
    qq.options.ensure_record("HOUSE").set_arg(2, gelev)
    exq = ModelView(qq).excavation()
    iface = set(exq["fsin"]) if exq["elements"] else set(soil.nodes)
    iface |= {n for g, e in soil.iter_elements() if g.type not in (SOLID, PLANE) for n in e.nodes if n}
    # matching tolerance (D-GEN-06, D-MDL-12 "coincidence = geometric tol"): the coarser of the two
    # models' tolerances.  An excavation generated by EXCAV with delta > 0 carries EDUOPT,GEOMTOL =
    # delta (its levels are cluster means of structure nodes scattered by up to delta, spec 09
    # section 4.1); that tolerance is passed on to the merged model, whose interface nodes then
    # scatter in the same way (INTGEN must still see the ground-surface face).
    tol0 = ModelView(out).tol
    soil_gt = eduopt_float(soil, "GEOMTOL", 0.0)
    if soil_gt > tol0:
        set_geomtol(out, soil_gt)
        notes.info(f"interface nodes matched within {soil_gt:g} (EDUOPT,GEOMTOL of the soil model, e.g. the EXCAV "
                   f"delta); EDUOPT,GEOMTOL,{fmt_num(soil_gt)} set on the merged model")
    v = ModelView(out)
    rad = max(v.tol, vq.tol)
    # interface pairs (spatial search)
    struct_ids = sorted(i for i in struct.nodes if i in v.P)
    used = {n for r in v.elems if r.group not in newg for n in r.elem.nodes if n}
    pairs: List[Tuple[int, int]] = []
    inner = 0
    if struct_ids:
        from scipy.spatial import cKDTree
        SP = np.array([v.P[i] for i in struct_ids])
        tree = cKDTree(SP)
        for old in sorted(soil.nodes):
            s = maps.node[old]
            p = v.P.get(s)
            if p is None or not np.all(np.isfinite(p)) or p[2] > gelev + v.tol:
                continue
            hit = tree.query_ball_point(p, r=rad)
            if hit and old not in iface:
                inner += 1
            elif hit:
                # used structure nodes first; among them the exactly coincident (tol0) lowest number,
                # otherwise the nearest one
                def key(k: int) -> Tuple[bool, bool, float, int]:
                    d = float(np.linalg.norm(SP[k] - p))
                    return struct_ids[k] not in used, d > tol0, 0.0 if d <= tol0 else d, struct_ids[k]
                pairs.append((s, struct_ids[min(hit, key=key)]))
    if inner:
        notes.info(f"{inner} soil nodes inside the excavation volume or on its ground-surface face coincide with "
                   f"structure nodes and are not joined (only the foundation-soil interface is, spec 01 rule 11)")
    # completeness of the interface (spec 01 rules 8 and 11): every foundation-soil interface node of
    # the excavation at or below grade should find a structure node, otherwise the structure is not
    # connected to the excavated soil there
    paired = {s for s, _ in pairs}
    need = sorted(maps.node[o] for o in exq["fsin"]
                  if maps.node[o] in v.P and v.P[maps.node[o]][2] <= gelev + v.tol)
    missing = [s for s in need if s not in paired]
    if missing:
        what = (f"none of the {len(need)} foundation-soil interface nodes of the soil model has" if not pairs
                else f"{len(missing)} of the {len(need)} foundation-soil interface nodes of the soil model have")
        text = (f"{what} a coincident structure node (tolerance {rad:.3g}): ")
        if mode == 0:
            notes.info(text + f"they stay unconnected anyway (Mode 0): {_short(missing)}")
        elif not pairs:
            notes.warn(text + "the structure is NOT connected to the excavated soil (spec 01 rules 8 and 11); check "
                       "the coordinates of both models, GROUNDELEV and EDUOPT,GEOMTOL (the EXCAV delta)")
        else:
            notes.warn(text + f"the structure is not connected to the excavated soil there (spec 01 rules 8 and "
                       f"11); merged-model nodes {_short(missing)}")
    node_map = {old: maps.node[old] for old in soil.nodes}
    res = MergeSoilResult(out, maps, pairs, node_map)
    if mode == 1:
        sub = {s: t for s, t in pairs}
        for gid in newg:
            for e in out.groups[gid].elements.values():
                e.nodes = [sub.get(n, n) if n else 0 for n in e.nodes]
        loaded = []
        for s, t in pairs:
            ns, nt = out.nodes[s], out.nodes[t]
            if 0 in ns.flags:
                nt.flags.add(0)
                ns.flags.discard(0)
            if s in out.tmass or s in out.rmass or s in out.forces or s in out.moments:
                loaded.append(s)
        if loaded:
            notes.warn(f"{len(loaded)} merged soil nodes carry masses or loads (they stay on the unused node)")
        for r in out.nout[len(struct.nout):]:
            r.nodes = [sub.get(n, n) for n in r.nodes]
        for r in out.rdnd[len(struct.rdnd):]:
            r.node = sub.get(r.node, r.node)
        inv = {v_: k for k, v_ in maps.node.items()}
        for s, t in pairs:
            node_map[inv[s]] = t
    elif mode in (2, 3):
        gid = _top(out.groups) + 1
        grp = Group(gid, SPRING, "MERGESOIL interface springs")
        p_hi = _top(out.springs) + 1
        out.springs[p_hi] = SpringProp(p_hi, stiff, stiff, stiff, 0.0, 0.0, 0.0, 0.0)
        p_lo = None
        if mode == 3:
            p_lo = p_hi + 1
            out.springs[p_lo] = SpringProp(p_lo, stiff2, stiff2, stiff2, 0.0, 0.0, 0.0, 0.0)
        for k, (s, t) in enumerate(pairs, start=1):
            prop = p_hi
            if mode == 3 and v.P[s][2] > seplevel + v.tol:   # z = SepLevel (within tol) counts as below (D-MDL-12)
                prop = p_lo
            grp.elements[k] = Element(k, [s, t], mat=1, prop=prop)
        if pairs:
            out.groups[gid] = grp
            res.spring_group = gid
        else:
            del out.springs[p_hi]
            if p_lo is not None:
                del out.springs[p_lo]
        res.springs = len(pairs)
    check_node_limit(out, notes)
    return res


def write_node_map(path, node_map: Dict[int, int], title: str = "") -> None:
    """Mapping file of MERGESOIL / NCOM: one ``old new`` line per node (spec 09 section 4.5 step 5)."""
    lines = [f"# {title}" if title else "# node map", "# old new"]
    lines += [f"{o} {n}" for o, n in sorted(node_map.items())]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


# ======================================================================================
# EXCAV (spec 09 section 4.1; D-MDL-19)
# ======================================================================================
def _structural_refs(v: ModelView) -> List[ElemRef]:
    return [r for r in v.elems if r.type in (SOLID, PLANE, SHELL, TSHELL, BEAMS) and not r.excavated]


#: faces within about 0.8 degrees of horizontal are template faces (EXCAV / SOILMESH)
HORIZONTAL_COS = 0.9999


def horizontal_faces(v: ModelView, refs: Sequence[ElemRef], allowed: Set[int]) -> List[List[int]]:
    """Nearly horizontal faces (Newell normal with ``|n_z| >= HORIZONTAL_COS``) of SHELL/TSHELL
    elements and of SOLID faces, all of whose nodes are in ``allowed``: ordered distinct nodes.

    The test is geometric, not by level, so a basemat whose nodes have small z variations (to be
    gathered by EXCAV ``delta``) is still recognised as the template grid."""
    from .check import newell_normal
    out = []
    for r in refs:
        ns = r.elem.nodes
        if r.type in (SHELL, TSHELL):
            polys = [_ordered_unique(ns)]
        elif r.type == SOLID:
            pad = list(ns) + [0] * (8 - len(ns))
            polys = [_ordered_unique(pad[k - 1] for k in f) for f in SOLID_FACES]
        else:
            continue
        for p in polys:
            if len(p) < 3 or not all(n in allowed for n in p):
                continue
            nv = newell_normal([v.P[n] for n in p])
            ln = float(np.linalg.norm(nv))
            if ln > 0 and abs(nv[2]) >= HORIZONTAL_COS * ln:
                out.append(p)
    return out


def lowest_faces(v: ModelView, faces: Sequence[List[int]], thr: float) -> List[List[int]]:
    """The lowest band of horizontal faces (the basemat grid): faces sorted by their lowest node
    elevation join the band while that elevation is within ``thr`` of the band's highest node."""
    spans = sorted((min(v.P[n][2] for n in f), max(v.P[n][2] for n in f), k) for k, f in enumerate(faces))
    top = spans[0][1]
    out = []
    for lo, hi, k in spans:
        if lo > top + thr:
            break
        out.append(faces[k])
        top = max(top, hi)
    return out


@dataclass
class TemplateGrid:
    """Plan template of EXCAV / SOILMESH: unique (x, y) points and CCW cells referencing them."""
    xy: np.ndarray
    source: List[int]                 # lowest source node at each point
    cells: List[Tuple[int, ...]]      # CCW (seen from above), 3 or 4 point indices


def template_grid(v: ModelView, faces: Sequence[List[int]], notes: Notes) -> TemplateGrid:
    """Unique plan points of the faces (coincident in plan within tol) ordered by their lowest source
    node number, and unique cells oriented counter-clockwise seen from above."""
    nodes = sorted({n for f in faces for n in f})
    xy = np.array([v.P[n][:2] for n in nodes])
    lab = unique_points(xy, v.tol)
    first: Dict[int, int] = {}
    for n, l in zip(nodes, lab):
        first.setdefault(int(l), n)
    order = sorted(first, key=lambda l: first[l])
    idx = {l: k for k, l in enumerate(order)}
    pidx = {n: idx[int(l)] for n, l in zip(nodes, lab)}
    pts = np.array([v.P[first[l]][:2] for l in order])
    cells: Dict[Tuple[int, ...], Tuple[int, ...]] = {}
    degenerate = 0
    for f in faces:
        c = list(dict.fromkeys(pidx[n] for n in f))
        if len(c) < 3 or len(c) > 4:
            degenerate += 1
            continue
        a = signed_area(pts[c])
        if abs(a) <= v.tol * v.tol:
            degenerate += 1
            continue
        if a < 0:
            c = c[::-1]
        k0 = c.index(min(c))                    # start at the lowest point index (deterministic)
        cells.setdefault(tuple(sorted(c)), tuple(c[k0:] + c[:k0]))
    if degenerate:
        notes.warn(f"{degenerate} template faces with zero plan area or more than 4 corners skipped")
    return TemplateGrid(pts, [first[l] for l in order], [cells[k] for k in sorted(cells)])


def _solid_nodes(bottom: Sequence[int], top: Sequence[int]) -> List[int]:
    """SOLID node list of a cell extruded between two levels: hexahedron for a quad, prism
    ``a b c c a' b' c' c'`` for a triangle (repeated nodes, spec 09 section 1.3)."""
    if len(bottom) == 4:
        return list(bottom) + list(top)
    a, b, c = bottom
    a2, b2, c2 = top
    return [a, b, c, c, a2, b2, c2, c2]


@dataclass
class ExcavResult:
    model: SSIModel
    levels: List[float]
    template_points: int
    cells: int
    elements: int
    nodes: int


def excav(src: SSIModel, delta: float = 0.0, notes: Optional[Notes] = None) -> ExcavResult:
    """EXCAV,<model>,[delta]: excavation volume of the active model's basement (spec 09 section 4.1,
    D-MDL-19).

    1. Nodes of structural SOLID/SHELL/TSHELL/BEAMS/PLANE elements with ``z <= gelev + thr`` are
       clustered into levels (``thr = delta`` if > 0, else the geometric tolerance; a level starts
       when ``z - z_start > thr``; level elevation = cluster mean).  The level containing ``gelev``
       is snapped to it, otherwise ``gelev`` is added as the top level.
    2. Template: the lowest band of (nearly) horizontal shell / solid faces, the basemat grid; the
       levels start at its lowest node level (parts of the basement that do not reach it get no
       excavation volume; with ``delta`` too small, the z scatter of the grid gives extra levels).
    3. Every template cell is extruded between consecutive levels into a SOLID (hexahedron, or prism
       for a triangle); nodes are numbered bottom-up, level by level, at the template (x, y) and the
       level elevation.
    4. One SOLID group, ETYPE 2; MSET = the TOPL soil layer containing the element mid-height
       (MSET 1 with a warning when no TOPL profile is defined).  The soil layers, TOPL and the
       HOUSE gravity and ground elevation are copied.
    5. ``delta`` > tol: the structure nodes lie up to ``delta`` off the (snapped) level elevations,
       so the excavation model gets ``EDUOPT,GEOMTOL = delta``; MERGESOIL matches the interface with
       that tolerance (D-GEN-06), so the structure is still joined to the excavation.
    """
    notes = _notes(notes)
    v = ModelView(src)
    gelev = v.gelev
    thr = float(delta) if delta and delta > 0 else v.tol
    refs = _structural_refs(v)
    cand = sorted({n for r in refs for n in r.dof_nodes if n in v.P and v.P[n][2] <= gelev + thr})
    if not cand:
        raise GenerationError(f"no structural node at or below the ground elevation {gelev:g}")
    z = np.array([v.P[n][2] for n in cand])
    means, lab = cluster_values(z, thr)
    level = {n: int(l) for n, l in zip(cand, lab)}
    faces = horizontal_faces(v, refs, set(cand))
    if not faces:
        raise GenerationError("no horizontal SHELL or SOLID face below grade: the basemat grid defines the "
                              "EXCAV template")
    base = lowest_faces(v, faces, thr)
    l0 = min(level[n] for f in base for n in f)
    if l0 > 0:
        notes.warn(f"{l0} node levels below the lowest basemat grid (z = {means[l0]:g}) ignored")
    levels = [float(x) for x in means[l0:]]
    if abs(levels[-1] - gelev) <= thr:
        levels[-1] = gelev
    elif levels[-1] < gelev:
        levels.append(gelev)
    if len(levels) < 2:
        raise GenerationError("the basemat lies at the ground elevation: there is no embedment")
    tg = template_grid(v, base, notes)
    nT = len(tg.xy)
    out = SSIModel()
    for k, zk in enumerate(levels):
        for i in range(nT):
            out.define_node(k * nT + i + 1, (float(tg.xy[i, 0]), float(tg.xy[i, 1]), zk), 0)
    grp = Group(1, SOLID, "EXCAV excavation volume")
    eid = 0
    no_topl = below = 0
    for k in range(len(levels) - 1):
        zmid = 0.5 * (levels[k] + levels[k + 1])
        li, deeper = _layer_index(v, gelev - zmid)
        if li is None:
            mset, no_topl = 1, no_topl + 1
        else:
            mset = src.topl[li]
            below += int(deeper)
        for c in tg.cells:
            eid += 1
            ns = _solid_nodes([k * nT + i + 1 for i in c], [(k + 1) * nT + i + 1 for i in c])
            grp.elements[eid] = Element(eid, ns, mat=mset, etype=2)
    out.groups[1] = grp
    out.group_active = 1
    for k in sorted(src.layers):
        s = src.layers[k]
        out.layers[k] = SoilLayer(k, s.thick, s.weight, s.vp, s.vs, s.pdamp, s.sdamp)
    out.topl = list(src.topl)
    rec = out.options.ensure_record("HOUSE")
    rec.set_arg(1, src.gravity)
    rec.set_arg(2, gelev)
    if delta and delta > v.tol:
        # the generated nodes sit on the level means, the structure nodes up to delta off them
        # (a little more at a cluster snapped to gelev): record that scatter as the model's geometric
        # tolerance so that MERGESOIL joins them (spec 09 section 4.1 step 5, "MERGESOIL later welds
        # coincident nodes to the structure")
        dev = max((abs(float(v.P[n][2]) - levels[level[n] - l0]) for n in cand if level[n] >= l0), default=0.0)
        gt = max(float(delta), dev * (1.0 + 1e-9))
        set_geomtol(out, gt)
        notes.info(f"EDUOPT,GEOMTOL,{fmt_num(gt)} set on the excavation model (the z scatter allowed by delta): "
                   f"MERGESOIL joins structure nodes within that distance of the generated levels")
    if no_topl:
        notes.warn("no TOPL soil-layer profile: MSET 1 used for every excavation element")
    if below:
        notes.warn(f"{below} element layers lie below the last TOPL layer: its layer number is used")
    off = _off_interface_levels(v, levels)
    if off:
        notes.warn(f"levels {', '.join(f'{z:g}' for z in off)} are not soil-layer interfaces (interaction nodes "
                   f"must lie on interfaces, EDU-01)")
    return ExcavResult(out, levels, nT, len(tg.cells), eid, nT * len(levels))


# ======================================================================================
# SOILMESH (spec 09 section 4.7; D-MDL-19)
# ======================================================================================
def perimeter(xy: np.ndarray, tol: float) -> List[int]:
    """Indices of the points on the convex hull boundary (collinear edge points included), in
    counter-clockwise order starting from the hull's first vertex."""
    p = np.asarray(xy, float)
    h = convex_hull(p)
    if len(h) < 3:
        raise GenerationError("the basement perimeter is degenerate (fewer than 3 hull vertices)")
    out: List[Tuple[float, int]] = []
    s0 = 0.0
    taken: Set[int] = set()
    for k in range(len(h)):
        a, b = p[h[k]], p[h[(k + 1) % len(h)]]
        ab = b - a
        L = float(np.hypot(*ab))
        t = ((p - a) @ ab) / (L * L)
        dist = np.abs(ab[0] * (p[:, 1] - a[1]) - ab[1] * (p[:, 0] - a[0])) / L
        on = np.where((dist <= tol) & (t >= -tol / L) & (t < 1.0 - tol / L))[0]
        for i in on:
            if int(i) not in taken:
                taken.add(int(i))
                out.append((s0 + float(t[i]) * L, int(i)))
        s0 += L
    return [i for _, i in sorted(out)]


@dataclass
class SoilMeshResult:
    model: SSIModel
    levels: List[float]
    perimeter_points: int
    footprint_cells: int
    elements: int
    new_nodes: int
    springs: int
    group: int
    spring_group: Optional[int] = None


def soilmesh(src: SSIModel, sx: float, sy: float, hori: int, vert: int, xadj: float = 0.0, yadj: float = 0.0,
             zdepth: float = 0.0, contact: int = 0, rnum: int = 1, notes: Optional[Notes] = None) -> SoilMeshResult:
    """SOILMESH: near-field soil mesh around the basement for soil-pressure models (spec 09
    section 4.7, D-MDL-19; reconstruction, the manual gives only the parameters).

    * Basement levels: nodes of structural elements at or below grade, clustered with tol.  The
      external wall is the convex-hull perimeter of each level, which must be the same at every
      level (uniform external wall); centroid ``c`` = perimeter mean + (xAdj, yAdj).
    * ``hori`` rings: perimeter point p of ring i lies at ``c + (p - c)(1 + i s/100)`` (linear growth
      with sX in x and sY in y) at every basement level; hexahedra join rings i-1 and i between
      consecutive levels.
    * ``vert`` layers of thickness ``Zdepth`` below the basemat cover the rings and the basemat
      footprint (template of the lowest level's faces, as EXCAV).
    * Soil SOLIDs: ETYPE 1 (near-field soil is part of the structure), material = a new M entry
      (type 3, Vp/Vs/weight/damping) of the TOPL layer at the element mid-height.
    * ``contact`` = 0: the soil shares the basement wall and mat nodes; otherwise duplicate soil-side
      nodes are tied to the basement by SPRING elements of property ``rNum`` (contact surface).
    * Numbering starts after the source maxima (nodes, groups, materials); the result model holds
      only the new entities plus copies of the referenced basement nodes, so it can be written
      (WRITE) and read over the original model (INP).
    """
    notes = _notes(notes)
    if hori < 0 or vert < 0 or (hori == 0 and vert == 0):
        raise GenerationError("<hori> and <vert> must be >= 0 and not both 0")
    if vert > 0 and zdepth <= 0:
        raise GenerationError("<Zdepth> must be > 0 when <vert> > 0")
    if 1 + hori * sx / 100.0 <= 0 or 1 + hori * sy / 100.0 <= 0:
        raise GenerationError("<sX>/<sY> make a ring collapse")
    v = ModelView(src)
    gelev, tol = v.gelev, v.tol
    depths = v.interfaces()
    if depths is None:
        raise GenerationError("SOILMESH needs the TOPL soil-layer profile (materials of the soil elements)")
    refs = _structural_refs(v)
    cand = sorted({n for r in refs for n in r.dof_nodes if n in v.P and v.P[n][2] <= gelev + tol})
    if not cand:
        raise GenerationError(f"no structural node at or below the ground elevation {gelev:g}")
    means, lab = cluster_values([v.P[n][2] for n in cand], tol)
    level = {n: int(l) for n, l in zip(cand, lab)}
    faces = horizontal_faces(v, refs, set(cand))
    if not faces:
        raise GenerationError("no basemat (horizontal SHELL or SOLID faces) below grade")
    bfaces = lowest_faces(v, faces, tol)
    l0 = min(level[n] for f in bfaces for n in f)
    levels = [float(x) for x in means[l0:]]
    nL = len(levels)
    if nL < 2 and hori > 0:
        notes.warn("the basement has a single level: the rings have no height")
    # perimeter of the bottom level (reference) and its nodes at every level
    at = {k: [n for n in cand if level[n] == k + l0] for k in range(nL)}
    xy0 = np.array([v.P[n][:2] for n in at[0]])
    per0 = perimeter(xy0, tol)
    P = len(per0)
    pxy = xy0[per0]
    wall: List[List[int]] = []
    for k in range(nL):
        ns = at[k]
        xyk = np.array([v.P[n][:2] for n in ns])
        row = []
        for p in pxy:
            d = np.hypot(*(xyk - p).T)
            j = int(np.argmin(d))
            if d[j] > max(tol, 1e-9):
                raise GenerationError(f"no basement node at ({p[0]:g}, {p[1]:g}) on level z = {levels[k]:g}: "
                                      f"SOILMESH needs a uniform external wall at all levels")
            row.append(min(n for n, dd in zip(ns, d) if dd <= max(tol, 1e-9)))
        wall.append(row)
        perk = perimeter(xyk, tol)
        if len(perk) != P:
            notes.warn(f"level z = {levels[k]:g}: the perimeter has {len(perk)} nodes, the basemat {P} "
                       f"(non-uniform external wall)")
    c = pxy.mean(axis=0) + np.array([xadj, yadj])
    out = SSIModel()
    nid = _top(src.nodes)
    copied: Set[int] = set()

    def base(n: int) -> int:
        if n not in copied:
            copied.add(n)
            p = v.P[n]
            out.define_node(n, (float(p[0]), float(p[1]), float(p[2])), 0)
        return n

    dup: Dict[int, int] = {}

    def soil_side(n: int) -> int:
        nonlocal nid
        if not contact:
            return base(n)
        if n not in dup:
            base(n)
            nid += 1
            p = v.P[n]
            out.define_node(nid, (float(p[0]), float(p[1]), float(p[2])), 0)
            dup[n] = nid
        return dup[n]

    def new_node(x: float, y: float, z: float) -> int:
        nonlocal nid
        nid += 1
        out.define_node(nid, (x, y, z), 0)
        return nid

    # ring nodes ring[i][k][j] (i = 0: the wall, soil side)
    ring: List[List[List[int]]] = [[[0] * P for _ in range(nL)] for _ in range(hori + 1)]
    for k in range(nL):
        for j in range(P):
            ring[0][k][j] = soil_side(wall[k][j]) if hori > 0 else 0
    rxy = [pxy]
    for i in range(1, hori + 1):
        f = np.array([1.0 + i * sx / 100.0, 1.0 + i * sy / 100.0])
        rxy.append(c + (pxy - c) * f)
        for k in range(nL):
            for j in range(P):
                ring[i][k][j] = new_node(float(rxy[i][j, 0]), float(rxy[i][j, 1]), levels[k])
    # materials: one M entry per TOPL layer used
    moff = _top(src.materials)
    lay_mat: Dict[int, int] = {}

    def material(zbot: float, ztop: float) -> int:
        li, _ = _layer_index(v, gelev - 0.5 * (zbot + ztop))
        lid = src.topl[li]
        if lid not in lay_mat:
            L = src.layers[lid]
            mid = moff + len(lay_mat) + 1
            out.materials[mid] = Material(mid, L.vp, L.vs, L.weight, L.pdamp, L.sdamp, 3)
            lay_mat[lid] = mid
        return lay_mat[lid]

    gid = _top(src.groups) + 1
    grp = Group(gid, SOLID, "SOILMESH near-field soil")
    eid = 0

    def add_solid(bottom: List[int], top: List[int], zb: float, zt: float, xy_b: np.ndarray) -> None:
        nonlocal eid
        if signed_area(xy_b) < 0:
            bottom, top = bottom[::-1], top[::-1]
        eid += 1
        grp.elements[eid] = Element(eid, _solid_nodes(bottom, top), mat=material(zb, zt), etype=1)

    for i in range(1, hori + 1):
        for k in range(nL - 1):
            for j in range(P):
                j2 = (j + 1) % P
                q = [(i - 1, j), (i, j), (i, j2), (i - 1, j2)]
                add_solid([ring[a][k][b] for a, b in q], [ring[a][k + 1][b] for a, b in q], levels[k], levels[k + 1],
                          np.array([rxy[a][b] for a, b in q]))
    fcells = 0
    if vert > 0:
        tg = template_grid(v, bfaces, notes)
        fcells = len(tg.cells)
        # bottom level nodes of the footprint (basemat nodes, soil side) and of the rings
        src_at_bottom = {}
        for t, n0 in enumerate(tg.source):
            p = tg.xy[t]
            ns = [n for n in at[0] if np.hypot(*(v.P[n][:2] - p)) <= max(tol, 1e-9)]
            src_at_bottom[t] = min(ns) if ns else n0
        zb = [levels[0] - j * zdepth for j in range(vert, 0, -1)] + [levels[0]]   # bottom-up
        col: Dict[Tuple[str, int, int], List[int]] = {}
        for t in range(len(tg.xy)):
            col[("f", 0, t)] = [new_node(float(tg.xy[t, 0]), float(tg.xy[t, 1]), z) for z in zb[:-1]] + \
                               [soil_side(src_at_bottom[t])]
        # ring 0 of the footprint = the template points of the perimeter
        tidx = {}
        for j in range(P):
            d = np.hypot(*(tg.xy - pxy[j]).T)
            tidx[j] = int(np.argmin(d))
        for i in range(1, hori + 1):
            for j in range(P):
                col[("r", i, j)] = [new_node(float(rxy[i][j, 0]), float(rxy[i][j, 1]), z) for z in zb[:-1]] + \
                                   [ring[i][0][j]]
        for j in range(P):
            col[("r", 0, j)] = col[("f", 0, tidx[j])]
        for k in range(vert):
            for cell in tg.cells:
                add_solid([col[("f", 0, t)][k] for t in cell], [col[("f", 0, t)][k + 1] for t in cell], zb[k],
                          zb[k + 1], tg.xy[list(cell)])
            for i in range(1, hori + 1):
                for j in range(P):
                    j2 = (j + 1) % P
                    q = [(i - 1, j), (i, j), (i, j2), (i - 1, j2)]
                    add_solid([col[("r", a, b)][k] for a, b in q], [col[("r", a, b)][k + 1] for a, b in q], zb[k],
                              zb[k + 1], np.array([rxy[a][b] for a, b in q]))
    if not grp.elements:
        raise GenerationError("no soil element generated")
    out.groups[gid] = grp
    out.group_active = gid
    res = SoilMeshResult(out, levels, P, fcells, eid, nid - _top(src.nodes), 0, gid)
    if contact and dup:
        sg = gid + 1
        sgrp = Group(sg, SPRING, "SOILMESH contact springs")
        for k, (b, s) in enumerate(sorted(dup.items()), start=1):
            sgrp.elements[k] = Element(k, [b, s], mat=1, prop=rnum)
        out.groups[sg] = sgrp
        res.springs, res.spring_group = len(dup), sg
        if rnum in src.springs:
            s = src.springs[rnum]
            out.springs[rnum] = SpringProp(rnum, s.scx, s.scy, s.scz, s.scxx, s.scyy, s.sczz, s.damp)
        else:
            notes.warn(f"spring property {rnum} (contact surface) is not defined: define SC,{rnum}")
    off = _off_interface_levels(v, levels)
    if off:
        notes.warn(f"basement levels {', '.join(f'{z:g}' for z in off)} are not soil-layer interfaces")
    return res
