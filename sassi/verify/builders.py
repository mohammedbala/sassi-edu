"""Builders of canonical SSI models as module decks, and runners of the SITE -> POINT -> HOUSE ->
(FORCE) -> ANALYS -> MOTION chain (shared helper library of the verification problems and examples).

Why this module exists
----------------------
A verification problem of the SSI chain needs the same steps every time: a layered site, the point
loads of POINT, a HOUSE model (structure, excavated soil, interaction nodes), possibly FORCE loads,
an ANALYS run and a MOTION run.  In SASSI-EDU every module reads one *deck* written by AFWRITE
(:mod:`sassi.io.decks`); this module writes those decks directly from small Python descriptions, so
a verification problem is a few lines and exercises exactly the deck -> module -> file path of a
real run (``sassi.modules.base.run_module``, the function behind RUNSITE, RUNANALYS ...).

Contents
--------
* Site and frequencies: :class:`SoilLayer`, :class:`Site` (``layered_site``, ``uniform_site``,
  ``Site.subdivided``), :class:`FrequencySet` (Fourier grid ``df = 1/(NFFT dt)`` or a harmonic step).
* Decks: :func:`site_deck`, :func:`point_deck`, :func:`force_deck`, :func:`analys_deck`,
  :func:`motion_deck` and the general HOUSE deck builder :class:`HouseBuilder`.
* Canonical models (each returns a :class:`Model`: the HOUSE builder plus named nodes and the
  recommended POINT parameters):

  - :func:`surface_rigid_mat` -- square or circular surface mat of interaction nodes made rigid by
    stiff massless BEAMS to a centre node (or by stiff SHELL elements); optional foundation mass;
  - :func:`embedded_box` -- a box excavation through the top TOPL layers (excavated SOLID elements
    with the layer properties) with the interaction sets of the methods FV, FI-FSIN, FI-EVBN, FFV
    and a structure that is 'none', 'soil' (structure = excavated soil, the zero-SSI identity of
    VP-16), 'solid' (an embedded block) or 'shell' (basement walls, base and roof slabs);
  - :func:`stick_on_mat` -- a lumped-mass BEAMS stick on a rigid surface mat;
  - :func:`sdof_on_node` -- a spring-mass oscillator on one surface interaction node (SPRING, or the
    equivalent GENERAL element in global or local axes), the "fixed-base" SDOF of VP-01 / VP-39 on
    a very stiff site;
  - 2D (plane strain, HOUSE ``<dim>`` = 1, POINT2): :func:`rigid_strip_2d` (rigid surface strip,
    VP-42) and :func:`embedded_plane_2d` (PLANE excavation with an optional PLANE structure and stick,
    VP-T3); pass ``dim=mdl.dim`` (1) to :func:`run_soil` / :func:`point_deck` for POINT2;
  - SYMM planes: :meth:`HouseBuilder.symmetry` writes the HOUSE deck ``symm`` rows (half / quarter
    models, D-ANL-12; see sassi/verify/problems/vp_twod_symm.py).

* Runners: :func:`run` (one module, raises :class:`ChainError` with the listing tail on failure),
  :func:`run_soil` (SITE + POINT), :func:`run_site_xyz` (three SITE runs SV x' / SH y' / P z' and the
  FILE1X/Y/Z copies of D-ANL-06), :func:`run_house`, :func:`run_force`, :func:`run_analys`,
  :func:`run_motion`; result readers :func:`read_file8`, :func:`tf`.

Conventions: consistent units, ``Site.gravity`` converts weights to masses (rho = weight / g); the
ground surface is at ``Site.gelev`` (default 0, z up); user interface i is the top of TOPL layer i;
node ids are assigned bottom-up so interaction nodes are numbered bottom-up (EDU-21).  The central
zone radius defaults to R1 3.4: R0 = 0.90 h (square meshes) or 0.85 h (circular meshes).

Example (a surface mat on a layered site; units m, t, kN, s)::

    from sassi.verify import builders as B
    site = B.layered_site([(2.0, 150.0, 300.0, 1.9, 0.05)], (400.0, 800.0, 2.1, 0.03))
    fs = B.FrequencySet.fourier(0.01, 1024, [4, 20, 60])          # df = 1/(NFFT dt)
    mdl = B.surface_rigid_mat(site, half_width=4.0, ndiv=4)
    B.run_soil(wd, "m", site, fs, layer=mdl.layer, rad=mdl.rad)   # SITE (vertical SV) + POINT
    B.run_house(wd, "m", mdl)                                      # HOUSE
    B.run_analys(wd, "m", fs, impe=2)                              # ANALYS + global impedance
    H = B.tf(B.read_file8(wd), mdl["centre"], 1)                   # centre-node ATF, x
    # rigid-foundation compliance instead: unit loads at the centre node in a vibration run
    B.run_force(wd, "m", [(mdl["centre"], 1, 1.0, 0.0)], fs)
    B.run_analys(wd, "m", fs, type=1)
"""
from __future__ import annotations

import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..io import decks
from ..io.container import Container, read_container
from ..io.deckfmt import Table
from ..modules.base import run_module

GRAVITY = 9.81
#: wave name -> (SITE WAVE type, control direction <cm>, wave family <wopt>)
WAVES = {"SV": (2, 0, 0), "SH": (4, 1, 1), "P": (3, 2, 0)}
#: default stiffness factor of "rigid" links relative to the largest soil shear modulus
RIGID_FACTOR = 1.0e5
PathLike = Union[str, Path]


class ChainError(RuntimeError):
    """A module of the chain failed; the message holds the end of its listing."""


# ======================================================================================
# Site and frequencies
# ======================================================================================
@dataclass
class SoilLayer:
    """One soil layer: thickness, velocities, mass density and damping ratios (fractions)."""

    thick: float
    vs: float
    vp: float
    rho: float
    beta_s: float = 0.05
    beta_p: Optional[float] = None

    @property
    def dp(self) -> float:
        return self.beta_s if self.beta_p is None else self.beta_p

    @property
    def G(self) -> float:
        return self.rho * self.vs ** 2

    @classmethod
    def from_nu(cls, thick: float, vs: float, rho: float, nu: float = 1.0 / 3.0, beta: float = 0.05,
                beta_p: Optional[float] = None) -> "SoilLayer":
        """Layer with Vp from Poisson's ratio: ``Vp = Vs sqrt(2(1 - nu)/(1 - 2 nu))``."""
        vp = vs * math.sqrt(2.0 * (1.0 - nu) / (1.0 - 2.0 * nu))
        return cls(thick, vs, vp, rho, beta, beta_p)


@dataclass
class Site:
    """TOPL layers (top first) on a half-space (``nl`` generated sublayers; ``nl = 0``: rigid base)."""

    layers: List[SoilLayer]
    halfspace: SoilLayer
    gravity: float = GRAVITY
    nl: int = 20
    cmodform: int = 0
    hslaw: Optional[str] = None
    gelev: float = 0.0

    @property
    def nI(self) -> int:
        """Number of user interfaces (TOPL layers + 1)."""
        return len(self.layers) + 1

    def depths(self) -> np.ndarray:
        """Depths of the user interfaces 1..nI (interface i = top of TOPL layer i)."""
        return np.concatenate([[0.0], np.cumsum([l.thick for l in self.layers])])

    def elevations(self) -> np.ndarray:
        return self.gelev - self.depths()

    def layer_rows(self, halfspace_thick: float = 0.0) -> List[List[float]]:
        """Rows ``[no, thick, weight, vp, vs, dp, ds]`` of the TOPL layers and the half-space."""
        rows = []
        for k, l in enumerate(self.layers + [self.halfspace]):
            t = l.thick if k < len(self.layers) else halfspace_thick
            rows.append([k + 1, float(t), l.rho * self.gravity, float(l.vp), float(l.vs), float(l.dp),
                         float(l.beta_s)])
        return rows

    def max_G(self) -> float:
        return max(l.G for l in self.layers + [self.halfspace])

    def subdivided(self, max_thick: float) -> "Site":
        """The same site with every layer split into equal sublayers no thicker than ``max_thick``."""
        out = []
        for l in self.layers:
            n = max(1, int(math.ceil(l.thick / max_thick - 1e-12)))
            out += [SoilLayer(l.thick / n, l.vs, l.vp, l.rho, l.beta_s, l.beta_p)] * n
        return Site(out, self.halfspace, self.gravity, self.nl, self.cmodform, self.hslaw, self.gelev)


def layered_site(layers: Sequence[Sequence[float]], halfspace: Sequence[float], gravity: float = GRAVITY,
                 nl: int = 20, cmodform: int = 0, hslaw: Optional[str] = None, gelev: float = 0.0) -> Site:
    """Site from tuples ``(thick, vs, vp, rho, beta_s[, beta_p])`` (TOPL layers, top first) and
    ``(vs, vp, rho, beta_s[, beta_p])`` for the half-space."""
    ls = [SoilLayer(*[float(v) for v in row]) for row in layers]
    hs = SoilLayer(0.0, *[float(v) for v in halfspace])
    return Site(ls, hs, gravity, nl, cmodform, hslaw, gelev)


def uniform_site(depth: float, nsub: int, vs: float, rho: float, nu: float = 1.0 / 3.0, beta: float = 0.05,
                 hs_vs: Optional[float] = None, hs_beta: Optional[float] = None, gravity: float = GRAVITY,
                 nl: int = 20, cmodform: int = 0) -> Site:
    """A uniform stratum of ``nsub`` equal TOPL layers on a half-space (by default the same soil:
    a uniform half-space)."""
    lay = SoilLayer.from_nu(depth / nsub, vs, rho, nu, beta)
    hs = SoilLayer.from_nu(0.0, hs_vs or vs, rho, nu, beta if hs_beta is None else hs_beta)
    return Site([lay] * nsub, hs, gravity, nl, cmodform)


@dataclass
class FrequencySet:
    """SSI frequencies ``f = n df`` (requirements 4.0.1).

    Build with :meth:`fourier` (time-history analyses: ``df = 1/(NFFT dt)``, needed by MOTION) or
    :meth:`harmonic` (a free frequency step, ``fstep``; SITE/POINT/FORCE/ANALYS only)."""

    fnum: List[int]
    df: float
    delt: float = 0.005
    nft: int = 4096
    fstep: float = 0.0

    @classmethod
    def fourier(cls, delt: float, nft: int, fnum: Iterable[int]) -> "FrequencySet":
        return cls(sorted(int(n) for n in fnum), 1.0 / (delt * nft), float(delt), int(nft), 0.0)

    @classmethod
    def harmonic(cls, df: float, fnum: Iterable[int], nft: Optional[int] = None) -> "FrequencySet":
        """Free frequency step ``df`` (SITE/FORCE ``fstep``).  ``nft`` defaults to the smallest power
        of 2 (>= 4096) with every frequency number below NFFT/2; ``delt = 1/(df NFFT)`` so that a
        MOTION run on the same grid is possible."""
        fn = sorted(int(n) for n in fnum)
        if nft is None:
            nft = 4096
            while fn and nft // 2 < fn[-1]:
                nft *= 2
        return cls(fn, float(df), 1.0 / (df * nft), int(nft), float(df))

    @property
    def freq(self) -> np.ndarray:
        return np.asarray(self.fnum, float) * self.df

    def subset(self, fnum: Iterable[int]) -> "FrequencySet":
        return FrequencySet(sorted(int(n) for n in fnum), self.df, self.delt, self.nft, self.fstep)

    def fill(self, d: decks.Deck) -> decks.Deck:
        """Write delt / nft / df (and fstep where the schema has it) and the frequency table."""
        d["delt"], d["nft"], d["df"] = self.delt, self.nft, self.df
        if "fstep" in d.params and d.module in ("SITE", "FORCE"):
            d["fstep"] = self.fstep
        d.tables["freqs"] = Table(columns=["number"])
        for n in self.fnum:
            d.table("freqs").append([int(n)])
        return d


# ======================================================================================
# Decks of SITE, POINT, FORCE, ANALYS, MOTION
# ======================================================================================
def _common(d: decks.Deck, model: str, title: str) -> decks.Deck:
    d["model"] = model
    d["title"] = title or f"{model} ({d.module})"
    return d


def site_deck(site: Site, fs: FrequencySet, wave: str = "SV", cl: int = 1, cm: Optional[int] = None,
              waves: Optional[Sequence[Sequence[float]]] = None, wopt: Optional[int] = None, mode1: int = 1,
              mode2: int = 1, model: str = "m", title: str = "") -> decks.Deck:
    """SITE deck: vertically incident ``wave`` ('SV' x', 'SH' y' or 'P' z') with the control point at
    the top of TOPL layer ``cl`` (within motion), or explicit WAVE records ``waves``
    ``(type, opt, ratio1, ratio2, angle)`` with ``cm``/``wopt``."""
    d = _common(decks.new("SITE"), model, title)
    d["gravity"] = site.gravity
    d["nl"] = site.nl
    d["cmodform"] = site.cmodform
    if site.hslaw is not None:
        d["hslaw"] = site.hslaw
    d["hs"] = len(site.layers) + 1
    d["mode1"], d["mode2"] = int(mode1), int(mode2)
    d["cl"] = int(cl)
    wt, wcm, wfam = WAVES[wave.upper()]
    d["cm"] = wcm if cm is None else int(cm)
    d["wopt"] = wfam if wopt is None else int(wopt)
    for row in site.layer_rows()[:-1]:
        d.table("layers").append(row)
    d.table("halfspace").append(site.layer_rows()[-1])
    for w in (waves if waves is not None else [(wt, 1, 1.0, 1.0, 0.0)]):
        d.table("waves").append(list(w))
    return fs.fill(d)


def point_deck(fs: FrequencySet, layer: int = 0, rad: float = 1.0, model: str = "m", title: str = "",
               dim: int = 2) -> decks.Deck:
    """POINT deck: point loads at the user interfaces 1..layer+1, central zone radius ``rad``;
    ``dim`` 2 = POINT3 (3D), 1 = POINT2 (2D line loads, D-PNT-01)."""
    d = _common(decks.new("POINT"), model, title)
    d["layer"], d["rad"], d["dim"] = int(layer), float(rad), int(dim)
    d["df"] = fs.df
    for n in fs.fnum:
        d.table("freqs").append([int(n)])
    return d


def force_deck(loads: Sequence[Sequence[float]], fs: FrequencySet, gravity: float = GRAVITY, mforce: int = 1,
               model: str = "m", title: str = "") -> decks.Deck:
    """FORCE deck: rows ``(node, dof, factor, arrival)`` (dof 1-3 forces, 4-6 moments)."""
    d = _common(decks.new("FORCE"), model, title)
    d["gravity"], d["mforce"] = gravity, int(mforce)
    for r in loads:
        d.table("loads").append([int(r[0]), int(r[1]), float(r[2]), float(r[3]) if len(r) > 3 else 0.0])
    return fs.fill(d)


def analys_deck(fs: FrequencySet, gravity: float = GRAVITY, model: str = "m", title: str = "",
                **params) -> decks.Deck:
    """ANALYS deck; ``params`` are deck parameters (type, mode, save, prnt, fopt, ang, xc, yc, zc,
    impe, simul, delrst ...)."""
    d = _common(decks.new("ANALYS"), model, title)
    d["gravity"] = gravity
    for k, v in params.items():
        if k not in d.params:
            raise KeyError(f"unknown ANALYS parameter {k}")
        d[k] = v
    return fs.fill(d)


def motion_deck(fs: FrequencySet, thfile: str = "", nout: Sequence[Sequence[int]] = (), damp: Sequence[float] = (),
                gravity: float = GRAVITY, model: str = "m", title: str = "", **params) -> decks.Deck:
    """MOTION deck on the Fourier grid of ``fs`` (requires ``df = 1/(NFFT dt)``); ``nout`` rows
    ``(node, dof, c1..c6)``; ``params`` are other deck parameters (type, resp, interp, file8 ...)."""
    d = _common(decks.new("MOTION"), model, title)
    d["delt"], d["nft"], d["df"] = fs.delt, fs.nft, 1.0 / (fs.delt * fs.nft)
    d["gravity"] = gravity
    d["thfile"] = thfile
    d["mult"], d["max"] = 1.0, 0.0
    for k, v in params.items():
        if k not in d.params:
            raise KeyError(f"unknown MOTION parameter {k}")
        d[k] = v
    for r in nout:
        d.table("nout").append(list(r))
    for z in damp:
        d.table("damp").append([float(z)])
    return d


def write_history(path: PathLike, values: Sequence[float], dt: float) -> Path:
    """Control-motion / load history file (MOTION <fopt> 0: dt on the first line), 17 digits."""
    lines = [f"{dt:.16g}"] + [f"{v:.16e}" for v in np.asarray(values, float)]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return Path(path)


# ======================================================================================
# HOUSE deck builder
# ======================================================================================
def _spring_pattern(k6: Sequence[float]) -> np.ndarray:
    """12x12 two-node spring stiffness pattern with constants ``k6`` on the six DOFs."""
    k6 = np.asarray(k6, float)
    K = np.zeros((12, 12))
    i = np.arange(6)
    K[i, i] = k6
    K[i + 6, i + 6] = k6
    K[i, i + 6] = -k6
    K[i + 6, i] = -k6
    return K


class HouseBuilder:
    """Accumulates nodes, elements, tables, masses and interaction nodes and writes the HOUSE deck.

    Node ids are assigned in creation order (create the lowest nodes first to get a bottom-up
    numbering).  Tables are filled exactly as AFWRITE does: global coordinates, fixity flags, the
    resolved ETYPE, the M/L/R/SC/MX tables, nodal masses (mass units by default) and the site table
    (``sitelayers``: TOPL rows + half-space row of thickness 0)."""

    def __init__(self, site: Site, dim: int = 2, imp: int = 0, incomp: int = 1, title: str = ""):
        self.site = site
        self.dim, self.imp, self.incomp, self.title = int(dim), int(imp), int(incomp), title
        self.nodes: Dict[int, List[float]] = {}
        self.fixity: Dict[int, List[int]] = {}
        self._coord: Dict[Tuple[int, int, int], int] = {}
        self.groups: List[List] = []
        self.elements: List[List] = []
        self.materials: List[List] = []
        self.sections: List[List] = []
        self.springs: List[List] = []
        self.matrices: List[List] = []
        self.masses: List[List] = []
        self.interaction: List[int] = []
        self.symm: List[List[int]] = []             # SYMM rows (no, type, n1, n2, n3)
        self._eid: Dict[int, int] = {}
        self._mat_cache: Dict[tuple, int] = {}
        self._sec_cache: Dict[tuple, int] = {}
        self._group_cache: Dict[tuple, int] = {}
        self._knodes: Dict[Tuple[int, int, int], int] = {}
        self._scale = 1e-6

    # ---- nodes ----------------------------------------------------------------------------
    def _key(self, x: float, y: float, z: float) -> Tuple[int, int, int]:
        return tuple(int(round(v / self._scale)) for v in (x, y, z))

    def add_node(self, x: float, y: float, z: float, fix: Sequence[int] = (0,) * 6, nid: Optional[int] = None,
                 register: bool = True) -> int:
        """New node (id = next free id unless ``nid``); ``fix`` = six 0/1 flags (UX..ROTZ).
        ``register=False`` keeps the node out of the coordinate lookup of :meth:`node` / :meth:`find`
        (orientation K nodes, which must never be picked up as structural nodes)."""
        if nid is None:
            nid = (max(self.nodes) + 1) if self.nodes else 1
        if nid in self.nodes:
            raise ValueError(f"node {nid} already defined")
        self.nodes[nid] = [float(x), float(y), float(z)]
        self.fixity[nid] = [int(v) for v in fix]
        if register:
            self._coord.setdefault(self._key(x, y, z), nid)
        return nid

    def orientation_node(self, x: float, y: float, z: float) -> int:
        """A K (orientation) node at (x, y, z): fixed, carries no DOF, shared by members that need
        the same point, never returned by :meth:`node`."""
        k = self._key(x, y, z)
        if k not in self._knodes:
            self._knodes[k] = self.add_node(x, y, z, fix=(1,) * 6, register=False)
        return self._knodes[k]

    def node(self, x: float, y: float, z: float, fix: Optional[Sequence[int]] = None) -> int:
        """The node at (x, y, z) (created if absent; coordinates matched to 1e-6)."""
        k = self._key(x, y, z)
        if k in self._coord:
            nid = self._coord[k]
            if fix is not None:
                self.fix(nid, [i + 1 for i, v in enumerate(fix) if v])
            return nid
        return self.add_node(x, y, z, fix if fix is not None else (0,) * 6)

    def find(self, x: float, y: float, z: float) -> Optional[int]:
        return self._coord.get(self._key(x, y, z))

    def fix(self, node: int, dofs: Iterable[int]) -> None:
        """Fix DOFs (1..6) of ``node`` (D command)."""
        for k in dofs:
            self.fixity[node][int(k) - 1] = 1

    def xyz(self, node: int) -> np.ndarray:
        return np.asarray(self.nodes[node])

    # ---- properties ---------------------------------------------------------------------------
    def material(self, mtype: int, val1: float, val2: float, weight: float, pdamp: float = 0.0,
                 sdamp: Optional[float] = None) -> int:
        """M-table material (type 1: E, nu; 2: M, G; 3: Vp, Vs); returns its id."""
        row = (int(mtype), float(val1), float(val2), float(weight), float(pdamp),
               float(pdamp if sdamp is None else sdamp))
        if row not in self._mat_cache:
            mid = len(self.materials) + 1
            self.materials.append([mid, *row])
            self._mat_cache[row] = mid
        return self._mat_cache[row]

    def elastic(self, E: float, nu: float, rho: float = 0.0, beta: float = 0.0) -> int:
        """Type-1 material from E, nu and the mass density rho (weight = rho g)."""
        return self.material(1, E, nu, rho * self.site.gravity, beta, beta)

    def soil_material(self, layer: int) -> int:
        """Type-3 material (Vp, Vs, weight, damping) identical to TOPL layer ``layer`` (1-based):
        a structure made of it has exactly the excavated-soil matrices (VP-16, D-ELM-12)."""
        l = self.site.layers[layer - 1]
        return self.material(3, l.vp, l.vs, l.rho * self.site.gravity, l.dp, l.beta_s)

    def section(self, A: float, As2: float = 0.0, As3: float = 0.0, J: float = 0.0, I2: float = 0.0,
                I3: float = 0.0) -> int:
        row = tuple(float(v) for v in (A, As2, As3, J, I2, I3))
        if row not in self._sec_cache:
            sid = len(self.sections) + 1
            self.sections.append([sid, *row])
            self._sec_cache[row] = sid
        return self._sec_cache[row]

    def group(self, code: int, title: str = "") -> int:
        """Element group of type ``code`` (1 SOLID, 2 BEAMS, 3 SHELL, 7 SPRING, 9 GENERAL); groups
        with the same (code, title) are shared."""
        key = (int(code), title)
        if key not in self._group_cache:
            gid = len(self.groups) + 1
            self.groups.append([gid, int(code), title or f"group {gid}"])
            self._group_cache[key] = gid
            self._eid[gid] = 0
        return self._group_cache[key]

    def _element(self, group: int, nodes: Sequence[int], etype: int = 1, mat: int = 0, prop: int = 0,
                 eint: int = 0, thick: float = 0.0, ki: str = "000000", kj: str = "000000") -> int:
        self._eid[group] += 1
        eid = self._eid[group]
        n = (list(int(v) for v in nodes) + [0] * 8)[:8]
        self.elements.append([group, eid, int(etype), int(mat), int(prop), int(eint), float(thick)] + n + [ki, kj])
        return eid

    # ---- elements --------------------------------------------------------------------------------
    def solid(self, nodes8: Sequence[int], mat: int, etype: int = 1, group: Optional[int] = None,
              eint: int = 0) -> int:
        """SOLID element (nodes 1-4 one face, 5-8 the opposite face); ``etype`` 1 structure (M table
        ``mat``), 2 excavated soil (``mat`` = L layer number, MSET)."""
        g = group or self.group(1, "excavated soil" if etype == 2 else "solid structure")
        return self._element(g, nodes8, etype=etype, mat=mat, eint=eint)

    def plane(self, nodes4: Sequence[int], mat: int, etype: int = 1, group: Optional[int] = None) -> int:
        """PLANE element (2D plane strain in the X-Z plane, nodes I J K L); ``etype`` 1 structure (M table
        ``mat``), 2 excavated soil (``mat`` = L layer number, MSET)."""
        g = group or self.group(4, "excavated soil" if etype == 2 else "plane structure")
        return self._element(g, nodes4, etype=etype, mat=mat)

    def beam(self, i: int, j: int, k: int, mat: int, section: int, group: Optional[int] = None,
             ki: str = "000000", kj: str = "000000") -> int:
        """BEAMS element I-J with orientation node K (axis 2 toward K)."""
        g = group or self.group(2, "beams")
        return self._element(g, (i, j, k), mat=mat, prop=section, ki=ki, kj=kj)

    def shell(self, nodes: Sequence[int], mat: int, thick: float, group: Optional[int] = None, etype: int = 1) -> int:
        g = group or self.group(3, "shells")
        return self._element(g, nodes, etype=etype, mat=mat, thick=thick)

    def spring(self, i: int, j: int, k6: Sequence[float], damp: float = 0.0, group: Optional[int] = None) -> int:
        """SPRING element: six uncoupled global constants ``k6`` with hysteretic damping ``damp``."""
        pid = len(self.springs) + 1
        self.springs.append([pid] + [float(v) for v in k6] + [float(damp)])
        g = group or self.group(7, "springs")
        return self._element(g, (i, j), prop=pid)

    def general(self, nodes: Sequence[int], KR: np.ndarray, KI: Optional[np.ndarray] = None,
                MM: Optional[np.ndarray] = None, group: Optional[int] = None) -> int:
        """GENERAL element: full symmetric 12x12 ``KR + i KI`` and mass ``MM`` (mass units unless
        the deck says otherwise); 2 nodes = global axes, 3 nodes (I, J, K) = local axes."""
        pid = len({r[0] for r in self.matrices}) + 1
        for kind, A in (("R", KR), ("I", KI), ("M", MM)):
            if A is None:
                continue
            A = np.asarray(A, float)
            for r in range(12):
                terms = list(A[r, r:]) + [0.0] * r
                self.matrices.append([pid, kind, r + 1] + [float(t) for t in terms])
        g = group or self.group(9, "general")
        return self._element(g, nodes, prop=pid)

    def mass(self, node: int, mx: float, my: Optional[float] = None, mz: Optional[float] = None,
             mxx: float = 0.0, myy: float = 0.0, mzz: float = 0.0, units: int = 0) -> None:
        """Nodal mass (MT/MR); ``units`` 0 mass, 1 weight."""
        my = mx if my is None else my
        mz = mx if mz is None else mz
        self.masses.append([int(node), float(mx), float(my), float(mz), float(mxx), float(myy), float(mzz),
                            int(units)])

    def set_interaction(self, nodes: Iterable[int]) -> None:
        for n in nodes:
            if int(n) not in self.interaction:
                self.interaction.append(int(n))

    def symmetry(self, no: int, stype: int, nodes: Sequence[int]) -> None:
        """SYMM,<no>,<type>,<n1>,<n2>[,<n3>]: symmetry (0) / antisymmetry (1) plane through the nodes
        (3D: three nodes of a plane parallel to XZ or YZ; 2D: two nodes of a line parallel to Z)."""
        ids = (list(int(v) for v in nodes) + [0, 0, 0])[:3]
        self.symm = [r for r in self.symm if r[0] != int(no)] + [[int(no), int(stype)] + ids]

    # ---- deck --------------------------------------------------------------------------------------
    def interaction_order(self) -> List[int]:
        """Interaction nodes bottom-up (ascending z, then node id; EDU-21)."""
        return sorted(self.interaction, key=lambda n: (round(self.nodes[n][2], 9), n))

    def deck(self, model: str = "m", title: str = "") -> decks.Deck:
        d = _common(decks.new("HOUSE"), model, title or self.title)
        s = self.site
        d["gravity"], d["gelev"], d["dim"], d["imp"] = s.gravity, s.gelev, self.dim, self.imp
        d["incomp"], d["cmodform"] = self.incomp, s.cmodform
        for nid in sorted(self.nodes):
            d.table("nodes").append([nid] + self.nodes[nid] + self.fixity[nid])
        for n in self.interaction_order():
            d.table("interaction").append([n])
        for g in self.groups:
            d.table("groups").append(list(g))
        for e in self.elements:
            d.table("elements").append(list(e))
        for m in self.materials:
            d.table("materials").append(list(m))
        rows = s.layer_rows()
        for r in rows:
            d.table("layers").append(list(r))
            d.table("sitelayers").append(list(r))
        for r in self.sections:
            d.table("beamprops").append(list(r))
        for r in self.springs:
            d.table("springprops").append(list(r))
        for r in self.matrices:
            d.table("matrices").append(list(r))
        for r in self.masses:
            d.table("masses").append(list(r))
        for r in sorted(self.symm):
            d.table("symm").append(list(r))
        return d


@dataclass
class Model:
    """A canonical HOUSE model: the builder, named nodes and the recommended POINT parameters."""

    house: HouseBuilder
    tags: Dict[str, object] = field(default_factory=dict)
    rad: float = 1.0          # POINT central-zone radius R0 (R1 3.4)
    layer: int = 0            # POINT <layer>: number of embedded TOPL layers
    description: str = ""
    dim: int = 2              # HOUSE <dim> / POINT <dim>: 2 3D (POINT3), 1 2D (POINT2)

    def __getitem__(self, name: str):
        return self.tags[name]

    def deck(self, model: str = "m") -> decks.Deck:
        return self.house.deck(model, self.description)


# ======================================================================================
# Canonical models
# ======================================================================================
def _rigid_link_props(hb: HouseBuilder, length: float, E: Optional[float] = None) -> Tuple[int, int]:
    """Massless stiff beam material and section: E = RIGID_FACTOR x max soil G, square section of
    side ``length`` (shear areas and inertias of that square)."""
    E = E if E is not None else RIGID_FACTOR * hb.site.max_G()
    mat = hb.material(1, E, 0.25, 0.0, 0.0, 0.0)
    a = float(length)
    sec = hb.section(A=a * a, As2=a * a * 5.0 / 6.0, As3=a * a * 5.0 / 6.0, J=0.1406 * a ** 4,
                     I2=a ** 4 / 12.0, I3=a ** 4 / 12.0)
    return mat, sec


def _orientation_node(hb: HouseBuilder, i: int, j: int) -> int:
    """A K node not collinear with I-J: above I for non-vertical members, beside I for vertical ones."""
    xi, xj = hb.xyz(i), hb.xyz(j)
    d = xj - xi
    L = float(np.linalg.norm(d))
    vertical = abs(d[2]) > 0.99 * L
    off = np.array([L, 0.0, 0.0]) if vertical else np.array([0.0, 0.0, L])
    return hb.orientation_node(*(xi + off))


def rigid_link(hb: HouseBuilder, i: int, j: int, size: float, E: Optional[float] = None,
               group: Optional[int] = None) -> int:
    """Stiff massless BEAMS between nodes i and j (a "rigid link")."""
    mat, sec = _rigid_link_props(hb, size, E)
    return hb.beam(i, j, _orientation_node(hb, i, j), mat, sec, group=group or hb.group(2, "rigid links"))


def _mat_points(half_width: float, ndiv: int, shape: str) -> Tuple[np.ndarray, float, float]:
    """Plan points of a surface mat: square grid (ndiv x ndiv cells) or a circular mesh of ``ndiv``
    rings with 6k nodes on ring k.  Returns (xy, element size h, R0)."""
    if shape == "square":
        h = 2.0 * half_width / ndiv
        g = np.linspace(-half_width, half_width, ndiv + 1)
        X, Y = np.meshgrid(g, g, indexing="xy")
        return np.column_stack([X.ravel(), Y.ravel()]), h, 0.90 * h
    if shape == "disk":
        h = half_width / ndiv
        pts = [(0.0, 0.0)]
        for k in range(1, ndiv + 1):
            r = k * h
            for m in range(6 * k):
                t = 2.0 * math.pi * m / (6 * k)
                pts.append((r * math.cos(t), r * math.sin(t)))
        return np.asarray(pts), h, 0.85 * h
    raise ValueError(f"mat shape must be 'square' or 'disk', got {shape!r}")


def surface_rigid_mat(site: Site, half_width: float = 5.0, ndiv: int = 4, shape: str = "square",
                      rigid: str = "beams", centre_height: float = 0.0, E_rigid: Optional[float] = None,
                      mass: float = 0.0, rotary: Sequence[float] = (0.0, 0.0, 0.0), thick: Optional[float] = None,
                      title: str = "") -> Model:
    """Rigid surface foundation: interaction nodes on the ground surface, made rigid by stiff massless
    BEAMS from a centre node C (``rigid='beams'``) or by stiff massless SHELL elements on the square
    grid (``rigid='shell'``, drilling rotations fixed).

    ``centre_height`` > 0 puts C above grade (a structural, non-interaction node); 0 uses the mat
    centre node (needs an even ``ndiv`` for a square).  ``mass`` / ``rotary`` = lumped foundation
    mass and rotary inertias (mass units) at C.  Tags: 'centre', 'mat' (interaction nodes), 'h'."""
    hb = HouseBuilder(site, title=title or f"surface rigid {shape} mat, half width {half_width:g}")
    z0 = site.gelev
    xy, h, R0 = _mat_points(half_width, ndiv, shape)
    mat_nodes = [hb.node(x, y, z0) for x, y in xy]
    hb.set_interaction(mat_nodes)
    if rigid == "beams":
        if centre_height > 0:
            c = hb.node(0.0, 0.0, z0 + centre_height)
        else:
            c = hb.find(0.0, 0.0, z0)
            if c is None:
                raise ValueError("the mat has no centre node: use an even ndiv or centre_height > 0")
        for n in mat_nodes:
            if n != c:
                rigid_link(hb, c, n, h, E_rigid)
    elif rigid == "shell":
        if shape != "square":
            raise ValueError("rigid='shell' needs a square mat")
        if centre_height > 0:
            raise ValueError("rigid='shell' uses the mat centre node (centre_height must be 0)")
        E = E_rigid if E_rigid is not None else RIGID_FACTOR * site.max_G()
        mat = hb.material(1, E, 0.2, 0.0, 0.0, 0.0)
        t = thick if thick is not None else h
        g = hb.group(3, "rigid mat")
        idx = {(i, j): mat_nodes[j * (ndiv + 1) + i] for i in range(ndiv + 1) for j in range(ndiv + 1)}
        for j in range(ndiv):
            for i in range(ndiv):
                hb.shell([idx[(i, j)], idx[(i + 1, j)], idx[(i + 1, j + 1)], idx[(i, j + 1)]], mat, t, group=g)
        for n in mat_nodes:
            hb.fix(n, [6])                       # flat mat: drilling rotation has no stiffness (FIXSHLROT)
        c = hb.find(0.0, 0.0, z0)
        if c is None:
            raise ValueError("the mat has no centre node: use an even ndiv")
    else:
        raise ValueError(f"rigid must be 'beams' or 'shell', got {rigid!r}")
    if mass > 0 or any(rotary):
        hb.mass(c, mass, mass, mass, *[float(v) for v in rotary])
    return Model(hb, {"centre": c, "mat": mat_nodes, "h": h}, rad=R0, layer=0,
                 description=hb.title)


@dataclass
class Stick:
    """Lumped-mass BEAMS stick: node elevations above its base, masses (mass units) and section."""

    heights: Sequence[float]
    masses: Sequence[float]
    E: float
    A: float
    I: float
    nu: float = 0.2
    beta: float = 0.05
    J: Optional[float] = None
    As: float = 0.0
    rotary: Sequence[float] = ()


def add_stick(hb: HouseBuilder, base: int, stick: Stick) -> List[int]:
    """Vertical BEAMS stick from node ``base`` (doubly symmetric section, massless beams, lumped
    masses at the stick nodes); returns the stick node ids (base excluded)."""
    x0 = hb.xyz(base)
    mat = hb.material(1, stick.E, stick.nu, 0.0, stick.beta, stick.beta)
    J = stick.J if stick.J is not None else 2.0 * stick.I
    sec = hb.section(stick.A, stick.As, stick.As, J, stick.I, stick.I)
    g = hb.group(2, "stick")
    nodes, prev = [], base
    for k, z in enumerate(stick.heights):
        n = hb.node(x0[0], x0[1], x0[2] + float(z))
        hb.beam(prev, n, _orientation_node(hb, prev, n), mat, sec, group=g)
        m = float(stick.masses[k]) if k < len(stick.masses) else 0.0
        r = float(stick.rotary[k]) if k < len(stick.rotary) else 0.0
        if m > 0 or r > 0:
            hb.mass(n, m, m, m, r, r, 0.0)
        nodes.append(n)
        prev = n
    return nodes


def stick_on_mat(site: Site, stick: Stick, half_width: float = 5.0, ndiv: int = 4, shape: str = "square",
                 rigid: str = "beams", mat_mass: float = 0.0, E_rigid: Optional[float] = None, title: str = "") -> Model:
    """A lumped-mass stick on a rigid surface mat (centre node at grade).  Tags as
    :func:`surface_rigid_mat` plus 'stick' (stick nodes bottom-up) and 'top'."""
    m = surface_rigid_mat(site, half_width, ndiv, shape, rigid, 0.0, E_rigid, mat_mass,
                          title=title or f"stick on a rigid {shape} mat")
    nodes = add_stick(m.house, m["centre"], stick)
    m.tags.update(stick=nodes, top=nodes[-1] if nodes else m["centre"])
    return m


def embedded_box(site: Site, half_width: float, n_emb: int, ndiv: int = 4, method: str = "FV",
                 structure: str = "soil", skip: int = 2, E: Optional[float] = None, nu: float = 0.2,
                 rho: float = 0.0, beta: float = 0.02, thick: Optional[float] = None,
                 stick: Optional[Stick] = None, title: str = "") -> Model:
    """Square box excavation (plan 2a x 2a, ``ndiv`` x ``ndiv`` elements) through the top ``n_emb``
    TOPL layers (one element per layer -- subdivide the site layers for a finer vertical mesh).

    ``method`` selects the interaction set (R1 4.5): 'FV' every excavation node; 'FSIN' (FI-FSIN)
    the lateral and bottom faces; 'EVBN' (FI-EVBN) all boundary nodes incl. the ground-surface face;
    'FFV' EVBN plus every ``skip``-th interior horizontal plane (counted from the surface).

    ``structure``: 'none' (an empty pit), 'soil' (structure = SOLIDs with type-3 materials equal to
    the layers: K*_s = K*_e, the zero-SSI identity of VP-16; FV only), 'solid' (an embedded SOLID
    block with E, nu, rho, beta; FV only, because every excavation node is a structure node --
    EXCSTRCHK) or 'shell' (basement: SHELL walls, base slab and roof slab of thickness ``thick`` on
    the excavation boundary; drilling rotations of flat-face nodes fixed).  ``stick``: optional
    lumped-mass stick on the top centre node, tied to its four neighbours on the top face by
    rigid links (a rigid "spider", so the stick is clamped to the foundation).  Tags: 'top_centre', 'bottom_centre', 'excavation' (all box nodes), 'levels' (node ids
    per level, top first), 'stick', 'top'."""
    method = method.upper()
    if method not in ("FV", "FSIN", "EVBN", "FFV"):
        raise ValueError(f"method must be FV, FSIN, EVBN or FFV, got {method!r}")
    if not 1 <= n_emb <= len(site.layers):
        raise ValueError(f"n_emb = {n_emb} must lie in 1..{len(site.layers)}")
    if structure in ("soil", "solid") and method != "FV":
        raise ValueError(f"structure={structure!r} uses every excavation node: only method 'FV' is valid "
                         "(EXCSTRCHK: interior excavation nodes shared with the structure must interact)")
    imp = {"FV": 0, "FFV": 1, "FSIN": 2, "EVBN": 2}[method]
    hb = HouseBuilder(site, imp=imp, incomp=1,
                      title=title or f"embedded box {2 * half_width:g} x {2 * half_width:g}, {n_emb} layers, {method}, "
                                     f"structure {structure}")
    z = site.elevations()[:n_emb + 1]                       # top first
    h = 2.0 * half_width / ndiv
    g = np.linspace(-half_width, half_width, ndiv + 1)
    ids: Dict[Tuple[int, int, int], int] = {}
    for lev in range(n_emb, -1, -1):                      # bottom level first: bottom-up numbering
        for j in range(ndiv + 1):
            for i in range(ndiv + 1):
                ids[(i, j, lev)] = hb.node(g[i], g[j], z[lev])
    gexc = hb.group(1, "excavated soil")
    gstr = hb.group(1, "structure") if structure in ("soil", "solid") else None
    smat = None
    if structure == "solid":
        Es = E if E is not None else RIGID_FACTOR * site.max_G() / 100.0
        smat = hb.elastic(Es, nu, rho, beta)
    for lev in range(n_emb):                              # element between level lev (top) and lev+1
        for j in range(ndiv):
            for i in range(ndiv):
                bot = [ids[(i, j, lev + 1)], ids[(i + 1, j, lev + 1)], ids[(i + 1, j + 1, lev + 1)], ids[(i, j + 1, lev + 1)]]
                top = [ids[(i, j, lev)], ids[(i + 1, j, lev)], ids[(i + 1, j + 1, lev)], ids[(i, j + 1, lev)]]
                hb.solid(bot + top, mat=lev + 1, etype=2, group=gexc)          # MSET = TOPL layer lev+1
                if structure == "soil":
                    hb.solid(bot + top, mat=hb.soil_material(lev + 1), etype=1, group=gstr)
                elif structure == "solid":
                    hb.solid(bot + top, mat=smat, etype=1, group=gstr)
    # ---- interaction set (R1 4.5) --------------------------------------------------------------------
    def lateral(i, j):
        return i in (0, ndiv) or j in (0, ndiv)

    inter = []
    for (i, j, lev), n in ids.items():
        on_bottom = lev == n_emb
        on_top = lev == 0
        if method == "FV":
            keep = True
        elif method == "FSIN":
            keep = lateral(i, j) or on_bottom
        elif method == "EVBN":
            keep = lateral(i, j) or on_bottom or on_top
        else:
            keep = lateral(i, j) or on_bottom or on_top or (lev % max(1, int(skip)) == 0)
        if keep:
            inter.append(n)
    hb.set_interaction(inter)
    # ---- shell basement ---------------------------------------------------------------------------
    if structure == "shell":
        Es = E if E is not None else RIGID_FACTOR * site.max_G() / 100.0
        smat = hb.elastic(Es, nu, rho, beta)
        t = thick if thick is not None else 0.1 * h
        gs = hb.group(3, "basement")
        faces = 0
        for lev in (0, n_emb):                              # roof and base slabs
            for j in range(ndiv):
                for i in range(ndiv):
                    hb.shell([ids[(i, j, lev)], ids[(i + 1, j, lev)], ids[(i + 1, j + 1, lev)], ids[(i, j + 1, lev)]],
                             smat, t, group=gs)
                    faces += 1
        for lev in range(n_emb):                            # walls x = +-a and y = +-a
            for k in range(ndiv):
                for i in (0, ndiv):
                    hb.shell([ids[(i, k, lev + 1)], ids[(i, k + 1, lev + 1)], ids[(i, k + 1, lev)], ids[(i, k, lev)]],
                             smat, t, group=gs)
                for j in (0, ndiv):
                    hb.shell([ids[(k, j, lev + 1)], ids[(k + 1, j, lev + 1)], ids[(k + 1, j, lev)], ids[(k, j, lev)]],
                             smat, t, group=gs)
        # drilling rotations: nodes on exactly one flat face (not on an edge) need the rotation about
        # the face normal fixed (FIXSHLROT, EDU-06)
        for (i, j, lev), n in ids.items():
            planes = []
            if i in (0, ndiv):
                planes.append(4)
            if j in (0, ndiv):
                planes.append(5)
            if lev in (0, n_emb):
                planes.append(6)
            if len(planes) == 1:
                hb.fix(n, planes)
    elif structure not in ("none", "soil", "solid"):
        raise ValueError(f"structure must be 'none', 'soil', 'solid' or 'shell', got {structure!r}")
    tags: Dict[str, object] = {"excavation": sorted(ids.values()),
                               "levels": [[ids[(i, j, lev)] for j in range(ndiv + 1) for i in range(ndiv + 1)]
                                          for lev in range(n_emb + 1)], "h": h}
    if ndiv % 2 == 0:
        tags["top_centre"] = ids[(ndiv // 2, ndiv // 2, 0)]
        tags["bottom_centre"] = ids[(ndiv // 2, ndiv // 2, n_emb)]
    if stick is not None:
        if "top_centre" not in tags:
            raise ValueError("a stick needs a top centre node: use an even ndiv")
        base = tags["top_centre"]
        c = ndiv // 2
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):        # rigid spider: fixed-base stick
            rigid_link(hb, base, ids[(c + di, c + dj, 0)], h)
        nodes = add_stick(hb, base, stick)
        tags.update(stick=nodes, top=nodes[-1] if nodes else base)
    return Model(hb, tags, rad=0.90 * h, layer=n_emb, description=hb.title)


def sdof_on_node(site: Site, k: float, m: float, beta: float, element: str = "spring", direction: int = 1,
                 height: float = 1.0, title: str = "") -> Model:
    """Spring-mass oscillator on one surface interaction node (the "fixed-base" SDOF of VP-01/VP-39
    when the site is very stiff).

    Base node B at (0, 0, gelev) (interaction node, rotations fixed); mass node M at height
    ``height`` above it with mass ``m`` (mass units) and only DOF ``direction`` (1 X, 2 Y, 3 Z) free;
    ``element`` 'spring' (SPRING, K* = k c(beta) with the site CMODFORM), 'general' (GENERAL, 2 nodes,
    global axes, MXR/MXI from c(beta)) or 'general3' (GENERAL, 3 nodes, local axes: axis 1 = B->M
    (vertical), axis 2 toward a K node offset along the global ``direction``, so local DOF 2 carries
    the spring).  Tags: 'base', 'mass'."""
    from ..conventions import cfactor
    hb = HouseBuilder(site, title=title or f"SDOF ({element}) k = {k:g}, m = {m:g}, beta = {beta:g}")
    z0 = site.gelev
    b = hb.add_node(0.0, 0.0, z0, fix=(0, 0, 0, 1, 1, 1))
    fixm = [1] * 6
    fixm[direction - 1] = 0
    mn = hb.add_node(0.0, 0.0, z0 + height, fix=fixm)
    hb.set_interaction([b])
    k6 = [0.0] * 6
    k6[direction - 1] = k
    c = complex(cfactor(beta, site.cmodform))
    if element == "spring":
        hb.spring(b, mn, k6, beta)
    elif element == "general":
        hb.general((b, mn), _spring_pattern(np.asarray(k6) * c.real), _spring_pattern(np.asarray(k6) * c.imag))
    elif element == "general3":
        if direction == 3:
            raise ValueError("'general3' needs a horizontal direction (local axis 1 is vertical)")
        off = np.zeros(3)
        off[direction - 1] = 1.0
        kn = hb.orientation_node(*(np.array([0.0, 0.0, z0]) + off))
        kl = [0.0, k, 0.0, 0.0, 0.0, 0.0]
        hb.general((b, mn, kn), _spring_pattern(np.asarray(kl) * c.real), _spring_pattern(np.asarray(kl) * c.imag))
    else:
        raise ValueError(f"element must be 'spring', 'general' or 'general3', got {element!r}")
    hb.mass(mn, m, m, m)
    return Model(hb, {"base": b, "mass": mn}, rad=1.0, layer=0, description=hb.title)


# ======================================================================================
# 2D (plane-strain) models: HOUSE <dim> = 1, POINT2 (requirements 1.5 PLANE, 4.3, 4.6 item 2)
# ======================================================================================
#: out-of-plane DOFs of a 2D model in the X-Z plane (UY, ROTX, ROTZ): fixed on 6-DOF nodes
OUT_OF_PLANE = (2, 4, 6)


def rigid_strip_2d(site: Site, half_width: float = 1.0, ndiv: int = 8, mass: float = 0.0, rotary: float = 0.0,
                   E_rigid: Optional[float] = None, title: str = "") -> Model:
    """Rigid surface strip of half-width B (2D, X-Z plane, per unit length): ``ndiv + 1`` interaction nodes
    at x = -B..B on the ground surface (``ndiv`` even), made rigid by stiff massless BEAMS from the centre
    node; the out-of-plane DOFs (UY, ROTX, ROTZ) are fixed.  ``mass`` / ``rotary`` = foundation mass and
    rotary inertia about Y at the centre (per unit length).  R0 = h (R1 3.4, 2D rule).  Tags: 'centre',
    'strip' (interaction nodes, x ascending), 'h'."""
    if ndiv < 2 or ndiv % 2:
        raise ValueError("rigid_strip_2d needs an even ndiv >= 2 (centre node)")
    hb = HouseBuilder(site, dim=1, title=title or f"2D rigid surface strip, half width {half_width:g}")
    z0 = site.gelev
    h = 2.0 * half_width / ndiv
    nodes = [hb.node(x, 0.0, z0) for x in np.linspace(-half_width, half_width, ndiv + 1)]
    hb.set_interaction(nodes)
    c = nodes[ndiv // 2]
    for n in nodes:
        hb.fix(n, OUT_OF_PLANE)
    for n in nodes:
        if n != c:                       # section of side B: links much stiffer than the soil for any ndiv
            rigid_link(hb, c, n, half_width, E_rigid)
    if mass > 0 or rotary > 0:
        hb.mass(c, mass, mass, mass, 0.0, rotary, 0.0)
    return Model(hb, {"centre": c, "strip": nodes, "h": h}, rad=h, layer=0, description=hb.title, dim=1)


def embedded_plane_2d(site: Site, half_width: float, n_emb: int, ndiv: int = 4, structure: str = "soil",
                      E: Optional[float] = None, nu: float = 0.2, rho: float = 0.0, beta: float = 0.02,
                      stick: Optional[Stick] = None, x0: float = 0.0, title: str = "") -> Model:
    """2D excavation of PLANE elements (plane strain, X-Z plane) of width 2a (``ndiv`` columns) through
    the top ``n_emb`` TOPL layers (one element per layer), flexible volume (every excavation node
    interacts), centred at x = ``x0``.  ``structure``: 'none' (empty pit), 'soil' (PLANE structure with
    type-3 materials equal to the layers: K*_s = K*_e, the 2D zero-SSI identity of VP-T3) or 'plane' (an
    embedded PLANE block with E, nu, rho, beta).  ``stick``: optional lumped-mass BEAMS stick on the top
    centre node (out-of-plane DOFs fixed).  Tags: 'top_centre', 'excavation', 'levels' (top first), 'h',
    'stick', 'top'."""
    if not 1 <= n_emb <= len(site.layers):
        raise ValueError(f"n_emb = {n_emb} must lie in 1..{len(site.layers)}")
    hb = HouseBuilder(site, dim=1, incomp=1,
                      title=title or f"2D embedded excavation {2 * half_width:g} wide, {n_emb} layers, structure {structure}")
    z = site.elevations()[:n_emb + 1]                       # top first
    h = 2.0 * half_width / ndiv
    xs = x0 + np.linspace(-half_width, half_width, ndiv + 1)
    ids: Dict[Tuple[int, int], int] = {}
    for lev in range(n_emb, -1, -1):                      # bottom level first: bottom-up numbering
        for i in range(ndiv + 1):
            ids[(i, lev)] = hb.node(xs[i], 0.0, z[lev])
    gexc = hb.group(4, "excavated soil")
    gstr = hb.group(4, "structure") if structure in ("soil", "plane") else None
    smat = None
    if structure == "plane":
        Es = E if E is not None else RIGID_FACTOR * site.max_G() / 100.0
        smat = hb.elastic(Es, nu, rho, beta)
    elif structure not in ("none", "soil"):
        raise ValueError(f"structure must be 'none', 'soil' or 'plane', got {structure!r}")
    for lev in range(n_emb):                              # element between level lev (top) and lev + 1
        for i in range(ndiv):
            quad = [ids[(i, lev + 1)], ids[(i + 1, lev + 1)], ids[(i + 1, lev)], ids[(i, lev)]]   # counter-clockwise
            hb.plane(quad, mat=lev + 1, etype=2, group=gexc)
            if structure == "soil":
                hb.plane(quad, mat=hb.soil_material(lev + 1), etype=1, group=gstr)
            elif structure == "plane":
                hb.plane(quad, mat=smat, etype=1, group=gstr)
    hb.set_interaction(list(ids.values()))
    tags: Dict[str, object] = {"excavation": sorted(ids.values()), "h": h,
                               "levels": [[ids[(i, lev)] for i in range(ndiv + 1)] for lev in range(n_emb + 1)]}
    if ndiv % 2 == 0:
        tags["top_centre"] = ids[(ndiv // 2, 0)]
        tags["bottom_centre"] = ids[(ndiv // 2, n_emb)]
    if stick is not None:
        if "top_centre" not in tags:
            raise ValueError("a stick needs a top centre node: use an even ndiv")
        base = tags["top_centre"]
        c = ndiv // 2
        spider = [ids[(c - 1, 0)], ids[(c + 1, 0)]]
        for n in spider:                                  # rigid links to the neighbours: the stick is clamped
            rigid_link(hb, base, n, h)                    # (a PLANE node carries no rotation, HINGED)
        nodes = add_stick(hb, base, stick)
        for n in [base] + spider + nodes:
            hb.fix(n, OUT_OF_PLANE)
        tags.update(stick=nodes, top=nodes[-1] if nodes else base)
    return Model(hb, tags, rad=h, layer=n_emb, description=hb.title, dim=1)


# ======================================================================================
# Runners
# ======================================================================================
def write_deck(workdir: PathLike, model: str, deck: decks.Deck) -> Path:
    """Write ``deck`` as ``<model>.<ext>`` in ``workdir`` (what AFWRITE does)."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    deck["model"] = model
    return decks.write(decks.deck_path(workdir, model, deck.module), deck)


def run(module: str, workdir: PathLike, model: str, check: bool = True) -> int:
    """Run one module like RUN<MODULE>; with ``check`` raise :class:`ChainError` on failure."""
    rc = run_module(module, model, Path(workdir))
    if rc != 0 and check:
        p = Path(workdir) / f"{model}_{module.upper()}.out"
        tail = p.read_text(encoding="utf-8")[-4000:] if p.exists() else ""
        raise ChainError(f"{module} failed (status {rc}):\n{tail}")
    return rc


def listing(workdir: PathLike, model: str, module: str) -> str:
    p = Path(workdir) / f"{model}_{module.upper()}.out"
    return p.read_text(encoding="utf-8") if p.exists() else ""


def copy_file(workdir: PathLike, src: str, dst: str) -> Path:
    """FCOPY,src,dst in the model directory (e.g. FILE1 -> FILE1X, FILE9 -> FILE9001)."""
    return Path(shutil.copyfile(Path(workdir) / src, Path(workdir) / dst))


def run_soil(workdir: PathLike, model: str, site: Site, fs: FrequencySet, layer: int = 0, rad: float = 1.0,
             wave: str = "SV", cl: int = 1, mode2: bool = True, point: bool = True,
             waves: Optional[Sequence[Sequence[float]]] = None, cm: Optional[int] = None,
             dim: int = 2) -> Dict[str, Container]:
    """SITE (Mode 1 + Mode 2 for ``wave``) and POINT (``dim`` 2 POINT3, 1 POINT2); returns {'FILE1',
    'FILE2', 'FILE3'}."""
    write_deck(workdir, model, site_deck(site, fs, wave=wave, cl=cl, mode2=int(mode2), waves=waves, cm=cm,
                                         model=model))
    run("SITE", workdir, model)
    out = {"FILE2": read_container(Path(workdir) / "FILE2", "FILE2")}
    if mode2:
        out["FILE1"] = read_container(Path(workdir) / "FILE1", "FILE1")
    if point:
        write_deck(workdir, model, point_deck(fs, layer, rad, model=model, dim=dim))
        run("POINT", workdir, model)
        out["FILE3"] = read_container(Path(workdir) / "FILE3", "FILE3")
    return out


def run_site_xyz(workdir: PathLike, model: str, site: Site, fs: FrequencySet, cl: int = 1,
                 keep: str = "X") -> List[str]:
    """Three SITE runs (vertical SV x', SH y', P z'; angle 0) copied to FILE1X, FILE1Y, FILE1Z
    (requirements 2.4, D-ANL-06).  FILE1 is left equal to FILE1<keep>.  The .sit deck finally
    written is the SV one (HOUSE reads its layer table)."""
    names = []
    for wave, c in (("SH", "Y"), ("P", "Z"), ("SV", "X")):
        write_deck(workdir, model, site_deck(site, fs, wave=wave, cl=cl, mode1=1, mode2=1, model=model))
        run("SITE", workdir, model)
        copy_file(workdir, "FILE1", f"FILE1{c}")
        names.append(f"FILE1{c}")
    copy_file(workdir, f"FILE1{keep}", "FILE1")
    return sorted(names)


def run_house(workdir: PathLike, model: str, mdl: Union[Model, HouseBuilder]) -> Container:
    """Write the HOUSE deck (the .sit deck must exist) and run HOUSE; returns FILE4."""
    hb = mdl.house if isinstance(mdl, Model) else mdl
    write_deck(workdir, model, hb.deck(model, getattr(mdl, "description", "")))
    run("HOUSE", workdir, model)
    return read_container(Path(workdir) / f"{model}.N4", "FILE4")


def run_force(workdir: PathLike, model: str, loads: Sequence[Sequence[float]], fs: FrequencySet,
              copy_to: Optional[str] = None, gravity: float = GRAVITY) -> Container:
    """FORCE run (FILE9), optionally copied to ``copy_to`` (e.g. 'FILE9002')."""
    write_deck(workdir, model, force_deck(loads, fs, gravity=gravity, model=model))
    run("FORCE", workdir, model)
    if copy_to:
        copy_file(workdir, "FILE9", copy_to)
    return read_container(Path(workdir) / "FILE9", "FILE9")


def run_analys(workdir: PathLike, model: str, fs: FrequencySet, check: bool = True, gravity: float = GRAVITY,
               **params) -> int:
    """ANALYS run with deck parameters ``params``; returns the status (raises on failure when ``check``)."""
    write_deck(workdir, model, analys_deck(fs, gravity=gravity, model=model, **params))
    return run("ANALYS", workdir, model, check=check)


def run_motion(workdir: PathLike, model: str, fs: FrequencySet, thfile: str, nout: Sequence[Sequence[int]],
               check: bool = True, **params) -> int:
    write_deck(workdir, model, motion_deck(fs, thfile, nout, model=model, **params))
    return run("MOTION", workdir, model, check=check)


def read_file8(workdir: PathLike, name: str = "FILE8") -> Container:
    return read_container(Path(workdir) / name, "FILE8")


def tf(file8: Container, node: int, dof: int) -> np.ndarray:
    """Transfer function (nF,) of (node, dof) in a FILE8 (zeros when it is not an equation: fixed DOF)."""
    hit = np.flatnonzero((np.asarray(file8["eq_node"]) == int(node)) & (np.asarray(file8["eq_dof"]) == int(dof)))
    if hit.size == 0:
        return np.zeros(len(file8["fnum"]), complex)
    return np.asarray(file8["H"])[:, int(hit[0])]
