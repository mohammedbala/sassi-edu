"""HOUSE module: frequency-independent model matrices of the structure and of the excavated soil
-> FILE4 (``<model>.N4``), COOSK, COOSM, DOFSMAP, FILE90, FILE91 (requirements 2.1, 4.1, 4.4).

What HOUSE does (for the structural engineer)
---------------------------------------------
The flexible-volume method (Lysmer et al. 1981; R1 section 4.4) solves, per frequency,

    [ (K*_s - w^2 M_s) - (K*_e - w^2 M_e) + X_ff ] U = X_ff U'_f

where ``K*_s, M_s`` belong to the structure (including any near-field soil modelled as structure),
``K*_e, M_e`` to the *excavated* soil -- the soil volume the basement displaces, with the free-field
layer properties -- and ``X_ff`` is the soil impedance at the interaction nodes (ANALYS).  None of
the four matrices depends on the frequency (hysteretic damping is a complex modulus,
requirements 4.0.2), so HOUSE assembles them once, on one global equation numbering, and ANALYS
combines them at every frequency.

Steps (requirements 4.4 items 1-4):

1. read the HOUSE deck; build the materials: structure from the M table (``material_from_M``),
   excavated SOLID/PLANE elements from the L layer given by their MSET index
   (``material_from_layer``, D-ELM-12), with the deck gravity (rho = weight / g) and CMODFORM;
2. classify SOLID/PLANE elements (ETYPE 1 structure, 2 excavated soil, 0 resolved by D-HOU-01);
3. build the DOF map (D-ELM-09: only DOFs some element defines are kept; interaction nodes always
   keep their translations) and assemble ``K*_s, M_s, K*_e, M_e`` with :mod:`sassi.elements`;
4. place the interaction nodes on the SITE user interfaces (depth = gelev - z, D-GEN-06 tolerance,
   EDU-01) and check the interaction-set rules (above grade, bottom-up order, G-14, EXCSTRCHK,
   EDU-22, Error 124).  The interaction DOFs are the translations UX, UY, UZ in 3D and the
   in-plane UX, UZ in 2D (FILE4 ``int_eq`` has -1 in the UY column of a 2D model, meta
   ``int_dofs``);
5. write FILE4 (DOF map, interaction nodes, topology and the STRESS recovery operators), COOSK,
   COOSM, DOFSMAP, FILE90 (restart hashes, D-W1-13), FILE91 (metadata) and a listing with the
   excavated-element layer assignment table that the manual asks the user to review (05b R7),
   the EDU-08/EDU-10 checks and the excavation numbering/grouping rules R5/R6 (warnings,
   D-HOU-04; the findings are stored in FILE4 meta ``excav_layering`` for STRESS).

2D models (<dim> = 1, P1; :mod:`sassi.core.ssi2d`): PLANE elements (plane strain, unit thickness) in the
X-Z plane, BEAMS / SPRING / GENERAL in the plane; the interaction DOFs are UX, UZ; the 6-DOF nodes keep
their out-of-plane DOFs (UY, ROTX, ROTZ: structural only, fix them with D unless something restrains
them).  SOLID / SHELL / TSHELL elements are errors (EDU-26, D-PNT-01).

Symmetry planes (SYMM, P1; D-ANL-12, :mod:`sassi.core.symmetry`): the planes of a half or quarter model
(3D: parallel to XZ or YZ, at most two and orthogonal; 2D: one line parallel to Z) are checked (manual
Errors 2-5, EDU-26 geometry; the model must lie on one side of every plane) and their boundary conditions
are added to the D fixities of the plane nodes -- symmetric loading (type 0): the normal translation and
the in-plane rotations, antisymmetric (type 1): the in-plane translations and the normal rotation.
FILE4 carries the planes (``x_symm``) and the constrained interaction translations (``x_int_symm``) for
ANALYS, which builds the soil impedance of the reduced set from image interaction nodes; FILE90
``int_hash`` includes the planes.  Incoherent motion, wave passage and multiple excitation are errors
with SYMM planes and in 2D (full 3D models only, manual 2.7 and 6.5.4, G-17).

Elements: SOLID, BEAMS, SHELL, PLANE, TSHELL (thick Mindlin-Reissner shell, EINT 0 reduced / 1 selective,
automatic drilling stiffness, D-ELM-08, :mod:`sassi.elements.tshell`), SPRING, GENERAL.  FILE4 carries the
TSHELL thickness, E and nu per element (``x_rec_TSHELL_*``) for the THSHLSTR face stresses of STRESS.

Node-numbering optimizer ("Optimize Model", HOUSEX <optimize> = 1; requirements 4.4 item 5, D-HOU-03,
D-FIL-09; :mod:`sassi.core.renumber`): the model is read, assembled and checked in its own numbering
(every listed error and warning quotes the node numbers of ``<model>.hou``); only a model without errors
is renumbered 1..N by reverse Cuthill-McKee, started from the interaction nodes, which keep their
bottom-up order, and only if the bandwidth and profile of the assembled ``|K_s| + |K_e|`` do not grow
(spec 05b test 8).  The renumbered deck is assembled again and FILE4, COOSK, COOSM are written from it --
so they and the FILE8 of ANALYS use the NEW node numbers -- together with ``<model>.hounew`` and the
``old new`` pairs ``<model>.map`` (a failed run writes neither).  As the manual warns, MOTION/RELDISP
output nodes must then be given with the new numbers; FORCE loads (FILE9) are in the model numbering and
ANALYS translates them to the new numbers (it stops when the pairs are missing).  FILE4 carries
``x_node_old_id`` and the meta ``x_optimized`` / ``x_node_map``.  The sparse solver of ANALYS orders the equations itself, so
the results do not change (VP-O1).

Incoherent motion, wave passage and multiple excitation (HOUSE <coh>, <wpass>, <me>; INCOH, WPASS, ME,
AMP; requirements 4.4 item 6, D-INC-*; :mod:`sassi.core.coherency`, :mod:`sassi.core.incoherency`): for
every SSI frequency and direction X, Y, Z HOUSE builds the coherency matrix of the interaction nodes
(horizontal projections; unlagged coherency model WPASS <cohf>), factorises it (lambda descending, sign
convention D-INC-03, trace check Eq. 6.1, the ``I N C O`` table when INCOH <ipr> = 1) and writes the
incoherency factors ``s`` of the interaction nodes -- deterministic AS sum or single mode: one FILE77;
stochastic simulation (HSeed/VSeed, RandPhz): FILE77001 .. FILE77Ns with HOUSEX <nsim> -- times the wave
passage delay ``exp(-i w tau)`` along Line D and the spectral amplification ratios of the ME zones.
ANALYS applies them to the free-field load (FFL) or motion (FFM).  A HOUSE run that writes FILE77
renames the incoherency files of the other kind left by an earlier run ``.prev``.  Option AA (ANSYS
model input, P2) stops the run with a clear message.

Non-linear soil SSI (HOUSEX <nlssi> = 1, requirements 4.4 item 7, :mod:`sassi.core.nlsoil`): the
SOLID (3D) / PLANE (2D) structure elements of the ``.pin`` groups get one internal material per
element -- iteration 0 ``G = GFAC G_layer``, ``beta = DFAC beta_layer`` of the free-field TOPL layer
at the element's mid-depth; later runs (``.liq`` = 1) the strain-compatible
``G = G_max (G/G_max)(gamma_eff)``, ``beta = D(gamma_eff)`` from the FILE74 strains and the FILE73
curves, with G_max the low-strain modulus of the element's M-table material -- and HOUSE writes
FILE78 (properties used), ``.liq`` = 1 and the convergence row of the previous iteration
(NLSOIL_CONVERGENCE.TXT, D-NLS-02).  A linear run renames a FILE78 left by an earlier non-linear run
FILE78.prev (it does not describe the new FILE4).

    python -m sassi.modules.house < HOUSE.inp      # batch use
"""
from __future__ import annotations

import time
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp

from .. import PRODUCT, __version__
from ..conventions import ELEMENT_COMPONENTS, ELEMENT_TYPE_NAMES, frequency_step, unit_system
from ..core import coherency as COH
from ..core import house_lib as hl
from ..core import incoherency as INC
from ..core import nlsoil as NL
from ..core import renumber as RN
from ..core import symmetry as SYM
from ..elements import (ELEMENTS, ElemRecord, ElementError, MaterialProps, assemble, build_dofmap,
                        material_from_layer, material_from_M, nodal_masses_from_table, unstiffened_rotations)
from ..elements.assemble import DofMap, element_dof_nodes, is_excavated_soil
from ..io import decks
from ..io.container import Container
from ..io.deckfmt import Table
from ..io.files import validate, write_container
from .base import ModuleContext, ModuleError, batch_main

NAME = "HOUSE"

#: the excavated-element layer assignment table lists at most this many rows
MAX_LAYER_ROWS = 5000
#: warnings repeated per node/element are listed individually up to this count
MAX_REPEAT = 25


def file4_name(model: str) -> str:
    """FILE4 is written as ``<model>.N4`` (requirements 1.7, 2.1)."""
    return f"{model}.N4"


def deck_model_hash(path) -> str:
    """The interpreter model hash AFWRITE writes as the deck comment ``* model_hash = <hash>``
    ('' when absent, e.g. for hand-written decks)."""
    try:
        with open(path, encoding="utf-8") as fh:
            for _, line in zip(range(20), fh):
                t = line.strip()
                if t.startswith("*") and "model_hash" in t and "=" in t:
                    return t.split("=", 1)[1].strip()
    except OSError:
        pass
    return ""


# ======================================================================================
# Model data read from the deck
# ======================================================================================
@dataclass
class HouseModel:
    """The analysis model of one HOUSE run (deck data resolved into element records).

    ``records[k]`` is the :class:`sassi.elements.ElemRecord` of deck element row ``rows[k]``;
    ``etype[k]`` its resolved ETYPE (D-HOU-01).  ``masses`` are in mass units."""

    gravity: float
    gelev: float
    dim: int
    form: int
    node_ids: np.ndarray
    xyz: np.ndarray
    fix: np.ndarray
    node_xyz: Dict[int, np.ndarray]
    groups: Dict[int, Tuple[int, str]]
    records: List[ElemRecord] = field(default_factory=list)
    rows: List[dict] = field(default_factory=list)
    etype: List[int] = field(default_factory=list)
    masses: Dict[int, np.ndarray] = field(default_factory=dict)
    mass_rows: List[dict] = field(default_factory=list)
    interaction: List[int] = field(default_factory=list)
    site: List[dict] = field(default_factory=list)
    materials: Dict[int, MaterialProps] = field(default_factory=dict)
    layers: Dict[int, MaterialProps] = field(default_factory=dict)
    tol: float = 1e-9
    imp: int = 0                       # HOUSE <imp>: 0 FV, 1 FFV, 2 FI (EDU-22)
    rigid_base: bool = False           # SITE <nl> = 0: the last interface is a rigid base
    n_etype0: int = 0                  # SOLID/PLANE elements with ETYPE 0 resolved here (D-HOU-01)
    n_etype0_exc: int = 0              # ... of which excavated
    excavation: Dict[str, set] = field(default_factory=dict)
    layering: Dict[str, List[str]] = field(default_factory=dict)   # R5/R6 findings (D-HOU-04)
    symm: List[SYM.SymmetryPlane] = field(default_factory=list)    # SYMM planes (D-ANL-12)
    symm_fix: Optional[np.ndarray] = None   # (nN, 6) DOFs fixed by the symmetry conditions
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    info: List[str] = field(default_factory=list)

    def z(self, n: int) -> float:
        return float(self.node_xyz[int(n)][2])

    @property
    def depths(self) -> np.ndarray:
        """Depths of the user interfaces 1..nI (the last ``site`` row is the half-space)."""
        if not self.site:
            return np.zeros(0)
        return hl.interface_depths([r["thick"] for r in self.site[:-1]])


def _err(hm: HouseModel, text: str) -> None:
    if text not in hm.errors:
        hm.errors.append(text)


def _warn(hm: HouseModel, text: str) -> None:
    if text not in hm.warnings:
        hm.warnings.append(text)


def _site_rows(d: decks.Deck, sit: Optional[decks.Deck], hm_warn: List[str]) -> List[dict]:
    """SITE layer table: the HOUSE deck ``sitelayers`` (TOPL layers top first + half-space row);
    cross-checked against the ``.sit`` deck (D-HOU-02), which is used when the HOUSE deck has none."""
    rows = d.rows("sitelayers")
    if sit is None:
        return rows
    srows = sit.rows("layers") + sit.rows("halfspace")
    if not rows:
        if srows:
            hm_warn.append("the HOUSE deck has no 'sitelayers' table: the layer table of the .sit deck is used")
        return [dict(r, thick=(r["thick"] if k < len(srows) - 1 else 0.0)) for k, r in enumerate(srows)]
    cols = ("no", "thick", "weight", "vp", "vs", "dp", "ds")
    if len(srows) != len(rows):
        hm_warn.append(f"the .sit deck has {len(srows)} layer rows (TOPL + half-space), the HOUSE deck "
                       f"{len(rows)}: the HOUSE deck table is used -- re-run AFWRITE")
    else:
        for k, (a, b) in enumerate(zip(rows, srows)):
            last = k == len(rows) - 1
            diff = [c for c in cols if not (last and c == "thick") and
                    abs(float(a[c]) - float(b[c])) > 1e-9 * max(1.0, abs(float(b[c])))]
            if diff:
                hm_warn.append(f"site layer row {k + 1} differs between the HOUSE and .sit decks ({', '.join(diff)}): "
                               "the HOUSE deck table is used -- re-run AFWRITE")
                break
    if abs(float(sit["gravity"]) - float(d["gravity"])) > 1e-9 * abs(float(d["gravity"])):
        hm_warn.append(f"gravity {d['gravity']} of the HOUSE deck differs from the .sit deck ({sit['gravity']})")
    if int(sit["cmodform"]) != int(d["cmodform"]):
        hm_warn.append("CMODFORM of the HOUSE deck differs from the .sit deck: excavated soil and free field "
                       "would use different complex moduli (D-CNV-04)")
    return rows


def build_house_model(d: decks.Deck, sit: Optional[decks.Deck] = None) -> HouseModel:
    """Resolve the HOUSE deck into element records, materials, nodal masses and the site table.

    Input errors are collected in ``HouseModel.errors`` (the run stops after listing them all)."""
    g = float(d["gravity"])
    form = int(d["cmodform"])
    nodes = d.rows("nodes")
    ids = np.array([int(r["id"]) for r in nodes], dtype=np.int64)
    xyz = np.array([[float(r["x"]), float(r["y"]), float(r["z"])] for r in nodes], dtype=float).reshape(-1, 3)
    fix = np.array([[int(r[c]) for c in ("fx", "fy", "fz", "frx", "fry", "frz")] for r in nodes],
                   dtype=np.int64).reshape(-1, 6)
    hm = HouseModel(gravity=g, gelev=float(d["gelev"]), dim=int(d["dim"]), form=form, node_ids=ids, xyz=xyz,
                    fix=(fix != 0).astype(np.int64), node_xyz={int(n): xyz[k] for k, n in enumerate(ids)},
                    groups={}, imp=int(d["imp"]), rigid_base=(sit is not None and int(sit["nl"]) == 0))
    if g <= 0:
        _err(hm, f"Error 1: acceleration of gravity {g} must be > 0")
        return hm
    if len(set(ids.tolist())) != ids.size:
        _err(hm, "duplicate node ids in the 'nodes' table")
    site_warn: List[str] = []
    hm.site = _site_rows(d, sit, site_warn)
    for w in site_warn:
        _warn(hm, w)
    for r in d.rows("groups"):
        hm.groups[int(r["id"])] = (int(r["type"]), str(r["title"]))

    # ---- property tables -------------------------------------------------------------------
    mtab = {int(r["id"]): r for r in d.rows("materials")}
    ltab = {int(r["no"]): r for r in d.rows("layers")}
    rtab = {int(r["id"]): r for r in d.rows("beamprops")}
    sctab = {int(r["id"]): r for r in d.rows("springprops")}
    try:
        mxtab = hl.general_matrix_tables(d.rows("matrices"))
    except ValueError as exc:
        _err(hm, f"GENERAL matrix table: {exc}")
        mxtab = {}

    def material(mid: int, where: str) -> Optional[MaterialProps]:
        if mid in hm.materials:
            return hm.materials[mid]
        r = mtab.get(mid)
        if r is None:
            _err(hm, f"Error 13: material {mid} of {where} is not defined (M table)")
            return None
        try:
            mat = material_from_M(int(r["type"]), float(r["val1"]), float(r["val2"]), float(r["weight"]),
                                  float(r["pdamp"]), float(r["sdamp"]), g, form)
        except ElementError as exc:
            _err(hm, f"material {mid}: {exc}")
            return None
        hm.materials[mid] = mat
        return mat

    def layer(lid: int, where: str) -> Optional[MaterialProps]:
        if lid in hm.layers:
            return hm.layers[lid]
        r = ltab.get(lid)
        if r is None:
            _err(hm, f"Error 19: soil layer {lid} (MSET of {where}) is not defined (L table)")
            return None
        try:
            mat = material_from_layer(float(r["vp"]), float(r["vs"]), float(r["weight"]), float(r["dp"]),
                                      float(r["ds"]), g, form)
        except ElementError as exc:
            _err(hm, f"soil layer {lid}: {exc}")
            return None
        hm.layers[lid] = mat
        return mat

    # ---- elements ----------------------------------------------------------------------------
    elems = d.rows("elements")
    used = sorted({int(r[f"n{k}"]) for r in elems for k in range(1, 9) if int(r[f"n{k}"]) > 0}
                  | {int(r["id"]) for r in d.rows("interaction")})
    act = np.array([hm.node_xyz[n] for n in used if n in hm.node_xyz]).reshape(-1, 3)
    hm.tol = hl.geometric_tolerance(act if act.size else xyz)
    incomp = int(d["incomp"]) == 0
    shell_damp_warned = set()
    seen = set()
    for r in elems:
        gid, eid = int(r["group"]), int(r["id"])
        where = f"element {eid} of group {gid}"
        if (gid, eid) in seen:
            _err(hm, f"{where} is defined twice")
            continue
        seen.add((gid, eid))
        if gid not in hm.groups:
            _err(hm, f"{where}: group {gid} is not defined in the 'groups' table")
            continue
        code = hm.groups[gid][0]
        if code not in ELEMENTS:
            _err(hm, f"{where}: unknown group type {code}")
            continue
        nodes = tuple(int(r[f"n{k}"]) for k in range(1, 9))
        live = [n for n in nodes if n]
        missing = [n for n in live if n not in hm.node_xyz]
        if missing:
            _err(hm, f"Error 41: {where} references undefined nodes {missing}")
            continue
        zs = [hm.z(n) for n in dict.fromkeys(live)]
        et = hl.resolve_etype(code, int(r["etype"]), zs, hm.gelev, hm.tol)
        if int(r["etype"]) == 0 and code in hl.SOIL_TYPES:
            hm.n_etype0 += 1
            hm.n_etype0_exc += int(et == 2)
        if int(r["etype"]) == 2 and code not in hl.SOIL_TYPES + hl.SHELL_TYPES:
            _warn(hm, f"{where}: ETYPE 2 has no effect on {ELEMENT_TYPE_NAMES[code]} elements (structure)")
        excav = et == 2 and code in hl.SOIL_TYPES
        props: dict = {}
        mat: Optional[MaterialProps] = None
        mid = int(r["mat"])
        pid = int(r["prop"])
        if code in hl.SOIL_TYPES:
            mat = layer(mid, where) if excav else material(mid, where)
            props["incompatible"] = bool(incomp and not excav)        # D-ELM-02 / D-ELM-03
            if code == 1:
                eint = int(r["eint"])
                if eint not in (0, 1, 2):
                    _err(hm, f"{where}: EINT {eint} must be 0, 1 or 2 (D-ELM-01)")
                    continue
                props["eint"] = eint
        elif code == 2:
            mat = material(mid, where)
            s = rtab.get(pid)
            if s is None:
                _err(hm, f"Error 26: real property {pid} of {where} is not defined (R table)")
                continue
            props["section"] = dict(A=float(s["axial"]), As2=float(s["shear2"]), As3=float(s["shear3"]),
                                    J=float(s["tors"]), I2=float(s["flex2"]), I3=float(s["flex3"]))
            props["ki"] = str(r["ki"]).strip()
            props["kj"] = str(r["kj"]).strip()
        elif code in (3, 5):                                          # SHELL, TSHELL
            mat = material(mid, where)
            t = float(r["thick"])
            if t <= 0:
                _err(hm, f"{where}: shell thickness {t} must be > 0 (THICK)")
                continue
            props["thick"] = t
            if code == 5:                                             # D-ELM-08, manual 9.4.9
                eint = int(r["eint"])
                if eint not in (0, 1):
                    _err(hm, f"{where}: TSHELL EINT {eint} must be 0 (reduced) or 1 (selective)")
                    continue
                props["eint"] = eint
            mr = mtab.get(mid)
            if mr is not None and float(mr["pdamp"]) != float(mr["sdamp"]) and mid not in shell_damp_warned:
                shell_damp_warned.add(mid)
                _warn(hm, f"material {mid} used by SHELL/TSHELL elements has pdamp != sdamp: the S-wave damping "
                          f"{float(mr['sdamp'])} is used (spec 08 5.20)")
        elif code == 7:
            s = sctab.get(pid)
            if s is None:
                _err(hm, f"Error 33: spring property {pid} of {where} is not defined (SC table)")
                continue
            props.update(k=tuple(float(s[c]) for c in ("scx", "scy", "scz", "scxx", "scyy", "sczz")),
                         damp=float(s["damp"]), form=form)
        elif code == 9:
            p = mxtab.get(pid)
            if p is None:
                _err(hm, f"Error 83: matrix property {pid} of {where} is not defined (MXR/MXI/MXM)")
                continue
            props.update(KR=p["R"], KI=p["I"], MM=p["M"], upper_rows=True, munits=int(d["gmunits"]), gravity=g)
        if code in (1, 2, 3, 4, 5) and mat is None:
            continue
        rec = ElemRecord(group=gid, id=eid, code=code, nodes=nodes, excavated=excav, mat=mat, props=props)
        try:
            element_dof_nodes(rec)
        except ElementError as exc:
            _err(hm, f"{where}: {exc}")
            continue
        hm.records.append(rec)
        hm.rows.append(r)
        hm.etype.append(et)

    # ---- nodal masses (MT/MR, MUNITS: 1 weight -> divided by g, D-MDL-08) -----------------------
    hm.mass_rows = d.rows("masses")
    bad = sorted({int(r["node"]) for r in hm.mass_rows if int(r["node"]) not in hm.node_xyz})
    if bad:
        _err(hm, f"Error 40: nodal masses on undefined nodes {bad[:20]}")
    rows = [[int(r["node"]), r["mx"], r["my"], r["mz"], r["mxx"], r["myy"], r["mzz"], int(r["units"])]
            for r in hm.mass_rows if int(r["node"]) in hm.node_xyz]
    for r in rows:
        if r[7] not in (0, 1):
            _err(hm, f"nodal mass of node {r[0]}: units flag {r[7]} must be 0 (mass) or 1 (weight)")
        if any(float(v) < 0 for v in r[1:7]):
            _err(hm, f"nodal mass of node {r[0]}: negative value")
    hm.masses = nodal_masses_from_table(rows, g)

    # ---- interaction nodes ------------------------------------------------------------------
    inter = [int(r["id"]) for r in d.rows("interaction")]
    dup = sorted(n for n, c in Counter(inter).items() if c > 1)
    if dup:
        _err(hm, f"interaction nodes listed twice: {dup[:20]}")
    und = [n for n in inter if n not in hm.node_xyz]
    if und:
        _err(hm, f"Error 41: interaction nodes {und[:20]} are not defined")
    hm.interaction = list(dict.fromkeys(n for n in inter if n in hm.node_xyz))
    if hm.interaction and not hm.site:
        _err(hm, "interaction nodes are defined but the SITE layer table is empty (TOPL / SITE <hs>)")
    apply_symmetry(hm, d.rows("symm"))
    return hm


def apply_symmetry(hm: HouseModel, rows: List[dict]) -> None:
    """SYMM planes of a half / quarter model (spec 07 9.2.39, D-ANL-12, :mod:`sassi.core.symmetry`):
    the manual's Errors 2-5 and the geometry rules (EDU-26), the model on one side of every plane, and
    the boundary conditions of the symmetry on the nodes of the planes -- symmetric (type 0): the normal
    translation and the in-plane rotations; antisymmetric (type 1): the in-plane translations and the
    normal rotation -- added to the D fixities (``hm.symm_fix`` records them for the listing and the
    interaction-set rules)."""
    if not rows:
        return
    planes, errs, warns = SYM.parse_planes(rows, hm.node_xyz, hm.dim, hm.tol)
    for e in errs:
        _err(hm, e)
    for w in warns:
        _warn(hm, w)
    if errs or not planes:
        return
    hm.symm = planes
    used = sorted({int(n) for r in hm.records for n in element_dof_nodes(r)[0]} | set(hm.interaction))
    P = np.array([hm.node_xyz[n] for n in used if n in hm.node_xyz], dtype=float).reshape(-1, 3)
    for no, neg, pos in SYM.side_violations(planes, P, hm.tol):
        p = next(q for q in planes if q.no == no)
        _err(hm, f"EDU-26: SYMM {p.describe(hm.dim)}: the model has nodes on both sides ({neg} below, {pos} above "
                 f"{SYM.AXIS_NAMES[p.axis].lower()} = {p.coord:.6g}); a half / quarter model lies on one side of "
                 "every symmetry plane")
    fix = SYM.symmetry_fixity(planes, hm.xyz, hm.tol)
    hm.symm_fix = fix
    hm.fix = ((hm.fix != 0) | fix).astype(np.int64)


# ======================================================================================
# Assembly and the interaction set
# ======================================================================================
@dataclass
class HouseResult:
    """Assembled matrices on the analysis numbering plus the interaction data."""

    dofmap: DofMap
    Ks: sp.csr_matrix
    Ms: sp.csr_matrix
    Ke: sp.csr_matrix
    Me: sp.csr_matrix
    recovery: Dict[str, Dict[str, np.ndarray]]
    messages: List[str]
    mass_total_s: np.ndarray          # (3,) structural mass incl. the mass at fixed DOFs
    mass_total_e: np.ndarray          # (3,) excavated soil mass
    mass_free_s: np.ndarray           # (3,) structural mass on the active DOFs
    ignored_masses: List[Tuple[int, int, float]]
    int_node: np.ndarray
    int_iface: np.ndarray
    int_eq: np.ndarray
    int_misfit: np.ndarray


def interaction_dofs(dim: int) -> Tuple[int, ...]:
    """Interaction translations (requirements 4.1): (UX, UZ) in 2D (<dim> = 1), (UX, UY, UZ) in 3D."""
    return (1, 3) if int(dim) == 1 else (1, 2, 3)


def assemble_house(hm: HouseModel) -> HouseResult:
    """DOF map and matrices (requirements 4.1 DOF management, 4.4 item 2).

    The elements are assembled once on a DOF map *without* the D fixities, from which the analysis
    matrices are extracted (rows/columns of the active DOFs): the result is identical to a direct
    assembly on the fixed DOF map, and the unconstrained matrices give the total masses for the
    listing (including the mass that sits on supports)."""
    tr = interaction_dofs(hm.dim)
    extra = {n: tr for n in hm.interaction}                  # interaction DOFs = translations
    dm = build_dofmap(hm.node_ids, hm.fix, hm.records, extra_dofs=extra)
    free = build_dofmap(hm.node_ids, np.zeros_like(hm.fix), hm.records, extra_dofs=extra)
    am = assemble(hm.node_xyz, hm.records, free, masses=hm.masses)
    mask = dm.eq_table >= 0
    keep = free.eq_table[mask]                                # fixed numbering order = free order
    to_fixed = np.full(free.neq, -1, dtype=np.int64)
    to_fixed[keep] = dm.eq_table[mask]

    def cut(A):
        A = sp.csr_matrix(A)[keep][:, keep]
        A.sort_indices()
        return A.tocsr()

    recovery = {}
    for name, rd in am.recovery.items():
        rd = dict(rd)
        eq = rd["eq"]
        rd["eq"] = np.where(eq >= 0, to_fixed[np.clip(eq, 0, None)], -1).astype(np.int64)
        recovery[name] = rd
    ignored = list(am.ignored_masses)
    for node, m6 in hm.masses.items():
        for k in range(6):
            if m6[k] != 0.0 and free.eq(node, k + 1) >= 0 and dm.eq(node, k + 1) < 0:
                ignored.append((int(node), k + 1, float(m6[k])))
    msgs = [m for m in am.messages if "nodal mass terms" not in m]
    Ks, Ms, Ke, Me = cut(am.Ks), cut(am.Ms), cut(am.Ke), cut(am.Me)
    return HouseResult(dofmap=dm, Ks=Ks, Ms=Ms, Ke=Ke, Me=Me, recovery=recovery, messages=msgs,
                       mass_total_s=hl.translational_mass(am.Ms, free.eq_dof),
                       mass_total_e=hl.translational_mass(am.Me, free.eq_dof),
                       mass_free_s=hl.translational_mass(Ms, dm.eq_dof),
                       ignored_masses=sorted(ignored), int_node=np.zeros(0, np.int64),
                       int_iface=np.zeros(0, np.int64), int_eq=np.zeros((0, 3), np.int64),
                       int_misfit=np.zeros(0))


def interaction_set(hm: HouseModel, res: HouseResult) -> None:
    """Interfaces and equation numbers of the interaction nodes and the rules of requirements 4.4
    item 3 (errors and warnings are added to ``hm``)."""
    inter = np.asarray(hm.interaction, dtype=np.int64)
    res.int_node = inter
    n = inter.size
    # The soil acts on the interaction nodes through their *interaction translations* only: UX, UY,
    # UZ in 3D, the in-plane UX, UZ in 2D (requirements 4.1; plane strain in the X-Z plane).  A 2D
    # node that also carries a 6-DOF element (BEAMS, SHELL, SPRING, GENERAL) has a UY equation, but
    # that equation is structural only and must not enter the interaction set (int_eq = -1).
    tr = interaction_dofs(hm.dim)
    res.int_eq = np.array([[res.dofmap.eq(int(k), d) if d in tr else -1 for d in (1, 2, 3)] for k in inter],
                          dtype=np.int64).reshape(n, 3)
    if n == 0:
        res.int_iface = np.zeros(0, np.int64)
        res.int_misfit = np.zeros(0)
        return
    depths = hm.depths
    z = np.array([hm.z(k) for k in inter])
    itol = hl.interface_tolerance(depths, hm.tol)
    iface, mis, above = hl.match_interfaces(z, hm.gelev, depths, itol)
    res.int_iface, res.int_misfit = iface, mis
    nI = depths.size
    for k in range(n):
        node = int(inter[k])
        if above[k]:
            _err(hm, f"EDU-01: interaction node {node} at z = {z[k]:.6g} is above the ground elevation "
                     f"{hm.gelev:.6g} (G-12)")
        elif mis[k] > itol:
            _err(hm, f"EDU-01: interaction node {node} at depth {hm.gelev - z[k]:.6g} is not on a soil-layer "
                     f"interface (nearest: interface {iface[k]} at depth {depths[iface[k] - 1]:.6g}, "
                     f"tolerance {itol:.3g})")
    # rigid base (SITE <nl> = 0): no point-load solution there
    if hm.rigid_base:
        on_base = [int(inter[k]) for k in range(n) if iface[k] == nI and mis[k] <= itol and not above[k]]
        if on_base:
            _err(hm, f"EDU-01: interaction nodes {on_base[:20]} lie on the rigid base (SITE <nl> = 0, "
                     f"interface {nI}): POINT has no solution there")
    # fixed translations (FIXEDINT, Error 124)
    active = [0, 2] if hm.dim == 1 else [0, 1, 2]
    no_dof = []
    for k in range(n):
        node = int(inter[k])
        idx = res.dofmap.node_index(node)
        sfix = hm.symm_fix[idx] if hm.symm_fix is not None else np.zeros(6, bool)
        fixed = [("UX", "UY", "UZ")[c] for c in active if hm.fix[idx, c] and not sfix[c]]
        nsym = sum(1 for c in active if sfix[c])
        if fixed and len(fixed) == len(active):
            _err(hm, f"Error 124: interaction node {node} has all translations fixed (D-CHK-06)")
        elif fixed:
            _warn(hm, f"FIXEDINT: interaction node {node} has fixed translations {fixed}")
        if nsym and len(fixed) + nsym >= len(active):
            no_dof.append(node)
    if no_dof:
        hm.info.append(f"SYMM: {len(no_dof)} interaction nodes on the symmetry planes have all their translations "
                       f"constrained (no interaction DOF; zero force and displacement by symmetry): "
                       f"{hl.first_items(no_dof)}")
    # bottom-up ascending order (EDU-21, G-13): warning for coherent analyses
    if n > 1:
        bad = np.flatnonzero(np.diff(z) < -hm.tol)
        if bad.size:
            k = int(bad[0])
            _warn(hm, f"EDU-21: interaction nodes are not numbered bottom-up: node {inter[k + 1]} "
                      f"(z = {z[k + 1]:.6g}) follows node {inter[k]} (z = {z[k]:.6g}); {bad.size} inversion(s) "
                      "(G-13; an error for incoherent analyses)")
    _structure_rules(hm, inter, z, itol)


def _structure_rules(hm: HouseModel, inter: np.ndarray, z: np.ndarray, itol: float) -> None:
    """G-14, EXCSTRCHK, EDU-22 and buried shells (requirements 4.4 item 3, 03 Item 13)."""
    iset = set(int(k) for k in inter)
    exc_elems = [(r.code, r.nodes) for r in hm.records if is_excavated_soil(r)]
    z_of = {int(k): float(v[2]) for k, v in hm.node_xyz.items()}
    ex = hl.excavation_sets(exc_elems, z_of, hm.gelev, hm.tol)
    struct_nodes: Dict[int, List[str]] = {}       # every non-excavated element node -> its elements
    exc_owners: Dict[int, List[str]] = {}         # the same without buried shells (EXCSTRCHK owners)
    exc_elems: Dict[int, List[Tuple[int, tuple]]] = {}   # EXCSTRCHK owners as (code, nodes); K node code 0
    buried: set = set()
    for r, et in zip(hm.records, hm.etype):
        if is_excavated_soil(r):
            continue
        dn, gn = element_dof_nodes(r)
        what = f"{ELEMENT_TYPE_NAMES[r.code]} element {r.id} group {r.group}"
        is_buried = r.code in hl.SHELL_TYPES and et == 2
        if is_buried:
            buried.update(dn)
        for k in dict.fromkeys(dn):
            struct_nodes.setdefault(int(k), []).append(what)
            if not is_buried:
                exc_owners.setdefault(int(k), []).append(what)
                exc_elems.setdefault(int(k), []).append((int(r.code), tuple(int(n) for n in dn)))
        for k in gn[len(dn):]:                                 # BEAMS / GENERAL K node
            struct_nodes.setdefault(int(k), []).append(what + " (K node)")
            exc_owners.setdefault(int(k), []).append(what + " (K node)")
            exc_elems.setdefault(int(k), []).append((0, (int(k),)))
    # EXCSTRCHK (G-15): interior excavation nodes connected to the structure.  Owners are the
    # elements with ETYPE != 2 and every BEAMS/SPRING/GENERAL element (spec 09 section 3.1, as
    # sassi/prep/check.py:excstrchk): a buried shell (ETYPE 2) lies *in* the soil and is covered by
    # the "buried shell nodes not interaction nodes" warning below.  Severity: hl.excstrchk_kind, the
    # rule CHECK applies (an error unless the node interacts and only SOLID/PLANE elements built on the
    # excavation mesh use it; final audit).
    n_ok = []
    for k in sorted(ex["interior"]):
        owners = exc_owners.get(k)
        if not owners:
            continue
        if hl.excstrchk_kind(exc_elems.get(k, []), ex["nodes"], k in iset) == "Warning":
            n_ok.append(k)
        elif k in iset:
            _err(hm, f"EXCSTRCHK: excavation interior node {k} is shared with the structure "
                     f"({'; '.join(owners[:3])}); as an interaction node it ties that structure to the soil "
                     f"continuum (residual of X_ff - excavated soil), which corrupts flexible subsystems "
                     f"(manual 1.5.1 rule 11): give the structure its own nodes inside the excavation")
        else:
            _err(hm, f"EXCSTRCHK: excavation interior node {k} is shared with the structure "
                     f"({'; '.join(owners[:3])}) and is not an interaction node")
    if n_ok:
        _warn(hm, f"EXCSTRCHK: {len(n_ok)} excavation interior nodes are shared with SOLID/PLANE structure "
                  f"elements built on the excavation mesh (interaction nodes): accepted -- the FV results are "
                  f"close for near-field soil or a solid basement and exact when the elements reproduce the "
                  f"excavated soil, but separate structure and excavation nodes are the correct model (manual "
                  f"1.5.1 rule 11): {hl.first_items(n_ok)}")
    # G-14: structural interaction nodes below grade must belong to the excavation (or a buried shell)
    g14 = [int(k) for k, zz in zip(inter, z) if int(k) in struct_nodes and int(k) not in ex["nodes"]
           and int(k) not in buried and hm.gelev - zz > itol]
    if g14:
        _warn(hm, f"G-14: interaction nodes below grade that belong to the structure but not to the "
                  f"excavated soil: {hl.first_items(g14)}")
    # interaction nodes that are only orientation (K) nodes
    kint = [int(k) for k in inter if int(k) in struct_nodes and all(o.endswith("(K node)")
                                                                   for o in struct_nodes[int(k)])]
    if kint:
        _warn(hm, f"KINT: interaction nodes used only as beam/GENERAL orientation nodes: {hl.first_items(kint)}")
    # buried shells: their nodes are interaction nodes (spec 08 4.6)
    nb = sorted(buried - iset)
    if nb:
        _warn(hm, f"buried shell (ETYPE 2) nodes that are not interaction nodes: {hl.first_items(nb)}")
    # EDU-22: excavation nodes that should be interaction nodes
    if ex["nodes"]:
        need = ex["nodes"] if hm.imp == 0 else ex["fsin"]
        miss = sorted(set(need) - iset)
        if miss:
            what = "FV needs every excavation node" if hm.imp == 0 else "lateral/bottom excavation surface (FSIN)"
            _warn(hm, f"EDU-22: {len(miss)} excavation nodes are not interaction nodes ({what}): "
                      f"{hl.first_items(miss)}")
    hm.excavation = ex


# ======================================================================================
# Node-numbering optimizer (HOUSEX <optimize> = 1; requirements 4.4 item 5, D-HOU-03, D-FIL-09)
# ======================================================================================
@dataclass
class NodeRenumbering:
    """Outcome of the HOUSE node optimizer: ``old_ids[k] -> new_ids[k]`` for every node of the deck."""

    old_ids: np.ndarray
    new_ids: np.ndarray
    result: RN.RenumberResult
    n_interaction: int
    files: List[str] = field(default_factory=list)

    @property
    def mapping(self) -> Dict[int, int]:
        return {int(o): int(n) for o, n in zip(self.old_ids, self.new_ids)}


#: node-number columns of the HOUSE deck tables (renumbered by the optimizer)
_NODE_COLUMNS = {"nodes": ("id",), "interaction": ("id",), "masses": ("node",), "symm": ("n1", "n2", "n3"),
                 "elements": tuple(f"n{k}" for k in range(1, 9))}


def renumbered_deck(d: decks.Deck, mapping: Dict[int, int]) -> decks.Deck:
    """Copy of a HOUSE deck with every node number replaced by ``mapping[old]`` (0 = unused slot is kept):
    the ``.hounew`` deck of the optimizer.  The ``nodes`` table is sorted by the new numbers (HOUSE numbers
    the equations in node-table order), the interaction table keeps its order, ``optimize`` is set to 0
    (the deck is already optimised).  A multiple-excitation zone (``me``: the interaction nodes with
    numbers ``nfirst..nlast``) is rewritten with the new numbers when its interaction nodes are exactly
    the interaction nodes of a range of new numbers; otherwise it keeps the original range (HOUSE resolves
    the zones on the original numbering before renumbering, so the analysis is not affected; only a
    later run of the ``.hounew`` deck would be -- :func:`_renumber_me_rows` lists such zones)."""
    out = decks.new("HOUSE")
    out.params = dict(d.params)
    out.params["optimize"] = 0
    for name, tab in d.tables.items():
        cols = list(tab.columns)
        rows = [list(r) for r in tab.rows]
        for c in _NODE_COLUMNS.get(name, ()):
            j = cols.index(c)
            for r in rows:
                v = int(r[j])
                if v > 0:
                    if v not in mapping:
                        raise ModuleError(f"node optimizer: node {v} of table '{name}' is not defined")
                    r[j] = mapping[v]
        if name == "nodes":
            rows.sort(key=lambda r: int(r[0]))
        if name == "me" and rows:
            inter = [int(r["id"]) for r in d.rows("interaction")]
            rows, _ = _renumber_me_rows(cols, rows, inter, mapping)
        out.tables[name] = Table(columns=cols, rows=rows)
    return out


def _renumber_me_rows(cols: List[str], rows: List[list], inter: List[int], mapping: Dict[int, int]
                      ) -> Tuple[List[list], List[int]]:
    """ME zone rows in the new numbering where the zone's interaction nodes form a range of new numbers
    containing no other interaction node; returns (rows, motion numbers kept in the old numbering)."""
    jn, j1, j2 = cols.index("no"), cols.index("nfirst"), cols.index("nlast")
    new_inter = {mapping.get(n, n) for n in inter}
    kept = []
    out = []
    for r in rows:
        r = list(r)
        n1, n2 = int(r[j1]), int(r[j2])
        members = sorted(mapping.get(n, n) for n in inter if n1 <= n <= n2)
        if members and {n for n in new_inter if members[0] <= n <= members[-1]} == set(members):
            r[j1], r[j2] = members[0], members[-1]
        else:
            kept.append(int(r[jn]))
        out.append(r)
    return out, kept


def optimize_model(hm: HouseModel, res: HouseResult) -> NodeRenumbering:
    """New node numbers of a checked and assembled model (reverse Cuthill-McKee from the interaction
    nodes, which keep their bottom-up order; :func:`sassi.core.renumber.optimize_numbering`).  The RCM
    order is compared with the original one on the assembled equation matrix ``|K_s| + |K_e|`` of the model
    (bandwidth and profile, spec 05b test 8)."""
    ids = np.asarray(hm.node_ids, dtype=np.int64)
    index = {int(n): k for k, n in enumerate(ids)}
    cliques = [[index[int(n)] for n in dict.fromkeys(element_dof_nodes(r)[0])] for r in hm.records]
    keep = [index[int(n)] for n in hm.interaction]
    dm = res.dofmap
    eq_node = np.array([index[int(n)] for n in dm.eq_node], dtype=np.int64)
    pattern = (abs(res.Ks) + abs(res.Ke)).tocsr()
    result = RN.optimize_numbering(ids.size, cliques, keep=keep, z=hm.xyz[:, 2] if ids.size else None,
                                   equations=(eq_node, np.asarray(dm.eq_dof, dtype=np.int64)), pattern=pattern)
    return NodeRenumbering(old_ids=ids, new_ids=result.new_number.astype(np.int64), result=result,
                           n_interaction=len(keep))


def renumbered_model(d: decks.Deck, sit: Optional[decks.Deck], hm: HouseModel,
                     renum: NodeRenumbering) -> Tuple[decks.Deck, HouseModel, HouseResult]:
    """The ``.hounew`` deck and its model, assembled: what FILE4, COOSK, COOSM are written from (manual
    section 3: "the .hounew model is what generates FILE4").  The model was checked in the original
    numbering; renumbering cannot create an error, so any error here is an internal one.  The element
    state set for this run on the original records (non-linear soil materials, incompatible modes) is
    carried over element by element (same order and numbers)."""
    d2 = renumbered_deck(d, renum.mapping)
    hm2 = build_house_model(d2, sit)
    if len(hm2.records) != len(hm.records) or any((a.group, a.id, a.code) != (b.group, b.id, b.code)
                                                  for a, b in zip(hm.records, hm2.records)):
        raise ModuleError("node optimizer: the renumbered model has different elements (internal error)")
    for a, b in zip(hm.records, hm2.records):
        b.mat = a.mat
        b.props = dict(a.props)
    try:
        res2 = assemble_house(hm2)
    except ElementError as exc:
        raise ModuleError(f"node optimizer: renumbered model: {exc} (internal error)") from None
    interaction_set(hm2, res2)
    if hm2.errors:
        raise ModuleError("node optimizer: the renumbered model fails the HOUSE checks (internal error): "
                          + "; ".join(hm2.errors[:3]))
    return d2, hm2, res2


def write_optimizer_files(ctx: ModuleContext, d2: decks.Deck, renum: NodeRenumbering) -> List[str]:
    """``<model>.hounew`` (the renumbered deck, ``optimize`` = 0) and ``<model>.map`` (``old new`` pairs,
    D-FIL-09), written together with FILE4 once every check has passed."""
    hn = ctx.path(f"{ctx.model}.hounew")
    decks.write(hn, d2, comments=[f"HOUSE node-numbering optimizer ({renum.result.method}): renumbered model of "
                                  f"{ctx.model}.hou", f"old -> new node numbers in {ctx.model}.map"])
    mp = RN.write_map(ctx.path(f"{ctx.model}.map"), renum.old_ids, renum.new_ids)
    renum.files = [hn.name, mp.name]
    return [f"{hn.name} (node-numbering optimizer: renumbered model)", f"{mp.name} (node-numbering optimizer: old new)"]


def _list_optimizer(lst, ctx: ModuleContext, hm: HouseModel, renum: NodeRenumbering, write: bool) -> None:
    """Listing of the node optimizer and the manual's warning on the new numbers."""
    r = renum.result
    lst.section("Node-numbering optimizer (HOUSEX <optimize> = 1, reverse Cuthill-McKee, D-HOU-03)")
    lst.write(f" Method                        : {r.method}")
    lst.write(f" Nodes renumbered              : {renum.old_ids.size} (new numbers 1..{renum.old_ids.size}; "
              f"{r.n_isolated} node(s) without element DOFs numbered last)")
    what = "equation matrix |Ks| + |Ke|" if r.level == "equation" else "node graph"
    lst.write(f" Bandwidth before/after        : {r.bandwidth_before} / {r.bandwidth_after}  ({what})")
    lst.write(f" Profile before/after          : {r.profile_before} / {r.profile_after}"
              + (f"  ({100.0 * (1.0 - r.profile_after / r.profile_before):.1f} % smaller)"
                 if r.profile_before > 0 and r.reduced else ""))
    lst.write(f" Interaction nodes             : {renum.n_interaction}, relative (bottom-up) order kept"
              + (" (last numbers of their part, ascending)" if r.method == "RCM" and renum.n_interaction else ""))
    for t in r.notes:
        lst.write(f" {t}")
    lst.write(f" This listing uses the node numbers of {ctx.model}.hou; FILE4, COOSK, COOSM, DOFSMAP and the FILE8 "
              f"of ANALYS use those of {ctx.model}.hounew")
    if write:
        lst.write(f" Renumbered model              : {ctx.model}.hounew;  old -> new pairs: {ctx.model}.map "
                  "(written with FILE4)")
    else:
        lst.write(" Data check only: .hounew and .map are not written")
    changed = int(np.sum(renum.old_ids != renum.new_ids))
    if changed:
        _warn(hm, f"node-numbering optimizer: {changed} node numbers changed -- FILE4, COOSK, COOSM and the FILE8 of "
                  f"ANALYS use the NEW numbers of {ctx.model}.hounew (old -> new pairs in {ctx.model}.map): select "
                  "MOTION/RELDISP output nodes with the new numbers (manual section 3); FORCE loads (FILE9) keep "
                  "the model numbering and ANALYS translates them through the old -> new pairs (FILE4 "
                  "x_node_old_id); STRESS element numbers do not change")
    else:
        lst.write(" Node numbers unchanged (the .map is the identity)")


def _retire_optimizer_files(ctx: ModuleContext) -> List[str]:
    """A run without the optimizer writes FILE4 in the original numbering: a ``<model>.map`` / ``.hounew``
    left by an earlier optimised run no longer describes it and is renamed ``.prev`` (as FILE78)."""
    out = []
    for ext in (".map", ".hounew"):
        p = ctx.path(f"{ctx.model}{ext}")
        if p.exists():
            q = p.with_name(p.name + ".prev")
            p.replace(q)
            out.append(f"{q.name} (renamed: {p.name} of an earlier optimised run does not describe this FILE4)")
    return out


# ======================================================================================
# Module entry point
# ======================================================================================
def _p1_options(d: decks.Deck, hm: HouseModel) -> None:
    """Options of later tiers (requirements 0.1): stop or warn.  Incoherency, wave passage and multiple
    excitation are implemented (:func:`_incoherency_plan`), SYMM planes (:func:`apply_symmetry`) and 2D
    models too; Option AA (P2) stops the run.  Incoherent motion, wave passage and multiple excitation
    need a full 3D model: they are errors with SYMM planes or HOUSE <dim> = 1 (manual 2.7 and 6.5.4,
    G-17; the incoherency settings check incoherent motion itself)."""
    if int(d["ansys"]) == 1:
        raise ModuleError("ANSYS model input (Option AA) not available in this build (tier P2)")
    keys = [k for k in ("coh", "wpass", "me") if int(d[k]) == 1]
    names = {"coh": "incoherent motion", "wpass": "wave passage", "me": "multiple excitation"}
    what = " and ".join(names[k] for k in keys) + f" (HOUSE <{'>/<'.join(keys)}>)"
    if keys and d.rows("symm"):
        _err(hm, f"EDU-26: {what} not allowed with SYMM planes: incoherent analysis, wave passage and multiple "
                 "excitation need the full model (manual 2.7 and 6.5.4, G-17)")
    elif keys and hm.dim == 1 and keys != ["coh"]:       # incoherent motion alone: the incoherency settings say it
        _err(hm, f"EDU-26: {what} need a 3D model (HOUSE <dim> = 2; manual 6.5.4, G-17)")


def run(ctx: ModuleContext) -> int:
    lst = ctx.listing
    d = decks.read(ctx.deck_path, NAME)
    # ---- D-HOU-02: the .sit deck is always required ------------------------------------------
    sit_path = ctx.path(f"{ctx.model}.sit")
    if not sit_path.exists():
        raise ModuleError(f"HOUSE needs {ctx.model}.sit (run AFWRITE with SITE enabled)")
    try:
        sit = decks.read(sit_path, "SITE")
    except (ValueError, KeyError) as exc:
        raise ModuleError(f"{sit_path.name}: {exc}") from None
    dim = int(d["dim"])
    if dim not in (1, 2):
        raise ModuleError(f"<dim> = {dim}: HOUSE needs 1 (2D) or 2 (3D); 1D analysis is not available")
    ctx.progress(0.05, "HOUSE: reading the deck")
    hm = build_house_model(d, sit)
    _p1_options(d, hm)
    optimize = int(d["optimize"]) == 1
    nlplan = _nonlinear_plan(ctx, hm) if int(d["nlssi"]) == 1 and not hm.errors else None
    g = hm.gravity
    if g > 0 and min(abs(g - 32.174) / 32.174, abs(g - 9.80665) / 9.80665) > 0.05:
        _warn(hm, f"acceleration of gravity {g:.6g} differs by more than 5 % from 32.174 ft/s2 and 9.80665 m/s2 "
                  "(D-CNV-08): check the unit system")
    if int(sit["soilmode"]) == 1 and any(is_excavated_soil(r) for r in hm.records):
        _warn(hm, "SITE uses the strain-compatible FILE88 properties (SITEX soil mode 1) but the excavated "
                  "elements use the L table: update the L layers with the SOIL results, otherwise the excavated "
                  "soil differs from the free field (D-ELM-12)")
    _list_input(lst, d, hm)
    if not hm.records and not hm.errors:
        _err(hm, "the HOUSE deck has no elements")
    if hm.errors:
        if optimize:
            _warn(hm, "node-numbering optimizer (HOUSEX <optimize> = 1) skipped: the input has errors")
        _fail(lst, hm)
    # ---- assembly and checks, in the numbering of the model ---------------------------------------
    ctx.progress(0.2, "HOUSE: element matrices and assembly")
    ctx.announce("HOUSE.elements", ne=len(hm.records), dim=dim)
    try:
        res = assemble_house(hm)
    except ElementError as exc:
        raise ModuleError(str(exc)) from None
    interaction_set(hm, res)
    ctx.announce("HOUSE.assembled", ne=len(hm.records), neq=int(res.Ks.shape[0]), nnz=int(res.Ks.nnz),
                 nnze=int(res.Ke.nnz), nint=int(len(res.int_node)))
    _dimension_checks(hm)
    _rotation_checks(hm, res)
    inco = _incoherency_plan(ctx, d, hm, res)          # incoherency / wave passage / ME (None when off)
    for m in res.messages:
        hm.info.append(m)
    if res.ignored_masses:
        _warn(hm, f"{len(res.ignored_masses)} nodal mass terms on fixed or undefined DOFs ignored (W5/W6): "
                  + ", ".join(f"node {n} dof {k}" for n, k, _ in res.ignored_masses[:10]))
    # ---- node-numbering optimizer (requirements 4.4 item 5): only a checked model is renumbered ----
    renum: Optional[NodeRenumbering] = None
    if optimize:
        if hm.errors:
            _warn(hm, "node-numbering optimizer (HOUSEX <optimize> = 1) skipped: the model has errors")
        else:
            ctx.progress(0.5, "HOUSE: node-numbering optimizer")
            renum = optimize_model(hm, res)
            _list_optimizer(lst, ctx, hm, renum, write=int(d["opmode"]) == 0)
    if inco is not None:
        _list_incoherency_input(lst, inco, renum is not None)
    _list_results(lst, d, hm, res)
    if nlplan is not None:
        _list_nonlinear(lst, ctx, hm, nlplan)
    if hm.errors:
        _fail(lst, hm)
    if int(d["opmode"]) == 1:
        lst.write("")
        lst.write(" Data check only (<opmode> = 1): no output files written")
        return 0
    paths: List[str] = []
    if renum is not None:
        if np.any(renum.old_ids != renum.new_ids):
            ctx.progress(0.6, "HOUSE: assembly of the renumbered model")
            d, hm, res = renumbered_model(d, sit, hm, renum)
        else:                                          # identity map: the checked model is the renumbered one
            d = renumbered_deck(d, renum.mapping)
    ctx.progress(0.8, "HOUSE: writing FILE4, COOSK, COOSM, DOFSMAP, FILE90, FILE91")
    paths += write_outputs(ctx, d, hm, res, renum)
    if renum is not None:
        paths += write_optimizer_files(ctx, d, renum)
    else:
        paths += _retire_optimizer_files(ctx)
    if nlplan is not None:
        paths += _nonlinear_outputs(ctx, d, nlplan)
    else:
        paths += _retire_file78(ctx)
    if inco is not None:
        ctx.progress(0.85, "HOUSE: incoherency factors (FILE77)")
        paths += _incoherency_outputs(ctx, hm, res, inco)
    lst.section("Files written")
    for p in paths:
        lst.write(f" {p}")
    ctx.progress(1.0, "HOUSE done")
    return 0


def _fail(lst, hm: HouseModel) -> None:
    for w in hm.warnings:
        lst.warning(w)
    hm.warnings.clear()
    for e in hm.errors:
        lst.error(e)
    raise ModuleError(f"{len(hm.errors)} error(s) in the HOUSE input (listed above); no output files written")


def _dimension_checks(hm: HouseModel) -> None:
    """EDU-26 (D-PNT-01): element types consistent with HOUSE <dim> -- an *error*, as in CHECK
    (sassi/prep/check.py), because <dim> also selects POINT2/POINT3: a PLANE (plane-strain) element
    has no meaning in a 3D model, and SOLID/SHELL/TSHELL elements cannot be coupled to a 2D
    (plane-strain) far field."""
    codes = {r.code for r in hm.records}
    if hm.dim == 2 and 4 in codes:
        _err(hm, "EDU-26: PLANE elements in a 3D model (HOUSE <dim> = 2)")
    if hm.dim == 1 and codes & {1, 3, 5}:
        _err(hm, "EDU-26: SOLID/SHELL/TSHELL elements in a 2D model (HOUSE <dim> = 1)")
    if hm.dim == 1 and codes & {2, 7, 9}:
        hm.info.append("2D model (plane strain, X-Z plane, per unit length): BEAMS/SPRING/GENERAL nodes keep six DOFs; "
                       "the out-of-plane DOFs UY, ROTX, ROTZ are structural only (the soil acts on UX, UZ) -- fix them "
                       "(D) unless something restrains them")


def _rotation_checks(hm: HouseModel, res: HouseResult) -> None:
    """EDU-06 (D-ELM-09): active rotational DOFs without stiffness (SHELL drilling) need
    FIXROT / FIXSHLROT."""
    if res.dofmap.neq == 0 or not np.any(res.dofmap.eq_dof >= 4):
        return
    flagged = unstiffened_rotations(res.Ks, res.dofmap)
    if flagged:
        nodes = sorted({n for n, _ in flagged})
        _warn(hm, f"EDU-06: {len(nodes)} nodes have an active rotation without stiffness (e.g. SHELL drilling; "
                  f"use FIXROT/FIXSHLROT, or D at nodes shared with SOLID/PLANE elements, which FIXROT does not "
                  f"treat): " + "; ".join(
                      f"node {n} axis ({a[0]:+.3f},{a[1]:+.3f},{a[2]:+.3f})" for n, a in flagged[:8])
                  + (" ..." if len(flagged) > 8 else ""))


# ======================================================================================
# Incoherency, wave passage and multiple excitation (requirements 4.4 item 6, D-INC-*)
# ======================================================================================
#: memory budget of the incoherency factors held at once (stochastic samples are computed in chunks)
FILE77_CHUNK_BYTES = 256 * 2 ** 20
#: the I N C O table lists at most this many modes per frequency and direction
INCO_MAX_ROWS = 100


@dataclass
class IncoherencyPlan:
    """Validated incoherency / wave-passage / multiple-excitation input of one HOUSE run."""

    settings: INC.IncoherencySettings
    spec: Optional[COH.CoherencySpec]
    zones: List[INC.MEZone]
    fnum: np.ndarray
    df: float
    embedded_levels: int
    notes: List[str] = field(default_factory=list)


def _incoherency_plan(ctx: ModuleContext, d: decks.Deck, hm: HouseModel, res: HouseResult) -> Optional[IncoherencyPlan]:
    """Settings, checks and coherency model of the incoherent / wave-passage / multiple-excitation input
    (requirements 4.4 item 6; spec 05b section 8): errors and warnings go to ``hm`` (listed with the
    other checks; an error stops the run before any file is written).  None when HOUSE <coh>, <wpass>
    and <me> are all 0 (coherent analysis: no FILE77)."""
    st = INC.IncoherencySettings.from_deck(d, INC.read_options(ctx.workdir))
    if not st.active:
        if d.rows("me") or d.rows("amp"):
            _warn(hm, "ME zones / AMP ratios in the deck are not used: HOUSE <me> = 0")
        return None
    notes: List[str] = []
    # ---- the SSI frequencies (the factors are computed at every frequency of the set)
    fn = [int(r["number"]) for r in d.rows("freqs")]
    df = house_frequency_step(d)
    if not fn:
        _err(hm, "Error 120: the frequency set is empty: incoherency, wave passage and multiple excitation are "
                 "computed at the SSI frequencies of the HOUSE deck")
    elif len(set(fn)) != len(fn) or min(fn) < 1:
        _err(hm, "EDU-23: the frequency set holds duplicate or non-positive frequency numbers")
    if not df > 0:
        _err(hm, "the frequency step is not defined (HOUSE deck df, or delt and NFFT)")
    fnum = np.array(sorted(set(fn)), dtype=np.int64)
    # ---- option combination (spec 05b section 8)
    iface = np.asarray(res.int_iface, np.int64)
    levels = int(np.unique(iface[iface > 1]).size) if iface.size else 0
    errs, warns = st.validate(len(hm.interaction), hm.dim, bool(d.rows("symm")), levels)
    for e in errs:
        _err(hm, e)
    for w in warns:
        _warn(hm, w)
    if st.coh:
        # R4: interaction nodes numbered bottom-up -- an error for incoherent analyses (D-HOU-05)
        for w in [w for w in hm.warnings if w.startswith("EDU-21")]:
            hm.warnings.remove(w)
            _err(hm, w.replace("(G-13; an error for incoherent analyses)", "(G-13): an error for incoherent analyses "
                               "(manual 6.5.4: interaction nodes numbered from the base up to the ground surface)"))
    # ---- coherency model
    spec: Optional[COH.CoherencySpec] = None
    xyz = np.array([hm.node_xyz[int(n)] for n in hm.interaction], dtype=float).reshape(-1, 3)
    if st.coh and 1 <= st.cohf <= 7:
        l2m = COH.length_to_metres(hm.gravity)
        user = None
        if st.cohf == 7:
            try:
                user = COH.UserCoherency.from_dir(ctx.workdir)
                notes += user.notes
            except COH.CoherencyError as exc:
                _err(hm, f"EDU-26: {exc}")
        spec = COH.CoherencySpec(model=st.cohf, gamma=st.gamma, alpha=st.alpha, angle_deg=st.wang, length_to_m=l2m,
                                 user=user)
        try:
            spec.validate()
        except COH.CoherencyError as exc:
            _err(hm, str(exc) if str(exc).startswith("Error") else f"EDU-26: {exc}")
            spec = None
        if spec is not None and st.cohf in (2, 3, 4, 5, 6):
            si = unit_system(hm.gravity) == "SI"
            if st.met != (1 if si else 0):
                _warn(hm, f"INCOH <met> = {st.met} ({'SI' if st.met else 'British'}) differs from the unit system of the "
                          f"model (gravity {hm.gravity:g}: {'SI, metres' if si else 'British, feet'}): the Abrahamson "
                          f"model distances are converted with the model units (x {l2m:g} to metres, D-CNV-08)")
            if xyz.shape[0] > 1:
                dmax = float(spec.distances(xyz[:, :2]).max()) * l2m
                mdl = COH.load_model(st.cohf)
                lo, hi = mdl.validated_range_m
                if dmax > hi:
                    _warn(hm, f"coherency model {st.cohf}: separations up to {dmax:.4g} m exceed the published range "
                              f"{lo:g}-{hi:g} m" + (f"; beyond {mdl.formula_max_m:g} m the coherency is held at its value "
                                                    f"at {mdl.formula_max_m:g} m (conservative)"
                                                    if dmax > mdl.formula_max_m else ""))
        if spec is not None and xyz.shape[0] > 1:
            pairs = INC.close_projections(xyz, hm.tol)
            if pairs:
                i, j = pairs[0]
                _warn(hm, f"EDU-17: Near-Coincident Horizontal Projections of Interaction Nodes: {len(pairs)} pair(s) on "
                          f"different levels, e.g. nodes {hm.interaction[i]} and {hm.interaction[j]} (D-INC-09: check "
                          "the trace of Eq. 6.1 in the I N C O table; EDUOPT,INCOHMERGE,1 or the per-level approach "
                          "with BUILDFILE77)")
    # ---- multiple excitation
    zones: List[INC.MEZone] = []
    if st.me:
        zones, e2, w2 = INC.me_zones(d.rows("me"), d.rows("amp"), hm.interaction, int(fnum.size), st.cmplxspec)
        for e in e2:
            _err(hm, e)
        for w in w2:
            _warn(hm, w)
    if st.wpass and st.appv >= 1.0e8:
        notes.append(f"apparent velocity {st.appv:.3g}: the wave-passage delays are negligible "
                     f"(largest delay about {_max_extent(xyz) / st.appv:.3g} s)")
    return IncoherencyPlan(settings=st, spec=spec, zones=zones, fnum=fnum, df=float(df), embedded_levels=levels,
                           notes=notes)


def _max_extent(xyz: np.ndarray) -> float:
    if xyz.shape[0] < 2:
        return 0.0
    p = xyz[:, :2]
    return float(np.linalg.norm(p.max(axis=0) - p.min(axis=0)))


def _list_incoherency_input(lst, plan: IncoherencyPlan, optimized: bool) -> None:
    """Listing of the incoherency / wave-passage / multiple-excitation input."""
    st = plan.settings
    lst.section("Incoherency, wave passage and multiple excitation (requirements 4.4 item 6, FILE77)")
    lst.write(f" Soil motion                   : {'incoherent' if st.coh else 'coherent'} (HOUSE <coh> = {st.coh})")
    if st.coh:
        lst.write(f" Coherency model               : " + (plan.spec.describe() if plan.spec else f"<cohf> = {st.cohf}"))
        if plan.spec is not None and st.cohf in (2, 3, 4, 5, 6):
            lst.write(f"   source                      : {COH.load_model(st.cohf).citation}")
        lst.write(f" Superposition                 : {st.describe_method()}")
        lst.write(f" Incoherent modes (<nmodes>)   : "
                  + ("all" if st.nmodes == 0 else (f"the first {st.nmodes}" if st.nmodes > 0 else f"mode {-st.nmodes} only")))
        lst.write(f" Mode sign convention          : {st.sign}" + (" (sum of every mode >= 0, D-INC-03)" if st.sign == "ADJUST"
                                                                   else " (eigensolver signs, EDUOPT,INCOHSIGN,RAW)"))
        lst.write(f" Plan positions                : " + ("merged within the geometric tolerance (EDUOPT,INCOHMERGE,1)"
                                                         if st.merge else "every interaction node (D-INC-09)"))
        lst.write(f" Embedded interaction levels   : {plan.embedded_levels} (INCOH <ngp> = {st.ngp})")
    lst.write(f" Wave passage                  : " + (f"V_app = {st.appv:.6g} along Line D at {st.wang:g} deg from X; delay "
                                                     f"measured from the control point ({st.xc:g}, {st.yc:g}) (D-INC-06)"
                                                     if st.wpass else "none"))
    if st.me:
        lst.write(f" Multiple excitation           : {len(plan.zones)} zone(s), "
                  f"{'complex' if st.cmplxspec else 'real'} spectral amplification ratios (D-INC-10)")
        for z in plan.zones:
            lst.write(f"   motion {z.no:>4d}: interaction nodes {z.nfirst}..{z.nlast} ({z.nodes.size} nodes); "
                      f"|SAR| {np.abs(z.sar).min():.4g} .. {np.abs(z.sar).max():.4g}")
        if optimized:
            lst.write("   (zones resolved on the node numbers of the .hou deck, before the node optimizer)")
    lst.write(f" SSI frequencies               : {plan.fnum.size} (df = {plan.df:.9g} Hz)")
    out = (f"FILE77001 .. {INC.file77_name(st.nsim)} (one per simulation; ANALYS <simul> = {st.nsim} with FILE1X/Y/Z)"
           if st.stochastic else "FILE77")
    lst.write(f" Output                        : {out}")
    for t in plan.notes:
        lst.write(f" {t}")


def _retire_file77(ctx: ModuleContext, keep: List[str]) -> List[str]:
    """Rename the incoherency files of an earlier HOUSE run that this run does not rewrite ``.prev``
    (a deterministic run retires FILE77001..; a stochastic run with Ns samples retires FILE77 and the
    samples above Ns), so that ANALYS cannot mix them with the new FILE77."""
    out = []
    names = ["FILE77"] + [INC.file77_name(k) for k in range(1, 1000)]
    for nm in names:
        if nm in keep:
            continue
        p = ctx.path(nm)
        if p.exists():
            q = p.with_name(nm + ".prev")
            p.replace(q)
            out.append(f"{q.name} (renamed: {nm} of an earlier HOUSE run is not part of this incoherency input)")
    return out


def _incoherency_outputs(ctx: ModuleContext, hm: HouseModel, res: HouseResult, plan: IncoherencyPlan) -> List[str]:
    """Compute and write FILE77 (deterministic, wave passage, ME) or FILE77001 .. FILE77Ns (stochastic) on
    the interaction nodes of the written FILE4 (the renumbered ones when the optimizer ran), with the
    decomposition summary and the I N C O tables in the listing (requirements 4.4 item 6)."""
    lst = ctx.listing
    st = plan.settings
    int_node = np.asarray(res.int_node, np.int64)
    xyz = np.array([hm.node_xyz[int(n)] for n in int_node], dtype=float).reshape(-1, 3)
    N = int(int_node.size)
    fnum = plan.fnum
    freq = fnum * plan.df
    nF = int(fnum.size)
    sar = INC.sar_factors(plan.zones, nF, N) if st.me else None
    run = INC.FactorRun(settings=st, spec=plan.spec, int_xyz=xyz, fnum=fnum, freq=freq, sar=sar, merge_tol=hm.tol)
    sims = list(range(st.nsim)) if st.stochastic else [0]
    names = [INC.file77_name(r + 1) for r in sims] if st.stochastic else [INC.file77_name()]
    out = _retire_file77(ctx, names)
    try:
        f90 = read_container_meta(ctx.path("FILE90"))
    except (OSError, ValueError):
        f90 = {}
    per_sample = max(1, nF * 3 * N * 16)
    chunk = max(1, int(FILE77_CHUNK_BYTES // per_sample))
    tables: List[str] = []
    printing = {"on": bool(st.ipr) and bool(st.coh)}

    def on_dec(q: int, d_: int, dec: INC.Decomposition) -> None:
        if printing["on"]:
            INC.inco_table(tables.append, q, int(fnum[q]), float(freq[q]), d_, dec, INCO_MAX_ROWS)

    stats = None
    info: Dict[str, object] = {}
    first_s = None
    for c0 in range(0, len(sims), chunk):
        part = sims[c0:c0 + chunk]
        try:
            S, st_q, inf = INC.compute_factors(run, part, on_dec)
        except (INC.IncoherencyError, COH.CoherencyError) as exc:
            raise ModuleError(f"incoherency factors: {exc}") from None
        if stats is None:
            stats, info, first_s = st_q, inf, S[0]
        printing["on"] = False                        # the I N C O tables once (first chunk)
        for j, r in enumerate(part):
            meta = _file77_meta(ctx, st, plan, info, f90, r + 1 if st.stochastic else None)
            extra = _file77_stats(stats) if st.coh else None
            INC.write_file77(ctx.path(names[c0 + j]), fnum, freq, int_node, xyz, S[j], meta, extra, module=NAME)
        if ctx.cancelled():
            raise ModuleError("run cancelled")
    _list_incoherency_results(lst, plan, stats, info, first_s, freq)
    for t in tables:
        lst.write(t)
    if st.stochastic:
        out.insert(0, f"{names[0]} .. {names[-1]} ({len(names)} stochastic incoherency samples, X/Y/Z factors)")
    else:
        out.insert(0, f"FILE77 ({st.method} factors of {N} interaction nodes, {nF} frequencies, X/Y/Z)")
    return out


def _file77_meta(ctx: ModuleContext, st: INC.IncoherencySettings, plan: IncoherencyPlan, info: dict, f90: dict,
                 sim: Optional[int]) -> Dict[str, object]:
    return {"df": plan.df, "method": st.method, "model": ctx.model, "coh": st.coh, "wpass": st.wpass, "me": st.me,
            "cohf": st.cohf if st.coh else 0, "coherency_model": plan.spec.describe() if plan.spec else "",
            "nmodes": st.nmodes, "mode": -st.nmodes if (st.coh and st.nmodes < 0) else 0, "supmode": st.supmode,
            "stochastic": int(st.stochastic), "sim": sim, "nsim": st.nsim if st.stochastic else 1,
            "hseed": st.hseed, "vseed": st.vseed, "randphz": st.randphz, "appv": st.appv, "wang": st.wang,
            "xc": st.xc, "yc": st.yc, "alpha": st.alpha, "gamma": list(st.gamma), "sign": st.sign, "merge": st.merge,
            "n_positions": int(info.get("n_positions", 0)) if info else 0,
            "int_hash": str(f90.get("int_hash", "")), "file4_hash": str(f90.get("file4_hash", "")),
            "me_zones": [[z.no, z.nfirst, z.nlast] for z in plan.zones],
            "content": "s[q, d, i]: factor of interaction node i (int_node), direction d (0 X, 1 Y, 2 Z) at frequency "
                       "q; incoherent free field u_i = s_i u_i,coherent (FFM) or load P_i = s_i P_i,coherent (FFL)"}


def _file77_stats(stats) -> Dict[str, np.ndarray]:
    """Decomposition summary arrays of FILE77 (SASSI-EDU additions, ``x_`` prefix)."""
    nF = len(stats)
    lam1 = np.array([[s.lam1_pct for s in row] for row in stats], float).reshape(nF, 3)
    n90 = np.array([[s.n90 for s in row] for row in stats], np.int64).reshape(nF, 3)
    tre = np.array([[s.trace_error for s in row] for row in stats], float).reshape(nF, 3)
    return {"x_lambda1_pct": lam1, "x_n90": n90, "x_trace_error": tre}


def _list_incoherency_results(lst, plan: IncoherencyPlan, stats, info: dict, s0: Optional[np.ndarray],
                              freq: np.ndarray) -> None:
    """Per-frequency summary of the spectral factorisation (lambda_1 share, modes for 90 % of the
    variance, Eq. 6.1 trace check) and the low-frequency factor check."""
    st = plan.settings
    lst.section("Incoherency factors: spectral factorisation summary (Eqs. 6.1-6.4)")
    if stats is None or not st.coh:
        lst.write(" coherent motion: no coherency decomposition (FILE77 holds the wave-passage / SAR factors)")
    else:
        lst.write(f" {int(info.get('n_positions', 0))} plan positions; lambda_1 % = 100 lambda_1 / N; n90 = modes for 90 % "
                  "of the variance; trace error = |sum(lambda) - N| / N")
        lst.write(f"{'q':>5s}{'number':>8s}{'f (Hz)':>11s}" + "".join(f"{c + ' lam1 %':>13s}{'n90':>6s}" for c in "XYZ")
                  + f"{'trace err':>11s}")
        for q in range(len(stats)):
            row = stats[q]
            lst.write(f"{q + 1:>5d}{int(plan.fnum[q]):>8d}{freq[q]:>11.5g}"
                      + "".join(f"{s.lam1_pct:>13.4f}{s.n90:>6d}" for s in row)
                      + f"{max(s.trace_error for s in row):>11.2e}")
        worst = float(info.get("max_trace_error", 0.0))
        lst.write(f" largest trace error {worst:.2e} (Eq. 6.1 check: < {INC.TRACE_TOL:g})")
        if worst > INC.TRACE_TOL:
            lst.warning(f"Eq. 6.1 check failed (|sum(lambda) - N| / N = {worst:.2e} > {INC.TRACE_TOL:g}): the coherency "
                        f"matrix is not positive semi-definite or ill-conditioned (most negative eigenvalue "
                        f"{float(info.get('min_raw', 0.0)):.3g} N, set to 0); see the I N C O table (INCOH <ipr> = 1)")
    if s0 is not None and s0.size:
        a = np.abs(s0[0])
        lst.write(f" factors at the lowest SSI frequency f = {freq[0]:.5g} Hz: |s| = {a.min():.6g} .. {a.max():.6g} "
                  "(X, Y, Z; deterministic AS -> 1 as f -> 0, the zero-frequency ATF check of the manual)")

def _nonlinear_plan(ctx: ModuleContext, hm: HouseModel) -> Optional[NL.NLPlan]:
    """Properties of the nonlinear soil elements for this run (:func:`sassi.core.nlsoil.plan_iteration`).

    Reads ``<model>.pin``, ``<model>.liq``, FILE73 (curves), and on later iterations FILE74 (effective
    strains) and the previous FILE78 (properties used, for the convergence measure).  The free field
    at every nonlinear element (the TOPL layer containing its mid-depth, ``sitelayers`` row with the
    L properties, as for the excavated soil, D-ELM-12) is the reference of the iteration-0 factors:
    ``G = GFAC G_layer``, ``beta = DFAC beta_layer`` (requirements 4.4 item 7); the element's M-table
    material is the low-strain soil (G_max of the curves).  The material of every nonlinear element is
    replaced by its own internal material (one per element); nonlinear soil elements are assembled
    without incompatible modes, like the excavated soil, so that equal properties give exactly the
    excavated-soil matrices (zero-SSI identity, VP-16).  Without a ``.pin`` file there are no nonlinear
    groups: a warning (CHECK EDU-42) and a linear run (a FILE78 of an earlier run is then retired,
    :func:`_retire_file78`)."""
    pin_path = ctx.path(f"{ctx.model}.pin")
    if not pin_path.exists():
        _warn(hm, f"non-linear soil SSI (HOUSEX <nlssi> = 1): {pin_path.name} not found -- the model is analysed "
                  "with its linear properties and FILE78 is not written (define the nonlinear groups with PINGRP / "
                  "PINMAT and AFWRITE, or edit the .pin with the HOUSE dialog Input Data)")
        return None
    try:
        pin = NL.read_pin(pin_path)
        liq = NL.read_liq(ctx.path(f"{ctx.model}.liq"))
        curves = NL.load_curves(ctx.path("FILE73")) if ctx.path("FILE73").exists() else None
        file74 = NL.read_nlfile(ctx.path("FILE74"), "FILE74") if liq == 1 and ctx.path("FILE74").exists() else None
    except (NL.NLSoilError, OSError) as exc:
        raise ModuleError(f"non-linear soil input: {exc}") from None
    prev78 = None
    if liq == 1 and ctx.path("FILE78").exists():
        try:
            prev78 = NL.read_nlfile(ctx.path("FILE78"), "FILE78")
        except NL.NLSoilError as exc:
            _warn(hm, f"previous FILE78 not readable ({exc}): no convergence measure")
    groups = {g.igrp for g in pin.groups}
    elements = []
    layer_cache: Dict[int, MaterialProps] = {}
    for k, (r, row) in enumerate(zip(hm.records, hm.rows)):
        e = NL.NLElementIn(k, r.group, r.id, r.code, is_excavated_soil(r), int(row["mat"]), r.mat)
        if r.group in groups and r.code in hl.SOIL_TYPES:
            e.layer, e.layer_no, e.depth = _free_field_layer(hm, r, layer_cache)
        elements.append(e)
    try:
        plan = NL.plan_iteration(pin, elements, hm.dim, liq, file74, curves, prev78)
    except NL.NLSoilError as exc:
        raise ModuleError(f"non-linear soil SSI: {exc}") from None
    n_inc = 0
    for e in plan.elements:
        rec = hm.records[e.index]
        rec.mat = e.material
        if rec.props.get("incompatible"):
            rec.props["incompatible"] = False
            n_inc += 1
    if n_inc:
        hm.info.append(f"non-linear soil SSI: incompatible modes switched off for {n_inc} near-field soil elements "
                       "(formulation of the excavated soil)")
    for w in plan.warnings:
        _warn(hm, w)
    if curves is None:
        _warn(hm, "FILE73 (soil curves) not found: STRESS <iter> = 1 and the next iterations need it -- run SOIL")
    return plan


def _free_field_layer(hm: HouseModel, r: ElemRecord,
                      cache: Dict[int, MaterialProps]) -> Tuple[Optional[MaterialProps], int, float]:
    """Free field at a near-field soil element (requirements 4.4 item 7, ``G_layer``/``beta_layer``):
    ``(material, L number, mid-depth)`` of the SITE layer table row (``sitelayers``: TOPL layers with
    their L properties + the half-space row) containing the mid-depth of the element (centroid of its
    distinct nodes, the assignment of the excavated soil check EDU-08).  Above grade: TOPL layer 1.
    ``(None, 0, depth)`` without a layer table.  ``cache`` keeps one material per table row."""
    zs = [hm.z(n) for n in dict.fromkeys(x for x in r.nodes if x)]
    depth = hm.gelev - float(np.mean(zs)) if zs else float("nan")
    if not hm.site or not np.isfinite(depth):
        return None, 0, depth
    idx = min(max(hl.layer_index_at_depth(hm.depths, depth), 1), len(hm.site))
    row = hm.site[idx - 1]
    if idx not in cache:
        try:
            cache[idx] = material_from_layer(float(row["vp"]), float(row["vs"]), float(row["weight"]),
                                             float(row["dp"]), float(row["ds"]), hm.gravity, hm.form)
        except ElementError as exc:
            raise ModuleError(f"non-linear soil SSI: free-field layer {int(row['no'])} at depth {depth:.5g}: {exc}") from None
    return cache[idx], int(row["no"]), depth


def _retire_file78(ctx: ModuleContext) -> List[str]:
    """A linear HOUSE run (<nlssi> = 0, or no .pin) writes a FILE4 that a FILE78 left from an earlier
    non-linear run does not describe: it is renamed FILE78.prev so that STRESS <iter> = 1 and the
    convergence checks cannot use it (the .liq iteration flag is not changed)."""
    return [f"{n} (FILE78 of an earlier non-linear run: it does not describe this FILE4)"
            for n in NL.retire(ctx.workdir, ("FILE78",))]


def _nonlinear_outputs(ctx: ModuleContext, d: decks.Deck, plan: NL.NLPlan) -> List[str]:
    """FILE78 (properties used, with the FILE4 hash and the deck and .pin digests), ``<model>.liq`` = 1
    and the convergence row."""
    try:
        file4 = str(read_container_meta(ctx.path("FILE90")).get("file4_hash", ""))
    except (OSError, ValueError):
        file4 = ""
    deck = NL.deck_digest(ctx.deck_path) if ctx.deck_path else ""
    NL.write_nlfile(ctx.path("FILE78"), plan.file78(model=ctx.model, file4=file4, deck=deck), title=str(d["title"]))
    out = ["FILE78"]
    liq = ctx.path(f"{ctx.model}.liq")
    if NL.read_liq(liq) != 1:
        NL.write_liq(liq, 1)
        out.append(f"{liq.name} (1: the next HOUSE run reads FILE74)")
    if plan.convergence is not None:
        NL.record_convergence(ctx.workdir, plan.convergence)
        out.append(NL.CONVERGENCE_FILE)
    return out


def read_container_meta(path) -> dict:
    """Meta data of a (small) container file, e.g. FILE90."""
    from ..io.container import read_container
    return dict(read_container(path).meta)


def _list_nonlinear(lst, ctx: ModuleContext, hm: HouseModel, plan: NL.NLPlan) -> None:
    """Listing of the non-linear soil iteration: .pin echo, element properties, group summary and the
    convergence check of the previous iteration (D-NLS-02)."""
    pin, curves = plan.pin, plan.curves
    lst.section(f"Non-linear soil SSI: iteration {plan.iteration} (HOUSEX <nlssi> = 1, {ctx.model}.pin, "
                "requirements 4.4 item 7)")
    lst.write(f" Properties from               : " + ("the .pin initial factors and the free field: G = GFAC G_layer, "
                                                      "beta = DFAC beta_layer (TOPL layer at the element mid-depth)"
                                                      if plan.source == "PIN" else
                                                      f"{plan.source} and the FILE73 curves: G = G_max (G/G_max)"
                                                      "(gamma_eff), beta = D(gamma_eff)"))
    lst.write(" G_max of the curves           : the low-strain modulus of the element's M-table material (column G_max)")
    lst.write(f" Effective strain factor ESF   : {pin.esf:g}")
    lst.write(f" Curves NCURV (FILE73)         : {pin.ncurv}" + (
        "  (" + ", ".join(f"{k} {cv.label}" for k, cv in sorted(curves.items())) + ")" if curves else "  (FILE73 not found)"))
    lst.write(f"{'group':>7s}{'type':>7s}{'elements':>10s}{'ISTR':>6s}  strain measure")
    for gi in plan.groups:
        lst.write(f"{gi.group:>7d}{gi.etype:>7s}{gi.nelem:>10d}{gi.istr:>6d}  {NL.STRAIN_FLAGS[gi.istr]}")
    lst.write(f"{'group':>7s}{'material':>10s}{'GFAC':>8s}{'DFAC':>8s}{'ICURVE':>8s}  {'curve':<10s}{'G_max':>13s}"
              f"{'beta_mat':>10s}")
    seen = set()
    for e in plan.elements:
        if (e.group, e.mat) in seen:
            continue
        seen.add((e.group, e.mat))
        lst.write(f"{e.group:>7d}{e.mat:>10d}{e.gfac:>8.4g}{e.dfac:>8.4g}{e.curve:>8d}  "
                  f"{NL.curve_label(curves, e.curve):<10s}{e.base.G0:>13.5g}{e.base.beta_s:>10.4f}")
    prev = plan.previous.by_key() if plan.previous is not None else {}
    lst.section(f"Non-linear soil elements: properties of iteration {plan.iteration} (one internal material per element)")
    first = plan.source == "PIN"
    if first:
        lst.write(" layer = L number of the free-field TOPL layer at the element mid-depth whose G and damping GFAC and "
                  "DFAC scale ('mat': no layer table, the element material)")
    lst.write(f"{'group':>6s}{'elem':>7s}{'mat':>5s}{'z centre':>11s}" + (f"{'layer':>6s}" if first else "")
              + f"{'g_eff %':>12s}{'G/Gmax':>9s}{'G':>13s}{'Vs':>10s}{'beta':>8s}"
              + (f"{'dG/G %':>9s}{'dbeta %':>9s}" if prev else ""))
    for k, e in enumerate(plan.elements):
        if k >= MAX_LAYER_ROWS:
            lst.write(f" ... {len(plan.elements) - MAX_LAYER_ROWS} more elements (FILE78 lists them all)")
            break
        r = hm.records[e.index]
        zc = float(np.mean([hm.z(n) for n in dict.fromkeys(x for x in r.nodes if x)]))
        m = e.material
        lay = (f"{(str(e.ref_layer) if e.ref_layer else 'mat'):>6s}" if first else "")
        line = (f"{e.group:>6d}{e.element:>7d}{e.mat:>5d}{zc:>11.5g}{lay}{e.gamma_eff:>12.5g}{e.ratio:>9.4f}"
                f"{m.G0:>13.5g}{m.vs:>10.5g}{m.beta_s:>8.4f}")
        p = prev.get((e.group, e.element))
        if p is not None:
            dG = (m.G0 - p.G) / m.G0 * 100.0
            line += f"{dG:>9.3f}{(m.beta_s - p.beta) * 100.0:>9.3f}"
        lst.write(line)
    lst.write(f"{'group':>7s}{'elements':>10s}{'G/Gmax min':>12s}{'mean':>8s}{'max':>8s}{'beta min':>10s}{'mean':>8s}"
              f"{'max':>8s}")
    for g, n, r0, r1, r2, b0, b1, b2 in NL.group_summary(plan.rows()):
        lst.write(f"{g:>7d}{n:>10d}{r0:>12.4f}{r1:>8.4f}{r2:>8.4f}{b0:>10.4f}{b1:>8.4f}{b2:>8.4f}")
    lst.section("Non-linear soil convergence (D-NLS-02: max |dG/G| < 2 %, max |d beta| < 0.5 % absolute)")
    if plan.convergence is not None:
        c = plan.convergence
        lst.write(f" Properties used in iteration {c.iteration} (FILE78) -> strain-compatible properties of its response "
                  f"(FILE74) = the properties of iteration {plan.iteration}:")
        lst.write(f"   {c.text()}")
        if c.converged:
            lst.write(f"   the results of iteration {c.iteration} are strain-compatible: the iterations have converged "
                      f"(this run repeats them with properties changed by less than the tolerances)")
    elif plan.iteration == 0:
        lst.write(" iteration 0: initial properties from the .pin (no strains yet); run ANALYS, STRESS (<iter> = 1) and, "
                  "for X/Y/Z input, COMBXYZSTRAIN, then HOUSE again")
    for t in plan.notes:
        lst.write(f" {t}")
    hist = NL.read_convergence(ctx.workdir)
    if plan.convergence is not None:
        hist = [h for h in hist if h["iteration"] != plan.convergence.iteration] + [
            {"iteration": plan.convergence.iteration, "max_dG_pct": plan.convergence.dG_pct,
             "max_dbeta_pct": plan.convergence.dbeta_pct, "mean_G_Gmax": plan.convergence.mean_ratio,
             "mean_beta": plan.convergence.mean_beta, "converged": int(plan.convergence.converged)}]
    hist = sorted(hist, key=lambda h: h["iteration"])
    if hist:
        lst.write(f"{'iteration':>10s}{'max dG/G %':>13s}{'max dbeta %':>13s}{'mean G/Gmax':>13s}{'mean beta':>11s}"
                  f"  converged")
        for h in hist:
            lst.write(f"{int(h['iteration']):>10d}{h['max_dG_pct']:>13.3f}{h['max_dbeta_pct']:>13.3f}"
                      f"{h['mean_G_Gmax']:>13.4f}{h['mean_beta']:>11.4f}  {'yes' if h['converged'] else 'no'}")


# ======================================================================================
# Output files
# ======================================================================================
def _site_array(hm: HouseModel) -> np.ndarray:
    cols = ("no", "thick", "weight", "vp", "vs", "dp", "ds")
    return np.array([[float(r[c]) for c in cols] for r in hm.site], dtype=float).reshape(-1, 7)


def file4_arrays(hm: HouseModel, res: HouseResult) -> Tuple[Dict[str, np.ndarray], Dict[str, list]]:
    """FILE4 arrays (files.SCHEMAS['FILE4']) and the recovery component lists.

    The ``elem_*`` arrays hold the output-capable elements (types with STRESS components) in deck
    order; ``rec_<T>_idx`` points into them."""
    dm = res.dofmap
    out_types = {name for name in res.recovery}
    elem_rows = []
    for r, row, et in zip(hm.records, hm.rows, hm.etype):
        name = ELEMENT_TYPE_NAMES[r.code]
        if name not in out_types:
            continue
        elem_rows.append((r, row, et))
    nE = len(elem_rows)
    pos = {(r.group, r.id): k for k, (r, _, _) in enumerate(elem_rows)}
    arr: Dict[str, np.ndarray] = {
        "node_id": hm.node_ids.astype(np.int64),
        "node_xyz": hm.xyz.astype(float),
        "eq_node": dm.eq_node.astype(np.int64),
        "eq_dof": dm.eq_dof.astype(np.int64),
        "int_node": res.int_node.astype(np.int64),
        "int_iface": res.int_iface.astype(np.int64),
        "int_eq": res.int_eq.astype(np.int64).reshape(-1, 3),
        "elem_group": np.array([r.group for r, _, _ in elem_rows], dtype=np.int64),
        "elem_id": np.array([r.id for r, _, _ in elem_rows], dtype=np.int64),
        "elem_type": np.array([r.code for r, _, _ in elem_rows], dtype=np.int64),
        "elem_nodes": np.array([list(r.nodes)[:8] for r, _, _ in elem_rows], dtype=np.int64).reshape(nE, 8),
        "elem_excav": np.array([int(is_excavated_soil(r)) for r, _, _ in elem_rows], dtype=np.int64),
        # SASSI-EDU additions (x_ prefix, files.py)
        "x_elem_etype": np.array([et for _, _, et in elem_rows], dtype=np.int64),
        "x_elem_mat": np.array([int(row["mat"]) for _, row, _ in elem_rows], dtype=np.int64),
        "x_elem_prop": np.array([int(row["prop"]) for _, row, _ in elem_rows], dtype=np.int64),
        "x_node_fix": hm.fix.astype(np.int64),
        "x_iface_depth": hm.depths.astype(float),
        "x_sitelayers": _site_array(hm),
        "x_int_xyz": np.array([hm.node_xyz[int(k)] for k in res.int_node], dtype=float).reshape(-1, 3),
    }
    if hm.symm:                         # SYMM planes and the interaction translations they constrain (D-ANL-12)
        arr["x_symm"] = SYM.planes_array(hm.symm)
        rows = [dm.node_index(int(k)) for k in res.int_node]
        sfix = hm.symm_fix if hm.symm_fix is not None else np.zeros((hm.node_ids.size, 6), bool)
        arr["x_int_symm"] = np.asarray(sfix[rows, :3] if rows else np.zeros((0, 3)), dtype=np.int64).reshape(-1, 3)
    comps: Dict[str, list] = {}
    for name, rd in res.recovery.items():
        idx = np.array([pos[(int(g), int(e))] for g, e in zip(rd["group"], rd["id"])], dtype=np.int64)
        arr[f"rec_{name}_S"] = rd["S"].astype(complex)
        arr[f"rec_{name}_eq"] = rd["eq"].astype(np.int64)
        arr[f"rec_{name}_idx"] = idx
        if "B" in rd:
            arr[f"x_rec_{name}_B"] = rd["B"].astype(complex)
        comps[name] = list(ELEMENT_COMPONENTS[name])
        if name == "TSHELL":
            # THSHLSTR face stresses/strains need t, E, nu per element (D-TSH-01; SASSI-EDU x_ arrays)
            recs = {(r.group, r.id): r for r in hm.records if r.code == 5}
            rows = [recs[(int(g), int(e))] for g, e in zip(rd["group"], rd["id"])]
            arr["x_rec_TSHELL_thick"] = np.array([float(r.props["thick"]) for r in rows], dtype=float)
            arr["x_rec_TSHELL_E"] = np.array([float(r.mat.E0) for r in rows], dtype=float)
            arr["x_rec_TSHELL_nu"] = np.array([float(r.mat.nu0) for r in rows], dtype=float)
            arr["x_rec_TSHELL_eint"] = np.array([int(r.props.get("eint", 0)) for r in rows], dtype=np.int64)
    return arr, comps


def restart_hashes(hm: HouseModel, res: HouseResult, file4: Dict[str, np.ndarray]) -> Dict[str, str]:
    """FILE90 keys (D-W1-13, requirements 2.5):

    * ``int_hash``   interaction-node ids, coordinates and interfaces (X_ff can be reused); for 2D models
      and SYMM half / quarter models also the dimension and the planes (they change X_ff);
    * ``layer_hash`` site layer table + gravity + complex-modulus form (FILE1/FILE3 consistency);
    * ``file4_hash`` FILE4 arrays + COOSK/COOSM matrices (factorised system can be reused)."""
    xyz = np.array([hm.node_xyz[int(k)] for k in res.int_node], dtype=float).reshape(-1, 3)
    items = [("int_node", res.int_node.astype(np.int64)), ("int_xyz", xyz), ("int_iface", res.int_iface.astype(np.int64))]
    if hm.dim == 1 or hm.symm:          # 3D models without SYMM keep their earlier hash
        items += [("dim", np.int64(hm.dim)), ("symm", SYM.planes_array(hm.symm))]
    ih = hl.stable_hash(items)
    lh = hl.stable_hash([("sitelayers", _site_array(hm)), ("gravity", np.float64(hm.gravity)),
                         ("cmodform", np.int64(hm.form))])
    items = [(k, file4[k]) for k in sorted(file4)]
    for name, A in (("Ks", res.Ks), ("Ke", res.Ke), ("Ms", res.Ms), ("Me", res.Me)):
        items += hl.sparse_items(name, A)
    return {"int_hash": ih, "layer_hash": lh, "file4_hash": hl.stable_hash(items)}


def write_outputs(ctx: ModuleContext, d: decks.Deck, hm: HouseModel, res: HouseResult,
                  renum: Optional[NodeRenumbering] = None) -> List[str]:
    """Write FILE4 (<model>.N4), COOSK, COOSM, DOFSMAP, FILE90 and FILE91; return the names.

    With the node optimizer FILE4 is in the new numbering and also carries ``x_node_old_id`` (the
    original number of every ``node_id``) and the meta keys ``x_optimized`` = 1 and ``x_node_map``
    (the ``.map`` file name), so that downstream modules can translate or check node requests."""
    arr, comps = file4_arrays(hm, res)
    if renum is not None:
        back = {int(n): int(o) for o, n in zip(renum.old_ids, renum.new_ids)}
        arr["x_node_old_id"] = np.array([back[int(n)] for n in arr["node_id"]], dtype=np.int64)
    hashes = restart_hashes(hm, res, arr)
    neq = res.dofmap.neq
    afw = deck_model_hash(ctx.deck_path) if ctx.deck_path else ""
    meta4 = {"gravity": hm.gravity, "gelev": hm.gelev, "dim": hm.dim, "nEq": neq, "cmodform": hm.form,
             "model_hash": hashes["file4_hash"], "afwrite_model_hash": afw, "components": comps, "model": ctx.model,
             "title": str(d["title"]), "units": unit_system(hm.gravity), "imp": int(d["imp"]),
             "incomp": int(d["incomp"]), "gmunits": int(d["gmunits"]), "nI": int(hm.depths.size),
             "tol": hm.tol, "int_dofs": list(interaction_dofs(hm.dim)),
             "excav_layering": {k: list(v) for k, v in hm.layering.items()} or {"R5": [], "R6": []},
             "int_eq": "int_eq columns UX, UY, UZ; only the int_dofs translations are interaction DOFs "
                       "(2D: UX, UZ; the UY column is -1)",
             "recovery":"component TF = rec_<T>_S[e] @ u[rec_<T>_eq[e]] (eq -1 = fixed, 0)",
             "elem_nodes": "node slots as in the E command (BEAMS/GENERAL: I, J, K; K is orientation only)",
             "shell_membrane": "stress (D-STR-10, D-W1-14)"}
    if renum is not None:
        meta4.update(x_optimized=1, x_node_map=f"{ctx.model}.map", x_renumber_method=renum.result.method)
    if hm.symm:
        meta4["x_symm_planes"] = [p.describe(hm.dim) for p in hm.symm]
        meta4["x_symm"] = ("x_symm rows: no, type (0 symmetry, 1 antisymmetry), normal axis (0 X, 1 Y), coordinate; "
                           "x_int_symm: interaction translations UX UY UZ constrained by the planes (D-ANL-12)")
    c4 = Container("FILE4", meta4, arr)
    probs = validate(c4)
    if probs:                                          # pragma: no cover - programming error
        raise ModuleError("FILE4: " + "; ".join(probs))
    name4 = file4_name(ctx.model)
    write_container(ctx.path(name4), "FILE4", arr, meta4, module=NAME)
    common = {"nEq": neq, "model_hash": hashes["file4_hash"], "model": ctx.model}
    write_container(ctx.path("COOSK"), "COOSK", {"Ks": res.Ks, "Ke": res.Ke},
                    dict(common, content="Ks: structure incl. near-field soil (complex K*); Ke: excavated soil"),
                    module=NAME)
    write_container(ctx.path("COOSM"), "COOSM", {"Ms": res.Ms, "Me": res.Me},
                    dict(common, content="Ms: structure + nodal masses; Me: excavated soil (real)"), module=NAME)
    write_container(ctx.path("DOFSMAP"), "DOFSMAP",
                    {"eq_node": arr["eq_node"], "eq_dof": arr["eq_dof"], "int_node": arr["int_node"],
                     "x_int_eq": arr["int_eq"], "x_int_iface": arr["int_iface"]}, common, module=NAME)
    write_container(ctx.path("FILE90"), "FILE90", {}, dict(hashes, model=ctx.model), module=NAME)
    counts = {ELEMENT_TYPE_NAMES[c]: int(sum(1 for r in hm.records if r.code == c))
              for c in sorted({r.code for r in hm.records})}
    meta91 = {"model": ctx.model, "title": str(d["title"]), "program": f"{PRODUCT} {__version__} HOUSE",
              "created": time.strftime("%Y-%m-%d %H:%M:%S"), "n_nodes": int(hm.node_ids.size),
              "n_elements": counts, "n_excavated": int(sum(1 for r in hm.records if is_excavated_soil(r))),
              "n_interaction": int(res.int_node.size), "nEq": neq, "gravity": hm.gravity, "gelev": hm.gelev,
              "dim": hm.dim, "imp": int(d["imp"]), "units": unit_system(hm.gravity),
              "mass_structure": res.mass_total_s.tolist(), "mass_excavated": res.mass_total_e.tolist(),
              "hashes": hashes, "afwrite_model_hash": afw}
    write_container(ctx.path("FILE91"), "FILE91", {}, meta91, module=NAME)
    return [name4, "COOSK", "COOSM", "DOFSMAP", "FILE90", "FILE91"]


# ======================================================================================
# Listing
# ======================================================================================
def _list_input(lst, d: decks.Deck, hm: HouseModel) -> None:
    g = hm.gravity
    lst.section("HOUSE input")
    lst.write(f" Title                         : {d['title']}")
    lst.write(f" Acceleration of gravity       : {g:.6g} ({'British' if unit_system(g) == 'BS' else 'SI'} units)"
              if g > 0 else f" Acceleration of gravity       : {g}")
    lst.write(f" Ground elevation              : {hm.gelev:.6g}")
    lst.write(f" Dimension                     : {'3D' if hm.dim == 2 else '2D'}  (<dim> = {hm.dim})")
    lst.write(f" Method                        : {('FV', 'FFV', 'FI')[min(max(int(d['imp']), 0), 2)]}"
              f"  (<imp> = {int(d['imp'])})")
    lst.write(f" Soil motion                   : {'coherent' if int(d['coh']) == 0 else 'incoherent'}")
    lst.write(f" Incompatible modes            : {'included' if int(d['incomp']) == 0 else 'suppressed'}"
              " (structural SOLID/PLANE only)")
    lst.write(f" GENERAL mass units (MOPT)     : {'weight (divided by g)' if int(d['gmunits']) == 1 else 'mass'}")
    lst.write(f" Complex-modulus form          : {'1 + 2 i beta' if hm.form == 1 else '1 - 2 beta^2 + 2 i beta sqrt(1 - beta^2)'}")
    lst.write(f" Geometric tolerance           : {hm.tol:.3g} (D-GEN-06)")
    lst.write(f" Nodes                         : {hm.node_ids.size}")
    # counts by type
    lst.section("Elements by type")
    lst.write(f"{'type':>10s}{'structure':>12s}{'excavated':>12s}{'buried':>10s}")
    for code in sorted({r.code for r in hm.records}):
        recs = [(r, et) for r, et in zip(hm.records, hm.etype) if r.code == code]
        nx = sum(1 for r, _ in recs if is_excavated_soil(r))
        nb = sum(1 for r, et in recs if code in hl.SHELL_TYPES and et == 2)
        lst.write(f"{ELEMENT_TYPE_NAMES[code]:>10s}{len(recs) - nx:>12d}{nx:>12d}{nb:>10d}")
    lst.section("Groups")
    lst.write(f"{'group':>7s}{'type':>9s}{'elements':>10s}{'excavated':>11s}  title")
    for gid, (code, title) in sorted(hm.groups.items()):
        recs = [r for r in hm.records if r.group == gid]
        nx = sum(1 for r in recs if is_excavated_soil(r))
        lst.write(f"{gid:>7d}{ELEMENT_TYPE_NAMES.get(code, str(code)):>9s}{len(recs):>10d}{nx:>11d}  {title}")
    if hm.materials:
        lst.section("Structural materials (M table; rho = weight / g, complex moduli D-CNV-04)")
        lst.write(f"{'mat':>5s}{'E0':>13s}{'nu0':>8s}{'Re G*':>13s}{'Im G*':>13s}{'Re M*':>13s}{'Im M*':>13s}"
                  f"{'rho':>12s}{'beta_p':>8s}{'beta_s':>8s}")
        for k, m in sorted(hm.materials.items()):
            lst.write(f"{k:>5d}{m.E0:>13.5g}{m.nu0:>8.4f}{m.G.real:>13.5g}{m.G.imag:>13.5g}{m.M.real:>13.5g}"
                      f"{m.M.imag:>13.5g}{m.rho:>12.5g}{m.beta_p:>8.4f}{m.beta_s:>8.4f}")
    if hm.layers:
        lst.section("Soil layers of the excavated elements (L table, D-ELM-12)")
        lst.write(f"{'layer':>6s}{'Vs':>11s}{'Vp':>11s}{'Re G*':>13s}{'Im G*':>13s}{'Re M*':>13s}{'Im M*':>13s}"
                  f"{'rho':>12s}{'beta_p':>8s}{'beta_s':>8s}")
        for k, m in sorted(hm.layers.items()):
            lst.write(f"{k:>6d}{m.vs:>11.5g}{m.vp:>11.5g}{m.G.real:>13.5g}{m.G.imag:>13.5g}{m.M.real:>13.5g}"
                      f"{m.M.imag:>13.5g}{m.rho:>12.5g}{m.beta_p:>8.4f}{m.beta_s:>8.4f}")
    if hm.site:
        lst.section("Site layer table (TOPL layers + half-space) and user interfaces")
        dep = hm.depths
        lst.write(f"{'iface':>6s}{'depth':>11s}{'elev.':>11s}{'layer':>7s}{'thick':>10s}{'Vs':>10s}{'Vp':>10s}"
                  f"{'weight':>10s}")
        for k, r in enumerate(hm.site):
            last = k == len(hm.site) - 1
            lst.write(f"{k + 1:>6d}{dep[k]:>11.5g}{hm.gelev - dep[k]:>11.5g}{int(r['no']):>7d}"
                      f"{'half-sp.' if last else format(float(r['thick']), '.5g'):>10s}{float(r['vs']):>10.5g}"
                      f"{float(r['vp']):>10.5g}{float(r['weight']):>10.5g}")
    if hm.mass_rows:
        tot = np.sum([m for m in hm.masses.values()], axis=0) if hm.masses else np.zeros(6)
        lst.section("Nodal masses (MT/MR; MUNITS 1 = weight, divided by g)")
        lst.write(f" {len(hm.mass_rows)} rows on {len(hm.masses)} nodes; sum of translational masses "
                  f"({tot[0]:.6g}, {tot[1]:.6g}, {tot[2]:.6g}), rotary ({tot[3]:.6g}, {tot[4]:.6g}, {tot[5]:.6g})")
    if hm.n_etype0:
        lst.write("")
        lst.write(f" ETYPE 0 resolved for {hm.n_etype0} SOLID/PLANE elements: {hm.n_etype0_exc} excavated, "
                  f"{hm.n_etype0 - hm.n_etype0_exc} structure (D-HOU-01)")


def _list_symmetry(lst, hm: HouseModel) -> None:
    """SYMM planes of a half / quarter model with the boundary conditions applied (D-ANL-12)."""
    if not hm.symm:
        return
    what = "half" if len(hm.symm) == 1 else "quarter"
    lst.section(f"Symmetry planes (SYMM, D-ANL-12): {what} model")
    lst.write(f"{'plane':>6s}{'loading':>14s}{'location':>16s}{'nodes on it':>13s}  DOFs fixed on the plane")
    on = SYM.on_planes(hm.symm, hm.xyz, hm.tol)
    for k, p in enumerate(hm.symm):
        loc = f"{SYM.AXIS_NAMES[p.axis].lower()} = {p.coord:.6g}"
        lst.write(f"{p.no:>6d}{SYM.TYPE_NAMES[p.type]:>14s}{loc:>16s}{int(on[:, k].sum()):>13d}  "
                  + " ".join(SYM.DOF_NAMES[dd - 1] for dd in p.fixed_dofs))
    lst.write(" The symmetry conditions are added to the D fixities of the plane nodes.  Elements, masses and loads "
              "lying IN a plane must carry 1/2 (1/4 on two planes) of their full-model values.  ANALYS forms the soil")
    lst.write(" impedance of the reduced interaction set with image interaction nodes; incoherent motion, wave passage, "
              "multiple excitation and the global impedance need the full model (manual 2.7, G-17).")


def _list_results(lst, d: decks.Deck, hm: HouseModel, res: HouseResult) -> None:
    dm = res.dofmap
    _list_symmetry(lst, hm)
    _layer_assignment(lst, d, hm)
    lst.section("Interaction nodes per user interface")
    if res.int_node.size == 0:
        lst.write(" none (fixed-base or vibration model without soil impedance)")
    else:
        dep = hm.depths
        lst.write(f"{'iface':>6s}{'depth':>11s}{'layer below':>13s}{'nodes':>7s}  node ids")
        for i in sorted(set(int(k) for k in res.int_iface)):
            sel = res.int_node[res.int_iface == i]
            lay = int(hm.site[i - 1]["no"]) if i - 1 < len(hm.site) else -1
            lst.write(f"{i:>6d}{dep[i - 1]:>11.5g}{lay:>13d}{sel.size:>7d}  {hl.first_items(sel)}")
        lst.write(f" {res.int_node.size} interaction nodes, {int((res.int_eq >= 0).sum())} interaction equations "
                  "(translations only, D-CNV-09)")
    lst.section("Model matrices")
    lst.write(f" Number of equations           : {dm.neq}")
    lst.write(f" DOFs fixed by D               : {int((dm.defined & dm.fixed).sum())}")
    lst.write(f" DOFs eliminated (no element)  : {dm.n_eliminated} (D-ELM-09, information)")
    lst.write(f" Non-zeros Ks / Ke             : {res.Ks.nnz} / {res.Ke.nnz}")
    lst.write(f" Non-zeros Ms / Me             : {res.Ms.nnz} / {res.Me.nnz}")
    ms, me, mf = res.mass_total_s, res.mass_total_e, res.mass_free_s
    lst.write(f" Total structural mass X Y Z   : {ms[0]:.6g} {ms[1]:.6g} {ms[2]:.6g}"
              f"  (weight {ms[2] * hm.gravity:.6g})")
    lst.write(f"   of which on active DOFs     : {mf[0]:.6g} {mf[1]:.6g} {mf[2]:.6g}")
    lst.write(f" Total excavated soil mass     : {me[0]:.6g} {me[1]:.6g} {me[2]:.6g}"
              f"  (weight {me[2] * hm.gravity:.6g})")
    for t in hm.info:
        lst.write(f" {t}")
    lst.section("Checks")
    if not hm.warnings and not hm.errors:
        lst.write(" no warnings")
    for w in hm.warnings:
        lst.warning(w)
    hm.warnings.clear()


def house_frequency_step(d: decks.Deck) -> float:
    """Frequency step of the HOUSE frequency set: the resolved deck ``df`` if > 0, otherwise
    ``1/(delt NFFT)`` (same rule as SITE/FORCE, conventions.frequency_step); 0 when neither is
    given (the frequency-dependent checks such as EDU-10 are then skipped)."""
    df = float(d["df"])
    if df > 0:
        return df
    delt, nft = float(d["delt"]), int(d["nft"])
    return frequency_step(delt, nft, 0.0) if delt > 0 and nft > 0 else 0.0


def _layer_assignment(lst, d: decks.Deck, hm: HouseModel) -> None:
    """Excavated-element layer assignment table (05b R7) with EDU-08 and EDU-10 checks."""
    exc = [(r, row) for r, row in zip(hm.records, hm.rows) if is_excavated_soil(r)]
    if not exc:
        return
    dep = hm.depths
    fn = [int(r["number"]) for r in d.rows("freqs")]
    df = house_frequency_step(d)
    fcut = max(fn) * df if fn and df > 0 else 0.0
    lst.section("Excavated soil elements: layer assignment (check against the site profile, 05b R7)")
    lst.write(f"{'group':>6s}{'elem':>7s}{'MSET':>6s}{'mid-depth':>11s}{'TOPL at depth':>15s}{'Vs':>10s}"
              f"{'rho':>11s}{'beta_s':>8s}{'height':>9s}{'f_pass':>9s}  check")
    n08 = n10 = 0
    placed: List[Tuple[int, int, int, int]] = []          # (group, element, TOPL row, layer no) for R5/R6
    for k, (r, row) in enumerate(exc):
        ns = list(dict.fromkeys(n for n in r.nodes if n))
        zs = np.array([hm.z(n) for n in ns])
        dmid = hm.gelev - float(zs.mean())
        h = float(zs.max() - zs.min())
        lid = int(row["mat"])
        m = hm.layers.get(lid)
        expect = None
        if dep.size:
            idx = hl.layer_index_at_depth(dep, dmid)
            if 1 <= idx <= len(hm.site):
                expect = int(hm.site[idx - 1]["no"])
                placed.append((r.group, r.id, idx, expect))
        flags = []
        if expect is not None and expect != lid:
            flags.append("EDU-08")
            n08 += 1
            if n08 <= MAX_REPEAT:
                _warn(hm, f"EDU-08: excavated element {r.id} of group {r.group} uses layer {lid}, the TOPL layer at "
                          f"its mid-depth {dmid:.5g} is {expect}")
        fpass = m.vs / (5.0 * h) if (m is not None and h > 0) else float("inf")
        if fcut > 0 and fpass < fcut * (1 - 1e-9):
            flags.append("EDU-10")
            n10 += 1
        if k < MAX_LAYER_ROWS:
            lst.write(f"{r.group:>6d}{r.id:>7d}{lid:>6d}{dmid:>11.5g}{('-' if expect is None else str(expect)):>15s}"
                      f"{(m.vs if m else 0.0):>10.5g}{(m.rho if m else 0.0):>11.5g}{(m.beta_s if m else 0.0):>8.4f}"
                      f"{h:>9.4g}{min(fpass, 9.99e8):>9.4g}  {' '.join(flags) if flags else 'ok'}")
    if len(exc) > MAX_LAYER_ROWS:
        lst.write(f" ... {len(exc) - MAX_LAYER_ROWS} more excavated elements not listed")
    if n08 > MAX_REPEAT:
        _warn(hm, f"EDU-08: {n08} excavated elements in total use a layer different from the TOPL layer at their depth")
    if n10:
        _warn(hm, f"EDU-10: {n10} excavated elements are taller than Vs/(5 f_cut) at f_cut = {fcut:.4g} Hz "
                  "(1/5-wavelength rule, G-06)")
    _excavation_groups(lst, hm, placed)


def _excavation_groups(lst, hm: HouseModel, placed: List[Tuple[int, int, int, int]]) -> None:
    """Excavation group summary and the numbering/grouping rules R5/R6 (spec 05b section 1.5,
    D-HOU-04): warnings here; the findings go to FILE4 meta ``excav_layering`` so that STRESS can
    stop when soil-pressure or nodal-contour output is requested (strict rule of the manual)."""
    hm.layering = hl.excavation_layering(placed)
    if not placed:
        return
    lst.section("Excavation groups (rules R5/R6: one group per embedment layer, numbered from the surface down)")
    lst.write(f"{'group':>6s}{'elements':>10s}{'element ids':>16s}  layers (top first)")
    by_group: Dict[int, List[Tuple[int, int, int]]] = {}
    for g, e, idx, lay in placed:
        by_group.setdefault(g, []).append((e, idx, lay))
    for g in sorted(by_group):
        rows = by_group[g]
        ids = [e for e, _, _ in rows]
        lays = [lay for _, lay in sorted({(i, lay) for _, i, lay in rows})]
        lst.write(f"{g:>6d}{len(rows):>10d}{f'{min(ids)}..{max(ids)}':>16s}  {' '.join(str(x) for x in lays)}")
    for rule in ("R5", "R6"):
        for m in hm.layering[rule][:MAX_REPEAT]:
            _warn(hm, f"{rule}: {m} (D-HOU-04: a warning for the SSI analysis; required for STRESS soil-pressure "
                      "and nodal-contour output)")
    if not hm.layering["R5"] and not hm.layering["R6"]:
        lst.write(" R5/R6 satisfied")


if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
