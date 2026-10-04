"""LOADGEN: ANSYS loads and boundary conditions from the SSI solution (Option A, the ACS SASSI-ANSYS
two-step approach).

Sources: manual section 6.4.15 (reference/acs-sassi.txt lines 4730-4870), spec 04 section 15.6 (the
dialogs "ANSYS Static Load Converter" and "ANSYS Dynamic Load Converter"), spec 01 (Option A),
requirements sections 1.1 (Option A, LOADGEN) and 5.3 (Modules menu: ANSYS Eq. Static Load / ANSYS
Dynamic Load).  The exact LOADGEN input format is in the separate "ACS SASSI-ANSYS Integration
Capability" manual, which is not available (spec 04 OQ 20): the deck, the file layouts and the APDL
written here are SASSI-EDU interpretations, labelled D-LGN-xx below and in docs/user/OPTION_A.md.

Why a two-step analysis (for the structural engineer)
-----------------------------------------------------
Step 1 is the SSI analysis (SITE -> POINT -> HOUSE -> ANALYS -> MOTION / RELDISP / STRESS): it gives
the motion of every node of the SSI model, soil-structure interaction included.  Step 2 is a refined
ANSYS model of the structure (finer mesh, other element types, nonlinear materials, contact ...)
loaded by that motion.  Option A assumes that the step-2 refinements do not change the motion of the
foundation-soil interface (manual, Introduction).  For a linear structure the decomposition is exact:
split the total displacement as ``u = u_r + iota u_g`` with ``iota u_g`` a rigid-body translation of
the reference (ground) motion.  Rigid-body motion produces no stiffness forces (``K iota = 0``), so the
equations of the structural DOFs s with the interface DOFs b prescribed become::

    M_ss u_r,s'' + C u_r,s' + K_ss u_r,s = -M_ss iota_s a_g(t) - K_sb u_r,b(t)

which is exactly an ANSYS transient run with ``ACEL = a_g(t)`` (the acceleration of the reference
frame, D-LGN-04) and the prescribed displacements ``D = u_r,b(t)`` of the interface nodes relative to
the reference motion -- the inputs of the "ANSYS Dynamic Load" dialog (ground acceleration file and the
relative-displacement files of RELDISP).  The only approximation of step 2 is the damping model: ANSYS
direct integration uses Rayleigh damping ``C = alpha M + beta K`` (``zeta(w) = alpha/(2 w) + beta w/2``)
while SASSI uses frequency-independent hysteretic damping; the manual calls this "a significant
limitation" of the dynamic second step.

The equivalent static loads of the "ANSYS Eq. Static Load" dialog freeze the same equation at a
critical time t* and drop the velocity terms: ``K_ss u_s = -M_ss a_s(t*) - K_sb u_b(t*)``, i.e. the
nodal inertia forces ``F = -m a(t*)`` of the *absolute* accelerations (MOTION ``.ACC``) and, at the
interface nodes, the displacements relative to the free field (RELDISP ``.THD``).  The relative and the
total interface displacements differ by a rigid translation, which produces no stress.  Without the
interface displacements ("Acceleration") the interface nodes are fixed: the classical fixed-base
equivalent static analysis.  By the dynamic equilibrium of the free body above the interface, the sum of
the inertia forces equals the SSI base shear (VP-LA1).

What LOADGEN does (D-LGN-01 ... D-LGN-08)
----------------------------------------
1. reads its deck ``<model>.lgn`` (written by RUNLOADGEN, :mod:`sassi.prep.commands.loadgen_cmds`), FILE4
   (``<model>.N4``: nodes, equations, interaction nodes, the HOUSE optimizer numbering) and the HOUSE
   deck (the dialog's "HOUSE Module Input");
2. maps the nodes: FILE4 numbering -> model numbering (``x_node_old_id`` of the HOUSE node optimizer,
   else identity) -> ANSYS numbering (identity, as the ``ANSYS`` command and ``CONVERT,ANSYS`` keep the
   node numbers; or a ``sassi ansys`` pairs file; or the ANSYS node at the same coordinates of an ANSYS
   ``.cdb``/APDL file, D-LGN-05).  2-D models (PLANE groups) use the rotation of the ANSYS export:
   SASSI X-Z plane -> ANSYS X-Y plane (``UZ -> UY``, ``ROTY -> -ROTZ``);
3. masses (static, D-LGN-03): **Lumped Mass** -- per node ``m_d = sum_j M[(i,d), (j,d)]`` (row sums of the
   HOUSE structure mass matrix COOSM ``Ms`` over the translations of the same direction: the total mass
   per direction is exact) and the rotary inertias ``M[(i,r), (i,r)]``; **Master Node Mass** -- load nodes
   (ANSYS numbering, with their own masses and coordinates) slaved to master nodes of the SSI model, whose
   rigid-body motion ``a_k = a_M + alpha_M x (x_k - x_M)`` drives them (a coarse stick -> a refined ANSYS
   model).  With "Generate Mass Data" the mass file is written (overwritten) from the HOUSE matrices,
   otherwise it is read;
4. histories (D-LGN-06): ``source`` RESULTS -- MOTION ``.ACC`` (absolute acceleration, g) and RELDISP
   ``.THD`` (relative displacement) per node and DOF, or, when a file is missing, the frames
   ``ACC/ACCR/THD/THDR``, and for a displacement MOTION's complex ``.TFI`` with RELDISP's free-field formula;
   ``source`` FILE8 -- the same histories computed from FILE8 and the control motion with the MOTION deck
   settings (the MOTION and RELDISP algorithms, so the results equal their files); every relative
   displacement is started at rest (its initial value, the zero-mean offset of the periodic FFT solution,
   is subtracted; deck ``rest``);
5. static (D-LGN-02, D-LGN-07): the critical times (the largest local maxima of the base shear / moment /
   a node response, or given times), then one ANSYS load step per critical time: ``F`` (inertia forces)
   and ``D`` (interface displacements, or 0 for the fixed base) -- one APDL file, or one file per
   critical time with "Use Multiple File List Inputs";
6. dynamic (D-LGN-04): ``*DIM`` TABLE arrays (time column filled with ``*VFILL ... RAMP``) of the
   reference acceleration (``ACEL``) and of the interface relative displacements (``D``), Rayleigh
   damping (``ALPHAD``, ``BETAD``) and the full transient solution commands;
7. writes the listing ``<model>_LOADGEN.out`` and, for static runs, the resultant histories
   ``<apdl>_res.txt`` (t, VX, VY, VZ, MX, MY, MZ of the inertia forces).

The APDL conventions (D-LGN-08): the model's consistent units (accelerations of the SSI files in g are
multiplied by the HOUSE gravity), global Cartesian axes of the model, values with ``digits`` (12)
significant digits, at most 10 values per array assignment, names ``LG...`` of at most 32 characters, a
``LG_SOLVE`` switch (0 = define the loads only) and a header that says what every block is and how to run
it.  Nothing of the user's ANSYS model (elements, materials, supports) is redefined, except the
interface supports of the fixed-base variants (``D = 0``).

Module entry: :func:`run` (``RUNLOADGEN``; LOADGEN is listed in :data:`sassi.modules.base.MODULE_NAMES`, so
:func:`sassi.modules.base.run_module` runs it; :func:`run_loadgen` is a convenience wrapper) and the batch protocol
``python -m sassi.modules.loadgen < LOADGEN.inp`` (three lines: model, deck, listing).
"""
from __future__ import annotations

import math
import re
import sys
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, TextIO, Tuple, Union

import numpy as np

from .. import PRODUCT, __version__
from .. import conventions as C
from ..core import interp as I
from ..core import signal as S
from ..io import decks, thfile
from ..io.container import read_container
from ..io.deckfmt import Table, read_raw, write_raw
from .base import Listing, ModuleContext, ModuleError

NAME = "LOADGEN"
DECK_EXT = ".lgn"
PathLike = Union[str, Path]

#: "Data to Add From ACS SASSI" radio of the static dialog (spec 04 section 15.6 (a))
DATA_NAMES = {1: "Displacement", 2: "Acceleration", 3: "Disp. and Accel.", 4: "Disp. for Soil Module"}
DATA_WORDS = {"DISP": 1, "DISPLACEMENT": 1, "D": 1, "ACC": 2, "ACCEL": 2, "ACCELERATION": 2, "A": 2,
              "DISPACC": 3, "BOTH": 3, "DA": 3, "SOILDISP": 4, "SOIL": 4}
MASS_NAMES = {1: "Lumped Mass", 2: "Master Node Mass"}
MASS_WORDS = {"LUMPED": 1, "LUMP": 1, "MASTER": 2}
SOURCES = ("RESULTS", "FILE8")
METHODS = ("REL", "ACC")
MAPMODES = ("IDENTITY", "PAIRS", "COORD")
#: critical-time criteria (D-LGN-07)
CRITERIA = ("V", "VX", "VY", "VZ", "MX", "MY", "MZ", "ACC", "DISP", "TIME", "STEP")
RESULTANT_NAMES = ("VX", "VY", "VZ", "MX", "MY", "MZ")

DISP_LABELS = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ")
FORCE_LABELS = ("FX", "FY", "FZ", "MX", "MY", "MZ")
#: SASSI X-Z plane -> ANSYS X-Y plane (inverse of the CONVERT,ANSYS rotation, sassi.io.apdl)
R_2D = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]])
VALUES_PER_LINE = 10          # APDL array assignment: at most 10 values per command (*SET)
APDL_MAX_LINE = 640           # ANSYS command-line length limit
APDL_NAME_MAX = 32
TABLE_WARN = 4000             # many TABLE parameters: ANSYS limits the number of parameters
DT_RTOL = 1e-6


def _quiet():
    """Silence the spurious floating-point flags of matmul with numpy 2.0 + macOS Accelerate (see
    :func:`sassi.elements.base.blas_quiet`); results are checked for finiteness where they are written."""
    return np.errstate(divide="ignore", over="ignore", invalid="ignore", under="ignore")


# =======================================================================================
# LOADGEN deck (D-LGN-01): SASSI-EDU deck format with its own schema
# =======================================================================================
@dataclass(frozen=True)
class Param:
    """One deck parameter: name, type (int, float, str or list of floats), default, meaning."""
    name: str
    type: type
    default: Any
    doc: str = ""


PARAMS: Tuple[Param, ...] = (
    Param("title", str, "", "model title"),
    Param("model", str, "", "model name"),
    Param("opmode", int, 0, "0 complete solution, 1 data check (no APDL written)"),
    Param("analysis", str, "STATIC", "STATIC (ANSYS Eq. Static Load) or DYNAMIC (ANSYS Dynamic Load)"),
    # ---- ANSYS Static Load Converter
    Param("data", int, 2, "static data to add: 1 Displacement, 2 Acceleration, 3 Disp. and Accel., "
                          "4 Disp. for Soil Module"),
    Param("multi", int, 0, "Use Multiple File List Inputs: one APDL file per critical time (0: one file, "
                           "the first critical time)"),
    Param("rotdisp", int, 0, "Rotational Disp.: prescribe the rotations too"),
    Param("rotacc", int, 0, "Rotational Accel.: rotary inertia moments / rotational accelerations too"),
    Param("masstype", int, 1, "1 Lumped Mass, 2 Master Node Mass"),
    Param("genmass", int, 0, "Generate Mass Data: (over)write the mass file from the HOUSE mass matrix"),
    Param("source", str, "RESULTS", "RESULTS (MOTION .ACC / RELDISP .THD files or frames) or FILE8"),
    # ---- ANSYS Dynamic Load Converter
    Param("method", str, "REL", "REL: ACEL reference motion + D relative displacements (manual); ACC: fixed "
                                "interface, ACEL = absolute acceleration of <refnode>"),
    Param("alpha", float, 0.0, "Rayleigh damping alpha (ALPHAD), C = alpha M + beta K"),
    Param("beta", float, 0.0, "Rayleigh damping beta (BETAD)"),
    Param("refnode", int, 0, "reference node of the ACEL motion (model numbering); 0 = control motion"),
    # ---- files and paths
    Param("ssipath", str, "", "SASSI Model and Results Input path ('' = model directory)"),
    Param("housefile", str, "", "HOUSE Module Input ('' = <model>.hou)"),
    Param("dispdir", str, "THD", "translational displacement frames (RELDISP)"),
    Param("disprotdir", str, "THDR", "rotational displacement frames (RELDISP)"),
    Param("accdir", str, "ACC", "translational acceleration frames (MOTION)"),
    Param("accrotdir", str, "ACCR", "rotational acceleration frames (MOTION)"),
    Param("ansyspath", str, "", "ANSYS Model and Data Input path ('' = model directory)"),
    Param("lumpfile", str, "", "lumped-mass file ('' = <model>.masl)"),
    Param("masterfile", str, "", "master-node mass file ('' = <model>.masm)"),
    Param("apdlfile", str, "", "APDL file ('' = <model>_LGS.inp static, <model>_LGD.inp dynamic)"),
    Param("groundfile", str, "", "Ground Acceleration File, g ('' = the control motion of the decks)"),
    Param("groundfopt", int, 0, "format of the ground file: 0 dt then values, 1 (t, a) pairs"),
    Param("groundmult", float, 1.0, "scale factor of the ground file"),
    # ---- node mapping SASSI-EDU -> ANSYS
    Param("mapmode", str, "IDENTITY", "IDENTITY, PAIRS (file of 'sassi ansys' pairs) or COORD (ANSYS file)"),
    Param("mapfile", str, "", "node map file (PAIRS) or ANSYS .cdb / APDL file with the nodes (COORD)"),
    Param("maptol", float, 0.0, "COORD matching tolerance, length (0 = 1e-6 x model size)"),
    # ---- critical times (static)
    Param("crit", str, "V", "criterion: V VX VY VZ MX MY MZ ACC DISP TIME STEP"),
    Param("ncrit", int, 1, "number of critical times (largest local maxima)"),
    Param("tsep", float, 0.0, "minimum separation of the critical times, s"),
    Param("critnode", int, 0, "node of the ACC / DISP criterion (model numbering)"),
    Param("critdof", int, 1, "DOF 1..6 of the ACC / DISP criterion"),
    Param("refpoint", list, [], "moment reference point x y z ([] = centroid of the interface nodes)"),
    Param("times", list, [], "critical times (TIME, s) or time steps (STEP, 1-based)"),
    Param("digits", int, 12, "significant digits of the APDL values (6..17)"),
    Param("rest", int, 1, "1: relative displacements start at rest (their initial value, the zero-mean offset of "
                          "the periodic SSI solution, is subtracted); 0: as computed"),
)
TABLES: Dict[str, Tuple[str, ...]] = {
    "dnodes": ("node",),        # nodes that receive D (interface; default: interaction nodes)
    "masters": ("node",),       # master nodes (Master Node Mass generation)
    "checknodes": ("node",),    # dynamic: nodes whose absolute accelerations are exported as tables
}
_PARAM = {p.name: p for p in PARAMS}


@dataclass
class LoadgenDeck:
    """Typed LOADGEN deck: ``params`` (name -> value) and node ``tables`` (name -> list of node ids)."""
    params: Dict[str, Any] = field(default_factory=lambda: {p.name: (list(p.default) if isinstance(p.default, list)
                                                                      else p.default) for p in PARAMS})
    tables: Dict[str, List[int]] = field(default_factory=lambda: {k: [] for k in TABLES})

    def __getitem__(self, name: str) -> Any:
        return self.params[name]

    def __setitem__(self, name: str, value: Any) -> None:
        if name not in _PARAM:
            raise KeyError(f"LOADGEN deck: unknown parameter {name!r}")
        self.params[name] = value


def new_deck() -> LoadgenDeck:
    """A LOADGEN deck with every parameter at its default."""
    return LoadgenDeck()


def deck_path(model_dir: PathLike, model: str) -> Path:
    return Path(model_dir) / f"{model}{DECK_EXT}"


def _coerce(v: Any, typ: type) -> Any:
    if typ is list:
        if isinstance(v, str):
            t = v.strip()
            if t.startswith("["):
                import json
                v = json.loads(t)
            else:
                v = [x for x in t.replace(",", " ").split() if x]
        return [float(str(x).replace("D", "E").replace("d", "e")) if isinstance(x, str) else float(x) for x in v]
    if typ is str:
        if isinstance(v, str):
            t = v.strip()
            if t.startswith('"'):
                import json
                return json.loads(t)
            return t
        return "" if v is None else str(v)
    if typ is int:
        f = float(str(v).strip().replace("D", "E").replace("d", "e") or 0) if isinstance(v, str) else float(v)
        if f != round(f):
            raise ValueError(f"expected an integer, got {v!r}")
        return int(round(f))
    if typ is float:
        return float(str(v).strip().replace("D", "E").replace("d", "e") or 0) if isinstance(v, str) else float(v)
    raise TypeError(typ)


def write_deck(path: PathLike, deck: LoadgenDeck, comments: Sequence[str] = ()) -> Path:
    """Write the LOADGEN deck (``SASSI-EDU LOADGEN DECK v1``, the format of :mod:`sassi.io.deckfmt`)."""
    unknown = set(deck.params) - set(_PARAM)
    if unknown:
        raise KeyError(f"LOADGEN deck: unknown parameters {sorted(unknown)}")
    params = {p.name: _coerce(deck.params.get(p.name, p.default), p.type) for p in PARAMS}
    tables = {}
    for name, cols in TABLES.items():
        t = Table(columns=list(cols))
        for n in deck.tables.get(name, []):
            t.rows.append([int(n)])
        tables[name] = t
    return write_raw(path, NAME, params, tables, version=1, comments=comments)


def read_deck(path: PathLike) -> Tuple[LoadgenDeck, List[str]]:
    """Read a LOADGEN deck; returns ``(deck, notes)`` (missing parameters take their defaults)."""
    raw = read_raw(path)
    if raw.module != NAME:
        raise ValueError(f"{Path(path).name}: expected a LOADGEN deck, found {raw.module}")
    d = new_deck()
    notes: List[str] = []
    for p in PARAMS:
        if p.name in raw.params:
            try:
                d.params[p.name] = _coerce(raw.params[p.name], p.type)
            except (ValueError, TypeError) as exc:
                raise ValueError(f"{Path(path).name}: parameter {p.name}: {exc}") from None
    for k in raw.params:
        if k not in _PARAM:
            notes.append(f"unknown deck parameter '{k}' ignored")
    for name, cols in TABLES.items():
        t = raw.tables.get(name)
        if t is None:
            continue
        if list(t.columns) != list(cols):
            raise ValueError(f"{Path(path).name}: table {name} has columns {t.columns}, expected {list(cols)}")
        d.tables[name] = [int(_coerce(r[0], int)) for r in t.rows]
    return d, notes


# =======================================================================================
# Model information from FILE4 (and the HOUSE deck)
# =======================================================================================
@dataclass
class ModelInfo:
    """Nodes, equations and numbering of the SSI model (FILE4 = ``<model>.N4``)."""
    model: str
    gravity: float
    dim: int
    two_d: bool                                   # PLANE elements: the ANSYS export rotates X-Z -> X-Y
    xyz: Dict[int, np.ndarray]                    # FILE4 node -> global coordinates
    model_of: Dict[int, int]                      # FILE4 node -> model (original) node
    file4_of: Dict[int, int]                      # model node -> FILE4 node
    dofs: Dict[int, List[int]]                    # FILE4 node -> DOFs with an equation
    eq_node: np.ndarray
    eq_dof: np.ndarray
    int_nodes: List[int]                          # FILE4 numbering, interaction order
    excav_nodes: List[int]                        # FILE4 numbering
    optimized: bool
    meta: dict
    title: str = ""
    house_masses: Dict[int, np.ndarray] = field(default_factory=dict)   # model node -> (6,) mass units

    @property
    def size(self) -> float:
        """Diagonal of the bounding box of the nodes (for default tolerances)."""
        if not self.xyz:
            return 1.0
        P = np.array(list(self.xyz.values()))
        return float(np.linalg.norm(P.max(axis=0) - P.min(axis=0))) or 1.0


def read_model_info(workdir: Path, model: str, house_path: Optional[Path] = None) -> ModelInfo:
    """FILE4 (required) and the HOUSE deck (title, PLANE groups and masses; optional)."""
    p4 = workdir / f"{model}.N4"
    if not p4.exists():
        raise ModuleError(f"{p4.name} (FILE4) missing -- run HOUSE first")
    try:
        f4 = read_container(p4, kind="FILE4")
    except (OSError, ValueError) as exc:
        raise ModuleError(f"{p4.name}: {exc}") from None
    ids = np.asarray(f4["node_id"], dtype=np.int64)
    xyz = np.asarray(f4["node_xyz"], dtype=float).reshape(-1, 3)
    old = np.asarray(f4.get("x_node_old_id", ids), dtype=np.int64)
    eq_node = np.asarray(f4["eq_node"], dtype=np.int64)
    eq_dof = np.asarray(f4["eq_dof"], dtype=np.int64)
    dofs: Dict[int, List[int]] = {}
    for n, k in zip(eq_node, eq_dof):
        dofs.setdefault(int(n), []).append(int(k))
    for n in dofs:
        dofs[n].sort()
    etype = np.asarray(f4["elem_type"], dtype=np.int64)
    excav = np.asarray(f4["elem_excav"], dtype=np.int64)
    enodes = np.asarray(f4["elem_nodes"], dtype=np.int64).reshape(len(etype), -1)
    excav_nodes = sorted({int(n) for row in enodes[excav == 1] for n in row if n > 0})
    meta = dict(f4.meta)
    info = ModelInfo(model=model, gravity=float(meta.get("gravity", 0.0)), dim=int(meta.get("dim", 2)),
                     two_d=bool(np.any(etype == 4)), xyz={int(n): p for n, p in zip(ids, xyz)},
                     model_of={int(n): int(o) for n, o in zip(ids, old)},
                     file4_of={int(o): int(n) for n, o in zip(ids, old)}, dofs=dofs, eq_node=eq_node,
                     eq_dof=eq_dof, int_nodes=[int(n) for n in f4["int_node"]], excav_nodes=excav_nodes,
                     optimized=bool(int(meta.get("x_optimized", 0) or 0)), meta=meta,
                     title=str(meta.get("title", "")))
    if house_path is not None and house_path.exists():
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                hd = decks.read(house_path, "HOUSE")
            except (OSError, ValueError, KeyError):
                hd = None
        if hd is not None:
            info.title = info.title or str(hd["title"])
            g = float(hd["gravity"]) or info.gravity
            types = {int(r["id"]): int(r["type"]) for r in hd.rows("groups")}
            used = {int(r["group"]) for r in hd.rows("elements")}
            info.two_d = info.two_d or any(types.get(gid) == 4 for gid in used)
            for r in hd.rows("masses"):
                v = np.array([r["mx"], r["my"], r["mz"], r["mxx"], r["myy"], r["mzz"]], dtype=float)
                if int(r["units"]) == 1 and g > 0:
                    v = v / g
                info.house_masses[int(r["node"])] = info.house_masses.get(int(r["node"]), 0.0) + v
    return info


def ansys_dof(dof: int, two_d: bool) -> Tuple[int, float]:
    """ANSYS DOF index (0..5) and sign of SASSI DOF ``dof`` (1..6) under the ANSYS export rotation.

    3-D models keep the axes.  2-D models (PLANE groups) are rotated from the SASSI X-Z plane to the
    ANSYS X-Y plane, ``x_A = R x_S`` with ``R = [[1,0,0],[0,0,1],[0,-1,0]]`` (sassi.io.apdl): UX -> UX,
    UY -> -UZ, UZ -> UY and, rotations being vectors of the same proper rotation, ROTY -> -ROTZ,
    ROTZ -> ROTY."""
    s = int(dof) - 1
    if not two_d:
        return s, 1.0
    col = R_2D[:, s % 3]
    a = int(np.argmax(np.abs(col)))
    return a + (3 if s >= 3 else 0), float(col[a])


def control_direction(cm: int, ang: float) -> np.ndarray:
    """Global unit vector of the control motion: SITE <cm> (0 x', 1 y', 2 z') rotated by ANALYS <ang>."""
    return np.array([I.rigid_body_anchor(k, int(cm), float(ang)) for k in (1, 2, 3)], dtype=float)


# =======================================================================================
# Node mapping SASSI-EDU -> ANSYS (D-LGN-05)
# =======================================================================================
@dataclass
class NodeMap:
    """Model node -> ANSYS node.  ``table`` is None for the identity."""
    mode: str = "IDENTITY"
    table: Optional[Dict[int, int]] = None
    source: str = ""
    notes: List[str] = field(default_factory=list)

    def ansys(self, model_node: int) -> Optional[int]:
        if self.table is None:
            return int(model_node)
        return self.table.get(int(model_node))


def read_pairs(path: Path) -> Dict[int, int]:
    """``sassi ansys`` node pairs, one per line; comments after ``!``, ``#`` or a leading ``*``."""
    out: Dict[int, int] = {}
    for k, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
        s = re.split(r"[!#]", line, maxsplit=1)[0].strip()
        if not s or s.startswith("*"):
            continue
        parts = s.replace(",", " ").split()
        if len(parts) < 2:
            raise ModuleError(f"{path.name} line {k}: expected 'sassi_node ansys_node', got {line.strip()!r}")
        try:
            a, b = int(float(parts[0])), int(float(parts[1]))
        except ValueError:
            raise ModuleError(f"{path.name} line {k}: node numbers expected, got {line.strip()!r}") from None
        if a in out and out[a] != b:
            raise ModuleError(f"{path.name} line {k}: SASSI node {a} mapped twice ({out[a]} and {b})")
        out[a] = b
    return out


def coordinate_map(info: ModelInfo, ansys_nodes: Dict[int, np.ndarray], tol: float) -> Tuple[Dict[int, int], int]:
    """Model node -> ANSYS node at the same position (after the 2-D rotation), within ``tol``.

    Returns the map and the number of SASSI nodes without a match.  Nodes of a refined ANSYS mesh that
    lie between SASSI nodes are not interpolated: give them their own load data (pairs file / master
    mass file) or keep the interface mesh of the SSI model."""
    if not ansys_nodes:
        return {}, len(info.xyz)
    from scipy.spatial import cKDTree
    aid = np.array(sorted(ansys_nodes), dtype=np.int64)
    A = np.array([ansys_nodes[int(n)] for n in aid], dtype=float).reshape(-1, 3)
    tree = cKDTree(A)
    sid = np.array(sorted(info.xyz), dtype=np.int64)
    P = np.array([info.xyz[int(n)] for n in sid], dtype=float).reshape(-1, 3)
    if info.two_d:
        with _quiet():
            P = P @ R_2D.T
    dist, idx = tree.query(P, k=1)
    out: Dict[int, int] = {}
    missing = 0
    for n, d_, j in zip(sid, dist, idx):
        if d_ <= tol:
            out[info.model_of[int(n)]] = int(aid[int(j)])
        else:
            missing += 1
    return out, missing


def build_node_map(d: LoadgenDeck, info: ModelInfo, ansys_dir: Path, workdir: Path, lst: Listing) -> NodeMap:
    mode = str(d["mapmode"] or "IDENTITY").strip().upper()
    if mode not in MAPMODES:
        raise ModuleError(f"node map mode {mode!r} must be one of {', '.join(MAPMODES)}")
    if mode == "IDENTITY":
        return NodeMap("IDENTITY", None, "ANSYS node = SASSI-EDU model node (ANSYS export / CONVERT,ANSYS numbering)")
    name = str(d["mapfile"] or "").strip()
    if not name:
        raise ModuleError(f"node map mode {mode} needs a file (LGMAP,{mode},...,<file>)")
    p = _resolve(name, workdir)
    if not p.exists():
        raise ModuleError(f"node map file {p} not found")
    if mode == "PAIRS":
        return NodeMap("PAIRS", read_pairs(p), f"pairs file {p.name}")
    from ..io.ansys_cdb import read_cdb
    try:
        cdb = read_cdb(p)
    except OSError as exc:
        raise ModuleError(f"cannot read {p}: {exc}") from None
    if not cdb.nodes:
        raise ModuleError(f"{p.name}: no ANSYS nodes found (NBLOCK or N commands)")
    tol = float(d["maptol"]) if float(d["maptol"]) > 0 else 1e-6 * info.size
    table, missing = coordinate_map(info, cdb.nodes, tol)
    nm = NodeMap("COORD", table, f"coordinates of {p.name} ({len(cdb.nodes)} ANSYS nodes, tolerance {tol:.3g})")
    if cdb.node_rot:
        nm.notes.append(f"{len(cdb.node_rot)} ANSYS nodes have rotated nodal coordinate systems: LOADGEN writes "
                        "global components (F and D act in the nodal systems in ANSYS)")
    nm.notes.append(f"{len(table)} SASSI nodes matched, {missing} without an ANSYS node at their position")
    if not int(d["opmode"]):
        out = ansys_dir / f"{info.model}_LG.map"
        lines = [f"! SASSI-EDU LOADGEN node map (model {info.model}): sassi_node ansys_node, matched by coordinates "
                 f"of {p.name} within {tol:.6g}"]
        lines += [f"{a} {b}" for a, b in sorted(table.items())]
        out.write_text("\n".join(lines) + "\n", encoding="utf-8")
        nm.notes.append(f"map written to {out.name} (usable as a PAIRS file)")
    return nm


# =======================================================================================
# Masses (D-LGN-03)
# =======================================================================================
@dataclass
class LoadMass:
    """Inertia of one load node: ``m`` = (mx, my, mz, mxx, myy, mzz) in mass units.

    Lumped Mass: ``node`` is a model node (its own acceleration).  Master Node Mass: ``node`` is an
    ANSYS node at ``xyz`` (SASSI-EDU global coordinates) driven by the rigid-body motion of model node
    ``master``."""
    node: int
    m: np.ndarray
    master: int = 0
    xyz: Optional[np.ndarray] = None


def lumped_from_matrix(info: ModelInfo, Ms) -> Dict[int, np.ndarray]:
    """Lumped nodal masses from the HOUSE structure mass matrix (COOSM ``Ms``, FILE4 equations).

    Translations: ``m_(i,d) = sum_j M[(i,d), (j,d)]`` over the translational equations j of the same
    direction d (row sums; for a lumped matrix the diagonal, for a consistent one the tributary mass, and
    the total ``iota_d^T M iota_d`` per direction is preserved).  Rotations: the diagonal rotary inertia
    ``M[(i,r), (i,r)]``.  Returns {model node: (6,)} for the nodes with a non-zero mass."""
    import scipy.sparse as sp
    M = sp.csr_matrix(Ms)
    neq = M.shape[0]
    if neq != len(info.eq_dof):
        raise ModuleError(f"COOSM has {neq} equations, FILE4 {len(info.eq_dof)}: run HOUSE again")
    lump = np.zeros(neq)
    for d in (1, 2, 3):
        mask = (info.eq_dof == d).astype(float)
        rs = M @ mask
        lump[info.eq_dof == d] = rs[info.eq_dof == d]
    diag = M.diagonal()
    rot = info.eq_dof >= 4
    lump[rot] = diag[rot]
    out: Dict[int, np.ndarray] = {}
    for i in np.flatnonzero(lump != 0.0):
        n = info.model_of[int(info.eq_node[i])]
        out.setdefault(n, np.zeros(6))[int(info.eq_dof[i]) - 1] += lump[i]
    return out


def _mass_header(kind: str, info: ModelInfo, how: str) -> List[str]:
    return [f"! SASSI-EDU {__version__} LOADGEN {kind} file (Option A), model {info.model}",
            f"! {how}",
            f"! units: mass units of the model (weight / g, g = {info.gravity:g}); rotary inertias in mass x length^2"]


def write_lumped(path: Path, masses: Dict[int, np.ndarray], info: ModelInfo, how: str) -> Path:
    """Lumped-mass file: ``node mx my mz mxx myy mzz`` (model numbering), one node per line."""
    lines = _mass_header("lumped-mass", info, how)
    lines.append(f"! {'node':>8s} {'mx':>22s} {'my':>22s} {'mz':>22s} {'mxx':>22s} {'myy':>22s} {'mzz':>22s}")
    for n in sorted(masses):
        lines.append(f"{n:10d} " + " ".join(f"{v:22.15e}" for v in masses[n]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _numbers_of(line: str) -> Optional[List[float]]:
    s = re.split(r"[!#]", line, maxsplit=1)[0].strip()
    if not s or s.startswith("*"):
        return None
    try:
        return [float(t.replace("D", "E").replace("d", "e")) for t in s.replace(",", " ").split()]
    except ValueError:
        return []


def read_lumped(path: Path) -> Dict[int, np.ndarray]:
    """Read a lumped-mass file: lines ``node mx my mz [mxx myy mzz]`` (a single value = mx = my = mz)."""
    out: Dict[int, np.ndarray] = {}
    for k, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
        v = _numbers_of(line)
        if v is None:
            continue
        if len(v) not in (2, 4, 7):
            raise ModuleError(f"{path.name} line {k}: expected 'node mx my mz [mxx myy mzz]', got {line.strip()!r}")
        m = np.zeros(6)
        if len(v) == 2:
            m[:3] = v[1]
        else:
            m[:len(v) - 1] = v[1:]
        if np.any(m < 0):
            raise ModuleError(f"{path.name} line {k}: negative mass")
        out[int(v[0])] = out.get(int(v[0]), 0.0) + m
    return out


def write_master(path: Path, masses: List[LoadMass], info: ModelInfo, how: str) -> Path:
    """Master-node mass file: ``node master x y z mx my mz mxx myy mzz`` (node = ANSYS load node,
    master = model node of the SSI model, x y z = load-node position in SASSI-EDU global axes)."""
    lines = _mass_header("master-node mass", info, how)
    lines.append("! load node (ANSYS numbering) driven by the rigid-body motion of its master node (SASSI-EDU model "
                 "numbering): a = a_M + alpha_M x (x - x_M)")
    lines.append(f"! {'node':>8s} {'master':>10s} {'x':>16s} {'y':>16s} {'z':>16s} {'mx':>22s} {'my':>22s} "
                 f"{'mz':>22s} {'mxx':>22s} {'myy':>22s} {'mzz':>22s}")
    for lm in masses:
        x = lm.xyz if lm.xyz is not None else np.zeros(3)
        lines.append(f"{lm.node:10d} {lm.master:10d} " + " ".join(f"{v:16.9e}" for v in x) + " "
                     + " ".join(f"{v:22.15e}" for v in lm.m))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def read_master(path: Path) -> List[LoadMass]:
    """Read a master-node mass file (lines ``node master x y z mx my mz [mxx myy mzz]``)."""
    out: List[LoadMass] = []
    for k, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
        v = _numbers_of(line)
        if v is None:
            continue
        if len(v) not in (8, 11):
            raise ModuleError(f"{path.name} line {k}: expected 'node master x y z mx my mz [mxx myy mzz]', got "
                              f"{line.strip()!r}")
        m = np.zeros(6)
        m[:len(v) - 5] = v[5:]
        if np.any(m < 0):
            raise ModuleError(f"{path.name} line {k}: negative mass")
        out.append(LoadMass(int(v[0]), m, master=int(v[1]), xyz=np.array(v[2:5], dtype=float)))
    return out


def master_assignment(info: ModelInfo, nodes: Iterable[int], masters: Sequence[int]) -> Dict[int, int]:
    """Master of each model node: the master closest in elevation (z), ties broken by the plan
    distance -- the floor master of a building model (D-LGN-03).  Model numbering in and out."""
    mxyz = []
    for mn in masters:
        f4 = info.file4_of.get(int(mn))
        if f4 is None:
            raise ModuleError(f"master node {mn} is not a node of the SSI model (FILE4)")
        mxyz.append(info.xyz[f4])
    mxyz = np.array(mxyz).reshape(-1, 3)
    out: Dict[int, int] = {}
    tol = 1e-6 * info.size
    for n in nodes:
        p = info.xyz[info.file4_of[int(n)]]
        dz = np.abs(mxyz[:, 2] - p[2])
        best = np.flatnonzero(dz <= dz.min() + tol)
        dxy = np.hypot(mxyz[best, 0] - p[0], mxyz[best, 1] - p[1])
        out[int(n)] = int(masters[int(best[int(np.argmin(dxy))])])
    return out


# =======================================================================================
# Histories (D-LGN-06)
# =======================================================================================
_FRAME_RE = re.compile(r"^(?P<tag>[A-Za-z]+)_(?P<val>-?\d+(?:\.\d*)?)_(?P<num>\d+)$")


def find_node_file(directory: Path, node: int, dof: int, ext: str) -> Optional[Path]:
    """MOTION/RELDISP file of (node, dof) with 5- or 6-digit node numbers (D-FIL-05)."""
    for nm in dict.fromkeys([C.nodal_result_name(node, dof, ext, 0), C.nodal_result_name(node, dof, ext, 100000)]):
        p = directory / nm
        if p.exists():
            return p
    return None


def read_history_file(path: Path) -> Tuple[np.ndarray, float]:
    """``.ACC``/``.THD`` history: ``#`` header lines, dt, one value per line (D-FIL-02)."""
    vals: List[float] = []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s[0] in "#*":
            continue
        for t in s.replace(",", " ").split():
            vals.append(float(t.replace("D", "E")))
    if not vals:
        raise ModuleError(f"{path.name}: empty history file")
    return np.asarray(vals[1:], dtype=float), float(vals[0])


def read_frames(directory: Path) -> Tuple[np.ndarray, Dict[int, int], np.ndarray]:
    """Frames ``<TAG>_<t>_<nnnnn>`` of a MOTION/RELDISP frame directory: (steps, {node: row},
    values (n_steps, n_nodes, 3)), steps in frame-number order (D-FIL-03)."""
    files = []
    for p in directory.iterdir():
        m = _FRAME_RE.match(p.name)
        if m and p.is_file():
            files.append((int(m.group("num")), p))
    if not files:
        raise ModuleError(f"no frame files in {directory}")
    files.sort()
    steps = np.array([k for k, _ in files], dtype=np.int64)
    if np.any(np.diff(steps) != 1) or steps[0] != 1:
        raise ModuleError(f"{directory.name}: frame numbers are not 1, 2, 3 ... (missing frames)")
    nodes: Optional[np.ndarray] = None
    vals = None
    for s, (k, p) in enumerate(files):
        rows = []
        lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
        for ln in lines[1:]:
            t = ln.split()
            if t:
                rows.append([float(x.replace("D", "E")) for x in t])
        A = np.asarray(rows, dtype=float)
        if nodes is None:
            nodes = A[:, 0].astype(np.int64)
            vals = np.zeros((len(files), len(nodes), A.shape[1] - 1))
        elif not np.array_equal(A[:, 0].astype(np.int64), nodes):
            raise ModuleError(f"{p.name}: the frame lists other nodes than the first frame of {directory.name}")
        vals[s] = A[:, 1:]
    return steps, {int(n): i for i, n in enumerate(nodes)}, vals


@dataclass
class HistorySet:
    """Histories in model units on a common time grid: accelerations (absolute; length/s^2, rad/s^2)
    and displacements (relative to the reference motion; length, rad), keyed by (FILE4 node, dof)."""
    dt: float
    n: int
    acc: Dict[Tuple[int, int], np.ndarray] = field(default_factory=dict)
    disp: Dict[Tuple[int, int], np.ndarray] = field(default_factory=dict)
    origin: Dict[str, str] = field(default_factory=dict)
    ref: Dict[int, np.ndarray] = field(default_factory=dict)      # source FILE8: reference-node acceleration per dof


def _deck_or_none(path: Path, module: str):
    import warnings
    if not path.exists():
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return decks.read(path, module)
        except (OSError, ValueError, KeyError):
            return None


class ResultsSource:
    """MOTION ``.ACC`` / RELDISP ``.THD`` files of the SSI results directory, their frames, or -- for the
    relative displacements -- MOTION's complex ``.TFI`` with RELDISP's free-field formula."""

    def __init__(self, d: LoadgenDeck, ssidir: Path, gravity: float, delt_hint: Optional[float], model: str = ""):
        self.d = d
        self.dir = ssidir
        self.g = gravity
        self.delt_hint = delt_hint
        self.model = model
        self._frames: Dict[str, Tuple[np.ndarray, Dict[int, int], np.ndarray]] = {}
        self._tfi: Optional[dict] = None
        self.used: Dict[str, int] = {}

    def _tfi_setup(self) -> Optional[dict]:
        """Control motion and grid of the MOTION deck for the ``.TFI`` fallback (None without a deck)."""
        if self._tfi is None:
            md = _deck_or_none(self.dir / f"{self.model}.mot", "MOTION")
            if md is None or not str(md["thfile"] or "").strip():
                self._tfi = {}
            else:
                nfft, dt = int(md["nft"]), float(md["delt"])
                f = S.fourier_grid(nfft, dt)
                _, A, name = control_motion_from_deck(md, self.dir)
                cm, ang = int(md["cm"]), float(md["ang"])
                p8 = self.dir / str(md["file8"] or "FILE8")
                if p8.exists():
                    try:
                        meta = read_container(p8).meta
                        cm, ang = int(meta.get("cm", cm)), float(meta.get("ang", ang))
                    except (OSError, ValueError):
                        pass
                self._tfi = dict(nfft=nfft, dt=dt, df=1.0 / (nfft * dt), nK=len(f), cm=cm, ang=ang,
                                 drive=S.displacement_spectrum(A, f, self.g),
                                 nout=S.output_length(nfft, dt, float(md["dur"])))
        return self._tfi or None

    def _from_tfi(self, node: int, dof: int) -> Tuple[Optional[np.ndarray], Optional[float]]:
        """Relative displacement (free-field reference) from the complex ``.TFI``: RELDISP's
        ``irfft((H - H_ff) U_g)`` (D-RDP-01/04), the same numbers as its ``.THD``."""
        p = find_node_file(self.dir, node, dof, "TFI")
        st = self._tfi_setup() if p is not None else None
        if st is None:
            return None, None
        from .reldisp import read_tf_on_grid
        H = read_tf_on_grid(p, st["nK"], st["df"])
        D = (H - I.rigid_body_anchor(dof, st["cm"], st["ang"])) * st["drive"]
        u = np.fft.irfft(S.hermitian_bins(D, st["nfft"]), n=st["nfft"])
        return u[:st["nout"]], st["dt"]

    def _frame(self, kind: str, dof: int) -> Optional[Tuple[np.ndarray, Dict[int, int], np.ndarray]]:
        key = {("acc", False): "accdir", ("acc", True): "accrotdir", ("disp", False): "dispdir",
               ("disp", True): "disprotdir"}[(kind, dof >= 4)]
        name = str(self.d[key] or "").strip()
        if not name:
            return None
        if name not in self._frames:
            p = _resolve(name, self.dir)
            if not p.is_dir():
                self._frames[name] = None            # type: ignore[assignment]
            else:
                self._frames[name] = read_frames(p)
        return self._frames[name]

    def get(self, kind: str, node: int, dof: int) -> Tuple[Optional[np.ndarray], Optional[float], str]:
        """(history in file units, dt or None for frames, origin) or (None, None, '')."""
        ext = "ACC" if kind == "acc" else "THD"
        p = find_node_file(self.dir, node, dof, ext)
        if p is not None:
            v, dt = read_history_file(p)
            self.used[ext] = self.used.get(ext, 0) + 1
            return v, dt, p.name
        fr = self._frame(kind, dof)
        if fr is not None:
            steps, rows, vals = fr
            r = rows.get(int(node))
            if r is not None:
                self.used["frames"] = self.used.get("frames", 0) + 1
                return vals[:, r, (dof - 1) % 3].copy(), None, "frames"
        if kind == "disp":
            try:
                v, dt = self._from_tfi(node, dof)
            except ModuleError:
                raise
            except (OSError, ValueError) as exc:
                raise ModuleError(f"relative displacement of node {node} {C.DOF_TAGS[dof]} from its .TFI: {exc}") from None
            if v is not None:
                self.used["TFI"] = self.used.get("TFI", 0) + 1
                return v, dt, "TFI"
        return None, None, ""


def _resolve(name: str, base: Path) -> Path:
    raw = str(name).strip().strip('"')
    if "\\" in raw and "/" not in raw:
        raw = raw.replace("\\", "/")
    p = Path(raw).expanduser()
    return p if p.is_absolute() else base / p


def load_results(d: LoadgenDeck, info: ModelInfo, ssidir: Path, acc_req: Iterable[Tuple[int, int]],
                 disp_req: Iterable[Tuple[int, int]], delt_hint: Optional[float], lst: Listing) -> HistorySet:
    """Histories of the SSI result files (``source`` RESULTS).  Node numbers of the files are the FILE4
    numbers (MOTION and RELDISP work on FILE8, which uses the HOUSE optimizer numbering)."""
    src = ResultsSource(d, ssidir, info.gravity, delt_hint, info.model)
    acc: Dict[Tuple[int, int], np.ndarray] = {}
    disp: Dict[Tuple[int, int], np.ndarray] = {}
    missing: List[str] = []
    dts: Dict[str, float] = {}
    for kind, req, store in (("acc", acc_req, acc), ("disp", disp_req, disp)):
        for node, dof in sorted(set(req)):
            v, dt, origin = src.get(kind, node, dof)
            if v is None:
                missing.append(f"{node} {C.DOF_TAGS[dof]} ({'.ACC' if kind == 'acc' else '.THD'})")
                continue
            if dt is not None:
                dts[origin] = dt
            store[(node, dof)] = v * (info.gravity if kind == "acc" else 1.0)
    if missing:
        raise ModuleError("SSI result histories missing for " + ", ".join(missing[:30]) + (" ..." if len(missing) > 30
                          else "") + f" in {ssidir}: run MOTION (NOUT flag 2 or Save ACC in All Points, Save Rotation "
                          "for Ansys; complex TFs for the displacements) and RELDISP (RDND, free-field reference) for "
                          "these nodes, or use LOADGEN source FILE8")
    if dts:
        vals = np.array(list(dts.values()))
        dt = float(vals[0])
        bad = [k for k, v in dts.items() if abs(v - dt) > DT_RTOL * dt]
        if bad:
            raise ModuleError(f"the result files have different time steps ({dt:g} s in "
                              f"{next(iter(dts))}, other values in {', '.join(bad[:5])})")
    elif delt_hint:
        dt = float(delt_hint)
    else:
        raise ModuleError("the time step of the frames is unknown (no MOTION/RELDISP deck): use .ACC/.THD files")
    lengths = {len(v) for v in list(acc.values()) + list(disp.values())}
    n = min(lengths) if lengths else 0
    if len(lengths) > 1:
        lst.warning(f"the result histories have different lengths ({min(lengths)} .. {max(lengths)} samples): the "
                    f"first {n} samples are used (MOTION and RELDISP <dur> should be equal)")
    hs = HistorySet(dt=dt, n=n, acc={k: v[:n] for k, v in acc.items()}, disp={k: v[:n] for k, v in disp.items()})
    hs.origin["histories"] = (", ".join(f"{v} {k}" for k, v in sorted(src.used.items())) + " files of MOTION / RELDISP"
                              + (f" in {ssidir.name}" if ssidir.name else "")
                              + ("; TFI = relative displacements computed from MOTION's complex .TFI (RELDISP formula, "
                                 "free field)" if "TFI" in src.used else ""))
    return hs


@dataclass
class MotionSetup:
    """Control motion and interpolation settings of the MOTION deck (``source`` FILE8)."""
    nfft: int
    dt: float
    nout: int
    gravity: float
    option: int
    smooth: float
    pzadj: int
    bl: bool
    cm: int
    ang: float
    A: np.ndarray            # rfft of the scaled, padded control motion (g)
    padded: np.ndarray       # control motion padded to NFFT (g)
    thfile: str


def control_motion_from_deck(dk, workdir: Path) -> Tuple[np.ndarray, np.ndarray, str]:
    """(padded control motion in g, its rfft, file name) of a MOTION / RELDISP / STRESS deck, scaled and
    checked exactly as those modules do (:func:`sassi.modules.motion.read_control_motion`)."""
    from .motion import read_control_motion
    name = str(dk["thfile"] or "").strip()
    if not name:
        raise ModuleError(f"no control-motion file (THFILE) in the {dk.module} deck")
    p = _resolve(name, workdir)
    cmot = read_control_motion(p, int(dk["fopt"]), int(dk["rec1"]), int(dk["rec2"]), float(dk["delt"]),
                               int(dk["nft"]), float(dk["mult"]), float(dk["max"]))
    return cmot.padded, cmot.A, p.name


def file8_histories(info: ModelInfo, workdir: Path, model: str, acc_req: Iterable[Tuple[int, int]],
                    disp_req: Iterable[Tuple[int, int]], refnode_f4: int, lst: Listing,
                    ref_baseline: bool = True) -> Tuple[HistorySet, MotionSetup]:
    """Absolute accelerations and relative displacements from FILE8 and the control motion (``source``
    FILE8), with the MOTION deck settings: ``a = irfft(H~ A)`` as MOTION (interpolation option, smoothing,
    phase adjustment, Nyquist rule, Hudson-Housner baseline when <bl> = 1) and ``u = irfft((H~ - H~_ref)
    U_g)`` with ``U_g = -g A / w^2`` as RELDISP (D-RDP-01, D-CNV-06).  ``H~_ref`` is the free field (the
    rigid-body anchor of the control direction, D-RDP-04) or, for a reference node, its interpolated TF
    (translations; rotations are absolute because the reference frame translates only).  The reference node's
    acceleration is baseline-corrected like MOTION's only with ``ref_baseline`` (method ACC); method REL needs
    the exact second derivative of the motion its relative displacements refer to."""
    from .motion import read_file8
    md = _deck_or_none(workdir / f"{model}.mot", "MOTION")
    if md is None:
        raise ModuleError(f"{model}.mot (MOTION deck) missing: source FILE8 takes the control motion and the "
                          "interpolation settings from it (run AFWRITE with MOTION enabled)")
    if int(md["type"]) != 0:
        raise ModuleError("LOADGEN needs a seismic analysis (MOTION <type> 0)")
    if int(md["srss"]):
        raise ModuleError("MOTION uses SRSS transfer functions (MOTIONX <srss> 1): use LOADGEN source RESULTS")
    nfft, dt = int(md["nft"]), float(md["delt"])
    df = 1.0 / (nfft * dt)
    f8name = str(md["file8"] or "FILE8")
    p8 = workdir / f8name
    if not p8.exists():
        raise ModuleError(f"{f8name} missing -- run ANALYS (source FILE8)")
    f8 = read_file8(p8, f8name)
    if abs(float(f8.meta["df"]) - df) > 1e-6 * df:
        raise ModuleError(f"{f8name} frequency step {float(f8.meta['df']):.9g} Hz differs from the MOTION deck "
                          f"1/(NFFT dt) = {df:.9g} Hz")
    if int(f8.meta["type"]) != 0:
        raise ModuleError(f"{f8name} is not a seismic solution (type {f8.meta['type']})")
    if str(f8.meta.get("model_hash", "")) and str(info.meta.get("model_hash", "")) and \
            str(f8.meta.get("model_hash")) != str(info.meta.get("model_hash")):
        lst.warning(f"{f8name} model_hash differs from FILE4: the solution may come from another HOUSE run")
    cm, ang = int(f8.meta.get("cm", md["cm"])), float(f8.meta.get("ang", md["ang"]))
    padded, A, thname = control_motion_from_deck(md, workdir)
    gravity = info.gravity
    option = int(md["interp"])
    setup = MotionSetup(nfft=nfft, dt=dt, nout=S.output_length(nfft, dt, float(md["dur"])), gravity=gravity,
                        option=option, smooth=0.0 if option == 6 else float(md["smo"]), pzadj=int(md["pzadj"]),
                        bl=bool(int(md["bl"])), cm=cm, ang=ang, A=A, padded=padded, thfile=thname)
    fnum = np.asarray(f8["fnum"], dtype=np.int64)
    f_ssi = fnum * df
    colmap = {(int(n), int(k)): i for i, (n, k) in enumerate(zip(f8["eq_node"], f8["eq_dof"]))}
    acc_req = sorted(set(acc_req))
    disp_req = sorted(set(disp_req))
    ref_cols: List[Tuple[int, int]] = []
    if refnode_f4:
        ref_cols = [(refnode_f4, k) for k in (1, 2, 3) if (refnode_f4, k) in colmap]
    need = sorted(set(acc_req) | set(disp_req) | set(ref_cols))
    absent = [q for q in need if q not in colmap]
    if absent:
        raise ModuleError("node DOFs without an equation in FILE8 (fixed DOFs or wrong node numbers): "
                          + ", ".join(f"{n} {C.DOF_TAGS[k]}" for n, k in absent[:20]))
    f_grid = S.fourier_grid(nfft, dt)
    fN = min(f_ssi[-1], f_grid[-1])
    kmax = int(np.searchsorted(f_grid, fN * (1 + 1e-12), side="right"))
    cols = np.array([colmap[q] for q in need], dtype=np.int64)
    H = np.asarray(f8["H"])[:, cols]
    h0 = np.array([I.rigid_body_anchor(k, cm, ang) for _, k in need], dtype=complex)
    Hg = np.zeros((len(f_grid), len(need)), dtype=complex)
    for b0 in range(0, len(need), 64):                     # MOTION processes 64 columns at a time
        sl = slice(b0, b0 + 64)
        Hg[:, sl] = I.interpolate_tf(f_ssi, H[:, sl], f_grid, option, smooth=setup.smooth, pzadj=setup.pzadj,
                                     h0=h0[sl], mode="seismic")
    Hg = S.hermitian_bins(Hg, nfft)
    Hg[kmax:] = 0.0                                        # .TFI rows stop at f_N (RELDISP reads them)
    pos = {q: j for j, q in enumerate(need)}
    hs = HistorySet(dt=dt, n=setup.nout)
    if acc_req:
        jj = [pos[q] for q in acc_req]
        a = np.fft.irfft(S.hermitian_bins(Hg[:, jj] * A[:, None], nfft), n=nfft, axis=0)
        for c, q in enumerate(acc_req):
            v = a[:, c]
            if setup.bl:
                v = S.baseline_correction(v, dt, scale=gravity).acc
            hs.acc[q] = v[:setup.nout] * gravity
    if disp_req:
        drive = S.displacement_spectrum(A, f_grid, gravity)
        Href = np.zeros((len(f_grid), len(disp_req)), dtype=complex)
        for c, (n_, k) in enumerate(disp_req):
            if k <= 3:
                Href[:, c] = Hg[:, pos[(refnode_f4, k)]] if refnode_f4 and (refnode_f4, k) in pos \
                    else I.rigid_body_anchor(k, cm, ang)
        jj = [pos[q] for q in disp_req]
        u = np.fft.irfft(S.hermitian_bins((Hg[:, jj] - Href) * drive[:, None], nfft), n=nfft, axis=0)
        for c, q in enumerate(disp_req):
            hs.disp[q] = u[:setup.nout, c].copy()
    if refnode_f4:
        jj = [pos[q] for q in ref_cols]
        a = np.fft.irfft(S.hermitian_bins(Hg[:, jj] * A[:, None], nfft), n=nfft, axis=0)
        for c, q in enumerate(ref_cols):
            v = a[:, c]
            if setup.bl and ref_baseline:
                v = S.baseline_correction(v, dt, scale=gravity).acc
            hs.ref[q[1]] = v[:setup.nout] * gravity
    hs.origin["histories"] = (f"{f8name} + control motion {thname} (MOTION deck: NFFT {nfft}, dt {dt:g} s, "
                              f"interpolation {option}, baseline {'on' if setup.bl else 'off'})")
    return hs, setup


# =======================================================================================
# Reference motion of the dynamic analysis (ACEL) and consistency of the D reference
# =======================================================================================
def reldisp_reference(workdir: Path, model: str) -> Tuple[Optional[str], Optional[int]]:
    """Reference of the RELDISP ``.THD`` files from the RELDISP deck: ('FREEFIELD', None), ('NODE', node in
    FILE8 numbering) or (None, None) when the deck is missing."""
    rd = _deck_or_none(workdir / f"{model}.rdi", "RELDISP")
    if rd is None:
        return None, None
    rel = str(rd["relfile"] or "").strip()
    if rel == "" or rel.upper() == "FREEFIELD":
        return "FREEFIELD", None
    m = re.match(r"^(\d{5,6})(TR_X|TR_Y|TR_Z|R_XX|R_YY|R_ZZ)\.TF[IU]$", Path(rel).name, re.IGNORECASE)
    return "NODE", (int(m.group(1)) if m else -1)


def ground_motion(d: LoadgenDeck, workdir: Path, model: str, dt: float, n: int,
                  setup: Optional[MotionSetup]) -> Tuple[np.ndarray, str]:
    """Reference (free-field) acceleration in g on the history grid: the Ground Acceleration File, or the
    control motion of the MOTION (source FILE8) / RELDISP / MOTION / STRESS deck."""
    gname = str(d["groundfile"] or "").strip()
    if gname:
        p = _resolve(gname, workdir)
        if not p.exists():
            raise ModuleError(f"ground acceleration file {p} not found")
        try:
            a, dtf = thfile.read_history(p, fopt=int(d["groundfopt"]))
        except ValueError as exc:
            raise ModuleError(f"{p.name}: {exc}") from None
        if dtf is None or abs(dtf - dt) > DT_RTOL * dt:
            raise ModuleError(f"{p.name}: time step {dtf} s differs from the SSI histories ({dt:g} s)")
        a = np.asarray(a, dtype=float) * float(d["groundmult"])
        out = np.zeros(n)
        out[:min(n, len(a))] = a[:n]
        return out, f"{p.name} x {float(d['groundmult']):g}"
    if setup is not None:
        return setup.padded[:n].copy(), f"control motion {setup.thfile} (MOTION deck scaling)"
    for mod, ext in (("RELDISP", ".rdi"), ("MOTION", ".mot"), ("STRESS", ".str")):
        dk = _deck_or_none(workdir / f"{model}{ext}", mod)
        if dk is None or not str(dk["thfile"] or "").strip():
            continue
        if abs(float(dk["delt"]) - dt) > DT_RTOL * dt:
            raise ModuleError(f"{mod} deck time step {float(dk['delt']):g} s differs from the SSI histories ({dt:g} s)")
        padded, _, name = control_motion_from_deck(dk, workdir)
        return padded[:n].copy() if len(padded) >= n else np.pad(padded, (0, n - len(padded))), \
            f"control motion {name} ({mod} deck scaling)"
    raise ModuleError("no ground acceleration: give the Ground Acceleration File (LGFILE,GROUND,<file>) or run "
                      "AFWRITE with MOTION/RELDISP enabled (their control motion is used)")


# =======================================================================================
# Critical times (D-LGN-07)
# =======================================================================================
def pick_peaks(y: np.ndarray, n: int, sep: int = 0) -> List[int]:
    """Indices of the ``n`` largest local maxima of ``|y|``, at least ``sep`` samples apart, by decreasing
    value.  A local maximum is a sample not smaller than its neighbours (the end samples count)."""
    a = np.abs(np.asarray(y, dtype=float))
    if a.size == 0:
        return []
    left = np.concatenate([[-np.inf], a[:-1]])
    right = np.concatenate([a[1:], [-np.inf]])
    cand = np.flatnonzero((a >= left) & (a >= right) & (a > 0))
    # plateaus: keep the first sample of a run of equal maxima
    cand = cand[np.concatenate([[True], np.diff(cand) > 1])] if cand.size else cand
    order = cand[np.argsort(-a[cand], kind="stable")]
    out: List[int] = []
    for k in order:
        if all(abs(int(k) - j) >= max(sep, 1) for j in out):
            out.append(int(k))
        if len(out) >= n:
            break
    return out


# =======================================================================================
# APDL helpers (D-LGN-08)
# =======================================================================================
def fmt_value(v: float, digits: int = 12) -> str:
    """APDL number with ``digits`` significant digits (``0`` exactly for zero)."""
    v = float(v)
    if not math.isfinite(v):
        raise ModuleError(f"non-finite value {v} in the APDL output (check the SSI results)")
    if v == 0.0:
        return "0"
    return f"{v:.{int(digits) - 1}E}"


def _exact(v: float) -> str:
    """Shortest text that reads back as the same double (time step, end time)."""
    return repr(float(v)).upper()


def table_lines(name: str, values: np.ndarray, digits: int, comment: str = "") -> List[str]:
    """``*DIM`` of a 1-column TABLE with primary variable TIME, its time column (``*VFILL ... RAMP,0,LG_DT``)
    and the values, at most :data:`VALUES_PER_LINE` per assignment ``NAME(i,1)=v1,...``."""
    v = np.asarray(values, dtype=float)
    out = [f"*DIM,{name},TABLE,{len(v)},1,1,TIME" + (f"   ! {comment}" if comment else ""),
           f"*VFILL,{name}(1,0),RAMP,0,LG_DT"]
    for i in range(0, len(v), VALUES_PER_LINE):
        out.append(f"{name}({i + 1},1)=" + ",".join(fmt_value(x, digits) for x in v[i:i + VALUES_PER_LINE]))
    return out


def apdl_name(prefix: str, node: int, label: str) -> str:
    name = f"{prefix}_{int(node)}_{label}"
    if len(name) > APDL_NAME_MAX:
        raise ModuleError(f"APDL parameter name {name} is longer than {APDL_NAME_MAX} characters")
    return name


def _comment_block(lines: Sequence[str], width: int = 96) -> List[str]:
    return ["! " + "=" * width] + [("! " + ln) if ln else "!" for ln in lines] + ["! " + "=" * width]


def _head(items: Sequence[Tuple[str, str]], width: int = 96) -> List[str]:
    """Header lines ``Label      : text`` wrapped at ``width``; the label ``>`` keeps the text verbatim and
    indented (APDL command lines of the instructions), ``^`` verbatim at the margin (title)."""
    import textwrap
    out: List[str] = []
    for label, text in items:
        if label.startswith(">"):
            out.append(" " * 13 + text)
            continue
        if label.startswith("^"):
            out.append(text)
            continue
        prefix = f"{label:<11s}: " if label else " " * 13
        parts = textwrap.wrap(text, width=width - 13, break_long_words=False, break_on_hyphens=False) or [""]
        out.append(prefix + parts[0])
        out += [" " * 13 + w for w in parts[1:]]
    return out


# =======================================================================================
# The run
# =======================================================================================
@dataclass
class Plan:
    """Resolved run data shared by the static and dynamic writers."""
    d: LoadgenDeck
    info: ModelInfo
    nmap: NodeMap
    ansys_dir: Path
    ssidir: Path
    dnodes: List[int]                      # FILE4 numbering
    dnode_src: str
    cm: int = 0
    ang: float = 0.0
    digits: int = 12
    files: List[str] = field(default_factory=list)


def _deck_flag(d: LoadgenDeck, name: str, allowed=(0, 1)) -> int:
    v = int(d[name])
    if v not in allowed:
        raise ModuleError(f"LOADGEN <{name}> = {v} must be one of {', '.join(map(str, allowed))}")
    return v


def validate_deck(d: LoadgenDeck) -> None:
    """Value checks of the LOADGEN deck (the commands repeat them as warnings)."""
    if str(d["analysis"]).upper() not in ("STATIC", "DYNAMIC"):
        raise ModuleError(f"<analysis> {d['analysis']!r} must be STATIC or DYNAMIC")
    _deck_flag(d, "data", (1, 2, 3, 4))
    for k in ("multi", "rotdisp", "rotacc", "genmass", "opmode"):
        _deck_flag(d, k)
    _deck_flag(d, "masstype", (1, 2))
    if str(d["source"]).upper() not in SOURCES:
        raise ModuleError(f"<source> {d['source']!r} must be RESULTS or FILE8")
    if str(d["method"]).upper() not in METHODS:
        raise ModuleError(f"<method> {d['method']!r} must be REL or ACC")
    if float(d["alpha"]) < 0 or float(d["beta"]) < 0:
        raise ModuleError("Rayleigh damping coefficients alpha and beta must be >= 0")
    if str(d["crit"]).upper() not in CRITERIA:
        raise ModuleError(f"critical-time criterion {d['crit']!r} must be one of {' '.join(CRITERIA)}")
    if int(d["ncrit"]) < 1:
        raise ModuleError("<ncrit> must be >= 1")
    if float(d["tsep"]) < 0:
        raise ModuleError("<tsep> must be >= 0")
    if not 1 <= int(d["critdof"]) <= 6:
        raise ModuleError("<critdof> must be 1..6")
    if not 6 <= int(d["digits"]) <= 17:
        raise ModuleError("<digits> must be 6..17")
    _deck_flag(d, "rest")
    if len(d["refpoint"]) not in (0, 3):
        raise ModuleError("<refpoint> must be empty or x y z")


def _f4(info: ModelInfo, model_node: int, what: str) -> int:
    n = info.file4_of.get(int(model_node))
    if n is None:
        raise ModuleError(f"{what} node {model_node} is not a node of the SSI model (FILE4)")
    return n


def _dnode_set(d: LoadgenDeck, info: ModelInfo, data: int, dynamic: bool) -> Tuple[List[int], str]:
    """Nodes that receive D (FILE4 numbering): the LGNODE,D list or the default set (interaction nodes;
    "Disp. for Soil Module": interaction nodes and the nodes of the excavated soil elements)."""
    if d.tables["dnodes"]:
        return [_f4(info, n, "interface") for n in d.tables["dnodes"]], "LGNODE,D list"
    if not dynamic and data == 4:
        nodes = sorted(set(info.int_nodes) | set(info.excav_nodes))
        return nodes, "interaction nodes and nodes of the excavated soil elements (Disp. for Soil Module)"
    return list(info.int_nodes), "interaction nodes (default)"


def _units_line(info: ModelInfo) -> str:
    sysname = "British-like (ft, s)" if info.gravity > 20.0 else "SI-like (m, s)"
    return (f"the SSI model's consistent units, {sysname}: g = {info.gravity:g} length/s^2; masses = weight / g; "
            "SSI accelerations in g are multiplied by g")


def run(ctx: ModuleContext) -> int:
    """LOADGEN module entry point (manual 6.4.15, spec 04 section 15.6)."""
    lst = ctx.listing
    if ctx.deck_path is None or not ctx.deck_path.exists():
        raise ModuleError(f"input deck {ctx.model}{DECK_EXT} not found -- run RUNLOADGEN (it writes the deck)")
    try:
        d, notes = read_deck(ctx.deck_path)
    except ValueError as exc:
        raise ModuleError(str(exc)) from None
    for t in notes:
        lst.warning(t)
    validate_deck(d)
    if d["title"]:
        lst.write(f" {d['title']}")
    dynamic = str(d["analysis"]).upper() == "DYNAMIC"
    ssidir = _resolve(str(d["ssipath"]), ctx.workdir) if str(d["ssipath"]).strip() else ctx.workdir
    ansys_dir = _resolve(str(d["ansyspath"]), ctx.workdir) if str(d["ansyspath"]).strip() else ctx.workdir
    hname = str(d["housefile"] or "").strip() or f"{ctx.model}.hou"
    hpath = _resolve(hname, ssidir)
    info = read_model_info(ssidir, ctx.model, hpath)
    if not hpath.exists():
        lst.warning(f"HOUSE Module Input {hpath.name} not found: PLANE groups and MT/MR masses are taken from FILE4 "
                    "and COOSM only")
    if info.gravity <= 0:
        raise ModuleError("FILE4 has no gravity (HOUSE <gravity>): cannot convert accelerations in g")
    ansys_dir.mkdir(parents=True, exist_ok=True)
    nmap = build_node_map(d, info, ansys_dir, ctx.workdir, lst)
    data = int(d["data"])
    dnodes, dsrc = _dnode_set(d, info, data, dynamic)
    plan = Plan(d=d, info=info, nmap=nmap, ansys_dir=ansys_dir, ssidir=ssidir, dnodes=dnodes, dnode_src=dsrc,
                digits=int(d["digits"]))

    lst.section("LOADGEN options (Option A, ACS SASSI-ANSYS two-step approach)")
    lst.write(f"   analysis                 : {'ANSYS Dynamic Load (transient, direct integration)' if dynamic else 'ANSYS Eq. Static Load'}")
    if not dynamic:
        lst.write(f"   data to add              : {data} {DATA_NAMES[data]}; multiple files {int(d['multi'])}")
        lst.write(f"   mass                     : {MASS_NAMES[int(d['masstype'])]}, generate {int(d['genmass'])}")
    else:
        lst.write(f"   method                   : {str(d['method']).upper()}; Rayleigh alpha {float(d['alpha']):g}, "
                  f"beta {float(d['beta']):g}; reference node {int(d['refnode']) or 'none (control motion)'}")
    lst.write(f"   rotations                : displacements {int(d['rotdisp'])}, accelerations {int(d['rotacc'])}")
    lst.write(f"   history source           : {str(d['source']).upper()}")
    lst.write(f"   SSI results path         : {ssidir}")
    lst.write(f"   ANSYS path               : {ansys_dir}")
    lst.write(f"   model                    : {info.model} ({info.title}); {len(info.xyz)} nodes, "
              f"{len(info.eq_dof)} equations, {'2-D (X-Z -> ANSYS X-Y)' if info.two_d else '3-D'}; g = {info.gravity:g}")
    if info.optimized:
        lst.write("   numbering                : HOUSE optimizer: FILE4/FILE8/MOTION/RELDISP use the new numbers; LOADGEN "
                  "translates them to the model numbers of the .hou (FILE4 x_node_old_id)")
    lst.write(f"   node map to ANSYS        : {nmap.mode}: {nmap.source}")
    for t in nmap.notes:
        lst.write(f"      {t}")
    lst.write(f"   interface (D) nodes      : {len(dnodes)}, {dsrc}")
    if not dnodes:
        if (not dynamic and data in (1, 4)) or (dynamic and str(d["method"]).upper() == "REL"):
            raise ModuleError("no interface nodes for the D commands (no interaction nodes in FILE4 and no LGNODE,D "
                              "list)")
        lst.warning("no interface nodes: the ANSYS model must provide its own supports")
    ctx.progress(0.1, "LOADGEN: inputs read")
    if dynamic:
        _run_dynamic(ctx, plan)
    else:
        _run_static(ctx, plan)
    ctx.progress(1.0, "LOADGEN: done")
    lst.section("Files written")
    for f in plan.files:
        lst.write(f"   {f}")
    if not plan.files:
        lst.write("   none (data check)" if int(d["opmode"]) else "   none")
    return 0


# ---------------------------------------------------------------------------------------
# helpers shared by static and dynamic
# ---------------------------------------------------------------------------------------
def _model_node(plan: Plan, f4: int) -> int:
    return plan.info.model_of[int(f4)]


def _ansys_node(plan: Plan, f4: int, missing: List[int]) -> Optional[int]:
    a = plan.nmap.ansys(_model_node(plan, f4))
    if a is None:
        missing.append(_model_node(plan, f4))
    return a


def _check_missing_map(plan: Plan, missing: List[int], what: str, lst: Listing) -> None:
    if missing:
        u = sorted(set(missing))
        raise ModuleError(f"{len(u)} {what} node(s) have no ANSYS node in the node map ({plan.nmap.mode}): "
                          + ", ".join(map(str, u[:20])) + (" ..." if len(u) > 20 else ""))


def _d_dofs(plan: Plan, f4: int, rot: bool) -> List[int]:
    return [k for k in plan.info.dofs.get(int(f4), []) if k <= 3 or rot]


def _apdl_path(plan: Plan, dynamic: bool) -> Path:
    name = str(plan.d["apdlfile"] or "").strip() or f"{plan.info.model}_{'LGD' if dynamic else 'LGS'}.inp"
    p = _resolve(name, plan.ansys_dir)
    if not p.suffix:
        p = p.with_suffix(".inp")
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def _direction(plan: Plan, workdir: Path, setup: Optional[MotionSetup]) -> Tuple[int, float, str]:
    """Control direction (cm, ang): FILE8 (via the MOTION setup), else the MOTION / RELDISP deck."""
    if setup is not None:
        return setup.cm, setup.ang, "FILE8"
    for mod, ext in (("MOTION", ".mot"), ("RELDISP", ".rdi"), ("STRESS", ".str")):
        dk = _deck_or_none(plan.ssidir / f"{plan.info.model}{ext}", mod)
        if dk is not None:
            cm, ang = int(dk["cm"]), float(dk["ang"])
            f8 = plan.ssidir / str(dk["file8"] or "FILE8")
            if f8.exists():
                try:
                    meta = read_container(f8).meta
                    return int(meta.get("cm", cm)), float(meta.get("ang", ang)), f8.name
                except (OSError, ValueError):
                    pass
            return cm, ang, f"{mod} deck"
    return 0, 0.0, "default (x)"


def _delt_hint(plan: Plan) -> Optional[float]:
    for mod, ext in (("MOTION", ".mot"), ("RELDISP", ".rdi")):
        dk = _deck_or_none(plan.ssidir / f"{plan.info.model}{ext}", mod)
        if dk is not None:
            return float(dk["delt"])
    return None


def _histories(ctx: ModuleContext, plan: Plan, acc_req, disp_req, refnode_f4: int = 0, ref_baseline: bool = True
               ) -> Tuple[HistorySet, Optional[MotionSetup]]:
    lst = ctx.listing
    if str(plan.d["source"]).upper() == "FILE8":
        return file8_histories(plan.info, plan.ssidir, plan.info.model, acc_req, disp_req, refnode_f4, lst,
                               ref_baseline=ref_baseline)
    if refnode_f4 and not ref_baseline:
        md = _deck_or_none(plan.ssidir / f"{plan.info.model}.mot", "MOTION")
        if md is not None and int(md["bl"]):
            lst.warning("MOTION applied the Hudson-Housner baseline correction (<bl> 1) to the .ACC of the reference "
                        "node, RELDISP does not correct the relative displacements: ACEL and D refer to slightly "
                        "different motions -- use source FILE8 (exact) or MOTION <bl> 0")
    acc_req = set(acc_req) | ({(refnode_f4, k) for k in (1, 2, 3) if k in plan.info.dofs.get(refnode_f4, [])}
                              if refnode_f4 else set())
    hs = load_results(plan.d, plan.info, plan.ssidir, acc_req, disp_req, _delt_hint(plan), lst)
    return hs, None


def start_at_rest(hs: HistorySet) -> float:
    """Subtract the initial value of every relative-displacement history (D-LGN-06); returns the largest
    removed value relative to the largest peak.

    MOTION/RELDISP work with the FFT, whose zero-frequency term of a displacement is undetermined and set
    to 0 (D-CNV-06): the periodic solution is the physical one minus its mean over the Fourier period.  The
    structure is at rest before the record starts, so the physical history is ``u(t) - u(0)`` as long as the
    quiet zone lets the response decay (requirements 4.0.1).  Without the correction the ANSYS transient
    would start with a jump (zero initial conditions) and the static interface displacements would carry a
    different constant at every node."""
    if not hs.disp:
        return 0.0
    peak = max(float(np.max(np.abs(v))) for v in hs.disp.values()) or 1.0
    off = max(abs(float(v[0])) for v in hs.disp.values())
    for k in list(hs.disp):
        hs.disp[k] = hs.disp[k] - hs.disp[k][0]
    return off / peak


def _rest(plan: Plan, hs: HistorySet, lst: Listing) -> None:
    if not hs.disp:
        return
    if int(plan.d["rest"]):
        r = start_at_rest(hs)
        lst.write(f"   start at rest            : initial values of the relative displacements subtracted (largest "
                  f"{r:.2%} of the peak; the zero-mean offset of the periodic SSI solution)")
    else:
        peak = max(float(np.max(np.abs(v))) for v in hs.disp.values()) or 1.0
        off = max(abs(float(v[0])) for v in hs.disp.values())
        if off > 1e-3 * peak:
            lst.warning(f"the relative displacements do not start at rest (|u(0)| up to {off / peak:.2%} of the peak) "
                        "and <rest> = 0: ANSYS starts from zero initial conditions (LGOPT,,,1 subtracts u(0))")


def _check_reference(plan: Plan, refnode_f4: int, needs_disp: bool, rot: bool, lst: Listing,
                     static: bool = False) -> None:
    """RESULTS source: the RELDISP ``.THD`` files must be relative to a rigid-body reference that agrees
    with LOADGEN's.  Dynamic: the ACEL motion (free field or ``refnode``) must be the RELDISP reference.
    Static: any translation reference is a rigid-body motion (no stress), but rotations relative to a
    node are not."""
    if str(plan.d["source"]).upper() != "RESULTS" or not needs_disp:
        return
    kind, node = reldisp_reference(plan.ssidir, plan.info.model)
    if kind is None:
        lst.warning("RELDISP deck not found: the .THD files are assumed to be relative to the free field")
        return
    if kind == "NODE" and rot:
        raise ModuleError("the RELDISP .THD rotations are relative to a reference node (RELFILE): that is not a "
                          "rigid-body reference -- use source FILE8 (rotations stay absolute) or the free field")
    if static:
        if kind == "NODE":
            lst.write(f"   RELDISP reference node {node}: the .THD translations differ from the total ones by a rigid "
                      "translation (no stress)")
        return
    if kind == "FREEFIELD" and refnode_f4:
        raise ModuleError("the RELDISP .THD files are relative to the free field but LOADGEN uses reference node "
                          f"{_model_node(plan, refnode_f4)}: set the reference node to 0 or use source FILE8")
    if kind == "NODE" and (not refnode_f4 or node != refnode_f4):
        raise ModuleError(f"the RELDISP .THD files are relative to node {node} (RELFILE) but LOADGEN uses "
                          + (f"reference node {_model_node(plan, refnode_f4)}" if refnode_f4 else "the free field")
                          + ": make them agree (LOADGENDYN <refnode>) or use source FILE8")


# ---------------------------------------------------------------------------------------
# static
# ---------------------------------------------------------------------------------------
@dataclass
class StaticStep:
    """One equivalent static load step: critical sample ``k`` (time ``t``), its rank and criterion value, the
    ANSYS ``F`` and ``D`` commands (ANSYS node, label, value) and the resultant of the inertia forces."""
    rank: int
    k: int
    t: float
    value: float
    forces: List[Tuple[int, str, float]] = field(default_factory=list)
    disps: List[Tuple[int, str, float]] = field(default_factory=list)
    resultant: np.ndarray = field(default_factory=lambda: np.zeros(6))


@dataclass
class InertiaLoad:
    """Inertia force history ``f`` (SASSI global axes, model units) of DOF ``dof`` (1..6) at ``node``
    (lumped: FILE4 node, applied at its ANSYS node; master: ANSYS load node ``ansys`` at ``xyz``)."""
    node: int
    dof: int
    f: np.ndarray
    xyz: np.ndarray
    ansys: Optional[int] = None


def _masses(ctx: ModuleContext, plan: Plan) -> Tuple[List[LoadMass], str]:
    """Load masses (D-LGN-03): generated from COOSM (or the HOUSE deck masses) or read from the mass file."""
    d, info, lst = plan.d, plan.info, ctx.listing
    mtype = int(d["masstype"])
    given = str(d["lumpfile"] if mtype == 1 else d["masterfile"]).strip()
    path = _resolve(given or f"{info.model}.{'masl' if mtype == 1 else 'masm'}", plan.ansys_dir)
    if not int(d["genmass"]):
        if not path.exists():
            raise ModuleError(f"mass file {path} not found: check Generate Mass Data (LOADGEN <genmass> = 1) or give "
                              f"the file (LGFILE,{'LUMPED' if mtype == 1 else 'MASTER'},<file>)")
        if mtype == 1:
            return [LoadMass(n, v) for n, v in sorted(read_lumped(path).items())], f"{path.name} (read)"
        return read_master(path), f"{path.name} (read)"
    pm = plan.ssidir / "COOSM"
    if pm.exists():
        try:
            Ms = read_container(pm, kind="COOSM").sparse("Ms")
        except (OSError, ValueError, KeyError) as exc:
            raise ModuleError(f"COOSM: {exc}") from None
        lumped = lumped_from_matrix(info, Ms)
        how = "generated from the HOUSE structure mass matrix COOSM Ms (row sums per direction, diagonal rotary)"
    else:
        if not info.house_masses:
            raise ModuleError("COOSM missing and no MT/MR masses in the HOUSE deck: run HOUSE before generating masses")
        lumped = {n: v.copy() for n, v in info.house_masses.items() if np.any(v)}
        how = "generated from the MT/MR masses of the HOUSE deck (COOSM missing: element masses not included)"
        lst.warning("COOSM missing: only the MT/MR nodal masses are used (element masses not included)")
    check = bool(int(d["opmode"]))
    if mtype == 1:
        if not check:
            write_lumped(path, lumped, info, how)
            plan.files.append(f"{path.name} (lumped masses, {len(lumped)} nodes)")
        return [LoadMass(n, v) for n, v in sorted(lumped.items())], f"{path.name} ({how})"
    masters = list(d.tables["masters"])
    if not masters:
        raise ModuleError("Master Node Mass generation needs the master nodes (LGNODE,M,<list>)")
    assign = master_assignment(info, lumped, masters)
    masses = []
    for n in sorted(lumped):
        a = plan.nmap.ansys(n)
        if a is None:
            raise ModuleError(f"node {n} has no ANSYS node in the node map: it cannot be written as a load node")
        masses.append(LoadMass(a, lumped[n], master=assign[n], xyz=info.xyz[info.file4_of[n]].copy()))
    if not check:
        write_master(path, masses, info, how + "; master = the master node closest in elevation, then in plan")
        plan.files.append(f"{path.name} (master-node masses, {len(masses)} load nodes, {len(masters)} masters)")
    return masses, f"{path.name} ({how})"


def _inertia_requests(plan: Plan, masses: List[LoadMass], skip: set, rotacc: bool
                      ) -> Tuple[List[Tuple[LoadMass, int, List[int]]], np.ndarray, float]:
    """(load mass, FILE4 node of the driving motion, DOFs) of every load node; the translational mass of the
    nodes in ``skip`` (their DOFs are prescribed: the inertia goes to the reactions); the rotary inertia
    left out because Rotational Accel. is off."""
    info = plan.info
    mtype = int(plan.d["masstype"])
    out: List[Tuple[LoadMass, int, List[int]]] = []
    excluded = np.zeros(3)
    skipped_rot = 0.0
    for lm in masses:
        if mtype == 1:
            drv = _f4(info, lm.node, "mass")
            have = info.dofs.get(drv, [])
            if drv in skip:
                excluded += lm.m[:3]
                continue
            dofs = [k for k in have if lm.m[k - 1] != 0.0 and (k <= 3 or rotacc)]
            if not rotacc:
                skipped_rot += float(sum(lm.m[k - 1] for k in have if k >= 4))
        else:
            drv = _f4(info, lm.master, "master")
            have = info.dofs.get(drv, [])
            dofs = [k for k in have if k <= 3 or rotacc]
        out.append((lm, drv, dofs))
    return out, excluded, skipped_rot


def inertia_histories(plan: Plan, loads: List[Tuple[LoadMass, int, List[int]]], hs: HistorySet,
                      rotacc: bool) -> List[InertiaLoad]:
    """Inertia force histories ``F = -m a`` (D-LGN-03).  Lumped Mass: ``a`` of the node itself.  Master Node
    Mass: rigid-body motion of the master, ``a_k = a_M + alpha_M x (x_k - x_M)`` (alpha_M = 0 without the
    rotational accelerations) and the moments ``-I alpha_M``."""
    info = plan.info
    mtype = int(plan.d["masstype"])
    n = hs.n
    out: List[InertiaLoad] = []
    for lm, drv, dofs in loads:
        if mtype == 1:
            for k in dofs:
                out.append(InertiaLoad(drv, k, -lm.m[k - 1] * hs.acc[(drv, k)], info.xyz[drv]))
            continue
        aM = np.stack([hs.acc.get((drv, k), np.zeros(n)) for k in (1, 2, 3)], axis=1)
        alM = np.stack([hs.acc.get((drv, k), np.zeros(n)) if rotacc else np.zeros(n) for k in (4, 5, 6)], axis=1)
        with _quiet():
            a = aM + np.cross(alM, (lm.xyz - info.xyz[drv])[None, :])
        for k in (1, 2, 3):
            if lm.m[k - 1] != 0.0 and not (info.two_d and k == 2):
                out.append(InertiaLoad(drv, k, -lm.m[k - 1] * a[:, k - 1], lm.xyz, ansys=lm.node))
        if rotacc:
            for k in (4, 5, 6):
                if lm.m[k - 1] != 0.0 and k in info.dofs.get(drv, []):
                    out.append(InertiaLoad(drv, k, -lm.m[k - 1] * alM[:, k - 4], lm.xyz, ansys=lm.node))
    return out


def resultants(loads: List[InertiaLoad], n: int, x0: np.ndarray) -> np.ndarray:
    """(n, 6) histories VX VY VZ MX MY MZ of the inertia forces about ``x0`` (SASSI global axes):
    ``V = sum F``, ``M = sum (x - x0) x F + sum M_node``."""
    R = np.zeros((n, 6))
    for ld in loads:
        if ld.dof <= 3:
            e = np.zeros(3)
            e[ld.dof - 1] = 1.0
            R[:, ld.dof - 1] += ld.f
            with _quiet():
                R[:, 3:] += np.outer(ld.f, np.cross(ld.xyz - x0, e))
        else:
            R[:, ld.dof - 1] += ld.f
    return R


def _critical_steps(plan: Plan, crit: str, R: np.ndarray, hs: HistorySet, crit_q, x0: np.ndarray,
                    e_cm: np.ndarray, lst: Listing) -> Tuple[List[int], Optional[np.ndarray], str]:
    """Critical samples (D-LGN-07): the largest local maxima of ``|criterion|`` (at least ``tsep`` apart) or
    the given times / steps."""
    d = plan.d
    dt, n = hs.dt, hs.n
    if crit in ("TIME", "STEP"):
        times = list(d["times"])
        if not times:
            raise ModuleError(f"criterion {crit} needs the times (LGTIME,{crit},<t1>,<t2>,...)")
        ks = []
        for t in times:
            k = int(round(t / dt)) if crit == "TIME" else int(round(t)) - 1
            if crit == "TIME" and abs(t - k * dt) > 1e-6 * dt:
                lst.warning(f"critical time {t:g} s is not on the time grid: the nearest step, t = {k * dt:g} s, is used")
            if not 0 <= k < n:
                raise ModuleError(f"critical {'time' if crit == 'TIME' else 'step'} {t:g} is outside the histories "
                                  f"(0 .. {(n - 1) * dt:g} s, steps 1 .. {n})")
            if k not in ks:
                ks.append(k)
        return ks, None, f"given {'times' if crit == 'TIME' else 'time steps'}"
    if crit == "V":
        with _quiet():
            y = R[:, :3] @ e_cm
        yname = "base shear along the control direction (sum of the inertia forces)"
    elif crit in RESULTANT_NAMES:
        y = R[:, RESULTANT_NAMES.index(crit)]
        yname = (f"{crit} of the inertia forces" + (f" about ({x0[0]:g}, {x0[1]:g}, {x0[2]:g})" if crit[0] == "M"
                                                    else ""))
    elif crit == "ACC":
        y = hs.acc[crit_q] / plan.info.gravity
        yname = f"absolute acceleration (g) of node {d['critnode']} {C.DOF_TAGS[crit_q[1]]}"
    else:
        y = hs.disp[crit_q]
        yname = f"relative displacement of node {d['critnode']} {C.DOF_TAGS[crit_q[1]]}"
    ks = pick_peaks(y, int(d["ncrit"]), int(round(float(d["tsep"]) / dt)))
    if not ks:
        raise ModuleError(f"the criterion ({yname}) is zero at every time step: no critical time")
    return ks, y, yname


def _run_static(ctx: ModuleContext, plan: Plan) -> None:
    """ANSYS Eq. Static Load (D-LGN-02): critical times and one load step each."""
    lst, d, info = ctx.listing, plan.d, plan.info
    data = int(d["data"])
    use_acc = data in (2, 3)
    use_disp = data in (1, 3, 4)
    rotacc, rotdisp = bool(int(d["rotacc"])), bool(int(d["rotdisp"]))
    crit = str(d["crit"]).upper()
    needs_mass = use_acc or crit in ("V",) + RESULTANT_NAMES
    masses: List[LoadMass] = []
    mass_src = ""
    if needs_mass:
        masses, mass_src = _masses(ctx, plan)
    # nodes with prescribed DOFs (relative displacements, or the fixed base of "Acceleration")
    skip = set(plan.dnodes) if data in (2, 3) else set()
    loads, excluded, skipped_rot = _inertia_requests(plan, masses, skip, rotacc)
    acc_req = [(drv, k) for _, drv, dofs in loads for k in dofs]
    if int(d["masstype"]) == 2:
        acc_req = [(drv, k) for _, drv, dofs in loads for k in plan.info.dofs.get(drv, []) if k <= 3 or rotacc]
    disp_req: List[Tuple[int, int]] = []
    if use_disp:
        disp_req = [(nd, k) for nd in plan.dnodes for k in _d_dofs(plan, nd, rotdisp)]
    crit_q = None
    if crit in ("ACC", "DISP"):
        cn = _f4(info, int(d["critnode"]), "criterion")
        if int(d["critdof"]) not in info.dofs.get(cn, []):
            raise ModuleError(f"criterion node {d['critnode']} has no equation for DOF {d['critdof']} (fixed or absent)")
        crit_q = (cn, int(d["critdof"]))
        (acc_req if crit == "ACC" else disp_req).append(crit_q)
    _check_reference(plan, 0, bool(disp_req), rotdisp, lst, static=True)
    hs, setup = _histories(ctx, plan, acc_req, disp_req)
    dt, n = hs.dt, hs.n
    if n == 0:
        raise ModuleError("empty SSI histories")
    cm, ang, csrc = _direction(plan, ctx.workdir, setup)
    plan.cm, plan.ang = cm, ang
    e_cm = control_direction(cm, ang)
    lst.section("Inputs")
    lst.write(f"   histories                : {hs.origin.get('histories', '')}")
    lst.write(f"   time grid                : dt {dt:g} s, {n} samples (0 .. {(n - 1) * dt:g} s)")
    lst.write(f"   control direction        : cm {cm}, ang {ang:g} deg -> ({e_cm[0]:.4g}, {e_cm[1]:.4g}, {e_cm[2]:.4g}) "
              f"from {csrc}")
    _rest(plan, hs, lst)
    if masses:
        tot = np.sum([lm.m for lm in masses], axis=0)
        lst.write(f"   masses                   : {mass_src}")
        lst.write(f"                              {len(masses)} load nodes; total mx {tot[0]:.6g}, my {tot[1]:.6g}, "
                  f"mz {tot[2]:.6g}; rotary {tot[3]:.4g}, {tot[4]:.4g}, {tot[5]:.4g}")
        if np.any(excluded):
            lst.write(f"                              mass of the nodes with prescribed DOFs (inertia to the reactions): "
                      f"mx {excluded[0]:.6g}, my {excluded[1]:.6g}, mz {excluded[2]:.6g}")
        if skipped_rot > 0:
            lst.warning(f"rotary inertia {skipped_rot:.4g} not loaded: Rotational Accel. is off (LOADGEN <rotacc> 0)")
        if int(d["masstype"]) == 2 and not rotacc:
            lst.warning("Master Node Mass without Rotational Accel.: the master rotations are taken as zero "
                        "(a = a_M only)")
    x0 = (np.array(d["refpoint"], dtype=float) if len(d["refpoint"]) == 3 else
          np.mean([info.xyz[k] for k in plan.dnodes], axis=0) if plan.dnodes else np.zeros(3))
    F = inertia_histories(plan, loads, hs, rotacc)
    R = resultants(F, n, x0)
    if crit in ("V",) + RESULTANT_NAMES and not F:
        raise ModuleError("no inertia forces (no masses at the load nodes): choose LGTIME,ACC / DISP / TIME / STEP")
    ks, y, yname = _critical_steps(plan, crit, R, hs, crit_q, x0, e_cm, lst)
    multi = bool(int(d["multi"]))
    if not multi and len(ks) > 1:
        lst.warning(f"{len(ks)} critical times requested but Use Multiple File List Inputs is off: only the first "
                    f"(t = {ks[0] * dt:g} s) is written (LOADGEN <multi> 1 writes one file per critical time)")
        ks = ks[:1]
    steps = [StaticStep(rank=i + 1, k=k, t=k * dt, value=float(y[k]) if y is not None else 0.0, resultant=R[k].copy())
             for i, k in enumerate(ks)]
    # ---- ANSYS commands of every step (ANSYS numbering and axes)
    missing: List[int] = []
    f_nodes = [(ld, ld.ansys if ld.ansys is not None else _ansys_node(plan, ld.node, missing)) for ld in F]
    d_nodes = [(nd, _ansys_node(plan, nd, missing)) for nd in plan.dnodes]
    _check_missing_map(plan, missing, "load or interface", lst)
    for st in steps:
        for ld, a in f_nodes:
            ia, sg = ansys_dof(ld.dof, info.two_d)
            st.forces.append((int(a), FORCE_LABELS[ia], sg * float(ld.f[st.k])))
        if use_disp:
            for nd, a in d_nodes:
                for k in _d_dofs(plan, nd, rotdisp):
                    ia, sg = ansys_dof(k, info.two_d)
                    st.disps.append((int(a), DISP_LABELS[ia], sg * float(hs.disp[(nd, k)][st.k])))
        elif data == 2:
            for nd, a in d_nodes:                       # fixed-base supports of the equivalent static model
                for k in _d_dofs(plan, nd, True):
                    ia, _ = ansys_dof(k, info.two_d)
                    st.disps.append((int(a), DISP_LABELS[ia], 0.0))
    if use_disp and not rotdisp and any(k >= 4 for nd in plan.dnodes for k in info.dofs.get(nd, [])):
        lst.warning("interface nodes have rotational DOFs whose rotations are not prescribed (Rotational Disp. off, "
                    "LOADGEN <rotdisp> 0): the ANSYS model leaves them free")
    lst.section("Critical times")
    lst.write(f"   criterion: {yname}" + ("" if y is None else f"; maximum |value| {np.max(np.abs(y)):.6g} at "
                                          f"t = {int(np.argmax(np.abs(y))) * dt:g} s"))
    lst.write(f"{'rank':>6s} {'step':>7s} {'time (s)':>10s} {'criterion':>14s} {'VX':>14s} {'VY':>14s} {'VZ':>14s} "
              f"{'MX':>14s} {'MY':>14s} {'MZ':>14s}")
    for st in steps:
        lst.write(f"{st.rank:6d} {st.k + 1:7d} {st.t:10.4f} {st.value:14.6e} " + " ".join(f"{v:14.6e}" for v in st.resultant))
    lst.write(f"   V and M: resultants of the inertia forces F = -m a (SASSI-EDU global axes), moments about "
              f"({x0[0]:g}, {x0[1]:g}, {x0[2]:g}); the SSI base shear of the structure above the interface is -V")
    if int(d["opmode"]):
        lst.section("Data check")
        lst.write("   input checked; data-check mode writes no APDL file")
        return
    base = _apdl_path(plan, False)
    if F:
        res = base.with_name(base.stem + "_res.txt")
        hdr = (f"SASSI-EDU LOADGEN resultants of the inertia forces F = -m a (model {info.model}, SASSI-EDU global axes); "
               f"moments about ({x0[0]:.10g}, {x0[1]:.10g}, {x0[2]:.10g}); dt {dt:.10g} s\nt VX VY VZ MX MY MZ")
        np.savetxt(res, np.column_stack([np.arange(n) * dt, R]), fmt="%.12e", header=hdr, comments="# ")
        plan.files.append(f"{res.name} (resultant histories of the inertia forces)")
    paths = [base] if len(steps) == 1 and not multi else \
        [base.with_name(f"{base.stem}_{i + 1:02d}{base.suffix}") for i in range(len(steps))]
    for i, (st, path) in enumerate(zip(steps, paths)):
        text = apdl_static(plan, st, i + 1, len(steps), mass_src, yname, data, dt, n, path)
        path.write_text(text, encoding="utf-8")
        plan.files.append(f"{path.name} (load step {i + 1}: t = {st.t:g} s, {len(st.forces)} F, {len(st.disps)} D)")
    lst.section("ANSYS loads")
    for st, path in zip(steps, paths):
        fsum = np.zeros(6)
        for (_, lab, v) in st.forces:
            fsum[FORCE_LABELS.index(lab)] += v
        lst.write(f"   {path.name}: t = {st.t:g} s, {len(st.forces)} F (ANSYS axes: sum FX {fsum[0]:.6e}, FY {fsum[1]:.6e}, "
                  f"FZ {fsum[2]:.6e}), {len(st.disps)} D")


def apdl_static(plan: Plan, st: StaticStep, ls: int, nls: int, mass_src: str, yname: str, data: int, dt: float,
                n: int, path: Path) -> str:
    """APDL of one equivalent static load step (D-LGN-02, D-LGN-08): header, ``LG_SOLVE`` switch,
    ``ANTYPE,STATIC``, ``TIME,<ls>``, the ``F`` and ``D`` commands and the guarded ``SOLVE``."""
    info = plan.info
    digits = plan.digits
    fx = np.zeros(6)
    for (_, lab, v) in st.forces:
        fx[FORCE_LABELS.index(lab)] += v
    items: List[Tuple[str, str]] = [
        ("^", f"{PRODUCT} {__version__}  LOADGEN (Option A): ANSYS equivalent static seismic load, file {ls} of {nls}"),
        ("Model", info.model + (f"  ({info.title})" if info.title else "")),
        ("Written", f"{time.strftime('%Y-%m-%d %H:%M:%S')}  (LOADGEN listing {info.model}_LOADGEN.out)"),
        ("Load step", f"SSI time t = {st.t:.6g} s (sample {st.k + 1} of {n}, dt {dt:g} s); critical time {st.rank} by "
                      f"{yname}"),
        ("Data", DATA_NAMES[data]),
    ]
    if st.forces:
        items.append(("", f"inertia forces F = -m a(t) of the ABSOLUTE accelerations: {len(st.forces)} components at "
                          f"{len({f[0] for f in st.forces})} nodes; masses: {mass_src}"))
    if data in (1, 3, 4) and st.disps:
        items.append(("", f"SSI displacements RELATIVE to the free field at {len({x[0] for x in st.disps})} nodes "
                          f"({plan.dnode_src}); they differ from the total ones by a rigid translation, which causes "
                          "no stress"))
    elif data == 2 and st.disps:
        items.append(("", f"fixed base: D = 0 at the {len({x[0] for x in st.disps})} interface nodes ({plan.dnode_src})"))
    items += [
        ("Resultant", f"sum FX {fx[0]:.6e}, FY {fx[1]:.6e}, FZ {fx[2]:.6e} (ANSYS axes); the reactions balance it"),
        ("Units", _units_line(info)),
        ("Axes", "global Cartesian axes of the SASSI-EDU model (ANSYS CSYS,0), Z up"
         + ("; 2-D model: SASSI X-Z plane -> ANSYS X-Y plane (Y = SASSI Z), as the ANSYS command" if info.two_d else "")),
        ("Nodes", f"{plan.nmap.mode}: {plan.nmap.source}"),
        ("Excluded", "dead load (combine with a gravity load case) and the damping forces of the SSI solution"),
        ("How to run", "read your structural model, then this file, then review the results:"),
        (">", "/INPUT,<model>,inp       ! the ANSYS model (e.g. the output of the SASSI-EDU ANSYS command)"),
        (">", f"/INPUT,{path.stem},{path.suffix.lstrip('.') or 'inp'}   ! this file: loads, then SOLVE (LG_SOLVE = 1)"),
        (">", "/POST1, then SET,LAST     ! stresses and reactions"),
        ("", "The files of one LOADGEN run can be read one after the other: each one re-applies every F and D of its "
             "time and solves a new load step."),
    ]
    head = _head(items)
    L = _comment_block(head)
    L += ["LG_SOLVE = 1                ! 1 = solve at the end of this file; 0 = apply the loads only (then add your "
          "options and SOLVE)",
          f"LG_TSSI = {_exact(st.t)}             ! SSI time of this load step (s)",
          "/SOLU",
          "ANTYPE,STATIC",
          f"TIME,{ls}                      ! load step {ls} <- SSI time {st.t:g} s"]
    if st.forces:
        L += ["!", f"! ---- inertia forces F = -m a(t) ({len(st.forces)} components)"]
        L += [f"F,{a},{lab},{fmt_value(v, digits)}" for a, lab, v in st.forces]
    if st.disps:
        what = ("SSI displacements relative to the free field" if data != 2 else
                "fixed-base supports of the equivalent static model (D = 0)")
        L += ["!", f"! ---- {what} ({len(st.disps)} components)"]
        L += [f"D,{a},{lab},{fmt_value(v, digits)}" for a, lab, v in st.disps]
    L += ["!", "*IF,LG_SOLVE,EQ,1,THEN", "SOLVE", "*ENDIF", "FINISH", ""]
    return "\n".join(L)


# ---------------------------------------------------------------------------------------
# dynamic
# ---------------------------------------------------------------------------------------
def _compare_ground(plan: Plan, ag: np.ndarray, dt: float, n: int, setup: Optional[MotionSetup], lst: Listing) -> None:
    """The D tables are relative to the control motion (free-field reference): a Ground Acceleration File
    that differs from it breaks the exact decomposition u = u_r + iota u_g (warning with the difference)."""
    saved = plan.d["groundfile"]
    try:
        plan.d.params["groundfile"] = ""
        ref, name = ground_motion(plan.d, plan.ssidir, plan.info.model, dt, n, setup)
    except ModuleError:
        return
    finally:
        plan.d.params["groundfile"] = saved
    peak = float(np.max(np.abs(ref))) or 1.0
    diff = float(np.max(np.abs(ag - ref))) / peak
    if diff > 1e-6:
        lst.warning(f"the Ground Acceleration File differs from the {name} by {diff:.2%} of its peak: the D tables are "
                    "relative to the control motion, so ACEL should be the same motion (for a kinematic SSI motion use "
                    "a reference node, LOADGENDYN <refnode>, with source FILE8)")
    else:
        lst.write(f"   ground acceleration file = {name} (max difference {diff:.1e} of the peak)")


def _run_dynamic(ctx: ModuleContext, plan: Plan) -> None:
    lst, d, info = ctx.listing, plan.d, plan.info
    method = str(d["method"]).upper()
    rotdisp = bool(int(d["rotdisp"]))
    rotacc = bool(int(d["rotacc"]))
    refnode = int(d["refnode"])
    ref_f4 = _f4(info, refnode, "reference") if refnode else 0
    if method == "ACC" and not ref_f4:
        raise ModuleError("method ACC needs the reference node whose absolute acceleration drives the fixed base "
                          "(LOADGENDYN,...,ACC,<refnode>)")
    disp_req: List[Tuple[int, int]] = []
    if method == "REL":
        for nd in plan.dnodes:
            disp_req += [(nd, k) for k in _d_dofs(plan, nd, rotdisp)]
    check = [_f4(info, nn, "check") for nn in d.tables["checknodes"]]
    acc_req: List[Tuple[int, int]] = []
    for nd in check:
        acc_req += [(nd, k) for k in info.dofs.get(nd, []) if k <= 3 or rotacc]
    if method == "ACC":
        # interface accelerations, only to report how rigid the foundation is (when they are available)
        file8 = str(d["source"]).upper() == "FILE8"
        for nd in plan.dnodes:
            for k in info.dofs.get(nd, []):
                if k <= 3 and (file8 or find_node_file(plan.ssidir, nd, k, "ACC") is not None):
                    acc_req.append((nd, k))
    _check_reference(plan, ref_f4, bool(disp_req), rotdisp, lst)
    hs, setup = _histories(ctx, plan, acc_req, disp_req, ref_f4, ref_baseline=(method == "ACC"))
    dt, n = hs.dt, hs.n
    cm, ang, csrc = _direction(plan, ctx.workdir, setup)
    plan.cm, plan.ang = cm, ang
    e_cm = control_direction(cm, ang)
    # ---- reference (frame) acceleration in SASSI global axes, model units
    if ref_f4:
        aref = np.stack([hs.ref[k] if k in hs.ref else hs.acc.get((ref_f4, k), np.zeros(n)) for k in (1, 2, 3)],
                        axis=1)
        ref_txt = f"absolute acceleration of node {refnode} (" + ("FILE8" if setup else "MOTION .ACC") + ")"
    else:
        ag, ref_txt = ground_motion(d, plan.ssidir, info.model, dt, n, setup)
        aref = np.outer(ag * info.gravity, e_cm)
        if str(d["groundfile"] or "").strip() and disp_req:
            _compare_ground(plan, ag, dt, n, setup, lst)
    with _quiet():
        A_ans = aref @ (R_2D.T if info.two_d else np.eye(3))   # ANSYS components of the frame acceleration
    lst.section("Inputs")
    lst.write(f"   histories                : {hs.origin.get('histories', '')}")
    lst.write(f"   time grid                : dt {dt:g} s, {n} samples ({(n - 1) * dt:g} s)")
    lst.write(f"   control direction        : cm {cm}, ang {ang:g} deg -> ({e_cm[0]:.4g}, {e_cm[1]:.4g}, {e_cm[2]:.4g}) "
              f"from {csrc}")
    lst.write(f"   ACEL (reference) motion  : {ref_txt}; peak {np.max(np.abs(aref)) / info.gravity if n else 0:.6g} g")
    missing: List[int] = []
    # ---- consistency checks
    if method == "ACC":
        devs = []
        for nd in plan.dnodes:
            for k in (1, 2, 3):
                if (nd, k) in hs.acc and k in info.dofs.get(nd, []):
                    devs.append(np.max(np.abs(hs.acc[(nd, k)] - aref[:, k - 1])))
        if devs:
            peak = float(np.max(np.abs(aref))) or 1.0
            lst.write(f"   interface nodes vs reference: max |a_i - a_ref| = {max(devs) / info.gravity:.4g} g, "
                      f"{max(devs) / peak:.2%} of the peak (the fixed base assumes a rigid foundation)")
    _rest(plan, hs, lst)
    # ---- tables
    digits = plan.digits
    tables: List[str] = []
    ntab = 0
    acel_fields: List[str] = []
    for ia, lab in enumerate(("X", "Y", "Z")):
        col = A_ans[:, ia]
        if np.any(col != 0.0):
            name = f"LG_AC{lab}"
            tables += table_lines(name, col, digits, f"reference-frame acceleration {lab} (length/s^2)")
            acel_fields.append(f"%{name}%")
            ntab += 1
        else:
            acel_fields.append("0")
    d_cmds: List[str] = []
    for nd in plan.dnodes:
        a = _ansys_node(plan, nd, missing)
        if a is None:
            continue
        if method == "REL":
            for k in _d_dofs(plan, nd, rotdisp):
                ia, sg = ansys_dof(k, info.two_d)
                name = apdl_name("LGD", a, DISP_LABELS[ia])
                tables += table_lines(name, sg * hs.disp[(nd, k)], digits,
                                      f"node {a} {DISP_LABELS[ia]} relative to the reference motion")
                d_cmds.append(f"D,{a},{DISP_LABELS[ia]},%{name}%")
                ntab += 1
        else:
            for k in _d_dofs(plan, nd, True):
                ia, _ = ansys_dof(k, info.two_d)
                d_cmds.append(f"D,{a},{DISP_LABELS[ia]},0")
    chk: List[str] = []
    for nd in check:
        a = _ansys_node(plan, nd, missing)
        if a is None:
            continue
        for k in info.dofs.get(nd, []):
            if k <= 3 or rotacc:
                ia, sg = ansys_dof(k, info.two_d)
                name = apdl_name("LGA", a, DISP_LABELS[ia])
                chk += table_lines(name, sg * hs.acc[(nd, k)], digits,
                                   f"SSI absolute acceleration of node {a} {DISP_LABELS[ia]} (check only, not a load)")
                ntab += 1
    _check_missing_map(plan, missing, "interface or check", lst)
    if ntab > TABLE_WARN:
        lst.warning(f"{ntab} TABLE parameters: ANSYS limits the number of parameters (about 5000 in many releases); "
                    "reduce the interface nodes or split the model")
    if method == "REL" and plan.dnodes and not rotdisp and any(k >= 4 for nd in plan.dnodes
                                                                for k in info.dofs.get(nd, [])):
        lst.warning("interface nodes have rotational DOFs whose rotations are not prescribed (LOADGENDYN <rotdisp> 0)")
    alpha, beta = float(d["alpha"]), float(d["beta"])
    if alpha == 0.0 and beta == 0.0:
        lst.warning("no Rayleigh damping (alpha = beta = 0): give ALPHAD/BETAD (LOADGENDYN) for a realistic transient")
    lst.section("ANSYS transient loads")
    nchk = sum(1 for c in chk if c.startswith("*DIM"))
    lst.write(f"   {ntab} TABLE arrays of {n} values; ACEL,{','.join(acel_fields)}; {len(d_cmds)} D commands"
              + (f"; {nchk} check tables" if nchk else ""))
    if disp_req:
        peak = max(hs.disp.items(), key=lambda kv: float(np.max(np.abs(kv[1]))))
        lst.write(f"   largest relative displacement {np.max(np.abs(peak[1])):.6e} at node {_model_node(plan, peak[0][0])} "
                  f"{C.DOF_TAGS[peak[0][1]]}")
    if alpha or beta:
        for f in (1.0, 2.0, 5.0, 10.0, 20.0):
            w = 2 * math.pi * f
            lst.write(f"   Rayleigh damping ratio at {f:5.1f} Hz: {alpha / (2 * w) + beta * w / 2:.4f}")
    if int(d["opmode"]):
        lst.section("Data check")
        lst.write("   input checked; data-check mode writes no APDL file")
        return
    path = _apdl_path(plan, True)
    text = apdl_dynamic(plan, method, dt, n, tables, chk, acel_fields, d_cmds, ref_txt, alpha, beta, hs)
    path.write_text(text, encoding="utf-8")
    plan.files.append(f"{path.name} (transient: {ntab} tables, {len(d_cmds)} D, ACEL {','.join(acel_fields)})")


def apdl_dynamic(plan: Plan, method: str, dt: float, n: int, tables: List[str], chk: List[str], acel: List[str],
                 d_cmds: List[str], ref_txt: str, alpha: float, beta: float, hs: HistorySet) -> str:
    """APDL of the dynamic second step (D-LGN-04, D-LGN-08)."""
    info = plan.info
    tend = (n - 1) * dt
    items: List[Tuple[str, str]] = [
        ("^", f"{PRODUCT} {__version__}  LOADGEN (Option A): ANSYS dynamic load, full transient (direct integration)"),
        ("Model", info.model + (f"  ({info.title})" if info.title else "")),
        ("Written", f"{time.strftime('%Y-%m-%d %H:%M:%S')}  (LOADGEN listing {info.model}_LOADGEN.out)"),
        ("Histories", hs.origin.get("histories", "")),
        ("Time", f"{n} samples, dt = {dt:g} s, 0 .. {tend:g} s (TABLE arrays, primary variable TIME)"),
    ]
    if method == "REL":
        items += [("Method", f"REL -- relative formulation (manual 6.4.15): the reference frame moves with the ground "
                             f"motion, ACEL = {ref_txt}; the interface nodes follow their SSI displacements RELATIVE to "
                             f"that motion (D tables; {plan.dnode_src}).  ANSYS displacements and accelerations are "
                             "relative to the frame (add the ACEL table for absolute accelerations); the stresses are "
                             "those of the total motion, a rigid translation producing none.")]
    else:
        items += [("Method", f"ACC -- fixed base driven by the SSI foundation motion: ACEL = {ref_txt}; interface "
                             f"nodes fixed (D = 0; {plan.dnode_src}).  Valid for a rigid foundation.")]
    items += [
        ("Damping", f"Rayleigh ALPHAD = {alpha:g}, BETAD = {beta:g} (zeta(w) = alpha/(2 w) + beta w/2).  The SSI "
                    "analysis used frequency-independent hysteretic damping: Rayleigh damping is the main "
                    "approximation of the dynamic 2nd step (manual warning)."),
        ("Units", _units_line(info)),
        ("Axes", "global Cartesian axes of the SASSI-EDU model (ANSYS CSYS,0), Z up"
         + ("; 2-D model: SASSI X-Z plane -> ANSYS X-Y plane (Y = SASSI Z)" if info.two_d else "")),
        ("Nodes", f"{plan.nmap.mode}: {plan.nmap.source}"),
        ("Mass", "the ANSYS model must carry the mass of the structure (DENS, MASS21): ACEL loads it with -M a"),
        ("How to run", "read your structural model, then this file:"),
        (">", "/INPUT,<model>,inp       ! the ANSYS model (e.g. the output of the SASSI-EDU ANSYS command)"),
        (">", f"/INPUT,{_apdl_path(plan, True).stem},inp   ! this file: tables, loads and the transient solution"),
        ("", "With LG_SOLVE = 0 (edit below) everything is defined but not solved: add your own options (NLGEOM, "
             "OUTRES, contact ...) and SOLVE."),
    ]
    head = _head(items)
    L = _comment_block(head)
    L += [f"LG_SOLVE = 1                ! 1 = SOLVE at the end of this file, 0 = define the loads only",
          f"LG_DT = {_exact(dt)}          ! time step of the SSI histories (s)",
          f"LG_NT = {n}                  ! samples per table",
          f"LG_TEND = {_exact(tend)}        ! end time (s) = (LG_NT - 1) LG_DT"]
    if tables:
        L += ["!", "! ---- TABLE arrays: column 0 = time (s), column 1 = value"] + tables
    if chk:
        L += ["!", "! ---- SSI absolute accelerations for comparison with the ANSYS response (not loads)"] + chk
    L += ["!", "/SOLU", "ANTYPE,TRANS", "TRNOPT,FULL", "TIMINT,ON", "AUTOTS,OFF", "KBC,0",
          "DELTIM,LG_DT                ! constant time step = SSI time step (tables are evaluated at each substep)",
          "TIME,LG_TEND",
          f"ALPHAD,{_exact(alpha)}", f"BETAD,{_exact(beta)}",
          "OUTRES,ALL,ALL              ! every substep (reduce for large models, e.g. OUTRES,ALL,10)",
          f"ACEL,{','.join(acel)}" + ("            ! acceleration of the reference frame = ground motion" if any(
              a != "0" for a in acel) else "")]
    if d_cmds:
        L += ["!", "! ---- interface nodes"] + d_cmds
    L += ["!", "*IF,LG_SOLVE,EQ,1,THEN", "SOLVE", "*ENDIF", "FINISH", ""]
    return "\n".join(L)


# =======================================================================================
# Runner (equivalent of sassi.modules.base.run_module) and batch protocol
# =======================================================================================
def run_loadgen(model: str, workdir: PathLike, deck_path: Optional[PathLike] = None,
                listing_path: Optional[PathLike] = None, echo: Optional[Callable[[str], None]] = None,
                progress: Optional[Callable[[float, str], None]] = None,
                cancelled: Optional[Callable[[], bool]] = None) -> int:
    """Run LOADGEN for ``model`` in ``workdir`` like :func:`sassi.modules.base.run_module` (listing
    ``<model>_LOADGEN.out``, status line, 0 on success, 1 module error, 2 unexpected failure)."""
    workdir = Path(workdir)
    deck_path = Path(deck_path) if deck_path is not None else workdir / f"{model}{DECK_EXT}"
    listing = Listing(Path(listing_path) if listing_path is not None else workdir / f"{model}_{NAME}.out", echo=echo)
    ctx = ModuleContext(module=NAME, model=model, workdir=workdir, deck_path=deck_path, listing=listing,
                        progress=progress or (lambda f, m: None), cancelled=cancelled or (lambda: False))
    t0 = time.time()
    rc = 2
    try:
        listing.header(NAME, "ACS SASSI-ANSYS integration, Option A (ANSYS Eq. Static Load / ANSYS Dynamic Load)")
        rc = int(run(ctx) or 0)
    except ModuleError as exc:
        listing.error(str(exc))
        rc = 1
    except Exception as exc:  # unexpected: keep the traceback in the listing
        listing.error(f"unexpected failure: {exc!r}")
        listing.write(traceback.format_exc())
        rc = 2
    finally:
        listing.write("")
        listing.write(f" {NAME} finished with status {'OK' if rc == 0 else 'FAILED'} in {time.time() - t0:.2f} s; "
                      f"{len(listing.warnings)} warning(s), {len(listing.errors)} error(s)")
        listing.close()
    return rc


def batch_main(stdin: Optional[TextIO] = None) -> int:
    """``python -m sassi.modules.loadgen < LOADGEN.inp`` (model, deck, listing; manual section 3.2)."""
    stdin = stdin or sys.stdin
    lines = [ln.strip() for ln in stdin.read().splitlines() if ln.strip()]
    if len(lines) < 3:
        print(f"{NAME}: expected three input lines (model, deck, listing)", file=sys.stderr)
        return 2
    model, deck, out = lines[:3]
    wd = Path.cwd()
    return run_loadgen(model, wd, deck_path=wd / deck, listing_path=wd / out, echo=print)


if __name__ == "__main__":
    raise SystemExit(batch_main())
