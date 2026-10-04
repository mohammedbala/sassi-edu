"""ANALYS module: soil impedance and flexible-volume SSI solution -> FILE8 (requirements 2.1, 2.4-2.6,
4.6; R1 section 4; decisions D-ANL-01..13, D-W1-13; spec 05c Part A).

What ANALYS does (for the structural engineer)
----------------------------------------------
For every SSI frequency f_q = n_q df (w = 2 pi f_q), in the order of the frequency list (q = 1..nF):

1. **Frequency survey** (D-ANL-08): the frequency numbers to solve (the deck frequency set, or all
   frequencies of FILE1/FILE9 when <fopt> = 1) must exist in the free-field / load files and in
   FILE3 (or in the restart files of a restart run); otherwise the run stops before any solve.
2. **Soil flexibility** ``F_ff`` (3 nInt x 3 nInt) at the interaction nodes from the POINT solution
   (FILE3; :func:`sassi.core.flexibility.flexibility_matrix`, R1 4.1, D-ANL-02).
3. **Soil impedance** ``X_ff = F_ff^-1`` (dense complex LU, D-ANL-04).
4. **System** ``C = (K*_s - w^2 M_s) - (K*_e - w^2 M_e) + A_f^T X_ff A_f`` (general assembly of the
   flexible-volume equation, R1 4.4, D-ANL-01).  The method variants FV / FI-FSIN / FI-EVBN / FFV
   differ *only* by the set of interaction nodes chosen in HOUSE (R1 4.5).
5. **Load** -- seismic: ``b = A_f^T X_ff U'_f`` with the free field ``U'_f`` of FILE1 at each
   interaction node (:func:`sassi.core.freefield.free_field_at_nodes`: SITE axes rotated by the
   coordinate transformation angle <ang>, phase ``exp(-i k x')`` of inclined/surface waves with the
   control point (xc, yc) as origin); foundation vibration: the FORCE load vector of FILE9.
   Incoherent motion, wave passage and multiple excitation (HOUSE <coh>, <wpass>, <me> = 1;
   requirements 4.6 item 5, spec 05c A.5.9): the HOUSE factors ``s`` of FILE77 (one complex factor
   per interaction node, frequency and direction; the factor of the control-motion direction multiplies
   the three components of a node) act on the free-field **load**, ``b = A_f^T (s .* X_ff U'_f)``
   (FFL, ANALYSX <ffm> = 0, the default validated by EPRI), or on the free-field **motion**,
   ``b = A_f^T X_ff (s .* U'_f)`` (FFM, <ffm> = 1; surface foundations only, D-INC-12)
   (:func:`sassi.core.ssi_solver.incoherent_seismic_load`).
6. **Solve** by the Schur complement on the interaction DOFs (sparse LU of the non-interaction
   block, dense LU of ``S + X_ff``; full sparse LU fallback with a warning), see
   :mod:`sassi.core.ssi_solver` (D-ANL-03).  All simultaneous load cases share the factorisations.
7. **Output**: the transfer functions of every equation of FILE4 to FILE8 (seismic: total motion
   per unit control motion; vibration: displacement per unit load factor; fixed DOFs are not
   equations, H = 0, D-ANL-13).
8. **Restart files** (<save> = 1): ``COOXqqq`` = X_ff, ``COOTKqqq`` = the factorised system, both
   carrying the frequency number, the frequency value and the FILE90 hashes; ``COOXI`` / ``COOTKI``
   (with the frequency step ``df``) index them by frequency number, so restart runs may use any subset
   of the saved frequencies (D-ANL-05).
9. **Global impedance** (<impe> 1/2): ``K_G = T^T X_ff T`` about the control point -> FOUNSTIF,
   FOUNDASH, FOUNDAMP, FOUNIMPD and FILE11 (D-ANL-07).

Modes of analysis (requirements 2.5, spec 05c A.5.3)
---------------------------------------------------
====  =========================  ===========  =============================================
deck  manual name                solver mode  reuses
====  =========================  ===========  =============================================
0     Initiation                 1            nothing (FILE3 needed)
1     New Structure              2            X_ff (COOXqqq); structure re-assembled
2     New Seismic Environment    3            X_ff and the factorised system (COOTKqqq)
3     New Dynamic Loading        3            X_ff and the factorised system (COOTKqqq)
====  =========================  ===========  =============================================

A restart record is accepted only when its FILE90 hashes match the current FILE90 of HOUSE:
interaction nodes and site layering for COOX (``int_hash``, ``layer_hash``) and, in addition, the
FILE4/COOSK/COOSM content for COOTK (``file4_hash``) -- and only when it belongs to the same
frequency: the frequency step of the index and the frequency stored in every record must equal
those of the run (relative 1e-6, as the D-ANL-08 step check), otherwise frequency number n would
denote another frequency (a new NFFT or time step in SITE/FORCE).  Every restart result equals a
fresh initiation run to round-off (VP-22).

Simultaneous cases (D-ANL-06)
-----------------------------
* seismic, <simul> = 0: FILE1 -> FILE8;  <simul> = 1: FILE1X, FILE1Y, FILE1Z -> FILE8X, FILE8Y,
  FILE8Z (X + Y + Z, the angle must be 0).  With a deterministic FILE77 (AS sum, single incoherent
  mode, wave passage or multiple excitation) the same files are produced, the FILE77 factors of the
  control direction (X, Y, Z) multiplying the load of each case;
* seismic stochastic simulation, <simul> = Ns (1..50, the HOUSEX <nsim> of HOUSE) with the samples
  FILE77001 .. FILE77Ns: FILE1X/Y/Z -> FILE8001 .. FILE8(3 Ns), simulation s and direction d
  (1 X, 2 Y, 3 Z) in ``FILE8{3(s-1)+d:03d}`` (requirements 1.8, spec 05c A.5.4).  <simul> = 1 is a
  stochastic run when only stochastic samples are present (FILE77001, no FILE77);
* vibration, <simul> = 0 or 1: FILE9 -> FILE8;  <simul> = Nl >= 2: FILE9001 .. FILE9Nl ->
  FILE8001 .. FILE8Nl (Nl <= 500).  (The single-case rule for <simul> = 1 is the rule of the RUN
  prerequisite check, sassi/prep/commands/modules_cmd.py.)

Listing
-------
Input summary, frequency survey table, resource check, a per-frequency summary (solver path, sizes,
condition estimates, time), the low-frequency rigid-body check (G-19 / EDU-18), the transfer
functions per <prnt> (amplitudes, or real and imaginary parts, of all six DOFs of at most
:data:`MAX_PRINT_NODES` nodes -- the lowest node numbers; use MOTION for the others) and the global
impedance tables.

2D analysis (HOUSE <dim> = 1; requirements 4.3, 4.6 item 2, D-PNT-01, D-W2-07; :mod:`sassi.core.ssi2d`)
-----------------------------------------------------------------------------------------------------
A plane-strain model in the X-Z plane (PLANE elements, BEAMS/SPRING/GENERAL in the plane), per unit
length.  FILE3 must hold POINT2 line-load solutions (an error otherwise); the soil flexibility is the
2x2 P-SV block per node pair with the odd terms times sgn(x_i - x_j) and the impedance is the inverse of
the in-plane (UX, UZ) block, stored in the 3D layout with zero UY rows and columns.  Seismic input must
be in-plane (control direction x' or z': SV, P, Rayleigh waves) with the angle 0 or 180 deg; the
free field is evaluated at (x, yc, z) (the y coordinate is ignored).  <simul> = 1 analyses the X and Z
cases (FILE1X, FILE1Z -> FILE8X, FILE8Z; the anti-plane Y case is skipped).  The global impedance
(<impe>) and incoherent / wave-passage / multiple-excitation input are 3D only (errors).

Symmetry and antisymmetry planes (SYMM, D-ANL-12; :mod:`sassi.core.symmetry`)
----------------------------------------------------------------------------
HOUSE stores the planes of a half / quarter model in FILE4 (``x_symm``) and fixes the symmetry DOFs
of the plane nodes (``x_int_symm`` flags the constrained interaction translations).  The soil
flexibility of the reduced interaction set is the image sum ``F_red(i, j) = sum_S s_S F(i, j_S) P_S``
over the 2^k combinations of the k planes; the constrained translations are removed before
``X_red = F_red^-1``.  The seismic free field must have the symmetry of every plane (``U'(image) =
s P U'``, checked at every frequency and case: an X input is antisymmetric about a plane normal to X
and symmetric about a plane normal to Y, a Z input symmetric about both) -- the run stops otherwise.
Incoherent motion, wave passage, multiple excitation and the global impedance are not allowed with
SYMM (manual 2.7 and 6.5.4, G-17).  Vibration loads on plane nodes are those of the reduced model
(one half / one quarter of the full-model load, like the masses).

Node-numbering optimizer (HOUSEX <optimize> = 1, D-HOU-03)
---------------------------------------------------------
FILE4 is then in the optimised numbering (meta ``x_optimized`` = 1) but FORCE writes FILE9 with the
node numbers of the model (the interpreter numbering of ``<model>.hou``).  ANALYS translates every
FILE9 load node through the old -> new pairs of FILE4 (``x_node_old_id``, or ``<model>.map``) and lists
the translation; without the pairs the run stops (applying the loads to whatever node now carries the
old number would be silently wrong).

    python -m sassi.modules.analys < ANALYS.inp      # batch use
"""
from __future__ import annotations

import contextlib
import contextvars
import os
import re
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..conventions import DOF_LABELS, frequency_step, incoherent_file8_name, loadcase_name, restart_names
from ..core import incoherency as INC
from ..core import renumber as RN
from ..core import ssi2d as S2
from ..core import ssi_solver as ss
from ..core import symmetry as SYM
from ..core.flexibility import flexibility_block, flexibility_matrix
from ..core.freefield import free_field_at_nodes
from ..core.interp import rigid_body_anchor
from ..io import decks
from ..io.container import Container, read_container
from ..io.files import validate, write_container
from .base import ModuleContext, ModuleError, batch_main

NAME = "ANALYS"

#: maximum number of SSI frequencies of one ANALYS run (requirements 4.0.1)
MAX_FREQ = 500
#: maximum number of simultaneous vibration load cases (spec 05c A.5.4)
MAX_LOAD_CASES = 500
#: maximum number of incoherent stochastic simulations (spec 05c A.5.4, requirements 4.4 item 6)
MAX_SIMULATIONS = INC.MAX_SIM
#: the listing prints the transfer functions of at most this many nodes (lowest node numbers)
MAX_PRINT_NODES = 50
#: relative tolerance of the frequency-step comparison between the deck and the input files (D-ANL-08)
DF_RTOL = 1e-6
#: refuse the run when the memory estimate exceeds this fraction of the physical RAM (D-ANL-10)
MEMLIMIT_FRACTION = 0.8
#: assumed fill-in factor of the sparse LU (memory estimate only)
SPARSE_FILL = 30
#: dense (3 nInt)^2 complex matrices of the memory estimate (D-ANL-10): X_ff and the LU of S + X_ff
#: (2) plus one for the flexibility coordinate tables and the restart-file buffers
DENSE_COPIES = 3.0
#: ... with <save> = 1: a record (1) plus its in-memory write buffer (WRITE_COPIES) may coexist
DENSE_COPIES_SAVE = 3.25
#: size of the byte buffer of a container write relative to the array written (buffer growth included)
WRITE_COPIES = 2.25
#: complex work arrays of the flexibility assembly blocks (flexibility_matrix, 3D)
FLEX_WORK_ARRAYS = 8
#: dense copies of a frequency on the full sparse LU fallback path (X_ff, the scattered block, the
#: assembled matrix and the LU of the dense block)
FALLBACK_COPIES = 5.0
#: node pairs per block of the flexibility assembly (flexibility_matrix ``block`` = this // nInt rows;
#: its default): about 8 complex work arrays of this size are alive while F_ff is assembled
FLEX_BLOCK_ELEMS = 2_000_000
#: FILE11 holds X_ff of every frequency only while that takes at most this many bytes
FILE11_MAX_X_BYTES = 256 * 2 ** 20
#: context-local override of FILE11_MAX_X_BYTES (:func:`file11_x_cap`); None = the module default
_FILE11_X_CAP = contextvars.ContextVar("analys_file11_x_cap", default=None)


def file11_x_limit() -> int:
    """Byte limit for X_ff in FILE11 in the current thread / context (D-ANL-07): the value set by an
    enclosing :func:`file11_x_cap`, else :data:`FILE11_MAX_X_BYTES`."""
    cap = _FILE11_X_CAP.get()
    return FILE11_MAX_X_BYTES if cap is None else int(cap)


@contextlib.contextmanager
def file11_x_cap(nbytes: int):
    """Inside the ``with`` block, ANALYS runs of the *calling* thread (context variable, not a module
    global) store X_ff in FILE11 only while it takes at most ``nbytes``; 0 = K_G and T only.  ANALYS
    runs in other threads keep the default (D-ANL-07: FILE11 holds X_ff and T)."""
    token = _FILE11_X_CAP.set(max(0, int(nbytes)))
    try:
        yield
    finally:
        _FILE11_X_CAP.reset(token)

#: low-frequency rigid-body check of the seismic TFs (G-19 / EDU-18)
LOWFREQ_TOL = 0.05

SOLVER_MODE = {0: 1, 1: 2, 2: 3, 3: 3}
MODE_NAMES = {0: "Initiation", 1: "New Structure", 2: "New Seismic Environment", 3: "New Dynamic Loading"}
METHOD_NAMES = {0: "FV (flexible volume)", 1: "FFV (fast flexible volume)", 2: "FI (flexible interface)"}
IMPE_NAMES = {0: "none", 1: "decoupled (diagonal) impedances", 2: "full 6x6 rigid-body impedance matrix"}
RIGID_DOFS = ("X", "Y", "Z", "XX", "YY", "ZZ")
PATH_CODE = {p: k for k, p in enumerate(ss.PATHS)}
#: restart file names: COOXqqq / COOTKqqq (3 or more digits: order numbers >= 1000 give 4 digits) and the
#: indexes COOXI / COOTKI -- not COOSK / COOSM of HOUSE
_RESTART_RE = re.compile(r"^COO(X|TK)(\d{3,}|I)$")


# ======================================================================================
# Input data
# ======================================================================================
@dataclass
class LoadCase:
    """One simultaneous case: input file (FILE1x / FILE9lll) -> output file (FILE8x / FILE8lll)."""

    label: object                  # FILE8 meta 'case': 0 single, 'X'/'Y'/'Z', load case number, 'S001X' ...
    infile: str
    outfile: str
    data: Optional[Container] = None
    rows: Dict[int, int] = field(default_factory=dict)       # frequency number -> row
    cm: int = 0
    # incoherent / wave-passage / multiple-excitation factors (requirements 4.6 item 5)
    f77: str = ""                  # FILE77 or FILE77sss ('' = none)
    f77_dir: int = -1              # direction row of the factors (0 X, 1 Y, 2 Z); -1 = control direction cm
    sim: int = 0                   # stochastic simulation (1-based), 0 = deterministic
    lowfreq_skip: str = ""         # reason to skip the G-19 rigid-body check ('' = check)


@dataclass
class Options:
    """Resolved and checked ANALYS deck options."""

    title: str
    opmode: int
    type: int
    mode: int
    save: int
    prnt: int
    fopt: int
    ang: float
    xc: float
    yc: float
    zc: float
    impe: int
    simul: int
    delrst: int
    df: float
    delt: float
    nft: int
    gravity: float
    freqs: List[int]
    freqset: int = 1
    coh: int = 0
    wpass: int = 0
    me: int = 0
    ffm: int = 0

    @property
    def seismic(self) -> bool:
        return self.type == 0

    @property
    def factored(self) -> bool:
        """Seismic input with HOUSE incoherency factors (FILE77): incoherent motion, wave passage or
        multiple excitation (requirements 4.6 item 5)."""
        return self.type == 0 and bool(self.coh or self.wpass or self.me)

    @property
    def solver_mode(self) -> int:
        return SOLVER_MODE[self.mode]


@dataclass
class ModelData:
    """FILE4 (``<model>.N4``) and the matrices / hashes of the HOUSE run."""

    file4: Container
    neq: int
    eq_node: np.ndarray
    eq_dof: np.ndarray
    int_node: np.ndarray
    int_iface: np.ndarray
    int_eq: np.ndarray
    int_xyz: np.ndarray
    part: ss.Partition
    model_hash: str
    method: int
    Ks: object = None
    Ms: object = None
    Ke: object = None
    Me: object = None
    hashes: Dict[str, str] = field(default_factory=dict)
    eqmap: Dict[Tuple[int, int], int] = field(default_factory=dict)     # (node, dof) -> equation
    dim: int = 2                                   # HOUSE <dim>: 1 2D (POINT2), 2 3D (POINT3)
    tol: float = 0.0                               # geometric tolerance of HOUSE (D-GEN-06)
    symm: List[SYM.SymmetryPlane] = field(default_factory=list)          # SYMM planes (FILE4 x_symm)
    int_symm: Optional[np.ndarray] = None          # (nInt, 3) translations constrained by the planes
    node_map: Optional[Dict[int, int]] = None      # optimizer: model (old) node number -> FILE4 number
    map_source: str = ""                           # where node_map comes from (listing)

    @property
    def nint(self) -> int:
        return int(self.int_node.size)

    @property
    def reduced(self) -> bool:
        """The impedance is computed on a subset of the 3 nInt translations (2D or SYMM)."""
        return self.dim == 1 or bool(self.symm)

    def keep_mask(self) -> np.ndarray:
        """(3 nInt,) interaction translations of the soil impedance (2D: UX, UZ; SYMM: minus the
        constrained translations of plane nodes)."""
        return SYM.keep_mask(self.dim, self.int_symm, self.nint)


@dataclass
class FreqRecord:
    """Per-frequency summary row of the listing."""

    q: int
    fnum: int
    f: float
    path: str = ""
    nf: int = 0
    nn: int = 0
    rcond_F: float = float("nan")
    rcond_nn: float = float("nan")
    rcond_ff: float = float("nan")
    seconds: float = 0.0
    note: str = ""


# ======================================================================================
# Deck options
# ======================================================================================
def read_options(d: decks.Deck, lst) -> Options:
    """Check the deck options (spec 05c A.5, A.7) and stop on tier P1/P2 requests."""
    def ival(name: str, allowed: Sequence[int]) -> int:
        v = int(d[name])
        if v not in allowed:
            raise ModuleError(f"<{name}> = {v} must be one of {list(allowed)}")
        return v

    opmode = ival("opmode", (0, 1))
    typ = ival("type", (0, 1))
    mode = int(d["mode"])
    if mode == 6:
        raise ModuleError("Mode 6 (New Load Vector) is not available in this build (D-ANL-09, out of scope)")
    if mode not in MODE_NAMES:
        raise ModuleError(f"<mode> = {mode} must be 0 (initiation), 1 (new structure), 2 (new seismic "
                          "environment) or 3 (new dynamic loading)")
    save = ival("save", (0, 1))
    prnt = ival("prnt", (0, 1))
    fopt = ival("fopt", (0, 1))
    impe = ival("impe", (0, 1, 2))
    simul = int(d["simul"])
    if simul < 0:
        raise ModuleError(f"<simul> = {simul} must be >= 0")
    coh, wpass, me = ival("coh", (0, 1)), ival("wpass", (0, 1)), ival("me", (0, 1))
    ffm = ival("ffm", (0, 1))
    if typ == 0 and simul > 1 and not coh:
        raise ModuleError(f"<simul> = {simul} for seismic analysis means {simul} incoherent stochastic simulations "
                          "(FILE77001 ...): it needs incoherent motion (HOUSE <coh> = 1); the coherent X/Y/Z cases "
                          "use <simul> = 1")
    if typ == 0 and coh and simul > MAX_SIMULATIONS:
        raise ModuleError(f"<simul> = {simul} incoherent simulations exceed the maximum {MAX_SIMULATIONS} "
                          "(spec 05c A.5.4)")
    if typ == 1 and simul > MAX_LOAD_CASES:
        raise ModuleError(f"<simul> = {simul} load cases exceed the maximum {MAX_LOAD_CASES} (spec 05c A.5.4)")
    if typ == 1 and (coh or wpass or me):
        lst.warning("incoherency, wave passage and multiple excitation (HOUSE <coh>/<wpass>/<me>) apply to seismic "
                    "input only: ignored for foundation vibration")
        coh = wpass = me = 0
    if ffm and not (typ == 0 and (coh or wpass or me)):
        lst.warning("free-field motion (FFM, ANALYSX <ffm> = 1) applies to incoherent, wave-passage or "
                    "multiple-excitation input only: ignored")
        ffm = 0
    ang = float(d["ang"])
    if not 0.0 <= ang < 360.0:
        a2 = ang % 360.0
        lst.warning(f"coordinate transformation angle {ang:g} deg normalised to {a2:g} deg (D-ANL-11)")
        ang = a2
    if typ == 0 and simul >= 1 and ang != 0.0:
        raise ModuleError(f"simultaneous X/Y/Z seismic cases need the coordinate transformation angle 0 "
                          f"(<ang> = {ang:g}; D-ANL-06, spec 05c A.5.4)")
    if typ == 1 and mode == 2:
        lst.warning("New Seismic Environment (<mode> = 2) with foundation vibration: solved as New Dynamic "
                    "Loading (same solver Mode 3)")
    if typ == 0 and mode == 3:
        lst.warning("New Dynamic Loading (<mode> = 3) with seismic analysis: solved as New Seismic "
                    "Environment (same solver Mode 3)")
    delrst = int(d["delrst"])
    delt, nft = float(d["delt"]), int(d["nft"])
    df = float(d["df"])
    if df <= 0:
        try:
            df = frequency_step(delt, nft, 0.0)
        except ValueError:
            df = 0.0                      # taken from the input files
    freqs = [int(r["number"]) for r in d.rows("freqs")]
    if len(set(freqs)) != len(freqs):
        raise ModuleError("the frequency set contains duplicate frequency numbers")
    if freqs and min(freqs) < 1:
        raise ModuleError("frequency numbers must be positive integers")
    return Options(title=str(d["title"]), opmode=opmode, type=typ, mode=mode, save=save, prnt=prnt, fopt=fopt,
                   ang=ang, xc=float(d["xc"]), yc=float(d["yc"]), zc=float(d["zc"]), impe=impe, simul=simul,
                   delrst=delrst, df=df, delt=delt, nft=nft, gravity=float(d["gravity"]), freqs=sorted(freqs),
                   freqset=int(d["freq"]), coh=coh, wpass=wpass, me=me, ffm=ffm)


# ======================================================================================
# Input files
# ======================================================================================
def _read(ctx: ModuleContext, name: str, kind: str, producer: str) -> Container:
    p = ctx.require(name, producer)
    try:
        c = read_container(p, kind=kind)
    except (ValueError, OSError) as exc:
        raise ModuleError(f"{name}: {exc}") from None
    probs = validate(c)
    if probs:
        raise ModuleError(f"{name} does not satisfy the {kind} schema: " + "; ".join(probs))
    return c


def read_model(ctx: ModuleContext, opt: Options, lst) -> ModelData:
    """FILE4 family (requirements 2.1): FILE4 always; COOSK/COOSM when the system is assembled
    (solver modes 1 and 2); FILE90 for saved or reused restart files; DOFSMAP and FILE91 for restarts."""
    name4 = f"{ctx.model}.N4"
    f4 = _read(ctx, name4, "FILE4", "HOUSE")
    dim = int(f4.meta["dim"])
    if dim not in (1, 2):
        raise ModuleError(f"{name4}: <dim> = {dim}: ANALYS needs a 2D (<dim> = 1) or 3D (<dim> = 2) model "
                          "(1D analysis is not available)")
    if dim == 1 and (("int_dofs" in f4.meta and [int(v) for v in f4.meta["int_dofs"]] != list(S2.IN_PLANE_DOFS))
                     or np.any(np.asarray(f4["int_eq"], np.int64).reshape(-1, 3)[:, 1] >= 0)):
        raise ModuleError(f"{name4}: a 2D model has the interaction DOFs UX, UZ (int_dofs [1, 3], UY column of int_eq "
                          f"-1; D-W2-07), the file says {f4.meta.get('int_dofs')} -- re-run HOUSE")
    neq = int(f4.meta["nEq"])
    eq_node = np.asarray(f4["eq_node"], np.int64)
    eq_dof = np.asarray(f4["eq_dof"], np.int64)
    if eq_node.size != neq:
        raise ModuleError(f"{name4}: {eq_node.size} equations listed, meta nEq = {neq}")
    int_node = np.asarray(f4["int_node"], np.int64)
    int_iface = np.asarray(f4["int_iface"], np.int64)
    int_eq = np.asarray(f4["int_eq"], np.int64).reshape(-1, 3)
    if "x_int_xyz" in f4:
        int_xyz = np.asarray(f4["x_int_xyz"], float).reshape(-1, 3)
    else:
        pos = {int(n): k for k, n in enumerate(np.asarray(f4["node_id"]))}
        xyz = np.asarray(f4["node_xyz"], float)
        int_xyz = np.array([xyz[pos[int(n)]] for n in int_node], float).reshape(-1, 3)
    try:
        part = ss.Partition.from_int_eq(int_eq, neq)
    except ValueError as exc:
        raise ModuleError(f"{name4}: {exc}") from None
    if int_node.size:
        key = S2.duplicate_key(int_xyz, dim)
        if np.unique(key, axis=0).shape[0] != int_node.size:
            raise ModuleError("two interaction nodes have the same coordinates" + (" x and z (2D: y is ignored)"
                                                                                   if dim == 1 else "")
                              + ": the soil flexibility matrix would be singular (merge the nodes)")
    md = ModelData(file4=f4, neq=neq, eq_node=eq_node, eq_dof=eq_dof, int_node=int_node, int_iface=int_iface,
                   int_eq=int_eq, int_xyz=int_xyz, part=part, model_hash=str(f4.meta["model_hash"]),
                   method=int(f4.meta.get("imp", 0)),
                   eqmap={(int(n), int(k)): i for i, (n, k) in enumerate(zip(eq_node, eq_dof))}, dim=dim)
    xyz_all = np.asarray(f4["node_xyz"], float).reshape(-1, 3)
    md.tol = float(f4.meta.get("tol", 0.0) or 0.0) or max(1e-9, 1e-6 * float(np.ptp(xyz_all, axis=0).max())
                                                          if xyz_all.size else 1e-9)
    _symmetry_data(md, f4, name4)
    _node_map(ctx, md, f4, name4)
    if opt.solver_mode in (1, 2):
        ck = _read(ctx, "COOSK", "COOSK", "HOUSE")
        cm = _read(ctx, "COOSM", "COOSM", "HOUSE")
        for c, nm in ((ck, "COOSK"), (cm, "COOSM")):
            if int(c.meta["nEq"]) != neq:
                raise ModuleError(f"{nm} has {c.meta['nEq']} equations, {name4} {neq} -- re-run HOUSE")
            if "model_hash" in c.meta and str(c.meta["model_hash"]) != md.model_hash:
                raise ModuleError(f"{nm} and {name4} come from different HOUSE runs -- re-run HOUSE")
        md.Ks, md.Ke = ck.sparse("Ks"), ck.sparse("Ke")
        md.Ms, md.Me = cm.sparse("Ms"), cm.sparse("Me")
        for nm, A in (("Ks", md.Ks), ("Ke", md.Ke), ("Ms", md.Ms), ("Me", md.Me)):
            if A.shape != (neq, neq):
                raise ModuleError(f"{nm} has shape {A.shape}, expected ({neq}, {neq}) -- re-run HOUSE")
    need90 = opt.save == 1 or opt.mode > 0
    if need90 or ctx.path("FILE90").exists():
        if need90:
            f90 = _read(ctx, "FILE90", "FILE90", "HOUSE")
        else:
            try:
                f90 = read_container(ctx.path("FILE90"), kind="FILE90")
            except (ValueError, OSError):
                f90 = None
        if f90 is not None:
            md.hashes = {k: str(f90.meta[k]) for k in ("int_hash", "layer_hash", "file4_hash") if k in f90.meta}
            if md.hashes.get("file4_hash", md.model_hash) != md.model_hash:
                msg = f"FILE90 and {name4} come from different HOUSE runs"
                if need90:
                    raise ModuleError(msg + " -- re-run HOUSE")
                lst.warning(msg)
    if opt.mode > 0:
        dm = _read(ctx, "DOFSMAP", "DOFSMAP", "HOUSE")
        if not (np.array_equal(np.asarray(dm["eq_node"]), eq_node) and np.array_equal(np.asarray(dm["eq_dof"]), eq_dof)
                and np.array_equal(np.asarray(dm["int_node"]), int_node)):
            raise ModuleError(f"DOFSMAP does not match {name4} -- re-run HOUSE")
        ctx.require("FILE91", "HOUSE")
    return md


def _symmetry_data(md: ModelData, f4: Container, name4: str) -> None:
    """SYMM planes of a half / quarter model (FILE4 ``x_symm``) and the interaction translations they
    constrain (``x_int_symm``, recomputed from the planes when absent); D-ANL-12."""
    if "x_symm" not in f4 or np.asarray(f4["x_symm"]).size == 0:
        return
    md.symm = SYM.planes_from_array(f4["x_symm"])
    if "x_int_symm" in f4:
        cons = np.asarray(f4["x_int_symm"]).reshape(-1, 3) != 0
    else:
        cons = SYM.constrained_translations(md.symm, md.int_xyz, md.tol)
    if cons.shape != md.int_eq.shape:
        raise ModuleError(f"{name4}: x_int_symm has shape {cons.shape}, the interaction set {md.int_eq.shape} -- "
                          "re-run HOUSE")
    bad = np.flatnonzero((cons & (md.int_eq >= 0)).any(axis=1))
    if bad.size:
        raise ModuleError(f"{name4}: interaction nodes {[int(md.int_node[k]) for k in bad[:10]]} lie on a SYMM plane but "
                          "the translations the symmetry constrains are active equations -- re-run HOUSE (it fixes "
                          "them, D-ANL-12)")
    md.int_symm = cons


def _node_map(ctx: ModuleContext, md: ModelData, f4: Container, name4: str) -> None:
    """Old -> new node numbers of an optimised FILE4 (D-HOU-03): FILE4 ``x_node_old_id`` or the
    ``<model>.map`` file named by the meta ``x_node_map``.  Left None when FILE4 is not optimised or
    the pairs are missing (a vibration run then stops, :func:`translate_load_nodes`)."""
    if int(f4.meta.get("x_optimized", 0) or 0) != 1:
        return
    if "x_node_old_id" in f4:
        old = np.asarray(f4["x_node_old_id"], np.int64).reshape(-1)
        new = np.asarray(f4["node_id"], np.int64).reshape(-1)
        if old.size == new.size:
            md.node_map = {int(o): int(n) for o, n in zip(old, new)}
            md.map_source = f"{name4} (x_node_old_id)"
            return
    mp = ctx.path(str(f4.meta.get("x_node_map") or f"{ctx.model}.map"))
    if mp.exists():
        try:
            md.node_map = RN.read_map(mp)
        except (ValueError, OSError) as exc:
            raise ModuleError(f"{mp.name}: {exc}") from None
        md.map_source = mp.name


def translate_load_nodes(md: ModelData, cases: List["LoadCase"], lst) -> None:
    """FORCE writes FILE9 in the model (interpreter) numbering; an optimised FILE4 uses the new numbers
    (D-HOU-03).  Replace every FILE9 load node by its new number (in memory) and list the translation;
    stop when the old -> new pairs are missing or a load node is not in them (Error 62)."""
    if int(md.file4.meta.get("x_optimized", 0) or 0) != 1:
        return
    if md.node_map is None:
        raise ModuleError(f"{md.file4.meta.get('model', 'the model')}.N4 was written by the HOUSE node-numbering "
                          "optimizer (HOUSEX <optimize> = 1) but the old -> new node pairs are missing (FILE4 "
                          "x_node_old_id and the .map file): the FILE9 loads (model node numbers) cannot be placed "
                          "on the renumbered model -- re-run HOUSE")
    done = set()
    for c in cases:
        if id(c.data) in done:
            continue
        done.add(id(c.data))
        nodes = np.asarray(c.data["load_node"], np.int64)
        miss = sorted({int(n) for n in nodes if int(n) not in md.node_map})
        if miss:
            raise ModuleError(f"{c.infile}: loads on nodes {miss[:10]}, which are not in the HOUSE model (Error 62; "
                              f"old -> new pairs of {md.map_source})")
        new = np.array([md.node_map[int(n)] for n in nodes], np.int64)
        c.data.arrays["load_node"] = new
        moved = int(np.sum(new != nodes))
        lst.write(f" {c.infile}: {nodes.size} loaded DOFs, node numbers translated to the optimised numbering of FILE4 "
                  f"({moved} renumbered; old -> new pairs of {md.map_source}, D-HOU-03)")


def load_cases(opt: Options) -> List[LoadCase]:
    """Input / output files of the simultaneous cases (D-ANL-06, spec 05c A.5.4)."""
    if opt.seismic:
        if opt.simul == 0:
            return [LoadCase(0, "FILE1", "FILE8")]
        return [LoadCase(c, f"FILE1{c}", f"FILE8{c}", cm=k) for k, c in enumerate("XYZ")]
    if opt.simul <= 1:
        return [LoadCase(0, "FILE9", "FILE8")]
    return [LoadCase(k, loadcase_name("FILE9", k), loadcase_name("FILE8", k)) for k in range(1, opt.simul + 1)]


# ======================================================================================
# Incoherency, wave passage and multiple excitation (FILE77; requirements 4.6 item 5)
# ======================================================================================
@dataclass
class FactorSet:
    """The FILE77 factors of a run: every file, its frequency rows and the position of every FILE4
    interaction node in it."""

    files: Dict[str, Container] = field(default_factory=dict)
    rows: Dict[str, Dict[int, int]] = field(default_factory=dict)
    perm: Dict[str, np.ndarray] = field(default_factory=dict)

    def factor(self, case: LoadCase, fnum: int) -> np.ndarray:
        """Factors ``s`` (nInt,) of ``case`` at frequency number ``fnum`` in FILE4 interaction order:
        the row of the case direction (X/Y/Z cases) or of the control direction of its FILE1."""
        c = self.files[case.f77]
        d = case.f77_dir if case.f77_dir >= 0 else case.cm
        return np.asarray(c["s"])[self.rows[case.f77][fnum], d][self.perm[case.f77]]


def factored_cases(ctx: ModuleContext, opt: Options, lst) -> List[LoadCase]:
    """Load cases of seismic input with FILE77 factors (D-ANL-06, spec 05c A.5.4).

    Stochastic simulation (incoherent motion with <simul> >= 2, or <simul> = 1 when only stochastic
    samples FILE77001 ... are present): simulation s = 1..<simul>, direction d -> FILE1X/Y/Z with the
    factors of FILE77sss -> ``FILE8{3(s-1)+d:03d}``.  Otherwise the deterministic FILE77 multiplies the
    single case (FILE1 -> FILE8, factors of the control direction) or the X/Y/Z cases (<simul> = 1)."""
    det = ctx.path(INC.file77_name()).exists()
    sto = ctx.path(INC.file77_name(1)).exists()
    if opt.coh and (opt.simul >= 2 or (opt.simul == 1 and sto and not det)):
        cases = [LoadCase(f"S{s:03d}{c}", f"FILE1{c}", incoherent_file8_name(s, k + 1), cm=k,
                          f77=INC.file77_name(s), f77_dir=k, sim=s)
                 for s in range(1, opt.simul + 1) for k, c in enumerate("XYZ")]
        if det:
            lst.warning(f"FILE77 (deterministic factors) is present but not used: this run analyses the {opt.simul} "
                        "stochastic samples FILE77001 ...")
        return cases
    if not det:
        if opt.coh and sto:
            raise ModuleError("FILE77 missing -- the directory holds the stochastic samples FILE77001 ...: set "
                              "<simul> to the number of simulations (with FILE1X/Y/Z) or re-run HOUSE")
        raise ModuleError("FILE77 missing -- run HOUSE (incoherent motion, wave passage or multiple excitation) first")
    if opt.coh and sto:
        lst.warning("stochastic samples FILE77001 ... are present but the deterministic FILE77 is used "
                    f"(<simul> = {opt.simul})")
    cases = load_cases(opt)
    for c in cases:
        c.f77 = INC.file77_name()
        c.f77_dir = c.cm if c.label in ("X", "Y", "Z") else -1
    return cases


def read_factors(ctx: ModuleContext, opt: Options, md: ModelData, cases: List[LoadCase], lst) -> FactorSet:
    """Read the FILE77 files of the cases (HOUSE): schema, frequency step (D-ANL-08), the factors of
    every FILE4 interaction node at its coordinates (a FILE77 of another model is refused), and the
    rigid-body exceptions of the low-frequency check."""
    fs = FactorSet()
    xyz = md.int_xyz
    tol = float(md.file4.meta.get("tol", 0.0)) or max(1e-9, 1e-6 * float(np.abs(xyz).max()) if xyz.size else 1e-9)
    for name in dict.fromkeys(c.f77 for c in cases):
        if not ctx.path(name).exists():
            raise ModuleError(f"{name} missing -- run HOUSE (stochastic simulation with HOUSEX <nsim> >= {opt.simul})")
        try:
            c77 = INC.read_file77(ctx.path(name))
            fs.perm[name] = INC.node_permutation(c77, md.int_node, xyz, max(tol, 1e-9), name)
        except INC.IncoherencyError as exc:
            raise ModuleError(str(exc)) from None
        _check_df(name, c77, opt)
        n77 = int(np.asarray(c77["int_node"]).size)
        if n77 > md.nint:
            lst.warning(f"{name} has factors for {n77 - md.nint} nodes that are not interaction nodes of FILE4 "
                        "(not used)")
        fs.files[name] = c77
        fs.rows[name] = {int(n): k for k, n in enumerate(np.asarray(c77["fnum"]))}
    for c in cases:
        meta = fs.files[c.f77].meta
        if bool(c.sim) != bool(int(meta.get("stochastic", 0) or 0)):
            lst.warning(f"{c.f77}: written by a {'stochastic' if int(meta.get('stochastic', 0) or 0) else 'deterministic'} "
                        f"HOUSE run, used as {'a stochastic sample' if c.sim else 'deterministic factors'}")
        if str(meta.get("method", "")) == "SINGLE" and int(meta.get("mode", 0) or 0) > 1:
            c.lowfreq_skip = f"single incoherent mode {int(meta['mode'])}: its rigid-body (f -> 0) limit is 0"
        elif int(meta.get("me", 0) or 0):
            c.lowfreq_skip = "multiple excitation: the low-frequency value is the SAR of the zone"
    return fs


def survey_factors(lst, fnums: List[int], fs: FactorSet) -> None:
    """Every SSI frequency must be in every FILE77 used (frequency survey, D-ANL-08)."""
    for name, rows in fs.rows.items():
        miss = [n for n in fnums if n not in rows]
        if miss:
            raise ModuleError(f"frequency survey failed: Frequency {', '.join(str(v) for v in miss[:10])}"
                              f"{' ...' if len(miss) > 10 else ''} not in {name} -- re-run HOUSE with the frequency "
                              "set of ANALYS")
    names = list(fs.rows)
    lst.write(f" all {len(fnums)} frequencies found in {names[0]}" + (f" .. {names[-1]} ({len(names)} files)"
                                                                         if len(names) > 1 else ""))


def _check_df(name: str, c: Container, opt: Options) -> None:
    dfc = float(c.meta["df"])
    if opt.df <= 0:
        opt.df = dfc
    elif abs(dfc - opt.df) > DF_RTOL * opt.df:
        raise ModuleError(f"frequency step of {name} ({dfc:.9g} Hz) differs from the ANALYS step {opt.df:.9g} Hz "
                          "(D-ANL-08): re-run the producing module with the same NFFT / time step / frequency step")


def _check_depths(name: str, c: Container, md: ModelData) -> None:
    """The site of FILE1/FILE3 must be the site of HOUSE (same user interfaces)."""
    if "x_iface_depth" not in md.file4 or "depth_user" not in c:
        return
    a = np.asarray(md.file4["x_iface_depth"], float)
    b = np.asarray(c["depth_user"], float)
    tol = 1e-6 * max(1.0, float(np.abs(a).max()) if a.size else 1.0)
    if a.size != b.size or (a.size and float(np.abs(a - b).max()) > tol):
        raise ModuleError(f"{name} was computed for another site layering (interface depths {b.tolist()[:8]} vs "
                          f"HOUSE {a.tolist()[:8]}) -- re-run SITE/POINT with the HOUSE site")


def read_load_cases(ctx: ModuleContext, opt: Options, md: ModelData, cases: List[LoadCase], lst) -> None:
    producer = "SITE (Mode 2)" if opt.seismic else "FORCE"
    kind = "FILE1" if opt.seismic else "FILE9"
    seen: Dict[str, Container] = {}               # stochastic cases share FILE1X/Y/Z: read each file once
    warned: set = set()
    for c in cases:
        if c.infile in seen:
            c.data = seen[c.infile]
        else:
            if not ctx.path(c.infile).exists():
                hint = ""
                if c.infile != kind:
                    hint = f" (run {producer.split()[0]} and copy {kind} to {c.infile}, e.g. FCOPY,{kind},{c.infile})"
                raise ModuleError(f"{c.infile} missing -- run {producer} first{hint}")
            c.data = seen[c.infile] = _read(ctx, c.infile, kind, producer)
            _check_df(c.infile, c.data, opt)
            if opt.seismic:
                _check_depths(c.infile, c.data, md)
        c.rows = {int(n): k for k, n in enumerate(np.asarray(c.data["fnum"]))}
        if opt.seismic:
            cm = int(c.data.meta["cm"])
            if (c.label in ("X", "Y", "Z") or c.sim) and cm != c.cm and c.infile not in warned:
                warned.add(c.infile)
                lst.warning(f"{c.infile}: control direction {('x', 'y', 'z')[cm]}' but the "
                            f"{'XYZ'[c.cm]} case needs {('x', 'y', 'z')[c.cm]}' (SITE runs: SV x' / SH y' / P z', "
                            "spec 05c A.5.4)")
            c.cm = cm


# ======================================================================================
# Frequency survey and restart records
# ======================================================================================
@dataclass
class RestartIndex:
    kind: str                          # 'COOX' or 'COOTK'
    order: Dict[int, int]              # frequency number -> order number qqq
    meta: Dict[str, object]


def _same_frequency(saved: float, f: float) -> bool:
    """Frequency values agree within the relative tolerance of the frequency-step check (D-ANL-08)."""
    return abs(float(saved) - float(f)) <= DF_RTOL * abs(float(f))


def read_index(ctx: ModuleContext, kind: str, md: ModelData, opt: Options) -> RestartIndex:
    """COOXI / COOTKI with the FILE90 hash validation and the frequency-step check (requirements
    2.5 and 4.6 item 8, D-ANL-05, D-ANL-08).

    The records are keyed by frequency *value*: a record saved as frequency number n at the step
    df_saved holds X_ff (or the factorised system) of f = n df_saved, which is a different frequency
    when this run uses another step (another NFFT, time step or deck frequency step).  The FILE90
    hashes do not contain the frequencies, so the step stored in the index (meta ``df``) and the
    frequency of every entry (array ``freq``) are compared with this run before any record is used."""
    name = f"{kind}I"
    p = ctx.path(name)
    if not p.exists():
        raise ModuleError(f"{name} missing -- run ANALYS (initiation) with Save Restart Files (<save> = 1)")
    try:
        c = read_container(p, kind=name)
    except (ValueError, OSError) as exc:
        raise ModuleError(f"{name}: {exc}") from None
    keys = ("int_hash", "layer_hash") + (("file4_hash",) if kind == "COOTK" else ())
    bad = [k for k in keys if str(c.meta.get(k)) != md.hashes.get(k)]
    if bad:
        what = {"int_hash": "interaction nodes", "layer_hash": "soil layers",
                "file4_hash": "structure model (FILE4/COOSK/COOSM)"}
        raise ModuleError(f"{name}: restart files were saved for different " + " and ".join(what[k] for k in bad)
                          + " (FILE90 hash mismatch) -- run an initiation (<mode> = 0) with <save> = 1")
    fnum = np.asarray(c["fnum"], np.int64).reshape(-1)
    df_saved = c.meta.get("df")
    if df_saved is None:
        raise ModuleError(f"{name} does not record the frequency step of its restart files -- run an initiation "
                          "(<mode> = 0) with <save> = 1")
    if not abs(float(df_saved) - opt.df) <= DF_RTOL * opt.df:
        raise ModuleError(f"{name}: the restart files were saved for the frequency step df = {float(df_saved):.9g} Hz, "
                          f"this run uses df = {opt.df:.9g} Hz (D-ANL-05/08): frequency number n is another "
                          "frequency now -- use the NFFT / time step / frequency step of the initiation, or run an "
                          "initiation (<mode> = 0) with <save> = 1")
    if "freq" in c:
        fsaved = np.asarray(c["freq"], float).reshape(-1)
        for n, fv in zip(fnum, fsaved):
            if not _same_frequency(fv, n * opt.df):
                raise ModuleError(f"{name}: frequency number {int(n)} was saved at f = {fv:.9g} Hz, this run has "
                                  f"f = {n * opt.df:.9g} Hz (D-ANL-05) -- run an initiation with <save> = 1")
    order = {int(n): int(o) for n, o in zip(fnum, np.asarray(c["order"]))}
    return RestartIndex(kind, order, dict(c.meta))


def frequency_list(opt: Options, cases: List[LoadCase], lst=None) -> List[int]:
    """Frequency numbers to solve: the deck frequency set (<fopt> = 0) or every frequency of the first
    FILE1/FILE9 (<fopt> = 1).  More than :data:`MAX_FREQ` gives the EDU-02 warning (manual limit, as
    in SITE); the run continues: D-GEN-07 makes the manual size limits warnings (hard errors only where
    the physics needs them).  Order numbers from 1000 on give 4-digit restart names (COOX1000), which
    :data:`_RESTART_RE` covers."""
    if opt.fopt == 1:
        fn = sorted(int(n) for n in np.asarray(cases[0].data["fnum"]))
    else:
        fn = list(opt.freqs)
    if not fn:
        raise ModuleError("Error 120: the frequency set is empty")
    if len(fn) > MAX_FREQ and lst is not None:
        lst.warning(f"EDU-02: {len(fn)} SSI frequencies exceed the manual limit of {MAX_FREQ} per ANALYS run "
                    "(requirements 4.0.1); split the set and merge the FILE8s with COMBIN to stay within it")
    return fn


def survey(lst, opt: Options, fnums: List[int], cases: List[LoadCase], file3: Optional[Container],
           coox: Optional[RestartIndex], cootk: Optional[RestartIndex]) -> Dict[int, int]:
    """Frequency survey table (requirements 2.6, D-ANL-08); returns {fnum: FILE3 row}."""
    rows3 = {int(n): k for k, n in enumerate(np.asarray(file3["fnum"]))} if file3 is not None else {}
    lst.section("Frequency survey")
    cases = list({c.infile: c for c in cases}.values())     # one column per input file (stochastic cases share them)
    cols = [c.infile for c in cases]
    if file3 is not None:
        cols.append("FILE3")
    if coox is not None:
        cols.append("COOXqqq")
    if cootk is not None:
        cols.append("COOTKqqq")
    lst.write(f"{'q':>5s}{'number':>9s}{'f (Hz)':>14s}" + "".join(f"{c:>10s}" for c in cols) + "  status")
    missing: Dict[str, List[int]] = {}
    for q, n in enumerate(fnums, start=1):
        cells, ok = [], True
        for c in cases:
            r = c.rows.get(n)
            cells.append("-" if r is None else str(r + 1))
            if r is None:
                ok = False
                missing.setdefault(c.infile, []).append(n)
        if file3 is not None:
            r = rows3.get(n)
            cells.append("-" if r is None else str(r + 1))
            if r is None:
                ok = False
                missing.setdefault("FILE3", []).append(n)
        for idx, nm in ((coox, "COOX"), (cootk, "COOTK")):
            if idx is None:
                continue
            o = idx.order.get(n)
            cells.append("-" if o is None else f"{o:03d}")
            if o is None:
                ok = False
                missing.setdefault(nm + "I", []).append(n)
        lst.write(f"{q:5d}{n:9d}{n * opt.df:14.6g}" + "".join(f"{v:>10s}" for v in cells)
                  + ("  ok" if ok else "  MISSING"))
    if missing:
        msgs = []
        for nm, lst_n in missing.items():
            sample = ", ".join(str(v) for v in lst_n[:10]) + (" ..." if len(lst_n) > 10 else "")
            prod = {"FILE3": "POINT", "COOXI": "ANALYS with <save> = 1", "COOTKI": "ANALYS with <save> = 1"}.get(
                nm, "SITE" if nm.startswith("FILE1") else "FORCE")
            msgs.append(f"Frequency {sample} not in {nm} -- re-run {prod}")
        raise ModuleError("frequency survey failed: " + "; ".join(msgs))
    lst.write(f" all {len(fnums)} frequencies found")
    return rows3


def _restart_record(ctx: ModuleContext, idx: RestartIndex, fnum: int, md: ModelData, opt: Options) -> Container:
    """Restart record COOXqqq / COOTKqqq of frequency number ``fnum``, validated by its frequency
    number, its frequency value ``f = fnum df`` and the FILE90 hashes (requirements 4.6 item 8,
    D-ANL-05)."""
    order = idx.order[fnum]
    name = restart_names(order)[0 if idx.kind == "COOX" else 1]
    p = ctx.path(name)
    if not p.exists():
        raise ModuleError(f"{name} (frequency number {fnum}) missing -- the restart files are incomplete")
    try:
        c = read_container(p, kind=idx.kind)
    except (ValueError, OSError) as exc:
        raise ModuleError(f"{name}: {exc}") from None
    if int(c.meta.get("fnum", -1)) != int(fnum):
        raise ModuleError(f"{name} holds frequency number {c.meta.get('fnum')}, the index {idx.kind}I says {fnum} "
                          "-- the restart files were overwritten; run an initiation with <save> = 1")
    f = fnum * opt.df
    fsaved = c.meta.get("freq")
    if fsaved is None or not _same_frequency(fsaved, f):
        raise ModuleError(f"{name} holds frequency f = {fsaved} Hz, frequency number {fnum} of this run is "
                          f"{f:.9g} Hz (D-ANL-05) -- the restart files were saved for another frequency step; "
                          "run an initiation with <save> = 1")
    keys = ("int_hash", "layer_hash") + (("file4_hash",) if idx.kind == "COOTK" else ())
    if any(str(c.meta.get(k)) != md.hashes.get(k) for k in keys):
        raise ModuleError(f"{name}: FILE90 hash mismatch (the model changed since the restart file was saved)")
    return c


# ======================================================================================
# Soil impedance
# ======================================================================================
def soil_impedance(file3: Container, row: int, xy: np.ndarray, iface: np.ndarray) -> Tuple[np.ndarray, float]:
    """``X_ff = F_ff^-1`` at one frequency with a single dense ``3 nInt x 3 nInt`` array (steps 2-3;
    D-ANL-02, D-ANL-04, D-ANL-10): ``F_ff`` is built by
    :func:`sassi.core.flexibility.flexibility_matrix`, symmetrised in place and inverted in its own
    buffer.  Same values as ``impedance_matrix(flexibility_matrix(...))`` up to round-off.  Returns
    ``(X_ff, rcond(F_ff))``."""
    n = int(np.asarray(xy).reshape(-1, 2).shape[0])
    F = flexibility_matrix(file3, row, xy, iface, symmetrize=False, block=max(1, FLEX_BLOCK_ELEMS // max(n, 1)))
    ss.symmetrize_inplace(F)                       # (F + F^T) / 2, D-ANL-02
    return ss.impedance_matrix(F, overwrite=True)


def reduced_soil_impedance(file3: Container, row: int, md: ModelData) -> Tuple[np.ndarray, float]:
    """``X_ff`` of a 2D model and / or a SYMM half / quarter model (requirements 4.6 items 2-3, D-ANL-12):

    * the flexibility of the interaction nodes -- with SYMM planes the image sum
      ``F_red = sum_S s_S F(nodes, images_S(nodes)) P_S`` (:func:`sassi.core.symmetry.reduced_flexibility`,
      blocks by :func:`sassi.core.flexibility.flexibility_block`); POINT2 files give the plane-strain blocks;
    * symmetrised, restricted to the active translations (2D: UX, UZ; SYMM: without the translations the
      symmetry constrains on the planes, whose rows and columns of ``F_red`` vanish) and inverted
      (:func:`sassi.core.symmetry.reduced_impedance`); the removed DOFs get zero rows and columns.

    Returns ``(X_ff (3 nInt x 3 nInt), rcond of the inverted block)``."""
    xy, iface = md.int_xyz[:, :2], md.int_iface
    blk = max(1, FLEX_BLOCK_ELEMS // max(md.nint, 1))
    if md.symm:
        def block(xyz_load: np.ndarray) -> np.ndarray:
            return flexibility_block(file3, row, xy, iface, np.asarray(xyz_load)[:, :2], iface, block=blk)
        F = SYM.reduced_flexibility(block, md.int_xyz, md.symm)
    else:
        F = flexibility_matrix(file3, row, xy, iface, symmetrize=False, block=blk)
    return SYM.reduced_impedance(F, md.keep_mask())


#: relative tolerance of the symmetry check of the seismic free field (D-ANL-12)
SYMM_FIELD_RTOL = 1e-6


def model_rules(opt: Options, md: ModelData) -> None:
    """Option combinations of 2D and SYMM models (requirements 4.4 item 6, 4.6 item 9, G-17, D-ANL-12):
    incoherent motion, wave passage and multiple excitation need a full 3D model, the global impedance
    a 3D model without SYMM planes."""
    if md.symm:
        planes = "; ".join(p.describe(md.dim) for p in md.symm)
        if opt.factored:
            raise ModuleError("incoherent motion, wave passage and multiple excitation (HOUSE <coh>/<wpass>/<me>, "
                              f"FILE77) are not allowed with SYMM planes ({planes}): they need the full model "
                              "(manual 2.7 and 6.5.4, G-17)")
        if opt.impe:
            raise ModuleError(f"global impedance (<impe> = {opt.impe}) is not allowed with SYMM planes ({planes}): "
                              "the impedance of a half / quarter interaction set is not the foundation impedance "
                              "(manual 2.7, G-17)")
    if md.dim == 1:
        if opt.factored:
            raise ModuleError("incoherent motion, wave passage and multiple excitation (HOUSE <coh>/<wpass>/<me>) need "
                              "a 3D model (HOUSE <dim> = 2; manual 6.5.4, G-17)")
        if opt.impe:
            raise ModuleError(f"global impedance (<impe> = {opt.impe}) is available for 3D models only (requirements "
                              "4.6 item 9); 2D model (HOUSE <dim> = 1)")


def check_symmetric_input(md: ModelData, case: LoadCase, fnum: int, opt: Options, xyz: np.ndarray,
                          U: np.ndarray) -> None:
    """The free field of a SYMM model must have the symmetry of every plane: ``U'(image) = s P U'`` with
    s = +1 (type 0) or -1 (type 1) (D-ANL-12).  A vertically incident horizontal input is antisymmetric
    about the plane normal to its direction and symmetric about the other one; vertical (P) input is
    symmetric; inclined or surface waves propagating across a plane have neither symmetry."""
    imgs = [free_field_at_nodes(case.data, case.rows[fnum], p.mirror(xyz), md.int_iface, opt.ang, opt.xc, opt.yc)
            for p in md.symm]
    comps = (0, 2) if md.dim == 1 else (0, 1, 2)
    errs = SYM.field_symmetry_error(md.symm, U, imgs, comps)
    for p, (e_sym, e_anti) in zip(md.symm, errs):
        e = e_sym if p.type == SYM.SYMMETRIC else e_anti
        if e <= SYMM_FIELD_RTOL:
            continue
        other = SYM.ANTISYMMETRIC if p.type == SYM.SYMMETRIC else SYM.SYMMETRIC
        fits = min(e_sym, e_anti) <= SYMM_FIELD_RTOL
        hint = (f"it is {SYM.TYPE_NAMES[other]}: use SYMM type {other} for this plane, or another model for "
                "this input" if fits else
                "it is neither symmetric nor antisymmetric (inclined or surface waves crossing the plane, or a "
                "rotated input): analyse the full model")
        raise ModuleError(f"{case.infile} (frequency number {fnum}): the seismic free field does not have the "
                          f"{SYM.TYPE_NAMES[p.type]} of SYMM {p.describe(md.dim)} (relative deviation {e:.2e}); "
                          f"{hint} (an X input is antisymmetric about planes normal to X and symmetric about planes "
                          "normal to Y; a Z input is symmetric about both; D-ANL-12)")


# ======================================================================================
# Resource check (D-ANL-10, G-22, EDU-20)
# ======================================================================================
def physical_memory() -> Optional[int]:
    try:
        return int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))
    except (ValueError, OSError, AttributeError):
        return None


@dataclass
class MemoryEstimate:
    """Peak-memory estimate of one ANALYS run in bytes (D-ANL-10, G-22); see :func:`memory_estimate`."""

    inputs: int
    dense: int
    dense_copies: float
    work: int
    sparse: int
    results: int
    file11: int
    fallback: int

    @property
    def total(self) -> int:
        return self.inputs + self.dense + self.work + self.sparse + self.results + self.file11

    @property
    def fallback_total(self) -> int:
        """Peak of a frequency on the full sparse LU fallback path (it replaces the Schur work)."""
        return self.inputs + self.fallback + self.sparse + self.results + self.file11


def memory_estimate(md: ModelData, nF: int, ncases: int, opt: Options, input_bytes: int = 0) -> MemoryEstimate:
    """Peak memory of the frequency loop (D-ANL-10; requirements 4.6 item 10).

    With ``D = (3 nInt)^2 x 16 B`` (one dense complex matrix on the interaction DOFs) one frequency
    holds at most (measured with tracemalloc, tests/unit/test_analys_memory.py):

    * soil impedance: ``F_ff`` built and inverted in place (1 D) + the flexibility work arrays
      (about 0.6 D of coordinate tables + the assembly blocks of :data:`FLEX_BLOCK_ELEMS`);
    * Schur path: ``X_ff`` + ``S + X_ff`` / its LU (2 D) + the bounded Schur work
      (:func:`sassi.core.ssi_solver.schur_work_bytes`);
    * restart files (<save> = 1): the record (1 D) + its in-memory write buffer, which may grow to
      :data:`WRITE_COPIES` D while it is written (sassi.io.container writes through a byte buffer).

    The estimate is ``inputs + DENSE_COPIES D (DENSE_COPIES_SAVE D with <save> = 1) + max(work) +
    sparse factors + transfer functions + FILE11``; the full sparse LU fallback (used only at
    frequencies where ``C_nn`` is ill-conditioned) needs about :data:`FALLBACK_COPIES` D plus the
    sparse factors instead of the dense and work terms and is reported separately
    (:attr:`MemoryEstimate.fallback_total`).  Not counted: the LAPACK getri work array (64 x 3 nInt
    complex numbers) and other arrays linear in the number of equations."""
    n = md.nint
    n3 = 3 * n
    D = n3 * n3 * 16
    if opt.solver_mode == 3:                      # X_ff (COOX) and the LU of S + X_ff (COOTK) are read
        copies = 2.0
    else:
        copies = DENSE_COPIES_SAVE if opt.save else DENSE_COPIES
        if md.reduced and opt.mode == 0:          # 2D / SYMM: the image block or the kept block is a copy
            copies += 1.0
    rows = min(n, max(1, FLEX_BLOCK_ELEMS // max(n, 1)))
    flex_work = FLEX_WORK_ARRAYS * 16 * rows * n if opt.mode == 0 else 0
    tile = int(np.sqrt(ss.DENSE_BLOCK_ELEMS))
    tile_work = 2 * 16 * min(n3, tile) ** 2
    work = max(flex_work, ss.schur_work_bytes(md.part.nn, md.part.nf), tile_work) if opt.solver_mode != 3 else 0
    nnz = sum(int(A.nnz) for A in (md.Ks, md.Ke) if A is not None) or 30 * md.neq
    sparse = SPARSE_FILL * nnz * 16
    per_case = nF * md.neq * 16
    results = per_case * ncases + int(WRITE_COPIES * per_case)
    xall = nF * D if (opt.impe and nF * D <= file11_x_limit()) else 0
    file11 = int((1 + WRITE_COPIES) * xall)
    fallback = int(FALLBACK_COPIES * D)
    return MemoryEstimate(inputs=int(input_bytes), dense=int(copies * D), dense_copies=copies, work=int(work),
                          sparse=sparse, results=results, file11=file11, fallback=fallback)


def resource_check(lst, md: ModelData, nF: int, ncases: int, opt: Options, input_bytes: int = 0) -> MemoryEstimate:
    """Listing of :func:`memory_estimate` and the EDU-20 guard (G-22): refuse the run when the
    estimate exceeds :data:`MEMLIMIT_FRACTION` of the physical memory; warn when only the full-LU
    fallback would exceed it."""
    est = memory_estimate(md, nF, ncases, opt, input_bytes)
    mb = 2.0 ** 20
    ram = physical_memory()
    lst.section("Resource check (D-ANL-10)")
    lst.write(f" input files held in memory                        : {est.inputs / mb:12.1f} MB")
    lst.write(f" dense matrices {est.dense_copies:g} x (3 x {md.nint})^2 x 16 B            : {est.dense / mb:12.1f} MB")
    lst.write(f" bounded work arrays (flexibility / Schur blocks)  : {est.work / mb:12.1f} MB")
    lst.write(f" sparse factors (estimate, fill {SPARSE_FILL})            : {est.sparse / mb:12.1f} MB")
    lst.write(f" transfer functions {nF} x {md.neq} x {ncases} case(s) + buffer : {est.results / mb:12.1f} MB")
    if est.file11:
        lst.write(f" FILE11 impedance matrices + write buffer          : {est.file11 / mb:12.1f} MB")
    lst.write(f" estimated peak memory                             : {est.total / mb:12.1f} MB")
    if md.part.nf and md.part.nn:
        lst.write(f" a frequency on the full sparse LU fallback path needs about {est.fallback_total / mb:.1f} MB "
                  f"({FALLBACK_COPIES:g} x (3 x {md.nint})^2 x 16 B instead of the dense matrices and work arrays)")
    if ram:
        lst.write(f" physical memory                                   : {ram / mb:12.1f} MB")
        if est.total > MEMLIMIT_FRACTION * ram:
            raise ModuleError(f"EDU-20 (G-22): the estimated memory {est.total / 2 ** 30:.2f} GB exceeds "
                              f"{100 * MEMLIMIT_FRACTION:.0f} % of the physical memory ({ram / 2 ** 30:.2f} GB): reduce "
                              "the interaction nodes (FFV/FI), the frequencies per run or the simultaneous cases")
        fb = est.fallback_total
        if md.part.nf and md.part.nn and fb > MEMLIMIT_FRACTION * ram:
            lst.warning(f"EDU-20 (G-22): a frequency that needs the full sparse LU fallback (ill-conditioned C_nn) "
                        f"would need about {fb / 2 ** 30:.2f} GB, more than {100 * MEMLIMIT_FRACTION:.0f} % of the "
                        "physical memory: the run may fail at such a frequency")
    else:
        lst.write(" physical memory unknown: check skipped")
    return est


# ======================================================================================
# Loads
# ======================================================================================
def vibration_loads(md: ModelData, cases: List[LoadCase], fnum: int, lst=None, warned=None) -> np.ndarray:
    """FILE9 load vectors of every case on the FILE4 equations, (neq, ncases)."""
    eqmap = md.eqmap
    b = np.zeros((md.neq, len(cases)), complex)
    for j, c in enumerate(cases):
        r = c.rows[fnum]
        P = np.asarray(c.data["P"])[r]
        for k, (n, dof) in enumerate(zip(np.asarray(c.data["load_node"]), np.asarray(c.data["load_dof"]))):
            e = eqmap.get((int(n), int(dof)))
            if e is None:
                if int(n) not in np.asarray(md.file4["node_id"]):
                    raise ModuleError(f"{c.infile}: load on node {int(n)}, which is not in the HOUSE model (Error 62)")
                if warned is not None and (c.infile, int(n), int(dof)) not in warned:
                    warned.add((c.infile, int(n), int(dof)))
                    if lst is not None:
                        lst.warning(f"{c.infile}: load on node {int(n)} {DOF_LABELS[int(dof)]} ignored -- the DOF is "
                                    "fixed or not active in the HOUSE model (W10/W11)")
                continue
            b[e, j] += P[k]
    return b


def _singular_hint(md: ModelData, C) -> str:
    """Diagnosis of a singular system: equations without any stiffness or mass (e.g. the rotations a
    SPRING or SHELL adds to a node that nothing else restrains, EDU-06)."""
    if C is None:
        return ""
    rows = np.flatnonzero(np.asarray(abs(C).sum(axis=1)).ravel() == 0.0)
    plane = (" (2D model: the out-of-plane DOFs UY, ROTX, ROTZ of BEAMS/SPRING/GENERAL nodes are structural only "
             "-- fix them with D unless something restrains them)") if md.dim == 1 else ""
    if rows.size:
        what = ", ".join(f"node {int(md.eq_node[i])} {DOF_LABELS[int(md.eq_dof[i])]}" for i in rows[:10])
        return (f" -- {rows.size} equation(s) without stiffness or mass: {what}{' ...' if rows.size > 10 else ''}; "
                "fix them (D, FIXROT/FIXSHLROT; see EDU-06 in the HOUSE listing)" + plane)
    return (" -- check the model for mechanisms (unrestrained rotations, EDU-06 in the HOUSE listing; hinges; "
            "undamped resonances of a structure without interaction nodes)" + plane)


# ======================================================================================
# Output files
# ======================================================================================
def _fmt(v: float) -> str:
    return f"{v:.16e}"


def write_impedance_files(ctx: ModuleContext, opt: Options, md: ModelData, fnums: List[int],
                          KG: np.ndarray, T: np.ndarray, Xall: Optional[np.ndarray]) -> List[str]:
    """FOUNSTIF, FOUNDASH, FOUNDAMP, FOUNIMPD (text, one row per frequency) and FILE11 (D-ANL-07)."""
    f = np.asarray(fnums, float) * opt.df
    w = 2.0 * np.pi * f
    re, im = KG.real, KG.imag
    with np.errstate(divide="ignore", invalid="ignore"):
        damp = np.where(re != 0.0, im / (2.0 * np.abs(re)), np.where(im == 0.0, 0.0, np.inf))
    quantities = {"FOUNSTIF": ("dynamic stiffness Re K", re),
                  "FOUNDASH": ("viscous damping coefficient Im K / w", im / w[:, None, None]),
                  "FOUNDAMP": ("effective damping ratio Im K / (2 |Re K|)", damp),
                  "FOUNIMPD": ("absolute value |K|", np.abs(KG))}
    if opt.impe == 1:
        labels = [f"K{d}{d}" for d in RIGID_DOFS]
    else:
        labels = [f"K{a}-{b}" for a in RIGID_DOFS for b in RIGID_DOFS]
    names = []
    for name, (what, V) in quantities.items():
        lines = [f"# SASSI-EDU ANALYS global (unconstrained) soil impedance K_G = T^T X_ff T (D-ANL-07)",
                 f"# {name}: {what}; reference point ({opt.xc:g}, {opt.yc:g}, {opt.zc:g}); "
                 f"{'diagonal terms' if opt.impe == 1 else '6x6 matrix row-major'}",
                 "# columns: f(Hz) " + " ".join(labels)]
        for q in range(len(fnums)):
            vals = np.diagonal(V[q]) if opt.impe == 1 else V[q].reshape(-1)
            lines.append(f"{f[q]:.10f} " + " ".join(_fmt(float(v)) for v in vals))
        ctx.path(name).write_text("\n".join(lines) + "\n", encoding="utf-8")
        names.append(name)
    arrays = {"fnum": np.asarray(fnums, np.int64), "freq": f, "KG": KG, "T": T, "int_node": md.int_node,
              "int_xyz": md.int_xyz}
    if Xall is not None:
        arrays["X"] = Xall
    meta = {"df": opt.df, "impe": opt.impe, "ref": [opt.xc, opt.yc, opt.zc], "dof_order": list(RIGID_DOFS),
            "x_stored": Xall is not None, "model_hash": md.model_hash,
            "note": "K_G = T^T X_ff T over the translational interaction DOFs (unconstrained impedance)"
                    + ("" if Xall is not None else f"; X_ff not stored (more than {file11_x_limit() / 2 ** 20:g} MB, "
                                                     "save the restart files COOXqqq instead)")}
    write_container(ctx.path("FILE11"), "FILE11", arrays, meta, module=NAME)
    names.append("FILE11")
    return names


def _delete_restart_files(ctx: ModuleContext) -> List[str]:
    gone = []
    for p in sorted(ctx.workdir.iterdir()):
        if p.is_file() and _RESTART_RE.match(p.name):
            p.unlink()
            gone.append(p.name)
    return gone


# ======================================================================================
# Listing helpers
# ======================================================================================
def _list_input(lst, opt: Options, md: ModelData, cases: List[LoadCase], ctx: ModuleContext) -> None:
    lst.section("ANALYS input")
    lst.write(f" Title                          : {opt.title}")
    lst.write(f" Operation mode                 : {'data check only' if opt.opmode else 'solution'}")
    lst.write(f" Type of analysis               : {'seismic' if opt.seismic else 'foundation vibration'}")
    lst.write(f" Mode of analysis               : {opt.mode} {MODE_NAMES[opt.mode]} (solver Mode {opt.solver_mode})")
    lst.write(f" Save restart files             : {'yes' if opt.save else 'no'}")
    lst.write(f" Print                          : {'amplitudes only' if opt.prnt else 'real and imaginary parts'}")
    lst.write(f" Frequencies                    : "
              f"{'all of ' + cases[0].infile if opt.fopt else f'frequency set {opt.freqset} ({len(opt.freqs)} numbers)'}")
    lst.write(f" Frequency step                 : {opt.df:.9g} Hz")
    if opt.seismic:
        lst.write(f" Coordinate transformation angle: {opt.ang:g} deg (x' -> x)")
    lst.write(f" Control / reference point      : ({opt.xc:g}, {opt.yc:g}, {opt.zc:g})")
    lst.write(f" Global impedance               : {IMPE_NAMES[opt.impe]}")
    items = [f"{c.infile}{(' x ' + c.f77) if c.f77 else ''} -> {c.outfile}" for c in cases]
    if len(items) > 9:
        items = items[:6] + [f"... ({len(cases)} cases)"] + items[-2:]
    lst.write(f" Simultaneous cases             : {opt.simul}: " + ", ".join(items))
    lst.section("Model (FILE4 = " + f"{ctx.model}.N4)")
    lst.write(f" dimension                      : {S2.dimension_name(md.dim)}"
              + ("; interaction DOFs UX, UZ, per unit length (D-W2-07)" if md.dim == 1 else ""))
    lst.write(f" equations                      : {md.neq}")
    if md.dim == 1:
        lst.write(f" interaction nodes              : {md.nint} ({md.part.nf} interaction DOFs UX, UZ, "
                  f"{2 * md.nint - md.part.nf} fixed interaction translations)")
    else:
        lst.write(f" interaction nodes              : {md.nint} ({md.part.nf} interaction DOFs, "
                  f"{3 * md.nint - md.part.nf} fixed interaction translations)")
    if md.symm:
        cons = md.int_symm[:, [0, 2]] if (md.int_symm is not None and md.dim == 1) else md.int_symm
        ncons = int(cons.sum()) if cons is not None else 0
        what = "half" if len(md.symm) == 1 else "quarter"
        lst.write(f" symmetry (SYMM, D-ANL-12)      : {what} model; "
                  + "; ".join(p.describe(md.dim) for p in md.symm))
        lst.write(f"                                  soil flexibility by {2 ** len(md.symm)} image terms "
                  f"F_red = sum_S s_S F(i, j_S) P_S; {ncons} interaction translations constrained on the "
                  "plane(s) (removed before X = F_red^-1)")
        if not opt.seismic:
            lst.write("                                  vibration: loads on plane nodes are those of the "
                      "reduced model (1/2 or 1/4 of the full-model load)")
    lst.write(f" non-interaction equations      : {md.part.nn}")
    lst.write(f" method (HOUSE <imp>)           : {METHOD_NAMES.get(md.method, str(md.method))}; the interaction "
              "set defines the method (R1 4.5)")
    if md.nint:
        ifs = np.unique(md.int_iface)
        lst.write(f" interaction interfaces         : {', '.join(str(int(i)) for i in ifs)}")


def _list_summary(lst, recs: List[FreqRecord], mode0: bool) -> None:
    lst.section("Per-frequency summary (solver path: schur = sparse LU of C_nn + dense LU of S + X_ff; "
                "dense = all equations interact; full = sparse LU of the whole system)")
    lst.write(f"{'q':>5s}{'number':>8s}{'f (Hz)':>12s}{'path':>8s}{'n_f':>7s}{'n_n':>8s}"
              + (f"{'rcond F_ff':>12s}" if mode0 else "") + f"{'rcond C_nn':>12s}{'rcond S+X':>12s}{'time (s)':>10s}")
    for r in recs:
        lst.write(f"{r.q:5d}{r.fnum:8d}{r.f:12.6g}{r.path:>8s}{r.nf:7d}{r.nn:8d}"
                  + (f"{r.rcond_F:12.3e}" if mode0 else "")
                  + f"{r.rcond_nn:12.3e}{r.rcond_ff:12.3e}{r.seconds:10.3f}")
    notes = [(r.q, r.note) for r in recs if r.note]
    for q, n in notes[:50]:
        lst.write(f"   q = {q}: {n}")


def _print_tfs(lst, opt: Options, md: ModelData, fnums: List[int], case: LoadCase, H: np.ndarray) -> None:
    nodes = np.unique(md.eq_node)
    shown = nodes[:MAX_PRINT_NODES]
    col = {(int(n), int(k)): i for i, (n, k) in enumerate(zip(md.eq_node, md.eq_dof))}
    what = "amplitudes" if opt.prnt else "real and imaginary parts"
    kind = "U/U_cp" if opt.seismic else "displacement per unit load factor"
    lst.section(f"Transfer functions {case.outfile} ({kind}, {what})")
    if nodes.size > shown.size:
        lst.write(f" the first {shown.size} of {nodes.size} nodes are printed (MAX_PRINT_NODES); use MOTION for others")
    head = f"{'node':>8s}{'':>4s}" + "".join(f"{DOF_LABELS[k]:>14s}" for k in range(1, 7))
    for q, n in enumerate(fnums):
        lst.write(f" frequency {q + 1}: number {n}, f = {n * opt.df:.6g} Hz")
        lst.write(head)
        for nd in shown:
            vals = [H[q, col[(int(nd), k)]] if (int(nd), k) in col else 0.0 for k in range(1, 7)]
            if opt.prnt:
                lst.write(f"{int(nd):8d}{'':>4s}" + "".join(f"{abs(v):14.6e}" for v in vals))
            else:
                lst.write(f"{int(nd):8d}{'Re':>4s}" + "".join(f"{np.real(v):14.6e}" for v in vals))
                lst.write(f"{'':8s}{'Im':>4s}" + "".join(f"{np.imag(v):14.6e}" for v in vals))


def _print_impedance(lst, opt: Options, fnums: List[int], KG: np.ndarray, embedded: bool) -> None:
    lst.section(f"Global soil impedance K_G = T^T X_ff T about ({opt.xc:g}, {opt.yc:g}, {opt.zc:g}) (D-ANL-07)")
    lst.write(" these are 'unconstrained' impedances: they equal the rigid-foundation impedances only for surface "
              "foundations (spec 05c A.5.10)")
    if embedded:
        lst.warning("global impedance of an embedded interaction set: not the rigid-foundation impedance "
                    "(use FI-FSIN for embedded foundations, spec 05c A.5.10)")
    if opt.impe == 1:
        lst.write(f"{'f (Hz)':>12s}" + "".join(f"{'Re K' + d:>14s}{'Im K' + d:>14s}" for d in RIGID_DOFS))
        for q, n in enumerate(fnums):
            dg = np.diagonal(KG[q])
            lst.write(f"{n * opt.df:12.6g}" + "".join(f"{v.real:14.6e}{v.imag:14.6e}" for v in dg))
        return
    for q, n in enumerate(fnums):
        lst.write(f" f = {n * opt.df:.6g} Hz (rows/columns {' '.join(RIGID_DOFS)})")
        for part, A in (("Re", KG[q].real), ("Im", KG[q].imag)):
            for i in range(6):
                lst.write(f"{part:>6s} {RIGID_DOFS[i]:>3s}" + "".join(f"{v:14.6e}" for v in A[i]))


def _lowfreq_check(lst, opt: Options, md: ModelData, cases: List[LoadCase], H: Dict[str, np.ndarray],
                   fnums: List[int]) -> None:
    """G-19 / EDU-18: at the first SSI frequency every translational ATF in the control direction should
    be close to its rigid-body value (|ATF| ~ 1)."""
    if not opt.seismic or not fnums:
        return
    lst.section(f"Low-frequency check (G-19): |ATF| at f_1 = {fnums[0] * opt.df:.6g} Hz vs the rigid-body value")
    skipped: Dict[str, List[str]] = {}
    for c in cases:
        why = getattr(c, "lowfreq_skip", "")
        if why:                                    # incoherent single modes k > 1, ME zones (FILE77)
            skipped.setdefault(why, []).append(c.outfile)
            continue
        h0 = np.array([rigid_body_anchor(int(k), c.cm, opt.ang) for k in md.eq_dof])
        big = (np.abs(h0) >= 0.1) & (md.eq_dof <= 3)
        if not big.any():
            continue
        dev = np.abs(np.abs(H[c.outfile][0]) - np.abs(h0))
        rel = np.where(big, dev / np.where(big, np.abs(h0), 1.0), 0.0)
        k = int(np.argmax(rel))
        lst.write(f" {c.outfile}: largest deviation {100 * rel[k]:.3f} % at node {int(md.eq_node[k])} "
                  f"{DOF_LABELS[int(md.eq_dof[k])]} ({big.sum()} DOFs checked)")
        bad = np.flatnonzero(rel > LOWFREQ_TOL)
        if bad.size:
            lst.warning(f"EDU-18 (G-19): {c.outfile}: |ATF| at the first SSI frequency deviates more than "
                        f"{100 * LOWFREQ_TOL:.0f} % from the rigid-body value at {bad.size} DOFs, e.g. "
                        + ", ".join(f"{int(md.eq_node[i])} {DOF_LABELS[int(md.eq_dof[i])]} ({abs(H[c.outfile][0, i]):.3f})"
                                    for i in bad[:10]) + "; check the model and the lowest SSI frequency")
    for why, names in skipped.items():
        lst.write(f" {', '.join(names[:6])}{' ...' if len(names) > 6 else ''}: not checked ({why})")


def _list_factors(lst, opt: Options, md: ModelData, cases: List[LoadCase], fs: FactorSet) -> None:
    """Listing of the FILE77 factors used (requirements 4.6 item 5) and the FFM recommendation (D-INC-12)."""
    lst.section("Incoherency, wave passage and multiple excitation (FILE77, requirements 4.6 item 5)")
    first = fs.files[cases[0].f77].meta
    what = {"AS": "deterministic algebraic sum of the incoherent modes (AS)",
            "SINGLE": f"single incoherent mode {int(first.get('mode', 0) or 0)}",
            "SS": "stochastic simulation samples", "COHERENT": "coherent motion with wave passage / ME factors",
            "BUILT": "per-level factors combined by BUILDFILE77"}.get(str(first.get("method", "")), str(first.get("method")))
    lst.write(f" HOUSE switches (deck)          : coherency {'incoherent' if opt.coh else 'coherent'}, wave passage "
              f"{'on' if opt.wpass else 'off'}, multiple excitation {'on' if opt.me else 'off'}")
    lst.write(f" Factors                        : {what}; {first.get('coherency_model') or 'no coherency model'}")
    lst.write(f" Application                    : " + ("FFM, free-field motion: b = A_f^T X_ff (s .* U'_f)" if opt.ffm
                                                       else "FFL, free-field load: b = A_f^T (s .* X_ff U'_f) (default)"))
    names = list(fs.files)
    lst.write(f" Files                          : {names[0]}" + (f" .. {names[-1]} ({len(names)} samples)"
                                                                 if len(names) > 1 else ""))
    if any(c.sim for c in cases):
        lst.write(f" Stochastic cases               : {len(cases)} = {opt.simul} simulations x 3 directions "
                  f"({cases[0].outfile} .. {cases[-1].outfile}); response statistics over the samples are formed "
                  "after MOTION (e.g. AVERAGE of the spectra)")
    if opt.ffm and md.nint and np.any(md.int_iface > 1):
        lst.warning("D-INC-12: free-field motion (FFM) with embedded interaction nodes -- FFM is recommended for "
                    "surface foundations only (spec 05c A.5.9: incoherent motion with the coherent impedance of "
                    "the excavation creates artificial scattering); use FFL")


# ======================================================================================
# Module entry point
# ======================================================================================
def run(ctx: ModuleContext) -> int:
    lst = ctx.listing
    d = decks.read(ctx.deck_path, NAME)
    opt = read_options(d, lst)
    md = read_model(ctx, opt, lst)
    if md.nint == 0:
        if opt.seismic:
            raise ModuleError("the model has no interaction nodes: seismic input acts through the interaction "
                              "nodes (INT) -- define them or run a foundation vibration analysis")
        if opt.impe:
            raise ModuleError("global impedance requested but the model has no interaction nodes")
        lst.warning("no interaction nodes: the structure is analysed without soil (fixed-base vibration)")
    model_rules(opt, md)
    cases = factored_cases(ctx, opt, lst) if opt.factored else load_cases(opt)
    skipped_y = ""
    if md.dim == 1 and opt.seismic and any(c.label == "Y" for c in cases):
        cases = [c for c in cases if c.label != "Y"]
        skipped_y = ("2D model: the anti-plane Y case (FILE1Y, SH) is not analysed; <simul> = 1 runs the in-plane "
                     "X and Z cases only (FILE1X -> FILE8X, FILE1Z -> FILE8Z)")
    read_load_cases(ctx, opt, md, cases, lst)
    if md.dim == 1 and opt.seismic:
        for c in cases:
            err = S2.inplane_input_error(c.cm, opt.ang, c.infile)
            if err:
                raise ModuleError(err)
    factors = read_factors(ctx, opt, md, cases, lst) if opt.factored else None
    if opt.df <= 0:
        raise ModuleError("the frequency step is not defined (deck df, delt/nft or input files)")
    _list_input(lst, opt, md, cases, ctx)
    if skipped_y:
        lst.write(f" {skipped_y}")
    if not opt.seismic:
        translate_load_nodes(md, cases, lst)
    if factors is not None:
        _list_factors(lst, opt, md, cases, factors)

    # ---- frequency survey -------------------------------------------------------------------------
    fnums = frequency_list(opt, cases, lst)
    file3 = None
    coox = cootk = None
    if opt.mode == 0 and md.nint:
        file3 = _read(ctx, "FILE3", "FILE3", "POINT")
        _check_df("FILE3", file3, opt)
        err = S2.point_dimension_error(md.dim, int(file3.meta.get("dim", 2)))
        if err:
            raise ModuleError(err)
        _check_depths("FILE3", file3, md)
        li = set(int(v) for v in np.asarray(file3["load_iface"]))
        off = sorted(set(int(v) for v in md.int_iface) - li)
        if off:
            raise ModuleError(f"interaction nodes on interfaces {off} but FILE3 holds point-load solutions for "
                              f"interfaces {sorted(li)} -- increase POINT <layer> and re-run POINT")
    elif opt.mode > 0:
        coox = read_index(ctx, "COOX", md, opt)
        if opt.solver_mode == 3:
            cootk = read_index(ctx, "COOTK", md, opt)
    rows3 = survey(lst, opt, fnums, cases, file3, coox, cootk)
    if factors is not None:
        survey_factors(lst, fnums, factors)
    nF = len(fnums)
    inputs = [f"{ctx.model}.N4"] + list(dict.fromkeys(c.infile for c in cases)) \
        + (["FILE3"] if file3 is not None else []) + (["COOSK", "COOSM"] if opt.solver_mode in (1, 2) else []) \
        + (list(factors.files) if factors is not None else [])
    input_bytes = sum(ctx.path(nm).stat().st_size for nm in inputs if ctx.path(nm).exists())
    resource_check(lst, md, nF, len(cases), opt, input_bytes)
    if opt.opmode == 1:
        lst.write("")
        lst.write(" Data check only (<opmode> = 1): inputs and frequency survey checked; nothing solved or written")
        return 0
    if opt.save and opt.solver_mode == 3:
        lst.write(" <save> = 1 in a solver Mode 3 restart: the restart files are reused, nothing new to save")
    if opt.save and not md.hashes:
        raise ModuleError("FILE90 missing -- run HOUSE (restart files need its hashes)")

    # ---- frequency loop -------------------------------------------------------------------------------
    part = md.part
    xy, iface = md.int_xyz[:, :2], md.int_iface
    xyz_ff = S2.plane_coordinates(md.int_xyz, opt.yc) if md.dim == 1 else md.int_xyz   # free-field points
    H = {c.outfile: np.zeros((nF, md.neq), complex) for c in cases}
    T = ss.rigid_body_transform(md.int_xyz, (opt.xc, opt.yc, opt.zc)) if opt.impe else None
    KG = np.zeros((nF, 6, 6), complex) if opt.impe else None
    n3 = 3 * md.nint
    Xall = np.zeros((nF, n3, n3), complex) if (opt.impe and nF * n3 * n3 * 16 <= file11_x_limit()) else None
    recs: List[FreqRecord] = []
    written_x: List[Tuple[int, int, float]] = []
    written_tk: List[Tuple[int, int, float]] = []
    warned: set = set()
    hashes = dict(md.hashes)
    t_all = time.perf_counter()
    for q, fn in enumerate(fnums):
        if ctx.cancelled():
            raise ModuleError("run cancelled")
        t0 = time.perf_counter()
        order = q + 1
        f = fn * opt.df
        w = 2.0 * np.pi * f
        rec = FreqRecord(q=order, fnum=fn, f=f, nf=part.nf, nn=part.nn)
        # -- impedance X_ff (steps 2-3)
        if opt.mode == 0:
            try:
                if md.nint and md.reduced:                 # 2D in-plane block and / or SYMM image sums
                    X, rec.rcond_F = reduced_soil_impedance(file3, rows3[fn], md)
                elif md.nint:
                    X, rec.rcond_F = soil_impedance(file3, rows3[fn], xy, iface)
                else:
                    X = np.zeros((0, 0), complex)
            except (ValueError, FloatingPointError, ss.SolverError) as exc:
                raise ModuleError(f"frequency {f:.6g} Hz: soil impedance failed: {exc}") from None
            if opt.save:
                name = restart_names(order)[0]
                write_container(ctx.path(name), "COOX", {"X": X},
                                dict(hashes, fnum=fn, freq=f, order=order, model=ctx.model, rcond_F=rec.rcond_F),
                                module=NAME)
                written_x.append((fn, order, f))
        else:
            X = np.asarray(_restart_record(ctx, coox, fn, md, opt)["X"], complex)
            if X.shape != (n3, n3):
                raise ModuleError(f"COOX record of frequency number {fn} has shape {X.shape}, expected ({n3}, {n3})")
        # -- system and factorisation (steps 4, 6)
        C = None
        try:
            if opt.solver_mode in (1, 2):
                C = ss.dynamic_matrix(md.Ks, md.Ms, md.Ke, md.Me, w)
                fac = ss.factorize(C, ss.restrict(X, part), part, rcond_min=ss.RCOND_MIN)
            else:
                fac = ss.SSIFactor.from_arrays(_restart_record(ctx, cootk, fn, md, opt), md.neq)
                if fac.part.nf != part.nf or not np.array_equal(fac.part.f_eq, part.f_eq):
                    raise ModuleError(f"COOTK record of frequency number {fn} has another equation partition "
                                      "than FILE4 -- run an initiation with <save> = 1")
                rec.note = "factorisation reused (COOTK)"
            rec.path, rec.rcond_nn, rec.rcond_ff = fac.path, fac.rcond_nn, fac.rcond_ff
            if fac.note:
                rec.note = fac.note
                if fac.path == "full" and part.nf:
                    lst.warning(f"frequency {f:.6g} Hz: {fac.note}")
            # -- loads and solution (steps 5, 6)
            if opt.seismic:
                Up = np.stack([free_field_at_nodes(c.data, c.rows[fn], xyz_ff, iface, opt.ang, opt.xc, opt.yc)
                               .reshape(-1) for c in cases], axis=1)                    # (3 nInt, ncases)
                if md.symm:                                # the input must have the symmetry of SYMM (D-ANL-12)
                    for j, c in enumerate(cases):
                        check_symmetric_input(md, c, fn, opt, xyz_ff, Up[:, j])
                if factors is not None:                    # FFL / FFM with the FILE77 factors (4.6 item 5)
                    s = np.stack([factors.factor(c, fn) for c in cases], axis=1)       # (nInt, ncases)
                    U = fac.solve(ss.incoherent_seismic_load(X, Up, s, part, ffm=bool(opt.ffm)))
                else:
                    U = fac.solve(ss.seismic_load(X, Up, part))
            else:
                b = vibration_loads(md, cases, fn, lst, warned)
                U = fac.solve(b[part.f_eq], b[part.n_eq])
        except ss.SolverError as exc:
            raise ModuleError(f"frequency {f:.6g} Hz (number {fn}): {exc}{_singular_hint(md, C)}") from None
        except ValueError as exc:
            raise ModuleError(f"frequency {f:.6g} Hz (number {fn}): {exc}") from None
        if not np.all(np.isfinite(U)):
            raise ModuleError(f"frequency {f:.6g} Hz: non-finite transfer functions")
        for j, c in enumerate(cases):
            H[c.outfile][q] = U[:, j]
        if opt.impe:
            KG[q] = ss.global_impedance(X, T)
            if Xall is not None:
                Xall[q] = X
        # X_ff is no longer needed: release it before the factorisation is written, so that at most
        # two dense matrices of this frequency are alive at any time (D-ANL-10)
        X = C = None
        if opt.save and opt.solver_mode in (1, 2):
            name = restart_names(order)[1]
            write_container(ctx.path(name), "COOTK", fac.arrays(),
                            dict(hashes, fnum=fn, freq=f, order=order, model=ctx.model, nEq=md.neq,
                                 path=fac.path, rcond_nn=fac.rcond_nn, rcond_ff=fac.rcond_ff),
                            module=NAME)
            written_tk.append((fn, order, f))
        fac = None                                   # the next frequency starts without dense arrays
        rec.seconds = time.perf_counter() - t0
        recs.append(rec)
        ctx.progress((q + 1) / nF, f"ANALYS: frequency {q + 1}/{nF} ({f:.4g} Hz, {rec.path})")
    t_all = time.perf_counter() - t_all

    # ---- outputs ------------------------------------------------------------------------------------
    _list_summary(lst, recs, opt.mode == 0)
    lst.write(f" total solution time {t_all:.3f} s for {nF} frequencies")
    files: List[str] = []
    fn_arr = np.asarray(fnums, np.int64)
    for c in cases:
        meta = {"df": opt.df, "type": opt.type, "case": c.label, "ang": opt.ang,
                "cm": c.cm if opt.seismic else 0, "nfft": opt.nft, "delt": opt.delt, "model_hash": md.model_hash,
                "mode": opt.mode, "simul": opt.simul, "source": c.infile, "model": ctx.model, "title": opt.title,
                "method": md.method, "xc": opt.xc, "yc": opt.yc, "zc": opt.zc,
                "content": ("seismic transfer functions U/U_cp (total motion per unit control motion)" if opt.seismic
                            else "displacement per unit load factor (vibration)"),
                "x_dim": md.dim}
        if md.symm:
            meta["x_symm"] = SYM.planes_array(md.symm).tolist()
        if factors is not None:
            m77 = factors.files[c.f77].meta
            meta.update(coh=opt.coh, wpass=opt.wpass, me=opt.me, ffm=opt.ffm, file77=c.f77, sim=c.sim,
                        x_incoherency=str(m77.get("method", "")), x_mode=int(m77.get("mode", 0) or 0),
                        x_stochastic=int(m77.get("stochastic", 0) or 0),
                        x_file77_direction="XYZ"[c.f77_dir if c.f77_dir >= 0 else c.cm])
        arrays = {"fnum": fn_arr, "freq": fn_arr * opt.df, "eq_node": md.eq_node, "eq_dof": md.eq_dof,
                  "H": H[c.outfile], "x_paths": np.array([PATH_CODE[r.path] for r in recs], np.int64)}
        probs = validate(Container("FILE8", meta, arrays))
        if probs:                                          # pragma: no cover - programming error
            raise ModuleError("FILE8: " + "; ".join(probs))
        write_container(ctx.path(c.outfile), "FILE8", arrays, meta, module=NAME)
        files.append(c.outfile)
    if written_x:
        a = np.array(written_x)
        write_container(ctx.path("COOXI"), "COOXI", {"fnum": a[:, 0].astype(np.int64),
                                                     "order": a[:, 1].astype(np.int64), "freq": a[:, 2]},
                        dict(hashes, df=opt.df, model=ctx.model), module=NAME)
        files += [restart_names(int(o))[0] for o in a[:, 1]] + ["COOXI"]
    if written_tk:
        a = np.array(written_tk)
        write_container(ctx.path("COOTKI"), "COOTKI", {"fnum": a[:, 0].astype(np.int64),
                                                       "order": a[:, 1].astype(np.int64), "freq": a[:, 2]},
                        dict(hashes, df=opt.df, model=ctx.model), module=NAME)
        files += [restart_names(int(o))[1] for o in a[:, 1]] + ["COOTKI"]
    if opt.impe:
        embedded = bool(md.nint and np.any(md.int_iface > 1))
        _print_impedance(lst, opt, fnums, KG, embedded)
        files += write_impedance_files(ctx, opt, md, fnums, KG, T, Xall)
    _lowfreq_check(lst, opt, md, cases, H, fnums)
    printed = [c for c in cases if c.sim <= 1]            # stochastic runs: the first simulation only
    for c in printed:
        _print_tfs(lst, opt, md, fnums, c, H[c.outfile])
    if len(printed) < len(cases):
        lst.write(f" transfer functions of simulations 2 .. {opt.simul} are not printed (FILE8 holds them all)")
    if opt.delrst:
        if opt.mode == 0:
            lst.warning("Delete Restart Files applies to restart runs only (<mode> > 0): ignored")
        elif opt.save:
            lst.warning("Delete Restart Files with Save Restart Files: the restart files are kept")
        else:
            gone = _delete_restart_files(ctx)
            lst.write(f" restart files deleted: {', '.join(gone) if gone else 'none found'}")
    lst.section("Files written")
    shown = files if len(files) <= 40 else files[:20] + ["..."] + files[-5:]
    for name in shown:
        lst.write(f" {name}")
    ctx.progress(1.0, "ANALYS done")
    return 0



if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
