"""STRESS: element stress, strain and force transfer functions, time histories and maxima.

What the module does (requirements section 4.10, spec 05d section 1, decisions D-STR-01 ... D-STR-12)
-------------------------------------------------------------------------------------------------
ANALYS solved the SSI problem at the SSI frequencies and stored the nodal transfer functions (TFs)
``H`` of every equation in FILE8; HOUSE stored, per output-capable element, a *recovery operator*
``S`` in FILE4 (``<model>.N4``, arrays ``rec_<T>_S``).  STRESS

1. computes the **stress transfer functions** at the SSI frequencies, ``STF_c(f_j) = S_c . U_e(f_j)``
   (:func:`sassi.core.stress_lib.element_stf`): SOLID/PLANE centroid stresses ``D* B u`` (and strains
   ``B u`` with ``EDUOPT,STRAINOUT,1``, D-STR-04) in global axes, SHELL membrane *stresses* and
   moments per unit length in the local axes x'y' (D-STR-10, D-W1-14), TSHELL membrane forces NXX NYY
   NXY and transverse shear forces QXZ QYZ per unit length and moments MXX MYY MXY per unit length in
   the local axes (spec 05d, D-STR-10; :mod:`sassi.elements.tshell`), BEAMS local end forces
   ``k_L T u`` = forces exerted *on the element* at I and J (D-STR-05), SPRING ``k* (u_J - u_I)``
   per global component.  Rigid-body motion produces no strain, so the total-motion TFs of a
   seismic FILE8 are used directly; the rigid-body part (the control motion projected on each
   element DOF) is subtracted first, ``S . (U_e - r_e)``, which is the same number in exact
   arithmetic but keeps the low-frequency STFs accurate (``H ~ 1 + O((f/f_n)^2)``);
2. **interpolates** each STF -- not the nodal TFs (D-STR-01) -- onto the Fourier grid with the
   MOTION schemes (``<interopt>`` 0-6, STRESSX smoothing and phase adjustment;
   :func:`sassi.core.interp.interpolate_tf`).  Below the first SSI frequency a seismic STF goes
   linearly to 0 at f = 0 (a rigid-body motion has no stress), above the last one it is 0;
3. **convolves** with the control *displacement* spectrum ``U_g = -g A / w^2`` (seismic: FILE8 holds
   ``U/U_cp``, the same ratio for displacement and acceleration, so a quantity linear in
   displacement needs the ground displacement; D-CNV-06, D-STR-02) or with the reference load
   spectrum ``F`` (foundation vibration), and returns to the time domain with ``irfft``;
4. computes the **derived** quantities in the time domain: the SOLID octahedral shear stress
   ``SOCT`` (all six stress components are computed and their maxima listed whenever it is
   requested) and, with strain output, the octahedral shear strain ``EOCT``;
5. writes, per element request (EOUT code digit per component: 0 none, 1 print maximum, 2 print
   maximum and save the history, D-STR-03): the listing tables of maxima with their times,
   ``<etype>_<ggg>_<eeeee>_<comp>.THS`` (``<save>`` = 1, every ``<skip>``-th sample, D-STR-07),
   ``.TFU`` / ``.TFI`` (``<itran>`` = 1, all element types, D-STR-06), FILE14 (formatted STFs at the
   SSI frequencies, ``<itran>`` = 1) and FILE15 (formatted histories of the code-2 components,
   written whatever ``<save>``: ``<save>`` controls the ``.THS`` files only);
6. all-element post-processing (STRESSX, SECDATAOPT; tier P1 basic versions):
   ``ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`` and ``<model>_ABS_MAX.sig/.tau/.bdsig/.bdtau``
   (``<savemax>``), all-element histories ``<model>.sig/.tau/.bdsig/.bdtau`` (``<saveth>``),
   ``NSTRESS/ESTRESS_nnnnn.ess`` frames (``<secdataopt>`` = 1, D-STR-12; nnnnn = 1-based time step,
   every ``<skip>``-th step, listed in ``NSTRESS/ESTRESS.lst`` for CALCSECTHIST), nodal stress frames
   ``NSTRESS/stress_<t>_<n>_<cls>`` (``<rstns>``, D-STR-08) and soil-pressure frames
   ``SOILPRES/pres_<t>_<n>_<ele|nod>`` with ``pres_max_ele``, ``pres_max_nod`` and
   ``STATIC_SOIL_PRESSURES.TXT`` (``<rstsp>``, D-STR-09); frame and soil-group selection from
   ``Frames.txt``.

7. non-linear soil SSI (``<iter>`` = 1, "Auto Computation of Strains in Soil El.", requirements
   4.10 item 5, :mod:`sassi.core.nlsoil`): for every nonlinear soil element of FILE78 (written by
   HOUSE), the strain-component histories (FILE4 strain operators through the same interpolation and
   convolution), the shear-strain measure gamma(t) per ISTR, ``gamma_eff = ESF max|gamma(t)|`` and
   the strain-compatible G and beta of the FILE73 curve -> FILE74 of this input direction (EOUT
   requests are not needed for these elements).  FILE8 must be the response to the FILE78
   properties (its FILE4 hash): a FILE8 of another HOUSE run stops the run.

8. TSHELL elements (manual sections 6 and 9.19): for every requested TSHELL element the maxima of the
   8 basic components are listed and written to ``TSHELL_ELEMENT_MAX.TXT`` (the input of THSHLSMH);
   with THSHLSTR = 1 the top/bottom face stresses and strains are computed from those maxima with the
   sign permutations of D-TSH-01 (:func:`sassi.elements.tshell.face_stresses`) -> listing and
   ``TSHELL_FACE_STRESSES.TXT``.  The flag is read from a ``thshlstr`` deck line when present (the deck
   schema has no such parameter yet), else from ``THSHLSTR.opt`` written by the THSHLSTR command
   (:func:`thshlstr_flag`).  As in the manual, no text frame files are generated for TSHELL elements.

Not available in this build: the binary database ``<model>_STRESS.bin`` (P2).

File layouts defined here (the manual does not give them, spec 05d OQ-S9; documented so that the
section-cut tools and the GUI read them back):

* ``<model>_ABS_MAX.<cls>``: ELEMENT_CENTER layout (D-FIL-06) with the class columns -- ``sig``
  normal stresses (SOLID SXX SYY SZZ; SHELL membrane FXX FYY 0; PLANE SXX SZZ 0), ``tau`` shear
  stresses (SOLID SXY SXZ SYZ; SHELL FXY 0 0; PLANE TXZ 0 0), ``bdsig`` SHELL bending moments per
  unit length (MXX MYY) and ``bdtau`` the SHELL twisting moment (MXY);
* ``<model>.<cls>``: ``#`` header lines (title, dt, column names) then one row per output time step
  ``t v1 v2 ...`` with one column per element component (``numpy.loadtxt`` reads it);
* element soil-pressure frames ``pres_*_ele``: frame header ``nrows ncols`` then rows
  ``group element p`` (D-FIL-03 with the element identified by its group and number);
* nodal stress frames ``NSTRESS/stress_*_<cls>`` (D-FIL-03 rows ``node v1 v2 ...``): for a structure of
  SHELL elements only, ``sig`` FXX FYY 0, ``tau`` FXY 0 0, ``bdsig`` MXX MYY, ``bdtau`` MXY in the element
  local axes; for a structure with SOLID elements, ``sig`` SXX SYY SZZ and ``tau`` SXY SXZ SYZ in global
  axes, the SHELL membrane stresses being rotated to global axes before the averaging (excavated soil
  elements are never part of the structure).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from .. import conventions as C
from ..core import interp as I
from ..core import nlsoil as NL
from ..core import signal as S
from ..core import stress_lib as SL
from ..elements import ELEMENTS
from ..elements import shell as SHELL_ELEMENT
from ..elements import tshell as TSHELL_ELEMENT
from ..io import decks
from ..io import library as LIB
from ..io.container import Container, read_container
from ..io.deckfmt import read_raw
from ..io.files import validate
from .base import ModuleContext, ModuleError, batch_main
from .house import file4_name
from .motion import (DF_RTOL, MAX_NFFT, MAX_SSI_FREQ, ControlMotion, find_file_nocase, frame_time,
                     read_control_motion, read_file8, write_frame, write_history_file, write_tf_file)

NAME = "STRESS"

BATCH_COLS = 384                 # STF columns (elements x components) interpolated/convolved together
MAX_LISTED = 20                  # elements listed in a warning
MEMORY_WARN = 2.0e9              # bytes of stored all-element histories that trigger a warning
DERIVED = {"SOCT": "octahedral shear stress"}
CODE_NAMES = {0: "no request", 1: "print maximum", 2: "print maximum and save time history"}

#: columns of the all-element component classes as indices into the ELEMENT_CENTER columns (-1 = 0)
CLASS_CENTER_IDX: Dict[str, Dict[str, Tuple[int, ...]]] = {
    "sig": {"SOLID": (0, 1, 2), "SHELL": (0, 1, -1), "PLANE": (0, 1, -1)},
    "tau": {"SOLID": (3, 4, 5), "SHELL": (2, -1, -1), "PLANE": (2, -1, -1)},
    "bdsig": {"SHELL": (3, 4)},
    "bdtau": {"SHELL": (5,)},
}
CLASS_TITLES = {"sig": "normal stresses (SOLID SXX SYY SZZ; SHELL membrane FXX FYY; PLANE SXX SZZ)",
                "tau": "shear stresses (SOLID SXY SXZ SYZ; SHELL membrane FXY; PLANE TXZ)",
                "bdsig": "SHELL bending moments per unit length (MXX MYY)",
                "bdtau": "SHELL twisting moment per unit length (MXY)"}


# =======================================================================================
# FILE4: element table and recovery operators
# =======================================================================================
@dataclass
class TypeRecovery:
    """Recovery data of one element type (FILE4 ``rec_<T>_*``, ARCHITECTURE section 4)."""
    etype: str
    comps: List[str]
    S: np.ndarray                    # (nE, nc, nd) complex
    eq: np.ndarray                   # (nE, nd) FILE4 equation numbers, -1 fixed / padding
    group: np.ndarray                # (nE,)
    elem: np.ndarray                 # (nE,)
    nodes: np.ndarray                # (nE, 8) node ids, 0 = unused
    excav: np.ndarray                # (nE,) 1 excavated soil
    B: Optional[np.ndarray] = None   # (nE, ns, nd) strain operator (SOLID, PLANE)
    props: Optional[Dict[str, np.ndarray]] = None   # FILE4 x_rec_<T>_<name> per-element data (TSHELL thick, E, nu)

    @property
    def n(self) -> int:
        return len(self.elem)


@dataclass
class File4Data:
    meta: dict
    dim: int
    gravity: float
    node_xyz: Dict[int, np.ndarray]
    eq_node: np.ndarray
    eq_dof: np.ndarray
    types: Dict[str, TypeRecovery]
    lookup: Dict[Tuple[int, int], Tuple[str, int]]       # (group, element) -> (type, row in rec_<T>)
    group_type: Dict[int, str]


def read_file4(path: Path) -> File4Data:
    """Read FILE4 (``<model>.N4``) and the recovery operators of every output-capable type."""
    try:
        c = read_container(path, kind="FILE4")
    except (ValueError, FileNotFoundError) as exc:
        raise ModuleError(f"{Path(path).name}: {exc}") from exc
    probs = validate(c)
    if probs:
        raise ModuleError(f"{Path(path).name} does not satisfy the FILE4 schema: " + "; ".join(probs))
    comps = dict(c.meta.get("components") or {})
    eg, ei = np.asarray(c["elem_group"], np.int64), np.asarray(c["elem_id"], np.int64)
    enodes = np.asarray(c["elem_nodes"], np.int64).reshape(len(eg), -1)
    excav = np.asarray(c["elem_excav"], np.int64)
    types: Dict[str, TypeRecovery] = {}
    lookup: Dict[Tuple[int, int], Tuple[str, int]] = {}
    group_type: Dict[int, str] = {}
    for name, clist in comps.items():
        if f"rec_{name}_S" not in c:
            continue
        Sx = np.asarray(c[f"rec_{name}_S"], dtype=complex)
        eq = np.asarray(c[f"rec_{name}_eq"], dtype=np.int64)
        idx = np.asarray(c[f"rec_{name}_idx"], dtype=np.int64)
        if Sx.ndim != 3 or eq.shape != (Sx.shape[0], Sx.shape[2]) or idx.shape != (Sx.shape[0],):
            raise ModuleError(f"{Path(path).name}: recovery arrays of {name} have inconsistent shapes")
        if Sx.shape[1] != len(clist):
            raise ModuleError(f"{Path(path).name}: {name} recovery has {Sx.shape[1]} rows for {len(clist)} "
                              "components")
        B = c.get(f"x_rec_{name}_B")
        pre = f"x_rec_{name}_"
        props = {k[len(pre):]: np.asarray(c[k]) for k in c.arrays if k.startswith(pre) and k != f"{pre}B"
                 and np.asarray(c[k]).shape[:1] == (Sx.shape[0],)}
        t = TypeRecovery(name, list(clist), Sx, eq, eg[idx], ei[idx], enodes[idx], excav[idx],
                         None if B is None else np.asarray(B, dtype=complex), props or None)
        types[name] = t
        for r, (g, e) in enumerate(zip(t.group, t.elem)):
            lookup[(int(g), int(e))] = (name, r)
            group_type[int(g)] = name
    xyz = np.asarray(c["node_xyz"], dtype=float).reshape(-1, 3)
    return File4Data(meta=dict(c.meta), dim=int(c.meta.get("dim", 2)), gravity=float(c.meta.get("gravity", 0.0)),
                     node_xyz={int(n): p for n, p in zip(c["node_id"], xyz)},
                     eq_node=np.asarray(c["eq_node"], np.int64), eq_dof=np.asarray(c["eq_dof"], np.int64),
                     types=types, lookup=lookup, group_type=group_type)


# =======================================================================================
# Element output requests (EOUT, D-STR-03; CHECK Errors 79-82)
# =======================================================================================
@dataclass
class ElemRequest:
    group: int
    elem: int
    etype: str
    row: int                         # row in rec_<T>
    codes: np.ndarray                # (n slots of the type,) 0/1/2


def parse_requests(rows: List[dict], f4: File4Data, lst) -> Dict[Tuple[int, int], ElemRequest]:
    """EOUT table rows -> requests keyed by (group, element).

    Errors 80 (group without output-capable elements in FILE4: unknown group, GENERAL) and
    81 (element not in its group) stop the run; Error 82 (element requested twice) is merged with
    the larger code per component and a warning (CHECK reports it before AFWRITE)."""
    out: Dict[Tuple[int, int], ElemRequest] = {}
    bad_g: List[int] = []
    bad_e: List[Tuple[int, int]] = []
    dups: List[Tuple[int, int]] = []
    extra: List[Tuple[int, int]] = []
    for r in rows:
        g, e = int(r["group"]), int(r["element"])
        codes = np.array([int(r[f"c{k}"]) for k in range(1, 13)], dtype=np.int64)
        if np.any((codes < 0) | (codes > 2)):
            raise ModuleError(f"EOUT group {g} element {e}: output codes must be 0, 1 or 2 (D-STR-03)")
        if g not in f4.group_type:
            if g not in bad_g:
                bad_g.append(g)
            continue
        key = (g, e)
        if key not in f4.lookup:
            bad_e.append(key)
            continue
        etype, row = f4.lookup[key]
        nc = len(f4.types[etype].comps)
        if np.any(codes[nc:] != 0):
            extra.append(key)
        q = ElemRequest(g, e, etype, row, codes[:nc].copy())
        if key in out:
            dups.append(key)
            out[key].codes = np.maximum(out[key].codes, q.codes)
        else:
            out[key] = q
    errs = []
    if bad_g:
        errs.append("Error 80: illegal group for output request: " + ", ".join(map(str, bad_g))
                    + " (no output-capable elements in FILE4: unknown group, or GENERAL matrices, which have no "
                      "STRESS output)")
    if bad_e:
        errs.append("Error 81: illegal element output request: " + ", ".join(
            f"element {e} group {g}" for g, e in bad_e[:MAX_LISTED]) + (" ..." if len(bad_e) > MAX_LISTED else ""))
    if errs:
        for m in errs:
            lst.error(m)
        raise ModuleError("; ".join(errs))
    if dups:
        lst.warning("Error 82: element output requests defined more than once, merged (larger code per "
                    "component): " + ", ".join(f"element {e} group {g}" for g, e in dups[:MAX_LISTED]))
    if extra:
        lst.warning("output codes beyond the components of the element type are ignored for "
                    + ", ".join(f"element {e} group {g}" for g, e in extra[:MAX_LISTED]))
    return out


# =======================================================================================
# Output components of one element type
# =======================================================================================
@dataclass
class OutComp:
    """One output quantity of an element type."""
    name: str                        # file-name code (SXX, MZJ, SOCT, EXX, EOCT ...)
    slot: int                        # EOUT digit (0-based) that requests it
    col: int                         # column in the history array (-1: derived)
    derived: str = ""                # 'oct_stress' / 'oct_strain'
    strain: bool = False


@dataclass
class TypePlan:
    """Columns processed for one element type: primary stress components (rows of S), strains
    (rows of B, D-STR-04) and the derived quantities."""
    etype: str
    rec: TypeRecovery
    stress_rows: np.ndarray          # rows of S that are primary components
    strain_rows: np.ndarray          # rows of B used (empty without STRAINOUT)
    outs: List[OutComp]
    center_cols: List[int]           # history columns of the six ELEMENT_CENTER columns (-1 zero)

    @property
    def ncols(self) -> int:
        return len(self.stress_rows) + len(self.strain_rows)


def type_plan(rec: TypeRecovery, strainout: bool) -> TypePlan:
    comps = rec.comps
    prim = [i for i, c in enumerate(comps) if c not in DERIVED]
    outs = [OutComp(comps[i], i, k) for k, i in enumerate(prim)]
    for i, c in enumerate(comps):
        if c == "SOCT":
            outs.append(OutComp(c, i, -1, derived="oct_stress"))
    strain_rows = np.zeros(0, dtype=np.int64)
    if strainout and rec.B is not None and rec.etype in ("SOLID", "PLANE"):
        ns = rec.B.shape[1]
        strain_rows = np.arange(ns)
        names = SL.strain_component_names(rec.etype, [comps[i] for i in prim[:ns]])
        outs += [OutComp(nm, prim[k], len(prim) + k, strain=True) for k, nm in enumerate(names)]
        if "SOCT" in comps:
            outs.append(OutComp("EOCT", comps.index("SOCT"), -1, derived="oct_strain", strain=True))
    center = []
    if rec.etype in SL.CENTER_COLUMNS:
        pos = {comps[i]: k for k, i in enumerate(prim)}
        center = [pos[c] if c is not None else -1 for c in SL.CENTER_COLUMNS[rec.etype]]
    return TypePlan(rec.etype, rec, np.asarray(prim, dtype=np.int64), strain_rows, outs, center)


# =======================================================================================
# Deck checks
# =======================================================================================
def _validate_deck(d, lst) -> None:
    """Runtime checks of the STRESS deck (spec 05d section 1.14; CHECK repeats them before AFWRITE)."""
    if d["interopt"] not in I.OPTION_NAMES:
        raise ModuleError(f"interpolation option {d['interopt']} must be 0..6")
    if d["pzadj"] not in (0, 1):
        raise ModuleError(f"phase adjustment {d['pzadj']} must be 0 or 1")
    if d["smo"] < 0:
        raise ModuleError("smoothing parameter must be >= 0")
    if d["skip"] < 0:
        raise ModuleError("Error 67: negative output time-history step (Skip Time History Steps)")
    if d["dur"] < 0:
        raise ModuleError("Error 68: negative total duration")
    if d["type"] not in (0, 1):
        raise ModuleError(f"type of analysis {d['type']} must be 0 (seismic) or 1 (vibration)")
    for k in ("iter", "save", "itran", "savemax", "saveth", "rstns", "rstsp", "secdataopt"):
        if int(d[k]) not in (0, 1):
            raise ModuleError(f"STRESS <{k}> = {d[k]} must be 0 or 1")
    nft = int(d["nft"])
    if nft <= 0:
        raise ModuleError("Error 50: illegal number of values for the Fourier transform")
    if not C.is_power_of_two(nft):
        raise ModuleError(f"NFFT = {nft} is not a power of 2 (AFWRITE rounds it, Warning 9)")
    if nft > MAX_NFFT:
        raise ModuleError(f"NFFT = {nft} exceeds the STRESS maximum {MAX_NFFT} (D-CNV-11)")
    if d["delt"] <= 0:
        raise ModuleError("Error 49: illegal time step of the control motion")
    if d["smo"] > 0 and d["interopt"] == 6:
        lst.warning("EDU-14: smoothing has no effect with interpolation option 6 (spline); ignored")
    if d["pzadj"] == 1 and d["interopt"] != 6 and d["smo"] == 0:
        lst.warning("phase adjustment with S = 0 leaves the phase unchanged (rho = 1); it needs S > 0")


def _strainout(deck_path: Path) -> bool:
    """``EDUOPT,STRAINOUT,1`` (D-STR-04).  The STRESS deck schema has no ``strainout`` parameter
    (lead-owned ``sassi.io.decks``); a ``strainout`` line in the deck is honoured when present."""
    try:
        raw = read_raw(deck_path)
    except Exception:  # pragma: no cover - the deck was read successfully just before
        return False
    v = str(raw.params.get("strainout", "0")).strip().strip('"')
    try:
        return int(float(v)) == 1
    except ValueError:
        return False


#: sidecar written by the THSHLSTR command in the model directory (sassi.prep.commands.thickshell)
THSHLSTR_FILE = "THSHLSTR.opt"
#: STRESS outputs of the TSHELL elements (requested elements; maxima of the 8 basic components / faces)
TSHELL_MAX_FILE = "TSHELL_ELEMENT_MAX.TXT"
TSHELL_FACE_FILE = "TSHELL_FACE_STRESSES.TXT"


def thshlstr_flag(deck_path: Optional[Path], workdir: Path) -> Tuple[int, str]:
    """THSHLSTR face-stress flag (manual 9.19.2, D-TSH-01) and where it was found.

    The STRESS deck carries ``thshlstr`` (AFWRITE writes 0 or 1 from the model's THSHLSTR record) and is
    then authoritative.  A deck value of -1 (not set: decks built without AFWRITE) or a deck without the
    parameter falls back to the ``THSHLSTR.opt`` file that the THSHLSTR command writes in the model
    directory, otherwise 0 (manual default)."""
    if deck_path is not None:
        try:
            raw = read_raw(deck_path)
            v = raw.params.get("thshlstr")
            # The deck (AFWRITE copies the model's THSHLSTR record) is authoritative.
            if v is not None and str(v).strip().strip('"') != "" and int(float(str(v).strip().strip('"'))) in (0, 1):
                return int(float(str(v).strip().strip('"'))), "STRESS deck"
        except Exception:  # pragma: no cover - the deck was read successfully before
            pass
    p = find_file_nocase(workdir, THSHLSTR_FILE)
    if p is not None:
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            t = line.split("!", 1)[0].strip()
            if t.upper().startswith("THSHLSTR"):
                parts = [x.strip() for x in t.split(",")]
                try:
                    return (1 if len(parts) > 1 and int(float(parts[1] or 0)) == 1 else 0), p.name
                except ValueError:
                    return 0, f"{p.name} (unreadable: {t!r}; 0 used)"
    return 0, "default"


# =======================================================================================
# The run
# =======================================================================================
@dataclass
class _Flags:
    savemax: bool
    saveth: bool
    rstns: bool
    rstsp: bool
    ess: bool

    @property
    def all_elements(self) -> bool:
        return self.savemax or self.saveth or self.ess

    @property
    def any(self) -> bool:
        return self.savemax or self.saveth or self.rstns or self.rstsp or self.ess


def run(ctx: ModuleContext) -> int:
    """STRESS module entry point (requirements section 4.10)."""
    lst = ctx.listing
    d = decks.read(ctx.deck_path, NAME)
    if d["title"]:
        lst.write(f" {d['title']}")
    _validate_deck(d, lst)
    strainout = _strainout(ctx.deck_path)
    seismic = int(d["type"]) == 0
    mode = "seismic" if seismic else "vibration"
    nfft, dt = int(d["nft"]), float(d["delt"])
    df = 1.0 / (nfft * dt)
    option, smooth, pzadj = int(d["interopt"]), float(d["smo"]), int(d["pzadj"])
    flags = _Flags(bool(d["savemax"]), bool(d["saveth"]), bool(d["rstns"]), bool(d["rstsp"]), bool(d["secdataopt"]))

    lst.section("STRESS options")
    lst.write(f"   operation mode            : {'data check' if d['opmode'] else 'complete solution'}")
    lst.write(f"   type of analysis          : {mode}")
    lst.write(f"   interpolation option      : {option} ({I.OPTION_NAMES[option]})")
    lst.write(f"   smoothing S, phase adj.   : {smooth:g}, {pzadj}")
    lst.write(f"   NFFT, dt, df              : {nfft}, {dt:g} s, {df:.9g} Hz (Fourier period {nfft * dt:g} s)")
    lst.write(f"   save histories (.THS)     : {int(d['save'])};  transfer functions (.TFU/.TFI, FILE14): "
              f"{int(d['itran'])};  skip {int(d['skip'])}")
    lst.write(f"   {'control motion' if seismic else 'load history  '}            : {d['thfile']}  (fopt {d['fopt']}, "
              f"records {d['rec1']}..{d['rec2'] or 'last'}, mult {d['mult']:g}, max {d['max']:g})")
    lst.write(f"   all-element options       : save max {int(flags.savemax)}, save histories {int(flags.saveth)}, "
              f"nodal stress frames {int(flags.rstns)}, soil pressures {int(flags.rstsp)}, .ess frames "
              f"{int(flags.ess)}")
    if strainout:
        lst.write("   strain output             : on (EDUOPT,STRAINOUT,1; SOLID/PLANE strains, files E..., D-STR-04)")
    if int(d["iter"]):
        lst.write("   non-linear soil strains   : on (<iter> = 1: FILE78 elements -> FILE74, requirements 4.10 item 5)")
    thshlstr, th_src = thshlstr_flag(ctx.deck_path, ctx.workdir)

    # ---------------------------------------------------------------- FILE4 and requests
    f4 = read_file4(ctx.require(file4_name(ctx.model), "HOUSE"))
    reqs = parse_requests(d.rows("eout"), f4, lst)
    nlin = _nonlinear_input(ctx, f4, lst) if int(d["iter"]) else None
    if not reqs and not flags.any and nlin is None:
        raise ModuleError("Error 79: no element output request (EOUT) and no all-element option")
    lst.section("FILE4 elements")
    for t, rec in f4.types.items():
        lst.write(f"   {t:<7s}: {rec.n} elements in groups {', '.join(map(str, sorted(set(rec.group.tolist()))))}"
                  f"  (components {' '.join(rec.comps)})")
    if reqs:
        _list_requests(lst, reqs, f4)
    if "TSHELL" in f4.types:
        lst.write(f"   TSHELL face stresses/strains (THSHLSTR): {thshlstr} ({th_src}; 0 = the 8 basic components only, "
                  "1 = also top/bottom face values from their maxima, D-TSH-01)")
        if flags.any:
            lst.warning("no text frame files are generated for the TSHELL elements (manual Ch. 3): they are not "
                        "part of the all-element outputs (ELEMENT_CENTER, .sig/.tau, ESTRESS, NSTRESS frames)")

    frames: List[int] = []
    pgroups: List[int] = []
    if flags.rstns or flags.rstsp:
        p = find_file_nocase(ctx.workdir, "Frames.txt")
        if p is None:
            raise ModuleError("restart for nodal stress / soil pressure contours requested but Frames.txt is "
                              "missing (Frame Selection, D-FIL-07)")
        try:
            frames, pgroups = SL.read_frames_file(p)
        except ValueError as exc:
            raise ModuleError(str(exc)) from None
        lst.write(f"   Frames.txt: {len(frames)} frames, soil-pressure groups {pgroups or 'none'}")

    # ---------------------------------------------------------------- control motion
    thpath = LIB.module_path(d["thfile"], ctx.workdir) if str(d["thfile"] or "").strip() else None
    if thpath is None:
        raise ModuleError("Error 73: no time-history file (THFILE) given")
    cmot = read_control_motion(thpath, int(d["fopt"]), int(d["rec1"]), int(d["rec2"]), dt, nfft,
                               float(d["mult"]), float(d["max"]))
    lst.section("Control motion" if seismic else "Reference load history")
    if LIB.is_library_name(d["thfile"]):
        lst.write(f"   {LIB.note(d['thfile'])}")
    imax = int(np.argmax(np.abs(cmot.acc))) if len(cmot.acc) else 0
    lst.write(f"   file {cmot.path.name}: {cmot.nrec} records, {len(cmot.acc)} used, duration "
              f"{len(cmot.acc) * dt:g} s, quiet zone {(nfft - len(cmot.acc)) * dt:g} s")
    lst.write(f"   peak value after scaling {np.max(np.abs(cmot.acc)) if len(cmot.acc) else 0:.6g} "
              f"{'g' if seismic else '(load factor)'} at t = {imax * dt:g} s")

    # ---------------------------------------------------------------- FILE8
    f8name = d["file8"] or "FILE8"
    f8 = read_file8(ctx.require(f8name, "ANALYS (or COMBIN)"), f8name)
    f_ssi = _check_file8(f8, f8name, f4, d, df, lst)
    if nlin is not None:
        _nonlinear_file8(f8, f8name, nlin[0], lst)
    if not np.all(np.isfinite(np.asarray(f8["H"]))):
        raise ModuleError(f"{f8name} holds non-finite transfer functions (NaN/Inf): check the ANALYS run")
    for t, rec in f4.types.items():
        if not np.all(np.isfinite(rec.S)):
            raise ModuleError(f"FILE4 recovery operators of {t} hold non-finite values: check the HOUSE run")
    gravity = float(d["gravity"])
    if f4.gravity > 0 and abs(f4.gravity - gravity) > 1e-9 * max(gravity, 1.0):
        lst.warning(f"STRESS gravity {gravity:g} differs from the HOUSE gravity {f4.gravity:g} in FILE4; the "
                    "deck value is used for U_g = -g A / w^2")
    if int(d["opmode"]) == 1:
        lst.section("Data check")
        lst.write("   input checked; data-check mode writes no output files")
        return 0

    cm = int(f8.meta.get("cm", d["cm"]))
    ang = float(f8.meta.get("ang", d["ang"]))
    if cm != int(d["cm"]) or abs(ang - float(d["ang"])) > 1e-9:
        lst.warning(f"control direction/angle of {f8name} (cm {cm}, ang {ang:g}) differ from the deck "
                    f"(cm {d['cm']}, ang {d['ang']:g}); the {f8name} values are used")
    sr = _StressRun(ctx, d, f4, np.asarray(f8["H"]), f_ssi, cmot, strainout, flags, frames, pgroups, mode, cm, ang,
                    thshlstr=bool(thshlstr))
    if nlin is not None:
        _nonlinear_strains(sr, nlin[0], nlin[1], cm, str(f8.meta.get("case", "")))
    sr.execute(reqs)
    return 0


def _check_file8(f8: Container, f8name: str, f4: File4Data, d, df: float, lst) -> np.ndarray:
    """Consistency of FILE8 with the STRESS deck and FILE4; returns the SSI frequencies ``fnum df``."""
    meta = f8.meta
    fnum = np.asarray(f8["fnum"], dtype=np.int64)
    df8 = float(meta["df"])
    if abs(df - df8) > DF_RTOL * df8:
        raise ModuleError(f"frequency step of STRESS 1/(NFFT dt) = {df:.9g} Hz differs from the {f8name} step "
                          f"{df8:.9g} Hz (D-MOT-14): use the SITE NFFT and time step")
    if int(meta["type"]) != int(d["type"]):
        raise ModuleError(f"type of analysis {d['type']} does not match {f8name} (type {meta['type']})")
    if (not np.array_equal(np.asarray(f8["eq_node"]), f4.eq_node)
            or not np.array_equal(np.asarray(f8["eq_dof"]), f4.eq_dof)):
        raise ModuleError(f"{f8name} and FILE4 have different equation maps (eq_node/eq_dof): HOUSE was rerun "
                          "after ANALYS -- rerun ANALYS")
    h4, h8 = f4.meta.get("model_hash"), meta.get("model_hash")
    if h4 and h8 and h4 != h8:
        lst.warning(f"{f8name} model_hash differs from FILE4 (the TFs may come from another HOUSE run)")
    if fnum[-1] > int(d["nft"]) // 2:
        lst.warning(f"G-02: frequency numbers above NFFT/2 = {int(d['nft']) // 2}; interpolation stops at the "
                    "Nyquist frequency (D-CNV-12)")
    if len(fnum) > MAX_SSI_FREQ:
        lst.warning(f"EDU-02: {len(fnum)} SSI frequencies exceed the STRESS input limit {MAX_SSI_FREQ} "
                    "(requirements 4.0.1); the run continues")
    f_ssi = fnum * df
    lst.section(f"Transfer functions ({f8name})")
    lst.write(f"   {len(fnum)} SSI frequencies {f_ssi[0]:.6g} .. {f_ssi[-1]:.6g} Hz, {len(f4.eq_node)} equations")
    option = int(d["interopt"])
    eff = I.effective_option(option, len(fnum))
    if eff != option:
        lst.warning(f"{len(fnum)} SSI frequencies: interpolation option {option} replaced by "
                    f"{'complex linear' if eff < 0 else 'option 6'}")
    return f_ssi


def _list_requests(lst, reqs: Dict[Tuple[int, int], ElemRequest], f4: File4Data) -> None:
    lst.section("Element output requests (EOUT; code 0 none, 1 maximum, 2 maximum + history)")
    by_group: Dict[int, List[ElemRequest]] = {}
    for q in reqs.values():
        by_group.setdefault(q.group, []).append(q)
    for g in sorted(by_group):
        qs = sorted(by_group[g], key=lambda q: q.elem)
        comps = f4.types[qs[0].etype].comps
        lst.write(f"   group {g} ({qs[0].etype}), {len(qs)} elements; components {' '.join(comps)}")
        codes: Dict[str, List[int]] = {}
        for q in qs:
            codes.setdefault("".join(str(int(c)) for c in q.codes), []).append(q.elem)
        for code, els in codes.items():
            lst.write(f"      code {code}: elements {_ranges(els)}")


def _ranges(ids: Sequence[int]) -> str:
    ids = sorted(set(int(i) for i in ids))
    out, i = [], 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        out.append(str(ids[i]) if i == j else f"{ids[i]}-{ids[j]}")
        i = j + 1
    return ", ".join(out)


# =======================================================================================
# Processing
# =======================================================================================
class _StressRun:
    """STF, interpolation, convolution and outputs of one STRESS run."""

    def __init__(self, ctx: ModuleContext, d, f4: File4Data, H: np.ndarray, f_ssi: np.ndarray, cmot: ControlMotion,
                 strainout: bool, flags: _Flags, frames: List[int], pgroups: List[int], mode: str, cm: int = 0,
                 ang: float = 0.0, thshlstr: bool = False):
        self.ctx, self.d, self.f4, self.H, self.f_ssi, self.cmot = ctx, d, f4, H, f_ssi, cmot
        self.thshlstr = bool(thshlstr)
        self.lst = ctx.listing
        self.wd = ctx.workdir
        self.flags, self.frames, self.pgroups = flags, frames, pgroups
        self.mode = mode
        self.seismic = mode == "seismic"
        self.nfft, self.dt = int(d["nft"]), float(d["delt"])
        self.option = int(d["interopt"])
        self.smooth = 0.0 if self.option == 6 else float(d["smo"])
        self.pzadj = int(d["pzadj"])
        self.skip = max(1, int(d["skip"]))
        self.save = bool(d["save"])
        self.itran = bool(d["itran"])
        self.gravity = float(d["gravity"])
        self.f_grid = S.fourier_grid(self.nfft, self.dt)
        fN = min(f_ssi[-1], self.f_grid[-1])
        self.kmax = int(np.searchsorted(self.f_grid, fN * (1 + 1e-12), side="right"))
        self.nout = S.output_length(self.nfft, self.dt, float(d["dur"]))
        self.steps = SL.output_steps(self.nout, self.skip)
        # convolution kernel: control displacement U_g = -g A / w^2 (seismic) or the load spectrum F
        self.drive = (S.displacement_spectrum(cmot.A, self.f_grid, self.gravity) if self.seismic else cmot.A)
        self.plans = {t: type_plan(rec, strainout) for t, rec in f4.types.items()}
        # seismic rigid-body motion of the element DOFs, subtracted before the recovery (exact, see
        # sassi.core.stress_lib.element_stf): avoids the low-frequency cancellation of total-motion TFs
        self.rigid: Dict[str, Optional[np.ndarray]] = {}
        for t, rec in f4.types.items():
            spec = ELEMENTS.get(C.ELEMENT_TYPE_CODES.get(t, -1))
            ok = self.seismic and spec is not None and spec.nnodes * len(spec.dofs) == rec.S.shape[2]
            self.rigid[t] = SL.rigid_body_dofs(spec.dofs, spec.nnodes, cm, ang) if ok else None
        self.maxima: Dict[Tuple[int, int, str], Tuple[float, float, bool]] = {}   # (g, e, comp) -> (max, t, aux)
        self.center_max: Dict[str, np.ndarray] = {}
        self.center_hist: Dict[str, np.ndarray] = {}
        self.processed: Dict[str, np.ndarray] = {}
        self.files: List[str] = []

    # ------------------------------------------------------------------ strain histories (non-linear soil)
    def strain_histories(self, rec: TypeRecovery, rows: np.ndarray) -> np.ndarray:
        """Strain-component histories ``(NFFT, len(rows), ns)`` of SOLID/PLANE elements ``rows``: the rows
        of the strain operator ``B`` (EXX EYY EZZ GXY GXZ GYZ / EXX EZZ GXZ) through the same STF ->
        interpolation -> convolution chain as the stresses (:meth:`_batch`), over the whole Fourier
        period (requirements 4.10 items 1-3 and 5)."""
        if rec.B is None:
            raise ModuleError(f"FILE4 holds no strain operators for {rec.etype} elements: re-run HOUSE")
        nF, b = len(self.f_ssi), len(rows)
        stf = SL.element_stf(rec.B[rows], rec.eq[rows], self.H, self.rigid.get(rec.etype))     # (nF, b, ns)
        ns = stf.shape[2]
        Hg = I.interpolate_tf(self.f_ssi, stf.reshape(nF, b * ns), self.f_grid, self.option, smooth=self.smooth,
                              pzadj=self.pzadj, h0=None, mode=self.mode)
        Hg = S.hermitian_bins(Hg, self.nfft) if self.seismic else _nyquist_real(Hg, self.nfft)
        R = S.hermitian_bins(Hg * self.drive[:, None], self.nfft)
        return np.fft.irfft(R, n=self.nfft, axis=0).reshape(self.nfft, b, ns)

    # ------------------------------------------------------------------ work sets
    def _work_rows(self, reqs: Dict[Tuple[int, int], ElemRequest]) -> Tuple[Dict[str, np.ndarray], Set[str]]:
        rows: Dict[str, Set[int]] = {t: set() for t in self.f4.types}
        for q in reqs.values():
            rows[q.etype].add(q.row)
        store: Set[str] = set()
        fl = self.flags
        for t, rec in self.f4.types.items():
            if t not in SL.CENTER_COLUMNS:
                continue
            if fl.all_elements:
                rows[t].update(range(rec.n))
            if fl.saveth or fl.ess:
                store.add(t)
            if fl.rstns and t in ("SOLID", "SHELL") and self.f4.dim == 2:
                rows[t].update(int(r) for r in np.where(rec.excav == 0)[0])
                store.add(t)
            if fl.rstsp and t == "SOLID":
                rows[t].update(int(r) for r in np.where(np.isin(rec.group, self.pgroups))[0])
                store.add(t)
        return {t: np.array(sorted(r), dtype=np.int64) for t, r in rows.items() if r}, store

    # ------------------------------------------------------------------ main
    def execute(self, reqs: Dict[Tuple[int, int], ElemRequest]) -> None:
        lst = self.lst
        work, store = self._work_rows(reqs)
        nbytes = sum(self.nout * self.f4.types[t].n * 6 * 8 for t in store)
        if nbytes > MEMORY_WARN:
            lst.warning(f"all-element histories need {nbytes / 1e9:.1f} GB of memory (reduce NFFT / dur or the "
                        "all-element options)")
        lst.section("Output")
        lst.write(f"   {sum(len(r) for r in work.values())} elements processed; history length {self.nout} samples "
                  f"({self.nout * self.dt:g} s), output every {self.skip} sample(s); TFI rows 0..{self.kmax - 1}")
        req_by_type: Dict[str, Dict[int, ElemRequest]] = {}
        for q in reqs.values():
            req_by_type.setdefault(q.etype, {})[q.row] = q
        fh14 = open(self.wd / "FILE14", "w", encoding="utf-8") if (self.itran and reqs) else None
        need15 = any(np.any(q.codes == 2) for q in reqs.values())
        fh15 = open(self.wd / "FILE15", "w", encoding="utf-8") if need15 else None
        try:
            unit = "per unit control displacement" if self.seismic else "per unit load factor"
            for fh, what in ((fh14, f"stress transfer functions at the SSI frequencies (STF = S . U_e, {unit})"),
                             (fh15, "time histories of the code-2 components")):
                if fh is not None:
                    fh.write(f"# SASSI-EDU STRESS {what}; model {self.ctx.model}\n")
            total = sum(len(r) for r in work.values())
            done = 0
            for t, rows in work.items():
                plan = self.plans[t]
                rec = plan.rec
                if t in SL.CENTER_COLUMNS:
                    self.center_max[t] = np.zeros((rec.n, 6))
                if t in store:
                    self.center_hist[t] = np.zeros((self.nout, rec.n, 6))
                self.processed[t] = rows
                per = max(1, BATCH_COLS // max(plan.ncols, 1))
                for b0 in range(0, len(rows), per):
                    if self.ctx.cancelled():
                        raise ModuleError("run cancelled")
                    rb = rows[b0:b0 + per]
                    self._batch(plan, rb, req_by_type.get(t, {}), fh14, fh15)
                    done += len(rb)
                    self.ctx.progress(min(1.0, done / max(total, 1)), f"STRESS: {done}/{total} elements")
        finally:
            for fh in (fh14, fh15):
                if fh is not None:
                    fh.close()
        if fh14 is not None:
            self.files.append("FILE14")
        if fh15 is not None:
            self.files.append("FILE15")
        self._list_maxima(reqs)
        self._tshell_outputs(reqs)
        if self.flags.savemax:
            self._write_center_max()
        if self.flags.saveth:
            self._write_class_histories()
        if self.flags.ess:
            self._write_ess()
        if self.flags.rstns:
            self._write_stress_frames()
        if self.flags.rstsp:
            self._write_soil_pressures()
        nth = sum(1 for f in self.files if f.endswith(".THS"))
        ntf = sum(1 for f in self.files if f.endswith(".TFU"))
        lst.section("Files written")
        lst.write(f"   {nth} .THS histories, {ntf} .TFU/.TFI pairs")
        for f in self.files:
            if not f.endswith((".THS", ".TFU", ".TFI")):
                lst.write(f"   {f}")

    # ------------------------------------------------------------------ one batch of elements
    def _batch(self, plan: TypePlan, rows: np.ndarray, reqs: Dict[int, ElemRequest], fh14, fh15) -> None:
        rec = plan.rec
        nF = len(self.f_ssi)
        b = len(rows)
        rig = self.rigid.get(plan.etype)
        parts = [SL.element_stf(rec.S[rows][:, plan.stress_rows, :], rec.eq[rows], self.H, rig)]
        if len(plan.strain_rows):
            parts.append(SL.element_stf(rec.B[rows][:, plan.strain_rows, :], rec.eq[rows], self.H, rig))
        stf = np.concatenate(parts, axis=2)                                  # (nF, b, nc)
        nc = stf.shape[2]
        Hg = I.interpolate_tf(self.f_ssi, stf.reshape(nF, b * nc), self.f_grid, self.option, smooth=self.smooth,
                              pzadj=self.pzadj, h0=None, mode=self.mode)
        Hg = S.hermitian_bins(Hg, self.nfft) if self.seismic else _nyquist_real(Hg, self.nfft)
        R = S.hermitian_bins(Hg * self.drive[:, None], self.nfft)
        hist = np.fft.irfft(R, n=self.nfft, axis=0)[:self.nout].reshape(self.nout, b, nc)
        Hg = Hg.reshape(len(self.f_grid), b, nc)
        # derived quantities in the time domain (requirements 4.10 item 3)
        derived: Dict[str, np.ndarray] = {}
        for oc in plan.outs:
            if oc.derived == "oct_stress":
                derived[oc.name] = SL.octahedral_shear_stress(*(hist[:, :, k] for k in range(6)))
            elif oc.derived == "oct_strain":
                s0 = len(plan.stress_rows)
                derived[oc.name] = SL.octahedral_shear_strain(*(hist[:, :, s0 + k] for k in range(6)))
        # all-element arrays
        if plan.center_cols and rec.etype in self.center_max:
            cc = np.array(plan.center_cols)
            vals = np.where(cc[None, None, :] >= 0, hist[:, :, np.clip(cc, 0, None)], 0.0)   # (nout, b, 6)
            self.center_max[rec.etype][rows] = np.max(np.abs(vals), axis=0)
            if rec.etype in self.center_hist:
                self.center_hist[rec.etype][:, rows, :] = vals
        # per-request outputs
        for j, r in enumerate(rows):
            q = reqs.get(int(r))
            if q is None:
                continue
            self._element_outputs(plan, q, j, stf, Hg, hist, derived, fh14, fh15)

    def _element_outputs(self, plan: TypePlan, q: ElemRequest, j: int, stf, Hg, hist, derived, fh14, fh15) -> None:
        wd = self.wd
        et = plan.etype
        need_aux = set()
        if et == "TSHELL" and np.any(q.codes > 0):
            # all 8 basic maxima: TSHELL_ELEMENT_MAX.TXT and the THSHLSTR face values (D-TSH-01)
            need_aux.update(o.name for o in plan.outs if o.col >= 0 and not o.strain)
        for oc in plan.outs:
            if oc.derived == "oct_stress" and q.codes[oc.slot] > 0:
                need_aux.update(o.name for o in plan.outs if o.col >= 0 and not o.strain)
            if oc.derived == "oct_strain" and q.codes[oc.slot] > 0:
                need_aux.update(o.name for o in plan.outs if o.col >= 0 and o.strain)
        for oc in plan.outs:
            code = int(q.codes[oc.slot])
            aux = code == 0 and oc.name in need_aux
            if code == 0 and not aux:
                continue
            r = derived[oc.name][:, j] if oc.col < 0 else hist[:, j, oc.col]
            k = int(np.argmax(np.abs(r)))
            self.maxima[(q.group, q.elem, oc.name)] = (float(abs(r[k])), k * self.dt, aux)
            if aux:
                continue
            if code == 2:
                rs = r[self.steps]
                if self.save:
                    nm = C.element_result_name(et, q.group, q.elem, oc.name, "THS")
                    write_history_file(wd / nm, rs, self.dt * self.skip)
                    self.files.append(nm)
                if fh15 is not None:
                    fh15.write(f"# {et} group {q.group} element {q.elem} {oc.name}: dt {self.dt * self.skip:g} s, "
                               f"{len(rs)} values; columns time value\n")
                    fh15.write("\n".join(f"{s * self.dt:.6f} {v:.10e}" for s, v in zip(self.steps, rs)) + "\n")
            if self.itran and oc.col >= 0:
                hdr = (f"{et} group {q.group} element {q.elem} {oc.name} "
                       + ("stress TF per unit control displacement" if self.seismic else "TF per unit load factor"))
                nm_u = C.element_result_name(et, q.group, q.elem, oc.name, "TFU")
                nm_i = C.element_result_name(et, q.group, q.elem, oc.name, "TFI")
                write_tf_file(wd / nm_u, self.f_ssi, stf[:, j, oc.col], True, hdr)
                write_tf_file(wd / nm_i, self.f_grid[:self.kmax], Hg[:self.kmax, j, oc.col], True,
                              hdr + f" interp={self.option}")
                self.files += [nm_u, nm_i]
                if fh14 is not None:
                    fh14.write(f"# {et} group {q.group} element {q.elem} {oc.name}: f(Hz) |STF| phase(deg) Re Im\n")
                    for fi, h in zip(self.f_ssi, stf[:, j, oc.col]):
                        fh14.write(f"{fi:.6f} {abs(h):.10e} {np.degrees(np.angle(h)):.6f} {h.real:.10e} "
                                   f"{h.imag:.10e}\n")

    # ------------------------------------------------------------------ listing of maxima
    def _list_maxima(self, reqs: Dict[Tuple[int, int], ElemRequest]) -> None:
        if not self.maxima:
            return
        lst = self.lst
        unit = "stress, force or moment" + (" (strain: dimensionless)" if any(
            o.strain for p in self.plans.values() for o in p.outs) else "")
        by_group: Dict[int, List[Tuple[int, str, float, float, bool]]] = {}
        order = {t: [o.name for o in p.outs] for t, p in self.plans.items()}
        for (g, e, comp), (v, t, aux) in self.maxima.items():
            by_group.setdefault(g, []).append((e, comp, v, t, aux))
        for g in sorted(by_group):
            et = self.f4.group_type[g]
            lst.section(f"Maximum absolute response, group {g} ({et}) -- {unit}")
            lst.write(f"{'element':>9s} {'comp':>6s} {'max |r|':>16s} {'time (s)':>10s}")
            for e, comp, v, t, aux in sorted(by_group[g], key=lambda r: (r[0], order[et].index(r[1]))):
                lst.write(f"{e:9d} {comp:>6s} {v:16.6e} {t:10.4f}" + (("   (computed for the TSHELL maxima/faces, not "
                                                                       "requested)" if et == "TSHELL" else
                                                                       "   (computed for SOCT/EOCT, not requested)")
                                                                      if aux else ""))
        big = max(((v, g, e, c, t) for (g, e, c), (v, t, aux) in self.maxima.items() if not aux), default=None)
        if big:
            lst.write("")
            lst.write(f"   largest requested value {big[0]:.6e}: group {big[1]} element {big[2]} {big[3]} at "
                      f"t = {big[4]:.4f} s")
        if self.itran and any(o.derived for p in self.plans.values() for o in p.outs):
            lst.write("   note: SOCT/EOCT are computed in the time domain and have no transfer function files")

    # ------------------------------------------------------------------ TSHELL maxima and THSHLSTR faces
    def _tshell_outputs(self, reqs: Dict[Tuple[int, int], ElemRequest]) -> None:
        """Requested TSHELL elements: ``TSHELL_ELEMENT_MAX.TXT`` (maxima of the 8 basic components, the
        input of THSHLSMH) and, with THSHLSTR = 1, the top/bottom face stresses and strains computed from
        those maxima (``TSHELL_FACE_STRESSES.TXT`` and the listing; D-TSH-01,
        :func:`sassi.elements.tshell.face_stresses`).

        ``TSHELL_ELEMENT_MAX.TXT``: ``#`` header lines, then one row ``group element NXX NYY NXY QXZ QYZ MXX
        MYY MXY`` (absolute maxima over the output duration, local element axes).
        ``TSHELL_FACE_STRESSES.TXT``: ``#`` header lines, then four rows per element, one per sign
        permutation ``(s_N, s_M)`` = ``++ -- +- -+``: ``group element perm SXX SYY TXY S1 S2 EXX EYY GXY E1 E2
        TXZ TYZ`` (stresses F/L^2, strains dimensionless; TXZ/TYZ = 1.5 Q/t at the mid-surface)."""
        req_t = sorted((q for q in reqs.values() if q.etype == "TSHELL"), key=lambda q: (q.group, q.elem))
        # an EOUT row whose 8 codes are all 0 requests nothing: no maxima, no row (no NaN)
        qs = [q for q in req_t if np.any(np.asarray(q.codes) > 0)]
        if len(qs) < len(req_t):
            skipped = [q for q in req_t if not np.any(np.asarray(q.codes) > 0)]
            self.lst.write(f"   {len(skipped)} TSHELL element request(s) with all output codes 0 (no output): not in "
                           f"{TSHELL_MAX_FILE}" + (f" / {TSHELL_FACE_FILE}" if self.thshlstr else "") + ", e.g. "
                           + ", ".join(f"group {q.group} element {q.elem}" for q in skipped[:5]))
        self._retire_tshell_files(self._tshell_files(qs) if qs else [])

    def _retire_tshell_files(self, written) -> None:
        """A TSHELL output file of an earlier run that this run does not rewrite (no TSHELL output requests,
        or THSHLSTR = 0 now) would be read by THSHLSMH as current: renamed ``<name>.prev`` (as FILE78)."""
        for name in (TSHELL_MAX_FILE, TSHELL_FACE_FILE):
            if name in written:
                continue
            p = find_file_nocase(self.wd, name)
            if p is not None:
                q = p.with_name(p.name + ".prev")
                p.replace(q)
                self.files.append(f"{q.name}  (renamed: {p.name} of an earlier run does not describe this run)")

    def _tshell_files(self, qs: List[ElemRequest]) -> List[str]:
        """Write ``TSHELL_ELEMENT_MAX.TXT`` and, with THSHLSTR = 1, ``TSHELL_FACE_STRESSES.TXT`` for the
        TSHELL requests ``qs`` (all with maxima); returns the names written."""
        comps = self.f4.types["TSHELL"].comps
        rows = []
        for q in qs:
            vals = [self.maxima.get((q.group, q.elem, c), (float("nan"), 0.0, True))[0] for c in comps]
            rows.append((q, np.array(vals, dtype=float)))
        lines = ["# SASSI-EDU STRESS: maxima of the TSHELL basic components (absolute values over the output "
                 f"duration, local element axes); model {self.ctx.model}",
                 "# N and Q per unit length (F/L), M per unit length (F.L/L)",
                 "# group element " + " ".join(comps)]
        lines += [f"{q.group} {q.elem} " + " ".join(f"{v:.10e}" for v in vals) for q, vals in rows]
        (self.wd / TSHELL_MAX_FILE).write_text("\n".join(lines) + "\n", encoding="utf-8")
        self.files.append(f"{TSHELL_MAX_FILE}  ({len(rows)} TSHELL elements; input of THSHLSMH)")
        if not self.thshlstr:
            return [TSHELL_MAX_FILE]
        rec = self.f4.types["TSHELL"]
        props = rec.props or {}
        if not all(k in props for k in ("thick", "E", "nu")):
            self.lst.warning("THSHLSTR = 1 but FILE4 has no TSHELL thickness/material data (x_rec_TSHELL_*): re-run "
                             "HOUSE; no face stresses written")
            return [TSHELL_MAX_FILE]
        lst = self.lst
        lst.section("TSHELL face stresses and strains (THSHLSTR,1; from the maxima of the basic components, D-TSH-01)")
        lst.write("   sigma = s_N N/t + s_M 6 M/t^2 for (NXX, MXX) and (NYY, MYY) with the sign permutations (s_N, s_M); "
                  "TXY = NXY/t + 6 MXY/t^2 (++);")
        lst.write("   top face z' = +t/2: (s_N, s_M); bottom face: (s_N, -s_M).  Principal values from each permutation; "
                  "strains by plane-stress Hooke's law.")
        hdr = f"{'group':>6s}{'elem':>7s}{'perm':>5s}" + "".join(f"{c:>13s}" for c in TSHELL_ELEMENT.FACE_COLUMNS)
        lst.write(hdr)
        out = ["# SASSI-EDU STRESS: TSHELL face stresses and strains from the maxima of the basic components "
               f"(THSHLSTR,1; D-TSH-01); model {self.ctx.model}",
               "# SXX = s_N NXX/t + s_M 6 MXX/t^2, SYY likewise, TXY = NXY/t + 6 MXY/t^2; S1/S2 principal stresses; "
               "EXX EYY GXY E1 E2 plane-stress strains; TXZ TYZ = 1.5 Q/t (mid-surface)",
               "# group element perm " + " ".join(TSHELL_ELEMENT.FACE_COLUMNS) + " TXZ TYZ"]
        for q, vals in rows:
            r = q.row
            fs = TSHELL_ELEMENT.face_stresses(vals, float(props["thick"][r]), float(props["E"][r]),
                                              float(props["nu"][r]))
            for name, _, _ in TSHELL_ELEMENT.FACE_PERMUTATIONS:
                v = fs["rows"][name]
                lst.write(f"{q.group:>6d}{q.elem:>7d}{name:>5s}" + "".join(f"{x:>13.5e}" for x in v))
                out.append(f"{q.group} {q.elem} {name} " + " ".join(f"{x:.10e}" for x in v)
                           + f" {fs['TXZ']:.10e} {fs['TYZ']:.10e}")
            lst.write(f"{'':>13s}transverse shear at the mid-surface: TXZ = {fs['TXZ']:.5e}, TYZ = {fs['TYZ']:.5e}")
        (self.wd / TSHELL_FACE_FILE).write_text("\n".join(out) + "\n", encoding="utf-8")
        self.files.append(f"{TSHELL_FACE_FILE}  ({len(rows)} TSHELL elements x 4 sign permutations)")
        return [TSHELL_MAX_FILE, TSHELL_FACE_FILE]

    # ------------------------------------------------------------------ all-element files
    def _center_blocks(self, values: Dict[str, np.ndarray], types: Sequence[str] = SL.CENTER_TYPES,
                       cols: Optional[Dict[str, Tuple[int, ...]]] = None) -> List[SL.CenterBlock]:
        """Group blocks (ascending group number) of the processed elements of ``types``."""
        groups = []
        for t in types:
            if t in values:
                rec = self.f4.types[t]
                groups += [(int(g), t) for g in set(rec.group[self.processed[t]].tolist())]
        ordered = SL.ordered_group_numbers(groups)
        blocks = []
        for g, t in sorted(groups):
            rec = self.f4.types[t]
            rows = self.processed[t][rec.group[self.processed[t]] == g]
            rows = rows[np.argsort(rec.elem[rows], kind="stable")]
            v = values[t][rows]
            if cols is not None:
                idx = np.array(cols[t])
                v = np.where(idx[None, :] >= 0, v[:, np.clip(idx, 0, None)], 0.0)
            blocks.append(SL.CenterBlock(t, g, ordered[g], rec.elem[rows], v))
        return blocks

    def _write_center_max(self) -> None:
        blocks = self._center_blocks(self.center_max)
        if not blocks:
            self.lst.warning("Save Max Value: no SOLID, SHELL, PLANE or SPRING elements in FILE4")
            return
        SL.write_element_center(self.wd / "ELEMENT_CENTER_ABS_MAX_STRESSES.TXT", blocks)
        self.files.append("ELEMENT_CENTER_ABS_MAX_STRESSES.TXT  (BEAMS excluded: end forces, spec 05d 1.9)")
        for cls, spec in CLASS_CENTER_IDX.items():
            bl = self._center_blocks(self.center_max, tuple(spec), spec)
            if bl:
                nm = f"{self.ctx.model}_ABS_MAX.{cls}"
                SL.write_element_center(self.wd / nm, bl)
                self.files.append(f"{nm}  ({CLASS_TITLES[cls]})")

    def _write_class_histories(self) -> None:
        for cls, spec in CLASS_CENTER_IDX.items():
            cols, names = [], []
            for t in spec:
                if t not in self.center_hist:
                    continue
                rec = self.f4.types[t]
                rows = self.processed[t]
                rows = rows[np.lexsort((rec.elem[rows], rec.group[rows]))]
                comps = SL.CENTER_COLUMNS[t]
                for r in rows:
                    for k in spec[t]:
                        if k < 0:
                            continue
                        cols.append(self.center_hist[t][self.steps, r, k])
                        names.append(f"{t}_{int(rec.group[r]):03d}_{int(rec.elem[r]):05d}_{comps[k]}")
            if not cols:
                continue
            nm = f"{self.ctx.model}.{cls}"
            M = np.column_stack([self.steps * self.dt] + cols)
            hdr = (f"SASSI-EDU STRESS all-element histories: {CLASS_TITLES[cls]}\n"
                   f"dt = {self.dt * self.skip:.10g} s, {len(self.steps)} steps, {len(cols)} element components\n"
                   + "time " + " ".join(names))
            np.savetxt(self.wd / nm, M, fmt="%.10e", header=hdr)
            self.files.append(f"{nm}  ({len(cols)} element components x {len(self.steps)} steps)")

    def _write_ess(self) -> None:
        if not self.center_hist:
            self.lst.warning("SECDATAOPT = 1: no SOLID, SHELL, PLANE or SPRING elements; no .ess frames")
            return
        out = self.wd / "NSTRESS"
        out.mkdir(exist_ok=True)
        names = []
        for s in self.steps:
            blocks = self._center_blocks({t: v[s] for t, v in self.center_hist.items()})
            nm = f"ESTRESS_{s + 1:05d}.ess"
            SL.write_element_center(out / nm, blocks, fmt="{:.10e}")
            names.append(nm)
        (out / "ESTRESS.lst").write_text(
            f"# ESTRESS frames of model {self.ctx.model}: dt = {self.dt:g} s, frame n is time (n-1) dt\n"
            + "\n".join(names) + "\n", encoding="utf-8")
        self.files.append(f"NSTRESS/ESTRESS_nnnnn.ess  ({len(names)} frames, list NSTRESS/ESTRESS.lst, D-STR-12)")

    def _frame_indices(self, what: str) -> List[int]:
        ok = [n for n in self.frames if 1 <= n <= self.nout]
        bad = [n for n in self.frames if not 1 <= n <= self.nout]
        if bad:
            self.lst.warning(f"{what}: frame numbers outside 1..{self.nout} ignored: {bad[:MAX_LISTED]}")
        return ok

    def _write_stress_frames(self) -> None:
        """Nodal stress frames (D-STR-08, spec 05d 1.9; plotting only): plain mean of the adjacent element-centre
        values of the *structural* (non-excavated) SOLID and SHELL elements of a 3D model.

        Excavated soil is not part of the structure, so the frame classes depend only on the structural element
        types actually averaged -- never on which all-element options made STRESS process the soil elements:

        * SHELL-only structure: membrane frames ``sig`` (FXX FYY 0) and ``tau`` (FXY 0 0) plus the bending frames
          ``bdsig`` (MXX MYY) and ``bdtau`` (MXY), in the element local axes x'y' (the ELEMENT_CENTER
          convention of the SHELL outputs);
        * structure with SOLID elements (with or without SHELL elements): membrane frames only, in **global
          axes**: ``sig`` = SXX SYY SZZ, ``tau`` = SXY SXZ SYZ.  The SHELL membrane stresses are first rotated
          to global axes (``sigma_g = Lam^T sigma_l Lam``, :func:`sassi.core.stress_lib.membrane_to_global`),
          so that the average at a node shared by a wall and a solid combines components of one coordinate
          system."""
        lst = self.lst
        if self.f4.dim != 2:
            lst.warning("restart for nodal stress contours: only SOLID and SHELL elements of 3D models (spec 05d "
                        "1.9); no frames written")
            return
        hists: List[Tuple[str, np.ndarray]] = []
        for t in ("SOLID", "SHELL"):
            if t in self.center_hist:
                rows = np.where(self.f4.types[t].excav == 0)[0]
                if len(rows):
                    hists.append((t, rows))
        if not hists:
            lst.warning("restart for nodal stress contours: no structural SOLID or SHELL elements")
            return
        struct = [t for t, _ in hists]
        elems = [self.f4.types[t].nodes[r] for t, rows in hists for r in rows]
        nodes, A = SL.averaging_matrix(elems)
        shell_only = struct == ["SHELL"]
        if shell_only:
            classes = ["sig", "tau", "bdsig", "bdtau"]
            rows = hists[0][1]
            v = self.center_hist["SHELL"][:, rows, :]
            E_cls = {}
            for cls in classes:
                idx = np.array(CLASS_CENTER_IDX[cls]["SHELL"])
                E_cls[cls] = np.where(idx[None, None, :] >= 0, v[:, :, np.clip(idx, 0, None)], 0.0)
            axes = "SHELL local axes x'y'"
        else:
            classes = ["sig", "tau"]
            G = np.concatenate([self.center_hist[t][:, rows, :] if t == "SOLID" else
                                SL.membrane_to_global(self.center_hist[t][:, rows, :3], self._shell_axes(rows))
                                for t, rows in hists], axis=1)                       # (nout, nE, 6) global
            E_cls = {"sig": G[:, :, :3], "tau": G[:, :, 3:]}
            axes = ("global axes; membrane stresses only (SOLID + SHELL, SHELL rotated to global)"
                    if len(struct) > 1 else "global axes")
        out = self.wd / "NSTRESS"
        out.mkdir(exist_ok=True)
        frames = self._frame_indices("Frames.txt (nodal stresses)")
        for cls in classes:
            Nn = SL.average_histories(A, E_cls[cls])
            for n in frames:
                write_frame(out / f"stress_{frame_time((n - 1) * self.dt)}_{n:05d}_{cls}", nodes, Nn[n - 1])
            write_frame(out / f"stress_ABS_MAX_{cls}", nodes, np.max(np.abs(Nn), axis=0))
        self.files.append(f"NSTRESS/stress_<t>_<n>_<cls>  ({len(frames)} frames x {', '.join(classes)}; "
                          f"{len(nodes)} nodes; {axes})")

    def _shell_axes(self, rows: np.ndarray) -> np.ndarray:
        """Local axes ``Lam (n, 3, 3)`` (rows x' y' z', spec 08 4.6) of the SHELL elements ``rows``, rebuilt from
        the FILE4 node slots and coordinates exactly as HOUSE built them (distinct corners in the order of
        their first appearance, so a triangle given as I J K K has the axes of I J K)."""
        rec = self.f4.types["SHELL"]
        Lam = np.empty((len(rows), 3, 3))
        for j, r in enumerate(rows):
            ns = list(dict.fromkeys(int(v) for v in rec.nodes[r][:4] if v > 0))
            Lam[j] = SHELL_ELEMENT.local_frame(np.array([self.f4.node_xyz[n] for n in ns]))[0]
        return Lam

    def _write_soil_pressures(self) -> None:
        """Soil pressures on the structure (D-STR-09): ``p = -n^T sigma n`` (compression positive) on the faces
        of the soil SOLID elements of the Frames.txt groups that are shared with the structure."""
        lst = self.lst
        rec = self.f4.types.get("SOLID")
        if rec is None or "SOLID" not in self.center_hist or not self.pgroups:
            lst.warning("restart for soil pressure contours: no soil-pressure SOLID groups (Frames.txt)")
            return
        bad = [g for g in self.pgroups if self.f4.group_type.get(g) != "SOLID"]
        if bad:
            lst.warning(f"soil-pressure groups that are not SOLID groups ignored: {bad}")
        soil = np.isin(rec.group, self.pgroups)
        # structure nodes: the DOF nodes of the non-excavated elements outside the soil-pressure groups.  Only
        # the first ELEMENTS[code].nnodes slots carry DOFs: the BEAMS K node (slot 3) only orients the cross
        # section and may be any soil or free-field node, so it must not make a soil face 'structural'.
        struct_nodes: Set[int] = set()
        for t, tr in self.f4.types.items():
            if t == "SPRING":
                continue
            spec = ELEMENTS.get(C.ELEMENT_TYPE_CODES.get(t, -1))
            nslot = spec.nnodes if spec is not None else tr.nodes.shape[1]
            mask = tr.excav == 0
            if t == "SOLID":
                mask &= ~soil
            for ns in tr.nodes[mask][:, :nslot]:
                struct_nodes.update(int(v) for v in ns if v > 0)
        rows, normals, face_nodes = [], [], []
        for r in np.where(soil)[0]:
            faces = SL.shared_faces(rec.nodes[r], struct_nodes)
            if not faces:
                continue
            rows.append(int(r))
            normals.append([SL.face_normal(np.array([self.f4.node_xyz[n] for n in f])) for f in faces])
            face_nodes.append(sorted(set(n for f in faces for n in f)))
        if not rows:
            lst.warning("soil-pressure groups have no SOLID face shared with the structure; no pressures")
            return
        sig = self.center_hist["SOLID"]
        p = np.zeros((self.nout, len(rows)))
        for j, r in enumerate(rows):
            p[:, j] = -np.mean([SL.normal_stress(sig[:, r, :], nv) for nv in normals[j]], axis=0)
        # static (geological) pressures, added algebraically
        groups = sorted(set(int(rec.group[r]) for r in rows))
        ordered = SL.ordered_group_numbers((g, "SOLID") for g in groups)
        static = np.zeros(len(rows))
        sp_path = find_file_nocase(self.wd, "STATIC_SOIL_PRESSURES.TXT")
        key = {(int(rec.group[r]), int(rec.elem[r])): j for j, r in enumerate(rows)}
        if sp_path is not None:
            try:
                for b in SL.read_element_center(sp_path):
                    for e, v in zip(b.elements, b.values):
                        j = key.get((b.group, int(e)))
                        if j is not None and v.size:
                            static[j] = float(v[0])
            except ValueError as exc:
                raise ModuleError(f"STATIC_SOIL_PRESSURES.TXT: {exc}") from None
            lst.write(f"   static soil pressures read from {sp_path.name} (max |p_static| {np.max(np.abs(static)):.6g})")
        else:
            SL.write_element_center(self.wd / "STATIC_SOIL_PRESSURES.TXT", self._pressure_blocks(rows, static, ordered))
            self.files.append("STATIC_SOIL_PRESSURES.TXT  (zeros; enter static pressures and rerun)")
        p = p + static[None, :]
        nodes, A = SL.averaging_matrix(face_nodes)
        out = self.wd / "SOILPRES"
        out.mkdir(exist_ok=True)
        grp = np.array([[rec.group[r], rec.elem[r]] for r in rows], dtype=np.int64)
        frames = self._frame_indices("Frames.txt (soil pressures)")
        for n in frames:
            tag = f"{frame_time((n - 1) * self.dt)}_{n:05d}"
            _write_element_frame(out / f"pres_{tag}_ele", grp, p[n - 1])
            write_frame(out / f"pres_{tag}_nod", nodes, A @ p[n - 1])
        pmax = np.max(np.abs(p), axis=0)
        _write_element_frame(out / "pres_ABS_MAX_ele", grp, pmax)
        write_frame(out / "pres_ABS_MAX_nod", nodes, A @ pmax)
        SL.write_element_center(self.wd / "pres_max_ele", self._pressure_blocks(rows, pmax, ordered))
        write_frame(self.wd / "pres_max_nod", nodes, A @ pmax)
        self.files.append(f"SOILPRES/pres_<t>_<n>_<ele|nod>  ({len(frames)} frames, {len(rows)} elements, "
                          f"{len(nodes)} nodes), pres_max_ele, pres_max_nod")
        lst.write(f"   soil pressure: {len(rows)} elements adjacent to the structure; max |p| {pmax.max():.6e} "
                  "(compression positive, static included)")

    def _pressure_blocks(self, rows: List[int], vals: np.ndarray, ordered: Dict[int, int]) -> List[SL.CenterBlock]:
        rec = self.f4.types["SOLID"]
        rows_a = np.asarray(rows)
        blocks = []
        for g in sorted(ordered):
            m = rec.group[rows_a] == g
            o = np.argsort(rec.elem[rows_a[m]], kind="stable")
            blocks.append(SL.CenterBlock("SOLID", g, ordered[g], rec.elem[rows_a[m]][o], vals[m][o][:, None]))
        return blocks


# =======================================================================================
# Non-linear soil SSI: effective strains -> FILE74 (requirements 4.10 item 5, D-NLS-03)
# =======================================================================================
def _nonlinear_input(ctx: ModuleContext, f4: File4Data, lst) -> Optional[Tuple[NL.NLFile, Dict[int, object]]]:
    """FILE78 (nonlinear elements, ESF, ISTR, curves, G_max -- written by HOUSE) and FILE73 (curves),
    checked against FILE4: every element must be a SOLID/PLANE element with a strain operator.  Without
    FILE78 (no non-linear SSI HOUSE run) there is nothing to compute: a warning, FILE74 is not written."""
    p78 = ctx.path("FILE78")
    if not p78.exists():
        lst.warning("<iter> = 1 (Auto Computation of Strains in Soil El.) but FILE78 does not exist (HOUSE was not run "
                    "with the Non-Linear SSI option and a .pin file): FILE74 is not written and the run continues")
        return None
    try:
        f78 = NL.read_nlfile(p78, "FILE78")
        curves = NL.load_curves(ctx.path("FILE73"))
    except NL.NLSoilError as exc:
        raise ModuleError(f"non-linear soil: {exc}") from None
    if not f78.rows:
        raise ModuleError("FILE78 lists no nonlinear soil elements")
    h4 = str(f4.meta.get("model_hash", ""))
    if f78.file4 and h4 and f78.file4 != h4:
        lst.warning("FILE78 was written by another HOUSE run than FILE4 (model hash): re-run HOUSE, ANALYS and STRESS "
                    "of the same iteration")
    bad = [r.key for r in f78.rows if r.key not in f4.lookup or f4.lookup[r.key][0] not in NL.STRAIN_COMPONENTS]
    if bad:
        raise ModuleError(f"FILE78 elements not found as SOLID/PLANE elements of FILE4 (group, element): {bad[:8]}")
    nocurve = sorted({r.curve for r in f78.rows if r.curve not in curves})
    if nocurve:
        raise ModuleError(f"FILE78 refers to curves {nocurve} that are not in FILE73 (curves {sorted(curves)})")
    try:
        NL.check_curves(curves, {r.curve for r in f78.rows})
    except NL.NLSoilError as exc:
        raise ModuleError(f"non-linear soil: {exc}") from None
    return f78, curves


def _nonlinear_file8(f8: Container, f8name: str, f78: NL.NLFile, lst) -> None:
    """<iter> = 1: the transfer functions must be the response to the FILE78 properties -- ANALYS run
    on the FILE4 HOUSE wrote with them (FILE8 ``model_hash`` = FILE78 FILE4 hash).  Otherwise the
    strains would be labelled with the FILE78 iteration although they belong to other properties and
    the next convergence measure would be false (e.g. HOUSE re-run, ANALYS skipped)."""
    h8 = str(f8.meta.get("model_hash", "") or "")
    if f78.file4 and h8:
        if h8 != f78.file4:
            raise ModuleError(f"non-linear soil (<iter> = 1): {f8name} is not the response to the properties of "
                              f"FILE78 (iteration {f78.iteration}): ANALYS computed it from another FILE4 than the one "
                              "HOUSE wrote with these properties -- run ANALYS (New Structure restart, <mode> 1) after "
                              "HOUSE, then STRESS")
    else:
        lst.warning(f"non-linear soil (<iter> = 1): FILE78 or {f8name} carries no FILE4 hash -- it cannot be checked "
                    f"that {f8name} is the response to the FILE78 properties (iteration {f78.iteration})")


def _nonlinear_strains(sr: "_StressRun", f78: NL.NLFile, curves, cm: int, case: str = "") -> None:
    """Effective strains of the nonlinear soil elements of FILE78 for this input direction -> FILE74.

    gamma(t) per ISTR of the element's group (:func:`sassi.core.nlsoil.shear_strain_history`) from
    the strain histories over the whole Fourier period; ``gamma_eff = ESF max|gamma(t)|``; G and beta
    from the FILE73 curve at gamma_eff (FILE74 also lists max|gamma(t)| and G_max)."""
    lst, f4 = sr.lst, sr.f4
    gmax_pct: Dict[Tuple[int, int], float] = {}
    work: Dict[Tuple[str, int], List[Tuple[Tuple[int, int], int]]] = {}
    for r in f78.rows:
        etype, row = f4.lookup[r.key]
        gi = f78.group_info(r.group)
        work.setdefault((etype, gi.istr if gi is not None else 0), []).append((r.key, row))
    for (etype, istr), items in work.items():
        rec = f4.types[etype]
        per = max(1, BATCH_COLS // len(NL.STRAIN_COMPONENTS[etype]))
        for b0 in range(0, len(items), per):
            if sr.ctx.cancelled():
                raise ModuleError("run cancelled")
            chunk = items[b0:b0 + per]
            hist = sr.strain_histories(rec, np.array([row for _, row in chunk], dtype=np.int64))
            gam, _ = NL.effective_strain(NL.shear_strain_history(hist, istr, etype), 1.0)
            for (key, _), gm in zip(chunk, gam):
                gmax_pct[key] = 100.0 * float(gm)
    if sr.seismic:
        direction = NL.DIRECTIONS[cm] if 0 <= int(cm) < 3 else str(cm)
    else:                                  # foundation vibration: the strains of a load case
        direction = f"LOAD{('-' + case) if case not in ('', '0') else ''}"
    try:
        f74 = NL.strain_file(f78, gmax_pct, curves, direction=direction)
    except NL.NLSoilError as exc:
        raise ModuleError(f"non-linear soil: {exc}") from None
    NL.write_nlfile(sr.wd / "FILE74", f74, title=str(sr.d["title"]))
    sr.files.append(f"FILE74  (effective strains of {len(f74.rows)} nonlinear soil elements, input {direction}, "
                    f"iteration {f74.iteration})")
    # ---- listing
    lst.section(f"Non-linear soil: effective strains of iteration {f78.iteration}, input {direction} (FILE74; "
                "requirements 4.10 item 5)")
    lst.write(f"   ESF = {f78.esf:g};  gamma_eff = ESF x max|gamma(t)| over the Fourier period ({sr.nfft} steps); "
              "strains in percent; G, beta from the FILE73 curves")
    for gi in f78.groups:
        lst.write(f"   group {gi.group} ({gi.etype}, {gi.nelem} elements): ISTR {gi.istr} = {NL.STRAIN_FLAGS[gi.istr]}")
    used = f78.by_key()
    lst.write(f"{'group':>6s}{'elem':>7s}{'g_max %':>12s}{'g_eff %':>12s}{'curve':>6s}{'G/Gmax':>9s}{'beta':>8s}"
              f"{'G used/Gmax':>13s}{'beta used':>11s}{'dG/G %':>9s}{'dbeta %':>9s}")
    for k, r in enumerate(f74.rows):
        if k >= 5000:
            lst.write(f" ... {len(f74.rows) - 5000} more elements (FILE74 lists them all)")
            break
        u = used[r.key]
        dG = (r.G - u.G) / r.G * 100.0 if r.G else 0.0
        lst.write(f"{r.group:>6d}{r.element:>7d}{r.extra:>12.5g}{r.gamma_eff:>12.5g}{r.curve:>6d}{r.ratio:>9.4f}"
                  f"{r.beta:>8.4f}{u.ratio:>13.4f}{u.beta:>11.4f}{dG:>9.3f}{(r.beta - u.beta) * 100.0:>9.3f}")
    for g, n, r0, r1, r2, b0, b1, b2 in NL.group_summary(f74.rows):
        lst.write(f"   group {g}: {n} elements, strain-compatible G/Gmax {r0:.4f} .. {r2:.4f} (mean {r1:.4f}), "
                  f"beta {b0:.4f} .. {b2:.4f} (mean {b1:.4f})")
    conv = NL.compare(f78.rows, f74.rows, f78.iteration)
    lst.write(f"   input {direction} alone: {conv.text()}")
    lst.write("   (for X/Y/Z input the decisive measure is that of the combined strains: COMBXYZSTRAIN, D-NLS-04)")


def _write_element_frame(path: Path, group_elem: np.ndarray, values: np.ndarray) -> Path:
    """Element frame (D-FIL-03 header ``nrows ncols``): rows ``group element value``."""
    lines = [f"{len(values)} 3"]
    lines += [f"{int(g)} {int(e)} {float(v):.10e}" for (g, e), v in zip(group_elem, values)]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return Path(path)


def _nyquist_real(Hg: np.ndarray, nfft: int) -> np.ndarray:
    """Nyquist bin real (D-MOT-03) without touching the vibration f = 0 value (H_1)."""
    H = np.array(Hg, dtype=complex, copy=True)
    if nfft % 2 == 0:
        H[-1] = H[-1].real
    return H


if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
