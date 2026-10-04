"""MOTION: SSI transfer functions, nodal time histories and in-structure response spectra.

What the module does (requirements section 4.8, spec 05c Part B, R1 section 5)
-----------------------------------------------------------------------------
ANALYS solved the SSI problem at the SSI frequencies only and stored, for every equation
(node, DOF), the complex transfer function (TF) ``H = U/U_cp`` in FILE8 (seismic analysis:
total motion per unit control motion, dimensionless and equal for displacement and
acceleration; vibration analysis: displacement per unit load factor).  MOTION

1. reads the control-motion acceleration history (THFILE, in g), selects records
   rec1..rec2, scales it by ``mult`` or to ``max`` (exactly one non-zero, Errors 77/78), checks
   it fits in the Fourier period and transforms it: ``A = rfft(a)``;
2. checks that the Fourier grid ``df = 1/(NFFT dt)`` reproduces the FILE8 frequency step
   (D-MOT-14);
3. interpolates the requested TFs from the SSI frequencies to every Fourier frequency
   (:func:`sassi.core.interp.interpolate_tf`, options 0-6, smoothing, phase adjustment);
4. convolves: ``r = irfft(H A)`` -- the response acceleration in g (seismic), or the
   displacement / velocity / acceleration ``H F, i w H F, -w^2 H F`` for a load history F
   (vibration, MOTIONX <resp>);
5. writes the results: ``.TFU`` (computed TF at the SSI frequencies), ``.TFI`` (interpolated TF
   on the Fourier grid 0..f_N), ``.ACC`` (history, length 1.2*dur), ``.RS`` per damping
   (Nigam-Jennings over the full record, SA in g), the listing with the maxima (ZPA), FILE13 /
   FILE12, CONTTRS ``.RSO`` conversions, SRSS TFs (SRSSTF.txt) and post-processing frames.

Incoherent input (requirements 2.4, 4.4 item 6, 4.8 items 5-6): a FILE8 that ANALYS computed with
HOUSE incoherency factors (FILE77) is described in the listing; the FILE8 of a single incoherent mode
k > 1 (HOUSE <nmodes> = -k) has the zero-frequency anchor H(0) = 0 (its factors vanish as f -> 0,
D-INC-03), also when it enters an SRSS TF (SRSSTF.txt), whose modal files are identified by the
mode ANALYS recorded.  Stochastic samples (FILE8001 ...) are processed one at a time; their response
spectra are averaged afterwards (Simulation Mean, line mathematics AVERAGE).

File names follow requirements section 1.8 (``00012TR_X.TFU``, ``00012TR_X01.RS`` ...).  Text
formats follow D-FIL-02/03; values are written with 17 significant digits so files read back
bit-faithfully (RELDISP works from the complex ``.TFI``).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .. import conventions as C
from ..core import interp as I
from ..core import signal as S
from ..core.spectra import log_frequencies, response_spectrum
from ..io import decks, textfiles, thfile
from ..io.container import Container, read_container
from ..io.files import validate
from .base import ModuleContext, ModuleError, batch_main

NAME = "MOTION"

MAX_NFFT = 65536          # D-CNV-11 (MOTION, RELDISP, STRESS)
MAX_SSI_FREQ = 1500       # requirements section 4.0.1 (MOTION input)
MAX_RS_DAMP = 5           # manual MOTION limit (D-MOT-15: more values -> warning, all computed)
DF_RTOL = 1e-6            # D-MOT-14
DT_RTOL = 1e-6            # D-MOT-12 (time step of the history file vs <delt>)
BATCH = 64                # columns processed together (bounded memory)
RESPONSE_NAMES = {0: "displacement", 1: "velocity", 2: "acceleration"}


# =======================================================================================
# Text writers (D-FIL-02 / D-FIL-03 layouts, full precision) -- also used by RELDISP
# =======================================================================================
def write_tf_file(path: Path, f: np.ndarray, H: np.ndarray, cplx: bool = True, header: str = "") -> Path:
    """Write ``.TFU/.TFI/.TFD``: optional ``#`` header, rows ``f |H| [phase_rad]`` (D-FIL-02).

    Same layout as :func:`sassi.io.textfiles.write_tf` but with 17 significant digits, so a
    complex TF survives the text round trip to ~1e-16 (RELDISP verification, VP-34).
    """
    path = Path(path)
    H = np.asarray(H, dtype=complex) + 0.0       # + 0.0 turns signed zeros (-0.0) into +0.0 (phase 0, not pi)
    amp = np.abs(H)
    lines = [f"# {header}"] if header else []
    if cplx:
        ph = np.angle(H)
        lines += [f"{fi:.10f} {a:.16e} {p:.16e}" for fi, a, p in zip(f, amp, ph)]
    else:
        lines += [f"{fi:.10f} {a:.16e}" for fi, a in zip(f, amp)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_history_file(path: Path, values: np.ndarray, dt: float, header: str = "") -> Path:
    """Write ``.ACC/.THD``: first line dt, then one value per line (D-FIL-02, THFILE fopt 0).

    The optional header is written as a ``#`` line before dt (readers skip it).
    """
    path = Path(path)
    lines = [f"# {header}"] if header else []
    lines.append(f"{float(dt):.16g}")
    lines += [f"{v:.16e}" for v in np.asarray(values, dtype=float)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_xy_file(path: Path, x: np.ndarray, ys: Sequence[np.ndarray], header: str = "") -> Path:
    """Write columns ``x y1 y2 ...`` (``.RS``/``.RSO``: f SA)."""
    path = Path(path)
    lines = [f"# {header}"] if header else []
    Y = np.column_stack([np.asarray(y, dtype=float) for y in ys]) if len(ys) else np.zeros((len(x), 0))
    for i in range(len(x)):
        lines.append(" ".join([f"{x[i]:.10f}"] + [f"{v:.16e}" for v in Y[i]]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_frame(path: Path, nodes: np.ndarray, values: np.ndarray) -> Path:
    """Frame file (D-FIL-03): header ``nrows ncols``, rows ``node v1 ... v(ncols-1)``."""
    path = Path(path)
    nodes = np.asarray(nodes)
    values = np.asarray(values, dtype=float).reshape(len(nodes), -1)
    lines = [f"{len(nodes)} {values.shape[1] + 1}"]
    for nd, row in zip(nodes, values):
        lines.append(" ".join([str(int(nd))] + [f"{v:.10e}" for v in row]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def frame_freq(f: float) -> str:
    """Frame frequency field ``fff.ff`` (requirements section 1.8)."""
    return f"{f:06.2f}"


def frame_time(t: float) -> str:
    """Frame time field ``tt.ttt`` (requirements section 1.8)."""
    return f"{t:06.3f}"


def find_file_nocase(workdir: Path, name: str) -> Optional[Path]:
    """Case-insensitive lookup of ``name`` in ``workdir`` (D-FIL-07)."""
    p = workdir / name
    if p.exists():
        return p
    low = name.lower()
    for q in workdir.iterdir():
        if q.name.lower() == low:
            return q
    return None


# =======================================================================================
# Control motion (requirements section 4.8 item 1, UT-07)
# =======================================================================================
def scale_control_motion(acc: np.ndarray, mult: float, maxval: float) -> np.ndarray:
    """Scale a control motion by ``mult`` or to the peak ``maxval`` (exactly one non-zero).

    ``a <- mult a`` (mult != 0, max = 0) or ``a <- a max/max|a|`` (max != 0, mult = 0).
    Both zero is Error 77, both non-zero Error 78 (spec 05c B.5.10).
    """
    mult = float(mult or 0.0)
    maxval = float(maxval or 0.0)
    if mult == 0.0 and maxval == 0.0:
        raise ModuleError("Error 77: multiplication factor and maximum value are both zero "
                          "(exactly one must be non-zero)")
    if mult != 0.0 and maxval != 0.0:
        raise ModuleError("Error 78: multiplication factor and maximum value are both non-zero "
                          "(use only one of them)")
    acc = np.asarray(acc, dtype=float)
    if maxval != 0.0 and (acc.size == 0 or np.max(np.abs(acc)) == 0.0):
        raise ModuleError("control motion is identically zero: it cannot be scaled to a maximum value")
    return thfile.scale_history(acc, mult=mult, maxval=maxval)


@dataclass
class ControlMotion:
    """Scaled control motion (or load history) and its Fourier transform."""
    acc: np.ndarray            # scaled record (n samples, g for seismic)
    padded: np.ndarray         # zero-padded to NFFT
    A: np.ndarray              # rfft of the padded record
    dt: float
    nrec: int                  # records in the file (before rec1/rec2 selection)
    path: Path


def read_control_motion(path: Path, fopt: int, rec1: int, rec2: int, delt: float, nfft: int,
                        mult: float, maxval: float) -> ControlMotion:
    """Read, select, scale, check and transform the control motion (requirements section 4.8 item 1).

    Records are 1-based samples after the format header (D-MOT-12).  The time step of the file
    (first value for fopt 0, the time column for fopt 1) must equal ``delt`` (relative 1e-6).
    """
    path = Path(path)
    if not path.exists():
        raise ModuleError(f"Error 73: time-history file {path} not found")
    try:
        full, dt_file = thfile.read_history(path, fopt=fopt, rec1=1, rec2=0)
    except ValueError as exc:
        raise ModuleError(f"time-history file {path.name}: {exc}") from exc
    nrec = len(full)
    r1 = int(rec1 or 1)
    r2 = int(rec2 or 0)
    if r1 < 1 or r1 > max(nrec, 1):
        raise ModuleError(f"Error 74: first record {r1} outside 1..{nrec} of {path.name}")
    if r2 < 0 or r2 > nrec:
        raise ModuleError(f"Error 75: last record {r2} outside 0..{nrec} of {path.name}")
    if r2 and r1 > r2:
        raise ModuleError(f"Error 76: first record {r1} is greater than last record {r2}")
    acc = full[r1 - 1:(r2 if r2 else nrec)].copy()
    if dt_file is not None and abs(dt_file - delt) > DT_RTOL * abs(delt):
        raise ModuleError(f"time step of {path.name} ({dt_file:g} s) differs from the MOTION time step "
                          f"{delt:g} s")
    acc = scale_control_motion(acc, mult, maxval)
    if len(acc) > nfft:
        raise ModuleError(f"EDU-03: time history has {len(acc)} records, longer than the Fourier period "
                          f"NFFT = {nfft}; use NFFT >= {1 << int(np.ceil(np.log2(len(acc))))}")
    padded = S.pad_record(acc, nfft)
    return ControlMotion(acc=acc, padded=padded, A=np.fft.rfft(padded), dt=float(delt), nrec=nrec, path=path)


# =======================================================================================
# FILE8 and output requests
# =======================================================================================
def read_file8(path: Path, name: str = "FILE8") -> Container:
    """Read and validate a FILE8-type container (ARCHITECTURE section 4)."""
    try:
        c = read_container(path, kind="FILE8")
    except (ValueError, FileNotFoundError) as exc:
        raise ModuleError(f"{name}: {exc}") from exc
    probs = validate(c)
    if probs:
        raise ModuleError(f"{name} does not satisfy the FILE8 schema: " + "; ".join(probs))
    fnum = np.asarray(c["fnum"])
    if np.any(np.diff(fnum) <= 0):
        raise ModuleError(f"{name}: frequency numbers are not strictly increasing")
    H = np.asarray(c["H"])
    if H.shape != (len(fnum), len(c["eq_node"])):
        raise ModuleError(f"{name}: H has shape {H.shape}, expected ({len(fnum)}, {len(c['eq_node'])})")
    return c


@dataclass
class Request:
    """One nodal output request (NOUT flags, spec 05c B.5.8, merged per (node, DOF))."""
    node: int
    dof: int
    tf: bool = False           # flag 1: printed plot of TF -> .TFU/.TFI
    th: bool = False           # flag 2: save time history (.ACC)
    plot: bool = False         # flag 3: plot time history (UI only)
    plot_rs: bool = False      # flag 4: deprecated
    rs: bool = False           # flag 5: save RS (.RS)
    mx: bool = False           # flag 6: print maximum
    listing_tf: bool = False   # print the TF in the listing (explicit NOUT request)
    explicit: bool = False     # requested by NOUT (not only by an all-points option)
    col: int = -1              # FILE8 column

    def merge(self, other: "Request") -> None:
        for k in ("tf", "th", "plot", "plot_rs", "rs", "mx", "listing_tf", "explicit"):
            setattr(self, k, getattr(self, k) or getattr(other, k))


def nout_requests(rows: List[Dict], listing) -> Dict[Tuple[int, int], Request]:
    """NOUT table rows -> requests keyed by (node, dof); duplicates merged with OR (D-MOT-10)."""
    out: Dict[Tuple[int, int], Request] = {}
    dups = []
    for r in rows:
        node, dof = int(r["node"]), int(r["dir"])
        if dof not in C.DOF_TAGS:
            raise ModuleError(f"NOUT direction {dof} for node {node} must be 1..6")
        flags = [bool(int(r[f"c{i}"])) for i in range(1, 7)]
        q = Request(node, dof, tf=flags[0], th=flags[1], plot=flags[2], plot_rs=flags[3], rs=flags[4],
                    mx=flags[5], listing_tf=flags[0], explicit=True)
        if (node, dof) in out:
            dups.append((node, dof))
            out[(node, dof)].merge(q)
        else:
            out[(node, dof)] = q
    if dups:
        listing.warning("EDU-16: duplicate output requests merged (flags OR'ed): "
                        + ", ".join(f"{n} {C.DOF_TAGS[d]}" for n, d in dups[:20])
                        + (" ..." if len(dups) > 20 else ""))
    return out


# =======================================================================================
# Transfer-function sources: plain FILE8 or SRSS of modal FILE8s
# =======================================================================================
class _TFSource:
    """Computed and interpolated TFs of FILE8 columns (requirements section 4.8 items 2-6)."""

    def __init__(self, f_ssi: np.ndarray, H: np.ndarray, eq_dof: np.ndarray, mode: str, option: int,
                 smooth: float, pzadj: int, cm: int, ang: float):
        self.f_ssi = f_ssi
        self.H = H
        self.eq_dof = eq_dof
        self.mode = mode
        self.option = option
        self.smooth = smooth
        self.pzadj = pzadj
        self.cm = cm
        self.ang = ang
        #: scale of the zero-frequency anchor: 0 for the FILE8 of a single incoherent mode k > 1
        #: (its rigid-body limit vanishes, D-INC-03; :func:`incoherent_anchor_scale`)
        self.anchor_scale = 1.0

    def anchors(self, cols: np.ndarray) -> np.ndarray:
        return self.anchor_scale * np.array([I.rigid_body_anchor(int(self.eq_dof[c]), self.cm, self.ang)
                                             for c in cols], dtype=complex)

    def computed(self, cols: np.ndarray) -> np.ndarray:
        return self.H[:, cols]

    def grid(self, cols: np.ndarray, f_grid: np.ndarray) -> np.ndarray:
        return I.interpolate_tf(self.f_ssi, self.H[:, cols], f_grid, self.option, smooth=self.smooth,
                                pzadj=self.pzadj, h0=self.anchors(cols), mode=self.mode)


class _SRSSSource(_TFSource):
    """SRSS TF of modal FILE8s (D-MOT-06, spec 05c B.5.4 / B.6.9).

    ``|H| = sqrt(sum_k |H_k|^2)``; phase 0 (p = 0) or the coherent phase (p = 1).  Each modal TF
    is interpolated first, then combined (default order of D-MOT-06).  Zero-frequency anchors:
    the rigid-body value for the coherent TF and for mode 1, 0 for the higher modes -- as f -> 0
    the coherency matrix tends to all ones, its first eigenvector carries the coherent motion
    and the other eigenvalues vanish (D-INC-03).  The mode of every modal file is the incoherent
    mode recorded by ANALYS in its FILE8 (``x_mode``, from the HOUSE FILE77 of ``<nmodes> = -k``);
    for files without that record the first listed file is taken as mode 1.
    """

    def __init__(self, modal: List[np.ndarray], coherent: Optional[np.ndarray],
                 modal_modes: Optional[List[int]] = None, **kw):
        super().__init__(H=modal[0], **kw)
        self.modal = modal
        self.coherent = coherent
        self.modal_modes = list(modal_modes) if modal_modes else [0] * len(modal)

    def _mode_one(self, k: int) -> bool:
        mk = self.modal_modes[k] if k < len(self.modal_modes) else 0
        return mk == 1 or (mk == 0 and k == 0)

    def computed(self, cols: np.ndarray) -> np.ndarray:
        amp = np.sqrt(sum(np.abs(Hk[:, cols]) ** 2 for Hk in self.modal))
        if self.coherent is None:
            return amp.astype(complex)
        return amp * np.exp(1j * np.angle(self.coherent[:, cols]))

    def grid(self, cols: np.ndarray, f_grid: np.ndarray) -> np.ndarray:
        h0 = self.anchors(cols)
        sq = 0.0
        for k, Hk in enumerate(self.modal):
            Hi = I.interpolate_tf(self.f_ssi, Hk[:, cols], f_grid, self.option, smooth=self.smooth, pzadj=0,
                                  h0=h0 if self._mode_one(k) else np.zeros_like(h0), mode=self.mode)
            sq = sq + np.abs(Hi) ** 2
        amp = np.sqrt(sq)
        if self.coherent is None:
            out = amp.astype(complex)
        else:
            Hc = I.interpolate_tf(self.f_ssi, self.coherent[:, cols], f_grid, self.option, smooth=self.smooth,
                                  pzadj=0, h0=h0, mode=self.mode)
            out = amp * np.exp(1j * np.angle(Hc))
        if self.pzadj:
            eff = I.effective_option(self.option, len(self.f_ssi))
            out = I.phase_adjust(out, I.phase_factor(eff, self.smooth))
        return out


def incoherency_note(meta: Dict) -> str:
    """Description of the incoherency of a FILE8 written by ANALYS from HOUSE FILE77 factors
    (requirements 4.6 item 5): '' for a coherent FILE8."""
    if not meta.get("file77"):
        return ""
    method = str(meta.get("x_incoherency", ""))
    what = {"AS": "incoherent, deterministic AS sum of the modes",
            "SINGLE": f"incoherent, single mode {int(meta.get('x_mode', 0) or 0)}",
            "SS": f"incoherent, stochastic sample {int(meta.get('sim', 0) or 0)}",
            "COHERENT": "coherent motion with wave passage / multiple excitation",
            "BUILT": "incoherent, per-level factors (BUILDFILE77)"}.get(method, method or "FILE77 factors")
    return (f"{what}; {'FFM' if int(meta.get('ffm', 0) or 0) else 'FFL'}; factors {meta.get('file77')} "
            f"direction {meta.get('x_file77_direction', '?')}")


def incoherent_anchor_scale(meta: Dict) -> float:
    """Scale of the zero-frequency anchor H(0) of a FILE8 (D-MOT-03, D-INC-03): 0 for a single
    incoherent mode k > 1 (``sqrt(lambda_k) phi_k -> 0`` as f -> 0), 1 otherwise."""
    if str(meta.get("x_incoherency", "")) == "SINGLE" and int(meta.get("x_mode", 0) or 0) > 1:
        return 0.0
    return 1.0


def _check_srss_modes(listing, names: List[str], modal: List[Container], coh: Optional[Container]) -> List[int]:
    """Incoherent mode of every modal FILE8 of SRSSTF.txt (ANALYS meta ``x_mode`` of a single-mode
    run, HOUSE <nmodes> = -k; 0 = unknown) with consistency warnings (SRSS TF, spec 05b 3.18)."""
    modes = []
    for nm, c in zip(names, modal):
        single = str(c.meta.get("x_incoherency", "")) == "SINGLE"
        modes.append(int(c.meta.get("x_mode", 0) or 0) if single else 0)
        if not single:
            listing.warning(f"SRSSTF.txt: {nm} is not the FILE8 of a single incoherent mode (HOUSE <nmodes> = -k, "
                            "ANALYS with FILE77): its zero-frequency anchor follows its position in the list")
    known = [m for m in modes if m]
    dup = sorted({m for m in known if known.count(m) > 1})
    if dup:
        listing.warning(f"SRSSTF.txt: incoherent modes {dup} are listed more than once")
    if known:
        listing.write("   incoherent modes: " + ", ".join(str(m) if m else "?" for m in modes))
    if coh is not None and coh.meta.get("file77") and str(coh.meta.get("x_incoherency", "")) != "COHERENT":
        listing.warning("SRSSTF.txt: the coherent-phase file (first name, phase option 1) is an incoherent FILE8")
    return modes


def read_srsstf(workdir: Path) -> Tuple[int, int, List[str]]:
    """Parse ``SRSSTF.txt``: ``n p`` then FILE8 names, the coherent one first when p = 1 (D-FIL-07)."""
    p = find_file_nocase(workdir, "SRSSTF.txt")
    if p is None:
        raise ModuleError("SRSS TF requested (MOTIONX srss = 1) but SRSSTF.txt is missing")
    lines = [ln.strip() for ln in p.read_text(encoding="utf-8", errors="replace").splitlines()]
    lines = [ln for ln in lines if ln and ln[0] not in "#*"]
    if not lines:
        raise ModuleError("SRSSTF.txt is empty")
    head = lines[0].replace(",", " ").split()
    try:
        n, ph = int(head[0]), int(head[1]) if len(head) > 1 else 0
    except ValueError as exc:
        raise ModuleError(f"SRSSTF.txt: first line must be '<# of modes> <phase option>': {lines[0]!r}") from exc
    if ph not in (0, 1) or n < 1:
        raise ModuleError("SRSSTF.txt: number of modes must be >= 1 and the phase option 0 or 1")
    names = lines[1:]
    need = n + (1 if ph == 1 else 0)
    if len(names) < need:
        raise ModuleError(f"SRSSTF.txt lists {len(names)} files, {need} needed ({n} modes, phase option {ph})")
    return n, ph, names[:need]


# =======================================================================================
# Listing helpers
# =======================================================================================
def _printer_plot(listing, f: np.ndarray, H: np.ndarray, width: int = 40) -> None:
    """Printed plot of a TF (manual "Printed Plot of Transfer Functions"): table + bar."""
    amp = np.abs(H)
    top = float(amp.max()) if amp.size and amp.max() > 0 else 1.0
    listing.write(f"{'freq (Hz)':>12s} {'|H|':>14s} {'phase (deg)':>12s}")
    for fi, a, h in zip(f, amp, H):
        bar = "*" * int(round(width * a / top))
        listing.write(f"{fi:12.4f} {a:14.6e} {np.degrees(np.angle(h)):12.3f} |{bar}")


def _rs_listing(listing, node: int, tag: str, rs_f: np.ndarray, damps: np.ndarray, rs: Dict[str, np.ndarray],
                unit: str, sv_scale: float, sv_unit: str, max_rows: int = 31) -> None:
    """Acceleration (SA, absolute) and velocity (SV, relative) spectra in the listing (D-MOT-07).

    ``sv_scale`` converts SV from (acceleration unit x s) to velocity units (gravity for SA in g).
    """
    listing.section(f"Response spectra node {node} {tag} (SA absolute acceleration in {unit}, "
                    f"SV relative velocity in {sv_unit})")
    step = max(1, int(np.ceil(len(rs_f) / max_rows)))
    rows = sorted(set(range(0, len(rs_f), step)) | {len(rs_f) - 1})
    for idd, z in enumerate(damps):
        sa, sv = rs["SA"][idd], rs["SV"][idd] * sv_scale
        listing.write(f"   damping {z:g}: peak SA {sa.max():.6e} at {rs_f[np.argmax(sa)]:.4g} Hz, peak SV "
                      f"{sv.max():.6e} at {rs_f[np.argmax(sv)]:.4g} Hz, SA at {rs_f[-1]:.4g} Hz {sa[-1]:.6e}")
        listing.write(f"{'freq (Hz)':>12s} {'SA':>14s} {'SV':>14s}")
        for k in rows:
            listing.write(f"{rs_f[k]:12.4f} {sa[k]:14.6e} {sv[k]:14.6e}")


# =======================================================================================
# Module run
# =======================================================================================
@dataclass
class _Plan:
    """Per-column work list."""
    cols: np.ndarray
    node: np.ndarray
    dof: np.ndarray
    tf: np.ndarray
    th: np.ndarray
    rs: np.ndarray
    mx: np.ndarray
    listing_tf: np.ndarray
    frame: np.ndarray          # needed for ACC / RS frames
    th_print: np.ndarray
    explicit: np.ndarray       # NOUT request (RS summary in the listing)


def _validate_deck(d, listing) -> None:
    """Runtime checks of the MOTION deck (spec 05c B.8; CHECK repeats them before AFWRITE)."""
    if d["interp"] not in I.OPTION_NAMES:
        raise ModuleError(f"interpolation option {d['interp']} must be 0..6")
    if d["pzadj"] not in (0, 1):
        raise ModuleError(f"phase adjustment {d['pzadj']} must be 0 or 1")
    if d["smo"] < 0:
        raise ModuleError("smoothing parameter must be >= 0")
    if d["step"] < 0:
        raise ModuleError("Error 67: negative output time-history step")
    if d["dur"] < 0:
        raise ModuleError("Error 68: negative total duration")
    if d["freq1"] < 0:
        raise ModuleError("Error 69: negative first response-spectrum frequency")
    if d["freq2"] < 0:
        raise ModuleError("Error 70: negative last response-spectrum frequency")
    if d["fstep"] < 0:
        raise ModuleError("Error 71: negative number of response-spectrum frequency steps")
    if d["type"] not in (0, 1):
        raise ModuleError(f"type of analysis {d['type']} must be 0 (seismic) or 1 (vibration)")
    if d["resp"] not in (0, 1, 2):
        raise ModuleError(f"MOTIONX <resp> {d['resp']} must be 0, 1 or 2")
    if d["bl"] not in (0, 1):
        raise ModuleError(f"baseline correction <bl> {d['bl']} must be 0 or 1 (D-MOT-08)")
    if d["f1213"] not in (0, 1, 2):
        raise ModuleError(f"MOTIONX <f1213> {d['f1213']} must be 0, 1 or 2")
    nft = int(d["nft"])
    if not C.is_power_of_two(nft):
        raise ModuleError(f"NFFT = {nft} is not a power of 2 (AFWRITE rounds it, Warning 9)")
    if nft > MAX_NFFT:
        raise ModuleError(f"NFFT = {nft} exceeds the MOTION maximum {MAX_NFFT} (D-CNV-11)")
    if d["delt"] <= 0:
        raise ModuleError("time step <delt> must be > 0")
    if d["smo"] > 0 and d["interp"] == 6:
        listing.warning("EDU-14: smoothing has no effect with interpolation option 6 (spline); ignored")
    if d["pzadj"] == 1 and d["interp"] != 6 and 0 < d["smo"] < 10:
        listing.warning("phase adjustment with 0 < S < 10 is weak (rho = 1/(1+S))")
    if d["pzadj"] == 1 and d["interp"] != 6 and d["smo"] == 0:
        listing.warning("phase adjustment with S = 0 leaves the phase unchanged (rho = 1); it needs S > 0")


def rs_requested(d, reqs: Dict[Tuple[int, int], Request]) -> bool:
    """True when the run computes response spectra: NOUT flag 5, Save RS in All Points, RS frames
    (full output only, ``<out>`` = 0) or the CONTTRS conversion ``<cnvrt>`` = 1."""
    full = not int(d["out"])
    nodal = any(q.rs for q in reqs.values()) or bool(int(d["savers"])) or bool(int(d["rstrs"]))
    return bool(int(d["cnvrt"])) or (full and nodal)


def _rs_settings(d, listing, required: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Damping list and log-spaced RS frequencies (requirements section 4.8 item 8, D-MOT-07).

    When response spectra are ``required`` the frequency data must satisfy the implementation
    decision of spec 05c B.8: ``0 < freq1`` (Error 69), ``freq1 < freq2`` (Error 70) and
    ``fstep >= 2`` (Error 71); otherwise the RS grid would be empty or a single frequency.
    """
    damps = np.asarray([float(r["value"]) for r in d.rows("damp")], dtype=float)
    if np.any((damps <= 0) | (damps >= 1)):
        raise ModuleError("Error 72: response-spectrum damping must be in (0, 1)")
    if len(damps) > MAX_RS_DAMP:
        listing.warning(f"{len(damps)} damping values: the manual MOTION limit is {MAX_RS_DAMP}; all are "
                        "computed (D-MOT-15)")
    f1, f2, n = float(d["freq1"]), float(d["freq2"]), int(d["fstep"])
    if required:
        if f1 <= 0:
            raise ModuleError(f"Error 69: first response-spectrum frequency {f1:g} Hz must be > 0")
        if f2 <= f1:
            raise ModuleError(f"Error 70: last response-spectrum frequency {f2:g} Hz must exceed the first "
                              f"({f1:g} Hz)")
        if n < 2:
            raise ModuleError(f"Error 71: number of response-spectrum frequency steps {n} must be >= 2")
    return damps, (np.zeros(0) if n < 2 or f1 <= 0 or f2 <= f1 else log_frequencies(f1, f2, n))


def _convert_external(ctx: ModuleContext, d, damps: np.ndarray, rs_f: np.ndarray, nfft: int) -> None:
    """CONTTRS.txt: response spectra of external ``.ACC`` histories -> ``.RSO`` (spec 05c B.5.11)."""
    lst = ctx.listing
    lst.section("Convert time histories to response spectra (CONTTRS.txt)")
    p = find_file_nocase(ctx.workdir, "CONTTRS.txt")
    if p is None:
        lst.warning("<cnvrt> = 1 but CONTTRS.txt is missing; no conversion")
        return
    if len(rs_f) == 0 or len(damps) == 0:
        lst.warning("no response-spectrum frequencies or damping values; no conversion")
        return
    names = [ln.strip() for ln in p.read_text(encoding="utf-8", errors="replace").splitlines()]
    names = [n for n in names if n and n[0] not in "#*"]
    dt = float(d["delt"])
    for nm in names:
        src = Path(nm) if Path(nm).is_absolute() else ctx.workdir / nm
        if src.suffix.upper() != ".ACC":
            lst.warning(f"{nm}: external histories must have the .ACC extension; skipped")
            continue
        if not src.exists():
            lst.warning(f"{nm}: file not found; skipped")
            continue
        a, dtf = textfiles.read_history(src)
        if abs(dtf - dt) > DT_RTOL * dt:
            lst.warning(f"{nm}: time step {dtf:g} s differs from {dt:g} s; skipped")
            continue
        if len(a) > nfft:
            lst.warning(f"{nm}: {len(a)} values exceed NFFT = {nfft}; skipped")
            continue
        rs = response_spectrum(S.pad_record(a, nfft), dt, rs_f, damps)
        out = src.with_suffix(".RSO")
        write_xy_file(out, rs_f, [rs["SA"][i] for i in range(len(damps))],
                      header="f(Hz) SA for damping " + " ".join(f"{z:g}" for z in damps))
        lst.write(f"   {nm}: {len(a)} values, PGA {np.max(np.abs(a)):.6g}, ZPA of RS {rs['SA'][0, -1]:.6g}"
                  f" -> {out.name}")


def run(ctx: ModuleContext) -> int:
    """MOTION module entry point (requirements section 4.8)."""
    lst = ctx.listing
    d = decks.read(ctx.deck_path, NAME)
    if d["title"]:
        lst.write(f" {d['title']}")
    _validate_deck(d, lst)
    seismic = int(d["type"]) == 0
    mode = "seismic" if seismic else "vibration"
    nfft, dt = int(d["nft"]), float(d["delt"])
    df = 1.0 / (nfft * dt)
    option, smooth, pzadj = int(d["interp"]), float(d["smo"]), int(d["pzadj"])
    gravity = float(d["gravity"])

    lst.section("MOTION options")
    lst.write(f"   operation mode            : {'data check' if d['opmode'] else 'complete solution'}")
    lst.write(f"   type of analysis          : {mode}")
    lst.write(f"   output                    : {'transfer functions only' if d['out'] else 'full'}")
    lst.write(f"   interpolation option      : {option} ({I.OPTION_NAMES[option]})")
    lst.write(f"   smoothing S, phase adj.   : {smooth:g}, {pzadj}")
    lst.write(f"   NFFT, dt, df              : {nfft}, {dt:g} s, {df:.9g} Hz (Fourier period {nfft * dt:g} s)")
    lst.write(f"   control motion            : {d['thfile']}  (fopt {d['fopt']}, records {d['rec1']}..{d['rec2'] or 'last'},"
              f" mult {d['mult']:g}, max {d['max']:g})")
    lst.write(f"   baseline correction       : {'Hudson-Housner' if d['bl'] else 'none'}")
    if not seismic:
        lst.write(f"   vibration response        : {RESPONSE_NAMES[int(d['resp'])]} (RS from acceleration)")
    reqs = nout_requests(d.rows("nout"), lst)
    damps, rs_f = _rs_settings(d, lst, required=rs_requested(d, reqs))
    if len(damps):
        lst.write(f"   RS dampings               : {', '.join(f'{z:g}' for z in damps)}")
        lst.write(f"   RS frequencies            : {d['freq1']:g} .. {d['freq2']:g} Hz, {int(d['fstep'])} points (log)")

    for q in reqs.values():
        if q.plot_rs:
            lst.write(" NOTE: NOUT flag 4 (plot acceleration and velocity RS) is deprecated and ignored")
            break
    allpts = any(int(d[k]) for k in ("savetf", "saveacc", "savers", "saverot", "rsttf", "rstacc", "rstrs"))
    if not reqs and not allpts and not int(d["cnvrt"]):
        raise ModuleError("Error 64: no nodal output request")

    # ---------------------------------------------------------------- control motion
    thpath = Path(d["thfile"]) if d["thfile"] else None
    if thpath is not None and not thpath.is_absolute():
        thpath = ctx.workdir / thpath
    need_motion = bool(reqs or allpts) and not int(d["out"])
    cmot: Optional[ControlMotion] = None
    if need_motion or (int(d["opmode"]) == 1 and thpath is not None):
        if thpath is None:
            raise ModuleError("Error 73: no time-history file (THFILE) given")
        cmot = read_control_motion(thpath, int(d["fopt"]), int(d["rec1"]), int(d["rec2"]), dt, nfft,
                                   float(d["mult"]), float(d["max"]))
        lst.section("Control motion" if seismic else "Reference load history")
        lst.write(f"   file {cmot.path.name}: {cmot.nrec} records, {len(cmot.acc)} used, "
                  f"duration {len(cmot.acc) * dt:g} s, quiet zone {(nfft - len(cmot.acc)) * dt:g} s")
        unit = "g" if seismic else "load factor"
        imax = int(np.argmax(np.abs(cmot.acc))) if len(cmot.acc) else 0
        lst.write(f"   peak value after scaling {np.max(np.abs(cmot.acc)) if len(cmot.acc) else 0:.6g} {unit} "
                  f"at t = {imax * dt:g} s")

    if int(d["opmode"]) == 1:
        lst.section("Data check")
        if reqs or allpts:
            f8 = ctx.require(d["file8"] or "FILE8", "ANALYS (or COMBIN)")
            read_file8(f8, f8.name)
        lst.write("   input checked; data-check mode writes no output files")
        return 0

    if int(d["cnvrt"]):
        _convert_external(ctx, d, damps, rs_f, nfft)
    if not reqs and not allpts:
        return 0

    # ---------------------------------------------------------------- FILE8 / SRSS
    if int(d["srss"]):
        nmod, ph, names = read_srsstf(ctx.workdir)
        conts = []
        for nm in names:
            p = find_file_nocase(ctx.workdir, nm)
            if p is None:
                raise ModuleError(f"SRSSTF.txt: modal file {nm} not found")
            conts.append(read_file8(p, nm))
        ref = conts[0]
        for c, nm in zip(conts[1:], names[1:]):
            if (not np.array_equal(c["fnum"], ref["fnum"]) or not np.array_equal(c["eq_node"], ref["eq_node"])
                    or not np.array_equal(c["eq_dof"], ref["eq_dof"])):
                raise ModuleError(f"SRSSTF.txt: {nm} has a different frequency list or DOF map than {names[0]}")
        coh = conts[0] if ph == 1 else None
        modal = conts[1:] if ph == 1 else conts
        f8 = ref
        f8name = names[0]
        lst.section("SRSS transfer functions (SRSSTF.txt)")
        lst.write(f"   {nmod} modes, phase option {ph} ({'coherent phase' if ph else 'zero phase'})")
        for nm, c in zip(names, conts):
            note = incoherency_note(c.meta)
            lst.write(f"   {nm}" + (f"   ({note})" if note else "   (coherent)"))
        if pzadj:
            lst.warning("phase adjustment combined with the SRSS TF (should be 0 with SRSS)")
        modal_modes = _check_srss_modes(lst, names[1:] if ph == 1 else names, modal, coh)
    else:
        f8name = d["file8"] or "FILE8"
        f8 = read_file8(ctx.require(f8name, "ANALYS (or COMBIN)"), f8name)
    meta = f8.meta
    fnum = np.asarray(f8["fnum"], dtype=np.int64)
    eq_node = np.asarray(f8["eq_node"], dtype=np.int64)
    eq_dof = np.asarray(f8["eq_dof"], dtype=np.int64)
    df8 = float(meta["df"])
    if abs(df - df8) > DF_RTOL * df8:
        raise ModuleError(f"frequency step of MOTION 1/(NFFT dt) = {df:.9g} Hz differs from the {f8name} step "
                          f"{df8:.9g} Hz (D-MOT-14): use the SITE NFFT and time step")
    if int(meta["type"]) != int(d["type"]):
        raise ModuleError(f"type of analysis {d['type']} does not match {f8name} (type {meta['type']})")
    # frequency numbers identify the SSI frequencies (D-ANL-08); placing them on the MOTION grid (df agrees
    # with the FILE8 step to 1e-6) makes every computed frequency coincide with a Fourier bin
    f_ssi = fnum * df
    freq8 = np.asarray(f8["freq"], dtype=float)
    if np.max(np.abs(freq8 - fnum * df8)) > 1e-6 * max(df8, 1e-12) * max(1, fnum[-1]):
        lst.warning(f"{f8name}: stored frequencies differ from fnum*df; fnum*df is used")
    if len(fnum) > MAX_SSI_FREQ:
        lst.warning(f"EDU-02: {len(fnum)} SSI frequencies exceed the MOTION input limit {MAX_SSI_FREQ}")
    if fnum[-1] > nfft // 2:
        lst.warning(f"G-02: frequency numbers above NFFT/2 = {nfft // 2}; interpolation stops at the Nyquist "
                    "frequency (D-CNV-12)")
    cm = int(meta.get("cm", d["cm"]))
    ang = float(meta.get("ang", d["ang"]))
    if cm != int(d["cm"]) or abs(ang - float(d["ang"])) > 1e-9:
        lst.warning(f"control direction/angle of {f8name} (cm {cm}, ang {ang:g}) differ from the deck "
                    f"(cm {d['cm']}, ang {d['ang']:g}); the {f8name} values are used")
    lst.section(f"Transfer functions ({f8name})")
    lst.write(f"   {len(fnum)} SSI frequencies {f_ssi[0]:.6g} .. {f_ssi[-1]:.6g} Hz, {len(eq_node)} equations, "
              f"type {mode}, control direction {('x', 'y', 'z')[cm]}' at angle {ang:g} deg")
    inc_note = "" if int(d["srss"]) else incoherency_note(meta)
    if inc_note:
        lst.write(f"   incoherent input: {inc_note}")
        if incoherent_anchor_scale(meta) == 0.0:
            lst.write("   zero-frequency anchor H(0) = 0: the factors of a single incoherent mode k > 1 vanish as "
                      "f -> 0 (D-INC-03)")
        if str(meta.get("x_incoherency", "")) == "SS" and not pzadj:
            lst.write("   stochastic sample: average the response spectra of all samples (AVERAGE) for the "
                      "Simulation Mean; phase adjustment (<pzadj> = 1) gives the EPRI approximate upper bound")
    eff = I.effective_option(option, len(fnum))
    if eff != option:
        lst.warning(f"{len(fnum)} SSI frequencies: interpolation option {option} replaced by "
                    f"{'complex linear' if eff < 0 else 'option 6'}"
                    + ("; smoothing is not applied" if option != 6 and smooth > 0 else ""))
        if pzadj and I.phase_factor(eff, smooth) != I.phase_factor(option, smooth):
            lst.warning(f"phase adjustment follows the scheme actually used: rho = "
                        f"{I.phase_factor(eff, smooth):g} instead of {I.phase_factor(option, smooth):g} "
                        "(D-MOT-05; S has no effect on the spline / linear schemes)")

    kw = dict(f_ssi=f_ssi, eq_dof=eq_dof, mode=mode, option=option, smooth=0.0 if option == 6 else smooth,
              pzadj=pzadj, cm=cm, ang=ang)
    if int(d["srss"]):
        src: _TFSource = _SRSSSource([np.asarray(c["H"]) for c in modal],
                                     None if coh is None else np.asarray(coh["H"]), modal_modes=modal_modes, **kw)
    else:
        src = _TFSource(H=np.asarray(f8["H"]), **kw)
        src.anchor_scale = incoherent_anchor_scale(meta)

    optimizer_numbering_note(ctx, [n for n, _ in reqs], lst)
    plan = _make_plan(d, reqs, eq_node, eq_dof, lst)
    if plan is None:
        lst.warning("no requested node/DOF exists in the TF file; nothing to do")
        return 0
    maxnode = max(int(d["maxnode"]), int(eq_node.max()) if eq_node.size else 0)
    _process(ctx, d, src, plan, cmot, f_ssi, damps, rs_f, maxnode, seismic, gravity)
    return 0


def optimizer_numbering_note(ctx: ModuleContext, nodes, lst, what: str = "output requests (NOUT)") -> None:
    """Requirements 2.6 "HOUSE optimizer used": when HOUSE renumbered the model (``<model>.map`` present;
    HOUSE renames a stale map ``.prev``), FILE8 and the MOTION/RELDISP files use the NEW node numbers and
    the manual has the user select the post-processing nodes from the ``.hounew``.  A request given in the
    model numbering then silently names another node, so every requested node whose model number differs
    is listed (``new = model``), and requests absent from the ``.map`` targets are listed too (final audit:
    example 1 with Optimize Model reported mat nodes 37, 78-81 for the requested 41, 82-85 without a word)."""
    p = ctx.path(f"{ctx.model}.map")
    if not p.exists():
        return
    try:
        from ..core.renumber import read_map
        new_old = {int(n): int(o) for o, n in read_map(p).items()}
    except (OSError, ValueError) as exc:
        lst.warning(f"HOUSE node optimizer map {p.name} unreadable ({exc}): node requests are read as new numbers")
        return
    req = sorted({int(n) for n in nodes})
    moved = [(n, new_old[n]) for n in req if n in new_old and new_old[n] != n]
    absent = [n for n in req if n not in new_old]
    if moved or absent:
        lst.warning(f"HOUSE node optimizer ({p.name}): FILE8 uses the new node numbers of {ctx.model}.hounew and the "
                    f"{what} are read as NEW numbers (manual: select the post-processing nodes from the .hounew)"
                    + ("; requested new node = model node: " + ", ".join(f"{n} = {o}" for n, o in moved[:20])
                       + (" ..." if len(moved) > 20 else "") if moved else "")
                    + (f"; not new node numbers of the .map: {absent[:20]}" if absent else ""))


def _make_plan(d, reqs: Dict[Tuple[int, int], Request], eq_node, eq_dof, lst) -> Optional[_Plan]:
    """Resolve requests (NOUT + all-points options) to FILE8 columns."""
    colmap = {(int(n), int(k)): i for i, (n, k) in enumerate(zip(eq_node, eq_dof))}
    nodes_present = set(int(n) for n in eq_node)
    work: Dict[int, Request] = {}
    absent, fixed = [], []
    for (node, dof), q in reqs.items():
        c = colmap.get((node, dof))
        if c is None:
            (fixed if node in nodes_present else absent).append((node, dof))
            continue
        q.col = c
        work[c] = q
    if fixed:
        lst.warning("output requests on constrained (fixed) DOFs ignored: "
                    + ", ".join(f"{n} {C.DOF_TAGS[k]}" for n, k in fixed[:20]) + (" ..." if len(fixed) > 20 else ""))
    if absent:
        lst.warning("Error 65: nodes not in the TF file (check node numbers / HOUSE optimizer .map): "
                    + ", ".join(f"{n} {C.DOF_TAGS[k]}" for n, k in absent[:20]) + (" ..." if len(absent) > 20 else ""))
    trans = np.where(eq_dof <= 3)[0]
    rot = np.where(eq_dof >= 4)[0]
    full = not int(d["out"])
    frames = bool(int(d["rstacc"]) or int(d["rstrs"])) and full

    def add(cols, **flags):
        for c in cols:
            q = work.get(int(c))
            if q is None:
                q = Request(int(eq_node[c]), int(eq_dof[c]), col=int(c))
                work[int(c)] = q
            for k, v in flags.items():
                setattr(q, k, getattr(q, k) or v)

    if int(d["savetf"]):
        add(trans, tf=True)
    if int(d["saveacc"]):
        add(trans, th=True)
    if int(d["savers"]):
        add(trans, rs=True)
    if int(d["saverot"]):
        add(rot, th=True)
    if int(d["rsttf"]):
        add(trans)                     # TFU frames hold every translational DOF
    frame_cols = set()
    if frames:
        frame_cols = set(int(c) for c in trans) | (set(int(c) for c in rot) if int(d["rstacc"]) else set())
        add(sorted(frame_cols))
    if not work:
        return None
    if not full:
        # B.5.3 "Output Only Transfer Functions": extract the TFs (TFU/TFI) at the selected nodes,
        # whatever NOUT flags (or Save ... in All Points options) selected them
        promoted = [q for q in work.values() if (q.explicit or q.th or q.rs or q.mx) and not q.tf]
        for q in promoted:
            q.tf = True
        if promoted:
            lst.write(f" NOTE: <out> = 1 (transfer functions only): .TFU/.TFI written for {len(promoted)} "
                      "requested node/DOF(s) without NOUT flag 1")
    cols = np.array(sorted(work))
    get = lambda k: np.array([bool(getattr(work[c], k)) for c in cols])
    th = get("th") & full
    return _Plan(cols=cols, node=eq_node[cols], dof=eq_dof[cols], tf=get("tf"), th=th, rs=get("rs") & full,
                 mx=get("mx") & full, listing_tf=get("listing_tf"),
                 frame=np.array([c in frame_cols for c in cols]),
                 th_print=th & (int(d["step"]) > 0), explicit=get("explicit"))


def _process(ctx: ModuleContext, d, src: _TFSource, plan: _Plan, cmot: Optional[ControlMotion],
             f_ssi: np.ndarray, damps: np.ndarray, rs_f: np.ndarray, maxnode: int, seismic: bool,
             gravity: float) -> None:
    """Interpolate, convolve and write every planned column (in batches of :data:`BATCH`)."""
    lst = ctx.listing
    wd = ctx.workdir
    nfft, dt = int(d["nft"]), float(d["delt"])
    f_grid = S.fourier_grid(nfft, dt)
    fN = min(f_ssi[-1], f_grid[-1])
    kmax = int(np.searchsorted(f_grid, fN * (1 + 1e-12), side="right"))   # TFI rows 0..f_N
    full = not int(d["out"])
    nout = S.output_length(nfft, dt, float(d["dur"]))
    resp = int(d["resp"])
    cplx = bool(int(d["cplx"]))
    bl = bool(int(d["bl"]))
    f1213 = int(d["f1213"])
    if f1213 == 1 and not bl:
        lst.warning("FILE13 requested (f1213 = 1) but baseline correction is off: FILE13 is not written")
    if bl and not seismic and resp != 2:
        lst.warning("baseline correction applies to acceleration histories only; ignored for <resp> 0/1")
    if full and cmot is None:
        raise ModuleError("Error 73: time histories requested but no control motion")
    need_rs = full and (plan.rs.any() or (int(d["rstrs"]) and plan.frame.any()))
    if need_rs and (len(rs_f) == 0 or len(damps) == 0):
        lst.warning("response spectra requested but no damping values (DAMP list empty); RS skipped")
        need_rs = False
    unit = "g" if seismic else "model units"
    sv_scale, sv_unit = (gravity, "length/s") if seismic else (1.0, "model units")
    rname = "acceleration" if seismic else RESPONSE_NAMES[resp]
    w = 2.0 * np.pi * f_grid
    lst.section("Output")
    lst.write(f"   {len(plan.cols)} node/DOF columns; TFI rows 0..{fN:.6g} Hz ({kmax} Fourier frequencies)")
    if full:
        lst.write(f"   history length {nout} samples ({nout * dt:g} s); requested response: {rname}")

    maxima = []
    fh12 = open(wd / "FILE12", "w", encoding="utf-8") if (full and f1213 == 2) else None
    fh13 = open(wd / "FILE13", "w", encoding="utf-8") if (full and f1213 == 1 and bl) else None
    # frame accumulators
    fr_acc: Dict[int, np.ndarray] = {}
    fr_rs: Dict[int, np.ndarray] = {}
    baseline_report = []
    edu18: List[Tuple[int, str, float]] = []
    try:
        ncol = len(plan.cols)
        for b0 in range(0, ncol, BATCH):
            if ctx.cancelled():
                raise ModuleError("run cancelled")
            sl = slice(b0, min(ncol, b0 + BATCH))
            cols = plan.cols[sl]
            need = (plan.th[sl] | plan.mx[sl] | plan.rs[sl] | plan.frame[sl]) if full else np.zeros(len(cols), bool)
            if not (plan.tf[sl].any() or plan.listing_tf[sl].any() or need.any()):
                continue                                               # e.g. columns kept only for TFU frames
            Hc = src.computed(cols)                                    # (nF, mb)
            if seismic:                                                # EDU-18 / G-19 low-frequency check
                h0 = src.anchors(cols)
                big = (np.abs(h0) >= 0.1) & (plan.dof[sl] <= 3)
                dev = np.abs(np.abs(Hc[0]) - np.abs(h0)) > 0.05 * np.abs(h0)
                edu18 += [(int(plan.node[b0 + j]), C.DOF_TAGS[int(plan.dof[b0 + j])], float(abs(Hc[0, j])))
                          for j in np.where(big & dev)[0]]
            Hg = src.grid(cols, f_grid)                                # (nK, mb)
            Hg = S.hermitian_bins(Hg, nfft) if seismic else _nyquist_real(Hg, nfft)
            for jj, c in enumerate(cols):
                i = b0 + jj
                node, dof = int(plan.node[i]), int(plan.dof[i])
                tag = C.DOF_TAGS[dof]
                if plan.tf[i]:
                    hdr = f"node {node} dof {tag} analysis={'seismic' if seismic else 'vibration'} complex={int(cplx)}"
                    write_tf_file(wd / C.nodal_result_name(node, dof, "TFU", maxnode), f_ssi, Hc[:, jj], cplx, hdr)
                    write_tf_file(wd / C.nodal_result_name(node, dof, "TFI", maxnode), f_grid[:kmax], Hg[:kmax, jj],
                                  cplx, hdr + " interp=%d" % int(d["interp"]))
                if plan.listing_tf[i]:
                    lst.section(f"Transfer function node {node} {tag} (computed frequencies)")
                    _printer_plot(lst, f_ssi, Hc[:, jj])
            if not full:
                continue
            # ---- convolution
            if not need.any():
                continue
            if seismic:
                Racc = Hg * cmot.A[:, None]                             # acceleration in g
                Rreq = Racc
            else:
                HF = Hg * cmot.A[:, None]                               # displacement
                Racc = -(w ** 2)[:, None] * HF
                Rreq = {0: HF, 1: (1j * w)[:, None] * HF, 2: Racc}[resp]
            acc_full = np.fft.irfft(S.hermitian_bins(Racc, nfft), n=nfft, axis=0)
            req_full = acc_full if Rreq is Racc else np.fft.irfft(S.hermitian_bins(Rreq, nfft), n=nfft, axis=0)
            for jj, c in enumerate(cols):
                i = b0 + jj
                if not need[jj]:
                    continue
                node, dof = int(plan.node[i]), int(plan.dof[i])
                tag = C.DOF_TAGS[dof]
                a_full = acc_full[:, jj]
                r_full = req_full[:, jj]
                if bl and (seismic or resp == 2):
                    br = S.baseline_correction(a_full, dt, scale=gravity if seismic else 1.0)
                    a_full = br.acc
                    r_full = a_full
                    vel_out = br.vel[:nout]
                    dis_out = br.dis[:nout]
                    baseline_report.append((node, tag, br.c1, br.c2, float(br.dis[nout - 1]), br.final_disp))
                    if fh13 is not None:
                        fh13.write(f"# node {node} {tag}: time(s) acc({'g' if seismic else 'L/s2'}) vel(L/s) disp(L)\n")
                        t = np.arange(nout) * dt
                        for row in zip(t, a_full[:nout], vel_out, dis_out):
                            fh13.write(" ".join(f"{v:.10e}" for v in row) + "\n")
                r_out = r_full[:nout]
                if plan.th[i]:
                    write_history_file(wd / C.nodal_result_name(node, dof, "ACC", maxnode), r_out, dt)
                    if fh12 is not None:
                        fh12.write(f"# node {node} {tag} {rname} history: dt {dt:g} n {nout}\n")
                        fh12.write("\n".join(f"{v:.10e}" for v in r_out) + "\n")
                if plan.th_print[i]:
                    lst.section(f"Time history node {node} {tag} ({rname}, every {int(d['step'])}th step)")
                    for k in range(0, nout, int(d["step"])):
                        lst.write(f"{k * dt:12.4f} {r_out[k]:16.6e}")
                if plan.mx[i] or plan.th[i]:
                    k = int(np.argmax(np.abs(r_out)))
                    maxima.append((node, tag, float(np.abs(r_out[k])), k * dt, bool(plan.mx[i])))
                if plan.frame[i] and int(d["rstacc"]):
                    fr_acc[int(c)] = r_out.copy()
                if need_rs and (plan.rs[i] or (plan.frame[i] and int(d["rstrs"]) and dof <= 3)):
                    rs = response_spectrum(a_full, dt, rs_f, damps)
                    if plan.rs[i]:
                        for idd in range(len(damps)):
                            write_xy_file(wd / C.nodal_rs_name(node, dof, idd + 1, maxnode), rs_f, [rs["SA"][idd]],
                                          header=f"node {node} {tag} damping {damps[idd]:g}: f(Hz) SA({unit})")
                        if fh12 is not None:
                            for idd in range(len(damps)):
                                fh12.write(f"# node {node} {tag} RS damping {damps[idd]:g}: f SA({unit}) SV({sv_unit})\n")
                                for fi, sa, sv in zip(rs_f, rs["SA"][idd], rs["SV"][idd] * sv_scale):
                                    fh12.write(f"{fi:.6f} {sa:.10e} {sv:.10e}\n")
                    if plan.frame[i]:
                        fr_rs[int(c)] = rs["SA"].copy()
                    if plan.rs[i] and plan.explicit[i]:
                        _rs_listing(lst, node, tag, rs_f, damps, rs, unit, sv_scale, sv_unit)
            ctx.progress(min(1.0, (b0 + len(cols)) / ncol), f"MOTION: {b0 + len(cols)}/{ncol} columns")
    finally:
        for fh in (fh12, fh13):
            if fh is not None:
                fh.close()

    if edu18:
        lst.warning(f"EDU-18 (G-19): |ATF| at the first SSI frequency {f_ssi[0]:.4g} Hz deviates more than 5 % from "
                    "the rigid-body value at " + ", ".join(f"{n} {t} ({a:.3f})" for n, t, a in edu18[:20])
                    + (" ..." if len(edu18) > 20 else "") + "; check the model and the lowest SSI frequency")
    if maxima:
        lst.section(f"Maximum requested response ({rname}, {unit})")
        lst.write(f"{'node':>8s} {'dof':>6s} {'max |r|':>16s} {'time (s)':>10s}")
        for node, tag, v, t, show in maxima:
            if show:
                lst.write(f"{node:8d} {tag:>6s} {v:16.6e} {t:10.4f}")
        zpa = max(maxima, key=lambda r: r[2])
        lst.write(f"   largest value {zpa[2]:.6e} at node {zpa[0]} {zpa[1]}, t = {zpa[3]:.4f} s")
    if baseline_report:
        lst.section("Baseline correction (Hudson-Housner): removed acceleration c1 + c2 t")
        lst.write(f"{'node':>8s} {'dof':>6s} {'c1':>14s} {'c2':>14s} {'d(T_out)':>14s} {'d(T_end)':>14s}")
        for node, tag, c1, c2, dout, dend in baseline_report:
            lst.write(f"{node:8d} {tag:>6s} {c1:14.6e} {c2:14.6e} {dout:14.6e} {dend:14.6e}")
        lst.write("   check that the final displacements are about zero; use RELDISP for relative displacements")
    _write_frames(ctx, d, src, plan, f_ssi, rs_f, damps, fr_acc, fr_rs, dt)


def _nyquist_real(Hg: np.ndarray, nfft: int) -> np.ndarray:
    """Nyquist bin real (D-MOT-03) without touching the vibration f = 0 value (H_1)."""
    H = np.array(Hg, dtype=complex, copy=True)
    if nfft % 2 == 0:
        H[-1] = H[-1].real
    return H


def _node_table(plan: _Plan, cols_subset, dofs=(1, 2, 3)):
    """Rows of nodes with any of ``dofs`` among the planned columns -> (nodes, {col: (row, k)})."""
    nodes = sorted(set(int(plan.node[i]) for i, c in enumerate(plan.cols) if int(c) in cols_subset
                       and int(plan.dof[i]) in dofs))
    pos = {n: r for r, n in enumerate(nodes)}
    where = {}
    for i, c in enumerate(plan.cols):
        if int(c) in cols_subset and int(plan.dof[i]) in dofs:
            where[int(c)] = (pos[int(plan.node[i])], dofs.index(int(plan.dof[i])))
    return np.array(nodes, dtype=np.int64), where


def _write_frames(ctx: ModuleContext, d, src: _TFSource, plan: _Plan, f_ssi, rs_f, damps, fr_acc, fr_rs,
                  dt: float) -> None:
    """Post-processing restart frames (spec 05c B.5.12, requirements section 1.7.4; D-FIL-03)."""
    wd = ctx.workdir
    lst = ctx.listing
    if int(d["rsttf"]):
        trans_cols = set(int(c) for i, c in enumerate(plan.cols) if plan.dof[i] <= 3)
        nodes, where = _node_table(plan, trans_cols)
        if len(nodes):
            cols = np.array(sorted(where))
            Hc = src.computed(cols)
            out = wd / "TFU"
            out.mkdir(exist_ok=True)
            for q, f in enumerate(f_ssi):
                vals = np.zeros((len(nodes), 6))
                for j, c in enumerate(cols):
                    r, k = where[int(c)]
                    vals[r, 2 * k] = Hc[q, j].real
                    vals[r, 2 * k + 1] = Hc[q, j].imag
                write_frame(out / f"TFU_{frame_freq(f)}_{q + 1:05d}", nodes, vals)
            lst.write(f"   TFU frames: {len(f_ssi)} files in TFU/ ({len(nodes)} nodes, Re/Im of X Y Z)")
    if int(d["rstacc"]) and fr_acc:
        for dofs, sub in (((1, 2, 3), "ACC"), ((4, 5, 6), "ACCR")):
            nodes, where = _node_table(plan, set(fr_acc), dofs)
            if not len(nodes):
                continue
            n = len(next(iter(fr_acc.values())))
            vals = np.zeros((n, len(nodes), 3))
            for c, (r, k) in where.items():
                vals[:, r, k] = fr_acc[c]
            out = wd / sub
            out.mkdir(exist_ok=True)
            for s in range(n):
                write_frame(out / f"{sub}_{frame_time(s * dt)}_{s + 1:05d}", nodes, vals[s])
            if sub == "ACC":
                write_frame(wd / "ACC_max.txt", nodes, np.max(np.abs(vals), axis=0))
            lst.write(f"   {sub} frames: {n} files in {sub}/ ({len(nodes)} nodes)"
                      + ("; ZPA frame ACC_max.txt" if sub == "ACC" else ""))
    if int(d["rstrs"]) and fr_rs:
        nodes, where = _node_table(plan, set(fr_rs))
        if len(nodes):
            out = wd / "RS"
            out.mkdir(exist_ok=True)
            for idd in range(len(damps)):
                vals = np.zeros((len(rs_f), len(nodes), 3))
                for c, (r, k) in where.items():
                    vals[:, r, k] = fr_rs[c][idd]
                for q, f in enumerate(rs_f):
                    write_frame(out / f"RS{idd + 1:02d}_{frame_freq(f)}_{q + 1:05d}", nodes, vals[q])
            lst.write(f"   RS frames: {len(rs_f) * len(damps)} files in RS/ ({len(nodes)} nodes)")


if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
