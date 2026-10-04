"""Frame and frequency-refinement tools (manual 9.7.4, 9.7.11, 9.7.12, 9.7.24; requirements 3.4.I;
spec 09 sections 2.4, 2.11, 2.12, 2.24; D-FIL-03).

* CRITFREQ,<tol>,<minfilter>,<TF>,<Var>: compare the interpolated ATF (``<TF>.TFI``, Fourier
  frequencies) with the computed one (``<TF>.TFU``, SSI frequencies).  An interpolated peak that is
  more than ``tol`` % above (or below) the larger computed value of the two bracketing SSI
  frequencies is not supported by the SSI solution: a frequency should be added there.  Peaks lower
  than ``(1 - minfilter/100)`` times the global TFI maximum are not considered.  The variable
  ``<Var>`` receives the **frequency numbers** of the flagged peaks (``round(f / df)``, df = TFI
  step = FILE8 df, D-MOT-14), so ``FREQ,<set>,@Var[1],...`` adds them directly; the Hz values are
  printed.  This is the frequency-refinement loop of requirements 2.4.
* FRAMECOMBIN,<op>,<num>,<InFile1>,...,<InFileX>,<Outfile>: 0 SRSS, 1 sum, 2 average of frame
  files, cell by cell (node ids copied from the first file).
* FRAMESEL,<tol>,<Acc>,<Var>: 1-based frame numbers of the local maxima / minima of an acceleration
  history with ``|a| >= tol/100 max|a|`` -> ``<Var>`` (frame k = time (k-1) dt, D-STR-08).
* MODFRAMES,<cols>,<framelist>: rewrite the second header number (columns) of every frame of a list.
* HARMFRAME,<Src>,<Freq>,<OutDir>,[NFrames],[Ref] (SASSI-EDU extension): the steady-state harmonic
  motion of every node at one computed (SSI) frequency, ``u(t) = Re(H(f) e^{iwt})`` over one period, as
  frames ``HARM_<wt deg>_<k>`` (node, X, Y, Z) for PROCFRAME and DEFORMPLOT.  Source: a FILE8-type
  file of ANALYS / COMBIN or the ``TFU`` restart frames of MOTION (Restart for TF).

(REMOVEFREQ, the third frequency-refinement tool, is implemented with the other extension
commands in :mod:`sassi.prep.commands.extensions`.)
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path

import numpy as np

from ...plotting.lines import critfreq, framesel, numeric_rows
from ...plotting.state import (HARM_TAG, combine_frames, frame_name_info, harmonic_frame_name, harmonic_frames,
                               modframes_text, read_frame, read_frame_list, write_frame)
from ..interpreter import Variable
from ..registry import CommandError
from .plotting import plot_command

_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _variable(c, k: int) -> str:
    name = c.raw(k)
    if name.startswith("@"):
        name = name[1:]
    if not name or not _NAME_RE.match(name):
        raise CommandError(f"argument {k}: '{c.raw(k)}' is not a valid variable name")
    return name


def _set_variable(c, name: str, items) -> Variable:
    """Write a result list into a variable (redefinition: counter reset to 0, like VAR)."""
    v = Variable(name, [str(x) for x in items], 0)
    c.interp.variables[name.lower()] = v
    return v


def _tf_file(c, base: str, ext: str) -> Path:
    for e in (ext, ext.lower()):
        try:
            return c.interp.resolve_path(base + "." + e, must_exist=True)
        except CommandError:
            continue
    raise CommandError(f"{base}.{ext} not found (CRITFREQ reads <TF>.TFU and <TF>.TFI written by MOTION)")


def _amplitudes(path: Path):
    rows = [r for r in numeric_rows(path) if len(r) >= 2]
    if len(rows) < 2:
        raise CommandError(f"{path.name}: fewer than 2 (frequency, amplitude) rows")
    A = np.asarray([r[:2] for r in rows], dtype=float)
    order = np.argsort(A[:, 0], kind="stable")
    return A[order, 0], np.abs(A[order, 1])


@plot_command("CRITFREQ", max_args=4)
def cmd_critfreq(c):
    """CRITFREQ,<tol>,<minfilter>,<TF>,<Var>: frequencies where interpolated TF peaks deviate from TFU."""
    tol = c.float(1, required=True, what="tol")
    minf = c.float(2, required=True, what="minfilter")
    if tol < 0 or not 0 <= minf <= 100:
        raise CommandError("<tol> must be >= 0 and <minfilter> within 0..100 %")
    base = c.str(3)
    if not base:
        raise CommandError("<TF> (transfer-function file without extension) is required")
    if base.upper().endswith((".TFU", ".TFI")):
        base = base[:-4]
    var = _variable(c, 4)
    fu, au = _amplitudes(_tf_file(c, base, "TFU"))
    fi, ai = _amplitudes(_tf_file(c, base, "TFI"))
    peaks, df = critfreq(fu, au, fi, ai, tol, minf)
    flagged = [p for p in peaks if p.flagged]
    numbers = sorted(set(p.number for p in flagged))
    _set_variable(c, var, numbers)
    for p in peaks:
        c.info(f"CRITFREQ: TFI peak {p.freq:g} Hz (frequency number {p.number}): |TFI| {p.tfi:.5g}, computed "
               f"neighbours max {p.ref:.5g}, difference {p.diff_pct:.1f} %" + ("  -> added" if p.flagged else ""))
    c.confirm(f"{len(peaks)} TFI peak(s) considered, {len(numbers)} flagged at tol {tol:g} %: variable {var} = "
              + (", ".join(str(n) for n in numbers) or "(empty)") + f" (frequency numbers, df {df:g} Hz)")


@plot_command("FRAMECOMBIN")
def cmd_framecombin(c):
    """FRAMECOMBIN,<op>,<num>,<InFile1>,...,<InFileX>,<Outfile>: combine frames (0 SRSS, 1 sum, 2 average)."""
    op = c.int(1, required=True, what="op")
    if op not in (0, 1, 2):
        raise CommandError("<op> must be 0 (SRSS), 1 (sum) or 2 (average)")
    num = c.int(2, required=True, what="num")
    if num < 1:
        raise CommandError("<num> must be >= 1")
    files = [c.str(k) for k in range(3, c.nargs + 1)]
    if any(not f for f in files):
        raise CommandError("blank file name")
    if len(files) != num + 1:
        raise CommandError(f"<num> = {num} needs {num} input files and one output file ({len(files)} names given)")
    frames = [read_frame(c.input_path(f)) for f in files[:-1]]
    out = combine_frames(frames, op)
    path = c.output_path(files[-1])
    try:
        write_frame(path, out.nodes, out.values)
    except OSError as exc:
        raise CommandError(f"cannot write {path}: {exc}") from None
    c.confirm(f"{path.name}: {('SRSS', 'sum', 'average')[op]} of {num} frame(s) ({len(out.nodes)} rows, "
              f"{out.ncols} columns)")


@plot_command("FRAMESEL", max_args=3)
def cmd_framesel(c):
    """FRAMESEL,<tol>,<Acc>,<Var>: critical frame numbers (local extrema >= tol % of max|a|) -> variable."""
    tol = c.float(1, required=True, what="tol")
    if not 0 <= tol <= 100:
        raise CommandError("<tol> must be within 0..100 %")
    if not c.given(2):
        raise CommandError("<Acc> (acceleration history file) is required")
    path = c.input_path(c.str(2))
    var = _variable(c, 3)
    rows = numeric_rows(path)
    if len(rows) > 1 and all(len(r) == 2 for r in rows):
        a = np.asarray([r[1] for r in rows], dtype=float)           # time / value pairs
    else:
        flat = [v for r in rows for v in r]
        if len(flat) < 2:
            raise CommandError(f"{path.name}: no acceleration values (MOTION .ACC: dt first, then the values)")
        a = np.asarray(flat[1:], dtype=float)                      # one-column: dt first
    frames = framesel(a, tol)
    _set_variable(c, var, frames.tolist())
    c.confirm(f"{len(frames)} critical frame(s) of {len(a)}: variable {var} = "
              + (", ".join(str(int(k)) for k in frames[:30]) + (" ..." if len(frames) > 30 else "")
                 if len(frames) else "(empty)"))


@plot_command("MODFRAMES", max_args=2)
def cmd_modframes(c):
    """MODFRAMES,<cols>,<framelist>: set the column count (second header number) of every listed frame."""
    cols = c.int(1, required=True, what="cols")
    if cols < 2:
        raise CommandError("<cols> must be >= 2 (the node id is column 1)")
    if not c.given(2):
        raise CommandError("<framelist> is required")
    files = read_frame_list(c.input_path(c.str(2)))
    changed = 0
    for f in files:
        if not f.exists():
            raise CommandError(f"frame file {f} not found")
        text = f.read_text(encoding="utf-8", errors="replace")
        new, ch = modframes_text(text, cols)
        if ch:
            f.write_text(new, encoding="utf-8")
            changed += 1
    c.confirm(f"{changed} of {len(files)} frame header(s) rewritten with {cols} columns")


# ======================================================================================
# HARMFRAME (SASSI-EDU extension): steady-state harmonic frames at one SSI frequency
# ======================================================================================
#: frames per period accepted by HARMFRAME
HARM_FRAMES_MIN, HARM_FRAMES_MAX, HARM_FRAMES_DEFAULT = 4, 360, 24


def _harm_from_file8(path: Path, freq: float):
    """(SSI frequencies, index q, description, nodes, H (n, 3) at f_q, meta) of a FILE8-type file."""
    from ...modules.base import ModuleError
    from ...modules.motion import read_file8
    try:
        f8 = read_file8(path, path.name)
    except ModuleError as exc:
        raise CommandError(str(exc)) from None
    fnum = np.asarray(f8["fnum"], dtype=np.int64)
    f_ssi = fnum * float(f8.meta["df"])
    q = int(np.argmin(np.abs(f_ssi - freq)))
    eq_node = np.asarray(f8["eq_node"], dtype=np.int64)
    eq_dof = np.asarray(f8["eq_dof"], dtype=np.int64)
    trans = np.nonzero(eq_dof <= 3)[0]
    nodes = np.array(sorted(set(int(n) for n in eq_node[trans])), dtype=np.int64)
    pos = {int(n): r for r, n in enumerate(nodes)}
    Hq = np.asarray(f8["H"])[q]
    H = np.zeros((len(nodes), 3), dtype=complex)
    for col in trans:
        H[pos[int(eq_node[col])], int(eq_dof[col]) - 1] = Hq[col]
    return f_ssi, q, f"{path.name} frequency number {int(fnum[q])}", nodes, H, f8.meta


def _harm_from_tfu(path: Path, freq: float):
    """The same from the TFU restart frames of MOTION (``TFU_<f>_<n>``: node, Re/Im of X, Y, Z)."""
    cands = [(frame_name_info(f.name)["value"], f) for f in read_frame_list(path)]
    cands = [(v, f) for v, f in cands if v is not None]
    if not cands:
        raise CommandError(f"{path.name}: no frame file named <tag>_<frequency>_<n> (MOTION Restart for TF writes "
                           f"TFU/TFU_<f>_<n>); a FILE8-type file is read directly")
    f_ssi = np.array([v for v, _ in cands], dtype=float)
    q = int(np.argmin(np.abs(f_ssi - freq)))
    fr = read_frame(cands[q][1])
    if fr.values.shape[1] != 6:
        raise CommandError(f"{cands[q][1].name}: {fr.values.shape[1]} data columns; TFU frames hold 6 (Re, Im of X, "
                           f"Y, Z)")
    H = fr.values[:, 0::2] + 1j * fr.values[:, 1::2]
    order = np.argsort(f_ssi, kind="stable")
    return f_ssi[order], int(np.nonzero(order == q)[0][0]), f"TFU frame {cands[q][1].name}", fr.nodes, H, None


def _model_numbering(c, src: Path, nodes: np.ndarray):
    """FILE8 node numbers -> model numbers when HOUSE renumbered the model (``<model>.map`` beside the
    source): DEFORMPLOT draws the frames on the model, which keeps its own numbering."""
    name = c.model.name
    folder = src.parent
    mp = folder / f"{name}.map" if name else None
    if mp is None or not mp.is_file():
        return nodes, ""
    from ...core.renumber import read_map
    try:
        new_old = {int(n): int(o) for o, n in read_map(mp).items()}
    except (OSError, ValueError) as exc:
        c.warn(f"HOUSE optimizer map {mp.name} unreadable ({exc}): frame nodes keep the FILE8 numbers")
        return nodes, ""
    out = np.array([new_old.get(int(n), int(n)) for n in nodes], dtype=np.int64)
    changed = int(np.sum(out != nodes))
    return out, (f"; {changed} node number(s) translated to the model numbering ({mp.name})" if changed else "")


@plot_command("HARMFRAME", max_args=5)
def cmd_harmframe(c):
    """HARMFRAME,<Src>,<Freq>,<OutDir>,[NFrames],[Ref]: steady-state harmonic frames u = Re(H e^{iwt}) of
    every node over one period at the computed frequency closest to <Freq> (SASSI-EDU extension; PROCFRAME
    and DEFORMPLOT animate them).

    ``<Src>``: a FILE8-type file (FILE8, FILE8X, FILE81 ...) or the TFU frame folder (or list file) of
    MOTION Restart for TF.  ``<Freq>`` in Hz: the closest SSI frequency is used and reported.
    ``<OutDir>``: folder of the frames (created; earlier HARM frames in it are replaced).  ``<NFrames>``:
    frames per period (default 24, 4..360).  ``<Ref>``: blank = total motion per unit control motion;
    0 = relative to the free field (the control motion: unit amplitude, zero phase, in the input
    direction, as RELDISP without RELFILE); n = relative to node n (its X, Y, Z motion subtracted).
    """
    if not c.given(1) or not c.given(3):
        raise CommandError("<Src> (FILE8-type file or TFU frame folder) and <OutDir> (frame folder) are required")
    src = c.input_path(c.str(1))
    freq = c.float(2, required=True, what="Freq")
    if not freq > 0:
        raise CommandError("<Freq> must be > 0 Hz")
    out = c.output_path(c.str(3))
    nfr = c.int(4, default=HARM_FRAMES_DEFAULT)
    if not HARM_FRAMES_MIN <= nfr <= HARM_FRAMES_MAX:
        raise CommandError(f"<NFrames> must be within {HARM_FRAMES_MIN}..{HARM_FRAMES_MAX} (frames per period)")
    ref = c.int(5) if c.given(5) else None
    if ref is not None and ref < 0:
        raise CommandError("<Ref> must be blank (total motion), 0 (free field) or a node number")
    file8 = src.is_file() and zipfile.is_zipfile(src)
    f_ssi, q, what, nodes, H, meta = (_harm_from_file8 if file8 else _harm_from_tfu)(src, float(freq))
    if not len(nodes):
        raise CommandError(f"{src.name}: no translational degree of freedom")
    nodes, renum = _model_numbering(c, src, nodes)
    fq = float(f_ssi[q])
    if freq < f_ssi.min() * (1 - 1e-9) or freq > f_ssi.max() * (1 + 1e-9):
        c.warn(f"{freq:g} Hz is outside the computed range {f_ssi.min():.4g}..{f_ssi.max():.4g} Hz: the nearest "
               f"computed frequency {fq:.4g} Hz is used")
    elif abs(fq - freq) > 0.05 * freq:
        c.warn(f"the closest computed frequency {fq:.4g} Hz is more than 5 % from {freq:g} Hz (add an SSI frequency "
               f"there to animate it)")
    # ---------------------------------------------------------------- reference motion
    if ref is None:
        refname = "total motion"
    elif ref == 0:
        if meta is None:
            raise CommandError("<Ref> = 0 (free field) needs a FILE8-type source (control direction); use a node "
                               "number or blank with TFU frames")
        if int(meta.get("type", 0)) != 0:
            refname = "total motion (foundation vibration: no free-field motion, the free-field reference is 0)"
        else:
            from ...core.interp import rigid_body_anchor
            cm, ang = int(meta.get("cm", 0)), float(meta.get("ang", 0.0))
            H = H - np.array([rigid_body_anchor(d, cm, ang) for d in (1, 2, 3)], dtype=complex)[None, :]
            refname = f"relative to the free field (unit control motion in {('x', 'y', 'z')[cm]}')"
    else:
        hit = np.nonzero(nodes == ref)[0]
        if not len(hit):
            raise CommandError(f"reference node {ref} has no translational DOF in {src.name}")
        H = H - H[hit[0]][None, :]
        refname = f"relative to node {ref}"
    # ---------------------------------------------------------------- frames
    phases, U = harmonic_frames(H, nfr)
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob(f"{HARM_TAG}_*"):
        if old.is_file():
            old.unlink()
    others = [p.name for p in out.iterdir() if p.is_file() and not p.name.startswith(".")]
    if others:
        c.warn(f"{out.name} holds other files ({others[0]} ...): PROCFRAME of the folder would mix them with the "
               f"harmonic frames; use a folder of its own")
    try:
        for k in range(nfr):
            write_frame(out / harmonic_frame_name(float(phases[k]), k + 1), nodes, U[k])
    except OSError as exc:
        raise CommandError(f"cannot write the frames in {out}: {exc}") from None
    amp = np.abs(H)
    r, d = np.unravel_index(int(np.argmax(amp)), amp.shape)
    c.info(f"HARMFRAME: {fq:.6g} Hz ({what}), the computed frequency closest to {freq:g} Hz among {len(f_ssi)}; "
           f"period {1.0 / fq:.4g} s; {refname}{renum}")
    try:
        shown = out.resolve().relative_to(Path(c.interp.cwd).resolve()).as_posix()
    except ValueError:
        shown = str(out)
    deg = int(round(float(np.degrees(np.angle(H[r, d]))))) + 0          # + 0: no "-0"
    c.confirm(f"{nfr} frames over one period (wt = 0..{phases[-1]:g} deg) in {shown} ({len(nodes)} nodes, X Y Z); "
              f"largest amplitude {amp[r, d]:.4g} at node {int(nodes[r])} {'XYZ'[d]} (phase {deg} deg): PROCFRAME "
              f"the folder, then DEFORMPLOT")
