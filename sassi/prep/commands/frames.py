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

(REMOVEFREQ, the third frequency-refinement tool, is implemented with the other extension
commands in :mod:`sassi.prep.commands.extensions`.)
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ...plotting.lines import critfreq, framesel, numeric_rows
from ...plotting.state import combine_frames, modframes_text, read_frame, read_frame_list, write_frame
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
