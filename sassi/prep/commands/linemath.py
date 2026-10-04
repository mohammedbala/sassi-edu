"""Line-object commands: file I/O and line mathematics (manual 9.14; requirements 3.4.L, 4.13;
spec 10 sections 4.1-4.2; D-LIN-01..04, D-FIL-02, D-UI-16).

* READSPEC,<SpecFile>,<numLines>,<Line1>,...: spectrum-format columns -> lines (the frequency
  column is not counted; column 1 is named after the file, column k > 1 ``<file>[k]``).
* READTH,<THFile>,<Pair>,<Num>: time history (Pair 0: dt first then values; 1: time/value pairs).
* WRITESPEC,<SpecFile>,<Num1>,...: up to 50 lines on their union grid (exact text, 17 digits).
* WRITETH,<THFile>,<Num>: one-column history, ``dt = x2 - x1``.
* ADDITION, SUBTRACTION, LINECOMBIN, AVERAGE, SRSS, BROADEN: results on the union of the source
  abscissas, linear inside a line, constant outside (requirements 4.13).  The destination may be
  a source (computed first, then overwritten).  An undefined source line is an error (D-LIN-02).
* LBINCORS: parsed; the manual specifies no algorithm (D-LIN-03) -> error.

The maths lives in :mod:`sassi.plotting.lines` (no interpreter dependency).
"""
from __future__ import annotations

from typing import List

from ...plotting.lines import (DEFAULT_NAMES, Line, addition, average, broaden, clean_xy,
                               linear_combination, ranges_differ, read_spec_file, read_th_file, srss,
                               subtraction, write_spec_file, write_th_file)
from ...plotting.state import plot_state
from ..options import eduopt
from ..registry import CommandError
from .plotting import plot_command, refresh_line_plots

MAX_SOURCES = 100


def _store(c, lines: List[Line]) -> None:
    st = plot_state(c.interp)
    for L in lines:
        st.set_line(L)
    nums = [L.number for L in lines]
    st.emit("lines", None, numbers=nums)
    refresh_line_plots(c, nums)


def _sources(c, k_from: int, k_to=None) -> List[int]:
    k_to = c.nargs if k_to is None else k_to
    nums = [c.int(k) for k in range(k_from, k_to + 1) if c.given(k)]
    if not nums:
        raise CommandError("at least one source line is required")
    if len(nums) > MAX_SOURCES:
        raise CommandError(f"at most {MAX_SOURCES} source lines")
    return nums


def _result(c, dest: int, op: str, sources: List[Line], x, y) -> Line:
    kind = sources[0].kind
    if kind == "history" and len(sources) > 1 and ranges_differ(sources):
        c.warn("the source histories cover different time ranges: outside its range a line is held constant at "
               "its last value (not zero-padded)")
    return Line(dest, DEFAULT_NAMES[op], x, y, kind, source=c.line.strip())


def _combine(c, op: str, fn) -> None:
    dest = c.int(1, required=True, what="dest")
    nums = _sources(c, 2)
    st = plot_state(c.interp)
    src = st.get_lines(nums)
    x, y = fn(src)
    L = _result(c, dest, op, src, x, y)
    _store(c, [L])
    c.confirm(f"line {dest} = {op} of line(s) {', '.join(str(n) for n in nums)} ({L.n} points, '{L.name}')")


# ======================================================================================
# File I/O
# ======================================================================================
@plot_command("READSPEC")
def cmd_readspec(c):
    """READSPEC,<SpecFile>,<numLines>,<Line1>,...,<LineN>: load spectrum columns (frequency column not counted)."""
    if not c.given(1):
        raise CommandError("<SpecFile> is required")
    path = c.input_path(c.str(1))
    n = c.int(2, required=True, what="numLines")
    if n < 1:
        raise CommandError("<numLines> must be >= 1 (the frequency column is not counted)")
    refs = [c.int(k) for k in range(3, c.nargs + 1) if c.given(k)]
    if not refs:
        raise CommandError("line reference numbers are required")
    if len(refs) < n:
        extra = [refs[-1] + i for i in range(1, n - len(refs) + 1)]
        c.info(f"READSPEC: {n} columns, {len(refs)} line numbers: columns {len(refs) + 1}..{n} go to lines "
               + ", ".join(str(v) for v in extra))
        refs += extra
    elif len(refs) > n:
        c.warn(f"{len(refs)} line numbers for {n} columns; the extra numbers are ignored")
        refs = refs[:n]
    if len(set(refs)) != len(refs):
        raise CommandError("the line numbers must be different")
    x, ys = read_spec_file(path, n)
    out = []
    notes = []
    for k, (num, y) in enumerate(zip(refs, ys), start=1):
        xx, yy, nt = clean_xy(x, y)
        notes = nt or notes
        name = path.name if k == 1 else f"{path.name}[{k}]"
        out.append(Line(num, name, xx, yy, "spectrum", source=str(path)))
    for nt in notes:
        c.warn(f"{path.name}: {nt}")
    _store(c, out)
    c.confirm(f"{path.name}: {len(out)} line(s) {', '.join(str(v) for v in refs)} loaded ({len(out[0].x)} points)")


@plot_command("READTH", max_args=3)
def cmd_readth(c):
    """READTH,<THFile>,<Pair>,<Num>: load a time history (Pair 0 one column with dt first, 1 time/value pairs)."""
    if not c.given(1):
        raise CommandError("<THFile> is required")
    path = c.input_path(c.str(1))
    pair = c.int(2, required=True, what="Pair")
    num = c.int(3, required=True, what="Num")
    t, a = read_th_file(path, pair)
    t, a, notes = clean_xy(t, a)
    for nt in notes:
        c.warn(f"{path.name}: {nt}")
    _store(c, [Line(num, path.name, t, a, "history", source=str(path))])
    dt = t[1] - t[0] if len(t) > 1 else 0.0
    c.confirm(f"{path.name}: line {num} loaded ({len(t)} points, dt {dt:g})")


@plot_command("WRITESPEC")
def cmd_writespec(c):
    """WRITESPEC,<SpecFile>,<Num1>,...,[Num50]: write lines on the union of their abscissas."""
    if not c.given(1):
        raise CommandError("<SpecFile> is required")
    nums = _sources(c, 2)
    if len(nums) > 50:
        raise CommandError("at most 50 lines per file")
    lines = plot_state(c.interp).get_lines(nums)
    path = c.output_path(c.str(1))
    try:
        out, X = write_spec_file(path, lines)
    except OSError as exc:
        raise CommandError(f"cannot write {path}: {exc}") from None
    c.confirm(f"{out.name}: {len(lines)} line(s) written on {len(X)} abscissas")


@plot_command("WRITETH", max_args=2)
def cmd_writeth(c):
    """WRITETH,<THFile>,<Num>: write a line as a one-column history (dt = x2 - x1)."""
    if not c.given(1):
        raise CommandError("<THFile> is required")
    num = c.int(2, required=True, what="Num")
    L = plot_state(c.interp).get_lines([num])[0]
    path = c.output_path(c.str(1))
    try:
        out, dt, notes = write_th_file(path, L)
    except OSError as exc:
        raise CommandError(f"cannot write {path}: {exc}") from None
    for nt in notes:
        c.warn(nt)
    c.confirm(f"{out.name}: line {num} written ({L.n} values, dt {dt:g})")


# ======================================================================================
# Line mathematics
# ======================================================================================
@plot_command("ADDITION")
def cmd_addition(c):
    """ADDITION,<dest>,<source1>,...,<source100>: sum of lines ('Linear Combin.')."""
    _combine(c, "ADDITION", addition)


@plot_command("SUBTRACTION")
def cmd_subtraction(c):
    """SUBTRACTION,<dest>,<source1>,...,<source100>: first line minus the others ('Linear Combin.')."""
    _combine(c, "SUBTRACTION", subtraction)


@plot_command("AVERAGE")
def cmd_average(c):
    """AVERAGE,<dest>,<source1>,...,<source100>: mean of lines ('Average Line')."""
    _combine(c, "AVERAGE", average)


@plot_command("SRSS")
def cmd_srss(c):
    """SRSS,<dest>,<source1>,...,<source100>: square root of the sum of squares ('SRSS Line')."""
    _combine(c, "SRSS", srss)


@plot_command("LINECOMBIN")
def cmd_linecombin(c):
    """LINECOMBIN,<Dest>,<Line1>,<Coeff1>,...,[Line100],[Coeff100]: sum of c_i y_i ('Linear Combin.')."""
    dest = c.int(1, required=True, what="Dest")
    n = c.nargs - 1
    if n < 2 or n % 2:
        raise CommandError("a coefficient is required for every line: LINECOMBIN,<Dest>,<Line1>,<Coeff1>,...")
    if n // 2 > MAX_SOURCES:
        raise CommandError(f"at most {MAX_SOURCES} lines")
    nums = [c.int(k, required=True) for k in range(2, c.nargs + 1, 2)]
    coeffs = [c.float(k, required=True) for k in range(3, c.nargs + 1, 2)]
    src = plot_state(c.interp).get_lines(nums)
    x, y = linear_combination(src, coeffs)
    L = _result(c, dest, "LINECOMBIN", src, x, y)
    _store(c, [L])
    terms = " + ".join(f"{a:g} * line {b}" for a, b in zip(coeffs, nums))
    c.confirm(f"line {dest} = {terms} ({L.n} points)")


@plot_command("BROADEN")
def cmd_broaden(c):
    """BROADEN,<Dest>,<Smooth1>,<Smooth2>,<source1>,...: envelope, +-Smooth2 % peak broadening, Smooth1 % bridging."""
    dest = c.int(1, required=True, what="Dest")
    s1 = c.float(2, required=True, what="Smooth1")
    s2 = c.float(3, required=True, what="Smooth2")
    nums = _sources(c, 4)
    src = plot_state(c.interp).get_lines(nums)
    grid = eduopt(c.model, "BROADENGRID")
    bridge = eduopt(c.model, "BRIDGE")
    x, y, notes = broaden(src, s1, s2, grid=grid, bridge=bridge)
    for nt in notes:
        c.warn(nt)
    L = _result(c, dest, "BROADEN", src, x, y)
    _store(c, [L])
    c.confirm(f"line {dest} = envelope of line(s) {', '.join(str(n) for n in nums)}, broadened +-{s2:g} %, "
              f"bridged {s1:g} % ({L.n} points, grid {grid.lower()}, bridging {bridge.lower()})")


@plot_command("LBINCORS", max_args=2)
def cmd_lbincors(c):
    """LBINCORS,<out>,<in>: lower-bound incoherent RS -- no algorithm in the manual (D-LIN-03)."""
    c.int(1, required=True, what="out")
    c.int(2, required=True, what="in")
    raise CommandError("the lower-bound incoherent response spectrum algorithm is not specified in the manual "
                       "(D-LIN-03); no line was created")
