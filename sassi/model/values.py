"""Number parsing and formatting rules of the ``.pre`` language (shared by model and interpreter).

Requirements section 3.1:

* L6 / D-PAR-07 -- numbers are free format; Fortran exponents ``E`` and ``D`` are accepted
  (``1.0D9``, ``1e+008``).  Integer fields given as reals are rounded to the nearest integer
  (the caller warns when the fractional part is not zero).
* D-PAR-17 -- WRITE output must be byte-stable and must read back to the *identical* state.
  :func:`fmt_num` therefore writes the shortest text that parses back to exactly the same
  binary float (Python's ``repr``), and writes integral values without a decimal point.

These helpers live in :mod:`sassi.model` (not in :mod:`sassi.prep`) because typed option
records convert their stored tokens with them, and the model package must not import the
interpreter.
"""
from __future__ import annotations

import math
import re
from typing import Optional, Union

Number = Union[int, float]

_NUM_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([EeDd][+-]?\d+)?$")


class NumberError(ValueError):
    """A token that should be a number is not one."""


def is_number(text: Optional[str]) -> bool:
    """True when ``text`` is a free-format number (L6)."""
    if text is None:
        return False
    t = text.strip()
    if not t:
        return False
    if _NUM_RE.match(t):
        return True
    return t.lower() in ("inf", "+inf", "-inf", "nan")


def parse_float(text: str) -> float:
    """Parse a free-format real number; ``D`` exponents are accepted (L6, D-PAR-07)."""
    t = text.strip()
    if not t:
        raise NumberError("empty numeric field")
    if not _NUM_RE.match(t):
        if t.lower() in ("inf", "+inf", "-inf", "nan"):
            return float(t)
        raise NumberError(f"'{text.strip()}' is not a number")
    return float(t.replace("D", "E").replace("d", "e"))


def round_half_away(x: float) -> int:
    """Round to the nearest integer, halves away from zero (Fortran NINT, not banker's rounding)."""
    if x >= 0:
        return int(math.floor(x + 0.5))
    return -int(math.floor(-x + 0.5))


def parse_int(text: str):
    """Parse an integer field.

    Returns ``(value, exact)`` where ``exact`` is False when the field was a real number with
    a non-zero fractional part (D-PAR-07: rounded, the caller issues a warning).
    """
    x = parse_float(text)
    if not math.isfinite(x):
        raise NumberError(f"'{text.strip()}' is not a finite integer")
    n = round_half_away(x)
    return n, (float(n) == x)


def fmt_num(x: Number) -> str:
    """Shortest exact text for a number (WRITE, D-PAR-17).

    Integral floats are written without a decimal point (``10.0 -> '10'``) as long as they are
    exactly representable as integers; other floats use ``repr`` (shortest round-trip form).
    """
    if isinstance(x, bool):
        return "1" if x else "0"
    if isinstance(x, int):
        return str(x)
    xf = float(x)
    if math.isfinite(xf) and xf == int(xf) and abs(xf) < 1e15:
        s = str(int(xf))
        return "-0" if (s == "0" and math.copysign(1.0, xf) < 0) else s
    return repr(xf)


def fmt_fixed(x: float, width: int = 13, prec: int = 5) -> str:
    """Fixed-width number for listings (NLIST, MLIST ...)."""
    if x == 0:
        return f"{0.0:>{width}.{prec}f}"
    ax = abs(x)
    if 1e-3 <= ax < 10 ** (width - prec - 3):
        return f"{x:>{width}.{prec}f}"
    return f"{x:>{width}.{prec - 1}E}"
