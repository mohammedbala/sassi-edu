"""Line objects and line mathematics (requirements section 4.13, D-LIN-01..04, D-FIL-02, D-UI-16).

A *line object* is the ACS SASSI name for one plotted curve: an ordered list of (x, y) points with
a name and a markers flag, stored in the session under an integer reference number (spec 06
section 5.1, spec 10 section 1.2).  Spectra and transfer functions have x = frequency (Hz); time
histories have x = time (s).

Every multi-line operation (ADDITION, SUBTRACTION, LINECOMBIN, AVERAGE, SRSS, BROADEN, WRITESPEC)
works on the **union of the abscissas** of its sources (requirements 4.13, manual 9.14):

1. the result grid is ``X = sort(unique(union of all x))`` (values closer than 1e-9 relative are
   merged, the first one kept);
2. inside a line's own range its value is the **linear** interpolation between its points (linear
   x and linear y, even when the plot axes are logarithmic);
3. outside its range the line is held **constant** at its first / last value
   (``numpy.interp(X, xs, ys, left=ys[0], right=ys[-1])``).

Rule 3 is what the manual specifies.  For time histories of different durations it means the
shorter record *repeats its last sample* -- not the same as zero padding; the commands warn when
history ranges differ.

Spectrum broadening (D-LIN-01, ASCE 4 / RG 1.122 practice) is done in three steps:

1. **envelope** ``E(f) = max_i y_i(f)``;
2. **peak broadening** by ``b = Smooth2/100``: every ordinate at ``f0`` is spread over
   ``[f0(1-b), f0(1+b)]``, i.e. ``B(f) = max{E(f') : f/(1+b) <= f' <= f/(1-b)}``.  A sharp peak at
   ``f0`` becomes a plateau over ``[f0(1-b), f0(1+b)]`` and its flanks shift by the factors
   ``(1-b)`` and ``(1+b)``;
3. **peak bridging** by ``p = Smooth1/100``: valleys between adjacent peaks P_k, P_k+1 are
   filled up to ``min(P_k, P_k+1)`` when ``(min(P_k, P_k+1) - valley)/min(P_k, P_k+1) <= p``
   (``EDUOPT,BRIDGE,AMPLITUDE``: when ``|P_k - P_k+1| / max(P_k, P_k+1) <= p``); repeated
   until nothing changes.

With the default ``EDUOPT,BROADENGRID,AUGMENTED`` the result is the *exact* piecewise-linear
function: the grid is augmented with the points ``f(1 +- b)`` (plateau edges), the crossing points
of the source lines (envelope corners) and the crossing points created by the broadening and
bridging steps.  ``EDUOPT,BROADENGRID,ACS`` reproduces the discrete form of the original program:
every step is evaluated on the source points only (hence the manual's advice to use >= 301
frequencies, EDU-15).

The module is pure numpy (no interpreter or plotting imports), so the GUI and tests use it directly.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

PathLike = Union[str, Path]

#: relative tolerance for merging near-equal abscissas of the union grid (spec 10 section 1.3)
MERGE_RTOL = 1e-9
#: minimum number of points recommended for broadening (manual line 6887, SRP 3.7.1, EDU-15)
BROADEN_MIN_POINTS = 301
#: default names given by the creating commands (spec 10 section 1.2, exact strings)
DEFAULT_NAMES = {"ADDITION": "Linear Combin.", "SUBTRACTION": "Linear Combin.", "LINECOMBIN": "Linear Combin.",
                 "AVERAGE": "Average Line", "SRSS": "SRSS Line", "BROADEN": "Envelope"}
KINDS = ("spectrum", "history")


class LineError(ValueError):
    """Invalid line data or arguments (reported by the command handlers as errors)."""


# ======================================================================================
# Line objects
# ======================================================================================
@dataclass
class Line:
    """One line object of the session line store (spec 10 section 1.2).

    ``x`` is strictly increasing.  ``kind`` is ``'spectrum'`` (x = frequency, SPECPLOT, WRITESPEC)
    or ``'history'`` (x = time, THPLOT, WRITETH); it only selects defaults, any line may be
    plotted in either plot.  ``name`` and ``markers`` are *global* line properties (LINENAME,
    MARKERS): every graph showing the line uses them.
    """
    number: int
    name: str
    x: np.ndarray
    y: np.ndarray
    kind: str = "spectrum"
    source: str = ""
    markers: bool = False

    def __post_init__(self):
        self.x = np.asarray(self.x, dtype=float).copy()
        self.y = np.asarray(self.y, dtype=float).copy()
        if self.x.ndim != 1 or self.x.shape != self.y.shape:
            raise LineError(f"line {self.number}: x and y must be 1-D arrays of equal length")
        if len(self.x) == 0:
            raise LineError(f"line {self.number}: no points")
        if len(self.x) > 1 and np.any(np.diff(self.x) <= 0):
            raise LineError(f"line {self.number}: abscissas must be strictly increasing")
        if self.kind not in KINDS:
            raise LineError(f"line {self.number}: kind must be one of {KINDS}")

    @property
    def n(self) -> int:
        return len(self.x)

    def at(self, X) -> np.ndarray:
        """Values on ``X``: linear inside the line's range, constant outside (requirements 4.13)."""
        return resample(self.x, self.y, X)

    def copy(self, number: Optional[int] = None) -> "Line":
        return Line(self.number if number is None else number, self.name, self.x, self.y, self.kind,
                    self.source, self.markers)

    def to_dict(self) -> Dict:
        """JSON-serialisable form (GUI / ``PlotState.to_dict``)."""
        return {"number": int(self.number), "name": self.name, "kind": self.kind, "source": self.source,
                "markers": bool(self.markers), "n": int(self.n),
                "x": [float(v) for v in self.x], "y": [float(v) for v in self.y]}

    @classmethod
    def from_dict(cls, d: Dict) -> "Line":
        return cls(int(d["number"]), d.get("name", ""), np.asarray(d["x"], float), np.asarray(d["y"], float),
                   d.get("kind", "spectrum"), d.get("source", ""), bool(d.get("markers", False)))


def clean_xy(x: Sequence[float], y: Sequence[float]) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Make ``x`` strictly increasing: stable sort, keep the first of duplicate abscissas.

    Returns ``(x, y, notes)``; ``notes`` describe what was changed (spec 06 section 5.4, spec 10
    section 1.3: sort on load, keep the first occurrence of a duplicate and warn).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    notes: List[str] = []
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise LineError("non-finite values in the data")
    if len(x) > 1 and np.any(np.diff(x) < 0):
        order = np.argsort(x, kind="stable")
        x, y = x[order], y[order]
        notes.append("abscissas were not ascending; the points were sorted")
    if len(x) > 1:
        keep = np.concatenate([[True], np.diff(x) > 0])
        if not np.all(keep):
            notes.append(f"{int((~keep).sum())} duplicate abscissa(s); the first occurrence kept")
            x, y = x[keep], y[keep]
    return x, y, notes


# ======================================================================================
# Union grid and resampling (requirements 4.13)
# ======================================================================================
def merge_close(X: np.ndarray, rtol: float = MERGE_RTOL) -> np.ndarray:
    """Sorted unique values with near-duplicates (relative ``rtol``) merged, the first one kept."""
    X = np.unique(np.asarray(X, dtype=float))
    if len(X) < 2:
        return X
    scale = np.maximum(np.abs(X[1:]), np.abs(X[:-1]))
    close = np.diff(X) <= rtol * np.where(scale > 0, scale, 1.0)
    if not np.any(close):
        return X
    out = [X[0]]
    for v, c in zip(X[1:], close):
        # compare with the last *kept* value so a chain of close values collapses to its first member
        if abs(v - out[-1]) <= rtol * max(abs(v), abs(out[-1]), 1e-300):
            continue
        out.append(v)
    return np.asarray(out)


def union_grid(xs: Iterable[np.ndarray], rtol: float = MERGE_RTOL) -> np.ndarray:
    """Union of the abscissas of several lines (spec 10 section 1.3, rule 1)."""
    parts = [np.asarray(x, dtype=float).ravel() for x in xs]
    if not parts:
        return np.zeros(0)
    return merge_close(np.concatenate(parts), rtol)


def resample(xs: np.ndarray, ys: np.ndarray, X) -> np.ndarray:
    """Piecewise-linear value of the line (xs, ys) at ``X`` with constant extrapolation (rules 2-3)."""
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    return np.interp(np.asarray(X, dtype=float), xs, ys, left=ys[0], right=ys[-1])


def on_union(lines: Sequence[Line]) -> Tuple[np.ndarray, np.ndarray]:
    """``(X, Y)`` with ``Y[i]`` = line i resampled on the union grid ``X`` of all lines."""
    if not lines:
        raise LineError("no source lines")
    X = union_grid([L.x for L in lines])
    Y = np.vstack([L.at(X) for L in lines])
    return X, Y


def ranges_differ(lines: Sequence[Line], rtol: float = 1e-9) -> bool:
    """True when the sources do not all span the same x range (constant extrapolation used)."""
    lo = [L.x[0] for L in lines]
    hi = [L.x[-1] for L in lines]
    span = max(max(hi) - min(lo), 1e-300)
    return (max(lo) - min(lo)) > rtol * span or (max(hi) - min(hi)) > rtol * span


# ======================================================================================
# Linear operations (requirements 4.13; spec 10 sections 4.1.1, 4.1.2, 4.1.5-4.1.7)
# ======================================================================================
def linear_combination(lines: Sequence[Line], coeffs: Sequence[float]) -> Tuple[np.ndarray, np.ndarray]:
    """LINECOMBIN: ``y = sum c_i y_i`` on the union grid."""
    if len(coeffs) != len(lines):
        raise LineError("one coefficient per line is required")
    X, Y = on_union(lines)
    c = np.asarray(coeffs, dtype=float)
    return X, np.sum(c[:, None] * Y, axis=0)          # element-wise (no BLAS: exact sums, no spurious flags)


def addition(lines: Sequence[Line]) -> Tuple[np.ndarray, np.ndarray]:
    """ADDITION: ``y = sum y_i`` (LINECOMBIN with all coefficients 1)."""
    return linear_combination(lines, [1.0] * len(lines))


def subtraction(lines: Sequence[Line]) -> Tuple[np.ndarray, np.ndarray]:
    """SUBTRACTION: ``y = y_1 - sum_{i>=2} y_i`` (coefficients 1, -1, -1, ...)."""
    return linear_combination(lines, [1.0] + [-1.0] * (len(lines) - 1))


def average(lines: Sequence[Line]) -> Tuple[np.ndarray, np.ndarray]:
    """AVERAGE: ``y = (1/n) sum y_i`` (e.g. mean ISRS of several input-motion sets)."""
    X, Y = on_union(lines)
    return X, Y.mean(axis=0)


def srss(lines: Sequence[Line]) -> Tuple[np.ndarray, np.ndarray]:
    """SRSS: ``y = sqrt(sum y_i^2)`` (e.g. combination of the X, Y, Z input-direction ISRS)."""
    X, Y = on_union(lines)
    return X, np.sqrt(np.sum(Y * Y, axis=0))


# ======================================================================================
# Envelope, broadening and bridging (D-LIN-01)
# ======================================================================================
def _tol(y: np.ndarray) -> float:
    return 1e-12 * max(1.0, float(np.max(np.abs(y))) if len(y) else 1.0)


def envelope(lines: Sequence[Line], exact: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """Pointwise maximum of the lines (step 1 of BROADEN, test T-B3).

    The values at the union-grid points are ``max_i y_i``.  With ``exact`` the crossing points
    of the lines between grid points are inserted, so the piecewise-linear result *is* the upper
    envelope everywhere (without them the chord would cut the corner where two lines cross).
    """
    X, Y = on_union(lines)
    if not exact or len(lines) == 1 or len(X) < 2:
        return X, Y.max(axis=0)
    tol = _tol(Y)
    for _ in range(64):                      # each pass resolves one switch per interval
        E = Y.max(axis=0)
        dx = np.diff(X)
        S = np.diff(Y, axis=1) / dx          # slopes (nL, nI)
        Yl, Yr = Y[:, :-1], Y[:, 1:]
        # line that is maximal just right of the left end / just left of the right end
        a = _argmax_tiebreak(Yl, S, +1, tol)
        b = _argmax_tiebreak(Yr, S, -1, tol)
        k = np.nonzero(a != b)[0]
        if len(k) == 0:
            break
        ia, ib = a[k], b[k]
        sa, sb = S[ia, k], S[ib, k]
        den = sa - sb
        ok = np.abs(den) > 0
        k, ia, ib, sa, sb, den = k[ok], ia[ok], ib[ok], sa[ok], sb[ok], den[ok]
        # line a: Y[a,k] + sa (x - X_k);  line b: Y[b,k] + sb (x - X_k)
        xs = X[k] + (Y[ib, k] - Y[ia, k]) / den
        inside = (xs > X[k] + 1e-12 * dx[k]) & (xs < X[k + 1] - 1e-12 * dx[k])
        xs = xs[inside]
        if len(xs) == 0:
            break
        Xn = np.unique(np.concatenate([X, xs]))
        Y = np.vstack([L.at(Xn) for L in lines])
        X = Xn
    return X, Y.max(axis=0)


def _argmax_tiebreak(V: np.ndarray, S: np.ndarray, sign: int, tol: float) -> np.ndarray:
    """Per column, the row of the largest value; ties (within ``tol``) broken by ``sign * slope``."""
    vmax = V.max(axis=0)
    cand = V >= vmax - tol
    key = np.where(cand, sign * S, -np.inf)
    return np.argmax(key, axis=0)


class _RangeMax:
    """Sparse table for O(1) range-maximum queries over a fixed array (vectorised)."""

    def __init__(self, a: np.ndarray):
        a = np.asarray(a, dtype=float)
        self.levels = [a]
        k = 1
        while (1 << k) <= len(a):
            prev = self.levels[-1]
            h = 1 << (k - 1)
            self.levels.append(np.maximum(prev[:-h], prev[h:]))
            k += 1

    def query(self, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
        """max a[lo:hi] (half-open); -inf for empty ranges."""
        lo = np.asarray(lo, dtype=np.int64)
        hi = np.asarray(hi, dtype=np.int64)
        out = np.full(lo.shape, -np.inf)
        n = hi - lo
        ok = n > 0
        if np.any(ok):
            k = np.floor(np.log2(n[ok])).astype(np.int64)
            lo_ok, hi_ok = lo[ok], hi[ok]
            vals = np.empty(len(k))
            for lev in np.unique(k):
                m = k == lev
                arr = self.levels[int(lev)]
                vals[m] = np.maximum(arr[lo_ok[m]], arr[hi_ok[m] - (1 << int(lev))])
            out[ok] = vals
        return out


def _window_max(X: np.ndarray, E: np.ndarray, f: np.ndarray, b: float, rm: _RangeMax) -> np.ndarray:
    """``max{E(f') : f/(1+b) <= f' <= f/(1-b)}`` for each f (E piecewise linear on X)."""
    lo = f / (1.0 + b)
    hi = f / (1.0 - b)
    v = np.maximum(resample(X, E, lo), resample(X, E, hi))
    i0 = np.searchsorted(X, lo, side="right")      # first vertex > lo
    i1 = np.searchsorted(X, hi, side="left")       # first vertex >= hi
    return np.maximum(v, rm.query(i0, i1))


def broaden_peaks(X: np.ndarray, E: np.ndarray, b: float, exact: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """Peak broadening by ``+-b`` (step 2 of BROADEN; test T-B1).

    ``exact``: evaluate on the augmented grid ``X U X(1-b) U X(1+b)`` (clipped to the data range)
    and insert the corner points of ``B`` between grid points, so the result is the exact
    broadened function.  Otherwise (``EDUOPT,BROADENGRID,ACS``) the discrete form
    ``B(f_k) = max{E(f_j) : f_j in [f_k/(1+b), f_k/(1-b)]}`` on the source points.
    """
    X = np.asarray(X, dtype=float)
    E = np.asarray(E, dtype=float)
    if not 0.0 <= b < 1.0:
        raise LineError(f"peak broadening must satisfy 0 <= Smooth2 < 100 % (got {100 * b:g} %)")
    if b == 0.0 or len(X) < 2:
        return X.copy(), E.copy()
    if X[0] < 0.0:
        raise LineError("peak broadening needs non-negative abscissas (frequencies)")
    rm = _RangeMax(E)
    if not exact:
        lo, hi = X / (1.0 + b), X / (1.0 - b)
        i0 = np.searchsorted(X, lo * (1 - 1e-12), side="left")
        i1 = np.searchsorted(X, hi * (1 + 1e-12), side="right")
        return X.copy(), rm.query(i0, i1)
    G = np.concatenate([X, X * (1.0 - b), X * (1.0 + b)])
    G = merge_close(G[(G >= X[0]) & (G <= X[-1])], 1e-12)
    B = _window_max(X, E, G, b, rm)
    # corners between grid points: on (p, q) the window ends stay in fixed segments of E and the
    # interior vertices are fixed, so B = max(L(f), H(f), C) with L, H linear and C constant
    p, q = G[:-1], G[1:]
    m = 0.5 * (p + q)
    lo_m, hi_m = m / (1.0 + b), m / (1.0 - b)
    C = rm.query(np.searchsorted(X, lo_m, side="right"), np.searchsorted(X, hi_m, side="left"))
    Lp, Lq = resample(X, E, p / (1.0 + b)), resample(X, E, q / (1.0 + b))
    Hp, Hq = resample(X, E, p / (1.0 - b)), resample(X, E, q / (1.0 - b))
    cands = []
    for (ap, aq), (cp, cq) in (((Lp, Lq), (Hp, Hq)), ((Lp, Lq), (C, C)), ((Hp, Hq), (C, C))):
        d0, d1 = ap - cp, aq - cq
        with np.errstate(invalid="ignore", divide="ignore"):
            t = d0 / (d0 - d1)
        ok = np.isfinite(t) & (t > 1e-12) & (t < 1 - 1e-12) & np.isfinite(cp)
        cands.append(p[ok] + t[ok] * (q[ok] - p[ok]))
    xs = np.concatenate(cands) if cands else np.zeros(0)
    if len(xs):
        Bx = _window_max(X, E, xs, b, rm)
        chord = resample(G, B, xs)
        keep = chord - Bx > _tol(E)
        if np.any(keep):
            G = np.concatenate([G, xs[keep]])
            order = np.argsort(G, kind="stable")
            G = G[order]
            B = np.concatenate([B, Bx[keep]])[order]
            G, uniq = np.unique(G, return_index=True)
            B = B[uniq]
    return G, B


def _peaks(y: np.ndarray, tol: float) -> List[Tuple[int, int]]:
    """Local maxima of ``y`` as index ranges ``(i0, i1)``; a flat plateau counts as one peak.

    A plateau is a peak when its neighbours on both sides are lower; a plateau touching an end of
    the range only needs its inner neighbour to be lower (one-sided maximum).
    """
    n = len(y)
    out: List[Tuple[int, int]] = []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(y[j + 1] - y[i]) <= tol:
            j += 1
        left_ok = i == 0 or y[i - 1] < y[i] - tol
        right_ok = j == n - 1 or y[j + 1] < y[i] - tol
        if left_ok and right_ok and not (i == 0 and j == n - 1):
            out.append((i, j))
        i = j + 1
    return out


def _raise_between(X: np.ndarray, Y: np.ndarray, i0: int, i1: int, level: float, exact: bool,
                   tol: float) -> Tuple[np.ndarray, np.ndarray, bool]:
    """Raise ``Y`` to ``level`` strictly between indices i0 < i1, inserting the crossings (exact).

    A crossing point is inserted only where a segment *clearly* crosses the fill level (one end
    more than ``tol`` above it, the other more than ``tol`` below) and only when it falls strictly
    between the two vertices in floating point.  A vertex within ``tol`` of the level (e.g. a
    broadened plateau vertex 1 ulp below ``A``) is simply raised by the ``np.maximum`` below: the
    crossing would lie on that vertex, and rounding could otherwise put a second point on the
    same abscissa (a Line needs strictly increasing x).  The result changes by at most ``tol``.
    """
    # both ends are peaks (>= level), so only interior points can lie below the fill level
    if not np.any(Y[i0 + 1:i1] < level - tol):
        return X, Y, False
    if exact:
        xs, ys = [], []
        for k in range(i0, i1):
            y0, y1 = Y[k], Y[k + 1]
            if (y0 > level + tol and y1 < level - tol) or (y0 < level - tol and y1 > level + tol):
                t = (level - y0) / (y1 - y0)
                x = X[k] + t * (X[k + 1] - X[k])
                if X[k] < x < X[k + 1]:                  # rounding may land on a vertex: then skip it
                    xs.append(x)
                    ys.append(level)
        if xs:
            Xn = np.concatenate([X, xs])
            Yn = np.concatenate([Y, ys])
            order = np.argsort(Xn, kind="stable")
            X, Y = Xn[order], Yn[order]
            # the crossings lie strictly inside (X[i0], X[i1]): i0 is unchanged, i1 moves
            i1 = i1 + len(xs)
    Y = Y.copy()
    inner = slice(i0 + 1, i1)
    changed = bool(np.any(Y[inner] < level - tol))
    Y[inner] = np.maximum(Y[inner], level)
    return X, Y, changed


def bridge_peaks(X: np.ndarray, Y: np.ndarray, p: float, criterion: str = "VALLEY",
                 exact: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """Peak bridging (step 3 of BROADEN; test T-B2).

    Peaks are found anew in every pass; in a pass every adjacent pair (P_k, P_k+1), from low to
    high frequency, whose valley qualifies is filled up to ``A = min(P_k, P_k+1)``:

    * ``VALLEY`` (default): ``(A - min valley) / A <= p``;
    * ``AMPLITUDE`` (``EDUOPT,BRIDGE,AMPLITUDE``): ``|P_k - P_k+1| / max(P_k, P_k+1) <= p``.

    Passes repeat until nothing changes.  ``exact`` inserts the points where the curve crosses
    the fill level (so the fill starts exactly where the flank reaches ``A``).
    """
    X = np.asarray(X, dtype=float).copy()
    Y = np.asarray(Y, dtype=float).copy()
    if not 0.0 <= p <= 1.0:
        raise LineError(f"peak bridging must satisfy 0 <= Smooth1 <= 100 % (got {100 * p:g} %)")
    criterion = criterion.upper()
    if criterion not in ("VALLEY", "AMPLITUDE"):
        raise LineError(f"unknown bridging criterion {criterion}")
    if p == 0.0 or len(X) < 3:
        return X, Y
    tol = _tol(Y)
    for _ in range(10 * len(X) + 10):
        peaks = _peaks(Y, tol)
        if len(peaks) < 2:
            break
        anchors = [(X[i0], X[i1], Y[i0]) for i0, i1 in peaks]     # positions survive insertions
        changed_pass = False
        for (xa0, xa1, va), (xb0, xb1, vb) in zip(anchors[:-1], anchors[1:]):
            i_end = int(np.searchsorted(X, xa1, side="left"))
            i_start = int(np.searchsorted(X, xb0, side="left"))
            if i_start <= i_end + 1:
                continue
            A = min(va, vb)
            if A <= 0:
                continue
            valley = float(np.min(Y[i_end:i_start + 1]))
            if criterion == "VALLEY":
                ok = (A - valley) / A <= p + 1e-12
            else:
                ok = abs(va - vb) / max(va, vb) <= p + 1e-12
            if not ok:
                continue
            X, Y, ch = _raise_between(X, Y, i_end, i_start, A, exact, tol)
            changed_pass = changed_pass or ch
        if not changed_pass:
            break
    return X, Y


def broaden(lines: Sequence[Line], smooth1: float, smooth2: float, grid: str = "AUGMENTED",
            bridge: str = "VALLEY") -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """BROADEN,Dest,Smooth1,Smooth2,sources (D-LIN-01): envelope -> broaden -> bridge.

    ``smooth1`` = peak bridging percentage ("Peak Difference(%)"), ``smooth2`` = peak broadening
    percentage ("Broaden(%)").  Returns ``(x, y, notes)``; ``notes`` holds warnings (EDU-15 when
    a source has fewer than 301 points).
    """
    if not lines:
        raise LineError("no source lines")
    grid = grid.upper()
    if grid not in ("AUGMENTED", "ACS"):
        raise LineError(f"unknown broadening grid {grid}")
    b = float(smooth2) / 100.0
    p = float(smooth1) / 100.0
    if not 0.0 <= b < 1.0:
        raise LineError(f"Smooth2 (peak broadening) must be in [0, 100) % (got {smooth2:g})")
    if not 0.0 <= p <= 1.0:
        raise LineError(f"Smooth1 (peak bridging) must be in [0, 100] % (got {smooth1:g})")
    notes: List[str] = []
    short = [L.number for L in lines if L.n < BROADEN_MIN_POINTS]
    if short:
        notes.append(f"EDU-15: line(s) {', '.join(str(n) for n in short)} have fewer than "
                     f"{BROADEN_MIN_POINTS} points; a reduced number of frequency steps affects the accuracy "
                     f"of spectrum broadening (manual: use 301 frequencies, 0.1-100 Hz)")
    exact = grid == "AUGMENTED"
    X, E = envelope(lines, exact=exact)
    X, B = broaden_peaks(X, E, b, exact=exact)
    X, B = bridge_peaks(X, B, p, criterion=bridge, exact=exact)
    return X, B, notes


# ======================================================================================
# Text files (D-FIL-02, D-LIN-04, spec 10 section 4.2)
# ======================================================================================
def _parse_numbers(line: str) -> Optional[List[float]]:
    """Numbers of a text line (blank, comma or tab separated, Fortran D exponents), or None when
    the line is not purely numeric (headers and comments are skipped by every reader, D-FIL-02)."""
    s = line.strip()
    if not s or s[0] in "#*!":
        return None
    out: List[float] = []
    for tok in s.replace(",", " ").split():
        try:
            out.append(float(tok.replace("D", "E").replace("d", "e")))
        except ValueError:
            return None
    return out or None


def numeric_rows(path: PathLike) -> List[List[float]]:
    """All numeric rows of a text file (non-numeric lines skipped)."""
    rows: List[List[float]] = []
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    for ln in text.splitlines():
        v = _parse_numbers(ln)
        if v is not None:
            rows.append(v)
    return rows


def read_spec_file(path: PathLike, nlines: int) -> Tuple[np.ndarray, List[np.ndarray]]:
    """READSPEC: column 1 = frequency (not counted), data columns 1..nlines (spec 10 4.2.1).

    Every numeric row must have at least ``1 + nlines`` numbers.  Returns ``(x, [y_1 .. y_n])``
    as read (unsorted; :func:`clean_xy` is applied by the caller).
    """
    if nlines < 1:
        raise LineError("numLines must be >= 1")
    rows = numeric_rows(path)
    if not rows:
        raise LineError(f"{Path(path).name}: no numeric rows")
    short = [k for k, r in enumerate(rows, start=1) if len(r) < 1 + nlines]
    if short:
        have = min(len(r) for r in rows) - 1
        raise LineError(f"{Path(path).name}: numeric row {short[0]} has {len(rows[short[0] - 1])} values; "
                        f"{1 + nlines} needed (the frequency column is not counted in numLines; the file has "
                        f"{max(have, 0)} data column(s))")
    A = np.asarray([r[:1 + nlines] for r in rows], dtype=float)
    return A[:, 0], [A[:, k] for k in range(1, nlines + 1)]


def read_th_file(path: PathLike, pair: int) -> Tuple[np.ndarray, np.ndarray]:
    """READTH: ``pair = 0`` one-column format (first number = dt, then values, several per line
    allowed, ``t_k = (k-1) dt``, D-LIN-04); ``pair = 1`` rows ``t a``."""
    rows = numeric_rows(path)
    name = Path(path).name
    if pair == 1:
        good = [r for r in rows if len(r) >= 2]
        if len(good) < 1:
            raise LineError(f"{name}: no (time, value) pairs")
        A = np.asarray([r[:2] for r in good], dtype=float)
        return A[:, 0], A[:, 1]
    if pair != 0:
        raise LineError("Pair must be 0 (one-column history with dt first) or 1 (time/value pairs)")
    flat = [v for r in rows for v in r]
    if len(flat) < 2:
        raise LineError(f"{name}: a one-column history needs the time step followed by at least one value")
    dt = float(flat[0])
    if not dt > 0:
        raise LineError(f"{name}: time step {dt:g} (first number of the file) must be > 0")
    vals = np.asarray(flat[1:], dtype=float)
    return np.arange(len(vals)) * dt, vals


def fmt_value(v: float) -> str:
    """Exact text of a double (17 significant digits): WRITESPEC/WRITETH round trips are exact."""
    return f"{float(v):24.16E}"


def write_spec_file(path: PathLike, lines: Sequence[Line], header: bool = True) -> Tuple[Path, np.ndarray]:
    """WRITESPEC: one row per union-grid point ``x y_1 ... y_n`` (spec 10 4.2.3).

    A ``#`` header line names the columns (readers skip it, D-FIL-02).  Values are written with 17
    significant digits so READSPEC gives back exactly the resampled values (VP-53).
    """
    if not lines:
        raise LineError("no lines to write")
    if len(lines) > 50:
        raise LineError("at most 50 lines per spectrum file")
    X, Y = on_union(lines)
    path = Path(path)
    out = []
    if header:
        out.append("# x " + " | ".join(f"{L.number}: {L.name}" for L in lines))
    for k in range(len(X)):
        out.append(" ".join([fmt_value(X[k])] + [fmt_value(v) for v in Y[:, k]]))
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return path, X


def write_th_file(path: PathLike, line: Line, rtol: float = 1e-6) -> Tuple[Path, float, List[str]]:
    """WRITETH: one-column history, ``dt = x[1] - x[0]`` on the first line then the values.

    Returns ``(path, dt, notes)``; notes warn about a non-constant step (> 1e-6 dt) or a start
    time other than 0 (the one-column format restarts at t = 0, D-LIN-04).
    """
    if line.n < 2:
        raise LineError(f"line {line.number}: a history needs at least two points to define dt")
    dt = float(line.x[1] - line.x[0])
    notes: List[str] = []
    steps = np.diff(line.x)
    bad = np.abs(steps - dt) > rtol * abs(dt)
    if np.any(bad):
        notes.append(f"line {line.number}: the time step is not constant ({int(bad.sum())} intervals differ from "
                     f"dt = {dt:g} by more than {rtol:g} dt); the values are written with dt = {dt:g}")
    if abs(line.x[0]) > rtol * abs(dt):
        notes.append(f"line {line.number} starts at x = {line.x[0]:g}; the one-column format starts at t = 0")
    path = Path(path)
    out = [fmt_value(dt).strip()] + [fmt_value(v) for v in line.y]
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return path, dt, notes


# ======================================================================================
# CRITFREQ and FRAMESEL algorithms (spec 09 sections 2.4 and 2.12)
# ======================================================================================
def local_maxima(a: np.ndarray) -> np.ndarray:
    """Indices j with ``a[j] >= a[j-1]`` and ``a[j] > a[j+1]`` (interior points; spec 09 2.4 step 3)."""
    a = np.asarray(a, dtype=float)
    if len(a) < 3:
        return np.zeros(0, dtype=np.int64)
    j = np.arange(1, len(a) - 1)
    return j[(a[j] >= a[j - 1]) & (a[j] > a[j + 1])]


@dataclass
class CritPeak:
    freq: float          # TFI peak frequency (Hz)
    number: int          # frequency number round(f / df_TFI)
    tfi: float           # interpolated amplitude at the peak
    ref: float           # max of the computed amplitudes at the bracketing SSI frequencies
    diff_pct: float      # 100 |tfi - ref| / ref
    flagged: bool


def critfreq(f_u: np.ndarray, a_u: np.ndarray, f_i: np.ndarray, a_i: np.ndarray, tol: float,
             minfilter: float) -> Tuple[List[CritPeak], float]:
    """CRITFREQ (spec 09 section 2.4): interpolated TF peaks not supported by computed values.

    ``f_u, a_u``: computed amplitudes at the SSI frequencies (``.TFU``); ``f_i, a_i``: interpolated
    amplitudes at the Fourier frequencies (``.TFI``).  Peaks of the TFI lower than
    ``(1 - minfilter/100) max|TFI|`` are not considered; a considered peak at ``f_p`` is compared
    with ``A_ref = max(A_U(F_m), A_U(F_m+1))`` of the bracketing SSI frequencies and flagged when
    ``100 |A_I(f_p) - A_ref| / A_ref > tol``.  Returns ``(peaks, df_TFI)``.
    """
    f_u = np.asarray(f_u, float)
    a_u = np.abs(np.asarray(a_u, float))
    f_i = np.asarray(f_i, float)
    a_i = np.abs(np.asarray(a_i, float))
    if len(f_u) < 2 or len(f_i) < 3:
        raise LineError("CRITFREQ needs at least 2 computed and 3 interpolated frequencies")
    steps = np.diff(f_i)
    df = float(np.median(steps)) if len(steps) else 0.0
    amax = float(np.max(a_i))
    thresh = (1.0 - float(minfilter) / 100.0) * amax
    out: List[CritPeak] = []
    for j in local_maxima(a_i):
        fp, ap = float(f_i[j]), float(a_i[j])
        if ap < thresh:
            continue
        m = int(np.searchsorted(f_u, fp, side="right")) - 1
        if m < 0 or m >= len(f_u) - 1:
            continue                                  # outside the computed range: nothing to compare
        ref = max(a_u[m], a_u[m + 1])
        if ref <= 0:
            continue
        d = 100.0 * abs(ap - ref) / ref
        number = int(round(fp / df)) if df > 0 else j
        out.append(CritPeak(fp, number, ap, float(ref), float(d), bool(d > float(tol))))
    return out, df


def framesel(a: np.ndarray, tol: float) -> np.ndarray:
    """FRAMESEL (spec 09 section 2.12): 1-based frame numbers of the critical local extrema.

    Local maxima (``a_j >= a_j-1`` and ``a_j >= a_j+1``) and minima (``<=``) whose absolute value
    is at least ``tol/100 max|a|``; the end samples count with their single neighbour; a run of
    equal values counts once (its first sample).  Frame k is time ``(k-1) dt`` (D-STR-08).
    """
    a = np.asarray(a, dtype=float)
    n = len(a)
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    amax = float(np.max(np.abs(a)))
    if n == 1:
        return np.array([1], dtype=np.int64) if amax > 0 else np.zeros(0, dtype=np.int64)
    prev = np.concatenate([[a[1]], a[:-1]])        # one-sided ends: compare with the only neighbour
    nxt = np.concatenate([a[1:], [a[-2]]])
    is_max = (a >= prev) & (a >= nxt)
    is_min = (a <= prev) & (a <= nxt)
    flat_run = np.concatenate([[False], a[1:] == a[:-1]])
    ext = (is_max | is_min) & ~flat_run & (np.abs(a) >= float(tol) / 100.0 * amax)
    if amax == 0:
        ext[:] = False
    return np.nonzero(ext)[0].astype(np.int64) + 1
