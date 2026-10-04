"""Transfer-function interpolation in frequency (MOTION and STRESS).

ANALYS solves the SSI problem only at a few tens to a few hundred *SSI frequencies*
``f_1 < ... < f_N``.  To convolve a transfer function (TF) with a control motion we need
it at every Fourier frequency ``f_k = k df`` of the FFT grid.  This module fills that
gap.  Normative sources: requirements section 4.8 items 2-5, decisions D-MOT-01 ...
D-MOT-05 and D-CNV-12, theory R1 section 5.1-5.2 (prototype ``R1_checks/test5*.py``).

The physics behind options 0-5 (Tajirian 1981, R1 section 5.1)
--------------------------------------------------------------
Over a short frequency band every SSI transfer function looks like the response of a
*two-degree-of-freedom system with hysteretic damping*.  Such a TF is a ratio of two
quadratics in ``w^2`` (no odd powers of ``w`` because hysteretic damping is
frequency independent)::

    H(w) = (C1 w^4 + C2 w^2 + C3) / (w^4 + C4 w^2 + C5)

Five complex constants are fixed by five computed values ``H_p`` (a *window*), from the
linear system ``[w_p^4, w_p^2, 1, -w_p^2 H_p, -H_p] . C = w_p^4 H_p`` (p = 1..5).  We use
the scaled variable ``x = (w / w_max)^2`` of the window (D-MOT-01: "w scaled by the
window maximum"), which leaves the interpolant unchanged and keeps the matrix well
conditioned.  When the window data come from fewer than two modes (an SDOF-like or flat
TF) the system is rank deficient; then *every* member of the solution family
interpolates, and the minimum-norm least-squares solution (SVD, rcond 1e-10) is used
(R1 section 5.1 "Rank deficiency").

The guard (D-MOT-01, R1 section 5.2): when the fitted denominator nearly vanishes on the
window span, ``min|D| < 1e-3 max|D|`` with ``D(x) = x^2 + C4 x + C5``, the window may hold a
*spurious* near-real pole and the complex cubic through the 4 nearest window points is used
instead.  The ratio test alone cannot tell a spurious pole from a genuine one: a lightly
damped mode in a wide window (typical of a uniform, coarse low-frequency grid) also has
``min|D| << max|D|``, and so does the cancelled pole/zero pair of a rank-deficient window.  The
ratio test is therefore only the *trigger*.  The cubic replaces the rational form when, in a
triggered window, a pole of D is none of the following:

1. *cancelled*: a numerator zero z is so close that the doublet factor ``(x - z)/(x - p)``
   differs from 1 by at most :data:`DOUBLET_TOL` anywhere on the span
   (``|p - z| <= DOUBLET_TOL dist(p, span)``);
2. *physical*: under the ``e^{+iwt}`` convention a hysteretic mode ``k c(beta) - m w^2 = 0``
   has its pole at ``x_p = (w0/w_max)^2 c(beta)``.  Since the SASSI factor is
   ``c(beta) = exp(i theta)`` with ``beta = sin(theta/2)`` (D-CNV-03), a pole with
   ``0 < arg x_p < pi`` is a damped mode with implied damping ``beta_p = sin(arg x_p / 2)``.
   Poles with ``beta_p >= DAMPING_FLOOR`` (1e-4) are accepted as genuine resonances;
3. *remote*: its own factor ratio ``dist(p, span) / max|x - p|`` is at least
   ``sqrt(GUARD_RATIO)``.  Since ``ratio(D) >= ratio(x - p1) ratio(x - p2)``, two remote poles
   can never trip the trigger: in a triggered window at least one pole is near the span, and a
   remote pole only adds a smooth factor varying by less than ``1/sqrt(GUARD_RATIO)`` ~ 32.

The cubic is also used when a fit does not reproduce its own 5 points (node error above
:data:`NODE_FAIL_TOL`).  Exact 2-DOF and SDOF hysteretic data therefore keep the exact rational
form (VP-28), while a near-real pole with negative or vanishing damping that no zero cancels
(a spurious spike from noisy or non-2-DOF data) is still caught.

Window schemes (requirements section 4.8 table; 1-based start index m, window m..m+4,
candidate windows of interval j: ``max(1, j-3) <= m <= min(j, N-4)``):

====== ============================================== ==================================
Option Windows                                        Combination
====== ============================================== ==================================
0      all candidate windows (SASSI2000)              weighted, w = max(1e-3, 1-|f-c|/h)
1      starts 1, 5, 9, ... (SASSI 1982); tail N-4     single window
2      all candidate windows                          arithmetic mean
3      the 3 candidates with the smallest |f - c_m|   arithmetic mean
4      starts 2, 6, 10, ...; interval 1 uses m = 1    single window
5      starts 3, 7, 11, ...; intervals 1-2 use m = 1  single window
6      no windows: not-a-knot cubic splines of Re H and Im H through the computed points
       plus the f = 0 anchor (D-MOT-02)
====== ============================================== ==================================

``N < 5`` turns options 0-5 into option 6, ``N < 3`` into complex linear interpolation.
Every scheme returns ``H_j`` exactly at the computed frequencies (enforced bit-exactly).

Outside the computed range (D-MOT-03): ``f > f_N`` gives 0 (the SSI cut-off acts as a
low-pass filter); below ``f_1`` a seismic TF is linear between the rigid-body anchor
``H(0) = h0`` and ``H_1``, a vibration TF is held constant at ``H_1``.

Post-processing: smoothing (D-MOT-04, options 0-5 only) and phase adjustment (D-MOT-05).
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
from scipy.interpolate import CubicSpline

__all__ = [
    "OPTION_NAMES", "RCOND", "RCOND_LADDER", "NODE_TOL", "NODE_FAIL_TOL", "GUARD_RATIO",
    "DOUBLET_TOL", "DAMPING_FLOOR", "WEIGHT_FLOOR",
    "interpolate_tf", "effective_option", "window_table", "fit_windows", "window_poles",
    "smooth_tf", "phase_adjust", "phase_factor", "rigid_body_anchor",
]

#: Interpolation option names (MOTION dialog / requirements section 4.8; spec 05c B.5.5).
OPTION_NAMES: Dict[int, str] = {
    0: "SASSI2000 dense overlapping windows, weighted averaging",
    1: "original SASSI 1982 non-overlapping windows",
    2: "dense overlapping windows, averaging",
    3: "only three overlapping windows, averaging",
    4: "non-overlapping windows with one position shift",
    5: "non-overlapping windows with two position shift",
    6: "complex bicubic spline (not-a-knot splines of Re and Im)",
}

RCOND = 1e-10          # relative singular-value cut of the 5x5 window solve (D-MOT-01)
RCOND_LADDER = (1e-12, 1e-14, 1e-16)   # lower cuts tried when the 1e-10 solution does not interpolate
NODE_TOL = 1e-12       # accepted relative error of a window fit at its own 5 points
NODE_FAIL_TOL = 1e-6   # a fit worse than this at its own 5 points is not used (guard -> cubic)
GUARD_RATIO = 1e-3     # min|D| < GUARD_RATIO * max|D| on the window span triggers the guard (D-MOT-01)
DOUBLET_TOL = 1e-3     # pole/zero pair cancelled when |p - z| <= DOUBLET_TOL * dist(p, span)
DAMPING_FLOOR = 1e-4   # poles with implied damping sin(arg x_p / 2) >= this are genuine modes
WEIGHT_FLOOR = 1e-3    # option 0 weight floor epsilon (requirements section 4.8 table)
_MATCH_RTOL = 1e-9     # an output frequency within 1e-9 f_N of f_j *is* f_j
_WORK_BUDGET = 4_000_000   # complex elements per working array (column chunking)


# ---------------------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------------------
def effective_option(option: int, n_ssi: int) -> int:
    """Scheme actually used for ``n_ssi`` computed points (requirements section 4.8 item 2).

    Returns the option itself, 6 when options 0-5 have fewer than 5 points, or -1 for complex
    linear interpolation (fewer than 3 points).
    """
    option = int(option)
    if option not in OPTION_NAMES:
        raise ValueError(f"interpolation option must be 0..6, got {option}")
    if n_ssi < 3:
        return -1
    if n_ssi < 5 and option != 6:
        return 6
    return option


def phase_factor(option: int, smooth: float) -> float:
    """Residual phase factor rho of the phase adjustment (D-MOT-05, spec 05c B.6.8).

    ``rho = 1/(1+S)`` for the window schemes 0-5 (S = smoothing parameter, which acts on these
    schemes only); ``rho = 0`` (complete zeroing of the phase) for option 6 and for the complex
    linear fallback ``-1`` of :func:`effective_option`, on which S has no effect.  Pass the
    *effective* option (the scheme actually used): with fewer than 5 SSI frequencies options
    0-5 run the spline and therefore get its rho = 0.
    """
    return 1.0 / (1.0 + float(smooth)) if int(option) in (0, 1, 2, 3, 4, 5) else 0.0


def rigid_body_anchor(dof: int, cm: int = 0, ang_deg: float = 0.0) -> float:
    """Zero-frequency value H(0) of a seismic TF for global output DOF ``dof`` (D-MOT-03).

    As w -> 0 the structure moves rigidly with the free field, so the total-motion TF tends to
    the projection of the control-motion direction on the output DOF.  The control direction
    ``cm`` (SITE <cm>: 0 x', 1 y', 2 z') is given in the SITE axes x'y'z', which are rotated by
    the ANALYS angle ``ang`` about Z: x' = (cos a, sin a, 0), y' = (-sin a, cos a, 0),
    z' = (0, 0, 1) in global axes.  Rotational DOFs (4..6) have H(0) = 0.
    """
    if dof not in (1, 2, 3):
        return 0.0
    a = np.deg2rad(float(ang_deg))
    e = {0: (np.cos(a), np.sin(a), 0.0), 1: (-np.sin(a), np.cos(a), 0.0), 2: (0.0, 0.0, 1.0)}
    if int(cm) not in e:
        raise ValueError(f"control direction cm must be 0, 1 or 2, got {cm}")
    v = float(e[int(cm)][dof - 1])
    return 0.0 if abs(v) < 1e-15 else v


def window_table(option: int, n_ssi: int, f_ssi: np.ndarray) -> Tuple[np.ndarray, Optional[str]]:
    """Candidate windows of every interval for options 0-5.

    Returns ``(win, rule)``: ``win`` is an int array (N-1, 4) of 0-based window starts per
    interval j (0-based interval [f_j, f_{j+1}]), -1 where unused; ``rule`` is ``'single'``,
    ``'mean'``, ``'weighted'`` or ``'nearest3'`` (the latter two need the target frequency
    and are resolved in :func:`_combine`).
    """
    N = int(n_ssi)
    if N < 5:
        raise ValueError("window schemes need at least 5 computed frequencies")
    nint = N - 1
    win = -np.ones((nint, 4), dtype=np.int64)
    j1 = np.arange(1, nint + 1)                       # 1-based interval numbers
    if option in (0, 2, 3):
        lo = np.maximum(1, j1 - 3)
        hi = np.minimum(j1, N - 4)
        for c in range(4):
            m = lo + c
            ok = m <= hi
            win[ok, c] = m[ok] - 1
        rule = {0: "weighted", 2: "mean", 3: "nearest3"}[option]
        return win, rule
    # single-window partitions: starts s0, s0+4, ...; leading intervals before s0 use m = 1
    s0 = {1: 1, 4: 2, 5: 3}[option]
    m = np.where(j1 < s0, 1, s0 + 4 * ((j1 - s0) // 4))
    m = np.minimum(m, N - 4)                          # tail window = last 5 points
    win[:, 0] = m - 1
    return win, "single"


# ---------------------------------------------------------------------------------------
# 5-point rational windows (Tajirian form)
# ---------------------------------------------------------------------------------------
def fit_windows(f_ssi: np.ndarray, H: np.ndarray, rcond: Optional[float] = None,
                adaptive: bool = True) -> Dict[str, np.ndarray]:
    """Fit the 2-DOF hysteretic rational form on every 5-point window (D-MOT-01, R1 section 5.1).

    ``H`` has shape (N, m).  Returns a dict with ``C`` (W, m, 5) complex constants in the scaled
    variable ``x = (w/w_max)^2``, ``scale`` (W,) = w_max of each window, ``x0`` (W,) the scaled
    lowest window frequency, ``rank`` (W, m), ``rcond`` (W, m) the singular-value cut finally
    used, ``node_err`` (W, m) the relative interpolation error at the 5 window points,
    ``trigger`` (W, m) the literal D-MOT-01 ratio test ``min|D| < 1e-3 max|D|`` and ``guard``
    (W, m) booleans (True = use the cubic fallback in that window; see the module docstring:
    the trigger refined by doublet cancellation and the physical-pole test).

    Rank handling: the minimum-norm solution with the prescribed cut ``rcond`` (default
    :data:`RCOND` = 1e-10, D-MOT-01) is accepted when it reproduces the 5 window values to
    :data:`NODE_TOL`.  That is the case for genuinely rank-deficient (SDOF-like or flat) data,
    whose dropped singular value is at round-off level.  For smooth 2-DOF data in a narrow
    window the smallest singular value can be small (1e-11 relative) but *informative*;
    truncating it destroys the interpolation property (errors ~1e-8).  With ``adaptive`` the cut
    is then lowered along :data:`RCOND_LADDER` until the window values are reproduced, which keeps
    options 0-5 exact for 2-DOF data (VP-28) without affecting the rank-deficient cases.
    """
    if rcond is None:
        rcond = RCOND
    f_ssi = np.asarray(f_ssi, dtype=float)
    H = np.asarray(H, dtype=complex)
    N, m = H.shape
    W = N - 4
    idx = np.arange(W)[:, None] + np.arange(5)[None, :]           # (W, 5)
    w = 2.0 * np.pi * f_ssi
    scale = w[idx[:, -1]]                                          # (W,)
    x = (w[idx] / scale[:, None]) ** 2                             # (W, 5) real, x[:, -1] = 1
    Hw = np.transpose(H[idx], (0, 2, 1))                           # (W, m, 5)
    xb = np.broadcast_to(x[:, None, :], Hw.shape)
    A = np.empty(Hw.shape + (5,), dtype=complex)                   # (W, m, 5 rows, 5 cols)
    A[..., 0] = xb * xb
    A[..., 1] = xb
    A[..., 2] = 1.0
    A[..., 3] = -xb * Hw
    A[..., 4] = -Hw
    rhs = xb * xb * Hw                                             # (W, m, 5)
    # minimum-norm least squares by SVD with a relative singular-value cut (numpy lstsq rule)
    u, sv, vh = np.linalg.svd(A)
    uhb = np.einsum("...pi,...p->...i", u.conj(), rhs)
    hmax = np.max(np.abs(Hw), axis=-1)
    hscale = np.where(hmax > 0, hmax, 1.0)

    def solve(rc: float):
        keep = sv > rc * sv[..., :1]
        inv = np.where(keep, 1.0 / np.where(keep, sv, 1.0), 0.0)
        C = np.einsum("...ij,...i->...j", vh.conj(), inv * uhb)    # (W, m, 5)
        # node error: row p of A C - b is N(x_p) - H_p D(x_p), so the value error is that / D(x_p)
        num = (C[..., 0:1] * xb + C[..., 1:2]) * xb + C[..., 2:3]
        den = (xb + C[..., 3:4]) * xb + C[..., 4:5]
        with np.errstate(divide="ignore", invalid="ignore"):
            e = np.abs(num - Hw * den) / np.abs(den)
        e = np.where(np.isfinite(e), e, np.inf).max(axis=-1) / hscale
        return C, keep.sum(axis=-1), e

    C, rank, err = solve(rcond)
    used = np.full(rank.shape, rcond)
    if adaptive:
        for rc in RCOND_LADDER:
            if rc >= rcond:
                continue
            bad = err > NODE_TOL
            if not np.any(bad):
                break
            C2, rank2, err2 = solve(rc)
            better = bad & (err2 < err)
            C = np.where(better[..., None], C2, C)
            rank = np.where(better, rank2, rank)
            err = np.where(better, err2, err)
            used = np.where(better, rc, used)
    x0 = x[:, 0]
    trigger = _denominator_ratio(C[..., 3], C[..., 4], x0) < GUARD_RATIO
    guard = _spurious_pole_guard(C, x0, trigger) | ~(err <= NODE_FAIL_TOL)
    return {"C": C, "scale": scale, "x0": x0, "rank": rank, "rcond": used, "node_err": err,
            "trigger": trigger, "guard": guard}


def _span_lo(x0: np.ndarray, shape: Tuple[int, ...]) -> np.ndarray:
    """Broadcast the per-window lower span end ``x0`` (W,) to ``shape`` = (W, m, ...)."""
    return np.broadcast_to(x0.reshape(x0.shape + (1,) * (len(shape) - 1)), shape)


def _denominator_ratio(c4: np.ndarray, c5: np.ndarray, x0: np.ndarray) -> np.ndarray:
    """``min|D| / max|D|`` over the window span [x0, 1] for D(x) = x^2 + c4 x + c5 (D-MOT-01).

    The extrema of |D(x)|^2 on real x are at the interval ends and at the real roots of the
    cubic ``2x^3 + 3p x^2 + (p^2 + 2r + q^2) x + (p r + q s) = 0`` (c4 = p + iq, c5 = r + is),
    so the minimum and maximum over [x0, 1] are found exactly (batched companion eigenvalues).
    """
    p, q = c4.real, c4.imag
    r, s = c5.real, c5.imag
    shape = p.shape
    comp = np.zeros(shape + (3, 3))
    comp[..., 0, 0] = -1.5 * p
    comp[..., 0, 1] = -0.5 * (p * p + 2.0 * r + q * q)
    comp[..., 0, 2] = -0.5 * (p * r + q * s)
    comp[..., 1, 0] = 1.0
    comp[..., 2, 1] = 1.0
    roots = np.linalg.eigvals(comp)                                # (..., 3) complex
    lo = _span_lo(x0, shape)
    real_ok = np.abs(roots.imag) <= 1e-9 * (1.0 + np.abs(roots.real))
    xr = roots.real
    inside = real_ok & (xr > lo[..., None]) & (xr < 1.0)
    cand = np.concatenate([lo[..., None], np.ones(shape + (1,)), np.where(inside, xr, lo[..., None])], axis=-1)
    D = np.abs(cand * cand + c4[..., None] * cand + c5[..., None])
    dmax = D.max(axis=-1)
    dmin = D.min(axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(dmax > 0, dmin / np.where(dmax > 0, dmax, 1.0), 0.0)


def _quadratic_roots(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Both roots of ``a x^2 + b x + c`` (complex, batched, cancellation-free formula).

    A missing root (degree < 2) is returned as ``inf``; an identically zero polynomial gives
    ``nan`` roots (it has no defined zeros).
    """
    disc = np.sqrt(b * b - 4.0 * a * c + 0j)
    disc = np.where((b.conj() * disc).real < 0, -disc, disc)     # |b + disc| >= |b - disc|
    q = -0.5 * (b + disc)
    with np.errstate(divide="ignore", invalid="ignore"):
        r1 = np.where(a != 0, q / np.where(a != 0, a, 1.0), np.inf)
        r2 = np.where(q != 0, c / np.where(q != 0, q, 1.0), np.where(a != 0, 0.0, np.inf))
    dead = (a == 0) & (b == 0) & (c == 0)
    return np.where(dead, np.nan, r1), np.where(dead, np.nan, r2)


def window_poles(C: np.ndarray) -> Dict[str, np.ndarray]:
    """Poles, zeros and implied damping of fitted windows (``C`` (..., 5), scaled x = (w/w_max)^2).

    ``poles`` (..., 2) are the roots of ``x^2 + C4 x + C5``; ``zeros`` (..., 2) the roots of
    ``C1 x^2 + C2 x + C3``; ``damping`` (..., 2) the implied hysteretic damping
    ``sin(arg x_p / 2)`` of each pole (negative for a pole below the real axis, i.e. a mode with
    negative damping under the e^{+iwt} convention).  A diagnostic for the guard and for teaching:
    the pole of a hysteretic SDOF is ``(w0/w_max)^2 c(beta)`` and its implied damping is beta.
    """
    one = np.ones(C.shape[:-1], dtype=complex)
    p1, p2 = _quadratic_roots(one, C[..., 3], C[..., 4])
    z1, z2 = _quadratic_roots(C[..., 0], C[..., 1], C[..., 2])
    poles = np.stack([p1, p2], axis=-1)
    return {"poles": poles, "zeros": np.stack([z1, z2], axis=-1), "damping": np.sin(0.5 * np.angle(poles))}


def _spurious_pole_guard(C: np.ndarray, x0: np.ndarray, trigger: np.ndarray) -> np.ndarray:
    """Refined D-MOT-01 guard: True where the rational form must be replaced by the cubic.

    Applied where the ratio test ``trigger`` fires (``min|D| < GUARD_RATIO max|D|`` on [x0, 1]).
    Write ``D = (x - p1)(x - p2)``.  A pole p is *harmless* when

    * it is **cancelled**: a numerator zero z (each zero used once) lies so close that the
      doublet factor ``(x - z)/(x - p)`` differs from 1 by at most ``DOUBLET_TOL`` on the span,
      ``|p - z| <= DOUBLET_TOL dist(p, [x0, 1])`` -- the pole/zero pair a rank-deficient (SDOF or
      flat) window carries;
    * it is **physical**: a damped mode, implied damping ``sin(arg p / 2) >= DAMPING_FLOOR``
      (a genuine, possibly sharp resonance between computed points -- the Tajirian premise);
    * it is **remote**: its own factor ratio ``dist(p) / max|x - p|`` on the span is at least
      ``sqrt(GUARD_RATIO)``.  Because ``ratio(D) >= ratio(x - p1) ratio(x - p2)``, a triggered
      window always has a pole below that bound, and only such a near pole can make a spike.

    The guard fires when a triggered window holds a pole that is none of these: a near-real
    pole with negative or vanishing damping that no zero cancels (a spurious spike).
    """
    guard = np.zeros(trigger.shape, dtype=bool)
    if not np.any(trigger):
        return guard
    lo = _span_lo(x0, trigger.shape)[trigger]                      # (nT,)
    Ct = C[trigger]                                                # (nT, 5)
    pz = window_poles(Ct)
    P, Z, damp = pz["poles"], pz["zeros"], pz["damping"]          # (nT, 2)
    # distance of each pole from the real segment [x0, 1] and its largest distance on it
    dist = np.abs(P - np.clip(P.real, lo[:, None], 1.0))
    reach = np.maximum(np.abs(P - lo[:, None]), np.abs(P - 1.0))
    # doublet cancellation: the two poles are matched to distinct zeros (try both pairings)
    with np.errstate(invalid="ignore"):
        close = np.abs(P[:, :, None] - Z[:, None, :]) <= DOUBLET_TOL * dist[:, :, None]   # (nT, pole, zero)
    close &= ~np.isnan(np.abs(Z))[:, None, :]
    pair_a = np.stack([close[:, 0, 0], close[:, 1, 1]], axis=1)    # p1<->z1, p2<->z2
    pair_b = np.stack([close[:, 0, 1], close[:, 1, 0]], axis=1)    # p1<->z2, p2<->z1
    use_b = pair_b.sum(axis=1) > pair_a.sum(axis=1)
    cancelled = np.where(use_b[:, None], pair_b, pair_a)
    cancelled |= np.all(Ct[:, :3] == 0, axis=1)[:, None]          # N == 0: the window is exactly 0
    physical = damp >= DAMPING_FLOOR
    remote = dist >= np.sqrt(GUARD_RATIO) * reach
    guard[trigger] = np.any(~(cancelled | physical | remote), axis=1)
    return guard


def _rational_eval(fit: Dict[str, np.ndarray], wstart: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Evaluate window ``wstart[t]`` at target frequency ``f[t]`` for every column -> (nT, m)."""
    C = fit["C"][wstart]                                           # (nT, m, 5)
    x = ((2.0 * np.pi * f) / fit["scale"][wstart]) ** 2            # (nT,)
    x = x[:, None]
    num = (C[..., 0] * x + C[..., 1]) * x + C[..., 2]
    den = (x + C[..., 3]) * x + C[..., 4]
    return num / den


def _lagrange4(f_ssi: np.ndarray, H: np.ndarray, wstart: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Complex cubic through the 4 window points nearest to each target (guard fallback).

    Of the 5 sorted window points the 4 nearest to f are contiguous: drop the farther end.
    """
    first = f_ssi[wstart]
    last = f_ssi[wstart + 4]
    s = np.where(np.abs(f - first) <= np.abs(f - last), wstart, wstart + 1)   # (nT,)
    pts = s[:, None] + np.arange(4)[None, :]                       # (nT, 4)
    xs = f_ssi[pts]                                                # (nT, 4)
    L = np.ones((len(f), 4))
    for a in range(4):
        for b in range(4):
            if a != b:
                L[:, a] *= (f - xs[:, b]) / (xs[:, a] - xs[:, b])
    return np.einsum("ta,tam->tm", L, H[pts])


def _window_interp(f_ssi: np.ndarray, H: np.ndarray, f: np.ndarray, option: int,
                   fit: Optional[Dict[str, np.ndarray]] = None) -> np.ndarray:
    """Options 0-5 at interior targets ``f`` (f_1 <= f <= f_N) -> (nT, m)."""
    N = len(f_ssi)
    if fit is None:
        fit = fit_windows(f_ssi, H)
    win, rule = window_table(option, N, f_ssi)
    j = np.clip(np.searchsorted(f_ssi, f, side="right") - 1, 0, N - 2)   # interval index
    cand = win[j]                                                  # (nT, 4)
    valid = cand >= 0
    cs = np.where(valid, cand, 0)
    centre = 0.5 * (f_ssi[cs] + f_ssi[cs + 4])
    half = 0.5 * (f_ssi[cs + 4] - f_ssi[cs])
    dist = np.abs(f[:, None] - centre)
    if rule == "single":
        wts = valid.astype(float)
    elif rule == "mean":
        wts = valid.astype(float)
    elif rule == "weighted":
        wts = np.where(valid, np.maximum(WEIGHT_FLOOR, 1.0 - dist / half), 0.0)
    else:  # nearest3: smallest |f - c_m|, ties -> smaller m (stable sort on the column order)
        key = np.where(valid, dist, np.inf)
        order = np.argsort(key, axis=1, kind="stable")
        rank = np.empty_like(order)
        np.put_along_axis(rank, order, np.arange(4)[None, :].repeat(len(f), 0), axis=1)
        wts = (valid & (rank < 3)).astype(float)
    m = H.shape[1]
    acc = np.zeros((len(f), m), dtype=complex)
    for c in range(4):
        sel = wts[:, c] > 0
        if not np.any(sel):
            continue
        ws = cs[sel, c]
        vals = _rational_eval(fit, ws, f[sel])
        g = fit["guard"][ws]                                       # (nS, m)
        if np.any(g):
            lag = _lagrange4(f_ssi, H, ws, f[sel])
            vals = np.where(g, lag, vals)
        acc[sel] += wts[sel, c][:, None] * vals
    return acc / wts.sum(axis=1)[:, None]


# ---------------------------------------------------------------------------------------
# Option 6 and linear
# ---------------------------------------------------------------------------------------
def _spline_interp(f_ssi: np.ndarray, H: np.ndarray, f: np.ndarray, anchor: Optional[np.ndarray]) -> np.ndarray:
    """Option 6: separate not-a-knot cubic splines of Re H and Im H (D-MOT-02).

    The f = 0 anchor (rigid-body value for seismic TFs, H_1 for vibration) is a spline knot
    when f_1 > 0, so the low-frequency end of the spline is physically constrained.
    """
    if anchor is not None and f_ssi[0] > 0.0:
        x = np.concatenate([[0.0], f_ssi])
        y = np.concatenate([np.asarray(anchor, dtype=complex)[None, :], H], axis=0)
    else:
        x, y = f_ssi, H
    sr = CubicSpline(x, y.real, axis=0, bc_type="not-a-knot")
    si = CubicSpline(x, y.imag, axis=0, bc_type="not-a-knot")
    return sr(f) + 1j * si(f)


def _linear_interp(f_ssi: np.ndarray, H: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Complex linear interpolation between computed points (N < 3 fallback)."""
    if len(f_ssi) == 1:
        return np.repeat(H[:1], len(f), axis=0)
    out = np.empty((len(f), H.shape[1]), dtype=complex)
    for c in range(H.shape[1]):
        out[:, c] = np.interp(f, f_ssi, H[:, c].real) + 1j * np.interp(f, f_ssi, H[:, c].imag)
    return out


# ---------------------------------------------------------------------------------------
# Smoothing and phase adjustment
# ---------------------------------------------------------------------------------------
def smooth_tf(f_ssi: np.ndarray, H: np.ndarray, f: np.ndarray, Hi: np.ndarray, S: float) -> np.ndarray:
    """Band filter of interpolated values between computed points (D-MOT-04, spec 05c B.6.7).

    For f in [f_j, f_{j+1}] with the complex linear reference L(f), the excursion of |H(f)|
    outside the band [A_min, A_max] = [min, max](|H_j|, |H_{j+1}|) is measured by
    ``r = max(0, |H|/A_max - 1) + max(0, 1 - |H|/max(A_min, 1e-12 A_max))`` and pulled back:
    ``H_s = L + (H - L)/(1 + S r)``.  S = 0 is the identity; values inside the band and at the
    computed points are unchanged.  Genuine resonant peaks between computed points are clipped,
    which is why S must be 0 for coherent analyses.
    """
    if S <= 0:
        return Hi
    N = len(f_ssi)
    j = np.clip(np.searchsorted(f_ssi, f, side="right") - 1, 0, N - 2)
    fj, fj1 = f_ssi[j], f_ssi[j + 1]
    t = ((f - fj) / (fj1 - fj))[:, None]
    Hj, Hj1 = H[j], H[j + 1]
    L = Hj + (Hj1 - Hj) * t
    a1, a2 = np.abs(Hj), np.abs(Hj1)
    amax = np.maximum(a1, a2)
    amin = np.minimum(a1, a2)
    amag = np.abs(Hi)
    with np.errstate(divide="ignore", invalid="ignore"):
        over = np.where(amax > 0, np.maximum(0.0, amag / np.where(amax > 0, amax, 1.0) - 1.0), 0.0)
        floor = np.maximum(amin, 1e-12 * amax)
        under = np.where(floor > 0, np.maximum(0.0, 1.0 - amag / np.where(floor > 0, floor, 1.0)), 0.0)
    r = over + under
    # both neighbours exactly zero: any excursion is unsupported -> fall back to L (= 0)
    dead = (amax == 0) & (amag > 0)
    out = L + (Hi - L) / (1.0 + S * r)
    return np.where(dead, L, out)


def phase_adjust(H: np.ndarray, rho: float) -> np.ndarray:
    """Phase adjustment referred to the zero-frequency phase (D-MOT-05, spec 05c B.6.8)::

        phi = unwrap(arg H)  (along the grid),   phi_0 = phase at f = 0
        H'  = |H| exp(i (phi_0 + rho (phi - phi_0)))

    ``H`` (nK,) or (nK, m) must be ordered by ascending frequency starting at f = 0.  Every
    differential phase between neighbouring Fourier components is scaled by rho, as B.6.8
    requires, while the zero-frequency value -- the rigid-body anchor of a seismic TF (D-MOT-03)
    -- is kept: rho = 0 gives ``|H|`` times the sign of the anchor.  For a positive anchor
    (phi_0 = 0) this is the B.6.8 formula ``|H| exp(i rho phi)`` and rho = 0 gives exactly |H|.
    Referring the phase to phi_0 makes the operation sign-equivariant, ``PA(-H) = -PA(H)``:
    reversing the control direction (SITE angle + 180 deg, anchor -1) reverses the response
    instead of turning -|H| into +|H|.  When H(0) = 0 (rotations, components normal to the
    input) phi_0 is the phase of the lowest-frequency non-zero value (the limit f -> 0+).
    rho = 1 returns H unchanged.
    """
    H = np.asarray(H, dtype=complex)
    rho = float(rho)
    if rho == 1.0:
        return H.copy()
    one_d = H.ndim == 1
    H2 = H.reshape(-1, 1) if one_d else H
    amp = np.abs(H2)
    phi = np.unwrap(np.angle(H2), axis=0)
    first = np.argmax(amp > 0, axis=0)                           # lowest non-zero row (0 if all zero)
    phi0 = np.take_along_axis(phi, first[None, :], axis=0)       # (1, m)
    out = amp * np.exp(1j * (phi0 + rho * (phi - phi0)))
    return out[:, 0] if one_d else out


# ---------------------------------------------------------------------------------------
# Public driver (binding signature, ARCHITECTURE section 6.4)
# ---------------------------------------------------------------------------------------
def interpolate_tf(f_ssi, H, f_out, option: int, smooth: float = 0.0, pzadj: int = 0,
                   h0=None, mode: str = "seismic") -> np.ndarray:
    """Interpolate computed TFs ``H`` (nF, m) at SSI frequencies ``f_ssi`` onto ``f_out`` (nK,).

    Implements requirements section 4.8 items 2-5 (D-MOT-01 ... D-MOT-05):

    * interior ``f_1 <= f <= f_N``: scheme ``option`` (0-6, see the module docstring), then
      smoothing with parameter ``smooth`` (options 0-5; ignored for option 6);
    * values at the computed frequencies are returned *exactly* (TFI = TFU, VP-28);
    * ``f > f_N``: 0;  ``0 <= f < f_1``: complex linear from ``h0`` (m,) to H_1 for
      ``mode='seismic'`` (h0 defaults to 0), constant H_1 for ``mode='vibration'``;
    * ``pzadj = 1``: phase adjustment over the whole output grid (:func:`phase_adjust`), the
      phase referred to the lowest output frequency (normally f = 0), with rho of the scheme
      actually used (:func:`phase_factor` of :func:`effective_option`): 1/(1+S) for options
      0-5, 0 for option 6 and for the N < 5 / N < 3 fallbacks, on which S has no effect.

    ``H`` may be 1-D (one TF); the result then is 1-D.  ``f_ssi`` must be strictly increasing
    and positive.  The Nyquist-bin rule of D-MOT-03 is applied by the caller, which knows the
    FFT grid (see :func:`sassi.core.signal.hermitian_bins`).
    """
    f_ssi = np.asarray(f_ssi, dtype=float).ravel()
    Hin = np.asarray(H, dtype=complex)
    one_d = Hin.ndim == 1
    H2 = Hin.reshape(-1, 1) if one_d else Hin
    if H2.shape[0] != len(f_ssi):
        raise ValueError(f"H has {H2.shape[0]} rows for {len(f_ssi)} SSI frequencies")
    if len(f_ssi) == 0:
        raise ValueError("no SSI frequencies")
    if np.any(np.diff(f_ssi) <= 0):
        raise ValueError("SSI frequencies must be strictly increasing (duplicates are an error)")
    if f_ssi[0] < 0:
        raise ValueError("SSI frequencies must be non-negative")
    mode = mode.lower()
    if mode not in ("seismic", "vibration"):
        raise ValueError("mode must be 'seismic' or 'vibration'")
    eff = effective_option(option, len(f_ssi))
    smooth = float(smooth or 0.0)
    if smooth < 0:
        raise ValueError("smoothing parameter must be >= 0")
    m = H2.shape[1]
    if mode == "seismic":
        anchor = np.zeros(m, dtype=complex) if h0 is None else np.broadcast_to(np.asarray(h0, dtype=complex), (m,)).copy()
    else:
        anchor = H2[0].copy()

    f_out = np.asarray(f_out, dtype=float).ravel()
    order = np.argsort(f_out, kind="stable")
    fo = f_out[order]
    out = np.zeros((len(fo), m), dtype=complex)

    f1, fN = f_ssi[0], f_ssi[-1]
    tol = _MATCH_RTOL * max(fN, 1e-300)
    inside = (fo >= f1 - tol) & (fo <= fN + tol)
    fi = np.clip(fo[inside], f1, fN)
    if fi.size:
        # column chunks keep the working arrays small (D-GEN-10)
        nW = max(len(f_ssi) - 4, 1)
        per_col = max(nW * 25, len(fi) * 5, 1)
        step = max(1, _WORK_BUDGET // per_col)
        vals = np.empty((len(fi), m), dtype=complex)
        for c0 in range(0, m, step):
            cols = slice(c0, min(m, c0 + step))
            Hc = H2[:, cols]
            if eff in (0, 1, 2, 3, 4, 5):
                v = _window_interp(f_ssi, Hc, fi, eff)
                v = smooth_tf(f_ssi, Hc, fi, v, smooth)
            elif eff == 6:
                v = _spline_interp(f_ssi, Hc, fi, anchor[cols])
            else:
                v = _linear_interp(f_ssi, Hc, fi)
            vals[:, cols] = v
        # exactness at the computed frequencies (bit-exact TFI = TFU)
        k = np.clip(np.searchsorted(f_ssi, fi), 0, len(f_ssi) - 1)
        km = np.clip(k - 1, 0, len(f_ssi) - 1)
        near = np.where(np.abs(f_ssi[k] - fi) <= np.abs(f_ssi[km] - fi), k, km)
        hit = np.abs(f_ssi[near] - fi) <= tol
        vals[hit] = H2[near[hit]]
        out[inside] = vals
    below = fo < f1 - tol
    if np.any(below):
        if mode == "seismic":
            t = (fo[below] / f1)[:, None] if f1 > 0 else np.zeros((int(below.sum()), 1))
            out[below] = anchor[None, :] + (H2[0] - anchor)[None, :] * t
        else:
            out[below] = H2[0][None, :]
    if pzadj:
        out = phase_adjust(out, phase_factor(eff, smooth))
    res = np.empty_like(out)
    res[order] = out
    return res[:, 0] if one_d else res
