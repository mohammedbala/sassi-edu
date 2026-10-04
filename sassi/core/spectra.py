"""Response spectra and ground-motion parameters shared by MOTION, SOIL, EQUAKE and the UI.

Response spectra use the exact piecewise-linear recurrence of Nigam & Jennings (1969)
for the SDOF equation (requirements section 4.8 item 8, D-MOT-07)::

    x'' + 2 zeta w x' + w^2 x = -a(t)       (x relative displacement, a ground acceleration)

with the ground acceleration linear between samples.  Outputs per oscillator:

* ``SA``  max |absolute acceleration| = max |-(2 zeta w x' + w^2 x)|   (written to .RS)
* ``SV``  max |relative velocity|
* ``SD``  max |relative displacement|;   ``PSA = w^2 SD``,  ``PSV = w SD``

For oscillator frequencies above ``0.1/dt`` the record is FFT-upsampled by 4 before
integration (D-MOT-07).  The recurrence is evaluated with ``scipy.signal.lfilter`` on
the equivalent discrete state-space system, which is exact for piecewise-linear input.
"""
from __future__ import annotations

from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np
from scipy import signal as _sig


def log_frequencies(f1: float, f2: float, n: int) -> np.ndarray:
    """``n`` log-spaced frequencies from f1 to f2 inclusive (MOTION <freq1>,<freq2>,<fstep>)."""
    if n <= 1:
        return np.asarray([f1], dtype=float)
    return np.logspace(np.log10(f1), np.log10(f2), int(n))


def _nj_coefficients(w: float, z: float, dt: float):
    """Nigam-Jennings A (2x2) and B0, B1 (2,) for state s = (x, v): s+ = A s + B0 a_k + B1 a_{k+1}."""
    if z >= 1.0:
        raise ValueError("damping ratio must be < 1")
    sq = np.sqrt(1.0 - z * z)
    wd = w * sq
    E = np.exp(-z * w * dt)
    S = np.sin(wd * dt)
    C = np.cos(wd * dt)
    A = np.array([[E * (z / sq * S + C), E * S / wd],
                  [-w / sq * E * S, E * (C - z / sq * S)]])
    w2, w3 = w * w, w * w * w
    t1 = (2.0 * z * z - 1.0) / (w2 * dt)
    t2 = 2.0 * z / (w3 * dt)
    b11 = E * ((t1 + z / w) * S / wd + (t2 + 1.0 / w2) * C) - t2
    b12 = -E * (t1 * S / wd + t2 * C) - 1.0 / w2 + t2
    b21 = E * ((t1 + z / w) * (C - z / sq * S) - (t2 + 1.0 / w2) * (wd * S + z * w * C)) + 1.0 / (w2 * dt)
    b22 = -E * (t1 * (C - z / sq * S) - t2 * (wd * S + z * w * C)) - 1.0 / (w2 * dt)
    return A, np.array([b11, b21]), np.array([b12, b22])


def sdof_response(acc: np.ndarray, dt: float, freq: float, zeta: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact piecewise-linear SDOF response (x, v, absolute acceleration) to ground acceleration ``acc``.

    Initial conditions at rest.  Same units as ``acc`` (x in acc-units * s^2).
    """
    return _filter_response(np.asarray(acc, dtype=float), dt, 2.0 * np.pi * freq, zeta)


def _sdof_exact_loop(acc: np.ndarray, dt: float, w: np.ndarray, z: float):
    """Vectorised exact recurrence over many oscillators (columns) -- returns max |SA|, |SV|, |SD|."""
    nw = len(w)
    A = np.zeros((nw, 2, 2))
    B0 = np.zeros((nw, 2))
    B1 = np.zeros((nw, 2))
    for i, wi in enumerate(w):
        A[i], B0[i], B1[i] = _nj_coefficients(wi, z, dt)
    x = np.zeros(nw)
    v = np.zeros(nw)
    sa = np.zeros(nw)
    sv = np.zeros(nw)
    sd = np.zeros(nw)
    a11, a12, a21, a22 = A[:, 0, 0], A[:, 0, 1], A[:, 1, 0], A[:, 1, 1]
    b11, b21 = B0[:, 0], B0[:, 1]
    b12, b22 = B1[:, 0], B1[:, 1]
    w2 = w * w
    c2 = 2.0 * z * w
    for k in range(len(acc) - 1):
        ak, ak1 = acc[k], acc[k + 1]
        x, v = a11 * x + a12 * v + b11 * ak + b12 * ak1, a21 * x + a22 * v + b21 * ak + b22 * ak1
        np.maximum(sd, np.abs(x), out=sd)
        np.maximum(sv, np.abs(v), out=sv)
        np.maximum(sa, np.abs(w2 * x + c2 * v), out=sa)
    return sa, sv, sd


def response_spectrum(acc: np.ndarray, dt: float, freqs: Sequence[float], dampings: Iterable[float],
                      upsample: bool = True, method: str = "filter", substep: bool = True) -> Dict[str, np.ndarray]:
    """Response spectra of ``acc`` (any units) for each damping ratio.

    ``upsample`` enables the 4x FFT up-sampling above 0.1/dt (D-MOT-07); ``substep`` sub-divides
    the (piecewise-linear) input so every oscillator period has >= ``STEPS_PER_PERIOD`` steps.

    Returns a dict with arrays of shape (n_damping, n_freq): ``SA`` (absolute acceleration),
    ``SV`` (relative velocity), ``SD`` (relative displacement), ``PSA``, ``PSV``; plus ``freq``
    and ``damping``.
    """
    acc = np.asarray(acc, dtype=float)
    freqs = np.asarray(freqs, dtype=float)
    damps = np.asarray(list(dampings), dtype=float)
    nd, nf = len(damps), len(freqs)
    SA = np.zeros((nd, nf))
    SV = np.zeros((nd, nf))
    SD = np.zeros((nd, nf))
    hi = freqs > 0.1 / dt if upsample else np.zeros(nf, bool)
    acc_up = None
    if np.any(hi):
        acc_up = _sig.resample(acc, 4 * len(acc))
    for sel, a, h in ((~hi, acc, dt), (hi, acc_up, dt / 4.0)):
        if not np.any(sel) or a is None:
            continue
        idx = np.flatnonzero(sel)
        # Sub-stepping (linear interpolation of the input, which the Nigam-Jennings recurrence
        # treats exactly) so that every oscillator period spans >= STEPS_PER_PERIOD steps: this
        # resolves response peaks that fall between the record samples.
        m_all = np.ones(len(idx), int)
        if upsample and substep:
            m_all = np.clip(np.ceil(STEPS_PER_PERIOD * freqs[idx] * h).astype(int), 1, MAX_SUBSTEPS)
        for m in np.unique(m_all):
            sub = idx[m_all == m]
            if m > 1:
                tt = np.arange((len(a) - 1) * m + 1) / m
                am = np.interp(tt, np.arange(len(a)), a)
                hm = h / m
            else:
                am, hm = a, h
            w = 2.0 * np.pi * freqs[sub]
            for j, z in enumerate(damps):
                if method == "loop":
                    sa, sv, sd = _sdof_exact_loop(am, hm, w, z)
                else:
                    sa = np.zeros(len(w))
                    sv = np.zeros(len(w))
                    sd = np.zeros(len(w))
                    for i, wi in enumerate(w):
                        xx, vv, aa = _filter_response(am, hm, wi, z)
                        sa[i], sv[i], sd[i] = np.max(np.abs(aa)), np.max(np.abs(vv)), np.max(np.abs(xx))
                SA[j, sub], SV[j, sub], SD[j, sub] = sa, sv, sd
    w = 2.0 * np.pi * freqs
    return {"freq": freqs, "damping": damps, "SA": SA, "SV": SV, "SD": SD, "PSA": SD * w * w, "PSV": SD * w}


#: minimum number of integration steps per oscillator period (response-peak sampling error
#: <= 1 - cos(pi/64) = 0.12 %); the input is sub-divided linearly when needed
STEPS_PER_PERIOD = 64
MAX_SUBSTEPS = 64


def _filter_response(acc: np.ndarray, dt: float, w: float, z: float):
    """Exact recurrence via lfilter: the state starts at rest (x0 = v0 = 0)."""
    A, B0, B1 = _nj_coefficients(w, z, dt)
    # Recurrence s_{k+1} = A s_k + B0 a_k + B1 a_{k+1}, s_0 = 0.  Write it as an IIR filter on
    # the input sequence: S(q) = (I - A q^-1)^-1 (B0 q^-1 + B1 q^-1 q) ... use the delayed form
    # s_k = sum A^(k-1-m) (B0 a_m + B1 a_{m+1}).  Let u_k = B0 a_{k-1} + B1 a_k for k >= 1 (u_0 = 0):
    # s_k = A s_{k-1} + u_k.  Each component is lfilter([.. ], [1, -tr, det]) of u.
    n = len(acc)
    u = np.zeros((2, n))
    u[:, 1:] = B0[:, None] * acc[None, :-1] + B1[:, None] * acc[None, 1:]
    tr = A[0, 0] + A[1, 1]
    det = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    den = np.array([1.0, -tr, det])
    # (I - A z^-1)^-1 = adj(I - A z^-1)/det(I - A z^-1);  adj = [[1 - a22 z^-1, a12 z^-1],[a21 z^-1, 1 - a11 z^-1]]
    x = _sig.lfilter([1.0, -A[1, 1]], den, u[0]) + _sig.lfilter([0.0, A[0, 1]], den, u[1])
    v = _sig.lfilter([0.0, A[1, 0]], den, u[0]) + _sig.lfilter([1.0, -A[0, 0]], den, u[1])
    absacc = -(w * w * x + 2.0 * z * w * v)
    return x, v, absacc


# ---------------------------------------------------------------------------------------
# Ground-motion parameters
# ---------------------------------------------------------------------------------------
def integrate(acc: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    """Trapezoidal velocity and displacement (zero initial conditions)."""
    acc = np.asarray(acc, dtype=float)
    vel = np.concatenate([[0.0], np.cumsum(0.5 * (acc[1:] + acc[:-1]) * dt)])
    dis = np.concatenate([[0.0], np.cumsum(0.5 * (vel[1:] + vel[:-1]) * dt)])
    return vel, dis


def arias_intensity(acc: np.ndarray, dt: float, gravity: float = 1.0) -> np.ndarray:
    """Cumulative Arias intensity  I(t) = pi/(2g) * int a^2 dt  (acc in units of gravity*... )."""
    acc = np.asarray(acc, dtype=float)
    return np.pi / (2.0 * gravity) * np.concatenate([[0.0], np.cumsum(0.5 * (acc[1:] ** 2 + acc[:-1] ** 2) * dt)])


def strong_motion_window(acc: np.ndarray, dt: float, lo: float = 0.05, hi: float = 0.75) -> Tuple[float, float]:
    """Times (t_lo, t_hi) at which the normalised Arias intensity reaches ``lo`` and ``hi``."""
    ia = arias_intensity(acc, dt)
    if ia[-1] <= 0:
        return 0.0, 0.0
    r = ia / ia[-1]
    t = np.arange(len(acc)) * dt
    return float(np.interp(lo, r, t)), float(np.interp(hi, r, t))


def band_averaged_psd(acc: np.ndarray, dt: float, band: float = 0.20, window: Optional[Tuple[float, float]] = None):
    """One-sided PSD of the strong-motion part, averaged over +/- ``band`` frequency intervals.

    ``S0(f) = 2 |F(f)|^2 / (2 pi T_D)`` with F the Fourier transform (dt * rfft) of the
    window of duration T_D (requirements section 4.12 item 7).  Returns (f, psd).
    """
    acc = np.asarray(acc, dtype=float)
    if window is None:
        window = strong_motion_window(acc, dt)
    i0, i1 = int(round(window[0] / dt)), int(round(window[1] / dt)) + 1
    seg = acc[i0:i1]
    TD = len(seg) * dt
    F = dt * np.fft.rfft(seg)
    f = np.fft.rfftfreq(len(seg), dt)
    raw = 2.0 * np.abs(F) ** 2 / (2.0 * np.pi * TD)
    out = np.empty_like(raw)
    for i, fi in enumerate(f):
        sel = (f >= fi * (1 - band)) & (f <= fi * (1 + band))
        out[i] = raw[sel].mean() if np.any(sel) else raw[i]
    return f, out
