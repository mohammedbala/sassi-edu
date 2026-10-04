"""Numerical library of the EQUAKE module: spectrum-compatible acceleration histories.

Normative sources: requirements section 4.12, decisions D-EQK-01 ... D-EQK-08, D-GEN-08;
reference criteria R2 sections G.2 (SRP 3.7.1 Rev. 4) and G.3 (RG 1.60).

The generation follows the classical engineering sequence:

1. **Target** design response spectrum (Hz, SA in g) at damping ``zeta``, interpolated
   log-log; below its first frequency the spectral displacement is held constant, above its
   last frequency SA is held constant (zero-period plateau).
2. **Initial motion** (D-EQK-02): a stationary Gaussian process with random phases (PCG64)
   and Fourier amplitudes from the SIMQKE power-spectral-density estimate of Gasparini &
   Vanmarcke (1976), multiplied by a Saragoni-Hart type envelope with a stationary strong
   part; or the Fourier phases of a *seed record* (``accopt`` = 1).
3. **Levy-Wilkinson** frequency-domain iterations: Fourier amplitudes are multiplied by
   ``aim * SA_target / SA_computed`` (spectra by Nigam-Jennings at 100 points/decade,
   :func:`sassi.core.spectra.response_spectrum`), phases kept.  With random phases each
   iteration also applies *iterative clipping and filtering* (clip the stationary process at
   the target zero-period acceleration, keep the new phases, restore the amplitudes); this
   lowers the crest factor so that the PGA can equal the target ZPA while SA keeps the
   spectral amplification of the target ([R], needed for RG 1.60 where SA/PGA = 3.13).
4. **Time-domain refinement** (P1, "AB" algorithm): Al Atik & Abrahamson (2010) improved
   tapered-cosine wavelets are added at the response-peak times of a set of oscillators so
   that the peaks move to ``aim * target``; the influence matrix is computed exactly with
   the same discrete oscillator as the checker and solved with Tikhonov damping.
5. **Drift control and baseline correction** (D-EQK-03, SRP 3.7.1 "no baseline drift"), part
   of *every* step so that the match is built from oscillatory content, never from a
   record-length displacement arc: a zero-phase high-pass (corner 0.5 x the lowest check
   frequency, at least 1/duration) on every frequency-domain update; then zero-mean
   acceleration (f = 0 term removed) and a displacement polynomial (degree 5, i.e. a cubic
   acceleration correction) *fitted by least squares* to the displacement under the
   constraints of zero final velocity and displacement.  The correction is a linear
   projection; the wavelets of step 4 are projected with it before use.
6. **PSD floor** [R] (only when a target PSD is given): deficient +/-20 % band averages of the
   strong-motion (Arias 5-75 %) segment are raised to the target, interleaved with step 4.
7. **Checks** (D-EQK-04): the manual's criteria and SRP 3.7.1 Rev. 4 Option 1 Approach 2,
   reported separately; :func:`drift_measures` quantifies long-period (drift) content.
"""
from __future__ import annotations

import math
from functools import lru_cache
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy import optimize as _opt
from scipy import signal as _sig

from . import spectra as SP

PathLike = Union[str, Path]
G_IN_S2 = 386.08858267716535        # standard gravity in in/s^2
G_CM_S2 = 980.665                   # standard gravity in cm/s^2


# =======================================================================================
# Target spectra
# =======================================================================================
def read_two_columns(path: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Read a two-column text file (frequency, value); '#'/'*' and non-numeric lines skipped."""
    xs, ys = [], []
    for ln in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s[0] in "#*":
            continue
        tok = s.replace(",", " ").split()
        try:
            vals = [float(t.replace("D", "E")) for t in tok[:2]]
        except ValueError:
            continue
        if len(vals) == 2:
            xs.append(vals[0])
            ys.append(vals[1])
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


@dataclass
class TargetSpectrum:
    """Design response spectrum (SA in g) at one damping ratio, log-log interpolated."""

    freq: np.ndarray
    sa: np.ndarray
    damping: float = 0.05

    def __post_init__(self):
        f = np.asarray(self.freq, dtype=float)
        s = np.asarray(self.sa, dtype=float)
        order = np.argsort(f)
        f, s = f[order], s[order]
        if f.size < 2 or np.any(f <= 0) or np.any(s <= 0):
            raise ValueError("target spectrum needs >= 2 points with positive frequency and amplitude")
        if np.any(np.diff(f) <= 0):
            raise ValueError("target spectrum frequencies must be distinct")
        self.freq, self.sa = f, s

    @property
    def fmin(self) -> float:
        return float(self.freq[0])

    @property
    def fmax(self) -> float:
        return float(self.freq[-1])

    def __call__(self, f) -> np.ndarray:
        """SA (g) at ``f``: log-log inside; constant SD below fmin; constant SA above fmax."""
        f = np.asarray(f, dtype=float)
        fl = np.maximum(f, 1e-12)
        out = 10.0 ** np.interp(np.log10(fl), np.log10(self.freq), np.log10(self.sa))
        out = np.where(fl < self.fmin, self.sa[0] * (fl / self.fmin) ** 2, out)
        return out

    @classmethod
    def from_file(cls, path: PathLike, damping: float = 0.05) -> "TargetSpectrum":
        f, s = read_two_columns(path)
        return cls(f, s, damping)


# RG 1.60 Rev. 2 (R2 G.3): damping % -> (A, B, C, D) amplification factors; D is a displacement factor.
RG160_H = {0.5: (1.0, 4.96, 5.95, 3.20), 2.0: (1.0, 3.54, 4.25, 2.50), 5.0: (1.0, 2.61, 3.13, 2.05),
           7.0: (1.0, 2.27, 2.72, 1.88), 10.0: (1.0, 1.90, 2.28, 1.70)}
RG160_V = {0.5: (1.0, 4.96, 5.67, 2.13), 2.0: (1.0, 3.54, 4.05, 1.67), 5.0: (1.0, 2.61, 2.98, 1.37),
           7.0: (1.0, 2.27, 2.59, 1.25), 10.0: (1.0, 1.90, 2.17, 1.13)}
RG160_PGD_IN = 36.0          # maximum ground displacement (in) for 1.0 g


def rg160_control_points(component: str = "H", damping: float = 0.05) -> Dict[str, float]:
    """RG 1.60 control points at ``damping`` (linear interpolation in damping, RG 1.60 C.1/C.2).

    Returns frequencies fA = 33, fB = 9, fC (2.5 H / 3.5 V), fD = 0.25 Hz, the acceleration
    factors A, B, C and the spectral displacement at D, ``SD_D`` (in, for 1 g).
    """
    tab = RG160_H if component.upper().startswith("H") else RG160_V
    d = sorted(tab)
    pct = 100.0 * damping
    if not (d[0] <= pct <= d[-1]):
        raise ValueError(f"RG 1.60 damping must be within {d[0]} - {d[-1]} %")
    vals = [float(np.interp(pct, d, [tab[k][i] for k in d])) for i in range(4)]
    return {"fA": 33.0, "fB": 9.0, "fC": 2.5 if tab is RG160_H else 3.5, "fD": 0.25,
            "A": vals[0], "B": vals[1], "C": vals[2], "SD_D": vals[3] * RG160_PGD_IN}


def rg160_spectrum(freqs, component: str = "H", damping: float = 0.05, pga: float = 1.0) -> np.ndarray:
    """RG 1.60 design spectrum SA (g) anchored to ``pga`` (pseudo-acceleration, tripartite lines).

    Straight lines on log-log axes between the control points D (constant displacement
    below 0.25 Hz), C, B and A; constant ZPA above 33 Hz.
    """
    cp = rg160_control_points(component, damping)
    f = np.asarray(freqs, dtype=float)
    sa_d = (2 * np.pi * cp["fD"]) ** 2 * cp["SD_D"] / G_IN_S2
    xs = np.log10([cp["fD"], cp["fC"], cp["fB"], cp["fA"]])
    ys = np.log10([sa_d, cp["C"], cp["B"], cp["A"]])
    out = 10.0 ** np.interp(np.log10(np.maximum(f, 1e-12)), xs, ys)
    out = np.where(f >= cp["fA"], cp["A"], out)
    out = np.where(f < cp["fD"], (2 * np.pi * f) ** 2 * cp["SD_D"] / G_IN_S2, out)
    return pga * out


def rg160_target_psd(freqs, pga: float = 1.0, units: str = "SI") -> np.ndarray:
    """SRP 3.7.1 Appendix A minimum PSD for the RG 1.60 *horizontal* spectrum (R2 G.2).

    One-sided ``S0(w) = 2|F(w)|^2/(2 pi T_D)``; returned in cm^2/s^3 (``units='SI'``) or
    in^2/s^3 (``'BS'``), scaled by ``pga**2``.
    """
    f = np.asarray(freqs, dtype=float)
    out = np.empty_like(f)
    m1 = f < 2.5
    m2 = (f >= 2.5) & (f < 9.0)
    m3 = (f >= 9.0) & (f < 16.0)
    m4 = f >= 16.0
    out[m1] = 4190.0 * (np.maximum(f[m1], 0.0) / 2.5) ** 0.2
    out[m2] = 4190.0 * (2.5 / f[m2]) ** 1.8
    out[m3] = 418.0 * (9.0 / f[m3]) ** 3
    out[m4] = 74.2 * (16.0 / f[m4]) ** 8
    out = out * pga * pga
    if units.upper() == "BS":
        out = out / 2.54 ** 2
    return out


def write_spectrum_file(path: PathLike, f, sa, header: str = "") -> Path:
    """Write a two-column spectrum file (Hz, SA g), e.g. an RSIN target ``.rsi``."""
    path = Path(path)
    lines = [f"# {header}"] if header else []
    lines += [f"{fi:.6f} {si:.8e}" for fi, si in zip(f, sa)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# =======================================================================================
# Acceptance criteria (D-EQK-04)
# =======================================================================================
def check_grid(f1: float, f2: float, per_decade: int = 100) -> np.ndarray:
    """Log-uniform frequencies from f1 to f2 with at least ``per_decade`` points per decade."""
    n = int(math.ceil(per_decade * math.log10(f2 / f1) - 1e-9)) + 1
    return SP.log_frequencies(f1, f2, max(n, 2))


@dataclass
class RSCheck:
    """Spectral-matching criteria on one frequency band (manual or SRP 3.7.1 Rev. 4 variant)."""

    name: str
    fmin: float
    fmax: float
    npts: int
    per_decade: float
    min_ratio: float
    max_ratio: float
    mean_ratio: float
    n_below_10: int
    n_above_30: int
    max_run_below: int
    f_min_ratio: float
    f_max_ratio: float
    req_fmin: float = 0.0             # band the criteria set requires (SRP: 0.1 Hz - min(50 Hz, Nyquist))
    req_fmax: float = 0.0

    @property
    def band_complete(self) -> bool:
        """The checked band covers the required band (always true for the manual variant)."""
        return self.fmin <= self.req_fmin * (1 + 1e-9) and self.fmax >= self.req_fmax * (1 - 1e-9)

    @property
    def ok_below(self) -> bool:
        return self.n_below_10 == 0

    @property
    def ok_above(self) -> bool:
        return self.n_above_30 == 0

    @property
    def ok_run(self) -> bool:
        return self.max_run_below <= 9

    @property
    def ok_density(self) -> bool:
        return self.per_decade >= 100 - 1e-9

    @property
    def passed(self) -> bool:
        return self.ok_below and self.ok_above and self.ok_run and self.ok_density


def rs_criteria(name: str, freqs: np.ndarray, sa: np.ndarray, target: np.ndarray) -> RSCheck:
    """Never > 10 % below, never > 30 % above, <= 9 adjacent points below the target."""
    r = np.asarray(sa) / np.asarray(target)
    run = best = 0
    for b in r < 1.0:
        run = run + 1 if b else 0
        best = max(best, run)
    dec = math.log10(freqs[-1] / freqs[0]) if len(freqs) > 1 else 1.0
    return RSCheck(name, float(freqs[0]), float(freqs[-1]), len(freqs), (len(freqs) - 1) / max(dec, 1e-12),
                   float(r.min()), float(r.max()), float(r.mean()), int(np.sum(r < 0.9 - 1e-12)),
                   int(np.sum(r > 1.3 + 1e-12)), int(best), float(freqs[int(np.argmin(r))]),
                   float(freqs[int(np.argmax(r))]))


def check_bands(target: TargetSpectrum, dt: float) -> Dict[str, Tuple[float, float]]:
    """Frequency bands of the two criteria sets (D-EQK-04).

    * manual: between the lowest and highest target frequency (limited to Nyquist);
    * SRP 3.7.1 Rev. 4 Approach 2 (b): 0.1 Hz to 50 Hz or the Nyquist frequency.  Above the
      last target frequency the target is extended at constant SA (zero-period plateau, as
      :class:`TargetSpectrum` does), so the band always reaches min(50 Hz, Nyquist); below the
      first target frequency the target is not extrapolated: a target starting above 0.1 Hz
      leaves the band incomplete (reported by :func:`srp_required_band` / ``RSCheck.band_complete``).
    """
    fny = 0.5 / dt
    man = (target.fmin, min(target.fmax, fny))
    srp = (max(0.1, target.fmin), min(50.0, fny))
    if srp[0] >= srp[1]:                 # target entirely above the SRP band: check its first point only
        srp = (min(target.fmin, fny), min(target.fmin, fny) * (1 + 1e-6))
    return {"manual": man, "srp": srp}


def srp_required_band(dt: float) -> Tuple[float, float]:
    """Band required by SRP 3.7.1 Rev. 4 Option 1 Approach 2 (b): 0.1 Hz to min(50 Hz, Nyquist)."""
    return 0.1, min(50.0, 0.5 / dt)


def evaluate_rs(acc_g: np.ndarray, dt: float, target: TargetSpectrum, zeta: float) -> Dict[str, RSCheck]:
    """Both RS criteria sets for a motion (spectra with :func:`sassi.core.spectra.response_spectrum`)."""
    out = {}
    for key, (f1, f2) in check_bands(target, dt).items():
        fc = check_grid(f1, f2)
        sa = SP.response_spectrum(acc_g, dt, fc, [zeta])["SA"][0]
        c = rs_criteria("manual (ACS SASSI)" if key == "manual" else "SRP 3.7.1 Rev. 4 Opt. 1 App. 2",
                        fc, sa, target(fc))
        if key == "srp":
            c.req_fmin, c.req_fmax = srp_required_band(dt)
        else:
            c.req_fmin, c.req_fmax = f1, f2
        out[key] = c
    return out


# =======================================================================================
# Ground-motion parameters, PSD, correlation
# =======================================================================================
def ground_motion_parameters(acc_g: np.ndarray, dt: float, g_len: float) -> Dict[str, float]:
    """PGA (g), PGV, PGD (length units of ``g_len`` = gravity in in/s^2 or cm/s^2), V/A (s),
    AD/V^2, Arias intensity and the 5-75 % strong-motion window (SRP 3.7.1)."""
    a = np.asarray(acc_g, dtype=float)
    vel, dis = SP.integrate(a * g_len, dt)
    pga = float(np.max(np.abs(a)))
    pgv = float(np.max(np.abs(vel)))
    pgd = float(np.max(np.abs(dis)))
    A = pga * g_len
    t5, t75 = SP.strong_motion_window(a, dt)
    ia = SP.arias_intensity(a * g_len, dt, g_len)
    return {"PGA": pga, "PGV": pgv, "PGD": pgd, "V/A": pgv / A if A > 0 else 0.0,
            "AD/V2": A * pgd / pgv ** 2 if pgv > 0 else 0.0, "t5": t5, "t75": t75, "D5-75": t75 - t5,
            "Arias": float(ia[-1]), "vel_end": float(vel[-1]), "dis_end": float(dis[-1])}


def correlation(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation coefficient of two records (zero-mean; the common length is used)."""
    m = min(len(x), len(y))
    x = np.asarray(x[:m], float)
    y = np.asarray(y[:m], float)
    x = x - np.mean(x)
    y = y - np.mean(y)
    d = math.sqrt(float(np.sum(x * x) * np.sum(y * y)))
    return float(np.sum(x * y) / d) if d > 0 else 0.0


def moving_correlation(x: np.ndarray, y: np.ndarray, dt: float, window: float = 2.0,
                       step: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
    """Non-stationary correlation in a moving window (manual: 2 s); returns (centre times, rho)."""
    n = min(len(x), len(y))
    w = max(2, int(round(window / dt)))
    s = max(1, int(round(step / dt)))
    t, r = [], []
    for i0 in range(0, n - w + 1, s):
        t.append((i0 + 0.5 * w) * dt)
        r.append(correlation(x[i0:i0 + w], y[i0:i0 + w]))
    return np.asarray(t), np.asarray(r)


def strong_motion_psd(acc_phys: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    """One-sided PSD of the 5-75 % Arias window, +/-20 % band averaged (req. 4.12 item 7)."""
    return SP.band_averaged_psd(acc_phys, dt, band=0.20)


def psd_check(acc_phys: np.ndarray, dt: float, target_psd: Callable[[np.ndarray], np.ndarray],
              f1: float = 0.3, f2: float = 24.0, level: float = 0.8) -> Dict[str, float]:
    """SRP 3.7.1 Appendix A: band-averaged PSD >= ``level`` x target between f1 and f2."""
    f, p = strong_motion_psd(acc_phys, dt)
    m = (f >= f1) & (f <= f2)
    if not np.any(m):
        return {"min_ratio": float("nan"), "f_min": float("nan"), "frac_below": 1.0, "passed": False, "n": 0}
    r = p[m] / target_psd(f[m])
    k = int(np.argmin(r))
    return {"min_ratio": float(r[k]), "f_min": float(f[m][k]), "frac_below": float(np.mean(r < level)),
            "passed": bool(np.all(r >= level)), "n": int(np.sum(m))}


# =======================================================================================
# Envelope, SIMQKE estimate, wavelets, baseline correction
# =======================================================================================
def saragoni_hart_envelope(t: np.ndarray, duration: float, rise: float = 0.10, stationary: float = 0.45,
                           end_level: float = 0.05, b: float = 2.0) -> np.ndarray:
    """Saragoni-Hart type envelope with a stationary strong-motion part (D-EQK-02).

    The Saragoni & Hart (1974) function ``(t/t1)^b exp(b (1 - t/t1))`` rises to 1 at
    ``t1 = rise * duration``; the envelope then stays at 1 for ``stationary * duration`` and
    decays along the descending branch of the same function to ``end_level`` at the end.
    With the defaults (20 s: 2 s rise, 9 s plateau) the 5-75 % Arias duration is ~9-10 s.
    """
    t = np.asarray(t, dtype=float)
    t1 = rise * duration
    t2 = t1 + stationary * duration
    tau_e = _opt.brentq(lambda x: b * math.log(x) + b * (1.0 - x) - math.log(end_level), 1.0 + 1e-9, 1e3)
    td = max(duration - t2, 1e-9) / (tau_e - 1.0)
    w = np.ones_like(t)
    m = t < t1
    x = np.maximum(t[m] / t1, 1e-300)
    w[m] = x ** b * np.exp(b * (1.0 - x))
    m = t > t2
    x = 1.0 + (t[m] - t2) / td
    w[m] = x ** b * np.exp(b * (1.0 - x))
    return w


def simqke_psd(omega: np.ndarray, sa: np.ndarray, zeta: float, ts: float, p: float = 0.5) -> np.ndarray:
    """SIMQKE (Gasparini & Vanmarcke 1976) one-sided PSD estimate G(w) from a target spectrum.

    Recursive form ``G(w_i) = 4 zeta_s / (pi w_i - 4 zeta_s w_{i-1}) * (SA_i^2 / r_i^2 - int_0^{w_{i-1}} G)``
    with the time-dependent damping ``zeta_s = zeta / (1 - exp(-2 zeta w T_s))`` and the
    Vanmarcke peak factor ``r`` for ``N = (w T_s / 2 pi)(-ln p)^-1`` half-cycles (median p).
    ``sa`` in acceleration units; G in (acceleration units)^2 per rad/s.  Initial estimate only.
    """
    G = np.zeros_like(omega, dtype=float)
    cum = 0.0
    for i in range(1, len(omega)):
        w = omega[i]
        zs = zeta / (1.0 - math.exp(-2.0 * zeta * w * ts))
        N = max(w * ts / (2.0 * math.pi) / (-math.log(p)), 1.01)
        arg = 2.0 * N * (1.0 - math.exp(-((4.0 * zs / math.pi) ** 0.6) * math.sqrt(math.pi * math.log(2.0 * N))))
        r2 = max(2.0 * math.log(max(arg, 1.0 + 1e-9)), 1.0)
        g = 4.0 * zs / max(w * math.pi - 4.0 * zs * omega[i - 1], 1e-12) * (sa[i] ** 2 / r2 - cum)
        G[i] = max(g, 1e-6 * sa[i] ** 2 / w)
        cum += G[i] * (omega[i] - omega[i - 1])
    return G


def aa2010_wavelet(t: np.ndarray, f: float, tj: float, zeta: float, gamma_max: Optional[float] = None) -> np.ndarray:
    """Improved tapered-cosine wavelet of Al Atik & Abrahamson (2010).

    ``f(t) = cos(w' (t - tj + dtj)) exp(-((t - tj + dtj)/gamma)^2)`` with ``w' = w sqrt(1-zeta^2)``,
    ``gamma = 1.178 f^-0.93`` (optionally capped) and ``dtj = atan(sqrt(1-zeta^2)/zeta)/w'``, which
    places the oscillator's response peak at ``tj``.  Its Gaussian taper makes the added
    velocity and displacement nearly zero (no drift).
    """
    wd = 2.0 * math.pi * f * math.sqrt(1.0 - zeta * zeta)
    g = 1.178 * f ** -0.93
    if gamma_max is not None:
        g = min(g, gamma_max)
    dtj = math.atan(math.sqrt(1.0 - zeta * zeta) / zeta) / wd
    tau = np.asarray(t) - tj + dtj
    return np.cos(wd * tau) * np.exp(-(tau / g) ** 2)


BASELINE_DEGREE = 5     # degree of the fitted displacement polynomial (acceleration correction: cubic)


def _integrate_rows(a: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    """Trapezoidal velocity and displacement of every row of ``a`` from rest.

    Same discrete operator as :func:`sassi.core.spectra.integrate` (applied along the last
    axis), so that end conditions imposed with it hold exactly for the checker's integration.
    """
    a = np.asarray(a, dtype=float)
    z = np.zeros(a.shape[:-1] + (1,))
    v = np.concatenate([z, np.cumsum(0.5 * (a[..., 1:] + a[..., :-1]) * dt, axis=-1)], axis=-1)
    d = np.concatenate([z, np.cumsum(0.5 * (v[..., 1:] + v[..., :-1]) * dt, axis=-1)], axis=-1)
    return v, d


@lru_cache(maxsize=32)
def _baseline_operator(n: int, dt: float, degree: int):
    """Pre-computed pieces of the constrained least-squares baseline fit (see :func:`baseline_correct`)."""
    T = (n - 1) * dt
    tau = np.arange(n) / (n - 1)
    ks = np.arange(2, degree + 1)
    # acceleration basis b_k = d^2/dt^2 (t/T)^k; its displacement Db_k is integrated with the
    # same discrete trapezoidal rule as the record (so it equals (t/T)^k up to O(dt^2))
    Bacc = np.stack([k * (k - 1) * tau ** (k - 2) / T ** 2 for k in ks])
    Vb, Db = _integrate_rows(Bacc, dt)
    C = np.stack([Vb[:, -1], Db[:, -1]])                 # end velocity and displacement, (2, m)
    with _quiet_blas():
        Cpinv = np.linalg.pinv(C)                        # minimum-norm particular solution
        _, _, vt = np.linalg.svd(C)
        N = vt[2:].T                                     # null space of the end conditions, (m, m-2)
        Mpinv = np.linalg.pinv(Db.T @ N) if N.shape[1] else np.zeros((0, n))
    if not (np.all(np.isfinite(Cpinv)) and np.all(np.isfinite(Mpinv))):
        raise FloatingPointError("baseline correction: non-finite fit operator")
    return Bacc, Db, Cpinv, N, Mpinv


def baseline_correct(acc: np.ndarray, dt: float, degree: int = BASELINE_DEGREE) -> np.ndarray:
    """Baseline correction of D-EQK-03 (any acceleration units; 1-D record or 2-D array of rows).

    1. complex-frequency step: the f = 0 Fourier term (the mean) of the acceleration is removed;
    2. time-domain step: a displacement polynomial ``p(t) = sum_{k=2..degree} c_k (t/T)^k``
       (zero initial displacement and velocity) is *fitted by least squares* to the
       displacement of the record, under the two constraints that the corrected record,
       integrated from rest (trapezoidal rule, :func:`sassi.core.spectra.integrate`), ends
       with zero velocity and zero displacement; ``p''(t)`` is subtracted from the acceleration.

    With ``degree`` = 3 the two end conditions fix c2 and c3 completely (no fit is left: this
    only closes the record).  The default degree 5 keeps two least-squares degrees of freedom,
    which remove a record-length displacement arc ("baseline drift", SRP 3.7.1); the
    acceleration correction ``p''`` is then a cubic, i.e. still a polynomial of order <= 3 in
    the acceleration (D-EQK-03).  The operation is linear and a projection (applying it twice
    changes nothing), so wavelets that are corrected with it can be added to a corrected
    record without disturbing its end conditions.  The result is compatible: velocity and
    displacement are the integrals of the acceleration.
    """
    a = np.asarray(acc, dtype=float)
    one = a.ndim == 1
    A = np.atleast_2d(a)
    n = A.shape[-1]
    A = A - np.mean(A, axis=-1, keepdims=True)
    if n < max(degree + 2, 5):
        return A[0] if one else A
    Bacc, Db, Cpinv, N, Mpinv = _baseline_operator(n, float(dt), int(degree))
    V, D = _integrate_rows(A, dt)
    E = np.stack([V[:, -1], D[:, -1]], axis=1)          # (rows, 2)
    with _quiet_blas():
        cp = E @ Cpinv.T                                 # satisfies the end conditions
        if N.shape[1]:
            z = (D - cp @ Db) @ Mpinv.T                  # least-squares fit in the null space
            c = cp + z @ N.T
        else:
            c = cp
        out = A - c @ Bacc
    return out[0] if one else out


def highpass_gain(f: np.ndarray, fc: float, order: int = 4) -> np.ndarray:
    """Zero-phase Butterworth high-pass magnitude ``1/sqrt(1 + (fc/f)^(2 order))`` (0 at f = 0)."""
    f = np.asarray(f, dtype=float)
    out = np.zeros_like(f)
    pos = f > 0
    if fc <= 0:
        out[pos] = 1.0
        return out
    out[pos] = 1.0 / np.sqrt(1.0 + (fc / f[pos]) ** (2 * order))
    return out


def drift_measures(acc_g: np.ndarray, dt: float, g_len: float = G_CM_S2) -> Dict[str, float]:
    """Long-period ("baseline drift") content of a record (SRP 3.7.1: no baseline drift).

    ``frac_long``: fraction of the Fourier energy of the displacement (record zero-padded to
    >= 8 times its length) at periods longer than the record duration T, i.e. f < 1/T.  A
    record whose displacement is a single arc spanning the record has most of its energy
    there; a record whose displacement oscillates at the periods of the target spectrum has
    little.  Also returns PGD (length units of ``g_len``).
    """
    a = np.asarray(acc_g, dtype=float) * g_len
    _, d = SP.integrate(a, dt)
    T = (len(a) - 1) * dt
    nf = 1 << int(math.ceil(math.log2(max(8 * len(d), 2))))
    D = np.abs(np.fft.rfft(d, nf)) ** 2
    fr = np.fft.rfftfreq(nf, dt)
    tot = float(D.sum())
    frac = float(D[fr < 1.0 / T].sum() / tot) if tot > 0 else 0.0
    return {"frac_long": frac, "PGD": float(np.max(np.abs(d))), "T": T}


def _quiet_blas():
    """Silence the spurious 'divide by zero/overflow in matmul' warnings that numpy 2.x emits
    with the macOS Accelerate BLAS for finite arrays (numpy issue 26748); results are checked
    for finiteness explicitly by the callers."""
    return np.errstate(divide="ignore", over="ignore", invalid="ignore")


# =======================================================================================
# Oscillator bank (exact discrete influence of added wavelets)
# =======================================================================================
class OscillatorBank:
    """SDOF oscillators used by the wavelet refinement.

    Uses :func:`sassi.core.spectra.sdof_response` (the same exact piecewise-linear
    recurrence as :func:`~sassi.core.spectra.response_spectrum`), so predicted and checked
    responses coincide.  The checker integrates oscillators with ``f <= 0.1/dt`` on the record
    itself and those above on the record up-sampled 4 times (FFT resampling,
    ``scipy.signal.resample``); a bank with ``upsample=4`` does the same, so its predictions are
    exact for the high-frequency oscillators too (resampling is linear).  The absolute-
    acceleration response at (fine) sample n to an input ``u`` is
    ``y[n] = sum_{m=1..n} h[n-m] u[m] + r0[n] u[0]`` (the oscillator starts at rest).
    """

    def __init__(self, freqs: np.ndarray, dt: float, n: int, zeta: float, upsample: int = 1):
        self.freqs = np.asarray(freqs, dtype=float)
        self.up = int(upsample)
        if self.up not in (1, 4):
            raise ValueError("OscillatorBank: upsample must be 1 or 4 (as in the response-spectrum checker)")
        lim = 0.1 / dt if self.up == 1 else 0.5 / dt
        if np.any(self.freqs > lim * (1 + 1e-9)):
            raise ValueError(f"OscillatorBank: frequencies must be <= {'0.1' if self.up == 1 else '0.5'}/dt")
        self.n_rec, self.dt_rec = n, dt
        self.dt, self.n, self.zeta = dt / self.up, n * self.up, zeta
        nf, m = len(self.freqs), self.n
        self.h = np.zeros((nf, m))
        self.r0 = np.zeros((nf, m))
        d0 = np.zeros(m)
        d0[0] = 1.0
        d1 = np.zeros(m)
        d1[1] = 1.0
        for i, f in enumerate(self.freqs):
            self.r0[i] = SP.sdof_response(d0, self.dt, f, zeta)[2]
            r1 = SP.sdof_response(d1, self.dt, f, zeta)[2]
            self.h[i, :-1] = r1[1:]

    def fine(self, x: np.ndarray) -> np.ndarray:
        """Record(s) on the oscillator grid (the checker's 4x FFT up-sampling when ``upsample=4``)."""
        x = np.asarray(x, dtype=float)
        return x if self.up == 1 else _sig.resample(x, self.n, axis=-1)

    def peaks(self, acc: np.ndarray):
        """(|peak|, fine sample index, sign) of the absolute-acceleration response of every oscillator."""
        a = self.fine(acc)
        nf = len(self.freqs)
        pk = np.zeros(nf)
        idx = np.zeros(nf, dtype=int)
        sg = np.zeros(nf)
        for i, f in enumerate(self.freqs):
            y = SP.sdof_response(a, self.dt, f, self.zeta)[2]
            k = int(np.argmax(np.abs(y)))
            pk[i], idx[i], sg[i] = abs(y[k]), k, (1.0 if y[k] >= 0 else -1.0)
        return pk, idx, sg

    def influence(self, F: np.ndarray, idx: np.ndarray) -> np.ndarray:
        """C[i, j] = response of oscillator i at its peak sample idx[i] to wavelet F[j] (record grid)."""
        Ff = self.fine(np.atleast_2d(F))
        K = np.zeros((len(self.freqs), self.n))
        for i in range(len(self.freqs)):
            k = int(idx[i])
            if k > 0:
                K[i, 1:k + 1] = self.h[i, :k][::-1]
            K[i, 0] = self.r0[i, k]
        with _quiet_blas():
            C = K @ Ff.T
        if not np.all(np.isfinite(C)):
            raise FloatingPointError("non-finite influence matrix")
        return C


# =======================================================================================
# Matching driver
# =======================================================================================
@dataclass
class MatchOptions:
    """Tuning of the generation ([R] engineering choices, documented in the listing)."""

    aim: float = 1.05               # match to aim x target (Approach 2: mean ratio slightly > 1)
    lw_iterations: int = 20         # Levy-Wilkinson iterations (req. 4.12 item 3: up to 20)
    clip_level: float = 1.03        # crest-factor clipping level / target ZPA (random phases)
    clip_passes: int = 3
    wavelets: bool = True           # AA2010 time-domain refinement (P1)
    wavelet_iterations: int = 15    # per pass
    wavelet_gain: float = 1.0
    tikhonov: float = 0.2           # regularisation / median |C_ii|
    gamma_cap: float = 1.0          # wavelet width cap as a fraction of the duration
    outer_passes: int = 8           # further refinement passes while a criterion fails (+ PSD floor if a target PSD)
    psd_level: float = 1.0          # PSD floor level (the SRP requirement is 0.8)
    envelope: Tuple[float, float, float] = (0.10, 0.45, 0.05)   # rise, stationary, end level
    # drift control (D-EQK-03, SRP 3.7.1 "no baseline drift"), applied inside the matching loop
    highpass: float = 0.5           # high-pass corner / lowest check frequency (also >= 1/duration)
    highpass_order: int = 4         # Butterworth order of the (zero-phase) high-pass magnitude
    pad: int = 1                    # FFT length = pad x next power of two >= number of samples
    baseline_degree: int = BASELINE_DEGREE
    wavelet_fmax: float = 0.8       # highest wavelet frequency / Nyquist (<= 0.2: record-grid oscillators only)


@dataclass
class MatchResult:
    acc: np.ndarray                 # acceleration (g), baseline corrected
    dt: float
    checks: Dict[str, RSCheck]
    psd: Optional[Dict[str, float]] = None
    log: List[str] = field(default_factory=list)

    @property
    def rs_passed(self) -> bool:
        return all(c.passed for c in self.checks.values())

    @property
    def max_deviation(self) -> float:
        """Largest |SA/target - 1| over the check bands (D-EQK-05 trial selection)."""
        return max(max(abs(c.min_ratio - 1.0), abs(c.max_ratio - 1.0)) for c in self.checks.values())

    def score(self) -> Tuple[int, float]:
        """(number of failed criteria sets, maximum deviation): the trial ranking of D-EQK-05."""
        fails = sum(0 if c.passed else 1 for c in self.checks.values())
        if self.psd is not None and not self.psd["passed"]:
            fails += 1
        return fails, self.max_deviation

    @property
    def violations(self) -> int:
        """Number of offending points: > 10 % below, > 30 % above, adjacent points beyond 9, PSD bands < 80 %."""
        v = sum(c.n_below_10 + c.n_above_30 + max(0, c.max_run_below - 9) for c in self.checks.values())
        if self.psd is not None and not self.psd["passed"]:
            v += max(1, int(round(self.psd.get("frac_below", 0.0) * self.psd.get("n", 1))))
        return int(v)

    def rank(self) -> Tuple[int, int, float]:
        """Ranking of the iterates of one trial: failed sets, then offending points, then max deviation."""
        fails, dev = self.score()
        return fails, self.violations, dev


def plateau_start(target: TargetSpectrum, tol: float = 0.02) -> float:
    """Lowest target frequency from which SA stays within ``tol`` of its last value (zero-period plateau)."""
    s = target.sa
    k = len(s) - 1
    while k > 0 and abs(s[k - 1] / s[-1] - 1.0) <= tol:
        k -= 1
    return float(target.freq[k])


def _summ(c: RSCheck) -> str:
    return (f"min {c.min_ratio:.3f} max {c.max_ratio:.3f} mean {c.mean_ratio:.3f} "
            f"<-10%: {c.n_below_10} >+30%: {c.n_above_30} run below: {c.max_run_below}")


def _interp_ratio(fk: np.ndarray, fc: np.ndarray, ratio: np.ndarray) -> np.ndarray:
    return 10.0 ** np.interp(np.log10(np.maximum(fk, 1e-9)), np.log10(fc), np.log10(ratio))


def _end_taper(n: int, dt: float, length: float = 0.5) -> np.ndarray:
    w = np.ones(n)
    m = min(max(int(round(length / dt)), 1), n // 4)
    w[n - m:] = 0.5 * (1.0 + np.cos(np.pi * np.arange(1, m + 1) / m))
    return w


def match_spectrum(target: TargetSpectrum, dt: float, duration: float, zeta: Optional[float] = None,
                   rng: Optional[np.random.Generator] = None, seed_acc: Optional[np.ndarray] = None,
                   target_psd: Optional[Callable[[np.ndarray], np.ndarray]] = None, g_len: float = G_CM_S2,
                   options: Optional[MatchOptions] = None,
                   progress: Optional[Callable[[str], None]] = None) -> MatchResult:
    """Generate one spectrum-compatible acceleration history (g) of ``duration`` seconds.

    ``rng`` gives random phases (``seed_acc`` None); with ``seed_acc`` (g, sampled at ``dt``) the
    Fourier phases of the seed record are kept in the frequency-domain stage.  ``target_psd``
    (physical units, cm^2/s^3 or in^2/s^3 with ``g_len`` = gravity in cm/s^2 or in/s^2) enables
    the PSD floor.  Deterministic for a given ``rng`` state.

    Drift control (D-EQK-03; SRP 3.7.1 "no baseline drift") is part of every step, so that the
    spectral match is built from real oscillatory content and never from a record-length
    displacement arc:

    * a zero-phase Butterworth high-pass of corner ``f_hp = max(highpass x f_low, 1/duration)``
      (``f_low`` = lowest check frequency) multiplies the Fourier amplitudes of the initial
      motion and of every Levy-Wilkinson update; the LW ratio is never applied to content far
      below the target band;
    * every record is closed by the constrained least-squares displacement fit of
      :func:`baseline_correct` (an exact linear projection);
    * the AA2010 wavelets are passed through the same projection before their influence
      matrix is formed, so the predicted oscillator peaks are those of the corrected record
      and a wavelet truncated by the record ends adds no drift.
    """
    opt = options or MatchOptions()
    zeta = target.damping if zeta is None else zeta
    log: List[str] = []

    def note(msg: str) -> None:
        log.append(msg)
        if progress:
            progress(msg)

    n = int(round(duration / dt)) + 1
    t = np.arange(n) * dt
    nfft = max(int(opt.pad), 1) * (1 << int(math.ceil(math.log2(max(n, 2)))))
    fk = np.fft.rfftfreq(nfft, dt)
    fny = 0.5 / dt
    bands = check_bands(target, dt)
    f_lo = min(b[0] for b in bands.values())
    f_hi = max(b[1] for b in bands.values())
    fc = check_grid(f_lo, f_hi)
    tg = target(fc)
    taper = _end_taper(n, dt)
    zpa_defined = target.fmax >= 25.0
    zpa = float(target(np.array([min(target.fmax, fny)]))[0])
    deg = int(opt.baseline_degree)
    f_hp = max(opt.highpass * f_lo, 1.0 / duration) if opt.highpass > 0 else 0.0
    # applied below the check band only, so that repeated updates do not erode the band itself
    Hhp = np.where(fk >= f_lo, 1.0, highpass_gain(fk, f_hp, opt.highpass_order))

    def bc(x: np.ndarray) -> np.ndarray:
        return baseline_correct(x, dt, deg)

    def rs(a):
        return SP.response_spectrum(a, dt, fc, [zeta])["SA"][0]

    best: Optional[MatchResult] = None

    def consider(a: np.ndarray, label: str) -> MatchResult:
        nonlocal best
        res = MatchResult(a.copy(), dt, evaluate_rs(a, dt, target, zeta))
        if target_psd is not None:
            res.psd = psd_check(a * g_len, dt, target_psd)
        if best is None or res.rank() < best.rank():
            best = res
        return res

    note(f"drift control: high-pass corner {f_hp:.4g} Hz (Butterworth order {opt.highpass_order}, zero phase) "
         f"on every frequency-domain update; constrained least-squares displacement fit of degree {deg} "
         f"(zero end velocity and displacement) after every step; wavelets projected the same way")
    # ---------------- stage 1/2: initial motion and Levy-Wilkinson iterations ----------
    if seed_acc is None:
        if rng is None:
            rng = np.random.default_rng(0)
        rise, stat, eps = opt.envelope
        # D-EQK-02: stationary strong part of at least ~6 s (SRP 3.7.1: 5-75 % Arias >= 6 s)
        stat = max(stat, min(7.5 / duration, 0.75))
        rise = min(rise, max(0.0, 0.95 - stat) / 2.0) if rise + stat > 0.95 else rise
        env = saragoni_hart_envelope(t, duration, rise, stat, eps) * taper
        om = 2.0 * np.pi * fk
        sa_k = target(np.maximum(fk, 1e-3)) * opt.aim
        G = simqke_psd(om, sa_k, zeta, stat * duration)
        phase = rng.uniform(0.0, 2.0 * np.pi, len(fk))
        Z = 0.5 * nfft * np.sqrt(2.0 * G * (om[1] - om[0])) * np.exp(1j * phase) * Hhp
        Z[0] = 0.0
        note(f"initial motion: random phases, SIMQKE amplitudes, Saragoni-Hart envelope "
             f"(rise {rise * duration:.2f} s, stationary {stat * duration:.2f} s); NFFT {nfft}")

        def synth(Zc):
            return bc(env * np.fft.irfft(Zc, nfft)[:n])
    else:
        s = np.zeros(n)
        m = min(n, len(seed_acc))
        s[:m] = seed_acc[:m]
        Z = np.fft.rfft(s * taper, nfft) * Hhp
        note(f"initial motion: seed record ({len(seed_acc)} values), Fourier phases kept; NFFT {nfft}")

        def synth(Zc):
            return bc(taper * np.fft.irfft(Zc, nfft)[:n])

    a = synth(Z)
    for it in range(1, opt.lw_iterations + 1):
        sa = rs(a)
        res = consider(a, f"LW {it}")
        if it == 1 or it % 5 == 0 or it == opt.lw_iterations:
            note(f"LW iteration {it:2d}: " + "; ".join(f"{k}: {_summ(c)}" for k, c in res.checks.items()))
        Z = Z * _interp_ratio(fk, fc, opt.aim * tg / sa) * Hhp
        if seed_acc is None and zpa_defined and opt.clip_passes > 0:
            amp = np.abs(Z)
            lim = opt.clip_level * zpa
            for _ in range(opt.clip_passes):
                sproc = np.clip(np.fft.irfft(Z, nfft), -lim, lim)
                Z = amp * np.exp(1j * np.angle(np.fft.rfft(sproc)))
        a = synth(Z)
    consider(a, "LW final")

    # ---------------- stage 3: AA2010 wavelet refinement ---------------------------------
    if opt.wavelets:
        # Wavelet band: up to wavelet_fmax x Nyquist, but not on the zero-period plateau, where
        # SA = PGA for every oscillator, their peaks coincide and the influence matrix is singular
        # (the LW polish and the crest-factor clipping of stage 2 control that band).  Oscillators
        # up to 0.1/dt run on the record grid, those above on the checker's 4x up-sampled grid, so
        # the influence matrix is exact in both cases.  Wavelets truncated by the record ends
        # (gamma ~ 10 s at 0.1 Hz) are harmless: they are projected by bc() before use.
        fw = min(opt.wavelet_fmax * fny, f_hi, plateau_start(target))
        sub = fc[fc <= fw * (1.0 + 1e-12)]
        banks = []
        if sub.size:
            for off in range(3):
                fm = np.unique(np.concatenate([sub[:2], sub[off::4], sub[-1:]]))
                lo_f, hi_f = fm[fm <= 0.1 / dt], fm[fm > 0.1 / dt]
                parts = []
                if lo_f.size:
                    parts.append(OscillatorBank(lo_f, dt, n, zeta))
                if hi_f.size:
                    parts.append(OscillatorBank(hi_f, dt, n, zeta, upsample=4))
                banks.append((fm, parts, target(fm)))
            wgt = np.clip((fk - 0.8 * fw) / (0.2 * fw), 0.0, 1.0)
            note(f"wavelet refinement: {sub.size} check frequencies {sub[0]:.3g} - {sub[-1]:.3g} Hz (oscillators "
                 f"above 0.1/dt = {0.1 / dt:g} Hz on the 4x up-sampled record, as in the checker)"
                 + (f"; Levy-Wilkinson polish above {0.8 * fw:.3g} Hz" if f_hi > fw else ""))
        else:
            note(f"wavelet refinement skipped: no check frequency in the wavelet band; "
                 f"Levy-Wilkinson polish of the whole band only")
            fw = 0.0
            wgt = np.where(fk > 0, 1.0, 0.0)
        gcap = opt.gamma_cap * duration

        def refine(a: np.ndarray, niter: int, tag: str, psd_each: bool = False) -> np.ndarray:
            res = None
            for it in range(niter):
                if banks:
                    fm, parts, tgm = banks[it % 3]
                    pks = [bk.peaks(a) for bk in parts]
                    pk = np.concatenate([q[0] for q in pks])
                    sg = np.concatenate([q[2] for q in pks])
                    tpk = np.concatenate([q[1] * bk.dt for q, bk in zip(pks, parts)])
                    dR = (opt.aim * tgm - pk) * sg
                    F = bc(np.stack([aa2010_wavelet(t, fm[j], tpk[j], zeta, gcap) for j in range(len(fm))]))
                    C = np.vstack([bk.influence(F, q[1]) for q, bk in zip(pks, parts)])
                    lam = opt.tikhonov * float(np.median(np.abs(np.diag(C))))
                    with _quiet_blas():
                        b = np.linalg.solve(C.T @ C + lam * lam * np.eye(len(fm)), C.T @ dR)
                        da = opt.wavelet_gain * (b @ F)
                    if not np.all(np.isfinite(da)):
                        note(f"{tag}: non-finite wavelet correction, refinement stopped")
                        return a
                    a = bc(a + da)
                if f_hi > fw:   # Levy-Wilkinson polish of the band above the wavelet range
                    sa = rs(a)
                    A = np.fft.rfft(a, nfft)
                    A = A * (1.0 + wgt * (_interp_ratio(fk, fc, opt.aim * tg / sa) - 1.0)) * Hhp
                    a = bc(np.fft.irfft(A, nfft)[:n])
                if psd_each and target_psd is not None:
                    a = bc(psd_floor(a, dt, target_psd, g_len, opt.psd_level))
                res = consider(a, f"{tag} {it + 1}")
                c = list(res.checks.values())
                if all(x.passed for x in c) and min(x.min_ratio for x in c) >= 0.95 and \
                        max(x.max_run_below for x in c) <= 5 and max(x.max_ratio for x in c) <= 1.25:
                    if res.psd is None or res.psd["passed"]:
                        note(f"{tag}: converged after {it + 1} iteration(s): " +
                             "; ".join(f"{k}: {_summ(x)}" for k, x in res.checks.items()))
                        return a
            if res is not None:
                note(f"{tag}: {niter} iteration(s): " + "; ".join(f"{k}: {_summ(x)}" for k, x in res.checks.items()))
            return a

        a = refine(best.acc.copy(), opt.wavelet_iterations, "wavelet pass 1")
        # ---------------- stage 4: further passes while a criterion fails ------------------
        # with a target PSD: PSD floor of the best record, then refinement interleaved with
        # the floor; without: the refinement continues from the last iterate
        for p in range(2, opt.outer_passes + 2):
            if best.score()[0] == 0:              # the kept record already meets every criterion
                break
            if target_psd is not None:
                a = best.acc.copy()
                cur = consider(a, "check")
                for _ in range(3):
                    a = psd_floor(a, dt, target_psd, g_len, opt.psd_level)
                a = bc(a)
                note(f"PSD floor pass {p - 1}: min PSD ratio before {cur.psd['min_ratio']:.3f} "
                     f"at {cur.psd['f_min']:.2f} Hz")
            a = refine(a, opt.wavelet_iterations, f"wavelet pass {p}", psd_each=target_psd is not None)
    final = best
    final.log = log
    dm = drift_measures(final.acc, dt, g_len)
    note(f"selected record: " + "; ".join(f"{k}: {_summ(c)}" for k, c in final.checks.items()) +
         (f"; PSD min ratio {final.psd['min_ratio']:.3f}" if final.psd else "") +
         f"; displacement energy at periods > duration {100 * dm['frac_long']:.1f} %")
    return final


def psd_floor(acc_g: np.ndarray, dt: float, target_psd: Callable[[np.ndarray], np.ndarray], g_len: float,
              level: float = 1.0, f1: float = 0.3, f2: float = 24.0, ramp: float = 0.5) -> np.ndarray:
    """Raise deficient +/-20 % band averages of the strong-motion segment to ``level`` x target [R].

    The Fourier amplitudes of the 5-75 % Arias segment are multiplied, bin by bin, by
    ``sqrt(level T / P_avg)`` where the band average is below the target; the change is
    blended into the record with ``ramp``-second cosine tapers.  Phases are unchanged.
    """
    acc = np.asarray(acc_g, dtype=float) * g_len
    t5, t75 = SP.strong_motion_window(acc, dt)
    i0, i1 = int(round(t5 / dt)), int(round(t75 / dt)) + 1
    seg = acc[i0:i1]
    m = len(seg)
    if m < 8:
        return np.asarray(acc_g, dtype=float).copy()
    f = np.fft.rfftfreq(m, dt)
    F = np.fft.rfft(seg)
    raw = 2.0 * np.abs(dt * F) ** 2 / (2.0 * np.pi * m * dt)
    c = np.concatenate([[0.0], np.cumsum(raw)])
    lo = np.searchsorted(f, 0.8 * f, side="left")
    hi = np.searchsorted(f, 1.2 * f, side="right")
    cnt = np.maximum(hi - lo, 1)
    avg = (c[hi] - c[lo]) / cnt
    tgt = np.where(f > 0, target_psd(np.maximum(f, 1e-6)), 0.0)
    sel = (f >= 0.8 * f1) & (f <= 1.2 * f2) & (avg < level * tgt)
    g = np.ones_like(f)
    g[sel] = np.sqrt(np.clip(level * tgt[sel] / np.maximum(avg[sel], 1e-300), 1.0, 9.0))
    new = np.fft.irfft(F * g, m)
    w = np.ones(m)
    k = min(max(int(round(ramp / dt)), 1), m // 4)
    w[:k] = 0.5 * (1.0 - np.cos(np.pi * np.arange(k) / k))
    w[m - k:] = w[:k][::-1]
    out = acc.copy()
    out[i0:i1] += w * (new - seg)
    return out / g_len


def correlate_pair(x: np.ndarray, y: np.ndarray, dt: float, times: Sequence[float], values: Sequence[float],
                   window: float = 2.0) -> np.ndarray:
    """Mix ``y`` with ``x`` so that the 2-s windows have the target correlation rho(t) (D-EQK-06, P1).

    In each window (Hann windows, 50 % overlap) ``y' = rho x_n + sqrt(1 - rho^2) y`` with ``x_n``
    = x scaled to the RMS of y in that window; rho(t) is piecewise linear in the CORR pairs.
    The result is to be re-matched with ``y'`` as a seed record.
    """
    n = min(len(x), len(y))
    x, y = np.asarray(x[:n], float), np.asarray(y[:n], float)
    w = max(4, int(round(window / dt)))
    hop = w // 2
    out = np.zeros(n)
    wsum = np.zeros(n)
    han = np.hanning(w + 2)[1:-1]
    tt = np.asarray(times, float)
    vv = np.asarray(values, float)
    for i0 in range(-hop, n, hop):
        i_a, i_b = max(i0, 0), min(i0 + w, n)
        if i_b - i_a < 2:
            continue
        h = han[i_a - i0:i_b - i0]
        tc = 0.5 * (i_a + i_b) * dt
        rho = float(np.clip(np.interp(tc, tt, vv), -1.0, 1.0)) if len(tt) else 0.0
        xs, ys = x[i_a:i_b], y[i_a:i_b]
        sx = math.sqrt(float(np.mean(xs * xs))) or 1.0
        sy = math.sqrt(float(np.mean(ys * ys)))
        mix = rho * xs * (sy / sx) + math.sqrt(max(1.0 - rho * rho, 0.0)) * ys
        out[i_a:i_b] += h * mix
        wsum[i_a:i_b] += h
    return np.where(wsum > 1e-12, out / np.maximum(wsum, 1e-12), y)
