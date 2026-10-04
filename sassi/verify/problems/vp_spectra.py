"""VP-30: response spectra by the exact piecewise-linear (Nigam-Jennings) recurrence (R2 G.1).

Closed-form checks of :func:`sassi.core.spectra.response_spectrum` (used by MOTION, SOIL and
EQUAKE): step, velocity impulse, harmonic excitation at resonance and the high-frequency
limit SA -> PGA.  Tolerances: 1e-5 relative (dt <= T/100) and 1e-3 for the harmonic case.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from ...core import spectra as SP
from .. import VPResult, problem


def step_psa_exact(zeta: float) -> float:
    """PSA/a0 for a step of ground acceleration from rest: 1 + exp(-zeta pi / sqrt(1 - zeta^2))."""
    return 1.0 + math.exp(-zeta * math.pi / math.sqrt(1.0 - zeta * zeta))


def impulse_psv_exact(zeta: float) -> float:
    """PSV/dV for a velocity impulse: exp(-zeta arccos(zeta) / sqrt(1 - zeta^2))."""
    return math.exp(-zeta * math.acos(zeta) / math.sqrt(1.0 - zeta * zeta))


@problem("VP-30", "Response spectra (Nigam-Jennings) closed forms", tier="P0", modules=["MOTION", "EQUAKE", "SOIL"],
         source="R2 G.1")
def vp30(workdir: Path) -> VPResult:
    r = VPResult()
    z = 0.05
    # ---- step input: a(t) = a0 for t >= 0 (exact for the piecewise-linear recurrence)
    r.check("closed form 1 + exp(-zeta pi/sqrt(1-zeta^2)), zeta = 5 %", step_psa_exact(z), 1.854468, atol=5e-7,
            note="R2 value printed to 6 decimals")
    a0 = 0.37
    for f in (1.0, 5.0, 10.0):
        T = 1.0 / f
        dt = T / 1000.0
        acc = np.full(int(round(2.0 * T / dt)) + 1, a0)
        rs = SP.response_spectrum(acc, dt, [f], [z])
        r.check(f"step: PSA/a0 at {f:g} Hz (dt = T/1000)", rs["PSA"][0, 0] / a0, step_psa_exact(z), rtol=1e-5)
        r.notes.append(f"step at {f:g} Hz: PSA/a0 = {rs['PSA'][0, 0] / a0:.6f}, true absolute SA/a0 = "
                       f"{rs['SA'][0, 0] / a0:.6f} (SA >= PSA; they differ by O(zeta))")
        # peak time aligned with a sample: exact to round-off
        wd = 2 * math.pi * f * math.sqrt(1 - z * z)
        dt2 = (math.pi / wd) / 200.0
        acc2 = np.full(401, a0)
        rs2 = SP.response_spectrum(acc2, dt2, [f], [z])
        r.check(f"step: PSA/a0 at {f:g} Hz (peak on a sample)", rs2["PSA"][0, 0] / a0, step_psa_exact(z), rtol=1e-10)
    # ---- velocity impulse: triangle of area dV over 2 dt (dt = T/2000)
    r.check("closed form exp(-zeta acos(zeta)/sqrt(1-zeta^2)), zeta = 5 %", impulse_psv_exact(z), 0.92669,
            atol=5e-6, note="R2 value printed to 5 decimals")
    dV = 1.3
    for f in (0.5, 2.0, 8.0):
        T = 1.0 / f
        dt = T / 2000.0
        acc = np.zeros(int(round(1.5 * T / dt)) + 1)
        acc[1] = dV / dt
        rs = SP.response_spectrum(acc, dt, [f], [z])
        r.check(f"impulse: PSV/dV at {f:g} Hz (dt = T/2000)", rs["PSV"][0, 0] / dV, impulse_psv_exact(z), rtol=1e-5)
    # ---- harmonic excitation at resonance, long duration
    A = 0.2
    for zz, dur, ref_psa in ((0.05, 80.0, 10.000), (0.02, 150.0, 25.000)):
        f = 1.0
        dt = 1.0 / (200.0 * f)
        t = np.arange(int(round(dur / dt)) + 1) * dt
        acc = A * np.sin(2 * math.pi * f * t)
        rs = SP.response_spectrum(acc, dt, [f], [zz])
        r.check(f"harmonic at resonance: PSA/A, zeta = {zz:g}", rs["PSA"][0, 0] / A, ref_psa, rtol=1e-3)
        r.check(f"harmonic at resonance: closed form 1/(2 zeta), zeta = {zz:g}", 1.0 / (2 * zz), ref_psa, rtol=1e-12)
        if zz == 0.05:
            r.check("harmonic at resonance: max |abs. acc.|/A = sqrt(1+4 zeta^2)/(2 zeta), 5 %", rs["SA"][0, 0] / A,
                    10.0499, rtol=1e-3)
    # ---- high-frequency limit SA -> PGA
    dt = 0.001
    n = 10000                       # exactly 10 periods: the FFT up-sampling of the checker is exact
    t = np.arange(n) * dt
    acc = 0.3 * np.sin(2 * math.pi * 1.0 * t)
    rs = SP.response_spectrum(acc, dt, [1000.0, 2000.0], [z])
    for k, f in enumerate((1000.0, 2000.0)):
        r.check(f"SA({f:g} Hz) / PGA, 1 Hz sinusoid (oscillator up-sampled branch)", rs["SA"][0, k] / 0.3, 1.0,
                rtol=1e-5)
    rng = np.random.default_rng(3)
    dt = 0.002
    n = 8192
    spec = np.fft.rfft(rng.standard_normal(n))
    fr = np.fft.rfftfreq(n, dt)
    spec[fr > 5.0] = 0.0
    rec = np.fft.irfft(spec, n)
    rec *= 0.25 / np.max(np.abs(rec))
    rs = SP.response_spectrum(rec, dt, [40.0, 500.0], [z])
    r.check("SA(500 Hz) / PGA, band-limited random record (5 Hz)", rs["SA"][0, 1] / 0.25, 1.0, rtol=1e-3,
            note="PGA of the sampled record; the up-sampled interpolant peaks slightly higher")
    r.notes.append("Spectra from sassi.core.spectra.response_spectrum (lead-owned); PSA = w^2 SD and the true "
                   "absolute-acceleration SA are both reported for the step input.")
    return r
