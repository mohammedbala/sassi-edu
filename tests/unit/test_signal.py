"""Unit tests of sassi.core.signal (FFT conventions D-CNV-02, D-CNV-06, D-MOT-08/09)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.core import signal as S


def test_fourier_grid():
    assert np.allclose(S.fourier_grid(8, 0.5), [0, 0.25, 0.5, 0.75, 1.0])
    assert S.fourier_grid(4096, 0.005)[1] == pytest.approx(0.048828125)
    assert np.allclose(S.omega_grid(8, 0.5), 2 * np.pi * S.fourier_grid(8, 0.5))
    with pytest.raises(ValueError):
        S.fourier_grid(1, 0.01)


def test_pad_record():
    a = np.arange(5.0)
    p = S.pad_record(a, 8)
    assert p.tolist() == [0, 1, 2, 3, 4, 0, 0, 0]
    with pytest.raises(ValueError):
        S.pad_record(np.ones(9), 8)


def test_convolve_identity_and_delay():
    rng = np.random.default_rng(1)
    n, dt = 256, 0.01
    a = rng.standard_normal(n)
    A = np.fft.rfft(a)
    assert np.allclose(S.convolve(np.ones(n // 2 + 1), A), a, atol=1e-14)
    # e^{-i w tau} delays by tau (exp(+i w t) convention): circular shift by 7 samples
    w = S.omega_grid(n, dt)
    r = S.convolve(np.exp(-1j * w * 7 * dt), A)
    assert np.allclose(r, np.roll(a, 7), atol=1e-12)
    two = S.convolve(np.ones((n // 2 + 1, 2)), A)
    assert two.shape == (n, 2) and np.allclose(two[:, 1], a)


def test_hermitian_bins():
    H = np.ones(5, dtype=complex) * (1 + 1j)
    out = S.hermitian_bins(H, 8)
    assert out[-1] == 1.0 and out[0] == 1.0 and out[2] == 1 + 1j


def test_displacement_spectrum_of_a_harmonic():
    n, dt, g, k0 = 1024, 0.01, 9.81, 20
    t = np.arange(n) * dt
    f = S.fourier_grid(n, dt)
    w0 = 2 * np.pi * f[k0]
    a = np.sin(w0 * t)                       # acceleration in g, periodic on the grid
    U = S.displacement_spectrum(np.fft.rfft(a), f, g)
    assert U[0] == 0
    assert np.allclose(np.fft.irfft(U, n), -g * a / w0 ** 2, atol=1e-12)
    assert np.allclose(S.derivative_factor(f, 2), -(2 * np.pi * f) ** 2)


def test_output_length():
    assert S.output_length(4096, 0.005, 0.0) == 4096
    assert S.output_length(4096, 0.005, 10.0) == 2400       # 1.2 * dur
    assert S.output_length(4096, 0.005, 50.0) == 4096       # clipped to NFFT


def test_integrate_trapz_exact_for_linear():
    dt = 0.1
    t = np.arange(11) * dt
    v, d = S.integrate_trapz(2.0 + 3.0 * t, dt)
    assert np.allclose(v, 2 * t + 1.5 * t ** 2)


def test_baseline_correction_removes_linear_baseline():
    dt = 0.01
    t = np.arange(2000) * dt
    a = 0.1 * np.sin(2 * np.pi * 1.3 * t) * np.exp(-0.2 * t)
    drift = 0.004 - 0.0003 * t
    b_clean = S.baseline_correction(a, dt, scale=9.81)
    b_drift = S.baseline_correction(a + drift, dt, scale=9.81)
    assert np.allclose(b_drift.acc, b_clean.acc, atol=1e-13)        # a linear baseline is removed exactly
    assert b_drift.c1 - b_clean.c1 == pytest.approx(0.004, rel=1e-9)
    assert b_drift.c2 - b_clean.c2 == pytest.approx(-0.0003, rel=1e-8)
    assert b_drift.final_disp == pytest.approx(b_clean.final_disp, abs=1e-10)
    assert np.allclose(S.hudson_baseline(a + drift, dt), b_drift.acc)
    # the corrected velocity has no residual trend (least-squares residual is orthogonal to t, t^2/2)
    B = np.stack([t, 0.5 * t * t], axis=1)
    assert np.allclose(B.T @ (b_drift.vel / 9.81), 0.0, atol=1e-9 * np.abs(B).sum())
