"""Unit tests of sassi.core.equake_lib (targets, criteria, envelope, wavelets, baseline, matching)."""
from __future__ import annotations


import numpy as np
import pytest

from sassi.core import equake_lib as EL
from sassi.core import spectra as SP


# ------------------------------------------------------------------ RG 1.60 and App. A
def test_rg160_control_points_and_shape():
    f = np.array([0.1, 0.25, 2.5, 9.0, 33.0, 50.0])
    h = EL.rg160_spectrum(f, "H")
    assert h[2] == pytest.approx(3.13) and h[3] == pytest.approx(2.61) and h[4] == 1.0 and h[5] == 1.0
    assert h[1] == pytest.approx(0.4716, abs=5e-5)
    assert h[0] * EL.G_IN_S2 / (2 * np.pi * 0.1) ** 2 == pytest.approx(2.05 * 36.0)   # constant SD below D
    v = EL.rg160_spectrum(np.array([0.25, 3.5, 9.0]), "V")
    assert v[0] == pytest.approx(0.3152, abs=5e-5) and v[1] == pytest.approx(2.98) and v[2] == pytest.approx(2.61)
    assert EL.rg160_spectrum(np.array([2.5]), "H", pga=0.3)[0] == pytest.approx(0.3 * 3.13)
    # continuity of the log-log lines and damping interpolation (RG 1.60: linear in damping)
    ff = np.logspace(-1.5, 2, 2000)
    lg = np.log(EL.rg160_spectrum(ff, "H"))
    assert np.max(np.abs(np.diff(lg))) < 0.05
    cp = EL.rg160_control_points("H", 0.035)
    assert cp["C"] == pytest.approx(0.5 * (4.25 + 3.13))
    with pytest.raises(ValueError):
        EL.rg160_control_points("H", 0.2)


def test_app_a_psd_continuity_and_units():
    for fb in (2.5, 9.0, 16.0):
        lo, hi = EL.rg160_target_psd(np.array([fb * (1 - 1e-9), fb]))
        # the SRP coefficients are rounded: 418 (9/16)^3 = 74.40 vs 74.2 at 16 Hz (0.26 %)
        assert lo == pytest.approx(hi, rel=3e-3)
    si = EL.rg160_target_psd(np.array([5.0]), 0.5, "SI")[0]
    bs = EL.rg160_target_psd(np.array([5.0]), 0.5, "BS")[0]
    assert si == pytest.approx(0.25 * 4190.0 * (2.5 / 5.0) ** 1.8) and bs == pytest.approx(si / 2.54 ** 2)


def test_target_spectrum_interpolation(tmp_path):
    f = np.array([0.2, 1.0, 10.0, 33.0])
    sa = np.array([0.1, 1.0, 2.0, 1.0])
    EL.write_spectrum_file(tmp_path / "t.rsi", f, sa, header="test")
    t = EL.TargetSpectrum.from_file(tmp_path / "t.rsi")
    assert t(np.array([1.0]))[0] == pytest.approx(1.0)
    assert t(np.array([np.sqrt(10.0)]))[0] == pytest.approx(np.sqrt(2.0))    # log-log midpoint
    assert t(np.array([0.1]))[0] == pytest.approx(0.1 * 0.25)               # constant SD below fmin
    assert t(np.array([80.0]))[0] == pytest.approx(1.0)                     # constant SA above fmax
    with pytest.raises(ValueError):
        EL.TargetSpectrum(np.array([1.0, 1.0]), np.array([1.0, 2.0]))


# ------------------------------------------------------------------ criteria
def test_rs_criteria_counts():
    f = EL.check_grid(1.0, 10.0)
    assert len(f) == 101
    tg = np.ones_like(f)
    sa = np.ones_like(f) * 1.05
    sa[10:19] = 0.95                     # 9 adjacent below: allowed
    c = EL.rs_criteria("x", f, sa, tg)
    assert c.passed and c.max_run_below == 9 and c.per_decade == pytest.approx(100.0)
    sa[19] = 0.99                        # 10 adjacent below: fails
    sa[50] = 1.31
    sa[60] = 0.89
    c = EL.rs_criteria("x", f, sa, tg)
    assert not c.ok_run and c.n_above_30 == 1 and c.n_below_10 == 1 and not c.passed
    assert c.f_min_ratio == pytest.approx(f[60])


def test_check_bands():
    t = EL.TargetSpectrum(np.array([0.05, 100.0]), np.array([0.1, 1.0]))
    b = EL.check_bands(t, 0.005)
    assert b["manual"] == (0.05, 100.0) and b["srp"] == (0.1, 50.0)
    b = EL.check_bands(t, 0.02)
    assert b["srp"][1] == 25.0 and b["manual"][1] == 25.0
    # SRP (b) reaches min(50 Hz, Nyquist) through the zero-period plateau even when the target
    # file ends at 33 Hz; a target starting above 0.1 Hz leaves the band incomplete
    t33 = EL.TargetSpectrum(np.array([0.1, 2.5, 33.0]), np.array([0.2, 3.13, 1.0]))
    assert EL.check_bands(t33, 0.005)["srp"] == (0.1, 50.0)
    assert EL.check_bands(t33, 0.005)["manual"] == (0.1, 33.0)
    assert t33(np.array([45.0]))[0] == 1.0
    a = np.random.default_rng(0).standard_normal(4001) * 0.2
    c = EL.evaluate_rs(a, 0.005, t33, 0.05)
    assert c["srp"].band_complete and c["srp"].fmax == pytest.approx(50.0) and c["manual"].band_complete
    t02 = EL.TargetSpectrum(np.array([0.2, 2.5, 33.0]), np.array([0.4, 3.13, 1.0]))
    c = EL.evaluate_rs(a, 0.005, t02, 0.05)
    assert not c["srp"].band_complete and (c["srp"].req_fmin, c["srp"].req_fmax) == (0.1, 50.0)
    assert EL.plateau_start(t33) == 33.0 and EL.plateau_start(t02) == 33.0


def test_correlation_tools():
    rng = np.random.default_rng(0)
    x, y = rng.standard_normal(4000), rng.standard_normal(4000)
    assert EL.correlation(x, x) == pytest.approx(1.0)
    assert EL.correlation(x, -2 * x + 3) == pytest.approx(-1.0)
    assert abs(EL.correlation(x, y)) < 0.06
    t, r = EL.moving_correlation(x, x, 0.005)
    np.testing.assert_allclose(r, 1.0)
    assert t[0] == pytest.approx(1.0)
    mixed = EL.correlate_pair(x, y, 0.005, [0.0, 30.0], [0.5, 0.5])
    assert EL.correlation(x, mixed) == pytest.approx(0.5, abs=0.06)


def test_ground_motion_parameters_harmonic():
    dt = 0.001
    t = np.arange(20001) * dt
    a = 0.2 * np.sin(2 * np.pi * t)                      # g; exactly 20 periods
    p = EL.ground_motion_parameters(a, dt, 981.0)
    assert p["PGA"] == pytest.approx(0.2)
    assert p["PGV"] == pytest.approx(2 * 0.2 * 981.0 / (2 * np.pi), rel=1e-4)   # v from rest: (A/w)(1 - cos)
    assert p["D5-75"] == pytest.approx(0.70 * 20.0, abs=0.05)


# ------------------------------------------------------------------ envelope, PSD estimate, wavelets
def test_saragoni_hart_envelope():
    t = np.linspace(0, 20, 4001)
    w = EL.saragoni_hart_envelope(t, 20.0)
    assert w[0] == 0.0 and w.max() == pytest.approx(1.0)
    assert np.all(w[(t >= 2.0) & (t <= 11.0)] == 1.0)
    assert w[-1] == pytest.approx(0.05, rel=1e-6)
    t5, t75 = SP.strong_motion_window(w, t[1] - t[0])
    assert t75 - t5 > 6.0                                 # strong part of the intensity envelope


def test_simqke_psd_positive_and_scales():
    om = 2 * np.pi * np.linspace(0, 50, 2001)
    sa = EL.rg160_spectrum(np.maximum(om / (2 * np.pi), 1e-3), "H") * 981.0
    g1 = EL.simqke_psd(om, sa, 0.05, 9.0)
    g2 = EL.simqke_psd(om, 2 * sa, 0.05, 9.0)
    assert np.all(g1[1:] > 0)
    np.testing.assert_allclose(g2[1:], 4 * g1[1:], rtol=1e-10)


def test_aa2010_wavelet_nearly_no_drift():
    """Residual velocity of the tapered cosine = gamma sqrt(pi) exp(-(w' gamma / 2)^2) (tiny)."""
    dt = 0.002
    t = np.arange(25001) * dt
    for f in (0.5, 2.0, 10.0):
        w = EL.aa2010_wavelet(t, f, 25.0, 0.05)
        v, d = SP.integrate(w, dt)
        g = 1.178 * f ** -0.93
        wd = 2 * np.pi * f * np.sqrt(1 - 0.05 ** 2)
        exact = g * np.sqrt(np.pi) * np.exp(-(wd * g / 2) ** 2)     # integral over the whole time axis
        assert abs(v[-1]) == pytest.approx(exact, rel=0.05, abs=1e-12)
        assert abs(v[-1]) < 1e-4 * np.max(np.abs(v))


def test_oscillator_bank_influence_is_exact():
    dt, n = 0.01, 1500
    t = np.arange(n) * dt
    rng = np.random.default_rng(4)
    a = np.convolve(rng.standard_normal(n), np.ones(4) / 4, "same") * np.hanning(n)
    fr = np.array([0.3, 1.0, 4.0, 9.5])
    bank = EL.OscillatorBank(fr, dt, n, 0.05)
    pk, idx, sg = bank.peaks(a)
    rs = SP.response_spectrum(a, dt, fr, [0.05], substep=False)["SA"][0]   # bank = checker without sub-steps
    np.testing.assert_allclose(pk, rs, rtol=1e-12)
    F = np.stack([EL.aa2010_wavelet(t, f, idx[j] * dt, 0.05) for j, f in enumerate(fr)])
    C = bank.influence(F, idx)
    for i, f in enumerate(fr):
        for j in range(len(fr)):
            y = SP.sdof_response(F[j], dt, f, 0.05)[2]
            assert C[i, j] == pytest.approx(y[idx[i]], abs=1e-10 * max(1.0, np.max(np.abs(y))))
    with pytest.raises(ValueError):
        EL.OscillatorBank(np.array([20.0]), dt, n, 0.05)


def test_upsampled_oscillator_bank_is_exact():
    """Above 0.1/dt the checker integrates the 4x FFT-resampled record; the bank does the same."""
    dt, n = 0.01, 1500
    t = np.arange(n) * dt
    rng = np.random.default_rng(4)
    a = rng.standard_normal(n) * np.hanning(n)
    fr = np.array([12.0, 20.0, 33.0, 45.0])
    bank = EL.OscillatorBank(fr, dt, n, 0.05, upsample=4)
    pk, idx, sg = bank.peaks(a)
    np.testing.assert_allclose(pk, SP.response_spectrum(a, dt, fr, [0.05], substep=False)["SA"][0], rtol=1e-12)
    F = EL.baseline_correct(np.stack([EL.aa2010_wavelet(t, f, idx[j] * bank.dt, 0.05) for j, f in enumerate(fr)]),
                            dt)
    C = bank.influence(F, idx)
    b = np.array([0.1, -0.2, 0.05, 0.3])
    with np.errstate(all="ignore"):           # spurious Accelerate-BLAS matmul warnings (numpy 26748)
        a1 = a + b @ F
    y0 = np.array([SP.sdof_response(bank.fine(a), bank.dt, f, 0.05)[2][idx[i]] for i, f in enumerate(fr)])
    y1 = np.array([SP.sdof_response(bank.fine(a1), bank.dt, f, 0.05)[2][idx[i]] for i, f in enumerate(fr)])
    np.testing.assert_allclose(y1 - y0, C @ b, atol=1e-12 * np.max(np.abs(y0)))
    with pytest.raises(ValueError):
        EL.OscillatorBank(np.array([60.0]), dt, n, 0.05, upsample=4)


# ------------------------------------------------------------------ baseline correction (VP-31-like)
def test_baseline_correction_constant_offset():
    dt = 0.01
    t = np.arange(3001) * dt
    a = 0.3 * np.sin(2 * np.pi * 1.3 * t) * np.exp(-0.1 * t) + 0.02          # constant offset
    ac = EL.baseline_correct(a, dt)
    v, d = SP.integrate(ac, dt)
    assert abs(v[-1]) < 1e-12 and abs(d[-1]) < 1e-10
    assert abs(np.mean(v)) < 0.05 * np.max(np.abs(v))                          # no velocity trend left
    # D-EQK-03: the correction is a polynomial of order <= 3 in acceleration (displacement degree 5)
    corr = a - ac
    np.testing.assert_allclose(np.diff(corr, 4), 0.0, atol=1e-12 * np.max(np.abs(corr)) + 1e-15)
    # degree 3: the end conditions alone fix the correction (order <= 1 in acceleration)
    a3 = EL.baseline_correct(a, dt, degree=3)
    np.testing.assert_allclose(np.diff(a - a3, 2), 0.0, atol=1e-12)
    v3, d3 = SP.integrate(a3, dt)
    assert abs(v3[-1]) < 1e-12 and abs(d3[-1]) < 1e-10


def test_baseline_fit_removes_record_length_displacement_arc():
    """D-EQK-03 / SRP 3.7.1 'no baseline drift': the constrained least-squares displacement fit
    removes an arc spanning the record; the end-condition-only correction (degree 3) cannot."""
    dt, n = 0.01, 2001
    t = np.arange(n) * dt
    T = t[-1]
    d_osc = np.exp(-((t - 10.0) / 2.5) ** 2) * np.sin(2 * np.pi * 1.1 * t)     # oscillatory, compact
    x = t / T
    d_arc = 40.0 * x ** 2 * (1 - x) ** 2                                        # at rest at both ends
    d = d_osc + d_arc
    a = np.gradient(np.gradient(d, dt), dt)
    for deg, removed in ((5, True), (3, False)):
        _, dc = SP.integrate(EL.baseline_correct(a, dt, degree=deg), dt)
        err = np.max(np.abs(dc - d_osc))
        if removed:
            assert err < 0.02 * np.max(np.abs(d_osc))
        else:
            assert err > 0.5 * np.max(d_arc)
    m = EL.drift_measures(a, dt, 1.0)
    mc = EL.drift_measures(EL.baseline_correct(a, dt), dt, 1.0)
    assert m["frac_long"] > 0.9 and mc["frac_long"] < 0.05


def test_baseline_correct_is_a_linear_projection_on_rows():
    dt, n = 0.005, 1201
    rng = np.random.default_rng(5)
    A = np.cumsum(rng.standard_normal((3, n)), axis=1) * 1e-2 + 0.1
    B = EL.baseline_correct(A, dt)
    for i in range(3):
        np.testing.assert_allclose(B[i], EL.baseline_correct(A[i], dt), atol=1e-13)
    np.testing.assert_allclose(EL.baseline_correct(B, dt), B, atol=1e-12 * np.max(np.abs(B)))     # idempotent
    np.testing.assert_allclose(EL.baseline_correct(A[0] + 2 * A[1], dt), B[0] + 2 * B[1], atol=1e-12)   # linear
    V, D = (np.stack(x) for x in zip(*(SP.integrate(b, dt) for b in B)))
    assert np.max(np.abs(V[:, -1])) < 1e-12 and np.max(np.abs(D[:, -1])) < 1e-12


def test_highpass_gain_and_drift_measure():
    f = np.array([0.0, 0.05, 0.1, 1.0])
    h = EL.highpass_gain(f, 0.1, 4)
    assert h[0] == 0.0 and h[2] == pytest.approx(2 ** -0.5) and h[3] == pytest.approx(1.0, abs=1e-7)
    assert h[1] == pytest.approx(1 / np.sqrt(1 + 2.0 ** 8))
    assert np.all(EL.highpass_gain(f[1:], 0.0) == 1.0)


# ------------------------------------------------------------------ matching (small, fast)
def _small_match(seed, **kw):
    f = EL.check_grid(0.3, 50.0, 20)
    tgt = EL.TargetSpectrum(f, EL.rg160_spectrum(f, "H", 0.05, 0.3), 0.05)
    opts = EL.MatchOptions(lw_iterations=8, wavelet_iterations=6, **kw)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed, spawn_key=(1,))))
    return EL.match_spectrum(tgt, 0.01, 20.0, 0.05, rng=rng, options=opts), tgt


def test_match_is_deterministic_and_close():
    r1, tgt = _small_match(7)
    r2, _ = _small_match(7)
    np.testing.assert_array_equal(r1.acc, r2.acc)                  # D-GEN-09 bit-identical
    r3, _ = _small_match(8)
    assert not np.array_equal(r1.acc, r3.acc)
    c = r1.checks["manual"]
    assert 0.8 < c.min_ratio and c.max_ratio < 1.4 and 0.95 < c.mean_ratio < 1.15
    v, d = SP.integrate(r1.acc, 0.01)
    assert abs(v[-1]) < 1e-10 and abs(d[-1]) < 1e-9
    assert r1.log and "selected record" in r1.log[-1]


def test_match_without_wavelet_band_uses_lw_polish():
    """No check frequency in the wavelet band (here: wavelets limited to 0.1/dt = 20 Hz, target
    25-100 Hz): the wavelet stage is skipped and logged, the LW polish runs alone."""
    tg = EL.TargetSpectrum(np.array([25.0, 50.0, 100.0]), np.array([1.5, 1.2, 1.0]), 0.05)
    res = EL.match_spectrum(tg, 0.005, 20.0, 0.05, rng=np.random.default_rng(1),
                            options=EL.MatchOptions(lw_iterations=3, wavelet_iterations=2, outer_passes=1,
                                                    wavelet_fmax=0.2))
    assert res.acc.size == 4001 and np.all(np.isfinite(res.acc))
    assert any("wavelet refinement skipped" in ln for ln in res.log)
    assert res.checks["manual"].fmin == pytest.approx(25.0)


def test_match_controls_long_period_drift():
    """SRP 3.7.1 'no baseline drift': the displacement of a matched record oscillates at the target
    periods; little of its energy lies at periods longer than the record (cf. the review test)."""
    res, tgt = _small_match(3)
    m = EL.drift_measures(res.acc, 0.01, 980.665)
    assert m["frac_long"] < 0.2
    assert m["PGD"] < 1.5 * 0.3 * 36.0 * 2.54
    assert any("drift control" in ln for ln in res.log)


def test_match_with_seed_record_keeps_character():
    f = EL.check_grid(0.3, 50.0, 20)
    tgt = EL.TargetSpectrum(f, EL.rg160_spectrum(f, "H", 0.05, 0.3), 0.05)
    rng = np.random.default_rng(1)
    dt = 0.01
    n = 2001
    t = np.arange(n) * dt
    seed = rng.standard_normal(n) * np.exp(-((t - 8.0) / 3.0) ** 2)
    res = EL.match_spectrum(tgt, dt, 20.0, 0.05, seed_acc=seed, options=EL.MatchOptions(lw_iterations=8,
                                                                                    wavelets=False))
    assert res.checks["manual"].mean_ratio == pytest.approx(1.05, abs=0.1)
    # energy stays where the seed had it (phases kept)
    t5, t75 = SP.strong_motion_window(res.acc, dt)
    assert 3.0 < t5 < 8.0 < t75 < 13.0


def test_psd_floor_raises_deficient_bands():
    dt = 0.01
    n = 2001
    rng = np.random.default_rng(3)
    t = np.arange(n) * dt
    spec = np.fft.rfft(rng.standard_normal(n) * EL.saragoni_hart_envelope(t, 20.0))
    fr = np.fft.rfftfreq(n, dt)
    spec[(fr > 2.0) & (fr < 3.0)] *= 0.1                          # a gap in the spectrum
    a = np.fft.irfft(spec, n) * 0.05

    def target(ff):
        return np.full_like(np.asarray(ff, float), 200.0)
    before = EL.psd_check(a * 981.0, dt, target, 1.0, 5.0)
    b = a
    for _ in range(3):
        b = EL.psd_floor(b, dt, target, 981.0, level=1.0, f1=1.0, f2=5.0)
    after = EL.psd_check(b * 981.0, dt, target, 1.0, 5.0)
    assert after["min_ratio"] > before["min_ratio"]
