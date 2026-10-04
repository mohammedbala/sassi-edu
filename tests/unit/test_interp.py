"""Unit tests of sassi.core.interp (requirements section 4.8 items 2-5, D-MOT-01 ... D-MOT-05)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.core import interp as I


def _data(n=12, seed=0, m=2):
    rng = np.random.default_rng(seed)
    f = np.cumsum(rng.uniform(0.2, 0.6, n)) + 0.3
    H = rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m)) + 2.0
    return f, H


# ------------------------------------------------------------------ window tables
@pytest.mark.parametrize("option,N,expected", [
    (1, 10, [0, 0, 0, 0, 4, 4, 4, 4, 5]),          # starts 1, 5; tail window = last 5 points
    (4, 12, [0, 1, 1, 1, 1, 5, 5, 5, 5, 7, 7]),    # starts 2, 6; interval 1 uses m = 1
    (5, 12, [0, 0, 2, 2, 2, 2, 6, 6, 6, 6, 7]),    # starts 3, 7; intervals 1-2 use m = 1
])
def test_single_window_partitions(option, N, expected):
    win, rule = I.window_table(option, N, np.arange(N, dtype=float))
    assert rule == "single"
    assert win[:, 0].tolist() == expected
    assert np.all(win[:, 1:] == -1)
    j = np.arange(N - 1)                           # every window contains its interval
    assert np.all((win[:, 0] <= j) & (j <= win[:, 0] + 3)) and np.all(win[:, 0] + 4 <= N - 1)


@pytest.mark.parametrize("option,rule", [(0, "weighted"), (2, "mean"), (3, "nearest3")])
def test_candidate_windows(option, rule):
    N = 10
    win, r = I.window_table(option, N, np.arange(N, dtype=float))
    assert r == rule
    # 1-based: max(1, j-3) <= m <= min(j, N-4)
    for j1 in range(1, N):
        expect = list(range(max(1, j1 - 3) - 1, min(j1, N - 4)))
        got = [m for m in win[j1 - 1] if m >= 0]
        assert got == expect


def test_effective_option_fallbacks():
    assert I.effective_option(1, 5) == 1
    assert I.effective_option(1, 4) == 6
    assert I.effective_option(6, 4) == 6
    assert I.effective_option(0, 2) == -1
    with pytest.raises(ValueError):
        I.effective_option(7, 10)


# ------------------------------------------------------------------ exactness and range ends
@pytest.mark.parametrize("option", range(7))
@pytest.mark.parametrize("smooth", [0.0, 25.0])
def test_exact_at_computed_frequencies(option, smooth):
    f, H = _data()
    fo = np.concatenate([np.linspace(0, f[-1] * 1.2, 301), f])
    Hi = I.interpolate_tf(f, H, fo, option, smooth=smooth, h0=np.array([1.0, 0.0]))
    assert np.array_equal(Hi[-len(f):], H)          # bit exact (TFI = TFU)


@pytest.mark.parametrize("option", [0, 3, 6])
def test_range_ends(option):
    f, H = _data()
    fo = np.array([0.0, 0.5 * f[0], f[-1] * 1.0001, 2 * f[-1]])
    h0 = np.array([1.0, 0.25])
    Hs = I.interpolate_tf(f, H, fo, option, h0=h0)
    assert np.allclose(Hs[0], h0) and np.allclose(Hs[1], 0.5 * (h0 + H[0]))     # seismic: linear from H(0)
    assert np.all(Hs[2:] == 0)                                                    # above f_N: 0
    Hv = I.interpolate_tf(f, H, fo, option, mode="vibration")
    assert np.allclose(Hv[:2], H[0])                                              # vibration: H_1
    assert np.all(Hv[2:] == 0)


def test_one_dimensional_and_unsorted_output():
    f, H = _data(m=1)
    fo = np.array([3.0, 0.1, 2.0, f[4]])
    a = I.interpolate_tf(f, H[:, 0], fo, 1, h0=1.0)
    b = I.interpolate_tf(f, H, np.sort(fo), 1, h0=1.0)[:, 0]
    assert a.shape == (4,)
    assert np.allclose(a, b[np.argsort(np.argsort(fo))])
    assert a[3] == H[4, 0]


def test_invalid_inputs():
    f, H = _data()
    with pytest.raises(ValueError):
        I.interpolate_tf(np.r_[f[:3], f[2], f[4:]], H, f, 1)   # duplicate frequency
    with pytest.raises(ValueError):
        I.interpolate_tf(f, H[:-1], f, 1)
    with pytest.raises(ValueError):
        I.interpolate_tf(f, H, f, 1, smooth=-1)
    with pytest.raises(ValueError):
        I.interpolate_tf(f, H, f, 1, mode="other")


def test_linear_fallback_and_spline_fallback():
    f = np.array([1.0, 2.0])
    H = np.array([1 + 1j, 3 - 1j])
    Hi = I.interpolate_tf(f, H, [1.5], 0, h0=1.0)
    assert np.allclose(Hi, [2.0])                                  # N < 3: complex linear
    f4 = np.array([1.0, 2.0, 3.0, 4.0])
    cub = lambda x: 1 + 2 * x - 0.3 * x ** 2 + 0.05j * x ** 3
    Hi = I.interpolate_tf(f4, cub(f4), [1.3, 2.7, 3.9], 2, h0=cub(0.0))   # N < 5: option 6 (exact for cubics)
    assert np.allclose(Hi, cub(np.array([1.3, 2.7, 3.9])), rtol=1e-12)


# ------------------------------------------------------------------ the Tajirian form
def _two_dof(f):
    w2 = (2 * np.pi * np.asarray(f)) ** 2
    k1, k2, m1, m2 = 900.0, 2000.0, 1.0, 1.5
    c = cfactor(0.05)
    K = np.array([[k1, -k1], [-k1, k1 + k2]]) * c
    out = []
    for x in w2:
        u = np.linalg.solve(K - x * np.diag([m1, m2]), x * np.array([m1, m2]))
        out.append(1 + u[0])
    return np.array(out)


@pytest.mark.parametrize("option", range(6))
def test_options_0_5_exact_for_2dof(option):
    f = np.array([0.6, 1.4, 2.5, 3.1, 4.2, 5.0, 6.3, 7.7, 9.0, 10.5, 12.0])
    fo = np.linspace(f[0], f[-1], 400)
    Hi = I.interpolate_tf(f, _two_dof(f), fo, option, h0=1.0)
    assert np.max(np.abs(Hi - _two_dof(fo)) / np.abs(_two_dof(fo))) < 1e-10


def test_rank_deficient_sdof_windows_are_exact():
    ks = (2 * np.pi * 4) ** 2 * cfactor(0.05)
    sd = lambda f: ks / (ks - (2 * np.pi * np.asarray(f)) ** 2)
    f = np.array([2.0, 3.0, 4.0, 5.0, 6.0])
    fit = I.fit_windows(f, sd(f)[:, None])
    assert fit["rank"].tolist() == [[4]]                         # R1 V5: rank 4
    fo = np.linspace(2, 6, 101)
    assert np.max(np.abs(I.interpolate_tf(f, sd(f), fo, 1, h0=1.0) - sd(fo)) / np.abs(sd(fo))) < 1e-12


def test_adaptive_rcond_keeps_informative_singular_values():
    # narrow windows of smooth 2-DOF data (spacing 2 df up to 24 Hz): the 1e-10 cut alone truncates a
    # small but informative singular value and loses the interpolation property (~1e-8)
    f = np.arange(4, 500, 2) * 0.048828125
    H = _two_dof(f)[:, None]
    fixed = I.fit_windows(f, H, adaptive=False)
    adapt = I.fit_windows(f, H)
    assert np.max(fixed["node_err"]) > 1e-9
    assert np.max(adapt["node_err"]) <= I.NODE_TOL
    assert np.any(adapt["rcond"] < I.RCOND)
    fo = np.linspace(f[0], f[-1], 3000)
    assert np.max(np.abs(I.interpolate_tf(f, H[:, 0], fo, 1, h0=1.0) - _two_dof(fo)) / np.abs(_two_dof(fo))) < 1e-10


def test_guard_switches_to_cubic_near_a_real_pole():
    # nearly undamped SDOF: the pole is 1e-6 off the real axis inside the window
    ks = (2 * np.pi * 4) ** 2 * (1 + 2e-6j)
    sd = lambda f: ks / (ks - (2 * np.pi * np.asarray(f)) ** 2)
    f = np.array([3.0, 3.6, 4.1, 4.7, 5.2])
    fit = I.fit_windows(f, sd(f)[:, None])
    assert bool(fit["guard"][0, 0])
    ft = np.array([3.3, 4.4])
    Hi = I.interpolate_tf(f, sd(f), ft, 1, h0=1.0)
    # 4 nearest points: 3.3 -> f[0:4]; 4.4 -> f[1:5]
    for x, pts in ((3.3, slice(0, 4)), (4.4, slice(1, 5))):
        xs, ys = f[pts], sd(f[pts])
        lag = sum(ys[a] * np.prod([(x - xs[b]) / (xs[a] - xs[b]) for b in range(4) if b != a]) for a in range(4))
        assert np.isclose(Hi[list(ft).index(x)], lag, rtol=1e-12)


def test_option0_weights_and_option3_selection():
    f, H = _data(n=12, m=1)
    fit = I.fit_windows(f, H)
    x = 0.5 * (f[5] + f[6])                      # interval j = 6 (1-based): candidates m = 3..6 (1-based)
    vals = {m: I._rational_eval(fit, np.array([m]), np.array([x]))[0, 0] for m in range(2, 6)}
    vals = {m: (I._lagrange4(f, H, np.array([m]), np.array([x]))[0, 0] if fit["guard"][m, 0] else v)
            for m, v in vals.items()}
    c = {m: 0.5 * (f[m] + f[m + 4]) for m in vals}
    h = {m: 0.5 * (f[m + 4] - f[m]) for m in vals}
    w = {m: max(I.WEIGHT_FLOOR, 1 - abs(x - c[m]) / h[m]) for m in vals}
    ref0 = sum(w[m] * vals[m] for m in vals) / sum(w.values())
    assert np.isclose(I.interpolate_tf(f, H, [x], 0)[0, 0], ref0, rtol=1e-12)
    ref2 = np.mean(list(vals.values()))
    assert np.isclose(I.interpolate_tf(f, H, [x], 2)[0, 0], ref2, rtol=1e-12)
    near = sorted(vals, key=lambda m: (abs(x - c[m]), m))[:3]
    ref3 = np.mean([vals[m] for m in near])
    assert np.isclose(I.interpolate_tf(f, H, [x], 3)[0, 0], ref3, rtol=1e-12)


# ------------------------------------------------------------------ smoothing, phase, anchors
def test_smoothing_band_filter():
    f = np.array([1.0, 2.0])
    H = np.array([[1.0 + 0j], [2.0 + 0j]])
    fo = np.array([1.5, 1.5, 1.5])
    Hi = np.array([[1.5 + 0j], [4.0 + 0j], [1.2 + 0.3j]])        # inside the band, overshoot, inside
    assert np.array_equal(I.smooth_tf(f, H, fo, Hi, 0.0), Hi)    # S = 0: identity
    out = I.smooth_tf(f, H, fo, Hi, 10.0)
    assert out[0, 0] == Hi[0, 0] and out[2, 0] == Hi[2, 0]         # |H| within [A_min, A_max]: unchanged
    L = 1.5
    r = 4.0 / 2.0 - 1.0
    assert np.isclose(out[1, 0], L + (4.0 - L) / (1 + 10 * r))
    big = I.smooth_tf(f, H, fo, Hi, 1e12)
    assert np.isclose(big[1, 0], L)


def test_phase_adjustment():
    f = np.linspace(0, 10, 201)
    H = (1 + 0.1 * f) * np.exp(-1j * 0.9 * f)                    # phase lag growing beyond -pi
    assert np.array_equal(I.phase_adjust(H, 0.0), np.abs(H).astype(complex))
    assert np.allclose(I.phase_adjust(H, 1.0), H, rtol=1e-13)
    half = I.phase_adjust(H, 0.5)
    assert np.allclose(np.unwrap(np.angle(half)), -0.45 * f, atol=1e-12)
    assert I.phase_factor(6, 100.0) == 0.0 and I.phase_factor(2, 9.0) == 0.1


def test_rigid_body_anchor():
    a = 30.0
    c, s = np.cos(np.radians(a)), np.sin(np.radians(a))
    assert np.isclose(I.rigid_body_anchor(1, 0, a), c) and np.isclose(I.rigid_body_anchor(2, 0, a), s)
    assert I.rigid_body_anchor(3, 0, a) == 0.0
    assert np.isclose(I.rigid_body_anchor(1, 1, a), -s) and np.isclose(I.rigid_body_anchor(2, 1, a), c)
    assert I.rigid_body_anchor(3, 2, a) == 1.0 and I.rigid_body_anchor(1, 2, a) == 0.0
    assert I.rigid_body_anchor(4, 0, 0.0) == 0.0
    assert I.rigid_body_anchor(1, 0, 0.0) == 1.0 and I.rigid_body_anchor(2, 0, 0.0) == 0.0


# ------------------------------------------------------------------ refined guard (D-MOT-01 trigger + pole tests)
DF4 = 1.0 / (4096 * 0.005)                                       # VP grid: NFFT 4096, dt 0.005 s


def _hyst_sdof(f0, beta):
    ks = (2 * np.pi * f0) ** 2 * cfactor(beta)
    return lambda f: ks / (ks - (2 * np.pi * np.asarray(f, dtype=float)) ** 2)


@pytest.mark.parametrize("f0,beta,spacing", [(0.5, 0.02, 8), (0.4, 0.01, 4), (0.2, 0.01, 2), (0.6, 0.01, 8),
                                             (1.0, 0.01, 8), (0.3, 0.05, 8)])
def test_guard_keeps_exact_lightly_damped_sdof_on_uniform_grids(f0, beta, spacing):
    """Uniform coarse grids: the true resonance pole fails the literal ratio test but is a damped
    (physical) pole, so the exact rational form is kept for every window option (VP-28)."""
    sd = _hyst_sdof(f0, beta)
    fs = np.arange(4, 513, spacing) * DF4
    fo = np.linspace(fs[0], fs[-1], 4001)
    fit = I.fit_windows(fs, sd(fs)[:, None])
    assert fit["trigger"].any() and not fit["guard"].any()       # literal ratio test fires, refined guard not
    for opt in range(6):
        Hi = I.interpolate_tf(fs, sd(fs), fo, opt, h0=1.0)
        assert np.max(np.abs(Hi - sd(fo)) / np.abs(sd(fo))) < 1e-10, opt


def test_window_poles_give_the_hysteretic_damping():
    """Teaching check: the fitted pole of a hysteretic SDOF is (w0/w_max)^2 c(beta), implied damping beta."""
    beta = 0.04
    sd = _hyst_sdof(4.0, beta)
    f = np.array([3.0, 3.5, 4.2, 4.8, 5.5])
    fit = I.fit_windows(f, sd(f)[:, None])
    pz = I.window_poles(fit["C"][0, 0])
    true_pole = (4.0 / 5.5) ** 2 * cfactor(beta)
    k = int(np.argmin(np.abs(pz["poles"] - true_pole)))
    assert abs(pz["poles"][k] - true_pole) < 1e-10
    assert abs(pz["damping"][k] - beta) < 1e-10


def test_guard_ignores_cancelled_doublet_of_noisy_sdof():
    """A rank-deficient SDOF window with 1e-7 noise carries a pole/zero doublet with negative damping
    near the span (seed 1, window 47).  The doublet cancels (|p - z| ~ 4e-9 << 1e-3 dist), so the
    rational form is kept and the error stays at the noise level (the cubic gave 1e-4)."""
    fn = np.unique(np.round(np.geomspace(4, 512, 64)).astype(np.int64))
    fs = fn * DF4
    sd = _hyst_sdof(4.0, 0.05)
    rng = np.random.default_rng(1)
    H = sd(fs) * (1 + 1e-7 * (rng.standard_normal(len(fs)) + 1j * rng.standard_normal(len(fs))))
    fit = I.fit_windows(fs, H[:, None])
    assert bool(fit["trigger"][47, 0]) and not fit["guard"].any()
    pz = I.window_poles(fit["C"][47, 0])
    assert np.min(pz["damping"]) < 0                              # the doublet pole is "non-physical" ...
    k = int(np.argmin(pz["damping"]))
    assert np.min(np.abs(pz["zeros"] - pz["poles"][k])) < 1e-6    # ... but cancelled by a zero
    fo = np.linspace(fs[0], fs[-1], 5001)
    for opt in range(6):
        Hi = I.interpolate_tf(fs, H, fo, opt, h0=1.0)
        assert np.max(np.abs(Hi - sd(fo)) / np.abs(sd(fo))) < 1e-5, opt


@pytest.mark.parametrize("damping_sign,guarded", [(-1.0, True), (+1.0, False)])
def test_guard_rejects_near_real_pole_with_negative_damping(damping_sign, guarded):
    """H = (x - z)/(x - p) with p 2e-4 off the real axis inside the span and no cancelling zero.  With
    negative damping (Im p < 0) the pole is spurious and the cubic is used; with the same positive
    damping (beta_p ~ 1.4e-4 > DAMPING_FLOOR) it is a genuine sharp resonance and the rational form is
    kept (exact)."""
    w_max = 2 * np.pi * 5.0
    p, z = 0.7 + damping_sign * 2e-4j, 0.82 + 0.05j
    H = lambda f: (((2 * np.pi * np.asarray(f)) / w_max) ** 2 - z) / (((2 * np.pi * np.asarray(f)) / w_max) ** 2 - p)
    f = np.array([3.0, 3.6, 4.0, 4.5, 5.0])
    fit = I.fit_windows(f, H(f)[:, None])
    assert bool(fit["trigger"][0, 0])
    assert bool(fit["guard"][0, 0]) is guarded
    fo = np.sort(np.r_[np.linspace(3.0, 5.0, 801), 5.0 * np.sqrt(0.7)])   # include the peak frequency
    Hi = I.interpolate_tf(f, H(f), fo, 1, h0=1.0)
    if guarded:
        assert np.max(np.abs(Hi)) < 0.2 * np.max(np.abs(H(fo)))  # the spike is not reproduced
    else:
        assert np.max(np.abs(Hi - H(fo)) / np.abs(H(fo))) < 1e-9


def test_guard_uses_cubic_when_the_fit_cannot_reproduce_its_points():
    """Data 0, 0, 0, 0, 1: the 4 zeros force N = 0, so no rational of the Tajirian form passes through
    the fifth point (node error 1 > NODE_FAIL_TOL); the cubic of the 4 nearest points is used."""
    f = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    H = np.array([0, 0, 0, 0, 1.0], dtype=complex)
    fit = I.fit_windows(f, H[:, None])
    assert fit["node_err"][0, 0] > I.NODE_FAIL_TOL and bool(fit["guard"][0, 0])
    ft = np.array([4.5])                                            # nearest 4 points: f[1:5]
    xs, ys = f[1:5], H[1:5]
    lag = sum(ys[a] * np.prod([(ft[0] - xs[b]) / (xs[a] - xs[b]) for b in range(4) if b != a]) for a in range(4))
    assert np.isclose(I.interpolate_tf(f, H, ft, 1)[0], lag, rtol=1e-12)


# ------------------------------------------------------------------ phase adjustment referred to f = 0
def _phase_case():
    f = np.linspace(0, 10, 201)
    return f, (1 + 0.1 * f) * np.exp(-1j * 0.9 * f)               # H(0) = 1


@pytest.mark.parametrize("rho", [0.0, 0.001, 0.09, 0.5])
def test_phase_adjustment_sign_equivariance_and_anchor(rho):
    """D-MOT-05 'reference f = 0': PA(-H) = -PA(H) and the zero-frequency value is kept."""
    f, H = _phase_case()
    pos, neg = I.phase_adjust(H, rho), I.phase_adjust(-H, rho)
    assert np.max(np.abs(neg + pos)) < 1e-12
    assert abs(neg[0] - (-1.0)) < 1e-15 and abs(pos[0] - 1.0) < 1e-15
    # differential phase scaled by rho (B.6.8) and amplitude unchanged
    dphi = np.diff(np.unwrap(np.angle(neg)))
    assert np.allclose(dphi, rho * np.diff(np.unwrap(np.angle(H))), atol=1e-12)
    assert np.allclose(np.abs(neg), np.abs(H), rtol=1e-14)


def test_phase_adjustment_zero_dc_uses_the_first_nonzero_phase():
    """H(0) = 0 (rotation / normal component): the reference phase is the limit f -> 0+."""
    f, H = _phase_case()
    G = H * (f / 10.0) * np.exp(1j * 2.5)                          # G(0) = 0, arg G(0+) -> 2.5 rad
    for rho in (0.0, 0.2):
        a, b = I.phase_adjust(G, rho), I.phase_adjust(-G, rho)
        assert a[0] == 0 and np.max(np.abs(a + b)) < 1e-12
        assert np.isclose(np.angle(a[1]), np.angle(G[1]), atol=1e-12)   # reference = first non-zero bin
    assert np.allclose(I.phase_adjust(G, 0.0), np.abs(G) * np.exp(1j * np.angle(G[1])), rtol=1e-14)


def test_phase_adjustment_columns_and_rho_one():
    f, H = _phase_case()
    M = np.column_stack([H, -2.0 * H, 0.0 * H])
    out = I.phase_adjust(M, 0.3)
    assert np.allclose(out[:, 0], I.phase_adjust(H, 0.3), rtol=0, atol=1e-15)
    assert np.allclose(out[:, 1], -2.0 * I.phase_adjust(H, 0.3), rtol=0, atol=1e-14)
    assert np.all(out[:, 2] == 0)
    assert np.array_equal(I.phase_adjust(M, 1.0), M)             # rho = 1: identity


def test_phase_factor_follows_the_scheme_actually_used():
    """N < 5: options 0-5 run the spline, on which S has no effect -> rho = 0 like option 6."""
    assert I.phase_factor(2, 9.0) == 0.1 and I.phase_factor(6, 9.0) == 0.0 and I.phase_factor(-1, 9.0) == 0.0
    f = np.array([1.0, 2.0, 3.0, 4.0])
    H = np.exp(-0.4j * f) * (1 + 0.2 * f)
    fo = np.linspace(0, 4, 81)
    out = I.interpolate_tf(f, H, fo, 2, smooth=10.0, pzadj=1, h0=1.0)
    spline = I.interpolate_tf(f, H, fo, 6, pzadj=0, h0=1.0)
    assert np.max(np.abs(np.angle(out))) == 0.0                     # zero phase (rho = 0) ...
    assert np.allclose(np.abs(out), np.abs(spline), rtol=1e-14)    # ... of the spline actually used
