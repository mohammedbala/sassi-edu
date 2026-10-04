"""Unit tests of the SOIL-NON core (sassi/core/soilnon.py): MKZ backbone, Masing damping, curve fit,
extended Masing rules, shear-column discretisation, damping matrices, linear transfer function of the
discrete model, Newmark integration, .nls side file and sublayer resolution (decisions SN-1 ... SN-12)."""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import integrate

from sassi.core import shake as SH
from sassi.core import soilnon as SN

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")


# ---------------------------------------------------------------------------------------------
# MKZ backbone and Masing damping
# ---------------------------------------------------------------------------------------------
def test_mkz_backbone_tangent_and_reference_strain():
    G0, gr, beta, s = 5.0e4, 4e-4, 1.3, 0.8
    g = np.linspace(-5e-3, 5e-3, 41)
    eps = 1e-9
    num = (SN.mkz_stress(g + eps, G0, gr, beta, s) - SN.mkz_stress(g - eps, G0, gr, beta, s)) / (2 * eps)
    np.testing.assert_allclose(SN.mkz_tangent(g[g != 0], G0, gr, beta, s), num[g != 0], rtol=1e-5)
    assert SN.mkz_stress(-1e-3, G0, gr, beta, s) == pytest.approx(-SN.mkz_stress(1e-3, G0, gr, beta, s))
    # G/G0 = 1 / (1 + beta) at gamma = gamma_r, any s; linear for beta = 0
    assert SN.mkz_modulus_ratio(0.04, 0.04, 1.3, 0.8) == pytest.approx(1.0 / 2.3)
    assert SN.mkz_modulus_ratio(1.0, 0.04, 0.0, 0.8) == 1.0
    assert SN.mkz_tangent(0.0, G0, gr, beta, s) == pytest.approx(G0)


def _masing_d_hyperbolic_r2(g, gr):
    """R2 H.2 / Darendeli (2001): D_Masing,a=1 (%) of the hyperbolic model (independent transcription)."""
    return (100.0 / math.pi) * (4.0 * (g - gr * math.log((g + gr) / gr)) / (g * g / (g + gr)) - 2.0)


@pytest.mark.parametrize("ratio", [0.01, 0.3, 1.0, 3.0, 30.0])
def test_masing_damping_hyperbolic_closed_form(ratio):
    gr = 0.0352                                        # % (Darendeli sand)
    d = SN.masing_damping(ratio * gr, gr, 1.0, 1.0)
    assert 100.0 * d == pytest.approx(_masing_d_hyperbolic_r2(ratio * gr, gr), rel=1e-9, abs=1e-12)


@pytest.mark.parametrize("beta,s", [(1.0, 0.7), (1.4, 0.9), (0.8, 0.5), (1.2, 1.0)])
@pytest.mark.parametrize("ga", [1e-5, 3e-4, 2e-3, 5e-2])
def test_masing_damping_matches_quadrature(beta, s, ga):
    gr, G0 = 5e-4, 1.0
    F = lambda x: G0 * x / (1.0 + beta * (x / gr) ** s)          # noqa: E731
    integral, _ = integrate.quad(F, 0.0, ga, epsabs=0, epsrel=1e-13, limit=200)
    ref = (2.0 / math.pi) * (2.0 * integral / (ga * F(ga)) - 1.0)
    assert SN.masing_damping(ga, gr, beta, s) == pytest.approx(ref, rel=1e-8, abs=1e-12)


def test_masing_damping_limits():
    assert SN.masing_damping(0.0, 1e-3, 1.0, 1.0) == 0.0
    assert SN.masing_damping(1e-3, 1e-3, 0.0, 1.0) == 0.0                  # linear
    d = SN.masing_damping(np.logspace(-6, 0, 40), 1e-3, 1.0, 0.9)
    assert np.all(np.diff(d) > 0)                                           # increases with strain
    assert SN.masing_damping(1e3, 1e-3, 1.0, 1.0) == pytest.approx(2.0 / math.pi, rel=2e-2)   # hyperbolic limit


# ---------------------------------------------------------------------------------------------
# Curve fit (spec 11 section 7 test 9)
# ---------------------------------------------------------------------------------------------
def test_curve_fit_recovers_synthetic_parameters():
    beta, s, gr = 1.3, 0.85, 0.04
    x = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0])
    y = 1.0 / (1.0 + beta * (x / gr) ** s)
    fit = SN.fit_mkz(x, y)
    # beta and gamma_r appear only through beta / gamma_r^s (SN-3): beta = 1 and gamma_r' = gamma_r beta^(-1/s)
    assert fit.beta == 1.0
    assert fit.s == pytest.approx(s, rel=1e-6)
    assert fit.gr_pct == pytest.approx(gr * beta ** (-1.0 / s), rel=1e-6)
    assert fit.rms < 1e-8 and fit.npts == 10
    np.testing.assert_allclose(SN.mkz_modulus_ratio(x, fit.gr_pct, fit.beta, fit.s), y, atol=1e-8)


def test_curve_fit_bounds_and_degenerate_curves():
    x = np.array([1e-4, 1e-3, 1e-2, 1e-1, 1.0])
    fit = SN.fit_mkz(x, 1.0 / (1.0 + (x / 0.05) ** 1.4))        # s > 1 is outside the admissible range
    assert fit.s == pytest.approx(SN.S_BOUNDS[1])
    lin = SN.fit_mkz(x, np.ones_like(x))
    assert lin.beta == 0.0 and lin.note
    with pytest.raises(ValueError):
        SN.fit_mkz([0.0, 0.0], [1.0, 0.5])


# ---------------------------------------------------------------------------------------------
# Extended Masing rules
# ---------------------------------------------------------------------------------------------
G0, GR, BETA, S = 100.0, 1e-3, 1.0, 0.9


def _F(x):
    return SN.mkz_stress(x, G0, GR, BETA, S)


def _path(points, n=200):
    out = [0.0]
    for a, b in zip(points[:-1], points[1:]):
        out.extend(np.linspace(a, b, n + 1)[1:])
    return np.asarray(out)


def test_extended_masing_rules_follow_the_expected_curves():
    a = 1e-3
    m = SN.MasingMKZ([G0], [GR], [BETA], [S])
    pts = [0.0, 2 * a, -a, a, 0.0, 1.5 * a, 2.5 * a, -3 * a]
    path = _path(pts)
    tau = m.drive(path)[:, 0]
    t_r0 = _F(2 * a)
    t_r1 = t_r0 + 2 * _F((-a - 2 * a) / 2)
    t_r2 = t_r1 + 2 * _F((a + a) / 2)
    t_r3 = t_r2 + 2 * _F((0.0 - a) / 2)
    t_r4 = _F(2.5 * a)

    def expected(k, g):
        seg = k // 200
        if seg == 0:                          # rule 1: virgin backbone
            return _F(g)
        if seg == 1:                          # rule 2: unloading from R0
            return t_r0 + 2 * _F((g - 2 * a) / 2)
        if seg == 2:                          # reloading from R1
            return t_r1 + 2 * _F((g + a) / 2)
        if seg == 3:                          # inner unloading from R2
            return t_r2 + 2 * _F((g - a) / 2)
        if seg == 4:                          # inner reloading from R3; beyond R2 the loop closes (rule 4)
            return t_r3 + 2 * _F(g / 2) if g <= a else t_r1 + 2 * _F((g + a) / 2)
        if seg == 5:                          # beyond the largest past strain: backbone (rule 3)
            return t_r1 + 2 * _F((g + a) / 2) if g <= 2 * a else _F(g)
        return t_r4 + 2 * _F((g - 2.5 * a) / 2) if g >= -2.5 * a else _F(g)   # back to the backbone

    exp = np.array([expected(k - 1, g) if k else 0.0 for k, g in enumerate(path)])
    np.testing.assert_allclose(tau, exp, rtol=1e-12, atol=1e-12 * G0 * a)
    assert m.top[0] == 0                      # on the backbone at the end


def test_single_large_step_closes_several_loops():
    a = 1e-3
    fine = SN.MasingMKZ([G0], [GR], [BETA], [S])
    coarse = SN.MasingMKZ([G0], [GR], [BETA], [S])
    pts = [0.0, 2 * a, -a, a, 0.0]
    fine.drive(_path(pts))
    coarse.drive(_path(pts))
    # one step from 0 to 3a passes R2 (a) and R0 (2a): two loops close, the backbone is reached
    t_fine = fine.drive(np.linspace(0.0, 3 * a, 301)[1:])[-1, 0]
    t_coarse = coarse.drive([3 * a])[-1, 0]
    assert t_coarse == pytest.approx(_F(3 * a), rel=1e-12)
    assert t_fine == pytest.approx(t_coarse, rel=1e-12)
    assert coarse.top[0] == 0


def test_tangent_is_consistent_with_the_trial_stress():
    m = SN.MasingMKZ([G0, G0], [GR, 2 * GR], [BETA, 1.4], [S, 0.7])
    m.drive([np.array([1e-3, 2e-3]), np.array([-5e-4, 1e-3])])       # one element on a branch, one too
    for g in (np.array([-1e-3, 5e-4]), np.array([0.0, 1.5e-3])):
        tr = m.trial(g)
        h = 1e-10
        num = (m.trial(g + h).tau - m.trial(g - h).tau) / (2 * h)
        np.testing.assert_allclose(tr.kt, num, rtol=1e-5)


@pytest.mark.parametrize("ga_over_gr", [0.1, 1.0, 10.0])
def test_masing_loop_secant_and_area(ga_over_gr):
    ga = ga_over_gr * GR
    m = SN.MasingMKZ([G0], [GR], [BETA], [S])
    n = 4000
    up = np.linspace(0.0, ga, n + 1)
    cycle = np.concatenate([np.linspace(ga, -ga, 2 * n + 1)[1:], np.linspace(-ga, ga, 2 * n + 1)[1:]])
    m.drive(up[1:])
    tau = m.drive(cycle)[:, 0]
    g = np.concatenate([[ga], cycle])
    t = np.concatenate([[_F(ga)], tau])
    assert tau[2 * n - 1] == pytest.approx(-_F(ga), rel=1e-12)       # the loop passes (-ga, -F(ga))
    assert tau[-1] == pytest.approx(_F(ga), rel=1e-12)               # and closes at (ga, F(ga))
    area = abs(np.trapezoid(t, g)) if hasattr(np, "trapezoid") else abs(np.trapz(t, g))
    d = area / (2 * math.pi * ga * _F(ga))
    assert d == pytest.approx(SN.masing_damping(ga, GR, BETA, S), rel=2e-6)


def test_reversal_memory_grows_beyond_the_initial_capacity():
    m = SN.MasingMKZ([G0], [GR], [BETA], [S], cap=4)
    amps = 2e-3 * 0.95 ** np.arange(60)                   # nested, decaying cycles never close
    pts = [0.0]
    for k, amp in enumerate(amps):
        pts.append(amp if k % 2 == 0 else -amp)
    m.drive(pts[1:])
    assert m.top[0] == len(amps) - 1 and m.rg.shape[1] >= m.top[0]


def test_linear_elements_keep_no_memory():
    m = SN.MasingMKZ([G0], [GR], [0.0], [S])
    tau = m.drive([1e-3, -1e-3, 2e-3, 0.0])[:, 0]
    np.testing.assert_allclose(tau, G0 * np.array([1e-3, -1e-3, 2e-3, 0.0]))
    assert m.top[0] == 0


# ---------------------------------------------------------------------------------------------
# Column, modes, damping
# ---------------------------------------------------------------------------------------------
def _column(thick=(4.0, 6.0, 10.0), rho=(1.8, 1.9, 2.0), vs=(150.0, 220.0, 350.0), xi=0.03, beta=0.0,
            eta=0.0, base=(2.2, 1200.0), fmax=25.0, npw=10):
    models = [SN.SublayerModel(beta, 1.0, 0.05 if beta else float("inf"), xi, eta) for _ in thick]
    nel = SN.elements_per_sublayer(thick, vs, fmax, npw)
    return SN.SoilNonColumn(list(thick), list(rho), list(vs), models, base[0], base[1], nel)


def test_discretisation_rule_mass_and_mid_elements():
    col = _column()
    assert np.all(col.nel % 2 == 1)
    assert np.all(col.thick / col.nel <= col.vs / (10 * 25.0) + 1e-12)
    assert np.sum(col.mass) == pytest.approx(np.sum(col.rho * col.thick))
    # the middle element of a sublayer is centred on its mid-height
    z_top = np.concatenate([[0.0], np.cumsum(col.h)])
    zc = 0.5 * (z_top[col.mid_elem] + z_top[col.mid_elem + 1])
    np.testing.assert_allclose(zc, col.depth_top[:-1] + 0.5 * col.thick)
    assert col.top_node[-1] == col.ne


def test_fixed_base_modes_of_a_uniform_column():
    H, vs = 30.0, 200.0
    col = _column(thick=(H,), rho=(2.0,), vs=(vs,), fmax=50.0, npw=20)
    w, phi = SN.fixed_base_modes(col)
    f = w / (2 * np.pi)
    exact = (2 * np.arange(1, 4) - 1) * vs / (4 * H)
    np.testing.assert_allclose(f[:3], exact, rtol=2e-3)
    np.testing.assert_allclose(phi.T @ (col.mass[:col.ne, None] * phi), np.eye(col.ne), atol=1e-9)


@pytest.mark.parametrize("bedint", [0, 1])
def test_modal_damping_gives_the_ratio_of_every_mode(bedint):
    col = _column(xi=0.04)
    damp = SN.build_damping(col, 1, bedint)
    w, phi = SN.fixed_base_modes(col)
    cf = damp.full[:col.ne, :col.ne]
    modal = phi.T @ cf @ phi
    np.testing.assert_allclose(np.diag(modal), 2 * 0.04 * w, rtol=1e-8)
    off = modal - np.diag(np.diag(modal))
    assert np.max(np.abs(off)) < 1e-8 * np.max(np.abs(modal))
    if bedint == 1:
        np.testing.assert_allclose(damp.full.sum(axis=1), 0.0, atol=1e-9 * np.max(np.abs(damp.full)))


def test_modal_damping_ratios_are_strain_energy_weighted():
    models = [SN.SublayerModel(0.0, 1.0, float("inf"), xi) for xi in (0.01, 0.05)]
    col = SN.SoilNonColumn([10.0, 10.0], [2.0, 2.0], [200.0, 200.0], models, 2.2, 1000.0, [5, 5])
    damp = SN.build_damping(col, 1, 0)
    assert np.all(damp.modal_xi > 0.01) and np.all(damp.modal_xi < 0.05)
    w, phi = SN.fixed_base_modes(col)
    dphi = np.vstack([phi[1:] - phi[:-1], -phi[-1:]])
    ek = col.k0[:, None] * dphi ** 2
    np.testing.assert_allclose(damp.modal_xi, (col.elem_xi[:, None] * ek).sum(0) / ek.sum(0))


def test_rayleigh_and_viscous_damping():
    col = _column(xi=0.02)
    d3 = SN.build_damping(col, 3, 0)                          # multipliers 0: matched at f1 and 5 f1
    w, _ = SN.fixed_base_modes(col)
    for wk in (w[0], 5 * w[0]):
        assert d3.alpha / (2 * wk) + d3.beta * wk / 2 == pytest.approx(0.02, rel=1e-9)
    assert d3.notes
    d3u = SN.build_damping(col, 3, 1, mmmult=0.5, smmult=0.001)
    assert d3u.alpha == 0.5 and d3u.beta == 0.001 and d3u.kind == "tri" and len(d3u.diag) == col.ne + 1
    col2 = _column(eta=20.0)
    d2 = SN.build_damping(col2, 2, 0)
    v = np.zeros(col2.ne)
    v[0] = 1.0
    cv = d2.matvec(v)
    assert cv[0] == pytest.approx(20.0 / col2.h[0]) and cv[1] == pytest.approx(-20.0 / col2.h[0])
    d0 = SN.build_damping(col, 0, 0)
    assert d0.type == 1 and "not a documented choice" in d0.notes[0]


# ---------------------------------------------------------------------------------------------
# Linear transfer of the discrete model vs the continuum
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("bedint", [0, 1])
def test_discrete_transfer_function_converges_to_the_continuum(bedint):
    H, vs, rho, eta = 20.0, 200.0, 2.0, 100.0
    col = _column(thick=(H,), rho=(rho,), vs=(vs,), eta=eta, fmax=50.0, npw=40, base=(2.2, 800.0))
    damp = SN.build_damping(col, 2, bedint)
    f = np.array([0.5, 1.0, 2.5, 4.0, 6.0])
    Hd = SN.linear_transfer(col, damp, bedint, f, nodes=[0])[:, 0]
    w = 2 * np.pi * f
    gs = rho * vs ** 2 + 1j * w * eta                          # Kelvin-Voigt complex modulus
    vstar = np.sqrt(gs / rho)
    kH = w / vstar * H
    if bedint == 0:
        Hc = 1.0 / np.cos(kH)                                   # surface / rigid base
    else:
        a = rho * vstar / (2.2 * 800.0)
        Hc = 1.0 / (np.cos(kH) + 1j * a * np.sin(kH))           # surface / outcrop
    np.testing.assert_allclose(Hd, Hc, rtol=3e-3)
    assert SN.linear_transfer(col, damp, bedint, [0.0])[0, 0] == 1.0


# ---------------------------------------------------------------------------------------------
# Time integration
# ---------------------------------------------------------------------------------------------
def _record(n=600, dt=0.01, fmax=12.0, seed=2):
    rng = np.random.default_rng(seed)
    a = rng.standard_normal(n) * np.sin(np.pi * np.arange(n) / n) ** 2
    A = np.fft.rfft(a, 4 * n)
    A[np.fft.rfftfreq(4 * n, dt) > fmax] = 0
    a = np.fft.irfft(A, 4 * n)[:n]
    return a / np.max(np.abs(a))


def test_band_limited_input_interpolates_the_samples():
    a = _record(200)
    a_dt, a_fine = SN.band_limited_input(a, 0.01, 512, up=8)
    np.testing.assert_allclose(a_fine[::8], a_dt, atol=1e-12)
    np.testing.assert_allclose(a_dt[:200], a, atol=1e-12)
    b_dt, _ = SN.band_limited_input(a, 0.01, 512, cutoff=5.0)
    B = np.fft.rfft(b_dt)
    assert np.max(np.abs(B[np.fft.rfftfreq(512, 0.01) > 5.0 + 1e-9])) < 1e-9
    with pytest.raises(ValueError):
        SN.band_limited_input(a, 0.01, 128)


@pytest.mark.parametrize("damptype,bedint", [(1, 1), (2, 0), (3, 1)])
def test_linear_time_integration_matches_the_discrete_transfer_function(damptype, bedint):
    col = _column(xi=0.05, eta=150.0)
    damp = SN.build_damping(col, damptype, bedint)
    dt, nfft = 0.01, 2048
    a = _record()
    errs = []
    for nsub in (4, 8):
        ctl = SN.Controls(nsub=nsub, flexible=False, dg_max=5e-5, tol_d=1e-9, tol_f=1e-9, maxit=10)
        res = SN.run_soilnon(a, dt, nfft, col, damp, bedint, ctl)
        f = np.fft.rfftfreq(nfft, dt)
        Hs = SN.linear_transfer(col, damp, bedint, f, nodes=[0])[:, 0]
        ref = np.fft.irfft(np.fft.rfft(res.input_motion) * Hs, nfft)
        td = res.acc_top[:, 0]
        assert np.max(np.abs(td)) == pytest.approx(np.max(np.abs(ref)), rel=3e-3)
        errs.append(np.sqrt(np.mean((td - ref) ** 2) / np.mean(ref ** 2)))
    assert errs[1] < 0.35 * errs[0]                    # second-order accurate (ratio 1/4)
    assert res.stats["bisections"] == 0


def test_nonlinear_run_is_in_equilibrium_and_follows_the_masing_model():
    """Without viscous damping the total shear stress (inertia of the soil above, SN-12) equals the
    hysteretic stress of the element at every output step."""
    col = _column(thick=(5.0, 5.0), rho=(1.8, 1.9), vs=(120.0, 180.0), beta=1.0, eta=0.0)
    damp = SN.build_damping(col, 2, 0)
    a = 3.0 * _record(800)
    ctl = SN.Controls(nsub=2, flexible=True, dg_max=SN.FLEX_STRAIN_INCREMENT, tol_d=1e-9, tol_f=1e-9, maxit=25)
    res = SN.run_soilnon(a, 0.01, 2048, col, damp, 0, ctl)
    smax = np.max(np.abs(res.tau_h_mid))
    np.testing.assert_allclose(res.stress_mid, res.tau_h_mid, atol=1e-6 * smax)
    assert np.max(res.gmax_mid) > 0.05e-2                      # beyond gamma_r = 0.05 %: nonlinear range
    assert res.stats["strain_rejections"] > 0                  # the flexible step was refined
    # absolute acceleration at the rigid base is the input motion
    np.testing.assert_allclose(res.acc_top[:, -1], res.input_motion, atol=1e-12)


def test_convergence_failure_is_reported():
    col = _column(thick=(5.0,), rho=(1.8,), vs=(120.0,), beta=1.0)
    damp = SN.build_damping(col, 2, 0)
    ctl = SN.Controls(nsub=1, flexible=False, dg_max=5e-5, tol_d=1e-30, tol_f=1e-30, maxit=1)
    with pytest.raises(SN.SoilNonError, match="equilibrium not reached"):
        SN.run_soilnon(3.0 * _record(100), 0.01, 256, col, damp, 0, ctl)


def test_controls_from_options():
    c0 = SN.Controls.from_options(SN.NLSoilOptions(opt=1), 0.02)
    assert c0.flexible and c0.nsub == 8 and c0.tol_d == SN.DEFAULT_DISPCONV and c0.maxit == SN.DEFAULT_EQUALIT
    c1 = SN.Controls.from_options(SN.NLSoilOptions(opt=1, nsub=3, dispconv=1e-4, forceconv=1e-5, equalit=7), 0.01)
    assert not c1.flexible and c1.nsub == 3 and c1.tol_d == 1e-4 and c1.tol_f == 1e-5 and c1.maxit == 7


# ---------------------------------------------------------------------------------------------
# .nls side file, NLSOIL / NLSLAYER data, sublayer resolution
# ---------------------------------------------------------------------------------------------
def test_nls_text_round_trip(tmp_path):
    o = SN.NLSoilOptions(opt=1, nsub=4, dispconv=1e-5, forceconv=2e-5, equalit=12, bedint=1, damptype=3,
                         mmmult=0.25, smmult=0.001)
    lays = [SN.NLLayerData(2, 0, 1.2, 0.85, 0.05, 0.0), SN.NLLayerData(1, 1)]
    p = SN.write_nls(tmp_path / "m.nls", SN.NlsData(o, {l.num: l for l in lays}, "m", "abc"))
    back = SN.read_nls(p)
    assert back.options == o and back.model == "m" and back.model_hash == "abc"
    assert back.layers[2] == lays[0] and back.layers[1] == lays[1]
    assert p.read_text().splitlines()[3] == "NLSOIL,1,4,1e-05,2e-05,12,1,3,0.25,0.001"
    with pytest.raises(ValueError, match="no NLSOIL"):
        SN.parse_nls("NLSLAYER,1,1\n")
    with pytest.raises(ValueError, match="unknown record"):
        SN.parse_nls("NLSOIL,1\nFOO,1\n")


def test_option_and_layer_problems():
    assert not SN.NLSoilOptions(opt=1, damptype=2).problems()
    bad = SN.NLSoilOptions(opt=2, nsub=-1, bedint=3, damptype=7, mmmult=-1.0)
    assert len(bad.problems()) == 5
    assert not SN.NLLayerData(1, 0, 1.0, 0.9, 0.03).problems()
    assert not SN.NLLayerData(1, 0, 0.0).problems()                      # beta 0: linear (SN-2)
    assert SN.NLLayerData(1, 0, 1.0, 0.0, 0.03).problems()               # S missing
    assert SN.NLLayerData(0, 1).problems()


def _curve(label="Sand"):
    x = np.array([1e-4, 1e-3, 1e-2, 0.1, 1.0])
    return SH.DynamicProperty(label, x, 1.0 / (1.0 + (x / 0.05) ** 0.9), x, np.array([0.5, 1.0, 3.0, 9.0, 20.0]))


def test_resolve_sublayer():
    m = SN.resolve_sublayer(1, SN.NLLayerData(1, 1, 9.0, 9.0, 9.0, 9.0), _curve(), 0.05)
    assert m.source == "fit" and m.beta == 1.0 and m.s == pytest.approx(0.9, rel=1e-6)
    assert m.xi_min == pytest.approx(0.005) and "DYNP" in m.xi_source             # SN-4: D at 1e-4 %
    u = SN.resolve_sublayer(2, SN.NLLayerData(2, 0, 1.2, 0.8, 0.04, 30.0), None, 0.03)
    assert (u.source, u.beta, u.s, u.gr_pct, u.eta, u.xi_min) == ("user", 1.2, 0.8, 0.04, 30.0, 0.03)
    lin = SN.resolve_sublayer(3, SN.NLLayerData(3), None, 0.02)
    assert lin.linear and lin.source == "linear"
    with pytest.raises(SN.SoilNonError, match="Error 125"):
        SN.resolve_sublayer(4, None, None, 0.02)
    with pytest.raises(SN.SoilNonError, match="no DYNP"):
        SN.resolve_sublayer(5, SN.NLLayerData(5, 1), None, 0.02)
    gg, dd = SN.equivalent_properties(m, 0.05)
    assert gg == pytest.approx(SN.mkz_modulus_ratio(0.05, m.gr_pct, 1.0, m.s))
    assert dd == pytest.approx(m.xi_min + SN.masing_damping(0.05, m.gr_pct, 1.0, m.s))
    assert SN.equivalent_properties(lin, 1.0) == (1.0, 0.02)
