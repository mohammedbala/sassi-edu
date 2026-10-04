"""Unit tests of the Option NON hysteresis models and equivalent-linear properties (sassi/core/hysteresis.py;
requirements 4.15, D-NON-03, D-NON-04, D-NON-06, VP-45 basis)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.core import hysteresis as H


def epp(xy=0.01, fy=100.0, xend=1.0):
    return H.Backbone([xy, xend], [fy, fy], 1)


def smooth(dc=1e-4, vc=1000.0, vu=2500.0, gy=0.004, n=20):
    g = [dc] + [dc + k * (gy - dc) / n for k in range(1, n + 1)] + [0.02]
    v = [vc] + [vc + (vu - vc) * (1 - (1 - (x - dc) / (gy - dc)) ** 2) for x in g[1:n + 1]] + [1.02 * vu]
    return H.Backbone(g, v, n + 1)


# ======================================================================================
# Backbone
# ======================================================================================
def test_backbone_evaluation_and_integral():
    bb = H.Backbone([1.0, 3.0, 5.0], [10.0, 20.0, 25.0], 2)
    assert bb.k_el == 10.0 and (bb.x_cr, bb.f_cr, bb.x_y, bb.f_y) == (1.0, 10.0, 3.0, 20.0)
    assert bb.force(0.5) == pytest.approx(5.0) and bb.force(-2.0) == pytest.approx(-15.0)
    assert bb.force(10.0) == pytest.approx(25.0)                    # constant beyond the last point
    np.testing.assert_allclose(bb.force(np.array([-4.0, 0.0, 4.0])), [-22.5, 0.0, 22.5])
    # exact integral of the polyline: 5 + (10+20)/2*2 + (20+22.5)/2*1 = 5 + 30 + 21.25
    assert bb.integral(4.0) == pytest.approx(56.25)
    assert bb.integral(6.0) == pytest.approx(5 + 30 + 45 + 25)
    assert bb.secant(0.5) == 10.0 and bb.secant(4.0) == pytest.approx(22.5 / 4.0)


@pytest.mark.parametrize("x,y,yi,msg", [
    ([1.0], [1.0], 1, "at least 2 points"),
    ([1.0, 1.0], [1.0, 2.0], 1, "strictly increasing"),
    ([0.0, 1.0], [1.0, 2.0], 1, "> 0"),
    ([1.0, 2.0], [2.0, 1.0], 1, "must not decrease"),
    ([1.0, 2.0], [1.0, 2.0], 3, "yield point number"),
])
def test_backbone_problems(x, y, yi, msg):
    with pytest.raises(H.HysteresisError, match=msg):
        H.Backbone(x, y, yi)


def test_strict_models_reject_flat_segments():
    with pytest.raises(H.HysteresisError, match="increase strictly"):
        H.ChengMertzShear(epp())
    with pytest.raises(H.HysteresisError, match="increase strictly"):
        H.Takeda(epp())
    with pytest.raises(H.HysteresisError, match="not included"):
        H.make_model(2, smooth())


# ======================================================================================
# GMR (Masing): closed form, memory rules
# ======================================================================================
@pytest.mark.parametrize("xa", [0.015, 0.02, 0.05, 0.1, 0.5])
def test_gmr_epp_loop_closed_form(xa):
    """VP-45 basis: elastic-perfectly-plastic BBC cycled at x > x_y gives a closed loop with
    xi_h = 2 (x - x_y)/(pi x) and K_sec = F_y/x."""
    bb = epp()
    lp = H.cyclic_loop(H.MasingGMR(bb), xa)
    assert lp.closed == 0.0
    assert lp.f_pos == pytest.approx(100.0, rel=1e-14) and lp.f_neg == pytest.approx(-100.0, rel=1e-14)
    assert lp.damping == pytest.approx(2 * (xa - 0.01) / (math.pi * xa), rel=1e-12)
    assert lp.secant == pytest.approx(100.0 / xa, rel=1e-14)
    assert lp.energy == pytest.approx(4 * 100.0 * (xa - 0.01), rel=1e-12)
    assert lp.energy == pytest.approx(H.masing_loop_energy(bb, xa), rel=1e-12)


def test_gmr_smooth_backbone_matches_masing_energy():
    bb = smooth()
    for xa in (5e-4, 2e-3, 3.5e-3, 0.01):
        lp = H.cyclic_loop(H.MasingGMR(bb), xa)
        assert lp.closed < 1e-9 * bb.f_y
        assert lp.energy == pytest.approx(H.masing_loop_energy(bb, xa), rel=1e-10)
        assert lp.secant == pytest.approx(bb.secant(xa), rel=1e-12)


def test_gmr_masing_branch_and_memory_rule():
    bb = epp(xy=1.0, fy=10.0, xend=100.0)
    m = H.MasingGMR(bb)
    assert m.step(3.0) == pytest.approx(10.0)                        # virgin: backbone
    # unloading branch F = F_r + 2 F_bb((x - x_r)/2): elastic over 2 x_y
    assert m.step(2.0) == pytest.approx(0.0)
    assert m.step(1.0) == pytest.approx(-10.0)
    assert m.step(0.0) == pytest.approx(-10.0)
    # inner loop 0 -> 1.5 -> 0: reloading branch from (0, -10), then the memory rule at the reversal point
    assert m.step(1.5) == pytest.approx(-10.0 + 2 * 10.0 * min(0.75, 1.0))
    f = m.step(0.0)
    assert f == pytest.approx(-10.0)                                 # loop closed at (0, -10)
    # continuing down: the branch from (3, 10) continues (memory), then the backbone beyond -3
    assert m.step(-3.0) == pytest.approx(-10.0)
    assert m.stack == []                                             # back on the backbone
    assert m.step(-5.0) == pytest.approx(-10.0)


def test_gmr_vertices_record_exact_kinks():
    bb = H.Backbone([1.0, 2.0], [10.0, 15.0], 1)
    m = H.MasingGMR(bb)
    m.step(4.0)
    m.step(-4.0)
    xs = [v[0] for v in m.vertices]
    assert 1.0 in xs and 2.0 in xs                                   # backbone kinks on the virgin curve
    assert 2.0 in xs and 0.0 in xs                                   # branch kinks at x_r - 2 x_k


# ======================================================================================
# CMS (Cheng-Mertz shear, HYST04)
# ======================================================================================
def test_cms_statement_functions():
    bb = smooth()
    m = H.ChengMertzShear(bb)
    X = 10 * bb.x_cr
    assert m.SI == pytest.approx(bb.k_el)
    assert m.FS1(X) == pytest.approx(m.SI * 1.4675 * 0.1 ** 0.345)
    assert m.FS2(X) == pytest.approx(m.SI * 0.7761 * 0.1 ** 0.5195)
    assert m.FS3(X) == pytest.approx(m.SI * 0.0707 * 0.1 ** 1.369)
    assert m.FSR(X) == pytest.approx(m.SI * 0.1 ** 1.02)
    assert m.FS1(bb.x_cr) == m.SI                                    # capped at the initial stiffness


def test_cms_virgin_loading_follows_backbone_and_elastic_unloading():
    bb = smooth()
    m = H.ChengMertzShear(bb)
    for x in (0.5e-4, 1e-4, 1e-3, 3e-3):
        assert m.step(x) == pytest.approx(bb.force(x), rel=1e-12)
    m.reset()
    m.step(0.8e-4)
    assert m.step(-0.5e-4) == pytest.approx(-0.5e-4 * bb.k_el)       # elastic before the first cracking


def test_cms_unloading_bands():
    """Rules 2-4: from the peak (DM, PM) the unloading slope is S1 down to PM - PC, then max(S2, SO2) down
    to PC/2, then max(S3, SO3) to zero force."""
    bb = smooth()
    m = H.ChengMertzShear(bb)
    pm = m.step(3e-3)
    s1 = m.FS1(3e-3)
    f = m.step(3e-3 - 0.5 * (bb.f_cr / s1))          # half way down the first band
    assert f == pytest.approx(pm - 0.5 * bb.f_cr, rel=1e-12)
    assert m.rule == 2.0
    m.step(-1e-4)
    xs = np.array(m.vertices)
    zero = xs[np.isclose(xs[:, 1], 0.0, atol=1e-9)]
    assert zero.size                                                 # passes exactly through zero force (rule 4 limit)


def test_cms_cyclic_degradation_and_pinching():
    """Repeated cycles at the same amplitude reach the backbone only at 1.04 DMAX (ALPHA): the peak force
    drops below the backbone and stabilises; a low cracking point gives pinched (S-shaped) loops."""
    bb = smooth()
    lp = H.cyclic_loop(H.ChengMertzShear(bb), 3e-3, cycles=4)
    assert lp.f_pos < bb.force(3e-3) and lp.f_pos > 0.95 * bb.force(3e-3)
    assert lp.closed < 1e-6 * bb.f_y                                 # stabilised cycle
    assert 0.05 < lp.damping < 0.3
    # pinching: with a low cracking point the reloading below 0.75 PC is softer than with a high one
    lo = H.Backbone(*_fig13_bbc(2e-5, 795.0))
    hi = H.Backbone(*_fig13_bbc(6.2e-5, 2470.0))
    assert H.hysteretic_damping(1, lo, 8e-4) < H.hysteretic_damping(4, lo, 8e-4)


def _fig13_bbc(dc, vc, vu=5000.0, gpk=3e-4, gend=1.2e-3, n=20):
    g = list(np.linspace(dc, gpk, n + 1))
    v = [vc + (vu - vc) * (1 - (1 - (x - dc) / (gpk - dc)) ** 2) for x in g]
    return g + [gend], v + [vu * 1.002], n + 1


def test_cms_small_loops_and_memory():
    bb = smooth()
    m = H.ChengMertzShear(bb)
    for x in (3e-3, 1e-3, 2e-3, -3e-3, 3.5e-3):
        m.step(x)
    assert m.f > 0 and m.rule in (1.0, 10.0)
    assert m.failures == 0 and not m.warnings


def test_cms_is_deterministic_and_rate_independent():
    bb = smooth()
    hist = 3e-3 * np.sin(np.linspace(0, 6 * np.pi, 400)) * np.linspace(0.3, 1.0, 400)
    f1 = H.ChengMertzShear(bb).run(hist)
    f2 = H.ChengMertzShear(bb).run(hist)
    np.testing.assert_array_equal(f1, f2)
    # a finer sampling of the same path gives the same forces at the common samples
    fine = np.interp(np.linspace(0, 399, 3991), np.arange(400), hist)
    f3 = H.ChengMertzShear(bb).run(fine)
    np.testing.assert_allclose(f3[::10], f1, rtol=1e-9, atol=1e-9 * bb.f_y)


# ======================================================================================
# TAK (Takeda)
# ======================================================================================
def test_takeda_rules():
    bb = H.Backbone([1.0, 4.0, 10.0], [100.0, 160.0, 170.0], 2)   # Dc 1, Pc 100, Dy 4, Py 160
    m = H.Takeda(bb)
    assert (m.K0, m.K1, m.Kcy) == (100.0, 20.0, 52.0)
    assert m.step(0.5) == 50.0 and m.step(-0.5) == -50.0              # elastic both ways
    assert m.step(3.0) == pytest.approx(140.0)                        # primary curve after cracking
    # unloading before yield toward the opposite cracking point (-1, -100): slope (140 + 100)/(3 + 1)
    assert m.step(1.0) == pytest.approx(140.0 - 60.0 * 2.0)
    assert m.step(-1.0) == pytest.approx(-100.0)                      # reaches (-Dc, -Pc) exactly
    assert m.step(-2.0) == pytest.approx(-120.0) and m.mode == "primary"
    m.step(6.0)                                                       # reload to the +peak, yield, flat
    assert m.f == pytest.approx(160.0) and m.yielded[1]
    # unloading after yield: K_u = K_cy (Dy/Dmax)^0.4
    ku = 52.0 * (4.0 / 6.0) ** 0.4
    assert m.step(5.0) == pytest.approx(160.0 - ku)
    # partial unloading reversed: reload toward the unloading point (6, 160), then the envelope
    assert m.step(6.0) == pytest.approx(160.0) and m.step(7.0) == pytest.approx(160.0)


def test_takeda_zero_crossing_targets_opposite_peak():
    bb = H.Backbone([1.0, 4.0, 10.0], [100.0, 160.0, 170.0], 2)
    m = H.Takeda(bb)
    m.step(6.0)
    m.step(-6.0)
    assert m.f == pytest.approx(-160.0) and m.yielded[-1]
    m.step(0.0)
    ku = 52.0 * (4.0 / 6.0) ** 0.4
    x0 = -6.0 + 160.0 / ku                                            # zero crossing of the unloading line
    # reloading line from (x0, 0) to the positive peak (6, 160)
    assert m.step(2.0) == pytest.approx(160.0 * (2.0 - x0) / (6.0 - x0))


def test_takeda_backbone_is_capped_at_yield():
    bb = smooth()
    m = H.Takeda(bb)
    assert m.backbone_force(0.01) == pytest.approx(bb.f_y)            # peak capacity = yield force (manual)
    assert H.secant_ratio(3, bb, 0.01) == pytest.approx(bb.f_y / 0.01 / bb.k_el)


# ======================================================================================
# Equivalent linearisation, damping rule, convergence
# ======================================================================================
def test_combine_damping_rule():
    """D-NON-04: xi = min(cutoff, scale xi_h + [xi_el]); scale 0 -> 1, cutoff 0 -> none; cutoff in %."""
    assert H.combine_damping(0.10, 0.04) == pytest.approx(0.10)
    assert H.combine_damping(0.10, 0.04, include_elastic=True) == pytest.approx(0.14)
    assert H.combine_damping(0.10, 0.04, cutoff_pct=7.0, include_elastic=True) == pytest.approx(0.07)
    assert H.combine_damping(0.10, 0.04, scale=0.5, include_elastic=True) == pytest.approx(0.09)
    assert H.combine_damping(0.10, 0.04, scale=0.0) == pytest.approx(0.10)


def test_equivalent_linear_harmonic_epp():
    bb = epp()
    t = np.linspace(0, 4, 4001)
    x = 0.03 * np.sin(2 * np.pi * t)
    r = H.equivalent_linear(4, bb, x, edf=0.8, xi_el=0.02, include_elastic=True)
    xeq = 0.8 * 0.03
    assert r.x_max == pytest.approx(0.03) and r.x_eq == pytest.approx(xeq)
    assert r.ratio == pytest.approx(100.0 / xeq / bb.k_el, rel=1e-12)
    assert r.xi_h == pytest.approx(2 * (xeq - 0.01) / (math.pi * xeq), rel=1e-9)
    assert r.xi == pytest.approx(r.xi_h + 0.02)
    assert r.ductility == pytest.approx(3.0) and r.f_mu == pytest.approx(10000.0 * 0.03 / 100.0, rel=1e-9)
    # elastic response: ratio 1, no hysteretic damping, F_mu 1
    e = H.equivalent_linear(4, bb, 0.5e-2 * np.sin(2 * np.pi * t), 0.8, 0.02)
    assert (e.ratio, e.xi_h, e.f_mu) == (1.0, 0.0, pytest.approx(1.0))


def test_takeda_uses_elastic_damping():
    bb = smooth()
    r = H.equivalent_linear(3, bb, 3e-3 * np.sin(np.linspace(0, 8 * np.pi, 800)), 0.8, 0.04, include_elastic=False)
    assert r.xi == 0.04 and any("TAK" in n for n in r.notes)


def test_convergence_measure():
    c = H.convergence([100.0, 200.0], [101.0, 195.0], [0.05, 0.06], [0.052, 0.06])
    assert c.max_de == pytest.approx(0.025) and c.max_dxi == pytest.approx(0.002) and not c.converged
    c = H.convergence([100.0], [101.5], [0.05], [0.054])
    assert c.converged


def test_crv_table():
    bb = smooth()
    a = H.crv_amplitudes(bb, 20)
    assert a[0] == pytest.approx(bb.x_cr / 10) and np.all(np.diff(a) > 0) and set(bb.x) <= set(a)
    tab = H.eql_curve(1, bb, a)
    assert tab.shape == (len(a), 5)
    below = tab[:, 0] <= bb.x_cr
    assert np.all(tab[below, 1] == 1.0) and np.all(tab[below, 2] == 0.0)
    assert np.all(np.diff(tab[~below, 1]) < 0)                         # secant ratio decreases with amplitude
