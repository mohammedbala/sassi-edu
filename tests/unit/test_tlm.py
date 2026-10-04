"""Unit tests of the thin-layer core (sassi.core.tlm, greens_tlm, freefield half-space stiffness).

UT-19 (half-space buffer) is :func:`test_ut19_halfspace_buffer`.
"""
from __future__ import annotations

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.core import freefield as ff
from sassi.core import greens_tlm, tlm


def uniform_column(n, h, vs, vp, rho=2.0, beta=0.0, base="rigid"):
    G, M = tlm.layer_moduli(rho, vs, vp, beta, beta)
    return tlm.Column(h=[h] * n, rho=[rho] * n, G=[complex(G)] * n, M=[complex(M)] * n, base=base)


def layered_column():
    layers = [(1.0, 200, 400, 2.0, 0.05, 0.05)] * 5 + [(2.0, 400, 800, 2.1, 0.03, 0.03)] * 5
    h = [l[0] for l in layers]
    rho = [l[3] for l in layers]
    G, M = tlm.layer_moduli(rho, [l[1] for l in layers], [l[2] for l in layers], [l[4] for l in layers],
                            [l[5] for l in layers])
    return tlm.Column(h=h, rho=rho, G=G, M=M)


# ------------------------------------------------------------------ layer matrices
def test_single_layer_matrices_match_kausel_table():
    h, rho, G, lam = 2.0, 1.8, 1000.0 + 50j, 1500.0 + 70j
    col = tlm.Column(h=[h], rho=[rho], G=[G], M=[lam + 2 * G], base="dashpot")
    m = tlm.column_matrices(col, tlm.MASS_CONSISTENT)
    lp2 = lam + 2 * G
    np.testing.assert_allclose(m.Ax, h / 6 * lp2 * np.array([[2, 1], [1, 2]]))
    np.testing.assert_allclose(m.Az, h / 6 * G * np.array([[2, 1], [1, 2]]))
    np.testing.assert_allclose(m.Gs, G / h * np.array([[1, -1], [-1, 1]]))
    np.testing.assert_allclose(m.Gp, lp2 / h * np.array([[1, -1], [-1, 1]]))
    np.testing.assert_allclose(m.Bxz, 0.5 * np.array([[lam - G, -(lam + G)], [lam + G, -(lam - G)]]))
    np.testing.assert_allclose(m.Mm, rho * h / 6 * np.array([[2, 1], [1, 2]]))
    mixed = tlm.column_matrices(col).Mm
    np.testing.assert_allclose(mixed, rho * h * np.array([[5, 1], [1, 5]]) / 12)


def test_inplane_stiffness_is_complex_symmetric():
    m = tlm.column_matrices(layered_column()).free()
    K = m.K_inplane(2 * np.pi * 7.0, 0.4 - 0.05j)
    np.testing.assert_allclose(K, K.T, atol=1e-9 * np.abs(K).max())


# ------------------------------------------------------------------ eigen-solutions
@pytest.mark.parametrize("mass", ["consistent", "mixed"])
@pytest.mark.parametrize("k", [0.3, 1.7 + 0.2j])
def test_modal_flexibility_equals_direct_inverse(mass, k):
    """Kausel Eq. 49-52 vs direct inverse (R1 V1)."""
    col = layered_column()
    om = 2 * np.pi * 10.0
    md = tlm.column_modes(col, om, mass)
    Fd, Fyd = tlm.direct_flexibility(col, om, k, mass)
    Fm, Fym = tlm.modal_flexibility(md, k)
    assert np.abs(Fd - Fm).max() / np.abs(Fd).max() < 1e-10
    assert np.abs(Fyd - Fym).max() / np.abs(Fyd).max() < 1e-10


def test_mode_orthonormality():
    col = layered_column()
    om = 2 * np.pi * 8.0
    m = tlm.column_matrices(col).free()
    md = tlm.column_modes(col, om)
    # Love: phi^T Ay phi = I
    P = md.phiy.T @ m.Ay @ md.phiy
    np.testing.assert_allclose(P, np.eye(P.shape[0]), atol=1e-8)
    # Rayleigh: Y_i^T Abar Z_j = k_j delta_ij (Kausel Eq. 21-22a)
    n = m.n_free
    Abar = np.block([[m.Ax, np.zeros((n, n))], [m.Bxz.T, m.Az]])
    Y = np.vstack([md.phix * md.kR, md.phiz])
    Z = np.vstack([md.phix, md.phiz * md.kR])
    with np.errstate(all="ignore"):                 # spurious Accelerate FPE flags in matmul
        O = Y.T @ Abar @ Z
    np.testing.assert_allclose(O, np.diag(md.kR), atol=1e-7 * np.abs(md.kR).max())


def test_outgoing_root_branch():
    k = tlm.outgoing_root(np.array([4.0, -4.0, 4.0 - 0.1j, -4.0 + 0.1j, 4.0 + 0.1j]))
    assert k[0] == 2.0
    assert k[1] == -2.0j
    assert np.all(k.imag <= 0)
    assert k[2].real > 0


def test_rayleigh_speed_uniform_stratum_r1_v2():
    """R1 V2: c_R/Vs = 0.9210 (consistent), 0.9194 (lumped), 0.9202 (mixed) at h = lambda/20."""
    vs, nu, f = 100.0, 0.25, 20.0
    vp = vs * np.sqrt(2 * (1 - nu) / (1 - 2 * nu))
    lam = vs / f
    col = uniform_column(80, lam / 20, vs, vp, rho=1.0)
    om = 2 * np.pi * f
    for mass, ref in (("consistent", 0.9210), ("lumped", 0.9194), ("mixed", 0.9202)):
        md = tlm.column_modes(col, om, mass)
        j = tlm.select_mode(md.kR, 1)
        assert om / md.kR[j].real / vs == pytest.approx(ref, abs=1e-4)


def test_love_consistent_mass_converges_from_above():
    """R2 D.5: thin-layer Love velocities approach the exact value from above (consistent mass)."""
    b1, b2, H, rho, f = 200.0, 400.0, 10.0, 2.0, 20.0
    G1, M1 = tlm.layer_moduli(rho, b1, 2 * b1, 0, 0)
    G2, M2 = tlm.layer_moduli(rho, b2, 2 * b2, 0, 0)
    hg, _ = tlm.halfspace_sublayers(b2, f, 20, law="uniform")
    cs = []
    for n in (10, 20, 40):
        col = tlm.build_column([H / n] * n, [rho] * n, [complex(G1)] * n, [complex(M1)] * n, rho, complex(G2),
                               complex(M2), hg)
        md = tlm.column_modes(col, 2 * np.pi * f, tlm.MASS_CONSISTENT)
        k = md.kL
        w = 2 * np.pi * f
        ok = (k.real > w / b2) & (np.abs(k.imag) < 0.05 * k.real)
        cs.append(np.sort(w / k[ok].real)[1])           # first higher mode
    ref = 279.222
    assert cs[0] > cs[1] > cs[2] > ref


def test_select_mode_rules():
    k = np.array([5.0 - 0.01j, 3.0 - 0.001j, 9.0 - 8.0j, -1.0 - 2.0j])
    assert tlm.select_mode(k, 1) == 0             # largest Re k among propagating modes
    assert tlm.select_mode(k, 2) == 1             # smallest |Im k| with Re k > 0


def test_least_decay_ties_break_to_largest_re_k():
    """D-SIT-09 / D-GEN-09: |Im k| equal within TIE_RTOL max|k| is a tie, broken by the largest Re k."""
    k = np.array([3.0 + 0j, 5.0 - 1e-15j, 4.0 + 0j, 1.0 - 5.0j])
    c = tlm.choose_mode(k, 2)
    assert (c.index, c.n_tied, c.n_candidates) == (1, 3, 4)
    # a genuine (physical) difference of the decay is still respected
    c = tlm.choose_mode(np.array([5.0 - 1e-3j, 3.0 - 1e-5j]), 2)
    assert (c.index, c.n_tied) == (1, 1)


def test_least_decay_on_undamped_rigid_base_column_is_reproducible():
    """Undamped stratum on a rigid base: all propagating modes have Im k = 0 exactly, so the
    least-decay choice is the tie-break (fundamental mode), independent of rounding noise."""
    picks = []
    for eps in (0.0, 1e-13, 1e-12, 1e-11):
        rho = 2.0 * (1 + eps)
        md = tlm.column_modes(uniform_column(60, 0.5, 100.0, 200.0, rho=rho), 2 * np.pi * 20.0)
        c = tlm.choose_mode(md.kR, 2, loss=0.0)
        assert c.n_tied >= 10
        assert c.index == tlm.select_mode(md.kR, 1, loss=0.0)
        picks.append(md.kR[c.index].real)
    assert max(picks) - min(picks) <= 1e-9 * max(picks)


@pytest.mark.parametrize("beta", [0.0, 0.05, 0.2, 0.45, 0.49])
def test_material_loss_and_propagating_ratio(beta):
    """tan(delta/2) = beta/sqrt(1-beta^2) (SASSI form) and tan(atan(2 beta)/2) (form 1 + 2i beta);
    the propagating sector is the elastic one, |arg k| <= atan(0.5), widened by atan(loss)."""
    G, M = tlm.layer_moduli(2.0, 100.0, 200.0, beta, beta)
    loss = tlm.material_loss(G, M)
    assert loss == pytest.approx(beta / np.sqrt(1 - beta ** 2), abs=1e-14)
    G1, M1 = tlm.layer_moduli(2.0, 100.0, 200.0, beta, beta, form=1)
    assert tlm.material_loss(G1, M1) == pytest.approx(np.tan(0.5 * np.arctan(2 * beta)), abs=1e-14)
    assert tlm.propagating_ratio(loss) == pytest.approx(np.tan(np.arctan(0.5) + np.arcsin(beta)), rel=1e-12)
    assert tlm.propagating_ratio(0.0) == 0.5


@pytest.mark.parametrize("beta", [0.05, 0.3, 0.46, 0.49])
def test_damped_uniform_column_keeps_its_fundamental_surface_wave(beta):
    """Heavy damping (EDU-04 admits beta < 0.5): the shortest-wavelength rule must still find the
    fundamental Rayleigh and Love modes.  For a uniform column the damped eigenvalues are the elastic
    ones at the complex frequency w/sqrt(c(beta)), i.e. close to k_el (sqrt(1-beta^2) - i beta)."""
    col0 = uniform_column(40, 0.25, 100.0, 200.0)
    col = uniform_column(40, 0.25, 100.0, 200.0, beta=beta)
    om = 2 * np.pi * 20.0
    m0, md = tlm.column_modes(col0, om), tlm.column_modes(col, om)
    rot = np.sqrt(1 - beta ** 2) - 1j * beta
    for k0, kd in ((m0.kR, md.kR), (m0.kL, md.kL)):
        k_el = k0[tlm.select_mode(k0, 1, 0.0)]
        for loss in (tlm.column_loss(col), None):          # column loss, or estimated from the spectrum
            kd_sel = kd[tlm.select_mode(kd, 1, loss)]
            assert abs(kd_sel - k_el * rot) < 0.01 * abs(k_el)


def test_heavily_damped_top_layer_keeps_surface_wave():
    """Soft heavily damped top (beta 0.46) on stiff lightly damped soil: with the column loss the
    fundamental Rayleigh mode (concentrated in the top layer at 20 Hz) remains propagating; the fixed
    elastic sector would select a deep mode instead."""
    G1, M1 = tlm.layer_moduli(2.0, 100.0, 200.0, 0.46, 0.46)
    G2, M2 = tlm.layer_moduli(2.0, 300.0, 600.0, 0.02, 0.02)
    n1, n2 = 24, 20
    col = tlm.Column(h=[0.25] * n1 + [0.75] * n2, rho=[2.0] * (n1 + n2), G=[complex(G1)] * n1 + [complex(G2)] * n2,
                     M=[complex(M1)] * n1 + [complex(M2)] * n2)
    md = tlm.column_modes(col, 2 * np.pi * 20.0)
    w = 2 * np.pi * 20.0
    # the elastic sector rejects every top-layer mode and silently falls on a deep mode (c ~ 270 m/s)
    assert w / md.kR[tlm.select_mode(md.kR, 1, loss=0.0)].real > 250.0
    k = md.kR[tlm.select_mode(md.kR, 1, tlm.column_loss(col))]
    assert 100.0 < w / k.real < 120.0 and abs(k.imag) / k.real == pytest.approx(0.52, abs=0.05)
    assert k.real == md.kR[md.kR.real > 0].real.max()      # the slowest mode of the column


# ------------------------------------------------------------------ half-space buffer (UT-19)
@pytest.mark.parametrize("law", ["geometric", "uniform", "linear"])
def test_ut19_halfspace_buffer(law):
    """UT-19: Vs_hs = 1500 m/s -> H = 450 m at 5 Hz and 75 m at 30 Hz; the generated thicknesses sum
    to H, are non-decreasing and each <= lambda/8 (nl = 20)."""
    for f, H in ((5.0, 450.0), (30.0, 75.0)):
        assert tlm.halfspace_depth(1500.0, f) == pytest.approx(H, rel=1e-14)
        for h_last, vs_last in ((2.0, 300.0), (0.5, 200.0), (10.0, 1000.0), (0.0, 0.0)):
            h, used = tlm.halfspace_sublayers(1500.0, f, 20, h_last=h_last, vs_last=vs_last, law=law)
            assert h.size == 20
            assert h.sum() == pytest.approx(H, rel=1e-12)
            assert np.all(np.diff(h) >= -1e-12 * H)
            assert np.all(h <= 1500.0 / f / 8 * (1 + 1e-12))
            assert used in ("geometric", "uniform", "linear")


def test_halfspace_buffer_geometric_series_and_fallback():
    h, used = tlm.halfspace_sublayers(1000.0, 8.0, 20, h_last=1.2, vs_last=200.0)
    assert used == "geometric"
    q = h[1:] / h[:-1]
    np.testing.assert_allclose(q, q[0], rtol=1e-10)          # constant ratio
    assert h[0] == pytest.approx(1.2 * 1000 / 200)            # h1 = h_last Vs_hs / Vs_last
    h2, used2 = tlm.halfspace_sublayers(1000.0, 0.5, 20, h_last=1.2, vs_last=200.0)
    assert used2 == "uniform"                                 # cap lambda/8 violated -> uniform
    np.testing.assert_allclose(h2, 3000.0 / 20)


# ------------------------------------------------------------------ 1-D columns
def test_vertical_response_rigid_base_closed_form():
    """Uniform damped layer on a rigid base: surface/base = 1/cos(k* H) (R2 A.1)."""
    H, vs, beta = 30.0, 200.0, 0.05
    col = uniform_column(60, H / 60, vs, 2 * vs, beta=beta)
    for f in (0.5, 1.6667, 5.0):
        w = 2 * np.pi * f
        u = tlm.vertical_response(col, w, "s")
        ex = 1 / np.cos(w * H / (vs * np.sqrt(cfactor(beta))))
        assert abs(u[0] / u[-1] - ex) / abs(ex) < 1e-4


def test_vertical_response_p_uses_constrained_modulus():
    H, vs, vp, beta = 20.0, 150.0, 400.0, 0.03
    col = uniform_column(80, H / 80, vs, vp, beta=beta)
    w = 2 * np.pi * 4.0
    u = tlm.vertical_response(col, w, "p")
    ex = 1 / np.cos(w * H / (vp * np.sqrt(cfactor(beta))))
    assert abs(u[0] / u[-1] - ex) / abs(ex) < 1e-5


def test_dashpot_base_absorbs_homogeneous_halfspace():
    """Column of half-space material on its own buffer: the surface motion is the outcrop motion 2E
    (with the incident phase) -- the dashpot with V* is the exact impedance (R1 §2.4)."""
    vs, rho, beta, f = 500.0, 2.0, 0.02, 3.0
    G, M = tlm.layer_moduli(rho, vs, 2 * vs, beta, beta)
    hg, _ = tlm.halfspace_sublayers(vs, f, 20, law="uniform")
    col = tlm.build_column([5.0] * 4, [rho] * 4, [complex(G)] * 4, [complex(M)] * 4, rho, complex(G), complex(M), hg)
    w = 2 * np.pi * f
    u = tlm.vertical_response(col, w, "s")
    k = w / (vs * np.sqrt(cfactor(beta)))
    ex = 2 * np.exp(-1j * k * col.h.sum())            # 2E at the surface, E = 1 at the dashpot level
    assert abs(u[0] - ex) / abs(ex) < 1e-2          # buffer sublayers are lambda/13: small phase error
    assert abs(tlm.outcrop_response(col, w, "s") - 2 * np.exp(-1j * k * hg.sum())) < 1e-2 * 2


# ------------------------------------------------------------------ exact half-space stiffness (D-SIT-08)
def test_halfspace_stiffness_vertical_limit_is_dashpot():
    rho, vs, vp, beta, w = 2.0, 300.0, 600.0, 0.03, 2 * np.pi * 5.0
    G, M = tlm.layer_moduli(rho, vs, vp, beta, beta)
    Kout, Kin = ff.halfspace_stiffness(rho, complex(G), complex(M), w, 0.0)
    np.testing.assert_allclose(Kout, np.diag([1j * w * np.sqrt(rho * G), 1j * w * np.sqrt(rho * M)]), rtol=1e-12)
    np.testing.assert_allclose(Kin, -Kout, rtol=1e-12)


@pytest.mark.parametrize("k", [0.0, 0.1, 0.25])
def test_halfspace_stiffness_matches_deep_damped_stratum(k):
    """The exact half-space stiffness equals the surface dynamic stiffness of a deep damped
    thin-layer stratum (independent check of the signs and of the i-scaling)."""
    rho, vs, vp, b = 2.0, 100.0, 200.0, 0.10       # 10 % damping: no reflection from the deep base
    G, M = tlm.layer_moduli(rho, vs, vp, b, b)
    w = 2 * np.pi * 5.0
    lam = vs / 5.0
    n = 500                                         # 12.5 lambda_s deep, h = lambda/40
    col = uniform_column(n, lam / 40, vs, vp, rho=rho, beta=b)
    m = tlm.column_matrices(col).free()
    F = np.linalg.inv(m.K_inplane(w, k))
    Ks = np.linalg.inv(F[np.ix_([0, m.n_free], [0, m.n_free])])
    Kout, _ = ff.halfspace_stiffness(rho, complex(G), complex(M), w, k)
    assert np.abs(Ks - Kout).max() / np.abs(Kout).max() < 5e-3
    Fa = np.linalg.inv(m.K_antiplane(w, k))[0, 0]
    Ko, _ = ff.halfspace_stiffness(rho, complex(G), complex(M), w, k, inplane=False)
    assert abs(1 / Fa - Ko[0, 0]) / abs(Ko[0, 0]) < 5e-3


# ------------------------------------------------------------------ exact Green functions
def test_point_load_greens_static_boussinesq_cerruti():
    """R1 V3: quasi-static Kausel point-load functions on a deep graded stratum."""
    vs, nu, rho = 100.0, 0.3, 2.0
    vp = vs * np.sqrt(2 * (1 - nu) / (1 - 2 * nu))
    hs, z, h = [], 0.0, 0.05
    while z < 400:
        hs.append(h)
        z += h
        h = min(h * 1.06, 10.0)
    n = len(hs)
    G, M = tlm.layer_moduli(rho, vs, vp, 0.005, 0.005)
    col = tlm.Column(h=hs, rho=[rho] * n, G=[complex(G)] * n, M=[complex(M)] * n)
    md = tlm.column_modes(col, 2 * np.pi * 0.01)
    mu = rho * vs ** 2
    for r in (1.0, 2.0, 5.0):
        g = greens_tlm.point_load(md, 0, 0, r)
        B = 1 / (2 * np.pi * mu * r)
        assert g["u"].real / B == pytest.approx(1.0, rel=0.02)
        assert g["v"].real / B == pytest.approx(1 - nu, rel=0.02)
        assert g["w"].real / B == pytest.approx(-(1 - 2 * nu) / 2, rel=0.02)
        assert g["p"].real / B == pytest.approx((1 - 2 * nu) / 2, rel=0.02)
        assert g["q"].real / B == pytest.approx(1 - nu, rel=0.02)


def test_disk_load_compliance_monotone():
    vs, nu, rho = 100.0, 0.25, 2.0
    vp = vs * np.sqrt(2 * (1 - nu) / (1 - 2 * nu))
    vals = []
    for N in (4, 8, 16):
        col = uniform_column(N, 10.0 / N, vs, vp, rho=rho, beta=0.05)
        md = tlm.column_modes(col, 2 * np.pi * 2.0)
        vals.append(abs(greens_tlm.disk_load_average(md, 0, 0, 1.0)["vertical"]))
    assert vals[0] < vals[1] < vals[2]
