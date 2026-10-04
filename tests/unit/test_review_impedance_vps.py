"""Adversarial review tests of the impedance verification problems
(sassi/verify/problems/vp_impedance.py: VP-10, VP-11, VP-13, VP-14, VP-17; requirements §6.3,
R2 §B, §C, §D.5, §F; R1 §3.4; D-ANL-07).

These tests do not repeat the implementer's own checks.  They probe the helpers and the evidence
behind the four xfail VPs from independent angles:

* the exact low-frequency radiation damping of a rigid disk, re-derived here with explicit
  real-axis branches (e^{+iwt}) and anchored to two classical results that are independent of
  the implementer's code: the Rayleigh speeds of R2 D.5 and the Miller & Pursey (1955) energy
  partition of a vertical surface source (67.4 % Rayleigh waves at nu = 1/4);
* the dispersion of a layer on a rigid base (continuum, P-SV determinant written here), which
  confirms the zero-group-velocity onset below Vp/(4H) used to explain the VP-13 vertical peak;
* an independent square-cell static BEM (own grid, own analytic self-integrals) for the VP-14
  rocking evidence, plus exact half-space results the BEM must meet (Mindlin's relaxed
  horizontal punch, nu = 1/2 welded = relaxed, Reissner-Sagoci torsion for any nu);
* the R2 F.2 3-DOF model re-derived in absolute coordinates with a coupled impedance;
* energy (Schur complement) bounds between welded and relaxed impedances;
* a small end-to-end SITE + POINT + HOUSE + ANALYS run (4 rings, ~1 s, ~1 MB) checked against
  an independent inversion of F_ff, complex symmetry, radiation-damping signs and the exact
  low-frequency damping (the VP-11 rocking-damping evidence reproduced on a coarse mesh).
"""
from __future__ import annotations

import math
import shutil

import numpy as np
import pytest
from scipy import integrate, optimize

from sassi.verify import builders as B
from sassi.verify.problems import vp_impedance as V


# ---------------------------------------------------------------------------------------
# independent exact low-frequency damping (real-axis integrals, explicit branches)
# ---------------------------------------------------------------------------------------
def _rayleigh_x(nu: float) -> float:
    """c_R / Vs from the secular equation of R2 D.5."""
    al = math.sqrt((1 - 2 * nu) / (2 * (1 - nu)))
    return optimize.brentq(lambda x: (2 - x * x) ** 2 - 4 * math.sqrt(1 - x * x) * math.sqrt(1 - al * al * x * x),
                           0.5, 0.999999)


def _lowfreq_independent(nu: float) -> dict:
    """Leading-order damping of a rigid disk (relaxed contact), k_s = 1, G = 1.

    Imaginary parts of the surface flexibility symbols on 0 < xi < 1 written out by hand for
    e^{+iwt} (vertical wavenumbers nu = i sqrt(k^2 - xi^2) for propagating waves); the Rayleigh
    pole, which the outgoing contour passes above, contributes -pi Res."""
    al = math.sqrt((1 - 2 * nu) / (2 * (1 - nu)))      # k_p / k_s
    xiR = 1.0 / _rayleigh_x(nu)

    def frayl(xi):                                      # real Rayleigh function, xi > 1
        return (2 * xi * xi - 1) ** 2 - 4 * xi * xi * math.sqrt(xi * xi - al * al) * math.sqrt(xi * xi - 1)
    d = 1e-6
    dF = (frayl(xiR + d) - frayl(xiR - d)) / (2 * d)

    def im_zz(xi):
        b = math.sqrt(1 - xi * xi)
        if xi < al:
            a = math.sqrt(al * al - xi * xi)
            return -a / ((2 * xi * xi - 1) ** 2 + 4 * xi * xi * a * b)
        c = math.sqrt(xi * xi - al * al)
        return -4 * xi * xi * c * c * b / ((2 * xi * xi - 1) ** 4 + 16 * xi ** 4 * c * c * b * b)

    def im_rr(xi):
        b = math.sqrt(1 - xi * xi)
        if xi < al:
            a = math.sqrt(al * al - xi * xi)
            return -b / ((2 * xi * xi - 1) ** 2 + 4 * xi * xi * a * b)
        c = math.sqrt(xi * xi - al * al)
        return (-1j * b / complex((2 * xi * xi - 1) ** 2, -4 * xi * xi * c * b)).imag

    pole_zz = -math.pi * (-math.sqrt(xiR ** 2 - al * al) / dF)
    pole_rr = -math.pi * (-math.sqrt(xiR ** 2 - 1) / dF)

    def body(f, p):
        return sum(integrate.quad(lambda s: s ** p * f(s), a, b, limit=400)[0] for a, b in ((0, al), (al, 1)))
    b1, b3, brr = body(im_zz, 1), body(im_zz, 3), body(im_rr, 1)
    btt = -1.0                                            # int_0^1 xi (-1/sqrt(1 - xi^2)) dxi
    p1, p3, prr = pole_zz * xiR, pole_zz * xiR ** 3, pole_rr * xiR
    return dict(c_v0=-2 / (math.pi * (1 - nu)) * (b1 + p1), A_r=-2 / (3 * math.pi * (1 - nu)) * (b3 + p3),
                c_h0=-2 / (math.pi * (2 - nu)) * (brr + prr + btt), rayleigh_fraction=p1 / (b1 + p1))


@pytest.mark.parametrize("nu, cr", [(0.25, 0.919402), (1.0 / 3.0, 0.932526)])
def test_rayleigh_speed_matches_r2_d5(nu, cr):
    assert _rayleigh_x(nu) == pytest.approx(cr, abs=2e-6)


def test_miller_pursey_partition_anchors_the_vertical_damping_integral():
    """Miller & Pursey (1955): a vertical surface source on a nu = 1/4 half-space radiates 67.4 %
    of its power as Rayleigh waves.  At low frequency the disk acts as that point source, so the
    pole share of the c_v(0) integral must be 0.674."""
    assert _lowfreq_independent(0.25)["rayleigh_fraction"] == pytest.approx(0.674, abs=0.002)


@pytest.mark.parametrize("nu", [0.25, 1.0 / 3.0, 0.4])
def test_lowfreq_damping_agrees_with_independent_integration(nu):
    mine = _lowfreq_independent(nu)
    theirs = V.lowfreq_damping(nu)
    for k in ("c_v0", "c_h0", "A_r"):
        assert theirs[k] == pytest.approx(mine[k], rel=1e-6), k


def test_veletsos_verbic_rocking_damping_is_below_the_exact_low_frequency_limit():
    """The VP-11 xfail rests on c_r(VV) ~ 0.1 a0^2 vs the exact 0.240 a0^2 (nu = 1/3)."""
    A_exact = _lowfreq_independent(1.0 / 3.0)["A_r"]
    assert A_exact == pytest.approx(0.2401, abs=5e-4)
    a0 = 0.05
    S = V._vv_impedance(a0, 1.0)
    cr_vv = S[1, 1].imag / (a0 * 8 * V.GMOD / (3 * (1 - V.NU)))
    assert cr_vv / a0 ** 2 == pytest.approx(0.1, rel=1e-3)
    assert A_exact / (cr_vv / a0 ** 2) > 2.3


@pytest.mark.parametrize("a0", [0.5, 1.0, 1.5])
def test_vv_impedance_helper_reproduces_r2_c2_table(a0):
    """_vv_impedance (used to re-derive the R2 F.1 reference) vs the R2 C.2 nu = 1/3 table."""
    R = 7.0
    S = V._vv_impedance(a0, R)
    Kx0, Kr0 = 8 * V.GMOD * R / (2 - V.NU), 8 * V.GMOD * R ** 3 / (3 * (1 - V.NU))
    kr, cr, _, _ = V.VP11_VV[a0]
    assert S[1, 1].real / Kr0 == pytest.approx(kr, abs=1e-4)
    assert S[1, 1].imag / (a0 * Kr0) == pytest.approx(cr, abs=1e-4)
    assert S[0, 0] / Kx0 == pytest.approx(1 + 0.65j * a0, rel=1e-12)
    assert S[0, 1] == 0 and S[1, 0] == 0


# ---------------------------------------------------------------------------------------
# VP-13: dispersion of the layer over a rigid base (continuum, independent of SITE)
# ---------------------------------------------------------------------------------------
def _layer_det(k: float, w: float, H: float = 3.0, vs: float = 1.0, vp: float = 2.0) -> float:
    """P-SV secular determinant of a layer, free top (z = 0), fixed bottom (z = H); potentials
    phi = A cos(pz) + B sin(pz)/p, psi = i (C cos(qz) + D sin(qz)/q): every entry is real."""
    mu, lam = vs ** 2, vp ** 2 - 2 * vs ** 2

    def basis(m2, z):
        m = np.sqrt(complex(m2))
        c = np.cos(m * z).real
        s = (np.sin(m * z) / m).real if abs(m) > 1e-12 else z
        return c, s, -m2 * s, c                       # c, s, c', s'
    p2, q2 = w * w / vp ** 2 - k * k, w * w / vs ** 2 - k * k
    cp0, sp0, dcp0, dsp0 = basis(p2, 0.0)
    cpH, spH, dcpH, dspH = basis(p2, H)
    cq0, sq0, dcq0, dsq0 = basis(q2, 0.0)
    cqH, sqH, dcqH, dsqH = basis(q2, H)
    a = -(lam * w * w / vp ** 2 + 2 * mu * p2)
    M = np.array([[a * cp0, a * sp0, 2 * mu * k * dcq0, 2 * mu * k * dsq0],
                  [-2 * mu * k * dcp0, -2 * mu * k * dsp0, mu * (q2 - k * k) * cq0, mu * (q2 - k * k) * sq0],
                  [-k * cpH, -k * spH, -dcqH, -dsqH],
                  [dcpH, dspH, k * cqH, k * sqH]])
    return float(np.linalg.det(M))


def _real_roots(w: float, kmax: float = 2.0, n: int = 2500) -> np.ndarray:
    ks = np.linspace(1e-6, kmax, n)
    d = np.array([_layer_det(k, w) for k in ks])
    return ks[np.flatnonzero(np.sign(d[1:]) != np.sign(d[:-1]))]


def test_layer_zero_group_velocity_onset_below_vertical_cutoff():
    """VP-13 evidence: the 2nd P-SV mode of the r = 0.5, H = 3, Vs = 1, Vp = 2 layer appears at a
    ZGV point with k > 0 *below* the 1-D cut-off Vp/(4H) (A0 = 0.5236)."""
    r = 0.5
    A0s = np.arange(0.495, 0.5241, 0.0025)
    onset = None
    for A0 in A0s:
        roots = _real_roots(A0 / r)
        if roots.size >= 3:
            onset = (A0, roots[0] * r)
            break
    assert onset is not None
    A0z, kRz = onset
    assert 0.50 < A0z < 0.515 < math.pi / 6
    assert 0.05 < kRz < 0.25                                  # finite wavenumber: ZGV, not a k = 0 cut-off


# ---------------------------------------------------------------------------------------
# static BEM: independent square BEM and exact half-space results
# ---------------------------------------------------------------------------------------
def _my_kernel(dx, dy, nu):
    r = np.hypot(dx, dy)
    out = np.zeros(np.shape(dx) + (3, 3))
    b, c = 1 / (2 * np.pi * r), (1 - 2 * nu) / (4 * np.pi * r * r)
    out[..., 0, 0] = b * (1 - nu + nu * dx * dx / r ** 2)
    out[..., 1, 1] = b * (1 - nu + nu * dy * dy / r ** 2)
    out[..., 0, 1] = out[..., 1, 0] = b * nu * dx * dy / r ** 2
    out[..., 2, 2] = b * (1 - nu)
    out[..., 0, 2], out[..., 1, 2] = c * dx, c * dy           # upward P_z pushes the surface outward
    out[..., 2, 0], out[..., 2, 1] = -c * dx, -c * dy         # reciprocity
    return out


def _my_square_bem(n: int, nu: float, welded: bool, ng: int = 6) -> np.ndarray:
    """Square 2 x 2 (B = 1), n x n square cells, uniform tractions, collocation at the centres;
    self cell: int dA / r = 4 s ln(1 + sqrt 2), and the cos^2 part is half of it by symmetry."""
    s = 2.0 / n
    xc = -1 + s * (np.arange(n) + 0.5)
    X, Y = np.meshgrid(xc, xc, indexing="ij")
    P = np.c_[X.ravel(), Y.ravel()]
    nc = len(P)
    g, w = np.polynomial.legendre.leggauss(ng)
    gx, gy = np.meshgrid(0.5 * s * g, 0.5 * s * g, indexing="ij")
    W = np.outer(w, w).ravel() * (0.5 * s) ** 2
    off = np.c_[gx.ravel(), gy.ravel()]
    I0 = 4 * s * math.log(1 + math.sqrt(2))
    A = np.zeros((nc, 3, nc, 3))
    for i in range(nc):
        d = P[i][None, None, :] - (P[:, None, :] + off[None, :, :])
        with np.errstate(all="ignore"):
            A[i] = np.einsum("q,jqab->ajb", W, _my_kernel(d[..., 0], d[..., 1], nu))
        A[i, :, i, :] = 0.0
        A[i, 0, i, 0] = A[i, 1, i, 1] = (1 - nu / 2) * I0 / (2 * np.pi)
        A[i, 2, i, 2] = (1 - nu) * I0 / (2 * np.pi)
    A = A.reshape(3 * nc, 3 * nc)
    xyz = np.c_[P, np.zeros(nc)]
    T = np.zeros((nc, 3, 6))
    T[:, 0, 0] = T[:, 1, 1] = T[:, 2, 2] = 1.0
    T[:, 1, 5], T[:, 0, 5] = xyz[:, 0], -xyz[:, 1]
    T[:, 2, 3], T[:, 2, 4] = xyz[:, 1], -xyz[:, 0]
    T = T.reshape(3 * nc, 6)
    with np.errstate(all="ignore"):                       # spurious FPE flags of Accelerate matmul
        if welded:
            return s * s * T.T @ np.linalg.solve(A, T)
        iz = 3 * np.arange(nc) + 2
        Tz = T[iz][:, [2, 3, 4]]
        K = np.zeros((6, 6))
        K[np.ix_([2, 3, 4], [2, 3, 4])] = s * s * Tz.T @ np.linalg.solve(A[np.ix_(iz, iz)], Tz)
        return K


def test_independent_square_bem_confirms_vp14_rocking_evidence():
    """VP-14: welded K_xx = 6.457 G B^3 and relaxed 6.238 (both above Pais & Kausel's 6.0)."""
    nu = 1.0 / 3.0
    for welded, ref in ((True, 6.4572), (False, 6.2367)):
        K = {n: _my_square_bem(n, nu, welded) for n in (16, 24)}
        ext = V.richardson([1 / 16, 1 / 24], [K[16], K[24]])
        assert ext[3, 3] == pytest.approx(ext[4, 4], rel=1e-10)
        assert ext[3, 3] == pytest.approx(ref, rel=0.01)          # 16/24 vs the implementer's 24/32
        assert ext[3, 3] > 1.035 * V.VP14_PK["Kxx"]
    Kw = {n: _my_square_bem(n, nu, True) for n in (16, 24)}
    Kimpl = {n: V.static_bem(V.square_cells(n), nu) for n in (16, 24)}
    for n in (16, 24):                                            # same discretisation -> same matrix
        np.testing.assert_allclose(np.diag(Kimpl[n]), np.diag(Kw[n]), rtol=2e-3)


@pytest.mark.parametrize("nu", [0.25, 1.0 / 3.0])
def test_static_bem_relaxed_horizontal_is_mindlin_exact(nu):
    """Relaxed horizontal punch (normal traction zero): traction q0 (1 - r^2/a^2)^-1/2 gives a
    uniform u_x = (2 - nu) Q / (8 G a), so 8 G a/(2 - nu) is exact for this contact model."""
    K = V.bem_reference("disk", (8, 16), nu, relaxed=True)
    assert K[0, 0] == pytest.approx(8 / (2 - nu), rel=3e-3)


def test_static_bem_incompressible_welded_equals_relaxed_and_torsion_independent_of_nu():
    Kw = V.bem_reference("disk", (8, 16), 0.5)
    Kr = V.bem_reference("disk", (8, 16), 0.5, relaxed=True)
    assert Kw[2, 2] == pytest.approx(8.0, rel=3e-3)                      # 4GR/(1-nu), nu = 1/2
    assert Kw[2, 2] == pytest.approx(Kr[2, 2], rel=1e-9)                 # no normal-shear coupling
    assert Kw[3, 3] == pytest.approx(16 / 3, rel=3e-3)
    assert abs(Kw[0, 4]) < 1e-9 * Kw[0, 0]
    for nu in (0.0, 0.45):
        assert V.bem_reference("disk", (8, 16), nu)[5, 5] == pytest.approx(16 / 3, rel=3e-3)


def test_surface_kernel_reciprocity_and_signs():
    rng = np.random.default_rng(3)
    d = rng.uniform(-2, 2, (20, 2))
    K = V.surface_kernel(d[:, 0], d[:, 1], 0.3)
    Km = V.surface_kernel(-d[:, 0], -d[:, 1], 0.3)
    np.testing.assert_allclose(K, np.swapaxes(Km, 1, 2), rtol=1e-13)   # Betti: G_ab(d) = G_ba(-d)
    k = V.surface_kernel(np.array([1.0]), np.array([0.0]), 0.3)[0]
    assert k[0, 2] > 0          # upward point force: surface moves outward (Boussinesq: inward when pressed)
    assert k[2, 0] < 0          # point ahead of a +x surface push moves down
    assert k[2, 2] == pytest.approx(0.7 / (2 * np.pi), rel=1e-14)


# ---------------------------------------------------------------------------------------
# R2 F.2: independent derivation in absolute coordinates with a coupled impedance
# ---------------------------------------------------------------------------------------
def _absolute_model(w, m, h, kstar, S):
    """Absolute DOFs (mass x_m, base x0, rocking th); ground displacement 1.  Massless rigid base:
    S [x0 - 1, th]^T = (shear, moment) of the structural spring k*(x_m - x0 - h th) at height h."""
    A = np.array([[kstar - w * w * m, -kstar, -kstar * h],
                  [-kstar, S[0, 0] + kstar, S[0, 1] + kstar * h],
                  [-kstar * h, S[1, 0] + kstar * h, S[1, 1] + kstar * h * h]], complex)
    b = np.array([0.0, S[0, 0], S[1, 0]], complex)
    return np.linalg.solve(A, b)[0]


def test_three_dof_matches_absolute_coordinate_model_with_coupling():
    rng = np.random.default_rng(7)
    m, h = 942.48, 10.0
    kstar = 1.4883e5 * (1 + 0.1j)
    for _ in range(5):
        Kx, Kt = rng.uniform(5e5, 2e6), rng.uniform(5e7, 2e8)
        Kxt = -rng.uniform(0, 0.1) * Kx * 10
        S = np.array([[Kx * (1 + 0.6j), Kxt * (1 - 0.2j)], [Kxt * (1 - 0.2j), Kt * (0.9 + 0.1j)]])
        w = 2 * np.pi * rng.uniform(0.5, 3.0)
        assert V._total_motion(w, m, h, kstar, S) == pytest.approx(_absolute_model(w, m, h, kstar, S), rel=1e-10)


def test_three_dof_coupled_static_springs_period():
    """Real coupled springs, no damping: the 3-DOF resonance is at (T~/T)^2 = 1 + k (C_xx + 2h C_xt +
    h^2 C_tt), C = S^-1 (a negative K_xt, as for a welded surface mat, adds flexibility)."""
    m, h, T = 942.48, 10.0, 0.5
    k = m * (2 * np.pi / T) ** 2
    S = np.array([[9.77e5, -7.3e5], [-7.3e5, 8.3e7]])
    C = np.linalg.inv(S)
    ratio = math.sqrt(1 + k * (C[0, 0] + 2 * h * C[0, 1] + h * h * C[1, 1]))
    assert ratio > math.sqrt(1 + k * (1 / S[0, 0] + h * h / S[1, 1]))     # coupling softens
    wr = 2 * np.pi / (T * ratio)
    w2 = wr * wr
    det = np.linalg.det(np.array([[k - w2 * m, -w2 * m, -w2 * m * h],
                                  [-w2 * m, S[0, 0] - w2 * m, S[0, 1] - w2 * m * h],
                                  [-w2 * m * h, S[1, 0] - w2 * m * h, S[1, 1] - w2 * m * h * h]]))
    assert abs(det) < 1e-9 * k * S[0, 0] * S[1, 1]
    near = [abs(V._total_motion(wr * f, m, h, k + 0j, S.astype(complex))) for f in (0.999, 1.001)]
    assert min(near) > 100


# ---------------------------------------------------------------------------------------
# relaxed vs welded impedance, extrapolation
# ---------------------------------------------------------------------------------------
def test_relaxed_impedance_is_softer_than_welded_and_equals_schur_complement():
    """Relaxed = displacement prescribed on one component group, tractions zero on the other:
    the energy minimum over the free components, i.e. the Schur complement of F^-1 = (F_gg)^-1,
    which can never exceed the welded energy on the vertical/rocking or horizontal/torsion subspace."""
    rng = np.random.default_rng(11)
    xy = rng.uniform(-1, 1, (9, 2))
    T = V.rigid_transform(xy)
    A = rng.normal(size=(27, 27))
    F = A @ A.T + 5 * np.eye(27)
    Kr = V.relaxed_impedance(F, T).real
    Kw = T.T @ np.linalg.solve(F, T)
    X = np.linalg.inv(F)
    iz = np.arange(2, 27, 3)
    ih = np.setdiff1d(np.arange(27), iz)
    for g, f, dofs in ((iz, ih, [2, 3, 4]), (ih, iz, [0, 1, 5])):
        schur = X[np.ix_(g, g)] - X[np.ix_(g, f)] @ np.linalg.solve(X[np.ix_(f, f)], X[np.ix_(f, g)])
        Tg = T[np.ix_(g, dofs)]
        np.testing.assert_allclose(Kr[np.ix_(dofs, dofs)], Tg.T @ schur @ Tg, rtol=1e-9, atol=1e-12)
        D = Kw[np.ix_(dofs, dofs)] - Kr[np.ix_(dofs, dofs)]
        assert np.linalg.eigvalsh(0.5 * (D + D.T)).min() > -1e-10
    assert np.allclose(Kr[np.ix_([2, 3, 4], [0, 1, 5])], 0.0)


def test_richardson_complex_arrays_and_order():
    h = [0.3, 0.2, 0.1]
    K0 = np.array([[1 + 2j, 3 - 1j], [0.5j, 4.0]])
    K1 = np.array([[0.1 - 0.2j, 1j], [2.0, -1.0]])
    for p in (1.0, 2.0):
        v = [K0 + K1 * x ** p for x in h]
        np.testing.assert_allclose(V.richardson(h, v, p), K0, rtol=1e-12, atol=1e-12)
        assert V.observed_order(h, [4.0 - 0.7 * x ** p for x in h]) == pytest.approx(p, rel=1e-6)
    # non-uniform refinement ratios (the VP meshes are 8/12/16 and 16/24/32)
    hs = [1 / 8, 1 / 12, 1 / 16]
    assert V.observed_order(hs, [2.0 + x ** 1.3 for x in hs]) == pytest.approx(1.3, rel=1e-6)
    assert V.richardson(hs, [2.0 + 5 * x for x in hs]) == pytest.approx(2.0, rel=1e-12)


# ---------------------------------------------------------------------------------------
# small end-to-end run: SITE + POINT + HOUSE + ANALYS (<impe> = 2) on a 4-ring disk
# ---------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def small_disk(tmp_path_factory):
    wd = tmp_path_factory.mktemp("review_imp_small")
    R, nd = 5.0, 4
    a0s = [0.05, 0.25, 0.5]
    site = V.soil_column(R / nd, min(2 * math.pi * R / max(a0s) / 8, 4 * R), 20 * R)
    fs, a0 = V.a0_frequencies(a0s, R)
    mdl = B.surface_rigid_mat(site, half_width=R, ndiv=nd, shape="disk")
    res = V.impedance_chain(wd, site, fs, mdl)
    from sassi.core.flexibility import flexibility_matrix
    F = [np.array(flexibility_matrix(res["FILE3"], q, res["xy"], np.ones(len(res["xy"]), int)))
         for q in range(len(a0))]
    out = dict(res, a0=a0, R=R, F=F)
    shutil.rmtree(wd, ignore_errors=True)                  # keep the nearly full disk clean
    return out


def test_small_chain_kg_independent_inverse_and_symmetry(small_disk):
    T = V.rigid_transform(small_disk["xy"])
    for q, F in enumerate(small_disk["F"]):
        with np.errstate(all="ignore"):
            Kmine = T.T @ np.linalg.inv(F) @ T             # independent of ss.impedance_matrix
        KG = small_disk["KG"][q]
        assert np.abs(KG - Kmine).max() / np.abs(KG).max() < 1e-9
        assert np.abs(KG - KG.T).max() / np.abs(KG).max() < 1e-9   # complex symmetric (reciprocity)
        assert np.all(np.diag(KG).imag > 0)                         # radiation + material damping, e^{+iwt}


def test_small_chain_welded_stiffer_than_relaxed_and_coupling_sign(small_disk):
    KG, KR = small_disk["KG"][0].real, small_disk["KR"][0].real
    for i in (0, 2, 3, 5):
        assert KG[i, i] >= KR[i, i] * (1 - 1e-12)
    R = small_disk["R"]
    assert -0.12 < KG[0, 4] / (KG[0, 0] * R) < -0.04               # welded sliding-rocking coupling < 0
    assert KG[1, 3] == pytest.approx(-KG[0, 4], rel=1e-9)            # theta_x couples to y with the other sign


def test_small_chain_radiation_damping_vs_exact_low_frequency_limits(small_disk):
    """End-to-end check of the VP-11 evidence on a coarse mesh: SASSI-EDU's relaxed damping
    (material part 2 beta/a0 removed) against the exact limits derived independently above."""
    ex = _lowfreq_independent(1.0 / 3.0)
    a0 = small_disk["a0"]
    K0 = small_disk["KR"][0].real
    beta = V.BETA_ELASTIC
    cv = small_disk["KR"][0][2, 2].imag / (a0[0] * K0[2, 2]) - 2 * beta / a0[0]
    assert cv == pytest.approx(ex["c_v0"], rel=0.05)
    q = list(np.round(a0, 2)).index(0.25)
    cr = small_disk["KR"][q][3, 3].imag / (a0[q] * K0[3, 3]) - 2 * beta / a0[q]
    assert cr == pytest.approx(V.VP11_BEM[0.25]["cr"], rel=0.10)
    assert cr == pytest.approx(ex["A_r"] * a0[q] ** 2, rel=0.10)
    assert cr > 2.0 * 0.0062                                          # R2 C.2 Veletsos-Verbic c_r(0.25)


# ---------------------------------------------------------------------------------------
# VP-17 evidence: the *dynamic* sliding-rocking coupling of the welded disk
# ---------------------------------------------------------------------------------------
def _coupling_kernel(a0: float, r, nu: float) -> np.ndarray:
    """Dynamic part of G_xz = Gc(r) dx/r (u_x from an upward unit P_z; G_zx = -Gc dx/r by
    reciprocity), units G = 1, length R:  Gc = (a0/2pi) int xi g_c(xi) J1(a0 xi r) dxi with the
    P-SV cross symbol g_c = -xi (2 xi^2 - 1 - 2 nu_p nu_s)/F (static limit (1 - 2nu)/(2 xi), i.e.
    Gc -> (1 - 2nu)/(4 pi r)).  Same outgoing contour as :func:`V.lamb_kernels`."""
    from scipy import special
    xi, dxi = V._wavenumber_contour()
    eta = math.sqrt(2 * (1 - nu) / (1 - 2 * nu))
    vp, vs = np.sqrt(xi * xi - 1 / eta ** 2 + 0j), np.sqrt(xi * xi - 1 + 0j)
    F = (2 * xi * xi - 1) ** 2 - 4 * xi * xi * vp * vs
    fc = (xi * (-xi * (2 * xi * xi - 1 - 2 * vp * vs) / F) - (1 - 2 * nu) / 2) * dxi
    r = np.atleast_1d(np.asarray(r, float))
    out = np.zeros(r.size, complex)
    with np.errstate(all="ignore"):
        for i0 in range(0, r.size, 32):
            out[i0:i0 + 32] = a0 / (2 * np.pi) * (special.jv(1, a0 * xi[None] * r[i0:i0 + 32, None]) @ fc)
    return out


def test_dynamic_coupling_kernel_reproduces_wong_ur0():
    """R2 D.3 column U_r0 (radial displacement, vertical force pressing down; mu r u_r / P), nu = 0.33."""
    nu = 0.33
    r = np.array([0.5, 1.0, 2.0, 3.0, 5.0])
    U = -r * ((1 - 2 * nu) / (4 * np.pi * r) + _coupling_kernel(1.0, r, nu))
    wong = np.array([-.032 + .007j, -.033 + .025j, .006 + .060j, .074 + .035j, .005 - .120j])
    assert np.abs(U - wong).max() < 0.002


def _welded_dynamic_bem(cells: np.ndarray, a0: float, nu: float = 1.0 / 3.0, nq: int = 4) -> np.ndarray:
    """Welded-contact dynamic impedance (G = 1, length R): the implementer's static singular part and
    zz/xx/yy/xy dynamic parts plus the independently derived xz/zx coupling of _coupling_kernel."""
    from scipy.interpolate import CubicSpline
    nc = cells.shape[0]
    with np.errstate(all="ignore"):
        A0, Pc, area = V.bem_influence(cells, nu)
        Q, W = V._cell_quadrature(cells, nq)
        rg = np.linspace(0.0, 2.2 * float(np.abs(cells).max()), 300)
        spl = [CubicSpline(rg, v) for v in V.lamb_kernels(a0, rg, nu)] + [CubicSpline(rg, _coupling_kernel(a0, rg, nu))]
        A = A0.reshape(nc, 3, nc, 3).astype(complex)
        for i0 in range(0, nc, 32):
            d = Pc[i0:i0 + 32, None, None, :] - Q[None]
            r = np.hypot(d[..., 0], d[..., 1])
            ux, uy = d[..., 0] / r, d[..., 1] / r
            psi = np.arctan2(d[..., 1], d[..., 0])
            z, a, b, c = (s(r) for s in spl)
            c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
            A[i0:i0 + 32, 2, :, 2] += np.einsum("jq,ijq->ij", W, z)
            A[i0:i0 + 32, 0, :, 0] += np.einsum("jq,ijq->ij", W, a + b * c2)
            A[i0:i0 + 32, 1, :, 1] += np.einsum("jq,ijq->ij", W, a - b * c2)
            bxy = np.einsum("jq,ijq->ij", W, b * s2)
            A[i0:i0 + 32, 0, :, 1] += bxy
            A[i0:i0 + 32, 1, :, 0] += bxy
            cx, cy = np.einsum("jq,ijq->ij", W, c * ux), np.einsum("jq,ijq->ij", W, c * uy)
            A[i0:i0 + 32, 0, :, 2] += cx
            A[i0:i0 + 32, 1, :, 2] += cy
            A[i0:i0 + 32, 2, :, 0] -= cx
            A[i0:i0 + 32, 2, :, 1] -= cy
        return V._bem_rigid(A.reshape(3 * nc, 3 * nc), Pc, area, relaxed=False)


def test_welded_dynamic_impedance_matches_independent_welded_bem(tmp_path):
    """SASSI-EDU's welded K_G (ANALYS <impe> = 2) at the VP-17 resonance a0 = 1.07, rings 6/8
    extrapolated, against the welded exact-kernel BEM (6/8 polar rings extrapolated).  This checks
    the dynamic sliding-rocking coupling (Re < 0, Im > 0) that lowers the VP-17 peak, and the
    welded rocking damping of VP-11."""
    R, a0t = 5.0, 1.07
    unit = np.array([[R ** (1 + (i > 2) + (j > 2)) for j in range(6)] for i in range(6)]) * V.GMOD
    K = {}
    for nd in (6, 8):
        wd = tmp_path / f"w{nd}"
        wd.mkdir()
        site = V.soil_column(R / nd, min(2 * math.pi * R / a0t / 8, 4 * R), 20 * R)
        fs, a0 = V.a0_frequencies([a0t], R)
        mdl = B.surface_rigid_mat(site, half_width=R, ndiv=nd, shape="disk")
        B.run_soil(wd, "m", site, fs, layer=0, rad=mdl.rad)
        B.run_house(wd, "m", mdl)
        K[nd] = V.analys_global_impedance(wd, "m", fs)[0] / unit
        shutil.rmtree(wd, ignore_errors=True)
    Ks = V.richardson([1 / 6, 1 / 8], [K[6], K[8]])
    Kb = V.richardson([1 / 6, 1 / 8], [_welded_dynamic_bem(V.disk_cells(n), a0t) for n in (6, 8)])
    assert Kb[0, 4].real < 0 < Kb[0, 4].imag                               # BEM: coupling sign pattern
    assert Ks[0, 4].real < 0 < Ks[0, 4].imag
    assert abs(Ks[0, 4] - Kb[0, 4]) < 0.05 * abs(Kb[0, 4])
    for i in (0, 2, 3):
        assert abs(Ks[i, i] - Kb[i, i]) < 0.04 * abs(Kb[i, i]), i
    assert Ks[3, 3].imag == pytest.approx(Kb[3, 3].imag, rel=0.06)        # welded rocking damping


# ---------------------------------------------------------------------------------------
# robustness of helpers (defects found in review)
# ---------------------------------------------------------------------------------------
def test_graded_thicknesses_caps_every_sublayer():
    """Review defect (minor), fixed: graded_thicknesses did not cap the first sublayer at tmax, so
    soil_column(t0 > lambda/8, ...) silently broke the h <= lambda/8 rule it documents."""
    assert max(V.graded_thicknesses(2.0, 1.0, 10.0)) <= 1.0


def test_relaxed_bem_kernel_table_covers_all_pair_distances_of_a_square():
    """Review defect (minor), fixed: relaxed_bem_impedance tabulated the dynamic kernels only up to
    r = 2.2 max|V|; for a square (max distance 2.83 B) the cubic spline extrapolated.  The table now
    runs to V.max_pair_distance of the evaluated points (the original assertion re-implemented the
    old 2.2 max|V| bound literally, so it could not see the fix)."""
    cells = V.square_cells(4)
    Pw, _ = V._cell_quadrature(cells, 4)
    pts = Pw.reshape(-1, 2)
    dmax = float(np.max(np.hypot(*(pts[:, None, :] - pts[None, :, :]).transpose(2, 0, 1))))
    assert V.max_pair_distance(pts, pts) >= dmax
