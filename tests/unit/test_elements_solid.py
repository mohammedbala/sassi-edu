"""Unit tests of the SOLID and PLANE continuum elements (requirements 4.1, D-ELM-01..03, D-CNV-04/05)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.elements import plane, solid
from sassi.elements.base import ElementError, material_from_E_nu, material_from_M

# numpy 2.0 + macOS Accelerate raises spurious FP flags in matmul (results are finite and exact)
pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")

MAT = material_from_E_nu(1000.0, 0.25, 2.0)


def _distorted_hex(seed=1):
    rng = np.random.default_rng(seed)
    return (solid.NODE_NAT + 1.0) / 2.0 * np.array([2.0, 1.0, 1.5]) + rng.uniform(-0.08, 0.08, (8, 3))


def _strain_field():
    eps = np.array([1e-3, 2e-3, -1e-3, 5e-4, 3e-4, -2e-4])
    E = np.array([[eps[0], eps[3] / 2, eps[4] / 2], [eps[3] / 2, eps[1], eps[5] / 2], [eps[4] / 2, eps[5] / 2, eps[2]]])
    return eps, E


@pytest.mark.parametrize("inc", [False, True])
@pytest.mark.parametrize("eint", [0, 1, 2])
def test_solid_symmetry_rigid_modes_and_mass(inc, eint):
    xyz = _distorted_hex()
    K, M = solid.matrices(xyz, MAT, incompatible=inc, eint=eint)
    assert K.shape == (24, 24) and M.shape == (24, 24)
    assert np.abs(K - K.T).max() < 1e-12 * np.abs(K).max()
    ev = np.linalg.eigvalsh(K.real)
    assert np.sum(np.abs(ev) < 1e-9 * ev.max()) == 6               # 3 translations + 3 rotations
    assert ev.min() > -1e-9 * ev.max()
    vol = M.sum() / 3.0 / MAT.rho
    # volume of the trilinear element by direct Gauss integration
    pts, w = solid.gauss_3d(3)
    _, dN = solid.hex8_shape(pts)
    assert vol == pytest.approx(np.sum(w * np.linalg.det(np.einsum("gai,aj->gij", dN, xyz))), rel=1e-12)


def test_solid_mixed_mass_is_half_lumped_half_consistent():
    xyz = _distorted_hex()
    _, Mc = solid.matrices(xyz, MAT, mass="consistent")
    _, Ml = solid.matrices(xyz, MAT, mass="lumped")
    _, Mm = solid.matrices(xyz, MAT)
    assert np.allclose(Ml, np.diag(Mc.sum(axis=1)))
    assert np.allclose(Mm, 0.5 * Ml + 0.5 * Mc)


@pytest.mark.parametrize("inc", [False, True])
def test_solid_constant_strain_recovery_exact(inc):
    xyz = _distorted_hex()
    eps, E = _strain_field()
    u = (xyz @ E.T).ravel()
    S = solid.recovery(xyz, MAT, incompatible=inc)
    assert S.shape == (7, 24)
    assert np.allclose(S[6], 0.0)                                  # SOCT row (time domain in STRESS)
    assert np.allclose(S[:6] @ u, solid.isotropic_D(MAT.lam, MAT.G) @ eps, rtol=1e-12, atol=1e-14)
    B = solid.recovery(xyz, MAT, incompatible=inc, quantity="strain")
    assert np.allclose(B @ u, eps, atol=1e-15)


def test_solid_rigid_rotation_gives_no_stress():
    xyz = _distorted_hex()
    W = np.array([[0, -0.01, 0.02], [0.01, 0, -0.03], [-0.02, 0.03, 0]])
    u = (xyz @ W.T).ravel()
    assert np.abs(solid.recovery(xyz, MAT, incompatible=True) @ u).max() < 1e-12


def test_solid_incompatible_modes_exact_pure_bending():
    # rectangular brick under pure bending: u_x = k x z, u_z = -k x^2/2 (nu = 0) is exact with the modes
    m = material_from_E_nu(1000.0, 0.0, 1.0)
    a, b, c = 2.0, 1.0, 1.0
    xyz = (solid.NODE_NAT + 1.0) / 2.0 * np.array([a, b, c]) - np.array([0, 0, c / 2])
    k = 1e-3
    u = np.column_stack([k * xyz[:, 0] * xyz[:, 2], 0 * xyz[:, 0], -0.5 * k * xyz[:, 0] ** 2]).ravel()
    for inc, exact in ((True, True), (False, False)):
        K, _ = solid.matrices(xyz, m, incompatible=inc)
        energy = 0.5 * u @ K.real @ u
        ref = 0.5 * 1000.0 * k * k * b * c ** 3 / 12.0 * a        # 1/2 E I k^2 L
        if exact:
            assert energy == pytest.approx(ref, rel=1e-10)
        else:
            assert energy > 1.5 * ref                              # locking of the compatible brick


def test_solid_degenerate_prism_and_pyramid():
    base = (solid.NODE_NAT + 1.0) / 2.0
    seven = base.copy()
    seven[7] = seven[6]                                            # 7 = 8 (7 distinct nodes)
    wedge = base.copy()
    wedge[4] = wedge[5] = (0.5, 0.0, 1.0)                          # 5 = 6 and 7 = 8: ridge at x = 0.5
    wedge[6] = wedge[7] = (0.5, 1.0, 1.0)
    pyr = base.copy()
    pyr[4:] = np.array([0.5, 0.5, 1.0])                            # 5 = 6 = 7 = 8
    pts, w = solid.gauss_3d(4)
    _, dN = solid.hex8_shape(pts)
    vol7 = np.sum(w * np.linalg.det(np.einsum("gai,aj->gij", dN, seven)))
    for xyz, vol in ((seven, vol7), (wedge, 0.5), (pyr, 1.0 / 3.0)):
        for inc in (False, True):
            K, M = solid.matrices(xyz, MAT, incompatible=inc)
            assert M.sum() / 3 / MAT.rho == pytest.approx(vol, rel=1e-12)
            assert np.all(np.isfinite(K))
            # collapse the repeated nodes: the reduced element still has exactly 6 rigid modes
            nodes = np.unique(xyz, axis=0, return_inverse=True)[1].ravel()
            T = np.zeros((24, 3 * (nodes.max() + 1)))
            for a, n in enumerate(nodes):
                T[3 * a:3 * a + 3, 3 * n:3 * n + 3] = np.eye(3)
            ev = np.linalg.eigvalsh(T.T @ K.real @ T)
            assert np.sum(np.abs(ev) < 1e-9 * ev.max()) == 6


def test_solid_negative_jacobian_edu05():
    xyz = (solid.NODE_NAT + 1.0) / 2.0
    bad = xyz[[4, 5, 6, 7, 0, 1, 2, 3]]                            # faces swapped: inverted element
    with pytest.raises(ElementError, match="EDU-05"):
        solid.matrices(bad, MAT)
    with pytest.raises(ElementError):
        solid.matrices(xyz, MAT, eint=5)


@pytest.mark.parametrize("inc", [False, True])
def test_solid_complex_stiffness_equals_K0_c_when_dampings_equal(inc):
    b = 0.07
    md = material_from_M(1, 1000.0, 0.25, 2.0, b, b, 1.0)
    K, _ = solid.matrices(_distorted_hex(), md, incompatible=inc)
    K0, _ = solid.matrices(_distorted_hex(), md.undamped(), incompatible=inc)
    assert np.abs(K - K0.real * cfactor(b)).max() < 1e-12 * np.abs(K0).max()
    # different beta_p and beta_s: K* = lam* K_lam + G* K_G is not proportional to K0
    md2 = material_from_M(1, 1000.0, 0.25, 2.0, 0.02, 0.08, 1.0)
    K2, _ = solid.matrices(_distorted_hex(), md2, incompatible=inc)
    assert np.abs(K2 - K0.real * cfactor(0.08)).max() > 1e-6 * np.abs(K0).max()


def test_solid_batch_equals_single():
    xyz = np.stack([_distorted_hex(s) for s in range(3)])
    mats = [MAT, material_from_E_nu(500.0, 0.3, 1.0), MAT]
    out = solid.batch(xyz, mats, incompatible=True, eint=1)
    for e in range(3):
        K, M = solid.matrices(xyz[e], mats[e], incompatible=True, eint=1)
        assert np.allclose(out["K"][e], K) and np.allclose(out["M"][e], M)
        assert np.allclose(out["S"][e], solid.recovery(xyz[e], mats[e], incompatible=True, eint=1))


# ---------------------------------------------------------------- PLANE
def _quad():
    return np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.1], [1.8, 0.0, 1.2], [0.1, 0.0, 0.9]])


@pytest.mark.parametrize("inc", [False, True])
def test_plane_rigid_modes_and_constant_strain(inc):
    xyz = _quad()
    K, M = plane.matrices(xyz, MAT, incompatible=inc)
    ev = np.linalg.eigvalsh(K.real)
    assert np.sum(np.abs(ev) < 1e-9 * ev.max()) == 3
    x, z = xyz[:, 0], xyz[:, 2]
    u = np.column_stack([1e-3 * x + 2e-4 * z, 3e-4 * x - 2e-3 * z]).ravel()
    eps = np.array([1e-3, -2e-3, 5e-4])
    S = plane.recovery(xyz, MAT, incompatible=inc)
    assert np.allclose(S @ u, plane.plane_strain_D(MAT.lam, MAT.G) @ eps, rtol=1e-12)
    area = 0.5 * abs(np.dot(x, np.roll(z, -1)) - np.dot(np.roll(x, -1), z))
    assert M.sum() / 2 / MAT.rho == pytest.approx(area, rel=1e-12)


def test_plane_orientation_is_normalised():
    xyz = _quad()
    rev = xyz[[0, 3, 2, 1]]                                        # clockwise numbering
    assert plane.orientation(xyz) == 1 and plane.orientation(rev) == -1
    K, M = plane.matrices(xyz, MAT, incompatible=True)
    Kr, Mr = plane.matrices(rev, MAT, incompatible=True)
    p = np.array([0, 3, 2, 1])
    d = np.column_stack([2 * p, 2 * p + 1]).ravel()
    assert np.allclose(Kr, K[np.ix_(d, d)])
    assert np.allclose(Mr, M[np.ix_(d, d)])
    assert np.allclose(plane.recovery(rev, MAT)[:, :], plane.recovery(xyz, MAT)[:, d])


def test_plane_triangle_by_repeated_node_and_three_nodes():
    tri = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    K3, M3 = plane.matrices(tri, MAT)
    K4, M4 = plane.matrices(np.vstack([tri, tri[2]]), MAT)
    assert K3.shape == (6, 6) and K4.shape == (8, 8)
    # folding the repeated node of the 4-node form gives the 3-node form
    F = np.zeros((8, 6))
    F[:6, :6] = np.eye(6)
    F[6:, 4:] = np.eye(2)
    assert np.allclose(F.T @ K4 @ F, K3) and np.allclose(F.T @ M4 @ F, M3)
    assert M3.sum() / 2 / MAT.rho == pytest.approx(0.5)


def test_plane_complex_stiffness_proportional():
    b = 0.05
    md = material_from_M(1, 1000.0, 0.3, 1.0, b, b, 1.0)
    K, _ = plane.matrices(_quad(), md, incompatible=True)
    K0, _ = plane.matrices(_quad(), md.undamped(), incompatible=True)
    assert np.abs(K - K0.real * cfactor(b)).max() < 1e-12 * np.abs(K0).max()
