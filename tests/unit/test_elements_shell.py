"""Unit tests of the SHELL element (requirements 4.1, D-ELM-06/07, spec 08 4.6)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.elements import shell
from sassi.elements.base import ElementError, material_from_E_nu, material_from_M

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")

MAT = material_from_E_nu(1.0e6, 0.3, 2.0)
T = 0.1
QUAD = np.array([[0.0, 0, 0], [2.0, 0.2, 0], [2.2, 1.5, 0], [-0.1, 1.2, 0]])


def _rot(seed=3):
    Q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    return Q if np.linalg.det(Q) > 0 else -Q


def test_local_axes_spec08():
    Lam, c, xy = shell.local_frame(QUAD)
    mLI = 0.5 * (QUAD[3] + QUAD[0])
    mJK = 0.5 * (QUAD[1] + QUAD[2])
    xp = (mJK - mLI) / np.linalg.norm(mJK - mLI)
    assert np.allclose(Lam[0], xp)
    assert np.allclose(Lam[2], [0, 0, 1])                       # counter-clockwise seen from +z
    assert np.allclose(Lam[1], np.cross(Lam[2], Lam[0]))
    # triangle: L = K, x' parallel to I-J
    Lt, _, _ = shell.local_frame(QUAD[:3])
    assert np.allclose(Lt[0], (QUAD[1] - QUAD[0]) / np.linalg.norm(QUAD[1] - QUAD[0]))


@pytest.mark.parametrize("nodes", [4, 3])
def test_rigid_modes_and_zero_drilling(nodes):
    xyz = (_rot() @ QUAD[:nodes].T).T
    K, M = shell.matrices(xyz, MAT, thick=T)
    n = 6 * nodes
    ev = np.linalg.eigvalsh(K.real)
    # 6 rigid-body modes + one zero drilling DOF per node
    assert np.sum(np.abs(ev) < 1e-9 * ev.max()) == 6 + nodes
    assert ev.min() > -1e-9 * ev.max()
    Lam, _, _ = shell.local_frame(xyz)
    for a in range(nodes):
        dr = np.zeros(n)
        dr[6 * a + 3:6 * a + 6] = Lam[2]                        # rotation about the normal
        assert np.abs(K @ dr).max() < 1e-9 * np.abs(K).max()


def test_lumped_mass():
    xyz = (_rot() @ QUAD.T).T
    _, M = shell.matrices(xyz, MAT, thick=T)
    _, _, xy = shell.local_frame(xyz)
    A = shell.polygon_area(xy)
    assert np.allclose(M, np.diag(np.diag(M)))
    for a in range(4):
        assert np.allclose(np.diag(M)[6 * a:6 * a + 3], MAT.rho * T * A / 4)
        assert np.allclose(np.diag(M)[6 * a + 3:6 * a + 6], 0.0)


def test_triangle_with_repeated_node_is_padded():
    tri = QUAD[:3]
    K3, M3 = shell.matrices(tri, MAT, thick=T)
    K4, M4 = shell.matrices(np.vstack([tri, tri[2]]), MAT, thick=T)
    assert K4.shape == (24, 24) and np.allclose(K4[:18, :18], K3) and np.allclose(K4[18:], 0)
    S4 = shell.recovery(np.vstack([tri, tri[2]]), MAT, thick=T)
    assert S4.shape == (6, 24) and np.allclose(S4[:, 18:], 0)


@pytest.mark.parametrize("nodes", [4, 3])
def test_constant_membrane_and_curvature_recovery(nodes):
    xyz = QUAD[:nodes]
    E, nu = MAT.E0, MAT.nu0
    # membrane: u = 1e-3 x, v = 2e-3 y + 1e-3 x  -> eps = (1e-3, 2e-3, 1e-3)
    Lam, c, _ = shell.local_frame(xyz)
    u = np.zeros(6 * nodes)
    for a in range(nodes):
        x, y = xyz[a, 0], xyz[a, 1]
        u[6 * a:6 * a + 3] = [1e-3 * x, 2e-3 * y + 1e-3 * x, 0.0]
    eps_g = np.array([[1e-3, 0.5e-3], [0.5e-3, 2e-3]])
    A2 = Lam[:2, :2]
    eps_l = A2 @ eps_g @ A2.T
    sig_l = E / (1 - nu ** 2) * np.array([eps_l[0, 0] + nu * eps_l[1, 1], eps_l[1, 1] + nu * eps_l[0, 0],
                                         0.5 * (1 - nu) * 2 * eps_l[0, 1]])
    S = shell.recovery(xyz, MAT, thick=T)
    assert np.allclose((S @ u)[:3].real, sig_l, rtol=1e-10)
    Sf = shell.recovery(xyz, MAT, thick=T, membrane_output="force")
    assert np.allclose((Sf @ u)[:3].real, T * sig_l, rtol=1e-10)
    # bending: w = 0.5 k1 x^2 + 0.5 k2 y^2 + k12 x y ; theta_x = w,y ; theta_y = -w,x
    k1, k2, k12 = 1e-2, -2e-2, 5e-3
    q = np.zeros(6 * nodes)
    for a in range(nodes):
        x, y = xyz[a, 0], xyz[a, 1]
        q[6 * a + 2] = 0.5 * k1 * x * x + 0.5 * k2 * y * y + k12 * x * y
        q[6 * a + 3] = k2 * y + k12 * x
        q[6 * a + 4] = -(k1 * x + k12 * y)
    D = E * T ** 3 / (12 * (1 - nu ** 2))
    Mg = -D * np.array([[k1 + nu * k2, (1 - nu) * k12], [(1 - nu) * k12, k2 + nu * k1]])
    Ml = A2 @ Mg @ A2.T
    m = (S @ q).real
    assert np.allclose(m[3:], [Ml[0, 0], Ml[1, 1], Ml[0, 1]], rtol=1e-9, atol=1e-12 * D)
    assert np.allclose(m[:3], 0.0, atol=1e-12 * E)


def test_complex_stiffness_proportional():
    b = 0.05
    md = material_from_M(1, 1.0e6, 0.3, 2.0, b, b, 1.0)
    K, M = shell.matrices(QUAD, md, thick=T)
    K0, M0 = shell.matrices(QUAD, md.undamped(), thick=T)
    assert np.abs(K - K0.real * cfactor(b)).max() < 1e-12 * np.abs(K0).max()
    assert np.allclose(M, M0)


def test_warped_quad_projected_and_errors():
    warped = QUAD.copy()
    warped[2, 2] = 0.01
    K, M = shell.matrices(warped, MAT, thick=T)
    assert np.all(np.isfinite(K))
    with pytest.raises(ElementError):
        shell.matrices(QUAD, MAT, thick=0.0)
    with pytest.raises(ElementError):
        shell.matrices(np.array([[0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]]), MAT, thick=T)


def test_dkq_cylindrical_bending_matches_hermite_beam():
    # a single rectangular DKQ under cylindrical bending behaves as the exact Hermite beam strip
    a, b = 2.0, 1.0
    m = material_from_E_nu(1.0e6, 0.0, 1.0)
    xy = np.array([[0, 0], [a, 0], [a, b], [0, b]], dtype=float)
    Db = 1.0e6 * T ** 3 / 12 * shell._Dplate(0.0)
    Kb, _ = shell.dkq(xy, Db)
    # DOFs (w, theta_x, theta_y) per node; impose w = theta_x = 0 on y-edges pattern: w(x) only
    Tm = np.zeros((12, 4))
    for n, (x, y) in enumerate(xy):
        end = 0 if x == 0 else 1
        Tm[3 * n, 2 * end] = 1.0                # w_i
        Tm[3 * n + 2, 2 * end + 1] = -1.0       # theta_y = -dw/dx  (beam slope DOF)
    kb = Tm.T @ Kb @ Tm
    EI = 1.0e6 * b * T ** 3 / 12
    ref = EI / a ** 3 * np.array([[12, 6 * a, -12, 6 * a], [6 * a, 4 * a * a, -6 * a, 2 * a * a],
                                  [-12, -6 * a, 12, -6 * a], [6 * a, 2 * a * a, -6 * a, 4 * a * a]])
    assert np.allclose(kb, ref, rtol=1e-10)


def test_membrane_modes_independent_of_solid_switch():
    K1, _ = shell.matrices(QUAD, MAT, thick=T)
    K2, _ = shell.matrices(QUAD, MAT, thick=T, incompatible=False)       # MOPT <incomp> key is ignored
    K3, _ = shell.matrices(QUAD, MAT, thick=T, membrane_incompatible=False)
    assert np.allclose(K1, K2)
    assert not np.allclose(K1, K3)
    # in-plane bending of a rectangle: the incompatible membrane is softer (no locking)
    assert np.trace(K3.real) > np.trace(K1.real)


@pytest.mark.parametrize("order, keep", [((0, 1, 2, 2), [0, 1, 2]), ((0, 1, 2, 0), [0, 1, 2]),
                                         ((0, 1, 1, 2), [0, 1, 3]), ((0, 0, 1, 2), [0, 2, 3])])
def test_triangle_from_any_adjacent_repeated_corner(order, keep):
    """One pair of adjacent coincident corners (not only L = K) is a triangle: the matrices are
    the 3-node ones placed at the distinct corners, zero for the repeated one (spec 08 4.6)."""
    tri = (_rot(5) @ QUAD[:3].T).T
    X4 = tri[list(order)]
    assert shell.corner_rows(X4) == keep and shell.is_triangle(X4)
    K3, M3 = shell.matrices(tri, MAT, thick=T)
    S3 = shell.recovery(tri, MAT, thick=T)
    K4, M4 = shell.matrices(X4, MAT, thick=T)
    S4 = shell.recovery(X4, MAT, thick=T)
    idx = np.concatenate([np.arange(6 * r, 6 * r + 6) for r in keep])
    rest = np.setdiff1d(np.arange(24), idx)
    assert np.allclose(K4[np.ix_(idx, idx)], K3) and np.allclose(K4[rest], 0) and np.allclose(K4[:, rest], 0)
    assert np.allclose(M4[np.ix_(idx, idx)], M3) and np.allclose(M4[rest], 0)
    assert np.allclose(S4[:, idx], S3) and np.allclose(S4[:, rest], 0)


def test_coincident_corners_raise_error_9():
    """Reviewer case: (I, J, K, I) used to raise ZeroDivisionError; it is now a triangle.  Opposite
    or multiple coincident corners and a zero-length projected side raise ElementError (Error 9)."""
    X = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 0.0]])
    K, _ = shell.matrices(X, MAT, thick=T)
    assert K.shape == (24, 24) and np.allclose(K[18:], 0)
    bad = [np.array([[0, 0, 0], [1, 0, 0], [0, 0, 0], [0, 1, 0.0]]),        # I = K (opposite)
           np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0], [1, 0, 0.0]]),        # I = J and K = L
           np.array([[0, 0, 0], [1, 0, 0], [1, 0, 0.0]])]                    # 3-node with J = K
    for Xb in bad:
        with pytest.raises(ElementError, match="Error 9"):
            shell.matrices(Xb, MAT, thick=T)
    with pytest.raises(ElementError, match="Error 9"):
        shell._dk_constraints(np.array([[0, 0], [1, 0], [1, 1], [0, 0.0]]), [(0, 1), (1, 2), (2, 3), (3, 0)])
