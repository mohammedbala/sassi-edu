"""Adversarial review tests of the HOUSE element library (work package 'elements').

These tests probe the element formulations from angles independent of the implementer's own
tests: equivalences between element types (PLANE = tied extruded SOLID, GENERAL with a beam's
local matrix = BEAMS), exact closed-form fields (pure bending with anticlastic curvature,
Timoshenko propped cantilever, DKQ cylindrical bending of a strip, Navier plate deflection),
objectivity under rigid rotations, rigid-body kinetic energies, damping split between bulk
and shear (beta_p != beta_s), and assembly against a hand scatter of the element matrices.

Requirements exercised: 1.5, 4.0.2 (complex moduli per element type), 4.1 (formulations and
DOF management), 4.4 item 1-2 (structure / excavated matrices), D-ELM-01..09, D-CNV-04/05,
D-STR-05, spec 08 sections 4.4-4.9, R2 I.1.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from sassi.conventions import cfactor, complex_lame
from sassi.elements import (ElemRecord, assemble, beam, build_dofmap, general, material_from_E_nu,
                            material_from_M, plane, shell, solid, spring)
from sassi.elements.base import ElementError, section_rectangle

# numpy 2.0 + macOS Accelerate raises spurious FP flags in matmul for finite products
pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")


def _rotation(seed: int) -> np.ndarray:
    Q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    return Q if np.linalg.det(Q) > 0 else -Q


def _distorted_hex(seed: int = 3, amp: float = 0.15) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return (solid.NODE_NAT + 1.0) / 2.0 * np.array([2.0, 1.0, 1.5]) + rng.uniform(-amp, amp, (8, 3))


def _collapse(X: np.ndarray, pattern: str) -> np.ndarray:
    """Degenerate SOLID shapes of spec 08 section 4.4 by repeated node coordinates."""
    X = X.copy()
    if pattern == "7=8":
        X[7] = X[6]
    elif pattern == "wedge":            # 5 = 6 and 7 = 8
        X[5] = X[4]
        X[7] = X[6]
    elif pattern == "pyramid":          # 5 = 6 = 7 = 8
        X[5] = X[4]
        X[6] = X[4]
        X[7] = X[4]
    return X


def _hex_volume(X: np.ndarray) -> float:
    pts, wts = solid.gauss_3d(3)
    _, dN = solid.hex8_shape(pts)
    J = np.einsum("gai,aj->gij", dN, X)
    return float((np.linalg.det(J) * wts).sum())


# ======================================================================================
# SOLID
# ======================================================================================
@pytest.mark.parametrize("pattern", ["hex", "7=8", "wedge", "pyramid"])
@pytest.mark.parametrize("eint", [0, 1, 2])
def test_solid_incompatible_modes_inert_for_linear_fields(pattern, eint):
    """Taylor-corrected modes satisfy integral(B_alpha dV) = 0, so for ANY shape (incl. the
    degenerate prisms/pyramid) a linear displacement field gives K_inc u = K_comp u, and the
    recovered stress is D eps exactly (D-ELM-02, patch test at the element level)."""
    mat = material_from_E_nu(1000.0, 0.3, 2.0)
    X = _distorted_hex() if pattern == "hex" else _collapse(_distorted_hex(), pattern)
    Eps = np.array([[1e-3, 2e-4, -3e-4], [2e-4, -5e-4, 4e-4], [-3e-4, 4e-4, 7e-4]])
    u = (X @ Eps.T).ravel()
    Kc, _ = solid.matrices(X, mat, incompatible=False, eint=eint)
    Ki, _ = solid.matrices(X, mat, incompatible=True, eint=eint)
    f = Kc @ u
    assert np.abs(Ki @ u - f).max() <= 1e-12 * np.abs(f).max()
    eps6 = np.array([Eps[0, 0], Eps[1, 1], Eps[2, 2], 2 * Eps[0, 1], 2 * Eps[0, 2], 2 * Eps[1, 2]])
    s_exact = solid.isotropic_D(mat.lam, mat.G) @ eps6
    S = solid.recovery(X, mat, incompatible=True, eint=eint)
    assert S.shape == (7, 24)
    assert np.allclose(S[:6] @ u, s_exact, rtol=0, atol=1e-12 * np.abs(s_exact).max())
    assert np.all(S[6] == 0)                                      # SOCT row (ARCHITECTURE 6.2)
    # six rigid-body modes, all other eigenvalues positive
    w = np.linalg.eigvalsh(Ki.real)
    assert int((np.abs(w) < 1e-9 * w.max()).sum()) == 6 and w[6] > 1e-6 * w.max()


@pytest.mark.parametrize("inc", [False, True])
def test_solid_eint_identical_for_parallelepiped(inc):
    """For a rectangular parallelepiped the integrands are at most quadratic per direction, so
    the 2x2x2 rule (EINT 0) is already exact: EINT 1/2 must give the same K (D-ELM-01); for a
    distorted brick they must differ (the rule really changes)."""
    mat = material_from_E_nu(1000.0, 0.25, 2.0)
    box = (solid.NODE_NAT + 1.0) / 2.0 * np.array([2.0, 0.7, 1.3]) + np.array([5.0, -1.0, 2.0])
    K0 = solid.matrices(box, mat, incompatible=inc, eint=0)[0]
    for eint in (1, 2):
        assert np.allclose(solid.matrices(box, mat, incompatible=inc, eint=eint)[0], K0, rtol=0,
                           atol=1e-12 * np.abs(K0).max())
    Xd = _distorted_hex(amp=0.2)
    Kd0 = solid.matrices(Xd, mat, incompatible=inc, eint=0)[0]
    Kd2 = solid.matrices(Xd, mat, incompatible=inc, eint=2)[0]
    assert np.abs(Kd2 - Kd0).max() > 1e-6 * np.abs(Kd0).max()


def test_solid_pure_bending_with_anticlastic_curvature_energy_exact():
    """Pure bending of a block with nu != 0: u = -k x z, v = nu k y z, w = k/2 (x^2 + nu (z^2 - y^2))
    has eps_xx = -k z, eps_yy = eps_zz = nu k z, no shear.  The 9 Wilson modes contain all three
    quadratic terms of w, so the condensed element energy equals the exact 1/2 E k^2 int z^2 dV;
    the compatible brick locks (energy much higher)."""
    E, nu = 1000.0, 0.3
    mat = material_from_E_nu(E, nu, 1.0)
    a, b, h = 2.0, 1.0, 0.8
    box = (solid.NODE_NAT + 1.0) / 2.0 * np.array([a, b, h])
    c = np.array([a, b, h]) / 2.0
    k = 1e-3

    def field(p):
        x, y, z = p - c
        return np.array([-k * x * z, nu * k * y * z, 0.5 * k * (x * x + nu * (z * z - y * y))])

    u = np.concatenate([field(p) for p in box])
    U_exact = 0.5 * E * k * k * (a * b * h ** 3 / 12.0)
    Ki = solid.matrices(box, mat, incompatible=True)[0].real
    Kc = solid.matrices(box, mat, incompatible=False)[0].real
    assert 0.5 * u @ Ki @ u == pytest.approx(U_exact, rel=1e-12)
    assert 0.5 * u @ Kc @ u > 1.5 * U_exact


def test_solid_objectivity_under_rigid_rotation():
    """K(R X + c) = T K(X) T^T with T = blockdiag(R): the element (incl. incompatible modes and
    the centroid-Jacobian correction) has no preferred global direction; M is invariant too."""
    mat = material_from_M(1, 1000.0, 0.3, 2.0, 0.03, 0.07, 1.0)
    X = _distorted_hex(5)
    R = _rotation(1)
    Xr = X @ R.T + np.array([10.0, -3.0, 4.0])
    T = np.kron(np.eye(8), R)
    for inc in (False, True):
        K, M = solid.matrices(X, mat, incompatible=inc, eint=1)
        Kr, Mr = solid.matrices(Xr, mat, incompatible=inc, eint=1)
        assert np.abs(Kr - T @ K @ T.T).max() <= 1e-12 * np.abs(K).max()
        assert np.abs(Mr - T @ M @ T.T).max() <= 1e-12 * np.abs(M).max()
        # stresses rotate as tensors: sigma' = R sigma R^T for u' = R u
        S, Sr = solid.recovery(X, mat, incompatible=inc), solid.recovery(Xr, mat, incompatible=inc)
        u = np.random.default_rng(2).normal(size=24)
        s, sr = S @ u, Sr @ (T @ u)

        def tens(v):
            return np.array([[v[0], v[3], v[4]], [v[3], v[1], v[5]], [v[4], v[5], v[2]]])

        assert np.allclose(tens(sr), R @ tens(s) @ R.T, atol=1e-12 * np.abs(s).max())


@pytest.mark.parametrize("pattern", ["hex", "7=8", "wedge", "pyramid"])
@pytest.mark.parametrize("scheme", ["mixed", "consistent", "lumped"])
def test_solid_mass_total_equals_rho_volume(pattern, scheme):
    """A rigid translation carries the full mass rho V in every scheme and for every degenerate
    shape (D-CNV-05); the volumes are computed independently (pyramid = V of the base x h / 3)."""
    rho = 2.5
    mat = material_from_E_nu(1000.0, 0.3, rho)
    if pattern == "pyramid":
        X = _collapse((solid.NODE_NAT + 1.0) / 2.0, "pyramid")    # unit square base, apex (0,0,1)
        V = 1.0 / 3.0
    elif pattern == "wedge":
        X = _collapse((solid.NODE_NAT + 1.0) / 2.0, "wedge")      # triangular prism-like wedge
        V = _hex_volume(X)
        assert V == pytest.approx(0.5)
    else:
        X = _distorted_hex() if pattern == "hex" else _collapse(_distorted_hex(), pattern)
        V = _hex_volume(X)
    M = solid.matrices(X, mat, mass=scheme)[1]
    for d in range(3):
        e = np.zeros(24)
        e[d::3] = 1.0
        assert e @ M @ e == pytest.approx(rho * V, rel=1e-12)


def test_solid_bulk_and_shear_damping_are_separate():
    """beta_p acts on M = lam + 2G and beta_s on G (D-CNV-04): a pure shear strain gives
    tau = G0 c(beta_s) gamma, a uniaxial strain e_xx gives sigma_xx = M0 c(beta_p) e_xx."""
    G0, M0, bs, bp = 400.0, 1500.0, 0.08, 0.02
    mat = material_from_M(2, M0, G0, 1.0, bp, bs, 1.0)
    X = _distorted_hex(7)
    S = solid.recovery(X, mat, incompatible=True)
    gam = 1e-3
    u_shear = np.concatenate([[gam * p[1], 0.0, 0.0] for p in X])   # gamma_xy = gam
    u_uni = np.concatenate([[gam * p[0], 0.0, 0.0] for p in X])     # eps_xx = gam
    s1, s2 = S @ u_shear, S @ u_uni
    assert s1[3] == pytest.approx(G0 * cfactor(bs) * gam, rel=1e-12)
    assert s2[0] == pytest.approx(M0 * cfactor(bp) * gam, rel=1e-12)
    Gs, Ms_, lam = complex_lame(G0, M0, bs, bp)
    assert s2[1] == pytest.approx(lam * gam, rel=1e-12)


# ======================================================================================
# PLANE
# ======================================================================================
def _extruded_brick(quad_xz: np.ndarray, t: float):
    """SOLID brick obtained by extruding a PLANE quad (x, z) by t along y, and the tie matrix
    T (24 x 8) that makes front/back nodes move together with UY = 0 (plane strain)."""
    qmap = {(-1, -1): 0, (1, -1): 1, (1, 1): 2, (-1, 1): 3}
    X = np.zeros((8, 3))
    T = np.zeros((24, 8))
    for a, (xi, eta, zeta) in enumerate(solid.NODE_NAT):
        q = qmap[(int(xi), int(zeta))]
        X[a] = (quad_xz[q, 0], 0.5 * (eta + 1.0) * t, quad_xz[q, 1])
        T[3 * a, 2 * q] = 1.0
        T[3 * a + 2, 2 * q + 1] = 1.0
    return X, T


@pytest.mark.parametrize("scheme", ["mixed", "consistent", "lumped"])
def test_plane_equals_tied_extruded_solid(scheme):
    """Plane strain per unit thickness (spec 08 4.7): a SOLID brick extruded by t along y with
    UY = 0 and front/back nodes tied reproduces t x the PLANE stiffness (complex, beta_p != beta_s),
    t x its mass in every scheme, and the PLANE centroid stresses SXX SZZ TXZ."""
    mat = material_from_M(1, 1000.0, 0.3, 2.0, 0.04, 0.08, 1.0)
    quad = np.array([[0.1, 0.0], [1.3, 0.2], [1.2, 1.1], [-0.1, 0.9]])
    t = 0.7
    X, T = _extruded_brick(quad, t)
    Ks, Ms = solid.matrices(X, mat, incompatible=False, eint=0, mass=scheme)
    xyz = np.column_stack([quad[:, 0], np.zeros(4), quad[:, 1]])
    Kp, Mp = plane.matrices(xyz, mat, incompatible=False, mass=scheme)
    assert np.abs(T.T @ Ks @ T - t * Kp).max() <= 1e-12 * np.abs(Kp).max()
    assert np.abs(T.T @ Ms @ T - t * Mp).max() <= 1e-12 * np.abs(Mp).max()
    Ss = solid.recovery(X, mat) @ T
    assert np.abs(Ss[[0, 2, 4]] - plane.recovery(xyz, mat)).max() <= 1e-12 * np.abs(Ss).max()


def test_plane_strain_pure_bending_energy_exact_with_incompatible_modes():
    """Plane-strain pure bending of a rectangle: u = -k x z, w = k x^2/2 + c z^2/2 with
    c = lam/(lam+2G) k (sigma_zz = 0).  QM6-type modes make the element energy exact,
    1/2 Ebar k^2 a h^3/12 with Ebar = M - lam^2/M (D-ELM-03)."""
    mat = material_from_E_nu(1000.0, 0.3, 1.0)
    lam, G = mat.lam0, mat.G0
    Mc = lam + 2 * G
    a, h, k = 2.0, 1.0, 1e-3
    rect = np.array([[0, 0, 0], [a, 0, 0], [a, 0, h], [0, 0, h]], dtype=float)
    c = lam / Mc * k

    def field(p):
        x, z = p[0] - a / 2, p[2] - h / 2
        return np.array([-k * x * z, 0.5 * k * x * x + 0.5 * c * z * z])

    u = np.concatenate([field(p) for p in rect])
    U = 0.5 * (Mc - lam * lam / Mc) * k * k * a * h ** 3 / 12.0
    Ki = plane.matrices(rect, mat, incompatible=True)[0].real
    assert 0.5 * u @ Ki @ u == pytest.approx(U, rel=1e-12)


@pytest.mark.parametrize("inc", [False, True])
def test_plane_clockwise_input_gives_same_physics(inc):
    """Either orientation is accepted (spec 08 Q12): a clockwise node list gives the matrices of
    the counter-clockwise element permuted to the caller's order."""
    mat = material_from_E_nu(1000.0, 0.25, 2.0)
    ccw = np.array([[0.0, 0, 0.0], [1.2, 0, 0.1], [1.0, 0, 1.0], [0.1, 0, 0.8]])
    order = [0, 3, 2, 1]
    cw = ccw[order]
    K1, M1 = plane.matrices(ccw, mat, incompatible=inc)
    K2, M2 = plane.matrices(cw, mat, incompatible=inc)
    p = np.array([[2 * n, 2 * n + 1] for n in order]).ravel()
    assert np.allclose(K2, K1[np.ix_(p, p)], atol=1e-12 * np.abs(K1).max())
    assert np.allclose(M2, M1[np.ix_(p, p)], atol=1e-14 * np.abs(M1).max())
    assert np.allclose(plane.recovery(cw, mat, incompatible=inc), plane.recovery(ccw, mat, incompatible=inc)[:, p])


# ======================================================================================
# BEAMS
# ======================================================================================
_E, _NU, _RHO = 2.0e5, 0.3, 7.8
_SEC = dict(A=0.02, As2=0.015, As3=0.012, J=3.0e-5, I2=4.0e-5, I3=1.0e-4)


def _bmat():
    return material_from_E_nu(_E, _NU, _RHO)


def test_beam_cantilever_flexibility_all_six_components():
    """Clamped at I, unit load at J in each local direction: the exact Timoshenko cantilever
    flexibility including the sign of the rotation about axis 2 (theta2 = -du3/dx) and the
    coupling terms, axial and torsion (spec 08 4.5, D-ELM-04)."""
    mat = _bmat()
    G = mat.G0
    L = 4.0
    xyz = np.array([[0, 0, 0], [L, 0, 0], [0, 1, 0]], dtype=float)    # local = global
    K = beam.matrices(xyz, mat, section=_SEC)[0].real
    Fl = np.linalg.inv(K[6:, 6:])
    s = _SEC
    exp = np.zeros((6, 6))
    exp[0, 0] = L / (_E * s["A"])
    exp[3, 3] = L / (G * s["J"])
    exp[1, 1] = L ** 3 / (3 * _E * s["I3"]) + L / (G * s["As2"])
    exp[5, 5] = L / (_E * s["I3"])
    exp[1, 5] = exp[5, 1] = L ** 2 / (2 * _E * s["I3"])               # +P2 -> +theta3
    exp[2, 2] = L ** 3 / (3 * _E * s["I2"]) + L / (G * s["As3"])
    exp[4, 4] = L / (_E * s["I2"])
    exp[2, 4] = exp[4, 2] = -L ** 2 / (2 * _E * s["I2"])              # +P3 -> -theta2
    assert np.allclose(Fl, exp, rtol=1e-10, atol=1e-14 * np.abs(exp).max())


def test_beam_rigid_body_modes_arbitrary_orientation():
    """The six rigid motions u = a + w x r (r = node position) are in the null space of the global
    matrix for an oblique beam with shear areas; the rest of the spectrum is positive."""
    mat = _bmat()
    R = _rotation(4)
    xyz = np.vstack([np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]) + R @ np.array([3.0, 0, 0]),
                     np.array([1.0, 2.0, 3.0]) + R @ np.array([0.5, 2.0, 0.3])])
    K = beam.matrices(xyz, mat, section=_SEC)[0].real
    for k in range(6):
        a = np.zeros(3)
        w = np.zeros(3)
        (a if k < 3 else w)[k % 3] = 1.0
        u = np.concatenate([np.concatenate([a + np.cross(w, xyz[n]), w]) for n in (0, 1)])
        assert np.abs(K @ u).max() <= 1e-10 * np.abs(K).max()
    ev = np.linalg.eigvalsh(K)
    assert int((np.abs(ev) < 1e-10 * ev.max()).sum()) == 6


@pytest.mark.parametrize("kj, idx, I, As", [("000001", 7, "I3", "As2"), ("000010", 8, "I2", "As3")])
def test_beam_release_with_shear_gives_propped_timoshenko_stiffness(kj, idx, I, As):
    """KJ releasing the bending moment: the lateral stiffness at J of a beam clamped at I is
    1 / (L^3/(3EI) + L/(G As)) - static condensation must keep the shear flexibility."""
    mat = _bmat()
    L = 4.0
    xyz = np.array([[0, 0, 0], [L, 0, 0], [0, 1, 0]], dtype=float)
    K = beam.matrices(xyz, mat, section=_SEC, kj=kj)[0].real
    ref = 1.0 / (L ** 3 / (3 * _E * _SEC[I]) + L / (mat.G0 * _SEC[As]))
    assert K[idx, idx] == pytest.approx(ref, rel=1e-12)


def test_beam_recovery_cantilever_end_forces_on_element():
    """Cantilever with tip load P along e2 (oblique beam): forces exerted on the element (D-STR-05)
    are FYJ = +P, FYI = -P, MZI = -P L, MZJ = 0; all other components zero."""
    mat = _bmat()
    L, P = 4.0, 3.0
    R = _rotation(8)
    xyz = np.vstack([np.zeros(3), R @ np.array([L, 0, 0]), R @ np.array([1.0, 1.0, 0])])
    K = beam.matrices(xyz, mat, section=_SEC)[0].real
    e2 = R[:, 1]
    uJ = np.linalg.solve(K[6:, 6:], np.concatenate([P * e2, np.zeros(3)]))
    f = (beam.recovery(xyz, mat, section=_SEC) @ np.concatenate([np.zeros(6), uJ])).real
    exp = np.zeros(12)
    exp[1], exp[7], exp[5] = -P, P, -P * L
    assert np.allclose(f, exp, atol=1e-9 * P * L)


def test_beam_consistent_mass_rigid_motions_with_shear():
    """Kinetic energies of rigid motions (D-ELM-05): translations rho A L, twist rho (I2+I3) L,
    rigid rotation about the transverse axes through I rho A L^3 / 3 (no rotary inertia) - the
    Timoshenko interpolation reproduces v = theta x exactly."""
    mat = _bmat()
    L = 4.0
    R = _rotation(9)
    xyz = np.vstack([np.zeros(3), R @ np.array([L, 0, 0]), R @ np.array([0.3, 1.0, 0])])
    M = beam.matrices(xyz, mat, section=_SEC)[1]
    T = np.kron(np.eye(4), R.T)                                        # local = R^T global
    ml = T @ M @ T.T                                                    # back to local axes
    rA = _RHO * _SEC["A"]
    for d in range(3):
        u = np.zeros(12)
        u[d] = u[6 + d] = 1.0
        assert u @ ml @ u == pytest.approx(rA * L, rel=1e-12)
    u = np.zeros(12)
    u[3] = u[9] = 1.0
    assert u @ ml @ u == pytest.approx(_RHO * (_SEC["I2"] + _SEC["I3"]) * L, rel=1e-12)
    u = np.zeros(12)
    u[5], u[7], u[11] = 1.0, L, 1.0                                     # rotation about e3
    assert u @ ml @ u == pytest.approx(rA * L ** 3 / 3.0, rel=1e-12)
    u = np.zeros(12)
    u[4], u[8], u[10] = 1.0, -L, 1.0                                    # rotation about e2
    assert u @ ml @ u == pytest.approx(rA * L ** 3 / 3.0, rel=1e-12)


def test_beam_complex_moduli_from_M_and_G_when_dampings_differ():
    """E*, G* from (M*, G*) (requirements 4.0.2): axial stiffness E* A / L with
    E* = 2 G* (1 + nu*), nu* = (M* - 2G*)/(2(M* - G*)); torsion G* J / L."""
    mat = material_from_M(1, 2.0e5, 0.3, 7.8, 0.02, 0.07, 1.0)
    Gs, Ms_, _ = complex_lame(mat.G0, mat.M0, 0.07, 0.02)
    nus = (Ms_ - 2 * Gs) / (2 * (Ms_ - Gs))
    Es = 2 * Gs * (1 + nus)
    L = 2.5
    xyz = np.array([[0, 0, 0], [L, 0, 0], [0, 0, 1.0]], dtype=float)
    K = beam.matrices(xyz, mat, section=_SEC)[0]
    assert K[0, 0] == pytest.approx(Es * _SEC["A"] / L, rel=1e-12)
    assert K[3, 3] == pytest.approx(Gs * _SEC["J"] / L, rel=1e-12)
    assert abs(Es / (mat.E0 * cfactor(0.02)) - 1) > 1e-3                # E* is not E c(beta_p)


def test_general_with_beam_local_matrices_equals_beams():
    """A 3-node GENERAL element given the BEAMS local 12x12 matrices (K_R, K_I, M) in the I-J-K
    axes must reproduce the BEAMS element (same local-axis construction and T, spec 08 4.9)."""
    mat = material_from_M(1, 2.0e5, 0.3, 7.8, 0.05, 0.05, 1.0)
    R = _rotation(12)
    xyz = np.vstack([np.array([1.0, 0, 2.0]), np.array([1.0, 0, 2.0]) + R @ np.array([3.0, 0, 0]),
                     np.array([1.0, 0, 2.0]) + R @ np.array([1.0, 2.0, 0.0])])
    sec = dict(_SEC)
    T, _, Lb = beam.transformation(xyz)
    kl = beam.local_stiffness(Lb, *mat.beam_moduli(), sec)
    ml = beam.local_mass(Lb, mat.rho, mat.E0, mat.G0, sec)
    Kg, Mg = general.matrices(xyz, None, KR=kl.real, KI=kl.imag, MM=ml)
    Kb, Mb = beam.matrices(xyz, mat, section=sec)
    assert np.abs(Kg - Kb).max() <= 1e-12 * np.abs(Kb).max()
    assert np.abs(Mg - Mb).max() <= 1e-12 * np.abs(Mb).max()


def test_rectangle_section_orientation_swap():
    """J of the rectangle is symmetric in (b, h) (Roark with b <= h, spec 08 Q16); I2 and I3 swap."""
    s1, s2 = section_rectangle(0.5, 1.0), section_rectangle(1.0, 0.5)
    assert s1["J"] == pytest.approx(s2["J"], rel=1e-14)
    assert s1["I2"] == pytest.approx(s2["I3"]) and s1["I3"] == pytest.approx(s2["I2"])
    assert s1["J"] == pytest.approx(0.028610, abs=5e-7)


# ======================================================================================
# SHELL
# ======================================================================================
def _strip_model(nx: int, L: float, b: float, mat, t: float):
    node_xyz = {}

    def nid(i, j):
        return 1 + i + j * (nx + 1)

    for j in range(2):
        for i in range(nx + 1):
            node_xyz[nid(i, j)] = (L * i / nx, b * j, 0.0)
    recs = [ElemRecord(1, i + 1, 3, (nid(i, 0), nid(i + 1, 0), nid(i + 1, 1), nid(i, 1)), False, mat,
                       dict(thick=t)) for i in range(nx)]
    return node_xyz, recs, nid


def test_dkq_cantilever_strip_exact_and_moment_sign():
    """nu = 0 strip clamped at x = 0 with a tip line load: w is cubic, which DKQ represents
    exactly (beta_n is linear on every side), so w_tip = P L^3/(3 E I) for any nx.  The recovered
    Mx'x' = -P (L - x_c)/b (positive = tension on the +z' face, the face compressed here)."""
    E, t = 1.0e4, 0.1
    mat = material_from_E_nu(E, 0.0, 2.0)
    L, b, nx, P = 10.0, 1.0, 4, 1.0
    node_xyz, recs, nid = _strip_model(nx, L, b, mat, t)
    ids = sorted(node_xyz)
    fix = np.array([[1] * 6 if node_xyz[n][0] == 0 else [1, 1, 0, 0, 0, 1] for n in ids])
    dm = build_dofmap(ids, fix, recs)
    am = assemble(node_xyz, recs, dm)
    F = np.zeros(dm.neq)
    for j in range(2):
        F[dm.eq(nid(nx, j), 3)] = P / 2
    u = spla.spsolve(sp.csc_matrix(am.Ks.real), F)
    I = b * t ** 3 / 12
    assert u[dm.eq(nid(nx, 0), 3)] == pytest.approx(P * L ** 3 / (3 * E * I), rel=1e-9)
    rec = am.recovery["SHELL"]
    ue = np.where(rec["eq"] >= 0, u[np.maximum(rec["eq"], 0)], 0.0)
    s = np.einsum("ecd,ed->ec", rec["S"].real, ue)
    xc = L * (np.arange(nx) + 0.5) / nx
    assert np.allclose(s[:, 3], -P * (L - xc) / b, rtol=1e-8)
    assert np.allclose(s[:, [0, 1, 2, 4, 5]], 0.0, atol=1e-9 * P * L)


def test_dkq_simply_supported_plate_navier_deflection():
    """Uniformly loaded simply supported square plate: w_c = 0.00406235 q a^4 / D (Navier); DKQ
    on a 16 x 16 mesh with lumped loads is within 0.05 % (converges, no residual bias)."""
    E, nu, t, a = 200e9, 0.3, 0.05, 10.0
    mat = material_from_E_nu(E, nu, 8000.0)
    D = E * t ** 3 / (12 * (1 - nu * nu))
    n = 16
    node_xyz = {}

    def nid(i, j):
        return 1 + i + j * (n + 1)

    for j in range(n + 1):
        for i in range(n + 1):
            node_xyz[nid(i, j)] = (a * i / n, a * j / n, 0.0)
    recs = [ElemRecord(1, i + j * n + 1, 3, (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)), False, mat,
                       dict(thick=t)) for j in range(n) for i in range(n)]
    ids = sorted(node_xyz)

    def on_edge(x, y):
        return min(abs(x), abs(x - a), abs(y), abs(y - a)) < 1e-9

    fix = np.array([[1, 1, 1 if on_edge(*node_xyz[k][:2]) else 0, 0, 0, 1] for k in ids])
    dm = build_dofmap(ids, fix, recs)
    am = assemble(node_xyz, recs, dm)
    F = np.zeros(dm.neq)
    h = a / n
    for k in ids:
        e = dm.eq(k, 3)
        if e >= 0:
            F[e] = h * h                                              # q = 1, interior nodes only
    u = spla.spsolve(sp.csc_matrix(am.Ks.real), F)
    wc = u[dm.eq(nid(n // 2, n // 2), 3)] * D / a ** 4
    assert wc == pytest.approx(0.00406235, rel=5e-4)


@pytest.mark.parametrize("tri", [False, True])
def test_shell_objectivity_null_space_and_mass(tri):
    """Oblique flat facet: K(R X) = T K(X) T^T; null space = 6 rigid motions + n drilling
    rotations (D-ELM-06, zero drilling stiffness); total translational mass rho t A (D-ELM-07)."""
    mat = material_from_M(1, 1.0e4, 0.3, 2.0, 0.05, 0.05, 1.0)
    t = 0.2
    flat = np.array([[0.0, 0, 0], [2.0, 0.2, 0], [1.8, 1.5, 0], [-0.1, 1.2, 0]])
    if tri:
        flat = flat[:3]
    n = flat.shape[0]
    A = abs(shell.polygon_area(flat[:, :2]))
    R = _rotation(17)
    X = flat @ R.T + np.array([1.0, 2.0, 3.0])
    K0, M0 = shell.matrices(flat, mat, thick=t)
    K1, M1 = shell.matrices(X, mat, thick=t)
    T = np.kron(np.eye(2 * n), R)
    assert np.abs(K1 - T @ K0 @ T.T).max() <= 1e-10 * np.abs(K0).max()
    assert np.allclose(M1, M0)
    ev = np.linalg.eigvalsh(K1.real / cfactor(0.05).real)
    assert int((np.abs(ev) < 1e-9 * ev.max()).sum()) == 6 + n
    for d in range(3):
        e = np.zeros(6 * n)
        e[d::6] = 1.0
        assert e @ M1 @ e == pytest.approx(mat.rho * t * A, rel=1e-12)
    # drilling: a pure rotation about the normal at one node carries no force
    znode = np.zeros(6 * n)
    znode[3:6] = R[:, 2]
    assert np.abs(K1 @ znode).max() <= 1e-10 * np.abs(K1).max()


def test_shell_membrane_inplane_bending_exact_for_rectangle():
    """Plane-stress in-plane bending of a rectangular facet: u = -k x y, v = k/2 (x^2 + nu y^2)
    (sigma_yy = 0).  The Q4 + incompatible modes membrane gives the exact energy
    1/2 E t k^2 a b^3/12 (D-ELM-06)."""
    E, nu, t = 1.0e4, 0.3, 0.1
    mat = material_from_E_nu(E, nu, 1.0)
    a, b, k = 2.0, 1.0, 1e-3
    X = np.array([[0, 0, 0], [a, 0, 0], [a, b, 0], [0, b, 0]], dtype=float)
    c = np.array([a / 2, b / 2])
    u = np.zeros(24)
    for i, p in enumerate(X):
        x, y = p[0] - c[0], p[1] - c[1]
        u[6 * i] = -k * x * y
        u[6 * i + 1] = 0.5 * k * (x * x + nu * y * y)
    K = shell.matrices(X, mat, thick=t)[0].real
    U = 0.5 * E * t * k * k * a * b ** 3 / 12.0
    assert 0.5 * u @ K @ u == pytest.approx(U, rel=1e-12)


# ======================================================================================
# SPRING, GENERAL
# ======================================================================================
def test_spring_zero_length_and_force_sign():
    """Zero-length springs are allowed (requirements 4.1); F = k* (u_J - u_I) (D-STR-05)."""
    k = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    K, M = spring.matrices(np.zeros((2, 3)), None, k=k, damp=0.1)
    assert np.all(M == 0)
    u = np.arange(12, dtype=float)
    f = spring.recovery(None, None, k=k, damp=0.1) @ u
    assert np.allclose(f, np.array(k) * cfactor(0.1) * 6.0)
    # nodal forces K u are -f at I and +f at J
    assert np.allclose(K @ u, np.concatenate([-f, f]))
    with pytest.raises(ElementError):
        spring.matrices(None, None, k=(1, 1, -1, 0, 0, 0))


# ======================================================================================
# Assembly
# ======================================================================================
def _mixed_model():
    mat = material_from_M(1, 1000.0, 0.3, 2.0, 0.03, 0.06, 1.0)
    mshell = material_from_M(1, 1000.0, 0.3, 2.0, 0.05, 0.05, 1.0)
    node_xyz = {}
    hexa = (solid.NODE_NAT + 1.0) / 2.0
    for a in range(8):
        node_xyz[10 + a] = tuple(hexa[a])
    node_xyz[30] = (0.5, 0.5, 2.0)          # beam end
    node_xyz[31] = (0.0, 3.0, 1.0)          # beam K node (orientation only)
    node_xyz[40] = (2.0, 0.0, 1.0)          # shell node
    node_xyz[41] = (2.0, 1.0, 1.0)          # shell node
    node_xyz[50] = (5.0, 5.0, 5.0)          # unused node
    node_xyz[60] = (1.0, 0.0, 1.0)          # spring partner of node 15 (coincident allowed elsewhere)
    recs = [
        ElemRecord(1, 1, 1, tuple(range(10, 18)), False, mat, dict(incompatible=True, eint=1)),
        ElemRecord(2, 1, 2, (16, 30, 31), False, mat, dict(section=_SEC, ki="000000", kj="000001")),
        ElemRecord(3, 1, 3, (15, 40, 41, 16), False, mshell, dict(thick=0.1)),
        ElemRecord(3, 2, 3, (40, 41, 30, 0), False, mshell, dict(thick=0.1)),   # triangle, L = 0
        ElemRecord(7, 1, 7, (15, 60), False, None, dict(k=(1, 2, 3, 4, 5, 6), damp=0.04)),
        ElemRecord(9, 1, 9, (60, 40, 31), False, None,
                   dict(KR=np.eye(12) * 7.0, KI=np.eye(12) * 0.5, MM=np.eye(12) * 0.25)),
    ]
    return node_xyz, recs


def test_assembly_equals_hand_scatter_of_element_matrices():
    """Ks/Ms equal a dense hand scatter of every element's (K, M) through DofMap.eq, with fixed
    DOFs dropped; the K node and the unused node get no equations (requirements 4.1)."""
    node_xyz, recs = _mixed_model()
    ids = sorted(node_xyz, reverse=True)                               # numbering follows node_ids
    fix = np.zeros((len(ids), 6), dtype=int)
    fix[ids.index(10)] = 1
    fix[ids.index(11), :3] = 1
    dm = build_dofmap(ids, fix, recs)
    assert dm.eq(31, 1) == -1 and dm.eq(50, 1) == -1                   # K node / unused node
    assert dm.eq(12, 4) == -1 and dm.eq(16, 4) >= 0                    # SOLID-only node: no rotations
    assert np.all(np.diff([dm.eq(n, d) for n in ids for d in range(1, 7) if dm.eq(n, d) >= 0]) == 1)
    am = assemble(node_xyz, recs, dm, masses={30: (1, 1, 1, 0, 0, 0)})
    Kd = np.zeros((dm.neq, dm.neq), dtype=complex)
    Md = np.zeros((dm.neq, dm.neq))
    from sassi.elements.assemble import element_dof_nodes
    from sassi.elements import ELEMENTS
    for r in recs:
        dn, gn = element_dof_nodes(r)
        spec = ELEMENTS[r.code]
        xyz = np.array([node_xyz[n] for n in gn], dtype=float)
        K, M = spec.matrices(xyz, r.mat, **r.props)
        eq = dm.eqs(dn, spec.dofs)
        for i, ei in enumerate(eq):
            for j, ej in enumerate(eq):
                if ei >= 0 and ej >= 0:
                    Kd[ei, ej] += K[i, j]
                    Md[ei, ej] += M[i, j]
    Md[dm.eq(30, 1), dm.eq(30, 1)] += 1
    Md[dm.eq(30, 2), dm.eq(30, 2)] += 1
    Md[dm.eq(30, 3), dm.eq(30, 3)] += 1
    assert np.abs(am.Ks.toarray() - Kd).max() <= 1e-12 * np.abs(Kd).max()
    assert np.abs(am.Ms.toarray() - Md).max() <= 1e-12 * np.abs(Md).max()
    assert am.Ke.nnz == 0 and am.Me.nnz == 0
    # complex-symmetric (not Hermitian) global stiffness
    assert np.abs(am.Ks - am.Ks.T).max() <= 1e-12 * np.abs(Kd).max()


def test_degenerate_solid_assembles_as_merged_nodes():
    """A wedge given with repeated node numbers (5 = 6, 7 = 8) assembles to T^T K8 T where T merges
    the repeated nodes; its rigid translation mass is rho V (spec 08 4.4)."""
    mat = material_from_E_nu(1000.0, 0.3, 2.0)
    X = (solid.NODE_NAT + 1.0) / 2.0
    X = _collapse(X, "wedge")
    nodes = (1, 2, 3, 4, 5, 5, 6, 6)
    node_xyz = {1: X[0], 2: X[1], 3: X[2], 4: X[3], 5: X[4], 6: X[6]}
    recs = [ElemRecord(1, 1, 1, nodes, False, mat, dict(incompatible=True))]
    dm = build_dofmap(sorted(node_xyz), None, recs)
    am = assemble(node_xyz, recs, dm)
    K8, M8 = solid.matrices(X, mat, incompatible=True)
    T = np.zeros((24, dm.neq))
    for a, n in enumerate(nodes):
        for d in range(3):
            T[3 * a + d, dm.eq(n, d + 1)] = 1.0
    assert np.allclose(am.Ks.toarray(), T.T @ K8 @ T, atol=1e-12 * np.abs(K8).max())
    e = np.zeros(dm.neq)
    e[0::3] = 1.0
    assert e @ am.Ms.toarray() @ e == pytest.approx(2.0 * 0.5, rel=1e-12)


def test_excavated_soil_goes_to_Ke_without_incompatible_modes():
    """Excavated SOLID/PLANE: matrices only in Ke/Me, compatible element even if the record asks
    for incompatible modes (D-ELM-02/03), structural copy only in Ks/Ms."""
    mat = material_from_E_nu(1000.0, 0.3, 2.0)
    X = _distorted_hex(11)
    node_xyz = {i + 1: X[i] for i in range(8)}
    rec_e = ElemRecord(1, 1, 1, tuple(range(1, 9)), True, mat, dict(incompatible=True))
    rec_s = ElemRecord(2, 1, 1, tuple(range(1, 9)), False, mat, dict(incompatible=True))
    dm = build_dofmap(sorted(node_xyz), None, [rec_e, rec_s])
    am = assemble(node_xyz, [rec_e, rec_s], dm)
    Kc, Mc = solid.matrices(X, mat, incompatible=False)
    Ki, _ = solid.matrices(X, mat, incompatible=True)
    assert np.allclose(am.Ke.toarray(), Kc, atol=1e-12 * np.abs(Kc).max())
    assert np.allclose(am.Me.toarray(), Mc)
    assert np.allclose(am.Ks.toarray(), Ki, atol=1e-12 * np.abs(Ki).max())
    assert list(am.recovery["SOLID"]["excavated"]) == [1, 0]


def test_excavated_flag_on_non_soil_elements_is_not_subtracted():
    """A buried SHELL (deck etype 2) or any non-continuum element flagged excavated must not end up
    in the excavated-soil matrices Ke/Me (they are subtracted in the flexible-volume impedance);
    the assembler should keep it in Ks or reject the record."""
    mat = material_from_E_nu(1000.0, 0.3, 2.0)
    node_xyz = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0)}
    recs = [ElemRecord(1, 1, 3, (1, 2, 3, 4), True, mat, dict(thick=0.1))]
    dm = build_dofmap(sorted(node_xyz), None, recs)
    try:
        am = assemble(node_xyz, recs, dm)
    except ElementError:
        return                                                         # rejecting is acceptable
    assert am.Ke.nnz == 0 and am.Ks.nnz > 0
