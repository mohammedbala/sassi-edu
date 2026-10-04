"""Unit tests of the BEAMS element (requirements 4.1, D-ELM-04/05, D-STR-05, spec 08 4.5)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.elements import beam
from sassi.elements.base import ElementError, material_from_E_nu, material_from_M

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")

MAT = material_from_E_nu(2.0e5, 0.3, 7.8)
SEC = dict(A=0.02, As2=0.015, As3=0.012, J=3e-5, I2=4e-5, I3=1e-4)
L = 10.0
XYZ = np.array([[0.0, 0, 0], [L, 0, 0], [0, 1.0, 0]])


def _rot(seed):
    Q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    return Q if np.linalg.det(Q) > 0 else -Q


def test_local_axes_follow_K_node():
    T, Lam, Lb = beam.transformation(np.array([[0, 0, 0], [0, 0, 3.0], [1.0, 5.0, 1.0]]))
    assert Lb == 3.0
    assert np.allclose(Lam[0], [0, 0, 1])          # e1 = I -> J
    assert np.allclose(Lam[1], [1, 5, 0] / np.sqrt(26))   # e2 toward K, normal to e1
    assert np.allclose(Lam[2], np.cross(Lam[0], Lam[1]))
    with pytest.raises(ElementError, match="Error 9"):
        beam.transformation(np.array([[0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]]))


def test_symmetry_and_rigid_modes():
    K, M = beam.matrices(XYZ, MAT, section=SEC)
    assert np.abs(K - K.T).max() < 1e-12 * np.abs(K).max()
    ev = np.linalg.eigvalsh(K.real)
    assert np.sum(np.abs(ev) < 1e-9 * ev.max()) == 6
    assert np.all(np.linalg.eigvalsh(M) > -1e-12)


def test_cantilever_tip_deflections_exact():
    E, G = MAT.E0, MAT.G0
    K, _ = beam.matrices(XYZ, MAT, section=SEC)
    Kjj = K.real[6:, 6:]
    assert np.linalg.solve(Kjj, np.eye(6)[1])[1] == pytest.approx(
        L ** 3 / (3 * E * SEC["I3"]) + L / (G * SEC["As2"]), rel=1e-12)
    assert np.linalg.solve(Kjj, np.eye(6)[2])[2] == pytest.approx(
        L ** 3 / (3 * E * SEC["I2"]) + L / (G * SEC["As3"]), rel=1e-12)
    assert np.linalg.solve(Kjj, np.eye(6)[0])[0] == pytest.approx(L / (E * SEC["A"]), rel=1e-12)
    assert np.linalg.solve(Kjj, np.eye(6)[3])[3] == pytest.approx(L / (G * SEC["J"]), rel=1e-12)
    eb = dict(SEC, As2=0.0, As3=0.0)
    K, _ = beam.matrices(XYZ, MAT, section=eb)
    assert np.linalg.solve(K.real[6:, 6:], np.eye(6)[1])[1] == pytest.approx(L ** 3 / (3 * E * SEC["I3"]), rel=1e-12)


def test_rotation_invariance_and_k_node_choice():
    R = _rot(4)
    K1, M1 = beam.matrices(XYZ, MAT, section=SEC)
    xyz = (R @ XYZ.T).T
    xyz2 = xyz.copy()
    xyz2[2] = R @ np.array([7.0, 3.0, 0.0])                   # another K in the same 1-2 half plane
    T = np.kron(np.eye(4), R)
    for x in (xyz, xyz2):
        K2, M2 = beam.matrices(x, MAT, section=SEC)
        assert np.allclose(T.T @ K2 @ T, K1, atol=1e-10 * np.abs(K1).max())
        assert np.allclose(T.T @ M2 @ T, M1, atol=1e-12 * np.abs(M1).max())


def test_consistent_mass_totals_and_euler_bernoulli_matrix():
    _, M = beam.matrices(XYZ, MAT, section=SEC)
    rAL = MAT.rho * SEC["A"] * L
    for d in range(3):
        assert M[d::6, d::6].sum() == pytest.approx(rAL, rel=1e-12)          # translational mass
    assert M[3::6, 3::6].sum() == pytest.approx(MAT.rho * (SEC["I2"] + SEC["I3"]) * L, rel=1e-12)
    # no rotary inertia: the bending-rotation terms come from the translational inertia only
    _, Meb = beam.matrices(XYZ, MAT, section=dict(SEC, As2=0.0, As3=0.0))
    assert Meb[5, 5] == pytest.approx(rAL * 4 * L * L / 420, rel=1e-12)
    assert Meb[5, 11] == pytest.approx(-rAL * 3 * L * L / 420, rel=1e-12)
    assert Meb[4, 4] == pytest.approx(rAL * 4 * L * L / 420, rel=1e-12)
    mb = beam._bending_m(420.0, 0.0, 1.0)
    assert np.allclose(mb, [[156, 22, 54, -13], [22, 4, 13, -3], [54, 13, 156, -22], [-13, -3, -22, 4]])
    # shear-deformable mass: Przemieniecki closed form (phi = 0.3)
    phi = 0.3
    mt = beam._bending_m(1.0, phi, 1.0)
    c = 1.0 / (1 + phi) ** 2
    assert mt[0, 0] == pytest.approx(c * (13 / 35 + 7 * phi / 10 + phi ** 2 / 3))
    assert mt[0, 3] == pytest.approx(-c * (13 / 420 + 3 * phi / 40 + phi ** 2 / 24))
    assert mt[1, 1] == pytest.approx(c * (1 / 105 + phi / 60 + phi ** 2 / 120))


def test_end_releases():
    eb = dict(SEC, As2=0.0, As3=0.0)
    E = MAT.E0
    K, _ = beam.matrices(XYZ, MAT, section=eb)
    assert K.real[7, 7] == pytest.approx(12 * E * SEC["I3"] / L ** 3, rel=1e-12)
    Kr, Mr = beam.matrices(XYZ, MAT, section=eb, kj="000001")
    assert np.allclose(Kr[11], 0.0) and np.allclose(Kr[:, 11], 0.0)   # released M3 at J
    assert np.allclose(Mr[11], 0.0)
    # lateral stiffness at J with rotation at J free (condensed): 3 E I / L^3
    keep = [7]
    Kjj = Kr.real[np.ix_(keep, keep)]
    assert Kjj[0, 0] == pytest.approx(3 * E * SEC["I3"] / L ** 3, rel=1e-12)
    # condensed mass keeps the translational total
    assert Mr[1::6, 1::6].sum() == pytest.approx(MAT.rho * SEC["A"] * L, rel=1e-12)
    with pytest.raises(ElementError, match="Error 10"):
        beam.matrices(XYZ, MAT, section=eb, ki=(0, 0, 0, 0, 0, 1), kj="000001")


def test_recovery_end_forces_equilibrium():
    R = _rot(9)
    xyz = (R @ XYZ.T).T
    S = beam.recovery(xyz, MAT, section=SEC)
    u = np.random.default_rng(0).normal(size=12)
    f = (S @ u).real
    Tm, Lam, Lb = beam.transformation(xyz)
    # equilibrium of the element in local axes: forces and moments about I sum to zero
    assert np.allclose(f[0:3] + f[6:9], 0.0, atol=1e-9 * np.abs(f).max())
    mom = f[3:6] + f[9:12] + np.cross([Lb, 0, 0], f[6:9])
    assert np.allclose(mom, 0.0, atol=1e-9 * np.abs(f).max())
    K, _ = beam.matrices(xyz, MAT, section=SEC)
    assert np.allclose(Tm @ (K @ u), S @ u)                      # f_local = T K_global u


def test_complex_stiffness_proportional_when_dampings_equal():
    b = 0.04
    md = material_from_M(1, 2.0e5, 0.3, 7.8, b, b, 1.0)
    K, M = beam.matrices(XYZ, md, section=SEC, ki="000011")
    K0, M0 = beam.matrices(XYZ, md.undamped(), section=SEC, ki="000011")
    assert np.abs(K - K0.real * cfactor(b)).max() < 1e-12 * np.abs(K0).max()
    assert np.isrealobj(M) and np.allclose(M, M0)


def test_section_input_validation():
    with pytest.raises(ElementError):
        beam.matrices(XYZ, MAT, section=dict(SEC, A=0.0))
    with pytest.raises(ElementError):
        beam.matrices(XYZ, MAT)
    K1, _ = beam.matrices(XYZ, MAT, section=[SEC[k] for k in ("A", "As2", "As3", "J", "I2", "I3")])
    K2, _ = beam.matrices(XYZ, MAT, section=SEC)
    assert np.allclose(K1, K2)


@pytest.mark.parametrize("As3", [0.0, 0.012])
@pytest.mark.parametrize("kj", ["000010", "001000", "001010"])
def test_release_in_a_bending_plane_without_stiffness(kj, As3):
    """I2 = 0 is valid input (only I2 < 0 is Error 31): the 1-3 bending plane then has no
    stiffness and a release of P3 or M2 in it is redundant.  The condensation must not fail
    (it used to raise numpy LinAlgError) and the matrices are the continuous limit I2 -> 0+."""
    xyz = (_rot(4) @ XYZ.T).T
    sec0 = dict(SEC, I2=0.0, As3=As3)
    K0, M0 = beam.matrices(xyz, MAT, section=sec0, kj=kj)
    S0 = beam.recovery(xyz, MAT, section=sec0, kj=kj)
    # stiffness: identical to the unreleased element (the released component carries no force)
    Kn, _ = beam.matrices(xyz, MAT, section=sec0)
    assert np.abs(K0 - Kn).max() <= 1e-12 * np.abs(Kn).max()
    rel = beam.releases_to_dofs(None, kj)
    assert np.allclose(S0[rel], 0.0)
    # continuity in I2: matrices for a tiny positive I2 differ by O(I2)
    for eps in (1e-8, 1e-11):
        Ke, Me = beam.matrices(xyz, MAT, section=dict(sec0, I2=eps * SEC["I2"]), kj=kj)
        assert np.abs(K0 - Ke).max() <= 20 * eps * np.abs(Ke).max()
        assert np.abs(M0 - Me).max() <= 20 * eps * np.abs(Me).max()
    # rigid translation carries the full mass rho A L in every direction
    for d in range(3):
        e = np.zeros(12)
        e[d::6] = 1.0
        assert e @ M0 @ e == pytest.approx(MAT.rho * SEC["A"] * L, rel=1e-12)
