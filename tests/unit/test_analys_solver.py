"""Unit tests of sassi.core.ssi_solver (ANALYS numerics: R1 4.2-4.6, requirements 4.6, D-ANL-03/04/05/07)."""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from sassi.core import ssi_solver as ss
from sassi.io.container import read_container, write_container


def _cs_sparse(n: int, seed: int, density: float = 0.08, shift: float = 8.0) -> sp.csr_matrix:
    """Random sparse complex *symmetric* (not Hermitian) matrix, diagonally dominated."""
    rng = np.random.default_rng(seed)
    A = sp.random(n, n, density=density, random_state=rng, format="csr")
    B = sp.random(n, n, density=density, random_state=rng, format="csr")
    C = (A + 1j * 0.1 * B)
    C = C + C.T + sp.diags(shift + 0.3j + rng.random(n))
    return sp.csr_matrix(C)


def _cs_dense(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n)) + 0.2j * rng.standard_normal((n, n))
    return A + A.T + n * np.eye(n)


def _system(neq=40, nint=6, seed=3):
    """C (sparse), X (3 nint x 3 nint dense) and a partition with the interaction DOFs scattered."""
    rng = np.random.default_rng(seed)
    eqs = rng.permutation(neq)[:3 * nint].reshape(nint, 3)
    eqs[1, 2] = -1                                  # one fixed interaction translation
    part = ss.Partition.from_int_eq(eqs, neq)
    C = _cs_sparse(neq, seed)
    X = _cs_dense(3 * nint, seed + 1)
    return C, X, part


def _dense_reference(C, X, part, b):
    A = C.toarray().astype(complex)
    A[np.ix_(part.f_eq, part.f_eq)] += X[np.ix_(part.x_pos, part.x_pos)]
    return np.linalg.solve(A, b)


def test_partition_from_int_eq():
    eqs = np.array([[4, 5, 6], [0, -1, 2]])
    p = ss.Partition.from_int_eq(eqs, 8)
    np.testing.assert_array_equal(p.f_eq, [4, 5, 6, 0, 2])
    np.testing.assert_array_equal(p.x_pos, [0, 1, 2, 3, 5])
    np.testing.assert_array_equal(p.n_eq, [1, 3, 7])
    with pytest.raises(ValueError):
        ss.Partition.from_int_eq(np.array([[1, 1, 2]]), 4)


@pytest.mark.parametrize("seismic", [True, False])
def test_schur_solution_equals_dense_solve(seismic):
    """D-ANL-03: Schur complement (sparse LU of C_nn, dense LU of S + X) = direct dense solve."""
    C, X, part = _system()
    fac = ss.factorize(C, X[np.ix_(part.x_pos, part.x_pos)], part)
    assert fac.path == "schur" and fac.rcond_nn > ss.RCOND_MIN and fac.rcond_ff > 0
    rng = np.random.default_rng(7)
    b = np.zeros((part.neq, 2), complex)
    if seismic:
        Up = rng.standard_normal((X.shape[0], 2)) + 1j * rng.standard_normal((X.shape[0], 2))
        bf = ss.seismic_load(X, Up, part)
        b[part.f_eq] = bf
        U = fac.solve(bf)
    else:
        b = rng.standard_normal((part.neq, 2)) + 1j * rng.standard_normal((part.neq, 2))
        U = fac.solve(b[part.f_eq], b[part.n_eq])
    ref = _dense_reference(C, X, part, b)
    np.testing.assert_allclose(U, ref, rtol=0, atol=1e-11 * np.abs(ref).max())


def test_dense_path_when_every_equation_interacts():
    nint = 5
    C = _cs_sparse(3 * nint, 2)
    X = _cs_dense(3 * nint, 4)
    part = ss.Partition.from_int_eq(np.arange(3 * nint).reshape(nint, 3), 3 * nint)
    fac = ss.factorize(C, X, part)
    assert fac.path == "dense"
    b = np.random.default_rng(1).standard_normal(3 * nint) + 0j
    np.testing.assert_allclose(fac.solve(b)[:, 0], np.linalg.solve(C.toarray() + X, b), rtol=1e-11, atol=1e-13)


def test_full_fallback_path_and_forced_threshold():
    """Requirements 4.6 item 6: when C_nn does not meet the rcond threshold the full system with the
    dense X block is factorised by a sparse LU; a DOF without any stiffness makes the full system
    singular too (SolverError)."""
    C, X, part = _system()
    C = C.tolil()
    k = int(part.n_eq[0])                      # a non-interaction DOF without stiffness or mass
    C[k, :] = 0
    C[:, k] = 0
    with pytest.raises(ss.SolverError, match="singular"):
        ss.factorize(C.tocsr(), X[np.ix_(part.x_pos, part.x_pos)], part)
    C2, X2, part2 = _system(seed=5)
    fac2 = ss.factorize(C2, X2[np.ix_(part2.x_pos, part2.x_pos)], part2, rcond_min=2.0)
    assert fac2.path == "full" and "full sparse LU" in fac2.note
    b2 = np.random.default_rng(2).standard_normal((part2.neq, 1)) + 0j
    ref = _dense_reference(C2, X2, part2, b2)
    np.testing.assert_allclose(fac2.solve(b2[part2.f_eq], b2[part2.n_eq]), ref, rtol=0, atol=1e-11 * np.abs(ref).max())


def test_singular_column_coupled_only_to_soil_goes_full():
    """A non-interaction DOF that is only stiffened through an interaction DOF is fine in the full
    path; C_nn of that DOF alone is singular."""
    C, X, part = _system(seed=9)
    C = C.tolil()
    k = int(part.n_eq[1])
    j = int(part.f_eq[0])
    C[k, :] = 0
    C[:, k] = 0
    C[k, k] = 1.0
    C[k, j] = C[j, k] = -1.0
    C[j, j] = C[j, j] + 1.0
    fac = ss.factorize(C.tocsr(), X[np.ix_(part.x_pos, part.x_pos)], part)
    b = np.random.default_rng(4).standard_normal((part.neq, 1)) + 0j
    ref = _dense_reference(C.tocsr(), X, part, b)
    np.testing.assert_allclose(fac.solve(b[part.f_eq], b[part.n_eq]), ref, rtol=0, atol=1e-10 * np.abs(ref).max())


def test_no_interaction_dofs():
    C = _cs_sparse(12, 1)
    part = ss.Partition.from_int_eq(np.zeros((0, 3), int), 12)
    fac = ss.factorize(C, np.zeros((0, 0)), part)
    b = np.arange(12.0) + 0j
    np.testing.assert_allclose(fac.solve(np.zeros((0, 1)), b)[:, 0], np.linalg.solve(C.toarray(), b), rtol=1e-12)


@pytest.mark.parametrize("seed,neq,nint", [(3, 40, 6), (11, 30, 10)])
def test_restart_record_round_trip(tmp_path, seed, neq, nint):
    """D-ANL-05: the factorisation stored in a COOTK container and restored solves identically
    (restored SuperLU factors are applied by triangular solves)."""
    C, X, part = _system(neq, nint, seed)
    fac = ss.factorize(C, X[np.ix_(part.x_pos, part.x_pos)], part)
    write_container(tmp_path / "COOTK001", "COOTK", fac.arrays(), {"rcond_nn": fac.rcond_nn, "rcond_ff": fac.rcond_ff})
    back = ss.SSIFactor.from_arrays(read_container(tmp_path / "COOTK001", "COOTK"), neq)
    rng = np.random.default_rng(seed)
    b = rng.standard_normal((neq, 3)) + 1j * rng.standard_normal((neq, 3))
    U1 = fac.solve(b[part.f_eq], b[part.n_eq])
    U2 = back.solve(b[part.f_eq], b[part.n_eq])
    np.testing.assert_allclose(U2, U1, rtol=0, atol=1e-12 * np.abs(U1).max())
    assert back.path == fac.path and back.rcond_nn == fac.rcond_nn


def test_restart_record_dense_path(tmp_path):
    nint = 4
    C = _cs_sparse(3 * nint, 6)
    X = _cs_dense(3 * nint, 7)
    part = ss.Partition.from_int_eq(np.arange(3 * nint).reshape(nint, 3), 3 * nint)
    fac = ss.factorize(C, X, part)
    write_container(tmp_path / "TK", "COOTK", fac.arrays(), {})
    back = ss.SSIFactor.from_arrays(read_container(tmp_path / "TK"), part.neq)
    b = np.arange(3 * nint) + 1j
    assert back.path == "dense" and np.array_equal(back.solve(b), fac.solve(b))


def test_restart_record_full_path(tmp_path):
    C, X, part = _system(seed=13)
    fac = ss.factorize(C, X[np.ix_(part.x_pos, part.x_pos)], part, rcond_min=2.0)
    write_container(tmp_path / "TK", "COOTK", fac.arrays(), {})
    back = ss.SSIFactor.from_arrays(read_container(tmp_path / "TK"), part.neq)
    b = np.random.default_rng(0).standard_normal((part.neq, 1)) + 0j
    np.testing.assert_allclose(back.solve(b[part.f_eq], b[part.n_eq]), fac.solve(b[part.f_eq], b[part.n_eq]),
                               rtol=1e-12, atol=1e-14)


def test_multiple_rhs_equal_single_rhs_bitwise():
    """Every simultaneous case gets bit-for-bit the single-case result (VP-23)."""
    C, X, part = _system(seed=21)
    fac = ss.factorize(C, X[np.ix_(part.x_pos, part.x_pos)], part)
    rng = np.random.default_rng(3)
    Up = rng.standard_normal((X.shape[0], 3)) + 1j * rng.standard_normal((X.shape[0], 3))
    U = fac.solve(ss.seismic_load(X, Up, part))
    for j in range(3):
        Uj = fac.solve(ss.seismic_load(X, Up[:, j], part))
        assert np.array_equal(U[:, j], Uj[:, 0])


def test_impedance_matrix_is_symmetric_inverse():
    F = _cs_dense(12, 8) * 1e-6
    F = 0.5 * (F + F.T)
    X, rc = ss.impedance_matrix(F)
    with np.errstate(all="ignore"):                       # spurious Accelerate matmul flags (numpy 2.0)
        np.testing.assert_allclose(X @ F, np.eye(12), atol=1e-10)
    assert np.array_equal(X, X.T) and 0 < rc <= 1
    with pytest.raises(ss.SolverError):
        ss.impedance_matrix(np.ones((3, 3)))
    with pytest.raises(ss.SolverError):
        ss.impedance_matrix(np.array([[np.nan]]))


def test_dynamic_matrix():
    Ks = sp.csr_matrix(np.array([[2.0 + 0.1j, -1.0], [-1.0, 1.0]]))
    Ms = sp.csr_matrix(np.diag([1.0, 2.0]))
    Ke = sp.csr_matrix(np.array([[0.5, 0.0], [0.0, 0.0]]))
    Me = sp.csr_matrix(np.diag([0.25, 0.0]))
    C = ss.dynamic_matrix(Ks, Ms, Ke, Me, 3.0).toarray()
    np.testing.assert_allclose(C, Ks.toarray() - Ke.toarray() - 9.0 * (Ms.toarray() - Me.toarray()))


def test_rcond_estimate_matches_exact_order():
    """Hager/Higham estimate of ||A^-1||_1 is a lower bound, usually within a small factor: the
    estimated rcond lies in [exact, 10 exact] for an ill-conditioned (n^2) tridiagonal matrix."""
    n = 200
    A = sp.diags([-np.ones(n - 1), (2.0 + 1e-6j) * np.ones(n), -np.ones(n - 1)], [-1, 0, 1], format="csc")
    A = A + sp.csc_matrix(_cs_sparse(n, 4, density=0.005, shift=0.0) * 1e-3)
    fac = ss.SparseFactor.factor(A)
    est = ss.sparse_rcond(A, fac)
    Ad = A.toarray()
    exact = 1.0 / (np.linalg.norm(Ad, 1) * np.linalg.norm(np.linalg.inv(Ad), 1))
    assert exact < 1e-3
    assert exact * (1 - 1e-12) <= est <= 10 * exact


def test_rigid_body_transform_and_global_impedance():
    """K_G = T^T X T (requirements 4.6 item 9): two points at x = +-d with isotropic point springs k
    give K_xx = K_yy = K_zz = 2k, K_rx = 0 (points on the x axis), K_ry = K_rz = 2 k d^2."""
    d, k = 1.5, 3.0
    xyz = np.array([[d, 0.0, 0.0], [-d, 0.0, 0.0]])
    T = ss.rigid_body_transform(xyz, (0.0, 0.0, 0.0))
    np.testing.assert_allclose(T[0:3], [[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, d], [0, 0, 1, 0, -d, 0]])
    KG = ss.global_impedance(k * np.eye(6), T)
    np.testing.assert_allclose(np.diag(KG).real, [2 * k, 2 * k, 2 * k, 0.0, 2 * k * d * d, 2 * k * d * d])
    np.testing.assert_allclose(KG, KG.T)
    # reference point below the mat: horizontal-rocking coupling K_x,ry = 2 k dz
    KG2 = ss.global_impedance(k * np.eye(6), ss.rigid_body_transform(xyz, (0.0, 0.0, -2.0)))
    assert abs(KG2[0, 4] - 2 * k * 2.0) < 1e-12 and abs(KG2[1, 3] + 2 * k * 2.0) < 1e-12
