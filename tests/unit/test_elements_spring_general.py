"""Unit tests of the SPRING and GENERAL elements (requirements 4.1, 4.0.2, spec 08 4.8-4.9, UT-15)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.conventions import CMOD_SIMPLE, cfactor
from sassi.elements import general, spring
from sassi.elements.base import ElementError

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")

KV = (1.0e3, 2.0e3, 3.0e3, 10.0, 20.0, 30.0)


def test_spring_matrix_and_damping():
    K, M = spring.matrices(None, None, k=KV, damp=0.05)
    c = cfactor(0.05)
    for d in range(6):
        assert K[d, d] == pytest.approx(KV[d] * c)
        assert K[d + 6, d + 6] == pytest.approx(KV[d] * c)
        assert K[d, d + 6] == pytest.approx(-KV[d] * c)
    assert np.count_nonzero(K) == 24 and np.all(M == 0)
    K1, _ = spring.matrices(None, None, k=KV, damp=0.05, form=CMOD_SIMPLE)
    assert K1[0, 0] == pytest.approx(KV[0] * (1 + 0.1j))
    with pytest.raises(ElementError):
        spring.matrices(None, None, k=(1, 2, 3), damp=0.0)
    with pytest.raises(ElementError):
        spring.matrices(None, None, k=(-1, 0, 0, 0, 0, 0))
    with pytest.raises(ElementError):
        spring.matrices(None, None, k=KV, damp=0.6)


def test_spring_recovery_relative_displacement():
    S = spring.recovery(None, None, k=KV, damp=0.0)
    u = np.arange(12, dtype=float)
    assert np.allclose(S @ u, np.asarray(KV) * (u[6:] - u[:6]))


def test_general_from_upper_rows():
    rows = {1: [1000, 0, 0, 0, 0, 0, -1000], 7: [1000]}
    A = general.from_upper_rows(rows)
    assert A[0, 0] == 1000 and A[0, 6] == -1000 and A[6, 0] == -1000 and A[6, 6] == 1000
    assert np.count_nonzero(A) == 4
    tab = np.zeros((12, 12))
    tab[0, [0, 6]] = [1000, -1000]
    tab[6, 0] = 1000
    assert np.allclose(general.from_upper_rows(tab), A)
    with pytest.raises(ElementError):
        general.from_upper_rows({12: [1.0, 2.0]})                  # row 12 holds one term only
    bad = tab.copy()
    bad[11, 5] = 1.0
    with pytest.raises(ElementError):
        general.from_upper_rows(bad)


def test_general_global_equals_spring():
    b = 0.05
    c = cfactor(b)
    KR = np.zeros((12, 12))
    d = np.arange(6)
    for M_, s in ((KR, 1.0),):
        M_[d, d] = KV
        M_[d + 6, d + 6] = KV
        M_[d, d + 6] = -np.asarray(KV)
        M_[d + 6, d] = -np.asarray(KV)
    Kg, Mg = general.matrices(np.zeros((2, 3)), None, KR=KR * c.real, KI=KR * c.imag)
    Ks, _ = spring.matrices(None, None, k=KV, damp=b)
    assert np.abs(Kg - Ks).max() < 1e-12 * max(KV)
    assert np.all(Mg == 0)
    assert general.recovery().shape == (0, 12)


def test_general_local_transformation():
    R, _ = np.linalg.qr(np.random.default_rng(1).normal(size=(3, 3)))
    if np.linalg.det(R) < 0:
        R = -R
    xyz = np.array([np.zeros(3), 2.0 * R[:, 0], R[:, 1] + 0.3 * R[:, 0]])
    A = np.random.default_rng(2).normal(size=(12, 12))
    A = A + A.T
    Kg, Mg = general.matrices(xyz, None, KR=A, MM=np.eye(12))
    T = np.kron(np.eye(4), R.T)                                   # rows of Lam = columns of R
    assert np.allclose(Kg, T.T @ A @ T)
    assert np.allclose(Mg, np.eye(12))                             # isotropic mass is invariant
    with pytest.raises(ElementError, match="Error 9"):
        general.matrices(np.array([[0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]]), None, KR=A)
    with pytest.raises(ElementError):
        general.matrices(None, None, KR=np.triu(A))                # not symmetric


def test_ut15_general_mass_weight_units():
    MM = np.zeros((12, 12))
    MM[6, 6] = MM[7, 7] = MM[8, 8] = 32.2
    _, Mw = general.matrices(None, None, MM=MM, munits=1, gravity=32.2)
    _, Mm = general.matrices(None, None, MM=MM / 32.2, munits=0)
    assert np.array_equal(Mw, Mm)
    with pytest.raises(ElementError):
        general.matrices(None, None, MM=MM, munits=1)


def test_general_deck_table_layout():
    tab = np.zeros((12, 12))
    tab[0, [0, 6]] = [1000.0, -1000.0]          # MXR,p,1,1000,0,0,0,0,0,-1000
    tab[6, 0] = 1000.0                          # MXR,p,7,1000
    K, _ = general.matrices(None, None, KR=tab, upper_rows=True)
    assert K[0, 0] == 1000 and K[6, 0] == -1000 and K[0, 6] == -1000 and K[6, 6] == 1000
    with pytest.raises(ElementError):
        general.matrices(None, None, KR=tab)    # a deck table is not a full symmetric matrix
