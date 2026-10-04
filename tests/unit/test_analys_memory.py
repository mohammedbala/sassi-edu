"""Memory behaviour of the ANALYS numerics and the resource check (D-ANL-10, G-22, EDU-20;
requirements 4.6 item 10).

The review of the ANALYS work package measured about 9 dense ``(3 nInt)^2`` complex matrices alive
per frequency while the resource check budgeted 3.  These tests pin the in-place behaviour that
brought the frequency loop down to two dense matrices plus bounded work arrays, and check with
``tracemalloc`` that the printed estimate bounds the measured peak of complete ANALYS runs.
(numpy registers its array allocations with tracemalloc; SuperLU's internal factor storage is not
traced, which is why the sparse factors have their own term in the estimate.)
"""
from __future__ import annotations

import re
import tracemalloc

import numpy as np
import pytest
import scipy.sparse as sp

import sassi.modules.analys as analys
from sassi.core import ssi_solver as ss
from sassi.verify import builders as B


def _peak_increment(fn, *args, **kwargs):
    """(result, peak traced bytes allocated during the call)."""
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        c0, _ = tracemalloc.get_traced_memory()
        out = fn(*args, **kwargs)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return out, peak - c0


def _cs_dense(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n)) + 0.2j * rng.standard_normal((n, n))
    return A + A.T + n * np.eye(n)


# ======================================================================================
# In-place dense helpers
# ======================================================================================
@pytest.mark.parametrize("block", [1, 7, 64, None])
def test_symmetrize_inplace_equals_the_full_expression_bit_for_bit(block):
    rng = np.random.default_rng(4)
    A = rng.standard_normal((53, 53)) + 1j * rng.standard_normal((53, 53))
    ref = 0.5 * (A + A.T)
    out = ss.symmetrize_inplace(A, block)
    assert out is A
    assert np.array_equal(A, ref)
    with pytest.raises(ValueError):
        ss.symmetrize_inplace(np.zeros((3, 4)))


def test_abs_column_sums_in_blocks(monkeypatch):
    monkeypatch.setattr(ss, "DENSE_BLOCK_ELEMS", 50)
    rng = np.random.default_rng(5)
    A = rng.standard_normal((41, 23)) + 1j * rng.standard_normal((41, 23))
    for M in (A, np.asfortranarray(A), A[::2, ::3]):
        np.testing.assert_allclose(ss.abs_column_sums(M), np.abs(M).sum(axis=0), rtol=1e-14)
    B_ = A.copy()
    B_[3, 7] = np.nan
    assert not np.isfinite(ss.abs_column_sums(B_)[7])


@pytest.mark.parametrize("order", ["C", "F"])
def test_impedance_matrix_inverts_in_place(monkeypatch, order):
    """``impedance_matrix(F, overwrite=True)`` returns X_ff in F's own buffer (getrf + getri; a
    C-ordered F through its Fortran-ordered transpose) and allocates only bounded work arrays."""
    monkeypatch.setattr(ss, "DENSE_BLOCK_ELEMS", 2 ** 12)
    n = 360
    F0 = _cs_dense(n, 6) * 1e-6
    F0 = 0.5 * (F0 + F0.T)
    ref, rc_ref = ss.impedance_matrix(F0)                 # copying variant: F0 untouched
    assert np.array_equal(F0, 0.5 * (_cs_dense(n, 6) * 1e-6 + (_cs_dense(n, 6) * 1e-6).T))
    F = np.array(F0, order=order)
    (X, rc), inc = _peak_increment(ss.impedance_matrix, F, overwrite=True)
    assert np.shares_memory(X, F)
    # bounded work only: tiles / column blocks of DENSE_BLOCK_ELEMS and the getri work array (~64 n entries)
    assert inc < 0.3 * F.nbytes, f"impedance_matrix allocated {inc / F.nbytes:.2f} dense copies"
    assert np.array_equal(X, X.T)
    np.testing.assert_allclose(X, ref, rtol=0, atol=1e-12 * np.abs(ref).max())
    with np.errstate(all="ignore"):                       # spurious Accelerate matmul flags (numpy 2.0)
        np.testing.assert_allclose(X @ F0, np.eye(n), rtol=0, atol=1e-9)
    assert rc == pytest.approx(rc_ref, rel=1e-6)


def test_impedance_rcond_is_the_one_norm_estimate_for_a_c_ordered_matrix():
    """The C-ordered path factorises F^T; gecon with the infinity norm then estimates rcond_1(F)
    (not rcond_inf(F)) -- checked on a matrix whose two condition numbers differ by 500x."""
    n = 60
    rng = np.random.default_rng(1)
    F = np.eye(n) * (1 + rng.random(n)) + 0.01j * rng.standard_normal((n, n))
    F[0, :] += 5.0                                        # one heavy row: ||F||_inf >> ||F||_1
    F = np.ascontiguousarray(F.astype(complex))
    inv = np.linalg.inv(F)
    r1 = 1.0 / (np.abs(F).sum(0).max() * np.abs(inv).sum(0).max())
    rinf = 1.0 / (np.abs(F).sum(1).max() * np.abs(inv).sum(1).max())
    assert r1 > 100 * rinf
    for M in (F.copy(), np.asfortranarray(F)):
        X, rc = ss.impedance_matrix(M, overwrite=True)
        assert r1 * (1 - 1e-12) <= rc <= 3.0 * r1
        np.testing.assert_allclose(X, 0.5 * (inv + inv.T), rtol=0, atol=1e-13 * np.abs(inv).max())


def _schur_system(nint: int, nn: int, seed: int = 7):
    rng = np.random.default_rng(seed)
    neq = 3 * nint + nn
    eqs = rng.permutation(neq)[:3 * nint].reshape(nint, 3)
    part = ss.Partition.from_int_eq(eqs, neq)
    A = sp.random(neq, neq, density=0.004, random_state=rng, format="csr")
    C = A + A.T + sp.diags(8.0 + 0.3j + rng.random(neq))
    X = _cs_dense(3 * nint, seed + 1) * 0.01
    return sp.csr_matrix(C, dtype=complex), X, part


def test_factorize_schur_path_holds_one_dense_matrix(monkeypatch):
    """Schur path: ``S`` is built in Fortran order, ``X_ff`` added in place and ``S + X_ff``
    factorised in place -- one dense matrix plus bounded blocks; ``X`` is not modified and is
    passed without a copy when no interaction translation is fixed (:func:`restrict`)."""
    monkeypatch.setattr(ss, "DENSE_BLOCK_ELEMS", 2 ** 12)
    C, X, part = _schur_system(nint=150, nn=200)
    D = X.nbytes
    X0 = X.copy()
    Xr = ss.restrict(X, part)
    assert Xr is X
    fac, inc = _peak_increment(ss.factorize, C, Xr, part, block_bytes=2 ** 16)
    assert fac.path == "schur"
    assert inc < 1.35 * D, f"factorize allocated {inc / D:.2f} dense copies"
    assert np.array_equal(X, X0)
    assert fac.lu.flags.f_contiguous
    rng = np.random.default_rng(3)
    b = rng.standard_normal(part.neq) + 1j * rng.standard_normal(part.neq)
    A = C.toarray()
    A[np.ix_(part.f_eq, part.f_eq)] += X
    U = fac.solve(b[part.f_eq], b[part.n_eq])[:, 0]
    np.testing.assert_allclose(U, np.linalg.solve(A, b), rtol=0, atol=1e-10 * np.abs(U).max())
    # one interaction translation fixed (dropped from the partition): restrict() copies the active block
    p2 = ss.Partition(part.neq, part.f_eq[:-1], part.x_pos[:-1], np.r_[part.n_eq, part.f_eq[-1]])
    R = ss.restrict(X, p2)
    assert R.shape == (part.nf - 1, part.nf - 1) and not np.shares_memory(R, X)
    np.testing.assert_array_equal(R, X[:-1, :-1])


@pytest.mark.parametrize("block", [5, 2 ** 21])
def test_scatter_dense_block_on_unsorted_equations(monkeypatch, block):
    """The full-LU fallback inserts ``A_f^T X A_f`` built directly in CSC form (sorted equations,
    column blocks): it must equal the dense scatter for any order of the interaction equations."""
    monkeypatch.setattr(ss, "DENSE_BLOCK_ELEMS", block)
    rng = np.random.default_rng(9)
    eqs = rng.permutation(30)[:12].reshape(4, 3)
    part = ss.Partition.from_int_eq(eqs, 30)
    X = rng.standard_normal((12, 12)) + 1j * rng.standard_normal((12, 12))
    Bm = ss._scatter_dense_block(part, X)
    assert Bm.format == "csc"
    ref = np.zeros((30, 30), complex)
    ref[np.ix_(part.f_eq, part.f_eq)] = X
    np.testing.assert_array_equal(Bm.toarray(), ref)
    C, X2, p3 = _schur_system(nint=10, nn=12, seed=11)
    fac = ss.factorize(C, X2, p3, rcond_min=2.0)          # force the full sparse LU path
    assert fac.path == "full"
    A = C.toarray()
    A[np.ix_(p3.f_eq, p3.f_eq)] += X2
    b = np.arange(p3.neq) + 1j
    np.testing.assert_allclose(fac.solve(b[p3.f_eq], b[p3.n_eq])[:, 0], np.linalg.solve(A, b), rtol=1e-10)


# ======================================================================================
# The resource check bounds the measured peak of complete ANALYS runs
# ======================================================================================
FS = B.FrequencySet.fourier(0.01, 1024, [20, 60])


@pytest.fixture(scope="module")
def mat_dir(tmp_path_factory):
    """An 81-node rigid surface mat (one dense matrix = 0.9 MB; small files: the test disk may be tight)."""
    wd = tmp_path_factory.mktemp("analys_mem")
    site = B.layered_site([(2.0, 150.0, 300.0, 1.9, 0.05)], (400.0, 800.0, 2.1, 0.03))
    mdl = B.surface_rigid_mat(site, half_width=10.0, ndiv=8, mass=10.0)
    B.run_soil(wd, "m", site, FS, layer=0, rad=mdl.rad)
    B.run_house(wd, "m", mdl)
    return wd


def _measured_run(wd, **params):
    B.write_deck(wd, "m", B.analys_deck(FS, model="m", **params))
    rc, peak = _peak_increment(B.run, "ANALYS", wd, "m", False)
    text = B.listing(wd, "m", "ANALYS")
    est = float(re.search(r"estimated peak memory\s+:\s+([\d.]+) MB", text).group(1)) * 2 ** 20
    return rc, peak, est, text


@pytest.mark.parametrize("small_blocks", [True, False])
def test_estimate_bounds_the_measured_peak(mat_dir, monkeypatch, small_blocks):
    """D-ANL-10: the printed estimate must not be below the traced peak of the run -- initiation
    with and without restart files, New Structure, New Seismic Environment and the global impedance
    (FILE11 holds X_ff of every frequency).  With small work blocks the dense matrices dominate (the
    case of large models); with the default blocks the bounded work arrays do."""
    if small_blocks:
        monkeypatch.setattr(ss, "SCHUR_BLOCK_BYTES", 2 ** 16)
        monkeypatch.setattr(ss, "DENSE_BLOCK_ELEMS", 2 ** 10)
        monkeypatch.setattr(analys, "FLEX_BLOCK_ELEMS", 2 ** 10)
    D = (3 * 81) ** 2 * 16
    for params in (dict(), dict(save=1), dict(mode=1), dict(mode=1, save=1), dict(mode=2), dict(impe=2)):
        rc, peak, est, text = _measured_run(mat_dir, **params)
        assert rc == 0, text[-2000:]
        assert peak <= est, f"{params}: traced peak {peak / 2 ** 20:.1f} MB > estimate {est / 2 ** 20:.1f} MB"
        if small_blocks and params in (dict(), dict(mode=1)):
            # two dense matrices per frequency + the coordinate tables of F_ff + inputs/outputs of this
            # small model (observed about 3.0 D; the code reviewed held about 9 D per frequency)
            assert peak < 3.5 * D, f"{params}: traced peak {peak / D:.2f} dense matrices"


def test_fallback_path_memory_is_reported(mat_dir, monkeypatch):
    """The full sparse LU fallback needs more memory than the Schur path: it is listed, and EDU-20
    warns (without refusing the run) when only the fallback would exceed the memory limit."""
    md_est = {}

    def capture(*args, **kwargs):
        est = real(*args, **kwargs)
        md_est["est"] = est
        return est

    real = analys.memory_estimate
    monkeypatch.setattr(analys, "memory_estimate", capture)
    monkeypatch.setattr(ss, "SCHUR_BLOCK_BYTES", 2 ** 16)       # dense-dominated, as for large models
    monkeypatch.setattr(analys, "FLEX_BLOCK_ELEMS", 2 ** 10)
    monkeypatch.setattr(ss, "DENSE_BLOCK_ELEMS", 2 ** 10)
    rc, _, est, text = _measured_run(mat_dir)
    assert rc == 0 and "a frequency on the full sparse LU fallback path needs about" in text
    e = md_est["est"]
    need_fb = e.fallback_total
    assert need_fb > e.total
    ram = int((e.total + need_fb) / 2 / analys.MEMLIMIT_FRACTION)        # between the two requirements
    assert e.total < analys.MEMLIMIT_FRACTION * ram < need_fb
    monkeypatch.setattr(analys, "physical_memory", lambda: ram)
    rc, _, _, text = _measured_run(mat_dir)
    assert rc == 0 and "EDU-20 (G-22): a frequency that needs the full sparse LU fallback" in text
    monkeypatch.setattr(analys, "physical_memory", lambda: int(e.total / analys.MEMLIMIT_FRACTION * 0.9))
    rc, _, _, text = _measured_run(mat_dir)
    assert rc == 1 and "EDU-20 (G-22): the estimated memory" in text
