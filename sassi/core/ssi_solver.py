"""Per-frequency flexible-volume SSI solver: impedance, assembly, Schur-complement solve and the
storable factorisation used by the ANALYS restarts (requirements 4.6 items 3-6 and 8, R1 4.2-4.6,
decisions D-ANL-01/03/04/05).

The flexible-volume equation (for the structural engineer)
---------------------------------------------------------
At one circular frequency w the soil-structure system is (R1 4.4, general assembly D-ANL-01)::

    C(w) U = b,      C = (K*_s - w^2 M_s) - (K*_e - w^2 M_e) + A_f^T X_ff A_f

* ``K*_s, M_s``: structure (with any near-field soil modelled as structure), HOUSE COOSK/COOSM;
* ``K*_e, M_e``: the *excavated* soil -- the soil volume the basement displaces, with the free-field
  layer properties; it is subtracted because the soil impedance ``X_ff`` below already contains it;
* ``X_ff = F_ff^-1``: the dynamic stiffness of the free-field site at the interaction DOFs (the
  translations of the interaction nodes), from the point-load flexibility ``F_ff`` (POINT, FILE3);
* ``A_f`` gathers the interaction DOFs out of the global equations.

The load is ``b = A_f^T X_ff U'_f`` for seismic input (``U'_f`` = free-field motion at the
interaction nodes, SITE FILE1) or the nodal load vector of FORCE (FILE9) for foundation vibration.
Incoherent, wave-passage and multiple-excitation input multiply either the free-field load
(FFL) or the free-field motion (FFM) by the HOUSE factors of FILE77
(:func:`incoherent_seismic_load`).
``U`` are *total* displacement amplitudes per unit control motion (or per unit load factor).

Why a Schur complement (D-ANL-03)
---------------------------------
``X_ff`` is a *full* complex matrix on the interaction DOFs, the rest of ``C`` is sparse.  The
equations are therefore split into the non-interaction set ``n`` and the interaction set ``f``::

    [ C_nn  C_nf      ] [U_n]   [b_n]
    [ C_fn  C_ff + X  ] [U_f] = [b_f]

``C_nn`` is factorised by a sparse LU (SuperLU), the dense Schur complement
``S = C_ff - C_fn C_nn^-1 C_nf`` is formed and ``(S + X)`` is factorised by a dense LU::

    U_f = (S + X)^-1 (b_f - C_fn C_nn^-1 b_n),       U_n = C_nn^-1 (b_n - C_nf U_f)

Every right-hand side (all simultaneous load cases, restarts) reuses both factorisations.  All
matrices are complex *symmetric* (not Hermitian): LU factorisations only, never Cholesky.

If ``C_nn`` alone is singular or ill-conditioned (estimated reciprocal condition number below
:data:`RCOND_MIN`, e.g. an undamped structure with its interaction DOFs clamped at one of its
fixed-base frequencies) the full system with the dense ``X`` block inserted is factorised by a
sparse LU instead (path ``'full'``).  Models without non-interaction equations (pure flexible
volume with a structure made only of interaction nodes) use the dense path ``'dense'``.

Memory (D-ANL-10, G-22)
-----------------------
The dense ``3 nInt x 3 nInt`` matrices dominate the memory of large models (576 MB each for 2,000
interaction nodes).  The helpers below therefore work *in place* and with bounded work arrays:
``F_ff`` is symmetrised in place and overwritten by its inverse ``X_ff``
(:func:`impedance_matrix` with ``overwrite=True``: LAPACK getrf + getri on the same buffer), the
Schur complement ``S`` is built in Fortran order, ``X_ff`` is added to it in place and ``S + X_ff``
is overwritten by its LU factors; norms, finiteness checks and the symmetrisation run over blocks of
at most :data:`DENSE_BLOCK_ELEMS` entries and the Schur complement is formed from column blocks of
at most :data:`SCHUR_BLOCK_BYTES`.  One frequency of the Schur path thus holds two dense matrices
(``X_ff`` and the LU of ``S + X_ff``) plus bounded work arrays (see :func:`schur_work_bytes`).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Callable, Dict, Optional, Tuple

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .tlm import quiet_fpe

#: reciprocal condition number below which C_nn is treated as singular (requirements 4.6 item 6)
RCOND_MIN = 1e-12
#: memory budget of one block of C_nn^-1 C_nf columns while forming the Schur complement
SCHUR_BLOCK_BYTES = 64 * 2 ** 20
#: number of block-sized work arrays alive at once while forming the Schur complement (the dense
#: block of C_nf, its SuperLU copy, Z = C_nn^-1 C_nf, C_fn Z and the update of S; D-ANL-10)
SCHUR_WORK_ARRAYS = 5
#: entries per block of the in-place dense helpers (norm, finiteness, symmetrisation): 32 MB complex
DENSE_BLOCK_ELEMS = 2 ** 21
#: path codes stored in the restart records (COOTKqqq)
PATHS = ("schur", "dense", "full")


class SolverError(RuntimeError):
    """A singular or non-finite system that cannot be solved at one frequency."""


# ---------------------------------------------------------------------------------------
# Dense helpers
# ---------------------------------------------------------------------------------------
@quiet_fpe
def abs_column_sums(A: np.ndarray) -> np.ndarray:
    """Column sums of ``|A|`` of a dense matrix, computed over blocks of at most
    :data:`DENSE_BLOCK_ELEMS` entries along the contiguous axis (no full ``|A|`` temporary).  A
    non-finite entry gives a non-finite sum."""
    A = np.asarray(A)
    n, m = A.shape
    if A.flags.f_contiguous and not A.flags.c_contiguous:
        out = np.empty(m)
        step = max(1, DENSE_BLOCK_ELEMS // max(n, 1))
        for j0 in range(0, m, step):
            out[j0:j0 + step] = np.abs(A[:, j0:j0 + step]).sum(axis=0)
        return out
    out = np.zeros(m)
    step = max(1, DENSE_BLOCK_ELEMS // max(m, 1))
    for i0 in range(0, n, step):
        out += np.abs(A[i0:i0 + step]).sum(axis=0)
    return out


def _norm1(A) -> float:
    """1-norm (maximum absolute column sum) of a dense or sparse matrix."""
    if sp.issparse(A):
        if A.shape[0] == 0:
            return 0.0
        return float(abs(A).sum(axis=0).max())
    A = np.asarray(A)
    return float(abs_column_sums(A).max()) if A.size else 0.0


def symmetrize_inplace(A: np.ndarray, block: Optional[int] = None) -> np.ndarray:
    """``A <- (A + A^T) / 2`` in place, tile by tile (D-ANL-02, D-ANL-10).

    Gives bit-for-bit the values of ``0.5 * (A + A.T)`` (addition and multiplication are
    commutative in IEEE arithmetic) without the two full temporaries of that expression."""
    n = A.shape[0]
    if A.ndim != 2 or A.shape[1] != n:
        raise ValueError(f"symmetrize_inplace needs a square matrix, got shape {A.shape}")
    b = int(block) if block else max(1, int(np.sqrt(DENSE_BLOCK_ELEMS)))
    for i0 in range(0, n, b):
        i1 = min(n, i0 + b)
        for j0 in range(i0, n, b):
            j1 = min(n, j0 + b)
            T = A[i0:i1, j0:j1] + A[j0:j1, i0:i1].T          # new array: old values of both tiles
            T *= 0.5
            A[i0:i1, j0:j1] = T
            A[j0:j1, i0:i1] = T.T
    return A


@quiet_fpe
def dense_lu(A: np.ndarray, what: str = "matrix", overwrite: bool = False) -> Tuple[np.ndarray, np.ndarray, float]:
    """Dense LU with partial pivoting (LAPACK getrf) and the LAPACK 1-norm rcond estimate (gecon).

    ``overwrite=True`` factorises a complex Fortran-ordered ``A`` in place (its content is lost; a
    C-ordered array is copied as LAPACK needs column-major storage).  Raises :class:`SolverError`
    for a singular or non-finite matrix."""
    A = np.asarray(A, dtype=complex)
    if A.size == 0:
        return np.zeros((0, 0), complex), np.zeros(0, np.int32), 1.0
    colsum = abs_column_sums(A)
    if not np.all(np.isfinite(colsum)):
        raise SolverError(f"{what} contains non-finite values")
    anorm = float(colsum.max())
    del colsum
    with warnings.catch_warnings():
        warnings.simplefilter("error", sla.LinAlgWarning)
        try:
            lu, piv = sla.lu_factor(A, overwrite_a=bool(overwrite and A.flags.f_contiguous), check_finite=False)
        except (sla.LinAlgWarning, sla.LinAlgError, ValueError) as exc:
            raise SolverError(f"{what} is singular ({exc})") from None
    gecon, = sla.get_lapack_funcs(("gecon",), (lu,))
    rc, info = gecon(lu, anorm, norm="1")
    rc = float(rc) if info == 0 else float("nan")
    if not rc > 0.0:
        raise SolverError(f"{what} is singular (rcond = {rc:.3g})")
    return lu, piv, rc


@quiet_fpe
def impedance_matrix(F: np.ndarray, overwrite: bool = False) -> Tuple[np.ndarray, float]:
    """``X_ff = F_ff^-1`` by a dense LU factorisation and inversion (LAPACK getrf + getri;
    D-ANL-04, R1 4.2).

    ``F_ff`` is complex symmetric (D-ANL-02); the inverse is symmetrised ``(X + X^T)/2`` to remove
    the round-off asymmetry of the inversion.  ``overwrite=True`` inverts ``F`` in its own buffer
    (no other dense matrix is allocated, D-ANL-10): a C-ordered array is factorised through its
    transpose, which is Fortran-ordered, because ``inv(F^T)^T = inv(F)``.  Returns
    ``(X, rcond(F))`` with the LAPACK 1-norm estimate of ``rcond(F)``."""
    F = np.asarray(F, dtype=complex)
    n = F.shape[0]
    if n == 0:
        return np.zeros((0, 0), complex), 1.0
    A = F if overwrite else F.copy(order="F")
    if A.flags.f_contiguous:
        V, trans = A, False
    elif A.flags.c_contiguous:
        V, trans = A.T, True                    # Fortran view: LU of F^T in the same buffer
    else:
        A = V = np.asfortranarray(A)
        trans = False
    colsum = abs_column_sums(A)
    if not np.all(np.isfinite(colsum)):
        raise SolverError("the soil flexibility matrix F_ff contains non-finite values")
    anorm = float(colsum.max())                 # ||F||_1 = ||F^T||_inf
    del colsum
    with warnings.catch_warnings():
        warnings.simplefilter("error", sla.LinAlgWarning)
        try:
            lu, piv = sla.lu_factor(V, overwrite_a=True, check_finite=False)
        except (sla.LinAlgWarning, sla.LinAlgError, ValueError) as exc:
            raise SolverError(f"the soil flexibility matrix F_ff is singular ({exc})") from None
    gecon, getri, getri_lwork = sla.get_lapack_funcs(("gecon", "getri", "getri_lwork"), (lu,))
    # rcond_1(F) = rcond_inf(F^T): the infinity-norm estimate of the transposed factorisation
    rc, info = gecon(lu, anorm, norm="I" if trans else "1")
    rc = float(rc) if info == 0 else float("nan")
    if not rc > 0.0:
        raise SolverError(f"the soil flexibility matrix F_ff is singular (rcond = {rc:.3g})")
    lw = getri_lwork(n)
    lwork = max(n, int(np.real(lw[0] if isinstance(lw, tuple) else lw)))
    inv, info = getri(lu, piv, lwork=lwork, overwrite_lu=1)
    if info != 0:
        raise SolverError(f"the soil flexibility matrix F_ff is singular (getri info = {info})")
    X = inv.T if trans else inv                 # same buffer as A, in A's memory order
    del lu, inv
    return symmetrize_inplace(X), rc


def dynamic_matrix(Ks, Ms, Ke, Me, omega: float) -> sp.csr_matrix:
    """``C = (K*_s - w^2 M_s) - (K*_e - w^2 M_e)`` on the FILE4 equations (requirements 4.6 item 4)."""
    w2 = float(omega) ** 2
    C = (sp.csr_matrix(Ks, dtype=complex) - sp.csr_matrix(Ke, dtype=complex)) \
        - w2 * (sp.csr_matrix(Ms, dtype=complex) - sp.csr_matrix(Me, dtype=complex))
    C = sp.csr_matrix(C)
    C.sum_duplicates()
    C.sort_indices()
    return C


# ---------------------------------------------------------------------------------------
# Equation partition
# ---------------------------------------------------------------------------------------
@dataclass
class Partition:
    """Interaction (``f``) and non-interaction (``n``) equations.

    ``f_eq`` are the equation numbers of the interaction translations in the order of the
    impedance matrix (node-major x, y, z; translations that are fixed are dropped), ``x_pos``
    their positions in the full ``3 nInt`` impedance matrix, ``n_eq`` the other equations."""

    neq: int
    f_eq: np.ndarray
    x_pos: np.ndarray
    n_eq: np.ndarray

    @classmethod
    def from_int_eq(cls, int_eq, neq: int) -> "Partition":
        flat = np.asarray(int_eq, dtype=np.int64).reshape(-1)
        x_pos = np.flatnonzero(flat >= 0).astype(np.int64)
        f_eq = flat[x_pos]
        if np.unique(f_eq).size != f_eq.size:
            raise ValueError("an equation appears twice among the interaction DOFs")
        if f_eq.size and (f_eq.min() < 0 or f_eq.max() >= neq):
            raise ValueError("interaction equation numbers outside the FILE4 numbering")
        mask = np.ones(int(neq), bool)
        mask[f_eq] = False
        return cls(int(neq), f_eq, x_pos, np.flatnonzero(mask).astype(np.int64))

    @property
    def nf(self) -> int:
        return int(self.f_eq.size)

    @property
    def nn(self) -> int:
        return int(self.n_eq.size)

    def scatter(self, U_f: np.ndarray, U_n: np.ndarray) -> np.ndarray:
        """Assemble ``U`` (neq, m) from its interaction and non-interaction parts."""
        m = U_f.shape[1] if U_f.ndim == 2 else (U_n.shape[1] if U_n.ndim == 2 else 1)
        U = np.zeros((self.neq, m), complex)
        U[self.f_eq] = np.asarray(U_f).reshape(self.nf, m)
        U[self.n_eq] = np.asarray(U_n).reshape(self.nn, m)
        return U


# ---------------------------------------------------------------------------------------
# Sparse LU that can be stored and re-used (COOTKqqq, D-ANL-05)
# ---------------------------------------------------------------------------------------
class SparseFactor:
    """LU factors of a sparse complex matrix in the SuperLU convention ``Pr A Pc = L U``.

    Built by :meth:`factor` (SuperLU, scipy ``splu``) or restored from a restart record by
    :meth:`from_arrays`; a restored factor solves by two sparse triangular solves::

        A x = b   <=>   L U (Pc^T x) = Pr b   ->   y[perm_r] = b,  z = U^-1 L^-1 y,  x = z[perm_c]
    """

    def __init__(self, n: int, lu=None, L=None, U=None, perm_r=None, perm_c=None):
        self.n = int(n)
        self._lu = lu
        self._L = L
        self._U = U
        self.perm_r = perm_r
        self.perm_c = perm_c

    @classmethod
    def factor(cls, A, what: str = "matrix") -> "SparseFactor":
        A = sp.csc_matrix(A, dtype=complex)
        if A.shape[0] == 0:
            return cls(0)
        if not np.all(np.isfinite(A.data)):
            raise SolverError(f"{what} contains non-finite values")
        try:
            lu = spla.splu(A)
        except RuntimeError as exc:          # "Factor is exactly singular"
            raise SolverError(f"{what} is singular ({exc})") from None
        return cls(A.shape[0], lu=lu, perm_r=np.asarray(lu.perm_r, np.int64),
                   perm_c=np.asarray(lu.perm_c, np.int64))

    def solve(self, B: np.ndarray) -> np.ndarray:
        B = np.asarray(B, dtype=complex)
        if self.n == 0:
            return np.zeros_like(B)
        if self._lu is not None:
            return self._lu.solve(B)
        y = np.empty_like(B)
        y[self.perm_r] = B
        z = spla.spsolve_triangular(self._L, y, lower=True, unit_diagonal=True)
        z = spla.spsolve_triangular(self._U, z, lower=False)
        return np.asarray(z)[self.perm_c]

    def solve_h(self, B: np.ndarray) -> np.ndarray:
        """``A^H x = b`` (only for a fresh SuperLU factor; used by the condition estimate)."""
        return self._lu.solve(np.asarray(B, dtype=complex), trans="H")

    def arrays(self, prefix: str) -> Dict[str, object]:
        """Arrays of a restart record: L, U (CSR) and the row/column permutations."""
        if self.n == 0:
            return {f"{prefix}n": np.int64(0)}
        L = self._L if self._L is not None else sp.csr_matrix(self._lu.L)
        U = self._U if self._U is not None else sp.csr_matrix(self._lu.U)
        return {f"{prefix}n": np.int64(self.n), f"{prefix}L": sp.csr_matrix(L), f"{prefix}U": sp.csr_matrix(U),
                f"{prefix}perm_r": self.perm_r, f"{prefix}perm_c": self.perm_c}

    @classmethod
    def from_arrays(cls, c, prefix: str) -> "SparseFactor":
        n = int(c[f"{prefix}n"])
        if n == 0:
            return cls(0)
        return cls(n, L=c.sparse(f"{prefix}L").tocsr(), U=c.sparse(f"{prefix}U").tocsr(),
                   perm_r=np.asarray(c[f"{prefix}perm_r"], np.int64), perm_c=np.asarray(c[f"{prefix}perm_c"], np.int64))


@quiet_fpe
def inverse_norm1_estimate(solve: Callable[[np.ndarray], np.ndarray], solve_h: Callable[[np.ndarray], np.ndarray],
                           n: int, itmax: int = 5) -> float:
    """Estimate ``||A^-1||_1`` from solves with A and A^H (Hager 1984 / Higham 1988, the algorithm of
    LAPACK xLACON).  Deterministic (fixed start vectors), a few solves only."""
    if n == 0:
        return 0.0
    x = np.full(n, 1.0 / n, dtype=complex)
    est = 0.0
    for it in range(itmax):
        y = solve(x)
        gamma = float(np.abs(y).sum())
        if it > 0 and gamma <= est:
            break
        est = gamma
        ay = np.abs(y)
        xi = np.where(ay > 0, y / np.where(ay > 0, ay, 1.0), 1.0)
        z = solve_h(xi)
        j = int(np.argmax(np.abs(z)))
        if it > 0 and np.abs(z[j]) <= np.real(np.vdot(z, x)):
            break
        x = np.zeros(n, dtype=complex)
        x[j] = 1.0
    # Higham's alternative vector guards against the rare failures of the power iteration
    i = np.arange(n)
    b = ((-1.0) ** i) * (1.0 + i / max(n - 1, 1))
    alt = 2.0 * float(np.abs(solve(b.astype(complex))).sum()) / (3.0 * n)
    return max(est, alt)


def sparse_rcond(A, fac: SparseFactor) -> float:
    """Estimated reciprocal 1-norm condition number ``1/(||A||_1 ||A^-1||_1)`` of a factorised
    sparse matrix (decides the Schur / full-LU path, requirements 4.6 item 6)."""
    n = A.shape[0]
    if n == 0:
        return 1.0
    inv = inverse_norm1_estimate(fac.solve, fac.solve_h, n)
    nrm = _norm1(A)
    if not np.isfinite(inv) or inv <= 0 or nrm <= 0:
        return 0.0
    return 1.0 / (nrm * inv)


# ---------------------------------------------------------------------------------------
# The factorised SSI system of one frequency
# ---------------------------------------------------------------------------------------
@dataclass
class SSIFactor:
    """Factorised SSI system at one frequency (path ``'schur'``, ``'dense'`` or ``'full'``).

    It is everything a "New Seismic Environment" / "New Dynamic Loading" restart (ANALYS solver
    Mode 3) needs besides the impedance: the sparse factors of ``C_nn``, the coupling blocks
    ``C_nf``/``C_fn`` and the dense LU of ``S + X`` -- or the sparse factors of the full matrix."""

    path: str
    part: Partition
    nn_fac: Optional[SparseFactor] = None
    Cnf: Optional[sp.csr_matrix] = None
    Cfn: Optional[sp.csr_matrix] = None
    lu: Optional[np.ndarray] = None
    piv: Optional[np.ndarray] = None
    full_fac: Optional[SparseFactor] = None
    rcond_nn: float = float("nan")
    rcond_ff: float = float("nan")
    note: str = ""

    @quiet_fpe
    def solve(self, b_f: np.ndarray, b_n: Optional[np.ndarray] = None) -> np.ndarray:
        """Solve for ``U`` (neq, m) given the interaction loads ``b_f`` (nf, m) and the
        non-interaction loads ``b_n`` (nn, m; None = zero, the seismic case).

        The right-hand sides are solved one at a time with the shared factorisations, so each
        load case gets bit-for-bit the result of a single-case run (VP-23)."""
        p = self.part
        b_f = np.asarray(b_f, dtype=complex)
        if p.nf == 0 and b_n is not None:
            m = int(np.asarray(b_n).reshape(p.nn, -1).shape[1])
        else:
            m = int(b_f.reshape(p.nf, -1).shape[1]) if p.nf else 1
        b_f = b_f.reshape(p.nf, m)
        b_n = np.zeros((p.nn, m), complex) if b_n is None else np.asarray(b_n, dtype=complex).reshape(p.nn, m)
        U = np.zeros((p.neq, m), complex)
        for j in range(m):
            U[:, j] = self._solve1(b_f[:, j:j + 1], b_n[:, j:j + 1])[:, 0]
        return U

    def _solve1(self, b_f: np.ndarray, b_n: np.ndarray) -> np.ndarray:
        p = self.part
        if self.path == "full":
            b = np.zeros((p.neq, 1), complex)
            b[p.f_eq] = b_f
            b[p.n_eq] = b_n
            return np.asarray(self.full_fac.solve(b)).reshape(p.neq, 1)
        if self.path == "dense":
            U_f = sla.lu_solve((self.lu, self.piv), b_f, check_finite=False)
            return p.scatter(U_f, np.zeros((0, 1), complex))
        # Schur complement path (R1 4.6)
        rhs = b_f - self.Cfn @ np.asarray(self.nn_fac.solve(b_n)).reshape(p.nn, 1) if np.any(b_n) else b_f
        U_f = sla.lu_solve((self.lu, self.piv), rhs, check_finite=False)
        U_n = self.nn_fac.solve(b_n - self.Cnf @ U_f)
        return p.scatter(U_f, np.asarray(U_n).reshape(p.nn, 1))

    # ---- restart record (COOTKqqq) ------------------------------------------------------
    def arrays(self) -> Dict[str, object]:
        p = self.part
        out: Dict[str, object] = {"f_eq": p.f_eq, "x_pos": p.x_pos, "n_eq": p.n_eq,
                                  "path": np.int64(PATHS.index(self.path))}
        if self.path == "full":
            out.update(self.full_fac.arrays("full_"))
        else:
            out["lu"] = self.lu
            out["piv"] = np.asarray(self.piv, dtype=np.int64)
        if self.path == "schur":
            out.update(self.nn_fac.arrays("nn_"))
            out["Cnf"] = self.Cnf
            out["Cfn"] = self.Cfn
        return out

    @classmethod
    def from_arrays(cls, c, neq: int) -> "SSIFactor":
        part = Partition(int(neq), np.asarray(c["f_eq"], np.int64), np.asarray(c["x_pos"], np.int64),
                         np.asarray(c["n_eq"], np.int64))
        path = PATHS[int(c["path"])]
        f = cls(path=path, part=part, rcond_nn=float(c.meta.get("rcond_nn", float("nan"))),
                rcond_ff=float(c.meta.get("rcond_ff", float("nan"))))
        if path == "full":
            f.full_fac = SparseFactor.from_arrays(c, "full_")
            return f
        f.lu = np.asarray(c["lu"], dtype=complex)
        f.piv = np.asarray(c["piv"], dtype=np.int32)
        if path == "schur":
            f.nn_fac = SparseFactor.from_arrays(c, "nn_")
            f.Cnf = c.sparse("Cnf").tocsr()
            f.Cfn = c.sparse("Cfn").tocsr()
        return f


def _scatter_dense_block(part: Partition, Xr: np.ndarray) -> sp.csc_matrix:
    """The dense interaction block ``A_f^T X A_f`` as a sparse (CSC) matrix on the full numbering.

    Built directly in compressed form (sorted equation numbers, 32-bit indices where possible):
    one gathered copy of ``Xr`` plus a quarter of it for the row indices (D-ANL-10)."""
    nf, neq = part.nf, part.neq
    if nf == 0:
        return sp.csc_matrix((neq, neq), dtype=complex)
    order = np.argsort(part.f_eq, kind="stable")
    eq = part.f_eq[order]
    itype = np.int32 if max(neq, nf * nf) < 2 ** 31 - 1 else np.int64
    Xr = np.asarray(Xr, complex)
    data = np.empty(nf * nf, complex)
    for j0 in range(0, nf, max(1, DENSE_BLOCK_ELEMS // nf)):       # column j of the block -> CSC column eq[j]
        cj = order[j0:j0 + max(1, DENSE_BLOCK_ELEMS // nf)]
        data[j0 * nf:(j0 + cj.size) * nf] = Xr[np.ix_(order, cj)].T.reshape(-1)
    indices = np.tile(eq.astype(itype), nf)
    indptr = np.zeros(neq + 1, itype)
    indptr[eq + 1] = nf
    np.cumsum(indptr, out=indptr)
    return sp.csc_matrix((data, indices, indptr), shape=(neq, neq))


@quiet_fpe
def factorize(C: sp.csr_matrix, Xr: np.ndarray, part: Partition, rcond_min: float = RCOND_MIN,
              block_bytes: Optional[int] = None) -> SSIFactor:
    """Factorise ``C + A_f^T Xr A_f`` (requirements 4.6 item 6, D-ANL-03).

    ``C`` is the dynamic matrix (:func:`dynamic_matrix`), ``Xr`` the impedance restricted to the
    active interaction DOFs (``X[x_pos][:, x_pos]``, see :func:`restrict`; it is not modified).
    ``block_bytes`` (default :data:`SCHUR_BLOCK_BYTES`) bounds one column block of the Schur
    complement work.  Raises :class:`SolverError` when the system is singular."""
    block_bytes = int(block_bytes or SCHUR_BLOCK_BYTES)
    C = sp.csr_matrix(C, dtype=complex)
    nf, nn = part.nf, part.nn
    Xr = np.asarray(Xr, dtype=complex).reshape(nf, nf)
    if nf == 0:                                 # no soil: plain sparse system (vibration, fixed base)
        fac = SparseFactor.factor(C, "the dynamic stiffness matrix")
        return SSIFactor(path="full", part=part, full_fac=fac, note="no interaction DOFs")
    if nn == 0:                                 # every equation is an interaction DOF
        A = C[part.f_eq][:, part.f_eq].toarray(order="F")
        A += Xr
        lu, piv, rc = dense_lu(A, "the SSI matrix C_ff + X_ff", overwrite=True)
        return SSIFactor(path="dense", part=part, lu=lu, piv=piv, rcond_ff=rc)
    Cnn = C[part.n_eq][:, part.n_eq].tocsc()
    fac: Optional[SparseFactor] = None
    rc_nn = 0.0
    why = ""
    try:
        fac = SparseFactor.factor(Cnn, "C_nn")
        rc_nn = sparse_rcond(Cnn, fac)
        if not rc_nn >= rcond_min:
            why = f"C_nn ill-conditioned (estimated rcond {rc_nn:.2e} < {rcond_min:.0e})"
    except SolverError as exc:
        why = str(exc)
    if why:
        del fac, Cnn
        B = _scatter_dense_block(part, Xr)
        full = C.tocsc() + B                         # CSC throughout: SuperLU takes it without a copy
        del B
        ffac = SparseFactor.factor(full, "the SSI system matrix")
        del full
        return SSIFactor(path="full", part=part, full_fac=ffac, rcond_nn=rc_nn,
                         note=f"{why}: full sparse LU with the dense X_ff block")
    Cnf = C[part.n_eq][:, part.f_eq].tocsc()
    Cfn = C[part.f_eq][:, part.n_eq].tocsr()
    S = C[part.f_eq][:, part.f_eq].toarray(order="F")   # C_ff; Fortran order: LU in place below
    cols = np.flatnonzero(np.diff(Cnf.indptr))       # interaction DOFs coupled to the n set
    step = max(1, int(block_bytes // (16 * max(nn, nf, 1))))
    for i0 in range(0, cols.size, step):
        cb = cols[i0:i0 + step]
        Z = fac.solve(Cnf[:, cb].toarray())          # C_nn^-1 C_nf (block of columns)
        W = Cfn @ Z
        del Z
        if cb[-1] - cb[0] + 1 == cb.size:            # contiguous columns: update the view in place
            S[:, cb[0]:cb[-1] + 1] -= W
        else:
            S[:, cb] -= W
        del W
    S += Xr                                          # S + X_ff, in place
    lu, piv, rc = dense_lu(S, "the condensed SSI matrix S + X_ff", overwrite=True)
    del S
    return SSIFactor(path="schur", part=part, nn_fac=fac, Cnf=Cnf.tocsr(), Cfn=Cfn, lu=lu, piv=piv,
                     rcond_nn=rc_nn, rcond_ff=rc)


def schur_work_bytes(nn: int, nf: int, block_bytes: Optional[int] = None) -> int:
    """Bound of the work arrays of :func:`factorize` besides ``X_ff`` and ``S``: at most
    :data:`SCHUR_WORK_ARRAYS` blocks of ``max(nn, nf) x step`` complex entries (D-ANL-10)."""
    if nn == 0 or nf == 0:
        return 0
    block_bytes = int(block_bytes or SCHUR_BLOCK_BYTES)
    step = max(1, int(block_bytes // (16 * max(nn, nf, 1))))
    return SCHUR_WORK_ARRAYS * 16 * max(nn, nf) * min(step, nf)


def restrict(X: np.ndarray, part: Partition) -> np.ndarray:
    """``X_ff`` restricted to the active interaction DOFs (``X[x_pos][:, x_pos]``); ``X`` itself
    (no copy) when no interaction translation is fixed (D-ANL-10)."""
    n3 = X.shape[0]
    if part.nf == n3 and np.array_equal(part.x_pos, np.arange(n3)):
        return X
    return X[np.ix_(part.x_pos, part.x_pos)]


# ---------------------------------------------------------------------------------------
# Loads and rigid-body impedance
# ---------------------------------------------------------------------------------------
@quiet_fpe
def seismic_load(X: np.ndarray, Up: np.ndarray, part: Partition) -> np.ndarray:
    """Interaction loads ``b_f = A_f^T X_ff U'_f`` (requirements 4.6 item 5).

    ``Up`` is the free field at the interaction nodes, (nInt, 3) or (nInt, 3, m) for m cases, or
    already flattened node-major (3 nInt,) / (3 nInt, m).  Returns (nf, m)."""
    Up = np.asarray(Up, dtype=complex)
    n3 = X.shape[0]
    Uf = Up.reshape(n3, -1)
    out = np.zeros((part.nf, Uf.shape[1]), complex)
    for j in range(Uf.shape[1]):                     # one case at a time: same rounding as a single run
        out[:, j] = (X @ Uf[:, j])[part.x_pos]
    return out


@quiet_fpe
def incoherent_seismic_load(X: np.ndarray, Up: np.ndarray, s: np.ndarray, part: Partition,
                            ffm: bool = False) -> np.ndarray:
    """Interaction loads of incoherent / wave-passage / multiple-excitation input (requirements 4.6
    item 5, spec 05c A.5.9; D-INC-12):

    * FFL (free-field load, default): ``b_f = A_f^T (s .* X_ff U'_f)`` (``.*`` element by element) -- the coherent free-field
      load ``X_ff U'_f`` of every interaction node is multiplied by its factor;
    * FFM (free-field motion): ``b_f = A_f^T X_ff (s .* U'_f)`` -- the free-field motion is multiplied
      by the factors before the coherent impedance acts on it.

    ``s`` (nInt,) or (nInt, m) holds one complex factor per interaction node (the factor of the
    control-motion direction multiplies the three components of the node); ``Up`` as in
    :func:`seismic_load`.  With ``s = 1`` both forms reproduce :func:`seismic_load` exactly (a product
    by 1 + 0i is exact in IEEE arithmetic).  Returns (nf, m)."""
    Up = np.asarray(Up, dtype=complex)
    n3 = X.shape[0]
    Uf = Up.reshape(n3, -1)
    m = Uf.shape[1]
    s = np.asarray(s, dtype=complex).reshape(n3 // 3, -1)
    if s.shape[1] not in (1, m):
        raise ValueError(f"incoherency factors for {s.shape[1]} cases, {m} load cases")
    out = np.zeros((part.nf, m), complex)
    for j in range(m):                               # one case at a time: same rounding as a single run
        s3 = np.repeat(s[:, j if s.shape[1] > 1 else 0], 3)
        if ffm:
            out[:, j] = (X @ (s3 * Uf[:, j]))[part.x_pos]
        else:
            out[:, j] = (s3 * (X @ Uf[:, j]))[part.x_pos]
    return out


def rigid_body_transform(xyz: np.ndarray, ref) -> np.ndarray:
    """``T`` (3n, 6): translations of n points for a rigid-body motion (u0, theta) about ``ref``,
    ``u_j = u0 + theta x d_j`` with ``d_j = r_j - ref`` (requirements 4.6 item 9)::

        T_j = [[1,0,0,  0, dz,-dy],
               [0,1,0,-dz,  0, dx],
               [0,0,1, dy,-dx,  0]]
    """
    d = np.asarray(xyz, float).reshape(-1, 3) - np.asarray(ref, float).reshape(1, 3)
    n = d.shape[0]
    T = np.zeros((n, 3, 6))
    T[:, 0, 0] = T[:, 1, 1] = T[:, 2, 2] = 1.0
    dx, dy, dz = d[:, 0], d[:, 1], d[:, 2]
    T[:, 0, 4], T[:, 0, 5] = dz, -dy
    T[:, 1, 3], T[:, 1, 5] = -dz, dx
    T[:, 2, 3], T[:, 2, 4] = dy, -dx
    return T.reshape(3 * n, 6)


@quiet_fpe
def global_impedance(X: np.ndarray, T: np.ndarray) -> np.ndarray:
    """"Unconstrained" foundation impedance ``K_G = T^T X_ff T`` (6x6, DOF order X, Y, Z, XX, YY,
    ZZ; requirements 4.6 item 9, D-ANL-07)."""
    return T.T @ np.asarray(X, complex) @ T
