"""GENERAL (GROUP type 9, "GM"): user 12x12 complex stiffness and mass between two nodes.

Normative sources: requirements 1.5 and 4.1 (GENERAL row: ``K* = K_R + i K_I``, mass divided by g
in weight units, 3 nodes -> local axes as BEAMS, ``T = blockdiag(Lam, Lam, Lam, Lam)``), spec 08
sections 4.9 and 5.26-5.29 (MXR/MXI/MXM upper-triangle rows; DOFs 1-6 node I, 7-12 node J;
global input with nodes I, J; local input with I, J, K), D-CNV-04 (GM damping as entered),
UT-15 (MOPT <matrix> = 1: mass entered as weight).

Each GENERAL element couples exactly two nodes; larger super-elements are assembled from several
GENERAL elements, one per node pair (the user splits shared diagonal blocks).
"""
from __future__ import annotations

from typing import Mapping, Optional, Sequence, Union

import numpy as np

from .base import ElementError, ElementSpec, blas_quiet, check_finite, frame_from_three_points, register

NAME = "GENERAL"
CODE = 9
DOFS = (1, 2, 3, 4, 5, 6)
NNODES = 2
COMPONENTS: list = []          # STRESS has no output components for GENERAL elements


def from_upper_rows(rows: Union[np.ndarray, Mapping[int, Sequence[float]], None]) -> np.ndarray:
    """Full symmetric 12x12 matrix from upper-triangle rows (MXR/MXI/MXM layout).

    ``rows`` is either a mapping ``{row (1..12): [t1, t2, ...]}`` with ``t1 = A[r, r]``,
    ``t2 = A[r, r+1]`` ... (at most 13 - r terms, missing trailing terms are 0), or a (12, 12) table
    whose row r holds those terms in its first 13 - r columns (HOUSE deck ``matrices`` table).
    """
    A = np.zeros((12, 12))
    if rows is None:
        return A
    if isinstance(rows, Mapping):
        items = [(int(r), list(v)) for r, v in rows.items()]
    else:
        tab = np.asarray(rows, dtype=float)
        if tab.shape != (12, 12):
            raise ElementError("upper-triangle table must be 12 x 12")
        items = [(r + 1, list(tab[r, :12 - r])) for r in range(12)]
        if any(np.any(tab[r, 12 - r:] != 0) for r in range(12)):
            raise ElementError("row r of an MX table may hold at most 13 - r terms")
    for r, terms in items:
        if not 1 <= r <= 12:
            raise ElementError(f"matrix row {r} outside 1..12")
        if len(terms) > 13 - r:
            raise ElementError(f"matrix row {r} has {len(terms)} terms (at most {13 - r})")
        A[r - 1, r - 1:r - 1 + len(terms)] = terms
    return np.triu(A) + np.triu(A, 1).T


def _full(Mx, name: str, upper_rows: bool = False) -> np.ndarray:
    if Mx is None:
        return np.zeros((12, 12))
    if isinstance(Mx, Mapping) or upper_rows:
        A = from_upper_rows(Mx)
    else:
        A = np.asarray(Mx, dtype=float)
    if A.shape != (12, 12):
        raise ElementError(f"GENERAL {name} must be 12 x 12")
    if not np.allclose(A, A.T, rtol=1e-12, atol=1e-12 * max(1.0, float(np.abs(A).max()))):
        raise ElementError(f"GENERAL {name} must be symmetric (enter the upper triangle with MX{name[-1]})")
    return 0.5 * (A + A.T)


def transformation(xyz) -> Optional[np.ndarray]:
    """T (12x12) for local input (3 nodes I, J, K) or None for global input (2 nodes)."""
    if xyz is None:
        return None
    X = np.asarray(xyz, dtype=float)
    if X.shape[0] == 2:
        return None
    if X.shape != (3, 3):
        raise ElementError("GENERAL element needs 2 (global) or 3 (local, I J K) nodes")
    Lam, _ = frame_from_three_points(X[0], X[1], X[2])
    return np.kron(np.eye(4), Lam)


def matrices(xyz=None, mat=None, KR=None, KI=None, MM=None, munits: int = 0, gravity: Optional[float] = None,
             upper_rows: bool = False, **_ignored):
    """(K* = K_R + i K_I, M) in global axes (12x12).

    ``KR``, ``KI``, ``MM``: full symmetric 12x12 arrays, upper-row mappings ``{row: terms}``, or -
    with ``upper_rows=True`` - (12, 12) tables in the MXR/MXI/MXM row layout of the HOUSE deck
    (see :func:`from_upper_rows`).  ``munits`` = MOPT <matrix>: 0 mass units, 1 weight units (the
    mass is divided by ``gravity``).  With 3 rows in ``xyz`` the matrices are in the local I-J-K
    axes and are transformed as ``T^T A T``.
    """
    K = _full(KR, "KR", upper_rows) + 1j * _full(KI, "KI", upper_rows)
    M = _full(MM, "MM", upper_rows)
    if int(munits) == 1:
        if gravity is None or gravity <= 0:
            raise ElementError("GENERAL mass in weight units needs gravity > 0")
        M = M / float(gravity)
    T = transformation(xyz)
    if T is not None:
        with blas_quiet():
            K = T.T @ K @ T
            M = T.T @ M @ T
    return check_finite(0.5 * (K + K.T), "GENERAL K"), 0.5 * (M + M.T)


def recovery(xyz=None, mat=None, **_ignored):
    """GENERAL elements have no STRESS output components: an empty (0, 12) operator."""
    return np.zeros((0, 12), dtype=complex)


SPEC = register(ElementSpec(code=CODE, name=NAME, nnodes=NNODES, dofs=DOFS, components=list(COMPONENTS),
                            matrices=matrices, recovery=recovery,
                            description="user 12x12 complex stiffness and mass (2 nodes global, 3 nodes local)"))
