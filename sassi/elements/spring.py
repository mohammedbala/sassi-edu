"""SPRING (GROUP type 7): six uncoupled global spring constants between nodes I and J.

Normative sources: requirements 1.5 and 4.1 (SPRING row), 4.0.2 (k* = k c(damp)), spec 08
sections 4.8 and 5.34 (SC table: kx ky kz kxx kyy kzz in global axes + one damping ratio; no mass;
zero length allowed), D-STR-05 (recovery F = k* (u_J - u_I) per global component).

For each global DOF d = UX..ROTZ the spring adds::

    K[d_I, d_I] += k_d*   K[d_J, d_J] += k_d*   K[d_I, d_J] -= k_d*   K[d_J, d_I] -= k_d*

with the hysteretic constant ``k_d* = k_d c(damp)``.  The spring has no orientation (the local
system and the node positions are irrelevant) and carries no mass.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from ..conventions import CMOD_SASSI, cfactor
from .base import ElementError, ElementSpec, register

NAME = "SPRING"
CODE = 7
DOFS = (1, 2, 3, 4, 5, 6)
NNODES = 2
COMPONENTS = ["FX", "FY", "FZ", "MXX", "MYY", "MZZ"]


def spring_constants(k: Sequence[float], damp: float = 0.0, form: int = CMOD_SASSI) -> np.ndarray:
    """Complex constants k* = k c(damp) (6,)."""
    kk = np.asarray(k, dtype=float).ravel()
    if kk.size != 6:
        raise ElementError("SPRING needs six constants (kx, ky, kz, kxx, kyy, kzz)")
    if np.any(kk < 0):
        raise ElementError("SPRING constants must be >= 0 (Errors 34-39)")
    try:
        c = cfactor(float(damp), form)
    except ValueError as exc:
        raise ElementError(str(exc)) from None
    return kk * c


def matrices(xyz=None, mat=None, k: Sequence[float] = None, damp: float = 0.0, form: int = CMOD_SASSI,
             **_ignored):
    """(K* (12,12) complex, M = 0 (12,12)) in global axes; ``xyz`` and ``mat`` are not used."""
    if k is None:
        raise ElementError("SPRING element needs its SC constants k=(kx, ky, kz, kxx, kyy, kzz)")
    ks = spring_constants(k, damp, form)
    K = np.zeros((12, 12), dtype=complex)
    d = np.arange(6)
    K[d, d] = ks
    K[d + 6, d + 6] = ks
    K[d, d + 6] = -ks
    K[d + 6, d] = -ks
    return K, np.zeros((12, 12))


def recovery(xyz=None, mat=None, k: Sequence[float] = None, damp: float = 0.0, form: int = CMOD_SASSI,
             **_ignored):
    """Spring forces (6, 12): F_d = k_d* (u_J,d - u_I,d), global components FX..MZZ (D-STR-05)."""
    ks = spring_constants(k, damp, form)
    S = np.zeros((6, 12), dtype=complex)
    d = np.arange(6)
    S[d, d] = -ks
    S[d, d + 6] = ks
    return S


SPEC = register(ElementSpec(code=CODE, name=NAME, nnodes=NNODES, dofs=DOFS, components=list(COMPONENTS),
                            matrices=matrices, recovery=recovery,
                            description="six uncoupled global springs with hysteretic damping, no mass"))
