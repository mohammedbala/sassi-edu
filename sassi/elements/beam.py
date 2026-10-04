"""BEAMS (GROUP type 2): 2-node 3D Timoshenko frame element with an orientation node K.

Normative sources: requirements 1.5 and 4.1 (BEAMS row), D-ELM-04 (Timoshenko with shear areas,
0 = no shear deformation; end releases by static condensation, the condensed mass uses the same
transformation), D-ELM-05 (consistent mass, torsional inertia rho (I2+I3), no rotary inertia),
D-CNV-04 (E*, G* from M*, G*), D-STR-05 (end forces exerted on the element, local axes),
spec 08 sections 4.5, 5.16-5.17, 5.30 (local axes, R table, KI/KJ).

Local axes (spec 08 4.5): e1 = I->J, e2 = component of (K - I) normal to e1 (axis 2 toward K),
e3 = e1 x e2.  Local DOFs per node: u1 u2 u3 (translations) and t1 t2 t3 (rotations).

Section (R table): A (axial), As2, As3 (shear areas for local axes 2 and 3), J (torsion), I2, I3
(inertias about local axes 2 and 3).  Bending in the 1-2 plane (deflection u2, rotation t3) uses
E I3 and As2; bending in the 1-3 plane (deflection u3, rotation t2) uses E I2 and As3.  The shear
deformation parameter is ``phi = 12 E I / (G As L^2)`` (phi = 0 when As = 0, Euler-Bernoulli).

Stiffness for one bending plane (deflection v, rotation t = dv/dx, Przemieniecki 1968)::

               E I         [ 12      6L        -12     6L      ]
    k = ---------------  * [ 6L   (4+phi)L^2  -6L  (2-phi)L^2  ]
         (1 + phi) L^3     [-12     -6L        12     -6L      ]
                           [ 6L   (2-phi)L^2  -6L  (4+phi)L^2  ]

In the 1-3 plane the rotation about axis 2 is t2 = -du3/dx, which flips the sign of the
translation-rotation coupling terms.  The single element is exact for end loads (cantilever tip
deflection ``P L^3/(3EI) + P L/(G As)``).

Mass (D-ELM-05): consistent translational mass from the Timoshenko displacement interpolation
(reduces to the classical 156/22/54/13 Euler-Bernoulli matrix for phi = 0), axial
``rho A L/6 [2 1; 1 2]``, torsional ``rho (I2+I3) L/6 [2 1; 1 2]``; no rotary inertia.  The mass is
real, so it is built with the undamped phi0.

End releases (KI/KJ, 1 = released): the released local DOFs r are condensed statically,
``K_c = K_cc - K_cr K_rr^-1 K_rc`` (complex K*), and the mass is transformed with the real
undamped transformation ``u = T u_c``, ``T = [I; -K0_rr^-1 K0_rc]`` (D-ELM-04).  Releasing the same
component at both ends is CHECK Error 10 (mechanism for P1/M1; not allowed by the manual).
A zero inertia (I2 or I3 = 0 is accepted, only < 0 is Error 31/32) leaves its bending plane
without stiffness, so a release in that plane is redundant (the end force is already zero) and
``K_rr`` would be singular: such components are left out of the stiffness condensation, and the
mass transformation uses the limit I -> 0+ (the Euler-Bernoulli Hermite shape of that plane,
because phi -> 0), so the matrices are continuous in I.

Global matrices: ``K = T^T K_l T`` with ``T = blockdiag(Lam, Lam, Lam, Lam)`` (rows of Lam = e1 e2 e3).
Recovery (D-STR-05): local end forces ``f = K_l,c T u`` = forces exerted *on the element* by the
nodes, ordered FXI FYI FZI MXI MYI MZI FXJ ... MZJ (local 1/2/3 written as X/Y/Z).
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np

from .base import (ElementError, ElementSpec, MaterialProps, as_bool6, blas_quiet, check_finite,
                   frame_from_three_points, gauss_legendre, register)

NAME = "BEAMS"
CODE = 2
DOFS = (1, 2, 3, 4, 5, 6)
NNODES = 2
COMPONENTS = ["FXI", "FYI", "FZI", "MXI", "MYI", "MZI", "FXJ", "FYJ", "FZJ", "MXJ", "MYJ", "MZJ"]

_SECTION_KEYS = ("A", "As2", "As3", "J", "I2", "I3")


def parse_section(section) -> Dict[str, float]:
    """Section dict(A, As2, As3, J, I2, I3) from a dict or a 6-sequence in R-table order."""
    if section is None:
        raise ElementError("BEAMS element needs a section (R table: A, As2, As3, J, I2, I3)")
    if isinstance(section, dict):
        try:
            s = {k: float(section[k]) for k in _SECTION_KEYS}
        except KeyError as exc:
            raise ElementError(f"beam section is missing {exc}") from None
    else:
        vals = list(section)
        if len(vals) != 6:
            raise ElementError("beam section needs 6 values A, As2, As3, J, I2, I3")
        s = dict(zip(_SECTION_KEYS, map(float, vals)))
    if s["A"] <= 0:
        raise ElementError("Error 27: beam axial area must be > 0")
    if s["J"] <= 0:
        raise ElementError("Error 30: beam torsion constant must be > 0")
    if min(s["As2"], s["As3"], s["I2"], s["I3"]) < 0:
        raise ElementError("beam shear areas and inertias must be >= 0 (Errors 28, 29, 31, 32)")
    return s


def parse_release(code) -> Tuple[bool, ...]:
    """KI/KJ code -> 6 bools for P1 P2 P3 M1 M2 M3 (1 = released); accepts '000011' or sequences."""
    return as_bool6(code)


def _phi(E, G, I, As, L):
    if As <= 0.0:
        return 0.0 * E / G                                     # Euler-Bernoulli (keeps the dtype)
    return 12.0 * E * I / (G * As * L * L)


def _bending_k(EI, phi, L):
    """4x4 stiffness for (v1, t1, v2, t2) with t = dv/dx."""
    c = EI / ((1.0 + phi) * L ** 3)
    L2 = L * L
    return c * np.array([[12.0, 6 * L, -12.0, 6 * L],
                         [6 * L, (4.0 + phi) * L2, -6 * L, (2.0 - phi) * L2],
                         [-12.0, -6 * L, 12.0, -6 * L],
                         [6 * L, (2.0 - phi) * L2, -6 * L, (4.0 + phi) * L2]])


def timoshenko_shape(x: np.ndarray, L: float, phi: float) -> np.ndarray:
    """Transverse-displacement shape functions (nx, 4) of the exact Timoshenko element for
    DOFs (v1, t1, v2, t2), t = dv/dx at the nodes is the *bending* rotation (Przemieniecki)."""
    s = np.asarray(x, float) / L
    f = 1.0 / (1.0 + phi)
    N1 = f * (1 - 3 * s ** 2 + 2 * s ** 3 + phi * (1 - s))
    N2 = f * L * (s - 2 * s ** 2 + s ** 3 + 0.5 * phi * (s - s ** 2))
    N3 = f * (3 * s ** 2 - 2 * s ** 3 + phi * s)
    N4 = f * L * (-s ** 2 + s ** 3 + 0.5 * phi * (s ** 2 - s))
    return np.column_stack([N1, N2, N3, N4])


def _bending_m(rhoA, phi, L):
    """4x4 consistent translational mass for (v1, t1, v2, t2) (no rotary inertia)."""
    xg, wg = gauss_legendre(4)                                 # integrand is a degree-6 polynomial
    x = 0.5 * L * (xg + 1.0)
    N = timoshenko_shape(x, L, phi)
    return rhoA * 0.5 * L * np.einsum("g,ga,gb->ab", wg, N, N)


# index sets of the two bending planes in the 12-DOF local vector
_PLANE12 = [1, 5, 7, 11]          # u2_I, t3_I, u2_J, t3_J   (t3 = +du2/dx)
_PLANE13 = [2, 4, 8, 10]          # u3_I, t2_I, u3_J, t2_J   (t2 = -du3/dx)
_FLIP13 = np.diag([1.0, -1.0, 1.0, -1.0])


def local_stiffness(L: float, E, G, sec: Dict[str, float]) -> np.ndarray:
    """12x12 local Timoshenko stiffness (complex if E, G are complex)."""
    dtype = complex if np.iscomplexobj(np.asarray([E, G])) else float
    k = np.zeros((12, 12), dtype=dtype)
    ka = E * sec["A"] / L
    kt = G * sec["J"] / L
    for (i, j, v) in ((0, 6, ka), (3, 9, kt)):
        k[i, i] += v
        k[j, j] += v
        k[i, j] -= v
        k[j, i] -= v
    phi2 = _phi(E, G, sec["I3"], sec["As2"], L)
    phi3 = _phi(E, G, sec["I2"], sec["As3"], L)
    k[np.ix_(_PLANE12, _PLANE12)] += _bending_k(E * sec["I3"], phi2, L)
    k[np.ix_(_PLANE13, _PLANE13)] += _FLIP13 @ _bending_k(E * sec["I2"], phi3, L) @ _FLIP13
    return k


def local_mass(L: float, rho: float, E0: float, G0: float, sec: Dict[str, float]) -> np.ndarray:
    """12x12 local consistent mass (real, D-ELM-05)."""
    m = np.zeros((12, 12))
    rA = rho * sec["A"]
    rJ = rho * (sec["I2"] + sec["I3"])                          # polar inertia from I2 + I3
    for (i, j, v) in ((0, 6, rA * L), (3, 9, rJ * L)):
        m[i, i] += v / 3.0
        m[j, j] += v / 3.0
        m[i, j] += v / 6.0
        m[j, i] += v / 6.0
    for idx, I, As in ((_PLANE12, sec["I3"], sec["As2"]), (_PLANE13, sec["I2"], sec["As3"])):
        phi0 = float(np.real(_phi(E0, G0, I, As, L)))
        mb = _bending_m(rA, phi0, L)
        if idx is _PLANE13:
            mb = _FLIP13 @ mb @ _FLIP13
        m[np.ix_(idx, idx)] += mb
    return m


def releases_to_dofs(ki, kj) -> np.ndarray:
    """Local DOF indices (0..11) released at I (KI) and J (KJ); Error 10 if both ends release
    the same component."""
    ri, rj = parse_release(ki), parse_release(kj)
    both = [n for n, (a, b) in enumerate(zip(ri, rj)) if a and b]
    if both:
        names = ["P1", "P2", "P3", "M1", "M2", "M3"]
        raise ElementError("Error 10: the same end-release component is released at both I and J ("
                           + ", ".join(names[c] for c in both) + ")")
    return np.array([c for c in range(6) if ri[c]] + [6 + c for c in range(6) if rj[c]], dtype=int)


def condense_releases(k: np.ndarray, m: Optional[np.ndarray], k0: np.ndarray, rel: np.ndarray,
                      k0_mass: Optional[np.ndarray] = None):
    """Static condensation of the released local DOFs (D-ELM-04).

    Returns 12x12 matrices with zero rows/columns at the released DOFs: the element transmits no
    force there, and the released end displacement is an internal (condensed) quantity.

    Released DOFs whose stiffness is identically zero (a bending plane with I = 0) carry no force
    anyway; they are excluded from the stiffness condensation (whose ``K_rr`` they would make
    singular).  The mass transformation ``T`` is computed from ``k0_mass`` (default ``k0``), in
    which the caller supplies a regular stiffness for such planes (see :func:`_local`).
    """
    if rel.size == 0:
        return k, m
    keep = np.setdiff1d(np.arange(12), rel)
    # k0 is positive semi-definite, so a zero diagonal means a zero row and column: such a released
    # DOF (bending plane with I = 0) is decoupled from all others and needs no condensation
    stiff = rel[np.abs(np.diag(k0))[rel] > 0.0]
    k0m = k0 if k0_mass is None else k0_mass
    try:
        with blas_quiet():
            kc = np.zeros_like(k)
            kc[np.ix_(keep, keep)] = k[np.ix_(keep, keep)]
            if stiff.size:
                X = np.linalg.solve(k[np.ix_(stiff, stiff)], k[np.ix_(stiff, keep)])
                kc[np.ix_(keep, keep)] -= k[np.ix_(keep, stiff)] @ X
            mc = None
            if m is not None:
                T = np.zeros((12, keep.size))
                T[keep, np.arange(keep.size)] = 1.0
                T[rel, :] = -np.linalg.solve(k0m[np.ix_(rel, rel)], k0m[np.ix_(rel, keep)])
                mc = np.zeros_like(m)
                mc[np.ix_(keep, keep)] = T.T @ m @ T
    except np.linalg.LinAlgError:
        raise ElementError("BEAMS end releases leave a singular released block (check KI/KJ against "
                           "the section properties)") from None
    return kc, mc


def transformation(xyz) -> Tuple[np.ndarray, np.ndarray, float]:
    """(T 12x12, Lam 3x3, L) from the I, J, K coordinates."""
    xyz = np.asarray(xyz, dtype=float)
    if xyz.shape != (3, 3):
        raise ElementError("BEAMS needs the coordinates of nodes I, J and K (3 x 3)")
    Lam, L = frame_from_three_points(xyz[0], xyz[1], xyz[2])
    T = np.kron(np.eye(4), Lam)
    return T, Lam, L


def _mass_condensation_section(sec: Dict[str, float], rel: np.ndarray) -> Optional[Dict[str, float]]:
    """Section for the release transformation of the mass when a released bending plane has no
    stiffness (I = 0): the limit I -> 0+ has phi -> 0, i.e. the Euler-Bernoulli condensation,
    which is independent of the size of EI, so I = 1 with As = 0 represents it exactly."""
    ref = dict(sec)
    changed = False
    for dofs, I, As in ((_PLANE12, "I3", "As2"), (_PLANE13, "I2", "As3")):
        if sec[I] == 0.0 and np.intersect1d(rel, dofs).size:
            ref[I], ref[As] = 1.0, 0.0
            changed = True
    return ref if changed else None


def _local(xyz, mat: MaterialProps, section, ki, kj, want_mass: bool):
    sec = parse_section(section)
    T, Lam, L = transformation(xyz)
    E, G = mat.beam_moduli()
    k = local_stiffness(L, E, G, sec)
    k0 = local_stiffness(L, mat.E0, mat.G0, sec)
    m = local_mass(L, mat.rho, mat.E0, mat.G0, sec) if want_mass else None
    rel = releases_to_dofs(ki, kj)
    ref = _mass_condensation_section(sec, rel) if want_mass else None
    k0_mass = local_stiffness(L, mat.E0, mat.G0, ref) if ref is not None else None
    k, m = condense_releases(k, m, k0, rel, k0_mass)
    return k, m, T


def matrices(xyz, mat: MaterialProps, section=None, ki=None, kj=None, **_ignored):
    """(K* (12,12) complex, M (12,12) real) in global axes; xyz rows = I, J, K."""
    k, m, T = _local(xyz, mat, section, ki, kj, True)
    with blas_quiet():
        K = T.T @ k.astype(complex) @ T
        M = T.T @ m @ T
    return check_finite(0.5 * (K + K.T), "BEAMS K"), check_finite(0.5 * (M + M.T), "BEAMS M")


def recovery(xyz, mat: MaterialProps, section=None, ki=None, kj=None, **_ignored):
    """Local end forces (12, 12): f = K_l,condensed T u (forces on the element, D-STR-05)."""
    k, _, T = _local(xyz, mat, section, ki, kj, False)
    with blas_quiet():
        S = k.astype(complex) @ T
    return check_finite(S, "BEAMS recovery")


SPEC = register(ElementSpec(code=CODE, name=NAME, nnodes=NNODES, dofs=DOFS, components=list(COMPONENTS),
                            matrices=matrices, recovery=recovery,
                            description="3D Timoshenko frame (I, J; K orientation), releases, consistent mass"))
