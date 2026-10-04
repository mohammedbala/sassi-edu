"""SHELL (GROUP type 3): flat-facet thin (Kirchhoff) shell with no drilling stiffness.

Normative sources: requirements 1.5 and 4.1 (SHELL row), D-ELM-06 (membrane = plane-stress Q4 with
incompatible modes, CST for triangles; bending = DKQ quad / DKT triangle; zero drilling stiffness),
D-ELM-07 (lumped mass rho t A / n on translations, no rotary inertia), D-CNV-04 (E* = E c(beta),
nu real), spec 08 section 4.6 (nodes I J K L counter-clockwise, L omitted = triangle; local axes),
D-STR-10 and spec 05d (recovery: membrane stresses Sx'x' Sy'y' Sx'y' and moments per unit length
Mx'x' My'y' Mx'y' at the centroid in local axes).

Local axes (spec 08 4.6)::

    x' = unit(m_JK - m_LI)        (m_AB = midpoint of side AB; triangles use L = K)
    z' = unit(x' x (m_KL - m_IJ))  (normal; I->J->K->L counter-clockwise seen from +z')
    y' = z' x x'

Warped quadrilaterals are projected on the mean plane through the centroid (the out-of-plane
offsets are ignored; CHECK reports W3/E12 by the D-CHK-05 thresholds).  A triangle is given with
3 nodes, or with 4 nodes of which one pair of adjacent corners coincides (L = K, also L = I, J = K,
...; :func:`corner_rows`); any other coincidence of corners is CHECK Error 9 (spec 08 5.5).

Membrane (plane stress, DOFs u', v'): ``N = D_m eps`` with
``D_m = E t/(1-nu^2) [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]``; quads use the bilinear Q4 with the
4 Wilson-Taylor incompatible modes (1-xi^2, 1-eta^2 for u' and v') condensed statically (exact in
pure in-plane bending of rectangles), triangles the constant-strain triangle.

Plate bending (DOFs w', theta_x', theta_y'): Discrete Kirchhoff elements (Batoz, Bathe & Ho 1980;
Batoz & Ben Tahar 1982).  The rotations of the normal ``beta_x = -dw/dx = theta_y``,
``beta_y = -dw/dy = -theta_x`` are interpolated quadratically (8-node serendipity for DKQ, 6-node
triangle for DKT).  The mid-side values are eliminated with the discrete Kirchhoff constraints:

* along each side w is cubic, so the tangential rotation at the mid-side is
  ``beta_s,k = 3/(2L) (w_i - w_j) - 1/4 (beta_s,i + beta_s,j)`` (i.e. beta_s = -dw/ds there);
* the normal rotation varies linearly: ``beta_n,k = (beta_n,i + beta_n,j)/2``;
* at the corners beta = -grad(w) (Kirchhoff).

The curvatures ``kappa = [beta_x,x, beta_y,y, beta_x,y + beta_y,x]`` and moments
``M = D_b kappa``, ``D_b = E t^3 / (12 (1 - nu^2)) [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]``, so
``Mx'x' = integral(sigma_x'x' z' dz')`` (positive = tension on the +z' face).  Transverse shear
deformation is neglected (thin plate).  The rotation about z' (drilling) has **no stiffness**: it must
be restrained with FIXROT / FIXSHLROT (EDU-06).

Global matrices ``K = T^T K_l T`` with ``T = blockdiag(Lam, ..., Lam)`` (one Lam per translation and
rotation triple, rows of Lam = x', y', z').  Because the formulation is linear in E with nu real,
``K* = c(beta_s) K0`` exactly.
"""
from __future__ import annotations

from typing import List, Tuple

import numpy as np

from .base import ElementError, ElementSpec, MaterialProps, blas_quiet, check_finite, gauss_legendre, register
from .plane import quad4_shape, strain_matrix_2d

NAME = "SHELL"
CODE = 3
DOFS = (1, 2, 3, 4, 5, 6)
NNODES = 4
COMPONENTS = ["FXX", "FYY", "FXY", "MXX", "MYY", "MXY"]

#: integration points per direction of the DKQ bending stiffness (Batoz & Ben Tahar: 2 x 2)
DKQ_GAUSS = 2


# ======================================================================================
# geometry
# ======================================================================================
#: two corners closer than this fraction of the element size are coincident
_COINCIDENT_RTOL = 1e-12


def corner_rows(xyz) -> List[int]:
    """Rows of ``xyz`` that are the distinct corners of the facet (CHECK Error 9 otherwise).

    * 3 rows: a triangle; any coincident pair is Error 9.
    * 4 rows, all distinct: a quadrilateral, ``[0, 1, 2, 3]``.
    * 4 rows with exactly one pair of *cyclically adjacent* coincident corners, e.g. L = K (the
      usual triangle input), L = I or J = K: the triangle of the three distinct corners in the
      order of their first appearance, e.g. ``[0, 1, 2]`` for (I, J, K, I) and ``[0, 1, 3]`` for
      (I, J, J, L).  The circulation (normal z') is unchanged.
    * any other coincidence (two pairs, or opposite corners I = K / J = L, which fold the facet
      onto itself) is Error 9.
    """
    X = np.asarray(xyz, dtype=float)
    if X.ndim != 2 or X.shape[1] != 3 or X.shape[0] not in (3, 4):
        raise ElementError("SHELL needs 3 or 4 nodes with 3 coordinates")
    n = X.shape[0]
    span = max(float(np.ptp(X, axis=0).max()), 1e-300)
    same = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2) <= _COINCIDENT_RTOL * span
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n) if same[a, b]]
    if not pairs:
        return list(range(n))
    if n == 4 and len(pairs) == 1 and (pairs[0][1] - pairs[0][0]) in (1, 3):
        return [a for a in range(4) if a != pairs[0][1]]           # drop the later of the two
    raise ElementError("Error 9: SHELL element with coincident corner nodes (corners "
                       + ", ".join(f"{'IJKL'[a]}={'IJKL'[b]}" for a, b in pairs)
                       + "); only one pair of adjacent corners may coincide, giving a triangle")


def is_triangle(xyz: np.ndarray) -> bool:
    """True if the facet has three distinct corners (3 rows, or 4 rows with one adjacent pair
    coincident, see :func:`corner_rows`)."""
    return len(corner_rows(xyz)) == 3


def local_frame(xyz) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(Lam (3,3) rows x' y' z', centroid (3,), xy (n,2) projected local coordinates).

    ``xy`` has one row per distinct corner (:func:`corner_rows`): 3 for a triangle, 4 for a quad.
    """
    X = np.asarray(xyz, dtype=float)
    keep = corner_rows(X)
    tri = len(keep) == 3
    X = X[keep]
    I, J, K = X[0], X[1], X[2]
    L = K if tri else X[3]
    xp = 0.5 * (J + K) - 0.5 * (L + I)
    nx = np.linalg.norm(xp)
    if nx <= 0:
        raise ElementError("degenerate SHELL element (zero x' axis)")
    xp = xp / nx
    zp = np.cross(xp, 0.5 * (K + L) - 0.5 * (I + J))
    nz = np.linalg.norm(zp)
    if nz <= 1e-12 * max(float(np.ptp(X, axis=0).max()), 1e-300):
        raise ElementError("Error 9: SHELL element with collinear nodes (degenerate facet)")
    zp = zp / nz
    yp = np.cross(zp, xp)
    Lam = np.vstack([xp, yp, zp])
    c = X.mean(axis=0)
    xy = (X - c) @ Lam[:2].T
    return Lam, c, xy


def polygon_area(xy: np.ndarray) -> float:
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def _Dplate(nu: float) -> np.ndarray:
    return np.array([[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, 0.5 * (1.0 - nu)]])


# ======================================================================================
# membrane
# ======================================================================================
def membrane_q4(xy: np.ndarray, Dm: np.ndarray, incompatible: bool = True):
    """Plane-stress Q4 membrane (8x8) with condensed incompatible modes, and the centroid strain
    operator (3x8).  ``Dm`` is the membrane rigidity (force per unit length)."""
    x, w = gauss_legendre(2)
    pts = np.array([[a, b] for a in x for b in x])
    wts = np.array([wa * wb for wa in w for wb in w])
    _, dN = quad4_shape(pts)
    J = np.einsum("gai,aj->gij", dN, xy)
    det = np.linalg.det(J)
    if np.any(det <= 0):
        raise ElementError("EDU-05: non-positive Jacobian in a SHELL membrane (re-entrant quadrilateral?)")
    Bu = strain_matrix_2d(np.einsum("gjk,gak->gaj", np.linalg.inv(J), dN))
    _, dN0 = quad4_shape(np.zeros((1, 2)))
    J0 = dN0[0].T @ xy
    det0 = np.linalg.det(J0)
    J0inv = np.linalg.inv(J0)
    B0 = strain_matrix_2d((J0inv @ dN0[0].T).T[None])[0]
    if incompatible:
        dNa = np.zeros((4, 2, 2))
        dNa[:, 0, 0] = -2.0 * pts[:, 0]
        dNa[:, 1, 1] = -2.0 * pts[:, 1]
        Ba = strain_matrix_2d(np.einsum("jk,gak->gaj", J0inv, dNa) * (det0 / det)[:, None, None])
        Bf = np.concatenate([Bu, Ba], axis=2)
    else:
        Bf = Bu
    Kf = np.einsum("g,gia,ij,gjb->ab", wts * det, Bf, Dm, Bf)
    if incompatible:
        X = np.linalg.solve(Kf[8:, 8:], Kf[8:, :8])
        K = Kf[:8, :8] - Kf[:8, 8:] @ X
        # recovery includes the condensed modes: eps = B0 u + Ba(0) alpha, alpha = -X u; the mode
        # derivatives (-2 xi, -2 eta) vanish at the centroid, so Ba(0) = 0 here
        Ba0 = np.zeros((3, 4))
        Brec = B0 - Ba0 @ X
    else:
        K, Brec = Kf, B0
    return 0.5 * (K + K.T), Brec


def membrane_cst(xy: np.ndarray, Dm: np.ndarray):
    """Constant-strain triangle membrane (6x6) and its strain operator (3x6)."""
    x, y = xy[:, 0], xy[:, 1]
    A = polygon_area(xy)
    if A <= 0:
        raise ElementError("degenerate SHELL triangle (zero or negative area)")
    b = np.array([y[1] - y[2], y[2] - y[0], y[0] - y[1]])
    c = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]])
    B = np.zeros((3, 6))
    B[0, 0::2] = b
    B[1, 1::2] = c
    B[2, 0::2] = c
    B[2, 1::2] = b
    B /= 2.0 * A
    return A * B.T @ Dm @ B, B


# ======================================================================================
# Discrete Kirchhoff plate bending
# ======================================================================================
def _dk_constraints(xy: np.ndarray, sides) -> np.ndarray:
    """Matrix A (2 nb, 3 nc) giving [beta_x(nodes), beta_y(nodes)] of the quadratic rotation field
    (corners then mid-sides) from the corner DOFs [w, theta_x, theta_y] per corner."""
    nc = xy.shape[0]
    nb = nc + len(sides)
    span = max(float(np.ptp(xy, axis=0).max()), 1e-300)
    A = np.zeros((2 * nb, 3 * nc))
    for i in range(nc):
        A[i, 3 * i + 2] = 1.0          # beta_x = theta_y
        A[nb + i, 3 * i + 1] = -1.0    # beta_y = -theta_x
    for k, (i, j) in enumerate(sides):
        m = nc + k
        d = xy[j] - xy[i]
        L = float(np.hypot(d[0], d[1]))
        if L <= _COINCIDENT_RTOL * span:
            raise ElementError("Error 9: SHELL element with a zero-length side (coincident nodes)")
        t = d / L
        n = np.array([-t[1], t[0]])
        bsum_x = A[i] + A[j]           # beta_x,i + beta_x,j (rows)
        bsum_y = A[nb + i] + A[nb + j]
        w_diff = np.zeros(3 * nc)
        w_diff[3 * i] = 1.0
        w_diff[3 * j] = -1.0
        beta_s = 1.5 / L * w_diff - 0.25 * (t[0] * bsum_x + t[1] * bsum_y)
        beta_n = 0.5 * (n[0] * bsum_x + n[1] * bsum_y)
        A[m] = t[0] * beta_s + n[0] * beta_n
        A[nb + m] = t[1] * beta_s + n[1] * beta_n
    return A


def _curvature_B(dNx: np.ndarray) -> np.ndarray:
    """Curvature operator (3, 2 nb) acting on [beta_x(nodes), beta_y(nodes)]."""
    nb = dNx.shape[0]
    B = np.zeros((3, 2 * nb))
    B[0, :nb] = dNx[:, 0]
    B[1, nb:] = dNx[:, 1]
    B[2, :nb] = dNx[:, 1]
    B[2, nb:] = dNx[:, 0]
    return B


_SER_NAT = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1], [0, -1], [1, 0], [0, 1], [-1, 0]], dtype=float)


def serendipity8_dshape(p: np.ndarray) -> np.ndarray:
    """Natural derivatives (8, 2) of the 8-node serendipity functions at point p = (xi, eta)."""
    xi, eta = float(p[0]), float(p[1])
    d = np.zeros((8, 2))
    for a in range(4):
        xa, ya = _SER_NAT[a]
        d[a, 0] = 0.25 * xa * (1 + eta * ya) * (2 * xi * xa + eta * ya)
        d[a, 1] = 0.25 * ya * (1 + xi * xa) * (xi * xa + 2 * eta * ya)
    for a in (4, 6):                          # xi_a = 0
        ya = _SER_NAT[a, 1]
        d[a, 0] = -xi * (1 + eta * ya)
        d[a, 1] = 0.5 * (1 - xi * xi) * ya
    for a in (5, 7):                          # eta_a = 0
        xa = _SER_NAT[a, 0]
        d[a, 0] = 0.5 * xa * (1 - eta * eta)
        d[a, 1] = -eta * (1 + xi * xa)
    return d


def dkq(xy: np.ndarray, Db: np.ndarray):
    """DKQ bending stiffness (12x12) and centroid curvature operator (3x12) of a quadrilateral."""
    A = _dk_constraints(xy, [(0, 1), (1, 2), (2, 3), (3, 0)])
    x, w = gauss_legendre(DKQ_GAUSS)
    K = np.zeros((12, 12))
    for a, wa in zip(x, w):
        for b, wb in zip(x, w):
            p = np.array([a, b])
            _, dN4 = quad4_shape(p[None])
            J = dN4[0].T @ xy
            det = np.linalg.det(J)
            if det <= 0:
                raise ElementError("EDU-05: non-positive Jacobian in a SHELL (DKQ) element")
            dNx = (np.linalg.inv(J) @ serendipity8_dshape(p).T).T
            BA = _curvature_B(dNx) @ A
            K += wa * wb * det * BA.T @ Db @ BA
    _, dN4 = quad4_shape(np.zeros((1, 2)))
    J0 = dN4[0].T @ xy
    B0 = _curvature_B((np.linalg.inv(J0) @ serendipity8_dshape(np.zeros(2)).T).T) @ A
    return 0.5 * (K + K.T), B0


def tri6_dshape(xi: float, eta: float) -> np.ndarray:
    """Derivatives (6, 2) of the quadratic triangle functions (corners 1-3, mid-sides 12, 23, 31)."""
    L1, L2, L3 = 1.0 - xi - eta, xi, eta
    dL = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])   # dL_i/d(xi, eta)
    Ls = (L1, L2, L3)
    d = np.zeros((6, 2))
    for i in range(3):
        d[i] = (4.0 * Ls[i] - 1.0) * dL[i]
    for k, (i, j) in enumerate(((0, 1), (1, 2), (2, 0))):
        d[3 + k] = 4.0 * (Ls[i] * dL[j] + Ls[j] * dL[i])
    return d


def dkt(xy: np.ndarray, Db: np.ndarray):
    """DKT bending stiffness (9x9) and centroid curvature operator (3x9) of a triangle."""
    A = _dk_constraints(xy, [(0, 1), (1, 2), (2, 0)])
    J = np.array([xy[1] - xy[0], xy[2] - xy[0]])            # d(x, y)/d(xi, eta), constant
    det = np.linalg.det(J)
    if det <= 0:
        raise ElementError("degenerate SHELL triangle (zero or negative area)")
    Jinv = np.linalg.inv(J)
    K = np.zeros((9, 9))
    for (a, b) in ((0.5, 0.0), (0.5, 0.5), (0.0, 0.5)):      # mid-side rule, exact for quadratics
        BA = _curvature_B((Jinv @ tri6_dshape(a, b).T).T) @ A
        K += (det / 6.0) * BA.T @ Db @ BA
    B0 = _curvature_B((Jinv @ tri6_dshape(1.0 / 3.0, 1.0 / 3.0).T).T) @ A
    return 0.5 * (K + K.T), B0


# ======================================================================================
# element
# ======================================================================================
def _local(xyz, mat: MaterialProps, thick: float, incompatible: bool = True):
    """Real local matrices: K0 (6n x 6n), recovery operators, area, Lam, corner rows, thickness.

    ``n`` is the number of distinct corners (:func:`corner_rows`).  ``incompatible`` is the SHELL
    membrane switch (on by default, D-ELM-06), not MOPT <incomp>.
    """
    if thick is None or float(thick) <= 0:
        raise ElementError("SHELL thickness must be > 0 (THICK)")
    t = float(thick)
    keep = corner_rows(xyz)
    Lam, _, xy = local_frame(xyz)
    tri = xy.shape[0] == 3
    n = xy.shape[0]
    E0, nu = mat.E0, mat.nu0
    Dm = E0 * t / (1.0 - nu * nu) * _Dplate(nu)
    Db = E0 * t ** 3 / (12.0 * (1.0 - nu * nu)) * _Dplate(nu)
    with blas_quiet():
        if tri:
            Km, Bm = membrane_cst(xy, Dm)
            Kb, Bb = dkt(xy, Db)
        else:
            Km, Bm = membrane_q4(xy, Dm, incompatible=incompatible)
            Kb, Bb = dkq(xy, Db)
    mem = np.array([[6 * a, 6 * a + 1] for a in range(n)]).ravel()
    ben = np.array([[6 * a + 2, 6 * a + 3, 6 * a + 4] for a in range(n)]).ravel()
    K = np.zeros((6 * n, 6 * n))
    K[np.ix_(mem, mem)] = Km
    K[np.ix_(ben, ben)] = Kb
    R = np.zeros((6, 6 * n))
    R[0:3, mem] = Dm @ Bm                 # membrane forces per unit length N = D_m eps
    R[3:6, ben] = Db @ Bb                 # moments per unit length M = D_b kappa
    return K, R, polygon_area(xy), Lam, keep, t


def _expand(a: np.ndarray, keep: List[int], n_rows: int, square: bool) -> np.ndarray:
    """Place the node blocks of a matrix built on the distinct corners ``keep`` into the layout
    of the ``n_rows`` input rows: rows and columns (``square``, K and M) or columns only (S).
    A repeated corner gets zero rows/columns."""
    if len(keep) == n_rows:
        return a
    idx = np.concatenate([np.arange(6 * r, 6 * r + 6) for r in keep])
    if square:
        out = np.zeros((6 * n_rows, 6 * n_rows), dtype=a.dtype)
        out[np.ix_(idx, idx)] = a
    else:
        out = np.zeros((a.shape[0], 6 * n_rows), dtype=a.dtype)
        out[:, idx] = a
    return out


def matrices(xyz, mat: MaterialProps, thick: float = None, membrane_incompatible: bool = True, **_ignored):
    """(K* complex, M real) in global axes, size 6 x (number of rows of xyz); a triangle given with
    4 rows (one adjacent corner repeated, e.g. L = K, see :func:`corner_rows`) returns zero
    rows/columns for the repeated corner.

    ``membrane_incompatible`` (default True, D-ELM-06) keeps the incompatible modes of the Q4
    membrane; it is independent of the SOLID/PLANE switch MOPT <incomp> (``incompatible`` is ignored).
    """
    X = np.asarray(xyz, dtype=float)
    K0, _, A, Lam, keep, t = _local(X, mat, thick, membrane_incompatible)
    n = len(keep)
    T = np.kron(np.eye(2 * n), Lam)
    cE = mat.shell_moduli()[0] / mat.E0                      # c(beta_s): E* = E0 c(beta)
    with blas_quiet():
        K = cE * (T.T @ K0 @ T)
    M = np.zeros((6 * n, 6 * n))
    mnode = mat.rho * t * A / n                               # D-ELM-07
    for a in range(n):
        for d in range(3):
            M[6 * a + d, 6 * a + d] = mnode
    K, M = _expand(K, keep, X.shape[0], True), _expand(M, keep, X.shape[0], True)
    return check_finite(0.5 * (K + K.T), "SHELL K"), M


def recovery(xyz, mat: MaterialProps, thick: float = None, membrane_incompatible: bool = True,
             membrane_output: str = "stress", **_ignored):
    """Centroid resultants (6, 6n) in local axes: FXX FYY FXY and MXX MYY MXY.

    ``membrane_output='stress'`` (default, D-STR-10 / spec 05d: membrane stresses F/L^2) or
    ``'force'`` (membrane forces per unit length F/L).  Moments are always per unit length.
    """
    X = np.asarray(xyz, dtype=float)
    _, R, _, Lam, keep, t = _local(X, mat, thick, membrane_incompatible)
    n = len(keep)
    if membrane_output == "stress":
        R[0:3] /= t
    elif membrane_output != "force":
        raise ElementError("membrane_output must be 'stress' or 'force'")
    T = np.kron(np.eye(2 * n), Lam)
    cE = mat.shell_moduli()[0] / mat.E0
    with blas_quiet():
        S = cE * (R @ T)
    S = _expand(S, keep, X.shape[0], False)
    return check_finite(S, "SHELL recovery")


SPEC = register(ElementSpec(code=CODE, name=NAME, nnodes=NNODES, dofs=DOFS, components=list(COMPONENTS),
                            matrices=matrices, recovery=recovery,
                            description="flat Kirchhoff shell: Q4+incompatible/CST membrane, DKQ/DKT bending, "
                                        "no drilling, lumped mass"))
