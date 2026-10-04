"""PLANE (GROUP type 4): 4-node plane-strain quadrilateral in the global X-Z plane.

Normative sources: requirements 1.5 and 4.1 (PLANE row), D-ELM-03 (same incompatible-mode switch
as SOLID for structural PLANE elements, never for excavated soil), D-CNV-04 (complex Lame
constants), D-CNV-05 (1/2 lumped + 1/2 consistent mass), spec 08 section 4.7 (nodes I J K L
counter-clockwise, triangle by omitting / repeating L, unit thickness, 2x2 Gauss, either
orientation accepted), ARCHITECTURE 6.2 (recovery SXX SZZ TXZ at the centroid).

Formulation
-----------
2D models live in the X-Z plane (Y ignored); the DOFs are UX and UZ (labels 1 and 3).  The
bilinear map uses natural coordinates (xi, zeta) with nodes at (-1,-1), (1,-1), (1,1), (-1,1).
Plane strain (eps_yy = 0) per unit thickness::

    [s_xx]   [lam+2G  lam    0] [e_xx ]
    [s_zz] = [lam     lam+2G 0] [e_zz ]          (complex lam*, G*)
    [t_xz]   [0       0      G] [g_xz ]

Incompatible modes (Wilson QM6 with the Taylor correction, 4 internal DOFs ``1-xi^2, 1-zeta^2``
for UX and UZ) are condensed statically exactly as for SOLID; see :mod:`sassi.elements.solid`.

Orientation: "counter-clockwise" is ambiguous in the X-Z plane (spec 08 Q12).  The element
accepts either orientation: if the Jacobian at the centroid is negative the node order is
reversed internally and the matrices are returned in the caller's node order.
:func:`sassi.elements.assemble.assemble` reports such elements in ``AssembledModel.messages``
(requirements 4.1: "node order normalised to positive Jacobian with a warning"); use
:func:`orientation` / :func:`centroid_jacobian` to detect them elsewhere.
"""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np

from .base import (ElementError, ElementSpec, MaterialProps, blas_quiet, check_finite, gauss_legendre,
                   mixed_mass, register)

NAME = "PLANE"
CODE = 4
DOFS = (1, 3)
NNODES = 4
COMPONENTS = ["SXX", "SZZ", "TXZ"]

NODE_NAT = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype=float)
GAUSS_POINTS = 2

_D_LAM = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
_D_G = np.diag([2.0, 2.0, 1.0])
_REVERSE = np.array([0, 3, 2, 1])


def plane_strain_D(lam: complex, G: complex) -> np.ndarray:
    """3x3 plane-strain elasticity matrix (xx, zz, xz) for (complex) Lame constants."""
    return lam * _D_LAM + G * _D_G


def quad4_shape(pts: np.ndarray):
    """Bilinear shape functions N (np, 4) and natural derivatives dN (np, 4, 2)."""
    pts = np.atleast_2d(pts)
    f = 1.0 + pts[:, None, :] * NODE_NAT[None, :, :]
    N = f.prod(axis=2) / 4.0
    dN = np.empty(f.shape)
    dN[..., 0] = NODE_NAT[None, :, 0] * f[..., 1] / 4.0
    dN[..., 1] = NODE_NAT[None, :, 1] * f[..., 0] / 4.0
    return N, dN


def gauss_2d(n: int):
    x, w = gauss_legendre(n)
    X, Y = np.meshgrid(x, x, indexing="ij")
    return np.column_stack([X.ravel(), Y.ravel()]), (w[:, None] * w[None, :]).ravel()


def strain_matrix_2d(dNx: np.ndarray) -> np.ndarray:
    """B (..., 3, 2n) for strains (e_xx, e_zz, g_xz) from physical derivatives dNx (..., n, 2)."""
    shp = dNx.shape[:-2]
    n = dNx.shape[-2]
    B = np.zeros(shp + (3, 2 * n), dtype=dNx.dtype)
    dx, dz = dNx[..., 0], dNx[..., 1]
    B[..., 0, 0::2] = dx
    B[..., 1, 1::2] = dz
    B[..., 2, 0::2] = dz
    B[..., 2, 1::2] = dx
    return B


def _xz(xyz: np.ndarray) -> np.ndarray:
    xyz = np.asarray(xyz, dtype=float)
    if xyz.shape[-1] == 3:
        return xyz[..., [0, 2]]
    if xyz.shape[-1] == 2:
        return xyz
    raise ElementError("PLANE coordinates must have 2 (x, z) or 3 (x, y, z) columns")


def _as_quads(xz: np.ndarray) -> np.ndarray:
    if xz.ndim == 2:
        xz = xz[None]
    if xz.shape[1] == 3:                                      # triangle: repeat node K as L
        xz = np.concatenate([xz, xz[:, 2:3, :]], axis=1)
    if xz.shape[1] != 4:
        raise ElementError("PLANE needs 4 nodes (or 3 for a triangle)")
    return xz


def centroid_jacobian(xyz) -> np.ndarray:
    """Jacobian determinant at the centroid, (nE,), for xyz (nE, 4|3, 3|2) or one element.

    It equals a quarter of the signed area of the quadrilateral in the (x right, z up) view:
    positive for counter-clockwise I-J-K-L, negative for clockwise input (spec 08 Q12).
    """
    xz = _as_quads(_xz(xyz))
    _, dN0 = quad4_shape(np.zeros((1, 2)))
    return np.linalg.det(np.einsum("ai,eaj->eij", dN0[0], xz))


def orientation(xyz) -> int:
    """+1 if I-J-K-L is counter-clockwise in the (x right, z up) view, -1 otherwise."""
    return 1 if centroid_jacobian(np.asarray(xyz, float)[None])[0] > 0 else -1


def plane_kernel(xyz, lam, G, rho, incompatible: bool = False, mass: str = "mixed",
                 want_mass: bool = True, want_recovery: bool = True) -> Dict[str, np.ndarray]:
    """Vectorised PLANE matrices for nE elements.

    xyz (nE, 4, 3) or (nE, 4, 2) [x, z]; lam, G complex (nE,); rho (nE,).  Returns K (nE,8,8)
    complex, M (nE,8,8) real, S (nE,3,8) complex centroid stresses and B (nE,3,8) strains, all in
    the caller's node order.
    """
    xz = _as_quads(_xz(xyz))
    nE = xz.shape[0]
    lam = np.broadcast_to(np.asarray(lam, dtype=complex), (nE,))
    G = np.broadcast_to(np.asarray(G, dtype=complex), (nE,))
    rho = np.broadcast_to(np.asarray(rho, dtype=float), (nE,))
    # orientation normalisation (spec 08 Q12); the assembler reports the flipped elements
    flip = centroid_jacobian(xz) < 0
    xzw = xz.copy()
    xzw[flip] = xz[flip][:, _REVERSE, :]
    with blas_quiet():
        out = _plane_kernel(xzw, lam, G, rho, incompatible, mass, want_mass, want_recovery)
    if np.any(flip):
        p2 = np.column_stack([2 * _REVERSE, 2 * _REVERSE + 1]).ravel()
        for k, v in out.items():
            w = v[flip]
            back = np.empty_like(w)
            if k in ("K", "M"):
                back[:, p2[:, None], p2[None, :]] = w
            else:
                back[:, :, p2] = w
            v[flip] = back
    for k, v in out.items():
        check_finite(v, f"PLANE {k}")
    return out


def _plane_kernel(xz, lam, G, rho, incompatible, mass, want_mass, want_recovery):
    nE = xz.shape[0]
    pts, wts = gauss_2d(GAUSS_POINTS)
    _, dN = quad4_shape(pts)
    J = np.einsum("gai,eaj->egij", dN, xz)
    det = np.linalg.det(J)
    span = np.linalg.norm(xz.max(axis=1) - xz.min(axis=1), axis=1)
    if np.any(det <= 1e-12 * (span ** 2)[:, None]):
        raise ElementError("EDU-05: non-positive Jacobian at a Gauss point of a PLANE element "
                           "(distorted or self-intersecting quadrilateral)")
    Jinv = np.linalg.inv(J)
    Bu = strain_matrix_2d(np.einsum("egjk,gak->egaj", Jinv, dN))     # (nE, ng, 3, 8)

    J0 = np.einsum("ai,eaj->eij", quad4_shape(np.zeros((1, 2)))[1][0], xz)
    det0 = np.linalg.det(J0)
    J0inv = np.linalg.inv(J0)
    Bu0 = strain_matrix_2d(np.einsum("ejk,ak->eaj", J0inv, quad4_shape(np.zeros((1, 2)))[1][0]))

    if incompatible:
        dNa = np.zeros((pts.shape[0], 2, 2))
        dNa[:, 0, 0] = -2.0 * pts[:, 0]                          # d(1-xi^2)/dxi
        dNa[:, 1, 1] = -2.0 * pts[:, 1]                          # d(1-zeta^2)/dzeta
        scale = det0[:, None] / det                              # Taylor correction
        Ba = strain_matrix_2d(np.einsum("ejk,gak->egaj", J0inv, dNa) * scale[:, :, None, None])
        Bf = np.concatenate([Bu, Ba], axis=3)
        Ba0 = np.zeros((nE, 3, 4))                               # mode derivatives vanish at the centroid
    else:
        Bf = Bu
    wdet = det * wts[None, :]
    div = Bf[:, :, 0, :] + Bf[:, :, 1, :]
    K_lam = np.einsum("eg,ega,egb->eab", wdet, div, div)
    K_G = np.einsum("eg,egia,i,egib->eab", wdet, Bf, np.diag(_D_G), Bf)
    Kf = lam[:, None, None] * K_lam + G[:, None, None] * K_G
    if incompatible:
        Kuu, Kua, Kau, Kaa = Kf[:, :8, :8], Kf[:, :8, 8:], Kf[:, 8:, :8], Kf[:, 8:, 8:]
        X = np.linalg.solve(Kaa, Kau)
        K = Kuu - Kua @ X
        Brec = Bu0.astype(complex) - Ba0 @ X
    else:
        K = Kf
        Brec = Bu0.astype(complex)
    out = {"K": 0.5 * (K + np.swapaxes(K, 1, 2))}
    if want_mass:
        N, _ = quad4_shape(pts)                                  # 2x2 Gauss is exact for the bilinear mass
        Mc4 = np.einsum("eg,ga,gb->eab", wdet * rho[:, None], N, N)
        Mc = np.einsum("eab,ij->eaibj", Mc4, np.eye(2)).reshape(nE, 8, 8)
        out["M"] = mixed_mass(Mc, mass)
    if want_recovery:
        D = lam[:, None, None] * _D_LAM[None] + G[:, None, None] * _D_G[None]
        out["S"] = D @ Brec
        out["B"] = Brec
    return out


def _fold_triangle(K: np.ndarray, n_in: int) -> np.ndarray:
    """Sum the repeated 4th node into node 3 when the caller gave only 3 nodes."""
    if n_in != 3:
        return K
    K = K.copy()
    if K.ndim == 2 and K.shape[0] == K.shape[1] == 8:
        K[4:6, :] += K[6:8, :]
        K[:, 4:6] += K[:, 6:8]
        return K[:6, :6]
    K[:, 4:6] += K[:, 6:8]
    return K[:, :6]


def matrices(xyz, mat: MaterialProps, incompatible: bool = False, mass: str = "mixed", **_ignored):
    """(K* (8,8) complex, M (8,8) real) of one PLANE element (6x6 for a 3-node triangle)."""
    n_in = np.asarray(xyz).shape[0]
    r = plane_kernel(np.asarray(xyz, float)[None], mat.lam, mat.G, mat.rho, incompatible=incompatible,
                     mass=mass, want_recovery=False)
    return _fold_triangle(r["K"][0], n_in), _fold_triangle(r["M"][0], n_in)


def recovery(xyz, mat: MaterialProps, incompatible: bool = False, quantity: str = "stress", **_ignored):
    """Centroid stresses SXX SZZ TXZ (3, 8) or strains EXX EZZ GXZ with ``quantity='strain'``."""
    n_in = np.asarray(xyz).shape[0]
    r = plane_kernel(np.asarray(xyz, float)[None], mat.lam, mat.G, mat.rho, incompatible=incompatible,
                     want_mass=False)
    S = r["B"][0] if quantity == "strain" else r["S"][0]
    return _fold_triangle(S, n_in)


def batch(xyz, mats: Sequence[MaterialProps], incompatible: bool = False, mass: str = "mixed", **_ignored):
    """dict(K, M, S, B) for nE PLANE elements (4 nodes each, repeated node for triangles)."""
    lam = np.array([m.lam for m in mats], dtype=complex)
    G = np.array([m.G for m in mats], dtype=complex)
    rho = np.array([m.rho for m in mats], dtype=float)
    return plane_kernel(xyz, lam, G, rho, incompatible=incompatible, mass=mass)


SPEC = register(ElementSpec(code=CODE, name=NAME, nnodes=NNODES, dofs=DOFS, components=list(COMPONENTS),
                            matrices=matrices, recovery=recovery, batch=batch,
                            description="4-node plane strain (X-Z), optional 4 incompatible modes, 1/2+1/2 mass"))
