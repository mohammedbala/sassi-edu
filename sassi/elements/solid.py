"""SOLID (GROUP type 1): 8-node isoparametric hexahedron with optional incompatible modes.

Normative sources: requirements 1.5 and 4.1 (SOLID row), D-ELM-01 (Gauss rule by EINT),
D-ELM-02 (Wilson-Taylor incompatible modes, structural SOLIDs only, condensed from K*),
D-CNV-04 (complex Lame constants), D-CNV-05 (1/2 lumped + 1/2 consistent mass), EDU-05 (positive
Jacobian), spec 08 section 4.4 (node numbering, degenerate shapes), ARCHITECTURE 6.2 (recovery:
centroid stresses SXX SYY SZZ SXY SXZ SYZ in global axes, SOCT row zero).

Formulation (for the student)
-----------------------------
Nodes 1-2-3-4 form one face and 5-6-7-8 the opposite face (node 4+i opposite node i).  Natural
coordinates of the nodes are (xi, eta, zeta) = (-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1),
(-1,-1,1), (1,-1,1), (1,1,1), (-1,1,1) and the trilinear shape functions are
``N_a = (1 + xi xi_a)(1 + eta eta_a)(1 + zeta zeta_a) / 8``.  Prisms and pyramids are obtained by
repeating node numbers (7 = 8; 5 = 6 and 7 = 8; 5 = 6 = 7 = 8): the repeated nodes simply carry
the same coordinates and assemble into the same global equations.

Stress-strain law with complex Lame constants (hysteretic damping, requirements 4.0.2)::

    sigma = lam* tr(eps) I + 2 G* eps   ->   D* = lam* D_lam + G* D_G

so the element stiffness is ``K* = lam* K_lam + G* K_G`` with two real matrices.  When
beta_p = beta_s the whole matrix is simply ``c(beta) K0``.

Incompatible modes (D-ELM-02).  Trilinear bricks are far too stiff in bending ("shear locking"):
a pure bending field needs a quadratic displacement that the element cannot represent.  Wilson
added the internal modes ``1 - xi^2, 1 - eta^2, 1 - zeta^2`` for each of the 3 displacement
components (9 internal DOFs alpha).  Their derivatives are evaluated with the *centroid*
Jacobian J0 and scaled by det J0 / det J (Taylor, Beresford & Wilson 1976), which makes
``integral(B_alpha dV) = 0`` for any shape: a constant stress state then does no work on the
internal modes, so the element passes the patch test.  The alpha DOFs are condensed statically
from the complex matrix::

    K = K_uu - K_ua K_aa^-1 K_au          alpha = -K_aa^-1 K_au u

and the stress recovery uses the full strain ``B_u u + B_a alpha``.  (At the centroid the
derivatives of 1 - xi^2 vanish, so the internal modes do not change the centroid strain operator,
but they change the displacements that are solved for.)

Mass (D-CNV-05): ``M = 1/2 diag(row sums of Mc) + 1/2 Mc`` with the consistent mass
``Mc = integral(rho N^T N dV)``; the ``mass`` keyword selects 'consistent' or 'lumped' for studies
(VP-03).
"""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np

from .base import (ElementError, ElementSpec, MaterialProps, blas_quiet, check_finite, gauss_legendre,
                   mixed_mass, register)

NAME = "SOLID"
CODE = 1
DOFS = (1, 2, 3)
NNODES = 8
COMPONENTS = ["SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ", "SOCT"]

#: Natural coordinates of the 8 nodes (spec 08 section 4.4).
NODE_NAT = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                     [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=float)

#: Gauss points per direction for EINT 0/1/2 (D-ELM-01).
EINT_POINTS = {0: 2, 1: 3, 2: 4}

# Voigt order of the strain/stress vector: xx, yy, zz, xy, xz, yz (engineering shear strains)
_D_LAM = np.zeros((6, 6))
_D_LAM[:3, :3] = 1.0
_D_G = np.diag([2.0, 2.0, 2.0, 1.0, 1.0, 1.0])


def isotropic_D(lam: complex, G: complex) -> np.ndarray:
    """6x6 isotropic elasticity matrix (Voigt xx yy zz xy xz yz) for (complex) Lame constants."""
    return lam * _D_LAM + G * _D_G


def hex8_shape(pts: np.ndarray):
    """Trilinear shape functions N (np, 8) and natural derivatives dN (np, 8, 3) at ``pts``."""
    pts = np.atleast_2d(pts)
    f = 1.0 + pts[:, None, :] * NODE_NAT[None, :, :]          # (np, 8, 3)
    N = f.prod(axis=2) / 8.0
    dN = np.empty(f.shape)
    dN[..., 0] = NODE_NAT[None, :, 0] * f[..., 1] * f[..., 2] / 8.0
    dN[..., 1] = NODE_NAT[None, :, 1] * f[..., 0] * f[..., 2] / 8.0
    dN[..., 2] = NODE_NAT[None, :, 2] * f[..., 0] * f[..., 1] / 8.0
    return N, dN


def gauss_3d(n: int):
    """Tensor-product n x n x n Gauss rule: points (n^3, 3), weights (n^3,)."""
    x, w = gauss_legendre(n)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    W = (w[:, None, None] * w[None, :, None] * w[None, None, :])
    return np.column_stack([X.ravel(), Y.ravel(), Z.ravel()]), W.ravel()


def _incompatible_dshape(pts: np.ndarray) -> np.ndarray:
    """Natural derivatives (np, 3 modes, 3) of the Wilson modes 1-xi^2, 1-eta^2, 1-zeta^2."""
    pts = np.atleast_2d(pts)
    d = np.zeros((pts.shape[0], 3, 3))
    for m in range(3):
        d[:, m, m] = -2.0 * pts[:, m]
    return d


def strain_matrix(dNx: np.ndarray) -> np.ndarray:
    """Strain-displacement matrices B (..., 6, 3n) from physical derivatives dNx (..., n, 3)."""
    shp = dNx.shape[:-2]
    n = dNx.shape[-2]
    B = np.zeros(shp + (6, 3 * n), dtype=dNx.dtype)
    dx, dy, dz = dNx[..., 0], dNx[..., 1], dNx[..., 2]
    B[..., 0, 0::3] = dx
    B[..., 1, 1::3] = dy
    B[..., 2, 2::3] = dz
    B[..., 3, 0::3] = dy
    B[..., 3, 1::3] = dx
    B[..., 4, 0::3] = dz
    B[..., 4, 2::3] = dx
    B[..., 5, 1::3] = dz
    B[..., 5, 2::3] = dy
    return B


def _jacobian(dN: np.ndarray, xyz: np.ndarray):
    """J[e, g, i, j] = d x_j / d xi_i, its determinant and inverse."""
    J = np.einsum("gai,eaj->egij", dN, xyz)
    det = np.linalg.det(J)
    return J, det


def _check_det(det: np.ndarray, xyz: np.ndarray, what: str) -> None:
    span = np.linalg.norm(xyz.max(axis=1) - xyz.min(axis=1), axis=1)      # (nE,)
    tol = 1e-12 * (span / 2.0) ** 3
    bad = np.where((det <= tol[:, None]).any(axis=1))[0]
    if bad.size:
        raise ElementError(f"EDU-05: non-positive Jacobian at {what} Gauss points of SOLID element "
                           f"(batch index {int(bad[0])}); check the node order (1-2-3-4 face, 5-8 opposite)")


def solid_kernel(xyz: np.ndarray, lam: np.ndarray, G: np.ndarray, rho: np.ndarray, eint: int = 0,
                 incompatible: bool = False, mass: str = "mixed", want_mass: bool = True,
                 want_recovery: bool = True) -> Dict[str, np.ndarray]:
    """Vectorised SOLID matrices for nE elements (see :func:`_solid_kernel` for the arguments)."""
    with blas_quiet():
        out = _solid_kernel(xyz, lam, G, rho, eint, incompatible, mass, want_mass, want_recovery)
    for k, v in out.items():
        check_finite(v, f"SOLID {k}")
    return out


def _solid_kernel(xyz, lam, G, rho, eint, incompatible, mass, want_mass, want_recovery):
    """Vectorised SOLID matrices for nE elements.

    Parameters
    ----------
    xyz : (nE, 8, 3) node coordinates (repeated rows for degenerate shapes)
    lam, G : (nE,) complex Lame constants; rho : (nE,) mass densities
    eint : 0/1/2 -> 2/3/4 Gauss points per direction (D-ELM-01)
    incompatible : add the 9 Wilson-Taylor modes (structural elements only, D-ELM-02)

    Returns dict with K (nE,24,24) complex, M (nE,24,24) real, S (nE,7,24) complex stress
    recovery at the centroid (SOCT row zero) and B (nE,6,24) strain recovery.
    """
    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim == 2:
        xyz = xyz[None]
    nE = xyz.shape[0]
    if xyz.shape[1:] != (8, 3):
        raise ElementError("SOLID needs 8 nodes with 3 coordinates (repeat nodes for prisms/pyramids)")
    lam = np.broadcast_to(np.asarray(lam, dtype=complex), (nE,))
    G = np.broadcast_to(np.asarray(G, dtype=complex), (nE,))
    rho = np.broadcast_to(np.asarray(rho, dtype=float), (nE,))
    try:
        ng = EINT_POINTS[int(eint)]
    except KeyError:
        raise ElementError(f"SOLID EINT must be 0, 1 or 2 (got {eint})") from None

    # ---- Gauss points, Jacobians and compatible strain matrices ----------------------------
    pts, wts = gauss_3d(ng)
    _, dN = hex8_shape(pts)
    J, det = _jacobian(dN, xyz)
    _check_det(det, xyz, "stiffness")
    Jinv = np.linalg.inv(J)
    dNx = np.einsum("egjk,gak->egaj", Jinv, dN)                 # dN/dx = J^-1 dN/dxi
    Bu = strain_matrix(dNx)                                      # (nE, ng, 6, 24)

    # centroid quantities (recovery, Taylor correction)
    _, dN0 = hex8_shape(np.zeros((1, 3)))
    J0, det0 = _jacobian(dN0, xyz)                               # (nE, 1, 3, 3)
    if np.any(det0[:, 0] <= 0):
        raise ElementError("EDU-05: non-positive Jacobian at the centroid of a SOLID element")
    J0inv = np.linalg.inv(J0[:, 0])                              # (nE, 3, 3)
    Bu0 = strain_matrix(np.einsum("ejk,ak->eaj", J0inv, dN0[0]))  # (nE, 6, 24)

    if incompatible:
        # Wilson modes with the Taylor centroid-Jacobian correction: dNa/dx = (detJ0/detJ) J0^-1 dNa/dxi
        dNa = _incompatible_dshape(pts)                          # (ng, 3, 3)
        scale = det0[:, 0][:, None] / det                        # (nE, ng)
        dNax = np.einsum("ejk,gak->egaj", J0inv, dNa) * scale[:, :, None, None]
        Ba = strain_matrix(dNax)                                 # (nE, ng, 6, 9)
        Bf = np.concatenate([Bu, Ba], axis=3)                    # (nE, ng, 6, 33)
        Ba0 = strain_matrix(np.einsum("ejk,ak->eaj", J0inv, _incompatible_dshape(np.zeros((1, 3)))[0]))
    else:
        Bf = Bu
    wdet = det * wts[None, :]                                    # (nE, ng)
    div = Bf[:, :, 0, :] + Bf[:, :, 1, :] + Bf[:, :, 2, :]       # m^T B (volumetric strain row)
    K_lam = np.einsum("eg,ega,egb->eab", wdet, div, div)
    K_G = np.einsum("eg,egia,i,egib->eab", wdet, Bf, np.diag(_D_G), Bf)
    Kf = lam[:, None, None] * K_lam + G[:, None, None] * K_G     # complex (nE, nf, nf)

    nu = 24
    if incompatible:
        Kuu, Kua = Kf[:, :nu, :nu], Kf[:, :nu, nu:]
        Kau, Kaa = Kf[:, nu:, :nu], Kf[:, nu:, nu:]
        X = np.linalg.solve(Kaa, Kau)                            # alpha = -X u
        K = Kuu - Kua @ X
        Brec = Bu0.astype(complex) - Ba0 @ X                     # strain of u + condensed modes
    else:
        K = Kf
        Brec = Bu0.astype(complex)
    K = 0.5 * (K + np.swapaxes(K, 1, 2))                         # remove round-off asymmetry
    out: Dict[str, np.ndarray] = {"K": K}

    if want_mass:
        ngm = max(ng, 3)                                         # exact for trilinear geometry
        ptm, wtm = gauss_3d(ngm)
        Nm, dNm = hex8_shape(ptm)
        _, detm = _jacobian(dNm, xyz)
        _check_det(detm, xyz, "mass")
        Mc8 = np.einsum("eg,ga,gb->eab", detm * wtm[None, :] * rho[:, None], Nm, Nm)
        Mc = np.einsum("eab,ij->eaibj", Mc8, np.eye(3)).reshape(nE, 24, 24)
        out["M"] = mixed_mass(Mc, mass)
    if want_recovery:
        D = lam[:, None, None] * _D_LAM[None] + G[:, None, None] * _D_G[None]
        S = np.zeros((nE, 7, 24), dtype=complex)
        S[:, :6, :] = D @ Brec
        out["S"] = S
        out["B"] = Brec
    return out


# ------------------------------------------------------------------------------------------
# Registry API (ARCHITECTURE 6.2)
# ------------------------------------------------------------------------------------------
def _mat_arrays(mats: Sequence[MaterialProps]):
    lam = np.array([m.lam for m in mats], dtype=complex)
    G = np.array([m.G for m in mats], dtype=complex)
    rho = np.array([m.rho for m in mats], dtype=float)
    return lam, G, rho


def matrices(xyz, mat: MaterialProps, incompatible: bool = False, eint: int = 0, mass: str = "mixed",
             **_ignored):
    """(K* (24,24) complex, M (24,24) real) of one SOLID element."""
    r = solid_kernel(np.asarray(xyz, float)[None], mat.lam, mat.G, mat.rho, eint=eint,
                     incompatible=incompatible, mass=mass, want_recovery=False)
    return r["K"][0], r["M"][0]


def recovery(xyz, mat: MaterialProps, incompatible: bool = False, eint: int = 0, quantity: str = "stress",
             **_ignored):
    """Centroid recovery operator: stresses (7, 24) [SXX SYY SZZ SXY SXZ SYZ, SOCT row = 0] or,
    with ``quantity='strain'``, engineering strains (6, 24) [EXX EYY EZZ GXY GXZ GYZ]."""
    r = solid_kernel(np.asarray(xyz, float)[None], mat.lam, mat.G, mat.rho, eint=eint,
                     incompatible=incompatible, want_mass=False)
    return r["B"][0] if quantity == "strain" else r["S"][0]


def batch(xyz, mats: Sequence[MaterialProps], incompatible: bool = False, eint: int = 0, mass: str = "mixed",
          **_ignored):
    """dict(K (nE,24,24), M (nE,24,24), S (nE,7,24), B (nE,6,24)) for elements sharing
    eint/incompatible/mass."""
    lam, G, rho = _mat_arrays(mats)
    return solid_kernel(xyz, lam, G, rho, eint=eint, incompatible=incompatible, mass=mass)


SPEC = register(ElementSpec(code=CODE, name=NAME, nnodes=NNODES, dofs=DOFS, components=list(COMPONENTS),
                            matrices=matrices, recovery=recovery, batch=batch,
                            description="8-node hexahedron, optional 9 incompatible modes, 1/2+1/2 mass"))
