"""TSHELL (GROUP type 5): flat-facet thick shell, Mindlin-Reissner plate + membrane + drilling.

Normative sources: requirements 1.5 (TSHELL row: 4 nodes, 6 DOF, EINT 0 reduced / 1 selective,
automatic small drilling stiffness, lumped mass), 4.1 (TSHELL row: MITC4; EINT 0 reduced -- 1-point
bending and shear, hourglass-stabilised; EINT 1 selective -- 2x2 bending, 1-point shear; drilling
stiffness 1e-4 x min membrane diagonal; THSHLSTR face stresses), decision D-ELM-08, D-CNV-04
(E* = E c(beta_s), nu real), D-CNV-05 (lumped mass), D-STR-10 (N, Q per unit length), D-TSH-01
(face stresses), spec 08 section 4.6 (nodes I J K L counter-clockwise, triangles, local axes), spec
05d (STRESS components NXX NYY NXY QXZ QYZ MXX MYY MXY), manual sections 6 (HOUSE remarks 1-2) and
9.4.9 (EINT).

What a structural engineer should know
--------------------------------------
The thin SHELL element (Kirchhoff) neglects transverse shear deformation.  For walls and slabs whose
thickness is not small compared with the span or the wavelength of the response (thick basemats,
shear walls of nuclear buildings, the higher modes of any plate) the Mindlin-Reissner theory is the
better model: the normal to the mid-surface stays straight but not normal, so the rotations
``beta`` of the normal are independent of the slope of ``w`` and the transverse shear strain
``gamma = grad w + beta`` carries the shear forces ``Q = kappa G t gamma`` (kappa = 5/6).

Its well-known numerical difficulty is **shear locking**: with low-order interpolation a thin plate
cannot make ``gamma = 0`` everywhere and the element becomes far too stiff (a cantilever plate 1000
times longer than thick, 8 x 2 elements, deflects about 6500 times too little with plain 2x2 integration
of the shear, VP-TS1).  The cure used here is
the MITC4 *assumed natural strain* field of Bathe & Dvorkin (1985): the transverse shear strains are
sampled at the four mid-side tying points and interpolated from there, which removes locking (the
thin-plate limit converges to the Kirchhoff solution, VP-TS1) and introduces no spurious mode.

Local axes, DOFs and kinematics (spec 08 4.6, as SHELL)
--------------------------------------------------------
``x' = unit(m_JK - m_LI)``, ``z' = unit(x' x (m_KL - m_IJ))``, ``y' = z' x x'`` (m_AB mid-point of
side AB; triangles use L = K).  Warped quadrilaterals are projected on the mean plane through the
centroid.  A triangle is given with 3 nodes or with 4 nodes of which one pair of adjacent corners
coincides (:func:`sassi.elements.shell.corner_rows`).  Per node the local DOFs are
``u', v', w', theta_x', theta_y', theta_z'`` (right-hand rotations) and::

    u(z) = u' + z beta_x,  v(z) = v' + z beta_y,  beta_x = theta_y',  beta_y = -theta_x'

    membrane  eps   = [u,x, v,y, u,y + v,x]            N = D_m eps,   D_m = E t/(1-nu^2) [1 nu 0; nu 1 0; 0 0 (1-nu)/2]
    bending   kappa = [beta_x,x, beta_y,y, beta_x,y + beta_y,x]   M = D_b kappa,  D_b = E t^3/(12(1-nu^2)) [...]
    shear     gamma = [w,x + beta_x, w,y + beta_y]      Q = D_s gamma,  D_s = kappa G t I,  kappa = 5/6

so ``M_xx = integral(sigma_xx z dz)`` is positive with tension on the +z' face and the plate
equilibrium reads ``Q_x = M_xx,x + M_xy,y`` (no distributed couples).  In the Kirchhoff limit
``gamma -> 0`` gives ``beta = -grad w``, the SHELL convention.

Quadrilateral (MITC4)
---------------------
* **Membrane** (both EINT values; the manual reduces only bending and shear): the SHELL membrane, a
  plane-stress Q4 with the 4 Wilson-Taylor incompatible modes condensed (2x2 Gauss), which is exact in
  in-plane bending of rectangles and passes the patch test (:func:`sassi.elements.shell.membrane_q4`).
* **Transverse shear** (both EINT values): MITC4.  The covariant shear strains
  ``g_xi = w,xi + beta . x,xi`` and ``g_eta = w,eta + beta . x,eta`` of the displacement field are
  tied at the edge mid-points A (0, 1), C (0, -1) and B (1, 0), D (-1, 0) and interpolated::

      g~_xi  = (1 + eta)/2 g_xi(A) + (1 - eta)/2 g_xi(C)
      g~_eta = (1 + xi)/2 g_eta(B) + (1 - xi)/2 g_eta(D)          gamma = J^-1 [g~_xi, g~_eta]

  integrated with 2x2 Gauss points.  This is how the manual's "1-point transverse shear" is realised
  for both EINT options: at the centre ``g~`` is the 1-point (tying-average) value, and the linear
  terms act as an *assumed-strain hourglass stabilisation* of that one-point shear (the one-point
  shell elements of Belytschko and co-workers stabilise their shear with the same Bathe-Dvorkin
  field).  A literal 1-point shear would leave the ``w = +-1`` hourglass mode without stiffness (a
  spurious mode that pollutes SSI responses); the MITC4 field has none and does not lock.
* **Bending**, EINT 1 *selective*: 2x2 Gauss integration (exact for the bilinear rotations).
* **Bending**, EINT 0 *reduced* (the manual's default): one-point integration at the centre plus a
  physical hourglass stabilisation::

      K_b = A B_0^T D_b B_0 + eps_hg sum_g w_g |J_g| (B_g - B_0)^T D_b (B_g - B_0)

  with ``eps_hg`` = :data:`HOURGLASS` (0.1, the customary stiffness hourglass coefficient of one-point
  shell elements).  ``B_g - B_0`` vanishes for every linear rotation field, so constant curvature
  states are exact for any quadrilateral (patch test) and the two rotation hourglass modes
  (``beta = xi eta`` patterns, which neither the one-point curvature nor the MITC4 shear sees) get
  ``eps_hg`` times their exact (2x2) stiffness.  ``eps_hg = 1`` reproduces EINT 1 on parallelograms.

Triangle (MITC3 shear + condensed rotation bubble + stabilised linear shear part)
-------------------------------------------------------------------------------
CST membrane.  Plate (:func:`tri_plate_operators`; EINT has no effect, every integral is exact):

* ``w`` linear; rotations linear in the corner values **plus a cubic bubble** ``f4 beta_b``,
  ``f4 = 27 r s (1 - r - s)`` (r, s area coordinates), whose two DOFs are condensed statically inside
  the element (no load acts on them: SASSI has no distributed couples);
* transverse shear: the MITC3 field of the corner DOFs (Lee & Bathe 2004; the rotated lowest-order
  Raviart-Thomas interpolant) ``gamma~ = a + c (-(y - y_c), x - x_c)``, whose tangential component along
  each side equals the side average of the displacement-based strain
  ``(w_j - w_i + (beta_i + beta_j).d/2) / |d|``; the bubble vanishes on the sides, so it enters the
  shear through its element mean: ``gamma_0 = a + (27/60) beta_b``;
* shear energy ``A gamma_0^T D_s gamma_0 + w_lin int (c rot)^T D_s (c rot) dA`` with the
  Lyly-Stenberg-Vihinen (1993) weight ``w_lin = t^2/(t^2 + alpha h^2)`` (alpha =
  :data:`TRI_STAB_ALPHA` = 0.2, h the longest side) on the linear ("curl") part only.

Why (VP-TS1 section 4): the plain MITC3 triangle has as many shear constraints as DOFs in the thin limit
(one per side) and locks -- a simply supported plate with t/a = 1/1000 deflects 46 % of Navier on a
16 x 16 mesh with all diagonals the same way.  The bubble lets every element satisfy the constant-shear
constraint at the cost of its bending energy (a discrete-Kirchhoff-like coupling of w and the
rotations, as the bubble of the MITC3+ element of Lee, Lee & Bathe 2014 or of Arnold & Falk 1989), and
the weight relaxes the curl constraint only where it would lock (``h >> t``); for ``h << t`` the
element tends to MITC3.  Unlike the MITC3+ tying (interior points where f4 = 1/2, curl part scaled by
d = 1e-4) this keeps the full curl stiffness of thick elements (no near-zero-energy element mode) and
the shear sees the bubble with exactly the weight of its work (mean 27/60), so the constant-shear
equilibrium patch test with distributed couples is exact.  Constant curvature: the bubble decouples
(``int grad f4 dA = 0``) and the MITC3 shear of a Kirchhoff state is zero, so the bending patch test is
exact; with the bubble at its exact value (0 for linear rotations) the strain operators reproduce
constant shear.  At the centroid ``grad f4 = 0``: the bubble does not change the recovered moments.

Drilling stiffness (manual HOUSE remark 1, D-ELM-08)
-----------------------------------------------------
"A small rotational stiffness is automatically added in HOUSE": the nodal drilling rotations
``theta_z'`` are tied to the in-plane rotation of the membrane at the centre,
``omega = (v,x - u,y)/2``, by the penalty ``k_d sum_i (theta_z',i - omega)^2`` with
``k_d = DRILL_FACTOR x min(diag K_m) x A`` (:data:`DRILL_FACTOR` = 1e-4; the membrane diagonal (F/L)
times the element area gives a rotational stiffness F.L, so the factor is unit independent).  The
penalty vanishes for a rigid in-plane rotation (``theta_z = omega``), so the free element keeps
its six rigid-body modes, and it is 1e-4 of the membrane stiffness, so it does not alter the
membrane response; it only makes FIXROT / FIXSHLROT unnecessary for TSHELL models.

Mass (D-CNV-05 "lumped", spec 08 4.6)
-------------------------------------
``rho t A / n`` on each translation and the Mindlin rotary inertia ``rho t^3/12 A / n`` on the two
bending rotations ``theta_x', theta_y'`` (spec 08 4.6, inferred option; required to reproduce the
Mindlin thick-plate frequencies of NAFEMS Test 21, which include rotary inertia); no inertia on the
drilling rotation (its stiffness is a numerical device, an inertia would create spurious
low-frequency drilling modes).  In global axes the rotational block is
``rho t^3/12 A/n (I - z' z'^T)``.

Damping (D-CNV-04): ``E* = E0 c(beta_s)`` with real nu, so ``K* = c(beta_s) K0`` and the recovery
operator is scaled by the same factor.

Recovery (ARCHITECTURE 6.2, spec 05d, D-STR-10): at the centre, local axes x'y'z'::

    NXX NYY NXY  membrane forces per unit length  N = D_m eps      (F/L)
    QXZ QYZ      transverse shear forces per unit length Q = D_s gamma~   (F/L; MITC field)
    MXX MYY MXY  moments per unit length          M = D_b kappa    (F.L/L)

THSHLSTR face stresses (D-TSH-01): :func:`face_stresses`.
"""
from __future__ import annotations

import re
from typing import Dict, List, Sequence, Tuple

import numpy as np

from ..conventions import ELEMENT_COMPONENTS
from .base import ElementError, ElementSpec, MaterialProps, blas_quiet, check_finite, gauss_legendre, register
from .plane import quad4_shape
from .shell import corner_rows, local_frame, membrane_cst, membrane_q4, polygon_area

__all__ = ["NAME", "CODE", "DOFS", "NNODES", "COMPONENTS", "SHEAR_FACTOR", "DRILL_FACTOR", "HOURGLASS",
           "EINT_REDUCED", "EINT_SELECTIVE", "TRI_STAB_ALPHA", "TRI_BUBBLE_MEAN", "matrices", "recovery",
           "local_matrices", "face_stresses", "FACE_PERMUTATIONS", "FACE_COLUMNS", "plate_rigidities",
           "mitc4_shear_operator", "tri_plate_operators", "SPEC"]

NAME = "TSHELL"
CODE = 5
DOFS = (1, 2, 3, 4, 5, 6)
NNODES = 4
COMPONENTS = list(ELEMENT_COMPONENTS[NAME])        # NXX NYY NXY QXZ QYZ MXX MYY MXY

#: transverse shear correction factor (Reissner 5/6; Mindlin's pi^2/12 differs by 1.3 %)
SHEAR_FACTOR = 5.0 / 6.0
#: drilling stiffness k_d = DRILL_FACTOR x min diag(K_membrane) x area (D-ELM-08)
DRILL_FACTOR = 1.0e-4
#: EINT 0 bending hourglass stabilisation coefficient (fraction of the exact hourglass stiffness)
HOURGLASS = 0.1
#: EINT values (manual 9.4.9): 0 reduced (default), 1 selective
EINT_REDUCED, EINT_SELECTIVE = 0, 1

#: natural coordinates of the MITC4 tying points: A (0, 1), C (0, -1) for g_xi; B (1, 0), D (-1, 0) for g_eta
_TIE_XI = np.array([[0.0, 1.0], [0.0, -1.0]])
_TIE_ETA = np.array([[1.0, 0.0], [-1.0, 0.0]])


# ======================================================================================
# material rigidities
# ======================================================================================
def _plate_matrix(nu: float) -> np.ndarray:
    return np.array([[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, 0.5 * (1.0 - nu)]])


def plate_rigidities(E: float, nu: float, t: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(D_m (3,3), D_b (3,3), D_s (2,2)): membrane, bending and transverse shear rigidities of a
    homogeneous plate of thickness ``t`` (per unit length)."""
    P = _plate_matrix(nu)
    G = E / (2.0 * (1.0 + nu))
    return (E * t / (1.0 - nu * nu) * P, E * t ** 3 / (12.0 * (1.0 - nu * nu)) * P,
            SHEAR_FACTOR * G * t * np.eye(2))


# ======================================================================================
# quadrilateral operators (bending DOFs per node: w, theta_x, theta_y)
# ======================================================================================
def _jacobian(xy: np.ndarray, p: np.ndarray):
    """(N (4,), dN natural (4, 2), J (2, 2) with rows x,xi and x,eta) at natural point p."""
    N, dN = quad4_shape(np.asarray(p, float)[None])
    J = dN[0].T @ xy
    return N[0], dN[0], J


def _bending_B(dNx: np.ndarray) -> np.ndarray:
    """Curvature operator (3, 12): kappa = [beta_x,x, beta_y,y, beta_x,y + beta_y,x] with
    beta_x = theta_y, beta_y = -theta_x (columns per node: w, theta_x, theta_y)."""
    B = np.zeros((3, 12))
    B[0, 2::3] = dNx[:, 0]
    B[1, 1::3] = -dNx[:, 1]
    B[2, 1::3] = -dNx[:, 0]
    B[2, 2::3] = dNx[:, 1]
    return B


def _covariant_shear(N: np.ndarray, dN: np.ndarray, J: np.ndarray) -> np.ndarray:
    """Displacement-based covariant shear strains (2, 12): rows g_xi = w,xi + beta . x,xi and
    g_eta = w,eta + beta . x,eta (beta = (theta_y, -theta_x))."""
    B = np.zeros((2, 12))
    for r in range(2):                       # r = 0: xi, 1: eta ; J[r] = (x,r, y,r)
        B[r, 0::3] = dN[:, r]
        B[r, 1::3] = -N * J[r, 1]
        B[r, 2::3] = N * J[r, 0]
    return B


def mitc4_shear_operator(xy: np.ndarray):
    """Function ``p -> (B_s (2, 12), det J)`` of the MITC4 assumed transverse shear strains
    ``gamma = [gamma_xz, gamma_yz]`` at natural point p (Bathe & Dvorkin 1985)."""
    tie = {}
    for name, pts, row in (("A", _TIE_XI[0], 0), ("C", _TIE_XI[1], 0), ("B", _TIE_ETA[0], 1),
                           ("D", _TIE_ETA[1], 1)):
        N, dN, J = _jacobian(xy, pts)
        tie[name] = _covariant_shear(N, dN, J)[row]

    def op(p):
        xi, eta = float(p[0]), float(p[1])
        _, _, J = _jacobian(xy, p)
        cov = np.vstack([0.5 * (1 + eta) * tie["A"] + 0.5 * (1 - eta) * tie["C"],
                         0.5 * (1 + xi) * tie["B"] + 0.5 * (1 - xi) * tie["D"]])
        det = float(np.linalg.det(J))
        return np.linalg.solve(J, cov), det
    return op


def _quad_plate(xy: np.ndarray, Db: np.ndarray, Ds: np.ndarray, eint: int, hourglass: float):
    """Bending + shear stiffness (12, 12) and centre operators (kappa (3, 12), gamma (2, 12))."""
    x, w = gauss_legendre(2)
    pts = [np.array([a, b]) for a in x for b in x]
    wts = [wa * wb for wa in w for wb in w]
    shear = mitc4_shear_operator(xy)
    _, dN0, J0 = _jacobian(xy, np.zeros(2))
    det0 = float(np.linalg.det(J0))
    B0 = _bending_B(np.linalg.solve(J0, dN0.T).T)
    Kb = np.zeros((12, 12))
    Khg = np.zeros((12, 12))
    Ks = np.zeros((12, 12))
    for p, wg in zip(pts, wts):
        _, dN, J = _jacobian(xy, p)
        det = float(np.linalg.det(J))
        if det <= 0.0:
            raise ElementError("EDU-05: non-positive Jacobian in a TSHELL element (re-entrant or "
                               "clockwise quadrilateral?)")
        Bg = _bending_B(np.linalg.solve(J, dN.T).T)
        Kb += wg * det * Bg.T @ Db @ Bg
        dB = Bg - B0
        Khg += wg * det * dB.T @ Db @ dB
        Bs, _ = shear(p)
        Ks += wg * det * Bs.T @ Ds @ Bs
    if int(eint) == EINT_REDUCED:
        Kb = 4.0 * det0 * B0.T @ Db @ B0 + float(hourglass) * Khg      # area = 4 det J(0) (bilinear map)
    Bs0, _ = shear(np.zeros(2))
    return Kb + Ks, B0, Bs0


# ======================================================================================
# triangle operators (MITC3 shear + condensed rotation bubble + stabilised linear shear part)
# ======================================================================================
#: Lyly-Stenberg-Vihinen weight t^2/(t^2 + alpha h^2) of the linear part of the triangle shear (h = longest side)
TRI_STAB_ALPHA = 0.2
#: element mean of the cubic bubble f4 = 27 r s (1 - r - s) (integral / area = 27/60)
TRI_BUBBLE_MEAN = 27.0 / 60.0
#: DOFs of the uncondensed triangle plate: (w, theta_x, theta_y) per corner, then the bubble rotation (beta_x, beta_y)
TRI_NODAL, TRI_BUBBLE = slice(0, 9), slice(9, 11)
#: 7-point degree-5 rule on the triangle (area coordinates r, s; weights sum to 1): exact for the bubble bending
_T7_A, _T7_B, _T7_C, _T7_D = 0.0597158717897698, 0.4701420641051151, 0.7974269853530873, 0.1012865073234563
_TRI7_RS = np.array([[1.0 / 3.0, 1.0 / 3.0], [_T7_B, _T7_B], [_T7_A, _T7_B], [_T7_B, _T7_A],
                     [_T7_D, _T7_D], [_T7_C, _T7_D], [_T7_D, _T7_C]])
_TRI7_W = np.array([0.225] + [0.1323941527885062] * 3 + [0.1259391805448271] * 3)


def _tri_derivatives(xy: np.ndarray) -> Tuple[np.ndarray, float]:
    """Cartesian derivatives (3, 2) of the linear shape functions and the area."""
    A = polygon_area(xy)
    if A <= 0.0:
        raise ElementError("degenerate TSHELL triangle (zero or negative area)")
    x, y = xy[:, 0], xy[:, 1]
    b = np.array([y[1] - y[2], y[2] - y[0], y[0] - y[1]]) / (2.0 * A)
    c = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]]) / (2.0 * A)
    return np.column_stack([b, c]), A


def _mitc3_coefficients(xy: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """``(T (3, 9), centroid)``: coefficients ``[a_x, a_y, c] = T u_b`` of the MITC3 field
    ``gamma~(x) = a + c (-(y - y_c), x - x_c)`` from the corner bending DOFs u_b (w, theta_x, theta_y per
    node) with linear rotations (Lee & Bathe 2004; the rotated lowest-order Raviart-Thomas interpolant).

    Each side k = (i, j), d = x_j - x_i: the tangential integral of the field equals that of the
    displacement-based strain, ``gamma~(m_k) . d = (w_j - w_i) + (beta_i + beta_j) . d / 2``."""
    xc = xy.mean(axis=0)
    Q = np.zeros((3, 3))
    R = np.zeros((3, 9))
    for k, (i, j) in enumerate(((0, 1), (1, 2), (2, 0))):
        d = xy[j] - xy[i]
        m = 0.5 * (xy[i] + xy[j]) - xc
        Q[k] = [d[0], d[1], m[0] * d[1] - m[1] * d[0]]
        R[k, 3 * j] += 1.0
        R[k, 3 * i] -= 1.0
        for n in (i, j):
            R[k, 3 * n + 2] += 0.5 * d[0]               # beta_x = theta_y
            R[k, 3 * n + 1] -= 0.5 * d[1]               # beta_y = -theta_x
    return np.linalg.solve(Q, R), xc


def _bubble_gradient(xy: np.ndarray, r: float, s: float) -> np.ndarray:
    """Cartesian gradient (2,) of the cubic bubble f4 = 27 r s (1 - r - s) at area coordinates (r, s)."""
    J = np.array([xy[1] - xy[0], xy[2] - xy[0]])          # rows x,r and x,s
    return np.linalg.solve(J, [27.0 * s * (1.0 - 2.0 * r - s), 27.0 * r * (1.0 - r - 2.0 * s)])


def tri_plate_operators(xy: np.ndarray, Db: np.ndarray, Ds: np.ndarray, t: float,
                        alpha: float = TRI_STAB_ALPHA) -> Dict[str, np.ndarray]:
    """Uncondensed plate operators of a TSHELL triangle with corners ``xy`` (3, 2) (counter-clockwise, local
    axes); 11 DOFs: ``(w, theta_x, theta_y)`` per corner, then the bubble rotation ``(beta_x, beta_y)``.

    Fields::

        w = sum h_i w_i                        beta = sum h_i beta_i + f4 beta_b     (f4 = 27 r s (1 - r - s))
        gamma~ = a + TRI_BUBBLE_MEAN beta_b + c (-(y - y_c), x - x_c)              (a, c: MITC3 tying)

    ``a`` and ``c`` are the MITC3 coefficients of the corner DOFs (:func:`_mitc3_coefficients`); the
    bubble, which vanishes on the sides and so is invisible to the side tying, enters the shear through its
    element mean.  Energy: bending ``int kappa^T D_b kappa`` (7-point rule, exact), shear
    ``A gamma_0^T D_s gamma_0`` of the constant part ``gamma_0 = a + 0.45 beta_b`` plus ``w_lin int
    (c rot)^T D_s (c rot)`` of the linear part with the weight ``w_lin = t^2 / (t^2 + alpha h^2)``.

    Returns ``K`` (11, 11), the centroid operators ``Bk`` (3, 11) (curvatures) and ``Bs`` (2, 11) (shear
    strains, the constant part), ``c`` (11,) and ``weight`` (w_lin)."""
    xy = np.asarray(xy, dtype=float)
    dNx, A = _tri_derivatives(xy)
    T, xc = _mitc3_coefficients(xy)
    Bk = np.zeros((3, 11))                                # curvatures of the linear corner rotations
    Bk[0, 2:9:3], Bk[1, 1:9:3] = dNx[:, 0], -dNx[:, 1]
    Bk[2, 1:9:3], Bk[2, 2:9:3] = -dNx[:, 0], dNx[:, 1]
    Bs = np.zeros((2, 11))
    Bs[:, TRI_NODAL] = T[:2]
    Bs[0, 9] = Bs[1, 10] = TRI_BUBBLE_MEAN
    c = np.zeros(11)
    c[TRI_NODAL] = T[2]
    h = max(float(np.linalg.norm(xy[j] - xy[i])) for i, j in ((0, 1), (1, 2), (2, 0)))
    weight = t * t / (t * t + float(alpha) * h * h)
    d = xy - xc                                           # second moments of area about the centroid
    Ixx, Iyy, Ixy = (A / 12.0 * float(np.sum(u * v)) for u, v in ((d[:, 0], d[:, 0]), (d[:, 1], d[:, 1]),
                                                                 (d[:, 0], d[:, 1])))
    rot_energy = Ds[0, 0] * Iyy - (Ds[0, 1] + Ds[1, 0]) * Ixy + Ds[1, 1] * Ixx   # int rot^T D_s rot dA
    K = A * (Bk.T @ Db @ Bk + Bs.T @ Ds @ Bs) + weight * rot_energy * np.outer(c, c)
    for (r, s), wg in zip(_TRI7_RS, _TRI7_W):              # bubble bending (int grad f4 = 0: no coupling)
        g = _bubble_gradient(xy, r, s)
        Bb = np.array([[g[0], 0.0], [0.0, g[1]], [g[1], g[0]]])
        K[TRI_BUBBLE, TRI_BUBBLE] += wg * A * Bb.T @ Db @ Bb
    return {"K": 0.5 * (K + K.T), "Bk": Bk, "Bs": Bs, "c": c, "weight": np.array(weight)}


def _tri_plate(xy: np.ndarray, Db: np.ndarray, Ds: np.ndarray, t: float):
    """Condensed bending + shear stiffness (9, 9) and centroid operators (kappa (3, 9), gamma (2, 9)) of a
    triangle: the bubble rotation is eliminated by static condensation, ``beta_b = -K_bb^-1 K_bn u``
    (no load acts on it in SASSI: there are no distributed couples)."""
    op = tri_plate_operators(xy, Db, Ds, t)
    K, n, b = op["K"], TRI_NODAL, TRI_BUBBLE
    X = np.linalg.solve(K[b, b], K[b, n])
    Kc = K[n, n] - K[n, b] @ X
    return 0.5 * (Kc + Kc.T), op["Bk"][:, n], op["Bs"][:, n] - op["Bs"][:, b] @ X


# ======================================================================================
# element (local axes)
# ======================================================================================
def _drilling(dNx0: np.ndarray, n: int, kd: float) -> np.ndarray:
    """Drilling penalty (3n, 3n) on [u', v', theta_z'] per node: kd sum_i (theta_z,i - omega)^2,
    omega = (v,x - u,y)/2 at the centre (derivatives ``dNx0`` (n, 2))."""
    b = np.zeros(3 * n)
    b[0::3] = -0.5 * dNx0[:, 1]
    b[1::3] = 0.5 * dNx0[:, 0]
    K = np.zeros((3 * n, 3 * n))
    for i in range(n):
        r = -b.copy()
        r[3 * i + 2] += 1.0
        K += kd * np.outer(r, r)
    return K


def local_matrices(xy: np.ndarray, E: float, nu: float, t: float, eint: int = EINT_REDUCED,
                   hourglass: float = HOURGLASS) -> Tuple[np.ndarray, np.ndarray]:
    """Real local stiffness (6n, 6n) and recovery operator (8, 6n) of a flat element with corner
    coordinates ``xy`` (n = 3 or 4, counter-clockwise, local axes), Young's modulus ``E``, Poisson's
    ratio ``nu`` and thickness ``t``.  DOFs per node: u', v', w', theta_x', theta_y', theta_z'."""
    n = xy.shape[0]
    if n not in (3, 4):
        raise ElementError("TSHELL needs 3 or 4 distinct corners")
    Dm, Db, Ds = plate_rigidities(E, nu, t)
    with blas_quiet():
        if n == 4:
            Km, Bm = membrane_q4(xy, Dm, incompatible=True)
            Kp, Bk, Bg = _quad_plate(xy, Db, Ds, eint, hourglass)
            _, dN0, J0 = _jacobian(xy, np.zeros(2))
            dNx0 = np.linalg.solve(J0, dN0.T).T
            area = 4.0 * float(np.linalg.det(J0))
        else:
            Km, Bm = membrane_cst(xy, Dm)
            Kp, Bk, Bg = _tri_plate(xy, Db, Ds, t)
            dNx0, area = _tri_derivatives(xy)
    mem = np.array([[6 * a, 6 * a + 1] for a in range(n)]).ravel()
    ben = np.array([[6 * a + 2, 6 * a + 3, 6 * a + 4] for a in range(n)]).ravel()
    drl = np.array([[6 * a, 6 * a + 1, 6 * a + 5] for a in range(n)]).ravel()
    K = np.zeros((6 * n, 6 * n))
    K[np.ix_(mem, mem)] += Km
    K[np.ix_(ben, ben)] += Kp
    kd = DRILL_FACTOR * float(np.min(np.diag(Km))) * area
    K[np.ix_(drl, drl)] += _drilling(dNx0, n, kd)
    R = np.zeros((8, 6 * n))
    R[0:3, mem] = Dm @ Bm                 # N  (F/L)
    R[3:5, ben] = Ds @ Bg                 # Q  (F/L)
    R[5:8, ben] = Db @ Bk                 # M  (F.L/L)
    return 0.5 * (K + K.T), R


def _thickness(thick) -> float:
    if thick is None or not np.isfinite(float(thick)) or float(thick) <= 0:
        raise ElementError("TSHELL thickness must be > 0 (THICK)")
    return float(thick)


def _eint(eint) -> int:
    e = 0 if eint is None else int(eint)
    if e not in (EINT_REDUCED, EINT_SELECTIVE):
        raise ElementError(f"TSHELL EINT must be 0 (reduced) or 1 (selective), got {eint}")
    return e


def _element(xyz, mat: MaterialProps, thick, eint, hourglass):
    X = np.asarray(xyz, dtype=float)
    t = _thickness(thick)
    e = _eint(eint)
    try:
        keep = corner_rows(X)
        Lam, _, xy = local_frame(X)
        K0, R0 = local_matrices(xy, mat.E0, mat.nu0, t, e, hourglass)
    except ElementError as exc:                       # helpers shared with SHELL: name the element type
        raise ElementError(re.sub(r"\bSHELL\b", "TSHELL", str(exc))) from None
    return X, keep, Lam, xy, t, K0, R0


def _expand(a: np.ndarray, keep: List[int], n_rows: int, square: bool) -> np.ndarray:
    """Matrices built on the distinct corners ``keep`` -> layout of the ``n_rows`` input rows (zero
    rows/columns for a repeated corner)."""
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


def matrices(xyz, mat: MaterialProps, thick: float = None, eint: int = EINT_REDUCED,
             hourglass: float = HOURGLASS, **_ignored):
    """(K* complex, M real) in global axes, size 6 x (rows of ``xyz``).

    ``eint`` 0 reduced (default, manual 9.4.9) or 1 selective; ``hourglass`` the EINT 0 bending
    stabilisation coefficient.  A triangle given with 4 rows (one adjacent corner repeated) gets
    zero rows/columns for the repeated corner."""
    X, keep, Lam, xy, t, K0, _ = _element(xyz, mat, thick, eint, hourglass)
    n = len(keep)
    T = np.kron(np.eye(2 * n), Lam)
    cE = mat.shell_moduli()[0] / mat.E0                          # c(beta_s): E* = E0 c(beta_s)
    with blas_quiet():
        K = cE * (T.T @ K0 @ T)
    A = polygon_area(xy)
    m_tr = mat.rho * t * A / n
    m_rot = mat.rho * t ** 3 / 12.0 * A / n
    P = np.eye(3) - np.outer(Lam[2], Lam[2])                     # bending rotations only (no drilling)
    M = np.zeros((6 * n, 6 * n))
    for a in range(n):
        M[6 * a:6 * a + 3, 6 * a:6 * a + 3] = m_tr * np.eye(3)
        M[6 * a + 3:6 * a + 6, 6 * a + 3:6 * a + 6] = m_rot * P
    K, M = _expand(K, keep, X.shape[0], True), _expand(M, keep, X.shape[0], True)
    return check_finite(0.5 * (K + K.T), "TSHELL K"), 0.5 * (M + M.T)


def recovery(xyz, mat: MaterialProps, thick: float = None, eint: int = EINT_REDUCED,
             hourglass: float = HOURGLASS, **_ignored):
    """Centre resultants (8, 6 x rows) in local axes: NXX NYY NXY QXZ QYZ (F/L), MXX MYY MXY (F.L/L)."""
    X, keep, Lam, _, _, _, R0 = _element(xyz, mat, thick, eint, hourglass)
    n = len(keep)
    T = np.kron(np.eye(2 * n), Lam)
    cE = mat.shell_moduli()[0] / mat.E0
    with blas_quiet():
        S = cE * (R0 @ T)
    return check_finite(_expand(S, keep, X.shape[0], False), "TSHELL recovery")


# ======================================================================================
# THSHLSTR face stresses and strains (D-TSH-01, manual section 6 STRESS, spec 11 section 4.1)
# ======================================================================================
#: sign permutations (s_N, s_M) of the maximum force and moment contributions (manual: ++, --, +-, -+)
FACE_PERMUTATIONS = (("++", 1.0, 1.0), ("--", -1.0, -1.0), ("+-", 1.0, -1.0), ("-+", -1.0, 1.0))
#: columns of one face-result row (:func:`face_stresses`)
FACE_COLUMNS = ("SXX", "SYY", "TXY", "S1", "S2", "EXX", "EYY", "GXY", "E1", "E2")


def face_stresses(maxima: Sequence[float], t: float, E: float, nu: float) -> Dict[str, object]:
    """Face (top/bottom) stresses and strains of a TSHELL from the **maximum values** of its eight basic
    components ``(NXX NYY NXY QXZ QYZ MXX MYY MXY)`` (THSHLSTR,1; D-TSH-01).

    Plate formulas (local axes; the top face is z' = +t/2)::

        sigma_face = N/t + s 6 M/t^2         (s = +1 top face, -1 bottom face)

    The maxima of N and M occur at different times with unknown signs, so the manual combines them
    with the four sign permutations (s_N, s_M) = ++, --, +-, -+ for (NXX, MXX) and (NYY, MYY), and with
    ``++`` only for the in-plane shear (NXY, MXY), "assuming that the in-plane shear stress is the
    largest".  Since ``(s_N, s_M)`` on the top face equals ``(s_N, -s_M)`` on the bottom face, the four
    permutations cover both faces.  For each permutation::

        SXX = s_N NXX/t + s_M 6 MXX/t^2,  SYY = s_N NYY/t + s_M 6 MYY/t^2,  TXY = NXY/t + 6 MXY/t^2
        S1,S2 = (SXX + SYY)/2 +- sqrt(((SXX - SYY)/2)^2 + TXY^2)              (principal stresses)
        EXX = (SXX - nu SYY)/E,  EYY = (SYY - nu SXX)/E,  GXY = TXY/G,  G = E/(2(1+nu))  (plane stress)
        E1,E2 = (EXX + EYY)/2 +- sqrt(((EXX - EYY)/2)^2 + (GXY/2)^2)         (principal strains)

    The transverse shear stresses are largest at the mid-surface (parabolic distribution, zero on the
    faces): ``TXZ = 1.5 QXZ/t``, ``TYZ = 1.5 QYZ/t``.

    Returns ``{'rows': {perm: (10,) array in FACE_COLUMNS order}, 'TXZ': float, 'TYZ': float}``.
    """
    v = np.abs(np.asarray(maxima, dtype=float).ravel())
    if v.size != 8:
        raise ElementError("face_stresses needs the 8 TSHELL components NXX NYY NXY QXZ QYZ MXX MYY MXY")
    t = _thickness(t)
    if not E > 0:
        raise ElementError("face strains need E > 0")
    Nxx, Nyy, Nxy, Qxz, Qyz, Mxx, Myy, Mxy = v
    G = E / (2.0 * (1.0 + nu))
    rows = {}
    for name, sn, sm in FACE_PERMUTATIONS:
        sxx = sn * Nxx / t + sm * 6.0 * Mxx / t ** 2
        syy = sn * Nyy / t + sm * 6.0 * Myy / t ** 2
        txy = Nxy / t + 6.0 * Mxy / t ** 2
        c, r = 0.5 * (sxx + syy), np.hypot(0.5 * (sxx - syy), txy)
        exx, eyy, gxy = (sxx - nu * syy) / E, (syy - nu * sxx) / E, txy / G
        ce, re = 0.5 * (exx + eyy), np.hypot(0.5 * (exx - eyy), 0.5 * gxy)
        rows[name] = np.array([sxx, syy, txy, c + r, c - r, exx, eyy, gxy, ce + re, ce - re])
    return {"rows": rows, "TXZ": 1.5 * Qxz / t, "TYZ": 1.5 * Qyz / t}


SPEC = register(ElementSpec(code=CODE, name=NAME, nnodes=NNODES, dofs=DOFS, components=list(COMPONENTS),
                            matrices=matrices, recovery=recovery,
                            description="flat Mindlin-Reissner shell: Q4+incompatible/CST membrane, MITC4 shear "
                                        "(quads) / MITC3 + condensed rotation bubble + stabilised linear shear "
                                        "(triangles), EINT 0 reduced (1-point + hourglass) / 1 selective (2x2) "
                                        "bending, drilling penalty, lumped mass with rotary inertia"))
