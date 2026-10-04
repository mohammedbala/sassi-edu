"""2D (plane-strain) SSI helpers: HOUSE ``<dim>`` = 1, POINT2, PLANE elements (requirements 1.5 PLANE,
4.3 POINT2, 4.6 item 2, D-PNT-01, D-W2-07; R1 2.6, 3.5, 4.1).

How a 2D SSI analysis works in SASSI-EDU (for the structural engineer)
---------------------------------------------------------------------
The model lies in the global X-Z plane (Z up; the y coordinates are ignored, manual 9.3.1) and is a
slice of unit thickness of an infinitely long structure and site: PLANE elements are plane-strain
elements of unit thickness, masses, stiffnesses and loads are per unit length out of plane.

* SITE: the free field and the Rayleigh/Love modes are the same for 2D and 3D (one FILE1/FILE2
  serves both, requirements 4.2).  The in-plane input is SV (x'), P (z') or Rayleigh waves; SH and
  Love waves are anti-plane (y') and are not analysed by the in-plane 2D model (below).
* POINT2 (``<dim>`` = 1, D-PNT-01): unit *line* loads on a plane-strain strip ``|x| <= R0`` with
  Waas-Lysmer transmitting boundaries (:mod:`sassi.core.strip2d`); the far field is the exact modal
  expansion ``sum_j alpha_j phi_j exp(-i k_j (|x| - R0))`` (R1 3.5).
* HOUSE: the interaction DOFs are the in-plane translations UX, UZ (D-W2-07, FILE4 ``int_dofs`` =
  [1, 3]); BEAMS / SPRING / GENERAL nodes keep their six DOFs (the out-of-plane DOFs UY, ROTX, ROTZ
  are structural only -- fix them unless an out-of-plane response is wanted).
* ANALYS: the soil flexibility is the 2x2 P-SV block per node pair with the odd coupling terms times
  ``sgn(x_i - x_j)`` (requirements 4.6 item 2); the impedance ``X_ff`` is the inverse of the
  in-plane block (:func:`sassi.core.symmetry.keep_mask`), ``2 nInt x 2 nInt``, stored in the 3D
  layout with zero UY rows and columns.  The free field ``U'_f`` is evaluated at ``(x, yc, z)``; the
  coordinate transformation angle must be 0 or 180 deg (the plane of propagation is the X-Z plane).
  The SSI equation, the Schur solver, restarts and the FILE8 layout are those of 3D.
* MOTION, RELDISP, STRESS (PLANE: SXX SZZ TXZ at the centroid) work on FILE8 and FILE4 unchanged.

Not available in 2D (documented restrictions): anti-plane SH / Love input with PLANE elements (they
have no UY DOF), incoherent motion (manual 6.5.4, G-17), the global impedance (requirements 4.6 item 9:
3D only; :func:`rigid_strip_impedance` gives the 2D rigid-strip impedance for verification).
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .tlm import quiet_fpe

#: interaction DOFs of a 2D model (D-W2-07)
IN_PLANE_DOFS = (1, 3)
#: rigid-body DOFs of a 2D foundation: X translation, Z translation, rotation about Y
RIGID_DOFS_2D = ("X", "Z", "YY")
#: allowed ANALYS coordinate transformation angles in 2D (x' along +X or -X)
ANGLES_2D = (0.0, 180.0)


def dimension_name(dim: int) -> str:
    return "2D (plane strain, POINT2)" if int(dim) == 1 else "3D (POINT3)"


def point_dimension_error(dim4: int, dim3: int) -> str:
    """D-PNT-01: FILE3 must come from POINT2 for a 2D HOUSE model and from POINT3 for a 3D one;
    returns the error text ('' when consistent)."""
    if int(dim4) == int(dim3):
        return ""
    return (f"FILE3 holds {'POINT2 (2D line-load)' if int(dim3) == 1 else 'POINT3 (3D point-load)'} solutions but "
            f"the HOUSE model is {dimension_name(dim4)} (HOUSE <dim> = {dim4}; D-PNT-01) -- re-run POINT with "
            f"<dim> = {dim4}")


def inplane_input_error(cm: int, ang: float, name: str = "FILE1") -> str:
    """The seismic input of a 2D (in-plane) model: control direction x' or z' and the angle 0 or 180 deg
    (x' in the X-Z plane); returns the error text ('' when admissible)."""
    if int(cm) == 1:
        return (f"{name}: control motion direction y' (SH / Love waves) is anti-plane: the in-plane 2D model "
                "(PLANE elements, interaction DOFs UX, UZ) cannot be excited by it -- use SV/P/Rayleigh "
                "input (x' or z') or a 3D model")
    a = float(ang) % 360.0
    if min(abs(a - v) for v in ANGLES_2D) > 1e-9:
        return (f"coordinate transformation angle {ang:g} deg: a 2D model lies in the X-Z plane, so x' must be "
                "along +X (0 deg) or -X (180 deg)")
    return ""


def plane_coordinates(xyz: np.ndarray, yc: float) -> np.ndarray:
    """Node coordinates of a 2D model for the free field: the y coordinate (ignored in 2D, manual
    9.3.1) is replaced by the control-point ``yc`` so that it never enters the phase ``exp(-i k x')``."""
    out = np.array(xyz, dtype=float, copy=True).reshape(-1, 3)
    out[:, 1] = float(yc)
    return out


def duplicate_key(xyz: np.ndarray, dim: int) -> np.ndarray:
    """Rounded coordinates that identify coincident interaction nodes (x, z in 2D; x, y, z in 3D)."""
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    cols = [0, 2] if int(dim) == 1 else [0, 1, 2]
    P = xyz[:, cols]
    scale = max(1e-9, float(np.abs(P).max()) * 1e-9) if P.size else 1e-9
    return np.round(P / scale).astype(np.int64)


def rigid_transform_2d(xyz: np.ndarray, ref) -> np.ndarray:
    """``T`` (3n, 3): node translations ``[ux, uy, uz]`` of a rigid-body motion ``(u0, w0, theta_y)``
    of a 2D foundation about ``ref = (xc, zc)`` (rotation about +Y: ``ux = u0 + theta (z - zc)``,
    ``uz = w0 - theta (x - xc)``); the UY rows are zero."""
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    n = xyz.shape[0]
    xc, zc = float(ref[0]), float(ref[1])
    T = np.zeros((n, 3, 3))
    T[:, 0, 0] = 1.0
    T[:, 2, 1] = 1.0
    T[:, 0, 2] = xyz[:, 2] - zc
    T[:, 2, 2] = -(xyz[:, 0] - xc)
    return T.reshape(3 * n, 3)


@quiet_fpe
def rigid_strip_impedance(X: np.ndarray, xyz: np.ndarray, ref) -> np.ndarray:
    """2D rigid-foundation impedance ``K = T^T X_ff T`` (3x3, DOFs X, Z, YY; per unit length) of the
    interaction nodes ``xyz`` -- the plane-strain counterpart of the global impedance of requirements
    4.6 item 9 (used by VP-42 as an independent check of the FORCE + ANALYS compliance)."""
    T = rigid_transform_2d(xyz, ref)
    return T.T @ np.asarray(X, complex) @ T


def strip_stiffness_reference(B: float, H: float, G: float, nu: float = 0.30) -> Tuple[float, float]:
    """Static stiffness per unit length of a rigid strip of half-width ``B`` on a layer of thickness
    ``H`` over a rigid base (Jakub & Roesset 1977, MIT R77-36, Tables 1/2; R2 B.4), valid
    1/8 <= B/H <= 1/2: ``K_x = G 1.175 (1 + 2.15 B/H)``, ``K_phi = G B^2 2.394 (1 + 0.17 B/H)`` for
    nu = 0.30 (the Poisson-specific rows of R2 B.4 for nu = 0, 0.15, 0.30, 0.45)."""
    table = {0.0: (0.904, 2.5, 1.891, 0.17), 0.15: (1.017, 2.37, 2.094, 0.17),
             0.30: (1.175, 2.15, 2.394, 0.17), 0.45: (1.419, 1.95, 2.907, 0.24)}
    key = min(table, key=lambda v: abs(v - nu))
    if abs(key - nu) > 1e-12:
        raise ValueError(f"no strip reference for nu = {nu} (tabulated: {sorted(table)})")
    ax, bx, ap, bp = table[key]
    r = B / H
    return G * ax * (1.0 + bx * r), G * B * B * ap * (1.0 + bp * r)


def stiffness_from_compliance(C: np.ndarray) -> np.ndarray:
    """Stiffness matrix from a (complex) compliance matrix of a vibration analysis (unit loads)."""
    return np.linalg.inv(np.asarray(C, complex))


def in_plane(eq_dof: np.ndarray) -> np.ndarray:
    """Mask of the in-plane DOFs (UX, UZ, ROTY) of a 2D FILE4/FILE8 equation list."""
    d = np.asarray(eq_dof)
    return np.isin(d, (1, 3, 5))


def optional_int(meta: dict, key: str, default: int = 0) -> int:
    v: Optional[object] = meta.get(key, default)
    try:
        return int(v)        # type: ignore[arg-type]
    except (TypeError, ValueError):
        return default
