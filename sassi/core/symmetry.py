"""SYMM: symmetry / antisymmetry planes of half and quarter models (requirements 3.4.F SYMM, D-ANL-12,
spec 07 section 9.2.39, manual sections 2.7 "Symmetry of the System" and 9.2.39).

What a symmetry plane means (for the structural engineer)
----------------------------------------------------------
A structure that is symmetric about a vertical plane, on a horizontally layered site, under a loading
that is *symmetric* or *antisymmetric* with respect to that plane, has a response with the same
property.  Only one half (one plane) or one quarter (two orthogonal planes) of the model is then
analysed, with the boundary conditions of the symmetry on the nodes of the plane ("with respect to the
loading", manual 9.2.39).  For a plane with unit normal ``n`` (along X or Y; the planes must be
parallel to the XZ or YZ planes, in 2D the line must be parallel to Z) and the reflection
``P = I - 2 n n^T`` (``diag(-1, 1, 1)`` for a plane normal to X):

* symmetric loading (``SYMM`` type 0):      ``u(x') =  P u(x)``  -> on the plane the normal translation
  and the two in-plane rotations are zero;
* antisymmetric loading (type 1):           ``u(x') = -P u(x)``  -> on the plane the two in-plane
  translations and the normal rotation are zero.

``x'`` is the mirror image of ``x``.  Rotations are pseudo-vectors (``theta(x') = -s P theta(x)`` with
``s = +1`` / ``-1`` for type 0 / 1), which gives the rotation conditions above (D-ANL-12).  A
horizontal X input (vertical SV with x' = X) is antisymmetric about a plane normal to X and symmetric
about a plane normal to Y; a vertical (P) input is symmetric about both.

The soil: image interaction nodes (D-ANL-12)
-------------------------------------------
HOUSE fixes the DOFs above on the nodes of the plane; ANALYS needs the soil flexibility of the *reduced*
interaction set.  The interaction forces of the full model have the symmetry of the loading:
``f(j') = s P f(j)`` for the image ``j'`` of node ``j``.  The free-field displacement at an interaction
node ``i`` of the reduced model is therefore (unit-load flexibility ``F(i, j)`` of POINT, R1 4.1)::

    (u - u')(i) = sum_j [ F(i, j) g(j) + F(i, j') s P g(j) ]  =  sum_j F_red(i, j) g(j)

    F_red(i, j) = F(i, j) + s F(i, j') P                       one plane
    F_red(i, j) = sum over the 2^k subsets S of the k planes of  s_S F(i, j_S) P_S     (k = 1, 2)

with ``j_S`` the image of ``j`` in the planes of ``S`` (``j_{}`` = j), ``s_S`` and ``P_S`` the products
of the signs and reflections of those planes.  (D-ANL-12 writes the image term as ``P F(i, j') P``;
for a layered site ``F(i, j') P = P F(i', j)``, so the reflection belongs on the load side as written
here.)  ``g`` is the soil force on the *reduced* structure: a node on a plane carries one half (one
quarter on two planes) of the full-model force, exactly as the structure, the excavated soil, the
masses and the loads of a half model carry half of the plane node's full-model values -- the image
sum gives the plane node's own term ``F(i, p)(I + s P)``, i.e. twice the force of the half model.
For the DOFs that the symmetry constrains on the plane, ``(I + s P)`` is zero: the rows and columns of
``F_red`` vanish (zero force, zero displacement) and they are removed *before* the inversion
``X_red = F_red^-1`` (:func:`reduced_impedance`).  With reciprocity and the reflection invariance of
the layered site ``F_red`` is symmetric; it is symmetrised like ``F_ff`` (D-ANL-02).

Modelling rules of a half / quarter model (manual 9.3.2: "constraining appropriate degrees of freedom
on the planes of symmetry"): the model lies on one side of every plane; elements and masses lying *in*
a plane (e.g. a stick on the axis) carry 1/2 (1/4 on two planes) of their full-model values; loads on
plane nodes likewise.  Incoherent motion, wave passage, multiple excitation and the global impedance
are not allowed with SYMM (manual 2.7 and 6.5.4, G-17).
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

SYMMETRIC, ANTISYMMETRIC = 0, 1
TYPE_NAMES = {SYMMETRIC: "symmetry", ANTISYMMETRIC: "antisymmetry"}
AXIS_NAMES = ("X", "Y")
DOF_NAMES = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ")
#: at most two planes (manual 9.2.39)
MAX_PLANES = 2

#: manual Chapter 10 errors of the SYMM command (spec 11 section 10.2)
ERROR_TEXT = {2: "Error 2: Illegal Node for Symmetry Plane/Line {i}",
              3: "Error 3: Illegal Number of Nodes for Symmetry Plane/Line {i}",
              4: "Error 4: Nodes from Symmetry Line {i} Are Equal",
              5: "Error 5: Nodes from Symmetry Plane {i} Are Collinear"}


@dataclass(frozen=True)
class SymmetryPlane:
    """One SYMM plane (3D, parallel to XZ or YZ) or line (2D, parallel to Z).

    ``axis`` is the direction of the plane normal: 0 = X (the plane ``x = coord``, parallel to YZ),
    1 = Y (the plane ``y = coord``, parallel to XZ).  ``type`` 0 symmetric, 1 antisymmetric loading."""

    no: int
    type: int
    axis: int
    coord: float
    nodes: Tuple[int, ...] = ()

    @property
    def sign(self) -> float:
        """``s = +1`` (symmetric) or ``-1`` (antisymmetric): ``u(image) = s P u``."""
        return 1.0 if self.type == SYMMETRIC else -1.0

    @property
    def reflection(self) -> np.ndarray:
        """Diagonal of the reflection ``P`` of translations: -1 on the normal component."""
        p = np.ones(3)
        p[self.axis] = -1.0
        return p

    @property
    def fixed_dofs(self) -> Tuple[int, ...]:
        """DOF labels (1..6) fixed on the nodes of the plane (D-ANL-12)."""
        a = self.axis
        b = 1 - a                                 # the other horizontal axis
        if self.type == SYMMETRIC:                # normal translation, in-plane rotations (about b and z)
            return tuple(sorted((a + 1, b + 4, 6)))
        return tuple(sorted((b + 1, 3, a + 4)))   # in-plane translations, normal rotation

    def mirror(self, xyz: np.ndarray) -> np.ndarray:
        """Mirror images of points ``xyz`` (n, 3) (or (n, 2) plan coordinates)."""
        out = np.array(xyz, dtype=float, copy=True)
        out[..., self.axis] = 2.0 * self.coord - out[..., self.axis]
        return out

    def distance(self, xyz: np.ndarray) -> np.ndarray:
        """Signed distance of points from the plane along its normal."""
        return np.asarray(xyz, float)[..., self.axis] - self.coord

    def describe(self, dim: int = 2) -> str:
        what = "line" if dim == 1 else "plane"
        return (f"{what} {self.no}: {AXIS_NAMES[self.axis].lower()} = {self.coord:.6g} "
                f"({TYPE_NAMES[self.type]}, type {self.type})")


# ======================================================================================
# Plane definitions (HOUSE deck table 'symm')
# ======================================================================================
def parse_planes(rows: Iterable[dict], node_xyz: Dict[int, Sequence[float]], dim: int, tol: float
                 ) -> Tuple[List[SymmetryPlane], List[str], List[str]]:
    """Planes of the HOUSE deck rows ``(no, type, n1, n2, n3)`` with the manual's Errors 2-5 and the
    geometry rules of 9.2.39 (EDU-26): at most two planes, 3D planes parallel to XZ or YZ, a 2D line
    parallel to Z, two planes orthogonal (parallel symmetry planes would need infinitely many images).

    ``tol`` is the geometric tolerance of the model (D-GEN-06).  Returns ``(planes, errors, warnings)``;
    rows with ``n1 = 0`` are resets (ignored with a warning)."""
    planes: List[SymmetryPlane] = []
    errors: List[str] = []
    warnings: List[str] = []
    seen = set()
    what = "line" if int(dim) == 1 else "plane"
    for r in rows:
        no, typ = int(r["no"]), int(r["type"])
        ids = [int(r[k]) for k in ("n1", "n2", "n3")]
        if ids[0] == 0:
            warnings.append(f"SYMM {what} {no}: <node1> = 0 resets the {what} (row ignored)")
            continue
        if no in seen:
            errors.append(f"EDU-26: SYMM {what} {no} is defined twice")
            continue
        seen.add(no)
        if no not in (1, 2):
            errors.append(f"EDU-26: SYMM {what} number {no} must be 1 or 2 (at most two planes, manual 9.2.39)")
        if typ not in (SYMMETRIC, ANTISYMMETRIC):
            errors.append(f"EDU-26: SYMM {what} {no}: <type> = {typ} must be 0 (symmetry) or 1 (antisymmetry)")
            continue
        if any(n < 0 or (n and n not in node_xyz) for n in ids):
            errors.append(ERROR_TEXT[2].format(i=no) + f" (nodes {ids})")
            continue
        given = [n for n in ids if n]
        if len(given) < 2:
            errors.append(ERROR_TEXT[3].format(i=no))
            continue
        P = np.array([node_xyz[n] for n in given], dtype=float).reshape(-1, 3)
        plane = _plane_from_points(no, typ, tuple(given), P, dim, tol, errors)
        if plane is not None:
            planes.append(plane)
    if len(planes) > MAX_PLANES:
        errors.append(f"EDU-26: {len(planes)} SYMM planes (at most {MAX_PLANES}, manual 9.2.39)")
    if int(dim) == 1 and len(planes) > 1:
        errors.append("EDU-26: a 2D model has at most one SYMM line (two lines parallel to Z would need "
                      "infinitely many images)")
    elif len(planes) == 2 and planes[0].axis == planes[1].axis:
        errors.append(f"EDU-26: SYMM planes {planes[0].no} and {planes[1].no} are parallel "
                      f"({AXIS_NAMES[planes[0].axis].lower()} = const): two symmetry planes must be orthogonal")
    return sorted(planes, key=lambda p: p.no), errors, warnings


def _plane_from_points(no: int, typ: int, ids: Tuple[int, ...], P: np.ndarray, dim: int, tol: float,
                       errors: List[str]) -> Optional[SymmetryPlane]:
    spread = np.ptp(P, axis=0)
    if int(dim) == 1:                                    # 2D: line parallel to Z through 2 (or 3) nodes
        if len(set(ids)) < len(ids) or (spread[0] <= tol and spread[2] <= tol):
            errors.append(ERROR_TEXT[4].format(i=no))
            return None
        if spread[0] > tol:
            errors.append(f"EDU-26: SYMM line {no} (nodes {list(ids)}) is not parallel to the Z axis "
                          "(2D models: the line of symmetry must be parallel to Z, manual 9.2.39)")
            return None
        return SymmetryPlane(no, typ, 0, float(P[:, 0].mean()), ids)
    if len(ids) == 3:
        a, b = P[1] - P[0], P[2] - P[0]
        if len(set(ids)) < 3 or np.linalg.norm(np.cross(a, b)) <= 1e-9 * np.linalg.norm(a) * max(np.linalg.norm(b), 1e-300):
            errors.append(ERROR_TEXT[5].format(i=no))
            return None
    elif len(set(ids)) < 2 or float(np.linalg.norm(P[1] - P[0])) <= tol:
        errors.append(ERROR_TEXT[4].format(i=no))
        return None
    flat_x, flat_y = spread[0] <= tol, spread[1] <= tol
    if flat_x and not flat_y:
        return SymmetryPlane(no, typ, 0, float(P[:, 0].mean()), ids)
    if flat_y and not flat_x:
        return SymmetryPlane(no, typ, 1, float(P[:, 1].mean()), ids)
    if flat_x and flat_y:
        errors.append(f"EDU-26: SYMM plane {no}: the nodes {list(ids)} lie on one vertical line and do not define "
                      "the plane -- give a node off that line")
    else:
        errors.append(f"EDU-26: SYMM plane {no} (nodes {list(ids)}) is not parallel to the XZ or YZ plane "
                      "(manual 9.2.39: 3D symmetry planes are vertical and parallel to XZ or YZ)")
    return None


def on_planes(planes: Sequence[SymmetryPlane], xyz: np.ndarray, tol: float) -> np.ndarray:
    """(n, nP) boolean: point k lies on plane p (within ``tol``)."""
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    if not planes:
        return np.zeros((xyz.shape[0], 0), bool)
    return np.stack([np.abs(p.distance(xyz)) <= tol for p in planes], axis=1)


def symmetry_fixity(planes: Sequence[SymmetryPlane], xyz: np.ndarray, tol: float) -> np.ndarray:
    """(n, 6) boolean: DOFs fixed on the nodes of the planes (D-ANL-12)."""
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    fix = np.zeros((xyz.shape[0], 6), bool)
    on = on_planes(planes, xyz, tol)
    for k, p in enumerate(planes):
        for d in p.fixed_dofs:
            fix[on[:, k], d - 1] = True
    return fix


def side_violations(planes: Sequence[SymmetryPlane], xyz: np.ndarray, tol: float) -> List[Tuple[int, int, int]]:
    """For every plane: ``(plane no, points on the negative side, points on the positive side)`` when the
    points lie on *both* sides (a half / quarter model must lie on one side of each plane)."""
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    out = []
    for p in planes:
        d = p.distance(xyz)
        neg, pos = int(np.sum(d < -tol)), int(np.sum(d > tol))
        if neg and pos:
            out.append((p.no, neg, pos))
    return out


# ======================================================================================
# Image sums (ANALYS)
# ======================================================================================
def image_terms(planes: Sequence[SymmetryPlane]) -> List[Tuple[Tuple[SymmetryPlane, ...], float, np.ndarray]]:
    """All ``2^k`` subsets ``S`` of the planes with ``(S, s_S, diag(P_S))``; the empty subset first."""
    out = []
    for k in range(len(planes) + 1):
        for S in combinations(planes, k):
            s = 1.0
            P = np.ones(3)
            for p in S:
                s *= p.sign
                P = P * p.reflection
            out.append((tuple(S), s, P))
    return out


def mirror_points(xyz: np.ndarray, S: Sequence[SymmetryPlane]) -> np.ndarray:
    """Images of ``xyz`` in the planes of ``S`` (the reflections of orthogonal planes commute)."""
    out = np.array(xyz, dtype=float, copy=True)
    for p in S:
        out = p.mirror(out)
    return out


def reduced_flexibility(block: Callable[[np.ndarray], np.ndarray], xyz: np.ndarray,
                        planes: Sequence[SymmetryPlane]) -> np.ndarray:
    """Flexibility of the reduced interaction set by image superposition (module docstring)::

        F_red = sum_S  s_S  F(obs = nodes, load = images_S(nodes))  P_S

    ``block(xyz_load)`` must return a *new* array: the free-field flexibility (3n, 3n), node-major
    ``[ux, uy, uz]``, of the observation points ``xyz`` (the reduced interaction nodes, fixed) under
    unit loads at the points ``xyz_load`` (the nodes or their images, same interfaces); it is scaled
    in place.  Not symmetrised."""
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    n = xyz.shape[0]
    F = None
    for S, s, P in image_terms(planes):                   # accumulated in place: at most two dense arrays
        B = np.asarray(block(mirror_points(xyz, S)), dtype=complex).reshape(n, 3, n, 3)
        if S:
            B *= (s * P)[None, None, None, :]             # columns: reflected load components
        if F is None:
            F = B
        else:
            F += B
        del B
    return np.zeros((0, 0), complex) if F is None else F.reshape(3 * n, 3 * n)


def constrained_translations(planes: Sequence[SymmetryPlane], xyz: np.ndarray, tol: float) -> np.ndarray:
    """(n, 3) boolean: translations of interaction nodes constrained by the symmetry (zero force and
    zero displacement: removed from F_red before the inversion)."""
    return symmetry_fixity(planes, xyz, tol)[:, :3]


def keep_mask(dim: int, constrained: Optional[np.ndarray], n: int) -> np.ndarray:
    """(3n,) boolean: interaction DOFs of the soil impedance (node-major x, y, z): in 2D the in-plane
    UX, UZ only (D-W2-07), minus the translations constrained by SYMM planes."""
    keep = np.ones((n, 3), bool)
    if int(dim) == 1:
        keep[:, 1] = False
    if constrained is not None and np.size(constrained):
        keep &= ~np.asarray(constrained, bool).reshape(n, 3)
    return keep.reshape(-1)


def reduced_impedance(F: np.ndarray, keep: np.ndarray) -> Tuple[np.ndarray, float]:
    """Soil impedance on the active interaction DOFs (requirements 4.6 items 2-3, D-ANL-02/04/12):
    ``F`` (3n x 3n, the free-field or image-reduced flexibility) is symmetrised in place, the DOFs
    with ``keep = False`` (2D: UY; SYMM: translations constrained on a plane) are removed, the rest
    is inverted and scattered back into a 3n x 3n ``X`` with zero rows and columns at the removed
    DOFs (the 3D layout of COOX / FILE8 is kept).  Returns ``(X, rcond(F_kept))``."""
    from .ssi_solver import impedance_matrix, symmetrize_inplace
    F = np.asarray(F, complex)
    n3 = F.shape[0]
    symmetrize_inplace(F)
    keep = np.asarray(keep, bool).reshape(-1)
    if keep.size != n3:
        raise ValueError(f"keep mask of length {keep.size} for a {n3} x {n3} flexibility")
    if keep.all():
        return impedance_matrix(F, overwrite=True)
    idx = np.flatnonzero(keep)
    X = np.zeros((n3, n3), complex)
    if idx.size == 0:
        return X, 1.0
    Fk = np.asfortranarray(F[np.ix_(idx, idx)])
    Xk, rc = impedance_matrix(Fk, overwrite=True)
    X[np.ix_(idx, idx)] = Xk
    return X, rc


def field_symmetry_error(planes: Sequence[SymmetryPlane], U: np.ndarray, U_images: Sequence[np.ndarray],
                         components: Sequence[int] = (0, 1, 2)) -> List[Tuple[float, float]]:
    """For every plane: the relative deviations of a field from symmetry and from antisymmetry,
    ``max |U(image) - s P U| / max |U|`` with s = +1 and s = -1 (``U_images[p]`` = the field at the
    images of the points in plane p).  The seismic free field must have the symmetry of the SYMM type
    of every plane (ANALYS check)."""
    U = np.asarray(U, complex).reshape(-1, 3)
    c = list(components)
    scale = max(float(np.abs(U[:, c]).max()) if U.size else 0.0, 1e-300)
    out = []
    for p, Ui in zip(planes, U_images):
        Ui = np.asarray(Ui, complex).reshape(-1, 3)
        PU = U * p.reflection[None, :]
        e_sym = float(np.abs(Ui[:, c] - PU[:, c]).max()) / scale if U.size else 0.0
        e_anti = float(np.abs(Ui[:, c] + PU[:, c]).max()) / scale if U.size else 0.0
        out.append((e_sym, e_anti))
    return out


def planes_array(planes: Sequence[SymmetryPlane]) -> np.ndarray:
    """FILE4 array ``x_symm`` (nP, 4): ``no, type, axis, coord`` of every plane."""
    return np.array([[p.no, p.type, p.axis, p.coord] for p in planes], dtype=float).reshape(-1, 4)


def planes_from_array(a) -> List[SymmetryPlane]:
    """Inverse of :func:`planes_array` (FILE4 ``x_symm``)."""
    a = np.asarray(a, float).reshape(-1, 4)
    return [SymmetryPlane(int(r[0]), int(r[1]), int(r[2]), float(r[3])) for r in a]


def required_type(axis: int, cm: int) -> int:
    """SYMM type consistent with a vertically incident input in global direction ``cm`` (0 X, 1 Y, 2 Z)
    for a plane normal to ``axis``: the input is antisymmetric about the plane normal to its own
    horizontal direction, symmetric about the other one and symmetric for vertical input."""
    return ANTISYMMETRIC if int(cm) == int(axis) else SYMMETRIC
