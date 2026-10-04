"""Coordinate-system geometry of the model database (spec 08 section 2, decision D-MDL-04).

Conventions
-----------
* The global system is right-handed Cartesian with **Z up**.
* A local Cartesian system ``s`` is stored as its origin ``O`` (global coordinates) and a
  rotation matrix ``R`` whose **columns are the local x, y, z unit vectors expressed in global
  axes**.  Then::

      X_global = O + R @ x_local          x_local = R.T @ (X_global - O)

* ``LOC,<ns>,0,<x0>,<y0>,<z0>,<txy>,<tyz>,<txz>`` (D-MDL-04, ANSYS ``LOCAL`` convention) applies
  three successive rotations about the *moving* axes, angles in degrees:

  1. ``txy`` about Z (turns X toward Y),
  2. ``tyz`` about the new X (turns Y toward Z),
  3. ``txz`` about the new Y (turns Z toward X),

  so that ``R = Rz(txy) @ Rx(tyz) @ Ry(txz)``.
* ``LOCAL,<ns>,0,<n1>,<n2>,<n3>``: origin at n1, x along n1->n2, z normal to the plane
  (n1, n2, n3) and ``y = z x x`` so that n3 lies at local y > 0 (spec 08 Q7: the manual's
  "X x Z" would be left-handed).
"""
from __future__ import annotations

import math
from typing import Sequence, Tuple

import numpy as np


def cosd(a: float) -> float:
    """cos of an angle in degrees, exact (0, +-1) at multiples of 90 degrees."""
    r = math.fmod(a, 360.0)
    if r % 90.0 == 0.0:
        q = int(round(r / 90.0)) % 4
        return (1.0, 0.0, -1.0, 0.0)[q]
    return math.cos(math.radians(a))


def sind(a: float) -> float:
    """sin of an angle in degrees, exact (0, +-1) at multiples of 90 degrees."""
    r = math.fmod(a, 360.0)
    if r % 90.0 == 0.0:
        q = int(round(r / 90.0)) % 4
        return (0.0, 1.0, 0.0, -1.0)[q]
    return math.sin(math.radians(a))


def rot_z(a: float) -> np.ndarray:
    c, s = cosd(a), sind(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rot_x(a: float) -> np.ndarray:
    c, s = cosd(a), sind(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def rot_y(a: float) -> np.ndarray:
    c, s = cosd(a), sind(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def loc_matrix(txy: float, tyz: float, txz: float) -> np.ndarray:
    """Rotation matrix of ``LOC``: ``R = Rz(txy) Rx(tyz) Ry(txz)`` (columns = local axes)."""
    return rot_z(txy) @ rot_x(tyz) @ rot_y(txz)


def local_from_points(p1: Sequence[float], p2: Sequence[float], p3: Sequence[float]) -> np.ndarray:
    """Rotation matrix of ``LOCAL`` from three global points (spec 08 section 3.9).

    ``ex = (p2-p1)/|p2-p1|``, ``ez = (a x b)/|a x b|`` with ``a = p2-p1``, ``b = p3-p1``, and
    ``ey = ez x ex``.  Raises ``ValueError`` when the points are coincident or collinear.
    """
    p1 = np.asarray(p1, float)
    a = np.asarray(p2, float) - p1
    b = np.asarray(p3, float) - p1
    la = float(np.linalg.norm(a))
    if la == 0.0:
        raise ValueError("nodes n1 and n2 coincide")
    n = np.cross(a, b)
    ln = float(np.linalg.norm(n))
    if ln <= 1e-12 * la * max(float(np.linalg.norm(b)), 1e-300):
        raise ValueError("nodes n1, n2, n3 are collinear or coincident")
    ex = a / la
    ez = n / ln
    ey = np.cross(ez, ex)
    return np.column_stack([ex, ey, ez])


def euler_from_matrix(R: np.ndarray) -> Tuple[float, float, float]:
    """Inverse of :func:`loc_matrix`: angles (txy, tyz, txz) in degrees.

    From ``R = Rz(a) Rx(b) Ry(c)``: ``R[2,1] = sin b``, ``R[0,1] = -sin a cos b``,
    ``R[1,1] = cos a cos b``, ``R[2,0] = -cos b sin c``, ``R[2,2] = cos b cos c``.
    At gimbal lock (cos b = 0) ``c`` is set to 0.
    """
    R = np.asarray(R, float)
    sb = max(-1.0, min(1.0, R[2, 1]))
    b = math.asin(sb)
    cb = math.cos(b)
    if abs(cb) > 1e-12:
        a = math.atan2(-R[0, 1], R[1, 1])
        c = math.atan2(-R[2, 0], R[2, 2])
    else:
        c = 0.0
        a = math.atan2(R[1, 0], R[0, 0])
    return math.degrees(a), math.degrees(b), math.degrees(c)


def to_global(origin: np.ndarray, R: np.ndarray, xyz_local: np.ndarray) -> np.ndarray:
    """``X = O + R x`` for one point (3,) or many points (n, 3)."""
    x = np.asarray(xyz_local, float)
    return np.asarray(origin, float) + x @ np.asarray(R, float).T


def to_local(origin: np.ndarray, R: np.ndarray, xyz_global: np.ndarray) -> np.ndarray:
    """``x = R^T (X - O)`` for one point (3,) or many points (n, 3)."""
    X = np.asarray(xyz_global, float) - np.asarray(origin, float)
    return X @ np.asarray(R, float)
