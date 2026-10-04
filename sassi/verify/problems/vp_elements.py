"""Verification problems of the HOUSE element library: VP-03, VP-37, VP-38, VP-39.

* VP-03  discrete soil-column dispersion of the SOLID/PLANE mass schemes (R2 A.3);
* VP-37  element closed forms: Euler-Bernoulli cantilever, Timoshenko tip deflections, end
         releases, simply supported Kirchhoff plate, axial bar, patch tests (R2 I.1, spec 08 sec. 11);
* VP-38  NAFEMS free-vibration benchmarks FV12, FV16 (SHELL bending) and FV32 (SHELL membrane)
         (R2 I.2);
* VP-39  spring-mass SDOF and GENERAL-element equivalence (requirements 6.3, spec 08 sec. 11).

All checks run directly on the element library (:mod:`sassi.elements`): the models are
assembled with :func:`sassi.elements.assemble` and the modal checks use the undamped real part
of the stiffness with :func:`sassi.elements.natural_frequencies` (requirements 6.3 type E/D/P).
"""
from __future__ import annotations

import math
from typing import Callable, Dict, List, Sequence

import numpy as np
import scipy.optimize as so
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from sassi.conventions import cfactor
from sassi.elements import (ElemRecord, ElementError, assemble, build_dofmap, material_from_E_nu,
                            material_from_M, natural_frequencies)
from sassi.elements import beam, plane, shell, solid
from sassi.elements.base import blas_quiet
from sassi.verify import VPResult, problem, worse


# ======================================================================================
# helpers
# ======================================================================================
def _model(node_xyz: Dict[int, Sequence[float]], recs: List[ElemRecord], fixfun: Callable, masses=None,
           undamped: bool = False):
    ids = sorted(node_xyz)
    fix = np.array([fixfun(n, *node_xyz[n]) for n in ids], dtype=int)
    dm = build_dofmap(ids, fix, recs)
    return dm, assemble(node_xyz, recs, dm, masses=masses, undamped=undamped)


def _solve(K, F):
    with blas_quiet():
        return spla.spsolve(sp.csc_matrix(K), F)


def _closed_ratio(kh: float, scheme: str) -> float:
    """R2 A.3 discrete-to-exact frequency ratio of a 1-D linear element chain."""
    c = math.cos(kh)
    if scheme == "consistent":
        return math.sqrt(6.0 * (1.0 - c) / (2.0 + c)) / kh
    if scheme == "lumped":
        return math.sqrt(2.0 * (1.0 - c)) / kh
    return math.sqrt(12.0 * (1.0 - c) / (5.0 + c)) / kh


# ---------------------------------------------------------------------------------------
# VP-03: shear column built from SOLID or PLANE elements (faces tied to move rigidly in x)
# ---------------------------------------------------------------------------------------
_VS = 100.0
_RHO = 2.0


def _column_material():
    G = _RHO * _VS ** 2
    nu = 1.0 / 3.0
    return material_from_E_nu(2.0 * G * (1.0 + nu), nu, _RHO)


def _column(kind: str, ne: int, h: float, scheme: str):
    """Assembled column (Ne elements of height h) and the tie matrix T (level x-DOFs -> nodes)."""
    mat = _column_material()
    node_xyz, recs, level = {}, [], {}
    if kind == "SOLID":
        corners = [(0, 0), (1, 0), (1, 1), (0, 1)]
        for lv in range(ne + 1):
            for c, (x, y) in enumerate(corners):
                n = 1 + 4 * lv + c
                node_xyz[n] = (x * h, y * h, lv * h)
                level[n] = lv
        for e in range(ne):
            b = 1 + 4 * e
            recs.append(ElemRecord(1, e + 1, 1, tuple(range(b, b + 4)) + tuple(range(b + 4, b + 8)), False, mat,
                                   dict(incompatible=False, eint=0, mass=scheme)))
        fixfun = lambda n, x, y, z: [1 if level[n] == 0 else 0, 1, 1, 1, 1, 1]
    else:
        for lv in range(ne + 1):
            for c, x in enumerate((0.0, 1.0)):
                n = 1 + 2 * lv + c
                node_xyz[n] = (x * h, 0.0, lv * h)
                level[n] = lv
        for e in range(ne):
            b = 1 + 2 * e
            recs.append(ElemRecord(1, e + 1, 4, (b, b + 1, b + 3, b + 2), False, mat,
                                   dict(incompatible=False, mass=scheme)))
        fixfun = lambda n, x, y, z: [1 if level[n] == 0 else 0, 1, 1, 1, 1, 1]
    dm, am = _model(node_xyz, recs, fixfun)
    T = np.zeros((dm.neq, ne))
    for e in range(dm.neq):
        T[e, level[int(dm.eq_node[e])] - 1] = 1.0
    return am, T, mat


def _element_dispersion(kind: str, scheme: str, kh: float) -> float:
    """omega_h/omega of the infinite periodic column from one element's tied 2x2 matrices."""
    mat = _column_material()
    h = 1.0
    if kind == "SOLID":
        xyz = (solid.NODE_NAT + 1.0) / 2.0 * h
        K, M = solid.matrices(xyz, mat, mass=scheme)
        T = np.zeros((24, 2))
        for a in range(8):
            T[3 * a, 0 if a < 4 else 1] = 1.0
    else:
        xyz = np.array([[0, 0, 0], [h, 0, 0], [h, 0, h], [0, 0, h]], dtype=float)
        K, M = plane.matrices(xyz, mat, mass=scheme)
        T = np.zeros((8, 2))
        for a in range(4):
            T[2 * a, 0 if a < 2 else 1] = 1.0
    with blas_quiet():
        k = T.T @ K.real @ T
        m = T.T @ M @ T
    # Bloch wave u_j = U exp(i j kh): each node sees k11 + k22 + 2 k12 cos(kh) (same for m)
    w2 = (k[0, 0] + k[1, 1] + 2.0 * k[0, 1] * math.cos(kh)) / (m[0, 0] + m[1, 1] + 2.0 * m[0, 1] * math.cos(kh))
    return math.sqrt(w2) / (_VS * kh / h)


@problem("VP-03", "Discrete soil column dispersion (SOLID/PLANE mass schemes)", tier="P0",
         modules=["HOUSE"], source="R2 A.3")
def vp03(workdir):
    """Plane shear wave in a column of linear elements: consistent, lumped and 50/50 mass.

    Each element face is tied to move rigidly in x (the column has rollers), so the element
    reduces exactly to the 1-D shear element ``G A/h [1 -1; -1 1]`` with its mass scheme.  The
    references of R2 A.3 are the closed-form discrete dispersion ratios; the tabulated values are
    those formulas rounded to 5 decimals, so they are checked with a half-unit-of-last-digit
    tolerance (5e-6) while the element results are checked against the full-precision formulas
    with the VP tolerance 1e-6.
    """
    r = VPResult()
    tab5 = {"consistent": 1.06632, "lumped": 0.93549, "mixed": 0.99451}
    kh = 2.0 * math.pi / 5.0
    for kind in ("SOLID", "PLANE"):
        for scheme, ref in tab5.items():
            val = _element_dispersion(kind, scheme, kh)
            r.check(f"{kind} {scheme} mass: w_h/w at N = 5 per wavelength vs closed form", val,
                    _closed_ratio(kh, scheme), rtol=1e-6)
            r.check(f"{kind} {scheme} mass: w_h/w at N = 5 vs R2 A.3 table (5 decimals)", val, ref, atol=5e-6)
    # fixed-free column eigenvalues (R2 A.3 second table)
    table = {  # Ne: (consistent w1, lumped w1, mixed w1, mixed w2)
        2: (1.02586, 0.97450, 0.99919, 0.92712),
        4: (1.00644, 0.99359, 0.99995, 0.99578),
        8: (1.00161, 0.99839, 1.00000, 0.99975),
    }
    H = 8.0
    for kind in ("SOLID", "PLANE"):
        for ne, refs in table.items():
            h = H / ne
            for scheme, ref1 in zip(("consistent", "lumped", "mixed"), refs[:3]):
                am, T, mat = _column(kind, ne, h, scheme)
                with blas_quiet():
                    Kr = T.T @ (am.Ks.real @ T)
                    Mr = T.T @ (am.Ms @ T)
                f = natural_frequencies(Kr, Mr)
                w1 = 2.0 * math.pi * f[0] / (math.pi * _VS / (2.0 * H))
                kh1 = math.pi / (2.0 * ne)
                r.check(f"{kind} column Ne={ne} {scheme}: w1/w1,exact vs closed form", w1,
                        _closed_ratio(kh1, scheme), rtol=1e-6)
                r.check(f"{kind} column Ne={ne} {scheme}: w1/w1,exact vs R2 A.3 table", w1, ref1, atol=5e-6)
                if scheme == "mixed":
                    w2 = 2.0 * math.pi * f[1] / (3.0 * math.pi * _VS / (2.0 * H))
                    kh2 = 3.0 * math.pi / (2.0 * ne)
                    r.check(f"{kind} column Ne={ne} mixed: w2/w2,exact vs closed form", w2,
                            _closed_ratio(kh2, "mixed"), rtol=1e-6)
                    r.check(f"{kind} column Ne={ne} mixed: w2/w2,exact vs R2 A.3 table", w2, refs[3], atol=5e-6)
    r.notes.append("Column H = 8, Vs = 100, rho = 2; faces tied in x, rollers (UY, UZ fixed), base fixed.")
    r.notes.append("R2 A.3 tabulated ratios are 5-decimal roundings of the closed forms; the element results "
                   "match the closed forms to ~1e-15 (tolerance 1e-6) and the tables to their printed precision.")
    return r


# ======================================================================================
# VP-37: element closed forms
# ======================================================================================
def _beam_line(ne: int, L: float, mat, sec, fixfun):
    node_xyz = {i + 1: (L * i / ne, 0.0, 0.0) for i in range(ne + 1)}
    node_xyz[9999] = (0.0, 1.0, 0.0)                       # K node: axis 2 = +Y, axis 3 = +Z
    recs = [ElemRecord(1, i + 1, 2, (i + 1, i + 2, 9999), False, mat, dict(section=sec)) for i in range(ne)]
    return _model(node_xyz, recs, fixfun)


def _plate_mesh(nx: int, ny: int, a: float, b: float, mat, t: float, origin=(0.0, 0.0), xy_of=None):
    node_xyz, recs = {}, []
    nid = lambda i, j: 1 + i + j * (nx + 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            x, y = origin[0] + a * i / nx, origin[1] + b * j / ny
            if xy_of is not None:
                x, y = xy_of(i, j)
            node_xyz[nid(i, j)] = (x, y, 0.0)
    e = 0
    for j in range(ny):
        for i in range(nx):
            e += 1
            recs.append(ElemRecord(1, e, 3, (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)), False, mat,
                                   dict(thick=t)))
    return node_xyz, recs


# MacNeal & Harder (1985) patch geometries
_MH2D_OUT = [(0.0, 0.0), (0.24, 0.0), (0.24, 0.12), (0.0, 0.12)]
_MH2D_IN = [(0.04, 0.02), (0.18, 0.03), (0.16, 0.08), (0.08, 0.08)]
_MH2D_ELEMS = [(0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7)]   # indices: out 0-3, in 4-7
_MH3D_IN = [(0.249, 0.342, 0.192), (0.826, 0.288, 0.288), (0.850, 0.649, 0.263), (0.273, 0.750, 0.230),
            (0.320, 0.186, 0.643), (0.677, 0.305, 0.683), (0.788, 0.693, 0.644), (0.165, 0.745, 0.702)]
_HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def _solid_patch(incompatible: bool, eint: int) -> float:
    """MacNeal-Harder 7-element distorted hexahedron patch; returns max relative stress error."""
    mat = material_from_E_nu(1.0e6, 0.25, 1.0)
    outer = (solid.NODE_NAT + 1.0) / 2.0
    xyz = np.vstack([outer, np.array(_MH3D_IN)])
    node_xyz = {i + 1: xyz[i] for i in range(16)}
    conn = [tuple(range(8, 16))]
    for f in _HEX_FACES:
        el = [f[k] for k in range(4)] + [8 + f[k] for k in range(4)]
        X = xyz[el]
        _, dN0 = solid.hex8_shape(np.zeros((1, 3)))
        if np.linalg.det(dN0[0].T @ X) < 0:                  # orient: outer face first or second
            el = el[4:] + el[:4]
        conn.append(tuple(el))
    recs = [ElemRecord(1, i + 1, 1, tuple(n + 1 for n in c), False, mat, dict(incompatible=incompatible, eint=eint))
            for i, c in enumerate(conn)]
    eps = 1e-3 * np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])     # u = 1e-3 (2x+y+z)/2 etc.
    Eps = 1e-3 * np.array([[1.0, 0.5, 0.5], [0.5, 1.0, 0.5], [0.5, 0.5, 1.0]])
    ufield = lambda p: Eps @ np.asarray(p)
    return _patch_solve(node_xyz, recs, lambda n: n <= 8, ufield, "SOLID",
                        (solid.isotropic_D(mat.lam, mat.G) @ eps).real)


def _patch_solve(node_xyz, recs, is_boundary, ufield, rec_name, s_exact, rot_field=None,
                 free_dofs=None) -> float:
    """Prescribe the exact field on boundary nodes, solve the interior, return the max error of the
    recovered constant field relative to its largest component."""
    ids = sorted(node_xyz)
    dm_all = build_dofmap(ids, np.zeros((len(ids), 6), dtype=int), recs)
    am = assemble(node_xyz, recs, dm_all)
    u_ex = np.zeros(dm_all.neq)
    for e in range(dm_all.neq):
        n, d = int(dm_all.eq_node[e]), int(dm_all.eq_dof[e])
        if d <= 3:
            u_ex[e] = ufield(node_xyz[n])[d - 1]
        elif rot_field is not None:
            u_ex[e] = rot_field(node_xyz[n])[d - 4]
    free = np.array([not is_boundary(int(dm_all.eq_node[e])) and
                     (free_dofs is None or int(dm_all.eq_dof[e]) in free_dofs) for e in range(dm_all.neq)])
    K = am.Ks.real.tocsr()
    u = u_ex.copy()
    if free.any():
        Kff = K[free][:, free]
        Kfb = K[free][:, ~free]
        u[free] = _solve(Kff, -(Kfb @ u_ex[~free]))
    rec = am.recovery[rec_name]
    ue = np.where(rec["eq"] >= 0, u[np.maximum(rec["eq"], 0)], 0.0)
    s = np.einsum("ecd,ed->ec", rec["S"].real, ue)
    ncomp = len(s_exact)
    err = np.abs(s[:, :ncomp] - np.asarray(s_exact)[None, :]).max()
    return float(err / np.abs(s_exact).max())


def _plane_patch(incompatible: bool) -> float:
    mat = material_from_E_nu(1.0e6, 0.25, 1.0)
    pts = _MH2D_OUT + _MH2D_IN
    node_xyz = {i + 1: (p[0], 0.0, p[1]) for i, p in enumerate(pts)}
    recs = [ElemRecord(1, i + 1, 4, tuple(n + 1 for n in c), False, mat, dict(incompatible=incompatible))
            for i, c in enumerate(_MH2D_ELEMS)]
    ufield = lambda p: 1e-3 * np.array([p[0] + 0.5 * p[2], 0.0, p[2] + 0.5 * p[0]])
    eps = 1e-3 * np.array([1.0, 1.0, 1.0])
    return _patch_solve(node_xyz, recs, lambda n: n <= 4, ufield, "PLANE",
                        (plane.plane_strain_D(mat.lam, mat.G) @ eps).real)


def _rotation(seed: int = 7) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return Q if np.linalg.det(Q) > 0 else -Q


def _shell_membrane_patch(triangles: bool) -> float:
    """MacNeal-Harder membrane patch rotated into an oblique plane (tests the transformations)."""
    mat = material_from_E_nu(1.0e6, 0.25, 1.0)
    t = 0.01
    R = _rotation()
    pts = [np.array([p[0], p[1], 0.0]) for p in _MH2D_OUT + _MH2D_IN]
    node_xyz = {i + 1: R @ p for i, p in enumerate(pts)}
    recs = []
    for i, c in enumerate(_MH2D_ELEMS):
        n = tuple(k + 1 for k in c)
        if triangles:
            recs.append(ElemRecord(1, 2 * i + 1, 3, (n[0], n[1], n[2]), False, mat, dict(thick=t)))
            recs.append(ElemRecord(1, 2 * i + 2, 3, (n[0], n[2], n[3]), False, mat, dict(thick=t)))
        else:
            recs.append(ElemRecord(1, i + 1, 3, n, False, mat, dict(thick=t)))
    # in-plane field in the patch coordinates (x, y): u = 1e-3 (x + y/2), v = 1e-3 (y + x/2)
    def ufield(p):
        q = R.T @ np.asarray(p)
        loc = 1e-3 * np.array([q[0] + 0.5 * q[1], q[1] + 0.5 * q[0], 0.0])
        return R @ loc
    E, nu = mat.E0, mat.nu0
    s_exact = E / (1 - nu ** 2) * np.array([1e-3 * (1 + nu), 1e-3 * (1 + nu), 0.5 * (1 - nu) * 1e-3])
    # each element's x' axis follows its node numbering: the stresses are rotated to the patch axes
    return _patch_solve_shell_membrane(node_xyz, recs, R, ufield, s_exact)


def _patch_solve_shell_membrane(node_xyz, recs, R, ufield, s_exact) -> float:
    ids = sorted(node_xyz)
    dm = build_dofmap(ids, np.zeros((len(ids), 6), dtype=int), recs)
    am = assemble(node_xyz, recs, dm)
    u_ex = np.zeros(dm.neq)
    for e in range(dm.neq):
        n, d = int(dm.eq_node[e]), int(dm.eq_dof[e])
        if d <= 3:
            u_ex[e] = ufield(node_xyz[n])[d - 1]
    # interior nodes: translations free, rotations fixed (exact rotations are zero)
    free = np.array([int(dm.eq_node[e]) > 4 and int(dm.eq_dof[e]) <= 3 for e in range(dm.neq)])
    K = am.Ks.real.tocsr()
    u = u_ex.copy()
    u[free] = _solve(K[free][:, free], -(K[free][:, ~free] @ u_ex[~free]))
    rec = am.recovery["SHELL"]
    ue = np.where(rec["eq"] >= 0, u[np.maximum(rec["eq"], 0)], 0.0)
    s = np.einsum("ecd,ed->ec", rec["S"].real, ue)
    # rotate each element's local membrane stress tensor to the patch axes
    errs = []
    for k, r_ in enumerate(recs):
        Lam, _, _ = shell.local_frame(np.array([node_xyz[n] for n in r_.nodes]))
        A2 = (Lam @ R)[:2, :2]                                  # local x'y' in patch xy coordinates
        sig_loc = np.array([[s[k, 0], s[k, 2]], [s[k, 2], s[k, 1]]])
        sig = A2.T @ sig_loc @ A2
        sx = np.array([sig[0, 0], sig[1, 1], sig[0, 1]])
        errs.append(np.abs(sx - s_exact).max())
        errs.append(np.abs(s[k, 3:]).max())                     # no bending
    return float(worse(*errs) / np.abs(s_exact).max())


def _shell_bending_patch(triangles: bool) -> float:
    """Constant-curvature plate patch (w quadratic) on the distorted MacNeal-Harder mesh."""
    mat = material_from_E_nu(1.0e6, 0.25, 1.0)
    t = 0.01
    pts = _MH2D_OUT + _MH2D_IN
    node_xyz = {i + 1: (p[0], p[1], 0.0) for i, p in enumerate(pts)}
    recs = []
    for i, c in enumerate(_MH2D_ELEMS):
        n = tuple(k + 1 for k in c)
        if triangles:
            recs.append(ElemRecord(1, 2 * i + 1, 3, (n[0], n[1], n[2]), False, mat, dict(thick=t)))
            recs.append(ElemRecord(1, 2 * i + 2, 3, (n[0], n[2], n[3]), False, mat, dict(thick=t)))
        else:
            recs.append(ElemRecord(1, i + 1, 3, n, False, mat, dict(thick=t)))
    k1, k2, k12 = 1e-3, 2e-3, 0.5e-3
    w = lambda p: 0.5 * k1 * p[0] ** 2 + 0.5 * k2 * p[1] ** 2 + k12 * p[0] * p[1]
    wx = lambda p: k1 * p[0] + k12 * p[1]
    wy = lambda p: k2 * p[1] + k12 * p[0]
    ufield = lambda p: np.array([0.0, 0.0, w(p)])
    rot = lambda p: np.array([wy(p), -wx(p), 0.0])             # theta_x = w,y ; theta_y = -w,x
    D = mat.E0 * t ** 3 / (12 * (1 - mat.nu0 ** 2))
    nu = mat.nu0
    m_exact = -D * np.array([k1 + nu * k2, k2 + nu * k1, (1 - nu) * k12])
    ids = sorted(node_xyz)
    dm = build_dofmap(ids, np.zeros((len(ids), 6), dtype=int), recs)
    am = assemble(node_xyz, recs, dm)
    u_ex = np.zeros(dm.neq)
    for e in range(dm.neq):
        n, d = int(dm.eq_node[e]), int(dm.eq_dof[e])
        u_ex[e] = ufield(node_xyz[n])[d - 1] if d <= 3 else rot(node_xyz[n])[d - 4]
    free = np.array([int(dm.eq_node[e]) > 4 and int(dm.eq_dof[e]) in (3, 4, 5) for e in range(dm.neq)])
    K = am.Ks.real.tocsr()
    u = u_ex.copy()
    u[free] = _solve(K[free][:, free], -(K[free][:, ~free] @ u_ex[~free]))
    rec = am.recovery["SHELL"]
    ue = np.where(rec["eq"] >= 0, u[np.maximum(rec["eq"], 0)], 0.0)
    s = np.einsum("ecd,ed->ec", rec["S"].real, ue)
    errs = []
    for k, r_ in enumerate(recs):
        Lam, _, _ = shell.local_frame(np.array([node_xyz[n] for n in r_.nodes]))
        A2 = Lam[:2, :2]
        m_loc = np.array([[s[k, 3], s[k, 5]], [s[k, 5], s[k, 4]]])
        mg = A2.T @ m_loc @ A2
        errs.append(np.abs(np.array([mg[0, 0], mg[1, 1], mg[0, 1]]) - m_exact).max())
    return float(worse(*errs) / np.abs(m_exact).max())


@problem("VP-37", "HOUSE element closed forms (beams, plates, bars, patch tests)", tier="P0",
         modules=["HOUSE"], source="R2 I.1, spec 08 sec. 11")
def vp37(workdir):
    r = VPResult()
    E, nu, rho = 2.0e5, 0.3, 7.8
    mat = material_from_E_nu(E, nu, rho)
    G = mat.G0
    L = 10.0

    # ---- 1. Euler-Bernoulli cantilever (planar bending in X-Y, As = 0), 20 elements ------------
    sec = dict(A=0.02, As2=0.0, As3=0.0, J=3.0e-5, I2=4.0e-5, I3=1.0e-4)
    ne = 20
    dm, am = _beam_line(ne, L, mat, sec, lambda n, x, y, z: [1] * 6 if n == 1 or n == 9999 else [1, 0, 1, 1, 1, 0])
    f = natural_frequencies(am.Ks, am.Ms, 5)
    mbar = rho * sec["A"]
    for k, ref in enumerate([1.875104, 4.694091, 7.854757, 10.995541, 14.137168]):
        bl = L * ((2 * math.pi * f[k]) ** 2 * mbar / (E * sec["I3"])) ** 0.25
        r.check(f"cantilever beta_{k + 1} L (20 BEAMS elements)", bl, ref, rtol=0.01)

    # ---- 2. Timoshenko tip deflections (single element) -----------------------------------------
    xyz = np.array([[0, 0, 0], [L, 0, 0], [0, 1, 0]], dtype=float)
    P = 1.0
    for As2, As3, label in ((0.0, 0.0, "no shear area"), (0.015, 0.012, "with shear areas")):
        s2 = dict(sec, As2=As2, As3=As3)
        K, _ = beam.matrices(xyz, mat, section=s2)
        Kjj = K.real[6:, 6:]
        u2 = np.linalg.solve(Kjj, P * np.eye(6)[1])[1]
        u3 = np.linalg.solve(Kjj, P * np.eye(6)[2])[2]
        ref2 = P * L ** 3 / (3 * E * sec["I3"]) + (P * L / (G * As2) if As2 > 0 else 0.0)
        ref3 = P * L ** 3 / (3 * E * sec["I2"]) + (P * L / (G * As3) if As3 > 0 else 0.0)
        r.check(f"tip deflection along axis 2, {label}", u2, ref2, rtol=1e-10)
        r.check(f"tip deflection along axis 3, {label}", u3, ref3, rtol=1e-10)
    # rotated beam with an arbitrary K node: the response along e2 is unchanged (spec 08 11.9)
    s2 = dict(sec, As2=0.015, As3=0.012)
    Rm = _rotation(11)
    K_node = Rm @ np.array([3.0, 2.5, 0.0])                    # any point in the rotated 1-2 plane, +e2 side
    xyzR = np.vstack([np.zeros(3), Rm @ np.array([L, 0, 0]), K_node])
    K, _ = beam.matrices(xyzR, mat, section=s2)
    e2 = Rm[:, 1]
    u = np.linalg.solve(K.real[6:, 6:], np.concatenate([P * e2, np.zeros(3)]))
    r.check("rotated beam (arbitrary K): tip deflection along e2", u[:3] @ e2,
            P * L ** 3 / (3 * E * sec["I3"]) + P * L / (G * 0.015), rtol=1e-10)

    # ---- 3. end release: lateral stiffness 12EI/L^3 -> 3EI/L^3 ----------------------------------
    for kj, ref, label in ((None, 12 * E * sec["I3"] / L ** 3, "no release"),
                           ("000001", 3 * E * sec["I3"] / L ** 3, "KJ releases M3")):
        node_xyz = {1: (0, 0, 0), 2: (L, 0, 0), 3: (0, 1, 0)}
        recs = [ElemRecord(1, 1, 2, (1, 2, 3), False, mat, dict(section=sec, kj=kj))]
        dm1, am1 = _model(node_xyz, recs, lambda n, x, y, z: [1] * 6 if n != 2 else [1, 0, 1, 1, 1, 1])
        r.check(f"lateral stiffness at J along +Y, {label}", am1.Ks.real.toarray()[0, 0], ref, rtol=1e-10)
    try:
        beam.matrices(xyz, mat, section=sec, ki="000001", kj="000001")
        raised = False
    except ElementError as exc:
        raised = "Error 10" in str(exc)
    r.require("M3 released at both I and J -> Error 10", raised)

    # ---- 4. axial bar, fixed-free, 50 elements ------------------------------------------------
    dmb, amb = _beam_line(50, L, mat, sec, lambda n, x, y, z: [1] * 6 if n == 1 or n == 9999 else [0, 1, 1, 1, 1, 1])
    fb = natural_frequencies(amb.Ks, amb.Ms, 3)
    for n in range(1, 4):
        r.check(f"bar f_{n} = (2n-1) sqrt(E/rho)/(4L) (50 BEAMS elements)", fb[n - 1],
                (2 * n - 1) * math.sqrt(E / rho) / (4 * L), rtol=0.01)

    # ---- 5. simply supported Kirchhoff plate (SHELL, 32 x 16 mesh on 2 x 1) ----------------------
    pm = material_from_E_nu(200e9, 0.3, 8000.0)
    t, a, b = 0.05, 2.0, 1.0
    node_xyz, recs = _plate_mesh(32, 16, a, b, pm, t)
    on_edge = lambda x, y: abs(x) < 1e-9 or abs(x - a) < 1e-9 or abs(y) < 1e-9 or abs(y - b) < 1e-9
    dmp, amp = _model(node_xyz, recs, lambda n, x, y, z: [1, 1, 1 if on_edge(x, y) else 0, 0, 0, 1])
    fp = natural_frequencies(amp.Ks, amp.Ms, 4)
    Dp = 200e9 * t ** 3 / (12 * (1 - 0.09))
    mn = sorted(((m / a) ** 2 + (n / b) ** 2, m, n) for m in range(1, 6) for n in range(1, 4))
    for k in range(4):
        s_, m_, n_ = mn[k]
        r.check(f"SS plate f_{m_}{n_} (SHELL 32x16)", fp[k], math.pi / 2 * s_ * math.sqrt(Dp / (8000.0 * t)),
                rtol=0.01)

    # ---- 6. patch tests: exact constant stress --------------------------------------------------
    for inc in (False, True):
        for eint in (0, 1, 2):
            r.check(f"SOLID patch (MacNeal-Harder, incompatible={inc}, EINT={eint}): max stress error",
                    _solid_patch(inc, eint), 0.0, atol=1e-10)
        r.check(f"PLANE patch (MacNeal-Harder, incompatible={inc}): max stress error", _plane_patch(inc), 0.0,
                atol=1e-10)
    r.check("SHELL membrane patch, oblique plane, Q4 + incompatible modes", _shell_membrane_patch(False), 0.0,
            atol=1e-10)
    r.check("SHELL membrane patch, oblique plane, CST triangles", _shell_membrane_patch(True), 0.0, atol=1e-10)
    r.check("SHELL bending patch (constant curvature), DKQ", _shell_bending_patch(False), 0.0, atol=1e-10)
    r.check("SHELL bending patch (constant curvature), DKT", _shell_bending_patch(True), 0.0, atol=1e-10)

    # ---- 7. complex stiffness: K* = c(beta) K0 when beta_p = beta_s (requirements 4.0.2) -----------
    bdamp = 0.06
    md = material_from_M(1, 1.0e4, 0.3, 2.0, bdamp, bdamp, 1.0)
    rng = np.random.default_rng(21)
    hexa = (solid.NODE_NAT + 1.0) / 2.0 + rng.uniform(-0.05, 0.05, (8, 3))
    quad = np.array([[0.0, 0, 0], [1.2, 0, 0.1], [1.1, 0, 1.0], [0.0, 0, 0.9]])
    cases = (("SOLID + incompatible modes", solid.matrices, hexa, dict(incompatible=True)),
             ("PLANE + incompatible modes", plane.matrices, quad, dict(incompatible=True)),
             ("BEAMS with shear areas and KI release", beam.matrices, np.array([[0, 0, 0], [3.0, 1, 0], [0, 0, 1.0]]),
              dict(section=dict(sec, As2=0.015, As3=0.012), ki="000011")),
             ("SHELL quad", shell.matrices, quad[:, [0, 2, 1]], dict(thick=0.1)))
    for label, fn, xyz_c, kw in cases:
        Kc, _ = fn(xyz_c, md, **kw)
        K0, _ = fn(xyz_c, md.undamped(), **kw)
        r.check(f"{label}: max|K* - c(beta) K0| / max|K0|", np.abs(Kc - cfactor(bdamp) * K0.real).max()
                / np.abs(K0).max(), 0.0, atol=1e-12)

    # ---- informative: SOLID cantilever with / without incompatible modes -------------------------
    r.notes.append(_solid_cantilever_note())
    r.notes.append("Beams: E = 2e5, nu = 0.3, rho = 7.8, L = 10, A = 0.02, I3 = 1e-4, I2 = 4e-5; "
                   "plate: E = 200 GPa, nu = 0.3, rho = 8000, t = 0.05, a x b = 2 x 1.")
    return r


def _solid_cantilever_note() -> str:
    """Tip deflection of a 10 x 1 x 1 SOLID cantilever (10 x 1 x 1 elements) / Timoshenko beam."""
    mat = material_from_E_nu(1.0e4, 0.0, 1.0)
    L, h = 10.0, 1.0
    nx = 10
    node_xyz, idx = {}, {}
    nid = 0
    for i in range(nx + 1):
        for j in range(2):
            for k in range(2):
                nid += 1
                node_xyz[nid] = (L * i / nx, h * j, h * k)
                idx[(i, j, k)] = nid
    res = {}
    for inc in (False, True):
        recs = []
        for i in range(nx):
            c = [idx[(i, 0, 0)], idx[(i + 1, 0, 0)], idx[(i + 1, 1, 0)], idx[(i, 1, 0)],
                 idx[(i, 0, 1)], idx[(i + 1, 0, 1)], idx[(i + 1, 1, 1)], idx[(i, 1, 1)]]
            recs.append(ElemRecord(1, i + 1, 1, tuple(c), False, mat, dict(incompatible=inc)))
        dm, am = _model(node_xyz, recs, lambda n, x, y, z: [1] * 6 if abs(x) < 1e-9 else [0, 0, 0, 1, 1, 1])
        F = np.zeros(dm.neq)
        tip = [idx[(nx, j, k)] for j in range(2) for k in range(2)]
        for n in tip:
            F[dm.eq(n, 3)] = 0.25
        u = _solve(am.Ks.real, F)
        res[inc] = np.mean([u[dm.eq(n, 3)] for n in tip])
    I = h ** 4 / 12.0
    ref = L ** 3 / (3 * mat.E0 * I) + L / (mat.G0 * h * h / 1.2)
    return ("SOLID cantilever 10x1x1 (10 elements, tip shear 1): tip deflection / Timoshenko beam = "
            f"{res[False] / ref:.4f} without and {res[True] / ref:.4f} with incompatible modes "
            "(spec 08 11.14: incompatible modes approach beam theory).")


# ======================================================================================
# VP-38: NAFEMS free vibration
# ======================================================================================
_NAFEMS_MAT = dict(E=200e9, nu=0.3, rho=8000.0)


def _nafems_material():
    return material_from_E_nu(_NAFEMS_MAT["E"], _NAFEMS_MAT["nu"], _NAFEMS_MAT["rho"])


def fv12(n: int = 8) -> np.ndarray:
    """FV12 free thin square plate 10 x 10 x 0.05 m, in-plane restrained (u = v = rot_z = 0)."""
    node_xyz, recs = _plate_mesh(n, n, 10.0, 10.0, _nafems_material(), 0.05)
    _, am = _model(node_xyz, recs, lambda k, x, y, z: [1, 1, 0, 0, 0, 1])
    return natural_frequencies(am.Ks, am.Ms, 9)[3:]


def fv16(n: int = 8) -> np.ndarray:
    """FV16 cantilevered thin square plate 10 x 10 x 0.05 m clamped along x = 0."""
    node_xyz, recs = _plate_mesh(n, n, 10.0, 10.0, _nafems_material(), 0.05)
    _, am = _model(node_xyz, recs, lambda k, x, y, z: [1] * 6 if abs(x) < 1e-9 else [1, 1, 0, 0, 0, 1])
    return natural_frequencies(am.Ks, am.Ms, 6)


def fv32(nx: int = 16, ny: int = 8) -> np.ndarray:
    """FV32 cantilevered tapered membrane: corners (0, -2.5), (10, -0.5), (10, 0.5), (0, 2.5) m,
    root x = 0 clamped, out-of-plane DOFs restrained (geometry of the Abaqus Benchmarks Guide
    input file nfv32i4f.inp, 16 x 8 CPS4I mesh)."""
    def xy(i, j):
        x = 10.0 * i / nx
        half = 2.5 - 0.2 * x
        return x, -half + 2.0 * half * j / ny
    node_xyz, recs = _plate_mesh(nx, ny, 10.0, 5.0, _nafems_material(), 0.05, xy_of=xy)
    _, am = _model(node_xyz, recs, lambda k, x, y, z: [1] * 6 if abs(x) < 1e-9 else [0, 0, 1, 1, 1, 1])
    return natural_frequencies(am.Ks, am.Ms, 6)


FV12_REF = [1.622, 2.360, 2.922, 4.233, 4.233, 7.416]
FV16_REF = [0.421, 1.029, 2.582, 3.306, 3.753, 6.555]
FV32_REF = [44.623, 130.03, 162.70, 246.05, 379.90, 391.44]


@problem("VP-38", "NAFEMS free vibration FV12, FV16, FV32 (SHELL)", tier="P0", modules=["HOUSE"],
         source="R2 I.2")
def vp38(workdir):
    """NAFEMS FV12 / FV16 (Kirchhoff bending, DKQ + lumped mass) and FV32 (membrane, Q4 with
    incompatible modes + lumped mass); 2 % on the first 4 elastic modes.

    Lead decision (wave-1 review): the manual's SHELL uses lumped translational mass without rotary
    inertia (D-ELM-07), which converges from below at O(h^2) and is 2.4-6.7 % low on the coarse NAFEMS
    8x8 mesh.  The element is therefore verified on a converged 32x32 mesh (FV12, FV16) together with
    an observed convergence order of about 2; the 8x8 and 16x16 results are reported in the notes.
    FV32 is checked on the NAFEMS 16x8 mesh."""
    r = VPResult()
    meshes = {n: (fv12(n), fv16(n)) for n in (8, 16, 32)}
    f12, f16 = meshes[32]
    f32 = fv32(16, 8)
    for name, f, ref, mesh in (("FV12", f12, FV12_REF, "32x32"), ("FV16", f16, FV16_REF, "32x32"),
                               ("FV32", f32, FV32_REF, "16x8")):
        for k in range(4):
            r.check(f"{name} mode {k + 1} (SHELL {mesh}, Hz)", f[k], ref[k], rtol=0.02)
    # observed order of convergence of mode 1 (Richardson, meshes 8/16/32)
    for j, name in ((0, "FV12"), (1, "FV16")):
        a, b, c = (meshes[n][j][0] for n in (8, 16, 32))
        order = np.log2(abs(b - a) / abs(c - b))
        r.check(f"{name} mode 1 observed convergence order (8/16/32 meshes)", order, 2.0, rtol=0.15)
    for n, (g12, g16) in meshes.items():
        r.notes.append(f"{n}x{n}: FV12 errors " + ", ".join(f"{100 * (v / q - 1):+.2f}%" for v, q in zip(g12[:6], FV12_REF))
                       + "; FV16 errors " + ", ".join(f"{100 * (v / q - 1):+.2f}%" for v, q in zip(g16[:6], FV16_REF)))
    r.notes.append("FV32 16x8: modes 1-6 = " + ", ".join(f"{v:.4f}" for v in f32[:6]) + " Hz")
    r.notes.append("DKQ (Batoz & Ben Tahar) with the lumped translational mass of D-ELM-07 converges from below at "
                   "O(h^2).  The FV16 fundamental converges to lambda1 = w a^2 sqrt(rho h / D) = 3.471 (0.4179 Hz), the "
                   "value of the accurate superposition solutions for the cantilever square plate; the NAFEMS target "
                   "0.421 Hz (lambda1 = 3.49, a Ritz upper bound) is 0.6 % higher.  FV32 uses the geometry of the "
                   "Abaqus input file nfv32i4f.inp (symmetric taper, 16 x 8 CPS4I mesh).")
    return r


# ======================================================================================
# VP-39: spring-mass SDOF and GENERAL equivalence
# ======================================================================================
def _spring_pattern(kvec) -> np.ndarray:
    K = np.zeros((12, 12))
    d = np.arange(6)
    K[d, d] = kvec
    K[d + 6, d + 6] = kvec
    K[d, d + 6] = -np.asarray(kvec)
    K[d + 6, d] = -np.asarray(kvec)
    return K


@problem("VP-39", "Spring-mass SDOF and GENERAL element equivalence", tier="P0", modules=["HOUSE"],
         source="requirements 6.3, spec 08 sec. 11 (C-28)")
def vp39(workdir):
    r = VPResult()
    k, beta, m = 1000.0, 0.05, 1.0
    # ---- SPRING SC,1,1000,0,0,0,0,0,0.05 between fixed node 1 and node 2 with MT = 1 (mass units)
    node_xyz = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0)}
    sc = dict(k=(k, 0, 0, 0, 0, 0), damp=beta)
    fixfun = lambda n, x, y, z: [1] * 6 if n == 1 else [0, 1, 1, 1, 1, 1]
    rec_s = [ElemRecord(1, 1, 7, (1, 2), False, None, sc)]
    dm, am = _model(node_xyz, rec_s, fixfun, masses={2: (m, m, m, 0, 0, 0)})
    _, am0 = _model(node_xyz, rec_s, fixfun, masses={2: (m, m, m, 0, 0, 0)}, undamped=True)
    f = natural_frequencies(am0.Ks, am0.Ms)                   # undamped K0 (Re K* = K0 (1 - 2 b^2))
    r.check("SPRING-mass undamped frequency vs sqrt(k/m)/(2 pi)", f[0], math.sqrt(k / m) / (2 * math.pi), rtol=1e-8)
    r.check("SPRING-mass undamped frequency vs 5.033 Hz (3 decimals)", f[0], 5.033, atol=5e-4)
    ks = am.Ks.toarray()[0, 0]
    r.check("damped K* / (k c(beta)) (real part)", (ks / (k * cfactor(beta))).real, 1.0, rtol=1e-12)
    r.check("damped K* / (k c(beta)) (imaginary part)", (ks / (k * cfactor(beta))).imag, 0.0, atol=1e-12)
    # fixed-base total-acceleration transfer function k*/(k* - m w^2): peak 1/(2 b sqrt(1-b^2))
    w0 = math.sqrt(k / m)
    H = lambda w: abs(ks / (ks - m * w * w))
    opt = so.minimize_scalar(lambda x: -H(x * w0), bounds=(0.9, 1.1), method="bounded",
                             options={"xatol": 1e-12})
    peak = 1.0 / (2 * beta * math.sqrt(1 - beta * beta))
    r.check("ATF peak |H|max vs 1/(2 beta sqrt(1-beta^2))", H(opt.x * w0), peak, rtol=1e-8)
    r.check("ATF peak |H|max vs 10.0125 (4 decimals)", H(opt.x * w0), 10.0125, atol=5e-5)
    r.check("ATF peak location w/w0 vs sqrt(1 - 2 beta^2)", opt.x, math.sqrt(1 - 2 * beta * beta), rtol=1e-6)

    # ---- GENERAL, 2 nodes, global input: MXR = k(1-2b^2), MXI = 2kb sqrt(1-b^2) in the spring pattern
    kr = k * (1 - 2 * beta * beta)
    ki = 2 * k * beta * math.sqrt(1 - beta * beta)
    rows_R = {1: [kr, 0, 0, 0, 0, 0, -kr], 7: [kr]}               # MXR,1,1,kr,0,0,0,0,0,-kr ; MXR,1,7,kr
    rows_I = {1: [ki, 0, 0, 0, 0, 0, -ki], 7: [ki]}
    rec_g = [ElemRecord(1, 1, 9, (1, 2), False, None, dict(KR=rows_R, KI=rows_I))]
    dmg, amg = _model(node_xyz, rec_g, fixfun, masses={2: (m, m, m, 0, 0, 0)})
    dK = abs(amg.Ks - am.Ks).max()
    r.check("GENERAL (2-node global) K* - SPRING K* (max abs / k)", dK / k, 0.0, atol=1e-8)
    # same damped response: |k*/(k* - m w^2)| at the SPRING peak
    Hg = abs(amg.Ks.toarray()[0, 0] / (amg.Ks.toarray()[0, 0] - m * (opt.x * w0) ** 2))
    r.check("GENERAL (2-node global) ATF peak vs SPRING", Hg, H(opt.x * w0), rtol=1e-8)
    # undamped GM (MXR = k, MXI = 0): same frequency as the undamped SPRING
    rows_R0 = {1: [k, 0, 0, 0, 0, 0, -k], 7: [k]}
    rec_g0 = [ElemRecord(1, 1, 9, (1, 2), False, None, dict(KR=rows_R0))]
    _, amg0 = _model(node_xyz, rec_g0, fixfun, masses={2: (m, m, m, 0, 0, 0)})
    fg = natural_frequencies(amg0.Ks, amg0.Ms)
    r.check("GENERAL (2-node global, MXR = k) frequency vs SPRING", fg[0], f[0], rtol=1e-8)
    # mass entered in the GM matrix instead of MT, in weight units (MOPT <matrix> = 1)
    g = 32.2
    MM = np.zeros((12, 12))
    MM[6, 6] = m * g
    rec_gm = [ElemRecord(1, 1, 9, (1, 2), False, None, dict(KR=rows_R0, MM=MM, munits=1, gravity=g))]
    _, amgm = _model(node_xyz, rec_gm, fixfun)
    r.check("GENERAL mass in weight units / g: frequency vs SPRING + MT", natural_frequencies(amgm.Ks, amgm.Ms)[0],
            f[0], rtol=1e-8)

    # ---- GENERAL, 3 nodes, local input, axes aligned with permuted global axes -----------------
    # I->J along +Y, K on +Z: e1 = Y, e2 = Z, e3 = X.  Local constants (k1..k6) map to global
    # kx = k3, ky = k1, kz = k2, kxx = k6, kyy = k4, kzz = k5.
    kloc = np.array([1000.0, 2000.0, 3000.0, 40.0, 50.0, 60.0])
    c = cfactor(beta)
    node3 = {1: (0.0, 0.0, 0.0), 2: (0.0, 1.0, 0.0), 3: (0.0, 0.0, 1.0)}
    rec_l = [ElemRecord(1, 1, 9, (1, 2, 3), False, None,
                        dict(KR=_spring_pattern(kloc * c.real), KI=_spring_pattern(kloc * c.imag)))]
    kglob = (kloc[2], kloc[0], kloc[1], kloc[5], kloc[3], kloc[4])
    rec_sp = [ElemRecord(1, 1, 7, (1, 2), False, None, dict(k=kglob, damp=beta))]
    free_all = lambda n, x, y, z: [1] * 6 if n == 1 else [0] * 6
    dml, aml = _model(node3, rec_l, free_all)
    dms, ams = _model(node3, rec_sp, free_all)
    r.check("GENERAL (3-node local, permuted axes) K* - equivalent SPRING K* (max abs / max k)",
            abs(aml.Ks - ams.Ks).max() / kloc.max(), 0.0, atol=1e-8)

    # ---- GENERAL, 3 nodes, oblique axial spring along I->J ------------------------------------
    R = _rotation(5)
    n_ax = R[:, 0]
    nodeo = {1: (0.0, 0.0, 0.0), 2: tuple(2.0 * n_ax), 3: tuple(R[:, 1])}
    kax = np.array([k, 0, 0, 0, 0, 0])
    rec_o = [ElemRecord(1, 1, 9, (1, 2, 3), False, None,
                        dict(KR=_spring_pattern(kax * c.real), KI=_spring_pattern(kax * c.imag)))]
    fixo = lambda nn, x, y, z: [1] * 6 if nn == 1 else [0, 0, 0, 1, 1, 1]
    dmo, amo = _model(nodeo, rec_o, fixo, masses={2: (m, m, m, 0, 0, 0)})
    Kt = amo.Ks.toarray()
    r.check("GENERAL oblique axial spring: K* translational block - k c n n^T (max abs / k)",
            np.abs(Kt - k * c * np.outer(n_ax, n_ax)).max() / k, 0.0, atol=1e-8)
    rec_o0 = [ElemRecord(1, 1, 9, (1, 2, 3), False, None, dict(KR=_spring_pattern(kax)))]   # undamped MXR = k
    _, amo0 = _model(nodeo, rec_o0, fixo, masses={2: (m, m, m, 0, 0, 0)})
    fo, modes = natural_frequencies(amo0.Ks, amo0.Ms, return_modes=True)
    r.check("GENERAL oblique axial spring: frequency along I->J", fo[-1], math.sqrt(k / m) / (2 * math.pi), rtol=1e-8)
    r.check("GENERAL oblique axial spring: mode direction |phi . n|", abs(modes[:, -1] @ n_ax), 1.0, rtol=1e-8)
    r.notes.append("Reference values 5.033 Hz and 10.0125 are 3- and 4-decimal roundings; the exact closed forms "
                   "are checked with the VP tolerance 1e-8.")
    return r
