"""Verification problems of the TSHELL thick-shell element and the HOUSE node optimizer (tier P2/P1):
VP-38T, VP-54T, VP-TS1, VP-O1.

* VP-38T  TSHELL part of VP-38 (R2 I.2): NAFEMS Test 21 simply supported thick plate 10 x 10 x 1 m
          (E = 200 GPa, nu = 0.3, rho = 8000 kg/m^3; hard simple supports w = 0 and the tangential
          rotation fixed on every edge, in-plane restrained): 45.897, 109.44 (x2), 167.89, 204.51 (x2),
          256.50 (x2) Hz.  The reference is the Mindlin closed form with rotary inertia and the shear
          factor pi^2/12 (reproduced here to 1e-4).  Criteria (lead decision D-W3-01, superseding the
          requirements 6.3 "8 x 8 NAFEMS mesh" row, which is now informative): 16 x 16 (first four) and
          32 x 32 (all eight, quads and both triangle patterns) within 2 % with the observed order 2.  Thin limit: FV12 and FV16
          (t = 0.05 m) with TSHELL within 2 % on the 32 x 32 mesh and within 0.5 % of the Kirchhoff
          SHELL (DKQ) on the same mesh; FV16 with triangles within 2 %.
* VP-54T  THSHLSTR part of VP-54 (spec 11 section 7 item 7, spec 05d test 6): HOUSE -> STRESS with a
          TSHELL plate under a pure membrane and a pure bending state (synthetic FILE8): the face
          stresses are N/t on both faces (membrane) and +-6 M/t^2 (bending), from the maxima of the
          basic components (D-TSH-01), and the face strains follow plane-stress Hooke's law.
* VP-TS1  TSHELL patch tests and locking: constant membrane, constant bending and constant
          transverse-shear states exact on distorted patches (EINT 0 and 1, quadrilaterals and
          triangles); the constant-shear equilibrium patch test with distributed couples exact on
          parallelogram meshes (MITC4) and on the distorted triangle patch (work-consistent loads incl.
          the condensed rotation bubble); no shear locking: cantilever plate with t/L = 1/1000 within
          1 % of Kirchhoff (quads and triangles), simply supported plate t/a = 1/1000 and 1/10 with both
          triangulation patterns within 2 % (16 x 16) / 0.5 % (32 x 32) of Navier at O(h^2) (the plain
          MITC3 triangle gives 46 %, informative row); exactly six rigid-body modes and no near-zero
          triangle mode; K* = c(beta) K0; distorted quad and triangle meshes converge.
* VP-O1   HOUSE "Optimize Model" (requirements 4.4 item 5, D-HOU-03): an embedded TSHELL basement on
          a layered site, deliberately badly numbered, run SITE -> POINT -> HOUSE -> ANALYS with and
          without the optimizer: the FILE8 transfer functions mapped back through ``<model>.map`` are
          identical (1e-10), the equation bandwidth and profile are reduced, the interaction nodes keep
          their relative (bottom-up) order, and ``<model>.hounew`` reproduces the optimised FILE4.

The NAFEMS 8 x 8 mesh (requirements 6.3 "<= 2 % first 4 modes (NAFEMS mesh)"): TSHELL misses it on mode 4,
the (2, 2) mode with four elements per half-wave: -5.42 % (EINT 0) and -4.09 % (EINT 1).  This is the O(h^2)
discretisation error of the element, not a mass-lumping bias that a converged mesh would hide: against the
element's own kappa = 5/6 Mindlin solution mode 4 is about -5.6 % / -4.3 % on 8 x 8, -1.41 % / -1.08 % on
16 x 16 and -0.35 % / -0.27 % on 32 x 32; a consistent mass gives +1.3 to +6.3 % on 8 x 8 (larger) and a
1/2-1/2 mixed mass still +2.5 % (EINT 1, modes 2-3).  Lead decision D-W3-01 (as D-W1-05 for SHELL) made the
8 x 8 rows informative; the criteria are 2 % on 16 x 16 and 32 x 32 with the observed order ~2.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from sassi.conventions import cfactor
from sassi.core import renumber as RN
from sassi.core import signal as S
from sassi.elements import (ElemRecord, ElementError, assemble, build_dofmap, material_from_E_nu,
                            material_from_M, natural_frequencies)
from sassi.elements import shell, tshell
from sassi.elements.base import blas_quiet, gauss_legendre
from sassi.elements.plane import quad4_shape
from sassi.io import decks
from sassi.io.files import read_container
from sassi.verify import VPResult, problem, worse

# ======================================================================================
# helpers
# ======================================================================================
NAFEMS = dict(E=200e9, nu=0.3, rho=8000.0)
TEST21_REF = [45.897, 109.44, 109.44, 167.89, 204.51, 204.51, 256.50, 256.50]
#: note of the requirements 6.3 rows on the NAFEMS 8 x 8 mesh (see the module docstring)
NAFEMS_8X8_NOTE = ("informative (lead decision D-W3-01): mode 4 (2, 2) misses 2 % on the NAFEMS 8 x 8 mesh by the "
                   "O(h^2) discretisation error of the lumped-mass element (D-CNV-05); 16 x 16 and 32 x 32 are criteria")


def plate_mesh(nx: int, ny: int, a: float, b: float, mat, t: float, code: int = 5, eint: int = 0,
               xy_of: Callable = None) -> Tuple[Dict[int, Tuple[float, float, float]], List[ElemRecord]]:
    """Rectangular nx x ny plate in the X-Y plane (node (i, j) = 1 + i + j (nx + 1)), SHELL (code 3)
    or TSHELL (code 5, ``eint``)."""
    nid = lambda i, j: 1 + i + j * (nx + 1)
    node_xyz = {}
    for j in range(ny + 1):
        for i in range(nx + 1):
            x, y = (a * i / nx, b * j / ny) if xy_of is None else xy_of(i, j)
            node_xyz[nid(i, j)] = (float(x), float(y), 0.0)
    recs = []
    for j in range(ny):
        for i in range(nx):
            props = dict(thick=t) if code == 3 else dict(thick=t, eint=eint)
            recs.append(ElemRecord(1, len(recs) + 1, code, (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)),
                                   False, mat, props))
    return node_xyz, recs


def model(node_xyz, recs, fixfun, undamped: bool = True):
    ids = sorted(node_xyz)
    fix = np.array([fixfun(n, *node_xyz[n]) for n in ids], dtype=int)
    dm = build_dofmap(ids, fix, recs)
    return dm, assemble(node_xyz, recs, dm, undamped=undamped)


def _solve(K, F):
    with blas_quiet():
        return spla.spsolve(sp.csc_matrix(K), F)


def mindlin_ss_frequencies(a: float, b: float, t: float, E: float, nu: float, rho: float,
                           kappa: float = math.pi ** 2 / 12.0, rotary: bool = True, nmodes: int = 8) -> np.ndarray:
    """Closed-form flexural frequencies (Hz) of a simply supported (hard: w = 0 and tangential rotation 0)
    Mindlin plate a x b x t (Navier solution, Mindlin 1951): for each (m, n), k^2 = (m pi/a)^2 + (n pi/b)^2,
    the lower root of ``(D k^2 + kGt - I w^2)(kGt k^2 - rho t w^2) - (kGt)^2 k^2 = 0`` with
    ``I = rho t^3/12`` (0 without rotary inertia)."""
    G = E / (2.0 * (1.0 + nu))
    D = E * t ** 3 / (12.0 * (1.0 - nu * nu))
    kGt = kappa * G * t
    rI = rho * t ** 3 / 12.0 if rotary else 0.0
    out = []
    for m in range(1, 8):
        for n in range(1, 8):
            k2 = (m * math.pi / a) ** 2 + (n * math.pi / b) ** 2
            A, B = rI * rho * t, -(rI * kGt * k2 + rho * t * (D * k2 + kGt))
            C = (D * k2 + kGt) * kGt * k2 - kGt ** 2 * k2
            w2 = -C / B if A == 0 else (-B - math.sqrt(B * B - 4 * A * C)) / (2 * A)
            out.append(math.sqrt(w2) / (2 * math.pi))
    return np.sort(out)[:nmodes]


def test21(n: int, eint: int) -> np.ndarray:
    """NAFEMS Test 21 with an n x n TSHELL mesh: hard simple supports, in-plane restrained."""
    mat = material_from_E_nu(**NAFEMS)
    a = 10.0
    node_xyz, recs = plate_mesh(n, n, a, a, mat, 1.0, eint=eint)

    def fixfun(k, x, y, z):
        ex = abs(x) < 1e-9 or abs(x - a) < 1e-9
        ey = abs(y) < 1e-9 or abs(y - a) < 1e-9
        return [1, 1, 1 if (ex or ey) else 0, 1 if ex else 0, 1 if ey else 0, 0]
    _, am = model(node_xyz, recs, fixfun)
    return natural_frequencies(am.Ks, am.Ms, 8)


def test21_triangles(n: int, diag: str = "right") -> np.ndarray:
    """NAFEMS Test 21 with n x n cells split into two TSHELL triangles each (``diag`` 'right' or 'cross')."""
    node_xyz, recs = triangle_plate_mesh(n, 10.0, material_from_E_nu(**NAFEMS), 1.0, diag)
    _, am = model(node_xyz, recs, _ss_fix(10.0))
    return natural_frequencies(am.Ks, am.Ms, 8)


def fv16_triangles(n: int, diag: str = "right") -> np.ndarray:
    """FV16 cantilevered thin plate (clamped along x = 0) with n x n cells of two TSHELL triangles each."""
    node_xyz, recs = triangle_plate_mesh(n, 10.0, material_from_E_nu(**NAFEMS), 0.05, diag)
    _, am = model(node_xyz, recs, lambda k, x, y, z: [1] * 6 if abs(x) < 1e-9 else [1, 1, 0, 0, 0, 1])
    return natural_frequencies(am.Ks, am.Ms, 6)


def ss_plate_distorted(n: int, eint: int, t: float = 0.05, dist: float = 0.25, seed: int = 3,
                       nmodes: int = 1) -> np.ndarray:
    """Hard simply supported 10 x 10 plate (NAFEMS material) on an n x n TSHELL mesh whose interior nodes are
    moved at random by up to ``dist`` times the element size in x and y (fixed seed)."""
    mat = material_from_E_nu(**NAFEMS)
    a = 10.0
    rng = np.random.default_rng(seed)
    h = a / n
    off = {(i, j): (dist * h * rng.uniform(-1, 1) if 0 < i < n else 0.0, dist * h * rng.uniform(-1, 1) if 0 < j < n
                    else 0.0) for j in range(n + 1) for i in range(n + 1)}
    node_xyz, recs = plate_mesh(n, n, a, a, mat, t, eint=eint,
                                xy_of=lambda i, j: (a * i / n + off[(i, j)][0], a * j / n + off[(i, j)][1]))

    def fixfun(k, x, y, z):
        ex = abs(x) < 1e-9 or abs(x - a) < 1e-9
        ey = abs(y) < 1e-9 or abs(y - a) < 1e-9
        return [1, 1, 1 if (ex or ey) else 0, 1 if ex else 0, 1 if ey else 0, 0]
    _, am = model(node_xyz, recs, fixfun)
    return natural_frequencies(am.Ks, am.Ms, nmodes)


def fv12_tshell(n: int, eint: int) -> np.ndarray:
    """FV12 free thin plate 10 x 10 x 0.05 m with TSHELL, in-plane restrained (u = v = rot_z = 0)."""
    node_xyz, recs = plate_mesh(n, n, 10.0, 10.0, material_from_E_nu(**NAFEMS), 0.05, eint=eint)
    _, am = model(node_xyz, recs, lambda k, x, y, z: [1, 1, 0, 0, 0, 1])
    return natural_frequencies(am.Ks, am.Ms, 9)[3:]


def fv16_tshell(n: int, eint: int) -> np.ndarray:
    """FV16 cantilevered thin plate 10 x 10 x 0.05 m with TSHELL, clamped along x = 0."""
    node_xyz, recs = plate_mesh(n, n, 10.0, 10.0, material_from_E_nu(**NAFEMS), 0.05, eint=eint)
    _, am = model(node_xyz, recs, lambda k, x, y, z: [1] * 6 if abs(x) < 1e-9 else [1, 1, 0, 0, 0, 1])
    return natural_frequencies(am.Ks, am.Ms, 6)


# ======================================================================================
# VP-38T
# ======================================================================================
@problem("VP-38T", "NAFEMS Test 21 thick plate and thin-limit FV12/FV16 (TSHELL)", tier="P2", modules=["HOUSE"],
         source="R2 I.2; requirements 4.1 TSHELL, D-ELM-08, D-W1-05")
def vp38t(workdir):
    from .vp_elements import FV12_REF, FV16_REF, fv12, fv16
    r = VPResult()
    E, nu, rho = NAFEMS["E"], NAFEMS["nu"], NAFEMS["rho"]
    # ---- the reference: Mindlin closed form (rotary inertia, kappa = pi^2/12, hard simple supports)
    cf = mindlin_ss_frequencies(10.0, 10.0, 1.0, E, nu, rho)
    for k, ref in enumerate(TEST21_REF):
        r.check(f"Test 21 closed form (Mindlin, kappa = pi^2/12, rotary inertia) mode {k + 1} vs NAFEMS (Hz)",
                cf[k], ref, rtol=2e-4)
    cf56 = mindlin_ss_frequencies(10.0, 10.0, 1.0, E, nu, rho, kappa=5.0 / 6.0)
    cfnr = mindlin_ss_frequencies(10.0, 10.0, 1.0, E, nu, rho, kappa=5.0 / 6.0, rotary=False)
    r.notes.append("Closed form with the TSHELL shear factor 5/6: " + ", ".join(
        f"{v:.3f} ({100 * (v / q - 1):+.2f} %)" for v, q in zip(cf56, TEST21_REF))
        + "; without rotary inertia: " + ", ".join(f"{100 * (v / q - 1):+.2f} %" for v, q in zip(cfnr, TEST21_REF))
        + " -- hence the Mindlin rotary inertia rho t^3/12 in the TSHELL lumped mass.")
    # ---- TSHELL meshes
    for eint in (0, 1):
        res = {n: test21(n, eint) for n in (8, 16, 32)}
        # D-W3-01 (extends D-W1-05 to TSHELL): the lumped mass of D-CNV-05 converges from below at O(h^2);
        # the coarse NAFEMS 8x8 mesh is reported (informative), the criteria are the 16x16/32x32 meshes and
        # the observed convergence order.
        for k in range(4):
            r.inform(f"Test 21 mode {k + 1} (TSHELL EINT {eint}, NAFEMS 8x8 mesh, Hz)", res[8][k], TEST21_REF[k],
                     note=NAFEMS_8X8_NOTE)
        for k in range(8):
            r.check(f"Test 21 mode {k + 1} (TSHELL EINT {eint}, 32x32, Hz)", res[32][k], TEST21_REF[k], rtol=0.02)
        for k in range(4):
            r.check(f"Test 21 mode {k + 1} (TSHELL EINT {eint}, 16x16, Hz)", res[16][k], TEST21_REF[k], rtol=0.02)
        for k in range(4):
            r.inform(f"Test 21 mode {k + 1} (TSHELL EINT {eint}, 16x16) vs the element's own Mindlin solution "
                     "(kappa = 5/6, Hz)", res[16][k], cf56[k])
        a, b, c = (res[n][0] for n in (8, 16, 32))
        r.check(f"Test 21 mode 1 observed convergence order (TSHELL EINT {eint}, 8/16/32)",
                math.log2(abs(b - a) / abs(c - b)), 2.0, rtol=0.15)
        for n, f in res.items():
            r.notes.append(f"Test 21 TSHELL EINT {eint} {n}x{n}: " + ", ".join(
                f"{100 * (v / q - 1):+.2f}%" for v, q in zip(f, TEST21_REF)))
    # ---- triangles (two per cell, both patterns): Test 21 and FV16 on 32 x 32 cells
    for diag in ("right", "cross"):
        f21 = test21_triangles(32, diag)
        for k in range(8):
            r.check(f"Test 21 mode {k + 1} (TSHELL triangles '{diag}', 32x32 cells, Hz)", f21[k], TEST21_REF[k],
                    rtol=0.02)
        f16 = fv16_triangles(32, diag)
        for k in range(4):
            r.check(f"FV16 mode {k + 1} (TSHELL triangles '{diag}', 32x32 cells, Hz)", f16[k], FV16_REF[k], rtol=0.02)
        r.notes.append(f"triangles '{diag}' 32x32: Test 21 " + ", ".join(f"{100 * (v / q - 1):+.2f}%" for v, q in
                                                                          zip(f21, TEST21_REF))
                       + "; FV16 " + ", ".join(f"{100 * (v / q - 1):+.2f}%" for v, q in zip(f16, FV16_REF)))
    # ---- thin limit: FV12 / FV16 (t/a = 1/200) with TSHELL vs the targets and vs the Kirchhoff SHELL
    s12, s16 = fv12(32), fv16(32)
    for eint in (0, 1):
        t12, t16 = fv12_tshell(32, eint), fv16_tshell(32, eint)
        for name, f, ref, sk in (("FV12", t12, FV12_REF, s12), ("FV16", t16, FV16_REF, s16)):
            for k in range(4):
                r.check(f"{name} mode {k + 1} (TSHELL EINT {eint}, 32x32, Hz)", f[k], ref[k], rtol=0.02)
                r.check(f"{name} mode {k + 1}: TSHELL EINT {eint} / Kirchhoff SHELL (DKQ), same 32x32 mesh", f[k],
                        sk[k], rtol=0.005)
            r.notes.append(f"{name} TSHELL EINT {eint} 32x32: " + ", ".join(
                f"{100 * (v / q - 1):+.2f}%" for v, q in zip(f[:6], ref)) + "; SHELL DKQ: " + ", ".join(
                f"{100 * (v / q - 1):+.2f}%" for v, q in zip(sk[:6], ref)))
    r.notes.append("TSHELL: MITC4 transverse shear (shear factor 5/6), EINT 0 one-point bending + hourglass "
                   f"stabilisation {tshell.HOURGLASS:g}, EINT 1 2x2 bending, lumped mass rho t A/4 + rotary inertia "
                   "rho t^3/12 A/4 on the bending rotations; hard simple supports (w = 0, tangential rotation 0).")
    return r


# ======================================================================================
# VP-TS1: patch tests, locking, rigid-body modes, damping
# ======================================================================================
_MH_OUT = [(0.0, 0.0), (0.24, 0.0), (0.24, 0.12), (0.0, 0.12)]
_MH_IN = [(0.04, 0.02), (0.18, 0.03), (0.16, 0.08), (0.08, 0.08)]
_MH_EL = [(0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7)]
_PATCH_MAT = dict(E=1.0e6, nu=0.25, rho=1.0)
_PATCH_T = 0.01


def _rotation(seed: int = 7) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return Q if np.linalg.det(Q) > 0 else -Q


def mh_patch(eint: int, triangles: bool = False, R: np.ndarray = None):
    """MacNeal-Harder 5-element distorted patch (outer nodes 1-4) of TSHELL quads (or 10 triangles),
    optionally rotated into an oblique plane by ``R``."""
    R = np.eye(3) if R is None else R
    mat = material_from_E_nu(**_PATCH_MAT)
    node_xyz = {i + 1: R @ np.array([p[0], p[1], 0.0]) for i, p in enumerate(_MH_OUT + _MH_IN)}
    recs = []
    for i, c in enumerate(_MH_EL):
        n = tuple(k + 1 for k in c)
        props = dict(thick=_PATCH_T, eint=eint)
        if triangles:
            recs.append(ElemRecord(1, len(recs) + 1, 5, (n[0], n[1], n[2]), False, mat, dict(props)))
            recs.append(ElemRecord(1, len(recs) + 1, 5, (n[0], n[2], n[3]), False, mat, dict(props)))
        else:
            recs.append(ElemRecord(1, len(recs) + 1, 5, n, False, mat, props))
    return node_xyz, recs, mat


def parallelogram_patch(eint: int, skew: float = 0.3, n: int = 3):
    """n x n mesh of equal parallelograms (TSHELL) with interior nodes; returns (node_xyz, recs, interior)."""
    mat = material_from_E_nu(**_PATCH_MAT)
    nid = lambda i, j: 1 + i + j * (n + 1)
    node_xyz = {nid(i, j): np.array([0.1 * i + skew * 0.1 * j, 0.07 * j, 0.0]) for j in range(n + 1) for i in range(n + 1)}
    recs = [ElemRecord(1, k + 1, 5, (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)), False, mat,
                       dict(thick=_PATCH_T, eint=eint)) for k, (j, i) in enumerate((j, i) for j in range(n) for i in range(n))]
    interior = {nid(i, j) for i in range(1, n) for j in range(1, n)}
    return node_xyz, recs, interior, mat


def patch_solve(node_xyz, recs, exact: Callable, free: Callable, loads: Callable = None):
    """Prescribe the exact nodal field on the non-free DOFs, solve the free ones (``free(node, dof)``),
    return (recovered TSHELL resultants (nE, 8) local axes, max relative nodal error)."""
    ids = sorted(node_xyz)
    dm = build_dofmap(ids, np.zeros((len(ids), 6), dtype=int), recs)
    am = assemble(node_xyz, recs, dm)
    u_ex = np.array([exact(np.asarray(node_xyz[int(n)]))[int(d) - 1] for n, d in zip(dm.eq_node, dm.eq_dof)])
    fr = np.array([bool(free(int(n), int(d))) for n, d in zip(dm.eq_node, dm.eq_dof)])
    K = am.Ks.real.tocsr()
    f = np.zeros(dm.neq) if loads is None else loads(dm)
    u = u_ex.copy()
    if fr.any():
        u[fr] = _solve(K[fr][:, fr], f[fr] - K[fr][:, ~fr] @ u_ex[~fr])
    rec = am.recovery["TSHELL"]
    ue = np.where(rec["eq"] >= 0, u[np.maximum(rec["eq"], 0)], 0.0)
    s = np.einsum("ecd,ed->ec", rec["S"].real, ue)
    return s, float(np.abs(u - u_ex).max() / max(np.abs(u_ex).max(), 1e-300))


def to_patch_axes(node_xyz, recs, s: np.ndarray, R: np.ndarray = None) -> np.ndarray:
    """Resultants rotated from each element's local axes x'y' to the patch axes: rows
    [Nxx Nyy Nxy Qx Qy Mxx Myy Mxy]."""
    R = np.eye(3) if R is None else R
    out = []
    for k, r in enumerate(recs):
        X = np.array([np.asarray(node_xyz[n]) for n in dict.fromkeys(r.nodes)])
        Lam, _, _ = shell.local_frame(X)
        A = (Lam @ R)[:2, :2]                         # local x', y' in patch coordinates
        N = A.T @ np.array([[s[k, 0], s[k, 2]], [s[k, 2], s[k, 1]]]) @ A
        M = A.T @ np.array([[s[k, 5], s[k, 7]], [s[k, 7], s[k, 6]]]) @ A
        Q = A.T @ s[k, 3:5]
        out.append([N[0, 0], N[1, 1], N[0, 1], Q[0], Q[1], M[0, 0], M[1, 1], M[0, 1]])
    return np.array(out)


def couple_loads(node_xyz, recs, m: np.ndarray) -> Callable:
    """Consistent nodal loads of a uniform distributed couple ``m = (m_x, m_y)`` (moment per unit area
    acting on beta_x, beta_y; theta_y = beta_x, theta_x = -beta_y) on TSHELL quadrilaterals in the X-Y plane:
    ``f = integral(N_i m) dA``.  Triangles also load their condensed bubble rotation: see
    :func:`triangle_shear_patch`."""
    def loads(dm):
        f = np.zeros(dm.neq)
        x, w = gauss_legendre(3)
        for r in recs:
            ns = list(dict.fromkeys(r.nodes))
            if len(ns) != 4:
                raise ValueError("couple_loads: quadrilaterals only (triangles: triangle_shear_patch)")
            X = np.array([np.asarray(node_xyz[n])[:2] for n in ns])
            for xa, wa in zip(x, w):
                for xb, wb in zip(x, w):
                    N, dN = quad4_shape(np.array([[xa, xb]]))
                    dA = wa * wb * np.linalg.det(dN[0].T @ X)
                    for i, nd in enumerate(ns):
                        f[dm.eq(nd, 5)] += N[0, i] * m[0] * dA
                        f[dm.eq(nd, 4)] -= N[0, i] * m[1] * dA
        return f
    return loads


def triangle_shear_patch(t: float, gam: Sequence[float] = (1.0e-3, -2.0e-3)) -> Tuple[float, float, float]:
    """Constant transverse shear ``gamma`` (w linear, beta = 0) on the MacNeal-Harder patch split into 10
    TSHELL triangles in the X-Y plane (E, nu of the patch material, thickness ``t``).

    Returns ``(op_err, q_err, u_err)``: the error of the triangle strain operators applied to the exact
    field (bubble rotation 0), and the Q and nodal errors of the equilibrium patch test: uniform distributed
    couple ``m = Q = D_s gamma`` (moment equilibrium Q = div M + m with M = 0), boundary nodes prescribed,
    interior nodes free, work-consistent element loads (``A/3 m`` on each corner rotation, ``27/60 A m`` on
    the bubble rotation, condensed with the bubble)."""
    E, nu = _PATCH_MAT["E"], _PATCH_MAT["nu"]
    _, Db, Ds = tshell.plate_rigidities(E, nu, t)
    gam = np.asarray(gam, dtype=float)
    m = Ds @ gam
    xy = np.array(_MH_OUT + _MH_IN)
    tris = []
    for c in _MH_EL:
        tris += [(c[0], c[1], c[2]), (c[0], c[2], c[3])]
    nn = len(xy)
    K = np.zeros((3 * nn, 3 * nn))
    F = np.zeros(3 * nn)
    u_ex = np.zeros(3 * nn)
    u_ex[0::3] = xy @ gam
    op_err, kept = 0.0, []
    for tri in tris:
        X = xy[list(tri)]
        op = tshell.tri_plate_operators(X, Db, Ds, t)
        dofs = np.array([3 * k + d for k in tri for d in range(3)])
        op_err = worse(op_err, float(np.abs(Ds @ op["Bs"][:, :9] @ u_ex[dofs] - m).max() / np.abs(m).max()))
        A = shell.polygon_area(X)
        f = np.zeros(11)
        f[2:9:3], f[1:9:3] = A / 3.0 * m[0], -A / 3.0 * m[1]          # theta_y = beta_x, theta_x = -beta_y
        f[9:] = tshell.TRI_BUBBLE_MEAN * A * m
        Kf = op["K"]
        Xc = np.linalg.solve(Kf[9:, 9:], np.column_stack([Kf[9:, :9], f[9:]]))   # bubble = Xc[:, 9] - Xc[:, :9] u
        K[np.ix_(dofs, dofs)] += Kf[:9, :9] - Kf[:9, 9:] @ Xc[:, :9]
        F[dofs] += f[:9] - Kf[:9, 9:] @ Xc[:, 9]
        kept.append((dofs, op, Xc))
    free = np.repeat(np.arange(nn) >= 4, 3)
    u = u_ex.copy()
    u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free] - K[np.ix_(free, ~free)] @ u_ex[~free])
    q_err = 0.0
    for dofs, op, Xc in kept:
        bubble = Xc[:, 9] - Xc[:, :9] @ u[dofs]
        Q = Ds @ (op["Bs"][:, :9] @ u[dofs] + op["Bs"][:, 9:] @ bubble)
        q_err = worse(q_err, float(np.abs(Q - m).max() / np.abs(m).max()))
    return op_err, q_err, float(np.abs(u - u_ex).max() / np.abs(u_ex).max())


def triangle_plate_mesh(n: int, a: float, mat, t: float, diag: str = "right", dist: float = 0.0, seed: int = 3):
    """n x n cells of the square [0, a]^2 in the X-Y plane, each split into two TSHELL triangles: ``diag``
    'right' (every diagonal from (i, j) to (i+1, j+1)) or 'cross' (alternating); interior nodes moved at
    random by up to ``dist`` x the cell size (fixed seed)."""
    rng = np.random.default_rng(seed)
    h = a / n
    nid = lambda i, j: 1 + i + j * (n + 1)
    node_xyz = {}
    for j in range(n + 1):
        for i in range(n + 1):
            dx = dist * h * rng.uniform(-1, 1) if 0 < i < n else 0.0
            dy = dist * h * rng.uniform(-1, 1) if 0 < j < n else 0.0
            node_xyz[nid(i, j)] = (a * i / n + dx, a * j / n + dy, 0.0)
    recs = []
    for j in range(n):
        for i in range(n):
            c = (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1))
            pair = ([(c[0], c[1], c[2]), (c[0], c[2], c[3])] if diag == "right" or (i + j) % 2 == 0
                    else [(c[0], c[1], c[3]), (c[1], c[2], c[3])])
            for nodes in pair:
                recs.append(ElemRecord(1, len(recs) + 1, 5, nodes, False, mat, dict(thick=t, eint=0)))
    return node_xyz, recs


def _ss_fix(a: float):
    """Hard simple supports of the square [0, a]^2: w = 0 and the tangential rotation fixed on the edges,
    in-plane restrained."""
    def fixfun(k, x, y, z):
        ex = abs(x) < 1e-9 or abs(x - a) < 1e-9
        ey = abs(y) < 1e-9 or abs(y - a) < 1e-9
        return [1, 1, 1 if (ex or ey) else 0, 1 if ex else 0, 1 if ey else 0, 0]
    return fixfun


def navier_ss_centre(a: float, t: float, E: float, nu: float, q: float = 1.0, kappa: float = 5.0 / 6.0,
                     nterms: int = 199) -> float:
    """Centre deflection of the hard simply supported square Mindlin plate under a uniform load q (Navier):
    ``w_mn = q_mn (1/(D k^4) + 1/(kappa G t k^2))``, ``q_mn = 16 q/(pi^2 m n)``, m, n odd.  The shear factor is
    the literal 5/6 of the Reissner-Mindlin theory (not the element constant: final audit)."""
    D = E * t ** 3 / (12.0 * (1.0 - nu * nu))
    G = E / (2.0 * (1.0 + nu))
    m = np.arange(1, nterms + 1, 2, dtype=float)
    M, N = np.meshgrid(m, m, indexing="ij")
    k2 = (M * math.pi / a) ** 2 + (N * math.pi / a) ** 2
    s = np.sin(M * math.pi / 2) * np.sin(N * math.pi / 2)
    return float(np.sum(16.0 * q / (math.pi ** 2 * M * N) * (1.0 / (D * k2 * k2) + 1.0 / (kappa * G * t * k2)) * s))


def ss_plate_centre(node_xyz, recs, a: float) -> float:
    """Centre deflection of a hard simply supported square plate mesh under a unit uniform load lumped to the
    nodes by tributary area (regular n x n node grid)."""
    dm, am = model(node_xyz, recs, _ss_fix(a))
    xs = sorted({round(p[0], 9) for p in node_xyz.values() if abs(p[1]) < 1e-9})
    h = xs[1] - xs[0]
    F = np.zeros(dm.neq)
    for k, (x, y, _) in node_xyz.items():
        e = dm.eq(k, 3)
        if e >= 0:
            F[e] += (h if 1e-9 < x < a - 1e-9 else h / 2) * (h if 1e-9 < y < a - 1e-9 else h / 2)
    u = _solve(am.Ks.real, F)
    centre = min(node_xyz, key=lambda k: (node_xyz[k][0] - a / 2) ** 2 + (node_xyz[k][1] - a / 2) ** 2)
    return float(u[dm.eq(centre, 3)])


def _locking_reference(nx: int, ny: int, L: float, b: float, t: float, mat) -> float:
    """Tip deflection of the cantilever plate with *plain* 2x2 Gauss integration of the displacement-based
    transverse shear (the locking element the MITC4 field replaces), for the informative comparison."""
    node_xyz, recs = plate_mesh(nx, ny, L, b, mat, t, eint=1)
    dm, am = model(node_xyz, recs, lambda k, x, y, z: [1] * 6 if abs(x) < 1e-9 else [0] * 6)
    K = am.Ks.real.tolil()
    E, nu = mat.E0, mat.nu0
    _, _, Ds = tshell.plate_rigidities(E, nu, t)
    gp, gw = gauss_legendre(2)
    for r in recs:
        xy = np.array([node_xyz[n][:2] for n in r.nodes])
        op = tshell.mitc4_shear_operator(xy)
        dK = np.zeros((12, 12))
        for a_, wa in zip(gp, gw):
            for b_, wb in zip(gp, gw):
                p = np.array([a_, b_])
                N, dN, J = tshell._jacobian(xy, p)
                Bd = np.linalg.solve(J, tshell._covariant_shear(N, dN, J))
                Bm, det = op(p)
                dK += wa * wb * det * (Bd.T @ Ds @ Bd - Bm.T @ Ds @ Bm)
        eq = dm.eqs(r.nodes, (3, 4, 5))
        for i in range(12):
            for j in range(12):
                if eq[i] >= 0 and eq[j] >= 0:
                    K[eq[i], eq[j]] += dK[i, j]
    tip = [n for n, p in node_xyz.items() if abs(p[0] - L) < 1e-9]
    F = _tip_load(dm, node_xyz, tip, b, ny)
    u = _solve(K.tocsr(), F)
    return float(np.mean([u[dm.eq(n, 3)] for n in tip]))


def _tip_load(dm, node_xyz, tip, b: float, ny: int) -> np.ndarray:
    """Unit total transverse load spread consistently along the free edge."""
    F = np.zeros(dm.neq)
    for n in tip:
        y = node_xyz[n][1]
        F[dm.eq(n, 3)] += (b / ny) * (0.5 if (abs(y) < 1e-9 or abs(y - b) < 1e-9) else 1.0) / b
    return F


def cantilever_tip(nx: int, ny: int, L: float, b: float, t: float, eint: int, mat) -> float:
    node_xyz, recs = plate_mesh(nx, ny, L, b, mat, t, eint=eint)
    dm, am = model(node_xyz, recs, lambda k, x, y, z: [1] * 6 if abs(x) < 1e-9 else [0] * 6)
    tip = [n for n, p in node_xyz.items() if abs(p[0] - L) < 1e-9]
    u = _solve(am.Ks.real, _tip_load(dm, node_xyz, tip, b, ny))
    return float(np.mean([u[dm.eq(n, 3)] for n in tip]))


def cantilever_tip_triangles(nx: int, ny: int, L: float, b: float, t: float, mat, diag: str = "right") -> float:
    """Tip deflection of the cantilever plate of :func:`cantilever_tip` with every cell split into two TSHELL
    triangles (``diag`` as :func:`triangle_plate_mesh`)."""
    node_xyz, recs = plate_mesh(nx, ny, L, b, mat, t)
    out = []
    for rq in recs:
        c = rq.nodes
        i = (rq.id - 1) % nx
        j = (rq.id - 1) // nx
        pair = ([(c[0], c[1], c[2]), (c[0], c[2], c[3])] if diag == "right" or (i + j) % 2 == 0
                else [(c[0], c[1], c[3]), (c[1], c[2], c[3])])
        out += [ElemRecord(1, len(out) + k + 1, 5, nodes, False, mat, dict(thick=t)) for k, nodes in enumerate(pair)]
    dm, am = model(node_xyz, out, lambda k, x, y, z: [1] * 6 if abs(x) < 1e-9 else [0] * 6)
    tip = [n for n, p in node_xyz.items() if abs(p[0] - L) < 1e-9]
    u = _solve(am.Ks.real, _tip_load(dm, node_xyz, tip, b, ny))
    return float(np.mean([u[dm.eq(n, 3)] for n in tip]))


def plain_mitc3_ss_centre(n: int, t: float, mat, a: float = 10.0) -> float:
    """Centre deflection of the simply supported plate of :func:`ss_plate_centre` ('right' triangles) with the
    *plain* MITC3 triangle: the bubble rotation fixed at 0 and the full weight on the linear shear part
    (``tri_plate_operators(..., alpha=0)``, corner block) -- the locking element, for the informative row."""
    node_xyz, recs = triangle_plate_mesh(n, a, mat, t, "right")
    ids = sorted(node_xyz)
    pos = {k: i for i, k in enumerate(ids)}
    _, Db, Ds = tshell.plate_rigidities(mat.E0, mat.nu0, t)
    K = sp.lil_matrix((3 * len(ids), 3 * len(ids)))
    for rq in recs:
        xy = np.array([node_xyz[k][:2] for k in rq.nodes])
        Ke = tshell.tri_plate_operators(xy, Db, Ds, t, alpha=0.0)["K"][:9, :9]
        dofs = [3 * pos[k] + d for k in rq.nodes for d in range(3)]
        for i, di in enumerate(dofs):
            for j, dj in enumerate(dofs):
                K[di, dj] += Ke[i, j]
    fix = _ss_fix(a)
    keep = np.array([not fix(k, *node_xyz[k])[2 + d] for k in ids for d in range(3)])
    h = a / n
    F = np.zeros(3 * len(ids))
    for k, (x, y, _) in node_xyz.items():
        F[3 * pos[k]] = (h if 1e-9 < x < a - 1e-9 else h / 2) * (h if 1e-9 < y < a - 1e-9 else h / 2)
    Kc = K.tocsr()[keep][:, keep]
    u = np.zeros(3 * len(ids))
    u[keep] = _solve(Kc, F[keep])
    centre = min(node_xyz, key=lambda k: (node_xyz[k][0] - a / 2) ** 2 + (node_xyz[k][1] - a / 2) ** 2)
    return float(u[3 * pos[centre]])


@problem("VP-TS1", "TSHELL patch tests, shear locking, rigid-body modes", tier="P2", modules=["HOUSE"],
         source="requirements 1.5/4.1 TSHELL, D-ELM-08; MacNeal & Harder (1985) patch")
def vpts1(workdir):
    r = VPResult()
    E, nu, t = _PATCH_MAT["E"], _PATCH_MAT["nu"], _PATCH_T
    Dm, Db, Ds = tshell.plate_rigidities(E, nu, t)
    R = _rotation()
    # ---- 1. constant membrane state with a rigid in-plane rotation (drilling consistency), oblique plane
    Om = 2.0e-4
    eps = 1.0e-3 * np.array([1.0, 1.0, 1.0])                   # u = 1e-3 (x + y/2), v = 1e-3 (y + x/2)
    n_exact = Dm @ eps

    def mem_field(p):
        q = R.T @ p
        loc = np.array([1e-3 * (q[0] + 0.5 * q[1]) - Om * q[1], 1e-3 * (q[1] + 0.5 * q[0]) + Om * q[0], 0.0])
        return np.concatenate([R @ loc, Om * R[:, 2]])          # theta = Om z' (drilling)
    for eint in (0, 1):
        for tri in (False, True):
            node_xyz, recs, _ = mh_patch(eint, tri, R)
            s, du = patch_solve(node_xyz, recs, mem_field, lambda n, d: n > 4)
            g = to_patch_axes(node_xyz, recs, s, R)
            what = "CST triangles" if tri else "Q4 + incompatible modes"
            r.check(f"membrane patch ({what}, EINT {eint}, oblique plane, drilling free): max N error / |N|",
                    np.abs(g[:, :3] - n_exact).max() / np.abs(n_exact).max(), 0.0, atol=1e-10)
            r.check(f"membrane patch ({what}, EINT {eint}): nodal field incl. theta_z = omega exact", du, 0.0,
                    atol=1e-8, note="all 6 DOFs of the interior nodes free: cond(K_ff) ~ 1e7-1e8 (drilling penalty "
                                    "1e-4, thin-plate bending), so nodal values carry ~cond x eps round-off")
            r.check(f"membrane patch ({what}, EINT {eint}): no bending or shear", np.abs(g[:, 3:]).max()
                    / np.abs(n_exact).max(), 0.0, atol=1e-10)
    # ---- 2. constant bending (w quadratic, beta = -grad w): exact M, Q = 0
    k1, k2, k12 = 1e-3, 2e-3, 0.5e-3
    m_exact = Db @ np.array([-k1, -k2, -2.0 * k12])

    def bend_field(p):
        x, y = p[0], p[1]
        wx, wy = k1 * x + k12 * y, k2 * y + k12 * x
        return np.array([0.0, 0.0, 0.5 * k1 * x * x + 0.5 * k2 * y * y + k12 * x * y, wy, -wx, 0.0])
    for eint in (0, 1):
        for tri in (False, True):
            node_xyz, recs, _ = mh_patch(eint, tri)
            s, du = patch_solve(node_xyz, recs, bend_field, lambda n, d: n > 4 and d in (3, 4, 5))
            g = to_patch_axes(node_xyz, recs, s)
            what = "triangles" if tri else "MITC4"
            r.check(f"bending patch ({what}, EINT {eint}): max M error / |M|", np.abs(g[:, 5:] - m_exact).max()
                    / np.abs(m_exact).max(), 0.0, atol=1e-10)
            r.check(f"bending patch ({what}, EINT {eint}): transverse shear Q = 0 (relative to |M|/patch size)",
                    np.abs(g[:, 3:5]).max() * 0.24 / np.abs(m_exact).max(), 0.0, atol=1e-10)
            r.check(f"bending patch ({what}, EINT {eint}): interior nodal values exact", du, 0.0, atol=1e-10)
    # ---- 3. constant transverse shear (w linear, beta = 0): Q = kappa G t gamma
    gam = np.array([1.0e-3, -2.0e-3])
    q_exact = Ds @ gam
    shear_field = lambda p: np.array([0.0, 0.0, gam[0] * p[0] + gam[1] * p[1], 0.0, 0.0, 0.0])
    for eint in (0, 1):
        node_xyz, recs, _ = mh_patch(eint, False)
        s, _ = patch_solve(node_xyz, recs, shear_field, lambda n, d: False)
        g = to_patch_axes(node_xyz, recs, s)
        r.check(f"constant shear state reproduced in every distorted element (MITC4, EINT {eint}): max Q error",
                np.abs(g[:, 3:5] - q_exact).max() / np.abs(q_exact).max(), 0.0, atol=1e-10)
        # equilibrium patch test: distributed couple m = Q balances the constant shear (M = 0)
        s, du = patch_solve(node_xyz, recs, shear_field, lambda n, d: n > 4 and d in (3, 4, 5),
                            couple_loads(node_xyz, recs, q_exact))
        g = to_patch_axes(node_xyz, recs, s)
        r.inform(f"shear equilibrium patch, distributed couples (MITC4, distorted MacNeal-Harder, EINT {eint}): "
                 "max Q error", np.abs(g[:, 3:5] - q_exact).max() / np.abs(q_exact).max(), 0.0,
                 note="the MITC4 field reproduces constant shear exactly but its element integral differs from "
                      "the displacement-based one on non-parallelograms")
        # ... and on parallelograms the equilibrium patch is exact (final audit: this block ran inside the
        # triangle loop below with the leaked loop variable, so EINT 1 was recorded three times and EINT 0
        # never; tests/unit/test_audit_vpsuite.py)
        node_xyz, recs, interior, _ = parallelogram_patch(eint)
        s, du = patch_solve(node_xyz, recs, shear_field, lambda n, d: n in interior and d in (3, 4, 5),
                            couple_loads(node_xyz, recs, q_exact))
        g = to_patch_axes(node_xyz, recs, s)
        r.check(f"shear equilibrium patch, distributed couples (MITC4, skewed parallelograms, EINT {eint}): max Q error",
                np.abs(g[:, 3:5] - q_exact).max() / np.abs(q_exact).max(), 0.0, atol=1e-10)
        r.check(f"shear equilibrium patch (MITC4, skewed parallelograms, EINT {eint}): interior nodal values exact",
                du, 0.0, atol=1e-10)
    # triangles (EINT has no effect): the condensed bubble rotation is an element DOF, so the constant-shear
    # state is checked with the exact field (bubble 0) and as an equilibrium state with work-consistent loads
    for tt in (t, 1.0e-4, 0.05):
        op_err, q_err, u_err = triangle_shear_patch(tt, gam)
        r.check(f"constant shear state reproduced by the triangle strain operators (distorted patch, t = {tt:g}, "
                "bubble rotation 0): max Q error", op_err, 0.0, atol=1e-10)
        r.check(f"shear equilibrium patch, distributed couples (triangles, distorted, t = {tt:g}): max Q error",
                q_err, 0.0, atol=1e-10)
        r.check(f"shear equilibrium patch (triangles, distorted, t = {tt:g}): interior nodal values exact", u_err,
                0.0, atol=1e-8, note="cond(K_ff) grows like (h/t)^2 (1e6 at t = 1e-4): nodal values carry ~cond x eps "
                                     "round-off, the resultants (row above) do not")
    # ---- 4. no shear locking: cantilever plate L = 10, b = 2, nu = 0 (Kirchhoff = beam), tip line load 1
    E0 = 200e9
    mat0 = material_from_E_nu(E0, 0.0, 8000.0)
    L, b = 10.0, 2.0
    for tL in (1e-3, 1e-1):
        tt = tL * L
        I = b * tt ** 3 / 12.0
        ref_k = L ** 3 / (3 * E0 * I)                                   # Kirchhoff (Euler-Bernoulli)
        ref_m = ref_k + L / (5.0 / 6.0 * E0 / 2.0 * b * tt)               # + Mindlin shear deflection
        for eint in (0, 1):
            d8 = cantilever_tip(8, 2, L, b, tt, eint, mat0)
            if tL == 1e-3:
                r.check(f"cantilever plate t/L = 1/1000 (TSHELL EINT {eint}, 8x2): tip deflection / Kirchhoff", d8 / ref_k,
                        1.0, rtol=0.01)
                d32 = cantilever_tip(32, 4, L, b, tt, eint, mat0)
                r.check(f"cantilever plate t/L = 1/1000 (TSHELL EINT {eint}, 32x4): tip deflection / Kirchhoff",
                        d32 / ref_k, 1.0, rtol=0.001)
            else:
                r.check(f"cantilever plate t/L = 1/10 (TSHELL EINT {eint}, 8x2): tip deflection / Timoshenko",
                        d8 / ref_m, 1.0, rtol=0.01)
            if tL == 1e-3:
                r.inform(f"cantilever t/L = 1/1000 (EINT {eint}, 8x2) vs the discrete result 1 - 1/(4 N^2) of one-point "
                         "shear Timoshenko elements (N = 8)", d8 / ref_k, 1.0 - 1.0 / 256.0)
    lock = _locking_reference(8, 2, L, b, 1e-3 * L, mat0)
    r.inform("plain 2x2 displacement-based shear (locking element), t/L = 1/1000, 8x2: tip deflection / Kirchhoff",
             lock / (L ** 3 / (3 * E0 * b * (1e-2) ** 3 / 12.0)), 1.0, note="shear locking: the MITC4 field removes it")
    # triangles: cantilever (16 x 4 cells, two triangles each) and the simply supported plate (Navier, kappa 5/6)
    tt = 1e-3 * L
    for diag in ("right", "cross"):
        d16 = cantilever_tip_triangles(16, 4, L, b, tt, mat0, diag)
        r.check(f"cantilever plate t/L = 1/1000 (TSHELL triangles '{diag}', 16x4 cells): tip deflection / Kirchhoff",
                d16 / (L ** 3 / (3 * E0 * b * tt ** 3 / 12.0)), 1.0, rtol=0.01)
    matn = material_from_E_nu(**NAFEMS)
    for ta in (1e-3, 0.1):
        ref = navier_ss_centre(10.0, ta * 10.0, NAFEMS["E"], NAFEMS["nu"])
        for diag in ("right", "cross"):
            w = {n: ss_plate_centre(*triangle_plate_mesh(n, 10.0, matn, ta * 10.0, diag), 10.0) / ref for n in (8, 16, 32)}
            r.check(f"SS plate t/a = {ta:g}, uniform load (triangles '{diag}', 16x16 cells): centre w / Navier", w[16],
                    1.0, rtol=0.02)
            r.check(f"SS plate t/a = {ta:g}, uniform load (triangles '{diag}', 32x32 cells): centre w / Navier", w[32],
                    1.0, rtol=0.005)
            e8, e16, e32 = (abs(w[n] - 1.0) for n in (8, 16, 32))
            r.check(f"SS plate t/a = {ta:g} (triangles '{diag}'): observed order 8/16/32", math.log2(e16 / e32), 2.0,
                    rtol=0.15)
            r.notes.append(f"SS plate t/a = {ta:g} triangles '{diag}' 8/16/32: "
                           + ", ".join(f"{w[n]:.5f}" for n in (8, 16, 32)) + " x Navier")
    ref = navier_ss_centre(10.0, 0.01, NAFEMS["E"], NAFEMS["nu"])
    r.inform("plain MITC3 triangle (no bubble, unweighted linear shear), SS plate t/a = 1/1000, 16x16 'right': centre "
             "w / Navier", plain_mitc3_ss_centre(16, 0.01, matn) / ref, 1.0,
             note="shear locking of the plain MITC3 triangle (as many shear constraints as DOFs): the condensed "
                  "rotation bubble and the weighted linear shear part remove it")
    # ---- 5. rigid-body modes: exactly 6 zero eigenvalues per free element
    mat = material_from_E_nu(**NAFEMS)
    quad = np.array([[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [1.0, 0.9, 0.0], [-0.1, 1.1, 0.0]]) @ R.T
    for label, xyz in (("distorted quad", quad), ("triangle", quad[:3]), ("triangle given as I J K K",
                                                                            quad[[0, 1, 2, 2]])):
        for eint in (0, 1):
            for tt in (0.01, 0.5):
                K, _ = tshell.matrices(xyz, mat, thick=tt, eint=eint)
                w = np.linalg.eigvalsh(0.5 * (K.real + K.real.T))
                nz = int(np.sum(np.abs(w) < 1e-9 * np.abs(w).max()))
                expect = 12 if xyz.shape[0] == 4 and len(set(map(tuple, np.round(xyz, 12)))) == 3 else 6
                r.check(f"{label}, t = {tt:g}, EINT {eint}: zero-energy modes (6 rigid-body"
                        + (" + 6 of the repeated node)" if expect == 12 else ")"), nz, expect, atol=0)
                if label == "triangle" and eint == 0:
                    ratio = float(w[6] / np.abs(w).max())
                    r.require(f"triangle, t = {tt:g}: no near-zero-energy mode (7th eigenvalue / max = {ratio:.2e} "
                              "> 1e-7, the drilling-penalty level)", ratio > 1e-7)
    # free 4x4 plate in an oblique plane: six rigid-body modes, no spurious mode
    node_xyz, recs = plate_mesh(4, 4, 2.0, 2.0, mat, 0.2)
    node_xyz = {k: tuple(R @ np.asarray(v)) for k, v in node_xyz.items()}
    for eint in (0, 1):
        for rr in recs:
            rr.props["eint"] = eint
        _, am = model(node_xyz, recs, lambda *a: [0] * 6)
        K = am.Ks.toarray().real
        w = np.linalg.eigvalsh(0.5 * (K + K.T))
        r.check(f"free 4x4 TSHELL plate (oblique plane, EINT {eint}): zero-energy modes", int(np.sum(np.abs(w) < 1e-9
                * np.abs(w).max())), 6, atol=0)
    # ---- 6. complex stiffness K* = c(beta) K0 (D-CNV-04) and recovery scaled alike
    bd = 0.06
    md = material_from_M(1, 3.0e7, 0.2, 24.0, bd, bd, 9.81)
    for eint in (0, 1):
        Kc, _ = tshell.matrices(quad, md, thick=0.3, eint=eint)
        K0, _ = tshell.matrices(quad, md.undamped(), thick=0.3, eint=eint)
        r.check(f"TSHELL EINT {eint}: max|K* - c(beta) K0| / max|K0|", np.abs(Kc - cfactor(bd) * K0.real).max()
                / np.abs(K0).max(), 0.0, atol=1e-12)
        Sc = tshell.recovery(quad, md, thick=0.3, eint=eint)
        S0 = tshell.recovery(quad, md.undamped(), thick=0.3, eint=eint)
        r.check(f"TSHELL EINT {eint}: max|S* - c(beta) S0| / max|S0|", np.abs(Sc - cfactor(bd) * S0.real).max()
                / np.abs(S0).max(), 0.0, atol=1e-12)
    # ---- 7. distorted meshes converge without locking: thin (t/a = 1/200) simply supported plate, interior nodes
    #         moved at random by up to 25 % of h; reference = Mindlin closed form with the element's kappa = 5/6
    ref1 = mindlin_ss_frequencies(10.0, 10.0, 0.05, NAFEMS["E"], NAFEMS["nu"], NAFEMS["rho"], kappa=5.0 / 6.0,
                                  nmodes=1)[0]
    for eint in (0, 1):
        f1 = [ss_plate_distorted(n, eint)[0] for n in (8, 16, 32)]
        errs = [abs(f / ref1 - 1.0) for f in f1]
        r.require(f"distorted thin SS plate (EINT {eint}): mode-1 error decreases 8 -> 16 -> 32 "
                  f"({', '.join(f'{100 * e:.3f} %' for e in errs)})", errs[0] > errs[1] > errs[2])
        r.check(f"distorted thin SS plate (EINT {eint}, 32x32, 25 % random distortion): mode 1 (Hz)", f1[2], ref1,
                rtol=0.005)
        r.inform(f"distorted thin SS plate (EINT {eint}): observed order 8/16/32 (random distortion per mesh)",
                 math.log2(errs[0] / errs[1]) if errs[1] > 0 else float("nan"), 2.0)
    f1 = []
    for n in (8, 16, 32):
        node_xyz, recs = triangle_plate_mesh(n, 10.0, mat, 0.05, "cross", dist=0.25)
        _, am = model(node_xyz, recs, _ss_fix(10.0))
        f1.append(natural_frequencies(am.Ks, am.Ms, 1)[0])
    errs = [abs(f / ref1 - 1.0) for f in f1]
    r.require(f"distorted thin SS plate (triangles 'cross'): mode-1 error decreases 8 -> 16 -> 32 "
              f"({', '.join(f'{100 * e:.3f} %' for e in errs)})", errs[0] > errs[1] > errs[2])
    r.check("distorted thin SS plate (triangles 'cross', 32x32 cells, 25 % random distortion): mode 1 (Hz)", f1[2],
            ref1, rtol=0.005)
    r.notes.append("Patch: MacNeal-Harder 5-element distorted patch (0.24 x 0.12), E = 1e6, nu = 0.25, t = 0.01; "
                   "membrane test rotated into an oblique plane with a rigid in-plane rotation (theta_z = omega).")
    r.notes.append("The MITC4 constant-shear equilibrium patch test with distributed couples is exact on parallelogram "
                   "meshes; on the distorted patch it is not (informative row), a known property of the assumed covariant "
                   "shear interpolation; constant shear itself is reproduced exactly in every distorted element and distorted "
                   "meshes converge without locking (section 7).")
    r.notes.append("Triangles: MITC3 shear of the corner DOFs + condensed cubic rotation bubble (seen by the shear through "
                   f"its element mean 27/60) + weight t^2/(t^2 + {tshell.TRI_STAB_ALPHA:g} h^2) on the linear shear part "
                   "(sassi.elements.tshell.tri_plate_operators); EINT has no effect on triangles.")
    return r


# ======================================================================================
# VP-54T: THSHLSTR face stresses through HOUSE -> STRESS
# ======================================================================================
NFFT54, DT54 = 512, 0.01
DF54 = 1.0 / (NFFT54 * DT54)
FNUM54 = np.array([2, 4, 6, 8, 10, 12])
K0_54, A0_54 = 7, 0.1                     # harmonic control motion: bin 7, 0.1 g
T54 = 0.4
MAT54 = (1, 1, 3.0e7, 0.2, 24.0, 0.0, 0.0)        # M row: undamped (real STF: N(t) = N0 u_g(t))


def tshell_plate_deck(eint: int = 0, nx: int = 2, ny: int = 2, a: float = 2.0, b: float = 1.5):
    """HOUSE deck of a flat TSHELL plate in an inclined plane (local x' = e1, y' = e2, z' = e1 x e2 for
    every element) and the frame (origin, e1, e2)."""
    from .vp_house import add_element, add_node, house_deck
    d = house_deck()
    o = np.array([0.3, -0.2, 0.5])
    e1 = np.array([2.0, 1.0, 2.0]) / 3.0
    v = np.array([-1.0, 2.0, 0.5])
    e2 = v - (v @ e1) * e1
    e2 /= np.linalg.norm(e2)
    nid = lambda i, j: 1 + i + j * (nx + 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            add_node(d, nid(i, j), *(o + a * i / nx * e1 + b * j / ny * e2))
    d.table("materials").append(list(MAT54))
    d.table("groups").append([1, 5, "thick plate"])
    for j in range(ny):
        for i in range(nx):
            add_element(d, 1, 1 + i + j * nx, [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)], etype=1,
                        mat=1, eint=eint, thick=T54)
    return d, (o, e1, e2)


def read_face_file(path: Path) -> Dict[Tuple[int, int, str], np.ndarray]:
    out = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split()
        out[(int(p[0]), int(p[1]), p[2])] = np.array([float(x) for x in p[3:]])
    return out


@problem("VP-54T", "THSHLSTR TSHELL face stresses: pure membrane N/t, pure bending +-6M/t^2", tier="P2",
         modules=["HOUSE", "STRESS"], source="spec 11 sections 4.1 and 7 item 7; spec 05d test 6; D-TSH-01")
def vp54t(workdir):
    from .vp_house import run_house
    from .vp_stress import add_eout, field_file8, harmonic_file, read_tfu, run_stress, stress_deck
    from ...modules.stress import TSHELL_FACE_FILE, TSHELL_MAX_FILE
    from ...prep.commands.thickshell import read_tshell_max, write_thshlstr_file
    r = VPResult()
    wd = Path(workdir)
    _, _, E, nu, w, _, _ = MAT54
    G_SI = 9.81
    Dm, Db, Ds = tshell.plate_rigidities(E, nu, T54)
    a = harmonic_file(wd / "eq.acc", NFFT54, DT54, K0_54, A0_54)
    ug = np.fft.irfft(S.displacement_spectrum(np.fft.rfft(a), S.fourier_grid(NFFT54, DT54), G_SI), NFFT54)
    ugmax = float(np.max(np.abs(ug)))
    eps = np.array([2.0e-4, -0.5e-4, 1.0e-4])         # membrane strains (exx, eyy, gxy)
    kap = np.array([3.0e-4, 1.0e-4, -0.5e-4])         # curvatures (beta_x,x, beta_y,y, beta_x,y + beta_y,x)
    for eint in (0, 1):
        d, (o, e1, e2) = tshell_plate_deck(eint)
        rc, out = run_house(wd, d)
        r.require(f"EINT {eint}: HOUSE run with TSHELL elements succeeds", rc == 0, out[-600:] if rc else "")
        if rc:
            continue
        n3 = np.cross(e1, e2)

        def loc(x):
            q = np.asarray(x) - o
            return q @ e1, q @ e2

        def membrane(node, x, f):
            X, Y = loc(x)
            u = eps[0] * X + 0.5 * eps[2] * Y
            v = eps[1] * Y + 0.5 * eps[2] * X
            return np.concatenate([u * e1 + v * e2, np.zeros(3)]).astype(complex)

        def bending(node, x, f):                           # beta = (kxx X + kxy/2 Y, kyy Y + kxy/2 X), beta = -grad w
            X, Y = loc(x)
            bx = kap[0] * X + 0.5 * kap[2] * Y
            by = kap[1] * Y + 0.5 * kap[2] * X
            wv = -(0.5 * kap[0] * X * X + 0.5 * kap[1] * Y * Y + 0.5 * kap[2] * X * Y)
            theta = by * (-e1) + bx * e2                   # theta_x' = -beta_y, theta_y' = beta_x
            return np.concatenate([wv * n3, theta]).astype(complex)
        write_thshlstr_file(wd, 1)                          # what THSHLSTR,1 writes for STRESS
        for name, field, n0, m0 in (("pure membrane", membrane, Dm @ eps, np.zeros(3)),
                                    ("pure bending", bending, np.zeros(3), Db @ kap)):
            field_file8(wd, FNUM54, DF54, field, NFFT54, DT54)
            dk = stress_deck(NFFT54, DT54, itran=1)
            add_eout(dk, 1, [1, 2, 3, 4], [1] * 8)
            rc, out = run_stress(wd, dk)
            r.require(f"EINT {eint}, {name}: STRESS run succeeds", rc == 0, out[-800:] if rc else "")
            if rc:
                continue
            r.require(f"EINT {eint}, {name}: listing reports THSHLSTR 1 from THSHLSTR.opt",
                      "(THSHLSTR): 1 (THSHLSTR.opt" in out)
            ref8 = np.concatenate([n0, np.zeros(2), m0])
            scale = max(np.abs(n0).max() / T54, 6.0 * np.abs(m0).max() / T54 ** 2)
            # STF of the basic components (frequency independent) = closed form
            err = 0.0
            for c, comp in enumerate(tshell.COMPONENTS):
                _, h = read_tfu(wd, "TSHELL", 1, 1, comp)
                err = worse(err, float(np.max(np.abs(h - ref8[c]))))
            r.check(f"EINT {eint}, {name}: STF NXX..MXY = D eps / D kappa (max error / max |N|, |M|)", err
                    / max(np.abs(ref8).max(), 1e-300), 0.0, atol=1e-10)
            comps, mx = read_tshell_max(wd / TSHELL_MAX_FILE)
            vals = mx[(1, 1)]
            r.check(f"EINT {eint}, {name}: maxima of the 8 components = |closed form| x max|u_g|",
                    float(np.max(np.abs(vals - np.abs(ref8) * ugmax))) / (np.abs(ref8).max() * ugmax), 0.0, atol=1e-9)
            faces = read_face_file(wd / TSHELL_FACE_FILE)
            Nh, Mh = vals[[0, 1, 2]], vals[[5, 6, 7]]
            top, bot = faces[(1, 1, "++")], faces[(1, 1, "+-")]
            if name == "pure membrane":
                r.check(f"EINT {eint}, pure membrane: sigma_xx top = bottom = Nxx/t (++ vs +-)",
                        abs(top[0] - bot[0]) / scale, 0.0, atol=1e-12)
                r.check(f"EINT {eint}, pure membrane: SXX = Nxx/t = |E/(1-nu^2)(exx + nu eyy)| max|u_g|", top[0],
                        abs(n0[0]) / T54 * ugmax, rtol=1e-9)
                r.check(f"EINT {eint}, pure membrane: SYY = Nyy/t", top[1], abs(n0[1]) / T54 * ugmax, rtol=1e-9)
            else:
                # N = 0 up to round-off (the seismic rigid-body subtraction of STRESS, |S| eps): 1e-9
                r.check(f"EINT {eint}, pure bending: SXX(++) = +6 Mxx/t^2", top[0], 6.0 * Mh[0] / T54 ** 2, rtol=1e-9)
                r.check(f"EINT {eint}, pure bending: SXX(+-) = -6 Mxx/t^2", bot[0], -6.0 * Mh[0] / T54 ** 2, rtol=1e-9)
                r.check(f"EINT {eint}, pure bending: 6 Mxx/t^2 = 6 |D (kxx + nu kyy)| max|u_g| / t^2", top[0],
                        6.0 * abs(m0[0]) * ugmax / T54 ** 2, rtol=1e-9)
                r.check(f"EINT {eint}, pure bending: SYY(--) = -6 Myy/t^2", faces[(1, 1, "--")][1],
                        -6.0 * Mh[1] / T54 ** 2, rtol=1e-9)
            # every face value is s_N N/t + s_M 6 M/t^2 of the listed maxima (exact arithmetic)
            ferr = 0.0
            for perm, sn, sm in tshell.FACE_PERMUTATIONS:
                v = faces[(1, 1, perm)]
                for j in (0, 1):
                    ferr = worse(ferr, abs(v[j] - (sn * Nh[j] / T54 + sm * 6.0 * Mh[j] / T54 ** 2)) / scale)
            r.check(f"EINT {eint}, {name}: SXX, SYY of the 4 permutations = s_N N/t + s_M 6 M/t^2 (D-TSH-01)", ferr, 0.0,
                    atol=1e-12)
            # shear (++ only), principal values and strains of every permutation
            txy = Nh[2] / T54 + 6.0 * Mh[2] / T54 ** 2
            Gm = E / (2 * (1 + nu))
            perr = 0.0
            for perm in ("++", "--", "+-", "-+"):
                v = faces[(1, 1, perm)]
                sxx, syy, t_ = v[0], v[1], v[2]
                c_, rr = 0.5 * (sxx + syy), math.hypot(0.5 * (sxx - syy), t_)
                exx, eyy = (sxx - nu * syy) / E, (syy - nu * sxx) / E
                perr = worse(perr, abs(t_ - txy) / scale, abs(v[3] - (c_ + rr)) / scale, abs(v[4] - (c_ - rr)) / scale,
                           abs(v[5] - exx) * E / scale, abs(v[6] - eyy) * E / scale, abs(v[7] - t_ / Gm) * E / scale)
            r.check(f"EINT {eint}, {name}: TXY (++), principal stresses and plane-stress strains of all 4 permutations",
                    perr, 0.0, atol=1e-12)
    r.notes.append(f"Plate 2 x 1.5 x {T54} (2 x 2 TSHELL elements) in an inclined plane, E = 3e7, nu = 0.2, undamped; "
                   f"harmonic control motion {A0_54} g at bin {K0_54} (f = {K0_54 * DF54:.4f} Hz), max|u_g| = {ugmax:.6g}.")
    return r


# ======================================================================================
# VP-O1: the HOUSE node optimizer
# ======================================================================================
def equation_profile(K) -> Tuple[int, int]:
    """(bandwidth, profile) of the sparsity pattern of a square matrix (lower envelope)."""
    A = sp.coo_matrix(K)
    if A.nnz == 0:
        return 0, 0
    n = A.shape[0]
    first = np.arange(n)
    np.minimum.at(first, A.row, A.col)
    return int(np.max(np.abs(A.row - A.col))), int(np.sum(np.arange(n) - first))


def optimizer_model(seed: int = 5):
    """Embedded TSHELL basement (walls, base and roof slabs) with a two-storey TSHELL superstructure
    (perimeter walls and floor slabs) on a 2-layer site, FI-FSIN interaction set (lateral and bottom
    faces of the excavation), with the node numbers scrambled at random (the interaction table keeps
    its bottom-up order).  Units kN, m, t."""
    from .. import builders as B
    from ...modules.house import renumbered_deck
    site = B.layered_site([(1.5, 150.0, 300.0, 1.9, 0.05), (1.5, 200.0, 400.0, 1.9, 0.05)], (350.0, 700.0, 2.0, 0.03),
                          nl=8)
    a, ndiv = 1.5, 3
    mdl = B.embedded_box(site, half_width=a, n_emb=2, ndiv=ndiv, method="FSIN", structure="shell", thick=0.25)
    hb = mdl.house
    for g in hb.groups:
        if g[1] == 3:
            g[1] = 5                                        # SHELL basement -> TSHELL
    mat = hb.elastic(3.0e7, 0.2, 2.4, 0.05)
    gsup = hb.group(5, "superstructure")
    gx = np.linspace(-a, a, ndiv + 1)
    levels = (0.0, 2.0, 4.0)
    for zb, zt in zip(levels[:-1], levels[1:]):             # perimeter walls, storey by storey
        for k in range(ndiv):
            for c in (gx[0], gx[-1]):
                hb.shell([hb.node(gx[k], c, zb), hb.node(gx[k + 1], c, zb), hb.node(gx[k + 1], c, zt),
                          hb.node(gx[k], c, zt)], mat, 0.25, group=gsup)
                hb.shell([hb.node(c, gx[k], zb), hb.node(c, gx[k + 1], zb), hb.node(c, gx[k + 1], zt),
                          hb.node(c, gx[k], zt)], mat, 0.25, group=gsup)
    for zt in levels[1:]:                                   # floor slabs
        for j in range(ndiv):
            for i in range(ndiv):
                hb.shell([hb.node(gx[i], gx[j], zt), hb.node(gx[i + 1], gx[j], zt), hb.node(gx[i + 1], gx[j + 1], zt),
                          hb.node(gx[i], gx[j + 1], zt)], mat, 0.3, group=gsup)
    d = hb.deck("m")
    ids = [int(row[0]) for row in d.table("nodes").rows]
    rng = np.random.default_rng(seed)
    scrambled = dict(zip(ids, (int(v) for v in rng.permutation(ids))))
    return site, mdl, renumbered_deck(d, scrambled)


@problem("VP-O1", "HOUSE node-numbering optimizer: identical TFs through the .map, smaller profile", tier="P1",
         modules=["HOUSE", "ANALYS"], source="requirements 4.4 item 5, D-HOU-03, D-FIL-09; spec 05b test 8")
def vpo1(workdir):
    from .. import builders as B
    r = VPResult()
    wd = Path(workdir)
    site, mdl, d0 = optimizer_model()
    fs = B.FrequencySet.harmonic(1.0, [3, 9])
    B.run_soil(wd, "m", site, fs, layer=mdl.layer, rad=mdl.rad)
    # ---- reference run: original (scrambled) numbering
    d0["optimize"] = 0
    B.write_deck(wd, "m", d0)
    B.run("HOUSE", wd, "m")
    B.run_analys(wd, "m", fs)
    f4a = read_container(wd / "m.N4", "FILE4")
    f8a = B.read_file8(wd)
    ka = read_container(wd / "COOSK", "COOSK")
    r.require("reference run writes no .map", not (wd / "m.map").exists())
    # ---- optimised run
    d1 = decks.read(wd / "m.hou", "HOUSE")
    d1["optimize"] = 1
    B.write_deck(wd, "m", d1)
    B.run("HOUSE", wd, "m")
    out = B.listing(wd, "m", "HOUSE")
    B.run_analys(wd, "m", fs)
    f4b = read_container(wd / "m.N4", "FILE4")
    f8b = B.read_file8(wd)
    kb = read_container(wd / "COOSK", "COOSK")
    mp = RN.read_map(wd / "m.map")
    old_ids = [int(v) for v in f4a["node_id"]]
    r.require(".map holds every node once (old -> new is a bijection onto 1..N)",
              sorted(mp) == sorted(old_ids) and sorted(mp.values()) == list(range(1, len(old_ids) + 1)))
    r.require("listing: node-numbering optimizer section and the new-number warning",
              "Node-numbering optimizer" in out and "use the NEW numbers" in out)
    # ---- identical transfer functions through the map
    eqb = {(int(n), int(k)): i for i, (n, k) in enumerate(zip(f8b["eq_node"], f8b["eq_dof"]))}
    cols = [eqb.get((mp[int(n)], int(k)), -1) for n, k in zip(f8a["eq_node"], f8a["eq_dof"])]
    r.require("every FILE8 equation of the reference run exists in the optimised run", min(cols) >= 0)
    Ha, Hb = np.asarray(f8a["H"]), np.asarray(f8b["H"])
    if min(cols) >= 0:
        r.check("FILE8 H mapped back through .map: max |H_opt - H_ref| / max |H_ref|",
                float(np.abs(Hb[:, cols] - Ha).max() / np.abs(Ha).max()), 0.0, atol=1e-10)
    r.check("number of equations unchanged", len(f8b["eq_node"]), len(f8a["eq_node"]), atol=0)
    # ---- bandwidth and profile of the structure + excavated-soil matrices
    pat = lambda c: (abs(c.sparse("Ks")) + abs(c.sparse("Ke"))).tocsr()
    bwa, pfa = equation_profile(pat(ka))
    bwb, pfb = equation_profile(pat(kb))
    r.require(f"equation bandwidth reduced ({bwa} -> {bwb})", bwb < bwa)
    r.require(f"equation profile reduced ({pfa} -> {pfb})", pfb < pfa)
    r.inform("equation profile after / before", pfb / pfa, 1.0)
    # ---- interaction nodes: same nodes, same order, ascending new numbers
    ia, ib = [int(v) for v in f4a["int_node"]], [int(v) for v in f4b["int_node"]]
    r.require("interaction nodes keep their relative (bottom-up) order", ib == [mp[n] for n in ia])
    r.require("interaction nodes numbered ascending in interaction order", all(np.diff(ib) > 0))
    r.require("FILE4 x_node_old_id maps the new numbers back", [mp[int(o)] for o in f4b["x_node_old_id"]]
              == [int(v) for v in f4b["node_id"]] and int(f4b.meta.get("x_optimized", 0)) == 1)
    # ---- the .hounew deck is the model HOUSE analysed: re-run it without the optimizer
    dn = decks.read(wd / "m.hounew", "HOUSE")
    r.require(".hounew has optimize = 0", int(dn["optimize"]) == 0)
    B.write_deck(wd, "m", dn)
    B.run("HOUSE", wd, "m")
    f4c = read_container(wd / "m.N4", "FILE4")
    kc = read_container(wd / "COOSK", "COOSK")
    same = all(np.array_equal(np.asarray(f4c[k]), np.asarray(f4b[k])) for k in
               ("node_id", "eq_node", "eq_dof", "int_node", "int_eq", "elem_nodes"))
    r.require(".hounew re-run (optimize 0) reproduces the optimised FILE4 numbering", same)
    r.check(".hounew re-run: max |Ks - Ks_opt| / max |Ks_opt|", float(abs(kc.sparse("Ks") - kb.sparse("Ks")).max()
                                                                      / abs(kb.sparse("Ks")).max()), 0.0, atol=1e-14)
    r.notes.append(f"{len(old_ids)} nodes, {len(ia)} interaction nodes (FI-FSIN), {len(f8a['eq_node'])} equations; "
                   f"equation bandwidth {bwa} -> {bwb}, profile {pfa} -> {pfb}; frequencies "
                   f"{', '.join(f'{float(v):g}' for v in fs.freq)} Hz.")
    return r
