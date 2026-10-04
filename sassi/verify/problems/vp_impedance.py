"""Foundation impedance and inertial SSI: the end-to-end physics checks of SITE + POINT + HOUSE (+ FORCE)
+ ANALYS (requirements §6.3 VP-10, VP-11, VP-13, VP-14, VP-17; R2 §B, §C, §F; R1 §3.4; D-ANL-07, D-PNT-05).

* VP-10  rigid disk on a half-space: static stiffness vs 4GR/(1-nu), 8GR/(2-nu), 8GR^3/(3(1-nu)), 16GR^3/3
* VP-11  rigid disk dynamic impedance vs Veletsos-Verbic (a0 <= 1.5) and the high-frequency dashpots
* VP-13  Tajirian & Tabatabaie (1985): disk on a layer over a rigid base -- resonances and statics
* VP-14  square rigid footing statics vs Pais & Kausel (1988)
* VP-17  inertial SSI: SDOF on a rigid massless disk = the exact 3-DOF model (R2 F.2) fed with the
         code's own impedances, and the Veletsos-Meek / Veletsos-Verbic worked example (R2 F.1)

How a SASSI model gives a foundation impedance (for the structural engineer)
---------------------------------------------------------------------------
The foundation is a set of *interaction nodes* on the ground surface: a ring mesh for a disk (ring k
carries 6k nodes, element size h = R/n) or a square grid.  POINT computes the soil flexibility F_ff
between the nodes from point-load solutions (R1 §3.3); the near-singular diagonal terms come from a
small axisymmetric core of radius R0 = 0.85 h (ring/triangular meshes) or 0.90 h (square meshes), the
manual's central-zone rules (R1 §3.4; for a square cell, 0.90 h is the RADIUS formula 0.9 sqrt(A_plan)
of D-PNT-05).  ANALYS inverts F_ff into the soil impedance X_ff and condenses
it onto the six rigid-body motions of the mat, ``K_G = T^T X_ff T`` (``<impe>`` = 2, D-ANL-07).  Every
translation of every node follows the rigid body, so K_G is the impedance of a rigid foundation in
**welded** contact.  The classical closed forms (Boussinesq punch, Veletsos & Verbic) assume
**relaxed** contact (no shear traction under a vertically moving / rocking disk).  For a like-for-like
comparison the relaxed impedance is also formed here from the same F_ff: vertical and rocking from
the vertical block only, ``K = T_z^T F_zz^-1 T_z`` (horizontal tractions zero), horizontal and torsion
from the horizontal block only.

Model rules used by every problem (documented with the mesh/R0 study in
docs/verification/impedance_study.md):

* mat: interaction nodes on the surface made rigid by stiff massless BEAMS to the centre node
  (:func:`sassi.verify.builders.surface_rigid_mat`); K_G does not depend on the mat (X_ff only);
* soil column: top sublayer thickness = the node spacing h (the setting the 0.85h / 0.90h rules are
  calibrated for: cube-like near-surface elements), growing by 12 % per sublayer up to
  min(lambda_min/8, 4R) (h <= lambda/8 layer rule at the highest frequency), down to 20-30 R, then
  20 generated half-space sublayers (UNIFORM law, D-W1-01) and base dashpots;
* near-elastic soil (beta = 1e-4) where elastic closed forms are the reference;
* "after mesh extrapolation": three meshes with the same R0/h and t0/h ratios (self-similar
  refinement); Richardson extrapolation with the formal order p = 1 (tributary-area and
  central-zone errors are O(h)) from the two finest meshes; the observed order of the three meshes
  is reported.

Independent references computed here (derived, type "D")
--------------------------------------------------------
The closed forms of R2 are approximations for some quantities (welded vs relaxed contact, curve
fits).  Where a VP fails, the evidence is an independent rigorous solution:

* :func:`static_bem` -- rigid foundation on an elastic half-space with the exact Boussinesq /
  Cerruti surface kernels (piecewise-constant tractions, collocation; welded or relaxed contact).
  It reproduces the exact disk values (relaxed 6.000 / 4.000, Mossakovskii's welded vertical
  4GR ln(3-4nu)/(1-2nu) = 6.130, Reissner-Sagoci torsion 5.333 GR^3) to 0.1 %.
* :func:`lamb_kernels` / :func:`relaxed_bem_impedance` / :func:`welded_bem_impedance` -- the same BEM
  with the exact dynamic surface Green's functions of the half-space (Lamb's problem, wavenumber
  integrals on a contour in the upper half plane; the welded version adds the vertical/horizontal
  coupling kernel); the kernels reproduce Wong's (1975) table (R2 D.3) to 0.001.  :data:`VP11_BEM`
  stores the mesh-extrapolated relaxed disk results (regenerate in ~1 min, see the unit tests);
  :func:`welded_disk_impedance` is the rigorous counterpart of K_G (VP-17 guard).
* :func:`layer_zgv_onset` -- continuum P-SV dispersion of a layer on a rigid base (propagator
  expm(A H), no thin-layer discretisation): the zero-group-velocity onset behind the VP-13 vertical peak.
* :func:`lowfreq_damping` -- the exact low-frequency radiation damping of the disk (first-order
  perturbation with the static punch tractions and the imaginary part of the Lamb kernel).

Guards of the known failures: every check that fails against an approximate R2 reference has a
passing check of the same quantity against a derived reference (VP-11 relaxed c_r vs the BEM, VP-13
ZGV onset and the lightly damped second peak, VP-14 K_xx vs the welded BEM, VP-17 peak vs F.2 with the
welded BEM impedance), so that an xfail cannot hide a regression.  Every K_G is checked against an
independent LU solve T^T F_ff^-1 T and for complex symmetry (:func:`kg_checks`).

Disk space: FILE11 normally stores X_ff of every frequency (up to 256 MB, D-ANL-07); the problems run
ANALYS with FILE11 holding K_G and T only (:func:`file11_without_matrices`, local to the calling
thread), and delete the binary intermediate files of every run once the results are in memory
(:func:`purge_work_files`; decks and listings are kept).

Units: m, t, kN, s (VP-13: the dimensionless units of Tajirian & Tabatabaie, r = 0.5, H = 3, Vs = rho = 1).
"""
from __future__ import annotations

import contextlib
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy import integrate, linalg, optimize, special
from scipy.interpolate import CubicSpline

from ...core import ssi_solver as ss
from ...core.flexibility import flexibility_matrix
from ...core.tlm import quiet_fpe
from ...io.container import read_container
from .. import VPResult, problem
from .. import builders as B

# ======================================================================================
# Soil and references
# ======================================================================================
VS, RHO, NU = 100.0, 2.0, 1.0 / 3.0            # m/s, t/m3; G = 20 000 kPa
VPP = VS * math.sqrt(2 * (1 - NU) / (1 - 2 * NU))   # 200 m/s
GMOD = RHO * VS ** 2
GROWTH = 1.12                                   # sublayer growth below the surface
BETA_ELASTIC = 1.0e-4                           # near-elastic soil

#: R2 B.1 (nu = 1/3), units G R and G R^3
VP10_REF = {"Kv": 6.0, "Kh": 4.8, "Kr": 4.0, "Kt": 16.0 / 3.0}
#: R2 C.2 Veletsos-Verbic, nu = 1/3: a0 -> (k_r, c_r, k_v, c_v); horizontal k = 1, c = alpha1 = 0.65
VP11_VV = {0.5: (0.9529, 0.0235, 0.9517, 0.7886), 1.0: (0.8400, 0.0800, 0.8634, 0.8593),
           1.5: (0.7120, 0.1440, 0.7934, 0.9152)}
VP11_ALPHA1 = 0.65
#: R2 C.1 normalised high-frequency dashpots c(inf) (= plane-wave dashpots rho V A, rho Vp I, ...)
VP11_CINF = {"h": math.pi * (2 - NU) / 8, "v": math.pi * (1 - NU) / 4 * 2.0, "r": 3 * math.pi / 32 * (1 - NU) * 2.0,
             "t": 3 * math.pi / 32}
#: Independent rigorous reference (derived here): relaxed-contact disk, exact Lamb kernels, BEM on
#: polar meshes n = 16 and 24 extrapolated (p = 1) -- see :func:`relaxed_bem_impedance`.
#: a0 -> {coefficient: value}; k = Re K / K0, c = Im K / (a0 K0) with the BEM static K0
#: (static: Kv 6.0011, Kh 4.8009, Kr 4.0023, Kt 5.3365 vs exact 6, 4.8 (approx.), 4, 5.3333).
VP11_BEM = {0.25: dict(kh=0.9974, ch=0.5815, kr=0.9826, cr=0.0143, kv=0.9921, cv=0.7854, kt=0.9879, ct=0.0085),
            0.5: dict(kh=0.9899, ch=0.5853, kr=0.9381, cr=0.0501, kv=0.9689, cv=0.7931, kt=0.9557, ct=0.0310),
            1.0: dict(kh=0.9652, ch=0.5991, kr=0.8253, cr=0.1353, kv=0.8826, cv=0.8239, kt=0.8667, ct=0.0913),
            1.5: dict(kh=0.9410, ch=0.6170, kr=0.7243, cr=0.2005, kv=0.7645, cv=0.8722, kt=0.7825, ct=0.1440)}
#: R2 B.2 Pais & Kausel square (L = B, nu = 1/3), units G B and G B^3
VP14_PK = {"Kz": 7.05, "Kx": 5.52, "Ky": 5.52, "Kzz": 8.31, "Kxx": 6.0, "Kyy": 6.0}
VP14_GAZETAS = {"Kz": 6.81, "Kx": 5.40, "Ky": 5.40, "Kzz": 8.3471, "Kxx": 5.3975, "Kyy": 5.5836}


# ======================================================================================
# Helpers: soil column, frequencies, chain
# ======================================================================================
def graded_thicknesses(t0: float, tmax: float, depth: float, growth: float = GROWTH) -> List[float]:
    """Sublayer thicknesses from the surface down to ``depth``: t0, t0 g, t0 g^2 ..., every one of them
    (the first included) capped at tmax -- tmax is the h <= lambda/8 layer rule of :func:`soil_column`,
    which must hold even when the requested top sublayer t0 is thicker."""
    if t0 <= 0 or tmax <= 0 or depth <= 0:
        raise ValueError("t0, tmax and depth must be positive")
    out, z, t = [], 0.0, min(float(t0), float(tmax))
    while z < depth * (1 - 1e-12):
        out.append(t)
        z += t
        t = min(t * growth, tmax)
    return out


def soil_column(t0: float, tmax: float, depth: float, vs: float = VS, rho: float = RHO, nu: float = NU,
                beta: float = BETA_ELASTIC, nl: int = 20, cmodform: int = 0, gravity: float = B.GRAVITY) -> B.Site:
    """Uniform soil as graded TOPL sublayers on the same soil as half-space (``nl`` generated
    sublayers, D-SIT-02 UNIFORM law, base dashpots); ``nl = 0``: rigid base."""
    lay = [B.SoilLayer.from_nu(t, vs, rho, nu, beta) for t in graded_thicknesses(t0, tmax, depth)]
    hs = B.SoilLayer.from_nu(0.0, vs, rho, nu, beta)
    return B.Site(lay, hs, gravity=gravity, nl=nl, cmodform=cmodform)


def a0_frequencies(a0s: Sequence[float], R: float, da: float = 0.01, vs: float = VS) -> Tuple[B.FrequencySet, np.ndarray]:
    """Harmonic frequency set with ``a0 = w R / Vs`` on a grid of step ``da``; returns the set and the
    exact a0 of each frequency number."""
    df = da * vs / (2 * math.pi * R)
    fnum = sorted({max(1, int(round(a / da))) for a in a0s})
    return B.FrequencySet.harmonic(df, fnum), np.asarray(fnum, float) * da


def rigid_transform(xy: np.ndarray) -> np.ndarray:
    xyz = np.c_[np.asarray(xy, float).reshape(-1, 2), np.zeros(len(xy))]
    return ss.rigid_body_transform(xyz, (0.0, 0.0, 0.0))


@quiet_fpe
def relaxed_impedance(F: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Relaxed-contact rigid impedance (6x6) from the flexibility F_ff (3n x 3n, node-major): vertical
    and rocking from the vertical block (horizontal tractions zero), horizontal and torsion from the
    horizontal block (vertical tractions zero); couplings between the two groups are zero."""
    n = F.shape[0] // 3
    K = np.zeros((6, 6), complex)
    iz = 3 * np.arange(n) + 2
    ih = np.sort(np.r_[3 * np.arange(n), 3 * np.arange(n) + 1])
    for idx, dofs in ((iz, [2, 3, 4]), (ih, [0, 1, 5])):
        Tr = T[np.ix_(idx, dofs)]
        K[np.ix_(dofs, dofs)] = Tr.T @ np.linalg.solve(F[np.ix_(idx, idx)], Tr)
    return K


@contextlib.contextmanager
def file11_without_matrices():
    """Run ANALYS with FILE11 holding K_G and T only.

    FILE11 normally also stores X_ff of every frequency while that takes at most
    ``analys.FILE11_MAX_X_BYTES`` (256 MB, D-ANL-07): up to 100 MB per frequency for the finest meshes
    here.  The impedance problems only need K_G, so the limit is set to 0 with
    :func:`sassi.modules.analys.file11_x_cap` -- a context variable, so only the ANALYS runs of this
    thread are affected (FILE11 meta ``x_stored`` = False); the module default is never touched."""
    from ...modules import analys as _analys
    with _analys.file11_x_cap(0):
        yield


def analys_global_impedance(wd: Path, model: str, fs: B.FrequencySet, **params) -> np.ndarray:
    """K_G (nF, 6, 6) of ANALYS ``<impe>`` = 2 about the origin (FILE11 'KG', D-ANL-07); ``params`` are
    further ANALYS deck parameters.  FILE8 of the run is written as usual."""
    with file11_without_matrices():
        B.run_analys(wd, model, fs, impe=2, **params)
    f11 = read_container(Path(wd) / "FILE11", "FILE11")
    if bool(f11.meta.get("x_stored")):                       # pragma: no cover - guarded by the context
        raise RuntimeError("FILE11 unexpectedly stores X_ff")
    return np.asarray(f11["KG"], complex)


#: keep every file of the VP runs (debugging); by default only the decks and listings survive a run
KEEP_WORK_FILES = False


def purge_work_files(wd: Path) -> int:
    """Delete the binary intermediate files of a finished run in ``wd`` (FILE1/2/3, FILE8, FILE11,
    HOUSE matrices, FOUN* tables ...) once the results are in memory; the input decks and the module
    listings are kept.  A VP writes 7-75 MB of such files, which the verification suite cannot afford
    on a nearly full disk.  Returns the number of bytes freed (0 when :data:`KEEP_WORK_FILES`)."""
    from ...modules.base import DECK_EXT
    if KEEP_WORK_FILES or not Path(wd).is_dir():
        return 0
    keep = {e for e in DECK_EXT.values() if e} | {".out"}
    freed = 0
    for p in Path(wd).iterdir():
        if p.is_file() and p.suffix not in keep:
            freed += p.stat().st_size
            p.unlink()
    return freed


def independent_global_impedance(F: np.ndarray, T: np.ndarray) -> np.ndarray:
    """K = T^T F^-1 T by an LU solve of F (``numpy.linalg.solve``): the reference for ANALYS's K_G,
    independent of the impedance inversion and condensation of ``ssi_solver`` (``impedance_matrix``,
    ``global_impedance``) that ANALYS uses.  F is complex symmetric, not Hermitian (LU, no Cholesky)."""
    with np.errstate(all="ignore"):
        return T.T @ np.linalg.solve(F, T)


def kg_checks(r: VPResult, label: str, KG: np.ndarray, KD: np.ndarray) -> None:
    """Two consistency checks of the ANALYS global impedance (D-ANL-07): K_G vs the independent
    T^T F_ff^-1 T, and complex symmetry K_G = K_G^T (reciprocity of the soil flexibility)."""
    KG, KD = np.asarray(KG), np.asarray(KD)
    err = float(np.abs(KG - KD).max() / np.abs(KD).max())
    asym = float(np.abs(KG - np.swapaxes(KG, -1, -2)).max() / np.abs(KG).max())
    r.check(f"{label}: ANALYS K_G = T^T X_ff T vs T^T F_ff^-1 T by an independent LU solve (max rel. difference)",
            err, 0.0, atol=1e-8)
    r.check(f"{label}: complex symmetry of K_G (max |K_G - K_G^T| / max |K_G|)", asym, 0.0, atol=1e-8)


def impedance_chain(wd: Path, site: B.Site, fs: B.FrequencySet, mdl: B.Model, model: str = "m",
                    relaxed: bool = True, cleanup: bool = True) -> Dict[str, object]:
    """SITE (vertical SV) -> POINT -> HOUSE -> ANALYS ``<impe>`` = 2 on a surface mat.

    Returns ``KG`` (ANALYS, welded), ``KD`` = T^T F_ff^-1 T from the binding flexibility API and an
    independent LU solve (:func:`independent_global_impedance`), ``KR`` (relaxed contact) and the
    FILE2 / FILE3 containers (in memory).  With ``cleanup`` the intermediate files are deleted
    afterwards (:func:`purge_work_files`)."""
    wd = Path(wd)
    wd.mkdir(parents=True, exist_ok=True)
    files = B.run_soil(wd, model, site, fs, layer=0, rad=mdl.rad)
    B.run_house(wd, model, mdl)
    nodes = list(mdl["mat"])
    KG = analys_global_impedance(wd, model, fs)
    if cleanup:
        purge_work_files(wd)
    xy = np.array([mdl.house.xyz(n)[:2] for n in nodes])
    T = rigid_transform(xy)
    KD, KR = [], []
    f3 = files["FILE3"]
    for q in range(len(fs.fnum)):
        F = flexibility_matrix(f3, q, xy, np.ones(len(nodes), int))
        KR.append(relaxed_impedance(F, T) if relaxed else np.zeros((6, 6), complex))
        KD.append(independent_global_impedance(F, T))
    return {"KG": KG, "KD": np.array(KD), "KR": np.array(KR), "FILE3": f3, "FILE2": files["FILE2"], "xy": xy}


def richardson(h: Sequence[float], v: Sequence, p: float = 1.0):
    """Richardson extrapolation to h -> 0 from the two finest values (``h`` descending)."""
    h1, h2 = float(h[-2]), float(h[-1])
    v1, v2 = np.asarray(v[-2]), np.asarray(v[-1])
    return v2 + (v2 - v1) * h2 ** p / (h1 ** p - h2 ** p)


def observed_order(h: Sequence[float], v: Sequence[float]) -> float:
    """Order p of ``v(h) = v0 + a h^p`` through three values (NaN when not monotone)."""
    h1, h2, h3 = (float(x) for x in h[-3:])
    v1, v2, v3 = (float(x) for x in v[-3:])
    d1, d2 = v1 - v2, v2 - v3
    if d1 == 0 or d2 == 0 or d1 * d2 < 0:
        return float("nan")
    ratio = d1 / d2

    def g(p):
        return (h1 ** p - h2 ** p) / (h2 ** p - h3 ** p) - ratio
    lo, hi = 0.05, 8.0
    if g(lo) * g(hi) > 0:
        return float("nan")
    return float(optimize.brentq(g, lo, hi))


def peak_near(x: np.ndarray, y: np.ndarray, centre: float, window: float = 0.10) -> float:
    """Location of the largest local maximum of y(x) within ``centre (1 +- window)``, refined by a
    parabola through the three grid points (NaN when there is none)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    idx = [i for i in range(1, len(x) - 1) if y[i] >= y[i - 1] and y[i] > y[i + 1]
           and abs(x[i] - centre) <= window * centre]
    if not idx:
        return float("nan")
    i = max(idx, key=lambda k: y[k])
    a, b, c = y[i - 1], y[i], y[i + 1]
    den = a - 2 * b + c
    d = 0.5 * (a - c) / den if den != 0 else 0.0
    return float(x[i] + d * (x[i + 1] - x[i]))


# ======================================================================================
# Independent reference 1: static rigid-foundation BEM (exact Boussinesq / Cerruti kernels)
# ======================================================================================
def surface_kernel(dx, dy, nu: float, G: float = 1.0) -> np.ndarray:
    """Static surface Green's function of a half-space, ``K[..., a, b]`` = u_a at (dx, dy) = obs - load
    due to a unit surface force P_b (z up): Cerruti for P_x, P_y, Boussinesq for P_z (R2 D.1)."""
    r = np.hypot(dx, dy)
    Bf = 1.0 / (2 * np.pi * G * r)
    C = (1 - 2 * nu) / (4 * np.pi * G * r * r)
    K = np.empty(np.shape(dx) + (3, 3))
    K[..., 0, 0] = Bf * ((1 - nu) + nu * dx * dx / r ** 2)
    K[..., 0, 1] = K[..., 1, 0] = Bf * nu * dx * dy / r ** 2
    K[..., 1, 1] = Bf * ((1 - nu) + nu * dy * dy / r ** 2)
    K[..., 0, 2] = C * dx
    K[..., 1, 2] = C * dy
    K[..., 2, 0] = -C * dx
    K[..., 2, 1] = -C * dy
    K[..., 2, 2] = Bf * (1 - nu)
    return K


def _cell_quadrature(V: np.ndarray, n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Gauss points (nc, n*n, 2) and weights (nc, n*n) of bilinear quadrilaterals V (nc, 4, 2)."""
    g, w = np.polynomial.legendre.leggauss(n)
    s, t = (a.ravel() for a in np.meshgrid(g, g, indexing="ij"))
    ww = np.outer(w, w).ravel()
    N = 0.25 * np.array([(1 - s) * (1 - t), (1 + s) * (1 - t), (1 + s) * (1 + t), (1 - s) * (1 + t)])
    dNs = 0.25 * np.array([-(1 - t), (1 - t), (1 + t), -(1 + t)])
    dNt = 0.25 * np.array([-(1 - s), -(1 + s), (1 + s), (1 - s)])
    P = np.einsum("aq,cad->cqd", N, V)
    Js = np.einsum("aq,cad->cqd", dNs, V)
    Jt = np.einsum("aq,cad->cqd", dNt, V)
    return P, ww[None, :] * (Js[..., 0] * Jt[..., 1] - Js[..., 1] * Jt[..., 0])


def _self_integral(V: np.ndarray, P: np.ndarray, nu: float, n: int = 24) -> np.ndarray:
    """Integral of the 1/r kernel over the cell V (4, 2) seen from P inside it: polar coordinates around
    P, ``Int g(phi)/rho dA = Int g(phi) rho_edge(phi) dphi`` over the four edge triangles."""
    g, w = np.polynomial.legendre.leggauss(n)
    out = np.zeros((3, 3))
    for a in range(4):
        pa, pb = V[a] - P, V[(a + 1) % 4] - P
        fa, fb = math.atan2(pa[1], pa[0]), math.atan2(pb[1], pb[0])
        dphi = (fb - fa + math.pi) % (2 * math.pi) - math.pi
        phi = fa + 0.5 * dphi * (g + 1)
        e = pb - pa
        nrm = np.array([e[1], -e[0]]) / math.hypot(e[0], e[1])
        d = abs(float(pa @ nrm))
        dirs = np.c_[np.cos(phi), np.sin(phi)]
        rho = d / np.abs(dirs @ nrm)
        out += np.einsum("q,qab->ab", 0.5 * dphi * w * rho, surface_kernel(-dirs[:, 0], -dirs[:, 1], nu))
    return out


def bem_influence(V: np.ndarray, nu: float, nfar: int = 3, nnear: int = 8, near: float = 2.0):
    """Influence matrix A (3nc x 3nc): displacement at the cell centroids per unit uniform traction on
    each cell; returns (A, centroids, areas)."""
    nc = V.shape[0]
    Pw, Ww = _cell_quadrature(V, 2)
    area = Ww.sum(axis=1)
    Pc = np.einsum("cq,cqd->cd", Ww, Pw) / area[:, None]
    size = np.sqrt(area)
    A = np.zeros((nc, 3, nc, 3))
    for nq, which in ((nfar, False), (nnear, True)):
        Q, W = _cell_quadrature(V, nq)
        for i0 in range(0, nc, 64):
            i1 = min(nc, i0 + 64)
            dist = np.hypot(Pc[i0:i1, None, 0] - Pc[None, :, 0], Pc[i0:i1, None, 1] - Pc[None, :, 1])
            mask = (dist < near * np.maximum(size[i0:i1, None], size[None, :])) == which
            mask[np.arange(i1 - i0), np.arange(i0, i1)] = False
            ii, jj = np.nonzero(mask)
            if ii.size:
                d = Pc[i0 + ii][:, None, :] - Q[jj]
                A[i0 + ii, :, jj, :] = np.einsum("pq,pqab->pab", W[jj], surface_kernel(d[..., 0], d[..., 1], nu))
    for i in range(nc):
        A[i, :, i, :] = _self_integral(V[i], Pc[i], nu)
    return A.reshape(3 * nc, 3 * nc), Pc, area


def _bem_rigid(A: np.ndarray, Pc: np.ndarray, area: np.ndarray, relaxed: bool) -> np.ndarray:
    nc = Pc.shape[0]
    T = rigid_transform(Pc)
    W = np.repeat(area, 3)
    K = np.zeros((6, 6), A.dtype)
    if not relaxed:
        return (T * W[:, None]).T @ np.linalg.solve(A, T)
    iz = 3 * np.arange(nc) + 2
    ih = np.sort(np.r_[3 * np.arange(nc), 3 * np.arange(nc) + 1])
    for idx, dofs in ((iz, [2, 3, 4]), (ih, [0, 1, 5])):
        Tr = T[np.ix_(idx, dofs)]
        K[np.ix_(dofs, dofs)] = (Tr * W[idx, None]).T @ np.linalg.solve(A[np.ix_(idx, idx)], Tr)
    return K


def static_bem(V: np.ndarray, nu: float = NU, relaxed: bool = False) -> np.ndarray:
    """Static 6x6 stiffness (G = 1) of a rigid surface foundation on the cells V (nc, 4, 2): welded
    (all tractions) or relaxed contact.  Piecewise-constant tractions, collocation at the centroids;
    the error is O(h) (extrapolate)."""
    with np.errstate(all="ignore"):
        A, Pc, area = bem_influence(V, nu)
        return _bem_rigid(A, Pc, area, relaxed)


def disk_cells(n: int, R: float = 1.0) -> np.ndarray:
    """Polar mesh of a disk: n rings of width R/n, ring k with round(2 pi (k + 1/2)) cells (4 quarter
    cells in the centre)."""
    V, dr = [], R / n
    for k in range(n):
        r0, r1 = k * dr, (k + 1) * dr
        m = 4 if k == 0 else int(round(2 * math.pi * (k + 0.5)))
        th = np.linspace(0, 2 * math.pi, m + 1)
        for a in range(m):
            t0, t1 = th[a], th[a + 1]
            if k == 0:
                tm = 0.5 * (t0 + t1)
                V.append([[0, 0], [r1 * math.cos(t0), r1 * math.sin(t0)], [r1 * math.cos(tm), r1 * math.sin(tm)],
                          [r1 * math.cos(t1), r1 * math.sin(t1)]])
            else:
                V.append([[r0 * math.cos(t0), r0 * math.sin(t0)], [r1 * math.cos(t0), r1 * math.sin(t0)],
                          [r1 * math.cos(t1), r1 * math.sin(t1)], [r0 * math.cos(t1), r0 * math.sin(t1)]])
    return np.asarray(V, float)


def square_cells(n: int, B_: float = 1.0) -> np.ndarray:
    g = np.linspace(-B_, B_, n + 1)
    return np.asarray([[[g[i], g[j]], [g[i + 1], g[j]], [g[i + 1], g[j + 1]], [g[i], g[j + 1]]]
                       for j in range(n) for i in range(n)], float)


def bem_reference(shape: str, ns: Tuple[int, int], nu: float = NU, relaxed: bool = False) -> np.ndarray:
    """Mesh-extrapolated (p = 1) static BEM stiffness in units of G a, G a^3 (a = R or half-width B)."""
    cells = disk_cells if shape == "disk" else square_cells
    K = [static_bem(cells(n), nu, relaxed) for n in ns]
    return richardson([1.0 / ns[0], 1.0 / ns[1]], K)


# ======================================================================================
# Independent reference 2: dynamic relaxed BEM with the exact Lamb kernels
# ======================================================================================
def _wavenumber_contour(t1: float = 2.0, delta: float = 0.25, X: float = 100.0, n1: int = 2400, dx: float = 0.05):
    """Integration path xi = t + i delta sin(pi t/t1) (t <= t1), then the real axis to X: passes above the
    branch points and the Rayleigh pole of the elastic half-space (time factor e^{+iwt}).  The
    integrands of :func:`lamb_kernels` decay like xi^-2 beyond the pole: truncation at X = 100 changes
    the kernels by less than 1e-4 (checked against X = 400)."""
    g, w = np.polynomial.legendre.leggauss(8)

    def panels(a, b, step):
        e = np.arange(a, b + 1e-12, step)
        t = (0.5 * np.diff(e)[:, None] * (g[None] + 1) + e[:-1, None]).ravel()
        return t, (0.5 * np.diff(e)[:, None] * w[None]).ravel()
    t, wt = panels(0.0, t1, t1 / (n1 // 8))
    xi1 = t + 1j * delta * np.sin(np.pi * t / t1)
    dxi1 = (1 + 1j * delta * (np.pi / t1) * np.cos(np.pi * t / t1)) * wt
    t2, w2 = panels(t1, X, dx)
    return np.concatenate([xi1, t2 + 0j]), np.concatenate([dxi1, w2 + 0j])


def _lamb_symbols(xi: np.ndarray, nu: float):
    """Wavenumber-domain surface flexibilities (G = 1, k_s = 1): vertical g_zz, in-plane horizontal
    g_rr (P-SV) and anti-plane g_tt (SH); principal square-root branch (Re >= 0)."""
    eta = math.sqrt(2 * (1 - nu) / (1 - 2 * nu))
    vp = np.sqrt(xi * xi - 1.0 / eta ** 2 + 0j)
    vs = np.sqrt(xi * xi - 1.0 + 0j)
    F = (2 * xi * xi - 1) ** 2 - 4 * xi * xi * vp * vs
    return -vp / F, -vs / F, 1.0 / vs


def _lamb_coupling_symbol(xi: np.ndarray, nu: float) -> np.ndarray:
    """P-SV cross symbol g_c (G = 1, k_s = 1): radial surface displacement per unit upward vertical
    surface traction J0(xi r).  Large xi: g_c -> (1 - 2 nu)/(2 xi), the Boussinesq inward/outward
    term (1 - 2 nu)/(4 pi G r) of :func:`surface_kernel`."""
    eta = math.sqrt(2 * (1 - nu) / (1 - 2 * nu))
    vp = np.sqrt(xi * xi - 1.0 / eta ** 2 + 0j)
    vs = np.sqrt(xi * xi - 1.0 + 0j)
    F = (2 * xi * xi - 1) ** 2 - 4 * xi * xi * vp * vs
    return -xi * (2 * xi * xi - 1 - 2 * vp * vs) / F


def _bessel_012(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """J0, J1, J2 of a real or complex argument.  Real arguments (the real-axis part of the
    wavenumber contour, ~90 % of the points) use the fast real routines and J2 = 2 J1/x - J0
    (series x^2/8 - x^4/96 near 0); complex arguments use the general routine."""
    if np.iscomplexobj(x):
        return special.jv(0, x), special.jv(1, x), special.jv(2, x)
    J0, J1 = special.j0(x), special.j1(x)
    small = np.abs(x) < 1e-3
    xs = np.where(small, 1.0, x)
    J2 = np.where(small, x * x / 8 - x ** 4 / 96, 2 * J1 / xs - J0)
    return J0, J1, J2


def lamb_kernels(a0: float, r, nu: float = NU, coupling: bool = False, **contour):
    """Dynamic parts (exact minus static) of the half-space surface Green's functions at distances r
    (units G = 1, length R, a0 = w R / Vs):

        G_zz = (1-nu)/(2 pi r) + Dz,   G_xx = A + B cos 2psi,  G_yy = A - B cos 2psi,  G_xy = B sin 2psi,
        A = (2-nu)/(4 pi r) + DA,      B = nu/(4 pi r) + DB,
        G_xz = Gc dx/r,  G_yz = Gc dy/r,  G_zx = -G_xz,  G_zy = -G_yz,  Gc = (1-2nu)/(4 pi r) + Dc

    (G_ab = u_a at the observation point per unit force P_b at the source, dx = obs - source, z up;
    the static parts are :func:`surface_kernel`).  Hankel transforms of the Lamb symbols minus their
    static limits, integrated on :func:`_wavenumber_contour` (all remainders are smooth in r).
    Returns (Dz, DA, DB), and Dc as a fourth array when ``coupling`` (welded contact needs it)."""
    xi, dxi = _wavenumber_contour(**contour)
    gzz, grr, gtt = _lamb_symbols(xi, nu)
    f = {"z": (xi * gzz - (1 - nu)) * dxi, "a": (xi * (grr + gtt) - (2 - nu)) * dxi,
         "b": (xi * (grr - gtt) + nu) * dxi}
    if coupling:
        f["c"] = (xi * _lamb_coupling_symbol(xi, nu) - (1 - 2 * nu) / 2) * dxi
    lifted = xi.imag != 0                         # complex part of the contour; the rest is the real axis
    parts = [(xi[lifted], {k: v[lifted] for k, v in f.items()}),
             (xi[~lifted].real, {k: v[~lifted] for k, v in f.items()})]
    r = np.asarray(r, float).reshape(-1)
    out = {k: np.zeros(r.size, complex) for k in f}
    with np.errstate(all="ignore"):
        for i0 in range(0, r.size, 32):
            for x, fk in parts:
                J0, J1, J2 = _bessel_012(a0 * x[None] * r[i0:i0 + 32, None])
                out["z"][i0:i0 + 32] += a0 / (2 * np.pi) * (J0 @ fk["z"])
                out["a"][i0:i0 + 32] += a0 / (4 * np.pi) * (J0 @ fk["a"])
                out["b"][i0:i0 + 32] += -a0 / (4 * np.pi) * (J2 @ fk["b"])
                if coupling:
                    out["c"][i0:i0 + 32] += a0 / (2 * np.pi) * (J1 @ fk["c"])
    return tuple(out[k] for k in ("z", "a", "b", "c") if k in out)


#: spacing of the distance table of the dynamic kernels, in units of the foundation size (~R/136)
KERNEL_TABLE_STEP = 2.2 / 299


def kernel_table(a0: float, rmax: float, nu: float = NU, coupling: bool = False) -> List[CubicSpline]:
    """Cubic splines of the dynamic kernel remainders of :func:`lamb_kernels` on 0 <= r <= rmax (the
    remainders are smooth, so a table plus spline is exact to ~1e-6 and far cheaper than evaluating
    the wavenumber integrals at every quadrature-point pair).  ``rmax`` must cover the largest
    source-observation distance: outside the table a spline extrapolates."""
    n = max(300, int(math.ceil(rmax / KERNEL_TABLE_STEP)) + 1)
    rg = np.linspace(0.0, float(rmax), n)
    return [CubicSpline(rg, v) for v in lamb_kernels(a0, rg, nu, coupling=coupling)]


def max_pair_distance(P: np.ndarray, Q: np.ndarray) -> float:
    """Upper bound of |p - q| over two plane point sets (|p - q| <= |p| + |q|, both about the origin),
    with a 0.1 % margin: the range of :func:`kernel_table` for a foundation centred at the origin
    (2 R for a disk, 2 sqrt(2) B for a square)."""
    rp = float(np.hypot(P[..., 0], P[..., 1]).max())
    rq = float(np.hypot(Q[..., 0], Q[..., 1]).max())
    return 1.001 * (rp + rq)


def dynamic_bem_impedance(V: np.ndarray, a0: float, nu: float = NU, nq: int = 4, relaxed: bool = True) -> np.ndarray:
    """Dynamic impedance (6x6, G = 1, length R) of a rigid surface foundation on the cells V at ``a0``,
    relaxed (default) or welded contact: static singular part as :func:`static_bem`, smooth dynamic
    remainders of :func:`lamb_kernels` by Gauss quadrature (``nq`` x ``nq`` points per cell).  Welded
    contact adds the vertical/horizontal coupling kernel G_xz, which couples sliding and rocking."""
    nc = V.shape[0]
    with np.errstate(all="ignore"):
        A0, Pc, area = bem_influence(V, nu)
        Q, W = _cell_quadrature(V, nq)
        spl = kernel_table(a0, max_pair_distance(Pc, Q), nu, coupling=not relaxed)
        A = A0.reshape(nc, 3, nc, 3).astype(complex)
        for i0 in range(0, nc, 32):
            d = Pc[i0:i0 + 32, None, None, :] - Q[None]
            r = np.hypot(d[..., 0], d[..., 1])
            psi = np.arctan2(d[..., 1], d[..., 0])
            z, a, b = (s(r) for s in spl[:3])
            c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
            A[i0:i0 + 32, 2, :, 2] += np.einsum("jq,ijq->ij", W, z)
            A[i0:i0 + 32, 0, :, 0] += np.einsum("jq,ijq->ij", W, a + b * c2)
            A[i0:i0 + 32, 1, :, 1] += np.einsum("jq,ijq->ij", W, a - b * c2)
            bxy = np.einsum("jq,ijq->ij", W, b * s2)
            A[i0:i0 + 32, 0, :, 1] += bxy
            A[i0:i0 + 32, 1, :, 0] += bxy
            if not relaxed:                      # G_xz = Gc cos(psi), G_zx = -G_xz (reciprocity)
                gc = spl[3](r)
                cx = np.einsum("jq,ijq->ij", W, gc * np.cos(psi))
                cy = np.einsum("jq,ijq->ij", W, gc * np.sin(psi))
                A[i0:i0 + 32, 0, :, 2] += cx
                A[i0:i0 + 32, 1, :, 2] += cy
                A[i0:i0 + 32, 2, :, 0] -= cx
                A[i0:i0 + 32, 2, :, 1] -= cy
        return _bem_rigid(A.reshape(3 * nc, 3 * nc), Pc, area, relaxed=relaxed)


def relaxed_bem_impedance(V: np.ndarray, a0: float, nu: float = NU, nq: int = 4) -> np.ndarray:
    """Relaxed-contact dynamic impedance (6x6, G = 1, length R) of a rigid surface foundation on the
    cells V at ``a0`` (:func:`dynamic_bem_impedance` with ``relaxed=True``)."""
    return dynamic_bem_impedance(V, a0, nu, nq, relaxed=True)


def welded_bem_impedance(V: np.ndarray, a0: float, nu: float = NU, nq: int = 4) -> np.ndarray:
    """Welded-contact dynamic impedance (6x6, G = 1, length R) of a rigid surface foundation on the
    cells V at ``a0`` (:func:`dynamic_bem_impedance` with ``relaxed=False``): the rigorous counterpart
    of SASSI's K_G = T^T X_ff T, including the dynamic sliding-rocking coupling."""
    return dynamic_bem_impedance(V, a0, nu, nq, relaxed=False)


def lowfreq_damping(nu: float = NU) -> Dict[str, float]:
    """Exact low-frequency radiation damping of a rigid disk (relaxed contact) on an elastic half-space:
    first-order perturbation ``Im C = (1/4 pi^2) Int Im G~(kappa) |t~(kappa)|^2 d2kappa`` with the
    static punch tractions (t~ -> P, i kappa_x M) and the exact Lamb symbols (body waves + Rayleigh
    pole).  Returns c_v(0), c_h(0) and A_r with c_r ~ A_r a0^2 (normalised as R2 C.2)."""
    eta = math.sqrt(2 * (1 - nu) / (1 - 2 * nu))

    def F(x):
        return (2 * x * x - 1) ** 2 - 4 * x * x * np.sqrt(x * x - 1 / eta ** 2 + 0j) * np.sqrt(x * x - 1 + 0j)

    def im_g(which, x):                          # imaginary parts on the real axis, x < 1 (outgoing branch)
        vp = 1j * np.sqrt(1 / eta ** 2 - x * x + 0j) if x < 1 / eta else np.sqrt(x * x - 1 / eta ** 2 + 0j)
        vs = 1j * np.sqrt(1 - x * x + 0j)
        Fx = (2 * x * x - 1) ** 2 - 4 * x * x * vp * vs
        return {"zz": -vp / Fx, "rr": -vs / Fx, "tt": 1.0 / vs}[which].imag

    xr = optimize.brentq(lambda x: F(x).real, 1.0 + 1e-6, 1.5)
    dF = (F(xr + 1e-7) - F(xr - 1e-7)).real / 2e-7
    res = {"zz": (-np.sqrt(xr * xr - 1 / eta ** 2) / dF), "rr": (-np.sqrt(xr * xr - 1) / dF)}

    def integral(which, power):
        body = sum(integrate.quad(lambda x: x ** power * im_g(which, x), a, b, limit=200)[0]
                   for a, b in ((0, 1 / eta), (1 / eta, 1.0)))
        if which == "tt":
            return body
        return body - math.pi * abs(res[which]) * xr ** power
    Iv = integral("zz", 1)
    Ir = integral("zz", 3)
    Ih = integral("rr", 1) + integral("tt", 1)
    return {"c_v0": -2.0 / (math.pi * (1 - nu)) * Iv, "c_h0": -2.0 / (math.pi * (2 - nu)) * Ih,
            "A_r": -2.0 / (3 * math.pi * (1 - nu)) * Ir}


# ======================================================================================
# R2 F.2: exact 3-DOF model of an SDOF on a rigid massless foundation
# ======================================================================================
def three_dof(omega: float, m: float, h: float, kstar: complex, S: np.ndarray, m0: float = 0.0,
              I0: float = 0.0) -> Tuple[complex, complex, complex]:
    """R2 F.2 generalised to a coupled impedance ``S = [[S_xx, S_xt], [S_tx, S_tt]]`` (translation and
    rocking theta about the horizontal axis, u_x(h) = u0 + h theta).  Unknowns relative to the free field
    for a unit control *displacement* (u_g = 1, ddu_g = -w^2): structural deformation u, foundation
    translation u0, rocking theta.  The total motion of the mass per unit control motion is
    ``1 + u0 + h theta + u`` (= the FILE8 transfer function)."""
    w2 = omega * omega
    A = np.array([[kstar - w2 * m, -w2 * m, -w2 * m * h],
                  [-w2 * m, S[0, 0] - w2 * (m + m0), S[0, 1] - w2 * m * h],
                  [-w2 * m * h, S[1, 0] - w2 * m * h, S[1, 1] - w2 * (m * h * h + I0)]], complex)
    b = w2 * np.array([m, m + m0, m * h], complex)          # -m ddu_g with ddu_g = -w^2
    u, u0, th = np.linalg.solve(A, b)
    return u, u0, th


def _total_motion(omega: float, m: float, h: float, kstar: complex, S: np.ndarray) -> complex:
    """Total motion of the SDOF mass per unit control motion, ``1 + u0 + h theta + u`` (R2 F.2)."""
    u, u0, th = three_dof(omega, m, h, kstar, S)
    return 1 + u0 + h * th + u


# ======================================================================================
# VP-10  rigid disk statics
# ======================================================================================
VP10_MESHES = (8, 12, 16)


def _disk_runs(workdir: Path, R: float, meshes: Sequence[int], a0s: Sequence[float], depth: float, tmax: float,
               shape: str = "disk", tag: str = "d") -> Dict[int, dict]:
    out = {}
    for nd in meshes:
        h = R / nd if shape == "disk" else 2 * R / nd
        site = soil_column(h, tmax, depth)
        fs, a0 = a0_frequencies(a0s, R)
        mdl = B.surface_rigid_mat(site, half_width=R, ndiv=nd, shape=shape)
        res = impedance_chain(Path(workdir) / f"{tag}{nd}", site, fs, mdl)
        res.update(h=h, a0=a0, nint=len(mdl["mat"]), nlayers=len(site.layers), R0=mdl.rad)
        out[nd] = res
    return out


def _check_kg_identity(r: VPResult, runs: Dict[int, dict], label: str) -> None:
    for nd, res in runs.items():
        kg_checks(r, f"{label} mesh {nd}", res["KG"], res["KD"])


@problem("VP-10", "Rigid disk static stiffness on a half-space (welded mat, mesh-extrapolated)", tier="P0",
         modules=["SITE", "POINT", "HOUSE", "ANALYS"], source="R2 B.1; R1 §3.4; D-ANL-07")
def vp10(workdir: Path) -> VPResult:
    """Disk R = 5 m on a uniform half-space (Vs = 100 m/s, rho = 2, nu = 1/3, beta = 1e-4) at
    a0 = 0.05; ring meshes with 8, 12 and 16 rings (h = R/n, R0 = 0.85 h), top sublayer h growing by 12 %
    to 4R, 30 R deep + 20 generated sublayers; K_G from ANALYS <impe> = 2; extrapolated (p = 1)."""
    r = VPResult()
    R = 5.0
    runs = _disk_runs(Path(workdir), R, VP10_MESHES, [0.05], depth=30 * R, tmax=4 * R)
    _check_kg_identity(r, runs, "disk")
    hs = [runs[n]["h"] for n in VP10_MESHES]
    unit = {"Kv": (2, GMOD * R), "Kh": (0, GMOD * R), "Kr": (3, GMOD * R ** 3), "Kt": (5, GMOD * R ** 3)}
    seq = {k: [runs[n]["KG"][0][i, i].real / u for n in VP10_MESHES] for k, (i, u) in unit.items()}
    rel = {k: [runs[n]["KR"][0][i, i].real / u for n in VP10_MESHES] for k, (i, u) in unit.items()}
    ext = {k: float(richardson(hs, v)) for k, v in seq.items()}
    ext_rel = {k: float(richardson(hs, v)) for k, v in rel.items()}
    for k, ref in VP10_REF.items():
        r.check(f"{k} (welded K_G, extrapolated) vs R2 B.1 relaxed closed form", ext[k], ref, rtol=0.05,
                note=f"meshes {', '.join(f'{v:.4f}' for v in seq[k])}; observed order "
                     f"{observed_order(hs, seq[k]):.2f}")
    for k in ("Kv", "Kr", "Kt"):
        r.check(f"{k} (relaxed contact from F_ff, extrapolated) vs the exact relaxed value", ext_rel[k], VP10_REF[k],
                rtol=0.05, note=f"meshes {', '.join(f'{v:.4f}' for v in rel[k])}")
    K0 = runs[VP10_MESHES[-1]]["KG"][0].real
    r.check("isotropy of the ring mesh: K_xx = K_yy (finest mesh)", abs(K0[0, 0] - K0[1, 1]) / K0[0, 0], 0.0,
            atol=1e-8)
    # independent welded reference: exact-kernel BEM (derived)
    Kb = bem_reference("disk", (16, 24))
    bem = {"Kv": Kb[2, 2], "Kh": Kb[0, 0], "Kr": Kb[3, 3], "Kt": Kb[5, 5]}
    for k in VP10_REF:
        r.check(f"{k} (welded, extrapolated) vs welded-contact BEM with exact kernels (derived reference)", ext[k],
                bem[k], rtol=0.05)
    mos = 4 * math.log(3 - 4 * NU) / (1 - 2 * NU)
    r.notes.append("Welded vs relaxed contact: K_G = T^T X_ff T ties every translation of the mat nodes to the rigid "
                   f"body (welded).  Welded references (BEM): Kv {bem['Kv']:.4f} (Mossakovskii exact {mos:.4f}), "
                   f"Kh {bem['Kh']:.4f}, Kr {bem['Kr']:.4f}, Kt {bem['Kt']:.4f}; relaxed: 6, 4.8, 4, 5.3333.")
    for k in VP10_REF:
        r.notes.append(f"{k}: meshes n = {VP10_MESHES}: {', '.join(f'{v:.4f}' for v in seq[k])} -> extrapolated "
                       f"{ext[k]:.4f} ({100 * (ext[k] / VP10_REF[k] - 1):+.2f} % vs R2 B.1, "
                       f"{100 * (ext[k] / bem[k] - 1):+.2f} % vs welded BEM); relaxed from F_ff "
                       f"{ext_rel[k]:.4f}")
    r.notes.append("Mesh: " + "; ".join(f"n = {n}: {runs[n]['nint']} nodes, h = {runs[n]['h']:.4g} m, R0 = "
                                         f"{runs[n]['R0']:.4g} m, {runs[n]['nlayers']} sublayers"
                                         for n in VP10_MESHES))
    return r


# ======================================================================================
# VP-11  rigid disk dynamic impedance
# ======================================================================================
VP11_MESH = 12
VP11_A0 = (0.05, 0.5, 1.0, 1.5)
VP11_A0_HIGH = (6.0, 7.0, 8.0)


def _coefficients(K: np.ndarray, K0: np.ndarray, a0: float) -> Dict[str, float]:
    """k = Re K / K0 and c = Im K / (a0 K0) of the diagonal terms (R2 §0 normalisation)."""
    out = {}
    for name, i in (("h", 0), ("v", 2), ("r", 3), ("t", 5)):
        out["k" + name] = float(K[i, i].real / K0[i, i])
        out["c" + name] = float(K[i, i].imag / (a0 * K0[i, i]))
    return out


@problem("VP-11", "Rigid disk dynamic impedance (Veletsos-Verbic a0 <= 1.5; high-frequency dashpots)", tier="P0",
         modules=["SITE", "POINT", "HOUSE", "ANALYS"], source="R2 C.1, C.2; D-ANL-07")
def vp11(workdir: Path) -> VPResult:
    """Disk R = 5 m, 12 rings (h = R/12, R0 = 0.85 h), near-elastic half-space (beta = 1e-4).
    Low band a0 = 0.5, 1.0, 1.5 normalised by the a0 = 0.05 stiffness of the same model; top sublayer h
    growing to lambda_s(a0 = 1.5)/8, 20 R deep.  High band a0 = 6, 7, 8: sublayers <= lambda_s(8)/8,
    3 R deep (the 20 generated sublayers carry the radiation); Im K / w vs the plane-wave dashpots."""
    r = VPResult()
    R = 5.0
    h = R / VP11_MESH
    wd = Path(workdir)
    lam = 2 * math.pi * R / max(VP11_A0)
    site = soil_column(h, min(lam / 8, 4 * R), 20 * R)
    fs, a0 = a0_frequencies(VP11_A0, R)
    mdl = B.surface_rigid_mat(site, half_width=R, ndiv=VP11_MESH, shape="disk")
    low = impedance_chain(wd / "low", site, fs, mdl)
    kg_checks(r, "low band", low["KG"], low["KD"])
    K0w, K0r = low["KG"][0].real, low["KR"][0].real
    worst_vv = {}
    for q in range(1, len(a0)):
        x = float(a0[q])
        cw = _coefficients(low["KG"][q], K0w, x)
        cr = _coefficients(low["KR"][q], K0r, x)
        kr_, cr_, kv_, cv_ = VP11_VV[round(x, 2)]
        for name, val, ref in (("k_h", cw["kh"], 1.0), ("c_h", cw["ch"], VP11_ALPHA1), ("k_r", cw["kr"], kr_),
                               ("c_r", cw["cr"], cr_), ("k_v", cw["kv"], kv_), ("c_v", cw["cv"], cv_)):
            if name == "c_r":    # D-W2-01: superseded by the exact-kernel relaxed BEM check below
                r.inform(f"{name} at a0 = {x:g} (K_G, welded) vs Veletsos-Verbic (R2 C.2)", val, ref,
                         note="informative: Veletsos-Verbic rocking damping is 28-53 % below the rigorous solution")
            else:
                r.check(f"{name} at a0 = {x:g} (K_G, welded) vs Veletsos-Verbic (R2 C.2)", val, ref, rtol=0.10)
            worst_vv[(name, x)] = val / ref - 1
        bem = VP11_BEM[round(x, 2)]
        for name in ("kh", "ch", "kr", "cr", "kv", "cv", "kt", "ct"):
            r.check(f"{name[0]}_{name[1]} at a0 = {x:g} (relaxed contact from F_ff) vs relaxed BEM with exact Lamb "
                    f"kernels (derived reference)", cr[name], bem[name], rtol=0.10)
        r.notes.append(f"a0 = {x:g}: welded k_h {cw['kh']:.4f} c_h {cw['ch']:.4f} k_r {cw['kr']:.4f} c_r {cw['cr']:.4f} "
                       f"k_v {cw['kv']:.4f} c_v {cw['cv']:.4f} k_t {cw['kt']:.4f} c_t {cw['ct']:.4f} | relaxed k_h "
                       f"{cr['kh']:.4f} c_h {cr['ch']:.4f} k_r {cr['kr']:.4f} c_r {cr['cr']:.4f} k_v {cr['kv']:.4f} "
                       f"c_v {cr['cv']:.4f} | BEM c_h {bem['ch']:.4f} c_r {bem['cr']:.4f} | VV c_h 0.65 c_r {cr_:.4f}")
    # high-frequency dashpots
    hf_h = R / VP11_MESH
    lam8 = 2 * math.pi * R / max(VP11_A0_HIGH)
    site_hf = soil_column(hf_h, lam8 / 8, 3 * R)
    fs_hf, a0_hf = a0_frequencies(VP11_A0_HIGH, R)
    mdl_hf = B.surface_rigid_mat(site_hf, half_width=R, ndiv=VP11_MESH, shape="disk")
    high = impedance_chain(wd / "high", site_hf, fs_hf, mdl_hf, relaxed=False)
    kg_checks(r, "high band", high["KG"], high["KD"])
    area, I, J = math.pi * R ** 2, math.pi * R ** 4 / 4, math.pi * R ** 4 / 2
    dash = {"h": (0, RHO * VS * area), "v": (2, RHO * VPP * area), "r": (3, RHO * VPP * I), "t": (5, RHO * VS * J)}
    for q, x in enumerate(a0_hf):
        w = x * VS / R
        vals = {k: high["KG"][q][i, i].imag / w / C for k, (i, C) in dash.items()}
        for k, v in vals.items():
            r.check(f"c_{k}(a0 = {x:g}) / c(inf): Im K / w vs the plane-wave dashpot (R2 C.1, c(inf) = "
                    f"{VP11_CINF[k]:.3f})", v, 1.0, rtol=0.05)
        r.notes.append(f"a0 = {x:g}: Im K/(w C_inf) h {vals['h']:.4f} v {vals['v']:.4f} r {vals['r']:.4f} "
                       f"t {vals['t']:.4f}")
    lf = lowfreq_damping(NU)
    vv_vs_bem = "; ".join(f"a0 = {x:g}: c_r {100 * (VP11_VV[x][1] / VP11_BEM[x]['cr'] - 1):+.0f} %, c_h "
                          f"{100 * (VP11_ALPHA1 / VP11_BEM[x]['ch'] - 1):+.1f} %, k_r "
                          f"{100 * (VP11_VV[x][0] / VP11_BEM[x]['kr'] - 1):+.1f} %, k_v "
                          f"{100 * (VP11_VV[x][2] / VP11_BEM[x]['kv'] - 1):+.1f} %, c_v "
                          f"{100 * (VP11_VV[x][3] / VP11_BEM[x]['cv'] - 1):+.1f} %" for x in sorted(VP11_VV))
    r.notes.append(
        "Reference consistency (evidence, docs/verification/impedance_study.md): the exact low-frequency radiation "
        f"damping of the relaxed disk (perturbation with the exact Lamb kernel) is c_h(0) = {lf['c_h0']:.3f}, c_v(0) = "
        f"{lf['c_v0']:.3f} and c_r = {lf['A_r']:.3f} a0^2; Veletsos-Verbic use c_h = 0.65 and c_r = 0.1 a0^2.  "
        "Veletsos-Verbic vs the rigorous relaxed BEM (VP11_BEM, exact Lamb kernels): " + vv_vs_bem + ".  SASSI-EDU "
        "(relaxed contact) agrees with the BEM within a few %.")
    bad = [f"{n} at a0 = {x:g}: {100 * e:+.1f} %" for (n, x), e in worst_vv.items() if abs(e) > 0.10]
    if bad:
        r.notes.append("Veletsos-Verbic checks outside 10 %: " + "; ".join(bad))
    return r


# ======================================================================================
# VP-13  Tajirian & Tabatabaie (1985): disk on a layer over a rigid base
# ======================================================================================
VP13_RES_V = (math.pi / 6, math.pi / 2, 5 * math.pi / 6)            # (2n-1) Vp/(4H) -> A0 = 0.524, 1.571, 2.618
VP13_RES_H = (math.pi / 12, math.pi / 4, 5 * math.pi / 12)          # (2n-1) Vs/(4H) -> A0 = 0.262, 0.785, 1.309
VP13_DA = 0.01
VP13_STATIC_MESHES = ((4, 24), (8, 48), (12, 72))


def _tt_site(nsub: int, beta: float) -> B.Site:
    lay = [B.SoilLayer(3.0 / nsub, 1.0, 2.0, 1.0, beta)] * nsub
    return B.Site(lay, B.SoilLayer(0.0, 1.0, 2.0, 1.0, beta), nl=0)


def _tt_run(wd: Path, nd: int, nsub: int, beta: float, a0max: float, nsample: int = 3) -> Dict[str, object]:
    """One SITE + POINT + HOUSE + ANALYS sweep A0 = 0.01 .. a0max of the Tajirian-Tabatabaie disk.

    Returns ``a0``, ``KG`` (nF, 6, 6) and ``KD``, ``KGs``: K_G and the independent T^T F_ff^-1 T at
    ``nsample`` frequencies spread over the sweep (first, middle, last) for :func:`kg_checks`; the
    intermediate files are deleted afterwards (:func:`purge_work_files`)."""
    R = 0.5
    site = _tt_site(nsub, beta)
    fs, a0 = a0_frequencies(np.arange(VP13_DA, a0max + 1e-9, VP13_DA), R, da=VP13_DA, vs=1.0)
    mdl = B.surface_rigid_mat(site, half_width=R, ndiv=nd, shape="disk")
    files = B.run_soil(wd, "m", site, fs, layer=0, rad=mdl.rad)
    B.run_house(wd, "m", mdl)
    KG = analys_global_impedance(wd, "m", fs)
    purge_work_files(wd)
    nodes = list(mdl["mat"])
    xy = np.array([mdl.house.xyz(n)[:2] for n in nodes])
    T = rigid_transform(xy)
    qs = sorted({int(round(x)) for x in np.linspace(0, len(a0) - 1, nsample)})
    KD = np.array([independent_global_impedance(flexibility_matrix(files["FILE3"], q, xy, np.ones(len(nodes), int)), T)
                   for q in qs])
    return {"a0": a0, "KG": KG, "KD": KD, "KGs": KG[qs]}


def _zgv_onset(wd: Path) -> Tuple[float, float]:
    """A0 and k R at which the second Rayleigh mode of the (lightly damped) layer starts to propagate
    near the thickness-stretch cut-off Vp/(4H) (SITE Mode 1 only, 24 sublayers, beta = 1e-4)."""
    site = _tt_site(24, 1e-4)
    da = 0.001
    fs, a0 = a0_frequencies(np.arange(0.45, 0.56, da), 0.5, da=da, vs=1.0)
    B.write_deck(wd, "m", B.site_deck(site, fs, mode2=0, model="m"))
    B.run("SITE", wd, "m")
    f2 = read_container(Path(wd) / "FILE2", "FILE2")
    purge_work_files(wd)
    for q in range(len(a0)):
        k = np.asarray(f2["kR"][q])
        prop = k[(k.real > 0) & (np.abs(k.imag) <= 3e-3 * np.abs(k.real))]
        if prop.size >= 2:
            return float(a0[q]), float(np.sort(prop.real)[0] * 0.5)
    return float("nan"), float("nan")


def psv_layer_system(k, w: float, vs: float, vp: float, rho: float = 1.0) -> np.ndarray:
    """System matrix A (..., 4, 4) of the P-SV equations y' = A y of a homogeneous elastic layer
    (continuum, independent of the thin-layer discretisation of SITE).  z down, plane waves
    e^{i(wt - kx)}; y = [U, W, T, S] with U and T the horizontal displacement and shear traction
    divided by i, W the vertical displacement and S the normal traction on horizontal planes, so that
    A is real for real k and w:

        U' = k W + T/mu,                 W' = -k lam/M U + S/M,
        T' = (k^2 (M - lam^2/M) - rho w^2) U + k lam/M S,    S' = -rho w^2 W - k T

    (mu = rho Vs^2, M = lam + 2 mu = rho Vp^2).  At k = 0 it splits into the shear and the
    compression column."""
    mu, M = rho * vs ** 2, rho * vp ** 2
    lam = M - 2 * mu
    k = np.asarray(k, float)
    A = np.zeros(k.shape + (4, 4))
    A[..., 0, 1], A[..., 0, 2] = k, 1.0 / mu
    A[..., 1, 0], A[..., 1, 3] = -k * lam / M, 1.0 / M
    A[..., 2, 0], A[..., 2, 3] = k * k * (M - lam * lam / M) - rho * w * w, k * lam / M
    A[..., 3, 1], A[..., 3, 2] = -rho * w * w, -k
    return A


def layer_dispersion_det(k, w: float, H: float, vs: float, vp: float, rho: float = 1.0) -> np.ndarray:
    """Secular function of the P-SV (Rayleigh-type) modes of a layer with a free top and a fixed
    base: the determinant of the displacement-to-displacement block of the propagator expm(A H) (a
    free top has T = S = 0, so [U, W](H) = P_uu [U, W](0) must vanish for a mode).  Real for real k;
    an entire function of k (no poles), so its sign changes are the real wavenumbers."""
    P = linalg.expm(psv_layer_system(k, w, vs, vp, rho) * H)
    return P[..., 0, 0] * P[..., 1, 1] - P[..., 0, 1] * P[..., 1, 0]


def layer_zgv_onset(H: float = 3.0, vs: float = 1.0, vp: float = 2.0, R: float = 0.5, a0_lo: float = 0.48,
                    kmax: float = 2.0, nk: int = 2500) -> Tuple[float, float]:
    """Continuum reference for the VP-13 vertical evidence: A0 = w R/Vs at which a pair of real
    wavenumbers of the second P-SV mode is born below the k = 0 thickness-stretch cut-off
    Vp/(4H) (zero-group-velocity point), and k R there.  Real roots counted by sign changes of
    :func:`layer_dispersion_det` on 0 < k <= kmax (the pair is resolved once it is 2 kmax/nk apart,
    ~1e-7 in A0 after the true onset); onset bracketed on an A0 grid from ``a0_lo`` (step 0.0025) and
    bisected.  VP-13 layer: A0 = 0.5094, k R = 0.170."""
    ks = np.linspace(1e-6, kmax, nk)

    def roots(A0):
        d = layer_dispersion_det(ks, A0 * vs / R, H, vs, vp)
        i = np.flatnonzero(np.sign(d[1:]) != np.sign(d[:-1]))
        return 0.5 * (ks[i] + ks[i + 1])
    cut = math.pi * vp / (2 * H) * R / vs
    grid = np.arange(a0_lo, cut, 0.0025)
    n0 = roots(grid[0]).size
    hit = [a for a in grid if roots(a).size > n0]
    if not hit:
        return float("nan"), float("nan")
    hi = float(hit[0])
    lo = hi - 0.0025
    for _ in range(20):                                     # 0.0025 / 2^20 ~ 2e-9
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if roots(mid).size > n0 else (mid, hi)
    kr = np.sort(roots(hi)) * R
    gaps = np.diff(kr)
    j = int(np.argmin(gaps))
    return hi, float(0.5 * (kr[j] + kr[j + 1]))


@problem("VP-13", "Tajirian & Tabatabaie (1985): rigid disk on a layer over a rigid base (resonances, statics)",
         tier="P0", modules=["SITE", "POINT", "HOUSE", "ANALYS"], source="R2 C.7, B.4")
def vp13(workdir: Path) -> VPResult:
    """r = 0.5, H = 3, Vs = 1, Vp = 2, rho = 1, beta = 5 % (checked), 15 % and 1 % (diagnostics and
    the vertical-peak-2 guard), 24 sublayers of 0.125 on a rigid
    base (SITE nl = 0); 4 rings (h = 0.125 = sublayer, R0 = 0.85 h); A0 = 0.01 .. 3.0 step 0.01.
    Resonance = local maximum of the compliance amplitude |C_ii|, C = K_G^-1 (within +-10 % of the
    expected value, parabolic refinement).  Statics: |K_G| at A0 = 0.01 (|c(beta)| = 1 in the SASSI
    damping form, so |K(0)| is the elastic stiffness) on 4/8/12 rings with h = sublayer (24/48/72
    sublayers), extrapolated (p = 1).

    Guards of the two vertical checks that fail against the 1-D rule (derived references, so that a
    regression of the very quantities that fail is still detected): the second vertical resonance of
    the lightly damped (beta = 1 %) layer vs (2n-1) Vp/(4H), and the first vertical peak (beta = 5 %)
    against the zero-group-velocity onset of the layer, computed both by SITE and from the continuum
    P-SV dispersion (:func:`layer_zgv_onset`)."""
    r = VPResult()
    wd = Path(workdir)
    curves = {}
    for beta in (0.05, 0.15, 0.01):
        run = _tt_run(wd / f"b{int(round(100 * beta))}", 4, 24, beta, 3.0)
        kg_checks(r, f"beta = {100 * beta:g} % sweep (3 sampled frequencies)", run["KGs"], run["KD"])
        curves[beta] = (run["a0"], run["KG"], np.linalg.inv(run["KG"]))
    a0, KG, C = curves[0.05]
    peaks5 = {}
    for name, i, refs in (("vertical", 2, VP13_RES_V), ("horizontal", 0, VP13_RES_H)):
        for n, ref in enumerate(refs, start=1):
            pk = peak_near(a0, np.abs(C[:, i, i]), ref)
            peaks5[(name, n)] = pk
            if name == "vertical":   # D-W2-02: 1-D rule superseded by the ZGV-onset and beta = 1 % guards
                r.inform(f"{name} resonance {n}: A0 of the |compliance| peak, beta = 5 %", pk, ref,
                         note="informative: the vertical peaks are not at (2n-1)Vp/(4H)" +
                              ("; no peak within +-10 %" if not math.isfinite(pk) else ""))
            else:
                r.check(f"{name} resonance {n}: A0 of the |compliance| peak, beta = 5 %", pk, ref, rtol=0.02,
                        note="no peak within +-10 %" if not math.isfinite(pk) else "")
    found = []
    a015, _, C15 = curves[0.15]
    for name, i, refs in (("vertical", 2, VP13_RES_V), ("horizontal", 0, VP13_RES_H)):
        for ref in refs:
            pk = peak_near(a015, np.abs(C15[:, i, i]), ref)
            found.append(f"{name} {ref:.3f}: {'none' if not math.isfinite(pk) else f'{pk:.3f}'}")
    r.notes.append("beta = 15 %: |compliance| peaks near the expected resonances: " + "; ".join(found))
    a01, _, C1 = curves[0.01]                                      # lightly damped layer
    for name, i in (("vertical", 2), ("horizontal", 0)):
        y = np.abs(C1[:, i, i])
        pk = [peak_near(a01, y, float(a01[j]), window=1e-9) for j in range(1, len(y) - 1)
              if y[j] >= y[j - 1] and y[j] > y[j + 1]]
        r.notes.append(f"beta = 1 % (diagnostic): all |C| peaks of the {name} compliance at A0 = "
                       + ", ".join(f"{v:.3f}" for v in pk))
    pk2 = peak_near(a01, np.abs(C1[:, 2, 2]), VP13_RES_V[1])
    r.check("guard (vertical peak 2): A0 of the |compliance| peak of the lightly damped layer (beta = 1 %) "
            "vs 3 Vp/(4H)", pk2, VP13_RES_V[1], rtol=0.02,
            note="the peak exists; 5 % damping hides it" if math.isfinite(pk2) else "no peak within +-10 %")
    # first vertical peak: the onset of radiation of the layer's second P-SV mode (ZGV point)
    a0z, kz = _zgv_onset(wd / "zgv")
    a0c, kc = layer_zgv_onset()
    r.check("ZGV onset of the layer's 2nd P-SV mode: SITE (24 sublayers) vs the continuum P-SV dispersion "
            "(derived reference)", a0z, a0c, rtol=0.01, note=f"k R: SITE {kz:.3f}, continuum {kc:.3f}")
    pk1 = peaks5[("vertical", 1)]
    r.require("guard (vertical peak 1): the beta = 5 % |compliance| peak lies below the continuum ZGV onset "
              f"(A0 {pk1:.4f} < {a0c:.4f})", math.isfinite(pk1) and pk1 < a0c)
    r.check("guard (vertical peak 1): A0 of the beta = 5 % |compliance| peak vs the continuum ZGV onset "
            "(derived reference: radiation sets in there; damping moves the peak below it)", pk1, a0c, rtol=0.06)
    r.notes.append(
        f"Vertical: the layer's second Rayleigh mode (thickness-stretch family) starts to propagate at A0 = {a0z:.3f} "
        f"(SITE; continuum P-SV dispersion {a0c:.4f}) with k R = {kz:.3f} (continuum {kc:.3f}) > 0, i.e. at a "
        "zero-group-velocity point below the k = 0 cut-off Vp/(4H) (A0 = 0.524): radiation sets in there and, with 5 % "
        f"damping, the vertical compliance peak lies further below ({pk1:.4f}, {100 * (pk1 / a0c - 1):+.1f} % of the "
        "onset).  The '(2n-1)Vp/(4H)' rule of R2 C.7 is a 1-D column frequency, not the compliance peak.")
    # statics
    hs, seq = [], {"Kh": [], "Kr": [], "Kv": [], "Kt": []}
    for nd, nsub in VP13_STATIC_MESHES:
        st = _tt_run(wd / f"s{nd}", nd, nsub, 0.05, VP13_DA, nsample=1)
        kg_checks(r, f"statics mesh {nd}/{nsub}", st["KGs"], st["KD"])
        K = np.abs(st["KG"][0])
        hs.append(0.5 / nd)
        seq["Kh"].append(K[0, 0] / 0.5)
        seq["Kr"].append(K[3, 3] / 0.5 ** 3)
        seq["Kv"].append(K[2, 2] / 0.5)
        seq["Kt"].append(K[5, 5] / 0.5 ** 3)
    ext = {k: float(richardson(hs, v)) for k, v in seq.items()}
    RH = 0.5 / 3.0
    refs = {"Kh": 8 / (2 - NU) * (1 + RH / 2), "Kr": 8 / (3 * (1 - NU)) * (1 + RH / 6)}
    for k, ref in refs.items():
        r.check(f"static {k} (|K_G| at A0 = 0.01, extrapolated) vs Kausel/Elsabee-Morray (R2 B.4)", ext[k], ref,
                rtol=0.05, note=f"meshes {', '.join(f'{v:.4f}' for v in seq[k])}")
    kv_ref, kt_ref = 4 / (1 - NU) * (1 + 1.28 * RH), 16 / 3
    r.notes.append(f"statics [VERIFY formulas, not checked]: Kv {ext['Kv']:.4f} vs 4GR/(1-nu)(1+1.28R/H) = {kv_ref:.4f}; "
                   f"Kt {ext['Kt']:.4f} vs 16GR^3/3 = {kt_ref:.4f}")
    r.notes.append("statics meshes (rings, sublayers) " + ", ".join(f"{a}/{b}" for a, b in VP13_STATIC_MESHES) + ": " +
                   "; ".join(f"{k} {', '.join(f'{v:.4f}' for v in vals)} -> {ext[k]:.4f}" for k, vals in seq.items()))
    return r


# ======================================================================================
# VP-14  square footing statics
# ======================================================================================
VP14_MESHES = (16, 24, 32)


@problem("VP-14", "Square rigid footing statics vs Pais & Kausel (welded mat, mesh-extrapolated)", tier="P0",
         modules=["SITE", "POINT", "HOUSE", "ANALYS"], source="R2 B.2; R1 §3.4")
def vp14(workdir: Path) -> VPResult:
    """Square 2B x 2B, B = 5 m, on the VP-10 half-space at a0 = w B / Vs = 0.05; grids of 16, 24 and 32
    cells per side (h = 2B/n, R0 = 0.90 h), top sublayer h growing 12 % to 4B, 30 B deep; K_G from
    ANALYS <impe> = 2, extrapolated (p = 1)."""
    r = VPResult()
    Bh = 5.0
    runs = _disk_runs(Path(workdir), Bh, VP14_MESHES, [0.05], depth=30 * Bh, tmax=4 * Bh, shape="square", tag="s")
    _check_kg_identity(r, runs, "square")
    hs = [runs[n]["h"] for n in VP14_MESHES]
    idx = {"Kz": (2, 1), "Kx": (0, 1), "Ky": (1, 1), "Kzz": (5, 3), "Kxx": (3, 3), "Kyy": (4, 3)}
    seq = {k: [runs[n]["KG"][0][i, i].real / (GMOD * Bh ** p) for n in VP14_MESHES] for k, (i, p) in idx.items()}
    ext = {k: float(richardson(hs, v)) for k, v in seq.items()}
    Kb = bem_reference("square", (24, 32))
    bem = {k: float(Kb[i, i]) for k, (i, p) in idx.items()}
    for k, ref in VP14_PK.items():
        note = f"meshes {', '.join(f'{v:.4f}' for v in seq[k])}; observed order {observed_order(hs, seq[k]):.2f}"
        if k in ("Kxx", "Kyy"):   # D-W2-03: Pais-Kausel rocking fit superseded by the welded BEM check below
            r.inform(f"{k} (welded K_G, extrapolated) vs Pais & Kausel (R2 B.2)", ext[k], ref,
                     note="informative: the fit is 7 % below the welded exact-kernel BEM; " + note)
        else:
            r.check(f"{k} (welded K_G, extrapolated) vs Pais & Kausel (R2 B.2)", ext[k], ref, rtol=0.05, note=note)
    K = runs[VP14_MESHES[-1]]["KG"][0].real
    r.check("square symmetry K_x = K_y (finest mesh)", abs(K[0, 0] - K[1, 1]) / K[0, 0], 0.0, atol=1e-8)
    r.check("square symmetry K_xx = K_yy (finest mesh)", abs(K[3, 3] - K[4, 4]) / K[3, 3], 0.0, atol=1e-8)
    for k in VP14_PK:
        r.check(f"{k} (welded, extrapolated) vs welded-contact BEM with exact kernels (derived reference)", ext[k],
                bem[k], rtol=0.05)
    for k in VP14_PK:
        r.notes.append(f"{k}: meshes {VP14_MESHES}: {', '.join(f'{v:.4f}' for v in seq[k])} -> {ext[k]:.4f} "
                       f"({100 * (ext[k] / VP14_PK[k] - 1):+.2f} % vs Pais-Kausel {VP14_PK[k]}, "
                       f"{100 * (ext[k] / bem[k] - 1):+.2f} % vs welded BEM {bem[k]:.4f}; Gazetas {VP14_GAZETAS[k]})")
    Kbr = bem_reference("square", (24, 32), relaxed=True)
    r.notes.append(f"Rocking evidence: the exact-kernel BEM gives K_xx = {bem['Kxx']:.3f} GB^3 (welded) and "
                   f"{Kbr[3, 3]:.3f} GB^3 (relaxed); Pais & Kausel 6.0 is below both, Gazetas 5.40/5.58 further below.  The "
                   "same BEM reproduces the exact disk values to 0.1 % (VP-10), and SASSI-EDU's square agrees with it "
                   "within ~1 %: the fitted Pais-Kausel rocking coefficient at L/B = 1 is ~7 % low for welded contact.")
    return r


# ======================================================================================
# VP-17  inertial SSI: SDOF on a rigid massless disk
# ======================================================================================
VP17 = dict(R=10.0, h=10.0, m=942.48, T=0.5, beta=0.05, nd=9, E=3.0e7, rigid=1.0e6)
VP17_REF = dict(Kx=9.60e5, Kt=8.00e7, ratio=1.15805, f_res=1.717, peak=6.71)
#: welded exact-kernel BEM reference of VP-17 (derived): polar meshes (rings) extrapolated with p = 1
#: (8/12 agree with 12/16 within 0.2 %) and the a0 nodes of its cubic interpolation (error < 1e-7)
VP17_BEM_MESHES = (8, 12)
VP17_BEM_A0 = (0.90, 0.96, 1.02, 1.08, 1.14, 1.20, 1.26)


def _vv_impedance(a0: float, R: float) -> np.ndarray:
    """Veletsos-Verbic (R2 C.2) horizontal and rocking impedances, uncoupled (nu = 1/3)."""
    Kx0, Kr0 = 8 * GMOD * R / (2 - NU), 8 * GMOD * R ** 3 / (3 * (1 - NU))
    x = (0.5 * a0) ** 2
    kr, cr = 1 - 0.8 * x / (1 + x), 0.8 * 0.5 * x / (1 + x)
    return np.array([[Kx0 * (1 + 1j * a0 * VP11_ALPHA1), 0], [0, Kr0 * (kr + 1j * a0 * cr)]])


def _bem_impedance(a0: float, R: float) -> np.ndarray:
    """Rigorous relaxed-contact horizontal and rocking impedances (uncoupled): exact static K0 and the
    coefficients of :data:`VP11_BEM`, interpolated linearly in a0 (a0 <= 1.5)."""
    Kx0, Kr0 = 8 * GMOD * R / (2 - NU), 8 * GMOD * R ** 3 / (3 * (1 - NU))
    a = [0.0] + sorted(VP11_BEM)

    def get(key, v0):
        return float(np.interp(a0, a, [v0] + [VP11_BEM[x][key] for x in sorted(VP11_BEM)]))
    ch0 = lowfreq_damping(NU)["c_h0"]
    return np.array([[Kx0 * (get("kh", 1.0) + 1j * a0 * get("ch", ch0)), 0],
                     [0, Kr0 * (get("kr", 1.0) + 1j * a0 * get("cr", 0.0))]])


def welded_disk_impedance(a0s: Sequence[float], meshes: Sequence[int] = VP17_BEM_MESHES) -> np.ndarray:
    """Rigorous welded-contact impedance of a rigid disk on an elastic half-space (len(a0s), 6, 6),
    units G = 1, R = 1: :func:`welded_bem_impedance` on the polar meshes ``meshes`` (rings),
    extrapolated with p = 1 (two meshes)."""
    K = [np.array([welded_bem_impedance(disk_cells(n), float(a)) for a in a0s]) for n in meshes]
    return richardson([1.0 / n for n in meshes], K)


def _resonance(f: np.ndarray, H: np.ndarray, centre: float) -> Tuple[float, float]:
    """Resonance frequency (parabolic refinement) and peak amplitude of |H(f)|."""
    a = np.abs(H)
    return peak_near(f, a, centre, window=0.15), float(a.max())


@problem("VP-17", "Inertial SSI: SDOF on a rigid massless disk vs the exact 3-DOF model and Veletsos-Meek/-Verbic",
         tier="P0", modules=["SITE", "POINT", "HOUSE", "FORCE", "ANALYS"], source="R2 F.1, F.2", slow=True)
def vp17(workdir: Path) -> VPResult:
    """R = 10 m disk (9 rings, R0 = 0.85 h) on the half-space Vs = 100 m/s, rho = 2, nu = 1/3 (beta = 1e-4);
    SDOF: massless Euler-Bernoulli BEAMS cantilever h = 10 m (3EI/h^3 = k) with m = 942.48 t at the top,
    T = 0.5 s, hysteretic beta = 5 % with CMODFORM 1 (k* = k(1 + 2i beta), as R2 F.2).

    1. Impedance S of the soil + mat at the centre node from FORCE unit loads (X, ROTY) on the mat alone
       (ANALYS vibration, 2 load cases): S = C^-1 (it includes the stiff-link flexibility of the mat).
    2. Seismic run (vertical SV) of the SDOF on the mat: FILE8 = the 3-DOF solution fed with S
       (identity, 1e-8).
    3. Worked example (R2 F.1): static K_x, K_theta (ANALYS K_G at 0.05 Hz), the Veletsos-Meek period
       ratio, and the resonance / peak of |u_t/u_g| at the mass vs the values of F.2 with
       Veletsos-Verbic impedances (5 %); the reference values themselves are re-derived.
    4. Guard of the failing peak (derived reference): resonance and peak of F.2 fed with the rigorous
       welded-contact impedance of the disk (exact-kernel BEM incl. the dynamic sliding-rocking
       coupling, :func:`welded_disk_impedance`), 5 %; plus K_G vs an independent T^T F_ff^-1 T and its
       complex symmetry."""
    r = VPResult()
    p = VP17
    R, h, m, beta = p["R"], p["h"], p["m"], p["beta"]
    k = m * (2 * math.pi / p["T"]) ** 2
    kstar = k * (1 + 2j * beta)
    wd = Path(workdir) / "ssi"
    wd.mkdir(parents=True, exist_ok=True)
    df = 0.01
    fnum = [5] + list(range(155, 186))                   # 0.05 Hz and 1.55 .. 1.85 Hz
    fs = B.FrequencySet.harmonic(df, fnum)
    freq = fs.freq
    lam = VS / freq[-1]
    site = soil_column(R / p["nd"], min(lam / 8, 4 * R), 20 * R, cmodform=1)
    # 1. mat alone: impedance at the centre node by unit loads
    E_rigid = p["rigid"] * GMOD                          # mat links 10 x stiffer than the builder default
    mat = B.surface_rigid_mat(site, half_width=R, ndiv=p["nd"], shape="disk", E_rigid=E_rigid)
    c = mat["centre"]
    files = B.run_soil(wd, "m", site, fs, layer=0, rad=mat.rad)
    B.run_house(wd, "m", mat)
    B.run_force(wd, "m", [(c, 1, 1.0, 0.0)], fs, copy_to="FILE9001")
    B.run_force(wd, "m", [(c, 5, 1.0, 0.0)], fs, copy_to="FILE9002")
    B.run_analys(wd, "m", fs, type=1, simul=2)
    f81, f82 = B.read_file8(wd, "FILE8001"), B.read_file8(wd, "FILE8002")
    Cc = np.stack([np.stack([B.tf(f81, c, 1), B.tf(f82, c, 1)], -1), np.stack([B.tf(f81, c, 5), B.tf(f82, c, 5)], -1)],
                  1)
    S = np.linalg.inv(Cc)                                 # (nF, 2, 2): [x, theta_y] at the centre node
    # 2. SDOF on the mat, seismic
    I = k * h ** 3 / (3 * p["E"])
    stick = B.Stick(heights=[h], masses=[m], E=p["E"], A=10.0, I=I, beta=beta, J=2 * I)
    mdl = B.stick_on_mat(site, stick, half_width=R, ndiv=p["nd"], shape="disk", E_rigid=E_rigid)
    r.require("the SDOF model reuses the mat node numbering", mdl["centre"] == c)
    B.run_house(wd, "m", mdl)
    KG = analys_global_impedance(wd, "m", fs)              # seismic run with the global impedance
    f8 = B.read_file8(wd)
    purge_work_files(wd)                                   # all results are in memory now
    top = mdl["top"]
    Hm, Hc, Ht = B.tf(f8, top, 1), B.tf(f8, c, 1), B.tf(f8, c, 5)
    F2 = np.array([three_dof(2 * math.pi * f, m, h, kstar, S[q]) for q, f in enumerate(freq)])
    H2m = 1 + F2[:, 1] + h * F2[:, 2] + F2[:, 0]
    r.check("3-DOF identity: max |H_mass(SASSI) - H_mass(F.2 with the code's impedance)| / |H|",
            float(np.max(np.abs(Hm - H2m) / np.abs(H2m))), 0.0, atol=1e-8)
    r.check("3-DOF identity: foundation translation 1 + u0", float(np.max(np.abs(Hc - (1 + F2[:, 1])) / np.abs(Hc))), 0.0,
            atol=1e-8)
    r.check("3-DOF identity: foundation rocking theta (x h, relative to |H_mass|)",
            float(np.max(np.abs(Ht - F2[:, 2]) * h / np.abs(H2m))), 0.0, atol=1e-8)
    sel = np.ix_([0, 4], [0, 4])
    dS = max(float(np.abs(S[q] - KG[q][sel]).max() / np.abs(KG[q][sel]).max()) for q in range(len(freq)))
    F2g = np.array([_total_motion(2 * math.pi * f, m, h, kstar, KG[q][sel]) for q, f in enumerate(freq)])
    r.notes.append(f"centre-node impedance S (FORCE runs) vs the rigid-body K_G of ANALYS <impe> = 2: max relative "
                   f"difference {dS:.2e} (flexibility of the stiff massless BEAMS links, E = {p['rigid']:.0e} G); F.2 "
                   f"fed with K_G instead of S: max |dH|/|H| = {float(np.max(np.abs(F2g - Hm) / np.abs(Hm))):.2e}")
    # 3. worked example
    Kx, Kt, Kxt = KG[0][0, 0].real, KG[0][4, 4].real, KG[0][0, 4].real
    r.check("static K_x (K_G at 0.05 Hz) vs 8GR/(2-nu) = 9.60e5 kN/m", Kx, VP17_REF["Kx"], rtol=0.05)
    r.check("static K_theta (K_G at 0.05 Hz) vs 8GR^3/(3(1-nu)) = 8.00e7 kN m", Kt, VP17_REF["Kt"], rtol=0.05)
    ratio = math.sqrt(1 + k / Kx + k * h * h / Kt)
    r.check("Veletsos-Meek T~/T with the code's static springs vs 1.15805", ratio, VP17_REF["ratio"], rtol=0.05)
    band = slice(1, None)
    fpk, amp = _resonance(freq[band], Hm[band], VP17_REF["f_res"])
    r.check("resonance of |u_t/u_g| at the mass vs 1.717 Hz (F.2 with Veletsos-Verbic)", fpk, VP17_REF["f_res"],
            rtol=0.05)
    # D-W2-04: the R2 F.1 peak (F.2 + Veletsos-Verbic impedances) is superseded by the guard against F.2
    # fed with the rigorous welded impedance (below)
    r.inform("peak |u_t/u_g| at the mass vs 6.71 (F.2 with Veletsos-Verbic)", amp, VP17_REF["peak"],
             note="informative: Veletsos-Verbic impedances underestimate rocking damping and omit sliding-rocking coupling")
    # the reference re-derived: F.2 with Veletsos-Verbic impedances (R2 F.1)
    ff = np.arange(1.5, 1.95, 0.0005)
    hv = np.array([_total_motion(2 * math.pi * f, m, h, kstar, _vv_impedance(2 * math.pi * f * R / VS, R)) for f in ff])
    fv, av = _resonance(ff, hv, VP17_REF["f_res"])
    r.check("R2 F.1 reference re-derived: F.2 + Veletsos-Verbic resonance (Hz)", fv, 1.717, atol=5e-4)
    r.check("R2 F.1 reference re-derived: F.2 + Veletsos-Verbic peak", av, 6.71, atol=5e-3)
    # evidence: relaxed contact, rigorous impedances
    xy = np.array([mat.house.xyz(n)[:2] for n in mat["mat"]])
    T = rigid_transform(xy)
    KR, KD = [], []
    for q in range(len(freq)):                             # one F_ff at a time (11 MB each)
        F = flexibility_matrix(files["FILE3"], q, xy, np.ones(len(xy), int))
        KR.append(relaxed_impedance(F, T))
        KD.append(independent_global_impedance(F, T))
    KR = np.array(KR)
    kg_checks(r, "seismic run with the SDOF", KG, np.array(KD))
    hr = np.array([_total_motion(2 * math.pi * f, m, h, kstar, KR[q][sel]) for q, f in enumerate(freq)])
    hb = np.array([_total_motion(2 * math.pi * f, m, h, kstar, _bem_impedance(2 * math.pi * f * R / VS, R)) for f in ff])
    hn = np.array([_total_motion(2 * math.pi * f, m, h, kstar, np.diag(np.diag(S[q]))) for q, f in enumerate(freq)])
    frr, arr = _resonance(freq[band], hr[band], VP17_REF["f_res"])
    fbb, abb = _resonance(ff, hb, VP17_REF["f_res"])
    fnn, ann = _resonance(freq[band], hn[band], VP17_REF["f_res"])
    r.check("F.2 with the code's relaxed-contact impedance (from F_ff) vs F.2 with the rigorous relaxed BEM impedance: "
            "resonance (derived reference)", frr, fbb, rtol=0.05)
    r.check("F.2 with the code's relaxed-contact impedance vs F.2 with the rigorous relaxed BEM impedance: peak "
            "(derived reference)", arr, abb, rtol=0.05)
    # guard of the failing peak: F.2 fed with the rigorous *welded* impedance (exact-kernel BEM, derived)
    Kw = welded_disk_impedance(VP17_BEM_A0)[:, [0, 4]][:, :, [0, 4]]
    spl = CubicSpline(np.asarray(VP17_BEM_A0), Kw, axis=0)
    unit = GMOD * np.array([[R, R ** 2], [R ** 2, R ** 3]])
    hw = np.array([_total_motion(2 * math.pi * f, m, h, kstar, spl(2 * math.pi * f * R / VS) * unit) for f in ff])
    fww, aww = _resonance(ff, hw, VP17_REF["f_res"])
    r.check("guard (peak): peak |u_t/u_g| at the mass (SASSI, welded mat) vs F.2 fed with the welded exact-kernel "
            "BEM impedance (derived reference)", amp, aww, rtol=0.05)
    r.check("guard (peak): resonance of |u_t/u_g| (SASSI, welded mat) vs F.2 fed with the welded exact-kernel BEM "
            "impedance (derived reference)", fpk, fww, rtol=0.05)
    Kw0 = spl(1.07)                                         # [x, theta_y] block, G = R = 1
    r.notes.append(
        f"k = {k:.1f} kN/m; code statics (welded) K_x = {Kx:.4e} kN/m, K_theta = {Kt:.4e} kN m, K_x-theta = {Kxt:.4e} kN "
        f"({Kxt / (Kx * R):+.3f} K_x R); T~/T = {ratio:.5f}.")
    r.notes.append(
        f"Peak of |u_t/u_g| (resonance): SASSI welded mat {amp:.3f} ({fpk:.3f} Hz); F.2 with the rigorous welded BEM "
        f"impedance {aww:.3f} ({fww:.3f} Hz; its dynamic coupling at a0 = 1.07: K_x-theta = "
        f"{Kw0[0, 1].real:+.3f}{Kw0[0, 1].imag:+.3f}i G R^2, {Kw0[0, 1].real / Kw0[0, 0].real:+.3f} Re K_x R); "
        f"same SASSI impedance without the "
        f"sliding-rocking coupling {ann:.3f} ({fnn:.3f} Hz); code's relaxed-contact impedance {arr:.3f} ({frr:.3f} Hz); "
        f"rigorous relaxed BEM impedance {abb:.3f} ({fbb:.3f} Hz); Veletsos-Verbic {av:.3f} ({fv:.4f} Hz).  The R2 value "
        "6.71 inherits the low Veletsos-Verbic rocking damping (VP-11); the welded mat adds the sliding-rocking "
        "coupling (K_x-theta ~ -0.08 K_x R, confirmed by the welded static BEM), which lowers the peak further.")
    return r
