"""Verification problems of 2D (plane-strain) SSI and of SYMM symmetry / antisymmetry planes
(requirements 4.3 POINT2, 4.6 item 2, D-PNT-01, D-ANL-12, 6.3; R1 2.6, 3.5, 4.1; R2 B.4).

* VP-42  rigid strip on a layer over a rigid base (2D, per unit length): the static horizontal and
         rocking stiffnesses through SITE -> POINT2 -> HOUSE (rigid BEAMS links) -> FORCE -> ANALYS
         (vibration, three load cases) vs Jakub & Roesset (1977): ``K_x/G = 1.175 (1 + 2.15 B/H)``,
         ``K_phi/(G B^2) = 2.394 (1 + 0.17 B/H)`` for nu = 0.30 and 1/8 <= B/H <= 1/2 (R2 B.4, 5 %),
         extrapolated to a vanishing mesh size; the chain equals the direct rigid-strip impedance
         ``T^T X_ff T`` of FILE3; ``K_x -> 0`` as H -> infinity with the slope of the 2D Cerruti
         line-load compliance ``d(1/K_x)/d ln H = (1 - nu)/(pi G)``
* VP-T1  POINT2 far field vs Kausel's exact line-load Green functions (R1 V7 configuration: f = 4 Hz,
         R0 = 1 m, load at interface 4, observation at the surface; ``|x| >= 2 R0``, 0.3 %)
* VP-T2  SYMM: symmetric 3D models -- a flexible surface mat with a stick, and an embedded basement
         (FV, excavated SOLIDs, SHELL walls and slabs) with a stick -- analysed as full, half and
         quarter models give the same transfer functions to 1e-8, for X input (antisymmetry about
         x = 0, symmetry about y = 0) and Z input (symmetry about both); the DOFs the symmetry
         constrains vanish in the full model
* VP-T3  2D zero-SSI identity: an embedded PLANE excavation whose structure equals the excavated soil
         gives U = U'_f at every node (vertical SV and P, and an inclined SV wave), 1e-8

Units: m, t, kN, s (g = 9.81 m/s2).
"""
from __future__ import annotations

import math
import shutil
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from ...core import greens_tlm, ssi2d, strip2d, tlm
from ...core import symmetry as SYM
from ...core.flexibility import flexibility_block, flexibility_matrix, frequency_row
from ...core.freefield import free_field_at_nodes
from ...io.container import read_container
from .. import VPResult, problem, worse
from .. import builders as B


def _fresh(path: Path) -> Path:
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    return path


# ======================================================================================
# VP-42  rigid strip on a layer over a rigid base (2D statics)
# ======================================================================================
#: soil of VP-42 (nu = 0.30 as in the R2 B.4 row, small hysteretic damping: |c(beta)| = 1, so |K| is the
#: elastic stiffness in the static limit -- correspondence principle with one damping ratio)
VP42_NU, VP42_VS, VP42_RHO, VP42_BETA = 0.30, 100.0, 2.0, 0.02
#: strip half-width and the static-limit frequency a0 = w B / Vs
VP42_B, VP42_A0 = 1.0, 0.01
#: meshes (elements across the strip); the stiffness converges as O(h) (edge singularity of the contact
#: stress), so the h -> 0 value is the Richardson extrapolation of the two finest meshes
VP42_NDIV = (8, 16, 32)
#: sublayer grading: thickness h = 2B/ndiv at the surface growing by this factor with depth
VP42_GROWTH = 1.15


def graded_layers(H: float, h0: float, growth: float = VP42_GROWTH) -> List[float]:
    """Sublayer thicknesses from h0 at the surface, growing geometrically, summing to H."""
    out: List[float] = []
    z, h = 0.0, h0
    while z < H - 1e-9 * H:
        hh = min(h, H - z)
        if H - z - hh < 0.5 * hh:
            hh = H - z
        out.append(hh)
        z += hh
        h *= growth
    return out


def strip_site(H: float, h0: float) -> B.Site:
    vp = VP42_VS * math.sqrt(2.0 * (1.0 - VP42_NU) / (1.0 - 2.0 * VP42_NU))
    rows = [(t, VP42_VS, vp, VP42_RHO, VP42_BETA) for t in graded_layers(H, h0)]
    return B.layered_site(rows, (VP42_VS, vp, VP42_RHO, VP42_BETA), nl=0)       # rigid base


def strip_stiffness(wd: Path, H: float, ndiv: int) -> Dict[str, object]:
    """Rigid strip (half-width VP42_B, ``ndiv`` elements) on a layer of thickness H over a rigid base:
    FORCE unit loads Fx, Fz, My at the centre -> ANALYS vibration (FILE9001..003) -> compliance C (3x3,
    X Z YY) -> K = C^-1; and the direct impedance ``T^T X_ff T`` of the same interaction nodes from
    FILE3 (2D in-plane block).  Returns the stiffness matrices and the soil G."""
    _fresh(wd)
    h = 2.0 * VP42_B / ndiv
    site = strip_site(H, h)
    G = VP42_RHO * VP42_VS ** 2
    df = VP42_A0 * VP42_VS / VP42_B / (2.0 * math.pi)
    fs = B.FrequencySet.harmonic(df, [1])
    mdl = B.rigid_strip_2d(site, half_width=VP42_B, ndiv=ndiv)
    B.run_soil(wd, "m", site, fs, layer=0, rad=mdl.rad, mode2=False, dim=1)
    f4 = B.run_house(wd, "m", mdl)
    c = mdl["centre"]
    for k, dof in enumerate((1, 3, 5), start=1):
        B.run_force(wd, "m", [(c, dof, 1.0, 0.0)], fs, copy_to=f"FILE9{k:03d}")
    B.run_analys(wd, "m", fs, type=1, simul=3)
    f8 = [B.read_file8(wd, f"FILE8{k:03d}") for k in (1, 2, 3)]
    C = np.array([[B.tf(f8[j], c, d)[0] for j in range(3)] for d in (1, 3, 5)])
    K_chain = np.linalg.inv(C)
    # independent path: X_ff of the interaction nodes from FILE3 and the 2D rigid-body transformation
    f3 = read_container(wd / "FILE3", "FILE3")
    xyz = np.asarray(f4["x_int_xyz"], float)
    F = flexibility_matrix(f3, 0, xyz[:, :2], np.asarray(f4["int_iface"]))
    n = xyz.shape[0]
    ip = np.array([[3 * i, 3 * i + 2] for i in range(n)]).ravel()
    X = np.zeros_like(F)
    X[np.ix_(ip, ip)] = np.linalg.inv(F[np.ix_(ip, ip)])
    xc = np.asarray(f4["node_xyz"], float)[list(np.asarray(f4["node_id"])).index(c)]
    K_direct = ssi2d.rigid_strip_impedance(X, xyz, (xc[0], xc[2]))
    shutil.rmtree(wd, ignore_errors=True)
    return {"K": K_chain, "K_direct": K_direct, "G": G}


def bonded_strip_rocking(nu: float, n: int = 800) -> float:
    """Rocking stiffness ``K_phi/(G B^2)`` of a rigid strip *bonded* (welded) to an elastic half-plane
    (plane strain): own boundary-element solution with piecewise-constant shear and normal tractions on
    n panels and the exact surface kernels of the Flamant / Cerruti line loads (``(1 - nu)/(pi G) ln r``
    and the coupling ``(1 - 2 nu)/(4 G) sgn x``); u_x = 0, u_z = theta x.  (The relaxed punch has the exact
    ``pi/(2 (1 - nu))``.)  Informative reference of VP-42 for H -> infinity."""
    G, Bw = 1.0, 1.0
    c = (1.0 - nu) / (math.pi * G)
    a = (1.0 - 2.0 * nu) / (4.0 * G)
    e = np.linspace(-Bw, Bw, n + 1)
    xm = 0.5 * (e[:-1] + e[1:])
    h = e[1] - e[0]

    def iln(lo, hi):
        f = lambda s: np.where(s == 0, 0.0, s * np.log(np.abs(np.where(s == 0, 1.0, s))) - s)
        return f(hi) - f(lo)

    lo, hi = xm[:, None] - e[None, 1:], xm[:, None] - e[None, :-1]
    L = -c * iln(lo, hi)
    S = a * (np.abs(hi) - np.abs(lo))
    N = 2 * n
    M = np.zeros((N + 2, N + 2))
    M[:N, :N] = np.block([[L, S], [-S, L]])
    M[:n, N] = -1.0
    M[n:N, N + 1] = -1.0
    M[N, :n] = h
    M[N + 1, n:N] = h
    rhs = np.zeros(N + 2)
    rhs[n:N] = xm
    sol = np.linalg.solve(M, rhs)
    return float(np.sum(sol[n:N] * xm * h))


@problem("VP-42", "2D rigid strip on a layer over a rigid base: static K_x and K_phi (POINT2, ANALYS 2D)", tier="P1",
         modules=["SITE", "POINT", "HOUSE", "FORCE", "ANALYS"], source="R2 B.4 (Jakub & Roesset 1977, MIT R77-36)")
def vp42(workdir: Path) -> VPResult:
    """Rigid massless strip of half-width B = 1 m (``ndiv`` interaction nodes linked to the centre node by
    stiff BEAMS; out-of-plane DOFs fixed) on the surface of a uniform layer (nu = 0.30) of thickness H on
    a rigid base, at a0 = w B / Vs = 0.01.  The layer is divided into sublayers of thickness h = 2B/ndiv
    at the surface growing by 15 % with depth; POINT2 central zone R0 = h (R1 3.4).  The stiffness matrix
    is the inverse of the compliance of three FORCE load cases (Fx, Fz, My at the centre) through ANALYS.
    The FE stiffness converges as O(h) (singular edge stresses of the rigid punch): the criterion uses the
    Richardson value 2 K(32) - K(16) with the observed ratio of the mesh differences checked (~2)."""
    r = VPResult()
    wd = Path(workdir)
    res: Dict[Tuple[float, int], Dict[str, object]] = {}
    ratios = (0.5, 0.25, 0.125, 1.0 / 16.0, 1.0 / 32.0)
    for bh in ratios:
        H = VP42_B / bh
        for ndiv in (VP42_NDIV if bh >= 0.125 else VP42_NDIV[1:]):
            res[(bh, ndiv)] = strip_stiffness(wd / f"s{ndiv}", H, ndiv)
    G = float(res[(0.5, 32)]["G"])

    def kx(bh, nd):
        return abs(res[(bh, nd)]["K"][0, 0]) / G

    def kp(bh, nd):
        return abs(res[(bh, nd)]["K"][2, 2]) / (G * VP42_B ** 2)

    worst = 0.0
    for key, v in res.items():
        K, Kd = v["K"], v["K_direct"]
        for i in (0, 1, 2):
            worst = worse(worst, abs(K[i, i] - Kd[i, i]) / abs(Kd[i, i]))
    r.check("FORCE + ANALYS stiffness = direct rigid-strip impedance T^T X_ff T (diagonal terms, all meshes)",
            worst, 0.0, atol=1e-4, note="the stiff links (E = 1e5 G) leave an O(1e-5) difference")
    ext: Dict[float, Tuple[float, float]] = {}
    for bh in ratios:
        a16, a32 = kx(bh, 16), kx(bh, 32)
        b16, b32 = kp(bh, 16), kp(bh, 32)
        ext[bh] = (2.0 * a32 - a16, 2.0 * b32 - b16)
        if bh >= 0.125:
            a8, b8 = kx(bh, 8), kp(bh, 8)
            for name, q8, q16, q32 in (("K_x", a8, a16, a32), ("K_phi", b8, b16, b32)):
                r.check(f"B/H = {bh:g}: O(h) convergence of {name}: (K8 - K16)/(K16 - K32)", (q8 - q16) / (q16 - q32), 2.0,
                        rtol=0.15, note="first-order convergence justifies the Richardson extrapolation")
            rx, rp = ssi2d.strip_stiffness_reference(VP42_B, VP42_B / bh, 1.0, VP42_NU)
            r.check(f"B/H = {bh:g}: K_x/G (h -> 0) vs 1.175 (1 + 2.15 B/H)", ext[bh][0], rx, rtol=0.05)
            r.check(f"B/H = {bh:g}: K_phi/(G B^2) (h -> 0) vs 2.394 (1 + 0.17 B/H)", ext[bh][1], rp, rtol=0.05)
            r.inform(f"B/H = {bh:g}: K_x/G finest mesh (32 elements)", a32, rx)
            r.inform(f"B/H = {bh:g}: K_phi/(G B^2) finest mesh (32 elements)", b32, rp)
    # K_x -> 0 as H -> infinity: the compliance grows like the 2D Cerruti line-load displacement
    kxs = [ext[bh][0] for bh in ratios]
    r.require("K_x decreases monotonically as H grows (B/H = 1/2 ... 1/32): K_x -> 0 for a half-space",
              all(b < a for a, b in zip(kxs, kxs[1:])))
    slope = (1.0 / ext[1.0 / 32.0][0] - 1.0 / ext[1.0 / 16.0][0]) / math.log(2.0)
    r.check("d(G/K_x)/d ln H between B/H = 1/16 and 1/32 vs the half-plane line-load value (1 - nu)/pi",
            slope, (1.0 - VP42_NU) / math.pi, rtol=0.05,
            note="u_x = (1 - nu)/(pi G) Q ln(1/r) (Cerruti, plane strain): the static horizontal compliance diverges "
                 "logarithmically with the layer depth")
    kb = bonded_strip_rocking(VP42_NU)
    r.inform("B/H = 1/32: K_phi/(G B^2) vs bonded half-plane BEM (own Flamant/Cerruti panel solution)", ext[1.0 / 32.0][1], kb)
    r.inform("B/H = 1/32: K_phi/(G B^2) vs relaxed (smooth) half-plane punch pi/(2(1 - nu))", ext[1.0 / 32.0][1],
             math.pi / (2.0 * (1.0 - VP42_NU)))
    r.notes.append("Extrapolated (h -> 0) K_x/G, K_phi/(G B^2): " + "; ".join(
        f"B/H = {bh:g}: {ext[bh][0]:.4f}, {ext[bh][1]:.4f}" for bh in ratios))
    r.notes.append("Jakub & Roesset's fits are FE results with consistent boundaries (R2 B.4, approximate); the K_x fit "
                   "tends to a constant for H -> infinity while the exact static K_x of a strip on a half-plane is 0 "
                   "(K_x/G drops below the fit for B/H < 1/8: informative rows above).  The rocking stiffness of the "
                   "welded strip tends to the bonded half-plane value (BEM), above the smooth-punch pi/(2(1 - nu)).")
    return r


# ======================================================================================
# VP-T1  POINT2 far field vs exact line-load Green functions
# ======================================================================================
#: R1 V7 (docs/spec/R1_checks/point2.py): 6 x 1 m (Vs 150, Vp 300, rho 1.9) + 6 x 2 m (Vs 300, Vp 600,
#: rho 2.0), 5 % damping, rigid base, f = 4 Hz, R0 = 1 m, unit line loads at interface 4 (depth 3 m)
T1_LAYERS = [(1.0, 150.0, 300.0, 1.9, 0.05)] * 6 + [(2.0, 300.0, 600.0, 2.0, 0.05)] * 6
T1_HS = (300.0, 600.0, 2.0, 0.05)
T1_R0 = 1.0
T1_F = 4.0
T1_X = (2.0, 3.3, 5.0, 12.0)                     # distances in R0 (R1 V7: 2, 5, 12)
T1_COMPONENTS = ("xx", "zx", "xz", "zz", "yy")
#: FE + transmitting-boundary values printed by docs/spec/R1_checks/point2.py (consistent mass) at x = 2, 5, 12 m,
#: surface observation, unit line loads at interface 4: u_x and u_z for the x load (xx, zx) and the z load (xz, zz)
R1V7_PRINTED = {"xx": (5.677e-06 - 1.2108e-05j, 4.952e-06 - 1.1952e-05j, 1.817e-06 - 1.1192e-05j),
                "zx": (1.478e-06 - 2.1e-07j, 1.284e-06 - 2.61e-07j, 6.97e-07 - 3.32e-07j),
                "xz": (1.586e-06 - 1.58e-07j, 1.159e-06 - 8.5e-08j, 8.8e-08 + 1.27e-07j),
                "zz": (4.681e-06 - 6.61e-07j, 1.706e-06 - 3.06e-07j, 1.9e-07 - 7.8e-08j)}


def point2_components(f3, iface_obs: int, iface_load: int, xs: Sequence[float]) -> Dict[str, np.ndarray]:
    """POINT2 values at abscissae ``xs`` (observation on ``iface_obs``) for unit line loads at x = 0 on
    ``iface_load``, from FILE3 through :func:`sassi.core.flexibility.flexibility_block`."""
    xs = np.asarray(xs, float)
    F = flexibility_block(f3, 0, np.c_[xs, np.zeros_like(xs)], [iface_obs] * xs.size, np.zeros((1, 2)), [iface_load])
    F = F.reshape(xs.size, 3, 1, 3)[:, :, 0, :]
    return {"xx": F[:, 0, 0], "zx": F[:, 2, 0], "xz": F[:, 0, 2], "zz": F[:, 2, 2], "yy": F[:, 1, 1]}


@problem("VP-T1", "POINT2 far field vs Kausel's exact line-load Green functions (R1 V7)", tier="P1",
         modules=["SITE", "POINT"], source="R1 2.6, 3.5, V7; Kausel (1981) Eq. 47-48")
def vpt1(workdir: Path) -> VPResult:
    """SITE Mode 1 and POINT (``<dim>`` = 1, POINT2: a two-element plane-strain strip |x| <= R0 with
    Waas-Lysmer boundaries) for the R1 V7 column; the FILE3 far field ``sum alpha_j phi_j exp(-i k_j
    (|x| - R0))`` is compared with the exact discrete-layer line-load Green functions of the same column
    (Kausel Eq. 47-48, :func:`sassi.core.greens_tlm.line_load`, an independent modal-sum implementation)."""
    r = VPResult()
    wd = _fresh(Path(workdir) / "t1")
    site = B.layered_site(T1_LAYERS, T1_HS, nl=0)
    fs = B.FrequencySet.harmonic(T1_F, [1])
    B.run_soil(wd, "m", site, fs, layer=3, rad=T1_R0, mode2=False, dim=1)
    f2 = read_container(wd / "FILE2", "FILE2")
    f3 = read_container(wd / "FILE3", "FILE3")
    r.require("FILE3 written by POINT2 (meta dim = 1, program POINT2; D-PNT-01)",
              int(f3.meta["dim"]) == 1 and str(f3.meta.get("program")) == "POINT2")
    md = tlm.modes_from_file2(f2, 0)
    iu = np.asarray(f2["iface_user"], int)
    xs = np.asarray(T1_X) * T1_R0
    got = point2_components(f3, 1, 4, xs)
    for k, x in enumerate(xs):
        g = greens_tlm.line_load(md, iu[0], iu[3], x)
        # D-W3-12: 0.33 % at 2 R0 (R1 V7's prototype prints 0.32 % there, rounded by R1 to 0.3 %), 0.3 % beyond
        tol = 0.0033 if abs(x / T1_R0 - 2.0) < 1e-9 else 0.003
        for c in T1_COMPONENTS:
            r.check(f"R1 V7: load at interface 4, surface, |x| = {x / T1_R0:g} R0: {c}", abs(got[c][k] - g[c]) / abs(g[c]),
                    0.0, atol=tol, note="relative difference")
    # all load / observation interface pairs at large distance; parity of the coupling terms
    worst = 0.0
    for lo in (1, 2, 3, 4):
        for ob in (1, 2, 3, 4):
            v = point2_components(f3, ob, lo, [12.0 * T1_R0, 20.0 * T1_R0])
            for k, x in enumerate((12.0 * T1_R0, 20.0 * T1_R0)):
                g = greens_tlm.line_load(md, iu[ob - 1], iu[lo - 1], x)
                for c in T1_COMPONENTS:
                    worst = worse(worst, abs(v[c][k] - g[c]) / abs(g[c]))
    r.check("all load/observation interfaces 1..4, |x| = 12 and 20 R0: max relative difference", worst, 0.0, atol=0.003)
    a = point2_components(f3, 1, 4, [5.0 * T1_R0])
    b = point2_components(f3, 1, 4, [-5.0 * T1_R0])
    odd = max(abs(a[c][0] + b[c][0]) / abs(a[c][0]) for c in ("zx", "xz"))
    even = max(abs(a[c][0] - b[c][0]) / abs(a[c][0]) for c in ("xx", "zz", "yy"))
    r.check("parity: u_z | P_x and u_x | P_z odd in x (|F(x) + F(-x)| / |F(x)|, x = 5 R0)", odd, 0.0, atol=1e-12)
    r.check("parity: u_x | P_x, u_z | P_z, u_y | P_y even in x", even, 0.0, atol=1e-12)
    for lo, ob in ((1, 1), (4, 4)):
        v = point2_components(f3, ob, lo, [2.0 * T1_R0])
        g = greens_tlm.line_load(md, iu[ob - 1], iu[lo - 1], 2.0 * T1_R0)
        for c in T1_COMPONENTS:
            r.inform(f"load at interface {lo}, observation at {ob}, |x| = 2 R0: {c} (POINT2 vs exact)", abs(v[c][0]), abs(g[c]))
    # the POINT2 core is the R1 V7 prototype: with the prototype's consistent mass it reproduces the values that
    # docs/spec/R1_checks/point2.py prints (9 decimals), so the 2 R0 accuracy is the prototype's own
    worst = 0.0
    col = tlm.column_from_file2(f2, 0)
    mdc = tlm.column_modes(col, 2.0 * math.pi * T1_F, tlm.MASS_CONSISTENT)
    sol = strip2d.solve_point2(col, mdc, 2.0 * math.pi * T1_F, T1_R0, [iu[3]], tlm.MASS_CONSISTENT)
    ext = strip2d.exterior_2d(mdc.kR, mdc.kL, mdc.phix[[iu[0]]], mdc.phiz[[iu[0]]], mdc.phiy[[iu[0]]], sol.alpha_x,
                              sol.alpha_z, sol.alpha_y, [2.0, 5.0, 12.0], T1_R0)
    for k in range(3):
        for c in ("xx", "zx", "xz", "zz"):
            worst = worse(worst, abs(complex(ext[c][k, 0, 0]) - R1V7_PRINTED[c][k]))
    r.check("POINT2 core with the prototype's consistent mass = R1 V7 printed FE+TB values (x = 2, 5, 12 m; 9 decimals)",
            worst, 0.0, atol=1e-9, note="max |difference| (m per unit load); rounding of the printout <= 7e-10")
    r.notes.append("R1 V7 (docs/spec/R1_checks/point2.py) prints, at |x| = 2 R0, u_z due to the vertical load "
                   "4.681e-06 - 6.61e-07j vs Eq. 47 4.666e-06 - 6.59e-07j, i.e. 0.32 %, summarised in R1 as "
                   "'within 0.3 %'; SASSI-EDU's POINT2 reproduces the prototype (0.325 %).")
    shutil.rmtree(wd, ignore_errors=True)
    return r


# ======================================================================================
# VP-T2  SYMM half and quarter models vs the full model
# ======================================================================================
T2_SITE_LAYERS = [(2.0, 150.0, 300.0, 1.9, 0.05), (3.0, 250.0, 500.0, 2.0, 0.04), (5.0, 350.0, 700.0, 2.1, 0.03)]
T2_SITE_HS = (600.0, 1200.0, 2.2, 0.02)
T2_FS = (4, 10, 20, 41)                     # frequency numbers on the 0.01 s x 1024 grid
#: (input, SITE wave, SYMM type of the plane x = 0, of the plane y = 0)
T2_INPUTS = (("X", "SV", SYM.ANTISYMMETRIC, SYM.SYMMETRIC), ("Z", "P", SYM.SYMMETRIC, SYM.SYMMETRIC))
PART_FACTOR = {"full": 1.0, "half": 0.5, "quarter": 0.25}


def symmetric_model(site: B.Site, part: str, types: Tuple[int, int], embedded: bool, a: float = 4.0,
                    ndiv: int = 4) -> B.Model:
    """A doubly symmetric model about x = 0 and y = 0 (``part`` 'full', 'half' = x >= 0 with SYMM plane 1
    at x = 0, 'quarter' = x, y >= 0 with planes 1 (x = 0) and 2 (y = 0); ``types`` = SYMM types of the
    two planes):

    * surface (``embedded`` False): a flexible concrete SHELL mat 2a x 2a (t = 0.8 m), every mat node an
      interaction node, drilling rotations fixed;
    * embedded: a 2a x 2a box through the two upper layers (FV: excavated SOLIDs, every excavation node
      interacts) with a concrete SHELL basement (walls, base and roof slabs, t = 0.4 m);

    plus a lumped-mass BEAMS stick on the centre of the top surface (on both planes).  The stick section,
    its masses and rotary inertias are multiplied by 1/2 (half) or 1/4 (quarter): elements lying in a
    symmetry plane carry that share of the full model (manual 9.3.2); the stick's orientation nodes are
    off the mesh."""
    f = PART_FACTOR[part]
    hb = B.HouseBuilder(site, title=f"symmetric {'embedded basement' if embedded else 'surface mat'} + stick ({part})")
    g = np.linspace(-a, a, ndiv + 1)
    xs = g[g >= -1e-9] if part != "full" else g
    ys = g[g >= -1e-9] if part == "quarter" else g
    n_emb = 2 if embedded else 0
    z = site.elevations()[:n_emb + 1]
    ids: Dict[Tuple[int, int, int], int] = {}
    for lev in range(n_emb, -1, -1):                       # bottom-up numbering
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                ids[(i, j, lev)] = hb.node(x, y, z[lev])
    nx, ny = len(xs) - 1, len(ys) - 1
    if embedded:
        gexc = hb.group(1, "excavated soil")
        for lev in range(n_emb):
            for j in range(ny):
                for i in range(nx):
                    bot = [ids[(i, j, lev + 1)], ids[(i + 1, j, lev + 1)], ids[(i + 1, j + 1, lev + 1)], ids[(i, j + 1, lev + 1)]]
                    top = [ids[(i, j, lev)], ids[(i + 1, j, lev)], ids[(i + 1, j + 1, lev)], ids[(i, j + 1, lev)]]
                    hb.solid(bot + top, mat=lev + 1, etype=2, group=gexc)
    hb.set_interaction(list(ids.values()))
    conc = hb.elastic(3.0e7, 0.2, 2.4, 0.05)
    gs = hb.group(3, "basement" if embedded else "mat")
    t = 0.4 if embedded else 0.8
    for lev in ((0, n_emb) if embedded else (0,)):
        for j in range(ny):
            for i in range(nx):
                hb.shell([ids[(i, j, lev)], ids[(i + 1, j, lev)], ids[(i + 1, j + 1, lev)], ids[(i, j + 1, lev)]], conc, t, group=gs)
    if embedded:
        for lev in range(n_emb):
            for k in range(ny):
                for i in range(nx + 1):
                    if abs(abs(xs[i]) - a) < 1e-9:
                        hb.shell([ids[(i, k, lev + 1)], ids[(i, k + 1, lev + 1)], ids[(i, k + 1, lev)], ids[(i, k, lev)]],
                                 conc, t, group=gs)
            for k in range(nx):
                for j in range(ny + 1):
                    if abs(abs(ys[j]) - a) < 1e-9:
                        hb.shell([ids[(k, j, lev + 1)], ids[(k + 1, j, lev + 1)], ids[(k + 1, j, lev)], ids[(k, j, lev)]],
                                 conc, t, group=gs)
    for (i, j, lev), n in ids.items():                     # rotations without stiffness (EDU-06)
        faces = [d for d, on in ((4, abs(abs(xs[i]) - a) < 1e-9), (5, abs(abs(ys[j]) - a) < 1e-9),
                                 (6, lev in (0, n_emb))) if on]
        if not embedded:
            hb.fix(n, [6])                                 # flat mat: drilling rotation
        elif len(faces) == 1:
            hb.fix(n, faces)                               # one flat face: rotation about its normal
        elif not faces:
            hb.fix(n, [4, 5, 6])                           # interior excavation node: SOLID only
    c = hb.find(0.0, 0.0, z[0])
    mat = hb.material(1, 3.0e7, 0.2, 0.0, 0.05, 0.05)
    sec = hb.section(2.0 * f, 0.0, 0.0, 2.4 * f, 1.2 * f, 1.2 * f)
    gst = hb.group(2, "stick")
    prev = c
    nodes = []
    for k, (hgt, m, rot) in enumerate(((3.0, 60.0, 20.0), (6.0, 50.0, 20.0), (9.0, 30.0, 10.0))):
        n = hb.node(0.0, 0.0, z[0] + hgt)
        kn = hb.orientation_node(0.37 * a, 0.0, z[0] + hgt - 1.5)     # off the mesh, beside the stick
        hb.beam(prev, n, kn, mat, sec, group=gst)
        hb.mass(n, m * f, m * f, m * f, rot * f, rot * f, 0.0)
        nodes.append(n)
        prev = n
    if part in ("half", "quarter"):
        hb.symmetry(1, types[0], [hb.find(0.0, 0.0, z[0]), hb.find(0.0, a, z[0]), nodes[-1]])
    if part == "quarter":
        hb.symmetry(2, types[1], [hb.find(0.0, 0.0, z[0]), hb.find(a, 0.0, z[0]), nodes[-1]])
    return B.Model(hb, {"centre": c, "stick": nodes, "top": nodes[-1]}, rad=0.9 * 2.0 * a / ndiv, layer=n_emb,
                   description=hb.title)


def compare_reduced(full: Tuple[object, Dict[int, np.ndarray]], red: Tuple[object, Dict[int, np.ndarray]]
                    ) -> Tuple[float, int]:
    """max |H_reduced - H_full| / max |H_full| over every equation of the reduced FILE8 (nodes matched
    by coordinates among the nodes with equations), and the number of DOFs compared."""
    f8f, xyzf = full
    f8r, xyzr = red
    key = lambda p: tuple(np.round(np.asarray(p, float), 6))
    active = set(int(v) for v in np.asarray(f8f["eq_node"]))
    by_xyz = {key(p): n for n, p in xyzf.items() if n in active}
    scale = float(np.abs(f8f["H"]).max())
    err = 0.0
    for k, (n, d) in enumerate(zip(np.asarray(f8r["eq_node"]), np.asarray(f8r["eq_dof"]))):
        nf = by_xyz[key(xyzr[int(n)])]
        err = worse(err, float(np.abs(np.asarray(f8r["H"])[:, k] - B.tf(f8f, nf, int(d))).max()))
    return err / scale, int(np.asarray(f8r["eq_node"]).size)


def constrained_in_full(full: Tuple[object, Dict[int, np.ndarray]], planes: Sequence[SYM.SymmetryPlane],
                        tol: float = 1e-6) -> float:
    """max |H_full| / max |H_full| of the DOFs that the symmetry constrains on the planes (they must vanish)."""
    f8f, xyzf = full
    scale = float(np.abs(f8f["H"]).max())
    worst = 0.0
    nodes = np.asarray(f8f["eq_node"])
    dofs = np.asarray(f8f["eq_dof"])
    H = np.asarray(f8f["H"])
    for k, (n, d) in enumerate(zip(nodes, dofs)):
        p = xyzf[int(n)]
        for pl in planes:
            if abs(pl.distance(p[None, :])[0]) <= tol and int(d) in pl.fixed_dofs:
                worst = worse(worst, float(np.abs(H[:, k]).max()))
    return worst / scale


@problem("VP-T2", "SYMM half and quarter models give the transfer functions of the full model", tier="P1",
         modules=["SITE", "POINT", "HOUSE", "ANALYS"], source="D-ANL-12, spec 07 9.2.39 (image method)")
def vpt2(workdir: Path) -> VPResult:
    """Doubly symmetric models (a flexible surface mat with a stick; an embedded FV basement with a stick)
    under vertically incident SV (X input) and P (Z input) are analysed as full, half (x >= 0, SYMM plane
    x = 0) and quarter (x, y >= 0, planes x = 0 and y = 0) models.  X input is antisymmetric about x = 0
    (type 1) and symmetric about y = 0 (type 0); Z input is symmetric about both.  HOUSE applies the
    symmetry boundary conditions on the plane nodes and ANALYS forms the soil impedance of the reduced
    interaction set from image nodes (D-ANL-12).  The reduced models must reproduce the full model at every
    DOF (1e-8), and the DOFs the symmetry fixes must vanish in the full model (1e-8)."""
    r = VPResult()
    site = B.layered_site(T2_SITE_LAYERS, T2_SITE_HS)
    fs = B.FrequencySet.fourier(0.01, 1024, T2_FS)
    for embedded in (False, True):
        what = "embedded FV basement + stick" if embedded else "surface mat + stick"
        for inp, wave, tx, ty in T2_INPUTS:
            res = {}
            for part in ("full", "half", "quarter"):
                wd = _fresh(Path(workdir) / f"{'emb' if embedded else 'mat'}_{inp}_{part}")
                mdl = symmetric_model(site, part, (tx, ty), embedded)
                B.run_soil(wd, "m", site, fs, layer=mdl.layer, rad=mdl.rad, wave=wave)
                f4 = B.run_house(wd, "m", mdl)
                B.run_analys(wd, "m", fs)
                xyz = {int(n): np.asarray(p, float) for n, p in zip(np.asarray(f4["node_id"]), np.asarray(f4["node_xyz"]))}
                res[part] = (B.read_file8(wd), xyz, int(np.asarray(f4["int_node"]).size))
                if part != "full":
                    r.require(f"{what}, {inp} input, {part} model: FILE4 carries the SYMM planes (x_symm) and ANALYS lists "
                              "the image terms", "x_symm" in f4 and "image terms" in B.listing(wd, "m", "ANALYS"))
                shutil.rmtree(wd, ignore_errors=True)
            for part in ("half", "quarter"):
                err, ndof = compare_reduced(res["full"][:2], res[part][:2])
                r.check(f"{what}, {inp} input: {part} model ({res[part][2]} of {res['full'][2]} interaction nodes, "
                        f"{ndof} DOFs) vs full model: max |H_red - H_full| / max |H_full|", err, 0.0, atol=1e-8)
            planes = [SYM.SymmetryPlane(1, tx, 0, 0.0), SYM.SymmetryPlane(2, ty, 1, 0.0)]
            r.check(f"{what}, {inp} input: full-model DOFs constrained by the symmetry types ({tx}, {ty}) vanish",
                    constrained_in_full(res["full"][:2], planes), 0.0, atol=1e-8)
    r.notes.append("Elements and masses lying in a symmetry plane (the stick) carry 1/2 (half) or 1/4 (quarter) of the "
                   "full-model section, masses and rotary inertias; everything else is the corresponding part of the "
                   "full model.")
    return r


# ======================================================================================
# VP-T3  2D zero-SSI identity
# ======================================================================================
T3_FS = (4, 10, 20, 41, 82, 123)


@problem("VP-T3", "2D zero-SSI identity: embedded PLANE excavation with structure = soil gives U = U'_f", tier="P1",
         modules=["SITE", "POINT", "HOUSE", "ANALYS"], source="R1 4.4 (built-in check), 02 Part 6")
def vpt3(workdir: Path) -> VPResult:
    """2D (plane strain) excavation 6 m wide through the two upper layers (PLANE elements, 4 columns, FV:
    every excavation node interacts) whose structure is made of PLANE elements with type-3 materials equal
    to the layers (K*_s = K*_e, M_s = M_e).  The flexible-volume equation reduces to X_ff (U - U'_f) = 0, so
    U = U'_f at every node, for vertical SV (X case) and P (Z case) input (<simul> = 1 runs the in-plane
    cases X and Z) and for an inclined SV wave (30 deg, phase exp(-i k x') along the model)."""
    r = VPResult()
    site = B.layered_site(T2_SITE_LAYERS, T2_SITE_HS)
    fs = B.FrequencySet.fourier(0.01, 1024, T3_FS)
    wd = _fresh(Path(workdir) / "t3")
    mdl = B.embedded_plane_2d(site, half_width=3.0, n_emb=2, ndiv=4, structure="soil")
    B.run_site_xyz(wd, "m", site, fs)
    B.write_deck(wd, "m", B.point_deck(fs, mdl.layer, mdl.rad, model="m", dim=1))
    B.run("POINT", wd, "m")
    f4 = B.run_house(wd, "m", mdl)
    r.require("2D FILE4: interaction DOFs UX, UZ (int_dofs [1, 3]) and no UY interaction equation",
              list(f4.meta["int_dofs"]) == [1, 3] and bool(np.all(np.asarray(f4["int_eq"])[:, 1] == -1)))
    B.run_analys(wd, "m", fs, simul=1)
    r.require("<simul> = 1 in 2D writes FILE8X and FILE8Z only (anti-plane Y case skipped)",
              (wd / "FILE8X").exists() and (wd / "FILE8Z").exists() and not (wd / "FILE8Y").exists())
    nodes = np.asarray(f4["int_node"])
    xyz = np.asarray(f4["x_int_xyz"], float)
    iface = np.asarray(f4["int_iface"])

    def identity(f8name: str, f1name: str, ang: float = 0.0) -> float:
        f8 = B.read_file8(wd, f8name)
        f1 = read_container(wd / f1name, "FILE1")
        err, scale = 0.0, 0.0
        for q, n in enumerate(fs.fnum):
            Up = free_field_at_nodes(f1, frequency_row(f1, n), xyz, iface, ang, 0.0, 0.0)
            H = np.array([[B.tf(f8, int(nd), d)[q] for d in (1, 3)] for nd in nodes])
            err = worse(err, float(np.abs(H - Up[:, [0, 2]]).max()))
            scale = worse(scale, float(np.abs(Up).max()))
        return err / scale

    for case, name in (("X", "vertical SV"), ("Z", "vertical P")):
        r.check(f"{name} ({case}): max |U - U'_f| / max |U'_f| over all nodes (UX, UZ) and frequencies",
                identity(f"FILE8{case}", f"FILE1{case}"), 0.0, atol=1e-8)
        Hc = B.tf(B.read_file8(wd, f"FILE8{case}"), mdl["top_centre"], 1 if case == "X" else 3)
        r.check(f"{name} ({case}): control point (surface centre) max |ATF - 1|", float(np.abs(Hc - 1.0).max()),
                0.0, atol=1e-8)
    # inclined SV wave in the X-Z plane (x' = X): a travelling free field exp(-i k x)
    B.write_deck(wd, "m", B.site_deck(site, fs, waves=[(2, 1, 1.0, 1.0, 30.0)], cm=0, wopt=0, model="m"))
    B.run("SITE", wd, "m")
    B.run_analys(wd, "m", fs)
    f1 = read_container(wd / "FILE1", "FILE1")
    kmax = float(np.abs(np.asarray(f1["k"])).max())
    r.require("inclined SV: the free field travels along X (k > 0)", kmax > 0.0)
    r.check("inclined SV (30 deg): max |U - U'_f| / max |U'_f| over all nodes and frequencies",
            identity("FILE8", "FILE1"), 0.0, atol=1e-8)
    shutil.rmtree(wd, ignore_errors=True)
    return r
