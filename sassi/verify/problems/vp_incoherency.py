"""Verification problems of incoherent seismic input, wave passage and multiple excitation (HOUSE
FILE77 and the ANALYS free-field load / motion; requirements 4.4 item 6, 4.6 item 5, 6.3).

* VP-26   FFL vs FFM for a rigid massless surface foundation on rock (05c C6): identical where the
          theory makes them identical -- exactly for coherent factors, and in the rigid limit for the
          generalised force in the input direction and for a decoupled foundation DOF (vertical input
          of a doubly symmetric mat); the other components are reported (informative)
* VP-27   coherency decomposition (05b section 11): trace check Eq. 6.1 for models 1 and 7, the
          f -> 0 limit (lambda_1 = N, AS factors = 1, zero-frequency ATF 1.00), ensemble mean of
          s s^H -> Sigma, directional distances sqrt(2 alpha) / sqrt(2 (1 - alpha)) through HOUSE, the
          wave-passage phase difference w L / V_app
* VP-I1   rigid massless surface foundation on a uniform half-space under vertical SH with Luco-Wong
          coherency (deterministic AS): ATF = 1 at f -> 0, decreasing with frequency, equal to the
          independent rigid-foundation average U = K_G^-1 T^T X_ff (s .* U'_f) (FFM) and
          K_G^-1 T^T (s .* X_ff U'_f) (FFL) evaluated from FILE3 / FILE1 / FILE77 (1e-8)
* VP-I2   stochastic simulation: the mean of s s^H over Ns samples converges to Sigma with an error
          proportional to 1/sqrt(Ns) (the exact variance of the random-phase estimator is known)

A *rigid* massless foundation is modelled, as in VP-15, by stiff massless BEAMS from the mat nodes to
a centre node.  For coherent input the uniform free field is an exact rigid-body solution, but the
incoherent load deforms the links by O(K_soil / K_link): the rigid limit is obtained by Richardson
extrapolation in the link stiffness (factors k, 2k, 4k; the response is an analytic function of 1/k,
so the 3-point extrapolation leaves an O(1/k^3) error), and the raw O(1/k) deviation is reported.
Units: m, t, kN, s (g = 9.81 m/s2).
"""
from __future__ import annotations

import math
import shutil
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from ...core import coherency as COH
from ...core import incoherency as INC
from ...core import ssi_solver as ss
from ...core.flexibility import flexibility_matrix, frequency_row
from ...core.freefield import free_field_at_nodes
from ...io.container import read_container
from .. import VPResult, problem, worse
from .. import builders as B

#: link stiffness factors of the rigid-limit extrapolation (x the largest soil shear modulus)
RIGID_FACTORS = (1.0e5, 2.0e5, 4.0e5)


def _fresh(path: Path) -> Path:
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    return path


def richardson3(h1: np.ndarray, h2: np.ndarray, h4: np.ndarray) -> np.ndarray:
    """Rigid limit from the results at link stiffnesses k, 2k, 4k: ``(8 h4 - 6 h2 + h1) / 3`` removes
    the O(1/k) and O(1/k^2) terms."""
    return (8.0 * np.asarray(h4) - 6.0 * np.asarray(h2) + np.asarray(h1)) / 3.0


def house_incoherent(wd: Path, mdl: B.Model, fs: B.FrequencySet, model: str = "m", **params) -> None:
    """HOUSE deck of ``mdl`` with the frequency set of ``fs`` and the incoherency deck parameters
    ``params`` (coh, wpass, me, cohf, alpha, gammax ...), then run HOUSE."""
    d = fs.fill(mdl.deck(model))
    for k, v in params.items():
        if k not in d.params:
            raise KeyError(f"unknown HOUSE parameter {k}")
        d[k] = v
    B.write_deck(wd, model, d)
    B.run("HOUSE", wd, model)


def luco_wong_sigma(xy: np.ndarray, f: float, gamma: float, vs: float) -> np.ndarray:
    """Independent Luco-Wong coherency matrix exp[-(gamma w D / Vs)^2] (own code, not sassi.core)."""
    d = xy[:, None, :] - xy[None, :, :]
    D = np.sqrt((d ** 2).sum(axis=2))
    return np.exp(-(gamma * 2.0 * math.pi * f * D / vs) ** 2)


def as_factors(S: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Independent AS synthesis of D-INC-03 (own code, written from the documented rules of
    sassi.core.incoherency): eigen-decomposition, lambda descending, negative round-off clipped; every
    group of equal eigenvalues (within 1e-9 lambda_1, above 1e-12 lambda_1) re-based on its eigenspace:
    first the normalised projection of the uniform field 1, then, in node order, the normalised
    projections of the unit vectors with the largest remaining norm (pivoted Gram-Schmidt); every mode
    oriented so that its sum is >= 0 (first significant component > 0 for a vanishing sum) -> (s, lambda)."""
    w, V = np.linalg.eigh(S)
    lam = np.clip(w[::-1], 0.0, None)
    V = V[:, ::-1].copy()
    n = lam.size
    i = 0
    while i < n:
        j = i + 1
        while j < n and abs(lam[j - 1] - lam[j]) <= 1e-9 * lam[0]:
            j += 1
        if j - i > 1 and lam[i] > 1e-12 * lam[0]:
            P = V[:, i:j] @ V[:, i:j].T                 # projector on the eigenspace
            basis = []
            u = P @ np.ones(n)
            if np.linalg.norm(u) > 1e-8 * math.sqrt(n):
                basis.append(u / np.linalg.norm(u))
            while len(basis) < j - i:
                R = P.copy()
                for b in basis:
                    R -= np.outer(b, b)
                norms = np.einsum("ij,ij->j", R, R)
                k = int(np.flatnonzero(norms >= norms.max() * (1.0 - 1e-9))[0])
                basis.append(R[:, k] / math.sqrt(norms[k]))
            V[:, i:j] = np.column_stack(basis)
        i = j
    for k in range(V.shape[1]):
        col = V[:, k]
        tot = col.sum()
        if abs(tot) > 1e-9 * np.abs(col).sum():
            sgn = 1.0 if tot > 0 else -1.0
        else:
            big = np.flatnonzero(np.abs(col) > 1e-6 * np.abs(col).max())
            sgn = -1.0 if col[big[0]] < 0 else 1.0
        V[:, k] *= sgn
    return V @ np.sqrt(lam), lam


def random_phase_variance(S: np.ndarray) -> float:
    """Exact sum over all entries of the variance of one sample of ``s s^H`` for uniform random modal
    phases on the full circle: ``Var_ij = Sigma_ii Sigma_jj - sum_k lambda_k^2 phi_ik^2 phi_jk^2``
    (VP-I2; the mean of Ns independent samples has Var/Ns)."""
    w, V = np.linalg.eigh(S)
    lam = np.clip(w, 0.0, None)
    a2 = (V ** 2) * lam[None, :]                 # a_ik^2 = lambda_k phi_ik^2
    diag = np.diag(S)
    return float(np.sum(np.outer(diag, diag)) - np.sum((a2 @ a2.T)))


def _rigid_reference(f1, f3, f77, xyz: np.ndarray, iface: np.ndarray, fnum: Sequence[int], ref: np.ndarray,
                     direction: int, ffm: bool) -> np.ndarray:
    """Independent rigid massless foundation response (6 rigid-body DOFs per frequency) from FILE3
    (flexibility -> X_ff), FILE1 (U'_f), FILE77 (s):  q = K_G^-1 T^T b,  K_G = T^T X_ff T,
    b = X_ff (s .* U'_f) (FFM) or s .* (X_ff U'_f) (FFL)."""
    T = ss.rigid_body_transform(xyz, ref)
    rows77 = {int(n): k for k, n in enumerate(np.asarray(f77["fnum"]))}
    rows1 = {int(n): k for k, n in enumerate(np.asarray(f1["fnum"]))}
    out = []
    with np.errstate(all="ignore"):              # spurious BLAS floating-point flags (numpy 2 / Accelerate)
        for n in fnum:
            F = flexibility_matrix(f3, frequency_row(f3, int(n)), xyz[:, :2], iface)
            X = np.linalg.inv(F)
            X = 0.5 * (X + X.T)
            Up = free_field_at_nodes(f1, rows1[int(n)], xyz, iface, 0.0, 0.0, 0.0).reshape(-1)
            s3 = np.repeat(np.asarray(f77["s"])[rows77[int(n)], direction], 3)
            b = X @ (s3 * Up) if ffm else s3 * (X @ Up)
            KG = T.T @ X @ T
            out.append(np.linalg.solve(KG, T.T @ b))
    return np.asarray(out)


# ======================================================================================
# VP-I1  rigid surface foundation, Luco-Wong, vertical SH
# ======================================================================================
I1_SOIL = dict(depth=16.0, nsub=4, vs=250.0, rho=2.0, nu=1.0 / 3.0, beta=0.02)
I1_GAMMA = 0.5
I1_FNUM = [1, 10, 20, 40, 60, 80, 100, 120, 150]


@problem("VP-I1", "Rigid surface foundation under incoherent vertical SH (Luco-Wong, AS): ATF vs the "
         "rigid-foundation formula K_G^-1 T^T X_ff (s .* U'_f) from FILE3/FILE1/FILE77", tier="P1",
         modules=["SITE", "POINT", "HOUSE", "ANALYS"], source="Luco & Wong (1986); requirements 4.4 item 6, 4.6 item 5")
def vp_i1(workdir: Path) -> VPResult:
    """Square 10 m x 10 m rigid massless mat (5 x 5 interaction nodes) on a uniform half-space (Vs = 250 m/s),
    vertically incident SH (y'), Luco-Wong coherency gamma = 0.5 with Vs = 250 m/s, deterministic AS factors
    (D-INC-03).  HOUSE's FILE77 is recomputed independently (own coherency matrix and eigen-decomposition);
    the ANALYS centre-node ATF (FFL and FFM, rigid limit) is compared with the rigid-foundation formula
    evaluated from FILE3 / FILE1 / FILE77.  That formula is built with the library functions ANALYS also
    uses for X_ff (``flexibility_matrix``) and U'_f (``free_field_at_nodes``), which are verified by VP-07 ...
    VP-16: the comparison checks the ANALYS incoherent load assembly, solve and rigid limit, not X_ff or
    U'_f themselves (final audit).  The f -> 0 and monotonicity criteria apply to the ANALYS result."""
    r = VPResult()
    site = B.uniform_site(I1_SOIL["depth"], I1_SOIL["nsub"], I1_SOIL["vs"], I1_SOIL["rho"], I1_SOIL["nu"],
                          I1_SOIL["beta"])
    fs = B.FrequencySet.fourier(0.01, 2048, I1_FNUM)
    wd = _fresh(Path(workdir) / "i1")
    H: Dict[Tuple[float, int], np.ndarray] = {}
    refs: Dict[int, np.ndarray] = {}
    for k, fac in enumerate(RIGID_FACTORS):
        mdl = B.surface_rigid_mat(site, half_width=5.0, ndiv=4, E_rigid=fac * site.max_G())
        if k == 0:
            B.run_soil(wd, "m", site, fs, layer=0, rad=mdl.rad, wave="SH")
        house_incoherent(wd, mdl, fs, coh=1, cohf=1, alpha=I1_SOIL["vs"], gammax=I1_GAMMA, gammay=I1_GAMMA,
                         gammaz=I1_GAMMA)
        c = mdl["centre"]
        for ffm in (0, 1):
            B.run_analys(wd, "m", fs, coh=1, ffm=ffm)
            f8 = B.read_file8(wd)
            H[(fac, ffm)] = B.tf(f8, c, 2)
        if k == 0:
            f1 = read_container(wd / "FILE1", "FILE1")
            f3 = read_container(wd / "FILE3", "FILE3")
            f4 = read_container(wd / "m.N4", "FILE4")
            f77 = read_container(wd / "FILE77", "FILE77")
            xyz = np.asarray(f4["x_int_xyz"], float)
            iface = np.asarray(f4["int_iface"])
            ctr = np.asarray(mdl.house.nodes[c], float)
            # (a) FILE77 = an independent AS synthesis of the Luco-Wong coherency matrix
            worst77 = 0.0
            for q, f in enumerate(fs.freq):
                s_ind, _ = as_factors(luco_wong_sigma(xyz[:, :2], f, I1_GAMMA, I1_SOIL["vs"]))
                worst77 = worse(worst77, float(np.abs(np.asarray(f77["s"])[q, 1] - s_ind).max()))
            r.check("FILE77 (Y) vs independent Luco-Wong AS factors, max |difference| over all frequencies",
                    worst77, 0.0, atol=1e-10)
            for ffm in (0, 1):
                refs[ffm] = _rigid_reference(f1, f3, f77, xyz, iface, fs.fnum, ctr, 1, bool(ffm))[:, 1]
    for ffm in (0, 1):
        tag = "FFM" if ffm else "FFL"
        ref = refs[ffm]
        rigid = richardson3(*(H[(fac, ffm)] for fac in RIGID_FACTORS))
        rel = np.abs(rigid - ref) / np.abs(ref)
        r.check(f"{tag}: rigid-limit centre ATF_y vs K_G^-1 T^T b (FILE3/FILE1/FILE77), max relative difference",
                float(rel.max()), 0.0, atol=1e-8)
        raw = np.abs(H[(RIGID_FACTORS[0], ffm)] - ref) / np.abs(ref)
        r.inform(f"{tag}: link stiffness {RIGID_FACTORS[0]:g} G (no extrapolation), max relative difference "
                 "(O(K_soil/K_link))", float(raw.max()), 0.0)
        e1 = np.abs(H[(RIGID_FACTORS[0], ffm)] - ref).max()
        e2 = np.abs(H[(RIGID_FACTORS[1], ffm)] - ref).max()
        r.check(f"{tag}: error ratio k / 2k of the raw results (1/k convergence)", e1 / e2, 2.0, rtol=0.05)
        a = np.abs(rigid)                        # the ANALYS result (final audit: was the VP's own ref)
        r.check(f"{tag}: |ATF_y| at f = {fs.freq[0]:.4g} Hz (f -> 0 limit)", float(a[0]), 1.0, atol=1e-4)
        r.require(f"{tag}: |ATF_y| decreases monotonically over {fs.freq[0]:.3g} .. {fs.freq[-1]:.3g} Hz "
                  f"({' > '.join(f'{v:.4f}' for v in a)})", bool(np.all(np.diff(a) < 0)))
    r.inform("FFL / FFM centre |ATF_y| at the highest frequency",
             float(abs(refs[0][-1])), float(abs(refs[1][-1])))
    r.notes.append(f"Rigid mat 10 m x 10 m, Vs = 250 m/s, gamma = {I1_GAMMA}: |ATF_y| from {abs(refs[0][0]):.6f} at "
                   f"{fs.freq[0]:.3g} Hz to {abs(refs[0][-1]):.4f} at {fs.freq[-1]:.3g} Hz (FFL).")
    return r


# ======================================================================================
# VP-I2  stochastic simulation: ensemble mean of s s^H
# ======================================================================================
I2_NS = (2, 5, 10, 20, 50)


def ensemble_errors(files: List, xy: np.ndarray, freq: np.ndarray, gamma: float, vs: float,
                    ns_list: Sequence[int]) -> Tuple[np.ndarray, np.ndarray]:
    """Mean squared Frobenius error of the Ns-sample mean of s s^H (samples 1..Ns) against the
    independent Luco-Wong Sigma, averaged over frequencies and directions, and its exact expectation
    (sum of Var_ij)/Ns (:func:`random_phase_variance`)."""
    S = np.stack([np.asarray(c["s"]) for c in files])         # (Ns, nF, 3, N)
    mse = np.zeros(len(ns_list))
    theo = np.zeros(len(ns_list))
    cnt = 0
    for q, f in enumerate(freq):
        Sig = luco_wong_sigma(xy, f, gamma, vs)
        var = random_phase_variance(Sig)
        for d in range(3):
            prod = np.einsum("ri,rj->rij", S[:, q, d], np.conj(S[:, q, d]))
            csum = np.cumsum(prod, axis=0)
            for j, ns in enumerate(ns_list):
                est = csum[ns - 1] / ns
                mse[j] += float(np.sum(np.abs(est - Sig) ** 2))
                theo[j] += var / ns
            cnt += 1
    return mse / cnt, theo / cnt


@problem("VP-I2", "Stochastic incoherency simulation: mean of s s^H over Ns samples -> coherency matrix, error "
         "proportional to 1/sqrt(Ns)", tier="P1", modules=["HOUSE"], source="requirements 4.4 item 6, D-INC-07; 05b 11 test 3")
def vp_i2(workdir: Path) -> VPResult:
    """4 x 4 surface grid of interaction nodes (spacing 4 m) with Luco-Wong coherency (gamma = 0.3,
    Vs = 200 m/s) at 10 frequencies, HOUSE stochastic simulation (HSeed 11, VSeed 17, RandPhz 180, 50
    samples -> FILE77001..FILE77050).  For uniform random phases on the full circle the Ns-sample mean of
    s s^H is unbiased with the exact variance sum_ij Var_ij / Ns, Var_ij = 1 - sum_k lambda_k^2
    phi_ik^2 phi_jk^2; the observed mean squared error is compared with it and its log-log slope with -1."""
    r = VPResult()
    site = B.uniform_site(10.0, 2, vs=200.0, rho=2.0, nu=1.0 / 3.0, beta=0.03)
    fs = B.FrequencySet.fourier(0.01, 1024, [10, 30, 50, 70, 90, 110, 130, 150, 170, 190])
    mdl = B.surface_rigid_mat(site, half_width=6.0, ndiv=3, centre_height=1.0)
    wd = _fresh(Path(workdir) / "i2")
    B.write_deck(wd, "m", B.site_deck(site, fs, model="m"))
    params = dict(coh=1, cohf=1, alpha=200.0, gammax=0.3, gammay=0.3, gammaz=0.3, hseed=11, vseed=17,
                  randphz=180.0, nsim=50)
    house_incoherent(wd, mdl, fs, **params)
    files = [read_container(wd / INC.file77_name(k), "FILE77") for k in range(1, 51)]
    r.require("HOUSE wrote the 50 stochastic samples FILE77001 .. FILE77050 (no deterministic FILE77)",
              all(int(c.meta.get("stochastic", 0)) == 1 for c in files) and not (wd / "FILE77").exists())
    f4 = read_container(wd / "m.N4", "FILE4")
    xy = np.asarray(f4["x_int_xyz"], float)[:, :2]
    mse, theo = ensemble_errors(files, xy, fs.freq, 0.3, 200.0, I2_NS)
    for ns, m_, t_ in zip(I2_NS, mse, theo):
        r.check(f"Ns = {ns}: mean squared error of mean(s s^H) / exact expectation (sum Var_ij)/Ns", m_ / t_, 1.0,
                rtol=0.25)
    slope = float(np.polyfit(np.log(I2_NS), np.log(np.sqrt(mse)), 1)[0])
    r.check("log-log slope of the RMS error of mean(s s^H) vs Ns (1/sqrt(Ns) -> -0.5)", slope, -0.5, atol=0.1)
    r.inform("RMS error of the 50-sample mean relative to ||Sigma||_F (averaged)",
             float(np.sqrt(mse[-1]) / xy.shape[0]), 0.0)
    # reproducibility (D-GEN-09, D-INC-07) and independence of the samples
    s1 = np.asarray(files[0]["s"]).copy()
    house_incoherent(wd, mdl, fs, **params)
    again = np.asarray(read_container(wd / INC.file77_name(1), "FILE77")["s"])
    r.check("re-run of HOUSE with the same seeds: FILE77001 identical (max |difference|)",
            float(np.abs(again - s1).max()), 0.0, atol=0.0)
    r.require("FILE77001 and FILE77002 are different samples", bool(np.abs(s1 - np.asarray(files[1]["s"])).max() > 1e-3))
    r.notes.append("MSE ratios " + ", ".join(f"Ns={n}: {m_ / t_:.3f}" for n, m_, t_ in zip(I2_NS, mse, theo))
                   + f"; slope {slope:.3f}")
    return r


# ======================================================================================
# VP-27  decomposition checks
# ======================================================================================
def _user_tables(wd: Path, freq: Sequence[float], dist: Sequence[float], fun) -> None:
    """Write FREQCOH, DISTCOH and COHXUSER/COHYUSER/COHZUSER with gamma = fun(f, D, direction)."""
    (wd / COH.USER_FREQ).write_text(" ".join(f"{v:.10g}" for v in freq) + "\n", encoding="utf-8")
    (wd / COH.USER_DIST).write_text(" ".join(f"{v:.10g}" for v in dist) + "\n", encoding="utf-8")
    for k, name in enumerate(COH.USER_TABLES):
        rows = [" ".join(f"{fun(f, D, k):.15e}" for D in dist) for f in freq]
        (wd / name).write_text("\n".join(rows) + "\n", encoding="utf-8")


def _pair_model(site: B.Site, p1, p2) -> B.Model:
    """Two surface interaction nodes joined by a stiff massless beam (the smallest HOUSE model)."""
    hb = B.HouseBuilder(site, title="two interaction nodes")
    a = hb.node(p1[0], p1[1], site.gelev)
    b = hb.node(p2[0], p2[1], site.gelev)
    B.rigid_link(hb, a, b, 1.0)
    hb.set_interaction([a, b])
    return B.Model(hb, {"a": a, "b": b}, rad=0.5, layer=0)


@problem("VP-27", "Incoherency decomposition: Eq. 6.1 trace, f -> 0 limit, ensemble mean, directional distances, "
         "wave-passage phase", tier="P1", modules=["HOUSE", "ANALYS"], source="05b section 11; requirements 4.4 item 6")
def vp27(workdir: Path) -> VPResult:
    r = VPResult()
    site = B.uniform_site(12.0, 3, vs=300.0, rho=2.0, nu=1.0 / 3.0, beta=0.03)
    # ---- (1) Eq. 6.1 trace check on a 5 x 5 grid, model 1 and model 7 ------------------------------
    fs = B.FrequencySet.fourier(0.01, 8192, [1, 2, 100, 400, 800, 1200, 1600])    # 0.0122 .. 19.5 Hz
    wd = _fresh(Path(workdir) / "grid")
    mdl = B.surface_rigid_mat(site, half_width=6.0, ndiv=4)
    B.write_deck(wd, "m", B.site_deck(site, fs, model="m"))
    house_incoherent(wd, mdl, fs, coh=1, cohf=1, alpha=300.0, gammax=0.3, gammay=0.3, gammaz=0.5, ipr=1)
    f77 = read_container(wd / "FILE77", "FILE77")
    f4 = read_container(wd / "m.N4", "FILE4")
    xy = np.asarray(f4["x_int_xyz"], float)[:, :2]
    N = xy.shape[0]
    r.check(f"model 1: Eq. 6.1 |sum(lambda) - N|/N, worst of {len(fs.fnum)} frequencies x 3 directions (FILE77)",
            float(np.asarray(f77["x_trace_error"]).max()), 0.0, atol=1e-8)
    worst = 0.0
    for f in fs.freq:
        for g in (0.3, 0.5):
            lam = np.linalg.eigvalsh(luco_wong_sigma(xy, f, g, 300.0))
            worst = worse(worst, abs(float(np.clip(lam, 0, None).sum()) - N) / N)
    r.check("model 1: Eq. 6.1 recomputed independently (numpy eigvalsh)", worst, 0.0, atol=1e-8)
    r.require("HOUSE listing has the I N C O tables (INCOH <ipr> = 1)", "I N C O" in B.listing(wd, "m", "HOUSE"))
    # f -> 0: Sigma -> 1 1^T, lambda_1 = N, AS factors -> 1 (|s - 1| = O(f): the higher modes enter with
    # sqrt(lambda_k) = O(f) since 1 - Sigma_ij = O(f^2))
    lam1 = float(np.asarray(f77["x_lambda1_pct"])[0].min())
    r.check(f"f = {fs.freq[0]:.4g} Hz: lambda_1 / N (%)", lam1, 100.0, atol=1e-3)
    s_all = np.asarray(f77["s"])
    dev1, dev2 = float(np.abs(s_all[0] - 1.0).max()), float(np.abs(s_all[1] - 1.0).max())
    r.check(f"f = {fs.freq[0]:.4g} Hz: AS factors, max |s - 1| (X, Y, Z)", dev1, 0.0, atol=0.01)
    r.check(f"AS factors -> 1 linearly as f -> 0: |s - 1| at {fs.freq[1]:.4g} Hz / at {fs.freq[0]:.4g} Hz",
            dev2 / dev1, 2.0, rtol=0.02)
    # zero-frequency ATF of the rigid mat (manual check: ~1.00 within 5 %)
    B.write_deck(wd, "m", B.point_deck(fs, 0, mdl.rad, model="m"))
    B.run("SITE", wd, "m")
    B.run("POINT", wd, "m")
    B.run_analys(wd, "m", fs, coh=1)
    f8 = B.read_file8(wd)
    a0 = np.array([abs(B.tf(f8, n, 1)[0]) for n in mdl["mat"]])
    r.check(f"zero-frequency ATF (f = {fs.freq[0]:.4g} Hz, X, all mat nodes): max ||ATF| - 1|",
            float(np.abs(a0 - 1.0).max()), 0.0, atol=0.05)
    # model 7: tabulated exp(-c f D) (a positive definite kernel); the same checks, and vs the formula
    cexp = 0.004
    tf_grid = np.linspace(0.0, 25.0, 51)
    td_grid = np.linspace(0.0, 20.0, 81)
    _user_tables(wd, tf_grid, td_grid, lambda f, D, k: math.exp(-cexp * (1.0 + 0.5 * k) * f * D))
    house_incoherent(wd, mdl, fs, coh=1, wpass=1, cohf=7, alpha=0.5, ipr=0)
    f77u = read_container(wd / "FILE77", "FILE77")
    r.check("model 7: Eq. 6.1 |sum(lambda) - N|/N, worst (FILE77)", float(np.asarray(f77u["x_trace_error"]).max()),
            0.0, atol=1e-8)
    # the eigenvalues are well conditioned (unlike the eigenvectors of near-equal eigenvalues): compare the
    # share of the first mode with the one of the tabulated function itself (bilinear interpolation error)
    dev = 0.0
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(axis=2))
    for q, f in enumerate(fs.freq):
        for k in range(3):
            lam = np.linalg.eigvalsh(np.exp(-cexp * (1.0 + 0.5 * k) * f * d))
            dev = worse(dev, abs(float(np.asarray(f77u["x_lambda1_pct"])[q, k]) - 100.0 * float(lam[-1]) / N))
    r.check("model 7 (bilinear tables) vs the tabulated function exp(-c f D): lambda_1 / N in %, max |difference| "
            "(X, Y, Z)", dev, 0.0, atol=0.05)
    # ---- (2) ensemble mean of s s^H -> Sigma (VP-I2 in detail) ------------------------------------
    fs2 = B.FrequencySet.fourier(0.01, 8192, [400, 800, 1200])
    house_incoherent(wd, mdl, fs2, coh=1, cohf=1, alpha=300.0, gammax=0.3, gammay=0.3, gammaz=0.3,
                     hseed=5, vseed=9, randphz=180.0, nsim=20)
    files = [read_container(wd / INC.file77_name(k), "FILE77") for k in range(1, 21)]
    mse, theo = ensemble_errors(files, xy, fs2.freq, 0.3, 300.0, (5, 20))
    r.check("ensemble of 20 samples: mean squared error of mean(s s^H) / exact expectation", mse[1] / theo[1], 1.0,
            rtol=0.3)
    r.check("error ratio Ns = 5 vs Ns = 20 (1/sqrt(Ns): ratio 2)", float(np.sqrt(mse[0] / mse[1])), 2.0, rtol=0.25)
    # ---- (3) directional distances through HOUSE (model 7, gamma = 1 - D/100) ----------------------
    fs3 = B.FrequencySet.fourier(0.01, 1024, [100])
    wd3 = _fresh(Path(workdir) / "pairs")
    _user_tables(wd3, [0.0, 50.0], [0.0, 100.0], lambda f, D, k: 1.0 - D / 100.0)
    B.write_deck(wd3, "m", B.site_deck(site, fs3, model="m"))
    L, alpha, ang = 10.0, 0.1, 30.0
    ca, sa = math.cos(math.radians(ang)), math.sin(math.radians(ang))
    dists = {}
    for name, (ux, uy) in (("X'", (ca, sa)), ("Y'", (-sa, ca))):
        pm = _pair_model(site, (0.0, 0.0), (L * ux, L * uy))
        house_incoherent(wd3, pm, fs3, coh=1, wpass=1, cohf=7, alpha=alpha, wang=ang)
        s = np.asarray(read_container(wd3 / "FILE77", "FILE77")["s"])[0, 0].real
        g = float(s[0] * s[1])                  # AS of [[1, g], [g, 1]]: s1 s2 = g
        dists[name] = 100.0 * (1.0 - g)
    r.check(f"directional distance along Line D (alpha = {alpha}, L = {L:g}, Line D at {ang:g} deg): D / L",
            dists["X'"] / L, math.sqrt(2.0 * alpha), rtol=1e-9)
    r.check("directional distance across Line D: D / L", dists["Y'"] / L, math.sqrt(2.0 * (1.0 - alpha)), rtol=1e-9)
    r.check("ratio of the two distances at alpha = 0.1 (Y' separations weighted 3 times)", dists["X'"] / dists["Y'"],
            1.0 / 3.0, rtol=1e-9)
    # ---- (4) wave-passage phase difference w L / V_app (coherent motion, FILE77 of wave passage) ---
    vapp = 400.0
    pm = _pair_model(site, (0.0, 0.0), (L * ca, L * sa))
    house_incoherent(wd3, pm, fs3, coh=0, wpass=1, appv=vapp, wang=ang)
    s = np.asarray(read_container(wd3 / "FILE77", "FILE77")["s"])[0, 0]
    w = 2.0 * math.pi * fs3.freq[0]
    dphi = float(np.angle(s[1] / s[0]))
    expect = -((w * L / vapp + math.pi) % (2.0 * math.pi) - math.pi)
    r.check(f"wave passage: phase of s_2 / s_1 for L = {L:g} m along Line D, V_app = {vapp:g} m/s "
            f"(-w L / V_app, wrapped)", dphi, expect, atol=1e-12)
    r.check("wave passage: |s| = 1 (pure delay)", float(np.abs(np.abs(s) - 1.0).max()), 0.0, atol=1e-14)
    return r


# ======================================================================================
# VP-26  FFL vs FFM, rigid surface foundation on rock
# ======================================================================================
ROCK = dict(depth=30.0, nsub=3, vs=2000.0, rho=2.5, nu=0.25, beta=0.01)
R26_FNUM = [2, 50, 100, 200, 300, 400]


@problem("VP-26", "FFL vs FFM for a rigid massless surface foundation on rock (hard-rock coherency, AS)", tier="P1",
         modules=["SITE", "POINT", "HOUSE", "ANALYS"], source="05c C6; requirements 4.6 item 5")
def vp26(workdir: Path) -> VPResult:
    """Square 20 m x 20 m rigid massless mat on a uniform rock half-space (Vs = 2000 m/s), X/Y/Z input
    (SV x', SH y', P z'), Abrahamson 2007 hard-rock coherency (isotropic, V_app = 1e9), deterministic AS.

    Reference (requirements VP-26, spec 05c C6, manual section 6 ANALYS options): FFL and FFM give
    identical results for a surface rigid foundation on rock.  Lead decision D-W3-05 refined it: the
    free-field load (FFL, ``s .* X U'``) and the free-field motion (FFM, ``X (s .* U')``) load vectors
    differ by the commutator ``[S, X] U'``.  For a rigid massless foundation the rigid-body response is
    ``q = K_G^-1 T^T b`` and ``T^T b_FFL - T^T b_FFM = T^T [S, X] T t_d``, an antisymmetric matrix times
    the input direction t_d: the generalised force in the input direction is identical (zero diagonal),
    and so is any foundation DOF that K_G does not couple to the others -- the vertical translation of a
    doubly symmetric mat -- and every DOF for coherent factors (s = 1).

    Criteria (D-W3-05): those identities (rigid limit, Richardson extrapolation in the link stiffness)
    and the exact relation ``K_G (q_FFL - q_FFM) = T^T [S, X] T t_d`` computed from FILE11 and FILE77.
    The horizontal ATF differences of FFL and FFM (about 1e-3 to 3e-2, coupled through K_G) are
    informative rows."""
    r = VPResult()
    rock = B.uniform_site(ROCK["depth"], ROCK["nsub"], ROCK["vs"], ROCK["rho"], ROCK["nu"], ROCK["beta"])
    fs = B.FrequencySet.fourier(0.005, 2048, R26_FNUM)
    wd = _fresh(Path(workdir) / "rock")
    B.run_site_xyz(wd, "m", rock, fs)
    res: Dict[Tuple[float, int, str], np.ndarray] = {}
    KG = None
    mdl = None
    for k, fac in enumerate(RIGID_FACTORS):
        mdl = B.surface_rigid_mat(rock, half_width=10.0, ndiv=4, E_rigid=fac * rock.max_G())
        if k == 0:
            B.write_deck(wd, "m", B.point_deck(fs, 0, mdl.rad, model="m"))
            B.run("POINT", wd, "m")
        house_incoherent(wd, mdl, fs, coh=1, wpass=1, cohf=5, alpha=0.5, appv=1.0e9, met=1)
        c = mdl["centre"]
        for ffm in (0, 1):
            B.run_analys(wd, "m", fs, coh=1, wpass=1, ffm=ffm, simul=1, impe=2)
            for case in "XYZ":
                f8 = B.read_file8(wd, f"FILE8{case}")
                res[(fac, ffm, case)] = np.array([B.tf(f8, c, dof) for dof in range(1, 7)]).T   # (nF, 6)
            if KG is None:
                KG = np.asarray(read_container(wd / "FILE11", "FILE11")["KG"])                  # (nF, 6, 6)
    rigid = {(ffm, case): richardson3(*(res[(fac, ffm, case)] for fac in RIGID_FACTORS))
             for ffm in (0, 1) for case in "XYZ"}
    # (1) generalised force in the input direction: (K_G q)_d identical
    for d, case in enumerate("XYZ"):
        gl = np.einsum("qij,qj->qi", KG, rigid[(0, case)])[:, d]
        gm = np.einsum("qij,qj->qi", KG, rigid[(1, case)])[:, d]
        r.check(f"{case} input: generalised force (K_G q)_{case} of FFL vs FFM (rigid limit), max relative difference",
                float(np.max(np.abs(gl - gm) / np.abs(gl))), 0.0, atol=1e-8)
    # (2) D-W3-05: the spec's "identical results" holds exactly only where the commutator vanishes.
    #     Criteria: (a) the decoupled vertical translation is identical for Z input; (b) for every input
    #     the FFL - FFM difference of the generalised forces equals T^T [S, X] T t_d computed independently
    #     from FILE11 (X_ff, T) and FILE77 (s) -- an exact identity of the two load paths.  The horizontal
    #     foundation ATF differences are reported (informative).
    f11 = read_container(wd / "FILE11", "FILE11")
    f77s = read_container(wd / "FILE77", "FILE77")
    if "X" not in f11.arrays:
        raise RuntimeError("FILE11 holds no X_ff (raise the FILE11 cap for this VP)")
    Xall, T = np.asarray(f11["X"]), np.asarray(f11["T"])
    if not np.array_equal(np.asarray(f11["int_node"]), np.asarray(f77s["int_node"])):
        raise RuntimeError("FILE11 and FILE77 interaction orders differ")
    sfac = np.asarray(f77s["s"])                                        # (nF, 3, nInt)
    diffs = {}
    for d, case in enumerate("XYZ"):
        L, M = rigid[(0, case)][:, d], rigid[(1, case)][:, d]
        diffs[case] = float(np.max(np.abs(L - M) / np.abs(L)))
        if case == "Z":
            r.check("Z input: foundation ATF_z of FFL vs FFM (rigid limit; vertical translation decoupled by K_G), "
                    "max relative difference", diffs[case], 0.0, atol=1e-8)
        else:
            r.inform(f"{case} input: foundation ATF_{case.lower()} of FFL vs FFM (rigid limit), max relative "
                     "difference (spec 05c C6 'identical'; differs through the sliding-rocking coupling of K_G)",
                     diffs[case], 0.0)
        td = np.zeros(6)
        td[d] = 1.0
        worst = 0.0
        for q in range(len(Xall)):
            S3 = np.repeat(sfac[q, d], 3)
            X = Xall[q]
            with np.errstate(all="ignore"):                              # spurious Accelerate-BLAS flags
                comm = S3[:, None] * X - X * S3[None, :]                 # [S, X]
                pred = T.T @ (comm @ (T @ td))
                got = KG[q] @ (rigid[(0, case)][q] - rigid[(1, case)][q])
                scale = np.linalg.norm(KG[q] @ rigid[(0, case)][q])
            e = float(np.linalg.norm(got - pred) / scale)
            worst = e if not np.isfinite(e) else max(worst, e)          # a NaN must fail, not vanish
            if not np.isfinite(worst):
                break
        r.check(f"{case} input: K_G (q_FFL - q_FFM) = T^T [S, X] T t_d (independent commutator from FILE11 and "
                "FILE77), max relative residual", worst, 0.0, atol=1e-6)
    raw = np.abs(res[(RIGID_FACTORS[0], 0, "Z")][:, 2] - res[(RIGID_FACTORS[0], 1, "Z")][:, 2])
    r.inform("Z input, link stiffness 1e5 G (no extrapolation): max relative FFL/FFM difference",
             float(np.max(raw / np.abs(res[(RIGID_FACTORS[0], 0, 'Z')][:, 2]))), 0.0)
    # rotations times the mat half-width (10 m): the edge displacement they cause, per unit control motion
    rot = max(float(np.max(np.abs(rigid[(0, case)][:, 3:] - rigid[(1, case)][:, 3:]))) for case in "XYZ")
    rot_ffl = max(float(np.max(np.abs(rigid[(0, case)][:, 3:]))) for case in "XYZ")
    r.inform("rotations x half-width 10 m, per unit control motion: max |FFL - FFM| of all inputs (reference: "
             "max |FFL| rotation x 10 m)", rot * 10.0, rot_ffl * 10.0)
    # (4) coherent factors s = 1: FFL = FFM = the coherent analysis, every DOF
    f77 = read_container(wd / "FILE77", "FILE77")
    meta = {"df": float(f77.meta["df"]), "method": "COHERENT", "stochastic": 0, "coh": 0, "wpass": 0, "me": 0}
    INC.write_file77(wd / "FILE77", f77["fnum"], f77["freq"], f77["int_node"], f77["x_int_xyz"],
                     np.ones_like(np.asarray(f77["s"])), meta, module="VP-26")
    out = {}
    for tag, params in (("FFL", dict(coh=1, ffm=0)), ("FFM", dict(coh=1, ffm=1)), ("coherent", dict(coh=0))):
        B.run_analys(wd, "m", fs, simul=1, **params)
        out[tag] = np.concatenate([np.asarray(B.read_file8(wd, f"FILE8{c}")["H"]) for c in "XYZ"], axis=1)
    ref = out["coherent"]
    den = float(np.abs(ref).max())
    for tag in ("FFL", "FFM"):
        r.check(f"factors s = 1: {tag} vs the coherent analysis, all DOFs, max |difference| / max |H|",
                float(np.abs(out[tag] - ref).max()) / den, 0.0, atol=1e-12)
    r.notes.append(f"D-W3-05: FFL vs FFM foundation ATF differs by {diffs['X']:.3g} "
                   f"(X input) and {diffs['Y']:.3g} (Y input), {diffs['Z']:.3g} for Z input. Independent evidence "
                   "that the manual's 'identical results' is approximate: T^T [S, X] T is antisymmetric, so FFL = FFM "
                   "exactly only for the generalised force in the input direction, K_G-decoupled DOFs (the vertical "
                   "translation) and coherent factors -- all verified above, together with the exact commutator "
                   "identity of the difference; the horizontal translations and rotations are coupled through K_G "
                   "and differ (informative rows).")
    return r
