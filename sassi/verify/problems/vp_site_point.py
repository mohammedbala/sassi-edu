"""Verification problems of SITE, POINT and the layered-soil core (requirements §6.3, R2 A/D, R1 §9).

Every problem writes its decks with :func:`sassi.io.decks.new`, runs the modules through
:func:`sassi.modules.base.run_module` and evaluates FILE1/FILE2/FILE3 through the binding APIs
(:func:`sassi.core.freefield.free_field_at_nodes`,
:func:`sassi.core.flexibility.flexibility_matrix`), so the whole chain
deck -> module -> file -> API is verified.

* VP-02b  SITE 1-D column vs the closed-form layer-on-half-space TFs (gates D-SIT-02)
* VP-05   Rayleigh phase velocity of a deep uniform stratum
* VP-06   Love-wave dispersion of a layer over a half-space
* VP-07   static Boussinesq/Cerruti point-load solutions and reciprocity of F_ff
* VP-08   dynamic surface Green functions of a half-space (Wong 1975 / Apsel 1979)
* VP-09   POINT3 far field vs the exact thin-layer point-load series (Kausel 1981)
* VP-43   (P1) thin-layer convergence of the surface vertical compliance
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from ...conventions import cfactor
from ...core import greens_tlm, tlm
from ...core.flexibility import flexibility_matrix, frequency_row
from ...core.freefield import free_field_at_nodes
from ...io import decks
from ...io.container import Container
from ...io.files import read_container
from ...modules.base import run_module
from .. import VPResult, problem, worse

GRAV = 9.81


# ---------------------------------------------------------------------------------------
# helpers: decks -> modules -> files
# ---------------------------------------------------------------------------------------
def _run(module: str, workdir: Path, model: str) -> None:
    rc = run_module(module, model, workdir)
    if rc != 0:
        listing = (workdir / f"{model}_{module}.out").read_text(encoding="utf-8")
        raise RuntimeError(f"{module} failed (status {rc}):\n{listing[-3000:]}")


def run_site(workdir: Path, model: str, layers: Sequence[Sequence[float]], halfspace: Sequence[float],
             fnum: Sequence[int], df: float, nl: int, mode2: bool = True, cl: int = 1, cm: int = 0,
             waves: Sequence[Sequence[float]] = ((2, 1, 1.0, 1.0, 0.0),), hslaw: Optional[str] = None,
             wopt: int = 0) -> Dict[str, Container]:
    """Write a SITE deck (layers: (thick, rho, vs, vp, beta_s, beta_p); unit weight = rho g)
    and run SITE.  ``hslaw=None`` keeps the deck default (D-SIT-02)."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    d = decks.new("SITE")
    d["title"] = f"{model} verification site"
    d["model"] = model
    d["gravity"] = GRAV
    d["nl"] = nl
    d["df"] = df
    d["mode2"] = int(mode2)
    d["cl"] = cl
    d["cm"] = cm
    d["wopt"] = wopt
    if hslaw is not None:
        d["hslaw"] = hslaw
    for i, (h, rho, vs, vp, bs, bp) in enumerate(layers):
        d.table("layers").append([i + 1, h, rho * GRAV, vp, vs, bp, bs])
    h, rho, vs, vp, bs, bp = halfspace
    d.table("halfspace").append([len(layers) + 1, 0.0, rho * GRAV, vp, vs, bp, bs])
    for w in waves:
        d.table("waves").append(list(w))
    for n in fnum:
        d.table("freqs").append([int(n)])
    decks.write(workdir / f"{model}.sit", d)
    _run("SITE", workdir, model)
    out = {"FILE2": read_container(workdir / "FILE2", "FILE2")}
    if mode2:
        out["FILE1"] = read_container(workdir / "FILE1", "FILE1")
    return out


def run_point(workdir: Path, model: str, layer: int, rad: float, df: float, fnum: Sequence[int] = (),
              dim: int = 2) -> Container:
    d = decks.new("POINT")
    d["title"] = f"{model} verification point loads"
    d["model"] = model
    d["layer"] = layer
    d["rad"] = rad
    d["dim"] = dim
    d["df"] = df
    for n in fnum:
        d.table("freqs").append([int(n)])
    decks.write(Path(workdir) / f"{model}.poi", d)
    _run("POINT", Path(workdir), model)
    return read_container(Path(workdir) / "FILE3", "FILE3")


def _vp(vs: float, nu: float) -> float:
    return vs * np.sqrt(2.0 * (1.0 - nu) / (1.0 - 2.0 * nu))


# ---------------------------------------------------------------------------------------
# VP-02b
# ---------------------------------------------------------------------------------------
VP02_F = np.array([0.5, 1.0, 1.5, 5.0 / 3.0, 2.0, 3.0, 5.0, 25.0 / 3.0])
VP02_WITHIN = np.array([1.1216, 1.6931, 5.7710, 12.7152, 3.1156, 1.0411, 4.2038, 2.4815])   # R2 A.2
VP02_OUTCROP = np.array([1.1150, 1.6207, 3.4116, 3.8327, 2.4115, 1.0087, 2.3554, 1.6702])


def _vp02_exact(f):
    """Closed form of R2 A.1 (SASSI damping form): |1/cos k*H| and |1/(cos k*H + i a* sin k*H)|."""
    H, vs, rho, xi, vr, rr, xr = 30.0, 200.0, 2.0, 0.05, 1000.0, 2.2, 0.01
    w = 2 * np.pi * np.asarray(f)
    vss = vs * np.sqrt(cfactor(xi))
    vrs = vr * np.sqrt(cfactor(xr))
    k = w / vss
    a = rho * vss / (rr * vrs)
    return np.abs(1 / np.cos(k * H)), np.abs(1 / (np.cos(k * H) + 1j * a * np.sin(k * H)))


def _vp02_site(workdir: Path, nsub: int, nl: int, law: str):
    """Uniform 30 m layer (Vs 200, rho 2.0, 5 %) on rock (Vr 1000, rho 2.2, 1 %): SITE with the
    control point at the top of the half-space (within motion), vertical SV."""
    H, vs, rho, xi = 30.0, 200.0, 2.0, 0.05
    layers = [(H / nsub, rho, vs, 2 * vs, xi, xi)] * nsub
    hs = (0.0, 2.2, 1000.0, 2000.0, 0.01, 0.01)
    df = 1.0 / 6.0
    fnum = np.round(VP02_F / df).astype(int)
    files = run_site(workdir, "vp02b", layers, hs, fnum, df, nl, cl=nsub + 1, cm=0, hslaw=law)
    f1 = files["FILE1"]
    within = np.array([abs(free_field_at_nodes(f1, q, [[0.0, 0.0, 0.0]], [1], 0.0, 0.0, 0.0)[0, 0])
                       for q in range(len(fnum))])
    # outcrop TF: raw surface motion / discrete outcrop motion at the half-space top (R1 §2.7a)
    outcrop = np.abs(f1["U"][0, :, 0, 0] * f1["x_ucp"][0] / f1["x_outcrop"][0])
    uniform_used = files["FILE2"]["x_hs_uniform"]
    return within, outcrop, uniform_used


@problem("VP-02b", "SITE uniform layer on a half-space vs closed form (TLM + variable-depth half-space)",
         tier="P0", modules=["SITE"], source="R2 A.2, R1 §2.4/V4, D-SIT-02")
def vp02b(workdir: Path) -> VPResult:
    r = VPResult()
    ex_w, ex_o = _vp02_exact(VP02_F)
    for f, a, b in zip(VP02_F, ex_w, VP02_WITHIN):
        r.check(f"closed form |surf/base-within| at {f:.4g} Hz (R2 A.2 value)", a, b, rtol=5e-5)
    for f, a, b in zip(VP02_F, ex_o, VP02_OUTCROP):
        r.check(f"closed form |surf/outcrop| at {f:.4g} Hz (R2 A.2 value)", a, b, rtol=5e-5)
    default_law = str(decks.new("SITE")["hslaw"]).lower()
    laws = ["uniform"] + ([default_law] if default_law != "uniform" else [])
    results = {}
    for law in laws:
        for nl in (20, 10):
            results[(law, nl)] = _vp02_site(Path(workdir) / f"{law}_{nl}", 25, nl, law)
    w20, o20, _ = results[("uniform", 20)]
    # within TF (independent of the half-space model, R1 §2.7a): h = 1.2 m <= lambda/20 at 8.33 Hz
    for f, a, b in zip(VP02_F, w20, ex_w):
        r.check(f"SITE within TF |surf/HS-top| at {f:.4g} Hz (h = lambda/20, 1 % provisional)", a, b, rtol=0.01)
    wc, _, _ = _vp02_site(Path(workdir) / "coarse", 13, 20, "uniform")
    e_fine = np.max(np.abs(w20 - ex_w) / ex_w)
    e_coarse = np.max(np.abs(wc - ex_w) / ex_w)
    r.require("within TF converges: max error at h = lambda/20 below max error at h = lambda/10.4",
              e_fine < e_coarse, note=f"{e_fine:.2e} vs {e_coarse:.2e}")
    tol = {20: 0.01, 10: 0.02}
    gate = {}
    for (law, nl), (_, out, uni) in results.items():
        err = np.abs(out - ex_o) / ex_o
        gate[(law, nl)] = float(err.max())
        for f, a, b in zip(VP02_F, out, ex_o):
            r.check(f"outcrop TF |surf/HS outcrop| at {f:.4g} Hz, nl = {nl}, law {law}", a, b, rtol=tol[nl],
                    note=f"uniform fallback used at {int(uni.sum())}/{uni.size} frequencies")
    for (law, nl), e in gate.items():
        r.notes.append(f"D-SIT-02 gate: law '{law}', nl = {nl}: max outcrop error {100 * e:.3f} % "
                       f"(tolerance {100 * tol[nl]:.0f} %) -> {'PASS' if e <= tol[nl] else 'FAIL'}")
    if default_law != "uniform" and any(gate[(default_law, nl)] > tol[nl] for nl in (20, 10)):
        r.notes.append(f"D-SIT-02: the deck default law '{default_law}' fails VP-02b; the decision makes "
                       "'uniform' the default -- the SITE deck schema default (sassi/io/decks.py, lead) "
                       "and EDUOPT,HSLAW must be switched to UNIFORM.")
    r.notes.append(f"within TF max error {100 * e_fine:.4f} % at h = lambda/20 (mixed mass)")
    return r


# ---------------------------------------------------------------------------------------
# VP-05
# ---------------------------------------------------------------------------------------
VP05_REF = {0.25: 0.919402, 1.0 / 3.0: 0.932526, 0.45: 0.948960, 0.49: 0.954074}   # R2 D.5


@problem("VP-05", "Rayleigh phase velocity of a deep uniform stratum (SITE Mode 1)", tier="P0",
         modules=["SITE"], source="R2 D.5, R1 V2")
def vp05(workdir: Path) -> VPResult:
    """Uniform stratum Vs = 100 m/s, 30 m deep on a rigid base, sublayers h = 0.25 m; at
    10 / 20 / 40 Hz this is h = lambda/40, lambda/20, lambda/10 (depth >= 3 lambda).  The
    fundamental Rayleigh mode is the propagating mode with the largest Re k."""
    r = VPResult()
    vs, rho, h, depth = 100.0, 2.0, 0.25, 30.0
    n = int(round(depth / h))
    fnum, df = [1, 2, 4], 10.0
    for nu, ref in VP05_REF.items():
        vp = _vp(vs, nu)
        files = run_site(Path(workdir) / f"nu{nu:.3f}", "vp05", [(h, rho, vs, vp, 0.0, 0.0)] * n,
                         (0.0, rho, vs, vp, 0.0, 0.0), fnum, df, nl=0, mode2=False)
        f2 = files["FILE2"]
        err = {}
        for q, N in zip(range(3), (40, 20, 10)):
            k = f2["kR"][q]
            j = tlm.select_mode(k, 1)
            c = 2 * np.pi * f2["freq"][q] / k[j].real / vs
            err[N] = (c - ref) / ref
            # Lead decision (wave-1 review): linear thin layers lock volumetrically for nearly
            # incompressible soil, so for nu >= 0.45 the 0.5 % criterion applies at 40 sublayers per
            # lambda (SASSI-EDU warns G-09/EDU-11 for nu > 0.47); nu <= 1/3 keeps 10 and 20 per lambda.
            if (nu < 0.4 and N in (20, 10)) or (nu >= 0.4 and N == 40):
                r.check(f"c_R/Vs, nu = {nu:.4g}, {N} sublayers per lambda", c, ref, rtol=0.005)
        p = np.log2(abs(err[20]) / abs(err[40]))
        if nu < 0.4:
            r.check(f"observed convergence order (lambda/20 -> lambda/40), nu = {nu:.4g}", p, 2.0, rtol=0.1)
        else:
            r.require(f"monotone convergence 10 -> 20 -> 40 sublayers per lambda, nu = {nu:.4g}",
                      abs(err[10]) > abs(err[20]) > abs(err[40]))
        r.notes.append(f"nu = {nu:.4g}: errors {100 * err[10]:+.3f} % (10/lambda), {100 * err[20]:+.3f} % "
                       f"(20/lambda), {100 * err[40]:+.3f} % (40/lambda); order {p:.2f}")
    r.notes.append("The thin-layer matrices are linear in z (Kausel TLM, D-SIT-01): for nearly incompressible "
                   "soil (nu >= 0.45) they lock volumetrically and the error is stiffness-driven (it is the "
                   "same with lumped mass). nu = 0.49 needs about 40 sublayers per lambda for 0.5 %.")
    return r


# ---------------------------------------------------------------------------------------
# VP-06
# ---------------------------------------------------------------------------------------
VP06_C0 = {2: 387.536, 5: 300.026, 10: 224.175, 20: 205.949}     # R2 D.5
VP06_C1_20 = 279.222
VP06_CUTOFF = (0.0, 11.547, 23.094)                               # R2 D.5: f_n = n b1 / (2 H sqrt(1 - b1^2/b2^2))
VP06_DF = 0.05                                                     # frequency step of the SITE runs (Hz)
VP06_SWEEPS = ((0.05, 1.5), (10.5, 13.5), (22.0, 26.0))            # cut-off sweeps (Hz) of modes 0, 1, 2


def _love_speeds(f2: Container, q: int, beta2: float):
    """Phase velocities ``w / Re k`` of the Love modes with |Im k| <= 0.05 Re k at row q, split into
    trapped (c < beta2) and leaky (c >= beta2) modes, slowest first."""
    k = f2["kL"][q]
    w = 2 * np.pi * f2["freq"][q]
    ok = (k.real > 0) & (np.abs(k.imag) <= 0.05 * k.real)
    c = np.sort(w / k[ok].real)
    return c[c < beta2], c[c >= beta2]


def _trapped_love(f2: Container, q: int, beta2: float) -> np.ndarray:
    """Phase velocities of the trapped Love modes (c < beta2), slowest first."""
    return _love_speeds(f2, q, beta2)[0]


def _love_cutoff(f2: Container, n: int, beta2: float, f_lo: float, f_hi: float):
    """Discrete cut-off frequency of Love mode ``n`` (0 = fundamental) on the sweep [f_lo, f_hi]:
    the frequency where its phase velocity crosses beta2, i.e. where the number of trapped modes
    grows from n to n + 1.  Between the two bracketing sweep frequencies the velocity is
    interpolated linearly (at the lower one the mode is the slowest leaky mode, c just above
    beta2).  Returns ``(f_c, note)``; ``f_c = f_lo`` when the mode is already trapped at f_lo
    (the cut-off lies below the sweep), NaN when it is never trapped on the sweep."""
    freq = np.asarray(f2["freq"])
    rows = np.flatnonzero((freq >= f_lo - 1e-9) & (freq <= f_hi + 1e-9))
    first_t, _ = _love_speeds(f2, int(rows[0]), beta2)
    if first_t.size > n:
        return float(freq[rows[0]]), f"mode {n} already trapped at {freq[rows[0]]:.4g} Hz (cut-off below the sweep)"
    for a, b in zip(rows[:-1], rows[1:]):
        ta, la = _love_speeds(f2, int(a), beta2)
        tb, _ = _love_speeds(f2, int(b), beta2)
        if ta.size <= n < tb.size:
            fa, fb, cb = float(freq[a]), float(freq[b]), float(tb[n])
            if la.size == 0:
                return 0.5 * (fa + fb), f"bracketed by {fa:.4g} and {fb:.4g} Hz"
            ca = float(la[0])
            return fa + (ca - beta2) * (fb - fa) / (ca - cb), f"c = {ca:.2f} at {fa:.4g} Hz, {cb:.2f} at {fb:.4g} Hz"
    return float("nan"), f"mode {n} not trapped anywhere on {f_lo:.4g}..{f_hi:.4g} Hz"


@problem("VP-06", "Love-wave dispersion of a layer over a half-space (SITE Mode 1, SH)", tier="P0",
         modules=["SITE"], source="R2 D.5")
def vp06(workdir: Path) -> VPResult:
    """H = 10 m, beta1 = 200, beta2 = 400 m/s, equal density, elastic; the half-space is SITE's
    variable-depth half-space (20 generated sublayers over 1.5 lambda + dashpots, D-SIT-02).  Layer
    split in 20 sublayers (lambda/20 at 20 Hz); 10 and 40 sublayers for the convergence study.

    Reference criteria (requirements §6.3 as amended by lead decision D-W1-03): cut-offs 0, 11.547 and
    23.094 Hz; fundamental and first higher mode phase velocities within 0.5 %.

    * cut-offs: SITE is run on fine frequency sweeps (df = 0.05 Hz) around each cut-off; the
      discrete cut-off of mode n is where its phase velocity crosses beta2.  The variable-depth
      half-space cannot represent the grazing field near a Love cut-off (D-W1-03 a), so the cut-offs
      are checked with the provisional 10 % (for the zero cut-off of the fundamental mode 10 % of the
      first non-zero cut-off, a relative tolerance on 0 being undefined);
    * convergence (D-W1-03 b, mixed mass of D-CNV-05): the 20 Hz phase velocities converge
      monotonically with the sublayer count (10, 20, 40), not necessarily from above.
    """
    r = VPResult()
    b1, b2, H, rho = 200.0, 400.0, 10.0, 2.0
    fref = {fn: int(round(fn / VP06_DF)) for fn in VP06_C0}
    sweep = sorted({int(round(f / VP06_DF)) for lo, hi in VP06_SWEEPS
                    for f in np.arange(lo, hi + 0.5 * VP06_DF, VP06_DF)})
    c0, c1 = {}, {}
    f2_20 = None
    for nsub in (10, 20, 40):
        fnum = sorted(set(fref.values()) | (set(sweep) if nsub == 20 else set()))
        files = run_site(Path(workdir) / f"n{nsub}", "vp06", [(H / nsub, rho, b1, 2 * b1, 0.0, 0.0)] * nsub,
                         (0.0, rho, b2, 2 * b2, 0.0, 0.0), fnum, VP06_DF, nl=20, mode2=False)
        f2 = files["FILE2"]
        cs = {fn: _trapped_love(f2, frequency_row(f2, n), b2) for fn, n in fref.items()}
        c0[nsub] = {fn: float(cs[fn][0]) for fn in cs}
        c1[nsub] = float(cs[20][1]) if cs[20].size > 1 else float("nan")
        if nsub == 20:
            f2_20 = f2
            for fn in VP06_C0:
                r.check(f"fundamental Love c at {fn} Hz", c0[20][fn], VP06_C0[fn], rtol=0.005)
            r.check("first higher Love mode c at 20 Hz", c1[20], VP06_C1_20, rtol=0.005)
            r.require("one trapped mode at 2, 5 and 10 Hz (below the 11.547 Hz cut-off)",
                      all(cs[fn].size == 1 for fn in (2, 5, 10)))
            r.require("two trapped modes at 20 Hz (between the 11.547 and 23.094 Hz cut-offs)", cs[20].size == 2)
    # cut-off frequencies (requirements §6.3: 0, 11.547, 23.094 Hz, 0.5 %)
    for n, (fex, (lo, hi)) in enumerate(zip(VP06_CUTOFF, VP06_SWEEPS)):
        fc, how = _love_cutoff(f2_20, n, b2, lo, hi)
        # Lead decision (wave-1 review): near a cut-off the mode penetrates deeply into the half-space;
        # the variable-depth half-space (1.5 lambda + dashpots, R1 §2.4) shifts the discrete cut-offs up,
        # a known limitation of the method.  Cut-offs are checked with a provisional 10 % tolerance.
        if fex == 0.0:
            r.check(f"cut-off of Love mode {n} (exact 0 Hz; tolerance 10 % of the 11.547 Hz cut-off, provisional)",
                    fc, fex, atol=0.10 * VP06_CUTOFF[1], note=how)
        else:
            r.check(f"cut-off of Love mode {n} (exact {fex:.5g} Hz; provisional 10 %)", fc, fex, rtol=0.10, note=how)
        r.notes.append(f"Love mode {n}: discrete cut-off {fc:.4f} Hz vs exact {fex:.5g} Hz ({how})")
    # convergence from above (R2 D.5)
    series = [(f"fundamental c at {fn} Hz", [c0[n][fn] for n in (10, 20, 40)], VP06_C0[fn]) for fn in VP06_C0]
    series.append(("first higher mode c at 20 Hz", [c1[n] for n in (10, 20, 40)], VP06_C1_20))
    # Lead decision: R2 D.5's 'convergence from above' holds for consistent mass; the mixed mass of
    # D-CNV-05 converges from below.  Required: monotone convergence of the layer-discretisation-driven
    # cases (20 Hz, where the error exceeds the half-space-model error).
    for name, (ca, cb, cc), ref in series:
        if "20 Hz" not in name:
            continue
        r.require(f"monotone convergence: {name}, |c - exact| decreasing for 10/20/40 sublayers",
                  abs(ca - ref) > abs(cb - ref) > abs(cc - ref),
                  note=f"{ca - ref:+.4f} / {cb - ref:+.4f} / {cc - ref:+.4f} m/s")
    sign = "below" if c1[40] < VP06_C1_20 else "above"
    r.notes.append(f"mode-1 c at 20 Hz: {c1[10]:.3f}, {c1[20]:.3f}, {c1[40]:.3f} m/s (10/20/40 sublayers): "
                   f"converges from {sign}. R2 D.5 states convergence from above for the thin-layer method; "
                   "that holds for consistent mass (tests/unit/test_tlm.py), the mixed mass of D-CNV-05 "
                   "converges from below.  At 2-10 Hz the fundamental-mode error (< 0.05 %) is dominated by the "
                   "variable-depth half-space, not by the layer discretisation, so it is not monotone in N.")
    r.notes.append("Cut-offs: near a cut-off the mode penetrates deeply into the half-space (decay exponent "
                   "-> 0); the 1.5-lambda buffer with Lysmer-Kuhlemeyer dashpots (R1 §2.4) cannot represent this "
                   "grazing field, so the discrete modes stay leaky above the exact cut-offs.  The error falls "
                   "only slowly when half-space material is added as user layers above the buffer (mode 1: "
                   "+2.5 % with 100 m, +1.1 % with 300 m added; scratch study during the review fix).")
    return r


# ---------------------------------------------------------------------------------------
# VP-07
# ---------------------------------------------------------------------------------------
def _graded_layers(h0: float, growth: float, hmax: float, depth: float) -> List[float]:
    hs, z, h = [], 0.0, h0
    while z < depth:
        hs.append(h)
        z += h
        h = min(h * growth, hmax)
    return hs


def _boussinesq_cerruti(x: float, y: float, mu: float, nu: float) -> np.ndarray:
    """Surface displacements (z up) at (x, y) due to unit surface forces at the origin (R2 D.1):
    columns P_x, P_y, P_z (upward).  Cerruti for horizontal forces, Boussinesq for the vertical."""
    r = np.hypot(x, y)
    B = 1.0 / (2 * np.pi * mu * r)
    C = (1 - 2 * nu) / (4 * np.pi * mu * r * r)
    return np.array([[B * ((1 - nu) + nu * x * x / r ** 2), B * nu * x * y / r ** 2, C * x],
                     [B * nu * x * y / r ** 2, B * ((1 - nu) + nu * y * y / r ** 2), C * y],
                     [-C * x, -C * y, B * (1 - nu)]])


@problem("VP-07", "Static point-load Green functions (Boussinesq/Cerruti) and reciprocity of F_ff",
         tier="P0", modules=["SITE", "POINT", "ANALYS"], source="R2 D.1, R1 V3/V10")
def vp07(workdir: Path) -> VPResult:
    """Quasi-static (0.01 Hz) deep graded stratum (115 sublayers, 0.05 m growing by 6 % up to 10 m,
    406 m on a rigid base), Vs = 100 m/s, nu = 0.3, beta = 0.5 %; POINT3 with R0 = 0.1 m; the
    flexibility matrix of surface nodes at r = 1-5 m from the origin is compared component by
    component with Boussinesq/Cerruti (real parts)."""
    r = VPResult()
    vs, nu, rho, beta = 100.0, 0.3, 2.0, 0.005
    vp = _vp(vs, nu)
    mu = rho * vs * vs
    hs = _graded_layers(0.05, 1.06, 10.0, 400.0)
    layers = [(h, rho, vs, vp, beta, beta) for h in hs]
    run_site(workdir, "vp07", layers, (0.0, rho, vs, vp, beta, beta), [1], 0.01, nl=0, mode2=False)
    f3 = run_point(workdir, "vp07", layer=0, rad=0.1, df=0.01, fnum=[1])
    pts = np.array([[0, 0], [1, 0], [2, 0], [5, 0], [0, 3], [3 / np.sqrt(2), 3 / np.sqrt(2)], [-4, 1.5]], float)
    F = flexibility_matrix(f3, 0, pts, np.ones(len(pts), int), symmetrize=False)
    names = ("x", "y", "z")
    for i in range(1, len(pts)):
        x, y = pts[i]
        ref = _boussinesq_cerruti(x, y, mu, nu)
        got = F[3 * i:3 * i + 3, 0:3].real
        big = np.abs(ref) > 1e-3 * np.abs(ref).max()
        for a in range(3):
            for b in range(3):
                if big[a, b]:
                    r.check(f"u_{names[a]} due to P_{names[b]} at ({x:.3g}, {y:.3g}), r = {np.hypot(x, y):.3g} m",
                            got[a, b], ref[a, b], rtol=0.02)
    asym = np.linalg.norm(F - F.T) / np.linalg.norm(F)
    r.check("reciprocity ||F - F^T|| / ||F|| before symmetrisation", asym, 0.0, atol=1e-3)
    Fs = flexibility_matrix(f3, 0, pts, np.ones(len(pts), int))
    r.require("symmetrised F_ff is exactly symmetric", np.array_equal(Fs, Fs.T))
    return r


# ---------------------------------------------------------------------------------------
# VP-08
# ---------------------------------------------------------------------------------------
# Wong (1975) via Apsel (1979) Table 5.1, nu = 0.33: mu r u / P (Re, Im) -- R2 D.3
VP08_TABLE = np.array([
    [0.0, -.027, .000, .106, .000, .159, .000, -.106, .000],
    [0.5, -.032, .007, .087, -.062, .146, -.058, -.089, .058],
    [1.0, -.033, .025, .037, -.102, .112, -.105, -.045, .099],
    [1.5, -.020, .047, -.029, -.108, .063, -.133, .015, .112],
    [2.0, .006, .060, -.087, -.077, .011, -.137, .075, .093],
    [2.5, .041, .058, -.120, -.017, -.034, -.120, .118, .045],
    [3.0, .074, .035, -.114, .053, -.064, -.090, .132, -.020],
    [3.5, .092, -.005, -.070, .110, -.076, -.054, .111, -.084],
    [4.0, .087, -.054, -.001, .134, -.075, -.024, .059, -.132],
    [4.5, .056, -.096, .072, .118, -.065, -.004, -.014, -.144],
    [5.0, .005, -.120, .127, .064, -.054, .004, -.088, -.129],
    [5.5, -.054, -.116, .144, -.013, -.048, .003, -.145, -.077]])
VP08_APSEL_R3 = {"U_r1": (-.068, -.090), "U_theta1": (.128, -.020)}    # Apsel's own integration at r0 = 3


def _wong_components(F: np.ndarray, i: int, r: float, mu: float) -> Dict[str, complex]:
    """Wong's normalised components from the F block of node i (on the +x axis, load at node 0):
    U_r0 = -p mu r (radial, downward vertical force), U_z0 = q mu r, U_r1 = u mu r, U_theta1 = -v mu r."""
    s = mu * r
    return {"U_r0": -F[3 * i, 2] * s, "U_z0": F[3 * i + 2, 2] * s, "U_r1": F[3 * i, 0] * s,
            "U_theta1": -F[3 * i + 1, 1] * s}


@problem("VP-08", "Dynamic surface Green functions of a half-space (Wong 1975)", tier="P0",
         modules=["SITE", "POINT"], source="R2 D.3")
def vp08(workdir: Path) -> VPResult:
    """Near-elastic half-space (beta = 1e-4), Vs = 100 m/s, nu = 0.33: 119 graded sublayers (0.05 m
    growing by 8 % to lambda/20 = 0.5 m, 50 m deep) + 20 generated half-space sublayers with base
    dashpots; POINT3 with R0 = 0.1 m.  f = 10 Hz gives r0 = w r/Vs = 0.5..5.5 at r = 0.8..8.8 m;
    the static row (r0 = 0) is taken at 0.1 Hz and r = 1 m (r0 = 0.006).  Tolerance +-0.005."""
    r = VPResult()
    vs, nu, rho, beta = 100.0, 0.33, 2.0, 1e-4
    vp = _vp(vs, nu)
    mu = rho * vs * vs
    f = 10.0
    df = f / 100.0
    hs = _graded_layers(0.05, 1.08, vs / f / 20.0, 50.0)
    layers = [(h, rho, vs, vp, beta, beta) for h in hs]
    run_site(workdir, "vp08", layers, (0.0, rho, vs, vp, beta, beta), [1, 100], df, nl=20, mode2=False)
    f3 = run_point(workdir, "vp08", layer=0, rad=0.1, df=df, fnum=[1, 100])
    keys = ("U_r0", "U_z0", "U_r1", "U_theta1")
    # static row
    qs = frequency_row(f3, 1)
    rs = 1.0
    F0 = flexibility_matrix(f3, qs, np.array([[0.0, 0.0], [rs, 0.0]]), [1, 1])
    c0 = _wong_components(F0, 1, rs, mu)
    for j, key in enumerate(keys):
        r.check(f"{key} Re at r0 = 0 (static limit)", c0[key].real, VP08_TABLE[0, 1 + 2 * j], atol=0.005)
        r.check(f"{key} Im at r0 = 0 (static limit)", c0[key].imag, VP08_TABLE[0, 2 + 2 * j], atol=0.005)
    # dynamic rows
    qd = frequency_row(f3, 100)
    kk = 2 * np.pi * f / vs
    rr = VP08_TABLE[1:, 0] / kk
    pts = np.concatenate([[[0.0, 0.0]], np.c_[rr, np.zeros_like(rr)]])
    F = flexibility_matrix(f3, qd, pts, np.ones(len(pts), int))
    worst = 0.0
    for i, row in enumerate(VP08_TABLE[1:], start=1):
        c = _wong_components(F, i, rr[i - 1], mu)
        for j, key in enumerate(keys):
            for part, val in (("Re", c[key].real), ("Im", c[key].imag)):
                ref = row[1 + 2 * j + (part == "Im")]
                r.check(f"{key} {part} at r0 = {row[0]:.1f}", val, ref, atol=0.005)
                worst = worse(worst, abs(val - ref))
        if row[0] == 3.0:
            for key, (re, im) in VP08_APSEL_R3.items():
                r.check(f"{key} Re at r0 = 3.0 vs Apsel wavenumber integration", c[key].real, re, atol=0.005)
                r.check(f"{key} Im at r0 = 3.0 vs Apsel wavenumber integration", c[key].imag, im, atol=0.005)
    r.notes.append(f"max |error| vs Wong's table over r0 = 0.5..5.5: {worst:.4f} (Wong and Apsel themselves "
                   "differ by up to 0.004, e.g. U_r1 at r0 = 3.0)")
    return r


# ---------------------------------------------------------------------------------------
# VP-09
# ---------------------------------------------------------------------------------------
VP09_LAYERS = ([(0.5, 1.9, 150.0, 300.0, 0.05, 0.05)] * 8 + [(1.0, 2.0, 250.0, 500.0, 0.04, 0.04)] * 10
               + [(2.0, 2.1, 400.0, 800.0, 0.03, 0.03)] * 10)


# distance (m) -> (tolerance, dominant components only?)  -- requirements §6.3 / R1 §3.3 (v), V8
VP09_TOL = {3.0: (0.0235, False), 6.0: (0.008, True), 12.0: (0.0025, False), 25.0: (0.0025, False)}


@problem("VP-09", "POINT3 far field vs the exact thin-layer point-load series", tier="P0",
         modules=["SITE", "POINT"], source="R1 §3.3 (v), V8")
def vp09(workdir: Path) -> VPResult:
    """R1 V8 setup: 28 layers on a rigid base, f = 5 Hz, R0 = 0.9 m; loads at interfaces 1 and 7,
    observations at interfaces 1 and 5, r = 3.0, 6.0, 12.0, 25.0 m (3.3, 6.7, 13.3, 27.8 R0).
    The reference is Kausel's point-load series for the *same* discrete medium (FILE2 modes), so the
    comparison isolates the central-zone approximation.

    Tolerances (requirements §6.3; 3.3 R0 restated as 2.35 % by lead decision: R1 V8's own prototype
    gives 2.310 %, which R1 had rounded to 2.3 %): <= 2.35 % at 3.3 R0 and <= 0.25 % beyond 13 R0 on *all*
    components; <= 0.8 % at 6.7 R0 on the dominant components only (|value| >= 10 % of the largest
    component of the same load, observation and harmonic) -- R1 V8 exempts at that distance two
    components that are 10-100 times below the dominant one."""
    r = VPResult()
    R0 = 0.9
    files = run_site(workdir, "vp09", VP09_LAYERS, (0.0, 2.1, 400.0, 800.0, 0.03, 0.03), [1], 5.0, nl=0,
                     mode2=False)
    f2 = files["FILE2"]
    f3 = run_point(workdir, "vp09", layer=6, rad=R0, df=5.0, fnum=[1])
    md = tlm.modes_from_file2(f2, 0)
    loads, obs = (1, 7), (1, 5)
    xy = [[0.0, 0.0], [0.0, 0.0]]
    ifc = list(loads)
    index = {}
    for rho_ in VP09_TOL:
        for m in obs:
            index[(rho_, m)] = len(xy)
            xy.append([rho_, 0.0])
            ifc.append(m)
    F = flexibility_matrix(f3, 0, np.asarray(xy), np.asarray(ifc), symmetrize=False)
    allmax = {}
    for rho_, (tol, dominant_only) in VP09_TOL.items():
        worst_all = (0.0, "")
        for jl, n in enumerate(loads):
            for m in obs:
                i = index[(rho_, m)]
                blk = F[3 * i:3 * i + 3, 3 * jl:3 * jl + 3]
                got = {"u": blk[0, 0], "v": blk[1, 1], "w": blk[2, 0], "p": blk[0, 2], "q": blk[2, 2]}
                g = greens_tlm.point_load(md, m - 1, n - 1, rho_)
                for group in (("u", "v", "w"), ("p", "q")):
                    gmax = max(abs(complex(g[c])) for c in group)
                    for c in group:
                        ref = complex(g[c])
                        e = abs(got[c] - ref) / abs(ref)
                        where = f"{c}~, load iface {n}, obs iface {m}"
                        if e > worst_all[0]:
                            worst_all = (e, where)
                        minor = abs(ref) < 0.1 * gmax
                        if dominant_only and minor:
                            continue
                        r.check(f"|{c}~ - exact|/|exact|, load iface {n}, obs iface {m}, r = {rho_ / R0:.1f} R0"
                                + (" (minor component)" if minor else ""), e, 0.0, atol=tol)
        allmax[rho_] = worst_all
    r.notes.append("all-component maxima: " + ", ".join(f"{rho_ / R0:.1f} R0: {100 * v[0]:.3f} % ({v[1]})"
                                                      for rho_, v in allmax.items()))
    e33 = allmax[3.0][0]
    if e33 > VP09_TOL[3.0][0]:
        r.notes.append(f"3.3 R0: {100 * e33:.3f} % ({allmax[3.0][1]}) exceeds the tolerance "
                       f"{100 * VP09_TOL[3.0][0]:.2f} % of lead decision D-W1-04 (R1 V8's own prototype, "
                       "docs/spec/R1_checks/point3.py, gives 2.310 % for the same component).")
    return r


# ---------------------------------------------------------------------------------------
# VP-43 (P1): monotone thin-layer convergence (Kausel-Peek disk loads)
# ---------------------------------------------------------------------------------------
@problem("VP-43", "Thin-layer monotone convergence of disk-load compliances (P1, Lamb pulse not included)",
         tier="P1", modules=["SITE"], source="R2 D.5 (Kausel-Peek 1981), R1 §3.2")
def vp43(workdir: Path) -> VPResult:
    """Homogeneous stratum H = 10 m, Vs = 100 m/s, nu = 1/4, beta = 5 %, rigid base, discretised
    with N = 4, 6, 12, 24 sublayers (N = 96 as the converged reference); average displacement of a
    uniform surface disk load of radius 1 m from the SITE Mode 1 modes (FILE2) with Kausel's
    closed-form disk integrals, at 2 Hz.  The discretisation stiffens the system: the compliance
    must increase monotonically with N towards the converged value (R2 D.5).  The Lamb-pulse part of
    VP-43 (inverse FFT of the surface vertical Green function) is not implemented."""
    r = VPResult()
    vs, nu, rho, beta, H, R = 100.0, 0.25, 2.0, 0.05, 10.0, 1.0
    vp = _vp(vs, nu)
    comp = {}
    for N in (4, 6, 12, 24, 96):
        files = run_site(Path(workdir) / f"N{N}", "vp43", [(H / N, rho, vs, vp, beta, beta)] * N,
                         (0.0, rho, vs, vp, beta, beta), [4], 0.5, nl=0, mode2=False)
        md = tlm.modes_from_file2(files["FILE2"], 0)
        d = greens_tlm.disk_load_average(md, 0, 0, R)
        comp[N] = (abs(d["vertical"]), abs(d["horizontal"]))
    for key, label in ((0, "vertical"), (1, "horizontal")):
        seq = [comp[N][key] for N in (4, 6, 12, 24, 96)]
        r.require(f"{label} disk compliance increases monotonically with N = 4, 6, 12, 24 towards N = 96",
                  all(a < b for a, b in zip(seq, seq[1:])),
                  note=", ".join(f"{v / seq[-1]:.4f}" for v in seq) + " (ratio to N = 96)")
    r.notes.append("Lamb's problem check values (time domain) are not evaluated in this build.")
    return r
