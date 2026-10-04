"""Verification problems of the non-linear (near-field) soil SSI iterations (requirements 4.4 item 7,
4.10 item 5, 2.4; decisions D-NLS-01 ... D-NLS-04; :mod:`sassi.core.nlsoil`).

* VP-N1  consistency with SOIL (SHAKE): a laterally uniform column of near-field soil, modelled with
  nonlinear SOLID elements (structure, ETYPE 1, low-strain material = G_max) inside a flexible-volume
  (FV) excavation of a uniform sand deposit on rock, driven by vertically propagating SV waves.  The
  free field (SITE) and the excavated soil carry the strain-compatible properties of a SOIL (SHAKE)
  run of the same column, motion and curves; the SSI control motion is the SOIL surface motion.
  Without a structure there is no secondary nonlinearity, so

  (a) the near-field iterations started from the LOW-STRAIN properties -- one material per layer
      with GFAC = G_max / G_layer and DFAC = beta_0 / beta_layer, because the .pin factors scale the
      free field (requirements 4.4 item 7) -- must converge, monotonically and within the iteration
      limit of D-NLS-02, to the SHAKE strain-compatible G/G_max and damping profile: 5 % on G, 10 % on
      damping per layer (the differences are the discretisation of the column: thin-layer and 8-node
      elements of 1 m against the exact continuous layers of SHAKE, and the interpolation of the
      transfer functions between the SSI frequencies);
  (b) started from the free field itself (GFAC = DFAC = 1, the manual's "same as in free-field") the
      near field is strain-compatible at once: converged at iteration 0 (D-NLS-02), with the same
      profile within the same tolerances.
* VP-N2  restart consistency (requirements 2.5, VP-22 applied to the iterations): (a) an iteration
  whose properties do not change (flat curves G/G_max = GFAC, D = DFAC beta; the free field is the
  low-strain soil of the near-field material, so GFAC G_layer = GFAC G_max) reproduces, through the
  "New Structure" restart (ANALYS <mode> 1, X_ff from COOXqqq), the FILE8 of the initiation run to
  1e-10; (b) repeating HOUSE with the same FILE74 strains gives the same properties (FILE78) and the
  restart reproduces the previous FILE8 and FILE74 to 1e-10.

Both problems write their decks with :mod:`sassi.verify.builders` and run the modules as RUNxxx does.
Disk use is kept small: one element per layer in plan and the factorisation records COOTKqqq (not
needed by a New Structure restart) are deleted after the initiation run.
"""
from __future__ import annotations

import math
import shutil
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .. import VPResult, problem
from ...core import nlsoil as NL
from ...core import shake as SH
from ...io import decks
from ...io.container import read_container
from ...io.thfile import write_history
from .. import builders as B

GRAV = 9.81
# ---- VP-N1 site: 10 m of sand (low-strain Vs 200 m/s, nu 1/3, rho 2 t/m3, 2 %) on rock
N1_H, N1_NS = 10.0, 10
N1_VS0, N1_RHO, N1_NU, N1_B0 = 200.0, 2.0, 1.0 / 3.0, 0.02
N1_ROCK = dict(vs=1000.0, vp=2000.0, rho=2.2, b=0.01)
N1_DT, N1_NFT = 0.01, 2048
N1_FCUT = 15.0                  # SSI cut-off = SOIL cut-off (EDUOPT,SOILCUTOFF), Hz
N1_STEP = 5                     # SSI frequencies every N1_STEP Fourier bins (0.244 Hz) up to the cut-off
N1_PGA = 0.25                   # rock outcrop peak acceleration, g
N1_ESF = 0.65                   # effective strain factor = SOIL <ratio>
N1_SOIL_ITER = 15               # SHAKE iterations of the reference (converged far below 0.1 %)
N1_TOL_G, N1_TOL_D = 0.05, 0.10
#: Seed & Idriss (1970) upper-range sand / Idriss (1990) damping (SHAKE91 Table B-1, sassi/data/dynp_library.pre)
SAND = dict(strain=[1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0, 10.0],
            g=[1.0, 1.0, 0.99, 0.96, 0.85, 0.64, 0.37, 0.18, 0.08, 0.05, 0.035],
            d=[0.24, 0.42, 0.80, 1.40, 2.80, 5.10, 9.80, 15.50, 21.0, 25.0, 28.0])


def synthetic_rock_motion(n: int = 1500, dt: float = N1_DT, pga: float = N1_PGA, seed: int = 20261002) -> np.ndarray:
    """Deterministic synthetic rock motion (g): Kanai-Tajimi filtered noise (f_g = 3 Hz, zeta_g = 0.6,
    4th-order high pass at 0.3 Hz, nothing above 20 Hz), trapezoidal-exponential envelope, scaled to ``pga``."""
    rng = np.random.default_rng(seed)
    t = np.arange(n) * dt
    nf = 4096
    W = np.fft.rfft(rng.standard_normal(n), nf)
    f = np.fft.rfftfreq(nf, dt)
    fg, zg = 3.0, 0.6
    kt = (fg ** 4 + (2 * zg * fg * f) ** 2) / ((fg ** 2 - f ** 2) ** 2 + (2 * zg * fg * f) ** 2 + 1e-12)
    hp = (f / 0.3) ** 4 / (1.0 + (f / 0.3) ** 4)
    a = np.fft.irfft(W * np.sqrt(kt * hp) * (f < 20.0), nf)[:n]
    env = np.where(t < 2.0, (t / 2.0) ** 2, np.where(t < 9.0, 1.0, np.exp(-0.5 * (t - 9.0))))
    a = a * env
    return a * (pga / np.max(np.abs(a)))


def soil_deck(model: str, layers: Sequence[Tuple[float, float, float, float, float, str]], rock: Dict[str, float],
              curves: Dict[str, Dict[str, List[float]]], thfile: str, nrval: int, dt: float, nft: int, ratio: float,
              iterations: int, cutoff: float) -> decks.Deck:
    """SOIL deck: rows (thick, vs, vp, rho, beta, label) top first, the rock half-space last (D-SOL-02), the
    motion an outcrop at the rock top (file with a dt header line), ACC001.TH (surface, within) saved."""
    d = decks.new("SOIL")
    d.params.update(dict(title=f"{model}: SOIL (SHAKE) reference", model=model, nrval=nrval, grav=GRAV, header=1,
                         outcrop=1, save=1, iter=iterations, ratio=ratio, gravmult=1.0, cof=0.0, soilcutoff=cutoff,
                         delt=dt, nft=nft, cl=len(layers) + 1, thfile=thfile, thtit="synthetic rock outcrop",
                         mult=1.0, max=0.0, indir=0, cmodform=0))
    for lab, cv in curves.items():
        for x, y in zip(cv["strain"], cv["g"]):
            d.table("dynp").append([lab, "G", x, y])
        for x, y in zip(cv["strain"], cv["d"]):
            d.table("dynp").append([lab, "D", x, y])
    for i, (h, vs, vp, rho, b, lab) in enumerate(layers, start=1):
        d.table("profile").append([i, h, rho * GRAV, vp, vs, b, b, lab])
    d.table("profile").append([len(layers) + 1, 0.0, rock["rho"] * GRAV, rock["vp"], rock["vs"], rock["b"], rock["b"], ""])
    d.table("sacc").append([1, 2, 0])
    return d


def near_field_column(site: B.Site, n_layers: int, mat_low: Tuple[float, float, float, float],
                      width: float = 1.0, per_layer: bool = False) -> Tuple[B.HouseBuilder, int, List[int]]:
    """FV excavation of the top ``n_layers`` TOPL layers, one element per layer in plan (``width`` x ``width``):
    excavated SOLIDs (ETYPE 2, the layer properties) and, at the same nodes, the near-field soil SOLIDs
    (structure, ETYPE 1) of the type-3 material ``mat_low`` = (vp, vs, rho, beta) -- the low-strain soil.
    ``per_layer``: one material (same properties) per layer, numbered from the top, so that the .pin can
    give every layer its own GFAC/DFAC line.  Every node is an interaction node (FV).  Returns (builder,
    near-field group, interaction nodes)."""
    hb = B.HouseBuilder(site, imp=0, incomp=1, title="near-field soil column in an FV excavation")
    z = site.elevations()
    h = 0.5 * width
    ids: Dict[Tuple[int, int, int], int] = {}
    for lev in range(n_layers, -1, -1):                       # bottom-up numbering (EDU-21)
        for j, y in enumerate((-h, h)):
            for i, x in enumerate((-h, h)):
                ids[(i, j, lev)] = hb.node(x, y, z[lev])
    gexc = hb.group(1, "excavated soil")
    gnl = hb.group(1, "near-field soil (nonlinear)")
    vp, vs, rho, beta = mat_low
    mid = hb.material(3, vp, vs, rho * site.gravity, beta, beta)
    mids = [mid] * n_layers
    if per_layer:                                             # distinct numbers (the builder merges equal rows)
        mids = [mid] + [len(hb.materials) + k for k in range(1, n_layers)]
        for m in mids[1:]:
            hb.materials.append([m, 3, float(vp), float(vs), float(rho * site.gravity), float(beta), float(beta)])
    for lev in range(n_layers):
        bot = [ids[(0, 0, lev + 1)], ids[(1, 0, lev + 1)], ids[(1, 1, lev + 1)], ids[(0, 1, lev + 1)]]
        top = [ids[(0, 0, lev)], ids[(1, 0, lev)], ids[(1, 1, lev)], ids[(0, 1, lev)]]
        hb.solid(bot + top, mat=lev + 1, etype=2, group=gexc)
        hb.solid(bot + top, mat=mids[lev], etype=1, group=gnl)
    inter = sorted(ids.values())
    hb.set_interaction(inter)
    return hb, gnl, inter


def stress_deck(model: str, thfile: str, dt: float, nft: int, interopt: int = 1) -> decks.Deck:
    d = decks.new("STRESS")
    d.params.update(dict(model=model, title=f"{model}: STRESS (non-linear soil strains)", iter=1, save=0,
                         thfile=thfile, mult=1.0, max=0.0, gravity=GRAV, delt=dt, nft=nft, df=1.0 / (dt * nft),
                         interopt=interopt))
    return d


def _drop_cootk(wd: Path) -> None:
    """The factorised systems COOTKqqq are not used by a New Structure restart (it re-assembles and
    re-factorises): remove them to keep the disk use of the problem small."""
    for p in wd.glob("COOTK*"):
        p.unlink()


def _iteration(wd: Path, model: str, fs: B.FrequencySet, mode: int) -> NL.Convergence:
    """HOUSE -> ANALYS (<mode>) -> STRESS; returns the convergence of this iteration (recorded)."""
    B.run("HOUSE", wd, model)
    B.run_analys(wd, model, fs, gravity=GRAV, mode=mode, save=1 if mode == 0 else 0)
    if mode == 0:
        _drop_cootk(wd)
    B.run("STRESS", wd, model)
    st = NL.status(wd)
    if st is None:
        raise B.ChainError("FILE78/FILE74 of the iteration missing")
    NL.record_convergence(wd, st)
    return st


# ======================================================================================
# VP-N1
# ======================================================================================
@problem("VP-N1", "Near-field soil iterations reproduce SOIL (SHAKE) for a uniform column without structure",
         tier="P1", modules=["SOIL", "SITE", "POINT", "HOUSE", "ANALYS", "STRESS"],
         source="requirements 4.4 item 7, 4.10 item 5, R1 section 6; D-NLS-02/03")
def vpn1(workdir: Path) -> VPResult:
    r = VPResult()
    wd = Path(workdir)
    m = "n1"
    vp0 = N1_VS0 * math.sqrt(2.0 * (1.0 - N1_NU) / (1.0 - 2.0 * N1_NU))
    acc = synthetic_rock_motion()
    write_history(wd / "rock.acc", acc, N1_DT)
    # ---- 1. SOIL (SHAKE): strain-compatible profile, FILE73 curves, surface motion ACC001.TH
    layers = [(N1_H / N1_NS, N1_VS0, vp0, N1_RHO, N1_B0, "Sand")] * N1_NS
    B.write_deck(wd, m, soil_deck(m, layers, N1_ROCK, {"Sand": SAND}, "rock.acc", len(acc), N1_DT, N1_NFT, N1_ESF,
                                  N1_SOIL_ITER, N1_FCUT))
    B.run("SOIL", wd, m)
    f88 = SH.read_file88(wd / "FILE88")
    gmax0 = N1_RHO * N1_VS0 ** 2
    ref_ratio = np.asarray(f88["G"]) / gmax0
    ref_beta = np.asarray(f88["beta_s"])
    # ---- 2. free field and excavated soil with the SOIL strain-compatible properties (SITEX,1 practice)
    site = B.layered_site([(f88["thick"][i], f88["Vs"][i], f88["Vp"][i], N1_RHO, f88["beta_s"][i], f88["beta_p"][i])
                           for i in range(N1_NS)], (N1_ROCK["vs"], N1_ROCK["vp"], N1_ROCK["rho"], N1_ROCK["b"]),
                          gravity=GRAV, nl=20)
    df = 1.0 / (N1_DT * N1_NFT)
    fn = sorted(set([1, 2, 3] + list(range(N1_STEP, int(N1_FCUT / df) + 1, N1_STEP))))
    fs = B.FrequencySet.fourier(N1_DT, N1_NFT, fn)
    hb, gnl, inter = near_field_column(site, N1_NS, (vp0, N1_VS0, N1_RHO, N1_B0), width=1.0, per_layer=True)
    B.run_soil(wd, m, site, fs, layer=N1_NS, rad=0.9)
    d = hb.deck(m)
    d["nlssi"] = 1
    B.write_deck(wd, m, d)
    B.write_deck(wd, m, stress_deck(m, "ACC001.TH", N1_DT, N1_NFT))

    def profile() -> Tuple[np.ndarray, np.ndarray]:
        rows = sorted(NL.update_rows(NL.read_nlfile(wd / "FILE74", "FILE74").rows, NL.load_curves(wd / "FILE73")),
                      key=lambda q: q.element)
        return np.array([q.ratio for q in rows]), np.array([q.beta for q in rows])

    # ---- 3a. near-field iterations from the low-strain properties: the .pin factors scale the free field
    #          (requirements 4.4 item 7), so layer i starts at GFAC_i G_layer_i = G_max, DFAC_i beta_i = beta_0
    gfac = [(N1_VS0 / float(f88["Vs"][i])) ** 2 for i in range(N1_NS)]
    dfac = [N1_B0 / float(f88["beta_s"][i]) for i in range(N1_NS)]
    lines = [NL.PinMaterial(gfac[i], dfac[i], 1) for i in range(N1_NS)]
    NL.write_pin(wd / f"{m}.pin", NL.PinData(N1_ESF, 1, [NL.PinGroup(gnl, N1_NS, N1_NS, 0, lines)]))
    hist: List[NL.Convergence] = [_iteration(wd, m, fs, mode=0)]
    start = NL.read_nlfile(wd / "FILE78", "FILE78")
    while not hist[-1].converged and len(hist) <= NL.MAX_ITERATIONS:
        hist.append(_iteration(wd, m, fs, mode=1))
    last = hist[-1]
    ratio, beta = profile()
    # ---- 3b. a new analysis from the free field itself (GFAC = DFAC = 1): strain-compatible at once
    NL.write_liq(wd / f"{m}.liq", 0)
    NL.write_pin(wd / f"{m}.pin", NL.PinData(N1_ESF, 1, [NL.PinGroup(gnl, N1_NS, N1_NS, 0,
                                                                      [NL.PinMaterial(1.0, 1.0, 1)] * N1_NS)]))
    ff = _iteration(wd, m, fs, mode=1)
    ff_start = NL.read_nlfile(wd / "FILE78", "FILE78")
    ratio_ff, beta_ff = profile()
    # ---- checks
    r.require("(a) iteration 0 starts from the low-strain properties: max |G_0/G_max - 1| = "
              f"{max(abs(q.ratio - 1.0) for q in start.rows):.2e}", all(abs(q.ratio - 1.0) < 1e-9 for q in start.rows))
    r.require(f"(a) converged (D-NLS-02) after {len(hist) - 1} restart iterations <= {NL.MAX_ITERATIONS}",
              last.converged)
    dG = [abs(c.dG_pct) for c in hist]
    r.require("(a) monotone convergence: max |dG/G| decreases at every iteration ("
              + ", ".join(f"{v:.3g}" for v in dG) + " %)", len(hist) >= 3 and all(b < a for a, b in zip(dG, dG[1:])))
    for i in range(N1_NS):
        r.check(f"(a) layer {i + 1}: G/Gmax near field vs SOIL", ratio[i], ref_ratio[i], rtol=N1_TOL_G)
    for i in range(N1_NS):
        r.check(f"(a) layer {i + 1}: damping near field vs SOIL", beta[i], ref_beta[i], rtol=N1_TOL_D)
    r.require("(b) iteration 0 = the free field: G_0/G_max = (Vs_layer/Vs_0)^2",
              all(abs(q.ratio - (float(f88["Vs"][k]) / N1_VS0) ** 2) < 1e-9
                  for k, q in enumerate(sorted(ff_start.rows, key=lambda q: q.element))))
    r.require(f"(b) started from the free field: converged at iteration 0 ({ff.text()})", ff.iteration == 0 and
              ff.converged)
    for i in range(N1_NS):
        r.check(f"(b) layer {i + 1}: G/Gmax near field vs SOIL", ratio_ff[i], ref_ratio[i], rtol=N1_TOL_G)
    for i in range(N1_NS):
        r.check(f"(b) layer {i + 1}: damping near field vs SOIL", beta_ff[i], ref_beta[i], rtol=N1_TOL_D)
    r.notes.append(f"SOIL (SHAKE, {N1_SOIL_ITER} iterations, ESF {N1_ESF}, cut-off {N1_FCUT:g} Hz): G/Gmax "
                   + " ".join(f"{v:.4f}" for v in ref_ratio) + "; damping " + " ".join(f"{v:.4f}" for v in ref_beta))
    r.notes.append("(a) near field from low strain: G/Gmax " + " ".join(f"{v:.4f}" for v in ratio) + "; damping "
                   + " ".join(f"{v:.4f}" for v in beta))
    r.notes.append("(a) iterations: " + "; ".join(c.text() for c in hist))
    r.notes.append("(b) near field from the free field: G/Gmax " + " ".join(f"{v:.4f}" for v in ratio_ff)
                   + "; damping " + " ".join(f"{v:.4f}" for v in beta_ff) + "; " + ff.text())
    r.notes.append(f"{len(fn)} SSI frequencies {fn[0] * df:.3g}-{fn[-1] * df:.3g} Hz, {len(inter)} interaction nodes, "
                   f"rock outcrop PGA {N1_PGA:g} g")
    return r


# ======================================================================================
# VP-N2
# ======================================================================================
N2_NS = 4
N2_STEP = 8


def _max_rel(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a), np.asarray(b)
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


@problem("VP-N2", "Restart consistency of the near-field soil iterations (New Structure restart)", tier="P1",
         modules=["HOUSE", "ANALYS", "STRESS"], source="requirements 2.5 (VP-22), 4.4 item 7; D-NLS-01")
def vpn2(workdir: Path) -> VPResult:
    r = VPResult()
    wd = Path(workdir)
    m = "n2"
    vp0 = N1_VS0 * math.sqrt(2.0 * (1.0 - N1_NU) / (1.0 - 2.0 * N1_NU))
    site = B.layered_site([(1.0, N1_VS0, vp0, N1_RHO, N1_B0)] * N2_NS,
                          (N1_ROCK["vs"], N1_ROCK["vp"], N1_ROCK["rho"], N1_ROCK["b"]), gravity=GRAV, nl=20)
    df = 1.0 / (N1_DT * N1_NFT)
    fn = sorted(set([1, 2] + list(range(N2_STEP, int(N1_FCUT / df) + 1, N2_STEP))))
    fs = B.FrequencySet.fourier(N1_DT, N1_NFT, fn)
    acc = synthetic_rock_motion(n=1200)
    write_history(wd / "ctrl.acc", acc, N1_DT)
    hb, gnl, _ = near_field_column(site, N2_NS, (vp0, N1_VS0, N1_RHO, N1_B0), width=1.0)
    B.run_soil(wd, m, site, fs, layer=N2_NS, rad=0.9)
    d = hb.deck(m)
    d["nlssi"] = 1
    B.write_deck(wd, m, d)
    B.write_deck(wd, m, stress_deck(m, "ctrl.acc", N1_DT, N1_NFT))
    # ---- (a) flat curves: G/Gmax = GFAC = 0.6 and D = DFAC beta = 4 % at every strain (the free-field layers
    #          are the near-field material here: iteration 0 = GFAC G_layer = 0.6 G_max, DFAC beta_layer = 4 %)
    gfac, dfac = 0.6, 2.0
    flat = {"Flat": dict(strain=[1e-4, 10.0], g=[gfac, gfac], d=[100 * dfac * N1_B0] * 2)}
    SH.write_file73(wd / "FILE73", {k: SH.DynamicProperty(k, np.array(v["strain"]), np.array(v["g"]),
                                                          np.array(v["strain"]), np.array(v["d"]))
                                    for k, v in flat.items()})
    NL.write_pin(wd / f"{m}.pin", NL.PinData(0.65, 1, [NL.PinGroup(gnl, 1, N2_NS, 1, [NL.PinMaterial(gfac, dfac, 1)])]))
    NL.write_liq(wd / f"{m}.liq", 0)
    B.run("HOUSE", wd, m)
    B.run_analys(wd, m, fs, gravity=GRAV, mode=0, save=1)
    _drop_cootk(wd)
    H0 = np.asarray(read_container(wd / "FILE8", "FILE8")["H"])
    B.run("STRESS", wd, m)
    p78_0 = NL.read_nlfile(wd / "FILE78", "FILE78")
    B.run("HOUSE", wd, m)                                   # iteration 1: properties from FILE74 + flat curves
    p78_1 = NL.read_nlfile(wd / "FILE78", "FILE78")
    B.run_analys(wd, m, fs, gravity=GRAV, mode=1)
    H1 = np.asarray(read_container(wd / "FILE8", "FILE8")["H"])
    r.check("(a) unchanged properties: max |G_1 - G_0| / G_0", _max_rel([q.G for q in p78_1.rows],
                                                                       [q.G for q in p78_0.rows]), 0.0, atol=1e-12)
    r.check("(a) unchanged properties: max |beta_1 - beta_0|", np.max(np.abs(np.array([q.beta for q in p78_1.rows])
                                                                            - np.array([q.beta for q in p78_0.rows]))),
            0.0, atol=1e-12)
    r.check("(a) FILE8 of the New Structure restart (iteration 1) vs the initiation FILE8", _max_rel(H1, H0), 0.0,
            atol=1e-10)
    r.require("(a) the iteration-1 run is a restart (FILE78 iteration 1, ANALYS <mode> 1)",
              p78_1.iteration == 1 and int(read_container(wd / "FILE8", "FILE8").meta.get("mode", -1)) == 1)
    # ---- (b) the Sand curves: iteration 1, then HOUSE again with the same strains (FILE74 of iteration 0)
    SH.write_file73(wd / "FILE73", {"Sand": SH.DynamicProperty("Sand", np.array(SAND["strain"]), np.array(SAND["g"]),
                                                               np.array(SAND["strain"]), np.array(SAND["d"]))})
    NL.write_pin(wd / f"{m}.pin", NL.PinData(0.65, 1, [NL.PinGroup(gnl, 1, N2_NS, 0, [NL.PinMaterial(1.0, 1.0, 1)])]))
    NL.write_liq(wd / f"{m}.liq", 0)
    B.run("HOUSE", wd, m)                                   # iteration 0 (low strain), FILE4 changes: COOX reused
    B.run_analys(wd, m, fs, gravity=GRAV, mode=1)
    B.run("STRESS", wd, m)
    shutil.copyfile(wd / "FILE74", wd / "FILE74_it0")
    B.run("HOUSE", wd, m)                                   # iteration 1 from FILE74 (iteration 0)
    B.run_analys(wd, m, fs, gravity=GRAV, mode=1)
    Ha = np.asarray(read_container(wd / "FILE8", "FILE8")["H"])
    f78a = NL.read_nlfile(wd / "FILE78", "FILE78")
    B.run("STRESS", wd, m)
    f74a = NL.read_nlfile(wd / "FILE74", "FILE74")
    shutil.copyfile(wd / "FILE74_it0", wd / "FILE74")      # the same strains again
    B.run("HOUSE", wd, m)
    f78b = NL.read_nlfile(wd / "FILE78", "FILE78")
    B.run_analys(wd, m, fs, gravity=GRAV, mode=1)
    Hb = np.asarray(read_container(wd / "FILE8", "FILE8")["H"])
    B.run("STRESS", wd, m)
    f74b = NL.read_nlfile(wd / "FILE74", "FILE74")
    r.require("(b) the properties changed in iteration 1 (the test is not trivial)",
              _max_rel([q.G for q in f78a.rows], [q.gmax for q in f78a.rows]) > 1e-3)
    r.check("(b) same strains -> same properties: max |G_b - G_a| / G_a", _max_rel([q.G for q in f78b.rows],
                                                                                  [q.G for q in f78a.rows]), 0.0,
            atol=1e-12)
    r.check("(b) restart reproduces the previous FILE8", _max_rel(Hb, Ha), 0.0, atol=1e-10)
    r.check("(b) and the effective strains of FILE74", _max_rel([q.gamma_eff for q in f74b.rows],
                                                                [q.gamma_eff for q in f74a.rows]), 0.0, atol=1e-10)
    r.notes.append(f"{len(fn)} SSI frequencies; iteration 1 G/Gmax "
                   + " ".join(f"{q.ratio:.4f}" for q in f78a.rows) + " (Sand curves, ESF 0.65)")
    return r
