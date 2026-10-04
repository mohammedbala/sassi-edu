"""Verification problems of SOIL-NON, the nonlinear time-domain option of the SOIL module (requirements 1.2
and 3.4.O, spec 05a section 6.6; numerical core :mod:`sassi.core.soilnon`, decisions SN-1 ... SN-12).

* VP-SN1  linear limit: with a very large reference strain (gamma_r = 1e6 %) the time-domain shear column
  (elastic half-space through the Joyner-Chen dashpot, outcrop input at bedrock, frequency-independent
  viscous damping) reproduces the frequency-domain linear SOIL solution (SHAKE recursion, 0 iterations)
  of the same profile and motion: surface PGA and 5 % spectral accelerations within 2 %.  The rigid
  base with a within input is checked on the PGA; Rayleigh and visco-elastic damping matched at the
  fundamental frequency are reported as informative comparisons.
* VP-SN2  single-element cyclic strain loop: the extended-Masing loop of the MKZ backbone has the secant
  modulus of the backbone at the strain amplitude and a loop-area damping equal to the analytical Masing
  damping (closed form for s = 1, quadrature of the backbone otherwise), 1e-3.
* VP-SN3  moderate shaking (SHAKE91 sample profile and motion, 0.1 g rock outcrop): surface PGA, spectra and
  strains of SOIL-NON (MKZ fitted to the same curves) and SOIL-EQL (8 iterations) -- informative, with the
  shear strain index of Kim et al. (2016) that predicts when the two methods differ.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy import integrate

from ...core import shake as SH
from ...core import soilnon as SN
from ...core import spectra as SP
from ...io import decks, textfiles, thfile
from ...modules.base import run_module
from ...modules.soil import write_deck
from .. import VPResult, problem
from .vp_soil import build_shake91_deck, parse_shake91_input, shake91_dir

PERIODS = (0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0)


def _write_nls(path: Path, options: SN.NLSoilOptions, layers: Sequence[SN.NLLayerData]) -> None:
    path.write_text(SN.nls_text(options, list(layers), path.stem), encoding="utf-8")


def _run(workdir: Path, model: str) -> Tuple[int, str]:
    rc = run_module("SOIL", model, workdir)
    return rc, (workdir / f"{model}_SOIL.out").read_text(encoding="utf-8")


# =======================================================================================
# VP-SN1  linear limit
# =======================================================================================
SN1_PROFILE = [(5.0, 18.0, 150.0), (5.0, 18.0, 150.0), (5.0, 19.0, 250.0), (5.0, 19.0, 250.0),
               (5.0, 20.0, 400.0), (5.0, 20.0, 400.0)]          # thickness m, unit weight kN/m3, Vs m/s
SN1_HALFSPACE = (22.0, 1200.0)
SN1_XI = 0.03
SN1_GRAV = 9.81


def _sn1_deck(workdir: Path, model: str, outcrop: int) -> None:
    nlay = len(SN1_PROFILE) + 1
    d = decks.new("SOIL")
    d.params.update(dict(title="VP-SN1 linear limit", model=model, nrval=1900, grav=SN1_GRAV, header=3,
                         outcrop=outcrop, save=1, iter=0, ratio=0.65, gravmult=1.0, cof=0.0, soilcutoff=0.0,
                         delt=0.02, nft=4096, cl=nlay, thfile=str(shake91_dir() / "DIAM.ACC"),
                         thtit="Loma Prieta 1989 Diamond Heights (SHAKE91 sample record)", mult=0.0, max=0.1,
                         indir=0, cmodform=0))
    for i, (h, w, vs) in enumerate(SN1_PROFILE, start=1):
        d.table("profile").append([i, h, w, 2.0 * vs, vs, SN1_XI, SN1_XI, ""])
    w, vs = SN1_HALFSPACE
    # the Joyner-Chen dashpot has no material damping (SN-7): the SHAKE half-space is undamped too
    d.table("profile").append([nlay, 0.0, w, 2.0 * vs, vs, 0.0, 0.0, ""])
    d.table("sacc").append([1, 2, 0])
    write_deck(workdir / f"{model}.soi", d)


def _surface(workdir: Path) -> Tuple[np.ndarray, float]:
    a, dt = textfiles.read_history(workdir / "ACC001.TH")
    return a, dt


def _spectrum(a: np.ndarray, dt: float) -> np.ndarray:
    return SP.response_spectrum(a, dt, 1.0 / np.asarray(PERIODS), [0.05])["SA"][0]


def _sn1_column() -> Tuple[SN.SoilNonColumn, float]:
    """The SOIL-NON column the module builds for the VP-SN1 profile (for the matched viscosity)."""
    thick = np.array([p[0] for p in SN1_PROFILE])
    rho = np.array([p[1] for p in SN1_PROFILE]) / SN1_GRAV
    vs = np.array([p[2] for p in SN1_PROFILE])
    models = [SN.SublayerModel(1.0, 1.0, 1e6, SN1_XI) for _ in SN1_PROFILE]
    nel = SN.elements_per_sublayer(thick, vs, min(SN.F_MAX_DISCRETISATION, 0.5 / 0.02))
    col = SN.SoilNonColumn(thick, rho, vs, models, SN1_HALFSPACE[0] / SN1_GRAV, SN1_HALFSPACE[1], nel)
    w, _ = SN.fixed_base_modes(col)
    return col, float(w[0])


@problem("VP-SN1", "SOIL-NON linear limit equals the linear SOIL (SHAKE) solution", tier="P2", modules=["SOIL"],
         source="package soil_non; R1 section 6 (SHAKE recursion); Joyner & Chen (1975); Phillips & Hashash (2009)",
         slow=True)
def vpsn1(workdir: Path) -> VPResult:
    """Six 5 m sublayers (Vs 150/250/400 m/s, 3 % damping) on an elastic half-space (1200 m/s, undamped),
    SHAKE91 sample record scaled to 0.1 g as rock outcrop motion; SOIL-NON with gamma_r = 1e6 % (NLSLAYER
    curvefit 0, B = S = 1), damping type 1, default sub-increments and tolerances; SOIL-EQL with 0
    iterations.  Criteria (package spec): surface acceleration within 2 % (PGA and 5 % SA at 0.05-2 s)."""
    r = VPResult()
    workdir = Path(workdir)
    nlay = len(SN1_PROFILE) + 1
    lin_layers = [SN.NLLayerData(k, 0, 1.0, 1.0, 1e6, 0.0) for k in range(1, nlay)]
    cases: Dict[str, Tuple[np.ndarray, float]] = {}

    def run_case(tag: str, outcrop: int, opts=None, layers=None) -> bool:
        wd = workdir / tag
        wd.mkdir(parents=True, exist_ok=True)
        _sn1_deck(wd, "sn1", outcrop)
        if opts is not None:
            _write_nls(wd / "sn1.nls", opts, layers)
        rc, out = _run(wd, "sn1")
        if not r.require(f"SOIL run '{tag}' status 0", rc == 0):
            r.notes.append(f"{tag}: " + out[-2500:])
            return False
        cases[tag] = _surface(wd)
        return True

    # the discrete column of SOIL-NON vs the continuum: small-strain amplification surface / outcrop
    col, _ = _sn1_column()
    f = np.arange(0.2, 8.0, 0.005)
    Hd = np.abs(SN.linear_transfer(col, SN.build_damping(col, 1, 1), 1, f, nodes=[0])[:, 0])
    w = np.array([p[1] for p in SN1_PROFILE] + [SN1_HALFSPACE[0]])
    v = np.array([p[2] for p in SN1_PROFILE] + [SN1_HALFSPACE[1]])
    csh = SH.SoilColumn(np.array([p[0] for p in SN1_PROFILE] + [0.0]), w, v,
                        np.array([SN1_XI] * len(SN1_PROFILE) + [0.0]), SN1_GRAV)
    Hc = np.abs(SH.column_transfer(f, csh, csh.gmax, csh.beta0, 1, False, nlay, True, 0))
    kd, kc = int(np.argmax(Hd)), int(np.argmax(Hc))
    r.check("small-strain amplification peak |surface / outcrop|, discrete SOIL-NON model vs continuum", Hd[kd],
            Hc[kc], rtol=0.01)
    r.check("frequency of the small-strain amplification peak (Hz), discrete vs continuum", f[kd], f[kc], rtol=0.01)
    ok = run_case("eql_elastic", 1)
    ok &= run_case("non_elastic", 1, SN.NLSoilOptions(opt=1, bedint=1, damptype=1), lin_layers)
    if ok:
        ae, dt = cases["eql_elastic"]
        an, _ = cases["non_elastic"]
        r.check("elastic base, outcrop input: surface PGA (g), SOIL-NON vs SOIL-EQL", float(np.max(np.abs(an))),
                float(np.max(np.abs(ae))), rtol=0.02)
        se, sn = _spectrum(ae, dt), _spectrum(an, dt)
        for T, x, y in zip(PERIODS, sn, se):
            r.check(f"elastic base: 5 % SA at T = {T:g} s (g), SOIL-NON vs SOIL-EQL", float(x), float(y), rtol=0.02)
        rms = float(np.sqrt(np.mean((an - ae) ** 2) / np.mean(ae ** 2)))
        r.inform("elastic base: rms difference of the surface histories / rms (informative)", rms, 0.0)
    ok = run_case("eql_rigid", 0)
    ok &= run_case("non_rigid", 0, SN.NLSoilOptions(opt=1, bedint=0, damptype=1), lin_layers)
    if ok:
        ae, dt = cases["eql_rigid"]
        an, _ = cases["non_rigid"]
        r.check("rigid base, within input at the base: surface PGA (g), SOIL-NON vs SOIL-EQL",
                float(np.max(np.abs(an))), float(np.max(np.abs(ae))), rtol=0.02)
        se, sn = _spectrum(ae, dt), _spectrum(an, dt)
        k = int(np.argmax(np.abs(sn / se - 1.0)))
        r.inform(f"rigid base: largest 5 % SA deviation (at T = {PERIODS[k]:g} s), SOIL-NON vs SOIL-EQL",
                 float(sn[k]), float(se[k]),
                 note="no radiation damping at a rigid base: sharp resonances; converges with more "
                      "sub-increments and elements")
    # informative: frequency-dependent viscous damping matched at the fundamental frequency
    if "eql_elastic" in cases:
        ae, dt = cases["eql_elastic"]
        se = _spectrum(ae, dt)
        col, w1 = _sn1_column()
        visc = [SN.NLLayerData(k, 0, 1.0, 1.0, 1e6, 2.0 * SN1_XI * col.rho[k - 1] * col.vs[k - 1] ** 2 / w1)
                for k in range(1, nlay)]
        for tag, opts, layers, label in (
                ("non_rayleigh", SN.NLSoilOptions(opt=1, bedint=1, damptype=3), lin_layers,
                 "Rayleigh damping matched at f1 and 5 f1"),
                ("non_visco", SN.NLSoilOptions(opt=1, bedint=1, damptype=2), visc,
                 "visco-elastic damping matched at f1")):
            if run_case(tag, 1, opts, layers):
                an, _ = cases[tag]
                r.inform(f"{label}: surface PGA (g) vs SOIL-EQL", float(np.max(np.abs(an))), float(np.max(np.abs(ae))))
                sn = _spectrum(an, dt)
                k = int(np.argmax(np.abs(sn / se - 1.0)))
                r.inform(f"{label}: largest 5 % SA deviation (T = {PERIODS[k]:g} s) vs SOIL-EQL", float(sn[k]),
                         float(se[k]))
    r.notes.append("SOIL-EQL is the SHAKE recursion with the SASSI complex modulus (cmodform 0) and the "
                   "frequency-independent damping 3 %; SOIL-NON integrates the lumped-mass column (odd number of "
                   "elements per sublayer, h <= Vs / 250 Hz) with Newmark average acceleration, band-limited "
                   "input, the modal frequency-independent damping (every fixed-base mode 3 %) and, for the "
                   "elastic base, the Joyner-Chen dashpot rho_r V_r.  Rayleigh (f1, 5 f1) and Kelvin-Voigt (f1) "
                   "damping are frequency dependent and differ by more than the criterion: informative only.")
    return r


# =======================================================================================
# VP-SN2  single-element Masing loop
# =======================================================================================
def _masing_reference(ga: float, gr: float, beta: float, s: float) -> float:
    """Analytical Masing damping of the MKZ backbone: the closed form of the hyperbolic model for s = 1
    (Darendeli 2001, R2 H.2, gamma_r -> gamma_r / beta) and quadrature of the backbone otherwise."""
    if s == 1.0:
        g, grr = ga, gr / beta
        return (1.0 / math.pi) * (4.0 * (g - grr * math.log((g + grr) / grr)) / (g * g / (g + grr)) - 2.0)
    F = lambda x: x / (1.0 + beta * (x / gr) ** s)           # noqa: E731  (G0 = 1)
    integral, _ = integrate.quad(F, 0.0, ga, epsabs=0.0, epsrel=1e-12, limit=200)
    return (2.0 / math.pi) * (2.0 * integral / (ga * F(ga)) - 1.0)


@problem("VP-SN2", "SOIL-NON single-element Masing loop: secant modulus and loop-area damping", tier="P2",
         modules=["SOIL"], source="package soil_non; Masing (1926); Kramer (1996) 6.4.3; R2 H.2 (D_Masing)")
def vpsn2(workdir: Path) -> VPResult:
    """An MKZ element (G0 = 1e5, gamma_r = 0.05 %) is loaded on the backbone to the amplitude gamma_a and
    cycled twice between -gamma_a and +gamma_a (2000 strain increments per half cycle), for nine parameter
    sets (beta, s) x gamma_a / gamma_r (one material point each, driven together).  The secant modulus of
    the loop equals the backbone secant at gamma_a and the loop area (trapezoidal rule over the second
    cycle) gives the Masing damping A / (2 pi gamma_a tau_a) (1e-3)."""
    r = VPResult()
    G0, gr = 1.0e5, 0.05e-2
    n = 2000
    sets = [(beta, s, ratio) for beta, s in ((1.0, 1.0), (1.3, 0.85), (1.0, 0.7)) for ratio in (0.1, 1.0, 10.0)]
    beta = np.array([p[0] for p in sets])
    s = np.array([p[1] for p in sets])
    ga = np.array([p[2] for p in sets]) * gr
    m = SN.MasingMKZ(np.full(len(sets), G0), np.full(len(sets), gr), beta, s)
    unit_up = np.linspace(0.0, 1.0, n + 1)[1:]
    unit_cycle = np.concatenate([np.linspace(1.0, -1.0, 2 * n + 1)[1:], np.linspace(-1.0, 1.0, 2 * n + 1)[1:]])
    m.drive(unit_up[:, None] * ga[None, :])                                # virgin loading
    m.drive(unit_cycle[:, None] * ga[None, :])                             # first cycle
    tau = m.drive(unit_cycle[:, None] * ga[None, :])                       # second cycle
    trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    for j, (bj, sj, ratio) in enumerate(sets):
        tau_a = float(SN.mkz_stress(ga[j], G0, gr, bj, sj))
        tag = f"beta {bj:g}, s {sj:g}, gamma_a = {ratio:g} gamma_r"
        secant = (tau[-1, j] - tau[2 * n - 1, j]) / (2.0 * ga[j])
        r.check(f"{tag}: loop secant modulus / backbone secant", secant / (tau_a / ga[j]), 1.0, rtol=1e-3)
        g = np.concatenate([[ga[j]], unit_cycle * ga[j]])
        t = np.concatenate([[tau_a], tau[:, j]])
        d_loop = abs(float(trap(t, g))) / (2.0 * math.pi * ga[j] * tau_a)
        d_ref = _masing_reference(float(ga[j]), gr, bj, sj)
        r.check(f"{tag}: loop-area damping vs analytical Masing damping", d_loop, d_ref, rtol=1e-3)
        r.check(f"{tag}: SOIL-NON masing_damping (FILE88) vs analytical", float(SN.masing_damping(ga[j], gr, bj, sj)),
                d_ref, rtol=1e-6)
    r.notes.append("Masing loop of amplitude ga: branches tau = +-tau_a + 2 F((g -+ ga)/2); damping "
                   "D = (2/pi)(2 int_0^ga F dg / (ga F(ga)) - 1).  The s = 1 reference is the hyperbolic closed form "
                   "printed in R2 H.2 (Darendeli 2001) with gamma_r / beta; other s by adaptive quadrature.")
    return r


# =======================================================================================
# VP-SN3  moderate shaking: SOIL-NON vs SOIL-EQL (informative)
# =======================================================================================
def _strain_index(acc_g: np.ndarray, dt: float, grav: float, thick: Sequence[float], vs: Sequence[float]) -> float:
    """Shear strain index I_gamma = PGV of the input / V_S30 (Kim et al. 2016), in percent."""
    v = np.cumsum(acc_g * grav) * dt
    v = v - np.mean(v)
    z30 = 30.0 / 0.3048 if grav > 20.0 else 30.0            # 30 m in the model units
    tt, z = 0.0, 0.0
    for h, c in zip(thick, vs):
        hh = min(h, z30 - z)
        if hh <= 0:
            break
        tt += hh / c
        z += hh
    return 100.0 * float(np.max(np.abs(v))) / (z / tt)


@problem("VP-SN3", "SOIL-NON vs SOIL-EQL at moderate shaking (SHAKE91 sample, informative)", tier="P2",
         modules=["SOIL"], source="package soil_non; R2 H.1 (SHAKE91 sample); Kim et al. (2016) Earthquake Spectra 32")
def vpsn3(workdir: Path) -> VPResult:
    """SHAKE91 manual example (16 sublayers, 0.1 g rock outcrop motion at the base): SOIL-EQL with 8
    iterations (as VP-04) and SOIL-NON with the MKZ backbone fitted to the same G/Gmax curves (NLSLAYER
    curvefit 1), elastic base, frequency-independent damping D_min of the curves.  The comparison is
    informative: equivalent-linear and nonlinear analyses are expected to agree within an engineering band
    of about +-25 % in PGA when the shear strain index is small (Kim et al. 2016: practically identical for
    I_gamma < 0.03 %, increasingly different above)."""
    r = VPResult()
    workdir = Path(workdir)
    src = shake91_dir()
    inp = parse_shake91_input(src / "INP.DAT")
    nl = len(inp["layers"])
    res: Dict[str, Dict] = {}
    for tag in ("eql", "non"):
        wd = workdir / tag
        wd.mkdir(parents=True, exist_ok=True)
        build_shake91_deck(inp, wd, "s91", str(src / inp["motion"]["file"]))
        if tag == "non":
            _write_nls(wd / "s91.nls", SN.NLSoilOptions(opt=1, bedint=1, damptype=1),
                       [SN.NLLayerData(k, 1) for k in range(1, nl)])
        rc, out = _run(wd, "s91")
        if not r.require(f"SOIL-{tag.upper()} run status 0", rc == 0):
            r.notes.append(out[-2500:])
            return r
        a1, dt = textfiles.read_history(wd / "ACC001.TH")
        prof = {}
        for lay in range(1, nl + 1):
            a, _ = textfiles.read_history(wd / f"ACC{lay:03d}.TH")
            prof[lay] = float(np.max(np.abs(a)))
        strains = {}
        for lay in range(1, nl):
            g, _ = textfiles.read_history(wd / f"SN{lay:03d}.TH")
            strains[lay] = float(np.max(np.abs(g)))
        res[tag] = dict(a1=a1, dt=dt, pga=prof, strain=strains, out=out)
    m = inp["motion"]
    acc = thfile.read_soil_history(src / m["file"], m["nv"], m["nhead"])
    acc = thfile.scale_history(acc, m["mult"], m["xmax"])
    thick = [ly["thick"] for ly in inp["layers"][:-1]]
    vs = [ly["vs"] for ly in inp["layers"][:-1]]
    ig = _strain_index(acc, m["dt"], 32.2, thick, vs)
    r.inform("shear strain index I_gamma = PGV_input / V_S30 (%)", ig, 0.03,
             note="Kim et al. (2016): EL and NL practically identical below about 0.03 %")
    e, n = res["eql"], res["non"]
    r.inform("surface PGA (g): SOIL-NON vs SOIL-EQL", n["pga"][1], e["pga"][1],
             note="documented engineering band +-25 % (informative)")
    depth = np.concatenate([[0.0], np.cumsum(thick)])
    for lay in (5, 9, 13):
        r.inform(f"PGA at {depth[lay - 1]:.0f} ft (g): SOIL-NON vs SOIL-EQL", n["pga"][lay], e["pga"][lay])
    for lay in (4, 8, 12):
        r.inform(f"max shear strain sublayer {lay} (%): SOIL-NON vs SOIL-EQL", n["strain"][lay], e["strain"][lay])
    T = (0.2, 0.5, 1.0)
    sa_e = SP.response_spectrum(e["a1"], e["dt"], 1.0 / np.asarray(T), [0.05])["SA"][0]
    sa_n = SP.response_spectrum(n["a1"], n["dt"], 1.0 / np.asarray(T), [0.05])["SA"][0]
    for Tk, x, y in zip(T, sa_n, sa_e):
        r.inform(f"surface 5 % SA at T = {Tk:g} s (g): SOIL-NON vs SOIL-EQL", float(x), float(y))
    ratio = n["pga"][1] / e["pga"][1]
    r.require("SOIL-NON listing reports the equivalent-linear properties of the nonlinear run",
              "Maximum strains and equivalent-linear properties" in n["out"])
    r.notes.append(f"Surface PGA SOIL-NON / SOIL-EQL = {ratio:.3f} "
                   f"({'inside' if abs(ratio - 1.0) <= 0.25 else 'outside'} the +-25 % band); shear strain index "
                   f"{ig:.3f} %.  Differences come from the methods: SOIL-EQL uses one strain-compatible modulus and "
                   "damping per sublayer for the whole record, SOIL-NON follows the MKZ backbone (fitted to the "
                   "same G/Gmax curves) cycle by cycle with Masing damping and the small-strain damping D_min.")
    return r
