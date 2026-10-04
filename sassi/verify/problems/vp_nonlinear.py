"""Verification problems of Option NON (NONLINEAR module, hysteresis models, SHEAR / BBCGEN; requirements 4.15,
6.3; spec 05d section 8 items 7-8; spec 10 section 3.11; spec 11 section 7 item 3).

* VP-45    Option NON element models: an elastic-perfectly-plastic BBC cycled at x > x_y with the General
           Masing Rule gives a closed loop, ``xi_h = 2 (x - x_y)/(pi x)`` and ``K_sec = F_y/x``; the NONLINEAR
           module reproduces them from a relative-displacement history, and its state machine works
           (no .NON -> elastic run + SPRING.NON = 1; .NON = 1 -> iteration; SPRING_EQL_Matl_Prop.txt used)
* VP-46    SHEAR / BBCGEN numbers of the spec 10 section 3.11 worked example (British and SI) and the
           22-point BBCGEN curve (yield index 21, point 22 = (0.02, 1.02 V_u))
* VP-NON1  SDOF with a nonlinear spring on a rigid base under harmonic input, iterated through the whole
           chain HOUSE -> ANALYS (New Structure restart) -> MOTION -> RELDISP -> NONLINEAR: the converged
           equivalent-linear stiffness and damping equal the secant stiffness and the loop-area damping of
           the hysteresis model at the converged amplitude, every SSI response equals the closed-form
           steady state of the equivalent-linear SDOF, the iteration follows an independent closed-form
           fixed-point iteration, and it converges within the limits of D-NON-06

The SSI models are written as module decks with :mod:`sassi.verify.builders` and run exactly as the RUN
commands do; NONLINEAR is run with :func:`sassi.modules.nonlinear.run_nonlinear` (RUNNONLINEAR).
"""
from __future__ import annotations

import math
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from ...core import hysteresis as HY
from ...core import panels as PN
from ...io import decks
from ...modules import nonlinear as NLM
from .. import VPResult, problem, worse
from .. import builders as B


# ======================================================================================
# VP-45
# ======================================================================================
VP45_XY, VP45_FY = 0.01, 100.0
VP45_AMPLITUDES = (0.0125, 0.015, 0.02, 0.05, 0.1, 0.5)


def _epp_closed_form(x: float, xy: float = VP45_XY) -> Tuple[float, float]:
    """(xi_h, K_sec / K_el) of an elastic-perfectly-plastic Masing loop, written out independently:
    loop area 4 F_y (x - x_y), strain energy F_y x / 2."""
    return 2.0 * (x - xy) / (math.pi * x), xy / x


@problem("VP-45", "Option NON element models: EPP backbone under the General Masing Rule, NONLINEAR state machine",
         tier="P2", modules=["NONLINEAR"], source="requirements 6.3 VP-45; spec 05d section 8 items 7-8")
def vp45(workdir: Path) -> VPResult:
    r = VPResult()
    bb = HY.Backbone([VP45_XY, 100.0], [VP45_FY, VP45_FY], 1)
    for x in VP45_AMPLITUDES:
        loop = HY.cyclic_loop(HY.MasingGMR(bb), x)
        xi, ks = _epp_closed_form(x)
        r.check(f"x = {x / VP45_XY:g} x_y: loop closed |F(end) - F(start)|/F_y", loop.closed / VP45_FY, 0.0, atol=1e-12)
        r.check(f"x = {x / VP45_XY:g} x_y: xi_h vs 2(x - x_y)/(pi x)", loop.damping, xi, rtol=1e-6)
        r.check(f"x = {x / VP45_XY:g} x_y: K_sec vs F_y/x", loop.secant, VP45_FY / x, rtol=1e-6)
        r.check(f"x = {x / VP45_XY:g} x_y: peak forces +-F_y", 0.5 * (loop.f_pos - loop.f_neg), VP45_FY, rtol=1e-12)
    # ---- the module: a spring whose deformation is a harmonic history of amplitude 3 x_y
    wd = Path(workdir) / "module"
    shutil.rmtree(wd, ignore_errors=True)
    wd.mkdir(parents=True)
    amp = 3.0 * VP45_XY
    nt, dt = 400, 0.01
    t = np.arange(nt) * dt
    x = amp * np.sin(2.0 * math.pi * 2.5 * t)          # 2.5 Hz: the samples hit the peaks exactly
    _vp45_model(wd, bb, x, dt)
    rc = NLM.run_nonlinear("m", wd)
    r.require("NONLINEAR run 1 (no SPRING.NON) finishes with status OK", rc == 0)
    lst = (wd / "m_NONLINEAR.out").read_text(encoding="utf-8")
    r.require("run 1 is an ELASTIC run", "ELASTIC run" in lst)
    r.require("run 1 creates SPRING.NON containing 1", (wd / "SPRING.NON").exists()
              and (wd / "SPRING.NON").read_text().strip() == "1")
    it, rows = NLM.read_props(wd / "SPRING_EQL_Matl_Prop.txt")
    xi, ks = _epp_closed_form(amp)
    r.check("run 1: properties written for SSI iteration 1", it, 1, atol=0)
    r.check("module: K_sec/K_el at x = 3 x_y vs x_y/x", rows[1].ratio, ks, rtol=1e-6)
    r.check("module: xi_h at x = 3 x_y vs 2(x - x_y)/(pi x)", rows[1].xi_h, xi, rtol=1e-6)
    new = decks.read(wd / "m_new.hou", "HOUSE")
    sc = new.rows("springprops")[0]
    r.check("m_new.hou: spring constant k_el K_sec/K_el", sc["scx"], 1.0e4 * ks, rtol=1e-6)
    r.check("m_new.hou: spring damping = xi_h (hysteretic only, EQL <ElasicD> 0)", sc["damp"], xi, rtol=1e-6)
    f, _ = NLM.read_history(wd / "SPRING0001.ths")
    xs, _ = NLM.read_history(wd / "SPRING0001.thd")
    k = int(np.argmax(xs))
    r.check("time domain: F at max|x| = F_y (plastic plateau)", f[k], VP45_FY, rtol=1e-12)
    r.check("time domain: F_mu = K_el max|x| / F = 3", (wd / "Spring.fmu").exists() and float(
        np.loadtxt(wd / "Spring.fmu")[-1]), 3.0, rtol=1e-7)
    # ---- run 2: SPRING.NON = 1 -> iteration; the properties of SPRING_EQL_Matl_Prop.txt are the 'used' ones
    marker = 0.111
    p = wd / "SPRING_EQL_Matl_Prop.txt"
    p.write_text(p.read_text(encoding="utf-8").replace(f"{rows[1].xi_new:>23.16e}", f"{marker:>23.16e}"),
                 encoding="utf-8")
    shutil.copyfile(wd / "m_new.hou", wd / "m.hou")
    rc = NLM.run_nonlinear("m", wd)
    lst = (wd / "m_NONLINEAR.out").read_text(encoding="utf-8")
    r.require("NONLINEAR run 2 (SPRING.NON = 1) finishes with status OK", rc == 0)
    r.require("run 2 is a nonlinear iteration processing SSI analysis 1",
              "nonlinear iteration" in lst and "SSI analysis processed: 1" in lst)
    conv = NLM.read_convergence(wd)[-1]
    r.check("run 2: the damping 'used' is the one of SPRING_EQL_Matl_Prop.txt (|d xi| vs |marker - xi_h|)",
            conv["max_dxi_pct"], 100.0 * abs(marker - xi), rtol=1e-4)
    it2, _ = NLM.read_props(p)
    r.check("run 2: properties written for SSI iteration 2", it2, 2, atol=0)
    # ---- deleting both state files restarts from the elastic analysis
    (wd / "SPRING.NON").unlink()
    p.unlink()
    rc = NLM.run_nonlinear("m", wd)
    r.require("after deleting SPRING.NON and SPRING_EQL_Matl_Prop.txt: an ELASTIC run again",
              rc == 0 and "ELASTIC run" in (wd / "m_NONLINEAR.out").read_text(encoding="utf-8"))
    r.notes.append("GMR = Masing rule with factor 2 and the extended memory rules (D-NON-03); the closed form is "
                   "written out independently of sassi.core.hysteresis: loop area 4 F_y (x - x_y).")
    return r


def _vp45_model(wd: Path, bb: HY.Backbone, x: np.ndarray, dt: float) -> None:
    site = B.layered_site([(10.0, 300.0, 600.0, 2.0, 0.05)], (800.0, 1600.0, 2.2, 0.02))
    hb = B.HouseBuilder(site)
    a = hb.add_node(0.0, 0.0, 0.0, fix=(0, 0, 0, 1, 1, 1))
    b = hb.add_node(0.0, 0.0, 1.0, fix=(0, 1, 1, 1, 1, 1))
    hb.set_interaction([a])
    hb.spring(a, b, (1.0e4, 0, 0, 0, 0, 0), 0.02)
    B.write_deck(wd, "m", hb.deck("m"))
    e = NLM.EqlDeck()
    e.params.update(model="m", edf=1.0, nonlinopts=NLM.OPT_SPRINGS, elasticd=0, gravity=site.gravity)
    for p, (xx, yy) in enumerate(zip(bb.x, bb.y), start=1):
        e.add("bbc", **{"bbc": 1, "type": 4, "yield": 1, "point": p, "x": float(xx), "y": float(yy)})
    e.add("springs", num=1, group=1, elem=1, bbc=1, disp=1, force=4, prop=1, node_i=a, node_j=b, scx=1.0e4, scy=0.0,
          scz=0.0, scxx=0.0, scyy=0.0, sczz=0.0, damp=0.02, k_el=1.0e4)
    NLM.write_eql(wd / "m.eql", e)
    from ...conventions import nodal_result_name
    NLM.write_history(wd / nodal_result_name(a, 1, "THD"), np.zeros_like(x), dt)
    NLM.write_history(wd / nodal_result_name(b, 1, "THD"), x, dt)


# ======================================================================================
# VP-46
# ======================================================================================
VP46_BRITISH = dict(aci=4424.82, upper=6109.40, wood_raw=648.00, wood_lower=3665.64, barda=4026.77, gw=2405.90)
VP46_GW_ABE = 2617.54
VP46_SI = dict(aci=19682.6, barda=17912.0, gw=10702.0)
KIPS_TO_KN = 4.44822162           # the conversion of spec 10 section 3.11 (lb -> kN = 0.00444822162)


@problem("VP-46", "SHEAR / BBCGEN: shear capacities of the spec 10 worked example and the 22-point backbone",
         tier="P2", modules=["SHEAR", "BBCGEN"], source="spec 10 section 3.11; spec 11 section 7 item 3")
def vp46(workdir: Path) -> VPResult:
    """h_W = 20 ft, l_W = 30 ft, t_W = 2 ft, f'c = 5 ksi, f_y = 60 ksi, rho = 0.005, N_U = 1000 kips, A_BE = 0;
    the SI version 6.096 m x 9.144 m x 0.6096 m, f'c = 34473.8 kN/m2, f_y = 413685 kN/m2, N_U = 4448.22 kN."""
    r = VPResult()
    c = PN.shear_capacities(20.0, 30.0, 2.0, 5.0, 60.0, 0.005, 1000.0, gravity=32.2)
    names = dict(aci="V_ACI 318-08 (raw = capped)", upper="upper bound 10 sqrt(f'c) A_W", wood_raw="V_Wood raw",
                 wood_lower="Wood lower bound 6 sqrt(f'c) A_W", barda="V_Barda 1977", gw="V_Gulec-Whittaker 2009")
    for key, ref in VP46_BRITISH.items():
        r.check(f"British: {names[key]} (kips)", getattr(c, key), ref, rtol=1e-6)
    c2 = PN.shear_capacities(20.0, 30.0, 2.0, 5.0, 60.0, 0.005, 1000.0, a_be=0.1, fy_be=60.0, gravity=32.2)
    r.check("British: V_GW with A_BE = 0.1 ft2, f_y,BE = 60 ksi (kips)", c2.gw, VP46_GW_ABE, rtol=1e-6)
    s = PN.shear_capacities(6.096, 9.144, 0.6096, 34473.8, 413685.0, 0.005, 4448.22, gravity=9.81)
    # D-W3-08: the SI references are printed to 0.1 kN, so they are compared within half their last digit
    # (0.05 kN); the 1e-6 relative check of the SI path is the British-conversion check below.
    for key, ref in VP46_SI.items():
        r.check(f"SI: {names[key]} (kN) vs the printed reference (printed to 0.1 kN)", getattr(s, key), ref,
                atol=0.05)
    # independent consistency of the SI path: the British results converted to kN
    for key in VP46_SI:
        r.check(f"SI: {names[key]} vs the British result x {KIPS_TO_KN} kN/kip (input rounding of the SI data)",
                getattr(s, key), getattr(c, key) * KIPS_TO_KN, rtol=1e-6)
    # ---- BBCGEN through the commands on a 30 ft x 20 ft x 2 ft panel model (British units)
    from ...prep import Interpreter
    from ...prep.commands import nonlinear_cmds as NC
    ui = Interpreter(cwd=workdir)
    ui.run_text("""GRAVITY,32.2
N,1,0,0,0
N,4,30,0,0
FILL,1,4
NGEN,2,4,1,4,1,0,0,10
M,1,580000,0.2,0.15,0.04,0.04,1
GROUP,1,SHELL
MACT,1
E,1,1,2,6,5
EGEN,2,1,1
EGEN,1,4,1,3
THICK,1,6,1,2
PNLGEN
P,2,1,2,1,1
BBCGEN,1,4,5,60,0.005,1000,0,0,0
BBCGEN,2,4,5,60,0.005,1000,0,0,0.3
""")
    cv = NC.curve(ui.model, 1)
    gw = c.gw
    r.check("BBCGEN: number of points", len(cv["x"]), 22, atol=0)
    r.check("BBCGEN: yield point index", cv["yield_"], 21, atol=0)
    r.check("BBCGEN: point 22 strain", cv["x"][21], 0.02, rtol=1e-12)
    r.check("BBCGEN: point 22 force = 1.02 V_u (Gulec-Whittaker)", cv["y"][21], 1.02 * gw, rtol=1e-9)
    r.check("BBCGEN: yield point (0.004, V_u): force", cv["y"][20], gw, rtol=1e-9)
    r.check("BBCGEN: yield point strain", cv["x"][20], 0.004, rtol=1e-12)
    G = 580000.0 / 2.4
    ga = G * 2.0 * 30.0
    # 3 sqrt(f'c [psi]) A_W [in2] in lb, written out (final audit: the reference was the program's own
    # PN.shear_capacities value, the one BBCGEN stores): f'c = 5000 psi, A_W = 360 in x 24 in
    cracking = 3.0 * math.sqrt(5000.0) * (30.0 * 12.0) * (2.0 * 12.0) / 1000.0
    r.check("BBCGEN: cracking force 3 sqrt(f'c) A_W (kips)", cv["y"][0], cracking, rtol=1e-9)
    r.check("BBCGEN: initial slope Y1/X1 = G A_W (EDU-09 consistent)", cv["y"][0] / cv["x"][0], ga, rtol=1e-9)
    steps = np.diff(cv["x"][:21])
    r.check("BBCGEN: points 2-21 equally spaced in strain (max deviation / step)",
            float(np.max(np.abs(steps - steps.mean())) / steps.mean()), 0.0, atol=1e-9)
    cv2 = NC.curve(ui.model, 2)
    r.check("BBCGEN CrackingForceLevel 0.3: Y1 = 0.3 V_u", cv2["y"][0], 0.3 * gw, rtol=1e-9)
    r.notes.append("The SI references of spec 10 section 3.11 are printed to 0.1 kN: the computed 19682.569, "
                   f"{s.barda:.3f} and {s.gw:.3f} kN round to the printed 19682.6, 17912.0 and 10702.0 kN, but the "
                   "rounding itself (1.4e-6 to 2.5e-6) exceeds the 1e-6 tolerance; the SI results agree with the "
                   "British results converted with 4.44822162 kN/kip to 7.5e-7 (the SI input data are themselves "
                   "rounded: f'c = 34473.8 for 34473.786 kN/m2, N_U = 4448.22 for 4448.2216 kN).")
    return r


# ======================================================================================
# VP-NON1
# ======================================================================================
NON1 = dict(k=1000.0, m=1.0, xi_el=0.02, xy=0.01, fy=10.0, hard=0.1, A=2.5, dt=0.005, nft=4096, n0=164)


def non1_backbone() -> HY.Backbone:
    p = NON1
    x2 = 20.0 * p["xy"]
    return HY.Backbone([p["xy"], x2], [p["fy"], p["fy"] + p["hard"] * p["k"] * (x2 - p["xy"])], 1)


def _bilinear_force(x: float) -> float:
    """Closed form of the bilinear BBC (independent of sassi.core.hysteresis)."""
    p = NON1
    if x <= p["xy"]:
        return p["k"] * x
    return p["fy"] + p["hard"] * p["k"] * (x - p["xy"])


def _bilinear_masing_damping(x: float) -> float:
    """Loop-area damping of the Masing loop on the bilinear BBC: E_D = 8 int_0^x F - 4 x F(x)."""
    p = NON1
    if x <= p["xy"]:
        return 0.0
    integ = 0.5 * p["fy"] * p["xy"] + p["fy"] * (x - p["xy"]) + 0.5 * p["hard"] * p["k"] * (x - p["xy"]) ** 2
    f = _bilinear_force(x)
    return (8.0 * integ - 4.0 * x * f) / (2.0 * math.pi * x * f)


def _complex_modulus(xi: float) -> complex:
    """SASSI complex modulus factor 1 - 2 xi^2 + 2 i xi sqrt(1 - xi^2) (D-CNV-03), written out."""
    return complex(1.0 - 2.0 * xi * xi, 2.0 * xi * math.sqrt(1.0 - xi * xi))


def _steady_amplitude(k: float, xi: float, w: float) -> float:
    """Amplitude of the relative displacement of the equivalent-linear SDOF (complex stiffness k c(xi)) under
    the base acceleration A g sin(w t): m A g / |k c(xi) - m w^2|."""
    p = NON1
    return p["m"] * p["A"] * 9.81 / abs(k * _complex_modulus(xi) - p["m"] * w * w)


def _sampled_amplitude(k: float, xi: float) -> float:
    """max |x(t_n)| of the periodic steady state sampled at t_n = n dt, n = 0 .. NFFT-1 (what RELDISP writes):
    x(t) = Re{C e^(i w t)}, C = i m A g / (k c(xi) - m w^2) for a(t) = A sin(w t) (e^(+i w t) convention)."""
    p = NON1
    w = 2.0 * math.pi * p["n0"] / (p["dt"] * p["nft"])
    c = 1j * p["m"] * p["A"] * 9.81 / (k * _complex_modulus(xi) - p["m"] * w * w)
    t = np.arange(p["nft"]) * p["dt"]
    return float(np.max(np.abs((c * np.exp(1j * w * t)).real)))


def non1_closed_form_iteration(maxit: int = HY.MAX_ITERATIONS) -> List[Dict[str, float]]:
    """The fixed-point iteration of the equivalent-linear SDOF written out independently: sampled steady-state
    amplitude, secant and Masing damping of the bilinear BBC (EDF 1, elastic damping added, no cut-off),
    convergence test of D-NON-06."""
    p = NON1
    k, xi = p["k"], p["xi_el"]
    out = []
    for it in range(maxit + 1):
        x = _sampled_amplitude(k, xi)
        kn = _bilinear_force(x) / x if x > p["xy"] else p["k"]
        xin = _bilinear_masing_damping(x) + p["xi_el"]
        out.append(dict(k=k, xi=xi, x=x, k_new=kn, xi_new=xin))
        conv = abs(kn - k) / k < HY.TOL_E and abs(xin - xi) < HY.TOL_XI
        k, xi = kn, xin
        if conv:
            break
    return out


@problem("VP-NON1", "Option NON SDOF: nonlinear spring on a rigid base under harmonic input, iterated through "
         "HOUSE-ANALYS-MOTION-RELDISP-NONLINEAR", tier="P2",
         modules=["HOUSE", "ANALYS", "MOTION", "RELDISP", "NONLINEAR"], source="requirements 4.15, D-NON-04, D-NON-06")
def vp_non1(workdir: Path) -> VPResult:
    """A 1 t mass on a spring (k_el = 1000 kN/m, bilinear GMR backbone: yield 10 kN at 0.01 m, 10 % hardening,
    2 % elastic damping) on one interaction node of a practically rigid site (as VP-01), driven by a harmonic
    control acceleration 2.5 g at f0 = 8.008 Hz (above the elastic frequency 5.03 Hz) with an integer number of
    cycles in the NFFT window: MOTION/RELDISP give the exact periodic steady state.  EDF = 1 (harmonic
    motion), damping = hysteretic + elastic (no cut-off)."""
    r = VPResult()
    p = NON1
    wd = Path(workdir) / "chain"
    shutil.rmtree(wd, ignore_errors=True)
    wd.mkdir(parents=True)
    dt, nft, n0 = p["dt"], p["nft"], p["n0"]
    df = 1.0 / (dt * nft)
    f0 = n0 * df
    w = 2.0 * math.pi * f0
    site = B.layered_site([(10.0, 1.0e5, 2.0e5, 2.0, 0.01)], (1.0e5, 2.0e5, 2.0, 0.01))
    fs = B.FrequencySet.fourier(dt, nft, [n0 - 1, n0, n0 + 1])
    B.run_soil(wd, "m", site, fs, layer=0, rad=1.0)
    mdl = B.sdof_on_node(site, p["k"], p["m"], p["xi_el"])
    base, mass = mdl["base"], mdl["mass"]
    B.run_house(wd, "m", mdl)
    B.run_analys(wd, "m", fs, save=1)
    t = np.arange(nft) * dt
    B.write_history(wd / "harm.acc", p["A"] * np.sin(w * t), dt)
    bb = non1_backbone()
    e = NLM.EqlDeck()
    e.params.update(model="m", edf=1.0, nonlinopts=NLM.OPT_SPRINGS, elasticd=1, gravity=site.gravity)
    for k, (xx, yy) in enumerate(zip(bb.x, bb.y), start=1):
        e.add("bbc", **{"bbc": 1, "type": 4, "yield": 1, "point": k, "x": float(xx), "y": float(yy)})
    e.add("springs", num=1, group=1, elem=1, bbc=1, disp=1, force=4, prop=1, node_i=base, node_j=mass, scx=p["k"],
          scy=0.0, scz=0.0, scxx=0.0, scyy=0.0, sczz=0.0, damp=p["xi_el"], k_el=p["k"])
    NLM.write_eql(wd / "m.eql", e)

    def response() -> None:
        B.run_motion(wd, "m", fs, "harm.acc", [(base, 1, 1, 0, 0, 0, 0, 0), (mass, 1, 1, 0, 0, 0, 0, 0)], out=1,
                     cplx=1, type=0, cm=0, ang=0.0, gravity=site.gravity)
        d = decks.new("RELDISP")
        d["model"], d["delt"], d["nft"], d["df"] = "m", dt, nft, df
        d["thfile"], d["mult"], d["max"], d["gravity"] = "harm.acc", 1.0, 0.0, site.gravity
        d["relfile"], d["reldisoutput"] = "", 1
        d.table("rdnd").append([base, 1, 0, 0, 0, 0, 0])
        d.table("rdnd").append([mass, 1, 0, 0, 0, 0, 0])
        B.write_deck(wd, "m", d)
        B.run("RELDISP", wd, "m")

    ref = non1_closed_form_iteration()
    history = []
    k_used, xi_used = p["k"], p["xi_el"]
    for it in range(HY.MAX_ITERATIONS + 1):
        response()
        rc = NLM.run_nonlinear("m", wd)
        if rc != 0:
            r.require(f"NONLINEAR run of SSI analysis {it} finishes with status OK", False,
                      (wd / "m_NONLINEAR.out").read_text(encoding="utf-8")[-1500:])
            return r
        _, rows = NLM.read_props(wd / "SPRING_EQL_Matl_Prop.txt")
        conv = NLM.read_convergence(wd)[-1]
        history.append(dict(k=k_used, xi=xi_used, x=rows[1].x_max, k_new=rows[1].k_new, xi_new=rows[1].xi_new,
                            converged=conv["converged"]))
        if conv["converged"]:
            break
        shutil.copyfile(wd / "m_new.hou", wd / "m.hou")
        B.run("HOUSE", wd, "m")
        B.run_analys(wd, "m", fs, mode=1)
        k_used, xi_used = rows[1].k_new, rows[1].xi_new
    n_it = len(history) - 1
    r.require(f"converged within {HY.MAX_ITERATIONS} iterations (D-NON-06): {n_it} iterations",
              bool(history[-1]["converged"]) and n_it <= HY.MAX_ITERATIONS)
    r.check("number of iterations equals the closed-form fixed-point iteration", n_it, len(ref) - 1, atol=0)
    # every SSI response is the closed-form steady state of the SDOF with the properties it was computed with
    worst = 0.0
    for h in history:
        worst = worse(worst, abs(h["x"] / _sampled_amplitude(h["k"], h["xi"]) - 1.0))
    r.check("max |x_SSI / x_closed-form - 1| over the iterations (sampled steady-state amplitude)", worst, 0.0,
            atol=1e-6)
    for i, (h, g) in enumerate(zip(history, ref)):
        r.check(f"iteration {i}: chain amplitude vs the independent fixed-point iteration", h["x"], g["x"], rtol=1e-6)
        r.check(f"iteration {i}: chain k_new vs the independent fixed-point iteration", h["k_new"], g["k_new"],
                rtol=1e-6)
    last = history[-1]
    xc = last["x"]
    r.check("converged: k_new = secant of the BBC at the converged amplitude, F(x)/x", last["k_new"],
            _bilinear_force(xc) / xc, rtol=1e-6)
    r.check("converged: xi_new = Masing loop-area damping at x + elastic damping", last["xi_new"],
            _bilinear_masing_damping(xc) + p["xi_el"], rtol=1e-6)
    r.check("converged: relative change of k between the last two iterations < 2 %",
            abs(last["k_new"] - last["k"]) / last["k"], 0.0, atol=HY.TOL_E)
    r.check("converged: change of damping between the last two iterations < 0.5 %",
            abs(last["xi_new"] - last["xi"]), 0.0, atol=HY.TOL_XI)
    # the time-domain hysteresis of the last iteration: steady loop on the BBC with the same energy
    x_t, _ = NLM.read_history(wd / "SPRING0001.thd")
    vmax = float(np.max(np.abs(x_t)))
    loop = HY.cyclic_loop(HY.MasingGMR(bb), vmax)
    r.check("time domain: stabilised loop damping at max|x| vs closed form", loop.damping,
            _bilinear_masing_damping(vmax), rtol=1e-9)
    r.inform("sampled peak max|x(t_n)| vs the continuous steady-state amplitude of the converged SDOF", xc,
             _steady_amplitude(last["k"], last["xi"], w),
             note="RELDISP samples the harmonic every dt: the sampled peak is below the amplitude by < (pi/1024)^2/2")
    r.notes.append(f"f0 = {f0:.4f} Hz, elastic frequency {math.sqrt(p['k'] / p['m']) / (2 * math.pi):.3f} Hz; "
                   f"iterations {n_it}: amplitudes " + ", ".join(f"{h['x']:.5f}" for h in history)
                   + "; k/k_el " + ", ".join(f"{h['k_new'] / p['k']:.4f}" for h in history)
                   + "; damping " + ", ".join(f"{h['xi_new']:.4f}" for h in history) + ".")
    return r
