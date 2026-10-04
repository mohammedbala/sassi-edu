"""Verification problems of the SOIL module (requirements section 6.3, R2 sections A.2 and H.1).

* VP-02a  uniform damped layer on an elastic half-space: SOIL transfer functions against the
  exact closed forms (both complex-modulus forms), sub-division invariance.
* VP-04   SHAKE91 sample problem (Idriss & Sun 1992): INP.DAT and DIAM.ACC are parsed here, the
  SOIL deck is built from them, SOIL is run and its output files are compared with Table B-2
  of the SHAKE91 manual and the SHAKE16 output (R2 H.1).
"""
from __future__ import annotations

import math

import re
import shutil
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from scipy import signal as _sig

from ...conventions import cfactor
from ...core import shake as SH
from ...core import spectra as SP
from ...io import decks, textfiles
from ...modules.base import run_module
from ...modules.soil import saf_tf_name, soil_rs_name, write_deck
from .. import VPResult, problem

_PKG_DATA = Path(__file__).resolve().parents[2] / "data" / "shake91_example"
_DOC_DATA = Path(__file__).resolve().parents[3] / "docs" / "spec" / "R2_benchmark_data" / "shake91_example"


def shake91_dir() -> Path:
    """Directory holding INP.DAT / DIAM.ACC (package data first, then docs/spec)."""
    for d in (_PKG_DATA, _DOC_DATA):
        if (d / "INP.DAT").exists() and (d / "DIAM.ACC").exists():
            return d
    raise FileNotFoundError("SHAKE91 sample data (INP.DAT, DIAM.ACC) not found")


# =======================================================================================
# VP-02a  uniform layer, closed form (R2 A.1/A.2)
# =======================================================================================
VP02A_FREQS = [0.5, 1.0, 1.5, 1.6667, 2.0, 3.0, 5.0, 8.3333]
VP02A_REF = {   # (form) -> (|surf / base-within|, |surf / rock outcrop|), R2 A.2 table (4 decimals)
    0: ([1.1216, 1.6931, 5.7710, 12.7152, 3.1156, 1.0411, 4.2038, 2.4815],
        [1.1150, 1.6207, 3.4116, 3.8327, 2.4115, 1.0087, 2.3554, 1.6702]),
    1: ([1.1209, 1.6878, 5.6735, 12.7633, 3.1590, 1.0436, 4.2202, 2.4918],
        [1.1143, 1.6161, 3.3916, 3.8344, 2.4313, 1.0110, 2.3610, 1.6756]),
}
VP02A_LAYER = dict(H=30.0, vs=200.0, rho=2.0, xi=0.05, vr=1000.0, rhor=2.2, xir=0.01)


def uniform_layer_closed_form(f, form: int = 0, H=30.0, vs=200.0, rho=2.0, xi=0.05, vr=1000.0, rhor=2.2,
                              xir=0.01):
    """Exact TFs of a uniform layer on an elastic half-space (R2 A.1).

    Returns (surface / base-within, surface / rock-outcrop) = (1/cos k*H, 1/(cos k*H + i a* sin k*H)).
    """
    w = 2 * np.pi * np.asarray(f, dtype=float)
    vss = vs * np.sqrt(cfactor(xi, form))
    vrs = vr * np.sqrt(cfactor(xir, form))
    kH = w / vss * H
    a = rho * vss / (rhor * vrs)
    return 1.0 / np.cos(kH), 1.0 / (np.cos(kH) + 1j * a * np.sin(kH))


def _synthetic_motion(path: Path, n: int, dt: float, seed: int = 7) -> None:
    """A tapered random acceleration record (g) for runs whose TFs do not depend on the motion."""
    rng = np.random.default_rng(seed)
    t = np.arange(n) * dt
    env = np.sin(np.pi * t / (n * dt)) ** 2
    a = 0.05 * env * rng.standard_normal(n)
    path.write_text("\n".join(f"{v:.8e}" for v in a) + "\n", encoding="utf-8")


def _vp02a_deck(workdir: Path, model: str, thicknesses: List[float], form: int, fstep: float, cutoff: float) -> Path:
    p = VP02A_LAYER
    grav = 9.81
    d = decks.new("SOIL")
    d.params.update(dict(title=f"VP-02a uniform layer (form {form})", model=model, nrval=512, grav=grav,
                         header=0, outcrop=1, save=1, iter=0, ratio=0.65, gravmult=1.0, cof=0.0,
                         soilcutoff=cutoff, delt=0.01, nft=1024, cl=len(thicknesses) + 1, thfile="vp02a.acc",
                         thtit="synthetic", mult=0.0, max=0.1, indir=0, cmodform=form))
    for i, h in enumerate(thicknesses, start=1):
        d.table("profile").append([i, h, p["rho"] * grav, 2.0 * p["vs"], p["vs"], p["xi"], p["xi"], ""])
    d.table("profile").append([len(thicknesses) + 1, 0.0, p["rhor"] * grav, 2.0 * p["vr"], p["vr"], p["xir"],
                               p["xir"], ""])
    nb = len(thicknesses) + 1
    d.table("ssaf").append([1, 1, 1, 0, nb, fstep, "surface / base within"])
    d.table("ssaf").append([1, 1, 1, 1, nb, fstep, "surface / rock outcrop"])
    path = workdir / f"{model}.soi"
    write_deck(path, d)
    return path


def _read_saf(path: Path, fstep: float):
    f, H, _ = textfiles.read_tf(path)
    k = np.rint(f / fstep).astype(int)
    return k, H


@problem("VP-02a", "Uniform layer TFs, closed form (SOIL, 0 iterations)", tier="P0", modules=["SOIL"],
         source="R2 A.1, A.2")
def vp02a(workdir: Path) -> VPResult:
    """SOIL transfer functions of a uniform damped layer (H = 30 m, Vs = 200 m/s, rho = 2.0, 5 %)
    on an elastic half-space (1000 m/s, 2.2, 1 %): module output vs exact closed form (1e-8),
    closed form vs the R2 table (printed precision), any sub-division identical."""
    r = VPResult()
    workdir = Path(workdir)
    fstep, cutoff = 1e-4, 10.0
    _synthetic_motion(workdir / "vp02a.acc", 512, 0.01)
    kref = np.rint(np.asarray(VP02A_FREQS) / fstep).astype(int)
    fgrid = kref * fstep
    nb_cases = {"1 sublayer": [30.0], "5 sublayers": [3.0, 7.0, 10.0, 4.0, 6.0]}
    for form in (0, 1):
        tag = "SASSI form" if form == 0 else "CMODFORM,1 (1+2i beta)"
        hw, ho = uniform_layer_closed_form(fgrid, form)
        # closed form vs the tabulated R2 values (4 decimals -> half a unit of the last digit)
        for j, f in enumerate(VP02A_FREQS):
            r.check(f"closed form |surf/base-within| {f} Hz, {tag} vs R2 table", abs(hw[j]), VP02A_REF[form][0][j],
                    atol=5e-5, note="R2 values are printed to 4 decimals")
            r.check(f"closed form |surf/outcrop| {f} Hz, {tag} vs R2 table", abs(ho[j]), VP02A_REF[form][1][j],
                    atol=5e-5, note="R2 values are printed to 4 decimals")
        results = {}
        for case, th in nb_cases.items():
            model = f"vp02a_f{form}_{len(th)}"
            _vp02a_deck(workdir, model, th, form, fstep, cutoff)
            rc = run_module("SOIL", model, workdir)
            if not r.require(f"SOIL run ({case}, {tag}) status 0", rc == 0):
                r.notes.append((workdir / f"{model}_SOIL.out").read_text()[-2000:])
                continue
            nb = len(th) + 1
            kw, Hw = _read_saf(workdir / saf_tf_name(1, True, nb, False), fstep)
            ko, Ho = _read_saf(workdir / saf_tf_name(1, True, nb, True), fstep)
            results[case] = (Hw, Ho)
            iw = np.searchsorted(kw, kref)
            io = np.searchsorted(ko, kref)
            for j, f in enumerate(VP02A_FREQS):
                r.check(f"SOIL |surf/base-within| {f} Hz ({case}, {tag})", abs(Hw[iw[j]]), abs(hw[j]), rtol=1e-8)
                r.check(f"SOIL |surf/outcrop| {f} Hz ({case}, {tag})", abs(Ho[io[j]]), abs(ho[j]), rtol=1e-8)
            # the whole grid against the closed form (complex values)
            fall = kw * fstep
            cw, co = uniform_layer_closed_form(fall, form)
            r.check(f"max rel. error of complex surf/base-within on 0-{cutoff:g} Hz ({case}, {tag})",
                    float(np.max(np.abs(Hw - cw) / np.abs(cw))), 0.0, atol=1e-8)
            r.check(f"max rel. error of complex surf/outcrop on 0-{cutoff:g} Hz ({case}, {tag})",
                    float(np.max(np.abs(Ho - co) / np.abs(co))), 0.0, atol=1e-8)
        if len(results) == 2:
            (a1, b1), (a5, b5) = results["1 sublayer"], results["5 sublayers"]
            r.check(f"sub-division invariance, within ({tag})", float(np.max(np.abs(a5 - a1) / np.abs(a1))), 0.0,
                    atol=1e-8)
            r.check(f"sub-division invariance, outcrop ({tag})", float(np.max(np.abs(b5 - b1) / np.abs(b1))), 0.0,
                    atol=1e-8)
    r.notes.append("Module outputs are the SSAF Fourier amplification files (exact frequencies k*1e-4 Hz); "
                   "the R2 table columns correspond to f = 1.6667 and 8.3333 Hz exactly as printed.")
    return r


# =======================================================================================
# VP-04  SHAKE91 sample problem (R2 H.1)
# =======================================================================================
def _fixed(line: str, start: int, width: int) -> Optional[float]:
    s = line[start:start + width].strip()
    return float(s) if s else None


def _numbers_until(lines: List[str], i: int, n: int):
    vals: List[float] = []
    while len(vals) < n:
        vals.extend(float(t) for t in lines[i].split())
        i += 1
    return vals[:n], i


def parse_shake91_input(path: Path) -> Dict:
    """Parse the SHAKE91 input file (options 1-5 and the output requests 6, 9, 10).

    Returns a dict with ``materials`` (list of dicts with G/D curves), ``layers`` (type, thick,
    damping, unit weight, Vs; last = base), ``motion`` (nv, nfft, dt, file, fmt, mult, xmax, fmax,
    nhead, nperline), ``input_layer``, ``outcrop`` (True when SHAKE91 code 0), ``iterations``,
    ``ratio``, ``rs`` and ``amplification`` requests.
    """
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    out: Dict = {}
    i = 0
    n = len(lines)

    def is_int_line(s, val):
        t = s.split()
        return len(t) == 1 and t[0] == str(val)

    while i < n:
        ln = lines[i]
        t = ln.split()
        if len(t) == 1 and t[0].isdigit():
            opt = int(t[0])
            i += 1
            if opt == 0:
                break
            if opt == 1:
                nmat = int(lines[i].split()[0])
                i += 1
                mats = []
                for _m in range(nmat):
                    curves = []
                    for _c in range(2):
                        head = lines[i]
                        npts = int(head.split()[0])
                        title = head[5:].strip() if len(head) > 5 else ""
                        xs, i = _numbers_until(lines, i + 1, npts)
                        ys, i = _numbers_until(lines, i, npts)
                        curves.append((title, xs, ys))
                    mats.append({"g_title": curves[0][0], "g_strain": curves[0][1], "g_ratio": curves[0][2],
                                 "d_title": curves[1][0], "d_strain": curves[1][1], "d_pct": curves[1][2]})
                out["materials"] = mats
                i += 1      # material-number list line
            elif opt == 2:
                head = lines[i].split()
                nl = int(head[1])
                out["profile_title"] = " ".join(head[2:])
                i += 1
                layers = []
                for _k in range(nl):
                    s = lines[i]
                    layers.append({"no": int(s[0:5]), "type": int(s[5:10]), "thick": _fixed(s, 15, 10) or 0.0,
                                   "gmax": _fixed(s, 25, 10), "damp": _fixed(s, 35, 10),
                                   "weight": _fixed(s, 45, 10), "vs": _fixed(s, 55, 10)})
                    i += 1
                out["layers"] = layers
            elif opt == 3:
                a = lines[i].split()
                b = lines[i + 1].split()
                m = {"nv": int(a[0]), "nfft": int(a[1]), "dt": float(a[2]), "file": a[3],
                     "fmt": a[4] if len(a) > 4 else ""}
                if len(b) == 4:
                    m.update(mult=0.0, xmax=float(b[0]), fmax=float(b[1]), nhead=int(b[2]), nperline=int(b[3]))
                else:
                    m.update(mult=float(b[0]), xmax=float(b[1]), fmax=float(b[2]), nhead=int(b[3]),
                             nperline=int(b[4]))
                out["motion"] = m
                i += 2
            elif opt == 4:
                a = lines[i].split()
                out["input_layer"] = int(a[0])
                out["outcrop"] = int(a[1]) == 0          # SHAKE91: 0 outcropping, 1 within
                i += 1
            elif opt == 5:
                a = lines[i].split()
                out["iterations"] = int(a[1])
                out["ratio"] = float(a[2])
                i += 1
            elif opt == 6:
                lay = [int(x) for x in lines[i].split()]
                typ = [int(x) for x in lines[i + 1].split()]
                sav = [int(x) for x in lines[i + 2].split()]
                out.setdefault("acc_requests", []).extend(zip(lay, typ, sav))
                i += 3
            elif opt == 7:
                out.setdefault("stress_requests", []).append(lines[i].split()[:5])
                out["stress_requests"].append(lines[i + 1].split()[:5])
                i += 2
            elif opt == 9:
                a = lines[i].split()
                b = lines[i + 1].split()
                c = lines[i + 2].split()
                out["rs"] = {"layer": int(a[0]), "outcrop": int(a[1]) == 0, "gravity": float(b[2]),
                             "damping": [float(x) for x in c]}
                i += 3
            elif opt == 10:
                a = lines[i].split()
                out["amplification"] = {"layer_in": int(a[0]), "outcrop_in": int(a[1]) == 0,
                                        "layer_out": int(a[2]), "outcrop_out": int(a[3]) == 0, "df": float(a[4])}
                i += 1
            continue
        i += 1
    return out


def build_shake91_deck(inp: Dict, workdir: Path, model: str = "shake91", motion_file: str = "DIAM.ACC") -> Path:
    """Translate the parsed SHAKE91 input into a SOIL deck (D-SOL-02: the base is the last SPRO row)."""
    grav = 32.2
    m = inp["motion"]
    nl = len(inp["layers"])
    d = decks.new("SOIL")
    d.params.update(dict(title=f"VP-04 SHAKE91 sample: {inp.get('profile_title', '')}", model=model,
                         nrval=m["nv"], grav=grav, header=m["nhead"], outcrop=1 if inp["outcrop"] else 0, save=1,
                         iter=inp["iterations"], ratio=inp["ratio"], gravmult=1.0, cof=0.0, soilcutoff=m["fmax"],
                         delt=m["dt"], nft=m["nfft"], cl=inp["input_layer"], thfile=motion_file,
                         thtit="Loma Prieta 1989 Diamond Heights H1_90", mult=m["mult"], max=m["xmax"], indir=0,
                         cmodform=0))
    for k, mat in enumerate(inp["materials"], start=1):
        lab = f"MAT{k}"
        for x, y in zip(mat["g_strain"], mat["g_ratio"]):
            d.table("dynp").append([lab, "G", x, y])
        for x, y in zip(mat["d_strain"], mat["d_pct"]):
            d.table("dynp").append([lab, "D", x, y])
    for ly in inp["layers"]:
        vs = ly["vs"]
        d.table("profile").append([ly["no"], ly["thick"], ly["weight"], 2.0 * vs, vs, ly["damp"], ly["damp"],
                                   f"MAT{ly['type']}"])
    for lay in range(1, nl + 1):
        d.table("sacc").append([lay, 2, 1 if lay == 1 else 0])
    for lay in range(1, nl):
        d.table("sstr").append([lay, 1, 1, 1, 1])
    rs = inp.get("rs", {"layer": 1, "outcrop": True, "damping": [0.05]})
    d.table("srs").append([rs["layer"], 1, 1 if rs["outcrop"] else 0])
    for z in rs["damping"]:
        d.table("damp").append([z])
    amp = inp.get("amplification")
    if amp:   # SHAKE91 option 10: output / input; SSAF: <layer> numerator, <layer2> denominator
        d.table("ssaf").append([amp["layer_out"], 1, 1 if amp["outcrop_out"] else 0, 1 if amp["outcrop_in"] else 0,
                                amp["layer_in"], amp["df"], "surface / rock outcrop (SHAKE91 option 10)"])
    path = workdir / f"{model}.soi"
    write_deck(path, d)
    return path


# Table B-2 of the SHAKE91 manual (R2 H.1): sublayer -> (gamma_eff %, damping, G ksf, G/Gmax, gamma_max %, tau_max psf)
TABLE_B2 = {
    1: (0.00077, 0.007, 3851.5, 0.992, 0.00154, 59.43), 2: (0.00295, 0.014, 3020.0, 0.960, 0.00591, 178.41),
    3: (0.00634, 0.023, 2803.8, 0.892, 0.01267, 355.31), 4: (0.00976, 0.028, 2985.8, 0.852, 0.01952, 582.73),
    5: (0.01099, 0.030, 3621.7, 0.933, 0.02197, 795.85), 6: (0.01403, 0.035, 3540.5, 0.912, 0.02806, 993.46),
    7: (0.01362, 0.034, 4296.0, 0.915, 0.02723, 1169.83), 8: (0.01566, 0.037, 4239.8, 0.903, 0.03132, 1327.93),
    9: (0.01356, 0.034, 5402.8, 0.792, 0.02711, 1464.73), 10: (0.01505, 0.037, 5266.1, 0.772, 0.03011, 1585.43),
    11: (0.01336, 0.034, 6288.3, 0.795, 0.02671, 1679.79), 12: (0.01413, 0.035, 6203.6, 0.784, 0.02825, 1752.67),
    13: (0.01233, 0.032, 7357.0, 0.810, 0.02467, 1814.84), 14: (0.01282, 0.033, 7290.6, 0.803, 0.02563, 1868.58),
    15: (0.01115, 0.030, 8570.2, 0.829, 0.02230, 1911.02), 16: (0.00865, 0.026, 11292.4, 0.863, 0.01729, 1952.92),
}
# peak acceleration profile (g), Table B-2 option 6: depth ft -> PGA (top of sublayer, within except surface)
PGA_PROFILE = {0: 0.19037, 5: 0.19006, 10: 0.18876, 20: 0.18258, 30: 0.17208, 40: 0.15947, 50: 0.14288,
               60: 0.12652, 70: 0.11050, 80: 0.09840, 90: 0.08999, 100: 0.08268, 110: 0.08559, 120: 0.08547,
               130: 0.08198, 140: 0.07769, 150: 0.07617}
OUTCROP_TF = {1.0: 1.3270, 2.0: 3.3181, 2.125: 3.4369, 3.0: 1.7386, 5.0: 2.0450, 10.0: 1.6690}
SA_REF = {0.10: 0.3977, 0.40: 0.7832, 1.0: 0.1351}


@problem("VP-04", "SHAKE91 sample problem (SOIL complete)", tier="P0", modules=["SOIL"], source="R2 H.1")
def vp04(workdir: Path) -> VPResult:
    """SHAKE91 manual example: 16 sublayers + rock, outcrop input 0.10 g at the base, 8 iterations,
    strain ratio 0.50, cut-off 25 Hz.  Tolerances (requirements 6.3): PGA 0.5 %; strain, G, D 1 %;
    SA 2 %.  Amplification values (no tolerance stated) use 1 %."""
    r = VPResult()
    workdir = Path(workdir)
    src = shake91_dir()
    inp = parse_shake91_input(src / "INP.DAT")
    shutil.copy(src / inp["motion"]["file"], workdir / inp["motion"]["file"])
    model = "shake91"
    build_shake91_deck(inp, workdir, model, inp["motion"]["file"])
    rc = run_module("SOIL", model, workdir)
    listing = (workdir / f"{model}_SOIL.out").read_text(encoding="utf-8")
    if not r.require("SOIL run status 0", rc == 0):
        r.notes.append(listing[-3000:])
        return r
    nl = len(inp["layers"])
    depth_top = np.concatenate([[0.0], np.cumsum([ly["thick"] for ly in inp["layers"]][:-1])])

    # ---- peak accelerations (0.5 %)
    for lay in range(1, nl + 1):
        a, dt = textfiles.read_history(workdir / f"ACC{lay:03d}.TH")
        dep = int(round(depth_top[lay - 1]))
        if dep in PGA_PROFILE:
            what = "surface (outcrop)" if lay == 1 else ("base within" if lay == nl else f"{dep} ft within")
            r.check(f"PGA {what} (g)", float(np.max(np.abs(a))), PGA_PROFILE[dep], rtol=5e-3)
    # ---- strain-compatible properties (FILE88) and maxima (TH files)
    f88 = SH.read_file88(workdir / "FILE88")
    rho = np.array([ly["weight"] for ly in inp["layers"]]) / 32.2
    gmax = rho * np.array([ly["vs"] for ly in inp["layers"]]) ** 2
    mats = inp["materials"]
    for i in range(nl - 1):
        lay = i + 1
        geff, D, G, GG, gmx, tau = TABLE_B2[lay]
        r.check(f"sublayer {lay}: effective strain (%)", f88["gamma_eff_pct"][i], geff, rtol=1e-2)
        r.check(f"sublayer {lay}: G (ksf)", f88["G"][i], G, rtol=1e-2)
        r.check(f"sublayer {lay}: G/Gmax", f88["G"][i] / gmax[i], GG, rtol=1e-2)
        r.check(f"sublayer {lay}: damping vs Table B-2 (3 decimals)", f88["beta_s"][i], D, atol=5e-4,
                note="reference printed to 3 decimals: tolerance = half a unit of the last digit")
        mat = mats[inp["layers"][i]["type"] - 1]
        # the SHAKE91 rule written out (linear in log10 strain, constant outside the table), independent of
        # sassi.core.shake.interp_log_strain, which SOIL itself uses (final audit)
        d_ref = float(np.interp(math.log10(geff), np.log10(np.asarray(mat["d_strain"], float)),
                                np.asarray(mat["d_pct"], float))) / 100.0
        r.check(f"sublayer {lay}: damping vs D(published effective strain)", f88["beta_s"][i], d_ref, rtol=1e-2)
        g_strain, _ = textfiles.read_history(workdir / f"SN{lay:03d}.TH")
        r.check(f"sublayer {lay}: max strain (%)", float(np.max(np.abs(g_strain))), gmx, rtol=1e-2)
        s_hist, _ = textfiles.read_history(workdir / f"SS{lay:03d}.TH")
        r.check(f"sublayer {lay}: max stress (psf)", 1000.0 * float(np.max(np.abs(s_hist))), tau, rtol=1e-2,
                note="SOIL: tau = G* gamma; SHAKE91 prints G gamma(t)")

    # ---- amplification surface / rock outcrop (SSAF file, df = 0.125 Hz)
    amp = inp["amplification"]
    fn = saf_tf_name(amp["layer_out"], amp["outcrop_out"], amp["layer_in"], amp["outcrop_in"])
    f, H, _ = textfiles.read_tf(workdir / fn)
    for fr, ref in OUTCROP_TF.items():
        j = int(np.argmin(np.abs(f - fr)))
        r.check(f"|surface/rock outcrop| at {fr} Hz", abs(H[j]), ref, rtol=1e-2)
    j = int(np.argmax(np.abs(H)))
    r.check("max outcrop amplification", abs(H[j]), 3.44, rtol=1e-2)
    r.check("frequency of max outcrop amplification (Hz)", f[j], 2.12, rtol=1e-2)

    # ---- final surface / within-base amplification (listing, FFT grid as SHAKE91)
    mm = re.search(r"surface / base WITHIN \(final\) = ([0-9.]+) at ([0-9.]+) Hz", listing)
    if r.require("listing reports the final surface/within-base amplification", mm is not None):
        r.check("final max amplification surface/within base", float(mm.group(1)), 20.47, rtol=1e-2)
        r.check("frequency of final max amplification (Hz)", float(mm.group(2)), 2.11, rtol=1e-2)

    # ---- response spectrum of the surface motion (5 %)
    a1, dt = textfiles.read_history(workdir / "ACC001.TH")
    T = np.array(sorted(SA_REF))
    rs = SP.response_spectrum(a1, dt, 1.0 / T, [0.05])
    # Lead decision (wave-1 review): the printed SHAKE91 spectral values are spectra of a record
    # corrupted by SHAKE91's DRCTSP unit conversion (see the diagnostic below and the note), so
    # they are not compared with the spectrum of SOIL's surface motion.  The pass criterion is
    # that applying the DRCTSP defect to SOIL's surface motion reproduces the printed values.
    rsf = textfiles.read_xy(workdir / soil_rs_name(1, 1))
    for fr in (1.0, 10.0):
        j = int(np.argmin(np.abs(rsf[:, 0] - fr)))
        k = int(np.argmin(np.abs(1.0 / T - fr)))
        r.check(f"RS file value at {fr:g} Hz = spectrum of ACC001.TH", rsf[j, 1], rs["SA"][0, k], rtol=1e-6)
    # ---- diagnostic: the published spectrum is that of a record corrupted by SHAKE91 DRCTSP
    ab = shake91_drctsp_record(a1[:int(inp["motion"]["nv"] * 1.1)], inp["rs"]["gravity"])
    rsb = SP.response_spectrum(ab, dt, 1.0 / T, [0.05], upsample=False)
    for k, Tk in enumerate(T):
        r.check(f"SA at T = {Tk:.2f} s of ACC001.TH with the SHAKE91 DRCTSP scaling defect (g) (D-W1-06)",
                rsb["SA"][0, k] / inp["rs"]["gravity"], SA_REF[Tk], rtol=2e-2,
                note="SHAKE91 printed value reproduced from SOIL's surface motion with the DRCTSP defect emulated")
    # converged spectrum of the SOIL surface motion: band-limited FFT up-sampling (x32) before the
    # exact Nigam-Jennings recurrence; the lead's checker up-samples only above 0.1/dt (5 Hz here)
    up = 32
    a_up = _sig.resample(a1, up * len(a1))
    sa_conv = [float(np.max(np.abs(SP.sdof_response(a_up, dt / up, 1.0 / Tk, 0.05)[2]))) for Tk in T]
    sa_chk = [float(v) for v in rs["SA"][0]]
    conv_txt = ", ".join(f"{v:.4f} g (T = {Tk:g} s)" for v, Tk in zip(sa_conv, T))
    chk_txt = "/".join(f"{v:.4f}" for v in sa_chk)
    dev_txt = ", ".join(f"{100 * (v / SA_REF[Tk] - 1):+.0f} %" for v, Tk in zip(sa_conv, T))
    r.notes.append(
        "The published SHAKE91/SHAKE16 spectral accelerations are not spectra of the published motion, and the "
        "cause is identified: SHAKE91 subroutine DRCTSP converts the record from g to cm/s^2 with "
        "'if (zmax .gt. ABS(A(K))) then A(K) = GGT*A(K) else zmax = ABS(A(K))', so every sample that sets a new "
        "running maximum (including the surface PGA at 11.28 s) is left in g, i.e. ~981 times too small, and "
        "the spectrum is computed for that corrupted record (the loop also never computes the last period, "
        "hence SA(10 s) = 0 in the printout).  SOIL reproduces the SHAKE91 motion itself (PGA profile, "
        "strains, moduli, damping, transfer functions, mean-square frequencies 2.42/2.52 Hz) to the printed "
        f"digits; the converged spectrum of that motion (x{up} FFT up-sampling, Nigam-Jennings) is {conv_txt} "
        f"vs published 0.3977, 0.7832, 0.1351 g ({dev_txt}); sassi.core.spectra.response_spectrum with its default "
        f"sampling (up-sampling only above 0.1/dt = {0.1 / dt:g} Hz) gives {chk_txt} g.  "
        "Applying the DRCTSP defect to ACC001.TH reproduces all 151 "
        "printed values within 0.6 %.  The published SA references (R2 H.1 / req. 6.3) "
        "are therefore not valid spectra of the published motion; by lead decision D-W1-06 the VP checks the "
        "three SA values of the requirement table on the emulated record (2 %) and reports the correct spectrum "
        "of the SOIL surface motion here.")
    return r


def shake91_drctsp_record(acc_g: np.ndarray, ggt: float) -> np.ndarray:
    """Reproduce the unit conversion of SHAKE91 subroutine DRCTSP (diagnostic only).

    DRCTSP multiplies a sample by GGT only when it is below the running maximum of |a|; samples
    that set a new running maximum stay in g units.  Returns the record (GGT units) that SHAKE91
    actually integrates for its response spectra.
    """
    out = np.asarray(acc_g, dtype=float).copy()
    zmax = 0.0
    for k in range(len(out)):
        if zmax > abs(out[k]):
            out[k] = ggt * out[k]
        else:
            zmax = abs(out[k])
    return out
