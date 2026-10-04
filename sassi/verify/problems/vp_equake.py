"""VP-36: RG 1.60 spectrum-compatible motions from EQUAKE (requirements 6.3; R2 G.2, G.3).

Two horizontal components and one vertical component are generated for the RG 1.60 spectra
(5 %, anchored to 1.0 g, dt = 0.005 s, 20 s) with the SRP 3.7.1 Appendix A minimum PSD as the
target PSD of the horizontal components.  The output files are read back and checked against
the SRP 3.7.1 Rev. 4 Option 1 Approach 2 criteria, the manual's criteria, strong-motion
duration, statistical independence (|rho| <= 0.16), the PSD requirement, compatibility of the
.acc/.vel/.dis files and the absence of baseline drift (displacement energy at periods longer
than the record < 50 %; PGD not above the target's spectral displacement at 0.1 Hz).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ...core import equake_lib as EL
from ...core import spectra as SP
from ...io import decks, textfiles
from ...modules.base import run_module
from .. import VPResult, problem

# Frequencies of the target files: all RG 1.60 control points (0.25, 2.5, 3.5, 9, 33 Hz) are
# included, so log-log interpolation of the file reproduces the tripartite lines exactly.
RG160_FILE_FREQS = [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0,
                    7.0, 8.0, 9.0, 12.0, 15.0, 20.0, 25.0, 33.0, 50.0, 100.0]
PSD_FILE_FREQS = sorted(set(np.round(np.logspace(np.log10(0.05), np.log10(100.0), 121), 6).tolist()
                            + [2.5, 9.0, 16.0]))


def write_rg160_inputs(workdir: Path, gravity: float = 9.80665):
    """Write RG 1.60 targets (H.rsi, V.rsi) and the App. A target PSD (H.tpsd) into ``workdir``."""
    units = "BS" if gravity > 20 else "SI"
    f = np.asarray(RG160_FILE_FREQS)
    # two numeric columns only (no header line): the record count must equal <nrfreq> (Error 89)
    EL.write_spectrum_file(workdir / "RG160H.rsi", f, EL.rg160_spectrum(f, "H", 0.05, 1.0))
    EL.write_spectrum_file(workdir / "RG160V.rsi", f, EL.rg160_spectrum(f, "V", 0.05, 1.0))
    fp = np.asarray(PSD_FILE_FREQS)
    textfiles.write_xy(workdir / "RG160H.tpsd", fp, EL.rg160_target_psd(fp, 1.0, units))


@problem("VP-36", "Spectrum-compatible motion (RG 1.60), SRP 3.7.1 criteria", tier="P0", modules=["EQUAKE"],
         source="R2 G.2, G.3", slow=True)
def vp36(workdir: Path) -> VPResult:
    r = VPResult()
    workdir = Path(workdir)
    # ---- target generator: RG 1.60 control points (R2 G.3)
    H = lambda f: float(EL.rg160_spectrum(np.array([f]), "H")[0])  # noqa: E731
    V = lambda f: float(EL.rg160_spectrum(np.array([f]), "V")[0])  # noqa: E731
    r.check("RG 1.60 H: A at 33 Hz (g)", H(33.0), 1.0, rtol=1e-12)
    r.check("RG 1.60 H: B at 9 Hz (g)", H(9.0), 2.61, rtol=1e-12)
    r.check("RG 1.60 H: C at 2.5 Hz (g)", H(2.5), 3.13, rtol=1e-12)
    r.check("RG 1.60 H: D at 0.25 Hz, SA from 2.05 x 36 in (g)", H(0.25), 0.4716, atol=5e-5,
            note="R2 value printed to 4 decimals")
    r.check("RG 1.60 V: C at 3.5 Hz (g)", V(3.5), 2.98, rtol=1e-12)
    r.check("RG 1.60 V: D at 0.25 Hz, SA from 1.37 x 36 in (g)", V(0.25), 0.3152, atol=5e-5,
            note="R2 value printed to 4 decimals")
    r.check("RG 1.60 H constant displacement below 0.25 Hz: SD(0.1 Hz) (in)",
            H(0.1) * EL.G_IN_S2 / (2 * np.pi * 0.1) ** 2, 2.05 * 36.0, rtol=1e-12)
    r.check("App. A PSD at 2.5 Hz (cm^2/s^3)", float(EL.rg160_target_psd(np.array([2.5]))[0]), 4190.0, rtol=1e-12)

    # ---- EQUAKE run
    gravity = 9.80665
    g_len = gravity * 100.0
    write_rg160_inputs(workdir, gravity)
    for name in ("RG160H.rsi", "RG160V.rsi"):
        tg = EL.TargetSpectrum.from_file(workdir / name)
        comp = "H" if "H" in name else "V"
        fchk = EL.check_grid(0.1, 100.0)
        err = float(np.max(np.abs(tg(fchk) / EL.rg160_spectrum(fchk, comp) - 1.0)))
        r.check(f"{name}: log-log interpolation reproduces RG 1.60 (max rel. error, 0.1-100 Hz)", err, 0.0,
                atol=1e-8, note="file values written with 9 significant digits")
    d = decks.new("EQUAKE")
    d.params.update(dict(title="VP-36 RG 1.60", model="vp36", opmode=0, accopt=0, nrfreq=len(RG160_FILE_FREQS),
                         rand=11975, damp=0.05, dur=20.0, corr=0, seeds=1, tpsd=1, delt=0.005, gravity=gravity,
                         eqtit="RG 1.60, 5 %, 1.0 g"))
    d.table("spectra").append([1, "RG160H.rsi", "H1.rso", "", "H1.acc", "RG160H.tpsd"])
    d.table("spectra").append([2, "RG160H.rsi", "H2.rso", "", "H2.acc", "RG160H.tpsd"])
    d.table("spectra").append([3, "RG160V.rsi", "V.rso", "", "V.acc", ""])
    decks.write(workdir / "vp36.equ", d)
    rc = run_module("EQUAKE", "vp36", workdir)
    listing = (workdir / "vp36_EQUAKE.out").read_text(encoding="utf-8")
    if not r.require("EQUAKE run status 0", rc == 0):
        r.notes.append(listing[-3000:])
        return r

    accs = {}
    for name, tfile, psd in (("H1", "RG160H.rsi", True), ("H2", "RG160H.rsi", True), ("V", "RG160V.rsi", False)):
        a, dt = textfiles.read_history(workdir / f"{name}.acc")
        accs[name] = a
        tg = EL.TargetSpectrum.from_file(workdir / tfile)
        r.require(f"{name}: time step <= 0.005 s and Nyquist >= 50 Hz (dt = {dt:g} s)", dt <= 0.005 + 1e-12)
        r.require(f"{name}: total duration >= 20 s ({(len(a) - 1) * dt:.3f} s)", (len(a) - 1) * dt >= 20.0 - 1e-9)
        checks = EL.evaluate_rs(a, dt, tg, 0.05)
        for key, c in checks.items():
            lab = "manual" if key == "manual" else "SRP Rev.4"
            r.require(f"{name} {lab}: >= 100 points/decade over {c.fmin:g}-{c.fmax:g} Hz ({c.per_decade:.1f})",
                      c.ok_density)
            r.require(f"{name} {lab}: no point > 10 % below (min SA/target {c.min_ratio:.3f})", c.ok_below)
            r.require(f"{name} {lab}: no point > 30 % above (max SA/target {c.max_ratio:.3f})", c.ok_above)
            r.require(f"{name} {lab}: <= 9 adjacent points below ({c.max_run_below})", c.ok_run)
        srp = checks["srp"]
        r.require(f"{name}: SRP band covers 0.1-50 Hz ({srp.fmin:g}-{srp.fmax:g})",
                  srp.fmin <= 0.1 + 1e-9 and srp.fmax >= 50.0 - 1e-9)
        pm = EL.ground_motion_parameters(a, dt, g_len)
        r.require(f"{name}: strong-motion duration >= 6 s ({pm['D5-75']:.2f} s)", pm["D5-75"] >= 6.0)
        # baseline: compatible, no drift (files .vel/.dis vs integration of .acc)
        vel, _ = textfiles.read_history(workdir / f"{name}.vel")
        dis, _ = textfiles.read_history(workdir / f"{name}.dis")
        v_i, d_i = SP.integrate(a * g_len, dt)
        r.check(f"{name}: .vel = integral of .acc (max diff / PGV)", float(np.max(np.abs(vel - v_i)) / pm["PGV"]),
                0.0, atol=1e-6)
        r.check(f"{name}: .dis = integral of .vel (max diff / PGD)", float(np.max(np.abs(dis - d_i)) / pm["PGD"]),
                0.0, atol=1e-6)
        r.check(f"{name}: final velocity / PGV", abs(v_i[-1]) / pm["PGV"], 0.0, atol=1e-6)
        r.check(f"{name}: final displacement / PGD", abs(d_i[-1]) / pm["PGD"], 0.0, atol=1e-6)
        # SRP 3.7.1 "no baseline drift" (D-EQK-03): the displacement must oscillate at the periods of
        # the target, not form an arc spanning the record (which puts most of its Fourier energy at
        # periods longer than the record); and the PGD cannot exceed the long-period spectral
        # displacement of the target itself (SD -> PGD as the oscillator period grows)
        dm = EL.drift_measures(a, dt, g_len)
        r.require(f"{name}: no baseline drift: displacement energy at periods > {dm['T']:g} s "
                  f"({100 * dm['frac_long']:.1f} %) < 50 %", dm["frac_long"] < 0.5)
        sd_low = float(tg(np.array([0.1]))[0]) * g_len / (2 * np.pi * 0.1) ** 2
        r.require(f"{name}: PGD {pm['PGD']:.1f} cm <= target SD at 0.1 Hz ({sd_low:.1f} cm)", pm["PGD"] <= sd_low)
        if psd:
            pc = EL.psd_check(a * g_len, dt, lambda f: EL.rg160_target_psd(f, 1.0, "SI"))
            r.require(f"{name}: PSD >= 80 % of SRP App. A target over 0.3-24 Hz (min ratio {pc['min_ratio']:.3f} "
                      f"at {pc['f_min']:.2f} Hz)", pc["passed"])
        r.notes.append(f"{name}: PGA {pm['PGA']:.3f} g, PGV {pm['PGV']:.1f} cm/s, PGD {pm['PGD']:.1f} cm, "
                       f"V/A {pm['V/A']:.3f} s, AD/V^2 {pm['AD/V2']:.2f}, D5-75 {pm['D5-75']:.2f} s, "
                       f"long-period displacement energy {100 * dm['frac_long']:.1f} % "
                       f"(RG 1.60 anchors: PGV 48 in/s = 121.9 cm/s, PGD 36 in = 91.4 cm per g, V/A 0.124 s, "
                       f"AD/V^2 6.0)")
    for x, y in (("H1", "H2"), ("H1", "V"), ("H2", "V")):
        rho = EL.correlation(accs[x], accs[y])
        r.require(f"|rho({x}, {y})| <= 0.16 ({rho:+.4f})", abs(rho) <= 0.16)
    r.require("listing reports the SRP 3.7.1 summary", "Acceptance criteria summary" in listing)
    r.notes.append("Vertical component: SRP 3.7.1 App. A defines the minimum PSD for the RG 1.60 horizontal "
                   "spectrum only, so no PSD criterion is checked for V.")
    r.notes.append("The SRP Approach 2 checks are made at 100 points per decade from 0.1 to 50 Hz; the manual "
                   "variant over the target range 0.1-100 Hz (Nyquist).")
    r.notes.append("Baseline drift: before the review fix the H components had PGD of about 410 cm with 91-95 % of "
                   "the displacement energy at periods longer than the record; drift is now controlled inside the "
                   "matching loop (high-pass + constrained least-squares displacement fit, projected wavelets).")
    return r
