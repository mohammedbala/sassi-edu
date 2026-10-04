"""Verification problems of ANALYS (impedance and flexible-volume SSI solution) and of the
SITE -> POINT -> HOUSE -> (FORCE) -> ANALYS -> MOTION chain (requirements 6.3, R2).

* VP-01   fixed-base hysteretic SDOF transfer-function peak (SPRING + mass on a very stiff site,
          SASSI and CMODFORM 1 damping forms), MOTION .TFU of the ANALYS FILE8
* VP-15   surface rigid massless mat under vertical SV and SH: H_u = 1, H_theta = 0 (simultaneous
          X/Y/Z cases, rigid BEAMS and rigid SHELL mats), MOTION .TFU/.TFI of FILE8X
* VP-16   zero-SSI identity: embedded FV model whose structure equals the excavated soil gives
          U = U'_f at every node for SV, SH and P input
* VP-22   restart equivalence: Mode 1 (save) / New Structure / New Seismic Environment / New
          Dynamic Loading / frequency subsets equal initiation runs
* VP-23   simultaneous cases: FILE8X/Y/Z = single runs with FILE1X/Y/Z; 3 vibration load cases
          FILE8001..003 = single runs
* VP-24   coordinate transformation angle on a rotationally symmetric model
* VP-39b  spring-mass ATF peak 10.0125 through ANALYS; GENERAL (global and local) identical
* VP-40   moving-load phase: FORCE factor 0.5 and arrival 0.1 s -> 0.5 x the reference response
          shifted by 0.1 s (ANALYS type 1 and MOTION)
* VP-41   low-frequency ATF ~ 1 at every node (G-19), MOTION on the ANALYS FILE8

The models are built with :mod:`sassi.verify.builders` as module decks and run exactly as the
RUN<MODULE> commands do, so every problem verifies the deck -> module -> file contracts as well
as the numbers.  Units: m, t, kN, s (g = 9.81 m/s2).
"""
from __future__ import annotations

import math
import shutil
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

from ...conventions import nodal_result_name
from ...core.flexibility import frequency_row
from ...core.freefield import free_field_at_nodes
from ...io.container import read_container
from .. import VPResult, problem, worse
from .. import builders as B

#: layered site of the SSI problems: (thick, Vs, Vp, rho, beta); Vp = 2 Vs (nu = 1/3)
SITE_LAYERS = [(2.0, 150.0, 300.0, 1.9, 0.05), (3.0, 250.0, 500.0, 2.0, 0.04), (5.0, 350.0, 700.0, 2.1, 0.03)]
SITE_HS = (600.0, 1200.0, 2.2, 0.02)
#: concrete (kPa, t/m3)
E_CONC, RHO_CONC = 3.0e7, 2.4


def site_a() -> B.Site:
    return B.layered_site(SITE_LAYERS, SITE_HS)


def stiff_site() -> B.Site:
    """A practically rigid site (Vs = 100 km/s): the soil impedance at a single interaction node is
    ~1e11 kN/m, so a 1-t oscillator on it is a fixed-base oscillator driven by the control motion."""
    return B.layered_site([(10.0, 1.0e5, 2.0e5, 2.0, 0.01)], (1.0e5, 2.0e5, 2.0, 0.01))


def _point(wd: Path, fs: B.FrequencySet, mdl: B.Model, model: str = "m") -> None:
    B.write_deck(wd, model, B.point_deck(fs, mdl.layer, mdl.rad, model=model))
    B.run("POINT", wd, model)


def _fresh(path: Path) -> Path:
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    return path


def _relmax(a, b) -> float:
    a, b = np.asarray(a), np.asarray(b)
    den = float(np.abs(b).max())
    return float(np.abs(a - b).max()) / (den if den > 0 else 1.0)


def read_tf_text(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """.TFU/.TFI written by MOTION: rows ``f |H| phase`` (D-FIL-02) -> (f, complex H)."""
    rows = [ln.split() for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    a = np.array([[float(v) for v in r] for r in rows])
    return a[:, 0], a[:, 1] * np.exp(1j * a[:, 2])


def read_history_text(path: Path) -> Tuple[float, np.ndarray]:
    """.ACC written by MOTION: first value dt, then one value per line."""
    vals = [float(ln) for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    return vals[0], np.asarray(vals[1:])


# ======================================================================================
# VP-01  fixed-base hysteretic SDOF
# ======================================================================================
VP01_BETAS = (0.02, 0.05, 0.10, 0.15, 0.20)
VP01_PEAK_SASSI = (25.0050, 10.0125, 5.0252, 3.3715, 2.5516)
VP01_RATIO = (0.99960, 0.99750, 0.98995, 0.97724, 0.95917)
VP01_PEAK_FORM1 = (25.0200, 10.0499, 5.0990, 3.4801, 2.6926)
VP01_DAMP = {0.05: 4.99, 0.10: 9.95, 0.15: 14.83, 0.20: 19.60}          # recovered 1/(2|H|max), %


@problem("VP-01", "Fixed-base hysteretic SDOF transfer-function peak (SPRING + mass through ANALYS)", tier="P0",
         modules=["SITE", "POINT", "HOUSE", "ANALYS", "MOTION"], source="R2 section 0 (Ostadan et al. 2004)")
def vp01(workdir: Path) -> VPResult:
    """SPRING (damping beta) + 1 t mass on one surface interaction node of a practically rigid site,
    vertical SV with the control point at the surface: the base moves with the control motion and
    the mass TF is ``H = k*/(k* - m w^2)``.  The stiffness is chosen so that the theoretical peak
    frequency (``w0 sqrt(1 - 2 beta^2)`` for the SASSI form, ``w0`` for CMODFORM 1) is the grid
    frequency f_p = 5 Hz of a fine grid (df = 1 mHz, +-5 points): the sample there must be the
    largest and equal the closed-form peak (frequency step permitting: the peak location is then
    resolved to df/f0 = 2e-4)."""
    r = VPResult()
    df, fp, m = 1.0e-3, 5.0, 1.0
    npk = int(round(fp / df))
    fs = B.FrequencySet.harmonic(df, range(npk - 5, npk + 6))
    worst_peak, worst_base = 0.0, 0.0
    for form in (0, 1):
        site = stiff_site()
        site.cmodform = form
        wd = _fresh(Path(workdir) / f"form{form}")
        B.run_soil(wd, "m", site, fs, layer=0, rad=1.0)
        for ib, beta in enumerate(VP01_BETAS):
            wp = 2.0 * math.pi * fp
            k = m * wp ** 2 / (1.0 - 2.0 * beta ** 2) if form == 0 else m * wp ** 2
            f0 = math.sqrt(k / m) / (2.0 * math.pi)
            mdl = B.sdof_on_node(site, k, m, beta)
            B.run_house(wd, "m", mdl)
            B.run_analys(wd, "m", fs)
            f8 = B.read_file8(wd)
            H = B.tf(f8, mdl["mass"], 1)
            Hb = B.tf(f8, mdl["base"], 1)
            imax = int(np.argmax(np.abs(H)))
            peak = float(np.abs(H[imax]))
            ratio = fs.fnum[imax] * df / f0
            if form == 0:
                exact = 1.0 / (2.0 * beta * math.sqrt(1.0 - beta ** 2))
                r.check(f"SASSI form beta = {beta:.2f}: |H|max vs 1/(2 b sqrt(1-b^2))", peak, exact, rtol=1e-6)
                r.check(f"SASSI form beta = {beta:.2f}: |H|max vs R2 table (4 decimals)", peak, VP01_PEAK_SASSI[ib],
                        atol=5e-5)
                r.check(f"SASSI form beta = {beta:.2f}: peak w/w0 vs sqrt(1-2b^2) (grid step df/f0)", ratio,
                        VP01_RATIO[ib], atol=df / f0 + 5e-6)
                if beta in VP01_DAMP:
                    r.check(f"SASSI form beta = {beta:.2f}: recovered damping 1/(2|H|max) in %", 100.0 / (2.0 * peak),
                            VP01_DAMP[beta], atol=0.005)
            else:
                exact = math.sqrt(1.0 + 4.0 * beta ** 2) / (2.0 * beta)
                r.check(f"CMODFORM 1 beta = {beta:.2f}: |H|max vs sqrt(1+4b^2)/(2b)", peak, exact, rtol=1e-6)
                r.check(f"CMODFORM 1 beta = {beta:.2f}: |H|max vs R2 table (4 decimals)", peak, VP01_PEAK_FORM1[ib],
                        atol=5e-5)
                r.check(f"CMODFORM 1 beta = {beta:.2f}: peak at w0 (grid step df/f0)", ratio, 1.0, atol=df / f0 + 5e-6)
            worst_peak = worse(worst_peak, abs(peak / exact - 1.0))
            worst_base = worse(worst_base, float(np.abs(Hb - 1.0).max()))
        if form == 0:
            # MOTION on the ANALYS FILE8 of the last SASSI-form run (beta = 0.20): TFU = FILE8 at the SSI
            # frequencies (the ANALYS -> MOTION contract)
            node = mdl["mass"]
            B.run_motion(wd, "m", fs, "", [(node, 1, 1, 0, 0, 0, 0, 0)], out=1, type=0, cm=0, ang=0.0)
            fu, Hu = read_tf_text(wd / nodal_result_name(node, 1, "TFU"))
            r.check("MOTION .TFU of the mass node vs ANALYS FILE8 (max rel. difference)", _relmax(Hu, H), 0.0,
                    atol=1e-12)
            r.check("MOTION .TFU frequencies vs FILE8 frequencies (max abs. difference, Hz)",
                    float(np.abs(fu - fs.freq).max()), 0.0, atol=1e-9)
    r.require("base node moves with the control motion: max |H_base - 1| < 1e-6", worst_base < 1e-6,
              f"observed {worst_base:.2e}")
    r.notes.append(f"Largest relative peak error {worst_peak:.2e} (soil flexibility of the stiff site); base motion "
                   f"deviates from the control motion by at most {worst_base:.2e}.  Grid df = {df:g} Hz around "
                   f"f_p = {fp:g} Hz: the peak location is resolved to df/f0 ~ 2e-4, enough to separate the SASSI "
                   "peak w0 sqrt(1-2b^2) from w0 (relative offset ~ b^2 >= 4e-4 for b >= 0.02).")
    return r


# ======================================================================================
# VP-15  surface rigid massless foundation
# ======================================================================================
@problem("VP-15", "Surface rigid massless foundation under vertical coherent SV/SH: H_u = 1, H_theta = 0",
         tier="P0", modules=["SITE", "POINT", "HOUSE", "ANALYS", "MOTION"], source="R2 E6")
def vp15(workdir: Path) -> VPResult:
    """Square 8 m x 8 m mat (4 x 4 cells) on the layered site, made rigid by stiff massless BEAMS to
    its centre node or by stiff massless SHELLs; SV (x'), SH (y') and P (z') fields run as the
    simultaneous X/Y/Z cases.  A uniform surface free field is a rigid-body motion of a massless
    foundation, so the foundation input motion equals the free field: H_u = 1, H_theta = 0."""
    r = VPResult()
    site = site_a()
    fs = B.FrequencySet.fourier(0.01, 1024, [5, 10, 20, 40, 80, 120, 160, 200])
    L = 4.0
    worst = 0.0
    for variant in ("beams", "shell"):
        wd = _fresh(Path(workdir) / variant)
        mdl = B.surface_rigid_mat(site, half_width=L, ndiv=4, rigid=variant)
        B.run_site_xyz(wd, "m", site, fs)
        _point(wd, fs, mdl)
        B.run_house(wd, "m", mdl)
        B.run_analys(wd, "m", fs, simul=1)
        c = mdl["centre"]
        for case, u, th, name in (("X", 1, 5, "SV"), ("Y", 2, 4, "SH")):
            f8 = B.read_file8(wd, f"FILE8{case}")
            Hu, Ht = B.tf(f8, c, u), B.tf(f8, c, th)
            e_u = float(np.abs(Hu - 1.0).max())
            e_t = float(np.abs(Ht).max()) * L
            worst = worse(worst, e_u, e_t)
            r.check(f"{variant} mat, vertical {name} ({case}): max |H_u - 1| over {len(fs.fnum)} frequencies",
                    e_u, 0.0, atol=0.01)
            r.check(f"{variant} mat, vertical {name} ({case}): max |H_theta| L", e_t, 0.0, atol=0.01)
            Hall = np.array([B.tf(f8, n, u) for n in mdl["mat"]])
            r.check(f"{variant} mat, {name}: every mat node, max |H_{'xy'[u - 1]} - 1|",
                    float(np.abs(Hall - 1.0).max()), 0.0, atol=0.01)
        if variant == "beams":
            # ANALYS -> MOTION: TFU = FILE8X at the SSI frequencies, |TFI| = 1 below the last SSI frequency
            B.run_motion(wd, "m", fs, "", [(c, 1, 1, 0, 0, 0, 0, 0)], out=1, file8="FILE8X", cm=0, ang=0.0)
            fu, Hu_m = read_tf_text(wd / nodal_result_name(c, 1, "TFU"))
            r.check("MOTION .TFU (centre UX, FILE8X) vs ANALYS (max rel. difference)", _relmax(Hu_m, B.tf(
                B.read_file8(wd, "FILE8X"), c, 1)), 0.0, atol=1e-12)
            fi, Hi = read_tf_text(wd / nodal_result_name(c, 1, "TFI"))
            band = (fi > 0) & (fi <= fs.freq[-1])
            r.check("MOTION .TFI (centre UX): max ||H| - 1| up to the last SSI frequency",
                    float(np.abs(np.abs(Hi[band]) - 1.0).max()), 0.0, atol=0.01)
    r.notes.append(f"Largest deviation {worst:.2e} (round-off: the uniform free field is an exact rigid-body "
                   "solution of the discrete system of a massless foundation).")
    return r


# ======================================================================================
# VP-16  zero-SSI identity
# ======================================================================================
@problem("VP-16", "Zero-SSI identity: embedded FV model whose structure equals the excavated soil", tier="P0",
         modules=["SITE", "POINT", "HOUSE", "ANALYS"], source="02 Part 6, R1 4.4")
def vp16(workdir: Path) -> VPResult:
    """6 m x 6 m box through the two upper layers (FV, every excavation node interacts); the
    structure is made of SOLIDs with M-table type-3 materials equal to the layers (K*_s = K*_e,
    M_s = M_e).  The flexible-volume equation then reduces to X_ff (U - U'_f) = 0: U = U'_f at every
    node, for SV, SH and P input (simultaneous X/Y/Z cases)."""
    r = VPResult()
    site = site_a()
    fs = B.FrequencySet.fourier(0.01, 1024, [4, 10, 20, 41, 82, 123, 164, 205])
    wd = _fresh(Path(workdir) / "box")
    mdl = B.embedded_box(site, half_width=3.0, n_emb=2, ndiv=4, method="FV", structure="soil")
    B.run_site_xyz(wd, "m", site, fs)
    _point(wd, fs, mdl)
    f4 = B.run_house(wd, "m", mdl)
    B.run_analys(wd, "m", fs, simul=1)
    nodes = np.asarray(f4["int_node"])
    xyz = np.asarray(f4["x_int_xyz"])
    iface = np.asarray(f4["int_iface"])
    r.require("every excavation node is an interaction node (FV)", len(nodes) == len(mdl["excavation"]))
    for case, comp, name in (("X", 0, "SV"), ("Y", 1, "SH"), ("Z", 2, "P")):
        f8 = B.read_file8(wd, f"FILE8{case}")
        f1 = read_container(wd / f"FILE1{case}", "FILE1")
        err, scale = 0.0, 0.0
        for q, n in enumerate(fs.fnum):
            Up = free_field_at_nodes(f1, frequency_row(f1, n), xyz, iface, 0.0, 0.0, 0.0)
            H = np.array([[B.tf(f8, int(nd), d)[q] for d in (1, 2, 3)] for nd in nodes])
            err = worse(err, float(np.abs(H - Up).max()))
            scale = worse(scale, float(np.abs(Up).max()))
        r.check(f"vertical {name} ({case}): max |U - U'_f| / max |U'_f| over all nodes and frequencies",
                err / scale, 0.0, atol=1e-8)
        Hc = B.tf(f8, mdl["top_centre"], comp + 1)
        r.check(f"vertical {name} ({case}): control point (surface centre) max |ATF - 1|",
                float(np.abs(Hc - 1.0).max()), 0.0, atol=1e-8)
    return r


# ======================================================================================
# VP-22  restart equivalence
# ======================================================================================
def _vp22_model(site: B.Site, E_stick: float) -> B.Model:
    stick = B.Stick(heights=[4.0, 8.0], masses=[200.0, 150.0], E=E_stick, A=4.0, I=3.0, beta=0.05)
    return B.embedded_box(site, half_width=3.0, n_emb=2, ndiv=4, method="FSIN", structure="shell", E=E_CONC,
                          rho=RHO_CONC, beta=0.05, thick=0.4, stick=stick)


@problem("VP-22", "Restart equivalence: ANALYS modes 1/2/3 with COOX/COOTK equal initiation runs", tier="P0",
         modules=["ANALYS"], source="05c C3")
def vp22(workdir: Path) -> VPResult:
    """Embedded FI-FSIN box with a shell basement and a stick (Schur path with excavated
    non-interaction nodes).  Initiation with Save Restart Files, then: New Structure with the same
    FILE4; New Seismic Environment with the same FILE1, also on a frequency subset; New Structure
    with a changed stick vs a fresh initiation of the changed model; New Seismic Environment with
    an SH field vs an initiation; New Dynamic Loading with a FORCE load vs an initiation.  Restart
    files of another structure are refused (FILE90 hashes)."""
    r = VPResult()
    site = site_a()
    fs = B.FrequencySet.fourier(0.01, 1024, [3, 8, 20, 41, 82, 123, 164, 205])
    wd = _fresh(Path(workdir) / "restart")
    mA = _vp22_model(site, E_CONC)
    B.run_soil(wd, "m", site, fs, layer=mA.layer, rad=mA.rad)
    B.run_house(wd, "m", mA)
    B.run_analys(wd, "m", fs, mode=0, save=1)
    ref = B.read_file8(wd)["H"].copy()
    r.require("COOX001..COOX008, COOTK001..COOTK008, COOXI and COOTKI written",
              all((wd / n).exists() for n in ["COOXI", "COOTKI"] + [f"COOX{q:03d}" for q in range(1, 9)]
                  + [f"COOTK{q:03d}" for q in range(1, 9)]))
    B.run_analys(wd, "m", fs, mode=1)
    r.check("Mode 1 (New Structure, same FILE4) vs initiation: max rel. difference of FILE8",
            _relmax(B.read_file8(wd)["H"], ref), 0.0, atol=1e-10)
    B.run_analys(wd, "m", fs, mode=2)
    r.check("Mode 2 (New Seismic Environment, same FILE1) vs initiation", _relmax(B.read_file8(wd)["H"], ref), 0.0,
            atol=1e-10)
    sub = [8, 82, 164]
    B.run_analys(wd, "m", fs.subset(sub), mode=2)
    rows = [fs.fnum.index(n) for n in sub]
    f8 = B.read_file8(wd)
    r.require("frequency subset: FILE8 holds the subset frequency numbers", list(f8["fnum"]) == sub)
    r.check("Mode 2 on a frequency subset (3 of 8) vs initiation rows", _relmax(f8["H"], ref[rows]), 0.0, atol=1e-10)
    # ---- New Structure: stiffer stick, then a fresh initiation of the changed model ------------------
    mB = _vp22_model(site, 3.0 * E_CONC)
    B.run_house(wd, "m", mB)
    rc = B.run_analys(wd, "m", fs, mode=2, check=False)
    msg = B.listing(wd, "m", "ANALYS")
    r.require("Mode 2 after a structure change is refused (FILE90 hash of FILE4)", rc != 0 and "hash" in msg)
    B.run_analys(wd, "m", fs, mode=1, save=1)
    H1 = B.read_file8(wd)["H"].copy()
    B.run_analys(wd, "m", fs, mode=0, save=0)
    H0 = B.read_file8(wd)["H"].copy()
    r.check("Mode 1 (New Structure: stick 3x stiffer) vs initiation of the new model", _relmax(H1, H0), 0.0, atol=1e-10)
    r.require("the structure change changes the TFs (> 1e-3)", _relmax(H1, ref) > 1e-3)
    # ---- New Seismic Environment: SH field (factorisation of the new structure reused) ---------------
    B.write_deck(wd, "m", B.site_deck(site, fs, wave="SH", model="m"))
    B.run("SITE", wd, "m")
    B.run_analys(wd, "m", fs, mode=2)
    Hs = B.read_file8(wd)["H"].copy()
    B.run_analys(wd, "m", fs, mode=0)
    r.check("Mode 2 (New Seismic Environment: SH instead of SV) vs initiation", _relmax(Hs, B.read_file8(wd)["H"]),
            0.0, atol=1e-10)
    # ---- New Dynamic Loading: FORCE loads at the stick top and the box ------------------------------------
    top = mB["top"]
    B.run_force(wd, "m", [(top, 1, 1.0, 0.0), (top, 3, 0.5, 0.05), (mB["stick"][0], 5, 2.0, 0.02)], fs)
    B.run_analys(wd, "m", fs, type=1, mode=3)
    Hv = B.read_file8(wd)["H"].copy()
    B.run_analys(wd, "m", fs, type=1, mode=0)
    r.check("Mode 3 (New Dynamic Loading) vs vibration initiation", _relmax(Hv, B.read_file8(wd)["H"]), 0.0,
            atol=1e-10)
    return r


# ======================================================================================
# VP-23  simultaneous cases
# ======================================================================================
def _stick_mat(site: B.Site) -> B.Model:
    stick = B.Stick(heights=[3.0, 6.0, 9.0], masses=[100.0, 100.0, 80.0], E=E_CONC, A=3.0, I=2.0, beta=0.05)
    return B.stick_on_mat(site, stick, half_width=4.0, ndiv=4, mat_mass=150.0)


@problem("VP-23", "Simultaneous cases: FILE1X/Y/Z and vibration load cases equal single runs", tier="P0",
         modules=["ANALYS"], source="05c C4")
def vp23(workdir: Path) -> VPResult:
    r = VPResult()
    site = site_a()
    fs = B.FrequencySet.fourier(0.01, 1024, [4, 12, 30, 61, 102, 153, 204])
    wd = _fresh(Path(workdir) / "simul")
    mdl = _stick_mat(site)
    B.run_site_xyz(wd, "m", site, fs)
    _point(wd, fs, mdl)
    B.run_house(wd, "m", mdl)
    B.run_analys(wd, "m", fs, simul=1)
    sim = {c: B.read_file8(wd, f"FILE8{c}") for c in "XYZ"}
    for k, c in enumerate("XYZ"):
        r.require(f"FILE8{c} meta: case {c}, control direction {k}", sim[c].meta["case"] == c and sim[c].meta["cm"] == k)
        B.copy_file(wd, f"FILE1{c}", "FILE1")
        B.run_analys(wd, "m", fs, simul=0)
        r.check(f"FILE8{c} (simul = 1) vs single run with FILE1 = FILE1{c}", _relmax(sim[c]["H"], B.read_file8(wd)["H"]),
                0.0, atol=1e-12)
    top, s1 = mdl["top"], mdl["stick"][0]
    cases = [[(top, 1, 1.0, 0.0)], [(top, 2, 0.7, 0.05)], [(s1, 3, 1.0, 0.0), (top, 4, 0.5, 0.1)]]
    for k, loads in enumerate(cases, start=1):
        B.run_force(wd, "m", loads, fs, copy_to=f"FILE9{k:03d}")
    B.run_analys(wd, "m", fs, type=1, simul=3)
    vib = {k: B.read_file8(wd, f"FILE8{k:03d}") for k in (1, 2, 3)}
    for k in (1, 2, 3):
        B.copy_file(wd, f"FILE9{k:03d}", "FILE9")
        B.run_analys(wd, "m", fs, type=1, simul=0)
        r.check(f"FILE8{k:03d} (3 simultaneous load cases) vs single run with FILE9 = FILE9{k:03d}",
                _relmax(vib[k]["H"], B.read_file8(wd)["H"]), 0.0, atol=1e-12)
        r.require(f"FILE8{k:03d} meta case = {k}", int(vib[k].meta["case"]) == k)
    return r


# ======================================================================================
# VP-24  coordinate transformation angle
# ======================================================================================
@problem("VP-24", "Coordinate transformation angle on a rotationally symmetric model", tier="P0",
         modules=["ANALYS"], source="05c C5")
def vp24(workdir: Path) -> VPResult:
    """Stick on a square rigid mat (symmetric under 90-degree rotations and reflections, isotropic
    stick section): the x' (SV) input at angle a gives, on the symmetry axis (mat centre and stick
    nodes), ATF_x(a) = cos a ATF_x(0), ATF_y(a) = sin a ATF_x(0) and the rotation vector rotated by a;
    at every node H(a) = cos a H(0) + sin a H(90 deg) (superposition of the rotated free field)."""
    r = VPResult()
    site = site_a()
    fs = B.FrequencySet.fourier(0.01, 1024, [4, 12, 30, 61, 102, 153, 204])
    wd = _fresh(Path(workdir) / "angle")
    mdl = _stick_mat(site)
    B.run_soil(wd, "m", site, fs, layer=0, rad=mdl.rad)
    B.run_house(wd, "m", mdl)
    H: Dict[float, object] = {}
    for a in (0.0, 90.0, 30.0, 135.0):
        B.run_analys(wd, "m", fs, ang=a)
        H[a] = B.read_file8(wd)
        r.require(f"FILE8 meta ang = {a:g}", float(H[a].meta["ang"]) == a)
    axis = [mdl["centre"]] + list(mdl["stick"])
    for a in (30.0, 135.0):
        ca, sa = math.cos(math.radians(a)), math.sin(math.radians(a))
        e_x = e_y = e_r = 0.0
        scale = max(float(np.abs(B.tf(H[0.0], n, 1)).max()) for n in axis)
        rscale = max(float(np.abs(B.tf(H[0.0], n, 5)).max()) for n in axis)
        for n in axis:
            hx0, ry0 = B.tf(H[0.0], n, 1), B.tf(H[0.0], n, 5)
            e_x = worse(e_x, float(np.abs(B.tf(H[a], n, 1) - ca * hx0).max()) / scale)
            e_y = worse(e_y, float(np.abs(B.tf(H[a], n, 2) - sa * hx0).max()) / scale)
            e_r = worse(e_r, float(np.abs(B.tf(H[a], n, 5) - ca * ry0).max()) / rscale,
                      float(np.abs(B.tf(H[a], n, 4) + sa * ry0).max()) / rscale)
        r.check(f"a = {a:g} deg, axis nodes: max |ATF_x(a) - cos a ATF_x(0)| / max|ATF_x(0)|", e_x, 0.0, atol=1e-8)
        r.check(f"a = {a:g} deg, axis nodes: max |ATF_y(a) - sin a ATF_x(0)| / max|ATF_x(0)|", e_y, 0.0, atol=1e-8)
        r.check(f"a = {a:g} deg, axis nodes: rotations (ROTX, ROTY) = R(a) (0, ROTY(0))", e_r, 0.0, atol=1e-8)
        comb = ca * np.asarray(H[0.0]["H"]) + sa * np.asarray(H[90.0]["H"])
        r.check(f"a = {a:g} deg, all nodes and DOFs: H(a) vs cos a H(0) + sin a H(90)", _relmax(H[a]["H"], comb), 0.0,
                atol=1e-8)
    return r


# ======================================================================================
# VP-39b  spring-mass ATF through ANALYS, GENERAL equivalence
# ======================================================================================
@problem("VP-39b", "Spring-mass SDOF ATF peak through ANALYS; GENERAL (global and local axes) identical",
         tier="P0", modules=["HOUSE", "ANALYS"], source="requirements 6.3 VP-39, spec 08 sec. 11 (C-28)")
def vp39b(workdir: Path) -> VPResult:
    """SC,1,1000,0,0,0,0,0,0.05 with MT = 1 on a practically rigid site (fixed base), and the same
    oscillator as a 2-node GENERAL element (MXR = k(1-2b^2), MXI = 2kb sqrt(1-b^2)) and a 3-node GENERAL
    element in local axes.  ATF peak on a grid *aligned* with the peak frequency w0 sqrt(1 - 2b^2)
    (df = f_peak / 25000 ~ 0.2 mHz, frequency numbers 24996..25004), so the closed-form peak is a grid
    sample and the requirements 6.3 tolerance 1e-8 applies to the ANALYS chain itself."""
    r = VPResult()
    k, m, beta = 1000.0, 1.0, 0.05
    f0 = math.sqrt(k / m) / (2.0 * math.pi)
    fpk = f0 * math.sqrt(1.0 - 2.0 * beta ** 2)
    npk = 25000
    df = fpk / npk
    fs = B.FrequencySet.harmonic(df, range(npk - 4, npk + 5))
    site = stiff_site()
    wd = _fresh(Path(workdir) / "sdof")
    B.run_soil(wd, "m", site, fs, layer=0, rad=1.0)
    out = {}
    for element in ("spring", "general", "general3"):
        mdl = B.sdof_on_node(site, k, m, beta, element=element)
        B.run_house(wd, "m", mdl)
        B.run_analys(wd, "m", fs)
        out[element] = B.tf(B.read_file8(wd), mdl["mass"], 1)
    H = out["spring"]
    peak = float(np.abs(H).max())
    exact = 1.0 / (2.0 * beta * math.sqrt(1.0 - beta ** 2))
    r.check("SPRING-mass ATF peak vs 1/(2 b sqrt(1-b^2)) (aligned grid, df = f_peak/25000)", peak, exact,
            rtol=1e-8)
    r.check("SPRING-mass ATF peak vs 10.0125 (4 decimals)", peak, 10.0125, atol=5e-5)
    r.check("undamped frequency sqrt(k/m)/(2 pi) vs 5.033 Hz (3 decimals)", f0, 5.033, atol=5e-4)
    imax = int(np.argmax(np.abs(H)))
    r.require("the largest grid sample is the one at f0 sqrt(1-2b^2) (frequency number 25000)",
              fs.fnum[imax] == npk, f"largest at number {fs.fnum[imax]}")
    for element in ("general", "general3"):
        r.check(f"{element} ATF vs SPRING ATF (max rel. difference)", _relmax(out[element], H), 0.0, atol=1e-8)
    r.notes.append(f"Peak sample |H| = {peak:.12f} at f = {fs.fnum[imax] * df:.9f} Hz; relative error "
                   f"{abs(peak / exact - 1.0):.2e} against 1/(2 b sqrt(1-b^2)) = {exact:.12f} (tolerance 1e-8). "
                   "The grid is aligned with the peak: on a non-aligned grid the sampling offset delta costs "
                   "(delta/f0)^2/(2 b^2) relative.")
    return r


# ======================================================================================
# VP-40  moving load phase
# ======================================================================================
def _ref_load(n: int, dt: float) -> np.ndarray:
    """Zero-mean reference load history: one cycle of a 2 Hz sine after 0.5 s and a smaller 5 Hz
    double cycle, minus the record mean (zero DC term, so the f = 0 value of H does not matter)."""
    t = np.arange(n) * dt
    f = np.where((t >= 0.5) & (t < 1.0), np.sin(2 * np.pi * 2.0 * (t - 0.5)), 0.0)
    f += np.where((t >= 1.2) & (t < 1.6), 0.4 * np.sin(2 * np.pi * 5.0 * (t - 1.2)), 0.0)
    return f - f.mean()


@problem("VP-40", "Moving load phase: FORCE factor 0.5, arrival 0.1 s through ANALYS and MOTION", tier="P0",
         modules=["FORCE", "ANALYS", "MOTION"], source="05b section 11, D-FRC-01")
def vp40(workdir: Path) -> VPResult:
    """Horizontal load at the top of a stick on a rigid mat.  Reference: factor 1, arrival 0; test:
    factor 0.5, arrival 0.1 s.  The SSI frequencies are every Fourier frequency (no interpolation),
    the reference load has zero mean; then the displacement response of the test is exactly
    0.5 x the reference response circularly shifted by 0.1 s (5 samples)."""
    r = VPResult()
    site = site_a()
    nft, dt = 256, 0.02
    fs = B.FrequencySet.fourier(dt, nft, range(1, nft // 2 + 1))
    wd = _fresh(Path(workdir) / "moving")
    mdl = _stick_mat(site)
    B.run_soil(wd, "m", site, fs, layer=0, rad=mdl.rad, mode2=False)
    B.run_house(wd, "m", mdl)
    top = mdl["top"]
    B.write_history(wd / "load.th", _ref_load(nft, dt), dt)
    res = {}
    for tag, fac, t0 in (("ref", 1.0, 0.0), ("test", 0.5, 0.1)):
        B.run_force(wd, "m", [(top, 1, fac, t0)], fs)
        B.run_analys(wd, "m", fs, type=1)
        B.copy_file(wd, "FILE8", f"FILE8_{tag}")
        B.run_motion(wd, "m", fs, "load.th", [(top, 1, 1, 1, 0, 0, 0, 0)], type=1, resp=0, file8=f"FILE8_{tag}")
        _, u = read_history_text(wd / nodal_result_name(top, 1, "ACC"))
        res[tag] = (B.read_file8(wd, f"FILE8_{tag}"), u)
    f8r, ur = res["ref"]
    f8t, ut = res["test"]
    w = 2 * np.pi * fs.freq
    expect = 0.5 * np.exp(-1j * w * 0.1)[:, None] * np.asarray(f8r["H"])
    r.check("FILE8 (test) vs 0.5 exp(-i w 0.1) FILE8 (ref), all DOFs", _relmax(f8t["H"], expect), 0.0, atol=1e-10)
    shift = int(round(0.1 / dt))
    r.check("displacement history of the loaded node: u_test(t) vs 0.5 u_ref(t - 0.1 s) (circular)",
            _relmax(ut, 0.5 * np.roll(ur, shift)), 0.0, atol=1e-10)
    r.require("the reference response is not trivial", float(np.abs(ur).max()) > 0)
    return r


# ======================================================================================
# VP-41  low-frequency ATF
# ======================================================================================
@problem("VP-41", "Low-frequency ATF ~ 1 at all nodes (G-19) and MOTION on the ANALYS FILE8", tier="P0",
         modules=["ANALYS", "MOTION"], source="03 section 7 (G-19)")
def vp41(workdir: Path) -> VPResult:
    """Embedded FFV box (three layers, shell basement) with a stick; vertical SV, x input.  At the
    first SSI frequency (0.098 Hz) every node moves rigidly with the free field: |ATF_x| within 5 % of
    1 at every node; MOTION reads the FILE8 (TF output only) and reports no EDU-18 deviation."""
    r = VPResult()
    site = site_a()
    fs = B.FrequencySet.fourier(0.01, 1024, [1, 4, 16, 41, 82, 164])
    wd = _fresh(Path(workdir) / "lowfreq")
    stick = B.Stick(heights=[4.0, 8.0], masses=[300.0, 200.0], E=E_CONC, A=4.0, I=3.0, beta=0.05)
    mdl = B.embedded_box(site, half_width=3.0, n_emb=3, ndiv=4, method="FFV", structure="shell", E=E_CONC,
                         rho=RHO_CONC, beta=0.05, thick=0.4, stick=stick)
    B.run_soil(wd, "m", site, fs, layer=mdl.layer, rad=mdl.rad)
    B.run_house(wd, "m", mdl)
    B.run_analys(wd, "m", fs)
    f8 = B.read_file8(wd)
    H = np.asarray(f8["H"])
    ux = np.asarray(f8["eq_dof"]) == 1
    dev = float(np.abs(np.abs(H[0, ux]) - 1.0).max())
    r.check(f"max ||ATF_x| - 1| at f_1 = {fs.freq[0]:.4g} Hz over {int(ux.sum())} nodes", dev, 0.0, atol=0.05)
    r.require("ANALYS listing: no EDU-18 warning", "EDU-18" not in B.listing(wd, "m", "ANALYS"))
    nodes = [mdl["top"], mdl["top_centre"], mdl["bottom_centre"]]
    B.run_motion(wd, "m", fs, "", [(n, 1, 1, 0, 0, 0, 0, 0) for n in nodes], out=1, cm=0, ang=0.0)
    for n in nodes:
        _, Hm = read_tf_text(wd / nodal_result_name(n, 1, "TFU"))
        r.check(f"MOTION .TFU node {n} UX vs ANALYS FILE8", _relmax(Hm, B.tf(f8, n, 1)), 0.0, atol=1e-12)
    r.require("MOTION listing: no EDU-18 warning", "EDU-18" not in B.listing(wd, "m", "MOTION"))
    r.notes.append(f"Largest deviation of |ATF_x| from 1 at f_1: {dev:.2e}")
    return r
