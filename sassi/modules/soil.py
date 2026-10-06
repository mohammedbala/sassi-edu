"""SOIL module: SHAKE-type equivalent-linear free-field site response (SOIL-EQL).

Requirements section 4.11; decisions D-SOL-01 ... D-SOL-13, D-CNV-03, D-FIL-02, D-FIL-10;
theory R1 section 6.  The numerical core is :mod:`sassi.core.shake`.

Input
-----
``<model>.soi`` (written by AFWRITE, schema ``sassi.io.decks.SOIL``):

* ``profile`` table: SPRO sublayers top first (thickness, unit weight, Vp, Vs, damping and
  DYNP label resolved from the L table); the **last row is the half-space** (D-SOL-02).
* ``dynp`` table: curve points (strain %, G/Gmax) and (strain %, damping %) per label.
* control motion ``thfile`` (accelerations in g): ``header`` lines skipped, ``nrval``
  values read, scaled by ``mult`` or to peak ``max`` (exactly one non-zero, D-SOL-11).
* the motion is applied at the top of sublayer ``cl`` as an outcrop (``outcrop`` = 1,
  amplitude 2E) or within (E+F) motion.

Algorithm
---------
SHAKE recursion with complex moduli (``cmodform``), strains at mid-height of every
sublayer, effective strain ``ratio * max|gamma(t)|``, curves interpolated linearly in
log10(strain) with constant extrapolation, exactly ``iter`` property updates (D-SOL-04),
the final response computed with the last updated properties.  ``indir`` = 1 (vertical
input, SOILX) uses Vp, the P-wave damping and no iterations.  Fourier components above
``soilcutoff`` (``EDUOPT,SOILCUTOFF``) are removed; otherwise the cut-off is the Nyquist
frequency (``<cof>`` is disabled, D-SOL-08).

Outputs (model directory)
-------------------------
* listing ``<model>_SOIL.out``: input echo, iteration tables, strain-compatible properties,
  maximum accelerations / strains / stresses, amplification and spectra summaries;
* ``ACCxxx.TH``  acceleration (g) at the top of sublayer xxx (SACC opt 2), within or outcrop;
* ``SNxxx.TH`` / ``SSxxx.TH`` strain (%) / stress (model units) at mid-height (SSTR);
* ``RSxxx_zz.RS``  response spectrum at the top of sublayer xxx for damping number zz
  (SRS save = 1): rows ``f SA`` with SA in g times ``gravmult`` (D-SOL-07);
* ``SAFxxxO_yyyW_zz.RS``  spectral amplification factor RS(xxx)/RS(yyy) for damping number zz
  sampled at f = k * freqstep, k >= 1, up to the cut-off (SSAF, spec 05a 6.5 / 07 9.2.35;
  O = outcrop, W = within motion), and, as a documented extra output, ``SAFxxxO_yyyW.TFU``
  the Fourier amplification X(xxx)/X(yyy) at the same frequencies (SHAKE91 option 10);
* ``FILE73`` soil material curves and ``FILE88`` strain-compatible properties (``save`` = 1),
  formats in :func:`sassi.core.shake.write_file73` / :func:`~sassi.core.shake.write_file88`.

FILE88 holds, per soil sublayer (half-space excluded, matched to TOPL by position,
D-SOL-12): thickness, effective strain, G, Vs, beta_s, Vp, beta_p.  Each row is a
strain-compatible set: the effective strain is that of the last iteration, from which G and
beta_s were read (SHAKE91 Table B-2 pairing).  In iterated sublayers Vp follows Vs at
constant Poisson's ratio and beta_p = beta_s (D-SOL-06 default policy); linear sublayers
(no DYNP label) keep their input Vp and beta_p.

SOIL-NON (requirements 1.2 and 3.4.O, spec 05a section 6.6, manual 9.17.19-9.17.20)
---------------------------------------------------------------------------------
When the side file ``<deck stem>.nls`` (written by AFWRITE from the NLSOIL / NLSLAYER commands, see
:func:`sassi.core.soilnon.nls_text`) holds ``NLSOIL <Opt>`` = 1, SOIL runs the nonlinear time-domain
analysis of :mod:`sassi.core.soilnon` instead of SHAKE: MKZ hyperbolic backbone per sublayer
(NLSLAYER parameters or curve fit to the DYNP G/Gmax curve), extended Masing unloading/reloading,
viscous small-strain damping (NLDampType 1 frequency independent, 2 visco-elastic, 3 Rayleigh), rigid
(BedInt 0, within input) or elastic (BedInt 1, Joyner-Chen dashpot, outcrop input) base, Newmark
average acceleration with sub-increments and Newton-Raphson equilibrium iterations.  The input must
be the motion at bedrock (control layer = the half-space, the last SPRO sublayer; manual warning).
The output files have the same names as in SOIL-EQL (manual): ACCxxx.TH, SNxxx.TH, SSxxx.TH, RS
files, SSAF files, FILE73 and FILE88 (equivalent-linear properties of the nonlinear run, decision SN-12
of :mod:`sassi.core.soilnon`).  Without the side file, or with ``<Opt>`` = 0, SOIL is SOIL-EQL.
"""
from __future__ import annotations

import math
import re
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .. import conventions as C
from ..core import shake as SH
from ..core import soilnon as SN
from ..core import spectra as SP
from ..io import deckfmt, decks, textfiles, thfile
from ..io import library as LIB
from .base import ModuleContext, ModuleError, batch_main

NAME = "SOIL"

RS_FREQS = SP.log_frequencies(0.1, 100.0, 301)    # 100 points per decade (EDU-15: >= 301 points)
MAX_SAF_POINTS = 200001
MAX_SAF_RS_POINTS = 20000          # SSAF response-spectrum ratio frequencies (cost: one SDOF per point)

# ---------------------------------------------------------------------------------------
# Deck I/O.  Local workaround for a schema defect (reported to the lead): the shared SOIL
# schema types ssaf.freqstep as int, but SSAF <freqstep> is a frequency step in Hz (manual
# example 0.1, SHAKE91 0.125).  These two functions apply the schema of sassi.io.decks with
# that single column read/written as float; they stay valid once the schema is corrected.
# ---------------------------------------------------------------------------------------
_FLOAT_COLUMNS = {("ssaf", "freqstep")}


def _table_types(name: str):
    spec = decks.SCHEMAS[NAME].table_spec(name)
    return [(c, float if (name, c) in _FLOAT_COLUMNS else t) for c, t in spec.columns]


def _coerce(v, typ):
    if typ is str:
        if isinstance(v, str):
            t = v.strip()
            return deckfmt.parse_value(t) if t.startswith('"') else v
        return str(v)
    if isinstance(v, str):
        t = v.strip().replace("D", "E").replace("d", "e")
        f = float(t) if t else 0.0
    else:
        f = float(v)
    if typ is int:
        if f != int(round(f)):
            raise ValueError(f"expected an integer, got {v!r}")
        return int(round(f))
    return f


def read_deck(path) -> decks.Deck:
    """Read a SOIL deck (``decks.read`` semantics: missing parameters default with a warning)."""
    raw = deckfmt.read_raw(path)
    if raw.module != NAME:
        raise ValueError(f"{Path(path).name}: expected a {NAME} deck, found {raw.module}")
    s = decks.SCHEMAS[NAME]
    d = decks.new(NAME)
    for p in s.params:
        if p.name in raw.params:
            d.params[p.name] = _coerce(deckfmt.parse_value(raw.params[p.name]) if p.type is not str
                                       else raw.params[p.name], p.type)
        else:
            warnings.warn(f"{NAME} deck: parameter '{p.name}' missing, default {p.default!r} used")
    for t in s.tables:
        rt = raw.tables.get(t.name)
        if rt is None:
            continue
        cols = _table_types(t.name)
        if rt.columns != [c for c, _ in cols]:
            raise ValueError(f"{NAME} deck table {t.name}: columns {rt.columns} != {[c for c, _ in cols]}")
        for row in rt.rows:
            d.tables[t.name].rows.append([_coerce(v, typ) for v, (_, typ) in zip(row, cols)])
    return d


def write_deck(path, d: decks.Deck) -> Path:
    """Write a SOIL deck (``decks.write`` semantics with ``ssaf.freqstep`` as float)."""
    s = decks.SCHEMAS[NAME]
    unknown = set(d.params) - {p.name for p in s.params}
    if unknown:
        raise KeyError(f"{NAME} deck: unknown parameters {sorted(unknown)}")
    params = {p.name: _coerce(d.params.get(p.name, p.default), p.type) for p in s.params}
    tables = {}
    for t in s.tables:
        cols = _table_types(t.name)
        out = deckfmt.Table(columns=[c for c, _ in cols])
        tab = d.tables.get(t.name)
        if tab is not None:
            if list(tab.columns) != out.columns:
                raise ValueError(f"{NAME} deck table {t.name}: columns {tab.columns} != {out.columns}")
            for row in tab.rows:
                out.rows.append([_coerce(v, typ) for v, (_, typ) in zip(row, cols)])
        tables[t.name] = out
    return deckfmt.write_raw(path, NAME, params, tables, version=decks.DECK_VERSION)


# ---------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------
def soil_rs_name(layer: int, idamp: int) -> str:
    """Response-spectrum file of SOIL: ``RS001_01.RS`` (layer 1, damping number 1)."""
    return f"RS{layer:03d}_{idamp:02d}.RS"


def _saf_stem(layer: int, outcrop1: bool, layer2: int, outcrop2: bool) -> str:
    return f"SAF{layer:03d}{'O' if outcrop1 else 'W'}_{layer2:03d}{'O' if outcrop2 else 'W'}"


def saf_tf_name(layer: int, outcrop1: bool, layer2: int, outcrop2: bool) -> str:
    """SSAF Fourier amplification file, e.g. ``SAF001O_017O.TFU`` (O outcrop, W within)."""
    return _saf_stem(layer, outcrop1, layer2, outcrop2) + ".TFU"


def saf_rs_name(layer: int, outcrop1: bool, layer2: int, outcrop2: bool, idamp: int) -> str:
    """SSAF response-spectrum ratio file, e.g. ``SAF001O_017O_01.RS``."""
    return _saf_stem(layer, outcrop1, layer2, outcrop2) + f"_{idamp:02d}.RS"


def _units(grav: float) -> Dict[str, str]:
    if C.unit_system(grav) == "BS":
        return {"len": "ft", "wt": "kcf", "mod": "ksf", "vel": "ft/s", "acc": "ft/s^2"}
    return {"len": "m", "wt": "kN/m3", "mod": "kN/m2", "vel": "m/s", "acc": "m/s^2"}


def _resolve(ctx: ModuleContext, name: str) -> Path:
    """Input file of the deck: ``@name`` from the built-in library (D-W5-01), else absolute or in the model
    directory."""
    return LIB.module_path(name, ctx.workdir)


def _history_header(d: decks.Deck, L) -> int:
    """Header lines skipped before the values.  A built-in record (``@`` name) holds its time step on line 1:
    SOIL skips it when ``<header>`` is 0 and checks that it equals ``<delt>`` (D-W5-06)."""
    header = int(d["header"])
    name = str(d["thfile"] or "")
    dt_lib = LIB.history_dt(name)
    if dt_lib is None:
        return header
    L.write(f" {LIB.note(name)}")
    dt = float(d["delt"])
    if dt > 0 and abs(dt_lib - dt) > 1e-6 * dt:
        raise ModuleError(f"time step of the built-in record {name.strip()} ({dt_lib:g} s) differs from the SOIL "
                          f"time step {dt:g} s")
    if header == 0:
        L.write(f" {name.strip()}: line 1 is the time step; SOIL skips it (<header> 0 read as 1, D-W5-06)")
        return 1
    return header


def _validate(d: decks.Deck, nlay: int) -> None:
    """Fatal input errors (CHECK catalogue numbers of Manual Ch. 10, repeated defensively)."""
    if d["grav"] <= 0:
        raise ModuleError("Error 1: gravity acceleration must be > 0")
    if d["header"] < 0:
        raise ModuleError("Error 103: number of header lines < 0")
    if d["iter"] < 0:
        raise ModuleError("Error 105: number of iterations < 0")
    if not (0.0 < d["ratio"] < 1.0):
        raise ModuleError("Error 106: strain ratio not between 0 and 1")
    if d["gravmult"] <= 0:
        raise ModuleError("Error 108: gravity multiplier must be > 0")
    if d["cof"] < 0:
        raise ModuleError("Error 101: cut-off frequency < 0")
    if nlay < 2:
        raise ModuleError("SOIL profile needs at least one soil sublayer and the half-space (D-SOL-02)")
    if not (1 <= d["cl"] <= nlay):
        raise ModuleError(f"Error 104: illegal control layer {d['cl']} (1..{nlay})")
    if d["mult"] == 0 and d["max"] == 0:
        raise ModuleError("Error 77: multiplication factor and maximum value are both zero")
    if d["mult"] != 0 and d["max"] != 0:
        raise ModuleError("Error 78: give either the multiplication factor or the maximum value, not both")


def _fmt_row(vals, widths, fmts) -> str:
    out = []
    for v, w, f in zip(vals, widths, fmts):
        out.append(f"{v:>{w}s}" if isinstance(v, str) else f"{v:>{w}{f}}")
    return "".join(out)


def _peak(x: np.ndarray, dt: float) -> Tuple[float, float]:
    k = int(np.argmax(np.abs(x)))
    return float(abs(x[k])), k * dt


# ---------------------------------------------------------------------------------------
def run(ctx: ModuleContext) -> int:
    d = read_deck(ctx.deck_path)
    L = ctx.listing
    prof = d.rows("profile")
    nlay = len(prof)
    _validate(d, nlay)
    # ---- SOIL-NON branch: NLSOIL <Opt> = 1 in the side file <deck>.nls (AFWRITE) ------------
    nls = _read_nls(ctx)
    if nls is not None and nls.options.opt == 1:
        return _run_soil_non(ctx, d, nls)
    grav = float(d["grav"])
    form = int(d["cmodform"])
    units = _units(grav)
    indir = int(d["indir"])
    dt = float(d["delt"])
    nfft = int(d["nft"])
    cl = int(d["cl"])
    outcrop = bool(d["outcrop"])
    ratio = float(d["ratio"])
    niter = int(d["iter"])
    cutoff = float(d["soilcutoff"]) if d["soilcutoff"] > 0 else 0.0
    if dt <= 0:
        raise ModuleError("time step <delt> must be > 0")
    fnyq = 0.5 / dt

    L.section("SOIL-EQL  equivalent-linear site response (SHAKE methodology)")
    if nls is not None:
        L.write(f" {_nls_path(ctx).name}: NLSOIL <Opt> = 0 -- equivalent-linear analysis (the nonlinear "
                f"time-domain option SOIL-NON is off)")
    if d["title"]:
        L.write(f" Title            : {d['title']}")
    L.write(f" Model            : {d['model'] or ctx.model}")
    L.write(f" Units            : {units['len']}, {units['wt']}, {units['mod']}, {units['vel']}; "
            f"gravity {grav:g} {units['acc']} ({C.unit_system(grav)})")
    L.write(f" Complex modulus  : {'1 + 2 i beta (CMODFORM,1)' if form else '1 - 2 beta^2 + 2 i beta sqrt(1 - beta^2)'}")
    L.write(f" Input direction  : {'1 vertical (Vp, P-wave damping, no iterations)' if indir else '0 horizontal (Vs)'}")
    L.write(f" Iterations       : {niter}    strain ratio (effective / maximum) = {ratio:g}")
    L.write(f" Control motion   : top of sublayer {cl}, {'OUTCROP (2E)' if outcrop else 'WITHIN (E+F)'}")
    L.write(f" NFFT = {nfft}, dt = {dt:g} s, df = {1.0 / (nfft * dt):.6g} Hz, Nyquist = {fnyq:g} Hz")
    if cutoff:
        L.write(f" Cut-off frequency: {cutoff:g} Hz (EDUOPT,SOILCUTOFF)")
    else:
        L.write(" Cut-off frequency: Nyquist (<cof> is disabled, D-SOL-08)")
    if d["cof"] > 0 and abs(d["cof"] - fnyq) > 1e-9 * fnyq:
        L.write(f" (SOIL <cof> = {d['cof']:g} ignored: the cut-off is fixed to the Nyquist frequency;"
                f" use EDUOPT,SOILCUTOFF)")
    if not C.is_power_of_two(nfft):
        L.warning(f"Warning 9: NFFT = {nfft} is not a power of 2 (nearest {C.nearest_power_of_two(nfft)})")

    # ---------------- curves -----------------------------------------------------------
    try:
        curves = SH.curves_from_dynp_rows(d.rows("dynp"))
    except ValueError as exc:
        raise ModuleError(str(exc))
    for lab, cv in curves.items():
        if len(cv.g_strain) == 0:
            raise ModuleError(f"Error 97: dynamic property {lab} has no modulus curve")
        if len(cv.d_strain) == 0:
            raise ModuleError(f"Error 99: dynamic property {lab} has no damping curve")
        if np.any(cv.g_ratio < 0) or np.any(cv.g_ratio > 1):
            raise ModuleError(f"Error 98: dynamic property {lab}: G/Gmax outside [0, 1]")
        if max(len(cv.g_strain), len(cv.d_strain)) > 11:
            L.warning(f"dynamic property {lab}: more than 11 points per curve (manual limit, D-SOL-05)")
    _echo_curves(L, curves)

    # ---------------- profile ----------------------------------------------------------
    labels = [str(r["dynprop"]).strip() for r in prof]
    used = sorted({lab for lab in labels[:-1] if lab})
    for lab in used:
        if lab not in curves:
            raise ModuleError(f"Error 97/99: dynamic property '{lab}' assigned by SPRO is not defined (DYNP)")
    if indir == 0 and not used:
        L.warning("Error 95 condition: no dynamic soil property assigned to the soil sublayers; "
                  "the analysis is linear")
    thick = np.array([r["thick"] for r in prof], dtype=float)
    weight = np.array([r["weight"] for r in prof], dtype=float)
    vs0 = np.array([r["vs"] for r in prof], dtype=float)
    vp0 = np.array([r["vp"] for r in prof], dtype=float)
    ds0 = np.array([r["ds"] for r in prof], dtype=float)
    dp0 = np.array([r["dp"] for r in prof], dtype=float)
    if np.any(thick[:-1] <= 0):
        raise ModuleError("soil sublayer thickness must be > 0")
    if np.any(weight <= 0):
        raise ModuleError("unit weight must be > 0")
    if indir == 1:
        vel0, beta0 = vp0, dp0
        if np.any(vp0 <= 0):
            raise ModuleError("vertical input (indir 1) needs Vp > 0 in every sublayer")
    else:
        vel0, beta0 = vs0, ds0
        if np.any(vs0 <= 0):
            raise ModuleError("Vs must be > 0 in every sublayer")
    if np.any(beta0 < 0) or np.any(beta0 >= 0.5):
        raise ModuleError("EDU-04: damping ratios must satisfy 0 <= beta < 0.5")
    col = SH.SoilColumn(thick, weight, vel0, beta0, grav, labels if indir == 0 else [""] * nlay)
    _echo_profile(L, col, labels, units, indir)

    if indir == 1 and niter > 0:
        L.warning(f"vertical input: {niter} iterations requested, none performed (manual: set iterations = 0)")

    # ---------------- control motion ---------------------------------------------------
    if not d["thfile"]:
        raise ModuleError("Error 73: no time-history file (THFILE)")
    hpath = _resolve(ctx, d["thfile"])
    if not hpath.exists():
        raise ModuleError(f"Error 73: time-history file {d['thfile']} not found")
    header = _history_header(d, L)
    try:
        acc_g = thfile.read_soil_history(hpath, int(d["nrval"]), header)
    except ValueError as exc:
        raise ModuleError(str(exc))
    if d["nrval"] <= 0:
        L.warning(f"Error 100 condition: number of values <= 0; all {len(acc_g)} values of the file are used")
    if acc_g.size == 0:
        raise ModuleError(f"{d['thfile']}: no acceleration values read")
    if header == 0 and abs(acc_g[0] - dt) <= 1e-9 * max(dt, 1e-12) and np.max(np.abs(acc_g)) > 0:
        L.warning("the first value of the history equals the time step: if the file has the THFILE/.acc "
                  "layout (dt on line 1), set Number of Header Lines = 1")
    if len(acc_g) > nfft:
        raise ModuleError(f"EDU-03: {len(acc_g)} values do not fit in NFFT = {nfft}; use NFFT = "
                          f"{1 << int(math.ceil(math.log2(len(acc_g))))}")
    pk0, tk0 = _peak(acc_g, dt)
    acc_g = thfile.scale_history(acc_g, float(d["mult"]), float(d["max"]))
    pk1, _ = _peak(acc_g, dt)
    L.section("Input motion")
    L.write(f" File               : {d['thfile']}   {d['thtit']}")
    L.write(f" Values read        : {len(acc_g)} after {header} header line(s); dt = {dt:g} s")
    L.write(f" Maximum acceleration = {pk0:.5f} g at t = {tk0:.2f} s")
    fac = pk1 / pk0 if pk0 > 0 else 0.0
    L.write(f" Multiplied by {fac:.5f} to give a maximum of {pk1:.5f} g "
            f"({'MULT' if d['mult'] else 'MAX'} scaling)")

    if int(d["opmode"]) == 1:
        L.write("")
        L.write(" Operation mode 1: data check only -- input read and checked, no analysis performed")
        return 0

    # ---------------- analysis ---------------------------------------------------------
    ctx.progress(0.05, "SOIL: equivalent-linear iterations")
    try:
        res = SH.run_shake(acc_g * grav, dt, nfft, col, curves, cl, outcrop, ratio, niter, form=form,
                           cutoff=cutoff, iterate=(indir == 0))
    except (ValueError, FloatingPointError) as exc:
        raise ModuleError(str(exc))
    a_in = np.fft.irfft(res.input_spectrum, n=nfft) / grav
    pk2, _ = _peak(a_in, dt)
    A = np.abs(res.input_spectrum) ** 2
    msf = float(np.sum(res.freqs * A) / np.sum(A)) if np.sum(A) > 0 else 0.0
    L.write(f" Maximum acceleration after the cut-off = {pk2:.5f} g;  mean square frequency = {msf:.2f} Hz")

    _initial_summary(L, col, form, cl, outcrop, res, cutoff or fnyq, units)
    _iteration_tables(L, col, labels, res, units)
    ctx.progress(0.6, "SOIL: output")

    # ---------------- final properties, strains, stresses --------------------------------
    G = res.G
    beta = res.beta
    vel = np.sqrt(G / col.rho)
    nsoil = nlay - 1
    L.section("Strain-compatible soil properties (final)")
    g_compat = res.gamma_eff_compatible
    hdr = ["NO", "LABEL", "DEPTH", "EFF.STRAIN%", "MAX.STRAIN%", "DAMPING", "G" if indir == 0 else "M",
           "G/GMAX" if indir == 0 else "M/M0", "VS" if indir == 0 else "VP", "MAX.STRESS", "TIME"]
    w = [4, 10, 9, 13, 13, 9, 12, 8, 10, 12, 8]
    L.write(_fmt_row(hdr, w, [""] * len(w)))
    stress_max = np.zeros(nsoil)
    stress_t = np.zeros(nsoil)
    for i in range(nsoil):
        tau = res.stress_history(i + 1, form)
        stress_max[i], stress_t[i] = _peak(tau, dt)
        L.write(_fmt_row([i + 1, labels[i] or "-", col.depth_mid[i], g_compat[i], res.gamma_max[i],
                          beta[i], G[i], G[i] / col.gmax[i], vel[i], stress_max[i], stress_t[i]],
                         w, ["d", "s", ".2f", ".5f", ".5f", ".4f", ".1f", ".3f", ".1f", ".4g", ".2f"]))
    L.write(f" (stress in {units['mod']}, computed as tau(w) = G*(w) gamma(w); depth to mid-height in {units['len']})")
    if res.iterations:
        L.write(" EFF.STRAIN% is the effective strain of the last iteration, from which DAMPING and G were read "
                "(strain-compatible pairs, as in FILE88);")
        L.write(" MAX.STRAIN%, MAX.STRESS and TIME are those of the final response computed with these properties")
        d_eff = np.max(np.abs(res.gamma_eff - g_compat) / np.where(g_compat > 0, g_compat, 1.0)) * 100.0
        L.write(f" (effective strain of the final response differs by at most {d_eff:.2f} % from EFF.STRAIN%)")
    v_avg = col.average_velocity(vel)
    h_tot = float(np.sum(thick[:-1]))
    L.write(f" Period = {4 * h_tot / v_avg:.2f} s from the average {'shear' if indir == 0 else 'P'} "
            f"velocity {v_avg:.1f} {units['vel']} (4H/V)")
    _amplification(L, col, G, beta, form, cutoff or fnyq, res, "final")

    # ---------------- accelerations ------------------------------------------------------
    files: List[str] = []
    sacc = d.rows("sacc")
    if sacc:
        L.section("Maximum accelerations (top of sublayer)")
        L.write(_fmt_row(["LAYER", "DEPTH", "MOTION", "MAX.ACC(g)", "TIME(s)", "MEAN.SQ.FR", "FILE"],
                         [6, 10, 9, 12, 9, 11, 12], [""] * 7))
        for r in sacc:
            lay, opt, oc = int(r["layer"]), int(r["opt"]), bool(r["outcrop"])
            if opt <= 0:
                continue
            if not (1 <= lay <= nlay):
                L.warning(f"SACC layer {lay} is not a sublayer (1..{nlay}); skipped")
                continue
            a = res.acc_history(lay, oc) / grav
            pk, tp = _peak(a, dt)
            S = np.abs(res.acc_spectrum(lay, oc)) ** 2
            ms = float(np.sum(res.freqs * S) / np.sum(S)) if np.sum(S) > 0 else 0.0
            fname = ""
            if opt >= 2:
                fname = C.layer_th_name("ACC", lay)
                textfiles.write_history(ctx.path(fname), a, dt)
                files.append(fname)
            L.write(_fmt_row([lay, col.depth_top[lay - 1], "OUTCROP" if oc else "WITHIN", pk, tp, ms, fname],
                             [6, 10, 9, 12, 9, 11, 12], ["d", ".2f", "s", ".5f", ".2f", ".2f", "s"]))

    # ---------------- stresses / strains ------------------------------------------------
    sstr = d.rows("sstr")
    if sstr:
        L.section("Stress and strain at mid-height of sublayers")
        for r in sstr:
            lay = int(r["layer"])
            if not (1 <= lay <= nsoil):
                L.warning(f"SSTR layer {lay} is not a soil sublayer (1..{nsoil}); skipped")
                continue
            if r["opt3"] or r["opt4"]:
                g = res.strain_history(lay) * 100.0
                pk, tp = _peak(g, dt)
                L.write(f" layer {lay:3d}: maximum strain = {pk:.5f} % at t = {tp:.2f} s")
                if r["opt4"]:
                    fname = C.layer_th_name("SN", lay)
                    textfiles.write_history(ctx.path(fname), g, dt)
                    files.append(fname)
            if r["opt1"] or r["opt2"]:
                tau = res.stress_history(lay, form)
                pk, tp = _peak(tau, dt)
                L.write(f" layer {lay:3d}: maximum stress = {pk:.6g} {units['mod']} at t = {tp:.2f} s")
                if r["opt2"]:
                    fname = C.layer_th_name("SS", lay)
                    textfiles.write_history(ctx.path(fname), tau, dt)
                    files.append(fname)

    # ---------------- response spectra --------------------------------------------------
    damps = [float(r["value"]) for r in d.rows("damp")]
    srs = d.rows("srs")
    saf = d.rows("ssaf")
    rs_cache: Dict[Tuple[int, bool], Dict[str, np.ndarray]] = {}

    def spectrum(lay: int, oc: bool) -> Dict[str, np.ndarray]:
        key = (lay, oc)
        if key not in rs_cache:
            rs_cache[key] = SP.response_spectrum(res.acc_history(lay, oc) / grav, dt, RS_FREQS, damps)
        return rs_cache[key]

    if srs and not damps:
        raise ModuleError("Error 107: no damping ratios (DAMP) for the response spectra")
    if srs and damps:
        L.section("Response spectra at the top of sublayers (absolute acceleration, g x gravmult)")
        L.write(f" {len(RS_FREQS)} frequencies, 0.1 - 100 Hz (100 per decade); dampings {damps}")
        for r in srs:
            lay, oc = int(r["layer"]), bool(r["outcrop"])
            if not (1 <= lay <= nlay):
                L.warning(f"SRS layer {lay} is not a sublayer; skipped")
                continue
            rs = spectrum(lay, oc)
            for j, z in enumerate(damps):
                sa = rs["SA"][j] * float(d["gravmult"])
                k = int(np.argmax(sa))
                fname = ""
                if int(r["save"]):
                    fname = soil_rs_name(lay, j + 1)
                    textfiles.write_xy(ctx.path(fname), RS_FREQS, sa,
                                       header=f"SOIL RS layer {lay} {'outcrop' if oc else 'within'} damping {z:g}"
                                              f" (f Hz, SA g x {d['gravmult']:g})")
                    files.append(fname)
                L.write(f" layer {lay:3d} {'OUTCROP' if oc else 'WITHIN '} damping {z:5.3f}: "
                        f"max SA = {sa[k]:.4f} at {RS_FREQS[k]:.3f} Hz (T = {1 / RS_FREQS[k]:.3f} s); "
                        f"SA(100 Hz) = {sa[-1]:.4f}  {fname}")

    # ---------------- spectral amplification factors ----------------------------------
    if saf:
        L.section("Spectral amplification factors (SSAF)")
        fmax = cutoff or fnyq
        for r in saf:
            lay, lay2 = int(r["layer"]), int(r["layer2"])
            oc1, oc2 = bool(r["outcrop1"]), bool(r["outcrop2"])
            fst = float(r["freqstep"])
            if not (1 <= lay <= nlay) or not (1 <= lay2 <= nlay):
                L.warning(f"Error 109: illegal SSAF layers {lay}/{lay2}; skipped")
                continue
            if fst <= 0:
                L.warning(f"Error 110: SSAF frequency step {fst:g} <= 0; skipped")
                continue
            npts = int(math.floor(fmax / fst + 1e-9)) + 1
            if npts > MAX_SAF_POINTS:
                L.warning(f"SSAF: {npts} frequencies requested, truncated to {MAX_SAF_POINTS}")
                npts = MAX_SAF_POINTS
            f = np.arange(npts) * fst
            H = SH.column_transfer(f, col, G, beta, lay, oc1, lay2, oc2, form)
            k = int(np.argmax(np.abs(H)))
            title = r["title"] or ""
            L.write(f" {title}")
            L.write(f"  X({lay}, {'outcrop' if oc1 else 'within'}) / X({lay2}, {'outcrop' if oc2 else 'within'}),"
                    f" df = {fst:g} Hz: max |amplification| = {abs(H[k]):.4f} at {f[k]:.4f} Hz")
            if int(r["save"]):
                fname = saf_tf_name(lay, oc1, lay2, oc2)
                textfiles.write_tf(ctx.path(fname), f, H, complex_=True,
                                   header=f"SOIL SSAF Fourier amplification {title}")
                files.append(fname)
            if damps:
                # spectral amplification factor RS(layer)/RS(layer2) sampled at f_k = k * freqstep,
                # k >= 1, up to the cut-off (spec 05a 6.5, command spec 9.2.35)
                fr = f[1:]
                if len(fr) > MAX_SAF_RS_POINTS:
                    L.warning(f"SSAF: {len(fr)} response-spectrum frequencies requested, truncated to "
                              f"{MAX_SAF_RS_POINTS} (up to {fr[MAX_SAF_RS_POINTS - 1]:g} Hz)")
                    fr = fr[:MAX_SAF_RS_POINTS]
                if fr.size == 0:
                    L.warning(f"SSAF: frequency step {fst:g} Hz exceeds the cut-off {fmax:g} Hz; no RS ratio")
                else:
                    a1 = res.acc_history(lay, oc1) / grav
                    a2 = res.acc_history(lay2, oc2) / grav
                    r1 = SP.response_spectrum(a1, dt, fr, damps)
                    r2 = SP.response_spectrum(a2, dt, fr, damps)
                    for j, z in enumerate(damps):
                        rat = r1["SA"][j] / np.where(r2["SA"][j] > 0, r2["SA"][j], np.inf)
                        kk = int(np.argmax(rat))
                        L.write(f"  RS ratio damping {z:5.3f}, f = k x {fst:g} Hz ({len(fr)} values): max "
                                f"{rat[kk]:.4f} at {fr[kk]:.3f} Hz")
                        if int(r["save"]):
                            fname = saf_rs_name(lay, oc1, lay2, oc2, j + 1)
                            textfiles.write_xy(ctx.path(fname), fr, rat,
                                               header=f"SOIL SSAF RS ratio {title} damping {z:g} "
                                                      f"(f = k x {fst:g} Hz, RS({lay})/RS({lay2}))")
                            files.append(fname)

    # ---------------- FILE73, FILE88 -----------------------------------------------------
    if curves:
        SH.write_file73(ctx.path("FILE73"), curves, title=d["title"])
        files.append("FILE73")
    if int(d["save"]):
        if indir == 0:
            vs = vel
            bs = beta.copy()
            # D-SOL-06 applies to the iterated sublayers only: Vp at constant Poisson's ratio and
            # beta_p = beta_s.  Linear sublayers (no DYNP label) and the half-space keep their input
            # Vp and P-wave damping.
            itr = np.array([bool(lab) for lab in labels[:-1]] + [False])
            if str(d["vppolicy"]).upper() == "VP":     # EDUOPT,VPPOLICY,VP: keep input Vp and beta_p
                itr = np.zeros_like(itr)
            vp = np.where(itr, vp0 * vs / np.where(vs0 > 0, vs0, 1.0), vp0)
            bp = np.where(itr, bs, dp0)
        else:
            vs, bs, vp, bp = vs0.copy(), ds0.copy(), vel, beta.copy()
        Gs = col.rho * vs ** 2
        hs_row = [thick[-1], weight[-1], vp[-1], vs[-1], bp[-1], bs[-1]]
        SH.write_file88(ctx.path("FILE88"), thick[:-1], g_compat, Gs[:-1], vs[:-1], bs[:-1], vp[:-1],
                        bp[:-1], units=f"{units['len']} {units['mod']} {units['vel']}; damping ratios; strain %",
                        halfspace=hs_row, title=d["title"])
        files.append("FILE88")
        L.write("")
        L.write(" FILE88 written: strain-compatible properties of the "
                f"{nsoil} soil sublayers (iterated sublayers: Vp at constant Poisson's ratio, beta_p = beta_s, "
                "D-SOL-06; linear sublayers keep their input Vp and beta_p)" if indir == 0 else
                " FILE88 written: P-wave properties of the vertical analysis; Vs and beta_s as input")
    L.section("Files written")
    for fn in files:
        L.write(f"  {fn}")
    ctx.progress(1.0, "SOIL: done")
    return 0


# ---------------------------------------------------------------------------------------
# listing pieces
# ---------------------------------------------------------------------------------------
def _echo_curves(L, curves: Dict[str, SH.DynamicProperty]) -> None:
    L.section("Dynamic soil properties (DYNP): strain %, G/Gmax, damping %")
    for k, (lab, cv) in enumerate(curves.items(), start=1):
        L.write(f" Curve {k}: {lab}")
        L.write("     STRAIN %   G/GMAX       STRAIN %  DAMPING %")
        n = max(len(cv.g_strain), len(cv.d_strain))
        for i in range(n):
            a = f"{cv.g_strain[i]:12.5g}{cv.g_ratio[i]:9.4f}" if i < len(cv.g_strain) else " " * 21
            b = f"{cv.d_strain[i]:15.5g}{cv.d_pct[i]:11.3f}" if i < len(cv.d_strain) else ""
            L.write(a + b)


def _echo_profile(L, col: SH.SoilColumn, labels, units, indir) -> None:
    L.section("Soil profile (SPRO): last entry is the half-space")
    vlab = "VS" if indir == 0 else "VP"
    hdr = ["NO", "LABEL", "THICK", "DEPTH.MID", "UNIT.WT", vlab, "GMAX" if indir == 0 else "M0", "DAMPING"]
    w = [4, 10, 10, 11, 10, 10, 12, 9]
    L.write(_fmt_row(hdr, w, [""] * len(w)))
    for i in range(col.n):
        if i < col.n - 1:
            L.write(_fmt_row([i + 1, labels[i] or "-", col.thick[i], col.depth_mid[i], col.weight[i], col.vel0[i],
                              col.gmax[i], col.beta0[i]], w,
                             ["d", "s", ".3f", ".3f", ".4f", ".1f", ".1f", ".4f"]))
        else:
            L.write(_fmt_row([i + 1, "BASE", "-", "-", col.weight[i], col.vel0[i], col.gmax[i], col.beta0[i]], w,
                             ["d", "s", "s", "s", ".4f", ".1f", ".1f", ".4f"]))
    L.write(f" Depth to the half-space: {float(np.sum(col.thick[:-1])):.3f} {units['len']}; "
            f"mass density = unit weight / gravity")


def _amplification(L, col, G, beta, form, fmax, res: SH.ShakeResult, tag: str) -> None:
    """SHAKE-style maximum amplification surface / within base, plus outcrop."""
    n = col.n
    f = res.freqs[(res.freqs > 0) & (res.freqs <= fmax + 1e-12)]
    Hw = np.abs(SH.column_transfer(f, col, G, beta, 1, True, n, False, form))
    k = int(np.argmax(Hw))
    a, fa = SH.max_amplification(col, G, beta, 1, True, n, False, fmax, form)
    Ho = np.abs(SH.column_transfer(f, col, G, beta, 1, True, n, True, form))
    ko = int(np.argmax(Ho))
    b, fb = SH.max_amplification(col, G, beta, 1, True, n, True, fmax, form)
    L.write(f" Maximum amplification surface / base WITHIN ({tag}) = {Hw[k]:.4f} at {f[k]:.4f} Hz "
            f"(T = {1 / f[k]:.3f} s) on the FFT grid; continuous peak {a:.4f} at {fa:.4f} Hz")
    L.write(f" Maximum amplification surface / base OUTCROP ({tag}) = {Ho[ko]:.4f} at {f[ko]:.4f} Hz "
            f"(T = {1 / f[ko]:.3f} s) on the FFT grid; continuous peak {b:.4f} at {fb:.4f} Hz")


def _initial_summary(L, col, form, cl, outcrop, res, fmax, units) -> None:
    L.section("Low-strain profile")
    v_avg = col.average_velocity()
    h = float(np.sum(col.thick[:-1]))
    L.write(f" Period = {4 * h / v_avg:.2f} s from the average velocity {v_avg:.1f} {units['vel']} (4H/V)")
    _amplification(L, col, col.gmax, col.beta0, form, fmax, res, "initial")


def _iteration_tables(L, col, labels, res: SH.ShakeResult, units) -> None:
    if not res.iterations:
        L.section("No equivalent-linear iteration (linear analysis with the low-strain properties)")
        return
    for rec in res.iterations:
        L.section(f"ITERATION NUMBER {rec.number:3d}")
        L.write("   NO  LABEL       DEPTH  EFF.STRAIN%  <------- DAMPING ------>  <--------- SHEAR MODULUS --------->  G/GMAX")
        L.write("                                       NEW     USED   ERROR%          NEW         USED   ERROR%")
        dG, dB = rec.dG_pct, rec.dbeta_pct
        for i in range(len(rec.gamma_eff)):
            L.write(f" {i + 1:4d}  {labels[i] or '-':<10s}{col.depth_mid[i]:7.2f} {rec.gamma_eff[i]:11.5f} "
                    f"{rec.beta_new[i]:8.4f} {rec.beta_used[i]:8.4f} {dB[i]:8.1f} "
                    f"{rec.G_new[i]:12.1f} {rec.G_used[i]:12.1f} {dG[i]:8.1f} {rec.G_new[i] / col.gmax[i]:7.3f}")
        L.write(f" max |error|: damping {np.max(np.abs(dB)):.2f} %, modulus {np.max(np.abs(dG)):.2f} %")
    last = res.iterations[-1]
    L.write("")
    L.write(f" {len(res.iterations)} iteration(s) performed; last changes: damping "
            f"{np.max(np.abs(last.dbeta_pct)):.2f} %, modulus {np.max(np.abs(last.dG_pct)):.2f} % "
            f"(the final response uses the properties of the last update)")


# =======================================================================================
# SOIL-NON: nonlinear time-domain site response (NLSOIL <Opt> = 1)
# =======================================================================================
def _nls_path(ctx: ModuleContext) -> Path:
    """Side file of the NLSOIL / NLSLAYER data: ``<deck stem>.nls`` next to the SOIL deck (AFWRITE)."""
    if ctx.deck_path is not None:
        return Path(ctx.deck_path).with_suffix(".nls")
    return ctx.path(f"{ctx.model}.nls")


def _read_nls(ctx: ModuleContext) -> Optional[SN.NlsData]:
    p = _nls_path(ctx)
    if not p.exists():
        return None
    try:
        return SN.read_nls(p)
    except (OSError, ValueError) as exc:
        raise ModuleError(f"{p.name}: {exc}")


def _deck_hash(path: Optional[Path]) -> str:
    if path is None:
        return ""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for _ in range(40):
                ln = fh.readline()
                if not ln:
                    break
                m = re.search(r"model_hash\s*=\s*(\S+)", ln)
                if m:
                    return m.group(1)
    except OSError:
        pass
    return ""


def _read_control_motion(ctx: ModuleContext, d: decks.Deck, L, dt: float, nfft: int) -> np.ndarray:
    """Control motion (g) read and scaled as in SOIL-EQL (THFILE, header, nrval, mult / max)."""
    if not d["thfile"]:
        raise ModuleError("Error 73: no time-history file (THFILE)")
    hpath = _resolve(ctx, d["thfile"])
    if not hpath.exists():
        raise ModuleError(f"Error 73: time-history file {d['thfile']} not found")
    header = _history_header(d, L)
    try:
        acc_g = thfile.read_soil_history(hpath, int(d["nrval"]), header)
    except ValueError as exc:
        raise ModuleError(str(exc))
    if d["nrval"] <= 0:
        L.warning(f"Error 100 condition: number of values <= 0; all {len(acc_g)} values of the file are used")
    if acc_g.size == 0:
        raise ModuleError(f"{d['thfile']}: no acceleration values read")
    if header == 0 and abs(acc_g[0] - dt) <= 1e-9 * max(dt, 1e-12) and np.max(np.abs(acc_g)) > 0:
        L.warning("the first value of the history equals the time step: if the file has the THFILE/.acc "
                  "layout (dt on line 1), set Number of Header Lines = 1")
    if len(acc_g) > nfft:
        raise ModuleError(f"EDU-03: {len(acc_g)} values do not fit in NFFT = {nfft}; use NFFT = "
                          f"{1 << int(math.ceil(math.log2(len(acc_g))))}")
    pk0, tk0 = _peak(acc_g, dt)
    acc_g = thfile.scale_history(acc_g, float(d["mult"]), float(d["max"]))
    pk1, _ = _peak(acc_g, dt)
    L.section("Input motion (at bedrock)")
    L.write(f" File               : {d['thfile']}   {d['thtit']}")
    L.write(f" Values read        : {len(acc_g)} after {header} header line(s); dt = {dt:g} s")
    L.write(f" Maximum acceleration = {pk0:.5f} g at t = {tk0:.2f} s")
    fac = pk1 / pk0 if pk0 > 0 else 0.0
    L.write(f" Multiplied by {fac:.5f} to give a maximum of {pk1:.5f} g "
            f"({'MULT' if d['mult'] else 'MAX'} scaling)")
    return acc_g


def _run_soil_non(ctx: ModuleContext, d: decks.Deck, nls: SN.NlsData) -> int:
    """SOIL-NON analysis (see the module docstring and :mod:`sassi.core.soilnon`, decisions SN-1 ... SN-12)."""
    L = ctx.listing
    o = nls.options
    prof = d.rows("profile")
    nlay = len(prof)
    nsoil = nlay - 1
    grav = float(d["grav"])
    units = _units(grav)
    dt = float(d["delt"])
    nfft = int(d["nft"])
    cl = int(d["cl"])
    outcrop = bool(d["outcrop"])
    ratio = float(d["ratio"])
    cutoff = float(d["soilcutoff"]) if d["soilcutoff"] > 0 else 0.0
    indir = int(d["indir"])
    if dt <= 0:
        raise ModuleError("time step <delt> must be > 0")
    fnyq = 0.5 / dt
    nls_name = _nls_path(ctx).name

    L.section("SOIL-NON  nonlinear time-domain site response (MKZ hyperbolic model, extended Masing rules)")
    L.write(" WARNING (manual): not for direct use in nuclear licensing site response unless fully justified;")
    L.write("                   a benchmarking tool for soft soils with strong nonlinearity (DEEPSOIL-like).")
    L.write(" WARNING (manual): valid only for nonlinear convolution with the input motion defined as the")
    L.write("                   compatible motion at bedrock.")
    if d["title"]:
        L.write(f" Title            : {d['title']}")
    L.write(f" Model            : {d['model'] or ctx.model}")
    L.write(f" Units            : {units['len']}, {units['wt']}, {units['mod']}, {units['vel']}; "
            f"gravity {grav:g} {units['acc']} ({C.unit_system(grav)})")
    L.write(f" Options file     : {nls_name} (NLSOIL / NLSLAYER written by AFWRITE)")
    h_deck, h_nls = _deck_hash(ctx.deck_path), nls.model_hash
    if h_deck and h_nls and h_deck != h_nls:
        L.warning(f"{nls_name} was not written together with {Path(ctx.deck_path).name} (different model_hash): "
                  f"run AFWRITE again")
    probs = o.problems()
    if probs:
        raise ModuleError("; ".join(probs))
    if indir != 0:
        raise ModuleError("SOIL-NON models vertically propagating shear waves: Input Direction must be 0 "
                          "(SOILX <indir>)")
    if cl != nlay:
        raise ModuleError(f"SOIL-NON needs the control motion at bedrock: the control layer ({cl}) must be the "
                          f"half-space, sublayer {nlay} = last SPRO entry (SOILX <cl>, D-SOL-13); manual: valid "
                          f"only with the input motion defined at bedrock")
    bed = int(o.bedint)
    if bed == 0 and outcrop:
        L.warning("BedInt 0 (rigid bedrock): the control motion is applied as the motion of the base itself "
                  "(within motion); SOIL <outcrop> = 1 is ignored -- use BedInt 1 for an outcrop motion")
    if bed == 1 and not outcrop:
        raise ModuleError("BedInt 1 (viscoelastic bedrock) needs an outcrop control motion (SOIL <outcrop> = 1); "
                          "for a within (borehole) motion use BedInt 0 (rigid bedrock)")
    ctl = SN.Controls.from_options(o, dt)
    dtype = int(o.damptype) if int(o.damptype) in (1, 2, 3) else 1

    L.section("Nonlinear soil options (NLSOIL)")
    L.write(f" NLSOIL,{','.join(SN.format_value(v) for v in o.values())}")
    if ctl.flexible:
        L.write(f" Sub-increments     : 0 = flexible -- {ctl.nsub} base sub-step(s) per time step (<= {SN.DT_MAX_FLEX:g} s),"
                f" refined when a strain increment exceeds {100 * ctl.dg_max:g} % (SN-9)")
    else:
        L.write(f" Sub-increments     : {ctl.nsub} per time step (fixed), sub-step {dt / ctl.nsub:g} s")
        if dt / ctl.nsub > SN.DT_MAX_FLEX * (1 + 1e-9):
            L.warning(f"sub-step {dt / ctl.nsub:g} s > {SN.DT_MAX_FLEX:g} s: the Newmark period elongation exceeds "
                      f"1.3 % below 25 Hz (use more sub-increments or 0 = flexible)")
    L.write(f" Convergence        : displacement {ctl.tol_d:g}{' (default)' if o.dispconv <= 0 else ''}, force "
            f"{ctl.tol_f:g}{' (default)' if o.forceconv <= 0 else ''} (relative 2-norms, both required, SN-10)")
    L.write(f" Equilibrium iter.  : {ctl.maxit}{' (default)' if o.equalit <= 0 else ''} per sub-step "
            f"(Newton-Raphson, consistent tangent; bisection when exceeded)")
    L.write(f" Bedrock interface  : {bed} {SN.BEDROCK_TYPES[bed]}")
    if bed == 1:
        L.write("                      (BedInt 1 is 'currently disabled' in ACS SASSI V3; SASSI-EDU implements it with "
                "the Joyner & Chen 1975 dashpot, SN-7)")
    L.write(f" Damping type       : {o.damptype} -> {dtype} {SN.DAMPING_TYPES[dtype]}")
    if dtype == 3:
        L.write(f" Rayleigh multipliers: mass {o.mmmult:g}, stiffness {o.smmult:g}")
    L.write(" Time integration   : Newmark average acceleration (gamma 1/2, beta 1/4), relative to the input frame")
    if int(d["iter"]) > 0:
        L.write(f" (SOIL <iter> = {d['iter']} is not used: no equivalent-linear iterations in SOIL-NON; "
                f"<ratio> = {ratio:g} only defines the effective strain of FILE88)")

    # ---------------- curves and profile (as SOIL-EQL) ----------------------------------
    try:
        curves = SH.curves_from_dynp_rows(d.rows("dynp"))
    except ValueError as exc:
        raise ModuleError(str(exc))
    for lab, cv in curves.items():
        if len(cv.g_strain) == 0:
            raise ModuleError(f"Error 97: dynamic property {lab} has no modulus curve")
        if len(cv.d_strain) == 0:
            raise ModuleError(f"Error 99: dynamic property {lab} has no damping curve")
        if np.any(cv.g_ratio < 0) or np.any(cv.g_ratio > 1):
            raise ModuleError(f"Error 98: dynamic property {lab}: G/Gmax outside [0, 1]")
    if curves:
        _echo_curves(L, curves)
    labels = [str(r["dynprop"]).strip() for r in prof]
    for lab in sorted({lab for lab in labels[:-1] if lab}):
        if lab not in curves:
            raise ModuleError(f"Error 97/99: dynamic property '{lab}' assigned by SPRO is not defined (DYNP)")
    thick = np.array([r["thick"] for r in prof], dtype=float)
    weight = np.array([r["weight"] for r in prof], dtype=float)
    vs0 = np.array([r["vs"] for r in prof], dtype=float)
    vp0 = np.array([r["vp"] for r in prof], dtype=float)
    ds0 = np.array([r["ds"] for r in prof], dtype=float)
    dp0 = np.array([r["dp"] for r in prof], dtype=float)
    if np.any(thick[:-1] <= 0):
        raise ModuleError("soil sublayer thickness must be > 0")
    if np.any(weight <= 0):
        raise ModuleError("unit weight must be > 0")
    if np.any(vs0 <= 0):
        raise ModuleError("Vs must be > 0 in every sublayer and in the half-space")
    if np.any(ds0 < 0) or np.any(ds0 >= 0.5):
        raise ModuleError("EDU-04: damping ratios must satisfy 0 <= beta < 0.5")
    col_sh = SH.SoilColumn(thick, weight, vs0, ds0, grav, labels)
    _echo_profile(L, col_sh, labels, units, 0)
    rho = weight / grav

    # ---------------- nonlinear model of every sublayer (SN-1 ... SN-5) -----------------
    models: List[SN.SublayerModel] = []
    for j in range(nsoil):
        lab = labels[j]
        try:
            models.append(SN.resolve_sublayer(j + 1, nls.layers.get(j + 1), curves.get(lab) if lab else None,
                                              ds0[j]))
        except (SN.SoilNonError, ValueError) as exc:
            raise ModuleError(str(exc))
    extra = sorted(k for k in nls.layers if k > nsoil)
    fdisc = min(SN.F_MAX_DISCRETISATION, fnyq)
    nel = SN.elements_per_sublayer(thick[:-1], vs0[:-1], fdisc)
    col = SN.SoilNonColumn(thick[:-1], rho[:-1], vs0[:-1], models, float(rho[-1]), float(vs0[-1]), nel)
    w_fb, _ = SN.fixed_base_modes(col)
    w1 = float(w_fb[0])
    if dtype == 2:
        for j, mdl in enumerate(models):
            if nls.layers[j + 1].curvefit == 1:
                mdl.eta = 2.0 * mdl.xi_min * rho[j] * vs0[j] ** 2 / w1
        col = SN.SoilNonColumn(thick[:-1], rho[:-1], vs0[:-1], models, float(rho[-1]), float(vs0[-1]), nel)
    damping = SN.build_damping(col, o.damptype, bed, o.mmmult, o.smmult)

    L.section("Nonlinear soil layers (NLSLAYER): tau = G0 g / (1 + beta (g/g_r)^s) [+ eta dg/dt]")
    hdr = ["NO", "LABEL", "MODEL", "BETA", "S", "G_REF %", "FIT RMS", "XI_MIN", "ETA", "ELEMS", "H_ELEM", "FMAX"]
    w = [4, 10, 8, 8, 8, 11, 9, 9, 11, 6, 9, 8]
    L.write(_fmt_row(hdr, w, [""] * len(w)))
    for j, mdl in enumerate(models):
        L.write(_fmt_row([j + 1, labels[j] or "-", mdl.source, mdl.beta, mdl.s,
                          mdl.gr_pct if math.isfinite(mdl.gr_pct) else float("inf"),
                          mdl.fit.rms if mdl.fit is not None else 0.0, mdl.xi_min, mdl.eta, int(nel[j]),
                          thick[j] / nel[j], col.deepsoil_fmax()[j]],
                         w, ["d", "s", "s", ".4f", ".4f", ".5g", ".4f", ".4f", ".4g", "d", ".3f", ".1f"]))
    L.write(" MODEL: user = NLSLAYER parameters, fit = MKZ fitted to the DYNP G/Gmax curve (beta = 1, SN-3), linear =")
    L.write("        beta 0; G_REF in %; FIT RMS = rms error of G/Gmax at the curve points; XI_MIN small-strain")
    L.write("        damping (SN-4); ETA viscosity (NLDampType 2); ELEMS / H_ELEM: shear elements of the sublayer")
    L.write(f"        (h <= Vs / ({SN.ELEMENTS_PER_WAVELENGTH} x {fdisc:g} Hz), SN-11); FMAX = Vs / (4 h) (DEEPSOIL check)")
    for j, mdl in enumerate(models):
        if mdl.source == "linear" and nls.layers[j + 1].curvefit == 0:
            L.warning(f"NLSLAYER {j + 1}: Beta = 0 -- linear elastic sublayer (SN-2)")
        if mdl.source == "user" and mdl.s > 1.0:
            L.warning(f"NLSLAYER {j + 1}: S = {mdl.s:g} > 1 -- the backbone stress decreases beyond "
                      f"{mdl.gr_pct * (1.0 / (mdl.beta * (mdl.s - 1.0))) ** (1.0 / mdl.s):.4g} % strain")
        if dtype == 2 and mdl.eta == 0.0:
            L.warning(f"NLDampType 2: sublayer {j + 1} has no viscosity (NLSLAYER Vis = 0): no small-strain damping")
    if extra:
        L.write(f" NLSLAYER set(s) {extra} beyond the {nsoil} soil sublayers (the half-space is the elastic base) "
                f"are not used (SN-1)")
    _soil_non_damping_listing(L, col, damping, w_fb, models, dtype)

    # ---------------- small-strain amplification: discrete model vs continuum ----------
    fmax_tf = min(fdisc, fnyq)
    ftf = np.arange(0.05, fmax_tf + 1e-9, 0.05)
    if ftf.size:
        Hd = np.abs(SN.linear_transfer(col, damping, bed, ftf, nodes=[0])[:, 0])
        xi_c = np.array([m.xi_min for m in models] + [0.0])
        Hc = np.abs(SH.column_transfer(ftf, col_sh, col_sh.gmax, xi_c, 1, False, nlay, bed == 1, 0))
        kd, kc = int(np.argmax(Hd)), int(np.argmax(Hc))
        L.section("Small-strain amplification surface / " + ("outcrop" if bed == 1 else "base (within)"))
        L.write(f" discrete SOIL-NON model (linear, damping type {dtype}): max {Hd[kd]:.4f} at {ftf[kd]:.2f} Hz")
        L.write(f" continuum (SHAKE recursion, damping XI_MIN, half-space undamped): max {Hc[kc]:.4f} at "
                f"{ftf[kc]:.2f} Hz  (0.05 Hz grid up to {fmax_tf:g} Hz)")
        L.write(f" fundamental frequency of the column on a rigid base f1 = {w1 / (2 * np.pi):.4f} Hz")

    acc_g = _read_control_motion(ctx, d, L, dt, nfft)
    if int(d["opmode"]) == 1:
        L.write("")
        L.write(" Operation mode 1: data check only -- input read and checked, no analysis performed")
        return 0

    # ---------------- time integration ---------------------------------------------------
    ctx.progress(0.05, "SOIL-NON: time integration")
    try:
        res = SN.run_soilnon(acc_g * grav, dt, nfft, col, damping, bed, ctl, cutoff=cutoff,
                             progress=lambda f: ctx.progress(0.05 + 0.85 * f, "SOIL-NON: time integration"),
                             cancelled=ctx.cancelled)
    except SN.SoilNonError as exc:
        raise ModuleError(f"SOIL-NON: {exc}")
    st = res.stats
    L.section("Time integration")
    L.write(f" {nfft} time steps of {dt:g} s ({nfft * dt:g} s: the record of {len(acc_g)} values is zero-padded "
            f"to NFFT, as SOIL-EQL)")
    if cutoff:
        L.write(f" Input low-pass at {cutoff:g} Hz (EDUOPT,SOILCUTOFF), band-limited interpolation for the sub-steps (SN-8)")
    else:
        L.write(" Input interpolated band-limited (FFT) for the sub-steps (SN-8)")
    L.write(f" sub-steps {int(st['substeps'])}, Newton iterations {int(st['iterations'])} "
            f"(average {st['iterations'] / max(st['substeps'], 1):.2f}, largest {int(st['max_iterations'])} per sub-step)")
    L.write(f" smallest sub-step {st['min_substep']:.4g} s; bisections after non-convergence {int(st['bisections'])}; "
            f"sub-steps repeated for the strain-increment limit {int(st['strain_rejections'])}")
    L.write(f" cpu {st['cpu_s']:.2f} s")
    ctx.progress(0.92, "SOIL-NON: output")
    return _soil_non_outputs(ctx, d, L, res, col, models, labels, col_sh, thick, weight, vs0, vp0, ds0, dp0,
                             curves, units, bed, ratio, cutoff or fnyq)


def _soil_non_damping_listing(L, col: SN.SoilNonColumn, damping: SN.DampingModel, w_fb: np.ndarray,
                              models: List[SN.SublayerModel], dtype: int) -> None:
    L.section("Small-strain viscous damping (NLDampType)")
    for note in damping.notes:
        L.write(f" {note}")
    f_fb = w_fb / (2 * np.pi)
    if damping.type == 1:
        L.write(" Modal damping of all fixed-base modes, strain-energy weighted ratios (SN-5):")
        L.write("   MODE   FREQ(Hz)   DAMPING")
        for k in range(min(8, len(f_fb))):
            L.write(f" {k + 1:6d} {f_fb[k]:10.4f} {damping.modal_xi[k]:9.4f}")
    elif damping.type == 3:
        a, b = damping.alpha, damping.beta
        L.write(f" C = alpha M + beta K0 with alpha = {a:.6g} 1/s, beta = {b:.6g} s")
        for k in range(min(5, len(w_fb))):
            xi = a / (2 * w_fb[k]) + b * w_fb[k] / 2
            L.write(f"   mode {k + 1}: f = {f_fb[k]:.4f} Hz, damping ratio {xi:.4f}")
    else:
        L.write(" Element dashpots eta / h (the manual's eta dgamma/dt term, SN-5); damping ratio of mode 1:")
        for j, mdl in enumerate(models):
            xi1 = mdl.eta * w_fb[0] / (2.0 * col.rho[j] * col.vs[j] ** 2)
            L.write(f"   sublayer {j + 1:3d}: eta = {mdl.eta:.5g}  ->  xi(f1) = {xi1:.4f} (XI_MIN {mdl.xi_min:.4f})")


def _soil_non_outputs(ctx, d, L, res: SN.SoilNonResult, col: SN.SoilNonColumn, models, labels, col_sh, thick,
                      weight, vs0, vp0, ds0, dp0, curves, units, bed: int, ratio: float, fmax: float) -> int:
    """Listing tables and output files of SOIL-NON (same file names as SOIL-EQL, SN-12)."""
    grav = float(d["grav"])
    dt = res.dt
    nlay = len(thick)
    nsoil = nlay - 1
    files: List[str] = []
    warned_outcrop = []

    def acc_hist(lay: int, oc: bool) -> np.ndarray:
        """Absolute acceleration (length units / s^2) at the top of ``lay``; outcrop only at the elastic base."""
        if oc and lay == nlay and bed == 1:
            return res.input_motion
        if oc and not warned_outcrop:
            L.warning("SOIL-NON: the outcrop motion exists only at the elastic base (BedInt 1), where it is the "
                      "input; within motions are written for the other outcrop requests (SN-12)")
            warned_outcrop.append(1)
        return res.acc_history(lay)

    # ---------------- per-sublayer results ---------------------------------------------------
    L.section("Maximum strains and equivalent-linear properties of the nonlinear response (SN-12)")
    hdr = ["NO", "LABEL", "DEPTH", "MAX.STRAIN%", "TIME", "EFF.STRAIN%", "G/GMAX", "DAMPING", "DYNP G/GMAX",
           "DYNP DAMP", "MAX.STRESS"]
    w = [4, 10, 9, 13, 8, 13, 8, 9, 12, 10, 12]
    L.write(_fmt_row(hdr, w, [""] * len(w)))
    geff = ratio * res.gmax_mid * 100.0
    gg = np.ones(nsoil)
    damp = np.zeros(nsoil)
    for j in range(nsoil):
        gg[j], damp[j] = SN.equivalent_properties(models[j], geff[j])
        hist = res.strain_history(j + 1)
        _, tg = _peak(hist, dt)
        smax, _ = _peak(res.stress_history(j + 1), dt)
        lab = labels[j]
        if lab and lab in curves:
            cv = curves[lab]
            dg, dd = f"{cv.g_over_gmax(geff[j]):.3f}", f"{cv.damping(geff[j]):.4f}"
        else:
            dg, dd = "-", "-"
        L.write(_fmt_row([j + 1, lab or "-", col_sh.depth_mid[j], res.gmax_mid[j] * 100.0, tg, geff[j], gg[j],
                          damp[j], dg, dd, smax], w,
                         ["d", "s", ".2f", ".5f", ".2f", ".5f", ".3f", ".4f", "s", "s", ".4g"]))
    L.write(f" MAX.STRAIN% at mid-height over all sub-steps; EFF.STRAIN% = {ratio:g} x MAX.STRAIN% (SOIL <ratio>);")
    L.write(" G/GMAX = secant of the backbone and DAMPING = XI_MIN + Masing damping at EFF.STRAIN% (written to FILE88);")
    L.write(" DYNP columns: the curves at the same strain, for comparison with SOIL-EQL; MAX.STRESS total shear")
    L.write(f" stress (hysteretic + viscous, from equilibrium) in {units['mod']}")
    L.write(f" largest element strain in the column {100 * float(np.max(res.gmax_elem)) if res.gmax_elem.size else 0.0:.5f} %")

    L.section("Peak acceleration profile (absolute, top of sublayers)")
    L.write(_fmt_row(["LAYER", "DEPTH", "PGA(g)", "TIME(s)"], [6, 10, 10, 9], [""] * 4))
    for lay in range(1, nlay + 1):
        pk, tp = _peak(res.acc_history(lay) / grav, dt)
        L.write(_fmt_row([lay, col_sh.depth_top[lay - 1], pk, tp], [6, 10, 10, 9], ["d", ".2f", ".5f", ".2f"]))
    pk_in, _ = _peak(res.input_motion / grav, dt)
    L.write(f" input ({'outcrop' if bed == 1 else 'base'}) motion after the cut-off: {pk_in:.5f} g")

    # ---------------- accelerations ------------------------------------------------------
    sacc = d.rows("sacc")
    if sacc:
        L.section("Maximum accelerations (top of sublayer)")
        L.write(_fmt_row(["LAYER", "DEPTH", "MOTION", "MAX.ACC(g)", "TIME(s)", "MEAN.SQ.FR", "FILE"],
                         [6, 10, 9, 12, 9, 11, 12], [""] * 7))
        for r in sacc:
            lay, opt, oc = int(r["layer"]), int(r["opt"]), bool(r["outcrop"])
            if opt <= 0:
                continue
            if not (1 <= lay <= nlay):
                L.warning(f"SACC layer {lay} is not a sublayer (1..{nlay}); skipped")
                continue
            a = acc_hist(lay, oc) / grav
            is_oc = oc and lay == nlay and bed == 1
            pk, tp = _peak(a, dt)
            spec = np.fft.rfft(a)
            S = np.abs(spec) ** 2
            fr = np.fft.rfftfreq(len(a), dt)
            ms = float(np.sum(fr * S) / np.sum(S)) if np.sum(S) > 0 else 0.0
            fname = ""
            if opt >= 2:
                fname = C.layer_th_name("ACC", lay)
                textfiles.write_history(ctx.path(fname), a, dt)
                files.append(fname)
            L.write(_fmt_row([lay, col_sh.depth_top[lay - 1], "OUTCROP" if is_oc else "WITHIN", pk, tp, ms, fname],
                             [6, 10, 9, 12, 9, 11, 12], ["d", ".2f", "s", ".5f", ".2f", ".2f", "s"]))

    # ---------------- stresses / strains --------------------------------------------------
    sstr = d.rows("sstr")
    if sstr:
        L.section("Stress and strain at mid-height of sublayers")
        for r in sstr:
            lay = int(r["layer"])
            if not (1 <= lay <= nsoil):
                L.warning(f"SSTR layer {lay} is not a soil sublayer (1..{nsoil}); skipped")
                continue
            if r["opt3"] or r["opt4"]:
                g = res.strain_history(lay) * 100.0
                pk, tp = _peak(g, dt)
                L.write(f" layer {lay:3d}: maximum strain = {pk:.5f} % at t = {tp:.2f} s")
                if r["opt4"]:
                    fname = C.layer_th_name("SN", lay)
                    textfiles.write_history(ctx.path(fname), g, dt)
                    files.append(fname)
            if r["opt1"] or r["opt2"]:
                tau = res.stress_history(lay)
                pk, tp = _peak(tau, dt)
                L.write(f" layer {lay:3d}: maximum stress = {pk:.6g} {units['mod']} at t = {tp:.2f} s "
                        f"(total: hysteretic + viscous; hysteretic alone {float(np.max(np.abs(res.tau_h_mid[:, lay - 1]))):.6g})")
                if r["opt2"]:
                    fname = C.layer_th_name("SS", lay)
                    textfiles.write_history(ctx.path(fname), tau, dt)
                    files.append(fname)

    # ---------------- response spectra --------------------------------------------------
    damps = [float(r["value"]) for r in d.rows("damp")]
    srs = d.rows("srs")
    if srs and not damps:
        raise ModuleError("Error 107: no damping ratios (DAMP) for the response spectra")
    if srs and damps:
        L.section("Response spectra at the top of sublayers (absolute acceleration, g x gravmult)")
        L.write(f" {len(RS_FREQS)} frequencies, 0.1 - 100 Hz (100 per decade); dampings {damps}")
        for r in srs:
            lay, oc = int(r["layer"]), bool(r["outcrop"])
            if not (1 <= lay <= nlay):
                L.warning(f"SRS layer {lay} is not a sublayer; skipped")
                continue
            rs = SP.response_spectrum(acc_hist(lay, oc) / grav, dt, RS_FREQS, damps)
            is_oc = oc and lay == nlay and bed == 1
            for j, z in enumerate(damps):
                sa = rs["SA"][j] * float(d["gravmult"])
                k = int(np.argmax(sa))
                fname = ""
                if int(r["save"]):
                    fname = soil_rs_name(lay, j + 1)
                    textfiles.write_xy(ctx.path(fname), RS_FREQS, sa,
                                       header=f"SOIL-NON RS layer {lay} {'outcrop' if is_oc else 'within'} damping "
                                              f"{z:g} (f Hz, SA g x {d['gravmult']:g})")
                    files.append(fname)
                L.write(f" layer {lay:3d} {'OUTCROP' if is_oc else 'WITHIN '} damping {z:5.3f}: "
                        f"max SA = {sa[k]:.4f} at {RS_FREQS[k]:.3f} Hz (T = {1 / RS_FREQS[k]:.3f} s); "
                        f"SA(100 Hz) = {sa[-1]:.4f}  {fname}")

    # ---------------- spectral amplification factors ----------------------------------
    saf = d.rows("ssaf")
    if saf:
        L.section("Spectral amplification factors (SSAF) of the nonlinear motions")
        for r in saf:
            lay, lay2 = int(r["layer"]), int(r["layer2"])
            oc1, oc2 = bool(r["outcrop1"]), bool(r["outcrop2"])
            fst = float(r["freqstep"])
            if not (1 <= lay <= nlay) or not (1 <= lay2 <= nlay):
                L.warning(f"Error 109: illegal SSAF layers {lay}/{lay2}; skipped")
                continue
            if fst <= 0:
                L.warning(f"Error 110: SSAF frequency step {fst:g} <= 0; skipped")
                continue
            a1 = acc_hist(lay, oc1)
            a2 = acc_hist(lay2, oc2)
            npts = min(int(math.floor(fmax / fst + 1e-9)) + 1, MAX_SAF_POINTS)
            f = np.arange(npts) * fst
            fg = np.fft.rfftfreq(len(a1), dt)
            X1, X2 = np.fft.rfft(a1), np.fft.rfft(a2)
            # ratio of the Fourier transforms where the denominator carries energy (>= 1e-3 of its peak)
            ok = np.abs(X2) >= 1e-3 * float(np.max(np.abs(X2))) if X2.size and np.max(np.abs(X2)) > 0 \
                else np.zeros(X2.shape, bool)
            ratio_fft = np.zeros(X2.shape, complex)
            ratio_fft[ok] = X1[ok] / X2[ok]
            H = np.interp(f, fg, ratio_fft.real) + 1j * np.interp(f, fg, ratio_fft.imag)
            title = r["title"] or ""
            L.write(f" {title}")
            k = int(np.argmax(np.abs(H)))
            L.write(f"  Fourier ratio X({lay})/X({lay2}) of the computed motions (no transfer function exists for a "
                    f"nonlinear analysis; 0 where |X({lay2})| < 1e-3 of its peak): max {abs(H[k]):.4f} at "
                    f"{f[k]:.4f} Hz")
            if int(r["save"]):
                fname = saf_tf_name(lay, oc1, lay2, oc2)
                textfiles.write_tf(ctx.path(fname), f, H, complex_=True,
                                   header=f"SOIL-NON SSAF Fourier ratio of the computed motions {title}")
                files.append(fname)
            if damps:
                fr = f[1:][:MAX_SAF_RS_POINTS]
                if fr.size:
                    r1 = SP.response_spectrum(a1 / grav, dt, fr, damps)
                    r2 = SP.response_spectrum(a2 / grav, dt, fr, damps)
                    for j, z in enumerate(damps):
                        rat = r1["SA"][j] / np.where(r2["SA"][j] > 0, r2["SA"][j], np.inf)
                        kk = int(np.argmax(rat))
                        L.write(f"  RS ratio damping {z:5.3f}, f = k x {fst:g} Hz ({len(fr)} values): max "
                                f"{rat[kk]:.4f} at {fr[kk]:.3f} Hz")
                        if int(r["save"]):
                            fname = saf_rs_name(lay, oc1, lay2, oc2, j + 1)
                            textfiles.write_xy(ctx.path(fname), fr, rat,
                                               header=f"SOIL-NON SSAF RS ratio {title} damping {z:g} "
                                                      f"(f = k x {fst:g} Hz, RS({lay})/RS({lay2}))")
                            files.append(fname)

    # ---------------- FILE73, FILE88 -----------------------------------------------------
    if curves:
        SH.write_file73(ctx.path("FILE73"), curves, title=d["title"])
        files.append("FILE73")
    if int(d["save"]):
        rho = weight / grav
        G = rho[:-1] * vs0[:-1] ** 2 * gg
        vs = np.sqrt(G / rho[:-1])
        nonlin = np.array([not m.linear for m in models])
        if str(d["vppolicy"]).upper() == "VP":
            nonlin = np.zeros_like(nonlin)
        vp = np.where(nonlin, vp0[:-1] * vs / vs0[:-1], vp0[:-1])
        bp = np.where(nonlin, damp, dp0[:-1])
        hs_row = [thick[-1], weight[-1], vp0[-1], vs0[-1], dp0[-1], ds0[-1]]
        SH.write_file88(ctx.path("FILE88"), thick[:-1], geff, G, vs, damp, vp, bp,
                        units=f"{units['len']} {units['mod']} {units['vel']}; damping ratios; strain %",
                        halfspace=hs_row, title=(d["title"] + " " if d["title"] else "") + "(SOIL-NON)")
        files.append("FILE88")
        L.write("")
        L.write(f" FILE88 written: equivalent-linear properties of the nonlinear response of the {nsoil} soil sublayers "
                "(secant G and XI_MIN + Masing damping at the effective strain; Vp of nonlinear sublayers at "
                "constant Poisson's ratio and beta_p = beta_s, D-SOL-06)")
    L.section("Files written")
    for fn in files:
        L.write(f"  {fn}")
    ctx.progress(1.0, "SOIL-NON: done")
    return 0


if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
