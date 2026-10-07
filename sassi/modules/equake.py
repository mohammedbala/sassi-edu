"""EQUAKE module: spectrum-compatible acceleration histories (requirements section 4.12).

Decisions D-EQK-01 ... D-EQK-08, D-GEN-08 (PCG64 random numbers), D-CNV-08 (units);
numerics in :mod:`sassi.core.equake_lib`.

Input: ``<model>.equ`` (AFWRITE; schema ``sassi.io.decks.EQUAKE``).  For each spectrum
number 1..3 of the ``spectra`` table: target spectrum ``rsin`` (Hz, SA in g, ``nrfreq``
records, damping ``damp``), output spectrum ``rsout``, optional input history ``accin``
(seed record for ``accopt`` = 1, external history for ``accopt`` = 2), output history
``accout`` and optional target PSD ``tpsd_file`` (used when ``tpsd`` = 1; cm^2/s^3 or
in^2/s^3).

``accopt`` (D-EQK-01): 0 random phases, 1 seed-record phases, 2 external history (RS, PSD,
FFT and checks only).  ``seeds`` random-seed trials (0 -> 1) use the seeds rand, rand+1, ...;
the trial with the fewest failed criteria and then the smallest maximum deviation from the
target is kept (D-EQK-05).  Component ``k`` uses the independent PCG64 stream
``SeedSequence(rand + trial, spawn_key=(k,))``.  ``corr`` = 1 (P1) mixes spectrum 2 (Y) with
spectrum 1 (X) to the CORR correlation pairs in 2-s windows and re-matches it as a seed record,
with its own target PSD (D-EQK-06).

Outputs (names from the deck; ``.vel``, ``.dis``, ``.psd``, ``.fft`` share the ACCOUT stem):

* ``.acc`` acceleration in g, ``.vel`` / ``.dis`` in in/s, in (British) or cm/s, cm (SI)
  (D-EQK-07); first line dt, then one value per line (D-EQK-08, re-usable as THFILE);
* ``.rso`` response spectrum (Hz, SA g) at ``damp``, 100 points per decade over the target range;
* ``.psd`` one-sided PSD of the 5-75 % Arias window, +/-20 % band averaged (in^2/s^3 or cm^2/s^3);
* ``.fft`` (Hz, Re, Im) of the strong-motion window, f >= 0 (in/s or cm/s);
* listing with PGA/PGV/PGD, V/A, AD/V^2, strong-motion duration, long-period (drift) content,
  stationary and 2-s moving window correlations, and the manual and SRP 3.7.1 Rev. 4 criteria
  with pass/fail (D-EQK-04); every failed criterion is also issued as a warning (listing and
  screen, requirements 4.12 item 6) and counted in the run summary.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np

from .. import conventions as C
from ..core import equake_lib as EL
from ..core import spectra as SP
from ..io import decks, textfiles, thfile
from ..io import library as LIB
from .base import ModuleContext, ModuleError, batch_main

NAME = "EQUAKE"
MAX_STEPS = 32768
RHO_MAX = 0.16


@dataclass
class Component:
    """One generated (or external) acceleration history and its evaluation."""

    no: int
    acc: np.ndarray                       # g
    dt: float
    target: Optional[EL.TargetSpectrum]
    checks: Dict[str, EL.RSCheck] = field(default_factory=dict)
    psd: Optional[Dict[str, float]] = None
    params: Dict[str, float] = field(default_factory=dict)
    log: List[str] = field(default_factory=list)
    trial: int = 0
    files: List[str] = field(default_factory=list)


def _resolve(ctx: ModuleContext, name: str) -> Path:
    """Input file of the deck: ``@name`` from the built-in library (D-W5-01), else absolute or in the model
    directory."""
    return LIB.module_path(name, ctx.workdir)


def _resolve_out(ctx: ModuleContext, name: str) -> Path:
    """Output file of the deck (the built-in library is read-only: Errors 86 / 87)."""
    if LIB.is_library_name(name):
        raise ModuleError(f"Error 86/87: output file {str(name).strip()} is a built-in library file (read-only)")
    p = Path(str(name).strip())
    return p if p.is_absolute() else ctx.workdir / p


def _with_ext(path: Path, ext: str) -> Path:
    return path if path.suffix else path.with_suffix(ext)


def _units(gravity: float):
    """Physical acceleration unit (D-EQK-07): in/s^2 for British (ft/s^2 gravity), cm/s^2 for SI."""
    if C.unit_system(gravity) == "BS":
        return gravity * 12.0, "in", "BS"
    return gravity * 100.0, "cm", "SI"


def _psd_function(path: Path) -> Callable[[np.ndarray], np.ndarray]:
    f, p = EL.read_two_columns(path)
    keep = (f > 0) & (p > 0)
    f, p = f[keep], p[keep]
    if len(f) < 2:
        raise ModuleError(f"target PSD file {path.name}: needs >= 2 positive points")
    order = np.argsort(f)
    lf, lp = np.log10(f[order]), np.log10(p[order])

    def fn(x):
        return 10.0 ** np.interp(np.log10(np.maximum(np.asarray(x, float), 1e-12)), lf, lp)
    return fn


def _read_input_history(path: Path, dt: float, L) -> np.ndarray:
    """ACCIN history (fopt 0: dt on the first line, then values in g, D-EQK-08)."""
    acc, dt_file = thfile.read_history(path, fopt=0)
    if dt_file is None or dt_file <= 0:
        raise ModuleError(f"{path.name}: the first value must be the time step (THFILE fopt 0 layout)")
    if abs(dt_file - dt) > 1e-9 * dt:
        L.warning(f"{path.name}: time step {dt_file:g} s differs from EQUAKE time step {dt:g} s; "
                  f"the record is linearly re-sampled")
        t_old = np.arange(len(acc)) * dt_file
        t_new = np.arange(int(math.floor(t_old[-1] / dt + 1e-9)) + 1) * dt
        acc = np.interp(t_new, t_old, acc)
    return acc


def _validate(d, rows, L) -> None:
    if d["gravity"] <= 0:
        raise ModuleError("Error 1: gravity acceleration must be > 0")
    if not rows:
        raise ModuleError("Error 84: all spectrum input files (RSIN) are blank")
    if d["rand"] <= 0:
        raise ModuleError("Error 90: initial random number must be > 0")
    if d["nrfreq"] <= 0 and d["accopt"] != 2:
        raise ModuleError("Error 91: number of frequencies must be > 0")
    if d["dur"] <= 0:
        raise ModuleError("Error 92: total duration must be > 0")
    if d["delt"] <= 0:
        raise ModuleError("time step must be > 0")
    if d["corr"]:
        pairs = d.rows("corr")
        if not pairs:
            raise ModuleError("Error 93: correlated option selected but no correlation factors (CORR)")
        if any(abs(p["val"]) > 1.0 for p in pairs):
            raise ModuleError("Error 94: a correlation factor exceeds 1")
    if not 0.0 < d["damp"] < 0.5:
        raise ModuleError("damping of the target spectrum must be in (0, 0.5)")


# ---------------------------------------------------------------------------------------
def run(ctx: ModuleContext) -> int:
    d = decks.read(ctx.deck_path, NAME)
    L = ctx.listing
    accopt = int(d["accopt"])
    all_rows = sorted(d.rows("spectra"), key=lambda r: int(r["no"]))
    rows = [r for r in all_rows if str(r["rsin"]).strip() or (accopt == 2 and str(r["accin"]).strip())]
    _validate(d, rows, L)
    dt = float(d["delt"])
    dur = float(d["dur"])
    zeta = float(d["damp"])
    g_len, ulen, usys = _units(float(d["gravity"]))
    n = int(round(dur / dt)) + 1
    if n > MAX_STEPS:
        raise ModuleError(f"{n} time steps exceed the EQUAKE limit of {MAX_STEPS}")
    seeds = max(1, int(d["seeds"]))

    L.section("EQUAKE  spectrum-compatible acceleration time histories")
    if d["eqtit"]:
        L.write(f" Spectra title    : {d['eqtit']}")
    L.write(f" Model            : {d['model'] or ctx.model}   {d['title']}")
    L.write(f" Option <accopt>  : {accopt} ({['random phases', 'seed-record phases', 'external history: RS/PSD/FFT only'][min(accopt, 2)]})")
    L.write(f" Time step        : {dt:g} s;  total duration {dur:g} s ({n} values); Nyquist {0.5 / dt:g} Hz")
    L.write(f" Target damping   : {zeta:g};  random seed {int(d['rand'])};  seed trials {seeds}")
    L.write(f" Units            : gravity {d['gravity']:g} ({usys}); acceleration in g; velocity {ulen}/s; "
            f"displacement {ulen}; PSD {ulen}^2/s^3")
    if dur < 20.0:
        L.warning(f"total duration {dur:g} s < 20 s (manual and SRP 3.7.1)")
    if dt > 0.005:
        L.warning(f"time step {dt:g} s > 0.005 s (manual criterion; Nyquist {0.5 / dt:g} Hz < 100 Hz)")
    if dt > 0.010:
        L.warning(f"Nyquist frequency {0.5 / dt:g} Hz < 50 Hz (SRP 3.7.1 Rev. 4 Approach 2 (a))")
    if abs(zeta - 0.05) > 1e-12:
        L.write(" Note: the SRP 3.7.1 matching criteria are defined for 5 % damping; they are checked here "
                f"at the target damping {zeta:g}")
    opts = EL.MatchOptions()
    L.write(f" Matching (LW + AB): aim {opts.aim:g} x target, {opts.lw_iterations} LW iterations with crest-factor "
            f"clipping, AA2010 wavelets up to the zero-period plateau (at most {opts.wavelet_fmax * 0.5 / dt:g} Hz; "
            f"oscillators above {0.1 / dt:g} Hz on the 4x up-sampled record as in the checker), LW polish above")
    L.write(f" Drift control (D-EQK-03, SRP 3.7.1 'no baseline drift') in every step: zero-phase high-pass at "
            f"{opts.highpass:g} x the lowest check frequency (>= 1/duration), constrained least-squares fit of a "
            f"degree-{opts.baseline_degree} displacement polynomial (zero final velocity and displacement)")

    if int(d["opmode"]) == 1:
        for r in rows:
            for key, err in (("rsin", 85), ("accin", 88), ("tpsd_file", 0)):
                name = str(r[key]).strip()
                if name and not _resolve(ctx, name).exists():
                    raise ModuleError(f"{'Error %d: ' % err if err else ''}file {name} (spectrum {r['no']}) "
                                      f"does not exist")
            if str(r["rsin"]).strip():
                f, _ = EL.read_two_columns(_resolve(ctx, r["rsin"]))
                if int(d["nrfreq"]) > 0 and len(f) != int(d["nrfreq"]):
                    raise ModuleError(f"Error 89: {r['rsin']} has {len(f)} records, Number of Frequencies = "
                                      f"{d['nrfreq']}")
        L.write("")
        L.write(" Operation mode 1: data check only -- input files checked, no motion generated")
        return 0

    comps: List[Component] = []
    for k, r in enumerate(rows):
        no = int(r["no"])
        ctx.progress(0.05 + 0.85 * k / max(len(rows), 1), f"EQUAKE: spectrum {no}")
        ctx.announce("EQUAKE.fit", no=no, k=k + 1, n=len(rows), dur=float(dur), dt=float(dt), zeta=float(zeta),
                     accopt=int(accopt))
        comp = _one_spectrum(ctx, d, r, accopt, dt, dur, n, zeta, g_len, ulen, seeds, opts, comps)
        comps.append(comp)

    # correlated X-Y components (P1)
    if d["corr"] and accopt != 2 and len(comps) >= 2:
        _correlate(ctx, d, rows, comps, dt, dur, zeta, g_len, ulen, opts)

    _correlations(L, comps, dt)
    _srp_summary(L, comps, dt, dur, ulen)
    ctx.progress(1.0, "EQUAKE: done")
    return 0


def _one_spectrum(ctx, d, r, accopt, dt, dur, n, zeta, g_len, ulen, seeds, opts, previous) -> Component:
    L = ctx.listing
    no = int(r["no"])
    L.section(f"Spectrum number {no}")
    target = None
    if str(r["rsin"]).strip():
        p = _resolve(ctx, r["rsin"])
        if not p.exists():
            raise ModuleError(f"Error 85: spectrum input file {r['rsin']} does not exist")
        f, sa = EL.read_two_columns(p)
        if int(d["nrfreq"]) > 0 and len(f) != int(d["nrfreq"]):
            raise ModuleError(f"Error 89: {r['rsin']} has {len(f)} records, Number of Frequencies = {d['nrfreq']}")
        try:
            target = EL.TargetSpectrum(f, sa, zeta)
        except ValueError as exc:
            raise ModuleError(f"{r['rsin']}: {exc}")
        L.write(f" Target spectrum  : {r['rsin']} ({len(f)} points, {target.fmin:g} - {target.fmax:g} Hz, "
                f"damping {zeta:g}); log-log interpolation")
        if LIB.is_library_name(r["rsin"]):
            L.write(f"   {LIB.note(r['rsin'])}")
        f_srp = min(50.0, 0.5 / dt)
        if target.fmax < f_srp - 1e-9:
            L.write(f"   above {target.fmax:g} Hz the target is extended at constant SA = {target.sa[-1]:g} g "
                    f"(zero-period plateau) up to {f_srp:g} Hz for the SRP 3.7.1 band")
        if target.fmin > 0.1 + 1e-9:
            L.write(f"   the target starts above 0.1 Hz: the SRP 3.7.1 band (b) 0.1 - {f_srp:g} Hz cannot be "
                    f"checked completely (it is not extrapolated below {target.fmin:g} Hz)")
        L.write("        FREQ (Hz)     SA (g)")
        for fi, si in zip(target.freq, target.sa):
            L.write(f"     {fi:12.4f} {si:10.4f}")
    tpsd_fn = None
    if int(d["tpsd"]) and str(r["tpsd_file"]).strip():
        pp = _resolve(ctx, r["tpsd_file"])
        if not pp.exists():
            raise ModuleError(f"target PSD file {r['tpsd_file']} does not exist")
        tpsd_fn = _psd_function(pp)
        L.write(f" Target PSD       : {r['tpsd_file']} ({ulen}^2/s^3); check >= 80 % over 0.3 - 24 Hz (SRP App. A)")
        if LIB.is_library_name(r["tpsd_file"]):
            L.write(f"   {LIB.note(r['tpsd_file'])}")
    elif int(d["tpsd"]):
        L.write(" Target PSD       : none for this spectrum (PSD criterion not checked)")

    if accopt == 2:
        if not str(r["accin"]).strip():
            raise ModuleError(f"Error 88: acceleration input file {no} is blank (External Accel.)")
        pin = _resolve(ctx, r["accin"])
        if not pin.exists():
            raise ModuleError(f"acceleration input file {r['accin']} does not exist")
        acc, dt_file = thfile.read_history(pin, fopt=0)
        if dt_file is None or dt_file <= 0:
            raise ModuleError(f"{r['accin']}: the first value must be the time step")
        if abs(dt_file - dt) > 1e-9 * dt or len(acc) != n:
            L.warning(f"{r['accin']}: dt = {dt_file:g} s and {len(acc)} values (EQUAKE: {dt:g} s, {n} values); "
                      "the file's own time step is used")
        comp = Component(no, acc, dt_file, target)
        L.write(f" External history : {r['accin']} ({len(acc)} values, dt {dt_file:g} s); no simulation")
        if LIB.is_library_name(r["accin"]):
            L.write(f"   {LIB.note(r['accin'])}")
    else:
        if not str(r["accout"]).strip():
            raise ModuleError(f"Error 87: acceleration output file {no} is blank")
        if target is None:
            raise ModuleError(f"Error 84/85: spectrum {no} has no target spectrum")
        seed_acc = None
        if accopt == 1:
            if not str(r["accin"]).strip():
                raise ModuleError(f"Error 88: acceleration input file {no} is blank (Accel. Record)")
            pin = _resolve(ctx, r["accin"])
            if not pin.exists():
                raise ModuleError(f"seed record {r['accin']} does not exist")
            seed_acc = _read_input_history(pin, dt, L)
            L.write(f" Seed record      : {r['accin']} ({len(seed_acc)} values); its Fourier phases are kept")
            if LIB.is_library_name(r["accin"]):
                L.write(f"   {LIB.note(r['accin'])}")
        ntrial = 1 if seed_acc is not None else seeds
        if seed_acc is not None and seeds > 1:
            L.write(" (seed-record phases are deterministic: one trial)")
        best: Optional[Component] = None
        best_key = None
        for tr in range(ntrial):
            rng = None
            if seed_acc is None:
                rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(int(d["rand"]) + tr,
                                                                                 spawn_key=(no,))))
            res = EL.match_spectrum(target, dt, dur, zeta, rng=rng, seed_acc=seed_acc, target_psd=tpsd_fn,
                                    g_len=g_len, options=opts)
            cand = Component(no, res.acc, dt, target, res.checks, res.psd, log=res.log, trial=tr)
            fails, dev = res.score()
            rho = max([abs(EL.correlation(cand.acc, c.acc)) for c in previous] + [0.0])
            if rho > RHO_MAX:
                fails += 1
            key = (fails, dev)
            L.write(f" trial {tr + 1} (seed {int(d['rand']) + tr}): failed criteria {fails}, max deviation "
                    f"{dev:.3f}, max |rho| with previous components {rho:.3f}")
            if best_key is None or key < best_key:
                best, best_key = cand, key
        comp = best
        L.write(f" Selected trial {comp.trial + 1}.  Generation log:")
        for ln in comp.log:
            L.write(f"   {ln}")
    _evaluate_and_write(ctx, d, r, comp, zeta, g_len, ulen, tpsd_fn, accopt)
    return comp


def _evaluate_and_write(ctx, d, r, comp: Component, zeta, g_len, ulen, tpsd_fn, accopt) -> None:
    L = ctx.listing
    dt = comp.dt
    a = comp.acc
    if comp.target is not None:
        comp.checks = EL.evaluate_rs(a, dt, comp.target, zeta)
    if tpsd_fn is not None:
        comp.psd = EL.psd_check(a * g_len, dt, tpsd_fn)
    comp.params = EL.ground_motion_parameters(a, dt, g_len)
    comp.params["frac_long"] = EL.drift_measures(a, dt, g_len)["frac_long"]
    pm = comp.params
    L.write("")
    L.write(" Ground-motion parameters")
    L.write(f"   PGA = {pm['PGA']:.4f} g   PGV = {pm['PGV']:.4f} {ulen}/s   PGD = {pm['PGD']:.4f} {ulen}")
    L.write(f"   V/A = {pm['V/A']:.4f} s   AD/V^2 = {pm['AD/V2']:.3f}   Arias intensity = {pm['Arias']:.4g} {ulen}/s")
    L.write(f"   strong-motion duration (Arias 5-75 %) = {pm['D5-75']:.2f} s ({pm['t5']:.2f} - {pm['t75']:.2f} s)")
    L.write(f"   final velocity {pm['vel_end']:.3e} {ulen}/s, final displacement {pm['dis_end']:.3e} {ulen}")
    L.write(f"   long-period (drift) content: {100 * pm['frac_long']:.1f} % of the displacement Fourier energy at "
            f"periods longer than the record ({(len(a) - 1) * dt:g} s)")
    for key, c in comp.checks.items():
        L.write(f" Spectral match, {c.name}: {c.fmin:g} - {c.fmax:g} Hz, {c.npts} points "
                f"({c.per_decade:.1f} per decade)")
        L.write(f"   SA/target min {c.min_ratio:.3f} at {c.f_min_ratio:.3f} Hz, max {c.max_ratio:.3f} at "
                f"{c.f_max_ratio:.3f} Hz, mean {c.mean_ratio:.3f}")
        L.write(f"   points > 10 % below: {c.n_below_10} [{'PASS' if c.ok_below else 'FAIL'}]   "
                f"points > 30 % above: {c.n_above_30} [{'PASS' if c.ok_above else 'FAIL'}]   "
                f"max adjacent points below: {c.max_run_below} [{'PASS' if c.ok_run else 'FAIL'}]")
    if comp.psd is not None:
        L.write(f" PSD check (SRP App. A): min ratio {comp.psd['min_ratio']:.3f} at {comp.psd['f_min']:.2f} Hz "
                f"over 0.3 - 24 Hz [{'PASS' if comp.psd['passed'] else 'FAIL'}]")

    # ---- files
    if accopt != 2:
        pacc = _with_ext(_resolve_out(ctx, r["accout"]), ".acc")
        textfiles.write_history(pacc, a, dt)
        vel, dis = SP.integrate(a * g_len, dt)
        textfiles.write_history(pacc.with_suffix(".vel"), vel, dt)
        textfiles.write_history(pacc.with_suffix(".dis"), dis, dt)
        stem = pacc
        comp.files += [pacc.name, pacc.with_suffix(".vel").name, pacc.with_suffix(".dis").name]
    else:
        stem = _resolve_out(ctx, r["accout"]) if str(r["accout"]).strip() else \
            _resolve_out(ctx, LIB.bare(r["accin"]) if LIB.is_library_name(r["accin"]) else r["accin"])
    if str(r["rsout"]).strip() and comp.target is not None:
        band = EL.check_bands(comp.target, dt)["manual"]
        fc = EL.check_grid(*band)
    else:
        fc = EL.check_grid(0.1, min(100.0, 0.5 / dt))
    if str(r["rsout"]).strip():
        prs = _with_ext(_resolve_out(ctx, r["rsout"]), ".rso")
        sa = SP.response_spectrum(a, dt, fc, [zeta])["SA"][0]
        textfiles.write_xy(prs, fc, sa, header=f"EQUAKE response spectrum, damping {zeta:g} (f Hz, SA g)")
        comp.files.append(prs.name)
    elif accopt != 2:
        L.warning(f"Error 86 condition: spectrum output file {comp.no} is blank; no .rso written")
    acc_phys = a * g_len
    fp, psd = EL.strong_motion_psd(acc_phys, dt)
    ppsd = stem.with_suffix(".psd")
    textfiles.write_xy(ppsd, fp, psd, header=f"EQUAKE PSD of the 5-75 % Arias window, +/-20 % averaged "
                                             f"(f Hz, {ulen}^2/s^3)")
    t5, t75 = comp.params["t5"], comp.params["t75"]
    i0, i1 = int(round(t5 / dt)), int(round(t75 / dt)) + 1
    seg = acc_phys[i0:i1]
    F = dt * np.fft.rfft(seg)
    pfft = stem.with_suffix(".fft")
    textfiles.write_xy(pfft, np.fft.rfftfreq(len(seg), dt), F.real, F.imag,
                       header=f"EQUAKE Fourier transform of the strong-motion window (f Hz, Re, Im in {ulen}/s)")
    comp.files += [ppsd.name, pfft.name]
    L.write(" Files: " + ", ".join(comp.files))


def _correlate(ctx, d, rows, comps, dt, dur, zeta, g_len, ulen, opts) -> None:
    """D-EQK-06 (P1): impose rho(t) between X (spectrum 1) and Y (spectrum 2), then re-match Y as a seed record.

    The CORR pairs define the X-Y correlation, so the components are chosen by their spectrum
    numbers (1 = X, 2 = Y), not by their position in the table.  Y is re-matched with its own
    target PSD (when ``tpsd`` = 1), so the PSD floor applies to the re-matched record as well.
    """
    L = ctx.listing
    pairs = sorted(d.rows("corr"), key=lambda p: float(p["time"]))
    times = [float(p["time"]) for p in pairs]
    vals = [float(p["val"]) for p in pairs]
    by_no = {c.no: c for c in comps}
    L.section("Correlated components (P1, D-EQK-06)")
    if 1 not in by_no or 2 not in by_no:
        L.warning("correlated option: spectra 1 (X) and 2 (Y) are both required; correlation not imposed "
                  f"(spectra given: {sorted(by_no)})")
        return
    x, y = by_no[1], by_no[2]
    L.write(" target correlation pairs (time s, rho): " + ", ".join(f"({t:g}, {v:g})" for t, v in zip(times, vals)))
    r = next(rr for rr in rows if int(rr["no"]) == y.no)
    tpsd_fn = None
    if int(d["tpsd"]) and str(r["tpsd_file"]).strip():
        tpsd_fn = _psd_function(_resolve(ctx, r["tpsd_file"]))
    mixed = EL.correlate_pair(x.acc, y.acc, dt, times, vals)
    res = EL.match_spectrum(y.target, dt, dur, zeta, seed_acc=mixed, target_psd=tpsd_fn, g_len=g_len, options=opts)
    y.acc, y.checks, y.psd, y.log = res.acc, res.checks, res.psd, res.log
    y.files = []
    L.write(f" spectrum {y.no} (Y) re-generated with the mixed record as seed (phases kept"
            f"{', target PSD floor applied' if tpsd_fn is not None else ''}); the correlation with spectrum "
            f"{x.no} (X) after re-matching is reported below (manual WARNING: re-matching changes it)")
    for ln in y.log:
        L.write(f"   {ln}")
    _evaluate_and_write(ctx, d, r, y, zeta, g_len, ulen, tpsd_fn, 1)


def _correlations(L, comps: List[Component], dt: float) -> None:
    if len(comps) < 2:
        return
    L.section("Correlation between components")
    L.write(" stationary correlation coefficients (entire duration)")
    for i in range(len(comps)):
        for j in range(i + 1, len(comps)):
            r = EL.correlation(comps[i].acc, comps[j].acc)
            L.write(f"   spectrum {comps[i].no} - {comps[j].no}: rho = {r:+.4f} "
                    f"[{'PASS' if abs(r) <= RHO_MAX else 'FAIL'} |rho| <= {RHO_MAX}]")
    L.write(" non-stationary correlation (2-s moving window, 0.5-s step)")
    for i in range(len(comps)):
        for j in range(i + 1, len(comps)):
            t, r = EL.moving_correlation(comps[i].acc, comps[j].acc, dt)
            if len(r):
                k = int(np.argmax(np.abs(r)))
                L.write(f"   spectrum {comps[i].no} - {comps[j].no}: max |rho| = {abs(r[k]):.3f} at t = {t[k]:.1f} s; "
                        f"mean rho = {np.mean(r):+.3f}")
                L.write("     " + " ".join(f"{v:+.2f}" for v in r[::4]))


def _srp_summary(L, comps: List[Component], dt_deck: float, dur_deck: float, ulen: str) -> None:
    """Acceptance criteria tables (D-EQK-04) with a listing (and screen) warning for every failed one.

    Requirements 4.12 item 6: "checks: warnings in output and screen".  Every failed criterion of
    every component, and the inter-component independence check, issues ``L.warning`` in addition
    to its FAIL row; the total time step / duration of the deck were already warned at the start
    of the run and are warned again only for an external history whose own values differ.
    """
    L.section("Acceptance criteria summary (D-EQK-04: manual and SRP 3.7.1 Rev. 4 reported separately)")
    nfail = [0]

    def row(label: str, ok: bool, warn: Optional[str] = None) -> None:
        L.write(f"     {label:<66s}{' ' if len(label) >= 66 else ''}{'PASS' if ok else 'FAIL'}")
        if not ok:
            nfail[0] += 1
            if warn:
                L.warning(warn)

    rho = [(abs(EL.correlation(comps[i].acc, comps[j].acc)), comps[i].no, comps[j].no) for i in range(len(comps))
           for j in range(i + 1, len(comps))]
    for c in comps:
        dt = c.dt
        dur = (len(c.acc) - 1) * c.dt
        own = abs(dt - dt_deck) > 1e-12 or abs(dur - dur_deck) > 1e-9     # external history with its own dt/length
        sp = f"spectrum {c.no}"
        L.write(f" Spectrum {c.no}")
        man, srp = c.checks.get("manual"), c.checks.get("srp")
        L.write("   ACS SASSI manual criteria")
        row(f"total duration >= 20 s ({dur:g} s)", dur >= 20.0 - 1e-9,
            f"{sp}: total duration {dur:g} s < 20 s (manual)" if own else None)
        row(f"time step <= 0.005 s ({dt:g} s)", dt <= 0.005 + 1e-12,
            f"{sp}: time step {dt:g} s > 0.005 s (manual)" if own else None)
        if man is not None:
            row(f">= 100 points per decade ({man.per_decade:.1f}, {man.fmin:g} - {man.fmax:g} Hz)", man.ok_density,
                f"{sp}: fewer than 100 points per decade ({man.per_decade:.1f}, manual)")
            row(f"no point > 10 % below the target (min ratio {man.min_ratio:.3f})", man.ok_below,
                f"{sp}: {man.n_below_10} RS point(s) more than 10 % below the target (min SA/target "
                f"{man.min_ratio:.3f} at {man.f_min_ratio:.3f} Hz; manual)")
            row(f"no point > 30 % above the target (max ratio {man.max_ratio:.3f})", man.ok_above,
                f"{sp}: {man.n_above_30} RS point(s) more than 30 % above the target (max SA/target "
                f"{man.max_ratio:.3f} at {man.f_max_ratio:.3f} Hz; manual)")
            row(f"<= 9 adjacent points below the target ({man.max_run_below})", man.ok_run,
                f"{sp}: {man.max_run_below} adjacent RS points below the target (> 9; manual)")
        L.write("   SRP 3.7.1 Rev. 4 (Option 1, Approach 2, and general requirements)")
        row(f"(a) Nyquist >= 50 Hz ({0.5 / dt:g} Hz) and total duration >= 20 s", dt <= 0.010 + 1e-12 and dur >= 20.0 - 1e-9,
            f"{sp}: SRP 3.7.1 (a) not met (Nyquist {0.5 / dt:g} Hz, duration {dur:g} s)" if own else None)
        if srp is not None:
            lo, hi = srp.req_fmin, srp.req_fmax
            full = srp.band_complete
            row(f"(b) 100 points/decade, {srp.fmin:g} - {srp.fmax:g} Hz"
                f"{'' if full else f' (incomplete: {lo:g} - {hi:g} Hz required)'}",
                srp.ok_density and full,
                f"{sp}: SRP 3.7.1 (b) not met: "
                + (f"{srp.per_decade:.1f} points per decade" if not srp.ok_density else
                   f"checked band {srp.fmin:g} - {srp.fmax:g} Hz does not cover {lo:g} - {hi:g} Hz "
                   f"(target defined from {c.target.fmin if c.target is not None else srp.fmin:g} Hz)"))
            row(f"(c) none > 10 % below, <= 9 adjacent below (min {srp.min_ratio:.3f}, run {srp.max_run_below})",
                srp.ok_below and srp.ok_run,
                f"{sp}: SRP 3.7.1 (c) not met ({srp.n_below_10} point(s) > 10 % below, min SA/target "
                f"{srp.min_ratio:.3f} at {srp.f_min_ratio:.3f} Hz; {srp.max_run_below} adjacent points below)")
            row(f"(d) none > 30 % above (max {srp.max_ratio:.3f})", srp.ok_above,
                f"{sp}: SRP 3.7.1 (d) not met ({srp.n_above_30} point(s) > 30 % above, max SA/target "
                f"{srp.max_ratio:.3f} at {srp.f_max_ratio:.3f} Hz)")
        d575 = c.params.get("D5-75", 0.0)
        row(f"strong-motion duration (Arias 5-75 %) >= 6 s ({d575:.2f} s)", d575 >= 6.0,
            f"{sp}: strong-motion duration (Arias 5-75 %) {d575:.2f} s < 6 s (SRP 3.7.1)")
        if c.psd is not None:
            row(f"PSD >= 80 % of the target PSD, 0.3 - 24 Hz (min ratio {c.psd['min_ratio']:.3f})", c.psd["passed"],
                f"{sp}: PSD below 80 % of the target PSD (min ratio {c.psd['min_ratio']:.3f} at "
                f"{c.psd['f_min']:.2f} Hz; SRP 3.7.1 App. A)")
        else:
            L.write("     PSD: no target PSD given (not checked)")
        L.write(f"     V/A = {c.params.get('V/A', 0.0):.4f} s, AD/V^2 = {c.params.get('AD/V2', 0.0):.3f}, "
                f"PGD = {c.params.get('PGD', 0.0):.3f} {ulen} (to be compared with the controlling events); "
                f"long-period displacement energy {100 * c.params.get('frac_long', 0.0):.1f} %")
    if rho:
        L.write(" All components")
        mx = max(r[0] for r in rho)
        row(f"statistical independence |rho| <= {RHO_MAX} (max {mx:.4f})", mx <= RHO_MAX)
        for r, i, j in rho:
            if r > RHO_MAX:
                L.warning(f"spectra {i} and {j} are not statistically independent: |rho| = {r:.4f} > {RHO_MAX} "
                          f"(SRP 3.7.1)")
    L.write("")
    L.write(f" {nfail[0]} criterion check(s) failed" + (" (see the warnings above)" if nfail[0] else ""))


if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
