"""SITE module: layered free field by the thin-layer method (requirements §4.2, R1 §2).

Mode 1 (FILE2): for every analysis frequency, generate the variable-depth half-space sublayers
(1.5 lambda_s, law D-SIT-02) with base dashpots, assemble the thin-layer matrices (Kausel form,
mixed mass) and solve the Rayleigh and Love eigenproblems.  FILE2 serves SITE Mode 2 and POINT,
2D and 3D.

Mode 2 (FILE1): free-field motion at the user interfaces for unit control motion (within motion
at the top of control layer ``cl``, direction ``cm``, D-SIT-07) for each wave field of the WAVE
records, with the participation ratios interpolated between Frequency 1 and 2.

    python -m sassi.modules.site < SITE.inp        # batch use (three-line .inp protocol)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import numpy as np

from ..conventions import cfactor, frequency_step
from ..core import freefield as ff
from ..core import tlm
from ..io import decks
from ..io.files import read_container, validate, write_container
from .base import ModuleContext, ModuleError, batch_main

NAME = "SITE"
MASS = tlm.MASS_MIXED


# ---------------------------------------------------------------------------------------
# Input model
# ---------------------------------------------------------------------------------------
@dataclass
class Wave:
    type: int
    opt: int
    ratio1: float
    ratio2: float
    angle: float

    @property
    def name(self) -> str:
        return ff.WAVE_NAMES[self.type]


@dataclass
class SiteInput:
    """Resolved SITE input (user layers + half-space, frequencies, wave environment)."""

    title: str
    layer_no: List[int]
    thick: np.ndarray            # (nI-1,)
    weight: np.ndarray           # (nI,) user layers + half-space
    rho: np.ndarray
    vs: np.ndarray
    vp: np.ndarray
    ds: np.ndarray
    dp: np.ndarray
    G: np.ndarray                # complex (nI,)
    M: np.ndarray
    nl: int
    hslaw: str
    df: float
    fnum: np.ndarray
    mode1: bool
    mode2: bool
    cl: int
    cm: int
    wopt: int
    waves: List[Wave] = field(default_factory=list)
    freq1: int = 1
    freq2: int = 2048
    cmodform: int = 0
    gravity: float = 32.2
    hs_no: int = 0

    @property
    def nI(self) -> int:
        return len(self.thick) + 1

    @property
    def freq(self) -> np.ndarray:
        return self.fnum * self.df

    @property
    def depth_user(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.thick)])

    @property
    def base(self) -> str:
        return tlm.BASE_DASHPOT if self.nl > 0 else tlm.BASE_RIGID


def _read_file88(path: Path, n_layers: int, lst) -> Dict[str, np.ndarray]:
    """Strain-compatible layer properties written by SOIL (D-FIL-10, requirements §4.11)::

        # columns: layer thick gamma_eff_pct G Vs beta_s Vp beta_p
        1  5.0  7.7e-04  3851.5  996.2  0.0071  1863.8  0.0071

    Lines starting with '#' or '*' are comments (the ``# columns:`` line names the columns);
    without it the positional layout ``[layer] thick geff G Vs beta_s Vp beta_p`` is assumed.
    Rows are matched to the TOPL layers by position (D-SOL-12); a count mismatch is EDU-07.
    """
    names = None
    rows = []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s:
            continue
        if s[0] in "#*":
            low = s.lstrip("#*").strip().lower()
            if low.startswith("columns:"):
                names = low[len("columns:"):].split()
            continue
        try:
            rows.append([float(t.replace("D", "E").replace("d", "e")) for t in s.replace(",", " ").split()])
        except ValueError:
            raise ModuleError(f"FILE88: cannot read the line '{s}'") from None
    if not rows:
        raise ModuleError("FILE88 holds no numeric rows")
    width = min(len(r) for r in rows)
    data = np.asarray([r[:width] for r in rows])
    if names is not None and len(names) == width:
        def col(*keys):
            for kk in keys:
                if kk in names:
                    return data[:, names.index(kk)]
            raise ModuleError(f"FILE88 '# columns:' header lacks one of {keys}")
        out = {"vs": col("vs"), "ds": col("beta_s", "ds", "bs"), "vp": col("vp"), "dp": col("beta_p", "dp", "bp")}
    else:
        off = width - 7
        if off not in (0, 1):
            raise ModuleError("FILE88 rows must hold [layer] thick geff G Vs beta_s Vp beta_p")
        out = {"vs": data[:, off + 3], "ds": data[:, off + 4], "vp": data[:, off + 5], "dp": data[:, off + 6]}
    if len(out["vs"]) != n_layers:
        raise ModuleError(f"EDU-07: FILE88 has {len(out['vs'])} layers, TOPL has {n_layers}")
    lst.write(f" FILE88 strain-compatible properties read for {n_layers} layers (SITEX soil mode 1)")
    return out


def read_input(d: decks.Deck, workdir: Path, lst) -> SiteInput:
    """Validate the deck and resolve the layer properties (CHECK rules of 05a §7.3)."""
    g = float(d["gravity"])
    if g <= 0:
        raise ModuleError("Error 1: gravity must be positive")
    mode1, mode2 = bool(d["mode1"]), bool(d["mode2"])
    if not (mode1 or mode2):
        raise ModuleError("Error 45: SITE Mode 1 and Mode 2 are both off")
    rows = d.rows("layers")
    if not rows:
        raise ModuleError("Error 46: no top layers (TOPL)")
    hsr = d.rows("halfspace")
    if not hsr:
        raise ModuleError("SITE deck has no half-space row")
    allr = rows + [hsr[0]]
    thick = np.array([r["thick"] for r in rows], float)
    weight = np.array([r["weight"] for r in allr], float)
    vs = np.array([r["vs"] for r in allr], float)
    vp = np.array([r["vp"] for r in allr], float)
    ds = np.array([r["ds"] for r in allr], float)
    dp = np.array([r["dp"] for r in allr], float)
    if int(d["soilmode"]) == 1:
        p88 = Path(workdir) / "FILE88"
        if not p88.exists():
            raise ModuleError("FILE88 missing -- run SOIL first (SITEX soil mode 1)")
        sc = _read_file88(p88, len(rows), lst)
        vs[:-1], ds[:-1], vp[:-1], dp[:-1] = sc["vs"], sc["ds"], sc["vp"], sc["dp"]
    if np.any(thick <= 0):
        raise ModuleError("Error 20: layer thickness must be positive")
    if np.any(weight <= 0):
        raise ModuleError("Error 21: layer unit weight must be positive")
    if np.any(vs <= 0) or np.any(vp <= 0):
        raise ModuleError("Errors 22/23: layer velocities must be positive")
    if np.any(vp <= vs * np.sqrt(4.0 / 3.0)):
        raise ModuleError("layer with Vp <= sqrt(4/3) Vs (negative bulk modulus)")
    if np.any(ds < 0) or np.any(dp < 0):
        raise ModuleError("Error 25: negative layer damping")
    if np.any(ds >= 0.5) or np.any(dp >= 0.5):
        raise ModuleError("EDU-04: damping ratio >= 0.5")
    nl = int(d["nl"])
    if nl != 0 and not 4 <= nl <= 20:
        raise ModuleError(f"Error 47: number of generated half-space layers {nl} must be 0 or 4..20")
    if 0 < nl < 10:
        lst.warning(f"D-SIT-03: {nl} generated half-space layers (10-20 recommended)")
    hslaw = str(d["hslaw"] or tlm.HS_GEOMETRIC).lower()
    if hslaw not in tlm.HS_LAWS:
        raise ModuleError(f"unknown half-space sublayer law '{hslaw}' (geometric, uniform, linear)")
    if float(d["fstep"]) < 0:
        raise ModuleError("Error 48: frequency step < 0")
    if float(d["delt"]) < 0:
        raise ModuleError("Error 49: time step < 0")
    if int(d["nft"]) < 0:
        raise ModuleError("Error 50: number of Fourier components < 0")
    df = float(d["df"])
    if df <= 0:
        try:
            df = frequency_step(float(d["delt"]), int(d["nft"]), float(d["fstep"]))
        except ValueError as exc:
            raise ModuleError(f"cannot resolve the frequency step: {exc}") from None
    fn = [int(r["number"]) for r in d.rows("freqs")]
    if not fn:
        raise ModuleError("Error 120: frequency set is empty")
    if min(fn) <= 0:
        raise ModuleError("frequency numbers must be positive integers")
    if len(set(fn)) != len(fn):
        raise ModuleError("frequency set contains duplicate frequency numbers")
    fnum = np.array(sorted(fn), int)
    if fnum.size > 500:
        lst.warning(f"EDU-02: {fnum.size} analysis frequencies (manual limit 500)")
    cl, cm, wopt = int(d["cl"]), int(d["cm"]), int(d["wopt"])
    nI = len(rows) + 1
    if not 1 <= cl <= nI:
        raise ModuleError(f"control point layer {cl} outside 1..{nI}")
    if cm not in (0, 1, 2):
        raise ModuleError(f"control motion direction {cm} must be 0 (x'), 1 (y') or 2 (z')")
    if wopt not in (0, 1):
        raise ModuleError(f"wave family <wopt> = {wopt} must be 0 (R/SV/P) or 1 (SH/L)")
    rho = weight / g
    form = int(d["cmodform"])
    G = rho * vs ** 2 * np.asarray(cfactor(ds, form))
    M = rho * vp ** 2 * np.asarray(cfactor(dp, form))
    waves: List[Wave] = []
    for r in d.rows("waves"):
        w = Wave(int(r["type"]), int(r["opt"]), float(r["ratio1"]), float(r["ratio2"]), float(r["angle"]))
        if w.type not in ff.WAVE_NAMES:
            raise ModuleError(f"WAVE type {w.type} must be 1 R, 2 SV, 3 P, 4 SH or 5 L")
        if w.opt == 0:
            continue
        if ff.WAVE_FAMILY[w.type] != wopt:
            lst.warning(f"{w.name} wave field ignored: not in the selected wave family <wopt> = {wopt}")
            continue
        waves.append(w)
    inp = SiteInput(title=str(d["title"]), layer_no=[int(r["no"]) for r in rows], thick=thick, weight=weight,
                    rho=rho, vs=vs, vp=vp, ds=ds, dp=dp, G=np.asarray(G, complex), M=np.asarray(M, complex),
                    nl=nl, hslaw=hslaw, df=df, fnum=fnum, mode1=mode1, mode2=mode2, cl=cl, cm=cm, wopt=wopt,
                    waves=waves, freq1=int(d["freq1"]), freq2=int(d["freq2"]), cmodform=form, gravity=g,
                    hs_no=int(d["hs"]))
    if mode2:
        _check_waves(inp)
    _guideline_warnings(inp, lst)
    return inp


def _check_waves(inp: SiteInput) -> None:
    """Wave-environment checks (Errors 51-55, D-SIT-07..09)."""
    if not inp.waves:
        raise ModuleError("Error 51: all wave fields are off")
    if inp.freq1 <= 0 or inp.freq2 <= 0:
        raise ModuleError("Error 53: wave-ratio frequencies must be positive")
    for w in inp.waves:
        for rr in (w.ratio1, w.ratio2):
            if not 0.0 < rr <= 1.0:
                raise ModuleError(f"Error 54: {w.name} wave ratio {rr} must lie in (0, 1]")
        if not 0.0 <= w.angle < 360.0:
            raise ModuleError(f"Error 52: incident angle {w.angle} outside [0, 360)")
        if w.type in (ff.WAVE_SV, ff.WAVE_P, ff.WAVE_SH) and np.cos(np.radians(w.angle)) <= 1e-12:
            raise ModuleError(f"{w.name} wave: incident angle {w.angle} deg does not come from below")
    for which in ("ratio1", "ratio2"):
        s = sum(getattr(w, which) for w in inp.waves)
        if abs(s - 1.0) > 1e-6:
            raise ModuleError(f"Error 55: wave ratios ({which}) sum to {s:.6g}, not 1")
    if inp.wopt == 1 and inp.cm != 1:
        raise ModuleError("SH/L waves move along y' only: the control direction must be y' (<cm> = 1)")
    if inp.wopt == 0 and inp.cm == 1:
        raise ModuleError("R/SV/P waves move in the x'z' plane: the control direction cannot be y'")


def p_wave_input(inp: SiteInput) -> bool:
    """True when a P wave field is active in Mode 2 (D-SIT-10: Vp is then also checked)."""
    return inp.mode2 and any(w.type == ff.WAVE_P for w in inp.waves)


def _guideline_warnings(inp: SiteInput, lst) -> None:
    """Modelling guidelines (requirements §4.2 rules).

    G-05 / D-SIT-10, 1/5-wavelength rule: a layer of thickness h passes frequencies up to
    ``V / (5 h)``.  The governing velocity is the (strain-compatible, FILE88) shear-wave velocity
    Vs; for P-wave input the compression-wave velocity Vp is checked as well.
    """
    fmax = float(inp.freq.max())
    checks = [("Vs", inp.vs)] + ([("Vp", inp.vp)] if p_wave_input(inp) else [])
    for name, vel in checks:
        passing = vel[:-1] / (5.0 * inp.thick)
        bad = np.flatnonzero(passing < fmax)
        if bad.size:
            what = " (P-wave input, D-SIT-10)" if name == "Vp" else ""
            lst.warning(f"G-05{what}: {bad.size} layer(s) thicker than {name}/(5 f_max) (f_max = {fmax:.4g} Hz); "
                        f"first: layer {bad[0] + 1}, passing frequency {passing[bad[0]]:.4g} Hz")
    nu = (inp.vp ** 2 - 2 * inp.vs ** 2) / (2 * (inp.vp ** 2 - inp.vs ** 2))
    if np.any(nu > 0.47):
        lst.warning("G-09: Poisson's ratio > 0.47 in layer(s) " + ", ".join(str(i + 1) for i in np.flatnonzero(nu > 0.47)))
    if inp.nI - 1 > 200:
        lst.warning(f"EDU-02: {inp.nI - 1} top layers (V3 limit 200)")
    if inp.nI - 1 <= 20:
        lst.warning(f"{inp.nI - 1} top layers: the manual recommends more than 20 soil layers for accurate "
                    "Rayleigh and Love modes (requirements §4.2 rules)")


# ---------------------------------------------------------------------------------------
# Mode 1
# ---------------------------------------------------------------------------------------
def generated_sublayers(inp: SiteInput, f: float):
    """Generated half-space sublayers at frequency f (D-SIT-02)."""
    if inp.nl == 0:
        return np.zeros(0), tlm.HS_UNIFORM
    return tlm.halfspace_sublayers(inp.vs[-1], f, inp.nl, h_last=inp.thick[-1], vs_last=inp.vs[-2], law=inp.hslaw)


def column(inp: SiteInput, h_gen) -> tlm.Column:
    return tlm.build_column(inp.thick, inp.rho[:-1], inp.G[:-1], inp.M[:-1], inp.rho[-1], inp.G[-1], inp.M[-1],
                            h_gen)


def run_mode1(inp: SiteInput, ctx: ModuleContext) -> Dict[str, np.ndarray]:
    """SITE Mode 1 (requirements §4.2): per frequency, generated half-space sublayers (D-SIT-02),
    TLM column with base dashpots (R1 §2.4) and Rayleigh/Love eigen-solutions (R1 §2.5) -> FILE2."""
    nF = inp.fnum.size
    nl = inp.nl
    h_gen = np.zeros((nF, nl))
    law_used = []
    cols, modes = [], []
    for q, f in enumerate(inp.freq):
        if ctx.cancelled():
            raise ModuleError("run cancelled")
        hg, law = generated_sublayers(inp, float(f))
        h_gen[q] = hg
        law_used.append(law)
        col = column(inp, hg)
        md = tlm.column_modes(col, 2 * np.pi * float(f), MASS)
        if not (np.all(np.isfinite(md.kR)) and np.all(np.isfinite(md.kL))):
            raise ModuleError(f"eigen-solution failed at frequency {f:.6g} Hz")
        cols.append(col)
        modes.append(md)
        ctx.progress(0.7 * (q + 1) / nF, f"Mode 1: frequency {q + 1}/{nF} ({f:.4g} Hz)")
    col0 = cols[0]
    arrays = {
        "fnum": inp.fnum, "freq": inp.freq, "depth_user": inp.depth_user, "layer_thick": inp.thick,
        "layer_rho": inp.rho, "layer_G": inp.G, "layer_M": inp.M, "h_gen": h_gen,
        "iface_user": col0.user_free_index(),
        "kR": np.stack([m.kR for m in modes]), "phix": np.stack([m.phix for m in modes]),
        "phiz": np.stack([m.phiz for m in modes]), "kL": np.stack([m.kL for m in modes]),
        "phiy": np.stack([m.phiy for m in modes]),
        "x_hs_uniform": np.array([1 if lw == tlm.HS_UNIFORM else 0 for lw in law_used], int),
        "x_layer_vs": inp.vs, "x_layer_vp": inp.vp, "x_layer_ds": inp.ds, "x_layer_dp": inp.dp,
    }
    meta = {"df": inp.df, "nl": nl, "base": inp.base, "mass": MASS, "cmodform": inp.cmodform,
            "hslaw": inp.hslaw, "gravity": inp.gravity, "title": inp.title, "layer_no": inp.layer_no,
            "hs_layer": inp.hs_no}
    write_container(ctx.path("FILE2"), "FILE2", arrays, meta, module=NAME)
    _list_mode1(ctx.listing, inp, h_gen, law_used, modes)
    return {"cols": cols, "modes": modes}


FILE2_RTOL = 1e-9


def _check_file2_matches(f2, inp: SiteInput, rows) -> None:
    """Mode 2 alone must use the eigen-solutions of *this* site (requirements §4.2).

    FILE2 stores the resolved layer properties (after the FILE88 substitution of soil mode 1) and
    the generated half-space sublayers; they are compared with the deck to a tight tolerance.
    A FILE2 of another site would otherwise be used silently for vertical and surface waves while
    inclined waves use the deck properties, i.e. FILE1 would mix two sites.
    """
    def differs(a, b) -> bool:
        a, b = np.asarray(a), np.asarray(b)
        if a.shape != b.shape:
            return True
        scale = float(np.max(np.abs(b))) if b.size else 0.0
        return not np.allclose(a, b, rtol=FILE2_RTOL, atol=FILE2_RTOL * scale)

    bad = [name for name, a, b in (("layer thicknesses", f2["layer_thick"], inp.thick),
                                   ("densities", f2["layer_rho"], inp.rho),
                                   ("shear moduli G* (Vs, beta_s)", f2["layer_G"], inp.G),
                                   ("constrained moduli M* (Vp, beta_p)", f2["layer_M"], inp.M))
           if differs(a, b)]
    if int(f2.meta["cmodform"]) != inp.cmodform:
        bad.append("complex modulus form")
    if str(f2.meta["mass"]) != MASS:
        bad.append("mass matrices")
    # the generated sublayers themselves (law D-SIT-02, nl, half-space Vs) -- compared in substance, so
    # a law name whose sublayers coincide (e.g. geometric with uniform fallback) is accepted
    if inp.nl > 0 and any(differs(f2["h_gen"][q], generated_sublayers(inp, float(f2["freq"][q]))[0]) for q in rows):
        bad.append(f"generated half-space sublayers (FILE2 law {f2.meta.get('hslaw', '?')}, deck law {inp.hslaw})")
    if bad:
        raise ModuleError("FILE2 was written for another site (" + ", ".join(bad) + " differ from the deck)"
                          " -- re-run SITE Mode 1")


def _load_file2(inp: SiteInput, ctx: ModuleContext):
    """Mode 2 alone: recover the columns and modes from an existing FILE2 (requirements §4.2)."""
    f2 = read_container(ctx.require("FILE2", "SITE Mode 1"), "FILE2")
    probs = validate(f2)
    if probs:
        raise ModuleError("; ".join(probs))
    if abs(float(f2.meta["df"]) - inp.df) > 1e-9 * inp.df:
        raise ModuleError(f"FILE2 frequency step {f2.meta['df']} differs from the deck ({inp.df})")
    if len(f2["layer_thick"]) != inp.nI - 1 or int(f2.meta["nl"]) != inp.nl:
        raise ModuleError("FILE2 layering differs from the deck -- re-run SITE Mode 1")
    rows = []
    for n in inp.fnum:
        idx = np.flatnonzero(f2["fnum"] == n)
        if idx.size == 0:
            raise ModuleError(f"frequency number {n} not in FILE2 -- re-run SITE Mode 1")
        rows.append(int(idx[0]))
    _check_file2_matches(f2, inp, rows)
    cols = [tlm.column_from_file2(f2, q) for q in rows]
    modes = [tlm.modes_from_file2(f2, q) for q in rows]
    ctx.listing.write(f" Mode 2 uses the eigen-solutions of the existing FILE2 ({len(rows)} frequencies)")
    return {"cols": cols, "modes": modes}


# ---------------------------------------------------------------------------------------
# Mode 2
# ---------------------------------------------------------------------------------------
def wave_field(inp: SiteInput, col: tlm.Column, md: tlm.Modes, w: Wave, f: float) -> ff.WaveField:
    """Raw field of wave ``w`` at frequency f (vertical body waves P0; inclined/surface P1)."""
    omega = 2 * np.pi * f
    if w.type in (ff.WAVE_R, ff.WAVE_L):
        fld, _ = ff.surface_wave(col, md, w.type, w.opt)
        return fld
    if abs(np.sin(np.radians(w.angle))) < 1e-12:
        return ff.vertical_body_wave(col, omega, w.type, MASS)
    user = tlm.build_column(inp.thick, inp.rho[:-1], inp.G[:-1], inp.M[:-1], inp.rho[-1], inp.G[-1], inp.M[-1])
    return ff.inclined_body_wave(user, omega, w.type, w.angle, inp.rho[-1], inp.G[-1], inp.M[-1],
                                 inp.vs[-1], inp.vp[-1], rigid_base=(inp.nl == 0), mass=MASS)


def run_mode2(inp: SiteInput, ctx: ModuleContext, sol) -> None:
    """SITE Mode 2 (requirements §4.2, R1 §2.7): every wave field normalised to unit control motion
    (D-SIT-07), wave ratios interpolated between Frequency 1 and 2 (D-SIT-05) -> FILE1."""
    nW, nF, nI = len(inp.waves), inp.fnum.size, inp.nI
    U = np.zeros((nW, nF, nI, 3), complex)
    k = np.zeros((nW, nF), complex)
    ratio = np.zeros((nW, nF))
    ucp = np.zeros((nW, nF), complex)
    outc = np.zeros((nW, nF), complex)
    choice = {}                      # (wave, frequency) -> tlm.ModeChoice of the surface waves
    if any(abs(np.sin(np.radians(w.angle))) > 1e-12 and w.type in (2, 3, 4) for w in inp.waves) and inp.nl > 0:
        ctx.listing.write(" D-SIT-08: inclined body waves use the exact half-space stiffness; the generated"
                          " half-space sublayers are not used for them")
    for iw, w in enumerate(inp.waves):
        ratio[iw] = ff.wave_ratio(inp.fnum, inp.freq1, inp.freq2, w.ratio1, w.ratio2)
        for q, f in enumerate(inp.freq):
            try:
                fld = wave_field(inp, sol["cols"][q], sol["modes"][q], w, float(f))
                U[iw, q], ucp[iw, q] = ff.normalise(fld, inp.cl, inp.cm)
                if w.type in (ff.WAVE_R, ff.WAVE_L):
                    choice[iw, q] = ff.surface_mode_choice(sol["cols"][q], sol["modes"][q], w.type, w.opt)
            except (ValueError, np.linalg.LinAlgError) as exc:
                raise ModuleError(f"{w.name} wave at {f:.6g} Hz: {exc}") from None
            k[iw, q] = fld.k
            outc[iw, q] = fld.outcrop[inp.cm]
            ctx.progress(0.7 + 0.3 * (iw * nF + q + 1) / (nW * nF), f"Mode 2: {w.name} wave, frequency {q + 1}/{nF}")
    arrays = {"fnum": inp.fnum, "freq": inp.freq, "depth_user": inp.depth_user, "U": U, "k": k, "ratio": ratio,
              "wave_type": np.array([w.type for w in inp.waves], int),
              "x_angle": np.array([w.angle for w in inp.waves]), "x_opt": np.array([w.opt for w in inp.waves], int),
              "x_ucp": ucp, "x_outcrop": outc}
    meta = {"df": inp.df, "cl": inp.cl, "cm": inp.cm, "wopt": inp.wopt, "nl": inp.nl, "base": inp.base,
            "title": inp.title, "normalisation": "within motion at the top of layer cl, direction cm",
            "x_ucp": "raw control-point motion per unit incident wave (body waves) or mode amplitude",
            "x_outcrop": "discrete outcrop motion at the half-space top (direction cm) per unit incident wave"}
    write_container(ctx.path("FILE1"), "FILE1", arrays, meta, module=NAME)
    _list_mode2(ctx.listing, inp, U, k, ratio, ucp, outc, choice)


# ---------------------------------------------------------------------------------------
# Listing
# ---------------------------------------------------------------------------------------
def _sample_rows(n: int, m: int = 3) -> List[int]:
    if n <= m:
        return list(range(n))
    return sorted(set([0, n // 2, n - 1]))


def _list_input(lst, inp: SiteInput) -> None:
    lst.section("SITE input")
    lst.write(f" Title                          : {inp.title}")
    lst.write(f" Mode 1 (FILE2) / Mode 2 (FILE1): {int(inp.mode1)} / {int(inp.mode2)}")
    lst.write(f" Frequency step df             : {inp.df:.9g} Hz")
    lst.write(f" Analysis frequencies          : {inp.fnum.size} (numbers {inp.fnum[0]} .. {inp.fnum[-1]}, "
              f"{inp.freq[0]:.6g} .. {inp.freq[-1]:.6g} Hz)")
    lst.write(f" Base                          : {'rigid (nl = 0)' if inp.nl == 0 else f'half-space, {inp.nl} generated sublayers, law {inp.hslaw}, base dashpots'}")
    lst.write(f" Mass matrices                 : 1/2 lumped + 1/2 consistent (D-CNV-05)")
    lst.write(f" Complex modulus form          : {'1 - 2b^2 + 2ib sqrt(1-b^2)' if inp.cmodform == 0 else '1 + 2ib'}")
    lst.section("Soil layers (TOPL) and half-space")
    cols = ["layer", "L no", "thick", "depth top", "weight", "Vs", "Vp", "beta_s", "beta_p", "nu", "f_pass"]
    rows = []
    depth = inp.depth_user
    for i in range(inp.nI):
        nu = (inp.vp[i] ** 2 - 2 * inp.vs[i] ** 2) / (2 * (inp.vp[i] ** 2 - inp.vs[i] ** 2))
        if i < inp.nI - 1:
            rows.append([str(i + 1), str(inp.layer_no[i]), inp.thick[i], depth[i], inp.weight[i], inp.vs[i], inp.vp[i],
                         inp.ds[i], inp.dp[i], nu, inp.vs[i] / (5 * inp.thick[i])])
        else:
            rows.append(["HS", str(inp.hs_no), "-", depth[i], inp.weight[i], inp.vs[i], inp.vp[i], inp.ds[i], inp.dp[i],
                         nu, "-"])
    if p_wave_input(inp):                       # D-SIT-10: Vp governs the P-wave input
        cols.append("f_pass P")
        for i, row in enumerate(rows):
            row.append(inp.vp[i] / (5 * inp.thick[i]) if i < inp.nI - 1 else "-")
    lst.table(cols, rows, fmt="{:>14.5g}")
    lst.write(" f_pass = Vs/(5 h): highest frequency passed by the layer (1/5-wavelength rule, G-05)")
    if p_wave_input(inp):
        lst.write(" f_pass P = Vp/(5 h): the same rule for the P-wave input (D-SIT-10)")
    if inp.mode2:
        lst.section("Wave environment (Mode 2)")
        lst.write(f" Control point: top of layer {inp.cl} (depth {depth[inp.cl - 1]:.6g}), direction "
                  f"{ff.DIRECTION_NAMES[inp.cm]}, within motion (D-SIT-07)")
        lst.write(f" Wave-ratio frequencies: numbers {inp.freq1} and {inp.freq2} "
                  f"({inp.freq1 * inp.df:.6g} and {inp.freq2 * inp.df:.6g} Hz)")
        lst.table(["wave", "opt", "ratio 1", "ratio 2", "angle"],
                  [[w.name, str(w.opt), w.ratio1, w.ratio2, w.angle] for w in inp.waves], fmt="{:>14.6g}")


def _list_mode1(lst, inp: SiteInput, h_gen, law_used, modes) -> None:
    if inp.nl > 0:
        lst.section("Generated half-space sublayers (variable-depth method, D-SIT-02)")
        rows = []
        for q in range(inp.fnum.size):
            hg = h_gen[q]
            rows.append([float(inp.freq[q]), float(hg.sum()), law_used[q], float(hg[0]), float(hg[-1]),
                         float(inp.vs[-1] / inp.freq[q] / 8.0)])
        lst.table(["f (Hz)", "total 1.5L", "law used", "h first", "h last", "cap L/8"], rows, fmt="{:>14.6g}")
    loss = tlm.material_loss(inp.G, inp.M)
    ratio = tlm.propagating_ratio(loss)
    lst.section(f"Mode 1: propagating Rayleigh and Love modes (|Im k| <= {ratio:.4g} Re k, slowest first)")
    lst.write(f" Propagating sector |arg k| <= atan(0.5) + atan({loss:.4g}): the elastic sector widened by the "
              "largest material loss angle of the column (D-SIT-09)")
    for q in _sample_rows(inp.fnum.size):
        f = float(inp.freq[q])
        w = 2 * np.pi * f
        md = modes[q]
        for name, k in (("Rayleigh", md.kR), ("Love", md.kL)):
            prop = np.flatnonzero(tlm.propagating(k, ratio))
            prop = prop[np.argsort(-k.real[prop])]
            lst.write(f" f = {f:.6g} Hz, {name}: {prop.size} propagating of {k.size} modes")
            rows = [[str(int(j) + 1), k[j].real, k[j].imag, w / k[j].real, 2 * np.pi / k[j].real]
                    for j in prop[:5]]
            if rows:
                lst.table(["mode", "Re k", "Im k", "c = w/Re k", "wavelength"], rows, fmt="{:>14.6g}")


def _list_surface_modes(lst, inp: SiteInput, choice) -> None:
    """Which mode each surface-wave field uses (D-SIT-09), with ties of the least-decay rule."""
    if not choice:
        return
    lst.section("Mode 2: surface-wave mode selection (D-SIT-09)")
    rows, ties = [], 0
    for (iw, q), c in sorted(choice.items()):
        w = inp.waves[iw]
        rule = "least decay" if c.rule == 2 else "shortest wavelength"
        rows.append([w.name, float(inp.freq[q]), rule, str(c.index + 1), str(c.n_candidates), str(c.n_tied)])
        ties += c.n_tied > 1
    lst.table(["wave", "f (Hz)", "rule", "mode", "candidates", "tied"], rows, fmt="{:>14.6g}")
    if ties:
        lst.write(f" {ties} selection(s) of the least-decay rule found several modes with the same |Im k| "
                  f"(within {tlm.TIE_RTOL:.0e} max|k|, e.g. undamped soil): the mode with the largest Re k "
                  "(shortest wavelength) is used")


def _list_mode2(lst, inp: SiteInput, U, k, ratio, ucp, outc, choice=None) -> None:
    _list_surface_modes(lst, inp, choice or {})
    lst.section("Mode 2: free-field amplitudes at the user interfaces (unit control motion)")
    tot = np.einsum("wf,wfic->fic", ratio, U)
    depth = inp.depth_user
    for q in range(inp.fnum.size):
        f = float(inp.freq[q])
        lst.write("")
        info = "  ".join(f"{w.name}: r={ratio[i, q]:.4g} k={k[i, q].real:.5g}{k[i, q].imag:+.3g}i"
                         for i, w in enumerate(inp.waves))
        lst.write(f" f = {f:.6g} Hz (number {inp.fnum[q]})   {info}")
        rows = [[str(i + 1), depth[i], abs(tot[q, i, 0]), abs(tot[q, i, 1]), abs(tot[q, i, 2])] for i in range(inp.nI)]
        lst.table(["interface", "depth", "|U x'|", "|U y'|", "|U z'|"], rows, fmt="{:>14.6g}")
        for i, w in enumerate(inp.waves):
            if w.type in (2, 3, 4) and abs(outc[i, q]) > 0:
                surf = U[i, q, 0, inp.cm] * ucp[i, q] / outc[i, q]
                lst.write(f"   {w.name}: |surface / outcrop at half-space top| = {abs(surf):.6g}")


# ---------------------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------------------
def run(ctx: ModuleContext) -> int:
    d = decks.read(ctx.deck_path, NAME)
    lst = ctx.listing
    inp = read_input(d, ctx.workdir, lst)
    _list_input(lst, inp)
    if int(d["opmode"]) == 1:
        lst.write("")
        lst.write(" Data check only (<opmode> = 1): no FILE1/FILE2 written")
        return 0
    if inp.mode1:
        sol = run_mode1(inp, ctx)
    else:
        sol = _load_file2(inp, ctx)
    if inp.mode2:
        run_mode2(inp, ctx, sol)
    ctx.progress(1.0, "SITE done")
    return 0


if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
