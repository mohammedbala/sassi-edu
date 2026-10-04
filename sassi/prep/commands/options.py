"""Module analysis-option commands (requirements sections 3.3 and 3.4.F; spec 07 section 4).

Record setters replace the module's option record (every positional argument stored as typed);
indexed setters create or overwrite the entry at their key; string setters set a string.  The
records are the typed classes of :mod:`sassi.prep.options` (argument names of the manual, new-model
defaults, dialog labels), so the GUI dialogs read and write the same storage by field name.

Every command validates its arguments with the CHECK rules of the fields it owns (UI-06): a
non-numeric value in a numeric field skips the command (error); a value that CHECK will reject is
stored and reported as a warning quoting the Chapter 10 message, so that the model can be
completed in any order (CHECK and AFWRITE are the authority).

The list and request commands (FREQ, DAMP, TOPL, AMP, NOUT, EOUT, RDND) are in
:mod:`sassi.prep.commands.session`; the extension X-commands in
:mod:`sassi.prep.commands.extensions`.
"""
from __future__ import annotations

from typing import List

from ...model.values import NumberError, parse_int
from ..check import ERRORS, EDU_TEXT, WARNINGS, format_title
from ..options import AOPT_MODULES, problem_fmt, record_class
from ..registry import command, indexed_key


# ======================================================================================
# helpers
# ======================================================================================
def problem_text(kind: str, num, detail: str = "", **fmt) -> str:
    """``Error 47 : Illegal Number of Layers for Halfspace Simulation (detail)`` (manual titles;
    placeholders without a value are kept as ``<i>``, ``<w>`` ...)."""
    if isinstance(num, int):
        title = (ERRORS if kind == "Error" else WARNINGS).get(num, "")
    else:
        title = EDU_TEXT.get(num, ("", ""))[1]
    return f"{kind} {num} : {format_title(title, fmt)}" + (f" ({detail})" if detail else "")


def report_problems(c, rec) -> None:
    """Warn about the field-local CHECK rules violated by a stored record (UI-06).  The title
    placeholders come from the rule itself (:class:`sassi.prep.options.Problem`), as in CHECK."""
    for p in rec.problems():
        kind, num, det = p
        c.warn(problem_text(kind, num, det, **problem_fmt(p)) + " -- CHECK will report it")


def _announce_tier(c) -> None:
    spec = c.spec
    if spec.tier != "P0" and not spec.silent:
        c.warn(f"{spec.name} is not available in this build (tier {spec.tier}): stored only")


def build_record(c, name: str, tokens: List[str]):
    """Typed record of ``tokens``; extra arguments dropped with a warning, bad numbers rejected."""
    cls = record_class(name)
    nf = len(cls.FIELDS)
    toks = list(tokens)
    n_given = len(toks)
    while n_given and not toks[n_given - 1].strip():
        n_given -= 1
    if n_given > nf:
        c.warn(f"arguments after argument {nf} ignored")
        toks = toks[:nf]
    rec = cls(name, toks)
    errs = rec.token_errors()
    if errs:
        c.fail(errs[0])
    for w in rec.rounded_ints():
        c.warn(w)
    return rec


def store_record(c, name: str):
    """Record setter: the command replaces the record of ``name``."""
    rec = build_record(c, name, c.tokens)
    c.model.options.set_record(rec)
    report_problems(c, rec)
    _announce_tier(c)
    c.confirm(f"{name} options set ({len(rec)} arguments)")
    return rec


def store_entry(c, name: str, delete_when_blank: int = 0):
    """Indexed setter: create or overwrite the entry at its key (key rule of the catalogue).

    ``delete_when_blank`` = k: a blank argument k removes the entry (file setters RSIN ... TPSD).
    """
    toks = list(c.tokens)
    key = indexed_key(c.spec, toks)
    if delete_when_blank and not (toks[delete_when_blank - 1].strip() if len(toks) >= delete_when_blank else ""):
        existed = c.model.options.delete_entry(name, key)
        c.confirm(f"{name} {key} cleared" if existed else f"{name} {key} not defined (nothing cleared)")
        return None
    rec = build_record(c, name, toks)
    c.model.options.set_entry(name, key, rec)
    report_problems(c, rec)
    _announce_tier(c)
    c.confirm(f"{name} {key} set")
    return rec


def store_string(c, name: str) -> None:
    value = c.text(1) if c.spec.text_from == 1 else c.str(1)
    c.model.options.set_string(name, value)
    _announce_tier(c)
    c.confirm(f"{name}: {value}" if value else f"{name} cleared")


# ======================================================================================
# EQUAKE tab
# ======================================================================================
@command("EQUAKE")
def cmd_equake(c):
    """EQUAKE,<accopt>,<nrfreq>,<rand>,<damp>,<dur>,<corr>,<seeds>,[tpsd]: EQUAKE options (D-EQK-01)."""
    rec = store_record(c, "EQUAKE")
    if rec.accopt not in (0, 1, 2):
        c.warn("<accopt> must be 0 (none), 1 (seed record) or 2 (external history), D-EQK-01")


def _spectrum_file(c, name: str) -> None:
    no = c.int(1, required=True, what="no")
    if not 1 <= no <= 3:
        c.warn(f"spectrum number {no}: the EQUAKE dialog has spectra 1..3 (one per direction)")
    store_entry(c, name, delete_when_blank=2)


@command("RSIN")
def cmd_rsin(c):
    """RSIN,<no>,<file>: target response-spectrum file of spectrum <no> (blank file clears it)."""
    _spectrum_file(c, "RSIN")


@command("RSOUT")
def cmd_rsout(c):
    """RSOUT,<no>,<file>: response-spectrum output file of spectrum <no>."""
    _spectrum_file(c, "RSOUT")


@command("ACCIN")
def cmd_accin(c):
    """ACCIN,<no>,<file>: seed / external acceleration input file of spectrum <no>."""
    _spectrum_file(c, "ACCIN")


@command("ACCOUT")
def cmd_accout(c):
    """ACCOUT,<no>,<file>: generated acceleration output file of spectrum <no>."""
    _spectrum_file(c, "ACCOUT")


@command("TPSD")
def cmd_tpsd(c):
    """TPSD,<num>,<file>: target PSD file of spectrum <num>."""
    _spectrum_file(c, "TPSD")


@command("CORR")
def cmd_corr(c):
    """CORR,<no>,<time>,<val>: X-Y correlation pair <no> (P1)."""
    store_entry(c, "CORR")


@command("EQTIT", text_from=1)
def cmd_eqtit(c):
    """EQTIT,<title>: spectra title (rest of the line)."""
    store_string(c, "EQTIT")


# ======================================================================================
# SOIL tab
# ======================================================================================
@command("SOIL", paren=True)
def cmd_soil(c):
    """SOIL,<nrval>,<grav>,<header>,<outcrop>,<save>,<iter>,<ratio>,<gravmult>,<cof>: SOIL options.

    A final parenthesised token is the legacy PREP format (rule L9, D-SOL-10): the arguments are
    stored with the format string but not mapped (CHECK Error 102).
    """
    legacy = c.lexed.legacy_paren
    if legacy:
        rec = record_class("SOIL")("SOIL", list(c.tokens), legacy=legacy)
        c.model.options.set_record(rec)
        c.warn(f"legacy (PREP) format {legacy} detected; the arguments are stored but not mapped to the "
               f"current argument order (Error 102)")
        return
    rec = store_record(c, "SOIL")
    if rec.given(9) and rec.cof != 0:
        c.warn("<cof> is disabled: the cut-off is always the Nyquist frequency 1/(2 dt)")


@command("SPRO", text_from=3)
def cmd_spro(c):
    """SPRO,<layer>,<prop>,<dynprop>: SOIL sublayer -> L property and DYNP label (last = half-space)."""
    layer = c.int(1, required=True, what="layer")
    if layer < 1:
        c.fail("sublayer numbers start at 1")
    store_entry(c, "SPRO")


@command("DYNP", text_from=6)
def cmd_dynp(c):
    """DYNP,<no>,<sg>,<g>,<sd>,<d>,<label>: curve point (strain %, G/Gmax; strain %, damping %)."""
    no = c.int(1, required=True, what="no")
    if no < 1:
        c.fail("point numbers start at 1")
    store_entry(c, "DYNP")


def _soil_request(c, name: str) -> None:
    layer = c.int(1, required=True, what="layer")
    if layer < 1:
        c.fail("sublayer numbers start at 1")
    store_entry(c, name)


@command("SACC")
def cmd_sacc(c):
    """SACC,<layer>,<opt>,<outcrop>: acceleration output at the top of a SOIL sublayer."""
    _soil_request(c, "SACC")


@command("SRS")
def cmd_srs(c):
    """SRS,<layer>,<save>,<outcrop>: response-spectrum output at the top of a SOIL sublayer."""
    _soil_request(c, "SRS")


@command("SSTR")
def cmd_sstr(c):
    """SSTR,<layer>,<opt1>,<opt2>,<opt3>,<opt4>: stress/strain output at the centre of a sublayer."""
    _soil_request(c, "SSTR")


@command("SSAF", text_from=7)
def cmd_ssaf(c):
    """SSAF,<layer>,<save>,<outcrop1>,<outcrop2>,<layer2>,<freqstep>,<title>: spectral amplification."""
    _soil_request(c, "SSAF")


@command("SFOU")
def cmd_sfou(c):
    """SFOU,<layer>,<out>,<save>,<outcrop>,<smooth>,<nrval>: parsed and stored only."""
    _soil_request(c, "SFOU")
    c.warn("SFOU is not usable in this version (stored only; compute Fourier spectra with EQUAKE)")


@command("THFILE", text_from=1)
def cmd_thfile(c):
    """THFILE,<file>: control-motion acceleration file (shared by SOIL, MOTION, STRESS, RELDISP)."""
    store_string(c, "THFILE")


@command("THTIT", text_from=1)
def cmd_thtit(c):
    """THTIT,<title>: control-motion title."""
    store_string(c, "THTIT")


# ======================================================================================
# SITE, POINT, HOUSE, FORCE, ANALYS
# ======================================================================================
@command("SITE")
def cmd_site(c):
    """SITE,<opmode>,<mode1>,<fstep>,<nl>,<hs>,<mode2>,<wopt>,<freq1>,<freq2>,<cl>,<cm>,<delt>,<nft>,<freq>."""
    rec = store_record(c, "SITE")
    if rec.cm not in (0, 1, 2):
        c.warn("<cm> must be 0 (x'), 1 (y') or 2 (z')")
    if rec.wopt not in (0, 1):
        c.warn("<wopt> must be 0 (R-, SV-, P-waves) or 1 (SH- and L-waves)")


@command("WAVE")
def cmd_wave(c):
    """WAVE,<type>,<opt>,<ratio1>,<ratio2>,<angle>: wave field (1 R, 2 SV, 3 P, 4 SH, 5 L)."""
    t = c.int(1, required=True, what="type")
    if t not in (1, 2, 3, 4, 5):
        c.fail("<type> must be 1 (R), 2 (SV), 3 (P), 4 (SH) or 5 (L)")
    opt = c.int(2, default=1)
    if opt == 2 and t != 1:
        c.warn("<opt> = 2 (least decay) applies to Rayleigh waves only")
    store_entry(c, "WAVE")


@command("POINT")
def cmd_point(c):
    """POINT,<opmode>,<layer>,<rad>: embedment layers (0 surface) and central-zone radius."""
    store_record(c, "POINT")


@command("HOUSE")
def cmd_house(c):
    """HOUSE,<gravity>,<gelev>,<opmode>,<dim>,<imp>,<coh>,<wpass>,<me>,<cmplxspec>: HOUSE options."""
    rec = store_record(c, "HOUSE")
    if rec.dim not in (1, 2):
        c.warn("<dim> must be 1 (2D, POINT2) or 2 (3D, POINT3); 1D is not available")
    if rec.coh or rec.wpass or rec.me:
        c.warn("incoherency, wave passage and multiple excitation are tier P1 features")


@command("INCOH")
def cmd_incoh(c):
    """INCOH,<gammax>,<gammay>,<gammaz>,<alpha>,<ngp>,<ipr>,<nmodes>,<met>,<HSeed>,<VSeed>,<RandPhz> (P1)."""
    store_record(c, "INCOH")


@command("WPASS")
def cmd_wpass(c):
    """WPASS,<appv>,<ang>,<cohf>: wave passage (P1); <cohf> = unlagged coherency model 1-7 (D-INC-01)."""
    store_record(c, "WPASS")


@command("ME")
def cmd_me(c):
    """ME,<no>,<nfirst>,<nlast>,<xc>,<yc>,<zc>: multiple-excitation zone (P1)."""
    no = c.int(1, required=True, what="no")
    if no > 10:
        c.warn(f"motion number {no} > 10 (manual range 1-10; accepted)")
    store_entry(c, "ME")


@command("SYMM")
def cmd_symm(c):
    """SYMM,<no>,[<type>],[<node1>],[<node2>],[<node3>]: symmetry plane; <node1> = 0 resets plane <no>."""
    toks = list(c.tokens)
    key = indexed_key(c.spec, toks)
    n1 = toks[2].strip() if len(toks) > 2 else ""
    try:
        reset = (not n1) or parse_int(n1)[0] == 0
    except NumberError:
        reset = False
    if reset:
        existed = c.model.options.delete_entry("SYMM", key)
        c.confirm(f"SYMM {key} reset" if existed else f"SYMM {key} not defined (nothing reset)")
        return
    if key not in (1, 2):
        c.warn("symmetry plane numbers are 1 or 2 (at most two planes)")
    store_entry(c, "SYMM")


@command("FORCE")
def cmd_force(c):
    """FORCE,<opmode>: FORCE options (gravity, df, delt, NFFT and the frequency set are shared)."""
    store_record(c, "FORCE")


@command("ANALYS")
def cmd_analys(c):
    """ANALYS,<opmode>,<type>,<mode>,<save>,<prnt>,<fopt>,<ang>,<xc>,<yc>,<zc>,<impe>,[simul]."""
    rec = store_record(c, "ANALYS")
    if rec.type not in (0, 1):
        c.warn("<type> must be 0 (seismic) or 1 (foundation vibration)")
    if rec.mode not in (0, 1, 2, 3):
        c.warn("<mode> must be 0 (initiation), 1 (new structure), 2 (new seismic environment) or "
               "3 (new dynamic loading)")


# ======================================================================================
# MOTION, STRESS, RELDISP
# ======================================================================================
@command("MOTION")
def cmd_motion(c):
    """MOTION,<opmode>,<out>,<step>,<dur>,<res>,<freq1>,<freq2>,<fstep>,<mult>,<max>,<rec1>,<rec2>,<fopt>,
    <bl>,<smo>,<cplx>,<cnvrt>,<pzadj>,<interp>: MOTION options (19 arguments)."""
    rec = store_record(c, "MOTION")
    if rec.interp not in range(7):
        c.warn("<interp> must be 0..6")
    if rec.smo != 0 and rec.interp == 6:
        c.warn("smoothing must be 0 with interpolation option 6 (EDU-14)")


@command("STRESS")
def cmd_stress(c):
    """STRESS,<opmode>,<iter>,<save>,<itran>,<interopt>: STRESS options."""
    rec = store_record(c, "STRESS")
    if rec.interopt not in range(7):
        c.warn("<interopt> must be 0..6")


@command("RELD")
def cmd_reld(c):
    """RELD,<RelDisOutput>,<RelDispSAll>,<RelDispNumFiles>: RELDISP options (D-RDP-03)."""
    rec = store_record(c, "RELD")
    if rec.reldisoutput & 2:
        c.warn("<RelDisOutput> bit 1: relative displacements from .TFU (amplitudes) are discouraged (D-RDP-03)")


@command("RELFILE", text_from=1)
def cmd_relfile(c):
    """RELFILE,<FileName>: reference-node complex .TFI file of RELDISP."""
    store_string(c, "RELFILE")


# ======================================================================================
# AFWRITE tab
# ======================================================================================
@command("AOPT")
def cmd_aopt(c):
    """AOPT,<EQUAKE>,<SOIL>,<DEP>,<SITE>,<POINT>,<HOUSE>,<DEP>,<FORCE>,<ANALYS>,<COMBIN>,<MOTION>,<STRESS>,
    <RELDISP>,<PANEL>: modules processed by CHECK and AFWRITE (a non-zero flag includes the module)."""
    rec = build_record(c, "AOPT", c.tokens)
    for k in (3, 7):
        if rec.integer(k, 0):
            c.warn(f"argument {k} (DEP: module not in this version) must be 0; ignored")
    c.model.options.set_record(rec)
    on = rec.enabled()
    c.confirm("AOPT: " + (", ".join(on) if on else "no module enabled"))


__all__ = ["store_record", "store_entry", "store_string", "build_record", "report_problems", "problem_text"]

# the AOPT module order (exported for the GUI AFWRITE tab)
AOPT_ORDER = [m for m in AOPT_MODULES if m is not None]
