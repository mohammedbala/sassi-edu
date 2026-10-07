"""Defaults of blank analysis inputs from the built-in library (requirements section 7.19, D-W5-03 ... D-W5-14).

ACS SASSI requires the input files of an analysis: a blank THFILE is Error 73, no RSIN file Error 84, blank
RSOUT / ACCOUT Errors 86 / 87, a SOIL profile without curves Errors 95 / 97 / 99.  SASSI-EDU fills these
blanks from the built-in library (:mod:`sassi.io.library`) so that a new model runs before its project
inputs exist -- and **reports every default it uses** (CHECK Warning EDU-29, AFWRITE notes, the RUN<MODULE>
messages and the module listing), so the engineer always sees what input was assumed:

=======================  ============================================================================
blank input              default
=======================  ============================================================================
THFILE (MOTION, STRESS,  ``@rg160h_030g.acc`` (RG 1.60 record, 0.30 g, 20 s) for a seismic analysis
RELDISP)                 (ANALYS ``<type>`` 0); ``@ricker_5hz.th`` (5 Hz Ricker pulse) for a forced-
                         vibration analysis (``<type>`` 1).  Only with MOTION ``<fopt>`` 0 and SITE
                         ``<delt>`` = 0.005 s (the files' layout and time step); otherwise Error 73
                         stays, with the reason (D-W5-03 ... D-W5-05)
THFILE (SOIL; SOILX      ``@rg160h_030g.acc`` (SOIL is a seismic site response), same time-step
file blank too)          condition (D-W5-06)
RSIN (no RSIN file)      ``RSIN,1,@rg160h_030g.rsi`` (RG 1.60 H, 0.30 g, 5 %, 27 rows); a blank
                         Number of Frequencies takes its 27 rows (D-W5-07)
RSOUT / ACCOUT i         ``<model>_eq<i>.rso`` / ``<model>_eq<i>.acc`` in the model folder for every
                         spectrum i that runs (D-W5-08)
DYNP label not defined   the library curve Clay, Sand or Rock of that label (a model DYNP of the same
                         label wins; labels are case-sensitive) (D-W5-09)
no SPRO entry at all     the SITE profile: SPRO k = TOPL layer k with Sand (Vs < 2493 ft/s, 760 m/s)
                         or Rock (Vs >= 2493 ft/s), the SITE half-space ``<hs>`` as the last, linear
                         sublayer (D-W5-10)
=======================  ============================================================================

Not defaulted (ACS behaviour kept): a file that is *given* but missing (Errors 73 / 85 / 88 ...), ACCIN of
EQUAKE ``<accopt>`` 1 / 2 (Error 88), SPRO entries given without curves (Error 95), the numeric options.
``EDUOPT,DEFAULTS,OFF`` switches every default off (the ACS behaviour, D-W5-12); the ``@`` names themselves
always work.  The functions here are the single policy used by CHECK and AFWRITE
(:class:`sassi.prep.check.Resolved`), the RUN commands, LIBRARY,DEFAULTS and the GUI dialogs.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from ..conventions import unit_system
from ..io import library as LIB
from ..model.values import fmt_num
from .options import eduopt, get_entries, get_record, record_class

#: the default files (D-W5-03, D-W5-04, D-W5-07)
SEISMIC_RECORD = "@rg160h_030g.acc"
VIBRATION_LOAD = "@ricker_5hz.th"
TARGET_SPECTRUM = "@rg160h_030g.rsi"
#: shear-wave velocity from which a default SOIL sublayer is Rock (D-W5-10): 760 m/s (site class B/C
#: boundary of ASCE 7), 760 / 0.3048 = 2493.4 ft/s in British units
ROCK_VS = {"SI": 760.0, "BS": 760.0 / 0.3048}
#: relative tolerance of the time-step match (MOTION's DT_RTOL)
DT_RTOL = 1e-6
#: the label of the CHECK message reporting a default (sassi.prep.check.EDU_TEXT)
EDU_DEFAULT = "EDU-29"

_SHORT = {SEISMIC_RECORD: "RG 1.60 record", VIBRATION_LOAD: "5 Hz Ricker pulse"}
_DETAIL = {SEISMIC_RECORD: "0.30 g, 20 s, dt 0.005 s", VIBRATION_LOAD: "peak 1, 2 s, dt 0.005 s"}


@dataclass(frozen=True)
class DefaultUse:
    """One default that an analysis uses: the blank ``item``, the ``value`` used and the report line."""
    item: str             # 'THFILE', 'RSIN 1', 'RSOUT 1', 'ACCOUT 1', 'DYNP Sand', 'SPRO', 'SOIL header'
    value: str            # '@rg160h_030g.acc', 'm_eq1.rso', ...
    text: str             # the report line (CHECK EDU-29 detail, AFWRITE note, RUN message, listing)


@dataclass(frozen=True)
class HistoryChoice:
    """The history file a module reads: ``name`` ('' when none); ``default`` when ``name`` is a default;
    ``why_not`` when the input is blank and no default applies (the detail of Error 73)."""
    name: str
    default: Optional[DefaultUse] = None
    why_not: str = ""


def enabled(model) -> bool:
    """False with ``EDUOPT,DEFAULTS,OFF`` (ACS behaviour: blank inputs are errors, D-W5-12)."""
    return eduopt(model, "DEFAULTS") != "OFF"


def _dt_mismatch(name: str, delt: float) -> str:
    dt = LIB.history_dt(name)
    if dt is None or delt <= 0:
        return ""
    if abs(dt - delt) > DT_RTOL * max(abs(delt), 1e-30):
        return (f"the built-in default {name} has dt = {fmt_num(dt)} s, the time step SITE <delt> is "
                f"{fmt_num(delt)} s: give THFILE")
    return ""


def control_history(model, given: str, seismic: bool, fopt: int, delt: float) -> HistoryChoice:
    """THFILE of MOTION, STRESS and RELDISP (D-W5-03 ... D-W5-05).

    A given name is used as it is.  A blank one is the RG 1.60 record (seismic) or the Ricker pulse
    (forced vibration, ANALYS ``<type>`` 1), provided MOTION ``<fopt>`` is 0 (the files hold dt on
    line 1) and the time step is the files' 0.005 s; otherwise no default (Error 73 with the reason)."""
    g = (given or "").strip()
    if g:
        return HistoryChoice(g)
    if not enabled(model):
        return HistoryChoice("", why_not="EDUOPT,DEFAULTS,OFF: no built-in default")
    name = SEISMIC_RECORD if seismic else VIBRATION_LOAD
    if int(fopt) != 0:
        return HistoryChoice("", why_not=f"the built-in default {name} has the <fopt> 0 layout (dt on line 1), "
                                         f"MOTION <fopt> = {int(fopt)}: give THFILE")
    why = _dt_mismatch(name, delt)
    if why:
        return HistoryChoice("", why_not=why)
    if seismic:
        text = f"THFILE blank: the built-in {_SHORT[name]} {name} ({_DETAIL[name]}) is used"
    else:
        text = (f"THFILE blank (forced vibration, ANALYS <type> 1): the built-in {_SHORT[name]} {name} "
                f"({_DETAIL[name]}) is the load history")
    return HistoryChoice(name, DefaultUse("THFILE", name, text))


def soil_history(model, soilx_file: str, thfile: str, delt: float) -> HistoryChoice:
    """SOIL input motion: SOILX ``<file>``, else THFILE, else the RG 1.60 record (D-W5-06)."""
    g = (soilx_file or "").strip() or (thfile or "").strip()
    if g:
        return HistoryChoice(g)
    if not enabled(model):
        return HistoryChoice("", why_not="EDUOPT,DEFAULTS,OFF: no built-in default")
    name = SEISMIC_RECORD
    why = _dt_mismatch(name, delt)
    if why:
        return HistoryChoice("", why_not=why)
    text = (f"THFILE blank: the built-in {_SHORT[name]} {name} ({_DETAIL[name]}) is the SOIL input motion")
    return HistoryChoice(name, DefaultUse("THFILE", name, text))


def soil_header_note(name: str, header: int) -> Optional[DefaultUse]:
    """A built-in record holds its time step on line 1; SOIL reads values after ``<header>`` lines, so with
    ``<header>`` = 0 it skips that line itself (D-W5-06).  The rule is reported like a default."""
    if not LIB.is_library_name(name) or LIB.history_dt(name) is None or int(header) != 0:
        return None
    return DefaultUse("SOIL header", "1", f"SOIL <header> = 0 with the built-in record {name.strip()}: its first "
                                         f"line is the time step, SOIL skips it (read as <header> = 1)")


# ======================================================================================
# EQUAKE
# ======================================================================================
FILE_NAMES = ("RSIN", "RSOUT", "ACCIN", "ACCOUT", "TPSD")


def equake_files(model, accopt: int) -> Tuple[Dict[str, Dict[int, str]], List[DefaultUse]]:
    """The EQUAKE spectrum files ``{name: {spectrum: file}}`` with the defaults applied (D-W5-07, D-W5-08),
    and the defaults used.

    * no RSIN file and ``<accopt>`` != 2: RSIN 1 = the RG 1.60 0.30 g spectrum;
    * every spectrum that runs (an RSIN file; for ``<accopt>`` 2 an ACCIN or RSOUT entry): a blank RSOUT
      is ``<model>_eq<i>.rso`` and (``<accopt>`` != 2) a blank ACCOUT ``<model>_eq<i>.acc``.
    ACCIN and TPSD have no default."""
    files: Dict[str, Dict[int, str]] = {}
    for name in FILE_NAMES:
        files[name] = {}
        for k, rec in get_entries(model, name):
            if isinstance(k, int):
                files[name][k] = (rec.file or "").strip()
    uses: List[DefaultUse] = []
    if not enabled(model):
        return files, uses
    accopt = int(accopt)
    if accopt != 2 and not any(v for v in files["RSIN"].values()):
        files["RSIN"][1] = TARGET_SPECTRUM
        e = LIB.entry(TARGET_SPECTRUM)
        uses.append(DefaultUse("RSIN 1", TARGET_SPECTRUM,
                               f"RSIN blank: the built-in RG 1.60 H spectrum {TARGET_SPECTRUM} (0.30 g, 5 %, "
                               f"{e.values if e else 27} frequencies) is the target of spectrum 1"))
    if accopt != 2:
        active = sorted(k for k, v in files["RSIN"].items() if v)
    else:
        active = sorted(set(files["ACCIN"]) | set(files["RSOUT"]))
    stem = (getattr(model, "name", "") or "equake").strip() or "equake"
    for i in active:
        if not files["RSOUT"].get(i):
            files["RSOUT"][i] = f"{stem}_eq{i}.rso"
            uses.append(DefaultUse(f"RSOUT {i}", files["RSOUT"][i],
                                   f"RSOUT {i} blank: the output spectrum is written to {files['RSOUT'][i]} "
                                   f"(model folder)"))
        if accopt != 2 and not files["ACCOUT"].get(i):
            files["ACCOUT"][i] = f"{stem}_eq{i}.acc"
            uses.append(DefaultUse(f"ACCOUT {i}", files["ACCOUT"][i],
                                   f"ACCOUT {i} blank: the generated record is written to {files['ACCOUT'][i]} "
                                   f"(.vel, .dis, .psd, .fft next to it; model folder)"))
    return files, uses


def default_output_name(model, kind: str, i: int) -> str:
    """``<model>_eq<i>.rso`` / ``.acc`` (D-W5-08)."""
    stem = (getattr(model, "name", "") or "equake").strip() or "equake"
    return f"{stem}_eq{i}.{'rso' if kind.upper() == 'RSOUT' else 'acc'}"


# ======================================================================================
# SOIL profile and curves
# ======================================================================================
def model_dynp_labels(model) -> List[str]:
    return list(dict.fromkeys(lab for (lab, _), _ in get_entries(model, "DYNP")))


def soil_unit_system(model) -> str:
    """Unit system of the SOIL column: SOIL ``<grav>`` when SOIL was given, else the HOUSE gravity."""
    if model.options.record("SOIL") is not None:
        return unit_system(float(get_record(model, "SOIL").grav))
    return unit_system(float(get_record(model, "HOUSE").gravity))


def soil_profile(model) -> Tuple[List[Tuple[int, object]], Optional[DefaultUse], str]:
    """SOIL sublayers ``[(k, SPRO record)]`` (top first, the last is the half-space), the default used and,
    when there is no profile, why the default does not apply.

    Stored SPRO entries are used as they are.  With no SPRO entry at all the SITE profile is the default
    (D-W5-10): sublayer k = TOPL layer k, labelled Rock when its Vs >= 2493 ft/s (760 m/s) and Sand
    otherwise, then the SITE half-space layer ``<hs>`` without a label (linear)."""
    stored = [(k, rec) for k, rec in get_entries(model, "SPRO") if isinstance(k, int)]
    if stored:
        return stored, None, ""
    if not enabled(model):
        return [], None, "EDUOPT,DEFAULTS,OFF: no default profile"
    topl = list(model.topl)
    if not topl:
        return [], None, "no SPRO and no TOPL layers for the default profile"
    missing = [l for l in dict.fromkeys(topl) if l not in model.layers]
    if missing:
        return [], None, f"no SPRO; TOPL layer {missing[0]} is not defined (default profile)"
    hs = int(get_record(model, "SITE").hs)
    if hs not in model.layers:
        return [], None, f"no SPRO; the SITE half-space layer <hs> = {hs} is not defined (default profile)"
    us = soil_unit_system(model)
    thr = ROCK_VS[us]
    cls = record_class("SPRO")
    prof: List[Tuple[int, object]] = []
    labels: List[str] = []
    for k, l in enumerate(topl, start=1):
        lab = "Rock" if float(model.layers[l].vs) >= thr - 1e-9 else "Sand"
        labels.append(lab)
        prof.append((k, cls("SPRO", [str(k), str(l), lab])))
    n = len(topl)
    prof.append((n + 1, cls("SPRO", [str(n + 1), str(hs)])))
    groups = []
    for lab in ("Sand", "Rock"):
        ks = [k for k, x in enumerate(labels, start=1) if x == lab]
        if ks:
            groups.append(f"{lab} {_ranges(ks)}")
    unit = "m/s" if us == "SI" else "ft/s"
    vs = fmt_num(round(thr, 1))
    text = (f"SPRO blank: the default SOIL profile is used -- sublayers 1-{n} = the TOPL layers "
            f"{_short_list(topl)} with the built-in curves Sand (Vs < {vs} {unit}) and Rock (Vs >= {vs} {unit}): "
            f"{'; '.join(groups)}; sublayer {n + 1} = the half-space L {hs} (linear)")
    return prof, DefaultUse("SPRO", f"{n + 1} sublayers", text), ""


def library_curves(model, labels_used: Sequence[str]) -> Tuple[Dict[str, List[LIB.DynpPoint]], List[DefaultUse]]:
    """Library curves for the used DYNP labels the model does not define (D-W5-09), and the defaults used.
    A model DYNP of the same label always wins."""
    out: Dict[str, List[LIB.DynpPoint]] = {}
    uses: List[DefaultUse] = []
    if not enabled(model):
        return out, uses
    have = set(model_dynp_labels(model))
    for lab in dict.fromkeys(labels_used):
        if not lab or lab in have:
            continue
        pts = LIB.dynp_points(lab)
        if not pts:
            continue
        out[lab] = pts
        uses.append(DefaultUse(f"DYNP {lab}", lab,
                               f"DYNP {lab} not defined in the model: the built-in curve {lab} "
                               f"({LIB.DYNP_LABELS.get(lab, 'SHAKE91')}; @{LIB.DYNP_FILE}) is used"))
    return out, uses


def dynp_table(model, labels_used: Sequence[str]) -> List[Tuple[str, int, float, float, float, float, bool, bool]]:
    """The DYNP points AFWRITE writes to the SOIL deck, in curve (FILE73) order: the used labels in profile
    order, then the other model labels; ``(label, no, sg, g, sd, d, has_g, has_d)``.  Labels the model does
    not define come from the library (D-W5-09)."""
    dyn = get_entries(model, "DYNP")
    used = list(dict.fromkeys(l for l in labels_used if l))
    labels = used + [lab for (lab, _), _ in dyn if lab not in used]
    lib, _ = library_curves(model, used)
    rows = []
    for lab in dict.fromkeys(labels):
        pts = sorted(((no, rec) for (l2, no), rec in dyn if l2 == lab), key=lambda t: t[0])
        if pts:
            for no, rec in pts:
                rows.append((lab, int(no), float(rec.sg), float(rec.g), float(rec.sd), float(rec.d),
                             rec.given(2) and rec.given(3), rec.given(4) and rec.given(5)))
        elif lab in lib:
            for pt in lib[lab]:
                rows.append((lab, pt.no, float(pt.sg), float(pt.g), float(pt.sd), float(pt.d), True, True))
    return rows


# ======================================================================================
# helpers
# ======================================================================================
def _ranges(ks: Sequence[int]) -> str:
    out, start, prev = [], None, None
    for k in list(ks) + [None]:
        if start is None:
            start = prev = k
            continue
        if k is not None and k == prev + 1:
            prev = k
            continue
        out.append(f"{start}" if start == prev else f"{start}-{prev}")
        start = prev = k
    return ", ".join(out)


def _short_list(vals: Sequence[int], n: int = 8) -> str:
    vals = list(vals)
    s = ", ".join(str(v) for v in vals[:n])
    return s + (" ..." if len(vals) > n else "")


__all__ = ["SEISMIC_RECORD", "VIBRATION_LOAD", "TARGET_SPECTRUM", "ROCK_VS", "EDU_DEFAULT", "DefaultUse",
           "HistoryChoice", "enabled", "control_history", "soil_history", "soil_header_note", "equake_files",
           "default_output_name", "model_dynp_labels", "soil_unit_system", "soil_profile", "library_curves",
           "dynp_table"]
