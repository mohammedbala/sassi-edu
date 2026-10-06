"""Non-linear (near-field) soil SSI commands of SASSI-EDU (requirements 4.4 item 7, 4.10 item 5, 2.4,
3.4.R; decisions D-NLS-01 ... D-NLS-04; :mod:`sassi.core.nlsoil`).

The iteration loop of the manual (spec 03 section 4.5) is::

    initiation  SITE -> POINT -> HOUSE (.pin; iteration 0) -> ANALYS (<mode> 0, <save> 1)
                -> STRESS (<iter> 1, per input direction) -> COMBXYZSTRAIN
    iteration   HOUSE (.liq = 1: FILE74 + FILE73 -> new properties, FILE78) -> ANALYS (<mode> 1,
                "New Structure", reuses COOXqqq) -> STRESS (per direction) -> COMBXYZSTRAIN -> converged?

Commands
--------
``.pin`` content (the HOUSE dialog "Non-Linear SSI / Input Data" file; AFWRITE writes
``<model>.pin`` from these records when HOUSEX <nlssi> = 1 and at least one PINGRP is defined --
without PINGRP a hand-edited ``.pin`` in the model directory is used as it is):

* ``PIN,<esf>,[<ncurv>]`` -- line 1 of the .pin: effective strain factor ESF (default: the SOIL
  ``<ratio>``, 0.65) and the number of curves NCURV (blank: the number of DYNP labels, i.e. the
  curves SOIL writes to FILE73).  NGRP is counted by AFWRITE.
* ``PINGRP,<igrp>,<istr>,[<gfac>],[<dfac>],[<icurve>]`` -- a nonlinear soil group: SOLID (3D) or
  PLANE (2D) elements of the *structure* model (ETYPE 1 -- an element left at ETYPE 0 below grade is
  excavated soil, D-HOU-01, and is rejected by CHECK, EDU-41) in the excavated volume; ISTR 0 maximum
  shear-strain component, 1 octahedral (3D) / maximum shear strain (2D); optional default factors
  and curve of all its materials.  NMAT and NELEM are counted by AFWRITE.
* ``PINMAT,<igrp>,<mat>,<gfac>,<dfac>,<icurve>`` -- the line of material ``<mat>`` (M-table number)
  of group ``<igrp>``: initial shear-modulus factor GFAC and initial damping factor DFAC, both
  relative to the FREE FIELD at the element (requirements 4.4 item 7; manual "1.00 indicates same
  shear modulus as in free-field"): iteration 0 has ``G = GFAC G_layer`` and ``beta = DFAC beta_layer``
  of the TOPL layer containing the element's mid-depth (its L properties, as for the excavated soil);
  and the FILE73 curve ICURVE, given as a curve number or as a DYNP label (AFWRITE converts a label
  to the number SOIL gives it in FILE73: the labels used by SPRO in profile order, then the other
  DYNP labels with at least one point).  Blank fields take the PINGRP defaults.
* ``PINDEL,[<igrp>]`` -- delete group ``<igrp>`` (its PINGRP and PINMAT records), all .pin records
  when blank.
* ``PINLIST`` -- list the .pin records and the .pin text AFWRITE would write.

What the near-field M-table material must hold: the LOW-STRAIN properties of the near-field soil --
density, Poisson's ratio (Vp/Vs) and G_max = rho Vs^2, the reference of the G/G_max curve in every
later iteration (``G = G_max (G/G_max)(gamma_eff)``).  Do not give it strain-compatible values (the
curve would soften them a second time); GFAC/DFAC only choose the starting point.  The curves the
.pin refers to must have G/G_max and damping points (CHECK EDU-41, as SOIL's Errors 97-99).

Run control:

* ``NLSSIRESET`` -- start a new nonlinear analysis: ``<model>.liq`` = 0 (the next HOUSE run takes
  the .pin initial properties), a new convergence history, and FILE78 / FILE74 of the previous
  analysis renamed ``FILE78.prev`` / ``FILE74.prev`` (no command or module uses them any more).
* ``COMBXYZSTRAIN,<FILE74x>,<FILE74y>,<FILE74z>,<out>`` -- COMB_XYZ_STRAIN (D-NLS-04): SRSS of the
  effective strains of the directional FILE74 files, element by element, ``out`` (default FILE74)
  with G and beta of the FILE73 curves at the combined strain.  Blank inputs are skipped (two
  directions, or one for the plain convergence check); file names are relative to the model
  directory.  The inputs must be of one iteration and one HOUSE run (an error otherwise).  Prints
  the convergence measure against FILE78 (D-NLS-02).
* ``NLSSIITER,<var>[+<var2>...],[<maxit>],[<tolG>],[<tolB>]`` -- iteration driver: runs the command
  lines stored in variable ``<var>`` (one command per item; quote items that contain commas; several
  variables joined with ``+`` run one after the other) once per iteration, at most ``<maxit>`` times
  (default 8, D-NLS-02), and stops as soon as the properties used in the last iteration (FILE78) and
  the strain-compatible properties of its response (FILE74) differ by less than ``<tolG>`` % in G
  (default 2) and ``<tolB>`` percentage points in damping (default 0.5).  Inside the body ``#`` is
  the pass number.  When the analysis in progress has already converged (``.liq`` = 1, FILE78 of the
  FILE4 in the directory, HOUSE deck and .pin unchanged since: :func:`sassi.core.nlsoil.current_status`)
  nothing is run; after NLSSIRESET or a change of the model the body always runs.  ``NLSSIITER``
  without a variable prints the convergence history.  Example (one input direction; ANALYS <mode> 1
  and the decks written before)::

      VAR,NLSTEP,RUNHOUSE,RUNANALYS,RUNSTRESS
      NLSSIITER,NLSTEP,8

  The same loop without the driver (always ``n`` passes, no convergence stop) is the plain FOREACH
  pattern of D-NLS-01: ``VAR,IT,1,2,...,n`` and ``FOREACH,IT,FOREACH,NLSTEP,@NLSTEP[#]``.

All the records survive WRITE -> INP (they are written under "Other stored options").
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ...core import nlsoil as NL
from ...model.options import Field, Record, register_record_type
from ...model.values import NumberError, parse_float, parse_int
from ..interpreter import LoopFrame
from ..options import get_entries, get_record
from ..registry import command
from .extensions import model_dir

TIER = "P1"


# ======================================================================================
# Records
# ======================================================================================
@register_record_type
class PinRecord(Record):
    """``PIN,<esf>,[<ncurv>]``: line 1 of the .pin (NGRP is counted by AFWRITE)."""
    COMMAND = "PIN"
    FIELDS = (Field("esf", float, 0.0, "effective strain factor (blank: SOIL <ratio>)"),
              Field("ncurv", int, 0, "number of FILE73 curves (blank: number of DYNP labels)"))


@register_record_type
class PingrpRecord(Record):
    """``PINGRP,<igrp>,<istr>,[<gfac>],[<dfac>],[<icurve>]``: a nonlinear soil group."""
    COMMAND = "PINGRP"
    FIELDS = (Field("igrp", int, 0, "group number"), Field("istr", int, 0, "strain flag 0/1"),
              Field("gfac", float, 1.0, "default initial G factor"), Field("dfac", float, 1.0, "default damping factor"),
              Field("icurve", str, "", "default curve number or DYNP label"))


@register_record_type
class PinmatRecord(Record):
    """``PINMAT,<igrp>,<mat>,<gfac>,<dfac>,<icurve>``: one material line of a nonlinear group."""
    COMMAND = "PINMAT"
    FIELDS = (Field("igrp", int, 0, "group number"), Field("mat", int, 0, "material number"),
              Field("gfac", float, 1.0, "initial G factor"), Field("dfac", float, 1.0, "initial damping factor"),
              Field("icurve", str, "", "curve number or DYNP label"))


def _check_numbers(c, toks: Sequence[str], fields: Sequence[Field]) -> None:
    for i, (t, f) in enumerate(zip(toks, fields), start=1):
        t = t.strip()
        if not t or f.type is str:
            continue
        try:
            v = parse_float(t)
        except NumberError:
            c.fail(f"argument {i} <{f.name}>: '{t}' is not a number")
        if f.type is int and v != round(v):
            c.fail(f"argument {i} <{f.name}>: {t} must be an integer")


@command("PIN", tier=TIER, max_args=2)
def cmd_pin(c):
    """PIN,<esf>,[<ncurv>]: .pin line 1 -- effective strain factor ESF (blank: SOIL <ratio>) and number of
    FILE73 curves NCURV (blank: number of DYNP labels); NGRP is counted by AFWRITE.  The groups follow with
    PINGRP / PINMAT (GFAC, DFAC relative to the free field; the materials hold the low-strain soil)."""
    toks = list(c.tokens[:2])
    _check_numbers(c, toks, PinRecord.FIELDS)
    rec = PinRecord("PIN", toks)
    esf = rec.number(1)
    if esf is not None and not 0.0 < esf <= 1.0:
        c.fail(f"ESF = {esf:g} must lie in (0, 1] (typically 0.6-0.7)")
    if rec.given(2) and rec.integer(2) < 1:
        c.fail("NCURV must be >= 1 (blank = number of DYNP labels)")
    c.model.options.set_record(rec)
    c.confirm(f"PIN: ESF {esf if esf is not None else '= SOIL <ratio>'}, NCURV "
              f"{rec.integer(2) if rec.given(2) else '= number of DYNP labels'}")


@command("PINGRP", tier=TIER, max_args=5)
def cmd_pingrp(c):
    """PINGRP,<igrp>,<istr>,[<gfac>],[<dfac>],[<icurve>]: nonlinear soil group (SOLID/PLANE structure
    elements, ETYPE 1, whose M-table materials hold the low-strain soil: G_max) with strain flag ISTR
    (0 maximum shear component, 1 octahedral / 2D maximum shear) and optional default GFAC, DFAC (initial
    factors on the free-field G and damping of the TOPL layer at each element) and ICURVE (number or DYNP
    label) of its materials."""
    toks = list(c.tokens[:5])
    _check_numbers(c, toks, PingrpRecord.FIELDS)
    rec = PingrpRecord("PINGRP", toks)
    g = rec.integer(1, 0)
    if g <= 0:
        c.fail("group number <igrp> >= 1 required")
    if rec.integer(2, 0) not in NL.STRAIN_FLAGS:
        c.fail("ISTR must be 0 (maximum shear-strain component) or 1 (octahedral / 2D maximum shear strain)")
    _factor_checks(c, rec.number(3), rec.number(4))
    c.model.options.set_entry("PINGRP", g, rec)
    grp = c.model.groups.get(g)
    if grp is None:
        c.warn(f"group {g} is not defined (yet)")
    elif int(grp.type) not in NL.NL_TYPES:
        c.warn(f"group {g} is not a SOLID or PLANE group: HOUSE will stop")
    c.confirm(f"PINGRP {g}: ISTR {rec.integer(2, 0)}")


@command("PINMAT", tier=TIER, max_args=5)
def cmd_pinmat(c):
    """PINMAT,<igrp>,<mat>,<gfac>,<dfac>,<icurve>: .pin line of material <mat> of nonlinear group <igrp>:
    iteration 0 has G = GFAC G_layer and beta = DFAC beta_layer of the free-field TOPL layer at the
    element (1.0 = the free field), later iterations G = G_max (G/G_max) with G_max of material <mat> (the
    low-strain soil); blank fields: the PINGRP defaults; ICURVE a FILE73 curve number or a DYNP label."""
    toks = list(c.tokens[:5])
    _check_numbers(c, toks, PinmatRecord.FIELDS)
    rec = PinmatRecord("PINMAT", toks)
    g, mid = rec.integer(1, 0), rec.integer(2, 0)
    if g <= 0 or mid <= 0:
        c.fail("group <igrp> and material <mat> numbers >= 1 required")
    _factor_checks(c, rec.number(3), rec.number(4))
    if c.model.options.entry("PINGRP", g) is None:
        c.warn(f"group {g} has no PINGRP record (yet): define it, otherwise the line is not written")
    c.model.options.set_entry("PINMAT", (g, mid), rec)
    c.confirm(f"PINMAT group {g} material {mid} set")


def _factor_checks(c, gfac: Optional[float], dfac: Optional[float]) -> None:
    if gfac is not None and gfac <= 0:
        c.fail(f"GFAC = {gfac:g} must be > 0")
    if dfac is not None and dfac < 0:
        c.fail(f"DFAC = {dfac:g} must be >= 0")


@command("PINDEL", tier=TIER, max_args=1)
def cmd_pindel(c):
    """PINDEL,[<igrp>]: delete nonlinear group <igrp> (PINGRP and its PINMAT records); blank = all .pin records."""
    opts = c.model.options
    if not c.given(1):
        n = len(opts.entries("PINGRP")) + len(opts.entries("PINMAT")) + (1 if opts.record("PIN") else 0)
        for k, _ in opts.entries("PINGRP"):
            opts.delete_entry("PINGRP", k)
        for k, _ in opts.entries("PINMAT"):
            opts.delete_entry("PINMAT", k)
        opts.delete_record("PIN")
        c.confirm(f"PINDEL: {n} .pin records deleted")
        return
    g = c.int(1)
    found = opts.delete_entry("PINGRP", g)
    for k, _ in opts.entries("PINMAT"):
        if isinstance(k, tuple) and int(k[0]) == g:
            opts.delete_entry("PINMAT", k)
            found = True
    c.confirm(f"PINDEL: group {g} deleted" if found else f"PINDEL: group {g} had no .pin records")


# ======================================================================================
# Model -> .pin (used by AFWRITE, CHECK and PINLIST)
# ======================================================================================
def _dynp_points(model, label: str) -> Tuple[List[Tuple[float, float]], List[Tuple[float, float]]]:
    """(G points, damping points) of a DYNP label, as AFWRITE writes them to the SOIL deck: the model's
    points, else the built-in curve of that label (D-W5-09; not with EDUOPT,DEFAULTS,OFF)."""
    from ..defaults import enabled
    from ...io import library as LIB
    g, dmp = [], []
    for (lab, _no), rec in get_entries(model, "DYNP"):
        if lab != label:
            continue
        if rec.given(2) and rec.given(3):
            g.append((float(rec.sg), float(rec.g)))
        if rec.given(4) and rec.given(5):
            dmp.append((float(rec.sd), float(rec.d)))
    if not g and not dmp and enabled(model) and label not in {lab for (lab, _), _ in get_entries(model, "DYNP")}:
        for pt in LIB.dynp_points(label):
            g.append((float(pt.sg), float(pt.g)))
            dmp.append((float(pt.sd), float(pt.d)))
    return g, dmp


def file73_labels(model) -> List[str]:
    """DYNP labels in FILE73 curve order: the labels used by SPRO in profile order, then the other DYNP
    labels in key order -- the order of the SOIL deck written by AFWRITE (sassi.prep.afwrite
    DeckBuilder.soil, :func:`sassi.prep.defaults.dynp_table`: the default profile and the built-in curves
    of labels the model does not define included), which SOIL keeps in FILE73 (curve number = position,
    1-based).  A label without any complete point is not in the deck (nor in FILE73) and is left out."""
    from ..defaults import dynp_table, soil_profile
    prof = soil_profile(model)[0]
    used = list(dict.fromkeys(rec.dynprop for _, rec in prof if rec.dynprop))
    rows = dynp_table(model, used)
    return list(dict.fromkeys(row[0] for row in rows if row[6] or row[7]))


def curve_label_problems(model, label: str) -> List[str]:
    """Why the DYNP curve ``label`` cannot be used by the property update (SOIL Errors 97-99 applied to
    a curve of the .pin, which SOIL does not check when no SPRO sublayer uses it)."""
    g, dmp = _dynp_points(model, label)
    out = []
    if not any(s > 0 for s, _ in g):
        out.append("no G/G_max point with a positive strain (cf. Error 97)")
    elif any(not 0.0 <= v <= 1.0 for s, v in g if s > 0):
        out.append("G/G_max values outside [0, 1] (cf. Error 98)")
    if not any(s > 0 for s, _ in dmp):
        out.append("no damping point with a positive strain (cf. Error 99)")
    elif any(v < 0 for s, v in dmp if s > 0):
        out.append("negative damping values")
    return out


def has_pin_data(model) -> bool:
    return bool(model.options.entries("PINGRP"))


def _curve_number(token: Optional[str], labels: Sequence[str]) -> Tuple[Optional[int], str]:
    """(curve number, problem text) of an ICURVE token: an integer, or a DYNP label."""
    if token is None or not str(token).strip():
        return None, "no ICURVE (give it on PINMAT or as the PINGRP default)"
    t = str(token).strip()
    try:
        v, exact = parse_int(t)
        if not exact or v < 1:
            return None, f"ICURVE {t} must be a positive integer or a DYNP label"
        return int(v), ""
    except NumberError:
        pass
    for k, lab in enumerate(labels, start=1):
        if lab == t:
            return k, ""
    for k, lab in enumerate(labels, start=1):
        if lab.upper() == t.upper():
            return k, ""
    return None, f"ICURVE label '{t}' is not a DYNP label ({', '.join(labels) or 'no DYNP data'})"


def _resolved_excavated(model, gid: int, elems: Dict[int, object], view=None) -> List[int]:
    """Elements of group ``gid`` that HOUSE treats as excavated soil: ETYPE 2, and ETYPE 0 SOLID/PLANE
    elements below grade (D-HOU-01, resolved like CHECK and AFWRITE, D-AFW-06).  ``view`` is a CHECK
    :class:`sassi.prep.check.ModelView` (built here when needed and not given)."""
    out = sorted(e for e, el in elems.items() if int(el.etype) == 2)
    if any(int(el.etype) == 0 for el in elems.values()):
        if view is None:
            from ..check import ModelView              # lazy: check.py imports this module lazily too
            view = ModelView(model)
        out += [r.id for r in view.elems if r.group == gid and int(r.elem.etype) == 0 and r.etype == 2]
    return sorted(set(out))


def pin_from_model(model, view=None) -> Tuple[Optional[NL.PinData], List[Tuple[int, str]], List[str],
                                             Dict[int, List[str]]]:
    """Build the .pin content from the PIN / PINGRP / PINMAT records and the model.

    Returns ``(pin, problems, notes, material_notes)``: ``problems`` are ``(group, text)`` pairs that
    make the .pin invalid (CHECK reports them as EDU-41, AFWRITE then writes no .pin); ``pin`` is None
    when there is no PINGRP record.  ``view`` (optional) is the CHECK model view used to resolve
    ETYPE 0 (D-HOU-01)."""
    groups = get_entries_raw(model, "PINGRP")
    if not groups:
        return None, [], [], {}
    labels = file73_labels(model)
    all_labels = list(dict.fromkeys(lab for (lab, _), _ in get_entries(model, "DYNP")))
    probs: List[Tuple[int, str]] = []
    notes: List[str] = []
    pinrec = model.options.record("PIN")
    esf = pinrec.number(1) if pinrec is not None else None
    if esf is None:
        esf = float(get_record(model, "SOIL").get("ratio"))
        notes.append(f"ESF = SOIL <ratio> = {esf:g} (PIN not given)")
    ncurv = pinrec.integer(2) if (pinrec is not None and pinrec.given(2)) else None
    if ncurv is None:
        ncurv = max(1, len(labels))
        if not labels:
            notes.append("no DYNP labels: NCURV = 1 (FILE73 must hold the curves)")
    mats = {k: r for k, r in get_entries_raw(model, "PINMAT") if isinstance(k, tuple)}
    pin = NL.PinData(esf=float(esf), ncurv=int(ncurv))
    mnotes: Dict[int, List[str]] = {}
    for g, grec in groups:
        g = int(g)
        grp = model.groups.get(g)
        if grp is None:
            probs.append((g, f"group {g} is not defined"))
            continue
        if int(grp.type) not in NL.NL_TYPES:
            probs.append((g, f"group {g} is not a SOLID or PLANE group"))
            continue
        elems = grp.elements
        if not elems:
            probs.append((g, f"group {g} has no elements"))
            continue
        exc = _resolved_excavated(model, g, elems, view)
        if exc:
            ids = ", ".join(str(e) for e in exc[:8]) + (" ..." if len(exc) > 8 else "")
            probs.append((g, f"group {g} has excavated soil elements ({ids}: ETYPE 2, or ETYPE 0 below grade, "
                             "D-HOU-01): the near-field soil of a nonlinear group is structure -- set ETYPE 1 "
                             "(ETYPE,<first>,<last>,<step>,1) -- with an M-table soil material, next to the excavated soil"))
            continue
        mids = sorted({int(el.mat) for el in elems.values()})
        istr = grec.integer(2, 0)
        pg = NL.PinGroup(igrp=g, nmat=len(mids), nelem=len(elems), istr=int(istr))
        mnotes[g] = []
        for mid in mids:
            mr = mats.get((g, mid))
            gfac = _first(mr.number(3) if mr else None, grec.number(3), 1.0)
            dfac = _first(mr.number(4) if mr else None, grec.number(4), 1.0)
            tok = (mr.arg(5) if (mr is not None and mr.given(5)) else None) or grec.arg(5)
            icurve, why = _curve_number(tok, labels)
            if icurve is None:
                empty = str(tok or "").strip()
                if empty and empty.upper() in {lab.upper() for lab in all_labels} - {lab.upper() for lab in labels}:
                    why = (f"DYNP label '{empty}' has no complete G/G_max or damping point (it is not written to "
                           "FILE73; cf. Errors 97, 99)")
                probs.append((g, f"group {g} material {mid}: {why}"))
                icurve = 1
            elif icurve > ncurv:
                probs.append((g, f"group {g} material {mid}: ICURVE {icurve} > NCURV {ncurv}"))
            else:
                lab = labels[icurve - 1] if icurve <= len(labels) else ""
                for t in (curve_label_problems(model, lab) if lab else []):
                    probs.append((g, f"group {g} material {mid}: curve {lab} (ICURVE {icurve}): {t}"))
            if mid not in model.materials:
                probs.append((g, f"group {g}: material {mid} is not defined (M table)"))
            pg.materials.append(NL.PinMaterial(float(gfac), float(dfac), int(icurve)))
            lab = labels[icurve - 1] if 1 <= icurve <= len(labels) else ""
            mnotes[g].append(f"material {mid}" + (f", curve {lab}" if lab else ""))
        for (gg, mid), _ in mats.items():
            if gg == g and mid not in mids:
                notes.append(f"PINMAT group {g} material {mid}: the group has no element of this material (ignored)")
        pin.groups.append(pg)
    for (gg, mid) in mats:
        if model.options.entry("PINGRP", gg) is None:
            notes.append(f"PINMAT group {gg} material {mid}: no PINGRP {gg} (ignored)")
    return pin, probs, notes, mnotes


def get_entries_raw(model, name: str):
    """Entries of an indexed record type registered here (typed through sassi.model.options)."""
    return list(model.options.entries(name))


def _first(*vals):
    for v in vals:
        if v is not None:
            return v
    return None


def pin_text(model, view=None) -> Tuple[Optional[str], List[Tuple[int, str]], List[str]]:
    """The .pin text AFWRITE writes (None without PINGRP records or with problems)."""
    pin, probs, notes, mnotes = pin_from_model(model, view=view)
    if pin is None or probs:
        return None, probs, notes
    head = [f"SASSI-EDU .pin of model '{model.name}': nonlinear soil SSI input (manual 6.5.3, spec 05b section 6)",
            "written by AFWRITE from the PIN / PINGRP / PINMAT commands (edit those, not this file)",
            "line 1: NGRP, ESF, NCURV;  per group: IGRP, NMAT, NELEM, ISTR;  then NMAT lines GFAC, DFAC, ICURVE",
            "(the material lines follow the group's material numbers in ascending order;",
            " GFAC, DFAC scale the free-field G and damping of the TOPL layer at each element, iteration 0;",
            " the materials hold the low-strain soil, G_max of the curves)"]
    return NL.format_pin(pin, head, mnotes), probs, notes


@command("PINLIST", tier=TIER, max_args=0)
def cmd_pinlist(c):
    """PINLIST: list the .pin records (PIN, PINGRP, PINMAT) and the .pin text AFWRITE writes."""
    m = c.model
    if not has_pin_data(m):
        c.info("no PINGRP records: AFWRITE writes no .pin (a hand-edited <model>.pin is used as it is)")
        return
    text, probs, notes = pin_text(m)
    for g, t in probs:
        c.warn(t)
    for t in notes:
        c.info(t)
    if text:
        for ln in text.rstrip("\n").splitlines():
            c.info(ln)


# ======================================================================================
# Run control
# ======================================================================================
def _model_name(c) -> str:
    if not c.model.name:
        c.fail("model name not defined -- use MDL,<Model>,<Path>")
    return c.model.name


@command("NLSSIRESET", tier=TIER, max_args=0)
def cmd_nlssireset(c):
    """NLSSIRESET: start a new nonlinear soil SSI analysis -- <model>.liq = 0 (the next HOUSE run takes the
    .pin initial properties), a new NLSOIL_CONVERGENCE.TXT history, and FILE78 / FILE74 of the previous
    analysis renamed FILE78.prev / FILE74.prev (STRESS, COMBXYZSTRAIN and NLSSIITER no longer see them)."""
    name = _model_name(c)
    d = model_dir(c)
    if not d.is_dir():
        c.fail(f"model directory {d} does not exist")
    NL.write_liq(d / f"{name}.liq", 0)
    NL.reset_convergence(d)
    old = NL.retire(d, ("FILE78", "FILE74"))
    c.confirm(f"{name}.liq = 0: the next HOUSE run starts the nonlinear soil iterations from the .pin (iteration 0)"
              + (f"; the files of the previous analysis are kept as {', '.join(old)}" if old else ""))


def _in_dir(c, name: str) -> Path:
    p = Path(name.strip().replace("\\", "/"))
    return p if p.is_absolute() else model_dir(c) / p


@command("COMBXYZSTRAIN", max_args=4)
def cmd_combxyzstrain(c):
    """COMBXYZSTRAIN,<FILE74x>,<FILE74y>,<FILE74z>,<out>: COMB_XYZ_STRAIN (D-NLS-04) -- SRSS of the directional
    effective strains per element -> <out> (default FILE74), with G and beta from FILE73; blank inputs skipped."""
    names = [c.str(k) for k in (1, 2, 3)]
    out_name = c.str(4) or "FILE74"
    given = [(k, n) for k, n in enumerate(names) if n]
    if not given:
        c.fail("at least one directional FILE74 is required")
    files, labels = [], []
    for k, n in given:
        p = _in_dir(c, n)
        try:
            files.append(NL.read_nlfile(p, "FILE74"))
        except NL.NLSoilError as exc:
            c.fail(str(exc))
        labels.append(files[-1].direction or "XYZ"[k])
    d = model_dir(c)
    curves = None
    try:
        curves = NL.load_curves(d / "FILE73") if (d / "FILE73").is_file() else None
    except NL.NLSoilError as exc:
        c.warn(str(exc))
    try:
        comb, warns = NL.combine_files(files, curves, labels)
    except NL.NLSoilError as exc:
        c.fail(str(exc))
    for w in warns:
        c.warn(w)
    outp = _in_dir(c, out_name)
    NL.write_nlfile(outp, comb, title="COMB_XYZ_STRAIN")
    c.confirm(f"{outp.name}: SRSS of {len(files)} directional effective strains for {len(comb.rows)} elements "
              f"(iteration {comb.iteration})")
    p78 = d / "FILE78"
    if curves and p78.is_file():
        try:
            f78 = NL.read_nlfile(p78, "FILE78")
            if f78.iteration == comb.iteration and NL.same_file4(f78.file4, comb.file4):
                c.info("COMBXYZSTRAIN: " + NL.compare(f78.rows, comb.rows, comb.iteration).text())
        except NL.NLSoilError:
            pass


def _stamp(path: Path) -> Optional[Tuple[int, str]]:
    """(modification time ns, content digest) of a file, None when absent: tells whether a pass rewrote
    FILE74 (a new iteration always changes its ITERATION line; the time covers identical rewrites)."""
    try:
        return (path.stat().st_mtime_ns, hashlib.sha1(path.read_bytes()).hexdigest())
    except OSError:
        return None


def _history_lines(rows) -> List[str]:
    out = [f"{'iteration':>10s}{'max dG/G %':>13s}{'max dbeta %':>13s}{'mean G/Gmax':>13s}{'mean beta':>11s}  converged"]
    for h in rows:
        out.append(f"{int(h['iteration']):>10d}{h['max_dG_pct']:>13.3f}{h['max_dbeta_pct']:>13.3f}"
                   f"{h['mean_G_Gmax']:>13.4f}{h['mean_beta']:>11.4f}  {'yes' if h['converged'] else 'no'}")
    return out


@command("NLSSIITER", tier=TIER, max_args=4)
def cmd_nlssiiter(c):
    """NLSSIITER,<var>[+<var2>...],[<maxit>],[<tolG>],[<tolB>]: run the iteration commands stored in the
    variable(s) (one per item) until the nonlinear soil properties converge (max |dG/G| < tolG %, max
    |d beta| < tolB %, D-NLS-02) or <maxit> passes (default 8); '#' in the body = pass number.  Nothing
    runs when the analysis in progress has already converged (.liq = 1, FILE78 of the FILE4 in the
    directory, HOUSE deck and .pin unchanged); after NLSSIRESET or a model change the body always runs.
    Every pass must write a new FILE74 (else it stops: no convergence row from stale files).  Without
    <var>: print the convergence history."""
    d = model_dir(c)
    tol_g = c.float(3, default=NL.TOL_G_PCT)
    tol_b = c.float(4, default=NL.TOL_BETA_PCT)
    if tol_g <= 0 or tol_b <= 0:
        c.fail("the tolerances must be > 0")
    if not c.given(1):
        rows = NL.read_convergence(d)
        if not rows:
            c.info("no convergence history (NLSOIL_CONVERGENCE.TXT) in the model directory")
        for ln in _history_lines(rows):
            c.info(ln)
        st = NL.status(d, tol_g, tol_b)
        if st is not None:
            cur, why = NL.current_status(d, c.model.name, tol_g, tol_b) if c.model.name else (None, "no model name")
            c.info("last complete iteration: " + st.text()
                   + ("" if cur is not None else f" -- not the analysis in progress: {why}"))
        return
    body: List[str] = []
    names = [t.strip().lstrip("@") for t in c.str(1).split("+")]
    for name in names:
        v = c.interp.variables.get(name.lower()) if name else None
        if v is None:
            c.fail(f"variable {name or '(blank)'} not defined (VAR,{name or 'NAME'},<command 1>,<command 2>,...)")
        body += list(v.items)
    var = c.interp.variables[names[0].lower()]
    if not body:
        c.fail(f"variable(s) {'+'.join(names)} hold no commands")
    maxit = c.int(2, default=NL.MAX_ITERATIONS)
    if maxit < 1:
        c.fail("<maxit> must be >= 1")
    loops = c.interp.loops
    if any(fr.var in {n.lower() for n in names} for fr in loops):
        c.fail(f"NLSSIITER inside a loop over {var.name}")
    # "already converged" only for the analysis in progress (never after NLSSIRESET or a model / .pin
    # change), and without writing a history row: only passes that ran are recorded
    cur, why = NL.current_status(d, _model_name(c), tol_g, tol_b)
    if cur is not None and cur.converged:
        c.info(f"NLSSIITER: already converged -- {cur.text()} (nothing run; NLSSIRESET starts a new analysis)")
        return
    old = NL.status(d, tol_g, tol_b)
    if cur is None and old is not None and old.converged:
        c.info(f"NLSSIITER: the converged iteration {old.iteration} in the model directory is not that of the "
               f"analysis in progress ({why}): the iterations are run")
    frame = LoopFrame(var.name.lower(), 0)
    loops.append(frame)
    passes = 0
    conv = None
    try:
        for k in range(1, maxit + 1):
            frame.index = k
            stamp = _stamp(d / "FILE74")
            for item in body:
                if not c.interp.execute(item):
                    c.fail(f"pass {k}: '{item}' failed; iterations stopped (fix the problem, then run NLSSIITER again)")
            passes = k
            if stamp is not None and _stamp(d / "FILE74") == stamp:
                c.fail(f"pass {k}: FILE74 was not written -- the body must run HOUSE (non-linear SSI), ANALYS and "
                       "STRESS (<iter> = 1, and COMBXYZSTRAIN for X/Y/Z input); no convergence row recorded")
            conv = NL.status(d, tol_g, tol_b)
            if conv is None:
                c.fail(f"pass {k}: no FILE78/FILE74 of one iteration in {d} -- the body must run HOUSE (non-linear SSI), "
                       "ANALYS and STRESS (<iter> = 1, and COMBXYZSTRAIN for X/Y/Z input)")
            NL.record_convergence(d, conv)
            c.info(f"NLSSIITER pass {k}: {conv.text()}")
            if conv.converged:
                break
    finally:
        loops.pop()
    if conv is not None and conv.converged:
        c.info(f"NLSSIITER: converged after {passes} pass(es): the properties of iteration {conv.iteration} "
               "(FILE78) are strain-compatible with its response (FILE74)")
    else:
        c.warn(f"not converged after {passes} pass(es) (limit {maxit}, D-NLS-02): run NLSSIITER again or check the model")
    for ln in _history_lines(NL.read_convergence(d)):
        c.info(ln)
