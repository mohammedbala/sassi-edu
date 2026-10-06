"""LIBRARY: the built-in input library and the defaults of blank inputs (requirements section 7.19,
D-W5-01 ... D-W5-14; :mod:`sassi.io.library`, :mod:`sassi.prep.defaults`).

``LIBRARY``                lists every built-in input file (``@`` name, kind, description);
``LIBRARY,<kind>``         only one kind: RECORD, LOAD, SPECTRUM, PSD or DYNP;
``LIBRARY,<@name>``        one file: description, units, source and its path in this installation;
``LIBRARY,DEFAULTS``       the defaults of blank inputs the active model uses, per AOPT-enabled module
                           (what CHECK reports as Warning EDU-29), or why a blank input has none.
"""
from __future__ import annotations

from ...io import library as LIB
from ..registry import command
from .checks import search_dirs

_KIND_WORDS = {"RECORD": "record", "RECORDS": "record", "ACC": "record", "LOAD": "load", "LOADS": "load",
               "SPECTRUM": "spectrum", "SPECTRA": "spectrum", "RSIN": "spectrum", "PSD": "psd", "TPSD": "psd",
               "DYNP": "dynp", "CURVES": "dynp", "SOIL": "dynp"}


@command("LIBRARY", max_args=1)
def cmd_library(c):
    """LIBRARY,[kind|@name|DEFAULTS]: list the built-in input library (files named @<file>, sassi/data/library),
    one kind (RECORD, LOAD, SPECTRUM, PSD, DYNP) or one file, or the defaults of blank inputs the active model
    uses (DEFAULTS)."""
    arg = c.str(1).strip()
    key = arg.upper()
    if key == "DEFAULTS":
        _defaults(c)
        return
    if arg and (LIB.is_library_name(arg) or (LIB.entry(arg) is not None and key not in _KIND_WORDS)):
        e = LIB.entry(arg)
        if e is None:
            c.fail(f"{arg}: not a built-in input ({', '.join(LIB.names())})")
        c.info(f"LIBRARY {e.ref}: {e.title}")
        c.info(f"  kind      {e.kind}: {LIB.KINDS[e.kind]}")
        c.info(f"  units     {e.units}")
        c.info(f"  source    {e.source}")
        if e.kind == "dynp":
            for lab in LIB.dynp_labels():
                c.info(f"  curve     {lab}: {LIB.DYNP_LABELS.get(lab, '')} ({len(LIB.dynp_points(lab))} points)")
        c.info(f"  file      {e.path}{'' if e.path.is_file() else '  (missing in this installation)'}")
        return
    kind = None
    if key:
        kind = _KIND_WORDS.get(key)
        if kind is None:
            c.fail(f"'{arg}': give RECORD, LOAD, SPECTRUM, PSD, DYNP, DEFAULTS or a built-in file name (@name)")
    rows = LIB.entries(kind)
    c.info(f"LIBRARY: built-in inputs{'' if kind is None else ' (' + LIB.KINDS[kind] + ')'}; use them as @<file>, "
           f"e.g. THFILE,{LIB.CATALOGUE[0].ref}")
    for e in rows:
        c.info(f"  {e.ref:<28s} {e.kind:<9s} {e.title}")
    if kind in (None, "dynp"):
        c.info(f"  DYNP labels {', '.join(LIB.dynp_labels())} resolve without INP when the model does not "
               f"define them (SOIL)")
    c.info(f"  folder {LIB.LIBRARY_DIR}; LIBRARY,<@name> describes one file, LIBRARY,DEFAULTS the defaults of "
           f"blank inputs of the active model (EDUOPT,DEFAULTS,OFF switches them off)")


def _defaults(c) -> None:
    """LIBRARY,DEFAULTS: the defaults the active model uses, module by module (AOPT order)."""
    from ..check import Resolved
    from ..defaults import enabled
    from ..options import get_record
    m = c.model
    if not enabled(m):
        c.info("LIBRARY,DEFAULTS: EDUOPT,DEFAULTS,OFF -- no defaults: blank inputs are errors (ACS behaviour)")
        return
    r = Resolved(m, search_dirs(c))
    mods = [x for x in get_record(m, "AOPT").enabled() if x in ("EQUAKE", "SOIL", "MOTION", "STRESS", "RELDISP")]
    c.info(f"LIBRARY,DEFAULTS: built-in defaults of blank inputs for model {m.name or '(unnamed)'} "
           f"(AOPT modules {', '.join(mods) or 'none of EQUAKE, SOIL, MOTION, STRESS, RELDISP'})")
    n = 0
    for mod in mods:
        uses = r.defaults_used(mod)
        for u in uses:
            c.info(f"  {mod}: {u.text}")
            n += 1
        if mod in ("MOTION", "STRESS", "RELDISP") and r.history_needed(mod) and not r.thfile \
                and r.history_choice.why_not:
            c.info(f"  {mod}: THFILE blank, no default: {r.history_choice.why_not} (CHECK Error 73)")
        if mod == "SOIL":
            hc = r.soil_history()
            if not hc.name and hc.why_not:
                c.info(f"  SOIL: THFILE blank, no default: {hc.why_not} (CHECK Error 73)")
            prof, _, why = r.soil_profile()
            if not prof and why:
                c.info(f"  SOIL: {why} (CHECK Error 95)")
    if not n:
        c.info("  none: every input these modules read is given")
