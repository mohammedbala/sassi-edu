"""CHECK and the model-checking commands (requirements sections 3.4.A and 3.4.H; manual sections 9.2.7
and 9.8; spec 09 section 3; D-CHK-01..09).

* CHECK lists the errors and warnings of the model and of every AOPT-enabled module (the Check
  Errors window) and writes ``<model>.err`` in the model directory (requirements section 5.5).
* EXCSTRCHK, FIXEDINT, FREESPRING, HINGED, INTCOUNT and KINT are read-only reports; USED fixes the
  unused nodes.  Listings are capped by the Options > Check break number (spec 09 section 3).

The Options > Check settings (Show Warnings, Show Errors, Suppress Error Window, Break Check at)
are session settings, not saved (UI-07): ``interp.session['check_options']``.  The last report is
kept in ``interp.session['check_report']`` for the GUI (View > Check Errors).
"""
from __future__ import annotations

from pathlib import Path
from typing import List

from ..check import (CheckOptions, ModelView, analys_memory_bytes, dense_matrix_bytes, excstrchk, excstrchk_severity,
                     fixed_interaction, free_springs, hinges, kint, physical_ram, run_check, unused_nodes,
                     write_err)
from ..options import eduopt_float
from ..registry import command


def check_options(interp) -> CheckOptions:
    """The session's Options > Check settings (created with the defaults on first use)."""
    opt = interp.session.get("check_options")
    if not isinstance(opt, CheckOptions):
        opt = CheckOptions()
        interp.session["check_options"] = opt
    return opt


def search_dirs(c) -> List[Path]:
    """Directories for relative input files: model directory, then the working directory (L15)."""
    out = []
    if c.model.path:
        out.append(Path(c.model.path))
    out.append(c.interp.cwd)
    return out


def _listing(c, title: str, lines: List[str], empty: str) -> None:
    """Print a report capped by the Check break number (totals always printed)."""
    opt = check_options(c.interp)
    c.info(title)
    if not lines:
        c.info(f"  {empty}")
        return
    cap = opt.break_at if opt.break_at and opt.break_at > 0 else len(lines)
    for ln in lines[:cap]:
        c.info("  " + ln)
    if len(lines) > cap:
        c.info(f"  ... {len(lines) - cap} more not shown (Break Check at {cap})")


@command("CHECK", max_args=0)
def cmd_check(c):
    """CHECK: check the model and the AOPT-enabled modules; write <model>.err (manual Chapter 10)."""
    m = c.model
    opt = check_options(c.interp)
    rep = run_check(m, dirs=search_dirs(c))
    c.interp.session["check_report"] = rep
    if m.name and m.path and Path(m.path).is_dir():
        p = write_err(rep, Path(m.path) / f"{m.name}.err", opt)
        c.interp.session["check_err_file"] = str(p)
    else:
        p = None
    if not opt.suppress_window:
        for ln in rep.format(opt).rstrip("\n").splitlines():
            c.info(ln)
    summary = rep.summary() + (f"; written to {p}" if p else " (.err not written: no MDL)")
    if rep.errors():
        c.warn(summary)
    else:
        c.confirm(summary)
        if opt.suppress_window:
            c.info(summary)


@command("EXCSTRCHK", max_args=0)
def cmd_excstrchk(c):
    """EXCSTRCHK: excavation interior nodes shared with structure, beam, spring or GM elements (G-15)."""
    v = ModelView(c.model)
    res = excstrchk(v)
    implicit = sum(1 for r in v.elems if r.type in (1, 4) and r.elem.etype == 0)
    if implicit:
        c.info(f"EXCSTRCHK: ETYPE 0 of {implicit} SOLID/PLANE elements resolved by the ground elevation "
               f"(as ETYPEGEN,0) for this check")
    iset, exc_nodes = set(v.interaction_nodes()), set(v.excavation()["nodes"])
    kinds = {n: excstrchk_severity(v, n, iset, exc_nodes) for n in res}
    lines = [f"node {n}: {'; '.join(own)} [CHECK {kinds[n].lower()}]" for n, own in sorted(res.items())]
    _listing(c, f"EXCSTRCHK: {len(res)} excavation interior nodes shared with the structure", lines,
             "none: structure and excavation share only interface nodes")
    if res:
        nerr = sum(1 for k in kinds.values() if k == "Error")
        c.warn(f"{len(res)} shared excavation interior nodes may cause incorrect SSI results ({nerr} are CHECK "
               f"errors; warnings: interaction nodes used only by SOLID/PLANE elements on the excavation mesh)")


@command("FIXEDINT", max_args=0)
def cmd_fixedint(c):
    """FIXEDINT: interaction nodes with any fixed translational DOF."""
    v = ModelView(c.model)
    res = fixed_interaction(v)
    lines = [f"node {n}: {' '.join(labs)} fixed" + (" (Error 124)" if len(labs) == 3 else "") for n, labs in res]
    _listing(c, f"FIXEDINT: {len(res)} fixed interaction nodes", lines, "no interaction node has a fixed translation")
    if res:
        c.warn(f"{len(res)} interaction nodes have fixed translations: correct the model (INT or D)")


@command("FREESPRING", max_args=0)
def cmd_freespring(c):
    """FREESPRING: unconstrained nodes connected only to springs."""
    v = ModelView(c.model)
    res = free_springs(v)
    lines = [f"node {n}" + (" (carries mass)" if mass else ": no mass -- singular stiffness") for n, mass in res]
    _listing(c, f"FREESPRING: {len(res)} free spring nodes", lines, "none")
    if any(not mass for _, mass in res):
        c.warn(f"{sum(1 for _, mass in res if not mass)} spring-only nodes without fixity or mass")


@command("HINGED", max_args=0)
def cmd_hinged(c):
    """HINGED: possible unintended hinges (6-DOF elements meeting SOLID/PLANE at one node; drilling joints)."""
    v = ModelView(c.model)
    res = hinges(v)
    lines = [f"node {n}: {txt}" for n, txt in res]
    _listing(c, f"HINGED: {len(res)} possible hinges", lines, "none found")
    if res:
        c.warn(f"{len(res)} possible unintended hinges: penetrate the solids with massless beams or shells")


@command("INTCOUNT", max_args=0)
def cmd_intcount(c):
    """INTCOUNT: number of interaction nodes and the ANALYS memory estimate (D-ANL-10, UT-18)."""
    v = ModelView(c.model)
    n = len(v.interaction_nodes())
    c.info(f"INTCOUNT: {n} interaction nodes ({3 * n} interaction DOFs)")
    if n:
        one = dense_matrix_bytes(n)
        tot = analys_memory_bytes(n)
        c.info(f"  one dense complex matrix (3N)^2 x 16 B = {one / 1e9:.4g} GB; ANALYS estimate 3 (3N)^2 x 16 B = "
               f"{tot / 1e9:.4g} GB (D-ANL-10)")
        ram = physical_ram()
        if ram:
            frac = eduopt_float(c.model, "MEMLIMIT", 0.8)          # same guard as CHECK EDU-20 (D-ANL-10)
            c.info(f"  physical memory {ram / 1e9:.4g} GB")
            if tot > frac * ram:
                c.warn(f"the ANALYS memory estimate exceeds {100 * frac:.4g} % of the physical memory "
                       f"(EDUOPT,MEMLIMIT = {frac:.4g}; EDU-20)")
        if n > 20000:
            c.warn("more than 20,000 interaction nodes (practical limit; use FI or FFV)")
    c.confirm(f"{n} interaction nodes")


@command("KINT", max_args=0)
def cmd_kint(c):
    """KINT: beam K nodes that are interaction nodes."""
    v = ModelView(c.model)
    res = kint(v)
    lines = [f"node {n}" + (" (used only as a K node)" if only else "") for n, only in res]
    _listing(c, f"KINT: {len(res)} beam K nodes are interaction nodes", lines, "none")
    if res:
        c.warn(f"{len(res)} K nodes are interaction nodes: results may be incorrect")


@command("USED", max_args=0)
def cmd_used(c):
    """USED: fix all DOFs of the nodes not used by any element (D,n,n,1,1,ALL); unused interaction
    nodes are reported, not fixed (fixing them would trip FIXEDINT)."""
    m = c.model
    v = ModelView(m)
    nodes = unused_nodes(v)
    fixed, inter = [], []
    for n in nodes:
        if 0 in m.nodes[n].flags:
            inter.append(n)
            continue
        m.nodes[n].fix = [1] * 6
        fixed.append(n)
    if inter:
        c.warn(f"{len(inter)} unused interaction nodes not fixed: {inter[:20]}")
    c.confirm(f"USED: {len(fixed)} unused nodes fixed")
