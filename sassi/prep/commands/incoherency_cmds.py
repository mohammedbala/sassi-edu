"""Commands of incoherent seismic input, wave passage and multiple excitation (requirements 3.4.F,
3.4.R, 4.4 item 6; spec 07 sections 9.2.4, 9.2.16, 9.2.17, 9.2.20, 9.2.46; spec 05b sections 3.8-3.19
and 7).

The option records are those of :mod:`sassi.prep.options` (HOUSE, INCOH, WPASS, ME and the AMP list,
HOUSEX <supmode>/<nsim>, ANALYSX <ffm>); AFWRITE copies them into the HOUSE and ANALYS decks.  Since
wave 3 HOUSE and ANALYS implement them (:mod:`sassi.core.coherency`, :mod:`sassi.core.incoherency`),
so the handlers below replace the earlier ones that stored the records with a "tier P1 not available"
note; they store exactly the same records:

* ``HOUSE,<gravity>,<gelev>,<opmode>,<dim>,<imp>,<coh>,<wpass>,<me>,<cmplxspec>`` -- <coh> 1 incoherent
  motion, <wpass> 1 wave passage along Line D, <me> 1 multiple excitation, <cmplxspec> 1 complex
  spectral amplification ratios;
* ``INCOH,<gammax>,<gammay>,<gammaz>,<alpha>,<ngp>,<ipr>,<nmodes>,<met>,<HSeed>,<VSeed>,<RandPhz>`` --
  Luco-Wong coherence parameters per component, directionality factor alpha (models 2-7) or mean Vs
  (model 1), embedded levels, I N C O table, incoherent modes (0 all, k first k, -k mode k only),
  units flag, seeds and random phase angle of stochastic simulation (stochastic iff a seed is non-zero
  and RandPhz > 0, D-INC-02);
* ``WPASS,<appv>,<ang>,<cohf>`` -- apparent velocity and angle of Line D, unlagged coherency model
  1-7 (D-INC-01; models 2 and 4 refuse to run: their coefficients are not available, D-INC-04);
* ``ME,<no>,<nfirst>,<nlast>,<xc>,<yc>,<zc>`` -- multiple-excitation zone: interaction nodes
  ``nfirst..nlast`` (the control-point coordinates are stored, not used);
* ``AMP,<no>,<a1>,...`` (list command of :mod:`sassi.prep.commands.session`) -- spectral
  amplification ratios of zone <no>, one per SSI frequency ((Re, Im) pairs with <cmplxspec> = 1);
* ``HOUSEX,<optimize>,<supmode>,<nsim>,<nlssi>,<ansys>`` -- <supmode> 0 Linear (AS) / 1 Quadratic
  (SRSS TF, one run per mode with <nmodes> = -k), <nsim> stochastic simulations (1..50);
* ``ANALYSX,<ffm>,<delrst>`` -- <ffm> 0 free-field load (FFL, default), 1 free-field motion (FFM).

Action command:

* ``BUILDFILE77,<out>,<in1>,...,<inN>`` -- Build_FILE77 (spec 05b section 7, the per-level approach of
  deeply embedded foundations): combines the FILE77s of separate HOUSE incoherent runs on the
  interaction nodes of one embedment level each into ``<out>`` covering all of them.  Inputs must have
  the same frequency numbers and step and disjoint node sets; when ``<model>.N4`` (FILE4 of the final
  HOUSE run with all interaction nodes) is in the model directory the nodes are arranged in its
  interaction order and every interaction node must be covered.  File names are relative to the model
  directory.  Typical sequence: for every level, HOUSE (incoherent, that level's interaction nodes
  only) and ``FCOPY,FILE77,FILE77_L<k>``; then the final HOUSE run with all interaction nodes; then
  ``BUILDFILE77,FILE77,FILE77_L1,FILE77_L2,...``; then ANALYS (<coh> = 1).  The accuracy rests on
  consistent mode signs at all levels (sign convention D-INC-03).

ACS SASSI runs Build_FILE77 as a separate DOS program; SASSI-EDU provides it as this extension
command (requirements 3.4.R).
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

# The option handlers registered by these modules are replaced below: import them first so that the
# registry sees their registration before ours (modules are imported in alphabetical order).
from . import extensions as _extensions  # noqa: F401
from . import options as _options  # noqa: F401
from ...core import incoherency as INC
from ...io.container import read_container, write_container
from ..registry import command, indexed_key
from .extensions import store_xrecord
from .options import build_record, report_problems


# ======================================================================================
# Option records (the same storage as before, without the "tier P1" notes)
# ======================================================================================
def _store_record(c, name: str):
    rec = build_record(c, name, c.tokens)
    c.model.options.set_record(rec)
    report_problems(c, rec)
    c.confirm(f"{name} options set ({len(rec)} arguments)")
    return rec


@command("HOUSE", replaces=True)
def cmd_house(c):
    """HOUSE,<gravity>,<gelev>,<opmode>,<dim>,<imp>,<coh>,<wpass>,<me>,<cmplxspec>: HOUSE options
    (<coh> 1 incoherent motion, <wpass> 1 wave passage, <me> 1 multiple excitation)."""
    rec = _store_record(c, "HOUSE")
    if rec.dim not in (1, 2):
        c.warn("<dim> must be 1 (2D, POINT2) or 2 (3D, POINT3); 1D is not available")
    if rec.coh and rec.dim != 2:
        c.warn("incoherent motion needs a 3D model (<dim> = 2, no SYMM planes; G-17) -- CHECK will report it")
    if rec.me and not rec.wpass:
        c.warn("multiple excitation needs wave passage (<wpass> = 1) -- CHECK will report it")


@command("INCOH", replaces=True)
def cmd_incoh(c):
    """INCOH,<gammax>,<gammay>,<gammaz>,<alpha>,<ngp>,<ipr>,<nmodes>,<met>,<HSeed>,<VSeed>,<RandPhz>:
    incoherency options of HOUSE (requirements 4.4 item 6)."""
    rec = _store_record(c, "INCOH")
    if (int(rec.hseed) != 0 or int(rec.vseed) != 0) and float(rec.randphz) <= 0:
        c.info("INCOH: seeds without a random phase angle give the deterministic input (RandPhz = 180 for "
               "stochastic simulation, D-INC-02)")


@command("WPASS", replaces=True)
def cmd_wpass(c):
    """WPASS,<appv>,<ang>,<cohf>: apparent velocity and angle of Line D; <cohf> = unlagged coherency model
    1-7 (D-INC-01)."""
    from ...core.coherency import model_available, unavailable_reason
    rec = _store_record(c, "WPASS")
    if 2 <= int(rec.cohf) <= 6 and not model_available(int(rec.cohf)):
        c.warn(unavailable_reason(int(rec.cohf)))


@command("ME", replaces=True)
def cmd_me(c):
    """ME,<no>,<nfirst>,<nlast>,<xc>,<yc>,<zc>: multiple-excitation zone <no> = the interaction nodes
    <nfirst>..<nlast> (control point stored, not used)."""
    no = c.int(1, required=True, what="no")
    if no > 10:
        c.warn(f"motion number {no} > 10 (manual range 1-10; accepted)")
    toks = list(c.tokens)
    key = indexed_key(c.spec, toks)
    rec = build_record(c, "ME", toks)
    c.model.options.set_entry("ME", key, rec)
    report_problems(c, rec)
    c.confirm(f"ME {key} set")


@command("HOUSEX", replaces=True)
def cmd_housex(c):
    """HOUSEX,<optimize>,<supmode>,<nsim>,<nlssi>,<ansys>: HOUSE dialog fields without an argument
    (node optimizer, Linear/Quadratic superposition, stochastic simulations, non-linear soil SSI,
    ANSYS model input)."""
    rec = store_xrecord(c, "HOUSEX")
    if not 1 <= rec.nsim <= INC.MAX_SIM:
        c.warn(f"<nsim> must be 1..{INC.MAX_SIM} (G-27)")
    if rec.ansys:
        c.warn("HOUSEX ANSYS model input (Option AA) is not available in this build (tier P2): stored only")


@command("ANALYSX", replaces=True)
def cmd_analysx(c):
    """ANALYSX,<ffm>,<delrst>: Free-Field Load (0) / Free-Field Motion (1) application of the incoherency
    factors; delete restart files after a successful restart run."""
    store_xrecord(c, "ANALYSX")


# ======================================================================================
# BUILDFILE77
# ======================================================================================
def _model_dir(c) -> Path:
    m = c.model
    return Path(m.path) if m.path else c.interp.cwd


def _in_dir(c, name: str) -> Path:
    raw = os.path.expanduser(name.strip())
    if "\\" in raw and os.sep == "/":
        raw = raw.replace("\\", "/")
    p = Path(raw)
    return p if p.is_absolute() else _model_dir(c) / p


@command("BUILDFILE77")
def cmd_buildfile77(c):
    """BUILDFILE77,<out>,<in1>,...,<inN>: Build_FILE77 -- combine the FILE77s of per-level HOUSE runs into
    one FILE77 (spec 05b section 7; nodes in the FILE4 interaction order when <model>.N4 exists)."""
    out_n = c.str(1)
    names = [c.str(k) for k in range(2, c.nargs + 1) if c.given(k)]
    if not out_n or not names:
        c.fail("BUILDFILE77,<out>,<in1>,...,<inN>: the output and at least one input file are required")
    inputs = []
    for nm in names:
        p = _in_dir(c, nm)
        if not p.is_file():
            c.fail(f"{p} not found")
        try:
            inputs.append(INC.read_file77(p))
        except INC.IncoherencyError as exc:
            c.fail(str(exc))
    order = None
    f4 = _model_dir(c) / f"{c.model.name}.N4" if c.model.name else None
    if f4 is not None and f4.is_file():
        try:
            order = np.asarray(read_container(f4, kind="FILE4")["int_node"], np.int64)
        except (OSError, ValueError, KeyError) as exc:
            c.warn(f"{f4.name} not readable ({exc}): the nodes are kept in input order")
    try:
        arrays, meta, notes = INC.build_file77(inputs, names, order)
    except INC.IncoherencyError as exc:
        c.fail(str(exc))
    out = _in_dir(c, out_n)
    if out.exists():
        c.warn(f"{out.name} overwritten")
    write_container(out, "FILE77", arrays, meta, module="BUILDFILE77")
    for t in notes:
        c.warn(t)
    c.confirm(f"BUILDFILE77: {out.name} with {arrays['int_node'].size} interaction nodes from {len(names)} file(s), "
              f"{arrays['fnum'].size} frequencies" + (f" (FILE4 order of {f4.name})" if order is not None else ""))
