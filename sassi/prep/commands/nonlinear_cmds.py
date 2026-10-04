"""Option NON commands: nonlinear wall panels and springs, backbone curves, panel-model construction and the
NONLINEAR run (manual section 9.17; requirements 3.4.O, 3.4.R, 4.15; decisions D-NON-01 ... D-NON-12).

Data commands (stored in ``model.options``; WRITE writes them back, INP restores the same state)
------------------------------------------------------------------------------------------------
* ``EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>`` -- NONLINEAR header: EDF (default 0.8),
  element types (bit mask 1 panels, 2 springs, 4 beams = not available; blank = the types defined by P
  and S records), damping cut-off in percent (0 = none), damping scale (0 = 1), include elastic damping.
* ``P,<num>,<group>,<bbc>,<disp>,<force>`` -- wall panel: a SHELL group in a vertical plane, its BBC,
  Disp. Opt 1 (shear strain; 2 bending, experimental) and Force Opt 1 (CMS; 3 TAK experimental, 2 CMB
  not included).  Errors 126/128 apply unless ``EDUOPT,NONEXT,1`` (D-NON-05).
* ``S,<num>,<group>,<elem>,<bbc>,<disp>,<force>`` -- nonlinear spring: an existing SPRING element, its BBC,
  the DOF (1 X, 2 Y, 3 Z) and Force Opt 4 (GMR; Error 127 otherwise).
* ``B,<num>,<group>,<spgroup>,<bbc>,<force>,<end1>,<end2>`` -- nonlinear beam: stored, not usable.
* Backbone curves (origin implied; point 1 = cracking point; ``yield`` = 1-based index of the yield
  point): ``BBC,<num>,<type>,<points>,<yield>,<file>`` (X Y pairs from a file), ``BBCI,<num>,<yield>,<type>``,
  ``BBCP,<num>,<point>,<X>,<Y>`` (one point; ``point`` = npoints + 1 appends, skipped rows are (0, 0)),
  ``BBCX,<num>,<points>,<yield>,<X1>,...`` and ``BBCY,...``.  A curve is stored as its BBCI (yield, type),
  BBCX and BBCY entries; BBC and BBCP update them (so WRITE emits BBCI/BBCX/BBCY).  ``type`` is used for
  titles only.
* Deletion and listing (``<start>,[<end>],[<stride>]``, spec 11 section 0.2): DELBBC, DELBM, DELSPR, PDEL,
  DELNLS (NLSLAYER sets), PLIST.

Capacity and curve generation
-----------------------------
* ``SHEAR,<panel>,[fc],[fy],[P],[Nu],[Abe],[fybe]`` -- shear capacities of a panel (0 = all panels) by
  ACI 318-08, Wood 1990 (upper and lower bound), Barda 1977 and Gulec-Whittaker 2009
  (:mod:`sassi.core.panels`; arguments 6-7 are A_BE and f_y,BE, D-NON-10; ``EDUOPT,SHEARFORCEARGS,1``
  reads them as the forces F_VW and F_BE).
* ``BBCGEN,<Panel>,<ShearModel>,[fc],[fy],[Pn],[Nu],[bre],[bys],[CrackingForceLevel]`` -- 22-point BBC of
  a panel (0 = all) from its shear capacity (D-NON-09), written to the panel's BBC number (the panel
  number when 0, and the panel record is updated).

Panel-model construction (they edit the active model; run them on a copy)
-------------------------------------------------------------------------
WALLFLR (delete non-shell elements, one group per plane with >= 5 shells, the rest in group 1),
PANELIZE (split groups along their intersections with other groups), EDGE / EDGEMODEL (split walls
along the edges of openings), UNIPNL (one group per element), DGRDFLR (scale E of floor materials),
MERGEPANEL (append the groups and materials of a panel model), PNLGEN / PANELGEN (one panel per
vertical shell group), SOLIDPILE (pile-soil interface springs), NONLINMOTDISP (corner nodes of the
panels and end nodes of the S springs added to the MOTION and RELDISP requests).

Runs
----
* ``RUNNONLINEAR,[model]`` -- writes ``<model>.eql`` (the NONLINEAR deck, resolved from the model; CHECK
  errors 121-128 and the EDU-09 slope warnings are reported here) and runs NONLINEAR in the model
  directory.  AFWRITE of this build does not write the .eql (requirements tier P2), so RUNNONLINEAR does.
* ``COMBXYZTHD,<inpfile>`` -- COMB_XYZ_THD (D-NON-08) on ``COMB_XYZ_THD.inp`` (default name).
* ``NONLINBAT,<Sel>`` -- writes ``<model>_NONLINBAT.pre``, a generic batch script of the whole Option NON
  analysis (0 one direction, 1 X/Y/Z with COMB_XYZ_THD), and for Sel 1 a ``COMB_XYZ_THD.inp`` listing
  the THD files of the panel corners and spring ends (when the file does not exist).

SASSI-EDU extension commands (full names only, requirements 3.4.R)
-----------------------------------------------------------------
* ``NONLINITER,<var>[+<var2>...],[<maxit>]`` -- iteration driver: runs the command lines stored in the
  variable(s) once per pass (``#`` = pass number) until NONLINEAR reports convergence (D-NON-06:
  |dE/E| < 2 %, |d xi| < 0.5 %) or ``<maxit>`` passes (default 10).  Every pass must run RUNNONLINEAR.
* ``NONLINSAVE,<tag>,[<file8>]`` -- keep copies of the NONLINEAR results of a pass under the names of
  manual Fig. 1.2 (``Panel0001_<tag>.thd``, ``Panel_EQL_Matl_Prop_<tag>.txt``, ``<model>_<tag>.hou`` ...;
  ``<file8>`` = 1 also copies FILE8 / FILE8X/Y/Z to ``FILE8..._<tag>``; D-NON-07).
* ``NONLINRESET`` -- delete PANEL.NON, SPRING.NON, the ``*_EQL_Matl_Prop.txt`` files and the convergence
  history: the next NONLINEAR run is an elastic run (manual section 6.5.4).
* ``NONLINTHD,<prefix>`` -- copy the ``.THD`` files NONLINEAR reads (panel corners X/Y/Z, spring ends) to
  ``<prefix>_<name>`` so that the X, Y and Z RELDISP runs can be combined by COMB_XYZ_THD.
"""
from __future__ import annotations

import math
import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from ...conventions import nodal_result_name
from ...core import hysteresis as HY
from ...core import panels as PN
from ...model import (Element, Field, Group, Material, NodalRequest, Record, RelDispRequest, SpringProp, fmt_num,
                      make_record, register_record_type)
from ...model.values import NumberError, parse_float
from ..interpreter import LoopFrame
from ..options import eduopt
from ..registry import CommandReported, command, register
from . import select_ids
from .extensions import model_dir

TIER = "P2"
SHELL, TSHELL, SPRING, SOLID = 3, 5, 7, 1
PANEL_TYPES = (SHELL, TSHELL)


# ======================================================================================
# Typed records
# ======================================================================================
@register_record_type
class EqlRecord(Record):
    """``EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>``: NONLINEAR header."""
    COMMAND = "EQL"
    FIELDS = (Field("disp", float, 0.8, "Equivalent-Linear Displacement Factor (EDF)"),
              Field("nonlinopts", int, 0, "Nonlinear element types: 1 panels, 2 springs, 4 beams (blank = P/S given)"),
              Field("dampcutoff", float, 0.0, "Damping Cutoff % (0 = none)"),
              Field("dampscale", float, 0.0, "Damping Scale Factor (0 = 1)"),
              Field("elasticd", int, 0, "Include Elastic Damping (0/1)"))


@register_record_type
class PanelRecord(Record):
    """``P,<num>,<group>,<bbc>,<disp>,<force>``: wall panel."""
    COMMAND = "P"
    FIELDS = (Field("num", int, 0), Field("group", int, 0), Field("bbc", int, 0), Field("disp", int, 1),
              Field("force", int, 1))


@register_record_type
class SpringRecord(Record):
    """``S,<num>,<group>,<elem>,<bbc>,<disp>,<force>``: nonlinear spring."""
    COMMAND = "S"
    FIELDS = (Field("num", int, 0), Field("group", int, 0), Field("elem", int, 0), Field("bbc", int, 0),
              Field("disp", int, 1), Field("force", int, 4))


@register_record_type
class BeamRecord(Record):
    """``B,<num>,<group>,<spgroup>,<bbc>,<force>,<end1>,<end2>``: nonlinear beam (not usable)."""
    COMMAND = "B"
    FIELDS = (Field("num", int, 0), Field("group", int, 0), Field("spgroup", int, 0), Field("bbc", int, 0),
              Field("force", int, 4), Field("end1", int, 0), Field("end2", int, 0))


@register_record_type
class BbciRecord(Record):
    """``BBCI,<num>,<yield>,<type>``: yield point number and type of a backbone curve."""
    COMMAND = "BBCI"
    FIELDS = (Field("num", int, 0), Field("yield", int, 0), Field("type", int, 0))


def _check_numbers(c, first: int, last: int) -> None:
    for k in range(first, last + 1):
        t = c.raw(k)
        if t:
            try:
                parse_float(t)
            except NumberError:
                c.fail(f"argument {k}: '{t}' is not a number")


def _tokens(c, n: int) -> List[str]:
    return [c.raw(k) for k in range(1, n + 1)]


# ======================================================================================
# EQL, P, S, B
# ======================================================================================
@command("EQL", tier=TIER, max_args=5)
def cmd_eql(c):
    """EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>: NONLINEAR header -- EDF (0.7-0.9, default
    0.8), element types (1 panels + 2 springs; 4 beams not available; blank = the types with P/S records),
    damping cut-off % (0 none), damping scale (0 = 1), include elastic damping (0/1) (D-NON-02, D-NON-04)."""
    _check_numbers(c, 1, 5)
    rec = EqlRecord("EQL", _tokens(c, 5))
    disp = rec.get("disp")
    if not 0.0 < disp <= 1.0:
        c.fail(f"<disp> = {disp:g}: the equivalent-linear displacement factor must lie in (0, 1] (typically 0.7-0.9)")
    if not 0.6 <= disp <= 1.0:
        c.warn(f"EDF {disp:g} is outside the usual 0.7-0.9 range (0.6 gives too stiff, 1.0 too soft elements)")
    opts = rec.get("nonlinopts")
    if rec.given(2) and not 0 <= opts <= 7:
        c.fail("<NonLinOpts> is a bit mask: 1 panels, 2 springs, 4 beams (0..7)")
    if rec.given(2) and opts & 4:
        c.warn("nonlinear beams (bit 4) are not available in this version: ignored")
    cut = rec.get("dampcutoff")
    if cut < 0 or rec.get("dampscale") < 0:
        c.fail("<dampCutoff> and <dampScale> must be >= 0")
    if 0 < cut < 1:
        c.warn(f"<dampCutoff> is in percent: {cut:g} means {cut:g} % (give 7 for 7 %)")
    if rec.get("elasticd") not in (0, 1):
        c.fail("<ElasicD> must be 0 (hysteretic damping only) or 1 (add the elastic damping)")
    c.model.options.set_record(rec)
    c.confirm(f"EQL: EDF {disp:g}, elements "
              + (_opts_text(opts) if rec.given(2) else "= the P/S records") + f", cut-off "
              + (f"{cut:g} %" if cut > 0 else "none") + f", scale {rec.get('dampscale') or 1:g}, elastic damping "
              + ("added" if rec.get("elasticd") else "not added"))


def _opts_text(opts: int) -> str:
    parts = [n for b, n in ((1, "panels"), (2, "springs"), (4, "beams")) if opts & b]
    return " + ".join(parts) if parts else "none"


@command("P", tier=TIER, max_args=5)
def cmd_p(c):
    """P,<num>,<group>,<bbc>,<disp>,<force>: wall panel <num> = SHELL group <group> (vertical, coplanar) with
    backbone curve <bbc>, Disp. Opt (1 shear strain; 2 bending) and Force Opt (1 CMS; 3 TAK; 2 CMB not
    included).  Only 1/1 pass CHECK (Errors 126/128) unless EDUOPT,NONEXT,1."""
    _check_numbers(c, 1, 5)
    rec = PanelRecord("P", _tokens(c, 5))
    num, grp = rec.get("num"), rec.get("group")
    if num < 1 or grp < 1:
        c.fail("panel number <num> and group <group> must be >= 1")
    if rec.get("bbc") < 0:
        c.fail("<bbc> must be >= 0 (0 = assigned later by BBCGEN)")
    disp, force = rec.get("disp"), rec.get("force")
    if disp not in (1, 2):
        c.fail("<disp> must be 1 (shear strain) or 2 (bending rotation)")
    if force not in (1, 2, 3):
        c.fail("<force> must be 1 (CMS), 2 (CMB) or 3 (TAK)")
    if force == 2:
        c.warn("Force Opt 2 (Cheng-Mertz Bending) is not included in this version: stored only")
    elif force == 3 or disp == 2:
        c.warn(f"Force Opt {force} / Disp. Opt {disp}: CHECK Errors 126/128 in this version (set 1/1); EDUOPT,NONEXT,1 "
               "unlocks them as experimental (D-NON-05)")
    g = c.model.groups.get(grp)
    if g is None:
        c.warn(f"group {grp} is not defined (yet)")
    elif g.type not in PANEL_TYPES:
        c.warn(f"group {grp} is a {g.type_name} group: a panel must be a SHELL group")
    if c.model.options.entry("P", num) is not None:
        c.info(f"panel {num} redefined")
    c.model.options.set_entry("P", num, rec)
    c.confirm(f"panel {num}: group {grp}, BBC {rec.get('bbc')}, disp {disp}, force {HY.MODEL_CODES.get(force, force)}")


@command("S", tier=TIER, max_args=6)
def cmd_s(c):
    """S,<num>,<group>,<elem>,<bbc>,<disp>,<force>: nonlinear spring <num> = SPRING element <elem> of group
    <group> with backbone curve <bbc> on DOF <disp> (1 X, 2 Y, 3 Z translation) and Force Opt 4 (GMR,
    Error 127 otherwise)."""
    _check_numbers(c, 1, 6)
    rec = SpringRecord("S", _tokens(c, 6))
    num, grp, el = rec.get("num"), rec.get("group"), rec.get("elem")
    if num < 1 or grp < 1 or el < 1:
        c.fail("<num>, <group> and <elem> must be >= 1")
    if rec.get("bbc") < 1:
        c.warn("<bbc> = 0: give the backbone curve number of the spring")
    if rec.get("disp") not in (1, 2, 3):
        c.fail("<disp> must be 1 (X), 2 (Y) or 3 (Z): nonlinear springs act on one translational DOF")
    if rec.get("force") != 4:
        c.warn(f"Force Opt {rec.get('force')}: springs use the General Masing Rule (4); CHECK Error 127")
    g = c.model.groups.get(grp)
    if g is None:
        c.warn(f"group {grp} is not defined (yet)")
    elif g.type != SPRING:
        c.warn(f"group {grp} is a {g.type_name} group: S needs a SPRING element")
    elif el not in g.elements:
        c.warn(f"group {grp} has no element {el} (yet)")
    c.model.options.set_entry("S", num, rec)
    c.confirm(f"nonlinear spring {num}: group {grp} element {el}, BBC {rec.get('bbc')}, DOF {'XYZ'[rec.get('disp') - 1]}")


@command("B", tier=TIER, max_args=7)
def cmd_b(c):
    """B,<num>,<group>,<spgroup>,<bbc>,<force>,<end1>,<end2>: nonlinear beam -- not usable in this version
    (stored only, as the manual)."""
    _check_numbers(c, 1, 7)
    rec = BeamRecord("B", _tokens(c, 7))
    if rec.get("num") < 1:
        c.fail("<num> must be >= 1")
    c.model.options.set_entry("B", rec.get("num"), rec)
    c.warn("nonlinear beams are not usable in this version (stored only)")


# ======================================================================================
# Backbone curves
# ======================================================================================
def curve(model, num: int) -> Optional[Dict[str, object]]:
    """The stored backbone curve ``num``: dict(type, yield, x, y) (lists may differ in length), or None."""
    opts = model.options
    bi, bx, by = opts.entry("BBCI", num), opts.entry("BBCX", num), opts.entry("BBCY", num)
    if bi is None and bx is None and by is None:
        return None

    def vec(rec) -> List[float]:
        if rec is None:
            return []
        n = rec.integer(2, 0) or 0
        return [rec.number(k, 0.0) for k in range(4, 4 + n)]

    yld = None
    for rec, k in ((bi, 2), (bx, 3), (by, 3)):
        if rec is not None and rec.given(k):
            yld = rec.integer(k)
            break
    return dict(type=(bi.integer(3, 0) if bi is not None else 0) or 0, yield_=yld or 0, x=vec(bx), y=vec(by))


def curve_numbers(model) -> List[int]:
    keys: Set[int] = set()
    for name in ("BBCI", "BBCX", "BBCY"):
        keys.update(int(k) for k, _ in model.options.entries(name) if isinstance(k, int))
    return sorted(keys)


def curve_problems(model, num: int, strict: bool = False) -> List[str]:
    cv = curve(model, num)
    if cv is None:
        return [f"BBC {num} is not defined"]
    out = []
    if not cv["x"] or not cv["y"]:
        out.append(f"BBC {num}: " + ("X" if not cv["x"] else "Y") + " values missing (BBCX/BBCY)")
        return out
    if len(cv["x"]) != len(cv["y"]):
        out.append(f"BBC {num}: {len(cv['x'])} X values but {len(cv['y'])} Y values")
        return out
    if not cv["yield_"]:
        out.append(f"BBC {num}: the yield point number is not defined (BBCI, BBCX or BBCY <yield>)")
        return out
    out += [f"BBC {num}: {p}" for p in HY.backbone_problems(cv["x"], cv["y"], cv["yield_"], strict=strict)]
    return out


def backbone(model, num: int, strict: bool = False) -> HY.Backbone:
    probs = curve_problems(model, num, strict)
    if probs:
        raise ValueError("; ".join(probs))
    cv = curve(model, num)
    return HY.Backbone(cv["x"], cv["y"], cv["yield_"])


def _set_vector(model, name: str, num: int, yld: Optional[int], values: Sequence[float], tokens: Optional[Sequence[str]] = None) -> None:
    toks = [str(num), str(len(values)), "" if yld is None else str(yld)]
    toks += list(tokens) if tokens is not None else [fmt_num(float(v)) for v in values]
    model.options.set_entry(name, num, make_record(name, toks))


def _sync_yield(c, num: int, yld: int, source: str) -> None:
    """Keep the yield number of BBCI / BBCX / BBCY of curve ``num`` equal (the last command wins)."""
    opts = c.model.options
    for name, k in (("BBCI", 2), ("BBCX", 3), ("BBCY", 3)):
        if name == source:
            continue
        rec = opts.entry(name, num)
        if rec is None:
            continue
        if rec.given(k) and rec.integer(k) != yld:
            c.warn(f"BBC {num}: yield point number {rec.integer(k)} of {name} replaced by {yld}")
        rec.set_arg(k, yld)


def _store_curve(c, num: int, typ: Optional[int], yld: int, x: Sequence[float], y: Sequence[float],
                 xt: Optional[Sequence[str]] = None, yt: Optional[Sequence[str]] = None) -> None:
    opts = c.model.options
    bi = opts.entry("BBCI", num)
    t = typ if typ is not None else (bi.integer(3) if bi is not None and bi.given(3) else None)
    opts.set_entry("BBCI", num, BbciRecord("BBCI", [str(num), str(yld), "" if t is None else str(t)]))
    _set_vector(c.model, "BBCX", num, yld, x, xt)
    _set_vector(c.model, "BBCY", num, yld, y, yt)


def _report_curve(c, num: int) -> None:
    probs = curve_problems(c.model, num)
    for p in probs:
        c.warn(p)
    cv = curve(c.model, num)
    if cv and not probs:
        x, y = cv["x"], cv["y"]
        c.confirm(f"BBC {num}: {len(x)} points, yield point {cv['yield_']} ({x[cv['yield_'] - 1]:g}, "
                  f"{y[cv['yield_'] - 1]:g}), initial slope Y1/X1 = {y[0] / x[0]:.6g}"
                  + (f", type {HY.MODEL_CODES.get(cv['type'], cv['type'])}" if cv["type"] else ""))


@command("BBC", tier=TIER, max_args=5)
def cmd_bbc(c):
    """BBC,<num>,<type>,<points>,<yield>,<file>: backbone curve <num> from a file of X Y pairs (one pair per
    line, the origin not given; exactly <points> pairs are read); <type> 1 CMS, 2 CMB, 3 TAK, 4 GMR
    (titles only); <yield> = number of the yield point."""
    num = c.int(1, required=True, what="num")
    typ = c.int(2, default=None)
    npts = c.int(3, required=True, what="points")
    yld = c.int(4, required=True, what="yield")
    fname = c.str(5)
    if num < 1 or npts < 1:
        c.fail("<num> and <points> must be >= 1")
    if typ is not None and typ not in (0, 1, 2, 3, 4):
        c.fail("<type> must be 1 CMS, 2 CMB, 3 TAK or 4 GMR")
    if not 1 <= yld <= npts:
        c.fail(f"<yield> = {yld} must lie in 1..{npts}")
    if not fname:
        c.fail("<file> required")
    p = c.input_path(fname)
    pairs: List[Tuple[float, float]] = []
    try:
        lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        c.fail(f"cannot read {p}: {exc.strerror or exc}")
    for ln in lines:
        t = ln.replace(",", " ").split()
        if len(t) < 2:
            continue
        try:
            pairs.append((parse_float(t[0]), parse_float(t[1])))
        except NumberError:
            continue
    if len(pairs) < npts:
        c.fail(f"{p.name} holds {len(pairs)} X Y pairs, {npts} expected: curve not created")
    if len(pairs) > npts:
        c.warn(f"{p.name}: {len(pairs) - npts} pairs after the first {npts} ignored")
    pairs = pairs[:npts]
    _store_curve(c, num, typ, yld, [a for a, _ in pairs], [b for _, b in pairs])
    _report_curve(c, num)


@command("BBCI", tier=TIER, max_args=3)
def cmd_bbci(c):
    """BBCI,<num>,<yield>,<type>: set the yield point number and the type (1 CMS, 2 CMB, 3 TAK, 4 GMR) of
    backbone curve <num>."""
    _check_numbers(c, 1, 3)
    rec = BbciRecord("BBCI", _tokens(c, 3))
    num = rec.get("num")
    if num < 1:
        c.fail("<num> must be >= 1")
    if rec.given(3) and rec.get("type") not in (0, 1, 2, 3, 4):
        c.fail("<type> must be 1 CMS, 2 CMB, 3 TAK or 4 GMR")
    if rec.given(2):
        if rec.get("yield") < 1:
            c.fail("<yield> must be >= 1")
        _sync_yield(c, num, rec.get("yield"), "BBCI")
    c.model.options.set_entry("BBCI", num, rec)
    c.confirm(f"BBC {num}: yield point {rec.get('yield') or '(not set)'}, type "
              f"{HY.MODEL_CODES.get(rec.get('type'), 'not set')}")


@command("BBCP", tier=TIER, max_args=4)
def cmd_bbcp(c):
    """BBCP,<num>,<point>,<X>,<Y>: set point <point> of backbone curve <num> to (X, Y), as one row of the BBC
    grid of the NONLINEAR dialog; <point> = number of points + 1 appends a point (a curve can be built point
    by point after BBCI).  Rows skipped by a larger <point> are created as (0, 0), the empty rows of the grid,
    with a warning: RUNNONLINEAR rejects the curve until they are defined."""
    num = c.int(1, required=True, what="num")
    pt = c.int(2, required=True, what="point")
    x = c.float(3, required=True, what="X")
    y = c.float(4, required=True, what="Y")
    if num < 1 or pt < 1:
        c.fail("<num> and <point> must be >= 1")
    cv = curve(c.model, num) or dict(type=0, yield_=0, x=[], y=[])
    xs, ys = list(cv["x"]), list(cv["y"])
    n = max(len(xs), len(ys))
    if pt > n + 1:
        c.warn(f"BBC {num} had {n} points: points {n + 1}..{pt - 1} are set to (0, 0) until BBCP defines them")
    m = max(n, pt)
    xs += [0.0] * (m - len(xs))
    ys += [0.0] * (m - len(ys))
    xs[pt - 1], ys[pt - 1] = x, y
    bi = c.model.options.entry("BBCI", num)
    yld = cv["yield_"] or None
    if bi is None:
        c.model.options.set_entry("BBCI", num, BbciRecord("BBCI", [str(num), "" if yld is None else str(yld)]))
    _set_vector(c.model, "BBCX", num, yld, xs)
    _set_vector(c.model, "BBCY", num, yld, ys)
    c.confirm(f"BBC {num}: point {pt} = ({x:g}, {y:g}); {len(xs)} points")


def _cmd_vector(c, name: str) -> None:
    num = c.int(1, required=True, what="num")
    npts = c.int(2, required=True, what="points")
    yld = c.int(3, default=None)
    if num < 1 or npts < 1:
        c.fail("<num> and <points> must be >= 1")
    if yld is not None and not 1 <= yld <= npts:
        c.fail(f"<yield> = {yld} must lie in 1..{npts}")
    toks = [c.raw(k) for k in range(4, c.nargs + 1)]
    if len(toks) != npts or any(not t for t in toks):
        c.fail(f"{npts} values expected after <yield>, {len([t for t in toks if t])} given")
    vals = []
    for t in toks:
        try:
            vals.append(parse_float(t))
        except NumberError:
            c.fail(f"'{t}' is not a number")
    other = c.model.options.entry("BBCY" if name == "BBCX" else "BBCX", num)
    if other is not None and other.integer(2, 0) != npts:
        c.warn(f"BBC {num}: {npts} {name[-1]} values but {other.integer(2, 0)} {'Y' if name == 'BBCX' else 'X'} values "
               "(both vectors must have <points> values)")
    _set_vector(c.model, name, num, yld, vals, toks)
    if yld is not None:
        _sync_yield(c, num, yld, name)
    if c.model.options.entry("BBCI", num) is None:
        c.model.options.set_entry("BBCI", num, BbciRecord("BBCI", [str(num), "" if yld is None else str(yld)]))
    if other is not None and other.integer(2, 0) == npts:
        _report_curve(c, num)
    else:
        c.confirm(f"BBC {num}: {npts} {name[-1]} values set")


@command("BBCX", tier=TIER)
def cmd_bbcx(c):
    """BBCX,<num>,<points>,<yield>,<X1>,...,<Xn>: X values (strain or displacement, the origin not given) of
    backbone curve <num>."""
    _cmd_vector(c, "BBCX")


@command("BBCY", tier=TIER)
def cmd_bbcy(c):
    """BBCY,<num>,<points>,<yield>,<Y1>,...,<Yn>: Y values (force or moment) of backbone curve <num>."""
    _cmd_vector(c, "BBCY")


def _range(c, keys: Sequence[int], what: str, start_default: Optional[int] = None) -> List[int]:
    """``<start>,[<end>],[<stride>]`` of spec 11 section 0.2 (stride blank or < 1 -> 1 with a warning)."""
    keys = sorted(int(k) for k in keys)
    if not c.given(1) and start_default is None:
        c.int(1, required=True, what="start")
        a = 0
    else:
        a = c.int(1, default=start_default)
    b = c.int(2, default=a if start_default is None else (max(keys) if keys else a))
    inc = c.int(3, default=None)
    if c.given(3) and (inc is None or inc < 1):
        c.warn("stride < 1: the default stride 1 is used")
        inc = 1
    return select_ids(set(keys), a, b, inc or 1)


@command("DELBBC", tier=TIER, max_args=3)
def cmd_delbbc(c):
    """DELBBC,<start>,[<end>],[<stride>]: delete backbone curves (all their BBCI/BBCX/BBCY data)."""
    ids = _range(c, curve_numbers(c.model), "curve")
    for k in ids:
        for name in ("BBCI", "BBCX", "BBCY"):
            c.model.options.delete_entry(name, k)
    c.confirm(f"{len(ids)} backbone curves deleted")


def _delete_entries(c, name: str, what: str) -> None:
    keys = [k for k, _ in c.model.options.entries(name) if isinstance(k, int)]
    ids = _range(c, keys, what)
    for k in ids:
        c.model.options.delete_entry(name, k)
    c.confirm(f"{len(ids)} {what} deleted")


@command("PDEL", tier=TIER, max_args=3)
def cmd_pdel(c):
    """PDEL,<start>,[<end>],[<stride>]: delete wall panel records (the elements stay)."""
    _delete_entries(c, "P", "panels")


@command("DELSPR", tier=TIER, max_args=3)
def cmd_delspr(c):
    """DELSPR,<start>,[<end>],[<stride>]: delete nonlinear spring records (the springs stay as linear elements)."""
    _delete_entries(c, "S", "nonlinear springs")


@command("DELBM", tier=TIER, max_args=3)
def cmd_delbm(c):
    """DELBM,<start>,[<end>],[<stride>]: delete nonlinear beam records (the beams stay)."""
    _delete_entries(c, "B", "nonlinear beams")


def cmd_delnls(c):
    """DELNLS,<start>,[<end>],[<stride>]: delete NLSLAYER sets (nonlinear soil layers of SOIL-NON); the L and
    SOIL data are not changed."""
    _delete_entries(c, "NLSLAYER", "NLSLAYER sets")


# DELNLS belongs to the SOIL-NON data; it is registered as a replaceable handler so that a SOIL-NON command
# module may provide its own without an import conflict.
register("DELNLS", cmd_delnls, tier=TIER, max_args=3, placeholder=True, module=__name__,
         summary="DELNLS,<start>,[<end>],[<stride>]: delete NLSLAYER sets")


# ======================================================================================
# Model views shared by the commands
# ======================================================================================
def coordinates(model) -> Dict[int, np.ndarray]:
    if not model.nodes:
        return {}
    ids, xyz = model.global_coordinates()
    return {int(i): xyz[k] for k, i in enumerate(ids)}


def material_users(model) -> Dict[int, Set[int]]:
    """``{material: groups}`` of the M-table materials used by structural elements."""
    out: Dict[int, Set[int]] = {}
    for g in model.groups.values():
        if g.type in (SPRING, 9):
            continue
        for e in g.elements.values():
            out.setdefault(int(e.mat), set()).add(g.id)
    return out


def spring_prop_users(model) -> Dict[int, List[Tuple[int, int]]]:
    out: Dict[int, List[Tuple[int, int]]] = {}
    for g in model.groups.values():
        if g.type != SPRING:
            continue
        for e in g.elements.values():
            out.setdefault(int(e.prop), []).append((g.id, e.id))
    return out


def panel_geometry(model, gid: int, xyz: Optional[Dict[int, np.ndarray]] = None) -> PN.PanelGeometry:
    g = model.groups.get(gid)
    if g is None:
        raise PN.PanelError(f"group {gid} is not defined")
    if g.type not in PANEL_TYPES:
        raise PN.PanelError(f"group {gid} is a {g.type_name} group, not SHELL")
    if not g.elements:
        raise PN.PanelError(f"group {gid} has no elements")
    xyz = coordinates(model) if xyz is None else xyz
    return PN.panel_geometry([e.nodes for e in g.sorted_elements()], xyz)


def panel_section(model, gid: int) -> Tuple[float, int, List[str]]:
    """(thickness, material, warnings) of a panel group: uniform thickness and material expected."""
    g = model.groups[gid]
    th = sorted({float(e.thick) for e in g.elements.values()})
    mats = sorted({int(e.mat) for e in g.elements.values()})
    warns = []
    if len(th) > 1:
        warns.append(f"group {gid}: shells of different thicknesses {th}: the mean is used")
    if len(mats) > 1:
        warns.append(f"group {gid}: shells of several materials {mats} (each panel needs its own material, D-NON-12)")
    return float(np.mean(th)), mats[0], warns


def panels(model) -> List[Tuple[int, Record]]:
    return [(int(k), r) for k, r in model.options.entries("P") if isinstance(k, int)]


def springs(model) -> List[Tuple[int, Record]]:
    return [(int(k), r) for k, r in model.options.entries("S") if isinstance(k, int)]


# ======================================================================================
# SHEAR and BBCGEN
# ======================================================================================
def _capacity_args(c, k_fc: int):
    fc = c.float(k_fc, default=0.0)
    fy = c.float(k_fc + 1, default=0.0)
    rho = c.float(k_fc + 2, default=0.0)
    nu = c.float(k_fc + 3, default=0.0)
    a6 = c.float(k_fc + 4, default=0.0)
    a7 = c.float(k_fc + 5, default=0.0)
    if fc <= 0:
        c.fail("f'c must be > 0")
    if fy < 0 or rho < 0 or a6 < 0 or a7 < 0:
        c.fail("fy, the reinforcement ratio and the boundary-element data must be >= 0")
    return fc, fy, rho, nu, a6, a7


def _panel_capacities(c, num: int, rec: Record, fc, fy, rho, nu, a6, a7, xyz) -> Tuple[PN.ShearCapacities, PN.PanelGeometry, float]:
    gid = rec.integer(2, 0)
    geo = panel_geometry(c.model, gid, xyz)
    th, _, warns = panel_section(c.model, gid)
    for w in warns + geo.warnings:
        c.warn(f"panel {num}: {w}")
    force_args = eduopt(c.model, "SHEARFORCEARGS") == "1"
    cap = PN.shear_capacities(geo.height_extent, geo.width_extent, th, fc, fy, rho, nu, a6, a7, c.model.gravity,
                              force_args=force_args)
    return cap, geo, th


def _selected_panels(c, num: int) -> List[Tuple[int, Record]]:
    allp = panels(c.model)
    if not allp:
        c.fail("no panels defined (P or PNLGEN)")
    if num == 0:
        return allp
    sel = [(k, r) for k, r in allp if k == num]
    if not sel:
        c.fail(f"panel {num} is not defined")
    return sel


@command("SHEAR", tier=TIER, max_args=7)
def cmd_shear(c):
    """SHEAR,<panel>,[fc],[fy],[P],[Nu],[Abe],[fybe]: shear capacities of wall panel <panel> (0 = all) by
    ACI 318-08, Wood 1990 (upper/lower bound), Barda 1977 and Gulec-Whittaker 2009; fc, fy in ksi (kN/m2),
    P = web reinforcement ratio, Nu axial force (kips, kN), Abe/fybe boundary-element area and yield
    stress (EDUOPT,SHEARFORCEARGS,1: forces F_VW, F_BE).  Results in kips (kN); the panel geometry comes
    from the model (h_W vertical extent, l_W horizontal extent, t_W shell thickness)."""
    num = c.int(1, default=0)
    fc, fy, rho, nu, a6, a7 = _capacity_args(c, 2)
    units, _, gwarn = PN.unit_system(c.model.gravity)
    if gwarn:
        c.warn(gwarn)
    sel = _selected_panels(c, num)
    xyz = coordinates(c.model)
    funit = "kips" if units == "british" else "kN"
    c.info(f"SHEAR ({'British: ft, ksi, kips' if units == 'british' else 'SI: m, kN/m2, kN'}): panel, ACI 318-08, "
           f"Wood 1990 upper, Wood 1990 lower, Barda 1977, Gulec-Whittaker 2009 [{funit}]")
    for k, rec in sel:
        try:
            cap, geo, th = _panel_capacities(c, k, rec, fc, fy, rho, nu, a6, a7, xyz)
        except PN.PanelError as exc:
            c.error(f"panel {k}: {exc}")
            continue
        c.info(f"{k:6d} {cap.aci:14.6g} {min(cap.wood_raw, cap.upper):14.6g} {cap.wood_lower:14.6g} "
               f"{cap.barda:14.6g} {cap.gw:14.6g}")
        c.info(f"       h_W {geo.height_extent:g}, l_W {geo.width_extent:g}, t_W {th:g}; alpha_c {cap.alpha_c:g}; raw "
               f"ACI {cap.aci_raw:.6g}, raw Wood {cap.wood_raw:.6g}; upper bound 10 sqrt(f'c) A_W {cap.upper:.6g}")
    c.confirm(f"SHEAR: {len(sel)} panel(s)")


@command("BBCGEN", tier=TIER, max_args=9)
def cmd_bbcgen(c):
    """BBCGEN,<Panel>,<ShearModel>,[fc],[fy],[Pn],[Nu],[bre],[bys],[CrackingForceLevel]: 22-point backbone
    curve of panel <Panel> (0 = all) from its ultimate shear by ShearModel 1 ACI 318-08, 2 Wood 1990, 3 Barda
    1977, 4 Gulec-Whittaker 2009 (D-NON-09): cracking point V_cr = 3 sqrt(f'c) A_W (CrackingForceLevel 0) or
    CFL V_u (0.10-0.50), gamma_cr = V_cr/(G A_W); 20 points to the yield point (0.004, V_u); failure point
    (0.02, 1.02 V_u); type 1 (CMS); written to the panel's BBC (the panel number when 0)."""
    num = c.int(1, default=0)
    model_no = c.int(2, required=True, what="ShearModel")
    if model_no not in PN.SHEAR_MODELS:
        c.fail("<ShearModel> must be 1 ACI 318-08, 2 Wood 1990, 3 Barda 1977 or 4 Gulec-Whittaker 2009")
    fc, fy, rho, nu, a6, a7 = _capacity_args(c, 3)
    cfl = c.float(9, default=0.0)
    if cfl != 0.0 and not 0.10 <= cfl <= 0.50:
        c.fail(f"CrackingForceLevel {cfl:g}: 0 (cracking stress 3 sqrt(f'c)) or a ratio in [0.10, 0.50]")
    _, _, gwarn = PN.unit_system(c.model.gravity)
    if gwarn:
        c.warn(gwarn)
    sel = _selected_panels(c, num)
    xyz = coordinates(c.model)
    m = c.model
    made = 0
    for k, rec in sel:
        try:
            cap, geo, th = _panel_capacities(c, k, rec, fc, fy, rho, nu, a6, a7, xyz)
            _, mid, _ = panel_section(m, rec.integer(2, 0))
            mat = m.materials.get(mid)
            if mat is None:
                raise PN.PanelError(f"material {mid} is not defined")
            G = mat.constants(m.gravity).G
            vu = cap.ultimate(model_no)
            vcr = cap.cracking if cfl == 0.0 else cfl * vu
            area = th * geo.width_extent
            g, v, yi = PN.bbcgen_curve(vu, vcr, G * area)
        except (PN.PanelError, ValueError) as exc:
            c.error(f"panel {k}: {exc}")
            continue
        bnum = rec.integer(3, 0) or k
        if not rec.integer(3, 0):
            rec.set_arg(3, bnum)
            c.info(f"panel {k}: BBC number set to {bnum}")
        _store_curve(c, bnum, 1, yi, g, v)
        s2 = (v[1] - v[0]) / (g[1] - g[0])
        if s2 > G * area * (1 + 1e-12):
            c.warn(f"panel {k}: the second BBC segment ({s2:.4g}) is steeper than G A_W ({G * area:.4g})")
        made += 1
        c.info(f"panel {k}: BBC {bnum} ({PN.SHEAR_MODELS[model_no]}): V_u {vu:.6g}, V_cr {vcr:.6g} at gamma_cr "
               f"{g[0]:.4e} (G A_W {G * area:.6g}), yield (0.004, V_u), failure (0.02, {1.02 * vu:.6g})")
    c.confirm(f"BBCGEN: {made} backbone curve(s) generated")


# ======================================================================================
# Panel-model construction
# ======================================================================================
def _model_size(xyz: Dict[int, np.ndarray]) -> float:
    if not xyz:
        return 1.0
    P = np.array(list(xyz.values()))
    return float(max(np.max(np.ptp(P, axis=0)), 1e-30))


def _warn_requests(c, what: str) -> None:
    m = c.model
    if m.eout:
        c.warn(f"{what}: the EOUT requests refer to the old groups -- check them (EOUT,0 clears)")
    if m.options.entries("P") or m.options.entries("S"):
        c.warn(f"{what}: panel/spring records refer to group numbers -- re-run PNLGEN or redefine P/S")
    if c.interp.session.get("cuts"):
        c.warn(f"{what}: the cuts refer to the old groups and elements")


def _next_group(m) -> int:
    return (max(m.groups) + 1) if m.groups else 1


def _orientation_title(n: np.ndarray, d: float) -> str:
    ax = np.abs(n)
    c1 = math.cos(math.radians(1.0))
    if ax[2] >= c1:
        return f"FLOOR z={fmt_num(round(d * np.sign(n[2]), 9))}"
    if ax[1] >= c1:
        return f"WALL along X (XZ plane) y={fmt_num(round(d * np.sign(n[1]), 9))}"
    if ax[0] >= c1:
        return f"WALL along Y (YZ plane) x={fmt_num(round(d * np.sign(n[0]), 9))}"
    if abs(n[2]) <= PN.VERTICAL_TOL:
        return "OBLIQUE WALL"
    return "OBLIQUE"


@command("WALLFLR", tier=TIER, max_args=0)
def cmd_wallflr(c):
    """WALLFLR: separate shell walls and floors -- deletes every non-shell element of the active model, puts
    each set of >= 5 coplanar shells in its own group (groups 2, 3 ...; titles by orientation) and the other
    shells in group 1.  Run it on a copy of the model (CPMODEL), then PANELIZE / EDGE / GROUPMAT."""
    m = c.model
    xyz = coordinates(m)
    tol = max(1e-6 * _model_size(xyz), float(eduopt(m, "GEOMTOL", "0") or 0))
    shells: List[Tuple[int, Element]] = []
    removed = 0
    for gid in sorted(m.groups):
        g = m.groups[gid]
        if g.type in PANEL_TYPES:
            shells += [(g.type, e) for e in g.sorted_elements()]
        else:
            removed += len(g.elements)
    if not shells:
        c.fail("the model has no shell elements")
    normals, offs = [], []
    for _, e in shells:
        nd = PN.element_nodes(e.nodes)
        try:
            n, d = PN.shell_plane(np.array([xyz[k] for k in nd]))
        except (KeyError, PN.PanelError) as exc:
            c.fail(f"shell element with nodes {nd}: {exc}")
        normals.append(n)
        offs.append(d)
    lab = PN.plane_clusters(np.array(normals), np.array(offs), angle_tol=1e-3, offset_tol=tol)
    groups: Dict[int, Group] = {}
    rest = [i for i in range(len(shells)) if np.sum(lab == lab[i]) < 5]
    nxt = 2
    for k in sorted(set(int(v) for v in lab), key=lambda v: int(np.flatnonzero(lab == v)[0])):
        idx = [int(i) for i in np.flatnonzero(lab == k)]
        if len(idx) < 5:
            continue
        typ = shells[idx[0]][0]
        g = Group(nxt, typ, _orientation_title(normals[idx[0]], offs[idx[0]]))
        for j, i in enumerate(idx, start=1):
            g.elements[j] = shells[i][1].copy(new_id=j)
        groups[nxt] = g
        nxt += 1
    if rest:
        g = Group(1, shells[rest[0]][0], "WALLFLR: shells of no wall or floor of >= 5 elements")
        for j, i in enumerate(rest, start=1):
            g.elements[j] = shells[i][1].copy(new_id=j)
        groups[1] = g
    _warn_requests(c, "WALLFLR")
    m.groups = dict(sorted(groups.items()))
    m.group_active = max(m.groups) if m.groups else None
    if removed:
        c.warn(f"{removed} non-shell elements deleted")
    c.confirm(f"WALLFLR: {len(groups) - (1 if rest else 0)} wall/floor groups"
              + (f", {len(rest)} other shells in group 1" if rest else ""))


@command("PANELIZE", tier=TIER, max_args=0)
def cmd_panelize(c):
    """PANELIZE: split every shell group of a panel model (after WALLFLR) along its intersections with the
    other shell groups (edges whose two nodes also belong to another group); each connected part becomes a
    group (the first keeps the number, the others are appended)."""
    m = c.model
    shell_groups = [gid for gid in sorted(m.groups) if m.groups[gid].type in PANEL_TYPES]
    if not shell_groups:
        c.fail("the model has no shell groups")
    owners: Dict[int, Set[int]] = {}
    for gid in shell_groups:
        for e in m.groups[gid].elements.values():
            for n in PN.element_nodes(e.nodes):
                owners.setdefault(n, set()).add(gid)
    made = 0
    for gid in shell_groups:
        g = m.groups[gid]
        els = g.sorted_elements()
        edge_map: Dict[Tuple[int, int], List[int]] = {}
        for i, e in enumerate(els):
            for ed in PN.element_edges(e.nodes):
                edge_map.setdefault(ed, []).append(i)
        pairs = []
        for (a, b), idx in edge_map.items():
            cut = bool((owners.get(a, set()) & owners.get(b, set())) - {gid})
            if not cut:
                for i in idx[1:]:
                    pairs.append((idx[0], i))
        comp = PN.connected_components(len(els), pairs)
        parts = sorted(set(int(v) for v in comp))
        if len(parts) <= 1:
            continue
        keep = [els[i] for i in range(len(els)) if comp[i] == parts[0]]
        g.elements = {e.id: e for e in keep}
        for p in parts[1:]:
            ng = Group(_next_group(m), g.type, f"{g.title} (part)".strip())
            for j, i in enumerate([i for i in range(len(els)) if comp[i] == p], start=1):
                ng.elements[j] = els[i].copy(new_id=j)
            m.groups[ng.id] = ng
            made += 1
    if made:
        _warn_requests(c, "PANELIZE")
    c.confirm(f"PANELIZE: {made} new groups (MERGEGROUP / GCOM combine them if needed)")


def _edge_split(c, gid: int, flags: Tuple[int, int, int], xyz) -> int:
    m = c.model
    g = m.groups.get(gid)
    if g is None:
        raise PN.PanelError(f"group {gid} is not defined")
    if g.type not in PANEL_TYPES:
        raise PN.PanelError(f"group {gid} is not a shell group")
    els = g.sorted_elements()
    nodes = sorted({n for e in els for n in PN.element_nodes(e.nodes)})
    P = np.array([xyz[n] for n in nodes])
    cen, nrm, _ = PN.plane_fit(P)
    if abs(nrm[2]) > PN.VERTICAL_TOL:
        raise PN.PanelError(f"group {gid} is not a wall (vertical plane)")
    eh = PN.horizontal_axis(nrm)
    count: Dict[Tuple[int, int], int] = {}
    for e in els:
        for ed in PN.element_edges(e.nodes):
            count[ed] = count.get(ed, 0) + 1
    bnd = [ed for ed, k in count.items() if k == 1]
    c1 = math.cos(math.radians(1.0))
    segs, skip = [], []
    for a, b in bnd:
        pa, pb = np.asarray(xyz[a], float), np.asarray(xyz[b], float)
        d = (pb - pa) / max(np.linalg.norm(pb - pa), 1e-300)
        cls = next((k for k in range(3) if abs(d[k]) >= c1), None)
        segs.append((pa, pb))
        skip.append(cls is not None and flags[cls] == 1)
    size = float(np.max(np.ptp(P, axis=0)))
    labels, _ = PN.edge_cells([np.array([xyz[n] for n in PN.element_nodes(e.nodes)]) for e in els], eh, segs, skip,
                              1e-6 * max(size, 1e-30))
    cells = sorted(set(int(v) for v in labels))
    if len(cells) <= 1:
        return 0
    cent = {}
    for k in cells:
        pts = np.array([np.mean([xyz[n] for n in PN.element_nodes(els[i].nodes)], axis=0)
                        for i in range(len(els)) if labels[i] == k])
        cc = pts.mean(axis=0)
        q = 1e-6 * max(size, 1e-30)
        cent[k] = (-round(float(cc[2]) / q), round(float(cc @ eh) / q))
    order = sorted(cells, key=lambda k: cent[k])
    keep = [els[i] for i in range(len(els)) if labels[i] == order[0]]
    g.elements = {}
    for j, e in enumerate(keep, start=1):
        g.elements[j] = e.copy(new_id=j)
    for k in order[1:]:
        ng = Group(_next_group(m), g.type, f"{g.title} (EDGE {gid})".strip())
        for j, i in enumerate([i for i in range(len(els)) if labels[i] == k], start=1):
            ng.elements[j] = els[i].copy(new_id=j)
        m.groups[ng.id] = ng
    return len(order) - 1


def _flags(c, k: int) -> Tuple[int, int, int]:
    out = []
    for j in range(3):
        v = c.int(k + j, default=0)
        if v not in (0, 1):
            c.fail("edge flags are 0 (use the edges parallel to the axis) or 1 (ignore them)")
        out.append(v)
    return tuple(out)


@command("EDGE", tier=TIER, max_args=4)
def cmd_edge(c):
    """EDGE,<panel>,[X],[Y],[Z]: split wall group <panel> along the lines of its outer and opening edges
    (flags: 0 use the edges parallel to that global axis, 1 ignore them; oblique edges are always used);
    the cells (top first, then left to right) become groups: the first keeps the number, the others are
    appended (manual Fig. 1.5)."""
    gid = c.int(1, required=True, what="panel group")
    flags = _flags(c, 2)
    try:
        n = _edge_split(c, gid, flags, coordinates(c.model))
    except (PN.PanelError, KeyError) as exc:
        c.fail(str(exc))
    if n:
        _warn_requests(c, "EDGE")
    c.confirm(f"EDGE: group {gid} split into {n + 1} groups" if n else f"EDGE: group {gid} not split (no interior edge line)")


@command("EDGEMODEL", tier=TIER, max_args=3)
def cmd_edgemodel(c):
    """EDGEMODEL,[x],[y],[z]: EDGE on every wall (vertical shell) group of the model with the same flags."""
    flags = _flags(c, 1)
    xyz = coordinates(c.model)
    total = 0
    for gid in sorted(c.model.groups):
        g = c.model.groups[gid]
        if g.type not in PANEL_TYPES or not g.elements:
            continue
        try:
            total += _edge_split(c, gid, flags, xyz)
        except PN.PanelError:
            continue                       # floors and oblique shells are not walls
    if total:
        _warn_requests(c, "EDGEMODEL")
    c.confirm(f"EDGEMODEL: {total} new groups")


@command("UNIPNL", tier=TIER, max_args=1)
def cmd_unipnl(c):
    """UNIPNL,<group>: one group per element of <group> (element 1 stays, the others go to new appended
    groups, each numbered element 1) -- one panel per shell, e.g. for curved walls."""
    gid = c.int(1, required=True, what="group")
    m = c.model
    g = m.groups.get(gid)
    if g is None:
        c.fail(f"group {gid} is not defined")
    els = g.sorted_elements()
    if len(els) <= 1:
        c.confirm(f"UNIPNL: group {gid} has {len(els)} element(s): nothing to do")
        return
    g.elements = {1: els[0].copy(new_id=1)}
    for e in els[1:]:
        ng = Group(_next_group(m), g.type, f"{g.title} element {e.id}".strip())
        ng.elements[1] = e.copy(new_id=1)
        m.groups[ng.id] = ng
    _warn_requests(c, "UNIPNL")
    c.confirm(f"UNIPNL: group {gid} split into {len(els)} groups (run GROUPMAT for one material per panel)")


def scale_modulus(mat: Material, s: float) -> None:
    if mat.mtype == 1:
        mat.val1 *= s
    elif mat.mtype == 2:
        mat.val1 *= s
        mat.val2 *= s
    else:
        mat.val1 *= math.sqrt(s)
        mat.val2 *= math.sqrt(s)


@command("DGRDFLR", tier=TIER, max_args=1)
def cmd_dgrdflr(c):
    """DGRDFLR,<scale>: multiply Young's modulus of the materials of the floor panels (horizontal shell
    groups) by <scale> (after WALLFLR and GROUPMAT); materials also used by other groups are skipped."""
    s = c.float(1, required=True, what="scale")
    if s <= 0:
        c.fail("<scale> must be > 0")
    m = c.model
    xyz = coordinates(m)
    floors: Set[int] = set()
    for gid, g in m.groups.items():
        if g.type not in PANEL_TYPES or not g.elements:
            continue
        nodes = sorted({n for e in g.elements.values() for n in PN.element_nodes(e.nodes)})
        if len(nodes) < 3:
            continue
        _, n, _ = PN.plane_fit(np.array([xyz[k] for k in nodes]))
        if abs(n[2]) >= math.cos(math.radians(1.0)):
            floors.add(gid)
    users = material_users(m)
    done = []
    for mid in sorted({int(e.mat) for gid in floors for e in m.groups[gid].elements.values()}):
        others = users.get(mid, set()) - floors
        if others:
            c.warn(f"material {mid} is also used by groups {sorted(others)} (not floors): not scaled (run GROUPMAT)")
            continue
        if mid not in m.materials:
            c.warn(f"material {mid} is not defined")
            continue
        scale_modulus(m.materials[mid], s)
        done.append(mid)
    c.confirm(f"DGRDFLR: E x {s:g} for {len(done)} floor materials ({len(floors)} floor groups)")


@command("MERGEPANEL", tier=TIER, max_args=1)
def cmd_mergepanel(c):
    """MERGEPANEL,<Panel>: append the groups and materials of panel model <Panel> to the active (original)
    model whose shell groups were deleted; nodes are matched by number (coordinates checked), missing nodes
    are added; element materials are renumbered after the active model's."""
    num = c.int(1, required=True, what="panel model")
    if num not in c.interp.models:
        c.fail(f"model {num} is not in memory")
    if num == c.interp.active_model:
        c.fail("<Panel> must be another model than the active one")
    src, m = c.interp.models[num], c.model
    if any(g.type in PANEL_TYPES for g in m.groups.values()):
        c.warn("the active model still has shell groups (delete them first: the panel model replaces them)")
    sxyz, mxyz = coordinates(src), coordinates(m)
    tol = 1e-6 * max(_model_size(sxyz), 1e-30)
    used = sorted({n for g in src.groups.values() for e in g.elements.values() for n in PN.element_nodes(e.nodes)})
    added = 0
    for n in used:
        if n in m.nodes:
            if np.linalg.norm(np.asarray(mxyz[n]) - np.asarray(sxyz[n])) > tol:
                c.fail(f"node {n} has other coordinates in the two models: nodes must match by number")
            continue
        m.define_node(n, tuple(float(v) for v in sxyz[n]), 0)
        added += 1
    mat_off = max(m.materials, default=0)
    mp = {}
    for mid in sorted(src.materials):
        s = src.materials[mid]
        mp[mid] = mat_off + mid
        m.materials[mp[mid]] = Material(mp[mid], s.val1, s.val2, s.weight, s.pdamp, s.sdamp, s.mtype)
    g0 = _next_group(m)
    for k, gid in enumerate(sorted(src.groups)):
        sg = src.groups[gid]
        ng = Group(g0 + k, sg.type, sg.title)
        for e in sg.sorted_elements():
            ne = e.copy()
            if sg.type not in (SPRING, 9):
                ne.mat = mp.get(e.mat, e.mat)
            ng.elements[ne.id] = ne
        m.groups[ng.id] = ng
    if added:
        c.warn(f"{added} nodes of the panel model were missing and have been added")
    c.confirm(f"MERGEPANEL: {len(src.groups)} groups appended as {g0}..{g0 + len(src.groups) - 1}, "
              f"{len(src.materials)} materials as {mat_off + 1}..{mat_off + max(src.materials, default=0)}")


@command("PNLGEN", tier=TIER, abbrev=("PANELGEN",), max_args=0)
def cmd_pnlgen(c):
    """PNLGEN (PANELGEN): one panel record per vertical shell group, numbered 1, 2 ... in group order, with
    BBC = panel number, Disp. Opt 1 and Force Opt 1 (replaces the existing P records)."""
    m = c.model
    xyz = coordinates(m)
    old = panels(m)
    for k, _ in old:
        m.options.delete_entry("P", k)
    k = 0
    for gid in sorted(m.groups):
        g = m.groups[gid]
        if g.type not in PANEL_TYPES or not g.elements:
            continue
        first = g.sorted_elements()[0]
        nd = PN.element_nodes(first.nodes)
        try:
            n, _ = PN.shell_plane(np.array([xyz[v] for v in nd]))
        except (KeyError, PN.PanelError):
            continue
        if abs(n[2]) > PN.VERTICAL_TOL:
            continue
        k += 1
        m.options.set_entry("P", k, PanelRecord("P", [str(k), str(gid), str(k), "1", "1"]))
    if old:
        c.warn(f"{len(old)} existing panel records replaced")
    c.confirm(f"PNLGEN: {k} panels (BBC numbers = panel numbers; define them with BBCGEN or BBCX/BBCY)")


@command("SOLIDPILE", tier=TIER, max_args=4)
def cmd_solidpile(c):
    """SOLIDPILE,<group>,[stiff],[soft],[stiff2]: separate the SOLID pile group <group> from the soil -- its
    interface nodes are duplicated, the pile elements reconnected to the duplicates and the two joined by
    springs (4 % damping) in 4 new groups: side-wall X and Y (stiff, 1e7), side-wall Z (soft, 10), tip Z
    (stiff2, 1e7); each spring has its own SC property (targets of S with GMR)."""
    gid = c.int(1, required=True, what="group")
    stiff = c.float(2, default=1.0e7)
    soft = c.float(3, default=10.0)
    stiff2 = c.float(4, default=1.0e7)
    m = c.model
    g = m.groups.get(gid)
    if g is None:
        c.fail(f"group {gid} is not defined")
    if g.type != SOLID:
        c.fail(f"group {gid} is a {g.type_name} group: SOLIDPILE needs SOLID piles")
    if min(stiff, soft, stiff2) < 0:
        c.fail("spring constants must be >= 0")
    pile_nodes = {n for e in g.elements.values() for n in PN.element_nodes(e.nodes)}
    outside: Set[int] = set()
    for og in m.groups.values():
        if og.id == gid:
            continue
        for e in og.elements.values():
            outside.update(PN.element_nodes(e.nodes))
    iface = sorted(pile_nodes & outside)
    if not iface:
        c.fail(f"group {gid} shares no node with other elements")
    xyz = coordinates(m)
    # pile tips: the lowest interface level of each connected pile
    els = g.sorted_elements()
    pairs = []
    first: Dict[int, int] = {}
    for i, e in enumerate(els):
        for n in PN.element_nodes(e.nodes):
            if n in first:
                pairs.append((first[n], i))
            else:
                first[n] = i
    comp = PN.connected_components(len(els), pairs)
    node_comp = {n: int(comp[first[n]]) for n in pile_nodes}
    size = _model_size(xyz)
    zmin: Dict[int, float] = {}
    for n in pile_nodes:
        zmin[node_comp[n]] = min(zmin.get(node_comp[n], float("inf")), float(xyz[n][2]))
    tip = {n for n in iface if abs(float(xyz[n][2]) - zmin[node_comp[n]]) <= 1e-6 * size}
    nxt = max(m.nodes) + 1
    dup: Dict[int, int] = {}
    for n in iface:
        nd = m.nodes[n]
        m.define_node(nxt, nd.xyz, nd.csys)
        m.nodes[nxt].fix = list(nd.fix)
        for k in (3, 4, 5):                    # solid nodes have no rotations; springs would add them
            m.nodes[nxt].fix[k] = 1
            nd.fix[k] = 1
        dup[n] = nxt
        nxt += 1
    for e in g.elements.values():
        e.nodes = [dup.get(n, n) for n in e.nodes]
    sc = max(m.springs, default=0) + 1
    made = []
    for title, sel, k6 in (("side wall X", [n for n in iface if n not in tip], (stiff, 0, 0)),
                           ("side wall Y", [n for n in iface if n not in tip], (0, stiff, 0)),
                           ("side wall Z", [n for n in iface if n not in tip], (0, 0, soft)),
                           ("pile tip Z", sorted(tip), (0, 0, stiff2))):
        ng = Group(_next_group(m), SPRING, f"SOLIDPILE group {gid}: {title}")
        for j, n in enumerate(sel, start=1):
            m.springs[sc] = SpringProp(sc, float(k6[0]), float(k6[1]), float(k6[2]), 0.0, 0.0, 0.0, 0.04)
            ng.elements[j] = Element(j, [n, dup[n]], mat=1, prop=sc)
            sc += 1
        m.groups[ng.id] = ng
        made.append((ng.id, len(sel)))
    c.confirm(f"SOLIDPILE: {len(iface)} interface nodes duplicated ({len(tip)} at the tips); spring groups "
              + ", ".join(f"{gg} ({k} springs)" for gg, k in made))


# ======================================================================================
# NONLINMOTDISP and PLIST
# ======================================================================================
def element_nodes_needed(model, xyz=None) -> Tuple[List[Tuple[int, int]], List[str]]:
    """(node, dof) pairs whose .THD NONLINEAR reads (panel corners X/Y/Z, spring ends on their DOF)."""
    xyz = coordinates(model) if xyz is None else xyz
    need: List[Tuple[int, int]] = []
    warns: List[str] = []
    for k, rec in panels(model):
        try:
            geo = panel_geometry(model, rec.integer(2, 0), xyz)
        except PN.PanelError as exc:
            warns.append(f"panel {k}: {exc}")
            continue
        for w in geo.warnings:
            warns.append(f"panel {k}: {w}")
        for n in geo.corners:
            for d in (1, 2, 3):
                need.append((int(n), d))
    for k, rec in springs(model):
        g = model.groups.get(rec.integer(2, 0))
        e = g.elements.get(rec.integer(3, 0)) if g is not None else None
        if e is None or len(PN.element_nodes(e.nodes)) < 2:
            warns.append(f"spring {k}: group {rec.integer(2, 0)} element {rec.integer(3, 0)} is not a 2-node element")
            continue
        d = rec.integer(5, 1)
        for n in PN.element_nodes(e.nodes)[:2]:
            need.append((int(n), int(d)))
    return sorted(set(need)), warns


@command("NONLINMOTDISP", tier=TIER, abbrev=("NONLINMODISP",), max_args=0)
def cmd_nonlinmotdisp(c):
    """NONLINMOTDISP: add the four corner nodes of every panel (node-connection counting) and the end nodes
    of every S spring to the MOTION (NOUT, transfer functions) and RELDISP (RDND) output requests, without
    duplicates.  Run it before AFWRITE."""
    m = c.model
    need, warns = element_nodes_needed(m)
    for w in warns:
        c.warn(w)
    if not need:
        c.fail("no panel corner or spring node found (define P / S first)")
    have = {(n, r.dir) for r in m.nout for n in r.nodes}
    add_nout: Dict[int, List[int]] = {}
    for n, d in need:
        if (n, d) not in have:
            add_nout.setdefault(d, []).append(n)
    for d in sorted(add_nout):
        m.nout.append(NodalRequest(d, [1, 0, 0, 0, 0, 0], sorted(add_nout[d])))
    rd = {r.node: r for r in m.rdnd}
    added_rd = 0
    for n, d in need:
        r = rd.get(n)
        if r is None:
            r = RelDispRequest(n, [0] * 6)
            m.rdnd.append(r)
            rd[n] = r
            added_rd += 1
        r.flags[d - 1] = 1
    nodes = sorted({n for n, _ in need})
    c.confirm(f"NONLINMOTDISP: {len(nodes)} nodes; NOUT +{sum(len(v) for v in add_nout.values())} node/DOF requests, "
              f"RDND +{added_rd} nodes (MOTION Save Complex TF must be on for RELDISP)")


@command("PLIST", tier=TIER, max_args=3)
def cmd_plist(c):
    """PLIST,[start],[end],[stride]: list / check the wall panels: group, elements, corners, width and height,
    coplanarity, material, BBC and the supported options."""
    m = c.model
    allp = dict(panels(m))
    if not allp:
        c.info("no panels defined")
        return
    ids = _range(c, list(allp), "panels", start_default=1)
    xyz = coordinates(m)
    nonext = eduopt(m, "NONEXT") == "1"
    c.info(f"{'panel':>6s} {'group':>6s} {'elem':>5s} {'mat':>5s} {'BBC':>5s} {'disp':>5s} {'force':>6s} "
           f"{'corners BL BR TR TL':>24s} {'L':>9s} {'H':>9s} {'t':>7s}  status")
    for k in ids:
        rec = allp[k]
        gid = rec.integer(2, 0)
        g = m.groups.get(gid)
        disp, force, bbc = rec.integer(4, 1), rec.integer(5, 1), rec.integer(3, 0)
        status = []
        if (force != 1 or disp != 1) and not nonext:
            status.append("Error 126/128")
        if bbc and curve(m, bbc) is None:
            status.append(f"BBC {bbc} missing")
        elif bbc:
            status += curve_problems(m, bbc)
        try:
            geo = panel_geometry(m, gid, xyz)
            th, mid, warns = panel_section(m, gid)
            status += warns + geo.warnings
            c.info(f"{k:>6d} {gid:>6d} {len(g.elements):>5d} {mid:>5d} {bbc:>5d} {disp:>5d} "
                   f"{HY.MODEL_CODES.get(force, str(force)):>6s} {' '.join(str(v) for v in geo.corners):>24s} "
                   f"{geo.length:>9.4g} {geo.height:>9.4g} {th:>7.4g}  " + ("; ".join(status) if status else "ok"))
        except PN.PanelError as exc:
            c.info(f"{k:>6d} {gid:>6d}  -- {exc}")


# ======================================================================================
# The NONLINEAR deck (.eql) from the model
# ======================================================================================
def build_eql(model) -> Tuple["object", List[str], List[str]]:
    """``(EqlDeck, errors, warnings)``: the resolved NONLINEAR input (D-NON-01) and the CHECK messages of
    Option NON (Errors 121-128, D-NON-12, EDU-09).  The deck is None when there are errors."""
    from ...modules.nonlinear import OPT_BEAMS, OPT_PANELS, OPT_SPRINGS, EqlDeck
    errors: List[str] = []
    warns: List[str] = []
    m = model
    nonext = eduopt(m, "NONEXT") == "1"
    rec = m.options.record("EQL")
    if rec is None:
        rec = EqlRecord("EQL")
        warns.append("EQL not given: EDF 0.8, no damping cut-off, scale 1, elastic damping not added")
    P, S = panels(m), springs(m)
    if rec.given(2):
        opts = rec.get("nonlinopts")
    else:
        opts = (OPT_PANELS if P else 0) | (OPT_SPRINGS if S else 0)
    if opts & OPT_BEAMS:
        warns.append("nonlinear beams are not available in this version: ignored")
    if m.options.entries("B"):
        warns.append("B records (nonlinear beams) are not usable in this version: ignored")
    if not opts & (OPT_PANELS | OPT_SPRINGS):
        errors.append("no nonlinear element type: EQL <NonLinOpts> 1 (panels) and/or 2 (springs), or P/S records")
    if opts & OPT_PANELS and not P:
        errors.append("Error 121: No Panels Specified")
    if opts & OPT_SPRINGS and not S:
        errors.append("no nonlinear springs specified (S)")
    d = EqlDeck()
    d.params.update(title=m.title, model=m.name, edf=rec.get("disp"), nonlinopts=opts & (OPT_PANELS | OPT_SPRINGS),
                    dampcutoff=rec.get("dampcutoff"), dampscale=rec.get("dampscale"), elasticd=rec.get("elasticd"),
                    gravity=m.gravity, nonext=1 if nonext else 0, maxnode=max(m.nodes) if m.nodes else 0)
    xyz = coordinates(m)
    mat_users = material_users(m)
    sc_users = spring_prop_users(m)
    used_bbc: Dict[int, Tuple[int, bool]] = {}
    if opts & OPT_PANELS:
        mats_seen: Dict[int, int] = {}
        for k, r in P:
            gid, bbc, disp, force = r.integer(2, 0), r.integer(3, 0), r.integer(4, 1), r.integer(5, 1)
            g = m.groups.get(gid)
            if g is None:
                errors.append(f"Error 123: Group {gid} referenced in panel {k} does not exist")
                continue
            if g.type not in PANEL_TYPES:
                errors.append(f"panel {k}: group {gid} is a {g.type_name} group, not SHELL")
                continue
            if force == 2:
                errors.append(f"panel {k}: Force Option 2 (Cheng-Mertz Bending) is not included in this version")
            elif force != 1 and not nonext:
                errors.append(f"Error 126: Force Option {force} not supported for panel {k} in this version")
            if disp != 1 and not nonext:
                errors.append(f"Error 128: Displacement Option {disp} not supported for panel {k}")
            if force == 3 or disp == 2:
                warns.append(f"panel {k}: experimental option (EDUOPT,NONEXT,1): Force Opt {force}, Disp. Opt {disp}")
            th, mid, swarn = panel_section(m, gid)
            warns += [f"panel {k}: {w}" for w in swarn if "thickness" in w]
            if len({int(e.mat) for e in g.elements.values()}) > 1:
                errors.append(f"panel {k}: group {gid} has several materials: each panel needs one material of its own "
                              "(GROUPMAT, D-NON-12)")
                continue
            mat = m.materials.get(mid)
            if mat is None:
                errors.append(f"Error 122: Material {mid} referenced in panel {k} does not exist")
                continue
            others = mat_users.get(mid, set()) - {gid}
            if others or mid in mats_seen:
                errors.append(f"panel {k}: material {mid} is also used by group(s) {sorted(others | ({mats_seen[mid]} if mid in mats_seen else set()))}"
                              ": each panel needs a material of its own (GROUPMAT, D-NON-12)")
                continue
            mats_seen[mid] = gid
            try:
                geo = panel_geometry(m, gid, xyz)
            except PN.PanelError as exc:
                errors.append(f"panel {k}: {exc}")
                continue
            warns += [f"panel {k}: {w}" for w in geo.warnings]
            try:
                ec = mat.constants(m.gravity)
            except ValueError as exc:
                errors.append(f"panel {k}: material {mid}: {exc}")
                continue
            if bbc < 1:
                errors.append(f"panel {k}: no backbone curve (P <bbc> = 0; BBCGEN or BBCX/BBCY)")
                continue
            used_bbc[bbc] = (used_bbc.get(bbc, (0, False))[0] + 1, used_bbc.get(bbc, (0, False))[1] or force in (1, 3))
            area = th * geo.length
            inertia = th * geo.length ** 3 / 12.0
            d.add("panels", num=k, group=gid, bbc=bbc, disp=disp, force=force, mat=mid, mtype=mat.mtype,
                  val1=mat.val1, val2=mat.val2, weight=mat.weight, pdamp=mat.pdamp, sdamp=mat.sdamp, e_el=ec.E,
                  nu_el=ec.nu, g_el=ec.G, thick=th, length=geo.length, height=geo.height, area=area, inertia=inertia,
                  bl=geo.corners[0], br=geo.corners[1], tr=geo.corners[2], tl=geo.corners[3], ex=float(geo.e_h[0]),
                  ey=float(geo.e_h[1]), nelem=len(g.elements))
            if mat.pdamp != mat.sdamp:
                warns.append(f"panel {k}: material {mid} has pdamp != sdamp: the S-wave damping {mat.sdamp:g} is the "
                             "elastic damping")
            cv = curve(m, bbc)
            if cv and cv["x"] and cv["y"]:
                ref = ec.G * area if disp == 1 else ec.E * inertia
                kb = cv["y"][0] / cv["x"][0]
                if ref > 0 and abs(kb / ref - 1.0) > 0.01:
                    warns.append(f"EDU-09: panel {k}: BBC {bbc} initial slope {kb:.6g} differs by {100 * (kb / ref - 1):.2f} % "
                                 f"from {'G A_shear' if disp == 1 else 'E I'} = {ref:.6g} (D-NON-11)")
    if opts & OPT_SPRINGS:
        for k, r in S:
            gid, eid, bbc, dof, force = (r.integer(2, 0), r.integer(3, 0), r.integer(4, 0), r.integer(5, 1),
                                         r.integer(6, 4))
            g = m.groups.get(gid)
            if g is None or g.type != SPRING:
                errors.append(f"spring {k}: group {gid} is not a SPRING group")
                continue
            e = g.elements.get(eid)
            if e is None:
                errors.append(f"spring {k}: group {gid} has no element {eid}")
                continue
            if force != 4 and not nonext:
                errors.append(f"Error 127: Force Option {force} not supported for spring {k}")
            if dof not in (1, 2, 3):
                errors.append(f"spring {k}: DOF {dof} must be 1 (X), 2 (Y) or 3 (Z)")
                continue
            sp = m.springs.get(int(e.prop))
            if sp is None:
                errors.append(f"spring {k}: SC property {e.prop} is not defined")
                continue
            users = sc_users.get(int(e.prop), [])
            if len(users) > 1:
                errors.append(f"spring {k}: SC property {e.prop} is shared by {len(users)} spring elements: each "
                              "nonlinear spring needs a property of its own")
                continue
            kel = sp.k[dof - 1]
            if kel <= 0:
                errors.append(f"spring {k}: the SC constant of DOF {'XYZ'[dof - 1]} is {kel:g} (must be > 0)")
                continue
            if any(v != 0 for j, v in enumerate(sp.k) if j != dof - 1):
                warns.append(f"spring {k}: SC property {e.prop} has other non-zero constants: they keep their values but "
                             "share the updated damping (split nonlinear springs into 1D springs, manual 1.5.4)")
            nd = PN.element_nodes(e.nodes)
            if len(nd) < 2:
                errors.append(f"spring {k}: element {eid} has fewer than 2 nodes")
                continue
            if bbc < 1:
                errors.append(f"spring {k}: no backbone curve")
                continue
            used_bbc[bbc] = (used_bbc.get(bbc, (0, False))[0] + 1, used_bbc.get(bbc, (0, False))[1] or force in (1, 3))
            d.add("springs", num=k, group=gid, elem=eid, bbc=bbc, disp=dof, force=force, prop=int(e.prop),
                  node_i=nd[0], node_j=nd[1], scx=sp.scx, scy=sp.scy, scz=sp.scz, scxx=sp.scxx, scyy=sp.scyy,
                  sczz=sp.sczz, damp=sp.damp, k_el=kel)
            cv = curve(m, bbc)
            if cv and cv["x"] and cv["y"]:
                kb = cv["y"][0] / cv["x"][0]
                if abs(kb / kel - 1.0) > 0.01:
                    warns.append(f"EDU-09: spring {k}: BBC {bbc} initial slope {kb:.6g} differs by "
                                 f"{100 * (kb / kel - 1):.2f} % from the spring constant {kel:.6g} (D-NON-11)")
    for b, (_, strict) in sorted(used_bbc.items()):
        probs = curve_problems(m, b, strict=strict)
        if probs:
            errors += probs
            continue
        cv = curve(m, b)
        for p, (x, y) in enumerate(zip(cv["x"], cv["y"]), start=1):
            d.add("bbc", **{"bbc": b, "type": int(cv["type"]), "yield": int(cv["yield_"]), "point": p, "x": float(x),
                            "y": float(y)})
        bbx = HY.Backbone(cv["x"], cv["y"], cv["yield_"])
        if bbx.max_slope_ratio() > 1.0 + 1e-9:
            warns.append(f"BBC {b}: a segment is steeper than the initial slope (the CMS rules cap stiffnesses at Y1/X1)")
    return (None if errors else d), errors, warns


def write_eql_deck(model, path: Path):
    """Write ``<model>.eql``; returns (path or None, errors, warnings)."""
    from ...modules.nonlinear import write_eql
    deck, errors, warns = build_eql(model)
    if deck is None:
        return None, errors, warns
    return write_eql(path, deck, comments=["NONLINEAR input deck written by RUNNONLINEAR from the EQL, P, S and BBC "
                                           "data of the model (D-NON-01)"]), errors, warns


# ======================================================================================
# Runs
# ======================================================================================
def _model_for_run(c):
    num = c.int(1, default=-1)
    if num is None or num < 0:
        num = c.interp.active_model
    if num not in c.interp.models:
        c.fail(f"model {num} is not in memory")
    m = c.interp.models[num]
    if not m.name or not m.path:
        c.fail("Model name/path not defined -- use MDL")
    return m, Path(m.path)


@command("RUNNONLINEAR", tier=TIER, max_args=1)
def cmd_runnonlinear(c):
    """RUNNONLINEAR,[model]: write <model>.eql from the EQL / P / S / BBC data (Errors 121-128 stop the run) and
    run the NONLINEAR module in the model directory (it reads <model>.hou and the .THD files of RELDISP)."""
    from ...modules.nonlinear import read_convergence, run_nonlinear
    m, mdir = _model_for_run(c)
    if not mdir.is_dir():
        c.fail(f"model directory {mdir} does not exist")
    path, errors, warns = write_eql_deck(m, mdir / f"{m.name}.eql")
    for w in warns:
        c.warn(w)
    if errors:
        for e in errors:
            c.error(e)
        raise CommandReported("NONLINEAR input errors")
    if not (mdir / f"{m.name}.hou").exists():
        c.fail(f"{m.name}.hou not found -- run AFWRITE with HOUSE enabled")
    c.info(f"RUNNONLINEAR: model {m.name} in {mdir} ({path.name} written)")
    rc = run_nonlinear(m.name, mdir, echo=lambda line: c.info(line))
    if rc != 0:
        c.fail(f"NONLINEAR finished with status FAILED ({rc}); see {m.name}_NONLINEAR.out")
    rows = read_convergence(mdir)
    if rows:
        r = rows[-1]
        c.info(f"NONLINEAR: SSI analysis {r['iteration']}: max |dE/E| {r['max_de_pct']:.3f} %, max |d xi| "
               f"{r['max_dxi_pct']:.3f} %" + (" -- CONVERGED" if r["converged"] else ""))
    c.confirm(f"NONLINEAR finished (listing {m.name}_NONLINEAR.out)")


@command("COMBXYZTHD", tier=TIER, max_args=1)
def cmd_combxyzthd(c):
    """COMBXYZTHD,<inpfile>: COMB_XYZ_THD (D-NON-08) -- for every line 'outfile fileX fileY fileZ' of <inpfile>
    (default COMB_XYZ_THD.inp; line 1 = number of lines) the output is the time-step-wise sum of the X, Y and
    Z relative-displacement histories (files in the model directory; '-' = no file)."""
    from ...modules.base import ModuleError
    from ...modules.nonlinear import COMB_INP, comb_xyz_thd
    d = model_dir(c)
    try:
        written, warns = comb_xyz_thd(d, c.str(1) or COMB_INP)
    except ModuleError as exc:
        c.fail(str(exc))
    for w in warns:
        c.warn(w)
    c.confirm(f"COMB_XYZ_THD: {len(written)} combined histories written")


def _thd_names(model) -> List[str]:
    need, _ = element_nodes_needed(model)
    mx = max(model.nodes) if model.nodes else 0
    return [nodal_result_name(n, d, "THD", mx) for n, d in need]


@command("NONLINTHD", tier=TIER, max_args=1)
def cmd_nonlinthd(c):
    """NONLINTHD,<prefix>: copy the .THD files NONLINEAR reads (panel corners X/Y/Z, spring ends) to
    <prefix>_<name> (e.g. X_00011TR_X.THD after the X-input RELDISP run) for COMB_XYZ_THD."""
    pre = c.str(1)
    if not pre:
        c.fail("<prefix> required (e.g. X, Y or Z)")
    d = model_dir(c)
    names = _thd_names(c.model)
    if not names:
        c.fail("no panel corner or spring node (define P / S)")
    missing = [n for n in names if not (d / n).exists()]
    if missing:
        c.fail(f"{len(missing)} .THD files missing (run RELDISP with the NONLINMOTDISP requests): "
               + ", ".join(missing[:8]) + (" ..." if len(missing) > 8 else ""))
    for n in names:
        shutil.copyfile(d / n, d / f"{pre}_{n}")
    c.confirm(f"NONLINTHD: {len(names)} .THD files copied to {pre}_*")


def comb_inp_text(model, prefixes: Sequence[str] = ("X", "Y", "Z")) -> str:
    names = _thd_names(model)
    lines = [f"{len(names)}"]
    for n in names:
        lines.append(" ".join([n] + [f"{p}_{n}" for p in prefixes]))
    return "\n".join(lines) + "\n"


@command("NONLINRESET", tier=TIER, max_args=0)
def cmd_nonlinreset(c):
    """NONLINRESET: delete PANEL.NON, SPRING.NON, the *_EQL_Matl_Prop.txt files and the convergence history in
    the model directory -- the next NONLINEAR run is an elastic run (re-run AFWRITE for the elastic HOUSE deck)."""
    from ...modules.nonlinear import CONV_FILE, NON_FILES, PROP_FILES
    d = model_dir(c)
    if not d.is_dir():
        c.fail(f"model directory {d} does not exist")
    targets = {n.lower() for n in list(NON_FILES.values()) + list(PROP_FILES.values()) + [CONV_FILE]}
    gone = []
    for p in sorted(d.iterdir()):
        if p.is_file() and p.name.lower() in targets:
            p.unlink()
            gone.append(p.name)
    c.confirm("NONLINRESET: " + (", ".join(gone) + " deleted" if gone else "no state file in the model directory")
              + "; the next NONLINEAR run is an elastic run")


@command("NONLINSAVE", tier=TIER, max_args=2)
def cmd_nonlinsave(c):
    """NONLINSAVE,<tag>,[<file8>]: copies of the NONLINEAR results under the Fig. 1.2 names -- PanelNNNN_<tag>
    .thd/.ths/.crv, Panel_<tag>.fmu, Panel_EQL_Matl_Prop_<tag>.txt (SPRING likewise), <model>_<tag>.hou;
    <file8> = 1 also copies FILE8 (FILE8X/Y/Z) to FILE8..._<tag> (D-NON-07)."""
    tag = c.str(1)
    if not tag:
        c.fail("<tag> required (e.g. elastic, It#)")
    d = model_dir(c)
    m = c.model
    pats = []
    for p in sorted(d.iterdir()):
        if not p.is_file():
            continue
        nm = p.name
        low = nm.lower()
        stem, ext = os.path.splitext(nm)
        if ((low.startswith("panel") or low.startswith("spring")) and ext.lower() in (".thd", ".ths", ".crv", ".fmu")
                and not stem.lower().endswith(f"_{tag.lower()}") and "_it" not in stem.lower()
                and "_elastic" not in stem.lower()):
            pats.append((p, f"{stem}_{tag}{ext}"))
        elif low in ("panel_eql_matl_prop.txt", "spring_eql_matl_prop.txt"):
            pats.append((p, f"{stem}_{tag}{ext}"))
        elif nm == f"{m.name}.hou":
            pats.append((p, f"{m.name}_{tag}.hou"))
        elif c.int(2, default=0) == 1 and nm in ("FILE8", "FILE8X", "FILE8Y", "FILE8Z"):
            pats.append((p, f"{nm}_{tag}"))
    for src, dst in pats:
        shutil.copy2(src, d / dst)
    c.confirm(f"NONLINSAVE: {len(pats)} files copied with the tag {tag}")


@command("NONLINITER", tier=TIER, max_args=2)
def cmd_nonliniter(c):
    """NONLINITER,<var>[+<var2>...],[<maxit>]: run the command lines stored in the variable(s) (one command per
    item; quote items with commas) once per pass, '#' = pass number, until NONLINEAR reports convergence
    (D-NON-06: |dE/E| < 2 %, |d xi| < 0.5 %) or <maxit> passes (default 10).  Every pass must run
    RUNNONLINEAR.  Nothing runs when the last NONLINEAR run in the model directory already converged."""
    from ...modules.nonlinear import read_convergence
    d = model_dir(c)
    names = [t.strip().lstrip("@") for t in c.str(1).split("+") if t.strip()]
    if not names:
        c.fail("variable name required (VAR,<name>,<command 1>,<command 2>,...)")
    body: List[str] = []
    for name in names:
        v = c.interp.variables.get(name.lower())
        if v is None:
            c.fail(f"variable {name} not defined")
        body += list(v.items)
    if not body:
        c.fail("the variable(s) hold no commands")
    maxit = c.int(2, default=HY.MAX_ITERATIONS)
    if maxit < 1:
        c.fail("<maxit> must be >= 1")
    loops = c.interp.loops
    if any(fr.var in {n.lower() for n in names} for fr in loops):
        c.fail("NONLINITER inside a loop over the same variable")
    rows = read_convergence(d)
    if rows and rows[-1]["converged"]:
        c.info(f"NONLINITER: already converged (SSI analysis {rows[-1]['iteration']}); nothing run (NONLINRESET starts "
               "a new analysis)")
        return
    frame = LoopFrame(names[0].lower(), 0)
    loops.append(frame)
    last = None
    try:
        for k in range(1, maxit + 1):
            frame.index = k
            n0 = len(read_convergence(d))
            for item in body:
                if not c.interp.execute(item):
                    c.fail(f"pass {k}: '{item}' failed; iterations stopped")
            rows = read_convergence(d)
            if len(rows) <= n0:
                c.fail(f"pass {k}: NONLINEAR did not run (the body must run RUNNONLINEAR)")
            last = rows[-1]
            c.info(f"NONLINITER pass {k}: SSI analysis {last['iteration']}: max |dE/E| {last['max_de_pct']:.3f} %, max "
                   f"|d xi| {last['max_dxi_pct']:.3f} %" + (" -- converged" if last["converged"] else ""))
            if last["converged"]:
                break
    finally:
        loops.pop()
    if last is not None and last["converged"]:
        c.confirm(f"NONLINITER: converged after {k} pass(es)")
    else:
        c.warn(f"not converged after {maxit} pass(es) (D-NON-06 limit {HY.MAX_ITERATIONS})")


# ======================================================================================
# NONLINBAT
# ======================================================================================
def nonlinbat_text(model, sel: int) -> str:
    """The generic Option NON batch script of NONLINBAT (a .pre file run with INP)."""
    name = model.name or "<model>"
    an = model.options.record("ANALYS")
    an_toks = an.to_tokens() if an is not None else []
    while len(an_toks) < 12:
        an_toks.append("")

    def analys(mode: int, save: int, simul: str) -> str:
        t = list(an_toks)
        t[2], t[3] = str(mode), str(save)
        if len(t) < 12:
            t += [""] * (12 - len(t))
        t[11] = simul
        return "ANALYS," + ",".join(t).rstrip(",")

    three = sel == 1
    L = ["*" * 78,
         f"* Option NON batch script for model {name}, written by NONLINBAT,{sel} (manual 9.17.21, D-NON-07)",
         "* Run it with INP after the model, its options and the SITE runs (FILE1, or FILE1X/Y/Z for X/Y/Z input)",
         "* are ready.  Edit it as needed (control motions per direction, MOTION scaling of the vertical input).",
         "* Prerequisites: NONLINMOTDISP before AFWRITE (panel corners / spring ends in NOUT and RDND), MOTION",
         "* Save Complex TF on, RELDISP with the free-field reference (RELFILE blank).",
         "*" * 78,
         "* ---- initial elastic SSI analysis (Fig. 1.2 step 1): restart files saved for the iterations",
         "NONLINRESET",
         analys(0, 1, "1" if three else (an_toks[11] or "0")),
         "AOPT,0,0,0,0,1,1,0,0,1,0,0,0,0,0",
         "AFWRITE",
         "RUNPOINT",
         "RUNHOUSE",
         "RUNANALYS"]
    # AFWRITE in the loop writes the MOTION and RELDISP decks only (AOPT first in each direction's variable,
    # so that the ANALYS deck written before NONLINITER is still current when RUNANALYS runs)
    mr = '"AOPT,0,0,0,0,0,0,0,0,0,0,1,0,1,0"'
    if three:
        L += ["* MOTION + RELDISP per input direction; the THD files are kept as X_*, Y_*, Z_* for COMB_XYZ_THD",
              f'VAR,NLDX,{mr},"EDUOPT,TFFILE,FILE8X",AFWRITE,RUNMOTION,RUNRELDISP,"NONLINTHD,X"',
              f'VAR,NLDY,{mr},"EDUOPT,TFFILE,FILE8Y",AFWRITE,RUNMOTION,RUNRELDISP,"NONLINTHD,Y"',
              f'VAR,NLDZ,{mr},"EDUOPT,TFFILE,FILE8Z",AFWRITE,RUNMOTION,RUNRELDISP,"NONLINTHD,Z"',
              'VAR,NLNON,"COMBXYZTHD,COMB_XYZ_THD.inp",RUNNONLINEAR',
              "FOREACH,NLDX,@NLDX[#]",
              "FOREACH,NLDY,@NLDY[#]",
              "FOREACH,NLDZ,@NLDZ[#]",
              "FOREACH,NLNON,@NLNON[#]"]
        body = "NLHOUSE+NLDX+NLDY+NLDZ+NLNON"
    else:
        L += ["* MOTION + RELDISP for the single input direction",
              f"VAR,NLD,{mr},AFWRITE,RUNMOTION,RUNRELDISP",
              "VAR,NLNON,RUNNONLINEAR",
              "FOREACH,NLD,@NLD[#]",
              "FOREACH,NLNON,@NLNON[#]"]
        body = "NLHOUSE+NLD+NLNON"
    L += ["* keep the elastic results (Fig. 1.2 step 3)",
          "NONLINSAVE,elastic",
          "* ---- iterations: the new HOUSE deck, HOUSE, ANALYS New Structure (mode 1), MOTION, RELDISP,",
          "*      COMB_XYZ_THD, NONLINEAR; stop when converged (D-NON-06) or after 10 passes",
          analys(1, 0, "1" if three else (an_toks[11] or "0")),
          "AOPT,0,0,0,0,0,0,0,0,1,0,0,0,0,0",
          "AFWRITE",
          f'VAR,NLHOUSE,"FCOPY,{name}_new.hou,{name}.hou",RUNHOUSE,RUNANALYS',
          f"NONLINITER,{body},10",
          "NONLINSAVE,final",
          "* the converged properties are in Panel_EQL_Matl_Prop.txt / SPRING_EQL_Matl_Prop.txt; FILE8 (FILE8X/Y/Z)",
          "* and the MOTION / RELDISP results are those of the converged model"]
    return "\n".join(L) + "\n"


@command("NONLINBAT", tier=TIER, max_args=1)
def cmd_nonlinbat(c):
    """NONLINBAT,<Sel>: write <model>_NONLINBAT.pre, a generic batch script of the whole Option NON analysis
    (elastic SSI, then iterations with NONLINITER): Sel 0 one input direction, 1 X/Y/Z input with
    COMB_XYZ_THD (a COMB_XYZ_THD.inp listing the THD files is written when the file does not exist)."""
    sel = c.int(1, default=0)
    if sel not in (0, 1):
        c.fail("<Sel> must be 0 (one direction) or 1 (three directions with COMB_XYZ_THD)")
    m = c.model
    if not m.name:
        c.fail("model name not defined -- use MDL")
    d = model_dir(c)
    if not d.is_dir():
        c.fail(f"model directory {d} does not exist")
    p = d / f"{m.name}_NONLINBAT.pre"
    p.write_text(nonlinbat_text(m, sel), encoding="utf-8")
    msg = f"NONLINBAT: {p.name} written"
    if sel == 1:
        q = d / "COMB_XYZ_THD.inp"
        if q.exists():
            msg += "; COMB_XYZ_THD.inp exists (kept)"
        else:
            q.write_text(comb_inp_text(m), encoding="utf-8")
            msg += "; COMB_XYZ_THD.inp written"
    c.confirm(msg)
