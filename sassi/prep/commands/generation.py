"""Model generation and combination commands (manual sections 9.7 and 9.9; requirements section 3.4.I;
spec 09 sections 2 and 4; spec 04 section 7).

ETYPEGEN INTGEN RADIUS GLB2LOC GROUPMAT RMVUNUSED NCOM GCOM WELD MERGE MERGESOIL MERGEGROUP ROTATE
TRANSLATE EXCAV SOILMESH.  The algorithms are in :mod:`sassi.prep.generation_lib`; the handlers
parse the arguments, call them and print their notes.  Every command validates its input first:
an error leaves the models unchanged.

Decisions applied (requirements section 7; spec 09 section 11 open questions):

* ETYPEGEN,0 makes the implicit classification explicit (SOLID/PLANE by D-HOU-01, shells 1);
  ``ETYPEGEN,0,RESET`` restores the implicit ETYPE 0 (spec 09 Q-G4).
* INTGEN (Q-G9): excavation elements = SOLID/PLANE with ETYPE 2, or ETYPE 0 below grade by D-HOU-01
  (with a warning that the manual requires an explicit ETYPE); FSIN = nodes of the boundary faces
  that are not ground-surface faces; FFV internal level j (0 = bottom) is selected when
  ``j mod (skip+1) == 0``; option 4 = SOLID/PLANE/SHELL/TSHELL nodes and BEAMS I/J nodes at grade;
  buried-shell nodes are added by every option; only code 0 (interaction) is set or cleared (D-MDL-07).
* RADIUS: ``r_e = Scale sqrt(A_plan)``, Scale default 0.9 (D-PNT-05).
* MERGE / MERGESOIL: D-MDL-12.  ROTATE: D-MDL-14.  WELD: D-MDL-15.  EXCAV / SOILMESH: D-MDL-19.
* EXCAV with ``delta`` > tol stores ``EDUOPT,GEOMTOL`` = delta (the z scatter of the structure nodes
  about the snapped levels) in the excavation model; MERGESOIL matches the interface with the larger
  of the two models' tolerances and passes that tolerance on to the merged model (D-GEN-06), so the
  excavation is still joined to a basement whose levels scatter (spec 09 section 4.1 step 5).
* MERGESOIL warns about foundation-soil interface nodes of the soil model without a structure
  partner (spec 01 rules 8 and 11); ROTATE warns about D fixities and MT/MR masses whose global-axis
  meaning changes under the rotation (re-run FIXROT after ROTATE).
* GCOM and MERGEGROUP update the session-global cuts (D-MDL-13) when they are stored as
  ``{cut: set of (group, element)}`` in ``interp.session['cuts']``.
"""
from __future__ import annotations

from typing import Callable, Optional, Tuple

import numpy as np

from ...model import SSIModel, fmt_num, is_number
from .. import generation_lib as gl
from ..generation_lib import GenerationError, Notes
from ..registry import command
from . import range_args


def _flush(c, notes: Notes) -> None:
    for w in notes.warnings:
        c.warn(w)
    for t in notes.infos:
        c.info(t)


def _run(c, fn: Callable, *args, **kw):
    """Call a generation function; a GenerationError aborts the command (model unchanged)."""
    notes = Notes()
    try:
        out = fn(*args, notes=notes, **kw)
    except GenerationError as exc:
        _flush(c, notes)
        c.fail(str(exc))
    _flush(c, notes)
    return out


def _model_number(c, k: int, what: str, must_exist: bool = True) -> int:
    n = c.int(k, required=True, what=what)
    if n < 0:
        c.fail(f"<{what}>: model numbers are >= 0")
    if must_exist and n not in c.interp.models:
        c.fail(f"model {n} is not in memory")
    return n


def _store_model(c, number: int, new: SSIModel) -> None:
    """Put a generated model into slot ``number``; an existing model keeps its name and path."""
    old = c.interp.models.get(number)
    if old is not None:
        new.name, new.path = old.name, old.path
        if not new.title:
            new.title = old.title
        if not old.is_empty() and number != c.interp.active_model:
            c.warn(f"model {number} overwritten")
    c.interp.models[number] = new


def _remap_cuts(c, fn: Callable[[int, int], Optional[Tuple[int, int]]]) -> None:
    """Remap the (group, element) pairs of the session-global cuts (D-MDL-13) after GCOM /
    MERGEGROUP; ``fn`` returns the new pair or None to keep it."""
    cuts = c.interp.session.get("cuts")
    if not cuts:
        return
    ok = isinstance(cuts, dict) and all(
        isinstance(s, (set, list, tuple)) and all(isinstance(p, tuple) and len(p) == 2 for p in s)
        for s in cuts.values())
    if not ok:
        c.warn("cuts are stored in an unknown form and were not updated: redefine them")
        return
    k = 0
    for cut, s in list(cuts.items()):
        new = []
        for p in s:
            q = fn(int(p[0]), int(p[1]))
            if q is not None and q != p:
                k += 1
                new.append(q)
            else:
                new.append(p)
        cuts[cut] = type(s)(new)
    if k:
        c.info(f"{k} cut entries updated to the new group/element numbers (cuts are session-global)")


# ======================================================================================
# ETYPEGEN, INTGEN, RADIUS
# ======================================================================================
@command("ETYPEGEN", max_args=2)
def cmd_etypegen(c):
    """ETYPEGEN,<type>,[RESET]: ETYPE of every SOLID/PLANE/SHELL/TSHELL element: 0 by location (explicit
    1/2), 1 structure, 2 excavated soil / buried shell; ``0,RESET`` restores the implicit ETYPE 0."""
    t = c.int(1, required=True, what="type")
    reset = c.given(2)
    if reset and (c.word(2) != "RESET" or t != 0):
        c.fail("argument 2 must be RESET (with <type> 0)")
    counts = _run(c, gl.etypegen, c.model, t, reset=reset)
    if reset:
        c.confirm(f"{counts[0]} elements set to the implicit ETYPE 0")
    else:
        c.confirm(f"ETYPE set: {counts[1]} structure, {counts[2]} excavated soil / buried shell")


@command("INTGEN", max_args=2)
def cmd_intgen(c):
    """INTGEN,<type>,[level skip]: generate interaction nodes: 0 clear, 1 FV, 2 FI-EVBN (MSM),
    3 FI-FSIN (SM), 4 surface, 5 FFV (internal levels every level skip + 1, default skip 1); 1-5 add."""
    m = c.model
    opt = c.int(1, required=True, what="type")
    skip = c.int(2, default=1)
    if c.given(2) and opt != 5:
        c.warn("[level skip] is used by option 5 (FFV) only; ignored")
    if opt == 0:
        before = sum(1 for nd in m.nodes.values() if 0 in nd.flags)
        gl.clear_interaction(m)
        c.confirm(f"INTGEN,0: {before} interaction flags cleared")
        return
    res = _run(c, gl.intgen, m, opt, skip)
    name = gl.INTGEN_METHODS[opt]
    if opt == 5:
        c.info(f"FFV: {len(res.levels)} excavation node levels; internal levels kept (0 = bottom): "
               f"{', '.join(str(j) for j in res.selected_levels) or 'none'} (level skip {skip})")
    if res.buried:
        c.info(f"{res.buried} nodes of buried shells (ETYPE 2) included")
    c.confirm(f"INTGEN,{opt} {name}: {len(res.nodes)} nodes in the set, {len(res.added)} added; "
              f"{res.total} interaction nodes in total")


@command("RADIUS", max_args=2)
def cmd_radius(c):
    """RADIUS,[Scale],[FileName]: POINT central-zone radius r_e = Scale sqrt(A_plan) of every excavation
    element, with min / average / max (D-PNT-05; Scale default 0.9)."""
    m = c.model
    scale = c.float(1, default=0.9)
    if scale <= 0:
        c.fail("<Scale> must be > 0")
    radii = gl.excavation_radii(m, scale)
    if not radii:
        c.fail("no excavation elements (SOLID/PLANE with ETYPE 2 or below grade)")
    r = np.array([x for _, _, x in radii])
    rmin, ravg, rmax = float(r.min()), float(r.mean()), float(r.max())
    fname = c.str(2)
    if fname:
        p = c.output_path(fname)
        lines = [f"# RADIUS: excavation element radii r_e = Scale*sqrt(A_plan) (D-PNT-05), Scale = {fmt_num(scale)}",
                 "# group element r_e"]
        lines += [f"{g} {e} {x:.8g}" for g, e, x in radii]
        lines += [f"MINIMUM {rmin:.8g}", f"AVERAGE {ravg:.8g}", f"MAXIMUM {rmax:.8g}"]
        try:
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        except OSError as exc:
            c.fail(f"cannot write {p}: {exc}")
        c.info(f"radii written to {p}")
    c.info(f"{len(radii)} excavation elements: radius min {rmin:.6g}, average {ravg:.6g}, max {rmax:.6g}")
    c.interp.session["radius"] = {"min": rmin, "average": ravg, "max": rmax, "scale": scale}
    c.confirm(f"average radius {ravg:.6g} (use it as the POINT radius of the central zone)")


# ======================================================================================
# Node and group numbering, coordinate systems
# ======================================================================================
@command("GLB2LOC", max_args=4)
def cmd_glb2loc(c):
    """GLB2LOC,<Start>,<End>,<Stride>,<Sysno>: store nodes in Cartesian system Sysno (0 = global)."""
    m = c.model
    s = c.int(4, required=True, what="Sysno")
    ids = range_args(c, 1, m.nodes, all_default=True)
    try:
        n = gl.glb2loc(m, ids, s)
    except GenerationError as exc:
        c.fail(str(exc))
    c.confirm(f"{n} nodes converted to system {s}")


@command("GROUPMAT", max_args=0)
def cmd_groupmat(c):
    """GROUPMAT: one new material per group (copy of the material of its first element); the original
    material table is replaced (Option NON preparation)."""
    mp = _run(c, gl.groupmat, c.model)
    c.confirm(f"{len(mp)} group materials created")


@command("RMVUNUSED", max_args=0)
def cmd_rmvunused(c):
    """RMVUNUSED: delete the nodes that are used by no element and are not interaction nodes."""
    removed = _run(c, gl.rmvunused, c.model)
    c.confirm(f"{len(removed)} unused nodes removed")


@command("NCOM", max_args=1)
def cmd_ncom(c):
    """NCOM,[MapFile]: renumber nodes 1..N without gaps (order kept) and remap every node reference."""
    m = c.model
    mp = _run(c, gl.ncom, m)
    changed = sum(1 for o, n in mp.items() if o != n)
    if c.given(1):
        p = c.output_path(c.str(1))
        try:
            gl.write_node_map(p, mp, "NCOM node map")
        except OSError as exc:
            c.warn(f"cannot write {p}: {exc}")
        else:
            c.info(f"node map written to {p}")
    c.confirm(f"{len(mp)} nodes numbered 1..{len(mp)} ({changed} renumbered)")


@command("GCOM", max_args=0)
def cmd_gcom(c):
    """GCOM: renumber groups 1..G without gaps (order kept); element numbers unchanged."""
    m = c.model
    mp = gl.gcom(m)
    changed = {o: n for o, n in mp.items() if o != n}
    if changed:
        _remap_cuts(c, lambda g, e: (changed[g], e) if g in changed else None)
    c.confirm(f"{len(mp)} groups numbered 1..{len(mp)} ({len(changed)} renumbered)")


@command("WELD", max_args=1)
def cmd_weld(c):
    """WELD,[FORCE]: merge coincident nodes in the connectivity (lowest number kept); the two ends of a
    zero-length spring are never welded together unless FORCE (other nodes at that point still are)."""
    force = c.given(1)
    if force and c.word(1) != "FORCE":
        c.fail("argument 1 must be FORCE")
    sub = _run(c, gl.weld, c.model, force=force)
    if sub:
        c.info("substituted nodes stay in the node list unattached (RMVUNUSED removes them)")
    c.confirm(f"{len(sub)} coincident nodes substituted")


# ======================================================================================
# Combination: MERGE, MERGESOIL, MERGEGROUP
# ======================================================================================
@command("MERGE", max_args=5)
def cmd_merge(c):
    """MERGE,<Mdl1>,<Mdl2>,<X>,<Y>,<Z>: active model = Mdl1 + Mdl2 (numbers offset, translated by X,Y,Z)."""
    a = _model_number(c, 1, "Mdl1")
    b = _model_number(c, 2, "Mdl2")
    d = [c.float(k, 0.0) for k in (3, 4, 5)]
    out, maps = _run(c, gl.merge_models, c.interp.models[a], c.interp.models[b], d)
    act = c.interp.active_model
    _store_model(c, act, out)
    c.info(f"model {b}: nodes offset by {maps.node_offset}, groups by {maps.group_offset}")
    c.confirm(f"models {a} and {b} merged into model {act}: {len(out.nodes)} nodes, {len(out.groups)} groups")


@command("MERGESOIL", max_args=7)
def cmd_mergesoil(c):
    """MERGESOIL,<Struct>,<Soil>,[Mode],[Stiff],[Stiff2],[SepLevel],[Mapping]: join a structure model and
    an excavation model into the active model (Mode 0 unbonded, 1 merged nodes, 2 stiff springs, 3 stiff
    below / soft above SepLevel).  A trailing non-numeric argument is the Mapping file (D-MDL-12)."""
    a = _model_number(c, 1, "Struct")
    b = _model_number(c, 2, "Soil")
    if a == b:
        c.fail("<Struct> and <Soil> must be different models")
    n = c.nargs
    mapping, last = None, n
    if n == 7:
        mapping, last = c.str(7), 6
    elif n >= 3 and c.given(n) and not is_number(c.raw(n)):
        # D-MDL-12 / spec 09 section 4.5: a non-numeric token in any trailing optional position,
        # [Mode] included (MERGESOIL,1,2,file), is the Mapping file
        mapping, last = c.str(n), n - 1
        c.info(f"argument {n} '{mapping}' taken as the Mapping file (canonical position 7)")
    mode = c.int(3, default=1) if last >= 3 else 1
    stiff = c.float(4, default=1.0e7) if last >= 4 else 1.0e7
    stiff2 = c.float(5, default=10.0) if last >= 5 else 10.0
    sep = c.float(6) if last >= 6 else None
    res = _run(c, gl.mergesoil, c.interp.models[a], c.interp.models[b], mode, stiff, stiff2, sep)
    act = c.interp.active_model
    _store_model(c, act, res.model)
    np_ = len(res.pairs)
    if mode == 1:
        c.info(f"{np_} interface node pairs merged; the {np_} higher-numbered soil nodes are unused (RMVUNUSED)")
    elif mode in (2, 3):
        c.info(f"{res.springs} interface springs in group {res.spring_group}")
    else:
        c.info(f"{np_} coincident interface node pairs left unbonded")
    if mapping:
        p = c.output_path(mapping)
        try:
            p.parent.mkdir(parents=True, exist_ok=True)
            gl.write_node_map(p, res.node_map, f"MERGESOIL excavation node map: soil model {b} node -> model {act}")
        except OSError as exc:
            c.warn(f"cannot write {p}: {exc}")
        else:
            c.info(f"mapping file written to {p}")
    c.confirm(f"models {a} (structure) and {b} (soil) merged into model {act}, mode {mode}: {np_} interface pairs")


@command("MERGEGROUP", max_args=11)
def cmd_mergegroup(c):
    """MERGEGROUP,<dest>,[G1],...,[G10]: merge same-type groups into dest (then GCOM)."""
    m = c.model
    dest = c.int(1, required=True, what="dest")
    others = [c.int(k) for k in range(2, c.nargs + 1) if c.given(k)]
    if not others:
        c.fail("no group to merge")
    emap = _run(c, gl.mergegroup, m, dest, others)
    if emap:
        _remap_cuts(c, lambda g, e: emap.get((g, e)))
    c.confirm(f"{len(emap)} elements merged into group {dest}; use GCOM to close the group numbering")


# ======================================================================================
# Geometry: ROTATE, TRANSLATE
# ======================================================================================
@command("ROTATE", max_args=6)
def cmd_rotate(c):
    """ROTATE,<x>,<y>,<z>,<rxy>,<ryz>,<rzx>: rotate the model about (x,y,z): rxy about Z, then ryz about X,
    then rzx about Y (degrees, right-hand rule, D-MDL-14)."""
    ctr = [c.float(k, 0.0) for k in (1, 2, 3)]
    ang = [c.float(k, 0.0) for k in (4, 5, 6)]
    k = _run(c, gl.rotate_model, c.model, ctr, *ang)
    c.confirm(f"{k} nodes rotated by rxy {fmt_num(ang[0])}, ryz {fmt_num(ang[1])}, rzx {fmt_num(ang[2])} deg")


@command("TRANSLATE", max_args=3)
def cmd_translate(c):
    """TRANSLATE,<x>,<y>,<z>: move every node (defaults 0)."""
    d = [c.float(k, 0.0) for k in (1, 2, 3)]
    k = _run(c, gl.translate_model, c.model, d)
    c.confirm(f"{k} nodes translated by ({fmt_num(d[0])}, {fmt_num(d[1])}, {fmt_num(d[2])})")


# ======================================================================================
# EXCAV, SOILMESH
# ======================================================================================
@command("EXCAV", max_args=2)
def cmd_excav(c):
    """EXCAV,<model>,[delta]: excavation volume of the active model's basement, stored in <model>."""
    dest = _model_number(c, 1, "model", must_exist=False)
    if dest == c.interp.active_model:
        c.fail("<model> must differ from the active model (the source)")
    delta = c.float(2, default=0.0)
    if delta < 0:
        c.warn("[delta] must be positive: the default 0 is used")
        delta = 0.0
    res = _run(c, gl.excav, c.model, delta)
    _store_model(c, dest, res.model)
    c.info(f"levels: {', '.join(fmt_num(round(z, 10)) for z in res.levels)}")
    c.confirm(f"EXCAV: model {dest} holds {res.elements} excavation SOLIDs ({res.cells} cells x "
              f"{len(res.levels) - 1} layers) and {res.nodes} nodes on {len(res.levels)} levels")


@command("SOILMESH", max_args=10)
def cmd_soilmesh(c):
    """SOILMESH,<dest>,<sX>,<sY>,<hori>,<vert>,<xAdj>,<yAdj>,<Zdepth>,<contact>,<rNum>: near-field soil mesh
    around the basement (soil-pressure models), stored in <dest> numbered after the active model."""
    dest = _model_number(c, 1, "dest", must_exist=False)
    if dest == c.interp.active_model:
        c.fail("<dest> must differ from the active model (the source)")
    sx, sy = c.float(2, 0.0), c.float(3, 0.0)
    hori, vert = c.int(4, 0), c.int(5, 0)
    xadj, yadj, zdepth = c.float(6, 0.0), c.float(7, 0.0), c.float(8, 0.0)
    contact, rnum = c.int(9, 0), c.int(10, 1)
    res = _run(c, gl.soilmesh, c.model, sx, sy, hori, vert, xadj, yadj, zdepth, contact, rnum)
    _store_model(c, dest, res.model)
    extra = f", {res.springs} contact springs (group {res.spring_group})" if res.springs else ""
    c.info("write the soil mesh (WRITE) and read it over the original model (INP); do not MERGE it")
    c.confirm(f"SOILMESH: model {dest} holds {res.elements} soil SOLIDs in group {res.group}, "
              f"{res.new_nodes} new nodes{extra}")
