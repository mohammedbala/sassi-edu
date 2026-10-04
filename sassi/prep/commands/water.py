"""Water modelling commands FILLPOOL, REFINEMODEL and LISTPOOLINTER (ACS SASSI manual section 9.16;
requirements section 3.4.N; spec 11 section 1; decisions D-WAT-01, D-CHK-07) and the SASSI-EDU
extensions MERGEPOOL and POOLDATA.  The algorithms are in :mod:`sassi.prep.water_lib`; the physics,
the decisions and the limits are explained there and in docs/user/WATER.md.

Workflow of the manual (section 9.16) with the SASSI-EDU commands (``*`` lines are comments)::

    * model 1 is the whole structure; sub-model 2 gets only the walls and floor of one pool
    ACTM,1
    CUTVOL,1,<xmin>,<xmax>,<ymin>,<ymax>,<zmin>,<zmax>
    CUT2SUB,1,2
    ACTM,2
    * optional: a finer wall/floor mesh gives a finer water mesh
    REFINEMODEL
    FILLPOOL,1E9,0,1,-1,<last node of model 1>,0
    ACTM,1
    * the interface nodes of the pool must be found in model 1 (same number and position)
    LISTPOOLINTER,2
    * SASSI-EDU: import the water and the springs into model 1
    MERGEPOOL,2

The water is modelled with SOLID elements of water (bulk modulus 2.2 GPa, shear modulus 1e-8 of it,
0.5 % damping): a nearly incompressible solid without shear stiffness behaves like an inviscid fluid,
so the water that is pushed by the walls acts as the **impulsive hydrodynamic mass** of Westergaard and
Housner (VP-W1).  **Sloshing (the convective mass) is not represented**, and the water needs the
incompatible-mode SOLIDs (``MOPT,0``), which FILLPOOL and MERGEPOOL switch on.
"""
from __future__ import annotations

from typing import Callable

from ...model import fmt_num, make_record
from .. import water_lib as wl
from .. import writer
from ..generation_lib import GenerationError, Notes
from ..registry import command


def _flush(c, notes: Notes) -> None:
    for w in notes.warnings:
        c.warn(w)
    for t in notes.infos:
        c.info(t)


def _run(c, fn: Callable, *args, **kw):
    """Call a water_lib function; a GenerationError aborts the command (model unchanged)."""
    notes = Notes()
    try:
        out = fn(*args, notes=notes, **kw)
    except GenerationError as exc:
        _flush(c, notes)
        c.fail(str(exc))
    _flush(c, notes)
    return out


def _pool_model(c, k: int = 1):
    n = c.int(k, required=True, what="Pool")
    if n not in c.interp.models:
        c.fail(f"model {n} is not in memory")
    if n == c.interp.active_model:
        c.fail(f"model {n} is the active model: activate the original model (ACTM) and give the pool model number")
    return n, c.interp.models[n]


def _range_text(ids) -> str:
    ids = sorted(ids)
    if not ids:
        return "none"
    return f"{ids[0]}-{ids[-1]}" if len(ids) > 1 else str(ids[0])


# ======================================================================================
# FILLPOOL
# ======================================================================================
@command("FILLPOOL", max_args=6)
def cmd_fillpool(c):
    """FILLPOOL,<Stiff>,<Sensitivity>,<EmptyLevels>,<ShellArea>,<offset>,<stiff2>: fill the pool sub-model with water SOLIDs, wall springs and optional area shells.

    The active model must hold only the walls and floor of one pool (SHELL/TSHELL or SOLID elements).
    Arguments (manual 9.16.1; defaults 1e6, 0, 0, -1, -1, 0):

    * ``Stiff`` (> 0): interface spring stiffness along the wall normal (force/length);
    * ``Sensitivity`` (>= 0): allowed Z variation of the nodes of one Z-level (as EXCAV ``delta``);
    * ``EmptyLevels`` (>= 0): Z-levels, counted down from the highest, left without water (freeboard);
    * ``ShellArea``: 1 adds a group of dummy SHELLs on the wetted wall faces (areas for ANSYS processing);
    * ``offset``: water nodes are numbered from offset + 1; <= 0 uses the largest node number of the pool
      model, a positive value below it is an error -- give the last node number of the original model to
      import the water back (MERGEPOOL);
    * ``stiff2`` (>= 0): interface stiffness along the wall (0: frictionless water; > 0 for a pool filled
      with a material that carries shear).

    Fill (the EXCAV algorithm): the floor (lowest horizontal shells or free upward SOLID faces) is the plan
    template; every template cell is extruded between the Z-levels of the wall nodes from the floor up to
    the level ``EmptyLevels`` below the top (hexahedra; prisms for triangles).  The new groups are appended:
    the water SOLIDs (ETYPE 1, a new type-3 M entry: Vp = sqrt(K/rho) with K = 2.2 GPa, Vs = 1e-4 Vp, unit
    weight 9.81 kN/m3 or 0.0624 kcf by the HOUSE gravity, damping 0.5 %), the interface SPRINGs (wall node
    I, water node J, zero length: ``Stiff`` along the wall normals, ``stiff2`` along the wall, no rotation
    stiffness, no damping; oblique walls use the diagonal of the coupled matrix) and the area SHELLs.  The
    rotations of the water nodes and of SOLID-wall interface nodes are fixed (the springs add them) and
    MOPT ``<incomp>`` is set to 0.

    Physics.  The water SOLIDs have the bulk modulus of water and a shear modulus 10^8 times smaller: such a
    solid carries only pressure, so it obeys the acoustic wave equation of an inviscid compressible fluid.
    Well below the compression frequencies of the water (c/4H, about 74 Hz for 5 m of water) the water is
    practically incompressible and the walls drive the potential flow of Westergaard (1933) and Housner
    (1963): the impulsive hydrodynamic mass (VP-W1).  Limits: no sloshing / convective mass (the solid has
    no free-surface gravity stiffness), spurious very-low-frequency shear modes (ignore the water response
    below ~0.5 Hz), and the water needs the incompatible modes (MOPT,0) not to lock.  The springs must be
    stiff compared with the water: FILLPOOL prints the frequency of the water mass on the springs, which
    should be far above the frequencies of interest.
    """
    stiff = c.float(1, default=wl.DEFAULT_STIFF)
    sens = c.float(2, default=0.0)
    if sens < 0:
        c.warn("<Sensitivity> must be >= 0: the default 0 is used")
        sens = 0.0
    empty = c.int(3, default=0)
    shell_area = c.int(4, default=-1)
    offset = c.int(5, default=-1)
    stiff2 = c.float(6, default=0.0)
    res = _run(c, wl.fill_pool, c.model, stiff, sens, empty, shell_area, offset, stiff2)
    w = res.water
    lev = ", ".join(fmt_num(round(z, 10)) for z in res.levels)
    c.info(f"Z-levels above the floor: {lev}; water from {fmt_num(round(res.z_floor, 10))} to "
           f"{fmt_num(round(res.z_surface, 10))} ({res.filled} of {len(res.levels) - 1} intervals, "
           f"{len(res.levels) - 1 - res.filled} empty)")
    c.info(f"water material M,{res.water_material}: type 3, Vp = {w.vp:.6g}, Vs = {w.vs:.4g} (K = {w.K:.6g}, "
           f"G = {w.G:.4g} = 1e-8 K), unit weight {fmt_num(w.weight)}, damping {fmt_num(w.damping)} "
           f"({w.label} units; redefine with M,{res.water_material},... for other units)")
    c.info(f"water: volume {res.volume:.6g}, mass {res.mass:.6g}, wetted area {res.wetted_area:.6g}; "
           f"{res.template_points} plan points x {res.filled + 1} levels")
    if res.spring_props:
        props = "; ".join(f"SC {k}: ({fmt_num(a)}, {fmt_num(b)}, {fmt_num(cc)})"
                          for k, (a, b, cc) in sorted(res.spring_props.items()))
        c.info(f"interface spring properties (kx, ky, kz): {props}")
        fx = [f"{d} {f:.4g} Hz" for d, f in res.interface_freq.items() if f > 0]
        if fx:
            c.info("water mass on the interface springs alone: " + ", ".join(fx)
                   + " -- keep these far above the frequencies of interest (increase <Stiff> otherwise)")
    if res.pairs:
        extra = ""
        if res.wall_rotations_fixed:
            extra += f", all rotations at {len(res.wall_rotations_fixed)} SOLID-wall interface nodes"
        if res.wall_drilling_fixed:
            extra += f", the drilling rotation at {len(res.wall_drilling_fixed)} flat SHELL-wall interface nodes"
        c.info(f"rotations fixed at the {len(res.pairs)} attached water nodes{extra} (the springs add rotational "
               f"DOFs without stiffness, EDU-06; FIXROT no longer sees these nodes as solid-only / shell-only)")
    if res.sloshing:
        c.info("not modelled: sloshing (convective mass); rectangular-tank estimate from the plan extents: "
               + "; ".join(f"{d}: first sloshing frequency {f:.3g} Hz, Housner convective mass {100 * mc:.3g} % "
                           f"of the water" for d, (_, f, mc) in sorted(res.sloshing.items())))
    if res.mopt_changed:
        c.warn("MOPT <incomp> set to 0 (incompatible modes included): the nearly incompressible water SOLIDs "
               "lock without them (D-ELM-02); this applies to every structural SOLID of the model")
    shells = f", {res.shells} interface-area shells (group {res.shell_group})" if res.shells else ""
    c.confirm(f"FILLPOOL: {res.elements} water SOLIDs (group {res.water_group}), {len(res.pairs)} interface "
              f"springs (group {res.spring_group}){shells}; water nodes {_range_text(res.water_nodes)}")


# ======================================================================================
# REFINEMODEL
# ======================================================================================
def _remap_cuts(c, maps) -> int:
    """Session-global cuts (D-MDL-13) follow the new element numbers of the refined groups."""
    cuts = c.interp.session.get("cuts")
    if not cuts:
        return 0
    ok = isinstance(cuts, dict) and all(
        isinstance(s, (set, list, tuple)) and all(isinstance(p, tuple) and len(p) == 2 for p in s)
        for s in cuts.values())
    if not ok:
        c.warn("cuts are stored in an unknown form and were not updated: redefine them")
        return 0
    k = 0
    for cut, s in list(cuts.items()):
        new = []
        for g, e in s:
            mp = maps.get(int(g))
            kids = mp.get(int(e)) if mp else None
            if kids is None:
                new.append((g, e))
            else:
                k += 1
                new.extend((int(g), int(x)) for x in kids)
        cuts[cut] = type(s)(dict.fromkeys(new)) if not isinstance(s, set) else set(new)
    return k


@command("REFINEMODEL", max_args=0)
def cmd_refinemodel(c):
    """REFINEMODEL: split every quadrilateral SHELL, TSHELL and PLANE element into 4 and every hexahedral SOLID into 8.

    Manual 9.16.2: the elements are split at the edge midpoints (with the face and body centres the split
    needs, all at the bilinear / trilinear centres, so area and volume are preserved); midpoints and face
    centres are shared between neighbours, so a conforming mesh stays conforming.  Triangles, prisms /
    pyramids (repeated nodes), BEAMS, SPRING and GENERAL elements are unchanged (hanging nodes where an
    unrefined triangle or prism borders a refined element are reported).  New nodes are numbered after the
    last node and inherit an INT code or a D fixity when every parent corner node has it (re-check the
    interaction nodes, e.g. with INTGEN, and refine the TOPL soil layers when new interaction nodes fall
    between layer interfaces).  The elements of a refined group are renumbered 1..n with the children of
    each parent consecutive; EOUT requests and the session cuts follow; loads and masses are not
    redistributed.  Refine a pool's walls before FILLPOOL (the interface springs are not regenerated).
    """
    m = c.model
    res = _run(c, wl.refine_model, m)
    k = _remap_cuts(c, res.groups)
    if k:
        c.info(f"{k} cut entries replaced by their children (cuts are session-global)")
    if any(r.group in res.groups for r in m.eout):
        c.info("EOUT element lists of the refined groups now list the children of the requested elements")
    if res.interaction_added:
        c.info(f"{res.interaction_added} new nodes are interaction nodes (every parent corner was one): "
               f"check the interaction set (INTLIST, INTGEN)")
    groups = ", ".join(str(g) for g in sorted(res.groups))
    c.confirm(f"REFINEMODEL: {res.quads} quadrilaterals split into 4 and {res.hexes} hexahedra into 8 "
              f"(groups {groups}); {len(res.new_nodes)} new nodes {_range_text(res.new_nodes)}; "
              f"{res.untouched} elements unchanged")


# ======================================================================================
# LISTPOOLINTER
# ======================================================================================
@command("LISTPOOLINTER", max_args=1)
def cmd_listpoolinter(c):
    """LISTPOOLINTER,<Pool>: list the interface nodes of FILLPOOL pool model <Pool> that the active (original) model has with the same number and position.

    Manual 9.16.3: the original model and the pool sub-model must both be in memory and the original model
    must be active (ACTM).  The interface nodes are the wall/floor nodes attached to the water by the
    FILLPOOL springs (it works only for a pool filled with FILLPOOL).  Interface nodes missing from the
    original model, or at another position there, are reported as warnings: they are the nodes to fix
    before the water can be imported (MERGEPOOL).
    """
    n, pool = _pool_model(c)
    try:
        res = wl.match_interface(c.model, pool)
    except GenerationError as exc:
        c.fail(f"model {n}: {exc}")
    total = len(res.matched) + len(res.missing) + len(res.moved)
    c.info(f"{'node':>8} {'x':>14} {'y':>14} {'z':>14}")
    for nid, p in res.matched:
        c.info(f"{nid:>8} {p[0]:>14.6g} {p[1]:>14.6g} {p[2]:>14.6g}")
    if res.missing:
        c.warn(f"{len(res.missing)} interface nodes of pool model {n} are not defined in the active model: "
               f"{wl._short(res.missing)}")
    if res.moved:
        c.warn(f"{len(res.moved)} interface nodes have another position in the active model: "
               + ", ".join(f"{a} (distance {d:.4g})" for a, d in res.moved[:10])
               + (" ..." if len(res.moved) > 10 else ""))
    c.confirm(f"LISTPOOLINTER: {len(res.matched)} of {total} interface nodes of pool model {n} have the same "
              f"number and position in model {c.interp.active_model}")


# ======================================================================================
# MERGEPOOL (SASSI-EDU extension)
# ======================================================================================
@command("MERGEPOOL", max_args=1, tier="P2")
def cmd_mergepool(c):
    """MERGEPOOL,<Pool>: (SASSI-EDU) import the water, interface springs and area shells of FILLPOOL pool model <Pool> into the active model.

    Completes the workflow of manual 9.16 ("import the water and spring group back into the original
    model") without renumbering the nodes, which MERGE always does (D-MDL-12): every pool interface node
    must exist in the active model with the same number and position (LISTPOOLINTER), and no water node
    number may be in use (fill the pool with ``offset`` >= the last node number of the original model).
    The FILLPOOL groups are appended after the last group with their element numbers; the water and shell
    materials and the spring properties get new numbers after the last ones; the water nodes keep their
    fixed rotations, interface wall nodes without a rotational DOF in the active model get theirs fixed,
    and MOPT ``<incomp>`` is set to 0 (incompatible modes for the water SOLIDs).
    """
    n, pool = _pool_model(c)
    res = _run(c, wl.merge_pool, c.model, pool)
    if res.wall_rotations_fixed:
        c.info(f"rotations fixed at {len(res.wall_rotations_fixed)} SOLID-wall interface nodes: "
               f"{wl._short(res.wall_rotations_fixed)}")
    if res.wall_drilling_fixed:
        c.info(f"drilling rotation fixed at {len(res.wall_drilling_fixed)} flat SHELL-wall interface nodes: "
               f"{wl._short(res.wall_drilling_fixed)}")
    if res.mopt_changed:
        c.warn("MOPT <incomp> set to 0 (incompatible modes included): needed by the water SOLIDs; this applies "
               "to every structural SOLID of the model")
    gm = ", ".join(f"{a}->{b}" for a, b in sorted(res.groups.items()))
    mm = ", ".join(f"{a}->{b}" for a, b in sorted(res.materials.items()))
    c.confirm(f"MERGEPOOL: pool model {n} imported into model {c.interp.active_model}: groups {gm}; materials "
              f"{mm or 'none'}; {len(res.nodes)} water nodes {_range_text(res.nodes)}")


# ======================================================================================
# POOLDATA (SASSI-EDU storage record written by WRITE)
# ======================================================================================
@command("POOLDATA", cls="record", tier="P2")
def cmd_pooldata(c):
    """POOLDATA,<water>,<springs>,<shells>,<stiff>,<sens>,<empty>,<shellarea>,<offset>,<stiff2>,<watermat>,<zfloor>,<zsurface>: (SASSI-EDU) FILLPOOL pool data.

    Written by WRITE after the groups of a model filled with FILLPOOL (the groups it created, its
    arguments, the water material and the floor and surface elevations) so that LISTPOOLINTER and
    MERGEPOOL work after WRITE -> INP; not meant to be typed.  ``POOLDATA`` without arguments deletes it.
    """
    vals = [c.raw(k) or None for k in range(1, c.nargs + 1)]
    if not vals:
        c.model.options.delete_record("POOLDATA")
        c.confirm("pool data deleted")
        return
    rec = make_record("POOLDATA", vals)
    if len(vals) > len(wl.PoolDataRecord.FIELDS):
        c.warn(f"arguments after argument {len(wl.PoolDataRecord.FIELDS)} ignored")
        rec = make_record("POOLDATA", vals[:len(wl.PoolDataRecord.FIELDS)])
    c.model.options.set_record(rec)
    c.confirm("pool data stored")


def _write_pooldata(m) -> list:
    rec = m.options.record("POOLDATA")
    if rec is None:
        return []
    return ["* FILLPOOL pool data (SASSI-EDU): water, spring and area-shell groups, FILLPOOL arguments, water "
            "material, floor and surface elevations", writer.record_line(rec)]


writer.register_writer("POOLDATA", _write_pooldata)
