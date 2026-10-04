"""Cut and submodeling commands (manual section 9.10) and calculation commands (manual section 9.13);
requirements section 3.4.J, spec 09 section 5, spec 10 sections 1.4 and 3, spec 04 section 8.

CUTADD CUTRMV CUTVOL SLICE CUTCLR CUT2SUB CSECT EXTRACTEXCAV SPLITGROUP TRANELEM TRANVOL READSTR
SECDATAOPT CALCPAR CALCMOI CALCC CALCM CALCSECTHIST.  (CALCSECTHISTDB needs the binary stress
database, tier P2: its catalogue placeholder prints the tier message.)

The algorithms are in :mod:`sassi.prep.cuts_lib`; the handlers parse the arguments, call them and
print their notes.  Every command validates its input first: an error leaves the models and the cuts
unchanged.

Typical section-cut workflow (manual 5.8.2, spec 04 section 8.3)::

    INP,model.pre
    READSTR,NSTRESS/ESTRESS_00101.ess      element stresses of one time step (STRESS with SECDATAOPT,1)
    CUTVOL,3,52.5,52.8,-320,320,2.53,45.22  select the elements of one wall
    CSECT,1,3,52.65,0,4,1,0,0              cross-section model in model 1
    ACTM,1
    CALCPAR,1,0,0,0,1,0,10                 area, centroid, inertias and the six resultants
    ACTM,0
    CALCSECTHIST,NSTRESS/ESTRESS.lst,3,52.65,0,4,1,0,0,0,1,0,10,0.01,wall.csv    history of the resultants

CALCPAR / CALCMOI on a model that is not a cross-section give the **base section** of the whole structure
(the plane through its end on the -n side; Q-11, spec 10 section 3.8): ``CALCPAR,0,0,1,1,0,0`` on a building
model gives its base forces and moments.  Any other plane needs a cut and CSECT.  BEAMS are drawn by CSECT
but never belong to the section (no ``.ess`` record); CALCPAR / CALCMOI list them in a warning and CALCC on a
cross-section model leaves them out.

Decisions: D-MDL-13 (session-global cuts, CUTVOL box test, CUT2SUB contents), D-SEC-01 ... D-SEC-06,
D-STR-12 (see :mod:`sassi.prep.cuts_lib`).
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from ...model import SSIModel, fmt_num, make_record
from .. import cuts_lib as cl
from ..check import SHELL
from ..cuts_lib import SectionError
from ..generation_lib import Notes
from ..registry import command
from . import select_ids


# ======================================================================================
# helpers
# ======================================================================================
def _flush(c, notes: Notes) -> None:
    for w in notes.warnings:
        c.warn(w)
    for t in notes.infos:
        c.info(t)


def _run(c, fn: Callable, *args, **kw):
    """Call a cuts_lib function with a Notes collector; a SectionError aborts the command."""
    notes = Notes()
    try:
        out = fn(*args, notes=notes, **kw)
    except SectionError as exc:
        _flush(c, notes)
        c.fail(str(exc))
    _flush(c, notes)
    return out


def _model_number(c, k: int, what: str) -> int:
    n = c.int(k, required=True, what=what)
    if n < 0:
        c.fail(f"<{what}>: model numbers are >= 0")
    return n


def _cut_number(c, k: int = 1) -> int:
    n = c.int(k, required=True, what="cutnum")
    if n < 1:
        c.fail("<cutnum>: cut numbers are >= 1")
    return n


def _store_model(c, number: int, new: SSIModel) -> None:
    """Put a generated model into slot ``number``; an existing model keeps its name, path and title."""
    old = c.interp.models.get(number)
    if old is not None:
        new.name, new.path = old.name, old.path
        if not new.title:
            new.title = old.title
        if not old.is_empty():
            c.warn(f"model {number} overwritten")
    c.interp.models[number] = new


def _vec(c, k: int, what: str, default: float = 0.0) -> np.ndarray:
    return np.array([c.float(k + i, default, what=f"{what}{'xyz'[i]}") for i in range(3)])


def _box(c, k: int) -> List[Optional[float]]:
    return [c.float(k + i) for i in range(6)]


def _elements_arg(c, k: int) -> List[int]:
    """Element numbers of ``<e1>,...,<eN>`` (ids and ``a-b`` ranges, L8) or ``RANGE,<start>,[end],[stride]``
    from argument k (spec 09 section 5.4)."""
    if c.word(k) == "RANGE":
        a = c.int(k + 1, required=True, what="elem start")
        b = c.int(k + 2, default=a)
        s = c.int(k + 3, default=1)
        if s <= 0:
            c.warn("[stride] must be >= 1: 1 used")
            s = 1
        if b < a:
            a, b = b, a
        return list(range(a, b + 1, s))
    ids = c.id_list(k)
    if not ids:
        c.fail("element list (or RANGE,<start>,[end],[stride]) required")
    return ids


def _remap_cuts(c, gid: int, moved: Dict[int, int], new: int) -> int:
    """Cut entries of elements moved by SPLITGROUP follow them (cuts are session-global, D-MDL-13)."""
    k = 0
    for cut, s in cl.cut_store(c.interp).items():
        hit = [(g, e) for g, e in s if g == gid and e in moved]
        for g, e in hit:
            s.discard((g, e))
            s.add((new, moved[e]))
            k += 1
    return k


# ======================================================================================
# Cuts: CUTADD, CUTRMV, CUTVOL, SLICE, CUTCLR (spec 09 sections 5.1, 5.4-5.7, 5.9)
# ======================================================================================
@command("CUTADD")
def cmd_cutadd(c):
    """CUTADD,<cutnum>,<group>,<e1>,...,<eN> | CUTADD,<cutnum>,<group>,RANGE,<start>,[end],[stride]: add
    elements of one group to a cut (created if needed; duplicates ignored)."""
    cut = _cut_number(c)
    gid = c.int(2, required=True, what="group")
    grp = c.model.groups.get(gid)
    if grp is None:
        c.fail(f"group {gid} is not defined in the active model")
    ids = _elements_arg(c, 3)
    have = [e for e in ids if e in grp.elements]
    miss = [e for e in ids if e not in grp.elements]
    if miss:
        c.warn(f"{len(miss)} elements are not in group {gid}: {cl._short(miss)}")
    store = cl.cut_store(c.interp)
    s = store.setdefault(cut, set())
    before = len(s)
    s.update((gid, e) for e in have)
    c.confirm(f"cut {cut}: {len(s) - before} elements added ({len(s)} in total)")


@command("CUTRMV")
def cmd_cutrmv(c):
    """CUTRMV,<cutnum>,<group>,<e1>,...,<eN> | CUTRMV,<cutnum>,<group>,RANGE,<start>,[end],[stride]: remove
    elements from a cut."""
    cut = _cut_number(c)
    gid = c.int(2, required=True, what="group")
    ids = _elements_arg(c, 3)
    try:
        s = cl.get_cut(c.interp, cut)
    except SectionError as exc:
        c.fail(str(exc))
    rm = [(gid, e) for e in ids if (gid, e) in s]
    for k in rm:
        s.discard(k)
    if len(rm) < len(ids):
        c.info(f"{len(ids) - len(rm)} of the listed elements were not in cut {cut}")
    c.confirm(f"cut {cut}: {len(rm)} elements removed ({len(s)} left)")


@command("CUTVOL", max_args=7)
def cmd_cutvol(c):
    """CUTVOL,<cutnum>,[Xmin],[Xmax],[Ymin],[Ymax],[Zmin],[Zmax]: add the elements whose nodes all lie in
    the closed box (blank bounds = the active model's extents; D-MDL-13)."""
    cut = _cut_number(c)
    m = c.model
    v = cl.model_view(m)
    try:
        b = cl.box_bounds(v, _box(c, 2))
    except SectionError as exc:
        c.fail(str(exc))
    keys = cl.elements_in_box(m, b, v)
    s = cl.cut_store(c.interp).setdefault(cut, set())
    before = len(s)
    s.update(keys)
    c.confirm(f"cut {cut}: {len(keys)} elements in the box, {len(s) - before} added ({len(s)} in total)")


@command("SLICE", max_args=7)
def cmd_slice(c):
    """SLICE,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>: add the elements that cross (or touch) the plane."""
    cut = _cut_number(c)
    P = _vec(c, 2, "p")
    n = _vec(c, 5, "n")
    try:
        keys = cl.elements_on_plane(c.model, P, n)
    except SectionError as exc:
        c.fail(str(exc))
    s = cl.cut_store(c.interp).setdefault(cut, set())
    before = len(s)
    s.update(keys)
    c.confirm(f"cut {cut}: {len(keys)} elements cross the plane, {len(s) - before} added ({len(s)} in total)")


@command("CUTCLR", max_args=3)
def cmd_cutclr(c):
    """CUTCLR,<first>,[last],[step]: delete cuts (they can be rebuilt by adding elements again)."""
    a = c.int(1, required=True, what="first cut")
    store = cl.cut_store(c.interp)
    ids = select_ids(store, a, c.int(2, default=a), c.int(3, default=1))
    for k in ids:
        del store[k]
    if not ids:
        c.warn("no cut defined in that range")
    else:
        c.confirm(f"{len(ids)} cuts deleted: {cl._short(ids)}")


# ======================================================================================
# Submodels: CUT2SUB, CSECT, EXTRACTEXCAV, SPLITGROUP, TRANELEM, TRANVOL (spec 09 section 5)
# ======================================================================================
def _dest(c, k: int, what: str = "dest") -> int:
    n = _model_number(c, k, what)
    if n == c.interp.active_model:
        c.fail(f"<{what}> must differ from the active model (the source)")
    return n


def _cut_keys(c, cut: int) -> List[Tuple[int, int]]:
    try:
        keys = cl.get_cut(c.interp, cut)
    except SectionError as exc:
        c.fail(str(exc))
    if not keys:
        c.fail(f"cut {cut} is empty")
    return sorted(keys)


@command("CUT2SUB", max_args=3)
def cmd_cut2sub(c):
    """CUT2SUB,<cutnum>,<dest>,[solid]: submodel of the cut's elements in model <dest> (same group and element
    numbers); solid >= 1 turns shells into thick-shell solids for plotting."""
    cut = _cut_number(c)
    dest = _dest(c, 2)
    solid = c.int(3, default=-1)
    keys = _cut_keys(c, cut)
    sub = _run(c, cl.submodel, c.model, keys, solid_shells=solid >= 1)
    _store_model(c, dest, sub)
    c.confirm(f"CUT2SUB: model {dest} holds {sub.n_elements()} elements and {len(sub.nodes)} nodes of cut {cut}")


@command("CSECT", max_args=8)
def cmd_csect(c):
    """CSECT,<dest>,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>: cross-section model (unit thickness along n, same
    group and element numbers) of the cut's elements in model <dest>."""
    dest = _dest(c, 1)
    cut = _cut_number(c, 2)
    P = _vec(c, 3, "p")
    n = _vec(c, 6, "n")
    keys = _cut_keys(c, cut)
    have, miss = cl.model_elements(c.model, keys)
    if miss:
        c.warn(f"{len(miss)} elements of cut {cut} are not in the active model: {cl._short(cl._ge(miss))}")
    res = _run(c, cl.csect_model, c.model, have, P, n, source_model=c.interp.active_model, cut=cut)
    _store_model(c, dest, res.model)
    c.confirm(f"CSECT: model {dest} holds {res.elements} section elements ({len(res.pieces.pieces)} pieces, "
              f"{len(res.pieces.beams)} beams) and {res.nodes} nodes")


@command("EXTRACTEXCAV", max_args=1)
def cmd_extractexcav(c):
    """EXTRACTEXCAV,<Model>: submodel of the excavation elements (explicit ETYPE 2 SOLID/PLANE) with their
    nodes, interaction flags and soil layers; the active model is unchanged."""
    dest = _dest(c, 1, "Model")
    sub, n = _run(c, cl.extract_excavation, c.model)
    _store_model(c, dest, sub)
    nint = sum(1 for nd in sub.nodes.values() if 0 in nd.flags)
    c.confirm(f"EXTRACTEXCAV: model {dest} holds {n} excavation elements, {len(sub.nodes)} nodes "
              f"({nint} interaction nodes)")


@command("SPLITGROUP", max_args=3)
def cmd_splitgroup(c):
    """SPLITGROUP,<group>,<split>,[dir]: split a group by the plane of the first shell of SHELL group <split>,
    or (dir = X, Y, Z) by the plane <axis> = <split>; elements on the + side move to a new group."""
    from ...elements.base import ElementError
    from ...elements.shell import local_frame
    m = c.model
    gid = c.int(1, required=True, what="group")
    if gid not in m.groups:
        c.fail(f"group {gid} is not defined")
    if c.given(3):
        d = c.word(3)
        if d not in ("X", "Y", "Z"):
            c.fail("[dir] must be X, Y or Z")
        x = c.float(2, required=True, what="split")
        n = np.zeros(3)
        n["XYZ".index(d)] = 1.0
        P = n * x
        what = f"plane {d} = {fmt_num(x)}"
    else:
        sg = c.int(2, required=True, what="split")
        grp = m.groups.get(sg)
        if grp is None or grp.type != SHELL:
            c.fail(f"<split> group {sg} must be a SHELL group (or give [dir] X, Y or Z with a coordinate)")
        if not grp.elements:
            c.fail(f"SHELL group {sg} has no element")
        first = grp.sorted_elements()[0]
        corners = cl._shell_corners(first.nodes)
        v = cl.model_view(m)
        if len(corners) < 3 or not v.defined(corners):
            c.fail(f"the first shell of group {sg} has undefined nodes")
        try:
            Lam, P, _ = local_frame(np.array([v.P[k] for k in corners]))
        except ElementError as exc:
            c.fail(f"the first shell of group {sg}: {exc}")
        n = Lam[2]
        what = f"the plane of shell {sg}/{first.id}"
        if sg == gid:
            c.warn("the splitting shell group is the group being split")
    res = _run(c, cl.split_group, m, gid, P, n)
    if res.new_group is None:
        return
    k = _remap_cuts(c, gid, res.moved, res.new_group)
    if k:
        c.info(f"{k} cut entries follow the moved elements (cuts are session-global)")
    c.confirm(f"SPLITGROUP: {len(res.moved)} elements of group {gid} on the + side of {what} moved to group "
              f"{res.new_group} ({res.kept} stay)")


def _transfer(c, dest: int, keys) -> None:
    new = dest not in c.interp.models
    target = SSIModel() if new else c.interp.models[dest]
    res = _run(c, cl.transfer_elements, c.model, target, keys)     # validates before it changes anything
    if new:
        c.interp.models[dest] = target
        c.info(f"model {dest} created")
    c.confirm(f"{c.name}: {res.elements} elements and {res.nodes} nodes written to model {dest}"
              + (f", {res.tables_added} table entries added" if res.tables_added else ""))


@command("TRANELEM", max_args=5)
def cmd_tranelem(c):
    """TRANELEM,<dest>,<group>,<begin>,<end>,<stride>: copy elements of a group into model <dest>
    (overwrites elements, group type and nodes there)."""
    dest = _dest(c, 1)
    gid = c.int(2, required=True, what="group")
    grp = c.model.groups.get(gid)
    if grp is None:
        c.fail(f"group {gid} is not defined in the active model")
    a = c.int(3, required=True, what="begin")
    ids = select_ids(grp.elements, a, c.int(4, default=a), c.int(5, default=1))
    if not ids:
        c.fail(f"no element of group {gid} in the range")
    _transfer(c, dest, [(gid, e) for e in ids])


@command("TRANVOL", max_args=7)
def cmd_tranvol(c):
    """TRANVOL,<dest>,[Xmin],[Xmax],[Ymin],[Ymax],[Zmin],[Zmax]: copy the elements inside the box (CUTVOL
    test) into model <dest>."""
    dest = _dest(c, 1)
    m = c.model
    v = cl.model_view(m)
    try:
        b = cl.box_bounds(v, _box(c, 2))
    except SectionError as exc:
        c.fail(str(exc))
    keys = cl.elements_in_box(m, b, v)
    if not keys:
        c.fail("no element inside the box")
    _transfer(c, dest, keys)


# ======================================================================================
# Element stresses: READSTR, SECDATAOPT (spec 10 sections 3.9, 3.10)
# ======================================================================================
@command("READSTR", max_args=2)
def cmd_readstr(c):
    """READSTR,<Filename>,[Dir]: attach the element stresses of an .ess file to the active model."""
    name = c.str(1)
    if not name:
        c.fail("<Filename> required")
    if c.given(2):
        base = Path(c.str(2)).expanduser()
        if not base.is_absolute():
            base = c.interp.cwd / base
        path = base / name
        if not path.exists():
            c.fail(f"file {path} not found")
    else:
        path = c.input_path(name)
    notes = Notes()
    try:
        ess = cl.read_ess(path)
        keys = cl.check_ess_against_model(c.model, ess, notes)
    except SectionError as exc:
        _flush(c, notes)
        c.fail(str(exc))
    _flush(c, notes)
    maxima = bool(cl.MAXIMA_NAME.search(path.name))
    if maxima:
        c.warn("the file holds per-component absolute maxima, which do not occur at the same time: section "
               "resultants computed from them are not physically consistent (D-SEC-06)")
    n = cl.attach_stresses(c.model, ess, keys, maxima=maxima)
    c.confirm(f"READSTR: stresses of {n} elements loaded from {path.name}")


@command("SECDATAOPT", max_args=1)
def cmd_secdataopt(c):
    """SECDATAOPT,<flag>: STRESS writes NSTRESS/ESTRESS_<n>.ess frames for the whole history (1) or not (0)."""
    flag = c.int(1, required=True, what="flag")
    if flag not in (0, 1):
        c.fail("<flag> must be 0 (no .ess files) or 1 (save the .ess frames)")
    c.model.options.set_record(make_record("SECDATAOPT", [str(flag)]))
    c.confirm("STRESS " + ("writes NSTRESS/ESTRESS_<n>.ess frames" if flag else "writes no .ess frames")
              + " (takes effect at the next AFWRITE / RUNSTRESS)")


# ======================================================================================
# Calculations: CALCPAR, CALCMOI, CALCC, CALCM (spec 10 sections 3.3, 3.6-3.8)
# ======================================================================================
def _fmt(x: float) -> str:
    return f"{x: .6E}"


def _section(c, n: np.ndarray, r: np.ndarray):
    """Pieces, axes and properties of the active model's section (CALCPAR / CALCMOI).  The BEAMS that
    cross the plane and the elements without a section piece are listed in a warning (spec 10 section
    3.2), also on a CSECT model, whose section never includes the beams it draws."""
    notes = Notes()
    try:
        ex, ey, ez = cl.section_axes(n, r)
        sec = cl.active_section(c.model, ez, notes)
        if sec.skipped is not None:
            cl.report_skipped(sec.skipped, notes)
        props = cl.section_properties(sec.pieces, sec.point, ex, ey, ez)
    except SectionError as exc:
        _flush(c, notes)
        c.fail(str(exc))
    return sec, props, notes


@command("CALCPAR", max_args=8)
def cmd_calcpar(c):
    """CALCPAR,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sysno>,[verbose]: area, centroid, moments of inertia and the six
    section resultants (D-SEC-02 order) of the active cross-section model (on any other model: the base
    section of the whole structure, Q-11)."""
    n, r = _vec(c, 1, "normal"), _vec(c, 4, "right")
    sysno = c.int(7, default=0)
    verbose = c.int(8, default=-1)
    m = c.model
    sec, props, notes = _section(c, n, r)
    table = cl.stress_table(m)
    R = np.zeros(6)
    if not table:
        unmatched = cl.stress_info(m).get("unmatched")
        if unmatched:          # CSECT found stresses on the source model, none of the cut's elements
            notes.warn(f"{unmatched}: the resultants are reported as 0")
        else:
            notes.warn("no element stresses are loaded (READSTR on the original model before CSECT): the resultants "
                       "are reported as 0")
    else:
        S, ok = cl.gather_stresses(sec.pieces, table, sec.parents)
        if not ok.any():
            mismatch = cl.stress_mismatch(table, [p.key for p in sec.pieces])
            notes.warn((mismatch + ": the resultants are reported as 0") if mismatch else
                       "no section piece has stress data: the resultants are reported as 0")
        elif not ok.all():
            miss = [p.key for p, k in zip(sec.pieces, ok) if not k]
            notes.warn(f"{len(miss)} section elements without stress data are left out of the resultants: "
                       f"{cl._short(cl._ge(miss))}")
        T = cl.resultant_operators(sec.pieces, props, cl.bending_on(m))
        R = cl.resultants(T, S, ok)
        if cl.stress_info(m).get("maxima"):
            notes.warn("the stresses are absolute maxima (not simultaneous): the resultants are not physically "
                       "consistent (D-SEC-06)")
    cl.store_cut_system(m, sysno, props.C, props.ex, props.ey, props.ez, notes)
    _flush(c, notes)
    vals = props.values() + [float(x) for x in R]
    c.interp.session["calcpar"] = dict(zip(cl.CALCPAR_NAMES, vals))
    if verbose == -1:
        c.info(f"CALCPAR: section of {len(sec.pieces)} pieces (local axes ex, ey in the plane, ez = n; "
               f"forces of the +n side on the -n side, Fz > 0 tension)")
        for name, x in zip(cl.CALCPAR_NAMES, vals):
            c.info(f"{name:>4s} = {_fmt(x)}")
    else:
        c.info(" ".join(_fmt(x) for x in vals))


@command("CALCMOI", max_args=7)
def cmd_calcmoi(c):
    """CALCMOI,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sysno>: section moments of inertia of the active cross-section
    model (on any other model: the base section of the whole structure, Q-11) about its area centroid, in the
    local cut axes."""
    n, r = _vec(c, 1, "normal"), _vec(c, 4, "right")
    sysno = c.int(7, default=0)
    sec, props, notes = _section(c, n, r)
    cl.store_cut_system(c.model, sysno, props.C, props.ex, props.ey, props.ez, notes)
    _flush(c, notes)
    names = cl.CALCPAR_NAMES[:8]
    vals = props.values()
    c.interp.session["calcmoi"] = dict(zip(names, vals))
    for name, x in zip(names, vals):
        c.info(f"{name:>4s} = {_fmt(x)}")
    I1, I2, ang = props.principal()
    c.info(f"principal moments I1 = {_fmt(I1)}, I2 = {_fmt(I2)} (I1 axis at {ang:.3f} deg from local x about n)")


@command("CALCC", max_args=0)
def cmd_calcc(c):
    """CALCC: volume centroid of the active model (SOLID, SHELL, PLANE, BEAMS elements; on a cross-section
    model the BEAMS, which are not section pieces, are left out so that it equals the area centroid)."""
    vols = _run(c, cl.centroid_volumes, c.model)
    try:
        V, C = cl.volume_centroid(vols)
    except SectionError as exc:
        c.fail(str(exc))
    c.interp.session["calcc"] = {"Volume": V, "Xc": float(C[0]), "Yc": float(C[1]), "Zc": float(C[2])}
    c.info(f"CALCC: Xc = {_fmt(C[0])}, Yc = {_fmt(C[1])}, Zc = {_fmt(C[2])}; volume = {_fmt(V)} "
           f"({len(vols)} elements)")


@command("CALCM", max_args=0)
def cmd_calcm(c):
    """CALCM: element mass (material weight x volume / g) and lumped MT/MR masses of the active model."""
    s = _run(c, cl.model_masses, c.model)
    tot = s.total
    c.interp.session["calcm"] = {"element_weight": s.weight, "element_mass": s.mass, "soil_weight": s.soil_weight,
                                 "soil_mass": s.soil_mass, "lumped": [float(x) for x in s.lumped],
                                 "total": [float(x) for x in tot]}
    c.info(f"CALCM (g = {fmt_num(s.gravity)}): element weight = {_fmt(s.weight)}, element mass = {_fmt(s.mass)}")
    if s.soil_weight:
        c.info(f"excavated soil (not included): weight = {_fmt(s.soil_weight)}, mass = {_fmt(s.soil_mass)}")
    L = s.lumped
    c.info(f"lumped masses: Mx = {_fmt(L[0])}, My = {_fmt(L[1])}, Mz = {_fmt(L[2])}; rotational inertia: "
           f"Mxx = {_fmt(L[3])}, Myy = {_fmt(L[4])}, Mzz = {_fmt(L[5])}")
    c.info(f"total mass: X = {_fmt(tot[0])}, Y = {_fmt(tot[1])}, Z = {_fmt(tot[2])}")


# ======================================================================================
# CALCSECTHIST (spec 10 section 3.4)
# ======================================================================================
@command("CALCSECTHIST", max_args=14)
def cmd_calcsecthist(c):
    """CALCSECTHIST,<infile>,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sysno>,[ts],<outfile>: section
    resultants of a cut of the active (original) model for every .ess frame of a list file (7-column CSV with
    a final MAX row)."""
    if not c.given(1) or not c.given(14):
        c.fail("<infile> and <outfile> required")
    infile = c.input_path(c.str(1))
    cut = _cut_number(c, 2)
    P = _vec(c, 3, "p")
    n, r = _vec(c, 6, "normal"), _vec(c, 9, "right")
    sysno = c.int(12, default=0)
    ts = c.float(13, default=0.0)
    if ts < 0:
        ts = 0.0
    out = c.output_path(c.str(14))
    m = c.model
    keys = _cut_keys(c, cut)
    notes = Notes()
    if isinstance(m.ui_state.get("csect"), dict):
        notes.warn("the active model is a CSECT cross-section model: CALCSECTHIST works on the original model whose "
                   "elements the .ess frames describe")
    try:
        ex, ey, ez = cl.section_axes(n, r)
        frames = cl.read_frame_list(infile, c.interp.cwd)
        sp = cl.section_pieces(m, keys, P, ez, notes=notes)
        cl.report_skipped(sp, notes)
        props = cl.section_properties(sp.pieces, P, ex, ey, ez)
        T = cl.resultant_operators(sp.pieces, props, cl.bending_on(m))
        R = cl.history_rows(T, sp.pieces, frames, notes)
    except SectionError as exc:
        _flush(c, notes)
        c.fail(str(exc))
    cl.store_cut_system(m, sysno, props.C, ex, ey, ez, notes)
    out.parent.mkdir(parents=True, exist_ok=True)
    cl.write_section_csv(out, R, ts)
    _flush(c, notes)
    mx = cl.signed_absmax(R)
    c.interp.session["calcsecthist"] = {"file": str(out), "frames": len(frames), "max": [float(x) for x in mx]}
    c.info("CALCSECTHIST: MAX " + ", ".join(f"{k} = {_fmt(x)}" for k, x in zip(cl.RESULTANT_NAMES, mx)))
    c.confirm(f"CALCSECTHIST: {len(frames)} frames, {len(sp.pieces)} section pieces (A = {_fmt(props.area)}) "
              f"-> {out.name}")
