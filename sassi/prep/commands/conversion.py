"""File-conversion commands (requirements 3.4.K, manual section 9.11; spec 09 section 6; spec 04
sections 10-12; decisions D-ANS-01..08, D-MDL-03):

* ``CONVERT,ANSYS,<model>,<filename>,<gravity>,[<prefile>],[<damp>]`` -- ANSYS ``.cdb`` -> model
  (:mod:`sassi.io.ansys_cdb`); ``<prefile>`` (the dialog's "Output .pre File Name") and ``<damp>``
  (damping ratio of materials without DMPR) are SASSI-EDU extensions of the manual's argument list.
  The element map ``<model>_cdb.map`` (D-ANS-02) is written in the model directory (MDL), else as
  ``<cdb name>_cdb.map`` in the working directory.
* ``CONVERT,SSI,<model>,<filename>,[<prefile>]`` -- SASSI-EDU HOUSE deck ``.hou`` (+ ``.sit``,
  ``.poi`` with the same base name) -> model, the round-trip helper of D-ANS-08 (legacy SASSI2000
  fixed-format decks are out of scope).
* ``CONVERT,STRUDL,...`` -- not applicable in this version (manual).
* ``ANSYS,[FileName],[Dir],[<dmap>]`` -- APDL export (:mod:`sassi.io.apdl`); default
  ``<model>.inp`` in the model directory, ``unnamed.inp`` in the working directory without MDL
  (D-MDL-03); ``<dmap>`` 0 ``DMPR = 2 beta`` (D-ANS-06), 1 exact complex-modulus mapping (R1 1.2),
  2 no damping; ``EDUOPT,ANSYSMODERN,1`` selects current ANSYS elements.
* ``ANSYSREFORMAT,<Org>,<Map>`` -- copy model ``<Org>`` into the (empty) active model, one BEAMS
  group per end-release pattern (ANSYS sets releases per element type), map file
  ``old_group old_elem new_group new_elem``.
* ``ANSYSMODELTYPE,<type>`` (P2, Option AA) -- stored; ``GENMATRIXDAMP`` -- message only (manual:
  not usable in this version).

The destination model of CONVERT (blank = active) is *replaced* by the converted model; its MDL
name and path are kept so that WRITE/AFWRITE go where the user set them.  The active model does
not change.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ...io import decks
from ...io.ansys_cdb import convert_cdb, read_cdb
from ...io.apdl import DAMPING_MODES, ApdlOptions, export_apdl
from ...model import SSIModel, fmt_num, make_record
from ...model.entities import (BeamSection, Element, ElementRequest, Group, Material, MatrixProp, SoilLayer,
                               SpringProp)
from ..options import eduopt
from ..registry import command

SOLID, BEAMS, PLANE, GENERAL = 1, 2, 4, 9


# ======================================================================================
# helpers
# ======================================================================================
def _target(c, k: int) -> int:
    """Destination model number of argument ``k`` (blank = the active model)."""
    n = c.int(k, default=c.interp.active_model, what="model")
    if n < 0:
        c.fail("model numbers are >= 0")
    return n


def _install(c, n: int, new: SSIModel) -> None:
    """Replace model ``n`` by ``new``; the MDL name and path of the old model are kept."""
    old = c.interp.models.get(n)
    if old is not None:
        new.name, new.path = old.name, old.path
        if not old.is_empty():
            c.warn(f"model {n} replaced by the converted model")
    c.interp.models[n] = new


def _write_pre(c, model: SSIModel, name: str) -> None:
    from ..writer import write_pre
    target = c.output_path(name if Path(name).suffix else name + ".pre")
    text, notes = write_pre(model, filename=target.name)
    try:
        target.write_text(text, encoding="utf-8")
    except OSError as exc:
        c.fail(f"cannot write {target}: {exc.strerror or exc}")
    for t in notes:
        c.warn(t)
    c.info(f"converted model written to {target}")


# ======================================================================================
# CONVERT
# ======================================================================================
@command("CONVERT")
def cmd_convert(c):
    """CONVERT,<ConSel>,<model>,<filename>,...: run a file converter (SSI, ANSYS; STRUDL not applicable)."""
    if c.nargs == 0:
        c.info("CONVERT: converters (the Model > Converters dialogs in the GUI):\n"
               "  CONVERT,ANSYS,<model>,<file.cdb>,<gravity>,[<prefile>],[<damp>]   ANSYS CDWRITE file\n"
               "  CONVERT,SSI,<model>,<file.hou>,[<prefile>]                       SASSI-EDU HOUSE deck "
               "(+ .sit, .poi)\n"
               "  <model> blank = active model; gravity 32.2 (ft/s^2), 386.4 (in/s^2) or 9.81 (m/s^2). The "
               "converters have had limited testing: check every converted model before analysis.")
        return
    sel = c.word(1)
    if sel == "STRUDL":
        c.warn("CONVERT,STRUDL: the GT-STRUDL converter is not applicable in this version")
        return
    if sel not in ("ANSYS", "SSI"):
        c.fail(f"<ConSel> must be SSI, ANSYS or STRUDL (got '{c.str(1)}')")
    n = _target(c, 2)
    fname = c.str(3)
    if not fname:
        c.fail("<filename> is required")
    path = c.input_path(fname)
    if sel == "ANSYS":
        _convert_ansys(c, n, path)
    else:
        _convert_ssi(c, n, path)


def _convert_ansys(c, n: int, path: Path) -> None:
    g = c.float(4, what="gravity")
    if g is None:
        c.fail("<gravity> is required: 32.2 (ft/s^2), 386.4 (in/s^2) or 9.81 (m/s^2) -- it converts DENS to weight")
    if g <= 0:
        c.fail("<gravity> must be > 0")
    damp = c.float(6, what="damp")
    if damp is not None and not 0.0 <= damp < 0.5:
        c.fail("<damp> must be a damping ratio in [0, 0.5)")
    try:
        cdb = read_cdb(path)
    except OSError as exc:
        c.fail(f"cannot read {path}: {exc.strerror or exc}")
    if not cdb.nodes and not cdb.elements:
        c.fail(f"{path.name}: no nodes or elements found (not an ANSYS CDWRITE file?)")
    res = convert_cdb(cdb, g, damp)
    model = res.model
    _install(c, n, model)
    for kind, text in res.messages:
        (c.warn if kind == "warning" else c.info)(text)
    # element map (D-ANS-02)
    if model.name and model.path:
        mp = Path(model.path) / f"{model.name}_cdb.map"
    else:
        mp = c.interp.cwd / f"{path.stem}_cdb.map"          # no MDL: working directory (as ANSYS, D-MDL-03)
    try:
        mp.write_text(res.map_text(path.name), encoding="utf-8")
        c.info(f"element map (ansys_elem group elem) written to {mp}")
    except OSError as exc:
        c.warn(f"cannot write the element map {mp}: {exc.strerror or exc}")
    if c.given(5):
        _write_pre(c, model, c.str(5))
    c.info("CONVERT,ANSYS: only a subset of ANSYS elements and options is converted -- check the converted "
           "model before analysis (warnings above)")
    c.confirm(f"{path.name} converted to model {n}")


# ======================================================================================
# CONVERT,SSI: SASSI-EDU decks -> model (D-ANS-08)
# ======================================================================================
def _rec(name: str, values: Sequence) -> "object":
    return make_record(name, [v if isinstance(v, str) else fmt_num(v) for v in values])


def _bits(code: str) -> List[int]:
    s = str(code).strip() or "0"
    s = s.zfill(6)[-6:]
    return [1 if ch not in "0" else 0 for ch in s]


def model_from_decks(house: decks.Deck, site: Optional[decks.Deck] = None,
                     point: Optional[decks.Deck] = None) -> Tuple[SSIModel, List[str]]:
    """Rebuild a model from the HOUSE deck (and optionally the SITE and POINT decks) that AFWRITE
    wrote (D-ANS-08).  Returns ``(model, notes)``.

    The decks carry *resolved* data: global node coordinates (local systems are not restored),
    explicit ETYPE values (an implicit SOLID/PLANE 0 comes back as the 1/2 it resolved to; for the
    other types, where ETYPE 1 is the implicit value, 1 comes back as 0; D-HOU-01), and the
    AFWRITE write-time fix-ups (D-AFW-06) -- gap nodes (at the origin, all DOFs fixed, unreferenced)
    are dropped and the all-fixed flags AFWRITE gives to nodes without DOFs are cleared (AFWRITE
    sets them again).  Options without a deck parameter (e.g. ANALYS, MOTION) are not restored.
    """
    m = SSIModel()
    notes: List[str] = []
    h = house
    m.title = h["title"]
    # ---------------------------------------------------------------- groups and elements
    gtype: Dict[int, int] = {}
    for r in h.rows("groups"):
        g = Group(int(r["id"]), int(r["type"]), r["title"])
        m.groups[g.id] = g
        gtype[g.id] = g.type
    dof_nodes = set()
    refs = set()
    for r in h.rows("elements"):
        gid = int(r["group"])
        if gid not in m.groups:
            notes.append(f"element {r['id']} of undefined group {gid} skipped")
            continue
        nodes = [int(r[f"n{k}"]) for k in range(1, 9)]
        while nodes and nodes[-1] == 0:
            nodes.pop()
        etype = int(r["etype"])
        if etype == 1 and gtype[gid] not in (SOLID, PLANE):
            etype = 0             # resolved 1 = the implicit value for BEAMS/SHELL/SPRING/GENERAL (D-HOU-01)
        el = Element(int(r["id"]), nodes, mat=int(r["mat"]), prop=int(r["prop"]), etype=etype,
                     eint=int(r["eint"]), thick=float(r["thick"]), ki=_bits(r["ki"]), kj=_bits(r["kj"]))
        m.groups[gid].elements[el.id] = el
        refs.update(x for x in nodes if x)
        dof_nodes.update(x for x in (nodes[:2] if gtype[gid] in (BEAMS, GENERAL) else nodes) if x)
    if m.groups:
        m.group_active = max(m.groups)
    inter = {int(r["id"]) for r in h.rows("interaction")}
    mass_nodes = {int(r["node"]) for r in h.rows("masses")}
    # ---------------------------------------------------------------- nodes
    gaps = cleared = 0
    for r in h.rows("nodes"):
        nid = int(r["id"])
        xyz = (float(r["x"]), float(r["y"]), float(r["z"]))
        fix = [int(r[k]) for k in ("fx", "fy", "fz", "frx", "fry", "frz")]
        if (xyz == (0.0, 0.0, 0.0) and all(fix) and nid not in refs and nid not in inter
                and nid not in mass_nodes):
            gaps += 1                                  # AFWRITE gap node (Warning 1)
            continue
        m.define_node(nid, xyz, 0)
        if nid not in dof_nodes and nid not in inter and nid not in mass_nodes and all(fix):
            fix = [0] * 6                              # AFWRITE fix-up of nodes without DOFs (Warning 4)
            cleared += 1
        m.nodes[nid].fix = fix
        if nid in inter:
            m.nodes[nid].flags.add(0)
    if gaps:
        notes.append(f"{gaps} AFWRITE gap nodes (origin, all DOFs fixed, unreferenced) not restored")
    if cleared:
        notes.append(f"all-DOF fixities of {cleared} nodes without DOFs not restored (AFWRITE sets them, Warning 4)")
    missing_int = sorted(inter - set(m.nodes))
    if missing_int:
        notes.append(f"interaction nodes {missing_int[:10]} are not in the node table")
    # ---------------------------------------------------------------- property tables
    for r in h.rows("materials"):
        m.materials[int(r["id"])] = Material(int(r["id"]), r["val1"], r["val2"], r["weight"], r["pdamp"],
                                             r["sdamp"], int(r["type"]))
    for r in h.rows("layers"):
        m.layers[int(r["no"])] = SoilLayer(int(r["no"]), r["thick"], r["weight"], r["vp"], r["vs"], r["dp"], r["ds"])
    for r in h.rows("beamprops"):
        m.sections[int(r["id"])] = BeamSection(int(r["id"]), r["axial"], r["shear2"], r["shear3"], r["tors"],
                                               r["flex2"], r["flex3"])
    for r in h.rows("springprops"):
        m.springs[int(r["id"])] = SpringProp(int(r["id"]), r["scx"], r["scy"], r["scz"], r["scxx"], r["scyy"],
                                             r["sczz"], r["damp"])
    for r in h.rows("matrices"):
        pid, row = int(r["prop"]), int(r["row"])
        p = m.matrices.setdefault(pid, MatrixProp(pid))
        p.set_row(str(r["kind"]).upper(), row, [float(r[f"t{k}"]) for k in range(1, 14 - row)])
    for r in h.rows("masses"):
        n = int(r["node"])
        t = [float(r["mx"]), float(r["my"]), float(r["mz"])]
        q = [float(r["mxx"]), float(r["myy"]), float(r["mzz"])]
        if any(t) or not any(q):
            m.tmass[n] = t
            m.tmass_history.touch(n)
        if any(q):
            m.rmass[n] = q
            m.rmass_history.touch(n)
        if int(r["units"]) != SSIModel.MUNITS_DEFAULT:
            m.mass_units[n] = int(r["units"])
    # ---------------------------------------------------------------- options
    opt = m.options
    opt.set_record(_rec("HOUSE", [h["gravity"], h["gelev"], h["opmode"], h["dim"], h["imp"], h["coh"],
                                  h["wpass"], h["me"], h["cmplxspec"]]))
    hx = [h["optimize"], h["supmode"], h["nsim"], h["nlssi"], h["ansys"]]
    if hx != [0, 0, 1, 0, 0]:
        opt.set_record(_rec("HOUSEX", hx))
    if int(h["incomp"]) != 1 or int(h["gmunits"]) != 0:
        opt.set_record(_rec("MOPT", [h["incomp"]] + ([h["gmunits"]] if int(h["gmunits"]) else [])))
    if int(h["cmodform"]) != 0:
        opt.set_record(_rec("CMODFORM", [h["cmodform"]]))
    if int(h["coh"]):
        opt.set_record(_rec("INCOH", [h[k] for k in ("gammax", "gammay", "gammaz", "alpha", "ngp", "ipr", "nmodes",
                                                     "met", "hseed", "vseed", "randphz")]))
    if int(h["wpass"]):
        opt.set_record(_rec("WPASS", [h["appv"], h["wang"], h["cohf"]]))
    for r in h.rows("me"):
        opt.set_entry("ME", int(r["no"]), _rec("ME", [r["no"], r["nfirst"], r["nlast"]]))
    for r in h.rows("symm"):
        opt.set_entry("SYMM", int(r["no"]), _rec("SYMM", [r["no"], r["type"], r["n1"], r["n2"], r["n3"]]))
    amp: Dict[int, List[float]] = {}
    for r in h.rows("amp"):
        vals = amp.setdefault(int(r["no"]), [])
        vals += [float(r["re"]), float(r["im"])] if int(h["cmplxspec"]) == 1 else [float(r["re"])]
    m.amp = amp
    # site layers from the HOUSE deck (TOPL and the half-space row)
    sl = h.rows("sitelayers")
    for r in sl:
        if int(r["no"]) not in m.layers:
            m.layers[int(r["no"])] = SoilLayer(int(r["no"]), r["thick"], r["weight"], r["vp"], r["vs"], r["dp"],
                                               r["ds"])
    topl = [int(r["no"]) for r in sl[:-1]]
    hs_house = int(sl[-1]["no"]) if sl else 0
    freq_set = 1
    if site is not None:
        s = site
        freq_set = int(s["freq"])
        opt.set_record(_rec("SITE", [s["opmode"], s["mode1"], s["fstep"], s["nl"], s["hs"], s["mode2"], s["wopt"],
                                     s["freq1"], s["freq2"], s["cl"], s["cm"], s["delt"], s["nft"], s["freq"]]))
        if int(s["soilmode"]):
            opt.set_record(_rec("SITEX", [s["soilmode"]]))
        if str(s["hslaw"]).upper() != "UNIFORM":
            opt.set_entry("EDUOPT", "HSLAW", _rec("EDUOPT", ["HSLAW", str(s["hslaw"]).upper()]))
        for tab in ("layers", "halfspace"):
            for r in s.rows(tab):
                if int(r["no"]) not in m.layers:
                    m.layers[int(r["no"])] = SoilLayer(int(r["no"]), r["thick"], r["weight"], r["vp"], r["vs"],
                                                       r["dp"], r["ds"])
        topl = [int(r["no"]) for r in s.rows("layers")]
        for r in s.rows("waves"):
            opt.set_entry("WAVE", int(r["type"]), _rec("WAVE", [r["type"], r["opt"], r["ratio1"], r["ratio2"],
                                                                r["angle"]]))
        fn = [int(r["number"]) for r in s.rows("freqs")]
    else:
        fn = [int(r["number"]) for r in h.rows("freqs")]
        rec = make_record("SITE")
        rec.set_arg(12, fmt_num(h["delt"]))
        rec.set_arg(13, fmt_num(h["nft"]))
        if h["delt"] > 0 and h["nft"] > 0 and abs(h["df"] - 1.0 / (h["delt"] * h["nft"])) > 1e-12 * max(h["df"], 1e-300):
            rec.set_arg(3, fmt_num(h["df"]))
        if hs_house:
            rec.set_arg(5, str(hs_house))
        opt.set_record(rec)
        notes.append("no .sit deck: SITE options restored only from the HOUSE deck (delt, nft, df, half-space layer)")
    m.topl = topl
    if fn:
        m.freq_sets[freq_set] = sorted(set(fn))
    if point is not None:
        p = point
        opt.set_record(_rec("POINT", [p["opmode"], p["layer"], p["rad"]]))
    return m, notes


def _read_deck(path: Path, module: str) -> decks.Deck:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return decks.read(path, module)


def _convert_ssi(c, n: int, path: Path) -> None:
    try:
        house = _read_deck(path, "HOUSE")
    except (ValueError, KeyError) as exc:
        c.fail(f"{path.name} is not a SASSI-EDU HOUSE deck ({exc}); legacy SASSI2000 fixed-format decks are "
               "not supported (D-ANS-08)")
    except OSError as exc:
        c.fail(f"cannot read {path}: {exc.strerror or exc}")
    site = point = None
    for ext, mod in ((".sit", "SITE"), (".poi", "POINT")):
        p = path.with_suffix(ext)
        if p.exists():
            try:
                d = _read_deck(p, mod)
            except (ValueError, KeyError, OSError) as exc:
                c.warn(f"{p.name} not read: {exc}")
                continue
            if mod == "SITE":
                site = d
            else:
                point = d
            c.info(f"{p.name} read ({mod} options)")
    model, notes = model_from_decks(house, site, point)
    _install(c, n, model)
    for t in notes:
        c.info(t)
    if c.given(4):
        _write_pre(c, model, c.str(4))
    c.confirm(f"{path.name} converted to model {n}: {len(model.nodes)} nodes, {model.n_elements()} elements")


# ======================================================================================
# ANSYS (export)
# ======================================================================================
@command("ANSYS")
def cmd_ansys(c):
    """ANSYS,[FileName],[Dir],[<dmap>]: write the active model as ANSYS APDL input (<model>.inp)."""
    m = c.model
    fname = c.str(1)
    if not fname:
        if m.name:
            fname = m.name + ".inp"
        else:
            fname = "unnamed.inp"
            c.warn("the model has no name (MDL): unnamed.inp written in the working directory (D-MDL-03)")
    if not Path(fname).suffix:
        fname += ".inp"
    if c.given(2):
        d = Path(c.str(2)).expanduser()
        target = d / fname if d.is_absolute() else c.interp.cwd / d / fname
    elif m.name and m.path:
        target = c.output_path(fname)
    else:
        target = c.interp.cwd / fname if not Path(fname).is_absolute() else Path(fname)
    dmap = c.int(3, default=0, what="dmap")
    if dmap not in DAMPING_MODES:
        c.fail("<dmap> must be 0 (DMPR = 2 beta), 1 (exact complex-modulus mapping, R1 1.2) or 2 (no damping)")
    modern = str(eduopt(m, "ANSYSMODERN")).strip() == "1"
    try:
        res = export_apdl(m, ApdlOptions(modern=modern, damping=DAMPING_MODES[dmap]))
    except ValueError as exc:
        c.fail(str(exc))
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(res.text, encoding="utf-8")
    except OSError as exc:
        c.fail(f"cannot write {target}: {exc.strerror or exc}")
    for t in res.warnings:
        c.warn(t)
    for t in res.notes:
        c.info(t)
    if not m.groups:
        c.warn("the model has no elements")
    c.info(f"ANSYS: {res.counts.get('nodes', 0)} nodes, {res.counts.get('elements', 0)} elements; "
           + ("modern" if modern else "legacy ANSYS V11-15") + " elements; read it in ANSYS with /INPUT")
    c.confirm(f"APDL written to {target}")


# ======================================================================================
# ANSYSREFORMAT
# ======================================================================================
def reformat_beams(model: SSIModel) -> List[Tuple[int, int, int, int]]:
    """Regroup BEAMS groups by end-release pattern (KI, KJ) in place.

    The pattern of a group's first element keeps the group number; every other pattern gets a new
    group (numbered after the largest group) with the same title, elements renumbered 1..k in their
    old order.  EOUT requests of the moved elements follow them.  Returns the map rows
    ``(old_group, old_elem, new_group, new_elem)`` of every beam element.
    """
    rows: List[Tuple[int, int, int, int]] = []
    nxt = (max(model.groups) if model.groups else 0) + 1
    moved: Dict[Tuple[int, int], Tuple[int, int]] = {}
    for gid in sorted(model.groups):
        g = model.groups[gid]
        if g.type != BEAMS or not g.elements:
            continue
        parts: Dict[Tuple, List[Element]] = {}
        for e in g.sorted_elements():
            key = (tuple(int(bool(v)) for v in e.ki), tuple(int(bool(v)) for v in e.kj))
            parts.setdefault(key, []).append(e)
        if len(parts) == 1:
            rows += [(gid, e.id, gid, e.id) for e in g.sorted_elements()]
            continue
        keys = list(parts)
        g.elements = {}
        for i, key in enumerate(keys):
            new_gid = gid if i == 0 else nxt
            if i > 0:
                nxt += 1
                model.groups[new_gid] = Group(new_gid, BEAMS, g.title)
            ng = model.groups[new_gid]
            for k, e in enumerate(parts[key], start=1):
                rows.append((gid, e.id, new_gid, k))
                moved[(gid, e.id)] = (new_gid, k)
                ng.elements[k] = e.copy(new_id=k)
    rows.sort()
    if moved:
        eout: List[ElementRequest] = []
        for r in model.eout:
            split: Dict[int, List[int]] = {}
            for e in r.elements:
                ng, ne = moved.get((r.group, e), (r.group, e))
                split.setdefault(ng, []).append(ne)
            eout += [ElementRequest(list(r.codes), ng, sorted(v)) for ng, v in sorted(split.items())] or [r]
        model.eout = eout
    return rows


@command("ANSYSREFORMAT")
def cmd_ansysreformat(c):
    """ANSYSREFORMAT,<Org>,<Map>: copy model <Org> into the empty active model with one BEAMS group per
    release pattern (ANSYS KEYOPT(7)/(8) are per element type) and write the group map."""
    org = c.int(1, required=True, what="Org")
    mapname = c.str(2)
    if org not in c.interp.models:
        c.fail(f"model {org} is not in memory")
    if org == c.interp.active_model:
        c.fail("run ANSYSREFORMAT from an empty active model (ACTM to a new model number first)")
    if not c.model.is_empty():
        c.fail(f"the active model {c.interp.active_model} is not empty")
    new = c.interp.models[org].copy()
    rows = reformat_beams(new)
    c.interp.models[c.interp.active_model] = new
    n_new = len({(r[2]) for r in rows if r[0] != r[2]})
    if new.name:
        c.warn(f"model {c.interp.active_model} has the same name and path as model {org} (use MDL before AFWRITE)")
    if mapname:
        target = c.output_path(mapname)
        lines = ["! ANSYSREFORMAT map: old_group old_elem new_group new_elem"] + \
                [f"{a:8d} {b:8d} {g:8d} {e:8d}" for a, b, g, e in rows]
        try:
            target.write_text("\n".join(lines) + "\n", encoding="utf-8")
        except OSError as exc:
            c.fail(f"cannot write {target}: {exc.strerror or exc}")
        c.info(f"beam group map written to {target}")
    else:
        c.warn("<Map> not given: no map file written")
    c.confirm(f"model {org} reformatted into model {c.interp.active_model}: {n_new} new BEAMS groups")


# ======================================================================================
# ANSYSMODELTYPE, GENMATRIXDAMP
# ======================================================================================
@command("ANSYSMODELTYPE", max_args=1)
def cmd_ansysmodeltype(c):
    """ANSYSMODELTYPE,<type>: Option AA operation mode 1 embedded, 2 surface (P2: stored only)."""
    t = c.int(1, required=True, what="type")
    if t not in (1, 2):
        c.fail("<type> must be 1 (embedded) or 2 (surface)")
    c.model.options.set_record(make_record("ANSYSMODELTYPE", [str(t)]))
    c.warn(f"ANSYSMODELTYPE {'embedded' if t == 1 else 'surface'} stored; Option AA (ANSYS model input with "
           "HOUSEFSA/ANALYSFA) is not available in this build (tier P2)")


@command("GENMATRIXDAMP")
def cmd_genmatrixdamp(c):
    """GENMATRIXDAMP,<begin>,[end],[stride],<damp>: not usable in this version (manual)."""
    b = c.int(1, what="begin")
    e = c.int(2, default=-1, what="end")
    s = c.int(3, default=1, what="stride")
    d = c.float(4, what="damp")
    if b is not None and (e is None or e < 0):
        e = b
    what = (f" (matrix properties {b}..{e} step {max(s or 1, 1)}, damping {fmt_num(d)}: MXI = 2 damp MXR "
            "would be generated)") if b is not None and d is not None else ""
    c.warn("GENMATRIXDAMP is not usable in this version" + what + "; MXI unchanged -- enter MXI rows directly")
