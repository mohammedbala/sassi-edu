"""WRITE: serialise a model as ``.pre`` command text that INP reads back to the same state.

Requirements: UT-03 (WRITE -> INP into an empty model gives a deep-equal state), D-PAR-17
(canonical order of spec 07 section 9.2.47, X-commands only when non-default, byte-stable
output), D-PAR-14 (list and request clears are written first), spec 04 section 4.3 (header banner
and section comments).

Order of the output::

    banner
    * Model identification       TIT, [MDL], MOPT
    * Coordinate systems and nodes   LOC, CSYS, N, LOCAL (in dependency order; temporary
                                 defining nodes when needed), CSYS restore
    * Boundary conditions        D
    * Interaction nodes          INT
    * Material table             M
    * Soil layer table           L
    * Real property table        R, SC, MXR/MXI/MXM
    * Groups and elements        MACT/RACT, GROUP, GTIT, E, ETYPE/EINT/THICK/KI/KJ, [MTYPE],
                                 active group
    * Masses                     MT, MR, MUNITS
    * Loads                      F, MM, SYMM
    * Frequencies                FREQ
    * Analysis options           DYNP ... AOPT (OPTION_ORDER), extension sections
    [AFWR]

Exactness: numbers are written with :func:`sassi.model.values.fmt_num` (shortest round-trip
text), so every stored float reads back bit-identical.  A LOCAL system is re-created by LOCAL
itself, from its three nodes when they still have the coordinates they had at definition time,
otherwise from temporary nodes placed at the stored definition points (restored or deleted right
after); either way the computation is repeated bit for bit.  Data kept by MTYPE that the current
group type does not use is written under a type that accepts it, followed by MTYPE.  A model
with no active group ends the group section with a temporary group that is deleted again.

Extension points for other packages:

* :func:`register_writer` replaces the emission of one option command (e.g. a typed record
  writer that omits default values);
* :func:`register_section` adds a section (e.g. data kept in ``model.extensions``) before the
  AFWRITE section.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .. import PRODUCT, __version__
from ..model import SSIModel, fmt_num
from ..model.entities import DOF_NAMES, MATRIX_KINDS, MAX_ELEMENT_NODES, type_name
from ..model.geometry import euler_from_matrix
from .lexer import MAX_LINE, compress_ids, format_token, join_command

LineWriter = Callable[[SSIModel], List[str]]

#: canonical order of the analysis-option commands (spec 07 section 9.2.47 items 5-10)
OPTION_ORDER: List[Tuple[str, List[str]]] = [
    ("Free-field soil (SOIL) profile and output", ["DYNP", "SPRO", "SACC", "SRS", "SSTR", "SSAF", "SFOU", "DAMP"]),
    ("EQUAKE options", ["EQTIT", "EQUAKE", "RSIN", "RSOUT", "ACCIN", "ACCOUT", "TPSD", "CORR"]),
    ("Control motion", ["THFILE", "THTIT"]),
    ("SOIL options", ["SOIL", "SOILX"]),
    ("SITE options", ["SITE", "SITEX", "TOPL", "WAVE"]),
    ("POINT options", ["POINT"]),
    ("HOUSE options", ["HOUSE", "HOUSEX", "INCOH", "WPASS", "ME", "AMP"]),
    ("FORCE options", ["FORCE"]),
    ("ANALYS options", ["ANALYS", "ANALYSX"]),
    ("MOTION options", ["MOTION", "MOTIONX", "NOUT"]),
    ("STRESS options", ["STRESS", "STRESSX", "EOUT", "SECDATAOPT", "THSHLSTR"]),
    ("RELDISP options", ["RELD", "RELDX", "RELFILE", "RDND"]),
    ("SASSI-EDU options", ["CMODFORM", "EDUOPT", "BINOUT", "ANSYSMODELTYPE"]),
    ("NONLINEAR options (Option NON)", ["EQL", "P", "S", "B", "BBC", "BBCI", "BBCP", "BBCX", "BBCY",
                                        "NLSOIL", "NLSLAYER"]),
]
AFWRITE_SECTION = ("AFWRITE module selection", ["AOPT"])

#: extension X-commands: written only when non-default (requirements section 3.4.R) -- applies to
#: typed records, whose defaults are known
XCOMMANDS = {"SITEX", "SOILX", "HOUSEX", "ANALYSX", "MOTIONX", "STRESSX", "RELDX"}

WRITERS: Dict[str, LineWriter] = {}
SECTIONS: List[Tuple[str, LineWriter]] = []
#: commands written by the dedicated sections above (never repeated as "other options")
_HANDLED = {"MOPT", "SYMM"}


def register_writer(command: str, func: LineWriter) -> None:
    """Use ``func(model) -> lines`` to write option command ``command`` (replaces the generic writer)."""
    WRITERS[command.upper()] = func


def register_section(title: str, func: LineWriter) -> None:
    """Add a section written before the AFWRITE section (``func(model) -> lines``)."""
    for i, (t, _) in enumerate(SECTIONS):
        if t == title:
            SECTIONS[i] = (title, func)
            return
    SECTIONS.append((title, func))


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def _n(x) -> str:
    return fmt_num(x)


def _cmd(name: str, *vals) -> str:
    return ",".join([name] + [v if isinstance(v, str) else _n(v) for v in vals])


def _spec(name: str):
    from .registry import CATALOGUE, REGISTRY
    return REGISTRY.get(name) or CATALOGUE.get(name)


def record_line(rec) -> str:
    """One option record as a command line (generic writer)."""
    toks = rec.to_tokens()
    spec = _spec(rec.command)
    if getattr(rec, "legacy", "") and toks:
        parts = [format_token(t) for t in toks[:-1]] + [toks[-1]]
        return ",".join([rec.command] + parts)
    text_last = bool(spec and spec.text_from and len(toks) == spec.text_from)
    return join_command(rec.command, toks, text_last=text_last)


def _runs(items: Sequence[int], value_of) -> List[Tuple[int, int, object]]:
    """Runs of consecutive entries of ``items`` (sorted ids) with equal ``value_of(id)``."""
    out: List[Tuple[int, int, object]] = []
    for i in items:
        v = value_of(i)
        if out and out[-1][2] == v:
            out[-1] = (out[-1][0], i, v)
        else:
            out.append((i, i, v))
    return out


def _chunks(vals: Sequence, size: int = 10, avoid_leading_zero: bool = False) -> List[List]:
    out: List[List] = []
    i = 0
    vals = list(vals)
    while i < len(vals):
        j = min(i + size, len(vals))
        if avoid_leading_zero:
            while j < len(vals) and vals[j] == 0:
                j += 1
        out.append(vals[i:j])
        i = j
    return out


def _d_labels(fix: Sequence[int]) -> List[str]:
    f = list(fix)
    if all(f):
        return ["ALL"]
    if f[:3] == [1, 1, 1] and not any(f[3:]):
        return ["DISP"]
    if f[3:] == [1, 1, 1] and not any(f[:3]):
        return ["ROT"]
    return [DOF_NAMES[i] for i in range(6) if f[i]]


# --------------------------------------------------------------------------------------
# sections
# --------------------------------------------------------------------------------------
def _identification(m: SSIModel, mdl: bool) -> List[str]:
    out: List[str] = []
    if m.title:
        out.append(f"TIT,{m.title}")
    if mdl and m.name:
        out.append(join_command("MDL", [m.name, m.path]))
    rec = m.options.record("MOPT")
    if rec is not None:
        out.append(record_line(rec))
    return out


def _nodes(m: SSIModel) -> Tuple[List[str], List[str]]:
    """Coordinate systems and nodes in dependency order; returns (lines, notes).

    1. LOC systems (they depend on nothing) and the nodes stored in them and in system 0;
    2. LOCAL systems, each followed by the nodes stored in it.  A LOCAL system is written as
       ``LOCAL,ns,type,n1,n2,n3`` from its defining nodes when they are already written and still
       have the global coordinates they had at definition time.  Otherwise (helper nodes deleted,
       moved or redefined, a very common pattern) the three defining points are written as
       temporary nodes in system 0, LOCAL is applied to them, and the nodes are then restored to
       their current data (or deleted with NDEL when they no longer exist).  LOCAL repeats the
       original computation on bit-identical global points, so the system -- and every node
       stored in it -- reads back bit for bit (UT-03).  This block comes before the D / INT /
       load sections, so the NDEL cascade (D-MDL-11) removes nothing that is written later.
    3. Only a LOCAL record without its defining nodes/points (not produced by the commands) is
       written as the equivalent LOC (Euler angles, exact to rounding only), with a note.
    """
    out: List[str] = []
    notes: List[str] = []
    cur = [0]

    def csys(s: int) -> None:
        if cur[0] != s:
            out.append(f"CSYS,{s}")
            cur[0] = s

    def emit_node(i: int) -> None:
        n = m.nodes[i]
        csys(n.csys)
        out.append(_cmd("N", i, n.x, n.y, n.z))

    def emit_nodes(s: int) -> None:
        for i in m.nodes_in_system(s):
            emit_node(i)

    def from_nodes(sid: int) -> bool:
        """True when LOCAL ``sid`` can be re-created from its (already written) defining nodes."""
        cs = m.csys[sid]
        if len(cs.nodes) != 3 or len(cs.points) != 3:
            return False
        if not all(n in m.nodes and m.nodes[n].csys in done for n in cs.nodes):
            return False
        try:
            pts = [m.node_global(n) for n in cs.nodes]
        except KeyError:
            return False
        return all(np.array_equal(np.asarray(p), np.asarray(q, float)) for p, q in zip(pts, cs.points))

    def unchanged(n: int, p) -> bool:
        """Defining node ``n`` is written (its system is done) and is still at point ``p``."""
        return (n in m.nodes and m.nodes[n].csys in done
                and np.array_equal(m.node_global(n), np.asarray(p, float)))

    def via_helper_nodes(sid: int) -> None:
        cs = m.csys[sid]
        moved = [(n, p) for n, p in dict(zip(cs.nodes, cs.points)).items() if not unchanged(n, p)]
        what = "node" if len(moved) == 1 else "nodes"
        out.append(f"* LOCAL system {sid}: defining {what} {', '.join(str(n) for n, _ in moved)} placed at "
                   f"the definition points, then restored")
        csys(0)
        for n, p in moved:
            out.append(_cmd("N", n, *[float(v) for v in p]))
        out.append(_cmd("LOCAL", sid, cs.type, *cs.nodes))
        for n, _ in moved:
            if n not in m.nodes:
                out.append(f"NDEL,{n}")
            elif m.nodes[n].csys in done:
                emit_node(n)          # its system exists: restore the current data now
            # else: the node is written with the nodes of its own (later) system

    done = {0}
    # 1. LOC systems depend on nothing
    for sid in sorted(m.csys):
        cs = m.csys[sid]
        if cs.kind == "LOC":
            p = list(cs.params) + [0.0] * (6 - len(cs.params))
            out.append(_cmd("LOC", sid, cs.type, *p[:6]))
            done.add(sid)
    for s in sorted(done):
        emit_nodes(s)
    # 2. LOCAL systems: from their nodes when possible, else through temporary helper nodes
    pending = [sid for sid in sorted(m.csys) if sid not in done]
    while pending:
        sid = next((s for s in pending if from_nodes(s)), None)
        if sid is not None:
            cs = m.csys[sid]
            out.append(_cmd("LOCAL", sid, cs.type, *cs.nodes))
        else:
            sid = pending[0]
            cs = m.csys[sid]
            if len(cs.nodes) == 3 and len(cs.points) == 3:
                via_helper_nodes(sid)
            else:
                # 3. no defining data: equivalent LOC
                a, b, c = euler_from_matrix(cs.R)
                out.append(f"* LOCAL system {sid} written as the equivalent LOC (no defining nodes stored)")
                out.append(_cmd("LOC", sid, cs.type, *[float(v) for v in cs.origin], a, b, c))
                notes.append(f"coordinate system {sid} (LOCAL) written as LOC; the round trip is exact to "
                             f"rounding only")
        pending.remove(sid)
        done.add(sid)
        emit_nodes(sid)
    csys(m.csys_active)
    return out, notes


def _fixities(m: SSIModel) -> List[str]:
    out = []
    for a, b, fix in _runs(m.node_ids(), lambda i: tuple(m.nodes[i].fix)):
        if any(fix):
            out.append(",".join(["D", str(a), str(b), "1", "1"] + _d_labels(fix)))
    return out


def _interaction(m: SSIModel) -> List[str]:
    out = []
    ids = m.node_ids()
    for code in (0, 1, 2, 3):
        for a, b, on in _runs(ids, lambda i: code in m.nodes[i].flags):
            if on:
                out.append(f"INT,{a},{b},1,1,{code}")
    return out


def _tables(m: SSIModel) -> Tuple[List[str], List[str], List[str]]:
    mat = [_cmd("M", k, v.val1, v.val2, v.weight, v.pdamp, v.sdamp, v.mtype) for k, v in sorted(m.materials.items())]
    lay = [_cmd("L", k, v.thick, v.weight, v.vp, v.vs, v.pdamp, v.sdamp) for k, v in sorted(m.layers.items())]
    real = [_cmd("R", k, v.axial, v.shear2, v.shear3, v.tors, v.flex2, v.flex3) for k, v in sorted(m.sections.items())]
    real += [_cmd("SC", k, v.scx, v.scy, v.scz, v.scxx, v.scyy, v.sczz, v.damp) for k, v in sorted(m.springs.items())]
    for k, p in sorted(m.matrices.items()):
        for kind in MATRIX_KINDS:
            for r, terms in sorted(p.rows.get(kind, {}).items()):
                t = list(terms)
                while len(t) > 1 and t[-1] == 0.0:
                    t.pop()
                real.append(_cmd("MX" + kind, k, r, *t))
    return mat, lay, real


# group type codes (sassi.conventions.ELEMENT_TYPE_NAMES)
_SOLID, _BEAMS, _SHELL, _PLANE, _TSHELL, _SPRING, _GENERAL = 1, 2, 3, 4, 5, 7, 9


def _accepts(cmd: str, t: int, values) -> bool:
    """Would INP accept attribute command ``cmd`` with these (non-default) values in a group of
    type ``t``?  Mirrors the checks of :mod:`sassi.prep.commands.elements`."""
    if cmd == "EINT":
        return (t == _SOLID and all(v in (0, 1, 2) for v in values)) or \
               (t == _TSHELL and all(v in (0, 1) for v in values))
    if cmd == "THICK":
        return t in (_SHELL, _TSHELL)
    if cmd in ("KI", "KJ"):
        return t == _BEAMS
    return True                      # ETYPE is stored for every type


#: type used to write an attribute the group's current type does not accept (first that accepts)
_CARRIERS = {"EINT": (_TSHELL, _SOLID), "THICK": (_SHELL, _TSHELL), "KI": (_BEAMS,), "KJ": (_BEAMS,)}


def _node_carrier(t: int, maxn: int) -> int:
    """Type in which E lines with ``maxn`` nodes are accepted (``t`` itself when possible)."""
    if maxn <= MAX_ELEMENT_NODES[t]:
        return t
    if maxn <= MAX_ELEMENT_NODES[_BEAMS]:
        return _BEAMS
    if maxn <= MAX_ELEMENT_NODES[_TSHELL]:
        return _TSHELL
    return _SOLID


def _group_lines(g, state: Dict[str, Optional[int]]) -> List[str]:
    """Lines of one group.

    MTYPE keeps the data the new type does not use (spec 08 section 5.24): elements with more
    nodes than the type allows, EINT of a former SOLID group, THICK of a former SHELL group, beam
    releases.  INP accepts such data only in a type that uses it, so the group is created with a
    type that accepts its E lines, switched with MTYPE for each attribute the current type does
    not accept, and finally set back to its current type with MTYPE.  A model without MTYPE
    history is written with one GROUP line and no MTYPE.
    """
    out: List[str] = []
    T = g.type
    els = g.sorted_elements()
    ids = [e.id for e in els]
    cur = _node_carrier(T, max((len(e.nodes) for e in els), default=0))
    out.append(f"GROUP,{g.id},{type_name(cur)}")
    if g.title:
        out.append(f"GTIT,{g.id},{g.title}")
    for e in els:
        if e.mat != state["mact"]:
            out.append(f"MACT,{e.mat}")
            state["mact"] = e.mat
        if e.prop != state["ract"]:
            out.append(f"RACT,{e.prop}")
            state["ract"] = e.prop
        out.append(",".join(["E", str(e.id)] + [str(n) for n in e.nodes]))
    attrs: List[Tuple[str, Callable]] = [("ETYPE", lambda e: e.etype), ("EINT", lambda e: e.eint),
                                         ("THICK", lambda e: e.thick), ("KI", lambda e: tuple(e.ki)),
                                         ("KJ", lambda e: tuple(e.kj))]
    for cmd, getter in attrs:
        runs = [(a, b, v) for a, b, v in _runs(ids, lambda i: getter(g.elements[i]))
                if (any(v) if isinstance(v, tuple) else v)]
        if not runs:
            continue
        values = [v for _, _, v in runs]
        if not _accepts(cmd, cur, values):
            cur = T if _accepts(cmd, T, values) else next(t for t in _CARRIERS[cmd] if _accepts(cmd, t, values))
            out.append(f"MTYPE,{g.id},{type_name(cur)}")
        for a, b, v in runs:
            if isinstance(v, tuple):
                out.append(",".join([cmd, str(a), str(b), "1"] + [str(x) for x in v]))
            else:
                out.append(_cmd(cmd, a, b, 1, v))
    if cur != T:
        out.append(f"MTYPE,{g.id},{type_name(T)}")
    if any(ln.startswith("MTYPE,") for ln in out):
        out.insert(1, f"* group {g.id} ({type_name(T)}) holds data of an earlier type, kept by MTYPE")
    return out


def _groups(m: SSIModel) -> List[str]:
    out: List[str] = []
    state: Dict[str, Optional[int]] = {"mact": None, "ract": None}
    for gid in sorted(m.groups):
        out += _group_lines(m.groups[gid], state)
    last_group = max(m.groups) if m.groups else None
    if m.group_active is not None and m.group_active != last_group and m.group_active in m.groups:
        out.append(f"GROUP,{m.group_active}")
    elif m.group_active is None and m.groups:
        # no active group (GDEL of the active group, spec 08 section 5.12): GROUP always activates,
        # and only GDEL of the active group deactivates -- so a temporary group is created and deleted
        tmp = last_group + 1
        out.append("* no group is active")
        out += [f"GROUP,{tmp},{type_name(_SOLID)}", f"GDEL,{tmp}"]
    cur_mact, cur_ract = state["mact"], state["ract"]
    if (cur_mact if cur_mact is not None else 1) != m.mact:
        out.append(f"MACT,{m.mact}")
    if (cur_ract if cur_ract is not None else 1) != m.ract:
        out.append(f"RACT,{m.ract}")
    return out


def _masses(m: SSIModel) -> List[str]:
    out = [_cmd("MT", k, *v) for k, v in sorted(m.tmass.items())]
    out += [_cmd("MR", k, *v) for k, v in sorted(m.rmass.items())]
    eligible = sorted(set(m.nodes) | set(m.tmass) | set(m.rmass))
    for a, b, u in _runs(eligible, m.mass_unit):
        if u != m.MUNITS_DEFAULT:
            out.append(f"MUNITS,{a},{b},1,{u}")
    return out


def _loads(m: SSIModel) -> List[str]:
    out = [_cmd("F", k, *v.factor, *v.arrival) for k, v in sorted(m.forces.items())]
    out += [_cmd("MM", k, *v.factor, *v.arrival) for k, v in sorted(m.moments.items())]
    for _, rec in m.options.entries("SYMM"):
        out.append(record_line(rec))
    return out


def _freqs(m: SSIModel) -> List[str]:
    out = []
    for s, nums in sorted(m.freq_sets.items()):
        out.append(f"FREQ,{s},0")
        for ch in _chunks(nums, 10):
            out.append(_cmd("FREQ", s, *ch))
    return out


# -- list and request commands of the option section
def _w_damp(m: SSIModel) -> List[str]:
    if not m.damp:
        return []
    return ["DAMP,0"] + [_cmd("DAMP", *ch) for ch in _chunks(m.damp, 10)]


def _w_topl(m: SSIModel) -> List[str]:
    if not m.topl:
        return []
    return ["TOPL,0"] + [_cmd("TOPL", *ch) for ch in _chunks(m.topl, 20)]


def _w_amp(m: SSIModel) -> List[str]:
    out = []
    for no, vals in sorted(m.amp.items()):
        out.append(f"AMP,{no},0")
        for ch in _chunks(vals, 100, avoid_leading_zero=True):
            out.append(_cmd("AMP", no, *ch))
    return out


def _list_line(head: List[str], ids: Sequence[int]) -> List[str]:
    """``head`` + compressed id list, split over several lines only if longer than L14 allows."""
    toks = compress_ids(ids)
    line = ",".join(head + toks)
    if len(line) <= MAX_LINE:
        return [line]
    out, cur = [], []
    for t in toks:
        if len(",".join(head + cur + [t])) > MAX_LINE and cur:
            out.append(",".join(head + cur))
            cur = []
        cur.append(t)
    if cur:
        out.append(",".join(head + cur))
    return out


def _w_nout(m: SSIModel) -> List[str]:
    if not m.nout:
        return []
    out = ["NOUT,0"]
    for r in m.nout:
        out += _list_line(["NOUT", str(r.dir)] + [str(c) for c in r.codes], r.nodes)
    return out


def _w_eout(m: SSIModel) -> List[str]:
    if not m.eout:
        return []
    out = ["EOUT,0"]
    for r in m.eout:
        out += _list_line(["EOUT"] + [str(c) for c in r.codes] + [str(r.group)], r.elements)
    return out


def _w_rdnd(m: SSIModel) -> List[str]:
    if not m.rdnd:
        return []
    return ["RDND,0"] + [",".join(["RDND", str(r.node)] + [str(f) for f in r.flags]) for r in m.rdnd]


_BUILTIN = {"DAMP": _w_damp, "TOPL": _w_topl, "AMP": _w_amp, "NOUT": _w_nout, "EOUT": _w_eout, "RDND": _w_rdnd}


def _option_lines(m: SSIModel, name: str) -> List[str]:
    if name in WRITERS:
        return list(WRITERS[name](m))
    if name in _BUILTIN:
        return _BUILTIN[name](m)
    out: List[str] = []
    rec = m.options.record(name)
    if rec is not None:
        if not (rec.FIELDS and name in XCOMMANDS and rec.is_default()):
            out.append(record_line(rec))
    for _, r in m.options.entries(name):
        out.append(record_line(r))
    if name in m.options.strings:
        spec = _spec(name)
        val = m.options.strings[name]
        out.append(f"{name},{val}" if (spec and spec.text_from == 1) else join_command(name, [val]))
    return out


# --------------------------------------------------------------------------------------
# public API
# --------------------------------------------------------------------------------------
def write_pre(model: SSIModel, *, mdl: bool = False, afwr: bool = False, filename: str = "") -> Tuple[str, List[str]]:
    """Return ``(text, notes)``: the ``.pre`` text of ``model`` and warnings about inexact parts.

    ``mdl`` writes the MDL line (Options > Write "MDL command in *.Pre"); ``afwr`` appends AFWR.
    """
    m = model
    lines: List[str] = ["*" * 72,
                        f"* Written by {PRODUCT} {__version__} (ACS SASSI V3 command language)",
                        f"* Reload with  INP,{filename or (m.name + '.pre' if m.name else '<this file>')}",
                        "*" * 72]

    def section(title: str, body: List[str]) -> None:
        if body:
            lines.append(f"* {title}")
            lines.extend(body)

    section("Model identification", _identification(m, mdl))
    node_lines, notes = _nodes(m)
    section("Coordinate systems and nodes", node_lines)
    section("Boundary conditions", _fixities(m))
    section("Interaction nodes", _interaction(m))
    mat, lay, real = _tables(m)
    section("Material table", mat)
    section("Soil layer table", lay)
    section("Real property table", real)
    section("Groups and elements", _groups(m))
    section("Masses", _masses(m))
    section("Loads", _loads(m))
    section("Frequencies", _freqs(m))
    covered = set(_HANDLED)
    for title, names in OPTION_ORDER:
        body: List[str] = []
        for name in names:
            body += _option_lines(m, name)
            covered.add(name)
        section(title, body)
    covered.update(AFWRITE_SECTION[1])
    others = [n for n in m.options.names() if n not in covered]
    body = []
    for n in others:
        body += _option_lines(m, n)
    section("Other stored options", body)
    for title, func in SECTIONS:
        section(title, list(func(m)))
    body = []
    for n in AFWRITE_SECTION[1]:
        body += _option_lines(m, n)
    section(AFWRITE_SECTION[0], body)
    if afwr:
        lines.append("AFWR")
    return "\n".join(lines) + "\n", notes
