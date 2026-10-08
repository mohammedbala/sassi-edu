"""The Properties panel (requirements 7.20, D-W6-16): the selected nodes and elements of a model, their
properties edited with command text, and the input file that built the model kept in step.

* :func:`selection` -- what the panel shows (``GET /api/selection``): the nodes of NODESEL and the elements
  of ELEMSEL with their properties, the materials, soil layers and sections to choose from, and the input
  file that built the model in this session (:class:`sassi.ui.runview.RunTracker` records it);
* :func:`commands` -- the command lines of a list of edits: ``N`` (coordinates), ``D`` (fixities), ``INT``
  (interaction flag), ``MT`` / ``MR`` (masses), ``MSET`` / ``RSET`` / ``THICK`` / ``ETYPE`` / ``EINT`` /
  ``KI`` / ``KJ`` (elements, the group activated with ``GROUP`` and the active group restored);
* :func:`patch_input` -- the input file changed so that it builds the edited model: the explicit ``N``,
  ``MT`` or ``MR`` line of a node rewritten in place when the file has exactly one; every other change in an
  *edits* section placed after the part of the file that last changed the model (so before its CHECK,
  AFWRITE and module runs, and before a CPMODEL copies it).

An edit is a dict::

    {"op": "xyz", "nodes": [...], "x": v | None, "y": ..., "z": ...}         global coordinates (None: kept)
    {"op": "fix", "nodes": [...], "dofs": {"UX": 1, "ROTZ": 0, ...}}        1 fixed, 0 free
    {"op": "interaction", "nodes": [...], "value": 1 | 0}
    {"op": "mass" | "rmass", "nodes": [...], "values": [a | None, b | None, c | None]}
    {"op": "elem", "group": g, "elements": [...], "attr": "mat" | "prop" | "thick" | "etype" | "eint", "value": v}
    {"op": "release", "group": g, "elements": [...], "end": "I" | "J", "codes": [k1, ..., k6]}
"""
from __future__ import annotations

import math
import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ..model.entities import DOF_NAMES, type_name
from ..model.values import fmt_num

#: the header of an edits section (one per model) in an input file
EDITS_HEADER = "* ---------------------------------------------------------------------- edits (Properties panel)"
#: commands before which an edits section goes when the file's history is unknown
STOP = {"CHECK", "AFWRITE", "AFWRBAT", "WRITE", "SAVE", "CPMODEL"}
ELEM_COMMANDS = {"mat": "MSET", "prop": "RSET", "thick": "THICK", "etype": "ETYPE", "eint": "EINT"}
MAX_LIST = 2000


class EditError(ValueError):
    """An edit the panel cannot turn into commands (bad value, unknown node ...)."""


# ====================================================================== what the panel shows
def _num(x: float) -> str:
    """A number as an input file writes it: at most 10 significant digits, no trailing zeros."""
    x = float(x)
    if x == int(x) and abs(x) < 1e15:
        return str(int(x))
    return fmt_num(float(f"{x:.10g}"))


def selection(interp, number: Optional[int] = None) -> Dict[str, Any]:
    """The selected nodes and elements of model ``number`` (default: the active model) with their properties."""
    num = interp.active_model if number is None else int(number)
    m = interp.models.get(num)
    out: Dict[str, Any] = {"model": num if m is not None else None, "name": (m.name or "") if m else "",
                           "nodes": [], "groups": [], "materials": [], "layers": [], "sections": [], "springs": []}
    if m is None:
        return out
    nsel = [int(n) for n in m.ui_state.get("nodesel", []) if int(n) in m.nodes][:MAX_LIST]
    esel = [(int(g), int(e)) for g, e in m.ui_state.get("elemsel", [])
            if int(g) in m.groups and int(e) in m.groups[int(g)].elements][:MAX_LIST]
    conn: Dict[int, List[List[int]]] = {n: [] for n in nsel}
    if nsel:
        for g, e in m.iter_elements():
            for n in e.nodes:
                if n in conn and [g.id, e.id] not in conn[n]:
                    conn[n].append([g.id, e.id])
        ids, xyz = m.global_coordinates(nsel)
        gx = {int(i): [float(v) for v in p] for i, p in zip(ids, xyz)}
        for n in nsel:
            nd = m.nodes[n]
            out["nodes"].append({"id": n, "xyz": gx[n], "csys": nd.csys, "fix": [int(v) for v in nd.fix],
                                 "interaction": 0 in nd.flags, "mass": _floats(m.tmass.get(n)),
                                 "rmass": _floats(m.rmass.get(n)), "munits": int(m.mass_units.get(n, 1)),
                                 "elements": conn[n][:50]})
    by_group: Dict[int, List[int]] = {}
    for g, e in esel:
        by_group.setdefault(g, []).append(e)
    for gid in sorted(by_group):
        g = m.groups[gid]
        els = []
        for eid in sorted(by_group[gid]):
            e = g.elements[eid]
            els.append({"id": eid, "nodes": [int(n) for n in e.nodes], "mat": e.mat, "prop": e.prop,
                        "thick": float(e.thick), "etype": e.etype, "eint": e.eint, "ki": list(e.ki), "kj": list(e.kj)})
        out["groups"].append({"id": gid, "type": g.type, "type_name": type_name(g.type), "title": g.title,
                              "elements": els})
    out["materials"] = [{"id": k, "text": f"E {_num(v.val1)}, nu {_num(v.val2)}, unit weight {_num(v.weight)}"
                         if v.mtype == 1 else f"{_num(v.val1)}, {_num(v.val2)}, unit weight {_num(v.weight)}"}
                        for k, v in sorted(m.materials.items())]
    out["layers"] = [{"id": k, "text": f"Vs {_num(v.vs)}, Vp {_num(v.vp)}, thick {_num(v.thick)}"}
                     for k, v in sorted(m.layers.items())]
    out["sections"] = [{"id": k, "text": f"A {_num(v.axial)}, I2 {_num(v.flex2)}, I3 {_num(v.flex3)}"}
                       for k, v in sorted(m.sections.items())]
    out["springs"] = [{"id": k, "text": f"Kx {_num(v.scx)}, Ky {_num(v.scy)}, Kz {_num(v.scz)}"}
                      for k, v in sorted(m.springs.items())]
    out["active_group"] = m.group_active
    out["csys_active"] = m.csys_active
    return out


def _floats(v) -> Optional[List[float]]:
    return None if v is None else [float(x) for x in v]


# ====================================================================== edits -> commands
def ranges(ids: Iterable[int]) -> List[Tuple[int, int]]:
    """Consecutive runs of the sorted ids: ``[1, 2, 3, 7]`` -> ``[(1, 3), (7, 7)]``."""
    s = sorted(set(int(i) for i in ids))
    out: List[Tuple[int, int]] = []
    for i in s:
        if out and i == out[-1][1] + 1:
            out[-1] = (out[-1][0], i)
        else:
            out.append((i, i))
    return out


def _value(v: Any, what: str, kind=float, positive: bool = False):
    try:
        x = kind(v)
    except (TypeError, ValueError):
        raise EditError(f"{what}: {v!r} is not a number") from None
    if kind is float and not math.isfinite(x):
        raise EditError(f"{what}: {v!r} is not a finite number")
    if positive and x <= 0:
        raise EditError(f"{what} must be > 0")
    return x


def _nodes(m, ids: Sequence[Any]) -> List[int]:
    out = []
    for n in ids:
        n = int(n)
        if n not in m.nodes:
            raise EditError(f"node {n} is not defined")
        out.append(n)
    return out


def node_lines(m, edit: Dict[str, Any]) -> List[str]:
    """The command lines of one node edit (no CSYS wrapping: see :func:`commands`)."""
    op = edit.get("op")
    nodes = _nodes(m, edit.get("nodes") or [])
    if not nodes:
        return []
    if op == "xyz":
        new = [None if edit.get(k) in (None, "") else _value(edit.get(k), k.upper()) for k in ("x", "y", "z")]
        if all(v is None for v in new):
            return []
        ids, xyz = m.global_coordinates(nodes)
        cur = {int(i): list(p) for i, p in zip(ids, xyz)}
        return [f"N,{n}," + ",".join(_num(new[k] if new[k] is not None else cur[n][k]) for k in range(3)) for n in nodes]
    if op == "fix":
        dofs = {str(k).upper(): int(v) for k, v in (edit.get("dofs") or {}).items()}
        bad = [k for k in dofs if k not in DOF_NAMES]
        if bad:
            raise EditError(f"unknown DOF label(s) {bad}")
        out = []
        for val in (1, 0):
            labels = [k for k in DOF_NAMES if dofs.get(k) == val]
            if labels:
                out += [f"D,{a},{b},1,{val}," + ",".join(labels) for a, b in ranges(nodes)]
        return out
    if op == "interaction":
        val = 1 if int(edit.get("value", 1)) else 0
        return [f"INT,{a},{b},1,{val},0" for a, b in ranges(nodes)]
    if op in ("mass", "rmass"):
        vals = edit.get("values") or [None, None, None]
        new = [None if v in (None, "") else _value(v, "mass") for v in list(vals)[:3]]
        if all(v is None for v in new):
            return []
        store, cmd = (m.tmass, "MT") if op == "mass" else (m.rmass, "MR")
        out = []
        for n in nodes:
            cur = list(store.get(n, [0.0, 0.0, 0.0])) + [0.0, 0.0, 0.0]
            out.append(f"{cmd},{n}," + ",".join(_num(new[k] if new[k] is not None else cur[k]) for k in range(3)))
        return out
    raise EditError(f"unknown node edit {op!r}")


def element_lines(m, edit: Dict[str, Any]) -> Tuple[int, List[str]]:
    """``(group, command lines)`` of one element edit (the group is activated by :func:`commands`)."""
    gid = int(edit.get("group", 0))
    g = m.groups.get(gid)
    if g is None:
        raise EditError(f"group {gid} does not exist")
    els = [int(e) for e in edit.get("elements") or []]
    missing = [e for e in els if e not in g.elements]
    if missing:
        raise EditError(f"element(s) {missing} of group {gid} are not defined")
    if not els:
        return gid, []
    if edit.get("op") == "elem":
        attr = edit.get("attr")
        if attr not in ELEM_COMMANDS:
            raise EditError(f"unknown element property {attr!r}")
        if attr == "thick":
            v = _num(_value(edit.get("value"), "thickness", positive=True))
        else:
            v = str(_value(edit.get("value"), attr, kind=int))
        return gid, [f"{ELEM_COMMANDS[attr]},{a},{b},1,{v}" for a, b in ranges(els)]
    if edit.get("op") == "release":
        end = str(edit.get("end", "I")).upper()
        codes = [1 if int(c) else 0 for c in (edit.get("codes") or [])][:6]
        if end not in ("I", "J") or len(codes) != 6:
            raise EditError("a release edit needs end I or J and six codes")
        return gid, [f"K{end},{a},{b},1," + ",".join(str(c) for c in codes) for a, b in ranges(els)]
    raise EditError(f"unknown element edit {edit.get('op')!r}")


def commands(m, edits: Sequence[Dict[str, Any]]) -> List[str]:
    """The command lines of ``edits`` on model ``m``: node edits first (N in global coordinates: wrapped in
    ``CSYS,0`` / ``CSYS,<active>`` when a local system is active), then element edits group by group (``GROUP,g``
    before them, the active group activated again after them)."""
    nlines: List[str] = []
    glines: Dict[int, List[str]] = {}
    for ed in edits:
        if ed.get("op") in ("elem", "release"):
            gid, ls = element_lines(m, ed)
            glines.setdefault(gid, []).extend(ls)
        else:
            nlines.extend(node_lines(m, ed))
    out: List[str] = []
    if any(ln.startswith("N,") for ln in nlines) and m.csys_active != 0:
        out.append("CSYS,0")
        out.extend(nlines)
        out.append(f"CSYS,{m.csys_active}")
    else:
        out.extend(nlines)
    last = m.group_active
    for gid, ls in glines.items():
        if not ls:
            continue
        if gid != last:
            out.append(f"GROUP,{gid}")
            last = gid
        out.extend(ls)
    if glines and m.group_active is not None and last != m.group_active and any(glines.values()):
        out.append(f"GROUP,{m.group_active}")
    return out


# ====================================================================== the input file
_PLAIN_NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


def _head_args(line: str) -> Tuple[str, List[str]]:
    """``(canonical command name, argument tokens)`` of a plain command line ("" for comments)."""
    from ..prep.lexer import is_comment, split_head
    from ..prep.registry import lookup
    s = line.strip()
    if not s or is_comment(s):
        return "", []
    name, rest = split_head(s)
    try:
        spec = lookup(name)
    except Exception:                            # noqa: BLE001
        spec = None
    return (spec.name if spec is not None else name.upper()), [t.strip() for t in rest.split(",")] if rest else []


def _is_break(line: str) -> bool:
    s = line.strip()
    return not s or s == "*"


def _single_lines(lines: Sequence[str], cmd: str) -> Dict[int, List[int]]:
    """Plain ``cmd,<id>,<numbers>`` lines of the file by node id (line indices)."""
    out: Dict[int, List[int]] = {}
    for i, ln in enumerate(lines):
        head, args = _head_args(ln)
        if head != cmd or not args or not args[0].isdigit():
            continue
        if not all(_PLAIN_NUM.match(a) for a in args[1:4] if a):
            continue
        out.setdefault(int(args[0]), []).append(i)
    return out


def _rewrite(line: str, values: Sequence[str]) -> str:
    """Replace arguments 2-4 (the coordinates or the three masses) of a plain command line."""
    head, rest = line.strip().split(",", 1)
    args = rest.split(",")
    while len(args) < 4:
        args.append("")
    for k, v in enumerate(values):
        args[1 + k] = v
    indent = line[:len(line) - len(line.lstrip())]
    return indent + head + "," + ",".join(a.strip() for a in args).rstrip(",")


def anchor_fallback(lines: Sequence[str]) -> int:
    """Where an edits section goes when the file's history is unknown: before the section of the first CHECK,
    AFWRITE, WRITE, SAVE, CPMODEL or RUN<MODULE> line; else at the end.  Returns the insertion index."""
    for i, ln in enumerate(lines):
        head, _ = _head_args(ln)
        if head in STOP or (head.startswith("RUN") and len(head) > 3):
            j = i
            while j > 0 and not _is_break(lines[j - 1]):
                j -= 1
            return j
    return len(lines)


def _section_end(lines: Sequence[str], i: int) -> int:
    """The index after the last line of the section that holds line ``i``."""
    j = i + 1
    while j < len(lines) and not _is_break(lines[j]):
        j += 1
    return j


def patch_input(text: str, m, edits: Sequence[Dict[str, Any]], anchor: Optional[int] = None) -> Tuple[str, Dict[str, Any]]:
    """The input file ``text`` changed so that it builds model ``m`` with ``edits`` (see the module docstring).

    ``anchor``: index of the line of the file that last changed the model when the file ran (the edits section
    goes after its section); None: :func:`anchor_fallback`.  Returns ``(new text, {"changed": [line indices of
    the new text], "inplace": k, "added": k})``.  Call it *before* the edits are applied to the model (the
    current coordinates and masses complete partial edits)."""
    nl = "\n" if not text.endswith("\r\n") else "\r\n"
    lines = text.splitlines()
    changed_in_place: List[int] = []
    block: List[str] = []
    rest: List[Dict[str, Any]] = []
    explicit = {"N": _single_lines(lines, "N"), "MT": _single_lines(lines, "MT"), "MR": _single_lines(lines, "MR")}
    for ed in edits:
        op = ed.get("op")
        if op in ("xyz", "mass", "rmass"):
            cmd = {"xyz": "N", "mass": "MT", "rmass": "MR"}[op]
            for ln in node_lines(m, ed):
                n = int(ln.split(",")[1])
                where = explicit[cmd].get(n, [])
                local = cmd == "N" and m.nodes[n].csys != 0
                if len(where) == 1 and not local:
                    i = where[0]
                    lines[i] = _rewrite(lines[i], ln.split(",")[2:5])
                    changed_in_place.append(i)
                else:
                    block.append(ln)
        else:
            rest.append(ed)
    if rest:
        block.extend(commands(m, rest))
    if any(ln.startswith("N,") for ln in block) and any(_head_args(x)[0] == "CSYS" for x in lines):
        block = ["CSYS,0"] + block
    added: List[int] = []
    if block:
        start = _find_block(lines, m.name)
        if start is not None:
            end = _section_end(lines, start)
            body = lines[start + 2:end]
            for ln in block:                     # a node's N / MT / MR line replaces its previous one
                key = ln.split(",")[:2] if ln.split(",")[0] in ("N", "MT", "MR") else None
                if key is not None:
                    body = [x for x in body if x.split(",")[:2] != key]
                body.append(ln)
            lines[start + 2:end] = body
            shift = len(body) - (end - start - 2)
            changed_in_place = [i + shift if i >= end else i for i in changed_in_place]
            added = list(range(start + 2 + len(body) - len(block), start + 2 + len(body)))
        else:
            at = anchor_fallback(lines) if anchor is None else _section_end(lines, max(0, min(anchor, len(lines) - 1)))
            sec = ["*", EDITS_HEADER, f"* model {m.name or '?'}: changes made in the GUI's Properties panel"] + block
            if at < len(lines) and not _is_break(lines[at]):
                sec.append("*")
            lines[at:at] = sec
            changed_in_place = [i + len(sec) if i >= at else i for i in changed_in_place]
            added = list(range(at + 3, at + 3 + len(block)))
    new = nl.join(lines) + (nl if text.endswith(("\n", "\r")) or not text else "")
    return new, {"changed": sorted(set(changed_in_place + added)), "inplace": len(changed_in_place), "added": len(block)}


def _find_block(lines: Sequence[str], name: str) -> Optional[int]:
    """Index of the header of the edits section of model ``name`` (the header and its "* model" line)."""
    for i in range(len(lines) - 1):
        if lines[i].strip() == EDITS_HEADER and lines[i + 1].strip().startswith(f"* model {name or '?'}:"):
            return i
    return None
