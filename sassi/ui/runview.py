"""The Run view (requirements 7.20, D-W6-15): what an INP run is doing, line by line, beside the model.

While an INP / macro file runs in the GUI session, :class:`RunTracker` (hooked to the interpreter's
``input_frame`` and ``line_done``) pushes ``run`` events:

* ``{event: "file", path, name, kind, lines, depth}`` when a file starts (its lines, for the listing in which
  the browser shades the section being processed and marks the line), ``{event: "end", path, depth}`` when it
  ends;
* ``{event: "line", path, line, total, cmd, name, focus, text}`` after a line that changed the active model:
  ``focus`` holds what it added or changed (:func:`diff`: nodes, elements, soil layers, interaction nodes,
  fixities, masses, output requests), ``text`` a one-line summary ("nodes 2-9 added (8)");
* ``{event: "scene", data}``: a snapshot of the active model drawn the way MODELPLOT draws it, with the
  free-field soil around it (SHOWSOIL) -- or its soil column (LAYERPLOT) while it has no elements --, at most
  every :data:`SCENE_INTERVAL` seconds while the geometry changes, and at the end of every file.

The visitor can ask to *watch* the build: ``pace`` seconds after every line that changed the model (the run
is otherwise as fast as without the view).  Module runs report their steps through the activity panel's
``progress`` / ``job`` events (``ctx.announce``); the browser maps a step to the part of the model it works on
with the lists of :func:`requests` (interaction nodes, output nodes and elements, loaded nodes).

Everything here is a display aid: a failure is swallowed and never changes or stops the run.
"""
from __future__ import annotations

import time
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from ..plotting.state import PlotError, jsonable

#: seconds between two model snapshots while the geometry changes (and between two line events in real time)
SCENE_INTERVAL = 0.3
LINE_INTERVAL = 0.08
#: watch speeds the browser offers (seconds of pause after a line that changed the model)
PACES = (0.0, 0.12, 0.4)
#: at most this many ids per list in an event (the summary text still counts them all)
MAX_IDS = 4000
#: interaction flags and fixities are compared node by node after every line up to this model size; above
#: it only after the commands that set them
FULL_NODE_DIFF = 4000
FLAG_COMMANDS = {"INT", "INTGEN", "D", "FIXROT", "EXCSTRCHK", "FILLPOOL", "MERGEPOOL", "CONVERT", "NCOM",
                 "RMVUNUSED", "FOREACH", "MACRO", "INP"}


# ====================================================================== what a line changed
class Footprint:
    """The state of the active model that :func:`diff` compares (cheap: key sets, and the node flags of
    models up to :data:`FULL_NODE_DIFF` nodes or after the commands that set them)."""

    def __init__(self, interp, flags: bool = True):
        self.model = interp.active_model
        m = interp.models.get(self.model)
        self.m = m
        if m is None:
            return
        self.nodes: Set[int] = set(m.nodes)
        self.groups: Dict[int, Set[int]] = {g.id: set(g.elements) for g in m.groups.values()}
        self.layers = set(m.layers)
        self.topl = list(m.topl)
        self.materials = set(m.materials)
        self.sections = set(m.sections)
        self.mass = set(m.tmass) | set(m.rmass)
        self.loads = set(m.forces) | set(m.moments)
        self.n_nout, self.n_eout, self.n_rdnd = len(m.nout), len(m.eout), len(m.rdnd)
        self.flags: Optional[Dict[int, Tuple[Tuple[int, ...], Tuple[int, ...]]]] = None
        if flags:
            self.flags = {k: (tuple(sorted(n.flags)), tuple(n.fix)) for k, n in m.nodes.items() if n.flags or any(n.fix)}


def _ids(xs: Iterable[int]) -> List[int]:
    return sorted(int(x) for x in xs)[:MAX_IDS]


def ranges(ids: Sequence[int], limit: int = 6) -> str:
    """``[1, 2, 3, 7, 9, 10]`` -> ``"1-3, 7, 9-10"`` (at most ``limit`` runs, then ``...``)."""
    s = sorted(set(int(i) for i in ids))
    out: List[str] = []
    k = 0
    while k < len(s):
        j = k
        while j + 1 < len(s) and s[j + 1] == s[j] + 1:
            j += 1
        out.append(str(s[k]) if j == k else f"{s[k]}-{s[j]}")
        k = j + 1
        if len(out) >= limit and k < len(s):
            out.append("...")
            break
    return ", ".join(out)


def diff(before: Footprint, after: Footprint) -> Dict[str, Any]:
    """What changed between two footprints of the same model: ``{model, nodes, elements [[g, e]], layers, topl,
    interaction, fixed, mass, loads, groups, requests, all}`` (empty lists left out; ``all`` when another model
    became active or the model was replaced: the browser then highlights nothing and redraws)."""
    if after.m is None:
        return {}
    if before.m is None or before.model != after.model or before.m is not after.m:
        return {"model": after.model, "all": True}
    m = after.m
    f: Dict[str, Any] = {"model": after.model}
    new_nodes = after.nodes - before.nodes
    if new_nodes:
        f["nodes"] = _ids(new_nodes)
    elems: List[List[int]] = []
    groups: List[int] = []
    for gid, keys in after.groups.items():
        old = before.groups.get(gid)
        if old is None:
            groups.append(gid)
            old = set()
        for e in sorted(keys - old):
            elems.append([gid, int(e)])
    new_mats = after.materials - before.materials
    changed_mats = set(new_mats)
    new_secs = after.sections - before.sections
    if changed_mats or new_secs:
        # a material or section defined after the elements that use it: those elements
        for g, e in m.iter_elements():
            if e.mat in changed_mats or e.prop in new_secs:
                elems.append([g.id, e.id])
    if elems:
        seen: Set[Tuple[int, int]] = set()
        f["elements"] = [x for x in elems if not (tuple(x) in seen or seen.add(tuple(x)))][:MAX_IDS]
    if groups:
        f["groups"] = sorted(groups)
    layers = after.layers - before.layers
    if after.topl != before.topl:
        f["topl"] = True
        layers |= set(after.topl[len(before.topl):]) if after.topl[:len(before.topl)] == before.topl else set(after.topl)
    if layers:
        f["layers"] = _ids(layers)
    if new_mats:
        f["materials"] = _ids(new_mats)
    if new_secs:
        f["sections"] = _ids(new_secs)
    if before.flags is not None and after.flags is not None:
        inter, fixed = [], []
        for k in set(before.flags) | set(after.flags):
            b, a = before.flags.get(k, ((), (0,) * 6)), after.flags.get(k, ((), (0,) * 6))
            if a[0] != b[0] and 0 in a[0] and 0 not in b[0]:
                inter.append(k)
            if a[1] != b[1]:
                fixed.append(k)
        if inter:
            f["interaction"] = _ids(inter)
        if fixed:
            f["fixed"] = _ids(fixed)
    mass = after.mass - before.mass
    if mass:
        f["mass"] = _ids(mass)
    loads = after.loads - before.loads
    if loads:
        f["loads"] = _ids(loads)
    req: Dict[str, Any] = {}
    if after.n_nout > before.n_nout:
        req["nout"] = _ids(n for r in m.nout[before.n_nout:] for n in r.nodes)
    if after.n_eout > before.n_eout:
        req["eout"] = _eout_elements(m, m.eout[before.n_eout:])[:MAX_IDS]
    if after.n_rdnd > before.n_rdnd:
        req["rdnd"] = _ids(r.node for r in m.rdnd[before.n_rdnd:])
    if req:
        f["requests"] = req
    return f if len(f) > 1 else {}


def _eout_elements(m, reqs) -> List[List[int]]:
    out: List[List[int]] = []
    for r in reqs:
        g = m.groups.get(int(r.group))
        if g is None:
            continue
        ids = [int(e) for e in r.elements if int(e) > 0]
        out.extend([g.id, e] for e in (ids or sorted(g.elements)) if e in g.elements)
    return out


def describe(f: Dict[str, Any]) -> str:
    """One-line summary of a focus: ``"nodes 2-9 added (8)"``, ``"interaction nodes 1-81 (81)"`` ..."""
    if not f:
        return ""
    if f.get("all"):
        return f"model {f.get('model')}"
    parts: List[str] = []

    def add(what: str, ids: Sequence[Any], verb: str = "") -> None:
        flat = [i[1] if isinstance(i, (list, tuple)) else i for i in ids]
        parts.append(f"{what} {ranges(flat)}{(' ' + verb) if verb else ''} ({len(ids)})")

    if f.get("nodes"):
        add("nodes", f["nodes"], "added")
    if f.get("elements"):
        gs = sorted(set(g for g, _ in f["elements"]))
        add(f"group {gs[0]} elements" if len(gs) == 1 else "elements", f["elements"])
    elif f.get("groups"):
        parts.append("group " + ranges(f["groups"]))
    if f.get("layers"):
        parts.append(("soil profile: layers " if f.get("topl") else "soil layers ") + ranges(f["layers"]))
    if f.get("interaction"):
        add("interaction nodes", f["interaction"])
    if f.get("fixed"):
        add("fixed DOFs at nodes", f["fixed"])
    if f.get("mass"):
        add("masses at nodes", f["mass"])
    if f.get("loads"):
        add("loads at nodes", f["loads"])
    if f.get("materials") and not f.get("elements"):
        parts.append("material " + ranges(f["materials"]))
    if f.get("sections") and not f.get("elements"):
        parts.append("section " + ranges(f["sections"]))
    req = f.get("requests") or {}
    if req.get("nout"):
        add("output nodes", req["nout"])
    if req.get("eout"):
        add("output elements", req["eout"])
    if req.get("rdnd"):
        add("relative-displacement nodes", req["rdnd"])
    return "; ".join(parts)


# ====================================================================== the snapshot drawn
def requests(m) -> Dict[str, Any]:
    """The lists the browser maps the module steps onto: interaction nodes, output nodes (NOUT), output
    elements (EOUT), relative-displacement nodes (RDND), loaded nodes (F, MM)."""
    return {"interaction": _ids(k for k, n in m.nodes.items() if 0 in n.flags),
            "nout": _ids(set(n for r in m.nout for n in r.nodes)),
            "eout": _eout_elements(m, m.eout)[:MAX_IDS],
            "rdnd": _ids(set(r.node for r in m.rdnd)),
            "loads": _ids(set(m.forces) | set(m.moments))}


def snapshot(interp, state) -> Dict[str, Any]:
    """The active model as the Run view draws it: the MODELPLOT data (:func:`model_scene`, the default view)
    with the free-field soil around it (SHOWSOIL) when it has a soil profile, its nodes that no element uses
    yet (orientation nodes, a mesh before its elements: ``scene.free`` marks them) appended to the scene; its
    soil column (LAYERPLOT data) while it has no node; ``{"empty": true}`` without anything to draw.  Plus
    ``run`` (the lists of :func:`requests`) and ``name``."""
    import numpy as np
    from ..plotting.state import (layer_table, model_scene, rotation_matrix, soil_island, view_frame)
    n = interp.active_model
    m = interp.models.get(n)
    if m is None:
        return {"empty": True, "model": None, "name": ""}
    soil_ok = bool(m.topl and m.layers)
    if not m.nodes:
        if not soil_ok:
            return {"empty": True, "model": n, "name": m.name or ""}
        try:
            table = layer_table(m)
        except PlotError:
            return {"empty": True, "model": n, "name": m.name or ""}
        return jsonable({"kind": "LAYERPLOT", "family": "layer", "model": n, "name": m.name or "", "table": table,
                         "show": {}, "title": "", "run": requests(m)})
    v = state.defaults3d.copy()
    scene = model_scene(m, color_by=v.color_by, show_dof=v.show_dof, show_mass=v.show_mass, palettes=state.palettes)
    ids = [int(x) for x in scene["node_id"]]
    have = set(ids)
    free = sorted(k for k in m.nodes if k not in have)
    if free:
        fid, fxyz = m.global_coordinates(free)
        k0 = len(ids)
        scene["node_id"] = np.concatenate([np.asarray(scene["node_id"], dtype=np.int64), np.asarray(fid, dtype=np.int64)])
        scene["xyz"] = np.vstack([np.asarray(scene["xyz"], dtype=float).reshape(-1, 3), np.asarray(fxyz, dtype=float).reshape(-1, 3)])
        scene["node_visible"] = np.concatenate([np.asarray(scene["node_visible"], dtype=bool), np.ones(len(fid), dtype=bool)])
        scene["fix"] = np.vstack([np.asarray(scene["fix"]).reshape(-1, 6),
                                  np.asarray([m.nodes[int(x)].fix for x in fid], dtype=np.int64).reshape(-1, 6)])
        extra_int = [k0 + j for j, x in enumerate(fid) if 0 in m.nodes[int(x)].flags]
        scene["interaction"] = np.concatenate([np.asarray(scene["interaction"], dtype=np.int64),
                                               np.asarray(extra_int, dtype=np.int64)])
        lo, hi = scene["xyz"].min(axis=0), scene["xyz"].max(axis=0)
        scene["bbox"] = [float(lo[0]), float(hi[0]), float(lo[1]), float(hi[1]), float(lo[2]), float(hi[2])]
        scene["center"] = [float(x) for x in 0.5 * (lo + hi)]
    scene["free"] = list(range(len(ids), len(ids) + len(free)))
    soil, note = None, ""
    if soil_ok and scene["n_elements"]:
        try:
            soil = soil_island(m, rotation_matrix(v.rx, v.ry, v.rz), cut=v.soil_cut, margin=v.soil_margin,
                               depth=v.soil_depth)
        except PlotError as exc:
            note = str(exc)
    if soil is not None:
        b, sb = scene["bbox"], soil["bbox"]
        b = [min(b[0], sb[0]), max(b[1], sb[1]), min(b[2], sb[2]), max(b[3], sb[3]), min(b[4], sb[4]), max(b[5], sb[5])]
        scene["bbox"] = [float(x) for x in b]
        scene["center"] = [0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3]), 0.5 * (b[4] + b[5])]
    v.show_soil = soil is not None
    d = {"kind": "MODELPLOT", "family": "3d", "model": n, "name": m.name or "", "title": "", "view": v.to_dict(),
         "shader": state.shader.to_dict(), "scene": scene, "soil": soil,
         "camera": view_frame(v, scene, extra=soil["corners"] if soil is not None else None)}
    if note:
        d["soil_note"] = note
    d = jsonable(d)
    d["run"] = requests(m)
    return d


# ====================================================================== the tracker
class RunTracker:
    """Hooks of the GUI session's interpreter (``input_frame``, ``line_done``) that push the ``run`` events.

    ``push(event_type, **data)`` is the session's event log; ``state()`` returns the plot state (palettes,
    defaults of the 3D view); ``enabled`` is set while a Run view is open in the browser (snapshots are made
    only then); ``pace`` (seconds) slows the lines that change the model down to watching speed.
    """

    def __init__(self, interp, push: Callable[..., Any], state: Callable[[], Any],
                 sleep: Callable[[float], None] = time.sleep, clock: Callable[[], float] = time.monotonic):
        self.interp = interp
        self.push = push
        self.state = state
        self.sleep = sleep
        self.clock = clock
        self.enabled = False
        self.pace = 0.0
        self._before: Optional[Footprint] = None
        self._scene_t = -1e9
        self._line_t = -1e9
        self._dirty = False
        self._cost = 0.0                         # seconds the last snapshot took
        self._last_scene = 0                     # sequence number of the newest snapshot event
        #: model number -> (path of the input file, index of its line) that last changed the model: the Properties
        #: panel writes its edits into that file, after that line's section (sassi/ui/properties.py)
        self.model_lines: Dict[int, Tuple[str, int]] = {}
        self.replace: Optional[Callable[..., Any]] = None   # EventLog.replace (the session sets it)
        self._pending: Optional[Dict[str, Any]] = None
        interp.input_frame = self.on_frame
        interp.line_done = self.on_line_done
        prev = interp.progress

        def progress(i: int, total: int, source: str, line: str) -> None:
            self.on_line_start(line)
            if prev is not None:
                prev(i, total, source, line)
        interp.progress = progress

    # ------------------------------------------------------------------ settings (POST /api/runview)
    def configure(self, enabled: Optional[bool] = None, pace: Optional[float] = None) -> Dict[str, Any]:
        if enabled is not None:
            self.enabled = bool(enabled)
        if pace is not None:
            self.pace = max(0.0, min(2.0, float(pace)))
        return {"enabled": self.enabled, "pace": self.pace}

    # ------------------------------------------------------------------ hooks
    def on_frame(self, event: str, frame, lines: Sequence[str]) -> None:
        path = str(frame.path) if frame.path is not None else ""
        depth = len(self.interp.inputs)
        if event == "start":
            self.push("run", event="file", path=path, name=frame.name, kind=frame.kind, depth=depth,
                      lines=[str(x) for x in lines])
            if depth == 1:
                self._scene_t = -1e9
                self._dirty = True
        else:
            self._flush_line()
            if self._dirty and self.enabled:
                self._scene(force=True)
            self.push("run", event="end", path=path, name=frame.name, depth=depth + 1)

    def on_line_start(self, line: str) -> None:
        head = _head(line)
        if not head:                             # a comment, a blank line or an unknown command
            self._before = None
            return
        m = self.interp.models.get(self.interp.active_model)
        flags = m is None or len(m.nodes) <= FULL_NODE_DIFF or head in FLAG_COMMANDS
        self._before = Footprint(self.interp, flags=flags)

    def on_line_done(self, i: int, total: int, source: str, line: str, ok: bool) -> None:
        before, self._before = self._before, None
        if before is None:
            return
        after = Footprint(self.interp, flags=before.flags is not None)
        f = diff(before, after)
        if not f:
            return
        outer = self.interp.inputs[0] if self.interp.inputs else None
        if outer is not None and outer.path is not None and f.get("model") is not None:
            self.model_lines[int(f["model"])] = (str(outer.path), outer.line - 1)
        frame = self.interp.inputs[-1] if self.interp.inputs else None
        ev = {"path": str(frame.path) if frame is not None and frame.path is not None else "",
              "line": i, "total": total, "cmd": line.strip()[:160], "name": _head(line), "focus": f,
              "text": describe(f), "ok": bool(ok)}
        self._dirty = True
        now = self.clock()
        if self.pace > 0 or now - self._line_t >= LINE_INTERVAL:
            self._pending = None
            self._line_t = now
            self.push("run", event="line", **ev)
        else:
            self._pending = ev                  # the last line of a fast burst is still shown
        # snapshots: at most every SCENE_INTERVAL s in real time, after every line when watching (but never
        # more often than 2.5 times their own cost, so that a large model is not slowed down by its pictures)
        gap = max(self.pace if self.pace > 0 else SCENE_INTERVAL, 2.5 * self._cost)
        if self.enabled and (f.get("all") or now - self._scene_t >= gap):
            self._scene()
        if self.pace > 0:
            self.sleep(self.pace)

    # ------------------------------------------------------------------ helpers
    def _flush_line(self) -> None:
        if self._pending is not None:
            ev, self._pending = self._pending, None
            self.push("run", event="line", **ev)

    def _scene(self, force: bool = False) -> None:
        if not (self._dirty or force):
            return
        self._flush_line()
        t0 = self.clock()
        self._dirty = False
        try:
            data = snapshot(self.interp, self.state())
        except Exception as exc:                 # noqa: BLE001 -- a display aid
            data = {"empty": True, "error": str(exc)}
        self._scene_t = self.clock()
        self._cost = self._scene_t - t0
        # the event log keeps thousands of events: only the newest snapshot keeps its data (a browser that
        # reads an older one in the same batch draws the newest)
        if self._last_scene and self.replace is not None:
            self.replace(self._last_scene, data={"stale": True})
        self._last_scene = self.push("run", event="scene", data=data) or 0


def _head(line: str) -> str:
    """The canonical name of the command of a line ("" for comments and unknown commands)."""
    from ..prep.lexer import is_comment, split_head
    from ..prep.registry import lookup
    text = (line or "").strip()
    if not text or is_comment(text):
        return ""
    try:
        spec = lookup(split_head(text)[0])
    except Exception:                            # noqa: BLE001
        return ""
    return spec.name if spec is not None else ""
