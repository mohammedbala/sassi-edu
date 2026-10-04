"""Plot commands (manual 9.14; requirements 3.4.L, 5.7, 5.8; spec 06; spec 10 section 4).

Plot creation: SPECPLOT, THPLOT, LAYERPLOT, SOILPROPPLOT, MODELPLOT, NODEPLOT, CUTPLOT, BUBBLEPLOT,
VECTORPLOT, CONTOURPLOT, DEFORMPLOT (+ PROCFRAME, which fills the frame store the animations
read).  Plot settings: AXES, PLOTRANGE, PLOTTITLE, XTITLE, YTITLE, YTITLE2, LINENAME, MARKERS,
STIPPLE, ELECOLOR, ELENUM, GROUPNUM, NODENUM, NODESEL, SHOWDOF, SHOWMASS, SHRINK, WIREFRAME,
CNGVIEW, RSTVIEW, CNGCENTER, RSTCENTER, PAUSE, DEBUG, SHADEROPTIONS, COLOR, WINDOWSETTINGS.
Output: CAPTUREPLOT, CLOSEPLOT.  Tabs: ACTIVATEPLOT (SASSI-EDU extension, rule L17).

Every command works on the session :class:`sassi.plotting.state.PlotState` (``plot_state``).  The
state notifies its listeners (the GUI) of every new or changed plot.  In console / batch mode
(``PlotState.auto_render``) each new or changed plot is also rendered headless to
``<model>_plot<NN>_<kind>.png`` in the model directory and the path is printed (UI-03).

Conventions (spec 10 section 4.0): a setting command acts on the **active** plot when that plot is
capable, is ignored with a warning otherwise, and becomes a default for new plots when no plot is
open (2D settings: new SPECPLOT / THPLOT plots only).  Model-level display requests (NODESEL,
WINDOWSETTINGS hide / show / VOLUME fields) edit the model of the active 3D plot, else the active
model.  Toggles: ELENUM, GROUPNUM, NODENUM, SHOWMASS, SHRINK, STIPPLE, WIREFRAME -1 toggle / 0 off /
1 on; DEBUG 2 toggle; PAUSE -1 toggle / 0 start / 1 stop.  Commands that open a dialog when given
without arguments (CUTPLOT, SOILPROPPLOT, PROCFRAME, the animations, SHOWDOF, SHADEROPTIONS,
WINDOWSETTINGS, SPECPLOT/THPLOT "Line Selection") fail in batch / console mode (D-UI-08).
"""
from __future__ import annotations

import functools
from pathlib import Path
from typing import List, Optional

from ...plotting.lines import LineError
from ...plotting.state import (ALL_KINDS, ANIMATIONS, LINE_PALETTES, PALETTES, SHOWDOF_LABELS, THREE_D, FrameStore,
                               PlotError, cut_elements, frame_sequence, layer_table, model_scene, plot_state,
                               process_frames, soil_property_curves)
from ..registry import CommandError, command


def plot_command(name: str, **opts):
    """``@command`` that reports :class:`PlotError` / :class:`LineError` as command errors."""
    def deco(fn):
        @functools.wraps(fn)
        def wrapped(c):
            try:
                return fn(c)
            except (PlotError, LineError) as exc:
                raise CommandError(str(exc)) from None
        command(name, **opts)(wrapped)
        return fn
    return deco


def _state(c):
    return plot_state(c.interp)


def needs_dialog(c, dialog: str) -> None:
    """A command given without its arguments: the GUI opens ``dialog``; batch mode fails (D-UI-08)."""
    if _state(c).dialog(dialog, command=c.name):
        c.info(f"{c.name}: '{dialog}' dialog opened")
        return
    raise CommandError(f"requires arguments in batch mode (the '{dialog}' dialog is available in the GUI only, "
                       f"D-UI-08)")


def publish(c, plot, opened: bool) -> Optional[Path]:
    """Notify the listeners and, in console / batch mode, render the plot to its image file."""
    st = _state(c)
    st.emit("open" if opened else "update", plot.id, plot=plot.to_dict())
    if not st.auto_render:
        if opened:
            c.confirm(f"{plot.caption}: plot {plot.id} opened")
        return None
    path = st.image_path(plot, c.interp)
    try:
        from ...plotting.render_mpl import render_plot
        out, notes = render_plot(st, plot, c.interp, path)
    except Exception as exc:          # the image is a by-product: never fail the command for it
        c.warn(f"plot {plot.id} could not be rendered: {exc}")
        return None
    plot.image = str(out)
    for n in notes:
        c.warn(n)
    if opened:
        c.info(f"{plot.caption}: plot {plot.id} written to {out}")
    else:
        c.confirm(f"{plot.caption}: plot {plot.id} updated ({out})")
    return out


def refresh_line_plots(c, numbers) -> None:
    """Re-publish the open plots that show the given lines (after READSPEC, LINENAME ...)."""
    for p in _state(c).plots_showing(numbers):
        publish(c, p, opened=False)


def toggle(c, k: int, current: bool, toggle_value: int = -1) -> bool:
    """Toggle argument: ``toggle_value`` (default) toggles, 0 off, 1 on."""
    v = c.int(k, default=toggle_value)
    if v == toggle_value:
        return not current
    if v in (0, 1):
        return bool(v)
    raise CommandError(f"argument {k} must be {toggle_value} (toggle), 0 (off) or 1 (on)")


def line_list(c, k_from: int = 1, limit: int = 50) -> List[int]:
    """Line numbers from argument ``k_from``; -1 ends the list (spec 10 section 4.3)."""
    nums: List[int] = []
    for k in range(k_from, c.nargs + 1):
        if not c.given(k):
            continue
        v = c.int(k)
        if v == -1:
            break
        nums.append(v)
    if len(nums) > limit:
        c.warn(f"at most {limit} lines; lines after the {limit}th ignored")
        nums = nums[:limit]
    return list(dict.fromkeys(nums))


def _target(c, command: str):
    """``(plot or None, ok)``: the active capable plot, or None for the session defaults."""
    st = _state(c)
    plot, where = st.target(command)
    if where == "ignored":
        c.warn(f"not applicable to the active plot '{plot.caption}'; ignored")
        return None, False
    if where == "none":
        c.warn("no active plot; ignored")
        return None, False
    return plot, True


def _done(c, plot, what: str, family: str) -> None:
    if plot is not None:
        publish(c, plot, opened=False)
    else:
        c.confirm(f"{what}: no plot open, stored as default for new {family} plots")


def _model(c, number: Optional[int] = None):
    n = c.interp.active_model if number is None or number == -1 else number
    m = c.interp.models.get(n)
    if m is None:
        raise CommandError(f"model {n} is not in memory")
    return n, m


def _display_model(c):
    """Model whose display state (``ui_state``: hide requests, display volume, NODESEL) a command
    edits: the model of the active 3D plot -- a CUTPLOT may show a model that is not the active
    one, and its dialog acts on what it shows -- otherwise the active model."""
    p = _state(c).active_plot
    if p is not None and p.family in ("3d", "anim") and p.model in c.interp.models:
        return p.model, c.interp.models[p.model]
    return c.interp.active_model, c.model


def _refresh_model_plots(c, number: int) -> None:
    """Re-publish the open 3D plots of model ``number`` (model-level display state changed)."""
    for p in list(_state(c).plots.values()):
        if p.model == number and p.family in ("3d", "anim"):
            publish(c, p, opened=False)


# ======================================================================================
# 2D line plots
# ======================================================================================
def _lineplot(c, kind: str) -> None:
    if c.nargs == 0:
        needs_dialog(c, "Line Selection")
        return
    st = _state(c)
    nums = line_list(c)
    unknown = [n for n in nums if n not in st.lines]
    if unknown:
        c.warn("line(s) " + ", ".join(str(n) for n in unknown) + " not defined; ignored")
    known = [n for n in nums if n in st.lines]
    if not known:
        raise CommandError("no defined line to plot (load lines with READSPEC / READTH)")
    plot = st.open_plot(kind, {"lines": known})
    publish(c, plot, opened=True)


@plot_command("SPECPLOT")
def cmd_specplot(c):
    """SPECPLOT,<Line1>,...,<Line50>: spectrum plot of line objects (-1 ends the list)."""
    _lineplot(c, "SPECPLOT")


@plot_command("THPLOT")
def cmd_thplot(c):
    """THPLOT,<Line1>,...,<Line50>: time-history plot of line objects (-1 ends the list)."""
    _lineplot(c, "THPLOT")


# ======================================================================================
# Soil plots
# ======================================================================================
LAYER_COLUMNS = ("THICK", "WEIGHT", "VP", "VS", "PDAMP", "SDAMP")


@plot_command("LAYERPLOT")
def cmd_layerplot(c):
    """LAYERPLOT: soil layer column and property table of the active model (TOPL / L, SITE half-space)."""
    n, m = _model(c)
    layer_table(m)                      # validates (error when no layer is defined)
    show = {"thick": False, "weight": True, "vp": True, "vs": True, "pdamp": True, "sdamp": True}
    plot = _state(c).open_plot("LAYERPLOT", {"start": 1, "end": -1, "show": show}, model=n)
    publish(c, plot, opened=True)


@plot_command("SOILPROPPLOT")
def cmd_soilpropplot(c):
    """SOILPROPPLOT,<PropName>: G/Gmax and damping vs shear strain % of a DYNP property (case-sensitive)."""
    name = c.raw(1) if c.nargs <= 1 else c.rest.strip()      # DYNP labels may contain commas (L7)
    if not name:
        needs_dialog(c, "Select Dynamic Soil Property")
        return
    n, m = _model(c)
    soil_property_curves(m, name)        # validates the name
    plot = _state(c).open_plot("SOILPROPPLOT", {"name": name, "show": {"modulus": True, "damping": True}},
                               model=n)
    publish(c, plot, opened=True)


# ======================================================================================
# 3D model plots
# ======================================================================================
def _require_elements(m, n: int) -> None:
    if not any(g.elements for g in m.groups.values()):
        raise CommandError(f"model {n} has no elements to plot")


@plot_command("MODELPLOT", max_args=0)
def cmd_modelplot(c):
    """MODELPLOT: element plot of the active model (colours by ELECOLOR)."""
    n, m = _model(c)
    _require_elements(m, n)
    plot = _state(c).open_plot("MODELPLOT", {}, model=n)
    publish(c, plot, opened=True)


@plot_command("NODEPLOT", max_args=0)
def cmd_nodeplot(c):
    """NODEPLOT: node plot of the active model (element-connected nodes; interaction nodes red)."""
    n, m = _model(c)
    _require_elements(m, n)
    plot = _state(c).open_plot("NODEPLOT", {}, model=n)
    publish(c, plot, opened=True)


@plot_command("CUTPLOT", max_args=2)
def cmd_cutplot(c):
    """CUTPLOT,[Cut],[Model]: wireframe of a model with the elements of a cut filled (any model)."""
    if c.nargs == 0:
        needs_dialog(c, "Select Cut to Display")
        return
    cut = c.int(1, default=1)
    n, m = _model(c, c.int(2, default=-1))
    elems = cut_elements(c.interp, cut)
    if elems is None:
        raise CommandError(f"cut {cut} is not defined (CUTADD / CUTVOL / SLICE)")
    _require_elements(m, n)
    keys = set((g, e) for g, e in elems)
    have = sum(1 for g in m.groups.values() for e in g.elements if (g.id, e) in keys)
    if have < len(keys):
        c.warn(f"{len(keys) - have} element(s) of cut {cut} are not in model {n}")
    plot = _state(c).open_plot("CUTPLOT", {"cut": cut, "elements": len(keys)}, model=n)
    publish(c, plot, opened=True)


# ======================================================================================
# Frame processing and animations
# ======================================================================================
@plot_command("PROCFRAME")
def cmd_procframe(c):
    """PROCFRAME,[AniFile],[BufferDIR],[Data],[anitype]: frame list -> frame store + SASSIani.xml."""
    if not c.given(1) or not c.given(2):
        needs_dialog(c, "Parse Frame Data")
        return
    st = _state(c)
    listfile = c.input_path(c.str(1))
    bufdir = c.output_path(c.str(2))
    desc = c.raw(3)
    anitype = 0
    if c.nargs >= 4:
        parts = c.rest.split(",")
        try:
            anitype = c.int(c.nargs)
            if c.nargs > 4 and len(parts) >= c.nargs:      # unquoted description containing commas
                desc = ",".join(parts[2:c.nargs - 1]).strip()
        except CommandError:                               # last token not a number: all description
            desc = ",".join(parts[2:]).strip() if len(parts) >= c.nargs else desc
    if anitype not in (0, 1, 2, 3):
        c.warn(f"anitype {anitype}: 0 bubble, 1 vector, 2 contour, 3 time history expected (sorting tag only)")

    def progress(frac: float, item: str) -> None:
        st.emit("progress", None, fraction=frac, item=item, command="PROCFRAME")

    idx = process_frames(listfile, bufdir, desc, anitype, progress=progress)
    for note in idx.get("notes", []):
        c.warn(note)
    entry ={"description": desc, "directory": str(Path(bufdir).resolve()), "type": anitype,
             "frames": idx["nframes"], "start": 1, "end": idx["nframes"], "stride": 1, "scale": 1.0,
             "cmin": idx["min"][0], "cmax": idx["max"][0], "listfile": str(listfile)}
    try:
        db = st.ani_db()
        db.upsert(entry)
        where = f"; {db.path.name} updated"
    except OSError as exc:
        c.warn(f"the animation database could not be written: {exc}")
        where = ""
    st.emit("animations", None, directory=entry["directory"], nframes=idx["nframes"])
    c.info(f"PROCFRAME: {idx['nframes']} frames ({idx['data_columns']} data columns, {idx['layout']}) stored in "
           f"{bufdir}{where}")


_ANIM_ARGS = {"BUBBLEPLOT": ("BufferDir", "MnF", "MxF", "ST", "MnR", "MxR", "Col"),
              "CONTOURPLOT": ("BufferDir", "MnF", "MxF", "ST", "MnR", "MxR", "Col"),
              "VECTORPLOT": ("BufferDir", "MnF", "MxF", "ST", "Scale"),
              "DEFORMPLOT": ("BufferDir", "MnF", "MxF", "ST", "Scale")}


def _num(entry, key, cast, default):
    try:
        return cast(entry[key]) if entry and entry.get(key) not in (None, "") else default
    except (TypeError, ValueError):
        return default


def _animation(c, kind: str) -> None:
    if c.nargs == 0:
        needs_dialog(c, "Load Frame Data")
        return
    st = _state(c)
    names = _ANIM_ARGS[kind]
    if not c.given(1):
        raise CommandError("<BufferDir> (processed frame directory) is required")
    bdir = c.interp.resolve_path(c.str(1), must_exist=True)
    store = FrameStore(bdir)
    nf = store.nframes
    try:
        entry = st.ani_db().find(bdir)
    except OSError:
        entry = None
    missing = [names[k - 1] for k in range(2, len(names) + 1) if not c.given(k)]
    if missing:
        c.warn("partial arguments: " + ", ".join(missing) + " taken from the animation defaults "
               "(SASSIani.xml / frame store; manual: partial entry may give inconsistent results)")
    start = c.int(2, default=_num(entry, "start", int, 1))
    end = c.int(3, default=min(_num(entry, "end", int, nf), nf))
    stride = c.int(4, default=_num(entry, "stride", int, 1))
    frames = frame_sequence(start, end, stride, nf)
    params = {"buffer_dir": str(bdir), "frames": frames, "start": start, "end": end, "stride": stride,
              "current": frames[0], "nframes": nf, "layout": store.layout, "points": st.shader.points}
    if kind in ("BUBBLEPLOT", "CONTOURPLOT"):
        col = c.int(7, default=1)
        if not 1 <= col <= store.data_columns:
            raise CommandError(f"<Col> {col} outside 1..{store.data_columns} (data columns of the frames; the node "
                               f"id is not a data column)")
        lo, hi = store.column_range(frames, col)
        vmin = c.float(5, default=_num(entry, "cmin", float, lo) if col == 1 else lo)
        vmax = c.float(6, default=_num(entry, "cmax", float, hi) if col == 1 else hi)
        if vmin > vmax:
            raise CommandError(f"<MnR> {vmin:g} > <MxR> {vmax:g}")
        if vmin == vmax:
            vmax = vmin + (abs(vmin) if vmin else 1.0)
            c.warn(f"constant data: colour-map range widened to {vmin:g}..{vmax:g}")
        params.update(col=col, vmin=float(vmin), vmax=float(vmax))
    else:
        scale = c.float(5, default=_num(entry, "scale", float, st.shader.scale))
        params.update(scale=float(scale))
    n, m = _model(c)
    sc = model_scene(m, apply_hide=False)          # connectivity check: hide requests / volume do not count
    if sc["n_elements"] == 0:
        raise CommandError(f"model {n} has no elements: activate the model whose results the frames hold")
    nodes = store.nodes(frames[0])
    known = set(int(v) for v in sc["node_id"])
    absent = sum(1 for v in nodes if int(v) not in known)
    if absent:
        c.warn(f"{absent} of {len(nodes)} frame nodes are not element-connected nodes of model {n} (not drawn)")
    plot = st.open_plot(kind, params, model=n)
    publish(c, plot, opened=True)


@plot_command("BUBBLEPLOT", max_args=7)
def cmd_bubbleplot(c):
    """BUBBLEPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>: bubble animation of nodal values."""
    _animation(c, "BUBBLEPLOT")


@plot_command("CONTOURPLOT", max_args=7)
def cmd_contourplot(c):
    """CONTOURPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>: contour animation on element faces."""
    _animation(c, "CONTOURPLOT")


@plot_command("VECTORPLOT", max_args=5)
def cmd_vectorplot(c):
    """VECTORPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>: animated complex-TF vectors (X red, Y green, Z blue)."""
    _animation(c, "VECTORPLOT")


@plot_command("DEFORMPLOT", max_args=5)
def cmd_deformplot(c):
    """DEFORMPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>: animated deformed shape x' = x + Scale u."""
    _animation(c, "DEFORMPLOT")


@plot_command("PAUSE", max_args=1)
def cmd_pause(c):
    """PAUSE,[pz]: animation -1 toggle start/stop (default), 0 start, 1 stop."""
    plot, ok = _target(c, "PAUSE")
    if not ok:
        return
    v = c.int(1, default=-1)
    if v == -1:
        plot.view.paused = not plot.view.paused
    elif v in (0, 1):
        plot.view.paused = bool(v)
    else:
        raise CommandError("argument 1 must be -1 (toggle), 0 (start) or 1 (stop)")
    _state(c).emit("update", plot.id, plot=plot.to_dict())
    c.confirm(f"animation {'stopped' if plot.view.paused else 'running'}")


# ======================================================================================
# 2D settings
# ======================================================================================
@plot_command("AXES", max_args=6)
def cmd_axes(c):
    """AXES,<MaxTickX>,<MaxTickY>,<MinTickX>,<MinTickY>,<LogX>,<LogY>: grids and log axes (0/1; blank unchanged)."""
    plot, ok = _target(c, "AXES")
    if not ok:
        return
    s = plot.settings if plot is not None else _state(c).defaults2d
    for k, attr in enumerate(("major_x", "major_y", "minor_x", "minor_y", "log_x", "log_y"), start=1):
        if c.given(k):
            v = c.int(k)
            if v not in (0, 1):
                raise CommandError(f"argument {k} must be 0 or 1")
            setattr(s, attr, bool(v))
    if s.log_x and s.xmin is not None and s.xmin <= 0:
        c.warn("log X axis with a non-positive PLOTRANGE minimum: the X range reset to the data extent")
        s.xmin = s.xmax = None
    if s.log_y and s.ymin is not None and s.ymin <= 0:
        c.warn("log Y axis with a non-positive PLOTRANGE minimum: the Y range reset to the data extent")
        s.ymin = s.ymax = None
    _done(c, plot, "AXES", "spectrum / time-history")


@plot_command("PLOTRANGE", max_args=4)
def cmd_plotrange(c):
    """PLOTRANGE,<Xmin>,<Xmax>,<Ymin>,<Ymax>: extent of the active 2D plot (blank = data extent)."""
    plot, ok = _target(c, "PLOTRANGE")
    if not ok:
        return
    s = plot.settings if plot is not None else _state(c).defaults2d
    x0, x1, y0, y1 = (c.float(k) for k in range(1, 5))
    for a, b, ax, log in ((x0, x1, "X", s.log_x), (y0, y1, "Y", s.log_y)):
        if a is not None and b is not None and a >= b:
            raise CommandError(f"{ax} minimum {a:g} must be smaller than the maximum {b:g}")
        if log and a is not None and a <= 0:
            raise CommandError(f"{ax} axis is logarithmic: the minimum must be > 0")
    s.xmin, s.xmax, s.ymin, s.ymax = x0, x1, y0, y1
    _done(c, plot, "PLOTRANGE", "spectrum / time-history")


@plot_command("PLOTTITLE")
def cmd_plottitle(c):
    """PLOTTITLE,<Title>: title of the active plot (2D and 3D; ignored by the soil layer plot)."""
    plot, ok = _target(c, "PLOTTITLE")
    if not ok:
        return
    if plot is None:
        c.warn("no active plot; ignored")
        return
    plot.title = c.text(1)
    publish(c, plot, opened=False)


def _axis_title(c, command: str, attr: str) -> None:
    plot, ok = _target(c, command)
    if not ok:
        return
    s = plot.settings if plot is not None else _state(c).defaults2d
    setattr(s, attr, c.text(1))
    _done(c, plot, command, "spectrum / time-history")


@plot_command("XTITLE")
def cmd_xtitle(c):
    """XTITLE,<label>: X-axis title of the active 2D plot."""
    _axis_title(c, "XTITLE", "xtitle")


@plot_command("YTITLE")
def cmd_ytitle(c):
    """YTITLE,<label>: left Y-axis title of the active 2D plot."""
    _axis_title(c, "YTITLE", "ytitle")


@plot_command("YTITLE2")
def cmd_ytitle2(c):
    """YTITLE2,<label>: right Y-axis title (soil property plot: shear modulus axis)."""
    _axis_title(c, "YTITLE2", "ytitle2")


@plot_command("STIPPLE", max_args=1, tier="P1")
def cmd_stipple(c):
    """STIPPLE,<switch>: dash patterns on the lines of the active plot (-1 toggle, 0 off, 1 on)."""
    plot, ok = _target(c, "STIPPLE")
    if not ok:
        return
    s = plot.settings if plot is not None else _state(c).defaults2d
    s.stipple = toggle(c, 1, s.stipple)
    _done(c, plot, "STIPPLE", "spectrum / time-history")


@plot_command("LINENAME", max_args=2)
def cmd_linename(c):
    """LINENAME,<Num>,<Name>: rename a line object (global: every graph shows the new name)."""
    num = c.int(1, required=True)
    st = _state(c)
    if num not in st.lines:
        raise CommandError(f"line {num} not defined")
    st.lines[num].name = c.text(2)
    st.emit("lines", None, numbers=[num])
    c.confirm(f"line {num}: name '{st.lines[num].name}'")
    refresh_line_plots(c, [num])


@plot_command("MARKERS")
def cmd_markers(c):
    """MARKERS,<Mark>,<Ln1>,...,<Ln50>: data-point markers off (0) / on (1) for lines (global)."""
    mark = c.int(1, required=True)
    if mark not in (0, 1):
        raise CommandError("<Mark> must be 0 (off) or 1 (on)")
    st = _state(c)
    nums = line_list(c, 2)
    unknown = [n for n in nums if n not in st.lines]
    if unknown:
        c.warn("line(s) " + ", ".join(str(n) for n in unknown) + " not defined; ignored")
    known = [n for n in nums if n in st.lines]
    for n in known:
        st.lines[n].markers = bool(mark)
    if known:
        st.emit("lines", None, numbers=known)
        c.confirm(f"markers {'on' if mark else 'off'} for line(s) " + ", ".join(str(n) for n in known))
        refresh_line_plots(c, known)


# ======================================================================================
# Output
# ======================================================================================
@plot_command("CAPTUREPLOT")
def cmd_captureplot(c):
    """CAPTUREPLOT,<FileName>: save the active plot (PNG if '.png' is in the name, else BMP)."""
    st = _state(c)
    plot = st.active_plot
    if plot is None:
        raise CommandError("no active plot to capture")
    name = c.rest.strip() if c.nargs > 1 else c.str(1)
    if not name:
        raise CommandError("<FileName> is required")
    path = c.output_path(name)
    from ...plotting.render_mpl import image_format, render_plot
    fmt = image_format(name)            # D-UI-10: decided by the <FileName> text, never by the directory
    try:
        out, notes = render_plot(st, plot, c.interp, path, fmt)
    except OSError as exc:
        raise CommandError(f"cannot write {path}: {exc}") from None
    except (PlotError, LineError, ValueError, KeyError) as exc:
        raise CommandError(f"the plot could not be rendered: {exc}") from None
    for n in notes:
        c.warn(n)
    st.emit("capture", plot.id, path=str(out), format=fmt)
    c.info(f"CAPTUREPLOT: {plot.caption} saved to {out} ({fmt.upper()})")


@plot_command("CLOSEPLOT", max_args=0)
def cmd_closeplot(c):
    """CLOSEPLOT: close the active plot (the most recent remaining plot becomes active)."""
    st = _state(c)
    p = st.close_plot()
    if p is None:
        c.warn("no plot is open")
        return
    st.emit("close", p.id, active=st.active)
    c.confirm(f"{p.caption} (plot {p.id}) closed" + (f"; plot {st.active} active" if st.active else ""))


@plot_command("ACTIVATEPLOT", max_args=1, tier="P1", cls="ui")
def cmd_activateplot(c):
    """ACTIVATEPLOT,<id>: make open plot <id> the active plot (SASSI-EDU extension, rule L17).

    Most plot commands act on the active plot (the tab in front).  Selecting a tab in the GUI is
    therefore a state change that a replayed session must reproduce: the GUI submits this
    command when the user brings a tab to the front.  Plot ids are the numbers printed when the
    plots open (deterministic in a replay).
    """
    st = _state(c)
    pid = c.int(1, required=True, what="id")
    if pid not in st.plots:
        opened = ", ".join(str(i) for i in st.order) or "none"
        raise CommandError(f"plot {pid} is not open (open plots: {opened})")
    p = st.activate(pid)
    st.emit("activate", p.id, plot=p.to_dict())
    c.confirm(f"{p.caption}: plot {p.id} active")


# ======================================================================================
# 3D display settings
# ======================================================================================
def _view_target(c, command: str):
    plot, ok = _target(c, command)
    if not ok:
        return None, None
    return plot, (plot.view if plot is not None else _state(c).defaults3d)


@plot_command("ELECOLOR", max_args=1)
def cmd_elecolor(c):
    """ELECOLOR,<val>: element colours by 1 group, 2 material, 3 property (ElemPalette, 128 colours)."""
    plot, v = _view_target(c, "ELECOLOR")
    if v is None:
        return
    val = c.int(1, required=True)
    if val not in (1, 2, 3):
        raise CommandError("<val> must be 1 (group), 2 (material) or 3 (property)")
    v.color_by = val
    _done(c, plot, "ELECOLOR", "3D")


def _label_toggle(c, command: str, attr: str, other: Optional[str] = None) -> None:
    plot, v = _view_target(c, command)
    if v is None:
        return
    on = toggle(c, 1, getattr(v, attr))
    setattr(v, attr, on)
    if on and other:
        setattr(v, other, False)          # ELENUM and GROUPNUM are mutually exclusive
    _done(c, plot, command, "3D")


@plot_command("ELENUM", max_args=1)
def cmd_elenum(c):
    """ELENUM,[opt]: element number labels (-1 toggle, 0 off, 1 on); turns GROUPNUM off."""
    _label_toggle(c, "ELENUM", "elem_labels", "group_labels")


@plot_command("GROUPNUM", max_args=1)
def cmd_groupnum(c):
    """GROUPNUM,[opt]: group number labels (-1 toggle, 0 off, 1 on); turns ELENUM off."""
    _label_toggle(c, "GROUPNUM", "group_labels", "elem_labels")


@plot_command("NODENUM", max_args=1)
def cmd_nodenum(c):
    """NODENUM,[opt]: node number labels (-1 toggle, 0 off, 1 on)."""
    _label_toggle(c, "NODENUM", "node_labels")


@plot_command("SHOWMASS", max_args=1)
def cmd_showmass(c):
    """SHOWMASS,[opt]: lumped-mass markers (-1 toggle, 0 off, 1 on); red X, green Y, blue Z."""
    _label_toggle(c, "SHOWMASS", "show_mass")


@plot_command("SHRINK", max_args=1)
def cmd_shrink(c):
    """SHRINK,[switch]: shrink the elements of the element plot (1 on, 0 off, -1 toggle)."""
    _label_toggle(c, "SHRINK", "shrink")


@plot_command("WIREFRAME", max_args=1)
def cmd_wireframe(c):
    """WIREFRAME,<Switch>: wireframe mode of the element plot (0 off, 1 on, -1 toggle)."""
    _label_toggle(c, "WIREFRAME", "wireframe")


@plot_command("DEBUG", max_args=1, tier="P1")
def cmd_debug(c):
    """DEBUG,[switch]: view values / animation information on 3D plots (0 off, 1 on, 2 toggle)."""
    plot, v = _view_target(c, "DEBUG")
    if v is None:
        return
    v.debug = toggle(c, 1, v.debug, toggle_value=2)
    _done(c, plot, "DEBUG", "3D")


@plot_command("SHOWDOF", max_args=6)
def cmd_showdof(c):
    """SHOWDOF,[label1],...,[label6]: mark nodes with fixed DOFs (X Y Z XX YY ZZ DISP ROT ALL; NONE off)."""
    labels = [c.word(k) for k in range(1, c.nargs + 1) if c.given(k)]
    if not labels:
        needs_dialog(c, "Boundary Conditions")
        return
    dofs = set()
    for lab in labels:
        if lab not in SHOWDOF_LABELS:
            raise CommandError(f"unknown label {lab} (X, Y, Z, XX, YY, ZZ, DISP, ROT, ALL, NONE)")
        dofs.update(SHOWDOF_LABELS[lab])
    plot, v = _view_target(c, "SHOWDOF")
    if v is None:
        return
    v.show_dof = sorted(dofs)
    _done(c, plot, "SHOWDOF", "3D")


@plot_command("CNGVIEW", max_args=6)
def cmd_cngview(c):
    """CNGVIEW,<rX>,<rY>,<rZ>,<px>,<py>,<zoom>: 3D view (degrees, pan in model units, zoom 1 = fit)."""
    plot, v = _view_target(c, "CNGVIEW")
    if v is None:
        return
    vals = [c.float(k) for k in range(1, 7)]
    if vals[5] is not None and not vals[5] > 0:
        raise CommandError("<zoom> must be > 0")
    for attr, val in zip(("rx", "ry", "rz", "px", "py", "zoom"), vals):
        if val is not None:
            setattr(v, attr, float(val))
    _done(c, plot, "CNGVIEW", "3D")


@plot_command("RSTVIEW", max_args=0)
def cmd_rstview(c):
    """RSTVIEW: reset the view (default isometric view, zoom and pan fitting the model)."""
    plot, v = _view_target(c, "RSTVIEW")
    if v is None:
        return
    v.reset_view()
    _done(c, plot, "RSTVIEW", "3D")


@plot_command("CNGCENTER", max_args=3)
def cmd_cngcenter(c):
    """CNGCENTER,<X>,<Y>,<Z>: centre of rotation of the 3D plot (global coordinates)."""
    plot, v = _view_target(c, "CNGCENTER")
    if v is None:
        return
    v.center = [c.float(k, required=True) for k in (1, 2, 3)]
    _done(c, plot, "CNGCENTER", "3D")


@plot_command("RSTCENTER", max_args=0)
def cmd_rstcenter(c):
    """RSTCENTER: rotation centre = centre of the bounding box of the element-connected nodes."""
    plot, v = _view_target(c, "RSTCENTER")
    if v is None:
        return
    v.center = None
    _done(c, plot, "RSTCENTER", "3D")


@plot_command("NODESEL", max_args=20)
def cmd_nodesel(c):
    """NODESEL,<N1>,...,<N20>: toggle the selection of nodes (model level: all current and future plots).

    The model is the one of the active 3D plot, else the active model (:func:`_display_model`).
    """
    num, m = _display_model(c)
    ids = [c.int(k) for k in range(1, c.nargs + 1) if c.given(k)]
    if not ids:
        raise CommandError("node numbers required")
    sel = list(m.ui_state.get("nodesel", []))
    added, removed, unknown = [], [], []
    for n in ids:
        if n not in m.nodes:
            unknown.append(n)
        elif n in sel:
            sel.remove(n)
            removed.append(n)
        else:
            sel.append(n)
            added.append(n)
    if unknown:
        c.warn("node(s) " + ", ".join(str(n) for n in unknown) + " not defined; ignored")
    m.ui_state["nodesel"] = sel
    st = _state(c)
    st.emit("nodesel", None, model=num, selected=list(sel))
    c.confirm(f"model {num}: {len(sel)} node(s) selected" + (f"; added {added}" if added else "") +
              (f"; deselected {removed}" if removed else ""))
    _refresh_model_plots(c, num)


@plot_command("SHADEROPTIONS", max_args=4, tier="P1")
def cmd_shaderoptions(c):
    """SHADEROPTIONS,[points],[linew],[shrink],[scale]: node size, outline, shrink fraction, scale factor."""
    if c.nargs == 0:
        needs_dialog(c, "Shader Options")
        return
    st = _state(c)
    sh = st.shader
    for k, attr in enumerate(("points", "linew", "shrink", "scale"), start=1):
        val = c.float(k)
        if val is None:
            continue
        if attr in ("linew", "shrink") and not 0.0 <= val < 0.5:
            raise CommandError(f"<{attr}> is a fraction of the element (0 <= value < 0.5)")
        if attr in ("points", "scale") and not val > 0:
            raise CommandError(f"<{attr}> must be > 0")
        setattr(sh, attr, float(val))
    plot = st.active_plot
    if plot is not None and plot.kind in ANIMATIONS:
        if c.given(4) and "scale" in plot.params:
            plot.params["scale"] = sh.scale
        if c.given(1):
            plot.params["points"] = sh.points
    if plot is not None and plot.family in ("3d", "anim"):
        publish(c, plot, opened=False)
    else:
        c.confirm("shader options set")


@plot_command("COLOR", max_args=5, tier="P1")
def cmd_color(c):
    """COLOR,<Palette>,<Num>,<R>,<G>,<B>: set colour Num (1-based) of a palette (0-255 values).

    SpecLines / THLines have one entry per line number, so any ``<Num> >= 1`` (the line number)
    is accepted; the other palettes have a fixed size.
    """
    name = c.word(1)
    if name not in PALETTES:
        raise CommandError(f"unknown palette {c.raw(1)} (" + ", ".join(sorted(PALETTES)) + ")")
    num = c.int(2, required=True)
    if name in LINE_PALETTES:
        if num < 1:
            raise CommandError(f"<Num> must be a line number >= 1 for palette {name}")
    elif not 1 <= num <= len(PALETTES[name]):
        raise CommandError(f"<Num> must be 1..{len(PALETTES[name])} for palette {name}")
    rgb = [c.int(k, required=True) for k in (3, 4, 5)]
    if any(not 0 <= v <= 255 for v in rgb):
        raise CommandError("<R>, <G>, <B> must be 0..255")
    st = _state(c)
    st.palettes.setdefault(name, {})[num] = "#{:02x}{:02x}{:02x}".format(*rgb)
    c.confirm(f"{name} colour {num} = {st.palettes[name][num]}")
    if st.active_plot is not None:
        publish(c, st.active_plot, opened=False)


# ======================================================================================
# WINDOWSETTINGS (dialog; SASSI-EDU field form for scripts and the GUI, L17)
# ======================================================================================
_WS_HELP = ("TITLE,<text> | START,<n> | END,<n> | SHOW,<column>,<0|1> | SCALE,<s> | RANGE,<min>,<max> | "
            "DIRECTION,<X|Y|Z|ALL> | UNDEFORMED,<0|1> | FRAME,<k> | FRAMEPAUSE,<ms> | "
            "VOLUME,[xmin],[xmax],[ymin],[ymax],[zmin],[zmax] | HIDEGROUP,<g> | SHOWGROUP,<g> | "
            "HIDEELEM,<g>,<ids> | SHOWELEM,<g>,<ids> | HIDENODE,<ids> | SHOWNODE,<ids> | SHOWALL")


def _need(plot, kinds, field_name):
    if plot is None or plot.kind not in kinds:
        raise CommandError(f"field {field_name} does not apply to the active plot"
                           + (f" '{plot.caption}'" if plot else " (no plot open)"))


@plot_command("WINDOWSETTINGS")
def cmd_windowsettings(c):
    """WINDOWSETTINGS,[<field>,<values>]: settings window of the active plot.

    Without arguments the GUI opens the dialog of the active plot (Graph Plot Options, Window
    Options, Soil Layer / Soil Properties Window Settings); in batch mode that is an error
    (D-UI-08).  The field form (SASSI-EDU extension suggested by spec 10 section 4.4.3) sets one
    dialog field so every dialog action has command text (L17).
    """
    st = _state(c)
    if c.nargs == 0:
        needs_dialog(c, "Window Settings")
        return
    field_name = c.word(1)
    plot = st.active_plot
    mnum, m = _display_model(c)          # display requests: the model shown by the active 3D plot
    ui = m.ui_state
    if field_name == "TITLE":
        _need(plot, ALL_KINDS - {"LAYERPLOT"}, field_name)
        plot.title = c.rest.split(",", 1)[1].strip() if "," in c.rest else ""
    elif field_name in ("START", "END"):
        _need(plot, {"LAYERPLOT"}, field_name)
        v = c.int(2, required=True)
        new = dict(plot.params, **{"start" if field_name == "START" else "end": v})
        layer_table(c.interp.models[plot.model], new["start"], new["end"])      # validate first
        plot.params.update(start=new["start"], end=new["end"])
    elif field_name == "SHOW":
        _need(plot, {"LAYERPLOT", "SOILPROPPLOT"}, field_name)
        col, on = c.word(2), c.int(3, required=True)
        keys = {"LAYERPLOT": {k: k.lower() for k in LAYER_COLUMNS},
                "SOILPROPPLOT": {"MODULUS": "modulus", "DAMPING": "damping"}}[plot.kind]
        if col not in keys:
            raise CommandError(f"SHOW column must be one of {', '.join(keys)}")
        plot.params.setdefault("show", {})[keys[col]] = bool(on)
    elif field_name == "SCALE":
        _need(plot, {"VECTORPLOT", "DEFORMPLOT"}, field_name)
        s = c.float(2, required=True)
        if not s > 0:
            raise CommandError("scale factor must be > 0")
        plot.params["scale"] = float(s)
    elif field_name == "RANGE":
        _need(plot, {"BUBBLEPLOT", "CONTOURPLOT"}, field_name)
        a, b = c.float(2, required=True), c.float(3, required=True)
        if not a < b:
            raise CommandError("colour-map minimum must be smaller than the maximum")
        plot.params["vmin"], plot.params["vmax"] = float(a), float(b)
    elif field_name == "DIRECTION":
        _need(plot, {"VECTORPLOT"}, field_name)
        d = c.word(2)
        if d not in ("X", "Y", "Z", "ALL"):
            raise CommandError("DIRECTION must be X, Y, Z or ALL")
        plot.view.direction = d
    elif field_name == "UNDEFORMED":
        _need(plot, {"DEFORMPLOT"}, field_name)
        plot.view.show_undeformed = bool(c.int(2, required=True))
    elif field_name == "FRAME":
        _need(plot, ANIMATIONS, field_name)
        k = c.int(2, required=True)
        if k not in plot.params["frames"]:
            raise CommandError(f"frame {k} is not in the displayed sequence {plot.params['start']}.."
                               f"{plot.params['end']} step {plot.params['stride']}")
        plot.params["current"] = k
    elif field_name == "FRAMEPAUSE":
        _need(plot, THREE_D, field_name)
        ms = c.int(2, required=True)
        if ms < 1:
            raise CommandError("frame pause must be >= 1 ms")
        plot.params["frame_pause"] = ms
    elif field_name == "VOLUME":
        vals = [c.float(k) for k in range(2, 8)]
        if all(v is None for v in vals):
            ui.pop("display_volume", None)
        else:
            sc = model_scene(m, apply_hide=False)
            bb = sc["bbox"]
            vol = [bb[i] if vals[i] is None else float(vals[i]) for i in range(6)]
            if vol[0] > vol[1] or vol[2] > vol[3] or vol[4] > vol[5]:
                raise CommandError("display volume: each minimum must not exceed its maximum")
            ui["display_volume"] = vol
    elif field_name in ("HIDEGROUP", "SHOWGROUP"):
        g = c.int(2, required=True)
        if g not in m.groups:
            raise CommandError(f"group {g} not defined")
        hg = set(ui.get("hide_groups", []))
        (hg.add if field_name == "HIDEGROUP" else hg.discard)(g)
        ui["hide_groups"] = sorted(hg)
    elif field_name in ("HIDEELEM", "SHOWELEM"):
        g = c.int(2, required=True)
        if g not in m.groups:
            raise CommandError(f"group {g} not defined")
        ids = c.id_list(3)
        he = set(tuple(x) for x in ui.get("hide_elements", []))
        for e in ids:
            if e in m.groups[g].elements:
                (he.add if field_name == "HIDEELEM" else he.discard)((g, e))
        ui["hide_elements"] = [list(x) for x in sorted(he)]
    elif field_name in ("HIDENODE", "SHOWNODE"):
        ids = c.id_list(2)
        hn = set(ui.get("hide_nodes", []))
        for n in ids:
            if n in m.nodes:
                (hn.add if field_name == "HIDENODE" else hn.discard)(n)
        ui["hide_nodes"] = sorted(hn)
    elif field_name == "SHOWALL":
        for k in ("hide_groups", "hide_elements", "hide_nodes", "display_volume"):
            ui.pop(k, None)
    else:
        raise CommandError(f"unknown field {field_name} ({_WS_HELP})")
    if field_name in ("VOLUME", "HIDEGROUP", "SHOWGROUP", "HIDEELEM", "SHOWELEM", "HIDENODE", "SHOWNODE",
                      "SHOWALL"):
        c.confirm(f"display requests of model {mnum} updated")
        _refresh_model_plots(c, mnum)
        return
    publish(c, plot, opened=False)

