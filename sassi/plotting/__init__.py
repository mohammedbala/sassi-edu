"""SASSI-EDU plotting: line objects, line mathematics, plot state and headless rendering.

Requirements 3.4.L, 4.13, 5.7, 5.8; decisions D-LIN-01..04, D-UI-06..16, D-FIL-02/03; spec 06
(plots, toolbars) and spec 10 section 4 (exact command semantics).  The commands themselves are
in ``sassi/prep/commands/plotting.py`` (plots, settings, animations), ``linemath.py`` (line I/O
and maths) and ``frames.py`` (CRITFREQ, FRAMECOMBIN, FRAMESEL, MODFRAMES).

Package layout
--------------
``lines.py``       :class:`~.lines.Line`, union-grid resampling, ADDITION ... SRSS, BROADEN
                   (envelope -> +-b peak broadening -> peak bridging, D-LIN-01), spectrum / history
                   text files, CRITFREQ and FRAMESEL algorithms.  Pure numpy.
``state.py``       :class:`~.state.PlotState` (session state), plot requests, model scenes, soil
                   plot data, frame files, frame store, ``SASSIani.xml``, animation maths.  No
                   matplotlib import.
``render_mpl.py``  headless matplotlib (Agg) figures drawn from the plot data; PNG / BMP output.

Contract for the GUI (and any other front end)
----------------------------------------------
**One state per interpreter.**  ``state = sassi.plotting.state.plot_state(interp)`` returns the
:class:`PlotState` stored in ``interp.session['plot']`` (created on first use).  Every plotting
command of the interpreter edits it; the GUI never edits it directly but submits command text
(rule L17), e.g. ``SPECPLOT,1,2``, ``ELECOLOR,2``, ``WINDOWSETTINGS,HIDEGROUP,3``.

**Listeners.**  ``state.subscribe(fn)``; ``fn(event)`` receives a :class:`~.state.PlotEvent`
(``event.kind``, ``event.plot_id``, ``event.data``; ``event.to_dict()`` is JSON):

=============  ==============================================================================
``open``       a new plot request (``data['plot']`` = ``Plot.to_dict()``); it is the active plot
``update``     settings, view or data of an open plot changed (``data['plot']``)
``close``      CLOSEPLOT (``data['active']`` = new active plot id or None)
``activate``   ACTIVATEPLOT,<id> made plot ``event.plot_id`` the active one (``data['plot']``)
``lines``      line objects changed (``data['numbers']``): READSPEC, READTH, line maths,
               LINENAME, MARKERS
``capture``    CAPTUREPLOT wrote ``data['path']`` (``data['format']`` 'png' / 'bmp': PNG when the
               <FileName> text contains '.png', D-UI-10; the directory never decides)
``progress``   PROCFRAME progress (``data['fraction']`` 0..1, ``data['item']``)
``animations`` a frame store / ``SASSIani.xml`` entry was written (``data['directory']``)
``nodesel``    NODESEL changed the selection of ``data['model']`` (``data['selected']``)
``dialog``     a command without arguments asked for ``data['dialog']`` (see below)
=============  ==============================================================================

**Tabs.**  Setting commands act on the active plot, so bringing a tab to the front is a state
change: the GUI submits ``ACTIVATEPLOT,<id>`` (SASSI-EDU extension; L17 replay) and draws the
``activate`` event.  Plot ids are the numbers printed when the plots open.

A GUI that draws the plots itself sets ``state.auto_render = False`` (otherwise every new or
changed plot is also rendered to ``<model name or 'sassi'>_plot<NN>_<kind>.png`` in the model
directory -- MDL path, else the working directory -- and the path is printed; this is the
console / batch behaviour, UI-03).  A GUI that can show the manual's dialogs sets
``state.dialog_handler = fn(name, context) -> bool``; commands given without their arguments
(SPECPLOT/THPLOT "Line Selection", CUTPLOT "Select Cut to Display", SOILPROPPLOT "Select Dynamic
Soil Property", PROCFRAME "Parse Frame Data", animations "Load Frame Data", SHOWDOF "Boundary
Conditions", SHADEROPTIONS "Shader Options", WINDOWSETTINGS "Window Settings") call it and emit a
``dialog`` event; without a handler they fail (D-UI-08).

**Pure data.**  ``state.to_dict()`` is a JSON snapshot: ``lines`` {number: Line.to_dict()},
``plots`` [Plot.to_dict() in tab order], ``active``, ``defaults2d``, ``defaults3d``, ``shader``,
``palettes`` (COLOR overrides), ``auto_render``.  ``state.plot_data(plot_id, interp)`` returns
everything needed to draw one plot (JSON types; ``raw=True`` keeps numpy arrays):

* common: ``id kind family caption title model params settings view shader``; ``family`` is
  ``2d`` (SPECPLOT, THPLOT), ``layer``, ``soilprop``, ``3d`` (MODELPLOT, NODEPLOT, CUTPLOT) or
  ``anim`` (BUBBLEPLOT, VECTORPLOT, CONTOURPLOT, DEFORMPLOT);
* ``2d``: ``lines`` [{number, name, x, y, markers, color, dash}] in request order (colour by line
  number from the SpecLines / THLines palette -- one entry per line number, never cyclic, so the
  up to 50 lines of a plot have distinct colours; COLOR,SpecLines,<line>,R,G,B overrides one;
  dash patterns when STIPPLE is on), ``extent`` [xmin, xmax, ymin, ymax] (PLOTRANGE or the data
  extent);
* ``layer``: ``table`` {layers [{index, layer, thick, weight, vp, vs, pdamp, sdamp, top}],
  halfspace {...} or None (its ``top`` = ``total_depth``, below *all* profile layers), source
  'TOPL' | 'L', start, end, n_layers, omitted_above, omitted_below (profile layers outside
  Start / End Layer, to be marked as not shown)} and ``show`` (table columns);
* ``soilprop``: ``curves`` {name, g_strain, g, d_strain, d} (strain %, G/Gmax, damping %);
* ``3d`` / ``anim``: ``scene`` (:func:`~.state.model_scene`: ``node_id``, ``xyz``, ``faces``
  (n, 4; -1 pads triangles), ``tri`` (triangles for plotly Mesh3d), ``face_elem``,
  ``face_boundary``, ``edges`` (beams / springs), ``elem_*`` attributes and colours (ELECOLOR),
  ``node_visible`` (False: hidden by HIDENODE or outside the display volume -- draw no dot,
  label, bubble or vector), ``interaction``, ``fixed`` (SHOWDOF), ``mass`` / ``mass_dir``
  (SHOWMASS), ``selected`` (NODESEL) -- marker lists of visible nodes only --, ``bbox``,
  ``center``) and ``camera`` {basis (3x3: rows = screen right, screen up,
  toward the viewer, in model axes), center, extent, half}.  Screen coordinates of a point p are
  ``basis @ (p - center)`` (orthographic, D-UI-06); for plotly use ``up = basis[1]`` and the eye
  direction ``basis[2]``;
* ``anim`` also: ``animation`` {buffer_dir, frames (sequence), current, col, vmin, vmax, scale,
  layout} and ``frame`` (current frame on the scene nodes: ``value``, ``t``, ``rgb``, ``size``,
  ``colorbar`` for bubble / contour; ``vectors`` (n, 3 directions, 3 components) for vector;
  ``xyz`` deformed coordinates for deformed).  Other frames: :func:`~.state.animation_frame`
  with a :class:`~.state.FrameStore`.

**Settings rules** (spec 10 section 4.0).  A setting command acts on the active plot when it is
capable (:data:`~.state.CAPABILITY`), is ignored with a warning otherwise, and with no plot open
becomes a default of new plots (2D, for new SPECPLOT / THPLOT plots only -- a SOILPROPPLOT starts
from its own defaults: AXES, PLOTRANGE, XTITLE, YTITLE, YTITLE2, STIPPLE; 3D:
ELECOLOR, SHRINK, WIREFRAME, ELENUM, GROUPNUM, NODENUM, SHOWDOF, SHOWMASS, CNGVIEW, RSTVIEW,
CNGCENTER, RSTCENTER, DEBUG).  Line names and markers are global line properties.  Model-level
display state lives in ``model.ui_state`` (saved by SAVE, not by WRITE, D-UI-11): ``nodesel``
[node ids], ``hide_groups`` [g], ``hide_elements`` [[g, e]], ``hide_nodes`` [n] and
``display_volume`` [xmin, xmax, ymin, ymax, zmin, zmax] -- edited with NODESEL and the
WINDOWSETTINGS field form, which act on the model of the active 3D plot (a CUTPLOT may show a
model that is not the active one), else on the active model.  Display volume (spec 06 section
1.4, D-UI-11): an element is hidden only when all its nodes are outside the box; nodes outside
the box are hidden, and nodes of hidden elements are not in the scene at all.

**WINDOWSETTINGS field form** (SASSI-EDU extension of the settings dialogs so that every dialog
action has command text): ``WINDOWSETTINGS,<field>,<values>`` with ``TITLE,<text>``; layer plot
``START,<n>``, ``END,<n>`` (-1 deepest), ``SHOW,<THICK|WEIGHT|VP|VS|PDAMP|SDAMP>,<0|1>``; soil
properties ``SHOW,<MODULUS|DAMPING>,<0|1>``; animations ``SCALE,<s>``, ``RANGE,<min>,<max>``,
``DIRECTION,<X|Y|Z|ALL>``, ``UNDEFORMED,<0|1>``, ``FRAME,<k>``, ``FRAMEPAUSE,<ms>``; display
requests (model of the active 3D plot, else the active model) ``VOLUME,[6 bounds]``,
``HIDEGROUP,<g>``, ``SHOWGROUP,<g>``, ``HIDEELEM,<g>,<ids>``, ``SHOWELEM,<g>,<ids>``,
``HIDENODE,<ids>``, ``SHOWNODE,<ids>``, ``SHOWALL``.

**Views** (D-UI-12).  ``CNGVIEW,<rX>,<rY>,<rZ>,<px>,<py>,<zoom>``: screen = Rx(rX) Ry(rY) Rz(rZ)
(p - centre), degrees; the default (-60, 0, -45) is isometric with Z up; ``px``, ``py`` pan the
window in model length units; ``zoom`` 1 fits the bounding box of the element-connected nodes.

**Frames and animations** (requirements 5.8, D-FIL-03, D-UI-13..15).  Frame files: header
``nrows ncols`` then rows ``node v1 .. v(ncols-1)``; TFU frames hold Re/Im of X, Y, Z
(layout ``complex_xyz``), ACC / THD / RS frames X, Y, Z (``xyz``).  PROCFRAME writes the frame
store (``frame_00001.npy`` float32, ``nodes.npy``, ``index.json``) and an entry of ``SASSIani.xml``
in the per-user settings directory (``SASSI_EDU_SETTINGS_DIR`` overrides it; ``state.ani_db_path``
overrides the file).  A frame *directory* (MOTION / RELDISP write no list file) is read in
frame-number order (:func:`~.state.frame_sort_key`; not name order, which breaks at t >= 100 s);
a directory mixing sequences (RS01 / RS02, stress components) gives a warning.  ``<Col>`` of
BUBBLEPLOT / CONTOURPLOT counts the *data* columns (1 = the first value after the node id).
Colours: jet, clamped (dark blue below MnR, dark red above MxR); bubble size
``S_max max(0.2, t)``; vectors ``V_X = s (a_X, b_X, b_X)``, ``V_Y = s (b_Y, a_Y, b_Y)``,
``V_Z = s (b_Z, b_Z, a_Z)``; deformed ``x' = x + Scale u``.
"""
from __future__ import annotations

from .lines import Line, LineError
from .state import PlotError, PlotEvent, PlotState, plot_state

__all__ = ["Line", "LineError", "PlotError", "PlotEvent", "PlotState", "plot_state"]
