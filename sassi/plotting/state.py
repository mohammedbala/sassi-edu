"""Session plot state: line store, plot tabs, view settings, model scenes, frames and animations.

This module is the *data level* of the plotting commands (requirements 3.4.L, 5.7, 5.8; spec 06;
spec 10 section 4).  It never imports matplotlib: the GUI (plotly) and the headless renderer
(:mod:`sassi.plotting.render_mpl`) both draw from the JSON-serialisable dictionaries built here,
so a plot looks the same whichever front end shows it.  The package docstring
(:mod:`sassi.plotting`) describes the contract for the GUI.

Contents
--------
* :class:`PlotState` -- one per interpreter session (``plot_state(interp)``): line objects,
  open plots ("tabs") with the active one, default settings, shader options, palette overrides,
  listeners, the auto-render switch and the animation database location.
* :class:`Plot`, :class:`Settings2D`, :class:`View3D` -- one plot request and its settings.
* :func:`model_scene` -- nodes, faces, lines and markers of a model for MODELPLOT / NODEPLOT /
  CUTPLOT and the animations (spec 06 sections 2-4).
* :func:`layer_table`, :func:`soil_property_curves` -- LAYERPLOT / SOILPROPPLOT data.
* Frame files (D-FIL-03), the frame store written by PROCFRAME (requirements 5.8: one ``.npy`` per
  frame + ``index.json``) and the ``SASSIani.xml`` animation database (D-UI-15).
* Animation maths: frame sequence, jet colour map with clamping, bubble size (D-UI-14), vector
  rule (manual 7.3.5, D-UI-13), deformed coordinates.
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

from .lines import Line, LineError

PathLike = Union[str, Path]

# ======================================================================================
# Plot kinds, families and command capabilities (spec 06 section 0.1, spec 10 section 4.0)
# ======================================================================================
#: kind -> (family, caption)
PLOT_KINDS: Dict[str, Tuple[str, str]] = {
    "SPECPLOT": ("2d", "Spectrum Plot"),
    "THPLOT": ("2d", "Time History Plot"),
    "LAYERPLOT": ("layer", "Soil Layer Plot"),
    "SOILPROPPLOT": ("soilprop", "Soil Property Plot"),
    "MODELPLOT": ("3d", "Model Plot"),
    "NODEPLOT": ("3d", "Node Plot"),
    "CUTPLOT": ("3d", "Cut Plot"),
    "BUBBLEPLOT": ("anim", "Bubble Plot"),
    "VECTORPLOT": ("anim", "Vector Plot"),
    "CONTOURPLOT": ("anim", "Contour Plot"),
    "DEFORMPLOT": ("anim", "Displacement Plot"),
}
LINE_PLOTS = frozenset({"SPECPLOT", "THPLOT"})
ANIMATIONS = frozenset({"BUBBLEPLOT", "VECTORPLOT", "CONTOURPLOT", "DEFORMPLOT"})
THREE_D = frozenset({"MODELPLOT", "NODEPLOT", "CUTPLOT"}) | ANIMATIONS
ELEMENT_PLOTS = frozenset({"MODELPLOT", "CUTPLOT", "CONTOURPLOT", "DEFORMPLOT"})
ALL_KINDS = frozenset(PLOT_KINDS)

#: plot kinds each setting command acts on ("if the plot is capable"; others: ignored + warning)
CAPABILITY: Dict[str, frozenset] = {
    "AXES": LINE_PLOTS | {"SOILPROPPLOT"},
    "PLOTRANGE": LINE_PLOTS | {"SOILPROPPLOT"},
    "XTITLE": LINE_PLOTS | {"SOILPROPPLOT"},
    "YTITLE": LINE_PLOTS | {"SOILPROPPLOT"},
    "YTITLE2": LINE_PLOTS | {"SOILPROPPLOT"},
    "STIPPLE": LINE_PLOTS,
    "PLOTTITLE": ALL_KINDS - {"LAYERPLOT"},
    "ELECOLOR": frozenset({"MODELPLOT"}),
    "SHRINK": frozenset({"MODELPLOT"}),
    "WIREFRAME": frozenset({"MODELPLOT"}),
    "ELENUM": ELEMENT_PLOTS,
    "GROUPNUM": ELEMENT_PLOTS,
    "NODENUM": THREE_D,
    "SHOWDOF": frozenset({"MODELPLOT", "NODEPLOT"}),
    "SHOWMASS": frozenset({"MODELPLOT", "NODEPLOT"}),
    "CNGVIEW": THREE_D, "RSTVIEW": THREE_D, "CNGCENTER": THREE_D, "RSTCENTER": THREE_D,
    "DEBUG": THREE_D,
    "PAUSE": ANIMATIONS,
}
#: setting commands that become session defaults for new plots when no plot is active
DEFAULTABLE_2D = frozenset({"AXES", "PLOTRANGE", "XTITLE", "YTITLE", "YTITLE2", "STIPPLE"})
DEFAULTABLE_3D = frozenset({"ELECOLOR", "SHRINK", "WIREFRAME", "ELENUM", "GROUPNUM", "NODENUM", "SHOWDOF",
                            "SHOWMASS", "CNGVIEW", "RSTVIEW", "CNGCENTER", "RSTCENTER", "DEBUG"})

#: SHOWDOF labels -> DOF indices 0..5 (UX UY UZ ROTX ROTY ROTZ); spec 10 section 4.6
SHOWDOF_LABELS: Dict[str, Tuple[int, ...]] = {
    "X": (0,), "Y": (1,), "Z": (2,), "XX": (3,), "YY": (4,), "ZZ": (5,),
    "DISP": (0, 1, 2), "ROT": (3, 4, 5), "ALL": (0, 1, 2, 3, 4, 5), "NONE": (),
}
#: default isometric view (D-UI-12): rotations about X, Y, Z in degrees, pan, zoom (1 = fit)
DEFAULT_VIEW = (-60.0, 0.0, -45.0, 0.0, 0.0, 1.0)
#: types drawn as faces and as lines (DOF nodes only: the BEAMS / GENERAL K node is orientation)
FACE_TYPES = frozenset({1, 3, 4, 5})
LINE_TYPES = frozenset({2, 7, 9})
DOF_NODE_COUNT = {2: 2, 9: 2}
#: hexahedron faces (nodes 1-4 one face, 5-8 the opposite one; spec 08 section 4.4)
HEX_FACES = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))
#: element types carrying a material / a property number (ELECOLOR 2 / 3; others: default colour)
HAS_MATERIAL = frozenset({1, 2, 3, 4, 5})
HAS_PROPERTY = frozenset({2, 7, 9})


class PlotError(ValueError):
    """Invalid plot request (reported by the command handlers as an error)."""


# ======================================================================================
# Palettes (Options > Colors, manual 6.5.6; COLOR command)
# ======================================================================================
def _hex(rgb: Sequence[float]) -> str:
    return "#{:02x}{:02x}{:02x}".format(*[int(round(max(0.0, min(1.0, v)) * 255)) for v in rgb[:3]])


def _elem_palette() -> List[str]:
    """ElemPalette: 128 distinct colours (index ``((n-1) mod 128) + 1``, D-UI-09)."""
    base = ["#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#bcf60c",
            "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000"]
    out = list(base)
    k = 0
    while len(out) < 128:
        h = (k * 0.618033988749895) % 1.0
        s = 0.55 + 0.35 * ((k // 7) % 2)
        v = 0.95 - 0.25 * ((k // 3) % 2)
        i = int(h * 6) % 6
        f = h * 6 - int(h * 6)
        p, q, t = v * (1 - s), v * (1 - f * s), v * (1 - (1 - f) * s)
        rgb = [(v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q)][i]
        out.append(_hex(rgb))
        k += 1
    return out


_LINE_BASE = ["#e00000", "#0000e0", "#00a000", "#c000c0", "#00a0a0", "#e08000", "#606060", "#8000ff",
              "#a05000", "#ff60a0", "#008060", "#000080", "#808000", "#00c0ff", "#ff0060", "#404000"]
#: entries listed for the line palettes (a SPECPLOT / THPLOT shows up to 50 lines); line numbers
#: beyond continue the same sequence (:func:`line_color`), never cyclically
LINE_PALETTE_SIZE = 50
LINE_PALETTES = frozenset({"SPECLINES", "THLINES"})


def line_color(num: int) -> str:
    """Default colour of line ``num`` (1-based) in the SpecLines / THLines palettes.

    The manual's line palettes are "populated as lines with different line numbers are plotted":
    every line number has its own entry.  Numbers 1-16 are hand-picked high-contrast colours;
    later numbers walk the hue circle by the golden ratio (well-spread, non-repeating hues) with
    alternating lightness and saturation, so no two line numbers share a colour.
    """
    n = int(num)
    if 1 <= n <= len(_LINE_BASE):
        return _LINE_BASE[n - 1]
    k = max(n - len(_LINE_BASE), 1)
    h = (0.11 + k * 0.618033988749895) % 1.0
    s = (0.95, 0.65, 0.8)[k % 3]
    v = (0.85, 0.6, 0.45, 0.72)[(k // 3) % 4]
    i = int(h * 6) % 6
    f = h * 6 - int(h * 6)
    p, q, t = v * (1 - s), v * (1 - f * s), v * (1 - (1 - f) * s)
    return _hex([(v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q)][i])


PALETTES: Dict[str, List[str]] = {
    "ELEMPALETTE": _elem_palette(),
    "SPECLINES": [line_color(n) for n in range(1, LINE_PALETTE_SIZE + 1)],
    "THLINES": [line_color(n) for n in range(1, LINE_PALETTE_SIZE + 1)],
    # Node plot: 1 ordinary, 2 interaction, 3 fixed-DOF border, 4 mass border, 5 selected border
    "NODE": ["#000000", "#e00000", "#00b000", "#9000c0", "#0060ff"],
    # Element plot: 1 background, 2 outline, 3 labels, 4 default fill, 5 beams/springs, 6 mass X, 7 Y, 8 Z,
    # 9 fixed-DOF marker, 10 interaction node, 11 selected node
    "ELEMENT": ["#ffffff", "#202020", "#000000", "#b0b0b0", "#008000", "#e00000", "#00b000", "#0000e0",
                "#00b000", "#e00000", "#0060ff"],
    "CUT": ["#ffffff", "#404040", "#e00000"],
    "BUBBLE": ["#ffffff", "#000000"], "VECTOR": ["#ffffff", "#000000", "#e00000", "#00b000", "#0000e0"],
    "CONTOUR": ["#ffffff", "#202020"], "DEFORMED": ["#ffffff", "#c00000", "#909090", "#c8c8c8"],
    "SPEC": ["#ffffff", "#000000", "#b0b0b0"], "TIMEHIST": ["#ffffff", "#000000", "#b0b0b0"],
    "SOILLAYER": ["#e8e8e8", "#c8c8c8", "#a8a8a8", "#888888", "#686868", "#d8d8d8", "#b8b8b8", "#989898"],
    "SOILPROP": ["#ffffff", "#000000", "#e00000", "#00a000"],
    "UI": ["#000000", "#0000c0", "#008000", "#800080", "#c08000", "#e00000"],
}
PALETTE_NAMES = tuple(PALETTES)


def palette_color(name: str, num: int, overrides: Optional[Dict[str, Dict[int, str]]] = None) -> str:
    """Colour ``num`` (1-based) of palette ``name`` with COLOR overrides applied.

    Fixed palettes are cyclic (ElemPalette: index ``((n-1) mod 128) + 1``, D-UI-09); the line
    palettes SpecLines / THLines have one entry per line number (:func:`line_color`), so the
    lines of one plot never share a colour.
    """
    key = name.upper()
    if key in LINE_PALETTES:
        idx = int(num)
        if overrides and key in overrides and idx in overrides[key]:
            return overrides[key][idx]
        return line_color(idx)
    pal = PALETTES[key]
    n = len(pal)
    idx = ((int(num) - 1) % n) + 1
    if overrides and key in overrides and idx in overrides[key]:
        return overrides[key][idx]
    return pal[idx - 1]


def jet(t: np.ndarray) -> np.ndarray:
    """Jet colour map (dark blue -> blue -> cyan -> yellow -> red -> dark red) of ``t`` in [0, 1].

    Returns RGB in [0, 1], shape ``t.shape + (3,)``; values are clamped first, so data below the
    range take the dark-blue end and above the dark-red end (spec 06 section 8.7).
    """
    t = np.clip(np.asarray(t, dtype=float), 0.0, 1.0)

    def ramp(x):
        return np.clip(1.5 - np.abs(4.0 * t - x), 0.0, 1.0)
    return np.stack([ramp(3.0), ramp(2.0), ramp(1.0)], axis=-1)


def colormap_position(v: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    """``t = clamp((v - MnR)/(MxR - MnR), 0, 1)`` (spec 06 section 8.7)."""
    v = np.asarray(v, dtype=float)
    span = vmax - vmin
    if not span > 0:
        return np.where(v > vmin, 1.0, 0.0) if np.ndim(v) else float(v > vmin)
    return np.clip((v - vmin) / span, 0.0, 1.0)


def bubble_size(t: np.ndarray, smax: float) -> np.ndarray:
    """Bubble size ``S_max max(0.2, t)`` (D-UI-14)."""
    return float(smax) * np.maximum(0.2, np.asarray(t, dtype=float))


def colorbar_labels(vmin: float, vmax: float, n: int = 5) -> List[str]:
    """Five labels evenly spaced from Max (top) to Min (bottom), 5 decimals (spec 06 section 8.7)."""
    return [f"{v:.5f}" for v in np.linspace(vmax, vmin, n)]


# ======================================================================================
# Settings of one plot
# ======================================================================================
@dataclass
class Settings2D:
    """Axes and titles of a 2D line / soil-property plot (AXES, PLOTRANGE, XTITLE ... STIPPLE)."""
    xtitle: str = ""
    ytitle: str = ""
    ytitle2: str = ""
    major_x: bool = True          # AXES <MaxTickX>: thick grid lines at the major X ticks
    major_y: bool = True
    minor_x: bool = False         # AXES <MinTickX>: thin grid lines at the minor ticks ("Show Ticks")
    minor_y: bool = False
    log_x: bool = False
    log_y: bool = False
    xmin: Optional[float] = None  # PLOTRANGE; None = data extent
    xmax: Optional[float] = None
    ymin: Optional[float] = None
    ymax: Optional[float] = None
    stipple: bool = False         # STIPPLE: dash patterns per line

    def copy(self) -> "Settings2D":
        return Settings2D(**asdict(self))

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class View3D:
    """3D view and display state (CNGVIEW ... WIREFRAME, PAUSE; spec 06 sections 1-2, 8)."""
    rx: float = DEFAULT_VIEW[0]   # rotation about the screen X axis (deg)
    ry: float = DEFAULT_VIEW[1]
    rz: float = DEFAULT_VIEW[2]   # rotation about the model Z axis (deg), applied first
    px: float = 0.0               # pan, model length units along screen right / up
    py: float = 0.0
    zoom: float = 1.0             # 1 = fit the bounding box of the element-connected nodes
    center: Optional[List[float]] = None   # CNGCENTER; None = bounding-box centre (RSTCENTER)
    color_by: int = 1             # ELECOLOR 1 group, 2 material, 3 property
    node_labels: bool = False     # NODENUM
    elem_labels: bool = False     # ELENUM  (exclusive with group_labels)
    group_labels: bool = False    # GROUPNUM
    shrink: bool = False          # SHRINK (element plot)
    wireframe: bool = False       # WIREFRAME (element plot)
    show_dof: List[int] = field(default_factory=list)   # SHOWDOF: DOF indices 0..5 marked
    show_mass: bool = False       # SHOWMASS
    debug: bool = False           # DEBUG overlay (view values, animation info)
    paused: bool = False          # PAUSE (animations): False = running
    direction: str = "X"          # vector plot Output Direction X / Y / Z / ALL (Window Options)
    show_undeformed: bool = False  # deformed plot overlay (Window Options)
    view_rev: int = 0             # +1 at every CNGVIEW / RSTVIEW / CNGCENTER / RSTCENTER: the browser then
                                  # replaces a view rotated or zoomed with the mouse by the commanded one

    def copy(self) -> "View3D":
        d = asdict(self)
        return View3D(**d)

    def reset_view(self) -> None:
        """RSTVIEW: default isometric view, no pan, zoom 1 (D-UI-12)."""
        self.rx, self.ry, self.rz, self.px, self.py, self.zoom = DEFAULT_VIEW

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ShaderOptions:
    """SHADEROPTIONS,[points],[linew],[shrink],[scale] (spec 06 section 1.5; fractions, not %)."""
    points: float = 10.0          # Node/Bubble Node Size (maximum point size, px)
    linew: float = 0.02           # Element Outline Thickness (fraction of the element)
    shrink: float = 0.06          # Element Shrink (fraction)
    scale: float = 1.0            # Vector/Displacement Scale Factor

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def rotation_matrix(rx: float, ry: float, rz: float) -> np.ndarray:
    """View rotation ``R = Rx(rx) Ry(ry) Rz(rz)`` (degrees): screen = R (p - centre).

    Rows of R in model coordinates: R[0] screen right, R[1] screen up, R[2] toward the viewer
    (orthographic projection, D-UI-06).  The default (-60, 0, -45) shows Z up on the screen.
    """
    a, b, c = (math.radians(v) for v in (rx, ry, rz))
    Rx = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    Ry = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
    Rz = np.array([[math.cos(c), -math.sin(c), 0], [math.sin(c), math.cos(c), 0], [0, 0, 1]])
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):   # spurious Accelerate flags
        return Rx @ Ry @ Rz


# ======================================================================================
# Plot requests and events
# ======================================================================================
@dataclass
class Plot:
    """One open plot ("tab"): what was requested and how it is displayed.

    ``params`` holds the request: ``lines`` (SPECPLOT/THPLOT numbers in order), ``name``
    (SOILPROPPLOT), ``cut`` (CUTPLOT), ``start``/``end``/``show`` (LAYERPLOT), the animation request
    (``buffer_dir``, ``frames``, ``col``, ``vmin``, ``vmax``, ``scale``, ``current``) ...
    """
    id: int
    kind: str
    params: Dict[str, Any] = field(default_factory=dict)
    model: Optional[int] = None
    title: str = ""
    settings: Settings2D = field(default_factory=Settings2D)
    view: View3D = field(default_factory=View3D)
    image: Optional[str] = None          # last image written by the headless renderer

    @property
    def family(self) -> str:
        return PLOT_KINDS[self.kind][0]

    @property
    def caption(self) -> str:
        name = PLOT_KINDS[self.kind][1]
        if self.kind == "CUTPLOT":
            return f"Model {self.model} - Cut # {self.params.get('cut', 1)} Plot"
        if self.kind in THREE_D:
            return f"Model {self.model} - {name}"
        return name

    def to_dict(self) -> Dict[str, Any]:
        return {"id": self.id, "kind": self.kind, "family": self.family, "caption": self.caption,
                "model": self.model, "title": self.title, "params": jsonable(self.params),
                "settings": self.settings.to_dict(), "view": self.view.to_dict(), "image": self.image}


@dataclass
class PlotEvent:
    """What listeners receive (``PlotState.subscribe``).

    ``kind``: ``'open'`` (new plot request), ``'update'`` (settings or data of a plot changed),
    ``'close'``, ``'activate'`` (ACTIVATEPLOT brought a plot to the front), ``'lines'`` (line
    store changed; ``data['numbers']``), ``'capture'`` (CAPTUREPLOT wrote ``data['path']``), ``'progress'`` (PROCFRAME,
    ``data['fraction']``), ``'animations'`` (frame store / SASSIani.xml changed), ``'nodesel'``,
    ``'dialog'`` (a command asked for a dialog; ``data['dialog']``).
    """
    kind: str
    plot_id: Optional[int] = None
    data: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {"event": self.kind, "plot_id": self.plot_id, "data": jsonable(self.data)}


Listener = Callable[[PlotEvent], None]


def jsonable(obj: Any) -> Any:
    """Recursively convert numpy arrays / scalars, tuples, sets and Paths to JSON types."""
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set, frozenset)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        v = float(obj)
        return v if math.isfinite(v) else None
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    return obj


def settings_dir() -> Path:
    """Per-user settings directory (requirements 5.10); ``SASSI_EDU_SETTINGS_DIR`` overrides it."""
    env = os.environ.get("SASSI_EDU_SETTINGS_DIR")
    if env:
        return Path(env).expanduser()
    from ..cli import settings_dir as _sd
    return _sd()


# ======================================================================================
# PlotState
# ======================================================================================
class PlotState:
    """Session-level plotting state held by the interpreter (``interp.session['plot']``).

    Line objects last for the session (spec 10 section 1.2); several plots may be open, one is
    active (``active``).  Setting commands act on the active plot when it is capable, otherwise
    they are ignored with a warning; with no active plot the 2D settings become the defaults of
    the next SPECPLOT / THPLOT plots and the 3D settings those of the next 3D plots (SASSI-EDU
    convenience for ``.pre`` scripts; a soil-property plot always starts from its own defaults).

    ``auto_render`` (default True): every new or changed plot is rendered headless to an image
    in the model directory (console and batch mode, requirements UI-03).  A GUI that draws the
    plots itself sets it to False after subscribing.  ``dialog_handler`` is set by a GUI that can
    open the dialogs of commands given without arguments (D-UI-08: otherwise an error).
    """

    def __init__(self):
        self.lines: Dict[int, Line] = {}
        self.plots: Dict[int, Plot] = {}
        self.order: List[int] = []                 # tab order
        self.active: Optional[int] = None
        self.next_id = 1
        self.defaults2d = Settings2D()
        self.defaults3d = View3D()
        self.shader = ShaderOptions()
        self.palettes: Dict[str, Dict[int, str]] = {}
        self.auto_render = True
        self.image_dir: Optional[Path] = None       # None: model directory (MDL path) or CWD
        self.image_dpi = 100
        self.ani_db_path: Optional[Path] = None     # None: settings_dir() / SASSIani.xml
        self.dialog_handler: Optional[Callable[[str, Dict[str, Any]], bool]] = None
        self._listeners: List[Listener] = []

    # ------------------------------------------------------------------ listeners
    def subscribe(self, fn: Listener) -> None:
        """Register a callback receiving every :class:`PlotEvent`."""
        if fn not in self._listeners:
            self._listeners.append(fn)

    def unsubscribe(self, fn: Listener) -> None:
        if fn in self._listeners:
            self._listeners.remove(fn)

    def emit(self, kind: str, plot_id: Optional[int] = None, **data) -> PlotEvent:
        ev = PlotEvent(kind, plot_id, data)
        for fn in list(self._listeners):
            fn(ev)
        return ev

    def dialog(self, name: str, **context) -> bool:
        """Ask the GUI to open dialog ``name``; False when no GUI can (batch / console)."""
        if self.dialog_handler is None:
            return False
        ok = bool(self.dialog_handler(name, context))
        if ok:
            self.emit("dialog", self.active, dialog=name, context=context)
        return ok

    # ------------------------------------------------------------------ lines
    def set_line(self, line: Line) -> Line:
        """Store a line; an existing number is overwritten without warning (spec 06 section 5.1)."""
        self.lines[line.number] = line
        return line

    def get_lines(self, numbers: Sequence[int]) -> List[Line]:
        """Lines by number; :class:`LineError` for an undefined one (line maths, D-LIN-02)."""
        missing = [n for n in numbers if n not in self.lines]
        if missing:
            raise LineError("line(s) " + ", ".join(str(n) for n in missing) + " not defined (D-LIN-02)")
        return [self.lines[n] for n in numbers]

    def plots_showing(self, numbers: Iterable[int]) -> List[Plot]:
        s = set(numbers)
        return [p for p in self.plots.values() if p.kind in LINE_PLOTS and s & set(p.params.get("lines", []))]

    # ------------------------------------------------------------------ plots
    @property
    def active_plot(self) -> Optional[Plot]:
        return self.plots.get(self.active) if self.active is not None else None

    def open_plot(self, kind: str, params: Optional[Dict[str, Any]] = None, model: Optional[int] = None) -> Plot:
        """Create a new plot (tab) from the current defaults and make it active."""
        if kind not in PLOT_KINDS:
            raise PlotError(f"unknown plot kind {kind}")
        p = Plot(self.next_id, kind, dict(params or {}), model)
        self.next_id += 1
        if kind in LINE_PLOTS:
            # 2D settings given with no plot open are defaults of new spectrum / time-history plots
            p.settings = self.defaults2d.copy()
        else:
            # never carried into a soil-property plot: a PLOTRANGE or XTITLE meant for spectra
            # would give it a frequency range and label; strain spans 1e-4..10 %: log axis
            # (spec 06 section 7.3)
            p.settings = Settings2D(log_x=kind == "SOILPROPPLOT")
        p.view = self.defaults3d.copy()
        if kind in ANIMATIONS:
            p.view.paused = False
        self.plots[p.id] = p
        self.order.append(p.id)
        self.active = p.id
        return p

    def close_plot(self, plot_id: Optional[int] = None) -> Optional[Plot]:
        """CLOSEPLOT: close the active (or given) plot; the last opened remaining one becomes active."""
        pid = self.active if plot_id is None else plot_id
        if pid is None or pid not in self.plots:
            return None
        p = self.plots.pop(pid)
        self.order.remove(pid)
        if self.active == pid:
            self.active = self.order[-1] if self.order else None
        return p

    def activate(self, plot_id: int) -> Plot:
        if plot_id not in self.plots:
            raise PlotError(f"plot {plot_id} is not open")
        self.active = plot_id
        return self.plots[plot_id]

    def target(self, command: str) -> Tuple[Optional[Plot], str]:
        """Where a setting command applies: ``(plot, 'plot' | 'ignored' | 'defaults' | 'none')``."""
        p = self.active_plot
        if p is None:
            if command in DEFAULTABLE_2D or command in DEFAULTABLE_3D:
                return None, "defaults"
            return None, "none"
        if p.kind in CAPABILITY.get(command, ALL_KINDS):
            return p, "plot"
        return p, "ignored"

    # ------------------------------------------------------------------ images
    def image_path(self, plot: Plot, interp=None) -> Path:
        """Auto-render image name: ``<model name or 'sassi'>_plot<NN>_<kind>.png`` in the model
        directory (MDL path), else the working directory (``image_dir`` overrides)."""
        base_dir = self.image_dir
        name = "sassi"
        if interp is not None:
            m = interp.models.get(plot.model if plot.model is not None else interp.active_model, interp.model)
            if m.name:
                name = m.name
            if base_dir is None:
                base_dir = Path(m.path) if m.path else Path(interp.cwd)
        if base_dir is None:
            base_dir = Path.cwd()
        return Path(base_dir) / f"{name}_plot{plot.id:02d}_{plot.kind.lower()}.png"

    def ani_db(self) -> "AnimationDB":
        return AnimationDB(self.ani_db_path or settings_dir() / "SASSIani.xml")

    # ------------------------------------------------------------------ serialisation
    def to_dict(self, include_lines: bool = True) -> Dict[str, Any]:
        """JSON-serialisable snapshot of the whole state (GUI)."""
        d: Dict[str, Any] = {
            "active": self.active,
            "plots": [self.plots[i].to_dict() for i in self.order],
            "defaults2d": self.defaults2d.to_dict(), "defaults3d": self.defaults3d.to_dict(),
            "shader": self.shader.to_dict(),
            "palettes": {k: {str(n): c for n, c in v.items()} for k, v in self.palettes.items()},
            "auto_render": self.auto_render,
        }
        if include_lines:
            d["lines"] = {str(n): self.lines[n].to_dict() for n in sorted(self.lines)}
        else:
            d["lines"] = {str(n): {"number": n, "name": self.lines[n].name, "kind": self.lines[n].kind,
                                   "n": self.lines[n].n, "markers": self.lines[n].markers}
                          for n in sorted(self.lines)}
        return d

    def plot_data(self, plot: Union[int, Plot], interp=None, raw: bool = False) -> Dict[str, Any]:
        """Everything needed to draw ``plot`` (JSON-serialisable): see :func:`plot_data`."""
        p = self.plots[plot] if isinstance(plot, int) else plot
        return plot_data(self, p, interp, raw=raw)


def plot_state(interp) -> PlotState:
    """The session :class:`PlotState` of an interpreter (created on first use)."""
    st = interp.session.get("plot")
    if not isinstance(st, PlotState):
        st = PlotState()
        interp.session["plot"] = st
    return st


# ======================================================================================
# Model scenes (MODELPLOT, NODEPLOT, CUTPLOT, animations)
# ======================================================================================
def _cycle_unique(ids: Sequence[int]) -> List[int]:
    """Drop repeated nodes of a degenerate face (prisms / pyramids repeat node numbers)."""
    out: List[int] = []
    for n in ids:
        if n not in out:
            out.append(n)
    return out


def element_nodes(gtype: int, nodes: Sequence[int]) -> List[int]:
    """Nodes that carry DOFs (the BEAMS / GENERAL K node is an orientation node only)."""
    k = DOF_NODE_COUNT.get(gtype)
    ids = [n for n in (nodes[:k] if k else nodes) if n]
    return ids


def model_scene(model, color_by: int = 1, show_dof: Sequence[int] = (), show_mass: bool = False,
                cut: Optional[Iterable[Tuple[int, int]]] = None, apply_hide: bool = True,
                palettes: Optional[Dict[str, Dict[int, str]]] = None) -> Dict[str, Any]:
    """Geometry and markers of a model for the 3D plots (spec 06 sections 2-4, 9).

    Returned arrays (numpy; :func:`jsonable` converts them):

    * ``node_id`` (n,), ``xyz`` (n, 3): element-connected nodes (DOF nodes; K nodes excluded) in
      global coordinates; every other index below refers to these rows;
    * ``faces`` (nf, 4) node indices (-1 pads triangles), ``face_elem`` (nf,), ``face_boundary``
      (nf,) -- False for a SOLID face shared by two elements (hidden unless shrink/wireframe);
      ``tri`` (nt, 3) / ``tri_face`` (nt,) the same faces split into triangles (plotly Mesh3d);
    * ``edges`` (ne, 2), ``edge_elem`` (ne,): BEAMS, SPRING and GENERAL elements as lines;
      ``points`` / ``point_elem``: zero-length springs;
    * ``elem_group``, ``elem_id``, ``elem_type``, ``elem_mat``, ``elem_prop``, ``elem_etype``,
      ``elem_key`` (ELECOLOR key: group, material or property number; 0 = default colour),
      ``elem_color`` (hex), ``elem_centroid`` (ne, 3), ``elem_in_cut`` (bool);
    * ``node_visible`` (n,) bool: False for nodes hidden with HIDENODE or outside the display
      volume (no dot, label, bubble, vector or marker is drawn for them);
    * ``interaction`` (indices of interaction nodes), ``fixed`` (indices of nodes with a fixed DOF
      among ``show_dof``), ``fix`` (n, 6) fixity codes, ``mass`` (indices), ``mass_dir`` (k, 3)
      bool (translational mass in X, Y, Z), ``mass_rot`` (k,) bool, ``selected`` (NODESEL); the
      marker lists hold visible nodes only;
    * ``bbox`` [xmin, xmax, ymin, ymax, zmin, zmax], ``center`` (bounding-box centre).

    Hide requests stored in the model (``model.ui_state``, D-UI-11) are applied when
    ``apply_hide``: ``hide_groups`` [g...], ``hide_elements`` [[g, e]...], ``hide_nodes`` [n...] and
    ``display_volume`` [xmin, xmax, ymin, ymax, zmin, zmax].  Display volume (spec 06 section 1.4,
    "geometry outside the box is hidden"): an element is clipped only when all its nodes are
    outside (D-UI-11); the node rows are the nodes of the remaining elements, and a node outside
    the box is not visible even when an element drawn across the box boundary uses it.
    """
    ui = model.ui_state if apply_hide else {}
    hide_g = set(int(g) for g in ui.get("hide_groups", []))
    hide_e = set((int(g), int(e)) for g, e in ui.get("hide_elements", []))
    hide_n = set(int(n) for n in ui.get("hide_nodes", []))
    vol = ui.get("display_volume")
    cut_set = set((int(g), int(e)) for g, e in cut) if cut is not None else None
    rows = []
    dangling = 0
    for g, e in model.iter_elements():
        if g.id in hide_g or (g.id, e.id) in hide_e:
            continue
        dof_nodes = element_nodes(g.type, e.nodes)
        if not dof_nodes:
            continue
        if any(n not in model.nodes for n in dof_nodes):
            dangling += 1
            continue
        rows.append((g, e, dof_nodes))

    def coordinates(rows_):
        used_ = sorted(set(n for _, _, nn_ in rows_ for n in nn_))
        if not used_:
            return np.zeros(0, dtype=np.int64), np.zeros((0, 3))
        ids_, xyz_ = model.global_coordinates(used_)
        return np.asarray(ids_, dtype=np.int64), np.asarray(xyz_, dtype=float).reshape(-1, 3)

    ids, xyz = coordinates(rows)
    in_box = np.ones(len(ids), dtype=bool)
    if vol is not None and len(rows):
        # Model Display Volume (spec 06 section 1.4, D-UI-11): an element is clipped only when all
        # its nodes lie outside the box; the node list is then rebuilt from the remaining elements,
        # so nodes of clipped elements disappear from every 3D plot (and from the camera fit)
        lo = np.array([vol[0], vol[2], vol[4]], float)
        hi = np.array([vol[1], vol[3], vol[5]], float)
        inside = np.all((xyz >= lo) & (xyz <= hi), axis=1)
        index = {int(n): k for k, n in enumerate(ids)}
        rows = [r for r in rows if any(inside[index[n]] for n in r[2])]
        ids, xyz = coordinates(rows)
        in_box = np.all((xyz >= lo) & (xyz <= hi), axis=1)
    index = {int(n): k for k, n in enumerate(ids)}
    faces: List[List[int]] = []
    face_elem: List[int] = []
    face_key: List[Tuple[int, ...]] = []
    face_sig: List[Tuple[int, ...]] = []
    edges: List[Tuple[int, int]] = []
    edge_elem: List[int] = []
    points: List[int] = []
    point_elem: List[int] = []
    eg, ei, et, em, ep, ee, ekey, ecent, ecut = [], [], [], [], [], [], [], [], []
    for k, (g, e, nn) in enumerate(rows):
        eg.append(g.id), ei.append(e.id), et.append(g.type), em.append(e.mat), ep.append(e.prop), ee.append(e.etype)
        if color_by == 2:
            key = e.mat if g.type in HAS_MATERIAL else 0
        elif color_by == 3:
            key = e.prop if g.type in HAS_PROPERTY else 0
        else:
            key = g.id
        ekey.append(int(key))
        idx = [index[n] for n in nn]
        ecent.append(xyz[_cycle_unique(idx)].mean(axis=0))
        ecut.append(cut_set is not None and (g.id, e.id) in cut_set)
        if g.type == 1:
            full = list(e.nodes) + [0] * (8 - len(e.nodes))
            if len(e.nodes) < 8 or any(n == 0 for n in full[:8]):
                continue
            sig = tuple(sorted(set(full[:8])))      # coincident solids (backfill on excavated soil) share it
            for fc in HEX_FACES:
                f = _cycle_unique([full[i] for i in fc])
                if len(f) >= 3:
                    faces.append([index[n] for n in f])
                    face_elem.append(k)
                    face_key.append(tuple(sorted(f)))
                    face_sig.append(sig)
        elif g.type in FACE_TYPES:
            f = _cycle_unique(nn[:4])
            if len(f) >= 3:
                faces.append([index[n] for n in f])
                face_elem.append(k)
                face_key.append(())          # plates are always drawn
                face_sig.append(())
            elif len(f) == 2:
                edges.append((index[f[0]], index[f[1]]))
                edge_elem.append(k)
        elif g.type in LINE_TYPES:
            if len(nn) >= 2 and nn[0] != nn[1]:
                edges.append((index[nn[0]], index[nn[1]]))
                edge_elem.append(k)
            else:
                points.append(index[nn[0]])
                point_elem.append(k)
    # a solid face is on the boundary unless solids with *different* nodes share it (two neighbours);
    # solids on the same nodes (e.g. backfill elements on the excavated-soil elements) are one solid for
    # the picture: their faces are drawn once, by the last of them (the backfill or structure element is
    # usually defined after the excavated soil, so its group colour shows)
    sigs: Dict[Tuple[int, ...], set] = {}
    last: Dict[Tuple[int, ...], int] = {}
    for r, fk in enumerate(face_key):
        if fk:
            sigs.setdefault(fk, set()).add(face_sig[r])
            last[fk] = r
    boundary = np.array([(not fk) or (len(sigs[fk]) == 1 and last[fk] == r) for r, fk in enumerate(face_key)],
                        dtype=bool)
    F = np.full((len(faces), 4), -1, dtype=np.int64)
    for r, f in enumerate(faces):
        F[r, :len(f)] = f[:4]
    tri, tri_face = [], []
    for r, f in enumerate(faces):
        for j in range(1, len(f) - 1):
            tri.append((f[0], f[j], f[j + 1]))
            tri_face.append(r)
    pal = palettes or {}
    colors = [palette_color("ELEMPALETTE", kk, pal) if kk > 0 else palette_color("ELEMENT", 4, pal) for kk in ekey]
    # nodes: visibility (HIDENODE requests, Model Display Volume), then the node markers --
    # interaction, fixities, masses, selection -- of the visible nodes only
    keep_nodes = in_box.copy()
    if hide_n:
        keep_nodes &= np.array([int(n) not in hide_n for n in ids], dtype=bool)
    nodes = model.nodes
    inter = [k for k, n in enumerate(ids) if keep_nodes[k] and 0 in nodes[int(n)].flags]
    fix = np.array([nodes[int(n)].fix for n in ids], dtype=np.int64).reshape(len(ids), 6)
    sd = sorted(set(int(d) for d in show_dof))
    fixed = np.nonzero(fix[:, sd].any(axis=1) & keep_nodes)[0] if sd and len(ids) else np.zeros(0, dtype=np.int64)
    mass_idx, mass_dir, mass_rot = [], [], []
    if show_mass:
        for k, n in enumerate(ids):
            if not keep_nodes[k]:
                continue
            t = model.tmass.get(int(n))
            r = model.rmass.get(int(n))
            tdir = [bool(t and abs(t[j]) > 0) for j in range(3)]
            rot = bool(r and any(abs(v) > 0 for v in r))
            if any(tdir) or rot:
                mass_idx.append(k)
                mass_dir.append(tdir)
                mass_rot.append(rot)
    sel = [index[int(n)] for n in model.ui_state.get("nodesel", [])
           if int(n) in index and keep_nodes[index[int(n)]]]
    if len(ids):
        lo, hi = xyz.min(axis=0), xyz.max(axis=0)
        bbox = [lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]]
        center = 0.5 * (lo + hi)
    else:
        bbox = [0.0] * 6
        center = np.zeros(3)
    return {
        "node_id": np.asarray(ids, dtype=np.int64), "xyz": np.asarray(xyz, dtype=float),
        "node_visible": keep_nodes,
        "faces": F, "face_elem": np.asarray(face_elem, dtype=np.int64), "face_boundary": boundary,
        "tri": np.asarray(tri, dtype=np.int64).reshape(-1, 3), "tri_face": np.asarray(tri_face, dtype=np.int64),
        "edges": np.asarray(edges, dtype=np.int64).reshape(-1, 2), "edge_elem": np.asarray(edge_elem, dtype=np.int64),
        "points": np.asarray(points, dtype=np.int64), "point_elem": np.asarray(point_elem, dtype=np.int64),
        "elem_group": np.asarray(eg, dtype=np.int64), "elem_id": np.asarray(ei, dtype=np.int64),
        "elem_type": np.asarray(et, dtype=np.int64), "elem_mat": np.asarray(em, dtype=np.int64),
        "elem_prop": np.asarray(ep, dtype=np.int64), "elem_etype": np.asarray(ee, dtype=np.int64),
        "elem_key": np.asarray(ekey, dtype=np.int64), "elem_color": colors,
        "elem_centroid": np.asarray(ecent, dtype=float).reshape(-1, 3),
        "elem_in_cut": np.asarray(ecut, dtype=bool),
        "interaction": np.asarray(inter, dtype=np.int64), "fix": fix, "fixed": np.asarray(fixed, dtype=np.int64),
        "show_dof": sd, "mass": np.asarray(mass_idx, dtype=np.int64),
        "mass_dir": np.asarray(mass_dir, dtype=bool).reshape(-1, 3), "mass_rot": np.asarray(mass_rot, dtype=bool),
        "selected": np.asarray(sel, dtype=np.int64),
        "bbox": [float(v) for v in bbox], "center": [float(v) for v in center],
        "n_elements": len(rows), "dangling": dangling,
    }


def to_screen(points: np.ndarray, basis: np.ndarray, center: Sequence[float]) -> np.ndarray:
    """Screen coordinates ``basis @ (p - center)`` of points (n, 3).

    Written with element-wise products instead of ``@``: numpy 2.0 on the macOS Accelerate BLAS
    raises spurious "divide by zero / overflow in matmul" warnings for such small products.
    """
    P = np.asarray(points, dtype=float).reshape(-1, 3) - np.asarray(center, dtype=float)
    R = np.asarray(basis, dtype=float)
    return P[:, 0:1] * R[:, 0] + P[:, 1:2] * R[:, 1] + P[:, 2:3] * R[:, 2]


def view_frame(view: View3D, scene: Dict[str, Any]) -> Dict[str, Any]:
    """Resolved camera of a 3D plot: rotation, centre, screen extent (orthographic).

    ``basis``: rows right / up / toward-viewer in model coordinates (plotly: ``up`` = basis[1],
    ``eye`` direction = basis[2]); ``extent`` [x0, x1, y0, y1] of the screen window in model length
    units (fit of the rotated bounding box, divided by ``zoom``, shifted by the pan).
    """
    R = rotation_matrix(view.rx, view.ry, view.rz)
    c = np.asarray(view.center if view.center is not None else scene["center"], dtype=float)
    xyz = scene["xyz"]
    if len(xyz):
        s = to_screen(xyz, R, c)
        half = max(float(np.max(np.abs(s[:, 0]))), float(np.max(np.abs(s[:, 1]))), 1e-9) * 1.08
    else:
        half = 1.0
    half /= max(float(view.zoom), 1e-9)
    return {"basis": R, "center": c, "extent": [view.px - half, view.px + half, view.py - half, view.py + half],
            "half": half}


# ======================================================================================
# Soil layers and soil properties (LAYERPLOT, SOILPROPPLOT)
# ======================================================================================
def layer_table(model, start: int = 1, end: int = -1) -> Dict[str, Any]:
    """Rows of the soil layer plot (spec 06 section 6): the TOPL layers (or every L layer when
    TOPL is empty) from ``start`` to ``end`` (-1 = deepest), with the SITE half-space layer.

    ``top`` is the depth of a layer top below the free surface (all layers above it counted, shown
    or not); the half-space ``top`` = ``total_depth`` = sum of all layer thicknesses.
    ``omitted_above`` / ``omitted_below``: number of profile layers before ``start`` / after
    ``end`` that the plot does not show (the figure marks them).
    """
    site = model.options.record("SITE")
    hs = site.integer(5) if site is not None else None
    topl = list(model.topl)
    source = "TOPL"
    if not topl:                       # no TOPL list: every L layer except the half-space layer
        topl = [l for l in sorted(model.layers) if l != hs]
        source = "L"
    missing = [l for l in topl if l not in model.layers]
    seq = [l for l in topl if l in model.layers]
    if not seq:
        raise PlotError("no soil layers defined (L / TOPL)")
    n = len(seq)
    s = max(1, int(start or 1))
    e = n if (end is None or int(end) < 0 or int(end) > n) else int(end)
    if s > e:
        raise PlotError(f"Start Layer {s} is after End Layer {e} ({n} layers)")
    rows = []
    depth = sum(model.layers[l].thick for l in seq[:s - 1])
    for pos in range(s, e + 1):
        L = model.layers[seq[pos - 1]]
        rows.append({"index": pos, "layer": L.id, "thick": L.thick, "weight": L.weight, "vp": L.vp, "vs": L.vs,
                     "pdamp": L.pdamp, "sdamp": L.sdamp, "top": depth})
        depth += L.thick
    # the half-space starts below *all* layers of the profile, whatever End Layer shows
    total = sum(model.layers[l].thick for l in seq)
    half = None
    if hs and hs in model.layers:
        H = model.layers[hs]
        half = {"layer": H.id, "thick": 0.0, "weight": H.weight, "vp": H.vp, "vs": H.vs, "pdamp": H.pdamp,
                "sdamp": H.sdamp, "top": total}
    return {"layers": rows, "halfspace": half, "source": source, "missing": missing, "start": s, "end": e,
            "n_layers": n, "halfspace_layer": hs or 0, "total_depth": total,
            "omitted_above": s - 1, "omitted_below": n - e}


def soil_property_curves(model, name: str) -> Dict[str, Any]:
    """DYNP curves of the dynamic soil property ``name`` (case-sensitive; spec 06 section 7).

    Strains in %, G/Gmax and damping (%) as entered with DYNP (SHAKE convention).
    """
    pts = []
    names = set()
    for key, rec in model.options.entries("DYNP"):
        label, no = key if isinstance(key, tuple) else (None, key)
        names.add(label)
        if label == name:
            pts.append((int(no), rec))
    if not pts:
        avail = ", ".join(sorted(str(n) for n in names if n)) or "none"
        raise PlotError(f"dynamic soil property '{name}' not defined (names are case-sensitive; defined: {avail})")
    pts.sort(key=lambda t: t[0])
    g = [(rec.number(2), rec.number(3)) for _, rec in pts if rec.given(2) and rec.given(3)]
    d = [(rec.number(4), rec.number(5)) for _, rec in pts if rec.given(4) and rec.given(5)]
    return {"name": name,
            "g_strain": [float(a) for a, _ in g], "g": [float(b) for _, b in g],
            "d_strain": [float(a) for a, _ in d], "d": [float(b) for _, b in d],
            "points": len(pts)}


def cut_elements(interp, cut: int) -> Optional[List[Tuple[int, int]]]:
    """Elements (group, element) of cut ``cut`` from the session cut store (D-MDL-13).

    The cut commands (another work package) keep session-global cuts in ``interp.session['cuts']``;
    accepted layouts: ``{cut: iterable of (group, element)}`` or ``{cut: obj}`` where ``obj`` has an
    ``elements`` attribute or key.  Returns None when the cut is not defined.
    """
    store = interp.session.get("cuts")
    if store is None:
        return None
    entry = None
    if isinstance(store, dict):
        entry = store.get(cut, store.get(str(cut)))
    elif hasattr(store, "get"):
        entry = store.get(cut)
    if entry is None:
        return None
    if hasattr(entry, "elements"):
        entry = entry.elements
    elif isinstance(entry, dict) and "elements" in entry:
        entry = entry["elements"]
    out = []
    for item in entry:
        g, e = item
        out.append((int(g), int(e)))
    return out


# ======================================================================================
# Frame files (D-FIL-03) and list files
# ======================================================================================
@dataclass
class Frame:
    """An ASCII frame: header ``nrows ncols`` then rows ``node v1 ... v(ncols-1)``."""
    nodes: np.ndarray            # (nrows,) int
    values: np.ndarray           # (nrows, ncols-1) float
    header: List[str] = field(default_factory=list)

    @property
    def ncols(self) -> int:
        return self.values.shape[1] + 1


def read_frame(path: PathLike) -> Frame:
    """Read a frame file (D-FIL-03).  The first line is the header ``nrows ncols``; the data rows
    must have ``ncols`` numbers (node id = column 1).  A header that does not match the data is
    an error ("legacy header": fix it with MODFRAMES)."""
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise PlotError(f"{path.name}: empty frame file")
    head = lines[0].replace(",", " ").split()
    try:
        nrows, ncols = int(float(head[0])), int(float(head[1]))
    except (IndexError, ValueError):
        raise PlotError(f"{path.name}: the first line must be the header 'nrows ncols' (legacy header? use "
                        f"MODFRAMES)") from None
    rows = []
    for ln in lines[1:]:
        s = ln.strip()
        if not s:
            continue
        try:
            rows.append([float(t.replace("D", "E")) for t in s.replace(",", " ").split()])
        except ValueError:
            continue
    if len(rows) != nrows:
        raise PlotError(f"{path.name}: header says {nrows} rows, {len(rows)} found")
    bad = [k for k, r in enumerate(rows, start=1) if len(r) != ncols]
    if bad:
        raise PlotError(f"{path.name}: row {bad[0]} has {len(rows[bad[0] - 1])} values, header says {ncols} columns "
                        f"(the node id is column 1; fix legacy headers with MODFRAMES)")
    A = np.asarray(rows, dtype=float).reshape(nrows, ncols)
    return Frame(A[:, 0].astype(np.int64), A[:, 1:].copy(), head)


def write_frame(path: PathLike, nodes: np.ndarray, values: np.ndarray) -> Path:
    """Write a frame file in the D-FIL-03 layout (same text as MOTION's frames, 11 digits)."""
    path = Path(path)
    nodes = np.asarray(nodes)
    values = np.asarray(values, dtype=float).reshape(len(nodes), -1)
    out = [f"{len(nodes)} {values.shape[1] + 1}"]
    for nd, row in zip(nodes, values):
        out.append(" ".join([str(int(nd))] + [f"{v:.10e}" for v in row]))
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return path


def combine_frames(frames: Sequence[Frame], op: int) -> Frame:
    """FRAMECOMBIN (spec 09 section 2.11): cell by cell 0 SRSS, 1 sum, 2 average.

    All frames must have the same rows x columns and the same node ids (column 1, copied)."""
    if not frames:
        raise PlotError("no input frames")
    if op not in (0, 1, 2):
        raise PlotError("op must be 0 (SRSS), 1 (sum) or 2 (average)")
    f0 = frames[0]
    for k, f in enumerate(frames[1:], start=2):
        if f.values.shape != f0.values.shape:
            raise PlotError(f"frame {k} has {f.values.shape[0]} x {f.ncols} values, frame 1 has "
                            f"{f0.values.shape[0]} x {f0.ncols}")
        if not np.array_equal(f.nodes, f0.nodes):
            raise PlotError(f"frame {k} lists different node ids than frame 1")
    V = np.stack([f.values for f in frames])
    if op == 0:
        out = np.sqrt(np.sum(V * V, axis=0))
    elif op == 1:
        out = np.sum(V, axis=0)
    else:
        out = np.mean(V, axis=0)
    return Frame(f0.nodes.copy(), out, list(f0.header))


def frame_sort_key(path: PathLike) -> Tuple:
    """Order of the frames of a directory: by tag, component, then **frame number** (D-STR-08).

    Plain name order is wrong: the time field ``tt.ttt`` grows to 7 characters at t >= 100 s
    (``THD_100.000_00021`` sorts between ``THD_10.995_..`` and ``THD_11.000_..``), and the
    frequency field likewise at f >= 1000 Hz.  Files that do not follow the Table 3.2 pattern
    come after the matching ones, in name order.
    """
    name = Path(path).name
    info = frame_name_info(name)
    if info["frame"] is None:
        return (1, "", "", 0, 0.0, name)
    return (0, info["tag"], info["comp"] or "", info["frame"], info["value"], name)


def frame_directory_notes(files: Sequence[PathLike]) -> List[str]:
    """Warnings for a frame *directory* that mixes sequences (several tags such as RS01 / RS02
    damping sets, several stress components, or files outside the frame-name pattern): the
    animation then plays one sequence after the other and store index k is not the frame number."""
    infos = [frame_name_info(Path(f).name) for f in files]
    notes = []
    seqs = sorted(set((i["tag"], i["comp"] or "") for i in infos if i["frame"] is not None))
    if len(seqs) > 1:
        shown = ", ".join(t + (f"/{c}" if c else "") for t, c in seqs[:6]) + (" ..." if len(seqs) > 6 else "")
        notes.append(f"the directory mixes {len(seqs)} frame sequences ({shown}); they are stored one after the "
                     f"other -- use a list file or one directory per sequence")
    odd = [Path(f).name for f, i in zip(files, infos) if i["frame"] is None]
    if odd and len(odd) < len(infos):
        notes.append(f"{len(odd)} file(s) do not follow the frame-name pattern <tag>_<value>_<n> (e.g. {odd[0]}); "
                     f"they are stored after the numbered frames")
    return notes


def read_frame_list(path: PathLike) -> List[Path]:
    """Animation frame list (``*.dispani``, ``*.tfiani`` ... D-FIL-03): line 1 ignored, then one frame
    path per line, relative to the list file's folder.  A *directory* is accepted as well
    (SASSI-EDU extension; MOTION/RELDISP write no list file): its frame files ordered by
    :func:`frame_sort_key` (tag, component, frame number), so store index k is frame k."""
    path = Path(path)
    if path.is_dir():
        files = sorted((p for p in path.iterdir() if p.is_file() and not p.name.startswith(".")
                        and p.suffix.lower() not in (".npy", ".json", ".xml", ".png", ".bmp")), key=frame_sort_key)
        if not files:
            raise PlotError(f"{path}: no frame files in the directory")
        return files
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    out = []
    for ln in lines[1:]:
        s = ln.strip().strip('"')
        if not s:
            continue
        s = s.replace("\\", "/") if ("\\" in s and os.sep == "/") else s
        p = Path(os.path.expanduser(s))
        out.append(p if p.is_absolute() else path.parent / p)
    if not out:
        raise PlotError(f"{path.name}: no frame files listed (line 1 is a header and is ignored)")
    return out


def modframes_text(text: str, cols: int) -> Tuple[str, bool]:
    """MODFRAMES: replace the second number of the frame header (line 1) by ``cols``.

    A one-number legacy header ``nrows`` gets ``cols`` appended.  Returns ``(text, changed)``."""
    lines = text.split("\n")
    if not lines or not lines[0].strip():
        raise PlotError("empty frame file")
    toks = lines[0].replace(",", " ").split()
    new = [toks[0], str(int(cols))] + toks[2:] if len(toks) >= 2 else [toks[0], str(int(cols))]
    newline = " ".join(new)
    changed = newline != lines[0].strip()
    lines[0] = newline
    return "\n".join(lines), changed


_FRAME_NAME = re.compile(r"^(?P<tag>[A-Za-z]+\d*)_(?P<val>-?\d+(?:\.\d*)?)_(?P<num>\d+)(?:_(?P<comp>\w+))?$")


def frame_name_info(name: str) -> Dict[str, Any]:
    """Tag, time / frequency and frame number from a Table 3.2 frame name
    (``ACC_00.000_00001``, ``TFU_000.02_00001``, ``RS01_000.10_00001``, ``stress_00.000_00001_sig``)."""
    stem = Path(name).name
    m = _FRAME_NAME.match(stem)
    if not m:
        return {"tag": stem.split("_")[0] if "_" in stem else stem, "value": None, "frame": None, "comp": None}
    return {"tag": m.group("tag"), "value": float(m.group("val")), "frame": int(m.group("num")),
            "comp": m.group("comp")}


#: tag of the steady-state harmonic frames written by HARMFRAME (value field = phase angle wt in degrees)
HARM_TAG = "HARM"


def harmonic_frame_name(phase_deg: float, k: int) -> str:
    """``HARM_<phase>_<k>`` (Table 3.2 pattern; the value field is the phase angle wt in degrees)."""
    return f"{HARM_TAG}_{phase_deg:05.1f}_{k:05d}"


def harmonic_frames(H: np.ndarray, nframes: int) -> Tuple[np.ndarray, np.ndarray]:
    """Steady-state harmonic motion over one period (HARMFRAME, SASSI-EDU extension).

    ``H`` (n, 3) complex: the transfer functions of X, Y, Z of n nodes at one frequency (per unit control
    motion, time factor ``e^{i w t}`` as in the convolution ``a(t) = IFFT[H A]``).  Frame k (1-based) is the
    displacement at ``wt = phi_k = 2 pi (k - 1) / N``:

        u_k = Re(H e^{i phi_k}) = Re(H) cos(phi_k) - Im(H) sin(phi_k)

    so a DOF with ``H = |H| e^{i theta}`` moves as ``|H| cos(wt + theta)`` while the control motion moves as
    ``cos(wt)``.  Returns ``(phases in degrees (N,), U (N, n, 3))``."""
    if int(nframes) < 1:
        raise PlotError("the number of frames must be >= 1")
    H = np.asarray(H, dtype=complex).reshape(-1, 3)
    phi = 2.0 * np.pi * np.arange(int(nframes)) / int(nframes)
    U = np.real(H[None, :, :] * np.exp(1j * phi)[:, None, None])
    return np.degrees(phi), U


# ======================================================================================
# Frame store (PROCFRAME; requirements 5.8) and SASSIani.xml (D-UI-15)
# ======================================================================================
STORE_INDEX = "index.json"
ANITYPES = {0: "Bubble", 1: "Vector", 2: "Contour", 3: "Time History"}


def frame_layout(ncols: int, tag: str) -> str:
    """Meaning of the data columns: 'complex_xyz' (Re/Im of X Y Z, TFU frames), 'xyz', 'scalar' or 'columns'."""
    nd = ncols - 1
    if nd == 6 and tag.upper().startswith("TFU"):
        return "complex_xyz"
    if nd == 3:
        return "xyz"
    if nd == 1:
        return "scalar"
    if nd == 6:
        return "complex_xyz"
    return "columns"


def process_frames(listfile: PathLike, buffer_dir: PathLike, description: str = "", anitype: int = 0,
                   progress: Optional[Callable[[float, str], None]] = None) -> Dict[str, Any]:
    """PROCFRAME: convert the frames of a list file into the frame store (one ``.npy`` per frame).

    Store layout (requirements 5.8): ``frame_00001.npy`` ... (float32, rows x data columns, the
    "GPU-ready" single-precision frames of the manual), ``nodes.npy`` (int64, when every frame has
    the same node list; otherwise ``nodes_00001.npy`` per frame) and ``index.json`` (description,
    type, list file, frame count, data columns, layout, per-frame source / tag / time or
    frequency / min / max, global min / max, ``notes``: warnings such as a directory that mixes
    frame sequences, :func:`frame_directory_notes`).  Users should not edit the store, only delete it.
    """
    files = read_frame_list(listfile)
    notes = frame_directory_notes(files) if Path(listfile).is_dir() else []
    out = Path(buffer_dir)
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("frame_*.npy"):
        old.unlink()
    for old in out.glob("nodes*.npy"):
        old.unlink()
    if (out / STORE_INDEX).exists():              # a failed run must not leave a stale index behind
        (out / STORE_INDEX).unlink()
    entries = []
    first_nodes = None
    shared = True
    ncols = None
    gmin = gmax = None
    tag0 = ""
    n = len(files)
    for k, fp in enumerate(files, start=1):
        if not fp.exists():
            raise PlotError(f"frame file {fp} (list entry {k}) not found")
        fr = read_frame(fp)
        if ncols is None:
            ncols = fr.ncols
            tag0 = frame_name_info(fp.name)["tag"]
        elif fr.ncols != ncols:
            raise PlotError(f"{fp.name}: {fr.ncols} columns, the first frame has {ncols}")
        if first_nodes is None:
            first_nodes = fr.nodes
        elif shared and not np.array_equal(fr.nodes, first_nodes):
            shared = False
            for j, e in enumerate(entries, start=1):         # write the node lists of the earlier frames
                np.save(out / f"nodes_{j:05d}.npy", first_nodes)
        np.save(out / f"frame_{k:05d}.npy", fr.values.astype(np.float32))
        if not shared:
            np.save(out / f"nodes_{k:05d}.npy", fr.nodes)
        fmin = fr.values.min(axis=0) if len(fr.nodes) else np.zeros(fr.ncols - 1)
        fmax = fr.values.max(axis=0) if len(fr.nodes) else np.zeros(fr.ncols - 1)
        gmin = fmin if gmin is None else np.minimum(gmin, fmin)
        gmax = fmax if gmax is None else np.maximum(gmax, fmax)
        info = frame_name_info(fp.name)
        entries.append({"k": k, "file": f"frame_{k:05d}.npy", "source": str(fp), "name": fp.name,
                        "tag": info["tag"], "value": info["value"], "frame": info["frame"], "comp": info["comp"],
                        "nrows": int(len(fr.nodes)), "min": fmin.tolist(), "max": fmax.tolist()})
        if progress is not None:
            progress(k / n, fp.name)
    if shared and first_nodes is not None:
        np.save(out / "nodes.npy", first_nodes)
    index = {"format": "SASSI-EDU frame store", "version": 1, "description": description,
             "anitype": int(anitype), "type": ANITYPES.get(int(anitype), str(anitype)),
             "listfile": str(listfile), "nframes": n, "ncols": int(ncols), "data_columns": int(ncols) - 1,
             "layout": frame_layout(int(ncols), tag0), "tag": tag0, "shared_nodes": bool(shared),
             "min": gmin.tolist(), "max": gmax.tolist(), "notes": notes, "frames": entries}
    (out / STORE_INDEX).write_text(json.dumps(index, indent=1), encoding="utf-8")
    return index


class FrameStore:
    """Read access to a processed frame store (BufferDir)."""

    def __init__(self, directory: PathLike):
        self.dir = Path(directory)
        p = self.dir / STORE_INDEX
        if not p.exists():
            raise PlotError(f"{self.dir} is not a processed frame directory (no {STORE_INDEX}; run PROCFRAME)")
        self.index = json.loads(p.read_text(encoding="utf-8"))
        self._nodes = None

    @property
    def nframes(self) -> int:
        return int(self.index["nframes"])

    @property
    def data_columns(self) -> int:
        return int(self.index["data_columns"])

    @property
    def layout(self) -> str:
        return self.index.get("layout", "columns")

    def nodes(self, k: int) -> np.ndarray:
        if self.index.get("shared_nodes", True):
            if self._nodes is None:
                self._nodes = np.load(self.dir / "nodes.npy")
            return self._nodes
        return np.load(self.dir / f"nodes_{k:05d}.npy")

    def values(self, k: int) -> np.ndarray:
        """Frame k (1-based) as float64 (rows x data columns)."""
        if not 1 <= k <= self.nframes:
            raise PlotError(f"frame {k} outside 1..{self.nframes}")
        return np.load(self.dir / self.index["frames"][k - 1]["file"]).astype(float)

    def column_range(self, frames: Sequence[int], col: int) -> Tuple[float, float]:
        """Min / max of data column ``col`` (1-based) over the given frames."""
        fr = self.index["frames"]
        lo = min(fr[k - 1]["min"][col - 1] for k in frames)
        hi = max(fr[k - 1]["max"][col - 1] for k in frames)
        return float(lo), float(hi)

    def frame_label(self, k: int) -> str:
        e = self.index["frames"][k - 1]
        tag, val = e.get("tag", ""), e.get("value")
        if val is None:
            return f"frame {k}: {e.get('name', '')}"
        if str(tag).upper() == HARM_TAG:
            return f"frame {k}: ωt = {val:g}°"
        unit = "Hz" if str(tag).upper().startswith(("TFU", "RS")) else "s"
        return f"frame {k}: {tag} {val:g} {unit}"


class AnimationDB:
    """``SASSIani.xml``: ``<SASSIani><Animation description directory type frames start end stride
    scale cmin cmax listfile/>...</SASSIani>`` (D-UI-15) in the per-user settings directory."""

    ATTRS = ("description", "directory", "type", "frames", "start", "end", "stride", "scale", "cmin", "cmax",
             "listfile")

    def __init__(self, path: PathLike):
        self.path = Path(path)

    def entries(self) -> List[Dict[str, str]]:
        if not self.path.exists():
            return []
        try:
            root = ET.parse(str(self.path)).getroot()
        except ET.ParseError:
            return []
        return [dict(el.attrib) for el in root.findall("Animation")]

    def find(self, directory: PathLike) -> Optional[Dict[str, str]]:
        key = str(Path(directory).resolve())
        for e in self.entries():
            if str(Path(e.get("directory", "")).resolve()) == key:
                return e
        return None

    def _write(self, entries: List[Dict[str, str]]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        root = ET.Element("SASSIani")
        for e in entries:
            ET.SubElement(root, "Animation", {k: str(e[k]) for k in self.ATTRS if k in e})
        tree = ET.ElementTree(root)
        tmp = self.path.with_name(self.path.name + ".tmp")
        tree.write(str(tmp), encoding="utf-8", xml_declaration=True)
        os.replace(tmp, self.path)

    def upsert(self, entry: Dict[str, Any]) -> None:
        """Add an animation, replacing an entry with the same directory."""
        key = str(Path(entry["directory"]).resolve())
        keep = [e for e in self.entries() if str(Path(e.get("directory", "")).resolve()) != key]
        keep.append({k: str(v) for k, v in entry.items() if k in self.ATTRS})
        self._write(keep)

    def remove(self, directory: PathLike, delete_files: bool = False) -> bool:
        """Remove Animation: delete the entry (and, with ``delete_files``, the processed frames)."""
        key = str(Path(directory).resolve())
        ents = self.entries()
        keep = [e for e in ents if str(Path(e.get("directory", "")).resolve()) != key]
        if len(keep) == len(ents):
            return False
        self._write(keep)
        if delete_files and Path(directory).is_dir():
            shutil.rmtree(directory)
        return True


def frame_sequence(start: int, end: int, stride: int, nframes: int) -> List[int]:
    """Frames ``start, start+stride, ...`` up to ``end`` (spec 06 section 8.5; 1-based).

    Validation: ``1 <= start <= end <= nframes`` and ``stride >= 1``."""
    if stride < 1:
        raise PlotError(f"stride {stride} must be >= 1")
    if not 1 <= start <= end <= nframes:
        raise PlotError(f"frame range {start}..{end} must satisfy 1 <= start <= end <= {nframes}")
    return list(range(int(start), int(end) + 1, int(stride)))


# ======================================================================================
# Animation maths (spec 06 sections 8.7-8.10)
# ======================================================================================
def vector_components(values: np.ndarray, layout: str) -> Tuple[np.ndarray, np.ndarray]:
    """Real and imaginary parts (n, 3) of the X, Y, Z data of a frame."""
    v = np.asarray(values, dtype=float)
    if layout == "complex_xyz":
        return v[:, [0, 2, 4]], v[:, [1, 3, 5]]
    re = np.zeros((len(v), 3))
    m = min(3, v.shape[1])
    re[:, :m] = v[:, :m]
    return re, np.zeros_like(re)


def vector_rule(re: np.ndarray, im: np.ndarray, scale: float = 1.0) -> np.ndarray:
    """Drawn vectors (n, 3 directions, 3 components) of complex nodal data (manual 7.3.5, D-UI-13).

    ``V_X = Scale (a_X, b_X, b_X)``, ``V_Y = Scale (b_Y, a_Y, b_Y)``, ``V_Z = Scale (b_Z, b_Z, a_Z)``:
    a real component lies along its own axis; the imaginary part tilts it equally into the other
    two axes (an X value 5 + 0.3i is drawn as (5, 0.3, 0.3))."""
    re = np.asarray(re, dtype=float)
    im = np.asarray(im, dtype=float)
    V = np.repeat(im[:, :, None], 3, axis=2)          # V[n, d, :] = b_d
    for d in range(3):
        V[:, d, d] = re[:, d]
    return float(scale) * V


def deformed_coordinates(xyz: np.ndarray, u: np.ndarray, scale: float = 1.0) -> np.ndarray:
    """``x' = x + Scale u`` (spec 06 section 8.10); nodes without data have u = 0."""
    return np.asarray(xyz, dtype=float) + float(scale) * np.asarray(u, dtype=float)


def frame_on_scene(store: FrameStore, k: int, node_ids: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Frame k aligned with the scene nodes: ``(values (n, ncol) with NaN for missing nodes, found)``."""
    nodes = store.nodes(k)
    vals = store.values(k)
    pos = {int(n): r for r, n in enumerate(nodes)}
    out = np.full((len(node_ids), vals.shape[1]), np.nan)
    found = np.zeros(len(node_ids), dtype=bool)
    for i, n in enumerate(node_ids):
        r = pos.get(int(n))
        if r is not None:
            out[i] = vals[r]
            found[i] = True
    return out, found


def animation_frame(plot: Plot, scene: Dict[str, Any], store: FrameStore, k: Optional[int] = None) -> Dict[str, Any]:
    """Data of one animation frame on the scene nodes (current frame when ``k`` is None)."""
    p = plot.params
    k = int(p["current"] if k is None else k)
    vals, found = frame_on_scene(store, k, scene["node_id"])
    out: Dict[str, Any] = {"frame": k, "label": store.frame_label(k), "found": found}
    if plot.kind in ("BUBBLEPLOT", "CONTOURPLOT"):
        col = int(p["col"])
        v = vals[:, col - 1]
        t = colormap_position(np.nan_to_num(v, nan=p["vmin"]), p["vmin"], p["vmax"])
        out.update(value=v, t=np.where(found, t, np.nan), rgb=jet(t),
                   colorbar=colorbar_labels(p["vmin"], p["vmax"]))
        if plot.kind == "BUBBLEPLOT":
            out["size"] = np.where(found, bubble_size(t, p.get("points", 10.0)), 0.0)
    elif plot.kind == "VECTORPLOT":
        re, im = vector_components(np.nan_to_num(vals), store.layout)
        out["vectors"] = vector_rule(re, im, p["scale"])
        out["direction"] = plot.view.direction
    elif plot.kind == "DEFORMPLOT":
        re, _ = vector_components(np.nan_to_num(vals), store.layout)
        out["xyz"] = deformed_coordinates(scene["xyz"], re, p["scale"])
    return out


# ======================================================================================
# Export Table (File > Export Table, D-UI-17)
# ======================================================================================
def export_table_csv(state: PlotState, plot: Plot, path: PathLike, interp=None) -> Path:
    """Write the data of a spectrum / time-history / soil-property plot as CSV (D-UI-17).

    Line plots: header ``x name, line names``, rows on the union grid of the plotted lines (same
    resampling as WRITESPEC).  Soil-property plots: the four DYNP columns
    ``Strain, Mod. Red., Strain, Damp``.
    """
    import csv
    from .lines import on_union
    path = Path(path)
    if plot.kind in LINE_PLOTS:
        lines = [state.lines[n] for n in plot.params.get("lines", []) if n in state.lines]
        if not lines:
            raise PlotError("the plot shows no line")
        X, Y = on_union(lines)
        xname = plot.settings.xtitle or ("Frequency" if plot.kind == "SPECPLOT" else "Time")
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow([xname] + [L.name or f"line {L.number}" for L in lines])
            for k in range(len(X)):
                w.writerow([repr(float(X[k]))] + [repr(float(v)) for v in Y[:, k]])
        return path
    if plot.kind == "SOILPROPPLOT":
        if interp is None:
            raise PlotError("soil property data needs the interpreter (models)")
        cv = soil_property_curves(interp.models[plot.model], plot.params["name"])
        n = max(len(cv["g"]), len(cv["d"]))
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["Strain", "Mod. Red.", "Strain", "Damp"])
            for k in range(n):
                row = [cv["g_strain"][k], cv["g"][k]] if k < len(cv["g"]) else ["", ""]
                row += [cv["d_strain"][k], cv["d"][k]] if k < len(cv["d"]) else ["", ""]
                w.writerow(row)
        return path
    raise PlotError(f"Export Table is available for spectrum, time-history and soil-property plots, not {plot.kind}")


# ======================================================================================
# Plot data (what a front end draws)
# ======================================================================================
def line_plot_data(state: PlotState, plot: Plot) -> Dict[str, Any]:
    """SPECPLOT / THPLOT data: lines with colour (SpecLines / THLines by line number), markers,
    dash (STIPPLE) and the resolved extent (PLOTRANGE or the data extent)."""
    pal = "SPECLINES" if plot.kind == "SPECPLOT" else "THLINES"
    dashes = ["solid", "dashed", "dotted", "dashdot"]
    out_lines = []
    xs, ys = [], []
    for k, n in enumerate(plot.params.get("lines", [])):
        L = state.lines.get(n)
        if L is None:
            continue             # deleted since: ignored like an unknown number
        out_lines.append({"number": n, "name": L.name, "x": L.x, "y": L.y, "markers": L.markers,
                          "color": palette_color(pal, n, state.palettes),
                          "dash": dashes[k % len(dashes)] if plot.settings.stipple else "solid"})
        xs.append(L.x)
        ys.append(L.y)
    s = plot.settings
    ext = [None, None, None, None]
    if xs:
        X = np.concatenate(xs)
        Y = np.concatenate(ys)
        if s.log_x:
            X = X[X > 0]
        if s.log_y:
            Y = Y[Y > 0]
        if len(X):
            ext[0], ext[1] = float(X.min()), float(X.max())
        if len(Y):
            ext[2], ext[3] = float(Y.min()), float(Y.max())
    for i, v in enumerate((s.xmin, s.xmax, s.ymin, s.ymax)):
        if v is not None:
            ext[i] = float(v)
    return {"lines": out_lines, "extent": ext, "settings": s.to_dict()}


def plot_data(state: PlotState, plot: Plot, interp=None, raw: bool = False) -> Dict[str, Any]:
    """JSON-serialisable description of ``plot`` (the contract of :mod:`sassi.plotting`).

    Keys: ``id kind family caption title model settings view shader`` plus, by family:
    ``lines extent`` (2d), ``table show`` (layer), ``curves`` (soilprop), ``scene camera`` (3d) and
    ``scene camera animation frame`` (animations).  3D and model-based plots need ``interp``.
    ``raw=True`` keeps numpy arrays (headless renderer); otherwise everything is converted by
    :func:`jsonable` (lists, None for NaN).
    """
    fin = (lambda x: x) if raw else jsonable
    d: Dict[str, Any] = {"id": plot.id, "kind": plot.kind, "family": plot.family, "caption": plot.caption,
                         "title": plot.title, "model": plot.model, "params": plot.params,
                         "settings": plot.settings.to_dict(), "view": plot.view.to_dict(),
                         "shader": state.shader.to_dict()}
    fam = plot.family
    if fam == "2d":
        d.update(line_plot_data(state, plot))
        return fin(d)
    if interp is None:
        raise PlotError(f"{plot.kind} data needs the interpreter (models)")
    model = interp.models.get(plot.model)
    if model is None:
        raise PlotError(f"model {plot.model} is no longer in memory")
    if fam == "layer":
        d["table"] = layer_table(model, plot.params.get("start", 1), plot.params.get("end", -1))
        d["show"] = plot.params.get("show", {})
        return fin(d)
    if fam == "soilprop":
        d["curves"] = soil_property_curves(model, plot.params["name"])
        return fin(d)
    cut = None
    if plot.kind == "CUTPLOT":
        cut = cut_elements(interp, int(plot.params.get("cut", 1))) or []
    v = plot.view
    scene = model_scene(model, color_by=v.color_by, show_dof=v.show_dof, show_mass=v.show_mass, cut=cut,
                        palettes=state.palettes)
    cam = view_frame(v, scene)
    d["scene"] = scene
    d["camera"] = cam
    if fam == "anim":
        store = FrameStore(plot.params["buffer_dir"])
        d["animation"] = {k: plot.params[k] for k in ("buffer_dir", "frames", "current", "col", "vmin", "vmax",
                                                         "scale", "layout") if k in plot.params}
        d["frame"] = animation_frame(plot, scene, store)
    return fin(d)
