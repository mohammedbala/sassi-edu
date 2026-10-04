"""Preview pictures of the examples gallery (Learn tab, Model > Open Example): one PNG per ``examples/*.pre``.

Usage (from the project root)::

    python -m sassi.ui.thumbnails                         # every example -> sassi/ui/static/examples/<name>.png
    python -m sassi.ui.thumbnails ex01_surface_stick      # only these examples
    python -m sassi.ui.thumbnails --out DIR               # another output folder

A picture shows what clicking it opens in the GUI.  The example is copied into a temporary folder
(:func:`sassi.ui.lessons.copy_example`, deleted afterwards) and its model part runs there
(:func:`sassi.ui.learn.example_model_lines`: the ``.pre`` without CHECK, AFWRITE, RUN<MODULE>, WRITE ...;
no module runs, a fraction of a second per example), then the view command
(:func:`sassi.ui.learn.example_view_command`): ``MODELPLOT``, or ``LAYERPLOT`` for an example without a
structure (example 4, the free field).  The plot data of that plot (:func:`sassi.plotting.state.plot_data`,
the data the GUI and CAPTUREPLOT draw) is drawn here without the axes, triad and title:

* **MODELPLOT** -- the default isometric view (CNGVIEW -60, 0, -45, Z up) with the element colours by group
  (ELECOLOR 1, ElemPalette), the faces flat-shaded and outlined, beams and springs as lines in their group
  colour, the interaction nodes as red dots (the colour of the 3D view) and the lumped masses (MT) as dark
  diamonds when there are at most :data:`MAX_MASS_MARKERS` of them.  Faces, lines and markers are drawn
  in one depth-sorted collection (painter's algorithm), so the far side of the model stays hidden;
* **LAYERPLOT** -- the soil column (heights proportional to the thicknesses, SoilLayer palette), the
  half-space band, and the shear-wave velocity of each run of identical layers.

Rendering uses the matplotlib object API on the Agg canvas only (no ``pyplot``, no backend switch;
:func:`sassi.plotting.render_mpl._figure`).  The output is deterministic: no time stamp or software tag
in the PNG, a fixed colour palette (Pillow, at most :data:`COLORS` colours, no dithering).
"""
from __future__ import annotations

import argparse
import io
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..plotting.render_mpl import _figure
from ..plotting.state import palette_color, plot_data, plot_state, to_screen
from . import learn
from .lessons import EXAMPLES_DIR, copy_example

#: where the pictures are written (served as static/examples/<name>.png, learn.THUMB_URL)
OUT_DIR = Path(__file__).resolve().parent / "static" / "examples"
#: picture size in pixels; the cards show it at 240 x 150 CSS pixels (crisp at 2x density)
WIDTH, HEIGHT = 480, 300
DPI = 100
#: background: the preview band of the cards (styles.css .xthumb)
BG = "#f3f6fa"
MARGIN = 16                       # pixels around the model
INTERACTION = "#e00000"           # interaction nodes (the colour of the 3D model view)
MASS = "#1d2733"                  # lumped masses
MAX_MASS_MARKERS = 12             # more mass nodes (e.g. every slab node) would clutter the picture
LIGHT = np.array([-0.35, 0.55, 0.76])    # screen coordinates: from the upper left, toward the viewer
COLORS = 128                      # palette size of the PNG


class ThumbnailError(RuntimeError):
    """An example's model part failed, or it has nothing to show."""


# ======================================================================================
# the example's model and plot data
# ======================================================================================
def example_plot(name: str, examples_dir: Path = EXAMPLES_DIR) -> Dict[str, Any]:
    """Run the model part of example ``name`` and its view command in a temporary copy (deleted
    afterwards); returns the plot data of that plot (:func:`sassi.plotting.state.plot_data`, numpy arrays)."""
    from ..prep import Interpreter, Kind
    with tempfile.TemporaryDirectory(prefix="sassi-thumb-") as tmp:
        pre = copy_example(name, tmp, examples_dir)
        lines = learn.example_model_lines(pre.read_text(encoding="utf-8", errors="replace"))
        view = learn.example_view_command(lines)
        ui = Interpreter(cwd=tmp, seed=1)
        state = plot_state(ui)
        state.auto_render = False                       # no image file next to the model
        # SHOWMASS,1: the lumped masses are in the plot data (drawn when there are few of them)
        for ln in lines + [view] + (["SHOWMASS,1"] if view == "MODELPLOT" else []):
            if not ui.execute(ln):
                errs = ui.sink.texts(Kind.ERROR)
                raise ThumbnailError(f"{name}: {ln}: {errs[-1] if errs else 'failed'}")
        plot = state.active_plot
        if plot is None or plot.kind != view:
            raise ThumbnailError(f"{name}: {view} opened no plot")
        return plot_data(state, plot, ui, raw=True)


# ======================================================================================
# drawing
# ======================================================================================
def _rgb(color: str) -> np.ndarray:
    from matplotlib.colors import to_rgb
    return np.asarray(to_rgb(color), float)


def _shade(color: str, normal: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Flat shading of a face (two-sided) and its darker outline colour."""
    n = np.linalg.norm(normal)
    lam = abs(float(normal @ LIGHT)) / (n * np.linalg.norm(LIGHT)) if n > 0 else 0.5
    base = _rgb(color)
    face = np.clip(base * (0.70 + 0.34 * lam), 0.0, 1.0)
    return face, np.clip(base * 0.42, 0.0, 1.0)


def _fit(ax, pts: np.ndarray) -> float:
    """Frame the screen points ``pts`` (k, 2) with :data:`MARGIN`; returns model units per pixel."""
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    span = np.maximum(hi - lo, 1e-9)
    upp = max(span[0] / (WIDTH - 2 * MARGIN), span[1] / (HEIGHT - 2 * MARGIN))
    mid = 0.5 * (lo + hi)
    ax.set_xlim(mid[0] - upp * WIDTH / 2, mid[0] + upp * WIDTH / 2)
    ax.set_ylim(mid[1] - upp * HEIGHT / 2, mid[1] + upp * HEIGHT / 2)
    return upp


def _blank():
    fig = _figure(WIDTH / DPI, HEIGHT / DPI, DPI)
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()
    ax.set_facecolor(BG)
    return fig, ax


def shown_faces(sc: Dict[str, Any]) -> List[int]:
    """The faces to draw: the scene's ``face_boundary`` (plates, and the faces of solids on the surface of the
    solid volume; coincident solids -- an excavated-soil element and the backfill or structure element on
    the same nodes, the usual flexible-volume layout of example 6 -- count once, the later element's face
    is drawn), the same rule as the GUI's 3D view (:func:`sassi.plotting.state.model_scene`)."""
    return [int(r) for r in np.nonzero(np.asarray(sc["face_boundary"], dtype=bool))[0]]


def model_figure(d: Dict[str, Any]):
    """MODELPLOT picture: faces, beams and node markers of the scene, depth-sorted, no axes."""
    from matplotlib.collections import PolyCollection
    sc, cam = d["scene"], d["camera"]
    P = to_screen(np.asarray(sc["xyz"], float).reshape(-1, 3), cam["basis"], cam["center"])  # z: toward viewer
    if not len(P):
        raise ThumbnailError("the model has no element to draw")
    faces = np.asarray(sc["faces"], dtype=np.int64).reshape(-1, 4)
    face_elem = np.asarray(sc["face_elem"], dtype=np.int64)
    colors = list(sc["elem_color"])
    edges = np.asarray(sc["edges"], dtype=np.int64).reshape(-1, 2)
    edge_elem = np.asarray(sc["edge_elem"], dtype=np.int64)
    points = np.asarray(sc["points"], dtype=np.int64)
    visible = np.asarray(sc["node_visible"], dtype=bool)
    shown = shown_faces(sc)
    fig, ax = _blank()
    upp = _fit(ax, P[visible, :2] if visible.any() else P[:, :2])     # the GUI's fit: every scene node
    # a node marker or a line must be drawn after the faces it lies on: its depth is at least the
    # depth of every face at its nodes (a node on the far side stays behind the faces in front of it)
    front = np.full(len(P), -np.inf)
    items: List[Tuple[float, int, Any, Any, Any, float]] = []     # depth, order, polygon, face, edge, width
    etype = np.asarray(sc["elem_type"], dtype=np.int64)
    eps = 1e-6 * float(np.ptp(P[:, 2]) + 1.0)
    for r in shown:
        f = faces[r][faces[r] >= 0]
        Q = P[f]
        normal = np.cross(Q[2] - Q[0], Q[-1] - Q[1] if len(f) == 4 else Q[1] - Q[0])
        face, edge = _shade(colors[face_elem[r]], normal)
        # a plate on the face of a solid (a basement wall on the excavated soil) is drawn over that face
        depth = float(Q[:, 2].mean()) + (eps if etype[face_elem[r]] != 1 else 0.0)
        front[f] = np.maximum(front[f], depth)
        items.append((depth, 0, Q[:, :2], face, edge, 0.55))
    w = 1.25 * upp                                                  # half width of a line, 2.5 px
    for k, (a, b) in enumerate(edges):
        u = P[b, :2] - P[a, :2]
        L = float(np.hypot(*u))
        if L < 1e-12:
            continue
        nrm = np.array([-u[1], u[0]]) / L * w
        quad = np.array([P[a, :2] + nrm, P[b, :2] + nrm, P[b, :2] - nrm, P[a, :2] - nrm])
        depth = max(float(P[[a, b], 2].mean()), float(max(front[a], front[b])) + 2 * eps)
        col = _rgb(colors[edge_elem[k]])
        items.append((depth, 1, quad, col, np.clip(col * 0.55, 0, 1), 0.4))

    def marker(i: int, radius: float, sides: int, rot: float) -> np.ndarray:
        t = rot + np.arange(sides) * 2 * np.pi / sides
        return P[i, :2] + radius * upp * np.c_[np.cos(t), np.sin(t)]

    def above(i: int) -> float:
        return max(float(P[i, 2]), float(front[i])) + 3 * eps

    for i in points:
        items.append((above(i), 2, marker(i, 3.0, 4, np.pi / 4), _rgb(colors[0]), _rgb("#202020"), 0.5))
    for i in np.asarray(sc["interaction"], dtype=np.int64):
        if visible[i]:
            items.append((above(i), 3, marker(i, 2.4, 16, 0.0), _rgb(INTERACTION), _rgb("#ffffff"), 0.5))
    mass = np.asarray(sc.get("mass", []), dtype=np.int64)
    if 0 < len(mass) <= MAX_MASS_MARKERS:
        for i in mass:
            items.append((above(i), 4, marker(i, 5.0, 4, np.pi / 2), _rgb(MASS), _rgb("#ffffff"), 0.7))
    items.sort(key=lambda it: (it[0], it[1]))
    ax.add_collection(PolyCollection([it[2] for it in items], facecolors=[it[3] for it in items],
                                     edgecolors=[it[4] for it in items], linewidths=[it[5] for it in items],
                                     joinstyle="round", antialiased=True))
    return fig


#: soil tones of the layer picture, softest (lowest Vs) first
SOIL_TONES = ("#efe3c4", "#dcc597", "#c4a771", "#a98a57", "#8c6f43", "#71582f")
INK, MUTED = "#1c222a", "#5d6774"


def layer_figure(d: Dict[str, Any]):
    """LAYERPLOT picture: the soil column (heights proportional to the thicknesses, a thin line between
    layers), the half-space, and the shear-wave velocity of each run of identical layers (darker = stiffer)."""
    from matplotlib.patches import Rectangle
    t = d["table"]
    rows = t["layers"]
    if not rows:
        raise ThumbnailError("no soil layer to draw")
    hs = t.get("halfspace")
    fig, ax = _blank()
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(HEIGHT, 0)                                  # pixel coordinates, y down
    total = float(sum(r["thick"] for r in rows)) or 1.0
    y0, y1 = float(MARGIN), float(HEIGHT - MARGIN)
    hs_px = 0.17 * (y1 - y0) if hs else 0.0
    ppu = (y1 - y0 - hs_px) / total                         # pixels per length unit
    x0, cw = 124.0, 104.0
    # runs of the same L layer (TOPL repeats a layer to subdivide it)
    runs: List[List[Any]] = []
    for r in rows:
        if runs and runs[-1][2]["layer"] == r["layer"]:
            runs[-1][1] += r["thick"]
        else:
            runs.append([sum(x[1] - x[0] for x in runs) if runs else 0.0, None, r])
            runs[-1][1] = runs[-1][0] + r["thick"]
    ranks = sorted({float(r["vs"]) for _, _, r in runs})
    step = min(2.0, (len(SOIL_TONES) - 1) / max(1, len(ranks) - 1))
    tone = {v: SOIL_TONES[int(round(k * step))] for k, v in enumerate(ranks)}
    y = y0
    for r in rows:
        h = r["thick"] * ppu
        c = _rgb(tone[float(r["vs"])])
        ax.add_patch(Rectangle((x0, y), cw, h, facecolor=c, edgecolor=np.clip(c * 0.82, 0, 1), linewidth=0.5))
        y += h
    for top, bot, r in runs:
        ax.add_patch(Rectangle((x0, y0 + top * ppu), cw, (bot - top) * ppu, facecolor="none", edgecolor=INK,
                               linewidth=0.9))
        ax.text(x0 + cw + 14, y0 + 0.5 * (top + bot) * ppu, f"Vs {r['vs']:g}", ha="left", va="center",
                fontsize=15, color=INK)
    if hs:
        ax.add_patch(Rectangle((x0, y1 - hs_px), cw, hs_px, facecolor="#dfe4ea", edgecolor=INK, hatch="////",
                               linewidth=0.9))
        ax.text(x0 + cw + 14, y1 - hs_px / 2, f"Vs {hs['vs']:g}  half-space", ha="left", va="center",
                fontsize=15, color=INK)
    ax.text(x0 - 10, y0, "0", ha="right", va="center", fontsize=12, color=MUTED)
    ax.text(x0 - 10, y1 - hs_px, f"{total:g}", ha="right", va="center", fontsize=12, color=MUTED)
    ax.text(x0 - 10, 0.5 * (y0 + y1 - hs_px), "depth", ha="right", va="center", fontsize=12, color=MUTED,
            rotation=90)
    return fig


def figure(d: Dict[str, Any]):
    """The picture of a plot-data dictionary (MODELPLOT or LAYERPLOT)."""
    if d["family"] == "layer":
        return layer_figure(d)
    if d["kind"] == "MODELPLOT":
        return model_figure(d)
    raise ThumbnailError(f"no picture for {d['kind']}")


def png_bytes(fig) -> bytes:
    """The figure as a small, deterministic PNG (palette of at most :data:`COLORS` colours, no metadata)."""
    from PIL import Image
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=DPI, facecolor=fig.get_facecolor(), metadata={"Software": None})
    buf.seek(0)
    im = Image.open(buf).convert("RGB")
    pal = im.quantize(colors=COLORS, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    out = io.BytesIO()
    pal.save(out, format="PNG", optimize=True)
    return out.getvalue()


def thumbnail(name: str, examples_dir: Path = EXAMPLES_DIR) -> bytes:
    """The PNG picture of example ``name``."""
    return png_bytes(figure(example_plot(name, examples_dir)))


def example_names(examples_dir: Path = EXAMPLES_DIR) -> List[str]:
    return [p.stem for p in sorted(Path(examples_dir).glob("*.pre")) if learn.EXAMPLE_RE.match(p.stem)]


def generate(names: Optional[Sequence[str]] = None, out: Path = OUT_DIR,
             examples_dir: Path = EXAMPLES_DIR) -> Dict[str, int]:
    """Write ``<out>/<name>.png`` for the examples ``names`` (default: all); returns {file name: bytes}."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    sizes: Dict[str, int] = {}
    for name in names or example_names(examples_dir):
        data = thumbnail(name, examples_dir)
        (out / f"{name}.png").write_bytes(data)
        sizes[f"{name}.png"] = len(data)
    return sizes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m sassi.ui.thumbnails",
                                 description="write the preview pictures of the examples gallery")
    ap.add_argument("names", nargs="*", help="example names (default: every examples/*.pre)")
    ap.add_argument("--out", default=str(OUT_DIR), help=f"output folder (default {OUT_DIR})")
    args = ap.parse_args(argv)
    known = example_names()
    bad = [n for n in args.names if n not in known]
    if bad:
        print(f"unknown example(s): {', '.join(bad)} (known: {', '.join(known)})", file=sys.stderr)
        return 2
    try:
        sizes = generate(args.names or None, Path(args.out))
    except ThumbnailError as exc:
        print(f"thumbnail failed: {exc}", file=sys.stderr)
        return 1
    for name, n in sizes.items():
        print(f"  {name:32s} {n / 1e3:6.1f} kB")
    print(f"{len(sizes)} picture(s), {sum(sizes.values()) / 1e3:.1f} kB, {WIDTH} x {HEIGHT} px, in {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
