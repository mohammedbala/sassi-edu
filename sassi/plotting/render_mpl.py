"""Headless matplotlib rendering of the plot data (requirements UI-03, 5.7; D-UI-06, D-UI-10).

Every function draws from the dictionaries of :func:`sassi.plotting.state.plot_data`, so the image
written in batch mode and the GUI (plotly) show the same content.  Only the object-oriented
matplotlib API with the Agg canvas is used (no ``pyplot``, no global backend switch), so rendering
works without a display and does not interfere with a GUI process.

* SPECPLOT / THPLOT: 2D lines (colours by line number, legend, AXES grids and log scales,
  PLOTRANGE, titles, STIPPLE dashes, MARKERS); tick labels below 1e-12 of the axis range print 0.
* LAYERPLOT: soil column with heights proportional to the thicknesses, half-space band, and the
  property table (Unit Weight, Vp, Vs, Dp, Ds; Thickness optional); PLOTTITLE is ignored.  Layers
  outside Start / End Layer are marked by a dashed "not shown" band and table row.
* SOILPROPPLOT: damping (left axis, red) and G/Gmax (right axis, green) vs shear strain %
  (log axis by default).
* MODELPLOT / NODEPLOT / CUTPLOT and the animations: mplot3d poly collections.  The model is
  rotated by the CNGVIEW matrix ``R`` (:func:`sassi.plotting.state.rotation_matrix`) and viewed
  from the top with an orthographic camera, so the screen coordinates of a node are exactly
  ``R (p - centre)``; mplot3d sorts the faces by depth (painter's algorithm).
"""
from __future__ import annotations

import io
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .state import (ANIMATIONS, PlotState, colorbar_labels, jet, palette_color, plot_data, to_screen)

PathLike = Union[str, Path]
DPI = 100
#: transparent fill (mplot3d collections fail to draw with the string 'none' in matplotlib 3.9)
NOFILL = (0.0, 0.0, 0.0, 0.0)


def _txt(s: Any) -> str:
    """User text for matplotlib: a '$' is literal (no mathtext), so any title or name renders."""
    return str(s or "").replace("$", r"\$")


def _figure(w: float, h: float, dpi: int = DPI):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fig = Figure(figsize=(w, h), dpi=dpi)
    FigureCanvasAgg(fig)
    return fig


def _clean_formatter(ax_axis):
    """Print values smaller than 1e-12 of the axis span as 0 (spec 06 section 5.7: no -2.7756e-17)."""
    from matplotlib.ticker import FuncFormatter, ScalarFormatter
    base = ScalarFormatter(useOffset=False)
    base.set_axis(ax_axis)

    def fmt(v, pos=None):
        lo, hi = ax_axis.get_view_interval()
        span = abs(hi - lo) or 1.0
        if abs(v) < 1e-12 * span:
            v = 0.0
        return f"{v:.6g}"
    return FuncFormatter(fmt)


# ======================================================================================
# 2D line plots
# ======================================================================================
def _apply_axes(ax, s: Dict[str, Any], extent: Sequence[Optional[float]]) -> None:
    if s.get("log_x"):
        ax.set_xscale("log")
    else:
        ax.xaxis.set_major_formatter(_clean_formatter(ax.xaxis))
    if s.get("log_y"):
        ax.set_yscale("log")
    else:
        ax.yaxis.set_major_formatter(_clean_formatter(ax.yaxis))
    x0, x1, y0, y1 = extent
    if x0 is not None and x1 is not None and x1 > x0:
        ax.set_xlim(x0, x1)
    elif x0 is not None and x1 is not None:
        pad = abs(x0) * 0.05 or 1.0
        ax.set_xlim(x0 - pad, x1 + pad)
    if y0 is not None and y1 is not None and y1 > y0:
        pad = 0.0 if s.get("ymin") is not None or s.get("ymax") is not None else 0.04 * (y1 - y0)
        if s.get("log_y"):
            ax.set_ylim(y0, y1)
        else:
            ax.set_ylim(y0 - pad if s.get("ymin") is None else y0, y1 + pad if s.get("ymax") is None else y1)
    elif y0 is not None and y1 is not None:
        pad = abs(y0) * 0.05 or 1.0
        ax.set_ylim(y0 - pad, y1 + pad)
    if s.get("minor_x") or s.get("minor_y"):
        ax.minorticks_on()
    ax.grid(bool(s.get("major_x")), which="major", axis="x", color="#a0a0a0", linewidth=1.0)
    ax.grid(bool(s.get("major_y")), which="major", axis="y", color="#a0a0a0", linewidth=1.0)
    if s.get("minor_x"):
        ax.grid(True, which="minor", axis="x", color="#d0d0d0", linewidth=0.5)
    if s.get("minor_y"):
        ax.grid(True, which="minor", axis="y", color="#d0d0d0", linewidth=0.5)


def line_figure(d: Dict[str, Any], dpi: int = DPI):
    """SPECPLOT / THPLOT figure."""
    fig = _figure(8.5, 5.2, dpi)
    ax = fig.add_subplot(111)
    s = d["settings"]
    for L in d["lines"]:
        ax.plot(np.asarray(L["x"], float), np.asarray(L["y"], float), color=L["color"],
                linestyle=L.get("dash", "solid"), linewidth=1.4, marker="o" if L.get("markers") else None,
                markersize=3.5, label=_txt(f"{L['name'] or 'line'} ({L['number']})"))
    _apply_axes(ax, s, d["extent"])
    if d["title"]:
        ax.set_title(_txt(d["title"]))
    ax.set_xlabel(_txt(s.get("xtitle", "")))
    ax.set_ylabel(_txt(s.get("ytitle", "")))
    if s.get("ytitle2"):
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yscale(ax.get_yscale())
        ax2.set_ylabel(_txt(s["ytitle2"]))
    if d["lines"]:
        ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    return fig


# ======================================================================================
# Soil layer plot and soil property plot
# ======================================================================================
def _fmt(v: Optional[float]) -> str:
    if v is None:
        return "-"
    return f"{v:.6g}"


def layer_figure(d: Dict[str, Any], dpi: int = DPI):
    """LAYERPLOT figure: layer column (heights proportional to thickness) and property table."""
    from matplotlib.patches import Rectangle
    t = d["table"]
    show = {"thick": False, "weight": True, "vp": True, "vs": True, "pdamp": True, "sdamp": True}
    show.update(d.get("show") or {})
    rows = t["layers"]
    above, below = int(t.get("omitted_above", 0)), int(t.get("omitted_below", 0))
    nrows = len(rows) + 1 + (above > 0) + (below > 0)
    fig = _figure(11.0, max(4.5, 0.32 * nrows + 2.0), dpi)
    axc = fig.add_axes([0.04, 0.06, 0.2, 0.84])
    total = sum(r["thick"] for r in rows) or 1.0
    hs_h = 0.12 * total
    gap_h = 0.08 * total

    def gap(y0: float, text: str) -> float:
        # layers of the profile that the plot does not show (Start / End Layer): a dashed band, so
        # the half-space is not drawn as if it began right under the last displayed layer
        axc.add_patch(Rectangle((0, -(y0 + gap_h)), 1, gap_h, facecolor="none", edgecolor="#808080",
                                linestyle="--", linewidth=0.8))
        axc.text(0.5, -(y0 + gap_h / 2), text, ha="center", va="center", fontsize=7, color="#404040")
        return y0 + gap_h

    y = 0.0
    if above:
        y = gap(y, f"layers 1-{above} not shown" if above > 1 else "layer 1 not shown")
    for k, r in enumerate(rows):
        h = r["thick"]
        axc.add_patch(Rectangle((0, -(y + h)), 1, h, facecolor=palette_color("SOILLAYER", k + 1),
                                edgecolor="#404040", linewidth=0.8))
        if h > 0.012 * total:
            axc.text(-0.05, -(y + h / 2), str(r["index"]), ha="right", va="center", fontsize=8)
        y += h
    if below:
        first, last = t["end"] + 1, t["end"] + below
        y = gap(y, f"layers {first}-{last} not shown" if below > 1 else f"layer {first} not shown")
    axc.add_patch(Rectangle((0, -(y + hs_h)), 1, hs_h, facecolor="#f4f4f4", edgecolor="#404040", hatch="//",
                            linewidth=0.8))
    hsr = t.get("halfspace")
    hs_label = "Halfspace" + (f"\n(depth {_fmt(hsr['top'])})" if hsr and (above or below) else "")
    axc.text(0.5, -(y + hs_h / 2), hs_label, ha="center", va="center", fontsize=8,
             bbox=dict(facecolor="#f4f4f4", edgecolor="none", pad=1.0))
    axc.set_xlim(-0.35, 1.05)
    axc.set_ylim(-(y + hs_h) * 1.01, 0.01 * total)
    axc.set_axis_off()
    cols = [("Layer", None)]
    labels = {"thick": "Thickness", "weight": "Unit Weight", "vp": "P-Wave Velocity", "vs": "S-Wave Velocity",
              "pdamp": "P-Wave Damping Ratio", "sdamp": "S-Wave Damping Ratio"}
    for key in ("thick", "weight", "vp", "vs", "pdamp", "sdamp"):
        if show.get(key):
            cols.append((labels[key], key))
    cells = []
    if above:
        cells.append([f"1-{above} (not shown)" if above > 1 else "1 (not shown)"] + ["..."] * (len(cols) - 1))
    for r in rows:
        cells.append([f"{r['index']} (L {r['layer']})"] + [_fmt(r[k]) for _, k in cols[1:]])
    if below:
        first, last = t["end"] + 1, t["end"] + below
        cells.append([f"{first}-{last} (not shown)" if below > 1 else f"{first} (not shown)"]
                     + ["..."] * (len(cols) - 1))
    cells.append(["Halfspace" + (f" (L {hsr['layer']})" if hsr else "")] +
                 [(_fmt(hsr[k]) if hsr and k != "thick" else "-") for _, k in cols[1:]])
    axt = fig.add_axes([0.27, 0.06, 0.71, 0.84])
    axt.set_axis_off()
    tab = axt.table(cellText=[[_txt(v) for v in row] for row in cells], colLabels=[c for c, _ in cols],
                    loc="upper center", cellLoc="center")
    tab.auto_set_font_size(False)
    tab.set_fontsize(8)
    tab.auto_set_column_width(list(range(len(cols))))
    tab.scale(1.0, 1.25)
    fig.suptitle("Soil Layer Plot", fontsize=11)
    return fig


def soilprop_figure(d: Dict[str, Any], dpi: int = DPI):
    """SOILPROPPLOT figure: damping (left, red) and G/Gmax (right, green) vs shear strain %."""
    c = d["curves"]
    s = d["settings"]
    show = d["params"].get("show", {}) or {}
    fig = _figure(8.5, 5.2, dpi)
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    hd = hg = None
    if show.get("damping", True) and c["d_strain"]:
        hd, = ax.plot(c["d_strain"], c["d"], color=palette_color("SOILPROP", 3), marker="o", markersize=4,
                      label="Damping Ratio/Shear Strain")
    if show.get("modulus", True) and c["g_strain"]:
        hg, = ax2.plot(c["g_strain"], c["g"], color=palette_color("SOILPROP", 4), marker="o", markersize=4,
                       label="Shear Modulus/Shear Strain")
    if s.get("log_x"):
        ax.set_xscale("log")
    if s.get("log_y"):
        ax.set_yscale("log")
    if s.get("xmin") is not None or s.get("xmax") is not None:
        ax.set_xlim(s.get("xmin"), s.get("xmax"))          # None keeps the automatic bound
    if s.get("ymin") is not None or s.get("ymax") is not None:
        ax.set_ylim(s.get("ymin"), s.get("ymax"))
    ax2.set_ylim(0.0, 1.05)
    ax.grid(bool(s.get("major_x")), which="major", axis="x", color="#c0c0c0")
    ax.grid(bool(s.get("major_y")), which="major", axis="y", color="#c0c0c0")
    if s.get("minor_x"):
        ax.grid(True, which="minor", axis="x", color="#e0e0e0", linewidth=0.5)
    ax.set_xlabel(_txt(s.get("xtitle") or "Shear Strain %"))
    ax.set_ylabel(_txt(s.get("ytitle") or "Damping Ratio (%)"), color=palette_color("SOILPROP", 3))
    ax2.set_ylabel(_txt(s.get("ytitle2") or "Shear Modulus (G/Gmax)"), color=palette_color("SOILPROP", 4))
    hs = [h for h in (hd, hg) if h is not None]
    if hs:
        ax.legend(hs, [h.get_label() for h in hs], loc="center left", fontsize=8)
    ax.set_title(_txt(d["title"] or c["name"]))
    fig.tight_layout()
    return fig


# ======================================================================================
# 3D plots
# ======================================================================================
def _triad(fig, R: np.ndarray) -> None:
    """Axis triad in the lower-left corner: X red, Y green, Z blue (spec 06 section 1.1)."""
    ax = fig.add_axes([0.01, 0.01, 0.12, 0.12])
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal")
    ax.set_axis_off()
    for k, (lab, col) in enumerate((("X", "#e00000"), ("Y", "#00a000"), ("Z", "#0000e0"))):
        v = R[:2, k]
        ax.plot([0, v[0]], [0, v[1]], color=col, linewidth=2)
        ax.text(1.15 * v[0], 1.15 * v[1], lab, color=col, fontsize=8, ha="center", va="center")


def _colorbar(fig, vmin: float, vmax: float) -> None:
    """Vertical jet colour bar on the right with five labels from Max (top) to Min (spec 06 8.7)."""
    ax = fig.add_axes([0.9, 0.2, 0.025, 0.6])
    grad = jet(np.linspace(0, 1, 256))[:, None, :]
    ax.imshow(grad, aspect="auto", origin="lower", extent=(0, 1, 0, 1))
    ax.set_xticks([])
    ax.yaxis.tick_right()
    ax.set_yticks(np.linspace(1, 0, 5))
    ax.set_yticklabels(colorbar_labels(vmin, vmax))
    ax.tick_params(labelsize=7)


def _screen(xyz: np.ndarray, cam: Dict[str, Any]) -> np.ndarray:
    return to_screen(xyz, cam["basis"], cam["center"])


def model_figure(d: Dict[str, Any], dpi: int = DPI):
    """MODELPLOT / NODEPLOT / CUTPLOT and animation figures (mplot3d, orthographic)."""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the projection)
    from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
    sc = d["scene"]
    cam = d["camera"]
    v = d["view"]
    kind = d["kind"]
    shader = d.get("shader", {})
    fig = _figure(8.5, 7.5, dpi)
    ax = fig.add_axes([0.0, 0.0, 0.88 if kind in ("BUBBLEPLOT", "CONTOURPLOT") else 1.0, 0.94], projection="3d")
    ax.set_proj_type("ortho")
    ax.view_init(elev=90, azim=-90)
    ax.set_axis_off()
    xyz = np.asarray(sc["xyz"], float).reshape(-1, 3)
    frame = d.get("frame") or {}
    if kind == "DEFORMPLOT" and "xyz" in frame:
        xyz_draw = np.asarray(frame["xyz"], float).reshape(-1, 3)
    else:
        xyz_draw = xyz
    P = _screen(xyz_draw, cam)
    P0 = _screen(xyz, cam)
    faces = np.asarray(sc["faces"], dtype=np.int64).reshape(-1, 4)
    face_elem = np.asarray(sc["face_elem"], dtype=np.int64)
    boundary = np.asarray(sc["face_boundary"], dtype=bool)
    ecolor = list(sc["elem_color"])
    cent = _screen(np.asarray(sc["elem_centroid"], float), cam) if len(sc["elem_centroid"]) else np.zeros((0, 3))
    shrink = float(shader.get("shrink", 0.06)) if v.get("shrink") and kind == "MODELPLOT" else 0.0
    wire = bool(v.get("wireframe")) and kind == "MODELPLOT"
    show_all = shrink > 0 or wire
    sel = np.nonzero(boundary | show_all)[0] if len(faces) else np.zeros(0, dtype=np.int64)

    def polys_of(idx, pts):
        out = []
        for r in idx:
            f = faces[r][faces[r] >= 0]
            poly = pts[f]
            if shrink > 0:
                c = cent[face_elem[r]] if kind == "MODELPLOT" else poly.mean(axis=0)
                poly = c + (1.0 - shrink) * (poly - c)
            out.append(poly)
        return out

    outline = palette_color("ELEMENT", 2)
    if kind == "NODEPLOT":
        pass
    elif kind == "CUTPLOT":
        incut = np.asarray(sc["elem_in_cut"], dtype=bool)
        cutf = [r for r in range(len(faces)) if incut[face_elem[r]]]
        allf = list(sel)
        if allf:
            ax.add_collection3d(Poly3DCollection(polys_of(allf, P), facecolors=NOFILL,
                                                 edgecolors=palette_color("CUT", 2), linewidths=0.4))
        if cutf:
            ax.add_collection3d(Poly3DCollection(polys_of(cutf, P), facecolors=palette_color("CUT", 3),
                                                 edgecolors=outline, linewidths=0.6))
    elif kind == "CONTOURPLOT" and len(sel):
        t = np.asarray(frame.get("t", []), float)
        fc = []
        for r in sel:
            f = faces[r][faces[r] >= 0]
            tv = t[f] if len(t) else np.array([np.nan])
            fc.append(jet(np.nanmean(tv))[()].tolist() if np.any(np.isfinite(tv)) else [0.75, 0.75, 0.75])
        ax.add_collection3d(Poly3DCollection(polys_of(sel, P), facecolors=fc, edgecolors=outline, linewidths=0.4))
    elif kind == "DEFORMPLOT" and len(sel):
        if v.get("show_undeformed"):
            ax.add_collection3d(Poly3DCollection(polys_of(sel, P0), facecolors=NOFILL,
                                                 edgecolors=palette_color("DEFORMED", 4), linewidths=0.4))
        ax.add_collection3d(Poly3DCollection(polys_of(sel, P), facecolors=palette_color("DEFORMED", 3),
                                             edgecolors=palette_color("DEFORMED", 2), linewidths=0.5))
    elif kind in ("BUBBLEPLOT", "VECTORPLOT") and len(sel):
        ax.add_collection3d(Poly3DCollection(polys_of(sel, P), facecolors=NOFILL, edgecolors="#c8c8c8",
                                             linewidths=0.3))
    elif len(sel):
        fc = NOFILL if wire else [ecolor[face_elem[r]] for r in sel]
        ec = [ecolor[face_elem[r]] for r in sel] if wire else outline
        ax.add_collection3d(Poly3DCollection(polys_of(sel, P), facecolors=fc, edgecolors=ec,
                                             linewidths=0.8 if wire else 0.4))
    # beams, springs, GENERAL elements as lines
    edges = np.asarray(sc["edges"], dtype=np.int64).reshape(-1, 2)
    if len(edges) and kind != "NODEPLOT":
        segs = [P[e] for e in edges]
        cols = [ecolor[k] for k in np.asarray(sc["edge_elem"], dtype=np.int64)]
        if kind == "CUTPLOT":
            incut = np.asarray(sc["elem_in_cut"], dtype=bool)
            cols = [palette_color("CUT", 3) if incut[k] else palette_color("CUT", 2)
                    for k in np.asarray(sc["edge_elem"], dtype=np.int64)]
        ax.add_collection3d(Line3DCollection(segs, colors=cols, linewidths=2.0))
    pts = np.asarray(sc["points"], dtype=np.int64)
    if len(pts) and kind != "NODEPLOT":
        ax.scatter(P[pts, 0], P[pts, 1], P[pts, 2], marker="x", color=palette_color("ELEMENT", 5), s=30)
    half = float(cam["half"])
    msize = float(shader.get("points", 10.0))
    vis = np.asarray(sc.get("node_visible", np.ones(len(P), bool)), dtype=bool)
    # nodes and markers
    if kind == "NODEPLOT" and len(P):
        inter = np.zeros(len(P), dtype=bool)
        inter[np.asarray(sc["interaction"], dtype=np.int64)] = True
        ordinary = vis & ~inter
        ax.scatter(P[ordinary, 0], P[ordinary, 1], P[ordinary, 2], s=msize * 1.2, color=palette_color("NODE", 1),
                   depthshade=False)
        ii = vis & inter
        ax.scatter(P[ii, 0], P[ii, 1], P[ii, 2], s=msize * 1.6, color=palette_color("NODE", 2), depthshade=False)
        fx = np.asarray(sc["fixed"], dtype=np.int64)
        if len(fx):
            ax.scatter(P[fx, 0], P[fx, 1], P[fx, 2], s=msize * 6, facecolors=NOFILL,
                       edgecolors=palette_color("NODE", 3), linewidths=1.2, depthshade=False)
        ms = np.asarray(sc["mass"], dtype=np.int64)
        if len(ms):
            ax.scatter(P[ms, 0], P[ms, 1], P[ms, 2], s=msize * 10, facecolors=NOFILL,
                       edgecolors=palette_color("NODE", 4), linewidths=1.2, depthshade=False)
    elif kind in ("MODELPLOT", "CUTPLOT") and len(P):
        ii = np.asarray(sc["interaction"], dtype=np.int64)
        if len(ii) and kind == "MODELPLOT":
            ax.scatter(P[ii, 0], P[ii, 1], P[ii, 2], s=msize * 0.8, color=palette_color("ELEMENT", 10),
                       depthshade=False)
        fx = np.asarray(sc["fixed"], dtype=np.int64)
        if len(fx):
            ax.scatter(P[fx, 0], P[fx, 1], P[fx, 2], s=msize * 4, marker="^", color=palette_color("ELEMENT", 9),
                       depthshade=False)
        ms = np.asarray(sc["mass"], dtype=np.int64)
        if len(ms):
            md = np.asarray(sc["mass_dir"], dtype=bool).reshape(-1, 3)
            R = np.asarray(cam["basis"], float)
            L = 0.04 * half
            for j, col in enumerate((palette_color("ELEMENT", 6), palette_color("ELEMENT", 7),
                                     palette_color("ELEMENT", 8))):
                rows = ms[md[:, j]]
                if len(rows):
                    dvec = R[:, j] * L
                    segs = [np.vstack([P[r], P[r] + dvec]) for r in rows]
                    ax.add_collection3d(Line3DCollection(segs, colors=col, linewidths=2.5))
    if kind == "BUBBLEPLOT" and len(P):
        rgb = np.asarray(frame.get("rgb", []), float).reshape(-1, 3)
        size = np.asarray(frame.get("size", []), float)
        found = np.asarray(frame.get("found", []), dtype=bool)
        if len(rgb):
            m = found & vis
            ax.scatter(P[m, 0], P[m, 1], P[m, 2], s=(size[m] * 1.5) ** 2, c=rgb[m], edgecolors="#202020",
                       linewidths=0.3, depthshade=False)
    if kind == "VECTORPLOT" and len(P):
        V = np.asarray(frame.get("vectors", []), float).reshape(-1, 3, 3)
        found = np.asarray(frame.get("found", []), dtype=bool)
        R = np.asarray(cam["basis"], float)
        ax.scatter(P[vis, 0], P[vis, 1], P[vis, 2], s=6, color="#000000", depthshade=False)
        dirs = {"X": [0], "Y": [1], "Z": [2]}.get(str(frame.get("direction", "X")).upper(), [0, 1, 2])
        for dd in range(3):                     # X red, Y green, Z blue (Vector palette 3..5)
            col = palette_color("VECTOR", 3 + dd)
            if dd not in dirs or not len(V):
                continue
            Vs = to_screen(V[:, dd, :], R, (0.0, 0.0, 0.0))
            segs = [np.vstack([P[i], P[i] + Vs[i]]) for i in range(len(P)) if found[i] and vis[i]]
            if segs:
                ax.add_collection3d(Line3DCollection(segs, colors=col, linewidths=1.6))
    sel_n = np.asarray(sc["selected"], dtype=np.int64)
    if len(sel_n):
        ax.scatter(P[sel_n, 0], P[sel_n, 1], P[sel_n, 2], s=msize * 9, marker="s", facecolors=NOFILL,
                   edgecolors=palette_color("NODE", 5), linewidths=1.5, depthshade=False)
    # labels
    if v.get("node_labels") and len(P):
        ids = np.asarray(sc["node_id"], dtype=np.int64)
        for i in np.nonzero(vis)[0][:5000]:
            ax.text(P[i, 0], P[i, 1], P[i, 2], f" {ids[i]}", fontsize=6, color="#000080")
    if (v.get("elem_labels") or v.get("group_labels")) and len(cent) and kind != "NODEPLOT":
        lab = np.asarray(sc["elem_id"] if v.get("elem_labels") else sc["elem_group"], dtype=np.int64)
        for k in range(min(len(cent), 5000)):
            ax.text(cent[k, 0], cent[k, 1], cent[k, 2], str(lab[k]), fontsize=6, color="#800000",
                    ha="center", va="center")
    x0, x1, y0, y1 = cam["extent"]
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    if len(P):
        z0, z1 = float(P[:, 2].min()), float(P[:, 2].max())
    else:
        z0, z1 = -1.0, 1.0
    if z1 - z0 < 1e-9 * max(1.0, x1 - x0):
        z0, z1 = z0 - 0.5, z1 + 0.5
    ax.set_zlim(z0, z1)
    ax.set_box_aspect((1.0, (y1 - y0) / (x1 - x0), max((z1 - z0) / (x1 - x0), 0.01)), zoom=1.3)
    _triad(fig, np.asarray(cam["basis"], float))
    title = d.get("title") or d.get("caption", "")
    if kind in ANIMATIONS and frame.get("label"):
        title = f"{title}   [{frame['label']}]"
    fig.suptitle(_txt(title), fontsize=10)
    if kind in ("BUBBLEPLOT", "CONTOURPLOT"):
        a = d.get("animation", {})
        _colorbar(fig, float(a.get("vmin", 0.0)), float(a.get("vmax", 1.0)))
        if kind == "CONTOURPLOT":
            fig.text(0.5, 0.955, "nodal averages for plotting only; design values are element-centre values",
                     ha="center", fontsize=7, color="#606060")
    if v.get("debug"):
        fig.text(0.01, 0.97, f"rX {v['rx']:g}  rY {v['ry']:g}  rZ {v['rz']:g}  px {v['px']:g}  py {v['py']:g}  "
                             f"zoom {v['zoom']:g}", fontsize=7, family="monospace")
    return fig


# ======================================================================================
# Entry points
# ======================================================================================
def figure(d: Dict[str, Any], dpi: int = DPI):
    """Matplotlib figure of a plot-data dictionary (:func:`sassi.plotting.state.plot_data`)."""
    fam = d["family"]
    if fam == "2d":
        return line_figure(d, dpi)
    if fam == "layer":
        return layer_figure(d, dpi)
    if fam == "soilprop":
        return soilprop_figure(d, dpi)
    return model_figure(d, dpi)


def image_format(file_name: PathLike) -> str:
    """CAPTUREPLOT rule (D-UI-10): PNG when ``.png`` appears anywhere in ``<FileName>`` (any case),
    else BMP.  Pass the user's FileName text, not the resolved path: the model or working
    directory (e.g. ``plots.png.d``) must not decide the format."""
    return "png" if ".png" in str(file_name).lower() else "bmp"


def save_figure(fig, path: PathLike, fmt: str = "png") -> Tuple[Path, List[str]]:
    """Write ``fig`` as PNG or BMP; BMP needs Pillow (else PNG is written with a note)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    notes: List[str] = []
    if fmt == "bmp":
        try:
            from PIL import Image
        except ImportError:  # pragma: no cover - Pillow ships with matplotlib
            notes.append("Pillow is not installed: the image was written in PNG format")
            fig.savefig(str(path), format="png")
            return path, notes
        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        buf.seek(0)
        Image.open(buf).convert("RGB").save(str(path), format="BMP")
        return path, notes
    fig.savefig(str(path), format="png")
    return path, notes


def render_plot(state: PlotState, plot, interp, path: PathLike, fmt: Optional[str] = None,
                dpi: Optional[int] = None) -> Tuple[Path, List[str]]:
    """Render ``plot`` (a :class:`~sassi.plotting.state.Plot` or its id) to ``path``."""
    p = state.plots[plot] if isinstance(plot, int) else plot
    d = plot_data(state, p, interp, raw=True)
    fig = figure(d, dpi or state.image_dpi)
    return save_figure(fig, path, fmt or "png")
