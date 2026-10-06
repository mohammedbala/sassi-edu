"""The soil island picture (SHOWSOIL; requirements 7.20, D-W6-01 to D-W6-04).

SASSI has no soil island: the layered site is horizontally infinite.  SHOWSOIL draws the free-field layers
around the foundation as a display aid.  Checked here: the geometry (layers at their depths, the opening the
model fills, the cut-away quarter, outward faces, extents), the command (toggle, fields, errors, defaults,
message), the plot data (box and camera fit), both renderers' inputs and the lesson action."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sassi.plotting import state as S
from sassi.plotting.render_mpl import _ribbons, _soil_hidden
from sassi.plotting.state import PlotError, plot_data, plot_state, soil_island
from sassi.prep import Interpreter, Kind

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"

#: a 4 x 4 x 2 m excavation (two 1 m embedment layers) in a profile of 2 + 3 layers on a half-space
SITE = """L,1,1.0,19,600,200,0.05,0.05
L,2,1.0,19,700,250,0.05,0.05
L,3,2.0,20,900,400,0.04,0.04
L,4,0,22,2000,1000,0.02,0.02
TOPL,1,2,3,3,3
SITE,0,1,0,10,4,1,0,1,2048,1,0,0.005,4096,1
"""


def _box(depth_levels):
    """Nodes of a 4 x 4 m plan, 3 x 3 nodes per level (2 m spacing), at the given z levels; excavated
    SOLID elements between the levels; interaction nodes on all of them (FV)."""
    lines, n = [], 0
    ids = {}
    for k, z in enumerate(depth_levels):
        for j in range(3):
            for i in range(3):
                n += 1
                ids[(i, j, k)] = n
                lines.append(f"N,{n},{2 * i - 2},{2 * j - 2},{z}")
    lines.append("GROUP,1,SOLID")
    e = 0
    for k in range(len(depth_levels) - 1):
        lines.append(f"MACT,{k + 1}")
        for j in range(2):
            for i in range(2):
                e += 1
                b = [ids[(i, j, k)], ids[(i + 1, j, k)], ids[(i + 1, j + 1, k)], ids[(i, j + 1, k)]]
                t = [ids[(i, j, k + 1)], ids[(i + 1, j, k + 1)], ids[(i + 1, j + 1, k + 1)], ids[(i, j + 1, k + 1)]]
                lines.append("E," + ",".join(str(x) for x in [e] + b + t))
    lines.append(f"ETYPE,1,{e},1,2")
    lines.append(f"INT,1,{n},1")
    return lines


def _model(tmp_path, embedded=True):
    ui = Interpreter(cwd=tmp_path)
    plot_state(ui).auto_render = False
    lines = SITE.splitlines()
    if embedded:
        lines += _box([-2.0, -1.0, 0.0])
    else:                                   # a surface mat: shells at z = 0
        for j in range(3):
            for i in range(3):
                lines.append(f"N,{1 + i + 3 * j},{2 * i - 2},{2 * j - 2},0")
        lines += ["GROUP,1,SHELL", "E,1,1,2,5,4", "EGEN,1,1,1", "EGEN,1,3,1,2", "INT,1,9,1"]
    for ln in lines:
        assert ui.execute(ln), (ln, ui.sink.texts(Kind.ERROR)[-1:])
    return ui


def _inside(p, so):
    """Is point p inside the drawn soil (island box, not in the opening, not in the cut-away quarter)?"""
    X0, X1, Y0, Y1, zb, zg = so["bbox"]
    if not (X0 < p[0] < X1 and Y0 < p[1] < Y1 and zb < p[2] < zg):
        return False
    h = so["hole"]
    if h is not None and h[0] < p[0] < h[1] and h[2] < p[1] < h[3] and p[2] > h[4]:
        return False
    if so["cut"] is not None:
        sx, sy = so["cut"]
        cx, cy = so["center"]
        if (p[0] - cx) * sx > 0 and (p[1] - cy) * sy > 0:
            return False
    return True


def _area(q):
    return float(np.linalg.norm(np.cross(q[1] - q[0], q[3] - q[0])))


# ------------------------------------------------------------------ geometry
def test_embedded_island_layers_opening_and_cut(tmp_path):
    ui = _model(tmp_path)
    so = soil_island(ui.model, S.rotation_matrix(*S.DEFAULT_VIEW[:3]))
    # bands: the five TOPL layers at their depths (top down), then the half-space; softer = lighter
    tops = [b["z_top"] for b in so["bands"]]
    assert tops == [0.0, -1.0, -2.0, -4.0, -6.0, -8.0]
    assert [b["kind"] for b in so["bands"]] == ["layer"] * 5 + ["halfspace"]
    assert so["profile_depth"] == 8.0 and so["embedment"] == 2.0 and so["hole"] == [-2.0, 2.0, -2.0, 2.0, -2.0]
    lum = lambda c: sum(int(c[k:k + 2], 16) for k in (1, 3, 5))
    assert lum(so["bands"][0]["color"]) > lum(so["bands"][4]["color"])
    # automatic extent: the whole profile and a half-space band of a quarter of it; margin max(B/2, 0.3 d, D)
    assert so["depth"] == pytest.approx(10.0) and so["margin"] == pytest.approx(3.0)
    assert so["bbox"] == pytest.approx([-5.0, 5.0, -5.0, 5.0, -10.0, 0.0])
    # embedded: the quarter facing the default view (+X, -Y) is cut away through the foundation centre
    assert so["cut"] == [1.0, -1.0] and so["center"] == [0.0, 0.0]
    Q = np.asarray(so["quads"])
    N = np.asarray(so["quad_normal"])
    assert len(Q) and len(Q) == len(N) == len(so["quad_band"])
    eps = 1e-3
    for q, n in zip(Q, N):
        c = q.mean(axis=0)
        # every face bounds the soil, with its normal pointing out of it
        assert _inside(c - eps * n, so), (c, n)
        assert not _inside(c + eps * n, so), (c, n)
        # no face against the opening (the model fills it)
        h = so["hole"]
        o = c + eps * n
        assert not (h[0] < o[0] < h[1] and h[2] < o[1] < h[3] and o[2] > h[4])
    # the cut faces show the layers: faces with normals +X and -Y at the cut planes
    assert any(abs(q[:, 0] - 0).max() < 1e-9 and n[0] == 1 for q, n in zip(Q, N))
    assert any(abs(q[:, 1] - 0).max() < 1e-9 and n[1] == -1 for q, n in zip(Q, N))
    # the ground surface (top faces) = island plan minus the opening minus the cut quarter
    top = sum(_area(q) for q, n in zip(Q, N) if n[2] == 1)
    assert top == pytest.approx(100.0 - 16.0 - 25.0 + 4.0)
    # the layer interfaces become seams, the outline and creases edges
    zs = {round(float(z), 9) for s_ in np.asarray(so["seams"]) for z in s_[:, 2]}
    assert {-1.0, -2.0, -4.0, -6.0, -8.0} <= zs
    assert len(so["edges"]) > 0
    assert "half-space" in so["summary"] and "cut away" in so["summary"] and "embedded 2" in so["summary"]
    assert "infinity" in so["note"]


def test_surface_island_closed_below_the_mat(tmp_path):
    ui = _model(tmp_path, embedded=False)
    so = soil_island(ui.model, S.rotation_matrix(*S.DEFAULT_VIEW[:3]))
    assert so["hole"] is None and so["cut"] is None and so["embedment"] == 0.0     # automatic: no cut
    Q, N = np.asarray(so["quads"]), np.asarray(so["quad_normal"])
    X0, X1, Y0, Y1, zb, zg = so["bbox"]
    plan = (X1 - X0) * (Y1 - Y0)
    # the ground under the mat is the mat's own face: not drawn; the bottom is whole
    assert sum(_area(q) for q, n in zip(Q, N) if n[2] == 1) == pytest.approx(plan - 16.0)
    assert sum(_area(q) for q, n in zip(Q, N) if n[2] == -1) == pytest.approx(plan)
    # the four sides: perimeter x depth
    side = sum(_area(q) for q, n in zip(Q, N) if n[2] == 0)
    assert side == pytest.approx(2 * ((X1 - X0) + (Y1 - Y0)) * (zg - zb))
    # cut forced on, another view (from -X, +Y): that quarter goes
    so2 = soil_island(ui.model, S.rotation_matrix(-60, 0, 135), cut=1)
    assert so2["cut"] == [-1.0, 1.0]


def test_given_extent_and_truncated_profile(tmp_path):
    ui = _model(tmp_path)
    so = soil_island(ui.model, None, cut=0, margin=10, depth=5)
    assert so["bbox"] == pytest.approx([-12.0, 12.0, -12.0, 12.0, -5.0, 0.0]) and so["cut"] is None
    # layer 4 (-4 to -6) is clipped at 5, layer 5 and the half-space are not drawn
    assert [b["z_bot"] for b in so["bands"]] == [-1.0, -2.0, -4.0, -5.0]
    assert so["bands"][-1]["clipped"] and so["truncated"] == 1
    assert "cut at depth 5" in so["summary"] and "not drawn" in so["summary"]
    # deeper than the profile: the half-space fills the rest
    so = soil_island(ui.model, None, depth=20)
    assert so["bands"][-1]["kind"] == "halfspace" and so["bands"][-1]["z_bot"] == -20.0


def test_no_soil_layers(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    for ln in ("N,1,0,0,0", "N,2,1,0,0", "GROUP,1,BEAMS", "E,1,1,2"):
        ui.execute(ln)
    with pytest.raises(PlotError, match="no soil layers"):
        soil_island(ui.model)


# ------------------------------------------------------------------ command and plot data
def test_showsoil_command(tmp_path):
    ui = _model(tmp_path)
    st = plot_state(ui)
    # no plot: stored as the default of new 3D plots
    assert ui.execute("SHOWSOIL,1") and st.defaults3d.show_soil
    assert any("stored as default" in t for t in ui.sink.texts(Kind.CONFIRM))
    assert ui.execute("SHOWSOIL,0") and not st.defaults3d.show_soil
    assert ui.execute("MODELPLOT")
    v = st.active_plot.view
    assert not v.show_soil
    assert ui.execute("SHOWSOIL") and v.show_soil                   # toggle
    info = [t for t in ui.sink.texts(Kind.INFO) if t.startswith("SHOWSOIL:")]
    assert info and "5 layers" in info[-1] and "display only" in info[-1]
    assert ui.execute("SHOWSOIL,1,0,6,12") and (v.soil_cut, v.soil_margin, v.soil_depth) == (0, 6.0, 12.0)
    assert ui.execute("SHOWSOIL,1") and (v.soil_cut, v.soil_margin, v.soil_depth) == (0, 6.0, 12.0)   # blank keeps
    # wrong fields: an error, nothing changed
    for bad in ("SHOWSOIL,1,2", "SHOWSOIL,1,,-1", "SHOWSOIL,3"):
        assert not ui.execute(bad)
    assert (v.show_soil, v.soil_cut, v.soil_margin, v.soil_depth) == (True, 0, 6.0, 12.0)
    assert ui.execute("SHOWSOIL,0") and not v.show_soil
    # the soil layer plot is not capable: ignored with a warning
    assert ui.execute("LAYERPLOT") and ui.execute("SHOWSOIL,1")
    assert any("not applicable" in t for t in ui.sink.texts(Kind.WARNING))


def test_plot_data_box_and_fit_include_the_soil(tmp_path):
    ui = _model(tmp_path)
    st = plot_state(ui)
    assert ui.execute("NODEPLOT")
    p = st.active_plot
    d0 = plot_data(st, p, ui)
    assert d0["soil"] is None and d0["scene"]["bbox"] == pytest.approx([-2, 2, -2, 2, -2, 0])
    assert ui.execute("SHOWSOIL,1")
    d1 = plot_data(st, p, ui)
    assert d1["soil"]["bbox"] == pytest.approx([-5, 5, -5, 5, -10, 0])
    assert d1["scene"]["bbox"] == pytest.approx([-5, 5, -5, 5, -10, 0])
    assert d1["scene"]["center"] == pytest.approx([0, 0, -5])
    assert d1["camera"]["half"] > 2 * d0["camera"]["half"]
    assert isinstance(d1["soil"]["quads"], list) and isinstance(d1["soil"]["bands"][0]["vs"], float)   # JSON
    # a model without soil layers: a note instead of the soil
    ui2 = Interpreter(cwd=tmp_path / "b")
    plot_state(ui2).auto_render = False
    for ln in ("N,1,0,0,0", "N,2,1,0,0", "GROUP,1,BEAMS", "E,1,1,2", "MODELPLOT", "SHOWSOIL,1"):
        assert ui2.execute(ln), ln
    assert any("no soil to draw" in t for t in ui2.sink.texts(Kind.WARNING))
    d2 = plot_data(plot_state(ui2), plot_state(ui2).active_plot, ui2)
    assert d2["soil"] is None and "no soil layers" in d2["soil_note"]


def test_png_rendering_with_the_soil(tmp_path):
    ui = _model(tmp_path)
    plot_state(ui).auto_render = True
    for ln in ("MODELPLOT", "SHOWSOIL,1", "CAPTUREPLOT,soil_model.png", "NODEPLOT", "SHOWSOIL,1",
               "CAPTUREPLOT,soil_nodes.png"):
        assert ui.execute(ln), (ln, ui.sink.texts(Kind.ERROR)[-1:])
    pngs = sorted(tmp_path.rglob("soil_*.png"))
    assert len(pngs) == 2 and all(p.stat().st_size > 5000 for p in pngs)


def test_renderer_helpers():
    # ribbons: width across the segment, split into pieces
    quads, owner = _ribbons([np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 1.0]])], 0.2, 2.6)
    assert len(quads) == 4 and owner == [0, 0, 0, 0]
    assert quads[0][0][1] == pytest.approx(0.1) and quads[0][3][1] == pytest.approx(-0.1)
    assert quads[-1][1][0] == pytest.approx(10.0)
    # markers inside the drawn soil are hidden, those in the cut-away quarter and above grade are not
    so = {"bbox": [-5, 5, -5, 5, -10, 0], "cut": [1.0, -1.0], "center": [0.0, 0.0]}
    pts = np.array([[-1.0, 1.0, -1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]])
    assert _soil_hidden(pts, so).tolist() == [True, False, False, False]


# ------------------------------------------------------------------ GUI and lessons
def test_front_end_and_lesson_action():
    js = (STATIC / "plots.js").read_text(encoding="utf-8")
    assert 'SHOWSOIL: ["MODELPLOT", "NODEPLOT"]' in js
    assert '"Show the soil around the foundation (SHOWSOIL)", () => S.command("SHOWSOIL")' in js
    assert "function soilTraces" in js and 'role: "soil"' in js and "soilNote(d)" in js
    from sassi.ui import learn
    assert learn.resolve_action("plot-soil", "", ROOT, None, [], None) == \
        {"kind": "commands", "lines": ["MODELPLOT", "SHOWSOIL,1"]}
    assert learn.action_label("plot-soil", "") == "Plot the model in its soil"
    for lesson in ("05_embedded.md", "11_embedded_building.md"):
        text = (ROOT / "sassi" / "ui" / "lessons" / lesson).read_text(encoding="utf-8")
        assert "\nplot-soil\n" in text and "SHOWSOIL" in text
