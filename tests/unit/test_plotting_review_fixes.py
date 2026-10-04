"""Regression tests of the review of the plot state and plot commands (spec 06 sections 1.4, 5.1, 6;
spec 10 section 4; D-UI-10, D-UI-11; D-STR-08; rule L17).

Covered defects: display volume and node visibility, half-space depth of the layer table and the
omitted-layer marks, frame-directory order, CAPTUREPLOT format, display requests of a CUTPLOT of a
non-active model, line palettes with one colour per line number, ACTIVATEPLOT, and 2D defaults
that must not reach a soil-property plot.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from sassi.plotting import state as S
from sassi.plotting.state import plot_state, write_frame
from sassi.prep import Interpreter, Kind


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


@pytest.fixture
def ui(tmp_path, monkeypatch):
    monkeypatch.setenv("SASSI_EDU_SETTINGS_DIR", str(tmp_path / "settings"))
    u = Interpreter(cwd=tmp_path)
    plot_state(u).auto_render = False
    return u


def write_xy(path, x, y):
    path.write_text("\n".join(f"{float(a)!r} {float(b)!r}" for a, b in zip(x, y)) + "\n")


def box_model(ui, n=2, extra=""):
    """(n+1)^3 nodes of a unit grid, n^3 hexes (group 1 SOLID), plus extra command text."""
    out = []
    nid = lambda i, j, k: 1 + i + (n + 1) * j + (n + 1) ** 2 * k
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                out.append(f"N,{nid(i, j, k)},{i},{j},{k}")
    out.append("GROUP,1,SOLID")
    e = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                c = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k)]
                c += [v + (n + 1) ** 2 for v in c]
                out.append(f"E,{e}," + ",".join(map(str, c)))
                e += 1
    ui.run_text("\n".join(out) + "\n" + extra)
    assert not errors(ui)


# ------------------------------------------------------------------ display volume (spec 06 1.4, D-UI-11)
def test_display_volume_hides_nodes_and_markers_outside_the_box(ui):
    # 4 x 4 x 4 box; interaction nodes 1 6 11 16 21 (x = 0) and 2 (x = 1), a fixity and a mass at
    # the far corner node 125 (x = 4), node 2 fixed too
    box_model(ui, 4, "INT,1,21,5,1\nINT,2,2,1,1\nD,125,125,1,1,UX\nD,2,2,1,1,UX\nMT,125,1,1,1\n"
              "NODESEL,1,2,125\n")
    sc_all = S.model_scene(ui.model, show_dof=[0], show_mass=True)
    assert len(sc_all["node_id"]) == 125
    assert ui.execute("MODELPLOT") and ui.execute("WINDOWSETTINGS,VOLUME,0,0.5,,,,")
    st = plot_state(ui)
    sc = S.model_scene(ui.model, show_dof=[0], show_mass=True)
    assert sc["n_elements"] == 16                             # elements with a node at x <= 0.5
    x = sc["xyz"][:, 0]
    # only the nodes of the remaining elements are in the scene (x = 0 and x = 1) ...
    assert sorted(set(x.tolist())) == [0.0, 1.0] and len(sc["node_id"]) == 50
    # ... and those at x = 1 lie outside the box: drawn as element corners, but not as nodes
    vis = np.asarray(sc["node_visible"], bool)
    assert np.all(vis == (x <= 0.5))
    ids = sc["node_id"]
    assert sorted(ids[sc["interaction"]].tolist()) == [1, 6, 11, 16, 21]   # node 2 (x = 1) not shown
    assert len(sc["fixed"]) == 0 and len(sc["mass"]) == 0     # node 2 outside, node 125 clipped
    assert ids[sc["selected"]].tolist() == [1]
    # the camera fits the remaining geometry (x 0..1), not the clipped nodes (x up to 4)
    d = st.plot_data(st.active, ui, raw=True)
    assert np.max(d["scene"]["xyz"][:, 0]) == 1.0
    full = S.view_frame(st.active_plot.view, sc_all)["half"]
    R = S.rotation_matrix(*S.DEFAULT_VIEW[:3])
    slab = np.array([[i, j, k] for i in (0, 1) for j in range(5) for k in range(5)], float)
    scr = S.to_screen(slab, R, [0.5, 2.0, 2.0])               # centre of the remaining 1 x 4 x 4 slab
    assert d["camera"]["half"] == pytest.approx(1.08 * np.max(np.abs(scr[:, :2])), rel=1e-12)
    assert d["camera"]["half"] < full
    json.dumps(st.plot_data(st.active, ui))
    # HIDENODE hides a node's markers too (node-based plots)
    assert ui.execute("WINDOWSETTINGS,SHOWALL") and ui.execute("WINDOWSETTINGS,HIDENODE,1")
    sc = S.model_scene(ui.model)
    k = int(np.nonzero(sc["node_id"] == 1)[0][0])
    assert not sc["node_visible"][k] and k not in sc["interaction"].tolist() and k not in sc["selected"].tolist()


def test_node_plot_with_display_volume_renders(ui, tmp_path):
    from sassi.plotting.render_mpl import render_plot
    box_model(ui, 3)
    assert ui.execute("NODEPLOT") and ui.execute("WINDOWSETTINGS,VOLUME,,,,,0,0.2")
    st = plot_state(ui)
    out, notes = render_plot(st, st.active_plot, ui, tmp_path / "np.png")
    assert out.stat().st_size > 2000 and not notes


# ------------------------------------------------------------------ layer table (spec 06 6)
LAYERS = ("L,1,10,120,2000,1000,0.02,0.03\nL,2,20,125,2500,1200,0.02,0.02\nL,3,5,128,3000,1500,0.02,0.02\n"
          "L,4,0,130,4000,2000,0.01,0.01\nTOPL,1,2,3\nSITE,0,1,0,20,4\n")


def test_layer_table_halfspace_depth_and_omitted_layers(ui):
    ui.run_text(LAYERS)
    full = S.layer_table(ui.model)
    assert full["halfspace"]["top"] == 35.0 == full["total_depth"]
    assert (full["omitted_above"], full["omitted_below"]) == (0, 0)
    mid = S.layer_table(ui.model, 2, 2)
    assert [r["layer"] for r in mid["layers"]] == [2] and mid["layers"][0]["top"] == 10.0
    assert mid["halfspace"]["top"] == 35.0
    assert (mid["omitted_above"], mid["omitted_below"]) == (1, 1)


def test_layer_figure_marks_layers_not_shown(ui, tmp_path):
    from sassi.plotting.render_mpl import figure
    ui.run_text(LAYERS + "LAYERPLOT\nWINDOWSETTINGS,START,2\nWINDOWSETTINGS,END,2\n")
    assert not errors(ui)
    st = plot_state(ui)
    fig = figure(st.plot_data(st.active, ui, raw=True))
    texts = [t.get_text() for ax in fig.axes for t in ax.texts]
    assert "layer 1 not shown" in texts and "layer 3 not shown" in texts
    assert any(t.startswith("Halfspace") and "35" in t for t in texts)
    cells = [c.get_text().get_text() for ax in fig.axes for tab in ax.tables for c in tab.get_celld().values()]
    assert "1 (not shown)" in cells and "3 (not shown)" in cells


# ------------------------------------------------------------------ frame directories (D-STR-08)
def test_procframe_directory_stores_frames_in_time_order_beyond_100_s(ui, tmp_path):
    d = tmp_path / "THD"
    d.mkdir()
    for s in range(25):                                       # t = 0 .. 120 s every 5 s
        write_frame(d / f"THD_{5.0 * s:06.3f}_{s + 1:05d}", [1, 2], [[float(s), 0.0, 0.0], [0.0, 0.0, 0.0]])
    assert ui.execute("PROCFRAME,THD,store,THD frames,0")
    idx = json.loads((tmp_path / "store" / "index.json").read_text())
    assert [e["frame"] for e in idx["frames"]] == list(range(1, 26))
    assert [e["value"] for e in idx["frames"]] == [5.0 * s for s in range(25)]
    st = S.FrameStore(tmp_path / "store")
    assert [float(st.values(k)[0, 0]) for k in range(1, 26)] == [float(s) for s in range(25)]
    assert idx["notes"] == [] and not warnings(ui)


def test_frame_directory_mixing_sequences_is_ordered_and_warned(ui, tmp_path):
    d = tmp_path / "RS"
    d.mkdir()
    for tag in ("RS02", "RS01"):
        for q, f in enumerate((0.5, 50.0, 999.0, 1000.0)):
            write_frame(d / f"{tag}_{f:06.2f}_{q + 1:05d}", [1], [[f, 0.0, 0.0]])
    names = [p.name for p in S.read_frame_list(d)]
    assert names[:4] == ["RS01_000.50_00001", "RS01_050.00_00002", "RS01_999.00_00003", "RS01_1000.00_00004"]
    assert names[4].startswith("RS02_000.50")
    notes = S.frame_directory_notes(S.read_frame_list(d))
    assert len(notes) == 1 and "RS01" in notes[0] and "RS02" in notes[0]
    assert ui.execute("PROCFRAME,RS,rsstore,RS,0")
    assert any("mixes 2 frame sequences" in w for w in warnings(ui))
    # files outside the frame-name pattern come last, with a note
    write_frame(tmp_path / "RS" / "zzz_extra", [1], [[0.0, 0.0, 0.0]])
    write_frame(tmp_path / "RS" / "aaa_extra", [1], [[0.0, 0.0, 0.0]])
    files = S.read_frame_list(d)
    assert [f.name for f in files[-2:]] == ["aaa_extra", "zzz_extra"]
    assert any("frame-name pattern" in n for n in S.frame_directory_notes(files))


# ------------------------------------------------------------------ display requests follow the plotted model
def test_windowsettings_and_nodesel_act_on_the_model_of_the_active_3d_plot(ui):
    box_model(ui, 2, "GROUP,2,SHELL\nE,1,1,2,5,4\n")
    assert ui.execute("CPMODEL,1")
    ui.session["cuts"] = {1: [(1, 1)]}
    assert ui.execute("CUTPLOT,1,1")                          # model 1 while model 0 is active
    st = plot_state(ui)
    assert st.active_plot.model == 1 and ui.active_model == 0
    assert ui.execute("WINDOWSETTINGS,HIDEGROUP,2")
    assert ui.models[1].ui_state.get("hide_groups") == [2]
    assert not ui.models[0].ui_state.get("hide_groups")
    assert ui.execute("NODESEL,5")
    assert ui.models[1].ui_state.get("nodesel") == [5] and not ui.models[0].ui_state.get("nodesel")
    assert S.model_scene(ui.models[1])["n_elements"] == 8     # group 2 hidden in model 1 only
    # with a 2D plot in front the active model is edited (no 3D plot to follow)
    write_xy(ui.cwd / "s.rs", [1, 2], [1, 2])
    ui.run_text("READSPEC,s.rs,1,1\nSPECPLOT,1\n")
    assert ui.execute("NODESEL,7")
    assert ui.models[0].ui_state.get("nodesel") == [7] and ui.models[1].ui_state.get("nodesel") == [5]


# ------------------------------------------------------------------ line palettes (spec 06 5.1)
def test_line_palettes_have_one_colour_per_line_number(ui, tmp_path):
    for pal in ("SPECLINES", "THLINES"):
        cols = [S.palette_color(pal, n) for n in range(1, 201)]
        assert len(set(cols)) == 200                          # never cyclic
        assert len(S.PALETTES[pal]) == S.LINE_PALETTE_SIZE == 50
        assert S.PALETTES[pal] == cols[:50]
    assert S.palette_color("SPECLINES", 17) != S.palette_color("SPECLINES", 1)
    assert S.palette_color("ElemPalette", 129) == S.palette_color("ElemPalette", 1)   # still cyclic (D-UI-09)
    # COLOR accepts any line number for the line palettes and the plot uses it
    for n in (1, 17, 77):
        write_xy(tmp_path / f"l{n}.rs", [1, 2], [n, n])
        assert ui.execute(f"READSPEC,l{n}.rs,1,{n}")
    assert ui.execute("COLOR,SpecLines,77,1,2,3")
    assert not ui.execute("COLOR,SpecLines,0,1,2,3")
    assert not ui.execute("COLOR,Node,6,1,2,3")               # fixed palettes keep their size
    assert ui.execute("SPECPLOT,1,17,77")
    st = plot_state(ui)
    colors = [L["color"] for L in st.plot_data(st.active)["lines"]]
    assert colors[2] == "#010203" and len(set(colors)) == 3


# ------------------------------------------------------------------ ACTIVATEPLOT (L17)
def test_activateplot_switches_the_target_of_setting_commands(ui, tmp_path):
    write_xy(tmp_path / "s.rs", [1, 2], [1, 2])
    ui.run_text("READSPEC,s.rs,1,1\nSPECPLOT,1\nTHPLOT,1\n")
    st = plot_state(ui)
    p1, p2 = st.order
    seen = []
    st.subscribe(seen.append)
    assert st.active == p2
    assert ui.execute(f"ACTIVATEPLOT,{p1}")
    assert st.active == p1
    assert seen[-1].kind == "activate" and seen[-1].plot_id == p1 and seen[-1].data["plot"]["kind"] == "SPECPLOT"
    json.dumps(seen[-1].to_dict())
    assert ui.execute("XTITLE,Frequency (Hz)")
    assert st.plots[p1].settings.xtitle == "Frequency (Hz)" and st.plots[p2].settings.xtitle == ""
    assert not ui.execute("ACTIVATEPLOT,99") and any("open plots: 1, 2" in e for e in errors(ui))
    assert not ui.execute("ACTIVATEPLOT") and st.active == p1
    assert not ui.execute("ACTI,1")                           # extension command: full name only (L11)


def test_session_with_tab_switches_replays_identically(tmp_path, monkeypatch):
    monkeypatch.setenv("SASSI_EDU_SETTINGS_DIR", str(tmp_path / "settings"))
    write_xy(tmp_path / "s.rs", [1, 2], [1, 2])
    script = "READSPEC,s.rs,1,1\nSPECPLOT,1\nTHPLOT,1\nACTIVATEPLOT,1\nYTITLE,Sa (g)\nACTIVATEPLOT,2\nYTITLE,a (g)\n"
    snaps = []
    for _ in range(2):
        u = Interpreter(cwd=tmp_path)
        plot_state(u).auto_render = False
        u.run_text(script)
        assert not errors(u)
        snaps.append(plot_state(u).to_dict())
    assert snaps[0] == snaps[1]
    assert [p["settings"]["ytitle"] for p in snaps[0]["plots"]] == ["Sa (g)", "a (g)"]


# ------------------------------------------------------------------ 2D defaults (spec 10 4.0)
def test_2d_defaults_apply_to_line_plots_only(ui, tmp_path):
    write_xy(tmp_path / "s.rs", [0.1, 1, 10], [1, 2, 3])
    ui.run_text("READSPEC,s.rs,1,1\nPLOTRANGE,0.1,100,0,5\nXTITLE,Frequency (Hz)\nAXES,1,1,1,0,0,0\n"
                "DYNP,1,0.0001,1.0,0.0001,0.5,Clay\nDYNP,2,0.01,0.6,0.01,5,Clay\n")
    assert not errors(ui)
    st = plot_state(ui)
    assert ui.execute("SOILPROPPLOT,Clay")
    s = st.active_plot.settings
    assert s.log_x and not s.log_y and s.xmin is None and s.xmax is None and s.xtitle == ""
    assert not s.minor_x
    assert ui.execute("SPECPLOT,1")
    s = st.active_plot.settings
    assert (s.xmin, s.xmax, s.ymin, s.ymax) == (0.1, 100.0, 0.0, 5.0) and s.xtitle == "Frequency (Hz)"
    assert s.minor_x and not s.log_x
