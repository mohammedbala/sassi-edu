"""Plot state, plot commands and headless rendering (requirements 3.4.L, 5.7, UI-03; spec 06; spec 10 4.0-4.6;
D-UI-06..12, D-UI-08)."""
from __future__ import annotations

import json

import numpy as np
import pytest

from sassi.plotting import state as S
from sassi.plotting.state import PlotState, plot_state
from sassi.prep import Interpreter, Kind


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


@pytest.fixture
def ui(tmp_path):
    u = Interpreter(cwd=tmp_path)
    plot_state(u).auto_render = False
    return u


MODEL = """N,1,0,0,0
N,2,1,0,0
N,3,1,1,0
N,4,0,1,0
N,5,0,0,1
N,6,1,0,1
N,7,1,1,1
N,8,0,1,1
N,9,2,0,0
N,10,2,1,0
N,11,2,0,1
N,12,2,1,1
N,13,1,0.5,3
N,14,5,5,5
N,15,3,0,0
N,16,0,0,5
GROUP,1,SOLID
E,1,1,2,3,4,5,6,7,8
MACT,2
E,2,2,9,10,3,6,11,12,7
MACT,1
GROUP,2,SHELL
E,1,5,6,7,8
GROUP,3,BEAMS
RACT,4
E,1,7,13,14
RACT,1
GROUP,4,SPRING
E,1,9,15
E,2,15,15
INT,1,4,1,1
D,1,4,1,1,UX,UY
MT,13,10,10,0
"""


@pytest.fixture
def model_ui(ui):
    ui.run_text(MODEL)
    assert not errors(ui)
    return ui


def write_xy(path, x, y):
    path.write_text("\n".join(f"{float(a)!r} {float(b)!r}" for a, b in zip(x, y)) + "\n")


# ------------------------------------------------------------------ state and events
def test_open_close_activate_and_events():
    st = PlotState()
    seen = []
    st.subscribe(seen.append)
    p1 = st.open_plot("SPECPLOT", {"lines": [1]})
    p2 = st.open_plot("MODELPLOT", {}, model=0)
    assert st.active == p2.id and p2.caption == "Model 0 - Model Plot" and p1.caption == "Spectrum Plot"
    st.emit("open", p2.id, plot=p2.to_dict())
    assert seen[-1].kind == "open" and seen[-1].to_dict()["plot_id"] == p2.id
    assert st.close_plot().id == p2.id and st.active == p1.id
    assert st.close_plot().id == p1.id and st.active is None and st.close_plot() is None
    st.unsubscribe(seen.append)
    json.dumps(st.to_dict())


def test_target_rules():
    st = PlotState()
    assert st.target("AXES") == (None, "defaults")
    assert st.target("PLOTTITLE") == (None, "none")
    p = st.open_plot("LAYERPLOT", {}, model=0)
    assert st.target("PLOTTITLE") == (p, "ignored")      # the soil layer plot ignores PLOTTITLE
    assert st.target("AXES")[1] == "ignored"
    q = st.open_plot("SOILPROPPLOT", {"name": "x"}, model=0)
    assert st.target("YTITLE2") == (q, "plot") and q.settings.log_x       # strain axis log by default
    m = st.open_plot("NODEPLOT", {}, model=0)
    assert st.target("SHRINK")[1] == "ignored" and st.target("NODENUM") == (m, "plot")


def test_rotation_matrix_default_view_is_isometric_with_z_up():
    R = S.rotation_matrix(*S.DEFAULT_VIEW[:3])
    assert np.allclose(R @ R.T, np.eye(3))
    z = R @ np.array([0.0, 0.0, 1.0])
    assert z[1] > 0.8 and z[2] > 0                          # Z points up on the screen, tilted to the viewer
    assert np.allclose(S.rotation_matrix(0, 0, 0), np.eye(3))
    assert np.allclose(S.rotation_matrix(0, 0, 90) @ [1, 0, 0], [0, 1, 0])


def test_palettes_jet_and_bubble_rules():
    assert S.palette_color("ElemPalette", 130) == S.palette_color("ElemPalette", 2)       # spec 06 12.12
    assert len(set(S.PALETTES["ELEMPALETTE"])) == 128
    assert S.palette_color("SPECLINES", 1, {"SPECLINES": {1: "#123456"}}) == "#123456"
    assert np.allclose(S.jet(-1.0), [0, 0, 0.5]) and np.allclose(S.jet(2.0), [0.5, 0, 0])
    t = S.colormap_position(np.array([-11.0, 0.0, 11.0]), -10.0, 10.0)
    assert t.tolist() == [0.0, 0.5, 1.0]
    mid = S.jet(0.5)
    assert mid[1] == pytest.approx(1.0) and mid[0] < 0.6 and mid[2] < 0.6   # green
    assert S.colorbar_labels(-0.47718, 0.43260) == ["0.43260", "0.20515", "-0.02229", "-0.24974", "-0.47718"]
    assert S.bubble_size(np.array([0.0, 0.5, 1.0]), 10).tolist() == [2.0, 5.0, 10.0]


def test_vector_and_deformed_rules():
    V = S.vector_rule(np.array([[5.0, 0.0, 0.0]]), np.array([[0.3, 0.0, 0.0]]), 1.0)
    assert np.allclose(V[0, 0], [5.0, 0.3, 0.3])            # manual example 5 + 0.3i
    V = S.vector_rule(np.array([[0.0, 2.0, 3.0]]), np.array([[0.0, 0.1, 0.2]]), 2.0)
    assert np.allclose(V[0, 1], [0.2, 4.0, 0.2]) and np.allclose(V[0, 2], [0.4, 0.4, 6.0])   # D-UI-13
    assert np.allclose(S.deformed_coordinates([[1.0, 2.0, 3.0]], [[0.01, 0.0, 0.0]], 10.0), [[1.1, 2.0, 3.0]])
    re, im = S.vector_components(np.array([[1, 2, 3, 4, 5, 6.0]]), "complex_xyz")
    assert re.tolist() == [[1, 3, 5]] and im.tolist() == [[2, 4, 6]]


def test_frame_sequence():
    assert len(S.frame_sequence(1, 3110, 1, 3110)) == 3110
    assert S.frame_sequence(1, 10, 3, 3110) == [1, 4, 7, 10]
    with pytest.raises(S.PlotError):
        S.frame_sequence(5, 4, 1, 10)
    with pytest.raises(S.PlotError):
        S.frame_sequence(1, 11, 1, 10)


# ------------------------------------------------------------------ model scene
def test_model_scene_geometry_and_markers(model_ui):
    m = model_ui.model
    sc = S.model_scene(m, color_by=1, show_dof=[0], show_mass=True)
    ids = sc["node_id"].tolist()
    assert 14 not in ids and 16 not in ids                  # beam K node and unconnected node excluded
    assert sc["n_elements"] == 6
    assert len(sc["faces"]) == 2 * 6 + 1                    # two hexes + one shell
    assert int((~sc["face_boundary"]).sum()) == 2           # the shared hex face, counted per element
    assert len(sc["edges"]) == 2 and len(sc["points"]) == 1  # beam + spring, zero-length spring
    assert sorted(sc["node_id"][sc["interaction"]].tolist()) == [1, 2, 3, 4]
    assert sorted(sc["node_id"][sc["fixed"]].tolist()) == [1, 2, 3, 4]
    assert sc["node_id"][sc["mass"]].tolist() == [13] and sc["mass_dir"].tolist() == [[True, True, False]]
    assert sc["bbox"] == [0.0, 3.0, 0.0, 1.0, 0.0, 3.0]
    assert len(sc["tri"]) == 2 * len(sc["faces"])
    by_mat = S.model_scene(m, color_by=2)
    keys = dict(zip(zip(by_mat["elem_group"].tolist(), by_mat["elem_id"].tolist()), by_mat["elem_key"].tolist()))
    assert keys[(1, 1)] == 1 and keys[(1, 2)] == 2 and keys[(4, 1)] == 0      # springs: default colour
    by_prop = S.model_scene(m, color_by=3)
    keys = dict(zip(zip(by_prop["elem_group"].tolist(), by_prop["elem_id"].tolist()), by_prop["elem_key"].tolist()))
    assert keys[(3, 1)] == 4 and keys[(1, 1)] == 0
    json.dumps(S.jsonable(sc))


def test_model_scene_hide_requests_and_degenerate_solid(ui):
    ui.run_text("N,1,0,0,0\nN,2,1,0,0\nN,3,0,1,0\nN,4,0,0,1\nN,5,1,0,1\nN,6,0,1,1\n"
                "GROUP,1,SOLID\nE,1,1,2,3,3,4,5,6,6\nGROUP,2,SHELL\nE,1,1,2,5,4\n")
    m = ui.model
    sc = S.model_scene(m)
    sizes = sorted(int((f >= 0).sum()) for f in sc["faces"][sc["face_elem"] == 0])
    assert sizes == [3, 3, 4, 4, 4]                         # prism: 2 triangles + 3 quads
    m.ui_state["hide_groups"] = [2]
    assert S.model_scene(m)["n_elements"] == 1
    m.ui_state["hide_elements"] = [[1, 1]]
    assert S.model_scene(m)["n_elements"] == 0
    m.ui_state.clear()
    m.ui_state["display_volume"] = [5, 6, 5, 6, 5, 6]
    assert S.model_scene(m)["n_elements"] == 0
    assert S.model_scene(m, apply_hide=False)["n_elements"] == 2


def test_layer_table_and_soil_curves(ui):
    ui.run_text("L,1,10,120,2000,1000,0.02,0.03\nL,2,20,125,2500,1200,0.02,0.02\nL,3,0,130,4000,2000,0.01,0.01\n"
                "TOPL,2,1\nSITE,0,1,0,20,3\nDYNP,2,0.001,0.9,0.002,1.5,Sand\nDYNP,1,0.0001,1.0,0.0001,0.5,Sand\n"
                "DYNP,3,0.01,0.7,,,Sand\n")
    t = S.layer_table(ui.model)
    assert [r["layer"] for r in t["layers"]] == [2, 1] and t["halfspace"]["layer"] == 3
    assert t["layers"][1]["top"] == 20.0 and t["source"] == "TOPL"
    assert [r["layer"] for r in S.layer_table(ui.model, 2, -1)["layers"]] == [1]
    with pytest.raises(S.PlotError):
        S.layer_table(ui.model, 3, 2)
    cv = S.soil_property_curves(ui.model, "Sand")
    assert cv["g_strain"] == [0.0001, 0.001, 0.01] and cv["g"] == [1.0, 0.9, 0.7]
    assert cv["d_strain"] == [0.0001, 0.002] and cv["d"] == [0.5, 1.5]
    with pytest.raises(S.PlotError):
        S.soil_property_curves(ui.model, "sand")          # case-sensitive


# ------------------------------------------------------------------ 2D plot commands
def test_specplot_line_list_rules(ui, tmp_path):
    write_xy(tmp_path / "a.rs", [1, 2, 3], [1, 4, 9])
    ui.run_text("READSPEC,a.rs,1,1\nREADSPEC,a.rs,1,2\nREADSPEC,a.rs,1,3\n")
    st = plot_state(ui)
    assert ui.execute("SPECPLOT,1,2,-1,3")                   # -1 ends the list
    assert st.active_plot.params["lines"] == [1, 2]
    assert ui.execute("THPLOT,3,77")                          # unknown numbers ignored with a warning
    assert st.active_plot.params["lines"] == [3] and any("77" in w for w in warnings(ui))
    assert not ui.execute("SPECPLOT,55")
    assert not ui.execute("SPECPLOT")                          # Line Selection dialog: batch error (D-UI-08)
    assert any("D-UI-08" in e for e in errors(ui))
    d = st.plot_data(st.active)
    assert d["family"] == "2d" and d["lines"][0]["color"] == S.palette_color("THLINES", 3)
    json.dumps(d)


def test_2d_settings_defaults_and_validation(ui, tmp_path):
    st = plot_state(ui)
    ui.execute("AXES,1,1,1,1,1,0")                            # no plot: becomes the default
    ui.execute("XTITLE,Frequency, Hz")
    write_xy(tmp_path / "a.rs", [0.1, 1, 10], [1, 2, 1])
    ui.run_text("READSPEC,a.rs,1,1\n")
    ui.execute("SPECPLOT,1")
    p = st.active_plot
    assert p.settings.log_x and p.settings.minor_x and p.settings.xtitle == "Frequency, Hz"
    assert not ui.execute("PLOTRANGE,10,1")
    assert not ui.execute("PLOTRANGE,0,10")                   # log axis: minimum must be > 0
    assert ui.execute("PLOTRANGE,0.5,5,,3")
    assert (p.settings.xmin, p.settings.xmax, p.settings.ymin, p.settings.ymax) == (0.5, 5.0, None, 3.0)
    assert st.plot_data(p.id)["extent"] == [0.5, 5.0, 1.0, 3.0]
    ui.execute("PLOTTITLE,ISRS, node 17")
    assert p.title == "ISRS, node 17"
    ui.execute("STIPPLE")
    assert p.settings.stipple and st.plot_data(p.id)["lines"][0]["dash"] == "solid"
    ui.execute("MARKERS,1,1,5")
    assert st.lines[1].markers and any("5" in w for w in warnings(ui))
    ui.execute("LINENAME,1,Node 17, SRSS")
    assert st.lines[1].name == "Node 17, SRSS"
    assert not ui.execute("LINENAME,9,x")
    assert not ui.execute("AXES,2")


def test_plottitle_ignored_by_layer_plot_and_3d_commands_ignored_on_2d(ui, tmp_path):
    ui.run_text("L,1,10,120,2000,1000,0.02,0.03\nTOPL,1\n")
    assert ui.execute("LAYERPLOT")
    ui.execute("PLOTTITLE,x")
    assert any("not applicable" in w for w in warnings(ui))
    assert plot_state(ui).active_plot.title == ""
    ui.execute("SHRINK,1")
    assert sum("not applicable" in w for w in warnings(ui)) == 2


def test_closeplot_and_captureplot_formats(ui, tmp_path):
    write_xy(tmp_path / "a.rs", [1, 2, 3], [1, 4, 9])
    ui.execute("READSPEC,a.rs,1,1")
    assert not ui.execute("CAPTUREPLOT,x.png")               # no active plot
    ui.execute("SPECPLOT,1")
    ui.execute("THPLOT,1")
    seen = []
    plot_state(ui).subscribe(seen.append)
    assert ui.execute("CAPTUREPLOT,cap.PNG.tmp")              # '.png' anywhere, any case -> PNG
    assert (tmp_path / "cap.PNG.tmp").read_bytes()[:4] == b"\x89PNG"
    assert ui.execute("CAPTUREPLOT,cap.bmp")
    assert (tmp_path / "cap.bmp").read_bytes()[:2] == b"BM"
    assert [e.kind for e in seen] == ["capture", "capture"] and seen[1].data["format"] == "bmp"
    assert ui.execute("CLOSEPLOT") and plot_state(ui).active_plot.kind == "SPECPLOT"
    assert ui.execute("CLOSEPLOT") and plot_state(ui).active is None
    assert ui.execute("CLOSEPLOT") and any("no plot is open" in w for w in warnings(ui))


def test_auto_render_writes_images_in_the_model_directory(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    (tmp_path / "mdl").mkdir()
    write_xy(tmp_path / "a.rs", [1, 2, 3], [1, 4, 9])
    ui.run_text("MDL,demo," + str(tmp_path / "mdl") + "\n")
    ui.execute("READSPEC," + str(tmp_path / "a.rs") + ",1,1")
    assert ui.execute("SPECPLOT,1")
    img = tmp_path / "mdl" / "demo_plot01_specplot.png"
    assert img.exists() and plot_state(ui).active_plot.image == str(img)
    assert any(str(img) in t for t in ui.sink.texts(Kind.INFO))
    st = plot_state(ui)
    st.auto_render = False
    st.dialog_handler = lambda name, ctx: True               # a GUI with dialogs
    events = []
    st.subscribe(events.append)
    assert ui.execute("THPLOT")                              # dialog requested instead of an error
    assert events[-1].kind == "dialog" and events[-1].data["dialog"] == "Line Selection"
    assert ui.execute("THPLOT,1") and events[-1].kind == "open"
    assert not (tmp_path / "mdl" / "demo_plot02_thplot.png").exists()


# ------------------------------------------------------------------ 3D plot commands
def test_modelplot_view_and_label_commands(model_ui):
    ui = model_ui
    st = plot_state(ui)
    ui.execute("ELENUM,1")                                   # no plot: default for new plots
    assert ui.execute("MODELPLOT")
    v = st.active_plot.view
    assert v.elem_labels
    ui.execute("GROUPNUM")                                   # toggle on -> ELENUM off
    assert v.group_labels and not v.elem_labels
    ui.execute("ELECOLOR,3")
    assert v.color_by == 3 and not ui.execute("ELECOLOR,4")
    ui.execute("SHRINK")
    ui.execute("WIREFRAME,1")
    assert v.shrink and v.wireframe
    ui.execute("DEBUG")
    assert v.debug
    ui.execute("DEBUG,2")
    assert not v.debug
    ui.execute("CNGVIEW,10,20,30,0.5,,2")
    assert (v.rx, v.ry, v.rz, v.px, v.py, v.zoom) == (10.0, 20.0, 30.0, 0.5, 0.0, 2.0)
    assert not ui.execute("CNGVIEW,,,,,,0")
    ui.execute("RSTVIEW")
    assert (v.rx, v.ry, v.rz, v.zoom) == (-60.0, 0.0, -45.0, 1.0)
    ui.execute("CNGCENTER,1,2,3")
    cam = S.view_frame(v, S.model_scene(ui.model))
    assert cam["center"].tolist() == [1.0, 2.0, 3.0]
    ui.execute("RSTCENTER")
    cam = S.view_frame(v, S.model_scene(ui.model))
    assert cam["center"].tolist() == [1.5, 0.5, 1.5]
    ui.execute("SHOWDOF,DISP")
    assert v.show_dof == [0, 1, 2]
    ui.execute("SHOWDOF,NONE")
    assert v.show_dof == []
    assert not ui.execute("SHOWDOF")                          # dialog in batch mode
    assert not ui.execute("SHOWDOF,Q")
    ui.execute("SHOWMASS,1")
    assert v.show_mass
    ui.execute("PAUSE")
    assert any("not applicable" in w for w in warnings(ui))
    d = st.plot_data(st.active, ui)
    json.dumps(d)
    assert d["family"] == "3d" and len(d["scene"]["elem_color"]) == 6


def test_nodeplot_nodesel_toggles_at_model_level(model_ui):
    ui = model_ui
    assert ui.execute("NODEPLOT")
    ui.execute("NODESEL,13,7,99")
    assert ui.model.ui_state["nodesel"] == [13, 7] and any("99" in w for w in warnings(ui))
    ui.execute("NODESEL,13")
    assert ui.model.ui_state["nodesel"] == [7]
    sc = S.model_scene(ui.model)
    assert sc["node_id"][sc["selected"]].tolist() == [7]
    assert not ui.execute("NODESEL")


def test_modelplot_needs_elements_and_cutplot(model_ui):
    ui = model_ui
    ui.execute("ACTM,1")
    assert not ui.execute("MODELPLOT")
    assert not ui.execute("CUTPLOT")                          # Select Cut to Display dialog: batch error
    assert not ui.execute("CUTPLOT,1,0")                      # no cut defined
    ui.session["cuts"] = {1: [(1, 2), (3, 1), (9, 9)]}
    assert ui.execute("CUTPLOT,1,0")                          # model 0 while model 1 is active
    st = plot_state(ui)
    p = st.active_plot
    assert p.model == 0 and p.caption == "Model 0 - Cut # 1 Plot"
    assert any("not in model 0" in w for w in warnings(ui))
    d = st.plot_data(p.id, ui)
    incut = d["scene"]["elem_in_cut"]
    assert sum(incut) == 2


def test_windowsettings_field_form(model_ui):
    ui = model_ui
    st = plot_state(ui)
    assert not ui.execute("WINDOWSETTINGS")                   # dialog: batch error
    ui.execute("MODELPLOT")
    assert ui.execute("WINDOWSETTINGS,HIDEGROUP,2")
    assert ui.model.ui_state["hide_groups"] == [2]
    assert ui.execute("WINDOWSETTINGS,HIDEELEM,1,1-2")
    assert ui.model.ui_state["hide_elements"] == [[1, 1], [1, 2]]
    assert S.model_scene(ui.model)["n_elements"] == 3
    assert ui.execute("WINDOWSETTINGS,SHOWALL")
    assert S.model_scene(ui.model)["n_elements"] == 6
    assert ui.execute("WINDOWSETTINGS,VOLUME,,0.5")
    assert ui.model.ui_state["display_volume"][:2] == [0.0, 0.5]
    assert ui.execute("WINDOWSETTINGS,TITLE,Model, top view")
    assert st.active_plot.title == "Model, top view"
    assert not ui.execute("WINDOWSETTINGS,SCALE,2")           # not an animation
    assert not ui.execute("WINDOWSETTINGS,NOPE,1")
    ui.run_text("L,1,10,120,2000,1000,0.02,0.03\nTOPL,1\nLAYERPLOT\n")
    assert ui.execute("WINDOWSETTINGS,SHOW,THICK,1") and st.active_plot.params["show"]["thick"]
    assert not ui.execute("WINDOWSETTINGS,START,5")


def test_shaderoptions_and_color(ui):
    st = plot_state(ui)
    assert not ui.execute("SHADEROPTIONS")
    assert ui.execute("SHADEROPTIONS,12,,0.1")
    assert (st.shader.points, st.shader.linew, st.shader.shrink, st.shader.scale) == (12.0, 0.02, 0.1, 1.0)
    assert not ui.execute("SHADEROPTIONS,,0.7")
    assert ui.execute("COLOR,SpecLines,1,255,0,0")
    assert st.palettes["SPECLINES"][1] == "#ff0000"
    assert not ui.execute("COLOR,Nope,1,0,0,0") and not ui.execute("COLOR,Node,1,0,0,300")


def test_soil_plot_commands(ui):
    ui.run_text("L,1,10,120,2000,1000,0.02,0.03\nL,2,0,130,4000,2000,0.01,0.01\nTOPL,1\nSITE,0,1,0,20,2\n"
                "DYNP,1,0.0001,1.0,0.0001,0.5,Clay\nDYNP,2,0.01,0.6,0.01,5,Clay\n")
    st = plot_state(ui)
    assert ui.execute("LAYERPLOT")
    d = st.plot_data(st.active, ui)
    assert d["table"]["halfspace"]["layer"] == 2
    assert not ui.execute("SOILPROPPLOT,clay") and ui.execute("SOILPROPPLOT,Clay")
    assert not ui.execute("SOILPROPPLOT")
    assert st.plot_data(st.active, ui)["curves"]["d"] == [0.5, 5.0]
    ui.execute("ACTM,3")
    assert not ui.execute("LAYERPLOT")


# ------------------------------------------------------------------ rendering (headless)
def test_render_every_plot_kind(model_ui, tmp_path):
    from sassi.plotting.render_mpl import figure, image_format, render_plot
    ui = model_ui
    st = plot_state(ui)
    write_xy(tmp_path / "a.rs", [0.1, 1, 10], [1, 2, 1])
    ui.run_text("READSPEC,a.rs,1,1\nL,1,10,120,2000,1000,0.02,0.03\nTOPL,1\nDYNP,1,0.0001,1.0,0.0001,0.5,Clay\n"
                "DYNP,2,0.01,0.6,0.01,5,Clay\n")
    ui.session["cuts"] = {1: [(1, 1)]}
    cmds = ["SPECPLOT,1", "THPLOT,1", "LAYERPLOT", "SOILPROPPLOT,Clay", "MODELPLOT", "NODEPLOT", "CUTPLOT,1"]
    for k, cmd in enumerate(cmds):
        assert ui.execute(cmd), errors(ui)
        if cmd == "MODELPLOT":
            ui.run_text("NODENUM,1\nELENUM,1\nSHOWDOF,ALL\nSHOWMASS,1\nSHRINK,1\nDEBUG,1\nNODESEL,7\n")
        out, notes = render_plot(st, st.active_plot, ui, tmp_path / f"k{k}.png")
        assert out.stat().st_size > 2000 and not notes
    assert image_format("a.PnG.x") == "png" and image_format("a.bmp") == "bmp" and image_format("a") == "bmp"
    assert not errors(ui)


def test_export_table_csv(ui, tmp_path):
    write_xy(tmp_path / "a.rs", [1, 2, 3], [1, 4, 9])
    write_xy(tmp_path / "b.rs", [1.5, 2.5], [10, 20])
    ui.run_text("READSPEC,a.rs,1,1\nREADSPEC,b.rs,1,2\nSPECPLOT,1,2\n"
                "DYNP,1,0.0001,1.0,0.0001,0.5,Clay\nDYNP,2,0.01,0.6,,,Clay\n")
    st = plot_state(ui)
    p = S.export_table_csv(st, st.active_plot, tmp_path / "t.csv")
    rows = p.read_text().splitlines()
    assert rows[0] == "Frequency,a.rs,b.rs" and len(rows) == 1 + 5
    assert rows[2].split(",") == ["1.5", "2.5", "10.0"]
    ui.execute("SOILPROPPLOT,Clay")
    p = S.export_table_csv(st, st.active_plot, tmp_path / "s.csv", ui)
    assert p.read_text().splitlines()[2] == "0.01,0.6,,"
    st.open_plot("MODELPLOT", {}, model=0)
    with pytest.raises(S.PlotError):
        S.export_table_csv(st, st.active_plot, tmp_path / "m.csv")


def test_layer_plot_end_layer_minus_one_and_layers_without_topl(ui):
    ui.run_text("\n".join(f"L,{k},{k},120,2000,1000,0.02,0.02" for k in range(1, 9)) + "\nSITE,0,1,0,20,8\n")
    t = S.layer_table(ui.model)                    # no TOPL: every L layer except the half-space layer 8
    assert [r["layer"] for r in t["layers"]] == [1, 2, 3, 4, 5, 6, 7] and t["halfspace"]["layer"] == 8
    assert t["source"] == "L" and t["end"] == 7     # EndLayer -1 = deepest layer (spec 06 12.14)
    ui.execute("LAYERPLOT")
    st = plot_state(ui)
    assert not ui.execute("WINDOWSETTINGS,START,9") and st.active_plot.params["start"] == 1
    assert ui.execute("WINDOWSETTINGS,END,3")
    assert [r["layer"] for r in st.plot_data(st.active, ui)["table"]["layers"]] == [1, 2, 3]


def test_headless_3d_projection_is_the_view_matrix(model_ui):
    """The renderer views the rotated model from the top (orthographic): screen x, y = R (p - c)."""
    from mpl_toolkits.mplot3d import proj3d
    from sassi.plotting.render_mpl import model_figure
    ui = model_ui
    ui.execute("MODELPLOT")
    ui.execute("CNGVIEW,-30,10,-120")
    st = plot_state(ui)
    d = st.plot_data(st.active, ui, raw=True)
    R = d["camera"]["basis"]
    P = S.to_screen(d["scene"]["xyz"], R, d["camera"]["center"])
    with np.errstate(all="ignore"):
        assert np.allclose(P, (d["scene"]["xyz"] - d["camera"]["center"]) @ R.T)
    fig = model_figure(d)
    fig.canvas.draw()
    ax = [a for a in fig.axes if a.name == "3d"][0]
    M = ax.get_proj()
    disp = []
    for p in P[:6]:
        x, y, _ = proj3d.proj_transform(p[0], p[1], p[2], M)
        disp.append(ax.transData.transform((x, y)))
    disp = np.asarray(disp)
    # display coordinates are an affine image of the screen coordinates (x right, y up, z ignored)
    A = np.column_stack([P[:6, :2], np.ones(6)])
    coef, res, *_ = np.linalg.lstsq(A, disp, rcond=None)
    assert np.allclose(A @ coef, disp, atol=1e-6)
    assert coef[0, 0] > 0 and abs(coef[1, 0]) < 1e-9 and coef[1, 1] > 0 and abs(coef[0, 1]) < 1e-9


def test_manual_graphing_macro_end_to_end(tmp_path):
    """Manual 5.6.2: READSPEC x3 -> SRSS -> WRITESPEC -> SPECPLOT -> CAPTUREPLOT -> CLOSEPLOT in batch."""
    f = np.linspace(0.1, 50, 301)
    for k, d in enumerate("xyz"):
        write_xy(tmp_path / f"Node17{d}.rs", f, 0.2 + 0.1 * (k + 1) * np.exp(-((f - 5 * (k + 1)) / 2) ** 2))
    (tmp_path / "SRSS-macro.pre").write_text("READSPEC,$1$,1,1\nREADSPEC,$2$,1,2\nREADSPEC,$3$,1,3\nSRSS,4,1,2,3\n"
                                             "WRITESPEC,$4$,4\nSPECPLOT,1,2,3,4\nCAPTUREPLOT,$4$.png\nCLOSEPLOT\n")
    ui = Interpreter(cwd=tmp_path)
    s = ui.run_text("LOADMACRO,SRSS,SRSS-macro.pre\nMACRO,SRSS,Node17x.rs,Node17y.rs,Node17z.rs,Node17srss.rs\n")
    assert s.errors == 0, errors(ui)
    assert (tmp_path / "Node17srss.rs.png").read_bytes()[:4] == b"\x89PNG"
    assert (tmp_path / "sassi_plot01_specplot.png").exists()          # batch auto-render
    st = plot_state(ui)
    assert st.active is None and st.lines[4].name == "SRSS Line"
    from sassi.plotting.lines import read_spec_file
    x, ys = read_spec_file(tmp_path / "Node17srss.rs", 1)
    assert np.allclose(ys[0], np.sqrt(sum(st.lines[k].y ** 2 for k in (1, 2, 3))))


def test_view_commands_count_view_revisions(model_ui):
    """Every view command increments view_rev (also RSTVIEW of an unchanged view), so the browser replaces a
    view rotated or zoomed with the mouse by the commanded one; display toggles do not."""
    ui = model_ui
    assert ui.execute("MODELPLOT")
    v = plot_state(ui).active_plot.view
    assert v.view_rev == 0
    for k, line in enumerate(["CNGVIEW,10,20,30", "RSTVIEW", "RSTVIEW", "CNGCENTER,1,2,3", "RSTCENTER"], start=1):
        assert ui.execute(line)
        assert v.view_rev == k, line
    for line in ("NODENUM,1", "ELENUM,1", "WIREFRAME,1", "SHRINK"):
        assert ui.execute(line)
    assert v.view_rev == 5 and v.to_dict()["view_rev"] == 5
