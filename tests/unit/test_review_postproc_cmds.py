"""Adversarial review tests of the post-processing / plotting work package (requirements 3.4.L, 4.13, 5.7,
5.8; D-LIN-01..04, D-UI-10..17, D-FIL-02/03; spec 06 section 12; spec 09 sections 2.4, 2.11, 2.12; spec 10
section 4).

The checks are derived independently of the implementation: line maths against ``numpy.interp`` on the
union grid, BROADEN against a dense brute-force evaluation of the D-LIN-01 definition, the ACS discrete
grid by hand, FRAMESEL / CRITFREQ against straightforward loops written from the spec text, and the
closed-form cases of spec 06 section 12.

Tests marked ``xfail(strict=True)`` document defects found in the review (see the reason string); they
turn into XPASS failures once the defect is fixed, so the marker must then be removed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from sassi.plotting import lines as LM
from sassi.plotting import state as S
from sassi.plotting.lines import Line, LineError
from sassi.plotting.state import plot_state, write_frame
from sassi.prep import Interpreter, Kind


# ======================================================================================
# helpers
# ======================================================================================
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


def write_xy(path: Path, x, y) -> None:
    path.write_text("\n".join(f"{float(a)!r} {float(b)!r}" for a, b in zip(x, y)) + "\n")


def union_ref(lines_xy):
    """Independent union grid (exact values; no near-duplicates in the test data) and resampled rows."""
    X = np.unique(np.concatenate([np.asarray(x, float) for x, _ in lines_xy]))
    Y = np.vstack([np.interp(X, x, y, left=y[0], right=y[-1]) for x, y in lines_xy])
    return X, Y


def brute_broaden(X, E, b, f):
    """max{E(f') : f/(1+b) <= f' <= f/(1-b)} with E piecewise linear on X (constant outside),
    evaluated exactly at the window ends and at every vertex strictly inside the window."""
    out = np.empty(len(f))
    for i, fi in enumerate(f):
        lo, hi = fi / (1.0 + b), fi / (1.0 - b)
        pts = np.concatenate([[lo, hi], X[(X > lo) & (X < hi)]])
        out[i] = np.max(np.interp(pts, X, E, left=E[0], right=E[-1]))
    return out


def box_model(ui, n=2):
    """(n+1)^3 nodes of a unit grid, n^3 hexes (group 1 SOLID)."""
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
    ui.run_text("\n".join(out))
    assert not errors(ui)


# ======================================================================================
# Line maths on the union grid (requirements 4.13; spec 10 1.3, 4.1)
# ======================================================================================
def test_line_maths_commands_match_independent_union_grid(ui, tmp_path):
    rng = np.random.default_rng(101)
    data = []
    for k in range(3):
        x = np.sort(rng.choice(np.arange(1, 400) * 0.0625, size=12 + 5 * k, replace=False))   # exact binary
        y = rng.normal(size=len(x))
        write_xy(tmp_path / f"l{k}.rs", x, y)
        data.append((x, y))
    ui.run_text("READSPEC,l0.rs,1,1\nREADSPEC,l1.rs,1,2\nREADSPEC,l2.rs,1,3\n"
                "ADDITION,10,1,2,3\nSUBTRACTION,11,1,2,3\nAVERAGE,12,1,2,3\nSRSS,13,1,2,3\n"
                "LINECOMBIN,14,1,1.0,2,0.4,3,-0.4\n")
    assert not errors(ui)
    X, Y = union_ref(data)
    L = plot_state(ui).lines
    refs = {10: Y.sum(0), 11: Y[0] - Y[1] - Y[2], 12: Y.mean(0), 13: np.sqrt((Y ** 2).sum(0)),
            14: Y[0] + 0.4 * Y[1] - 0.4 * Y[2]}
    for num, ref in refs.items():
        assert np.array_equal(L[num].x, X)
        np.testing.assert_allclose(L[num].y, ref, rtol=0, atol=1e-13)
    assert L[10].name == L[11].name == L[14].name == "Linear Combin."
    assert L[12].name == "Average Line" and L[13].name == "SRSS Line"


def test_destination_may_be_a_source_and_undefined_source_is_an_error(ui, tmp_path):
    write_xy(tmp_path / "a.rs", [0, 1, 2], [0, 1, 2])
    write_xy(tmp_path / "b.rs", [0.5, 1.5], [10, 20])
    ui.run_text("READSPEC,a.rs,1,1\nREADSPEC,b.rs,1,2\n")
    assert ui.execute("ADDITION,1,1,2")                       # computed first, then overwritten
    assert plot_state(ui).lines[1].y.tolist() == [10, 10.5, 16, 21.5, 22]
    assert not ui.execute("ADDITION,5,1,99")                  # D-LIN-02: undefined source -> error
    assert 5 not in plot_state(ui).lines
    assert not ui.execute("LINECOMBIN,6,1,2.0,2")             # odd count after Dest
    assert not ui.execute("SRSS,7," + ",".join(["1"] * 101))  # more than 100 sources


# ======================================================================================
# BROADEN (D-LIN-01; spec 10 4.1.3)
# ======================================================================================
def test_broaden_multiline_exact_against_definition():
    """Envelope of three lines with different grids, +-12 % broadening: the piecewise-linear result
    equals the definition at 6000 points (between grid points too).  Independent evaluation: the
    maximum of the upper envelope over a window is reached at a window end or at a vertex of one of
    the lines (a crossing of two lines is a convex corner, never a maximum), so
    ``B(f) = max_i max{L_i(p) : p in {lo, hi} U (vertices inside (lo, hi))}``."""
    rng = np.random.default_rng(7)
    srcs = []
    for k in range(3):
        x = np.unique(np.round(rng.uniform(0.2, 40.0, 25 + 7 * k), 3))
        srcs.append(Line(k + 1, "", x, rng.uniform(0.1, 3.0, len(x))))
    x, y, notes = LM.broaden(srcs, 0.0, 12.0)
    assert notes and "EDU-15" in notes[0]                     # fewer than 301 points
    Xall = np.unique(np.concatenate([L.x for L in srcs]))
    b = 0.12
    f = np.unique(np.concatenate([np.linspace(Xall[0], Xall[-1], 6000), x]))
    ref = np.empty(len(f))
    for i, fi in enumerate(f):
        lo, hi = fi / (1 + b), fi / (1 - b)
        pts = np.concatenate([[lo, hi], Xall[(Xall > lo) & (Xall < hi)]])
        ref[i] = max(float(np.max(L.at(pts))) for L in srcs)
    got = np.interp(f, x, y)
    assert np.max(np.abs(got - ref)) < 1e-9


def test_broaden_acs_grid_is_the_discrete_form_on_the_source_points(ui, tmp_path):
    """EDUOPT,BROADENGRID,ACS: B(f_k) = max{E(f_j): f_j in [f_k/(1+b), f_k/(1-b)]} on the source
    grid only (T-B1 data, by hand: 9 -> {9,10}, 10 -> {9,10,11}, 11 -> {10,11} all reach 5)."""
    write_xy(tmp_path / "tb1.rs", [1, 5, 9, 10, 11, 15, 20], [1, 1, 1, 5, 1, 1, 1])
    ui.run_text("READSPEC,tb1.rs,1,1\nEDUOPT,BROADENGRID,ACS\nBROADEN,2,0,15,1\n")
    assert not errors(ui)
    B = plot_state(ui).lines[2]
    assert B.x.tolist() == [1, 5, 9, 10, 11, 15, 20]
    assert B.y.tolist() == [1, 1, 5, 5, 5, 1, 1]
    assert any("301" in w for w in warnings(ui))              # EDU-15


def test_tb2_amplitude_criterion_bridges_at_10_percent(ui, tmp_path):
    """Spec 10 T-B2: with Smooth1 = 10 the valley criterion leaves the line unchanged but the
    peak-amplitude criterion (EDUOPT,BRIDGE,AMPLITUDE) bridges: flat 3.8 from 5.2 to 8."""
    xb = [1, 4, 5, 6, 6.5, 7, 8, 9, 12, 20]
    yb = [0.5, 2, 4, 3, 3.5, 3.0, 3.8, 2, 1, 0.5]
    write_xy(tmp_path / "tb2.rs", xb, yb)
    ui.run_text("READSPEC,tb2.rs,1,1\nEDUOPT,BRIDGE,AMPLITUDE\nBROADEN,2,10,0,1\n")
    assert not errors(ui)
    B = plot_state(ui).lines[2]
    for xv in (5.2, 5.6, 6.0, 6.5, 7.0, 7.9, 8.0):
        assert B.at([xv])[0] == pytest.approx(3.8, abs=1e-12)
    for xv, yv in ((1, 0.5), (4, 2.0), (5, 4.0), (5.1, 3.9), (9, 2.0), (20, 0.5)):
        assert B.at([xv])[0] == pytest.approx(yv, abs=1e-12)


def test_bridging_is_idempotent_and_never_lowers_the_broadened_line():
    rng = np.random.default_rng(3)
    for _ in range(150):
        x = np.unique(np.round(rng.uniform(0.5, 30.0, 25), 2))
        y = rng.uniform(0.5, 4.0, len(x))                     # continuous values (no exact ties)
        L = Line(1, "", x, y)
        xe, e = LM.envelope([L])
        xg, g = LM.broaden_peaks(xe, e, 0.1)
        xr, r = LM.bridge_peaks(xg, g, 0.2)
        assert np.all(np.diff(xr) > 0)
        f = np.linspace(x[0], x[-1], 3001)
        assert np.min(np.interp(f, xr, r) - np.interp(f, xg, g)) > -1e-12
        x2, r2 = LM.bridge_peaks(xr, r, 0.2)
        assert np.max(np.abs(np.interp(f, x2, r2) - np.interp(f, xr, r))) < 1e-12


def test_broaden_with_equal_peak_values_does_not_fail(ui, tmp_path):
    x = [5.27, 6.03, 9.59, 10.11, 10.55, 12.97, 13.06, 15.94, 19.46, 21.24, 24.69]
    y = [1.0, 5.0, 1.0, 1.0, 1.0, 3.0, 5.0, 4.0, 3.0, 1.0, 4.0]
    write_xy(tmp_path / "eq.rs", x, y)
    ui.run_text("READSPEC,eq.rs,1,1\n")
    ok = ui.execute("BROADEN,2,15,15,1")
    assert ok, errors(ui)
    B = plot_state(ui).lines[2]
    assert np.all(np.diff(B.x) > 0)
    assert np.min(B.at(np.asarray(x)) - np.asarray(y)) >= -1e-12


# ======================================================================================
# Line files (D-FIL-02, D-LIN-04; spec 10 4.2; spec 06 12.5-12.6)
# ======================================================================================
def test_readth_pair0_several_values_per_line_and_writeth_warning(ui, tmp_path):
    (tmp_path / "h.acc").write_text("# header\n0.01\n1 2 3\n4, 5\n")
    assert ui.execute("READTH,h.acc,0,1")
    L = plot_state(ui).lines[1]
    np.testing.assert_allclose(L.x, [0.0, 0.01, 0.02, 0.03, 0.04], rtol=0, atol=1e-15)
    assert L.y.tolist() == [1, 2, 3, 4, 5] and L.kind == "history"
    write_xy(tmp_path / "nu.th", [0.0, 0.01, 0.025, 0.03], [1, 2, 3, 4])
    ui.run_text("READTH,nu.th,1,2\n")
    n0 = len(warnings(ui))
    assert ui.execute("WRITETH,nu_out.acc,2")
    assert len(warnings(ui)) > n0                             # non-uniform step -> warning
    first = (tmp_path / "nu_out.acc").read_text().split()[0]
    assert float(first) == pytest.approx(0.01)                # dt = x2 - x1


def test_readspec_numlines_counting_and_names(ui, tmp_path):
    (tmp_path / "c.tfu").write_text("# f a p\n1 10 0.1\n2 20 0.2\n3 30 0.3\n")
    assert ui.execute("READSPEC,c.tfu,1,4")                   # frequency column not counted
    assert plot_state(ui).lines[4].y.tolist() == [10, 20, 30] and plot_state(ui).lines[4].name == "c.tfu"
    assert ui.execute("READSPEC,c.tfu,2,5,6")
    assert plot_state(ui).lines[6].y.tolist() == [0.1, 0.2, 0.3] and plot_state(ui).lines[6].name == "c.tfu[2]"
    assert not ui.execute("READSPEC,c.tfu,3,7,8,9")           # only 2 data columns


# ======================================================================================
# 2D plot requests (spec 06 12.13; spec 10 4.3)
# ======================================================================================
def test_specplot_minus_one_ends_the_list_and_unknown_lines_are_ignored(ui, tmp_path):
    for k in (1, 2, 3):
        write_xy(tmp_path / f"s{k}.rs", [1, 2], [k, k])
        ui.execute(f"READSPEC,s{k}.rs,1,{k}")
    assert ui.execute("SPECPLOT,1,2,-1,3")
    assert plot_state(ui).active_plot.params["lines"] == [1, 2]
    assert ui.execute("SPECPLOT,1,77,3")
    assert plot_state(ui).active_plot.params["lines"] == [1, 3]
    assert any("77" in w for w in warnings(ui))


def test_plotrange_validation_on_log_axis(ui, tmp_path):
    write_xy(tmp_path / "s.rs", [0.1, 1, 10], [1, 2, 3])
    ui.run_text("READSPEC,s.rs,1,1\nSPECPLOT,1\nAXES,1,1,0,0,1,0\n")
    assert not ui.execute("PLOTRANGE,0,10,0,5")               # log X: minimum must be > 0
    assert not ui.execute("PLOTRANGE,5,1,0,5")                # min >= max
    assert ui.execute("PLOTRANGE,0.1,100,0,5")
    d = plot_state(ui).plot_data(plot_state(ui).active, ui)
    assert d["extent"] == [0.1, 100.0, 0.0, 5.0]


def test_captureplot_format_decided_by_the_file_name_only(tmp_path, monkeypatch):
    monkeypatch.setenv("SASSI_EDU_SETTINGS_DIR", str(tmp_path / "settings"))
    wd = tmp_path / "plots.png.d"
    wd.mkdir()
    u = Interpreter(cwd=wd)
    plot_state(u).auto_render = False
    write_xy(wd / "s.rs", [1, 2, 3], [1, 2, 1])
    u.run_text("READSPEC,s.rs,1,1\nSPECPLOT,1\n")
    assert u.execute("CAPTUREPLOT,out.bmp")
    assert (wd / "out.bmp").read_bytes()[:2] == b"BM"


# ======================================================================================
# 3D scenes (spec 06 sections 1-4, 12)
# ======================================================================================
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
N,16,0,0,5
GROUP,1,SOLID
E,1,1,2,3,4,5,6,7,8
MACT,2
E,2,2,9,10,3,6,11,12,7
GROUP,2,BEAMS
E,1,7,13,14
INT,1,4,1,1
D,2,2,1,1,UZ
MT,13,0,0,5
"""


def test_model_scene_connectivity_faces_markers(ui):
    ui.run_text(MODEL)
    assert not errors(ui)
    sc = S.model_scene(ui.model, color_by=2, show_dof=[2], show_mass=True)
    ids = sc["node_id"].tolist()
    # element-connected DOF nodes only: K node 14 (beam orientation) and orphan 16 are excluded
    assert ids == list(range(1, 14))
    # two hexes sharing one face: 12 faces, 2 interior (shared), 10 boundary
    assert len(sc["faces"]) == 12 and int(sc["face_boundary"].sum()) == 10
    assert len(sc["edges"]) == 1
    # ELECOLOR 2: key = material number of the solids (MACT 2 for element 2)
    key = dict(zip(zip(sc["elem_group"].tolist(), sc["elem_id"].tolist()), sc["elem_key"].tolist()))
    assert key[(1, 1)] == 1 and key[(1, 2)] == 2
    assert [ids[i] for i in sc["interaction"]] == [1, 2, 3, 4]
    assert [ids[i] for i in sc["fixed"]] == [2]               # UZ fixed at node 2 only
    assert [ids[i] for i in sc["mass"]] == [13] and sc["mass_dir"].tolist() == [[False, False, True]]
    # RSTCENTER: centre of the bounding box of element-connected nodes
    np.testing.assert_allclose(sc["center"], [1.0, 0.5, 1.5])


def test_display_volume_hides_nodes_of_clipped_elements_on_node_plot(ui):
    box_model(ui, 4)
    assert ui.execute("WINDOWSETTINGS,VOLUME,0,0.5,,,,")
    assert ui.execute("NODEPLOT")
    d = plot_state(ui).plot_data(plot_state(ui).active, ui, raw=True)
    sc = d["scene"]
    vis = np.asarray(sc["node_visible"], bool)
    x = sc["xyz"][:, 0]
    assert sc["n_elements"] == 16                             # elements with a node inside x <= 0.5
    # nodes at x >= 2 belong to no displayed element and lie outside the box: they must not be shown
    assert int(np.sum(vis & (x >= 2.0))) == 0


def test_camera_projection_and_default_view(ui):
    ui.run_text(MODEL)
    assert ui.execute("MODELPLOT")
    d = plot_state(ui).plot_data(plot_state(ui).active, ui, raw=True)
    R = np.asarray(d["camera"]["basis"])
    # D-UI-12 default isometric (rX -60, rZ -45): Z points up on the screen, R orthonormal
    np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-14)
    assert (R @ [0, 0, 1])[1] == pytest.approx(np.sin(np.radians(60)), abs=1e-14)
    ui.execute("CNGVIEW,0,0,90,,,2")
    d = plot_state(ui).plot_data(plot_state(ui).active, ui, raw=True)
    R = np.asarray(d["camera"]["basis"])
    np.testing.assert_allclose(R @ [1, 0, 0], [0, 1, 0], atol=1e-14)     # Rz(90): model X -> screen up
    half = d["camera"]["half"]
    s = S.to_screen(d["scene"]["xyz"], R, d["camera"]["center"])
    assert half == pytest.approx(1.08 * np.max(np.abs(s[:, :2])) / 2.0)  # zoom 2 halves the window
    ui.execute("RSTVIEW")
    assert plot_state(ui).active_plot.view.zoom == 1.0


# ======================================================================================
# Layer table (spec 06 6; spec 06 12.14)
# ======================================================================================
def test_layer_table_halfspace_depth_independent_of_end_layer(ui):
    ui.run_text("L,1,10,120,2000,1000,0.02,0.03\nL,2,20,125,2500,1200,0.02,0.02\nL,3,0,130,4000,2000,0.01,0.01\n"
                "TOPL,1,2\nSITE,0,1,0,20,3\n")
    full = S.layer_table(ui.model)
    assert full["halfspace"]["top"] == 30.0
    part = S.layer_table(ui.model, 1, 1)
    assert [r["layer"] for r in part["layers"]] == [1]
    assert part["halfspace"]["top"] == 30.0


def test_layer_table_end_minus_one_lists_all_layers_plus_halfspace(ui):
    ui.run_text("\n".join(f"L,{k},{k},120,2000,1000,0.02,0.02" for k in range(1, 9)) +
                "\nTOPL,1,2,3,4,5,6,7\nSITE,0,1,0,20,8\n")
    t = S.layer_table(ui.model, 1, -1)
    assert [r["layer"] for r in t["layers"]] == [1, 2, 3, 4, 5, 6, 7]
    assert t["halfspace"]["layer"] == 8
    assert [r["top"] for r in t["layers"]] == [0, 1, 3, 6, 10, 15, 21]


# ======================================================================================
# Animation maths (spec 06 12.8-12.11; D-UI-13, D-UI-14)
# ======================================================================================
def test_animation_closed_forms():
    # spec 06 12.9 vector rule (X), D-UI-13 for Y and Z
    re = np.array([[5.0, 2.0, -1.0]])
    im = np.array([[0.3, 0.7, 0.2]])
    V = S.vector_rule(re, im, 1.0)
    np.testing.assert_allclose(V[0, 0], [5.0, 0.3, 0.3])
    np.testing.assert_allclose(V[0, 1], [0.7, 2.0, 0.7])
    np.testing.assert_allclose(V[0, 2], [0.2, 0.2, -1.0])
    # spec 06 12.10 deformed shape
    np.testing.assert_allclose(S.deformed_coordinates([[1.0, 2.0, 3.0]], [[0.01, 0.0, 0.0]], 10.0), [[1.1, 2, 3]])
    # spec 06 12.11 frame selection
    assert len(S.frame_sequence(1, 3110, 1, 3110)) == 3110
    assert S.frame_sequence(1, 10, 3, 10) == [1, 4, 7, 10]
    # spec 06 12.8 colour map: below / above the range clamp, midpoint green; D-UI-14 bubble size
    t = S.colormap_position(np.array([-11.0, 0.0, 11.0]), -10.0, 10.0)
    np.testing.assert_allclose(t, [0.0, 0.5, 1.0])
    np.testing.assert_allclose(S.jet(0.5), [0.5, 1.0, 0.5])
    np.testing.assert_allclose(S.bubble_size(np.array([0.0, 0.1, 0.6]), 10.0), [2.0, 2.0, 6.0])
    assert S.colorbar_labels(-10.0, 10.0) == ["10.00000", "5.00000", "0.00000", "-5.00000", "-10.00000"]


def test_frame_directory_is_ordered_by_frame_number(tmp_path):
    d = tmp_path / "THD"
    d.mkdir()
    dt = 5.0
    for s in range(25):                                       # t = 0 .. 120 s
        write_frame(d / f"THD_{s * dt:06.3f}_{s + 1:05d}", [1], [[float(s), 0.0, 0.0]])
    files = S.read_frame_list(d)
    numbers = [S.frame_name_info(f.name)["frame"] for f in files]
    assert numbers == list(range(1, 26))


# ======================================================================================
# FRAMESEL, CRITFREQ, FRAMECOMBIN (spec 09 2.4, 2.11, 2.12)
# ======================================================================================
def test_framesel_matches_a_direct_loop(ui, tmp_path):
    rng = np.random.default_rng(12)
    a = np.round(rng.normal(size=300), 6)
    a[0] = 0.0
    tol = 35.0
    amax = np.max(np.abs(a))
    ref = []
    for j in range(len(a)):
        nb = [a[j - 1]] if j == len(a) - 1 else ([a[j + 1]] if j == 0 else [a[j - 1], a[j + 1]])
        is_ext = all(a[j] >= v for v in nb) or all(a[j] <= v for v in nb)
        if is_ext and abs(a[j]) >= tol / 100 * amax:
            ref.append(j + 1)
    (tmp_path / "r.acc").write_text("0.01\n" + "\n".join(repr(float(v)) for v in a) + "\n")
    assert ui.execute(f"FRAMESEL,{tol},r.acc,SEL")
    assert ui.variables["sel"].items == [str(k) for k in ref]


def test_critfreq_motion_format_minfilter_and_frequency_numbers(ui, tmp_path):
    """TFU/TFI written by sassi.io.textfiles.write_tf (MOTION format, 6-decimal frequencies) with
    df = 1/(4096 * 0.005).  Two interpolation overshoots of 40 % over the larger bracketing computed
    value: a high one at fnum 1550 (3.5 vs 2.5) and a low one at fnum 350 (1.4 vs 1.0, i.e. 40 % of
    the TFI maximum).  minfilter 50 considers only the high one; minfilter 100 both."""
    from sassi.io.textfiles import write_tf
    df = 1.0 / (4096 * 0.005)
    fnum_u = np.arange(100, 2001, 100)
    fu = fnum_u * df
    au = np.full(len(fu), 1.0)
    au[fnum_u == 1500] = 2.5
    k = np.arange(0, 2049)
    fi = k * df
    ai = np.interp(fi, fu, au, left=0.0, right=0.0)
    for k0, h in ((1550, 1.75), (350, 0.4)):
        ai = ai + h * np.clip(1.0 - np.abs(k - k0) / 20.0, 0.0, None)
    write_tf(tmp_path / "T.TFU", fu, au.astype(complex), complex_=True)
    write_tf(tmp_path / "T.TFI", fi, ai.astype(complex), complex_=True)
    assert int(np.argmax(ai)) == 1550 and ai[1550] == pytest.approx(3.5)
    assert ai[350] == pytest.approx(1.4)
    assert ui.execute("CRITFREQ,20,50,T,HI")
    assert ui.variables["hi"].items == ["1550"]
    assert ui.execute("CRITFREQ,20,100,T,ALL")               # minfilter 100: every peak considered
    assert ui.variables["all"].items == ["350", "1550"]
    assert ui.execute("CRITFREQ,45,100,T,NONE")
    assert ui.variables["none"].items == []
    assert not errors(ui)


def test_framecombin_three_frames_against_numpy(ui, tmp_path):
    rng = np.random.default_rng(4)
    nodes = np.array([4, 9, 17])
    vals = [np.round(rng.normal(size=(3, 6)), 4) for _ in range(3)]
    for k, v in enumerate(vals):
        write_frame(tmp_path / f"TFU_{k}", nodes, v)
    V = np.stack(vals)
    for op, ref in ((0, np.sqrt((V ** 2).sum(0))), (1, V.sum(0)), (2, V.mean(0))):
        assert ui.execute(f"FRAMECOMBIN,{op},3,TFU_0,TFU_1,TFU_2,out{op}")
        fr = S.read_frame(tmp_path / f"out{op}")
        assert fr.nodes.tolist() == nodes.tolist()
        np.testing.assert_allclose(fr.values, ref, rtol=1e-10, atol=1e-12)
