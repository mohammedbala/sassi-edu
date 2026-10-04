"""ETYPEGEN, INTGEN and RADIUS (requirements 3.4.I; spec 09 sections 2.6, 2.19, 2.25; D-HOU-01,
D-PNT-05, D-MDL-07).  The closed-form INTGEN counts on boxes are VP-50."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.prep import Interpreter, Kind
from sassi.prep import generation_lib as gl
from sassi.prep.writer import write_pre
from sassi.verify.problems.vp_generation import box_lines, box_node, box_reference, plate_lines


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def run(ui, lines):
    ui.run_text("\n".join(lines) if isinstance(lines, list) else lines)
    return ui.model


def n_int(m):
    return sum(1 for nd in m.nodes.values() if 0 in nd.flags)


def int_set(m):
    return {i for i, nd in m.nodes.items() if 0 in nd.flags}


def roundtrip(m):
    text, _ = write_pre(m)
    b = Interpreter()
    b.run_text(text)
    assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)
    return b.model


# ======================================================================================
# ETYPEGEN
# ======================================================================================
def _mixed(ui):
    """2 SOLIDs below grade, 1 SOLID straddling grade, 1 SOLID above, 1 SHELL (buried), 1 BEAM."""
    lines = []
    for k, z in enumerate((-10, -5, 0, 5, 10)):
        for j in range(2):
            for i in range(2):
                lines.append(f"N,{1 + i + 2 * j + 4 * k},{5 * i},{5 * j},{z}")
    lines += ["N,100,20,0,0", "N,101,25,0,0", "N,102,20,5,0", "GROUP,1,SOLID"]
    for e in range(4):
        b = 1 + 4 * e
        lines.append(f"E,{e + 1},{b},{b + 1},{b + 3},{b + 2},{b + 4},{b + 5},{b + 7},{b + 6}")
    lines += ["GROUP,2,SHELL", "E,1,1,2,4,3", "ETYPE,1,1,1,2", "GROUP,3,BEAMS", "E,1,100,101,102",
              "ETYPE,1,1,1,1", "GROUNDELEV,0"]
    return run(ui, lines)


def test_etypegen_1_and_2_set_solid_and_shell_types_but_not_beams(ui):
    m = _mixed(ui)
    m.groups[3].elements[1].etype = 0
    ui.execute("ETYPEGEN,1")
    assert all(e.etype == 1 for e in m.groups[1].elements.values())
    assert m.groups[2].elements[1].etype == 1
    assert m.groups[3].elements[1].etype == 0                     # beams untouched
    assert any("buried shells (ETYPE 2) are now structural" in w for w in ui.sink.texts(Kind.WARNING))
    ui.execute("ETYPEGEN,2")
    assert all(e.etype == 2 for e in m.groups[1].elements.values()) and m.groups[2].elements[1].etype == 2
    assert any("now buried shells" in w for w in ui.sink.texts(Kind.WARNING))


def test_etypegen_0_resolves_by_location_and_reset(ui):
    m = _mixed(ui)
    ui.execute("ETYPEGEN,0")
    # D-HOU-01: elements 1-2 below grade (all nodes <= gelev, centroid below); 3-4 above or straddling
    assert [m.groups[1].elements[e].etype for e in (1, 2, 3, 4)] == [2, 2, 1, 1]
    assert m.groups[2].elements[1].etype == 1                     # implicit shells are structure
    ui.execute("ETYPEGEN,0,RESET")
    assert all(e.etype == 0 for g in (1, 2) for e in m.groups[g].elements.values())
    assert not ui.execute("ETYPEGEN,3")
    assert not ui.execute("ETYPEGEN,1,RESET")


# ======================================================================================
# INTGEN
# ======================================================================================
def test_ffv_level_selection_rule(ui):
    m = run(ui, box_lines(3, 3, 6))
    v = gl.ModelView(m)
    nodes = v.excavation()["nodes"]
    for skip, sel in ((0, [1, 2, 3, 4, 5]), (1, [2, 4]), (2, [3]), (5, [])):
        means, level, s = gl.ffv_levels(v, nodes, skip)
        assert len(means) == 7 and s == sel
    ui.execute("INTGEN,5,1")
    sel = {box_node(3, 3, i, j, k) for k in (2, 4) for j in range(4) for i in range(4)}
    evbn = {box_node(3, 3, i, j, k) for k in range(7) for j in range(4) for i in range(4)
            if i in (0, 3) or j in (0, 3) or k in (0, 6)}
    assert int_set(m) == sel | evbn


def test_intgen_sets_the_same_flags_as_int_and_round_trips(ui):
    m = run(ui, box_lines(2, 2, 2))
    ui.execute("INTGEN,3")
    ref = Interpreter()
    ref.run_text("\n".join(box_lines(2, 2, 2) + [f"INT,{n},{n},1,1,0" for n in sorted(int_set(m))]))
    assert m.same_state(ref.model)
    assert roundtrip(m).same_state(m)
    # other INT codes are untouched by options 1-5 and by option 0
    ui.execute("INT,14,14,1,1,3")
    ui.execute("INTGEN,0")
    assert n_int(m) == 0 and m.nodes[14].flags == {3}


def test_intgen_2d_plane_excavation(ui):
    nx, nz = 4, 3
    lines = [f"N,{1 + i + (nx + 1) * k},{2.0 * i},0,{-2.0 * nz + 2.0 * k}" for k in range(nz + 1) for i in range(nx + 1)]
    lines.append("GROUP,1,PLANE")
    e = 0
    for k in range(nz):
        for i in range(nx):
            e += 1
            a = 1 + i + (nx + 1) * k
            lines.append(f"E,{e},{a},{a + 1},{a + nx + 2},{a + nx + 1}")
    lines += [f"ETYPE,1,{e},1,2", "GROUNDELEV,0"]
    m = run(ui, lines)
    fv = (nx + 1) * (nz + 1)
    evbn = fv - (nx - 1) * (nz - 1)
    for opt, ref in ((1, fv), (2, evbn), (3, evbn - (nx - 1)), (4, nx + 1)):
        ui.execute("INTGEN,0")
        assert ui.execute(f"INTGEN,{opt}")
        assert n_int(m) == ref, opt


def test_intgen_diagnostics_and_argument_checks(ui):
    m = run(ui, box_lines(2, 2, 2) + ["L,1,4,0.12,1500,800,0.05,0.05", "L,2,20,0.12,1500,800,0.05,0.05", "TOPL,1,2",
                                      "D,1,1,1,1,UX"])
    ui.execute("INTGEN,1")
    w = ui.sink.texts(Kind.WARNING)
    # levels at depth 0, 5, 10; interfaces at 0, 4, 24 -> the 18 nodes at depth 5 and 10 are off-interface
    assert any("18 interaction nodes are not on a soil-layer interface" in t for t in w)
    assert any("fixed translations" in t and "1" in t for t in w)
    assert not ui.execute("INTGEN,5,-1")
    assert not ui.execute("INTGEN,7")
    ui.execute("INTGEN,2,3")
    assert any("[level skip] is used by option 5" in t for t in ui.sink.texts(Kind.WARNING))


def test_intgen_numbering_warning_and_implicit_etype(ui):
    # nodes numbered top-down: warning about bottom-up numbering
    m = run(ui, box_lines(1, 1, 1, etype=0))
    ui.execute("INTGEN,1")
    w = ui.sink.texts(Kind.WARNING)
    assert any("implicit ETYPE 0" in t for t in w)
    assert not any("bottom-up" in t for t in w)
    ui2 = Interpreter()
    ui2.run_text("N,1,0,0,0\nN,2,1,0,0\nN,3,1,1,0\nN,4,0,1,0\nN,5,0,0,-1\nN,6,1,0,-1\nN,7,1,1,-1\nN,8,0,1,-1\n"
                 "GROUP,1,SOLID\nE,1,5,6,7,8,1,2,3,4\nETYPE,1,1,1,2\nINTGEN,1")
    assert n_int(ui2.model) == 8
    assert any("bottom-up" in t for t in ui2.sink.texts(Kind.WARNING))


def test_intgen_surface_needs_nodes_at_grade(ui):
    run(ui, plate_lines(2, 2, z=-3.0) + ["GROUNDELEV,0"])
    assert not ui.execute("INTGEN,4")
    assert n_int(ui.model) == 0


def test_intgen_surface_beam_grid_and_fixrot_ground_nodes(ui):
    # rigid beam grid at grade (K nodes above grade), a column, and a FIXROT plate: K nodes, the column top
    # and the fixed ground nodes of the drilling springs are not interaction nodes
    m = run(ui, ["N,1,0,0,0", "N,2,5,0,0", "N,3,0,5,0", "N,4,0,0,0.0", "N,10,0,0,9", "N,11,1,0,9", "GROUP,1,BEAMS",
                 "E,1,1,2,10", "E,2,1,3,10", "E,3,1,10,11", "N,20,10,0,0", "N,21,15,0,0", "N,22,15,5,0",
                 "N,23,10,5,0", "GROUP,2,SHELL", "E,1,20,21,22,23", "ROTATE,12.5,2.5,0,0,0,0", "GROUNDELEV,0"])
    ui.execute("FIXSHLROT")
    ui.execute("INTGEN,4")
    assert int_set(m) == {1, 2, 3, 20, 21, 22, 23}


# ======================================================================================
# RADIUS
# ======================================================================================
def test_radius_uniform_mesh_is_0_9_h(ui, tmp_path):
    m = run(ui, box_lines(3, 2, 2, h=6.0))
    ui.execute("RADIUS,,radii.txt")
    txt = (tmp_path / "radii.txt").read_text().splitlines()
    rows = [ln.split() for ln in txt if ln and not ln.startswith("#")]
    data = [r for r in rows if r[0].isdigit()]
    assert len(data) == 12 and all(math.isclose(float(r[2]), 0.9 * 6.0) for r in data)
    assert rows[-2] == ["AVERAGE", "5.4"]
    assert math.isclose(ui.session["radius"]["average"], 5.4)


def test_radius_non_uniform_and_scale(ui):
    # 2 elements of 4 x 4 and 2 x 2 in plan
    lines = ["N,1,0,0,-2", "N,2,4,0,-2", "N,3,4,4,-2", "N,4,0,4,-2", "N,5,0,0,0", "N,6,4,0,0", "N,7,4,4,0", "N,8,0,4,0",
             "N,11,10,0,-2", "N,12,12,0,-2", "N,13,12,2,-2", "N,14,10,2,-2", "N,15,10,0,0", "N,16,12,0,0",
             "N,17,12,2,0", "N,18,10,2,0", "GROUP,1,SOLID", "E,1,1,2,3,4,5,6,7,8", "E,2,11,12,13,14,15,16,17,18",
             "ETYPE,1,2,1,2"]
    m = run(ui, lines)
    r = gl.excavation_radii(m, 1.0)
    assert [x for _, _, x in r] == [4.0, 2.0]
    ui.execute("RADIUS,0.85")
    assert ui.session["radius"]["min"] == pytest.approx(1.7) and ui.session["radius"]["max"] == pytest.approx(3.4)
    assert ui.session["radius"]["average"] == pytest.approx(2.55)
    assert not ui.execute("RADIUS,-1")


def test_radius_prism_uses_plan_hull_and_plane_uses_width(ui):
    m = run(ui, ["N,1,0,0,-1", "N,2,2,0,-1", "N,3,0,2,-1", "N,5,0,0,0", "N,6,2,0,0", "N,7,0,2,0", "GROUP,1,SOLID",
                 "E,1,1,2,3,3,5,6,7,7", "ETYPE,1,1,1,2"])
    assert gl.excavation_radii(m, 1.0)[0][2] == pytest.approx(math.sqrt(2.0))
    ui2 = Interpreter()
    ui2.run_text("N,1,0,0,-2\nN,2,3,0,-2\nN,3,3,0,0\nN,4,0,0,0\nGROUP,1,PLANE\nE,1,1,2,3,4\nETYPE,1,1,1,2")
    assert gl.excavation_radii(ui2.model, 1.0)[0][2] == pytest.approx(3.0)


def test_radius_without_excavation_fails(ui):
    run(ui, box_lines(1, 1, 1, etype=1))
    assert not ui.execute("RADIUS")
