"""Node and coordinate-system commands (spec 08 sections 2-3): UT-12, UT-13 (node part), generation,
deletion cascade, fixities, INT flags, listings."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.prep import Interpreter, Kind


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def run(ui, text):
    ui.run_text(text)
    return ui.model


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


def xyz(m, n):
    return np.array(m.nodes[n].xyz)


# ------------------------------------------------------------------ UT-12 coordinate systems
def test_ut12_loc_csys_global(ui):
    m = run(ui, "LOC,1,0,10,0,0,45,0,0\nCSYS,1\nN,100,5,0,0\n")
    assert m.nodes[100].csys == 1 and m.nodes[100].xyz == (5.0, 0.0, 0.0)   # stored in the local system
    g = m.node_global(100)
    assert g == pytest.approx([10 + 5 / math.sqrt(2), 5 / math.sqrt(2), 0.0], abs=1e-12)
    run(ui, "GLOBAL,100,100,1")
    assert m.nodes[100].csys == 0
    assert xyz(m, 100) == pytest.approx([13.5355339, 3.5355339, 0.0], abs=1e-6)
    assert m.csys_active == 1                         # GLOBAL does not deactivate the system
    run(ui, "N,101,1,0,0")
    assert m.nodes[101].csys == 1


def test_ut12_local_equals_loc_90(ui):
    m = run(ui, "N,1,0,0,0\nN,2,0,1,0\nN,3,-1,0,0\nLOCAL,2,0,1,2,3\nLOC,3,0,0,0,0,90,0,0\n")
    assert np.max(np.abs(m.csys[2].R - m.csys[3].R)) < 1e-12
    assert np.max(np.abs(m.csys[2].origin - m.csys[3].origin)) < 1e-12
    assert m.csys[2].R[:, 0] == pytest.approx([0, 1, 0])     # ex = +Y
    assert m.csys[2].R[:, 1] == pytest.approx([-1, 0, 0])    # ey = -X


def test_ut12_euler_order_local_y_to_global_z(ui):
    m = run(ui, "LOC,4,0,0,0,0,0,90,0\n")
    assert m.csys[4].R[:, 1] == pytest.approx([0, 0, 1], abs=1e-15)
    # full order check against R = Rz Rx Ry for arbitrary angles
    run(ui, "LOC,5,0,1,2,3,30,20,10\n")
    a, b, c = np.radians([30, 20, 10])
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(b), -np.sin(b)], [0, np.sin(b), np.cos(b)]])
    Ry = np.array([[np.cos(c), 0, np.sin(c)], [0, 1, 0], [-np.sin(c), 0, np.cos(c)]])
    assert np.allclose(m.csys[5].R, Rz @ Rx @ Ry, atol=1e-14)


def test_coordinate_system_errors(ui):
    run(ui, "CSYS,3\nLOC,1,1,0,0,0,0,0,0\nLOC,0,0,0,0,0,0,0,0\nN,1,0,0,0\nN,2,1,0,0\nN,3,2,0,0\n"
            "LOCAL,2,0,1,2,3\nLOCAL,2,0,1,2,99\n")
    e = errors(ui)
    assert "CSYS: coordinate system 3 is not defined" in e
    assert any("cylindrical" in t for t in e)
    assert any("system numbers are >= 1" in t for t in e)
    assert any("collinear" in t for t in e)
    assert any("nodes [99] are not defined" in t for t in e)


def test_loc_redefinition_keeps_local_coordinates(ui):
    m = run(ui, "LOC,1,0,0,0,0,0,0,0\nCSYS,1\nN,1,1,0,0\nLOC,1,0,5,0,0,0,0,0\n")
    assert any("keep their local coordinates" in w for w in warnings(ui))
    assert m.nodes[1].xyz == (1.0, 0.0, 0.0)
    assert m.node_global(1) == pytest.approx([6, 0, 0])


def test_sdel_converts_nodes_and_resets_active(ui):
    m = run(ui, "LOC,1,0,10,0,0,90,0,0\nCSYS,1\nN,1,1,0,0\nSDEL,1\n")
    assert 1 not in m.csys and m.csys_active == 0
    assert m.nodes[1].csys == 0 and xyz(m, 1) == pytest.approx([10, 1, 0], abs=1e-12)
    run(ui, "SDEL,0")
    assert any("cannot be deleted" in e for e in errors(ui))


def test_slist_lists(ui):
    run(ui, "N,1,0,0,0\nN,2,1,0,0\nN,3,0,1,0\nLOCAL,1,0,1,2,3\nLOC,2,0,0,0,0,10,0,0\nSLIST")
    info = ui.sink.texts(Kind.INFO)
    assert any("system 1" in t and "from nodes 1, 2, 3" in t for t in info)
    assert "2 systems listed" in info


# ------------------------------------------------------------------ UT-13 generation
def test_ut13_fill_and_ngen(ui):
    m = run(ui, "N,1,0,0,0\nN,5,40,0,0\nFILL\n")
    assert [xyz(m, n)[0] for n in (2, 3, 4)] == [10.0, 20.0, 30.0]
    run(ui, "NGEN,3,5,1,5,1,0,0,10\n")
    for k, z in enumerate((10.0, 20.0, 30.0), start=1):
        for n in range(1, 6):
            assert xyz(m, n + 5 * k) == pytest.approx([10.0 * (n - 1), 0, z], abs=0)
    assert sorted(m.nodes) == list(range(1, 21))
    # FILL without arguments right after NGEN uses the last two defined nodes (19, 20): nothing between
    run(ui, "FILL")
    assert "FILL: no nodes to generate" in ui.sink.texts(Kind.INFO)


def test_fill_numbering_rules(ui):
    m = run(ui, "N,10,0,0,0\nN,20,10,0,0\nFILL,10,20,4\n")      # d = 10/5 = 2 -> 12, 14, 16, 18
    assert [n for n in sorted(m.nodes)] == [10, 12, 14, 16, 18, 20]
    assert xyz(m, 14)[0] == pytest.approx(4.0)
    run(ui, "N,30,0,0,0\nN,33,3,0,0\nFILL,30,33,1\n")            # d = 1.5 not integer -> n1+k, warning
    assert 31 in m.nodes and any("not an integer" in w for w in warnings(ui))
    run(ui, "FILL,10,20,4")                                         # would collide with existing nodes
    assert any("already exist" in e for e in errors(ui))


def test_ngen_defaults_and_overwrite(ui):
    m = run(ui, "N,1,0,0,0\nN,2,1,0,0\nNGEN,2,,,,,0,1,0\n")     # pattern 1..2, step 2
    assert sorted(m.nodes) == [1, 2, 3, 4, 5, 6]
    assert xyz(m, 6) == pytest.approx([1, 2, 0])
    run(ui, "NGEN,1,1,1,2,1,0,0,5")                                 # node 2 overwritten (target of 1)
    assert any("overwritten" in w for w in warnings(ui))


def test_generation_in_active_system(ui):
    m = run(ui, "LOC,1,0,0,0,0,90,0,0\nN,1,1,0,0\nCSYS,1\nNGEN,1,1,1,1,1,1,0,0\n")
    # node 1 is (1,0,0) global = (0,-1,0) in system 1; +1 along local x
    assert m.nodes[2].csys == 1
    assert xyz(m, 2) == pytest.approx([1, -1, 0], abs=1e-12)
    assert m.node_global(2) == pytest.approx([1, 1, 0], abs=1e-12)


def test_nmed_lmove_nmove(ui):
    m = run(ui, "N,1,0,0,0\nN,2,2,0,0\nN,7,2,2,0\nN,6,0,2,0\nNMED,500,1,2,7,6\n")
    assert xyz(m, 500) == pytest.approx([1, 1, 0])
    run(ui, "N,3,1,1,1\nN,4,2,2,2\nLMOVE,0,0,10,101,1,2,3,4\n")
    assert [xyz(m, 101 + i)[2] for i in range(4)] == [10, 10, 11, 12]
    run(ui, "NMOVE,2,,,201,1,2,3\n")                               # blank factors = 1
    assert xyz(m, 202) == pytest.approx([4, 0, 0]) and xyz(m, 203) == pytest.approx([2, 1, 1])
    run(ui, "NMOVE,1,0,1,301,3\n")                                  # explicit 0 kept with a warning
    assert xyz(m, 301) == pytest.approx([1, 0, 1])
    assert any("collapses" in w for w in warnings(ui))
    run(ui, "LMOVE,0,0,1,400\nLMOVE,0,0,1,400,999\n")
    assert any("at least one node" in e for e in errors(ui))
    assert any("[999] are not defined" in e for e in errors(ui))


def test_ut13_nscale_zero_rule(ui):
    m = run(ui, "N,1,1,1,1\nN,2,2,2,2\nNSCALE,1,1,,0,2,0\n")
    assert xyz(m, 1) == pytest.approx([1, 2, 1])                   # y only
    assert xyz(m, 2) == pytest.approx([2, 2, 2])
    run(ui, "NSCALE,1,2,1,0.3048,0.3048,0.3048")
    assert xyz(m, 2) == pytest.approx([0.6096] * 3)
    # scaling does not change the definition order
    assert m.node_history == [1, 2]


# ------------------------------------------------------------------ fixities and flags
def test_ut13_d_labels(ui):
    m = run(ui, "N,1,0,0,0\nN,2,0,0,0\nN,5,0,0,0\nD,1,1,,1,DISP\n")
    assert m.nodes[1].fix == [1, 1, 1, 0, 0, 0]
    run(ui, "D,1,5,1,1,ALL\n")
    assert all(m.nodes[n].fix == [1] * 6 for n in (1, 2, 5))
    run(ui, "D,1,1,,,ALL\n")                                       # val defaults to 0: frees all
    assert m.nodes[1].fix == [0] * 6
    run(ui, "D,2,2,,0,UX,ROTZ\nD,5\nD,5,5,1,1,UQ\nD,5,5,1,2,UX\n")
    assert m.nodes[2].fix == [0, 1, 1, 1, 1, 0]
    e = errors(ui)
    assert any("at least one DOF label" in t for t in e)
    assert any("unknown DOF label UQ" in t for t in e)
    assert any("must be 0 (free) or 1 (fixed)" in t for t in e)


def test_int_flags_and_intlist(ui):
    m = run(ui, "NGEN,1,1,1,1\n")                                   # no nodes -> error
    assert any("no nodes defined" in e for e in errors(ui))
    run(ui, "N,1,0,0,0\nN,2,1,0,0\nN,3,2,0,0\nINT,1,3,,1\nINT,2,2,1,0,0\nINT,3,3,1,1,2\nINT,1,1,1,1,5\n")
    assert m.nodes[1].flags == {0} and m.nodes[2].flags == set() and m.nodes[3].flags == {0, 2}
    assert any("<code> must be" in e for e in errors(ui))
    ui.sink.clear()
    run(ui, "INTLIST,,,,1")
    assert "2 nodes listed" in ui.sink.texts(Kind.INFO)


def test_redefinition_keeps_fixities_and_flags(ui):
    m = run(ui, "N,1,0,0,0\nD,1,1,1,1,UX\nINT,1,1,1,1\nN,1,5,5,5\n")
    assert m.nodes[1].fix[0] == 1 and m.nodes[1].flags == {0}
    assert xyz(m, 1) == pytest.approx([5, 5, 5])


# ------------------------------------------------------------------ deletion cascade (D-MDL-11)
def test_ndel_cascade(ui):
    m = run(ui, "N,1,0,0,0\nN,2,1,0,0\nN,3,2,0,0\nGROUP,1,SPRING\nE,1,1,2\nF,2,1\nMM,2,0,1\nMT,2,1,1,1\n"
                "MR,2,1,1,1\nMUNITS,2,2,1,0\nNDEL,2\n")
    assert 2 not in m.nodes
    assert not m.forces and not m.moments and not m.tmass and not m.rmass and not m.mass_units
    assert m.groups[1].elements[1].nodes == [1, 2]                  # dangling reference kept (Error 41)
    assert any("Error 41" in w for w in warnings(ui))
    assert m.node_history == [1, 3]
    run(ui, "NDEL,1,10,2")                                          # 1, 3 (5, 7, 9 do not exist)
    assert not m.nodes


def test_nlist_columns(ui):
    run(ui, "N,1,1.5,2,3\nD,1,1,1,1,UZ\nINT,1,1,1,1,0\nINT,1,1,1,1,3\nNLIST")
    row = [t for t in ui.sink.texts(Kind.INFO) if t.strip().startswith("1 ")][0]
    assert "1.5" in row and "0  0  1  0  0  0" in row and row.rstrip().endswith("IN")


# ------------------------------------------------------------------ generated ids must be positive
def test_ngen_refuses_nonpositive_ids_and_changes_nothing(ui):
    m = run(ui, "N,1\nN,2,1")
    ui.execute("NGEN,1,-5,1,2")                      # would create nodes -4 and -3
    assert sorted(m.nodes) == [1, 2]
    assert any("must be positive" in e for e in errors(ui))
    ui.execute("NGEN,3,-1,1,2")                      # k = 1 gives 0..1, k = 3 gives -2..-1
    assert sorted(m.nodes) == [1, 2]
    ui.sink.clear()
    m = run(ui, "N,21,0,0,0\nN,22,1,0,0\nNGEN,2,-10,21,22,1,0,1,0")      # descending numbering is fine
    assert not errors(ui)
    assert {11, 12, 1, 2} <= set(m.nodes) and m.nodes[1].xyz == (0.0, 2.0, 0.0)


@pytest.mark.parametrize("cmd", ["LMOVE", "NMOVE"])
def test_lmove_nmove_refuse_nonpositive_nd(ui, cmd):
    m = run(ui, "N,1\nN,2,1")
    ui.execute(f"{cmd},0,0,1,-3,1,2")
    ui.execute(f"{cmd},0,0,1,0,1,2")
    assert sorted(m.nodes) == [1, 2]
    assert sum("must be positive" in e for e in errors(ui)) == 2
