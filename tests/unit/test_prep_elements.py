"""Group, element and attribute commands (spec 08 sections 4-5): UT-13 (EGEN grid), ECOMPR, GDEL,
attributes and their validity per group type."""
from __future__ import annotations

import pytest

from sassi.prep import Interpreter, Kind

GRID = "\n".join(f"N,{5 * r + i + 1},{10 * i},0,{10 * r}" for r in range(4) for i in range(5))   # 5 x 4 nodes


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


def test_ut13_egen_grid_example(ui):
    m = run(ui, GRID + "\nGROUP,1,PLANE\nE,1,1,2,7,6\nEGEN,3,1,1\nEGEN,2,5,1,4\n")
    g = m.groups[1]
    assert sorted(g.elements) == list(range(1, 13))
    assert g.elements[4].nodes == [4, 5, 10, 9]
    assert g.elements[5].nodes == [6, 7, 12, 11]
    assert g.elements[12].nodes == [14, 15, 20, 19]


def test_egen_increments_k_node_and_copies_attributes(ui):
    m = run(ui, "GROUP,2,BEAMS\nMACT,3\nRACT,4\nE,1,1,2,99\nKI,1,1,1,0,0,0,0,1,1\nETYPE,1,1,1,1\n"
                "MACT,1\nRACT,1\nEGEN,2,10,1\n")
    g = m.groups[2]
    assert g.elements[2].nodes == [11, 12, 109] and g.elements[3].nodes == [21, 22, 119]
    for e in (2, 3):
        el = g.elements[e]
        assert (el.mat, el.prop, el.ki, el.etype) == (3, 4, [0, 0, 0, 0, 1, 1], 1)
    run(ui, "EGEN,1,100,1,1,1,10")                   # explicit first new element number
    assert g.elements[10].nodes == [101, 102, 199]


def test_element_defaults_and_overwrite(ui):
    m = run(ui, "MACT,2\nRACT,5\nGROUP,1,SOLID\nE,1,1,2,3,4,5,6,7,8\nE,1,1,2,3,4,5,6,7,7\n")
    e = m.groups[1].elements[1]
    assert (e.mat, e.prop, e.etype, e.eint, e.thick, e.ki, e.kj) == (2, 5, 0, 0, 0.0, [0] * 6, [0] * 6)
    assert e.nodes == [1, 2, 3, 4, 5, 6, 7, 7]
    assert any("redefined" in w for w in warnings(ui))
    run(ui, "GROUP,3,SHELL\nE,1,1,2,3\nE,2,1,2,3,4,5\nE,3,0\n")
    assert m.groups[3].elements[1].nodes == [1, 2, 3]            # triangle
    assert any("at most 4 nodes" in t for t in errors(ui))
    assert any("node numbers required" in t for t in errors(ui))


def test_group_rules(ui):
    m = run(ui, "E,1,1,2\nGROUP,1\nGROUP,1,solid\nGROUP,2,5\nGROUP,7,Beam\nGROUP,8,gm\nGROUP,1,SHELL\n"
                "GROUP,9,WALL\nGROUP,10,6\nGROUP,1\n")
    e = errors(ui)
    assert "E: no active group (use GROUP)" in e
    assert any("group 1 does not exist" in t for t in e)
    assert any("use MTYPE" in t for t in e)
    assert any("unknown group type WALL" in t for t in e)
    assert any("is not 1, 2, 3, 4, 5, 7 or 9" in t for t in e)
    assert m.groups[1].type == 1 and m.groups[2].type == 5 and m.groups[7].type == 2 and m.groups[8].type == 9
    assert m.group_active == 1


def test_mtype_and_gtit(ui):
    m = run(ui, "GROUP,1,SOLID\nE,1,1,2,3,4,5,6,7,8\nMTYPE,SHELL\n")
    assert m.groups[1].type == 3
    assert any("Error 7" in w for w in warnings(ui))
    run(ui, "GROUP,2,BEAMS\nMTYPE,1,PLANE\n")
    assert m.groups[1].type == 4 and m.groups[2].type == 2
    run(ui, "GTIT,1,Reactor building walls, east\nGTIT,Beam frame, level 2\nGTIT,,Third title\n")
    assert m.groups[1].title == "Reactor building walls, east"
    assert m.groups[2].title == "Third title"
    run(ui, "GTIT,9,x")
    assert any("group 9 is not defined" in e for e in errors(ui))


def test_edel_ecompr_remaps_eout(ui):
    m = run(ui, GRID + "\nGROUP,1,PLANE\nE,1,1,2,7,6\nEGEN,3,1,1\nEGEN,2,5,1,4\n"
                       "EOUT,1,1,1,0,0,0,0,0,0,0,0,0,1,5-12\nEDEL,5\nECOMPR\n")
    g = m.groups[1]
    assert sorted(g.elements) == list(range(1, 12))
    assert g.elements[5].nodes == [7, 8, 13, 12]                  # old element 6, order kept
    assert m.eout[0].elements == list(range(5, 12))               # old 6..12 -> 5..11, deleted 5 dropped
    assert any("EOUT references" in w for w in warnings(ui))


def test_gdel_drops_eout(ui):
    m = run(ui, "GROUP,1,SPRING\nE,1,1,2\nGROUP,2,SPRING\nE,1,2,3\n"
                "EOUT,1,0,0,0,0,0,0,0,0,0,0,0,1,1\nEOUT,1,0,0,0,0,0,0,0,0,0,0,0,2,1\nGDEL,2\n")
    assert sorted(m.groups) == [1] and m.group_active is None
    assert [r.group for r in m.eout] == [1]
    assert any("EOUT requests" in w for w in warnings(ui))


def test_attributes_by_group_type(ui):
    m = run(ui, GRID + "\nGROUP,1,SOLID\nE,1,1,2,7,6,1,2,7,6\nE,2,2,3,8,7,2,3,8,7\nETYPE,1,2,1,2\n"
                       "EINT,2,2,1,1\nMSET,1,2,1,4\nTHICK,1,1,1,1.0\nKI,1,1,1,1\n")
    s = m.groups[1].elements
    assert (s[1].etype, s[2].etype, s[1].eint, s[2].eint, s[1].mat) == (2, 2, 0, 1, 4)
    e = errors(ui)
    assert any("THICK needs a SHELL or TSHELL group" in t for t in e)
    assert any("KI needs a BEAMS group" in t for t in e)
    run(ui, "GROUP,2,TSHELL\nE,1,1,2,7,6\nTHICK,1,1,1,0.5\nEINT,1,1,1,1\nEINT,1,1,1,2\nTHICK,1,1,1,0\n")
    t = m.groups[2].elements[1]
    assert (t.thick, t.eint) == (0.5, 1)
    assert any("TSHELL integration must be 0" in x for x in errors(ui))
    assert any("thickness must be > 0" in x for x in errors(ui))
    run(ui, "GROUP,3,SPRING\nE,1,1,2\nMSET,1,1,1,3\nRSET,1,1,1,6\nEINT,1,1,1,1\nETYPE,1,1,1,1\n")
    sp = m.groups[3].elements[1]
    assert (sp.mat, sp.prop, sp.eint) == (1, 6, 0)
    w = warnings(ui)
    assert any("MSET is not used by SPRING" in x for x in w)
    assert any("EINT applies to SOLID and TSHELL" in x for x in w)
    assert any("ETYPE has no effect on SPRING" in x for x in w)
    run(ui, "GROUP,4,PLANE\nE,1,1,2,7,6\nRSET,1,1,1,2\nETYPE,1,1,1,3\n")
    assert m.groups[4].elements[1].prop == 1
    assert any("RSET is not used by PLANE" in x for x in warnings(ui))
    assert any("<type> must be 0, 1 or 2" in x for x in errors(ui))


def test_beam_releases(ui):
    m = run(ui, "GROUP,1,BEAMS\nE,1,1,2,3\nE,2,2,4,3\nKI,1,2,1,0,0,0,0,1,1\nKJ,2,,,1,0,0,0,0,0\nKJ,1,1,1,2\n")
    g = m.groups[1].elements
    assert g[1].ki == [0, 0, 0, 0, 1, 1] and g[2].ki == [0, 0, 0, 0, 1, 1]
    assert g[2].kj == [1, 0, 0, 0, 0, 0] and g[1].kj == [0] * 6
    assert any("release codes are 0" in e for e in errors(ui))


def test_mact_ract_are_global(ui):
    m = run(ui, "MACT,3\nRACT,2\nGROUP,1,BEAMS\nE,1,1,2,3\nGROUP,2,SPRING\nE,1,1,2\n")
    assert m.groups[1].elements[1].mat == 3 and m.groups[2].elements[1].mat == 3
    assert m.groups[2].elements[1].prop == 2


def test_lists(ui):
    run(ui, "GROUP,1,BEAMS\nGTIT,1,Frame\nE,1,1,2,3\nGROUP,2,SOLID\nGLIST\nGROUP,1\nELIST")
    info = ui.sink.texts(Kind.INFO)
    assert "2 groups listed" in info and "1 elements listed" in info
    assert any("BEAMS" in t and "Frame" in t for t in info)


SPEC08_EXAMPLE = """GRAVITY,32.2
N,1,0,0,0
N,5,40,0,0
FILL
NGEN,4,5,1,5,1,0,10,0
INT,1,25,1,1,0
D,1,12,1,1,ROTZ
D,14,25,1,1,ROTZ
M,1,519120,0.17,0.150,0.04,0.04,1
GROUP,1,SHELL
GTIT,1,Basemat
MACT,1
E,1,1,2,7,6
EGEN,3,1,1
EGEN,3,5,1,4
THICK,1,16,1,2.0
N,100,20,20,30
NMED,101,1,5,21,25
NDEL,101
N,102,-10,20,0
D,102,102,1,1,ALL
GROUP,2,BEAMS
R,1,25,20.83,20.83,88.02,52.08,52.08
MACT,1
RACT,1
E,1,13,100,102
MUNITS,100,100,1,1
MT,100,500,500,500
L,1,20,0.120,2000,1000,0.02,0.02
"""


def test_spec08_worked_example(ui):
    """The illustrative surface-mat model of spec 08 section 10 loads without warnings or errors."""
    m = run(ui, SPEC08_EXAMPLE)
    assert not errors(ui) and not warnings(ui)
    assert sorted(m.nodes) == list(range(1, 26)) + [100, 102]
    assert m.nodes[25].xyz == (40.0, 40.0, 0.0) and m.nodes[13].xyz == (20.0, 20.0, 0.0)
    assert sum(1 for n in m.nodes.values() if 0 in n.flags) == 25
    assert m.nodes[13].fix == [0] * 6 and m.nodes[12].fix == [0, 0, 0, 0, 0, 1]
    shells = m.groups[1].elements
    assert len(shells) == 16 and shells[16].nodes == [19, 20, 25, 24] and shells[16].thick == 2.0
    assert m.groups[2].elements[1].nodes == [13, 100, 102]
    assert m.nodal_masses()[100][:3].tolist() == pytest.approx([500 / 32.2] * 3)


# ------------------------------------------------------------------ EGEN: valid ids only
def test_egen_refuses_nonpositive_element_and_node_numbers(ui):
    run(ui, "GROUP,1,BEAMS\nE,1,1,2,3")
    g = ui.model.groups[1]
    ui.execute("EGEN,1,1,1,1,1,-5")                  # ee < 1
    ui.execute("EGEN,1,1,1,1,1,0")
    ui.execute("EGEN,1,-2,1")                        # nodes 1, 2, 3 -> -1, 0, 1
    ui.execute("EGEN,2,-1,1")                        # second copy reaches node 0
    assert sorted(g.elements) == [1]
    assert sum("must be positive" in e or "start at 1" in e for e in errors(ui)) == 4
    ui.sink.clear()
    run(ui, "E,2,11,12,13\nEGEN,1,-10,2")            # 11, 12, 13 -> 1, 2, 3: valid
    assert not errors(ui) and g.elements[3].nodes == [1, 2, 3]


def test_egen_keeps_blank_node_slots(ui):
    run(ui, "GROUP,1,SOLID\nE,1,5,6,0,8,9,10,11,12\nEGEN,1,-4,1")
    assert not errors(ui)
    assert ui.model.groups[1].elements[2].nodes == [1, 2, 0, 4, 5, 6, 7, 8]
