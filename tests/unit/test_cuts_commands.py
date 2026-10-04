"""Cut and submodel commands: CUTADD CUTRMV CUTVOL SLICE CUTCLR CUT2SUB EXTRACTEXCAV SPLITGROUP TRANELEM TRANVOL
(requirements 3.4.J; spec 09 sections 5.1-5.12; D-MDL-13)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.prep import Interpreter, Kind
from sassi.prep import cuts_lib as cl
from sassi.prep.writer import write_pre
from sassi.verify.problems.vp_generation import box_lines


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def run(ui, lines):
    ui.run_text("\n".join(lines) if isinstance(lines, list) else lines)
    return ui.model


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


def infos(ui):
    return ui.sink.texts(Kind.INFO)


def roundtrip(m):
    text, _ = write_pre(m)
    b = Interpreter()
    b.run_text(text)
    assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)
    return b.model


def grid_lines(nx, ny, nz, h=1.0, z0=0.0, group=1, first_node=1, extra=()):
    """nx x ny x nz unit-ish hexahedra (size h), nodes numbered x fastest from ``first_node``."""
    nid = lambda i, j, k: first_node + i + (nx + 1) * j + (nx + 1) * (ny + 1) * k     # noqa: E731
    out = [f"N,{nid(i, j, k)},{h * i},{h * j},{z0 + h * k}"
           for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)]
    out.append(f"GROUP,{group},SOLID")
    e = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                e += 1
                ids = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                       nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                out.append(f"E,{e}," + ",".join(str(x) for x in ids))
    return out + list(extra)


def cut(ui, n):
    return sorted(ui.session["cuts"][n])


# ======================================================================================
# Cut sets (D-MDL-13: session-global sets of (group, element))
# ======================================================================================
def test_cutadd_list_and_range_forms(ui):
    run(ui, grid_lines(2, 2, 2, 0.5))
    assert ui.execute("CUTADD,3,1,1,2,2-4")
    assert cut(ui, 3) == [(1, 1), (1, 2), (1, 3), (1, 4)]
    assert isinstance(ui.session["cuts"][3], set)
    assert ui.execute("CUTADD,3,1,RANGE,5,9,2")            # 5, 7 (9 does not exist)
    assert cut(ui, 3) == [(1, k) for k in (1, 2, 3, 4, 5, 7)]
    assert any("not in group 1: 9" in w for w in warnings(ui))
    assert ui.execute("cutadd,3,1,range,8")                   # end defaults to start; case-insensitive
    assert (1, 8) in ui.session["cuts"][3]
    assert not ui.execute("CUTADD,3,7,1")                     # group 7 not defined
    assert len(ui.session["cuts"][3]) == 7
    assert not ui.execute("CUTADD,0,1,1")                     # cut numbers >= 1


def test_cutrmv_and_undefined_cut(ui):
    run(ui, grid_lines(2, 2, 2, 0.5) + ["CUTADD,1,1,RANGE,1,8"])
    assert ui.execute("CUTRMV,1,1,RANGE,1,3")
    assert cut(ui, 1) == [(1, k) for k in range(4, 9)]
    assert ui.execute("CUTRMV,1,1,4 6;99")
    assert cut(ui, 1) == [(1, 5), (1, 7), (1, 8)]
    assert not ui.execute("CUTRMV,9,1,1")
    assert "cut 9 is not defined" in errors(ui)[-1]


def test_cutvol_manual_example_blank_bounds(ui):
    """Manual 5.8.1: ``cutvol,3,,,,,2.53`` selects the elements above the ground surface at 2.53."""
    lines = []
    zs = [0.0, 1.265, 2.53, 3.53, 4.53]
    for k, z in enumerate(zs):
        for j in range(4):
            x, y = [(0, 0), (1, 0), (1, 1), (0, 1)][j]
            lines.append(f"N,{4 * k + j + 1},{x},{y},{z}")
    lines.append("GROUP,1,SOLID")
    for k in range(4):
        b = [4 * k + j + 1 for j in range(4)]
        lines.append(f"E,{k + 1}," + ",".join(str(x) for x in b + [n + 4 for n in b]))
    run(ui, lines)
    assert ui.execute("cutvol,3,,,,,2.53")
    assert cut(ui, 3) == [(1, 3), (1, 4)]
    # nodes 1e-7 below the bound are inside: the tolerance is 1e-6 x the model size (D-GEN-06), 4.7e-6 here
    assert ui.execute("CUTVOL,4,,,,,2.5300001")
    assert cut(ui, 4) == [(1, 3), (1, 4)]
    assert ui.execute("CUTVOL,5,,,,,2.6")
    assert cut(ui, 5) == [(1, 4)]
    assert not ui.execute("CUTVOL,6,2,1")                     # Xmin > Xmax


def test_slice_spec09_cube_of_eight(ui):
    """Spec 09 section 8 item 11: SLICE through x = 0.5 of a 1 x 1 x 1 cube of 2 x 2 x 2 solids selects all 8."""
    run(ui, grid_lines(2, 2, 2, 0.5))
    assert ui.execute("SLICE,1,0.5,0,0,1,0,0")
    assert cut(ui, 1) == [(1, k) for k in range(1, 9)]
    assert ui.execute("SLICE,2,0.25,0,0,1,0,0")               # crosses the first column only
    assert cut(ui, 2) == [(1, 1), (1, 3), (1, 5), (1, 7)]
    assert ui.execute("SLICE,2,0,0,0.75,0,0,1")               # adds (does not replace)
    assert cut(ui, 2) == [(1, 1), (1, 3), (1, 5), (1, 6), (1, 7), (1, 8)]
    assert not ui.execute("SLICE,3,0,0,0,0,0,0")              # zero normal


def test_cutclr_ranges_and_cleared_cut_is_an_error(ui):
    run(ui, grid_lines(2, 1, 1) + [f"CUTADD,{k},1,1" for k in (1, 2, 3, 5)])
    assert ui.execute("CUTCLR,1,3")
    assert sorted(ui.session["cuts"]) == [5]
    assert ui.execute("CUTCLR,2,6,3")                         # 2, 5
    assert ui.session["cuts"] == {}
    assert ui.execute("CUTCLR,7")
    assert any("no cut defined" in w for w in warnings(ui))
    assert not ui.execute("CSECT,4,1,0.5,0,0,1,0,0")
    assert "cut 1 is not defined" in errors(ui)[-1]


def test_cuts_are_read_by_cutplot_and_follow_gcom(ui):
    from sassi.plotting.state import cut_elements
    run(ui, grid_lines(2, 1, 1, group=4) + ["CUTADD,2,4,2"])
    assert cut_elements(ui, 2) == [(4, 2)]
    assert ui.execute("GCOM")                                  # group 4 -> 1: the cut follows
    assert cut(ui, 2) == [(1, 2)]


# ======================================================================================
# CUT2SUB
# ======================================================================================
MIXED = grid_lines(2, 1, 1) + [
    "N,20,0,0,3", "N,21,2,0,3", "N,22,0,1,5", "N,23,5,5,5", "N,24,6,5,5",
    "GROUP,2,SHELL", "E,1,1,2,5,4", "THICK,1,1,1,0.25", "MSET,1,1,1,2",
    "GROUP,3,BEAMS", "E,1,20,21,22", "RSET,1,1,1,2",
    "GROUP,4,SPRING", "E,1,23,24", "RSET,1,1,1,1",
    "M,1,3000,0.2,0.15,0.05,0.05,1", "M,2,4000,0.2,0.15,0.05,0.05,1", "M,3,5000,0.2,0.15,0.05,0.05,1",
    "R,1,1,1,1,1,1,1", "R,2,2,2,2,2,2,2", "SC,1,1,1,1,0,0,0,0.05",
    "MT,1,10,10,10", "MT,12,7,7,7", "MUNITS,12,12,1,0", "D,20,20,,1,ALL", "INT,4,4,1,1,0",
    "GRAVITY,9.81", "GROUNDELEV,-1",
]


def test_cut2sub_copies_elements_nodes_tables_masses(ui):
    m0 = run(ui, MIXED)
    assert not errors(ui), errors(ui)
    run(ui, ["CUTADD,7,1,2", "CUTADD,7,2,1", "CUTADD,7,3,1", "CUT2SUB,7,3"])
    assert not errors(ui), errors(ui)
    assert ui.active_model == 0 and ui.model is m0
    s = ui.models[3]
    assert sorted((g.id, e.id) for g, e in s.iter_elements()) == [(1, 2), (2, 1), (3, 1)]
    assert s.groups[1].elements[2].nodes == m0.groups[1].elements[2].nodes
    assert sorted(s.nodes) == sorted({2, 3, 5, 6, 8, 9, 11, 12, 1, 4, 20, 21, 22})
    assert s.nodes[20].fix == [1] * 6 and s.nodes[4].flags == {0}
    assert sorted(s.materials) == [1, 2] and sorted(s.sections) == [2] and not s.springs
    assert sorted(s.tmass) == [1, 12] and s.mass_units == {12: 0}
    assert s.gravity == 9.81 and s.ground_elevation == -1
    assert s.groups[2].elements[1].thick == 0.25
    b = roundtrip(s)
    assert b.same_state(s)


def test_cut2sub_solid_turns_shells_into_thick_solids(ui):
    run(ui, ["N,1,0,0,0", "N,2,2,0,0", "N,3,2,1,0", "N,4,0,1,0", "N,5,3,0,0",
             "GROUP,1,SHELL", "E,1,1,2,3,4", "E,2,2,5,3", "THICK,1,2,1,0.2", "M,1,1,0.2,1,0,0,1",
             "CUTADD,1,1,1-2", "Cut2Sub,1,2,1"])
    assert not errors(ui), errors(ui)
    s = ui.models[2]
    assert s.groups[1].type == 1 and len(s.groups[1].elements) == 2
    zs = sorted({round(n.z, 12) for n in s.nodes.values()})
    assert zs == [-0.1, 0.1]
    assert min(s.nodes) > 5                                    # new nodes after the source numbering
    ui.execute("ACTM,2")
    assert ui.execute("CALCC")
    assert ui.session["calcc"]["Volume"] == pytest.approx((2.0 + 0.5) * 0.2)
    assert any("thick-shell SOLIDs" in t for t in infos(ui))


def test_cut2sub_argument_errors(ui):
    run(ui, grid_lines(1, 1, 1) + ["CUTADD,1,1,1"])
    assert not ui.execute("CUT2SUB,1,0")                      # dest = active model
    assert not ui.execute("CUT2SUB,2,1")                      # undefined cut
    ui.session["cuts"][4] = set()
    assert not ui.execute("CUT2SUB,4,1")                      # empty cut


# ======================================================================================
# EXTRACTEXCAV
# ======================================================================================
def test_extractexcav_copies_soil_elements_flags_and_layers(ui):
    lines = box_lines(2, 2, 1, h=5.0, etype=2) + [
        "L,1,5,0.12,2000,1000,0.05,0.05", "L,2,10,0.13,3000,1500,0.05,0.05", "TOPL,1,2",
        "N,100,0,0,0", "N,101,5,0,0", "N,102,5,5,0", "N,103,0,5,0",
        "GROUP,2,SHELL", "E,1,100,101,102,103", "THICK,1,1,1,0.5", "M,1,3000,0.2,0.15,0.05,0.05,1",
        "INTGEN,1"]
    m = run(ui, lines)
    assert not errors(ui), errors(ui)
    assert ui.execute("EXTRACTEXCAV,2")
    s = ui.models[2]
    assert sorted(s.groups) == [1] and len(s.groups[1].elements) == 4
    assert all(e.etype == 2 for e in s.groups[1].elements.values())
    assert sorted(s.layers) == [1, 2] and s.topl == [1, 2] and not s.materials
    assert {n for n, nd in s.nodes.items() if 0 in nd.flags} == {n for n in s.nodes if 0 in m.nodes[n].flags}
    assert ui.active_model == 0


def test_extractexcav_needs_explicit_etype(ui):
    run(ui, box_lines(1, 1, 1, h=5.0, etype=0))               # implicit ETYPE 0 below grade
    assert not ui.execute("EXTRACTEXCAV,1")
    assert "ETYPEGEN" in errors(ui)[-1]
    run(ui, ["ETYPEGEN,0"])
    assert ui.execute("EXTRACTEXCAV,1")
    run(ui, ["ACTM,5"] + grid_lines(1, 1, 1, z0=1.0))       # above grade, no excavation
    assert not ui.execute("EXTRACTEXCAV,6")


# ======================================================================================
# SPLITGROUP
# ======================================================================================
def test_splitgroup_by_axis_plane_moves_positive_side(ui):
    m = run(ui, grid_lines(1, 1, 4) + ["GTIT,1,Column", "CUTADD,1,1,2-4", "EOUT,1,1,1,1,1,1,0,0,0,0,0,0,1,1-4"])
    assert ui.execute("SPLITGROUP,1,2,Z")
    assert sorted(m.groups[1].elements) == [1, 2]
    g2 = m.groups[2]
    assert g2.type == 1 and g2.title == "Column_B" and sorted(g2.elements) == [1, 2]
    assert g2.elements[1].nodes == [9, 10, 12, 11, 13, 14, 16, 15]          # old element 3
    assert cut(ui, 1) == [(1, 2), (2, 1), (2, 2)]
    assert [(r.group, r.elements) for r in m.eout] == [(1, [1, 2]), (2, [1, 2])]
    assert not ui.execute("SPLITGROUP,1,2,W")
    assert ui.execute("SPLITGROUP,1,10,X")                      # nothing on the + side
    assert any("nothing split" in w for w in warnings(ui))


def test_splitgroup_by_the_plane_of_a_shell(ui):
    m = run(ui, grid_lines(1, 1, 4) + ["N,50,0,0,1", "N,51,1,0,1", "N,52,1,1,1", "N,53,0,1,1",
                                        "GROUP,5,SHELL", "E,1,50,51,52,53", "THICK,1,1,1,0.1"])
    assert ui.execute("SPLITGROUP,1,5")                          # plane z = 1, normal +Z
    assert sorted(m.groups[1].elements) == [1] and sorted(m.groups[6].elements) == [1, 2, 3]
    assert not ui.execute("SPLITGROUP,1,6")                      # <split> must be a SHELL group


# ======================================================================================
# TRANELEM / TRANVOL
# ======================================================================================
def test_tranelem_overwrite_semantics(ui):
    run(ui, ["ACTM,1", "N,1,9,9,9", "N,30,0,0,0", "N,31,1,0,0", "N,32,1,1,0",
             "GROUP,1,SHELL", "E,1,30,31,32", "E,5,30,31,32", "M,1,1,0.3,1,0,0,1", "ACTM,0"])
    m0 = run(ui, grid_lines(2, 1, 1) + ["M,1,3000,0.2,0.15,0.05,0.05,1", "M,2,1,0.2,0.1,0,0,1", "MSET,2,2,1,2",
                                        "D,1,1,,1,UX"])
    assert ui.execute("TRANELEM,1,1,1,2,1")
    d = ui.models[1]
    assert ui.active_model == 0 and ui.model is m0
    assert d.groups[1].type == 1                                  # type changed to SOLID (warning)
    assert any("changed from SHELL to SOLID" in w for w in warnings(ui))
    assert d.groups[1].elements[1].nodes == m0.groups[1].elements[1].nodes
    assert 5 in d.groups[1].elements                             # elements outside the range are kept
    assert tuple(d.nodes[1].xyz) == (0.0, 0.0, 0.0) and d.nodes[1].fix[0] == 1
    assert d.materials[1].val1 == 1 and d.materials[2].val1 == 1  # 1 kept (warning), 2 added
    assert any("kept as defined in the destination" in w for w in warnings(ui))
    assert not ui.execute("TRANELEM,0,1,1,2,1")                   # dest = active
    assert not ui.execute("TRANELEM,1,9,1,1,1")                   # undefined group


def test_tranvol_box_into_new_model(ui):
    run(ui, grid_lines(2, 2, 1))
    assert ui.execute("TRANVOL,4,,1")
    d = ui.models[4]
    assert sorted(d.groups[1].elements) == [1, 3]
    assert sorted(d.nodes) == sorted({n for e in (1, 3) for n in ui.model.groups[1].elements[e].nodes})
    assert any("model 4 created" in t for t in infos(ui))
    assert not ui.execute("TRANVOL,5,10,11")


# ======================================================================================
# CSECT model structure (spec 09 section 5.2)
# ======================================================================================
def test_csect_spec09_unit_thickness_patch_keeps_element_numbers(ui):
    """Spec 09 section 8 item 11: CSECT of the 8-solid cube at x = 0.5 gives a 2 x 2 patch of unit-thickness
    solids with the original element numbers (the face on the plane counted once, -n side, D-SEC-04)."""
    run(ui, grid_lines(2, 2, 2, 0.5) + ["SLICE,1,0.5,0,0,1,0,0", "CSECT,4,1,0.5,0,0,1,0,0"])
    assert not errors(ui), errors(ui)
    s = ui.models[4]
    assert sorted(s.groups[1].elements) == [1, 3, 5, 7]
    assert sorted({round(n.x, 12) for n in s.nodes.values()}) == [0.0, 1.0]
    assert len(s.nodes) == 18                                     # 3 x 3 grid at x = 0 and x = 1
    assert any("counted once" in t for t in infos(ui))
    assert s.ui_state["csect"]["cut"] == 1 and s.ui_state["csect"]["normal"] == [1.0, 0.0, 0.0]
    b = roundtrip(s)
    assert b.same_state(s)


def test_csect_hexagon_is_split_into_two_solids(ui):
    run(ui, grid_lines(1, 1, 1, extra=["CUTADD,1,1,1", "CSECT,2,1,0.5,0.5,0.5,1,1,1"]))
    s = ui.models[2]
    assert sorted(s.groups[1].elements) == [1, 2]                  # sub-element numbered after the max
    assert s.ui_state["csect"]["parents"] == [[1, 2, 1]]
    assert len(s.ui_state["csect"]["pieces"][0]["pts"]) == 6


def test_csect_shells_beams_and_errors(ui):
    run(ui, ["N,1,0,0,0", "N,2,2,0,0", "N,3,2,0,2", "N,4,0,0,2", "N,5,3,1,0", "N,6,3,1,2", "N,7,4,1,1",
             "GROUP,1,SHELL", "E,1,1,2,3,4", "THICK,1,1,1,0.5",
             "GROUP,2,BEAMS", "E,1,5,6,7", "GROUP,3,SPRING", "E,1,1,6",
             "SLICE,1,0,0,1,0,0.6,0.8"])
    assert ui.execute("CSECT,3,1,0,0,1,0,0.6,0.8")
    s = ui.models[3]
    sh = s.groups[1].elements[1]
    assert sh.thick == pytest.approx(0.5 / 0.8)                   # t / |m x n|, m = Y
    bm = s.groups[2].elements[1]
    a, b = s.node_global(bm.nodes[0]), s.node_global(bm.nodes[1])
    assert np.allclose(b - a, [0, 0.6, 0.8])                      # unit length along n
    assert np.allclose(s.node_global(bm.nodes[2]), [4, 1, 1])     # K node position kept
    assert 3 not in s.groups                                       # springs are not section elements
    assert any("SPRING" in w for w in warnings(ui))
    assert not ui.execute("CSECT,0,1,0,0,1,0,0,1")                 # dest = active
    assert not ui.execute("CSECT,4,1,0,0,9,0,0,1")                 # the plane misses every element


def test_cut2sub_with_local_systems_round_trips(ui):
    """Nodes defined in LOC / LOCAL systems keep their systems; a LOCAL system whose defining nodes are not
    in the submodel is still written exactly (WRITE re-creates it from its stored points)."""
    run(ui, ["N,100,0,0,0", "N,101,1,0,0", "N,102,0,1,0", "LOCAL,2,0,100,101,102", "LOC,3,0,5,0,0,30,0,0",
             "CSYS,3"] + grid_lines(1, 1, 1) + ["CSYS,2", "N,50,0,0,4", "N,51,1,0,4", "N,52,1,1,4", "CSYS,0",
                                                "GROUP,2,SHELL", "E,1,50,51,52", "THICK,1,1,1,0.1",
                                                "CUTADD,1,1,1", "CUTADD,1,2,1", "CUT2SUB,1,1"])
    assert not errors(ui), errors(ui)
    s = ui.models[1]
    assert s.nodes[1].csys == 3 and s.nodes[50].csys == 2 and 100 not in s.nodes
    assert np.allclose(s.node_global(2), ui.model.node_global(2))
    assert roundtrip(s).same_state(s)


def test_manual_example_submodel_above_ground(ui, tmp_path):
    """Manual 5.8.1 / spec 04 section 8.2 example 1: the elements above the ground surface at 2.53 written as a
    submodel (the file name contains spaces; arguments are comma-delimited)."""
    (tmp_path / "out dir").mkdir()
    lines = grid_lines(1, 1, 3, h=2.53 / 2, z0=0.0) + ["M,1,3000,0.2,0.15,0.05,0.05,1", "GROUNDELEV,2.53",
                                                     "cutvol,3,,,,,2.53", "cut2sub,3,1", "Actm,1",
                                                     f"write,model above ground.pre,{tmp_path / 'out dir'}"]
    run(ui, lines)
    assert not errors(ui), errors(ui)
    assert sorted(ui.model.groups[1].elements) == [3]
    b = Interpreter(cwd=tmp_path)
    assert b.execute(f"INP,{tmp_path / 'out dir' / 'model above ground.pre'}")
    assert b.model.same_state(ui.model)
