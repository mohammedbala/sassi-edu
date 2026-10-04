"""GLB2LOC, GROUPMAT, RMVUNUSED, NCOM, GCOM, WELD, MERGEGROUP, ROTATE, TRANSLATE (requirements 3.4.I;
spec 09 sections 2.13-2.26, 4.3, 4.6, 4.8, 4.9; D-MDL-14, D-MDL-15)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.model.geometry import rot_x, rot_z
from sassi.prep import Interpreter, Kind
from sassi.prep import generation_lib as gl
from sassi.prep.writer import write_pre
from sassi.verify.problems.vp_generation import plate_lines


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def run(ui, lines):
    ui.run_text("\n".join(lines) if isinstance(lines, list) else lines)
    return ui.model


def roundtrip(m):
    text, _ = write_pre(m)
    b = Interpreter()
    b.run_text(text)
    assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)
    return b.model


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


CUBE = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def cube(first, x0=0.0, group=1):
    out = [f"N,{first + k},{x0 + p[0]},{p[1]},{p[2]}" for k, p in enumerate(CUBE)]
    return out + [f"GROUP,{group},SOLID", "E,1," + ",".join(str(first + k) for k in range(8))]


# ======================================================================================
# GLB2LOC
# ======================================================================================
def test_glb2loc_matches_the_analytic_transformation(ui):
    """Spec 09 section 8 item 7: node in system 1 (origin (1,2,3), 30 deg about Z) -> system 2."""
    m = run(ui, ["LOC,1,0,1,2,3,30,0,0", "LOC,2,0,-4,5,0.5,0,90,0", "CSYS,1", "N,7,2,1,-1", "CSYS,0", "N,8,3,3,3"])
    g7 = np.array([1, 2, 3]) + rot_z(30) @ np.array([2, 1, -1])
    R2 = rot_x(90)                                     # LOC tyz = 90: local axes = columns of Rx(90)
    ref7 = R2.T @ (g7 - np.array([-4, 5, 0.5]))
    ui.execute("GLB2LOC,7,8,1,2")
    assert m.nodes[7].csys == 2 and m.nodes[8].csys == 2
    assert np.allclose(m.nodes[7].xyz, ref7, atol=1e-12)
    assert np.allclose(m.node_global(7), g7, atol=1e-12) and np.allclose(m.node_global(8), [3, 3, 3], atol=1e-12)
    ui.execute("GLB2LOC,8,8,1,0")
    assert m.nodes[8].csys == 0 and np.allclose(m.nodes[8].xyz, [3, 3, 3], atol=1e-12)
    assert not ui.execute("GLB2LOC,7,7,1,9")


# ======================================================================================
# GROUPMAT
# ======================================================================================
def test_groupmat_one_material_per_group(ui):
    m = run(ui, plate_lines(1, 1, group=1) + ["M,2,3e6,0.3,0.1,0.02,0.02", "N,10,0,0,-1", "N,11,1,0,-1", "N,12,0,1,-1",
                                              "GROUP,2,SHELL", "MACT,2", "E,1,1,2,4", "E,2,10,11,12", "MACT,1",
                                              "E,3,1,10,12", "THICK,1,3,1,0.3",
                                              "GROUP,3,SPRING", "E,1,1,10", "SC,1,1,1,1,0,0,0,0"] + cube(20, 5.0, 4)
           + ["ETYPE,1,1,1,2", "MACT,1", "L,1,10,0.12,1500,800,0.05,0.05"])
    m.groups[4].elements[1].mat = 1
    ui.execute("GROUPMAT")
    assert sorted(m.materials) == [1, 2]
    assert m.materials[2].val1 == 3e6 and m.materials[1].val1 == 4.32e6
    assert {e.mat for e in m.groups[1].elements.values()} == {1}
    assert {e.mat for e in m.groups[2].elements.values()} == {2}
    assert m.groups[4].elements[1].mat == 1                       # excavated: L table index untouched
    assert any("group 2: 1 elements had another material" in w for w in warnings(ui))


def test_groupmat_undefined_first_material_fails_without_change(ui):
    m = run(ui, plate_lines(1, 1) + ["MACT,5", "GROUP,2,SHELL", "E,1,1,2,3"])
    before = m.canonical()
    assert not ui.execute("GROUPMAT")
    assert m.canonical() == before


# ======================================================================================
# RMVUNUSED, NCOM
# ======================================================================================
def test_rmvunused_keeps_used_k_and_interaction_nodes(ui):
    m = run(ui, ["N,1,0,0,0", "N,2,1,0,0", "N,3,0,1,0", "N,4,5,5,5", "N,5,6,6,6", "N,6,7,7,7", "N,7,8,8,8",
                 "GROUP,1,BEAMS", "E,1,1,2,3", "INT,5,5,1,1", "MT,4,1,1,1", "NOUT,1,1,0,0,0,0,0,4,6", "RDND,6,1,1,1,0,0,0",
                 "INT,7,7,1,1,2"])
    ui.execute("RMVUNUSED")
    assert sorted(m.nodes) == [1, 2, 3, 5]                        # K node 3 and interaction node 5 kept
    assert 4 not in m.tmass and m.nout == [] and m.rdnd == []
    w = warnings(ui)
    assert any("carried masses or loads" in t for t in w) and any("output requests" in t for t in w)


def test_ncom_remaps_every_reference(ui, tmp_path):
    m = run(ui, ["N,10,0,0,0", "N,20,1,0,0", "N,30,0,1,0", "N,40,0,0,1", "N,50,5,0,0", "N,60,5,1,0",
                 "LOCAL,1,0,20,30,40", "GROUP,1,BEAMS", "E,1,10,20,30", "GROUP,2,SPRING", "E,1,40,50",
                 "SC,1,1,1,1,1,1,1,0", "F,50,1,0,0", "MM,60,0,1,0", "MT,20,5,5,5", "MR,30,1,1,1", "MUNITS,30,30,1,0",
                 "D,40,40,1,1,ALL", "INT,50,50,1,1", "NOUT,1,1,0,0,0,0,0,20,60", "RDND,30,1,1,1,0,0,0",
                 "ME,1,15,45,0,0,0", "SYMM,1,0,10,20,30"])
    ui.execute("NCOM,map.txt")
    assert sorted(m.nodes) == [1, 2, 3, 4, 5, 6]
    assert m.groups[1].elements[1].nodes == [1, 2, 3] and m.groups[2].elements[1].nodes == [4, 5]
    assert m.csys[1].nodes == (2, 3, 4)
    assert list(m.forces) == [5] and m.forces[5].node == 5 and list(m.moments) == [6]
    assert list(m.tmass) == [2] and list(m.rmass) == [3] and m.mass_units == {3: 0}
    assert m.nodes[4].fix == [1] * 6 and m.nodes[5].flags == {0}
    assert m.nout[0].nodes == [2, 6] and m.rdnd[0].node == 3
    me = m.options.entry("ME", 1)
    assert (me.integer(2), me.integer(3)) == (2, 4)               # range 15..45 -> nodes 20..40 -> 2..4
    symm = m.options.entry("SYMM", 1)
    assert [symm.integer(k) for k in (3, 4, 5)] == [1, 2, 3]
    assert list(m.node_history) == [1, 2, 3, 4, 5, 6]
    lines = [ln for ln in (tmp_path / "map.txt").read_text().splitlines() if not ln.startswith("#")]
    assert lines[0] == "10 1" and lines[-1] == "60 6"
    assert roundtrip(m).same_state(m)


def test_ncom_refuses_dangling_references(ui):
    m = run(ui, ["N,5,0,0,0", "N,9,1,0,0", "GROUP,1,SPRING", "E,1,5,7"])
    before = m.canonical()
    assert not ui.execute("NCOM")
    assert m.canonical() == before


# ======================================================================================
# GCOM, MERGEGROUP
# ======================================================================================
def test_gcom_compresses_groups_and_requests(ui):
    m = run(ui, cube(1, 0.0, 2) + cube(11, 2.0, 5) + cube(21, 4.0, 9) + ["GTIT,5,middle", "EOUT,1,0,0,0,0,0,0,0,0,0,0,0,9,1",
                                                                       "GROUP,5"])
    ui.session["cuts"] = {1: {(5, 1), (9, 1)}, 2: {(2, 1)}}
    ui.execute("GCOM")
    assert sorted(m.groups) == [1, 2, 3] and m.groups[2].title == "middle"
    assert m.group_active == 2 and m.eout[0].group == 3
    assert ui.session["cuts"] == {1: {(2, 1), (3, 1)}, 2: {(1, 1)}}
    assert m.groups[3].elements[1].nodes == list(range(21, 29))


def test_mergegroup_same_type_only(ui):
    m = run(ui, cube(1, 0.0, 1) + cube(11, 2.0, 2) + cube(21, 4.0, 3) + ["GROUP,4,SPRING", "E,1,1,11",
                                                                       "EOUT,1,0,0,0,0,0,0,0,0,0,0,0,3,1"])
    before = m.canonical()
    assert not ui.execute("MERGEGROUP,1,2,4")
    assert m.canonical() == before
    ui.execute("MERGEGROUP,1,3,2")
    assert sorted(m.groups) == [1, 4]
    assert m.groups[1].elements[2].nodes == list(range(21, 29))   # argument order: group 3 first
    assert m.groups[1].elements[3].nodes == list(range(11, 19))
    assert m.eout[0].group == 1 and m.eout[0].elements == [2]
    ui.execute("GCOM")
    assert sorted(m.groups) == [1, 2] and m.groups[2].type == 7


# ======================================================================================
# WELD
# ======================================================================================
def test_weld_two_cubes(ui):
    m = run(ui, cube(1, 0.0, 1) + cube(9, 1.0, 2) + ["MT,13,1,1,1"])
    ui.execute("WELD")
    assert m.groups[2].elements[1].nodes == [2, 10, 11, 3, 6, 14, 15, 7]
    assert len(m.used_nodes()) == 12 and len(m.nodes) == 16
    assert 13 in m.tmass                                           # masses are not transferred
    assert any("not transferred" in w for w in warnings(ui))
    ui.execute("RMVUNUSED")
    assert len(m.nodes) == 12


def test_weld_protects_zero_length_springs_unless_force(ui):
    m = run(ui, plate_lines(1, 1) + ["FIXSHLROT", "N,50,0,0,0", "GROUP,9,SPRING", "E,1,50,4"])
    sp = [g for g in m.groups.values() if g.type == 7 and g.id != 9][0]
    pairs = [list(e.nodes) for e in sp.elements.values()]
    ui.execute("WELD")
    assert [list(e.nodes) for e in sp.elements.values()] == pairs      # drilling springs intact
    # node 50 coincides with the shell node 1 (a drilling-spring end) and its ground node: it is welded
    # to node 1, only the spring's own two ends are kept apart (D-MDL-15; review defect)
    assert m.groups[9].elements[1].nodes == [1, 4]
    assert any("zero-length SPRING" in t for t in ui.sink.texts(Kind.INFO))
    ui.execute("WELD,FORCE")
    assert all(len(set(e.nodes)) == 1 for e in sp.elements.values())
    assert any("connect a node to itself" in w for w in warnings(ui))
    assert not ui.execute("WELD,MAYBE")


def test_weld_cannot_link_keeps_only_the_spring_ends_apart(ui):
    """A node coincident with both ends of a zero-length spring (here lower-numbered than both) is welded
    to one of them; the spring is neither collapsed nor detached (D-MDL-15, union-find with cannot-link)."""
    m = run(ui, ["N,1,0,0,0", "N,2,1,0,0", "N,3,1,1,0", "N,4,0,1,0",          # plate B, corner 1 at the origin
                 "N,5,0,0,0", "N,6,-1,0,0", "N,7,-1,-1,0", "N,8,0,-1,0",      # plate A, corner 5 at the origin
                 "N,9,0,0,0", "M,1,4e6,0.2,0.15,0.05,0.05", "GROUP,1,SHELL", "E,1,1,2,3,4", "E,2,5,6,7,8",
                 "THICK,1,2,1,0.5", "SC,1,0,0,0,10,10,10,0", "GROUP,2,SPRING", "RACT,1", "E,1,5,9",
                 "D,9,9,1,1,ALL"])
    sub = gl.weld(m)
    assert sub == {5: 1}                                   # 9 would join the cluster of its spring partner 5
    assert m.groups[1].elements[2].nodes == [1, 6, 7, 8]   # plate A now shares the corner with plate B
    assert m.groups[2].elements[1].nodes == [1, 9]         # spring end 5 -> 1: still two distinct nodes
    sub = gl.weld(m, force=True)                           # FORCE: the spring collapses
    assert sub == {5: 1, 9: 1} and m.groups[2].elements[1].nodes == [1, 1]   # 5 is unused but still listed


def test_coincident_groups_unconstrained_clusters_are_unchanged():
    xyz = np.array([[0, 0, 0], [0, 0, 0], [5, 5, 5], [5, 5, 5], [5, 5, 5], [9, 9, 9.0]])
    ids = [10, 3, 7, 1, 4, 2]
    assert gl.coincident_groups(ids, xyz, 1e-9) == {10: 3, 4: 1, 7: 1}
    # cannot-link 1-7: 4 joins 1 (pair (1,4) first), 7 stays alone; the other cluster is not affected
    assert gl.coincident_groups(ids, xyz, 1e-9, [(1, 7)]) == {10: 3, 4: 1}


# ======================================================================================
# ROTATE, TRANSLATE
# ======================================================================================
@pytest.mark.parametrize("fix, rotation, warned", [
    ("ROTZ", "ROTATE,0,0,0,30,0,0", False),         # drilling restraint about Z: a rotation about Z keeps it
    ("ROTZ", "ROTATE,0,0,0,0,30,0", True),          # ... about X it becomes a bending restraint
    ("UX", "ROTATE,0,0,0,90,0,0", True),            # UX restrains global X; the rotated model needs Y
    ("UX", "ROTATE,0,0,0,180,0,0", False),          # x -> -x: same restrained line
    ("UZ", "ROTATE,0,0,0,45,0,0", False),
    ("ALL", "ROTATE,0,0,0,10,20,30", False),        # all-or-none is invariant
])
def test_rotate_warns_about_global_axis_fixities(ui, fix, rotation, warned):
    run(ui, ["N,1,0,0,0", "N,2,1,0,0", f"D,1,1,1,1,{fix}"])
    ui.execute(rotation)
    assert any("D fixities" in w for w in warnings(ui)) == warned, warnings(ui)


def test_rotate_fixrot_then_rotate_and_anisotropic_masses(ui):
    run(ui, plate_lines(2, 2) + ["FIXROT", "MT,1,5,5,5", "MT,2,5,5,9", "MR,3,1,2,3"])
    ui.execute("ROTATE,0,0,0,0,30,0")
    w = warnings(ui)
    fx = [t for t in w if "D fixities" in t]
    assert len(fx) == 1 and "9 nodes" in fx[0] and "FIXROT" in fx[0]
    assert any("2 nodes have direction-dependent MT/MR masses" in t and t.endswith(": 2, 3") for t in w)
    fixed, masses = gl.global_axis_data_changed(ui.model, np.eye(3))
    assert fixed == [] and masses == []                    # identity: nothing changes meaning
def test_rotate_order_and_local_systems(ui):
    m = run(ui, ["N,1,1,0,0", "LOC,1,0,10,0,0,0,0,0", "CSYS,1", "N,2,1,0,0", "CSYS,0", "N,3,0,0,0", "N,4,0,1,0",
                 "N,5,0,0,1", "LOCAL,2,0,3,4,5"])
    R2 = np.array(m.csys[2].R)
    ui.execute("ROTATE,0,0,0,90,90,0")
    # Rz(90) first: (1,0,0) -> (0,1,0); then Rx(90): (0,1,0) -> (0,0,1)
    assert np.array_equal(m.node_global(1), [0.0, 0.0, 1.0])
    T = rot_x(90) @ rot_z(90)
    assert np.allclose(m.node_global(2), T @ np.array([11.0, 0, 0]), atol=1e-12)   # system rotates too
    assert m.nodes[2].csys == 1 and m.nodes[2].xyz == (1.0, 0.0, 0.0)
    assert np.allclose(m.csys[2].R, T @ R2, atol=1e-12)
    assert np.allclose(m.csys[1].R, T, atol=1e-12) and np.allclose(m.csys[1].origin, T @ [10.0, 0, 0], atol=1e-12)
    assert roundtrip(m).same_state(m)


def test_rotate_about_a_point_and_warnings(ui):
    m = run(ui, ["N,1,2,1,0", "N,2,3,1,0", "GROUP,1,SPRING", "E,1,1,2", "SC,1,10,20,30,0,0,0,0", "F,1,1,0,0"])
    ui.execute("ROTATE,1,1,0,180,0,0")
    assert np.allclose(m.node_global(1), [0, 1, 0], atol=1e-15)
    w = warnings(ui)
    assert any("anisotropic" in t for t in w) and any("F/MM" in t for t in w)


def test_translate_moves_global_and_local_nodes(ui):
    m = run(ui, ["N,1,1,2,3", "LOC,1,0,5,5,5,90,0,0", "CSYS,1", "N,2,1,0,0"])
    g2 = m.node_global(2).copy()
    ui.execute("TRANSLATE,10,0,-1")
    assert m.nodes[1].xyz == (11.0, 2.0, 2.0)
    assert np.allclose(m.node_global(2), g2 + [10, 0, -1], atol=1e-12)
    assert np.allclose(m.csys[1].origin, [5, 5, 5])
    assert any("ground elevation is not shifted" in t for t in warnings(ui))
