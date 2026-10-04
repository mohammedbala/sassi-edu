"""REFINEMODEL (manual section 9.16.2; spec 11 sections 1.2 and 7 item 1; requirements VP-54): quadrilaterals
1 -> 4, hexahedra 1 -> 8, shared midpoints, inherited attributes, flags and fixities, renumbering."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.prep import Interpreter, Kind
from sassi.prep import water_lib as wl
from sassi.prep.writer import write_pre
from sassi.verify.problems.vp_water import element_xyz, hex_volume, min_hex_jacobian, quad_area, shell_pool_lines


def make_ui(tmp_path, lines):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text("\n".join(lines))
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    ui.sink.clear()
    return ui


def grid_lines(n=2, typ="SHELL", z=0.0):
    lines = [f"N,{1 + i + (n + 1) * j},{i:.1f},{j:.1f},{z}" for j in range(n + 1) for i in range(n + 1)]
    lines += ["M,1,3e7,0.2,24,0.05,0.05", f"GROUP,1,{typ}"]
    lines += [f"E,{1 + i + n * j},{1 + i + (n + 1) * j},{2 + i + (n + 1) * j},{n + 3 + i + (n + 1) * j},"
              f"{n + 2 + i + (n + 1) * j}" for j in range(n) for i in range(n)]
    if typ in ("SHELL", "TSHELL"):
        lines.append(f"THICK,1,{n * n},1,0.2")
    return lines


def hex_lines(nx=1, ny=1, nz=1):
    ids = {}
    lines = []

    def nid(i, j, k):
        if (i, j, k) not in ids:
            ids[(i, j, k)] = len(ids) + 1
            lines.append(f"N,{ids[(i, j, k)]},{i},{j},{k}")
        return ids[(i, j, k)]
    els = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                els.append([nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                            nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)])
    lines += ["M,1,3e7,0.2,24,0.05,0.05", "GROUP,1,SOLID"]
    lines += ["E," + ",".join(map(str, [e + 1] + h)) for e, h in enumerate(els)]
    return lines


def test_quad_children_node_order_and_positions(tmp_path):
    ui = make_ui(tmp_path, ["N,1,0,0,0", "N,2,2,0,0", "N,3,2,2,0", "N,4,0,2,0", "M,1,3e7,0.2,24,0.05,0.05",
                            "GROUP,1,SHELL", "E,1,1,2,3,4", "THICK,1,1,1,0.3"])
    ui.execute("REFINEMODEL")
    m = ui.model
    # new nodes: m12, m23, m34, m41 (in creation order) and the centre
    pos = {n: tuple(m.node_global(n)) for n in m.nodes if n > 4}
    assert pos == {5: (1, 0, 0), 6: (2, 1, 0), 7: (1, 2, 0), 8: (0, 1, 0), 9: (1, 1, 0)}
    assert [list(e.nodes) for e in m.groups[1].sorted_elements()] == \
        [[1, 5, 9, 8], [5, 2, 6, 9], [9, 6, 3, 7], [8, 9, 7, 4]]
    conf = ui.sink.texts(Kind.CONFIRM)[-1]
    assert "1 quadrilaterals split into 4" in conf and "5 new nodes 5-9" in conf


def test_grid_refinement_counts_and_contiguous_numbers(tmp_path):
    ui = make_ui(tmp_path, grid_lines(2))
    ui.execute("REFINEMODEL")
    m = ui.model
    assert len(m.groups[1].elements) == 16 and len(m.nodes) == 25
    assert sorted(m.groups[1].elements) == list(range(1, 17))
    assert sum(quad_area(element_xyz(m, 1, e))[0] for e in m.groups[1].elements) == pytest.approx(4.0, rel=1e-14)
    # second refinement: 64 elements, 81 nodes
    ui.execute("REFINEMODEL")
    assert len(m.groups[1].elements) == 64 and len(m.nodes) == 81


@pytest.mark.parametrize("typ", ["TSHELL", "PLANE"])
def test_tshell_and_plane_quadrilaterals_are_refined(tmp_path, typ):
    if typ == "PLANE":           # 2D model in the X-Z plane
        lines = [f"N,{1 + i + 3 * j},{i:.1f},0,{j:.1f}" for j in range(3) for i in range(3)]
        lines += ["M,1,3e7,0.2,24,0.05,0.05", "GROUP,1,PLANE"]
        lines += [f"E,{1 + i + 2 * j},{1 + i + 3 * j},{2 + i + 3 * j},{5 + i + 3 * j},{4 + i + 3 * j}"
                  for j in range(2) for i in range(2)]
    else:
        lines = grid_lines(2, "TSHELL")
    ui = make_ui(tmp_path, lines)
    ui.execute("REFINEMODEL")
    m = ui.model
    assert len(m.groups[1].elements) == 16 and len(m.nodes) == 25
    if typ == "PLANE":
        assert all(abs(m.node_global(n)[1]) == 0.0 for n in m.nodes)


def test_hexahedra_share_face_and_edge_nodes(tmp_path):
    ui = make_ui(tmp_path, hex_lines(2, 1, 1) + ["ETYPE,1,2,1,1", "EINT,1,2,1,1"])
    ui.execute("REFINEMODEL")
    m = ui.model
    g = m.groups[1]
    assert len(g.elements) == 16 and len(m.nodes) == 5 * 3 * 3       # 45-node lattice, shared face
    assert sum(hex_volume(element_xyz(m, 1, e)) for e in g.elements) == pytest.approx(2.0, rel=1e-14)
    assert min(min_hex_jacobian(element_xyz(m, 1, e)) for e in g.elements) > 0
    assert all(e.etype == 1 and e.eint == 1 and e.mat == 1 for e in g.elements.values())


def test_shell_on_solid_face_stays_conforming(tmp_path):
    lines = hex_lines(1, 1, 1) + ["GROUP,2,SHELL", "MACT,1", "E,1,5,6,7,8", "THICK,1,1,1,0.1"]
    ui = make_ui(tmp_path, lines)
    ui.execute("REFINEMODEL")
    m = ui.model
    assert len(m.nodes) == 27                       # the shell reuses the hexahedron's top face nodes
    top = {n for e in m.groups[2].elements.values() for n in e.nodes}
    solid = {n for e in m.groups[1].elements.values() for n in e.nodes}
    assert top <= solid and len(top) == 9


def test_untouched_elements_and_hanging_nodes(tmp_path):
    lines = ["N,1,0,0,0", "N,2,1,0,0", "N,3,1,1,0", "N,4,0,1,0", "N,5,2,0,0", "M,1,3e7,0.2,24,0.05,0.05",
             "R,1,1,1,1,1,1,1", "SC,1,1,1,1,0,0,0,0", "GROUP,1,SHELL", "E,1,1,2,3,4", "E,2,2,5,3,3", "THICK,1,2,1,0.1",
             "N,6,0,0,1", "GROUP,2,BEAMS", "E,1,1,2,6", "GROUP,3,SPRING", "RACT,1", "E,1,1,5"]
    ui = make_ui(tmp_path, lines)
    ui.execute("REFINEMODEL")
    m = ui.model
    assert len(m.groups[1].elements) == 5
    assert list(m.groups[1].elements[5].nodes) == [2, 5, 3, 3]        # triangle unchanged, renumbered last
    assert list(m.groups[2].elements[1].nodes) == [1, 2, 6]
    assert list(m.groups[3].elements[1].nodes) == [1, 5]
    assert any("hanging nodes" in w for w in ui.sink.texts(Kind.WARNING))
    assert "3 elements unchanged" in ui.sink.texts(Kind.CONFIRM)[-1]


def test_flags_and_fixities_need_every_parent_corner(tmp_path):
    ui = make_ui(tmp_path, grid_lines(1) + ["INT,1,4,1,1", "D,1,2,1,1,UX,UZ", "D,1,1,1,1,UY", "INT,1,2,1,1,2"])
    m = ui.model
    ui.execute("REFINEMODEL")
    by = {tuple(np.round(m.node_global(n), 9)): m.nodes[n] for n in m.nodes}
    edge12 = by[(0.5, 0.0, 0.0)]           # between nodes 1 and 2: both interaction, both UX/UZ fixed
    assert 0 in edge12.flags and 2 in edge12.flags
    assert edge12.fix[:3] == [1, 0, 1]
    centre = by[(0.5, 0.5, 0.0)]           # all four corners are interaction nodes; only two are fixed
    assert centre.flags == {0} and centre.fix == [0] * 6
    assert any("5 new nodes are interaction nodes" in t for t in ui.sink.texts(Kind.INFO))


def test_eout_requests_and_cuts_follow_the_children(tmp_path):
    ui = make_ui(tmp_path, grid_lines(2) + ["EOUT,1,1,1,1,1,1,0,0,0,0,0,0,1,2,4", "CUTADD,3,1,3"])
    ui.execute("REFINEMODEL")
    m = ui.model
    assert m.eout[0].elements == [5, 6, 7, 8, 13, 14, 15, 16]
    assert ui.session["cuts"][3] == {(1, 9), (1, 10), (1, 11), (1, 12)}


def test_loads_masses_and_pool_warnings(tmp_path):
    ui = make_ui(tmp_path, grid_lines(1) + ["MT,3,10,10,10"])
    ui.execute("REFINEMODEL")
    assert any("not redistributed" in w for w in ui.sink.texts(Kind.WARNING))
    ui = make_ui(tmp_path, ["GRAVITY,9.81"] + shell_pool_lines(2.0, 2.0, 1.0, 2, 2, 1) + ["FILLPOOL,1e9"])
    ui.execute("REFINEMODEL")
    assert any("refine the pool walls before FILLPOOL" in w for w in ui.sink.texts(Kind.WARNING))


def test_nothing_to_refine_is_an_error(tmp_path):
    ui = make_ui(tmp_path, ["N,1,0,0,0", "N,2,1,0,0", "N,3,0,1,0", "M,1,3e7,0.2,24,0.05,0.05", "GROUP,1,SHELL",
                            "E,1,1,2,3", "THICK,1,1,1,0.1"])
    before = ui.model.canonical()
    ui.execute("REFINEMODEL")
    assert any("no quadrilateral" in e for e in ui.sink.texts(Kind.ERROR))
    assert ui.model.canonical() == before
    ui.execute("REFINEMODEL,1")
    assert any("ignored" in w for w in ui.sink.texts(Kind.WARNING))


def test_refined_model_write_inp_round_trip(tmp_path):
    ui = make_ui(tmp_path, hex_lines(1, 1, 2) + ["ETYPE,1,2,1,1", "INT,1,4,1,1"])
    ui.execute("REFINEMODEL")
    text, notes = write_pre(ui.model)
    assert not notes
    ui2 = Interpreter(cwd=tmp_path)
    ui2.run_text(text)
    assert not ui2.sink.texts(Kind.ERROR) and ui.model.same_state(ui2.model)


def test_refine_pool_then_fill(tmp_path):
    """The manual's use: refine the pool walls (the water mesh follows), then FILLPOOL."""
    ui = make_ui(tmp_path, ["GRAVITY,9.81"] + shell_pool_lines(2.0, 2.0, 2.0, 2, 2, 2))
    ui.execute("REFINEMODEL")
    ui.execute("FILLPOOL,1e9,0,1")
    assert not ui.sink.texts(Kind.ERROR)
    m = ui.model
    rec = wl.pool_record(m)
    water = m.groups[int(rec.water)]
    assert len(water.elements) == 4 * 4 * 3                     # 0.5 m mesh, 1.5 m of water
    assert sum(hex_volume(element_xyz(m, int(rec.water), e)) for e in water.elements) == pytest.approx(6.0, rel=1e-12)
