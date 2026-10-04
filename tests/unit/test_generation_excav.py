"""EXCAV and SOILMESH (requirements 3.4.I; spec 09 sections 4.1, 4.7; D-MDL-19)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.prep import Interpreter, Kind
from sassi.prep.generation_lib import signed_area
from sassi.prep.writer import write_pre
from sassi.verify.problems.vp_generation import basement_lines


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


SOIL = ["L,1,5,0.12,1500,800,0.05,0.05", "L,2,5,0.12,1600,850,0.05,0.05", "L,3,20,0.13,2000,1000,0.02,0.02",
        "TOPL,1,2,3"]


def well_oriented(m, gid):
    """Every SOLID: bottom face counter-clockwise seen from above, top face above the bottom face."""
    for e in m.groups[gid].elements.values():
        b = list(dict.fromkeys(e.nodes[:4]))
        t = list(dict.fromkeys(e.nodes[4:]))
        pb = np.array([m.node_global(n) for n in b])
        pt = np.array([m.node_global(n) for n in t])
        if signed_area(pb[:, :2]) <= 0 or not np.all(pt[:, 2] > pb[:, 2]) or not np.allclose(pt[:, :2], pb[:, :2]):
            return False
    return True


# ======================================================================================
# EXCAV
# ======================================================================================
def test_excav_basement_box(ui):
    run(ui, basement_lines(2, 5.0, (-10.0, -5.0, 0.0)) + SOIL + ["GRAVITY,9.81"])
    ui.execute("EXCAV,4")
    x = ui.models[4]
    assert ui.active_model == 0
    assert sorted(x.nodes) == list(range(1, 28)) and sorted(x.groups) == [1] and len(x.groups[1].elements) == 8
    for k, z in enumerate((-10.0, -5.0, 0.0)):                     # numbered bottom-up, level by level
        for j in range(3):
            for i in range(3):
                assert x.nodes[1 + i + 3 * j + 9 * k].xyz == (5.0 * i, 5.0 * j, z)
    els = x.groups[1].elements
    assert els[1].nodes == [1, 2, 5, 4, 10, 11, 14, 13]
    assert all(e.etype == 2 for e in els.values())
    assert [els[e].mat for e in (1, 4, 5, 8)] == [2, 2, 1, 1]      # TOPL layer at mid-height
    assert well_oriented(x, 1)
    assert x.ground_elevation == 0.0 and x.gravity == 9.81 and x.topl == [1, 2, 3] and sorted(x.layers) == [1, 2, 3]
    assert roundtrip(x).same_state(x)
    ui.execute("ACTM,4")
    ui.execute("INTGEN,1")
    assert sum(1 for n in x.nodes.values() if 0 in n.flags) == 27


def test_excav_adds_ground_level_and_warns_off_interface(ui):
    run(ui, basement_lines(2, 5.0, (-10.0, -4.0)) + SOIL + ["GROUNDELEV,0"])
    ui.execute("EXCAV,1")
    x = ui.models[1]
    assert sorted({n.z for n in x.nodes.values()}) == [-10.0, -4.0, 0.0]
    assert len(x.groups[1].elements) == 8
    assert any("levels -4 are not soil-layer interfaces" in w for w in ui.sink.texts(Kind.WARNING))


def test_excav_triangles_give_prisms_and_solid_basemat(ui):
    run(ui, ["N,1,0,0,-6", "N,2,4,0,-6", "N,3,0,4,-6", "N,4,0,0,-3", "N,5,4,0,-3", "N,6,0,4,-3", "GROUP,1,SHELL",
             "E,1,1,3,2", "E,2,4,5,6", "GROUNDELEV,0"])
    ui.execute("EXCAV,1")
    x = ui.models[1]
    assert len(x.nodes) == 9 and len(x.groups[1].elements) == 2
    assert x.groups[1].elements[1].nodes == [1, 2, 3, 3, 4, 5, 6, 6]   # CCW prism, nodes 3 = 4 and 7 = 8
    assert any("no TOPL" in w for w in ui.sink.texts(Kind.WARNING))
    # SOLID basemat: template from the bottom faces of the solids
    ui2 = Interpreter()
    lines = []
    for k, z in enumerate((-11.0, -10.0)):
        for j in range(3):
            for i in range(3):
                lines.append(f"N,{1 + i + 3 * j + 9 * k},{5 * i},{5 * j},{z}")
    lines.append("GROUP,1,SOLID")
    e = 0
    for j in range(2):
        for i in range(2):
            e += 1
            a = 1 + i + 3 * j
            lines.append(f"E,{e},{a},{a + 1},{a + 4},{a + 3},{a + 9},{a + 10},{a + 13},{a + 12}")
    lines += ["ETYPE,1,4,1,1", "EXCAV,2"]
    ui2.run_text("\n".join(lines))
    x = ui2.models[2]
    assert sorted({n.z for n in x.nodes.values()}) == [-11.0, -10.0, 0.0]
    assert len(x.groups[1].elements) == 8 and well_oriented(x, 1)


def test_excav_delta_clustering_and_errors(ui):
    run(ui, basement_lines(2, 5.0, (-10.0, -5.0, 0.0), floors=True, jitter=0.002) + SOIL)
    ui.execute("EXCAV,1,0.01")
    x = ui.models[1]
    assert len({n.z for n in x.nodes.values()}) == 3 and len(x.groups[1].elements) == 8
    ui.execute("EXCAV,1")
    assert len({round(n.z, 9) for n in ui.models[1].nodes.values()}) > 3
    ui.execute("EXCAV,1,-1")
    assert any("[delta] must be positive" in w for w in ui.sink.texts(Kind.WARNING))
    assert not ui.execute("EXCAV,0")                                 # the active model itself
    ui3 = Interpreter()
    ui3.run_text("N,1,0,0,-5\nN,2,5,0,-5\nN,3,0,5,-5\nGROUP,1,BEAMS\nE,1,1,2,3")
    assert not ui3.execute("EXCAV,1") and 1 not in ui3.models


def test_excav_delta_records_the_level_scatter_as_geomtol(ui):
    """delta > 0: the generated nodes sit on the level means, the structure nodes up to delta off them;
    EXCAV stores EDUOPT,GEOMTOL so that MERGESOIL joins them (spec 09 section 4.1 step 5; review defect)."""
    run(ui, basement_lines(2, 5.0, (-10.0, -5.0, 0.0), jitter=0.001) + SOIL)
    ui.execute("EXCAV,1,0.01")
    x = ui.models[1]
    assert x.options.entry("EDUOPT", "GEOMTOL").arg(2) == "0.01"
    assert len({n.z for n in x.nodes.values()}) == 3                  # geometry still snapped to the levels
    assert roundtrip(x).same_state(x)                                  # WRITE writes the tolerance
    ui.execute("EXCAV,1")                                              # delta 0: nothing recorded
    assert ui.models[1].options.entry("EDUOPT", "GEOMTOL") is None
    # a level snapped to gelev can be further than delta from its structure nodes: the larger scatter is used
    run(ui, ["ACTM,5"] + basement_lines(1, 5.0, (-5.0, -0.007), jitter=0.004) + SOIL + ["GROUNDELEV,0", "EXCAV,6,0.01"])
    x = ui.models[6]
    assert sorted({n.z for n in x.nodes.values()}) == [-5.0, 0.0]
    gt = float(x.options.entry("EDUOPT", "GEOMTOL").arg(2))
    assert 0.011 < gt < 0.011 * (1 + 1e-6)
    ui.execute("ACTM,6")
    ui.execute("ETYPEGEN,2")
    ui.execute("ACTM,7")
    ui.execute("MERGESOIL,5,6,1")
    sg = [g for g in ui.model.groups.values() if g.type == 1][0]
    assert len({n for e in sg.elements.values() for n in e.nodes} & set(ui.models[5].nodes)) == 8


# ======================================================================================
# SOILMESH
# ======================================================================================
def test_soilmesh_shared_nodes(ui):
    run(ui, basement_lines(2, 5.0, (-10.0, -5.0, 0.0)) + SOIL + ["M,1,4e6,0.2,0.15,0.05,0.05"])
    src = ui.model
    top = max(src.nodes)
    ui.execute("SOILMESH,5,50,50,2,1,0,0,5,0,1")
    s = ui.models[5]
    g = s.groups[max(src.groups) + 1]
    P, nL, F = 8, 3, 4
    assert len(g.elements) == (nL - 1) * P * 2 + 1 * (F + P * 2)     # 32 ring + 20 below = 52
    new = [n for n in s.nodes if n > top]
    assert len(new) == 2 * P * nL + (9 + 2 * P) and min(new) == top + 1
    assert set(s.nodes) - set(new) == {n for n in src.nodes if src.nodes[n].z <= 0}   # shared basement nodes
    assert well_oriented(s, g.id)
    assert all(e.etype == 1 for e in g.elements.values())
    # materials: type 3 copies of the TOPL layers, numbered after the source materials
    assert sorted(s.materials) == [2, 3, 4]
    assert {(mt.val1, mt.val2, mt.mtype) for mt in s.materials.values()} == {(1500, 800, 3), (1600, 850, 3),
                                                                              (2000, 1000, 3)}
    # ring 1 of the corner (0, 0): c = (5, 5), factor 1.5 -> (-2.5, -2.5); ring 2 -> (-5, -5)
    xy = {(n.x, n.y) for n in s.nodes.values()}
    assert (-2.5, -2.5) in xy and (-5.0, -5.0) in xy and (15.0, 15.0) in xy
    assert min(n.z for n in s.nodes.values()) == -15.0


def test_soilmesh_contact_springs_and_errors(ui):
    run(ui, basement_lines(2, 5.0, (-10.0, -5.0, 0.0)) + SOIL + ["SC,7,1e5,1e5,1e5,0,0,0,0.02"])
    ui.execute("SOILMESH,5,50,50,2,1,0,0,5,1,7")
    s = ui.models[5]
    sp = [g for g in s.groups.values() if g.type == 7]
    assert len(sp) == 1 and len(sp[0].elements) == 8 * 3 + 1        # wall nodes on 3 levels + mat interior
    assert s.springs[7].scx == 1e5 and all(e.prop == 7 for e in sp[0].elements.values())
    soil = [g for g in s.groups.values() if g.type == 1][0]
    base = set(ui.model.nodes)
    assert not ({n for e in soil.elements.values() for n in e.nodes} & base)   # no shared basement node
    assert roundtrip(s).same_state(s)
    ui2 = Interpreter()
    ui2.run_text("\n".join(basement_lines(2, 5.0, (-10.0, -5.0, 0.0))))
    assert not ui2.execute("SOILMESH,5,50,50,2,1,0,0,5,0,1")             # no TOPL profile
    ui.execute("SOILMESH,6,50,50,0,0,0,0,5,0,1")
    assert 6 not in ui.models
    ui.execute("SOILMESH,6,50,50,1,0,0,0,0,1,9")
    assert any("spring property 9" in w for w in ui.sink.texts(Kind.WARNING))
