"""MERGE and MERGESOIL (requirements 3.4.I; spec 09 sections 4.2, 4.5; spec 04 section 7; D-MDL-12)."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from sassi.model.materials import elastic_constants
from sassi.prep import Interpreter, Kind
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


CUBE = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def model_a():
    """Nodes 1-8 (cube) + 9-11 (beam), groups 1 (SOLID) and 3 (BEAMS), M 1-2, R 1, SC 4, LOC 2."""
    out = [f"N,{1 + k},{p[0]},{p[1]},{p[2]}" for k, p in enumerate(CUBE)]
    out += ["N,9,3,0,0", "N,10,4,0,0", "N,11,3,1,0", "M,1,4e6,0.2,0.15,0.05,0.05", "M,2,3e6,0.2,0.15,0.05,0.05",
            "R,1,1,1,1,1,1,1", "SC,4,1,1,1,0,0,0,0", "LOC,2,0,0,0,0,0,0,0", "L,1,10,0.12,1500,800,0.05,0.05",
            "GROUP,1,SOLID", "E,1,1,2,3,4,5,6,7,8", "GROUP,3,BEAMS", "MACT,2", "E,1,9,10,11", "FREQ,1,1,2,3"]
    return out


def model_b():
    """Cube 1-8 in LOC 1 (origin (10,0,0)), spring group 2 with SC 1, beam group 1 with R 2, GENERAL with MX 1."""
    out = ["LOC,1,0,10,0,0,90,0,0", "CSYS,1"]
    out += [f"N,{1 + k},{p[0]},{p[1]},{p[2]}" for k, p in enumerate(CUBE)]
    out += ["CSYS,0", "N,20,0,0,0", "N,21,0,0,-1", "M,1,1e6,0.3,0.1,0.02,0.02", "R,2,2,2,2,2,2,2",
            "SC,1,9,9,9,0,0,0,0", "MXR,1,1,1", "L,1,10,0.12,1500,800,0.05,0.05", "L,2,5,0.13,1600,900,0.05,0.05",
            "GROUP,1,BEAMS", "RACT,2", "E,1,1,2,20", "GROUP,2,SPRING", "RACT,1", "E,1,20,21",
            "GROUP,4,GENERAL", "E,1,20,21", "GROUP,5,SOLID", "MACT,1", "E,1,1,2,3,4,5,6,7,8",
            "D,21,21,1,1,ALL", "INT,20,20,1,1", "F,20,1,0,0", "MT,20,2,2,2", "NOUT,1,1,0,0,0,0,0,20",
            "EOUT,1,0,0,0,0,0,0,0,0,0,0,0,5,1", "RDND,20,1,0,0,0,0,0", "FREQ,1,9", "GROUNDELEV,-7"]
    return out


def test_merge_offsets_and_translation(ui):
    run(ui, ["ACTM,1"] + model_a() + ["ACTM,2"] + model_b() + ["ACTM,0", "MDL,merged,.", "MERGE,1,2,0,100,-5"])
    m = ui.model
    b = ui.models[2]
    assert m.name == "merged"
    assert sorted(m.nodes) == list(range(1, 12)) + list(range(12, 20)) + [31, 32]
    assert sorted(m.groups) == [1, 3, 4, 5, 7, 8]                  # Mdl2 groups 1, 2, 4, 5 -> +3
    assert m.groups[4].type == 2 and m.groups[4].elements[1].nodes == [12, 13, 31]
    assert m.groups[4].elements[1].prop == 3                         # R 2 -> 3 (offset 1)
    assert m.groups[5].elements[1].prop == 5 and m.springs[5].scx == 9   # SC 1 -> 5 (offset 4)
    assert m.groups[7].elements[1].prop == 1 and 1 in m.matrices         # MX 1 -> 1 (offset 0)
    assert m.groups[8].elements[1].mat == 3 and m.materials[3].val1 == 1e6   # M 1 -> 3 (offset 2)
    assert m.groups[4].elements[1].mat == 4 - 1                      # beam MACT 1 of model 2 -> 3
    assert sorted(m.layers) == [1, 2]                                # L shared, missing layer 2 added
    assert sorted(m.csys) == [2, 3] and m.nodes[12].csys == 3        # LOC 1 -> 3 (offset 2)
    for old in range(1, 9):
        assert np.allclose(m.node_global(old + 11), b.node_global(old) + [0, 100, -5], atol=1e-12)
    assert np.allclose(m.node_global(31), [0, 100, -5])
    assert m.nodes[32].fix == [1] * 6 and m.nodes[31].flags == {0}
    assert list(m.forces) == [31] and list(m.tmass) == [31]
    assert m.nout[-1].nodes == [31] and m.eout[-1].group == 8 and m.rdnd[-1].node == 31
    assert m.freq_sets == {1: [1, 2, 3]} and m.ground_elevation == 0.0     # options from Mdl1
    assert roundtrip(m).same_state(m)


def test_merge_same_model_twice_and_layer_conflict(ui):
    run(ui, ["ACTM,1"] + model_a() + ["ACTM,2", "L,1,10,0.12,1700,800,0.05,0.05", "N,1,0,0,0", "ACTM,3",
                                      "MERGE,1,1,50,0,0", "MERGE,1,2,0,0,0"])
    m = ui.models[3]
    assert m.layers[1].vp == 1500
    assert any("soil layer 1 differs" in w for w in ui.sink.texts(Kind.WARNING))
    ui.execute("MERGE,1,1,50,0,0")
    assert len(ui.model.nodes) == 22 and len(ui.model.groups) == 4
    assert not ui.execute("MERGE,1,9,0,0,0")


# ======================================================================================
# MERGESOIL
# ======================================================================================
SOIL = ["L,1,5,0.12,1500,800,0.05,0.05", "L,2,5,0.12,1500,800,0.05,0.05", "L,3,20,0.13,2000,1000,0.02,0.02",
        "TOPL,1,2,3"]


def _structure_and_soil(ui):
    run(ui, ["ACTM,1"] + basement_lines(2, 5.0, (-10.0, -5.0, 0.0)) + SOIL + ["EXCAV,2", "ACTM,3"])
    return ui.models[1], ui.models[2]


def test_mergesoil_mode1_merges_and_moves_flags(ui, tmp_path):
    s, q = _structure_and_soil(ui)
    q.nodes[1].flags.add(0)                                   # soil node 1 at the mat corner (node 1 of S)
    ui.execute("MERGESOIL,1,2,1,,,,soil.map")
    m = ui.model
    off = max(s.nodes)
    assert 0 in m.nodes[1].flags and 0 not in m.nodes[1 + off].flags
    soil_g = [g for g in m.groups.values() if g.type == 1][0]
    used = {n for e in soil_g.elements.values() for n in e.nodes}
    assert len(used & set(s.nodes)) == 25 and len(used) == 27
    mp = [ln.split() for ln in (tmp_path / "soil.map").read_text().splitlines() if not ln.startswith("#")]
    assert len(mp) == 27 and ["1", "1"] in mp and ["14", str(14 + off)] in mp   # 14: interior node at -5
    assert all(e.etype == 2 for e in soil_g.elements.values())
    assert roundtrip(m).same_state(m)


def test_mergesoil_modes_0_2_3_and_argument_checks(ui):
    _structure_and_soil(ui)
    ui.execute("MERGESOIL,1,2,0")
    assert len(ui.model.used_nodes()) == len(ui.model.nodes)
    assert not any(g.type == 7 for g in ui.model.groups.values())
    ui.execute("MERGESOIL,1,2,2,5e8")
    sg = [g for g in ui.model.groups.values() if g.type == 7][0]
    assert len(sg.elements) == 25 and {ui.model.springs[e.prop].scx for e in sg.elements.values()} == {5e8}
    assert sg.elements[1].nodes[1] <= max(ui.models[1].nodes) < sg.elements[1].nodes[0]   # (soil s, structure t)
    assert not ui.execute("MERGESOIL,1,2,3")                       # SepLevel required
    for sep, stiff in ((-5, 17), (0, 25)):                          # z = SepLevel counts as below
        ui.execute(f"MERGESOIL,1,2,3,,20,{sep}")
        m = ui.model
        sg = [g for g in m.groups.values() if g.type == 7][0]
        ks = [m.springs[e.prop].scx for e in sg.elements.values()]
        assert ks.count(1e7) == stiff and ks.count(20.0) == 25 - stiff
    assert not ui.execute("MERGESOIL,1,1")
    assert not ui.execute("MERGESOIL,1,2,4")


def test_mergesoil_converts_soil_materials_to_layers(ui):
    # soil model built with structural materials (e.g. converted from ANSYS): M 2 -> L 2, existing L 1 kept
    lines = ["ACTM,1"] + basement_lines(1, 5.0, (-5.0, 0.0)) + ["L,1,5,0.12,1500,800,0.05,0.05", "TOPL,1"]
    lines += ["ACTM,2", "M,1,9e5,0.3,0.11,0.04,0.03", "M,2,8e5,0.25,0.12,0.05,0.04"]
    lines += [f"N,{1 + k},{5 * p[0]},{5 * p[1]},{5 * p[2] - 5}" for k, p in enumerate(
        [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])]
    lines += ["GROUP,1,SOLID", "MACT,2", "E,1,1,2,3,4,5,6,7,8", "ETYPE,1,1,1,1", "GROUP,2,SOLID", "MACT,1",
              "E,1,1,2,3,4,5,6,7,8", "ETYPE,1,1,1,1", "ACTM,3", "MERGESOIL,1,2"]
    run(ui, lines)
    m = ui.model
    assert sorted(m.layers) == [1, 2] and m.layers[1].vp == 1500          # existing layer 1 kept
    assert any("soil material 1 differs" in w for w in ui.sink.texts(Kind.WARNING))
    k = elastic_constants(1, 8e5, 0.25, 0.12, 32.2)
    L2 = m.layers[2]
    assert math.isclose(L2.vs, k.Vs) and math.isclose(L2.vp, k.Vp) and (L2.pdamp, L2.sdamp) == (0.05, 0.04)
    assert L2.thick == 5.0 and L2.weight == 0.12
    assert sorted(m.materials) == [1]                                      # soil materials not copied
    soil = [g for gid, g in m.groups.items() if gid > 1]
    assert {e.mat for g in soil for e in g.elements.values()} == {1, 2}
    assert all(e.etype == 2 for g in soil for e in g.elements.values())


def test_mergesoil_joins_only_the_foundation_soil_interface(ui):
    """Basement with floor slabs at every level: every excavation node coincides with a structure node, but
    the interior node at -5 and the interior ground-surface node stay soil nodes (spec 01 rules 4, 11)."""
    run(ui, ["ACTM,1"] + basement_lines(2, 5.0, (-10.0, -5.0, 0.0), floors=True) + SOIL + ["EXCAV,2", "ACTM,3",
                                                                                         "MERGESOIL,1,2,1"])
    m = ui.model
    soil_g = [g for g in m.groups.values() if g.type == 1][0]
    used = {n for e in soil_g.elements.values() for n in e.nodes}
    assert len(used & set(ui.models[1].nodes)) == 25 and len(used) == 27
    assert any("2 soil nodes inside the excavation volume" in t for t in ui.sink.texts(Kind.INFO))
    ui.execute("EXCSTRCHK")
    assert not any("shared excavation interior nodes" in w for w in ui.sink.texts(Kind.WARNING))
    ui.execute("MERGESOIL,1,2,2")
    sg = [g for g in ui.model.groups.values() if g.type == 7][0]
    assert len(sg.elements) == 25


# ======================================================================================
# MERGESOIL: review defects (manual workflow, interface completeness, tolerances, parsing)
# ======================================================================================
def _stacked_soil(mats=(1, 2)):
    """Two 5 m soil hexes stacked from -10 to 0 (nodes 1-12), MSET 1 (lower) and 2 (upper)."""
    out = []
    for k, z in enumerate((-10.0, -5.0, 0.0)):
        out += [f"N,{1 + 4 * k + q},{5 * p[0]},{5 * p[1]},{z}" for q, p in enumerate(((0, 0), (1, 0), (1, 1), (0, 1)))]
    out += ["GROUP,1,SOLID", f"MACT,{mats[0]}", "E,1,1,2,3,4,5,6,7,8", f"MACT,{mats[1]}", "E,2,5,6,7,8,9,10,11,12"]
    return out


def test_mergesoil_manual_sequence_converts_the_m_table_of_an_etypegen2_soil_model(ui, tmp_path):
    """Spec 04 section 7.2: INP soil (M table, e.g. from ANSYS), ETYPEGEN,2, GROUNDELEV, MERGESOIL."""
    (tmp_path / "soil.pre").write_text("\n".join(_stacked_soil() + ["M,1,9e5,0.3,0.11,0.04,0.03",
                                                                    "M,2,8e5,0.25,0.12,0.05,0.04"]) + "\n")
    run(ui, ["ACTM,1"] + basement_lines(1, 5.0, (-10.0, -5.0, 0.0)) + ["ETYPEGEN,1", "ACTM,2", "INP,soil.pre",
                                                                        "ETYPEGEN,2", "GROUNDELEV,0", "ACTM,3",
                                                                        "MERGESOIL,1,2,1"])
    assert len(ui.models[2].nodes) == 12 and sorted(ui.models[2].materials) == [1, 2]
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    m = ui.model
    for lid, (E, nu, w, dp, ds) in ((1, (9e5, 0.3, 0.11, 0.04, 0.03)), (2, (8e5, 0.25, 0.12, 0.05, 0.04))):
        k = elastic_constants(1, E, nu, w, 32.2)
        L = m.layers[lid]
        assert math.isclose(L.vs, k.Vs) and math.isclose(L.vp, k.Vp) and (L.weight, L.pdamp, L.sdamp) == (w, dp, ds)
        assert L.thick == 5.0
    assert sorted(m.materials) == [1]                      # only the structure's material
    sg = [g for g in m.groups.values() if g.type == 1][0]
    assert [e.mat for e in sg.elements.values()] == [1, 2] and all(e.etype == 2 for e in sg.elements.values())
    assert len({n for e in sg.elements.values() for n in e.nodes} & set(ui.models[1].nodes)) == 12
    assert not ui.sink.texts(Kind.WARNING)


def test_mergesoil_soil_layer_wins_over_a_different_material_and_missing_mset_fails(ui):
    run(ui, ["ACTM,1"] + basement_lines(1, 5.0, (-10.0, -5.0, 0.0)) + ["ETYPEGEN,1", "ACTM,2"] + _stacked_soil()
        + ["L,1,5,0.12,1500,800,0.05,0.05", "M,1,9e5,0.3,0.11,0.04,0.03", "M,2,8e5,0.25,0.12,0.05,0.04",
           "ETYPEGEN,2", "ACTM,3", "MERGESOIL,1,2,1"])
    m = ui.model
    assert m.layers[1].vp == 1500 and m.layers[1].vs == 800            # L 1 of the soil model used
    assert m.layers[2].weight == 0.12                                   # M 2 converted
    assert any("MSET 1 is both soil layer 1 and material 1" in w for w in ui.sink.texts(Kind.WARNING))
    # MSET 3 is neither an L nor an M entry: error, active model unchanged
    before = json.dumps(ui.model.to_json(), sort_keys=True)
    ui.models[2].groups[1].elements[2].mat = 3
    assert not ui.execute("MERGESOIL,1,2,1")
    assert any("MSET 3 is neither a soil layer (L) nor a material (M)" in t for t in ui.sink.texts(Kind.ERROR))
    assert json.dumps(ui.model.to_json(), sort_keys=True) == before


def test_mergesoil_reports_interface_nodes_without_a_structure_partner(ui):
    _structure_and_soil(ui)
    s = ui.models[1]
    s.nodes[10].x = 0.3                                     # wall node at (0, 0, -5) moved: its soil node is alone
    ui.execute("MERGESOIL,1,2,1")
    off = max(s.nodes)
    w = [t for t in ui.sink.texts(Kind.WARNING) if "foundation-soil interface nodes" in t]
    assert len(w) == 1 and "1 of the 25" in w[0] and w[0].endswith(f"merged-model nodes {10 + off}")
    soil_g = [g for g in ui.model.groups.values() if g.type == 1][0]
    assert len({n for e in soil_g.elements.values() for n in e.nodes} & set(s.nodes)) == 24
    n0 = len(ui.sink.texts(Kind.WARNING))
    ui.execute("MERGESOIL,1,2,0")                           # unbonded: reported as information only
    assert len(ui.sink.texts(Kind.WARNING)) == n0
    assert any("1 of the 25" in t and "Mode 0" in t for t in ui.sink.texts(Kind.INFO))
    ui.execute("ACTM,1")
    ui.execute("TRANSLATE,0.01,0,0")
    ui.execute("ACTM,3")
    ui.execute("MERGESOIL,1,2,2")                           # no pair at all: the structure is not connected
    assert any("none of the 25" in t and "NOT connected" in t for t in ui.sink.texts(Kind.WARNING))
    assert not any(g.type == 7 for g in ui.model.groups.values())


def test_mergesoil_excav_delta_joins_scattered_levels_in_every_mode(ui):
    """EXCAV,2,0.01 on a basement whose node levels scatter by +-0.001: the excavation carries
    EDUOPT,GEOMTOL = 0.01, MERGESOIL joins all 25 interface nodes and passes the tolerance on, so INTGEN,3
    still recognises the ground-surface face of the merged model."""
    run(ui, ["ACTM,1"] + basement_lines(2, 5.0, (-10.0, -5.0, 0.0), jitter=0.001) + SOIL
        + ["EXCAV,2,0.01", "ACTM,2", "ETYPEGEN,2", "ACTM,3"])
    s = ui.models[1]
    for mode in ("1", "2", "3,,,-5"):
        ui.execute(f"MERGESOIL,1,2,{mode}")
        m = ui.model
        assert not any("interface nodes" in t for t in ui.sink.texts(Kind.WARNING))
        assert m.options.entry("EDUOPT", "GEOMTOL").arg(2) == "0.01"
        if mode == "1":
            sg = [g for g in m.groups.values() if g.type == 1][0]
            assert len({n for e in sg.elements.values() for n in e.nodes} & set(s.nodes)) == 25
            assert roundtrip(m).same_state(m)
        else:
            sp = [g for g in m.groups.values() if g.type == 7][0]
            assert len(sp.elements) == 25
            if mode.startswith("3"):                       # 9 mat + 8 at -5 (+-0.001, within tol) stiff
                assert sum(1 for e in sp.elements.values() if m.springs[e.prop].scx == 1e7) == 17
        ui.execute("INTGEN,0")
        ui.execute("INTGEN,3")
        assert sum(1 for n in m.nodes.values() if 0 in n.flags) == 25


def test_mergesoil_mapping_file_in_any_trailing_optional_slot(ui, tmp_path):
    _structure_and_soil(ui)
    for args, mode in (("a.map", 1), ("2,b.map", 2), (",c.map", 1)):
        assert ui.execute(f"MERGESOIL,1,2,{args}"), ui.sink.texts(Kind.ERROR)[-1:]
        name = args.split(",")[-1]
        assert (tmp_path / name).exists()
        assert any(g.type == 7 for g in ui.model.groups.values()) == (mode == 2)
