"""Water modelling commands FILLPOOL, LISTPOOLINTER, MERGEPOOL, POOLDATA and the analytical references of
:mod:`sassi.prep.water_lib` (manual section 9.16; requirements 3.4.N; spec 11 sections 1 and 7; D-WAT-01)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.prep import Interpreter, Kind
from sassi.prep import registry
from sassi.prep import water_lib as wl
from sassi.prep.writer import write_pre
from sassi.verify.problems.vp_water import element_xyz, hex_volume, shell_pool_lines, solid_pool_lines


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def make_ui(tmp_path, lines=(), gravity=9.81):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text("\n".join([f"GRAVITY,{gravity}"] + list(lines)))
    ui.sink.clear()
    return ui


def pool_ui(tmp_path, a=4.0, b=2.0, H=3.0, nx=4, ny=2, nz=3, **kw):
    return make_ui(tmp_path, shell_pool_lines(a, b, H, nx, ny, nz, **kw))


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


def water_volume(m, rec):
    return sum(hex_volume(element_xyz(m, int(rec.water), e)) for e in m.groups[int(rec.water)].elements)


def spring_at(m, rec, xyz):
    for e in m.groups[int(rec.springs)].elements.values():
        if np.linalg.norm(m.node_global(e.nodes[1]) - np.asarray(xyz, float)) < 1e-9:
            return m.springs[e.prop].k
    return None


# --------------------------------------------------------------------------------------
# registry
# --------------------------------------------------------------------------------------
def test_commands_are_registered_handlers():
    registry.load_commands()
    for name in ("FILLPOOL", "REFINEMODEL", "LISTPOOLINTER", "MERGEPOOL", "POOLDATA"):
        spec = registry.REGISTRY[name]
        assert not spec.placeholder and spec.module.endswith("commands.water")
        assert spec.tier == "P2"
    assert registry.REGISTRY["POOLDATA"].cls == "record"


# --------------------------------------------------------------------------------------
# FILLPOOL
# --------------------------------------------------------------------------------------
def test_fillpool_defaults(tmp_path):
    ui = pool_ui(tmp_path)
    m = ui.model
    nmax = max(m.nodes)
    ui.execute("FILLPOOL")
    assert not errors(ui)
    rec = wl.pool_record(m)
    assert (rec.stiff, rec.sensitivity, rec.empty, rec.shellarea, rec.offset, rec.stiff2) == (1.0e6, 0.0, 0, -1, nmax, 0.0)
    assert int(rec.shells) == 0 and len(m.groups) == 3                     # no area shells by default
    assert min(m.groups[int(rec.water)].elements[1].nodes) == nmax + 1      # offset <= 0 -> largest node
    assert spring_at(m, rec, (0.0, 1.0, 1.0))[:3] == (1.0e6, 0.0, 0.0)     # frictionless x-wall
    assert water_volume(m, rec) == pytest.approx(4.0 * 2.0 * 3.0, rel=1e-12)
    assert float(rec.zfloor) == 0.0 and float(rec.zsurface) == 3.0
    assert m.groups[int(rec.water)].title == wl.TITLES["water"]
    assert all(e.etype == 1 for e in m.groups[int(rec.water)].elements.values())


def test_fillpool_messages(tmp_path):
    ui = pool_ui(tmp_path)
    ui.execute("FILLPOOL,1e9,0,1,1,0,0")
    infos = " ".join(ui.sink.texts(Kind.INFO))
    assert "Vp = 1483.24" in infos and "1e-8 K" in infos
    assert "not modelled: sloshing" in infos
    assert "water mass on the interface springs alone" in infos
    assert any("MOPT <incomp> set to 0" in w for w in warnings(ui))
    conf = ui.sink.texts(Kind.CONFIRM)[-1]
    assert "16 water SOLIDs" in conf and "39 interface springs" in conf and "32 interface-area shells" in conf


def test_fillpool_water_material_si_and_british(tmp_path):
    ui = pool_ui(tmp_path)
    ui.execute("FILLPOOL,1e9")
    m = ui.model
    mat = m.materials[int(wl.pool_record(m).watermat)]
    c = mat.constants(m.gravity)
    assert mat.mtype == 3 and mat.weight == 9.81 and mat.pdamp == mat.sdamp == 0.005
    assert c.M - 4 * c.G / 3 == pytest.approx(2.2e6, rel=1e-8)            # kPa
    assert c.G / (c.M - 4 * c.G / 3) == pytest.approx(1e-8, rel=1e-6)
    assert c.rho == pytest.approx(1.0)
    # British (kip, ft): 62.4 pcf = 0.0624 kcf, K = 2.2 GPa in ksf
    wp = wl.water_properties(32.2)
    assert wp.units == "BS" and wp.weight == 0.0624 and wp.recognised
    assert wp.K == pytest.approx(2.2e9 / 47880.259, rel=1e-7)
    assert wp.vp == pytest.approx(1483.24 / 0.3048, rel=2e-3)                # ft/s
    wsi = wl.water_properties(9.81)
    assert wsi.vs / wsi.vp == pytest.approx(1e-4, rel=1e-6)
    assert not wl.water_properties(1.0).recognised
    with pytest.raises(wl.GenerationError):
        wl.water_properties(0.0)


def test_fillpool_unrecognised_gravity_warns(tmp_path):
    ui = make_ui(tmp_path, shell_pool_lines(2.0, 2.0, 1.0, 2, 2, 1), gravity=1.0)
    ui.execute("FILLPOOL")
    assert any("neither ~9.81" in w for w in warnings(ui))


def test_fillpool_spring_constants_normals(tmp_path):
    ui = pool_ui(tmp_path, a=4.0, b=3.0, H=3.0, nx=4, ny=3, nz=3)
    ui.execute("FILLPOOL,5e8,0,0,-1,0,7e2")
    m = ui.model
    rec = wl.pool_record(m)
    S, s = 5e8, 7e2
    assert spring_at(m, rec, (0, 1, 1))[:3] == (S, s, s)          # x-wall
    assert spring_at(m, rec, (2, 3, 2))[:3] == (s, S, s)          # y-wall
    assert spring_at(m, rec, (2, 1, 0))[:3] == (s, s, S)          # floor
    assert spring_at(m, rec, (4, 3, 1))[:3] == (S, S, s)          # vertical corner
    assert spring_at(m, rec, (4, 1, 0))[:3] == (S, s, S)          # x-wall foot
    assert spring_at(m, rec, (0, 0, 0))[:3] == (S, S, S)          # bottom corner
    assert spring_at(m, rec, (0, 1, 3))[:3] == (S, s, s)          # free surface: no Z restraint
    assert all(k[3:] == (0, 0, 0) for k in (m.springs[p].k for p in m.springs))
    assert all(m.springs[p].damp == 0.0 for p in m.springs)        # OQ-3
    # one spring per coincident pair, wall node I / water node J
    for e in m.groups[int(rec.springs)].elements.values():
        assert np.linalg.norm(m.node_global(e.nodes[0]) - m.node_global(e.nodes[1])) < 1e-12
        assert e.nodes[1] > int(rec.offset)


def test_fillpool_empty_levels_and_volume(tmp_path):
    for empty in (0, 1, 2):
        ui = pool_ui(tmp_path, a=3.0, b=2.0, H=1.5, nx=3, ny=2, nz=3)
        ui.execute(f"FILLPOOL,1e9,0,{empty}")
        m = ui.model
        rec = wl.pool_record(m)
        assert water_volume(m, rec) == pytest.approx(3.0 * 2.0 * (1.5 - empty * 0.5), rel=1e-12)
        assert float(rec.zsurface) == pytest.approx(1.5 - empty * 0.5)


@pytest.mark.parametrize("cmd, text", [
    ("FILLPOOL,0", "<Stiff> must be > 0"),
    ("FILLPOOL,-5", "<Stiff> must be > 0"),
    ("FILLPOOL,1e6,0,0,-1,-1,-1", "<stiff2> must be >= 0"),
    ("FILLPOOL,1e6,0,-1", "<EmptyLevels> must be >= 0"),
    ("FILLPOOL,1e6,0,3", "leaves no water"),
    ("FILLPOOL,1e6,0,0,-1,5", "is smaller than the largest node number"),
])
def test_fillpool_argument_errors_leave_model_unchanged(tmp_path, cmd, text):
    ui = pool_ui(tmp_path)
    before = ui.model.canonical()
    ui.execute(cmd)
    assert any(text in e for e in errors(ui)), errors(ui)
    assert ui.model.canonical() == before


def test_fillpool_negative_sensitivity_warns(tmp_path):
    ui = pool_ui(tmp_path)
    ui.execute("FILLPOOL,1e9,-0.5")
    assert any("<Sensitivity> must be >= 0" in w for w in warnings(ui)) and not errors(ui)


def test_fillpool_model_errors(tmp_path):
    ui = make_ui(tmp_path)
    ui.execute("FILLPOOL")
    assert any("no nodes" in e for e in errors(ui))
    # a BEAMS element in the pool model
    ui = pool_ui(tmp_path)
    ui.run_text("N,900,9,9,9\nN,901,9,9,10\nN,902,10,9,9\nR,1,1,1,1,1,1,1\nGROUP,5,BEAMS\nE,1,900,901,902")
    ui.sink.clear()
    ui.execute("FILLPOOL")
    assert any("only the walls and floor" in e and "BEAMS" in e for e in errors(ui))
    # a vertical wall only: no floor
    ui = make_ui(tmp_path, ["N,1,0,0,0", "N,2,1,0,0", "N,3,1,0,1", "N,4,0,0,1", "M,1,3e10,0.2,0,0,0,1",
                            "GROUP,1,SHELL", "E,1,1,2,3,4", "THICK,1,1,1,0.5"])
    ui.execute("FILLPOOL")
    assert any("no horizontal floor" in e for e in errors(ui)), errors(ui)
    # excavated-soil SOLIDs
    ui = make_ui(tmp_path, solid_pool_lines(2.0, 2.0, 1.0, 2, 2, 2) + ["ETYPE,1,40,1,2"])
    ui.execute("FILLPOOL")
    assert any("ETYPE 2" in e for e in errors(ui)), errors(ui)


def test_fillpool_twice_is_refused_until_groups_deleted(tmp_path):
    ui = pool_ui(tmp_path)
    ui.execute("FILLPOOL,1e9")
    ui.sink.clear()
    ui.execute("FILLPOOL,1e9")
    assert any("SPRING" in e for e in errors(ui))
    rec = wl.pool_record(ui.model)
    ui.execute(f"GDEL,{rec.springs}")
    ui.sink.clear()
    ui.execute("FILLPOOL,1e9")
    assert any("already filled" in e for e in errors(ui))
    ui.execute(f"GDEL,{rec.water}")
    ui.sink.clear()
    ui.execute("FILLPOOL,1e9")
    assert not errors(ui)


def test_fillpool_floor_not_level_needs_sensitivity(tmp_path):
    ui = make_ui(tmp_path, shell_pool_lines(3.0, 2.0, 2.0, 3, 2, 2))
    m = ui.model
    n = 5                                    # floor node (0, 1, 0) lifted by 1 cm
    p = m.node_global(n)
    assert np.allclose(p, (0.0, 1.0, 0.0))
    ui.execute(f"N,{n},{p[0]},{p[1]},0.01")
    ui.sink.clear()
    ui.execute("FILLPOOL,1e9")
    assert any("not level" in e for e in errors(ui)), errors(ui)
    ui.sink.clear()
    ui.execute("FILLPOOL,1e9,0.05")
    assert not errors(ui)
    rec = wl.pool_record(m)
    # the lifted floor node is still attached (within Sensitivity)
    assert any(e.nodes[0] == n for e in m.groups[int(rec.springs)].elements.values())


def test_fillpool_sensitivity_merges_jittered_wall_levels(tmp_path):
    lines = shell_pool_lines(2.0, 2.0, 2.0, 2, 2, 2)
    ui = make_ui(tmp_path, lines)
    m = ui.model
    rng = np.random.default_rng(3)
    for n in sorted(m.nodes):
        p = m.node_global(n)
        if p[2] > 0:
            ui.execute(f"N,{n},{p[0]},{p[1]},{p[2] + rng.uniform(-0.004, 0.004)}")
    ui.sink.clear()
    ui.execute("FILLPOOL,1e9,0.02")
    assert not errors(ui)
    rec = wl.pool_record(m)
    assert water_volume(m, rec) == pytest.approx(2.0 * 2.0 * float(rec.zsurface), rel=1e-12)
    assert abs(float(rec.zsurface) - 2.0) < 0.004
    # every wall node at the water levels has its spring: floor 9 + 8 perimeter x 2 levels
    assert len(m.groups[int(rec.springs)].elements) == 9 + 8 * 2


def test_fillpool_oblique_walls_use_projected_constants(tmp_path):
    ui = pool_ui(tmp_path, a=2.0, b=2.0, H=1.0, nx=2, ny=2, nz=1)
    ui.execute("ROTATE,0,0,0,30,0,0")
    ui.sink.clear()
    ui.execute("FILLPOOL,1e6,0,0,-1,0,0")
    assert not errors(ui)
    assert any("oblique walls" in w for w in warnings(ui))
    m = ui.model
    rec = wl.pool_record(m)
    c, s = math.cos(math.radians(30)), math.sin(math.radians(30))
    # mid node of the wall x' = 0 (normal -x' rotated by 30 deg): diagonal of Stiff n n^T
    p = np.array([0.0, 1.0, 1.0])
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    k = spring_at(m, rec, R @ p)
    assert k[0] == pytest.approx(1e6 * c * c, rel=1e-9) and k[1] == pytest.approx(1e6 * s * s, rel=1e-9)
    assert k[2] == 0.0
    # floor springs stay exact (normal Z)
    assert spring_at(m, rec, R @ np.array([1.0, 1.0, 0.0]))[:3] == (0.0, 0.0, 1e6)


def test_fillpool_triangular_floor_gives_prisms(tmp_path):
    # a right-triangle pool: floor = 4 triangles (repeated 4th node), walls split at the edge midpoints
    pts = [(0, 0), (1, 0), (2, 0), (1, 1), (0, 2), (0, 1)]          # A D B E C F
    lines = [f"N,{k + 1},{x},{y},0" for k, (x, y) in enumerate(pts)]
    lines += [f"N,{k + 7},{x},{y},1" for k, (x, y) in enumerate(pts)]
    lines += ["M,1,3e10,0.2,0,0,0,1", "GROUP,1,SHELL", "E,1,1,2,6,6", "E,2,2,3,4,4", "E,3,6,4,5,5", "E,4,2,4,6,6"]
    lines += [f"E,{5 + k},{k + 1},{(k + 1) % 6 + 1},{(k + 1) % 6 + 7},{k + 7}" for k in range(6)]
    lines += ["THICK,1,10,1,0.3"]
    ui = make_ui(tmp_path, lines)
    ui.execute("FILLPOOL,1e6")
    assert not errors(ui), errors(ui)
    m = ui.model
    rec = wl.pool_record(m)
    els = m.groups[int(rec.water)].elements
    assert len(els) == 4 and all(len(set(e.nodes)) == 6 for e in els.values())   # prisms by repeated nodes
    assert water_volume(m, rec) == pytest.approx(2.0, rel=1e-12)
    # the hypotenuse mid node E (normal (1,1)/sqrt 2) is oblique: diagonal of Stiff n n^T; the corners
    # (two normals spanning the horizontal plane) and the floor are exact
    assert any("oblique" in w for w in warnings(ui))
    assert spring_at(m, rec, (1, 1, 1))[:3] == pytest.approx((0.5e6, 0.5e6, 0.0))
    assert spring_at(m, rec, (2, 0, 1))[:3] == (1e6, 1e6, 0.0)
    assert spring_at(m, rec, (1, 0, 1))[:3] == (0.0, 1e6, 0.0)


def test_fillpool_unconnected_water_nodes_warn(tmp_path):
    # the y = 1 wall is one element from x = 0 to 2: the water node (1, 1, 1) has no wall node
    lines = ["N,1,0,0,0", "N,2,1,0,0", "N,3,2,0,0", "N,4,0,1,0", "N,5,1,1,0", "N,6,2,1,0",
             "N,7,0,0,1", "N,8,1,0,1", "N,9,2,0,1", "N,10,0,1,1", "N,11,2,1,1",
             "M,1,3e10,0.2,0,0,0,1", "GROUP,1,SHELL", "E,1,1,2,5,4", "E,2,2,3,6,5", "E,3,1,2,8,7", "E,4,2,3,9,8",
             "E,5,4,6,11,10", "E,6,1,4,10,7", "E,7,3,6,11,9", "THICK,1,7,1,0.3"]
    ui = make_ui(tmp_path, lines)
    ui.execute("FILLPOOL,1e9")
    assert not errors(ui), errors(ui)
    w = [t for t in warnings(ui) if "no coincident wall/floor node" in t]
    assert w and w[0].startswith("FILLPOOL: 1 water nodes"), warnings(ui)
    m = ui.model
    rec = wl.pool_record(m)
    assert len(m.groups[int(rec.springs)].elements) == 6 + 5     # floor 6 + perimeter 6 - 1 at z = 1


def test_fillpool_solid_pool(tmp_path):
    ui = make_ui(tmp_path, solid_pool_lines(3.0, 2.0, 2.0, 3, 2, 4))
    ui.execute("FILLPOOL,1e9,0,1,1")
    assert not errors(ui)
    m = ui.model
    rec = wl.pool_record(m)
    assert water_volume(m, rec) == pytest.approx(9.0, rel=1e-12)
    walls = {e.nodes[0] for e in m.groups[int(rec.springs)].elements.values()}
    assert all(m.nodes[w].fix[3:6] == [1, 1, 1] for w in walls)       # SOLID wall nodes: rotations fixed
    shells = m.groups[int(rec.shells)]
    area = sum(_area(m, int(rec.shells), e) for e in shells.elements)
    assert area == pytest.approx(3 * 2 + 2 * (3 + 2) * 1.5, rel=1e-12)


def _area(m, g, e):
    X = element_xyz(m, g, e)
    return 0.5 * np.linalg.norm(np.cross(X[2] - X[0], X[3] - X[1]))


def test_fillpool_check_after_fixrot_has_no_unrestrained_rotation(tmp_path):
    ui = pool_ui(tmp_path)
    ui.run_text("MDL,pool,.\nFILLPOOL,1e9,0,1,1\nFIXROT\nINT,1,15,1,1\nCHECK")
    err = (tmp_path / "pool.err").read_text()
    assert "EDU-06" not in err and "MODEL: 0 errors" in err


def test_fillpool_mopt_kept_when_already_zero(tmp_path):
    ui = pool_ui(tmp_path)
    ui.execute("MOPT,0,1,0,1")
    ui.sink.clear()
    ui.execute("FILLPOOL")
    assert not any("MOPT" in w for w in warnings(ui))
    mo = ui.model.mopt
    assert (mo.incomp, mo.get("matrix"), mo.get("mass"), mo.get("force")) == (0, 1, 0, 1)


def test_fillpool_write_inp_round_trip_and_pooldata(tmp_path):
    ui = pool_ui(tmp_path)
    ui.execute("FILLPOOL,2e8,0,1,1,500,10")
    m = ui.model
    text, notes = write_pre(m)
    assert not notes
    assert "* FILLPOOL pool data (SASSI-EDU)" in text and "\nPOOLDATA,2,3,4,200000000,0,1,1,500,10,2,0,2\n" in text
    ui2 = Interpreter(cwd=tmp_path)
    ui2.run_text(text)
    assert not ui2.sink.texts(Kind.ERROR) and m.same_state(ui2.model)
    assert write_pre(ui2.model)[0] == text
    ui2.execute("POOLDATA")
    assert wl.pool_record(ui2.model) is None
    ui2.execute("POOLDATA,2,3,4")
    assert int(wl.pool_record(ui2.model).springs) == 3


# --------------------------------------------------------------------------------------
# LISTPOOLINTER and MERGEPOOL
# --------------------------------------------------------------------------------------
def original_and_pool(tmp_path, solid=False, offset=None, extra=500):
    """Original model (model 0): a pool (group 1) and a column (BEAMS, nodes extra..extra+2); pool sub-model 2
    (CUT2SUB of group 1) filled by FILLPOOL with ``offset`` (default: the last node of the original)."""
    lines = solid_pool_lines(2.0, 2.0, 1.0, 2, 2, 2) if solid else shell_pool_lines(4.0, 2.0, 2.0, 4, 2, 2)
    a, b, k = extra, extra + 1, extra + 2
    lines += [f"N,{a},6,0,0", f"N,{b},6,0,3", f"N,{k},7,0,1.5", "M,2,3e7,0.2,24,0.05,0.05",
              "R,1,0.25,0.2,0.2,0.01,0.005,0.005", "GROUP,2,BEAMS", "MACT,2", "RACT,1", f"E,1,{a},{b},{k}"]
    ui = make_ui(tmp_path, lines)
    last = max(ui.model.nodes)
    n1 = len(ui.model.groups[1].elements)
    off = last if offset is None else offset
    ui.run_text(f"CUTADD,1,1,RANGE,1,{n1}\nCUT2SUB,1,2\nACTM,2\nFILLPOOL,1e9,0,0,1,{off},0\nACTM,0")
    assert not errors(ui), errors(ui)
    ui.sink.clear()
    return ui, last


def test_listpoolinter_lists_matching_nodes(tmp_path):
    ui, _ = original_and_pool(tmp_path)
    ui.execute("LISTPOOLINTER,2")
    pool = ui.models[2]
    walls = sorted({w for w, _ in wl.pool_interface(pool)})
    rows = [t for t in ui.sink.texts(Kind.INFO) if t.split() and t.split()[0].isdigit()]
    assert [int(t.split()[0]) for t in rows] == walls
    assert f"{len(walls)} of {len(walls)} interface nodes" in ui.sink.texts(Kind.CONFIRM)[-1]
    assert not warnings(ui)


def test_listpoolinter_reports_moved_and_missing_nodes(tmp_path):
    ui, _ = original_and_pool(tmp_path)
    walls = sorted({w for w, _ in wl.pool_interface(ui.models[2])})
    p = ui.model.node_global(walls[0])
    ui.execute(f"N,{walls[0]},{p[0] + 0.5},{p[1]},{p[2]}")
    ui.execute(f"NDEL,{walls[1]}")
    ui.sink.clear()
    ui.execute("LISTPOOLINTER,2")
    w = " ".join(warnings(ui))
    assert "another position" in w and f"{walls[0]} (distance 0.5)" in w
    assert "not defined in the active model" in w and str(walls[1]) in w
    assert f"{len(walls) - 2} of {len(walls)} interface nodes" in ui.sink.texts(Kind.CONFIRM)[-1]
    ui.sink.clear()
    ui.execute("MERGEPOOL,2")
    assert any("LISTPOOLINTER" in e for e in errors(ui))


def test_listpoolinter_errors(tmp_path):
    ui, _ = original_and_pool(tmp_path)
    ui.execute("LISTPOOLINTER,7")
    assert any("not in memory" in e for e in errors(ui))
    ui.sink.clear()
    ui.execute("LISTPOOLINTER,0")
    assert any("is the active model" in e for e in errors(ui))
    ui.sink.clear()
    ui.execute("CPMODEL,3")                    # an unfilled copy of the original
    ui.sink.clear()
    ui.execute("LISTPOOLINTER,3")
    assert any("not filled with FILLPOOL" in e for e in errors(ui))
    ui.sink.clear()
    pool = ui.models[2]
    rec = wl.pool_record(pool)
    del pool.groups[int(rec.springs)]
    ui.execute("LISTPOOLINTER,2")
    assert any("no longer in the pool model" in e for e in errors(ui))


def test_mergepool_imports_with_node_numbers_kept(tmp_path):
    ui, last = original_and_pool(tmp_path)
    m = ui.model
    pool = ui.models[2]
    prec = wl.pool_record(pool)
    g0, mat0, sc0 = max(m.groups), max(m.materials), max(m.springs) if m.springs else 0
    ui.execute("MERGEPOOL,2")
    assert not errors(ui), errors(ui)
    assert any("MOPT <incomp> set to 0" in w for w in warnings(ui))
    assert int(m.mopt.get("incomp")) == 0
    new = sorted(g for g in m.groups if g > g0)
    assert [m.groups[g].type for g in new] == [1, 7, 3]
    assert [m.groups[g].title for g in new] == [wl.TITLES["water"], wl.TITLES["springs"], wl.TITLES["shells"]]
    # same element numbers and nodes as in the pool model
    for g_pool, g_new in zip((int(prec.water), int(prec.springs), int(prec.shells)), new):
        assert {k: list(e.nodes) for k, e in pool.groups[g_pool].elements.items()} == \
               {k: list(e.nodes) for k, e in m.groups[g_new].elements.items()}
    water = m.groups[new[0]]
    assert {e.mat for e in water.elements.values()} == {mat0 + 1}
    assert m.materials[mat0 + 1].val1 == pool.materials[int(prec.watermat)].val1
    assert min(e.prop for e in m.groups[new[1]].elements.values()) == sc0 + 1
    wn = sorted({n for e in water.elements.values() for n in e.nodes})
    assert wn[0] == last + 1
    assert all(np.allclose(m.node_global(n), pool.node_global(n)) for n in wn)
    attached = {e.nodes[1] for e in m.groups[new[1]].elements.values()}
    assert all(m.nodes[n].fix[3:6] == [1, 1, 1] for n in attached)
    # the original model is still consistent: CHECK of the model data has no error
    ui.sink.clear()
    ui.run_text("MDL,orig,.\nFIXROT\nCHECK")
    assert "MODEL: 0 errors" in (tmp_path / "orig.err").read_text()
    # the merged model survives WRITE -> INP (UT-03)
    text, notes = write_pre(m)
    assert not notes and "POOLDATA" not in text
    ui3 = Interpreter(cwd=tmp_path)
    ui3.run_text(text)
    assert not ui3.sink.texts(Kind.ERROR) and m.same_state(ui3.model)
    # importing again is refused
    ui.sink.clear()
    ui.execute("MERGEPOOL,2")
    assert any("already used" in e for e in errors(ui))


def test_mergepool_requires_offset_beyond_original_nodes(tmp_path):
    # the column nodes 40-42 follow the 39 pool nodes; FILLPOOL offset 0 numbers the water from 40 too
    ui, last = original_and_pool(tmp_path, offset=0, extra=40)
    assert last == 42
    ui.execute("MERGEPOOL,2")
    assert any("already used" in e and f">= {last}" in e for e in errors(ui)), errors(ui)
    # without an actual clash the import works whatever the offset
    ui, last = original_and_pool(tmp_path, offset=0, extra=500)
    ui.execute("MERGEPOOL,2")
    assert not errors(ui)


def test_mergepool_fixes_rotations_of_solid_wall_nodes(tmp_path):
    ui, _ = original_and_pool(tmp_path, solid=True)
    m = ui.model
    walls = sorted({w for w, _ in wl.pool_interface(ui.models[2])})
    assert not any(m.nodes[w].fix[3:6] == [1, 1, 1] for w in walls)
    ui.execute("MERGEPOOL,2")
    assert not errors(ui)
    assert all(m.nodes[w].fix[3:6] == [1, 1, 1] for w in walls)


# --------------------------------------------------------------------------------------
# water_lib helpers and analytical references
# --------------------------------------------------------------------------------------
def test_impulsive_ratio_exact_limits_and_identity():
    assert wl.impulsive_ratio_exact(1.0, 1.0) == pytest.approx(0.5, abs=1e-14)
    for r in (0.05, 0.3, 0.7, 2.0, 7.5):
        assert wl.impulsive_ratio_exact(r, 1.0) + wl.impulsive_ratio_exact(1.0, r) == pytest.approx(1.0, abs=1e-12)
    # brute-force series
    for L, H in ((2.0, 1.0), (0.4, 1.3)):
        n = np.arange(400000)
        lam = (2 * n + 1) * np.pi / (2 * H)
        brute = 2 / (L * H * H) * np.sum(np.tanh(lam * L) / lam ** 3)
        assert wl.impulsive_ratio_exact(L, H) == pytest.approx(brute, rel=1e-9)
    # narrow tank -> all water impulsive; long tank -> Westergaard's exact series 0.5428 H/L
    assert wl.impulsive_ratio_exact(1e-3, 1.0) == pytest.approx(1.0, abs=1e-3)
    k = 28 * 1.2020569031595942 / math.pi ** 3 / 2          # 0.54275...
    assert wl.impulsive_ratio_exact(100.0, 1.0) * 100.0 == pytest.approx(k, rel=1e-12)
    with pytest.raises(ValueError):
        wl.impulsive_ratio_exact(0.0, 1.0)


def test_housner_westergaard_and_heights():
    assert wl.housner_impulsive_ratio(1.0, 1.0) == pytest.approx(math.tanh(math.sqrt(3)) / math.sqrt(3))
    # Housner is a conservative approximation of the exact series (within ~11 %)
    for r in (0.25, 0.5, 1.0, 2.0, 4.0):
        e, h = wl.impulsive_ratio_exact(r, 1.0), wl.housner_impulsive_ratio(r, 1.0)
        assert e < h < 1.12 * e
    assert wl.westergaard_ratio(10.0, 1.0) == pytest.approx(7 / 120)
    assert wl.impulsive_height_exact(1.0, 1.0) == pytest.approx(0.40467, abs=2e-5)
    assert wl.impulsive_height_exact(50.0, 1.0) == pytest.approx(0.40142, abs=2e-5)  # long reservoir
    assert wl.HOUSNER_HEIGHT_EBP == 0.375


def test_housner_convective_and_sloshing():
    g = 9.81
    mc, fc, hc = wl.housner_convective(5.0, 5.0, g)
    assert mc == pytest.approx(math.sqrt(2.5) / 3 * math.tanh(math.sqrt(2.5)), rel=1e-12)
    assert fc == pytest.approx(wl.sloshing_frequency_exact(5.0, 5.0, g), rel=0.01)
    assert 0.5 < hc < 0.6
    # exact impulsive mass + first convective mode ~ the whole water (the first sloshing mode dominates)
    for r in (0.5, 1.0, 2.0):
        assert wl.impulsive_ratio_exact(r, 1.0) + wl.housner_convective(r, 1.0, g)[0] == pytest.approx(1.0, abs=0.05)
    assert wl.sloshing_frequency_exact(5.0, 5.0, g) == pytest.approx(
        math.sqrt(math.pi * g / 10 * math.tanh(math.pi / 2)) / (2 * math.pi), rel=1e-12)


def test_normal_projector_and_interface_constants():
    P = wl.normal_projector([np.array([1.0, 0, 0]), np.array([0, 0, -1.0])])
    assert np.allclose(P, np.diag([1, 0, 1]))
    k, exact = wl.interface_constants(P, 10.0, 2.0)
    assert exact and k == (10.0, 2.0, 10.0)
    n = np.array([math.cos(0.3), math.sin(0.3), 0.0])
    k, exact = wl.interface_constants(wl.normal_projector([n]), 10.0, 0.0)
    assert not exact and k == pytest.approx((10 * n[0] ** 2, 10 * n[1] ** 2, 0.0))
    assert np.allclose(wl.normal_projector([n, n, -n]), np.outer(n, n))       # repeated normals: rank 1
    assert np.allclose(wl.normal_projector([]), 0.0)


def test_boundary_normals_corners_and_smooth_walls():
    xy = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    nrm, edges = wl._boundary_normals(xy, [(0, 1, 2, 3)])
    assert len(edges) == 4
    got = sorted(tuple(np.round(v, 12)) for v in nrm[0])
    assert got == [(-1.0, 0.0, 0.0), (0.0, -1.0, 0.0)]                         # corner: two normals
    # a 16-gon: neighbouring edges differ by 22.5 deg -> one averaged (radial) normal per vertex
    t = np.linspace(0, 2 * np.pi, 17)[:-1]
    xy = np.column_stack([np.cos(t), np.sin(t)])
    xy = np.vstack([xy, [0.0, 0.0]])
    cells = [(16, k, (k + 1) % 16) for k in range(16)]
    nrm, _ = wl._boundary_normals(xy, cells)
    for k in range(16):
        assert len(nrm[k]) == 1 and np.allclose(nrm[k][0][:2], xy[k], atol=1e-12)


def test_resolve_offset(tmp_path):
    ui = pool_ui(tmp_path)
    m = ui.model
    nmax = max(m.nodes)
    assert wl.resolve_offset(m, -1) == nmax and wl.resolve_offset(m, 0) == nmax
    assert wl.resolve_offset(m, nmax) == nmax and wl.resolve_offset(m, nmax + 7) == nmax + 7
    with pytest.raises(wl.GenerationError):
        wl.resolve_offset(m, nmax - 1)


def test_fillpool_roof_touching_the_water_is_reported(tmp_path):
    # a closed box: a roof of SHELLs at the top level; full (EmptyLevels = 0) -> the surface nodes coincide
    # with roof nodes but stay free (reported); one empty level -> nothing to report
    lines = shell_pool_lines(2.0, 2.0, 2.0, 2, 2, 2)
    roof = ["GROUP,2,SHELL", "MACT,1"]
    ui = make_ui(tmp_path, lines)
    m = ui.model
    top = {tuple(np.round(m.node_global(n)[:2], 9)): n for n in m.nodes if abs(m.node_global(n)[2] - 2.0) < 1e-9}
    k = 0
    ids = []
    for j in range(2):
        for i in range(2):
            corners = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
            if all(c in top for c in corners):
                k += 1
                ids.append(f"E,{k}," + ",".join(str(top[c]) for c in corners))
    centre_missing = (1.0, 1.0) not in top
    if centre_missing:                       # the walls have no node at the roof centre: add it
        ui.execute(f"N,{max(m.nodes) + 1},1,1,2")
        top[(1.0, 1.0)] = max(m.nodes)
        ids = []
        for j in range(2):
            for i in range(2):
                corners = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
                ids.append(f"E,{len(ids) + 1}," + ",".join(str(top[c]) for c in corners))
    ui.run_text("\n".join(roof + ids + [f"THICK,1,{len(ids)},1,0.3"]))
    ui.sink.clear()
    ui.execute("FILLPOOL,1e9,0,0")
    assert not errors(ui), errors(ui)
    assert any("nodes of the water surface coincide with wall/roof nodes" in w for w in warnings(ui)), warnings(ui)
    rec = wl.pool_record(m)
    ui.execute(f"GDEL,{rec.water},{rec.springs}")
    ui.execute("POOLDATA")
    ui.sink.clear()
    ui.execute("FILLPOOL,1e9,0,1")
    assert not errors(ui) and not any("water surface" in w for w in warnings(ui))
