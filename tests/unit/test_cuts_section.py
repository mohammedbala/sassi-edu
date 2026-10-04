"""Section-cut geometry and resultants: pieces, section properties, resultant operators, CALCPAR, CALCMOI,
CALCC, CALCM (requirements 4.14, UT-21; spec 10 sections 1.4, 3.2-3.8; D-SEC-01 ... D-SEC-04)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.model.geometry import loc_matrix
from sassi.prep import Interpreter, Kind
from sassi.prep import cuts_lib as cl
from sassi.prep.writer import write_pre
from sassi.verify.problems.vp_cuts import shell_wall_lines, two_hex_lines

from .test_cuts_commands import grid_lines


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def run(ui, lines):
    ui.run_text("\n".join(lines) if isinstance(lines, list) else lines)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    return ui.model


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


def all_keys(m):
    return [(g.id, e.id) for g, e in m.iter_elements()]


def section(m, P, n, r, table=None, bending=True):
    """Properties and local resultants of every element of ``m`` cut by the plane (lib level)."""
    ex, ey, ez = cl.section_axes(n, r)
    sp = cl.section_pieces(m, all_keys(m), P, ez)
    props = cl.section_properties(sp.pieces, P, ex, ey, ez)
    R = np.zeros(6)
    if table is not None:
        T = cl.resultant_operators(sp.pieces, props, bending)
        S, ok = cl.gather_stresses(sp.pieces, table)
        R = cl.resultants(T, S, ok)
    return sp, props, R


def voigt(S):
    return [S[0, 0], S[1, 1], S[2, 2], S[0, 1], S[0, 2], S[1, 2]]


def tensor(s):
    return np.array([[s[0], s[3], s[4]], [s[3], s[1], s[5]], [s[4], s[5], s[2]]], float)


# ======================================================================================
# Local axes (D-SEC-01, spec 10 section 1.4)
# ======================================================================================
def test_section_axes():
    ex, ey, ez = cl.section_axes([0, 0, 2], [1, 0, 1])
    assert np.allclose([ex, ey, ez], np.eye(3))
    ex, ey, ez = cl.section_axes([1, 0, 0], [0, 1, 0])
    assert np.allclose(ex, [0, 1, 0]) and np.allclose(ey, [0, 0, 1]) and np.allclose(np.cross(ex, ey), ez)
    with pytest.raises(cl.SectionError, match="parallel"):
        cl.section_axes([0, 0, 1], [0, 0, 3])
    with pytest.raises(cl.SectionError):
        cl.section_axes([0, 0, 0], [1, 0, 0])


# ======================================================================================
# Pieces and section properties
# ======================================================================================
def test_regular_hexagon_section_of_a_cube(ui):
    """Plane x + y + z = 1.5 through the unit cube: regular hexagon of side sqrt(2)/2,
    A = 3 sqrt(3)/4, Ixx = Iyy = 5 sqrt(3)/64 (isotropic), centroid (0.5, 0.5, 0.5)."""
    m = run(ui, grid_lines(1, 1, 1))
    sp, props, _ = section(m, [0.5, 0.5, 0.5], [1, 1, 1], [1, -1, 0])
    assert len(sp.pieces) == 1 and len(sp.pieces[0].pts) == 6
    assert props.area == pytest.approx(3 * math.sqrt(3) / 4, rel=1e-12)
    assert np.allclose(props.C, 0.5, atol=1e-12)
    s4 = 0.25
    assert props.Ixx == pytest.approx(5 * math.sqrt(3) / 16 * s4, rel=1e-12)
    assert props.Iyy == pytest.approx(props.Ixx, rel=1e-12) and abs(props.Ixy) < 1e-14


def test_mesh_refinement_and_uniform_stress_equilibrium(ui):
    """A uniform stress state transmits F = A S.n through any plane and no moment about the centroid; the
    section of a 2 x 2 x 2 mesh equals the section of the single cube."""
    S = np.array([[3.0, 1.5, -2.0], [1.5, -4.0, 0.7], [-2.0, 0.7, 6.0]])
    n, r, P = np.array([1.0, 2.0, 3.0]), [0, 0, 1], [0.45, 0.5, 0.52]
    m1 = run(ui, grid_lines(1, 1, 1))
    _, p1, _ = section(m1, P, n, r)
    ui2 = Interpreter()
    run(ui2, grid_lines(2, 2, 2, 0.5))
    m8 = ui2.model
    table = {k: np.array(voigt(S)) for k in all_keys(m8)}
    sp, p8, R = section(m8, P, n, r, table)
    assert len(sp.pieces) > 1
    assert p8.area == pytest.approx(p1.area, rel=1e-12)
    assert np.allclose(p8.C, p1.C, atol=1e-12)
    assert np.allclose([p8.Ixx, p8.Iyy, p8.Ixy], [p1.Ixx, p1.Iyy, p1.Ixy], atol=1e-12)
    E = np.vstack([p8.ex, p8.ey, p8.ez])
    assert np.allclose(R[:3], E @ (p8.area * S @ p8.ez), atol=1e-12)
    assert np.allclose(R[3:], 0.0, atol=1e-12)


def test_t_s1_is_invariant_under_a_rotation_of_everything(ui):
    """Rotating the model, the stress tensors, n and r leaves the local properties and resultants unchanged."""
    m = run(ui, two_hex_lines())
    st = {(1, 1): np.array([1.0, -2, 10, 0.5, 3, 1]), (1, 2): np.array([2.0, 1, -10, -0.5, 3, -1])}
    _, p0, R0 = section(m, [0, 0, 0.5], [0, 0, 1], [1, 0, 0], st)
    Q = loc_matrix(30.0, 20.0, -50.0)
    ui2 = Interpreter()
    run(ui2, ["LOC,1,0,7,-3,2,30,20,-50", "CSYS,1"] + two_hex_lines())
    m2 = ui2.model
    st2 = {k: np.array(voigt(Q @ tensor(v) @ Q.T)) for k, v in st.items()}
    P2 = np.array([7.0, -3, 2]) + Q @ np.array([0, 0, 0.5])
    _, p2, R2 = section(m2, P2, Q @ [0, 0, 1], Q @ [1, 0, 0], st2)
    assert np.allclose([p2.area, p2.Ixx, p2.Iyy, p2.Ixy], [p0.area, p0.Ixx, p0.Iyy, p0.Ixy], atol=1e-12)
    assert np.allclose(R2, R0, atol=1e-11)
    assert np.allclose(R0, [6.0, 0.0, 0.0, 0.0, 10.0, -1.0], atol=1e-12)


def test_coincident_faces_counted_once(ui):
    """D-SEC-04: a plane on the face shared by two stacked solids counts it once (the -n side element)."""
    m = run(ui, grid_lines(1, 1, 2))
    ex, ey, ez = cl.section_axes([0, 0, 1], [1, 0, 0])
    sp = cl.section_pieces(m, [(1, 1), (1, 2)], [0, 0, 1], ez)
    assert [p.key for p in sp.pieces] == [(1, 1)] and sp.coincident == 1
    sp = cl.section_pieces(m, [(1, 2)], [0, 0, 1], ez)
    assert [p.key for p in sp.pieces] == [(1, 2)]
    sp = cl.section_pieces(m, [(1, 1), (1, 2)], [0, 0, 1], -ez)
    assert [p.key for p in sp.pieces] == [(1, 2)]
    # shells: an edge shared by two stacked walls
    ui2 = Interpreter()
    run(ui2, ["N,1,0,0,0", "N,2,1,0,0", "N,3,1,0,1", "N,4,0,0,1", "N,5,1,0,2", "N,6,0,0,2",
              "GROUP,1,SHELL", "E,1,1,2,3,4", "E,2,4,3,5,6", "THICK,1,2,1,0.2"])
    sp = cl.section_pieces(ui2.model, [(1, 1), (1, 2)], [0, 0, 1], ez)
    assert [p.key for p in sp.pieces] == [(1, 1)]
    props = cl.section_properties(sp.pieces, [0, 0, 1], ex, ey, ez)
    assert props.area == pytest.approx(0.2)


def test_oblique_cut_of_a_wall_transmits_the_wall_force(ui):
    """Wall y = 0 with uniform vertical stress sigma (Sy'y', y' = Z): a plane tilted towards the wall normal
    cuts a width t/cos(a) but the transmitted force is still sigma t L along Z."""
    m = run(ui, shell_wall_lines(0.5))
    sig, a = 100.0, math.radians(35.0)
    n = [0.0, math.sin(a), math.cos(a)]
    table = {(1, 1): np.array([0, sig, 0, 0, 0, 0.0]), (1, 2): np.array([0, sig, 0, 0, 0, 0.0])}
    sp, props, R = section(m, [0, 0, 1], n, [1, 0, 0], table)
    assert props.area == pytest.approx(4 * 0.5 / math.cos(a), rel=1e-12)
    F = np.vstack([props.ex, props.ey, props.ez]).T @ R[:3]
    assert np.allclose(F, [0, 0, sig * 0.5 * 4], atol=1e-10)


def test_shell_bending_couple_matches_a_layered_solid_and_sectbend_switch(ui):
    """Plate z = 0 (normal +Z) with Mx'x' = M: the couple across x = 1 equals the moment of the
    equivalent two-layer SOLID plate with sigma_xx = +-s0 (M = s0 t^2/4); EDUOPT,SECTBEND,0 removes it."""
    t, s0, W = 0.4, 50.0, 1.0
    M = s0 * t * t / 4
    m = run(ui, ["N,1,0,0,0", "N,2,2,0,0", "N,3,2,1,0", "N,4,0,1,0", "GROUP,1,SHELL", "E,1,1,2,3,4",
                 f"THICK,1,1,1,{t}"])
    table = {(1, 1): np.array([0, 0, 0, M, 0, 0.0])}
    _, ps, Rs = section(m, [1, 0, 0], [1, 0, 0], [0, 1, 0], table)
    ui2 = Interpreter()
    run(ui2, _two_layer_plate(t))
    st = {(1, 1): np.array([-s0, 0, 0, 0, 0, 0.0]), (1, 2): np.array([s0, 0, 0, 0, 0, 0.0])}
    _, pv, Rv = section(ui2.model, [1, 0, 0], [1, 0, 0], [0, 1, 0], st)
    assert ps.area == pytest.approx(pv.area) == pytest.approx(t * W)
    assert np.allclose(Rs, Rv, atol=1e-12)
    assert Rs[3] == pytest.approx(M * W)                       # Mx (ex = Y): global moment +Y
    _, _, R0 = section(m, [1, 0, 0], [1, 0, 0], [0, 1, 0], table, bending=False)
    assert np.allclose(R0, 0.0)
    ui.execute("EDUOPT,SECTBEND,0")
    assert not cl.bending_on(ui.model)


def _two_layer_plate(t):
    """2 x 1 plate of thickness t as two SOLID layers (z in [-t/2, 0] and [0, t/2])."""
    zs = [-t / 2, 0.0, t / 2]
    out = []
    for k, z in enumerate(zs):
        for j, (x, y) in enumerate([(0, 0), (2, 0), (2, 1), (0, 1)]):
            out.append(f"N,{4 * k + j + 1},{x},{y},{z}")
    out.append("GROUP,1,SOLID")
    for k in range(2):
        b = [4 * k + j + 1 for j in range(4)]
        out.append(f"E,{k + 1}," + ",".join(str(x) for x in b + [q + 4 for q in b]))
    return out


def test_plane_strain_strip(ui):
    m = run(ui, ["N,1,0,0,0", "N,2,2,0,0", "N,3,2,0,1", "N,4,0,0,1", "GROUP,1,PLANE", "E,1,1,2,3,4"])
    table = {(1, 1): np.array([7.0, 5.0, 2.0, 0, 0, 0])}       # Sxx, Szz, Txz
    sp, props, R = section(m, [0, 0, 0.5], [0, 0, 1], [1, 0, 0], table)
    assert props.area == pytest.approx(2.0)                     # unit plane-strain thickness
    assert np.allclose(R, [4.0, 0, 10.0, 0, 0, 0])
    sp = cl.section_pieces(m, [(1, 1)], [0, 0, 0.5], [0, 0.6, 0.8])
    assert not sp.pieces and "PLANE (the plane normal must lie in the X-Z plane)" in sp.skipped


# ======================================================================================
# CALCPAR / CALCMOI on CSECT models
# ======================================================================================
def _ts1(ui, tmp_path, vals=None):
    run(ui, two_hex_lines())
    vals = vals or {(1, 1): [0, 0, 10, 0, 3, 0], (1, 2): [0, 0, -10, 0, 3, 0]}
    cl.write_ess(tmp_path / "s.ess", {1: "SOLID"}, vals)
    run(ui, ["READSTR,s.ess", "CUTADD,1,1,1-2", "CSECT,5,1,0,0,0.5,0,0,1", "ACTM,5"])


def test_calcpar_verbose_one_line_and_sysno(ui, tmp_path):
    _ts1(ui, tmp_path)
    assert ui.execute("CALCPAR,0,0,1,1,0,0,11,1")
    line = ui.sink.texts(Kind.INFO)[-1].split()
    assert len(line) == 14 and float(line[8]) == pytest.approx(6.0) and float(line[12]) == pytest.approx(10.0)
    cs = ui.model.csys[11]
    assert cs.kind == "LOC" and np.allclose(cs.origin, [1, 0.5, 0.5]) and np.allclose(cs.R, np.eye(3))
    assert ui.model.csys_active == 0
    text, _ = write_pre(ui.model)
    b = Interpreter()
    b.run_text(text)
    assert b.model.same_state(ui.model)


def test_calcpar_normal_rules_on_a_csect_model(ui, tmp_path):
    _ts1(ui, tmp_path)
    assert not ui.execute("CALCPAR,1,0,0,0,1,0")                # not the CSECT normal
    assert "CSECT plane normal" in errors(ui)[-1]
    assert ui.execute("CALCPAR,0,0,-1,1,0,0")                   # the other sense: view from the other side
    r = ui.session["calcpar"]
    assert r["Fz"] == pytest.approx(0.0, abs=1e-12) and r["My"] == pytest.approx(10.0)
    assert r["Fx"] == pytest.approx(-6.0)                       # in-plane shear seen from the other side


def test_calcpar_without_stresses_and_missing_records(ui, tmp_path):
    run(ui, two_hex_lines() + ["CUTADD,1,1,1-2", "CSECT,5,1,0,0,0.5,0,0,1", "ACTM,5", "CALCPAR,0,0,1,1,0,0"])
    assert any("no element stresses are loaded" in w for w in warnings(ui))
    assert ui.session["calcpar"]["Area"] == pytest.approx(2.0) and ui.session["calcpar"]["My"] == 0.0
    ui2 = Interpreter(cwd=tmp_path)
    _ts1(ui2, tmp_path, {(1, 1): [0, 0, 10, 0, 0, 0], (1, 2): [0, 0, 0, 0, 0, 0]})
    ui2.model.ui_state["element_stress"]["rows"].pop()            # element 2 without data
    assert ui2.execute("CALCPAR,0,0,1,1,0,0")
    assert any("without stress data" in w for w in ui2.sink.texts(Kind.WARNING))
    assert ui2.session["calcpar"]["Fz"] == pytest.approx(10.0)


def test_calcmoi_after_write_inp_uses_the_extruded_geometry(ui, tmp_path):
    """Without the stored pieces (WRITE -> INP), CALCMOI intersects the unit-thickness elements with their
    mid-plane: same properties for solids (hexagon split in two), shells (t_eff) and the T-S2 wall."""
    run(ui, grid_lines(1, 1, 1) + ["CUTADD,1,1,1", "CSECT,2,1,0.5,0.5,0.5,1,1,1", "ACTM,2",
                                   "CALCMOI,1,1,1,1,-1,0"])
    ref = dict(ui.session["calcmoi"])
    assert ref["Area"] == pytest.approx(3 * math.sqrt(3) / 4)
    text, _ = write_pre(ui.model)
    b = Interpreter()
    b.run_text(text)
    assert "csect" not in b.model.ui_state
    assert b.execute("CALCMOI,1,1,1,1,-1,0")
    for k, v in ref.items():
        assert b.session["calcmoi"][k] == pytest.approx(v, abs=1e-12), k
    # oblique wall: the extruded shell carries t_eff
    ui3 = Interpreter()
    run(ui3, shell_wall_lines(0.5) + ["SLICE,1,0,0,1,0,0.6,0.8", "CSECT,2,1,0,0,1,0,0.6,0.8", "ACTM,2",
                                      "CALCMOI,0,0.6,0.8,1,0,0"])
    ref = dict(ui3.session["calcmoi"])
    assert ref["Area"] == pytest.approx(4 * 0.5 / 0.8)
    b = Interpreter()
    b.run_text(write_pre(ui3.model)[0])
    assert b.execute("CALCMOI,0,0.6,0.8,1,0,0")
    for k, v in ref.items():
        assert b.session["calcmoi"][k] == pytest.approx(v, abs=1e-12), k


def test_calcpar_on_a_unit_thickness_model_uses_its_mid_plane(ui, tmp_path):
    """A model whose DOF nodes lie on two planes one unit apart along n (here T-S1, or a CSECT model after
    WRITE -> INP) is a unit-thickness extrusion: CALCPAR cuts it at its mid-plane (spec 10 section 3.7)."""
    run(ui, two_hex_lines())
    cl.write_ess(tmp_path / "s.ess", {1: "SOLID"}, {(1, 1): [0, 0, 10, 0, 0, 1], (1, 2): [0, 0, -10, 0, 0, -1]})
    run(ui, ["READSTR,s.ess", "CALCPAR,0,0,1,1,0,0"])
    r = ui.session["calcpar"]
    assert r["Zc"] == pytest.approx(0.5) and r["My"] == pytest.approx(10.0) and r["Mz"] == pytest.approx(-1.0)
    assert any("mid-plane n.x = 0.5 of the unit-thickness extrusion" in t for t in ui.sink.texts(Kind.INFO))
    assert cl.active_section(ui.model, [0, 0, 1]).how == "extruded"


def _building_with_excavation(tmp_path):
    """2 x 1 x 2 SOLID structure (group 1, z in [0, 2]) on a 2 x 1 x 1 explicitly excavated soil layer
    (group 2, ETYPE 2, z in [-1, 0]); stresses: bottom row Szz = -10 / -30 with Sxz = 4, top row Szz = 100,
    soil Szz = 1000 (must never enter the structure's base section)."""
    lines = grid_lines(2, 1, 2) + grid_lines(2, 1, 1, z0=-1.0, group=2, first_node=101) + ["ETYPE,1,2,1,2"]
    st = {(1, 1): [0, 0, -10, 0, 4, 0], (1, 2): [0, 0, -30, 0, 4, 0], (1, 3): [0, 0, 100, 0, 0, 0],
          (1, 4): [0, 0, 100, 0, 0, 0], (2, 1): [0, 0, 1000, 0, 0, 0], (2, 2): [0, 0, 1000, 0, 0, 0]}
    cl.write_ess(tmp_path / "b.ess", {1: "SOLID", 2: "SOLID"}, st)
    return lines + ["READSTR,b.ess"]


def test_calcpar_on_a_whole_structure_gives_the_base_section(ui, tmp_path):
    """Review defect 6 / Q-11 (spec 10 section 3.8: "the global base forces and moments of a whole building"):
    on a model that is not a cross-section, CALCPAR cuts the whole structure (excavated soil left out) at its
    end on the -n side; the faces on that plane are the bottom row of elements (D-SEC-04).
    Bottom row: f1 = (4, 0, -10) at x = 0.5, f2 = (4, 0, -30) at x = 1.5, C = (1, 0.5, 0):
    F = (8, 0, -40), M = (-0.5 X) x f1 + (0.5 X) x f2 = (0, -5, 0) + (0, 15, 0) = (0, 10, 0)."""
    run(ui, _building_with_excavation(tmp_path) + ["CALCPAR,0,0,1,1,0,0"])
    r = ui.session["calcpar"]
    assert r["Area"] == pytest.approx(2.0)
    assert [r["Xc"], r["Yc"], r["Zc"]] == pytest.approx([1.0, 0.5, 0.0], abs=1e-12)
    assert [r["Fx"], r["Fy"], r["Fz"], r["Mx"], r["My"], r["Mz"]] == pytest.approx([8, 0, -40, 0, 10, 0], abs=1e-12)
    infos = ui.sink.texts(Kind.INFO)
    assert any("end of the structure on the -n side" in t and "n.x = 0" in t for t in infos)
    assert any("2 excavated-soil elements (ETYPE 2)" in t for t in infos)
    assert any("D-SEC-04" in t for t in infos)
    assert cl.active_section(ui.model, [0, 0, 1]).how == "base"
    # n = -Z: the end on the -n side is the top (n.x = -2); the top row's upward traction, seen with ez = -Z,
    # is a tension Fz = +200
    assert ui.execute("CALCPAR,0,0,-1,1,0,0")
    assert ui.session["calcpar"]["Fz"] == pytest.approx(200.0) and ui.session["calcpar"]["Zc"] == pytest.approx(2.0)
    # an oblique n only touches the base at a corner: error with guidance
    assert not ui.execute("CALCMOI,1,1,1,1,-1,0")
    assert "CSECT" in errors(ui)[-1] and "base plane" in errors(ui)[-1]


def test_csect_and_calcpar_list_the_beams_and_springs_left_out(ui, tmp_path):
    """Review defect 1 / spec 10 section 3.2: the beams CSECT draws and the springs it cannot cut are never
    section pieces; CSECT says so (information) and CALCPAR / CALCMOI repeat the warning listing them, also
    after SAVE -> RESUME; a beam deleted from the cross-section model is no longer listed."""
    run(ui, two_hex_lines() + ["N,20,5,5,0", "N,21,5,5,1", "N,22,6,5,4", "GROUP,2,BEAMS", "E,1,20,21,22",
                               "R,1,0.1,0,0,0,0,0", "N,30,7,5,0", "N,31,7,5,1", "GROUP,3,SPRING", "E,1,30,31",
                               "CUTADD,1,1,1-2", "CUTADD,1,2,1", "CUTADD,1,3,1", "CSECT,3,1,0,0,0.5,0,0,1"])
    assert any("BEAMS elements cross the plane: drawn" in t and "2/1" in t for t in ui.sink.texts(Kind.INFO))
    assert any("SPRING" in w and "3/1" in w for w in warnings(ui))
    sec = ui.models[3].ui_state["csect"]
    assert [(d["g"], d["e"]) for d in sec["beams"]] == [(2, 1)] and np.allclose(sec["beams"][0]["q"], [5, 5, 0.5])
    run(ui, ["ACTM,3", f"MDL,cs,{tmp_path}", "SAVE"])
    for interp in (ui, Interpreter(cwd=tmp_path)):
        if interp is not ui:
            run(interp, [f"MDL,cs,{tmp_path}", "RESUME"])
        n0 = len(interp.sink.texts(Kind.WARNING))
        for cmd in ("CALCPAR,0,0,1,1,0,0", "CALCMOI,0,0,1,1,0,0"):
            assert interp.execute(cmd)
            new = interp.sink.texts(Kind.WARNING)[n0:]
            assert any("BEAMS" in w and "2/1" in w for w in new), (cmd, new)
            assert any("SPRING" in w and "3/1" in w for w in new), (cmd, new)
            n0 = len(interp.sink.texts(Kind.WARNING))
        assert interp.session["calcpar"]["Area"] == pytest.approx(2.0)
    run(ui, ["GROUP,2", "EDEL,1"])
    n0 = len(warnings(ui))
    assert ui.execute("CALCPAR,0,0,1,1,0,0")
    assert not any("BEAMS" in w for w in warnings(ui)[n0:])


def test_cross_section_model_with_a_beam_after_write_inp(ui, tmp_path):
    """Review defects 2 and 4: after WRITE -> INP the extrusion is recognised from the nodes carrying DOFs
    only (a beam's K node 3.5 off the plane and a free node far away do not move the mid-plane), and CALCC
    leaves the beam out on the CSECT model and on its read-back copy, so it equals the area centroid."""
    run(ui, two_hex_lines() + ["N,20,5,5,0", "N,21,5,5,1", "N,22,6,5,4", "GROUP,2,BEAMS", "E,1,20,21,22",
                               "R,1,0.1,0,0,0,0,0", "CUTADD,1,1,1-2", "CUTADD,1,2,1", "CSECT,3,1,0,0,0.5,0,0,1",
                               "ACTM,3", "CALCMOI,0,0,1,1,0,0", "CALCC"])
    mo, cc = dict(ui.session["calcmoi"]), dict(ui.session["calcc"])
    assert mo["Area"] == pytest.approx(2.0) and cc["Volume"] == pytest.approx(2.0)
    assert [cc["Xc"], cc["Yc"], cc["Zc"]] == pytest.approx([mo["Xc"], mo["Yc"], mo["Zc"]], abs=1e-12)
    assert any("BEAMS elements are drawn but are not section pieces" in t for t in ui.sink.texts(Kind.INFO))
    run(ui, ["N,999,0,0,50"])                                    # a node of no element (e.g. a lumped mass)
    b = Interpreter()
    b.run_text(write_pre(ui.model)[0])
    assert not cl.is_csect_model(b.model)
    assert np.allclose(cl.extrusion_normal(b.model), [0, 0, 1])
    assert b.execute("CALCMOI,0,0,1,1,0,0") and b.execute("CALCC"), b.sink.texts(Kind.ERROR)
    for k, v in mo.items():
        assert b.session["calcmoi"][k] == pytest.approx(v, abs=1e-12), k
    assert [b.session["calcc"][k] for k in ("Xc", "Yc", "Zc")] == pytest.approx([mo["Xc"], mo["Yc"], mo["Zc"]],
                                                                               abs=1e-12)
    # an ordinary model with the same beam is not a cross-section: CALCC counts the beam (length 1 x 0.1)
    ui2 = Interpreter()
    run(ui2, two_hex_lines() + ["N,20,5,5,0", "N,21,5,5,3", "N,22,6,5,4", "GROUP,2,BEAMS", "E,1,20,21,22",
                                "R,1,0.1,0,0,0,0,0", "CALCC"])
    assert cl.extrusion_normal(ui2.model) is None
    assert ui2.session["calcc"]["Volume"] == pytest.approx(2.0 + 0.3)


def test_cut_system_never_redefines_the_active_system(ui):
    """Review defect 3 / spec 10 section 1.4: <sysno> equal to the active system (CSYS) is not stored, even
    without nodes (later N commands would move), for CALCMOI, CALCPAR and CALCSECTHIST alike."""
    m = run(ui, two_hex_lines() + ["LOC,5,0,0,0,0,0,0", "CSYS,5"])
    notes = cl.Notes()
    assert not cl.store_cut_system(m, 5, np.zeros(3), np.eye(3)[1], np.eye(3)[2], np.eye(3)[0], notes)
    assert "active coordinate system" in notes.warnings[0]
    assert np.allclose(m.csys[5].R, np.eye(3))
    assert ui.execute("CALCPAR,0,0,1,0,1,0,5")
    assert any("system 5 is the active coordinate system" in w for w in warnings(ui))
    assert ui.execute("CALCPAR,0,0,1,0,1,0,6") and 6 in ui.model.csys     # another number is stored
    run(ui, ["N,100,1,0,0"])
    assert np.allclose(ui.model.node_global(100), [1.0, 0.0, 0.0])


def test_stresses_loaded_under_old_group_numbers_are_reported(ui, tmp_path):
    """Review defect 5: READSTR then GCOM leaves the records under the old group number (GCOM belongs to
    generation.py); CSECT, and CALCPAR on the CSECT model or on the original model, say that the loaded
    records match no section element instead of 'no stresses are loaded'; READSTR again fixes it."""
    run(ui, grid_lines(2, 1, 1, group=4))
    cl.write_ess(tmp_path / "g.ess", {4: "SOLID"}, {(4, 1): [0, 0, 10, 0, 0, 0], (4, 2): [0, 0, 10, 0, 0, 0]})
    run(ui, ["READSTR,g.ess", "GCOM", "CALCPAR,0,0,1,1,0,0"])
    w = warnings(ui)[-1]
    assert "none matches an element of the section" in w and "stress table: 4" in w and "section: 1" in w
    run(ui, ["CUTADD,1,1,1-2", "CSECT,3,1,0,0,0.5,0,0,1"])
    assert "none matches" in warnings(ui)[-1] and "CSECT again" in warnings(ui)[-1]
    run(ui, ["ACTM,3", "CALCPAR,0,0,1,1,0,0"])
    assert "none matches" in warnings(ui)[-1] and ui.session["calcpar"]["Fz"] == 0.0
    assert not any("no element stresses are loaded" in w for w in warnings(ui))
    # CALCSECTHIST with frames written before the renumbering says why every record is missing
    (tmp_path / "g.lst").write_text("frames\ng.ess\n")
    run(ui, ["ACTM,0", "CALCSECTHIST,g.lst,1,0,0,0.5,0,0,1,1,0,0,0,0,h.csv"])
    assert any("are in frame g.ess but none matches" in w and "stress table: 4" in w for w in warnings(ui))
    cl.write_ess(tmp_path / "g1.ess", {1: "SOLID"}, {(1, 1): [0, 0, 10, 0, 0, 0], (1, 2): [0, 0, 10, 0, 0, 0]})
    run(ui, ["ACTM,0", "READSTR,g1.ess", "CSECT,3,1,0,0,0.5,0,0,1", "ACTM,3", "CALCPAR,0,0,1,1,0,0"])
    assert ui.session["calcpar"]["Fz"] == pytest.approx(20.0)


# ======================================================================================
# CALCC / CALCM
# ======================================================================================
def test_calcc_exact_volumes(ui):
    """Trilinear volume of a distorted hexahedron (a frustum-like block), a prism with repeated nodes,
    shell area x thickness and beam length x A."""
    run(ui, ["N,1,0,0,0", "N,2,2,0,0", "N,3,2,2,0", "N,4,0,2,0", "N,5,0,0,1", "N,6,1,0,1", "N,7,1,1,1", "N,8,0,1,1",
             "GROUP,1,SOLID", "E,1,1,2,3,4,5,6,7,8", "CALCC"])
    # volume of the trilinear block = int_0^1 (2 - z)^2 dz = 7/3; centroid z = int z (2-z)^2 / V = 5/12 / (7/3)
    r = ui.session["calcc"]
    assert r["Volume"] == pytest.approx(7.0 / 3.0, rel=1e-12)
    assert r["Zc"] == pytest.approx((4 / 2 - 4 / 3 + 1 / 4) / (7 / 3), rel=1e-12)
    ui2 = Interpreter()
    run(ui2, ["N,1,0,0,0", "N,2,1,0,0", "N,3,0,1,0", "N,4,0,0,2", "N,5,1,0,2", "N,6,0,1,2",
              "GROUP,1,SOLID", "E,1,1,2,3,3,4,5,6,6",
              "N,10,5,0,0", "N,11,7,0,0", "N,12,7,3,0", "N,13,5,3,0", "GROUP,2,SHELL", "E,1,10,11,12,13",
              "THICK,1,1,1,0.5", "N,20,0,0,5", "N,21,0,0,9", "N,22,1,0,5", "GROUP,3,BEAMS", "E,1,20,21,22",
              "R,1,0.25,0,0,0,0,0", "CALCC"])
    r = ui2.session["calcc"]
    assert r["Volume"] == pytest.approx(1.0 + 3.0 + 1.0, rel=1e-12)
    xc = (1.0 * (1 / 3) + 3.0 * 6.0 + 1.0 * 0.0) / 5.0
    zc = (1.0 * 1.0 + 3.0 * 0.0 + 1.0 * 7.0) / 5.0
    assert r["Xc"] == pytest.approx(xc, rel=1e-12) and r["Zc"] == pytest.approx(zc, rel=1e-12)


def test_calcm_element_and_lumped_masses(ui):
    run(ui, grid_lines(2, 1, 1) + ["M,1,3000,0.2,0.15,0.05,0.05,1", "GRAVITY,32.2",
                                   "MT,1,3.22,3.22,6.44", "MT,2,1,2,3", "MUNITS,2,2,1,0", "MR,3,10,20,30"])
    assert ui.execute("CALCM")
    r = ui.session["calcm"]
    assert r["element_weight"] == pytest.approx(0.3) and r["element_mass"] == pytest.approx(0.3 / 32.2)
    assert np.allclose(r["lumped"], [0.1 + 1, 0.1 + 2, 0.2 + 3, 10 / 32.2, 20 / 32.2, 30 / 32.2])
    assert np.allclose(r["total"], np.array(r["lumped"][:3]) + 0.3 / 32.2)
    # excavated soil is reported separately
    ui2 = Interpreter()
    run(ui2, grid_lines(1, 1, 1, z0=-1.0) + ["ETYPE,1,1,1,2", "L,1,1,0.12,1000,500,0.05,0.05", "CALCM"])
    r = ui2.session["calcm"]
    assert r["element_weight"] == 0.0 and r["soil_weight"] == pytest.approx(0.12)


def test_manual_wall_section_workflow(ui, tmp_path):
    """Manual 5.8.2 workflow pattern (spec 04 section 8.3 example 2, with the CSECT plane moved inside the
    CUTVOL slab, spec 09 errata): READSTR, CUTVOL of one wall, CSECT, ACTM, CALCPAR without sysno, CALCM."""
    zs = [0.0, 2.53, 20.0, 45.22]
    lines = []
    for k, z in enumerate(zs):
        for j, (x, y) in enumerate([(52.5, 0), (52.8, 0), (52.8, 10), (52.5, 10)]):
            lines.append(f"N,{4 * k + j + 1},{x},{y},{z}")
    lines.append("GROUP,1,SOLID")
    for k in range(3):
        b = [4 * k + j + 1 for j in range(4)]
        lines.append(f"E,{k + 1}," + ",".join(str(x) for x in b + [q + 4 for q in b]))
    lines += ["M,1,3000,0.2,0.15,0.05,0.05,1"]
    run(ui, lines)
    cl.write_ess(tmp_path / "stressfile.txt", {1: "SOLID"},
                 {(1, e): [100.0, 0, 0, 5.0, 0, 0] for e in (1, 2, 3)})
    run(ui, ["readstr,stressfile.txt", "cutvol,3,52.5,52.8,-320,320,2.53,45.22",
             "csect,1,3,52.65,0,4,1,0,0", "actm,1", "calcpar,1,0,0,0,1,0", "calcm"])
    assert sorted(ui.model.groups[1].elements) == [2, 3]
    r = ui.session["calcpar"]
    A = 10.0 * (45.22 - 2.53)
    assert r["Area"] == pytest.approx(A) and r["Xc"] == pytest.approx(52.65)
    assert r["Fz"] == pytest.approx(100.0 * A) and r["Fx"] == pytest.approx(5.0 * A)   # ex = Y: Sxy acts along Y
    # local axes: ez = X, ex = Y (r), ey = ez x ex = Z: Iyy = int x'^2 dA over the 10 m width
    assert r["Iyy"] == pytest.approx((45.22 - 2.53) * 10.0 ** 3 / 12)
    assert r["Ixx"] == pytest.approx(10.0 * (45.22 - 2.53) ** 3 / 12)
    assert any("not stored" in t for t in ui.sink.texts(Kind.INFO))
    m = ui.session["calcm"]
    assert m["element_weight"] == pytest.approx(0.15 * A * 1.0)                        # unit thickness
