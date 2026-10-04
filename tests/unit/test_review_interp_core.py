"""Adversarial review tests of work package A1 (interpreter core and model database).

These tests probe the implementation from angles other than the implementer's own tests:

* independent derivations (LOC rotation built from Rodrigues rotations about the *moving* axes,
  LOCAL invariants, generation in rotated systems, material conversion chains, section
  formulas in limiting cases, lazy weight -> mass conversion);
* round trips (WRITE -> INP -> WRITE on randomised models, SAVE/RESUME exactness);
* substitution order and loop/macro semantics (L13, D-PAR-15/16, spec 04 sections 5-6);
* list / request / "last defined" semantics (UT-04, D-MDL-05, D-MDL-08, D-MDL-11).

The ``test_defect_*`` tests document defects found by the review.  They were marked
``xfail(strict=True)`` until the defects were fixed; the markers have been removed since.
"""
from __future__ import annotations

import json
import math
import random

import numpy as np
import pytest

from sassi.model.materials import elastic_constants, section_rectangle, section_circle
from sassi.prep import Interpreter, Kind
from sassi.prep.writer import write_pre


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def run(lines, **kw):
    ui = Interpreter(seed=1, **kw)
    ui.run_text("\n".join(lines))
    return ui


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


def roundtrip(ui):
    """WRITE the active model, INP it into a fresh interpreter, return the new interpreter."""
    text, _ = write_pre(ui.model, filename="rt.pre")
    ui2 = Interpreter(seed=1)
    ui2.run_text(text)
    return ui2, text


def rodrigues(axis, deg):
    """Right-hand rotation matrix about a unit ``axis`` by ``deg`` degrees (independent of LOC code)."""
    k = np.asarray(axis, float)
    k = k / np.linalg.norm(k)
    a = math.radians(deg)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) * math.cos(a) + math.sin(a) * K + (1 - math.cos(a)) * np.outer(k, k)


# ======================================================================================
# Coordinate systems (D-MDL-04, spec 08 sections 2-3, UT-12)
# ======================================================================================
@pytest.mark.parametrize("angles", [(45, 0, 0), (0, 90, 0), (30, 20, 10), (-75, 33, 140), (10, 89.9, -45)])
def test_loc_equals_successive_rotations_about_moving_axes(angles):
    """LOC: txy about Z, then tyz about the *new* X, then txz about the *newest* Y (spec 08 3.8).

    Built here from Rodrigues rotations about the current axis vectors, without Rz Rx Ry.
    """
    txy, tyz, txz = angles
    A = np.eye(3)                                   # columns: current local axes in global
    A = rodrigues(A[:, 2], txy) @ A
    A = rodrigues(A[:, 0], tyz) @ A
    A = rodrigues(A[:, 1], txz) @ A
    ui = run([f"LOC,7,0,1.5,-2,3,{txy},{tyz},{txz}", "CSYS,7", "N,1,0.3,-0.7,2.1"])
    cs = ui.model.csys[7]
    assert np.allclose(cs.R, A, atol=1e-13)
    assert np.isclose(np.linalg.det(cs.R), 1.0)
    X = ui.model.node_global(1)
    assert np.allclose(X, np.array([1.5, -2, 3]) + 0.3 * A[:, 0] - 0.7 * A[:, 1] + 2.1 * A[:, 2], atol=1e-13)


def test_ut12_reference_values_and_global_keeps_position():
    ui = run(["LOC,1,0,10,0,0,45,0,0", "CSYS,1", "N,100,5,0,0"])
    before = ui.model.node_global(100)
    ui.execute("GLOBAL,100,100,1")
    n = ui.model.nodes[100]
    assert n.csys == 0 and ui.model.csys_active == 1
    assert np.allclose(n.xyz, (10 + 5 / math.sqrt(2), 5 / math.sqrt(2), 0.0), atol=1e-12)
    assert np.allclose(n.xyz, before, atol=0)


@pytest.mark.parametrize("seed", range(6))
def test_local_invariants_random_points(seed):
    """LOCAL: origin n1, n2 on +x, n3 in the x-y plane with y > 0, right-handed (spec 08 3.9, Q7)."""
    rng = np.random.default_rng(seed)
    p = rng.uniform(-10, 10, size=(3, 3))
    lines = [f"N,{i + 1},{float(p[i, 0])!r},{float(p[i, 1])!r},{float(p[i, 2])!r}" for i in range(3)] + ["LOCAL,4,0,1,2,3"]
    ui = run(lines)
    m = ui.model
    cs = m.csys[4]
    assert np.allclose(cs.R.T @ cs.R, np.eye(3), atol=1e-12)
    assert np.isclose(np.linalg.det(cs.R), 1.0)
    l1, l2, l3 = (m.to_system(m.node_global(i), 4) for i in (1, 2, 3))
    assert np.allclose(l1, 0, atol=1e-11)
    assert np.isclose(l2[0], np.linalg.norm(p[1] - p[0])) and np.allclose(l2[1:], 0, atol=1e-11)
    assert l3[1] > 0 and abs(l3[2]) < 1e-11


def test_ngen_increments_follow_active_axes():
    """NGEN increments are applied along the active system axes (spec 08 Q3)."""
    ui = run(["LOC,1,0,5,5,0,90,0,0", "N,1,1,0,0", "CSYS,1", "NGEN,2,10,1,1,1,1,0,0"])
    m = ui.model
    R = m.csys[1].R
    X1 = m.node_global(1)
    for k, nid in ((1, 11), (2, 21)):
        assert np.allclose(m.node_global(nid), X1 + k * R @ np.array([1.0, 0, 0]), atol=1e-12)
        assert m.nodes[nid].csys == 1


def test_fill_and_nmed_are_frame_independent():
    """FILL and NMED give the same geometry in any Cartesian active system (spec 08 2.2)."""
    pts = {1: (0.0, 0.0, 0.0), 5: (40.0, 8.0, -4.0), 7: (3.0, 1.0, 2.0)}
    base = [f"N,{k},{v[0]},{v[1]},{v[2]}" for k, v in pts.items()]
    for csys in (["CSYS,0"], ["LOC,2,0,3,-1,7,30,40,50", "CSYS,2"]):
        ui = run(base + csys + ["FILL,1,5,3", "NMED,50,1,5,7"])
        m = ui.model
        for k, nid in enumerate((2, 3, 4), start=1):
            assert np.allclose(m.node_global(nid), np.array(pts[1]) + k / 4 * (np.array(pts[5]) - pts[1]), atol=1e-12)
        assert np.allclose(m.node_global(50), np.mean([pts[1], pts[5], pts[7]], axis=0), atol=1e-12)


def test_sdel_keeps_global_positions_and_resets_active():
    ui = run(["LOC,1,0,1,2,3,10,20,30", "CSYS,1", "N,1,1,1,1", "N,2,-2,0,4"])
    g = {i: ui.model.node_global(i) for i in (1, 2)}
    ui.execute("SDEL,1")
    m = ui.model
    assert m.csys_active == 0 and 1 not in m.csys
    for i in (1, 2):
        assert m.nodes[i].csys == 0 and np.allclose(m.nodes[i].xyz, g[i], atol=0)


def test_redefining_loc_moves_stored_nodes_with_warning():
    """D-MDL-04: a redefined system keeps the local numbers of its nodes (they move), with a warning."""
    ui = run(["LOC,1,0,0,0,0,0,0,0", "CSYS,1", "N,1,1,0,0", "LOC,1,0,10,0,0,0,0,0"])
    assert np.allclose(ui.model.node_global(1), (11, 0, 0))
    assert any("redefined" in w for w in warnings(ui))


# ======================================================================================
# Node and element generation (UT-13, spec 08 sections 3 and 5)
# ======================================================================================
def test_fill_with_explicit_count_uses_integer_increment():
    ui = run(["N,1,0,0,0", "N,9,8,0,0", "FILL,1,9,3"])
    assert sorted(ui.model.nodes) == [1, 3, 5, 7, 9]
    assert [ui.model.nodes[i].x for i in (3, 5, 7)] == [2.0, 4.0, 6.0]


def test_fill_refuses_existing_target():
    ui = run(["N,1", "N,3", "N,5,4", "FILL,1,5"])
    assert any("already exist" in e for e in errors(ui))
    assert sorted(ui.model.nodes) == [1, 3, 5]


def test_nmove_explicit_zero_kept_and_blank_is_one():
    ui = run(["N,1,2,3,4", "NMOVE,,0,2,10,1"])
    assert ui.model.nodes[10].xyz == (2.0, 0.0, 8.0)
    assert any("collapses" in w for w in warnings(ui))


def test_lmove_translates_list_in_order():
    ui = run(["N,3,1,0,0", "N,1,0,0,0", "LMOVE,0,0,10,101,3,1"])
    assert ui.model.nodes[101].xyz == (1.0, 0.0, 10.0)
    assert ui.model.nodes[102].xyz == (0.0, 0.0, 10.0)


def test_egen_increments_k_node_and_copies_attributes():
    """D-MDL-09: EGEN increments every node number (K included) and copies the pattern attributes."""
    ui = run(["GROUP,1,BEAMS", "MACT,3", "RACT,4", "E,1,1,2,99", "KI,1,1,1,0,0,0,0,1,1", "ETYPE,1,1,1,1",
              "MACT,1", "RACT,1", "EGEN,2,10,1"])
    g = ui.model.groups[1]
    assert sorted(g.elements) == [1, 2, 3]
    for k, eid in ((1, 2), (2, 3)):
        e = g.elements[eid]
        assert e.nodes == [1 + 10 * k, 2 + 10 * k, 99 + 10 * k]
        assert (e.mat, e.prop, e.etype, e.ki) == (3, 4, 1, [0, 0, 0, 0, 1, 1])


def test_ecompr_keeps_order_and_remaps_eout():
    ui = run(["GROUP,2,SOLID", "E,1,1,2,3,4,5,6,7,8", "E,3,1,2,3,4,5,6,7,8", "E,7,1,2,3,4,5,6,7,8",
              "EOUT,1,0,0,0,0,0,0,0,0,0,0,0,2,3 7", "ECOMPR"])
    assert sorted(ui.model.groups[2].elements) == [1, 2, 3]
    assert ui.model.eout[0].elements == [2, 3]


def test_ndel_cascade_and_dangling_count():
    """D-MDL-11: loads, masses and units of deleted nodes go; element references stay (Error 41)."""
    ui = run(["N,1", "N,2", "N,3", "MT,2,1,1,1", "MR,2,1,1,1", "MUNITS,2,2,1,0", "F,2,1,0,0", "MM,2,0,1,0",
              "GROUP,1,SPRING", "E,1,1,2", "E,2,2,3", "NDEL,2"])
    m = ui.model
    assert 2 not in m.nodes and 2 not in m.tmass and 2 not in m.rmass and 2 not in m.forces
    assert 2 not in m.moments and 2 not in m.mass_units
    assert m.groups[1].elements[2].nodes == [2, 3]
    assert any("2 element node references" in w for w in warnings(ui))


def test_d_expansions_and_default_free():
    ui = run(["N,1", "N,2", "N,3", "D,1,3,1,1,ALL", "D,2,2,,,ROT", "D,3,3,1,0,DISP"])
    fx = {i: ui.model.nodes[i].fix for i in (1, 2, 3)}
    assert fx == {1: [1] * 6, 2: [1, 1, 1, 0, 0, 0], 3: [0, 0, 0, 1, 1, 1]}


# ======================================================================================
# Materials, sections, masses (UT-14, UT-15, spec 08 section 6)
# ======================================================================================
@pytest.mark.parametrize("E,nu", [(30000.0, 0.25), (519120.0, 0.17), (1.0e7, 0.45), (2.0e5, 0.01)])
def test_material_type_chain_1_2_3(E, nu):
    """Type 1 -> (M, G) as type 2 -> (Vp, Vs) as type 3 returns E and nu (independent Lame formulas)."""
    gamma, g = 0.15, 32.2
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    G = E / (2 * (1 + nu))
    M = lam + 2 * G
    c1 = elastic_constants(1, E, nu, gamma, g)
    assert np.isclose(c1.G, G, rtol=1e-13) and np.isclose(c1.M, M, rtol=1e-13) and np.isclose(c1.lam, lam, rtol=1e-12)
    c2 = elastic_constants(2, c1.M, c1.G, gamma, g)
    c3 = elastic_constants(3, c1.Vp, c1.Vs, gamma, g)
    for c in (c2, c3):
        assert np.isclose(c.E, E, rtol=1e-12) and np.isclose(c.nu, nu, rtol=1e-11)
    rho = gamma / g
    assert np.isclose(c1.Vs, math.sqrt(G / rho), rtol=1e-13) and np.isclose(c1.Vp, math.sqrt(M / rho), rtol=1e-13)


def test_material_command_uses_gravity_lazily():
    """Spec 08 1.6: raw values are stored; GRAVITY given after M changes the derived density."""
    ui = run(["M,1,4000,2000,0.12,0.05,0.05,3", "GRAVITY,9.81"])
    k = ui.model.materials[1].constants(ui.model.gravity)
    assert np.isclose(k.rho, 0.12 / 9.81) and np.isclose(k.G, 0.12 / 9.81 * 2000 ** 2)


def test_section_formulas_limits():
    # J symmetric in (b, h) by the swap rule; square coefficient of the Roark formula
    assert np.isclose(section_rectangle(0.5, 1.0)["tors"], section_rectangle(1.0, 0.5)["tors"])
    assert np.isclose(section_rectangle(1.0, 1.0)["tors"], 1 / 3 - 0.21 * (1 - 1 / 12))
    # thin strip: J -> h b^3 / 3
    b, h = 1e-3, 1.0
    assert np.isclose(section_rectangle(b, h)["tors"], h * b ** 3 / 3, rtol=1e-3)
    # I2 about axis 2 (b along 3), I3 about axis 3
    s = section_rectangle(0.5, 1.0)
    assert np.isclose(s["flex2"], 1.0 * 0.5 ** 3 / 12) and np.isclose(s["flex3"], 0.5 * 1.0 ** 3 / 12)
    c = section_circle(2.0)
    assert np.isclose(c["tors"], c["flex2"] + c["flex3"])        # polar = I2 + I3 for a circle
    assert np.isclose(c["shear2"], 0.9 * c["axial"])


def test_mass_units_lazy_gravity_and_rotational():
    """UT-15 / D-MDL-08: default MUNITS 1 (weight / g, also for MR), MUNITS 0 = mass; g read lazily."""
    ui = run(["MT,1,32.2,64.4,0", "MR,1,3.22,0,0", "MT,2,5,5,5", "MUNITS,2,2,1,0"])
    m = ui.model
    assert np.allclose(m.nodal_masses()[1], [1, 2, 0, 0.1, 0, 0])
    assert np.allclose(m.nodal_masses()[2], [5, 5, 5, 0, 0, 0])
    ui.execute("GRAVITY,9.81")
    assert np.allclose(m.nodal_masses()[1][:2], [32.2 / 9.81, 64.4 / 9.81])


def test_general_matrix_symmetric_fill_and_weight_units():
    ui = run(["MXR,1,1,1,2,3,4,5,6,7,8,9,10,11,12", "MXR,1,12,99", "MXM,1,3,32.2", "MOPT,,1"])
    m = ui.model
    K, Mm = m.general_matrices(1)
    assert np.allclose(K.real, K.real.T)
    assert K[0, 11] == 12 and K[11, 0] == 12 and K[11, 11] == 99
    assert np.isclose(Mm[2, 2], 1.0)          # weight 32.2 / g 32.2


# ======================================================================================
# Loads and "last defined" (D-MDL-05, D-MDL-08, spec 08 section 7)
# ======================================================================================
def test_fscale_default_range_uses_definition_order():
    ui = run(["F,10,1,0,0", "F,30,1,0,0", "F,20,1,0,0", "F,25,1,0,0", "F,22,1,0,0", "FSCALE,,,,3"])
    f = {k: v.factor[0] for k, v in ui.model.forces.items()}
    # last two defined: 25 and 22 -> range 22..25 -> 22 and 25 scaled
    assert f == {10: 1.0, 30: 1.0, 20: 1.0, 25: 3.0, 22: 3.0}


def test_mopt_add_mode_for_masses_and_forces():
    ui = run(["MOPT,1,0,0,0", "MT,1,1,2,3", "MT,1,1,1,1", "F,5,1,0,0,0.1", "F,5,2,0,0,0.3"])
    m = ui.model
    assert m.tmass[1] == [2.0, 3.0, 4.0]
    assert m.forces[5].factor[0] == 3.0 and m.forces[5].arrival[0] == 0.3
    ui.execute("MOPT")
    ui.execute("MT,1,7,7,7")
    assert m.tmass[1] == [7.0, 7.0, 7.0]


def test_mtgen_copies_units_and_increments():
    ui = run(["MT,100,1,2,3", "MUNITS,100,100,1,0", "MTGEN,3,10,100,100,1,1,0,0"])
    m = ui.model
    assert [m.tmass[n][0] for n in (110, 120, 130)] == [2.0, 3.0, 4.0]
    assert all(m.mass_unit(n) == 0 for n in (110, 120, 130))


# ======================================================================================
# Substitution, loops, macros (L13, D-PAR-15/16, spec 04 sections 5-6, UT-09/10)
# ======================================================================================
def test_substitution_left_to_right_once_each():
    ui = run(["VAR,K", "N,@K++,@K,@K++"])
    assert ui.model.nodes[1].xyz == (1.0, 2.0, 0.0)
    assert ui.variables["k"].counter == 2


def test_loop_reversed_nesting_numbering():
    """ForEach,X,ForEach,Y,ForEach,Z: X outermost, so node 25(ix-1)+5(iy-1)+iz at (ix, iy, iz)."""
    ui = run(["Var,NNUM", "Var,X,1,2,3,4,5", "Var,Y,1,2,3,4,5", "Var,Z,1,2,3,4,5",
              "ForEach,X,ForEach,Y,ForEach,Z,N,@NNUM++,@X[#],@Y[#],@Z[#]"])
    m = ui.model
    assert len(m.nodes) == 125
    for ix in (1, 3, 5):
        for iy in (1, 4):
            for iz in (2, 5):
                assert m.nodes[25 * (ix - 1) + 5 * (iy - 1) + iz].xyz == (ix, iy, iz)


def test_hash_binds_to_its_own_variable_and_bare_hash_is_innermost():
    ui = run(["VAR,A,10,20", "VAR,B,1,2,3", "VAR,C",
              "FOREACH,A,FOREACH,B,N,@C++,@A[#],@B[#],#"])
    m = ui.model
    assert len(m.nodes) == 6
    assert m.nodes[1].xyz == (10, 1, 1) and m.nodes[3].xyz == (10, 3, 3) and m.nodes[6].xyz == (20, 3, 3)


def test_index_out_of_range_skips_command():
    ui = run(["VAR,L,1,2", "N,5,@L[3]"])
    assert 5 not in ui.model.nodes and errors(ui)


def test_nested_macros_with_substring_substitution(tmp_path):
    (tmp_path / "inner.pre").write_text("N,$1$,$2$,0,0\n")
    (tmp_path / "outer.pre").write_text("MACRO,INNER,$1$0,$1$\n")
    ui = Interpreter(cwd=tmp_path)
    ui.run_text("LOADMACRO,inner,inner.pre\nLOADMACRO,Outer,outer.pre\nMACRO,outer,7\nMACRO,OUTER,3")
    assert ui.model.nodes[70].x == 7.0 and ui.model.nodes[30].x == 3.0
    assert not errors(ui)


def test_macro_name_does_not_override_command(tmp_path):
    (tmp_path / "n.pre").write_text("N,99\n")
    ui = Interpreter(cwd=tmp_path)
    ui.run_text("LOADMACRO,N,n.pre\nN,1,1,1,1")
    assert sorted(ui.model.nodes) == [1]


def test_inp_summary_and_single_eof_message(tmp_path):
    (tmp_path / "b.pre").write_text("N,2\nNOPE,1\n")
    (tmp_path / "a.pre").write_text("* c\nN,1\nINP,b.pre\nDAMP,0.05,0.02\n")
    ui = Interpreter(cwd=tmp_path)
    ui.execute("INP,a.pre")
    info = ui.sink.texts(Kind.INFO)
    assert info.count("INPUT FILE REACHED EOF, INPUT SWITCHED TO KEYBOARD") == 1
    assert "INP b.pre: 2 lines, 2 commands, 0 warnings, 1 errors" in info
    assert "INP a.pre: 4 lines, 5 commands, 0 warnings, 1 errors" in info
    assert not any(m.kind == Kind.ECHO for m in list(ui.sink.messages)[1:])   # only the INP line echoed


# ======================================================================================
# Lists and frequency sets (UT-04, UT-16 subset)
# ======================================================================================
def test_freq_sets_and_hz_from_site_shared_variables():
    ui = run(["SITE,0,1,,,,,,,,,,0.005,4096,1", "FREQ,1,512,8,1", "FREQ,2,3", "FREQ,2,0"])
    m = ui.model
    assert np.isclose(m.frequency_step(), 0.048828125, rtol=0, atol=1e-15)
    nums, hz = m.frequencies(1)
    assert nums == [1, 8, 512] and np.isclose(hz[-1], 25.0)
    assert 2 not in m.freq_sets
    ui.execute("SITE,0,1,0.9765625,,,,,,,,,0.005,4096,1")
    assert m.frequency_step() == 0.9765625


def test_requests_append_and_clear_written_first():
    ui = run(["NOUT,1,1,0,0,0,0,1,1-3 7;9", "NOUT,3,0,0,0,0,0,1,4", "RDND,15,1,0,1,0,1,0"])
    m = ui.model
    assert [r.nodes for r in m.nout] == [[1, 2, 3, 7, 9], [4]]
    text, _ = write_pre(m)
    lines = text.splitlines()
    assert lines.index("NOUT,0") < next(i for i, l in enumerate(lines) if l.startswith("NOUT,1,"))
    ui.execute("NOUT,0")
    assert m.nout == []


# ======================================================================================
# WRITE -> INP (UT-03, D-PAR-17) and SAVE/RESUME (D-MDL-02)
# ======================================================================================
def _random_model(seed):
    rng = random.Random(seed)
    nn = rng.randint(5, 25)
    L = [f"N,{i},{rng.uniform(-9, 9)!r},{rng.uniform(-9, 9)!r},{rng.uniform(-9, 9)!r}" for i in range(1, nn + 1)]
    for _ in range(rng.randint(15, 50)):
        a = rng.randint(1, nn)
        b = rng.randint(a, nn)
        L.append(rng.choice([
            f"LOC,{rng.randint(1, 3)},0,{rng.uniform(-3, 3)!r},0,1,{rng.uniform(0, 90)!r},{rng.uniform(0, 90)!r},7",
            f"CSYS,{rng.randint(0, 3)}",
            f"N,{rng.randint(1, nn + 5)},{rng.uniform(-9, 9)!r},0.5,{rng.uniform(-9, 9)!r}",
            f"D,{a},{b},1,1,{rng.choice(['UX', 'ROTZ', 'DISP', 'ROT', 'ALL'])}",
            f"INT,{a},{b},1,1,{rng.randint(0, 3)}",
            f"M,{rng.randint(1, 3)},{rng.uniform(1e3, 1e6)!r},0.3,0.15,0.04,0.04,1",
            f"GROUP,{rng.randint(1, 3)},{rng.choice(['SOLID', 'BEAMS', 'SHELL', 'SPRING'])}",
            "E," + ",".join(str(rng.randint(1, nn)) for _ in range(rng.randint(2, 4))),
            f"MACT,{rng.randint(1, 3)}",
            f"MT,{a},{rng.uniform(0, 5)!r},1,2",
            f"MUNITS,{a},{b},1,{rng.randint(0, 1)}",
            f"F,{a},{rng.uniform(-1, 1)!r},0,1,{rng.uniform(0, 1)!r}",
            f"FREQ,1,{rng.randint(1, 50)},{rng.randint(1, 50)}",
            f"DAMP,{rng.choice([0, 0.02, 0.05])}",
            f"TOPL,{rng.randint(0, 3)},2",
            f"THICK,1,3,1,{rng.uniform(0.1, 2)!r}",
            f"MXR,1,{rng.randint(1, 12)},{rng.uniform(-1, 1)!r},0,3",
            f"NDEL,{a}",
            f"SITE,0,1,,,,,,,,,,0.005,4096,{rng.randint(1, 2)}",
            f"DYNP,{rng.randint(1, 3)},0.001,0.9,0.001,2.5,Sand, dense",
            "SOIL,1,32.2,(5X,F10.4)",
        ]))
    return L


@pytest.mark.parametrize("seed", range(25))
def test_write_inp_write_idempotent_random_models(seed):
    ui = run(_random_model(seed))
    ui2, text = roundtrip(ui)
    assert ui.model.same_state(ui2.model)
    assert write_pre(ui2.model, filename="rt.pre")[0] == text


def test_numbers_round_trip_bit_exact():
    vals = [0.1 + 0.2, 1e-300, -1.0e300, 2.0 ** -1074, 123456789012345.6, -0.0, 1 / 3]
    ui = run([f"N,{i + 1},{v!r},0,0" for i, v in enumerate(vals)])
    ui2, _ = roundtrip(ui)
    for i, v in enumerate(vals):
        x = ui2.model.nodes[i + 1].x
        assert x == v and math.copysign(1, x) == math.copysign(1, v)


def test_save_resume_restores_full_state_including_histories(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(f"MDL,demo,{tmp_path}\nN,5\nN,1\nN,3\nF,3,1,0,0\nF,1,1,0,0\nLOC,2,0,1,1,1,5,6,7\nSAVE")
    saved = json.dumps(ui.model.to_json(), sort_keys=True)
    ui.run_text("N,77\nNDEL,5\nRESUME")
    assert json.dumps(ui.model.to_json(), sort_keys=True) == saved
    ui.execute("FILL")          # last two defined nodes are 1 and 3 (definition order restored)
    assert 2 in ui.model.nodes


# ======================================================================================
# Defects found by the review.  They were strict xfails and are fixed now (the markers were
# removed as the module docstring prescribes); each test states the requirement it checks.
# ======================================================================================
def test_defect_gdel_active_group_round_trip():
    """UT-03: 'no active group' after GDEL of the active group (spec 08 5.12) survives WRITE -> INP."""
    ui = run(["GROUP,1,SOLID", "GROUP,2,SHELL", "GDEL,2"])
    assert ui.model.group_active is None
    ui2, _ = roundtrip(ui)
    assert ui.model.same_state(ui2.model)


def test_defect_mtype_kept_attribute_round_trip():
    """UT-03 / spec 08 5.24: data kept after MTYPE (THICK of a former SHELL group) is written and read back."""
    ui = run(["GROUP,1,SHELL", "E,1,1,2,3,4", "THICK,1,1,1,0.5", "MTYPE,1,PLANE"])
    ui2, _ = roundtrip(ui)
    assert ui.model.same_state(ui2.model)


def test_defect_mtype_elements_lost_on_round_trip():
    """UT-03: after MTYPE BEAMS -> SPRING the 3-node elements survive WRITE -> INP."""
    ui = run(["GROUP,1,BEAMS", "E,1,1,2,3", "MTYPE,1,SPRING"])
    ui2, _ = roundtrip(ui)
    assert sorted(ui2.model.groups[1].elements) == [1]


def test_defect_local_with_deleted_helper_nodes_round_trip_exact():
    """UT-03: a LOCAL system whose helper nodes were deleted is re-created bit for bit."""
    ui = run(["N,901,0.3,0.1,0", "N,902,1.7,0.9,0.2", "N,903,-0.4,1.3,0.5", "LOCAL,1,0,901,902,903",
              "NDEL,901,903", "CSYS,1", "N,1,1.1,2.2,3.3"])
    ui2, _ = roundtrip(ui)
    assert np.array_equal(ui2.model.csys[1].R, ui.model.csys[1].R)
    assert np.array_equal(ui2.model.node_global(1), ui.model.node_global(1))


def test_defect_inp_utf8_bom(tmp_path):
    """Spec 04 4.1 (old .pre files load unchanged): a UTF-8 byte-order mark is ignored."""
    (tmp_path / "bom.pre").write_bytes("N,1,1,2,3\r\nN,2\r\n".encode("utf-8-sig"))
    ui = Interpreter(cwd=tmp_path)
    ui.execute("INP,bom.pre")
    assert sorted(ui.model.nodes) == [1, 2] and not errors(ui)


def test_defect_nmed_ninth_node_used_despite_warning():
    """NMED takes 1-8 nodes (spec 08 3.14): a 9th node reported as ignored is really ignored."""
    lines = [f"N,{i},{4 * (i - 1)},0,0" for i in range(1, 9)] + ["N,9,1000,0,0", "NMED,500,1,2,3,4,5,6,7,8,9"]
    ui = run(lines)
    if 500 in ui.model.nodes:          # accepted: then the 9th node must really be ignored
        assert ui.model.nodes[500].x == pytest.approx(14.0)


def test_defect_lfreq_with_unresolvable_df():
    """LFREQ lists the frequency numbers when the frequency step cannot be resolved (no internal error)."""
    ui = run(["SITE,0,1,0,,,,,,,,,0,4096,1", "FREQ,1,5", "LFREQ"])
    assert not any("internal error" in e for e in errors(ui))


def test_defect_status_mtlist_robustness():
    """STATUS / MTLIST report unavailable derived values (SITE delt = 0, gravity 0) instead of failing."""
    ui = run(["SITE,0,1,0,,,,,,,,,0,4096,1", "FREQ,1,5", "STATUS", "GRAVITY,0", "MT,1,1,1,1", "MTLIST"])
    assert not any("internal error" in e for e in errors(ui))


def test_defect_generators_create_nonpositive_ids():
    """Node and element numbers must be positive (as N/E check): NGEN, LMOVE and EGEN refuse ids <= 0."""
    ui = run(["N,1", "N,2,1", "NGEN,1,-5,1,2", "LMOVE,0,0,1,-3,1,2",
              "GROUP,1,BEAMS", "E,1,1,2,3", "EGEN,1,1,1,1,1,-5"])
    assert min(ui.model.nodes) > 0
    assert min(ui.model.groups[1].elements) > 0


def test_defect_inp_missing_file_reports_failure(tmp_path):
    """INP of a missing file fails: execute() returns False (sassi -c exit status, replay history)."""
    ui = Interpreter(cwd=tmp_path)
    assert ui.execute("INP,missing.pre") is False


def test_defect_text_argument_with_unbalanced_quote():
    """L7 / D-PAR-06: a rest-of-line text argument is taken verbatim, an unbalanced quote included."""
    ui = run(['GROUP,1,SOLID', 'GTIT,1,Walls, "north', 'THTIT,acc X, "raw'])
    assert ui.model.groups[1].title == 'Walls, "north'
    assert ui.model.options.string("THTIT") == 'acc X, "raw'
