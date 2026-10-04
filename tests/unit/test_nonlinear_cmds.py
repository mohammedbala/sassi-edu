"""Unit tests of the Option NON commands (sassi/prep/commands/nonlinear_cmds.py; manual 9.17, requirements
3.4.O, 3.4.R, D-NON-01 ... D-NON-12): records and WRITE -> INP, backbone-curve commands, deletions,
panel-model construction (WALLFLR, PANELIZE, EDGE, UNIPNL, DGRDFLR, MERGEPANEL, PNLGEN, SOLIDPILE),
NONLINMOTDISP, SHEAR / BBCGEN, the .eql deck (CHECK errors 121-128, EDU-09), NONLINBAT, COMBXYZTHD and
the run-control extension commands."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.prep import Interpreter
from sassi.prep.commands import nonlinear_cmds as NC
from sassi.prep.messages import Kind
from sassi.prep.writer import write_pre


def run(ui, text):
    """Run command text with INP semantics; the message sink only holds this run's messages."""
    ui.sink.clear()
    ui.run_text(text)
    return ui


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


def roundtrip(ui):
    text, _ = write_pre(ui.model)
    b = Interpreter()
    b.run_text(text)
    assert not errors(b), errors(b)[:5]
    assert b.model.same_state(ui.model)
    return b


WALL = """
* a 12 m x 4 m wall in the XZ plane (4 x 2 shells) and a spring
N,1,0,0,0
N,5,12,0,0
FILL,1,5
NGEN,2,5,1,5,1,0,0,2
M,1,3.0E7,0.2,24.0,0.04,0.04,1
M,2,3.0E7,0.2,24.0,0.04,0.04,1
GROUP,1,SHELL
MACT,1
E,1,1,2,7,6
EGEN,3,1,1
EGEN,1,5,1,4
THICK,1,8,1,0.3
N,20,0,0,-1
SC,1,1000,0,0,0,0,0,0.02
GROUP,2,SPRING
RACT,1
E,1,20,1
GRAVITY,9.81
"""


@pytest.fixture
def wall(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(WALL)
    assert not errors(ui), errors(ui)
    return ui


# ======================================================================================
# Records and WRITE -> INP
# ======================================================================================
def test_records_round_trip(wall):
    run(wall, """EQL,0.8,3,7,0,1
P,1,1,1,1,1
S,1,2,1,2,1,4
B,1,3,4,2,4,1,2
BBCX,1,3,2,1E-4,2E-4,3E-3
BBCY,1,3,2,4500,5000,6000
BBCI,1,2,1
BBCP,2,1,0.01,10
BBCP,2,2,0.2,29
BBCI,2,1,4
""")
    assert not errors(wall), errors(wall)
    m = wall.model
    assert m.options.record("EQL").get("disp") == 0.8 and m.options.record("EQL").get("dampcutoff") == 7.0
    assert m.options.entry("P", 1).get("group") == 1 and m.options.entry("S", 1).get("disp") == 1
    cv = NC.curve(m, 1)
    assert cv["x"] == [1e-4, 2e-4, 3e-3] and cv["y"] == [4500, 5000, 6000] and cv["yield_"] == 2 and cv["type"] == 1
    assert NC.curve(m, 2) == dict(type=4, yield_=1, x=[0.01, 0.2], y=[10.0, 29.0])
    b = roundtrip(wall)
    assert NC.curve(b.model, 2)["x"] == [0.01, 0.2]


def test_eql_validation(wall):
    run(wall, "EQL,1.5")
    assert errors(wall)
    run(wall, "EQL,0.8,8")
    assert any("bit mask" in e for e in errors(wall))
    run(wall, "EQL,0.8,1,0.07")
    assert any("percent" in w for w in warnings(wall))


def test_panel_and_spring_messages(wall):
    run(wall, "P,1,1,1,1,3")
    assert any("126" in w for w in warnings(wall))
    run(wall, "P,2,2,1,1,1")
    assert any("SHELL" in w for w in warnings(wall))
    run(wall, "S,1,2,1,1,1,1")
    assert any("127" in w for w in warnings(wall))
    run(wall, "S,2,2,1,1,4,4")
    assert errors(wall)
    run(wall, "B,1,1,2,1,4,1,2")
    assert any("not usable" in w for w in warnings(wall))


def test_bbc_from_file(wall, tmp_path):
    p = tmp_path / "bbc.txt"
    p.write_text("* shear strain, force\n1.0E-4 4500\n2.0E-4 5000\n3.0E-3 6000\n0.02 6100\n", encoding="utf-8")
    run(wall, f"BBC,3,1,3,2,{p}")
    assert not errors(wall) and any("ignored" in w for w in warnings(wall))
    assert NC.curve(wall.model, 3)["x"] == [1e-4, 2e-4, 3e-3]
    run(wall, f"BBC,4,1,5,2,{p}")
    assert errors(wall) and NC.curve(wall.model, 4) is None


def test_bbcx_bbcy_consistency_and_bbcp_append(wall):
    run(wall, "BBCX,5,3,2,1,2,3")
    run(wall, "BBCY,5,2,2,10,20")
    assert any("both vectors" in w for w in warnings(wall))
    assert any("Y values" in p for p in NC.curve_problems(wall.model, 5))
    run(wall, "BBCP,5,3,3,30")
    assert NC.curve(wall.model, 5)["y"] == [10.0, 20.0, 30.0]
    run(wall, "BBCP,5,5,9,40")                        # skips point 4: an empty grid row (0, 0) until defined
    assert not errors(wall) and any("points 4..4 are set to (0, 0)" in w for w in warnings(wall))
    assert NC.curve(wall.model, 5)["x"] == [1.0, 2.0, 3.0, 0.0, 9.0]
    run(wall, "BBCI,5,3,4")
    assert any("increas" in p for p in NC.curve_problems(wall.model, 5))
    run(wall, "BBCP,5,4,6,35")
    cv = NC.curve(wall.model, 5)
    assert cv["yield_"] == 3 and wall.model.options.entry("BBCX", 5).integer(3) == 3
    assert cv["y"] == [10.0, 20.0, 30.0, 35.0, 40.0] and not NC.curve_problems(wall.model, 5)


def test_range_deletions(wall):
    run(wall, "P,1,1,1\nP,2,1,1\nP,3,1,1\nPDEL,1,3,2")
    assert sorted(k for k, _ in NC.panels(wall.model)) == [2]
    run(wall, "S,1,2,1,1,1,4\nS,2,2,1,1,1,4\nDELSPR,1,2,0")
    assert any("stride" in w for w in warnings(wall)) and not NC.springs(wall.model)
    run(wall, "BBCX,1,2,1,1,2\nBBCY,1,2,1,1,2\nBBCX,2,2,1,1,2\nDELBBC,1,2")
    assert NC.curve_numbers(wall.model) == []
    run(wall, "B,1,1,1,1\nDELBM,1")
    assert not wall.model.options.entries("B")


# ======================================================================================
# Panels: PNLGEN, PLIST, NONLINMOTDISP, SHEAR, BBCGEN
# ======================================================================================
def test_pnlgen_plist_nonlinmotdisp(wall):
    run(wall, "PNLGEN")
    assert [(k, r.get("group"), r.get("bbc")) for k, r in NC.panels(wall.model)] == [(1, 1, 1)]
    run(wall, "PLIST")
    assert any("1 5 15 11" in t for t in wall.sink.texts(Kind.INFO))
    run(wall, "S,1,2,1,1,3,4\nNONLINMOTDISP")
    m = wall.model
    nout = {(n, r.dir) for r in m.nout for n in r.nodes}
    assert {(n, d) for n in (1, 5, 11, 15) for d in (1, 2, 3)} <= nout
    assert (20, 3) in nout and (1, 3) in nout
    rd = {r.node: r.flags for r in m.rdnd}
    assert rd[1][:3] == [1, 1, 1] and rd[20] == [0, 0, 1, 0, 0, 0]
    run(wall, "NONLINMOTDISP")                          # second call: no duplicates (Error 66)
    assert len({r.node for r in m.rdnd}) == len(m.rdnd)
    assert sum(len(r.nodes) for r in m.nout) == len(nout)


def test_shear_and_bbcgen_commands(wall):
    run(wall, "PNLGEN\nSHEAR,1,30000,420000,0.005,0")
    info = " ".join(wall.sink.texts(Kind.INFO))
    assert "SI" in info and "7527.35" in info                          # Gulec-Whittaker capacity, kN
    run(wall, "BBCGEN,1,4,30000,420000,0.005,0,0,0,0")
    assert not errors(wall), errors(wall)
    cv = NC.curve(wall.model, 1)
    assert len(cv["x"]) == 22 and cv["yield_"] == 21 and cv["type"] == 1
    assert cv["x"][-1] == 0.02 and cv["y"][-1] == pytest.approx(1.02 * cv["y"][20])
    G = 3.0e7 / 2.4
    assert cv["y"][0] / cv["x"][0] == pytest.approx(G * 0.3 * 12.0)          # initial slope = G A_W
    run(wall, "BBCGEN,1,4,30000,420000,0.005,0,0,0,0.6")
    assert errors(wall)                                                  # CrackingForceLevel outside [0.10, 0.50]
    run(wall, "P,2,1,0,1,1\nBBCGEN,2,1,30000,420000,0.005,0,0,0,0.3")
    assert wall.model.options.entry("P", 2).get("bbc") == 2              # BBC number = panel number
    cv2 = NC.curve(wall.model, 2)
    assert cv2["y"][0] == pytest.approx(0.3 * cv2["y"][20])
    roundtrip(wall)


# ======================================================================================
# The .eql deck
# ======================================================================================
def test_build_eql_errors(wall):
    run(wall, "EQL,0.8,1")
    deck, errs, _ = NC.build_eql(wall.model)
    assert deck is None and any("Error 121" in e for e in errs)
    run(wall, "P,1,9,1,1,1")
    deck, errs, _ = NC.build_eql(wall.model)
    assert any("Error 123" in e for e in errs)
    run(wall, "P,1,1,1,2,3")
    deck, errs, _ = NC.build_eql(wall.model)
    assert any("Error 126" in e for e in errs) and any("Error 128" in e for e in errs)
    run(wall, "EDUOPT,NONEXT,1")
    deck, errs, _ = NC.build_eql(wall.model)
    assert not any("126" in e or "128" in e for e in errs)


def test_build_eql_material_and_slope_checks(wall):
    run(wall, "P,1,1,1,1,1\nBBCX,1,2,1,1E-4,1E-3\nBBCY,1,2,1,4500,5000\nEQL,0.8,1")
    deck, errs, warns = NC.build_eql(wall.model)
    assert deck is not None and not errs
    row = deck.rows("panels")[0]
    assert (row["bl"], row["br"], row["tr"], row["tl"]) == (1, 5, 15, 11)
    assert row["length"] == pytest.approx(12.0) and row["height"] == pytest.approx(4.0)
    assert row["e_el"] == pytest.approx(3e7) and row["g_el"] == pytest.approx(1.25e7)
    assert not any("EDU-09" in w for w in warns)                        # Y1/X1 = 4.5e7 = G A_shear
    run(wall, "BBCY,1,2,1,4000,5000")
    _, _, warns = NC.build_eql(wall.model)
    assert any("EDU-09" in w and "-11.11 %" in w for w in warns)          # 4.0e7 vs 4.5e7
    # a second group using the same material: D-NON-12
    run(wall, """GROUP,3,SHELL
MACT,1
E,1,1,2,7,6
THICK,1,1,1,0.3""")
    deck, errs, _ = NC.build_eql(wall.model)
    assert deck is None and any("D-NON-12" in e for e in errs)


def test_build_eql_springs(wall):
    run(wall, "S,1,2,1,1,3,4\nBBCX,1,2,1,0.01,0.2\nBBCY,1,2,1,10,29\nEQL,1.0,2")
    deck, errs, warns = NC.build_eql(wall.model)
    assert any("constant of DOF Z is 0" in e for e in errs)
    run(wall, "S,1,2,1,1,1,4")
    deck, errs, warns = NC.build_eql(wall.model)
    assert deck is not None and not errs, errs
    r = deck.rows("springs")[0]
    assert (r["node_i"], r["node_j"], r["k_el"], r["damp"]) == (20, 1, 1000.0, 0.02)
    assert not any("EDU-09" in w for w in warns)                        # slope 10/0.01 = 1000 = k


# ======================================================================================
# Panel-model construction
# ======================================================================================
def box_model(ui):
    """A 2 x 2 x 1 box of shells: 4 walls (2 elements each) and a roof (4 elements), a SOLID below."""
    ui.run_text("""
N,1,0,0,0
N,3,2,0,0
FILL,1,3
NGEN,2,3,1,3,1,0,1,0
NGEN,1,9,1,9,1,0,0,1
M,1,1000,0.2,1,0,0,1
GROUP,1,SHELL
MACT,1
E,1,1,2,11,10
EGEN,1,1,1
E,3,7,8,17,16
EGEN,1,1,3
E,5,1,4,13,10
EGEN,1,3,5
E,7,3,6,15,12
EGEN,1,3,7
E,9,10,11,14,13
EGEN,1,1,9
EGEN,1,3,9,10
THICK,1,12,1,0.1
GROUP,2,SOLID
E,1,1,2,5,4,10,11,14,13
""")
    assert not errors(ui), errors(ui)


def test_wallflr_and_panelize(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    box_model(ui)
    run(ui, "WALLFLR")
    assert not errors(ui), errors(ui)
    m = ui.model
    assert 2 not in m.groups or m.groups[2].type == 3
    assert all(g.type == 3 for g in m.groups.values())
    # the four walls (2 shells each) cannot form groups of >= 5: group 1; the roof (4 shells) too
    assert sorted(m.groups) == [1] and len(m.groups[1].elements) == 12
    assert any("non-shell" in w for w in warnings(ui))


def wall_with_openings(ui):
    """Manual Fig. 1.5: a 6 x 3 shell wall (x = 0..6, z = 0..3) with openings in the middle row at
    columns 1 and 4 (x = 1..2 and 4..5)."""
    ui.run_text("""
N,1,0,0,0
N,7,6,0,0
FILL,1,7
NGEN,3,7,1,7,1,0,0,1
M,1,1000,0.2,1,0,0,1
GROUP,1,SHELL
MACT,1
""")
    k = 0
    for row in range(3):
        for col in range(6):
            if row == 1 and col in (1, 4):
                continue
            k += 1
            n1 = row * 7 + col + 1
            ui.execute(f"E,{k},{n1},{n1 + 1},{n1 + 8},{n1 + 7}")
    ui.execute("THICK,1,16,1,0.2")
    assert not errors(ui), errors(ui)


def test_edge_reproduces_manual_figure_1_5(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    wall_with_openings(ui)
    run(ui, "EDGE,1,0,0,1")
    m = ui.model
    assert sorted(m.groups) == [1, 2, 3]
    zc = {g: np.mean([np.mean([ui.model.node_global(n)[2] for n in e.nodes]) for e in m.groups[g].elements.values()])
          for g in m.groups}
    assert zc[1] > zc[2] > zc[3]                                         # top band keeps group 1
    assert [len(m.groups[g].elements) for g in (1, 2, 3)] == [6, 4, 6]
    run(ui, "EDGE,2")
    assert sorted(m.groups) == [1, 2, 3, 4, 5]
    xc = {g: np.mean([np.mean([ui.model.node_global(n)[0] for n in e.nodes]) for e in m.groups[g].elements.values()])
          for g in (2, 4, 5)}
    assert xc[2] < xc[4] < xc[5]                                         # piers left to right: 2, 4, 5
    assert sum(len(g.elements) for g in m.groups.values()) == 16         # element count conserved
    run(ui, "EDGEMODEL")
    assert not errors(ui)


def test_unipnl_dgrdflr_groupmat(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    box_model(ui)
    run(ui, "GDEL,2\nUNIPNL,1")
    m = ui.model
    assert len(m.groups) == 12 and all(len(g.elements) == 1 for g in m.groups.values())
    run(ui, "GROUPMAT\nDGRDFLR,0.5")
    floors = [g for g in m.groups.values()
              if np.ptp([m.node_global(n)[2] for n in g.elements[1].nodes]) == 0]
    assert len(floors) == 4
    for g in floors:
        assert m.materials[g.elements[1].mat].val1 == pytest.approx(500.0)
    walls = [g for g in m.groups.values() if g not in floors]
    assert all(m.materials[g.elements[1].mat].val1 == 1000.0 for g in walls)


def test_mergepanel(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(WALL)
    ui.run_text("CPMODEL,1\nACTM,1\nGDEL,2")             # panel model: the wall only
    ui.run_text("ACTM,0\nGDEL,1")                        # original: delete its shells
    run(ui, "MERGEPANEL,1")
    assert not errors(ui), errors(ui)
    m = ui.model
    assert sorted(m.groups) == [2, 3] and m.groups[3].type == 3 and len(m.groups[3].elements) == 8
    assert m.groups[3].elements[1].mat == 3                              # materials appended after 1, 2


def test_solidpile(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text("""
N,1,0,0,0
N,3,2,0,0
FILL,1,3
NGEN,2,3,1,3,1,0,1,0
NGEN,2,9,1,9,1,0,0,1
M,1,1000,0.3,1,0.05,0.05,1
GROUP,1,SOLID
MACT,1
E,1,1,2,5,4,10,11,14,13
EGEN,1,1,1
EGEN,1,3,1,2
EGEN,1,9,1,4
GROUP,2,SOLID
E,1,1,2,5,4,10,11,14,13
""")
    # group 2 (the 'pile') coincides with element 1 of the soil: interface = its 8 nodes
    run(ui, "SOLIDPILE,2")
    assert not errors(ui), errors(ui)
    m = ui.model
    springs = [g for g in m.groups.values() if g.type == 7]
    assert len(springs) == 4
    nsp = sum(len(g.elements) for g in springs)
    tip = [g for g in springs if "tip" in g.title][0]
    assert len(tip.elements) == 4 and nsp == 4 * 3 + 4                    # 4 tip nodes, 4 side nodes x 3
    assert len(m.springs) == nsp and all(sp.damp == 0.04 for sp in m.springs.values())
    pile_nodes = {n for e in m.groups[2].elements.values() for n in e.nodes}
    assert pile_nodes.isdisjoint(set(range(1, 28)))                       # reconnected to the duplicates


# ======================================================================================
# Batch, COMB_XYZ_THD and run control
# ======================================================================================
def test_nonlinbat_and_comb_inp(wall, tmp_path):
    (tmp_path / "w").mkdir()
    run(wall, f"MDL,w,{tmp_path / 'w'}\nANALYS,0,0,0,0,1,0,0,0,0,0,0,1\nPNLGEN\nNONLINBAT,1")
    assert not errors(wall), errors(wall)
    text = (tmp_path / "w" / "w_NONLINBAT.pre").read_text()
    assert "NONLINITER,NLHOUSE+NLDX+NLDY+NLDZ+NLNON,10" in text and "ANALYS,0,0,1,0,1" in text
    assert "ANALYS,0,0,0,1,1" in text and "COMBXYZTHD" in text
    # the loop's AFWRITE writes MOTION/RELDISP only: AOPT is the first command of each direction's variable,
    # and nothing changes the model between the AFWRITE of the ANALYS deck and NONLINITER
    assert 'VAR,NLDX,"AOPT,0,0,0,0,0,0,0,0,0,0,1,0,1,0","EDUOPT,TFFILE,FILE8X",AFWRITE' in text
    lines = text.splitlines()
    k = lines.index("NONLINITER,NLHOUSE+NLDX+NLDY+NLDZ+NLNON,10")
    assert lines[k - 2] == "AFWRITE" and lines[k - 1].startswith("VAR,NLHOUSE,")
    inp = (tmp_path / "w" / "COMB_XYZ_THD.inp").read_text().split("\n")
    assert inp[0] == "12" and inp[1] == "00001TR_X.THD X_00001TR_X.THD Y_00001TR_X.THD Z_00001TR_X.THD"
    run(wall, "NONLINBAT,0")
    assert "NLHOUSE+NLD+NLNON" in (tmp_path / "w" / "w_NONLINBAT.pre").read_text()


def _hist(path, vals, dt=0.01):
    path.write_text("\n".join([f"{dt}"] + [f"{v:.16e}" for v in vals]) + "\n")


def test_combxyzthd_and_nonlinthd(wall, tmp_path):
    d = tmp_path / "w"
    d.mkdir()
    run(wall, f"MDL,w,{d}\nPNLGEN")
    names = NC._thd_names(wall.model)
    rng = np.random.default_rng(1)
    data = {n: rng.normal(size=(3, 50)) for n in names}
    # one RELDISP run per direction writes the plain names; NONLINTHD keeps them as X_*, Y_*, Z_*
    for k, p in enumerate("XYZ"):
        for n in names:
            _hist(d / n, data[n][k])
        run(wall, f"NONLINTHD,{p}")
        assert not errors(wall), errors(wall)
    (d / "COMB_XYZ_THD.inp").write_text(NC.comb_inp_text(wall.model))
    run(wall, "COMBXYZTHD")
    assert not errors(wall), errors(wall)
    from sassi.modules.nonlinear import read_history
    for n in names:
        v, dt = read_history(d / n)
        assert dt == 0.01
        np.testing.assert_allclose(v, data[n].sum(axis=0), rtol=1e-15, atol=1e-15)
    (d / f"Y_{names[0]}").unlink()
    run(wall, "COMBXYZTHD")
    assert errors(wall)                                                  # a missing input is an error


def test_nonlinreset_and_nonlinsave(wall, tmp_path):
    d = tmp_path / "w"
    d.mkdir()
    run(wall, f"MDL,w,{d}")
    for nm in ("PANEL.NON", "Panel_EQL_Matl_Prop.txt", "NONLINEAR_CONVERGENCE.TXT", "Panel0001.thd", "Panel.fmu",
               "w.hou", "FILE8X"):
        (d / nm).write_text("1\n")
    run(wall, "NONLINSAVE,It3,1")
    for nm in ("Panel0001_It3.thd", "Panel_It3.fmu", "Panel_EQL_Matl_Prop_It3.txt", "w_It3.hou", "FILE8X_It3"):
        assert (d / nm).exists(), nm
    run(wall, "NONLINRESET")
    assert not (d / "PANEL.NON").exists() and not (d / "Panel_EQL_Matl_Prop.txt").exists()
    assert not (d / "NONLINEAR_CONVERGENCE.TXT").exists() and (d / "Panel0001.thd").exists()


def test_nonliniter_requires_nonlinear_runs(wall, tmp_path):
    d = tmp_path / "w"
    d.mkdir()
    run(wall, f'MDL,w,{d}\nVAR,NLB,"TIT,pass"\nNONLINITER,NLB,2')
    assert any("did not run" in e for e in errors(wall))
    run(wall, 'VAR,NLC,"TIT,pass",NOSUCHCOMMAND\nNONLINITER,NLC,2')
    assert any("failed; iterations stopped" in e for e in errors(wall))
    (d / "NONLINEAR_CONVERGENCE.TXT").write_text("# h\n   0  1  0.1  0.0  1\n")
    run(wall, "NONLINITER,NLB,2")
    assert any("already converged" in t for t in wall.sink.texts(Kind.INFO))
