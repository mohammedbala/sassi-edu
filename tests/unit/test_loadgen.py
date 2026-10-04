"""LOADGEN module (Option A, sassi/modules/loadgen.py): deck, numbering, masses, histories, critical times,
APDL writers and complete static / dynamic runs on a small stick model (D-LGN-01 ... D-LGN-08).

The run tests share one model directory built once per module (a three-mass stick on a 3 x 3 rigid mat
on soil springs, FILE8 exact at every Fourier frequency: sassi.verify.problems.vp_loadgen.spring_model),
copied into each test's tmp_path."""
from __future__ import annotations

import io
import shutil
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from sassi.conventions import nodal_result_name
from sassi.io import textfiles
from sassi.io.container import read_container
from sassi.modules import loadgen as LG
from sassi.modules.base import ModuleError
from sassi.modules.motion import write_frame
from sassi.verify import builders as B
from sassi.verify.problems import vp_loadgen as V

G = 9.81
NFFT, DT = 512, 0.02
STICK_NODES = (12, 14, 16)            # stick_on_mat numbering: mat 1-9 (centre 5), stick 12, 14, 16
MAT_NODES = tuple(range(1, 10))


# ======================================================================================
# fixtures
# ======================================================================================
def _build(wd: Path, optimize: bool = False):
    fs, mdl = V.spring_model(wd, B, nfft=NFFT, dt=DT, optimize=optimize)
    V.control_motion(wd, n=300, peak=0.3, seed=3)
    if not optimize:
        nodes = [(n, (1, 2, 3, 4, 5, 6)) for n in (5,) + STICK_NODES] + [(n, (1, 2, 3, 4, 5, 6)) for n in MAT_NODES
                                                                         if n != 5]
        rc1, rc2 = V.run_motion_reldisp(wd, B, fs, nodes, MAT_NODES, rd_dofs=(1, 2, 3, 4, 5, 6), savetf=0)
        assert rc1 == 0 and rc2 == 0
    return fs, mdl


@pytest.fixture(scope="module")
def model_dir(tmp_path_factory):
    wd = tmp_path_factory.mktemp("lgmodel")
    _build(wd)
    yield wd
    shutil.rmtree(wd, ignore_errors=True)


@pytest.fixture
def wd(model_dir, tmp_path):
    dst = tmp_path / "m"
    shutil.copytree(model_dir, dst)
    return dst


def run(wd: Path, **params):
    rc, lst = V.run_loadgen(wd, **params)
    return rc, lst


def hist(wd: Path, node: int, dof: int, ext: str) -> np.ndarray:
    return textfiles.read_history(wd / nodal_result_name(node, dof, ext))[0]


# ======================================================================================
# deck
# ======================================================================================
def test_deck_round_trip(tmp_path):
    d = LG.new_deck()
    d["model"], d["analysis"], d["data"], d["multi"] = "m", "DYNAMIC", 3, 1
    d["alpha"], d["beta"], d["refpoint"], d["times"] = 0.25, 1.5e-3, [1.0, 2.0, -3.5], [1.5, 2.25]
    d["apdlfile"], d["groundfile"], d["crit"] = "out dir/x y.inp", "C:/a b/eq.acc", "TIME"
    d.tables["dnodes"] = [3, 1, 2]
    d.tables["checknodes"] = [12]
    p = LG.write_deck(tmp_path / "m.lgn", d, comments=["test"])
    assert p.read_text().startswith("SASSI-EDU LOADGEN DECK v1")
    e, notes = LG.read_deck(p)
    assert notes == []
    assert e.params == d.params and e.tables == {"dnodes": [3, 1, 2], "masters": [], "checknodes": [12]}


def test_deck_errors(tmp_path):
    d = LG.new_deck()
    with pytest.raises(KeyError):
        d["nonsense"] = 1
    p = LG.write_deck(tmp_path / "m.lgn", d)
    text = p.read_text().replace("[params]", "[params]\nextra = 1")
    p.write_text(text)
    _, notes = LG.read_deck(p)
    assert any("extra" in n for n in notes)
    p.write_text(text.replace("LOADGEN DECK", "MOTION DECK"))
    with pytest.raises(ValueError, match="expected a LOADGEN deck"):
        LG.read_deck(p)
    bad = LG.new_deck()
    bad["data"] = 7
    with pytest.raises(ModuleError, match="data"):
        LG.validate_deck(bad)
    bad = LG.new_deck()
    bad["crit"] = "XX"
    with pytest.raises(ModuleError, match="criterion"):
        LG.validate_deck(bad)


# ======================================================================================
# axes, numbering, node maps
# ======================================================================================
def test_ansys_dof_mapping():
    assert [LG.ansys_dof(k, False) for k in range(1, 7)] == [(k - 1, 1.0) for k in range(1, 7)]
    # 2-D: x_A = R x_S with R = [[1,0,0],[0,0,1],[0,-1,0]]
    got = {k: LG.ansys_dof(k, True) for k in range(1, 7)}
    assert got == {1: (0, 1.0), 2: (2, -1.0), 3: (1, 1.0), 4: (3, 1.0), 5: (5, -1.0), 6: (4, 1.0)}
    assert np.allclose(LG.control_direction(0, 30.0), [np.cos(np.pi / 6), np.sin(np.pi / 6), 0])
    assert np.allclose(LG.control_direction(2, 0.0), [0, 0, 1])


def test_two_d_rotation_matches_the_ansys_export():
    """The 2-D rotation of LOADGEN is the one of the ANSYS command (sassi.io.apdl): node coordinates of
    the exported APDL = R_2D x_SASSI."""
    from sassi.io.apdl import export_apdl
    from sassi.prep import Interpreter
    ui = Interpreter()
    ui.run_text("N,1,0,0,0\nN,2,1,0,0\nN,3,1,0,2\nN,4,0,0,2\nM,1,3e7,0.2,24,0.05,0.05,1\nGROUP,1,PLANE\nMACT,1\n"
                "E,1,1,2,3,4\nETYPE,1,1,1,1\n")
    text = export_apdl(ui.model).text
    nodes = {int(t[1]): np.array([float(x) for x in t[2:5]]) for t in
             (ln.split(",") for ln in text.splitlines() if ln.startswith("N,"))}
    for n in (1, 2, 3, 4):
        assert np.allclose(nodes[n], LG.R_2D @ ui.model.node_global(n))


def _info(xyz, two_d=False, old=None):
    ids = sorted(xyz)
    old = old or {n: n for n in ids}
    return LG.ModelInfo(model="m", gravity=G, dim=1 if two_d else 2, two_d=two_d,
                        xyz={n: np.asarray(xyz[n], float) for n in ids}, model_of=dict(old),
                        file4_of={o: n for n, o in old.items()}, dofs={n: [1, 2, 3] for n in ids},
                        eq_node=np.repeat(ids, 3), eq_dof=np.tile([1, 2, 3], len(ids)), int_nodes=[ids[0]],
                        excav_nodes=[], optimized=False, meta={})


def test_pairs_file(tmp_path):
    p = tmp_path / "map.txt"
    p.write_text("! sassi ansys\n* comment\n1 101\n2, 102   # trailing\n\n3 103\n")
    assert LG.read_pairs(p) == {1: 101, 2: 102, 3: 103}
    p.write_text("1 101\n1 102\n")
    with pytest.raises(ModuleError, match="mapped twice"):
        LG.read_pairs(p)
    p.write_text("1\n")
    with pytest.raises(ModuleError, match="expected"):
        LG.read_pairs(p)


def test_coordinate_map_3d_and_2d():
    info = _info({1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 0, 2)})
    ans = {101: np.array([0.0, 0, 0]), 102: np.array([1.0, 0, 1e-9]), 7: np.array([5.0, 5, 5])}
    table, missing = LG.coordinate_map(info, ans, tol=1e-6)
    assert table == {1: 101, 2: 102} and missing == 1
    info2 = _info({1: (0, 0, 0), 2: (1, 0, 2)}, two_d=True)
    ans2 = {11: np.array([0.0, 0, 0]), 12: np.array([1.0, 2.0, 0.0])}       # ANSYS X-Y plane: Y = SASSI Z
    assert LG.coordinate_map(info2, ans2, tol=1e-6) == ({1: 11, 2: 12}, 0)


def test_master_assignment_by_floor():
    info = _info({1: (0, 0, 0), 2: (5, 0, 4), 3: (0, 0, 4), 4: (0, 0, 8), 5: (9, 9, 8.0000001), 6: (3, 0, 4)})
    a = LG.master_assignment(info, [2, 5, 6], [3, 4, 1])
    assert a == {2: 3, 5: 4, 6: 3}
    with pytest.raises(ModuleError, match="not a node"):
        LG.master_assignment(info, [2], [99])


# ======================================================================================
# masses
# ======================================================================================
def test_lumped_from_consistent_beam_matrix():
    """Row sums per direction of a consistent beam mass give the tributary masses and the exact total."""
    L, rhoA = 2.0, 3.0
    m = rhoA * L / 420.0
    # 2-node planar beam (X-Z bending), DOFs per node: UZ (3), ROTY (5); axial UX (1) lumped
    Mb = m * np.array([[156, 22 * L, 54, -13 * L], [22 * L, 4 * L * L, 13 * L, -3 * L * L],
                       [54, 13 * L, 156, -22 * L], [-13 * L, -3 * L * L, -22 * L, 4 * L * L]])
    eq_node = np.array([1, 1, 1, 2, 2, 2])
    eq_dof = np.array([1, 3, 5, 1, 3, 5])
    M = np.zeros((6, 6))
    M[np.ix_([1, 2, 4, 5], [1, 2, 4, 5])] = Mb
    M[0, 0] = M[3, 3] = rhoA * L / 2
    info = _info({1: (0, 0, 0), 2: (2, 0, 0)})
    info.eq_node, info.eq_dof = eq_node, eq_dof
    lm = LG.lumped_from_matrix(info, sp.csr_matrix(M))
    assert np.isclose(lm[1][2] + lm[2][2], rhoA * L)                              # total vertical mass
    assert np.isclose(lm[1][2], rhoA * L / 2) and np.isclose(lm[1][0], rhoA * L / 2)
    assert np.isclose(lm[1][4], 4 * L * L * m)                                    # rotary = diagonal


def test_lumped_and_master_files(tmp_path):
    info = _info({1: (0, 0, 0)})
    p = LG.write_lumped(tmp_path / "a.masl", {3: np.array([1.0, 2, 3, 0, 0, 0.5])}, info, "test")
    back = LG.read_lumped(p)
    assert list(back) == [3] and np.allclose(back[3], [1, 2, 3, 0, 0, 0.5])
    p.write_text("! comment\n7 2.5\n8 1 2 3\n")
    got = LG.read_lumped(p)
    assert np.allclose(got[7], [2.5, 2.5, 2.5, 0, 0, 0]) and np.allclose(got[8], [1, 2, 3, 0, 0, 0])
    p.write_text("7 1 2\n")
    with pytest.raises(ModuleError, match="expected"):
        LG.read_lumped(p)
    p.write_text("7 -1 1 1\n")
    with pytest.raises(ModuleError, match="negative"):
        LG.read_lumped(p)
    ms = [LG.LoadMass(1001, np.array([1.0, 1, 1, 0, 0, 0]), master=5, xyz=np.array([1.0, 2, 3]))]
    q = LG.write_master(tmp_path / "a.masm", ms, info, "test")
    back = LG.read_master(q)
    assert back[0].node == 1001 and back[0].master == 5 and np.allclose(back[0].xyz, [1, 2, 3])
    q.write_text("1 2 3\n")
    with pytest.raises(ModuleError, match="expected"):
        LG.read_master(q)


# ======================================================================================
# histories, peaks, APDL helpers
# ======================================================================================
def test_pick_peaks():
    y = np.array([0.0, 1.0, 0.5, -3.0, 0.2, 2.0, 2.0, 0.1, -2.5, 0.0])
    assert LG.pick_peaks(y, 1) == [3]
    assert LG.pick_peaks(y, 3) == [3, 8, 5]                     # the plateau keeps its first sample
    assert LG.pick_peaks(y, 3, sep=3) == [3, 8]                 # 5 is too close to 3 and 8
    assert LG.pick_peaks(np.zeros(5), 2) == []
    assert LG.pick_peaks(np.array([5.0, 1.0, 0.0]), 1) == [0]   # an end sample can be the maximum


def test_start_at_rest():
    hs = LG.HistorySet(dt=0.1, n=4, disp={(1, 1): np.array([0.2, 1.2, -0.8, 0.2])})
    assert LG.start_at_rest(hs) == pytest.approx(0.2 / 1.2)
    assert np.allclose(hs.disp[(1, 1)], [0.0, 1.0, -1.0, 0.0])


def test_fmt_value_and_table_lines_parse_back():
    assert LG.fmt_value(0.0) == "0" and LG.fmt_value(-1.5, 6) == "-1.50000E+00"
    with pytest.raises(ModuleError, match="non-finite"):
        LG.fmt_value(np.nan)
    rng = np.random.default_rng(1)
    v = rng.standard_normal(23) * 1e-3
    lines = ["LG_DT = 0.02"] + LG.table_lines("LGD_7_UX", v, 12, "test")
    assert lines[1] == "*DIM,LGD_7_UX,TABLE,23,1,1,TIME   ! test"
    assert lines[2] == "*VFILL,LGD_7_UX(1,0),RAMP,0,LG_DT"
    assert len(lines) == 3 + 3 and lines[5].startswith("LGD_7_UX(21,1)=")
    prog = V.parse_apdl("\n".join(lines))
    assert prog.errors == []
    t = prog.table("LGD_7_UX")
    assert np.allclose(t.column(1), v, rtol=1e-11, atol=0) and np.array_equal(t.index, np.arange(23) * 0.02)
    with pytest.raises(ModuleError, match="longer"):
        LG.apdl_name("LGD", 123456, "X" * 30)


def test_test_parser_reports_problems():
    text = "\n".join(["*DIM,T,TABLE,3,1,1,TIME", "T(1,1)=1,2", "D,1,UQ,%T%", "D,1,UX,%NOPE%", "FOO,1",
                      "*IF,A,EQ,1,THEN", "x" * 700])
    errs = V.parse_apdl(text).errors
    assert any("not filled" in e for e in errs)
    assert any("invalid D label" in e for e in errs)
    assert any("NOPE" in e for e in errs)
    assert any("unexpected command FOO" in e for e in errs)
    assert any("700 characters" in e for e in errs)
    assert any("unterminated" in e for e in errs)
    prog = V.parse_apdl("A = 1\n*IF,A,EQ,1,THEN\nSOLVE\n*ENDIF\n*IF,A,EQ,0,THEN\nFINISH\n*ENDIF")
    assert [c.name for c in prog.commands] == ["SOLVE"] and prog.errors == []
    prog = V.parse_apdl("*DIM,T,TABLE,2,1,1,TIME\nT(1,1)=" + ",".join(["1"] * 11))
    assert any("11 values" in e for e in prog.errors)


def test_read_frames(tmp_path):
    d = tmp_path / "ACC"
    d.mkdir()
    for s in range(3):
        write_frame(d / f"ACC_{s * 0.01:06.3f}_{s + 1:05d}", np.array([4, 7]), np.array([[s, 0, 1], [2 * s, 1, 0]]))
    steps, rows, vals = LG.read_frames(d)
    assert list(steps) == [1, 2, 3] and rows == {4: 0, 7: 1}
    assert np.allclose(vals[:, 1, 0], [0, 2, 4])
    (d / "ACC_00.010_00002").unlink()
    with pytest.raises(ModuleError, match="missing frames"):
        LG.read_frames(d)


# ======================================================================================
# static runs
# ======================================================================================
def test_static_acceleration_fixed_base(wd):
    rc, lst = run(wd, analysis="STATIC", data=2, genmass=1)
    assert rc == 0, lst
    prog = V.parse_file(wd / "m_LGS.inp")
    assert prog.errors == []
    d = prog.loads("D")
    assert {n for n, _, _ in d} == set(MAT_NODES) and all(v == 0.0 for _, _, v in d)
    assert {lab for _, lab, _ in d} == {"UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"}       # clamped rigid mat
    f = prog.loads("F")
    assert {n for n, _, _ in f} == set(STICK_NODES)
    res = np.loadtxt(wd / "m_LGS_res.txt")
    k = int(np.argmax(np.abs(res[:, 1])))
    assert sum(v for _, lab, v in f if lab == "FX") == pytest.approx(res[k, 1], rel=1e-10)
    # F = -m a with the MOTION .ACC of the masses (mass 2 t at node 12)
    a12 = hist(wd, 12, 1, "ACC")[k] * G
    assert [v for n, lab, v in f if n == 12 and lab == "FX"][0] == pytest.approx(-2.0 * a12, rel=1e-10)
    assert prog.scalars["LG_SOLVE"] == 1 and len(prog.find("SOLVE")) == 1
    assert [c.fields[0] for c in prog.find("ANTYPE")] == ["STATIC"]
    assert "Acceleration" in lst and (wd / "m.masl").exists()


def test_static_displacement_only_with_disp_criterion(wd):
    rc, lst = run(wd, analysis="STATIC", data=1, crit="DISP", critnode=1, critdof=3)
    assert rc == 0, lst
    prog = V.parse_file(wd / "m_LGS.inp")
    assert prog.errors == [] and prog.loads("F") == []
    thd = hist(wd, 1, 3, "THD")
    u = thd - thd[0]
    k = int(np.argmax(np.abs(u)))
    d = {(n, lab): v for n, lab, v in prog.loads("D")}
    assert len(d) == 27                                             # translations only (rotdisp 0)
    assert d[(1, "UZ")] == pytest.approx(u[k], rel=1e-11)
    assert "rotations are not prescribed" in lst
    rc, lst = run(wd, analysis="STATIC", data=1, crit="DISP", critnode=1, critdof=3, rotdisp=1, apdlfile="r.inp")
    assert rc == 0, lst
    d = {(n, lab): v for n, lab, v in V.parse_file(wd / "r.inp").loads("D")}
    ry = hist(wd, 1, 5, "THD")
    assert len(d) == 54 and d[(1, "ROTY")] == pytest.approx((ry - ry[0])[k], rel=1e-10)


def test_static_without_start_at_rest(wd):
    rc, lst = run(wd, analysis="STATIC", data=1, crit="TIME", times=[2.0], rest=0)
    assert rc == 0 and "do not start at rest" in lst
    d = {(n, lab): v for n, lab, v in V.parse_file(wd / "m_LGS.inp").loads("D")}
    assert d[(1, "UZ")] == pytest.approx(hist(wd, 1, 3, "THD")[100], rel=1e-11)       # the RELDISP value itself


def test_static_soil_module_default_nodes(wd):
    rc, lst = run(wd, analysis="STATIC", data=4, crit="TIME", times=[1.0])
    assert rc == 0, lst
    assert "Disp. for Soil Module" in lst
    assert {n for n, _, _ in V.parse_file(wd / "m_LGS.inp").loads("D")} == set(MAT_NODES)


def test_static_given_times_and_multiple_files(wd):
    rc, lst = run(wd, analysis="STATIC", data=3, genmass=1, crit="TIME", times=[2.0, 1.005], multi=1)
    assert rc == 0, lst
    assert "not on the time grid" in lst
    p1, p2 = wd / "m_LGS_01.inp", wd / "m_LGS_02.inp"
    assert p1.exists() and p2.exists() and not (wd / "m_LGS.inp").exists()
    a, b = V.parse_file(p1), V.parse_file(p2)
    assert a.errors == b.errors == []
    assert a.scalars["LG_TSSI"] == pytest.approx(2.0) and b.scalars["LG_TSSI"] == pytest.approx(1.0)
    assert [c.fields[0] for c in a.find("TIME")] == ["1"] and [c.fields[0] for c in b.find("TIME")] == ["2"]
    rc, lst = run(wd, analysis="STATIC", data=3, genmass=1, crit="STEP", times=[5, 9], multi=0, apdlfile="one.inp")
    assert rc == 0 and "only the first" in lst
    assert V.parse_file(wd / "one.inp").scalars["LG_TSSI"] == pytest.approx(4 * DT)
    rc, lst = run(wd, analysis="STATIC", data=3, genmass=1, crit="TIME", times=[99.0])
    assert rc == 1 and "outside the histories" in lst


def test_static_peaks_criteria(wd):
    rc, lst = run(wd, analysis="STATIC", data=2, genmass=1, crit="MY", ncrit=3, tsep=0.5, multi=1,
                  refpoint=[0.0, 0.0, 0.0])
    assert rc == 0, lst
    res = np.loadtxt(wd / "m_LGS_res.txt")
    my = res[:, 5]
    ts = sorted(V.parse_file(wd / f"m_LGS_{i:02d}.inp").scalars["LG_TSSI"] for i in (1, 2, 3))
    ks = [int(round(t / DT)) for t in ts]
    assert int(np.argmax(np.abs(my))) in ks and min(np.diff(sorted(ks))) * DT >= 0.5
    rc, lst = run(wd, analysis="STATIC", data=2, genmass=1, crit="ACC", critnode=16, critdof=1)
    assert rc == 0, lst
    a = hist(wd, 16, 1, "ACC")
    assert V.parse_file(wd / "m_LGS.inp").scalars["LG_TSSI"] == pytest.approx(int(np.argmax(np.abs(a))) * DT)


def test_static_data_check_writes_nothing(wd):
    rc, lst = run(wd, analysis="STATIC", data=3, genmass=1, opmode=1)
    assert rc == 0 and "data-check mode" in lst
    assert not (wd / "m.masl").exists() and not list(wd.glob("m_LGS*"))


def test_static_mass_file_needed(wd):
    rc, lst = run(wd, analysis="STATIC", data=2, genmass=0)
    assert rc == 1 and "Generate Mass Data" in lst
    (wd / "user.masl").write_text("12 10.0\n16 1 2 3\n")
    rc, lst = run(wd, analysis="STATIC", data=2, genmass=0, lumpfile="user.masl")
    assert rc == 0, lst
    f = {(n, lab): v for n, lab, v in V.parse_file(wd / "m_LGS.inp").loads("F")}
    k = int(round(V.parse_file(wd / "m_LGS.inp").scalars["LG_TSSI"] / DT))
    assert f[(12, "FX")] == pytest.approx(-10.0 * G * hist(wd, 12, 1, "ACC")[k], rel=1e-10)
    assert f[(16, "FZ")] == pytest.approx(-3.0 * G * hist(wd, 16, 3, "ACC")[k], rel=1e-10)
    (wd / "user.masl").write_text("99 1.0\n")
    rc, lst = run(wd, analysis="STATIC", data=2, genmass=0, lumpfile="user.masl")
    assert rc == 1 and "mass node 99" in lst


def test_master_node_mass_rigid_body_motion(wd):
    """A load node 1 m beside the top master node follows its rigid-body motion a = a_M + alpha x r."""
    (wd / "user.masm").write_text("! node master x y z mx my mz\n5001 16 1.0 0.0 12.0 2.0 2.0 2.0\n")
    rc, lst = run(wd, analysis="STATIC", data=2, masstype=2, genmass=0, masterfile="user.masm", rotacc=1,
                  crit="TIME", times=[2.0])
    assert rc == 0, lst
    f = {(n, lab): v for n, lab, v in V.parse_file(wd / "m_LGS.inp").loads("F")}
    k = int(round(2.0 / DT))
    aM = np.array([hist(wd, 16, d, "ACC")[k] for d in (1, 2, 3)]) * G
    al = np.array([hist(wd, 16, d, "ACC")[k] for d in (4, 5, 6)]) * G
    a = aM + np.cross(al, [1.0, 0.0, 0.0])
    for i, lab in enumerate(("FX", "FY", "FZ")):
        assert f[(5001, lab)] == pytest.approx(-2.0 * a[i], rel=1e-9, abs=1e-12)
    assert abs(al[1]) > 0                                     # rocking: the rotational term is exercised
    rc, lst = run(wd, analysis="STATIC", data=2, masstype=2, genmass=0, masterfile="user.masm", crit="TIME",
                  times=[2.0], apdlfile="norot.inp")
    assert rc == 0 and "master rotations are taken as zero" in lst


def test_master_node_mass_generation(wd):
    rc, lst = run(wd, analysis="STATIC", data=2, masstype=2, genmass=1, tables={"masters": [5, 16]})
    assert rc == 0, lst
    ms = LG.read_master(wd / "m.masm")
    assert {m.node: m.master for m in ms} == {12: 5, 14: 16, 16: 16}        # closest in elevation (z 0 / 12)
    rc, lst = run(wd, analysis="STATIC", data=2, masstype=2, genmass=1)
    assert rc == 1 and "LGNODE,M" in lst


def test_frames_fallback(wd):
    """Without the .ACC files the MOTION frames ACC/ACCR give the same loads (11 digits)."""
    rc, _ = run(wd, analysis="STATIC", data=2, genmass=1, apdlfile="files.inp")
    assert rc == 0
    ref = {(n, lab): v for n, lab, v in V.parse_file(wd / "files.inp").loads("F")}
    acc = {n: [hist(wd, n, d, "ACC") for d in (1, 2, 3)] for n in STICK_NODES}
    (wd / "ACC").mkdir()
    for s in range(len(acc[12][0])):
        write_frame(wd / "ACC" / f"ACC_{s * DT:06.3f}_{s + 1:05d}", np.array(STICK_NODES),
                    np.array([[acc[n][d][s] for d in range(3)] for n in STICK_NODES]))
    for n in STICK_NODES:
        for d in (1, 2, 3):
            (wd / nodal_result_name(n, d, "ACC")).unlink()
    rc, lst = run(wd, analysis="STATIC", data=2, genmass=1, apdlfile="frames.inp")
    assert rc == 0, lst
    assert "frames" in lst
    got = {(n, lab): v for n, lab, v in V.parse_file(wd / "frames.inp").loads("F")}
    peak = max(abs(v) for v in ref.values())
    assert max(abs(got[k] - ref[k]) for k in ref) <= 1e-9 * peak


def test_tfi_fallback_equals_reldisp(wd):
    """Without RELDISP's .THD the relative displacements come from MOTION's complex .TFI (RELDISP formula)."""
    rc, _ = run(wd, analysis="DYNAMIC", rotdisp=1, apdlfile="thd.inp")
    assert rc == 0
    ref = V.parse_file(wd / "thd.inp")
    for p in wd.glob("*.THD"):
        p.unlink()
    rc, lst = run(wd, analysis="DYNAMIC", rotdisp=1, apdlfile="tfi.inp")
    assert rc == 0, lst
    assert "computed from MOTION's complex .TFI" in lst
    got = V.parse_file(wd / "tfi.inp")
    names = [nm for (_, _, nm) in ref.loads("D")]
    assert len(names) == 54
    scale = max(float(np.max(np.abs(ref.table(nm).column(1)))) for nm in names)
    assert max(float(np.max(np.abs(got.table(nm).column(1) - ref.table(nm).column(1)))) for nm in names) <= 1e-10 * scale


def test_results_missing_files_message(wd):
    for p in wd.glob("*R_YY.ACC"):
        p.unlink()
    rc, lst = run(wd, analysis="STATIC", data=2, masstype=2, genmass=0, masterfile="x.masm")
    assert rc == 1                                                  # mass file missing first
    (wd / "x.masm").write_text("5001 16 1 0 12 1 1 1\n")
    rc, lst = run(wd, analysis="STATIC", data=2, masstype=2, masterfile="x.masm", rotacc=1)
    assert rc == 1 and "16 R_YY (.ACC)" in lst and "source FILE8" in lst
    rc, lst = run(wd, analysis="STATIC", data=2, masstype=2, masterfile="x.masm", rotacc=1, source="FILE8")
    assert rc == 0, lst


def test_node_maps(wd):
    (wd / "pairs.txt").write_text("\n".join(f"{n} {n + 1000}" for n in MAT_NODES + STICK_NODES if n != 9) + "\n")
    rc, lst = run(wd, analysis="STATIC", data=3, genmass=1, mapmode="PAIRS", mapfile="pairs.txt")
    assert rc == 1 and "no ANSYS node in the node map" in lst and ": 9" in lst
    (wd / "pairs.txt").write_text("\n".join(f"{n} {n + 1000}" for n in MAT_NODES + STICK_NODES) + "\n")
    rc, lst = run(wd, analysis="STATIC", data=3, genmass=1, mapmode="PAIRS", mapfile="pairs.txt")
    assert rc == 0, lst
    prog = V.parse_file(wd / "m_LGS.inp")
    assert {n for n, _, _ in prog.loads("D")} == {n + 1000 for n in MAT_NODES}
    assert {n for n, _, _ in prog.loads("F")} == {n + 1000 for n in STICK_NODES}
    # COORD: an APDL model with other node numbers at the same places
    f4 = read_container(wd / "m.N4", "FILE4")
    lines = [f"N,{int(n) + 500},{float(x[0])!r},{float(x[1])!r},{float(x[2])!r}"
             for n, x in zip(f4["node_id"], f4["node_xyz"])]
    (wd / "ansys.inp").write_text("/PREP7\n" + "\n".join(lines) + "\nN,9999,50,50,50\nFINISH\n")
    rc, lst = run(wd, analysis="STATIC", data=3, genmass=1, mapmode="COORD", mapfile="ansys.inp")
    assert rc == 0, lst
    prog = V.parse_file(wd / "m_LGS.inp")
    assert {n for n, _, _ in prog.loads("D")} == {n + 500 for n in MAT_NODES}
    assert LG.read_pairs(wd / "m_LG.map")[12] == 512


# ======================================================================================
# dynamic runs
# ======================================================================================
def test_dynamic_rel_reference_node_file8(wd):
    rc, lst = run(wd, analysis="DYNAMIC", method="REL", refnode=5, source="FILE8", rotdisp=1, alpha=0.3, beta=1e-3)
    assert rc == 0, lst
    prog = V.parse_file(wd / "m_LGD.inp")
    assert prog.errors == []
    for lab, d in (("X", 1), ("Y", 2), ("Z", 3)):              # ACEL = absolute acceleration of the mat centre
        a = hist(wd, 5, d, "ACC")
        if np.max(np.abs(a)) > 1e-12:
            assert V._relmax(prog.table(f"LG_AC{lab}").column(1) / G, a) < 1e-9
    assert np.max(np.abs(prog.table("LGD_5_UX").column(1))) == 0.0          # the reference itself
    # translations relative to the centre, rotations absolute (= free-field THD rotations, reference rotation 0)
    for n in (1, 9):
        for lab, dof in (("UX", 1), ("UZ", 3)):                 # UX: rigid mat (tiny); UZ: rocking
            u = hist(wd, n, dof, "THD") - hist(wd, 5, dof, "THD")
            scale = np.max(np.abs(hist(wd, n, dof, "THD")))
            err = np.max(np.abs(prog.table(f"LGD_{n}_{lab}").column(1) - (u - u[0])))
            assert err <= 1e-10 * scale
        ry = hist(wd, n, 5, "THD")
        assert V._relmax(prog.table(f"LGD_{n}_ROTY").column(1), ry - ry[0]) < 1e-9
    assert np.max(np.abs(prog.table("LGD_1_UZ").column(1))) > 1e3 * np.max(np.abs(prog.table("LGD_1_UX").column(1)))


def test_reference_acceleration_and_baseline_correction(wd):
    """MOTION <bl> 1: method ACC follows MOTION (corrected reference acceleration); method REL keeps the exact
    second derivative of the reference motion its D tables refer to (uncorrected)."""
    from sassi.io import decks
    md = decks.read(wd / "m.mot", "MOTION")
    md["bl"] = 1
    decks.write(wd / "m.mot", md)                    # the existing .ACC files are those of <bl> 0
    a5 = hist(wd, 5, 1, "ACC")
    rc, lst = run(wd, analysis="DYNAMIC", method="REL", refnode=5, source="FILE8", apdlfile="rel.inp")
    assert rc == 0, lst
    assert V._relmax(V.parse_file(wd / "rel.inp").table("LG_ACX").column(1) / G, a5) < 1e-9
    rc, lst = run(wd, analysis="DYNAMIC", method="ACC", refnode=5, source="FILE8", apdlfile="acc.inp")
    assert rc == 0, lst
    assert V._relmax(V.parse_file(wd / "acc.inp").table("LG_ACX").column(1) / G, a5) > 1e-6


def test_dynamic_reference_consistency_errors(wd):
    rc, lst = run(wd, analysis="DYNAMIC", method="REL", refnode=5)          # RESULTS: RELDISP used the free field
    assert rc == 1 and "relative to the free field" in lst
    rc, lst = run(wd, analysis="DYNAMIC", method="ACC")
    assert rc == 1 and "reference node" in lst


def test_dynamic_ground_file(wd):
    rc, lst = run(wd, analysis="DYNAMIC", groundfile="eq.acc")
    assert rc == 0 and "ground acceleration file = control motion" in lst
    rc, lst = run(wd, analysis="DYNAMIC", groundfile="eq.acc", groundmult=1.2)
    assert rc == 0 and "differs from the control motion" in lst
    ax = V.parse_file(wd / "m_LGD.inp").table("LG_ACX").column(1)
    acc = textfiles.read_history(wd / "eq.acc")[0]
    assert ax[:len(acc)] == pytest.approx(1.2 * G * acc, rel=1e-11)
    rc, lst = run(wd, analysis="DYNAMIC", groundfile="nothere.acc")
    assert rc == 1 and "not found" in lst


def test_dynamic_acc_method_and_rigidity_report(wd):
    rc, lst = run(wd, analysis="DYNAMIC", method="ACC", refnode=5, apdlfile="acc.inp")
    assert rc == 0, lst
    assert "interface nodes vs reference" in lst
    prog = V.parse_file(wd / "acc.inp")
    assert prog.errors == [] and all(v == 0.0 for _, _, v in prog.loads("D"))
    assert {lab for _, lab, _ in prog.loads("D")} == {"UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"}


# ======================================================================================
# numbering of the HOUSE optimizer, 2-D models, runner and batch protocol
# ======================================================================================
def test_optimizer_numbering_is_translated(tmp_path):
    wd = tmp_path / "opt"
    wd.mkdir()
    _build(wd, optimize=True)
    f4 = read_container(wd / "m.N4", "FILE4")
    new_of = dict(zip(f4["x_node_old_id"].tolist(), f4["node_id"].tolist()))
    assert new_of[16] != 16                                         # the optimizer renumbered the model
    md = B.motion_deck(B.FrequencySet.fourier(DT, NFFT, np.arange(1, NFFT // 2)), "eq.acc", [], gravity=G, model="m")
    B.write_deck(wd, "m", md)                                       # MOTION deck: source FILE8 settings
    rc, lst = run(wd, analysis="DYNAMIC", source="FILE8", tables={"checknodes": [16]})
    assert rc == 0, lst
    assert "HOUSE optimizer" in lst
    prog = V.parse_file(wd / "m_LGD.inp")
    assert "LGA_16_UX" in prog.arrays and {n for n, _, _ in prog.loads("D")} == set(MAT_NODES)
    rc, lst = run(wd, analysis="STATIC", data=2, genmass=1, source="FILE8")
    assert rc == 0, lst
    assert set(LG.read_lumped(wd / "m.masl")) == set(STICK_NODES)
    assert {n for n, _, _ in V.parse_file(wd / "m_LGS.inp").loads("F")} == set(STICK_NODES)


def test_two_d_plane_model(tmp_path):
    from sassi.verify.problems.vp_stress import stiff_site_file8
    wd = tmp_path / "w"
    wd.mkdir()
    site = B.uniform_site(depth=20.0, nsub=4, vs=300.0, rho=2.0, beta=0.05, gravity=G)
    fnum = np.arange(1, NFFT // 2)
    fs = B.FrequencySet.fourier(DT, NFFT, fnum)
    hb = B.HouseBuilder(site, dim=1)
    mat = hb.elastic(3e7, 0.2, rho=0.0, beta=0.05)
    n = {(i, j): hb.node(x, 0.0, z) for i, x in enumerate((0.0, 1.0, 2.0)) for j, z in enumerate((0.0, 2.0))}
    g = hb.group(4, "wall")
    for i in range(2):
        hb._element(g, (n[(i, 0)], n[(i + 1, 0)], n[(i + 1, 1)], n[(i, 1)]), etype=1, mat=mat)
    hb.set_interaction([n[(i, 0)] for i in range(3)])
    for i in range(3):
        hb.mass(n[(i, 1)], 1.0)
    B.write_deck(wd, "w", B.site_deck(site, fs, model="w"))
    B.write_deck(wd, "w", hb.deck("w"))
    B.run("HOUSE", wd, "w")
    stiff_site_file8(wd, fnum, 1.0 / (NFFT * DT), 5e3, NFFT, DT, cm=0, model="w")
    V.control_motion(wd, n=300, seed=3)
    B.write_deck(wd, "w", B.motion_deck(fs, "eq.acc", [], gravity=G, model="w"))
    rc, lst = V.run_loadgen(wd, model="w", analysis="STATIC", data=3, genmass=1, source="FILE8")
    assert rc == 0, lst
    assert "2-D (X-Z -> ANSYS X-Y)" in lst
    prog = V.parse_file(wd / "w_LGS.inp")
    assert {lab for _, lab, _ in prog.loads("F")} == {"FX", "FY"}       # SASSI Z -> ANSYS Y
    assert {lab for _, lab, _ in prog.loads("D")} == {"UX", "UY"}
    res = np.loadtxt(wd / "w_LGS_res.txt")
    k = int(round(prog.scalars["LG_TSSI"] / DT))
    assert sum(v for _, lab, v in prog.loads("F") if lab == "FY") == pytest.approx(res[k, 3], abs=1e-9)  # VZ -> FY
    rc, lst = V.run_loadgen(wd, model="w", analysis="DYNAMIC", source="FILE8")
    assert rc == 0, lst
    prog = V.parse_file(wd / "w_LGD.inp")
    assert {c.fields[0] for c in prog.find("ACEL")} == {"%LG_ACX%"}
    assert {lab for _, lab, _ in prog.loads("D")} == {"UX", "UY"}


def test_missing_inputs(tmp_path):
    rc, lst = V.run_loadgen(tmp_path, analysis="STATIC")
    assert rc == 1 and "FILE4" in lst
    assert LG.run_loadgen("nodeck", tmp_path) == 1
    assert "not found" in (tmp_path / "nodeck_LOADGEN.out").read_text()


def test_batch_protocol(wd, monkeypatch):
    d = LG.new_deck()
    d["model"], d["data"], d["genmass"] = "m", 2, 1
    LG.write_deck(wd / "m.lgn", d)
    monkeypatch.chdir(wd)
    assert LG.batch_main(io.StringIO("m\nm.lgn\nm_LOADGEN.out\n")) == 0
    assert (wd / "m_LGS.inp").exists() and "status OK" in (wd / "m_LOADGEN.out").read_text()
    assert LG.batch_main(io.StringIO("m\n")) == 2
