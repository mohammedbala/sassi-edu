"""2D (plane-strain) SSI end to end: SITE -> POINT2 -> HOUSE (<dim> = 1, PLANE, BEAMS in the plane) -> ANALYS
-> MOTION / STRESS, the 2D impedance, input rules and errors, simultaneous X/Z cases, a 2D SYMM line and
the CHECK rules (requirements 4.3 POINT2, 4.6 item 2, D-PNT-01, D-W2-07, D-ANL-12)."""
from __future__ import annotations

import math
import shutil

import numpy as np
import pytest

from sassi.conventions import nodal_result_name
from sassi.core import symmetry as SYM
from sassi.core.flexibility import flexibility_matrix, frequency_row
from sassi.io import decks
from sassi.io.container import read_container
from sassi.prep import Interpreter
from sassi.prep.check import run_check
from sassi.verify import builders as B
from sassi.verify.problems.vp_stress import add_eout, run_stress, stress_deck

SITE_LAYERS = [(2.0, 150.0, 300.0, 1.9, 0.05), (3.0, 250.0, 500.0, 2.0, 0.04), (5.0, 350.0, 700.0, 2.1, 0.03)]
SITE_HS = (600.0, 1200.0, 2.2, 0.02)
FS = B.FrequencySet.fourier(0.01, 1024, [1, 4, 10, 20, 41, 82])
STICK = B.Stick(heights=[3.0, 6.0], masses=[40.0, 30.0], E=3.0e7, A=1.0, I=0.5, beta=0.05)


def _site() -> B.Site:
    return B.layered_site(SITE_LAYERS, SITE_HS)


def _model(**kw) -> B.Model:
    args = dict(half_width=3.0, n_emb=2, ndiv=4, structure="plane", E=3.0e7, nu=0.2, rho=2.4, beta=0.05, stick=STICK)
    args.update(kw)
    return B.embedded_plane_2d(_site(), **args)


@pytest.fixture(scope="module")
def chain2d(tmp_path_factory):
    wd = tmp_path_factory.mktemp("chain2d")
    mdl = _model()
    B.run_soil(wd, "m", _site(), FS, layer=mdl.layer, rad=mdl.rad, dim=1)
    f4 = B.run_house(wd, "m", mdl)
    B.run_analys(wd, "m", FS, save=1)
    return wd, mdl, f4


def test_2d_files_and_listing(chain2d):
    wd, mdl, f4 = chain2d
    assert int(f4.meta["dim"]) == 1 and list(f4.meta["int_dofs"]) == [1, 3]
    assert np.all(np.asarray(f4["int_eq"])[:, 1] == -1)
    f8 = B.read_file8(wd)
    assert int(f8.meta["x_dim"]) == 1 and set(np.asarray(f8["eq_dof"]).tolist()) == {1, 3, 5}
    out = B.listing(wd, "m", "ANALYS")
    assert "2D (plane strain, POINT2)" in out and "interaction DOFs UX, UZ" in out
    assert "out-of-plane DOFs UY, ROTX, ROTZ" in B.listing(wd, "m", "HOUSE")
    # the out-of-plane DOFs of the stick are fixed by the builder: no UY equation off the interaction set
    stick_dofs = {int(d) for n, d in zip(f8["eq_node"], f8["eq_dof"]) if int(n) in mdl["stick"]}
    assert stick_dofs == {1, 3, 5}


def test_2d_low_frequency_rigid_body_and_resonance(chain2d):
    """G-19: at the first SSI frequency (0.098 Hz) every node follows the free field (|ATF_x| ~ 1); the stick
    amplifies the motion at higher frequency."""
    wd, mdl, _ = chain2d
    f8 = B.read_file8(wd)
    H = np.asarray(f8["H"])
    ux = np.asarray(f8["eq_dof"]) == 1
    assert float(np.abs(np.abs(H[0, ux]) - 1.0).max()) < 0.01
    assert "EDU-18" not in B.listing(wd, "m", "ANALYS")
    assert np.abs(B.tf(f8, mdl["top"], 1)).max() > 1.05


def test_2d_impedance_is_the_inverse_of_the_inplane_block(chain2d):
    """COOX holds X_ff = (F_ff in-plane block)^-1 with zero UY rows and columns (requirements 4.6 item 2)."""
    wd, _, f4 = chain2d
    f3 = read_container(wd / "FILE3", "FILE3")
    xyz = np.asarray(f4["x_int_xyz"])
    F = flexibility_matrix(f3, frequency_row(f3, FS.fnum[2]), xyz[:, :2], np.asarray(f4["int_iface"]))
    n = xyz.shape[0]
    ip = np.array([[3 * i, 3 * i + 2] for i in range(n)]).ravel()
    order = int(np.flatnonzero(np.asarray(read_container(wd / "COOXI")["fnum"]) == FS.fnum[2])[0]) + 1
    X = np.asarray(read_container(wd / f"COOX{order:03d}")["X"])
    np.testing.assert_allclose(X[np.ix_(ip, ip)], np.linalg.inv(F[np.ix_(ip, ip)]), rtol=1e-9,
                               atol=1e-12 * np.abs(X).max())
    uy = np.arange(1, 3 * n, 3)
    assert np.all(X[uy] == 0) and np.all(X[:, uy] == 0)


def test_2d_restart_equals_initiation(chain2d, tmp_path):
    wd, _, _ = chain2d
    dst = tmp_path / "m"
    shutil.copytree(wd, dst)
    H0 = np.asarray(B.read_file8(dst)["H"]).copy()
    assert B.run_analys(dst, "m", FS, mode=2) == 0
    np.testing.assert_allclose(np.asarray(B.read_file8(dst)["H"]), H0, rtol=0, atol=1e-10 * np.abs(H0).max())


def test_2d_motion_and_stress(chain2d, tmp_path):
    """MOTION (TF output equal to FILE8; acceleration history) and STRESS (PLANE SXX SZZ TXZ histories) run on
    the 2D FILE8."""
    from sassi.conventions import element_result_name
    from sassi.io import textfiles
    from sassi.verify.problems.vp_analys import read_tf_text
    src, mdl, _ = chain2d
    wd = tmp_path / "m"
    shutil.copytree(src, wd)
    top = mdl["top"]
    f8 = B.read_file8(wd)
    assert B.run_motion(wd, "m", FS, "", [(top, 1, 1, 0, 0, 0, 0, 0)], out=1, cm=0, ang=0.0) == 0
    _, Hm = read_tf_text(wd / nodal_result_name(top, 1, "TFU"))
    np.testing.assert_allclose(Hm, B.tf(f8, top, 1), rtol=0, atol=1e-12 * np.abs(B.tf(f8, top, 1)).max())
    t = np.arange(FS.nft) * FS.delt
    B.write_history(wd / "eq.acc", 0.1 * np.sin(2 * math.pi * 2.0 * t) * np.exp(-0.5 * t), FS.delt)
    assert B.run_motion(wd, "m", FS, "eq.acc", [(top, 1, 1, 1, 0, 0, 0, 1)], cm=0, ang=0.0) == 0
    a, _ = textfiles.read_history(wd / nodal_result_name(top, 1, "ACC"))
    assert np.all(np.isfinite(a)) and np.abs(a).max() > 0
    d = stress_deck(FS.nft, FS.delt, "eq.acc", gravity=9.81)
    add_eout(d, 2, [1, 2], [2, 2, 2])                    # group 2: structural PLANE elements
    rc, out = run_stress(wd, d)
    assert rc == 0, out[-2000:]
    for comp in ("SXX", "SZZ", "TXZ"):
        v, _ = textfiles.read_history(wd / element_result_name("PLANE", 2, 1, comp, "THS"))
        assert np.all(np.isfinite(v)) and np.abs(v).max() > 0


# ---------------------------------------------------------------------------------------------- errors
def test_2d_needs_point2(tmp_path):
    mdl = _model(stick=None)
    B.run_soil(tmp_path, "m", _site(), FS, layer=mdl.layer, rad=mdl.rad, dim=2)        # POINT3 FILE3
    B.run_house(tmp_path, "m", mdl)
    assert B.run_analys(tmp_path, "m", FS, check=False) == 1
    assert "FILE3 holds POINT3 (3D point-load) solutions but the HOUSE model is 2D" in B.listing(tmp_path, "m", "ANALYS")


@pytest.fixture(scope="module")
def base2d(tmp_path_factory):
    wd = tmp_path_factory.mktemp("base2d")
    mdl = _model(stick=None, structure="none")
    B.run_soil(wd, "m", _site(), FS, layer=mdl.layer, rad=mdl.rad, dim=1)
    B.run_house(wd, "m", mdl)
    return wd


@pytest.mark.parametrize("params,msg", [(dict(ang=30.0), "x' must be along +X (0 deg) or -X (180 deg)"),
                                        (dict(impe=1), "available for 3D models only")])
def test_2d_option_errors(base2d, tmp_path, params, msg):
    wd = tmp_path / "m"
    shutil.copytree(base2d, wd)
    assert B.run_analys(wd, "m", FS, check=False, **params) == 1
    assert msg in B.listing(wd, "m", "ANALYS")


def test_2d_refuses_antiplane_input(base2d, tmp_path):
    wd = tmp_path / "m"
    shutil.copytree(base2d, wd)
    B.write_deck(wd, "m", B.site_deck(_site(), FS, wave="SH", model="m"))
    B.run("SITE", wd, "m")
    assert B.run_analys(wd, "m", FS, check=False) == 1
    assert "control motion direction y' (SH / Love waves) is anti-plane" in B.listing(wd, "m", "ANALYS")


def test_2d_angle_180_mirrors_the_input(base2d, tmp_path):
    """x' along -X: the in-plane input changes sign (vertical SV), ATF_x(180) = -ATF_x(0)."""
    wd = tmp_path / "m"
    shutil.copytree(base2d, wd)
    B.run_analys(wd, "m", FS)
    H0 = np.asarray(B.read_file8(wd)["H"]).copy()
    B.run_analys(wd, "m", FS, ang=180.0)
    np.testing.assert_allclose(np.asarray(B.read_file8(wd)["H"]), -H0, rtol=1e-10, atol=1e-12)


def test_2d_simultaneous_cases_are_x_and_z(base2d, tmp_path):
    wd = tmp_path / "m"
    shutil.copytree(base2d, wd)
    B.run_site_xyz(wd, "m", _site(), FS)
    assert B.run_analys(wd, "m", FS, simul=1) == 0
    assert (wd / "FILE8X").exists() and (wd / "FILE8Z").exists() and not (wd / "FILE8Y").exists()
    assert "anti-plane Y case (FILE1Y, SH) is not analysed" in B.listing(wd, "m", "ANALYS")
    hz = B.tf(B.read_file8(wd, "FILE8Z"), int(read_container(wd / "m.N4")["int_node"][-1]), 3)
    assert abs(abs(hz[0]) - 1.0) < 0.01                    # vertical P: rigid-body UZ at low frequency


def test_2d_house_refuses_wave_passage(tmp_path):
    mdl = _model(stick=None)
    B.run_soil(tmp_path, "m", _site(), FS, layer=mdl.layer, rad=mdl.rad, dim=1)
    d = FS.fill(mdl.deck("m"))
    d["wpass"] = 1
    B.write_deck(tmp_path, "m", d)
    assert B.run("HOUSE", tmp_path, "m", check=False) == 1
    assert "wave passage (HOUSE <wpass>) need a 3D model" in B.listing(tmp_path, "m", "HOUSE")


# ------------------------------------------------------------------------------------- 2D SYMM line
def _symm2d(part: str, stype: int) -> B.Model:
    """2D embedded PLANE block (concrete) through two layers with a stick on its axis x = 0; 'half' keeps
    x >= 0 with a SYMM line at x = 0 and half the stick section and masses."""
    site = _site()
    f = 0.5 if part == "half" else 1.0
    hb = B.HouseBuilder(site, dim=1, title=f"2D symmetric block ({part})")
    a, nd = 3.0, 4
    g = np.linspace(-a, a, nd + 1)
    xs = g[g >= -1e-9] if part == "half" else g
    z = site.elevations()[:3]
    ids = {}
    for lev in (2, 1, 0):
        for i, x in enumerate(xs):
            ids[(i, lev)] = hb.node(x, 0.0, z[lev])
    gexc, gstr = hb.group(4, "excavated soil"), hb.group(4, "structure")
    conc = hb.elastic(3.0e7, 0.2, 2.4, 0.05)
    for lev in range(2):
        for i in range(len(xs) - 1):
            quad = [ids[(i, lev + 1)], ids[(i + 1, lev + 1)], ids[(i + 1, lev)], ids[(i, lev)]]
            hb.plane(quad, mat=lev + 1, etype=2, group=gexc)
            hb.plane(quad, mat=conc, etype=1, group=gstr)
    hb.set_interaction(list(ids.values()))
    c = hb.find(0.0, 0.0, z[0])
    spider = [n for n in (hb.find(-1.5, 0.0, z[0]), hb.find(1.5, 0.0, z[0])) if n is not None]
    for n in spider:                                      # rigid links clamp the stick (no hinge at a PLANE node)
        B.rigid_link(hb, c, n, 1.5)
    mat = hb.material(1, 3.0e7, 0.2, 0.0, 0.05, 0.05)
    sec = hb.section(1.0 * f, 0.0, 0.0, 1.0 * f, 0.5 * f, 0.5 * f)
    prev, nodes = c, []
    for hgt, m in ((3.0, 40.0), (6.0, 30.0)):
        n = hb.node(0.0, 0.0, z[0] + hgt)
        hb.beam(prev, n, hb.orientation_node(1.1, 0.0, z[0] + hgt - 1.0), mat, sec, group=hb.group(2, "stick"))
        hb.mass(n, m * f, m * f, m * f, 0.0, 5.0 * f, 0.0)
        nodes.append(n)
        prev = n
    for n in [c] + spider + nodes:
        hb.fix(n, B.OUT_OF_PLANE)
    if part == "half":
        hb.symmetry(1, stype, [hb.find(0.0, 0.0, z[0]), hb.find(0.0, 0.0, z[2])])
    return B.Model(hb, {"top": nodes[-1]}, rad=1.5, layer=2, dim=1)


@pytest.mark.parametrize("wave,stype", [("SV", SYM.ANTISYMMETRIC), ("P", SYM.SYMMETRIC)])
def test_2d_symm_line_half_equals_full(tmp_path, wave, stype):
    from sassi.verify.problems.vp_twod_symm import compare_reduced
    res = {}
    for part in ("full", "half"):
        wd = tmp_path / part
        wd.mkdir()
        mdl = _symm2d(part, stype)
        B.run_soil(wd, "m", _site(), FS, layer=2, rad=mdl.rad, wave=wave, dim=1)
        f4 = B.run_house(wd, "m", mdl)
        B.run_analys(wd, "m", FS)
        res[part] = (B.read_file8(wd), {int(n): np.asarray(p) for n, p in zip(f4["node_id"], f4["node_xyz"])})
    err, ndof = compare_reduced(res["full"], res["half"])
    assert err < 1e-8 and ndof > 10


# ---------------------------------------------------------------------------------------------- CHECK
BASE2D = """
N,1,0,0,-2
N,2,2,0,-2
N,3,2,0,0
N,4,0,0,0
INT,1,4,1,1
L,1,2,1.9,300,150,0.05,0.05
L,2,10,2.0,500,250,0.04,0.04
TOPL,1,1
GROUP,1,PLANE
E,1,1,2,3,4
ETYPE,1,1,1,2
MSET,1,1,1,1
HOUSE,9.81,0,0,1,0,0,0,0,0
"""


def test_check_2d_rules():
    ui = Interpreter()
    ui.run_text(BASE2D + "ANALYS,0,0,0,0,1,0,30,0,0,0,1\nWPASS,1e9,0,1\nHOUSE,9.81,0,0,1,0,0,1,0,0\n")
    rep = run_check(ui.model, modules=["HOUSE", "ANALYS"])
    errs = [m.line() for m in rep.errors()]
    assert any("wave passage not allowed in a 2D model" in e for e in errs)
    assert any("global impedance is available for 3D models only" in e for e in errs)
    assert any("2D model: the coordinate transformation angle must be 0 or 180 deg" in e for e in errs)
    ui = Interpreter()
    ui.run_text(BASE2D + "SITE,0,1,0,20,3,1,1,1,4096,1,1,0.01,1024,1\nANALYS,0,0,0,0,1,0,0,0,0,0,0\n")
    errs = [m.line() for m in run_check(ui.model, modules=["ANALYS"]).errors("ANALYS")]
    assert any("SITE control direction y' (SH / Love) is anti-plane" in e for e in errs)
    ui = Interpreter()
    ui.run_text(BASE2D + "SYMM,1,1,4,3\n")                         # a horizontal 'line' in 2D
    errs = [m.line() for m in run_check(ui.model, modules=["HOUSE"]).errors("HOUSE")]
    assert any("not parallel to the Z axis" in e for e in errs)
    ui = Interpreter()
    ui.run_text(BASE2D + "SYMM,1,1,1,4\n")
    assert not [m for m in run_check(ui.model, modules=["HOUSE"]).errors("HOUSE") if "SYMM" in m.line()]


# ----------------------------------------------------------------------------------- SPRING in 2D
def test_2d_spring_mass_oscillator_peak(tmp_path):
    """A SPRING + mass oscillator on one surface interaction node of a practically rigid 2D site: the ATF peak
    is the hysteretic SDOF value 1/(2 beta sqrt(1 - beta^2)) = 10.0125 at f0 sqrt(1 - 2 beta^2) (R2 section 0)."""
    site = B.layered_site([(10.0, 1.0e5, 2.0e5, 2.0, 0.01)], (1.0e5, 2.0e5, 2.0, 0.01))
    f0, m, beta = 5.0, 1.0, 0.05
    k = (2.0 * math.pi * f0) ** 2 * m
    hb = B.HouseBuilder(site, dim=1, title="2D SDOF")
    b = hb.add_node(0.0, 0.0, 0.0, fix=(0, 1, 0, 1, 1, 1))
    mn = hb.add_node(0.0, 0.0, 1.0, fix=(0, 1, 1, 1, 1, 1))
    hb.set_interaction([b])
    hb.spring(b, mn, [k, 0, 0, 0, 0, 0], beta)
    hb.mass(mn, m, m, m)
    fs = B.FrequencySet.harmonic(0.001, [4986, 4987, 4988, 4989])
    B.run_soil(tmp_path, "m", site, fs, layer=0, rad=1.0, dim=1)
    B.run_house(tmp_path, "m", B.Model(hb, {}, rad=1.0, dim=1))
    B.run_analys(tmp_path, "m", fs)
    H = B.tf(B.read_file8(tmp_path), mn, 1)
    assert np.abs(H).max() == pytest.approx(1.0 / (2.0 * beta * math.sqrt(1.0 - beta ** 2)), rel=2e-5)


def test_2d_strip_compliance_is_reciprocal(tmp_path):
    """Vibration of a rigid 2D strip: the compliance of FORCE unit loads Fx / My is reciprocal (C_x,phi =
    C_phi,x), and rocking and sliding couple (welded contact)."""
    site = B.layered_site([(1.0, 100.0, 187.08, 2.0, 0.02)] * 4, (100.0, 187.08, 2.0, 0.02), nl=0)
    fs = B.FrequencySet.harmonic(0.5, [1, 4])
    mdl = B.rigid_strip_2d(site, half_width=1.0, ndiv=8)
    B.run_soil(tmp_path, "m", site, fs, layer=0, rad=mdl.rad, mode2=False, dim=1)
    B.run_house(tmp_path, "m", mdl)
    c = mdl["centre"]
    B.run_force(tmp_path, "m", [(c, 1, 1.0, 0.0)], fs, copy_to="FILE9001")
    B.run_force(tmp_path, "m", [(c, 5, 1.0, 0.0)], fs, copy_to="FILE9002")
    B.run_analys(tmp_path, "m", fs, type=1, simul=2)
    f1, f2 = B.read_file8(tmp_path, "FILE8001"), B.read_file8(tmp_path, "FILE8002")
    cxp, cpx = B.tf(f2, c, 1), B.tf(f1, c, 5)
    np.testing.assert_allclose(cxp, cpx, rtol=1e-6)
    assert np.all(np.abs(cxp) > 1e-3 * np.abs(B.tf(f1, c, 1)))


# --------------------------------------------------------------------------- commands -> AFWRITE -> RUN
PRE2D = """TIT,2D embedded block with a stick (half model, SYMM line x = 0)
L,1,2.0,19.0,300,150,0.05,0.05
L,2,3.0,20.0,500,250,0.04,0.04
L,3,1.0,21.0,700,350,0.03,0.03
TOPL,1,2
FREQ,1,4,10,20,41
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.01,1024,1
WAVE,2,1,1,1,0
N,1,0,0,-5
N,3,3,0,-5
FILL,1,3
N,4,0,0,-2
N,6,3,0,-2
FILL,4,6
N,7,0,0,0
N,9,3,0,0
FILL,7,9
N,10,0,0,3
N,11,0,0,6
N,12,1,0,1.5
N,13,1,0,4.5
N,14,0.75,0,1
D,12,14,1,1,ALL
GROUP,1,PLANE
E,1,1,2,5,4
E,2,2,3,6,5
E,3,4,5,8,7
E,4,5,6,9,8
ETYPE,1,4,1,2
MSET,1,2,1,2
MSET,3,4,1,1
M,1,3.0E7,0.2,24.0,0.05,0.05,1
M,2,3.0E7,0.2,0.0,0.05,0.05,1
R,1,0.5,0,0,0.5,0.25,0.25
R,3,0.4,0.33,0.33,0.01,0.005,0.005
GROUP,2,BEAMS
MACT,1
RACT,3
E,1,7,8,12
E,2,8,9,12
E,3,3,6,12
E,4,6,9,12
E,5,1,2,12
E,6,2,3,12
GROUP,3,BEAMS
MACT,2
RACT,1
E,1,7,10,12
E,2,10,11,13
D,1,11,1,1,UY
D,1,11,1,1,ROTX
D,1,11,1,1,ROTZ
MT,10,100,0,100
MT,11,75,0,75
INT,1,9,1,1
SYMM,1,1,7,1
POINT,0,2,1.5
HOUSE,9.81,0,0,1,0,0,0,0,0
ANALYS,0,0,0,0,1,0,0,0,0,0,0
AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0
CHECK
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
"""


def test_2d_half_model_from_commands(tmp_path):
    """The command path: a 2D half model (PLANE excavation, BEAMS basement and stick, SYMM line x = 0) is
    checked, written by AFWRITE (POINT <dim> = 1 from HOUSE <dim>, the SYMM rows in the HOUSE deck) and run
    SITE -> POINT2 -> HOUSE -> ANALYS without errors."""
    from sassi.prep import Kind
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(f"MDL,t2d,{tmp_path}\n" + PRE2D)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    assert int(decks.read(tmp_path / "t2d.poi", "POINT")["dim"]) == 1
    assert decks.read(tmp_path / "t2d.hou", "HOUSE").rows("symm")[0]["type"] == 1
    f8 = read_container(tmp_path / "FILE8", "FILE8")
    assert int(f8.meta["x_dim"]) == 1 and f8.meta["x_symm"] == [[1.0, 1.0, 0.0, 0.0]]
    out = (tmp_path / "t2d_ANALYS.out").read_text()
    assert "3 interaction translations constrained" in out
