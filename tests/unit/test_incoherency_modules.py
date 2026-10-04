"""Incoherency, wave passage and multiple excitation through the modules and the interpreter: HOUSE
(FILE77, I N C O table, checks, ME zones, stale-file handling), ANALYS (FFL / FFM, X/Y/Z rows, stochastic
cases FILE8001 ..., errors), MOTION (single-mode anchor, SRSS modes), the commands (INCOH, WPASS, ME,
HOUSE, HOUSEX, ANALYSX, BUILDFILE77), CHECK and AFWRITE (requirements 2.4, 4.4 item 6, 4.6 item 5,
4.8 items 5-6)."""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from sassi.core import incoherency as I
from sassi.io.container import read_container
from sassi.prep import Interpreter, Kind
from sassi.prep.check import run_check
from sassi.prep.options import eduopt
from sassi.verify import builders as B
from sassi.verify.problems.vp_incoherency import house_incoherent
from tests.unit.test_check_catalogue import BASE, _files

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")

SITE = B.uniform_site(12.0, 3, vs=250.0, rho=2.0, nu=1.0 / 3.0, beta=0.03)
FS = B.FrequencySet.fourier(0.01, 1024, [5, 20, 60])
LW = dict(coh=1, cohf=1, alpha=250.0, gammax=0.3, gammay=0.3, gammaz=0.3)


def _mat():
    return B.surface_rigid_mat(SITE, half_width=4.0, ndiv=2)


@pytest.fixture(scope="module")
def soil(tmp_path_factory):
    """SITE (SV x', SH y', P z') and POINT once for the module."""
    wd = tmp_path_factory.mktemp("incoh_soil")
    B.run_site_xyz(wd, "m", SITE, FS)
    B.write_deck(wd, "m", B.point_deck(FS, 0, _mat().rad, model="m"))
    B.run("POINT", wd, "m")
    return wd


def _wd(soil, tmp_path) -> Path:
    for nm in ("FILE1", "FILE1X", "FILE1Y", "FILE1Z", "FILE3", "m.sit"):
        shutil.copyfile(soil / nm, tmp_path / nm)
    return tmp_path


def _house_rc(wd, mdl, **params):
    d = FS.fill(mdl.deck("m"))
    for k, v in params.items():
        d[k] = v
    B.write_deck(wd, "m", d)
    return B.run("HOUSE", wd, "m", check=False), B.listing(wd, "m", "HOUSE")


# ======================================================================================
# HOUSE
# ======================================================================================
def test_house_writes_deterministic_file77_and_the_inco_table(soil, tmp_path):
    wd = _wd(soil, tmp_path)
    mdl = _mat()
    house_incoherent(wd, mdl, FS, ipr=1, **LW)
    c = read_container(wd / "FILE77", "FILE77")
    f4 = read_container(wd / "m.N4", "FILE4")
    assert np.array_equal(c["int_node"], f4["int_node"]) and np.asarray(c["s"]).shape == (3, 3, 9)
    assert c.meta["method"] == "AS" and c.meta["stochastic"] == 0 and c.meta["int_hash"]
    assert float(np.asarray(c["x_trace_error"]).max()) < 1e-12
    out = B.listing(wd, "m", "HOUSE")
    assert "I N C O" in out and "spectral factorisation summary" in out and "FILE77 (AS factors of 9" in out
    # X and Y share the Luco-Wong matrix, so they share the decomposition and the AS factors
    assert np.array_equal(np.asarray(c["s"])[:, 0], np.asarray(c["s"])[:, 1])


def test_house_stochastic_samples_and_stale_files(soil, tmp_path):
    wd = _wd(soil, tmp_path)
    mdl = _mat()
    house_incoherent(wd, mdl, FS, **LW)
    house_incoherent(wd, mdl, FS, hseed=3, vseed=4, randphz=180.0, nsim=3, **LW)
    assert all((wd / I.file77_name(k)).exists() for k in (1, 2, 3))
    assert not (wd / "FILE77").exists() and (wd / "FILE77.prev").exists()
    c2 = read_container(wd / "FILE77002", "FILE77")
    assert c2.meta["method"] == "SS" and c2.meta["sim"] == 2 and c2.meta["nsim"] == 3
    house_incoherent(wd, mdl, FS, **LW)
    assert (wd / "FILE77").exists() and not (wd / "FILE77001").exists() and (wd / "FILE77003.prev").exists()


@pytest.mark.parametrize("params,msg", [
    (dict(coh=1, cohf=2, wpass=1), "coefficients not available"),
    (dict(coh=1, cohf=5), "applied only with wave passage"),
    (dict(coh=1, cohf=7, wpass=1, alpha=0.5), "needs the files"),
    (dict(coh=1, cohf=5, wpass=1, alpha=2.0), "Error 114"),
    (dict(me=1, wpass=1), "Error 115"),
    (dict(coh=1, cohf=1, alpha=250.0, hseed=1, vseed=2, randphz=180.0, nsim=51), "1..50"),
    (dict(coh=1, cohf=1, alpha=250.0, gammax=0.05), "Error 58"),
    (dict(coh=1, cohf=1, alpha=250.0, nmodes=-12), "mode 12 does not exist"),
])
def test_house_incoherency_errors(soil, tmp_path, params, msg):
    wd = _wd(soil, tmp_path)
    rc, out = _house_rc(wd, _mat(), **params)
    assert rc == 1 and msg in out and not (wd / "FILE77").exists(), out[-2000:]


def test_house_bottom_up_order_is_an_error_when_incoherent(soil, tmp_path):
    wd = _wd(soil, tmp_path)
    hb = B.HouseBuilder(SITE, title="two levels")
    a = hb.node(0.0, 0.0, -4.0)
    b = hb.node(0.0, 0.0, 0.0)
    c = hb.node(4.0, 0.0, 0.0)
    B.rigid_link(hb, a, b, 1.0)
    B.rigid_link(hb, b, c, 1.0)
    hb.set_interaction([a, b, c])
    d = FS.fill(hb.deck("m"))
    d.tables["interaction"].rows = d.tables["interaction"].rows[::-1]       # top-down: not bottom-up
    B.write_deck(wd, "m", d)
    assert B.run("HOUSE", wd, "m", check=False) == 0                         # coherent: a warning
    assert "EDU-21" in B.listing(wd, "m", "HOUSE")
    for k, v in LW.items():
        d[k] = v
    B.write_deck(wd, "m", d)
    assert B.run("HOUSE", wd, "m", check=False) == 1
    assert "ERROR: EDU-21" in B.listing(wd, "m", "HOUSE")


def test_house_multiple_excitation_factors(soil, tmp_path):
    wd = _wd(soil, tmp_path)
    mdl = _mat()
    d = FS.fill(mdl.deck("m"))
    d["me"], d["wpass"], d["cmplxspec"], d["appv"] = 1, 1, 1, 1.0e12
    ids = sorted(mdl["mat"])
    d.table("me").append([1, ids[0], ids[2]])
    for k in (1, 2, 3):
        d.table("amp").append([1, k, 0.5 * k, 1.0])
    B.write_deck(wd, "m", d)
    B.run("HOUSE", wd, "m")
    c = read_container(wd / "FILE77", "FILE77")
    s = np.asarray(c["s"])
    pos = {int(n): k for k, n in enumerate(c["int_node"])}
    zone = [pos[n] for n in ids[:3]]
    rest = [pos[n] for n in ids[3:]]
    for q, k in enumerate((1, 2, 3)):
        assert np.allclose(s[q, :, zone], 0.5 * k + 1.0j, atol=1e-6)
        assert np.allclose(s[q, :, rest], 1.0, atol=1e-6)
    assert c.meta["method"] == "COHERENT" and c.meta["me_zones"] == [[1, ids[0], ids[2]]]


def test_house_optimizer_keeps_factors_on_the_renumbered_nodes(soil, tmp_path):
    wd = _wd(soil, tmp_path)
    mdl = _mat()
    house_incoherent(wd, mdl, FS, **LW)
    ref = read_container(wd / "FILE77", "FILE77")
    house_incoherent(wd, mdl, FS, optimize=1, **LW)
    new = read_container(wd / "FILE77", "FILE77")
    f4 = read_container(wd / "m.N4", "FILE4")
    assert np.array_equal(new["int_node"], f4["int_node"])
    pairs = dict(tuple(int(v) for v in ln.split()) for ln in (wd / "m.map").read_text().splitlines() if ln.strip()
                 and not ln.startswith("#"))
    old_pos = {int(n): k for k, n in enumerate(ref["int_node"])}
    back = {v: k for k, v in pairs.items()}
    for k, n in enumerate(new["int_node"]):
        assert np.allclose(np.asarray(new["s"])[:, :, k], np.asarray(ref["s"])[:, :, old_pos[back[int(n)]]])


def test_house_data_check_writes_no_file77(soil, tmp_path):
    wd = _wd(soil, tmp_path)
    rc, out = _house_rc(wd, _mat(), opmode=1, **LW)
    assert rc == 0 and "Incoherency, wave passage and multiple excitation" in out and not (wd / "FILE77").exists()


# ======================================================================================
# ANALYS
# ======================================================================================
def _factors77(wd, rows):
    """Replace FILE77 by uniform factors ``rows`` (X, Y, Z) at every node and frequency."""
    c = read_container(wd / "FILE77", "FILE77")
    s = np.empty_like(np.asarray(c["s"]))
    for d_, v in enumerate(rows):
        s[:, d_, :] = v
    I.write_file77(wd / "FILE77", c["fnum"], c["freq"], c["int_node"], c["x_int_xyz"], s,
                   {"df": c.meta["df"], "method": "AS", "stochastic": 0})


@pytest.mark.parametrize("ffm", [0, 1])
def test_analys_applies_the_direction_rows(soil, tmp_path, ffm):
    wd = _wd(soil, tmp_path)
    house_incoherent(wd, _mat(), FS, **LW)
    B.run_analys(wd, "m", FS, simul=1)
    coh = {c: np.asarray(B.read_file8(wd, f"FILE8{c}")["H"]) for c in "XYZ"}
    _factors77(wd, (1.0, 2.0, 3.0 - 1.0j))
    B.run_analys(wd, "m", FS, simul=1, coh=1, ffm=ffm)
    for c, v in zip("XYZ", (1.0, 2.0, 3.0 - 1.0j)):
        f8 = B.read_file8(wd, f"FILE8{c}")
        assert np.allclose(np.asarray(f8["H"]), v * coh[c], rtol=1e-12, atol=1e-14)
        assert f8.meta["file77"] == "FILE77" and f8.meta["x_file77_direction"] == c and f8.meta["ffm"] == ffm
    out = B.listing(wd, "m", "ANALYS")
    assert ("FFM, free-field motion" if ffm else "FFL, free-field load") in out


def test_analys_stochastic_cases(soil, tmp_path):
    wd = _wd(soil, tmp_path)
    house_incoherent(wd, _mat(), FS, hseed=3, vseed=4, randphz=180.0, nsim=2, **LW)
    B.run_analys(wd, "m", FS, simul=2, coh=1)
    names = [f"FILE8{k:03d}" for k in range(1, 7)]
    assert all((wd / n).exists() for n in names)
    for k, n in enumerate(names):
        f8 = B.read_file8(wd, n)
        sim, d_ = k // 3 + 1, k % 3
        assert f8.meta["sim"] == sim and f8.meta["cm"] == d_ and f8.meta["file77"] == f"FILE77{sim:03d}"
        assert f8.meta["case"] == f"S{sim:03d}{'XYZ'[d_]}" and f8.meta["x_incoherency"] == "SS"
    assert not np.allclose(B.read_file8(wd, "FILE8001")["H"], B.read_file8(wd, "FILE8004")["H"])
    # <simul> = 1 with stochastic samples only (no FILE77): the first simulation, FILE8001 .. FILE8003
    h1 = np.asarray(B.read_file8(wd, "FILE8001")["H"])
    for n in names[3:]:
        (wd / n).unlink()
    B.run_analys(wd, "m", FS, simul=1, coh=1)
    assert np.array_equal(np.asarray(B.read_file8(wd, "FILE8001")["H"]), h1)
    assert all((wd / n).exists() for n in names[:3]) and not any((wd / n).exists() for n in names[3:])
    assert B.run_analys(wd, "m", FS, simul=3, coh=1, check=False) == 1
    assert "FILE77003 missing" in B.listing(wd, "m", "ANALYS")
    assert B.run_analys(wd, "m", FS, simul=0, coh=1, check=False) == 1
    assert "stochastic samples FILE77001" in B.listing(wd, "m", "ANALYS")


def test_analys_refuses_incompatible_file77(soil, tmp_path):
    wd = _wd(soil, tmp_path)
    house_incoherent(wd, _mat(), FS.subset([5, 20]), **LW)            # FILE77 without frequency number 60
    B.run_analys(wd, "m", FS.subset([5, 20]), coh=1)
    assert B.run_analys(wd, "m", FS, coh=1, check=False) == 1
    assert "Frequency 60 not in FILE77" in B.listing(wd, "m", "ANALYS")
    c = read_container(wd / "FILE77", "FILE77")
    I.write_file77(wd / "FILE77", c["fnum"], c["freq"], np.asarray(c["int_node"]) + 1000, c["x_int_xyz"],
                   c["s"], {"df": c.meta["df"], "method": "AS", "stochastic": 0})
    assert B.run_analys(wd, "m", FS.subset([5, 20]), coh=1, check=False) == 1
    assert "has no factors for interaction nodes" in B.listing(wd, "m", "ANALYS")


def test_analys_single_mode_skips_the_rigid_body_check(soil, tmp_path):
    wd = _wd(soil, tmp_path)
    house_incoherent(wd, _mat(), FS, nmodes=-2, **LW)
    B.run_analys(wd, "m", FS, coh=1)
    out = B.listing(wd, "m", "ANALYS")
    assert "single incoherent mode 2" in out and "EDU-18" not in out
    f8 = B.read_file8(wd)
    assert f8.meta["x_incoherency"] == "SINGLE" and f8.meta["x_mode"] == 2


# ======================================================================================
# MOTION
# ======================================================================================
def _tfi(path):
    rows = [ln.split() for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    a = np.array([[float(v) for v in r] for r in rows])
    return a[:, 0], a[:, 1] * np.exp(1j * a[:, 2])


def test_motion_single_mode_anchor_and_srss_modes(soil, tmp_path):
    from sassi.conventions import nodal_result_name
    wd = _wd(soil, tmp_path)
    mdl = _mat()
    c = mdl["centre"]
    for k in (1, 2):
        house_incoherent(wd, mdl, FS, nmodes=-k, **LW)
        B.run_analys(wd, "m", FS, coh=1)
        shutil.copyfile(wd / "FILE8", wd / f"FILE8_M{k}")
    B.run_motion(wd, "m", FS, "", [(c, 1, 1, 0, 0, 0, 0, 0)], out=1, file8="FILE8_M2", interp=6)
    f, H = _tfi(wd / nodal_result_name(c, 1, "TFI"))
    assert f[0] == 0.0 and abs(H[0]) == 0.0                     # mode 2: H(0) = 0 (D-INC-03)
    assert "single mode 2" in B.listing(wd, "m", "MOTION")
    B.run_motion(wd, "m", FS, "", [(c, 1, 1, 0, 0, 0, 0, 0)], out=1, file8="FILE8_M1", interp=6)
    f, H = _tfi(wd / nodal_result_name(c, 1, "TFI"))
    assert abs(H[0]) == pytest.approx(1.0)                      # mode 1 carries the coherent motion
    (wd / "SRSSTF.txt").write_text("2 0\nFILE8_M2\nFILE8_M1\n")
    B.run_motion(wd, "m", FS, "", [(c, 1, 1, 0, 0, 0, 0, 0)], out=1, srss=1, interp=6)
    out = B.listing(wd, "m", "MOTION")
    assert "incoherent modes: 2, 1" in out
    f, H = _tfi(wd / nodal_result_name(c, 1, "TFI"))
    assert abs(H[0]) == pytest.approx(1.0)                      # the anchor follows mode 1, listed second


# ======================================================================================
# interpreter: commands, CHECK, AFWRITE
# ======================================================================================
def _ui(tmp_path, extra=""):
    ui = Interpreter(cwd=tmp_path)
    mdir = tmp_path / "m"
    mdir.mkdir(exist_ok=True)
    _files(mdir)
    ui.run_text(BASE.format(dir=mdir) + "\n" + extra)
    return ui, mdir


def test_commands_no_longer_report_tier_p1(tmp_path):
    ui, _ = _ui(tmp_path)
    n0 = len(ui.sink.texts(Kind.WARNING))
    for line in ("HOUSE,32.2,0,0,2,0,1,1,1,1", "INCOH,0.2,0.2,0.3,0.5,1,1,0,1,5,7,180", "WPASS,1e9,30,5",
                 "ME,1,1,4,0,0,0", "AMP,1,1,0.5,1,0.5,1,0.5,1,0.5", "HOUSEX,0,0,20", "ANALYSX,1"):
        assert ui.execute(line)
    new = ui.sink.texts(Kind.WARNING)[n0:]
    assert not any("not available in this build" in w or "tier P1" in w for w in new), new
    assert ui.model.options.record("WPASS").cohf == 5 and ui.model.options.record("HOUSEX").nsim == 20
    ui.execute("WPASS,1e9,0,4")
    assert any("coefficients not available" in w for w in ui.sink.texts(Kind.WARNING))


def test_write_inp_round_trip_of_the_incoherency_records(tmp_path):
    ui, mdir = _ui(tmp_path, "HOUSE,32.2,0,0,2,0,1,1,1,1\nINCOH,0.2,0.2,0.3,0.5,1,1,-2,1,0,0,0\n"
                             "WPASS,1500,30,5\nME,1,1,2,0,0,0\nME,2,3,4,0,0,0\nAMP,1,1,0,2,0.5,1,-1,1,1\n"
                             "AMP,2,1,1,1,1,1,1,1,1\nHOUSEX,0,1,1\nANALYSX,1\nEDUOPT,INCOHSIGN,RAW\n"
                             "EDUOPT,INCOHMERGE,1")
    assert ui.execute(f"WRITE,{tmp_path / 'rt.pre'}")
    ui2 = Interpreter(cwd=tmp_path)
    assert ui2.execute(f"INP,{tmp_path / 'rt.pre'}")
    for name in ("HOUSE", "INCOH", "WPASS", "HOUSEX", "ANALYSX"):
        assert ui2.model.options.record(name).to_tokens() == ui.model.options.record(name).to_tokens(), name
    assert ui2.model.amp == ui.model.amp
    for no in (1, 2):
        assert ui2.model.options.entry("ME", no).to_tokens() == ui.model.options.entry("ME", no).to_tokens()
    assert eduopt(ui2.model, "INCOHSIGN") == "RAW" and eduopt(ui2.model, "INCOHMERGE") == "1"


def test_buildfile77_command(tmp_path):
    ui, mdir = _ui(tmp_path)
    xyz = np.array([[0, 0, -4.0], [1, 0, -4.0], [0, 0, 0.0]])
    m = {"df": 0.5, "method": "AS", "stochastic": 0}
    I.write_file77(mdir / "F77_L1", [3, 5], [1.5, 2.5], [1, 2], xyz[:2], np.ones((2, 3, 2)), m)
    I.write_file77(mdir / "F77_L2", [3, 5], [1.5, 2.5], [3], xyz[2:], 2 * np.ones((2, 3, 1)), m)
    assert ui.execute("BUILDFILE77,FILE77,F77_L1,F77_L2")
    c = I.read_file77(mdir / "FILE77")
    assert c["int_node"].tolist() == [1, 2, 3] and np.allclose(np.asarray(c["s"])[0, 0], [1, 1, 2])
    assert c.meta["method"] == "BUILT" and c.meta["built_from"] == ["F77_L1", "F77_L2"]
    assert not ui.execute("BUILDFILE77,FILE77")
    assert not ui.execute("BUILDFILE77,FILE77,F77_L1,F77_L1")


@pytest.mark.parametrize("extra,kind,number,module", [
    ("HOUSE,32.2,0,0,2,0,1,1,0,0\nWPASS,1e9,0,2", "Error", "EDU-26", "HOUSE"),
    ("HOUSE,32.2,0,0,2,0,1,1,0,0\nWPASS,1e9,0,5\nINCOH,0.1,0.1,0.2,1.5", "Error", 114, "HOUSE"),
    ("HOUSE,32.2,0,0,2,0,1,1,0,0\nWPASS,1e9,0,7\nINCOH,0.1,0.1,0.2,0.5", "Error", "EDU-26", "HOUSE"),
    ("HOUSE,32.2,0,0,2,0,1,0,0,0\nINCOH,0.1,0.1,0.2,300,1,0,0,0,5,7,180\nHOUSEX,0,0,3\n"
     "ANALYS,0,0,0,0,1,0,0,0,0,0,0,2", "Error", "EDU-26", "ANALYS"),
    ("ANALYS,0,0,0,0,1,0,0,0,0,0,0,2", "Error", "EDU-26", "ANALYS"),
    ("ANALYSX,1", "Warning", "EDU-28", "ANALYS"),
])
def test_check_incoherency_rules(tmp_path, extra, kind, number, module):
    ui, mdir = _ui(tmp_path, extra)
    rep = run_check(ui.model, dirs=[mdir, tmp_path])
    assert rep.has(kind, number, module), rep.format()


def test_check_complex_sar_is_checked_on_the_modulus(tmp_path):
    ui, mdir = _ui(tmp_path, "HOUSE,32.2,0,0,2,0,0,1,1,1\nME,1,1,4\nAMP,1,-1,0.5,0,2,-3,0,1,1")
    rep = run_check(ui.model, dirs=[mdir, tmp_path])
    assert not rep.has("Error", 118, "HOUSE") and not rep.has("Error", 119, "HOUSE"), rep.format()
    ui.execute("AMP,1,0")
    ui.execute("AMP,1,11,0,1,1,1,1,1,1")
    assert run_check(ui.model, dirs=[mdir, tmp_path]).has("Error", 118, "HOUSE")


def test_afwrite_writes_the_house_option_file(tmp_path):
    ui, mdir = _ui(tmp_path, "HOUSE,32.2,0,0,2,0,1,0,0,0\nINCOH,0.1,0.1,0.2,300\nEDUOPT,INCOHSIGN,RAW\n"
                             "EDUOPT,INCOHMERGE,1")
    assert ui.execute("AFWRITE")
    assert (mdir / "base.hou").exists()
    assert I.read_options(mdir) == {"INCOHSIGN": "RAW", "INCOHMERGE": "1"}
