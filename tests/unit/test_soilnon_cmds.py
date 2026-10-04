"""SOIL-NON through the interpreter and the SOIL module: NLSOIL / NLSLAYER / DELNLS commands, WRITE -> INP,
AFWRITE (<model>.nls side file), the SOIL-NON branch of the SOIL module (errors, data check, outputs,
FILE88), requirements 3.4.O and spec 05a section 6.6."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sassi.core import shake as SH
from sassi.core import soilnon as SN
from sassi.io import decks, textfiles
from sassi.modules import soil
from sassi.modules.base import run_module
from sassi.prep import Interpreter, Kind

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")

LIB = Path(soil.__file__).resolve().parents[1] / "data" / "dynp_library.pre"


def _msgs(ui, kind):
    return ui.sink.texts(kind)


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


# ---------------------------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------------------------
def test_nlsoil_record_and_write_inp_round_trip(tmp_path, ui):
    ui.execute("NLSOIL,1,2,1e-5,,20,1,3,0.1,0.002")
    for line in ("NLSLAYER,1,1", "NLSLAYER,2,0,1.2,0.85,0.05,0", "NLSLAYER,3,0,0"):
        ui.execute(line)
    rec = ui.model.options.record("NLSOIL")
    assert rec.values == ["1", "2", "1e-5", None, "20", "1", "3", "0.1", "0.002"]
    assert rec.opt == 1 and rec.nsub == 2 and rec.dispconv == 1e-5 and rec.forceconv == 0.0 and rec.damptype == 3
    assert ui.model.options.entry("NLSLAYER", 2).refstrain == 0.05
    assert any("Joyner" in w for w in _msgs(ui, Kind.WARNING))          # BedInt 1 is an extension of V3
    ui.execute(f"MDL,rt,{tmp_path}")
    ui.execute("WRITE,rt.pre")
    text = (tmp_path / "rt.pre").read_text()
    assert "NLSOIL,1,2,1e-5,,20,1,3,0.1,0.002" in text and "NLSLAYER,2,0,1.2,0.85,0.05,0" in text
    ui2 = Interpreter(cwd=tmp_path)
    ui2.execute("INP,rt.pre")
    assert ui2.model.options.record("NLSOIL").values == rec.values
    assert [k for k, _ in ui2.model.options.entries("NLSLAYER")] == [1, 2, 3]
    assert ui2.model.options.entry("NLSLAYER", 2).values == ui.model.options.entry("NLSLAYER", 2).values


@pytest.mark.parametrize("line,text", [("NLSOIL,1,0,0,0,0,2", "BedInt"), ("NLSOIL,1,0,0,0,0,0,7", "NLDampType"),
                                       ("NLSOIL,1,0,0,0,0,0,3,-1", "multipliers"), ("NLSOIL,x", "not a number"),
                                       ("NLSOIL,1,1.5", "integer"), ("NLSOIL,1,-2", "NSTimeSunInc")])
def test_nlsoil_rejects_invalid_values(ui, line, text):
    ui.execute("NLSOIL,1")
    assert not ui.execute(line)
    assert any(text in e for e in _msgs(ui, Kind.ERROR))
    assert ui.model.options.record("NLSOIL").values == ["1"]              # previous record kept


def test_nlslayer_storage_and_messages(ui):
    ui.execute("NLSLAYER,3,1,0.5")
    assert ui.model.options.entry("NLSLAYER", 3).values == ["3", "1", "0.5"]
    assert any("ignored" in t for t in _msgs(ui, Kind.INFO))
    ui.execute("NLSLAYER,2,0,0")
    assert any("linear elastic" in w for w in _msgs(ui, Kind.WARNING))
    ui.execute("NLSLAYER,4,0,1.2,1.3,0.05")
    assert any("S = 1.3 > 1" in w for w in _msgs(ui, Kind.WARNING))
    ui.execute("NLSLAYER,5,0,1.2,0,0.05")
    assert any("needs Beta, S exponent and Reference Strain" in w for w in _msgs(ui, Kind.WARNING))
    for bad in ("NLSLAYER,0,1", "NLSLAYER,6,2", "NLSLAYER,7,0,-1", "NLSLAYER"):
        assert not ui.execute(bad)
    assert [k for k, _ in ui.model.options.entries("NLSLAYER")] == [2, 3, 4, 5]


def test_delnls_range_stride_and_linear_data_untouched(ui):
    ui.execute("L,1,5.0,18.0,400,200,0.02,0.02")
    ui.execute("SPRO,1,1,Sand")
    for k in range(1, 8):
        ui.execute(f"NLSLAYER,{k},1")
    ui.sink.clear()
    ui.execute("DELNLS,2,6,2")
    assert [k for k, _ in ui.model.options.entries("NLSLAYER")] == [1, 3, 5, 7]
    assert not _msgs(ui, Kind.WARNING)
    ui.execute("DELNLS,5,7")
    assert any("default stride" in w for w in _msgs(ui, Kind.WARNING))
    assert [k for k, _ in ui.model.options.entries("NLSLAYER")] == [1, 3]
    ui.execute("DELNLS,3,3,0")
    assert [k for k, _ in ui.model.options.entries("NLSLAYER")] == [1]
    ui.execute("DELNLS,20")
    assert any("no NLSLAYER" in t for t in _msgs(ui, Kind.INFO))
    assert 1 in ui.model.layers and ui.model.options.entry("SPRO", 1) is not None
    assert not ui.execute("DELNLS")


def _model_text(wd: Path, nls: str = "", n: int = 400) -> str:
    rng = np.random.default_rng(1)
    a = 0.1 * rng.standard_normal(n) * np.hanning(n)
    (wd / "in.acc").write_text("\n".join(f"{v:.6e}" for v in a) + "\n")
    return f"""MDL,m,{wd}
L,1,4.0,18.0,400,200,0.02,0.02
L,2,1.0,22.0,2000,1000,0.01,0.01
TOPL,1,1
SITE,0,1,0,20,2,1,0,1,512,1,0,0.01,1024,1
INP,{LIB}
SPRO,1,1,Sand
SPRO,2,1,Sand
SPRO,3,2
SOIL,{n},9.81,0,1,1,8,0.65,1,0
SOILX,0,0,0.15,3
THFILE,in.acc
SACC,1,2,0
DAMP,0.05
{nls}AOPT,0,1,0,0,0,0,0,0,0,0,0,0,0,0
AFWRITE
"""


def test_afwrite_writes_the_nls_file_and_soil_runs_soil_non(tmp_path, ui):
    nls = "NLSOIL,1,0,0,0,0,1,1\nNLSLAYER,1,1\nNLSLAYER,2,0,1.0,0.9,0.05\nNLSLAYER,3,1\n"
    ui.run_text(_model_text(tmp_path, nls), name="m.pre")
    assert not _msgs(ui, Kind.ERROR), _msgs(ui, Kind.ERROR)
    data = SN.read_nls(tmp_path / "m.nls")
    assert data.options.opt == 1 and data.options.bedint == 1 and sorted(data.layers) == [1, 2, 3]
    assert data.model_hash and data.model_hash in (tmp_path / "m.soi").read_text()
    ui.execute("RUNSOIL")
    out = (tmp_path / "m_SOIL.out").read_text()
    assert "SOIL-NON" in out and "finished with status OK" in out
    assert (tmp_path / "ACC001.TH").exists() and (tmp_path / "FILE88").exists()
    # without NLSOIL a stale .nls is retired (SOIL is SOIL-EQL again)
    ui2 = Interpreter(cwd=tmp_path)
    ui2.run_text(_model_text(tmp_path, ""), name="m.pre")
    assert not (tmp_path / "m.nls").exists() and (tmp_path / "m.nls.bak").exists()
    ui2.execute("RUNSOIL")
    out = (tmp_path / "m_SOIL.out").read_text()
    assert "SOIL-EQL" in out and "SOIL-NON" not in out


def test_nlsoil_opt_0_keeps_soil_eql(tmp_path, ui):
    ui.run_text(_model_text(tmp_path, "NLSOIL,0\n"), name="m.pre")
    assert SN.read_nls(tmp_path / "m.nls").options.opt == 0
    ui.execute("RUNSOIL")
    out = (tmp_path / "m_SOIL.out").read_text()
    assert "NLSOIL <Opt> = 0" in out and "ITERATION NUMBER" in out


# ---------------------------------------------------------------------------------------------
# SOIL module, SOIL-NON branch (deck + side file written directly)
# ---------------------------------------------------------------------------------------------
def _deck(tmp_path: Path, n=300, **kw) -> decks.Deck:
    rng = np.random.default_rng(4)
    a = 0.15 * rng.standard_normal(n) * np.hanning(n)
    (tmp_path / "in.acc").write_text("\n".join(f"{v:.8e}" for v in a) + "\n")
    d = decks.new("SOIL")
    d.params.update(dict(title="unit", model="m", nrval=n, grav=9.81, header=0, outcrop=1, save=1, iter=8,
                         ratio=0.65, gravmult=1.0, cof=0.0, soilcutoff=0.0, delt=0.01, nft=1024, cl=3,
                         thfile="in.acc", thtit="t", mult=0.0, max=0.2, indir=0, cmodform=0))
    d.params.update(kw)
    for i, (h, vs, lab) in enumerate([(4.0, 150.0, "Sand"), (6.0, 250.0, "Clay"), (0.0, 900.0, "")], start=1):
        d.table("profile").append([i, h, 19.0, 2.0 * vs, vs, 0.03 if h else 0.0, 0.03 if h else 0.0, lab])
    for lab, gr in (("Sand", 0.04), ("Clay", 0.1)):
        for x in (1e-4, 1e-3, 1e-2, 1e-1, 1.0):
            d.table("dynp").append([lab, "G", x, 1.0 / (1.0 + (x / gr) ** 0.9)])
        for x, dd in ((1e-4, 0.6), (1e-3, 1.0), (1e-2, 3.0), (1e-1, 9.0), (1.0, 20.0)):
            d.table("dynp").append([lab, "D", x, dd])
    d.table("sacc").append([1, 2, 0])
    d.table("sacc").append([3, 2, 1])
    d.table("srs").append([1, 1, 0])
    d.table("sstr").append([1, 1, 1, 1, 1])
    d.table("sstr").append([2, 1, 1, 1, 1])
    d.table("ssaf").append([1, 1, 0, 1, 3, 0.25, "surface / outcrop"])
    d.table("damp").append([0.05])
    return d


def _nls(tmp_path, opts=None, layers=None, name="m"):
    o = opts or SN.NLSoilOptions(opt=1, bedint=1, damptype=1)
    lays = layers if layers is not None else [SN.NLLayerData(1, 1), SN.NLLayerData(2, 1)]
    (tmp_path / f"{name}.nls").write_text(SN.nls_text(o, lays, name))


def _run(tmp_path, d, **nls):
    soil.write_deck(tmp_path / "m.soi", d)
    _nls(tmp_path, **nls)
    rc = run_module("SOIL", "m", tmp_path)
    return rc, (tmp_path / "m_SOIL.out").read_text()


def test_soil_non_outputs_and_file88(tmp_path):
    rc, out = _run(tmp_path, _deck(tmp_path))
    assert rc == 0, out[-3000:]
    for fn in ("ACC001.TH", "ACC003.TH", "SN001.TH", "SS001.TH", "SN002.TH", "SS002.TH", "RS001_01.RS",
               "SAF001W_003O.TFU", "SAF001W_003O_01.RS", "FILE73", "FILE88"):
        assert (tmp_path / fn).exists(), fn
    for text in ("Nonlinear soil layers (NLSLAYER)", "Small-strain viscous damping", "Time integration",
                 "Maximum strains and equivalent-linear properties", "Peak acceleration profile",
                 "Small-strain amplification"):
        assert text in out, text
    a3, dt = textfiles.read_history(tmp_path / "ACC003.TH")
    assert dt == pytest.approx(0.01) and len(a3) == 1024
    assert np.max(np.abs(a3)) == pytest.approx(0.2, rel=1e-9)          # outcrop at the elastic base = input
    f88 = SH.read_file88(tmp_path / "FILE88")
    sn1, _ = textfiles.read_history(tmp_path / "SN001.TH")
    geff = f88["gamma_eff_pct"]
    assert geff[0] >= 0.65 * np.max(np.abs(sn1)) * (1 - 1e-6)            # sub-step maxima >= sampled maxima
    assert geff[0] <= 0.65 * np.max(np.abs(sn1)) * 1.05
    curves = SH.curves_from_dynp_rows(_deck(tmp_path).rows("dynp"))
    for i, lab in enumerate(("Sand", "Clay")):
        fit = SN.fit_mkz(curves[lab].g_strain, curves[lab].g_ratio)
        gmax = 19.0 / 9.81 * (150.0, 250.0)[i] ** 2
        assert f88["G"][i] / gmax == pytest.approx(SN.mkz_modulus_ratio(geff[i], fit.gr_pct, 1.0, fit.s), rel=1e-6)
        xi_min = curves[lab].damping(1e-4)
        assert f88["beta_s"][i] == pytest.approx(xi_min + SN.masing_damping(geff[i], fit.gr_pct, 1.0, fit.s),
                                                 rel=1e-6)
    np.testing.assert_allclose(f88["Vp"] / f88["Vs"], 2.0, rtol=1e-7)                 # D-SOL-06
    rs = textfiles.read_xy(tmp_path / "RS001_01.RS")
    assert rs.shape == (301, 2) and np.all(rs[:, 1] > 0)


def test_soil_non_rigid_base_applies_the_input_at_the_base(tmp_path):
    d = _deck(tmp_path)
    rc, out = _run(tmp_path, d, opts=SN.NLSoilOptions(opt=1, bedint=0, damptype=3, nsub=2))
    assert rc == 0, out[-3000:]
    assert "<outcrop> = 1 is ignored" in out
    a3, _ = textfiles.read_history(tmp_path / "ACC003.TH")
    assert np.max(np.abs(a3)) == pytest.approx(0.2, rel=1e-9)
    assert "Rayleigh multipliers 0" in out


@pytest.mark.parametrize("change,message", [
    (dict(cl=1), "at bedrock"),
    (dict(indir=1), "Input Direction"),
    (dict(outcrop=0), "needs an outcrop"),
])
def test_soil_non_input_errors(tmp_path, change, message):
    rc, out = _run(tmp_path, _deck(tmp_path, **change))
    assert rc != 0 and message in out


def test_soil_non_missing_layer_and_bad_options(tmp_path):
    rc, out = _run(tmp_path, _deck(tmp_path), layers=[SN.NLLayerData(1, 1)])
    assert rc != 0 and "Error 125" in out
    rc, out = _run(tmp_path, _deck(tmp_path), opts=SN.NLSoilOptions(opt=1, damptype=7))
    assert rc != 0 and "NLDampType" in out
    rc, out = _run(tmp_path, _deck(tmp_path), layers=[SN.NLLayerData(1, 0, 1.0, 0.0, 0.05), SN.NLLayerData(2, 1)])
    assert rc != 0 and "must be > 0" in out
    (tmp_path / "m.nls").write_text("NLSOIL,1\nBOGUS,2\n")
    assert run_module("SOIL", "m", tmp_path) != 0
    assert "unknown record" in (tmp_path / "m_SOIL.out").read_text()


def test_soil_non_data_check_and_visco_elastic_damping(tmp_path):
    d = _deck(tmp_path, opmode=1)
    lays = [SN.NLLayerData(1, 0, 1.2, 0.85, 0.05, 15.0), SN.NLLayerData(2, 1)]
    rc, out = _run(tmp_path, d, opts=SN.NLSoilOptions(opt=1, bedint=1, damptype=2), layers=lays)
    assert rc == 0 and "data check only" in out
    assert not (tmp_path / "ACC001.TH").exists()
    assert "eta = 15" in out                                    # user viscosity of sublayer 1
    assert "xi(f1)" in out


def test_soil_non_side_file_hash_mismatch_warns(tmp_path):
    soil.write_deck(tmp_path / "m.soi", _deck(tmp_path, opmode=1))
    head, rest = (tmp_path / "m.soi").read_text().split("\n", 1)
    (tmp_path / "m.soi").write_text(head + "\n* model_hash = aaa\n" + rest)
    o = SN.NLSoilOptions(opt=1, bedint=1, damptype=1)
    SN.write_nls(tmp_path / "m.nls", SN.NlsData(o, {1: SN.NLLayerData(1, 1), 2: SN.NLLayerData(2, 1)}, "m", "bbb"))
    assert run_module("SOIL", "m", tmp_path) == 0
    assert "different model_hash" in (tmp_path / "m_SOIL.out").read_text()
