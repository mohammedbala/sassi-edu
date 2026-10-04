"""Unit tests of the SOIL module (deck I/O, outputs, vertical input, errors, curve library)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.core import shake as SH
from sassi.io import decks, textfiles
from sassi.modules import soil
from sassi.modules.base import run_module

DATA = Path(soil.__file__).resolve().parents[1] / "data"


def _motion(path: Path, n=600, dt=0.01, seed=0, header=0):
    rng = np.random.default_rng(seed)
    a = 0.1 * rng.standard_normal(n) * np.hanning(n)
    lines = ["header line"] * header + [f"{v:.8e}" for v in a]
    path.write_text("\n".join(lines) + "\n")
    return a


def _deck(tmp_path: Path, **kw) -> decks.Deck:
    d = decks.new("SOIL")
    d.params.update(dict(title="unit", model="m", nrval=600, grav=9.81, header=0, outcrop=1, save=1, iter=3,
                         ratio=0.65, gravmult=1.0, cof=0.0, soilcutoff=0.0, delt=0.01, nft=1024, cl=3,
                         thfile="in.acc", thtit="t", mult=0.0, max=0.1, indir=0, cmodform=0))
    d.params.update(kw)
    for i, (h, vs, lab) in enumerate([(8.0, 180.0, "Sand"), (12.0, 300.0, "Clay"), (0.0, 900.0, "")], start=1):
        d.table("profile").append([i, h, 19.0, 2.0 * vs, vs, 0.05 if h else 0.01, 0.05 if h else 0.01, lab])
    for lab in ("Sand", "Clay"):
        for x, gg, dd in ((1e-4, 1.0, 0.5), (1e-3, 0.95, 1.0), (1e-2, 0.75, 3.0), (1e-1, 0.4, 9.0), (1.0, 0.1, 20.0)):
            d.table("dynp").append([lab, "G", x, gg])
            d.table("dynp").append([lab, "D", x, dd])
    d.table("sacc").append([1, 2, 1])
    d.table("sacc").append([3, 2, 0])
    d.table("srs").append([1, 1, 1])
    d.table("sstr").append([1, 1, 1, 1, 1])
    d.table("sstr").append([2, 1, 1, 1, 1])
    d.table("ssaf").append([1, 1, 1, 1, 3, 0.25, "surface / base outcrop"])
    d.table("damp").append([0.05])
    d.table("damp").append([0.02])
    return d


def test_deck_round_trip_with_real_freqstep(tmp_path):
    d = _deck(tmp_path)
    soil.write_deck(tmp_path / "m.soi", d)
    back = soil.read_deck(tmp_path / "m.soi")
    assert back.rows("ssaf")[0]["freqstep"] == 0.25
    assert back["ratio"] == 0.65 and back["thfile"] == "in.acc"
    assert back.rows("profile")[2]["dynprop"] == ""
    assert len(back.rows("dynp")) == 20


def test_soil_run_outputs(tmp_path):
    _motion(tmp_path / "in.acc")
    soil.write_deck(tmp_path / "m.soi", _deck(tmp_path))
    assert run_module("SOIL", "m", tmp_path) == 0
    for fn in ("ACC001.TH", "ACC003.TH", "SN001.TH", "SS002.TH", "RS001_01.RS", "RS001_02.RS",
               "SAF001O_003O.TFU", "SAF001O_003O_01.RS", "FILE73", "FILE88"):
        assert (tmp_path / fn).exists(), fn
    out = (tmp_path / "m_SOIL.out").read_text()
    assert out.count("ITERATION NUMBER") == 3
    f88 = SH.read_file88(tmp_path / "FILE88")
    assert len(f88["layer"]) == 2                                   # half-space excluded (D-SOL-12)
    # constant Poisson's ratio policy (D-SOL-06): Vp/Vs unchanged, beta_p = beta_s
    np.testing.assert_allclose(f88["Vp"] / f88["Vs"], 2.0, rtol=1e-7)     # 9-digit file values
    np.testing.assert_allclose(f88["beta_p"], f88["beta_s"])
    np.testing.assert_allclose(f88["G"], 19.0 / 9.81 * f88["Vs"] ** 2, rtol=1e-8)
    a, dt = textfiles.read_history(tmp_path / "ACC003.TH")
    assert dt == pytest.approx(0.01) and len(a) == 1024
    # input is an outcrop motion at the top of the half-space: 0.1 g after MAX scaling
    rs = textfiles.read_xy(tmp_path / "RS001_01.RS")
    assert rs.shape == (301, 2) and rs[-1, 1] > 0.1                 # amplified surface ZPA
    c = SH.read_file73(tmp_path / "FILE73")
    assert [v.label for v in c.values()] == ["Sand", "Clay"]


def test_vertical_input_closed_form(tmp_path):
    """indir 1: P waves with Vp and beta_p, no iteration; surface/base-within = 1/cos(w H/Vp*)."""
    _motion(tmp_path / "in.acc")
    d = _deck(tmp_path, indir=1, iter=4)
    d.tables["profile"].rows = [[1, 25.0, 19.0, 600.0, 300.0, 0.04, 0.05, "Sand"],
                                [2, 0.0, 21.0, 4000.0, 2000.0, 0.01, 0.01, ""]]
    d.tables["sacc"].rows = []
    d.tables["sstr"].rows = []
    d.tables["srs"].rows = []
    d.tables["ssaf"].rows = [[1, 1, 1, 0, 2, 0.5, "surface / base within"]]
    d["cl"] = 2
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 0
    f, H, _ = textfiles.read_tf(tmp_path / "SAF001O_002W.TFU")
    k = 2 * np.pi * f / (600.0 * np.sqrt(cfactor(0.04)))
    np.testing.assert_allclose(H, 1.0 / np.cos(k * 25.0), rtol=1e-7)
    out = (tmp_path / "m_SOIL.out").read_text()
    assert "ITERATION NUMBER" not in out and "none performed" in out


def test_header_lines_and_mult(tmp_path):
    a = _motion(tmp_path / "in.acc", header=2)
    d = _deck(tmp_path, header=2, mult=2.0, max=0.0, iter=0)
    d.tables["sacc"].rows = [[3, 2, 1]]
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 0
    out, _ = textfiles.read_history(tmp_path / "ACC003.TH")
    np.testing.assert_allclose(out[:600], 2.0 * a, atol=1e-7)       # outcrop control motion reproduced


@pytest.mark.parametrize("kw,msg", [(dict(mult=1.0, max=0.1), "Error 78"), (dict(mult=0.0, max=0.0), "Error 77"),
                                    (dict(ratio=1.2), "Error 106"), (dict(cl=9), "Error 104"),
                                    (dict(thfile="missing.acc"), "Error 73"), (dict(nft=256), "EDU-03")])
def test_input_errors(tmp_path, kw, msg):
    _motion(tmp_path / "in.acc")
    soil.write_deck(tmp_path / "m.soi", _deck(tmp_path, **kw))
    assert run_module("SOIL", "m", tmp_path) == 1
    assert msg in (tmp_path / "m_SOIL.out").read_text()


def test_undefined_dynamic_property(tmp_path):
    _motion(tmp_path / "in.acc")
    d = _deck(tmp_path)
    d.tables["profile"].rows[0][-1] = "Gravel"
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 1
    assert "Gravel" in (tmp_path / "m_SOIL.out").read_text()


def _parse_dynp_pre(path: Path):
    curves = {}
    for ln in path.read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith("*"):
            continue
        tok = s.split(",")
        assert tok[0] == "DYNP" and len(tok) == 7
        no, sg, g, sd, dd, lab = int(tok[1]), *map(float, tok[2:6]), tok[6]
        c = curves.setdefault(lab, {"G": [], "D": []})
        c["G"].append((sg, g))
        c["D"].append((sd, dd))
    return curves


def test_dynp_library_matches_shake91_curves():
    """sassi/data/dynp_library.pre reproduces the INP.DAT curves under log-strain interpolation."""
    from sassi.verify.problems.vp_soil import parse_shake91_input
    lib = _parse_dynp_pre(DATA / "dynp_library.pre")
    assert sorted(lib) == ["Clay", "Rock", "Sand"]
    inp = parse_shake91_input(DATA / "shake91_example" / "INP.DAT")
    mats = dict(zip(["Clay", "Sand", "Rock"], inp["materials"]))
    probe = np.logspace(-5, 1.5, 400)
    for lab, m in mats.items():
        g = np.array(lib[lab]["G"])
        dd = np.array(lib[lab]["D"])
        np.testing.assert_allclose(SH.interp_log_strain(probe, g[:, 0], g[:, 1]),
                                   SH.interp_log_strain(probe, m["g_strain"], m["g_ratio"]), atol=1e-12)
        np.testing.assert_allclose(SH.interp_log_strain(probe, dd[:, 0], dd[:, 1]),
                                   SH.interp_log_strain(probe, m["d_strain"], m["d_pct"]), atol=1e-6)


def test_shake91_parser():
    from sassi.verify.problems.vp_soil import parse_shake91_input
    inp = parse_shake91_input(DATA / "shake91_example" / "INP.DAT")
    assert len(inp["materials"]) == 3 and len(inp["layers"]) == 17
    assert inp["layers"][0] == dict(no=1, type=2, thick=5.0, gmax=None, damp=0.05, weight=0.125, vs=1000.0)
    assert inp["layers"][-1]["thick"] == 0.0 and inp["layers"][-1]["vs"] == 4000.0
    m = inp["motion"]
    assert (m["nv"], m["nfft"], m["dt"], m["xmax"], m["fmax"], m["nhead"]) == (1900, 4096, 0.02, 0.10, 25.0, 3)
    assert inp["input_layer"] == 17 and inp["outcrop"] is True
    assert inp["iterations"] == 8 and inp["ratio"] == 0.5
    assert inp["amplification"]["df"] == 0.125


def test_data_check_mode(tmp_path):
    _motion(tmp_path / "in.acc")
    soil.write_deck(tmp_path / "m.soi", _deck(tmp_path, opmode=1))
    assert run_module("SOIL", "m", tmp_path) == 0
    assert "data check only" in (tmp_path / "m_SOIL.out").read_text()
    assert not (tmp_path / "FILE88").exists()


def test_file88_linear_sublayer_keeps_input_p_wave_properties(tmp_path):
    """D-SOL-06 (Vp at constant Poisson's ratio, beta_p = beta_s) applies to iterated sublayers
    only; an unlabelled (linear) sublayer keeps its input Vp and P-wave damping."""
    _motion(tmp_path / "in.acc")
    d = _deck(tmp_path, max=0.3)
    d.tables["profile"].rows = [[1, 8.0, 19.0, 360.0, 180.0, 0.05, 0.05, "Sand"],
                                [2, 12.0, 19.0, 900.0, 300.0, 0.02, 0.06, ""],
                                [3, 0.0, 19.0, 1800.0, 900.0, 0.015, 0.01, ""]]
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 0
    f88 = SH.read_file88(tmp_path / "FILE88")
    assert f88["Vs"][1] == pytest.approx(300.0) and f88["beta_s"][1] == pytest.approx(0.06)
    assert f88["Vp"][1] == pytest.approx(900.0) and f88["beta_p"][1] == pytest.approx(0.02)
    # iterated sublayer: softened, Vp/Vs kept, beta_p = beta_s
    assert f88["Vs"][0] < 180.0
    assert f88["Vp"][0] / f88["Vs"][0] == pytest.approx(2.0, rel=1e-7)
    assert f88["beta_p"][0] == pytest.approx(f88["beta_s"][0])
    out = (tmp_path / "m_SOIL.out").read_text()
    assert "linear sublayers keep their input Vp and beta_p" in out


def test_final_table_and_file88_report_the_same_compatible_strain(tmp_path):
    _motion(tmp_path / "in.acc")
    d = _deck(tmp_path, iter=1, max=0.3)
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 0
    f88 = SH.read_file88(tmp_path / "FILE88")
    out = (tmp_path / "m_SOIL.out").read_text()
    sec = out.split("Strain-compatible soil properties (final)")[1].splitlines()[3:5]
    eff = [float(ln.split()[3]) for ln in sec]
    np.testing.assert_allclose(eff, f88["gamma_eff_pct"], atol=5.1e-6)          # listing: 5 decimals
    it1 = out.split("ITERATION NUMBER   1")[1].splitlines()[4:6]
    np.testing.assert_allclose([float(ln.split()[3]) for ln in it1], f88["gamma_eff_pct"], atol=5.1e-6)


def test_ssaf_rs_ratio_values(tmp_path):
    """SAF RS file = RS(layer)/RS(layer2) of the layer motions at f = k * freqstep (k >= 1)."""
    from sassi.core import spectra as SP
    _motion(tmp_path / "in.acc")
    d = _deck(tmp_path, iter=0)
    d.tables["ssaf"].rows = [[1, 1, 1, 0, 3, 2.5, "surface / base within"]]
    d.tables["sacc"].rows = [[1, 2, 1], [3, 2, 0]]
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 0
    rat = textfiles.read_xy(tmp_path / soil.saf_rs_name(1, True, 3, False, 1))
    np.testing.assert_allclose(rat[:, 0], 2.5 * np.arange(1, 21), atol=1e-6)       # up to Nyquist 50 Hz
    a1, dt = textfiles.read_history(tmp_path / "ACC001.TH")
    a3, _ = textfiles.read_history(tmp_path / "ACC003.TH")
    r1 = SP.response_spectrum(a1, dt, rat[:, 0], [0.05])["SA"][0]
    r3 = SP.response_spectrum(a3, dt, rat[:, 0], [0.05])["SA"][0]
    np.testing.assert_allclose(rat[:, 1], r1 / r3, rtol=1e-5)
