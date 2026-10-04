"""Unit tests of the lead-owned shared code: conventions, decks, containers, text files, spectra."""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from sassi import conventions as C
from sassi.core import spectra as S
from sassi.io import decks, textfiles, thfile
from sassi.io.container import read_container, write_container
from sassi.io.files import validate


# ---------------------------------------------------------------- conventions
def test_cfactor_forms():
    b = 0.05
    c = C.cfactor(b)
    assert abs(abs(c) - 1.0) < 1e-14                      # |c| = 1 for the SASSI form
    assert abs(c - (1 - 2 * b * b + 2j * b * np.sqrt(1 - b * b))) < 1e-15
    assert C.cfactor(b, C.CMOD_SIMPLE) == 1 + 0.1j
    with pytest.raises(ValueError):
        C.cfactor(0.5)


def test_sdof_peak_sassi_form():
    # hysteretic SDOF |H|max = 1/(2b sqrt(1-b^2)) for the SASSI form (VP-01 reference, C-28)
    b = 0.05
    k = C.cfactor(b)
    r = np.linspace(0.9, 1.1, 200001)
    H = k / (k - r ** 2)                                  # total/ground displacement ratio
    assert abs(np.abs(H).max() - 1 / (2 * b * np.sqrt(1 - b * b))) < 1e-4


def test_frequency_step_and_power_of_two():
    assert C.frequency_step(0.005, 4096) == pytest.approx(0.048828125)
    assert C.frequency_step(0.005, 4096, fstep=0.5) == 0.5
    assert C.nearest_power_of_two(3000) == 2048
    assert C.nearest_power_of_two(3072) == 4096           # tie goes up
    assert C.is_power_of_two(4096) and not C.is_power_of_two(3000)


def test_naming_formulas():
    assert C.nodal_result_name(12, 1, "TFU") == "00012TR_X.TFU"
    assert C.nodal_rs_name(12, 1, 1) == "00012TR_X01.RS"
    assert C.nodal_result_name(123456, 6, "ACC") == "123456R_ZZ.ACC"
    assert C.element_result_name("BEAMS", 3, 45, "MXJ", "THS") == "BEAMS_003_00045_MXJ.THS"
    assert C.incoherent_file8_name(17, 2) == "FILE8050"
    assert C.incoherent_file8_name(50, 3) == "FILE8150"
    assert C.restart_names(1) == ("COOX001", "COOTK001")
    assert C.layer_th_name("ACC", 7) == "ACC007.TH"


def test_moduli_conversions():
    G, M, lam = C.moduli_from_E_nu(30000.0, 0.25)
    rho = 0.15 / 32.2
    vs, vp = np.sqrt(G / rho), np.sqrt(M / rho)
    G2, M2, lam2, nu2, E2 = C.moduli_from_velocities(vp, vs, rho)
    assert G2 == pytest.approx(G, rel=1e-12) and nu2 == pytest.approx(0.25, rel=1e-12)
    assert E2 == pytest.approx(30000.0, rel=1e-12)


# ---------------------------------------------------------------- decks
def test_deck_roundtrip_all_schemas(tmp_path):
    for name, schema in decks.SCHEMAS.items():
        d = decks.new(name)
        d["title"] = 'Title, with "quotes" = and spaces'
        for t in schema.tables:
            row = []
            for col, typ in t.columns:
                row.append({int: 3, float: 1.25, str: "a b"}[typ])
            d.table(t.name).append(row)
        p = tmp_path / f"m{schema.ext or '.x'}"
        decks.write(p, d)
        d2 = decks.read(p, name)
        assert d2.params == {k: d.params[k] for k in d2.params}
        for t in schema.tables:
            assert d2.table(t.name).rows == d.table(t.name).rows


def test_deck_rejects_unknown_param(tmp_path):
    d = decks.new("SITE")
    d["bogus"] = 1
    with pytest.raises(KeyError):
        decks.write(tmp_path / "m.sit", d)


def test_deck_integer_coercion(tmp_path):
    d = decks.new("SITE")
    d["nl"] = 20.0
    decks.write(tmp_path / "m.sit", d)
    assert decks.read(tmp_path / "m.sit")["nl"] == 20
    d["nl"] = 20.5
    with pytest.raises(ValueError):
        decks.write(tmp_path / "m.sit", d)


# ---------------------------------------------------------------- containers
def test_container_roundtrip_and_validate(tmp_path):
    K = sp.random(6, 6, density=0.3, random_state=1, format="csr") + 1j * sp.eye(6)
    arrays = {"fnum": np.array([2, 4]), "freq": np.array([0.1, 0.2]), "eq_node": np.arange(3),
              "eq_dof": np.ones(3, int), "H": np.ones((2, 3), complex), "K": K}
    meta = {k: 0 for k in ["df", "type", "case", "ang", "cm", "nfft", "delt", "model_hash"]}
    write_container(tmp_path / "FILE8", "FILE8", arrays, meta, module="ANALYS")
    c = read_container(tmp_path / "FILE8", "FILE8")
    assert validate(c) == []
    assert abs(c.sparse("K") - K).max() == 0
    with pytest.raises(ValueError):
        read_container(tmp_path / "FILE8", "FILE1")
    del c.arrays["H"]
    assert any("missing array 'H'" in p for p in validate(c))


# ---------------------------------------------------------------- text files
def test_textfiles_roundtrip(tmp_path):
    f = np.linspace(0, 10, 11)
    H = np.exp(1j * f) * (1 + f)
    textfiles.write_tf(tmp_path / "a.TFI", f, H, complex_=True, header="node 1")
    f2, H2, cplx = textfiles.read_tf(tmp_path / "a.TFI")
    assert cplx and np.allclose(f2, f) and np.allclose(H2, H, rtol=1e-7)
    textfiles.write_history(tmp_path / "a.ACC", np.sin(f), 0.01)
    v, dt = textfiles.read_history(tmp_path / "a.ACC")
    assert dt == 0.01 and np.allclose(v, np.sin(f))


def test_thfile_formats_and_scaling(tmp_path):
    p0 = tmp_path / "h0.acc"
    p0.write_text("0.01\n1.0 2.0\n-3.0\n4.0D0\n")
    a, dt = thfile.read_history(p0, fopt=0)
    assert dt == 0.01 and list(a) == [1.0, 2.0, -3.0, 4.0]
    a2, _ = thfile.read_history(p0, fopt=0, rec1=2, rec2=3)
    assert list(a2) == [2.0, -3.0]
    p1 = tmp_path / "h1.acc"
    p1.write_text("0.0 1\n0.02 2\n0.04 3\n")
    a, dt = thfile.read_history(p1, fopt=1)
    assert dt == pytest.approx(0.02) and list(a) == [1, 2, 3]
    s = thfile.scale_history(np.array([1.0, -2.0]), mult=0.0, maxval=0.5)
    assert np.allclose(s, [0.25, -0.5])                   # UT-07 rule a*max/max|a|
    assert np.allclose(thfile.scale_history(np.array([1.0, -2.0]), mult=2.0, maxval=0.5), [2, -4])


# ---------------------------------------------------------------- spectra (VP-30 values, R2 G.1)
def test_spectra_step_and_harmonic():
    z = 0.05
    dt = 0.0005
    r = S.response_spectrum(np.ones(int(4 / dt)), dt, [2.0], [z], upsample=False)
    assert r["PSA"][0, 0] == pytest.approx(1 + np.exp(-z * np.pi / np.sqrt(1 - z * z)), rel=2e-5)
    t = np.arange(int(300 / 0.002)) * 0.002
    r = S.response_spectrum(np.sin(2 * np.pi * t), 0.002, [1.0], [z], upsample=False)
    assert r["PSA"][0, 0] == pytest.approx(10.0, rel=1e-3)
    assert r["SA"][0, 0] == pytest.approx(np.sqrt(1 + 4 * z * z) / (2 * z), rel=1e-3)


def test_spectra_high_frequency_limit_is_pga():
    rng = np.random.default_rng(3)
    a = rng.standard_normal(2048)
    r = S.response_spectrum(a, 0.01, [300.0], [0.05])
    assert r["SA"][0, 0] == pytest.approx(np.abs(a).max(), rel=0.02)


def test_spectra_filter_matches_loop():
    rng = np.random.default_rng(5)
    a = rng.standard_normal(1500)
    fr = S.log_frequencies(0.2, 40.0, 25)
    r1 = S.response_spectrum(a, 0.01, fr, [0.02, 0.07], upsample=False)
    r2 = S.response_spectrum(a, 0.01, fr, [0.02, 0.07], upsample=False, method="loop")
    for k in ("SA", "SV", "SD"):
        assert np.allclose(r1[k], r2[k], rtol=1e-9)


def test_arias_and_psd():
    dt = 0.01
    t = np.arange(2000) * dt
    a = np.where((t > 5) & (t < 15), 1.0, 0.0)
    t0, t1 = S.strong_motion_window(a, dt)
    assert t0 == pytest.approx(5.5, abs=0.02) and t1 == pytest.approx(12.5, abs=0.02)
    f, p = S.band_averaged_psd(a, dt)
    assert np.all(p >= 0) and len(f) == len(p)
