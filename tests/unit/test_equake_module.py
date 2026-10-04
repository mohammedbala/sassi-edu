"""Unit tests of the EQUAKE module (external history, generation, units, errors, correlation)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sassi.core import equake_lib as EL
from sassi.core import spectra as SP
from sassi.io import decks, textfiles
from sassi.modules.base import run_module

FREQS = [0.2, 0.25, 0.5, 1.0, 2.5, 5.0, 9.0, 20.0, 33.0, 50.0]


def _target(tmp_path: Path, name="t.rsi", pga=0.3, comp="H", freqs=FREQS):
    f = np.asarray(freqs)
    EL.write_spectrum_file(tmp_path / name, f, EL.rg160_spectrum(f, comp, 0.05, pga))
    return name


def _deck(tmp_path: Path, rows, **kw):
    d = decks.new("EQUAKE")
    d.params.update(dict(title="unit", model="e", accopt=0, nrfreq=len(FREQS), rand=11975, damp=0.05, dur=20.0,
                         corr=0, seeds=1, tpsd=0, delt=0.01, gravity=32.2, eqtit="unit test"))
    d.params.update(kw)
    for r in rows:
        d.table("spectra").append(r)
    return d


def test_external_history_rs_psd_fft(tmp_path):
    rng = np.random.default_rng(0)
    dt = 0.01
    t = np.arange(2001) * dt
    a = 0.2 * rng.standard_normal(2001) * EL.saragoni_hart_envelope(t, 20.0)
    textfiles.write_history(tmp_path / "rec.acc", a, dt)
    _target(tmp_path)
    decks.write(tmp_path / "e.equ", _deck(tmp_path, [[1, "t.rsi", "rec.rso", "rec.acc", "", ""]], accopt=2))
    assert run_module("EQUAKE", "e", tmp_path) == 0
    rso = textfiles.read_xy(tmp_path / "rec.rso")
    a_back, _ = textfiles.read_history(tmp_path / "rec.acc")
    grid = EL.check_grid(0.2, 50.0)                          # manual band, 100 points per decade
    np.testing.assert_allclose(rso[:, 0], grid, atol=5e-7)   # frequencies printed with 6 decimals
    sa = SP.response_spectrum(a_back, dt, grid, [0.05])["SA"][0]
    np.testing.assert_allclose(rso[:, 1], sa, rtol=1e-7)
    assert rso[0, 0] == pytest.approx(0.2) and rso[-1, 0] == pytest.approx(50.0)
    assert (tmp_path / "rec.psd").exists() and (tmp_path / "rec.fft").exists()
    fft = textfiles.read_xy(tmp_path / "rec.fft")
    assert fft.shape[1] == 3 and fft[0, 0] == 0.0
    out = (tmp_path / "e_EQUAKE.out").read_text()
    assert "Acceptance criteria summary" in out and "no simulation" in out


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    wd = tmp_path_factory.mktemp("eqk")
    _target(wd)
    decks.write(wd / "e.equ", _deck(wd, [[1, "t.rsi", "H1.rso", "", "H1", ""]]))
    rc = run_module("EQUAKE", "e", wd)
    return wd, rc


def test_generation_outputs_and_units(generated):
    wd, rc = generated
    assert rc == 0
    a, dt = textfiles.read_history(wd / "H1.acc")            # ACCOUT without extension -> .acc
    assert dt == pytest.approx(0.01) and len(a) == 2001
    vel, _ = textfiles.read_history(wd / "H1.vel")
    dis, _ = textfiles.read_history(wd / "H1.dis")
    v, d = SP.integrate(a * 32.2 * 12.0, dt)                 # British: in/s and in (D-EQK-07)
    np.testing.assert_allclose(vel, v, atol=1e-6 * np.max(np.abs(v)))
    np.testing.assert_allclose(dis, d, atol=1e-6 * np.max(np.abs(d)))
    assert abs(v[-1]) < 1e-6 * np.max(np.abs(v)) and abs(d[-1]) < 1e-6 * np.max(np.abs(d))
    tg = EL.TargetSpectrum.from_file(wd / "t.rsi")
    chk = EL.evaluate_rs(a, dt, tg, 0.05)
    assert 0.85 < chk["manual"].min_ratio and chk["manual"].max_ratio < 1.35
    out = (wd / "e_EQUAKE.out").read_text()
    for key in ("PGA =", "V/A =", "AD/V^2", "strong-motion duration", "SRP 3.7.1 Rev. 4", "ACS SASSI manual",
                "in/s", "Generation log"):
        assert key in out, key
    assert "time step 0.01 s > 0.005 s" in out               # manual warning


def test_generation_is_deterministic(generated, tmp_path):
    wd, _ = generated
    _target(tmp_path)
    decks.write(tmp_path / "e.equ", _deck(tmp_path, [[1, "t.rsi", "H1.rso", "", "H1", ""]]))
    assert run_module("EQUAKE", "e", tmp_path) == 0
    assert (tmp_path / "H1.acc").read_text() == (wd / "H1.acc").read_text()


def test_seed_record_and_errors(tmp_path):
    rng = np.random.default_rng(2)
    dt = 0.01
    t = np.arange(2001) * dt
    seed = 0.1 * rng.standard_normal(2001) * np.exp(-((t - 8) / 4) ** 2)
    textfiles.write_history(tmp_path / "seed.acc", seed, dt)
    _target(tmp_path)
    decks.write(tmp_path / "e.equ", _deck(tmp_path, [[1, "t.rsi", "S.rso", "seed.acc", "S.acc", ""]], accopt=1))
    assert run_module("EQUAKE", "e", tmp_path) == 0
    assert "seed record" in (tmp_path / "e_EQUAKE.out").read_text().lower()
    # Error 89: record count differs from <nrfreq>
    decks.write(tmp_path / "e.equ", _deck(tmp_path, [[1, "t.rsi", "S.rso", "", "S.acc", ""]], nrfreq=5))
    assert run_module("EQUAKE", "e", tmp_path) == 1
    assert "Error 89" in (tmp_path / "e_EQUAKE.out").read_text()
    # Error 84: no spectrum input file
    decks.write(tmp_path / "e.equ", _deck(tmp_path, [[1, "", "S.rso", "", "S.acc", ""]]))
    assert run_module("EQUAKE", "e", tmp_path) == 1
    assert "Error 84" in (tmp_path / "e_EQUAKE.out").read_text()
    # Error 93: correlated without pairs
    decks.write(tmp_path / "e.equ", _deck(tmp_path, [[1, "t.rsi", "S.rso", "", "S.acc", ""]], corr=1))
    assert run_module("EQUAKE", "e", tmp_path) == 1
    assert "Error 93" in (tmp_path / "e_EQUAKE.out").read_text()


def test_correlated_components(tmp_path):
    _target(tmp_path)
    d = _deck(tmp_path, [[1, "t.rsi", "X.rso", "", "X.acc", ""], [2, "t.rsi", "Y.rso", "", "Y.acc", ""]], corr=1)
    d.table("corr").append([1, 0.0, 0.4])
    d.table("corr").append([2, 20.0, 0.4])
    decks.write(tmp_path / "e.equ", d)
    assert run_module("EQUAKE", "e", tmp_path) == 0
    x, _ = textfiles.read_history(tmp_path / "X.acc")
    y, _ = textfiles.read_history(tmp_path / "Y.acc")
    out = (tmp_path / "e_EQUAKE.out").read_text()
    assert "Correlated components" in out and "2-s moving window" in out
    assert EL.correlation(x, y) > 0.15                       # correlation imposed (re-matching changes it)


def test_data_check_mode(tmp_path):
    _target(tmp_path)
    decks.write(tmp_path / "e.equ", _deck(tmp_path, [[1, "t.rsi", "S.rso", "", "S.acc", ""]], opmode=1))
    assert run_module("EQUAKE", "e", tmp_path) == 0
    assert "data check only" in (tmp_path / "e_EQUAKE.out").read_text()
    assert not (tmp_path / "S.acc").exists()


def test_failed_criteria_are_warned_in_listing(generated):
    """Req. 4.12 item 6: a failed criterion is a listing (and screen) warning, not only a FAIL row.
    The 0.2-50 Hz target cannot demonstrate SRP 3.7.1 (b) (0.1-50 Hz): incomplete band."""
    wd, rc = generated
    out = (wd / "e_EQUAKE.out").read_text()
    assert "(incomplete: 0.1 - 50 Hz required)" in out
    assert "*** WARNING: spectrum 1: SRP 3.7.1 (b) not met" in out
    assert "criterion check(s) failed" in out
    n_warn = out.count("*** WARNING")
    assert f"{n_warn} warning(s)" in out


def test_external_history_failures_warned_with_echo(tmp_path):
    """Every failed criterion of an external history is warned; warnings reach the screen echo."""
    rng = np.random.default_rng(9)
    dt = 0.005
    t = np.arange(4001) * dt
    a = 0.05 * rng.standard_normal(4001) * EL.saragoni_hart_envelope(t, 20.0)
    textfiles.write_history(tmp_path / "rec.acc", a, dt)
    f = np.asarray([0.1] + FREQS)
    EL.write_spectrum_file(tmp_path / "t.rsi", f, EL.rg160_spectrum(f, "H", 0.05, 0.5))
    decks.write(tmp_path / "e.equ", _deck(tmp_path, [[1, "t.rsi", "rec.rso", "rec.acc", "", ""]], accopt=2,
                                          delt=dt, gravity=9.81, nrfreq=len(f)))
    screen = []
    assert run_module("EQUAKE", "e", tmp_path, echo=screen.append) == 0
    warned = [ln for ln in screen if "*** WARNING" in ln]
    assert any("more than 10 % below the target" in ln and "manual" in ln for ln in warned)
    assert any("SRP 3.7.1 (c) not met" in ln for ln in warned)
    assert not any("(b) not met" in ln for ln in warned)          # 0.1-50 Hz covered


def test_correlation_requires_spectra_1_and_2(tmp_path):
    """D-EQK-06: the CORR pairs relate X (spectrum 1) and Y (spectrum 2), chosen by number."""
    _target(tmp_path)
    d = _deck(tmp_path, [[1, "t.rsi", "X.rso", "", "X.acc", ""], [3, "t.rsi", "Z.rso", "", "Z.acc", ""]], corr=1)
    d.table("corr").append([1, 0.0, 0.5])
    d.table("corr").append([2, 20.0, 0.5])
    decks.write(tmp_path / "e.equ", d)
    assert run_module("EQUAKE", "e", tmp_path) == 0
    out = (tmp_path / "e_EQUAKE.out").read_text()
    assert "spectra 1 (X) and 2 (Y) are both required" in out and "re-generated" not in out


def test_correlated_rematch_uses_component_target_psd(tmp_path, monkeypatch):
    """The re-matched Y component is generated with its own target PSD (tpsd = 1)."""
    calls = []
    orig = EL.match_spectrum

    def spy(*args, **kw):
        calls.append(kw)
        kw["options"] = EL.MatchOptions(lw_iterations=3, wavelet_iterations=2, outer_passes=1)
        return orig(*args, **kw)
    monkeypatch.setattr(EL, "match_spectrum", spy)
    _target(tmp_path)
    fp = np.logspace(-1, 2, 61)
    textfiles.write_xy(tmp_path / "t.tpsd", fp, EL.rg160_target_psd(fp, 0.3, "BS"))
    d = _deck(tmp_path, [[2, "t.rsi", "Y.rso", "", "Y.acc", "t.tpsd"], [1, "t.rsi", "X.rso", "", "X.acc", "t.tpsd"]],
              corr=1, tpsd=1)
    d.table("corr").append([1, 0.0, 0.4])
    d.table("corr").append([2, 20.0, 0.4])
    decks.write(tmp_path / "e.equ", d)
    assert run_module("EQUAKE", "e", tmp_path) == 0
    assert len(calls) == 3
    assert calls[-1]["seed_acc"] is not None and calls[-1]["target_psd"] is not None
    out = (tmp_path / "e_EQUAKE.out").read_text()
    assert "spectrum 2 (Y) re-generated" in out and "target PSD floor applied" in out
