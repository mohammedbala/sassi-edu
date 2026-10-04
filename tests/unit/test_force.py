"""Unit tests of the FORCE module (requirements 4.5, D-FRC-01..03, D-MDL-08)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.io import decks
from sassi.io.files import SCHEMAS, read_container, validate
from sassi.modules.base import run_module
from sassi.modules.force import combine_loads, load_amplitudes


def run_force(tmp_path, loads, fnum=(1, 2, 5), mforce=1, **params):
    d = decks.new("FORCE")
    d["model"] = "m"
    d["mforce"] = mforce
    d["delt"] = 0.01
    d["nft"] = 1024
    for k, v in params.items():
        d[k] = v
    for row in loads:
        d.table("loads").append(list(row))
    for n in fnum:
        d.table("freqs").append([n])
    decks.write(tmp_path / "m.frc", d)
    rc = run_module("FORCE", "m", tmp_path)
    return rc, (tmp_path / "m_FORCE.out").read_text()


def test_file9_values_and_phases(tmp_path):
    loads = [(12, 1, 0.5, 0.1), (12, 3, -2.0, 0.0), (7, 5, 1.5, 0.25)]
    rc, out = run_force(tmp_path, loads, fnum=(5, 1, 2))
    assert rc == 0, out
    f9 = read_container(tmp_path / "FILE9", "FILE9")
    assert validate(f9) == []
    arrays, meta = SCHEMAS["FILE9"]
    assert set(arrays) <= set(f9.arrays) and set(meta) <= set(f9.meta)
    df = 1.0 / (0.01 * 1024)
    assert f9.meta["df"] == pytest.approx(df, rel=1e-15)
    np.testing.assert_array_equal(f9["fnum"], [1, 2, 5])                  # sorted ascending
    np.testing.assert_allclose(f9["freq"], np.array([1, 2, 5]) * df, rtol=1e-15)
    # loaded DOFs sorted by (node, dof)
    np.testing.assert_array_equal(f9["load_node"], [7, 12, 12])
    np.testing.assert_array_equal(f9["load_dof"], [5, 1, 3])
    w = 2 * np.pi * f9["freq"]
    P = f9["P"]
    assert P.shape == (3, 3)
    np.testing.assert_allclose(P[:, 0], 1.5 * np.exp(-1j * w * 0.25), rtol=1e-14)
    np.testing.assert_allclose(P[:, 1], 0.5 * np.exp(-1j * w * 0.1), rtol=1e-14)
    np.testing.assert_allclose(P[:, 2], -2.0, rtol=0)
    # phase lag grows linearly with frequency: arg P = -w t (mod 2 pi)
    np.testing.assert_allclose(np.angle(P[:, 1] / 0.5), np.angle(np.exp(-1j * w * 0.1)), atol=1e-14)
    assert "Loads: P(f) = factor * exp(-i 2 pi f t_arrival)" in out and "MY" in out


def test_overwrite_mode_last_definition_wins(tmp_path):
    loads = [(3, 2, 1.0, 0.0), (3, 2, 4.0, 0.2)]
    rc, out = run_force(tmp_path, loads, mforce=1)
    assert rc == 0, out
    f9 = read_container(tmp_path / "FILE9")
    np.testing.assert_allclose(f9["x_factor"], [4.0])
    np.testing.assert_allclose(f9["x_arrival"], [0.2])
    np.testing.assert_allclose(f9["P"][:, 0], 4.0 * np.exp(-2j * np.pi * f9["freq"] * 0.2), rtol=1e-14)
    assert "replaces the earlier one" in out


def test_add_mode_adds_factors_and_keeps_latest_arrival(tmp_path):
    """D-FRC-02 / D-MDL-08 (spec 08 section 9.5): factors add; the arrival time of the most recent
    definition with a non-zero factor replaces the earlier one, with a warning."""
    loads = [(3, 2, 1.0, 0.3), (3, 2, 4.0, 0.2), (3, 2, -0.5, 0.1)]
    rc, out = run_force(tmp_path, loads, mforce=0)
    assert rc == 0, out
    f9 = read_container(tmp_path / "FILE9")
    np.testing.assert_allclose(f9["x_factor"], [4.5])
    np.testing.assert_allclose(f9["x_arrival"], [0.1])
    assert "arrival time 0.2 s replaced by 0.1 s" in out and "D-FRC-02" in out and "WARNING" in out


def test_add_mode_zero_factor_rows_never_change_the_arrival():
    """A zero-factor row adds no load, so it must not move the phase of an existing load; a
    zero-factor placeholder followed by a real load takes the real load's arrival time."""
    r = [dict(node=1, dof=1, factor=1.0, arrival=0.1), dict(node=1, dof=1, factor=0.0, arrival=0.5),
         dict(node=2, dof=3, factor=0.0, arrival=0.9), dict(node=2, dof=3, factor=2.0, arrival=0.2)]
    loads, warns = combine_loads(r, 0)
    assert loads == {(1, 1): (1.0, 0.1), (2, 3): (2.0, 0.2)}
    assert "does not change the arrival time 0.1 s" in warns[0] and "replaced" not in warns[1]


def test_add_mode_matches_the_interpreter_rule(tmp_path):
    """FORCE add mode on repeated deck rows gives the same (factor, arrival) as the F command in
    MOPT <force> = 0 (sassi/prep/commands/loads.py) for the same sequence."""
    from sassi.prep import Interpreter
    seq = [(1.0, 0.3), (0.0, 0.7), (4.0, 0.2), (-0.5, 0.1)]
    ui = Interpreter(cwd=tmp_path)
    ui.run_text("MDL,u,{d}\nN,3,0,0,0\nMOPT,1,0,1,0\n".format(d=tmp_path)
                + "".join(f"F,3,{a},0,0,{t},0,0\n" for a, t in seq))
    ld = ui.model.forces[3]
    loads, _ = combine_loads([dict(node=3, dof=1, factor=a, arrival=t) for a, t in seq], 0)
    assert loads[(3, 1)] == (pytest.approx(ld.factor[0]), pytest.approx(ld.arrival[0]))


def test_zero_total_factor_dropped(tmp_path):
    rc, out = run_force(tmp_path, [(1, 1, 1.0, 0.0), (1, 1, -1.0, 0.0), (2, 3, 1.0, 0.0)], mforce=0)
    assert rc == 0, out
    f9 = read_container(tmp_path / "FILE9")
    np.testing.assert_array_equal(f9["load_node"], [2])
    assert "zero total factor" in out


def test_no_loads_is_error_61(tmp_path):
    rc, out = run_force(tmp_path, [])
    assert rc == 1 and "Error 61" in out and not (tmp_path / "FILE9").exists()


def test_bad_dof_and_empty_frequency_set(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    rc, out = run_force(tmp_path / "a", [(1, 7, 1.0, 0.0)])
    assert rc == 1 and "dof 7 must be 1..6" in out
    rc, out = run_force(tmp_path / "b", [(1, 1, 1.0, 0.0)], fnum=())
    assert rc == 1 and "Error 120" in out


def test_explicit_frequency_step(tmp_path):
    """Single-harmonic runs give df directly (fstep > 0); delt and NFFT are not used."""
    rc, out = run_force(tmp_path, [(1, 1, 1.0, 0.05)], fnum=(3,), fstep=0.75, delt=0.0, nft=0)
    assert rc == 0, out
    f9 = read_container(tmp_path / "FILE9")
    assert f9.meta["df"] == 0.75 and f9["freq"][0] == pytest.approx(2.25)


def test_data_check_mode(tmp_path):
    rc, out = run_force(tmp_path, [(1, 1, 1.0, 0.0)], opmode=1)
    assert rc == 0 and "Data check only" in out and not (tmp_path / "FILE9").exists()


def test_arrival_time_is_a_circular_shift():
    """D-FRC-01 / VP-40 at unit level: irfft(P(w) F_ref(w)) = a f_ref(t - t0) (circular shift)
    when t0 is a multiple of the time step."""
    n, dt = 256, 0.01
    t = np.arange(n) * dt
    f_ref = np.exp(-((t - 0.4) / 0.05) ** 2) * np.sin(2 * np.pi * 8 * t)
    freq = np.fft.rfftfreq(n, dt)
    a, t0 = 0.5, 0.1
    P = load_amplitudes(freq, np.array([a]), np.array([t0]))[:, 0]
    y = np.fft.irfft(P * np.fft.rfft(f_ref), n)
    np.testing.assert_allclose(y, a * np.roll(f_ref, int(round(t0 / dt))), atol=1e-12)
