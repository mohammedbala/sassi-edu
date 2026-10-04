"""Line objects and line mathematics (requirements 4.13, D-LIN-01..04, D-FIL-02, D-UI-16; spec 10 4.1-4.2):
the algorithms of sassi.plotting.lines and the READSPEC ... BROADEN, LBINCORS commands."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.plotting import lines as LM
from sassi.plotting.lines import Line, LineError
from sassi.plotting.state import plot_state
from sassi.prep import Interpreter, Kind


@pytest.fixture
def ui(tmp_path):
    u = Interpreter(cwd=tmp_path)
    plot_state(u).auto_render = False
    return u


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


def write_xy(path, x, y):
    path.write_text("\n".join(f"{float(a)!r} {float(b)!r}" for a, b in zip(x, y)) + "\n")


def lines(ui):
    return plot_state(ui).lines


# ------------------------------------------------------------------ union grid and resampling
def test_union_grid_merges_near_duplicates_and_keeps_first():
    X = LM.union_grid([np.array([0.0, 1.0, 2.0]), np.array([1.0 + 1e-12, 3.0])])
    assert X.tolist() == [0.0, 1.0, 2.0, 3.0]
    assert LM.merge_close(np.array([5.0, 5.0 + 4e-9, 5.0 + 8e-9])).tolist() == [5.0, 5.0 + 8e-9] or \
        len(LM.merge_close(np.array([5.0, 5.0 + 4e-9, 5.0 + 8e-9]))) <= 2


def test_resample_is_linear_inside_and_constant_outside():
    L = Line(1, "A", [1.0, 3.0], [1.0, 3.0])
    assert L.at([0.0, 1.0, 2.0, 3.0, 4.0]).tolist() == [1.0, 1.0, 2.0, 3.0, 3.0]


def test_line_validation_and_clean_xy():
    with pytest.raises(LineError):
        Line(1, "x", [0.0, 0.0], [1.0, 2.0])
    with pytest.raises(LineError):
        Line(1, "x", [], [])
    x, y, notes = LM.clean_xy([2.0, 1.0, 1.0, 3.0], [20.0, 10.0, 11.0, 30.0])
    assert x.tolist() == [1.0, 2.0, 3.0] and y.tolist() == [10.0, 20.0, 30.0]
    assert len(notes) == 2


# ------------------------------------------------------------------ T-L1 (spec 10 4.1.7)
def test_t_l1_all_operations():
    A = Line(1, "A", [0, 1, 2], [0, 1, 2])
    B = Line(2, "B", [0.5, 1.5], [10, 20])
    X, y = LM.addition([A, B])
    assert X.tolist() == [0, 0.5, 1, 1.5, 2]
    assert y.tolist() == [10, 10.5, 16, 21.5, 22]
    assert LM.average([A, B])[1].tolist() == [5, 5.25, 8, 10.75, 11]
    assert LM.subtraction([A, B])[1].tolist() == [-10, -9.5, -14, -18.5, -18]
    assert LM.linear_combination([A, B], [2, 0.5])[1].tolist() == [5, 6, 9.5, 13, 14]
    np.testing.assert_allclose(LM.srss([A, B])[1], [10, 10.01249, 15.03330, 20.05617, 20.09975], atol=5e-6)
    with pytest.raises(LineError):
        LM.linear_combination([A, B], [1.0])


def test_spec06_resampling_case_and_closed_forms():
    L1 = Line(1, "", [1, 3], [1, 3])
    L2 = Line(2, "", [2, 4], [10, 20])
    assert LM.addition([L1, L2])[1].tolist() == [11, 12, 18, 23]
    one = Line(3, "", [0, 1], [1, 1])
    assert np.allclose(LM.srss([one, one, one])[1], math.sqrt(3))
    assert np.allclose(LM.average([one, Line(4, "", [0, 1], [3, 3])])[1], 2.0)


# ------------------------------------------------------------------ BROADEN (D-LIN-01)
def test_t_b1_peak_broadening():
    L = Line(1, "", [1, 5, 9, 10, 11, 15, 20], [1, 1, 1, 5, 1, 1, 1])
    X, B, notes = LM.broaden([L], 0, 15)
    f = lambda v: float(np.interp(v, X, B))
    assert all(abs(f(v) - 5) < 1e-12 for v in (8.5, 9.0, 10.0, 11.0, 11.5))
    assert abs(f(7.65) - 1) < 1e-12 and abs(f(12.65) - 1) < 1e-12
    assert abs(f(8.075) - 3) < 1e-12 and abs(f(12.075) - 3) < 1e-12
    for g in (7.65, 8.5, 9.35, 10.35, 11.5, 12.65):
        assert np.any(np.abs(X - g) < 1e-12)
    assert notes and "EDU-15" in notes[0]          # 7 < 301 points


def test_t_b2_bridging_criteria():
    L = Line(1, "", [1, 4, 5, 6, 6.5, 7, 8, 9, 12, 20], [0.5, 2, 4, 3, 3.5, 3.0, 3.8, 2, 1, 0.5])
    X, B, _ = LM.broaden([L], 15, 0)
    f = lambda v: float(np.interp(v, X, B))
    assert abs(f(5) - 4) < 1e-12 and abs(f(5.2) - 3.8) < 1e-12 and abs(f(5.1) - 3.9) < 1e-12
    for v in (5.2, 6, 6.5, 7, 7.625, 8):
        assert abs(f(v) - 3.8) < 1e-12
    assert abs(f(9) - 2) < 1e-12 and abs(f(4) - 2) < 1e-12
    # Smooth1 = 10: the valley-depth criterion leaves it unchanged, the amplitude criterion bridges
    X2, B2, _ = LM.broaden([L], 10, 0)
    assert np.allclose(np.interp(L.x, X2, B2), L.y)
    X3, B3, _ = LM.broaden([L], 10, 0, bridge="AMPLITUDE")
    assert abs(float(np.interp(7.0, X3, B3)) - 3.8) < 1e-12


def test_t_b3_envelope_inserts_crossings():
    a = Line(1, "", [0, 2], [0, 2])
    b = Line(2, "", [0, 2], [2, 0])
    X, E = LM.envelope([a, b])
    assert np.any(np.abs(X - 1.0) < 1e-12)          # the lines cross at x = 1, E(1) = 1
    assert abs(float(np.interp(1.0, X, E)) - 1.0) < 1e-12
    X0, E0 = LM.envelope([a, b], exact=False)       # union grid only: the chord cuts the corner
    assert X0.tolist() == [0, 2] and E0.tolist() == [2, 2]


def test_broadening_matches_the_definition_on_random_data():
    rng = np.random.default_rng(1)
    for trial in range(4):
        x = np.sort(rng.uniform(0.3, 25, 40))
        y = np.abs(rng.normal(1, 0.7, 40)) + 0.1
        b = (0.05, 0.1, 0.15, 0.3)[trial]
        X, B = LM.broaden_peaks(x, y, b)
        f = np.linspace(x[0], x[-1], 700)
        ref = np.empty(len(f))
        for i, fi in enumerate(f):
            lo, hi = fi / (1 + b), fi / (1 - b)
            fp = np.concatenate([np.linspace(lo, hi, 3001), x[(x > lo) & (x < hi)]])
            ref[i] = np.interp(fp, x, y).max()
        assert np.max(np.interp(f, X, B) - ref) < 1e-10           # never above the definition
        assert np.max(np.abs(np.interp(f, X, B) - ref)) < 1e-10   # exact (corner points inserted)


def test_broaden_acs_grid_is_discrete_on_source_points():
    L = Line(1, "", [1, 5, 9, 10, 11, 15, 20], [1, 1, 1, 5, 1, 1, 1])
    X, B, _ = LM.broaden([L], 0, 15, grid="ACS")
    assert X.tolist() == L.x.tolist()
    assert B.tolist() == [1, 1, 5, 5, 5, 1, 1]       # 9 and 11 see the peak at 10 within +-15 %


def test_broaden_limits_and_identity():
    L = Line(1, "", [1, 2, 3], [1, 3, 2])
    X, B, _ = LM.broaden([L], 0, 0)
    assert X.tolist() == L.x.tolist() and B.tolist() == L.y.tolist()
    with pytest.raises(LineError):
        LM.broaden([L], 0, 100)
    with pytest.raises(LineError):
        LM.broaden([L], 101, 0)
    with pytest.raises(LineError):
        LM.broaden_peaks(np.array([-1.0, 1.0]), np.array([1.0, 2.0]), 0.1)


# ------------------------------------------------------------------ files
def test_spec_file_round_trip_is_exact(tmp_path):
    rng = np.random.default_rng(2)
    A = Line(1, "a", np.sort(rng.uniform(0, 10, 13)), rng.normal(size=13))
    B = Line(2, "b", np.sort(rng.uniform(0, 12, 9)), rng.normal(size=9))
    p, X = LM.write_spec_file(tmp_path / "x.rs", [A, B])
    x, ys = LM.read_spec_file(p, 2)
    assert np.array_equal(x, X)
    assert np.array_equal(ys[0], A.at(X)) and np.array_equal(ys[1], B.at(X))
    with pytest.raises(LineError):
        LM.read_spec_file(p, 3)                    # 1 + 2 columns only


def test_th_file_formats(tmp_path):
    (tmp_path / "a.acc").write_text("0.01\n1.0 2.0 3.0\n4.0D0\n")       # several values per line, D exponent
    t, a = LM.read_th_file(tmp_path / "a.acc", 0)
    assert t.tolist() == [0.0, 0.01, 0.02, 0.03] and a.tolist() == [1, 2, 3, 4]
    (tmp_path / "p.th").write_text("time acc\n0 1\n0.5 2\n1.0 3\n")
    t, a = LM.read_th_file(tmp_path / "p.th", 1)
    assert t.tolist() == [0, 0.5, 1.0] and a.tolist() == [1, 2, 3]
    L = Line(5, "h", [0.0, 0.1, 0.2, 0.31], [1.0, 2.0, 3.0, 4.0], "history")
    p, dt, notes = LM.write_th_file(tmp_path / "o.acc", L)
    assert dt == pytest.approx(0.1) and notes and "not constant" in notes[0]
    with pytest.raises(LineError):
        LM.read_th_file(tmp_path / "p.th", 2)


# ------------------------------------------------------------------ CRITFREQ / FRAMESEL algorithms
def test_critfreq_overshoot_and_minfilter():
    F = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    AU = np.array([1.0, 2.0, 2.5, 2.0, 1.0])
    fi = np.arange(1.0, 5.0 + 1e-9, 0.1)
    ai = np.interp(fi, F, AU) + 1.25 * np.clip(1 - np.abs(fi - 2.5) / 0.3, 0, None)   # 2.25 + 1.25 at 2.5 Hz
    ai += 1.2 * np.clip(1 - np.abs(fi - 4.5) / 0.2, 0, None)                         # 1.5 + 1.2 at 4.5 Hz
    peaks, df = LM.critfreq(F, AU, fi, ai, 20, 10)       # minfilter 10 %: peaks >= 0.9 max only
    assert df == pytest.approx(0.1)
    flagged = [p for p in peaks if p.flagged]
    assert [p.number for p in flagged] == [25]
    assert flagged[0].diff_pct == pytest.approx(100 * (3.5 - 2.5) / 2.5)     # 40 % above max(A_U(2), A_U(3))
    peaks_all, _ = LM.critfreq(F, AU, fi, ai, 20, 100)  # minfilter 100 %: every peak considered
    assert sorted(p.number for p in peaks_all if p.flagged) == [25, 45]
    assert not [p for p in LM.critfreq(F, AU, fi, ai, 60, 100)[0] if p.flagged]
    coincident = [p for p in peaks_all if p.number == 30]  # TFI = TFU at a computed frequency
    assert coincident and coincident[0].diff_pct < 1e-9 and not coincident[0].flagged


def test_framesel_extrema_above_tolerance():
    a = np.array([0.0, 1.0, 0.0, -0.5, 0.0, 0.2, 0.0, -2.0, 0.0])
    assert LM.framesel(a, 40).tolist() == [2, 8]           # |1| and |-2| >= 0.8; -0.5, 0.2 below
    assert LM.framesel(a, 0).tolist() == [1, 2, 4, 6, 8, 9]  # monotone samples (index 3, 5, 7) are not extrema
    assert LM.framesel(np.zeros(4), 10).tolist() == []


# ------------------------------------------------------------------ commands
def test_readspec_column_counting_and_names(ui, tmp_path):
    (tmp_path / "R.tfu").write_text("# header\n1 10 20 30\n2 11 21 31\n3 12 22 32\n")
    assert ui.execute("READSPEC,R.tfu,3,5")              # one number: consecutive lines 5, 6, 7
    L = lines(ui)
    assert [L[n].name for n in (5, 6, 7)] == ["R.tfu", "R.tfu[2]", "R.tfu[3]"]
    assert L[7].y.tolist() == [30, 31, 32] and L[5].kind == "spectrum"
    assert not ui.execute("READSPEC,R.tfu,4,1")          # frequency column not counted
    assert any("not counted" in e for e in errors(ui))
    assert ui.execute("READSPEC,R.tfu,1,9,10")            # extra numbers ignored with a warning
    assert 10 not in L and any("ignored" in w for w in warnings(ui))
    assert not ui.execute("READSPEC,missing.rs,1,1")


def test_line_math_commands(ui, tmp_path):
    write_xy(tmp_path / "A.txt", [0, 1, 2], [0, 1, 2])
    write_xy(tmp_path / "B.txt", [0.5, 1.5], [10, 20])
    ui.run_text("READSPEC,A.txt,1,1\nREADSPEC,B.txt,1,2\nADDITION,3,1,2\nLINECOMBIN,4,1,2,2,0.5\n"
                "SRSS,1,1,2\n")                            # destination = source: computed first
    L = lines(ui)
    assert L[3].y.tolist() == [10, 10.5, 16, 21.5, 22] and L[3].name == "Linear Combin."
    assert L[4].y.tolist() == [5, 6, 9.5, 13, 14]
    assert L[1].name == "SRSS Line" and L[1].n == 5
    assert not errors(ui)
    assert not ui.execute("ADDITION,5,1,99")              # undefined source: error (D-LIN-02)
    assert any("not defined" in e for e in errors(ui))
    assert not ui.execute("LINECOMBIN,5,1,2,2")           # odd count: a coefficient is missing
    assert not ui.execute("AVERAGE,5")
    assert not ui.execute("LBINCORS,5,1")
    assert any("D-LIN-03" in e for e in errors(ui))


def test_history_commands_and_range_warning(ui, tmp_path):
    (tmp_path / "h1.acc").write_text("0.01\n" + "\n".join(str(v) for v in [0, 1, 0, -1, 0]))
    (tmp_path / "h2.acc").write_text("0.01\n" + "\n".join(str(v) for v in [1, 1, 1]))
    ui.run_text("READTH,h1.acc,0,1\nREADTH,h2.acc,0,2\nADDITION,3,1,2\nWRITETH,o.acc,3\nREADTH,o.acc,0,4\n")
    L = lines(ui)
    assert L[1].kind == "history" and L[3].kind == "history"
    assert L[3].y.tolist() == [1, 2, 1, 0, 1]               # h2 held at 1 after its end
    assert any("different time ranges" in w for w in warnings(ui))
    assert np.array_equal(L[4].x, L[3].x) and np.array_equal(L[4].y, L[3].y)
    assert not ui.execute("READTH,h1.acc,3,9")


def test_broaden_command_uses_eduopt_and_warns_on_short_lines(ui, tmp_path):
    write_xy(tmp_path / "s.rs", [1, 5, 9, 10, 11, 15, 20], [1, 1, 1, 5, 1, 1, 1])
    ui.run_text("READSPEC,s.rs,1,1\nBROADEN,2,0,15,1\nEDUOPT,BROADENGRID,ACS\nBROADEN,3,0,15,1\n")
    L = lines(ui)
    assert L[2].name == "Envelope" and L[2].n > L[1].n
    assert L[3].x.tolist() == L[1].x.tolist()
    assert any("EDU-15" in w for w in warnings(ui))
    assert not ui.execute("BROADEN,4,0,120,1")


def test_writespec_needs_defined_lines_and_limits(ui, tmp_path):
    write_xy(tmp_path / "s.rs", [1, 2], [1, 2])
    ui.execute("READSPEC,s.rs,1,1")
    assert not ui.execute("WRITESPEC,o.rs,1,2")
    assert ui.execute("WRITESPEC,o.rs,1")
    x, ys = LM.read_spec_file(tmp_path / "o.rs", 1)
    assert x.tolist() == [1, 2] and ys[0].tolist() == [1, 2]
    assert not ui.execute("WRITETH,o.acc,7")
