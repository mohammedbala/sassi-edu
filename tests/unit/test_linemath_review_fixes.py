"""Regression tests of the review of the line maths (D-LIN-01 BROADEN, peak bridging).

Defect: ``bridge_peaks`` / ``_raise_between`` could insert a fill-level crossing onto an existing
abscissa when a grid value lies within rounding of the fill level A (a broadened plateau vertex
1 ulp below A); the BROADEN result then had a duplicate x and the command failed.
"""
from __future__ import annotations

import numpy as np
import pytest

from sassi.plotting import lines as LM
from sassi.plotting.lines import Line


def test_raise_between_never_duplicates_an_abscissa_near_the_fill_level():
    level = 4.0
    below = np.nextafter(level, 0.0)                      # 3.9999999999999996: 1 ulp below A
    X = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    Y = np.array([5.0, below, 1.0, below, 4.0])
    Xn, Yn, changed = LM._raise_between(X, Y, 0, 4, level, True, LM._tol(Y))
    assert changed and np.all(np.diff(Xn) > 0)
    # nothing is inserted next to the 1-ulp vertices: they are raised to the level themselves
    assert Xn.tolist() == X.tolist()
    assert Yn.tolist() == [5.0, 4.0, 4.0, 4.0, 4.0]
    # a clear crossing is still inserted exactly where the flank reaches the level
    Y2 = np.array([6.0, 2.0, 1.0, 2.0, 4.0])
    Xn, Yn, _ = LM._raise_between(X, Y2, 0, 4, level, True, LM._tol(Y2))
    assert Xn.tolist() == [0.0, 0.5, 1.0, 2.0, 3.0, 4.0]
    assert Yn.tolist() == [6.0, 4.0, 4.0, 4.0, 4.0, 4.0]


def test_broaden_fuzz_coarse_spectra_with_equal_peaks():
    """Coarse spectra with integer ordinates (many equal peaks): BROADEN never fails, the result is
    strictly increasing in x, never below the broadened envelope and never below the source."""
    rng = np.random.default_rng(2024)
    for trial in range(600):
        x = np.unique(np.round(rng.uniform(1.0, 30.0, int(rng.integers(5, 15))), 2))
        if len(x) < 3:
            continue
        y = rng.integers(1, 6, len(x)).astype(float)
        p = float(rng.choice([5.0, 10.0, 15.0, 20.0, 50.0]))
        b = float(rng.choice([0.0, 5.0, 10.0, 15.0]))
        crit = "AMPLITUDE" if trial % 2 else "VALLEY"
        src = Line(1, "", x, y)
        X, Y, _ = LM.broaden([src], p, b, bridge=crit)
        out = Line(2, "", X, Y)                           # strictly increasing abscissas
        xe, e = LM.envelope([src])
        xb, B = LM.broaden_peaks(xe, e, b / 100.0)
        f = np.linspace(x[0], x[-1], 1501)
        assert np.min(out.at(f) - np.interp(f, xb, B)) > -1e-12
        assert np.min(out.at(x) - y) > -1e-12


def test_broaden_command_reviewer_case(tmp_path):
    from sassi.plotting.state import plot_state
    from sassi.prep import Interpreter
    ui = Interpreter(cwd=tmp_path)
    plot_state(ui).auto_render = False
    x = [5.27, 6.03, 9.59, 10.11, 10.55, 12.97, 13.06, 15.94, 19.46, 21.24, 24.69]
    y = [1.0, 5.0, 1.0, 1.0, 1.0, 3.0, 5.0, 4.0, 3.0, 1.0, 4.0]
    (tmp_path / "eq.rs").write_text("\n".join(f"{a!r} {b!r}" for a, b in zip(x, y)) + "\n")
    assert ui.execute("READSPEC,eq.rs,1,1")
    assert ui.execute("BROADEN,2,15,15,1")
    B = plot_state(ui).lines[2]
    # the trough before the end value 4.0 is filled at A = 4.0 from 18.331 on, the broadened
    # vertex that interpolated to 1 ulp below A (where the duplicate abscissa used to appear)
    assert B.at([18.331, 20.0, 24.69]).tolist() == pytest.approx([4.0, 4.0, 4.0], abs=1e-12)
    assert np.all(np.diff(B.x) > 0)
