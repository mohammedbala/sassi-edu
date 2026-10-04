"""VP-30: response spectra by the exact piecewise-linear recurrence (closed forms, R2 G.1)."""
from __future__ import annotations

from sassi.verify import run_problem


def test_vp30_response_spectra(tmp_path):
    res = run_problem("VP-30", tmp_path)
    assert res.passed, "\n".join(f"{c.quantity}: {c.computed:.8g} vs {c.reference:.8g} (err {c.error:.3g})"
                                 for c in res.checks if not c.passed)
