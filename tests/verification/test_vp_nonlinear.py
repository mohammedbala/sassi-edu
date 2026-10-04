"""Verification problems of Option NON (sassi/verify/problems/vp_nonlinear.py): VP-45, VP-46, VP-NON1.

Each test runs the registered problem through ``sassi.verify.run_problem`` (the code the ``VERIFY`` command
runs) and asserts that every check passes; the failure message lists every check with computed value,
reference, error and tolerance.
"""
from __future__ import annotations

import pytest

from sassi.verify import run_problem


def _report(res) -> str:
    lines = [f"{'ok ' if c.passed else 'FAIL'} {c.quantity}: computed {c.computed:.8g}, reference {c.reference:.8g}, "
             f"error {c.error:.3e} ({c.kind}), tolerance {c.tolerance:.1e} {c.note}" for c in res.checks]
    return "\n".join(lines + [f"note: {n}" for n in res.notes])


@pytest.mark.parametrize("vp", ["VP-45", "VP-NON1"])
def test_vp_nonlinear(vp, tmp_path):
    res = run_problem(vp, tmp_path)
    assert res.passed, _report(res)


def test_vp46_full(tmp_path):
    """VP-46: British capacities and BBCGEN at 1e-6; SI path vs the British results x 4.44822162 at 1e-6;
    SI printed references (0.1 kN) within half their last digit (D-W3-08)."""
    res = run_problem("VP-46", tmp_path)
    assert res.passed, _report(res)
