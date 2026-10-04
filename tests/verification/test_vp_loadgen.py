"""Verification problems of LOADGEN, Option A (sassi/verify/problems/vp_loadgen.py).

Each test runs the registered problem through ``sassi.verify.run_problem`` (the same code the ``VERIFY``
command runs) and asserts that every check passes; the failure message lists every check with computed
value, reference, error and tolerance.
"""
from __future__ import annotations

import pytest

from sassi.verify import run_problem


def _report(res) -> str:
    lines = [f"{'ok ' if c.passed else 'FAIL'} {c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
             f"error {c.error:.3e} ({c.kind}), tolerance {c.tolerance:.1e} {c.note}" for c in res.checks]
    return "\n".join(lines + [f"note: {n}" for n in res.notes])


@pytest.mark.parametrize("vp", ["VP-LA1", "VP-LA2"])
def test_vp_loadgen(vp, tmp_path):
    res = run_problem(vp, tmp_path)
    assert res.passed, _report(res)
