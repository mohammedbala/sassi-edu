"""Verification problems of 2D SSI (POINT2, PLANE, ANALYS <dim> = 1) and SYMM half / quarter models
(sassi/verify/problems/vp_twod_symm.py): VP-42, VP-T1, VP-T2, VP-T3 (requirements 6.3)."""
from __future__ import annotations

import pytest

from sassi.verify import run_problem


def _report(res) -> str:
    lines = [f"{'ok ' if c.passed else 'FAIL'} {c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
             f"error {c.error:.3e} ({c.kind}), tolerance {c.tolerance:.1e} {c.note}" for c in res.checks]
    return "\n".join(lines + [f"note: {n}" for n in res.notes])


@pytest.mark.p1
@pytest.mark.parametrize("vp", ["VP-42", "VP-T2", "VP-T3"])
def test_vp_twod_symm(vp, tmp_path):
    res = run_problem(vp, tmp_path)
    assert res.passed, _report(res)


@pytest.mark.p1
def test_vp_t1(tmp_path):
    res = run_problem("VP-T1", tmp_path)
    assert res.passed, _report(res)
