"""Verification problems of ANALYS and the SSI chain (sassi/verify/problems/vp_analys.py):
VP-01, VP-15, VP-16, VP-22, VP-23, VP-24, VP-39b, VP-40, VP-41 (requirements 6.3)."""
from __future__ import annotations

import pytest

from sassi.verify import run_problem


def _assert_passed(vp_id, tmp_path):
    res = run_problem(vp_id, tmp_path)
    failed = [c for c in res.checks if not c.passed]
    msg = "\n".join(f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
                    f"error {c.error:.3g} > {c.tolerance:.3g} {c.note}" for c in failed[:20])
    assert res.passed, msg + "\n" + "\n".join(res.notes)


@pytest.mark.parametrize("vp_id", ["VP-01", "VP-15", "VP-16", "VP-22", "VP-23", "VP-24", "VP-39b", "VP-40",
                                   "VP-41"])
def test_vp_analys(vp_id, tmp_path):
    _assert_passed(vp_id, tmp_path)
