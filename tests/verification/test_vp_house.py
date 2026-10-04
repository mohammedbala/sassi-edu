"""Verification problems of the HOUSE module (sassi/verify/problems/vp_house.py)."""
from __future__ import annotations

from sassi.verify import run_problem


def _assert_passed(vp_id, tmp_path):
    res = run_problem(vp_id, tmp_path)
    failed = [c for c in res.checks if not c.passed]
    msg = "\n".join(f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
                    f"error {c.error:.3g} > {c.tolerance:.3g} {c.note}" for c in failed[:20])
    assert res.passed, msg + "\n" + "\n".join(res.notes)


def test_vp_h1(tmp_path):
    _assert_passed("VP-H1", tmp_path)


def test_vp_h2(tmp_path):
    _assert_passed("VP-H2", tmp_path)
