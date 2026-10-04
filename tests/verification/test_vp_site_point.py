"""Verification problems of SITE / POINT (sassi/verify/problems/vp_site_point.py)."""
from __future__ import annotations

import pytest

from sassi.verify import run_problem


def _assert_passed(vp_id, tmp_path):
    res = run_problem(vp_id, tmp_path)
    failed = [c for c in res.checks if not c.passed]
    msg = "\n".join(f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
                    f"error {c.error:.3g} > {c.tolerance:.3g} {c.note}" for c in failed[:20])
    assert res.passed, msg + "\n" + "\n".join(res.notes)


def test_vp02b(tmp_path):
    _assert_passed("VP-02b", tmp_path)


def test_vp05(tmp_path):
    _assert_passed("VP-05", tmp_path)


def test_vp06(tmp_path):
    _assert_passed("VP-06", tmp_path)


def test_vp07(tmp_path):
    _assert_passed("VP-07", tmp_path)


def test_vp08(tmp_path):
    _assert_passed("VP-08", tmp_path)


def test_vp09(tmp_path):
    _assert_passed("VP-09", tmp_path)


@pytest.mark.p1
def test_vp43(tmp_path):
    _assert_passed("VP-43", tmp_path)
