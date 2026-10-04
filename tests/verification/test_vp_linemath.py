"""Verification problems of line mathematics, broadening and file formats (sassi/verify/problems/vp_linemath.py)."""
from __future__ import annotations

import pytest

from sassi.verify import run_problem


def _assert_passed(vp_id, tmp_path):
    res = run_problem(vp_id, tmp_path)
    failed = [c for c in res.checks if not c.passed]
    msg = "\n".join(f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
                    f"error {c.error:.3g} > {c.tolerance:.3g} {c.note}" for c in failed[:20])
    assert res.passed, msg + "\n" + "\n".join(res.notes)


@pytest.mark.p1
def test_vp47(tmp_path):
    _assert_passed("VP-47", tmp_path)


@pytest.mark.p1
def test_vp48(tmp_path):
    _assert_passed("VP-48", tmp_path)


@pytest.mark.p1
def test_vp53(tmp_path):
    _assert_passed("VP-53", tmp_path)
