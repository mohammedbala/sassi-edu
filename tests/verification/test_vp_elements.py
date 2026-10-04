"""Verification problems of the HOUSE element library (requirements 6.3): VP-03, VP-37, VP-38, VP-39.

Each test runs the registered problem through ``sassi.verify.run_problem`` and asserts that every
check (computed vs reference within tolerance) passed; the failing checks are printed.
"""
from __future__ import annotations

import pytest

from sassi.verify import run_problem

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")


def _assert_passed(vp_id, tmp_path):
    res = run_problem(vp_id, tmp_path)
    failed = [f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
              f"error {c.error:.3g} > {c.kind} tol {c.tolerance:g}" for c in res.checks if not c.passed]
    assert res.passed, f"{vp_id} failed checks:\n" + "\n".join(failed) + "\nnotes:\n" + "\n".join(res.notes)


def test_vp03_column_dispersion(tmp_path):
    _assert_passed("VP-03", tmp_path)


def test_vp37_element_closed_forms(tmp_path):
    _assert_passed("VP-37", tmp_path)


def test_vp38_nafems_free_vibration(tmp_path):
    _assert_passed("VP-38", tmp_path)


def test_vp39_spring_mass_and_general(tmp_path):
    _assert_passed("VP-39", tmp_path)
