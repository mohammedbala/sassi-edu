"""Verification problems of the TSHELL element and the HOUSE node optimizer
(sassi/verify/problems/vp_tshell.py): VP-38T, VP-54T, VP-TS1, VP-O1."""
from __future__ import annotations

import functools

import pytest

from sassi.verify import run_problem

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")


def _assert_passed(vp_id, tmp_path):
    res = run_problem(vp_id, tmp_path)
    failed = [f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, error {c.error:.3g} > "
              f"{c.kind} tol {c.tolerance:g} {c.note}" for c in res.checks if not c.passed]
    assert res.passed, f"{vp_id} failed checks:\n" + "\n".join(failed) + "\nnotes:\n" + "\n".join(res.notes)


@functools.lru_cache(maxsize=None)
def _vp38t():
    return run_problem("VP-38T")          # temporary directory, removed afterwards (D-W2-11)


def test_vp38t_nafems_test21_and_thin_limit():
    res = _vp38t()
    failed = [f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, error {c.error:.3g} > "
              f"{c.kind} tol {c.tolerance:g}" for c in res.checks if not c.passed]
    assert res.passed, "VP-38T failed checks:\n" + "\n".join(failed)


def test_vp38t_nafems_8x8_rows_are_informative():
    """D-W3-01: the NAFEMS 8x8 rows are reported as informative comparisons (lumped mass, O(h^2) from below);
    mode 4 stays the coarsest-mesh outlier (about -5.4 % / -4.1 % for EINT 0 / 1)."""
    res = _vp38t()
    rows = {c.quantity: c for c in res.checks if "NAFEMS 8x8 mesh" in c.quantity}
    assert len(rows) == 8 and all(c.kind == "info" for c in rows.values())
    for e in (0, 1):
        c = rows[f"Test 21 mode 4 (TSHELL EINT {e}, NAFEMS 8x8 mesh, Hz)"]
        assert -0.07 < (c.computed - c.reference) / c.reference < -0.03


def test_vp54t_thshlstr_face_stresses(tmp_path):
    _assert_passed("VP-54T", tmp_path)


def test_vpts1_patch_tests_and_locking(tmp_path):
    _assert_passed("VP-TS1", tmp_path)


def test_vpo1_node_optimizer(tmp_path):
    _assert_passed("VP-O1", tmp_path)
