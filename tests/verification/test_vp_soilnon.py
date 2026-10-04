"""Verification problems of SOIL-NON (sassi/verify/problems/vp_soilnon.py): VP-SN1 linear limit vs SHAKE,
VP-SN2 Masing loop, VP-SN3 nonlinear vs equivalent-linear at moderate shaking (informative)."""
from __future__ import annotations

import pytest

from sassi.verify import run_problem

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")


def _report(res) -> str:
    lines = [f"{'ok ' if c.passed else 'FAIL'} {c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
             f"error {c.error:.3e} ({c.kind}), tolerance {c.tolerance:.1e} {c.note}" for c in res.checks]
    return "\n".join(lines + [f"note: {n}" for n in res.notes])


@pytest.mark.parametrize("vp", [pytest.param("VP-SN1", marks=pytest.mark.slow), "VP-SN2", "VP-SN3"])
def test_vp_soilnon(vp, tmp_path):
    res = run_problem(vp, tmp_path)
    assert res.passed, _report(res)
