"""Verification problems of SOIL: VP-02a (uniform layer closed form) and VP-04 (SHAKE91 sample)."""
from __future__ import annotations

import pytest

from sassi.verify import run_problem


def _report(res):
    return "\n".join(f"{c.quantity}: computed {c.computed:.6g} ref {c.reference:.6g} err {c.error:.3g} "
                     f"tol {c.tolerance:.3g}" for c in res.checks if not c.passed)


def test_vp02a_uniform_layer(tmp_path):
    res = run_problem("VP-02a", tmp_path)
    assert res.passed, _report(res)


def test_vp04_shake91(tmp_path):
    res = run_problem("VP-04", tmp_path)
    assert res.passed, _report(res)


def test_vp04_checks_the_three_emulated_sa_values(tmp_path):
    """D-W1-06: the three SA values of the requirement table are checked (2 %) on SOIL's surface motion with
    the SHAKE91 DRCTSP defect emulated (final audit: these rows were labelled 'diagnostic' although they are
    the VP-04 spectrum criterion)."""
    res = run_problem("VP-04", tmp_path)
    sa = [c for c in res.checks if c.quantity.startswith("SA at T =") and "(D-W1-06)" in c.quantity]
    assert len(sa) == 3 and all(c.passed and c.kind == "rel" and c.tolerance == 0.02 for c in sa), _report(res)
    assert not any(c.quantity.startswith("diagnostic") for c in res.checks)
