"""VP-E1: tutorial example 1 run through the interpreter (CHECK, AFWRITE, RUN<MODULE>) gives the same
decks and the same FILE8 (1e-12) as the model written by hand as module decks
(sassi/verify/problems/vp_examples.py)."""
from __future__ import annotations

import pytest

from sassi.verify import run_problem
from sassi.verify.problems.vp_examples import EXAMPLES_DIR


def _assert_passed(vp_id, tmp_path):
    res = run_problem(vp_id, tmp_path)
    failed = [c for c in res.checks if not c.passed]
    msg = "\n".join(f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
                    f"error {c.error:.3g} > {c.tolerance:.3g} {c.note}" for c in failed[:20])
    assert res.passed, msg + "\n" + "\n".join(res.notes)


@pytest.mark.skipif(not (EXAMPLES_DIR / "ex01_surface_stick.pre").exists(),
                    reason="the tutorial examples ship with the source tree only")
def test_vp_e1(tmp_path):
    _assert_passed("VP-E1", tmp_path)


def test_vp_e1_reports_a_missing_example(tmp_path, monkeypatch):
    """Without the source-tree examples VP-E1 fails with an explicit check, not with an exception."""
    import sassi.verify.problems.vp_examples as V
    monkeypatch.setattr(V, "EXAMPLES_DIR", tmp_path / "no_examples")
    res = V.vp_e1(tmp_path)
    assert not res.passed
    assert res.checks[0].quantity.startswith("tutorial example") and not res.checks[0].passed
