"""Verification problem of the ANSYS interfaces (sassi/verify/problems/vp_ansys.py): VP-A1."""
from __future__ import annotations

from sassi.verify import run_problem


def _assert_passed(vp_id, tmp_path):
    res = run_problem(vp_id, tmp_path)
    failed = [c for c in res.checks if not c.passed]
    msg = "\n".join(f"{c.quantity}: computed {c.computed:.10g}, reference {c.reference:.10g}, "
                    f"error {c.error:.3g} > {c.tolerance:.3g} {c.note}" for c in failed[:20])
    assert res.passed, msg + "\n" + "\n".join(res.notes)


def test_vpa1_beam188_cantilever(tmp_path):
    """BEAM188 cantilever (.cdb) = native model (1e-8) = Euler-Bernoulli (1 %); export round trip."""
    _assert_passed("VP-A1", tmp_path)
