"""VP-52: model checks with planted faults (requirements section 6.4)."""
from sassi.verify import run_problem


def test_vp52_model_checks(tmp_path):
    res = run_problem("VP-52", tmp_path)
    failed = [(c.quantity, c.computed, c.reference) for c in res.checks if not c.passed]
    assert res.passed, failed
