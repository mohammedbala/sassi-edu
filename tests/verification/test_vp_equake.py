"""VP-36: RG 1.60 spectrum-compatible motions meeting SRP 3.7.1 (EQUAKE)."""
from __future__ import annotations

import pytest

from sassi.verify import run_problem


@pytest.mark.slow
def test_vp36_rg160_srp(tmp_path):
    res = run_problem("VP-36", tmp_path)
    assert res.passed, "\n".join(f"{c.quantity}: {c.computed:.6g} vs {c.reference:.6g}"
                                 for c in res.checks if not c.passed)
