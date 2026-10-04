"""Verification problems of incoherent seismic input, wave passage and multiple excitation
(sassi/verify/problems/vp_incoherency.py): VP-26, VP-27, VP-I1, VP-I2 (requirements 6.3)."""
from __future__ import annotations

import pytest

from sassi.verify import run_problem

#: VP-26 (D-W3-05): FFL and FFM differ by the antisymmetric commutator T^T [S, X] T, so the spec's
#: "identical results" holds exactly only for the input-direction generalised force, K_G-decoupled DOFs
#: (vertical translation) and coherent factors; those are criteria, together with the exact commutator
#: identity of the FFL - FFM difference.  The horizontal ATF differences are informative.


def _report(res) -> str:
    lines = [f"{'ok ' if c.passed else 'FAIL'} {c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
             f"error {c.error:.3e} ({c.kind}), tolerance {c.tolerance:.1e} {c.note}" for c in res.checks]
    return "\n".join(lines + [f"note: {n}" for n in res.notes])


@pytest.mark.p1
@pytest.mark.parametrize("vp", ["VP-26", "VP-27", "VP-I1", "VP-I2"])
def test_vp_incoherency(vp, tmp_path):
    res = run_problem(vp, tmp_path)
    assert res.passed, _report(res)
