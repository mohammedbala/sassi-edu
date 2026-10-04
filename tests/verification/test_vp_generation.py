"""VP-50 (INTGEN counts) and VP-51 (generation tools): requirements section 6.3, spec 09 section 8."""
import pytest

from sassi.verify import run_problem

#: VP-51 check whose published reference (48) was internally inconsistent; D-W2-05 set it to 32 hexahedra
EXCAV_HEXES = "EXCAV delta 0.01: hexes"


def _failed(res, skip=()):
    return [(c.quantity, c.computed, c.reference) for c in res.checks if not c.passed and c.quantity not in skip]


def test_vp50_intgen_counts(tmp_path):
    res = run_problem("VP-50", tmp_path)
    assert res.passed, _failed(res)


def test_vp51_generation_tools(tmp_path):
    """Every VP-51 check except the EXCAV hexahedron count (asserted by the test below, D-W2-05)."""
    res = run_problem("VP-51", tmp_path)
    assert any(c.quantity == EXCAV_HEXES for c in res.checks)
    assert not _failed(res, skip=(EXCAV_HEXES,)), _failed(res, skip=(EXCAV_HEXES,))


def test_vp51_excav_reference_32_hexes(tmp_path):
    res = run_problem("VP-51", tmp_path)
    assert res.passed, _failed(res)
