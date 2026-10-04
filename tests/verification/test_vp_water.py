"""Verification problems of the water modelling commands (sassi/verify/problems/vp_water.py): VP-54W
(REFINEMODEL and FILLPOOL properties of VP-54) and VP-W1 (impulsive mass of a rigid tank, full chain)."""
from __future__ import annotations

from sassi.verify import run_problem


def _assert_passed(vp_id, tmp_path):
    res = run_problem(vp_id, tmp_path)
    failed = [f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, error {c.error:.3g} > "
              f"{c.kind} tol {c.tolerance:g} {c.note}" for c in res.checks if not c.passed]
    assert res.passed, f"{vp_id} failed checks:\n" + "\n".join(failed) + "\nnotes:\n" + "\n".join(res.notes)
    return res


def test_vp54w_refinemodel_and_fillpool(tmp_path):
    res = _assert_passed("VP-54W", tmp_path)
    assert len(res.checks) >= 50


def test_vpw1_rigid_tank_impulsive_mass(tmp_path):
    res = _assert_passed("VP-W1", tmp_path)
    crit = [c for c in res.checks if "exact potential-flow impulsive mass" in c.quantity]
    assert len(crit) == 6 and all(c.kind == "rel" and c.tolerance == 0.03 for c in crit)     # 2 tanks x 3 Hz (D-W3-07)
    # the FE water is a little heavier than the exact potential flow (mesh), lighter than Housner's formula
    for c in crit:
        assert 0.0 < (c.computed - c.reference) / c.reference < 0.03
    info = [c for c in res.checks if "Housner (1963) tanh" in c.quantity]
    assert info and all(c.kind == "info" and c.computed < c.reference for c in info)
    # the SHELL and SOLID tanks give the same water response (same water mesh)
    by = {c.quantity.replace("solid walls", "shell walls"): c.computed for c in crit if c.quantity.startswith("solid")}
    assert all(abs(by[c.quantity] - c.computed) < 1e-4 for c in crit if c.quantity.startswith("shell"))
    # disk hygiene: the large binary module files are deleted after each variant
    total = sum(p.stat().st_size for p in tmp_path.rglob("*") if p.is_file())
    assert total < 5e6
