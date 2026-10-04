"""Final audit of the verification suite (docs/verification/AUDIT_REPORT.md, section "Verification-suite
audit"): regression tests of the VP corrections.  Each test shows that a corrected criterion now detects
the defect class its label claims (mutation tests) or that a mislabelled / stale statement is gone.

* VP-TS1: the MITC4 parallelogram equilibrium patch ran inside the triangle loop with the leaked loop
  variable ``eint``: EINT 1 was recorded three times and EINT 0 never.
* VP-H2: "G* = rho Vs^2 c(beta_s)" compared |G*| of the VP's own material object with rho Vs^2; since
  |c(beta)| = 1 it could not see the damping.  The moduli are now recovered from HOUSE's Ke with affine
  fields and compared in complex form.
* VP-04: the damping reference used sassi.core.shake.interp_log_strain, the function SOIL itself uses.
* VP-46: the BBCGEN cracking-force reference was the program's own ShearCapacities.cracking.
* VP-E1: a NaN in a deck table was skipped by the row comparison (``e > worst`` is False for NaN).
* stale or overstated texts that the verification manual prints (VP-26, VP-I1, VP-H1, VP-06, VP-38T).
* VP-TS1 / VP-54T take their reference resultants from tshell.plate_rigidities: an independent check of the
  TSHELL membrane, in-plane shear and transverse-shear rigidities against the plate formulas.
* VP-48: the brute-force criterion was 2e-3 max(y) although the brute force is exact (tightened).
"""
from __future__ import annotations

import dataclasses
import inspect

import numpy as np
import pytest

from sassi.verify import VPResult, load_all, run_problem

REG = load_all()


def _check(res: VPResult, text: str):
    found = [c for c in res.checks if text in c.quantity]
    assert found, text
    return found


def test_vpts1_records_each_eint_parallelogram_patch_once():
    res = run_problem("VP-TS1")
    assert res.passed
    for eint in (0, 1):
        rows = [c for c in res.checks if "skewed parallelograms" in c.quantity and f"EINT {eint}" in c.quantity]
        assert len(rows) == 2, (eint, [c.quantity for c in rows])


def test_vph2_detects_conjugated_excavated_soil_damping(monkeypatch):
    import sassi.modules.house as H
    base = run_problem("VP-H2")
    assert base.passed
    orig = H.material_from_layer

    def conjugated(*a, **k):
        m = orig(*a, **k)
        return dataclasses.replace(m, G=np.conj(m.G), M=np.conj(m.M))

    monkeypatch.setattr(H, "material_from_layer", conjugated)
    res = run_problem("VP-H2")
    rows = _check(res, "HOUSE Ke: sum of")
    assert len(rows) == 2 and not any(c.passed for c in rows)


def test_vp04_damping_reference_is_independent_of_soil(monkeypatch):
    import sassi.core.shake as SH
    monkeypatch.setattr(SH, "interp_log_strain", lambda s, xs, ys: np.interp(s, np.sort(np.asarray(xs, float)),
                                                                            np.asarray(ys, float)[np.argsort(xs)]))
    res = run_problem("VP-04")                      # SOIL now interpolates linearly in strain (a defect)
    rows = _check(res, "damping vs D(published effective strain)")
    assert any(not c.passed for c in rows)


def test_vp46_cracking_reference_is_written_out(monkeypatch):
    import sassi.core.panels as PN
    orig = PN.shear_capacities

    def wrong(*a, **k):
        c = orig(*a, **k)
        return dataclasses.replace(c, cracking=1.1 * c.cracking)

    monkeypatch.setattr(PN, "shear_capacities", wrong)
    res = run_problem("VP-46")
    rows = _check(res, "cracking force 3 sqrt(f'c) A_W")
    assert not rows[0].passed


def test_vpe1_row_comparison_fails_on_nan():
    from sassi.verify.problems.vp_examples import compare_rows
    r = VPResult()
    compare_rows(r, "t", [[1.0, float("nan")], [2.0, 3.0]], [[1.0, 5.0], [2.0, 3.0]])
    assert not r.passed
    r = VPResult()
    compare_rows(r, "t", [[1.0, 5.0]], [[1.0, 5.0]])
    assert r.passed


@pytest.mark.parametrize("vp, stale", [
    ("VP-26", "the X and Y criteria fail"),
    ("VP-I1", "the independent rigid-foundation average"),
    ("VP-06", "tolerance 0.5 % with convergence from above"),
])
def test_vp_descriptions_are_current(vp, stale):
    p = REG[vp]
    text = p.title + "\n" + (inspect.getdoc(p.func) or "")
    assert stale not in text


def test_module_descriptions_are_current():
    from sassi.verify.problems import vp_generation, vp_house, vp_tshell
    assert "independent element-library assembly" not in vp_house.__doc__
    assert "VP-38T fails" not in vp_tshell.__doc__ and "lead decision requested" not in vp_tshell.NAFEMS_8X8_NOTE
    assert "reported as failing" not in vp_generation.__doc__


def test_tshell_rigidities_against_hand_formulas():
    """VP-TS1 / VP-54T take their reference resultants from tshell.plate_rigidities, the constants the
    element itself uses (audit finding: the membrane stiffness had no independent check).  Here the
    element's recovered resultants under constant strain states of a distorted quad are compared with the
    plate formulas written out: N = E t/(1 - nu^2) [eps_x + nu eps_y, eps_y + nu eps_x, (1 - nu)/2 gamma_xy],
    Q = (5/6) G t gamma_z, in the element's local axes x'y' (only the frame comes from the element code)."""
    from sassi.elements import material_from_E_nu, shell, tshell
    E, nu, t = 2.0e8, 0.25, 0.3
    G = E / (2 * (1 + nu))
    mat = material_from_E_nu(E, nu, 2.5)
    xyz = np.array([[0.0, 0.0, 0.0], [2.0, 0.2, 0.0], [2.3, 1.4, 0.0], [0.1, 1.1, 0.0]])   # distorted quad
    A = shell.local_frame(xyz)[0][:2, :2]                  # rows x', y' in global x, y

    def field(f):
        u = np.zeros(24)
        for k, p in enumerate(xyz):
            u[6 * k:6 * k + 6] = f(p)
        return u

    def n_hand(eps_g):                                     # global strain tensor -> local resultants
        el = A @ eps_g @ A.T
        c = E * t / (1 - nu * nu)
        return np.array([c * (el[0, 0] + nu * el[1, 1]), c * (el[1, 1] + nu * el[0, 0]),
                         c * 0.5 * (1 - nu) * 2 * el[0, 1]])
    e, g, gz = 1e-3, 2e-3, 3e-3
    for eint in (0, 1):
        S = np.real(tshell.recovery(xyz, mat, thick=t, eint=eint))
        for u, eps in ((field(lambda p: [e * p[0], 0, 0, 0, 0, 0]), np.array([[e, 0], [0, 0]])),
                       (field(lambda p: [0.5 * g * p[1], 0.5 * g * p[0], 0, 0, 0, 0]),
                        np.array([[0, 0.5 * g], [0.5 * g, 0]]))):
            ref = n_hand(eps)
            assert np.abs(S @ u - np.concatenate([ref, np.zeros(5)])).max() < 1e-10 * np.abs(ref).max(), eint
        q = S @ field(lambda p: [0, 0, gz * p[0], 0, 0, 0])                   # w = gz x, rotations 0
        qref = 5.0 / 6.0 * G * t * (A @ np.array([gz, 0.0]))
        assert np.abs(np.abs(q[3:5]) - np.abs(qref)).max() < 1e-10 * np.abs(qref).max(), eint


def test_vp48_brute_force_criterion_is_tight():
    """VP-48: the brute force evaluates every vertex and both window ends, so it is exact; the criterion was
    2e-3 max(y) with a note claiming the brute force misses maxima (observed 1.2e-14)."""
    res = run_problem("VP-48")
    row = _check(res, "random spectrum: max |B_code - B_definition|")[0]
    assert row.passed and row.kind == "abs" and row.tolerance < 1e-8
