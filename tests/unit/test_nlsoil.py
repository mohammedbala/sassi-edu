"""Unit tests of the non-linear soil SSI core (sassi.core.nlsoil): .pin, .liq, FILE74/FILE78, curves,
strain measures, COMB_XYZ_STRAIN (UT-20), convergence (D-NLS-02) and the HOUSE iteration plan
(requirements 4.4 item 7, 4.10 item 5)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.core import nlsoil as NL
from sassi.core import shake as SH
from sassi.elements.base import material_from_M

SAND = SH.DynamicProperty("Sand", np.array([1e-4, 1e-3, 1e-2, 0.1, 1.0]), np.array([1.0, 0.99, 0.85, 0.37, 0.08]),
                          np.array([1e-4, 1e-3, 1e-2, 0.1, 1.0]), np.array([0.24, 0.8, 2.8, 9.8, 21.0]))
PIN_EXAMPLE = """1, 0.60, 2
2, 5, 180, 1
1.0, 1.0, 1
1.0, 1.0, 1
1.0, 1.0, 1
1.0, 1.0, 1
1.0, 1.0, 1
"""


# ======================================================================================
# UT-20 and the SRSS combination
# ======================================================================================
def test_ut20_srss_of_directional_strains():
    """UT-20: 0.10, 0.08, 0.02 % -> 0.1296 % (requirements 6.2, D-NLS-04)."""
    g = NL.srss(0.10, 0.08, 0.02)
    assert round(float(g), 4) == 0.1296
    assert float(g) == pytest.approx(math.sqrt(0.0168), rel=1e-15)
    np.testing.assert_allclose(NL.srss([0.1, 0.3], [0.08, 0.4], [0.02, 0.0]), [0.12961481396815722, 0.5], rtol=1e-15)
    assert float(NL.srss(0.05)) == 0.05
    with pytest.raises(ValueError):
        NL.srss()


def _f74(direction, gammas, iteration=3, esf=0.65, gmax=1.0e5):
    nl = NL.NLFile("FILE74", iteration, esf, 1, direction=direction, groups=[NL.NLGroupInfo(2, "SOLID", len(gammas), 0)])
    for k, g in enumerate(gammas, start=1):
        nl.rows.append(NL.NLRow(2, k, g, 0.0, 0.0, 1, gmax, g / esf))
    return nl


def test_combine_files_ut20_and_curve_update():
    fx, fy, fz = _f74("X", [0.10, 0.20]), _f74("Y", [0.08, 0.0]), _f74("Z", [0.02, 0.0])
    out, warns = NL.combine_files([fx, fy, fz], {1: SAND}, ["X", "Y", "Z"])
    assert not warns
    assert out.kind == "FILE74" and out.iteration == 3 and out.direction == "SRSS X+Y+Z"
    assert round(out.rows[0].gamma_eff, 4) == 0.1296
    assert out.rows[1].gamma_eff == pytest.approx(0.20)
    G, b = NL.strain_compatible(SAND, out.rows[0].gamma_eff, 1.0e5)
    assert out.rows[0].G == pytest.approx(float(G[0])) and out.rows[0].beta == pytest.approx(float(b[0]))
    # without curves: strains only and a warning
    out2, warns2 = NL.combine_files([fx, fy], None)
    assert out2.rows[0].G == 0.0 and any("FILE73" in w for w in warns2)
    # different element sets, iterations, ESF or HOUSE runs (FILE4 hash) -> errors: a stale directional
    # file must not be mixed into the iteration (review of the nonlinear_soil package)
    with pytest.raises(NL.NLSoilError, match="different elements"):
        NL.combine_files([fx, _f74("Y", [0.1])], {1: SAND})
    with pytest.raises(NL.NLSoilError, match="different iterations"):
        NL.combine_files([fx, _f74("Y", [0.1, 0.1], iteration=2)], {1: SAND})
    with pytest.raises(NL.NLSoilError, match="different ESF"):
        NL.combine_files([fx, _f74("Y", [0.1, 0.1], esf=0.6)], {1: SAND})
    fx.file4, fy.file4 = "h1", "h2"
    with pytest.raises(NL.NLSoilError, match="different HOUSE runs"):
        NL.combine_files([fx, fy], {1: SAND})
    fy.file4 = "h1"
    fz.file4 = ""                                                          # a hand-made file carries no hash
    out3, _ = NL.combine_files([fx, fy, fz], {1: SAND})
    assert out3.file4 == "h1"


# ======================================================================================
# .pin
# ======================================================================================
def test_parse_pin_manual_example():
    pin = NL.parse_pin(PIN_EXAMPLE)
    assert (pin.ngrp, pin.esf, pin.ncurv) == (1, 0.60, 2)
    g = pin.group(2)
    assert (g.nmat, g.nelem, g.istr) == (5, 180, 1)
    assert all((m.gfac, m.dfac, m.icurve) == (1.0, 1.0, 1) for m in g.materials)
    assert pin.group(7) is None


def test_pin_comments_blanks_and_round_trip():
    text = "* comment\n# another\n2 0.65 3   ! two groups\n4 2 10 0\n0.5, 2.0, 2\n1, 1, 3\n\n5,1,6,1\n0.8 1.5 1 ! tail\n"
    pin = NL.parse_pin(text)
    assert [g.igrp for g in pin.groups] == [4, 5]
    assert pin.group(4).materials[0].gfac == 0.5 and pin.group(4).materials[1].icurve == 3
    again = NL.parse_pin(NL.format_pin(pin, ["header"], {4: ["material 1", "material 2"]}))
    assert again == pin
    out = NL.format_pin(pin)
    assert out.splitlines()[0] == "2, 0.65, 3" and "0.5, 2.0, 2" in out


@pytest.mark.parametrize("text,msg", [
    ("", "empty"),
    ("1, 0.6\n", "3 items"),
    ("0, 0.6, 1\n", "NGRP"),
    ("1, 1.5, 1\n2,1,1,0\n1,1,1\n", "ESF"),
    ("1, 0.6, 1\n2,1,1,2\n1,1,1\n", "ISTR"),
    ("1, 0.6, 1\n2,1,1,0\n1,1,2\n", "ICURVE"),
    ("1, 0.6, 1\n2,2,1,0\n1,1,1\n", "only 1 lines"),
    ("2, 0.6, 1\n2,1,1,0\n1,1,1\n", "2 groups announced"),
    ("1, 0.6, 1\n2,1,1,0\n1,1,1\n9,9,9\n", "unexpected data"),
    ("1, 0.6, 1\n2,1,1,0\n0,1,1\n", "GFAC"),
    ("1, 0.6, 1\n2,1,1,0\n1,-1,1\n", "DFAC"),
    ("1, 0.6, 1\n2,1,1,0\n1,x,1\n", "not a number"),
    ("2, 0.6, 1\n2,1,1,0\n1,1,1\n2,1,1,0\n1,1,1\n", "twice"),
])
def test_parse_pin_errors(text, msg):
    with pytest.raises(NL.NLSoilError, match=msg):
        NL.parse_pin(text)


def test_liq_flag(tmp_path):
    p = tmp_path / "m.liq"
    assert NL.read_liq(p) == 0
    NL.write_liq(p, 1)
    assert p.read_text().strip() == "1" and NL.read_liq(p) == 1
    NL.write_liq(p, 0)
    assert NL.read_liq(p) == 0
    p.write_text("junk\n")
    assert NL.read_liq(p) == 0
    p.write_text("")
    assert NL.read_liq(p) == 0


# ======================================================================================
# FILE74 / FILE78 (D-FIL-08)
# ======================================================================================
def test_nlfile_round_trip(tmp_path):
    nl = NL.NLFile("FILE78", 4, 0.65, 2, source="FILE74 iteration 3 (SRSS X+Y+Z)", file4="abc123", model="m",
                   groups=[NL.NLGroupInfo(3, "SOLID", 2, 1)], deck="d1", pin="p1",
                   rows=[NL.NLRow(3, 1, 0.0123456789, 12345.678901234, 0.0712345, 2, 50000.0, 7.0),
                         NL.NLRow(3, 2, 1e-5, 49999.5, 0.0024, 2, 50000.0, 7.0)])
    p = NL.write_nlfile(tmp_path / "FILE78", nl, title="test")
    txt = p.read_text()
    assert txt.startswith("# SASSI-EDU FILE78 v1") and "ITERATION 4" in txt and "ESF 0.65" in txt
    assert "GROUP 3 SOLID 2 1" in txt and "# group element gamma_eff_pct G beta curve Gmax material" in txt
    back = NL.read_nlfile(p, "FILE78")
    assert (back.kind, back.iteration, back.esf, back.ncurv, back.source, back.file4, back.model) == \
           ("FILE78", 4, 0.65, 2, nl.source, "abc123", "m")
    assert (back.deck, back.pin) == ("d1", "p1") and "DECK d1" in txt and "PIN p1" in txt
    assert back.groups == nl.groups
    for a, b in zip(back.rows, nl.rows):
        assert (a.group, a.element, a.curve, int(a.extra)) == (b.group, b.element, b.curve, int(b.extra))
        for f in ("gamma_eff", "G", "beta", "gmax"):
            assert getattr(a, f) == pytest.approx(getattr(b, f), rel=1e-10)
    assert back.rows[0].ratio == pytest.approx(12345.678901234 / 50000.0)
    with pytest.raises(NL.NLSoilError, match="is a FILE78 file"):
        NL.read_nlfile(p, "FILE74")
    arr = back.arrays()
    assert arr["element"].tolist() == [1, 2] and arr["curve"].dtype == np.int64


def test_nlfile_errors(tmp_path):
    with pytest.raises(NL.NLSoilError, match="not found"):
        NL.read_nlfile(tmp_path / "FILE74")
    p = tmp_path / "FILE74"
    p.write_text("# nothing\n2 1 0.1 1 0.1 1 1\n")
    with pytest.raises(NL.NLSoilError, match="ITERATION"):
        NL.read_nlfile(p)
    p.write_text("ITERATION 0\nESF 0.6\n2 1 0.1 1 0.1 1 1\n2 1 0.2 1 0.1 1 1\n")
    with pytest.raises(NL.NLSoilError, match="twice"):
        NL.read_nlfile(p)
    p.write_text("ITERATION 0\nESF 0.6\n2 1 0.1\n")
    with pytest.raises(NL.NLSoilError, match="line 3"):
        NL.read_nlfile(p)


# ======================================================================================
# Curves and the strain measures
# ======================================================================================
def test_strain_compatible_is_the_soil_rule():
    g = np.array([1e-5, 1e-4, 3e-3, 0.05, 2.0])
    G, b = NL.strain_compatible(SAND, g, 2.0e5)
    np.testing.assert_allclose(G, 2.0e5 * SH.interp_log_strain(g, SAND.g_strain, SAND.g_ratio), rtol=1e-15)
    np.testing.assert_allclose(b, SH.interp_log_strain(g, SAND.d_strain, SAND.d_pct) / 100.0, rtol=1e-15)
    assert G[0] == 2.0e5 and G[-1] == pytest.approx(2.0e5 * 0.08)          # constant outside the table
    # log-linear between two points: at the geometric mean the value is the arithmetic mean
    Gm, bm = NL.strain_compatible(SAND, math.sqrt(1e-2 * 0.1), 1.0)
    assert Gm[0] == pytest.approx(0.5 * (0.85 + 0.37), rel=1e-12) and bm[0] == pytest.approx(0.5 * (2.8 + 9.8) / 100)


def test_unusable_curves_are_input_errors():
    """A curve without damping (or G/G_max) points must give an NLSoilError naming the curve, not a bare
    ValueError traceback (review finding: STRESS <iter> = 1 crashed in interp_log_strain)."""
    no_d = SH.DynamicProperty("Back", np.array([1e-4, 1.0]), np.array([1.0, 0.3]), np.zeros(0), np.zeros(0))
    zero_d = SH.DynamicProperty("Zero", np.array([1e-4, 1.0]), np.array([1.0, 0.3]), np.array([0.0]), np.array([5.0]))
    bad_g = SH.DynamicProperty("Big", np.array([1e-4, 1.0]), np.array([1.2, 0.3]), np.array([1e-4]), np.array([5.0]))
    assert NL.curve_problems(SAND) == []
    assert any("damping" in t for t in NL.curve_problems(no_d))
    assert any("damping" in t for t in NL.curve_problems(zero_d))
    assert any("Error 98" in t for t in NL.curve_problems(bad_g))
    for cv in (no_d, zero_d):
        with pytest.raises(NL.NLSoilError, match=f"curve {cv.label}"):
            NL.strain_compatible(cv, 0.01, 1.0e5)
        with pytest.raises(NL.NLSoilError):
            NL.update_rows([NL.NLRow(1, 1, 0.01, 0, 0, 1, 1.0e5)], {1: cv})
    with pytest.raises(NL.NLSoilError, match="curve 2 is not in FILE73"):
        NL.check_curves({1: SAND}, [1, 2])
    with pytest.raises(NL.NLSoilError, match="Back"):
        NL.check_curves({1: SAND, 2: no_d}, [1, 2])
    NL.check_curves({1: SAND, 2: no_d}, [1])
    # HOUSE (plan) refuses an unusable curve of the .pin as soon as FILE73 is there, even in iteration 0
    b = _base()
    pin = NL.PinData(0.65, 2, [NL.PinGroup(3, 1, 2, 0, [NL.PinMaterial(1.0, 1.0, 2)])])
    with pytest.raises(NL.NLSoilError, match="Back"):
        NL.plan_iteration(pin, _elements(b, mats=(7, 7)), dim=2, liq=0, curves={1: SAND, 2: no_d})


def test_digests_and_retire(tmp_path):
    pin = NL.parse_pin(PIN_EXAMPLE)
    assert NL.pin_digest(pin) == NL.pin_digest(NL.parse_pin("* c\n" + PIN_EXAMPLE.replace(",", " ")))
    pin.groups[0].materials[0].gfac = 0.9
    assert NL.pin_digest(pin) != NL.pin_digest(NL.parse_pin(PIN_EXAMPLE))
    d = tmp_path / "m.hou"
    d.write_text("SASSI-EDU HOUSE DECK v1\n* model_hash = 1\n[params]\ngravity = 9.81\n")
    h = NL.deck_digest(d)
    d.write_text("SASSI-EDU HOUSE DECK v1\n# other comment\n\n[params]\n  gravity = 9.81\n")
    assert NL.deck_digest(d) == h and NL.deck_digest(tmp_path / "none.hou") == ""
    (tmp_path / "FILE78").write_text("new")
    (tmp_path / "FILE78.prev").write_text("old")
    assert NL.retire(tmp_path, ("FILE78", "FILE74")) == ["FILE78.prev"]
    assert (tmp_path / "FILE78.prev").read_text() == "new" and not (tmp_path / "FILE78").exists()


def test_shear_strain_measures_solid():
    gam = 1.0e-3
    pure = np.zeros((4, 6))
    pure[:, 4] = [gam, -2 * gam, 0.5 * gam, 0.0]                           # simple shear g_xz
    assert NL.shear_strain_history(pure, 0).tolist() == pytest.approx([gam, 2 * gam, 0.5 * gam, 0.0])
    np.testing.assert_allclose(NL.shear_strain_history(pure, 1), math.sqrt(2.0 / 3.0) * np.abs(pure[:, 4]), rtol=1e-14)
    mix = np.array([[0.0, 0.0, 0.0, 1e-3, -3e-3, 2e-3]])
    assert NL.shear_strain_history(mix, 0)[0] == pytest.approx(3e-3)       # max component
    uni = np.array([[1e-3, 0, 0, 0, 0, 0]])                                 # uniaxial strain e_x
    assert NL.shear_strain_history(uni, 0)[0] == 0.0
    assert NL.shear_strain_history(uni, 1)[0] == pytest.approx((2.0 / 3.0) * math.sqrt(2.0) * 1e-3)
    gmax, geff = NL.effective_strain(NL.shear_strain_history(pure, 0)[:, None], 0.65)
    assert gmax[0] == pytest.approx(2 * gam) and geff[0] == pytest.approx(0.65 * 2 * gam)
    with pytest.raises(ValueError):
        NL.shear_strain_history(np.zeros((2, 3)), 0, "SOLID")
    with pytest.raises(ValueError):
        NL.shear_strain_history(np.zeros((2, 6)), 0, "BEAMS")


def test_shear_strain_measures_plane():
    e = np.array([[2e-4, -1e-4, 4e-4]])                                     # EXX EZZ GXZ
    assert NL.shear_strain_history(e, 0, "PLANE")[0] == pytest.approx(4e-4)
    assert NL.shear_strain_history(e, 1, "PLANE")[0] == pytest.approx(math.hypot(3e-4, 4e-4))


# ======================================================================================
# Convergence (D-NLS-02)
# ======================================================================================
def _rows(Gs, betas, gmax=100.0, gam=0.01):
    return [NL.NLRow(1, k + 1, gam, G, b, 1, gmax) for k, (G, b) in enumerate(zip(Gs, betas))]


def test_compare_measures_and_tolerances():
    used = _rows([100.0, 50.0], [0.02, 0.05])
    new = _rows([98.1, 50.0], [0.024, 0.05])
    c = NL.compare(used, new, 5)
    assert c.iteration == 5 and c.n == 2
    assert c.dG_pct == pytest.approx((98.1 - 100.0) / 98.1 * 100) and c.dG_elem == (1, 1)
    assert c.dbeta_pct == pytest.approx(0.4) and c.dbeta_elem == (1, 1)
    assert c.converged and "CONVERGED" in c.text()
    assert not NL.compare(used, _rows([97.9, 50.0], [0.02, 0.05]), 5).converged          # 2.1 % in G
    assert not NL.compare(used, _rows([100.0, 50.0], [0.02, 0.0551]), 5).converged      # 0.51 % in beta
    assert NL.compare(used, _rows([97.9, 50.0], [0.02, 0.05]), 5, tol_g=3.0).converged
    assert c.mean_ratio == pytest.approx(np.mean([0.981, 0.5])) and c.mean_beta == pytest.approx(0.037)
    with pytest.raises(NL.NLSoilError):
        NL.compare(used, [NL.NLRow(9, 9, 0, 1, 0, 1, 1)], 0)


def test_convergence_history_file(tmp_path):
    used = _rows([100.0], [0.02])
    for k, G in enumerate([60.0, 70.0, 69.5]):
        NL.record_convergence(tmp_path, NL.compare(used, _rows([G], [0.05]), k))
    NL.record_convergence(tmp_path, NL.compare(used, _rows([69.9], [0.05]), 2))          # replaces row 2
    rows = NL.read_convergence(tmp_path)
    assert [r["iteration"] for r in rows] == [0, 1, 2]
    assert rows[2]["max_dG_pct"] == pytest.approx((69.9 - 100.0) / 69.9 * 100, rel=1e-6)
    assert all(r["converged"] == 0 for r in rows)
    NL.reset_convergence(tmp_path)
    assert NL.read_convergence(tmp_path) == []


def test_status_needs_one_complete_iteration(tmp_path):
    assert NL.status(tmp_path) is None
    SH.write_file73(tmp_path / "FILE73", {"Sand": SAND})
    f78 = NL.NLFile("FILE78", 2, 0.65, 1, rows=[NL.NLRow(1, 1, 0.01, 85000.0, 0.028, 1, 1.0e5, 1)])
    NL.write_nlfile(tmp_path / "FILE78", f78)
    f74 = NL.NLFile("FILE74", 1, 0.65, 1, rows=[NL.NLRow(1, 1, 0.01, 0.0, 0.0, 1, 1.0e5, 0.02)])
    NL.write_nlfile(tmp_path / "FILE74", f74)
    assert NL.status(tmp_path) is None                                      # FILE74 of another iteration
    f74.iteration = 2
    NL.write_nlfile(tmp_path / "FILE74", f74)
    st = NL.status(tmp_path)
    assert st is not None and st.iteration == 2 and st.converged            # 0.01 % -> G/Gmax = 0.85 exactly
    assert abs(st.dG_pct) < 1e-9
    # the strains must be the response to the FILE78 properties (same FILE4 hash when both carry one)
    f78.file4, f74.file4 = "h-iteration-2", "h-other"
    NL.write_nlfile(tmp_path / "FILE78", f78)
    NL.write_nlfile(tmp_path / "FILE74", f74)
    assert NL.status(tmp_path) is None
    f74.file4 = f78.file4
    NL.write_nlfile(tmp_path / "FILE74", f74)
    assert NL.status(tmp_path) is not None


def test_current_status_requires_the_analysis_in_progress(tmp_path):
    """NLSSIITER's 'already converged' (review finding): only with .liq = 1, FILE78 of the FILE4 in the
    directory (FILE90 hash) and the HOUSE deck and .pin HOUSE read for FILE78."""
    from sassi.io.files import write_container
    SH.write_file73(tmp_path / "FILE73", {"Sand": SAND})
    pin = NL.PinData(0.65, 1, [NL.PinGroup(1, 1, 1, 0, [NL.PinMaterial(1.0, 1.0, 1)])])
    NL.write_pin(tmp_path / "m.pin", pin, ["a comment"])
    (tmp_path / "m.hou").write_text("SASSI-EDU HOUSE DECK v1\n* model_hash = aaa\n[params]\ngravity = 9.81\n")
    write_container(tmp_path / "FILE90", "FILE90", {}, {"int_hash": "i", "layer_hash": "l", "file4_hash": "h4"})
    f78 = NL.NLFile("FILE78", 2, 0.65, 1, file4="h4", deck=NL.deck_digest(tmp_path / "m.hou"),
                    pin=NL.pin_digest(pin), rows=[NL.NLRow(1, 1, 0.01, 85000.0, 0.028, 1, 1.0e5, 1)])
    NL.write_nlfile(tmp_path / "FILE78", f78)
    NL.write_nlfile(tmp_path / "FILE74", NL.NLFile("FILE74", 2, 0.65, 1, file4="h4",
                                                    rows=[NL.NLRow(1, 1, 0.01, 0.0, 0.0, 1, 1.0e5, 0.02)]))
    st, why = NL.current_status(tmp_path, "m")
    assert st is None and ".liq" in why
    NL.write_liq(tmp_path / "m.liq", 1)
    st, why = NL.current_status(tmp_path, "m")
    assert st is not None and st.converged and why == ""
    # the model-hash comment AFWRITE rewrites does not count, a data change does
    (tmp_path / "m.hou").write_text("SASSI-EDU HOUSE DECK v1\n* model_hash = bbb\n[params]\ngravity = 9.81\n")
    assert NL.current_status(tmp_path, "m")[0] is not None
    (tmp_path / "m.hou").write_text("SASSI-EDU HOUSE DECK v1\n[params]\ngravity = 9.80665\n")
    assert "changed" in NL.current_status(tmp_path, "m")[1]
    (tmp_path / "m.hou").write_text("SASSI-EDU HOUSE DECK v1\n[params]\ngravity = 9.81\n")
    # .pin comments do not count, its content does
    NL.write_pin(tmp_path / "m.pin", pin, ["another comment"])
    assert NL.current_status(tmp_path, "m")[0] is not None
    pin.esf = 0.6
    NL.write_pin(tmp_path / "m.pin", pin)
    assert "m.pin has changed" in NL.current_status(tmp_path, "m")[1]
    pin.esf = 0.65
    NL.write_pin(tmp_path / "m.pin", pin)
    # another FILE4 in the directory (HOUSE re-run)
    write_container(tmp_path / "FILE90", "FILE90", {}, {"int_hash": "i", "layer_hash": "l", "file4_hash": "h5"})
    assert "FILE90" in NL.current_status(tmp_path, "m")[1]
    write_container(tmp_path / "FILE90", "FILE90", {}, {"int_hash": "i", "layer_hash": "l", "file4_hash": "h4"})
    assert NL.current_status(tmp_path, "m")[0] is not None
    # NLSSIRESET retires the files of the analysis
    assert NL.retire(tmp_path) == ["FILE78.prev", "FILE74.prev"]
    assert NL.status(tmp_path) is None and NL.current_status(tmp_path, "m")[0] is None
    assert NL.retire(tmp_path) == []


# ======================================================================================
# Materials and the HOUSE plan
# ======================================================================================
def _base(beta=0.02):
    return material_from_M(3, 400.0, 200.0, 19.62, beta, beta, 9.81)        # Vp 400, Vs 200, rho 2


def test_nonlinear_material_keeps_poisson_ratio():
    b = _base()
    m = NL.nonlinear_material(b, 0.4 * b.G0, 0.08)
    assert m.G0 == pytest.approx(0.4 * b.G0) and m.M0 == pytest.approx(0.4 * b.M0)
    assert m.nu0 == pytest.approx(b.nu0, rel=1e-14) and m.rho == b.rho
    assert m.beta_s == m.beta_p == 0.08
    assert m.G == pytest.approx(m.G0 * complex(1 - 2 * 0.08 ** 2, 2 * 0.08 * math.sqrt(1 - 0.08 ** 2)))
    m2 = NL.nonlinear_material(b, b.G0, 0.05, 0.03)
    assert (m2.beta_s, m2.beta_p) == (0.05, 0.03)
    with pytest.raises(NL.NLSoilError, match="EDU-04"):
        NL.nonlinear_material(b, b.G0, 0.5)
    with pytest.raises(NL.NLSoilError):
        NL.nonlinear_material(b, 0.0, 0.05)


def _elements(base, mats=(7, 7, 8), group=3, code=1, excavated=False):
    return [NL.NLElementIn(k, group, k + 1, code, excavated, mid, base) for k, mid in enumerate(mats)]


def test_plan_initial_iteration_scales_the_free_field_layer():
    """Requirements 4.4 item 7: iteration 0 has G = GFAC G_layer, beta = DFAC beta_layer of the free field at
    the element; the element material (low-strain soil) keeps its density and Poisson's ratio and is G_max."""
    b = _base()                                                             # Vs 200, rho 2, nu 1/3, 2 %
    lay1 = material_from_M(3, 300.0, 150.0, 19.62, 0.05, 0.04, 9.81)        # strain-compatible free field
    lay2 = material_from_M(3, 360.0, 180.0, 19.62, 0.03, 0.03, 9.81)
    pin = NL.PinData(0.65, 2, [NL.PinGroup(3, 2, 3, 0, [NL.PinMaterial(0.5, 2.0, 1), NL.PinMaterial(1.0, 1.0, 2)])])
    els = _elements(b)
    for e, lay, no, dep in zip(els, (lay1, lay2, lay2), (11, 12, 12), (0.5, 1.5, 1.5)):
        e.layer, e.layer_no, e.depth = lay, no, dep
    others = [NL.NLElementIn(10, 1, 1, 1, True, 1, b)]                     # an excavated element of another group
    plan = NL.plan_iteration(pin, els + others, dim=2, liq=0)
    assert plan.iteration == 0 and plan.source == "PIN" and plan.convergence is None and not plan.warnings
    assert [e.element for e in plan.elements] == [1, 2, 3] and [e.index for e in plan.elements] == [0, 1, 2]
    e1, e2, e3 = plan.elements
    assert e1.material.G0 == pytest.approx(0.5 * lay1.G0, rel=1e-14) and e1.ref_layer == 11
    assert e1.material.beta_s == pytest.approx(2.0 * 0.04) and e1.material.beta_p == pytest.approx(2.0 * 0.05)
    assert e1.material.rho == b.rho and e1.material.nu0 == pytest.approx(b.nu0, rel=1e-13)
    assert e1.curve == 1 and e1.gfac == 0.5 and e1.ratio == pytest.approx(0.5 * 150.0 ** 2 / 200.0 ** 2)
    assert e2.material.G0 == pytest.approx(0.5 * lay2.G0) and e2.material.beta_s == pytest.approx(0.06)  # line 1
    assert e3.material.G0 == pytest.approx(lay2.G0) and e3.material.beta_s == pytest.approx(0.03)        # line 2
    assert e3.ref_layer == 12 and e3.curve == 2
    f78 = plan.file78("m", "hash", deck="d")
    assert f78.iteration == 0 and f78.groups[0].istr == 0 and [r.extra for r in f78.rows] == [7.0, 7.0, 8.0]
    assert f78.rows[0].gmax == b.G0 and f78.rows[0].gamma_eff == 0.0          # G_max = the element material
    assert f78.rows[0].G == pytest.approx(0.5 * lay1.G0) and f78.rows[0].beta == pytest.approx(0.08)
    assert (f78.file4, f78.deck, f78.pin) == ("hash", "d", NL.pin_digest(pin))
    # a free field equal to the material and GFAC = DFAC = 1: the material object itself (bit-exact zero SSI)
    els[2].layer = b
    plan = NL.plan_iteration(pin, els, dim=2, liq=0)
    assert plan.elements[2].material is b
    # a free field stiffer than G_max / GFAC: a warning (the material must hold the low-strain soil)
    els[0].layer = material_from_M(3, 600.0, 300.0, 19.62, 0.02, 0.02, 9.81)
    plan = NL.plan_iteration(NL.PinData(0.65, 2, [NL.PinGroup(3, 2, 3, 0, [NL.PinMaterial(1.0, 1.0, 1)] * 2)]), els,
                             dim=2, liq=0)
    assert any("start stiffer than the low-strain modulus" in t for t in plan.warnings)
    # above the ground surface: TOPL layer 1 (a note)
    els[0].depth = -0.5
    plan = NL.plan_iteration(pin, els, dim=2, liq=0)
    assert any("above the ground surface" in t for t in plan.notes)


def test_plan_initial_iteration_without_layer_table_scales_the_material():
    b = _base()
    pin = NL.PinData(0.65, 2, [NL.PinGroup(3, 2, 3, 0, [NL.PinMaterial(0.5, 2.0, 1), NL.PinMaterial(1.0, 1.0, 2)])])
    plan = NL.plan_iteration(pin, _elements(b), dim=2, liq=0)
    assert any("no free-field layer table" in t for t in plan.warnings)
    e1, e3 = plan.elements[0], plan.elements[2]
    assert e1.material.G0 == pytest.approx(0.5 * b.G0) and e1.material.beta_s == pytest.approx(0.04)
    assert e1.ref_layer == 0 and e3.material is b


def test_plan_later_iteration_uses_file74_and_curves():
    b = _base()
    pin = NL.PinData(0.65, 1, [NL.PinGroup(3, 1, 2, 1, [NL.PinMaterial(1.0, 1.0, 1)])])
    els = _elements(b, mats=(7, 7))
    f74 = NL.NLFile("FILE74", 2, 0.65, 1, direction="SRSS X+Y+Z",
                    rows=[NL.NLRow(3, 1, 0.01, 0, 0, 1, b.G0, 0.0154), NL.NLRow(3, 2, 0.1, 0, 0, 1, b.G0, 0.154)])
    prev = NL.NLFile("FILE78", 2, 0.65, 1, rows=[NL.NLRow(3, 1, 0.005, 0.9 * b.G0, 0.02, 1, b.G0, 7),
                                                  NL.NLRow(3, 2, 0.05, 0.37 * b.G0, 0.098, 1, b.G0, 7)])
    plan = NL.plan_iteration(pin, els, dim=2, liq=1, file74=f74, curves={1: SAND}, previous78=prev)
    assert plan.iteration == 3 and plan.source.startswith("FILE74 iteration 2")
    e1, e2 = plan.elements
    assert e1.material.G0 == pytest.approx(0.85 * b.G0) and e1.material.beta_s == pytest.approx(0.028)
    assert e2.material.G0 == pytest.approx(0.37 * b.G0) and e2.material.beta_p == pytest.approx(0.098)
    assert e1.gamma_eff == 0.01 and e1.istr == 1
    c = plan.convergence
    assert c is not None and c.iteration == 2
    assert c.dG_pct == pytest.approx((0.85 - 0.9) / 0.85 * 100) and c.dbeta_pct == pytest.approx(0.8)
    # a FILE78 of another iteration gives no convergence measure (only a note)
    prev.iteration = 1
    plan2 = NL.plan_iteration(pin, els, dim=2, liq=1, file74=f74, curves={1: SAND}, previous78=prev)
    assert plan2.convergence is None and any("no convergence measure" in t for t in plan2.notes)
    # ... nor strains that are not the response to the FILE78 properties (another FILE4 hash)
    prev.iteration, prev.file4, f74.file4 = 2, "h-a", "h-b"
    plan3 = NL.plan_iteration(pin, els, dim=2, liq=1, file74=f74, curves={1: SAND}, previous78=prev)
    assert plan3.convergence is None and any("not the response" in t for t in plan3.notes)
    f74.file4 = "h-a"
    assert NL.plan_iteration(pin, els, dim=2, liq=1, file74=f74, curves={1: SAND}, previous78=prev).convergence


@pytest.mark.parametrize("kwargs,msg", [
    (dict(elements_group=5), "has no elements"),
    (dict(code=4), "must hold SOLID elements"),
    (dict(excavated=True), "excavated soil elements"),
    (dict(nmat=2), "NMAT = 2"),
])
def test_plan_errors(kwargs, msg):
    b = _base()
    els = _elements(b, mats=(7, 7), group=kwargs.get("elements_group", 3), code=kwargs.get("code", 1),
                    excavated=kwargs.get("excavated", False))
    nmat = kwargs.get("nmat", 1)
    pin = NL.PinData(0.65, 1, [NL.PinGroup(3, nmat, 2, 0, [NL.PinMaterial(1.0, 1.0, 1)] * nmat)])
    with pytest.raises(NL.NLSoilError, match=msg):
        NL.plan_iteration(pin, els, dim=2, liq=0)


def test_plan_later_iteration_errors_and_warnings():
    b = _base()
    pin = NL.PinData(0.6, 2, [NL.PinGroup(3, 1, 5, 0, [NL.PinMaterial(1.0, 1.0, 1)])])
    els = _elements(b, mats=(7, 7))
    with pytest.raises(NL.NLSoilError, match="FILE74 missing"):
        NL.plan_iteration(pin, els, dim=2, liq=1, file74=None, curves={1: SAND})
    f74 = NL.NLFile("FILE74", 0, 0.65, 1, rows=[NL.NLRow(3, 1, 0.01, 0, 0, 1, b.G0, 0.015)])
    with pytest.raises(NL.NLSoilError, match="FILE73"):
        NL.plan_iteration(pin, els, dim=2, liq=1, file74=f74, curves=None)
    with pytest.raises(NL.NLSoilError, match="no strain for 1 nonlinear elements"):
        NL.plan_iteration(pin, els, dim=2, liq=1, file74=f74, curves={1: SAND})
    f74.rows.append(NL.NLRow(3, 2, 0.01, 0, 0, 1, b.G0, 0.015))
    plan = NL.plan_iteration(pin, els, dim=2, liq=1, file74=f74, curves={1: SAND})
    text = " ".join(plan.warnings)
    assert "NELEM = 5" in text and "ESF = 0.65" in text and "NCURV = 2" in text
    with pytest.raises(NL.NLSoilError, match="ICURVE = 3"):
        NL.plan_iteration(NL.PinData(0.6, 3, [NL.PinGroup(3, 1, 2, 0, [NL.PinMaterial(1.0, 1.0, 3)])]), els, dim=2,
                          liq=0, curves={1: SAND})


def test_plan_2d_needs_plane_elements():
    b = _base()
    pin = NL.PinData(0.65, 1, [NL.PinGroup(3, 1, 2, 1, [NL.PinMaterial(1.0, 1.0, 1)])])
    plan = NL.plan_iteration(pin, _elements(b, mats=(7, 7), code=4), dim=1, liq=0)
    assert plan.groups[0].etype == "PLANE" and any("2D" in t for t in plan.notes)
    with pytest.raises(NL.NLSoilError, match="PLANE elements in a 2D"):
        NL.plan_iteration(pin, _elements(b, mats=(7, 7), code=1), dim=1, liq=0)


def test_strain_file_and_group_summary():
    f78 = NL.NLFile("FILE78", 1, 0.5, 1, groups=[NL.NLGroupInfo(3, "SOLID", 2, 0)],
                    rows=[NL.NLRow(3, 1, 0.0, 9.0e4, 0.02, 1, 1.0e5, 7), NL.NLRow(3, 2, 0.0, 9.0e4, 0.02, 1, 1.0e5, 7)])
    f74 = NL.strain_file(f78, {(3, 1): 0.02, (3, 2): 0.2}, {1: SAND}, direction="Y")
    assert (f74.kind, f74.iteration, f74.esf, f74.direction) == ("FILE74", 1, 0.5, "Y")
    assert [r.gamma_eff for r in f74.rows] == pytest.approx([0.01, 0.1]) and [r.extra for r in f74.rows] == [0.02, 0.2]
    assert f74.rows[0].G == pytest.approx(0.85e5) and f74.rows[1].beta == pytest.approx(0.098)
    (g, n, r0, r1, r2, b0, b1, b2), = NL.group_summary(f74.rows)
    assert (g, n) == (3, 2) and r0 == pytest.approx(0.37) and r2 == pytest.approx(0.85) and b2 == pytest.approx(0.098)
