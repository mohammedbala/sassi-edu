"""HOUSE and STRESS with the non-linear soil SSI option (requirements 4.4 item 7, 4.10 item 5): the .pin /
.liq / FILE74 / FILE78 flow on a one-element near-field soil column in an FV excavation (the model of
VP-N1/VP-N2, sassi/verify/problems/vp_nlsoil.py), run as RUNHOUSE / RUNANALYS / RUNSTRESS do."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from sassi.core import nlsoil as NL
from sassi.core import shake as SH
from sassi.io import decks
from sassi.io.container import read_container
from sassi.io.thfile import write_history
from sassi.verify import builders as B
from sassi.verify.problems.vp_nlsoil import SAND, near_field_column, stress_deck, synthetic_rock_motion

G = 9.81
VS0, RHO, NU, B0 = 200.0, 2.0, 1.0 / 3.0, 0.02
VP0 = VS0 * math.sqrt(2.0 * (1.0 - NU) / (1.0 - 2.0 * NU))
GMAX = RHO * VS0 ** 2
DT, NFT = 0.02, 512


def _sand():
    return SH.DynamicProperty("Sand", np.array(SAND["strain"]), np.array(SAND["g"]), np.array(SAND["strain"]),
                              np.array(SAND["d"]))


def prepare(wd: Path, nlssi: int = 1):
    """SITE + POINT for a 2-layer sand column on rock, the HOUSE deck (near-field group 2), FILE73, a control
    motion and the STRESS deck; returns (house builder, near-field group, frequency set)."""
    site = B.layered_site([(1.0, VS0, VP0, RHO, B0)] * 2, (1000.0, 2000.0, 2.2, 0.01), gravity=G, nl=10)
    fs = B.FrequencySet.fourier(DT, NFT, [1, 3, 5, 8, 11, 14, 18, 22, 27, 32, 38, 44, 50, 60, 70])
    hb, gnl, _ = near_field_column(site, 2, (VP0, VS0, RHO, B0))
    B.run_soil(wd, "m", site, fs, layer=2, rad=0.9)
    d = hb.deck("m")
    d["nlssi"] = nlssi
    B.write_deck(wd, "m", d)
    SH.write_file73(wd / "FILE73", {"Sand": _sand()})
    write_history(wd / "ctrl.acc", synthetic_rock_motion(n=400, dt=DT, pga=0.3), DT)
    B.write_deck(wd, "m", stress_deck("m", "ctrl.acc", DT, NFT))
    return hb, gnl, fs


def _pin(wd, gnl, gfac=1.0, dfac=1.0, istr=0, nmat=1, nelem=2, group=None, esf=0.65):
    NL.write_pin(wd / "m.pin", NL.PinData(esf, 1, [NL.PinGroup(group or gnl, nmat, nelem, istr,
                                                                [NL.PinMaterial(gfac, dfac, 1)] * nmat)]))


def _house(wd):
    rc = B.run("HOUSE", wd, "m", check=False)
    return rc, B.listing(wd, "m", "HOUSE")


def test_house_without_pin_runs_linear_with_a_warning(tmp_path):
    prepare(tmp_path)
    rc, out = _house(tmp_path)
    assert rc == 0 and "m.pin not found" in out and "non-linear soil" in out
    assert not (tmp_path / "FILE78").exists() and not (tmp_path / "m.liq").exists()


def test_house_initial_iteration_and_zero_change_identity(tmp_path):
    _, gnl, _ = prepare(tmp_path, nlssi=0)
    assert _house(tmp_path)[0] == 0
    Ks_lin = read_container(tmp_path / "COOSK").sparse("Ks").toarray()
    # GFAC = DFAC = 1 and near-field material = free-field layers: the same matrices as the linear run
    prepare(tmp_path)
    _pin(tmp_path, gnl)
    rc, out = _house(tmp_path)
    assert rc == 0, out
    np.testing.assert_array_equal(read_container(tmp_path / "COOSK").sparse("Ks").toarray(), Ks_lin)
    assert NL.read_liq(tmp_path / "m.liq") == 1
    for t in ("Non-linear soil SSI: iteration 0", "Non-linear soil elements", "iteration 0: initial properties"):
        assert t in out
    # GFAC 0.5, DFAC 2: half the modulus, twice the damping, FILE78 of iteration 0 with the FILE4 hash
    NL.write_liq(tmp_path / "m.liq", 0)
    _pin(tmp_path, gnl, gfac=0.5, dfac=2.0)
    assert _house(tmp_path)[0] == 0
    f78 = NL.read_nlfile(tmp_path / "FILE78", "FILE78")
    assert f78.iteration == 0 and f78.source == "PIN" and [r.element for r in f78.rows] == [1, 2]
    assert all(r.G == pytest.approx(0.5 * GMAX) and r.gmax == pytest.approx(GMAX) and r.beta == pytest.approx(2 * B0)
               and r.gamma_eff == 0.0 for r in f78.rows)
    assert f78.file4 == read_container(tmp_path / "FILE90").meta["file4_hash"]
    Ks = read_container(tmp_path / "COOSK").sparse("Ks").toarray()
    assert np.max(np.abs(Ks - Ks_lin)) > 1e-3 * np.max(np.abs(Ks_lin))


@pytest.mark.parametrize("kind,msg", [("group", "has no elements"), ("excavated", "excavated soil elements"),
                                      ("nmat", "NMAT = 2"), ("syntax", ".pin")])
def test_house_input_errors(tmp_path, kind, msg):
    _, gnl, _ = prepare(tmp_path)
    if kind == "group":
        _pin(tmp_path, gnl, group=9)
    elif kind == "excavated":
        _pin(tmp_path, gnl, group=1)                                   # group 1 is the excavated soil
    elif kind == "nmat":
        _pin(tmp_path, gnl, nmat=2)
    else:
        (tmp_path / "m.pin").write_text("1, 0.6\n")
    rc, out = _house(tmp_path)
    assert rc != 0 and msg in out


def test_house_later_iteration_from_file74(tmp_path):
    _, gnl, _ = prepare(tmp_path)
    _pin(tmp_path, gnl)
    assert _house(tmp_path)[0] == 0
    # .liq = 1 but no FILE74: stop with a message
    rc, out = _house(tmp_path)
    assert rc != 0 and "FILE74 missing" in out
    f74 = NL.NLFile("FILE74", 0, 0.65, 1, direction="X",
                    rows=[NL.NLRow(gnl, 1, 0.01, 0, 0, 1, GMAX, 0.0154), NL.NLRow(gnl, 2, 0.1, 0, 0, 1, GMAX, 0.154)])
    NL.write_nlfile(tmp_path / "FILE74", f74)
    rc, out = _house(tmp_path)
    assert rc == 0, out
    f78 = NL.read_nlfile(tmp_path / "FILE78", "FILE78")
    assert f78.iteration == 1 and f78.source.startswith("FILE74 iteration 0")
    assert f78.rows[0].G == pytest.approx(0.85 * GMAX) and f78.rows[0].beta == pytest.approx(0.028)
    assert f78.rows[1].G == pytest.approx(0.37 * GMAX) and f78.rows[1].beta == pytest.approx(0.098)
    rows = NL.read_convergence(tmp_path)
    assert [r["iteration"] for r in rows] == [0] and rows[0]["converged"] == 0
    assert rows[0]["max_dG_pct"] == pytest.approx((0.37 - 1.0) / 0.37 * 100, rel=1e-6)
    assert "Non-linear soil convergence" in out and "not converged" in out


def test_stress_effective_strains_file74(tmp_path):
    _, gnl, fs = prepare(tmp_path)
    _pin(tmp_path, gnl, istr=0)
    assert _house(tmp_path)[0] == 0
    B.run_analys(tmp_path, "m", fs, gravity=G, save=1)
    rc = B.run("STRESS", tmp_path, "m", check=False)
    out = B.listing(tmp_path, "m", "STRESS")
    assert rc == 0, out[-2000:]                                        # no EOUT request: <iter> = 1 is the output
    assert "Non-linear soil: effective strains of iteration 0, input X" in out
    f74 = NL.read_nlfile(tmp_path / "FILE74", "FILE74")
    assert (f74.iteration, f74.esf, f74.direction) == (0, 0.65, "X")
    sand = _sand()
    for r in f74.rows:
        assert r.gamma_eff == pytest.approx(0.65 * r.extra, rel=1e-12)
        Gc, bc = NL.strain_compatible(sand, r.gamma_eff, GMAX)
        assert r.G == pytest.approx(Gc[0], rel=1e-9) and r.beta == pytest.approx(bc[0], rel=1e-9)
        assert 1e-4 < r.extra < 1.0                                    # a strong-motion strain in percent
    # element 2 (the lower layer) strains more: shear strain grows with depth in a uniform column
    assert f74.rows[1].extra > f74.rows[0].extra
    # ISTR 1 (octahedral): for the shear-dominated column gamma_oct ~ sqrt(2/3) gamma_xz
    f78 = NL.read_nlfile(tmp_path / "FILE78", "FILE78")
    f78.groups[0].istr = 1
    NL.write_nlfile(tmp_path / "FILE78", f78)
    assert B.run("STRESS", tmp_path, "m", check=False) == 0
    oct74 = NL.read_nlfile(tmp_path / "FILE74", "FILE74")
    for a, b in zip(oct74.rows, f74.rows):
        assert a.extra / b.extra == pytest.approx(math.sqrt(2.0 / 3.0), rel=0.03)
    # a FILE78 element that is not in FILE4 stops the run
    f78.rows[0].element = 99
    NL.write_nlfile(tmp_path / "FILE78", f78)
    assert B.run("STRESS", tmp_path, "m", check=False) != 0
    assert "not found as SOLID/PLANE elements of FILE4" in B.listing(tmp_path, "m", "STRESS")


def test_stress_iter_without_file78_warns(tmp_path):
    _, gnl, fs = prepare(tmp_path, nlssi=0)
    assert _house(tmp_path)[0] == 0
    B.run_analys(tmp_path, "m", fs, gravity=G)
    d = decks.read(tmp_path / "m.str", "STRESS")
    d.table("eout").append([gnl, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    B.write_deck(tmp_path, "m", d)
    assert B.run("STRESS", tmp_path, "m", check=False) == 0
    out = B.listing(tmp_path, "m", "STRESS")
    assert "FILE78 does not exist" in out and not (tmp_path / "FILE74").exists()


# ======================================================================================
# Review round: free-field reference of iteration 0, stale-file protection, curve errors
# ======================================================================================
def test_house_iteration0_scales_the_free_field_layer(tmp_path):
    """Requirements 4.4 item 7: G = GFAC G_layer, beta = DFAC beta_layer of the TOPL layer at the element's
    mid-depth; the M-table material (here a softer backfill, Vs 150 m/s, 3 %) is G_max of the later updates."""
    site = B.layered_site([(1.0, 180.0, 360.0, RHO, 0.03, 0.04), (1.0, 220.0, 440.0, RHO, 0.05, 0.06)],
                          (1000.0, 2000.0, 2.2, 0.01), gravity=G, nl=10)
    fs = B.FrequencySet.fourier(DT, NFT, [1, 3, 5])
    hb, gnl, _ = near_field_column(site, 2, (300.0, 150.0, RHO, 0.03))
    B.write_deck(tmp_path, "m", B.site_deck(site, fs, model="m"))
    d = hb.deck("m")
    d["nlssi"] = 1
    B.write_deck(tmp_path, "m", d)
    _pin(tmp_path, gnl, gfac=0.5, dfac=2.0)
    rc, out = _house(tmp_path)
    assert rc == 0, out
    rows = {r.element: r for r in NL.read_nlfile(tmp_path / "FILE78", "FILE78").rows}
    for e, (vs, ds) in zip((1, 2), ((180.0, 0.03), (220.0, 0.05))):
        assert rows[e].G == pytest.approx(0.5 * RHO * vs ** 2, rel=1e-12)
        assert rows[e].beta == pytest.approx(2.0 * ds, rel=1e-12)
        assert rows[e].gmax == pytest.approx(RHO * 150.0 ** 2, rel=1e-12)
    assert "G = GFAC G_layer" in out and "layer = L number of the free-field TOPL layer" in out
    # GFAC 1: the free field (Vs 180/220) is stiffer than the backfill's G_max (Vs 150): a warning
    NL.write_liq(tmp_path / "m.liq", 0)
    _pin(tmp_path, gnl)
    rc, out = _house(tmp_path)
    assert rc == 0 and "start stiffer than the low-strain modulus" in out


def test_linear_house_run_retires_a_stale_file78(tmp_path):
    """A FILE78 describes the FILE4 HOUSE wrote with it: a later linear run (<nlssi> = 0, or 1 without a .pin)
    renames it FILE78.prev so that STRESS <iter> = 1 cannot use it (review finding on EDU-42)."""
    _, gnl, fs = prepare(tmp_path)
    _pin(tmp_path, gnl)
    assert _house(tmp_path)[0] == 0 and (tmp_path / "FILE78").exists()
    (tmp_path / "m.pin").unlink()                                    # <nlssi> = 1 but no .pin: linear run
    rc, out = _house(tmp_path)
    assert rc == 0 and "m.pin not found" in out
    assert not (tmp_path / "FILE78").exists() and (tmp_path / "FILE78.prev").exists()
    assert NL.read_liq(tmp_path / "m.liq") == 1                       # the iteration flag is not touched
    B.run_analys(tmp_path, "m", fs, gravity=G, save=0)
    assert B.run("STRESS", tmp_path, "m", check=False) != 0             # no FILE78, no EOUT: Error 79
    assert "Error 79" in B.listing(tmp_path, "m", "STRESS")
    assert not (tmp_path / "FILE74").exists()


def test_stress_curve_errors_are_input_errors(tmp_path):
    """An unusable FILE73 curve stops STRESS <iter> = 1 with a message, not a ValueError traceback."""
    _, gnl, fs = prepare(tmp_path)
    _pin(tmp_path, gnl)
    assert _house(tmp_path)[0] == 0
    B.run_analys(tmp_path, "m", fs, gravity=G, save=0)
    SH.write_file73(tmp_path / "FILE73", {"Sand": SH.DynamicProperty("Sand", np.array([1e-4, 1.0]),
                                                                     np.array([1.0, 0.3]), np.zeros(0), np.zeros(0))})
    assert B.run("STRESS", tmp_path, "m", check=False) != 0
    out = B.listing(tmp_path, "m", "STRESS")
    assert "no damping point" in out and "Traceback" not in out and "ValueError" not in out
    # HOUSE (next iteration) refuses the curve as well
    rc, out = _house(tmp_path)
    assert rc != 0 and "not usable" in out and "Traceback" not in out


def test_stress_needs_the_file8_of_the_file78_properties(tmp_path):
    """<iter> = 1: FILE8 must be the response to the FILE78 properties (FILE4 hash); FILE74 carries the hash."""
    _, gnl, fs = prepare(tmp_path)
    _pin(tmp_path, gnl)
    assert _house(tmp_path)[0] == 0
    B.run_analys(tmp_path, "m", fs, gravity=G, save=1)
    assert B.run("STRESS", tmp_path, "m", check=False) == 0
    f78 = NL.read_nlfile(tmp_path / "FILE78", "FILE78")
    f74 = NL.read_nlfile(tmp_path / "FILE74", "FILE74")
    assert f74.file4 == f78.file4 == read_container(tmp_path / "FILE90").meta["file4_hash"]
    assert f78.deck == NL.deck_digest(tmp_path / "m.hou") and f78.pin == NL.pin_digest(NL.read_pin(tmp_path / "m.pin"))
    assert _house(tmp_path)[0] == 0                                   # iteration 1: new properties, new FILE4
    assert B.run("STRESS", tmp_path, "m", check=False) != 0            # ANALYS skipped
    assert "is not the response to the properties of FILE78" in B.listing(tmp_path, "m", "STRESS")
    assert NL.read_nlfile(tmp_path / "FILE74", "FILE74").iteration == 0     # the old FILE74 is untouched
    B.run_analys(tmp_path, "m", fs, gravity=G, mode=1)
    assert B.run("STRESS", tmp_path, "m", check=False) == 0
    assert NL.read_nlfile(tmp_path / "FILE74", "FILE74").iteration == 1
