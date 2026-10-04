"""Unit tests of the SITE module (sassi.modules.site) and the free field (sassi.core.freefield)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.core import freefield as ff
from sassi.core import tlm
from sassi.io import decks
from sassi.io.container import Container
from sassi.io.files import SCHEMAS, read_container, validate
from sassi.modules.base import run_module

G = 9.81


def site_deck(n=10, h=2.0, vs=200.0, vp=400.0, rho=2.0, beta=0.05, hs=(800.0, 1600.0, 2.2, 0.02), nl=20,
              fnum=(2, 4, 10), df=0.5, waves=((2, 1, 1.0, 1.0, 0.0),), cl=1, cm=0, wopt=0, **params):
    d = decks.new("SITE")
    d["gravity"] = G
    d["nl"] = nl
    d["df"] = df
    d["cl"] = cl
    d["cm"] = cm
    d["wopt"] = wopt
    for k, v in params.items():
        d[k] = v
    for i in range(n):
        d.table("layers").append([i + 1, h, rho * G, vp, vs, beta, beta])
    vsh, vph, rhoh, bh = hs
    d.table("halfspace").append([n + 1, 0.0, rhoh * G, vph, vsh, bh, bh])
    for w in waves:
        d.table("waves").append(list(w))
    for f in fnum:
        d.table("freqs").append([f])
    return d


def run_site(tmp_path: Path, d, model="m"):
    decks.write(tmp_path / f"{model}.sit", d)
    rc = run_module("SITE", model, tmp_path)
    return rc, (tmp_path / f"{model}_SITE.out").read_text()


def test_site_writes_files_satisfying_schemas(tmp_path):
    rc, out = run_site(tmp_path, site_deck())
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1", "FILE1")
    f2 = read_container(tmp_path / "FILE2", "FILE2")
    assert validate(f1) == [] and validate(f2) == []
    for kind, c in (("FILE1", f1), ("FILE2", f2)):
        arrays, meta = SCHEMAS[kind]
        assert set(arrays) <= set(c.arrays)
        assert set(meta) <= set(c.meta)
    nF, nI, Nf = 3, 11, 11 + 20
    assert f2["kR"].shape == (nF, 2 * Nf) and f2["phix"].shape == (nF, Nf, 2 * Nf)
    assert f2["kL"].shape == (nF, Nf) and f2["phiy"].shape == (nF, Nf, Nf)
    assert f2["h_gen"].shape == (nF, 20)
    np.testing.assert_array_equal(f2["iface_user"], np.arange(nI))
    assert f1["U"].shape == (1, nF, nI, 3)
    np.testing.assert_allclose(f1["U"][0, :, 0, 0], 1.0)            # unit control motion
    assert f1.meta["cl"] == 1 and f1.meta["cm"] == 0
    np.testing.assert_allclose(f2["freq"], np.array([2, 4, 10]) * 0.5)
    assert "Soil layers" in out and "Generated half-space sublayers" in out
    assert "propagating Rayleigh and Love modes" in out and "free-field amplitudes" in out


def test_rigid_base_column_matches_closed_form(tmp_path):
    """Uniform layer on a rigid base, control point at the base: surface = 1/cos(k* H)."""
    H, vs, beta = 30.0, 200.0, 0.05
    d = site_deck(n=60, h=H / 60, vs=vs, beta=beta, nl=0, cl=61, fnum=(1, 3, 10), df=1 / 6)
    rc, out = run_site(tmp_path, d)
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    for q, f in enumerate(f1["freq"]):
        ex = 1 / np.cos(2 * np.pi * f * H / (vs * np.sqrt(cfactor(beta))))
        got = ff.free_field_at_nodes(f1, q, [[0, 0, 0]], [1], 0.0, 0.0, 0.0)[0, 0]
        assert abs(got - ex) / abs(ex) < 1e-4
    f2 = read_container(tmp_path / "FILE2")
    assert f2["iface_user"][-1] == -1 and f2.meta["base"] == "rigid"


def test_p_wave_column_uses_constrained_modulus(tmp_path):
    H, vs, vp, beta = 20.0, 150.0, 450.0, 0.03
    d = site_deck(n=80, h=H / 80, vs=vs, vp=vp, beta=beta, nl=0, cl=81, cm=2, waves=((3, 1, 1.0, 1.0, 0.0),),
                  fnum=(8,), df=0.5)
    rc, out = run_site(tmp_path, d)
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    ex = 1 / np.cos(2 * np.pi * 4.0 * H / (vp * np.sqrt(cfactor(beta))))
    assert abs(f1["U"][0, 0, 0, 2] - ex) / abs(ex) < 1e-4
    assert np.all(f1["U"][0, 0, :, :2] == 0)


def test_sh_family_and_control_direction(tmp_path):
    d = site_deck(waves=((4, 1, 1.0, 1.0, 0.0),), wopt=1, cm=1, cl=3)
    rc, out = run_site(tmp_path, d)
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    np.testing.assert_allclose(f1["U"][0, :, 2, 1], 1.0)
    assert np.all(f1["U"][0, :, :, [0, 2]] == 0)
    assert list(f1["wave_type"]) == [4]


def test_wave_ratio_interpolation(tmp_path):
    r = ff.wave_ratio([1, 2, 6, 10, 12], 2, 10, 0.3, 0.5)
    np.testing.assert_allclose(r, [0.3, 0.3, 0.4, 0.5, 0.5])
    d = site_deck(waves=((1, 1, 0.3, 0.5, 0.0), (2, 1, 0.7, 0.5, 0.0)), fnum=(1, 2, 6, 10, 12), freq1=2, freq2=10)
    rc, out = run_site(tmp_path, d)
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    np.testing.assert_allclose(f1["ratio"].sum(0), 1.0)
    np.testing.assert_allclose(f1["ratio"][0], [0.3, 0.3, 0.4, 0.5, 0.5])


def test_rayleigh_surface_wave_field(tmp_path):
    d = site_deck(waves=((1, 1, 1.0, 1.0, 0.0),), fnum=(10, 20))
    rc, out = run_site(tmp_path, d)
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    f2 = read_container(tmp_path / "FILE2")
    loss = tlm.material_loss(f2["layer_G"], f2["layer_M"])       # column loss (D-SIT-09 sector)
    for q in range(2):
        kR = f2["kR"][q]
        j = tlm.select_mode(kR, 1, loss)
        assert f1["k"][0, q] == kR[j]
        assert f1["U"][0, q, 0, 0] == pytest.approx(1.0)
        idx = f2["iface_user"]
        ratio = f1["U"][0, q, :, 2] / f1["U"][0, q, :, 0]
        np.testing.assert_allclose(ratio, -1j * f2["phiz"][q][idx, j] / f2["phix"][q][idx, j], rtol=1e-10)
    # x' phase of a surface wave: U(x') = U(0) exp(-i k x')
    xyz = np.array([[0.0, 0.0, 0.0], [7.0, 0.0, 0.0]])
    u = ff.free_field_at_nodes(f1, 0, xyz, [1, 1], 0.0, 0.0, 0.0)
    assert u[1, 0] == pytest.approx(u[0, 0] * np.exp(-1j * f1["k"][0, 0] * 7.0))


@pytest.mark.parametrize("wave, wopt, cm", [((1, 1, 1.0, 1.0, 0.0), 0, 0), ((5, 1, 1.0, 1.0, 0.0), 1, 1)])
def test_heavily_damped_site_keeps_surface_waves(tmp_path, wave, wopt, cm):
    """beta = 0.46 is valid input (EDU-04: beta < 0.5); material damping alone gives
    |Im k|/Re k = 0.52 > 0.5, so the propagating sector must be widened by the loss angle (D-SIT-09)."""
    d = site_deck(n=24, h=0.5, vs=100.0, vp=200.0, beta=0.46, hs=(300.0, 600.0, 2.2, 0.02), fnum=(20, 40),
                  waves=(wave,), wopt=wopt, cm=cm)
    rc, out = run_site(tmp_path, d)
    assert rc == 0, out
    assert "surface-wave mode selection" in out
    f1 = read_container(tmp_path / "FILE1")
    f2 = read_container(tmp_path / "FILE2")
    np.testing.assert_allclose(f1["U"][0, :, 0, cm], 1.0)
    key = "kR" if wave[0] == 1 else "kL"
    for q in range(2):
        k = f1["k"][0, q]
        w = 2 * np.pi * f1["freq"][q]
        assert k.real == f2[key][q][f2[key][q].real > 0].real.max()    # the slowest (fundamental) mode
        assert 90.0 < w / k.real < 125.0 and abs(k.imag) / k.real > 0.5


def test_least_decay_ties_reported_and_resolved_to_fundamental(tmp_path):
    """Undamped stratum on a rigid base: every propagating mode has Im k = 0, the least-decay rule
    (opt 2) ties and takes the largest Re k; the listing reports the tie (D-SIT-09, D-GEN-09)."""
    d = site_deck(n=30, h=0.5, vs=100.0, vp=200.0, beta=0.0, hs=(300.0, 600.0, 2.2, 0.0), nl=0, fnum=(40,),
                  waves=((1, 2, 1.0, 1.0, 0.0),))
    rc, out = run_site(tmp_path, d)
    assert rc == 0, out
    assert "found several modes with the same |Im k|" in out
    f1 = read_container(tmp_path / "FILE1")
    f2 = read_container(tmp_path / "FILE2")
    assert f1["k"][0, 0] == f2["kR"][0][tlm.select_mode(f2["kR"][0], 1, 0.0)]


def test_free_field_at_nodes_rotation_and_phase():
    """ANALYS mapping (R1 §2.7d): phase exp(-i k x'), x' = (x-xc) cos a + (y-yc) sin a; R_z(a)."""
    U = np.zeros((2, 1, 2, 3), complex)
    U[0, 0, :, 0] = [1.0, 0.5]          # wave 1: x' motion
    U[1, 0, :, 1] = [2.0, 1.0j]         # wave 2: y' motion
    c = Container(kind="FILE1", meta={"df": 1.0, "cl": 1, "cm": 0, "wopt": 0},
                  arrays={"fnum": np.array([1]), "freq": np.array([1.0]), "depth_user": np.array([0.0, 5.0]),
                          "U": U, "k": np.array([[0.2 - 0.01j], [0.0]]), "ratio": np.array([[0.6], [0.4]]),
                          "wave_type": np.array([1, 4])})
    a, xc, yc = 30.0, 1.0, -2.0
    xyz = np.array([[4.0, 2.0, 0.0], [1.0, -2.0, -5.0]])
    out = ff.free_field_at_nodes(c, 0, xyz, [1, 2], a, xc, yc)
    ar = np.radians(a)
    for n, it in ((0, 0), (1, 1)):
        xp = (xyz[n, 0] - xc) * np.cos(ar) + (xyz[n, 1] - yc) * np.sin(ar)
        ux = 0.6 * U[0, 0, it, 0] * np.exp(-1j * (0.2 - 0.01j) * xp)
        uy = 0.4 * U[1, 0, it, 1]
        np.testing.assert_allclose(out[n], [np.cos(ar) * ux - np.sin(ar) * uy, np.sin(ar) * ux + np.cos(ar) * uy, 0],
                                   atol=1e-14)
    with pytest.raises(ValueError):
        ff.free_field_at_nodes(c, 0, xyz, [1, 3], a, xc, yc)


def test_inclined_sh_wave_matches_closed_form(tmp_path):
    """Inclined SH (theta = 30 deg) through a layer on the exact half-space (D-SIT-08):
    surface/outcrop = 1/[cos(nu1 H) + i (G1 nu1)/(G2 nu2) sin(nu1 H)]."""
    H, v1, v2, r1, r2, b1, b2, th = 20.0, 200.0, 600.0, 1.9, 2.2, 0.05, 0.02, 30.0
    d = site_deck(n=80, h=H / 80, vs=v1, vp=2 * v1, rho=r1, beta=b1, hs=(v2, 2 * v2, r2, b2), nl=20,
                  waves=((4, 1, 1.0, 1.0, th),), wopt=1, cm=1, cl=1, fnum=(4, 9), df=0.5)
    rc, out = run_site(tmp_path, d)
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    for q, f in enumerate(f1["freq"]):
        w = 2 * np.pi * f
        k = w * np.sin(np.radians(th)) / v2
        G1, G2 = r1 * v1 ** 2 * cfactor(b1), r2 * v2 ** 2 * cfactor(b2)
        nu1 = np.sqrt(r1 * w ** 2 / G1 - k ** 2)
        nu2 = np.sqrt(r2 * w ** 2 / G2 - k ** 2)
        nu2 = -nu2 if nu2.imag > 0 else nu2
        ex = 1 / (np.cos(nu1 * H) + 1j * G1 * nu1 / (G2 * nu2) * np.sin(nu1 * H))
        got = f1["x_ucp"][0, q] / f1["x_outcrop"][0, q]
        assert abs(got - ex) / abs(ex) < 1e-3
        assert f1["k"][0, q] == pytest.approx(k)
    assert "D-SIT-08" in out


@pytest.mark.parametrize("wave,cm", [(2, 0), (3, 2), (3, 0)])
def test_inclined_inplane_waves_on_homogeneous_halfspace(tmp_path, wave, cm):
    """User layers (15 m) of half-space material: the surface motion is the free-surface motion of
    the homogeneous half-space, i.e. the outcrop motion times the up-going phase exp(-i nu H) of the
    incident wave (checks the TLM B-term and k^2 terms against the exact half-space stiffness)."""
    v, rho, beta, H, f, th = 300.0, 2.0, 0.02, 15.0, 10.0, 30.0
    d = site_deck(n=60, h=H / 60, vs=v, vp=2 * v, rho=rho, beta=beta, hs=(v, 2 * v, rho, beta), nl=20,
                  waves=((wave, 1, 1.0, 1.0, th),), wopt=0, cm=cm, cl=1, fnum=(10,), df=1.0)
    rc, out = run_site(tmp_path, d)
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    w = 2 * np.pi * f
    vw = 2 * v if wave == 3 else v
    k = w * np.sin(np.radians(th)) / vw
    nu = np.sqrt(w ** 2 / (vw ** 2 * cfactor(beta)) - k ** 2)
    nu = -nu if nu.imag > 0 else nu
    assert abs(f1["x_ucp"][0, 0] / f1["x_outcrop"][0, 0] - np.exp(-1j * nu * H)) < 2e-3


def test_inclined_small_angle_matches_vertical_within_tf(tmp_path):
    """theta -> 0: the within TF (independent of the half-space model) of the inclined-wave solver
    equals the vertical-wave solution."""
    res = []
    for ang in (0.0, 1e-4):
        wd = tmp_path / f"a{ang}"
        wd.mkdir()
        d = site_deck(cl=11, waves=((2, 1, 1.0, 1.0, ang),), fnum=(3, 7))
        rc, out = run_site(wd, d)
        assert rc == 0, out
        res.append(read_container(wd / "FILE1")["U"][0, :, 0, 0])
    np.testing.assert_allclose(res[1], res[0], rtol=1e-6)


@pytest.mark.parametrize("change,message", [
    (dict(mode1=0, mode2=0), "Error 45"),
    (dict(nl=3), "Error 47"),
    (dict(waves=()), "Error 51"),
    (dict(waves=((2, 1, 0.0, 1.0, 0.0),)), "Error 54"),
    (dict(waves=((1, 1, 0.5, 0.5, 0.0), (2, 1, 0.4, 0.5, 0.0))), "Error 55"),
    (dict(waves=((2, 1, 1.0, 1.0, 400.0),)), "Error 52"),
    (dict(cm=1), "control direction cannot be y'"),
    (dict(cl=40), "control point layer"),
    (dict(fnum=(2, 2, 4)), "duplicate"),
    (dict(fnum=()), "Error 120"),
    (dict(beta=0.6), "EDU-04"),
])
def test_site_input_errors(tmp_path, change, message):
    d = site_deck(**change)
    rc, out = run_site(tmp_path, d)
    assert rc == 1
    assert message in out


def test_no_layers_is_error_46(tmp_path):
    d = site_deck(n=0)
    rc, out = run_site(tmp_path, d)
    assert rc == 1 and "Error 46" in out


def test_vertical_wave_without_control_motion_is_error(tmp_path):
    rc, out = run_site(tmp_path, site_deck(cm=2))     # vertical SV has no z' motion
    assert rc == 1 and "no motion in the control direction" in out


def test_guideline_warnings(tmp_path):
    rc, out = run_site(tmp_path, site_deck(h=20.0, n=2, nl=5, fnum=(40,)))
    assert rc == 0
    assert "G-05" in out and "D-SIT-03" in out


@pytest.mark.parametrize("wave, fnum, expect_vp", [((3, 1, 1.0, 1.0, 0.0), 100, True),    # 50 Hz: Vs and Vp fail
                                                   ((3, 1, 1.0, 1.0, 0.0), 60, False),    # 30 Hz: only Vs fails
                                                   ((2, 1, 1.0, 1.0, 0.0), 100, False)])  # SV input: Vs only
def test_one_fifth_wavelength_rule_uses_vp_for_p_waves(tmp_path, wave, fnum, expect_vp):
    """D-SIT-10 / G-05: h = 2 m passes Vs/(5h) = 20 Hz and Vp/(5h) = 40 Hz; Vp is checked for P-wave
    input only."""
    cm = 2 if wave[0] == 3 else 0
    rc, out = run_site(tmp_path, site_deck(fnum=(fnum,), waves=(wave,), cm=cm))
    assert rc == 0, out
    assert "thicker than Vs/(5 f_max)" in out
    assert ("G-05 (P-wave input, D-SIT-10)" in out and "thicker than Vp/(5 f_max)" in out) == expect_vp


def test_opmode_data_check_writes_nothing(tmp_path):
    rc, out = run_site(tmp_path, site_deck(opmode=1))
    assert rc == 0 and "Data check only" in out
    assert not (tmp_path / "FILE1").exists() and not (tmp_path / "FILE2").exists()


def test_mode2_alone_reuses_file2(tmp_path):
    rc, out = run_site(tmp_path, site_deck(waves=((1, 1, 0.4, 0.4, 0.0), (2, 1, 0.6, 0.6, 0.0))))
    assert rc == 0, out
    U1 = read_container(tmp_path / "FILE1")["U"]
    (tmp_path / "FILE1").unlink()
    rc, out = run_site(tmp_path, site_deck(mode1=0, waves=((1, 1, 0.4, 0.4, 0.0), (2, 1, 0.6, 0.6, 0.0))))
    assert rc == 0, out
    assert "existing FILE2" in out
    np.testing.assert_allclose(read_container(tmp_path / "FILE1")["U"], U1, rtol=1e-12)
    rc, out = run_site(tmp_path, site_deck(mode1=0, fnum=(2, 5)))
    assert rc == 1 and "not in FILE2" in out


@pytest.mark.parametrize("change", [dict(vs=210.0), dict(beta=0.04), dict(vp=420.0), dict(h=2.1),
                                    dict(hs=(850.0, 1600.0, 2.2, 0.02)), dict(hs=(800.0, 1600.0, 2.3, 0.02)),
                                    dict(hslaw="uniform"), dict(cmodform=1)])
def test_mode2_alone_rejects_file2_of_another_site(tmp_path, change):
    """Mode 2 alone uses FILE2 only when it was written for the same resolved site: layer properties,
    half-space, generated sublayers (law D-SIT-02) and modulus form (requirements §4.2)."""
    base = dict(n=10, h=2.0, vs=200.0, vp=400.0, beta=0.05, hs=(800.0, 1600.0, 2.2, 0.02), hslaw="geometric")
    rc, out = run_site(tmp_path, site_deck(mode2=0, **base))
    assert rc == 0, out
    rc, out = run_site(tmp_path, site_deck(mode1=0, **{**base, **change}))
    assert rc == 1
    assert "FILE2 was written for another site" in out and "re-run SITE Mode 1" in out
    rc, out = run_site(tmp_path, site_deck(mode1=0, **base))             # the same site is accepted
    assert rc == 0, out


def test_mode2_without_file2_is_error(tmp_path):
    rc, out = run_site(tmp_path, site_deck(mode1=0))
    assert rc == 1 and "FILE2 missing" in out


def test_soilmode_reads_file88(tmp_path):
    lines = ["# SASSI-EDU FILE88 v1  strain-compatible soil properties written by SOIL", "# nlayers = 10",
             "# columns: layer thick gamma_eff_pct G Vs beta_s Vp beta_p"]          # SOIL format (D-FIL-10)
    for i in range(10):
        lines.append(f"{i + 1} 2.0 0.01 7.0e4 {180 + i} 0.06 {390 + i} 0.06")
    lines.append("# halfspace (not iterated): 0.0 21.6 1600.0 800.0 0.02 0.02")
    (tmp_path / "FILE88").write_text("\n".join(lines) + "\n")
    rc, out = run_site(tmp_path, site_deck(soilmode=1))
    assert rc == 0, out
    f2 = read_container(tmp_path / "FILE2")
    np.testing.assert_allclose(f2["x_layer_vs"][:-1], 180 + np.arange(10))
    np.testing.assert_allclose(f2["x_layer_ds"][:-1], 0.06)
    np.testing.assert_allclose(f2["x_layer_vp"][:-1], 390 + np.arange(10))
    (tmp_path / "FILE88").write_text("\n".join(lines[:8]) + "\n")
    rc, out = run_site(tmp_path, site_deck(soilmode=1))
    assert rc == 1 and "EDU-07" in out


def test_halfspace_sublayers_in_file2_follow_ut19(tmp_path):
    d = site_deck(vs=300.0, hs=(1500.0, 3000.0, 2.3, 0.01), fnum=(10, 60), df=0.5, hslaw="uniform")
    rc, out = run_site(tmp_path, d)
    assert rc == 0, out
    f2 = read_container(tmp_path / "FILE2")
    np.testing.assert_allclose(f2["h_gen"].sum(1), [450.0, 75.0])
    assert np.all(f2["h_gen"] <= (1500.0 / f2["freq"] / 8)[:, None] * (1 + 1e-12))


def test_vp06_checks_every_listed_criterion(tmp_path):
    """VP-06 (requirements §6.3 as amended by D-W1-03) must check all its criteria: the cut-offs
    0 / 11.547 / 23.094 Hz (provisional 10 %, a known limitation of the variable-depth half-space),
    the five phase velocities (0.5 %) and monotone convergence of the 20 Hz phase velocities."""
    from sassi.verify import run_problem
    res = run_problem("VP-06", tmp_path)
    cut = [c for c in res.checks if c.quantity.startswith("cut-off of Love mode")]
    assert [c.reference for c in cut] == [0.0, 11.547, 23.094]
    assert (cut[0].kind, cut[0].tolerance) == ("abs", pytest.approx(0.10 * 11.547))
    assert all(c.kind == "rel" and c.tolerance == 0.10 for c in cut[1:])
    conv = [c for c in res.checks if c.quantity.startswith("monotone convergence")]
    assert len(conv) == 2
    speeds = [c for c in res.checks if " c at " in c.quantity and c.kind == "rel"]
    assert len(speeds) == 5 and all(c.passed for c in speeds)
    assert res.passed
