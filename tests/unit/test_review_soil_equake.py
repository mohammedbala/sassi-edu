"""Adversarial review tests of the SOIL (SHAKE) and EQUAKE work package.

These tests probe the implementation from angles that differ from the implementer's tests:

* an *independent* SHAKE implementation based on the layer transfer matrix of the state
  vector (u, tau) instead of the up/down-going wave recursion (R1 section 6.1-6.3), compared
  with ``sassi.core.shake.run_shake`` after several equivalent-linear iterations;
* physical identities: continuity of displacement and shear stress at every interface,
  free-surface stress, within motion ``cos(k* z)`` in a uniform layer (R2 A.1), mid-height
  location of the strain, low-frequency rigid-body limit, deconvolution round trip;
* module-level semantics: SSTR option mapping (command spec 9.2.36), SACC opt 1, the
  ``gravmult`` scaling (D-SOL-07), FILE88 contents for vertical input;
* EQUAKE: PSD normalisation by Parseval/white noise, baseline-correction idempotence,
  external-history mode in British units (D-EQK-07), errors 90/92, and long-period drift
  of generated motions (SRP 3.7.1 "no baseline drift", D-EQK-03).

Tests that documented a defect found in the review were marked ``xfail(strict=True)``; the
defects (FILE88 pairing, SSAF sampling, long-period drift, missing warnings, crash above the
wavelet band) have been fixed and the markers removed, so these tests now guard the fixes.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.core import equake_lib as EL
from sassi.core import shake as SH
from sassi.core import spectra as SP
from sassi.io import decks, textfiles
from sassi.modules import soil
from sassi.modules.base import run_module

G_SI = 9.81


# =======================================================================================
# helpers
# =======================================================================================
def _col(thick, vs, rho, beta, labels=None, g=G_SI):
    return SH.SoilColumn(np.asarray(thick, float), np.asarray(rho, float) * g, np.asarray(vs, float),
                         np.asarray(beta, float), g, list(labels) if labels else [])


THICK = [6.0, 9.0, 14.0, 21.0, 0.0]
VS = [140.0, 210.0, 320.0, 480.0, 1100.0]
RHO = [1.75, 1.85, 1.95, 2.05, 2.3]
BETA = [0.06, 0.05, 0.04, 0.03, 0.01]


def _transfer_matrix_shake(freqs, thick, rho, Gc):
    """Independent SHAKE solution with the layer transfer matrix of (u, tau).

    ``Gc`` complex moduli per sublayer (last = half-space).  Starting from the free surface
    (u = 1, tau = 0), the state is propagated down each sublayer::

        u(z)   =  u0 cos kz + tau0 sin kz / (G* k)
        tau(z) = -u0 G* k sin kz + tau0 cos kz

    In the half-space u = E e^{ikz} + F e^{-ikz}, so the up-going amplitude is
    E = (u_b + tau_b/(i k G*))/2 and the outcrop motion is 2E (R1 section 6.1).
    Returns (outcrop at the base 2E, within at the base u_b, list of mid-height
    'acceleration gradients' du/dz per soil sublayer), all per unit surface motion.
    """
    w = 2 * np.pi * np.asarray(freqs, float)
    nf = len(w)
    u = np.ones(nf, complex)
    tau = np.zeros(nf, complex)
    mid = []
    for m in range(len(thick) - 1):
        k = w * np.sqrt(rho[m] / Gc[m])
        gk = Gc[m] * k
        h2 = 0.5 * thick[m]
        # mid-height
        um = u * np.cos(k * h2) + np.where(gk != 0, tau * np.sin(k * h2) / np.where(gk != 0, gk, 1), tau * h2 / Gc[m])
        tm = -u * gk * np.sin(k * h2) + tau * np.cos(k * h2)
        mid.append(tm / Gc[m])
        h = thick[m]
        un = u * np.cos(k * h) + np.where(gk != 0, tau * np.sin(k * h) / np.where(gk != 0, gk, 1), tau * h / Gc[m])
        tn = -u * gk * np.sin(k * h) + tau * np.cos(k * h)
        u, tau = un, tn
        del um
    khs = w * np.sqrt(rho[-1] / Gc[-1])
    E = 0.5 * (u + np.where(khs != 0, tau / (1j * np.where(khs != 0, khs, 1) * Gc[-1]), 0.0))
    return 2 * E, u, mid


def _independent_equivalent_linear(acc, dt, nfft, thick, rho, gmax, beta0, curves, labels, ratio, niter):
    """Independent equivalent-linear loop (outcrop input at the top of the half-space)."""
    f = np.fft.rfftfreq(nfft, dt)
    w = 2 * np.pi * f
    A = np.fft.rfft(acc, nfft)
    G = np.array(gmax, float)
    b = np.array(beta0, float)
    for _ in range(niter):
        Gc = G * cfactor(b)
        outc, _, mid = _transfer_matrix_shake(f, thick, rho, Gc)
        gmx = []
        for m in range(len(thick) - 1):
            # strain = d(displacement)/dz = -(1/w^2) d(acc)/dz
            gam = np.zeros_like(A)
            pos = w > 0
            gam[pos] = -mid[m][pos] / w[pos] ** 2 / outc[pos] * A[pos]
            gmx.append(np.max(np.abs(np.fft.irfft(gam, nfft))) * 100.0)
        geff = ratio * np.asarray(gmx)
        for m, lab in enumerate(labels[:-1]):
            if not lab:
                continue
            gs, gr, ds, dd = curves[lab]
            G[m] = gmax[m] * np.interp(np.log10(geff[m]), np.log10(gs), gr)
            b[m] = np.interp(np.log10(geff[m]), np.log10(ds), dd) / 100.0
    Gc = G * cfactor(b)
    outc, _, _ = _transfer_matrix_shake(f, thick, rho, Gc)
    surf = np.zeros_like(A)
    surf[1:] = A[1:] / outc[1:]
    surf[0] = A[0]
    return G, b, np.fft.irfft(surf, nfft)


# =======================================================================================
# SHAKE core: physics identities
# =======================================================================================
@pytest.mark.parametrize("form", [0, 1])
def test_interface_continuity_of_displacement_and_stress(form):
    """E/F recursion must satisfy u and tau = G* du/dz continuity at every interface."""
    col = _col(THICK, VS, RHO, BETA)
    f = np.linspace(0.05, 40.0, 400)
    w = 2 * np.pi * f
    vstar = SH.complex_velocity(col.gmax, col.rho, col.beta0, form)
    E, F = SH.wave_amplitudes(f, col.thick, col.rho, vstar)
    Gs = col.gmax * cfactor(col.beta0, form)
    for m in range(col.n - 1):
        k1, k2, h = w / vstar[m], w / vstar[m + 1], col.thick[m]
        u_bot = E[:, m] * np.exp(1j * k1 * h) + F[:, m] * np.exp(-1j * k1 * h)
        u_top = E[:, m + 1] + F[:, m + 1]
        t_bot = Gs[m] * 1j * k1 * (E[:, m] * np.exp(1j * k1 * h) - F[:, m] * np.exp(-1j * k1 * h))
        t_top = Gs[m + 1] * 1j * k2 * (E[:, m + 1] - F[:, m + 1])
        np.testing.assert_allclose(u_bot, u_top, rtol=1e-11)
        np.testing.assert_allclose(t_bot, t_top, rtol=1e-11)
    # free surface: zero shear stress (E1 = F1)
    np.testing.assert_allclose(E[:, 0] - F[:, 0], 0.0, atol=1e-14)


@pytest.mark.parametrize("form", [0, 1])
def test_within_motion_is_cos_kz_in_uniform_layer(form):
    """R2 A.1: within motion at depth z of a uniform layer = surface motion x cos(k* z)."""
    H, nsub = 30.0, 6
    dz = H / nsub
    col = _col([dz] * nsub + [0.0], [200.0] * nsub + [1000.0], [2.0] * nsub + [2.2],
               [0.05] * nsub + [0.01])
    f = np.linspace(0.1, 12.0, 240)
    vss = 200.0 * np.sqrt(cfactor(0.05, form))
    for j in range(1, nsub + 1):
        z = (j - 1) * dz
        H_j = SH.column_transfer(f, col, col.gmax, col.beta0, j, False, 1, False, form)   # within(z)/surface
        np.testing.assert_allclose(H_j, np.cos(2 * np.pi * f * z / vss), rtol=1e-10, atol=1e-12)


def test_mid_height_strain_location_and_closed_form():
    """Strain at mid-height: (i) of sublayer 1 equals 2 sin(k h/2)/(w V*) for E1 = F1 = 1;
    (ii) splitting a sublayer in three equal parts leaves the mid-height strain unchanged."""
    col = _col(THICK, VS, RHO, BETA)
    f = np.linspace(0.2, 30.0, 150)
    w = 2 * np.pi * f
    vstar = SH.complex_velocity(col.gmax, col.rho, col.beta0)
    E, F = SH.wave_amplitudes(f, col.thick, col.rho, vstar)
    g = SH.mid_layer_strain_factor(f, E, F, col.thick, vstar)
    k = w / vstar[0]
    np.testing.assert_allclose(g[:, 0], 2 * np.sin(k * col.thick[0] / 2) / (w * vstar[0]), rtol=1e-11)
    # split sublayer 3 (14 m) in 3
    t2 = THICK[:2] + [14.0 / 3] * 3 + THICK[3:]
    v2 = VS[:2] + [VS[2]] * 3 + VS[3:]
    r2 = RHO[:2] + [RHO[2]] * 3 + RHO[3:]
    b2 = BETA[:2] + [BETA[2]] * 3 + BETA[3:]
    c2 = _col(t2, v2, r2, b2)
    vs2 = SH.complex_velocity(c2.gmax, c2.rho, c2.beta0)
    E2, F2 = SH.wave_amplitudes(f, c2.thick, c2.rho, vs2)
    g2 = SH.mid_layer_strain_factor(f, E2, F2, c2.thick, vs2)
    np.testing.assert_allclose(g2[:, 3], g[:, 2], rtol=1e-10)


def test_low_frequency_rigid_body_limit():
    """All motion ratios tend to 1 as f -> 0 (rigid-body motion); outcrop/within differ at O(f)."""
    col = _col(THICK, VS, RHO, BETA)
    f = np.array([1e-4, 1e-5, 1e-6])
    for lo, oo in ((1, True), (3, False), (5, True), (5, False)):
        for li, oi in ((5, True), (5, False), (1, False)):
            err = np.abs(SH.column_transfer(f, col, col.gmax, col.beta0, lo, oo, li, oi) - 1.0)
            assert err[-1] < 1e-5
            if err[0] > 1e-12:
                assert err[1] / err[0] == pytest.approx(0.1, rel=0.05) or err[1] / err[0] < 0.011


def test_deconvolution_round_trip():
    """Surface motion computed from a base outcrop motion, applied back as within motion at
    the surface, must return the base outcrop motion (linear analysis, cut-off below Nyquist)."""
    col = _col(THICK, VS, RHO, BETA)
    rng = np.random.default_rng(11)
    n, dt, nfft = 1500, 0.01, 2048
    acc = rng.standard_normal(n) * np.hanning(n)
    r1 = SH.run_shake(acc, dt, nfft, col, {}, 5, True, 0.6, 0, cutoff=45.0)
    surf = r1.acc_history(1, True)
    r2 = SH.run_shake(surf, dt, nfft, col, {}, 1, False, 0.6, 0, cutoff=45.0)
    ref = np.fft.irfft(r1.input_spectrum, nfft)
    np.testing.assert_allclose(r2.acc_history(5, True), ref, atol=1e-11 * np.max(np.abs(ref)))


def test_stress_is_complex_modulus_times_strain():
    col = _col(THICK, VS, RHO, BETA)
    rng = np.random.default_rng(3)
    acc = rng.standard_normal(700) * np.hanning(700)
    res = SH.run_shake(acc, 0.01, 1024, col, {}, 5, True, 0.6, 0)
    for lay in (1, 3, 4):
        Gs = col.gmax[lay - 1] * cfactor(col.beta0[lay - 1])
        tau = np.fft.rfft(res.stress_history(lay))
        gam = np.fft.rfft(res.strain_history(lay))
        np.testing.assert_allclose(tau[1:-1], Gs * gam[1:-1], rtol=1e-9, atol=1e-12 * np.max(np.abs(tau)))
        # SASSI form: |G*| = G, so |tau(w)| = G |gamma(w)|
        np.testing.assert_allclose(np.abs(tau[1:-1]), col.gmax[lay - 1] * np.abs(gam[1:-1]), rtol=1e-9,
                                   atol=1e-12 * np.max(np.abs(tau)))


# =======================================================================================
# SHAKE core: equivalent-linear iterations vs an independent implementation
# =======================================================================================
def _curves():
    gs = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0])
    gr = np.array([1.0, 0.99, 0.96, 0.88, 0.72, 0.52, 0.31, 0.17, 0.08])
    ds = np.array([1e-4, 1e-3, 1e-2, 0.1, 1.0])
    dd = np.array([0.6, 1.1, 3.0, 9.5, 20.0])
    return {"S": (gs, gr, ds, dd)}


@pytest.mark.parametrize("niter", [1, 3, 6])
def test_iterations_match_independent_transfer_matrix_implementation(niter):
    rng = np.random.default_rng(21)
    n, dt, nfft = 1200, 0.01, 2048
    t = np.arange(n) * dt
    acc = 2.5 * rng.standard_normal(n) * np.sin(np.pi * t / t[-1]) ** 2     # m/s^2, ~0.3 g peak
    labels = ["S", "S", "S", "", ""]
    col = _col(THICK, VS, RHO, BETA, labels)
    cv = {k: SH.DynamicProperty(k, *v) for k, v in _curves().items()}
    res = SH.run_shake(acc, dt, nfft, col, cv, 5, True, 0.65, niter)
    G, b, surf = _independent_equivalent_linear(acc, dt, nfft, col.thick, col.rho, col.gmax, col.beta0,
                                                _curves(), labels, 0.65, niter)
    np.testing.assert_allclose(res.G, G, rtol=1e-9)
    np.testing.assert_allclose(res.beta, b, rtol=1e-9)
    # strong non-linearity indeed happened
    assert np.min(res.G[:3] / col.gmax[:3]) < 0.9
    np.testing.assert_allclose(res.acc_history(1, True), surf, atol=1e-9 * np.max(np.abs(surf)))


def test_flat_curves_are_a_fixed_point():
    """Curves equal to the low-strain properties: iterations change nothing, errors are 0."""
    col = _col(THICK, VS, RHO, [0.02] * 5, ["F", "F", "F", "F", ""])
    flat = SH.DynamicProperty("F", np.array([1e-4, 1.0]), np.array([1.0, 1.0]), np.array([1e-4, 1.0]),
                              np.array([2.0, 2.0]))
    rng = np.random.default_rng(1)
    acc = rng.standard_normal(800) * np.hanning(800)
    res = SH.run_shake(acc, 0.01, 1024, col, {"F": flat}, 5, True, 0.65, 4)
    assert len(res.iterations) == 4
    np.testing.assert_allclose(res.G, col.gmax, rtol=1e-14)
    np.testing.assert_allclose(res.beta, 0.02, rtol=1e-14)
    for rec in res.iterations:
        assert np.max(np.abs(rec.dG_pct)) < 1e-10 and np.max(np.abs(rec.dbeta_pct)) < 1e-10
    # effective strain is ratio x max|gamma(t)| (percent) of the final response
    gm = np.array([np.max(np.abs(res.strain_history(i + 1))) * 100 for i in range(4)])
    np.testing.assert_allclose(res.gamma_max, gm, rtol=1e-12)
    np.testing.assert_allclose(res.gamma_eff, 0.65 * gm, rtol=1e-12)


def test_rigid_base_resonances_of_undamped_layer():
    """Undamped layer on a (nearly) rigid base: surface/base-within peaks at f_n = (2n-1) Vs / 4H."""
    H, vs = 20.0, 250.0
    col = _col([H, 0.0], [vs, 1e6], [2.0, 2.0], [1e-6, 0.0])
    for n in (1, 2, 3):
        fn = (2 * n - 1) * vs / (4 * H)
        f = fn * np.array([0.97, 1.0, 1.03])
        a = np.abs(SH.column_transfer(f, col, col.gmax, col.beta0, 1, True, 2, False))
        assert a[1] > 1e3 * max(a[0], a[2]) / 1e2 and a[1] > 1e4


# =======================================================================================
# SOIL module
# =======================================================================================
def _write_motion(path: Path, n=600, seed=0, amp=0.1):
    rng = np.random.default_rng(seed)
    a = amp * rng.standard_normal(n) * np.hanning(n)
    path.write_text("\n".join(f"{v:.8e}" for v in a) + "\n")
    return a


def _soil_deck(**kw) -> decks.Deck:
    d = decks.new("SOIL")
    d.params.update(dict(title="review", model="m", nrval=600, grav=G_SI, header=0, outcrop=1, save=1, iter=2,
                         ratio=0.65, gravmult=1.0, cof=0.0, soilcutoff=0.0, delt=0.01, nft=1024, cl=3,
                         thfile="in.acc", thtit="t", mult=0.0, max=0.2, indir=0, cmodform=0))
    d.params.update(kw)
    # profile columns: layer thick weight vp vs dp ds dynprop
    d.table("profile").append([1, 8.0, 18.5, 400.0, 180.0, 0.05, 0.05, "Soil"])
    d.table("profile").append([2, 12.0, 19.5, 650.0, 300.0, 0.04, 0.04, "Soil"])
    d.table("profile").append([3, 0.0, 22.0, 2000.0, 900.0, 0.01, 0.01, ""])
    gs, gr, ds, dd = _curves()["S"]
    for x, y in zip(gs, gr):
        d.table("dynp").append(["Soil", "G", x, y])
    for x, y in zip(ds, dd):
        d.table("dynp").append(["Soil", "D", x, y])
    d.table("damp").append([0.05])
    return d


def test_sstr_option_mapping(tmp_path):
    """SSTR,<layer>,<opt1 compute stress>,<opt2 save SS>,<opt3 compute strain>,<opt4 save SN>."""
    _write_motion(tmp_path / "in.acc")
    d = _soil_deck()
    d.table("sstr").append([1, 1, 0, 1, 0])      # compute both, save none
    d.table("sstr").append([2, 0, 1, 0, 1])      # save both
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 0
    assert not (tmp_path / "SS001.TH").exists() and not (tmp_path / "SN001.TH").exists()
    assert (tmp_path / "SS002.TH").exists() and (tmp_path / "SN002.TH").exists()
    out = (tmp_path / "m_SOIL.out").read_text()
    assert "layer   1: maximum strain" in out and "layer   1: maximum stress" in out


def test_sacc_opt1_lists_maximum_without_file_and_surface_outcrop_equals_within(tmp_path):
    _write_motion(tmp_path / "in.acc")
    d = _soil_deck()
    d.table("sacc").append([2, 1, 0])
    d.table("sacc").append([1, 2, 1])
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 0
    assert not (tmp_path / "ACC002.TH").exists()
    a_out, _ = textfiles.read_history(tmp_path / "ACC001.TH")
    d.tables["sacc"].rows = [[1, 2, 0]]
    sub = tmp_path / "w"
    sub.mkdir()
    (sub / "in.acc").write_text((tmp_path / "in.acc").read_text())
    soil.write_deck(sub / "m.soi", d)
    assert run_module("SOIL", "m", sub) == 0
    a_in, _ = textfiles.read_history(sub / "ACC001.TH")
    np.testing.assert_allclose(a_in, a_out, rtol=0, atol=1e-12)        # 2E1 = E1 + F1 at a free surface


def test_gravmult_scales_rs_ordinates(tmp_path):
    """D-SOL-07: RS ordinates are SA in g multiplied by <gravmult>."""
    res = {}
    for gm in (1.0, 2.5):
        wd = tmp_path / f"g{gm}"
        wd.mkdir()
        _write_motion(wd / "in.acc")
        d = _soil_deck(gravmult=gm)
        d.table("srs").append([1, 1, 1])
        soil.write_deck(wd / "m.soi", d)
        assert run_module("SOIL", "m", wd) == 0
        res[gm] = textfiles.read_xy(wd / soil.soil_rs_name(1, 1))
    np.testing.assert_allclose(res[2.5][:, 1], 2.5 * res[1.0][:, 1], rtol=1e-7)
    # and the values are the spectrum of ACC001 in g (independent of gravity units)
    assert 0.2 < res[1.0][-1, 1] < 1.0


def test_vertical_input_file88(tmp_path):
    """indir 1: Vp/beta_p from the P-wave analysis (no iteration); Vs, beta_s copied from the input."""
    _write_motion(tmp_path / "in.acc")
    d = _soil_deck(indir=1, iter=0)
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 0
    f88 = SH.read_file88(tmp_path / "FILE88")
    np.testing.assert_allclose(f88["Vs"], [180.0, 300.0])
    np.testing.assert_allclose(f88["Vp"], [400.0, 650.0])
    np.testing.assert_allclose(f88["beta_s"], [0.05, 0.04])
    np.testing.assert_allclose(f88["beta_p"], [0.05, 0.04])
    np.testing.assert_allclose(f88["G"], np.array([18.5, 19.5]) / G_SI * f88["Vs"] ** 2, rtol=1e-8)


def test_file88_rows_are_strain_compatible(tmp_path):
    """FILE88 row: G = Gmax (G/Gmax)(gamma_eff) and beta_s = D(gamma_eff) (SHAKE91 Table B-2 pairs)."""
    _write_motion(tmp_path / "in.acc", amp=0.3)
    d = _soil_deck(iter=1, max=0.4)
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 0
    f88 = SH.read_file88(tmp_path / "FILE88")
    cv = SH.read_file73(tmp_path / "FILE73")[1]
    gmax = np.array([18.5, 19.5]) / G_SI * np.array([180.0, 300.0]) ** 2
    np.testing.assert_allclose(f88["G"], gmax * cv.g_over_gmax(f88["gamma_eff_pct"]), rtol=1e-3)
    np.testing.assert_allclose(f88["beta_s"], cv.damping(f88["gamma_eff_pct"]), rtol=1e-3)


def test_ssaf_rs_ratio_sampled_at_freqstep(tmp_path):
    _write_motion(tmp_path / "in.acc")
    d = _soil_deck(iter=0)
    d.table("ssaf").append([1, 1, 1, 1, 3, 0.5, "surface / base outcrop"])
    soil.write_deck(tmp_path / "m.soi", d)
    assert run_module("SOIL", "m", tmp_path) == 0
    rat = textfiles.read_xy(tmp_path / soil.saf_rs_name(1, True, 3, True, 1))
    k = rat[:, 0] / 0.5
    np.testing.assert_allclose(k, np.round(k), atol=1e-6)


# =======================================================================================
# EQUAKE library
# =======================================================================================
def test_psd_normalisation_white_noise():
    """One-sided S0(w) = 2|F|^2/(2 pi T_D): white noise of variance s^2 at dt has S0 = s^2 dt / pi."""
    rng = np.random.default_rng(0)
    dt, n = 0.005, 40001
    s = 37.0
    a = s * rng.standard_normal(n)
    f, p = SP.band_averaged_psd(a, dt, window=(0.0, (n - 1) * dt))
    m = (f > 1.0) & (f < 80.0)
    assert np.mean(p[m]) == pytest.approx(s * s * dt / np.pi, rel=0.02)
    # Parseval on the raw (unaveraged) definition: int_0^inf S0 dw = mean square
    F = dt * np.fft.rfft(a)
    raw = 2 * np.abs(F) ** 2 / (2 * np.pi * n * dt)
    dw = 2 * np.pi / (n * dt)
    w8 = np.ones_like(raw)
    w8[0] = 0.5
    if n % 2 == 0:
        w8[-1] = 0.5
    assert np.sum(raw * w8) * dw == pytest.approx(np.mean(a * a), rel=1e-10)


def test_app_a_target_rms_is_plausible_for_1g():
    """Integral of the App. A PSD over w gives the strong-motion mean square: ~0.3-0.4 g RMS at 1 g."""
    f = np.linspace(1e-4, 200.0, 400001)
    p = EL.rg160_target_psd(f)
    ms = float(np.sum(0.5 * (p[1:] + p[:-1]) * np.diff(2 * np.pi * f)))     # trapezoidal rule
    rms_g = np.sqrt(ms) / 980.665
    assert 0.3 < rms_g < 0.4


def test_baseline_correct_is_idempotent_and_exact_on_compatible_records():
    dt = 0.01
    t = np.arange(2001) * dt
    # displacement with compact support -> acceleration already at rest at both ends
    d = np.exp(-((t - 10.0) / 2.0) ** 2) * np.sin(2 * np.pi * 0.7 * t)
    a = np.gradient(np.gradient(d, dt), dt)
    a1 = EL.baseline_correct(a, dt)
    a2 = EL.baseline_correct(a1, dt)
    np.testing.assert_allclose(a2, a1, atol=1e-12 * np.max(np.abs(a1)))
    assert np.max(np.abs(a1 - a)) < 1e-3 * np.max(np.abs(a))
    v, dd = SP.integrate(a1, dt)
    assert abs(v[-1]) < 1e-12 and abs(dd[-1]) < 1e-12


def test_target_spectrum_matches_rg160_extrapolations():
    """TargetSpectrum built from RG 1.60 points reproduces the constant-SD branch below fD and the ZPA."""
    f = np.array([0.25, 2.5, 9.0, 33.0])
    t = EL.TargetSpectrum(f, EL.rg160_spectrum(f, "H"), 0.05)
    ff = np.array([0.05, 0.1, 0.17, 0.6, 1.3, 5.0, 20.0, 40.0, 90.0])
    np.testing.assert_allclose(t(ff), EL.rg160_spectrum(ff, "H"), rtol=1e-12)


def test_sa_ratio_between_tripartite_points():
    """Inside B-C the RG 1.60 line is log-log straight: SA at the geometric mean is the geometric mean."""
    fm = np.sqrt(2.5 * 9.0)
    assert EL.rg160_spectrum(np.array([fm]), "H")[0] == pytest.approx(np.sqrt(3.13 * 2.61), rel=1e-12)


def test_generated_motion_has_no_long_period_drift():
    f = EL.check_grid(0.3, 50.0, 20)
    tgt = EL.TargetSpectrum(f, EL.rg160_spectrum(f, "H", 0.05, 0.3), 0.05)
    dt, dur = 0.01, 20.0
    for seed in (7, 8):
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed, spawn_key=(1,))))
        res = EL.match_spectrum(tgt, dt, dur, 0.05, rng=rng,
                                options=EL.MatchOptions(lw_iterations=8, wavelet_iterations=6))
        _, d = SP.integrate(res.acc * 980.665, dt)
        D = np.abs(np.fft.rfft(d, 1 << 14)) ** 2
        fr = np.fft.rfftfreq(1 << 14, dt)
        frac = D[fr < 1.0 / dur].sum() / D.sum()      # periods longer than the record (target starts at 0.3 Hz)
        assert frac < 0.5, f"seed {seed}: {frac:.2f} of the displacement energy below {1 / dur:g} Hz"
        # RG 1.60 is anchored to PGD = 36 in per g; allow a factor 1.5
        assert np.max(np.abs(d)) < 1.5 * 0.3 * 36.0 * 2.54


# =======================================================================================
# EQUAKE module
# =======================================================================================
FREQS = [0.2, 0.25, 0.5, 1.0, 2.5, 5.0, 9.0, 20.0, 33.0, 50.0]


def _equake_deck(rows, **kw):
    d = decks.new("EQUAKE")
    d.params.update(dict(title="review", model="e", accopt=0, nrfreq=len(FREQS), rand=11975, damp=0.05, dur=20.0,
                         corr=0, seeds=1, tpsd=0, delt=0.01, gravity=32.2, eqtit="review"))
    d.params.update(kw)
    for r in rows:
        d.table("spectra").append(r)
    return d


def test_external_history_british_units_psd_fft(tmp_path):
    """accopt 2 with British gravity: .psd in in^2/s^3, .fft in in/s of the 5-75 % Arias window."""
    rng = np.random.default_rng(5)
    dt = 0.01
    t = np.arange(2001) * dt
    a = 0.15 * rng.standard_normal(2001) * EL.saragoni_hart_envelope(t, 20.0)
    textfiles.write_history(tmp_path / "rec.acc", a, dt)
    f = np.asarray(FREQS)
    EL.write_spectrum_file(tmp_path / "t.rsi", f, EL.rg160_spectrum(f, "H", 0.05, 0.3))
    decks.write(tmp_path / "e.equ", _equake_deck([[1, "t.rsi", "rec.rso", "rec.acc", "", ""]], accopt=2))
    assert run_module("EQUAKE", "e", tmp_path) == 0
    a_back, _ = textfiles.read_history(tmp_path / "rec.acc")
    g_in = 32.2 * 12.0
    fp, psd = SP.band_averaged_psd(a_back * g_in, dt)
    got = textfiles.read_xy(tmp_path / "rec.psd")
    np.testing.assert_allclose(got[:, 0], fp, atol=1e-6)
    np.testing.assert_allclose(got[:, 1], psd, rtol=1e-7)
    t5, t75 = SP.strong_motion_window(a_back, dt)
    seg = a_back[int(round(t5 / dt)):int(round(t75 / dt)) + 1] * g_in
    F = dt * np.fft.rfft(seg)
    fft = textfiles.read_xy(tmp_path / "rec.fft")
    np.testing.assert_allclose(fft[:, 1] + 1j * fft[:, 2], F, rtol=1e-6, atol=1e-7 * np.max(np.abs(F)))
    # no .vel/.dis for an external history
    assert not (tmp_path / "rec.vel").exists()


@pytest.mark.parametrize("kw,msg", [(dict(rand=0), "Error 90"), (dict(dur=0.0), "Error 92"),
                                    (dict(nrfreq=0), "Error 91")])
def test_equake_input_errors(tmp_path, kw, msg):
    f = np.asarray(FREQS)
    EL.write_spectrum_file(tmp_path / "t.rsi", f, EL.rg160_spectrum(f, "H", 0.05, 0.3))
    decks.write(tmp_path / "e.equ", _equake_deck([[1, "t.rsi", "o.rso", "", "o.acc", ""]], **kw))
    assert run_module("EQUAKE", "e", tmp_path) == 1
    assert msg in (tmp_path / "e_EQUAKE.out").read_text()


def test_failed_matching_criteria_are_warned(tmp_path):
    rng = np.random.default_rng(9)
    dt = 0.005
    t = np.arange(4001) * dt
    a = 0.05 * rng.standard_normal(4001) * EL.saragoni_hart_envelope(t, 20.0)     # far below the target
    textfiles.write_history(tmp_path / "rec.acc", a, dt)
    f = np.asarray(FREQS)
    EL.write_spectrum_file(tmp_path / "t.rsi", f, EL.rg160_spectrum(f, "H", 0.05, 0.5))
    decks.write(tmp_path / "e.equ", _equake_deck([[1, "t.rsi", "rec.rso", "rec.acc", "", ""]], accopt=2, delt=dt,
                                                 gravity=9.81))
    assert run_module("EQUAKE", "e", tmp_path) == 0
    out = (tmp_path / "e_EQUAKE.out").read_text()
    assert "[FAIL]" in out                                   # the criteria do fail ...
    assert "*** WARNING" in out                              # ... and must be warned (output and screen)


def test_match_spectrum_target_above_wavelet_band():
    tg = EL.TargetSpectrum(np.array([25.0, 50.0, 100.0]), np.array([1.5, 1.2, 1.0]), 0.05)
    res = EL.match_spectrum(tg, 0.005, 20.0, 0.05, rng=np.random.default_rng(1),
                            options=EL.MatchOptions(lw_iterations=2, wavelet_iterations=2))
    assert res.acc.size == 4001
