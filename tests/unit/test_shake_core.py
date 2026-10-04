"""Unit tests of sassi.core.shake (SHAKE recursion, strains, curves, FILE73/FILE88; R1 section 6)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.core import shake as SH


def _column(thick=(10.0, 15.0, 20.0, 0.0), vs=(150.0, 250.0, 400.0, 800.0), rho=(1.8, 1.9, 2.0, 2.2),
            beta=(0.05, 0.04, 0.03, 0.01), labels=None, g=9.81):
    w = np.asarray(rho) * g
    return SH.SoilColumn(np.asarray(thick), w, np.asarray(vs), np.asarray(beta), g,
                         list(labels) if labels else [])


# ------------------------------------------------------------------ curves
def test_interp_log_strain_rules():
    xs, ys = [1e-4, 1e-3, 1e-2], [1.0, 0.8, 0.4]
    assert SH.interp_log_strain(1e-5, xs, ys) == 1.0              # constant below
    assert SH.interp_log_strain(1.0, xs, ys) == 0.4               # constant above
    assert SH.interp_log_strain(10 ** -3.5, xs, ys) == pytest.approx(0.9, abs=1e-14)   # linear in log10
    # repeated abscissa (padded curve) and unsorted input are accepted
    assert SH.interp_log_strain(3e-3, [1e-2, 1e-4, 1e-3, 1e-2], [0.4, 1.0, 0.8, 0.4]) == pytest.approx(
        0.8 - 0.4 * np.log10(3.0), abs=1e-12)


def test_curves_from_dynp_rows_keep_order_and_units():
    rows = [dict(label="Sand", kind="G", strain=1e-4, value=1.0), dict(label="Sand", kind="D", strain=1e-4, value=0.5),
            dict(label="Clay", kind="g", strain=1e-4, value=1.0), dict(label="Clay", kind="d", strain=1.0, value=20.0),
            dict(label="Sand", kind="G", strain=1.0, value=0.1), dict(label="Sand", kind="D", strain=1.0, value=20.0)]
    cv = SH.curves_from_dynp_rows(rows)
    assert list(cv) == ["Sand", "Clay"]
    assert cv["Sand"].damping(1e-4) == pytest.approx(0.005)       # damping percent -> fraction
    assert cv["Sand"].g_over_gmax(1.0) == pytest.approx(0.1)
    with pytest.raises(ValueError):
        SH.curves_from_dynp_rows([dict(label="X", kind="Q", strain=1.0, value=1.0)])


# ------------------------------------------------------------------ recursion
@pytest.mark.parametrize("form", [0, 1])
def test_uniform_layer_closed_forms(form):
    """R2 A.1: surface/base-within = 1/cos(k*H); surface/rock-outcrop = 1/(cos k*H + i a* sin k*H)."""
    H, vs, rho, xi, vr, rr, xr = 30.0, 200.0, 2.0, 0.05, 1000.0, 2.2, 0.01
    col = _column(thick=(H, 0.0), vs=(vs, vr), rho=(rho, rr), beta=(xi, xr))
    f = np.linspace(0.05, 12.0, 300)
    k = 2 * np.pi * f / (vs * np.sqrt(cfactor(xi, form)))
    a = rho * vs * np.sqrt(cfactor(xi, form)) / (rr * vr * np.sqrt(cfactor(xr, form)))
    hw = SH.column_transfer(f, col, col.gmax, col.beta0, 1, True, 2, False, form)
    ho = SH.column_transfer(f, col, col.gmax, col.beta0, 1, True, 2, True, form)
    np.testing.assert_allclose(hw, 1 / np.cos(k * H), rtol=1e-12)
    np.testing.assert_allclose(ho, 1 / (np.cos(k * H) + 1j * a * np.sin(k * H)), rtol=1e-12)


def test_subdivision_is_exact():
    """Splitting a uniform sublayer does not change any transfer function (R2 A.2)."""
    c1 = _column(thick=(10.0, 30.0, 0.0), vs=(150.0, 300.0, 900.0), rho=(1.8, 2.0, 2.2), beta=(0.05, 0.03, 0.01))
    c2 = _column(thick=(10.0, 12.0, 18.0, 0.0), vs=(150.0, 300.0, 300.0, 900.0), rho=(1.8, 2.0, 2.0, 2.2),
                 beta=(0.05, 0.03, 0.03, 0.01))
    f = np.linspace(0.0, 25.0, 200)
    h1 = SH.column_transfer(f, c1, c1.gmax, c1.beta0, 1, True, 3, True)
    h2 = SH.column_transfer(f, c2, c2.gmax, c2.beta0, 1, True, 4, True)
    np.testing.assert_allclose(h1, h2, rtol=1e-11)
    # within motion at the top of the split layer (interface 2) is identical too
    np.testing.assert_allclose(SH.column_transfer(f, c1, c1.gmax, c1.beta0, 2, False, 3, True),
                               SH.column_transfer(f, c2, c2.gmax, c2.beta0, 2, False, 4, True), rtol=1e-11)


def test_free_surface_and_static_limit():
    col = _column()
    f = np.array([0.0, 1.0, 3.0])
    vstar = SH.complex_velocity(col.gmax, col.rho, col.beta0)
    E, F = SH.wave_amplitudes(f, col.thick, col.rho, vstar)
    np.testing.assert_allclose(E[:, 0], F[:, 0])                      # E1 = F1 (free surface)
    np.testing.assert_allclose(E[0], 1.0)                              # f = 0: rigid-body motion
    np.testing.assert_allclose(F[0], 1.0)


def test_mid_layer_strain_matches_finite_difference():
    """gamma = du/dz at mid-height, with u from the E/F waves (acceleration amplitudes / -w^2)."""
    col = _column()
    f = np.array([0.7, 2.3, 6.1])
    w = 2 * np.pi * f
    vstar = SH.complex_velocity(col.gmax, col.rho, col.beta0)
    E, F = SH.wave_amplitudes(f, col.thick, col.rho, vstar)
    g = SH.mid_layer_strain_factor(f, E, F, col.thick, vstar)
    for m in range(col.n - 1):
        k = w / vstar[m]
        z0, eps = 0.5 * col.thick[m], 1e-4

        def u(z):
            return -(E[:, m] * np.exp(1j * k * z) + F[:, m] * np.exp(-1j * k * z)) / w ** 2
        fd = (u(z0 + eps) - u(z0 - eps)) / (2 * eps)
        np.testing.assert_allclose(g[:, m], fd, rtol=1e-6)


@pytest.mark.parametrize("outcrop", [True, False])
@pytest.mark.parametrize("cl", [1, 3, 4])
def test_control_motion_is_reproduced(outcrop, cl):
    col = _column()
    rng = np.random.default_rng(5)
    dt, nfft = 0.01, 1024
    acc = rng.standard_normal(700) * np.hanning(700)
    res = SH.run_shake(acc, dt, nfft, col, {}, cl, outcrop, 0.65, 0)
    np.testing.assert_allclose(res.acc_history(cl, outcrop)[:700], acc, atol=1e-12 * np.max(np.abs(acc)))
    assert res.iterations == []


def test_cutoff_removes_components():
    col = _column()
    rng = np.random.default_rng(1)
    acc = rng.standard_normal(512)
    res = SH.run_shake(acc, 0.01, 512, col, {}, 4, True, 0.65, 0, cutoff=10.0)
    S = np.fft.rfft(res.acc_history(1, True))
    f = np.fft.rfftfreq(512, 0.01)
    assert np.max(np.abs(S[f > 10.0 + 1e-9])) < 1e-10 * np.max(np.abs(S))


def test_iterations_exact_count_and_properties():
    sand = SH.DynamicProperty("Sand", np.array([1e-4, 1e-3, 1e-2, 1e-1, 1.0]), np.array([1.0, 0.97, 0.8, 0.4, 0.1]),
                              np.array([1e-4, 1e-3, 1e-2, 1e-1, 1.0]), np.array([0.5, 1.0, 3.0, 9.0, 20.0]))
    col = _column(labels=["Sand", "Sand", "", ""])
    rng = np.random.default_rng(2)
    acc = 3.0 * rng.standard_normal(800) * np.hanning(800)          # m/s^2
    res = SH.run_shake(acc, 0.01, 1024, col, {"Sand": sand}, 4, True, 0.65, 5)
    assert [it.number for it in res.iterations] == [1, 2, 3, 4, 5]
    last = res.iterations[-1]
    # final properties are the last update; strain-compatible on the iterated sublayers
    np.testing.assert_allclose(res.G[:2], last.G_new[:2])
    for i in range(2):
        assert res.G[i] == pytest.approx(col.gmax[i] * sand.g_over_gmax(last.gamma_eff[i]))
        assert res.beta[i] == pytest.approx(sand.damping(last.gamma_eff[i]))
    # the unlabelled sublayer and the half-space keep their low-strain values
    assert res.G[2] == col.gmax[2] and res.G[3] == col.gmax[3]
    assert res.beta[2] == col.beta0[2]
    # iteration 1 uses the low-strain properties
    np.testing.assert_allclose(res.iterations[0].G_used, col.gmax[:-1])
    # vertical input: no update
    res_v = SH.run_shake(acc, 0.01, 1024, col, {"Sand": sand}, 4, True, 0.65, 5, iterate=False)
    assert res_v.iterations == [] and np.all(res_v.G == col.gmax)
    np.testing.assert_array_equal(res_v.gamma_eff_compatible, res_v.gamma_eff)


def test_strain_compatible_pairs_after_one_iteration():
    """SHAKE91 Table B-2 pairing (FILE88): the reported effective strain is the one the final G and
    beta were read from (last iteration), not that of the extra response with the final properties."""
    sand = SH.DynamicProperty("Sand", np.array([1e-4, 1e-3, 1e-2, 1e-1, 1.0]), np.array([1.0, 0.97, 0.8, 0.4, 0.1]),
                              np.array([1e-4, 1e-3, 1e-2, 1e-1, 1.0]), np.array([0.5, 1.0, 3.0, 9.0, 20.0]))
    col = _column(labels=["Sand", "Sand", "", ""])
    rng = np.random.default_rng(2)
    acc = 4.0 * rng.standard_normal(800) * np.hanning(800)
    res = SH.run_shake(acc, 0.01, 1024, col, {"Sand": sand}, 4, True, 0.65, 1)
    g = res.gamma_eff_compatible
    np.testing.assert_array_equal(g, res.iterations[-1].gamma_eff)
    for i in range(2):
        assert res.G[i] == pytest.approx(col.gmax[i] * sand.g_over_gmax(g[i]), rel=1e-12)
        assert res.beta[i] == pytest.approx(sand.damping(g[i]), rel=1e-12)
    # not converged after one update: the final response's strain differs (why the pairing matters)
    assert np.max(np.abs(res.gamma_eff[:2] / g[:2] - 1)) > 0.01


def test_max_amplification_rigid_like_base():
    col = _column(thick=(30.0, 0.0), vs=(200.0, 1e5), rho=(2.0, 2.0), beta=(0.05, 0.0))
    amp, fpk = SH.max_amplification(col, col.gmax, col.beta0, 1, True, 2, False, 5.0)
    assert fpk == pytest.approx(200.0 / 120.0, rel=1e-2)          # near f1 = Vs/(4H)
    assert amp == pytest.approx(1 / abs(np.cos(2 * np.pi * fpk * 30 / (200 * np.sqrt(cfactor(0.05))))), rel=1e-9)


# ------------------------------------------------------------------ files
def test_file88_round_trip(tmp_path):
    n = 3
    arr = [np.arange(1, n + 1) * k for k in (1.0, 1e-3, 1e3, 2e2, 1e-2, 4e2, 1e-2)]
    SH.write_file88(tmp_path / "FILE88", *arr, units="m kN/m2 m/s", halfspace=[0, 22, 2000, 1000, 0.01, 0.01])
    d = SH.read_file88(tmp_path / "FILE88")
    assert list(d["layer"]) == [1, 2, 3]
    for c, a in zip(SH.FILE88_COLUMNS[1:], arr):
        np.testing.assert_allclose(d[c], a, rtol=1e-8)


def test_file73_round_trip(tmp_path):
    cv = {"Sand": SH.DynamicProperty("Sand", np.array([1e-4, 1.0]), np.array([1.0, 0.1]), np.array([1e-4, 1e-2, 1.0]),
                                     np.array([0.5, 3.0, 20.0])),
          "Clay A": SH.DynamicProperty("Clay A", np.array([1e-4]), np.array([1.0]), np.array([1e-4]), np.array([1.0]))}
    SH.write_file73(tmp_path / "FILE73", cv)
    back = SH.read_file73(tmp_path / "FILE73")
    assert list(back) == [1, 2]
    assert back[1].label == "Sand" and back[2].label == "Clay A"
    np.testing.assert_allclose(back[1].d_pct, [0.5, 3.0, 20.0])
