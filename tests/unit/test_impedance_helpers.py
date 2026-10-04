"""Helpers of the impedance verification problems (sassi/verify/problems/vp_impedance.py): soil
column, extrapolation, peak location, the independent half-space BEM references and the R2 F.2
3-DOF model."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.verify.problems import vp_impedance as V


# ---------------------------------------------------------------------------------------
# column, frequencies, extrapolation, peaks
# ---------------------------------------------------------------------------------------
def test_graded_thicknesses():
    t = V.graded_thicknesses(0.5, 4.0, 60.0)
    assert t[0] == 0.5
    assert sum(t) >= 60.0 and sum(t[:-1]) < 60.0
    assert max(t) == 4.0
    assert all(b >= a for a, b in zip(t, t[1:]))
    assert t[1] == pytest.approx(0.5 * V.GROWTH)
    with pytest.raises(ValueError):
        V.graded_thicknesses(0.0, 1.0, 1.0)


def test_soil_column_is_uniform_halfspace():
    s = V.soil_column(1.0, 4.0, 30.0)
    assert s.nl == 20 and s.hslaw is None                        # deck default law (UNIFORM, D-W1-01)
    for lay in s.layers + [s.halfspace]:
        assert lay.vs == V.VS and lay.vp == pytest.approx(V.VPP) and lay.rho == V.RHO
        assert lay.beta_s == V.BETA_ELASTIC


def test_a0_frequencies():
    fs, a0 = V.a0_frequencies([0.05, 0.5, 1.5], 5.0)
    assert list(fs.fnum) == [5, 50, 150]
    np.testing.assert_allclose(2 * np.pi * fs.freq * 5.0 / V.VS, a0, rtol=1e-12)
    assert fs.nft // 2 >= max(fs.fnum)


def test_richardson_and_observed_order():
    h = [1.0, 0.5, 0.25]
    v = [3.0 + 0.7 * x ** 1.5 for x in h]
    assert V.observed_order(h, v) == pytest.approx(1.5, rel=1e-8)
    lin = [2.0 + 0.3 * x for x in h]
    assert V.richardson(h, lin) == pytest.approx(2.0, rel=1e-12)
    assert V.richardson(h, v, p=1.5) == pytest.approx(3.0, rel=1e-12)
    assert math.isnan(V.observed_order(h, [1.0, 2.0, 1.5]))            # not monotone


def test_peak_near_parabola():
    x = np.arange(0.0, 2.0, 0.01)
    y = 5.0 - (x - 1.2345) ** 2
    assert V.peak_near(x, y, 1.2) == pytest.approx(1.2345, abs=1e-10)
    assert math.isnan(V.peak_near(x, y, 0.5, window=0.1))              # no maximum in the window


def test_relaxed_impedance_blocks():
    """Relaxed contact = welded contact when the flexibility has no z-x coupling."""
    rng = np.random.default_rng(1)
    xy = rng.uniform(-1, 1, (6, 2))
    T = V.rigid_transform(xy)
    A = rng.normal(size=(18, 18))
    F = A @ A.T + 18 * np.eye(18)
    iz = np.arange(2, 18, 3)
    ih = np.setdiff1d(np.arange(18), iz)
    F[np.ix_(iz, ih)] = 0.0
    F[np.ix_(ih, iz)] = 0.0
    K = V.relaxed_impedance(F, T)
    Kw = T.T @ np.linalg.solve(F, T)
    np.testing.assert_allclose(K, Kw, atol=1e-10 * np.abs(Kw).max())


# ---------------------------------------------------------------------------------------
# independent references
# ---------------------------------------------------------------------------------------
def test_static_bem_reproduces_exact_disk_values():
    """Exact-kernel BEM, n = 8 and 16 rings extrapolated: relaxed punch 4GR/(1-nu), 8GR^3/(3(1-nu)),
    Mossakovskii's welded vertical 4GR ln(3-4nu)/(1-2nu) and Reissner-Sagoci torsion 16GR^3/3."""
    nu = 1.0 / 3.0
    Kr = V.bem_reference("disk", (8, 16), nu, relaxed=True)
    Kw = V.bem_reference("disk", (8, 16), nu)
    assert Kr[2, 2] == pytest.approx(4 / (1 - nu), rel=3e-3)
    assert Kr[3, 3] == pytest.approx(8 / (3 * (1 - nu)), rel=3e-3)
    assert Kw[2, 2] == pytest.approx(4 * math.log(3 - 4 * nu) / (1 - 2 * nu), rel=3e-3)
    assert Kw[5, 5] == pytest.approx(16 / 3, rel=3e-3)
    assert Kw[0, 0] == pytest.approx(Kw[1, 1], rel=1e-6)                # isotropic
    assert Kw[0, 4] < 0                                                   # welded sliding-rocking coupling


def test_lamb_kernels_reproduce_wong_table():
    """Exact dynamic surface Green's functions vs Wong (1975) / Apsel (1979), nu = 0.33 (R2 D.3)."""
    nu = 0.33
    r0 = np.array([0.5, 1.0, 2.0])
    Dz, DA, DB = V.lamb_kernels(1.0, r0, nu)
    Uz0 = (1 - nu) / (2 * np.pi) + r0 * Dz
    Ur1 = r0 * ((2 - nu) / (4 * np.pi * r0) + DA + nu / (4 * np.pi * r0) + DB)
    Ut1 = -r0 * ((2 - nu) / (4 * np.pi * r0) + DA - nu / (4 * np.pi * r0) - DB)
    wong = {"Uz0": [(.087, -.062), (.037, -.102), (-.087, -.077)], "Ur1": [(.146, -.058), (.112, -.105), (.011, -.137)],
            "Ut1": [(-.089, .058), (-.045, .099), (.075, .093)]}
    for name, got in (("Uz0", Uz0), ("Ur1", Ur1), ("Ut1", Ut1)):
        ref = np.array([complex(*v) for v in wong[name]])
        assert np.abs(got - ref).max() < 0.003, (name, got, ref)


def test_dynamic_bem_static_limit_and_damping_sign():
    cells = V.disk_cells(4)
    K0 = V.static_bem(cells, relaxed=True)
    Kd = V.relaxed_bem_impedance(cells, 1e-4)
    np.testing.assert_allclose(Kd.real, K0, rtol=1e-4, atol=1e-8)
    K1 = V.relaxed_bem_impedance(cells, 1.0)
    assert all(K1[i, i].imag > 0 for i in range(6))                     # radiation damping (e^{+iwt})
    assert K1[3, 3].real < K0[3, 3]                                       # rocking softens with a0


def test_lowfreq_damping_limits():
    """Exact low-frequency damping of the disk; consistent with the stored BEM table at a0 = 0.25."""
    d = V.lowfreq_damping(1.0 / 3.0)
    assert d["c_h0"] == pytest.approx(0.580, abs=0.003)
    assert d["c_v0"] == pytest.approx(0.782, abs=0.003)
    assert d["A_r"] == pytest.approx(0.240, abs=0.003)
    # the BEM at a0 = 0.25 is within the O(a0^2) correction of the low-frequency limit (Veletsos-Verbic: 0.1 a0^2)
    assert 0.9 < V.VP11_BEM[0.25]["cr"] / (d["A_r"] * 0.25 ** 2) < 1.0
    assert V.VP11_BEM[0.25]["ch"] == pytest.approx(d["c_h0"], rel=0.01)


# ---------------------------------------------------------------------------------------
# R2 F.2
# ---------------------------------------------------------------------------------------
def test_three_dof_static_springs_give_veletsos_meek_period():
    """Massless foundation, real constant springs, no damping: the resonance of the 3-DOF model is
    at T~ with (T~/T)^2 = 1 + k/K_x + k h^2/K_theta (R2 F.1 worked example: 1.15805)."""
    m, h, T = 942.48, 10.0, 0.5
    k = m * (2 * np.pi / T) ** 2
    S = np.diag([9.6e5, 8.0e7]).astype(complex)
    ratio = math.sqrt(1 + k / 9.6e5 + k * h * h / 8.0e7)
    assert ratio == pytest.approx(1.15805, abs=5e-6)
    w = 2 * np.pi / (T * ratio)
    det = np.linalg.det(np.array([[k - w * w * m, -w * w * m, -w * w * m * h],
                                  [-w * w * m, S[0, 0] - w * w * m, -w * w * m * h],
                                  [-w * w * m * h, -w * w * m * h, S[1, 1] - w * w * m * h * h]]))
    assert abs(det) < 1e-9 * k * S[0, 0].real * S[1, 1].real
    # the total motion of the mass is the fixed-base TF when the springs are rigid
    Hr = V._total_motion(0.9 * 2 * np.pi / T, m, h, k * (1 + 0.1j), np.diag([1e15, 1e18]).astype(complex))
    r_ = 0.9
    assert Hr == pytest.approx((1 + 0.1j) / (1 + 0.1j - r_ ** 2), rel=1e-6)


def test_three_dof_reproduces_r2_worked_example():
    """R2 F.1: F.2 with Veletsos-Verbic impedances gives the resonance 1.717 Hz and peak 6.71."""
    m, h, T, R = 942.48, 10.0, 0.5, 10.0
    k = m * (2 * np.pi / T) ** 2
    f = np.arange(1.6, 1.85, 0.0005)
    H = np.array([abs(V._total_motion(2 * np.pi * x, m, h, k * (1 + 0.1j), V._vv_impedance(2 * np.pi * x * R / V.VS, R)))
                  for x in f])
    assert f[np.argmax(H)] == pytest.approx(1.717, abs=1e-3)
    assert H.max() == pytest.approx(6.71, abs=5e-3)


# ---------------------------------------------------------------------------------------
# review fixes: column cap, kernel table range, welded dynamic BEM, layer dispersion,
# independent K_G check, work-file cleanup, context-local FILE11 switch
# ---------------------------------------------------------------------------------------
def test_graded_thicknesses_caps_the_first_sublayer():
    """t0 > tmax: every sublayer (the first included) obeys the h <= lambda/8 cap."""
    t = V.graded_thicknesses(2.0, 1.0, 10.0)
    assert t[0] == 1.0 and max(t) == 1.0 and sum(t) == pytest.approx(10.0)
    assert V.graded_thicknesses(0.5, 4.0, 60.0)[0] == 0.5                 # t0 <= tmax unchanged


def test_kernel_table_covers_every_pair_distance_of_a_square():
    """The dynamic-kernel table of the BEM reaches the largest centroid-to-Gauss-point distance of a
    square (2.43 B for 4 x 4 cells, beyond the former 2.2 max|V| = 2.2 B) and reproduces the
    wavenumber integrals there."""
    cells = V.square_cells(4)
    _, Pc, _ = V.bem_influence(cells, 1.0 / 3.0)
    Q, _ = V._cell_quadrature(cells, 4)
    d = np.hypot(*(Pc[:, None, None, :] - Q[None]).transpose(3, 0, 1, 2))
    rmax = V.max_pair_distance(Pc, Q)
    assert d.max() > 2.2 and rmax >= d.max()
    r = np.array([d.max(), 0.5 * d.max(), 0.123])
    spl = V.kernel_table(1.0, rmax, 1.0 / 3.0, coupling=True)
    exact = V.lamb_kernels(1.0, r, 1.0 / 3.0, coupling=True)
    for s, e in zip(spl, exact):
        np.testing.assert_allclose(s(r), e, atol=1e-7)


def test_coupling_kernel_reproduces_wong_ur0_column():
    """G_xz (radial displacement under a vertical force) vs R2 D.3 column U_r0 (Wong 1975, nu = 0.33;
    mu r u_r / P with the force pressing down, hence the minus sign for an upward unit P_z)."""
    nu = 0.33
    r = np.arange(0.5, 5.51, 0.5)
    _, _, _, Dc = V.lamb_kernels(1.0, r, nu, coupling=True)
    U = -r * ((1 - 2 * nu) / (4 * np.pi * r) + Dc)
    wong = np.array([-.032 + .007j, -.033 + .025j, -.020 + .047j, .006 + .060j, .041 + .058j, .074 + .035j,
                     .092 - .005j, .087 - .054j, .056 - .096j, .005 - .120j, -.054 - .116j])
    assert np.abs(U - wong).max() < 0.002
    np.testing.assert_allclose(V.lamb_kernels(1.0, r, nu)[0], V.lamb_kernels(1.0, r, nu, coupling=True)[0])


def test_welded_dynamic_bem_static_limit_symmetry_and_coupling():
    cells = V.disk_cells(4)
    Ks = V.static_bem(cells)                                                # welded statics
    K0 = V.welded_bem_impedance(cells, 1e-4)
    np.testing.assert_allclose(K0.real, Ks, rtol=1e-4, atol=1e-6)
    K1 = V.welded_bem_impedance(cells, 1.07)
    assert np.abs(K1 - K1.T).max() < 2e-3 * np.abs(K1).max()               # collocation: nearly symmetric
    assert K1[0, 4].real < 0 < K1[0, 4].imag                                # dynamic sliding-rocking coupling
    assert K1[1, 3] == pytest.approx(-K1[0, 4], rel=1e-4)                  # polar mesh: not exactly 4-fold
    Kr = V.relaxed_bem_impedance(cells, 1.07)
    assert Kr[0, 4] == 0 and K1[3, 3].real > Kr[3, 3].real                  # welded is stiffer


def test_welded_disk_impedance_at_vp17_resonance():
    """Mesh-extrapolated welded disk impedance (VP-17 guard reference) at a0 = 1.07: the coupling
    K_x-theta is -0.11 Re K_x R (statics: -0.08) and 8/12 agrees with 12/16 within 0.3 %."""
    K = V.welded_disk_impedance([1.07])[0]
    K2 = V.welded_disk_impedance([1.07], meshes=(12, 16))[0]
    for i, j in ((0, 0), (0, 4), (3, 3)):
        assert abs(K[i, j] - K2[i, j]) < 3e-3 * abs(K2[i, j])
    assert -0.13 < K[0, 4].real / K[0, 0].real < -0.09


def test_layer_dispersion_cutoffs_and_zgv_onset():
    """Continuum P-SV dispersion of the VP-13 layer (H = 3, Vs = 1, Vp = 2): at k = 0 the secular
    function vanishes at the shear and compression column frequencies (2n-1) pi V/(2H); the 2nd mode
    is born at a ZGV point A0 = 0.509 (k R = 0.17) below the cut-off Vp/(4H) (A0 = 0.524)."""
    for w in (np.pi / 6, np.pi / 3, np.pi / 2):                           # shear 1, compression 1, shear 2
        assert abs(V.layer_dispersion_det(np.array([0.0]), w, 3.0, 1.0, 2.0)[0]) < 1e-12
    assert abs(V.layer_dispersion_det(np.array([0.0]), 0.8, 3.0, 1.0, 2.0)[0]) > 0.1
    a0z, kr = V.layer_zgv_onset()
    assert 0.505 < a0z < 0.513 < np.pi / 6
    assert 0.12 < kr < 0.22


def test_independent_global_impedance_and_kg_checks():
    rng = np.random.default_rng(5)
    xy = rng.uniform(-1, 1, (7, 2))
    T = V.rigid_transform(xy)
    A = rng.normal(size=(21, 21)) + 1j * rng.normal(size=(21, 21))
    F = A @ A.T + 10 * np.eye(21)                                             # complex symmetric
    KD = V.independent_global_impedance(F, T)
    KG = T.T @ np.linalg.inv(F) @ T
    r = V.VPResult()
    V.kg_checks(r, "random", KG, KD)
    assert r.passed and len(r.checks) == 2
    bad = KG.copy()
    bad[0, 4] *= 1.001
    r2 = V.VPResult()
    V.kg_checks(r2, "perturbed", bad, KD)
    assert [c.passed for c in r2.checks] == [False, False]


def test_purge_work_files_keeps_decks_and_listings(tmp_path, monkeypatch):
    names = ["m.sit", "m.poi", "m.hou", "m.anl", "m.frc", "m_ANALYS.out", "FILE1", "FILE2", "FILE3", "FILE8",
             "FILE11", "COOSK", "FOUNSTIF", "m.N4", "FILE9001"]
    for n in names:
        (tmp_path / n).write_bytes(b"x" * 10)
    monkeypatch.setattr(V, "KEEP_WORK_FILES", True)
    assert V.purge_work_files(tmp_path) == 0 and len(list(tmp_path.iterdir())) == len(names)
    monkeypatch.setattr(V, "KEEP_WORK_FILES", False)
    assert V.purge_work_files(tmp_path) == 90
    assert sorted(p.name for p in tmp_path.iterdir()) == ["m.anl", "m.frc", "m.hou", "m.poi", "m.sit", "m_ANALYS.out"]


def test_file11_switch_is_local_to_the_calling_thread(tmp_path):
    """file11_without_matrices no longer patches the analys module global: an ANALYS run in another
    thread during the ``with`` block still stores X_ff in FILE11 (D-ANL-07)."""
    import threading

    from sassi.io.container import read_container
    from sassi.modules import analys as A
    from sassi.verify import builders as B

    R = 5.0
    site = V.soil_column(R / 2, 4 * R, 10 * R)
    fs, _ = V.a0_frequencies([0.5], R)
    mdl = B.surface_rigid_mat(site, half_width=R, ndiv=2, shape="disk")
    B.run_soil(tmp_path, "m", site, fs, layer=0, rad=mdl.rad)
    B.run_house(tmp_path, "m", mdl)

    def stored() -> bool:
        return bool(read_container(tmp_path / "FILE11", "FILE11").meta["x_stored"])
    seen = {}
    with V.file11_without_matrices():
        worker = threading.Thread(target=lambda: (seen.setdefault("limit", A.file11_x_limit()),
                                                  B.run_analys(tmp_path, "m", fs, impe=2)))
        worker.start()
        worker.join()
        seen["worker_stored"] = stored()
        seen["inside"] = A.file11_x_limit()
        B.run_analys(tmp_path, "m", fs, impe=2)
        seen["caller_stored"] = stored()
    assert seen == {"limit": A.FILE11_MAX_X_BYTES, "worker_stored": True, "inside": 0, "caller_stored": False}
    assert A.file11_x_limit() == A.FILE11_MAX_X_BYTES
    V.purge_work_files(tmp_path)
