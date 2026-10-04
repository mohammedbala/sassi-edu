"""Adversarial review tests of the ANALYS work package (sassi/modules/analys.py, sassi/core/ssi_solver.py,
sassi/verify/builders.py; requirements 2.5, 4.6, D-ANL-01..13, R1 4.1-4.6, spec 05c A.5).

The tests probe ANALYS from angles the package's own tests do not take:

* independent dense re-derivations of the flexible-volume equation with *fixed* interaction
  translations (the gather ``A_f`` drops rows/columns, the seismic load keeps the full ``X U'``);
* reciprocity (Maxwell-Betti) of the complex-symmetric SSI system: vibration TFs are symmetric
  between structural and interaction DOFs, and a seismic TF equals the vibration influence line
  times the free-field load ``X_ff U'_f`` (inclined wave, rotated axes, control point off origin);
* an inclined single wave: moving the control point (xc, yc) multiplies every TF by the phase
  ``exp(+i k dx')`` (and only the component along x' = R(ang) x matters);
* rotated SV input equals the superposition cos(a) H_SV + sin(a) H_SH of independent SITE runs;
* the low-frequency rigid-body limit with a rotated input (G-19);
* global impedance: ``K_G`` (FILE11) is the inverse of the compliance of a stiff massless mat
  obtained from 6 simultaneous vibration load cases at its centre; the FOUN* text files match
  FILE11; ``K_G = T^T X T`` with ``T`` built from the cross product ``u = u0 + theta x d``;
* restart equivalence on the *full sparse LU* fallback path (restored triangular factors);
* restart records must belong to the analysis frequency step (D-ANL-05: records keyed by
  frequency value, every record stores f).

Tests marked ``xfail(strict=True)`` document defects found by the review.  (The restart
frequency-step test was such a test; the defect is fixed -- read_index / _restart_record compare
df and f -- and the marker was removed.)
"""
from __future__ import annotations

import math
import shutil

import numpy as np
import pytest

from sassi.core import ssi_solver as ss
from sassi.core.flexibility import flexibility_matrix, frequency_row
from sassi.core.freefield import free_field_at_nodes
from sassi.io.container import read_container
from sassi.verify import builders as B

FS = B.FrequencySet.fourier(0.01, 1024, [1, 4, 20, 60])          # df = 0.09766 Hz


def _site() -> B.Site:
    return B.layered_site([(2.0, 150.0, 300.0, 1.9, 0.05), (3.0, 250.0, 500.0, 2.0, 0.04)],
                          (400.0, 800.0, 2.1, 0.03))


def _stick_model() -> B.Model:
    stick = B.Stick(heights=[3.0, 6.0], masses=[40.0, 30.0], E=3.0e7, A=1.0, I=0.5, beta=0.05)
    return B.stick_on_mat(_site(), stick, half_width=2.0, ndiv=2, mat_mass=20.0)


def _relmax(a, b) -> float:
    a, b = np.asarray(a), np.asarray(b)
    return float(np.abs(a - b).max()) / max(float(np.abs(b).max()), 1e-300)


def _rc(wd, fs=FS, **params) -> int:
    return B.run_analys(wd, "m", fs, check=False, **params)


def _listing(wd) -> str:
    return B.listing(wd, "m", "ANALYS")


@pytest.fixture(scope="module")
def base_dir(tmp_path_factory):
    """SITE (vertical SV) + POINT + HOUSE of a small stick on a rigid surface mat."""
    wd = tmp_path_factory.mktemp("review_analys_base")
    mdl = _stick_model()
    B.run_soil(wd, "m", _site(), FS, layer=0, rad=mdl.rad)
    B.run_house(wd, "m", mdl)
    return wd, mdl


@pytest.fixture
def wd(base_dir, tmp_path):
    src, mdl = base_dir
    dst = tmp_path / "m"
    shutil.copytree(src, dst)
    return dst, mdl


def _inclined_sv(wd, angle: float = 25.0) -> None:
    """Replace FILE1 by a single inclined SV field (control direction x')."""
    B.write_deck(wd, "m", B.site_deck(_site(), FS, waves=[(2, 1, 1.0, 1.0, angle)], cm=0, model="m"))
    B.run("SITE", wd, "m")


# ======================================================================================
# 1. Impedance and the dense flexible-volume equation with fixed interaction translations
# ======================================================================================
def test_coox_holds_the_symmetric_inverse_of_the_flexibility(wd):
    """D-ANL-04: COOXqqq = X_ff = F_ff^-1 (complex symmetric, node-major x, y, z)."""
    wd, _ = wd
    assert _rc(wd, save=1) == 0
    f3, f4 = read_container(wd / "FILE3"), read_container(wd / "m.N4")
    xyz, iface = np.asarray(f4["x_int_xyz"]), np.asarray(f4["int_iface"])
    idx = read_container(wd / "COOXI")
    for n, o in zip(np.asarray(idx["fnum"]), np.asarray(idx["order"])):
        X = np.asarray(read_container(wd / f"COOX{int(o):03d}")["X"])
        F = flexibility_matrix(f3, frequency_row(f3, int(n)), xyz[:, :2], iface)
        with np.errstate(all="ignore"):                     # spurious Accelerate matmul flags (numpy 2.0)
            XF = X @ F
        np.testing.assert_allclose(XF, np.eye(F.shape[0]), rtol=0, atol=1e-10)
        assert np.abs(X - X.T).max() <= 1e-14 * np.abs(X).max()


def test_fixed_interaction_translations_match_an_independent_dense_solve(tmp_path):
    """Interaction nodes with some translations fixed (D on UZ of one node, UX/UY of another):
    the remaining DOFs satisfy (C_s - C_e + X_aa) U_a = (X U')_a, where X = F^-1 of *all* interaction
    translations (a fixed DOF has zero total motion).  FILE8 holds every FILE4 equation, rotations of
    6-DOF nodes included, and no fixed DOF (D-ANL-13, FILE8 schema)."""
    mdl = _stick_model()
    corners = [n for n in mdl["mat"] if n != mdl["centre"]]
    mdl.house.fix(corners[0], [3])
    mdl.house.fix(corners[1], [1, 2])
    B.run_soil(tmp_path, "m", _site(), FS, layer=0, rad=mdl.rad)
    f4 = B.run_house(tmp_path, "m", mdl)
    assert _rc(tmp_path) == 0
    f8 = B.read_file8(tmp_path)
    np.testing.assert_array_equal(f8["eq_node"], f4["eq_node"])
    np.testing.assert_array_equal(f8["eq_dof"], f4["eq_dof"])
    eqset = set(zip(np.asarray(f8["eq_node"]).tolist(), np.asarray(f8["eq_dof"]).tolist()))
    assert (corners[0], 3) not in eqset and (corners[1], 1) not in eqset and (corners[1], 2) not in eqset
    assert all((mdl["top"], d) in eqset for d in (1, 2, 3, 4, 5))          # beam node rotations are output
    f1, f3 = read_container(tmp_path / "FILE1"), read_container(tmp_path / "FILE3")
    ck, cm = read_container(tmp_path / "COOSK"), read_container(tmp_path / "COOSM")
    xyz, iface = np.asarray(f4["x_int_xyz"]), np.asarray(f4["int_iface"])
    eq = np.asarray(f4["int_eq"]).reshape(-1)
    act = eq >= 0
    assert (~act).sum() == 3
    for q, n in enumerate(FS.fnum):
        w = 2 * math.pi * n * FS.df
        X = np.linalg.inv(flexibility_matrix(f3, frequency_row(f3, n), xyz[:, :2], iface))
        A = ((ck.sparse("Ks") - ck.sparse("Ke")).toarray()
             - w * w * (cm.sparse("Ms") - cm.sparse("Me")).toarray()).astype(complex)
        A[np.ix_(eq[act], eq[act])] += X[np.ix_(act, act)]
        Up = free_field_at_nodes(f1, frequency_row(f1, n), xyz, iface, 0.0, 0.0, 0.0).reshape(-1)
        b = np.zeros(A.shape[0], complex)
        with np.errstate(all="ignore"):                     # spurious Accelerate matmul flags
            b[eq[act]] = (X @ Up)[act]
        U = np.linalg.solve(A, b)
        assert _relmax(f8["H"][q], U) < 1e-9


# ======================================================================================
# 2. Reciprocity of the complex-symmetric SSI system
# ======================================================================================
def test_vibration_transfer_functions_are_reciprocal(wd):
    """Maxwell-Betti: the displacement at DOF a for a unit load at b equals the displacement at b for
    a unit load at a -- between a structural DOF (stick top) and an interaction DOF (mat corner), and
    between a moment and a force (rotation per unit force = displacement per unit moment)."""
    wd, mdl = wd
    top, mid = mdl["top"], mdl["stick"][0]
    corner = next(n for n in mdl["mat"] if n != mdl["centre"])
    pairs = [((top, 1), (corner, 3)), ((top, 5), (mid, 1)), ((corner, 2), (mid, 4))]
    k = 0
    for a, b in pairs:
        for (node, dof) in (a, b):
            k += 1
            B.run_force(wd, "m", [(node, dof, 1.0, 0.0)], FS, copy_to=f"FILE9{k:03d}")
    assert _rc(wd, type=1, simul=k) == 0
    k = 0
    for a, b in pairs:
        Ha = B.read_file8(wd, f"FILE8{k + 1:03d}")
        Hb = B.read_file8(wd, f"FILE8{k + 2:03d}")
        k += 2
        g_ba = B.tf(Ha, *b)                  # response at b, load at a
        g_ab = B.tf(Hb, *a)                  # response at a, load at b
        assert np.abs(g_ba).max() > 0
        assert _relmax(g_ba, g_ab) < 1e-8, (a, b)


def test_seismic_tf_equals_vibration_influence_times_free_field_load(wd):
    """Betti between the seismic and the vibration problem of the same system C:
    U_seis[k] = e_k^T C^-1 A^T X U' = sum_f g_k[f] (X U')_f with g_k = C^-1 e_k (C symmetric).
    Inclined SV (k != 0), rotated axes (ang = 30) and a control point off the origin, so the
    rotation, the oblique phase and the load assembly of ANALYS are all checked against the
    vibration solution and the COOX impedance."""
    wd, mdl = wd
    _inclined_sv(wd)
    ang, xc, yc = 30.0, 0.7, -0.4
    assert _rc(wd, save=1, ang=ang, xc=xc, yc=yc) == 0
    Hs = B.read_file8(wd)
    top, mid = mdl["top"], mdl["stick"][0]
    probes = [(top, 1), (top, 2), (mid, 5), (mdl["centre"], 3)]
    for k, (node, dof) in enumerate(probes, start=1):
        B.run_force(wd, "m", [(node, dof, 1.0, 0.0)], FS, copy_to=f"FILE9{k:03d}")
    assert _rc(wd, type=1, simul=len(probes)) == 0
    f4, f1 = read_container(wd / "m.N4"), read_container(wd / "FILE1")
    xyz, iface = np.asarray(f4["x_int_xyz"]), np.asarray(f4["int_iface"])
    eq = np.asarray(f4["int_eq"]).reshape(-1)
    idx = read_container(wd / "COOXI")
    order = {int(n): int(o) for n, o in zip(idx["fnum"], idx["order"])}
    for k, (node, dof) in enumerate(probes, start=1):
        g = np.asarray(B.read_file8(wd, f"FILE8{k:03d}")["H"])
        ref = B.tf(Hs, node, dof)
        val = np.zeros(len(FS.fnum), complex)
        for q, n in enumerate(FS.fnum):
            X = np.asarray(read_container(wd / f"COOX{order[n]:03d}")["X"])
            Up = free_field_at_nodes(f1, frequency_row(f1, n), xyz, iface, ang, xc, yc).reshape(-1)
            with np.errstate(all="ignore"):
                P = X @ Up
            val[q] = np.sum(g[q, eq] * P)
        assert _relmax(val, ref) < 1e-7, (node, dof)


# ======================================================================================
# 3. Control point, angle and superposition
# ======================================================================================
def test_control_point_shift_is_a_pure_phase_for_one_inclined_wave(wd):
    """FILE1 with one inclined SV wave of wavenumber k: moving (xc, yc) by d changes the reference
    phase only, H(d) = exp(+i k d.x') H(0) with x' = (cos a, sin a); the component of d along y'
    does not matter (A.5.8)."""
    wd, _ = wd
    _inclined_sv(wd)
    k = np.asarray(read_container(wd / "FILE1")["k"])
    assert k.shape[0] == 1 and np.abs(k).min() > 0
    for ang in (0.0, 90.0):
        assert _rc(wd, ang=ang) == 0
        H0 = np.asarray(B.read_file8(wd)["H"]).copy()
        for xc, yc in ((3.0, 0.0), (0.0, 3.0), (-1.5, 2.0)):
            assert _rc(wd, ang=ang, xc=xc, yc=yc) == 0
            H = np.asarray(B.read_file8(wd)["H"])
            a = math.radians(ang)
            d = xc * math.cos(a) + yc * math.sin(a)
            expect = np.exp(1j * k[0] * d)[:, None] * H0
            assert _relmax(H, expect) < 1e-9, (ang, xc, yc)


def test_rotated_sv_equals_superposition_of_independent_sv_and_sh_runs(tmp_path):
    """For vertical incidence the SH (y') field of a horizontally layered site is the SV (x') field
    rotated by 90 deg, so ANALYS with FILE1 = SV and angle a must equal cos a FILE8X + sin a FILE8Y of a
    simultaneous run with FILE1X = SV and FILE1Y = SH (independent SITE runs)."""
    mdl = _stick_model()
    B.run_site_xyz(tmp_path, "m", _site(), FS)
    B.write_deck(tmp_path, "m", B.point_deck(FS, 0, mdl.rad, model="m"))
    B.run("POINT", tmp_path, "m")
    B.run_house(tmp_path, "m", mdl)
    assert _rc(tmp_path, simul=1) == 0
    HX = np.asarray(B.read_file8(tmp_path, "FILE8X")["H"]).copy()
    HY = np.asarray(B.read_file8(tmp_path, "FILE8Y")["H"]).copy()
    B.copy_file(tmp_path, "FILE1X", "FILE1")
    for a in (30.0, 200.0):
        assert _rc(tmp_path, ang=a) == 0
        Ha = B.read_file8(tmp_path)
        assert float(Ha.meta["ang"]) == a
        ca, sa = math.cos(math.radians(a)), math.sin(math.radians(a))
        assert _relmax(Ha["H"], ca * HX + sa * HY) < 1e-9


def test_low_frequency_limit_is_the_rotated_rigid_body_motion(wd):
    """G-19: at f_1 = 0.098 Hz the whole model moves with the free field; with ang = 30 deg the
    translational ATFs tend to (cos a, sin a, 0) and the rotations to 0.  A negative angle is
    normalised (D-ANL-11) and gives the same FILE8."""
    wd, _ = wd
    a = 30.0
    assert _rc(wd, ang=a) == 0
    f8 = B.read_file8(wd)
    assert "EDU-18" not in _listing(wd)
    H1 = np.asarray(f8["H"])[0]
    dof = np.asarray(f8["eq_dof"])
    h0 = {1: math.cos(math.radians(a)), 2: math.sin(math.radians(a)), 3: 0.0}
    for d, v in h0.items():
        assert np.abs(H1[dof == d] - v).max() < 0.02, d
    ref = np.asarray(f8["H"]).copy()
    assert _rc(wd, ang=a - 360.0) == 0
    assert "normalised" in _listing(wd)
    np.testing.assert_array_equal(np.asarray(B.read_file8(wd)["H"]), ref)


# ======================================================================================
# 4. Global impedance
# ======================================================================================
def _rigid_T(xyz: np.ndarray, ref) -> np.ndarray:
    """Rigid-body kinematics from the cross product u = u0 + theta x (r - ref), built column by column."""
    d = np.asarray(xyz, float) - np.asarray(ref, float)
    n = d.shape[0]
    T = np.zeros((3 * n, 6))
    for j in range(n):
        T[3 * j:3 * j + 3, :3] = np.eye(3)
        for k in range(3):
            T[3 * j:3 * j + 3, 3 + k] = np.cross(np.eye(3)[k], d[j])
    return T


def _foun_rows(path):
    rows = [ln.split() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    a = np.array([[float(v) for v in r] for r in rows])
    return a[:, 0], a[:, 1:]


def test_global_impedance_is_the_inverse_compliance_of_a_rigid_massless_mat(tmp_path):
    """A massless mat made rigid by stiff BEAMS (1e5 G) to its centre: the 6x6 compliance at the
    centre from 6 simultaneous vibration load cases equals K_G^-1 (FILE11) up to the finite link
    stiffness; K_G is symmetric and equals T^T X T with an independently built T; the FOUN* files
    hold Re K, Im K / w, Im K / (2 |Re K|) and |K| (D-ANL-07)."""
    site = _site()
    mdl = B.surface_rigid_mat(site, half_width=2.0, ndiv=2)
    fs = B.FrequencySet.fourier(0.01, 1024, [4, 20, 60])
    B.run_soil(tmp_path, "m", site, fs, layer=0, rad=mdl.rad)
    B.run_house(tmp_path, "m", mdl)
    ref = (0.3, -0.2, 0.1)                                 # reference point off the centre
    assert _rc(tmp_path, fs=fs, impe=2, xc=ref[0], yc=ref[1], zc=ref[2]) == 0
    f11 = read_container(tmp_path / "FILE11")
    KG, X = np.asarray(f11["KG"]), np.asarray(f11["X"])
    T = _rigid_T(np.asarray(f11["int_xyz"]), ref)
    np.testing.assert_allclose(np.asarray(f11["T"]), T, rtol=0, atol=1e-14)
    for q in range(len(fs.fnum)):
        assert _relmax(T.T @ X[q] @ T, KG[q]) < 1e-12
        assert np.abs(KG[q] - KG[q].T).max() < 1e-12 * np.abs(KG[q]).max()
    # FOUN* text files
    w = 2 * math.pi * fs.freq
    f, stif = _foun_rows(tmp_path / "FOUNSTIF")
    np.testing.assert_allclose(f, fs.freq, rtol=1e-9)
    np.testing.assert_allclose(stif, KG.real.reshape(len(f), 36), rtol=1e-12, atol=1e-12 * np.abs(KG).max())
    _, dash = _foun_rows(tmp_path / "FOUNDASH")
    np.testing.assert_allclose(dash, (KG.imag / w[:, None, None]).reshape(len(f), 36), rtol=1e-12,
                               atol=1e-12 * np.abs(KG).max())
    _, impd = _foun_rows(tmp_path / "FOUNIMPD")
    np.testing.assert_allclose(impd, np.abs(KG).reshape(len(f), 36), rtol=1e-12)
    _, damp = _foun_rows(tmp_path / "FOUNDAMP")
    dg = np.arange(6) * 7                                  # diagonal positions in the row-major 6x6
    kd = np.diagonal(KG, axis1=1, axis2=2)
    np.testing.assert_allclose(damp[:, dg], kd.imag / (2 * np.abs(kd.real)), rtol=1e-12)
    # compliance of the rigid mat at its centre: reference point = centre for this comparison
    assert _rc(tmp_path, fs=fs, impe=2) == 0
    KG0 = np.asarray(read_container(tmp_path / "FILE11")["KG"])
    c = mdl["centre"]
    for d in range(1, 7):
        B.run_force(tmp_path, "m", [(c, d, 1.0, 0.0)], fs, copy_to=f"FILE9{d:03d}")
    assert _rc(tmp_path, fs=fs, type=1, simul=6) == 0
    G = np.zeros((len(fs.fnum), 6, 6), complex)
    for d in range(1, 7):
        f8 = B.read_file8(tmp_path, f"FILE8{d:03d}")
        for e in range(1, 7):
            G[:, e - 1, d - 1] = B.tf(f8, c, e)
    for q in range(len(fs.fnum)):
        Kinv = np.linalg.inv(KG0[q])
        assert _relmax(G[q], Kinv) < 1e-4            # finite stiffness of the rigid links (1e5 G)


# ======================================================================================
# 5. Restarts
# ======================================================================================
def test_restart_from_a_full_sparse_lu_record_equals_the_full_path_initiation(wd, monkeypatch):
    """Requirements 4.6 item 6: the fallback path stores the full sparse LU in COOTKqqq; a New
    Seismic Environment / New Dynamic Loading restart solves with the *restored* triangular factors
    and must reproduce the initiation that produced them (VP-22 tolerance 1e-10)."""
    wd, mdl = wd
    monkeypatch.setattr(ss, "RCOND_MIN", 2.0)                   # every C_nn 'ill-conditioned'
    assert _rc(wd, save=1) == 0
    f8 = B.read_file8(wd)
    assert np.all(np.asarray(f8["x_paths"]) == ss.PATHS.index("full"))
    ref = np.asarray(f8["H"]).copy()
    assert _rc(wd, mode=2) == 0
    assert _relmax(B.read_file8(wd)["H"], ref) < 1e-10
    B.run_force(wd, "m", [(mdl["top"], 1, 1.0, 0.0), (mdl["centre"], 3, 0.5, 0.05)], FS)
    assert _rc(wd, type=1, mode=3) == 0
    Hv = np.asarray(B.read_file8(wd)["H"]).copy()
    assert _rc(wd, type=1, mode=0) == 0
    assert _relmax(Hv, B.read_file8(wd)["H"]) < 1e-10


@pytest.mark.parametrize("case", ["new_structure", "new_seismic_environment", "new_dynamic_loading"])
def test_restart_records_of_another_frequency_step_are_not_used(wd, case):
    """Initiation with <save> = 1 on df = 1/(1024 x 0.01); then the free-field / load files are re-made
    with NFFT = 2048 (same frequency numbers, half the frequencies).  A restart must either refuse the
    records or equal a fresh initiation at the new frequencies."""
    wd, mdl = wd
    fs2 = B.FrequencySet.fourier(0.01, 2048, FS.fnum)
    typ = 1 if case == "new_dynamic_loading" else 0
    mode = {"new_structure": 1, "new_seismic_environment": 2, "new_dynamic_loading": 3}[case]
    loads = [(mdl["top"], 1, 1.0, 0.0)]
    if typ:
        B.run_force(wd, "m", loads, FS)
    assert _rc(wd, type=typ, save=1) == 0
    if typ:
        B.run_force(wd, "m", loads, fs2)
    else:
        B.run_soil(wd, "m", _site(), fs2, point=False)              # FILE1 on the new grid
    rc = _rc(wd, fs=fs2, type=typ, mode=mode)
    if rc != 0:
        assert "freq" in _listing(wd).lower()
        return
    H = np.asarray(B.read_file8(wd)["H"]).copy()
    B.run_soil(wd, "m", _site(), fs2, layer=0, rad=mdl.rad)        # SITE + POINT on the new grid
    assert _rc(wd, fs=fs2, type=typ, mode=0) == 0
    assert _relmax(H, B.read_file8(wd)["H"]) < 1e-8


# ======================================================================================
# 6. VP-39 accuracy through ANALYS
# ======================================================================================
def test_spring_mass_atf_peak_meets_the_vp39_tolerance_on_an_aligned_grid(tmp_path):
    """Requirements 6.3 VP-39: SC,1,1000,...,0.05 with MT = 1, ATF peak 1/(2b sqrt(1-b^2)) to 1e-8.
    With the SSI grid aligned on the peak frequency f0 sqrt(1 - 2b^2) (df = f_peak / 25000) the
    ANALYS chain on a practically rigid site (Vs = 100 km/s) reaches it; the sample there is the
    largest of its neighbours (the SASSI damping form puts the peak below f0)."""
    k, m, beta = 1000.0, 1.0, 0.05
    f0 = math.sqrt(k / m) / (2.0 * math.pi)
    fpk = f0 * math.sqrt(1.0 - 2.0 * beta ** 2)
    N = 25000
    fs = B.FrequencySet.harmonic(fpk / N, range(N - 3, N + 4))
    site = B.layered_site([(10.0, 1.0e5, 2.0e5, 2.0, 0.01)], (1.0e5, 2.0e5, 2.0, 0.01))
    B.run_soil(tmp_path, "m", site, fs, layer=0, rad=1.0)
    mdl = B.sdof_on_node(site, k, m, beta, element="spring")
    B.run_house(tmp_path, "m", mdl)
    assert _rc(tmp_path, fs=fs) == 0
    H = np.abs(B.tf(B.read_file8(tmp_path), mdl["mass"], 1))
    i = int(np.argmax(H))
    assert fs.fnum[i] == N
    exact = 1.0 / (2.0 * beta * math.sqrt(1.0 - beta ** 2))
    assert abs(H[i] / exact - 1.0) < 1e-8
