"""Adversarial review tests of the 'motion' work package (MOTION, RELDISP, COMBIN, interp, signal).

These tests probe the implementation from angles that differ from the implementer's own tests:
independent derivations (time-domain SDOF integration, explicit 5x5 Tajirian solves), invariance
properties (scaling, frequency scaling, column independence, window locality), symmetry
(RELDISP antisymmetry / transitivity, sign equivariance), round trips (COMBIN -> MOTION) and
limiting cases.  Tests marked ``xfail(strict=True)`` document defects found in the review; they
turn into XPASS failures once the defect is fixed, so the marker must then be removed.

Requirements cited: section 4.8 (MOTION), 4.9 (RELDISP), 4.7 (COMBIN), D-MOT-01..05, D-MOT-08/09,
D-CNV-06, D-RDP-01, VP-28 statements of requirements section 6.3, R1 section 5.1.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest
from scipy import signal as sg

from sassi.conventions import cfactor, nodal_rs_name, nodal_result_name
from sassi.core import interp as I
from sassi.core import signal as S
from sassi.io import textfiles
from sassi.io.container import read_container
from sassi.modules.base import run_module
from sassi.verify.problems.vp_motion import (THREE_DOF, TWO_DOF, motion_deck, reldisp_deck, run_deck,
                                             shear_chain_tf, synthetic_motion, write_file8, write_motion_file)

# grid of the VPs: NFFT 4096, dt 0.005 s -> df = 0.048828125 Hz
NFFT4, DT4 = 4096, 0.005
DF4 = 1.0 / (NFFT4 * DT4)
LOGSET = np.unique(np.round(np.geomspace(4, 512, 64)).astype(np.int64))


def _hyst_sdof(f0: float, beta: float):
    ks = (2 * np.pi * f0) ** 2 * cfactor(beta)
    return lambda f: ks / (ks - (2 * np.pi * np.asarray(f, dtype=float)) ** 2)


def _prel(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.max(np.abs(a - b) / np.abs(b)))


def _rel(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


# =======================================================================================
# Interpolation core: independent derivations and invariances
# =======================================================================================
def _tajirian_window(fw, Hw, f):
    """Independent Tajirian fit (R1 section 5.1, minus sign in column 5) by plain LU of the 5x5 system.

    Solved in x = (w / w_max)^2, which leaves the interpolant unchanged (D-MOT-01).
    """
    x = (fw / fw[-1]) ** 2
    A = np.column_stack([x ** 2, x, np.ones(5), -x * Hw, -Hw])
    c = np.linalg.solve(A, x ** 2 * Hw)
    xx = (np.asarray(f) / fw[-1]) ** 2
    return (c[0] * xx ** 2 + c[1] * xx + c[2]) / (xx ** 2 + c[3] * xx + c[4])


def test_dense_options_equal_explicit_window_combinations():
    """Options 2, 0, 3 equal mean / weighted mean / nearest-3 mean of independently fitted windows."""
    fn = LOGSET
    fs = fn * DF4
    H = shear_chain_tf(fs, **THREE_DOF)[:, 2]                  # 3-DOF: not of the 2-DOF form, full rank
    N = len(fs)
    rng = np.random.default_rng(1)
    j = rng.integers(0, N - 1, 25)                              # 0-based intervals
    f = fs[j] + rng.uniform(0.1, 0.9, 25) * (fs[j + 1] - fs[j])
    ref = {0: [], 2: [], 3: []}
    for jj, ff in zip(j, f):
        j1 = jj + 1
        ms = range(max(1, j1 - 3), min(j1, N - 4) + 1)
        vals = np.array([_tajirian_window(fs[m - 1:m + 4], H[m - 1:m + 4], ff) for m in ms])
        c = np.array([0.5 * (fs[m - 1] + fs[m + 3]) for m in ms])
        h = np.array([0.5 * (fs[m + 3] - fs[m - 1]) for m in ms])
        w = np.maximum(1e-3, 1.0 - np.abs(ff - c) / h)
        ref[2].append(vals.mean())
        ref[0].append(np.sum(w * vals) / np.sum(w))
        near = np.argsort(np.abs(ff - c), kind="stable")[:3]
        ref[3].append(vals[near].mean())
    for opt in (0, 2, 3):
        got = I.interpolate_tf(fs, H, f, opt, h0=1.0)
        assert _prel(got, np.array(ref[opt])) < 1e-9, opt


@pytest.mark.parametrize("option", [1, 4, 5])
def test_single_window_locality(option):
    """Options 1/4/5: the value in interval j depends only on the points of its window (section 4.8 table)."""
    fn = LOGSET
    fs = fn * DF4
    H = shear_chain_tf(fs, **THREE_DOF)[:, 1]
    N = len(fs)
    s0 = {1: 1, 4: 2, 5: 3}[option]
    for j1 in (1, 2, 3, 7, 22, N - 3, N - 1):
        m = 1 if j1 < s0 else s0 + 4 * ((j1 - s0) // 4)
        m = min(m, N - 4)
        assert m <= j1 <= m + 3                                 # the window contains interval j1
        f = fs[j1 - 1] + 0.37 * (fs[j1] - fs[j1 - 1])
        base = I.interpolate_tf(fs, H, [f], option, h0=1.0)[0]
        Hp = H.copy()
        outside = np.ones(N, bool)
        outside[m - 1:m + 4] = False
        Hp[outside] *= (1.3 - 0.4j)                             # perturb every point outside the window
        pert = I.interpolate_tf(fs, Hp, [f], option, h0=1.0)[0]
        assert abs(pert - base) <= 1e-12 * abs(base), (option, j1)
        alone = _tajirian_window(fs[m - 1:m + 4], H[m - 1:m + 4], f)
        assert abs(alone - base) <= 1e-9 * abs(base), (option, j1)


def test_scale_frequency_and_column_invariance():
    """Linearity in H, invariance under f -> c f, and independence of the columns (all options)."""
    fs = LOGSET * DF4
    H = shear_chain_tf(fs, **THREE_DOF)                          # (nF, 3)
    fk = S.fourier_grid(NFFT4, DT4)
    alpha = 3.7 - 2.2j
    for opt in range(7):
        a = I.interpolate_tf(fs, H, fk, opt, h0=np.ones(3))
        b = I.interpolate_tf(fs, alpha * H, fk, opt, h0=alpha * np.ones(3))
        assert _rel(b, alpha * a) < 1e-11, opt
        c = I.interpolate_tf(2.0 * fs, H, 2.0 * fk, opt, h0=np.ones(3))
        assert _rel(c, a) < 1e-11, opt
        for k in range(3):
            single = I.interpolate_tf(fs, H[:, k], fk, opt, h0=1.0)
            assert np.allclose(single, a[:, k], rtol=0, atol=1e-13 * np.max(np.abs(a[:, k]))), (opt, k)


def test_options_0_5_exact_for_another_lightly_damped_2dof():
    """VP-28 statement with a different 2-DOF system (closer modes, 2 % damping) on the VP log grid."""
    chain = dict(masses=(2.0, 1.0), stiffs=(3.0 * (2 * np.pi * 2.0) ** 2, (2 * np.pi * 2.6) ** 2), beta=0.02)
    fs = LOGSET * DF4
    fk = S.fourier_grid(NFFT4, DT4)
    ins = (fk >= fs[0]) & (fk <= fs[-1])
    H = shear_chain_tf(fs, **chain)
    ex = shear_chain_tf(fk[ins], **chain)
    for opt in range(6):
        Hi = I.interpolate_tf(fs, H, fk, opt, h0=np.ones(2))
        assert _prel(Hi[ins], ex) < 1e-10, opt


def test_guard_root_cause_rational_fit_is_exact_in_guarded_window():
    """Documents the root cause of the guard defect: the 5-point rational fit of an exact hysteretic SDOF
    is exact (1e-13) in the window that the literal D-MOT-01 ratio test flags (true resonance pole).

    Fix (implementer): the ratio test is kept as the *trigger* (``fit["trigger"]``); the guard itself
    (``fit["guard"]``) no longer fires on a physically admissible (damped) pole."""
    sd = _hyst_sdof(0.5, 0.02)
    fs = np.arange(4, 513, 8) * DF4                               # uniform 0.39 Hz SSI spacing from 0.195 Hz
    fit = I.fit_windows(fs, sd(fs)[:, None])
    assert bool(fit["trigger"][0, 0])                             # the ratio test fires in window 1 ...
    assert not bool(fit["guard"][0, 0])                           # ... but the refined guard keeps the rational form
    C = fit["C"][0, 0]
    poles = np.roots([1.0, C[3], C[4]])
    true_pole = (0.5 / fs[4]) ** 2 * cfactor(0.02)                # ... on the genuine resonance pole
    assert np.min(np.abs(poles - true_pole)) < 1e-10 * abs(true_pole)
    fk = S.fourier_grid(NFFT4, DT4)
    sel = (fk >= fs[0]) & (fk <= fs[4])
    raw = I._rational_eval(fit, np.zeros(int(sel.sum()), dtype=np.int64), fk[sel])[:, 0]
    assert _prel(raw, sd(fk[sel])) < 1e-12


def test_guard_does_not_destroy_exact_lightly_damped_sdof():
    """VP-28: 'hysteretic SDOF windows rank-deficient but exact' -- also on a uniform SSI grid."""
    sd = _hyst_sdof(0.5, 0.02)
    fs = np.arange(4, 513, 8) * DF4
    fk = S.fourier_grid(NFFT4, DT4)
    ins = (fk >= fs[0]) & (fk <= fs[-1])
    for opt in range(6):
        Hi = I.interpolate_tf(fs, sd(fs), fk, opt, h0=1.0)
        assert _prel(Hi[ins], sd(fk[ins])) < 1e-10, opt


@pytest.mark.parametrize("option,smooth", [(6, 0.0), (2, 1.0), (2, 1000.0)])
def test_phase_adjustment_is_sign_equivariant(option, smooth):
    """Reversing the input direction (ang = 180 deg, H -> -H, anchor -1) must reverse the adjusted TF."""
    fs = LOGSET * DF4
    fk = S.fourier_grid(NFFT4, DT4)
    H = shear_chain_tf(fs, **TWO_DOF)[:, 1]
    assert I.rigid_body_anchor(1, cm=0, ang_deg=180.0) == -1.0
    pos = I.interpolate_tf(fs, H, fk, option, smooth=smooth, pzadj=1, h0=1.0)
    neg = I.interpolate_tf(fs, -H, fk, option, smooth=smooth, pzadj=1, h0=-1.0)
    assert abs(neg[0] - (-1.0)) < 1e-12                           # rigid-body limit kept
    assert _rel(neg, -pos) < 1e-12


def test_smoothing_zero_is_bitwise_identity_and_large_s_bounds_overshoot():
    """C14: S = 0 identical to no smoothing; S -> large keeps |H| inside the neighbour band (3-DOF data)."""
    fs = LOGSET[::2] * DF4                                         # coarse grid -> overshoots exist
    H = shear_chain_tf(fs, **THREE_DOF)[:, 2]
    fk = S.fourier_grid(NFFT4, DT4)
    for opt in range(6):
        a = I.interpolate_tf(fs, H, fk, opt, h0=1.0)
        b = I.interpolate_tf(fs, H, fk, opt, smooth=0.0, h0=1.0)
        assert np.array_equal(a, b)
    ins = (fk >= fs[0]) & (fk <= fs[-1])
    f = fk[ins]
    Hi = I.interpolate_tf(fs, H, f, 2, h0=1.0)
    Hs = I.interpolate_tf(fs, H, f, 2, smooth=1e9, h0=1.0)
    j = np.clip(np.searchsorted(fs, f, side="right") - 1, 0, len(fs) - 2)
    amax = np.maximum(np.abs(H[j]), np.abs(H[j + 1]))
    amin = np.minimum(np.abs(H[j]), np.abs(H[j + 1]))
    t = (f - fs[j]) / (fs[j + 1] - fs[j])
    L = H[j] + (H[j + 1] - H[j]) * t                               # complex linear reference of D-MOT-04
    inside = (np.abs(Hi) >= amin * (1 + 1e-9)) & (np.abs(Hi) <= amax * (1 - 1e-9))
    outside = (np.abs(Hi) > amax * (1 + 1e-6)) | (np.abs(Hi) < amin * (1 - 1e-6))
    assert inside.any() and outside.any()
    assert np.allclose(Hs[inside], Hi[inside], rtol=1e-14, atol=0)    # inside the band: untouched (to round-off)
    assert np.all(np.abs(Hs[outside] - L[outside]) <= 1e-3 * np.abs(Hi[outside] - L[outside]))   # pulled onto L


# =======================================================================================
# Signal core
# =======================================================================================
def test_hudson_housner_residual_orthogonality():
    """D-MOT-08: the corrected velocity is the least-squares residual -> orthogonal to t and t^2."""
    dt = 0.01
    a = 0.2 * synthetic_motion(1500, dt, seed=17) + 0.003 + 0.0004 * np.arange(1500) * dt
    br = S.baseline_correction(a, dt, scale=9.81)
    t = np.arange(len(a)) * dt
    v = br.vel
    for basis in (t, 0.5 * t * t):
        assert abs(np.dot(v, basis)) <= 1e-10 * np.linalg.norm(v) * np.linalg.norm(basis)
    # the removed baseline is exactly the difference of the accelerations
    assert np.allclose(a - br.acc, br.c1 + br.c2 * t, rtol=0, atol=1e-15)


def test_displacement_spectrum_double_integration_of_harmonic_bin():
    """D-CNV-06: U_g = -g A / w^2 recovers the displacement of a periodic harmonic acceleration."""
    nfft, dt, g = 512, 0.02, 9.81
    t = np.arange(nfft) * dt
    k = 7
    w = 2 * np.pi * k / (nfft * dt)
    d_true = 0.013 * np.sin(w * t + 0.4)                          # length units
    acc_g = -w ** 2 * d_true / g
    U = S.displacement_spectrum(np.fft.rfft(acc_g), S.fourier_grid(nfft, dt), g)
    assert _rel(np.fft.irfft(U, nfft), d_true) < 1e-12


# =======================================================================================
# MOTION through the module: independent time-domain physics
# =======================================================================================
def _foh_response(A, B, C, D, u, h):
    """Exact first-order-hold discretisation of x' = Ax + Bu, y = Cx + Du (independent of sassi)."""
    Ad, Bd, Cd, Dd, _ = sg.cont2discrete((A, B, C, D), h, method="foh")
    b, a = sg.ss2tf(Ad, Bd, Cd, Dd)
    return sg.lfilter(b[0], a, u)


def test_motion_matches_time_domain_viscous_sdof(tmp_path):
    """MOTION .ACC (e^{+iwt} TF, convolution with A in g) equals a direct time-domain solution.

    Viscous SDOF on a moving base: x'' + 2 z w0 x' + w0^2 x = -a_g, total acceleration
    a = -(2 z w0 x' + w0^2 x), i.e. H = (w0^2 + 2 i z w0 w) / (w0^2 - w^2 + 2 i z w0 w) under e^{+iwt}.
    """
    nfft, dt = 4096, 0.0025
    df = 1.0 / (nfft * dt)
    f0, z = 3.0, 0.05
    w0 = 2 * np.pi * f0
    fn = np.arange(1, nfft // 2 + 1)
    w = 2 * np.pi * fn * df
    H = (w0 ** 2 + 2j * z * w0 * w) / (w0 ** 2 - w ** 2 + 2j * z * w0 * w)
    write_file8(tmp_path / "FILE8", fn, df, [7], [1], H[:, None], nfft=nfft, delt=dt)
    ag = 0.3 * synthetic_motion(1600, dt, seed=2, fmax=12.0)
    write_motion_file(tmp_path / "eq.acc", ag, dt)
    d = motion_deck(nout=[(7, 1, 0, 1, 0, 0, 0, 1)], thfile="eq.acc", nft=nfft, delt=dt, interp=1)
    assert run_deck(tmp_path, "m", d) == 0
    acc, dto = textfiles.read_history(tmp_path / "00007TR_X.ACC")
    assert dto == dt and len(acc) == nfft
    up = 8                                                        # band-limited input, 8x finer time step
    u = sg.resample(S.pad_record(ag, nfft), nfft * up)
    A = np.array([[0.0, 1.0], [-w0 ** 2, -2 * z * w0]])
    y = _foh_response(A, np.array([[0.0], [-1.0]]), np.array([[-w0 ** 2, -2 * z * w0]]), np.zeros((1, 1)),
                      u, dt / up)[::up]
    assert _rel(acc, y) < 1e-3
    # the opposite time convention (e^{-iwt}) would be grossly wrong: sanity of the test itself
    yc = np.fft.irfft(np.r_[1.0, np.conj(H)] * np.fft.rfft(S.pad_record(ag, nfft)), nfft)
    assert _rel(yc, y) > 0.3


def test_motion_rigid_rotated_input_reproduces_projected_control_motion(tmp_path):
    """H = cos(a) at all frequencies (rigid body, ang = 150 deg): .ACC = cos(a) * input (D-MOT-03 anchor)."""
    nfft, dt = 1024, 0.01
    df = 1.0 / (nfft * dt)
    ang = 150.0
    ca, sa = np.cos(np.radians(ang)), np.sin(np.radians(ang))
    fn = np.arange(1, nfft // 2 + 1)
    H = np.column_stack([np.full(len(fn), ca), np.full(len(fn), sa)])
    write_file8(tmp_path / "FILE8", fn, df, [4, 4], [1, 2], H, nfft=nfft, delt=dt, ang=ang)
    ag = 0.2 * synthetic_motion(700, dt, seed=12, fmax=30.0)
    write_motion_file(tmp_path / "eq.acc", ag, dt)
    d = motion_deck(nout=[(4, 1, 1, 1, 0, 0, 0, 1), (4, 2, 1, 1, 0, 0, 0, 1)], thfile="eq.acc", nft=nfft, delt=dt,
                    ang=ang, cplx=1)
    assert run_deck(tmp_path, "m", d) == 0
    ax, _ = textfiles.read_history(tmp_path / "00004TR_X.ACC")
    ay, _ = textfiles.read_history(tmp_path / "00004TR_Y.ACC")
    a = S.pad_record(ag, nfft)
    assert _rel(ax, ca * a) < 1e-12 and _rel(ay, sa * a) < 1e-12
    f, hx, _ = textfiles.read_tf(tmp_path / "00004TR_X.TFI")
    assert f[0] == 0.0 and abs(hx[0] - ca) < 1e-14                  # TFI f = 0 row = rigid-body anchor


def test_rs_uses_full_record_and_tends_to_zpa(tmp_path):
    """D-MOT-09: RS over the full NFFT record whatever dur; SA at 100 Hz ~ PGA of the response (ZPA)."""
    nfft, dt = 2048, 0.005
    df = 1.0 / (nfft * dt)
    fn = np.unique(np.round(np.geomspace(2, 1024, 80)).astype(int))
    H = np.column_stack([np.ones(len(fn)), shear_chain_tf(fn * df, **TWO_DOF)])
    write_file8(tmp_path / "FILE8", fn, df, [1, 2, 3], [1, 1, 1], H, nfft=nfft, delt=dt)
    write_motion_file(tmp_path / "eq.acc", 0.25 * synthetic_motion(1500, dt, seed=31, fmax=20.0), dt)
    rs = {}
    for dur in (0.0, 1.0):
        sub = tmp_path / f"d{dur:g}"
        sub.mkdir()
        for nm in ("FILE8", "eq.acc"):
            shutil.copy(tmp_path / nm, sub / nm)
        d = motion_deck(nout=[(3, 1, 0, 1, 0, 0, 1, 1)], damp=[0.05], thfile="eq.acc", nft=nfft, delt=dt,
                        freq1=0.1, freq2=100.0, fstep=301, dur=dur)
        assert run_deck(sub, "m", d) == 0
        rs[dur] = textfiles.read_xy(sub / nodal_rs_name(3, 1, 1))[:, 1]
    acc_full, _ = textfiles.read_history(tmp_path / "d0" / "00003TR_X.ACC")
    acc_short, _ = textfiles.read_history(tmp_path / "d1" / "00003TR_X.ACC")
    assert len(acc_short) == S.output_length(nfft, dt, 1.0) < len(acc_full)
    assert np.array_equal(rs[0.0], rs[1.0])
    assert abs(rs[0.0][-1] / np.max(np.abs(acc_full)) - 1.0) < 0.02


def test_vibration_velocity_and_acceleration_are_time_derivatives(tmp_path):
    """Vibration (type 1): resp 1 = d/dt resp 0 and resp 2 = d/dt resp 1 (central differences)."""
    nfft, dt = 2048, 0.005
    df = 1.0 / (nfft * dt)
    fn = np.unique(np.round(np.geomspace(1, 1024, 90)).astype(int))
    w = 2 * np.pi * fn * df
    k, m, c = 4.0e3, 1.0, 2 * 0.05 * np.sqrt(4.0e3)
    H = 1.0 / (k - m * w ** 2 + 1j * c * w)                       # displacement per unit load
    write_file8(tmp_path / "FILE8", fn, df, [9], [3], H[:, None], nfft=nfft, delt=dt, type=1)
    write_motion_file(tmp_path / "load.acc", synthetic_motion(800, dt, seed=5, fmax=8.0), dt)
    out = {}
    for resp in (0, 1, 2):
        sub = tmp_path / f"r{resp}"
        sub.mkdir()
        for nm in ("FILE8", "load.acc"):
            shutil.copy(tmp_path / nm, sub / nm)
        d = motion_deck(nout=[(9, 3, 0, 1, 0, 0, 0, 1)], thfile="load.acc", nft=nfft, delt=dt, type=1, resp=resp)
        assert run_deck(sub, "v", d) == 0
        out[resp], _ = textfiles.read_history(sub / "00009TR_Z.ACC")
    for lo, hi in ((0, 1), (1, 2)):
        deriv = (np.roll(out[lo], -1) - np.roll(out[lo], 1)) / (2 * dt)   # periodic record
        assert _rel(deriv, out[hi]) < 2e-2


def test_srss_tfu_is_root_sum_square_of_modal_tfs(tmp_path):
    """D-MOT-06 / B.5.4: TFU amplitude sqrt(sum |H_k|^2), zero phase (p = 0), at the computed frequencies."""
    nfft, dt = 1024, 0.01
    df = 1.0 / (nfft * dt)
    fn = np.unique(np.round(np.geomspace(2, 300, 40)).astype(int))
    H1 = shear_chain_tf(fn * df, **TWO_DOF)
    H2 = 0.3 * (1 - 0.5j) * shear_chain_tf(fn * df, **TWO_DOF)[:, ::-1]
    for nm, Hm in (("FILE8_01", H1), ("FILE8_02", H2)):
        write_file8(tmp_path / nm, fn, df, [2, 3], [1, 1], Hm, nfft=nfft, delt=dt)
    write_motion_file(tmp_path / "eq.acc", 0.2 * synthetic_motion(500, dt, seed=4), dt)
    (tmp_path / "srsstf.TXT").write_text("2 0\nFILE8_01\nFILE8_02\n", encoding="utf-8")   # case-insensitive
    d = motion_deck(nout=[(3, 1, 1, 0, 0, 0, 0, 0)], thfile="eq.acc", nft=nfft, delt=dt, srss=1, cplx=1)
    assert run_deck(tmp_path, "s", d) == 0
    f, Hu, cplx = textfiles.read_tf(tmp_path / "00003TR_X.TFU")
    assert cplx and np.allclose(f, fn * df, rtol=0, atol=1e-9)
    ref = np.sqrt(np.abs(H1[:, 1]) ** 2 + np.abs(H2[:, 1]) ** 2)
    assert _rel(Hu, ref) < 1e-14 and np.all(np.angle(Hu) == 0.0)


# =======================================================================================
# RELDISP symmetry and COMBIN round trip
# =======================================================================================
def _motion_three_nodes(wd: Path, nfft=1024, dt=0.01, g=9.81):
    df = 1.0 / (nfft * dt)
    fn = np.unique(np.concatenate([np.arange(1, 6), np.round(np.geomspace(6, 400, 50))]).astype(int))
    Hc = shear_chain_tf(fn * df, **TWO_DOF)
    H = np.column_stack([np.ones(len(fn)), Hc[:, 0], Hc[:, 1]])
    write_file8(wd / "FILE8", fn, df, [1, 2, 3], [1, 1, 1], H, nfft=nfft, delt=dt)
    write_motion_file(wd / "eq.acc", 0.2 * synthetic_motion(600, dt, seed=8, fmax=15.0), dt)
    d = motion_deck(nout=[(n, 1, 1, 1, 0, 0, 0, 0) for n in (1, 2, 3)], thfile="eq.acc", nft=nfft, delt=dt, cplx=1,
                    gravity=g)
    assert run_deck(wd, "m", d) == 0
    return dict(thfile="eq.acc", nft=nfft, delt=dt, gravity=g)


def test_reldisp_antisymmetry_and_transitivity(tmp_path):
    """d(3 rel 2) = -d(2 rel 3) and d(3 rel FF) - d(2 rel FF) = d(3 rel 2) (linearity of D-RDP-01)."""
    common = _motion_three_nodes(tmp_path)
    runs = {}
    for name, ref, node in (("a", "00002TR_X.TFI", 3), ("b", "00003TR_X.TFI", 2), ("c", "FREEFIELD", 3),
                            ("d", "FREEFIELD", 2)):
        sub = tmp_path / name
        sub.mkdir()
        for p in tmp_path.iterdir():
            if p.is_file():
                shutil.copy(p, sub / p.name)
        assert run_deck(sub, "r", reldisp_deck(ref, [(node, 1, 0, 0, 0, 0, 0)], **common)) == 0
        runs[name], _ = textfiles.read_history(sub / nodal_result_name(node, 1, "THD"))
    scale = np.max(np.abs(runs["a"]))
    assert np.max(np.abs(runs["a"] + runs["b"])) <= 1e-12 * scale
    assert np.max(np.abs((runs["c"] - runs["d"]) - runs["a"])) <= 1e-10 * scale


def test_reldisp_rigid_node_with_rotated_free_field_reference_is_zero(tmp_path):
    """ang = 180 deg: a node moving rigidly with the ground (H = -1 for X up to Nyquist) has zero
    relative displacement with respect to the free-field reference (anchor projection -1)."""
    nfft, dt, g = 1024, 0.01, 9.81
    df = 1.0 / (nfft * dt)
    fn = np.arange(1, nfft // 2 + 1)
    write_file8(tmp_path / "FILE8", fn, df, [5], [1], -np.ones((len(fn), 1)), nfft=nfft, delt=dt, ang=180.0)
    write_motion_file(tmp_path / "eq.acc", 0.2 * synthetic_motion(600, dt, seed=8), dt)
    d = motion_deck(nout=[(5, 1, 1, 0, 0, 0, 0, 0)], thfile="eq.acc", nft=nfft, delt=dt, cplx=1, ang=180.0, gravity=g)
    assert run_deck(tmp_path, "m", d) == 0
    rd = reldisp_deck("FREEFIELD", [(5, 1, 0, 0, 0, 0, 0)], thfile="eq.acc", nft=nfft, delt=dt, gravity=g, ang=180.0)
    assert run_deck(tmp_path, "r", rd) == 0
    dd, _ = textfiles.read_history(tmp_path / "00005TR_X.THD")
    assert np.max(np.abs(dd)) < 1e-14


def test_combin_then_motion_equals_full_run(tmp_path):
    """COMBIN (odd + even frequency numbers) followed by MOTION reproduces MOTION on the full FILE8."""
    nfft, dt = 1024, 0.01
    df = 1.0 / (nfft * dt)
    fn = np.arange(2, 301, 3)
    H = np.column_stack([np.ones(len(fn)), shear_chain_tf(fn * df, **TWO_DOF)])
    acc = 0.2 * synthetic_motion(600, dt, seed=3)
    full, comb = tmp_path / "full", tmp_path / "comb"
    full.mkdir()
    comb.mkdir()
    write_file8(full / "FILE8", fn, df, [1, 2, 3], [1, 1, 1], H, nfft=nfft, delt=dt)
    odd = (np.arange(len(fn)) % 2) == 1
    write_file8(comb / "FILE82", fn[odd], df, [1, 2, 3], [1, 1, 1], H[odd], nfft=nfft, delt=dt)   # 82 first in f
    write_file8(comb / "FILE81", fn[~odd], df, [1, 2, 3], [1, 1, 1], H[~odd], nfft=nfft, delt=dt)
    assert run_module("COMBIN", "c", comb) == 0
    merged = read_container(comb / "FILE8", "FILE8")
    assert np.array_equal(merged["fnum"], fn) and np.array_equal(merged["H"], H)
    for wd in (full, comb):
        write_motion_file(wd / "eq.acc", acc, dt)
        d = motion_deck(nout=[(3, 1, 1, 1, 0, 0, 0, 1)], thfile="eq.acc", nft=nfft, delt=dt, interp=0, cplx=1)
        assert run_deck(wd, "m", d) == 0
    for ext in ("TFI", "ACC"):
        assert (full / f"00003TR_X.{ext}").read_text() == (comb / f"00003TR_X.{ext}").read_text()
