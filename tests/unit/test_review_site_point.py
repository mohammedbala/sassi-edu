"""Adversarial review tests of the SITE / POINT work package (layered-soil core).

These tests probe the implementation from angles that differ from the implementer's own tests:

* independent re-derivations (Waas-form layer matrices, R1 §2.3; Kausel-Roesset closed-form
  half-space stiffness, D-SIT-08; a Thomson-Haskell type potential propagator for inclined P-SV
  waves through layered soil, R1 §2.7b; an independent 1-D column and companion-form TLM solver
  for the VP-02b / VP-05 gate numbers);
* structural properties (rigid-body null spaces and mass of the POINT core elements, complex
  symmetry of transmitting boundaries, passivity Im(F_ff) <= 0 under exp(+i w t), mirror /
  translation / permutation covariance of F_ff, POINT2 reciprocity);
* physical limiting cases (Rayleigh-wave ellipticity and retrograde motion, within TF
  independent of the half-space model, unit control motion of a mixed wave field through the
  binding API);
* robustness (frequency-subset consistency of FILE3, input checks).

Tests marked ``xfail(strict=True)`` document defects found by the review.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy import linalg as sla
from scipy.optimize import brentq

from sassi.conventions import cfactor
from sassi.core import axisym, strip2d, tlm
from sassi.core import freefield as ff
from sassi.core.flexibility import flexibility_matrix, frequency_row
from sassi.io import decks
from sassi.io.files import read_container
from sassi.modules.base import run_module

GRAV = 9.81


@pytest.fixture(autouse=True)
def _quiet_fpe():
    # numpy 2.0 + Accelerate raises spurious FP flags in matmul (see the implementer's report)
    with np.errstate(all="ignore"):
        yield


# ---------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------
def _site(wd: Path, layers, hs, fnum, df, nl=20, waves=((2, 1, 1.0, 1.0, 0.0),), cl=1, cm=0, wopt=0,
          mode2=True, hslaw=None, **params):
    """layers / hs: (thick, rho, vs, vp, beta_s, beta_p)."""
    wd.mkdir(parents=True, exist_ok=True)
    d = decks.new("SITE")
    d["gravity"] = GRAV
    d["nl"] = nl
    d["df"] = df
    d["cl"] = cl
    d["cm"] = cm
    d["wopt"] = wopt
    d["mode2"] = int(mode2)
    if hslaw is not None:
        d["hslaw"] = hslaw
    for k, v in params.items():
        d[k] = v
    for i, (h, rho, vs, vp, bs, bp) in enumerate(layers):
        d.table("layers").append([i + 1, h, rho * GRAV, vp, vs, bp, bs])
    _, rho, vs, vp, bs, bp = hs
    d.table("halfspace").append([len(layers) + 1, 0.0, rho * GRAV, vp, vs, bp, bs])
    for w in waves:
        d.table("waves").append(list(w))
    for n in fnum:
        d.table("freqs").append([int(n)])
    decks.write(wd / "m.sit", d)
    rc = run_module("SITE", "m", wd)
    return rc, (wd / "m_SITE.out").read_text()


def _point(wd: Path, layer, rad, df, fnum=(), dim=2):
    d = decks.new("POINT")
    d["layer"] = layer
    d["rad"] = rad
    d["dim"] = dim
    d["df"] = df
    for n in fnum:
        d.table("freqs").append([int(n)])
    decks.write(wd / "m.poi", d)
    rc = run_module("POINT", "m", wd)
    return rc, (wd / "m_POINT.out").read_text()


def _vp(vs, nu):
    return vs * np.sqrt(2.0 * (1.0 - nu) / (1.0 - 2.0 * nu))


def _test_column(base_dashpot=True):
    h = np.array([0.5, 0.7, 1.0, 1.3, 2.0])
    rho = np.array([1.9, 2.0, 2.1, 2.2, 2.2])
    G, M = tlm.layer_moduli(rho, [150, 200, 250, 300, 350], [300, 420, 500, 640, 700], [0.05] * 5, [0.03] * 5)
    hG, hM = tlm.layer_moduli(2.3, 500.0, 1000.0, 0.02, 0.02)
    hg = [3.0, 4.0, 5.0, 6.0] if base_dashpot else None
    return tlm.build_column(h, rho, G, M, 2.3, complex(hG), complex(hM), hg)


# ---------------------------------------------------------------------------------------
# 1. Thin-layer eigen-solutions vs an independent (Waas-form) assembly, R1 §2.2-2.3
# ---------------------------------------------------------------------------------------
def _waas_inplane(col: tlm.Column, omega: float, k: complex, mass=(5 / 12, 1 / 12)) -> np.ndarray:
    """Waas/SASSI form (A k^2 + i B_W k + G - w^2 M) on {u, v_down} per interface (R1 §2.3),
    assembled here from scratch (interleaved order), with the base dashpots."""
    n = col.n_iface
    K = np.zeros((2 * n, 2 * n), complex)
    lam = col.M - 2 * col.G
    for j in range(col.n_sub):
        L, G, h, r = lam[j], col.G[j], col.h[j], col.rho[j]
        lp = L + 2 * G
        A = h / 6 * np.array([[2 * lp, 0, lp, 0], [0, 2 * G, 0, G], [lp, 0, 2 * lp, 0], [0, G, 0, 2 * G]])
        BW = 0.5 * np.array([[0, -(L - G), 0, L + G], [L - G, 0, L + G, 0],
                             [0, -(L + G), 0, L - G], [-(L + G), 0, -(L - G), 0]])
        Gg = 1 / h * np.array([[G, 0, -G, 0], [0, lp, 0, -lp], [-G, 0, G, 0], [0, -lp, 0, lp]])
        a, b = mass
        Mm = r * h * np.array([[a, 0, b, 0], [0, a, 0, b], [b, 0, a, 0], [0, b, 0, a]])
        idx = np.arange(2 * j, 2 * j + 4)
        K[np.ix_(idx, idx)] += A * k * k + 1j * BW * k + Gg - omega ** 2 * Mm
    if col.base == tlm.BASE_DASHPOT:
        K[2 * n - 2, 2 * n - 2] += 1j * omega * col.cs
        K[2 * n - 1, 2 * n - 1] += 1j * omega * col.cp
        return K
    return K[:2 * (n - 1), :2 * (n - 1)]


def test_rayleigh_modes_satisfy_independent_waas_quadratic_problem():
    """Every one of the 2Nf Rayleigh modes of tlm.column_modes, mapped back to physical
    Waas variables (v_down = -u_z = i phi_z), must annihilate the independently assembled Waas
    quadratic pencil; the opposite vertical sign must not (checks B sign and u_z = -i phi_z)."""
    col = _test_column()
    om = 2 * np.pi * 7.0
    md = tlm.column_modes(col, om)
    worst, wrong = 0.0, np.inf
    for j, k in enumerate(md.kR):
        V = np.zeros(2 * col.n_free, complex)
        V[0::2] = md.phix[:, j]
        V[1::2] = 1j * md.phiz[:, j]
        Kw = _waas_inplane(col, om, k)
        worst = max(worst, np.linalg.norm(Kw @ V) / (np.linalg.norm(Kw) * np.linalg.norm(V)))
        V[1::2] *= -1
        wrong = min(wrong, np.linalg.norm(Kw @ V) / (np.linalg.norm(Kw) * np.linalg.norm(V)))
    assert worst < 1e-11
    assert wrong > 1e-4


def test_love_modes_satisfy_independent_sh_problem():
    col = _test_column(base_dashpot=False)
    om = 2 * np.pi * 9.0
    md = tlm.column_modes(col, om)
    n = col.n_iface
    for j, k in enumerate(md.kL):
        K = np.zeros((n, n), complex)
        for s in range(col.n_sub):
            G, h, r = col.G[s], col.h[s], col.rho[s]
            K[s:s + 2, s:s + 2] += (G * h / 6 * np.array([[2, 1], [1, 2]]) * k * k
                                    + G / h * np.array([[1, -1], [-1, 1]])
                                    - om ** 2 * r * h * np.array([[5 / 12, 1 / 12], [1 / 12, 5 / 12]]))
        K = K[:-1, :-1]                         # rigid base
        v = md.phiy[:, j]
        assert np.linalg.norm(K @ v) <= 1e-10 * np.linalg.norm(K) * np.linalg.norm(v)
    # outgoing branch for every mode
    assert np.all(md.kL.imag <= 1e-12 * np.abs(md.kL))
    assert np.all(md.kR.imag <= 1e-12 * np.abs(md.kR))


# ---------------------------------------------------------------------------------------
# 2. Exact half-space stiffness (D-SIT-08) vs the Kausel-Roesset (1981) closed form
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("k", [0.0, 0.02, 0.05, 0.07, 0.1, 0.3])
def test_halfspace_stiffness_equals_kausel_roesset_closed_form(k):
    """K = 2kG[(1-s^2)/(2(1-rs)) [[r,1],[1,s]] - [[0,1],[1,0]]], r = sqrt(1-(kp/k)^2),
    s = sqrt(1-(ks/k)^2) in Kausel variables {u_x, i u_z}; k = 0 is the dashpot i w rho V*."""
    rho, vs, vp, b = 2.1, 400.0, 800.0, 0.03
    G, M = rho * vs * vs * cfactor(b), rho * vp * vp * cfactor(b)
    om = 2 * np.pi * 5.0
    Ko, Ki = ff.halfspace_stiffness(rho, G, M, om, k)
    if k == 0.0:
        ref = 1j * om * np.diag([np.sqrt(rho * G), np.sqrt(rho * M)])
    else:
        kp, ks = om / np.sqrt(M / rho), om / np.sqrt(G / rho)
        r, s = np.sqrt(1 - (kp / k) ** 2 + 0j), np.sqrt(1 - (ks / k) ** 2 + 0j)
        ref = 2 * k * G * ((1 - s * s) / (2 * (1 - r * s)) * np.array([[r, 1], [1, s]])
                           - np.array([[0, 1], [1, 0]]))
    np.testing.assert_allclose(Ko, ref, rtol=1e-10, atol=1e-10 * np.abs(ref).max())
    # incident (up-going) stiffness: the same with nu -> -nu; vertical limit -i w rho V*
    if k == 0.0:
        np.testing.assert_allclose(Ki, -ref, rtol=1e-12)


def test_transmitting_boundaries_are_complex_symmetric():
    """Reciprocity of the exterior regions: the cylindrical boundary R_mu (R1 §3.3 iii), the
    Waas P-SV and SH boundaries (R1 §2.6) and the exact half-space stiffness are complex
    symmetric (not Hermitian)."""
    col = _test_column()
    om = 2 * np.pi * 5.0
    md = tlm.column_modes(col, om)
    for mu in (0, 1):
        R, _ = axisym.boundary_stiffness(col, md, mu, 0.9)
        assert np.linalg.norm(R - R.T) <= 1e-10 * np.linalg.norm(R)
        assert np.linalg.norm(R - R.conj().T) > 1e-3 * np.linalg.norm(R)
    for side in (1, -1):
        R, _ = strip2d.boundary_psv(col, md, side)
        assert np.linalg.norm(R - R.T) <= 1e-10 * np.linalg.norm(R)
    R, _ = strip2d.boundary_sh(col, md)
    assert np.linalg.norm(R - R.T) <= 1e-10 * np.linalg.norm(R)
    for k in (0.01, 0.06, 0.2):
        Ko, Ki = ff.halfspace_stiffness(2.3, col.G[-1], col.M[-1], om, k)
        assert abs(Ko[0, 1] - Ko[1, 0]) <= 1e-12 * np.abs(Ko).max()
        assert abs(Ki[0, 1] - Ki[1, 0]) <= 1e-12 * np.abs(Ki).max()


# ---------------------------------------------------------------------------------------
# 3. POINT core elements: rigid-body null spaces and mass (R1 §3.3 i, §3.5)
# ---------------------------------------------------------------------------------------
def _free_column():
    """Column with a 'dashpot' base of zero strength: all interfaces free, no base terms."""
    h = np.array([0.5, 0.7, 1.0, 1.3])
    rho = np.array([1.9, 2.0, 2.1, 2.2])
    G, M = tlm.layer_moduli(rho, [150, 200, 250, 300], [300, 420, 500, 640], [0.05] * 4, [0.05] * 4)
    return tlm.Column(h=h, rho=rho, G=G, M=M, base=tlm.BASE_DASHPOT, cs=0j, cp=0j)


@pytest.mark.parametrize("mu", [0, 1])
def test_axisym_core_rigid_body_null_space_and_mass(mu):
    col = _free_column()
    ni, R0 = col.n_iface, 0.8
    z = -col.depth                                    # z up
    K = axisym.core_dynamic_stiffness(col, mu, R0, 0.0)

    def vec(fun):
        v = np.zeros(6 * ni, complex)
        for a, r in ((0, 0.0), (1, R0)):
            comps = fun(r)
            for c in range(3):
                v[(a * 3 + c) * ni:(a * 3 + c + 1) * ni] = comps[c]
        return v

    one, zero = np.ones(ni), np.zeros(ni)
    if mu == 1:
        modes = [vec(lambda r: (one, one, zero)),               # translation along x
                 vec(lambda r: (z, z, -r * one))]               # rotation about y: u = (z, 0, -x)
        mass_ref = 2 * np.sum(col.rho * col.h) * R0 ** 2 / 2
    else:
        modes = [vec(lambda r: (zero, zero, one)),              # translation along z
                 vec(lambda r: (zero, r * one, zero))]          # rigid torsion u_theta = rho
        mass_ref = np.sum(col.rho * col.h) * R0 ** 2 / 2
    for v in modes:
        assert np.abs(K @ v).max() <= 1e-12 * np.abs(K).max() * np.abs(v).max()
    assert np.abs(K - K.T).max() <= 1e-14 * np.abs(K).max()
    w = 3.0
    Mc = (K - axisym.core_dynamic_stiffness(col, mu, R0, w)) / w ** 2
    v = modes[0]
    assert (v @ Mc @ v).real == pytest.approx(mass_ref, rel=1e-10)


def test_strip_core_rigid_body_null_space():
    col = _free_column()
    ni, R0 = col.n_iface, 0.8
    z = -col.depth
    K = strip2d.core_psv(col, R0, 0.0)
    v = np.zeros(K.shape[0], complex)
    v[:3 * ni] = 1.0                                       # translation along x
    assert np.abs(K @ v).max() <= 1e-12 * np.abs(K).max()
    v = np.zeros(K.shape[0], complex)
    for ix, x in enumerate((-R0, 0.0, R0)):                # rotation about y
        v[ix * ni:(ix + 1) * ni] = z
        v[3 * ni + ix * ni:3 * ni + (ix + 1) * ni] = -x
    assert np.abs(K @ v).max() <= 1e-12 * np.abs(K).max() * np.abs(v).max()
    Ks = strip2d.core_sh(col, R0, 0.0)
    assert np.abs(Ks @ np.ones(Ks.shape[0])).max() <= 1e-12 * np.abs(Ks).max()


# ---------------------------------------------------------------------------------------
# 4. Free field (SITE Mode 2) vs an independent potential propagator (R1 §2.7a,b)
# ---------------------------------------------------------------------------------------
def _vnum(kb2, k):
    nu = np.sqrt(complex(kb2 - k * k))
    if nu.imag > 0:
        nu = -nu
    if abs(nu.imag) < 1e-14 * abs(nu) and nu.real < 0:
        nu = -nu
    return nu


def _E(rho, G, M, om, k, z):
    """State vector [u_x, u_z, s_xz, s_zz] (z up) of P-up, P-down, S-up, S-down potentials."""
    nup, nus = _vnum(rho * om ** 2 / M, k), _vnum(rho * om ** 2 / G, k)
    a = rho * om ** 2 / G - 2 * k * k
    cols = [np.array([-1j * k, 1j * s * nup, 2 * G * k * s * nup, -G * a]) * np.exp(1j * s * nup * z)
            for s in (-1, 1)]
    cols += [np.array([-1j * s * nus, -1j * k, G * a, 2 * G * k * s * nus]) * np.exp(1j * s * nus * z)
             for s in (-1, 1)]
    return np.array(cols).T


def _propagator(layers, hs, om, k, wave):
    """Displacements (u_x, u_z) at the layer interfaces (surface first) for a unit incident P or
    S wave from the half-space; free surface; exact continuous layers."""
    Ehs = _E(*hs, om, k, 0.0)

    def states(c):
        st = Ehs @ c
        out = [st]
        for (h, r, G, M) in reversed(layers):
            st = _E(r, G, M, om, k, h) @ np.linalg.solve(_E(r, G, M, om, k, 0.0), st)
            out.append(st)
        return out[::-1]

    inc = np.zeros(4, complex)
    inc[0 if wave == "P" else 2] = 1.0
    e1, e2 = np.eye(4)[1].astype(complex), np.eye(4)[3].astype(complex)
    s0, s1, s2 = states(inc)[0], states(e1)[0], states(e2)[0]
    x = np.linalg.solve(np.array([[s1[2], s2[2]], [s1[3], s2[3]]]), -np.array([s0[2], s0[3]]))
    return np.array([s[:2] for s in states(inc + x[0] * e1 + x[1] * e2)])


_PROF = [(6.0, 1.9, 180.0, 400.0, 0.04), (8.0, 2.0, 260.0, 560.0, 0.03), (10.0, 2.1, 350.0, 750.0, 0.02)]
_HS = (2.2, 600.0, 1200.0, 0.01)


@pytest.mark.parametrize("wave,wtype,angle", [("SV", 2, 0.0), ("SV", 2, 20.0), ("SV", 2, 45.0), ("P", 3, 30.0),
                                              ("P", 3, 0.0)])
def test_psv_free_field_matches_independent_propagator(tmp_path, wave, wtype, angle):
    """Three damped layers over a damped half-space, 5 Hz, 40 sublayers per layer (h <= lambda/14);
    FILE1 normalised at the top of the half-space (within motion) must agree with a continuous
    potential propagator solution at the four user interfaces of the profile (the within TF does
    not depend on the half-space model, so the vertical case is also covered).  SV at 45 deg is
    beyond the P critical angle (Vp/Vs = 2): evanescent reflected P in the half-space."""
    sub, f = 40, 5.0
    om = 2 * np.pi * f
    layers, props = [], []
    for (H, r, v, p, b) in _PROF:
        layers += [(H / sub, r, v, p, b, b)] * sub
        props.append((H, r, r * v * v * cfactor(b), r * p * p * cfactor(b)))
    rh, vh, ph, bh = _HS
    cm = 0 if wave == "SV" else 2
    rc, out = _site(tmp_path, layers, (0.0, rh, vh, ph, bh, bh), [int(f)], 1.0, nl=20,
                    waves=((wtype, 1, 1.0, 1.0, angle),), cl=sub * len(_PROF) + 1, cm=cm)
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    iu = [0] + list(np.cumsum([sub] * len(_PROF)))
    got = f1["U"][0, 0][iu][:, [0, 2]]
    k = om * np.sin(np.radians(angle)) / (vh if wave == "SV" else ph)
    assert f1["k"][0, 0] == pytest.approx(k, abs=1e-15)
    ref = _propagator(props, (rh, rh * vh * vh * cfactor(bh), rh * ph * ph * cfactor(bh)), om, k,
                      "S" if wave == "SV" else "P")
    ref = ref / ref[-1, 0 if wave == "SV" else 1]
    np.testing.assert_allclose(got, ref, atol=2e-3 * np.abs(ref).max())


def test_inclined_sv_mirror_symmetry(tmp_path):
    """Angle theta and 360 - theta are mirror images in x': k -> -k, U_x' unchanged and U_z'
    reversed after normalisation to the x' control motion."""
    layers = [(1.0, 1.9, 200.0, 400.0, 0.05, 0.05)] * 12
    hs = (0.0, 2.2, 700.0, 1400.0, 0.02, 0.02)
    res = {}
    for ang in (25.0, 335.0):
        wd = tmp_path / f"a{int(ang)}"
        rc, out = _site(wd, layers, hs, [3, 8], 1.0, waves=((2, 1, 1.0, 1.0, ang),))
        assert rc == 0, out
        res[ang] = read_container(wd / "FILE1")
    a, b = res[25.0], res[335.0]
    np.testing.assert_allclose(b["k"], -a["k"], rtol=1e-12)
    np.testing.assert_allclose(b["U"][..., 0], a["U"][..., 0], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(b["U"][..., 2], -a["U"][..., 2], rtol=1e-9, atol=1e-12)


def test_rayleigh_wave_surface_ellipticity_and_retrograde_motion(tmp_path):
    """Homogeneous elastic half-space, nu = 1/4: the Rayleigh wave exp(i(wt - kx)) has
    u_z/u_x = -1.4679 i at the surface (retrograde ellipse, z up) and u_x changes sign at
    0.1925 lambda_R depth (closed form from the Rayleigh potentials)."""
    vs, nu = 100.0, 0.25
    vp = _vp(vs, nu)
    rc, out = _site(tmp_path, [(0.125, 2.0, vs, vp, 0.0, 0.0)] * 240, (0.0, 2.0, vs, vp, 0.0, 0.0), [2], 10.0,
                    nl=20, waves=((1, 1, 1.0, 1.0, 0.0),), hslaw="uniform")
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    U = f1["U"][0, 0]
    k = f1["k"][0, 0]
    assert 2 * np.pi * 20.0 / k.real / vs == pytest.approx(0.919402, rel=3e-3)
    ratio = U[0, 2] / U[0, 0]
    assert abs(ratio.real) < 1e-6 and ratio.imag == pytest.approx(-1.46789, rel=5e-3)
    z = f1["depth_user"]
    ux = U[:, 0].real
    i = np.flatnonzero(np.sign(ux[:-1]) != np.sign(ux[1:]))[0]
    z0 = z[i] - ux[i] * (z[i + 1] - z[i]) / (ux[i + 1] - ux[i])
    assert z0 / (2 * np.pi / k.real) == pytest.approx(0.19251, rel=1e-2)


def test_within_tf_independent_of_halfspace_model(tmp_path):
    """Vertical SV with the control point at the top of the half-space: the field in the user
    layers is fixed by the motion there, so nl = 0 (rigid), 10 and 20 (dashpot buffer) must give
    identical FILE1 amplitudes (R1 §2.7a [V])."""
    layers = [(1.5, 1.9, 180.0, 360.0, 0.05, 0.05)] * 8 + [(2.0, 2.0, 300.0, 600.0, 0.03, 0.03)] * 5
    hs = (0.0, 2.2, 800.0, 1600.0, 0.02, 0.02)
    out = []
    for nl in (0, 10, 20):
        rc, txt = _site(tmp_path / f"nl{nl}", layers, hs, [2, 7, 15], 0.5, nl=nl, cl=14)
        assert rc == 0, txt
        out.append(read_container(tmp_path / f"nl{nl}" / "FILE1")["U"])
    np.testing.assert_allclose(out[1], out[0], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(out[2], out[0], rtol=1e-9, atol=1e-12)


def test_mixed_wave_field_unit_control_motion_through_api(tmp_path):
    """Rayleigh + inclined SV with ratios interpolated between frequency numbers 2 and 12: at the
    control point (xc, yc) on interface cl the horizontal motion returned by the binding API is
    exactly R_z(ang) e_x' (unit control motion in x') at every frequency, and the ratios follow the
    linear rule (constant outside)."""
    layers = [(1.0, 1.9, 200.0, 400.0, 0.05, 0.05)] * 15
    hs = (0.0, 2.2, 700.0, 1400.0, 0.02, 0.02)
    rc, out = _site(tmp_path, layers, hs, [1, 4, 8, 12, 20], 1.0, cl=3, cm=0,
                    waves=((1, 1, 0.3, 0.7, 0.0), (2, 1, 0.7, 0.3, 35.0)), freq1=2, freq2=12)
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    t = np.clip((np.array([1, 4, 8, 12, 20]) - 2) / 10.0, 0, 1)
    np.testing.assert_allclose(f1["ratio"][0], 0.3 + 0.4 * t, rtol=1e-12)
    np.testing.assert_allclose(f1["ratio"].sum(0), 1.0, rtol=1e-12)
    ang, xc, yc = 40.0, 3.0, -2.0
    a = np.radians(ang)
    for q in range(5):
        u = ff.free_field_at_nodes(f1, q, [[xc, yc, -2.0]], [3], ang, xc, yc)[0]
        np.testing.assert_allclose(u[:2], [np.cos(a), np.sin(a)], atol=1e-12)   # x' control motion = 1
        # one wavelength along x' later the Rayleigh part has rotated by exp(-i k L)
        kR = f1["k"][0, q]
        d = 5.0
        u2 = ff.free_field_at_nodes(f1, q, [[xc + d * np.cos(a), yc + d * np.sin(a), 0.0]], [3], ang, xc, yc)[0]
        Up = f1["ratio"][0, q] * f1["U"][0, q, 2] * np.exp(-1j * kR * d) + \
            f1["ratio"][1, q] * f1["U"][1, q, 2] * np.exp(-1j * f1["k"][1, q] * d)
        np.testing.assert_allclose(u2, [np.cos(a) * Up[0], np.sin(a) * Up[0], Up[2]], rtol=1e-12, atol=1e-14)


# ---------------------------------------------------------------------------------------
# 5. Independent confirmation of the two failing VPs (gate numbers are not artefacts)
# ---------------------------------------------------------------------------------------
def _col1d(h, G, rho, om, c):
    n = len(h) + 1
    K = np.zeros((n, n), complex)
    for j in range(len(h)):
        K[j:j + 2, j:j + 2] += (G[j] / h[j] * np.array([[1, -1], [-1, 1]])
                                - om ** 2 * rho[j] * h[j] * np.array([[5 / 12, 1 / 12], [1 / 12, 5 / 12]]))
    K[-1, -1] += 1j * om * c
    f = np.zeros(n, complex)
    f[-1] = 2j * om * c
    return np.linalg.solve(K, f)


def test_vp02b_geometric_gate_reproduced_independently(tmp_path):
    """VP-02b site with the deck-default geometric law (nl = 20): SITE's outcrop TF equals an
    independent 1-D FE column with my own geometric sublayer generator to 1e-8, and its error vs the
    closed form exceeds 1 % at 8.33 Hz -- the D-SIT-02 gate outcome is genuine."""
    H, vs, rho, xi, vr, rr, xr = 30.0, 200.0, 2.0, 0.05, 1000.0, 2.2, 0.01
    df = 1.0 / 6.0
    fn = [3, 10, 50]                                      # 0.5, 1.667, 8.333 Hz
    rc, out = _site(tmp_path, [(1.2, rho, vs, 2 * vs, xi, xi)] * 25, (0.0, rr, vr, 2 * vr, xr, xr), fn, df,
                    nl=20, cl=26, hslaw="geometric")
    assert rc == 0, out
    f1 = read_container(tmp_path / "FILE1")
    site_tf = np.abs(f1["U"][0, :, 0, 0] * f1["x_ucp"][0] / f1["x_outcrop"][0])
    Gl, Gr = rho * vs ** 2 * cfactor(xi), rr * vr ** 2 * cfactor(xr)
    c = np.sqrt(rr * Gr)
    for q, n in enumerate(fn):
        f = n * df
        om, lam = 2 * np.pi * f, vr / f
        Hh = 1.5 * lam
        h1 = min(1.2 * vr / vs, Hh / 20)
        i = np.arange(20)
        qq = brentq(lambda x: h1 * np.sum(x ** i) - Hh, 1.0 + 1e-14, 10.0) if h1 < Hh / 20 else 1.0
        hg = h1 * qq ** i
        hg *= Hh / hg.sum()
        if hg.max() > lam / 8 * (1 + 1e-12):
            hg = np.full(20, Hh / 20)
        hh = np.r_[np.full(25, 1.2), hg]
        u = _col1d(hh, np.r_[np.full(25, Gl), np.full(20, Gr)], np.r_[np.full(25, rho), np.full(20, rr)], om, c)
        uo = _col1d(hg, np.full(20, Gr), np.full(20, rr), om, c)[0]
        assert site_tf[q] == pytest.approx(abs(u[0] / uo), rel=1e-8)
        np.testing.assert_allclose(read_container(tmp_path / "FILE2")["h_gen"][q], hg, rtol=1e-10)
    vss, vrs = vs * np.sqrt(cfactor(xi)), vr * np.sqrt(cfactor(xr))
    om = 2 * np.pi * 50 * df
    ex = abs(1 / (np.cos(om / vss * H) + 1j * rho * vss / (rr * vrs) * np.sin(om / vss * H)))
    assert abs(site_tf[2] - ex) / ex > 0.01


def _cr_companion(nu, n_per_lambda, depth_lambdas=3):
    """Fundamental c_R/Vs of a uniform stratum on a rigid base from the Waas quadratic pencil,
    solved by a companion linearisation (independent of tlm.rayleigh_modes)."""
    lp = 2 * (1 - nu) / (1 - 2 * nu)
    om = 2 * np.pi
    h = 1.0 / n_per_lambda
    n_sub = int(round(depth_lambdas * n_per_lambda))
    G, M = tlm.layer_moduli(np.ones(n_sub), np.ones(n_sub), np.full(n_sub, np.sqrt(lp)), np.zeros(n_sub),
                            np.zeros(n_sub))
    col = tlm.Column(h=np.full(n_sub, h), rho=np.ones(n_sub), G=G, M=M)
    K0 = _waas_inplane(col, om, 0.0)
    K1 = (_waas_inplane(col, om, 1.0) - _waas_inplane(col, om, -1.0)) / 2.0      # i B_W
    K2 = (_waas_inplane(col, om, 1.0) + _waas_inplane(col, om, -1.0)) / 2.0 - K0  # A
    m = K0.shape[0]
    Z, I = np.zeros((m, m)), np.eye(m)
    k = sla.eig(np.block([[Z, I], [-K0, -K1]]), np.block([[I, Z], [Z, K2]]), right=False)
    k = k[(k.real > 0) & (np.abs(k.imag) < 1e-8 * np.abs(k))]
    return om / k.real.max()


def test_vp05_locking_reproduced_independently(tmp_path):
    """VP-05 at nu = 0.45, 10 sublayers per shear wavelength: SITE's FILE2 c_R equals an independent
    companion-form solution of the same discrete pencil, and both are +1.24 % off -- the VP-05
    failure is a property of the linear thin layers (D-SIT-01), not a coding error."""
    nu, ref = 0.45, 0.948960
    vs = 100.0
    vp = _vp(vs, nu)
    rc, out = _site(tmp_path, [(1.0, 2.0, vs, vp, 0.0, 0.0)] * 30, (0.0, 2.0, vs, vp, 0.0, 0.0), [1], 10.0,
                    nl=0, mode2=False)
    assert rc == 0, out
    f2 = read_container(tmp_path / "FILE2")
    k = f2["kR"][0]
    c_site = 2 * np.pi * 10.0 / k[tlm.select_mode(k, 1)].real / vs
    c_ind = _cr_companion(nu, 10)
    assert c_site == pytest.approx(c_ind, rel=1e-9)
    assert (c_site - ref) / ref == pytest.approx(0.01241, abs=2e-4)


# ---------------------------------------------------------------------------------------
# 6. Flexibility matrix: passivity and covariance (R1 §4.1, D-ANL-02)
# ---------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def point_files(tmp_path_factory):
    out = {}
    for dim in (2, 1):
        wd = tmp_path_factory.mktemp(f"pf{dim}")
        layers = [(0.5, 1.9, 150.0, 300.0, 0.05, 0.05)] * 6 + [(1.0, 1.9, 250.0, 500.0, 0.05, 0.05)] * 6
        rc, txt = _site(wd, layers, (0.0, 2.1, 400.0, 800.0, 0.03, 0.03), [1, 4], 2.0, nl=20, mode2=False,
                        hslaw="uniform")
        assert rc == 0, txt
        rc, txt = _point(wd, 2, 0.9, 2.0, dim=dim)
        assert rc == 0, txt
        out[dim] = (wd, read_container(wd / "FILE3"))
    return out


def _grid(dim):
    xs = np.arange(-3, 4) * 1.0
    X, Y = np.meshgrid(xs, xs if dim == 2 else [0.0])
    xy = np.c_[X.ravel(), Y.ravel()]
    xy = np.vstack([xy, xy[:5] + [0.3, 0.0]])
    iface = np.r_[np.ones(len(xy) - 5, int), [2, 2, 3, 3, 3]]
    return xy, iface


@pytest.mark.parametrize("dim", [2, 1])
def test_flexibility_is_passive(point_files, dim):
    """Under exp(+i w t) the power absorbed by the soil for any load vector P is
    -w/2 Im(P^H F P) >= 0, so Im(F_ff) (F symmetric) must be negative semi-definite."""
    _, f3 = point_files[dim]
    xy, iface = _grid(dim)
    for q in range(2):
        F = flexibility_matrix(f3, q, xy, iface)
        ev = np.linalg.eigvalsh(F.imag)
        assert ev.max() < 0.0, (q, ev.max(), ev.min())


def test_flexibility_mirror_translation_permutation_covariance(point_files):
    _, f3 = point_files[2]
    xy, iface = _grid(2)
    rng = np.random.default_rng(3)
    xy = xy + rng.uniform(-0.2, 0.2, xy.shape)
    n = len(xy)
    F = flexibility_matrix(f3, 1, xy, iface, symmetrize=False)
    # translation
    np.testing.assert_allclose(flexibility_matrix(f3, 1, xy + [17.3, -4.1], iface, symmetrize=False), F,
                               rtol=1e-9, atol=1e-12 * np.abs(F).max())
    # mirror y -> -y: F' = D F D, D = diag(1, -1, 1) per node
    D = np.tile([1.0, -1.0, 1.0], n)
    Fm = flexibility_matrix(f3, 1, xy * [1.0, -1.0], iface, symmetrize=False)
    np.testing.assert_allclose(Fm, D[:, None] * F * D[None, :], rtol=1e-9, atol=1e-12 * np.abs(F).max())
    # permutation of the nodes
    p = rng.permutation(n)
    Fp = flexibility_matrix(f3, 1, xy[p], iface[p], symmetrize=False)
    idx = (3 * p[:, None] + np.arange(3)[None, :]).ravel()
    np.testing.assert_allclose(Fp, F[np.ix_(idx, idx)], rtol=1e-9, atol=1e-12 * np.abs(F).max())


def test_point2_flexibility_reciprocity(point_files):
    """POINT2 (P1): odd coupling terms x sgn(x_i - x_j) must make the unsymmetrised plane-strain
    flexibility reciprocal: u_z(i) / P_x(j) = u_x(j) / P_z(i) across interfaces."""
    _, f3 = point_files[1]
    xy = np.array([[0.0, 0.0], [2.3, 0.0], [-3.7, 0.0], [5.1, 0.0], [1.0, 0.0]])
    iface = np.array([1, 2, 3, 1, 3])
    F = flexibility_matrix(f3, 1, xy, iface, symmetrize=False)
    assert np.linalg.norm(F - F.T) <= 2e-3 * np.linalg.norm(F)
    # sign: ahead of a horizontal surface push the surface moves down (low-frequency Lamb sense)
    Fs = flexibility_matrix(f3, 0, np.array([[0.0, 0.0], [3.0, 0.0], [-3.0, 0.0]]), [1, 1, 1], symmetrize=False)
    assert Fs[3 * 1 + 2, 0] == pytest.approx(-Fs[3 * 2 + 2, 0], rel=1e-12)
    assert Fs[3 * 1 + 2, 0].real < 0.0 and Fs[3 * 1, 2].real > 0.0


def test_point3_frequency_subset_consistency(point_files, tmp_path):
    """A FILE3 written for one frequency number equals the matching row of a FILE3 written for all
    frequencies (ANALYS frequency survey looks rows up by number)."""
    wd, f3_all = point_files[2]
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "FILE2").write_bytes((wd / "FILE2").read_bytes())
    rc, txt = _point(sub, 2, 0.9, 2.0, fnum=[4])
    assert rc == 0, txt
    f3_one = read_container(sub / "FILE3")
    xy, iface = _grid(2)
    F1 = flexibility_matrix(f3_one, frequency_row(f3_one, 4), xy, iface)
    F2 = flexibility_matrix(f3_all, frequency_row(f3_all, 4), xy, iface)
    np.testing.assert_allclose(F1, F2, rtol=1e-12, atol=1e-15)


# ---------------------------------------------------------------------------------------
# 7. Input checks not covered elsewhere
# ---------------------------------------------------------------------------------------
def test_site_rejects_wave_from_above_and_bad_law(tmp_path):
    layers = [(2.0, 1.9, 200.0, 400.0, 0.05, 0.05)] * 5
    hs = (0.0, 2.2, 700.0, 1400.0, 0.02, 0.02)
    rc, out = _site(tmp_path / "a", layers, hs, [2], 1.0, waves=((2, 1, 1.0, 1.0, 120.0),))
    assert rc != 0 and "does not come from below" in out
    rc, out = _site(tmp_path / "b", layers, hs, [2], 1.0, hslaw="parabolic")
    assert rc != 0 and "half-space sublayer law" in out
    rc, out = _site(tmp_path / "c", [(2.0, 1.9, 200.0, 220.0, 0.05, 0.05)] * 5, hs, [2], 1.0)
    assert rc != 0 and "bulk modulus" in out


def test_point_rejects_frequency_step_mismatch_and_site_mode2_rejects_stale_file2(tmp_path):
    layers = [(2.0, 1.9, 200.0, 400.0, 0.05, 0.05)] * 5
    hs = (0.0, 2.2, 700.0, 1400.0, 0.02, 0.02)
    rc, out = _site(tmp_path, layers, hs, [2, 3], 1.0, mode2=False)
    assert rc == 0, out
    rc, out = _point(tmp_path, 1, 1.0, 1.5)
    assert rc != 0 and "frequency step" in out
    rc, out = _site(tmp_path, layers, hs, [2, 3], 1.0, nl=10, mode1=0)
    assert rc != 0 and "re-run SITE Mode 1" in out


# ---------------------------------------------------------------------------------------
# 8. Defects found by the review (documented as strict xfails)
# ---------------------------------------------------------------------------------------
# REVIEW DEFECT (fixed): D-SIT-09 'least decay' (opt 2) was decided by rounding noise when several
# propagating modes have Im k = 0 (undamped site on a rigid base); tlm.choose_mode now treats |Im k|
# within TIE_RTOL max|k| as equal and breaks ties by the largest Re k.
def test_least_decay_selection_is_stable_for_undamped_site():
    picks = []
    for eps in (0.0, 1e-13, 1e-12):
        n = 60
        rho = np.full(n, 2.0) * (1 + eps)
        G, M = tlm.layer_moduli(rho, np.full(n, 100.0), np.full(n, 200.0), np.zeros(n), np.zeros(n))
        md = tlm.column_modes(tlm.Column(h=np.full(n, 0.5), rho=rho, G=G, M=M), 2 * np.pi * 20.0)
        picks.append(md.kR[tlm.select_mode(md.kR, 2)].real)
    assert max(picks) - min(picks) <= 1e-6 * max(picks)


@pytest.mark.xfail(strict=True, reason="Known limitation (D-W1-03): the variable-depth half-space (1.5 lambda "
                   "buffer + dashpots) puts the first Love cut-off at about 12.40 Hz instead of 11.547 Hz (+7.4 %); "
                   "VP-06 checks it with the provisional 10 % tolerance, this test documents that 0.5 % is not met")
def test_vp06_first_love_cutoff_frequency(tmp_path):
    b1, b2, H, rho = 200.0, 400.0, 10.0, 2.0
    fnum = list(range(110, 136, 2))                    # 11.0 .. 13.4 Hz with df = 0.1
    rc, out = _site(tmp_path, [(H / 20, rho, b1, 2 * b1, 0.0, 0.0)] * 20, (0.0, rho, b2, 2 * b2, 0.0, 0.0),
                    fnum, 0.1, nl=20, mode2=False)
    assert rc == 0, out
    f2 = read_container(tmp_path / "FILE2")
    c1 = []
    for q, f in enumerate(f2["freq"]):
        k = f2["kL"][q]
        w = 2 * np.pi * f
        ok = (k.real > 0) & (np.abs(k.imag) <= 0.05 * k.real) & (w / k.real < 1.3 * b2)
        c = np.sort(w / k[ok].real)
        c1.append(c[1])
    c1 = np.asarray(c1)
    i = np.flatnonzero((c1[:-1] - b2) * (c1[1:] - b2) <= 0)[0]
    fc = f2["freq"][i] + (b2 - c1[i]) * (f2["freq"][i + 1] - f2["freq"][i]) / (c1[i + 1] - c1[i])
    assert fc == pytest.approx(11.547, rel=0.005)


# REVIEW DEFECT (fixed): the fixed 'propagating' threshold |Im k| <= 0.5 Re k rejected every mode once
# beta >= ~0.447 (material attenuation tan(delta/2) = beta/sqrt(1-beta^2) > 0.5), although EDU-04 accepts
# beta < 0.5; the propagating sector is now widened by the material loss angle (tlm.propagating_ratio).
def test_surface_wave_selection_for_heavily_damped_soil():
    n = 40
    beta = 0.46
    rho = np.full(n, 2.0)
    G, M = tlm.layer_moduli(rho, np.full(n, 100.0), np.full(n, 200.0), np.full(n, beta), np.full(n, beta))
    md = tlm.column_modes(tlm.Column(h=np.full(n, 0.25), rho=rho, G=G, M=M), 2 * np.pi * 20.0)
    j = tlm.select_mode(md.kR, 1)
    # the fundamental mode of a uniform stratum travels at about 0.92 Vs in elastic soil.  Material
    # damping (SASSI form) multiplies the wavenumbers by (sqrt(1 - beta^2) - i beta), so the phase
    # velocity w / Re k of the damped fundamental mode is 0.92 Vs / sqrt(1 - beta^2).
    # [implementer's correction: the review version compared w / Re k with 0.92 (rel 0.1), which the
    # fundamental mode cannot meet at beta = 0.46 (it gives 1.050 = 0.934 / 0.889, the elastic
    # fundamental of this column continued to complex frequency w / sqrt(c(beta)))]
    assert 2 * np.pi * 20.0 / md.kR[j].real / 100.0 == pytest.approx(0.92 / np.sqrt(1 - beta ** 2), rel=0.1)


# REVIEW DEFECT (fixed): SITE Mode 2 alone (_load_file2) only compared layer counts and nl with FILE2; a
# FILE2 built from different layer properties was used silently for vertical/surface waves while inclined
# waves used the deck properties.  _check_file2_matches now compares properties and generated sublayers.
def test_site_mode2_alone_rejects_file2_with_other_properties(tmp_path):
    hs = (0.0, 2.2, 700.0, 1400.0, 0.02, 0.02)
    rc, out = _site(tmp_path, [(2.0, 1.9, 200.0, 400.0, 0.05, 0.05)] * 5, hs, [2, 3], 1.0, mode2=False)
    assert rc == 0, out
    rc, out = _site(tmp_path, [(2.0, 1.9, 300.0, 600.0, 0.05, 0.05)] * 5, hs, [2, 3], 1.0, mode1=0)
    assert rc != 0
