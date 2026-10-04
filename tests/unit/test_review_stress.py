"""Adversarial review tests of the 'stress' work package (STRESS module, sassi.core.stress_lib).

The implementer's tests compare STRESS with references built from the same FILE4 recovery operators
(``S . U_e``).  These tests check the module from other angles:

* closed-form physics that does not use the FILE4 operators: complex moduli written out by hand
  (D-CNV-03/04: G* = rho Vs^2 c(beta_s), M* = rho Vp^2 c(beta_p), lam* = M* - 2G*, shells
  E* = E c(beta) with nu real), uniform-strain and constant-curvature patch fields, spring and beam
  end-force signs (D-STR-05), SHELL membrane outputs as stresses (D-STR-10, D-W1-14);
* the convolution against time-domain closed forms: a single harmonic for seismic input
  (``U_g = -g A / w^2``, e^{+iwt}, D-CNV-06, D-STR-02) and a band-limited load for vibration input;
* invariances: linearity in gravity and in the motion scaling, independence of the element batching,
  any SITE direction / angle in the rigid-body subtraction;
* the interpolation range rules (D-MOT-03 applied to STFs), the listing time of the maximum, the
  ELEMENT_CENTER 'ordered group #' (D-FIL-06), the Frames.txt example of the manual (spec 05d 8 item 9);
* error paths (Error 80 for GENERAL groups, data-check mode).

Tests marked ``xfail(strict=True)`` document defects found in the review; they turn into XPASS failures
once the defect is fixed, so the marker must then be removed.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Callable, Tuple

import numpy as np
import pytest

from sassi.conventions import cfactor, element_result_name
from sassi.core import stress_lib as SL
from sassi.io import decks, textfiles
from sassi.io.container import read_container
from sassi.modules.base import run_module
from sassi.verify.problems.vp_house import G_SI, add_element, add_node, house_deck, run_house
from sassi.verify.problems.vp_motion import synthetic_motion, write_file8, write_motion_file

NFFT, DT = 256, 0.02
DF = 1.0 / (NFFT * DT)                     # 0.1953125 Hz
FNUM = np.arange(2, 101)                   # SSI frequencies 0.39 .. 19.5 Hz


# =============================================================================================
# small independent helpers (deck writing, FILE8 synthesis, file reading)
# =============================================================================================
def _stress_deck(gravity: float = G_SI, **params) -> decks.Deck:
    d = decks.new("STRESS")
    d["model"] = "m"
    d["thfile"] = "eq.acc"
    d["nft"], d["delt"], d["df"] = NFFT, DT, DF
    d["gravity"] = gravity
    d["mult"], d["max"] = 1.0, 0.0
    for k, v in params.items():
        d[k] = v
    return d


def _eout(d: decks.Deck, group: int, elems, codes) -> None:
    codes = (list(codes) + [0] * 12)[:12]
    for e in elems:
        d.table("eout").append([group, e] + codes)


def _run(wd: Path, d: decks.Deck, strainout: bool = False) -> Tuple[int, str]:
    p = decks.write(Path(wd) / "m.str", d)
    if strainout:                                           # EDUOPT,STRAINOUT,1 (D-STR-04), see contract issue
        p.write_text(p.read_text().replace("[params]\n", "[params]\nstrainout = 1\n", 1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rc = run_module("STRESS", "m", wd)
    return rc, (Path(wd) / "m_STRESS.out").read_text(encoding="utf-8")


def _house(wd: Path, d: decks.Deck) -> None:
    rc, out = run_house(wd, d)
    assert rc == 0, out[-2000:]


def _file8(wd: Path, field: Callable, fnum=FNUM, **meta) -> None:
    """FILE8 whose TF at equation (node, dof) is ``field(node, xyz, f)[dof - 1]``."""
    f4 = read_container(Path(wd) / "m.N4", "FILE4")
    xyz = {int(n): np.asarray(p, float) for n, p in zip(f4["node_id"], f4["node_xyz"])}
    fnum = np.asarray(fnum, dtype=np.int64)
    H = np.zeros((len(fnum), len(f4["eq_node"])), dtype=complex)
    for q, n in enumerate(fnum):
        for i, (node, dof) in enumerate(zip(f4["eq_node"], f4["eq_dof"])):
            H[q, i] = field(int(node), xyz[int(node)], n * DF)[int(dof) - 1]
    meta.setdefault("model_hash", str(f4.meta.get("model_hash", "")))
    write_file8(Path(wd) / "FILE8", fnum, DF, f4["eq_node"], f4["eq_dof"], H, nfft=NFFT, delt=DT, **meta)


def _tfu(wd, et, g, e, comp, ext="TFU"):
    f, H, cplx = textfiles.read_tf(Path(wd) / element_result_name(et, g, e, comp, ext))
    assert cplx
    return f, H


def _ths(wd, et, g, e, comp):
    return textfiles.read_history(Path(wd) / element_result_name(et, g, e, comp, "THS"))


def _d(f):
    """An arbitrary smooth complex frequency dependence of the deformation (not used by STRESS)."""
    r2 = (np.asarray(f) / 6.0) ** 2
    return (0.3 + r2) / (1.0 - r2 + 0.1j)


def _harmonic(wd, k0=7, a0=0.3, kind="cos"):
    t = np.arange(NFFT) * DT
    w0 = 2 * np.pi * k0 * DF
    a = a0 * (np.cos(w0 * t) if kind == "cos" else np.sin(w0 * t))
    write_motion_file(Path(wd) / "eq.acc", a, DT)
    return t, w0, a


# =============================================================================================
# models
# =============================================================================================
SPRING_K, SPRING_DAMP = (1.0e5, 2.0e5, 3.0e5, 10.0, 20.0, 30.0), 0.03


def _spring_model(wd, damp=SPRING_DAMP):
    d = house_deck()
    add_node(d, 1, 0.0, 0.0, 0.0)
    add_node(d, 2, 1.0, 0.5, 0.0)
    d.table("springprops").append([1, *SPRING_K, damp])
    d.table("groups").append([1, 7, "spring"])
    add_element(d, 1, 1, [1, 2], prop=1)
    _house(wd, d)


def _solid_model(wd, mat=(1, 3, 900.0, 400.0, 20.0, 0.06, 0.03), n_elem=1):
    """A row of distorted hexahedra along X (group 1); type-3 material with beta_p != beta_s."""
    d = house_deck()
    nid = {}
    n = 0
    for ix in range(n_elem + 1):
        for iy, y in enumerate((0.0, 1.0)):
            for iz, z in enumerate((0.0, 1.0)):
                n += 1
                nid[(ix, iy, iz)] = n
                add_node(d, n, ix * 1.0 + 0.07 * y * z + 0.03 * z, y + 0.05 * z * (ix % 2), z + 0.04 * ix * y)
    d.table("materials").append(list(mat))
    d.table("groups").append([1, 1, "block"])
    for e in range(n_elem):
        ns = [nid[(e, 0, 0)], nid[(e + 1, 0, 0)], nid[(e + 1, 1, 0)], nid[(e, 1, 0)],
              nid[(e, 0, 1)], nid[(e + 1, 0, 1)], nid[(e + 1, 1, 1)], nid[(e, 1, 1)]]
        add_element(d, 1, e + 1, ns, etype=1, mat=1)
    _house(wd, d)


def _moduli_type3(vp, vs, weight, dp, ds):
    """Independent complex moduli (D-CNV-03/04): G* = rho Vs^2 c(ds), M* = rho Vp^2 c(dp), lam* = M* - 2 G*."""
    rho = weight / G_SI
    c = lambda b: 1.0 - 2.0 * b * b + 2j * b * np.sqrt(1.0 - b * b)
    Gs = rho * vs ** 2 * c(ds)
    Ms = rho * vp ** 2 * c(dp)
    return Gs, Ms, Ms - 2.0 * Gs


# =============================================================================================
# 1. SPRING: sign, complex stiffness and the seismic convolution in closed form
# =============================================================================================
def test_spring_force_sign_complex_stiffness_and_seismic_harmonic(tmp_path):
    """D-STR-05: F = k*(u_J - u_I) per global component with k* = k c(damp) (D-CNV-04); a harmonic
    control acceleration a0 cos(w0 t) (in g) gives u_g = -g a0 cos(w0 t)/w0^2 (D-CNV-06), so the force
    history is Re(STF u_hat e^{+i w0 t}) (e^{+iwt}, D-CNV-01).  FILE8 holds total-motion TFs (H = 1 + ...)."""
    _spring_model(tmp_path)
    c = np.array([0.37, -0.21, 0.05, 2e-3, -1e-3, 4e-3])
    _file8(tmp_path, lambda n, x, f: np.r_[1.0, 0, 0, 0, 0, 0] + (c * _d(f) if n == 2 else 0.0))
    t, w0, _ = _harmonic(tmp_path)
    d = _stress_deck(itran=1)
    _eout(d, 1, [1], [2] * 6)
    rc, out = _run(tmp_path, d)
    assert rc == 0, out[-2000:]
    k0 = 7
    for j, comp in enumerate(("FX", "FY", "FZ", "MXX", "MYY", "MZZ")):
        ks = SPRING_K[j] * cfactor(SPRING_DAMP)
        f, tfu = _tfu(tmp_path, "SPRING", 1, 1, comp)
        stf = ks * c[j] * _d(f)
        assert np.allclose(tfu, stf, rtol=1e-12, atol=1e-14 * np.max(np.abs(stf))), comp
        th, dt = _ths(tmp_path, "SPRING", 1, 1, comp)
        u_hat = -G_SI * 0.3 / w0 ** 2                      # u_g(t) = Re(u_hat e^{i w0 t})
        ref = np.real(ks * c[j] * _d(k0 * DF) * u_hat * np.exp(1j * w0 * t))
        assert dt == DT
        assert np.max(np.abs(th - ref)) <= 1e-11 * np.max(np.abs(ref)), comp


def test_spring_vibration_band_limited_load_is_a_scaled_copy(tmp_path):
    """Vibration (type 1): the STF is used directly (no rigid body) and held at STF_1 below f_1
    (D-MOT-03); a real constant STF s on [0, f_N] and a load band-limited below f_N give s * load(t)."""
    _spring_model(tmp_path, damp=0.0)
    c = 0.37
    _file8(tmp_path, lambda n, x, f: np.r_[(c if n == 2 else 0.0) - 0.1, 0, 0, 0, 0, 0], type=1)
    rng = np.random.default_rng(3)
    F = np.zeros(NFFT // 2 + 1, complex)
    F[0] = 0.5
    F[1:90] = rng.standard_normal(89) + 1j * rng.standard_normal(89)    # bins < f_N bin (100)
    load = np.fft.irfft(F, NFFT)
    write_motion_file(tmp_path / "eq.acc", load, DT)
    d = _stress_deck(type=1, itran=1)
    _eout(d, 1, [1], [2] * 6)
    rc, out = _run(tmp_path, d)
    assert rc == 0, out[-2000:]
    fx, _ = _ths(tmp_path, "SPRING", 1, 1, "FX")
    assert np.max(np.abs(fx - SPRING_K[0] * c * load)) <= 1e-12 * np.max(np.abs(SPRING_K[0] * c * load))
    # TFI: f = 0 row = STF_1 (vibration), above f_N zero, SSI frequencies exact
    fi, tfi = _tfu(tmp_path, "SPRING", 1, 1, "FX", "TFI")
    assert fi[0] == 0.0 and tfi[0] == pytest.approx(SPRING_K[0] * c, rel=1e-12)
    assert len(fi) == FNUM[-1] + 1


def test_seismic_tfi_range_rules(tmp_path):
    """D-MOT-03 applied to the STF (D-STR-01): seismic STF linear from 0 at f = 0 to STF_1, exact at the SSI
    frequencies, rows up to f_N only."""
    _spring_model(tmp_path)
    fnum = np.arange(5, 61, 5)
    _file8(tmp_path, lambda n, x, f: np.r_[1.0 + (0.2 * _d(f) if n == 2 else 0.0), 0, 0, 0, 0, 0], fnum=fnum)
    _harmonic(tmp_path)
    d = _stress_deck(itran=1, interopt=6)
    _eout(d, 1, [1], [1])
    assert _run(tmp_path, d)[0] == 0
    f, tfu = _tfu(tmp_path, "SPRING", 1, 1, "FX")
    fi, tfi = _tfu(tmp_path, "SPRING", 1, 1, "FX", "TFI")
    k = np.rint(fi / DF).astype(int)
    assert np.array_equal(k, np.arange(fnum[-1] + 1))
    assert tfi[0] == 0.0
    assert np.allclose(tfi[1:5], tfu[0] * np.arange(1, 5) / 5.0, rtol=1e-12)       # linear 0 -> STF_1
    assert np.allclose(tfi[fnum], tfu, rtol=1e-12)


# =============================================================================================
# 2. continuum elements: complex moduli and strains from closed-form fields
# =============================================================================================
def test_solid_uniaxial_strain_and_simple_shear_use_the_complex_moduli(tmp_path):
    """SOLID centroid stresses for linear fields on a distorted hexahedron (patch test):
    u_x = eps x  ->  SXX = M* eps, SYY = SZZ = lam* eps;   u_x = gam z  ->  SXZ = G* gam.
    Moduli written out independently (type 3: Vp, Vs; beta_p = 0.06 on M, beta_s = 0.03 on G)."""
    _solid_model(tmp_path)
    Gs, Ms, Ls = _moduli_type3(900.0, 400.0, 20.0, 0.06, 0.03)
    eps, gam = 2.0e-4, 3.0e-4
    _harmonic(tmp_path)
    comps = ("SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ")
    cases = {
        "uniaxial strain": (lambda n, x, f: np.r_[1.0 + _d(f) * eps * x[0], 0, 0, 0, 0, 0],
                            {"SXX": Ms * eps, "SYY": Ls * eps, "SZZ": Ls * eps}, {"EXX": eps}),
        "simple shear xz": (lambda n, x, f: np.r_[1.0 + _d(f) * gam * x[2], 0, 0, 0, 0, 0],
                            {"SXZ": Gs * gam}, {"EXZ": gam}),
    }
    for name, (field, sref, eref) in cases.items():
        _file8(tmp_path, field)
        d = _stress_deck(itran=1)
        _eout(d, 1, [1], [1] * 7)
        rc, out = _run(tmp_path, d, strainout=True)
        assert rc == 0, out[-2000:]
        scale = max(abs(v) for v in sref.values())
        for c in comps:
            f, h = _tfu(tmp_path, "SOLID", 1, 1, c)
            ref = sref.get(c, 0.0) * _d(f)
            assert np.max(np.abs(h - ref)) <= 1e-9 * scale * np.max(np.abs(_d(f))), (name, c)
        escale = max(abs(v) for v in eref.values())
        for c in ("EXX", "EYY", "EZZ", "EXY", "EXZ", "EYZ"):          # engineering strains (D-STR-04)
            f, h = _tfu(tmp_path, "SOLID", 1, 1, c)
            ref = eref.get(c, 0.0) * _d(f)
            assert np.max(np.abs(h - ref)) <= 1e-9 * escale * np.max(np.abs(_d(f))), (name, c)


def test_solid_eoct_of_simple_shear(tmp_path):
    """Octahedral shear strain (requirements 4.10 item 5, engineering form) of simple shear gamma:
    (2/3) sqrt(1.5) gamma = (sqrt 6/3) gamma at every time step; SOCT = (sqrt 6/3) |SXZ|."""
    _solid_model(tmp_path)
    gam = 3.0e-4
    _file8(tmp_path, lambda n, x, f: np.r_[1.0 + _d(f) * gam * x[2], 0, 0, 0, 0, 0])
    write_motion_file(tmp_path / "eq.acc", 0.2 * synthetic_motion(200, DT, seed=2, fmax=15.0), DT)
    d = _stress_deck()
    _eout(d, 1, [1], [0, 0, 0, 0, 2, 0, 2])
    rc, out = _run(tmp_path, d, strainout=True)
    assert rc == 0, out[-2000:]
    exz, _ = _ths(tmp_path, "SOLID", 1, 1, "EXZ")
    eoct, _ = _ths(tmp_path, "SOLID", 1, 1, "EOCT")
    sxz, _ = _ths(tmp_path, "SOLID", 1, 1, "SXZ")
    soct, _ = _ths(tmp_path, "SOLID", 1, 1, "SOCT")
    r6 = np.sqrt(6.0) / 3.0
    assert np.max(np.abs(eoct - r6 * np.abs(exz))) <= 1e-9 * np.max(np.abs(exz))
    assert np.max(np.abs(soct - r6 * np.abs(sxz))) <= 1e-9 * np.max(np.abs(sxz))


def test_plane_uniform_strain_plane_strain_moduli(tmp_path):
    """PLANE (X-Z, plane strain): u_x = a x + b z, u_z = c x + e z ->
    SXX = M* a + lam* e, SZZ = lam* a + M* e, TXZ = G* (b + c); strains EXX, EZZ, EXZ (engineering)."""
    d = house_deck(dim=1)
    pts = {1: (0.0, 0.0), 2: (1.2, 0.1), 3: (1.1, 1.0), 4: (-0.1, 0.9)}
    for n, (x, z) in pts.items():
        add_node(d, n, x, 0.0, z)
    d.table("materials").append([1, 3, 700.0, 300.0, 19.0, 0.07, 0.02])
    d.table("groups").append([1, 4, "plane"])
    add_element(d, 1, 1, [1, 2, 3, 4], etype=1, mat=1)
    _house(tmp_path, d)
    Gs, Ms, Ls = _moduli_type3(700.0, 300.0, 19.0, 0.07, 0.02)
    a, b, c, e = 1.0e-4, 2.0e-4, -0.5e-4, -3.0e-4
    _file8(tmp_path, lambda n, x, f: np.r_[1.0 + _d(f) * (a * x[0] + b * x[2]), 0, _d(f) * (c * x[0] + e * x[2]),
                                           0, 0, 0])
    _harmonic(tmp_path)
    dk = _stress_deck(itran=1)
    _eout(dk, 1, [1], [1, 1, 1])
    rc, out = _run(tmp_path, dk, strainout=True)
    assert rc == 0, out[-2000:]
    ref = {"SXX": Ms * a + Ls * e, "SZZ": Ls * a + Ms * e, "TXZ": Gs * (b + c), "EXX": a, "EZZ": e, "EXZ": b + c}
    for comp, v in ref.items():
        f, h = _tfu(tmp_path, "PLANE", 1, 1, comp)
        assert np.allclose(h, v * _d(f), rtol=1e-9, atol=1e-12 * abs(v)), comp


# =============================================================================================
# 3. SHELL: membrane outputs are stresses (D-STR-10 / D-W1-14), moments per unit length
# =============================================================================================
def _shell_model(wd, E=3.0e7, nu=0.2, beta=0.05, t=0.3, a=2.0, b=1.5):
    d = house_deck()
    for n, (x, y) in enumerate(((0, 0), (a, 0), (a, b), (0, b)), start=1):
        add_node(d, n, x, y, 0.0, fix=(0, 0, 0, 0, 0, 1))
    d.table("materials").append([1, 1, E, nu, 24.0, beta, beta])
    d.table("groups").append([1, 3, "plate"])
    add_element(d, 1, 1, [1, 2, 3, 4], mat=1, thick=t)
    _house(wd, d)


def test_shell_membrane_stresses_and_bending_moments_closed_form(tmp_path):
    """Rectangle in the X-Y plane (local x'y'z' = XYZ, spec 08 4.6).  Uniform membrane strain:
    FXX = E*/(1-nu^2)(ex + nu ey) etc. -- a *stress*, independent of the thickness (D-W1-14);
    constant curvature w = -(kx x^2 + ky y^2 + kt x y)/2 (theta_x = dw/dy, theta_y = -dw/dx):
    MXX = D (kx + nu ky), MYY = D (ky + nu kx), MXY = D (1-nu)/2 kt, D = E* t^3/(12 (1-nu^2))."""
    E, nu, beta, t = 3.0e7, 0.2, 0.05, 0.3
    _shell_model(tmp_path, E, nu, beta, t)
    Es = E * cfactor(beta)
    ex, ey, gxy = 1.0e-4, -0.4e-4, 2.5e-4
    kx, ky, kt = 2.0e-4, 1.0e-4, -3.0e-4

    def field(n, x, f):
        X, Y = x[0], x[1]
        w = -(kx * X * X + ky * Y * Y + kt * X * Y) / 2.0
        thx = -(ky * Y + kt * X / 2.0)                     # dw/dy
        thy = kx * X + kt * Y / 2.0                        # -dw/dx
        return np.r_[1.0 + _d(f) * (ex * X + gxy / 2 * Y), _d(f) * (gxy / 2 * X + ey * Y), _d(f) * w,
                     _d(f) * thx, _d(f) * thy, 0.0]

    _file8(tmp_path, field)
    _harmonic(tmp_path)
    dk = _stress_deck(itran=1)
    _eout(dk, 1, [1], [1] * 6)
    rc, out = _run(tmp_path, dk)
    assert rc == 0, out[-2000:]
    Dm = Es / (1 - nu ** 2)
    Db = Es * t ** 3 / (12 * (1 - nu ** 2))
    ref = {"FXX": Dm * (ex + nu * ey), "FYY": Dm * (ey + nu * ex), "FXY": Dm * (1 - nu) / 2 * gxy,
           "MXX": Db * (kx + nu * ky), "MYY": Db * (ky + nu * kx), "MXY": Db * (1 - nu) / 2 * kt}
    for comp, v in ref.items():
        f, h = _tfu(tmp_path, "SHELL", 1, 1, comp)
        assert np.allclose(h, v * _d(f), rtol=1e-8, atol=1e-10 * abs(v)), comp


# =============================================================================================
# 4. BEAMS: end forces exerted on the element (D-STR-05)
# =============================================================================================
def test_beam_axial_and_torsion_end_force_signs(tmp_path):
    """Beam along +X, K on +Y: extension delta gives FXI = -EA delta/L, FXJ = +EA delta/L (forces exerted on
    the element: tension pulls J forward, I backward); a twist phi of J gives MXI = -GJ phi/L, MXJ = +GJ phi/L.
    E* = E c(beta), G* = E/(2(1+nu)) c(beta) (beta_p = beta_s, D-CNV-04)."""
    E, nu, beta, L, A, J = 2.0e8, 0.3, 0.04, 3.0, 0.25, 0.012
    d = house_deck()
    add_node(d, 1, 0.0, 0.0, 0.0)
    add_node(d, 2, L, 0.0, 0.0)
    add_node(d, 3, 0.0, 1.0, 0.0, fix=(1,) * 6)
    d.table("materials").append([1, 1, E, nu, 78.0, beta, beta])
    d.table("beamprops").append([1, A, 0.2, 0.18, J, 0.004, 0.006])
    d.table("groups").append([1, 2, "member"])
    add_element(d, 1, 1, [1, 2, 3], mat=1)
    _house(tmp_path, d)
    delta, phi = 1.0e-4, 2.0e-4
    _file8(tmp_path, lambda n, x, f: np.r_[1.0 + (_d(f) * delta if n == 2 else 0.0), 0, 0,
                                           (_d(f) * phi if n == 2 else 0.0), 0, 0])
    _harmonic(tmp_path)
    dk = _stress_deck(itran=1)
    _eout(dk, 1, [1], [1] * 12)
    rc, out = _run(tmp_path, dk)
    assert rc == 0, out[-2000:]
    Es, Gs = E * cfactor(beta), E / (2 * (1 + nu)) * cfactor(beta)
    ref = {"FXI": -Es * A * delta / L, "FXJ": Es * A * delta / L, "MXI": -Gs * J * phi / L, "MXJ": Gs * J * phi / L,
           "FYI": 0.0, "FZI": 0.0, "MYI": 0.0, "MZI": 0.0, "FYJ": 0.0, "FZJ": 0.0, "MYJ": 0.0, "MZJ": 0.0}
    scale = Es * A * delta / L
    for comp, v in ref.items():
        f, h = _tfu(tmp_path, "BEAMS", 1, 1, comp)
        assert np.max(np.abs(h - v * _d(f))) <= 1e-9 * abs(scale) * np.max(np.abs(_d(f))), comp


# =============================================================================================
# 5. invariances
# =============================================================================================
def test_histories_are_linear_in_gravity_and_motion_scaling(tmp_path):
    """U_g = -g A/w^2: doubling the deck gravity doubles the seismic stresses; mult = 2.5 scales them by 2.5;
    scaling to a peak 'max' equals mult = max / max|a|."""
    _spring_model(tmp_path)
    _file8(tmp_path, lambda n, x, f: np.r_[1.0 + (0.3 * _d(f) if n == 2 else 0.0), 0, 0, 0, 0, 0])
    acc = 0.2 * synthetic_motion(200, DT, seed=6, fmax=15.0)
    write_motion_file(tmp_path / "eq.acc", acc, DT)

    def fx(**kw):
        d = _stress_deck(**kw)
        _eout(d, 1, [1], [2])
        rc, out = _run(tmp_path, d)
        assert rc == 0, out[-1500:]
        return _ths(tmp_path, "SPRING", 1, 1, "FX")[0]

    base = fx()
    assert np.allclose(fx(gravity=2 * G_SI), 2 * base, rtol=1e-12, atol=1e-14 * np.max(np.abs(base)))
    assert np.allclose(fx(mult=2.5), 2.5 * base, rtol=1e-12, atol=1e-14 * np.max(np.abs(base)))
    pk = np.max(np.abs(acc))
    assert np.allclose(fx(mult=0.0, max=0.5), (0.5 / pk) * base, rtol=1e-10, atol=1e-13 * np.max(np.abs(base)))


def test_results_do_not_depend_on_element_batching(tmp_path, monkeypatch):
    """Elements are processed in column batches (BATCH_COLS); one element per batch must give bit-identical
    histories to the default batching (the reshapes (nF, b, nc) <-> (nF, b*nc) are consistent)."""
    from sassi.modules import stress
    _solid_model(tmp_path, n_elem=5)
    _file8(tmp_path, lambda n, x, f: np.r_[1.0 + _d(f) * (1e-4 * x[0] * x[2] + 2e-4 * x[1]),
                                           _d(f) * 1e-4 * x[0] ** 2, _d(f) * 3e-4 * x[0] * x[1], 0, 0, 0])
    write_motion_file(tmp_path / "eq.acc", 0.2 * synthetic_motion(200, DT, seed=8, fmax=15.0), DT)
    d = _stress_deck(savemax=1)
    _eout(d, 1, range(1, 6), [2] * 7)
    assert _run(tmp_path, d)[0] == 0
    ref = {(e, c): _ths(tmp_path, "SOLID", 1, e, c)[0] for e in range(1, 6) for c in ("SXX", "SYZ", "SOCT")}
    emax = (tmp_path / "ELEMENT_CENTER_ABS_MAX_STRESSES.TXT").read_text()
    monkeypatch.setattr(stress, "BATCH_COLS", 1)
    assert _run(tmp_path, d)[0] == 0
    for (e, c), v in ref.items():
        assert np.array_equal(_ths(tmp_path, "SOLID", 1, e, c)[0], v), (e, c)
    assert (tmp_path / "ELEMENT_CENTER_ABS_MAX_STRESSES.TXT").read_text() == emax


@pytest.mark.parametrize("cm,ang", [(1, 30.0), (2, 0.0), (0, 125.0)])
def test_rigid_body_subtraction_is_exact_for_any_control_direction(tmp_path, cm, ang):
    """The module subtracts the rigid body of the FILE8 control direction (cm, ang) before the recovery;
    since S annihilates translations the STF must equal the plain product S . U_e of the FILE4 operators
    for any direction (compared here with an explicit loop over the element DOFs)."""
    _solid_model(tmp_path, n_elem=2)
    a = np.deg2rad(ang)
    e = {0: (np.cos(a), np.sin(a), 0.0), 1: (-np.sin(a), np.cos(a), 0.0), 2: (0.0, 0.0, 1.0)}[cm]

    def field(n, x, f):
        return np.r_[np.asarray(e) + _d(f) * np.array([1e-3 * x[2], 2e-3 * x[0] * x[1], -1e-3 * x[1]]), 0, 0, 0]

    _file8(tmp_path, field, cm=cm, ang=ang)
    _harmonic(tmp_path)
    dk = _stress_deck(itran=1, cm=cm, ang=ang)
    _eout(dk, 1, [1, 2], [1] * 7)
    rc, out = _run(tmp_path, dk)
    assert rc == 0, out[-2000:]
    f4 = read_container(tmp_path / "m.N4", "FILE4")
    f8 = read_container(tmp_path / "FILE8", "FILE8")
    H = np.asarray(f8["H"])
    for row in range(2):
        Sx, eq = f4["rec_SOLID_S"][row], f4["rec_SOLID_eq"][row]
        g = int(f4["elem_group"][f4["rec_SOLID_idx"][row]])
        el = int(f4["elem_id"][f4["rec_SOLID_idx"][row]])
        for c, comp in enumerate(("SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ")):
            _, h = _tfu(tmp_path, "SOLID", g, el, comp)
            plain = np.zeros(len(FNUM), dtype=complex)        # explicit loop: sum_d S_cd U_d
            big = np.zeros(len(FNUM))
            for k, e_k in enumerate(eq):
                if e_k >= 0:
                    plain += Sx[c, k] * H[:, e_k]
                    big += np.abs(Sx[c, k] * H[:, e_k])
            assert np.all(np.abs(h - plain) <= 1e-9 * big), (row, comp)


# =============================================================================================
# 6. listing, ELEMENT_CENTER, Frames.txt, errors
# =============================================================================================
def test_listing_time_of_maximum_of_a_harmonic(tmp_path):
    """Spring force under a0 sin(w0 t) with a real STF s: F = s u_g = (s g a0/w0^2) sin(w0 t); the peaks of
    |F| are at t = T0/4 + n T0/2, on the 0.02 s grid here (equal to round-off, so any of them may be listed)."""
    _spring_model(tmp_path, damp=0.0)
    _file8(tmp_path, lambda n, x, f: np.r_[1.0 + (0.5 if n == 2 else 0.0), 0, 0, 0, 0, 0])
    k0 = 8                                                      # T0 = NFFT dt / 8 = 0.64 s, T0/4 = 0.16 s
    _harmonic(tmp_path, k0=k0, a0=0.1, kind="sin")
    d = _stress_deck()
    _eout(d, 1, [1], [1])
    rc, out = _run(tmp_path, d)
    assert rc == 0, out[-1500:]
    row = [ln for ln in out.splitlines() if ln.split()[:2] == ["1", "FX"]]
    assert len(row) == 1
    v, t = float(row[0].split()[2]), float(row[0].split()[3])
    w0 = 2 * np.pi * k0 * DF
    assert v == pytest.approx(SPRING_K[0] * 0.5 * G_SI * 0.1 / w0 ** 2, rel=1e-6)
    n = (t - 0.16) / 0.32
    assert abs(n - round(n)) < 1e-6, t


def test_element_center_ordered_group_numbers_and_layout(tmp_path):
    """D-FIL-06: 'ordered group #' = 1-based order among the groups of the same element type; group blocks in
    ascending group number (manual layout, spec 05d 1.9)."""
    assert SL.ordered_group_numbers([(1, "SOLID"), (2, "SHELL"), (3, "SOLID"), (5, "SPRING"), (4, "SOLID")]) == \
        {1: 1, 2: 1, 3: 2, 4: 3, 5: 1}
    d = house_deck()
    nid = {}
    n = 0
    for ix in range(3):
        for iy in range(2):
            for iz in range(2):
                n += 1
                nid[(ix, iy, iz)] = n
                add_node(d, n, float(ix), float(iy), float(iz))
    d.table("materials").append([1, 3, 900.0, 400.0, 20.0, 0.05, 0.05])
    d.table("groups").append([7, 1, "b"])
    d.table("groups").append([3, 1, "a"])
    for g, ix in ((7, 1), (3, 0)):
        add_element(d, g, 1, [nid[(ix, 0, 0)], nid[(ix + 1, 0, 0)], nid[(ix + 1, 1, 0)], nid[(ix, 1, 0)],
                              nid[(ix, 0, 1)], nid[(ix + 1, 0, 1)], nid[(ix + 1, 1, 1)], nid[(ix, 1, 1)]], mat=1)
    _house(tmp_path, d)
    _file8(tmp_path, lambda n, x, f: np.r_[1.0 + _d(f) * 1e-4 * x[2] * (1 + x[0]), 0, 0, 0, 0, 0])
    _harmonic(tmp_path)
    dk = _stress_deck(savemax=1)
    _eout(dk, 3, [1], [1])
    rc, out = _run(tmp_path, dk)
    assert rc == 0, out[-1500:]
    lines = (tmp_path / "ELEMENT_CENTER_ABS_MAX_STRESSES.TXT").read_text().split("\n")
    assert lines[0].strip() == "2"
    assert lines[1].split() == ["SOLID", "3", "1", "1"] and lines[2].split() == ["SOLID", "7", "2", "1"]
    blocks = SL.read_element_center(tmp_path / "ELEMENT_CENTER_ABS_MAX_STRESSES.TXT")
    assert all(b.values.shape == (1, 6) and np.all(b.values >= 0) for b in blocks)


def test_frames_file_manual_example(tmp_path):
    """spec 05d section 8 item 9: the manual example gives frames 1..10 and soil groups 8, 9, 10."""
    p = tmp_path / "Frames.txt"
    p.write_text("10\n" + "\n".join(str(i) for i in range(1, 11)) + "\n3\n8 9 10\n")
    assert SL.read_frames_file(p) == (list(range(1, 11)), [8, 9, 10])


def test_error_80_for_a_general_matrix_group_and_data_check_writes_nothing(tmp_path):
    """GENERAL (GM) elements have no STRESS components: a request on their group is Error 80; the data-check
    mode (opmode 1) writes no output file at all (no .THS, .TFU, FILE14, FILE15)."""
    from sassi.verify.problems.vp_house import add_matrix_property
    d = house_deck()
    add_node(d, 1, 0.0, 0.0, 0.0)
    add_node(d, 2, 1.0, 0.0, 0.0)
    rng = np.random.default_rng(1)
    B = rng.standard_normal((12, 6))
    with np.errstate(all="ignore"):                          # spurious Accelerate flags (blas_quiet)
        KB = B @ B.T
    add_matrix_property(d, 1, "R", 1e4 * KB)
    add_matrix_property(d, 1, "I", 4e2 * KB)
    add_matrix_property(d, 1, "M", np.eye(12))
    d.table("springprops").append([1, *SPRING_K, SPRING_DAMP])
    d.table("groups").append([1, 9, "gm"])
    d.table("groups").append([2, 7, "spring"])
    add_element(d, 1, 1, [1, 2], prop=1)
    add_element(d, 2, 1, [1, 2], prop=1)
    rc, out = run_house(tmp_path, d)
    if rc != 0:
        pytest.skip("HOUSE deck for a GENERAL element not accepted by this builder: " + out[-300:])
    _file8(tmp_path, lambda n, x, f: np.r_[1.0 + (0.1 * _d(f) if n == 2 else 0.0), 0, 0, 0, 0, 0])
    _harmonic(tmp_path)
    dk = _stress_deck()
    _eout(dk, 1, [1], [1])
    rc, out = _run(tmp_path, dk)
    assert rc == 1 and "Error 80" in out
    before = set(p.name for p in tmp_path.iterdir())
    dk = _stress_deck(opmode=1, itran=1)
    _eout(dk, 2, [1], [2] * 6)
    rc, out = _run(tmp_path, dk)
    assert rc == 0
    after = set(p.name for p in tmp_path.iterdir())
    assert after - before <= {"m.str", "m_STRESS.out"}


# =============================================================================================
# 7. defects found by the review (fixed; the strict-xfail markers were removed after the fix)
# =============================================================================================
def _shell_on_excavated_soil(wd):
    """SHELL basemat (group 1, z = 0) on two excavated SOLID soil elements (group 2, ETYPE 2) down to the
    first layer interface (z = -2): a structure made of SHELL elements only."""
    d = house_deck()
    nid = {}
    n = 0
    for iz, z in enumerate((0.0, -2.0)):
        for iy, y in enumerate((0.0, 1.0)):
            for ix, x in enumerate((0.0, 1.0, 2.0)):
                n += 1
                nid[(ix, iy, iz)] = n
                add_node(d, n, x, y, z)
    for k in nid.values():
        d.table("interaction").append([k])
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("groups").append([1, 3, "mat"])
    d.table("groups").append([2, 1, "excavated soil"])
    for e, ix in ((1, 0), (2, 1)):
        add_element(d, 1, e, [nid[(ix, 0, 0)], nid[(ix + 1, 0, 0)], nid[(ix + 1, 1, 0)], nid[(ix, 1, 0)]], mat=1,
                    thick=0.5)
        add_element(d, 2, e, [nid[(ix, 0, 1)], nid[(ix + 1, 0, 1)], nid[(ix + 1, 1, 1)], nid[(ix, 1, 1)],
                              nid[(ix, 0, 0)], nid[(ix + 1, 0, 0)], nid[(ix + 1, 1, 0)], nid[(ix, 1, 0)]],
                    etype=2, mat=1)
    _house(wd, d)
    _file8(wd, lambda n, x, f: np.r_[1.0 + 0.01 * x[2] * x[0] * _d(f), 0.002 * x[1] * x[0] * _d(f),
                                     0.003 * x[0] ** 2 * _d(f), 0.001 * x[1] * _d(f), 0.002 * x[0] * _d(f), 0.0])
    write_motion_file(wd / "eq.acc", 0.2 * synthetic_motion(200, DT, seed=4, fmax=15.0), DT)
    (wd / "Frames.txt").write_text("1\n10\n0\n")


@pytest.mark.parametrize("extra", [{"savemax": 1}, {"secdataopt": 1}])
def test_defect_shell_bending_frames_do_not_depend_on_unrelated_options(tmp_path, extra):
    """spec 05d 1.9 / D-STR-08: 'If the model has SHELL elements only, separate bending frames (bd files) are
    also generated'.  Excavated soil is not part of the structure (it is left out of the nodal frames), so a
    SHELL basemat on excavated soil must get its bdsig/bdtau frames whatever the all-element options."""
    _shell_on_excavated_soil(tmp_path)
    dk = _stress_deck(rstns=1)
    _eout(dk, 1, [1, 2], [1] * 6)
    assert _run(tmp_path, dk)[0] == 0
    assert (tmp_path / "NSTRESS" / "stress_ABS_MAX_bdsig").exists()            # rstns alone: written
    for p in (tmp_path / "NSTRESS").iterdir():
        p.unlink()
    dk = _stress_deck(rstns=1, **extra)
    _eout(dk, 1, [1, 2], [1] * 6)
    rc, out = _run(tmp_path, dk)
    assert rc == 0, out[-1500:]
    assert (tmp_path / "NSTRESS" / "stress_ABS_MAX_bdsig").exists()
    assert (tmp_path / "NSTRESS" / "stress_ABS_MAX_bdtau").exists()


def test_defect_ssi_frequency_input_limit_is_reported(tmp_path):
    """requirements 4.0.1: 'Maximum SSI frequencies: ... 1,500 in COMBIN output and MOTION/RELDISP/STRESS
    input [M]'.  MOTION reports EDU-02 above the limit; STRESS accepts 1,501 frequencies silently."""
    _spring_model(tmp_path)
    f4 = read_container(tmp_path / "m.N4", "FILE4")
    nfft, dt = 4096, DT
    df = 1.0 / (nfft * dt)
    fnum = np.arange(1, 1502)
    H = np.zeros((len(fnum), len(f4["eq_node"])), complex)
    H[:, :] = 1.0
    write_file8(tmp_path / "FILE8", fnum, df, f4["eq_node"], f4["eq_dof"], H, nfft=nfft, delt=dt,
                model_hash=str(f4.meta.get("model_hash", "")))
    write_motion_file(tmp_path / "eq.acc", 0.1 * synthetic_motion(300, dt, seed=1), dt)
    d = _stress_deck()
    d["nft"], d["df"] = nfft, df
    _eout(d, 1, [1], [1])
    rc, out = _run(tmp_path, d)
    assert rc == 0
    assert "1500" in out or "1,500" in out or "EDU-02" in out
