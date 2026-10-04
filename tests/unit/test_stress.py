"""Unit tests of the STRESS module (requirements 4.10, D-STR-01 ... D-STR-12, D-FIL-06).

The models are written as HOUSE decks (FILE4 with the recovery operators) and the nodal transfer
functions are synthetic FILE8 containers, so every test exercises deck -> HOUSE -> STRESS -> files.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from sassi.core import interp as I
from sassi.core import signal as S
from sassi.core import stress_lib as SL
from sassi.elements.base import blas_quiet
from sassi.io.files import read_container
from sassi.verify.problems.vp_house import G_SI, add_element, add_node, house_deck, run_house
from sassi.verify.problems.vp_motion import synthetic_motion, write_file8, write_motion_file
from sassi.verify.problems.vp_stress import (add_eout, field_file8, file4_maps, mixed_model_deck, read_ths, read_tfu,
                                             run_stress, stress_deck)

NFFT, DT = 512, 0.01
DF = 1.0 / (NFFT * DT)
FNUM = np.unique(np.concatenate([np.arange(1, 8), np.arange(8, 120, 4)])).astype(int)
FN, BETA = 4.0, 0.05


def _d(f):
    """Relative deformation factor of a damped mode: (f/fn)^2 / (1 - (f/fn)^2 + 2 i beta)."""
    r2 = (f / FN) ** 2
    return r2 / (1.0 - r2 + 2j * BETA)


def deform_field(node, x, f, rigid=True):
    """Total-motion TF: unit rigid X motion plus a smooth non-rigid field scaled by _d(f)."""
    phi = np.array([0.01 * x[2] ** 2 + 0.004 * x[0] * x[1], 0.003 * x[0] * x[2], 0.002 * x[0] ** 2 - 0.001 * x[1] * x[2],
                    0.002 * x[1], -0.003 * x[0], 0.001 * x[2]])
    return (np.array([1.0, 0, 0, 0, 0, 0]) if rigid else 0.0) + _d(f) * phi


@pytest.fixture
def mixed(tmp_path):
    rc, out = run_house(tmp_path, mixed_model_deck())
    assert rc == 0, out[-2000:]
    H = field_file8(tmp_path, FNUM, DF, deform_field, NFFT, DT)
    acc = 0.25 * synthetic_motion(380, DT, seed=4, fmax=20.0)
    write_motion_file(tmp_path / "eq.acc", acc, DT)
    return tmp_path, H, acc


def _listing(wd):
    return (Path(wd) / "m_STRESS.out").read_text(encoding="utf-8")


def _reference_history(wd, H, acc, etype, group, elem, comp, mode="seismic", option=1, nout=NFFT):
    """irfft(interp(STF) * U_g) computed independently of the module from FILE4 / FILE8."""
    f4, _, _ = file4_maps(wd)
    t = f4.meta["components"][etype]
    idx = [k for k, g in enumerate(f4[f"rec_{etype}_idx"])
           if int(f4["elem_group"][g]) == group and int(f4["elem_id"][g]) == elem][0]
    U = SL.element_dof_tf(H, f4[f"rec_{etype}_eq"][idx:idx + 1])[:, 0, :]
    with blas_quiet():
        stf = U @ f4[f"rec_{etype}_S"][idx].T                    # (nF, nc)
    c = t.index(comp)
    f_grid = S.fourier_grid(NFFT, DT)
    Hg = I.interpolate_tf(FNUM * DF, stf[:, c], f_grid, option, mode=mode)
    A = np.fft.rfft(S.pad_record(acc, NFFT))
    if mode == "seismic":
        R = S.hermitian_bins(Hg, NFFT) * S.displacement_spectrum(A, f_grid, G_SI)
    else:
        Hg = np.array(Hg)
        Hg[-1] = Hg[-1].real
        R = Hg * A
    return np.fft.irfft(S.hermitian_bins(R, NFFT), NFFT)[:nout], stf[:, c], Hg


# =============================================================================================
# convolution and transfer-function files
# =============================================================================================
def test_seismic_histories_and_tf_files(mixed):
    wd, H, acc = mixed
    d = stress_deck(NFFT, DT, itran=1)
    add_eout(d, 1, [1], [2] * 7)
    add_eout(d, 2, [1], [2] * 12)
    add_eout(d, 3, [1], [2] * 6)
    add_eout(d, 4, [1], [2] * 6)
    rc, out = run_stress(wd, d)
    assert rc == 0, out[-3000:]
    for et, g, comp in (("SOLID", 1, "SXZ"), ("BEAMS", 2, "MZI"), ("SHELL", 3, "MXX"), ("SPRING", 4, "FY")):
        ref, stf, Hg = _reference_history(wd, H, acc, et, g, 1, comp)
        th, dt = read_ths(wd, et, g, 1, comp)
        assert dt == DT and len(th) == NFFT
        assert np.max(np.abs(th - ref)) <= 1e-9 * np.max(np.abs(ref)), (et, comp)
        f, tfu = read_tfu(wd, et, g, 1, comp)
        assert np.allclose(f, FNUM * DF) and np.allclose(tfu, stf, rtol=1e-9, atol=1e-12 * np.max(np.abs(stf)))
        fi, tfi = read_tfu(wd, et, g, 1, comp, "TFI")
        k = np.rint(fi / DF).astype(int)
        assert k[0] == 0 and k[-1] == FNUM[-1]                    # TFI rows 0 .. f_N
        assert np.allclose(tfi[FNUM], tfu, rtol=1e-12, atol=0)    # TFI = TFU at the SSI frequencies
        assert abs(tfi[0]) == 0.0                                 # rigid body: no stress at f = 0
    # SOCT is computed from the six histories in the time domain
    six = [read_ths(wd, "SOLID", 1, 1, c)[0] for c in ("SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ")]
    soct, _ = read_ths(wd, "SOLID", 1, 1, "SOCT")
    assert np.allclose(soct, SL.octahedral_shear_stress(*six), rtol=1e-12, atol=0)
    assert not (wd / "SOLID_001_00001_SOCT.TFU").exists()
    txt = _listing(wd)
    assert "Maximum absolute response, group 2 (BEAMS)" in txt and "time (s)" in txt
    f14 = (wd / "FILE14").read_text()
    assert "BEAMS group 2 element 1 MZI" in f14
    assert "SOLID group 1 element 1 SOCT" in (wd / "FILE15").read_text()


def test_vibration_convolution(tmp_path):
    rc, _ = run_house(tmp_path, mixed_model_deck())
    assert rc == 0
    H = field_file8(tmp_path, FNUM, DF, lambda n, x, f: deform_field(n, x, f, rigid=False), NFFT, DT, type=1)
    load = synthetic_motion(300, DT, seed=9, fmax=20.0)
    write_motion_file(tmp_path / "eq.acc", load, DT)
    d = stress_deck(NFFT, DT, type=1)
    add_eout(d, 2, [1], [2] * 12)
    rc, out = run_stress(tmp_path, d)
    assert rc == 0, out[-2000:]
    ref, _, _ = _reference_history(tmp_path, H, load, "BEAMS", 2, 1, "FYJ", mode="vibration")
    th, _ = read_ths(tmp_path, "BEAMS", 2, 1, "FYJ")
    assert np.max(np.abs(th - ref)) <= 1e-10 * np.max(np.abs(ref))


def test_interpolation_options_and_phase_adjustment(mixed):
    wd, H, acc = mixed
    for opt in (0, 2, 6):
        d = stress_deck(NFFT, DT, interopt=opt, itran=1)
        add_eout(d, 1, [1], [0, 0, 0, 0, 2])
        assert run_stress(wd, d)[0] == 0
        ref, _, _ = _reference_history(wd, H, acc, "SOLID", 1, 1, "SXZ", option=opt)
        th, _ = read_ths(wd, "SOLID", 1, 1, "SXZ")
        assert np.max(np.abs(th - ref)) <= 1e-9 * np.max(np.abs(ref)), opt
    d = stress_deck(NFFT, DT, interopt=2, pzadj=1, smo=1000.0, itran=1)
    add_eout(d, 1, [1], [0, 0, 0, 0, 1])
    assert run_stress(wd, d)[0] == 0
    _, tfi = read_tfu(wd, "SOLID", 1, 1, "SXZ", "TFI")
    _, tfu = read_tfu(wd, "SOLID", 1, 1, "SXZ")
    assert np.allclose(np.abs(tfi[FNUM]), np.abs(tfu), rtol=1e-12)       # phase adjustment keeps |STF|


def test_skip_save_itran_and_duration(mixed):
    wd, H, acc = mixed
    d = stress_deck(NFFT, DT)
    add_eout(d, 2, [1], [2] * 12)
    assert run_stress(wd, d)[0] == 0
    full, _ = read_ths(wd, "BEAMS", 2, 1, "MYJ")
    assert not (wd / "BEAMS_002_00001_MYJ.TFU").exists() and not (wd / "FILE14").exists()
    d = stress_deck(NFFT, DT, skip=3, dur=2.0)
    add_eout(d, 2, [1], [2] * 12)
    assert run_stress(wd, d)[0] == 0
    th, dt = read_ths(wd, "BEAMS", 2, 1, "MYJ")
    nout = int(round(1.2 * 2.0 / DT))
    assert dt == pytest.approx(3 * DT) and np.allclose(th, full[:nout:3], rtol=0, atol=0)
    for p in wd.glob("*.THS"):
        p.unlink()
    d = stress_deck(NFFT, DT, save=0)
    add_eout(d, 2, [1], [0, 2])
    assert run_stress(wd, d)[0] == 0
    assert not list(wd.glob("*.THS"))
    assert "BEAMS group 2 element 1 FYI" in (wd / "FILE15").read_text()


def test_mixed_group_beams_identical(mixed):
    """D-STR-11: BEAMS results are the same whether requested alone or together with other groups."""
    wd, H, acc = mixed
    d = stress_deck(NFFT, DT)
    add_eout(d, 2, [1], [2] * 12)
    assert run_stress(wd, d)[0] == 0
    alone = {c: read_ths(wd, "BEAMS", 2, 1, c)[0] for c in ("FXI", "FXJ", "MZI")}
    d = stress_deck(NFFT, DT, savemax=1)
    add_eout(d, 1, [1], [2] * 7)
    add_eout(d, 2, [1], [2] * 12)
    add_eout(d, 3, [1], [1] * 6)
    assert run_stress(wd, d)[0] == 0
    for c, v in alone.items():
        assert np.array_equal(read_ths(wd, "BEAMS", 2, 1, c)[0], v)


def test_aux_components_for_soct_and_largest_value(mixed):
    wd, H, acc = mixed
    d = stress_deck(NFFT, DT)
    add_eout(d, 1, [1], [0, 0, 0, 0, 0, 0, 1])
    assert run_stress(wd, d)[0] == 0
    txt = _listing(wd)
    assert txt.count("(computed for SOCT/EOCT, not requested)") == 6
    assert "largest requested value" in txt and "SOCT" in txt
    assert not list(wd.glob("*.THS"))


# =============================================================================================
# errors and checks
# =============================================================================================
def _expect_error(wd, d, text):
    rc, out = run_stress(wd, d)
    assert rc == 1 and text in out, out[-1500:]


def test_request_errors(mixed):
    wd, _, _ = mixed
    _expect_error(wd, stress_deck(NFFT, DT), "Error 79")
    d = stress_deck(NFFT, DT)
    add_eout(d, 99, [1], [1])
    _expect_error(wd, d, "Error 80")
    d = stress_deck(NFFT, DT)
    add_eout(d, 1, [5], [1])
    _expect_error(wd, d, "Error 81")
    d = stress_deck(NFFT, DT)
    add_eout(d, 1, [1], [3])
    _expect_error(wd, d, "0, 1 or 2")
    d = stress_deck(NFFT, DT)
    add_eout(d, 1, [1], [1, 0])
    add_eout(d, 1, [1], [0, 2])
    add_eout(d, 4, [1], [1] * 8)
    rc, out = run_stress(wd, d)
    assert rc == 0 and "Error 82" in out and "beyond the components" in out
    assert (wd / "SOLID_001_00001_SYY.THS").exists()                   # merged code 2 for SYY
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        d = stress_deck(NFFT, DT, skip=-1)
    add_eout(d, 1, [1], [1])
    _expect_error(wd, d, "Error 67")


def test_file_consistency_errors(mixed):
    wd, H, acc = mixed
    f4 = read_container(wd / "m.N4", "FILE4")
    d = stress_deck(NFFT, DT)
    add_eout(d, 1, [1], [1])
    write_file8(wd / "FILE8", FNUM, 2 * DF, f4["eq_node"], f4["eq_dof"], H, nfft=NFFT, delt=DT)
    _expect_error(wd, d, "D-MOT-14")
    write_file8(wd / "FILE8", FNUM, DF, f4["eq_node"], f4["eq_dof"], H, nfft=NFFT, delt=DT, type=1)
    _expect_error(wd, d, "does not match")
    perm = np.roll(np.arange(len(f4["eq_node"])), 1)
    write_file8(wd / "FILE8", FNUM, DF, f4["eq_node"][perm], f4["eq_dof"][perm], H, nfft=NFFT, delt=DT)
    _expect_error(wd, d, "different equation maps")
    (wd / "FILE8").unlink()
    _expect_error(wd, d, "run ANALYS")
    (wd / "m.N4").unlink()
    _expect_error(wd, d, "run HOUSE")


def test_data_check_mode_and_iter_warning(mixed):
    wd, _, _ = mixed
    d = stress_deck(NFFT, DT, opmode=1)
    add_eout(d, 1, [1], [2] * 7)
    rc, out = run_stress(wd, d)
    assert rc == 0 and "data-check mode writes no output files" in out and not list(wd.glob("*.THS"))
    d = stress_deck(NFFT, DT, iter=1)
    add_eout(d, 1, [1], [1])
    rc, out = run_stress(wd, d)
    assert rc == 0 and "FILE74" in out and not (wd / "FILE74").exists()


def test_control_motion_errors(mixed):
    wd, _, _ = mixed
    d = stress_deck(NFFT, DT, thfile="missing.acc")
    add_eout(d, 1, [1], [1])
    _expect_error(wd, d, "Error 73")
    d = stress_deck(NFFT, DT)
    d["mult"] = 0.0
    add_eout(d, 1, [1], [1])
    _expect_error(wd, d, "Error 77")


# =============================================================================================
# strain output (D-STR-04)
# =============================================================================================
def test_strain_output_with_strainout(mixed):
    wd, H, acc = mixed
    d = stress_deck(NFFT, DT, itran=1)
    add_eout(d, 1, [1], [2] * 7)
    run_stress(wd, d)                                               # writes the deck
    p = wd / "m.str"
    p.write_text(p.read_text().replace("[params]\n", "[params]\nstrainout = 1\n", 1))
    from sassi.modules.base import run_module
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert run_module("STRESS", "m", wd) == 0
    f4, _, _ = file4_maps(wd)
    B = f4["x_rec_SOLID_B"][0]
    U = SL.element_dof_tf(H, f4["rec_SOLID_eq"][:1])[:, 0, :]
    with blas_quiet():
        eps = U @ B.T                                               # (nF, 6) engineering strains
    for c, name in enumerate(("EXX", "EYY", "EZZ", "EXY", "EXZ", "EYZ")):
        _, tfu = read_tfu(wd, "SOLID", 1, 1, name)
        assert np.allclose(tfu, eps[:, c], rtol=1e-9, atol=1e-12 * np.max(np.abs(eps)))
    six = [read_ths(wd, "SOLID", 1, 1, n)[0] for n in ("EXX", "EYY", "EZZ", "EXY", "EXZ", "EYZ")]
    eoct, _ = read_ths(wd, "SOLID", 1, 1, "EOCT")
    assert np.allclose(eoct, SL.octahedral_shear_strain(*six), rtol=1e-12, atol=0)
    # stress = D* strain (SOLID): SXZ = G* EXZ
    from sassi.verify.problems.vp_stress import MAT_PD, material
    _, sxz = read_tfu(wd, "SOLID", 1, 1, "SXZ")
    _, exz = read_tfu(wd, "SOLID", 1, 1, "EXZ")
    assert np.allclose(sxz, material(MAT_PD).G * exz, rtol=1e-9)


# =============================================================================================
# all-element outputs
# =============================================================================================
def test_savemax_saveth_and_ess(mixed):
    wd, H, acc = mixed
    d = stress_deck(NFFT, DT, savemax=1, saveth=1, secdataopt=1, skip=64)
    add_eout(d, 1, [1], [2] * 7)
    add_eout(d, 3, [1], [2] * 6)
    add_eout(d, 4, [1], [2] * 6)
    rc, out = run_stress(wd, d)
    assert rc == 0, out[-2000:]
    assert "BEAMS excluded" in out and "m_ABS_MAX.bdsig" in out
    full = {}
    d2 = stress_deck(NFFT, DT)
    add_eout(d2, 1, [1], [2] * 7)
    add_eout(d2, 3, [1], [2] * 6)
    add_eout(d2, 4, [1], [2] * 6)
    assert run_stress(wd, d2)[0] == 0                              # full-length histories (skip 0)
    for et, g, comps in (("SOLID", 1, SL.CENTER_COLUMNS["SOLID"]), ("SHELL", 3, SL.CENTER_COLUMNS["SHELL"]),
                         ("SPRING", 4, SL.CENTER_COLUMNS["SPRING"])):
        full[et] = np.column_stack([read_ths(wd, et, g, 1, c)[0] for c in comps])
    blocks = SL.read_element_center(wd / "ELEMENT_CENTER_ABS_MAX_STRESSES.TXT")
    assert [(b.etype, b.group, b.ordered) for b in blocks] == [("SOLID", 1, 1), ("SHELL", 3, 1), ("SPRING", 4, 1)]
    for b in blocks:
        assert np.allclose(b.values[0], np.max(np.abs(full[b.etype]), axis=0), rtol=1e-9)
    sig = SL.read_element_center(wd / "m_ABS_MAX.sig")
    assert [(b.etype, b.values.shape[1]) for b in sig] == [("SOLID", 3), ("SHELL", 3)]
    assert sig[1].values[0, 2] == 0.0 and np.allclose(sig[1].values[0, :2], np.max(np.abs(full["SHELL"][:, :2]), 0))
    bd = SL.read_element_center(wd / "m_ABS_MAX.bdtau")
    assert bd[0].etype == "SHELL" and np.allclose(bd[0].values[0, 0], np.max(np.abs(full["SHELL"][:, 5])))
    # all-element histories, every 64th step
    M = np.loadtxt(wd / "m.tau")
    steps = np.arange(0, NFFT, 64)
    assert M.shape == (len(steps), 1 + 3 + 1)                       # time, SOLID SXY SXZ SYZ, SHELL FXY
    assert np.allclose(M[:, 0], steps * DT) and np.allclose(M[:, 1:4], full["SOLID"][steps, 3:6], rtol=1e-9,
                                                               atol=1e-12)
    assert "SHELL_003_00001_FXY" in (wd / "m.tau").read_text().splitlines()[2]
    # ESTRESS frames (signed values, 1-based step numbers)
    lstf = (wd / "NSTRESS" / "ESTRESS.lst").read_text().splitlines()
    assert lstf[1:] == [f"ESTRESS_{s + 1:05d}.ess" for s in steps]
    fr = SL.read_element_center(wd / "NSTRESS" / "ESTRESS_00065.ess")
    assert np.allclose(fr[0].values[0], full["SOLID"][64], rtol=1e-9, atol=1e-12 * np.max(np.abs(full["SOLID"])))
    assert np.allclose(fr[2].values[0], full["SPRING"][64], rtol=1e-9, atol=1e-12 * np.max(np.abs(full["SPRING"])))


def _two_block_model(wd, shell_wall=False):
    """Structural SOLID (group 1, x 0..1) next to a near-field soil SOLID (group 2, x 1..2)."""
    d = house_deck()
    nid = {}
    n = 0
    for ix, x in enumerate((0.0, 1.0, 2.0)):
        for iy, y in enumerate((0.0, 1.0)):
            for iz, z in enumerate((0.0, 1.0)):
                n += 1
                nid[(ix, iy, iz)] = n
                add_node(d, n, x + 0.05 * y * z, y, z)
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("materials").append([2, 3, 500.0, 250.0, 19.0, 0.04, 0.04])
    d.table("groups").append([1, 1, "wall block"])
    d.table("groups").append([2, 1, "near-field soil"])
    for g, ix, mat in ((1, 0, 1), (2, 1, 2)):
        ns = [nid[(ix, 0, 0)], nid[(ix + 1, 0, 0)], nid[(ix + 1, 1, 0)], nid[(ix, 1, 0)],
              nid[(ix, 0, 1)], nid[(ix + 1, 0, 1)], nid[(ix + 1, 1, 1)], nid[(ix, 1, 1)]]
        add_element(d, g, 1, ns, etype=1, mat=mat)
    rc, out = run_house(wd, d)
    assert rc == 0, out[-1500:]
    return nid


def test_nodal_stress_and_soil_pressure_frames(tmp_path):
    wd = tmp_path
    nid = _two_block_model(wd)
    H = field_file8(wd, FNUM, DF, deform_field, NFFT, DT)
    write_motion_file(wd / "eq.acc", 0.25 * synthetic_motion(380, DT, seed=4, fmax=20.0), DT)
    d = stress_deck(NFFT, DT, rstns=1, rstsp=1)
    add_eout(d, 1, [1], [2] * 6)
    add_eout(d, 2, [1], [2] * 6)
    _expect_error(wd, d, "Frames.txt")
    (wd / "Frames.txt").write_text("3\n1\n101\n9999\n1\n2\n")
    rc, out = run_stress(wd, d)
    assert rc == 0, out[-2000:]
    assert "outside" in out                                            # frame 9999 ignored
    s1 = np.column_stack([read_ths(wd, "SOLID", 1, 1, c)[0] for c in SL.CENTER_COLUMNS["SOLID"]])
    s2 = np.column_stack([read_ths(wd, "SOLID", 2, 1, c)[0] for c in SL.CENTER_COLUMNS["SOLID"]])
    # nodal stress frames: plain mean of the adjacent element centres (D-STR-08)
    fr = np.loadtxt(wd / "NSTRESS" / f"stress_{1.0:06.3f}_00101_sig", skiprows=1)
    rows = {int(r[0]): r[1:] for r in fr}
    shared, only1, only2 = nid[(1, 0, 0)], nid[(0, 0, 0)], nid[(2, 1, 1)]
    assert np.allclose(rows[shared], 0.5 * (s1[100, :3] + s2[100, :3]), rtol=1e-9)
    assert np.allclose(rows[only1], s1[100, :3], rtol=1e-9) and np.allclose(rows[only2], s2[100, :3], rtol=1e-9)
    mx = np.loadtxt(wd / "NSTRESS" / "stress_ABS_MAX_tau", skiprows=1)
    assert np.allclose({int(r[0]): r[1:] for r in mx}[only1], np.max(np.abs(s1[:, 3:6]), axis=0), rtol=1e-9)
    assert not (wd / "NSTRESS" / "stress_ABS_MAX_bdsig").exists()     # SOLID model: no bending frames
    # soil pressure on the shared face x = 1 (normal X): p = -SXX of the soil element (D-STR-09)
    ele = np.loadtxt(wd / "SOILPRES" / f"pres_{0.0:06.3f}_00001_ele", skiprows=1, ndmin=2)
    assert ele.shape == (1, 3) and ele[0, 0] == 2 and ele[0, 1] == 1
    n = SL.face_normal(np.array([[1.0, 0, 0], [1.0, 1, 0], [1.05, 1, 1], [1.0, 0, 1]]))
    p_ref = -SL.normal_stress(s2, n)
    ele = np.loadtxt(wd / "SOILPRES" / f"pres_{1.0:06.3f}_00101_ele", skiprows=1, ndmin=2)
    assert ele[0, 2] == pytest.approx(p_ref[100], rel=1e-9)
    nod = np.loadtxt(wd / "SOILPRES" / f"pres_{1.0:06.3f}_00101_nod", skiprows=1)
    assert sorted(int(v) for v in nod[:, 0]) == sorted(nid[(1, iy, iz)] for iy in (0, 1) for iz in (0, 1))
    assert np.allclose(nod[:, 1], p_ref[100], rtol=1e-9)
    stat = SL.read_element_center(wd / "STATIC_SOIL_PRESSURES.TXT")
    assert stat[0].group == 2 and stat[0].values.shape == (1, 1) and stat[0].values[0, 0] == 0.0
    pmax = SL.read_element_center(wd / "pres_max_ele")
    assert pmax[0].values[0, 0] == pytest.approx(np.max(np.abs(p_ref)), rel=1e-9)
    # static pressure added algebraically on the next run
    stat[0].values[0, 0] = 10.0
    SL.write_element_center(wd / "STATIC_SOIL_PRESSURES.TXT", stat)
    assert run_stress(wd, d)[0] == 0
    ele = np.loadtxt(wd / "SOILPRES" / f"pres_{1.0:06.3f}_00101_ele", skiprows=1, ndmin=2)
    assert ele[0, 2] == pytest.approx(10.0 + p_ref[100], rel=1e-9)
    assert SL.read_element_center(wd / "pres_max_ele")[0].values[0, 0] == pytest.approx(np.max(np.abs(10 + p_ref)))


def test_shell_only_bending_frames(tmp_path):
    d = house_deck()
    for n, (x, y) in enumerate(((0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)), start=1):
        add_node(d, n, x, y, 0.0, fix=(0, 0, 0, 0, 0, 1))
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("groups").append([1, 3, "slab"])
    add_element(d, 1, 1, [1, 2, 5, 4], thick=0.3)
    add_element(d, 1, 2, [2, 3, 6, 5], thick=0.3)
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out[-1500:]
    field_file8(tmp_path, FNUM, DF, deform_field, NFFT, DT)
    write_motion_file(tmp_path / "eq.acc", 0.25 * synthetic_motion(380, DT, seed=4, fmax=20.0), DT)
    (tmp_path / "Frames.txt").write_text("1\n50\n0\n")
    dk = stress_deck(NFFT, DT, rstns=1)
    add_eout(dk, 1, [1, 2], [2] * 6)
    assert run_stress(tmp_path, dk)[0] == 0
    for cls in ("sig", "tau", "bdsig", "bdtau"):
        assert (tmp_path / "NSTRESS" / f"stress_{0.49:06.3f}_00050_{cls}").exists()
    fr = np.loadtxt(tmp_path / "NSTRESS" / f"stress_{0.49:06.3f}_00050_bdtau", skiprows=1)
    m1, m2 = read_ths(tmp_path, "SHELL", 1, 1, "MXY")[0], read_ths(tmp_path, "SHELL", 1, 2, "MXY")[0]
    rows = {int(r[0]): r[1] for r in fr}
    assert rows[2] == pytest.approx(0.5 * (m1[49] + m2[49]), rel=1e-9) and rows[1] == pytest.approx(m1[49], rel=1e-9)


def test_batch_protocol(mixed, monkeypatch):
    import io
    from sassi.modules import stress
    wd, _, _ = mixed
    d = stress_deck(NFFT, DT)
    add_eout(d, 1, [1], [1])
    from sassi.io import decks
    decks.write(wd / "m.str", d)
    monkeypatch.chdir(wd)
    rc = stress.batch_main("STRESS", io.StringIO("m\nm.str\nm_STRESS.out\n"))
    assert rc == 0 and "SXX" in (wd / "m_STRESS.out").read_text()


def test_model_hash_warning_and_2d_frames(mixed, tmp_path_factory):
    wd, H, _ = mixed
    f4 = read_container(wd / "m.N4", "FILE4")
    d = stress_deck(NFFT, DT)
    add_eout(d, 1, [1], [1])
    write_file8(wd / "FILE8", FNUM, DF, f4["eq_node"], f4["eq_dof"], H, nfft=NFFT, delt=DT, model_hash="other")
    rc, out = run_stress(wd, d)
    assert rc == 0 and "model_hash differs from FILE4" in out
    # nodal stress frames need SOLID/SHELL elements of a 3D model (spec 05d 1.9)
    from sassi.verify.problems.vp_stress import plane_model_deck
    w2 = tmp_path_factory.mktemp("plane")
    assert run_house(w2, plane_model_deck())[0] == 0
    field_file8(w2, FNUM, DF, lambda n, x, f: deform_field(n, x, f), NFFT, DT)
    write_motion_file(w2 / "eq.acc", 0.25 * synthetic_motion(380, DT, seed=4, fmax=20.0), DT)
    (w2 / "Frames.txt").write_text("1\n10\n0\n")
    d = stress_deck(NFFT, DT, rstns=1, savemax=1)
    add_eout(d, 1, [1, 2], [2, 2, 2])
    rc, out = run_stress(w2, d)
    assert rc == 0 and "only SOLID and SHELL elements of 3D models" in out and not (w2 / "NSTRESS").exists()
    blocks = SL.read_element_center(w2 / "ELEMENT_CENTER_ABS_MAX_STRESSES.TXT")
    assert blocks[0].etype == "PLANE" and np.all(blocks[0].values[:, 3:] == 0.0)          # D-FIL-06: Sxx Szz Txz 0 0 0
    sxx = read_ths(w2, "PLANE", 1, 2, "SXX")[0]
    assert blocks[0].values[1, 0] == pytest.approx(np.max(np.abs(sxx)), rel=1e-9)


# =============================================================================================
# nodal frames: structure-only classes, SHELL membrane in global axes, BEAMS K nodes (review fixes)
# =============================================================================================
WALL_ANGLE = np.radians(30.0)


def _solid_and_rotated_wall(wd, excavate_solid=False):
    """SOLID block (group 1, unit cube) and a SHELL wall (group 2) hinged on the cube edge x = 1, y = 0 and
    turned 30 deg about Z.  Wall nodes I (1,0,0), J (1+c,s,0), K (1+c,s,1), L (1,0,1): by spec 08 4.6
    x' = (c, s, 0), z' = (s, -c, 0), y' = (0, 0, 1)."""
    d = house_deck()
    nid = {}
    n = 0
    for iz, z in enumerate((0.0, 1.0)):
        for iy, y in enumerate((0.0, 1.0)):
            for ix, x in enumerate((0.0, 1.0)):
                n += 1
                nid[(ix, iy, iz)] = n
                add_node(d, n, x, y, z)
    c, s = np.cos(WALL_ANGLE), np.sin(WALL_ANGLE)
    add_node(d, 9, 1.0 + c, s, 0.0)
    add_node(d, 10, 1.0 + c, s, 1.0)
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("groups").append([1, 1, "block"])
    d.table("groups").append([2, 3, "wall"])
    ns = [nid[(0, 0, 0)], nid[(1, 0, 0)], nid[(1, 1, 0)], nid[(0, 1, 0)],
          nid[(0, 0, 1)], nid[(1, 0, 1)], nid[(1, 1, 1)], nid[(0, 1, 1)]]
    add_element(d, 1, 1, ns, etype=2 if excavate_solid else 1, mat=1)
    add_element(d, 2, 1, [nid[(1, 0, 0)], 9, 10, nid[(1, 0, 1)]], mat=1, thick=0.3)
    rc, out = run_house(wd, d)
    assert rc == 0, out[-1500:]
    field_file8(wd, FNUM, DF, deform_field, NFFT, DT)
    write_motion_file(wd / "eq.acc", 0.25 * synthetic_motion(380, DT, seed=4, fmax=20.0), DT)
    (wd / "Frames.txt").write_text("1\n101\n0\n")
    return nid


def test_mixed_solid_shell_frames_rotate_shell_membrane_to_global_axes(tmp_path):
    """spec 05d 1.9: SOLID + SHELL -> membrane frames only.  The SHELL membrane stresses (local x'y') are
    rotated to global axes before averaging with the SOLID stresses (review finding: no mixing of
    coordinate systems).  Closed form for x' = (c, s, 0), y' = (0, 0, 1):
    SXX = FXX c^2, SYY = FXX s^2, SZZ = FYY, SXY = FXX c s, SXZ = FXY c, SYZ = FXY s."""
    nid = _solid_and_rotated_wall(tmp_path)
    d = stress_deck(NFFT, DT, rstns=1)
    add_eout(d, 1, [1], [2] * 6)
    add_eout(d, 2, [1], [2] * 6)
    rc, out = run_stress(tmp_path, d)
    assert rc == 0, out[-2000:]
    assert "global axes; membrane stresses only (SOLID + SHELL" in out
    sol = np.column_stack([read_ths(tmp_path, "SOLID", 1, 1, k)[0] for k in SL.CENTER_COLUMNS["SOLID"]])
    fxx, fyy, fxy = (read_ths(tmp_path, "SHELL", 2, 1, k)[0] for k in ("FXX", "FYY", "FXY"))
    c, s = np.cos(WALL_ANGLE), np.sin(WALL_ANGLE)
    wall = np.column_stack([fxx * c * c, fxx * s * s, fyy, fxx * c * s, fxy * c, fxy * s])
    assert np.max(np.abs(wall[:, 4])) > 0 and np.max(np.abs(wall[:, 0])) > 0          # non-trivial check
    k = 100
    for cls, cols in (("sig", slice(0, 3)), ("tau", slice(3, 6))):
        fr = np.loadtxt(tmp_path / "NSTRESS" / f"stress_{1.0:06.3f}_00101_{cls}", skiprows=1)
        rows = {int(r[0]): r[1:] for r in fr}
        scale = np.max(np.abs(np.r_[sol[k], wall[k]]))
        for node in (9, 10):                                                        # wall only
            assert np.allclose(rows[node], wall[k, cols], rtol=1e-9, atol=1e-12 * scale), (cls, node)
        for node in (nid[(1, 0, 0)], nid[(1, 0, 1)]):                               # wall and block
            assert np.allclose(rows[node], 0.5 * (sol[k, cols] + wall[k, cols]), rtol=1e-9, atol=1e-12 * scale)
        assert np.allclose(rows[nid[(0, 1, 1)]], sol[k, cols], rtol=1e-9, atol=1e-12 * scale)   # block only
    assert not (tmp_path / "NSTRESS" / "stress_ABS_MAX_bdsig").exists()


@pytest.mark.parametrize("extra", [{}, {"savemax": 1}, {"saveth": 1}, {"secdataopt": 1}])
def test_wall_on_excavated_soil_is_a_shell_only_structure(tmp_path, extra):
    """D-STR-08 / spec 05d 1.9: excavated soil is not structure, so the frame classes (bending frames, local
    axes) do not depend on the all-element options that make STRESS process the soil element too."""
    _solid_and_rotated_wall(tmp_path, excavate_solid=True)
    d = stress_deck(NFFT, DT, rstns=1, **extra)
    add_eout(d, 2, [1], [2] * 6)
    rc, out = run_stress(tmp_path, d)
    assert rc == 0, out[-2000:]
    assert "SHELL local axes" in out and "SOLID + SHELL" not in out
    fxy = read_ths(tmp_path, "SHELL", 2, 1, "FXY")[0]
    mxy = read_ths(tmp_path, "SHELL", 2, 1, "MXY")[0]
    tau = np.loadtxt(tmp_path / "NSTRESS" / f"stress_{1.0:06.3f}_00101_tau", skiprows=1)
    bdt = np.loadtxt(tmp_path / "NSTRESS" / f"stress_{1.0:06.3f}_00101_bdtau", skiprows=1)
    assert sorted(int(v) for v in tau[:, 0]) == [2, 6, 9, 10]                       # wall nodes only
    assert np.allclose(tau[:, 1], fxy[100], rtol=1e-12) and np.all(tau[:, 2:] == 0.0)
    assert np.allclose(bdt[:, 1], mxy[100], rtol=1e-12)


def test_soil_pressure_ignores_beam_orientation_nodes(tmp_path):
    """D-STR-09: the structure nodes are the DOF nodes of the structural elements.  A BEAMS K node (slot 3,
    orientation only) on the far face x = 2 of the soil element must not make that face 'shared with the
    structure': only the true contact face x = 1 gets a pressure."""
    wd = tmp_path
    nid = _two_block_model(wd)
    f4 = read_container(wd / "m.N4", "FILE4")
    d = house_deck()
    for n, p in zip(f4["node_id"], np.asarray(f4["node_xyz"]).reshape(-1, 3)):
        add_node(d, int(n), *p)
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("materials").append([2, 3, 500.0, 250.0, 19.0, 0.04, 0.04])
    d.table("beamprops").append([1, 0.25, 0.2, 0.18, 0.012, 0.004, 0.006])
    d.table("groups").append([1, 1, "wall block"])
    d.table("groups").append([2, 1, "near-field soil"])
    d.table("groups").append([3, 2, "posts"])
    for g, ix, mat in ((1, 0, 1), (2, 1, 2)):
        ns = [nid[(ix, 0, 0)], nid[(ix + 1, 0, 0)], nid[(ix + 1, 1, 0)], nid[(ix, 1, 0)],
              nid[(ix, 0, 1)], nid[(ix + 1, 0, 1)], nid[(ix + 1, 1, 1)], nid[(ix, 1, 1)]]
        add_element(d, g, 1, ns, etype=1, mat=mat)
    far = [nid[(2, iy, iz)] for iy in (0, 1) for iz in (0, 1)]
    for e, k in enumerate(far, start=1):           # posts on the structure edge x = 0, K on the soil far face
        add_element(d, 3, e, [nid[(0, 0, 0)], nid[(0, 1, 0)], k], mat=1)
    rc, out = run_house(wd, d)
    assert rc == 0, out[-1500:]
    field_file8(wd, FNUM, DF, deform_field, NFFT, DT)
    write_motion_file(wd / "eq.acc", 0.25 * synthetic_motion(380, DT, seed=4, fmax=20.0), DT)
    (wd / "Frames.txt").write_text("1\n101\n1\n2\n")
    dk = stress_deck(NFFT, DT, rstsp=1)
    add_eout(dk, 2, [1], [2] * 6)
    rc, out = run_stress(wd, dk)
    assert rc == 0, out[-2000:]
    nod = np.loadtxt(wd / "SOILPRES" / f"pres_{1.0:06.3f}_00101_nod", skiprows=1)
    assert sorted(int(v) for v in nod[:, 0]) == sorted(nid[(1, iy, iz)] for iy in (0, 1) for iz in (0, 1))
