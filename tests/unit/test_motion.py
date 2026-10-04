"""Unit tests of the MOTION module (requirements section 4.8; UT-07 MOTION scaling)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sassi.conventions import nodal_rs_name, nodal_result_name
from sassi.core import signal as S
from sassi.core.spectra import log_frequencies, response_spectrum
from sassi.io import textfiles
from sassi.modules.base import ModuleError
from sassi.modules.motion import read_control_motion, scale_control_motion
from sassi.verify.problems.vp_motion import (TWO_DOF, listing_text, motion_deck, run_deck, shear_chain_tf,
                                             synthetic_motion, write_file8, write_motion_file)

NFFT, DT = 1024, 0.01
DF = 1.0 / (NFFT * DT)
FNUM = np.unique(np.round(np.geomspace(2, 300, 40)).astype(int))


def _setup(wd: Path, n=600, eq_node=(1, 2, 3), eq_dof=(1, 1, 1), ang=0.0, type=0, H=None, fnum=FNUM):
    if H is None:
        H = np.column_stack([np.ones(len(fnum)), shear_chain_tf(fnum * DF, **TWO_DOF)])
    write_file8(wd / "FILE8", fnum, DF, eq_node, eq_dof, H, nfft=NFFT, delt=DT, ang=ang, type=type)
    acc = 0.2 * synthetic_motion(n, DT, seed=4, fmax=20.0)
    write_motion_file(wd / "eq.acc", acc, DT)
    return acc


def _deck(**kw):
    base = dict(thfile="eq.acc", nft=NFFT, delt=DT)
    base.update(kw)
    return motion_deck(**base)


# ------------------------------------------------------------------ UT-07 scaling
def test_ut07_scaling_rules():
    a = np.array([0.1, -0.5, 0.25, 0.0, 0.4])
    assert np.array_equal(scale_control_motion(a, 2.5, 0.0), a * 2.5)
    s = scale_control_motion(a, 0.0, 0.3)
    assert np.array_equal(s, a * (0.3 / 0.5))                      # a * max / max|a|
    assert np.max(np.abs(s)) == pytest.approx(0.3, rel=1e-15)
    with pytest.raises(ModuleError, match="Error 77"):
        scale_control_motion(a, 0.0, 0.0)
    with pytest.raises(ModuleError, match="Error 78"):
        scale_control_motion(a, 1.0, 0.3)


@pytest.mark.parametrize("mult,mx,err", [(0.0, 0.0, "Error 77"), (1.0, 0.3, "Error 78")])
def test_ut07_module_reports_scaling_errors(tmp_path, mult, mx, err):
    _setup(tmp_path)
    d = _deck(nout=[(3, 1, 1, 1, 0, 0, 0, 1)])
    d["mult"], d["max"] = mult, mx
    assert run_deck(tmp_path, "m", d) == 1
    assert err in listing_text(tmp_path, "m", "MOTION")


def test_ut07_module_scales_to_max(tmp_path):
    acc = _setup(tmp_path, eq_node=(1,), eq_dof=(1,), H=np.ones((NFFT // 2, 1)), fnum=np.arange(1, NFFT // 2 + 1))
    d = _deck(nout=[(1, 1, 0, 1, 0, 0, 0, 1)], interp=1)
    d["mult"], d["max"] = 0.0, 0.35
    assert run_deck(tmp_path, "m", d) == 0
    out, _ = textfiles.read_history(tmp_path / nodal_result_name(1, 1, "ACC"))
    assert np.allclose(out[:len(acc)], acc * 0.35 / np.max(np.abs(acc)), atol=1e-14)


# ------------------------------------------------------------------ control motion reading
def test_record_selection_formats_and_checks(tmp_path):
    a = np.arange(1, 11) * 0.01
    write_motion_file(tmp_path / "a0.acc", a, DT)
    cm = read_control_motion(tmp_path / "a0.acc", 0, 3, 7, DT, 64, 1.0, 0.0)
    assert np.allclose(cm.acc, a[2:7]) and cm.nrec == 10 and len(cm.padded) == 64
    (tmp_path / "a1.acc").write_text("\n".join(f"{i * DT:.4f} {v}" for i, v in enumerate(a)))
    cm1 = read_control_motion(tmp_path / "a1.acc", 1, 1, 0, DT, 64, 2.0, 0.0)
    assert np.allclose(cm1.acc, 2 * a)
    with pytest.raises(ModuleError, match="time step"):
        read_control_motion(tmp_path / "a0.acc", 0, 1, 0, 0.02, 64, 1.0, 0.0)
    with pytest.raises(ModuleError, match="EDU-03"):
        read_control_motion(tmp_path / "a0.acc", 0, 1, 0, DT, 8, 1.0, 0.0)
    with pytest.raises(ModuleError, match="Error 74"):
        read_control_motion(tmp_path / "a0.acc", 0, 11, 0, DT, 64, 1.0, 0.0)
    with pytest.raises(ModuleError, match="Error 75"):
        read_control_motion(tmp_path / "a0.acc", 0, 1, 12, DT, 64, 1.0, 0.0)
    with pytest.raises(ModuleError, match="Error 76"):
        read_control_motion(tmp_path / "a0.acc", 0, 6, 4, DT, 64, 1.0, 0.0)
    with pytest.raises(ModuleError, match="Error 73"):
        read_control_motion(tmp_path / "none.acc", 0, 1, 0, DT, 64, 1.0, 0.0)


# ------------------------------------------------------------------ consistency errors
def test_df_type_and_file_errors(tmp_path):
    _setup(tmp_path)
    req = [(3, 1, 1, 1, 0, 0, 0, 1)]
    write_motion_file(tmp_path / "eq005.acc", np.ones(10), 0.005)
    assert run_deck(tmp_path, "a", _deck(nout=req, delt=0.005, thfile="eq005.acc")) == 1  # df mismatch
    assert "D-MOT-14" in listing_text(tmp_path, "a", "MOTION")
    assert run_deck(tmp_path, "b", _deck(nout=req, type=1)) == 1                         # type mismatch
    assert "does not match" in listing_text(tmp_path, "b", "MOTION")
    assert run_deck(tmp_path, "c", _deck(nout=req, file8="FILE8X")) == 1                 # missing FILE8X
    assert "FILE8X missing -- run ANALYS" in listing_text(tmp_path, "c", "MOTION")
    assert run_deck(tmp_path, "e", _deck()) == 1                                         # no request
    assert "Error 64" in listing_text(tmp_path, "e", "MOTION")
    assert run_deck(tmp_path, "f", _deck(nout=req, nft=1000)) == 1                       # NFFT not 2^n
    assert run_deck(tmp_path, "g", _deck(nout=req, interp=7)) == 1


def test_data_check_and_tf_only(tmp_path):
    _setup(tmp_path)
    req = [(3, 1, 1, 1, 0, 0, 1, 1)]
    assert run_deck(tmp_path, "dc", _deck(nout=req, opmode=1, damp=[0.05])) == 0
    assert not list(tmp_path.glob("*.TF?")) and not list(tmp_path.glob("*.ACC"))
    assert run_deck(tmp_path, "tf", _deck(nout=req, out=1, damp=[0.05])) == 0
    assert (tmp_path / "00003TR_X.TFU").exists() and (tmp_path / "00003TR_X.TFI").exists()
    assert not list(tmp_path.glob("*.ACC")) and not list(tmp_path.glob("*.RS"))


def test_amplitude_only_tf_files(tmp_path):
    _setup(tmp_path)
    assert run_deck(tmp_path, "m", _deck(nout=[(3, 1, 1, 0, 0, 0, 0, 0)], cplx=0)) == 0
    f, H, cplx = textfiles.read_tf(tmp_path / "00003TR_X.TFI")
    assert not cplx and np.all(H >= 0)
    assert f[0] == 0.0 and f[-1] <= FNUM[-1] * DF + 1e-9 and f[-1] > FNUM[-1] * DF - DF


def test_requests_merge_and_fixed_dofs(tmp_path):
    _setup(tmp_path)
    nout = [(3, 1, 1, 0, 0, 0, 0, 0), (3, 1, 0, 1, 0, 0, 0, 0), (3, 2, 1, 0, 0, 0, 0, 0), (99, 1, 1, 0, 0, 0, 0, 0)]
    assert run_deck(tmp_path, "m", _deck(nout=nout)) == 0
    txt = listing_text(tmp_path, "m", "MOTION")
    assert "EDU-16" in txt and "constrained (fixed) DOFs ignored: 3 TR_Y" in txt and "Error 65" in txt
    assert (tmp_path / "00003TR_X.TFU").exists() and (tmp_path / "00003TR_X.ACC").exists()


def test_rs_file_equals_response_spectrum_of_the_history(tmp_path):
    _setup(tmp_path)
    d = _deck(nout=[(3, 1, 0, 1, 0, 0, 1, 1)], damp=[0.02, 0.05], freq1=0.5, freq2=50.0, fstep=31)
    assert run_deck(tmp_path, "m", d) == 0
    acc, _ = textfiles.read_history(tmp_path / "00003TR_X.ACC")
    rs = response_spectrum(acc, DT, log_frequencies(0.5, 50.0, 31), [0.02, 0.05])
    for i in (1, 2):
        xy = textfiles.read_xy(tmp_path / nodal_rs_name(3, 1, i))
        assert np.allclose(xy[:, 0], log_frequencies(0.5, 50.0, 31), rtol=1e-9)
        assert np.allclose(xy[:, 1], rs["SA"][i - 1], rtol=1e-12)
    assert "Maximum requested response" in listing_text(tmp_path, "m", "MOTION")


def test_rotated_input_zero_frequency_anchor(tmp_path):
    fn = FNUM
    H = np.ones((len(fn), 2), dtype=complex)
    _setup(tmp_path, eq_node=(1, 1), eq_dof=(1, 2), ang=30.0, H=H * [np.cos(np.pi / 6), np.sin(np.pi / 6)])
    d = _deck(nout=[(1, 1, 1, 0, 0, 0, 0, 0), (1, 2, 1, 0, 0, 0, 0, 0)], ang=30.0)
    assert run_deck(tmp_path, "m", d) == 0
    _, hx, _ = textfiles.read_tf(tmp_path / "00001TR_X.TFI")
    _, hy, _ = textfiles.read_tf(tmp_path / "00001TR_Y.TFI")
    assert np.isclose(hx[0], np.cos(np.pi / 6), rtol=1e-14) and np.isclose(hy[0], 0.5, rtol=1e-14)


# ------------------------------------------------------------------ vibration
def test_vibration_response_types(tmp_path):
    fn = np.arange(1, NFFT // 2 + 1)
    c = 2e-3 * (1 - 0.1j)                                         # constant compliance, displacement per unit load
    write_file8(tmp_path / "FILE8", fn, DF, [5], [3], np.full((len(fn), 1), c), nfft=NFFT, delt=DT, type=1)
    load = synthetic_motion(500, DT, seed=9, fmax=10.0)
    write_motion_file(tmp_path / "load.acc", load, DT)
    F = np.fft.rfft(S.pad_record(load, NFFT))
    w = S.omega_grid(NFFT, DT)
    Hg = np.full(len(F), c)
    Hg[-1] = Hg[-1].real
    expect = {0: Hg * F, 1: 1j * w * Hg * F, 2: -w ** 2 * Hg * F}
    for resp in (0, 1, 2):
        sub = tmp_path / f"r{resp}"
        sub.mkdir()
        for nm in ("FILE8", "load.acc"):
            (sub / nm).write_bytes((tmp_path / nm).read_bytes())
        d = motion_deck(nout=[(5, 3, 1, 1, 0, 0, 0, 1)], thfile="load.acc", nft=NFFT, delt=DT, type=1, resp=resp)
        assert run_deck(sub, "v", d) == 0
        out, _ = textfiles.read_history(sub / "00005TR_Z.ACC")
        ref = np.fft.irfft(S.hermitian_bins(expect[resp], NFFT), NFFT)
        assert np.allclose(out, ref, atol=1e-12 * np.max(np.abs(ref)))


# ------------------------------------------------------------------ auxiliary outputs
def test_conttrs_conversion(tmp_path):
    a = 0.2 * synthetic_motion(400, DT, seed=6)
    textfiles.write_history(tmp_path / "ext.ACC", a, DT)
    (tmp_path / "CONTTRS.txt").write_text("ext.ACC\nbad.txt\n")
    d = _deck(cnvrt=1, damp=[0.05], freq1=0.5, freq2=50.0, fstep=21)
    assert run_deck(tmp_path, "m", d) == 0
    xy = textfiles.read_xy(tmp_path / "ext.RSO")
    a8, _ = textfiles.read_history(tmp_path / "ext.ACC")
    rs = response_spectrum(S.pad_record(a8, NFFT), DT, log_frequencies(0.5, 50.0, 21), [0.05])
    assert np.allclose(xy[:, 1], rs["SA"][0], rtol=1e-12)
    assert ".ACC extension" in listing_text(tmp_path, "m", "MOTION")


def test_all_points_frames_and_file12(tmp_path):
    fn = FNUM
    Hc = shear_chain_tf(fn * DF, **TWO_DOF)
    eq_node = [1, 1, 2, 2, 2]
    eq_dof = [1, 2, 1, 2, 5]
    H = np.column_stack([np.ones(len(fn)), np.zeros(len(fn)), Hc[:, 1], 0.01 * Hc[:, 0], 1e-3 * Hc[:, 1]])
    _setup(tmp_path, eq_node=eq_node, eq_dof=eq_dof, H=H)
    d = _deck(damp=[0.05], freq1=0.5, freq2=20.0, fstep=5, dur=2.0, savetf=1, saveacc=1, savers=1, saverot=1,
              rsttf=1, rstacc=1, rstrs=1, f1213=2)
    assert run_deck(tmp_path, "m", d) == 0
    for n, tag in ((1, "TR_X"), (1, "TR_Y"), (2, "TR_X"), (2, "TR_Y")):
        for ext in ("TFU", "TFI", "ACC"):
            assert (tmp_path / f"0000{n}{tag}.{ext}").exists()
        assert (tmp_path / f"0000{n}{tag}01.RS").exists()
    assert (tmp_path / "00002R_YY.ACC").exists() and not (tmp_path / "00002R_YY.TFU").exists()
    tfu = sorted((tmp_path / "TFU").iterdir())
    assert len(tfu) == len(fn) and tfu[0].name == f"TFU_{fn[0] * DF:06.2f}_00001"
    head = tfu[0].read_text().splitlines()
    assert head[0] == "2 7" and head[1].split()[0] == "1"
    nacc = S.output_length(NFFT, DT, 2.0)
    assert len(list((tmp_path / "ACC").iterdir())) == nacc and len(list((tmp_path / "ACCR").iterdir())) == nacc
    assert (tmp_path / "ACC" / "ACC_00.000_00001").exists()
    assert len(list((tmp_path / "RS").iterdir())) == 5
    lines = (tmp_path / "ACC_max.txt").read_text().splitlines()
    assert lines[0] == "2 4"                                     # frame header: nrows ncols (D-FIL-03)
    zpa = np.array([[float(v) for v in ln.split()] for ln in lines[1:]])
    ax, _ = textfiles.read_history(tmp_path / "00002TR_X.ACC")
    assert zpa[1, 0] == 2 and np.isclose(zpa[1, 1], np.max(np.abs(ax)), rtol=1e-9)
    txt = (tmp_path / "FILE12").read_text()
    assert "node 2 TR_X acceleration history" in txt and "RS damping 0.05" in txt


def test_srss_input_errors(tmp_path):
    _setup(tmp_path)
    req = [(3, 1, 1, 0, 0, 0, 0, 0)]
    assert run_deck(tmp_path, "a", _deck(nout=req, srss=1)) == 1
    assert "SRSSTF.txt is missing" in listing_text(tmp_path, "a", "MOTION")
    write_file8(tmp_path / "FILE8_02", FNUM[:-1], DF, [1, 2, 3], [1, 1, 1], np.ones((len(FNUM) - 1, 3)),
                nfft=NFFT, delt=DT)
    (tmp_path / "srsstf.TXT").write_text("2 0\nFILE8\nFILE8_02\n")          # case-insensitive lookup
    assert run_deck(tmp_path, "b", _deck(nout=req, srss=1)) == 1
    assert "different frequency list" in listing_text(tmp_path, "b", "MOTION")


def test_low_frequency_atf_check_and_rs_listing(tmp_path):
    H = np.column_stack([np.ones(len(FNUM)), shear_chain_tf(FNUM * DF, **TWO_DOF)])
    H[0, 2] *= 1.2                                                # |ATF(f_1)| 20 % above the rigid-body value
    _setup(tmp_path, H=H)
    d = _deck(nout=[(3, 1, 0, 1, 0, 0, 1, 1), (2, 1, 0, 0, 0, 0, 0, 1)], damp=[0.05], freq1=0.5, freq2=40.0,
              fstep=11, gravity=9.81)
    assert run_deck(tmp_path, "m", d) == 0
    txt = listing_text(tmp_path, "m", "MOTION")
    assert "EDU-18" in txt and "3 TR_X (1.2" in txt and "2 TR_X (" not in txt.split("EDU-18")[1].splitlines()[0]
    assert "Response spectra node 3 TR_X" in txt and "SV relative velocity in length/s" in txt
    acc, _ = textfiles.read_history(tmp_path / "00003TR_X.ACC")
    rs = response_spectrum(acc, DT, log_frequencies(0.5, 40.0, 11), [0.05])
    assert f"peak SV {rs['SV'][0].max() * 9.81:.6e}" in txt


# ------------------------------------------------------------------ review fixes
@pytest.mark.parametrize("freq1,freq2,fstep,err", [(0.0, 50.0, 31, "Error 69"), (10.0, 5.0, 31, "Error 70"),
                                                    (10.0, 10.0, 31, "Error 70"), (0.1, 50.0, 1, "Error 71"),
                                                    (0.1, 50.0, 0, "Error 71"), (-1.0, 50.0, 31, "Error 69")])
def test_rs_frequency_data_validated_when_rs_requested(tmp_path, freq1, freq2, fstep, err):
    """Spec 05c B.8 (impl. decision): RS output needs 0 < freq1 < freq2 and fstep >= 2."""
    _setup(tmp_path)
    rs_req = dict(freq1=freq1, freq2=freq2, fstep=fstep, damp=[0.05])
    assert run_deck(tmp_path, "a", _deck(nout=[(3, 1, 0, 1, 0, 0, 1, 1)], **rs_req)) == 1       # NOUT flag 5
    assert err in listing_text(tmp_path, "a", "MOTION")
    assert run_deck(tmp_path, "b", _deck(nout=[(3, 1, 0, 1, 0, 0, 0, 1)], savers=1, **rs_req)) == 1
    assert run_deck(tmp_path, "c", _deck(cnvrt=1, **rs_req)) == 1                                # CONTTRS conversion
    assert not list(tmp_path.glob("*.RS"))


def test_rs_frequency_data_not_required_without_rs(tmp_path):
    """Without RS output (or with <out> = 1) the RS frequency fields are not used and not checked."""
    _setup(tmp_path)
    bad = dict(freq1=10.0, freq2=5.0, fstep=1)
    assert run_deck(tmp_path, "a", _deck(nout=[(3, 1, 0, 1, 0, 0, 0, 1)], **bad)) == 0
    assert run_deck(tmp_path, "b", _deck(nout=[(3, 1, 0, 1, 0, 0, 1, 1)], out=1, damp=[0.05], **bad)) == 0


def test_tf_only_output_writes_tfs_for_every_request(tmp_path):
    """B.5.3 <out> = 1: the TFs of the selected nodes are extracted whatever NOUT flags selected them."""
    _setup(tmp_path)
    assert run_deck(tmp_path, "m", _deck(nout=[(3, 1, 0, 1, 0, 0, 1, 1), (2, 1, 0, 0, 0, 0, 0, 1)], out=1)) == 0
    for n in (2, 3):
        assert (tmp_path / f"0000{n}TR_X.TFU").exists() and (tmp_path / f"0000{n}TR_X.TFI").exists()
    assert not list(tmp_path.glob("*.ACC")) and not list(tmp_path.glob("*.RS"))
    assert "<out> = 1 (transfer functions only)" in listing_text(tmp_path, "m", "MOTION")
    f, H, _ = textfiles.read_tf(tmp_path / "00003TR_X.TFU")
    assert np.allclose(H, shear_chain_tf(FNUM * DF, **TWO_DOF)[:, 1], rtol=1e-14, atol=0)


def test_phase_adjustment_fallback_warns_and_zeroes_the_phase(tmp_path):
    """N = 4 SSI frequencies: option 2 runs the spline; rho follows the spline (0), with a warning."""
    fn = np.array([5, 20, 60, 150])
    H = np.column_stack([np.ones(len(fn)), shear_chain_tf(fn * DF, **TWO_DOF)])
    _setup(tmp_path, H=H, fnum=fn)
    assert run_deck(tmp_path, "m", _deck(nout=[(3, 1, 1, 0, 0, 0, 0, 0)], interp=2, smo=10.0, pzadj=1, cplx=1)) == 0
    txt = listing_text(tmp_path, "m", "MOTION")
    assert "replaced by option 6; smoothing is not applied" in txt
    assert "phase adjustment follows the scheme actually used: rho = 0" in txt
    _, Hi, _ = textfiles.read_tf(tmp_path / "00003TR_X.TFI")
    assert np.max(np.abs(np.angle(Hi))) == 0.0


def test_phase_adjustment_keeps_a_negative_rigid_body_anchor(tmp_path):
    """ang = 180 deg: H = -1 (rigid body, X).  With phase adjustment (option 6, rho = 0) the response
    stays the reversed control motion; the old absolute-phase rule turned it into +input."""
    fn = np.arange(1, NFFT // 2 + 1)
    acc = _setup(tmp_path, eq_node=(4,), eq_dof=(1,), ang=180.0, H=-np.ones((len(fn), 1)), fnum=fn)
    d = _deck(nout=[(4, 1, 1, 1, 0, 0, 0, 1)], ang=180.0, interp=6, pzadj=1, cplx=1)
    assert run_deck(tmp_path, "m", d) == 0
    _, Hi, _ = textfiles.read_tf(tmp_path / "00004TR_X.TFI")
    assert np.allclose(Hi, -1.0, rtol=0, atol=1e-14)
    a, _ = textfiles.read_history(tmp_path / "00004TR_X.ACC")
    assert np.max(np.abs(a + S.pad_record(acc, NFFT))) < 1e-12
