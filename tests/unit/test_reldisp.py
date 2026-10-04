"""Unit tests of the RELDISP module (requirements section 4.9, D-RDP-01 ... D-RDP-04)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.core import signal as S
from sassi.io import textfiles
from sassi.modules.reldisp import dof_from_name
from sassi.verify.problems.vp_motion import (TWO_DOF, listing_text, motion_deck, reldisp_deck, run_deck,
                                             shear_chain_tf, synthetic_motion, write_file8, write_motion_file)

NFFT, DT, G = 1024, 0.01, 32.2
DF = 1.0 / (NFFT * DT)
FNUM = np.unique(np.concatenate([np.arange(1, 6), np.round(np.geomspace(6, 300, 40))]).astype(int))


def _motion(wd, type=0, cplx=1):
    H = np.column_stack([np.ones(len(FNUM)), shear_chain_tf(FNUM * DF, **TWO_DOF), 0.1 * shear_chain_tf(FNUM * DF, **TWO_DOF)[:, 1]])
    write_file8(wd / "FILE8", FNUM, DF, [1, 2, 3, 3], [1, 1, 1, 3], H, nfft=NFFT, delt=DT, type=type)
    acc = 0.2 * synthetic_motion(500, DT, seed=8, fmax=15.0)
    write_motion_file(wd / "eq.acc", acc, DT)
    d = motion_deck(nout=[(1, 1, 1, 1, 0, 0, 0, 0), (2, 1, 1, 1, 0, 0, 0, 0), (3, 1, 1, 1, 0, 0, 0, 0),
                          (3, 3, 1, 0, 0, 0, 0, 0)], thfile="eq.acc", nft=NFFT, delt=DT, cplx=cplx, type=type, gravity=G)
    assert run_deck(wd, "m", d) == 0
    return acc


def _rd(relfile, rdnd, **kw):
    base = dict(thfile="eq.acc", nft=NFFT, delt=DT, gravity=G)
    base.update(kw)
    return reldisp_deck(relfile, rdnd, **base)


def _tfi(wd, name):
    f, H, _ = textfiles.read_tf(wd / name)
    out = np.zeros(NFFT // 2 + 1, dtype=complex)
    out[np.rint(f / DF).astype(int)] = H
    return out


def test_dof_tags():
    assert dof_from_name("00415TR_X.TFI") == 1 and dof_from_name("C:/x/00415R_ZZ.TFI") == 6
    assert dof_from_name("123456TR_Z.TFU") == 3 and dof_from_name("ref.TFI") is None


def test_relative_to_another_node(tmp_path):
    acc = _motion(tmp_path)
    assert run_deck(tmp_path, "r", _rd("00002TR_X.TFI", [(3, 1, 0, 1, 0, 0, 0)])) == 0
    txt = listing_text(tmp_path, "r", "RELDISP")
    assert "one DOF per run" in txt and "3 TR_Z" in txt                 # Z flag ignored for a TR_X reference
    assert not (tmp_path / "00003TR_Z.THD").exists()
    d, dt = textfiles.read_history(tmp_path / "00003TR_X.THD")
    f = S.fourier_grid(NFFT, DT)
    U = S.displacement_spectrum(np.fft.rfft(S.pad_record(acc, NFFT)), f, G)
    ref = np.fft.irfft((_tfi(tmp_path, "00003TR_X.TFI") - _tfi(tmp_path, "00002TR_X.TFI")) * U, NFFT)
    assert dt == DT and np.allclose(d, ref, atol=1e-12 * np.max(np.abs(ref)))


def test_free_field_reference_per_direction(tmp_path):
    acc = _motion(tmp_path)
    assert run_deck(tmp_path, "r", _rd("", [(3, 1, 0, 1, 0, 0, 0)], dur=3.0)) == 0
    f = S.fourier_grid(NFFT, DT)
    U = S.displacement_spectrum(np.fft.rfft(S.pad_record(acc, NFFT)), f, G)
    dx, _ = textfiles.read_history(tmp_path / "00003TR_X.THD")
    dz, _ = textfiles.read_history(tmp_path / "00003TR_Z.THD")
    assert len(dx) == S.output_length(NFFT, DT, 3.0)
    rx = np.fft.irfft((_tfi(tmp_path, "00003TR_X.TFI") - 1.0) * U, NFFT)[:len(dx)]
    rz = np.fft.irfft(_tfi(tmp_path, "00003TR_Z.TFI") * U, NFFT)[:len(dz)]    # reference 0 off the input direction
    assert np.allclose(dx, rx, atol=1e-12 * np.max(np.abs(rx)))
    assert np.allclose(dz, rz, atol=1e-12 * np.max(np.abs(rz)))


def test_tfd_amplitude_only_and_frames(tmp_path):
    _motion(tmp_path)
    d = _rd("FREEFIELD", [(2, 1, 0, 0, 0, 0, 0)], reldisoutput=0, rstframes=1, dur=0.5)
    assert run_deck(tmp_path, "r", d) == 0
    f, H, cplx = textfiles.read_tf(tmp_path / "00002TR_X.TFD")
    assert not cplx and len(f) == NFFT // 2 + 1 and H[0] == 0
    frames = sorted((tmp_path / "THD").iterdir())
    assert len(frames) == S.output_length(NFFT, DT, 0.5) and frames[0].name == "THD_00.000_00001"
    assert frames[0].read_text().splitlines()[0] == "1 4"


def test_all_nodes_option(tmp_path):
    _motion(tmp_path)
    assert run_deck(tmp_path, "r", _rd("FREEFIELD", [], reldispsall=1, dur=0.2)) == 0
    for n in (1, 2, 3):
        assert (tmp_path / f"0000{n}TR_X.THD").exists()
    d1, _ = textfiles.read_history(tmp_path / "00001TR_X.THD")
    d3, _ = textfiles.read_history(tmp_path / "00003TR_X.THD")
    # the control point moves with the free field up to the SSI cut-off; only the ground content above
    # f_N (absent from MOTION's TFI, D-MOT-03) remains
    assert np.max(np.abs(d1)) < 1e-4 * np.max(np.abs(d3))
    assert (tmp_path / "THD").is_dir()


def test_missing_or_amplitude_only_tfi(tmp_path):
    _motion(tmp_path, cplx=0)
    assert run_deck(tmp_path, "a", _rd("FREEFIELD", [(3, 1, 0, 0, 0, 0, 0)])) == 1
    assert "Save Complex TF" in listing_text(tmp_path, "a", "RELDISP")
    assert run_deck(tmp_path, "b", _rd("FREEFIELD", [(7, 1, 0, 0, 0, 0, 0)])) == 1
    assert "Run MOTION with Save Complex TF for nodes 7 TR_X" in listing_text(tmp_path, "b", "RELDISP")
    assert run_deck(tmp_path, "c", _rd("00099TR_X.TFI", [(3, 1, 0, 0, 0, 0, 0)])) == 1


def test_tfi_off_grid_is_an_error(tmp_path):
    _motion(tmp_path)
    textfiles.write_tf(tmp_path / "00003TR_X.TFI", np.array([0.0, 0.037, 0.08]), np.ones(3, complex))
    assert run_deck(tmp_path, "r", _rd("FREEFIELD", [(3, 1, 0, 0, 0, 0, 0)])) == 1
    assert "not on the Fourier grid" in listing_text(tmp_path, "r", "RELDISP")


def test_vibration_relative_displacement(tmp_path):
    load = _motion(tmp_path, type=1)
    assert run_deck(tmp_path, "r", _rd("00002TR_X.TFI", [(3, 1, 0, 0, 0, 0, 0)], type=1)) == 0
    d, _ = textfiles.read_history(tmp_path / "00003TR_X.THD")
    F = np.fft.rfft(S.pad_record(load, NFFT))
    ref = np.fft.irfft(S.hermitian_bins((_tfi(tmp_path, "00003TR_X.TFI") - _tfi(tmp_path, "00002TR_X.TFI")) * F, NFFT),
                       NFFT)
    assert np.allclose(d, ref, atol=1e-12 * np.max(np.abs(ref)))
    f, H, _ = textfiles.read_tf(tmp_path / "00003TR_X.TFD")                 # vibration: TFD = H_rel
    k = np.rint(f / DF).astype(int)
    assert np.allclose(H, (_tfi(tmp_path, "00003TR_X.TFI") - _tfi(tmp_path, "00002TR_X.TFI"))[k], atol=1e-15)


def test_tfu_input_option_warns(tmp_path):
    _motion(tmp_path)
    assert run_deck(tmp_path, "r", _rd("00002TR_X.TFU", [(3, 1, 0, 0, 0, 0, 0)], reldisoutput=3)) == 0
    assert "discouraged" in listing_text(tmp_path, "r", "RELDISP")


@pytest.mark.parametrize("bad", [dict(nft=1000), dict(delt=0.0)])
def test_deck_checks(tmp_path, bad):
    _motion(tmp_path)
    assert run_deck(tmp_path, "r", _rd("FREEFIELD", [(3, 1, 0, 0, 0, 0, 0)], **bad)) == 1


def test_six_digit_node_files_with_default_maxnode(tmp_path):
    """D-FIL-05: a model node >= 100000 makes MOTION write 6-digit names (000002TR_X.TFI); RELDISP with
    the default deck <maxnode> = 0 finds them and keeps the width for .TFD/.THD."""
    H = np.column_stack([np.ones(len(FNUM)), shear_chain_tf(FNUM * DF, **TWO_DOF)])
    write_file8(tmp_path / "FILE8", FNUM, DF, [1, 2, 123456], [1, 1, 1], H, nfft=NFFT, delt=DT)
    acc = 0.2 * synthetic_motion(500, DT, seed=8, fmax=15.0)
    write_motion_file(tmp_path / "eq.acc", acc, DT)
    d = motion_deck(nout=[(2, 1, 1, 0, 0, 0, 0, 0), (123456, 1, 1, 0, 0, 0, 0, 0)], thfile="eq.acc", nft=NFFT,
                    delt=DT, cplx=1, gravity=G)
    assert run_deck(tmp_path, "m", d) == 0
    assert (tmp_path / "000002TR_X.TFI").exists() and not (tmp_path / "00002TR_X.TFI").exists()
    assert run_deck(tmp_path, "r", _rd("FREEFIELD", [(2, 1, 0, 0, 0, 0, 0), (123456, 1, 0, 0, 0, 0, 0)])) == 0
    assert (tmp_path / "000002TR_X.THD").exists() and (tmp_path / "000002TR_X.TFD").exists()
    assert (tmp_path / "123456TR_X.THD").exists()
    d2, _ = textfiles.read_history(tmp_path / "000002TR_X.THD")
    U = S.displacement_spectrum(np.fft.rfft(S.pad_record(acc, NFFT)), S.fourier_grid(NFFT, DT), G)
    f, Hi, _ = textfiles.read_tf(tmp_path / "000002TR_X.TFI")
    Hg = np.zeros(NFFT // 2 + 1, dtype=complex)
    Hg[np.rint(f / DF).astype(int)] = Hi
    ref = np.fft.irfft((Hg - 1.0) * U, NFFT)
    assert np.allclose(d2, ref, atol=1e-12 * np.max(np.abs(ref)))
    # the all-nodes mode resolves the same files
    assert run_deck(tmp_path, "ra", _rd("FREEFIELD", [], reldispsall=1)) == 0
