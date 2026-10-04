"""Regression tests of the integration fixes made while building the tutorial examples
(interpreter -> AFWRITE -> modules; see the work-package report).

* AFWRITE ``DeckBuilder.environment()``: the control direction written to the MOTION, STRESS and
  RELDISP decks for EDUOPT,TFFILE = FILE8X / FILE8Y / FILE8Z is that of the simultaneous case
  (D-ANL-06) only when ANALYS ``<simul>`` = 1; Windows paths are accepted; any other file keeps SITE
  ``<cm>`` and ANALYS ``<ang>``.
* RELDISP takes ``<cm>`` / ``<ang>`` from the solution file (deck ``<file8>``) as MOTION and STRESS do,
  so the free-field reference (D-RDP-04) is in the direction of the case even with a stale deck.
* CHECK warns (EDU-28) when an L number is repeated among the embedment layers of TOPL (manual: an L
  number may be repeated for identical layers, "but not for embedment layers").
"""
from __future__ import annotations

import numpy as np
import pytest

from sassi.io import textfiles
from sassi.prep import Interpreter, Kind
from sassi.prep.afwrite import DeckBuilder
from sassi.verify.problems.vp_motion import (listing_text, motion_deck, reldisp_deck, run_deck, synthetic_motion,
                                             write_file8, write_motion_file)


# ======================================================================================
# AFWRITE: control direction of the MOTION / STRESS / RELDISP decks
# ======================================================================================
def _ui(tmp_path, site_cm: int, ang: float, simul: int) -> Interpreter:
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(f"MDL,env,{tmp_path / 'env'}\n"
                "L,1,1.0,19.0,400,200,0.02,0.02\nTOPL,1,1\nFREQ,1,4,8\n"
                f"SITE,0,1,0,20,1,1,0,1,4096,1,{site_cm},0.005,8192,1\n"
                f"ANALYS,0,0,0,0,1,0,{ang:g},0,0,0,0,{simul}\n")
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    return ui


def _builder(ui: Interpreter, tffile: str) -> DeckBuilder:
    ui.execute(f"EDUOPT,TFFILE,{tffile}")
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    return DeckBuilder(ui.model, [ui.model.path])


@pytest.mark.parametrize("tffile, cm", [("FILE8X", 0), ("FILE8Y", 1), ("file8z", 2), ("../ex05/FILE8Y", 1),
                                        ("C:\\SSI\\M1\\FILE8Y", 1), ("FILE8", 2)])
def test_simultaneous_case_files_set_the_direction(tmp_path, tffile, cm):
    """SITE <cm> = 2 (the last SITE run was the P wave) and ANALYS <simul> = 1: FILE8Y is the y' case."""
    b = _builder(_ui(tmp_path, site_cm=2, ang=0.0, simul=1), tffile)
    assert b.environment() == (cm, 0.0)
    for deck in (b.motion(), b.stress(), b.reldisp()):
        assert (deck["cm"], deck["ang"]) == (cm, 0.0), deck.module


@pytest.mark.parametrize("tffile", ["FILE8X", "FILE8Y", "FILE8"])
def test_single_case_file_keeps_site_direction_and_angle(tmp_path, tffile):
    """ANALYS <simul> = 0: a FILE8 copied to FILE8X by the user is still the single SITE/ANALYS case
    (direction y', angle 30 degrees); the file name must not override them."""
    b = _builder(_ui(tmp_path, site_cm=1, ang=30.0, simul=0), tffile)
    assert b.environment() == (1, 30.0)
    assert (b.reldisp()["cm"], b.reldisp()["ang"]) == (1, 30.0)


# ======================================================================================
# RELDISP: direction from the solution file
# ======================================================================================
NFFT, DT, G = 1024, 0.01, 9.81
DF = 1.0 / (NFFT * DT)


def _rigid_y_case(wd):
    """FILE8Y of a simultaneous run (control direction y', cm 1): node 1 moves rigidly with the ground
    in Y at every frequency up to Nyquist; MOTION writes its complex .TFI."""
    fnum = np.arange(1, NFFT // 2 + 1)
    write_file8(wd / "FILE8Y", fnum, DF, [1], [2], np.ones((len(fnum), 1), complex), cm=1, nfft=NFFT, delt=DT)
    write_motion_file(wd / "eq.acc", 0.2 * synthetic_motion(500, DT, seed=8, fmax=15.0), DT)
    md = motion_deck(nout=[(1, 2, 1, 0, 0, 0, 0, 0)], thfile="eq.acc", nft=NFFT, delt=DT, cplx=1, gravity=G,
                     file8="FILE8Y", cm=1)
    assert run_deck(wd, "m", md) == 0


def test_reldisp_uses_the_direction_stored_in_the_solution_file(tmp_path):
    _rigid_y_case(tmp_path)
    # a stale deck: written after the P-wave SITE run (cm 2), before the FILE8Y direction was known
    rd = reldisp_deck("", [(1, 0, 1, 0, 0, 0, 0)], thfile="eq.acc", nft=NFFT, delt=DT, gravity=G, file8="FILE8Y",
                      cm=2, ang=0.0)
    assert run_deck(tmp_path, "r", rd) == 0
    d, _ = textfiles.read_history(tmp_path / "00001TR_Y.THD")
    # rigid with the free field: no relative motion (with the deck's z' direction the reference in Y
    # would be 0 and the result the absolute ground displacement, centimetres)
    assert np.max(np.abs(d)) < 1e-9
    txt = listing_text(tmp_path, "r", "RELDISP")
    assert "control direction/angle of FILE8Y (cm 1, ang 0) differ from the deck (cm 2, ang 0)" in txt


def test_reldisp_keeps_the_deck_direction_without_a_solution_file(tmp_path):
    _rigid_y_case(tmp_path)
    (tmp_path / "FILE8Y").unlink()
    rd = reldisp_deck("", [(1, 0, 1, 0, 0, 0, 0)], thfile="eq.acc", nft=NFFT, delt=DT, gravity=G, file8="FILE8Y",
                      cm=1, ang=0.0)
    assert run_deck(tmp_path, "r", rd) == 0
    d, _ = textfiles.read_history(tmp_path / "00001TR_Y.THD")
    assert np.max(np.abs(d)) < 1e-9
    assert "differ from the deck" not in listing_text(tmp_path, "r", "RELDISP")


# ======================================================================================
# CHECK: L numbers of the embedment layers
# ======================================================================================
EMBEDDED = """
L,1,1.0,19.0,500,250,0.05,0.05
L,2,1.0,19.0,500,250,0.05,0.05
L,3,1.0,20.0,800,400,0.04,0.04
TOPL,{topl}
N,1,0,0,-2
N,2,2,0,-2
N,3,0,2,-2
N,4,2,2,-2
NGEN,2,4,1,4,1,0,0,1
GROUP,1,SOLID
MACT,{bottom}
E,1,1,2,4,3,5,6,8,7
ETYPE,1,1,1,2
GROUP,2,SOLID
MACT,1
E,1,5,6,8,7,9,10,12,11
ETYPE,1,1,1,2
INT,1,12,1,1
FREQ,1,4,8
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.005,8192,1
POINT,0,2,0.9
HOUSE,9.81,0,0,2,0,0,0,0,0
AOPT,0,0,0,1,1,1,0,0,0,0,0,0,0,0
CHECK
"""


@pytest.mark.parametrize("topl, bottom, warned", [("1,1,3", 1, True), ("1,2,3", 2, False)])
def test_check_warns_on_repeated_embedment_layers(tmp_path, topl, bottom, warned):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(f"MDL,emb,{tmp_path / 'emb'}\n" + EMBEDDED.format(topl=topl, bottom=bottom))
    rep = ui.session["check_report"]
    msgs = [m for m in rep.messages if m.number == "EDU-28" and "embedment" in m.detail]
    assert bool(msgs) is warned, rep.format()
    if warned:
        assert msgs[0].kind == "Warning" and "L 1 repeated" in msgs[0].detail
    assert not [m for m in rep.messages if m.number == "EDU-08"], rep.format()      # MSET = TOPL layer
