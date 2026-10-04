"""Interpreter commands of the non-linear soil SSI (sassi/prep/commands/nlsoil_cmds.py): the .pin records
PIN / PINGRP / PINMAT / PINDEL / PINLIST, AFWRITE of <model>.pin, the CHECK rules EDU-41/42/43 (and the
Error 79 / EXCSTRCHK adjustments), COMBXYZSTRAIN (UT-20), NLSSIRESET and the NLSSIITER driver, and the
WRITE -> INP round trip (requirements 3.4.R, 4.4 item 7, D-NLS-01/02/04)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.core import nlsoil as NL
from sassi.core import shake as SH
from sassi.prep import Interpreter
from sassi.prep import registry
from sassi.prep.messages import Kind
from sassi.prep.writer import write_pre
from sassi.io.files import write_container
from pathlib import Path

#: a one-element column of excavated soil (groups 1, 2) and near-field soil (group 3, two materials)
COLUMN = """MDL,m,m
TIT,nonlinear column
L,1,1.0,19.62,400,200,0.02,0.02
L,2,1.0,19.62,400,200,0.02,0.02
L,3,1.0,21.58,2000,1000,0.01,0.01
TOPL,1,2
FREQ,1,4,20,41
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.01,1024,1
WAVE,2,1,1,1,0
HOUSE,9.81,0,0,2,0,0,0,0,0
DYNP,1,0.0001,1.0,0.0001,0.24,Clay
DYNP,2,1.0,0.5,1.0,10.0,Clay
DYNP,1,0.0001,1.0,0.0001,0.24,Sand
DYNP,2,1.0,0.1,1.0,20.0,Sand
SPRO,1,1,Sand
SPRO,2,2,Sand
SPRO,3,3
N,1,0,0,-2
N,2,1,0,-2
N,3,1,1,-2
N,4,0,1,-2
NGEN,2,4,1,4,1,0,0,1
GROUP,1,SOLID
MACT,1
E,1,5,6,7,8,9,10,11,12
ETYPE,1,1,1,2
GROUP,2,SOLID
MACT,2
E,1,1,2,3,4,5,6,7,8
ETYPE,1,1,1,2
M,1,400,200,19.62,0.02,0.02,3
M,2,500,250,19.62,0.02,0.02,3
GROUP,3,SOLID
MACT,1
E,1,1,2,3,4,5,6,7,8
MACT,2
E,2,5,6,7,8,9,10,11,12
ETYPE,1,2,1,1
INT,1,12,1,1
POINT,0,2,0.9
"""


def _ui(tmp_path, extra=""):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(COLUMN + extra)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    return ui


def _checks(ui):
    """{label: [details]} of the CHECK messages of the last CHECK/AFWRITE (from the .err file)."""
    return (ui.model.path and (lambda p: p.read_text() if p.exists() else "")(
        __import__("pathlib").Path(ui.model.path) / f"{ui.model.name}.err")) or ""


# ======================================================================================
# .pin records
# ======================================================================================
def test_pin_records_and_validation(tmp_path):
    ui = _ui(tmp_path)
    opts = ui.model.options
    assert ui.execute("PIN,0.6,2") and opts.record("PIN").to_tokens() == ["0.6", "2"]
    assert not ui.execute("PIN,1.5")                                   # ESF outside (0, 1]
    assert not ui.execute("PIN,0.6,0")                                 # NCURV >= 1
    assert ui.execute("PINGRP,3,1,0.8,1.5,Sand")
    assert opts.entry("PINGRP", 3).to_tokens() == ["3", "1", "0.8", "1.5", "Sand"]
    assert not ui.execute("PINGRP,3,2")                                # ISTR 0/1
    assert not ui.execute("PINGRP,0,0")
    assert not ui.execute("PINGRP,3,0,0")                              # GFAC > 0
    assert not ui.execute("PINGRP,3,0,1,-1")                           # DFAC >= 0
    assert not ui.execute("PINGRP,3,0.5")                              # integer ISTR
    assert ui.execute("PINGRP,1,0")                                    # stored, but a warning (excavated SOLID ok here)
    assert ui.execute("PINGRP,9,0") and any("not defined" in t for t in ui.sink.texts(Kind.WARNING))
    assert ui.execute("PINMAT,3,2,0.5,2,Clay")
    assert opts.entry("PINMAT", (3, 2)).to_tokens() == ["3", "2", "0.5", "2", "Clay"]
    assert ui.execute("PINMAT,7,1,1,1,1") and any("no PINGRP" in t for t in ui.sink.texts(Kind.WARNING))
    assert not ui.execute("PINMAT,3,0,1,1,1")
    # PINDEL of one group removes its PINMAT lines, blank removes everything
    assert ui.execute("PINDEL,3")
    assert opts.entry("PINGRP", 3) is None and opts.entry("PINMAT", (3, 2)) is None
    assert opts.entry("PINMAT", (7, 1)) is not None
    assert ui.execute("PINDEL")
    assert not opts.entries("PINGRP") and not opts.entries("PINMAT") and opts.record("PIN") is None


def test_pin_from_model_curve_labels_and_defaults(tmp_path):
    from sassi.prep.commands.nlsoil_cmds import file73_labels, pin_from_model, pin_text
    ui = _ui(tmp_path, "PINGRP,3,1,1.0,1.0,Sand\nPINMAT,3,2,0.5,2.0,Clay\n")
    m = ui.model
    assert file73_labels(m) == ["Sand", "Clay"]                        # SPRO labels first, then the others
    pin, probs, notes, mnotes = pin_from_model(m)
    assert not probs
    assert (pin.esf, pin.ncurv) == (0.65, 2) and any("SOIL <ratio>" in t for t in notes)     # PIN not given
    g = pin.group(3)
    assert (g.nmat, g.nelem, g.istr) == (2, 2, 1)
    assert [(x.gfac, x.dfac, x.icurve) for x in g.materials] == [(1.0, 1.0, 1), (0.5, 2.0, 2)]
    assert mnotes[3] == ["material 1, curve Sand", "material 2, curve Clay"]
    text, _, _ = pin_text(m)
    assert NL.parse_pin(text) == pin
    ui.execute("PIN,0.7,5")
    ui.execute("PINMAT,3,2,0.5,2.0,2")                                 # a curve number is kept as given
    pin, probs, _, _ = pin_from_model(m)
    assert not probs and (pin.esf, pin.ncurv) == (0.7, 5) and pin.group(3).materials[1].icurve == 2


@pytest.mark.parametrize("extra,msg", [
    ("PINGRP,3,0\n", "no ICURVE"),
    ("PINGRP,3,0,1,1,Gravel\n", "not a DYNP label"),
    ("PINGRP,3,0,1,1,7\nPIN,0.6,2\n", "ICURVE 7 > NCURV 2"),
    ("PINGRP,1,0,1,1,Sand\n", "excavated soil elements"),
    ("PINGRP,8,0,1,1,Sand\n", "not defined"),
    # review findings: ETYPE 0 below grade is excavated soil (D-HOU-01); the curves of the .pin are checked
    # like SOIL's Errors 97-99 even when no SPRO sublayer uses them
    ("ETYPE,1,2,1,0\nPINGRP,3,0,1,1,Sand\n", "ETYPE 0 below grade"),
    ("DYNP,1,0.0001,1.0,,,Back\nDYNP,2,1.0,0.3,,,Back\nPINGRP,3,0,1,1,Back\n", "curve Back (ICURVE 2): no damping"),
    ("DYNP,1,0.0001,1.0,0.0,5.0,Back\nPINGRP,3,0,1,1,Back\n", "no damping point with a positive strain"),
    ("DYNP,1,,,0.0001,5.0,Back\nPINGRP,3,0,1,1,Back\n", "no G/G_max point"),
    ("DYNP,1,0.0001,1.5,0.0001,5.0,Back\nPINGRP,3,0,1,1,Back\n", "outside [0, 1]"),
    ("DYNP,1,,,,,Void\nPINGRP,3,0,1,1,Void\n", "no complete G/G_max or damping point"),
])
def test_pin_from_model_problems(tmp_path, extra, msg):
    from sassi.prep.commands.nlsoil_cmds import pin_from_model
    ui = _ui(tmp_path, extra)
    _, probs, _, _ = pin_from_model(ui.model)
    assert any(msg in t for _, t in probs), probs


def test_file73_labels_leave_out_labels_without_points(tmp_path):
    """AFWRITE writes no SOIL deck row for a DYNP label without any complete point, so SOIL does not number
    it in FILE73: the label -> ICURVE conversion must skip it too."""
    from sassi.prep.commands.nlsoil_cmds import file73_labels, pin_from_model
    ui = _ui(tmp_path, "DYNP,1,,,,,Void\nDYNP,1,0.0001,1.0,0.0001,1.0,Back\nDYNP,2,1.0,0.3,1.0,12.0,Back\n"
                       "PINGRP,3,0,1,1,Back\n")
    assert file73_labels(ui.model) == ["Sand", "Back", "Clay"]        # SPRO labels, then the others in key order
    pin, probs, _, _ = pin_from_model(ui.model)
    assert not probs and pin.group(3).materials[0].icurve == 2


def test_check_reports_the_resolved_etype_and_the_curves(tmp_path):
    """CHECK (EDU-41) and AFWRITE: near-field SOLIDs left at ETYPE 0 below grade, and a curve of the .pin
    without damping points, block the HOUSE deck instead of failing in HOUSE / STRESS."""
    ui = _ui(tmp_path, "ETYPE,1,2,1,0\nPINGRP,3,0,1,1,Sand\nHOUSEX,0,0,1,1\nSTRESS,0,1,0,0,1\n"
                       "AOPT,0,0,0,1,0,1,0,0,0,0,0,0,0,0\n")
    ui.execute("AFWRITE")
    err = _checks(ui)
    assert "EDU-41" in err and "ETYPE 0 below grade" in err and "set ETYPE 1" in err
    assert not (tmp_path / "m" / "m.hou").exists() and not (tmp_path / "m" / "m.pin").exists()
    ui.execute("ETYPE,1,2,1,1")
    ui.execute("DYNP,1,0.0001,1.0,,,Back")
    ui.execute("PINGRP,3,0,1,1,Back")
    ui.execute("AFWRITE")
    err = _checks(ui)
    assert "EDU-41" in err and "curve Back" in err and "damping" in err
    ui.execute("DYNP,1,0.0001,1.0,0.0001,2.0,Back")
    ui.execute("AFWRITE")
    assert "EDU-41" not in _checks(ui) and (tmp_path / "m" / "m.pin").exists()


def test_write_inp_round_trip_of_the_pin_records(tmp_path):
    ui = _ui(tmp_path, "PIN,0.6\nPINGRP,3,1,1.0,1.0,Sand\nPINMAT,3,2,0.5,2.0,Clay\nHOUSEX,0,0,1,1\nSTRESS,0,1,0,0,1\n")
    text, _ = write_pre(ui.model)
    assert "PIN,0.6" in text and "PINGRP,3,1,1.0,1.0,Sand" in text and "PINMAT,3,2,0.5,2.0,Clay" in text
    b = Interpreter(cwd=tmp_path)
    b.run_text(text)
    assert not b.sink.texts(Kind.ERROR)
    assert b.model.same_state(ui.model)
    body = [ln for ln in text.splitlines() if not ln.startswith("*")]
    assert [ln for ln in write_pre(b.model)[0].splitlines() if not ln.startswith("*")] == body


# ======================================================================================
# AFWRITE and CHECK
# ======================================================================================
def test_afwrite_writes_the_pin_file(tmp_path):
    ui = _ui(tmp_path, "PIN,0.6\nPINGRP,3,1,1.0,1.0,Sand\nPINMAT,3,2,0.5,2.0,Clay\nHOUSEX,0,0,1,1\n"
                       "STRESS,0,1,0,0,1\nAOPT,0,0,0,1,0,1,0,0,0,0,0,0,0,0\n")
    assert ui.execute("AFWRITE")
    p = tmp_path / "m" / "m.pin"
    pin = NL.read_pin(p)
    assert (pin.ngrp, pin.esf, pin.ncurv) == (1, 0.6, 2)
    assert [(x.gfac, x.dfac, x.icurve) for x in pin.group(3).materials] == [(1.0, 1.0, 1), (0.5, 2.0, 2)]
    first = p.read_text()
    assert ui.execute("AFWRITE") and p.read_text() == first and not (tmp_path / "m" / "m.pin.bak").exists()
    ui.execute("PINMAT,3,2,0.4,2.0,Clay")
    assert ui.execute("AFWRITE")
    assert (tmp_path / "m" / "m.pin.bak").read_text() == first
    assert NL.read_pin(p).group(3).materials[1].gfac == 0.4
    # the .liq iteration state is never touched by AFWRITE
    assert not (tmp_path / "m" / "m.liq").exists()


def test_afwrite_keeps_a_hand_edited_pin(tmp_path):
    ui = _ui(tmp_path, "HOUSEX,0,0,1,1\nSTRESS,0,1,0,0,1\nAOPT,0,0,0,1,0,1,0,0,0,0,0,0,0,0\n")
    (tmp_path / "m").mkdir(exist_ok=True)
    (tmp_path / "m" / "m.pin").write_text("1, 0.6, 1\n3, 2, 2, 0\n1, 1, 1\n1, 1, 1\n")
    assert ui.execute("AFWRITE")
    assert (tmp_path / "m" / "m.pin").read_text() == "1, 0.6, 1\n3, 2, 2, 0\n1, 1, 1\n1, 1, 1\n"
    assert "EDU-42" not in _checks(ui)


def test_check_rules_of_the_nonlinear_soil(tmp_path):
    ui = _ui(tmp_path, "HOUSEX,0,0,1,1\nAOPT,0,0,0,1,0,1,0,0,0,0,0,1,0,0\nTHFILE,acc.th\n")
    (tmp_path / "m").mkdir(exist_ok=True)
    (tmp_path / "m" / "acc.th").write_text("0.01\n0.0\n0.1\n0.0\n")
    ui.execute("CHECK")
    err = _checks(ui)
    assert "EDU-42" in err                                             # no PINGRP, no .pin
    assert "EDU-43" in err                                             # STRESS <iter> = 0
    assert "Error 79" in err
    ui.execute("STRESS,0,1,0,0,1")
    ui.execute("PINGRP,2,0,1,1,Sand")                                   # group 2 is excavated soil
    ui.execute("CHECK")
    err = _checks(ui)
    assert "Error 79" not in err and "EDU-43" not in err and "EDU-42" not in err
    assert "EDU-41" in err and "excavated soil elements" in err
    ui.execute("PINDEL")
    ui.execute("PINGRP,3,0,1,1,Sand")
    ui.execute("CHECK")
    err = _checks(ui)
    assert "EDU-41" not in err and "EDU-42" not in err
    # STRESS <iter> = 1 without the non-linear option: EDU-43 (STRESS)
    ui.execute("HOUSEX,0,0,1,0")
    ui.execute("CHECK")
    assert "EDU-43" in _checks(ui)


BLOCK = """MDL,b,b
L,1,1.0,19.62,400,200,0.02,0.02
L,2,1.0,19.62,400,200,0.02,0.02
L,3,1.0,21.58,2000,1000,0.01,0.01
TOPL,1,2
FREQ,1,4,20
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.01,1024,1
WAVE,2,1,1,1,0
HOUSE,9.81,0,0,2,0,0,0,0,0
DYNP,1,0.0001,1.0,0.0001,0.24,Sand
DYNP,2,1.0,0.1,1.0,20.0,Sand
N,1,0,0,-2
N,3,2,0,-2
FILL,1,3
NGEN,2,3,1,3,1,0,1,0
NGEN,2,9,1,9,1,0,0,1
GROUP,1,SOLID
MACT,1
E,1,10,11,14,13,19,20,23,22
EGEN,1,1,1
EGEN,1,3,1,2
ETYPE,1,4,1,2
GROUP,2,SOLID
MACT,2
E,1,1,2,5,4,10,11,14,13
EGEN,1,1,1
EGEN,1,3,1,2
ETYPE,1,4,1,2
M,1,400,200,19.62,0.02,0.02,3
GROUP,3,SOLID
MACT,1
E,1,1,2,5,4,10,11,14,13
EGEN,1,1,1
EGEN,1,3,1,2
EGEN,1,9,1,4
ETYPE,1,8,1,1
INT,1,27,1,1
POINT,0,2,0.9
AOPT,0,0,0,1,0,1,0,0,0,0,0,0,0,0
"""


def test_excstrchk_accepts_the_near_field_soil_of_a_nonlinear_group(tmp_path):
    """A 2 x 2 x 2 block: its centre node 14 is an interior excavation node shared with the near-field soil
    (group 3).  Without the non-linear option CHECK reports EXCSTRCHK (an error that blocks the decks); with
    PINGRP,3 and HOUSEX <nlssi> = 1 the near-field soil replaces the excavated soil there (FV, interaction
    node): no EXCSTRCHK message."""
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(BLOCK)
    ui.execute("CHECK")
    err = (tmp_path / "b" / "b.err").read_text()
    assert "EXCSTRCHK" in err and "Node 14" in err
    ui.execute("PINGRP,3,0,1,1,Sand")
    ui.execute("HOUSEX,0,0,1,1")
    ui.execute("CHECK")
    err = (tmp_path / "b" / "b.err").read_text()
    assert "EXCSTRCHK" not in err
    ui.execute("INT,14,14,1,0")                                        # not an interaction node: reported again
    ui.execute("CHECK")
    assert "EXCSTRCHK" in (tmp_path / "b" / "b.err").read_text()


# ======================================================================================
# COMBXYZSTRAIN (UT-20), NLSSIRESET
# ======================================================================================
def _write74(path, direction, gammas, iteration=2):
    nl = NL.NLFile("FILE74", iteration, 0.65, 1, direction=direction, groups=[NL.NLGroupInfo(3, "SOLID", 2, 0)])
    for k, g in enumerate(gammas, start=1):
        nl.rows.append(NL.NLRow(3, k, g, 0.0, 0.0, 1, 8.0e4, g / 0.65))
    NL.write_nlfile(path, nl)


SAND = SH.DynamicProperty("Sand", np.array([1e-4, 1.0]), np.array([1.0, 0.1]), np.array([1e-4, 1.0]),
                          np.array([0.24, 20.0]))


def test_combxyzstrain_ut20(tmp_path):
    ui = _ui(tmp_path)
    d = tmp_path / "m"
    d.mkdir(exist_ok=True)
    _write74(d / "FILE74X", "X", [0.10, 0.02])
    _write74(d / "FILE74Y", "Y", [0.08, 0.02])
    _write74(d / "FILE74Z", "Z", [0.02, 0.01])
    SH.write_file73(d / "FILE73", {"Sand": SAND})
    assert ui.execute("COMBXYZSTRAIN,FILE74X,FILE74Y,FILE74Z,FILE74")
    out = NL.read_nlfile(d / "FILE74", "FILE74")
    assert round(out.rows[0].gamma_eff, 4) == 0.1296                    # UT-20
    assert out.rows[1].gamma_eff == pytest.approx(0.03)
    assert out.direction == "SRSS X+Y+Z" and out.iteration == 2
    G, b = NL.strain_compatible(SAND, out.rows[0].gamma_eff, 8.0e4)
    assert out.rows[0].G == pytest.approx(G[0]) and out.rows[0].beta == pytest.approx(b[0])
    assert out.rows[0].extra == pytest.approx(NL.srss(0.10, 0.08, 0.02) / 0.65)
    # blank inputs are skipped, the default output is FILE74
    assert ui.execute("COMBXYZSTRAIN,FILE74X,,FILE74Z")
    assert NL.read_nlfile(d / "FILE74").rows[0].gamma_eff == pytest.approx(NL.srss(0.10, 0.02))
    # with the FILE78 of the same iteration the convergence measure is printed
    NL.write_nlfile(d / "FILE78", NL.NLFile("FILE78", 2, 0.65, 1, rows=[
        NL.NLRow(3, 1, 0.0, 8.0e4, 0.02, 1, 8.0e4, 1), NL.NLRow(3, 2, 0.0, 8.0e4, 0.02, 1, 8.0e4, 1)]))
    ui.execute("COMBXYZSTRAIN,FILE74X,FILE74Y,FILE74Z")
    assert any("not converged" in t for t in ui.sink.texts(Kind.INFO))
    # errors: nothing given, a missing file, different elements
    assert not ui.execute("COMBXYZSTRAIN")
    assert not ui.execute("COMBXYZSTRAIN,NOFILE")
    _write74(d / "FILE74W", "W", [0.1])
    assert not ui.execute("COMBXYZSTRAIN,FILE74X,FILE74W")


def test_nlssireset(tmp_path):
    ui = _ui(tmp_path)
    d = tmp_path / "m"
    d.mkdir(exist_ok=True)
    NL.write_liq(d / "m.liq", 1)
    (d / NL.CONVERGENCE_FILE).write_text("# old\n")
    _write74(d / "FILE74", "X", [0.1, 0.02])
    NL.write_nlfile(d / "FILE78", NL.NLFile("FILE78", 2, 0.65, 1, rows=[NL.NLRow(3, 1, 0.0, 8.0e4, 0.02, 1, 8.0e4, 1)]))
    assert ui.execute("NLSSIRESET")
    assert NL.read_liq(d / "m.liq") == 0 and not (d / NL.CONVERGENCE_FILE).exists()
    # the iteration files of the previous analysis are retired (kept as .prev, invisible to STRESS/NLSSIITER)
    assert not (d / "FILE74").exists() and not (d / "FILE78").exists()
    assert NL.read_nlfile(d / "FILE78.prev", "FILE78").iteration == 2 and (d / "FILE74.prev").exists()
    assert any("FILE78.prev" in t for t in ui.sink.texts(Kind.CONFIRM))
    assert ui.execute("NLSSIRESET")                                     # nothing to retire: still fine
    nomdl = Interpreter(cwd=tmp_path)
    assert not nomdl.execute("NLSSIRESET")                              # needs MDL


# ======================================================================================
# NLSSIITER
# ======================================================================================
@pytest.fixture
def fake_step():
    """ZZTESTNLSTEP,<k>: emulates the HOUSE/ANALYS/STRESS pass of iteration k -- .liq = 1, FILE90 with a new
    FILE4 hash, FILE78 carrying that hash and the digests of m.hou and m.pin (as HOUSE writes it) and FILE74
    of the response -- whose change in G is 60 % x 0.5^k (flat curve G/Gmax = 0.5, D = 5 %), so the
    iterations converge (< 2 %) at k = 5; k < 0 fails."""
    calls = []

    @registry.command("ZZTESTNLSTEP", tier="P0")
    def _step(c):
        k = c.int(1)
        if k < 0:
            c.fail("deliberate failure")
        calls.append(k)
        d = Path(c.model.path)
        h4 = f"file4-of-iteration-{k}"
        write_container(d / "FILE90", "FILE90", {}, {"int_hash": "i", "layer_hash": "l", "file4_hash": h4})
        NL.write_liq(d / "m.liq", 1)
        pin = NL.pin_digest(NL.read_pin(d / "m.pin")) if (d / "m.pin").exists() else ""
        rows78 = [NL.NLRow(3, 1, 0.01, 100.0 * (0.5 + 0.3 * 0.5 ** k), 0.05, 1, 100.0, 1)]
        NL.write_nlfile(d / "FILE78", NL.NLFile("FILE78", k, 0.65, 1, rows=rows78, file4=h4, pin=pin,
                                                deck=NL.deck_digest(d / "m.hou")))
        NL.write_nlfile(d / "FILE74", NL.NLFile("FILE74", k, 0.65, 1, file4=h4,
                                                rows=[NL.NLRow(3, 1, 0.02, 0, 0, 1, 100.0, 0.03)]))

    yield calls
    registry.REGISTRY.pop("ZZTESTNLSTEP", None)


def _flat(d):
    SH.write_file73(d / "FILE73", {"Flat": SH.DynamicProperty("Flat", np.array([1e-4, 10.0]), np.array([0.5, 0.5]),
                                                             np.array([1e-4, 10.0]), np.array([5.0, 5.0]))})
    NL.write_pin(d / "m.pin", NL.PinData(0.65, 1, [NL.PinGroup(3, 1, 1, 0, [NL.PinMaterial(1.0, 1.0, 1)])]))
    (d / "m.hou").write_text("SASSI-EDU HOUSE DECK v1\n* model_hash = 1\n[params]\ngravity = 9.81\n")


def test_nlssiiter_converges_and_stops(tmp_path, fake_step):
    ui = _ui(tmp_path)
    d = tmp_path / "m"
    d.mkdir(exist_ok=True)
    _flat(d)
    assert ui.execute('VAR,BODY,"ZZTESTNLSTEP,#"')
    assert ui.execute("NLSSIITER,BODY,8")
    assert fake_step == [1, 2, 3, 4, 5]                                 # '#' = pass number; stop at < 2 %
    rows = NL.read_convergence(d)
    assert [r["iteration"] for r in rows] == [1, 2, 3, 4, 5] and rows[-1]["converged"] == 1
    assert rows[0]["max_dG_pct"] == pytest.approx(-30.0) and rows[-1]["max_dG_pct"] == pytest.approx(-1.875)
    assert any("converged after 5 pass(es)" in t for t in ui.sink.texts(Kind.INFO))
    # already converged (the analysis in progress): no pass at all and no history row written
    fake_step.clear()
    history = (d / NL.CONVERGENCE_FILE).read_text()
    (d / "m.hou").write_text("SASSI-EDU HOUSE DECK v1\n* model_hash = 2\n[params]\ngravity = 9.81\n")  # AFWRITE again
    assert ui.execute("NLSSIITER,BODY,8") and fake_step == []
    assert any("already converged" in t for t in ui.sink.texts(Kind.INFO))
    assert (d / NL.CONVERGENCE_FILE).read_text() == history
    # NLSSIITER without a variable lists the history
    n = len(ui.sink.texts(Kind.INFO))
    assert ui.execute("NLSSIITER")
    assert any("max dG/G" in t for t in ui.sink.texts(Kind.INFO)[n:])


@pytest.mark.parametrize("change,why", [
    ("pin", "m.pin has changed"), ("deck", "m.hou has changed"), ("file90", "FILE90"), ("liq", ".liq is not 1"),
    ("reset", ".liq is not 1")])
def test_nlssiiter_reruns_when_the_converged_iteration_is_not_of_the_current_analysis(tmp_path, fake_step, change,
                                                                                       why):
    """Review finding: 'already converged' must not be reported for the files of another analysis -- after
    NLSSIRESET, a new .pin (ESF, curves, factors), a changed model (HOUSE deck data) or another FILE4."""
    ui = _ui(tmp_path)
    d = tmp_path / "m"
    d.mkdir(exist_ok=True)
    _flat(d)
    ui.execute('VAR,BODY,"ZZTESTNLSTEP,#"')
    assert ui.execute("NLSSIITER,BODY,8") and fake_step[-1] == 5
    fake_step.clear()
    if change == "pin":
        NL.write_pin(d / "m.pin", NL.PinData(0.6, 1, [NL.PinGroup(3, 1, 1, 0, [NL.PinMaterial(1.0, 1.0, 1)])]))
    elif change == "deck":
        (d / "m.hou").write_text("SASSI-EDU HOUSE DECK v1\n[params]\ngravity = 9.80665\n")
    elif change == "file90":
        write_container(d / "FILE90", "FILE90", {}, {"int_hash": "i", "layer_hash": "l", "file4_hash": "other"})
    elif change == "liq":
        NL.write_liq(d / "m.liq", 0)
    else:
        assert ui.execute("NLSSIRESET")
    n = len(ui.sink.texts(Kind.INFO))
    assert ui.execute("NLSSIITER,BODY,8") and fake_step == [1, 2, 3, 4, 5]
    infos = ui.sink.texts(Kind.INFO)[n:]
    assert not any("already converged" in t for t in infos)
    if change != "reset":                                               # after NLSSIRESET there is no old iteration
        assert any("not that of the analysis in progress" in t and why in t for t in infos), infos
    assert [r["iteration"] for r in NL.read_convergence(d)][-1] == 5


def test_nlssiiter_limit_tolerances_and_concatenation(tmp_path, fake_step):
    ui = _ui(tmp_path)
    d = tmp_path / "m"
    d.mkdir(exist_ok=True)
    _flat(d)
    ui.execute('VAR,B1,"ZZTESTNLSTEP,#"')
    ui.execute('VAR,B2,"SHOWVAR,B1"')
    assert ui.execute("NLSSIITER,B1+B2,3")
    assert fake_step == [1, 2, 3]
    assert any("not converged after 3 pass(es)" in t for t in ui.sink.texts(Kind.WARNING))
    NL.reset_convergence(d)
    (d / "FILE74").unlink()
    fake_step.clear()
    assert ui.execute("NLSSIITER,B1,8,20")                               # tolG 20 %: converged at k = 2 (15 %)
    assert fake_step == [1, 2]


def test_nlssiiter_errors(tmp_path, fake_step):
    ui = _ui(tmp_path)
    d = tmp_path / "m"
    d.mkdir(exist_ok=True)
    _flat(d)
    assert not ui.execute("NLSSIITER,NOVAR")
    ui.execute("VAR,EMPTYV")
    assert not ui.execute("NLSSIITER,EMPTYV")
    ui.execute('VAR,BAD,"ZZTESTNLSTEP,-1"')
    assert not ui.execute("NLSSIITER,BAD")
    assert any("failed; iterations stopped" in t for t in ui.sink.texts(Kind.ERROR))
    ui.execute('VAR,NOFILES,"SHOWVAR,BAD"')                             # a body that writes no FILE78/FILE74
    for nm in ("FILE78", "FILE74"):
        if (d / nm).exists():
            (d / nm).unlink()
    assert not ui.execute("NLSSIITER,NOFILES,2")
    assert any("no FILE78/FILE74" in t for t in ui.sink.texts(Kind.ERROR))
    ui.execute('VAR,OK,"ZZTESTNLSTEP,#"')
    assert not ui.execute("NLSSIITER,OK,0")
    assert not ui.execute("NLSSIITER,OK,8,0")
    assert not ui.loops                                                  # the loop frame was removed
    # a pass that writes no new FILE74 is not an iteration: stop, no convergence row from stale files
    assert ui.execute("NLSSIITER,OK,2")                                  # FILE78/FILE74 of iteration 2 now exist
    rows = NL.read_convergence(d)
    assert not ui.execute("NLSSIITER,NOFILES,3")
    assert any("FILE74 was not written" in t for t in ui.sink.texts(Kind.ERROR))
    assert NL.read_convergence(d) == rows and not ui.loops
