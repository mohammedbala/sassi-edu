"""UT-08: the CHECK catalogue -- one minimal model per Error 1-128 and Warning 1-11 that violates that
rule (spec 11 section 5), the fatal stop on Errors 42/43, no Error 8 for zero-length springs, and the
break count (requirements section 6.2, D-CHK-01..09).

Every case starts from a base model that enables every module and passes CHECK without errors; the
case appends the command lines (or applies a model edit) that break one rule.  The test asserts
that the message number appears and that no *other* error appears except the listed companions.
"""
from __future__ import annotations

import pytest

from sassi.prep import Interpreter, Kind
from sassi.prep.check import ERRORS, MODEL, WARNINGS, CheckOptions, run_check


def _files(d):
    (d / "H1.rsi").write_text("0.5 0.2\n1 0.4\n5 0.8\n10 0.6\n33 0.3\n")
    acc = "\n".join(f"{0.01 * ((-1) ** k) * (k % 7)}" for k in range(200))
    (d / "H1.acc").write_text("0.005\n" + acc + "\n")


BASE = """
MDL,base,{dir}
TIT,CHECK catalogue base model
N,1,0,0,0
N,2,10,0,0
N,3,10,10,0
N,4,0,10,0
N,5,0,0,10
N,6,5,0,0
N,7,0,0,12
D,7,7,1,1,ALL
D,2,4,1,1,ROTZ
INT,1,4,1,1
M,1,4.32e5,0.25,0.15,0.05,0.05
R,1,1,0.8,0.8,0.2,0.1,0.1
SC,1,100,100,100,0,0,0,0.02
MXR,1,1,100
GROUP,1,SHELL
E,1,1,2,3,4
THICK,1,1,1,0.5
GROUP,2,BEAMS
E,1,1,5,6
GROUP,3,SPRING
E,1,5,7
GROUP,4,GENERAL
E,1,5,7
MT,5,1,1,1
F,5,1,0,0
L,1,5,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
TOPL,1,1
FREQ,1,2,4,8,16
SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,4096,1
WAVE,2,1,1,1,0
POINT,0,0,4.5
HOUSE,32.2,0,0,2,0,0,0,0,0
FORCE,0
ANALYS,0,0,0,0,1,0,0,0,0,0,0
THFILE,H1.acc
DAMP,0.02,0.05
MOTION,0,0,0,0,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1
NOUT,1,1,1,0,0,1,1,1-5
STRESS,0,0,1,0,1
EOUT,1,1,1,1,1,1,0,0,0,0,0,0,1,1
RELD,1,0,1
RDND,5,1,0,0,0,0,0
EQUAKE,0,5,11975,0.05,20,0,1
RSIN,1,H1.rsi
RSOUT,1,H1.rso
ACCOUT,1,H1out.acc
DYNP,1,0.0001,1,0.0001,0.5,Sand
DYNP,2,0.1,0.5,0.1,10,Sand
SPRO,1,1,Sand
SPRO,2,1,Sand
SPRO,3,2,Sand
SOIL,100,32.2,0,1,1,8,0.65,1,0
SRS,1,1,0
AOPT,1,1,0,1,1,1,0,1,1,1,1,1,1,0
"""

_SIXTEEN = "\n".join(["SPRO,0"] + [f"DYNP,1,0.0001,1,0.0001,0.5,P{k}\nDYNP,2,0.1,0.5,0.1,10,P{k}\nSPRO,{k},1,P{k}"
                                   for k in range(1, 18)])

#: number -> (extra command lines, companion errors allowed, model edit)
ERROR_CASES = {
    1: ("GRAVITY,-1", set(), None),
    2: ("SYMM,1,0,1,2,99", set(), None),
    3: ("SYMM,1,0,1", set(), None),
    4: ("SYMM,1,0,1,1", set(), None),
    5: ("N,8,20,0,0\nSYMM,1,0,1,2,8", set(), None),
    6: ("GROUP,1\nE,3,1,2,3,4", set(), None),
    7: ("GROUP,1\nE,2,1,2", set(), None),
    8: ("GROUP,2\nE,2,1,1,6", set(), None),
    9: ("N,8,20,0,0\nGROUP,1\nE,2,1,2,8", set(), None),
    10: ("GROUP,2\nKI,1,1,1,1,0,0,0,0,0\nKJ,1,1,1,1,0,0,0,0,0", set(), None),
    11: ("N,8,20,10,0\nGROUP,1\nE,2,1,2,3,8", set(), None),
    12: ("N,8,0,10,5\nGROUP,1\nE,2,1,2,3,8", set(), None),
    13: ("GROUP,1\nMSET,1,1,1,5", set(), None),
    14: ("M,1,-4.32e5,0.25,0.15,0.05,0.05", set(), None),
    15: ("M,1,4.32e5,-0.1,0.15,0.05,0.05", set(), None),
    16: ("M,1,4.32e5,0.25,-0.15,0.05,0.05", set(), None),
    17: ("M,1,4.32e5,0.25,0.15,-0.05,0.05", set(), None),
    18: ("M,1,4.32e5,0.25,0.15,0.05,-0.05", set(), None),
    19: ("TOPL,0\nTOPL,1,3", set(), None),
    20: ("L,1,0,0.12,1500,800,0.05,0.05", set(), None),
    21: ("L,1,5,-0.12,1500,800,0.05,0.05", set(), None),
    22: ("L,1,5,0.12,-1500,800,0.05,0.05", set(), None),
    23: ("L,1,5,0.12,1500,-800,0.05,0.05", set(), None),
    24: ("L,1,5,0.12,1500,800,-0.05,0.05", set(), None),
    25: ("L,1,5,0.12,1500,800,0.05,-0.05", set(), None),
    26: ("GROUP,2\nRSET,1,1,1,9", set(), None),
    27: ("R,1,0,0.8,0.8,0.2,0.1,0.1", set(), None),
    28: ("R,1,1,-0.8,0.8,0.2,0.1,0.1", set(), None),
    29: ("R,1,1,0.8,-0.8,0.2,0.1,0.1", set(), None),
    30: ("R,1,1,0.8,0.8,0,0.1,0.1", set(), None),
    31: ("R,1,1,0.8,0.8,0.2,-0.1,0.1", set(), None),
    32: ("R,1,1,0.8,0.8,0.2,0.1,-0.1", set(), None),
    33: ("GROUP,3\nRSET,1,1,1,9", set(), None),
    34: ("SC,1,-100,100,100,0,0,0,0.02", set(), None),
    35: ("SC,1,100,-100,100,0,0,0,0.02", set(), None),
    36: ("SC,1,100,100,-100,0,0,0,0.02", set(), None),
    37: ("SC,1,100,100,100,-1,0,0,0.02", set(), None),
    38: ("SC,1,100,100,100,0,-1,0,0.02", set(), None),
    39: ("SC,1,100,100,100,0,0,-1,0.02", set(), None),
    40: ("MT,99,1,1,1", set(), None),
    41: ("GROUP,1\nE,2,1,2,3,99", set(), None),
    44: ("SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,4096,2", set(), None),
    45: ("SITE,0,0,0,20,2,0,0,1,2048,1,0,0.005,4096,1", set(), None),
    46: ("TOPL,0", set(), None),
    47: ("SITE,0,1,0,3,2,1,0,1,2048,1,0,0.005,4096,1", set(), None),
    48: ("SITE,0,1,-1,20,2,1,0,1,2048,1,0,0.005,4096,1", set(), None),
    49: ("SITE,0,1,0,20,2,1,0,1,2048,1,0,-0.005,4096,1", {69}, None),
    50: ("SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,-4096,1", set(), None),
    51: ("WAVE,2,0,1,1,0", set(), None),
    52: ("WAVE,2,1,1,1,400", set(), None),
    53: ("SITE,0,1,0,20,2,1,0,0,2048,1,0,0.005,4096,1", set(), None),
    54: ("WAVE,2,1,1.5,1,0", {55}, None),
    55: ("WAVE,2,1,0.5,0.5,0", set(), None),
    56: ("POINT,0,-1,4.5", set(), None),
    57: ("POINT,0,0,0", set(), None),
    58: ("HOUSE,32.2,0,0,2,0,1,0,0,0\nINCOH,0.05", set(), None),
    59: ("HOUSE,32.2,0,0,2,0,1,0,0,0\nINCOH,0.1,0.1,0.2,0", set(), None),
    60: ("N,8,10,10,-5\nINT,8,8,1,1\nPOINT,0,1,4.5\nHOUSE,32.2,0,0,2,0,1,0,0,0\nINCOH,0.1,0.1,0.2,0.5,0",
         {"EDU-21"}, None),
    61: ("FDEL,5", set(), None),
    62: ("F,99,1,0,0", set(), None),
    63: ("ANALYS,0,0,0,0,1,0,400,0,0,0,0", set(), None),
    64: ("NOUT,0", set(), None),
    65: ("NOUT,1,1,1,0,0,1,1,99", set(), None),
    66: ("RDND,5,1,0,0,0,0,0", set(), None),
    67: ("MOTION,0,0,-1,0,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1", set(), None),
    68: ("MOTION,0,0,0,-1,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1", set(), None),
    69: ("MOTION,0,0,0,0,0,-0.1,100,301,1,0,1,0,0,0,0,1,0,0,1", set(), None),
    70: ("MOTION,0,0,0,0,0,0.1,-100,301,1,0,1,0,0,0,0,1,0,0,1", set(), None),
    71: ("MOTION,0,0,0,0,0,0.1,100,-1,1,0,1,0,0,0,0,1,0,0,1", set(), None),
    72: ("DAMP,0\nDAMP,1.5", set(), None),
    73: ("THFILE,missing.acc", set(), None),
    74: ("MOTION,0,0,0,0,0,0.1,100,301,1,0,500,0,0,0,0,1,0,0,1", set(), None),
    75: ("MOTION,0,0,0,0,0,0.1,100,301,1,0,1,-1,0,0,0,1,0,0,1", set(), None),
    76: ("MOTION,0,0,0,0,0,0.1,100,301,1,0,50,10,0,0,0,1,0,0,1", set(), None),
    77: ("MOTION,0,0,0,0,0,0.1,100,301,0,0,1,0,0,0,0,1,0,0,1", set(), None),
    78: ("MOTION,0,0,0,0,0,0.1,100,301,1,0.5,1,0,0,0,0,1,0,0,1", set(), None),
    79: ("EOUT,0", set(), None),
    80: ("EOUT,1,1,1,1,1,1,0,0,0,0,0,0,9,1", set(), None),
    81: ("EOUT,1,1,1,1,1,1,0,0,0,0,0,0,1,5", set(), None),
    82: ("EOUT,1,1,1,1,1,1,0,0,0,0,0,0,1,1", set(), None),
    83: ("GROUP,4\nRSET,1,1,1,9", set(), None),
    84: ("RSIN,1,", set(), None),
    85: ("RSIN,1,missing.rsi", set(), None),
    86: ("RSOUT,1,", set(), None),
    87: ("ACCOUT,1,", set(), None),
    88: ("EQUAKE,1,5,11975,0.05,20,0,1", set(), None),
    89: ("EQUAKE,0,4,11975,0.05,20,0,1", set(), None),
    90: ("EQUAKE,0,5,0,0.05,20,0,1", set(), None),
    91: ("EQUAKE,0,0,11975,0.05,20,0,1", set(), None),
    92: ("EQUAKE,0,5,11975,0.05,0,0,1", set(), None),
    93: ("EQUAKE,0,5,11975,0.05,20,1,1", set(), None),
    94: ("EQUAKE,0,5,11975,0.05,20,1,1\nCORR,1,0,1.5", set(), None),
    95: ("SPRO,1,1,\nSPRO,2,1,\nSPRO,3,2,", set(), None),
    96: ("EDUOPT,LIMITS,PREP\n" + _SIXTEEN, set(), None),
    97: ("DYNP,1,,,0.001,1,Clay\nSPRO,1,1,Clay", set(), None),
    98: ("DYNP,3,1,1.5,1,15,Sand", set(), None),
    99: ("DYNP,1,0.001,1,,,Clay\nSPRO,1,1,Clay", set(), None),
    100: ("SOIL,0,32.2,0,1,1,8,0.65,1,0", set(), None),
    101: ("SOIL,100,32.2,0,1,1,8,0.65,1,-1", set(), None),
    102: ("SOIL,4000,32.2,0,1,1,8,0.65,1.0,(5X,F10.4)", set(), None),
    103: ("SOIL,100,32.2,-1,1,1,8,0.65,1,0", set(), None),
    104: ("SOILX,0,0,0.1,9", set(), None),
    105: ("SOIL,100,32.2,0,1,1,-1,0.65,1,0", set(), None),
    106: ("SOIL,100,32.2,0,1,1,8,1.5,1,0", set(), None),
    107: ("DAMP,0", set(), None),
    108: ("SOIL,100,32.2,0,1,1,8,0.65,0,0", set(), None),
    109: ("SSAF,1,1,0,0,9,0.1,amp", set(), None),
    110: ("SSAF,1,1,0,0,2,0,amp", set(), None),
    111: ("SFOU,1,0,0,0,-1,0", set(), None),
    112: ("SFOU,1,0,0,0,0,-1", set(), None),
    113: ("HOUSE,32.2,0,0,2,0,0,1,0,0\nWPASS,0,0,1", set(), None),
    114: ("HOUSE,32.2,0,0,2,0,0,1,0,0\nWPASS,1e9,0,9", set(), None),
    115: ("HOUSE,32.2,0,0,2,0,0,1,1,0", set(), None),
    116: ("HOUSE,32.2,0,0,2,0,0,1,1,0\nME,1,99,4\nAMP,1,1,1,1,1", set(), None),
    117: ("HOUSE,32.2,0,0,2,0,0,1,1,0\nME,1,1,99\nAMP,1,1,1,1,1", set(), None),
    118: ("HOUSE,32.2,0,0,2,0,0,1,1,0\nME,1,1,4\nAMP,1,11,1,1,1", set(), None),
    119: ("HOUSE,32.2,0,0,2,0,0,1,1,0\nME,1,1,4\nAMP,1,1,1,1", set(), None),
    120: ("", set(), lambda m: m.freq_sets.__setitem__(1, [])),
    121: ("AOPT,1,1,0,1,1,1,0,1,1,1,1,1,1,1", set(), None),
    122: ("GROUP,5,SHELL\nMACT,9\nE,1,1,2,3,4\nMACT,1\nAOPT,1,1,0,1,1,1,0,1,1,1,1,1,1,1\nP,1,5,1,1,1", {13}, None),
    123: ("AOPT,1,1,0,1,1,1,0,1,1,1,1,1,1,1\nP,1,9,1,1,1", set(), None),
    124: ("D,1,1,1,1,DISP", set(), None),
    125: ("NLSOIL,1\nNLSLAYER,1", set(), None),
    126: ("AOPT,1,1,0,1,1,1,0,1,1,1,1,1,1,1\nP,1,1,1,1,2", {"EDU-44"}, None),   # base model shares a material between panel groups (D-NON-12)
    127: ("AOPT,1,1,0,1,1,1,0,1,1,1,1,1,1,1\nP,1,1,1,1,1\nS,1,3,1,1,1,3", {"EDU-44"}, None),   # base model shares a material between panel groups (D-NON-12)
    128: ("AOPT,1,1,0,1,1,1,0,1,1,1,1,1,1,1\nP,1,1,1,2,1", {"EDU-44"}, None),   # base model shares a material between panel groups (D-NON-12)
}

WARNING_CASES = {
    1: "N,9,30,0,0",
    2: "N,8,100,5,0\nGROUP,1\nE,2,2,3,8",
    3: "N,8,0,10,0.1\nGROUP,1\nE,2,1,2,3,8",
    4: "N,8,50,50,50",
    5: "D,5,5,1,1,UZ",
    6: "MR,5,1,1,1\nD,5,5,1,1,ROTX",
    7: "GROUP,5,SHELL",
    8: "EDUOPT,LIMITS,PREP\nTOPL,0\n" + "\n".join("TOPL," + ",".join(["1"] * 20) for _ in range(5)) + "\nTOPL,1",
    9: "SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,3000,1",
    10: "D,5,5,1,1,UX",
    11: "MM,5,1,0,0\nD,5,5,1,1,ROTX",
}


def _run(tmp_path, extra: str = "", edit=None, modules=None):
    ui = Interpreter(cwd=tmp_path)
    mdir = tmp_path / "m"
    mdir.mkdir(exist_ok=True)
    _files(mdir)
    ui.run_text(BASE.format(dir=mdir) + "\n" + extra)
    if edit is not None:
        edit(ui.model)
    return ui, run_check(ui.model, modules=modules, dirs=[mdir, tmp_path])


def test_base_model_has_no_errors(tmp_path):
    ui, rep = _run(tmp_path)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    assert rep.errors() == [], rep.format()
    assert rep.modules == ["EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "COMBIN", "MOTION",
                           "STRESS", "RELDISP"]
    assert rep.sections[0] == MODEL


@pytest.mark.parametrize("number", sorted(ERROR_CASES))
def test_error_catalogue(tmp_path, number):
    extra, allowed, edit = ERROR_CASES[number]
    ui, rep = _run(tmp_path, extra, edit)
    nums = rep.numbers("Error")
    assert number in nums, rep.format()
    unexpected = nums - {number} - set(allowed)
    assert not unexpected, f"Error {number}: unexpected {unexpected}\n{rep.format()}"
    text = rep.format()
    title = ERRORS[number].split("{")[0].strip()
    assert f"Error {number} : {title}" in text


@pytest.mark.parametrize("number", sorted(WARNING_CASES))
def test_warning_catalogue(tmp_path, number):
    ui, rep = _run(tmp_path, WARNING_CASES[number])
    assert number in rep.numbers("Warning"), rep.format()
    assert rep.errors() == [], rep.format()
    title = WARNINGS[number].split("{")[0].strip()
    assert f"Warning {number} : {title}" in rep.format()


def test_every_catalogue_number_is_covered():
    assert set(ERROR_CASES) | {42, 43} == set(ERRORS)
    assert set(WARNING_CASES) == set(WARNINGS)


def test_fatal_errors_stop_check(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    rep = run_check(ui.model)                       # default AOPT: HOUSE enabled, no nodes
    assert rep.fatal and rep.numbers("Error") == {42} and rep.sections == [MODEL]
    assert all(rep.blocked(m) for m in rep.modules)
    ui.run_text("N,1,0,0,0\nN,2,1,0,0")
    rep = run_check(ui.model)
    assert rep.fatal and rep.numbers("Error") == {43} and rep.sections == [MODEL]
    assert "CHECK stopped" in rep.format()


def test_no_error_8_for_zero_length_springs(tmp_path):
    ui, rep = _run(tmp_path, "N,8,0,0,10\nGROUP,3\nE,2,5,8")
    assert 8 not in rep.numbers("Error") and rep.errors() == []
    ui, rep = _run(tmp_path, "N,8,0,0,10\nGROUP,4\nE,2,5,8")         # 2-node GENERAL: Error 8
    assert 8 in rep.numbers("Error")


def test_break_count_per_message_type(tmp_path):
    """Break Check at N: at most N errors and N warnings per module (spec 05a section 3, OQ19,
    D-CHK-02) -- the count is per type (errors / warnings), not per message number."""
    ui, rep = _run(tmp_path, "N,160,1,1,1")          # gaps 8..159 -> 152 Warning 1, plus W4 for 160
    w1 = [m for m in rep.warnings(MODEL) if m.number == 1]
    assert len(w1) == 152
    n_w = len(rep.warnings(MODEL))
    assert n_w > 152 and rep.warnings(MODEL)[-1].number == 4
    txt = rep.format(CheckOptions(break_at=100))
    sec = txt.split("Errors and Warnings for MODEL")[1].split("Errors and Warnings for")[0]
    assert sum(1 for ln in sec.splitlines() if ln.startswith("Warning ")) == 100
    assert sum(1 for ln in sec.splitlines() if ln.startswith("Warning 1 :")) == 100
    assert "Warning 4 : Unused Node 160" not in sec                  # beyond the 100 warnings shown
    assert f"{n_w - 100} more warnings not shown (Break Check at 100)" in sec
    assert f"MODEL: 0 errors, {n_w} warnings" in sec
    txt = rep.format(CheckOptions(break_at=200))
    assert "Warning 4 : Unused Node 160" in txt and "not shown" not in txt
    txt = rep.format(CheckOptions(show_warnings=False))
    assert "Warning 1 :" not in txt


def test_model_errors_gate_house_and_file4_modules(tmp_path):
    ui, rep = _run(tmp_path, "M,1,-4.32e5,0.25,0.15,0.05,0.05")        # Error 14 under MODEL
    assert rep.errors(MODEL) and [m.module for m in rep.errors()] == [MODEL]
    assert rep.blocked("HOUSE") and rep.blocked("ANALYS") and rep.blocked("STRESS")
    assert not rep.blocked("SITE") and not rep.blocked("MOTION") and not rep.blocked("POINT")


def test_shared_variable_errors_reported_under_each_module(tmp_path):
    ui, rep = _run(tmp_path, "GRAVITY,-1")
    mods = {m.module for m in rep.errors() if m.number == 1}
    assert mods == {"EQUAKE", "SITE", "HOUSE", "FORCE", "ANALYS", "MOTION", "STRESS", "RELDISP"}
    ui, rep = _run(tmp_path, "THFILE,missing.acc")
    assert {m.module for m in rep.errors() if m.number == 73} == {"SOIL", "MOTION", "STRESS", "RELDISP"}


def test_err_format_headers_in_run_order(tmp_path):
    ui, rep = _run(tmp_path, "RSIN,1,missing.rsi\nTHFILE,missing.acc")
    lines = rep.format().splitlines()
    assert lines[0] == "CHECK: Errors and Warning for - base"
    heads = [ln[len("Errors and Warnings for "):] for ln in lines if ln.startswith("Errors and Warnings for ")]
    assert heads == ["MODEL", "EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "COMBIN", "MOTION",
                     "STRESS", "RELDISP"]
    assert "Error 85 : RS Input File 1 Does Not Exist  [missing.rsi]" in lines
    assert "EQUAKE: 1 errors, 0 warnings" in lines
