"""Typed option records, dataclass views for the GUI, shared variables, X-commands, list semantics
(UT-04) and the WRITE -> INP round trip of every option command (UT-03 extended to options)."""
from __future__ import annotations

import dataclasses

import pytest

from sassi.model import make_record
from sassi.prep import Interpreter, Kind
from sassi.prep import options as O
from sassi.prep.writer import write_pre


def test_typed_records_have_new_model_defaults():
    m = Interpreter().model
    site = O.get_record(m, "SITE")
    assert (site.mode1, site.mode2, site.nl, site.delt, site.nft, site.freq, site.freq2) == (1, 1, 20, 0.005, 4096, 1, 2048)
    mo = O.get_record(m, "MOTION")
    assert (mo.freq1, mo.freq2, mo.fstep, mo.mult, mo.get("max"), mo.cplx, mo.interp) == (0.1, 100.0, 301, 1.0, 0.0, 1, 1)
    assert O.get_record(m, "ANALYS").prnt == 1                              # Print Amplitude Only on (D-UI-03)
    assert O.get_record(m, "SOIL").ratio == 0.65 and O.get_record(m, "EQUAKE").dur == 20.0
    assert O.get_record(m, "STRESS").save == 1 and O.get_record(m, "RELD").reldisoutput == 1
    assert O.aopt_modules(m) == ["SITE", "POINT", "HOUSE", "ANALYS", "MOTION"]        # D-AFW-05


def test_records_keep_the_typed_tokens_and_type_values():
    ui = Interpreter()
    ui.execute("WPASS,1.0E9,0.0,3")
    ui.execute("SITE,0,1,,20,3,1,0,1,2048,1,0,5D-3,4096,1")
    w = ui.model.options.record("WPASS")
    assert type(w) is O.WpassRecord and w.to_tokens() == ["1.0E9", "0.0", "3"] and w.appv == 1e9
    s = ui.model.options.record("SITE")
    assert s.fstep == 0.0 and s.delt == 0.005 and s.arg(3) is None
    assert make_record("SITE", ["0"]).problems() == []


def test_handlers_validate_with_check_messages():
    ui = Interpreter()
    ui.execute("SITE,0,1,0,3")
    assert any("Error 47 : Illegal Number of Layers for Halfspace Simulation" in w for w in ui.sink.texts(Kind.WARNING))
    assert ui.model.options.record("SITE").nl == 3                      # stored: CHECK is the authority
    assert not ui.execute("POINT,0,x,1")                                 # non-numeric: command skipped
    assert ui.model.options.record("POINT") is None
    ui.execute("MOTION,0,0,0,0,0,0.1,100,301,0,0")
    assert any("Error 77" in w for w in ui.sink.texts(Kind.WARNING))
    ui.execute("WAVE,2,1,1.5,1,0")
    assert any("Error 54 : Illegal Value for Wave 2 Ratio at Frequency 1" in w for w in ui.sink.texts(Kind.WARNING))
    ui.execute("SITE,0,1,0,20,3,1,0,1,2048,1,0,0.005,4096,1,99")
    assert any("arguments after argument 14 ignored" in w for w in ui.sink.texts(Kind.WARNING))
    assert len(ui.model.options.record("SITE")) == 14


def test_from_args_and_to_command_for_the_gui():
    o = O.from_args("SITE,0,1,0,20,3")
    assert dataclasses.is_dataclass(o) and o.hs == 3 and o.nft == 4096 and o.delt == 0.005
    f = {x.name: x for x in dataclasses.fields(o)}
    assert f["nl"].metadata["label"] == "Number of Generated Layers" and f["nl"].metadata["position"] == 4
    o.nl = 10
    assert O.to_command(o) == "SITE,0,1,0,10,3,1,0,1,2048,1,0,0.005,4096,1"
    assert O.to_command(O.from_args(["0", "24"], command="EQUAKE")) == "EQUAKE,0,24,11975,0.05,20,0,1"
    assert O.to_command(O.from_args("ANALYS,0,0,0,1")) == "ANALYS,0,0,0,1,1,0,0,0,0,0,0"
    assert O.to_command(O.from_args("HOUSEX,0,1,5")) == "HOUSEX,0,1,5"           # trailing defaults omitted
    assert O.to_command(O.from_args("SSAF,1,1,0,1,3,0.1,SAF top, vs rock")) == "SSAF,1,1,0,1,3,0.1,SAF top, vs rock"
    rec = make_record("SITE", ["0", "1", "", "20"])
    assert O.to_command(rec) == "SITE,0,1,,20"
    back = O.record_from_options(O.options_from_record(rec))
    assert back.nl == 20 and back.fstep == 0.0
    with pytest.raises(ValueError):
        O.from_args("SITE,a")
    names = [d["name"] for d in O.dialog_fields("MOTION")]
    assert names[:4] == ["opmode", "out", "step", "dur"] and len(names) == 19


def test_gui_command_round_trips_through_the_interpreter():
    ui = Interpreter()
    o = O.from_args("MOTION")
    o.dur, o.interp, o.mult, o.max = 25.0, 6, 0.0, 0.3
    ui.execute(O.to_command(o))
    got = O.options_from_record(ui.model.options.record("MOTION"))
    assert got == o


def test_shared_variables_have_one_storage_location():
    ui = Interpreter()
    assert O.SHARED_VARIABLES["delt"].owner == "SITE" and "MOTION" in O.SHARED_VARIABLES["delt"].tabs
    ui.execute("SITE,0,1,0,20,3,1,0,1,2048,1,0,0.005,4096,1")
    for line in O.shared_command(ui.model, "delt", 0.01):
        ui.execute(line)
    assert O.shared_value(ui.model, "delt") == 0.01 and ui.model.options.record("SITE").hs == 3
    assert O.shared_command(ui.model, "gravity", 9.81) == ["GRAVITY,9.81"]
    for line in O.shared_command(ui.model, "damp", [0.02, 0.05]):
        ui.execute(line)
    assert ui.model.damp == [0.02, 0.05]
    for line in O.shared_command(ui.model, "rec2", 4096):
        ui.execute(line)
    assert O.get_record(ui.model, "MOTION").rec2 == 4096


def test_x_commands_not_stored_when_default():
    ui = Interpreter()
    ui.execute("HOUSEX,0,0,1")
    assert ui.model.options.record("HOUSEX") is None
    ui.execute("HOUSEX,1")
    assert ui.model.options.record("HOUSEX").optimize == 1
    ui.execute("HOUSEX")
    assert ui.model.options.record("HOUSEX") is None
    ui.execute("SOILX,1,0,0.2,3,/data/rock, outcrop.acc")
    assert ui.model.options.record("SOILX").file == "/data/rock, outcrop.acc"


def test_eduopt_keys_and_reset():
    ui = Interpreter()
    ui.execute("EDUOPT,hslaw,GEOMETRIC")
    assert O.eduopt(ui.model, "HSLAW") == "GEOMETRIC"
    ui.execute("EDUOPT,NFFTROUND,SIDEWAYS")
    assert any("expected one of NEAREST, UP" in w for w in ui.sink.texts(Kind.WARNING))
    assert not ui.execute("EDUOPT,GEOMTOL,abc")
    ui.execute("EDUOPT,HSLAW")
    assert O.eduopt(ui.model, "HSLAW") == "UNIFORM"          # default per D-W1-01


def test_ut04_list_semantics():
    ui = Interpreter()
    ui.run_text("DAMP,0.02\nDAMP,0.05\nFREQ,1,2,4\nFREQ,1,6\nFREQ,2,9\nTOPL,1,2\nTOPL,3\nAMP,1,1.0,0,1.2\nAMP,1,1.1")
    m = ui.model
    assert m.damp == [0.02, 0.05] and m.freq_sets == {1: [2, 4, 6], 2: [9]} and m.topl == [1, 2, 3]
    assert m.amp == {1: [1.0, 0.0, 1.2, 1.1]}
    ui.run_text("DAMP,0\nFREQ,1,0\nTOPL,0\nAMP,1,0")
    assert m.damp == [] and m.freq_sets == {2: [9]} and m.topl == [] and m.amp == {}


OPTIONS_PRE = """N,1,0,0,0
N,2,1,0,0
INT,1,2,1,1
DYNP,1,0.0001,1.0,0.0001,0.4,Rock, hard
DYNP,2,0.001,1.0,0.0003,0.8,Rock, hard
SPRO,1,1,Rock, hard
SPRO,2,2,Sand
SACC,1,2,0
SRS,1,1,0
SSTR,2,1,0,1,0
SSAF,1,1,0,1,3,0.125,SAF top, vs rock
SFOU,1,0,0,0,0,0
DAMP,0.02,0.05
EQTIT,Design spectra 0.3g, 5% damping
EQUAKE,1,24,11975,0.05,20.0,1,5,1
RSIN,1,/tmp/spectra, set 1/H1.rsi
RSOUT,1,H1.rso
ACCIN,1,seed.acc
ACCOUT,1,H1.acc
TPSD,1,H1.psd
CORR,1,0.0,0.0
CORR,2,5.0,0.1
THFILE,/tmp/H1,x.acc
THTIT,acc X, 8192
SOIL,4000,32.2,0,1,1,8,0.65,1.0,0
SOILX,1,0,0.2,17,rock.acc
SITE,0,1,0,20,3,1,0,1,2048,1,0,0.005,4096,1
SITEX,1
TOPL,1,2,2,2
WAVE,2,1,0.5,0.5,0.0
WAVE,3,1,0.5,0.5,10.0
POINT,0,1,2.7
HOUSE,32.2,-5,0,2,1,1,1,1,1
HOUSEX,1,1,5,1,0
INCOH,0.1,0.1,0.2,0.5,1,1,0,1,7,9,180
WPASS,1.0E9,30,3
ME,1,1,2,0,0,0
AMP,1,1.0,0,0.5,0
SYMM,1,0,1,2,5
FORCE,1
ANALYS,0,0,1,1,1,0,30,1,2,3,2,1
ANALYSX,1,1
MOTION,0,0,0,25.0,0,0.1,100.0,301,1.0,0.0,1,4096,0,1,10,1,0,1,3
MOTIONX,1,2,0,1,0,1,0,0,1,0
NOUT,1,1,1,0,0,1,1,1-2
STRESS,0,1,1,1,4
STRESSX,1,0.5,2,1,1,0,0
EOUT,1,1,1,1,1,1,0,0,0,0,0,0,1,1
RELFILE,C:\\SSI\\00101TR_X.TFI
RDND,1,1,0,0,0,0,0
RELD,1,1,1
RELDX,1,1
CMODFORM,1
EDUOPT,HSLAW,UNIFORM
EDUOPT,LIMITS,PREP
AOPT,1,1,0,1,1,1,0,1,1,1,1,1,1,0
"""


def test_ut03_write_inp_round_trip_of_every_option_command(tmp_path):
    a = Interpreter(cwd=tmp_path)
    a.run_text(OPTIONS_PRE)
    assert not a.sink.texts(Kind.ERROR), a.sink.texts(Kind.ERROR)
    a.execute("WRITE,opts.pre")
    b = Interpreter(cwd=tmp_path)
    b.execute("INP,opts.pre")
    assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)
    assert a.model.same_state(b.model)
    text = (tmp_path / "opts.pre").read_text()
    for cmd in ("SITEX,1", "SOILX,1,0,0.2,17,rock.acc", "HOUSEX,1,1,5,1,0", "ANALYSX,1,1",
                "MOTIONX,1,2,0,1,0,1,0,0,1,0", "STRESSX,1,0.5,2,1,1,0,0", "RELDX,1,1", "CMODFORM,1",
                "WPASS,1.0E9,30,3", "SSAF,1,1,0,1,3,0.125,SAF top, vs rock"):
        assert cmd in text, cmd
    t1, _ = write_pre(b.model)
    t2, _ = write_pre(a.model)
    assert t1 == t2
    rb = b.model.options
    assert type(rb.record("MOTION")) is O.MotionRecord and rb.record("MOTION").interp == 3
    assert type(rb.entry("DYNP", ("Rock, hard", 2))) is O.DynpRecord


def test_default_x_commands_are_not_written(tmp_path):
    a = Interpreter(cwd=tmp_path)
    a.run_text("SITE,0,1\nSITEX,0\nMOTIONX,0,2\n")
    text, _ = write_pre(a.model)
    assert "SITEX" not in text and "MOTIONX" not in text and "SITE,0,1" in text


# ======================================================================================
# Catalogue placeholders come from the rule that detected the problem (spec 11 section 5)
# ======================================================================================
def test_problems_carry_the_catalogue_placeholder_values():
    site = make_record("SITE", ["0", "1", "0", "20", "2", "1", "0", "0", "-4"])
    p53 = [p for p in site.problems() if p[1] == 53]
    assert p53 == [("Error", 53, "<freq1> = 0"), ("Error", 53, "<freq2> = -4")]   # still plain 3-tuples
    assert [O.problem_fmt(p) for p in p53] == [{"i": 1}, {"i": 2}]
    wave = make_record("WAVE", ["3", "1", "1", "0", "400"])
    assert {(p[1], tuple(sorted(O.problem_fmt(p).items()))) for p in wave.problems()} == \
        {(52, (("w", 3),)), (54, (("i", 2), ("w", 3)))}
    assert O.problem_fmt(make_record("SSAF", ["4", "1", "0", "0", "2", "0"]).problems()[0]) == {"i": 4}
    sfou = make_record("SFOU", ["5", "1", "0", "0", "-1", "-2"]).problems()
    assert [(p[1], O.problem_fmt(p)) for p in sfou] == [(111, {"i": 5}), (112, {"i": 5})]
    assert O.problem_fmt(make_record("DYNP", ["2", "0.1", "1.5", "0.1", "1", "Clay"]).problems()[0]) == {"p": "Clay"}
    kind, num, det = site.problems()[0]                                     # unpacks like a tuple
    assert (kind, num) == ("Error", 53)


def test_command_time_warnings_name_the_frequency_index():
    ui = Interpreter()
    ui.execute("SITE,0,1,0,20,2,1,0,0,-4,1,0,0.01,4096,1")
    warns = ui.sink.texts(Kind.WARNING)
    assert any(w.startswith("SITE: Error 53 : Illegal Value for Frequency 1 (<freq1> = 0)") for w in warns), warns
    assert any(w.startswith("SITE: Error 53 : Illegal Value for Frequency 2 (<freq2> = -4)") for w in warns), warns
    ui.execute("DYNP,2,0.1,1.5,0.1,1,Clay")
    assert any("Error 98 : Dynamic property Clay has illegal shear modulus values" in w
               for w in ui.sink.texts(Kind.WARNING))


def test_check_names_the_layer_of_soil_request_errors():
    from sassi.prep.check import run_check
    ui = Interpreter()
    ui.run_text("SFOU,3,1,0,0,-1,0\nSSAF,2,1,0,0,1,0\n")
    rep = run_check(ui.model, modules=["SOIL"])
    texts = {m.text for m in rep.errors("SOIL")}
    assert "Illegal Number of Smoothings for Layer 3" in texts
    assert "Illegal Frequency Step for Layer 2" in texts


# ======================================================================================
# GUI parse = console parse (rule L17): the registered command specification is used
# ======================================================================================
def test_from_args_and_to_command_use_the_registered_specification():
    ui = Interpreter()
    line = "SOILX,0,0,0.1,0,a, b.acc"
    ui.execute(line)
    assert O.get_record(ui.model, "SOILX").file == "a, b.acc"
    o = O.from_args(line)
    assert o.file == "a, b.acc" and o.cl == 0
    assert O.to_command(o) == line
    assert O.to_command(O.get_record(ui.model, "SOILX")) == line
    assert O.from_args("soilx,1").indir == 1                              # name case-insensitive (L3)
