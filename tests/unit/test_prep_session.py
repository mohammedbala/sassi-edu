"""Session commands (requirements 3.4.A/B): models in memory, MDL/MDLNAME, CD/MKDIR, MOPT,
GRAVITY/GROUNDELEV, FREQ/LFREQ (UT-04), list and request commands (UT-04, D-PAR-14), STATUS."""
from __future__ import annotations

from pathlib import Path

import pytest

from sassi.prep import Interpreter, Kind


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def run(ui, text):
    ui.run_text(text)
    return ui.model


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


# ------------------------------------------------------------------ models in memory (D-MDL-01)
def test_actm_cpmodel_dmodel_modellist(ui):
    run(ui, "N,1,0,0,0\nACTM,2\nN,5,0,0,0\n")
    assert ui.active_model == 2 and set(ui.models) == {0, 2}
    assert set(ui.models[0].nodes) == {1} and set(ui.models[2].nodes) == {5}
    run(ui, "ACTM,0\nMDL,base,.\nCPMODEL,3\n")
    assert set(ui.models[3].nodes) == {1} and ui.models[3].name == "base"
    assert any("same name and path" in w for w in warnings(ui))
    ui.models[3].nodes[1].x = 99.0
    assert ui.models[0].nodes[1].x == 0.0                    # deep copy
    run(ui, "DMODEL,3\nDMODEL,0\nDMODEL,7\nCPMODEL,0\nACTM,-1\n")
    assert 3 not in ui.models and ui.models[0].is_empty() and ui.active_model == 0
    e = errors(ui)
    assert any("model 7 is not in memory" in t for t in e)
    assert any("destination is the active model" in t for t in e)
    assert any("model numbers are >= 0" in t for t in e)
    ui.sink.clear()
    run(ui, "MODELLIST")
    assert any(t.startswith("*    0") for t in ui.sink.texts(Kind.INFO))
    # Model > New: lowest unused number (D-MDL-01)
    assert ui.new_model() == 1 and ui.active_model == 1 and ui.new_model() == 3


def test_mdl_mdlname_title_and_cd_mkdir(ui, tmp_path):
    run(ui, "MKDIR,proj/a\nMKDIR,proj/a\nCD,proj\nCD,nothere\n")
    assert (tmp_path / "proj" / "a").is_dir()
    assert ui.cwd == (tmp_path / "proj").resolve()
    assert any("already exists" in w for w in warnings(ui))
    assert any("does not exist" in e for e in errors(ui))
    run(ui, "MDL,Demo5,a\n")
    m = ui.model
    assert m.name == "Demo5" and Path(m.path) == (tmp_path / "proj" / "a").resolve()
    assert ui.cwd == Path(m.path)                            # MDL also sets the working directory
    run(ui, "MDL,Demo6,new/dir\nMDLNAME,Renamed\nTIT,A title, with commas\n")
    assert Path(m.path) == (tmp_path / "proj" / "a" / "new" / "dir").resolve() and m.name == "Renamed"
    assert m.title == "A title, with commas"


# ------------------------------------------------------------------ global options
def test_mopt_record_and_validation(ui):
    m = run(ui, "MOPT,0,1,0,0\n")
    assert (m.mopt.incomp, m.mopt.matrix, m.mopt.mass, m.mopt.force) == (0, 1, 0, 0)
    run(ui, "MOPT,,,1\n")                                     # blank -> defaults (record replaced)
    assert (m.mopt.incomp, m.mopt.matrix, m.mopt.mass, m.mopt.force) == (1, 0, 1, 1)
    run(ui, "MOPT,2\n")
    assert any("must be 0 or 1" in e for e in errors(ui))


def test_gravity_groundelev_touch_only_their_house_argument(ui):
    m = run(ui, "HOUSE,32.2,0.0,0,2,1,0,0,0,0\nGRAVITY,9.81\nGROUNDELEV,-10\n")
    assert m.options.record("HOUSE").to_tokens() == ["9.81", "-10", "0", "2", "1", "0", "0", "0", "0"]
    assert m.gravity == 9.81 and m.ground_elevation == -10.0
    a = Interpreter()
    a.execute("GRAVITY,9.81")
    assert a.model.options.record("HOUSE").to_tokens() == ["9.81"]
    assert a.model.ground_elevation == 0.0


# ------------------------------------------------------------------ frequency sets (UT-04, UT-16 listing)
def test_ut04_freq_list_semantics(ui):
    m = run(ui, "FREQ,1,2,4,6\nFREQ,1,8,10\nFREQ,2,5\n")
    assert m.freq_sets == {1: [2, 4, 6, 8, 10], 2: [5]}
    run(ui, "FREQ,1,0\n")
    assert m.freq_sets == {2: [5]}
    run(ui, "FREQ,2,5\nFREQ,3,-4\n")
    assert any("duplicate frequency numbers in set 2: [5]" in w for w in warnings(ui))
    assert any("positive integers" in e for e in errors(ui))


def test_lfreq_lists_numbers_and_hz(ui):
    m = run(ui, "SITE,0,1,0.9765625\nFREQ,1,1,3,5,7,9,11,13,15,16,18\nLFREQ\n")
    nums, hz = m.frequencies(1)
    ref = [0.977, 2.930, 4.883, 6.836, 8.789, 10.742, 12.695, 14.648, 15.625, 17.578]
    assert nums == [1, 3, 5, 7, 9, 11, 13, 15, 16, 18]
    assert hz == pytest.approx(ref, abs=5e-4)
    info = ui.sink.texts(Kind.INFO)
    assert any("Frequency set 1: 10 frequencies, df = 0.9765625 Hz" in t for t in info)
    assert any("17.5781" in t for t in info)
    b = Interpreter()
    b.run_text("SITE,0,1,0,20,3,1,0,1,2048,1,0,0.005,4096,1\nFREQ,1,512")
    assert b.model.frequency_step() == 0.048828125 and b.model.frequencies(1)[1] == [25.0]


# ------------------------------------------------------------------ list commands (UT-04)
def test_ut04_damp_topl_amp(ui):
    m = run(ui, "DAMP,0.02\nDAMP,0.05\n")
    assert m.damp == [0.02, 0.05]
    run(ui, "DAMP,0\n")
    assert m.damp == []
    run(ui, "DAMP,0.02,0.05,0,0.07\nTOPL,1,2,3\nTOPL,4\n")
    assert m.damp == [0.02, 0.05, 0.07] and m.topl == [1, 2, 3, 4]
    run(ui, "TOPL,0\nTOPL,5,6\nAMP,1,1.0,0,0.5\nAMP,1,0.25\nAMP,2,2\n")
    assert m.topl == [5, 6]
    assert m.amp == {1: [1.0, 0.0, 0.5, 0.25], 2: [2.0]}      # AMP keeps zeros after the first value
    run(ui, "AMP,1,0\n")
    assert m.amp == {2: [2.0]}


# ------------------------------------------------------------------ requests (D-PAR-14)
def test_output_requests(ui):
    m = run(ui, "NOUT,1,1,1,0,0,1,1,1-5 7;9\nNOUT,3,1,1,0,0,1,1,1-3\n"
                "EOUT,1,1,1,1,1,1,2,0,0,0,0,0,3,1-20\nRDND,15,1,0,1,0,1,0\n")
    assert [(r.dir, r.nodes) for r in m.nout] == [(1, [1, 2, 3, 4, 5, 7, 9]), (3, [1, 2, 3])]
    assert m.nout[0].codes == [1, 1, 0, 0, 1, 1]
    assert m.eout[0].group == 3 and m.eout[0].elements == list(range(1, 21)) and m.eout[0].codes[6] == 2
    assert m.rdnd[0].node == 15 and m.rdnd[0].flags == [1, 0, 1, 0, 1, 0]
    run(ui, "EOUT,0,0,0,0,0,0,1,0,0,0,0,0,3,4\n")              # a code of 0 first is not a clear
    assert len(m.eout) == 2
    run(ui, "NOUT,0\nEOUT,0\nRDND,0\n")
    assert not m.nout and not m.eout and not m.rdnd
    run(ui, "NOUT,7,1,1,1,1,1,1,1\nNOUT,1,1,1,1,1,1,1\nNOUT,1,0,0,0,0,0,0,9-3\nEOUT,3,0,0,0,0,0,0,0,0,0,0,0,1,1\n")
    e = errors(ui)
    assert any("<dir> must be 1..6" in t for t in e)
    assert any("node list is empty" in t for t in e)
    assert any("descending range" in t for t in e)
    assert any("output codes are 0" in t for t in e)


def test_status(ui):
    run(ui, "MDL,st,.\nTIT,Status test\nN,1,0,0,0\nINT,1,1,1,1\nFREQ,1,2,4\nSYMM,1,0,1,2,3\nSITE,0\nSTATUS\n")
    info = "\n".join(ui.sink.texts(Kind.INFO))
    assert "Model 0: name 'st'" in info and "Title: Status test" in info
    assert "Nodes 1 (interaction 1)" in info and "Frequency sets (count): 1 (2)" in info
    assert "SYMM 1: 1,0,1,2,3" in info and "Stored options: SITE" in info


# ------------------------------------------------------------------ listings with unresolvable values
def test_lfreq_and_status_without_frequency_step(ui):
    """SITE delt = 0 and fstep = 0: df cannot be resolved; LFREQ/STATUS say so instead of failing."""
    run(ui, "SITE,0,1,0,,,,,,,,,0,4096,1\nFREQ,1,5,3\nLFREQ\nSTATUS")
    assert not any("internal error" in e for e in errors(ui))
    info = ui.sink.texts(Kind.INFO)
    assert "Frequency set 1: 2 frequencies, df unavailable" in info
    assert any(t.split() == ["3", "5"] for t in info)            # the numbers are still listed
    assert any(t.startswith("Frequency sets (count): 1 (2); df = unavailable") for t in info)
    assert any("frequency step unknown" in w for w in warnings(ui))


def test_lfreq_with_frequency_step_lists_hz(ui):
    run(ui, "SITE,0,1,0.5\nFREQ,1,4,2\nLFREQ")
    info = ui.sink.texts(Kind.INFO)
    assert "Frequency set 1: 2 frequencies, df = 0.5 Hz" in info
    assert any(t.split() == ["2", "1.0000", "4", "2.0000"] for t in info)


# ------------------------------------------------------------------ AMP keeps zeros after the first value
def test_amp_keeps_zero_ratios_and_round_trips(ui):
    """AMP: zeros after the first value are data (D-INC-10 (Re, Im) pairs, Error 118/119 count rule)."""
    from sassi.prep.writer import write_pre
    vals = [1.0] + [0.0] * 120 + [0.5, 0.0]          # a WRITE chunk of 100 values would start with 0
    m = run(ui, "AMP,1," + ",".join(str(v) for v in vals))
    assert m.amp[1] == vals
    run(ui, "AMP,1,0.7,0")                            # a later line appends (its first value is not 0)
    vals += [0.7, 0.0]
    assert m.amp[1] == vals
    run(ui, "AMP,2,1.0,0,1.02,0\nDAMP,0.02,0,0.05")
    assert m.amp[2] == [1.0, 0.0, 1.02, 0.0] and m.damp == [0.02, 0.05]       # DAMP drops zeros
    text, _ = write_pre(m)
    amp_lines = [ln.split(",") for ln in text.splitlines() if ln.startswith("AMP,")]
    clears = [t for t in amp_lines if len(t) == 3 and float(t[2]) == 0]
    data = [t for t in amp_lines if t not in clears]
    assert len(clears) == 2 and all(float(t[2]) != 0 for t in data)   # no data line starts with 0
    ui2 = Interpreter()
    ui2.run_text(text)
    assert ui2.model.amp == m.amp
    run(ui, "AMP,1,0")
    assert 1 not in m.amp
