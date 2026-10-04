"""RUN<MODULE> (requirements sections 2.6, 3.4.G; D-RUN-01..05), FCOPY / FMOVE / REMOVEFREQ / VERIFY
(section 3.4.R), INTCOUNT and the memory estimate (UT-18), and the naming formulas (UT-17)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi import conventions as C
from sassi.io.container import read_container, write_container
from sassi.prep import Interpreter, Kind
from sassi.prep.check import dense_matrix_bytes
from sassi.prep.commands.modules_cmd import prerequisites
from sassi.prep.registry import lookup

SURFACE = """
MDL,run,{dir}
N,1,0,0,0
N,2,1,0,0
N,3,1,1,0
N,4,0,1,0
INT,1,4,1,1
D,1,4,1,1,ROTZ
M,1,3e6,0.25,0.15,0.05,0.05
GROUP,1,SHELL
E,1,1,2,3,4
THICK,1,1,1,0.5
L,1,10,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
TOPL,1
SITE,0,1,0,10,2,1,0,1,2048,1,0,0.005,4096,1
WAVE,2,1,1,1,0
FREQ,1,4,8
POINT,0,0,0.45
"""


def _ui(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(SURFACE.format(dir=tmp_path / "run"))
    return ui


def test_run_commands_are_registered():
    for mod in ("EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "COMBIN", "MOTION", "RELDISP",
                "STRESS"):
        spec = lookup(f"RUN{mod}")
        assert spec is not None and not spec.placeholder, mod
    assert not lookup("RUNNONLINEAR").placeholder                     # P2, implemented (Option NON)


def test_run_prerequisites(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    assert not ui.execute("RUNSITE")
    assert any("Model name/path not defined -- use MDL" in e for e in ui.sink.texts(Kind.ERROR))
    ui = _ui(tmp_path)
    assert not ui.execute("RUNSITE")
    assert any("Run AFWRITE first" in e for e in ui.sink.texts(Kind.ERROR))
    ui.execute("AOPT,0,0,0,0,1,0,0,0,0,0,0,0,0,0")
    ui.execute("AFWRITE")
    ui.sink.clear()
    assert not ui.execute("RUNPOINT")
    assert any("FILE2 missing -- run SITE Mode 1" in e for e in ui.sink.texts(Kind.ERROR))
    assert not ui.execute("RUNPOINT,7")
    assert any("model 7 is not in memory" in e for e in ui.sink.texts(Kind.ERROR))


def test_prerequisite_table(tmp_path):
    from sassi.io import decks
    d = decks.new("ANALYS")
    d["simul"] = 1
    msgs = prerequisites("ANALYS", "m", tmp_path, d)
    assert any("FILE1X missing" in s for s in msgs) and any("m.N4 missing -- run HOUSE" in s for s in msgs)
    d["type"], d["simul"] = 1, 3
    assert any("FILE9003 missing" in s for s in prerequisites("ANALYS", "m", tmp_path, d))
    assert prerequisites("HOUSE", "m", tmp_path, None) == ["HOUSE needs m.sit (run AFWRITE with SITE enabled)"]
    r = decks.new("RELDISP")
    r.table("rdnd").append([12, 1, 0, 0, 0, 0, 0])
    assert prerequisites("RELDISP", "m", tmp_path, r) == ["Run MOTION with Save Complex TF for nodes 12 TR_X"]
    assert len(prerequisites("COMBIN", "m", tmp_path, None)) == 2


def test_missing_module_reports_a_clear_error(tmp_path, monkeypatch):
    import sassi.prep.commands.modules_cmd as mc
    ui = _ui(tmp_path)
    monkeypatch.setattr(mc, "_module_available", lambda module: module != "POINT")
    assert not ui.execute("RUNPOINT")
    assert any("module POINT is not available in this build yet" in e for e in ui.sink.texts(Kind.ERROR))


def test_runsite_and_runpoint_end_to_end(tmp_path):
    ui = _ui(tmp_path)
    ui.execute("AOPT,0,0,0,1,1,0,0,0,0,0,0,0,0,0")
    ui.execute("AFWRITE")
    ui.sink.clear()
    assert ui.execute("RUNSITE"), ui.sink.texts(Kind.ERROR)
    assert ui.execute("RUNPOINT"), ui.sink.texts(Kind.ERROR)
    mdir = tmp_path / "run"
    for f in ("FILE1", "FILE2", "FILE3", "run_SITE.out", "run_POINT.out"):
        assert (mdir / f).exists(), f
    infos = ui.sink.texts(Kind.INFO)
    assert any("SITE finished with status OK" in t for t in infos)       # listing streamed (D-RUN-05)
    assert np.allclose(read_container(mdir / "FILE1")["fnum"], [4, 8])
    ui.execute("TIT,changed")                                             # title: not part of the hash
    ui.execute("POINT,0,0,0.5")
    ui.sink.clear()
    ui.execute("RUNPOINT")
    assert any("model changed since the last AFWRITE" in w for w in ui.sink.texts(Kind.WARNING))


def test_run_blocked_module(tmp_path):
    ui = _ui(tmp_path)
    ui.execute("AOPT,0,0,0,1,1,0,0,0,0,0,0,0,0,0")
    ui.execute("AFWRITE")
    ui.execute("SITE,0,1,0,3,2,1,0,1,2048,1,0,0.005,4096,1")
    ui.execute("AFWRITE")
    ui.sink.clear()
    assert not ui.execute("RUNSITE")
    assert any("run.sit not found" in e for e in ui.sink.texts(Kind.ERROR))   # renamed .bak (D-AFW-01)


def test_fcopy_fmove(tmp_path):
    ui = _ui(tmp_path)
    mdir = tmp_path / "run"
    (mdir / "FILE1").write_text("x")
    assert ui.execute("FCOPY,FILE1,FILE1X")
    assert (mdir / "FILE1X").read_text() == "x" and (mdir / "FILE1").exists()
    assert ui.execute("FMOVE,FILE1,FILE81")
    assert not (mdir / "FILE1").exists() and (mdir / "FILE81").exists()
    assert not ui.execute("FCOPY,nothere,FILE1Y")
    ui.execute("FCOPY,FILE1X,FILE81")
    assert any("FILE81 overwritten" in w for w in ui.sink.texts(Kind.WARNING))


def test_removefreq(tmp_path):
    ui = _ui(tmp_path)
    mdir = tmp_path / "run"
    fnum = np.array([2, 4, 6, 8])
    H = np.arange(8, dtype=complex).reshape(4, 2)
    write_container(mdir / "FILE8", "FILE8", {"fnum": fnum, "freq": fnum * 0.5, "eq_node": [1, 1],
                                              "eq_dof": [1, 2], "H": H},
                    meta={"df": 0.5, "type": 0, "case": "", "ang": 0, "cm": 0, "nfft": 4096, "delt": 0.005,
                          "model_hash": "x"}, module="ANALYS")
    assert ui.execute("REMOVEFREQ,FILE8,FILE8,4,8")
    c = read_container(mdir / "FILE8")
    assert c["fnum"].tolist() == [2, 6] and np.array_equal(c["H"], H[[0, 2]]) and c.meta["df"] == 0.5
    assert (mdir / "FILE8.bak").exists() and c["eq_node"].tolist() == [1, 1]


def test_verify_list_and_unknown(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    assert ui.execute("VERIFY,LIST")
    assert any(t.startswith("VP-52") for t in ui.sink.texts(Kind.INFO))
    assert not ui.execute("VERIFY,VP-999")


def test_verify_runs_a_problem(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    assert ui.execute("VERIFY,VP-52")
    infos = ui.sink.texts(Kind.INFO)
    assert any("VP-52 PASSED" in t for t in infos), infos[-20:]


def test_ut18_memory_estimate_and_intcount(tmp_path):
    assert dense_matrix_bytes(10000) == pytest.approx(14.4e9)
    assert dense_matrix_bytes(20000) == pytest.approx(57.6e9)
    ui = _ui(tmp_path)
    ui.execute("INTCOUNT")
    infos = ui.sink.texts(Kind.INFO)
    assert any(t.startswith("INTCOUNT: 4 interaction nodes (12 interaction DOFs)") for t in infos)


def test_ut17_naming_formulas():
    assert C.incoherent_file8_name(17, 2) == "FILE8050"
    assert C.incoherent_file8_name(50, 3) == "FILE8150"
    assert C.nodal_rs_name(12, 1, 1) == "00012TR_X01.RS"
    assert C.nodal_result_name(123456, 3, "TFI") == "123456TR_Z.TFI"
    assert C.nodal_result_name(12, 1, "TFU", max_node=100000) == "000012TR_X.TFU"
    assert C.element_result_name("BEAMS", 3, 45, "MXJ", "THS") == "BEAMS_003_00045_MXJ.THS"


# ======================================================================================
# ANALYS prerequisites: FILE4 family incl. DOFSMAP (section 2.1) and the frequency survey (2.6)
# ======================================================================================
def _fake(path, kind, fnum):
    write_container(path, kind, {"fnum": np.asarray(fnum), "freq": np.asarray(fnum) * 0.5}, meta={"df": 0.5})


def _analys_dir(tmp_path, f1, f3, files=("m.N4", "COOSK", "COOSM", "DOFSMAP")):
    for f in files:
        (tmp_path / f).write_text("x")
    _fake(tmp_path / "FILE1", "FILE1", f1)
    _fake(tmp_path / "FILE3", "FILE3", f3)


def _analys_deck(freqs, fopt=0):
    from sassi.io import decks
    d = decks.new("ANALYS")
    d["fopt"] = fopt
    for n in freqs:
        d.table("freqs").append([n])
    return d


def test_analys_needs_dofsmap(tmp_path):
    _analys_dir(tmp_path, [2, 4, 8], [2, 4, 8], files=("m.N4", "COOSK", "COOSM"))
    assert prerequisites("ANALYS", "m", tmp_path, _analys_deck([2, 4, 8])) == ["DOFSMAP missing -- run HOUSE"]


def test_analys_frequency_survey(tmp_path):
    _analys_dir(tmp_path, [2, 4, 8], [2, 4, 8])
    assert prerequisites("ANALYS", "m", tmp_path, _analys_deck([2, 4, 8])) == []
    assert prerequisites("ANALYS", "m", tmp_path, _analys_deck([4])) == []          # subset is fine
    msgs = prerequisites("ANALYS", "m", tmp_path, _analys_deck([2, 4, 6, 10]))
    assert len(msgs) == 2
    assert msgs[0].startswith("Frequency 6, 10 not in FILE1 -- re-run SITE")
    assert msgs[1].startswith("Frequency 6, 10 not in FILE3 -- re-run POINT")
    _fake(tmp_path / "FILE3", "FILE3", [2, 4])
    msgs = prerequisites("ANALYS", "m", tmp_path, _analys_deck([2, 4, 8]))
    assert msgs == [m for m in msgs if m.startswith("Frequency 8 not in FILE3 -- re-run POINT")] and msgs
    # <fopt> = 1: all FILE1 frequencies are solved, so all must be in FILE3
    msgs = prerequisites("ANALYS", "m", tmp_path, _analys_deck([], fopt=1))
    assert len(msgs) == 1 and msgs[0].startswith("Frequency 8 not in FILE3")


def test_analys_frequency_survey_vibration_and_simultaneous_cases(tmp_path):
    from sassi.io import decks
    _analys_dir(tmp_path, [2, 4], [2, 4, 8])
    for k, f in ((1, [2, 4, 8]), (2, [2, 4])):
        _fake(tmp_path / f"FILE9{k:03d}", "FILE9", f)
    d = _analys_deck([2, 4, 8])
    d["type"], d["simul"] = 1, 2
    msgs = prerequisites("ANALYS", "m", tmp_path, d)
    assert msgs == [m for m in msgs if "not in FILE9002 -- re-run FORCE" in m] and len(msgs) == 1
    assert msgs[0].startswith("Frequency 8 not in FILE9002")
    d2 = _analys_deck([2, 4, 8])
    d2["simul"] = 1                                                     # seismic X/Y/Z: FILE1X/Y/Z
    for name, f in (("FILE1X", [2, 4, 8]), ("FILE1Y", [2, 4, 8]), ("FILE1Z", [2, 8])):
        _fake(tmp_path / name, "FILE1", f)
    msgs = prerequisites("ANALYS", "m", tmp_path, d2)
    assert len(msgs) == 1 and msgs[0].startswith("Frequency 4 not in FILE1Z -- re-run SITE")
    assert isinstance(decks.new("ANALYS"), decks.Deck)


def test_runanalys_stops_on_the_frequency_survey(tmp_path, monkeypatch):
    """RUNANALYS checks the prerequisites before launching (the ANALYS module is not needed)."""
    import sassi.prep.commands.modules_cmd as mc
    ui = _ui(tmp_path)
    ui.execute("AOPT,0,0,0,0,0,0,0,0,1,0,0,0,0,0")                    # ANALYS only
    assert ui.execute("AFWRITE")
    mdir = tmp_path / "run"
    assert (mdir / "run.anl").exists()
    _analys_dir(mdir, [4, 8], [4], files=("run.N4", "COOSK", "COOSM", "DOFSMAP"))
    monkeypatch.setattr(mc, "_module_available", lambda module: True)
    ui.sink.clear()
    assert not ui.execute("RUNANALYS")
    assert any("RUNANALYS: Frequency 8 not in FILE3 -- re-run POINT" in e for e in ui.sink.texts(Kind.ERROR))


# ======================================================================================
# INTCOUNT honours EDUOPT,MEMLIMIT like CHECK EDU-20 (D-ANL-10)
# ======================================================================================
def test_intcount_memory_guard_uses_memlimit(tmp_path, monkeypatch):
    import sassi.prep.commands.checks as ck
    from sassi.prep.check import analys_memory_bytes
    need = analys_memory_bytes(4)
    monkeypatch.setattr(ck, "physical_ram", lambda: int(need / 0.5))       # estimate = 50 % of the RAM
    ui = _ui(tmp_path)
    ui.execute("INTCOUNT")
    assert not any("memory estimate exceeds" in w for w in ui.sink.texts(Kind.WARNING))   # default 0.8
    ui.execute("EDUOPT,MEMLIMIT,4D-1")
    ui.sink.clear()
    ui.execute("INTCOUNT")
    assert any("exceeds 40 % of the physical memory (EDUOPT,MEMLIMIT = 0.4; EDU-20)" in w
               for w in ui.sink.texts(Kind.WARNING)), ui.sink.texts(Kind.WARNING)


def test_verify_runs_registered_problems_when_another_vp_module_is_broken(tmp_path, monkeypatch):
    import sassi.verify as V
    import sassi.verify.problems.vp_checks  # noqa: F401

    def broken():
        raise ImportError("simulated")

    monkeypatch.setattr(V, "load_all", broken)
    monkeypatch.setattr(V, "run_problem", lambda *a, **k: (_ for _ in ()).throw(AssertionError("not used")))
    ui = Interpreter(cwd=tmp_path)
    assert ui.execute("VERIFY,VP-52")
    assert any("VP-52 PASSED" in t for t in ui.sink.texts(Kind.INFO))
