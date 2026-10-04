"""Option A commands (sassi/prep/commands/loadgen_cmds.py): LOADGEN, LOADGENDYN, LGFILE, LGNODE, LGTIME,
LGMAP, LGOPT, LGLIST, RUNLOADGEN -- validation, the deck they produce, WRITE -> INP (UT-03) and a run of
the module through the interpreter (D-LGN-01)."""
from __future__ import annotations

import shutil

import numpy as np
import pytest

from sassi.modules import loadgen as LG
from sassi.prep import Interpreter
from sassi.prep.commands import loadgen_cmds as LC
from sassi.prep.messages import Kind
from sassi.prep.writer import write_pre
from sassi.verify import builders as B
from sassi.verify.problems import vp_loadgen as V

NODES = "N,1,0,0,0\nN,2,1,0,0\nN,3,2,0,0\nN,10,0,0,4\nN,11,0,0,8\nN,12,0,0,12\n"


def _ui(tmp_path, text=NODES):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(text)
    return ui


def _ok(ui, line):
    ok = ui.execute(line)
    return ok, ui.sink.texts(Kind.ERROR)[-1:] if not ok else []


# ======================================================================================
# records and validation
# ======================================================================================
def test_loadgen_record(tmp_path):
    ui = _ui(tmp_path)
    assert ui.execute("LOADGEN,dispacc,1,0,1,master,1,file8")
    rec = ui.model.options.record("LOADGEN")
    assert rec.to_tokens() == ["DISPACC", "1", "0", "1", "MASTER", "1", "FILE8"]
    assert rec.data == "DISPACC" and rec.multi == 1
    for bad in ("LOADGEN,9", "LOADGEN,2,2", "LOADGEN,2,0,0,0,HEAVY", "LOADGEN,2,0,0,0,1,0,DISK", "LOADGEN,2,x"):
        assert not ui.execute(bad), bad
    assert ui.model.options.record("LOADGEN").to_tokens()[0] == "DISPACC"          # failed commands store nothing


def test_loadgendyn_record_and_rayleigh(tmp_path):
    ui = _ui(tmp_path)
    assert ui.execute("LOADGENDYN,0.5,0.003")
    assert ui.model.options.record("LOADGENDYN").to_tokens() == ["0.5", "0.003"]
    assert ui.execute("LOADGENDYN,,,rel,0,results,1,0,0.05,1,10")
    assert any("Rayleigh alpha" in t for t in ui.sink.texts(Kind.INFO))
    a, b = LC.rayleigh(0.05, 1.0, 10.0)
    w1, w2 = 2 * np.pi, 20 * np.pi
    assert a / (2 * w1) + b * w1 / 2 == pytest.approx(0.05) and a / (2 * w2) + b * w2 / 2 == pytest.approx(0.05)
    d = LC.build_deck(ui.model, "DYNAMIC")
    assert (d["alpha"], d["beta"]) == pytest.approx((a, b)) and d["rotdisp"] == 1 and d["method"] == "REL"
    for bad in ("LOADGENDYN,-1", "LOADGENDYN,0,0,ACC", "LOADGENDYN,0,0,XYZ", "LOADGENDYN,0,0,REL,0,RESULTS,0,0,0.05,5,1",
                "LOADGENDYN,0,0,REL,0,RESULTS,0,0,2", "LOADGENDYN,0,0,REL,0,RESULTS,0,0,0,0,0,3"):
        assert not ui.execute(bad), bad
    assert ui.execute("LOADGENDYN,0,0,ACC,5")
    assert LC.build_deck(ui.model, "DYNAMIC")["refnode"] == 5


def test_lgfile_lgmap_lgopt(tmp_path):
    ui = _ui(tmp_path)
    assert ui.execute("LGFILE,apdl,out dir/my model, static.inp")
    assert ui.model.options.entry("LGFILE", "APDL").arg(2) == "out dir/my model, static.inp"
    assert ui.execute("LGFILE,APDLDYN,dyn.inp") and ui.execute("LGFILE,GROUND,kin.acc")
    assert LC.build_deck(ui.model, "STATIC")["apdlfile"] == "out dir/my model, static.inp"
    assert LC.build_deck(ui.model, "DYNAMIC")["apdlfile"] == "dyn.inp"
    assert ui.execute("LGFILE,APDL,")
    assert ui.model.options.entry("LGFILE", "APDL") is None
    assert not ui.execute("LGFILE,NOPE,x")
    assert not ui.execute("LGMAP,PAIRS")
    assert ui.execute("LGMAP,COORD,0.001,ansys model.cdb")
    d = LC.build_deck(ui.model, "STATIC")
    assert (d["mapmode"], d["maptol"], d["mapfile"]) == ("COORD", 0.001, "ansys model.cdb")
    assert ui.execute("LGMAP") and ui.model.options.record("LGMAP") is None
    assert ui.execute("LGOPT,12,0") and ui.model.options.record("LGOPT") is None        # defaults: not stored
    assert ui.execute("LGOPT,8,1")
    d = LC.build_deck(ui.model, "STATIC")
    assert d["digits"] == 8 and d["opmode"] == 1
    assert not ui.execute("LGOPT,3") and not ui.execute("LGOPT,12,2") and not ui.execute("LGOPT,12,0,2")
    assert ui.execute("LGOPT,12,0,0") and LC.build_deck(ui.model, "STATIC")["rest"] == 0
    assert ui.execute("LGOPT,12,0,1") and ui.model.options.record("LGOPT") is None


def test_lgnode_lists(tmp_path):
    ui = _ui(tmp_path)
    assert ui.execute("LGNODE,D,1-3") and ui.execute("LGNODE,d,3,10")
    assert LC.lgnode_lists(ui.model) == {"D": [1, 2, 3, 10]}
    assert ui.execute("LGNODE,M,10,11") and ui.execute("LGNODE,A,12")
    d = LC.build_deck(ui.model, "STATIC")
    assert d.tables == {"dnodes": [1, 2, 3, 10], "masters": [10, 11], "checknodes": [12]}
    assert ui.execute("LGNODE,D,99")
    assert any("not in the model" in t for t in ui.sink.texts(Kind.WARNING))
    assert ui.execute("LGNODE,D,0") and "D" not in LC.lgnode_lists(ui.model)
    for bad in ("LGNODE,X,1", "LGNODE,D", "LGNODE,D,0,1"):
        assert not ui.execute(bad), bad


def test_lgtime_variants(tmp_path):
    ui = _ui(tmp_path)
    good = {"LGTIME,vx,3,0.5": ("VX", 3, 0.5), "LGTIME,V": ("V", 1, 0.0), "LGTIME,MY,2,0,1,2,3": ("MY", 2, 0.0),
            "LGTIME,ACC,1,0,12,1": ("ACC", 1, 0.0), "LGTIME,DISP,2,1,1,3": ("DISP", 2, 1.0)}
    for line, (crit, n, tsep) in good.items():
        assert ui.execute(line), line
        d = LC.build_deck(ui.model, "STATIC")
        assert (d["crit"], d["ncrit"], d["tsep"]) == (crit, n, tsep), line
    assert d["critnode"] == 1 and d["critdof"] == 3
    assert ui.execute("LGTIME,MY,2,0,1,2,3") and LC.build_deck(ui.model, "STATIC")["refpoint"] == [1.0, 2.0, 3.0]
    assert ui.execute("LGTIME,TIME,1.5,2.25") and LC.build_deck(ui.model, "STATIC")["times"] == [1.5, 2.25]
    assert ui.execute("LGTIME,STEP,10,20") and LC.build_deck(ui.model, "STATIC")["times"] == [10.0, 20.0]
    for bad in ("LGTIME,Q", "LGTIME,TIME", "LGTIME,TIME,-1", "LGTIME,STEP,0", "LGTIME,V,0", "LGTIME,V,1,-1",
                "LGTIME,ACC,1,0,12", "LGTIME,ACC,1,0,12,7", "LGTIME,MY,1,0,1,2", "LGTIME,V,1,0,5"):
        assert not ui.execute(bad), bad


def test_lglist(tmp_path):
    ui = _ui(tmp_path, "MDL,m,m\n" + NODES)
    ui.run_text("LOADGEN,3,1\nLGTIME,V,2,0.5\nLGNODE,A,12\nLGFILE,HOUSE,other.hou\nLOADGENDYN,0.2,0.001,REL,0,FILE8")
    assert ui.execute("LGLIST")
    text = ui.sink.texts(Kind.INFO)[-1]
    assert "Disp. and Accel." in text and "2 largest peak(s)" in text and "housefile  = other.hou" in text
    assert "acceleration check nodes: 12" in text and "default (interaction nodes)" in text
    assert ui.execute("LGLIST,DYNAMIC")
    text = ui.sink.texts(Kind.INFO)[-1]
    assert "method REL" in text and "source FILE8" in text and "alpha 0.2" in text


# ======================================================================================
# WRITE -> INP
# ======================================================================================
def test_write_inp_round_trip(tmp_path):
    ui = _ui(tmp_path, "MDL,m,m\n" + NODES)
    long = ",".join(str(n) for n in range(1, 200, 2))
    ui.run_text("N,199,5,5,5\n" + "\n".join(f"N,{n},{n},1,0" for n in range(5, 199, 2)) + "\n")
    lines = ["LOADGEN,ACC,0,1,1,LUMPED,1,FILE8", "LOADGENDYN,,,ACC,10,RESULTS,1,1,0.03,0.5,20,1,1.5",
             "LGFILE,LUMPED,masses dir/my.masl", "LGFILE,SSIPATH,../results", "LGNODE,D,1-3",
             f"LGNODE,A,{long}", "LGNODE,M,10-12", "LGTIME,MX,2,0.25,0,0,-1", "LGMAP,PAIRS,0,map file.txt",
             "LGOPT,10,1"]
    for line in lines:
        assert ui.execute(line), (line, ui.sink.texts(Kind.ERROR)[-1:])
    text, notes = write_pre(ui.model)
    lg = [ln for ln in text.splitlines() if ln.startswith(("LOADGEN", "LG"))]
    assert "LGNODE,A,0" in lg and sum(ln.startswith("LGNODE,A,") for ln in lg) >= 4     # long list split
    assert all(len(ln) < 300 for ln in lg)
    ui2 = Interpreter(cwd=tmp_path)
    ui2.run_text(text)
    assert not ui2.sink.texts(Kind.ERROR)
    assert ui.model.same_state(ui2.model)
    assert write_pre(ui2.model)[0].splitlines()[4:] == text.splitlines()[4:]           # byte-stable (banner apart)
    a, b = LC.build_deck(ui2.model, "DYNAMIC").params, LC.build_deck(ui.model, "DYNAMIC").params
    assert {k: v for k, v in a.items() if k != "model"} == {k: v for k, v in b.items() if k != "model"}


# ======================================================================================
# RUNLOADGEN
# ======================================================================================
@pytest.fixture(scope="module")
def results(tmp_path_factory):
    """A model directory with HOUSE, FILE8, MOTION and RELDISP results (the vp_loadgen spring model)."""
    wd = tmp_path_factory.mktemp("lgcmd")
    fs, mdl = V.spring_model(wd, B, nfft=512, dt=0.02)
    V.control_motion(wd, n=300, seed=3)
    nodes = [(n, (1, 2, 3)) for n in (5, 12, 14, 16)]
    rc = V.run_motion_reldisp(wd, B, fs, nodes, range(1, 10))
    assert rc == (0, 0)
    yield wd
    shutil.rmtree(wd, ignore_errors=True)


def test_runloadgen_errors(tmp_path):
    ui = _ui(tmp_path)
    assert not ui.execute("RUNLOADGEN")
    assert "MDL" in ui.sink.texts(Kind.ERROR)[-1]
    ui.run_text("MDL,m,m")
    assert not ui.execute("RUNLOADGEN,STATIC")
    assert "FILE4" in " ".join(ui.sink.texts(Kind.ERROR)[-2:])
    assert not ui.execute("RUNLOADGEN,SOMETIMES")
    assert not ui.execute("RUNLOADGEN,STATIC,7")
    assert not ui.execute("RUNLOADGEN,7,STATIC") and "model 7" in ui.sink.texts(Kind.ERROR)[-1]
    assert not ui.execute("RUNLOADGEN,0,SOMETIMES") and "argument 2" in ui.sink.texts(Kind.ERROR)[-1]


def test_runloadgen_static_and_dynamic(results, tmp_path):
    md = tmp_path / "m"
    shutil.copytree(results, md)
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(f"MDL,m,{md}\nLOADGEN,DISPACC,0,0,0,LUMPED,1\nLGTIME,V,1\nLGNODE,A,12,14,16\n"
                "LOADGENDYN,0.4,0.002\nLGFILE,APDLDYN,dyn.inp")
    assert ui.execute("RUNLOADGEN"), ui.sink.texts(Kind.ERROR)
    assert (md / "m.lgn").exists() and (md / "m_LGS.inp").exists() and (md / "m.masl").exists()
    d, _ = LG.read_deck(md / "m.lgn")
    assert d["analysis"] == "STATIC" and d["data"] == 3 and d["genmass"] == 1
    assert any("LOADGEN finished with status OK" in t for t in ui.sink.texts(Kind.INFO))
    prog = V.parse_file(md / "m_LGS.inp")
    assert prog.errors == [] and len(prog.loads("F")) == 9 and len(prog.loads("D")) == 27
    assert ui.execute("RUNLOADGEN,0,DYNAMIC"), ui.sink.texts(Kind.ERROR)       # RUN<MODULE> argument order
    prog = V.parse_file(md / "dyn.inp")
    assert prog.errors == [] and {"LGA_12_UX", "LGA_14_UX", "LGA_16_UX", "LG_ACX"} <= set(prog.arrays)
    assert prog.value(prog.find("ALPHAD")[0].fields[0]) == 0.4
    # a failing run reports FAILED and the listing name
    ui.run_text("LGFILE,GROUND,missing.acc")
    assert not ui.execute("RUNLOADGEN,DYNAMIC")
    assert "m_LOADGEN.out" in ui.sink.texts(Kind.ERROR)[-1]
    ui.run_text("LGFILE,GROUND,\nLOADGEN,ACC,0,0,0,1,0,FILE8")
    (md / "m.mot").unlink()
    assert not ui.execute("RUNLOADGEN,STATIC")
    assert "m.mot" in " ".join(ui.sink.texts(Kind.ERROR)[-2:])
