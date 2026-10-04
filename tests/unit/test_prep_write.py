"""WRITE -> INP round trip (UT-03, D-PAR-17) and SAVE/RESUME (D-MDL-02, D-MDL-03)."""
from __future__ import annotations

import json

import numpy as np
import pytest

from sassi.io.container import write_container
from sassi.prep import Interpreter, Kind
from sassi.prep.writer import register_section, register_writer, write_pre, SECTIONS, WRITERS

RICH = r"""* rich model exercising every kind of stored datum
TIT,Demo surface mat, coherent - BS units
MOPT,1,1,0,1
GRAVITY,32.2
N,1,0,0,0
N,5,40,0,0
FILL
NGEN,4,5,1,5,1,0,10,0
INT,1,25,1,1,0
INT,3,7,2,1,2
D,1,12,1,1,ROTZ
D,14,25,1,1,ROTZ
D,3,3,1,1,UX,UZ,ROTY
LOC,1,0,10,0,0,45,0,0
LOC,3,0,1,2,3,10,20,30
N,30,0,0,5
N,31,0,1,5
N,32,-1,0,5
LOCAL,2,0,30,31,32
CSYS,2
N,40,1.1,2.2,3.3
CSYS,1
N,50,3,3,3
CSYS,3
N,60,0.1,0.2,0.3
M,1,519120,0.17,0.150,0.04,0.04,1
M,2,4000,1500,0.12,0.05,0.05,3
M,3,36000,12000,0.15,0.02,0.02,2
L,1,20,0.120,2000,1000,0.02,0.02
L,2,30,0.125,2500,1250,0.03,0.03
R,1,25,20.83,20.83,88.02,52.08,52.08
SC,1,1.0E5,1.0E5,2.0E5,0,0,0,0.05
MXR,1,1,1000,0,0,0,0,0,-1000
MXR,1,7,1000
MXI,1,1,40,0,0,0,0,0,-40
MXM,1,1,2.5
MXM,1,12,0
GROUP,1,SHELL
GTIT,1,Basemat, east part
MACT,1
E,1,1,2,7,6
EGEN,3,1,1
EGEN,3,5,1,4
THICK,1,8,1,2.0
THICK,9,16,1,1.5
ETYPE,5,8,1,2
MSET,10,11,1,3
N,100,20,20,30
N,102,-10,20,0
D,102,102,1,1,ALL
GROUP,2,BEAMS
MACT,1
RACT,1
E,1,13,100,102
E,2,13,102,100
KI,1,1,1,0,0,0,0,1,1
KJ,2,2,1,0,0,0,1,0,0
GROUP,3,SOLID
MACT,2
E,1,1,2,7,6,30,31,40,50
E,2,1,2,7,6,30,31,40,40
EINT,1,1,1,2
ETYPE,1,2,1,1
GROUP,4,SPRING
RACT,1
E,1,100,102
GROUP,5,GENERAL
E,1,13,100
GROUP,1
MACT,7
RACT,9
MUNITS,100,100,1,0
MT,100,500,500,500
MR,100,10,10,10
MT,7,1,2,3
F,50,1.0,0,0
F,51,0.5,0,0,0.1
MM,13,0,0,1,0,0,0.2
SYMM,1,0,1,2,5
SYMM,2,1,1,6,11
FREQ,1,0
FREQ,1,2,4,6,8,10,12,14,16,18,20
FREQ,1,24,28,32
FREQ,3,7
DYNP,1,0.0001,1.0,0.0001,0.4,Rock, hard
DYNP,2,0.001,1.0,0.0003,0.8,Rock, hard
DYNP,1,0.0001,1.0,0.0001,0.5,Sand
SPRO,1,1,Rock, hard
SPRO,2,2,Sand
SACC,1,2,0
SRS,1,1,0
SSTR,2,1,0,1,0
SSAF,1,1,0,1,3,10,SAF top, vs rock
DAMP,0
DAMP,0.02,0.05,0.07
EQTIT,Design spectra 0.3g, 5% damping
EQUAKE,0,24,11975,0.05,20.0,0,5,0
RSIN,1,/tmp/spectra, set 1/H1.rsi
RSOUT,1,/tmp/H1.rso
CORR,1,0.0,0.0
THFILE,/tmp/H1,x.acc
THTIT,acc X, 8192
SOIL,4000,32.2,0,1,1,8,0.65,1.0,(5X,F10.4)
SOILX,1,0,0.2
SITE,0,1,0,20,3,1,0,1,2048,1,0,0.005,4096,1
SITEX,1
TOPL,0
TOPL,1,2,2,2
WAVE,2,1,1.0,1.0,0.0
WAVE,4,1,1.0,1.0,0.0
POINT,0,0,2.7
HOUSE,32.2,0.0,0,2,0,0,0,0,0
GROUNDELEV,-5
HOUSEX,0,1,5
INCOH,0.1,0.1,0.2,0.5,1,1,0,1,0,0,0
WPASS,1.0E9,0.0,3
ME,1,1,69,0,0,0
AMP,1,1.0,0,0.5
AMP,2,0,1
ANALYS,0,0,0,1,1,0,0.0,0.0,0.0,0.0,0,0
ANALYSX,1
MOTION,0,0,0,25.0,0,0.1,100.0,301,1.0,0.0,1,4096,0,0,0,1,0,0,0
MOTIONX,0,2,0,1
NOUT,1,1,1,0,0,1,1,1-5 7;9
NOUT,3,1,1,0,0,1,1,1-5
STRESS,0,0,1,0,0
STRESSX,0,0.5
EOUT,1,1,1,1,1,1,1,0,0,0,0,0,1,1-16
EOUT,1,1,1,1,1,1,0,0,0,0,0,0,3,2
RELFILE,C:\SSI\00101TR_X.TFI
RDND,100,1,0,0,0,0,0
RDND,102,0,0,1,0,0,1
RELD,0,0,2
RELDX,1
EDUOPT,HSLAW,UNIFORM
EDUOPT,GEOMTOL,1e-5
CMODFORM,1
BINOUT,1,,0
SECDATAOPT,1
EQL,1,2,3
NLSLAYER,3,1,0.5
BBCP,2,3,0.1,5
AOPT,0,1,0,1,1,1,0,0,1,0,1,1,1,0
GROUP,2
CSYS,2
"""


def _build(tmp_path, text=RICH):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(text, name="rich.pre")
    return ui


def _replay(tmp_path, pre_name):
    ui = Interpreter(cwd=tmp_path)
    ui.execute(f"INP,{pre_name}")
    return ui


def _diff(a, b):
    ca, cb = a.canonical(), b.canonical()
    return {k: (ca[k], cb.get(k)) for k in ca if ca[k] != cb.get(k)}


def test_ut03_write_inp_round_trip_is_identical(tmp_path):
    a = _build(tmp_path)
    assert not a.sink.texts(Kind.ERROR), a.sink.texts(Kind.ERROR)
    a.execute("WRITE,rich_out.pre")
    b = _replay(tmp_path, "rich_out.pre")
    assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)
    assert _diff(a.model, b.model) == {}
    assert a.model.same_state(b.model)
    assert a.model.model_hash() == b.model.model_hash()
    m = b.model
    # spot checks of every UT-03 item
    assert m.title == "Demo surface mat, coherent - BS units"
    assert m.groups[1].title == "Basemat, east part"
    assert m.options.entry("DYNP", ("Rock, hard", 2)) is not None
    assert m.options.entry("SPRO", 1).values[-1] == "Rock, hard"
    assert m.options.entry("SSAF", 1).values[-1] == "SAF top, vs rock"
    assert m.options.string("THFILE") == "/tmp/H1,x.acc"
    assert m.options.string("EQTIT") == "Design spectra 0.3g, 5% damping"
    assert m.options.entry("RSIN", 1).values[-1] == "/tmp/spectra, set 1/H1.rsi"
    assert m.options.record("SOIL").legacy == "(5X,F10.4)"
    assert m.options.record("HOUSE").values[:2] == ["32.2", "-5"]
    assert m.options.record("SITEX").values == ["1"]
    assert sorted(k for k, _ in m.options.entries("SYMM")) == [1, 2]
    assert m.freq_sets == {1: [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32], 3: [7]}
    assert m.damp == [0.02, 0.05, 0.07] and m.topl == [1, 2, 2, 2]
    assert m.amp == {1: [1.0, 0.0, 0.5]}
    assert [r.nodes for r in m.nout] == [[1, 2, 3, 4, 5, 7, 9], [1, 2, 3, 4, 5]]
    assert [(r.group, r.elements) for r in m.eout] == [(1, list(range(1, 17))), (3, [2])]
    assert [r.node for r in m.rdnd] == [100, 102]
    assert (m.csys_active, m.group_active, m.mact, m.ract) == (2, 2, 7, 9)
    assert m.csys[2].kind == "LOCAL" and m.csys[3].kind == "LOC"
    assert m.groups[1].elements[10].mat == 3 and m.groups[1].elements[9].thick == 1.5
    assert m.mopt.matrix == 1 and m.mopt.mass == 0


def test_write_is_byte_stable(tmp_path):
    a = _build(tmp_path)
    t1, _ = write_pre(a.model, filename="x.pre")
    t2, _ = write_pre(a.model, filename="x.pre")
    assert t1 == t2
    (tmp_path / "x.pre").write_text(t1)
    b = _replay(tmp_path, "x.pre")
    t3, _ = write_pre(b.model, filename="x.pre")
    assert t3 == t1                                     # WRITE of the replayed model is identical


def test_write_order_follows_spec(tmp_path):
    a = _build(tmp_path)
    text, _ = write_pre(a.model)
    lines = [ln for ln in text.splitlines() if not ln.startswith("*")]
    first = {}
    for i, ln in enumerate(lines):
        first.setdefault(ln.split(",")[0], i)
    order = ["TIT", "MOPT", "LOC", "N", "LOCAL", "D", "INT", "M", "L", "R", "SC", "MXR", "GROUP", "E", "MT",
             "MUNITS", "F", "SYMM", "FREQ", "DYNP", "DAMP", "EQTIT", "EQUAKE", "THFILE", "SOIL", "SOILX", "SITE",
             "SITEX", "TOPL", "POINT", "HOUSE", "HOUSEX", "AMP", "ANALYS", "MOTION", "NOUT", "STRESS", "EOUT",
             "RELD", "RDND", "CMODFORM", "EQL", "AOPT"]
    idx = [first[k] for k in order]
    assert idx == sorted(idx)
    assert "NOUT,0" in lines and "EOUT,0" in lines and "RDND,0" in lines and "FREQ,1,0" in lines
    assert lines[-1].startswith("AOPT")


def test_write_defaults_and_options(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.execute("N,1,0,0,0")
    ui.execute("WRITE")
    assert any("no file name given" in e for e in ui.sink.texts(Kind.ERROR))
    ui.execute("MDL,Demo,mdir")
    ui.execute("WRITE")
    assert (tmp_path / "mdir" / "Demo.pre").exists()
    (tmp_path / "other").mkdir()
    ui.execute(f"WRITE,arch.pre,{tmp_path / 'other'}")
    assert (tmp_path / "other" / "arch.pre").exists()
    ui.execute("WRITE,x.pre,/nonexistent/dir/zz")
    assert any("does not exist" in e for e in ui.sink.texts(Kind.ERROR))
    ui.write_options["mdl"] = True
    ui.write_options["afwr"] = True
    text, _ = write_pre(ui.model, mdl=True, afwr=True)
    assert f"MDL,Demo,{ui.model.path}" in text and text.rstrip().endswith("AFWR")


def test_local_system_exact_when_defining_nodes_moved(tmp_path):
    """UT-03: a LOCAL system whose defining node moved is re-created by LOCAL itself (bit for bit)."""
    a = Interpreter(cwd=tmp_path)
    a.run_text("N,1,0.3,0.1,0\nN,2,0.1,1.7,0.2\nN,3,-1,0.4,0.5\nLOCAL,2,0,1,2,3\nCSYS,2\nN,10,1,2,3\n"
               "CSYS,0\nN,2,5,5,5\nD,2,2,1,1,UX\nF,2,1,0,0\n")
    text, notes = write_pre(a.model)
    assert notes == []
    lines = text.splitlines()
    assert "LOC,2" not in text and "LOCAL,2,0,1,2,3" in lines
    i = lines.index("LOCAL,2,0,1,2,3")
    assert lines[i + 1] == "N,2,5,5,5"                 # node 2 restored right after LOCAL
    (tmp_path / "f.pre").write_text(text)
    b = _replay(tmp_path, "f.pre")
    assert not b.sink.texts(Kind.ERROR)
    assert a.model.same_state(b.model)
    cs_a, cs_b = a.model.csys[2], b.model.csys[2]
    assert cs_b.kind == "LOCAL" and np.array_equal(cs_b.R, cs_a.R) and np.array_equal(cs_b.origin, cs_a.origin)
    assert np.array_equal(b.model.node_global(10), a.model.node_global(10))
    assert b.model.nodes[2].fix[0] == 1 and b.model.forces[2].factor[0] == 1.0
    assert write_pre(b.model)[0] == text


def test_local_system_helper_nodes_deleted_redefined_and_chained(tmp_path):
    """UT-03: helper nodes deleted, re-used in another system, and a LOCAL built on a LOCAL."""
    a = Interpreter(cwd=tmp_path)
    a.run_text("\n".join([
        "LOC,5,0,1,2,3,10,20,30",
        "N,901,0.3,0.1,0", "N,902,1.7,0.9,0.2", "N,903,-0.4,1.3,0.5", "LOCAL,1,0,901,902,903",
        "NDEL,901", "CSYS,5", "N,902,7,8,9",            # 901 deleted, 902 moved into system 5
        "CSYS,1", "N,903,0.5,0.5,0.5",                   # 903 redefined in the system it defined
        "N,1,1.1,2.2,3.3", "N,2,2.2,0.1,0", "N,3,0.2,3.1,0.4", "LOCAL,3,0,1,2,3",
        "CSYS,3", "N,50,0.25,-0.5,1.5", "NDEL,2", "MT,903,1,2,3", "D,903,903,1,1,ALL",
        "CSYS,1",
    ]))
    assert not a.sink.texts(Kind.ERROR)
    text, notes = write_pre(a.model)
    assert notes == [] and "LOC,1" not in text and "LOC,3" not in text
    (tmp_path / "g.pre").write_text(text)
    b = _replay(tmp_path, "g.pre")
    assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)
    assert a.model.same_state(b.model), _diff(a.model, b.model)
    for s in (1, 3, 5):
        assert np.array_equal(a.model.csys[s].R, b.model.csys[s].R)
        assert np.array_equal(a.model.csys[s].origin, b.model.csys[s].origin)
    ia, xa = a.model.global_coordinates()
    ib, xb = b.model.global_coordinates()
    assert np.array_equal(ia, ib) and np.array_equal(xa, xb)
    assert 901 not in b.model.nodes and 2 not in b.model.nodes
    assert b.model.tmass[903] == [1.0, 2.0, 3.0] and b.model.nodes[903].fix == [1] * 6
    assert write_pre(b.model)[0] == text


def test_local_record_without_defining_data_written_as_loc(tmp_path):
    """A LOCAL record without stored nodes/points (not produced by commands) falls back to LOC."""
    a = Interpreter(cwd=tmp_path)
    a.run_text("N,1,0,0,0\nN,2,0,1,0\nN,3,-1,0,0\nLOCAL,2,0,1,2,3\nCSYS,2\nN,10,1,2,3\n")
    cs = a.model.csys[2]
    cs.nodes, cs.points = (), []
    text, notes = write_pre(a.model)
    assert notes and "LOC" in notes[0]
    (tmp_path / "f.pre").write_text(text)
    b = _replay(tmp_path, "f.pre")
    assert np.allclose(b.model.csys[2].R, a.model.csys[2].R, atol=1e-14)
    assert np.allclose(b.model.node_global(10), a.model.node_global(10), atol=1e-12)
    a.execute("WRITE,f2.pre")
    assert any("exact to rounding" in w for w in a.sink.texts(Kind.WARNING))


def test_writer_hooks(tmp_path):
    a = Interpreter(cwd=tmp_path)
    a.run_text("SITE,0,1\n")
    a.model.extensions["demo"] = {"x": 1}
    register_writer("SITE", lambda m: ["SITE,typed"])
    register_section("Demo extension", lambda m: [f"* demo {m.extensions['demo']['x']}"])
    try:
        text, _ = write_pre(a.model)
        assert "SITE,typed" in text and "* Demo extension\n* demo 1" in text
    finally:
        WRITERS.pop("SITE", None)
        SECTIONS[:] = [s for s in SECTIONS if s[0] != "Demo extension"]


# ------------------------------------------------------------------ SAVE / RESUME
def test_save_resume_round_trip(tmp_path):
    ui = _build(tmp_path)
    ui.execute("SAVE")
    assert any("use MDL first" in e for e in ui.sink.texts(Kind.ERROR))
    ui.execute("MDL,rich,model")
    ui.model.ui_state["hidden_elements"] = [[1, 2]]
    before = json.dumps(ui.model.to_json(), sort_keys=True)
    ui.execute("SAVE")
    assert (tmp_path / "model" / "rich.sdb").exists()
    ui.run_text("NDEL,1,100\nGDEL,1,5\nTIT,changed\n")
    ui.execute("RESUME")
    assert not ui.sink.texts(Kind.ERROR)[1:]
    assert json.dumps(ui.model.to_json(), sort_keys=True) == before   # histories and GUI state included


def test_resume_refuses_unknown_versions(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.execute("MDL,v,.")
    ui.execute("RESUME")
    assert any("SAVE first" in e for e in ui.sink.texts(Kind.ERROR))
    write_container(tmp_path / "v.sdb", "SDB", {"model": np.asarray("{}")}, meta={"sdb_version": "9.0"})
    ui.execute("RESUME")
    assert any("database version 9.0" in e for e in ui.sink.texts(Kind.ERROR))


# ------------------------------------------------------------------ groups: MTYPE history, active group
def _rt(tmp_path, text):
    a = Interpreter(cwd=tmp_path)
    a.run_text(text)
    assert not a.sink.texts(Kind.ERROR), a.sink.texts(Kind.ERROR)
    out, _ = write_pre(a.model)
    b = Interpreter(cwd=tmp_path)
    b.run_text(out)
    return a, b, out


def test_mtype_kept_data_round_trips(tmp_path):
    """Spec 08 5.24: data MTYPE keeps (but the new type does not use) is written and read back."""
    a, b, out = _rt(tmp_path, "\n".join([
        "GROUP,1,SOLID", "E,1,1,2,3,4,5,6,7,8", "E,2,1,2,3,4", "EINT,1,1,1,2", "ETYPE,1,2,1,1",
        "MTYPE,1,TSHELL",                                    # 8-node element and EINT 2 kept
        "GROUP,2,SHELL", "E,1,1,2,3,4", "THICK,1,1,1,0.5", "MTYPE,2,PLANE",
        "GROUP,3,BEAMS", "E,1,1,2,3", "KI,1,1,1,1,0,0,0,0,1", "KJ,1,1,1,0,0,0,1,1,0", "MTYPE,3,SPRING",
        "GROUP,4,TSHELL", "E,1,1,2,3,4", "EINT,1,1,1,1", "THICK,1,1,1,0.25",   # no MTYPE: no switch
        "GROUP,2",
    ]))
    assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)
    assert a.model.same_state(b.model), _diff(a.model, b.model)
    assert [b.model.groups[g].type for g in (1, 2, 3, 4)] == [5, 4, 7, 5]
    assert b.model.groups[1].elements[1].nodes == [1, 2, 3, 4, 5, 6, 7, 8]
    assert b.model.groups[1].elements[1].eint == 2 and b.model.groups[2].elements[1].thick == 0.5
    assert b.model.groups[3].elements[1].ki == [1, 0, 0, 0, 0, 1]
    lines = out.splitlines()
    assert lines.index("GROUP,1,SOLID") < lines.index("EINT,1,1,1,2") < lines.index("MTYPE,1,TSHELL")
    i = lines.index("GROUP,2,PLANE")                 # E lines fit PLANE; THICK needs SHELL
    assert lines[i:].index("MTYPE,2,SHELL") < lines[i:].index("THICK,1,1,1,0.5") < lines[i:].index("MTYPE,2,PLANE")
    assert "GROUP,3,BEAMS" in lines and "MTYPE,3,SPRING" in lines
    assert "GROUP,4,TSHELL" in lines and not any(ln.startswith("MTYPE,4") for ln in lines)
    assert write_pre(b.model)[0] == out


def test_model_without_mtype_history_has_no_mtype_lines(tmp_path):
    a = _build(tmp_path)
    text, _ = write_pre(a.model)
    assert "MTYPE" not in text


def test_no_active_group_round_trips(tmp_path):
    """Spec 08 5.12: after GDEL of the active group no group is active, also after WRITE -> INP."""
    a, b, out = _rt(tmp_path, "GROUP,1,SOLID\nGROUP,3,SHELL\nE,1,1,2,3\nEOUT,1,0,0,0,0,0,0,0,0,0,0,0,4,1\n"
                              "GROUP,2,BEAMS\nGDEL,2")
    assert a.model.group_active is None and b.model.group_active is None
    assert sorted(b.model.groups) == [1, 3]
    assert [r.group for r in b.model.eout] == [4]      # the temporary group 4 did not drop this request
    assert a.model.same_state(b.model)
    assert "GROUP,4,SOLID\nGDEL,4" in out
    a2, b2, out2 = _rt(tmp_path, "GROUP,1,SOLID\nGDEL,1")
    assert "GDEL" not in out2 and b2.model.group_active is None and not b2.model.groups


def _random_lines(seed):
    """Random command sequences mixing coordinate systems, deletions, MTYPE and groups."""
    import random
    rng = random.Random(seed)
    nn = rng.randint(4, 12)
    out = [f"N,{i},{rng.uniform(-9, 9)!r},{rng.uniform(-9, 9)!r},{rng.uniform(-9, 9)!r}" for i in range(1, nn + 1)]
    types = ["SOLID", "BEAMS", "SHELL", "SPRING", "PLANE", "TSHELL", "GENERAL"]
    for _ in range(rng.randint(20, 60)):
        a = rng.randint(1, nn)
        b = rng.randint(a, nn)
        n3 = rng.sample(range(1, nn + 4), 3)
        out.append(rng.choice([
            f"LOC,{rng.randint(1, 4)},0,{rng.uniform(-3, 3)!r},0,1,{rng.uniform(0, 90)!r},{rng.uniform(0, 90)!r},7",
            f"LOCAL,{rng.randint(1, 4)},0,{n3[0]},{n3[1]},{n3[2]}", f"CSYS,{rng.randint(0, 4)}",
            f"SDEL,{rng.randint(1, 4)}", f"N,{rng.randint(1, nn + 3)},{rng.uniform(-9, 9)!r},0.5,{rng.uniform(-9, 9)!r}",
            f"NDEL,{a},{b}", f"GLOBAL,{a},{b}", f"D,{a},{b},1,1,ALL", f"INT,{a},{b},1,1,0", f"MT,{a},1,2,3",
            f"F,{a},{rng.uniform(-1, 1)!r},0,1,0.25", f"MUNITS,{a},{b},1,0",
            f"GROUP,{rng.randint(1, 3)},{rng.choice(types)}", f"GROUP,{rng.randint(1, 3)}",
            f"MTYPE,{rng.randint(1, 3)},{rng.choice(types)}", f"GDEL,{rng.randint(1, 3)}",
            "E," + ",".join(str(rng.randint(1, nn)) for _ in range(rng.randint(2, 9))),
            f"EINT,1,3,1,{rng.randint(0, 2)}", f"THICK,1,3,1,{rng.uniform(0.1, 2)!r}", f"ETYPE,1,3,1,{rng.randint(0, 2)}",
            f"KI,1,3,1,1,0,0,0,0,{rng.randint(0, 1)}", f"MACT,{rng.randint(1, 3)}", "ECOMPR",
        ]))
    return out


@pytest.mark.parametrize("seed", range(40))
def test_ut03_random_models_round_trip_bit_exact(tmp_path, seed):
    """UT-03 on random models: deep-equal state, bit-identical global coordinates, stable text."""
    a = Interpreter(cwd=tmp_path)
    a.run_text("\n".join(_random_lines(seed)))
    text, notes = write_pre(a.model)
    b = Interpreter(cwd=tmp_path)
    b.run_text(text)
    assert notes == [] and not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)
    assert a.model.same_state(b.model), _diff(a.model, b.model)
    ia, xa = a.model.global_coordinates()
    ib, xb = b.model.global_coordinates()
    assert np.array_equal(ia, ib) and np.array_equal(xa, xb)
    assert write_pre(b.model)[0] == text
