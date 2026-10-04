"""File-conversion commands: CONVERT (ANSYS, SSI, STRUDL), ANSYS, ANSYSREFORMAT, ANSYSMODELTYPE,
GENMATRIXDAMP (requirements 3.4.K; D-ANS-01..08, D-MDL-03; spec 09 section 6)."""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from sassi.io import decks
from sassi.prep import Interpreter, Kind

DATA = Path(__file__).resolve().parents[1] / "data" / "ansys"


def _ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def _msgs(ui, kind):
    return ui.sink.texts(kind)


def _any(ui, kind, frag):
    return any(frag in t for t in _msgs(ui, kind))


# ======================================================================================
# CONVERT,ANSYS
# ======================================================================================
def test_convert_ansys_into_another_model(tmp_path):
    ui = _ui(tmp_path)
    ui.execute("N,1,0,0,0")
    ui.execute("ACTM,2")
    ui.execute(f"MDL,cant,{tmp_path}")
    ui.execute("ACTM,0")
    assert ui.execute(f"CONVERT,ANSYS,2,{DATA / 'beam188_cantilever.cdb'},9.81,cant_conv")
    assert ui.active_model == 0 and list(ui.model.nodes) == [1]           # active model unchanged
    m2 = ui.models[2]
    assert (m2.name, m2.path) == ("cant", str(tmp_path))                  # MDL identity kept
    assert len(m2.nodes) == 22 and m2.groups[1].type == 2 and m2.gravity == 9.81
    mp = tmp_path / "cant_cdb.map"                                        # D-ANS-02 map file
    rows = [ln.split() for ln in mp.read_text().splitlines() if not ln.startswith("!")]
    assert rows[0] == ["1", "1", "1"] and len(rows) == 20
    # the optional .pre output reads back to the same model state (UT-03 style)
    pre = tmp_path / "cant_conv.pre"
    ui2 = _ui(tmp_path)
    ui2.run_file(pre)
    assert not ui2.sink.texts(Kind.ERROR)
    assert ui2.model.same_state(m2)
    assert _any(ui, Kind.INFO, "only a subset of ANSYS elements")


def test_convert_ansys_replaces_a_non_empty_model_and_default_damping(tmp_path):
    ui = _ui(tmp_path)
    ui.execute("N,99,1,2,3")
    assert ui.execute(f"CONVERT,ANSYS,,{DATA / 'beam188_default_orientation.cdb'},9.81,,0.04")
    assert 99 not in ui.model.nodes and _any(ui, Kind.WARNING, "model 0 replaced")
    assert ui.model.materials[1].sdamp == 0.04 and ui.model.materials[1].pdamp == 0.04
    assert (tmp_path / "beam188_default_orientation_cdb.map").exists()     # no MDL: working directory
    assert not list(DATA.glob("*.map"))                                     # nothing written next to the .cdb


def test_convert_argument_errors(tmp_path):
    ui = _ui(tmp_path)
    cdb = DATA / "beam188_cantilever.cdb"
    assert not ui.execute(f"CONVERT,ANSYS,,{cdb}")                         # gravity required
    assert _any(ui, Kind.ERROR, "<gravity> is required")
    assert not ui.execute(f"CONVERT,ANSYS,,{cdb},0")
    assert not ui.execute(f"CONVERT,ANSYS,,{cdb},9.81,,0.7")               # damping ratio range
    assert not ui.execute(f"CONVERT,ANSYS,,{tmp_path / 'missing.cdb'},9.81")
    assert not ui.execute("CONVERT,NASTRAN,,x.bdf")
    assert not ui.execute(f"CONVERT,ANSYS,-1,{cdb},9.81")
    (tmp_path / "empty.cdb").write_text("/PREP7\nFINISH\n")
    assert not ui.execute(f"CONVERT,ANSYS,,{tmp_path / 'empty.cdb'},9.81")
    assert ui.execute("CONVERT,STRUDL,1,joints.txt") and _any(ui, Kind.WARNING, "not applicable in this version")
    assert ui.execute("CONVERT") and _any(ui, Kind.INFO, "CONVERT,ANSYS,<model>,<file.cdb>,<gravity>")


# ======================================================================================
# CONVERT,SSI (D-ANS-08)
# ======================================================================================
SSI_MODEL = """MDL,ssi,{dir}
TIT,CONVERT SSI round trip
LOC,1,0,100,0,0,0,0,0
CSYS,1
N,1,0,0,0
N,2,10,0,0
N,3,10,10,0
N,4,0,10,0
N,5,0,0,-5
N,6,10,0,-5
N,7,10,10,-5
N,8,0,10,-5
CSYS,0
N,10,100,0,5
N,11,100,0,10
N,12,105,0,5
N,13,110,0,0
N,14,110,0,-1
N,15,105,0,10
N,16,105,5,10
N,17,105,5,15
INT,1,8,1,1
M,1,4.32e5,0.25,0.15,0.05,0.05,1
M,2,3e6,0.2,0.16,0.07,0.07,1
R,1,1,0.8,0.8,0.2,0.1,0.1
SC,1,100,0,200,0,0,5,0.02
MXR,1,1,1000,0,0,0,0,0,-1000
MXR,1,7,1000
MXI,1,1,20
MXM,1,1,2
L,1,5,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
L,3,0,0.14,3000,1500,0.01,0.01
TOPL,1,2
GROUP,1,SOLID
E,1,5,6,7,8,1,2,3,4
ETYPE,1,1,1,2
GROUP,2,BEAMS
E,1,10,11,12
E,2,11,15,12
KJ,2,2,1,0,0,0,0,1,0
GROUP,3,SHELL
MACT,2
E,1,1,2,3,4
THICK,1,1,1,0.4
GROUP,4,SPRING
E,1,13,14
GROUP,5,GENERAL
E,1,15,16,17
MT,11,32.2,32.2,32.2
MR,15,1,2,3
MUNITS,15,15,1,0
D,14,14,1,1,ALL
FREQ,1,2,4,8,16
SITE,0,1,0,20,3,1,0,1,2048,1,0,0.005,4096,1
WAVE,2,1,0.5,0.5,0
WAVE,1,1,0.5,0.5,0
POINT,0,1,4.5
HOUSE,32.2,0,0,2,0,0,0,0,0
MOPT,0
AOPT,0,0,0,1,1,1,0,0,0,0,0,0,0,0
AFWRITE
"""


def _read(path, module):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        return decks.read(path, module)


def test_convert_ssi_round_trip(tmp_path):
    d1 = tmp_path / "orig"
    d1.mkdir()
    ui = _ui(tmp_path)
    ui.run_text(SSI_MODEL.format(dir=d1))
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    assert {p.name for p in d1.iterdir()} >= {"ssi.hou", "ssi.sit", "ssi.poi"}
    m0 = ui.models[0]
    assert ui.execute(f"CONVERT,SSI,1,{d1 / 'ssi.hou'},{tmp_path / 'back.pre'}")
    assert _any(ui, Kind.INFO, "ssi.sit read") and _any(ui, Kind.INFO, "ssi.poi read")
    m1 = ui.models[1]
    converted = m1.copy()
    # nodes: global coordinates (the local system is resolved), fixities, interaction flags
    assert sorted(m1.nodes) == sorted(m0.nodes)
    for n in m0.nodes:
        np.testing.assert_allclose(m1.node_global(n), m0.node_global(n), atol=1e-12)
        assert m1.nodes[n].fix == m0.nodes[n].fix, n
        assert m1.nodes[n].flags == m0.nodes[n].flags
    # elements and tables
    for g, e in m0.iter_elements():
        e1 = m1.groups[g.id].elements[e.id]
        assert (e1.nodes, e1.mat, e1.prop, e1.eint, e1.thick, e1.ki, e1.kj) == \
               (e.nodes, e.mat, e.prop, e.eint, e.thick, e.ki, e.kj)
        assert m1.groups[g.id].type == g.type
    assert m1.groups[1].elements[1].etype == 2 and m1.groups[2].elements[1].etype == 0
    for tab in ("materials", "layers", "sections", "springs"):
        assert {k: v.to_json() for k, v in getattr(m1, tab).items()} == \
               {k: v.to_json() for k, v in getattr(m0, tab).items()}, tab
    assert m1.matrices[1].to_json() == m0.matrices[1].to_json()
    assert (m1.tmass, m1.rmass, m1.mass_units) == (m0.tmass, m0.rmass, m0.mass_units)
    # options
    for name in ("SITE", "POINT", "HOUSE", "MOPT"):
        a, b = m0.options.record(name), m1.options.record(name)
        assert [float(x) for x in a.to_tokens() if x] == [float(x) for x in b.to_tokens() if x], name
    assert m1.topl == m0.topl and m1.freq_sets == m0.freq_sets
    assert m1.options.entries("WAVE") == m0.options.entries("WAVE")
    assert m1.title == m0.title
    # a second AFWRITE from the converted model gives the same HOUSE, SITE and POINT decks
    d2 = tmp_path / "back"
    d2.mkdir()
    ui.execute(f"ACTM,1")
    ui.execute(f"MDL,ssi,{d2}")
    ui.execute("AOPT,0,0,0,1,1,1,0,0,0,0,0,0,0,0")
    ui.execute("AFWRITE")
    for mod in ("HOUSE", "SITE", "POINT"):
        a = _read(decks.deck_path(d1, "ssi", mod), mod)
        b = _read(decks.deck_path(d2, "ssi", mod), mod)
        assert a.params == b.params, mod
        assert {k: t.rows for k, t in a.tables.items()} == {k: t.rows for k, t in b.tables.items()}, mod
    # the .pre written by CONVERT reads back to the converted model
    ui3 = _ui(tmp_path)
    ui3.run_file(tmp_path / "back.pre")
    assert not ui3.sink.texts(Kind.ERROR)
    assert ui3.model.same_state(converted)


def test_convert_ssi_without_site_and_legacy_files(tmp_path):
    d1 = tmp_path / "orig"
    d1.mkdir()
    ui = _ui(tmp_path)
    ui.run_text(SSI_MODEL.format(dir=d1))
    (d1 / "ssi.sit").unlink()
    (d1 / "ssi.poi").unlink()
    assert ui.execute(f"CONVERT,SSI,3,{d1 / 'ssi.hou'}")
    m3 = ui.models[3]
    rec = m3.options.record("SITE")
    assert rec.number(12) == 0.005 and rec.integer(13) == 4096 and rec.integer(5) == 3
    assert m3.topl == [1, 2] and m3.freq_sets == {1: [2, 4, 8, 16]}
    assert _any(ui, Kind.INFO, "no .sit deck")
    # a legacy SASSI2000 fixed-format file is refused (D-ANS-08)
    (tmp_path / "old.hou").write_text("  TEST HOUSE\n    1    2    3\n")
    assert not ui.execute(f"CONVERT,SSI,4,{tmp_path / 'old.hou'}")
    assert _any(ui, Kind.ERROR, "legacy SASSI2000")


# ======================================================================================
# ANSYS (export command)
# ======================================================================================
def test_ansys_command_file_names(tmp_path):
    ui = _ui(tmp_path)
    ui.execute(f"CONVERT,ANSYS,,{DATA / 'beam188_cantilever.cdb'},9.81")
    assert ui.execute("ANSYS")                                             # no MDL: unnamed.inp in the CWD
    assert (tmp_path / "unnamed.inp").exists() and _any(ui, Kind.WARNING, "D-MDL-03")
    mdir = tmp_path / "model"
    mdir.mkdir()
    ui.execute(f"MDL,cant,{mdir}")
    assert ui.execute("ANSYS") and (mdir / "cant.inp").exists()           # <model>.inp in the model directory
    assert ui.execute("ANSYS,other,sub") and (mdir / "sub" / "other.inp").exists()   # Dir relative to the CWD (MDL)
    assert ui.execute("ANSYS,exact.apdl,,1")
    text = (mdir / "exact.apdl").read_text()
    from sassi.model.values import fmt_num
    b = 0.02
    assert f"MP,EX,1,{fmt_num(2e11 * (1 - 2 * b * b))}" in text
    assert not ui.execute("ANSYS,x,,3")
    ui.execute("EDUOPT,ANSYSMODERN,1")
    assert ui.execute("ANSYS,modern")
    assert "ET,1,BEAM188" in (mdir / "modern.inp").read_text()


def test_ansys_export_runs_on_the_samples(tmp_path):
    """Every converted sample exports (legacy and modern) and the export reads back."""
    from sassi.io.ansys_cdb import convert_cdb, read_cdb, read_cdb_text
    from sassi.io.apdl import ApdlOptions, export_apdl
    for p in sorted(DATA.glob("*.cdb")):
        m = convert_cdb(read_cdb(p), 9.81, 0.05).model
        for g, e in m.iter_elements():
            if g.type in (1, 4):
                e.etype = 1                  # structure (the 2-D sample lies below the ground elevation)
        for modern in (False, True):
            res = export_apdl(m, ApdlOptions(modern=modern))
            back = convert_cdb(read_cdb_text(res.text), 9.81, 0.05).model
            assert back.n_elements() >= 1 or m.n_elements() == 0, p.name


@pytest.mark.parametrize("name", sorted(p.name for p in DATA.glob("*.cdb")))
def test_converted_samples_write_and_read_back_identically(name):
    """A converted model is plain SASSI-EDU data: WRITE -> INP gives the same state (UT-03)."""
    from sassi.io.ansys_cdb import convert_cdb, read_cdb
    from sassi.prep.writer import write_pre
    m = convert_cdb(read_cdb(DATA / name), 9.81, 0.05).model
    text, notes = write_pre(m)
    ui = Interpreter()
    ui.run_text(text)
    assert not ui.sink.texts(Kind.ERROR) and not notes
    assert ui.model.same_state(m)


# ======================================================================================
# ANSYSREFORMAT
# ======================================================================================
REFORMAT = """N,1,0,0,0
N,2,1,0,0
N,3,2,0,0
N,4,3,0,0
N,5,0,1,0
M,1,1e7,0.3,0.1,0.05,0.05,1
R,1,1,0,0,0.1,0.1,0.1
GROUP,1,BEAMS
E,1,1,2,5
E,2,2,3,5
E,3,3,4,5
KI,2,2,1,0,0,0,0,1,1
GROUP,2,BEAMS
E,1,1,3,5
EOUT,1,1,1,1,1,1,1,1,1,1,1,1,1,1-3
"""


def test_ansysreformat(tmp_path):
    ui = _ui(tmp_path)
    ui.run_text(REFORMAT)
    ui.execute("ACTM,1")
    assert ui.execute("ANSYSREFORMAT,0,beams.map")
    m = ui.model
    assert sorted(m.groups) == [1, 2, 3]
    assert [e.nodes for e in m.groups[1].sorted_elements()] == [[1, 2, 5], [3, 4, 5]]
    assert [e.nodes for e in m.groups[3].sorted_elements()] == [[2, 3, 5]]
    assert m.groups[3].elements[1].ki == [0, 0, 0, 0, 1, 1]
    assert [(r.group, r.elements) for r in m.eout] == [(1, [1, 2]), (3, [1])]
    rows = [ln.split() for ln in (tmp_path / "beams.map").read_text().splitlines() if not ln.startswith("!")]
    assert rows == [["1", "1", "1", "1"], ["1", "2", "3", "1"], ["1", "3", "1", "2"], ["2", "1", "2", "1"]]
    assert len(ui.models[0].groups) == 2                                   # original untouched
    # errors
    assert not ui.execute("ANSYSREFORMAT,0,again.map")                     # active model not empty
    ui.execute("ACTM,5")
    assert not ui.execute("ANSYSREFORMAT,7,x.map")                         # no such model
    assert not ui.execute("ANSYSREFORMAT,5,x.map")                         # Org = active


# ======================================================================================
# ANSYSMODELTYPE, GENMATRIXDAMP
# ======================================================================================
def test_ansysmodeltype_and_genmatrixdamp(tmp_path):
    from sassi.prep.writer import write_pre
    ui = _ui(tmp_path)
    assert ui.execute("ANSYSMODELTYPE,2")
    assert ui.model.options.record("ANSYSMODELTYPE").integer(1) == 2
    assert _any(ui, Kind.WARNING, "tier P2")
    assert "ANSYSMODELTYPE,2" in write_pre(ui.model)[0]
    assert not ui.execute("ANSYSMODELTYPE,3")
    ui.run_text("MXR,1,1,100\n")
    before = ui.model.matrices[1].to_json()
    assert ui.execute("GENMATRIXDAMP,1,,,0.05")
    assert ui.model.matrices[1].to_json() == before
    assert _any(ui, Kind.WARNING, "GENMATRIXDAMP is not usable in this version")
