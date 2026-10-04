"""Element-stress files and section histories: .ess reader (D-STR-12), READSTR, SECDATAOPT, CALCSECTHIST
(spec 10 sections 3.1, 3.4, 3.9, 3.10; D-SEC-02, D-SEC-05, D-SEC-06)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.prep import Interpreter, Kind
from sassi.prep import cuts_lib as cl
from sassi.prep.writer import write_pre
from sassi.verify.problems.vp_cuts import shell_wall_lines, two_hex_lines


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def run(ui, lines):
    ui.run_text("\n".join(lines) if isinstance(lines, list) else lines)
    return ui.model


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


# ======================================================================================
# .ess reader
# ======================================================================================
def test_read_ess_tolerant_layout(tmp_path):
    p = tmp_path / "a.ess"
    p.write_text("# comment\n2\nSOLID 1 1 2\nSHELL 3 1 1\n\nSOLID 1 1\n  1 1.0D+01 2 3 4 5 6\n"
                 "2 -1 -2 -3 -4 -5 -6\n* another comment\nSHELL 3 1\n 7 1.5 2.5\n")
    d = cl.read_ess(p)
    assert d.groups == {1: "SOLID", 3: "SHELL"}
    assert np.allclose(d.values[(1, 1)], [10, 2, 3, 4, 5, 6]) and np.allclose(d.values[(1, 2)], -np.arange(1, 7))
    assert np.allclose(d.values[(3, 7)], [1.5, 2.5, 0, 0, 0, 0])            # short row padded
    p.write_text("1\nSOLID 1 1 1\nSHELL 1 1\n1 0 0 0 0 0 0\n")
    with pytest.raises(cl.SectionError, match="group table says SOLID"):
        cl.read_ess(p)
    p.write_text("1\n1 0 0 0 0 0 0\n")
    with pytest.raises(cl.SectionError, match="before the first block header"):
        cl.read_ess(p)
    p.write_text("")
    with pytest.raises(cl.SectionError, match="empty"):
        cl.read_ess(p)


def test_read_ess_reads_the_stress_module_writer(tmp_path):
    """File-format contract with STRESS (D-STR-12): frames written by sassi.core.stress_lib read back exactly."""
    SL = pytest.importorskip("sassi.core.stress_lib")
    if not hasattr(SL, "write_element_center"):
        pytest.skip("STRESS writer not available")
    rng = np.random.default_rng(3)
    v1, v2 = rng.normal(size=(3, 6)) * 1e3, rng.normal(size=(2, 6))
    blocks = [SL.CenterBlock("SOLID", 2, 1, [1, 2, 5], v1), SL.CenterBlock("SHELL", 4, 1, [3, 4], v2)]
    SL.write_element_center(tmp_path / "ESTRESS_00001.ess", blocks, fmt="{:.10e}")
    d = cl.read_ess(tmp_path / "ESTRESS_00001.ess")
    assert d.groups == {2: "SOLID", 4: "SHELL"}
    for k, e in enumerate([1, 2, 5]):
        assert np.allclose(d.values[(2, e)], v1[k], rtol=1e-10)
    for k, e in enumerate([3, 4]):
        assert np.allclose(d.values[(4, e)], v2[k], rtol=1e-10)


# ======================================================================================
# READSTR
# ======================================================================================
def test_readstr_attaches_and_checks(ui, tmp_path):
    m = run(ui, two_hex_lines())
    h0 = m.model_hash()
    cl.write_ess(tmp_path / "s.ess", {1: "SOLID"}, {(1, 1): [1, 2, 3, 4, 5, 6], (1, 7): [0] * 6})
    assert ui.execute("READSTR,s.ess")
    assert any("1 records of elements that are not in the model ignored" in w for w in warnings(ui))
    assert any("have no stress record: 1/2" in w for w in warnings(ui))
    tab = cl.stress_table(m)
    assert list(tab) == [(1, 1)] and np.allclose(tab[(1, 1)], [1, 2, 3, 4, 5, 6])
    assert m.model_hash() == h0                                       # results are not model content
    assert "element_stress" not in write_pre(m)[0]
    assert cl.stress_table(m.copy()) .keys() == tab.keys()            # SAVE / CPMODEL keep them
    # Dir argument and the maxima warning (D-SEC-06)
    (tmp_path / "res").mkdir()
    cl.write_ess(tmp_path / "res" / "ELEMENT_CENTER_ABS_MAX_STRESSES.TXT", {1: "SOLID"},
                 {(1, 1): [1] * 6, (1, 2): [2] * 6})
    assert ui.execute("READSTR,ELEMENT_CENTER_ABS_MAX_STRESSES.TXT,res")
    assert any("not occur at the same time" in w for w in warnings(ui))
    assert cl.stress_info(m)["maxima"] is True and len(cl.stress_table(m)) == 2


def test_readstr_rejects_files_of_another_model(ui, tmp_path):
    m = run(ui, two_hex_lines())
    cl.write_ess(tmp_path / "s.ess", {1: "SOLID"}, {(1, 1): [1] * 6})
    assert ui.execute("READSTR,s.ess")
    cl.write_ess(tmp_path / "b.ess", {1: "SHELL"}, {(1, 1): [9] * 6})
    assert not ui.execute("READSTR,b.ess")
    assert "group 1 is SHELL in the file but SOLID in the model" in errors(ui)[-1]
    cl.write_ess(tmp_path / "c.ess", {2: "SOLID"}, {(2, 1): [9] * 6})
    assert not ui.execute("READSTR,c.ess")
    assert "group 2 (SOLID) is not in the model" in errors(ui)[-1]
    assert np.allclose(cl.stress_table(m)[(1, 1)], 1.0)               # unchanged
    assert not ui.execute("READSTR,missing.ess")
    assert not ui.execute("READSTR,s.ess,nodir")


# ======================================================================================
# SECDATAOPT
# ======================================================================================
def test_secdataopt_record(ui):
    m = ui.model
    assert ui.execute("SECDATAOPT,1")
    assert m.options.record("SECDATAOPT").integer(1) == 1
    assert not warnings(ui)                                            # no tier placeholder message
    assert "SECDATAOPT,1" in write_pre(m)[0]
    assert not ui.execute("SECDATAOPT,2")
    assert m.options.record("SECDATAOPT").integer(1) == 1
    assert ui.execute("SECDATAOPT,0")
    assert m.options.record("SECDATAOPT").integer(1) == 0


def test_calcsecthistdb_is_a_p2_message(ui):
    ui.execute("CALCSECTHISTDB,1,0,0,0,0,0,1,1,0,0,1,0.01,,,,out.csv")
    assert any("CALCSECTHISTDB is not available in this build (tier P2)" in w for w in warnings(ui))


# ======================================================================================
# CALCSECTHIST
# ======================================================================================
def _frames(tmp_path, folder="NSTRESS"):
    d = tmp_path / folder
    d.mkdir(exist_ok=True)
    states = [{(1, 1): [0, 100, 0, 0, 0, 0], (1, 2): [0, -100, 0, 0, 0, 0]},
              {(1, 1): [0, 100, 0, 0, 0, 0], (1, 2): [0, 100, 0, 0, 0, 0]},
              {(1, 1): [0, -300, 0, 0, 0, 0], (1, 2): [0, 50, 0, 0, 0, 0]}]
    names = []
    for k, s in enumerate(states, start=1):
        cl.write_ess(d / f"ESTRESS_{k:05d}.ess", {1: "SHELL"}, s)
        names.append(f"ESTRESS_{k:05d}.ess")
    # the list STRESS writes: a '#' header line, names relative to the list's folder
    (d / "ESTRESS.lst").write_text("# ESTRESS frames of model x: dt = 0.01 s\n" + "\n".join(names) + "\n")
    return d


def test_calcsecthist_csv_and_max_row(ui, tmp_path):
    run(ui, shell_wall_lines() + ["SLICE,2,0,0,1,0,0,1"])
    _frames(tmp_path)
    assert ui.execute("CALCSECTHIST,NSTRESS/ESTRESS.lst,2,0,0,1,0,0,1,1,0,0,4,0.01,out/wall.csv")
    rows = [ln.split(",") for ln in (tmp_path / "out" / "wall.csv").read_text().splitlines()]
    assert rows[0] == ["Time", "Fx", "Fy", "Fz", "Mx", "My", "Mz"]
    assert [r[0] for r in rows[1:]] == ["0.000000E+00", "1.000000E-02", "2.000000E-02", "MAX"]
    fz = [float(r[3]) for r in rows[1:]]
    my = [float(r[5]) for r in rows[1:]]
    assert fz == pytest.approx([0.0, 200.0, -250.0, -250.0])          # signed abs max
    assert my == pytest.approx([200.0, 0.0, -350.0, -350.0])
    assert 4 in ui.model.csys                                          # sysno stored
    assert ui.session["calcsecthist"]["frames"] == 3


def test_calcsecthist_steps_ties_and_lists(ui, tmp_path):
    run(ui, shell_wall_lines() + ["SLICE,2,0,0,1,0,0,1"])
    d = _frames(tmp_path)
    (d / "plain.lst").write_text("ESTRESS_00002.ess\nESTRESS_00001.ess\n")    # no header line
    assert ui.execute(f"CALCSECTHIST,{d / 'plain.lst'},2,0,0,1,0,0,1,1,0,0,,,steps.csv")
    rows = [ln.split(",") for ln in (tmp_path / "steps.csv").read_text().splitlines()]
    assert rows[0][0] == "Step" and [r[0] for r in rows[1:]] == ["1", "2", "MAX"]
    assert float(rows[-1][5]) == pytest.approx(200.0)
    assert cl.signed_absmax(np.array([[-5.0, 1.0], [5.0, -1.0]])).tolist() == [-5.0, 1.0]   # ties: first
    (d / "bad.lst").write_text("header\nESTRESS_00009.ess\n")
    assert not ui.execute(f"CALCSECTHIST,{d / 'bad.lst'},2,0,0,1,0,0,1,1,0,0,0,0.01,x.csv")
    assert "not found" in errors(ui)[-1]
    assert not ui.execute("CALCSECTHIST,NSTRESS/ESTRESS.lst,9,0,0,1,0,0,1,1,0,0,0,0.01,x.csv")   # no cut 9
    assert not ui.execute("CALCSECTHIST,NSTRESS/ESTRESS.lst,2,0,0,1,0,0,1,0,0,1,0,0.01,x.csv")   # r // n


def test_calcsecthist_missing_records_warn(ui, tmp_path):
    run(ui, two_hex_lines() + ["CUTADD,1,1,1-2"])
    cl.write_ess(tmp_path / "f1.ess", {1: "SOLID"}, {(1, 1): [0, 0, 10, 0, 0, 0], (1, 2): [0, 0, 10, 0, 0, 0]})
    cl.write_ess(tmp_path / "f2.ess", {1: "SOLID"}, {(1, 1): [0, 0, 10, 0, 0, 0]})
    (tmp_path / "l.lst").write_text("frames\nf1.ess\nf2.ess\n")
    assert ui.execute("CALCSECTHIST,l.lst,1,0,0,0.5,0,0,1,1,0,0,0,1,h.csv")
    assert any("no record in some frames" in w and "1/2 (1 frames)" in w for w in warnings(ui))
    rows = [ln.split(",") for ln in (tmp_path / "h.csv").read_text().splitlines()]
    assert float(rows[1][3]) == pytest.approx(20.0) and float(rows[2][3]) == pytest.approx(10.0)


def test_calcsecthist_rejects_frames_of_another_model(ui, tmp_path):
    run(ui, two_hex_lines() + ["CUTADD,1,1,1-2"])
    cl.write_ess(tmp_path / "f1.ess", {1: "SHELL"}, {(1, 1): [1] * 6, (1, 2): [1] * 6})
    (tmp_path / "l.lst").write_text("frames\nf1.ess\n")
    assert not ui.execute("CALCSECTHIST,l.lst,1,0,0,0.5,0,0,1,1,0,0,0,1,h.csv")
    assert "group 1 is SHELL in the frame but SOLID in the model" in errors(ui)[-1]
    assert not (tmp_path / "h.csv").exists()
