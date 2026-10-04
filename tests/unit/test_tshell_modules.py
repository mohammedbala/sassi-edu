"""TSHELL through HOUSE, STRESS, CHECK and the thick-shell commands THSHLSTR / THSHLSMH
(requirements 1.5, 4.1, 4.10, 3.4.Q; D-ELM-08, D-TSH-01; spec 11 section 4)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sassi.elements import ElemRecord, assemble, build_dofmap, material_from_M
from sassi.io import decks
from sassi.io.files import read_container
from sassi.modules.stress import TSHELL_FACE_FILE, TSHELL_MAX_FILE, thshlstr_flag
from sassi.prep import Interpreter, Kind
from sassi.prep.commands import thickshell as TH
from sassi.verify.problems.vp_house import G_SI, run_house
from sassi.verify.problems.vp_stress import add_eout, field_file8, harmonic_file, read_tfu, run_stress, stress_deck
from sassi.verify.problems.vp_tshell import (DF54, DT54, FNUM54, MAT54, NFFT54, T54, read_face_file,
                                              tshell_plate_deck)

pytestmark = [pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning"),
              pytest.mark.filterwarnings("ignore:STRESS deck. unknown parameter 'thshlstr' ignored:UserWarning")]


def _msgs(ui, kind=None):
    return [m.text for m in ui.sink.messages if kind is None or m.kind == kind]


# ======================================================================================
# HOUSE
# ======================================================================================
def test_house_assembles_tshell_and_writes_recovery_data(tmp_path):
    d, _ = tshell_plate_deck(eint=1)
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out[-2000:]
    assert "EDU-06" not in out                                  # automatic drilling stiffness
    f4 = read_container(tmp_path / "m.N4", "FILE4")
    assert f4.meta["components"]["TSHELL"] == ["NXX", "NYY", "NXY", "QXZ", "QYZ", "MXX", "MYY", "MXY"]
    assert f4["rec_TSHELL_S"].shape == (4, 8, 24)
    assert np.allclose(f4["x_rec_TSHELL_thick"], T54) and np.allclose(f4["x_rec_TSHELL_E"], MAT54[2])
    assert np.allclose(f4["x_rec_TSHELL_nu"], MAT54[3]) and list(f4["x_rec_TSHELL_eint"]) == [1] * 4
    # Ks = element-library assembly of the same records
    mat = material_from_M(*MAT54[1:], G_SI)
    nodes = {int(r[0]): np.array(r[1:4], float) for r in d.table("nodes").rows}
    recs = [ElemRecord(int(r[0]), int(r[1]), 5, tuple(int(v) for v in r[7:11]), False, mat,
                       dict(thick=float(r[6]), eint=int(r[5]))) for r in d.table("elements").rows]
    dm = build_dofmap(sorted(nodes), None, recs)
    am = assemble(nodes, recs, dm)
    Ks = read_container(tmp_path / "COOSK", "COOSK").sparse("Ks")
    assert abs(Ks - am.Ks).max() <= 1e-12 * abs(am.Ks).max()
    assert list(f4["eq_node"]) == list(dm.eq_node)


@pytest.mark.parametrize("field,value,msg", [("eint", 2, "TSHELL EINT 2 must be 0 (reduced) or 1 (selective)"),
                                             ("thick", 0.0, "shell thickness 0.0 must be > 0")])
def test_house_rejects_bad_tshell_data(tmp_path, field, value, msg):
    d, _ = tshell_plate_deck()
    col = d.table("elements").columns.index(field)
    d.table("elements").rows[0][col] = value
    rc, out = run_house(tmp_path, d)
    assert rc == 1 and msg in out


def test_buried_tshell_is_structure(tmp_path):
    d, _ = tshell_plate_deck()
    col = d.table("elements").columns.index("etype")
    for r in d.table("elements").rows:
        r[col] = 2
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out[-1500:]
    row = [ln.split() for ln in out.splitlines() if ln.strip().startswith("TSHELL")][0]
    assert row[1:4] == ["4", "0", "4"]                       # structure, excavated, buried
    f4 = read_container(tmp_path / "m.N4", "FILE4")
    assert list(f4["elem_excav"]) == [0] * 4 and list(f4["x_elem_etype"]) == [2] * 4


# ======================================================================================
# STRESS
# ======================================================================================
EPS = np.array([2.0e-4, -0.5e-4, 1.0e-4])


@pytest.fixture()
def plate(tmp_path):
    """HOUSE run of the 2 x 2 TSHELL plate + a FILE8 of a constant membrane strain field."""
    d, (o, e1, e2) = tshell_plate_deck()
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out[-1500:]
    harmonic_file(tmp_path / "eq.acc", NFFT54, DT54, 7, 0.1)

    def field(node, x, f):
        q = np.asarray(x) - o
        X, Y = q @ e1, q @ e2
        u = EPS[0] * X + 0.5 * EPS[2] * Y
        v = EPS[1] * Y + 0.5 * EPS[2] * X
        return np.concatenate([u * e1 + v * e2, np.zeros(3)]).astype(complex)
    field_file8(tmp_path, FNUM54, DF54, field, NFFT54, DT54)
    return tmp_path


def _deck(**kw):
    d = stress_deck(NFFT54, DT54, **kw)
    add_eout(d, 1, [1, 2], [2] * 8)
    return d


def _run_with_param(wd, d, value):
    d["thshlstr"] = value                                  # STRESS deck parameter (lead schema, wave 3)
    p = decks.write(Path(wd) / "m.str", d)
    from sassi.modules.base import run_module
    rc = run_module("STRESS", "m", wd)
    return rc, (Path(wd) / "m_STRESS.out").read_text()


def test_stress_tshell_basic_components_without_faces(plate):
    rc, out = run_stress(plate, _deck(itran=1))
    assert rc == 0, out[-2000:]
    assert "(THSHLSTR): 0 (default" in out
    assert (plate / TSHELL_MAX_FILE).exists() and not (plate / TSHELL_FACE_FILE).exists()
    for comp in ("NXX", "QYZ", "MXY"):
        assert (plate / f"TSHELL_001_00001_{comp}.THS").exists()
    _, h = read_tfu(plate, "TSHELL", 1, 1, "NXX")
    E, nu = MAT54[2], MAT54[3]
    assert np.allclose(h, E * T54 / (1 - nu * nu) * (EPS[0] + nu * EPS[1]), rtol=1e-9)
    comps, rows = TH.read_tshell_max(plate / TSHELL_MAX_FILE)
    assert comps == ["NXX", "NYY", "NXY", "QXZ", "QYZ", "MXX", "MYY", "MXY"] and sorted(rows) == [(1, 1), (1, 2)]


def test_stress_faces_from_deck_parameter(plate):
    TH.write_thshlstr_file(plate, 0)                     # the deck line wins over THSHLSTR.opt
    rc, out = _run_with_param(plate, _deck(), 1)
    assert rc == 0, out[-2000:]
    assert "(THSHLSTR): 1 (STRESS deck" in out and "TSHELL face stresses and strains" in out
    faces = read_face_file(plate / TSHELL_FACE_FILE)
    assert sorted({k[2] for k in faces}) == ["++", "+-", "-+", "--"] and len(faces) == 8
    _, rows = TH.read_tshell_max(plate / TSHELL_MAX_FILE)
    n = rows[(1, 1)]
    assert faces[(1, 1, "++")][0] == pytest.approx(n[0] / T54, rel=1e-9)      # pure membrane: N/t
    assert faces[(1, 1, "+-")][0] == pytest.approx(n[0] / T54, rel=1e-9)
    assert faces[(1, 1, "--")][0] == pytest.approx(-n[0] / T54, rel=1e-9)
    (plate / TSHELL_FACE_FILE).unlink()
    TH.write_thshlstr_file(plate, 1)
    rc, out = _run_with_param(plate, _deck(), 0)
    assert rc == 0 and "(THSHLSTR): 0 (STRESS deck" in out and not (plate / TSHELL_FACE_FILE).exists()


def test_stress_faces_from_command_file_and_frame_warning(plate):
    TH.write_thshlstr_file(plate, 1)
    rc, out = run_stress(plate, _deck(savemax=1))
    assert rc == 0, out[-2000:]
    assert "(THSHLSTR): 1 (THSHLSTR.opt" in out and (plate / TSHELL_FACE_FILE).exists()
    assert "no text frame files are generated for the TSHELL elements" in out


def test_error_80_text_without_tier_note(plate):
    d = stress_deck(NFFT54, DT54)
    add_eout(d, 9, [1], [1])
    rc, out = run_stress(plate, d)
    assert rc == 1 and "Error 80" in out and "TSHELL is tier P2" not in out


def test_thshlstr_flag_sources(tmp_path):
    p = decks.write(tmp_path / "m.str", stress_deck(NFFT54, DT54))
    assert thshlstr_flag(p, tmp_path) == (0, "default")
    TH.write_thshlstr_file(tmp_path, 1)
    assert thshlstr_flag(p, tmp_path) == (1, "THSHLSTR.opt")
    p.write_text(p.read_text().replace("thshlstr = -1", "thshlstr = 0", 1))   # as AFWRITE writes it
    assert thshlstr_flag(p, tmp_path) == (0, "STRESS deck")


# ======================================================================================
# commands
# ======================================================================================
def test_thshlstr_command_record_file_and_round_trip(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.execute("THSHLSTR,1")
    assert ui.model.options.record("THSHLSTR").values == ["1"]
    assert any("model path is not defined" in t for t in _msgs(ui, Kind.INFO))
    ui.execute(f"MDL,m,{tmp_path}")
    ui.execute("THSHLSTR,1")
    assert (tmp_path / "THSHLSTR.opt").read_text() == "THSHLSTR,1\n"
    assert not ui.execute("THSHLSTR,2")
    assert any("<flag> must be 0" in t for t in _msgs(ui, Kind.ERROR))
    ui.execute("THSHLSTR,")
    assert ui.model.options.record("THSHLSTR").values == ["0"] and "THSHLSTR,0" in (tmp_path / "THSHLSTR.opt").read_text()
    ui.execute("THSHLSTR,1")
    ui.execute(f"WRITE,{tmp_path / 'w.pre'}")
    assert "THSHLSTR,1" in (tmp_path / "w.pre").read_text()
    ui2 = Interpreter(cwd=tmp_path)
    ui2.execute(f"INP,{tmp_path / 'w.pre'}")
    assert ui2.model.options.record("THSHLSTR").values == ["1"]


def test_parzen_weights_and_rings():
    assert np.allclose(TH.parzen_weights(0), [1.0])
    assert np.allclose(TH.parzen_weights(1), [1.0, 0.25])
    assert np.allclose(TH.parzen_weights(3), [1.0, 0.71875, 0.25, 0.03125])
    nb = {1: [2], 2: [1, 3], 3: [2, 4], 4: [3]}
    assert TH.ring_distances(nb, 1, 2) == {1: 0, 2: 1, 3: 2}


def test_smooth_shear_chain_and_groups():
    vals = {(1, 1): np.array([0.0, 0.0]), (1, 2): np.array([3.0, 6.0]), (1, 3): np.array([0.0, 0.0]),
            (2, 1): np.array([9.0, 9.0])}
    nodes = {(1, 1): [1, 2, 6, 5], (1, 2): [2, 3, 7, 6], (1, 3): [3, 4, 8, 7], (2, 1): [3, 4, 9, 10]}
    avg = TH.smooth_shear(vals, nodes, 1, typ=1)
    assert np.allclose(avg[(1, 2)], [1.0, 2.0]) and np.allclose(avg[(1, 1)], [1.5, 3.0])
    assert np.allclose(avg[(2, 1)], [9.0, 9.0])                  # other group: not mixed
    par = TH.smooth_shear(vals, nodes, 1, typ=0)
    assert np.allclose(par[(1, 2)], [2.0, 4.0]) and np.allclose(par[(1, 1)], [0.6, 1.2])
    assert np.allclose(TH.smooth_shear(vals, nodes, 0, typ=0)[(1, 2)], [3.0, 6.0])


def test_thshlsmh_command(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    cmds = [f"MDL,m,{tmp_path}", "N,1,0,0,0", "N,2,1,0,0", "N,3,2,0,0", "N,4,0,1,0", "N,5,1,1,0", "N,6,2,1,0",
            "M,1,3e7,0.2,24,0.05,0.05,1", "GROUP,1,TSHELL", "E,1,1,2,5,4", "E,2,2,3,6,5", "THICK,1,2,1,0.3"]
    for c in cmds:
        assert ui.execute(c), c
    assert not ui.execute("THSHLSMH,1")
    assert any("TSHELL_ELEMENT_MAX.TXT not found" in t for t in _msgs(ui, Kind.ERROR))
    (tmp_path / TSHELL_MAX_FILE).write_text("# header\n# group element NXX NYY NXY QXZ QYZ MXX MYY MXY\n"
                                            "1 1 1 1 1 4.0 2.0 1 1 1\n1 2 1 1 1 0.0 6.0 1 1 1\n9 9 0 0 0 1 1 0 0 0\n")
    assert ui.execute("THSHLSMH,1,1,out")
    assert any("NON VALIDATED" in t for t in _msgs(ui, Kind.WARNING))
    assert any("not TSHELL elements of the active model" in t for t in _msgs(ui, Kind.WARNING))
    lines = (tmp_path / "out" / TH.SMOOTHED_FILE).read_text().split("\n")
    assert lines[0] == "2 4"
    row = [float(v) for v in lines[1].split()]
    assert row[:2] == [1, 1] and row[2:] == pytest.approx([2.0, 4.0])
    assert not ui.execute("THSHLSMH,-1")
    assert not ui.execute("THSHLSMH,1,2")


# ======================================================================================
# CHECK
# ======================================================================================
def test_check_tshell_eint_and_thickness(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    for c in (f"MDL,m,{tmp_path}", "N,1,0,0,0", "N,2,1,0,0", "N,3,1,1,0", "N,4,0,1,0", "M,1,3e7,0.2,24,0.05,0.05,1",
              "GROUP,1,TSHELL", "E,1,1,2,3,4", "THICK,1,1,1,0.3", "CHECK"):
        ui.execute(c)
    assert not [t for t in _msgs(ui) if "EDU-26" in t and "TSHELL" in t]
    ui.model.groups[1].elements[1].eint = 2
    ui.model.groups[1].elements[1].thick = 0.0
    ui.execute("CHECK")
    msgs = [t for t in _msgs(ui) if "EDU-26" in t]
    assert any("EINT must be 0 (reduced) or 1 (selective)" in t for t in msgs)
    assert any("thickness (THICK) must be > 0" in t for t in msgs)


@pytest.mark.parametrize("scenario", ["mdl_change", "inp_then_mdl"])
def test_check_writes_thshlstr_file_from_the_record(tmp_path, scenario):
    """THSHLSTR,1 given in another directory, or restored by INP before MDL: the THSHLSTR.opt that STRESS reads
    is missing in the model directory until CHECK (run by AFWRITE) rewrites it from the model record."""
    from sassi.prep.check import Checker
    ui = Interpreter(cwd=tmp_path)
    model = ["N,1,0,0,0", "N,2,1,0,0", "N,3,1,1,0", "N,4,0,1,0", "M,1,3e7,0.2,24,0.05,0.05,1", "GROUP,1,TSHELL",
             "E,1,1,2,3,4", "THICK,1,1,1,0.3"]
    assert ui.execute(f"MDL,m,{tmp_path / 'A'}") and ui.execute("THSHLSTR,1")
    if scenario == "mdl_change":
        assert ui.execute(f"MDL,m,{tmp_path / 'B'}")
    else:
        assert ui.execute(f"WRITE,{tmp_path / 'w.pre'}")
        ui = Interpreter(cwd=tmp_path)
        ui.execute(f"INP,{tmp_path / 'w.pre'}")
        assert ui.execute(f"MDL,m,{tmp_path / 'B'}")
    for c in model:
        assert ui.execute(c), c
    assert thshlstr_flag(None, tmp_path / "B") == (0, "default")
    Checker(ui.model, modules=["STRESS"]).run()
    assert thshlstr_flag(None, tmp_path / "B") == (1, "THSHLSTR.opt")


def test_check_rewrites_a_stale_thshlstr_file_and_warns_without_directory(tmp_path):
    from sassi.prep.check import Checker
    ui = Interpreter(cwd=tmp_path)
    for c in (f"MDL,m,{tmp_path}", "N,1,0,0,0", "N,2,1,0,0", "N,3,1,1,0", "M,1,3e7,0.2,24,0.05,0.05,1",
              "GROUP,1,TSHELL", "E,1,1,2,3", "THICK,1,1,1,0.3"):
        assert ui.execute(c), c
    TH.write_thshlstr_file(tmp_path, 1)                     # left by another model: this one has no THSHLSTR
    Checker(ui.model, modules=["STRESS"]).run()
    assert thshlstr_flag(None, tmp_path) == (0, "THSHLSTR.opt")
    Checker(ui.model, modules=["STRESS"]).run()             # consistent: nothing written
    assert TH.read_thshlstr_file(tmp_path) == 0
    ui.model.path = None
    ui.execute("THSHLSTR,1")
    rep = Checker(ui.model, modules=["STRESS"]).run()
    assert any(m.number == "EDU-28" and "THSHLSTR,1" in m.detail for m in rep.messages)


def test_stress_skips_tshell_requests_without_codes_and_retires_stale_files(plate):
    """An EOUT row with all 8 codes 0 requests nothing: no row (no NaN) in the TSHELL files; a later run
    without TSHELL output requests renames the files of the earlier run (.prev) so that THSHLSMH cannot
    smooth stale maxima."""
    TH.write_thshlstr_file(plate, 1)
    d = stress_deck(NFFT54, DT54)
    add_eout(d, 1, [1], [1] * 8)
    add_eout(d, 1, [2], [0] * 8)
    rc, out = run_stress(plate, d)
    assert rc == 0, out[-1500:]
    _, rows = TH.read_tshell_max(plate / TSHELL_MAX_FILE)
    assert list(rows) == [(1, 1)] and np.all(np.isfinite(rows[(1, 1)]))
    assert "nan" not in (plate / "TSHELL_FACE_STRESSES.TXT").read_text().lower()
    assert "with all output codes 0" in out
    TH.write_thshlstr_file(plate, 0)
    rc, out = run_stress(plate, d)
    assert rc == 0 and (plate / TSHELL_MAX_FILE).exists() and not (plate / "TSHELL_FACE_STRESSES.TXT").exists()
    assert (plate / "TSHELL_FACE_STRESSES.TXT.prev").exists()
    d = stress_deck(NFFT54, DT54)
    add_eout(d, 1, [2], [0] * 8)
    rc, out = run_stress(plate, d)
    assert rc == 0 and not (plate / TSHELL_MAX_FILE).exists() and (plate / (TSHELL_MAX_FILE + ".prev")).exists()
