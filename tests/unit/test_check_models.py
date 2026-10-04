"""Model-level checks and the model-checking / conditioning commands (requirements 3.4.H, 3.4.I; spec 09
sections 2.7-2.10 and 3; D-HOU-01, D-CHK-05/06/09) and the MOTION scaling rules (UT-07)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.io import thfile
from sassi.prep import Interpreter, Kind
from sassi.prep.check import MODEL, ModelView, run_check

CUBE = """
N,1,0,0,{z0}
N,2,1,0,{z0}
N,3,1,1,{z0}
N,4,0,1,{z0}
N,5,0,0,{z1}
N,6,1,0,{z1}
N,7,1,1,{z1}
N,8,0,1,{z1}
GROUP,1,SOLID
E,1,1,2,3,4,5,6,7,8
"""


def _ui(text):
    ui = Interpreter()
    ui.run_text(text)
    return ui


def test_etype_zero_resolution_d_hou_01():
    ui = _ui(CUBE.format(z0=-1, z1=0))                      # touching grade from below: excavated
    v = ModelView(ui.model)
    assert v.elems[0].etype == 2
    ui = _ui(CUBE.format(z0=-0.5, z1=0.5))                  # crossing grade: structure
    assert ModelView(ui.model).elems[0].etype == 1
    ui = _ui(CUBE.format(z0=-1, z1=0) + "GROUNDELEV,-1")    # above the ground elevation
    assert ModelView(ui.model).elems[0].etype == 1
    ui = _ui(CUBE.format(z0=-1, z1=0) + "ETYPE,1,1,1,1")    # explicit ETYPE wins
    assert ModelView(ui.model).elems[0].etype == 1


def test_inverted_solid_gives_edu05_and_distorted_face_w2():
    ui = _ui(CUBE.format(z0=0, z1=-1).replace("E,1,1,2,3,4,5,6,7,8", "E,1,1,2,3,4,5,6,7,8"))
    rep = run_check(ui.model, modules=["HOUSE"])
    assert rep.has("Error", "EDU-05", MODEL)               # top face below the bottom face: det J < 0
    ui = _ui(CUBE.format(z0=0, z1=1) + "N,9,0.05,0.5,0\nGROUP,2,SHELL\nE,1,1,9,4\nM,1,1e5,0.2,0.1,0.05,0.05")
    rep = run_check(ui.model, modules=["HOUSE"])
    assert any(m.number == 2 and "Group 2" in m.text for m in rep.warnings(MODEL))


def test_interaction_node_interfaces_edu01_and_bottom_up_edu21():
    base = CUBE.format(z0=-5, z1=0) + "L,1,5,0.12,1500,800,0.05,0.05\nL,2,10,0.13,2400,1200,0.02,0.02\nTOPL,1\n" \
                                      "SITE,0,1,0,20,2\nPOINT,0,1,1\nINT,1,8,1,1\n"
    rep = run_check(_ui(base).model, modules=["HOUSE"])
    assert not rep.has("Error", "EDU-01") and "EDU-21" not in rep.numbers()                   # bottom-up: fine
    inv = base.replace("N,1,0,0,-5", "N,1,0,0,0").replace("N,2,1,0,-5", "N,2,1,0,0") \
        .replace("N,3,1,1,-5", "N,3,1,1,0").replace("N,4,0,1,-5", "N,4,0,1,0") \
        .replace("N,5,0,0,0", "N,5,0,0,-5").replace("N,6,1,0,0", "N,6,1,0,-5") \
        .replace("N,7,1,1,0", "N,7,1,1,-5").replace("N,8,0,1,0", "N,8,0,1,-5") \
        .replace("E,1,1,2,3,4,5,6,7,8", "E,1,5,6,7,8,1,2,3,4")                               # top nodes first
    rep = run_check(_ui(inv).model, modules=["HOUSE"])
    assert rep.errors(MODEL) == [] and any(m.number == "EDU-21" for m in rep.warnings(MODEL))  # coherent: warning
    rep = run_check(_ui(inv + "HOUSE,32.2,0,0,2,0,1").model, modules=["HOUSE"])
    assert rep.has("Error", "EDU-21", MODEL)                                                  # incoherent: error
    rep = run_check(_ui(base + "N,9,0,0,-2\nINT,9,9,1,1").model, modules=["HOUSE"])
    assert any(m.number == "EDU-01" and "Node 9" in m.text for m in rep.errors(MODEL))
    rep = run_check(_ui(base + "POINT,0,0,1").model, modules=["HOUSE"])                       # below POINT zone
    assert rep.has("Error", "EDU-01", MODEL)


def test_edu22_excavation_boundary_nodes_and_edu08():
    base = CUBE.format(z0=-5, z1=0) + "L,1,5,0.12,1500,800,0.05,0.05\nL,2,10,0.13,2400,1200,0.02,0.02\n" \
                                      "TOPL,1\nSITE,0,1,0,20,2\nPOINT,0,1,1\n"
    rep = run_check(_ui(base + "INT,1,4,1,1").model, modules=["HOUSE"])
    nodes = sorted(int(m.text.split()[3]) for m in rep.warnings(MODEL) if m.number == "EDU-22")
    assert nodes == [5, 6, 7, 8]
    rep = run_check(_ui(base + "INT,1,8,1,1\nMSET,1,1,1,2").model, modules=["HOUSE"])
    assert any(m.number == "EDU-08" for m in rep.warnings(MODEL))


def test_check_and_err_file(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(f"MDL,chk,{tmp_path}\n" + CUBE.format(z0=0, z1=1) + "M,1,1e5,0.2,0.1,0.05,0.05\n"
                "AOPT,0,0,0,0,0,1,0,0,0,0,0,0,0,0\nCHECK")
    err = (tmp_path / "chk.err").read_text()
    assert err.startswith("CHECK: Errors and Warning for - chk") and "Errors and Warnings for HOUSE" in err
    assert ui.session["check_report"].modules == ["HOUSE"]


def test_used_freespring_kint_intcount_commands():
    ui = _ui(CUBE.format(z0=0, z1=1) + "N,20,5,5,5\nN,21,6,5,5\nINT,21,21,1,1\nN,22,5,6,5\nGROUP,2,SPRING\n"
                                       "E,1,20,22\nSC,1,1,1,1,0,0,0,0\nUSED\nFREESPRING\nINTCOUNT")
    m = ui.model
    assert m.nodes[21].fix == [0] * 6                                       # unused interaction node: reported
    assert any("unused interaction nodes not fixed" in w for w in ui.sink.texts(Kind.WARNING))
    infos = ui.sink.texts(Kind.INFO)
    assert any("FREESPRING: 2 free spring nodes" in t for t in infos)
    assert any(t.startswith("INTCOUNT: 1 interaction nodes") for t in infos)


def test_fixsldrot_fixsprrot():
    ui = _ui(CUBE.format(z0=0, z1=1) + "N,20,0,0,3\nGROUP,2,SPRING\nE,1,8,20\nSC,1,10,0,10,0,5,0,0\n"
                                       "MT,20,0,1,0\nFIXSLDROT\nFIXSPRROT")
    m = ui.model
    assert m.nodes[1].fix == [0, 0, 0, 1, 1, 1]
    assert m.nodes[8].fix == [0, 0, 0, 1, 0, 1]                  # solid + spring: ROTY has spring stiffness
    assert m.nodes[20].fix == [0, 0, 0, 1, 0, 1]                 # spring-only: UY has mass, ROTY stiffness


def test_fixshlrot_and_fixrot_oblique():
    c, s = math.cos(math.radians(30)), math.sin(math.radians(30))
    text = "".join(f"N,{k + 1},{x},{y * c},{y * s}\n" for k, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]))
    text += "GROUP,1,SHELL\nE,1,1,2,3,4\nM,1,1e5,0.2,0.1,0.05,0.05\n"
    ui = _ui(text + "FIXSHLROT,5")
    m = ui.model
    g = m.groups[2]
    assert g.type == 7 and len(g.elements) == 4 and len(m.nodes) == 8
    sc = m.springs[1]
    assert (sc.scx, sc.scxx) == (0.0, 0.0) and sc.scyy == pytest.approx(5 * s * s) and sc.sczz == pytest.approx(5 * c * c)
    assert m.nodes[5].fix == [1] * 6 and np.allclose(m.node_global(5), m.node_global(1))
    ui.execute("FIXSHLROT")                                       # second call: nodes already restrained
    assert len(m.groups[2].elements) == 4 and len(m.groups) == 2
    rep = run_check(m, modules=["HOUSE"])
    assert not any(x.number == "EDU-06" for x in rep.warnings(MODEL))


def test_ut07_motion_scaling_rules():
    ui = _ui("MOTION,0,0,0,0,0,0.1,100,301,0,0")
    from sassi.prep.options import get_record
    rec = get_record(ui.model, "MOTION")
    assert ("Error", 77, "") in rec.problems()
    ui.execute("MOTION,0,0,0,0,0,0.1,100,301,2,0.5")
    assert any(p[1] == 78 for p in get_record(ui.model, "MOTION").problems())
    a = np.array([0.1, -0.4, 0.25, 0.05])
    assert np.array_equal(thfile.scale_history(a, 0.0, 0.3), a * 0.3 / 0.4)      # a max / max|a|
    assert np.array_equal(thfile.scale_history(a, 2.0, 0.0), 2.0 * a)
    ui.execute("SOILX,0,1,0.1")
    rep = run_check(ui.model, modules=["SOIL"])
    assert rep.has("Error", 78, "SOIL")


def test_edu19_quiet_zone(tmp_path):
    (tmp_path / "a.acc").write_text("0.005\n" + "\n".join("0.01" for _ in range(200)) + "\n")
    text = ("L,1,5,0.12,1500,800,0.05,0.05\nTOPL,1,1\nFREQ,1,2,4\nTHFILE,{f}\nDAMP,0.05\n"
            "NOUT,1,1,1,0,0,1,1,1\nN,1,0,0,0\nSITE,0,1,0,20,1,1,0,1,2048,1,0,0.005,{n},1")
    f = tmp_path / "a.acc"
    rep = run_check(_ui(text.format(f=f, n=4096)).model, modules=["MOTION"])
    assert "EDU-19" not in rep.numbers()            # 19.5 s >= ln(100)/(2 pi 20 Hz 0.05) = 0.73 s
    rep = run_check(_ui(text.format(f=f, n=256)).model, modules=["MOTION"])
    assert rep.has("Warning", "EDU-19", "MOTION")   # 0.28 s quiet zone


def test_solid_face_distortion_and_prism():
    text = CUBE.format(z0=0, z1=1).replace("N,3,1,1,0", "N,3,0.5,0.02,0") + "M,1,1e5,0.2,0.1,0.05,0.05"
    rep = run_check(_ui(text).model, modules=["HOUSE"])
    assert any(m.number == 2 and m.text.endswith("Face 1") for m in rep.warnings(MODEL))
    prism = CUBE.format(z0=0, z1=1).replace("E,1,1,2,3,4,5,6,7,8", "E,1,1,2,3,3,5,6,7,7") + "M,1,1e5,0.2,0.1,0.05,0.05"
    rep = run_check(_ui(prism).model, modules=["HOUSE"])
    assert rep.errors(MODEL) == [] and not any(m.number == "EDU-05" for m in rep.messages)
