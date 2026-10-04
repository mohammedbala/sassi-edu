"""Animations of computed results in the guided course: the HARMFRAME command (steady-state harmonic frames
at one SSI frequency, SASSI-EDU extension) and the lesson action ``animate`` (PROCFRAME + DEFORMPLOT /
VECTORPLOT / BUBBLEPLOT / CONTOURPLOT as command text, rule L17).

HARMFRAME is verified on a hysteretic SDOF on a practically rigid site (the model of VP-01): the base moves
with the control motion (H = 1) and the mass has the closed-form transfer function H = k*/(k* - m w^2), so
every frame must equal Re(H e^{i wt}) with the known amplitude and phase, from the FILE8 and from the TFU
restart frames of MOTION alike."""
from __future__ import annotations

import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.plotting.state import (FrameStore, harmonic_frame_name, harmonic_frames, plot_state, read_frame,
                                  read_frame_list, write_frame)
from sassi.prep import Interpreter, Kind
from sassi.ui import learn
from sassi.ui import lessons as L
from sassi.verify import builders as B


@pytest.fixture(autouse=True)
def _settings(tmp_path, monkeypatch):
    """PROCFRAME writes SASSIani.xml in the settings folder: keep it in the test folder."""
    monkeypatch.setenv("SASSI_EDU_SETTINGS_DIR", str(tmp_path / "settings"))


def _run(ui, line):
    """Execute one line; returns (ok, messages of that line)."""
    n0 = len(list(ui.sink.messages))
    ok = ui.execute(line)
    return ok, [(m.kind, m.text) for m in list(ui.sink.messages)[n0:]]


def _errors(msgs):
    return [t for k, t in msgs if k == Kind.ERROR]


# ====================================================================== the harmonic maths
def test_harmonic_frames_amplitude_and_phase():
    H = np.array([[2.0 * np.exp(-1j * np.pi / 3), 0.5j, 0.0]])
    phases, U = harmonic_frames(H, 12)
    assert U.shape == (12, 1, 3) and np.allclose(phases, np.arange(12) * 30.0)
    phi = np.radians(phases)
    assert np.allclose(U[:, 0, 0], 2.0 * np.cos(phi - np.pi / 3))      # |H| cos(wt + theta), theta = -60 deg
    assert np.allclose(U[:, 0, 1], -0.5 * np.sin(phi))                 # 0.5 i: a quarter period ahead
    assert np.allclose(U[:, 0, 2], 0.0)
    assert int(np.argmax(U[:, 0, 0])) == 2                              # the peak at wt = 60 deg
    assert harmonic_frame_name(15.0, 2) == "HARM_015.0_00002"


# ====================================================================== HARMFRAME on an SDOF (VP-01 model)
BETA, FP, MASS = 0.05, 5.0, 1.0
FS = B.FrequencySet.harmonic(0.25, [4, 8, 12, 16, 18, 19, 20, 21, 22, 24, 28, 32, 40])     # 1 ... 10 Hz


def _h_sdof(f: float) -> complex:
    """Closed-form TF of the mass (SASSI complex modulus): k* / (k* - m w^2), k = m wp^2 / (1 - 2 beta^2)."""
    k = MASS * (2 * math.pi * FP) ** 2 / (1.0 - 2.0 * BETA ** 2)
    ks = k * complex(cfactor(BETA, 0))
    return ks / (ks - MASS * (2 * math.pi * f) ** 2)


@pytest.fixture(scope="module")
def sdof(tmp_path_factory):
    """FILE8 of a 1-t SDOF (beta = 5 %, SASSI-form peak at 5 Hz) on a practically rigid site, and the TFU
    restart frames of MOTION (Restart for TF) on the same FILE8."""
    from sassi.verify.problems.vp_analys import stiff_site
    wd = tmp_path_factory.mktemp("sdof")
    site = stiff_site()
    k = MASS * (2 * math.pi * FP) ** 2 / (1.0 - 2.0 * BETA ** 2)
    B.run_soil(wd, "m", site, FS, layer=0, rad=1.0)
    mdl = B.sdof_on_node(site, k, MASS, BETA)
    B.run_house(wd, "m", mdl)
    B.run_analys(wd, "m", FS)
    B.run_motion(wd, "m", FS, "", [(mdl["mass"], 1, 1, 0, 0, 0, 0, 0)], out=1, type=0, cm=0, ang=0.0, rsttf=1)
    return wd, mdl["base"], mdl["mass"]


def _frames(folder: Path):
    files = read_frame_list(folder)
    frs = [read_frame(f) for f in files]
    return files, frs


def _ui(wd: Path) -> Interpreter:
    ui = Interpreter(cwd=wd)
    plot_state(ui).auto_render = False
    return ui


def test_harmframe_sdof_known_amplitude_and_phase(sdof, tmp_path):
    wd, base, mass = sdof
    ui = _ui(wd)
    ok, msgs = _run(ui, "HARMFRAME,FILE8,4.9,HARM,24")
    assert ok, msgs
    info = " ".join(t for _, t in msgs)
    # 5.0 Hz (frequency number 20) is the computed frequency closest to 4.9 Hz, and the message says so
    assert "5 Hz (FILE8 frequency number 20)" in info and "closest to 4.9 Hz among 13" in info
    assert "largest amplitude 10.01 at node 2 X" in info and "total motion" in info
    files, frs = _frames(wd / "HARM")
    assert [f.name for f in files[:3]] == ["HARM_000.0_00001", "HARM_015.0_00002", "HARM_030.0_00003"]
    assert len(files) == 24
    H = _h_sdof(5.0)
    assert abs(abs(H) - 1.0 / (2 * BETA * math.sqrt(1 - BETA ** 2))) < 1e-9          # 10.0125 at the peak
    phi = np.radians(np.arange(24) * 15.0)
    for k, fr in enumerate(frs):
        assert list(fr.nodes) == [base, mass]
        row = {int(n): v for n, v in zip(fr.nodes, fr.values)}
        # the base follows the unit control motion cos(wt); the mass moves as |H| cos(wt + theta)
        assert np.allclose(row[base], [math.cos(phi[k]), 0.0, 0.0], atol=1e-5)
        assert row[mass][0] == pytest.approx(abs(H) * math.cos(phi[k] + np.angle(H)), abs=1e-4 * abs(H))
        assert row[mass][1] == 0.0 and row[mass][2] == 0.0
    # at resonance the mass lags the ground by about a quarter period (theta = -84.3 deg for beta = 5 %)
    assert math.degrees(np.angle(H)) == pytest.approx(-84.27, abs=0.01)
    xs = np.array([dict(zip(fr.nodes, fr.values[:, 0]))[mass] for fr in frs])
    assert int(np.argmax(xs)) == 6                          # wt = 90 deg, the frame nearest to -theta


def test_harmframe_references(sdof):
    wd, base, mass = sdof
    ui = _ui(wd)
    H = _h_sdof(5.0)
    phi = np.radians(np.arange(8) * 45.0)
    for ref, label in (("0", "relative to the free field (unit control motion in x')"),
                       (str(base), f"relative to node {base}")):
        ok, msgs = _run(ui, f"HARMFRAME,FILE8,5,REL{ref},8,{ref}")
        assert ok and label in " ".join(t for _, t in msgs), msgs
        _, frs = _frames(wd / f"REL{ref}")
        for k, fr in enumerate(frs):
            row = {int(n): v for n, v in zip(fr.nodes, fr.values)}
            assert np.allclose(row[base], 0.0, atol=1e-5)                       # the base moves with the ground
            assert row[mass][0] == pytest.approx(((H - 1) * np.exp(1j * phi[k])).real, abs=1e-4 * abs(H))


def test_harmframe_from_tfu_restart_frames_equals_file8(sdof):
    wd, base, mass = sdof
    assert (wd / "TFU").is_dir() and len(read_frame_list(wd / "TFU")) == len(FS.fnum)
    ui = _ui(wd)
    ok, msgs = _run(ui, "HARMFRAME,FILE8,5,A8,24")
    assert ok, msgs
    ok, msgs = _run(ui, "HARMFRAME,TFU,5,ATFU,24")
    assert ok, msgs
    assert "TFU frame TFU_005.00_00007" in " ".join(t for _, t in msgs)
    (_, a), (_, b) = _frames(wd / "A8"), _frames(wd / "ATFU")
    assert len(a) == len(b) == 24
    for fa, fb in zip(a, b):
        assert list(fa.nodes) == list(fb.nodes)
        assert np.allclose(fa.values, fb.values, rtol=1e-9, atol=1e-9)
    # the free-field reference needs the control direction of a FILE8
    ok, msgs = _run(ui, "HARMFRAME,TFU,5,ATFU,24,0")
    assert not ok and "needs a FILE8-type source" in _errors(msgs)[0]


def test_harmframe_errors_and_warnings(sdof, tmp_path):
    wd, base, mass = sdof
    ui = _ui(wd)
    for line, text in (("HARMFRAME,FILE8,5", "are required"), ("HARMFRAME,FILE8,0,H", "must be > 0 Hz"),
                       ("HARMFRAME,FILE8,5,H,3", "within 4..360"), ("HARMFRAME,FILE8,5,H,24,77", "reference node 77"),
                       ("HARMFRAME,NOFILE,5,H", "not found"), ("HARMFRAME,FILE8,5,H,24,-1", "<Ref> must be")):
        ok, msgs = _run(ui, line)
        assert not ok and text in " ".join(_errors(msgs)), (line, msgs)
    (wd / "notframes").mkdir(exist_ok=True)
    (wd / "notframes" / "x.txt").write_text("1 2\n1 0.5\n", encoding="utf-8")
    ok, msgs = _run(ui, "HARMFRAME,notframes,5,H")
    assert not ok and "no frame file named" in _errors(msgs)[0]
    ok, msgs = _run(ui, "HARMFRAME,FILE8,30,H30")
    assert ok and any(k == Kind.WARNING and "outside the computed range" in t for k, t in msgs)
    ok, msgs = _run(ui, "HARMFRAME,FILE8,5.6,H56")               # 5.5 Hz used: 1.8 % away, no warning
    assert ok and not any(k == Kind.WARNING for k, t in msgs)
    # a second run replaces the earlier HARM frames; other files in the folder are reported
    ok, _ = _run(ui, "HARMFRAME,FILE8,5,H56,6")
    assert ok and len(read_frame_list(wd / "H56")) == 6
    (wd / "H56" / "note.txt").write_text("x\n", encoding="utf-8")
    ok, msgs = _run(ui, "HARMFRAME,FILE8,5,H56,6")
    assert ok and any(k == Kind.WARNING and "holds other files" in t for k, t in msgs)


def test_harmframe_translates_optimizer_numbering(sdof, tmp_path):
    """HOUSE optimizer: FILE8 has the new node numbers; the frames are written in the model numbering
    (``<model>.map`` old -> new beside the FILE8), which DEFORMPLOT draws."""
    wd, base, mass = sdof
    work = tmp_path / "opt"
    work.mkdir()
    shutil.copy2(wd / "FILE8", work / "FILE8")
    (work / "m.map").write_text(f"# old new\n101 {base}\n102 {mass}\n", encoding="utf-8")
    ui = _ui(tmp_path)
    assert _run(ui, f"MDL,m,{work}")[0]
    ok, msgs = _run(ui, "HARMFRAME,FILE8,5,HARM,8")
    assert ok and "2 node number(s) translated to the model numbering (m.map)" in " ".join(t for _, t in msgs)
    _, frs = _frames(work / "HARM")
    assert list(frs[0].nodes) == [101, 102]


# ====================================================================== the lesson action
def test_parse_animate_and_labels():
    s = L.parse_animate("ex01/HARM | deformed 0.25 2 | SSI at 3.49 Hz")
    assert (s.path, s.kind, s.scale, s.stride, s.title) == ("ex01/HARM", "deformed", 0.25, 2, "SSI at 3.49 Hz")
    s = L.parse_animate("THD")
    assert (s.kind, s.scale, s.stride, s.title) == ("deformed", None, 1, "")
    s = L.parse_animate("RS | bubble 2")
    assert (s.kind, s.col, s.scale) == ("bubble", 2, None)
    assert L.parse_animate("X | vector auto 3").scale is None
    s = L.parse_animate("X | deformed 0.5 front | t")
    assert (s.scale, s.stride, s.view, s.title) == (0.5, 1, "front", "t")
    assert L.parse_animate("X | deformed top 0.5 2").view == "top" and L.parse_animate("X").view == ""
    for bad in ("", "a, b", "a | spin", "a | deformed 0", "a | deformed x", "a | bubble 0", "a | deformed 1 0",
                "a | deformed 1 2 3", "a | deformed | t | u", "a | deformed front top", "a | deformed 1 2 3 side"):
        with pytest.raises(ValueError):
            L.parse_animate(bad)
    assert learn.action_label("animate", "ex01/HARM | deformed | SSI at 3.49 Hz") == "Animate: SSI at 3.49 Hz"
    assert learn.action_label("animate", "ex01/HARM | vector") == "Animate ex01/HARM (vector)"
    assert "animate" in L.ACTION_VERBS


def _frame_folder(ws: Path) -> Path:
    """A hand-made frame folder: nodes 1-4 of a 1 m square shell, X = cos(wt), Z = 0.5 sin(wt)."""
    d = ws / "res" / "HARM"
    d.mkdir(parents=True)
    H = np.array([[1.0, 0.0, -0.5j]] * 4)
    phases, U = harmonic_frames(H, 8)
    for k in range(8):
        write_frame(d / harmonic_frame_name(phases[k], k + 1), np.arange(1, 5), U[k])
    return d


SQUARE = ["N,1,0,0,0", "N,2,10,0,0", "N,3,10,10,0", "N,4,0,10,0", "M,1,3e7,0.2,24,0.05,0.05", "GROUP,1,SHELL",
          "E,1,1,2,3,4", "THICK,1,1,1,0.5"]


def test_animate_resolves_to_command_text(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    folder = _frame_folder(ws)
    ui = _ui(ws)
    for ln in SQUARE:
        assert ui.execute(ln), ln
    lines = learn.resolve_action("animate", "res/HARM | deformed | Square at 5 Hz", ws, ui, [], None)["lines"]
    # automatic scale: the largest displacement (1) drawn as 15 % of the model size (10 m)
    assert lines == ["PROCFRAME,res/HARM,res/HARM_ani,Square at 5 Hz,3", "DEFORMPLOT,res/HARM_ani,1,8,1,1.5",
                     "WINDOWSETTINGS,UNDEFORMED,1", "WINDOWSETTINGS,TITLE,Square at 5 Hz"]
    for ln in lines:
        ok, msgs = _run(ui, ln)
        assert ok and not [t for k, t in msgs if k in (Kind.ERROR, Kind.WARNING)], (ln, msgs)
    st = plot_state(ui)
    p = st.active_plot
    assert p.kind == "DEFORMPLOT" and p.title == "Square at 5 Hz" and p.view.show_undeformed
    d = st.plot_data(p, ui)
    i2 = list(d["scene"]["node_id"]).index(2)
    assert d["frame"]["label"] == "frame 1: ωt = 0°" and d["frame"]["xyz"][i2] == pytest.approx([11.5, 0.0, 0.0])
    # the store itself: no second PROCFRAME; an explicit scale and stride
    lines = learn.resolve_action("animate", "res/HARM_ani | deformed 2 2 front", ws, ui, [], None)["lines"]
    assert lines == ["DEFORMPLOT,res/HARM_ani,1,8,2,2", "WINDOWSETTINGS,UNDEFORMED,1", "CNGVIEW,-90,0,0"]
    assert all(ui.execute(ln) for ln in lines)
    v = plot_state(ui).active_plot.view
    assert (v.rx, v.ry, v.rz) == (-90.0, 0.0, 0.0)              # front view: X to the right, Z up
    # bubble and contour: the column and its colour range over the frames
    lines = learn.resolve_action("animate", "res/HARM | bubble 3", ws, ui, [], None)["lines"]
    assert lines[0] == "PROCFRAME,res/HARM,res/HARM_ani,res/HARM,0" and lines[1] == "BUBBLEPLOT,res/HARM_ani,1,8,1,-0.5,0.5,3"
    assert all(ui.execute(ln) for ln in lines)
    assert learn.resolve_action("animate", "res/HARM | vector", ws, ui, [], None)["lines"][1] == \
        "VECTORPLOT,res/HARM_ani,1,8,1,1.5"
    # paths of another model folder are written relative to it
    assert ui.execute(f"MDL,other,{ws / 'm2'}")
    lines = learn.resolve_action("animate", "res/HARM", ws, ui, [], None)["lines"]
    assert lines[0] == "PROCFRAME,../res/HARM,../res/HARM_ani,res/HARM,3"
    assert all(ui.execute(ln) for ln in lines)


@pytest.mark.parametrize("args,status", [("res/none", 404), ("../x", 400), ("/tmp", 400), ("res/HARM | spin", 400),
                                         ("res/HARM | bubble 9", 400), ("res/empty", 404)])
def test_animate_errors(tmp_path, args, status):
    ws = tmp_path / "ws"
    ws.mkdir()
    _frame_folder(ws)
    (ws / "res" / "empty").mkdir()
    ui = _ui(ws)
    with pytest.raises(learn.LearnError) as e:
        learn.resolve_action("animate", args, ws, ui, [], None)
    assert e.value.status == status


def test_animate_through_the_gui_api(tmp_path):
    """POST /api/lessons/<id>/action, the command text through /api/command, then the frames of the plot."""
    from sassi.ui.api import GuiSession
    lessons = tmp_path / "lessons"
    lessons.mkdir()
    (lessons / "01_anim.md").write_text(ANIM_LESSON, encoding="utf-8")
    (tmp_path / "start").mkdir()
    S = GuiSession(cwd=str(tmp_path / "start"), settings_dir=str(tmp_path / "settings"))
    S.lesson_dir = lessons
    st, r = S.route("POST", "/api/lessons/97-animate/open", {}, {"confirm": True})
    assert st == 200, r
    ws = Path(r["workspace"])
    for ln in r["commands"] + SQUARE:
        assert S.route("POST", "/api/command", {}, {"lines": [ln]})[1]["ok"], ln
    _frame_folder(ws)
    st, a = S.route("POST", "/api/lessons/97-animate/action", {}, {"verb": "animate", "args": "res/HARM | deformed 2"})
    assert st == 200 and a["kind"] == "commands" and a["lines"][0].startswith("PROCFRAME,res/HARM,res/HARM_ani"), a
    st, res = S.route("POST", "/api/command", {}, {"lines": a["lines"]})
    assert res["ok"], res
    pid = S.plots.active
    st, d = S.route("GET", f"/api/plot/{pid}", {"frame": ["3"]}, {})
    assert st == 200 and d["label"] == "frame 3: ωt = 90°"
    hist = S.route("GET", "/api/state", {}, {})[1]["history"]
    assert any(h.startswith("DEFORMPLOT,res/HARM_ani") for h in hist)          # L17: replayable
    st, r = S.route("GET", "/api/lessons/97-animate", {}, {})
    assert r["steps"][0]["actions"][0]["label"] == "Animate: Mat at 0.39 Hz" and r["steps"][0]["actions"][0]["known"]


# ====================================================================== lesson validation and headless run
ANIM_LESSON = """---
id: 97-animate
title: Animation sample
part: Fundamentals
order: 97
minutes: 2
summary: a surface mat solved at two frequencies and animated
objectives: [Animate]
prerequisites: []
---
Intro.

## Solve and write harmonic frames

```sassi
MDL,run,run
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
AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
HARMFRAME,FILE8,0.39,HARM,8
```

```action
animate: run/HARM | deformed | Mat at 0.39 Hz
animate: run/HARM | bubble 1 2
```
"""


def test_headless_run_checks_animate_actions(tmp_path):
    les = L.parse_lesson(ANIM_LESSON, "anim.md")
    assert L.validate_lesson(les) == []
    rep = L.run_lesson_headless(les, tmp_path / "root")
    assert rep.ok, rep
    ws = Path(rep.workspace)
    assert FrameStore(ws / "run" / "HARM_ani").nframes == 8
    assert (ws / "SASSIani.xml").is_file()                     # the animation database stays in the workspace
    assert not list((ws / "run").glob("*.png"))                # no image rendered
    # a folder the step does not write is reported as missing; a malformed action fails validation
    les.steps[0].blocks[-1].text = "animate: run/NOPE\nanimate: run/HARM | spin"
    assert any("plot kind 'spin'" in p for p in L.validate_lesson(les))
    rep = L.run_lesson_headless(les, tmp_path / "root")
    assert ("step 1 (Solve and write harmonic frames)", "run/NOPE") in rep.missing and rep.errors
