"""The Run view (sassi/ui/runview.py, static/runview.js; requirements 7.20, D-W6-15).

While an INP file runs in the GUI session, the server reports the file's lines, what every line added to the
model (nodes, elements, soil layers, interaction nodes, fixities, masses, output requests) and live snapshots
of the model; the browser shows the listing with the section and the line being processed beside the model,
the part being processed highlighted (module steps: the part each step works on), and marks the same line in
the File Editor.  Checked on example 1's model part."""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from sassi.ui import runview as RV

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"
NODE = shutil.which("node")


@pytest.fixture(autouse=True)
def _settings(tmp_path, monkeypatch):
    monkeypatch.setenv("SASSI_EDU_SETTINGS_DIR", str(tmp_path / "settings"))


def _session(tmp_path):
    from sassi.ui.api import GuiSession
    from sassi.ui.learn import example_model_lines
    from sassi.ui.lessons import EXAMPLES_DIR, copy_example
    pre = copy_example("ex01_surface_stick", tmp_path / "ex", EXAMPLES_DIR)
    (pre.parent / "model.pre").write_text("\n".join(example_model_lines(pre.read_text())) + "\n", encoding="utf-8")
    return GuiSession(cwd=str(pre.parent), settings_dir=str(tmp_path / "settings"), web=True)


def _run(S, line="INP,model.pre"):
    seq0 = S.events.last
    st, r = S.route("POST", "/api/command", {}, {"lines": [line]})
    assert st == 200 and r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    evs, _, _ = S.events.since(seq0)
    return [e for e in evs if e["type"] == "run"]


def test_watching_example_1_being_built(tmp_path):
    S = _session(tmp_path)
    slept = []
    S.runview.sleep = slept.append                        # the pace, without waiting
    assert S.route("POST", "/api/runview", {}, {"enabled": True, "pace": 0.12}) == (200, {"enabled": True, "pace": 0.12})
    runs = _run(S)
    kinds = [e["event"] for e in runs]
    assert kinds[0] == "file" and kinds[-1] == "end" and "scene" in kinds
    f = runs[0]
    assert f["name"] == "model.pre" and f["depth"] == 1 and f["path"].endswith("model.pre")
    assert f["lines"] == Path(f["path"]).read_text().splitlines()                 # (MDL then moved into ex01/)
    lines = {e["cmd"]: e for e in runs if e["event"] == "line"}
    # every line that changed the model, with what it changed (watching: none is skipped)
    assert len(slept) == len(lines) and all(s == 0.12 for s in slept)
    assert lines["L,1,1.6,0.120,2000,1000,0.05,0.05"]["focus"]["layers"] == [1]
    assert lines["TOPL,1,1,1,1,1,1,1,1,1,1"]["focus"]["topl"] is True
    assert lines["FILL,1,9"]["focus"]["nodes"] == list(range(2, 9)) and lines["FILL,1,9"]["text"] == "nodes 2-8 added (7)"
    assert lines["NGEN,8,9,1,9,1,0,8,0"]["focus"]["nodes"] == list(range(10, 82))
    assert lines["EGEN,7,9,1,8"]["focus"]["elements"] == [[1, e] for e in range(9, 65)]
    assert lines["EGEN,7,9,1,8"]["text"] == "group 1 elements 9-64 (56)"
    spider = next(e for c, e in lines.items() if c.startswith("FOREACH,SPIDER"))
    assert spider["focus"]["elements"] == [[3, e] for e in range(1, 9)]
    assert lines["FIXROT"]["focus"]["fixed"][:3] == [1, 2, 3] and len(lines["FIXROT"]["focus"]["fixed"]) == 72
    assert lines["INT,1,81,1,1"]["focus"]["interaction"] == list(range(1, 82))
    masses = next(e for c, e in lines.items() if c.startswith("FOREACH,FLOOR"))
    assert masses["focus"]["mass"] == [82, 83, 84, 85]
    assert lines["NOUT,1,1,1,0,0,1,1,41,82-85"]["focus"]["requests"]["nout"] == [41, 82, 83, 84, 85]
    assert all(not c.startswith("*") for c in lines)       # comments change nothing
    # only the newest snapshot keeps its data in the event log
    scenes = [e for e in runs if e["event"] == "scene"]
    assert all(e["data"].get("stale") for e in scenes[:-1]) and not scenes[-1]["data"].get("stale")
    d = scenes[-1]["data"]
    assert d["kind"] == "MODELPLOT" and d["soil"] and d["scene"]["n_elements"] == 76
    assert len(d["run"]["interaction"]) == 81 and d["run"]["nout"] == [37, 41, 45, 82, 83, 84, 85]
    assert d["run"]["eout"] == [[2, 1], [2, 2], [2, 3], [2, 4]]
    json.dumps(d)


def test_full_speed_and_without_a_view(tmp_path):
    S = _session(tmp_path)
    S.runview.sleep = lambda s: pytest.fail("no pause at full speed")
    runs = _run(S)                                         # the view is not open: no snapshot
    assert [e["event"] for e in runs if e["event"] != "line"] == ["file", "end"]
    assert 1 <= sum(e["event"] == "line" for e in runs) < 45   # a fast burst: the last line of it is still shown


def test_snapshots_through_the_build(tmp_path):
    from sassi.prep import Interpreter
    from sassi.prep.commands.plotting import plot_state
    ui = Interpreter(cwd=tmp_path)
    st = plot_state(ui)
    assert RV.snapshot(ui, st)["empty"]
    for ln in ["MDL,t,t", "L,1,1.6,0.120,2000,1000,0.05,0.05", "L,2,3,0.13,6600,3300,0.02,0.02", "TOPL,1,1,1"]:
        assert ui.execute(ln)
    d = RV.snapshot(ui, st)
    assert d["kind"] == "LAYERPLOT" and len(d["table"]["layers"]) == 3        # soil only: its column
    for ln in ["N,1,0,0,0", "N,3,20,0,0", "FILL,1,3", "NGEN,2,3,1,3,1,0,10,0", "INT,1,9,1,1", "N,20,0,0,10"]:
        assert ui.execute(ln)
    d = RV.snapshot(ui, st)                                 # nodes, no element yet: all of them drawn
    sc = d["scene"]
    assert d["kind"] == "MODELPLOT" and sc["n_elements"] == 0 and len(sc["node_id"]) == 10 and len(sc["free"]) == 10
    assert len(sc["interaction"]) == 9 and sc["bbox"][:2] == [0.0, 20.0] and not d["soil"]
    for ln in ["GROUP,1,SHELL", "E,1,1,2,5,4", "THICK,1,1,1,1"]:
        assert ui.execute(ln)
    d = RV.snapshot(ui, st)                                 # elements: with the soil around them
    assert d["scene"]["n_elements"] == 1 and len(d["scene"]["free"]) == 6 and d["soil"]


def test_routes_settings_and_a_failing_hook(tmp_path):
    S = _session(tmp_path)
    assert S.route("GET", "/api/runview", {}, {})[1]["empty"]          # before any model
    assert S.route("POST", "/api/runview", {}, {"pace": 9})[1] == {"enabled": False, "pace": 2.0}
    assert S.route("POST", "/api/runview", {}, {"pace": "x"})[0] == 400
    S.runview.push = lambda *a, **k: (_ for _ in ()).throw(RuntimeError("display"))
    st, r = S.route("POST", "/api/command", {}, {"lines": ["INP,model.pre"]})
    assert r["ok"] and len(S.interp.model.nodes) == 90               # a failing display never stops the run
    st, d = S.route("GET", "/api/runview", {}, {})
    assert st == 200 and d["kind"] == "MODELPLOT" and d["name"] == "ex01"


def test_describe_and_ranges():
    assert RV.ranges([1, 2, 3, 7, 9, 10]) == "1-3, 7, 9-10"
    assert RV.ranges(range(1, 100, 2), limit=3).endswith("...")
    assert RV.describe({"model": 1, "nodes": [2, 3], "interaction": [5]}) == "nodes 2-3 added (2); interaction nodes 5 (1)"
    assert RV.describe({"model": 2, "all": True}) == "model 2"
    assert RV.describe({}) == ""


def test_event_log_replace():
    from sassi.ui.events import EventLog
    log = EventLog()
    s1 = log.push("run", event="scene", data={"big": 1})
    log.push("run", event="line")
    assert log.replace(s1, data={"stale": True}, other=1)
    ev = log.since(0)[0][0]
    assert ev["data"] == {"stale": True} and "other" not in ev
    assert not log.replace(999, data={})


@pytest.mark.skipif(NODE is None, reason="Node.js is not installed")
def test_sections_and_step_focus_in_the_browser_code():
    deck = (ROOT / "examples" / "ex01_surface_stick.pre").read_text(encoding="utf-8").splitlines()
    code = """
      globalThis.SASSI = {el: () => ({}), selectTab: () => {}};
      require(process.argv[1]);
      const R = SASSI.RunView, lines = JSON.parse(require("fs").readFileSync(0, "utf8"));
      const secs = R.sections(lines);
      const at = (i) => secs.find((s) => i >= s.start && i <= s.end) || null;
      const ngen = lines.findIndex((l) => l.startsWith("NGEN,8,9"));
      const fixrot = lines.indexOf("FIXROT");
      const r = {interaction: [1, 2], nout: [41], eout: [[2, 1]], loads: []};
      console.log(JSON.stringify({ngen: at(ngen), fixrot: at(fixrot), blank: at(lines.indexOf("*")), n: secs.length,
        impedance: R.stepFocus("ANALYS.impedance", {}, r), spectra: R.stepFocus("MOTION.spectra", {node: 85}, r),
        site: R.stepFocus("SITE.modes", {}, r), stress: R.stepFocus("STRESS.elements", {}, r), none: R.stepFocus("X.y", {}, r)}));"""
    p = subprocess.run([NODE, "-e", code, str(STATIC / "runview.js")], input=json.dumps(deck), capture_output=True,
                       text=True, timeout=60)
    assert p.returncode == 0, p.stderr[-2000:]
    res = json.loads(p.stdout)
    assert res["ngen"]["part"] == "structure" and res["ngen"]["title"].startswith("Mat nodes 1..81")
    assert res["fixrot"]["title"].startswith("A flat shell has no stiffness for the drilling rotation")
    assert res["blank"] is None and res["n"] > 15
    assert res["impedance"] == {"nodes": [1, 2]} and res["spectra"] == {"nodes": [85]}
    assert res["site"] == {"soil": True} and res["stress"] == {"elements": [[2, 1]]} and res["none"] == {}


def test_front_end_wiring():
    page = (STATIC / "index.html").read_text(encoding="utf-8")
    assert page.index('src="static/plots.js"') < page.index('src="static/runview.js"') < page.index('src="static/main.js"')
    app = (STATIC / "app.js").read_text(encoding="utf-8")
    assert 'SASSI.RunView.onEvent(ev); else SASSI.RunView.onProgress(ev);' in app and "SASSI.RunView.onJob(job)" in app
    assert 't.kind === "runview") ? "view" : "main"' in app                       # the right-hand group
    assert '{label: "Run View", action: () => SASSI.RunView && SASSI.RunView.open()' in app
    assert "S.editorHighlight = function (path, m)" in app and 'wrap: "off"' in app
    plots = (STATIC / "plots.js").read_text(encoding="utf-8")
    assert "P.drawModel = function (t, d, extra, opts)" in plots and "P.sceneHelpers = {faceList, meshTrace, outlineTrace, markerTrace};" in plots
    js = (STATIC / "runview.js").read_text(encoding="utf-8")
    assert 'S.post("/api/runview", {enabled: true, pace: R.pace()})' in js and 'S.get("/api/runview")' in js
    assert "SASSI.Learn.isLessonOpen()" in js                                    # not opened over a lesson
    css = (STATIC / "styles.css").read_text(encoding="utf-8")
    assert ".rv-ln.cur" in css and ".editor-hl .hl-cur" in css and "body.runview-shown .activity" in css
