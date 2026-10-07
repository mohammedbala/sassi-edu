"""The run summary window (sassi/ui/runsummary.py, static/summary.js; requirements 7.20, D-W6-12).

After a run the GUI shows the key inputs (from the decks the modules read), the key outputs (from the result files
they wrote) and their graphs.  Checked on a real run of example 1 (a 4-mass stick on a surface mat, ft, kip, s):
the SSI resonance of its transfer functions is 3.47 Hz at the top of the stick (lesson 1, examples/README.md)."""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from sassi.plotting.state import jsonable
from sassi.ui import runsummary as RS

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"


@pytest.fixture(scope="module")
def ex01(tmp_path_factory):
    """Example 1 run in the browser version's session (module runs inline: run_pending)."""
    import os
    from sassi.ui.api import GuiSession
    from sassi.ui.lessons import EXAMPLES_DIR, copy_example
    tmp = tmp_path_factory.mktemp("rs")
    before = os.environ.get("SASSI_EDU_SETTINGS_DIR")
    os.environ["SASSI_EDU_SETTINGS_DIR"] = str(tmp / "settings")      # PROCFRAME / settings files: not the user's
    pre = copy_example("ex01_surface_stick", tmp / "ex", EXAMPLES_DIR)
    S = GuiSession(cwd=str(pre.parent), settings_dir=str(tmp / "settings"), web=True)
    st, r = S.route("POST", "/api/command", {}, {"lines": [f"INP,{pre.name}"]})
    assert st == 200 and r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    while S.jobs.run_pending():
        pass
    yield S
    if before is None:
        os.environ.pop("SASSI_EDU_SETTINGS_DIR", None)
    else:
        os.environ["SASSI_EDU_SETTINGS_DIR"] = before
    shutil.rmtree(tmp, ignore_errors=True)


def _rows(sections, title_part):
    sec = next(s for s in sections if title_part in s["title"])
    return {k: v for k, v in sec["rows"]}


def test_summary_of_example_1(ex01):
    s = RS.summary(ex01.interp, "ex01")
    assert s["model"]["name"] == "ex01" and s["units"] == "ft" and not s["note"]
    assert [m["module"] for m in s["modules"]] == ["SITE", "POINT", "HOUSE", "ANALYS", "MOTION", "STRESS", "RELDISP"]
    assert all(m["ok"] and m["errors"] == 0 for m in s["modules"])
    # key inputs, from the decks
    model = _rows(s["inputs"], "Model")
    assert model["Nodes / elements"].startswith("90 / 76") and "81 (243 DOFs); surface foundation" == model["Interaction nodes"]
    assert model["Structure mass X Y Z"].startswith("368.7")          # 11,872 kips / 32.2
    soil = _rows(s["inputs"], "Soil profile")
    assert soil["Layers"] == "22 to depth 55 ft" and soil["Half-space"].startswith("Vs 3300 ft/s")
    motion = _rows(s["inputs"], "input motion")
    assert motion["Input motion"].startswith("rg160h_030g.acc") and motion["Peak, duration"].startswith("0.3239 g, 20 s")
    assert motion["Spectrum damping"] == "2 %, 5 %"
    # key outputs, from the result files: the SSI resonance at the top of the stick
    tf = _rows(s["outputs"], "Transfer functions")
    assert tf["Largest amplitude"].startswith("13.11 at 3.467 Hz (node 85 X)")
    acc = _rows(s["outputs"], "Peak accelerations")
    assert acc["Input motion (peak)"] == "0.3239 g" and acc["node 85 X"].startswith("1.278 g")
    rs = next(x for x in s["outputs"] if x["title"].startswith("In-structure response spectra"))
    assert "(5 % damping)" in rs["title"] and rs["rows"][0][0] == "node 85 X"            # the design damping
    rel = _rows(s["outputs"], "Relative displacements")
    assert "node 41 X" not in rel                                                          # the reference node: 0
    # graphs: transfer functions, acceleration, spectra with the input motion's, the soil profile
    ids = [c["id"] for c in s["charts"]]
    assert ids == ["tf", "acc", "rs", "soil"]
    ch = {c["id"]: c for c in s["charts"]}
    assert ch["tf"]["files"][0] == "00085TR_X.TFI" and len(ch["tf"]["lines"]) <= RS.MAX_LINES
    assert [ln["name"] for ln in ch["acc"]["lines"]] == ["node 85 X (largest)", "input motion (rg160h_030g.acc)"]
    assert max(len(ln["x"]) for ln in ch["acc"]["lines"]) <= RS.MAX_POINTS + 2
    assert ch["rs"]["logx"] and ch["rs"]["lines"][-1]["name"] == "input motion" and ch["rs"]["files"][0].endswith("02.RS")
    assert ch["soil"]["reversey"] and not ch["soil"]["files"]
    jsonable(s)                                                                            # the API answer


def test_through_the_api_and_open_as_plot(ex01):
    S = ex01
    st, s = S.route("GET", "/api/run_summary", {"model": ["ex01"]}, {})
    assert st == 200 and s["model"]["name"] == "ex01" and [c["id"] for c in s["charts"]][0] == "tf"
    st, s0 = S.route("GET", "/api/run_summary", {}, {})                                   # the active model
    assert st == 200 and s0["model"]["number"] == S.interp.active_model
    st, p = S.route("POST", "/api/run_summary/plot", {}, {"model": "ex01", "chart": "rs"})
    assert st == 200 and p["lines"][0].startswith("READSPEC,") and p["lines"][-1] == "AXES,,,,,1"
    assert any(ln.startswith("SPECPLOT,") for ln in p["lines"])
    st, r = S.route("POST", "/api/command", {}, {"lines": p["lines"]})
    assert r["ok"] and S.plots.plots[S.plots.active].kind == "SPECPLOT"
    st, p = S.route("POST", "/api/run_summary/plot", {}, {"model": "ex01", "chart": "acc"})
    assert st == 200 and any(ln.startswith("THPLOT,") for ln in p["lines"])
    assert S.route("POST", "/api/run_summary/plot", {}, {"model": "ex01", "chart": "soil"})[0] == 400
    assert S.route("GET", "/api/run_summary", {"model": ["nope"]}, {})[0] == 404


def test_a_model_without_a_run(tmp_path):
    from sassi.prep import Interpreter
    ui = Interpreter(cwd=tmp_path)
    s = RS.summary(ui)
    assert s["modules"] == [] and "nothing has run" in s["note"]
    ui.execute("MDL,fresh,fresh")
    s = RS.summary(ui, "fresh")
    assert s["modules"] == [] and "no module has run" in s["note"] and s["charts"] == []


def test_front_end_wiring():
    page = (STATIC / "index.html").read_text(encoding="utf-8")
    assert page.index('src="static/activity.js"') < page.index('src="static/summary.js"') < page.index('src="static/main.js"')
    act = (STATIC / "activity.js").read_text(encoding="utf-8")
    assert "S.onRunFinished({kind: ended.kind" in act and "run.models.push(mm[1])" in act
    js = (STATIC / "summary.js").read_text(encoding="utf-8")
    assert 'S.get("/api/run_summary"' in js and 'S.post("/api/run_summary/plot"' in js
    # not while a lesson is open, not after an intermediate module, not when switched off
    assert "SASSI.Learn.isLessonOpen()" in js and "if (!R.auto() || S.keepFocus" in js
    assert 'const RESULT_MODULES = ["MOTION", "STRESS", "RELDISP", "SOIL", "EQUAKE", "COMBIN", "NONLINEAR"];' in js
    learn = (STATIC / "learn.js").read_text(encoding="utf-8")
    assert "L.isLessonOpen = () => !!(P.lesson && dock() && !dock().hidden);" in learn
    app = (STATIC / "app.js").read_text(encoding="utf-8")
    assert '{label: "Run Summary", action: () => S.openRunSummary && S.openRunSummary()' in app
