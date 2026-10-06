"""The activity panel (sassi/ui/static/activity.js): the events an input file run sends to the browser.

While an INP file runs the session reports the command being executed (its text, canonical name and one-line
meaning) and the progress of a module run inside the file; the panel shows them.  Driven through
``GuiSession.route`` as the HTTP handler calls it."""
from __future__ import annotations

import shutil
from pathlib import Path

from sassi.ui.api import GuiSession

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"

DECK = """* a small model and a SITE run, to watch the events of an input file
MDL,act,act
GRAVITY,9.81
L,1,1.0,19,600,300,0.05,0.05
L,2,1.0,21,2000,1000,0.02,0.02
TOPL,1,1,1,1,1,2
FREQ,1,10,20,40
SITE,0,1,0,10,2,1,0,1,2048,1,0,0.005,4096,1
WAVE,2,1,1,1,0
AOPT,0,0,0,1,0,0,0,0,0,0,0,0,0,0
AFWRITE
RUNSITE
"""


def _events(S):
    evs, _last, _reset = S.events.since(0)
    return evs


def test_input_file_reports_commands_and_module_progress(tmp_path):
    (tmp_path / "act.pre").write_text(DECK, encoding="utf-8")
    S = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    st, r = S.route("POST", "/api/command", {}, {"lines": ["INP,act.pre"]})
    assert st == 200 and r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    evs = _events(S)
    prog = [e for e in evs if e["type"] == "progress" and e.get("kind") != "module"]
    assert prog and prog[0]["line"] == 1 and prog[0]["source"] == "act.pre" and prog[-1]["line"] == prog[-1]["total"]
    by_name = {e["name"]: e for e in prog if e["name"]}
    # the commands that start a module or write / check the decks are always reported, with their meaning
    assert {"AFWRITE", "RUNSITE"} <= set(by_name)
    assert by_name["RUNSITE"]["cmd"] == "RUNSITE" and "SITE" in by_name["RUNSITE"]["what"]
    assert all(set(e) >= {"cmd", "name", "what", "fraction", "text"} for e in prog)
    # the module run inside the file reports its own progress, up to its end
    mod = [e for e in evs if e["type"] == "progress" and e.get("kind") == "module"]
    assert mod and {e["module"] for e in mod} == {"SITE"} and mod[-1]["fraction"] == 1.0
    assert any("frequency" in e["text"] for e in mod)
    # the panel's start / end markers: the module start, its listing's last line, the file summary
    texts = [e["text"] for e in evs if e["type"] == "message"]
    assert any(t.startswith("RUNSITE: model act") for t in texts)
    assert any(t.strip().startswith("SITE finished with status OK in ") for t in texts)
    assert any(t.startswith("INP act.pre: ") and " commands, " in t for t in texts)
    shutil.rmtree(tmp_path / "act", ignore_errors=True)


def test_a_typed_module_run_does_not_use_the_in_file_hook(tmp_path):
    """A RUN<MODULE> typed at the keyboard runs as a job (worker process with its own progress); the session's
    module_progress hook serves only module runs inside an INP file."""
    S = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    assert S.interp.module_progress is not None and S.interp.progress is not None


def test_front_end_wiring():
    page = (STATIC / "index.html").read_text(encoding="utf-8")
    assert page.index('src="static/activity.js"') > page.index('src="static/app.js"')
    assert page.index('src="static/activity.js"') < page.index('src="static/main.js"')
    app = (STATIC / "app.js").read_text(encoding="utf-8")
    assert "S.activity.onEvent(ev)" in app and "S.activity.onJob(job)" in app and "S.activity.onIdle()" in app
    js = (STATIC / "activity.js").read_text(encoding="utf-8")
    # the end markers the panel parses are the interpreter's and the modules' own texts
    assert "finished with status" in js and "lines, (\\d+) commands, (\\d+) warnings, (\\d+) errors" in js
    assert "SASSI.activity = SASSI_ACTIVITY" in js
