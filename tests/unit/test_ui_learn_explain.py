"""Learn: the command explainer (GET /api/explain; sassi.ui.explain and sassi/ui/command_args.json) and the
front-end guards of the Learn features (static/learn.js, the Learn menu, the lesson dock, Command Entry "?").
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from sassi.prep.registry import load_commands, lookup
from sassi.ui import explain as X
from sassi.ui.api import GuiSession

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"


@pytest.fixture
def S(tmp_path):
    return GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))


def ex(S, line, cursor=None):
    q = {"line": [line]}
    if cursor is not None:
        q["cursor"] = [str(cursor)]
    st, r = S.route("GET", "/api/explain", q, {})
    assert st == 200, r
    return r


def arg(r, pos):
    return next(a for a in r["args"] if a["pos"] == pos)


# ====================================================================== the endpoint
def test_model_command_with_curated_arguments(S):
    r = ex(S, "NGEN,3,5,1,5,1,0,0,10")
    assert r["kind"] == "command" and r["name"] == "NGEN" and not r["abbreviated"]
    assert r["syntax"].startswith("NGEN,") and r["tier"] == "P0" and r["status"] == "yes"
    assert r["status_text"] == "implemented in this build" and r["category"] == "3.4.C"
    assert r["category_title"] == "Nodes and coordinate systems" and r["configures"] == "Model: nodes"
    assert r["args_from"] == "curated" and len(r["args"]) == 8
    a1, a8 = arg(r, 1), arg(r, 8)
    assert (a1["name"], a1["value"], a1["given"]) == ("itim", "3", True) and a1["default"] == "1" and a1["meaning"]
    assert (a8["name"], a8["value"]) == ("dz", "10")
    assert r["help"] == {"doc": "docs/reference/COMMAND_REFERENCE.md", "anchor": "34c-nodes-and-coordinate-systems"}
    assert r["dialog"] is None


def test_abbreviation_and_coded_values(S):
    r = ex(S, "GROU,1,2")
    assert r["name"] == "GROUP" and r["typed"] == "GROU" and r["abbreviated"] and "GROU" in r["abbrev"]
    t = arg(r, 2)
    assert t["name"] == "type" and t["value"] == "2" and t["value_meaning"].startswith("BEAMS")
    r = ex(S, "group,1,shell")                                       # names are case-insensitive
    assert r["name"] == "GROUP" and arg(r, 2)["value"] == "shell"
    r = ex(S, "D,1,4,1,1,UX,ROTZ")
    a5, a6 = arg(r, 5), arg(r, 6)
    assert a5["value_meaning"] == "translation along X" and a6["repeat"] and a6["value_meaning"] == "rotation about Z"
    assert arg(r, 4)["value_meaning"] == "fixed"


def test_module_option_record_arguments(S):
    """Module option commands use the typed record fields (names = dialog / WRITE names) plus the curated
    meanings; a blank argument shows its default; the dialog that edits the command is named."""
    from sassi.prep.options import OPTION_SPECS
    r = ex(S, "SITE,0,1,0,12,3,1,0,1,,1,0,0.01,1024,1")
    assert r["args_from"] == "record" and r["configures"]
    assert [a["name"] for a in r["args"]] == [f.name for f in OPTION_SPECS["SITE"].record.FIELDS]
    f2 = arg(r, 9)
    assert f2["name"] == "freq2" and f2["given"] is False and f2["value"] == ""
    assert arg(r, 13)["value"] == "1024" and "NFFT" in (arg(r, 13)["meaning"] + arg(r, 13)["label"])
    assert arg(r, 1)["label"] == "Operation Mode (0 solution, 1 data check)" and arg(r, 1)["value_meaning"]
    assert r["dialog"] == {"name": "ANALYSIS", "tab": "SITE", "label": "Options > Analysis > SITE tab"}
    # a command record of sassi.ui.cmdrecords (Option NON)
    r = ex(S, "EQL,0.8,,0,1,1")
    assert r["args_from"] == "record" and arg(r, 1)["name"] == "disp" and r["dialog"]["tab"] == "NONLINEAR"


def test_variable_lists_text_arguments_and_extra_values(S):
    r = ex(S, "TOPL,1,1,2,2,3")
    assert [a["value"] for a in r["args"]] == ["1", "1", "2", "2", "3"] and all(a["name"] for a in r["args"])
    r = ex(S, "TIT,Ex01 - stick, surface mat")                    # rest-of-line text (L7): one argument
    assert len(r["args"]) == 1 and r["args"][0]["value"] == "Ex01 - stick, surface mat"
    r = ex(S, "N,1,0,0,0,9,9,9,9")                                # 5-6: legacy fields, 7-8: not arguments
    extra = [a for a in r["args"] if a.get("extra")]
    assert [a["pos"] for a in extra] == [7, 8] and arg(r, 5)["value"] == "9" and not arg(r, 5).get("extra")
    r = ex(S, 'MDL,m,"/Users/me/Smith, J"')
    assert arg(r, 2)["value"] == "/Users/me/Smith, J"


def test_run_commands_and_comments(S):
    r = ex(S, "RUNS")
    assert r["kind"] == "unknown" and "RUNSITE" in r["suggestions"]       # no prefix matching (L10)
    r = ex(S, "RUNSITE")
    assert r["name"] == "RUNSITE" and r["category"] == "3.4.G" and "SITE" in r["configures"]
    r = ex(S, "* a comment, with commas")
    assert r["kind"] == "comment" and r["known"] and "not executed" in r["summary"]
    r = ex(S, "   ")
    assert r["kind"] == "blank"


def test_unknown_commands(S):
    r = ex(S, "NOD,1,0,0,0")
    assert r["kind"] == "unknown" and not r["known"] and r["typed"] == "NOD"
    assert "not a SASSI-EDU command" in r["message"] and isinstance(r["suggestions"], list)
    r = ex(S, "<script>alert(1)</script>")
    assert r["kind"] == "unknown"
    st, r = S.route("GET", "/api/explain", {"line": ["N,1"], "cursor": ["x"]}, {})
    assert st == 400


def test_cursor_argument_for_the_live_hint(S):
    line = "NGEN,3,5,1"
    assert ex(S, line, 2)["cursor_arg"] == 0                              # in the name
    assert ex(S, line, 5)["cursor_arg"] == 1                              # right after the first comma
    assert ex(S, line, 7)["cursor_arg"] == 2
    assert ex(S, line, len(line))["cursor_arg"] == 3
    assert ex(S, "N,", 2)["cursor_arg"] == 1
    assert ex(S, 'MDL,m,"a, b', 11)["cursor_arg"] == 2                   # a comma inside quotes
    assert ex(S, "TIT,a, b, c", 11)["cursor_arg"] == 1                   # the title is one argument
    assert ex(S, "XYZ,1", 4)["cursor_arg"] == 0


def test_status_of_placeholder_commands(S):
    infos = X.infos()
    stored = next((n for n, i in infos.items() if i.status == "stored only"), None)
    if stored:
        r = ex(S, stored)
        assert r["status"] == "stored only" and "stored" in r["status_text"]


# ====================================================================== the curated table
def test_command_args_table_schema_and_names():
    data = json.loads((ROOT / "sassi" / "ui" / "command_args.json").read_text(encoding="utf-8"))
    load_commands()
    entries = {k: v for k, v in data.items() if not k.startswith("_")}
    assert len(entries) >= 200
    for name, v in entries.items():
        spec = lookup(name)
        assert spec is not None and spec.name == name, name
        for key in ("summary", "configures", "source", "args"):
            assert v.get(key) is not None, (name, key)
        # plain text (the browser shows it with textContent): no HTML or Markdown markup
        assert v["summary"].strip() and not re.search(r"</?(a|b|i|em|strong|code|br|p|span|div|script)\b|\*\*|`", v["summary"]), name
        pos = [a["pos"] for a in v["args"]]
        assert pos == list(range(1, len(pos) + 1)), (name, pos)
        for a in v["args"]:
            assert a["name"] and isinstance(a.get("meaning", ""), str), (name, a)
            assert not a.get("repeat") or a is v["args"][-1], (name, "repeat only on the last argument")
            assert isinstance(a.get("values", {}), dict)
        fields = X.record_fields(name)
        if fields:                       # record commands: the same names in the same order
            for a in v["args"]:
                if a["pos"] <= len(fields):
                    assert a["name"] == fields[a["pos"] - 1].name, (name, a["pos"])


def _example_commands():
    load_commands()
    used = set()
    for f in sorted((ROOT / "examples").glob("*.pre")):
        for ln in f.read_text(encoding="utf-8").splitlines():
            s = ln.strip()
            if not s or s.startswith("*"):
                continue
            spec = lookup(re.split(r"[,\s]", s, maxsplit=1)[0])
            assert spec is not None, (f.name, s)
            used.add(spec.name)
    return sorted(used)


def test_every_command_of_the_examples_is_explained():
    used = _example_commands()
    assert len(used) > 80
    assert X.coverage(used) == []
    for name in used:                    # every argument of these commands has a meaning
        r = X.explain(name)
        assert r["summary"], name
        for a in r["args"]:
            assert (a["meaning"] or a["label"]), (name, a["pos"])


def test_commands_named_in_the_course_brief_are_explained():
    names = ("N NGEN FILL E EGEN GROUP M L R SC D INT MT MR MUNITS F MM FREQ DAMP TOPL WAVE SITE SOIL SPRO DYNP POINT "
             "HOUSE FORCE ANALYS MOTION NOUT STRESS EOUT RELD RELFILE RDND THFILE AOPT MDL TIT CD INP FCOPY FMOVE EQUAKE "
             "RSIN RSOUT ACCIN ACCOUT CMODFORM EDUOPT SITEX HOUSEX MOTIONX INTGEN ETYPE ETYPEGEN RUNEQUAKE RUNSOIL RUNSITE "
             "RUNPOINT RUNHOUSE RUNFORCE RUNANALYS RUNMOTION RUNSTRESS RUNRELDISP READSPEC SPECPLOT AXES READTH THPLOT "
             "MODELPLOT NODEPLOT LAYERPLOT SOILPROPPLOT CHECK AFWRITE ACTM DMODEL").split()
    assert X.coverage(list(names)) == []
    assert all(n in X.curated() for n in names)


# ====================================================================== front end
def _js(name):
    return (STATIC / name).read_text(encoding="utf-8")


def test_front_end_wiring():
    index, app, plots, learn, css = _js("index.html"), _js("app.js"), _js("plots.js"), _js("learn.js"), _js("styles.css")
    # learn.js loads after the modules it uses and before main.js; the dock and the "?" exist
    order = [index.index(f'"static/{f}"') for f in ("app.js", "dialogs.js", "plots.js", "learn.js", "main.js")]
    assert order == sorted(order)
    for needle in ('id="lessondock"', 'id="dock-resizer"', 'id="workarea"', 'id="cmd-help"', 'id="cmd-hint"'):
        assert needle in index
    # the Learn menu, Model > Open Example..., the toolbar button, start-up and clickable history lines
    assert '{title: "Learn", items: (SASSI.Learn && SASSI.Learn.menuItems()) || []}' in app
    assert 'label: "Open Example..."' in app and "SASSI.Learn.examplesDialog()" in app
    assert "S.rebuildMenus = buildMenubar" in app and "SASSI.Learn.init(" in app and "div.dataset.line = text" in app
    assert '["learn", "Learn: guided course' in plots and "learn: '<path" in plots
    for f in ("Start Page", "Open Example...", "Explain a Command...", "Course Workspace Folder...", "Reset Course Progress"):
        assert f'label: "{f}"' in learn, f
    # every contract action verb and every section key reaches the player
    from sassi.ui.lessons import ACTION_VERBS, SECTION_KEYS
    needs_run = re.search(r"const NEEDS_RUN = \[([^\]]*)\]", learn).group(1)
    assert all(f'"{v}"' in needs_run for v, files in ACTION_VERBS.items() if files)    # results need the step to have run
    for k in SECTION_KEYS.values():
        assert f".lpanel-{k}" in css, k
    # the runner goes through S.command (L17) and awaits module jobs; progress persists in try/catch storage
    assert "await S.command(line)" in learn and "L.waitJob(res.job.id)" in learn
    assert learn.count("window.localStorage") == 2 and "catch (e) { return dflt; }" in learn
    # Show answer, explain on lesson lines, Run step / Run steps 1-k, Reset, Previous / Next
    for needle in ("Show answer", 'L.explain(ln, row)', "▶ Run step", "Run steps 1–", "Reset lesson", "◀ Previous", "Next ▶"):
        assert needle in learn, needle
    # no native prompt / confirm / alert (they block the page); server HTML only through the renderer
    assert not re.search(r"(?<![\w.])(prompt|confirm|alert)\(", learn)
    assert learn.count(".innerHTML = markup") == 1


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_learn_js_syntax():
    for f in ("learn.js", "app.js", "plots.js"):
        r = subprocess.run(["node", "--check", str(STATIC / f)], capture_output=True, text=True, timeout=60)
        assert r.returncode == 0, (f, r.stderr)
