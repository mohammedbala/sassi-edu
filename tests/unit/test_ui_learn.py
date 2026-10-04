"""Learn: the guided course in the GUI (sassi.ui.learn and the /api/lessons, /api/examples routes).

The session is driven through ``GuiSession.route`` exactly as the HTTP handler calls it, and the lesson
commands are run the way the browser's lesson player runs them: one ``POST /api/command`` per line, a
``RUN<MODULE>`` line starts a worker job whose end is awaited before the next line (rule L17, UI-04).
The lessons are test fixtures (the real course files are validated by tests/unit/test_lessons.py)."""
from __future__ import annotations

import time
from pathlib import Path

import pytest

from sassi.ui import learn
from sassi.ui import lessons as L
from sassi.ui.api import GuiSession
from sassi.ui.settings import IniSettings

SAMPLE = """---
id: 99-sample
title: Sample lesson
part: Fundamentals
order: 1
minutes: 3
summary: a tiny lesson used by the GUI tests
objectives: [Define nodes, Make a group]
prerequisites: []
---
Intro text with a [theory link](docs/theory/THEORY_MANUAL.md#3-flexible-volume-substructuring) and an
[internal one](docs/internal/lesson_format.md).

```sassi-setup
TIT,sample
```

## Define two nodes
Narrative of step 1.

```sassi
* two nodes
N,1,0,0,0
N,2,1,0,0
```

### What this does
`N,<nd>,<x>,<y>,<z>` defines a node.

```sassi-show
N,<nd>,<x>,<y>,<z>
```

### Why it matters
Nodes carry the masses.

### Check yourself
How many nodes?
Answer: two.

```action
plot-nodes
explain: N,1,0,0,0
```

## Group
```sassi
GROUP,1,2
```

### Something else
```text
## not a heading
```
"""

EXAMPLE_LESSON = """---
id: 98-example
title: Lesson on an example
part: Design applications
order: 2
minutes: 5
example: ex03_forced_vibration
summary: uses example 3
objectives: [Load it]
prerequisites: [99-sample]
---
Intro.

## Read it
```sassi
TIT,read
```
"""

# a surface foundation that runs SITE in about a second (tests/unit/test_ui_jobs.py)
RUN_LESSON = """---
id: 97-run
title: A lesson with a module run
part: Advanced
order: 3
minutes: 2
summary: runs SITE
objectives: [Run a module]
prerequisites: []
---
Intro.

## Build
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
AOPT,0,0,0,1,1,0,0,0,0,0,0,0,0,0
AFWRITE
```

## Run SITE
```sassi
RUNSITE
* the free field is in FILE1
TIT,after the run
```

```action
open-listing: SITE
plot-layers
```
"""


@pytest.fixture
def S(tmp_path):
    start = tmp_path / "start"
    start.mkdir()
    lessons = tmp_path / "lessons"
    lessons.mkdir()
    (lessons / "01_sample.md").write_text(SAMPLE, encoding="utf-8")
    (lessons / "02_example.md").write_text(EXAMPLE_LESSON, encoding="utf-8")
    (lessons / "03_run.md").write_text(RUN_LESSON, encoding="utf-8")
    s = GuiSession(cwd=str(start), settings_dir=str(tmp_path / "settings"))
    s.lesson_dir = lessons
    yield s
    run = s.jobs.running()
    if run is not None:
        s.jobs.cancel(run.id)
        _wait(s, run.id, 20)


def call(S, method, path, body=None, query=None):
    q = {k: [str(v)] for k, v in (query or {}).items()}
    return S.route(method, path, q, body or {})


def _wait(S, jid, timeout=120.0):
    t0 = time.time()
    while S.jobs.get(jid).active:
        if time.time() - t0 > timeout:
            raise AssertionError(f"job {jid} still running after {timeout} s")
        time.sleep(0.05)
    return S.jobs.get(jid).to_dict()


def run_like_the_player(S, lines):
    """The lesson player's runner: one command per request; a RUN<MODULE> job is awaited."""
    out = []
    for ln in lines:
        st, r = call(S, "POST", "/api/command", {"lines": [ln]})
        assert st == 200, r
        if r.get("job"):
            d = _wait(S, r["job"]["id"])
            assert d["ok"], [m["text"] for m in d["messages"]]
            # the commands after the RUN line wait for it: none were deferred by the one-line request
            assert r["deferred"] == []
            out.append(("job", d))
        else:
            assert r["ok"], (ln, [m["text"] for m in r["messages"] if m["kind"] == "ERROR"])
            out.append(("ok", r))
        # a RUN line of a later request must not start while a job runs (the player waits)
        while S.jobs.running() is not None:
            time.sleep(0.05)
    return out


def open_lesson(S, lid, confirm=True):
    st, r = call(S, "POST", f"/api/lessons/{lid}/open", {"confirm": confirm})
    assert st == 200, r
    run_like_the_player(S, r["commands"])
    return r


# ====================================================================== outline and rendering
def test_outline_by_part_and_parse_errors(S, tmp_path):
    (S.lesson_dir / "04_broken.md").write_text("no front matter\n## step\n", encoding="utf-8")
    st, r = call(S, "GET", "/api/lessons")
    assert st == 200
    assert [p["part"] for p in r["parts"]] == ["Fundamentals", "Design applications", "Advanced"]
    les = r["parts"][0]["lessons"][0]
    assert les["id"] == "99-sample" and les["nsteps"] == 2 and les["minutes"] == 3
    assert les["objectives"] == ["Define nodes", "Make a group"] and les["step_titles"] == ["Define two nodes", "Group"]
    assert les["workspace_exists"] is False
    assert r["count"] == 3 and r["errors"] and r["errors"][0]["file"] == "04_broken.md"
    assert r["root"] == str((tmp_path / "start" / "sassi-course").resolve())


def test_lesson_detail_rendering(S):
    st, d = call(S, "GET", "/api/lessons/99-sample")
    assert st == 200
    assert d["setup"] == ["TIT,sample"]
    # the setup block is shown in the introduction; documents link to the Help tab, docs/internal does not
    assert '<div class="lesson-block" data-block="0"></div>' in d["intro_html"]
    assert d["intro_blocks"] == [{"k": 0, "kind": "sassi-setup", "lang": "sassi-setup", "lines": ["TIT,sample"]}]
    assert 'data-doc="docs/theory/THEORY_MANUAL.md"' in d["intro_html"]
    assert 'data-anchor="3-flexible-volume-substructuring"' in d["intro_html"]
    assert 'data-doc="docs/internal' not in d["intro_html"] and 'class="md-nolink"' in d["intro_html"]
    s1, s2 = d["steps"]
    assert s1["commands"] == ["* two nodes", "N,1,0,0,0", "N,2,1,0,0"]
    assert '<div class="lesson-block" data-block="0"></div>' in s1["narrative_html"]
    assert [(b["k"], b["kind"]) for b in s1["blocks"]] == [(0, "sassi"), (1, "sassi-show")]
    keys = [x["key"] for x in s1["sections"]]
    assert keys == ["what", "why", "check"]
    what, why, check = s1["sections"]
    assert what["title"] == "What this does" and 'data-block="1"' in what["html"]
    assert "How many nodes?" in check["html"] and "two." in check["answer_html"] and "two." not in check["html"]
    assert "plot-nodes" not in check["html"]                     # action blocks become buttons
    assert s1["actions"] == [{"verb": "plot-nodes", "args": "", "label": "Plot the nodes", "known": True},
                             {"verb": "explain", "args": "N,1,0,0,0", "label": "Explain N,1,0,0,0", "known": True}]
    assert s2["sections"][0]["key"] == "other" and "## not a heading" in s2["sections"][0]["html"]
    assert d["session"] is None and d["workspace"].endswith("99-sample")


def test_block_numbering_matches_the_parser():
    """The placeholders number the fenced blocks in document order exactly as lessons.Step.blocks does."""
    for text in (SAMPLE, RUN_LESSON, EXAMPLE_LESSON):
        les = L.parse_lesson(text, "x.md")
        for st in les.steps:
            k, kinds = 0, []
            for md in [st.narrative] + [x.markdown for x in st.sections]:
                _, bl = learn._fence_pass(md, k)
                k += len(bl)
                kinds += [b["kind"] for b in bl]
            assert kinds == [b.kind for b in st.blocks]


def test_unknown_and_invalid_lesson_ids(S):
    assert call(S, "GET", "/api/lessons/00-nope")[0] == 404
    assert call(S, "GET", "/api/lessons/../etc")[0] == 404          # the route only accepts [0-9a-z-] ids
    assert call(S, "POST", "/api/lessons/00-nope/open", {"confirm": True})[0] == 404


# ====================================================================== workspaces
def test_open_lesson_prepares_workspace_and_runs_setup(S, tmp_path):
    # an unsaved model outside the course: the learner is asked first, nothing is prepared
    st, r = call(S, "POST", "/api/command", {"lines": ["MDL,mine,mine", "N,1,0,0,0", "ACTM,3", "N,7,1,1,1"]})
    assert st == 200 and r["ok"]
    st, r = call(S, "POST", "/api/lessons/98-example/open", {})
    assert st == 200 and r["needs_confirm"] and {u["number"] for u in r["unsaved"]} == {0, 3}
    root = tmp_path / "start" / "sassi-course"
    assert not (root / "98-example").exists()
    r = open_lesson(S, "98-example")
    ws = Path(r["workspace"])
    assert ws == (root / "98-example").resolve()
    assert r["commands"][:4] == [f"CD,{ws}", "ACTM,0", "DMODEL,3", "DMODEL,0"] and r["setup"] == []
    # the example, its data and the marker are there; the session is in the lesson workspace with one empty model
    assert (ws / "ex03_forced_vibration.pre").is_file() and (ws / "data" / "ricker_5hz.th").is_file()
    assert (ws / learn.MARKER).is_file()
    assert Path(S.interp.cwd) == ws and list(S.interp.models) == [0] and S.interp.model.is_empty()
    st, state = call(S, "GET", "/api/state")
    assert state["lesson"] == {"id": "98-example", "workspace": str(ws), "ran": [], "opened": state["lesson"]["opened"]}
    st, d = call(S, "GET", "/api/lessons/98-example")
    assert d["session"]["workspace"] == str(ws) and d["workspace_exists"]
    # a model of a course workspace is rebuilt by its lesson: opening another lesson does not ask again
    assert call(S, "POST", "/api/command", {"lines": ["MDL,les,les", "N,1,0,0,0"]})[1]["ok"]
    st, r2 = call(S, "POST", "/api/lessons/99-sample/open", {})
    assert st == 200 and not r2.get("needs_confirm") and r2["setup"] == ["TIT,sample"]
    run_like_the_player(S, r2["commands"])
    assert S.interp.model.title == "sample"
    # the files of the workspace are readable through the file API (the course root is an allowed root)
    st, f = call(S, "GET", "/api/file", query={"name": str(ws / "ex03_forced_vibration.pre")})
    assert st == 200 and "MDL,ex03,ex03" in f["text"]


def test_reopen_recreates_only_course_folders(S, tmp_path):
    r = open_lesson(S, "99-sample")
    ws = Path(r["workspace"])
    (ws / "scratch.txt").write_text("x")
    open_lesson(S, "99-sample")                               # fresh: the learner's scratch file is gone
    assert ws.is_dir() and not (ws / "scratch.txt").exists() and (ws / learn.MARKER).is_file()
    # a folder of the same name that the course did not create is never deleted
    other = tmp_path / "start" / "sassi-course" / "98-example"
    other.mkdir(parents=True)
    (other / "precious.dat").write_text("keep")
    st, r = call(S, "POST", "/api/lessons/98-example/open", {"confirm": True})
    assert st == 409 and "not created by the course" in r["error"]
    assert (other / "precious.dat").read_text() == "keep"


def test_progress_is_recorded_for_the_open_lesson(S):
    open_lesson(S, "99-sample")
    st, r = call(S, "POST", "/api/lessons/99-sample/progress", {"step": 2})
    assert st == 200 and r["session"]["ran"] == [2]
    st, r = call(S, "POST", "/api/lessons/99-sample/progress", {"step": 2, "clear": True})
    assert r["session"]["ran"] == []
    assert call(S, "POST", "/api/lessons/98-example/progress", {"step": 1})[0] == 409
    assert call(S, "POST", "/api/lessons/99-sample/progress", {"step": "x"})[0] == 400


def test_step_run_streams_a_module_run(S, tmp_path):
    """Run step: the lines go one by one through /api/command; RUNSITE becomes a worker job whose listing is
    streamed into the Command History; the next line runs after the job; the actions then find the results."""
    r = open_lesson(S, "97-run")
    ws = Path(r["workspace"])
    st, d = call(S, "GET", "/api/lessons/97-run")
    s1, s2 = d["steps"]
    run_like_the_player(S, s1["commands"])
    before = S.events.last
    res = run_like_the_player(S, s2["commands"])
    kinds = [k for k, _ in res]
    assert kinds == ["job", "ok", "ok"]
    job = res[0][1]
    assert job["line"] == "RUNSITE" and job["state"] == "done"
    assert any("SITE finished with status OK" in m["text"] for m in job["messages"])
    evs = S.events.since(before)[0]
    streamed = [e for e in evs if e["type"] == "message" and str(e.get("source", "")).startswith("job ")]
    assert any("SITE finished" in e["text"] for e in streamed)          # the listing reached the history
    assert any(e["type"] == "job" and e["job"]["state"] == "done" for e in evs)
    assert S.interp.model.title == "after the run"                      # the line after RUNSITE ran after it
    assert "RUNSITE" in S.interp.history and S.interp.history.index("RUNSITE") < S.interp.history.index("TIT,after the run")
    assert (ws / "run" / "FILE1").is_file()
    # actions of the step
    st, a = call(S, "POST", "/api/lessons/97-run/action", {"verb": "open-listing", "args": "SITE"})
    assert st == 200 and a["kind"] == "file" and Path(a["path"]).is_file() and Path(a["path"]).name.lower() == "run_site.out"
    st, a = call(S, "POST", "/api/lessons/97-run/action", {"verb": "plot-layers", "args": ""})
    assert a == {"kind": "commands", "lines": ["LAYERPLOT"]}
    assert call(S, "POST", "/api/command", {"lines": a["lines"]})[1]["ok"]


# ====================================================================== actions
def _ws_with_results(S):
    r = open_lesson(S, "99-sample")
    ws = Path(r["workspace"])
    (ws / "res").mkdir()
    (ws / "res" / "a01.RS").write_text("0.5 0.10\n1.0 0.30\n5.0 0.80\n20.0 0.35\n", encoding="utf-8")
    (ws / "res" / "b01.RS").write_text("0.5 0.12\n1.0 0.25\n5.0 0.60\n20.0 0.30\n", encoding="utf-8")
    (ws / "res" / "a.ACC").write_text("0.01\n0.0\n0.1\n-0.05\n0.02\n0.0\n", encoding="utf-8")
    return ws


def test_plot_actions_build_command_text(S):
    ws = _ws_with_results(S)
    assert call(S, "POST", "/api/command", {"lines": [f"CD,{ws / 'res'}"]})[1]["ok"]
    st, a = call(S, "POST", "/api/lessons/99-sample/action", {"verb": "plot-spectrum", "args": "res/a01.RS, res/b01.RS | log"})
    assert st == 200 and a["kind"] == "commands"
    # relative to the working directory when that names the same file; consecutive free line numbers
    assert a["lines"] == ["READSPEC,a01.RS,1,1", "READSPEC,b01.RS,1,2", "SPECPLOT,1,2", "AXES,,,,,1"]
    st, r = call(S, "POST", "/api/command", {"lines": a["lines"]})
    assert r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    p = S.plots.active_plot
    assert p is not None and p.params["lines"] == [1, 2]
    st, a = call(S, "POST", "/api/lessons/99-sample/action", {"verb": "plot-history", "args": "res/a.ACC"})
    assert a["lines"] == ["READTH,a.ACC,0,3", "THPLOT,3"]               # pair 0: dt first (ACS SASSI format)
    assert call(S, "POST", "/api/command", {"lines": a["lines"]})[1]["ok"]
    st, a = call(S, "POST", "/api/lessons/99-sample/action", {"verb": "plot-soilprops", "args": "Sand"})
    assert a["lines"] == ["SOILPROPPLOT,Sand"]
    for verb, cmd in (("plot-model", "MODELPLOT"), ("plot-nodes", "NODEPLOT"), ("plot-layers", "LAYERPLOT")):
        assert call(S, "POST", "/api/lessons/99-sample/action", {"verb": verb})[1] == {"kind": "commands", "lines": [cmd]}
    # two files of the same name: each line is named by its workspace-relative path (legend, button label)
    (ws / "res2").mkdir()
    (ws / "res2" / "a01.RS").write_text("0.5 0.2\n20.0 0.4\n", encoding="utf-8")
    st, a = call(S, "POST", "/api/lessons/99-sample/action", {"verb": "plot-spectrum", "args": "res/a01.RS, res2/a01.RS"})
    assert a["lines"][1] == "LINENAME,4,res/a01.RS" and a["lines"][3] == "LINENAME,5,res2/a01.RS"
    assert call(S, "POST", "/api/command", {"lines": a["lines"]})[1]["ok"]
    assert S.plots.lines[5].name == "res2/a01.RS"
    assert learn.action_label("plot-spectrum", "res/a01.RS, res2/a01.RS | log") == "Plot res/a01.RS, res2/a01.RS (log f)"
    assert learn.action_label("plot-spectrum", "res/a01.RS, res/b01.RS") == "Plot a01.RS, b01.RS"
    # a file outside the working directory is named by its absolute path
    assert call(S, "POST", "/api/command", {"lines": [f"CD,{ws}", "MDL,m,mdir"]})[1]["ok"]
    st, a = call(S, "POST", "/api/lessons/99-sample/action", {"verb": "plot-spectrum", "args": "res/a01.RS"})
    assert a["lines"][0] == f"READSPEC,{ws / 'res' / 'a01.RS'},1,6"


def test_gui_request_actions(S):
    ws = _ws_with_results(S)
    st, a = call(S, "POST", "/api/lessons/99-sample/action", {"verb": "open-file", "args": "res/a01.RS"})
    assert a == {"kind": "file", "path": str(ws / "res" / "a01.RS")}
    st, a = call(S, "POST", "/api/lessons/99-sample/action", {"verb": "open-dialog", "args": "ANALYSIS/SITE"})
    assert a == {"kind": "dialog", "name": "ANALYSIS", "tab": "SITE"}
    assert call(S, "POST", "/api/lessons/99-sample/action", {"verb": "open-dialog", "args": "MODEL"})[1] == \
        {"kind": "dialog", "name": "MODEL", "tab": ""}
    assert call(S, "POST", "/api/lessons/99-sample/action", {"verb": "open-dialog", "args": "MOTION"})[1] == \
        {"kind": "dialog", "name": "ANALYSIS", "tab": "MOTION"}
    st, a = call(S, "POST", "/api/lessons/99-sample/action", {"verb": "open-doc", "args": "docs/theory/THEORY_MANUAL.md#2-complex-modulus-damping"})
    assert a == {"kind": "doc", "doc": "docs/theory/THEORY_MANUAL.md", "anchor": "2-complex-modulus-damping"}
    st, a = call(S, "POST", "/api/lessons/99-sample/action", {"verb": "explain", "args": "NGEN,3,5"})
    assert a == {"kind": "explain", "line": "NGEN,3,5"}


@pytest.mark.parametrize("verb,args,status", [
    ("open-file", "../../etc/passwd", 400), ("open-file", "/etc/passwd", 400), ("open-file", "C:/x.txt", 400),
    ("open-file", "res/../../x", 400), ("plot-spectrum", "../x.RS", 400), ("plot-spectrum", "", 400),
    ("plot-spectrum", "res/missing.RS", 404), ("open-file", "nothing.txt", 404), ("open-listing", "../SITE", 400),
    ("open-listing", "SITE", 404), ("open-dialog", "NOPE", 400), ("open-dialog", "ANALYSIS/NOPE", 400),
    ("open-doc", "docs/internal/lesson_format.md", 404), ("open-doc", "../README.md", 404), ("open-doc", "/etc/passwd", 404),
    ("plot-soilprops", "", 400), ("explain", "", 400), ("format-disk", "", 400)])
def test_action_path_safety_and_errors(S, verb, args, status):
    _ws_with_results(S)
    st, r = call(S, "POST", "/api/lessons/99-sample/action", {"verb": verb, "args": args})
    assert st == status, r


def test_every_contract_verb_is_resolved(S):
    _ws_with_results(S)
    args = {"plot-soilprops": "Sand", "plot-spectrum": "res/a01.RS", "plot-history": "res/a.ACC", "open-file": "res/a01.RS",
            "open-dialog": "ANALYSIS", "open-doc": "docs/index.md", "explain": "N,1", "open-listing": None,
            "animate": None}                 # needs frames: tests/unit/test_lesson_animate.py
    for verb in L.ACTION_VERBS:
        if args.get(verb, "") is None:
            continue
        st, r = call(S, "POST", "/api/lessons/99-sample/action", {"verb": verb, "args": args.get(verb, "")})
        assert st == 200 and r["kind"] in ("commands", "file", "dialog", "doc", "explain"), (verb, r)


# ====================================================================== examples
def test_examples_gallery(S):
    st, r = call(S, "GET", "/api/examples")
    assert st == 200
    ex = {x["name"]: x for x in r["examples"]}
    assert set(ex) >= {f"ex0{k}_{n}" for k, n in ((1, "surface_stick"), (2, "embedded_box"), (3, "forced_vibration"),
                                                  (4, "site_response"), (5, "xyz_simultaneous"), (6, "nonlinear_soil"),
                                                  (7, "option_non"))}
    e1 = ex["ex01_surface_stick"]
    assert e1["number"] == 1 and e1["title"] == "lumped-mass stick on a rigid surface mat"
    assert e1["learn"] and e1["learn"][0].startswith("the complete ACS SASSI chain SITE -> POINT")
    assert "ISRS" in " ".join(e1["learn"]) and e1["physics"].startswith("A vertically incident SV wave")
    assert e1["runtime"] == "~4 s" and "SITE POINT HOUSE" in e1["modules"] and e1["results"] == "ex01"
    e3 = ex["ex03_forced_vibration"]
    assert e3["title"].startswith("forced vibration of a surface foundation") and "global impedance option" in e3["title"]
    assert e3["lessons"] == [{"id": "98-example", "title": "Lesson on an example", "part": "Design applications"}]
    assert ex["ex04_site_response"]["learn"][0].startswith("EQUAKE: generate an acceleration history")
    for x in r["examples"]:
        assert x["title"] and x["learn"] and x["topic"] and x["runtime"], x["name"]


def test_parse_example_header_variants():
    h = learn.parse_example_header("* ===\n* SASSI-EDU tutorial example 9: a\n* two-line title\n* ===\n* What you learn\n"
                                   "*   1. first item\n*      continues\n*   2. second\n*\n* The model: a box\n*   of soil.\n"
                                   "* ===\nMDL,m9,out9\n")
    assert h == {"number": 9, "title": "a two-line title", "learn": ["first item continues", "second"],
                 "physics": "a box of soil.", "model": "m9", "results": "out9"}
    assert learn.readme_rows("| `exA.pre` | topic | SITE | ~1 s |\n| other |") == {"exA": {"topic": "topic", "modules": "SITE", "runtime": "~1 s"}}


def test_load_and_run_example_workspaces(S, tmp_path):
    st, r = call(S, "POST", "/api/examples/ex04_site_response/prepare", {"run": False})
    assert st == 200
    ws = Path(r["workspace"])
    assert ws == (tmp_path / "start" / "sassi-course" / "examples" / "ex04_site_response").resolve()
    assert r["pre"] == str(ws / "ex04_site_response.pre") and Path(r["pre"]).is_file()
    # data/ and the repository support file read through ../sassi/data/ (INP,../sassi/data/dynp_library.pre)
    assert (ws / "data" / "rg160h_030g.rsi").is_file() and (ws / "sassi" / "data" / "dynp_library.pre").is_file()
    assert r["commands"] == [f"CD,{ws}", "DMODEL,0"]
    run_like_the_player(S, r["commands"])
    assert Path(S.interp.cwd) == ws
    st, r = call(S, "POST", "/api/examples/ex03_forced_vibration/prepare", {"run": True})
    assert r["commands"][-1] == "INP,ex03_forced_vibration.pre"
    for bad in ("nosuch", "..", "ex01_surface_stick.pre"):
        assert call(S, "POST", f"/api/examples/{bad}/prepare", {"confirm": True})[0] == 404


def test_example_support_files():
    text = "INP,../sassi/data/dynp_library.pre\nTHFILE,../data/x.acc\nINP,../sassi/data/../../etc/passwd\n"
    assert L.example_support_files(text) == ["sassi/data/dynp_library.pre"]


@pytest.mark.slow
def test_run_all_of_an_example(S):
    """Run all: the copied example 3 runs with INP in its workspace without an error (about 1 s)."""
    st, r = call(S, "POST", "/api/examples/ex03_forced_vibration/prepare", {"run": True})
    st, res = call(S, "POST", "/api/command", {"lines": r["commands"]})
    assert res["ok"], [m["text"] for m in res["messages"] if m["kind"] == "ERROR"][:5]
    assert (Path(r["workspace"]) / "ex03" / "FOUNSTIF").exists()


def test_open_refused_while_a_module_runs(S):
    open_lesson(S, "97-run")
    st, d = call(S, "GET", "/api/lessons/97-run")
    run_like_the_player(S, d["steps"][0]["commands"])
    st, r = call(S, "POST", "/api/command", {"lines": ["RUNSITE"]})
    assert r.get("job")
    st, o = call(S, "POST", "/api/lessons/99-sample/open", {"confirm": True})
    st2, e = call(S, "POST", "/api/examples/ex01_surface_stick/prepare", {"confirm": True})
    _wait(S, r["job"]["id"])
    if st == 409:                                   # refused while the job was running (it may already be done)
        assert "still running" in o["error"]
    if st2 == 409:
        assert "still running" in e["error"]


# ====================================================================== settings
def test_course_root_setting(S, tmp_path):
    other = tmp_path / "elsewhere"
    other.mkdir()
    st, r = call(S, "POST", "/api/settings/course", {"course": {"root": str(other / "course")}})
    assert st == 200 and r["settings"]["course"]["root"] == str((other / "course").resolve())
    assert IniSettings(tmp_path / "settings").course["root"] == str((other / "course").resolve())   # SASSIini.xml
    st, state = call(S, "GET", "/api/state")
    assert state["course_root"] == str((other / "course").resolve())
    ws = Path(open_lesson(S, "99-sample")["workspace"])
    assert ws.parent == (other / "course").resolve()
    (ws / "notes.txt").write_text("hello")
    assert call(S, "GET", "/api/file", query={"name": str(ws / "notes.txt")})[1]["text"] == "hello"
    # refused: a file, a folder whose parent does not exist, the filesystem root, the home folder
    f = tmp_path / "afile"
    f.write_text("x")
    for bad in (str(f), str(tmp_path / "no" / "such" / "dir"), "/", str(Path.home())):
        st, r = call(S, "POST", "/api/settings/course", {"course": {"root": bad}})
        assert st == 400, bad
    # blank: back to the default <start>/sassi-course
    st, r = call(S, "POST", "/api/settings/course", {"course": {"root": ""}})
    assert r["settings"]["course"]["root"] == "" and S.course_root() == (tmp_path / "start" / "sassi-course").resolve()


def test_fresh_model_lines():
    assert learn.fresh_model_lines([0], 0) == ["DMODEL,0"]
    assert learn.fresh_model_lines([0, 2, 5], 2) == ["ACTM,0", "DMODEL,2", "DMODEL,5", "DMODEL,0"]
    assert learn.fresh_model_lines([3], 3) == ["ACTM,0", "DMODEL,3", "DMODEL,0"]
    assert learn.cd_line(Path("/a, b/c")) == 'CD,"/a, b/c"'


# ====================================================================== Learn > Free Disk Space
def test_free_disk_space_deletes_only_unused_course_workspaces(S, tmp_path):
    root = tmp_path / "start" / "sassi-course"
    open_lesson(S, "98-example")                                   # two lesson workspaces, then an example
    open_lesson(S, "99-sample")                                    # (the open lesson of the session)
    st, r = call(S, "POST", "/api/examples/ex01_surface_stick/prepare", {"confirm": True})
    assert st == 200
    run_like_the_player(S, r["commands"])                          # the session now works in the example
    # a folder under the course root that the course did not create is neither listed nor deleted
    mine = root / "my-notes"
    mine.mkdir()
    (mine / "keep.txt").write_text("keep")
    st, d = call(S, "GET", "/api/course/workspaces")
    assert st == 200 and d["root"] == str(root.resolve())
    by = {w["name"]: w for w in d["workspaces"]}
    assert set(by) == {"98-example", "99-sample", "examples/ex01_surface_stick"}
    assert by["98-example"]["bytes"] > 0 and not by["98-example"]["in_use"]
    assert by["99-sample"]["in_use"]                               # the open lesson
    assert by["examples/ex01_surface_stick"]["in_use"]             # the working directory
    assert d["bytes"] == sum(w["bytes"] for w in d["workspaces"])
    st, r = call(S, "POST", "/api/course/workspaces/delete", {"names": ["98-example", "examples/ex01_surface_stick",
                                                                      "my-notes", "../start"]})
    assert st == 200 and r["deleted"] == ["98-example"] and r["freed"] == by["98-example"]["bytes"]
    assert {s["name"] for s in r["skipped"]} == {"examples/ex01_surface_stick", "my-notes", "../start"}
    assert not (root / "98-example").exists() and (root / "examples" / "ex01_surface_stick").is_dir()
    assert (mine / "keep.txt").read_text() == "keep"
    assert sorted(w["name"] for w in r["workspaces"]) == ["99-sample", "examples/ex01_surface_stick"]
    # with no names: every workspace not in use; a deleted lesson workspace is recreated when opened again
    assert call(S, "POST", "/api/course/workspaces/delete", {"names": "98-example"})[0] == 400
    open_lesson(S, "98-example")
    st, r = call(S, "POST", "/api/course/workspaces/delete", {})
    assert sorted(r["deleted"]) == ["99-sample", "examples/ex01_surface_stick"] and (root / "98-example").is_dir()
    assert (root / "98-example" / "ex03_forced_vibration.pre").is_file()
