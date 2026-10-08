"""The Properties panel (sassi/ui/properties.py, static/props.js; requirements 7.20, D-W6-16).

Nodes and elements are selected with command text (NODESEL with ranges, ELEMSEL, SELCLR -- a click in a plot
submits them); the panel shows their properties and edits them with command text; the input file that built
the model is changed so that it builds the edited model (a node's own N / MT / MR line rewritten in place,
every other change in an edits section after the part of the file that last changed the model).  The edited
example 1 rebuilt from its changed input file is the edited model."""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from sassi.prep import Interpreter
from sassi.ui import properties as PR

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"
NODE = shutil.which("node")


@pytest.fixture(autouse=True)
def _settings(tmp_path, monkeypatch):
    monkeypatch.setenv("SASSI_EDU_SETTINGS_DIR", str(tmp_path / "settings"))


def _session(tmp_path, name="ex01_surface_stick"):
    from sassi.ui.api import GuiSession
    from sassi.ui.learn import example_model_lines
    from sassi.ui.lessons import EXAMPLES_DIR, copy_example
    pre = copy_example(name, tmp_path / "ex", EXAMPLES_DIR)
    (pre.parent / "model.pre").write_text("\n".join(example_model_lines(pre.read_text())) + "\n", encoding="utf-8")
    S = GuiSession(cwd=str(pre.parent), settings_dir=str(tmp_path / "settings"), web=True)
    st, r = S.route("POST", "/api/command", {}, {"lines": ["INP,model.pre"]})
    assert r["ok"]
    return S


def _model(S):
    m = S.interp.model
    return (m.nodes[82].z, list(m.tmass[82]), list(m.nodes[1].fix), sorted(m.nodes[2].flags), m.groups[1].elements[5].mat,
            m.groups[1].elements[6].thick, list(m.groups[2].elements[1].kj), list(m.rmass.get(85, [])), m.group_active)


def test_select_edit_and_rebuild_from_the_changed_input(tmp_path):
    S = _session(tmp_path)
    st, r = S.route("POST", "/api/command", {}, {"lines": ["NODESEL,82,1-3", "ELEMSEL,1,5-6", "ELEMSEL,2,1"]})
    assert r["ok"]
    st, sel = S.route("GET", "/api/selection", {}, {})
    assert st == 200 and [n["id"] for n in sel["nodes"]] == [82, 1, 2, 3]
    n82 = sel["nodes"][0]
    assert n82["xyz"] == [0.0, 0.0, 16.0] and n82["mass"] == [2200.0, 2200.0, 2200.0] and n82["elements"] == [[2, 1], [2, 2]]
    assert [(g["id"], g["type_name"], [e["id"] for e in g["elements"]]) for g in sel["groups"]] == [(1, "SHELL", [5, 6]), (2, "BEAMS", [1])]
    assert sel["groups"][0]["elements"][0]["thick"] == 5.0 and [m["id"] for m in sel["materials"]] == [1, 2, 3]
    assert sel["input"]["name"] == "model.pre" and sel["active_group"] == 3
    edits = [{"op": "xyz", "nodes": [82], "z": 17},
             {"op": "mass", "nodes": [82], "values": [2500, 2500, None]},
             {"op": "rmass", "nodes": [85], "values": [10, None, None]},
             {"op": "fix", "nodes": [1, 2, 3], "dofs": {"UX": 1}},
             {"op": "interaction", "nodes": [2, 3], "value": 0},
             {"op": "elem", "group": 1, "elements": [5, 6], "attr": "mat", "value": 2},
             {"op": "elem", "group": 1, "elements": [5, 6], "attr": "thick", "value": 4.5},
             {"op": "release", "group": 2, "elements": [1], "end": "J", "codes": [0, 0, 0, 0, 1, 1]}]
    st, res = S.route("POST", "/api/properties", {}, {"edits": edits, "input": {"path": sel["input"]["path"]}})
    assert st == 200 and res["ok"], res
    assert res["lines"] == ["N,82,0,0,17", "MT,82,2500,2500,2200", "MR,85,10,0,0", "D,1,3,1,1,UX", "INT,2,3,1,0,0",
                            "GROUP,1", "MSET,5,6,1,2", "THICK,5,6,1,4.5", "GROUP,2", "KJ,1,1,1,0,0,0,0,1,1", "GROUP,3"]
    inp = res["input"]
    assert inp["saved"] and inp["inplace"] == 1 and inp["added"] == 10
    text = Path(sel["input"]["path"]).read_text()
    assert text == inp["text"]
    lines = text.splitlines()
    assert "N,82,0,0,17" in lines and "N,82,0,0,16" not in lines                # rewritten in place
    k = lines.index(PR.EDITS_HEADER)
    assert lines[k - 1] == "*" and lines[k + 1].startswith("* model ex01:") and lines[k + 2] == "MT,82,2500,2500,2200"
    assert [lines[i] for i in inp["changed"]][0] == "N,82,0,0,17"
    want = _model(S)
    assert want == (17.0, [2500.0, 2500.0, 2200.0], [1, 0, 0, 0, 0, 1], [], 2, 4.5, [0, 0, 0, 0, 1, 1], [10.0, 0.0, 0.0], 3)
    # the changed input file builds the edited model
    from sassi.ui.api import GuiSession
    S2 = GuiSession(cwd=str(tmp_path / "ex"), settings_dir=str(tmp_path / "settings2"), web=True)
    st, r2 = S2.route("POST", "/api/command", {}, {"lines": [f"INP,{sel['input']['path']}"]})
    assert r2["ok"] and _model(S2) == want
    # a second edit goes into the same section; a node's MT line replaces its previous one
    st, res = S.route("POST", "/api/properties", {}, {"edits": [{"op": "mass", "nodes": [82], "values": [None, None, 3000]}],
                                                       "input": {"path": sel["input"]["path"]}})
    lines = res["input"]["text"].splitlines()
    assert lines.count(PR.EDITS_HEADER) == 1 and "MT,82,2500,2500,3000" in lines and "MT,82,2500,2500,2200" not in lines


def test_a_file_editor_with_unsaved_edits_gets_the_change_in_its_buffer(tmp_path):
    S = _session(tmp_path)
    st, sel = S.route("GET", "/api/selection", {}, {})
    path = Path(sel["input"]["path"])
    before = path.read_text()
    buffer = before + "* an unsaved note\n"
    st, res = S.route("POST", "/api/properties", {}, {"edits": [{"op": "xyz", "nodes": [90], "x": 1}],
                                                       "input": {"path": str(path), "text": buffer}})
    assert st == 200 and not res["input"]["saved"] and path.read_text() == before                 # not written
    assert "N,90,1,0,8" in res["input"]["text"] and res["input"]["text"].endswith("* an unsaved note\n")


def test_errors_and_the_model_only(tmp_path):
    S = _session(tmp_path)
    assert S.route("POST", "/api/properties", {}, {"edits": []})[0] == 400
    assert S.route("POST", "/api/properties", {}, {"edits": [{"op": "xyz", "nodes": [999], "x": 1}]})[0] == 400
    assert S.route("POST", "/api/properties", {}, {"edits": [{"op": "elem", "group": 1, "elements": [1], "attr": "thick", "value": -1}]})[0] == 400
    st, res = S.route("POST", "/api/properties", {}, {"edits": [{"op": "fix", "nodes": [5], "dofs": {"UZ": 1}}]})
    assert st == 200 and res["ok"] and "input" not in res and S.interp.model.nodes[5].fix[2] == 1    # no input: the model only


def test_patch_input_placement_and_local_systems():
    ui = Interpreter()
    for ln in ["MDL,m,m", "N,1,0,0,0", "N,2,10,0,0", "FILL,1,2", "GROUP,1,BEAMS", "E,1,1,2", "MT,2,5,5,5"]:
        assert ui.execute(ln)
    m = ui.model
    text = "MDL,m,m\nN,1,0,0,0\nN,2,10,0,0\nGROUP,1,BEAMS\nE,1,1,2\nMT,2,5,5,5\n*\nCHECK\nAFWRITE\nRUNSITE\n"
    new, info = PR.patch_input(text, m, [{"op": "xyz", "nodes": [2], "y": 3}, {"op": "fix", "nodes": [1], "dofs": {"UZ": 1}}])
    lines = new.splitlines()
    assert lines[2] == "N,2,10,3,0" and info["inplace"] == 1                                   # in place
    k = lines.index(PR.EDITS_HEADER)                                                          # before CHECK's section
    assert lines[k + 2] == "D,1,1,1,1,UZ" and lines.index("CHECK") > k and lines[k - 1] == "*"
    # a node defined twice, or in a local system: the change goes into the section (N in the global system)
    text2 = text.replace("MT,2,5,5,5", "N,2,10,0,0\nMT,2,5,5,5")
    new2, _ = PR.patch_input(text2, m, [{"op": "xyz", "nodes": [2], "z": 1}])
    assert new2.count("N,2,10,0,0") == 2 and "N,2,10,0,1" in new2.splitlines()[new2.splitlines().index(PR.EDITS_HEADER):]
    assert PR.ranges([3, 1, 2, 7]) == [(1, 3), (7, 7)]
    assert PR.anchor_fallback(["N,1,0,0,0", "*", "* run", "RUNSITE"]) == 2 and PR.anchor_fallback(["N,1,0,0,0"]) == 1


def test_selection_commands_and_the_scene():
    from sassi.plotting.state import model_scene
    ui = Interpreter()
    for ln in ["MDL,m,m", "N,1,0,0,0", "N,4,3,0,0", "FILL,1,4", "NGEN,2,4,1,4,1,0,1,0", "GROUP,1,SHELL", "E,1,1,2,6,5",
               "EGEN,3,1,1"]:
        assert ui.execute(ln)
    m = ui.model
    assert ui.execute("NODESEL,1-3,8") and m.ui_state["nodesel"] == [1, 2, 3, 8]
    assert ui.execute("NODESEL,2") and m.ui_state["nodesel"] == [1, 3, 8]                     # toggles
    assert ui.execute("ELEMSEL,1,1-2") and m.ui_state["elemsel"] == [[1, 1], [1, 2]]
    assert not ui.execute("ELEMSEL,9,1")                                                     # no group 9
    sc = model_scene(m)
    assert [int(sc["elem_id"][k]) for k in sc["elem_selected"]] == [1, 2]
    assert ui.execute("SELCLR,ELEM") and m.ui_state["elemsel"] == [] and m.ui_state["nodesel"] == [1, 3, 8]
    assert ui.execute("SELCLR") and m.ui_state["nodesel"] == []
    assert not ui.execute("SELCLR,X")
    explain = json.loads((ROOT / "sassi" / "ui" / "command_args.json").read_text(encoding="utf-8"))
    assert {"NODESEL", "ELEMSEL", "SELCLR"} <= set(explain)


@pytest.mark.skipif(NODE is None, reason="Node.js is not installed")
def test_browser_helpers():
    code = """
      globalThis.SASSI = {el: () => ({}), editors: {}};
      require(process.argv[1]);
      const R = SASSI.Props;
      console.log(JSON.stringify({r: R.ranges([7, 1, 2, 3, 9, 10]), p: R.parseIds("1-3, 7 9;12"), bad: R.parseIds("1-x"),
        lines: R.selectLines("ELEMSEL,1", [1, 2, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41], 19)}));"""
    p = subprocess.run([NODE, "-e", code, str(STATIC / "props.js")], capture_output=True, text=True, timeout=60)
    assert p.returncode == 0, p.stderr[-2000:]
    res = json.loads(p.stdout)
    assert res["r"] == "1-3, 7, 9-10" and res["p"] == [1, 2, 3, 7, 9, 12] and res["bad"] is None
    assert res["lines"] == ["ELEMSEL,1,1-3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39", "ELEMSEL,1,41"]   # 19 tokens


def test_front_end_wiring():
    page = (STATIC / "index.html").read_text(encoding="utf-8")
    assert '<aside id="propdock" hidden aria-label="Properties"></aside>' in page
    assert page.index('src="static/plots.js"') < page.index('src="static/props.js"') < page.index('src="static/main.js"')
    app = (STATIC / "app.js").read_text(encoding="utf-8")
    assert "if (SASSI.Props) SASSI.Props.refresh();" in app and '{label: "Properties", check: () =>' in app
    assert "setText: (text, saved) =>" in app and "flash: (lines) =>" in app
    plots = (STATIC / "plots.js").read_text(encoding="utf-8")
    assert 'S.command(only ? ["SELCLR"] : ["SELCLR", `NODESEL,${id}`]);' in plots and "S.command(`ELEMSEL,${g},${e}`)" in plots
    assert "function pickAt(div, d, cx, cy, nodes)" in plots and "selectedTraces(sc, true)" in plots
    js = (STATIC / "props.js").read_text(encoding="utf-8")
    assert 'S.post("/api/properties", {model: d.model, edits, input})' in js and 'S.get("/api/selection")' in js
    css = (STATIC / "styles.css").read_text(encoding="utf-8")
    assert "#propdock" in css and ".editor-hl .hl-chg" in css
