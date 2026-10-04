"""End-to-end GUI flows through the JSON API (requirements 5.3-5.9; package gui_polish).

Every flow below is the command text / API sequence the browser front end submits for a menu item,
toolbar button or dialog -- the browser walk-through of docs/user/GUI.md was also done by hand
against a running ``sassi-gui`` -- so the GUI behaviour is pinned by tests that need no browser:

* Options > Analysis: every tab round-trips (open, change a field, OK, reopen shows the new value,
  WRITE writes the command); Options > Model / Write / Check;
* Plot > Spectrum TFU-TFI / Time History (Line Selection): lines added with Add Line(s) are those
  READSPEC wrote, the dialog pre-fills from the active plot (``GET /api/plots``) and its command text
  reproduces that plot;
* Model > Converters with the optional output .pre file; Model > Exit;
* Process Animation Frame List on a module frame folder, the four animations, Pause keeping the frame;
* menus and toolbars of requirements 5.3 / 5.9 present in the front end, and static regression
  guards for the front-end fixes (no native prompt/confirm, responsive dialogs, Help pages).
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

from sassi.ui import api as apimod
from sassi.ui.api import GuiSession

STATIC = Path(__file__).resolve().parents[2] / "sassi" / "ui" / "static"
ROOT = Path(__file__).resolve().parents[2]

MODEL = """MDL,box,{dir}
N,1,0,0,0
N,2,2,0,0
N,3,2,2,0
N,4,0,2,0
N,5,0,0,-2
N,6,2,0,-2
N,7,2,2,-2
N,8,0,2,-2
INT,1,8,1,1
D,5,8,1,1,ALL
M,1,3e6,0.25,0.15,0.05,0.05
GROUP,1,SOLID
E,1,5,6,7,8,1,2,3,4
GROUP,2,BEAMS
R,1,1,1,1,1,1,1
E,1,1,3,2
MT,3,10,10,10
L,1,10,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
TOPL,1
"""


@pytest.fixture
def S(tmp_path):
    s = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    yield s
    run = s.jobs.running()
    if run is not None:
        s.jobs.cancel(run.id)


def call(S, method, path, body=None, query=None):
    q = {k: [str(v)] for k, v in (query or {}).items()}
    return S.route(method, path, q, body or {})


def cmd(S, *lines):
    st, r = call(S, "POST", "/api/command", {"lines": list(lines)})
    assert st == 200 and r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    return r


def load_model(S, tmp_path):
    cmd(S, *[ln for ln in MODEL.format(dir=tmp_path / "box").splitlines() if ln.strip()])
    return tmp_path / "box"


def written(S, tmp_path, name="check.pre"):
    cmd(S, f"WRITE,{name},{tmp_path}")
    return (tmp_path / name).read_text(encoding="utf-8").splitlines()


# ================================================================== Options dialogs (5.4)
#: one field per tab of Options > Analysis: (tab, record payload, record, field, value)
TAB_EDITS = [
    ("EQUAKE", {"EQUAKE": {"damp": 0.07}}, "EQUAKE", "damp", 0.07),
    ("SOIL", {"SOIL": {"nrval": 4000, "iter": 5}}, "SOIL", "iter", 5),
    ("SITE", {"SITE": {"nl": 12, "hs": 2}}, "SITE", "nl", 12),
    ("POINT", {"POINT": {"rad": 1.5}}, "POINT", "rad", 1.5),
    ("HOUSE", {"HOUSE": {"gelev": -1.0}}, "HOUSE", "gelev", -1.0),
    ("FORCE", {"FORCE": {"opmode": 1}}, "FORCE", "opmode", 1),
    ("ANALYS", {"ANALYS": {"impe": 2}}, "ANALYS", "impe", 2),
    ("MOTION", {"MOTION": {"dur": 5.0}}, "MOTION", "dur", 5.0),
    ("STRESS", {"STRESS": {"itran": 1}}, "STRESS", "itran", 1),
    ("RELDISP", {"RELD": {"reldispsall": 1}}, "RELD", "reldispsall", 1),
    ("AFWRITE", {"AOPT": {"equake": 1}}, "AOPT", "equake", 1),
]


@pytest.mark.parametrize("tab,records,rec,field,value", TAB_EDITS, ids=[t[0] for t in TAB_EDITS])
def test_every_analysis_tab_round_trips_and_writes(S, tmp_path, tab, records, rec, field, value):
    """Open (GET), change a field, OK (POST), reopen shows the new value, WRITE writes the command."""
    load_model(S, tmp_path)
    st, before = call(S, "GET", f"/api/options/{tab}")
    assert st == 200 and [t["name"] for t in before["form"]["tabs"]] == [tab]
    assert before["values"]["records"][rec][field] != value
    st, r = call(S, "POST", "/api/options/ANALYSIS", {"records": records})
    assert st == 200 and r["ok"], r
    line = [c for c in r["commands"] if c.startswith(rec + ",")]
    assert len(line) == 1, r["commands"]
    st, after = call(S, "GET", f"/api/options/{tab}")
    assert after["values"]["records"][rec][field] == value
    assert line[0] in written(S, tmp_path)                      # WRITE shows the command (L17, UT-03)


def test_lists_strings_requests_round_trip(S, tmp_path):
    load_model(S, tmp_path)
    payload = {"lists": {"DAMP": "0.02 0.05 0.10", "TOPL": "1 1 2"},
               "strings": {"THFILE": "motion.acc", "THTIT": "control, motion"},
               "indexed": {"WAVE": {"3": {"type": 3, "opt": 1, "ratio1": 1, "ratio2": 1, "angle": 0}}},
               "requests": {"NOUT": [{"dir": 1, "c1": 1, "c2": 1, "c3": 0, "c4": 0, "c5": 1, "c6": 1, "nodes": "1-3 5"}]}}
    st, r = call(S, "POST", "/api/options/ANALYSIS", payload)
    assert st == 200 and r["ok"], r
    st, v = call(S, "GET", "/api/options/ANALYSIS")
    v = v["values"]
    assert v["lists"]["DAMP"] == "0.02 0.05 0.1" and v["lists"]["TOPL"] == "1 1 2"
    assert v["strings"]["THTIT"] == "control, motion" and v["indexed"]["WAVE"]["3"]["opt"] == 1
    assert v["requests"]["NOUT"][0]["nodes"] == "1-3 5"
    # the first WAVE page also stores the implicit vertical SV field (shown = used, GUI.md section 4)
    assert r["commands"] == ["WAVE,2,1,1,1,0", "WAVE,3,1,1,1,0", "THFILE,motion.acc", "THTIT,control, motion", "DAMP,0",
                             "DAMP,0.02,0.05,0.1", "TOPL,0", "TOPL,1,1,2", "NOUT,0", "NOUT,1,1,1,0,0,1,1,1-3,5"]
    lines = written(S, tmp_path)
    assert all(c in lines for c in r["commands"])                    # WRITE writes every committed command


#: the command-record parts of the Analysis window (SOIL Nonlinear Soil, NONLINEAR tab) and the Modules-menu
#: dialogs: (dialog GET, dialog POST, payload part, family, field, value, command it emits)
X_EDITS = [
    ("SOIL", "ANALYSIS", "xrecords", "NLSOIL", "opt", 1, "NLSOIL,1,0,0,0,0,0,0,0,0"),
    ("SOIL", "ANALYSIS", "xtables", "NLSLAYER", None, [{"num": 1, "curvefit": 1}], "NLSLAYER,1,1,0,0,0,0"),
    ("NONLINEAR", "ANALYSIS", "xrecords", "EQL", "disp", 0.75, "EQL,0.75,,0,0,0"),
    ("NONLINEAR", "ANALYSIS", "xtables", "P", None, [{"num": 1, "group": 1, "bbc": 1, "disp": 1, "force": 1}],
     "P,1,1,1,1,1"),
    ("NONLINEAR", "ANALYSIS", "xindexed", "BBC", "1", {"type": 4, "yield": 1, "points": [{"x": 0.01, "y": 10}]},
     "BBCX,1,1,1,0.01"),
    ("LOADGEN", "LOADGEN", "xrecords", "LOADGEN", "data", 3, "LOADGEN,3,0,0,0,1,0,RESULTS"),
    ("LOADGEN", "LOADGEN", "xrecords", "LGTIME", "n", 3, "LGTIME,V,3,0"),
    ("LOADGENDYN", "LOADGENDYN", "xrecords", "LOADGENDYN", "alpha", 0.5, "LOADGENDYN,0.5,0,REL,0,RESULTS,0,0,0,0,0,0,1"),
    ("LOADGENDYN", "LOADGENDYN", "xrecords", "LGFILE", "GROUND", "kin.acc", "LGFILE,GROUND,kin.acc"),
]


@pytest.mark.parametrize("tab,dialog,part,fam,field,value,command", X_EDITS,
                         ids=[f"{t[0]}-{t[3]}" for t in X_EDITS])
def test_command_record_dialogs_round_trip_and_write(S, tmp_path, tab, dialog, part, fam, field, value, command):
    """The new bindings through the API as the browser drives them: open (GET), change, OK (POST: the command
    text), reopen shows the value, WRITE writes the command (L17, UT-03)."""
    load_model(S, tmp_path)
    st, before = call(S, "GET", f"/api/options/{tab}")
    assert st == 200 and before["form"]["tabs"][0]["name"] == tab
    v = before["values"]
    if part == "xrecords":
        body = {part: {fam: dict(v[part][fam], **{field: value})}}
    elif part == "xtables":
        body = {part: {fam: value}}
    else:
        body = {part: {fam: {field: value}}}
    st, r = call(S, "POST", f"/api/options/{dialog}", body)
    assert st == 200 and r["ok"], r
    assert command in r["commands"], r["commands"]
    after = call(S, "GET", f"/api/options/{tab}")[1]["values"]
    if part == "xrecords":
        assert after[part][fam][field] == value
    elif part == "xtables":
        assert [row["num"] for row in after[part][fam]] == [row["num"] for row in value]
    else:
        assert after[part][fam][field]["type"] == value["type"]
    assert command in written(S, tmp_path)
    # posting the reopened values back changes nothing (OK without a change emits no command)
    keys = ("xrecords", "xtables", "xindexed")
    st, r = call(S, "POST", f"/api/options/{dialog}", {k: after[k] for k in keys if k in after})
    assert st == 200 and r["commands"] == []


def test_command_record_commit_refused_through_the_api(S, tmp_path):
    load_model(S, tmp_path)
    st, r = call(S, "POST", "/api/options/ANALYSIS", {"xrecords": {"EQL": {"disp": 3}}, "records": {"SITE": {"nl": 12}}})
    assert st == 422 and any("equivalent-linear displacement factor" in p for p in r["problems"])
    assert call(S, "GET", "/api/options/SITE")[1]["values"]["records"]["SITE"]["nl"] != 12      # nothing applied
    st, r = call(S, "POST", "/api/options/LOADGEN", {"xrecords": {"LGNODE": {"D": "a-b"}}, "dry_run": True})
    assert st == 422 and r["problems"][0].startswith("LGNODE: D nodes")
    st, r = call(S, "POST", "/api/options/LOADGEN", {"xrecords": {"LGOPT": {"digits": 9}}, "dry_run": True})
    assert st == 200 and r["commands"] == ["LGOPT,9,0,1"] and S.interp.model.options.record("LGOPT") is None
    assert call(S, "POST", "/api/options/NOSUCH", {})[0] == 404


def test_model_write_check_dialogs_round_trip(S, tmp_path):
    load_model(S, tmp_path)
    st, r = call(S, "POST", "/api/options/MODEL", {"records": {"MOPT": {"matrix": 1}}})
    assert st == 200 and r["commands"] == ["MOPT,1,1,1,1"]
    assert call(S, "GET", "/api/options/MODEL")[1]["values"]["records"]["MOPT"]["matrix"] == 1
    st, r = call(S, "POST", "/api/options/WRITE", {"records": {"WRITE": {"mdl": 1}}})
    assert st == 200 and call(S, "GET", "/api/options/WRITE")[1]["values"]["records"]["WRITE"]["mdl"] == 1
    lines = written(S, tmp_path)
    assert any(ln.startswith("MDL,box,") for ln in lines)                  # "MDL command in *.Pre"
    assert "MOPT,1,1,1,1" in lines
    st, r = call(S, "POST", "/api/options/CHECK", {"records": {"CHECK": {"break_at": 50}}})
    assert call(S, "GET", "/api/options/CHECK")[1]["values"]["records"]["CHECK"]["break_at"] == 50


# ================================================================== Line Selection (5.7, spec 06 5.2)
def _spec_file(path: Path, scale: float) -> Path:
    path.write_text("# f SA\n" + "".join(f"{f} {scale * f / (1 + f)}\n" for f in (0.1, 0.5, 1, 2, 5, 10, 20)))
    return path


def test_line_selection_adds_check_and_prefill(S, tmp_path):
    """The dialog's command sequence: Add Line(s) = READSPEC into the Starting Number (the lines a
    later 'lines' refresh shows are exactly those READSPEC wrote), Ok = SPECPLOT + PLOTTITLE + XTITLE +
    AXES + PLOTRANGE; reopening pre-fills from GET /api/plots (lines, titles, log axes, ranges)."""
    a, b = _spec_file(tmp_path / "a.RS", 1), _spec_file(tmp_path / "b.RS", 2)
    r = cmd(S, f"READSPEC,{a},1,1")
    r = cmd(S, f"READSPEC,{b},1,2")                                   # the second file: line 2
    assert sorted(int(k) for k in call(S, "GET", "/api/lines")[1]["lines"]) == [1, 2]
    cmd(S, "SPECPLOT,1,2", "PLOTTITLE,Node 41, X direction", "XTITLE,Frequency (Hz)", "AXES,1,1,0,0,1,1",
        "PLOTRANGE,0.1,50,,")
    st, ps = call(S, "GET", "/api/plots")
    act = next(p for p in ps["plots"] if p["id"] == ps["active"])
    assert act["kind"] == "SPECPLOT" and act["params"]["lines"] == [1, 2]
    assert act["title"] == "Node 41, X direction"
    s = act["settings"]
    assert (s["xtitle"], s["log_x"], s["log_y"], s["xmin"], s["xmax"], s["ymin"]) == ("Frequency (Hz)", True, True,
                                                                                       0.1, 50.0, None)
    d2 = ps["defaults2d"]                             # what a new plot starts from: the dialog submits the differences
    assert (d2["log_x"], d2["xtitle"], d2["xmin"]) == (False, "", None)
    # reopen + Add Line(s) + Replace: the same settings on the replacing plot, with the new line
    _spec_file(tmp_path / "c.RS", 3)
    cmd(S, f"READSPEC,{tmp_path / 'c.RS'},1,3", "CLOSEPLOT", "SPECPLOT,1,2,3", "PLOTTITLE,Node 41, X direction",
        "XTITLE,Frequency (Hz)", "AXES,1,1,0,0,1,1", "PLOTRANGE,0.1,50,,")
    st, ps = call(S, "GET", "/api/plots")
    assert len(ps["plots"]) == 1
    p = ps["plots"][0]
    assert p["params"]["lines"] == [1, 2, 3] and p["settings"] == s and p["title"] == act["title"]
    st, d = call(S, "GET", f"/api/plot/{p['id']}")
    assert [L["number"] for L in d["lines"]] == [1, 2, 3]


def test_time_history_line_selection(S, tmp_path):
    h = tmp_path / "x.ACC"
    h.write_text("0.01\n" + "".join(f"{v}\n" for v in (0, 0.1, -0.2, 0.3, 0)))
    st, fi = call(S, "GET", "/api/fileinfo", query={"name": "x.ACC"})
    assert fi["plot"] == "th" and fi["pair"] == 0 and fi["points"] == 5
    cmd(S, f"READTH,{h},0,4", "THPLOT,4", "YTITLE,acc (g)")
    st, d = call(S, "GET", "/api/plot/1")
    assert d["kind"] == "THPLOT" and d["lines"][0]["x"][-1] == pytest.approx(0.04) and d["settings"]["ytitle"] == "acc (g)"


# ================================================================== Model menu
def test_converters_with_output_pre_file(S, tmp_path):
    """Model > Converters (spec 09 6.4): blank model number = the active model; Output .pre File Name.
    The decks come from tutorial example 1 (INP without its RUN lines, then AFWRITE)."""
    src = (ROOT / "examples" / "ex01_surface_stick.pre").read_text(encoding="utf-8").splitlines()
    (tmp_path / "ex01_decks.pre").write_text("\n".join(ln for ln in src if not ln.upper().startswith(("RUN", "WRITE"))))
    cmd(S, "INP,ex01_decks.pre")
    hou = tmp_path / "ex01" / "ex01.hou"
    assert hou.is_file()
    cmd(S, f"CONVERT,SSI,2,{hou},{tmp_path / 'from_hou.pre'}")
    assert (tmp_path / "from_hou.pre").is_file()
    assert S.interp.models[2].counts()["nodes"] == S.interp.models[0].counts()["nodes"] == 90
    cdb = ROOT / "tests" / "data" / "ansys" / "beam188_cantilever.cdb"
    if cdb.is_file():
        cmd(S, "ACTM,3", f"CONVERT,ANSYS,,{cdb},9.81,{tmp_path / 'beam.pre'}")         # blank = active (3)
        assert S.interp.models[3].counts()["elements"] > 0 and (tmp_path / "beam.pre").is_file()
        cmd(S, "ACTM,4", f"INP,{tmp_path / 'beam.pre'}")                               # Model > Input of it
        assert S.interp.models[4].counts()["elements"] == S.interp.models[3].counts()["elements"]


def test_model_open_save_and_exit(S, tmp_path):
    load_model(S, tmp_path)
    assert call(S, "GET", "/api/exit")[1]["unsaved"] == [{"number": 0, "name": "box"}]
    cmd(S, "SAVE")
    call(S, "POST", "/api/db", {"action": "add_group", "group": "Tutorials"})
    call(S, "POST", "/api/db", {"action": "add_model", "group": "Tutorials", "name": "box",
                                "path": str(tmp_path / "box"), "title": "Box, test"})
    cmd(S, "ACTM,1")
    st, r = call(S, "POST", "/api/db", {"action": "open", "name": "box", "path": str(tmp_path / "box")})
    assert st == 200 and r["ok"] and S.interp.history[-2:] == [f"MDL,box,{tmp_path / 'box'}", "RESUME"]
    assert S.interp.model.counts()["nodes"] == 8
    assert call(S, "GET", "/api/exit")[1]["unsaved"] == []
    st, r = call(S, "POST", "/api/exit")
    assert r["ok"] and (tmp_path / "settings" / "SASSIini.xml").is_file()


# ================================================================== animations (5.8)
def _frames(folder: Path, n: int = 3) -> Path:
    folder.mkdir()
    for k in range(1, n + 1):
        rows = "".join(f"{node} {0.01 * k * node} 0 {-0.002 * k}\n" for node in (1, 2, 3, 4))
        (folder / f"THD_{k * 0.01:06.3f}_{k:05d}").write_text(f"4 4\n{rows}")
    return folder


def test_procframe_on_a_module_frame_folder_and_animations(S, tmp_path):
    """Process Animation Frame List with a frame folder (MOTION / RELDISP restart frames), then the
    four animations; Pause keeps the frame shown (PAUSE + WINDOWSETTINGS,FRAME,<k>)."""
    load_model(S, tmp_path)
    folder = _frames(tmp_path / "box" / "THD")
    cmd(S, f"PROCFRAME,{folder},ani_thd,box THD frames,0")
    st, a = call(S, "GET", "/api/animations")
    ent = [e for e in a["entries"] if e["description"] == "box THD frames"]
    assert ent and ent[0]["frames"] == "3"
    store = ent[0]["directory"]
    for line, kind in ((f"DEFORMPLOT,{store},1,3,1,50", "DEFORMPLOT"), (f"BUBBLEPLOT,{store},1,3,1,,,1", "BUBBLEPLOT"),
                       (f"CONTOURPLOT,{store},1,3,1,,,1", "CONTOURPLOT"), (f"VECTORPLOT,{store},1,3,1,2", "VECTORPLOT")):
        cmd(S, line)
        st, ps = call(S, "GET", "/api/plots")
        act = next(p for p in ps["plots"] if p["id"] == ps["active"])
        assert act["kind"] == kind and act["family"] == "anim"
        st, fr = call(S, "GET", f"/api/plot/{act['id']}", query={"frame": 2})
        assert st == 200, fr
    cmd(S, "PAUSE", "WINDOWSETTINGS,FRAME,2")
    st, ps = call(S, "GET", "/api/plots")
    act = next(p for p in ps["plots"] if p["id"] == ps["active"])
    assert act["view"]["paused"] is True and act["params"]["current"] == 2
    cmd(S, "WINDOWSETTINGS,DIRECTION,ALL")
    assert next(p for p in call(S, "GET", "/api/plots")[1]["plots"] if p["id"] == act["id"])["view"]["direction"] == "ALL"


# ================================================================== Modules menu
def test_optional_modules_follow_the_command_catalogue(monkeypatch):
    """Modules > NONLINEAR (Option NON, P2) and the Run of the LOADGEN dialogs (Option A, P2) run in worker jobs
    once RUNNONLINEAR / RUNLOADGEN are implemented (both are), and not while their command is a placeholder."""
    assert apimod.run_modules()[-2:] == ("NONLINEAR", "LOADGEN")
    assert not apimod.lookup("RUNNONLINEAR").placeholder and not apimod.lookup("RUNLOADGEN").placeholder

    class Spec:
        name, placeholder, handler = "RUNNONLINEAR", True, None

    real = apimod.lookup
    monkeypatch.setattr(apimod, "lookup", lambda n: Spec() if n.upper() == "RUNNONLINEAR" else real(n))
    assert "NONLINEAR" not in apimod.run_modules() and apimod.run_modules()[-1] == "LOADGEN"


# ================================================================== front end: menus, toolbars, regressions
def _js(name: str) -> str:
    return (STATIC / name).read_text(encoding="utf-8")


def test_menus_and_toolbars_of_requirements_5_3_and_5_9():
    app, plots = _js("app.js"), _js("plots.js")
    labels = set(re.findall(r'label: "([^"]+)"', app))
    menu_5_3 = ["New", "Open", "Save", "Input", "Converters", "SASSI .hou", "ANSYS .cdb", "GT-STRUDL Database", "Output",
                "Export to ANSYS", "Export to STRUDL", "Exit", "Export Image", "Export Table", "Model", "Elements", "Nodes",
                "Cuts", "Spectrum TFU-TFI", "Time History", "Soil Layers", "Soil Properties", "Non Uniform Soil Field",
                "Process Animation Frame List", "Bubble", "Vector", "Contour", "Deformed Shape", "Location", "Extension",
                "NONLINEAR", "ANSYS Eq. Static Load", "ANSYS Dynamic Load", "ANSYS Super Element Utilities", "Write",
                "Check", "Analysis", "Windows Settings", "Colors", "Font", "Shader Options", "Reset Plot", "Check Errors",
                "Command Window", "Command Display", "Command Echo", "Output Confirmation", "Comments",
                "Warnings & Errors", "Toolbars", "Main Toolbar", "Plot Toolbar", "Help", "About", "Verification"]
    assert not [m for m in menu_5_3 if m not in labels]
    for m in ("EQUAKE", "SOIL", "LIQUEF", "SITE", "POINT", "HOUSE", "PINT", "FORCE", "ANALYS", "COMBIN", "MOTION",
              "STRESS", "RELDISP"):
        assert f'"{m}"' in app
    tools_5_9 = ["(ACTM)", "Model Database", "(SAVE)", "(WRITE)", "(ANSYS)", "(CONVERT,SSI)", "(CONVERT,ANSYS)", "STRUDL",
                 "(CAPTUREPLOT)", "(MODELPLOT)", "(NODEPLOT)", "(CUTPLOT)", "SPECPLOT", "THPLOT", "(LAYERPLOT)",
                 "(SOILPROPPLOT)", "(PROCFRAME)", "(BUBBLEPLOT)", "(VECTORPLOT)", "(CONTOURPLOT)", "(DEFORMPLOT)",
                 "Analysis Options", "(AFWRITE)", "(CNGVIEW)", "(RSTVIEW)", "(CNGCENTER)", "(RSTCENTER)", "(WIREFRAME)",
                 "(SHRINK)", "(ELECOLOR,1)", "(ELECOLOR,2)", "(ELECOLOR,3)", "(NODENUM)", "(ELENUM)", "(GROUPNUM)",
                 "(SHOWDOF)", "(SHOWMASS)", "(PAUSE)", "(DEBUG)"]
    titles = re.findall(r'\["\w+", "([^"]+)"', plots)
    assert not [t for t in tools_5_9 if not any(t in x for x in titles)]


def test_front_end_binds_the_command_record_dialogs():
    """Modules > ANSYS Eq. Static Load / ANSYS Dynamic Load open their dialogs (no direct RUNLOADGEN), whose Run
    posts the dialog with "run"; the options form edits %-paths (xrecords / xtables / xindexed), and the 2D-safe
    plot code is unchanged."""
    app, dlg, css = _js("app.js"), _js("dialogs.js"), _js("styles.css")
    assert 'optionsDialog("LOADGEN")' in app and 'optionsDialog("LOADGENDYN")' in app
    assert 'S.command("RUNLOADGEN' not in app
    assert "payload.run = true" in dlg and "S.openJobTab(r.job)" in dlg and "form.run" in dlg
    assert "xrecords: {}, xtables: {}, xindexed: {}" in dlg                  # posted by OptionsForm.payload
    for kind in ("xrecords", "xtables", "xindexed"):
        assert f'kind: "{kind}"' in dlg
    assert "newkey:" in dlg and "delentry:" in dlg and "deleteEntry" in dlg and 'c.w === "index"' in dlg
    assert "choiceValue(c.choices, s.value)" in dlg and "choiceValue(it.choices, s.value)" in dlg
    assert ".opt-full" in css and "g.full" in dlg


def test_front_end_regression_guards():
    dlg, app, plots, css = _js("dialogs.js"), _js("app.js"), _js("plots.js"), _js("styles.css")
    # (1) Line Selection: the target line numbers are read BEFORE the READSPEC request (its 'lines'
    #     event resets the Starting Number while the request is in flight); new lines are checked
    body = dlg[dlg.index("D.lineSelection = async function"):dlg.index("D.onLinesChanged")]
    add = body[body.index("const addLines = async"):body.index("const axisBox")]
    assert add.index("const s0 = Number(start.value)") < add.index("await S.command(line)")
    assert "for (let k = 0; k < n; k++) checked.add(s0 + k)" in add
    # (2) pre-filled from the active plot of the same kind (lines, titles, log axes, minor ticks, ranges)
    assert 'S.get("/api/plots")' in body and "active.kind === kind" in body and "plots.active" in body
    for f in ("src.params.lines", "st.xtitle", "st.ytitle", "st.log_x", "st.log_y", "st.minor_x", "st.xmin", "st.ymax"):
        assert f in body, f
    # (3) dialogs fit narrow windows: responsive grids, wrapping rows, fieldsets that shrink
    assert "max-width: calc(100vw - 16px)" in css and "repeat(auto-fit, minmax(270px, 1fr))" in css
    assert re.search(r"fieldset \{[^}]*min-width: 0", css) and re.search(r"\.row \{[^}]*flex-wrap: wrap", css)
    assert "repeat(auto-fit, minmax(250px, 1fr))" in css and 'width: "640px"}' not in body
    assert 'class: "minmax"' in body and "grid-template-columns: auto minmax(60px, 1fr) auto minmax(60px, 1fr)" in css
    # no native prompt / confirm / alert (they block the page): in-page D.ask / D.confirm
    for src in (dlg, app, plots):
        assert not re.search(r"(?<![\w.])(prompt|confirm|alert)\(", src)
    # (4) Help renders the documents (no raw Markdown text dump)
    assert "/api/help/doc?name=" in dlg and 'getAttribute("data-doc")' in dlg and "h.guide" not in dlg
    # File Editor: its own Run (INP) is not appended to the connected file
    assert "skipEditor: id" in app and "S.connectedEditor.tabId !== opts.skipEditor" in app
    # plots: no Plotly cloud upload button; Pause keeps the frame shown
    assert '"sendChartToCloud"' in plots and "WINDOWSETTINGS,FRAME,${frames[t.anim.k]}" in plots


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_javascript_syntax():
    for f in ("app.js", "dialogs.js", "plots.js", "main.js"):
        r = subprocess.run(["node", "--check", str(STATIC / f)], capture_output=True, text=True, timeout=60)
        assert r.returncode == 0, (f, r.stderr)
