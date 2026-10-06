"""GUI JSON API on an in-process session (ARCHITECTURE section 9; requirements 5.1-5.10).

The session (:class:`sassi.ui.api.GuiSession`) is exercised through ``route`` exactly as the HTTP
handler calls it: command round trip (L17, UI-01), state, model JSON for 3D plotting, plot data,
options GET/POST, file listing and path safety, Check Errors window, settings and the Load Model
database."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from sassi.ui.api import GuiSession

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


def load_model(S, tmp_path):
    lines = [ln for ln in MODEL.format(dir=tmp_path / "box").splitlines() if ln.strip()]
    st, r = call(S, "POST", "/api/command", {"lines": lines})
    assert st == 200 and r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    return r


# ------------------------------------------------------------------ commands (L17)
def test_command_round_trip_returns_messages_and_history(S, tmp_path):
    st, r = call(S, "POST", "/api/command", {"line": "N,1,1.5,2,3"})
    assert st == 200 and r["ok"] and r["results"] == [{"line": "N,1,1.5,2,3", "ok": True}]
    kinds = [m["kind"] for m in r["messages"]]
    assert kinds[0] == "ECHO" and r["messages"][0]["text"] == "N,1,1.5,2,3"
    assert S.interp.model.nodes[1].xyz == (1.5, 2.0, 3.0)
    st, r = call(S, "POST", "/api/command", {"line": "NOSUCHCMD,1"})
    assert st == 200 and not r["ok"]
    assert any("NOSUCHCMD Command not found" in m["text"] for m in r["messages"] if m["kind"] == "ERROR")
    # the replay history holds the accepted keyboard commands only (GUI session = .pre file)
    assert S.interp.history == ["N,1,1.5,2,3"]
    st, state = call(S, "GET", "/api/state")
    assert state["history"] == ["N,1,1.5,2,3"] and state["active"] == 0
    assert state["title"] == "SASSI-EDU User Interface (ACS SASSI V3 methodology)"


def test_command_body_validation(S):
    assert call(S, "POST", "/api/command", {})[0] == 400
    assert call(S, "POST", "/api/command", {"line": "N,1\nN,2"})[0] == 400
    assert call(S, "GET", "/api/command")[0] == 405
    assert call(S, "GET", "/api/nothing")[0] == 404


def test_events_stream_messages_state_and_plots(S, tmp_path):
    load_model(S, tmp_path)
    call(S, "POST", "/api/command", {"line": "MODELPLOT"})
    st, ev = call(S, "GET", "/api/events", query={"since": 0})
    types = [e["type"] for e in ev["events"]]
    assert {"message", "state", "plot"} <= set(types)
    plot = [e for e in ev["events"] if e["type"] == "plot"]
    assert plot[0]["event"] == "open" and plot[0]["data"]["plot"]["kind"] == "MODELPLOT"
    last = ev["last"]
    st, ev2 = call(S, "GET", "/api/events", query={"since": last, "timeout": 0})
    assert ev2["events"] == [] and ev2["last"] == last


def test_inp_progress_events(S, tmp_path):
    pre = tmp_path / "m.pre"
    pre.write_text("\n".join(f"N,{i},{i},0,0" for i in range(1, 40)) + "\n", encoding="utf-8")
    st, r = call(S, "POST", "/api/command", {"line": f"INP,{pre}"})
    assert r["ok"] and len(S.interp.model.nodes) == 39
    st, ev = call(S, "GET", "/api/events", query={"since": 0})
    prog = [e for e in ev["events"] if e["type"] == "progress"]
    assert prog and prog[-1]["line"] == prog[-1]["total"] == 39     # current line / total lines (spec 04 4.2)


# ------------------------------------------------------------------ model JSON
def test_model_json_for_3d_plotting(S, tmp_path):
    load_model(S, tmp_path)
    st, m = call(S, "GET", "/api/model")
    assert st == 200 and m["name"] == "box" and m["counts"]["nodes"] == 8
    assert m["nodes"]["id"] == list(range(1, 9)) and m["nodes"]["xyz"][1] == [2.0, 0.0, 0.0]
    assert m["interaction"] == list(range(1, 9))
    assert m["nodes"]["fix"][4] == [1, 1, 1, 1, 1, 1] and m["nodes"]["fix"][0] == [0] * 6
    g = {x["id"]: x for x in m["groups"]}
    assert g[1]["type_name"] == "SOLID" and g[1]["elements"]["nodes"][0] == [5, 6, 7, 8, 1, 2, 3, 4]
    assert g[2]["type_name"] == "BEAMS" and g[2]["elements"]["prop"] == [1]
    assert m["masses"]["node"] == [3] and m["masses"]["tmass"] == [[10.0, 10.0, 10.0]]
    assert m["bbox"] == [0.0, 2.0, 0.0, 2.0, -2.0, 0.0]
    json.dumps(m)                                        # JSON-serialisable
    assert call(S, "GET", "/api/model", query={"number": 7})[0] == 404


def test_plot_data_and_activate_close_by_command(S, tmp_path):
    load_model(S, tmp_path)
    call(S, "POST", "/api/command", {"lines": ["MODELPLOT", "NODEPLOT", "LAYERPLOT"]})
    st, ps = call(S, "GET", "/api/plots")
    assert [p["kind"] for p in ps["plots"]] == ["MODELPLOT", "NODEPLOT", "LAYERPLOT"] and ps["active"] == 3
    assert not ps["auto_render"]                         # the browser draws the plots
    st, d = call(S, "GET", "/api/plot/1")
    assert d["family"] == "3d" and len(d["scene"]["tri"]) == 12 and len(d["scene"]["edges"]) == 1
    assert d["camera"]["basis"] and len(d["scene"]["interaction"]) == 8
    st, d = call(S, "GET", "/api/plot/3")
    assert d["family"] == "layer" and d["table"]["layers"][0]["vs"] == 800
    call(S, "POST", "/api/command", {"lines": ["ACTIVATEPLOT,1", "WIREFRAME,1", "ELECOLOR,2"]})
    st, d = call(S, "GET", "/api/plot/1")
    assert d["view"]["wireframe"] and d["view"]["color_by"] == 2
    call(S, "POST", "/api/command", {"line": "CLOSEPLOT"})
    assert [p["kind"] for p in call(S, "GET", "/api/plots")[1]["plots"]] == ["NODEPLOT", "LAYERPLOT"]
    assert call(S, "GET", "/api/plot/1")[0] == 404


#: 2D model (HOUSE <dim> 1 style): PLANE elements in the X-Z plane (an embedded block of 4 x 2 elements) and a
#: BEAMS stick on its top centre node
PLANE_2D = """N,1,-4,0,-4
N,5,4,0,-4
FILL,1,5
NGEN,2,5,1,5,1,0,0,2
M,1,3e6,0.25,0.15,0.05,0.05
GROUP,1,PLANE
MACT,1
E,1,1,2,7,6
EGEN,3,1,1
EGEN,1,5,1,4
N,20,0,0,4
N,21,0,0,8
GROUP,2,BEAMS
R,1,1,1,1,1,1,1
E,1,13,20,1
E,2,20,21,1
INT,1,15,1,1"""


def test_2d_plane_model_in_the_3d_model_view(S):
    """A 2D model (PLANE elements in the X-Z plane, every y = 0) draws in the 3D model view: one face per PLANE
    element (all drawn), the beams as lines, a flat bounding box (y range 0), a finite camera fit for the
    default view and for the view along Y; NODEPLOT and CUTPLOT draw it too."""
    import math
    r = call(S, "POST", "/api/command", {"lines": PLANE_2D.splitlines() + ["MODELPLOT"]})[1]
    assert r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    st, d = call(S, "GET", "/api/plot/1")
    sc = d["scene"]
    assert st == 200 and d["family"] == "3d" and d["kind"] == "MODELPLOT"
    assert len(sc["faces"]) == 8 and all(sc["face_boundary"]) and len(sc["tri"]) == 16 and len(sc["edges"]) == 2
    assert sc["bbox"] == [-4.0, 4.0, 0.0, 0.0, -4.0, 8.0] and all(p[1] == 0.0 for p in sc["xyz"])
    assert sorted(set(sc["elem_type"])) == [2, 4] and len(sc["interaction"]) == 15
    for view in (None, "CNGVIEW,-90,0,0,0,0,1"):
        if view:
            assert call(S, "POST", "/api/command", {"line": view})[1]["ok"]
            d = call(S, "GET", "/api/plot/1")[1]
        cam = d["camera"]
        assert cam["half"] > 0 and all(math.isfinite(v) for v in cam["extent"] + cam["center"])
        assert cam["center"] == [0.0, 0.0, 2.0]                     # the box centre
    # face-on (rX = -90: looking along +Y): the X-Z plane fills the window, the 12 m height decides the fit
    assert abs(d["camera"]["half"] - 6.0 * 1.08) < 1e-9
    assert [round(v, 9) for v in d["camera"]["basis"][2]] in ([0.0, -1.0, 0.0], [0.0, 1.0, 0.0])
    call(S, "POST", "/api/command", {"lines": ["NODEPLOT", "CUTADD,1,1,1,1-4", "CUTPLOT,1"]})
    for pid in (2, 3):
        st, d = call(S, "GET", f"/api/plot/{pid}")
        assert st == 200 and d["family"] == "3d" and len(d["scene"]["node_id"]) == 17, (pid, d.get("error"))


def test_commands_without_arguments_open_gui_dialogs(S, tmp_path):
    """D-UI-08: in the GUI a plot command without arguments opens its dialog (dialog event); inside a
    .pre file it fails as in batch mode."""
    st, r = call(S, "POST", "/api/command", {"line": "SPECPLOT"})
    assert r["ok"]
    ev = [e for e in call(S, "GET", "/api/events", query={"since": 0})[1]["events"] if e["type"] == "plot"]
    assert ev[-1]["event"] == "dialog" and ev[-1]["data"]["dialog"] == "Line Selection"
    pre = tmp_path / "d.pre"
    pre.write_text("SPECPLOT\n", encoding="utf-8")
    st, r = call(S, "POST", "/api/command", {"line": f"INP,{pre}"})
    assert any("requires arguments in batch mode" in m["text"] for m in r["messages"] if m["kind"] == "ERROR")


def test_line_plot_from_result_file(S, tmp_path):
    f = tmp_path / "00001TR_X.TFU"
    f.write_text("# test\n0.5 1.0 0.0\n1.0 2.0 0.1\n2.0 1.5 0.2\n", encoding="utf-8")
    st, info = call(S, "GET", "/api/fileinfo", query={"name": f.name})
    assert info["plot"] == "spec" and info["columns"] == 2 and info["points"] == 3
    call(S, "POST", "/api/command", {"lines": [f"READSPEC,{f.name},1,1", "SPECPLOT,1", "AXES,1,1,0,0,1,0"]})
    st, lines = call(S, "GET", "/api/lines")
    assert lines["lines"]["1"]["name"] == "00001TR_X.TFU" and lines["lines"]["1"]["n"] == 3
    st, d = call(S, "GET", "/api/plot/1")
    assert d["family"] == "2d" and d["lines"][0]["y"] == [1.0, 2.0, 1.5] and d["settings"]["log_x"]
    st, r = call(S, "POST", "/api/export_table", {"name": "tab.csv"})
    assert st == 200 and (tmp_path / "tab.csv").read_text().splitlines()[0].startswith("Frequency")


# ------------------------------------------------------------------ options dialogs
def test_options_get_and_post_emit_command_text(S, tmp_path):
    st, d = call(S, "GET", "/api/options/ANALYSIS")
    assert st == 200 and [t["name"] for t in d["form"]["tabs"]][0] == "EQUAKE"
    assert d["values"]["records"]["SITE"]["nft"] == 4096
    st, d = call(S, "GET", "/api/options/SITE")
    assert [t["name"] for t in d["form"]["tabs"]] == ["SITE"]
    st, r = call(S, "POST", "/api/options/ANALYSIS", {"records": {"SITE": {"nl": "10", "delt": "0.01"}},
                                                       "lists": {"TOPL": "1, 2;2"}, "strings": {"THFILE": "acc x.acc"}})
    assert st == 200 and r["ok"]
    # <freq2> stays blank: NFFT/2 (requirements 5.4), not frozen at 2048
    assert r["commands"] == ["SITE,0,1,0,10,0,1,0,1,,1,0,0.01,4096,1", "THFILE,acc x.acc", "TOPL,0", "TOPL,1,2,2"]
    assert S.interp.model.topl == [1, 2, 2] and S.interp.model.options.string("THFILE") == "acc x.acc"
    assert S.interp.history[-4:] == r["commands"]        # replayable (L17)
    # an unchanged commit emits nothing
    st, r = call(S, "POST", "/api/options/ANALYSIS", {"records": {"SITE": {"nl": 10}}})
    assert r["commands"] == []
    # dry run: command text only
    st, r = call(S, "POST", "/api/options/ANALYSIS", {"records": {"POINT": {"rad": 12.5}}, "dry_run": True})
    assert r["commands"] == ["POINT,0,0,12.5"] and S.interp.model.options.record("POINT") is None


def test_options_commit_refused_with_chapter10_text(S):
    """UI-06: invalid values are refused with the Chapter 10 message; nothing is executed."""
    st, r = call(S, "POST", "/api/options/ANALYSIS", {"records": {"SITE": {"nl": 3}}})
    assert st == 422 and r["problems"] == ["SITE: Error 47 : Illegal Number of Layers for Halfspace Simulation  "
                                           "[<nl> = 3 (0 or 4..20)]"]
    assert S.interp.model.options.record("SITE") is None
    st, r = call(S, "POST", "/api/options/ANALYSIS", {"records": {"SITE": {"delt": "abc"}}})
    assert st == 422 and "not a number" in r["problems"][0]


def test_model_write_check_dialogs(S):
    st, r = call(S, "POST", "/api/options/MODEL", {"records": {"MOPT": {"incomp": 0}}})
    assert r["commands"] == ["MOPT,0,0,1,1"] and S.interp.model.mopt.get("incomp") == 0
    st, r = call(S, "POST", "/api/options/WRITE", {"records": {"WRITE": {"mdl": 1, "sim": 1, "sim_location": "pre"}}})
    assert st == 200 and S.interp.write_options["mdl"] and S.interp.write_options["sim_location"] == "pre"
    st, r = call(S, "POST", "/api/options/CHECK", {"records": {"CHECK": {"break_at": 5, "show_warnings": 0}}})
    co = S.interp.session["check_options"]
    assert st == 200 and co.break_at == 5 and not co.show_warnings
    st, d = call(S, "GET", "/api/options/CHECK")
    assert d["values"]["records"]["CHECK"]["break_at"] == 5
    assert call(S, "POST", "/api/options/CHECK", {"records": {"CHECK": {"break_at": -1}}})[0] == 422


def test_point_radius_helper(S, tmp_path):
    """POINT tab "From mesh": RADIUS on the excavation elements (the SOLID below grade, D-PNT-05)."""
    load_model(S, tmp_path)
    st, r = call(S, "POST", "/api/options/POINT/radius")
    assert st == 200 and r["ok"] and r["radius"]["average"] == pytest.approx(0.9 * 2.0)
    assert S.interp.history[-1] == "RADIUS"


# ------------------------------------------------------------------ files and path safety
def test_files_listing_and_path_safety(S, tmp_path):
    load_model(S, tmp_path)
    mdir = tmp_path / "box"
    (mdir / "00012TR_X.TFU").write_text("1 2\n", encoding="utf-8")
    (mdir / "00012TR_X.ACC").write_text("0.01\n0.1\n0.2\n", encoding="utf-8")
    (mdir / "FILE1").write_bytes(b"\x00\x01binary")
    (mdir / "sub").mkdir()
    st, lst = call(S, "GET", "/api/files")
    assert st == 200 and lst["dir"] == str(mdir.resolve())
    kinds = {f["name"]: (f["kind"], f["plot"]) for f in lst["files"]}
    assert kinds["00012TR_X.TFU"] == ("spectrum", "spec") and kinds["00012TR_X.ACC"] == ("history", "th")
    assert kinds["FILE1"][0] == "binary" and [d["name"] for d in lst["dirs"]] == ["sub"]
    st, f = call(S, "GET", "/api/file", query={"name": "00012TR_X.ACC"})
    assert st == 200 and f["text"].startswith("0.01")
    st, info = call(S, "GET", "/api/fileinfo", query={"name": "00012TR_X.ACC"})
    assert info["plot"] == "th" and info["pair"] == 0 and info["points"] == 2
    assert call(S, "GET", "/api/file", query={"name": "FILE1"})[0] in (403, 404)       # binary refused
    for bad in ("../../../etc/passwd", "/etc/passwd", "sub/../../../etc/hosts", "~/.ssh/id_rsa"):
        st, r = call(S, "GET", "/api/file", query={"name": bad})
        assert st in (403, 404), bad
    outside = tmp_path.parent / "outside_ui_test.txt"
    st, r = call(S, "POST", "/api/file", {"name": str(outside), "text": "x"})
    assert st == 403 and not outside.exists()
    assert call(S, "GET", "/api/files", query={"dir": "/"})[0] == 403
    # File > Open of a new name creates the file (spec 04 section 4.5); save writes it
    st, r = call(S, "POST", "/api/file", {"name": "notes.pre"})
    assert st == 200 and r["created"] and (mdir / "notes.pre").exists()
    st, r = call(S, "POST", "/api/file", {"name": "notes.pre", "text": "N,9,0,0,0\n"})
    assert (mdir / "notes.pre").read_text() == "N,9,0,0,0\n"


def test_symlink_out_of_root_is_refused(S, tmp_path):
    target = tmp_path.parent / "ui_secret_target"
    target.mkdir(exist_ok=True)
    (target / "s.txt").write_text("secret", encoding="utf-8")
    link = tmp_path / "link"
    try:
        link.symlink_to(target, target_is_directory=True)
    except OSError:
        pytest.skip("symlinks not available")
    st, r = call(S, "GET", "/api/file", query={"name": "link/s.txt"})
    assert st == 403


# ------------------------------------------------------------------ Check Errors window
def test_check_window_shows_err_file(S, tmp_path):
    load_model(S, tmp_path)
    call(S, "POST", "/api/command", {"line": "CHECK"})
    st, d = call(S, "GET", "/api/check")
    assert d["title"] == "CHECK: Errors and Warning for - box"
    assert d["text"].startswith("CHECK: Errors and Warning for - box")
    assert "Errors and Warnings for" in d["text"] and d["source"].endswith("box.err")


# ------------------------------------------------------------------ settings and database
def test_settings_persist_in_sassiini(S, tmp_path):
    st, r = call(S, "POST", "/api/settings/display", {"display": {"echo": False}})
    assert r["ok"] and not r["settings"]["display"]["echo"]
    st, r = call(S, "POST", "/api/settings/locations", {"locations": {"SITE": "/no/such/SITE.exe", "POINT": "built-in"}})
    assert not r["ok"] and "SITE" in r["errors"][0] and r["settings"]["locations"]["SITE"] == "built-in"
    exe = tmp_path / "SITE.exe"
    exe.write_text("#!/bin/sh\n", encoding="utf-8")
    st, r = call(S, "POST", "/api/settings/locations", {"locations": {"SITE": str(exe)}})
    assert r["ok"] and r["settings"]["locations"]["SITE"] == str(exe)
    st, r = call(S, "POST", "/api/settings/extensions", {"extensions": {"SITE": [".site", "_s.out"]}})
    xml = (tmp_path / "settings" / "SASSIini.xml").read_text()
    assert 'name="SITE"' in xml and ".site" in xml and 'echo="0"' in xml
    S2 = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    assert S2.settings.extensions["SITE"] == [".site", "_s.out"] and not S2.settings.display["echo"]
    assert S2.settings.locations["SITE"] == str(exe)
    assert call(S, "POST", "/api/settings/fonts", {})[0] == 404


def test_load_model_database(S, tmp_path):
    mdir = tmp_path / "proj" / "m1"
    st, db = call(S, "POST", "/api/db", {"action": "add_group", "group": "Trial"})
    assert db["groups"] == [{"name": "Trial", "models": []}]
    st, db = call(S, "POST", "/api/db", {"action": "add_model", "group": "Trial", "name": "m1", "path": str(mdir),
                                         "title": "first"})
    assert mdir.is_dir() and db["groups"][0]["models"][0]["name"] == "m1"
    assert call(S, "POST", "/api/db", {"action": "add_model", "group": "Nope", "name": "x", "path": str(mdir)})[0] == 400
    # Open: MDL (+ RESUME when an .sdb exists); a never-saved model opens without error
    st, r = call(S, "POST", "/api/db", {"action": "open", "name": "m1", "path": str(mdir), "title": "first"})
    assert r["ok"] and S.interp.model.name == "m1" and S.interp.model.title == "first"
    call(S, "POST", "/api/command", {"lines": ["N,1,0,0,0", "SAVE"]})
    assert (mdir / "m1.sdb").exists()
    call(S, "POST", "/api/command", {"line": "N,2,0,0,0"})
    st, r = call(S, "POST", "/api/db", {"action": "open", "name": "m1", "path": str(mdir)})
    assert sorted(S.interp.model.nodes) == [1]                         # RESUME
    (mdir / "other.txt").write_text("keep", encoding="utf-8")
    st, db = call(S, "POST", "/api/db", {"action": "remove_model", "group": "Trial", "name": "m1", "path": str(mdir),
                                         "delete_files": True})
    assert not (mdir / "m1.sdb").exists() and (mdir / "other.txt").exists()      # only the model's files
    assert db["groups"][0]["models"] == []
    xml = (tmp_path / "settings" / "SASSIdb.xml").read_text()
    assert '<Group name="Trial"' in xml


def test_exit_reports_unsaved_models(S, tmp_path):
    load_model(S, tmp_path)
    st, r = call(S, "GET", "/api/exit")
    assert r["unsaved"] == [{"number": 0, "name": "box"}]
    call(S, "POST", "/api/command", {"line": "SAVE"})
    assert call(S, "GET", "/api/exit")[1]["unsaved"] == []
    st, r = call(S, "POST", "/api/exit")
    assert r["ok"] and Path(r["settings"]).name == "SASSIini.xml"


def test_about_and_help(S):
    st, a = call(S, "GET", "/api/about")
    assert a["product"] == "SASSI-EDU" and a["python"] and a["numpy"] and a["scipy"]
    st, h = call(S, "GET", "/api/help")
    names = {c["name"] for c in h["commands"]}
    assert {"SITE", "RUNSITE", "MODELPLOT", "VERIFY", "AFWRITE"} <= names
    ids = [d["id"] for d in h["docs"]]
    assert h["home"] == "docs/index.md" and "docs/user/GUI.md" in ids and "examples/README.md" in ids
    assert h["home_doc"]["id"] == "docs/index.md" and h["home_doc"]["html"].startswith("<h1 ")
    st, d = call(S, "GET", "/api/help/doc", query={"name": "docs/user/GUI.md"})
    assert st == 200 and "sassi-gui" in d["html"] and d["title"].startswith("SASSI-EDU graphical user interface")
    assert any(t["slug"].startswith("3-menus") for t in d["toc"])
    st, d = call(S, "GET", "/api/help/doc", query={"name": "../README.md"})
    assert st == 404
    st, r = call(S, "GET", "/api/help/search", query={"q": "Command Entry"})
    assert st == 200 and r["results"]


# ------------------------------------------------------------------ Check Errors pop-up, shader persistence, quoting
def test_check_and_afwrite_push_the_check_window_event(S, tmp_path):
    """Spec 05a section 3 / requirements 5.5: CHECK (and the CHECK of AFWRITE) with messages pops up
    the Check Errors window unless Options > Check > Suppress Error Window is set."""
    load_model(S, tmp_path)
    seq = S.events.last
    call(S, "POST", "/api/command", {"line": "CHECK"})
    evs = [e for e in S.events.since(seq)[0] if e["type"] == "check"]
    assert len(evs) == 1 and evs[0]["suppress"] is False and evs[0]["summary"]
    call(S, "POST", "/api/options/CHECK", {"records": {"CHECK": {"suppress_window": 1}}})
    seq = S.events.last
    call(S, "POST", "/api/command", {"line": "AFWRITE"})
    evs = [e for e in S.events.since(seq)[0] if e["type"] == "check"]
    assert len(evs) == 1 and evs[0]["suppress"] is True
    seq = S.events.last
    call(S, "POST", "/api/command", {"line": "N,99,1,1,1"})               # no CHECK: no event
    assert not [e for e in S.events.since(seq)[0] if e["type"] == "check"]


def test_shader_options_saved_on_change_and_restored(S, tmp_path):
    from sassi.plotting.state import plot_state
    assert call(S, "POST", "/api/command", {"line": "SHADEROPTIONS,14,0.04,0.12,3"})[1]["ok"]
    xml = (tmp_path / "settings" / "SASSIini.xml").read_text()               # written at OK, not only at Exit
    assert '<Shader points="14.0" linew="0.04" shrink="0.12" scale="3.0"' in xml
    S2 = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    sh = plot_state(S2.interp).shader
    assert (sh.points, sh.linew, sh.shrink, sh.scale) == (14.0, 0.04, 0.12, 3.0)
    # a damaged value never stops the GUI: the default is kept
    (tmp_path / "settings" / "SASSIini.xml").write_text(xml.replace('shrink="0.12"', 'shrink="0.9"').replace(
        'points="14.0"', 'points="x"'), encoding="utf-8")
    sh = plot_state(GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings")).interp).shader
    assert (sh.points, sh.linew, sh.shrink, sh.scale) == (10.0, 0.04, 0.06, 3.0)


def test_load_model_open_quotes_paths_and_titles(S, tmp_path):
    d = tmp_path / "Smith, J"
    d.mkdir()
    lines = S.db_open_lines("m", str(d), "Plant, unit 2")
    assert lines == [f'MDL,m,"{d}"', "TIT,Plant, unit 2"]
    r = S.execute(lines)
    assert r["ok"] and Path(S.interp.model.path) == d and S.interp.model.title == "Plant, unit 2"
    assert not [m for m in r["messages"] if m["kind"] == "WARNING"]
    assert S.interp.history[-2] == lines[0]                                  # replays the same folder (L17)


def test_front_end_quotes_path_tokens_and_handles_new_events():
    """The browser builds path tokens with q() (D-PAR-06), opens job tabs for typed runs and the
    Check Errors window on the 'check' event, and the 3D camera uses the centre and the pan."""
    static = Path(__file__).resolve().parents[2] / "sassi" / "ui" / "static"
    app = (static / "app.js").read_text(encoding="utf-8")
    dlg = (static / "dialogs.js").read_text(encoding="utf-8")
    plots = (static / "plots.js").read_text(encoding="utf-8")
    assert "`INP,${S.q(p)}`" in app and "`INP,${S.q(d.path)}`" in app
    # Export to ANSYS submits a bare ANSYS (<model>.inp in the model folder, no path tokens) and downloads it
    assert 'S.command("ANSYS")' in dlg and "APDL written to " in dlg
    for tok in ("`CD,${q(cwd.value)}`", "`WRITE,${q(file.value)},${q(dir.value)}`",
                "`CAPTUREPLOT,${q(name.value)}`", "`READSPEC,${q(fname)}", "`READTH,${q(fname)}", "`SOILPROPPLOT,${q(cur)}`",
                "`PROCFRAME,${q(lf.value)},${q(dir.value)},${q(desc.value)}", "`${kind},${q(sel.directory)}",
                '["SSI", model.value.trim(), q(file.value), q(pre.value)]'):
        assert tok in dlg, tok
    assert 'case "check":' in app and "D().checkErrors" in app and "!S.jobs[job.id]" in app
    assert "d.camera.center" in plots and "v.px" in plots and "v.py" in plots and 'aspectmode: "manual"' in plots
    assert "INCOH" in dlg and "randphz = 180" in dlg and "request_defaults" in dlg and "wave_defaults_by_wopt" in dlg


def test_loadgen_run_arguments():
    """RUNLOADGEN,[STATIC|DYNAMIC],[model] and the RUN<MODULE> order RUNLOADGEN,<model>,[STATIC|DYNAMIC] (as the
    command reads them) for the worker job of a typed RUNLOADGEN."""
    from sassi.ui.api import loadgen_args
    assert loadgen_args([]) == ("STATIC", None)
    assert loadgen_args(["dynamic"]) == ("DYNAMIC", None)
    assert loadgen_args(["DYNAMIC", "2"]) == ("DYNAMIC", 2)
    assert loadgen_args(["1"]) == ("STATIC", 1)
    assert loadgen_args(["0", "Dynamic"]) == ("DYNAMIC", 0)
    assert loadgen_args(["", "3"]) == ("STATIC", 3)
    for bad, text in ((["SOMETIMES"], "argument 1 must be STATIC or DYNAMIC"), (["0", "X"], "argument 2 must be"),
                      (["STATIC", "x"], "argument 2 <model> = 'x' is not an integer")):
        with pytest.raises(ValueError) as ei:
            loadgen_args(bad)
        assert text in str(ei.value)


def test_options_get_of_the_loadgen_dialogs(S, tmp_path):
    """GET /api/options/LOADGEN | LOADGENDYN: the dialog layout (one page, Run button) and its values."""
    load_model(S, tmp_path)
    for name, analysis in (("LOADGEN", "STATIC"), ("LOADGENDYN", "DYNAMIC")):
        st, r = call(S, "GET", f"/api/options/{name}")
        assert st == 200 and r["form"]["name"] == name and r["form"]["run"]["line"] == f"RUNLOADGEN,{analysis}"
        assert len(r["form"]["tabs"]) == 1 and r["values"]["context"]["analysis"] == analysis
        assert r["values"]["context"]["model"] == "box" and "LGFILE" in r["values"]["xrecords"]
        json.dumps(r)
    st, r = call(S, "GET", "/api/options/NONLINEAR")
    assert st == 200 and r["values"]["xrecords"]["EQL"]["disp"] == 0.8 and "BBC" in r["values"]["xindexed_defaults"]
