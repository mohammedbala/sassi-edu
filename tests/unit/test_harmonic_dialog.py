"""Deformed shapes from the analysis results (Load Frame Data > animate the analysis results; sassi/ui/harmonic.py).

No example writes animation frames, so after an analysis the Deformed Shape dialog had nothing to list.  It now
offers the steady-state motion at a computed frequency from FILE8: GET /api/harmonic lists the FILE8-type files
of the active model with the frequency of the largest deformation, POST /api/harmonic/plan gives the HARMFRAME
line, POST /api/harmonic/show the PROCFRAME + DEFORMPLOT lines, all run as command text (rule L17)."""
from __future__ import annotations

from pathlib import Path

import pytest

from sassi.ui import harmonic

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"

#: a surface mat on a two-layer site, solved at two frequencies (the model of the lesson-animate tests)
MODEL = """MDL,run,run
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
RUNANALYS""".splitlines()


@pytest.fixture(autouse=True)
def _settings(tmp_path, monkeypatch):
    """PROCFRAME writes SASSIani.xml in the settings folder: keep it in the test folder."""
    monkeypatch.setenv("SASSI_EDU_SETTINGS_DIR", str(tmp_path / "settings"))


@pytest.fixture
def session(tmp_path):
    """A session of the browser version: module runs are queued jobs that run_pending runs in-process."""
    from sassi.ui.api import GuiSession
    return GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"), web=True)


def _run(S, lines):
    st, r = S.route("POST", "/api/command", {}, {"lines": lines})
    assert st == 200 and r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    while S.jobs.run_pending():                        # RUN<MODULE> jobs, then the commands deferred after them
        pass
    return r


def _solve(S):
    _run(S, MODEL)
    assert all(j.state == "done" for j in S.jobs.jobs.values()) and (Path(S.interp.model.path) / "FILE8").is_file()


def test_sources_plan_show_through_the_api(session, tmp_path):
    S = session
    st, r = S.route("GET", "/api/harmonic", {}, {})
    assert st == 200 and r["sources"] == []                   # nothing solved yet
    _solve(S)
    st, r = S.route("GET", "/api/harmonic", {}, {})
    assert st == 200 and r["model"] == "run" and [s["file"] for s in r["sources"]] == ["FILE8"], r
    src = r["sources"][0]
    assert len(src["freqs"]) == 2 and src["peak"] in src["freqs"] and src["seismic"] and src["direction"] == "X"
    assert src["peak_dir"] == "X" and src["peak_node"] in (1, 2, 3, 4)
    # the HARMFRAME line of a choice, in a frame folder of its own; relative: ",,0"
    st, p = S.route("POST", "/api/harmonic/plan", {}, {"file": "FILE8", "freq": src["peak"]})
    assert st == 200 and p["folder"] == harmonic.folder_name(src["peak"], False)
    assert p["lines"] == [f"HARMFRAME,FILE8,{src['peak']:.6g},{p['folder']}"]
    st, pr = S.route("POST", "/api/harmonic/plan", {}, {"file": "FILE8", "freq": src["peak"], "relative": True})
    assert pr["folder"].endswith("R") and pr["lines"][0].endswith(",,0")
    _run(S, p["lines"])
    st, sh = S.route("POST", "/api/harmonic/show", {}, {"folder": p["folder"], "title": "run at 0.39 Hz"})
    assert st == 200 and sh["lines"][0].startswith(f"PROCFRAME,{p['folder']},{p['folder']}_ani")
    assert any(ln.startswith(f"DEFORMPLOT,{p['folder']}_ani,1,24,1,") for ln in sh["lines"])
    _run(S, sh["lines"])
    pid = S.plots.active
    assert S.plots.plots[pid].kind == "DEFORMPLOT"
    st, d = S.route("GET", f"/api/plot/{pid}", {"frame": ["2"]}, {})
    assert st == 200 and "frame 2" in d["label"]
    hist = S.route("GET", "/api/state", {}, {})[1]["history"]
    assert any(h.startswith("HARMFRAME,FILE8,") for h in hist) and any(h.startswith("DEFORMPLOT,") for h in hist)


def test_bad_requests(session):
    S = session
    st, r = S.route("POST", "/api/harmonic/plan", {}, {"file": "FILE8", "freq": 1})
    assert st == 404                                          # no model folder yet
    _run(S, ["MDL,m,m"])
    for body, status in (({"file": "FILE9", "freq": 1}, 404), ({"file": "../FILE8", "freq": 1}, 404)):
        st, r = S.route("POST", "/api/harmonic/plan", {}, body)
        assert st == status, (body, r)
    (Path(S.interp.model.path) / "FILE8").write_bytes(b"not a container")
    st, r = S.route("GET", "/api/harmonic", {}, {})
    assert st == 200 and r["sources"] == []                   # not a FILE8-type file: not listed
    for f in (0, -1, "x"):
        st, r = S.route("POST", "/api/harmonic/plan", {}, {"file": "FILE8", "freq": f})
        assert st == 400, f
    for folder in ("../x", "OTHER", "HARM_1p5/../..", ""):
        st, r = S.route("POST", "/api/harmonic/show", {}, {"folder": folder})
        assert st == 400, folder
    st, r = S.route("POST", "/api/harmonic/show", {}, {"folder": "HARM_1p5"})
    assert st == 404                                          # HARMFRAME has not written it


def test_folder_names():
    assert harmonic.folder_name(3.49121, False) == "HARM_3p491"
    assert harmonic.folder_name(12.0, True) == "HARM_12R"
    assert harmonic.FILE8_NAME.fullmatch("FILE8X") and harmonic.FILE8_NAME.fullmatch("file81")
    assert not harmonic.FILE8_NAME.fullmatch("FILE1") and not harmonic.FILE8_NAME.fullmatch("../FILE8")


def test_front_end_wiring():
    js = (STATIC / "dialogs.js").read_text(encoding="utf-8")
    assert 'S.get("/api/harmonic")' in js and 'S.post("/api/harmonic/plan"' in js and 'S.post("/api/harmonic/show"' in js
    assert "function fromResults(hs, getApi)" in js and "largest deformation" in js
    assert 'if (kind === "DEFORMPLOT") body.appendChild(fromResults(hs, () => api));' in js
