"""The browser version of the GUI (GitHub Pages, docs/internal/web_version.md).

Python side, in CPython: the inline jobs (:class:`sassi.ui.jobs.InlineJobManager`), the bridge of the Web
Worker (:mod:`sassi.web.bridge`) with a fake ``push`` -- a typed RUN<MODULE> queued by its request, run by
``run_pending``, its listing and the events pushed, the deferred commands and the replay history --,
``/api/events`` that never waits, the relative ``docs/`` links; the local session is unchanged.  The site:
``web/build.py`` (what the bundle holds and must never hold, every file the page references), the front-end
wiring of the web mode, and -- when ``web/node_modules/pyodide`` exists (``npm install`` in ``web/``) -- the
built bundle driven in Pyodide under Node (``web/test_pyodide.mjs --quick``)."""
from __future__ import annotations

import importlib.util
import json
import re
import shutil
import subprocess
import time
import zipfile
from pathlib import Path

import pytest

from sassi.ui.api import GuiSession
from sassi.ui.jobs import InlineJobManager, JobBusy, JobManager
from sassi.web import bridge as bridge_mod
from sassi.web.bridge import Bridge

ROOT = Path(__file__).resolve().parents[2]
WEB = ROOT / "web"
STATIC = ROOT / "sassi" / "ui" / "static"

SURFACE = """MDL,run,{dir}
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
"""


class Page:
    """The page side of the protocol (web/boot.js + static/app.js): requests as text, pushes collected, the
    job view de-duplicated by message index as app.js applyJob does."""

    def __init__(self, tmp_path):
        self.pushes = []
        self.bridge = Bridge(lambda text: self.pushes.append(json.loads(text)), cwd=str(tmp_path / "work"),
                             settings_dir=str(tmp_path / "settings"))
        self.shown = {}

    def call(self, method, path, body=None, query=""):
        text = self.bridge.handle(method, path, query, "" if body is None else json.dumps(body))
        d = json.loads(text)
        return d["status"], d["payload"]

    def apply_job(self, d):
        j = self.shown.setdefault(d["id"], {"next": 0, "lines": []})
        first = d["next"] - len(d["messages"])
        j["lines"] += [m["text"] for m in d["messages"][max(0, j["next"] - first):]]
        j["next"] = max(j["next"], d["next"])


@pytest.fixture
def page(tmp_path, monkeypatch):
    monkeypatch.setattr(bridge_mod, "PUSH_INTERVAL", 0.0)          # push at every event: deterministic
    return Page(tmp_path)


def _model(page, tmp_path):
    lines = [ln for ln in SURFACE.format(dir=tmp_path / "run").splitlines() if ln.strip()] + ["AFWRITE"]
    st, r = page.call("POST", "/api/command", {"lines": lines})
    assert st == 200 and r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]


# ====================================================================== the session in web mode
def test_web_session_uses_inline_jobs_and_never_long_polls(tmp_path):
    web = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "s"), web=True)
    local = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "s2"))
    assert isinstance(web.jobs, InlineJobManager) and isinstance(local.jobs, JobManager)
    assert web.state()["web"] is True and local.state()["web"] is False
    t0 = time.time()
    st, r = web.route("GET", "/api/events", {"since": ["1000"], "timeout": ["20"]}, {})
    assert st == 200 and r["events"] == [] and time.time() - t0 < 1.0           # the timeout is ignored
    t0 = time.time()
    local.route("GET", "/api/events", {"since": ["1000"], "timeout": ["0.3"]}, {})
    assert time.time() - t0 >= 0.25                                              # the local server still waits
    a = web.about()
    assert a["web"] and "browser" in a["build"] and not local.about()["web"]


def test_inline_job_manager(tmp_path):
    seen = {"msg": [], "progress": 0, "finished": []}
    m = InlineJobManager(lambda j, msg: seen["msg"].append(msg["text"]), lambda j: None,
                         lambda j: seen["finished"].append((j.id, j.state)))
    spec = {"cwd": str(tmp_path), "lines": ["VERIFY,LIST"], "models": {}, "afwrite": {}}
    job = m.start("verify", "VERIFY,LIST", spec)
    assert job.state == "starting" and m.pending and m.running() is job      # queued, not run
    with pytest.raises(JobBusy):
        m.start("verify", "VERIFY,LIST", spec)
    assert m.run_pending() is False and not m.pending
    assert job.state == "done" and job.ok and seen["finished"] == [(1, "done")]
    assert any(t.startswith("VP-") for t in seen["msg"])
    # cancel before the start: dropped from the queue; a running inline job is not interruptible
    j2 = m.start("verify", "VERIFY,LIST", spec)
    assert m.cancel(j2.id).state == "cancelled" and not m.pending and m.run_pending() is False
    # Modules > Location: an executable cannot run in a browser tab
    j3 = m.start("module", "RUNSITE", {"external": {"exe": "/bin/site", "dir": str(tmp_path), "stdin": ""}}, "SITE")
    m.run_pending()
    assert j3.state == "failed" and "cannot run in the browser" in j3.messages[-1]["text"]


def test_run_spec_restores_the_module_runners(tmp_path):
    """run_spec wraps run_module for the progress only while it runs (an inline job shares the process)."""
    import sassi.modules.base as base
    from sassi.ui.worker import run_spec
    original = base.run_module
    sent = []
    assert run_spec({"cwd": str(tmp_path), "lines": ["VERIFY,LIST"]}, sent.append)
    assert base.run_module is original and sent[-1] == {"t": "done", "ok": True}


# ====================================================================== the bridge
def test_bridge_typed_runs_inline_with_pushes_and_deferred_commands(page, tmp_path):
    st, s = page.call("GET", "/api/state")
    assert st == 200 and s["web"] and s["cwd"] == str((tmp_path / "work").resolve())
    _model(page, tmp_path)
    st, ev = page.call("GET", "/api/events", query="since=0&timeout=20")       # what the page polls
    page.pushes.clear()
    st, r = page.call("POST", "/api/command", {"lines": ["runsite", "N,50,9,9,9", "RUNPOINT,0"]})
    assert st == 200 and r["ok"] and r["job"]["state"] == "starting" and r["deferred"] == ["N,50,9,9,9", "RUNPOINT,0"]
    assert page.bridge.pending and 50 not in page.bridge.session.interp.model.nodes   # answered before the run
    assert page.bridge.run_pending() is True             # the deferred RUNPOINT,0 queued a second job
    assert page.bridge.run_pending() is False and not page.bridge.pending
    S = page.bridge.session
    jobs = sorted(S.jobs.jobs.values(), key=lambda j: j.id)
    assert [j.line for j in jobs] == ["RUNSITE", "RUNPOINT,0"] and all(j.state == "done" for j in jobs)
    assert 50 in S.interp.model.nodes
    assert S.interp.history[-3:] == ["runsite", "N,50,9,9,9", "RUNPOINT,0"]      # typed text, typed order (L17)
    assert (tmp_path / "run" / "FILE3").exists()
    # pushes: the job listing while it ran (de-duplicated with the later answer), the progress, the events
    for p in page.pushes:
        if "job" in p:
            page.apply_job(p["job"])
    for j in jobs:
        st, d = page.call("GET", f"/api/jobs/{j.id}", query="since=0")
        page.apply_job(d)
        assert page.shown[j.id]["lines"] == [m["text"] for m in d["messages"]]   # every line once, in order
    assert any("SITE finished with status OK" in t for t in page.shown[jobs[0].id]["lines"])
    seqs = [e["seq"] for p in page.pushes for e in p.get("events", [])]
    assert seqs == sorted(set(seqs)) and seqs[0] == ev["last"] + 1                 # no gap, no repeat
    evs = [e for p in page.pushes for e in p.get("events", [])]
    assert any(e["type"] == "job" and "frequency" in e["job"]["progress_text"] for e in evs)   # status bar
    assert any(e["type"] == "message" and "SITE finished with status OK" in e["text"] for e in evs)
    st, ev2 = page.call("GET", "/api/events", query=f"since={seqs[-1]}")
    assert ev2["events"] == []                                                       # nothing left behind


def test_bridge_requests_and_cancel(page, tmp_path):
    assert page.call("GET", "/static/app.js")[0] == 404
    assert page.call("POST", "/api/command", None)[0] == 400                        # no line
    assert json.loads(page.bridge.handle("POST", "/api/command", "", "{nope"))["status"] == 400   # not JSON
    assert json.loads(page.bridge.handle("POST", "/api/command", "", "[1]"))["status"] == 400     # not an object
    _model(page, tmp_path)
    st, r = page.call("POST", "/api/run/SITE", {})
    assert st == 200 and r["state"] == "starting"
    st, c = page.call("POST", f"/api/jobs/{r['id']}/cancel", {})
    assert st == 200 and c["state"] == "cancelled" and not page.bridge.pending
    st, ev = page.call("GET", "/api/events", query="since=0")
    assert any("cancelled by the user (before it started)" in e.get("text", "") for e in ev["events"])
    assert "RUNSITE" not in page.bridge.session.interp.history


def test_bridge_files_help_and_lessons(page, tmp_path):
    work = (tmp_path / "work").resolve()
    st, r = page.call("POST", "/api/file", {"name": f"{work}/up.pre", "text": "N,1,0,0,0\n"})     # Upload
    assert st == 200 and r["ok"]
    st, d = page.call("GET", "/api/file", query=f"name={work}/up.pre")                           # Download
    assert st == 200 and d["text"] == "N,1,0,0,0\n"
    assert page.call("POST", "/api/file", {"name": "/etc/x.pre", "text": "x"})[0] == 403
    st, d = page.call("GET", "/api/help/doc", query="name=docs/verification/VERIFICATION_MANUAL.md")
    assert re.search(r'<img src="docs/verification/figures/[^"]+\.png"', d["html"])            # relative URL
    st, o = page.call("POST", "/api/lessons/01-why-ssi/open", {})
    assert st == 200 and o["workspace"].startswith(str(work))                                   # in the workspace


# ====================================================================== the site
def _load_build():
    spec = importlib.util.spec_from_file_location("sassi_web_build", WEB / "build.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def site(tmp_path_factory):
    pytest.importorskip("plotly")
    out = tmp_path_factory.mktemp("web") / "site"
    summary = _load_build().build(out)
    yield out, summary
    shutil.rmtree(out, ignore_errors=True)


FORBIDDEN = ("reference/", "docs/spec/", "docs/internal/", "tests/", "examples/sassi-course/", ".venv/")


def test_bundle_contents(site):
    out, summary = site
    names = zipfile.ZipFile(out / "web" / "sassi-edu.zip").namelist()
    for bad in FORBIDDEN:
        assert not [n for n in names if n.startswith(bad)], bad
        assert not [f for f in summary["files"] if f.startswith(bad)], bad
    assert not [n for n in names if "__pycache__" in n or n.endswith(".pyc") or n.startswith("sassi/ui/static/")]
    assert not [n for n in names if n.startswith(("docs/ARCHITECTURE", "pyproject", "web/", "sassi_edu.egg"))]
    for need in ("README.md", "sassi/__init__.py", "sassi/ui/api.py", "sassi/web/bridge.py", "sassi/ui/worker.py",
                 "sassi/ui/command_args.json", "sassi/ui/lessons/01_why_ssi.md", "examples/README.md",
                 "examples/ex01_surface_stick.pre", "examples/data/ricker_5hz.th", "docs/index.md",
                 "docs/user/GUI.md", "docs/theory/THEORY_MANUAL.md", "docs/reference/COMMAND_REFERENCE.md",
                 "docs/verification/VERIFICATION_MANUAL.md", "docs/verification/figures/vp_margins.png"):
        assert need in names, need
    assert summary["zip_bytes"] < 5_000_000                                          # a few MB
    assert _load_build().build(out)["build"] == summary["build"]                     # reproducible


def test_site_page_references_exist(site):
    out, summary = site
    html = (out / "index.html").read_text(encoding="utf-8")
    assert "{{TOKEN}}" not in html and "sassi-token" not in html and "Content-Security-Policy" in html
    refs = re.findall(r'(?:src|href)="([^"]+)"', html)
    assert refs and not [r for r in refs if r.startswith(("/", "http"))]                # relative URLs only
    for r in refs:
        assert (out / r.split("?")[0]).is_file(), r
        assert r.endswith(f"?v={summary['build']}"), r                                  # cache busting
    assert html.index("web/boot.js") < html.index("static/plotly.min.js") < html.index("static/app.js")
    boot = (out / "web" / "boot.js").read_text(encoding="utf-8")
    worker = (out / "web" / "worker.js").read_text(encoding="utf-8")
    assert "__SASSI_BUILD__" not in boot + worker and summary["build"] in boot and summary["build"] in worker
    assert 'new Worker(`web/worker.js?v=${BUILD}`, {type: "module"})' in boot
    assert "sassi-edu.zip?v=${BUILD}" in worker and "cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}" in worker
    assert 'PYODIDE_VERSION = "314.0.7"' in worker and 'PYODIDE_VERSION = "314.0.7"' in boot
    css = (out / "static" / "katex" / "katex.min.css").read_text(encoding="utf-8")
    for font in set(re.findall(r"url\((fonts/[^)]+)\)", css)):
        assert (out / "static" / "katex" / font).is_file(), font
    assert (out / ".nojekyll").is_file() and (out / "docs/verification/figures/vp_margins.png").is_file()


def test_site_ships_the_example_pictures(site):
    """The gallery pictures (sassi/ui/static/examples/*.png, a subfolder of the front end) are files of the site
    at the URL the gallery gives them (learn.THUMB_URL, relative to the page), byte for byte."""
    from sassi.ui import learn
    out, summary = site
    pics = sorted((STATIC / "examples").glob("*.png"))
    assert len(pics) == len(list((ROOT / "examples").glob("*.pre")))
    for p in pics:
        name = f"{learn.THUMB_URL}/{p.name}"
        assert name in summary["files"] and (out / name).read_bytes() == p.read_bytes(), name
    assert not [f for f in summary["files"] if f.startswith("static/examples/") and not f.endswith(".png")]


def test_site_ships_the_explainer_videos(site):
    """Every file of sassi/ui/static/videos (pages, scripts, data, the list index.html and the recorded
    narration audio/NN/*.mp3) is a file of the site, byte for byte; only the GUI's own static/index.html is
    replaced by the built page."""
    out, summary = site
    videos = sorted(p for p in (STATIC / "videos").rglob("*") if p.is_file() and p.name != ".DS_Store")
    assert videos and any(p.suffix == ".mp3" for p in videos)
    for p in videos:
        name = "static/" + p.relative_to(STATIC).as_posix()
        assert name in summary["files"] and (out / name).read_bytes() == p.read_bytes(), name
    assert "static/videos/index.html" in summary["files"] and "static/index.html" not in summary["files"]


def test_site_service_worker_keeps_the_app_on_the_computer(site):
    """sw.js at the site root (its scope is the site folder), with this build and the Pyodide version of
    worker.js: the Pyodide cache is named by the version only (kept across deployments), the site cache by the
    build (replaced by the next build); boot.js registers it before starting Python."""
    out, summary = site
    sw = (out / "sw.js").read_text(encoding="utf-8")
    assert "__SASSI_BUILD__" not in sw and "__PYODIDE_VERSION__" not in sw
    assert f'const BUILD = "{summary["build"]}";' in sw
    version = re.search(r'const PYODIDE_VERSION = "([0-9.]+)"', (WEB / "worker.js").read_text(encoding="utf-8")).group(1)
    assert f'const PYODIDE_VERSION = "{version}";' in sw
    assert "sassi-edu-pyodide-${PYODIDE_VERSION}" in sw and "sassi-edu-site-${BUILD}" in sw
    boot = (out / "web" / "boot.js").read_text(encoding="utf-8")
    assert 'navigator.serviceWorker.register("sw.js", {scope: "./"})' in boot
    assert f'const PYODIDE_VERSION = "{version}";' in boot          # boot.js reads the same Pyodide cache name


def test_build_refuses_to_empty_other_folders(tmp_path):
    b = _load_build()
    (tmp_path / "mine.txt").write_text("keep")
    with pytest.raises(b.BuildError):
        b.build(tmp_path)
    with pytest.raises(b.BuildError):
        b.build(ROOT)
    assert (tmp_path / "mine.txt").read_text() == "keep"


def test_web_front_end_wiring():
    app = (STATIC / "app.js").read_text(encoding="utf-8")
    dialogs = (STATIC / "dialogs.js").read_text(encoding="utf-8")
    # the transport replaces fetch; the poll never long-polls in web mode; pushes are handled by number
    assert "S.transport = window.SASSI_TRANSPORT || null;" in app and "if (S.transport) return webApi(" in app
    assert "timeout=${S.web ? 0 : 20}" in app and "if (S.web) await webWait(250);" in app
    assert "S.push = function (p)" in app and "if (ev.seq <= S.lastEvent) continue;" in app
    assert "const msgs = d.messages.slice(Math.max(0, j.next - first));" in app
    # what cannot work in a browser is hidden; Upload / Download only there
    assert '{label: "Exit", action: () => D().exit(), web: false}' in app
    assert '{label: "Location", action: () => D().moduleLocation(), web: false}' in app
    assert 'label: "Upload to Workspace...", action: () => D().uploadFiles(), webOnly: true' in app
    assert 'label: "Download...", action: () => D().downloadFile(), webOnly: true' in app
    assert "(S.web && it.web === false) || (!S.web && it.webOnly)" in app
    assert "D.uploadFiles = function" in dialogs and "D.downloadFile = function" in dialogs
    assert "not qualified for design or licensing work" in dialogs


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_web_js_syntax():
    for f in (WEB / "boot.js", WEB / "worker.js", WEB / "sw.js", WEB / "test_pyodide.mjs", STATIC / "app.js", STATIC / "dialogs.js"):
        r = subprocess.run(["node", "--check", str(f)], capture_output=True, text=True, timeout=60)
        assert r.returncode == 0, (f.name, r.stderr)


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("node") is None or not (WEB / "node_modules" / "pyodide").is_dir(),
                    reason="Pyodide test: run npm install in web/ first")
def test_bundle_in_pyodide(site):
    """The built bundle in Pyodide (Node), driven like the page: lesson 01 with its module runs, plots,
    Help, upload / download (web/test_pyodide.mjs --quick)."""
    out, _ = site
    r = subprocess.run(["node", str(WEB / "test_pyodide.mjs"), "--quick", "--site", str(out)], cwd=str(ROOT),
                       capture_output=True, text=True, timeout=900)
    line = next((ln for ln in r.stdout.splitlines() if ln.startswith("RESULT ")), None)
    assert line is not None, (r.stdout[-3000:], r.stderr[-3000:])
    res = json.loads(line[len("RESULT "):])
    assert r.returncode == 0 and res["ok"], res["failures"]
    assert res["plots"] > 0 and [j["line"] for j in res["jobs"]][:5] == ["RUNSITE", "RUNPOINT", "RUNHOUSE",
                                                                         "RUNANALYS", "RUNMOTION"]
