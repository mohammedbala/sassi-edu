"""Background jobs of the GUI (requirements UI-04, 5.3 Modules menu, Help > Verification).

RUN<MODULE> runs in a worker process: the listing is streamed into the job and the Command History,
progress reaches the status bar, Cancel terminates the process, and a successful run is appended to
the replay history (L17)."""
from __future__ import annotations

import os
import stat
import sys
import time

import pytest

from sassi.ui.api import GuiSession
from sassi.ui.jobs import parse_verify

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


@pytest.fixture
def S(tmp_path):
    s = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    yield s
    run = s.jobs.running()
    if run is not None:
        s.jobs.cancel(run.id)
        _wait(s, run.id, 20)


def _wait(S, jid, timeout=120.0):
    t0 = time.time()
    while S.jobs.get(jid).active:
        if time.time() - t0 > timeout:
            raise AssertionError(f"job {jid} still running after {timeout} s")
        time.sleep(0.05)
    return S.jobs.get(jid).to_dict()


def _model(S, tmp_path):
    r = S.execute([ln for ln in SURFACE.format(dir=tmp_path / "run").splitlines() if ln.strip()])
    assert r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]


def test_run_site_and_point_in_worker(S, tmp_path):
    _model(S, tmp_path)
    assert S.execute(["AFWRITE"])["ok"]
    st, job = S.route("POST", "/api/run/SITE", {}, {})
    assert st == 200 and job["line"] == "RUNSITE" and job["kind"] == "module"
    d = _wait(S, job["id"])
    assert d["state"] == "done" and d["ok"] and d["progress"] == 1.0
    assert d["listing"] == "run_SITE.out"
    texts = [m["text"] for m in d["messages"]]
    assert any("SITE finished with status OK" in t for t in texts)          # listing streamed (D-RUN-05)
    st, job2 = S.route("POST", "/api/run/POINT", {}, {})
    d2 = _wait(S, job2["id"])
    assert d2["ok"], [m["text"] for m in d2["messages"]]
    for f in ("FILE1", "FILE2", "FILE3", "run_SITE.out", "run_POINT.out"):
        assert (tmp_path / "run" / f).exists(), f
    assert S.interp.history[-2:] == ["RUNSITE", "RUNPOINT"]                  # replayable (L17)
    evs = S.events.since(0)[0]
    hist = [e["text"] for e in evs if e["type"] == "message"]
    assert "RUNSITE" in hist and any("SITE finished with status OK" in t for t in hist)
    progress = [e["job"]["progress_text"] for e in evs if e["type"] == "job"]
    assert any("frequency" in p for p in progress)                          # ctx.progress -> status bar
    # incremental polling of the listing
    st, part = S.route("GET", f"/api/jobs/{job['id']}", {"since": ["5"]}, {})
    assert part["next"] == len(d["messages"]) and len(part["messages"]) == len(d["messages"]) - 5
    st, lst = S.route("GET", "/api/jobs", {}, {})
    assert [j["line"] for j in lst["jobs"]] == ["RUNSITE", "RUNPOINT"]


def test_run_reports_missing_prerequisites(S, tmp_path):
    _model(S, tmp_path)
    st, job = S.route("POST", "/api/run/SITE", {}, {})       # no AFWRITE yet
    d = _wait(S, job["id"])
    assert d["state"] == "failed" and not d["ok"]
    assert any("Run AFWRITE first" in m["text"] for m in d["messages"] if m["kind"] == "ERROR")
    assert "RUNSITE" not in S.interp.history
    assert S.route("POST", "/api/run/NOPE", {}, {})[0] == 404


def _sleeper(S, tmp_path, seconds=60, **kw):
    """A long job without disk output: the worker runs an 'external module' (Python reading its
    program from standard input) that sleeps."""
    pidfile = tmp_path / "child.pid"
    prog = (f"import os, time\nopen({str(pidfile)!r}, 'w').write(str(os.getpid()))\n"
            f"print('sleeping', flush=True)\ntime.sleep({seconds})\n")
    spec = {"cwd": str(tmp_path), "lines": [], "models": {}, "afwrite": {},
            "external": {"exe": sys.executable, "dir": str(tmp_path), "stdin": prog}}
    return S.jobs.start("module", "RUNSLEEP", spec, module="SLEEP", **kw), pidfile


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


@pytest.mark.skipif(os.name != "posix", reason="process-group termination is POSIX")
def test_cancel_terminates_the_worker_and_its_children(S, tmp_path):
    job, pidfile = _sleeper(S, tmp_path)
    t0 = time.time()
    while not (pidfile.exists() and pidfile.read_text()) and time.time() - t0 < 20:
        time.sleep(0.05)
    child = int(pidfile.read_text())
    assert _alive(child)
    assert S.route("POST", f"/api/jobs/{job.id}/cancel", {}, {})[0] == 200
    d = _wait(S, job.id, 20)
    assert d["state"] == "cancelled" and d["ok"] is False
    assert S.jobs.get(job.id).proc.poll() is not None                       # the worker is gone
    t0 = time.time()
    while _alive(child) and time.time() - t0 < 10:                           # and the module it started
        time.sleep(0.05)
    assert not _alive(child)
    evs = S.events.since(0)[0]
    assert any("cancelled by the user" in e.get("text", "") for e in evs if e["type"] == "message")


def test_cancel_while_starting(S, tmp_path):
    job, _ = _sleeper(S, tmp_path)
    S.jobs.cancel(job.id)                    # before or while the worker process starts
    d = _wait(S, job.id, 20)
    assert d["state"] == "cancelled"


def test_one_job_at_a_time(S, tmp_path):
    job, _ = _sleeper(S, tmp_path)
    st2, r = S.route("POST", "/api/run/SITE", {}, {})
    assert st2 == 409 and "still running" in r["error"]
    S.jobs.cancel(job.id)
    _wait(S, job.id, 20)


def test_verify_job_table(S):
    st, job = S.route("POST", "/api/verify", {}, {"select": "VP-30"})
    d = _wait(S, job["id"])
    t = d["table"]
    assert d["ok"] and [p["id"] for p in t["problems"]] == ["VP-30"] and t["problems"][0]["status"] == "PASSED"
    c = t["problems"][0]["checks"][0]
    assert c["passed"] and c["reference"] is not None and c["tolerance"] is not None
    assert t["summary"].startswith("VERIFY VP-30: 1 passed, 0 failed")
    st, lst = S.route("GET", "/api/verify/list", {}, {})
    assert any(p["id"] == "VP-30" for p in lst["listed"])
    assert S.route("POST", "/api/verify", {}, {"select": "VP;rm"})[0] == 400


def test_parse_verify_messages():
    msgs = [{"kind": "ECHO", "text": "VERIFY,P0"},
            {"kind": "INFO", "text": "VP-99: Test problem"},
            {"kind": "INFO", "text": "   pass  peak: computed 1.5, reference 1.5, error 0 (rel), tolerance 0.01"},
            {"kind": "INFO", "text": "   FAIL  shape"},
            {"kind": "INFO", "text": "   note: provisional"},
            {"kind": "INFO", "text": "   VP-99 FAILED (0.25 s)"},
            {"kind": "ERROR", "text": "VERIFY: VP-98 failed to run: RuntimeError('x')"},
            {"kind": "WARNING", "text": "VERIFY P0: 0 passed, 2 failed (work files in /tmp/x)"}]
    t = parse_verify(msgs)
    p = t["problems"][0]
    assert p["id"] == "VP-99" and p["status"] == "FAILED" and p["elapsed"] == 0.25 and p["notes"] == ["provisional"]
    assert p["checks"][0]["computed"] == 1.5 and p["checks"][1] == {"passed": False, "quantity": "shape", "computed": None,
                                                                    "reference": None, "error": None, "kind": "bool",
                                                                    "tolerance": None}
    assert t["problems"][1]["status"] == "ERROR"


@pytest.mark.skipif(sys.platform.startswith("win"), reason="POSIX executable script")
def test_external_module_location(S, tmp_path):
    """Modules > Location: an external executable gets the three-line batch protocol (no shell)."""
    _model(S, tmp_path)
    exe = tmp_path / "fake_site"
    exe.write_text(f"#!{sys.executable}\nimport sys\nlines = sys.stdin.read().split()\nprint('FAKE SITE got', lines)\n",
                   encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    st, r = S.route("POST", "/api/settings/locations", {}, {"locations": {"SITE": str(exe)}})
    assert r["ok"]
    assert S.execute(["AFWRITE"])["ok"]                       # the deck the executable reads
    st, job = S.route("POST", "/api/run/SITE", {}, {})
    d = _wait(S, job["id"])
    assert d["ok"] and d["external"] == str(exe)
    assert any("FAKE SITE got ['run', 'run.sit', 'run_site.out']" in m["text"] for m in d["messages"])
    os.remove(exe)


@pytest.mark.skipif(sys.platform.startswith("win"), reason="POSIX executable script")
def test_external_module_preconditions(S, tmp_path):
    """Spec 04 section 15.3: an external module runs only with the deck written by AFWRITE and no
    CHECK error for the module in the last AFWRITE (the messages of RUN<MODULE>)."""
    marker = tmp_path / "ran.txt"
    exe = tmp_path / "fake_site"
    exe.write_text(f"#!{sys.executable}\nopen({str(marker)!r}, 'w').write('x')\n", encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    S.route("POST", "/api/settings/locations", {}, {"locations": {"SITE": str(exe), "POINT": str(exe)}})
    st, r = S.route("POST", "/api/run/SITE", {}, {})          # no MDL
    assert st == 400 and "Model name/path not defined" in r["error"] and r["reported"]
    _model(S, tmp_path)
    st, r = S.route("POST", "/api/run/SITE", {}, {})          # no AFWRITE yet
    assert st == 400 and "run.sit not found -- Run AFWRITE first" in r["error"]
    texts = [e["text"] for e in S.events.since(0)[0] if e["type"] == "message"]
    assert "RUNSITE" in texts and any("Run AFWRITE first" in t for t in texts)   # in the Command History
    S.execute(["AFWRITE"])
    S.interp.session["afwrite"][0]["blocked"]["POINT"] = "Error 61"   # a module with CHECK errors
    (tmp_path / "run" / "run.poi").write_text("stale deck", encoding="utf-8")
    st, r = S.route("POST", "/api/run/POINT", {}, {})
    assert st == 400 and "POINT has CHECK errors in the last AFWRITE (Error 61)" in r["error"]
    st, r = S.route("POST", "/api/command", {}, {"line": "RUNSITE,7"})  # typed, unknown model
    assert st == 200 and not r["ok"] and "model 7 is not in memory" in r["messages"][-1]["text"]
    assert not marker.exists() and not S.jobs.jobs
    os.remove(exe)


def test_typed_run_uses_a_worker_job_and_defers_the_next_commands(S, tmp_path):
    """UI-04: RUN<MODULE> typed in Command Entry runs in the worker (Cancel, progress); the commands
    after it in the same request wait for the job, in order; the replay history keeps the order."""
    _model(S, tmp_path)
    S.execute(["AFWRITE"])
    st, r = S.route("POST", "/api/command", {}, {"lines": ["runsite", "N,50,9,9,9", "RUNPOINT,0"]})
    assert st == 200 and r["ok"] and r["job"]["line"] == "RUNSITE" and r["deferred"] == ["N,50,9,9,9", "RUNPOINT,0"]
    assert [x["line"] for x in r["results"]] == ["runsite"]
    assert 50 not in S.interp.model.nodes                     # not yet: waits for the run
    _wait(S, r["job"]["id"])
    t0 = time.time()
    while (len(S.jobs.jobs) < 2 or S.jobs.running() is not None) and time.time() - t0 < 120:
        time.sleep(0.05)
    jobs = sorted(S.jobs.jobs.values(), key=lambda j: j.id)
    assert [j.line for j in jobs] == ["RUNSITE", "RUNPOINT,0"] and all(j.ok for j in jobs)
    assert 50 in S.interp.model.nodes
    assert S.interp.history[-3:] == ["runsite", "N,50,9,9,9", "RUNPOINT,0"]   # typed text, typed order (L17)
    assert (tmp_path / "run" / "FILE3").exists()


def test_run_history_position_is_the_start_of_the_run(S, tmp_path):
    """A command typed while a module runs is recorded after the run in the replay history."""
    job, _ = _sleeper(S, tmp_path, seconds=1, history="RUNSLEEP", history_index=len(S.interp.history))
    assert S.execute(["N,1,0,0,0"])["ok"]                      # typed during the run
    _wait(S, job.id, 30)
    assert S.interp.history[-2:] == ["RUNSLEEP", "N,1,0,0,0"]


def test_run_inside_inp_stays_synchronous(S, tmp_path):
    """UI-04: inside a .pre file RUN<MODULE> runs in the session's interpreter (no job)."""
    _model(S, tmp_path)
    S.execute(["AFWRITE"])
    pre = tmp_path / "run.pre"
    pre.write_text("RUNSITE\n", encoding="utf-8")
    r = S.execute([f"INP,{pre}"])
    assert r["ok"] and not S.jobs.jobs and (tmp_path / "run" / "FILE2").exists()


@pytest.mark.skipif(os.name != "posix", reason="process-group termination is POSIX")
def test_jobs_shutdown_terminates_the_worker(S, tmp_path):
    """Server stop (Ctrl-C): JobManager.shutdown leaves no worker or module process behind."""
    job, pidfile = _sleeper(S, tmp_path)
    t0 = time.time()
    while not (pidfile.exists() and pidfile.read_text()) and time.time() - t0 < 20:
        time.sleep(0.05)
    child = int(pidfile.read_text())
    S.jobs.shutdown()
    assert S.jobs.get(job.id).proc.poll() is not None
    t0 = time.time()
    while _alive(child) and time.time() - t0 < 10:
        time.sleep(0.05)
    assert not _alive(child)


# ================================================================== tier-P2 module runs: NONLINEAR and LOADGEN
#: one nonlinear spring (GMR backbone) between interaction node 1 and node 2; AFWRITE writes <model>.hou
NL_SPRING = """MDL,nl,{dir}
N,1,0,0,0
N,2,0,0,1
D,1,1,1,1,ROTX,ROTY,ROTZ
D,2,2,1,1,UY,UZ,ROTX,ROTY,ROTZ
INT,1,1,1,1
SC,1,1000,0,0,0,0,0,0.02
GROUP,1,SPRING
RACT,1
E,1,1,2
GRAVITY,9.81
L,1,10,20,600,300,0.05,0.05
L,2,10,22,1200,600,0.02,0.02
TOPL,1
SITE,0,1,0,10,2,1,0,1,2048,1,0,0.005,4096,1
S,1,1,1,1,1,4
BBCX,1,2,1,0.01,0.2
BBCY,1,2,1,10,29
EQL,1.0,2
AOPT,0,0,0,0,0,1,0,0,0,0,0,0,0,0
AFWRITE
"""


def _nl_model(S, tmp_path):
    import numpy as np
    from sassi.conventions import nodal_result_name
    from sassi.modules import nonlinear as NLM
    r = S.execute([ln for ln in NL_SPRING.format(dir=tmp_path / "nl").splitlines() if ln.strip()])
    assert r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    t = np.arange(400) * 0.01                     # the RELDISP histories of the spring ends
    NLM.write_history(tmp_path / "nl" / nodal_result_name(1, 1, "THD"), np.zeros_like(t), 0.01)
    NLM.write_history(tmp_path / "nl" / nodal_result_name(2, 1, "THD"), 0.03 * np.sin(2 * np.pi * 2.5 * t), 0.01)


def test_modules_nonlinear_streams_its_listing(S, tmp_path):
    """Modules > NONLINEAR runs RUNNONLINEAR in a worker like the other modules: the listing streamed into the
    job (output tab) and the Command History, ctx.progress in the status bar, the listing name, the replay
    history (L17); the dialog values of the NONLINEAR tab are what RUNNONLINEAR writes to the .eql."""
    _nl_model(S, tmp_path)
    st, r = S.route("POST", "/api/options/ANALYSIS", {}, {"xrecords": {"EQL": {"disp": 0.9}}})
    assert st == 200 and r["commands"] == ["EQL,0.9,2,0,0,0"]
    st, job = S.route("POST", "/api/run/NONLINEAR", {}, {})
    assert st == 200 and job["line"] == "RUNNONLINEAR" and job["module"] == "NONLINEAR"
    d = _wait(S, job["id"])
    texts = [m["text"] for m in d["messages"]]
    assert d["ok"] and d["state"] == "done", texts[-5:]
    assert any("NONLINEAR finished with status OK" in t for t in texts)        # the module listing, line by line
    assert any("SPRING0001" in t for t in texts) and d["listing"] == "nl_NONLINEAR.out"
    progress = [e["job"]["progress_text"] for e in S.events.since(0)[0] if e["type"] == "job"]
    assert any(p.startswith("NONLINEAR ") for p in progress), progress   # ctx.progress -> status bar
    hist = [e["text"] for e in S.events.since(0)[0] if e["type"] == "message"]
    assert "RUNNONLINEAR" in hist and any("NONLINEAR finished with status OK" in t for t in hist)
    assert S.interp.history[-1] == "RUNNONLINEAR"
    eql = (tmp_path / "nl" / "nl.eql").read_text(encoding="utf-8").splitlines()
    assert "edf = 0.9" in eql and "nonlinopts = 2" in eql                   # the tab's EQL in the deck


@pytest.fixture(scope="module")
def lg_results(tmp_path_factory):
    """HOUSE, FILE8, MOTION and RELDISP results of the vp_loadgen spring model (as tests/unit/test_loadgen_cmds.py)."""
    import shutil
    from sassi.verify import builders as B
    from sassi.verify.problems import vp_loadgen as V
    wd = tmp_path_factory.mktemp("lgui")
    fs, _ = V.spring_model(wd, B, nfft=512, dt=0.02)
    V.control_motion(wd, n=300, seed=3)
    assert V.run_motion_reldisp(wd, B, fs, [(n, (1, 2, 3)) for n in (5, 12, 14, 16)], range(1, 10)) == (0, 0)
    yield wd
    shutil.rmtree(wd, ignore_errors=True)


def _lg_model(S, tmp_path, lg_results):
    import shutil
    md = tmp_path / "m"
    shutil.copytree(lg_results, md)
    assert S.execute([f"MDL,m,{md}"])["ok"]
    return md


def test_loadgen_dialog_run_streams_its_listing(S, tmp_path, lg_results):
    """ANSYS Eq. Static Load > Run: the dialog's commands, then RUNLOADGEN,STATIC in a worker job (listing
    streamed, progress, listing name, replay history: the commands before the run)."""
    md = _lg_model(S, tmp_path, lg_results)
    body = {"xrecords": {"LOADGEN": {"data": 3, "genmass": 1}, "LGTIME": {"crit": "V", "n": 1}}, "run": True}
    st, r = S.route("POST", "/api/options/LOADGEN", {}, body)
    assert st == 200 and r["ok"] and r["commands"] == ["LOADGEN,3,0,0,0,1,1,RESULTS"], r
    job = r["job"]
    assert job["line"] == "RUNLOADGEN,STATIC" and job["module"] == "LOADGEN"
    d = _wait(S, job["id"])
    texts = [m["text"] for m in d["messages"]]
    assert d["ok"], texts[-5:]
    assert any("LOADGEN finished with status OK" in t for t in texts) and d["listing"] == "m_LOADGEN.out"
    assert (md / "m_LGS.inp").is_file() and (md / "m.masl").is_file()
    progress = [e["job"]["progress_text"] for e in S.events.since(0)[0] if e["type"] == "job"]
    assert any(p.startswith("LOADGEN") for p in progress), progress
    assert S.interp.history[-2:] == ["LOADGEN,3,0,0,0,1,1,RESULTS", "RUNLOADGEN,STATIC"]
    # Ok only (no run): the settings are stored, no job; Run without a change still runs
    n = len(S.jobs.jobs)
    st, r = S.route("POST", "/api/options/LOADGENDYN", {}, {"xrecords": {"LOADGENDYN": {"alpha": 0.4, "beta": 0.002}}})
    assert st == 200 and "job" not in r and len(S.jobs.jobs) == n
    st, r = S.route("POST", "/api/options/LOADGENDYN", {}, {"run": True})
    assert st == 200 and r["commands"] == [] and r["job"]["line"] == "RUNLOADGEN,DYNAMIC"
    d = _wait(S, r["job"]["id"])
    assert d["ok"] and (md / "m_LGD.inp").is_file(), [m["text"] for m in d["messages"]][-5:]


def test_loadgen_run_refusals(S, tmp_path):
    """A Run that cannot start reports it (no MDL: the job fails with the command's message; a running job: 409)."""
    st, r = S.route("POST", "/api/options/LOADGEN", {}, {"run": True})
    assert st == 200 and r["job"]["line"] == "RUNLOADGEN,STATIC"
    d = _wait(S, r["job"]["id"])
    assert not d["ok"] and any("MDL" in m["text"] for m in d["messages"] if m["kind"] == "ERROR")
    assert S.route("POST", "/api/run/LOADGEN", {}, {"args": ["SOMETIMES"]})[0] == 400
    assert S.route("POST", "/api/run/SITE", {}, {"args": ["STATIC"]})[0] == 400
    assert S.route("POST", "/api/run/LOADGEN", {}, {"args": "STATIC"})[0] == 400
    job, _ = _sleeper(S, tmp_path, seconds=30)
    try:
        st, r = S.route("POST", "/api/options/LOADGEN", {}, {"xrecords": {"LGOPT": {"digits": 9}}, "run": True})
        assert st == 200 and not r["ok"] and "still running" in r["run_error"]
        assert r["commands"] == ["LGOPT,9,0,1"]                        # Ok's part is done; the run is refused
    finally:
        S.jobs.cancel(job.id)
        _wait(S, job.id, 20)


def test_typed_runloadgen_uses_a_worker_job(S, tmp_path, lg_results):
    """RUNLOADGEN typed in Command Entry runs like the menu (UI-04), in either argument order of the command."""
    _lg_model(S, tmp_path, lg_results)
    st, r = S.route("POST", "/api/command", {}, {"line": "RUNLOADGEN,0,dynamic"})
    assert st == 200 and r["job"]["line"] == "RUNLOADGEN,DYNAMIC,0"
    d = _wait(S, r["job"]["id"])
    assert d["ok"], [m["text"] for m in d["messages"]][-5:]
    assert S.interp.history[-1] == "RUNLOADGEN,0,dynamic"                    # the typed text replays (L17)
    st, r = S.route("POST", "/api/command", {}, {"line": "RUNLOADGEN,SOMETIMES"})
    assert st == 200 and not r["ok"] and "job" not in r
    assert any("argument 1 must be STATIC or DYNAMIC" in m["text"] for m in S.events.since(0)[0]
               if m["type"] == "message" and m["kind"] == "ERROR")
    # a module error: the job fails, the listing (with its error) is streamed and the run is not replayed
    st, job = S.route("POST", "/api/run/LOADGEN", {}, {"args": ["static"], "model": 0})
    assert st == 200 and job["line"] == "RUNLOADGEN,STATIC,0"
    d = _wait(S, job["id"])
    assert d["state"] == "failed" and any("m.masl not found" in m["text"] for m in d["messages"])
    assert any("status FAILED" in m["text"] and "m_LOADGEN.out" in m["text"] for m in d["messages"] if m["kind"] == "ERROR")
    assert S.interp.history[-1] == "RUNLOADGEN,0,dynamic"
    assert S.execute(["LOADGEN,2,0,0,0,1,1"])["ok"]                       # Generate Mass Data
    st, job = S.route("POST", "/api/run/LOADGEN", {}, {"args": ["static"]})
    d = _wait(S, job["id"])
    assert d["ok"], [m["text"] for m in d["messages"] if m["kind"] in ("ERROR", "WARNING")]
    assert S.interp.history[-1] == "RUNLOADGEN,STATIC"
