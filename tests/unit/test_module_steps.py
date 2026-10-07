"""The computation steps the modules report while they run (ctx.announce; requirements 7.20, D-W6-13).

Each module reports the step it enters -- ``<MODULE>.<step>`` with its sizes (the matrices being formed, inverted
or factorised and their dimensions) -- through ``run_module(step=...)``.  The GUI's activity panel shows the
step's equation (sassi/ui/static/activity.js STEPS, the Theory Manual's notation, with a link to its section)
and the numbers.  Checked here: the keys and sizes of a real run, the GUI session's events and job state, the
worker wrapper, and that the panel knows every key the modules report."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from sassi.prep import Interpreter, Kind

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"

#: a surface mat on a two-layer site, solved at two frequencies (the model of the deformed-shape tests)
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
    monkeypatch.setenv("SASSI_EDU_SETTINGS_DIR", str(tmp_path / "settings"))


def test_a_run_reports_its_steps_with_sizes(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    steps = []
    ui.module_step = lambda module, key, data: steps.append((module, key, dict(data)))
    for ln in MODEL:
        assert ui.execute(ln), (ln, ui.sink.texts(Kind.ERROR)[-1:])
    keys = [k for _, k, _ in steps]
    # every module's steps, in order; ANALYS four per frequency (two frequencies)
    assert keys[:4] == ["SITE.modes", "SITE.modes", "SITE.field", "SITE.field"]
    assert keys.count("POINT.green") == 2 and keys.index("HOUSE.elements") < keys.index("HOUSE.assembled")
    analys = [k for k in keys if k.startswith("ANALYS.")]
    assert analys == ["ANALYS.impedance", "ANALYS.dynamic", "ANALYS.factor", "ANALYS.solve"] * 2
    assert all(m == k.split(".")[0] for m, k, _ in steps)
    data = {k: d for _, k, d in steps}
    # the sizes are the model's: 4 interaction nodes -> 12 interaction DOFs, the impedance 12 x 12
    h = data["HOUSE.assembled"]
    assert h["nint"] == 4 and h["neq"] > 0 and h["nnz"] > 0
    assert data["ANALYS.impedance"]["n"] == 12 and data["ANALYS.factor"]["nf"] == 12
    assert data["ANALYS.factor"]["nn"] + data["ANALYS.factor"]["nf"] == h["neq"] == data["ANALYS.dynamic"]["neq"]
    assert data["ANALYS.solve"]["seismic"] == 1 and data["ANALYS.solve"]["q"] == 2 and data["ANALYS.solve"]["nF"] == 2
    assert data["SITE.modes"]["n"] > 10 and data["POINT.green"]["nR"] > 0


def test_a_failing_step_callback_never_fails_the_run(tmp_path):
    from sassi.modules.base import ModuleContext, Listing
    ctx = ModuleContext(module="SITE", model="m", workdir=tmp_path, deck_path=None,
                        listing=Listing(tmp_path / "x.out"),
                        step=lambda key, data: (_ for _ in ()).throw(RuntimeError("boom")))
    ctx.announce("SITE.modes", q=1)                         # swallowed: a display aid


def test_gui_session_events_and_job_state(tmp_path):
    from sassi.ui.api import GuiSession
    S = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"), web=True)
    st, r = S.route("POST", "/api/command", {}, {"lines": MODEL[:-1]})     # up to HOUSE (inline jobs below)
    assert st == 200 and r["ok"]
    while S.jobs.run_pending():
        pass
    st, r = S.route("POST", "/api/command", {}, {"lines": ["RUNANALYS"]})
    assert r["ok"] and r["job"]["state"] == "starting"
    while S.jobs.run_pending():
        pass
    job = max(S.jobs.jobs.values(), key=lambda j: j.id)
    d = job.to_dict()
    assert d["state"] == "done" and d["step"]["key"] == "ANALYS.solve" and d["step"]["data"]["nf"] == 12
    # an INP run in the session itself: progress events of kind "step"
    (Path(S.interp.cwd) / "again.pre").write_text("RUNANALYS\n", encoding="utf-8")   # MDL made run/ the cwd
    seq0 = S.events.since(0)[1]
    st, r = S.route("POST", "/api/command", {}, {"lines": ["INP,again.pre"]})
    assert r["ok"]
    evs, _, _ = S.events.since(seq0)
    steps = [e for e in evs if e["type"] == "progress" and e.get("kind") == "step"]
    assert steps and {e["module"] for e in steps} == {"ANALYS"}
    assert {e["key"] for e in steps} >= {"ANALYS.impedance", "ANALYS.solve"}


def test_worker_wrapper_passes_steps_only_to_runners_that_take_them():
    from sassi.ui.worker import _with_progress
    sent = []

    def with_step(*a, progress=None, step=None):
        step("ANALYS.solve", {"q": 1})
        return 0

    def without_step(*a, progress=None):
        return 0

    assert _with_progress(with_step, sent.append)() == 0 and sent == [{"t": "step", "key": "ANALYS.solve", "data": {"q": 1}}]
    assert _with_progress(without_step, sent.append)() == 0          # no TypeError for an old runner


def test_the_panel_knows_every_step_and_its_theory_section():
    from sassi.ui.markdown import slugify
    keys = set()
    for p in (ROOT / "sassi" / "modules").glob("*.py"):
        keys |= set(re.findall(r'ctx\.announce\("([A-Z]+\.[a-z]+)"', p.read_text(encoding="utf-8")))
    assert len(keys) >= 18
    js = (STATIC / "activity.js").read_text(encoding="utf-8")
    steps = dict(re.findall(r'"([A-Z]+\.[a-z]+)": \{short: "[^"]+", title: "[^"]+", ref: \["[\d.]+", "([\w-]+)"\]', js))
    assert keys <= set(steps), keys - set(steps)
    strip = set(re.findall(r'"([A-Z]+\.[a-z]+)"', js[js.index("const MODULE_STEPS"):js.index("A.STEPS = STEPS;")]))
    assert keys <= strip
    theory = (ROOT / "docs" / "theory" / "THEORY_MANUAL.md").read_text(encoding="utf-8")
    anchors = {slugify(h) for h in re.findall(r"^#{2,3} (.+)$", theory, flags=re.M)}
    assert set(steps.values()) <= anchors, set(steps.values()) - anchors
    assert "S.D.helpTab({doc: TH, anchor: st.ref[1]})" in js


def test_front_end_scripts_hold_no_control_characters():
    """A TeX command such as \\frac or \\right written through a non-raw string becomes a control character
    (form feed, carriage return ...): the equations then do not render."""
    for p in list(STATIC.glob("*.js")) + [ROOT / "web" / "boot.js", ROOT / "web" / "worker.js", ROOT / "web" / "sw.js"]:
        bad = [c for c in p.read_bytes() if c < 0x20 and c not in (0x0A,)]
        assert not bad, f"{p.name}: control bytes {sorted(set(bad))}"
