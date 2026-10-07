"""The JSON API of the SASSI-EDU GUI (ARCHITECTURE section 9) on top of one interpreter.

:class:`GuiSession` holds **one** :class:`sassi.prep.Interpreter` per server session (UI-01): the
Command Entry, the menus, the dialogs and the toolbars all submit command text to it (rule L17),
so the Command History of a GUI session replays as a ``.pre`` file.  The session also subscribes to

* the interpreter's message sink (the six message classes of the Command History, 5.6),
* the plot state (:mod:`sassi.plotting.state`; the GUI draws the plots itself, so
  ``auto_render`` is switched off, and commands given without arguments open their dialog in the
  browser instead of failing, D-UI-08),
* the INP progress hook (status-bar progress: current line / total lines, spec 04 section 4.2),

and appends what it receives to an :class:`~sassi.ui.events.EventLog` that the browser polls.

The same session serves the browser version of the GUI (GitHub Pages, :mod:`sassi.web.bridge`), where
Python runs in the visitor's browser (Pyodide) without processes or threads: ``GuiSession(web=True)``
runs the module jobs inline (:class:`~sassi.ui.jobs.InlineJobManager`) and never blocks in
``/api/events`` (the timeout is ignored; the bridge pushes events to the page instead).

Routes (``route(method, path, query, body)``; all answers are JSON)::

    GET  /api/state                     models, active model, cwd, plots, settings, history
    GET  /api/events?since=N&timeout=T  long poll of the event log
    POST /api/command  {"line"|"lines"} execute command text; returns the messages
    GET  /api/model?number=N            model for 3D plotting (nodes, elements by group, flags ...)
    GET  /api/options/<NAME>            dialog layout + current values (ANALYSIS, MODEL, WRITE, CHECK,
                                        one tab EQUAKE ... AFWRITE, or the Modules-menu dialogs LOADGEN /
                                        LOADGENDYN of ANSYS Eq. Static Load / ANSYS Dynamic Load)
    POST /api/options/<NAME>            dialog commit: validated command text, executed (dry_run: not);
                                        LOADGEN / LOADGENDYN with "run": true then start RUNLOADGEN
    POST /api/options/POINT/radius      the POINT tab "From mesh" helper (RADIUS)
    POST /api/run/<MODULE> {"args"}     RUN<MODULE> in a worker process; returns the job (a RUN<MODULE>
                                        typed in Command Entry runs the same way, UI-04); LOADGEN takes
                                        "args": ["STATIC" | "DYNAMIC"]
    GET  /api/jobs, /api/jobs/<id>?since=k, POST /api/jobs/<id>/cancel
    GET  /api/files?dir=, /api/file?name=, /api/fileinfo?name=; POST /api/file {name, text}
    GET  /api/plots, /api/plot/<id>[?frame=k], /api/lines, /api/dynp, /api/animations, /api/cuts
    GET  /api/harmonic                  FILE8-type files of the active model: frequencies, largest response
    POST /api/harmonic/plan             {"file", "freq", "relative"} -> {"lines": [HARMFRAME ...], "folder"}
    POST /api/harmonic/show             {"folder", "title"} -> {"lines": [PROCFRAME ..., DEFORMPLOT ...]}
    GET  /api/run_summary[?model=]      key inputs, key outputs and charts of a model's run (sassi/ui/runsummary.py)
    POST /api/run_summary/plot          {"model", "chart"} -> {"lines": [READSPEC ..., SPECPLOT ...]}
                                        (/api/dynp: the model's DYNP properties and the built-in library
                                        curves; /api/file and /api/fileinfo also read built-in @ names)
    POST /api/export_table {name}       File > Export Table (CSV of the active 2D plot, D-UI-17)
    GET  /api/check                     the Check Errors window text (<model>.err, 5.5); a CHECK or
                                        AFWRITE with messages pushes a 'check' event (the window pops
                                        up unless Options > Check > Suppress Error Window)
    POST /api/verify {"select"}         Help > Verification (VERIFY in a worker); GET /api/verify/list
    GET  /api/about
    GET  /api/help                      Help > Help: the documents (docs/, examples/README.md) and the
                                        live command index; /api/help/doc?name=<id> one document
                                        rendered to HTML; /api/help/search?q= sections containing q
    GET  /api/settings; POST /api/settings/<locations|extensions|display|colours>
    GET  /api/db; POST /api/db {"action": ...}          Load Model dialog (SASSIdb.xml)
    POST /api/exit                      Model > Exit (saves SASSIini.xml, lists unsaved models)

  Learn (the guided course, :mod:`sassi.ui.learn`; the browser submits the returned command text, L17):

    GET  /api/lessons                   course outline by part (+ files that do not parse), course root
    GET  /api/lessons/<id>              one lesson rendered for the lesson player
    POST /api/lessons/<id>/open         {"confirm"} fresh workspace <course root>/<id>/; returns the commands
                                        (CD, fresh model, sassi-setup) -- or "needs_confirm" and the unsaved
                                        models when the learner must be asked first
    POST /api/lessons/<id>/action       {"verb", "args"} a step action -> command text or a GUI request
    POST /api/lessons/<id>/progress     {"step"} a step has run in the current workspace
    GET  /api/examples                  the examples gallery (examples/*.pre headers + README rows)
    POST /api/examples/<name>/prepare   {"run", "model", "confirm"} fresh <course root>/examples/<name>/; returns
                                        the commands (CD, fresh model[, INP]) and the copied .pre; "model":
                                        the model part of the example and MODELPLOT / LAYERPLOT ("view")
    GET  /api/explain?line=&cursor=     the command explainer (name, syntax, tier, status, arguments)
    GET  /api/course/workspaces         the course workspaces on disk with their size
    POST /api/course/workspaces/delete  {"names"?} delete course workspaces not in use (Learn > Free Disk Space)
    POST /api/settings/course           {"course": {"root"}} Learn > Course Workspace Folder
"""
from __future__ import annotations

import json
import platform
import re
import secrets
import sys
import threading
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from .. import PRODUCT, __version__
from ..io import decks
from ..io import library as LIB
from ..model.values import NumberError, parse_int
from ..plotting.state import FrameStore, PlotError, animation_frame, jsonable, plot_state
from ..prep import Interpreter
from ..prep.lexer import LexError, is_comment, join_command, split_args, split_head
from ..prep.messages import Kind
from ..prep.registry import CommandError, lookup
from . import dialogs, files, harmonic, learn, modeldata, runsummary
from . import explain as explainer
from . import lessons as lessonfiles
from .events import EventLog
from .helpdocs import HOME, HelpDocs, HelpError
from .jobs import InlineJobManager, JobBusy, JobManager, parse_verify
from .settings import SHADER_FIELDS, IniSettings, ModelDatabase, model_files

TITLE = "SASSI-EDU User Interface (ACS SASSI V3 methodology)"
RUN_MODULES = ("EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "COMBIN", "MOTION", "RELDISP",
               "STRESS")
#: Modules-menu entries of later tiers: enabled as soon as their RUN<MODULE> command is implemented
#: (not a placeholder of the command catalogue): NONLINEAR (Option NON) and LOADGEN (Option A, the Run of
#: the ANSYS Eq. Static Load / ANSYS Dynamic Load dialogs), both tier P2
OPTIONAL_RUN_MODULES = ("NONLINEAR", "LOADGEN")
#: RUNLOADGEN,[STATIC|DYNAMIC],[model] (also RUNLOADGEN,<model>,[STATIC|DYNAMIC]): the analysis word
LOADGEN_ANALYSES = ("STATIC", "DYNAMIC")
LOCK_TIMEOUT = 60.0


def run_modules() -> Tuple[str, ...]:
    """The modules the Modules menu (and a typed RUN<MODULE>) can run in a worker process."""
    extra = []
    for m in OPTIONAL_RUN_MODULES:
        spec = lookup(f"RUN{m}")
        if spec is not None and not spec.placeholder and spec.handler is not None:
            extra.append(m)
    return RUN_MODULES + tuple(extra)


def loadgen_args(tokens: List[str]) -> Tuple[str, Optional[int]]:
    """``(analysis, model)`` of RUNLOADGEN's arguments: ``[STATIC|DYNAMIC],[model]`` or, in the RUN<MODULE>
    order, ``<model>,[STATIC|DYNAMIC]`` (as the command reads them); raises ValueError with the message."""
    toks = [str(t).strip() for t in tokens[:2]] + ["", ""]
    try:
        model_first = bool(toks[0]) and parse_int(toks[0]) is not None
    except NumberError:
        model_first = False
    word, num = (toks[1], toks[0]) if model_first else (toks[0], toks[1])
    k_word, k_num = (2, 1) if model_first else (1, 2)
    analysis = (word or "STATIC").upper()
    if analysis not in LOADGEN_ANALYSES:
        raise ValueError(f"argument {k_word} must be STATIC or DYNAMIC (got '{word}')")
    if not num:
        return analysis, None
    try:
        return analysis, parse_int(num)[0]
    except NumberError:
        raise ValueError(f"argument {k_num} <model> = '{num}' is not an integer") from None


class ApiError(Exception):
    """An API request that cannot be served: HTTP ``status`` and a message."""

    def __init__(self, status: int, message: str, **extra: Any):
        super().__init__(message)
        self.status = status
        self.message = message
        self.extra = extra


def _q(query: Dict[str, List[str]], key: str, default: Optional[str] = None) -> Optional[str]:
    v = query.get(key)
    return v[0] if v else default


def _int(v: Any, what: str) -> int:
    try:
        return int(v)
    except (TypeError, ValueError):
        raise ApiError(400, f"{what} must be an integer") from None


class GuiSession:
    """One GUI session: interpreter, event log, jobs and settings.

    ``web``: the browser version (module docstring) -- inline jobs, no long poll.
    """

    def __init__(self, cwd: Optional[str] = None, settings_dir: Optional[str] = None, token: Optional[str] = None,
                 web: bool = False):
        self.web = bool(web)
        self.lock = threading.RLock()
        self.events = EventLog()
        sdir = Path(settings_dir) if settings_dir else None
        self.settings = IniSettings(sdir)
        self.db = ModelDatabase(sdir)
        self.interp = Interpreter(cwd=cwd)
        self.start_dir = Path(self.interp.cwd)      # sassi-gui --model-dir: stays an allowed file root
        self.interp.sink.subscribe(self._on_message)
        self.interp.progress = self._on_progress
        self.interp.module_progress = self._on_module_progress
        self.interp.module_step = self._on_module_step
        self._last_step: Tuple[str, float] = ("", 0.0)
        self.plots = plot_state(self.interp)
        self.plots.auto_render = False            # the browser draws the plots (Plotly)
        self.plots.dialog_handler = self._dialog_handler
        self.plots.subscribe(self._on_plot)
        if sdir is not None:
            self.plots.ani_db_path = sdir / "SASSIani.xml"
        for k in SHADER_FIELDS:                   # shader options persist (UI-07, D-UI-07, 5.10)
            setattr(self.plots.shader, k, float(self.settings.shader[k]))
        manager = InlineJobManager if self.web else JobManager
        self.jobs = manager(self._job_message, self._job_progress, self._job_finished)
        self.token = token or secrets.token_hex(16)
        self.revision = 0
        self.started = time.time()
        self._collect: Dict[int, List[Dict[str, Any]]] = {}
        self._last_progress = 0.0
        self._last_module_progress = 0.0
        self.shutdown: Optional[Callable[[], None]] = None
        self.helpdocs = HelpDocs()                # Help > Help: the project documentation (5.3)
        # Learn: the guided course (lesson files, examples) and the lesson of this session
        self.lesson_dir = lessonfiles.LESSON_DIR
        self.examples_dir = lessonfiles.EXAMPLES_DIR
        self.lesson_session: Optional[Dict[str, Any]] = None
        self._routes = self._build_routes()

    # ================================================================== listeners
    def _on_message(self, msg) -> None:
        ev = {"kind": msg.kind.name, "text": msg.text, "source": msg.source, "command": msg.command}
        self.events.push("message", **ev)
        lst = self._collect.get(threading.get_ident())
        if lst is not None:
            lst.append(ev)

    def _on_progress(self, line: int, total: int, source: str, cmd: str = "") -> None:
        """An INP / macro file is about to execute ``line`` of ``total`` (``cmd``): a progress event (at most
        every 0.15 s, always for the first and last line and for the commands that start a module or write and
        check the decks) with the command, its canonical name and its one-line meaning (the command
        explainer), for the status bar and the activity panel (static/activity.js)."""
        now = time.time()
        text = (cmd or "").strip()
        name = ""
        if text and not is_comment(text):
            try:
                spec = lookup(split_head(text)[0])
            except Exception:                    # noqa: BLE001 - a line the interpreter will report itself
                spec = None
            name = spec.name if spec is not None else ""
        notable = name.startswith("RUN") or name in ("AFWRITE", "CHECK", "INP", "VERIFY")
        if line in (1, total) or notable or now - self._last_progress > 0.15:
            self._last_progress = now
            what = ""
            if name:
                try:
                    what = str(explainer.explain(name).get("meaning") or "")
                except Exception:                # noqa: BLE001 - the panel then shows the command only
                    what = ""
            self.events.push("progress", line=line, total=total, source=source,
                             fraction=(line / total) if total else 1.0, text=f"{source}: line {line} / {total}",
                             cmd=text[:160], name=name, what=what[:200])

    def _on_module_progress(self, module: str, fraction: float, text: str) -> None:
        """A module run in this interpreter (RUN<MODULE> inside an INP file) reports its progress (ANALYS:
        frequency k/n ...): a progress event of kind "module", at most every 0.15 s and at its end."""
        now = time.time()
        if fraction >= 1.0 or now - self._last_module_progress > 0.15:
            self._last_module_progress = now
            self.events.push("progress", kind="module", module=str(module).upper(),
                             fraction=max(0.0, min(1.0, float(fraction))), text=str(text))

    def _on_module_step(self, module: str, key: str, data: Dict[str, Any]) -> None:
        """A module run in this interpreter enters a computation step (``ctx.announce``: ANALYS.impedance ...): a
        progress event of kind "step" with its sizes -- at once when the step changes, else at most every 0.12 s
        (a step repeated at every frequency)."""
        now = time.time()
        last_key, last_t = self._last_step
        if key != last_key or now - last_t > 0.12:
            self._last_step = (key, now)
            self.events.push("progress", kind="step", module=str(module).upper(), key=str(key), data=jsonable(data))

    def _on_plot(self, ev) -> None:
        self.events.push("plot", **ev.to_dict())

    def _dialog_handler(self, name: str, context: Dict[str, Any]) -> bool:
        """Commands given without arguments open their dialog in the browser -- except while a
        ``.pre`` file, macro or loop runs, where they fail as in batch mode (D-UI-08)."""
        return not self.interp.batch

    # ================================================================== locking
    @contextmanager
    def locked(self, timeout: Optional[float] = LOCK_TIMEOUT):
        """Hold the interpreter lock; ``timeout`` None waits as long as needed (background work)."""
        if not self.lock.acquire(timeout=-1 if timeout is None else timeout):
            raise ApiError(503, "the interpreter is busy (a command is still running); try again")
        try:
            yield
        finally:
            self.lock.release()

    # ================================================================== commands
    def execute(self, lines: List[str], timeout: Optional[float] = LOCK_TIMEOUT) -> Dict[str, Any]:
        """Execute command lines as typed in Command Entry (keyboard semantics); collect messages.

        UI-04: a ``RUN<MODULE>[,model]`` typed at the keyboard does not run inside the server's
        interpreter (that would block every other request for the whole run, without Cancel or
        progress): it starts the same worker job as the Modules menu.  The lines after it are
        *deferred* -- they run when the job has ended, so the order of the command text is kept
        (a module run and the commands after it never overlap).  Inside INP, macros and FOREACH the
        interpreter is in batch mode and RUN<MODULE> stays synchronous (UI-04).
        """
        tid = threading.get_ident()
        job: Optional[Dict[str, Any]] = None
        deferred: List[str] = []
        with self.locked(timeout):
            self._collect[tid] = []
            results = []
            report0 = self.interp.session.get("check_report")
            try:
                for i, ln in enumerate(lines):
                    run = self._typed_run(ln)
                    if run is None:
                        ok = self.interp.execute(ln)
                        results.append({"line": ln, "ok": ok})
                        continue
                    started = self._start_typed_run(ln, run, list(lines[i + 1:]))
                    results.append({"line": ln, "ok": started is not None, "job": started})
                    if started is not None:
                        job, deferred = started, list(lines[i + 1:])
                        break
            finally:
                msgs = self._collect.pop(tid, [])
            self._after_commands(report0)
        out = {"ok": all(r["ok"] for r in results), "results": results, "messages": msgs,
               "revision": self.revision}
        if job is not None:
            out["job"] = job
            out["deferred"] = deferred
        return out

    def _after_commands(self, report0: Any) -> None:
        """Bookkeeping after commands (lock held): state event, Check Errors pop-up, shader settings."""
        self.revision += 1
        self.events.push("state", revision=self.revision, active=self.interp.active_model,
                         models=modeldata.models_summary(self.interp), cwd=str(self.interp.cwd))
        rep = self.interp.session.get("check_report")
        if rep is not None and rep is not report0 and rep.messages:
            # requirements 5.5 / spec 05a section 3: the Check Errors window pops up after CHECK (and
            # the CHECK of AFWRITE) unless Options > Check > Suppress Error Window is set
            from ..prep.commands.checks import check_options
            self.events.push("check", suppress=bool(check_options(self.interp).suppress_window),
                             summary=rep.summary(), errors=len(rep.errors()))
        self._sync_shader()

    def _sync_shader(self) -> None:
        """Shader options changed (OK of the Shader dialog submits SHADEROPTIONS): SASSIini.xml."""
        cur = {k: float(getattr(self.plots.shader, k)) for k in SHADER_FIELDS}
        if cur != self.settings.shader:
            self.settings.shader = cur
            try:
                self.settings.save()
            except OSError as exc:
                self.interp.sink.emit(Kind.WARNING, f"SASSIini.xml not written: {exc}", source="gui")

    def save_settings(self) -> Path:
        """Write SASSIini.xml with the current shader options (Model > Exit, server stop)."""
        cur = {k: float(getattr(self.plots.shader, k)) for k in SHADER_FIELDS}
        self.settings.shader = cur
        return self.settings.save()

    # ------------------------------------------------------------------ typed RUN<MODULE> (UI-04)
    def _typed_run(self, line: str) -> Optional[Tuple[str, str]]:
        """``(module, head)`` when ``line`` is RUN<MODULE>[,model] typed at the keyboard (full name or
        documented abbreviation, rule L10); None for every other line and in batch mode."""
        if self.interp.batch:
            return None
        text = line.strip()
        if not text or is_comment(text):
            return None
        head, _ = split_head(text)
        if "@" in head or "#" in head:
            return None                          # a substituted command name: the interpreter's path
        spec = lookup(head)
        if spec is None or not spec.name.startswith("RUN") or spec.name[3:] not in run_modules():
            return None
        return spec.name[3:], head

    def _start_typed_run(self, line: str, run: Tuple[str, str], pending: List[str]) -> Optional[Dict[str, Any]]:
        """Start the worker job of a typed RUN<MODULE>[,model]; errors are reported in the Command
        History like the interpreter's (None returned)."""
        module, head = run
        text = line.strip()
        try:
            sub = self.interp.substitute(text)       # variables substituted once, as INP would (L13)
            _, rest = split_head(sub)
            toks = split_args(rest).tokens
        except (CommandError, LexError) as exc:
            self.interp.sink.emit(Kind.ECHO, text, source="keyboard")
            self.interp.sink.emit(Kind.ERROR, f"RUN{module}: {exc}", source="keyboard", command=f"RUN{module}")
            return None
        number: Optional[int] = None
        args: List[str] = []
        try:
            if module == "LOADGEN":
                analysis, number = loadgen_args(toks)
                args = [analysis]
            else:
                arg = toks[0].strip() if toks else ""
                if arg:
                    try:
                        number = parse_int(arg)[0]
                    except NumberError:
                        raise ValueError(f"argument 1 <model> = '{arg}' is not an integer") from None
        except ValueError as exc:
            self.interp.sink.emit(Kind.ECHO, text, source="keyboard")
            self.interp.sink.emit(Kind.ERROR, f"RUN{module}: {exc}", source="keyboard", command=f"RUN{module}")
            return None
        canon = ",".join([f"RUN{module}"] + args + ([str(number)] if number is not None else []))
        try:
            return self.run_module(module, number, args=args, worker_line=sub, echo=text, history=text,
                                   pending=pending, line=canon)
        except ApiError:
            return None                              # already reported in the Command History

    # ================================================================== jobs
    def _job_message(self, job, msg: Dict[str, Any]) -> None:
        try:
            kind = Kind[msg["kind"]]
        except KeyError:
            kind = Kind.INFO
        # the listing streamed into the Command History (D-RUN-05); emitted on the sink directly
        # so that the interpreter's INP counters are not affected
        self.interp.sink.emit(kind, msg["text"], source=f"job {job.id}", command=msg.get("command", ""))

    def _job_progress(self, job) -> None:
        self.events.push("job", job=job.to_dict(with_messages=False))

    def _job_finished(self, job) -> None:
        if job.state == "cancelled":
            how = "before it started" if self.web else "worker process terminated"
            self.interp.sink.emit(Kind.WARNING, f"{job.line}: cancelled by the user ({how})", source=f"job {job.id}")
        if job.ok and job.kind == "module" and job.history:
            with self.lock:
                # replay history (L17): the run goes where it was *started*, before the commands
                # typed while it ran (an AFWRITE typed during the run must not precede it in a .pre)
                hist = self.interp.history
                idx = len(hist) if job.history_index is None else min(max(job.history_index, 0), len(hist))
                hist.insert(idx, job.history)
        self.events.push("job", job=job.to_dict(with_messages=False))
        pending = list(job.pending or [])
        if not pending:
            return
        if job.state == "cancelled":
            self.interp.sink.emit(Kind.WARNING, f"{len(pending)} command(s) after {job.line} not executed (run "
                                                f"cancelled): {' | '.join(pending)[:300]}", source=f"job {job.id}")
            return
        try:
            self.execute(pending, timeout=None)  # the commands typed after RUN<MODULE>, in their order
        except Exception as exc:              # the job thread must not die silently
            self.interp.sink.emit(Kind.ERROR, f"commands after {job.line} not executed: {exc!r}", source=f"job {job.id}")

    def _worker_spec(self, lines: List[str], with_models: bool = True) -> Dict[str, Any]:
        with self.locked():
            spec: Dict[str, Any] = {"cwd": str(self.interp.cwd), "active": self.interp.active_model,
                                    "lines": list(lines), "models": {}, "afwrite": {}}
            if with_models:
                spec["models"] = {str(n): m.to_json() for n, m in self.interp.models.items()}
                afw = self.interp.session.get("afwrite", {})
                spec["afwrite"] = {str(n): {"hash": st.get("hash"), "blocked": dict(st.get("blocked", {}))}
                                   for n, st in afw.items()}
            return json.loads(json.dumps(spec))

    def _refuse(self, echo: str, text: str, status: int) -> None:
        """A run that cannot start: echo + error in the Command History, then the API error."""
        self.interp.sink.emit(Kind.ECHO, echo, source="keyboard")
        self.interp.sink.emit(Kind.ERROR, text, source="keyboard")
        raise ApiError(status, text, reported=True)

    def _external_target(self, module: str, number: Optional[int]):
        """Model run by an external module executable, and the reason it cannot run (or None).

        The RUN<MODULE> preconditions of spec 04 section 15.3 / requirements 2.6, as in
        :func:`sassi.prep.commands.modules_cmd.run_command`: model name and path, the input deck
        written by AFWRITE, no CHECK error for the module in the last AFWRITE.
        """
        n = self.interp.active_model if number is None or number < 0 else int(number)
        if n not in self.interp.models:
            return None, f"model {n} is not in memory"
        m = self.interp.models[n]
        if not m.name or not m.path:
            return m, "Model name/path not defined -- use MDL"
        if module != "COMBIN":
            ie = self.settings.extensions.get(module, [""])[0]
            dp = Path(m.path) / f"{m.name}{ie}" if ie else decks.deck_path(m.path, m.name, module)
            if not dp.exists():
                return m, f"{dp.name} not found -- Run AFWRITE first (with {module} enabled in AOPT)"
            blocked = self.interp.session.get("afwrite", {}).get(n, {}).get("blocked", {})
            if module in blocked:
                return m, (f"{module} has CHECK errors in the last AFWRITE ({blocked[module]}); correct them and "
                           f"run AFWRITE again")
        return m, None

    def run_module(self, module: str, number: Optional[int] = None, *, args: Optional[List[str]] = None,
                   line: Optional[str] = None, worker_line: Optional[str] = None, echo: Optional[str] = None,
                   history: Optional[str] = None, pending: Optional[List[str]] = None) -> Dict[str, Any]:
        """Start ``RUN<MODULE>[,args][,number]`` in a worker process (Modules menu, typed command; UI-04).

        ``args``: the arguments before the model number -- the analysis of RUNLOADGEN (``STATIC`` or
        ``DYNAMIC``; no other module takes any); ``line``: the command of the job (``RUNSITE``,
        ``RUNSITE,1``, ``RUNLOADGEN,DYNAMIC``); ``worker_line``: the text the worker's interpreter executes
        (the typed text after variable substitution); ``echo`` and ``history``: the text shown in the
        Command History and appended to the replay history when the run succeeds (L17); ``pending``:
        commands to execute once the job has ended.
        """
        module = module.upper()
        if module not in run_modules():
            raise ApiError(404, f"unknown module {module} (Modules menu: {', '.join(run_modules())})")
        args = [str(a).strip().upper() for a in (args or [])]
        if module == "LOADGEN":
            if len(args) > 1 or any(a not in LOADGEN_ANALYSES for a in args):
                raise ApiError(400, "RUNLOADGEN: the analysis is STATIC or DYNAMIC")
        elif args:
            raise ApiError(400, f"RUN{module} takes no arguments before the model number")
        line = line or ",".join([f"RUN{module}"] + args + ([str(number)] if number is not None else []))
        echo = echo or line
        with self.locked():
            spec = self._worker_spec([worker_line or line])
            ext = self.settings.external_module(module)
            if ext:
                m, problem = self._external_target(module, number)
                if problem:
                    self._refuse(echo, f"RUN{module}: {problem}", 400)
                ie, oe = self.settings.extensions.get(module, ["", f"_{module.lower()}.out"])
                spec["external"] = {"exe": ext, "dir": m.path,
                                    "stdin": f"{m.name}\n{m.name}{ie}\n{m.name}{oe}\n"}
            return self._start_job("module", line, spec, module, echo=echo, history=history or line,
                                   pending=pending)

    def _start_job(self, kind: str, line: str, spec: Dict[str, Any], module: str = "", echo: Optional[str] = None,
                   history: Optional[str] = None, pending: Optional[List[str]] = None) -> Dict[str, Any]:
        """Echo the command text (it precedes the streamed listing), then start the worker."""
        echo = echo or line
        with self.locked():
            run = self.jobs.running()
            if run is not None:
                self._refuse(echo, f"{line}: job {run.id} ({run.line}) is still running", 409)
            self.interp.sink.emit(Kind.ECHO, echo, source="keyboard")
            try:
                job = self.jobs.start(kind, line, spec, module=module, history=history,
                                      history_index=len(self.interp.history), pending=list(pending or []))
            except JobBusy as exc:
                self.interp.sink.emit(Kind.ERROR, f"{line}: {exc}", source="keyboard")
                raise ApiError(409, str(exc), reported=True) from None
        self.events.push("job", job=job.to_dict(with_messages=False))
        return job.to_dict()

    def verify(self, select: str) -> Dict[str, Any]:
        sel = (select or "P0").strip()
        if not re.fullmatch(r"[A-Za-z0-9_.\-]+", sel):
            raise ApiError(400, "selection: P0, P1, P2, ALL or a VP id")
        line = f"VERIFY,{sel}"
        return self._start_job("verify", line, self._worker_spec([line], with_models=False))

    # ================================================================== state
    def state(self) -> Dict[str, Any]:
        with self.locked():
            m = self.interp.model
            run = self.jobs.running()
            from ..prep.commands.checks import check_options
            co = check_options(self.interp)
            return {
                "product": PRODUCT, "version": __version__, "title": TITLE,
                "models": modeldata.models_summary(self.interp), "active": self.interp.active_model,
                "model": {"name": m.name, "path": m.path, "title": m.title},
                "cwd": str(self.interp.cwd), "revision": self.revision, "events": self.events.last,
                "job": run.to_dict(with_messages=False) if run else None,
                "plots": self.plots.to_dict(include_lines=False),
                "settings": self.settings.to_dict(),
                "history": list(self.interp.history[-1000:]),
                "write_options": dict(self.interp.write_options),
                "check_options": {"show_warnings": co.show_warnings, "show_errors": co.show_errors,
                                  "suppress_window": co.suppress_window, "break_at": co.break_at},
                "run_modules": list(run_modules()),
                "lesson": self._lesson_state(),
                "course_root": str(self.course_root()),
                "web": self.web,
            }

    # ================================================================== options dialogs
    def options_get(self, name: str) -> Dict[str, Any]:
        name = name.upper()
        with self.locked():
            if name in dialogs.TABS:
                fm = dialogs.form("ANALYSIS")
                fm["tabs"] = [t for t in fm["tabs"] if t["name"] == name]
                return {"form": fm, "values": dialogs.values(self.interp, "ANALYSIS")}
            try:
                fm = dialogs.form(name)
            except KeyError:
                raise ApiError(404, f"unknown dialog {name}") from None
            return {"form": fm, "values": dialogs.values(self.interp, name)}

    def options_post(self, name: str, body: Dict[str, Any]) -> Dict[str, Any]:
        name = name.upper()
        if name in ("WRITE", "CHECK"):
            with self.locked():
                try:
                    notes = dialogs.apply_session_dialog(self.interp, name, body)
                except dialogs.DialogError as exc:
                    raise ApiError(422, "the dialog values were refused", problems=exc.problems) from None
                for n in notes:
                    self.interp.sink.emit(Kind.CONFIRM, n, source="dialog")
                return {"ok": True, "commands": [], "notes": notes, "messages": []}
        if name in dialogs.MODULE_DIALOGS or name == "MODEL":
            kind = name
        else:
            kind = "ANALYSIS"
            if name not in ("ANALYSIS",) + tuple(dialogs.TABS):
                raise ApiError(404, f"unknown dialog {name}")
        with self.locked():
            try:
                lines, notes = dialogs.commit_commands(self.interp, body, kind)
            except dialogs.DialogError as exc:
                raise ApiError(422, "the dialog values were refused (UI-06)", problems=exc.problems) from None
            if body.get("dry_run"):
                return {"ok": True, "commands": lines, "notes": notes, "messages": []}
            for n in notes:
                self.interp.sink.emit(Kind.INFO, n, source="dialog")
            res = self.execute(lines) if lines else {"ok": True, "results": [], "messages": []}
        out = {"ok": res["ok"], "commands": lines, "notes": notes, "messages": res["messages"],
               "results": res["results"]}
        if body.get("run") and kind in dialogs.MODULE_DIALOGS and res["ok"]:
            # Run of ANSYS Eq. Static Load / ANSYS Dynamic Load: the dialog's Ok, then RUNLOADGEN in a worker
            try:
                out["job"] = self.run_module("LOADGEN", args=[dialogs.MODULE_DIALOGS[kind][0]])
            except ApiError as exc:
                out["ok"] = False
                out["run_error"] = exc.message
        return out

    def radius(self) -> Dict[str, Any]:
        res = self.execute(["RADIUS"])
        r = self.interp.session.get("radius") if res["ok"] else None
        return {"ok": res["ok"], "radius": r, "messages": res["messages"]}

    # ================================================================== model
    def model(self, number: Optional[str]) -> Dict[str, Any]:
        with self.locked():
            try:
                return modeldata.model_json(self.interp, -1 if number in (None, "") else _int(number, "number"))
            except KeyError as exc:
                raise ApiError(404, str(exc.args[0])) from None

    # ================================================================== files
    def _roots(self) -> List[Path]:
        extra = [self.start_dir]
        root = self.course_root()
        if root.is_dir():
            extra.append(root)                   # the lesson and example workspaces (Learn)
        return files.allowed_roots(self.interp, extra=extra)

    def files(self, directory: Optional[str]) -> Dict[str, Any]:
        with self.locked():
            roots = self._roots()
            d = directory or str(roots[0])
        try:
            return files.listing(Path(d), roots)
        except files.PathError as exc:
            raise ApiError(403 if "refused" in str(exc) else 404, str(exc)) from None

    def _library_file(self, name: Optional[str]) -> Optional[Path]:
        """A built-in input (``@name``, requirements 7.19): read-only, outside the model roots."""
        if not LIB.is_library_name(name):
            return None
        p = LIB.library_path(name)
        if p is None:
            raise ApiError(404, f"{str(name).strip()}: not a built-in input (LIBRARY lists them)")
        return p

    def file(self, name: Optional[str]) -> Dict[str, Any]:
        lib = self._library_file(name)
        if lib is not None:
            return {"name": LIB.PREFIX + lib.name, "path": str(lib), "text": files.read_text(lib),
                    "info": files.classify(lib), "readonly": True}
        with self.locked():
            roots = self._roots()
        try:
            p = files.safe_path(name or "", roots, must_exist=True)
            text = files.read_text(p)
        except files.PathError as exc:
            raise ApiError(403 if "refused" in str(exc) else 404, str(exc)) from None
        return {"name": p.name, "path": str(p), "text": text, "info": files.classify(p)}

    def file_write(self, body: Dict[str, Any]) -> Dict[str, Any]:
        if LIB.is_library_name(str(body.get("name", ""))):
            if body.get("text") is None:              # File > Open / the dialogs' Edit: show it, read-only
                lib = self._library_file(str(body.get("name")))
                return {"ok": True, "path": LIB.PREFIX + lib.name, "created": False, "text": files.read_text(lib),
                        "readonly": True}
            raise ApiError(403, f"{str(body.get('name')).strip()}: built-in library files are read-only; save a "
                                f"copy under another name")
        with self.locked():
            roots = self._roots()
        try:
            p = files.safe_path(str(body.get("name", "")), roots, must_exist=False)
            if p.exists() and not p.is_file():
                raise files.PathError(f"{p} is not a file")
            existed = p.exists()
            if body.get("text") is None:              # File > Open of a new name: create it (spec 04 4.5)
                if not existed:
                    files.write_text(p, "")
                return {"ok": True, "path": str(p), "created": not existed, "text": files.read_text(p)}
            files.write_text(p, str(body["text"]))
        except files.PathError as exc:
            raise ApiError(403 if "refused" in str(exc) else 400, str(exc)) from None
        return {"ok": True, "path": str(p), "created": not existed}

    def fileinfo(self, name: Optional[str]) -> Dict[str, Any]:
        lib = self._library_file(name)
        if lib is not None:
            return files.file_info(lib)
        with self.locked():
            roots = self._roots()
        try:
            return files.file_info(files.safe_path(name or "", roots, must_exist=True))
        except files.PathError as exc:
            raise ApiError(403 if "refused" in str(exc) else 404, str(exc)) from None

    # ================================================================== plots
    def plot(self, pid: str, frame: Optional[str]) -> Dict[str, Any]:
        with self.locked():
            p = self.plots.plots.get(_int(pid, "plot id"))
            if p is None:
                raise ApiError(404, f"plot {pid} is not open")
            try:
                if frame not in (None, ""):
                    from ..plotting.state import model_scene
                    model = self.interp.models.get(p.model)
                    if model is None or p.family != "anim":
                        raise ApiError(400, "frames are available for animation plots")
                    v = p.view
                    scene = model_scene(model, color_by=v.color_by, show_dof=v.show_dof, show_mass=v.show_mass,
                                        palettes=self.plots.palettes)
                    data = animation_frame(p, scene, FrameStore(p.params["buffer_dir"]), _int(frame, "frame"))
                    return jsonable(data)
                return self.plots.plot_data(p, self.interp)
            except (PlotError, KeyError, ValueError, OSError) as exc:
                if isinstance(exc, ApiError):
                    raise
                raise ApiError(422, f"plot {pid}: {exc}") from None

    def lines(self) -> Dict[str, Any]:
        with self.locked():
            return {str(n): {"number": n, "name": L.name, "kind": L.kind, "n": L.n, "markers": L.markers}
                    for n, L in sorted(self.plots.lines.items())}

    def dynp(self) -> Dict[str, Any]:
        with self.locked():
            props: Dict[str, List[Dict[str, Any]]] = {}
            for key, rec in self.interp.model.options.entries("DYNP"):
                if not isinstance(key, tuple):
                    continue
                label, no = key
                vals = dialogs.record_values(dialogs.typed(rec))
                props.setdefault(label, []).append(vals)
            return {"properties": {k: sorted(v, key=lambda r: r.get("no", 0)) for k, v in props.items()},
                    # built-in curves (Clay, Sand, Rock; requirements 7.19, D-W5-09): used by their label
                    # while the model does not define it
                    "library": {k: v for k, v in LIB.dynp_curves().items() if k not in props}}

    def animations(self) -> Dict[str, Any]:
        try:
            entries = self.plots.ani_db().entries()
            path = str(self.plots.ani_db().path)
        except OSError as exc:
            raise ApiError(500, f"animation database: {exc}") from None
        return {"path": path, "entries": entries}

    def harmonic(self, action: str, body: Dict[str, Any]) -> Dict[str, Any]:
        """Deformed shape from the analysis results (Load Frame Data, :mod:`sassi.ui.harmonic`): the FILE8-type
        files of the active model, then the command text of a choice (run by the page as typed, L17)."""
        try:
            with self.locked():
                if action == "sources":
                    return harmonic.sources(self.interp)
                if action == "plan":
                    return harmonic.plan(self.interp, str(body.get("file", "")), body.get("freq"),
                                         bool(body.get("relative")))
                return harmonic.show(self.interp, str(body.get("folder", "")), str(body.get("title", "")))
        except harmonic.HarmonicError as exc:
            raise ApiError(exc.status, str(exc)) from None

    def run_summary(self, query: Dict[str, Any], body: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """The summary window after a run (:mod:`sassi.ui.runsummary`); with ``body`` the command text that opens
        one of its charts as a plot."""
        try:
            with self.locked():
                if body is None:
                    model = _q(query, "model")
                    return jsonable(runsummary.summary(self.interp, model))
                lines = runsummary.plot_lines(self.interp, body.get("model"), str(body.get("chart", "")),
                                              self.helpdocs, list(self.plots.lines))
                return {"lines": lines}
        except runsummary.SummaryError as exc:
            raise ApiError(exc.status, str(exc)) from None
        except learn.LearnError as exc:
            raise ApiError(exc.status, str(exc)) from None

    def animation_remove(self, body: Dict[str, Any]) -> Dict[str, Any]:
        """Load Frame Data > Remove Animation: drop the SASSIani.xml entry and (optionally) the frame
        store files PROCFRAME wrote (``frame_*.npy``, ``nodes*.npy``, ``index.json``) -- other files of
        the directory are never touched (spec 06 section 8.5)."""
        d = str(body.get("directory") or "")
        if not d:
            raise ApiError(400, "directory required")
        try:
            removed = self.plots.ani_db().remove(d, delete_files=False)
        except OSError as exc:
            raise ApiError(500, f"animation database: {exc}") from None
        deleted: List[str] = []
        p = Path(d)
        if removed and body.get("delete_files") and (p / "index.json").is_file():
            for f in sorted(p.iterdir()):
                if f.is_file() and (f.name == "index.json" or re.fullmatch(r"(frame|nodes)(_\d+)?\.npy", f.name)):
                    f.unlink()
                    deleted.append(f.name)
            try:
                p.rmdir()
            except OSError:
                pass
        self.events.push("plot", event="animations", plot_id=None, data={"directory": d, "removed": removed})
        return {"ok": removed, "deleted": len(deleted), "entries": self.plots.ani_db().entries()}

    def cuts(self) -> Dict[str, Any]:
        with self.locked():
            store = self.interp.session.get("cuts") or {}
            out = []
            if isinstance(store, dict):
                for k in sorted(store, key=lambda x: int(x) if str(x).isdigit() else 0):
                    try:
                        from ..plotting.state import cut_elements
                        n = len(cut_elements(self.interp, int(k)) or [])
                    except (TypeError, ValueError):
                        n = 0
                    out.append({"cut": int(k) if str(k).isdigit() else k, "elements": n})
            return {"cuts": out}

    def export_table(self, body: Dict[str, Any]) -> Dict[str, Any]:
        from ..plotting.state import export_table_csv
        with self.locked():
            p = self.plots.active_plot
            if p is None:
                raise ApiError(400, "no active plot")
            try:
                path = files.safe_path(str(body.get("name") or f"plot{p.id:02d}.csv"), self._roots())
                out = export_table_csv(self.plots, p, path, self.interp)
            except files.PathError as exc:
                raise ApiError(403, str(exc)) from None
            except PlotError as exc:
                raise ApiError(422, str(exc)) from None
            self.interp.sink.emit(Kind.CONFIRM, f"Export Table: {p.caption} written to {out}", source="dialog")
            return {"ok": True, "path": str(out)}

    # ================================================================== check window
    def check(self) -> Dict[str, Any]:
        from ..prep.commands.checks import check_options
        with self.locked():
            m = self.interp.model
            title = f"CHECK: Errors and Warning for - {m.name or '(unnamed model)'}"
            rep = self.interp.session.get("check_report")
            errfile = Path(m.path) / f"{m.name}.err" if m.name and m.path else None
            if errfile is not None and errfile.is_file():
                text, source = errfile.read_text(encoding="utf-8", errors="replace"), str(errfile)
            elif rep is not None:
                text, source = rep.format(check_options(self.interp)), "last CHECK (no .err file: MDL not set)"
            else:
                text, source = "", ""
            summary = rep.summary() if rep is not None else ""
            return {"title": title, "text": text, "source": source, "summary": summary,
                    "suppress": check_options(self.interp).suppress_window}

    # ================================================================== help, about
    def about(self) -> Dict[str, Any]:
        import numpy
        import scipy
        try:
            import plotly
            pv = plotly.__version__
        except ImportError:                      # pragma: no cover
            pv = "not installed"
        try:
            import matplotlib
            mv = matplotlib.__version__
        except ImportError:                      # pragma: no cover
            mv = "not installed"
        out = {"product": PRODUCT, "version": __version__, "title": TITLE,
               "build": "pure Python (no compiled extensions)",
               "methodology": "flexible-volume SSI after Lysmer et al. (1981); ACS SASSI V3 nomenclature",
               "python": sys.version.split()[0], "numpy": numpy.__version__, "scipy": scipy.__version__,
               "plotly": pv, "matplotlib": mv, "platform": platform.platform(),
               "settings": str(self.settings.path), "web": self.web}
        if self.web:
            # the browser version: Python in the page (Pyodide); Plotly.js is loaded by the page itself
            pyodide = sys.modules.get("pyodide")
            out.update(build=f"pure Python in this browser tab (Pyodide {getattr(pyodide, '__version__', '')})".strip(),
                       plotly="Plotly.js of the page", matplotlib="not used (plots drawn by the page)",
                       settings=f"{self.settings.path} (in this browser tab only)")
        return out

    def help(self) -> Dict[str, Any]:
        """Help > Help (requirements 5.3, F1): the documents of the help set (rendered on request by
        :meth:`help_doc`) and the live command index of the interpreter (abbreviation, tier,
        availability and the first docstring line of every command)."""
        from ..prep.registry import CATALOGUE, REGISTRY
        rows = []
        for name in sorted(set(CATALOGUE) | set(REGISTRY)):
            spec = REGISTRY.get(name) or CATALOGUE.get(name)
            doc = ""
            if spec is not None and spec.handler is not None and spec.handler.__doc__:
                doc = spec.handler.__doc__.strip().splitlines()[0]
            rows.append({"name": name, "abbrev": list(getattr(spec, "abbrev", ()) or ()),
                         "tier": getattr(spec, "tier", ""), "class": getattr(spec, "cls", ""),
                         "available": bool(spec is not None and not spec.placeholder),
                         "summary": doc or getattr(spec, "summary", "")})
        docs = self.helpdocs.catalogue() if self.helpdocs.available else []
        home = HOME if any(d["id"] == HOME for d in docs) else (docs[0]["id"] if docs else "")
        try:
            home_doc = self.helpdocs.document(home) if home else None     # the first page, rendered
        except (HelpError, OSError):
            home_doc = None
        return {"docs": docs, "home": home, "home_doc": home_doc, "root": str(self.helpdocs.root), "commands": rows}

    def help_doc(self, name: Optional[str]) -> Dict[str, Any]:
        """One help document rendered to HTML (:mod:`sassi.ui.markdown`: no raw HTML passes)."""
        try:
            return self.helpdocs.document(name or HOME)
        except HelpError as exc:
            raise ApiError(exc.status, str(exc)) from None
        except OSError as exc:
            raise ApiError(404, f"{name}: {exc}") from None

    def help_search(self, query: Optional[str]) -> Dict[str, Any]:
        return self.helpdocs.search(query or "")

    def verify_list(self) -> Dict[str, Any]:
        res = self.execute(["VERIFY,LIST"])
        return {"ok": res["ok"], "listed": parse_verify(res["messages"])["listed"]}

    # ================================================================== settings
    def settings_post(self, kind: str, body: Dict[str, Any]) -> Dict[str, Any]:
        kind = kind.lower()
        errs: List[str] = []
        if kind == "locations":
            errs = self.settings.set_locations(dict(body.get("locations") or {}))
        elif kind == "extensions":
            self.settings.set_extensions(dict(body.get("extensions") or {}))
        elif kind == "display":
            for g, v in (body.get("display") or {}).items():
                if g in self.settings.display:
                    self.settings.display[g] = bool(v)
        elif kind == "course":
            try:
                root = learn.check_root(str((body.get("course") or {}).get("root", "")), self.start_dir)
            except learn.LearnError as exc:
                raise ApiError(exc.status, str(exc)) from None
            given = str((body.get("course") or {}).get("root", "")).strip()
            self.settings.course["root"] = str(root) if given else ""
        elif kind == "colours":
            for k, v in (body.get("colours") or {}).items():
                if k.upper() in self.settings.colours and re.fullmatch(r"#[0-9a-fA-F]{6}", str(v)):
                    self.settings.colours[k.upper()] = str(v)
        else:
            raise ApiError(404, f"unknown settings {kind}")
        try:
            path = self.settings.save()
        except OSError as exc:
            raise ApiError(500, f"cannot write SASSIini.xml: {exc}") from None
        return {"ok": not errs, "errors": errs, "path": str(path), "settings": self.settings.to_dict()}

    # ================================================================== Load Model database
    def db_get(self) -> Dict[str, Any]:
        self.db.load()
        return {"path": str(self.db.path), "groups": self.db.groups}

    def db_post(self, body: Dict[str, Any]) -> Dict[str, Any]:
        act = str(body.get("action", ""))
        deleted: List[str] = []
        try:
            if act == "add_group":
                self.db.add_group(str(body.get("group", "")))
            elif act == "remove_group":
                g = self.db.remove_group(str(body.get("group", "")))
                if body.get("delete_files"):
                    for m in g["models"]:
                        deleted += self._delete_model_files(m)
            elif act == "add_model":
                self.db.add_model(str(body.get("group", "")), str(body.get("name", "")), str(body.get("path", "")),
                                  str(body.get("title", "")))
            elif act == "remove_model":
                m = self.db.remove_model(str(body.get("group", "")), str(body.get("name", "")),
                                         str(body.get("path", "")))
                if body.get("delete_files"):
                    deleted += self._delete_model_files(m)
            else:
                raise ApiError(400, f"unknown action {act}")
        except ValueError as exc:
            raise ApiError(400, str(exc)) from None
        except OSError as exc:
            raise ApiError(500, str(exc)) from None
        out = self.db_get()
        out["deleted"] = deleted
        return out

    def _delete_model_files(self, m: Dict[str, Any]) -> List[str]:
        out = []
        for p in model_files(m["name"], m["path"]):
            try:
                p.unlink()
                out.append(str(p))
            except OSError:
                pass
        return out

    def db_open_lines(self, name: str, path: str, title: str = "") -> List[str]:
        """Load Model > Open: ``MDL`` + ``RESUME`` (requirements 5.3); a model never saved has no .sdb
        (expected and harmless, spec 04 section 3.3), then only MDL (and TIT) are submitted.

        Tokens holding a comma (a folder ``Smith, J``) are double-quoted (L6/L7, D-PAR-06): the
        command text must name the same directory when the history is replayed.
        """
        lines = [join_command("MDL", [name, path])]
        if (Path(path) / f"{name}.sdb").is_file():
            lines.append("RESUME")
        elif title:
            lines.append(join_command("TIT", [title], text_last=True))
        return lines

    # ================================================================== Learn: guided course
    def course_root(self) -> Path:
        """Folder of the lesson and example workspaces (Learn > Course Workspace Folder)."""
        return learn.course_root(self.settings.course.get("root", ""), self.start_dir)

    def _lessons(self) -> Tuple[List[Any], List[Dict[str, str]]]:
        return learn.load_all(self.lesson_dir)

    def _find_lesson(self, lesson_id: str):
        if not lessonfiles.ID_RE.match(lesson_id or ""):
            raise ApiError(400, f"invalid lesson id {lesson_id!r}")
        for les in self._lessons()[0]:
            if les.id == lesson_id:
                return les
        raise ApiError(404, f"no lesson {lesson_id}")

    def _lesson_state(self) -> Optional[Dict[str, Any]]:
        ls = self.lesson_session
        if not ls:
            return None
        return {"id": ls["id"], "workspace": ls["workspace"], "ran": sorted(ls["ran"]), "opened": ls["opened"]}

    def lessons(self) -> Dict[str, Any]:
        """``GET /api/lessons``: the course outline (start page, Learn menu)."""
        good, bad = self._lessons()
        out = learn.outline(good)
        root = self.course_root()
        for part in out["parts"]:
            for d in part["lessons"]:
                d["workspace_exists"] = (root / d["id"]).is_dir()
        out.update(errors=bad, root=str(root), session=self._lesson_state())
        return out

    def lesson(self, lesson_id: str) -> Dict[str, Any]:
        """``GET /api/lessons/<id>``: the lesson rendered for the lesson player."""
        les = self._find_lesson(lesson_id)
        d = learn.render_lesson(les, self.helpdocs)
        ws = self.course_root() / les.id
        d["workspace"] = str(ws)
        d["workspace_exists"] = ws.is_dir()
        ls = self.lesson_session
        d["session"] = self._lesson_state() if ls and ls["id"] == les.id else None
        return d

    def _unsaved_outside_course(self) -> List[Dict[str, Any]]:
        """Models changed since their last SAVE, except those of the course workspaces (a lesson or example
        rebuilds them)."""
        root = self.course_root()
        out = []
        with self.locked():
            paths = {n: m.path for n, m in self.interp.models.items()}
        for u in self.exit_info()["unsaved"]:
            p = paths.get(u["number"]) or ""
            if p:
                try:
                    Path(p).resolve().relative_to(root)
                    continue                      # a model of a lesson / example workspace
                except (ValueError, OSError):
                    pass
            out.append(u)
        return out

    def _prepare_root(self) -> Path:
        root = self.course_root()
        try:
            root.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise ApiError(500, f"cannot create the course workspace folder {root}: {exc}") from None
        return root

    def _refuse_while_running(self, what: str) -> None:
        run = self.jobs.running()
        if run is not None:
            raise ApiError(409, f"{what}: {run.line} is still running -- wait for it to finish or Cancel it")

    def lesson_open(self, lesson_id: str, body: Dict[str, Any]) -> Dict[str, Any]:
        """``POST /api/lessons/<id>/open``: prepare a fresh workspace and return the commands that the browser
        runs through the interpreter: ``CD`` there, a fresh model, the lesson's ``sassi-setup`` block."""
        les = self._find_lesson(lesson_id)
        self._refuse_while_running(f"open lesson {les.id}")
        unsaved = self._unsaved_outside_course()
        if unsaved and not body.get("confirm"):
            return {"needs_confirm": True, "unsaved": unsaved, "lesson": les.id}
        root = self._prepare_root()
        try:
            ws = learn.prepare_lesson(les, root, Path(self.examples_dir))
        except learn.LearnError as exc:
            raise ApiError(exc.status, str(exc)) from None
        with self.locked():
            lines = [learn.cd_line(ws)] + learn.fresh_model_lines(list(self.interp.models), self.interp.active_model)
        lines += list(les.setup)
        self.lesson_session = {"id": les.id, "workspace": str(ws), "ran": set(), "opened": time.time()}
        return {"ok": True, "lesson": les.id, "workspace": str(ws), "commands": lines, "setup": list(les.setup)}

    def lesson_progress(self, lesson_id: str, body: Dict[str, Any]) -> Dict[str, Any]:
        """``POST /api/lessons/<id>/progress {"step": k}``: step k has run in the current workspace."""
        ls = self.lesson_session
        if not ls or ls["id"] != lesson_id:
            raise ApiError(409, f"lesson {lesson_id} is not the open lesson of this session")
        step = _int(body.get("step"), "step")
        if body.get("clear"):
            ls["ran"].discard(step)
        else:
            ls["ran"].add(step)
        return {"ok": True, "session": self._lesson_state()}

    def lesson_action(self, lesson_id: str, body: Dict[str, Any]) -> Dict[str, Any]:
        """``POST /api/lessons/<id>/action {"verb", "args"}``: resolve a step action (paths relative to the
        lesson workspace)."""
        les = self._find_lesson(lesson_id)
        ls = self.lesson_session
        ws = Path(ls["workspace"]) if ls and ls["id"] == les.id else self.course_root() / les.id
        try:
            with self.locked():
                lines_mem = list(self.plots.lines)
                return learn.resolve_action(str(body.get("verb", "")), str(body.get("args", "")), ws, self.interp,
                                            lines_mem, self.helpdocs)
        except learn.LearnError as exc:
            raise ApiError(exc.status, str(exc)) from None

    def examples(self) -> Dict[str, Any]:
        """``GET /api/examples``: the examples gallery."""
        good, _ = self._lessons()
        return {"examples": learn.examples(Path(self.examples_dir), good), "root": str(self.course_root())}

    def example_prepare(self, name: str, body: Dict[str, Any]) -> Dict[str, Any]:
        """``POST /api/examples/<name>/prepare {"run", "model", "confirm"}``: copy the example into a fresh
        workspace; returns the commands (CD there, fresh model, and ``INP`` when ``run``; with ``model`` the
        model part of the example and its view command, :func:`sassi.ui.learn.example_model_lines`) and the
        copied ``.pre``."""
        self._refuse_while_running(f"load example {name}")
        unsaved = self._unsaved_outside_course()
        if unsaved and not body.get("confirm"):
            return {"needs_confirm": True, "unsaved": unsaved, "example": name}
        root = self._prepare_root()
        try:
            ws, pre = learn.prepare_example(name, root, Path(self.examples_dir))
        except learn.LearnError as exc:
            raise ApiError(exc.status, str(exc)) from None
        with self.locked():
            lines = [learn.cd_line(ws)] + learn.fresh_model_lines(list(self.interp.models), self.interp.active_model)
        out: Dict[str, Any] = {"ok": True, "example": name, "workspace": str(ws), "pre": str(pre)}
        if body.get("run"):
            lines.append(join_command("INP", [pre.name]))
        elif body.get("model"):
            # the picture of the gallery card: build the model without module runs, then show it
            try:
                model = learn.example_model_lines(pre.read_text(encoding="utf-8", errors="replace"))
            except OSError as exc:
                raise ApiError(500, f"cannot read {pre}: {exc}") from None
            out["view"] = learn.example_view_command(model)
            lines += model + [out["view"]]
        out["commands"] = lines
        return out

    def _in_use_paths(self) -> List[Optional[Path]]:
        with self.locked():
            paths: List[Optional[Path]] = [Path(self.interp.cwd)]
            paths += [Path(m.path) for m in self.interp.models.values() if m.path]
        ls = self.lesson_session
        if ls:
            paths.append(Path(ls["workspace"]))
        return paths

    def course_workspaces(self) -> Dict[str, Any]:
        """``GET /api/course/workspaces``: the course workspaces on disk and their size (Learn > Free Disk
        Space)."""
        root = self.course_root()
        ws = learn.course_workspaces(root, self._in_use_paths())
        return {"root": str(root), "workspaces": ws, "bytes": sum(w["bytes"] for w in ws)}

    def course_workspaces_delete(self, body: Dict[str, Any]) -> Dict[str, Any]:
        """``POST /api/course/workspaces/delete {"names"?}``: delete course workspaces (all that are not in
        use, or the named ones); only folders the course created are ever deleted."""
        self._refuse_while_running("delete course workspaces")
        names = body.get("names")
        if names is not None and (not isinstance(names, list) or not all(isinstance(n, str) for n in names)):
            raise ApiError(400, "names must be a list of workspace names")
        out = learn.delete_workspaces(self.course_root(), names, self._in_use_paths())
        out.update(self.course_workspaces())
        return out

    def explain(self, line: Optional[str], cursor: Optional[str]) -> Dict[str, Any]:
        """``GET /api/explain?line=&cursor=``: the command explainer (:mod:`sassi.ui.explain`)."""
        c = None if cursor in (None, "") else _int(cursor, "cursor")
        return explainer.explain(line or "", c)

    # ================================================================== exit
    def exit_info(self) -> Dict[str, Any]:
        """Models that changed since their last SAVE (Model > Exit prompt, extension)."""
        out = []
        with self.locked():
            for n, m in sorted(self.interp.models.items()):
                if m.is_empty():
                    continue
                saved = None
                if m.name and m.path and (Path(m.path) / f"{m.name}.sdb").is_file():
                    try:
                        from ..io.container import read_container
                        saved = read_container(Path(m.path) / f"{m.name}.sdb").meta.get("model_hash")
                    except (OSError, ValueError):
                        saved = None
                if saved != m.model_hash():
                    out.append({"number": n, "name": m.name})
        return {"unsaved": out}

    def exit(self) -> Dict[str, Any]:
        try:
            path = str(self.save_settings())          # SASSIini.xml at Exit (5.10), shader options included
        except OSError as exc:
            path = f"not written: {exc}"
        run = self.jobs.running()
        if run is not None:
            self.jobs.cancel(run.id)
        if self.shutdown is not None:
            threading.Thread(target=self._delayed_shutdown, daemon=True).start()
        return {"ok": True, "settings": path}

    def _delayed_shutdown(self) -> None:
        time.sleep(0.3)
        if self.shutdown is not None:
            self.shutdown()

    # ================================================================== routing
    def _build_routes(self) -> List[Tuple[str, "re.Pattern", Callable]]:
        R: List[Tuple[str, "re.Pattern", Callable]] = []

        def add(method: str, pattern: str, fn: Callable) -> None:
            R.append((method, re.compile("^" + pattern + "$"), fn))

        add("GET", r"/api/state", lambda m, q, b: self.state())
        add("GET", r"/api/events", lambda m, q, b: self._events(q))
        add("POST", r"/api/command", lambda m, q, b: self._command(b))
        add("GET", r"/api/model", lambda m, q, b: self.model(_q(q, "number")))
        add("GET", r"/api/models", lambda m, q, b: {"models": modeldata.models_summary(self.interp),
                                                    "active": self.interp.active_model})
        add("POST", r"/api/options/POINT/radius", lambda m, q, b: self.radius())
        add("GET", r"/api/options/(\w+)", lambda m, q, b: self.options_get(m.group(1)))
        add("POST", r"/api/options/(\w+)", lambda m, q, b: self.options_post(m.group(1), b))
        add("POST", r"/api/run/(\w+)", lambda m, q, b: self._run(m.group(1), b))
        add("GET", r"/api/jobs", lambda m, q, b: {"jobs": [j.to_dict(with_messages=False)
                                                          for j in self.jobs.jobs.values()]})
        add("GET", r"/api/jobs/(\d+)", lambda m, q, b: self._job(m.group(1), _q(q, "since", "0")))
        add("POST", r"/api/jobs/(\d+)/cancel", lambda m, q, b: self._cancel(m.group(1)))
        add("GET", r"/api/files", lambda m, q, b: self.files(_q(q, "dir")))
        add("GET", r"/api/file", lambda m, q, b: self.file(_q(q, "name")))
        add("POST", r"/api/file", lambda m, q, b: self.file_write(b))
        add("GET", r"/api/fileinfo", lambda m, q, b: self.fileinfo(_q(q, "name")))
        add("GET", r"/api/plots", lambda m, q, b: self.plots.to_dict(include_lines=False))
        add("GET", r"/api/plot/(\d+)", lambda m, q, b: self.plot(m.group(1), _q(q, "frame")))
        add("GET", r"/api/lines", lambda m, q, b: {"lines": self.lines()})
        add("GET", r"/api/dynp", lambda m, q, b: self.dynp())
        add("GET", r"/api/animations", lambda m, q, b: self.animations())
        add("POST", r"/api/animations/remove", lambda m, q, b: self.animation_remove(b))
        add("GET", r"/api/harmonic", lambda m, q, b: self.harmonic("sources", {}))
        add("GET", r"/api/run_summary", lambda m, q, b: self.run_summary(q))
        add("POST", r"/api/run_summary/plot", lambda m, q, b: self.run_summary(q, b or {}))
        add("POST", r"/api/harmonic/plan", lambda m, q, b: self.harmonic("plan", b))
        add("POST", r"/api/harmonic/show", lambda m, q, b: self.harmonic("show", b))
        add("GET", r"/api/cuts", lambda m, q, b: self.cuts())
        add("POST", r"/api/export_table", lambda m, q, b: self.export_table(b))
        add("GET", r"/api/check", lambda m, q, b: self.check())
        add("POST", r"/api/verify", lambda m, q, b: self.verify(str(b.get("select", "P0"))))
        add("GET", r"/api/verify/list", lambda m, q, b: self.verify_list())
        add("GET", r"/api/about", lambda m, q, b: self.about())
        add("GET", r"/api/help", lambda m, q, b: self.help())
        add("GET", r"/api/help/doc", lambda m, q, b: self.help_doc(_q(q, "name")))
        add("GET", r"/api/help/search", lambda m, q, b: self.help_search(_q(q, "q")))
        add("GET", r"/api/settings", lambda m, q, b: self.settings.to_dict())
        add("POST", r"/api/settings/(\w+)", lambda m, q, b: self.settings_post(m.group(1), b))
        add("GET", r"/api/db", lambda m, q, b: self.db_get())
        add("POST", r"/api/db", lambda m, q, b: self._db(b))
        add("GET", r"/api/lessons", lambda m, q, b: self.lessons())
        add("GET", r"/api/lessons/([0-9a-z][0-9a-z-]*)", lambda m, q, b: self.lesson(m.group(1)))
        add("POST", r"/api/lessons/([0-9a-z][0-9a-z-]*)/open", lambda m, q, b: self.lesson_open(m.group(1), b))
        add("POST", r"/api/lessons/([0-9a-z][0-9a-z-]*)/action", lambda m, q, b: self.lesson_action(m.group(1), b))
        add("POST", r"/api/lessons/([0-9a-z][0-9a-z-]*)/progress", lambda m, q, b: self.lesson_progress(m.group(1), b))
        add("GET", r"/api/examples", lambda m, q, b: self.examples())
        add("POST", r"/api/examples/([A-Za-z0-9][A-Za-z0-9_-]*)/prepare", lambda m, q, b: self.example_prepare(m.group(1), b))
        add("GET", r"/api/explain", lambda m, q, b: self.explain(_q(q, "line", ""), _q(q, "cursor")))
        add("GET", r"/api/course/workspaces", lambda m, q, b: self.course_workspaces())
        add("POST", r"/api/course/workspaces/delete", lambda m, q, b: self.course_workspaces_delete(b))
        add("GET", r"/api/exit", lambda m, q, b: self.exit_info())
        add("POST", r"/api/exit", lambda m, q, b: self.exit())
        return R

    def route(self, method: str, path: str, query: Dict[str, List[str]], body: Dict[str, Any]) -> Tuple[int, Any]:
        """Dispatch one API request; returns ``(status, JSON payload)``."""
        allowed = False
        for meth, pat, fn in self._routes:
            mo = pat.match(path)
            if not mo:
                continue
            allowed = True
            if meth != method:
                continue
            try:
                return 200, fn(mo, query, body or {})
            except ApiError as exc:
                d = {"error": exc.message}
                d.update(exc.extra)
                return exc.status, d
            except Exception as exc:          # a bug must not kill the server; report it
                import traceback
                return 500, {"error": f"internal error: {exc!r}", "traceback": traceback.format_exc()}
        if allowed:
            return 405, {"error": f"method {method} not allowed for {path}"}
        return 404, {"error": f"no API route {path}"}

    # ------------------------------------------------------------------ small handlers
    def _events(self, q: Dict[str, List[str]]) -> Dict[str, Any]:
        since = _int(_q(q, "since", "0"), "since")
        timeout = min(max(float(_q(q, "timeout", "0") or 0), 0.0), 30.0)
        if self.web:
            timeout = 0.0                     # one Python engine: a long poll would block it (bridge pushes)
        evs, last, reset = self.events.since(since, timeout=timeout)
        return {"events": evs, "last": last, "reset": reset}

    def _command(self, b: Dict[str, Any]) -> Dict[str, Any]:
        if "lines" in b:
            lines = [str(x) for x in (b.get("lines") or [])]
        elif "line" in b:
            lines = [str(b.get("line"))]
        else:
            raise ApiError(400, "body must hold 'line' or 'lines'")
        if any("\n" in ln or "\r" in ln for ln in lines):
            raise ApiError(400, "one command per line (use 'lines' for several)")
        return self.execute(lines)

    def _run(self, module: str, b: Dict[str, Any]) -> Dict[str, Any]:
        """``POST /api/run/<MODULE> {"args": [...], "model": n}`` (both optional)."""
        args = b.get("args") or []
        if not isinstance(args, list) or not all(isinstance(a, str) for a in args):
            raise ApiError(400, "args must be a list of words")
        num = b.get("model")
        return self.run_module(module, None if num in (None, "") else _int(num, "model"), args=args)

    def _job(self, jid: str, since: Optional[str]) -> Dict[str, Any]:
        j = self.jobs.get(_int(jid, "job id"))
        if j is None:
            raise ApiError(404, f"job {jid} not found")
        return j.to_dict(since=_int(since or "0", "since"))

    def _cancel(self, jid: str) -> Dict[str, Any]:
        j = self.jobs.get(_int(jid, "job id"))
        if j is None:
            raise ApiError(404, f"job {jid} not found")
        return self.jobs.cancel(j.id).to_dict(with_messages=False)

    def _db(self, b: Dict[str, Any]) -> Dict[str, Any]:
        if b.get("action") == "open":
            name, path = str(b.get("name", "")), str(b.get("path", ""))
            if not name or not path:
                raise ApiError(400, "select a model")
            res = self.execute(self.db_open_lines(name, path, str(b.get("title", ""))))
            return {"ok": res["ok"], "messages": res["messages"]}
        return self.db_post(b)
