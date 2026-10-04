"""Background jobs of the GUI: module runs and verification (requirements UI-04, 5.3 Modules menu).

``POST /api/run/<MODULE>`` (Modules > EQUAKE ... RELDISP) starts a :class:`Job` that executes the
command text ``RUN<MODULE>`` in a **worker process** (:mod:`sassi.ui.worker`) driven by a
background thread of the server:

* the GUI stays responsive (the interpreter of the session is not blocked: the worker holds a copy
  of the models -- a module run only writes files in the model directory);
* the module listing is streamed line by line into the job (``GET /api/jobs/<id>``), into the
  Command History and into the module-output tab; ``ctx.progress`` drives the status-bar progress;
* **Cancel** terminates the worker process (``POST /api/jobs/<id>/cancel``);
* a successful run is inserted into the interpreter's replay history at the position where it
  was *started* (commands typed during the run come after it), so the session still replays as a
  ``.pre`` file (rule L17);
* commands submitted after a typed ``RUN<MODULE>`` in the same request wait in ``Job.pending`` and
  are executed by the session when the job has ended.

Help > Verification starts ``VERIFY,<selection>`` the same way; :func:`parse_verify` turns its
messages (computed, reference, error, tolerance, pass/fail per check; requirements 6.1) into the
result table.

Only one job runs at a time (module runs of one model directory must not overlap).  Inside a
``.pre`` file RUNxxx stays synchronous (UI-04): INP runs in the session's own interpreter.

**Browser version** (GitHub Pages, :mod:`sassi.web.bridge`): Pyodide has neither processes nor threads,
so :class:`InlineJobManager` (same interface) queues the job in :meth:`~InlineJobManager.start` -- the
request that started it answers at once -- and :meth:`~InlineJobManager.run_pending` runs it afterwards
*in process* with :func:`sassi.ui.worker.run_spec` (a fresh interpreter holding a copy of the models,
as in the worker process).  The messages, the progress and the end reach the session through the
same callbacks, so the Command History, the replay history and the deferred commands behave as with
a worker process.  A running inline job cannot be cancelled (the Python engine is busy with it); a
job that has not started yet can.
"""
from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

PACKAGE_ROOT = Path(__file__).resolve().parents[2]      # directory containing the 'sassi' package
TERMINATE_GRACE = 3.0


class JobBusy(RuntimeError):
    """Another job is still running."""


class Job:
    """One background command (module run or verification)."""

    def __init__(self, jid: int, kind: str, line: str, module: str = ""):
        self.id = jid
        self.kind = kind                 # 'module' | 'verify'
        self.line = line
        self.module = module
        self.state = "starting"          # starting, running, done, failed, cancelled
        self.ok: Optional[bool] = None
        self.messages: List[Dict[str, Any]] = []
        self.progress = 0.0
        self.progress_text = ""
        self.started = time.time()
        self.finished: Optional[float] = None
        self.returncode: Optional[int] = None
        self.listing = ""
        self.external = ""
        self.proc: Optional[subprocess.Popen] = None
        self.cancel_requested = False
        self.history: Optional[str] = None          # replay-history text of a module run (L17)
        self.history_index: Optional[int] = None    # its position: the history length at the start
        self.pending: List[str] = []                # commands to execute after the job (typed RUN)
        self._lock = threading.Lock()

    @property
    def active(self) -> bool:
        return self.state in ("starting", "running")

    def add(self, msg: Dict[str, Any]) -> None:
        with self._lock:
            self.messages.append(msg)

    def to_dict(self, since: int = 0, with_messages: bool = True) -> Dict[str, Any]:
        with self._lock:
            msgs = self.messages[since:] if with_messages else []
            n = len(self.messages)
        d = {"id": self.id, "kind": self.kind, "line": self.line, "module": self.module, "state": self.state,
             "ok": self.ok, "progress": self.progress, "progress_text": self.progress_text,
             "started": self.started, "finished": self.finished, "returncode": self.returncode,
             "listing": self.listing, "external": self.external, "messages": msgs, "next": n,
             "elapsed": (self.finished or time.time()) - self.started, "pending": list(self.pending)}
        if self.kind == "verify":
            d["table"] = parse_verify(self.messages)
        return d


class _Jobs:
    """What both job managers share: the job table, creating a job, messages and the end of a job.

    ``on_message(job, msg)``, ``on_progress(job)`` and ``on_finished(job)`` are the session's
    callbacks (Command History, status bar, replay history).
    """

    def __init__(self, on_message: Callable[[Job, Dict[str, Any]], None],
                 on_progress: Callable[[Job], None], on_finished: Callable[[Job], None]):
        self.jobs: Dict[int, Job] = {}
        self._next = 1
        self._lock = threading.Lock()
        self._on_message = on_message
        self._on_progress = on_progress
        self._on_finished = on_finished

    def running(self) -> Optional[Job]:
        for j in self.jobs.values():
            if j.active:
                return j
        return None

    def get(self, jid: int) -> Optional[Job]:
        return self.jobs.get(jid)

    def _new_job(self, kind: str, line: str, spec: Dict[str, Any], module: str, history: Optional[str],
                 history_index: Optional[int], pending: Optional[List[str]]) -> Job:
        """A new job in state ``starting``; raises :class:`JobBusy` while another one is active.

        ``history`` / ``history_index`` / ``pending`` are set before the job runs (a job that
        fails at once still finds them in the finished callback).
        """
        with self._lock:
            cur = self.running()
            if cur is not None:
                raise JobBusy(f"job {cur.id} ({cur.line}) is still running")
            job = Job(self._next, kind, line, module)
            job.history, job.history_index, job.pending = history, history_index, list(pending or [])
            self._next += 1
            self.jobs[job.id] = job
        if spec.get("external"):
            job.external = spec["external"]["exe"]
        return job

    def _message(self, job: Job, obj: Dict[str, Any]) -> None:
        msg = {"kind": str(obj.get("kind", "INFO")).upper(), "text": str(obj.get("text", "")),
               "source": str(obj.get("source", "")), "command": str(obj.get("command", ""))}
        job.add(msg)
        m = re.search(r"see (\S+\.out)\b|\(listing (\S+\.out)\)", msg["text"])
        if m:
            job.listing = m.group(1) or m.group(2)
        self._on_message(job, msg)

    def _finish(self, job: Job, ok: bool, state: Optional[str]) -> None:
        job.ok = ok
        job.state = state or ("done" if ok else "failed")
        job.finished = time.time()
        if ok:
            job.progress = 1.0
        self._on_finished(job)


class JobManager(_Jobs):
    """Starts, follows and cancels the background jobs of a :class:`sassi.ui.api.GuiSession` (a worker
    process per job, followed by a thread of the server)."""

    def start(self, kind: str, line: str, spec: Dict[str, Any], module: str = "", history: Optional[str] = None,
              history_index: Optional[int] = None, pending: Optional[List[str]] = None) -> Job:
        """Start a worker for ``spec`` (see :mod:`sassi.ui.worker`); raises :class:`JobBusy`.

        ``history`` / ``history_index`` / ``pending`` are set before the worker starts (a job that
        fails at once still finds them in the finished callback).
        """
        job = self._new_job(kind, line, spec, module, history, history_index, pending)
        t = threading.Thread(target=self._run, args=(job, spec), name=f"sassi-job-{job.id}", daemon=True)
        t.start()
        return job

    def cancel(self, jid: int) -> Job:
        job = self.jobs[jid]
        if not job.active:
            return job
        job.cancel_requested = True
        proc = job.proc
        if proc is not None and proc.poll() is None:
            _terminate(proc)
        return job

    def shutdown(self, timeout: float = TERMINATE_GRACE) -> None:
        """Server stop (Ctrl-C, Model > Exit): cancel the running job and make sure its worker process
        group is gone before the server process exits -- the worker leads its own session, so it does
        not receive the terminal's SIGINT, and the delayed SIGKILL thread of :func:`_terminate` dies
        with the server."""
        job = self.running()
        if job is None:
            return
        self.cancel(job.id)
        t0 = time.time()
        while job.proc is None and job.active and time.time() - t0 < timeout:
            time.sleep(0.02)                     # the worker is still being created
        if job.proc is not None and job.proc.poll() is None:
            _terminate(job.proc)                 # (again: a Popen created after the first cancel)
            try:
                job.proc.wait(timeout)
            except subprocess.TimeoutExpired:
                _signal_group(job.proc, signal.SIGKILL if hasattr(signal, "SIGKILL") else signal.SIGTERM)
                try:
                    job.proc.wait(timeout)
                except subprocess.TimeoutExpired:
                    pass
        if job.proc is not None:
            _signal_group(job.proc, signal.SIGKILL if hasattr(signal, "SIGKILL") else signal.SIGTERM)

    # ------------------------------------------------------------------ worker thread
    def _run(self, job: Job, spec: Dict[str, Any]) -> None:
        env = dict(os.environ)
        env["PYTHONPATH"] = str(PACKAGE_ROOT) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
        env["MPLBACKEND"] = "Agg"
        env["PYTHONUNBUFFERED"] = "1"
        cwd = spec.get("cwd") or os.getcwd()
        last_progress = 0.0
        try:
            proc = subprocess.Popen([sys.executable, "-m", "sassi.ui.worker"], stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd, env=env,
                                    text=True, encoding="utf-8", errors="replace", bufsize=1,
                                    start_new_session=(os.name == "posix"))
        except OSError as exc:
            self._message(job, {"t": "msg", "kind": "ERROR", "text": f"cannot start the worker process: {exc}",
                                "source": f"job {job.id}", "command": ""})
            self._finish(job, False, None)
            return
        job.proc = proc
        job.state = "running"
        if job.cancel_requested:              # Cancel pressed while the process was starting
            _terminate(proc)
        try:
            proc.stdin.write(json.dumps(spec))
            proc.stdin.close()
        except (BrokenPipeError, OSError):
            pass
        ok: Optional[bool] = None
        for raw in proc.stdout:
            raw = raw.rstrip("\n")
            if not raw:
                continue
            try:
                obj = json.loads(raw)
            except ValueError:
                obj = None
            if not isinstance(obj, dict) or "t" not in obj:
                # stray output (warnings printed by a library): show it as information
                self._message(job, {"t": "msg", "kind": "INFO", "text": raw, "source": f"job {job.id}", "command": ""})
                continue
            if obj["t"] == "msg":
                self._message(job, obj)
            elif obj["t"] == "progress":
                job.progress = max(0.0, min(1.0, float(obj.get("fraction", 0.0))))
                job.progress_text = str(obj.get("text", ""))
                now = time.time()
                if now - last_progress > 0.2 or job.progress >= 1.0:
                    last_progress = now
                    self._on_progress(job)
            elif obj["t"] == "done":
                ok = bool(obj.get("ok"))
        rc = proc.wait()
        job.returncode = rc
        if job.cancel_requested:
            self._finish(job, False, "cancelled")
            return
        if ok is None:
            self._message(job, {"t": "msg", "kind": "ERROR", "text": f"worker process ended with status {rc}",
                                "source": f"job {job.id}", "command": ""})
            ok = False
        self._finish(job, ok, None)


class InlineJobManager(_Jobs):
    """The jobs of the browser version (module docstring): no process, no thread.

    :meth:`start` only queues the job (state ``starting``) and returns, so the request that started it
    answers at once; :meth:`run_pending` -- called by the web worker after it has posted that answer --
    runs it synchronously with :func:`sassi.ui.worker.run_spec`, with the semantics of the worker
    process (a fresh interpreter holding a copy of the session's models).
    """

    #: seconds between two progress callbacks (status bar), as for a worker process
    PROGRESS_INTERVAL = 0.2

    def __init__(self, on_message: Callable[[Job, Dict[str, Any]], None],
                 on_progress: Callable[[Job], None], on_finished: Callable[[Job], None]):
        super().__init__(on_message, on_progress, on_finished)
        self._queue: List[Tuple[Job, Dict[str, Any]]] = []

    @property
    def pending(self) -> bool:
        """A job waits for :meth:`run_pending`."""
        return bool(self._queue)

    def start(self, kind: str, line: str, spec: Dict[str, Any], module: str = "", history: Optional[str] = None,
              history_index: Optional[int] = None, pending: Optional[List[str]] = None) -> Job:
        """Queue a job for ``spec`` (see :mod:`sassi.ui.worker`); raises :class:`JobBusy`."""
        job = self._new_job(kind, line, spec, module, history, history_index, pending)
        self._queue.append((job, spec))
        return job

    def run_pending(self) -> bool:
        """Run the next queued job to its end; returns True when another job is queued (the commands
        deferred after a ``RUN<MODULE>`` may start one)."""
        if self._queue:
            job, spec = self._queue.pop(0)
            self._run(job, spec)
        return bool(self._queue)

    def cancel(self, jid: int) -> Job:
        """Cancel a job that has not started; a running inline job ends by itself (not interruptible)."""
        job = self.jobs[jid]
        if job.state != "starting":
            return job
        job.cancel_requested = True
        self._queue = [(j, s) for j, s in self._queue if j is not job]
        self._finish(job, False, "cancelled")
        return job

    def shutdown(self, timeout: float = TERMINATE_GRACE) -> None:
        """Drop the queued jobs (nothing runs in the background)."""
        for job, _ in list(self._queue):
            self.cancel(job.id)

    def _run(self, job: Job, spec: Dict[str, Any]) -> None:
        from .worker import run_spec
        job.state = "running"
        if spec.get("external"):
            # Modules > Location: an executable cannot be started from a browser tab
            self._message(job, {"kind": "ERROR", "source": f"job {job.id}", "text": f"external module "
                                f"{job.external} (Modules > Location) cannot run in the browser version"})
            self._finish(job, False, None)
            return
        last_progress = [0.0]
        result: Dict[str, Optional[bool]] = {"ok": None}

        def send(obj: Dict[str, Any]) -> None:
            if obj.get("t") == "msg":
                self._message(job, obj)
            elif obj.get("t") == "progress":
                job.progress = max(0.0, min(1.0, float(obj.get("fraction", 0.0))))
                job.progress_text = str(obj.get("text", ""))
                now = time.time()
                if now - last_progress[0] > self.PROGRESS_INTERVAL or job.progress >= 1.0:
                    last_progress[0] = now
                    self._on_progress(job)
            elif obj.get("t") == "done":
                result["ok"] = bool(obj.get("ok"))

        try:
            run_spec(spec, send)
        except Exception as exc:              # a failing job must not take the session down
            self._message(job, {"kind": "ERROR", "text": f"the job failed: {exc!r}", "source": f"job {job.id}"})
            result["ok"] = False
        job.returncode = 0 if result["ok"] else 1
        self._finish(job, bool(result["ok"]), None)


def _signal_group(proc: subprocess.Popen, sig: int) -> None:
    """Send ``sig`` to the worker's process group (POSIX), else to the worker."""
    if os.name == "posix":
        try:
            os.killpg(proc.pid, sig)
            return
        except (ProcessLookupError, PermissionError, OSError):
            pass
    try:
        proc.send_signal(sig)
    except (ProcessLookupError, OSError):
        pass


def _terminate(proc: subprocess.Popen) -> None:
    """Cancel (UI-04): terminate the worker and everything it started (an external module
    executable runs as a child of the worker).  On POSIX the worker leads its own process group,
    which receives SIGTERM, then SIGKILL after a grace period."""
    _signal_group(proc, signal.SIGTERM)

    def _kill() -> None:
        time.sleep(TERMINATE_GRACE)
        if proc.poll() is None:
            _signal_group(proc, signal.SIGKILL if hasattr(signal, "SIGKILL") else signal.SIGTERM)
    threading.Thread(target=_kill, daemon=True).start()


# ======================================================================================
# VERIFY output -> table (Help > Verification)
# ======================================================================================
_VP_HEAD = re.compile(r"^(VP-[\w.]+): (.*)$")
_VP_CHECK = re.compile(r"^\s+(pass|FAIL)\s+(.*?): computed (\S+), reference (\S+), error (\S+) \((\w*)\), "
                       r"tolerance (\S+)\s*$")
_VP_BOOL = re.compile(r"^\s+(pass|FAIL)\s+(.*)$")
_VP_END = re.compile(r"^\s+(VP-[\w.]+) (PASSED|FAILED) \(([\d.]+) s\)")
_VP_NOTE = re.compile(r"^\s+note: (.*)$")
_VP_LIST = re.compile(r"^(VP-[\w.]+)\s+(P\d|OOS)\s+(.*?)\s+\[(.*)\]\s*$")


def _num(t: str) -> Optional[float]:
    try:
        return float(t)
    except ValueError:
        return None


def parse_verify(messages: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Rows of the verification table from the messages of ``VERIFY`` (sassi.prep.commands.extensions).

    Returns ``{"problems": [{id, title, status, elapsed, checks: [{quantity, computed, reference,
    error, kind, tolerance, passed}], notes}], "summary": str, "listed": [...]}``; ``listed`` holds the
    rows of ``VERIFY,LIST``.
    """
    probs: List[Dict[str, Any]] = []
    listed: List[Dict[str, Any]] = []
    summary = ""
    cur: Optional[Dict[str, Any]] = None
    for m in messages:
        text = m.get("text", "")
        kind = m.get("kind", "")
        if kind == "ERROR":
            mm = re.search(r"(VP-[\w.]+) failed to run: (.*)", text)
            if mm:
                probs.append({"id": mm.group(1), "title": "", "status": "ERROR", "elapsed": None, "checks": [],
                              "notes": [mm.group(2)]})
                cur = None
            continue
        for line in text.splitlines():
            ml = _VP_LIST.match(line)
            if ml and kind == "INFO" and cur is None and "computed" not in line:
                listed.append({"id": ml.group(1), "tier": ml.group(2), "title": ml.group(3),
                               "modules": ml.group(4)})
                continue
            mh = _VP_HEAD.match(line)
            if mh and not line.startswith(" "):
                cur = {"id": mh.group(1), "title": mh.group(2), "status": "RUNNING", "elapsed": None,
                       "checks": [], "notes": []}
                probs.append(cur)
                continue
            if line.startswith("VERIFY "):
                summary = line
                continue
            if cur is None:
                continue
            me = _VP_END.match(line)
            if me:
                cur["status"] = me.group(2)
                cur["elapsed"] = float(me.group(3))
                cur = None
                continue
            mc = _VP_CHECK.match(line)
            if mc:
                cur["checks"].append({"passed": mc.group(1) == "pass", "quantity": mc.group(2),
                                      "computed": _num(mc.group(3)), "reference": _num(mc.group(4)),
                                      "error": _num(mc.group(5)), "kind": mc.group(6), "tolerance": _num(mc.group(7))})
                continue
            mn = _VP_NOTE.match(line)
            if mn:
                cur["notes"].append(mn.group(1))
                continue
            mb = _VP_BOOL.match(line)
            if mb:
                cur["checks"].append({"passed": mb.group(1) == "pass", "quantity": mb.group(2), "computed": None,
                                      "reference": None, "error": None, "kind": "bool", "tolerance": None})
    return {"problems": probs, "summary": summary, "listed": listed}
