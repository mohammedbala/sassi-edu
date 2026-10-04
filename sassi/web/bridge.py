"""The bridge between the page and the Python session of the browser version (GitHub Pages).

``web/worker.js`` -- a Web Worker running Pyodide -- creates one :class:`Bridge` and hands it the requests
of the page (``window.SASSI_TRANSPORT`` of ``web/boot.js``, carried by ``postMessage``)::

    bridge = Bridge(push)                              # push(text): post a JSON push to the page
    answer = bridge.handle(method, path, query, body)  # JSON text {"status": 200, "payload": {...}}
    while bridge.run_pending(): ...                    # after the answer has been posted: queued module runs

* **Requests** are the routes of :meth:`sassi.ui.api.GuiSession.route` (the API of the local server,
  ARCHITECTURE section 9), answered by one session in web mode: ``RUN<MODULE>`` and ``VERIFY`` jobs are
  queued by the request that starts them (:class:`sassi.ui.jobs.InlineJobManager`) and run by
  :meth:`Bridge.run_pending` once the answer is on its way; ``/api/events`` never waits.  No token is
  needed: the page and the worker share one origin and nothing else can reach the worker.
* **Pushes.**  A Web Worker answers one message at a time, so while a long command (``INP`` of a lesson
  step) or a module run keeps Python busy no request is answered.  Instead, the new events of the session
  (Command History, status-bar progress, job states, plot events) and the new listing lines of the
  running job are pushed to the page at most every :data:`PUSH_INTERVAL` seconds::

      {"events": [...], "last": 1234, "reset": false, "job": <GET /api/jobs/<id>?since=k answer>}

  The page skips what it already has (events by sequence number, listing lines by index), so a push and a
  later answer of ``/api/events`` or ``/api/jobs/<id>`` never show a line twice.
* **Files** live in the worker's in-memory file system: the workspace :data:`WORK_DIR` (the start
  directory; the course workspaces are under it) and the settings in :data:`SETTINGS_DIR`.  A reload of
  the page starts afresh.

Nothing here needs the browser: the tests drive a :class:`Bridge` in CPython with a fake ``push``
(``tests/unit/test_web.py``) and in Pyodide under Node (``web/test_pyodide.mjs``).
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple
from urllib.parse import parse_qs

from ..ui.api import GuiSession

#: the workspace of the browser version (Pyodide's home directory is /home/pyodide)
WORK_DIR = "/home/pyodide/work"
#: SASSIini.xml, SASSIdb.xml and SASSIani.xml of the browser version (in memory, like the workspace)
SETTINGS_DIR = "/home/pyodide/.sassi-edu"
#: seconds between two pushes while Python is busy
PUSH_INTERVAL = 0.1
#: largest request body accepted (as the local server)
MAX_BODY = 32 * 1024 * 1024


def _dumps(obj: Any) -> str:
    """JSON text as the local server writes it (no NaN; values JSON does not know as strings)."""
    return json.dumps(obj, allow_nan=False, default=str)


class Bridge:
    """One GUI session of the browser version and its request / push protocol (module docstring).

    ``push(text)`` posts a push (JSON text) to the page (None: pushes are only counted in :attr:`pushed`).
    ``cwd`` / ``settings_dir`` default to :data:`WORK_DIR` / :data:`SETTINGS_DIR`.
    """

    def __init__(self, push: Optional[Callable[[str], Any]] = None, cwd: Optional[str] = None,
                 settings_dir: Optional[str] = None):
        work = Path(cwd or WORK_DIR)
        work.mkdir(parents=True, exist_ok=True)
        self.session = GuiSession(cwd=str(work), settings_dir=settings_dir or SETTINGS_DIR, web=True)
        self._push = push
        self.pushed = 0                               # number of pushes sent
        self._sent_seq = 0                            # newest event the page has been sent
        self._sent_msgs: Dict[int, int] = {}          # job id -> listing lines the page has been sent
        self._next_push = 0.0
        self.session.events.subscribe(self._on_event)

    # ------------------------------------------------------------------ requests
    def handle(self, method: str, path: str, query: str = "", body: str = "") -> str:
        """Answer one request of the page: ``{"status": <HTTP status>, "payload": <JSON answer>}`` as text."""
        self._next_push = time.time() + PUSH_INTERVAL     # a quick request is not interrupted by a push
        status, payload = self._route(str(method or "GET").upper(), str(path or ""), str(query or ""), body)
        if path == "/api/events" and status == 200:
            self._sent_seq = max(self._sent_seq, int(payload.get("last", 0)))
        try:
            return _dumps({"status": status, "payload": payload})
        except (ValueError, TypeError) as exc:            # NaN or a value that cannot be encoded
            return _dumps({"status": 500, "payload": {"error": f"response could not be encoded: {exc}"}})

    def _route(self, method: str, path: str, query: str, body: Any) -> Tuple[int, Any]:
        if not path.startswith("/api/"):
            return 404, {"error": f"no API route {path}"}
        data: Dict[str, Any] = {}
        if method == "POST" and body:
            if len(body) > MAX_BODY:
                return 400, {"error": "request body too large"}
            try:
                data = json.loads(body)
            except ValueError:
                return 400, {"error": "request body is not JSON"}
            if not isinstance(data, dict):
                return 400, {"error": "request body must be a JSON object"}
        return self.session.route(method, path, parse_qs(query.lstrip("?")), data)

    # ------------------------------------------------------------------ jobs
    @property
    def pending(self) -> bool:
        """A module run or verification waits for :meth:`run_pending`."""
        return self.session.jobs.pending

    def run_pending(self) -> bool:
        """Run the next queued job to its end (its output pushed on the way, its final state at the end);
        returns True when another job is queued (deferred commands may start one)."""
        jobs = self.session.jobs
        job = jobs.running()                              # the queued job (state 'starting')
        if job is None or not jobs.pending:
            return jobs.pending
        self._next_push = time.time() + PUSH_INTERVAL
        more = jobs.run_pending()
        self.flush(job)
        return more

    # ------------------------------------------------------------------ pushes
    def _on_event(self, ev: Dict[str, Any]) -> None:
        now = time.time()
        if now >= self._next_push:
            self._next_push = now + PUSH_INTERVAL
            self.flush(self.session.jobs.running())

    def flush(self, job=None) -> None:
        """Push the events the page has not been sent and the new listing lines of ``job`` (if any)."""
        out: Dict[str, Any] = {}
        evs, last, reset = self.session.events.since(self._sent_seq)
        if evs:
            out.update(events=evs, last=last, reset=reset)
        if job is not None:
            d = job.to_dict(since=self._sent_msgs.get(job.id, 0))
            out["job"] = d
        if not out:
            return
        try:
            text = _dumps(out)
        except (ValueError, TypeError):
            return                                        # the next poll reports it
        if evs:
            self._sent_seq = last
        if job is not None:
            self._sent_msgs[job.id] = out["job"]["next"]
        self.pushed += 1
        if self._push is not None:
            try:
                self._push(text)
            except Exception:                             # the page is gone: nothing to tell
                pass
