"""Worker of the GUI (requirements UI-04): runs ``RUN<MODULE>`` or ``VERIFY`` off the session's interpreter.

``python -m sassi.ui.worker`` reads one JSON job description from standard input::

    {"cwd": "...", "active": 0, "models": {"0": <SSIModel.to_json()>, ...},
     "afwrite": {"0": {"hash": "...", "blocked": {...}}},     # AFWRITE state of the session
     "lines": ["RUNSITE"],                                     # command text (rule L17)
     "external": null | {"exe": "/path/SITE.exe", "stdin": "model\\nmodel.sit\\nmodel_SITE.out\\n",
                         "dir": "/model/dir"}}

and writes one JSON object per line to standard output:

* ``{"t": "msg", "kind": "INFO", "text": ..., "source": ..., "command": ...}`` -- every message of
  the worker's interpreter (the module listing is streamed as INFO messages by RUN<MODULE>,
  D-RUN-05); the echo of the command itself is not repeated (the GUI has echoed it);
* ``{"t": "progress", "fraction": 0.4, "text": "SITE Mode 1: frequency 3/8"}``;
* ``{"t": "done", "ok": true}`` at the end.

The commands run in an interpreter of their own holding a copy of the session's models, so a
module run neither blocks nor changes the GUI's interpreter (a run only writes files in the model
directory).  Cancel terminates this process (UI-04).

:func:`run_spec` does the work for one job description and reports through a ``send`` callback: the
worker process passes one that writes the JSON lines above; the browser version of the GUI (no
processes there, :class:`sassi.ui.jobs.InlineJobManager`) calls it in its own process with a callback
that feeds the job directly.

Progress: ``RUN<MODULE>`` (:mod:`sassi.prep.commands.modules_cmd`) calls
:func:`sassi.modules.base.run_module` without a progress callback, so :func:`run_spec` -- and only
while it runs -- wraps ``run_module`` to pass one that reports ``ctx.progress`` to the GUI
status bar.  RUNNONLINEAR and RUNLOADGEN (tier P2) run their modules through
:func:`sassi.modules.nonlinear.run_nonlinear` / :func:`sassi.modules.loadgen.run_loadgen`, which are
wrapped the same way when the job runs them.  The original functions are restored when the job ends
(an inline job shares its process with the GUI session).

``external``: Modules > Location may name an external executable (spec 04 section 15.1, e.g. the
commercial module for a cross-check).  It is run directly (no shell) in the model directory with
the three-line batch protocol on its standard input, and its output is streamed.
"""
from __future__ import annotations

import functools
import importlib
import json
import subprocess
import sys
from contextlib import contextmanager
from typing import Any, Callable, Dict, Iterator, List, Sequence, Tuple

#: the reporting callback of :func:`run_spec`: one message dict per call (``{"t": "msg" | "progress" | "done", ...}``)
Send = Callable[[Dict[str, Any]], None]


def _send(obj: Dict[str, Any]) -> None:
    sys.stdout.write(json.dumps(obj) + "\n")
    sys.stdout.flush()


def _with_progress(original, send: Send = _send):
    """``original`` (a module runner with a ``progress`` keyword) reporting its progress to the GUI."""
    def run(*args, **kw):
        if kw.get("progress") is None:
            kw["progress"] = lambda frac, msg: send({"t": "progress", "fraction": float(frac), "text": str(msg)})
        return original(*args, **kw)

    functools.update_wrapper(run, original)
    return run


#: runners of the tier-P2 modules that their RUN command calls (module, function), by command name
_P2_RUNNERS = {"RUNNONLINEAR": ("sassi.modules.nonlinear", "run_nonlinear"),
               "RUNLOADGEN": ("sassi.modules.loadgen", "run_loadgen")}


@contextmanager
def progress_reporting(send: Send = _send, lines: Sequence[str] = ()) -> Iterator[None]:
    """Pass a progress callback to every run_module call while the job runs (see module docstring), and to
    the NONLINEAR / LOADGEN runners when ``lines`` run them (they are imported only then); the original
    functions are put back afterwards."""
    import sassi.modules.base as base
    patched: List[Tuple[Any, str, Any]] = [(base, "run_module", base.run_module)]
    base.run_module = _with_progress(base.run_module, send)
    try:
        from sassi.prep.lexer import split_head
        from sassi.prep.registry import lookup
        for line in lines:
            spec = lookup(split_head(str(line).strip())[0]) if str(line).strip() else None
            target = _P2_RUNNERS.get(spec.name) if spec is not None else None
            if target is None:
                continue
            try:
                mod = importlib.import_module(target[0])
            except Exception:                      # the run reports the missing module itself
                continue
            fn = getattr(mod, target[1])
            if not getattr(fn, "__wrapped__", None):
                patched.append((mod, target[1], fn))
                setattr(mod, target[1], _with_progress(fn, send))
        yield
    finally:
        for mod, name, fn in reversed(patched):
            setattr(mod, name, fn)


def _run_external(ext: Dict[str, Any], send: Send = _send) -> bool:
    exe = ext["exe"]
    send({"t": "msg", "kind": "INFO", "text": f"external module {exe} (Modules > Location) in {ext['dir']}",
          "source": "worker", "command": ""})
    try:
        proc = subprocess.Popen([exe], cwd=ext["dir"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, errors="replace")
    except OSError as exc:
        send({"t": "msg", "kind": "ERROR", "text": f"cannot start {exe}: {exc}", "source": "worker", "command": ""})
        return False
    try:
        proc.stdin.write(ext.get("stdin", ""))
        proc.stdin.close()
    except OSError:
        pass
    for line in proc.stdout:
        send({"t": "msg", "kind": "INFO", "text": line.rstrip("\n"), "source": "external", "command": ""})
    rc = proc.wait()
    if rc != 0:
        send({"t": "msg", "kind": "ERROR", "text": f"{exe} finished with status {rc}", "source": "worker",
              "command": ""})
    return rc == 0


def run_spec(spec: Dict[str, Any], send: Send = _send) -> bool:
    """Run one job description (module docstring) and report through ``send``: the messages of the
    job's interpreter, its progress and ``{"t": "done", "ok": ...}`` at the end.  Returns ``ok``."""
    if spec.get("external"):
        ok = _run_external(spec["external"], send)
        send({"t": "done", "ok": ok})
        return ok

    from sassi.model import SSIModel
    from sassi.prep import Interpreter
    from sassi.prep.messages import Kind, MessageSink

    sink = MessageSink(keep=2000)

    def forward(msg) -> None:
        if msg.kind == Kind.ECHO:
            return                      # the GUI has already echoed the command line
        send({"t": "msg", "kind": msg.kind.name, "text": msg.text, "source": msg.source, "command": msg.command})

    sink.subscribe(forward)
    interp = Interpreter(cwd=spec.get("cwd") or None, sink=sink)
    models = spec.get("models") or {}
    if models:
        interp.models = {int(k): SSIModel.from_json(v) for k, v in models.items()}
        interp.active_model = int(spec.get("active", min(interp.models)))
    afw = spec.get("afwrite") or {}
    if afw:
        interp.session["afwrite"] = {int(k): dict(v) for k, v in afw.items()}
    lines = list(spec.get("lines", []))
    ok = True
    with progress_reporting(send, lines):
        for line in lines:
            ok = interp.execute(line) and ok
    send({"t": "done", "ok": ok})
    return ok


def main() -> int:
    spec = json.loads(sys.stdin.read() or "{}")
    return 0 if run_spec(spec, _send) else 1


if __name__ == "__main__":
    raise SystemExit(main())
