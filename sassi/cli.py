"""``sassi`` text console (requirements UI-02 b).

Usage::

    sassi                         interactive console (readline history, Up/Down recall)
    sassi run model.pre [--quiet] execute a .pre file (INP) and exit
    sassi -c "N,1,0,0,0" ...      execute command lines and exit
    sassi --version

The console and the batch mode drive the same :class:`sassi.prep.Interpreter` as the GUI (UI-01).
Exit status of ``sassi run``: 0 no errors, 1 errors were reported, 2 the file could not be read.
In the console, ``exit`` / ``quit`` or Ctrl-D leave the program (Model > Exit).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import List, Optional

from . import PRODUCT, __version__


def settings_dir() -> Path:
    """Per-user settings directory (requirements section 5.10, D-UI-07)."""
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "SASSI-EDU"
    if os.name == "nt":
        return Path(os.environ.get("APPDATA", str(Path.home()))) / "SASSI-EDU"
    return Path(os.environ.get("XDG_CONFIG_HOME", str(Path.home() / ".config"))) / "sassi-edu"


def _setup_readline():
    try:
        import readline
    except ImportError:  # pragma: no cover - platform dependent
        return None
    hist = settings_dir() / "console_history"
    try:
        hist.parent.mkdir(parents=True, exist_ok=True)
        if hist.exists():
            readline.read_history_file(str(hist))
        readline.set_history_length(2000)
    except OSError:
        pass
    return readline, hist


def _save_history(rl) -> None:
    if not rl:
        return
    readline, hist = rl
    try:
        readline.write_history_file(str(hist))
    except OSError:
        pass


def repl(interp, stdin=None) -> int:
    """Interactive loop: one command line per input line."""
    from .prep.messages import Kind
    rl = _setup_readline() if stdin is None else None
    interp.sink.set_visible(Kind.ECHO, False)      # the terminal already shows the typed line
    print(f"{PRODUCT} {__version__} command console (ACS SASSI V3 command language). "
          f"Type exit to quit.")
    try:
        while True:
            try:
                if stdin is None:
                    line = input(f"Model {interp.active_model}> ")
                else:
                    line = stdin.readline()
                    if not line:
                        break
            except EOFError:
                print()
                break
            except KeyboardInterrupt:
                print()
                continue
            if line.strip().lower() in ("exit", "quit"):
                break
            interp.execute(line)
    finally:
        _save_history(rl)
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    from .prep import Interpreter
    from .prep.messages import console_sink

    ap = argparse.ArgumentParser(prog="sassi", description=f"{PRODUCT} command interpreter")
    ap.add_argument("--version", action="version", version=f"{PRODUCT} {__version__}")
    ap.add_argument("-c", dest="commands", action="append", metavar="LINE",
                    help="execute a command line (repeatable) and exit")
    ap.add_argument("--cwd", default=None, help="initial working directory")
    sub = ap.add_subparsers(dest="action")
    rp = sub.add_parser("run", help="execute a .pre file and exit")
    rp.add_argument("file", help=".pre command file")
    rp.add_argument("--quiet", action="store_true", help="print only warnings, errors and the summary")
    args = ap.parse_args(argv)

    quiet = bool(getattr(args, "quiet", False))
    ui = Interpreter(cwd=args.cwd, sink=console_sink(quiet=quiet))
    if args.action == "run":
        p = Path(args.file).expanduser()
        if not p.is_absolute():
            p = (ui.cwd / p)
        summary = ui.run_file(str(p), resolve=False) if p.exists() else ui.run_file(args.file)
        if quiet:
            print(summary.text())
        if not summary.ok:
            return 2
        return 1 if summary.errors else 0
    if args.commands:
        failures = ui.execute_lines(args.commands)
        return 1 if failures else 0
    return repl(ui)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
