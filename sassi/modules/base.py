"""Common infrastructure of the SSI modules (EQUAKE, SOIL, SITE, POINT, HOUSE, FORCE,
ANALYS, COMBIN, MOTION, RELDISP, STRESS).

Every module is a Python module ``sassi.modules.<name>`` exposing::

    NAME = "SITE"                      # upper-case module name
    def run(ctx: ModuleContext) -> int # 0 success, non-zero failure

and can be executed:

* from the interpreter (``RUNSITE``) through :func:`run_module`;
* in batch exactly like the original executables, reading the three-line ``.inp``
  protocol from standard input (manual section 3.2)::

      python -m sassi.modules.site < SITE.inp
      # SITE.inp:
      # modelname
      # modelname.sit
      # modelname_SITE.out

All files are read/written in the working directory (the model directory).  The module
writes a plain-text *listing* (``<model>_<MODULE>.out``) through ``ctx.listing``.
"""
from __future__ import annotations

import importlib
import io
import sys
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, List, Optional, Sequence, TextIO, Union

from .. import PRODUCT, __version__

MODULE_NAMES = ("EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "COMBIN",
                "MOTION", "RELDISP", "STRESS", "NONLINEAR", "LOADGEN")
DECK_EXT = {"EQUAKE": ".equ", "SOIL": ".soi", "SITE": ".sit", "POINT": ".poi", "HOUSE": ".hou",
            "FORCE": ".frc", "ANALYS": ".anl", "COMBIN": "", "MOTION": ".mot", "RELDISP": ".rdi",
            "STRESS": ".str", "NONLINEAR": ".eql", "LOADGEN": ".lgn"}


class ModuleError(RuntimeError):
    """Fatal module error: the message is written to the listing and the run stops."""


class Listing:
    """Plain-text output listing of a module run (also echoed to a callback)."""

    def __init__(self, path: Optional[Path], echo: Optional[Callable[[str], None]] = None):
        self.path = path
        self._fh: Optional[TextIO] = open(path, "w", encoding="utf-8") if path else None
        self._echo = echo
        self.warnings: List[str] = []
        self.errors: List[str] = []

    def write(self, text: str = "") -> None:
        if self._fh:
            self._fh.write(text + "\n")
        if self._echo:
            self._echo(text)

    def header(self, module: str, title: str = "") -> None:
        self.write("=" * 78)
        self.write(f" {PRODUCT} {__version__}   module {module}")
        if title:
            self.write(f" {title}")
        self.write(f" run started {time.strftime('%Y-%m-%d %H:%M:%S')}")
        self.write("=" * 78)

    def section(self, title: str) -> None:
        self.write("")
        self.write(title)
        self.write("-" * len(title))

    def warning(self, text: str) -> None:
        self.warnings.append(text)
        self.write(f" *** WARNING: {text}")

    def error(self, text: str) -> None:
        self.errors.append(text)
        self.write(f" *** ERROR: {text}")

    def table(self, columns: List[str], rows, fmt: str = "{:>14.6g}") -> None:
        self.write("".join(f"{c:>14s}" for c in columns))
        for r in rows:
            self.write("".join(fmt.format(v) if not isinstance(v, str) else f"{v:>14s}" for v in r))

    def close(self) -> None:
        if self._fh:
            self._fh.close()
            self._fh = None


@dataclass
class ModuleContext:
    """Everything a module needs to run."""

    module: str
    model: str                      # model name
    workdir: Path                   # model directory (all files relative to it)
    deck_path: Optional[Path]       # <model>.<ext> (None for COMBIN)
    listing: Listing
    progress: Callable[[float, str], None] = field(default=lambda frac, msg: None)
    cancelled: Callable[[], bool] = field(default=lambda: False)

    def path(self, name: str) -> Path:
        return self.workdir / name

    def require(self, name: str, producer: str) -> Path:
        """Return the path of a required input file or raise ModuleError naming its producer."""
        p = self.path(name)
        if not p.exists():
            raise ModuleError(f"{name} missing -- run {producer} first")
        return p


def run_module(module: str, model: str, workdir: Union[str, Path], deck_path: Optional[Union[str, Path]] = None,
               listing_path: Optional[Union[str, Path]] = None, echo: Optional[Callable[[str], None]] = None,
               progress: Optional[Callable[[float, str], None]] = None,
               cancelled: Optional[Callable[[], bool]] = None, notes: Optional[Sequence[str]] = None) -> int:
    """Run SSI module ``module`` for ``model`` in ``workdir``; return 0 on success.

    ``notes``: the built-in defaults of blank inputs the run uses (RUN<MODULE> passes them, requirements
    section 7.19, D-W5-11); they are written to the listing after its header."""
    module = module.upper()
    if module not in MODULE_NAMES:
        raise ValueError(f"unknown module {module}")
    workdir = Path(workdir)
    if deck_path is None and DECK_EXT[module]:
        deck_path = workdir / f"{model}{DECK_EXT[module]}"
    if listing_path is None:
        listing_path = workdir / f"{model}_{module}.out"
    listing = Listing(Path(listing_path), echo=echo)
    ctx = ModuleContext(module=module, model=model, workdir=workdir,
                        deck_path=Path(deck_path) if deck_path else None, listing=listing,
                        progress=progress or (lambda f, m: None), cancelled=cancelled or (lambda: False))
    t0 = time.time()
    try:
        mod = importlib.import_module(f"sassi.modules.{module.lower()}")
        listing.header(module)
        if notes:
            listing.section("Built-in defaults of blank inputs used (CHECK Warning EDU-29; "
                            "EDUOPT,DEFAULTS,OFF switches them off)")
            for text in notes:
                listing.write(f"   {text}")
        if ctx.deck_path is not None and not ctx.deck_path.exists():
            raise ModuleError(f"input deck {ctx.deck_path.name} not found -- run AFWRITE first")
        rc = int(mod.run(ctx) or 0)
    except ModuleError as exc:
        listing.error(str(exc))
        rc = 1
    except Exception as exc:  # unexpected: keep the traceback in the listing
        listing.error(f"unexpected failure: {exc!r}")
        listing.write(traceback.format_exc())
        rc = 2
    finally:
        listing.write("")
        listing.write(f" {module} finished with status {'OK' if 'rc' in locals() and rc == 0 else 'FAILED'}"
                      f" in {time.time() - t0:.2f} s; {len(listing.warnings)} warning(s), {len(listing.errors)} error(s)")
        listing.close()
    return rc


def batch_main(module: str, stdin: Optional[TextIO] = None) -> int:
    """Entry point for ``python -m sassi.modules.<name> < <MODULE>.inp``."""
    stdin = stdin or sys.stdin
    lines = [ln.strip() for ln in stdin.read().splitlines() if ln.strip()]
    if len(lines) < 3:
        print(f"{module}: expected three input lines (model, deck, listing)", file=sys.stderr)
        return 2
    model, deck, out = lines[:3]
    workdir = Path.cwd()
    return run_module(module, model, workdir, deck_path=workdir / deck if deck else None,
                      listing_path=workdir / out, echo=print)
