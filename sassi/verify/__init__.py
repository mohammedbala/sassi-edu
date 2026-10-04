"""Verification framework (VERIFY command, pytest verification suite, verification report).

A verification problem (VP) is a function registered with :func:`problem` that returns a
:class:`VPResult`.  Each VP compares computed quantities with reference values from
docs/spec/R2_benchmarks.md (published, exact, approximate or derived) and states the
tolerance.  The same functions are run by pytest (tests/verification) and by the
``VERIFY`` command, and they produce the verification report (docs/verification).

Example::

    from sassi.verify import problem, VPResult

    @problem("VP-30", "Response spectra (Nigam-Jennings)", tier="P0", modules=["MOTION"],
             source="R2 G.1")
    def vp30(workdir):
        r = VPResult()
        r.check("PSA/a0 step, zeta=5%", computed, 1.854468, rtol=1e-5)
        return r
"""
from __future__ import annotations

import importlib
import math
import pkgutil
import shutil
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional


@dataclass
class Check:
    quantity: str
    computed: float
    reference: float
    error: float           # relative (or absolute when atol used)
    tolerance: float
    kind: str              # 'rel' or 'abs'
    passed: bool
    note: str = ""


@dataclass
class VPResult:
    checks: List[Check] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    figures: List[str] = field(default_factory=list)   # paths of PNG figures written by the VP
    elapsed: float = 0.0

    def check(self, quantity: str, computed: float, reference: float, rtol: Optional[float] = None,
              atol: Optional[float] = None, note: str = "") -> bool:
        computed = float(computed)
        reference = float(reference)
        if atol is not None:
            err = abs(computed - reference)
            ok = err <= atol and math.isfinite(computed)
            self.checks.append(Check(quantity, computed, reference, err, atol, "abs", ok, note))
        else:
            tol = 1e-6 if rtol is None else rtol
            denom = abs(reference) if reference != 0 else 1.0
            err = abs(computed - reference) / denom
            ok = err <= tol and math.isfinite(computed)
            self.checks.append(Check(quantity, computed, reference, err, tol, "rel", ok, note))
        return self.checks[-1].passed

    def inform(self, quantity: str, computed: float, reference: float, note: str = "") -> None:
        """Record an *informative* comparison (kind 'info'): reported with its relative difference but
        not a pass/fail criterion.  Used where a published approximate reference has been superseded
        by a rigorous one (the superseding check must be an ordinary ``check``)."""
        computed = float(computed)
        reference = float(reference)
        denom = abs(reference) if reference != 0 else 1.0
        err = abs(computed - reference) / denom if math.isfinite(computed) else float("nan")
        self.checks.append(Check(quantity, computed, reference, err, float("nan"), "info", True, note))

    def require(self, quantity: str, condition: bool, note: str = "") -> bool:
        self.checks.append(Check(quantity, float(bool(condition)), 1.0, 0.0 if condition else 1.0, 0.0, "bool",
                                 bool(condition), note))
        return bool(condition)

    @property
    def passed(self) -> bool:
        return bool(self.checks) and all(c.passed for c in self.checks)


def worse(*values: float) -> float:
    """``max`` for error accumulation that propagates NaN.

    Python's built-in ``max(0.0, nan)`` returns 0.0 (NaN compares false), which would silently hide a
    failed computation from a VP criterion.  ``worse`` returns NaN if any argument is NaN, else the
    largest value.  Use it for every ``worst = worse(worst, err)`` accumulation in verification code."""
    out = -math.inf
    for v in values:
        v = float(v)
        if math.isnan(v):
            return math.nan
        out = v if v > out else out
    return out


@dataclass
class Problem:
    id: str
    title: str
    func: Callable[[Path], VPResult]
    tier: str = "P0"
    modules: List[str] = field(default_factory=list)
    source: str = ""
    slow: bool = False


REGISTRY: Dict[str, Problem] = {}


def problem(id: str, title: str, tier: str = "P0", modules: Optional[List[str]] = None, source: str = "",
            slow: bool = False):
    def deco(func: Callable[[Path], VPResult]):
        REGISTRY[id] = Problem(id, title, func, tier, list(modules or []), source, slow)
        return func
    return deco


LOAD_ERRORS: Dict[str, str] = {}


def load_all() -> Dict[str, Problem]:
    """Import every module of sassi.verify.problems (each registers its VPs).

    A problem module that fails to import is recorded in :data:`LOAD_ERRORS` and skipped, so
    one broken module does not disable the other verification problems.
    """
    from . import problems
    for m in pkgutil.iter_modules(problems.__path__):
        name = f"{problems.__name__}.{m.name}"
        try:
            importlib.import_module(name)
        except Exception as exc:  # pragma: no cover - reported by VERIFY
            LOAD_ERRORS[name] = repr(exc)
    return REGISTRY


def run_problem(id: str, workdir: Optional[Path] = None, keep: bool = False) -> VPResult:
    """Run verification problem ``id``.

    With ``workdir`` the problem's files are written there and kept.  Without it a temporary
    directory is used and removed afterwards unless ``keep`` is True (then its path is added to
    the result notes).
    """
    if id not in REGISTRY:
        load_all()
    if id not in REGISTRY:
        raise KeyError(f"verification problem {id} is not registered (import errors: {LOAD_ERRORS})")
    p = REGISTRY[id]
    temporary = workdir is None
    wd = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix=f"sassi_{id}_"))
    wd.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    try:
        res = p.func(wd)
    finally:
        if temporary and not keep:
            shutil.rmtree(wd, ignore_errors=True)
    res.elapsed = time.time() - t0
    if temporary and keep:
        res.notes.append(f"files kept in {wd}")
    return res
