"""Verification manual and command reference of SASSI-EDU (the V&V report generator).

This module turns the registered verification problems (VPs, :mod:`sassi.verify`) into the
**Verification Manual** ``docs/verification/VERIFICATION_MANUAL.md`` and the interpreter's command
catalogue (:mod:`sassi.prep.registry`) into the **Command Reference**
``docs/reference/COMMAND_REFERENCE.md``.

What the verification manual contains (requirements section 6.1 item 4: "each VP report states the
inputs, the computed and reference values, the error and the tolerance"):

* a summary table of every VP (tier, modules, criteria passed, informative comparisons, largest
  error/tolerance ratio, run time, status);
* the coverage matrix modules x VPs (requirements section 6.5);
* the lead decisions of requirements section 7 (the subsections titled "Lead decisions", 7.16 to
  7.18), which explain why some published approximate references are reported as *informative*
  comparisons (``VPResult.inform``, kind ``info``) and which rigorous references supersede them;
* one section per VP: purpose and model (requirements section 6.3/6.4 row and the VP's own
  description), modules verified, reference sources (R1/R2 sections resolved to their titles,
  decisions cited), the table computed / reference / error / tolerance / result, and the notes;
* optional figures (``--figures``): small PNGs drawn with matplotlib's Agg backend (no display).

Tolerances and reference values are never changed here: the report prints what the VPs return.
A VP that raises is reported as ERROR with the last lines of its traceback, never hidden.

Usage (from the project root)::

    .venv/bin/python -m sassi.verify.report                     # run every VP, write the manual
    .venv/bin/python -m sassi.verify.report --skip-slow         # skip VPs registered with slow=True
    .venv/bin/python -m sassi.verify.report --only VP-01,VP-30  # selected VPs -> ./VERIFICATION_REPORT_partial.md
    .venv/bin/python -m sassi.verify.report --figures           # add the summary figure (PNG)
    .venv/bin/python -m sassi.verify.report --commands-only     # only the command reference

The exit status is 0 when every VP that ran passed, 1 otherwise.  The interpreter command
``VERIFYREPORT`` (:mod:`sassi.prep.commands.verifyreport`) calls :func:`generate_manual`.

Disk use (requirements D-W2-11): every VP runs in its own temporary directory, which is deleted
as soon as the VP has finished (``keep=True`` keeps it for inspection).
"""
from __future__ import annotations

import argparse
import datetime
import inspect
import math
import platform
import re
import shutil
import sys
import tempfile
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

from . import LOAD_ERRORS, Check, Problem, VPResult, load_all, run_problem

# --------------------------------------------------------------------------------------
# Locations
# --------------------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]          # project root (sassi/verify/report.py)
DOCS = ROOT / "docs"
REQUIREMENTS_PATH = DOCS / "spec" / "00_requirements.md"
R1_PATH = DOCS / "spec" / "R1_sassi_theory.md"
R2_PATH = DOCS / "spec" / "R2_benchmarks.md"
DEFAULT_MANUAL = DOCS / "verification" / "VERIFICATION_MANUAL.md"
DEFAULT_FIGURE_DIR = DOCS / "verification" / "figures"
DEFAULT_COMMAND_REFERENCE = DOCS / "reference" / "COMMAND_REFERENCE.md"
#: default output of a partial run (``--only`` / ``--tiers``), in the working directory
PARTIAL_REPORT = "VERIFICATION_REPORT_partial.md"

#: the module chain of manual Fig. 1.1, in the order of the coverage matrix
MODULE_COLUMNS: Tuple[str, ...] = ("EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS",
                                   "COMBIN", "MOTION", "RELDISP", "STRESS", "NONLINEAR", "LOADGEN")
#: column of the coverage matrix for everything that is not one of the analysis modules
OTHER_COLUMN = "UI/tools"

#: largest PNG copied from a VP work directory into the figure directory (bytes)
MAX_FIGURE_BYTES = 2_000_000

STATUS_PASS, STATUS_FAIL, STATUS_ERROR, STATUS_SKIPPED = "PASS", "FAIL", "ERROR", "SKIPPED"
STATUS_LETTER = {STATUS_PASS: "P", STATUS_FAIL: "F", STATUS_ERROR: "E", STATUS_SKIPPED: "S"}


# --------------------------------------------------------------------------------------
# Results of one run
# --------------------------------------------------------------------------------------
@dataclass
class VPRun:
    """One verification problem as run (or skipped) by :func:`run_vps`."""
    id: str
    title: str
    tier: str = "P0"
    modules: List[str] = field(default_factory=list)
    source: str = ""
    slow: bool = False
    description: str = ""              # Markdown: the VP's own description (docstring)
    result: Optional[VPResult] = None
    error: str = ""                    # traceback text when the VP raised
    skipped: str = ""                  # reason when the VP was not run
    elapsed: float = 0.0
    figures: List[str] = field(default_factory=list)   # PNG paths copied to the figure directory

    @classmethod
    def from_problem(cls, p: Problem) -> "VPRun":
        return cls(id=p.id, title=p.title, tier=p.tier, modules=list(p.modules), source=p.source,
                   slow=bool(p.slow), description=describe_problem(p))

    @property
    def status(self) -> str:
        if self.skipped:
            return STATUS_SKIPPED
        if self.error or self.result is None:
            return STATUS_ERROR
        return STATUS_PASS if self.result.passed else STATUS_FAIL

    @property
    def checks(self) -> List[Check]:
        return list(self.result.checks) if self.result is not None else []

    @property
    def criteria(self) -> List[Check]:
        """Pass/fail checks (every kind except the informative ``info``)."""
        return [c for c in self.checks if c.kind != "info"]

    @property
    def informative(self) -> List[Check]:
        return [c for c in self.checks if c.kind == "info"]

    @property
    def notes(self) -> List[str]:
        return list(self.result.notes) if self.result is not None else []

    def margin(self) -> Optional[float]:
        """Largest error/tolerance ratio over the numerical criteria (None if there is none).

        A ratio below 1 passes; ``inf`` marks a failed check with zero tolerance or a non-finite
        computed value.  Boolean conditions are excluded (they have no tolerance)."""
        ratios = [check_ratio(c) for c in self.criteria if c.kind in ("rel", "abs")]
        ratios = [r for r in ratios if r is not None]
        return max(ratios) if ratios else None


def check_ratio(c: Check) -> Optional[float]:
    """error/tolerance of one numerical check (``inf`` for a failure with zero tolerance)."""
    if c.kind not in ("rel", "abs"):
        return None
    err, tol = c.error, c.tolerance
    if err is None or (isinstance(err, float) and math.isnan(err)) or not math.isfinite(c.computed):
        return math.inf
    if tol is None or not math.isfinite(tol):
        return None
    if tol <= 0:
        return 0.0 if err == 0 else math.inf
    return float(err) / float(tol)


# --------------------------------------------------------------------------------------
# Selection and execution
# --------------------------------------------------------------------------------------
_ID_RE = re.compile(r"^VP-(\d+)([A-Za-z]*)$", re.IGNORECASE)


def vp_sort_key(vp_id: str):
    """Natural order: VP-01, VP-02a, VP-02b, ..., VP-38, VP-38T, VP-39, ..., then VP-A1, VP-E1 ..."""
    m = _ID_RE.match(vp_id.strip())
    if m:
        return (0, int(m.group(1)), m.group(2).lower(), "")
    return (1, 0, "", vp_id.upper())


def normalise_id(token: str) -> str:
    """'1' -> 'VP-01', 'vp-2a' -> 'VP-02a', 'VP-A1' unchanged (case of the suffix kept)."""
    t = token.strip()
    if not t:
        return t
    if t.isdigit():
        return f"VP-{int(t):02d}"
    m = _ID_RE.match(t)
    if m:
        return f"VP-{int(m.group(1)):02d}{m.group(2)}"
    if not t.upper().startswith("VP-"):
        return "VP-" + t.upper()
    return "VP-" + t[3:].upper()


def resolve_ids(tokens: Iterable[str], registry: Dict[str, Problem]) -> List[str]:
    """Registered ids for user tokens (case-insensitive); unknown tokens raise ``KeyError``."""
    by_upper = {k.upper(): k for k in registry}
    out: List[str] = []
    for tok in tokens:
        if not tok.strip():
            continue
        cand = normalise_id(tok)
        key = by_upper.get(cand.upper()) or by_upper.get(tok.strip().upper())
        if key is None:
            raise KeyError(f"verification problem {tok.strip()} is not registered")
        if key not in out:
            out.append(key)
    return out


def select_problems(registry: Dict[str, Problem], only: Optional[Sequence[str]] = None,
                    tiers: Optional[Sequence[str]] = None, skip_slow: bool = False
                    ) -> List[Tuple[Problem, str]]:
    """``(problem, skip_reason)`` in natural id order.

    ``only`` restricts the run to the given ids; ``tiers`` to the given tiers (P0/P1/P2).  With
    ``skip_slow`` a VP registered with ``slow=True`` is listed with a skip reason instead of being
    run, so the manual still shows it.
    """
    ids = sorted(registry, key=vp_sort_key)
    if only:
        wanted = set(resolve_ids(only, registry))
        ids = [i for i in ids if i in wanted]
    if tiers:
        tset = {t.strip().upper() for t in tiers if t.strip()}
        ids = [i for i in ids if registry[i].tier.upper() in tset]
    out = []
    for i in ids:
        p = registry[i]
        reason = "slow (skipped with --skip-slow)" if (skip_slow and p.slow) else ""
        out.append((p, reason))
    return out


def _safe_name(vp_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", vp_id)


def _collect_figures(paths: Sequence[str], workdir: Path, figdir: Path, vp_id: str) -> List[str]:
    """Copy the small PNG figures a VP wrote (``VPResult.figures``) to ``figdir``."""
    out: List[str] = []
    for k, p in enumerate(paths, start=1):
        src = Path(p)
        if not src.is_absolute():
            src = workdir / src
        try:
            if src.suffix.lower() != ".png" or not src.is_file() or src.stat().st_size > MAX_FIGURE_BYTES:
                continue
            figdir.mkdir(parents=True, exist_ok=True)
            dst = figdir / f"{_safe_name(vp_id)}_{k}_{src.name}"
            shutil.copyfile(src, dst)
            out.append(str(dst))
        except OSError:
            continue
    return out


def _run_one(p: Problem, wd: Path) -> VPResult:
    """Run one problem in ``wd`` (through :func:`sassi.verify.run_problem` when it is the registered
    problem of that id, directly otherwise, e.g. for test problems that are not registered)."""
    from . import REGISTRY
    if REGISTRY.get(p.id) is p:
        return run_problem(p.id, wd)
    wd.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    res = p.func(wd)
    if isinstance(res, VPResult):
        res.elapsed = time.time() - t0
    return res


def _default_progress(k: int, n: int, run: VPRun) -> None:
    extra = f" ({run.skipped})" if run.skipped else f" ({run.elapsed:.1f} s)"
    print(f"[{k:3d}/{n}] {run.id:8s} {run.status}{extra}", file=sys.stderr, flush=True)


def run_vps(selection: Sequence[Tuple[Problem, str]], keep: bool = False, workroot: Optional[Path] = None,
            figdir: Optional[Path] = None,
            progress: Optional[Callable[[int, int, VPRun], None]] = _default_progress) -> List[VPRun]:
    """Run the selected problems, each in its own work directory.

    The work directories live under ``workroot`` (default: a new temporary directory) and are
    deleted after each VP unless ``keep`` is True.  An exception inside a VP is caught and
    reported (status ERROR) so that one broken VP does not stop the report.
    """
    temporary_root = workroot is None
    root = Path(tempfile.mkdtemp(prefix="sassi_vreport_")) if temporary_root else Path(workroot)
    root.mkdir(parents=True, exist_ok=True)
    runs: List[VPRun] = []
    try:
        n = len(selection)
        for k, (p, skip) in enumerate(selection, start=1):
            run = VPRun.from_problem(p)
            if skip:
                run.skipped = skip
            else:
                wd = root / _safe_name(p.id)
                t0 = time.time()
                try:
                    res = _run_one(p, wd)
                    if not isinstance(res, VPResult):
                        raise TypeError(f"{p.id} returned {type(res).__name__}, not VPResult")
                    run.result = res
                except Exception:     # a VP must never stop the report (reported as ERROR)
                    run.error = traceback.format_exc()
                run.elapsed = time.time() - t0
                if figdir is not None and run.result is not None and run.result.figures:
                    run.figures = _collect_figures(run.result.figures, wd, Path(figdir), p.id)
                if not keep:
                    shutil.rmtree(wd, ignore_errors=True)
            runs.append(run)
            if progress is not None:
                progress(k, n, run)
    finally:
        if temporary_root and not keep:
            shutil.rmtree(root, ignore_errors=True)
    return runs


# --------------------------------------------------------------------------------------
# Markdown helpers
# --------------------------------------------------------------------------------------
_CODE_SPAN = re.compile(r"(`+)(.+?)\1", re.DOTALL)
_ROLE = re.compile(r":(?:py:)?(?:func|mod|class|meth|data|attr|obj|exc|const|ref|term):`~?([^`]+)`")


#: reST/Markdown emphasis around words (``*by hand*``): kept as emphasis, not escaped
_EMPHASIS = re.compile(r"(?<![\w*\\])\*([A-Za-z](?:[^*\n]*?[A-Za-z.)])?)\*(?![\w*])")


def _escape_plain(text: str) -> str:
    """Escape Markdown specials in text outside code spans (HTML brackets, stray emphasis markers);
    word emphasis ``*like this*`` is kept."""
    text = text.replace("\\", "\\\\")
    text = re.sub(r"<(?=[A-Za-z/!?])", "&lt;", text)       # would start an HTML tag
    if text.startswith(">"):
        text = "\\" + text                                  # would start a block quote
    keep: List[str] = []

    def _protect(m: "re.Match") -> str:
        keep.append(m.group(1))
        return f"\x00{len(keep) - 1}\x00"

    text = _EMPHASIS.sub(_protect, text)
    text = text.replace("*", "\\*")
    text = re.sub(r"(?<![0-9A-Za-z])_|_(?![0-9A-Za-z])", r"\\_", text)
    return re.sub(r"\x00(\d+)\x00", lambda m: "*" + keep[int(m.group(1))] + "*", text)


def md_escape(text: str) -> str:
    """Escape plain text for Markdown, leaving ``code spans`` untouched."""
    out, pos = [], 0
    for m in _CODE_SPAN.finditer(text):
        out.append(_escape_plain(text[pos:m.start()]))
        out.append(m.group(0))
        pos = m.end()
    out.append(_escape_plain(text[pos:]))
    return "".join(out)


def md_cell(text: object, escape: bool = True) -> str:
    """One table cell: single line, ``|`` escaped (also inside code spans, as GFM requires)."""
    s = "" if text is None else str(text)
    s = " ".join(s.split())
    if escape:
        s = md_escape(s)
    return s.replace("|", "\\|") if escape else re.sub(r"(?<!\\)\|", "\\|", s)


def md_table(headers: Sequence[str], rows: Iterable[Sequence[object]], align: Optional[Sequence[str]] = None,
             escape: bool = True) -> str:
    """GitHub-flavoured Markdown table; ``align`` items are 'l', 'r' or 'c'."""
    hdr = "| " + " | ".join(md_cell(h, escape=False) for h in headers) + " |"
    marks = []
    for k in range(len(headers)):
        a = (align[k] if align and k < len(align) else "l")
        marks.append({"l": "---", "r": "---:", "c": ":---:"}.get(a, "---"))
    lines = [hdr, "|" + "|".join(marks) + "|"]
    for r in rows:
        cells = [md_cell(c, escape=escape) for c in r]
        cells += [""] * (len(headers) - len(cells))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def slug(text: str) -> str:
    """GitHub-style heading anchor (lower case, punctuation removed, blanks to hyphens)."""
    s = text.strip().lower()
    s = re.sub(r"[^\w\- ]", "", s)
    return s.replace(" ", "-")


_BULLET = re.compile(r"^(\s*)([*+-]|\d+[.)]|\(\d+\)|[a-z][.)])\s+(.*)$")


def rst_to_md(text: str) -> str:
    """Convert a reST-flavoured docstring to Markdown paragraphs and lists.

    Roles (``:func:`x```) become code spans, ``double backticks`` single ones, paragraphs are
    re-flowed onto one line each, bullet/enumerated items become list items and literal blocks
    (after ``::``) fenced code.  Text outside code spans is escaped (:func:`md_escape`).
    """
    if not text:
        return ""
    text = _ROLE.sub(lambda m: "``" + m.group(1) + "``", text)
    text = re.sub(r"``([^`]+?)``", r"`\1`", text)
    blocks: List[List[str]] = []
    cur: List[str] = []
    for line in text.splitlines():
        if line.strip():
            cur.append(line.rstrip())
        elif cur:
            blocks.append(cur)
            cur = []
    if cur:
        blocks.append(cur)
    out: List[str] = []
    literal_next = False
    for block in blocks:
        if literal_next:
            indent = min(len(b) - len(b.lstrip()) for b in block)
            out.append("```text\n" + "\n".join(b[indent:] for b in block) + "\n```")
            literal_next = False
            continue
        first = block[0]
        if _BULLET.match(first):
            items: List[Tuple[str, str]] = []
            base_indent = len(first) - len(first.lstrip())
            for line in block:
                m = _BULLET.match(line)
                ind = len(line) - len(line.lstrip())
                if m and ind <= base_indent + 1:
                    marker = m.group(2)
                    kind = "1." if marker[0].isdigit() or marker.startswith("(") else "-"
                    items.append((kind, m.group(3).strip()))
                elif items:
                    items[-1] = (items[-1][0], items[-1][1] + " " + line.strip())
                else:
                    items.append(("-", line.strip()))
            lines = []
            for kind, body in items:
                body = " ".join(body.split())
                lit = body.endswith("::")
                if lit:
                    body = body[:-1]
                    literal_next = True
                lines.append(f"{kind} {md_escape(body)}")
            out.append("\n".join(lines))
            continue
        para = " ".join(" ".join(b.split()) for b in block)
        if para.endswith("::"):
            para = para[:-1]
            literal_next = True
        para = md_escape(para)
        if para.startswith("#"):
            para = "\\" + para
        out.append(para)
    return "\n\n".join(out)


def fmt_num(x: Optional[float]) -> str:
    """Compact number: 6 significant digits, scientific outside [1e-4, 1e6)."""
    if x is None:
        return ""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    if x == 0:
        return "0"
    ax = abs(x)
    if 1e-4 <= ax < 1e6:
        return f"{x:.6g}"
    return f"{x:.4e}"


def fmt_small(x: Optional[float]) -> str:
    """Errors and ratios: 3 significant digits."""
    if x is None:
        return ""
    if isinstance(x, float) and math.isnan(x):
        return "nan"
    if isinstance(x, float) and math.isinf(x):
        return "inf"
    if x == 0:
        return "0"
    return f"{x:.3g}"


def check_row(k: int, c: Check, with_note: bool) -> List[str]:
    """Row of the per-VP check table."""
    if c.kind == "bool":
        row = [str(k), c.quantity, "true" if c.passed else "false", "true", "", "condition",
               "pass" if c.passed else "FAIL"]
    elif c.kind == "info":
        row = [str(k), c.quantity, fmt_num(c.computed), fmt_num(c.reference), fmt_small(c.error) + " rel",
               "informative", "info"]
    else:
        row = [str(k), c.quantity, fmt_num(c.computed), fmt_num(c.reference),
               fmt_small(c.error) + (" rel" if c.kind == "rel" else " abs"),
               fmt_small(c.tolerance) + (" rel" if c.kind == "rel" else " abs"),
               "pass" if c.passed else "FAIL"]
    if with_note:
        row.append(c.note)
    return row


# --------------------------------------------------------------------------------------
# Specification lookups (requirements tables, R1/R2 section titles, decisions)
# --------------------------------------------------------------------------------------
def _read(path: Path) -> str:
    try:
        return Path(path).read_text(encoding="utf-8")
    except OSError:
        return ""


def split_row(line: str) -> List[str]:
    """Cells of a Markdown table row; escaped pipes (``\\|``) stay inside their cell."""
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|") and not s.endswith("\\|"):
        s = s[:-1]
    cells, cur, i = [], [], 0
    while i < len(s):
        ch = s[i]
        if ch == "\\" and i + 1 < len(s) and s[i + 1] == "|":
            cur.append("\\|")
            i += 2
            continue
        if ch == "|":
            cells.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
        i += 1
    cells.append("".join(cur).strip())
    return cells


def _section(text: str, heading_re: str) -> str:
    """Text of the section whose heading matches ``heading_re`` (up to the next heading of the same
    or a higher level)."""
    lines = text.splitlines()
    start, level = None, 0
    for i, line in enumerate(lines):
        m = re.match(r"^(#+)\s+(.*)$", line)
        if not m:
            continue
        if start is None:
            if re.search(heading_re, m.group(2)):
                start, level = i, len(m.group(1))
        elif len(m.group(1)) <= level:
            return "\n".join(lines[start:i])
    return "\n".join(lines[start:]) if start is not None else ""


def parse_tables(text: str) -> List[Tuple[List[str], List[List[str]]]]:
    """Every Markdown table of ``text`` as (header cells, rows of cells)."""
    tables, lines, i = [], text.splitlines(), 0
    while i < len(lines):
        if lines[i].lstrip().startswith("|") and i + 1 < len(lines) and re.match(r"^\s*\|[\s:|-]+\|\s*$",
                                                                                 lines[i + 1]):
            header = split_row(lines[i])
            rows = []
            j = i + 2
            while j < len(lines) and lines[j].lstrip().startswith("|"):
                rows.append(split_row(lines[j]))
                j += 1
            tables.append((header, rows))
            i = j
        else:
            i += 1
    return tables


def requirement_rows(path: Path = REQUIREMENTS_PATH) -> Dict[str, Dict[str, str]]:
    """Rows of the verification tables of requirements sections 6.2-6.4, keyed by id.

    Returns ``{id: {column: cell}}`` with the original Markdown of each cell, plus ``'section'``.
    """
    text = _read(path)
    out: Dict[str, Dict[str, str]] = {}
    for sec_re, label in ((r"^6\.3\b", "6.3"), (r"^6\.4\b", "6.4"), (r"^6\.2\b", "6.2")):
        sec = _section(text, sec_re)
        for header, rows in parse_tables(sec):
            if not header or header[0].strip().upper() != "ID":
                continue
            for r in rows:
                if not r or not r[0].strip():
                    continue
                d = {h.strip(): (r[k] if k < len(r) else "") for k, h in enumerate(header)}
                d["section"] = label
                out[r[0].strip()] = d
    return out


def requirement_row_for(vp_id: str, rows: Dict[str, Dict[str, str]]) -> Tuple[Optional[Dict[str, str]], str]:
    """The requirements row of ``vp_id`` or of its base id (VP-38T -> VP-38, VP-39b -> VP-39).

    Returns (row, relation) with relation '' (same id), 'extends' (base id) or 'none'."""
    if vp_id in rows:
        return rows[vp_id], ""
    m = _ID_RE.match(vp_id)
    if m and m.group(2):
        base = f"VP-{int(m.group(1)):02d}"
        if base in rows:
            return rows[base], "extends"
        for k in rows:          # VP-02a / VP-02b style ids in the table
            if k.upper() == vp_id.upper():
                return rows[k], ""
    return None, "none"


@dataclass
class Decision:
    id: str
    subject: str
    decision: str
    evidence: str = ""
    section: str = ""


def parse_decisions(path: Path = REQUIREMENTS_PATH) -> Dict[str, Decision]:
    """Every decision row ``D-...`` of requirements section 7, keyed by its id."""
    text = _read(path)
    sec7 = _section(text, r"^7\.\s")
    out: Dict[str, Decision] = {}
    current = ""
    for line in sec7.splitlines():
        m = re.match(r"^###\s+(7\.\d+)\s", line)
        if m:
            current = m.group(1)
            continue
        if not line.lstrip().startswith("|"):
            continue
        cells = split_row(line)
        m_id = re.match(r"^(D-[A-Z0-9]+-\d+)", cells[0].strip()) if cells else None
        if m_id is None:
            continue
        did = m_id.group(1)
        if len(cells) >= 4 and current in ("7.16", "7.17"):
            out[did] = Decision(did, cells[1], cells[2], cells[3], current)
        elif len(cells) >= 3:
            # (ID, Question, Decision, Rationale)
            out[did] = Decision(did, cells[1], cells[2], cells[3] if len(cells) > 3 else "", current)
    return out


def lead_section_numbers(text: str) -> List[str]:
    """Numbers of the subsections of requirements section 7 titled "Lead decisions ..." (7.16, 7.17,
    7.18 ...), in document order; 7.16 and 7.17 when there are none."""
    nums = re.findall(r"^###\s+(7\.\d+)\s+Lead decisions\b", text, flags=re.M)
    return nums or ["7.16", "7.17"]


def lead_sections_text(nums: Sequence[str]) -> str:
    """'7.16 and 7.17' or '7.16 to 7.18' for the manual's prose."""
    nums = list(nums) or ["7.16", "7.17"]
    if len(nums) == 1:
        return nums[0]
    return f"{nums[0]} and {nums[1]}" if len(nums) == 2 else f"{nums[0]} to {nums[-1]}"


def lead_decision_sections(path: Path = REQUIREMENTS_PATH) -> List[Tuple[str, str, List[Decision]]]:
    """(heading, preamble Markdown, decisions) of the lead-decision sections of requirements section 7
    (7.16, 7.17, 7.18 ..., see :func:`lead_section_numbers`)."""
    text = _read(path)
    out = []
    decisions = parse_decisions(path)
    for num in lead_section_numbers(text):
        sec = _section(text, r"^" + re.escape(num) + r"\b")
        if not sec:
            continue
        lines = sec.splitlines()
        heading = re.sub(r"^#+\s*", "", lines[0]).strip()
        pre = []
        for line in lines[1:]:
            if line.lstrip().startswith("|"):
                break
            pre.append(line)
        rows = [d for d in decisions.values() if d.section == num]
        out.append((heading, "\n".join(pre).strip(), rows))
    return out


def _heading_titles(text: str) -> Dict[str, str]:
    """Map section numbers to heading titles: '2.4' / 'C.1' / 'A' / '0' -> full heading text."""
    out: Dict[str, str] = {}
    for line in text.splitlines():
        m = re.match(r"^#{2,4}\s+((?:[A-Z]|\d+)(?:\.\d+)*)\.?\s+(.*)$", line)
        if m:
            out.setdefault(m.group(1), f"{m.group(1)} {m.group(2).strip()}")
    return out


def _id_rows(text: str, pattern: str) -> Dict[str, str]:
    """First-column ids of table rows matching ``pattern`` -> 'id: second/third column'."""
    out: Dict[str, str] = {}
    for header, rows in parse_tables(text):
        for r in rows:
            if r and re.fullmatch(pattern, r[0].strip()):
                rid = r[0].strip()
                desc = r[2] if len(r) > 2 and header and header[1].strip().lower() == "script" else (
                    r[1] if len(r) > 1 else "")
                out.setdefault(rid, desc)
    return out


@dataclass
class SourceIndex:
    """Titles of the R1/R2 sections and items cited by the ``source`` strings of the VPs."""
    r1: Dict[str, str] = field(default_factory=dict)
    r2: Dict[str, str] = field(default_factory=dict)
    r1_checks: Dict[str, str] = field(default_factory=dict)
    r2_items: Dict[str, str] = field(default_factory=dict)
    decisions: Dict[str, Decision] = field(default_factory=dict)

    @classmethod
    def load(cls, r1: Path = R1_PATH, r2: Path = R2_PATH, req: Path = REQUIREMENTS_PATH) -> "SourceIndex":
        t1, t2 = _read(r1), _read(r2)
        return cls(r1=_heading_titles(t1), r2=_heading_titles(t2), r1_checks=_id_rows(t1, r"V\d+"),
                   r2_items=_id_rows(t2, r"[A-L]\d+"), decisions=parse_decisions(req))

    def cite(self, source: str) -> List[Tuple[str, str]]:
        """(label, description) for every R1/R2 section and decision cited in ``source``."""
        out: List[Tuple[str, str]] = []
        seen = set()

        def add(label: str, desc: str) -> None:
            if label not in seen:
                seen.add(label)
                out.append((label, desc))

        for part in re.split(r";", source):
            # walk the tokens; an R1 / R2 prefix applies to the ids that follow it
            prefix = None
            for tok in re.findall(r"R[12]\b|spec\b|requirements\b|D-[A-Z0-9]+-\d+|§?\s*[A-L]\.\d+|"
                                  r"§?\s*\d+(?:\.\d+)*|V\d+|[A-L]\d+|section\s+\d+", part):
                tok = tok.strip()
                if tok in ("R1", "R2"):
                    prefix = tok
                    continue
                if tok in ("spec", "requirements"):
                    prefix = None
                    continue
                if tok.startswith("D-"):
                    d = self.decisions.get(tok)
                    add(tok, _plain(d.subject) if d else "decision of requirements section 7")
                    continue
                if prefix is None:
                    continue
                key = re.sub(r"^(§|section)\s*", "", tok).strip()
                if prefix == "R2":
                    if re.fullmatch(r"[A-L]\d+", key):
                        desc = self.r2_items.get(key, "")
                        sec = self.r2.get(key[0], "")
                        add(f"R2 {key}", "; ".join(x for x in (_plain(desc), _plain(sec)) if x))
                    elif key in self.r2:
                        add(f"R2 §{key}" if key.isdigit() else f"R2 {key}", _plain(self.r2[key]))
                else:
                    if re.fullmatch(r"V\d+", key):
                        add(f"R1 {key}", _plain(self.r1_checks.get(key, "numerical check of R1 section 9")))
                    elif key in self.r1:
                        add(f"R1 §{key}", _plain(self.r1[key]))
        return out


def _plain(md: str) -> str:
    """Strip Markdown emphasis and LaTeX dollars from a heading/cell for use in a table."""
    s = re.sub(r"\*\*([^*]+)\*\*", r"\1", md or "")
    s = s.replace("$", "")
    return " ".join(s.split())


# --------------------------------------------------------------------------------------
# VP descriptions
# --------------------------------------------------------------------------------------
def _module_paragraph(doc: str, vp_id: str) -> str:
    """The paragraph of a module docstring that describes ``vp_id`` (bullet or leading id)."""
    if not doc:
        return ""
    lines = doc.splitlines()
    pat = re.compile(r"^\s*(?:[*-]\s+)?" + re.escape(vp_id) + r"(?![\w,])\s*(?::|--|-(?=\s))?\s*(.*)$")
    start = None
    first = ""
    other = re.compile(r"^\s*(?:[*-]\s+)?VP-[0-9A-Za-z]+\b")
    for i, line in enumerate(lines):
        m = pat.match(line)
        if m and not re.search(r",\s*VP-[0-9A-Za-z]+", m.group(1)[:60]):   # not a list "VP-25 x, VP-28 y"
            start = i
            first = m.group(1)
            break
    if start is None:
        # a bullet item that names the problem in bold or as a label: "* defines **VP-E1**: ..."
        idre = re.escape(vp_id)
        pat2 = re.compile(r"^\s*[*-]\s+.*?(?:\*\*" + idre + r"\*\*\s*:?|(?<![\w-])" + idre + r":)\s*(.*)$")
        for i, line in enumerate(lines):
            m = pat2.match(line)
            if m:
                start = i
                first = m.group(1)
                other = re.compile(r"^\s*[*-]\s+")
                break
    if start is None:
        return ""
    para = [first] if first else []
    for line in lines[start + 1:]:
        if not line.strip() or other.match(line):
            break
        para.append(line.strip())
    text = " ".join(p for p in para if p).strip()
    return text


def describe_problem(p: Problem) -> str:
    """Markdown description of a VP: its function docstring, else its paragraph in the docstring of
    its module."""
    doc = inspect.getdoc(p.func) or ""
    if doc.strip():
        return rst_to_md(doc)
    mod = sys.modules.get(getattr(p.func, "__module__", ""), None)
    mdoc = inspect.getdoc(mod) if mod is not None else ""
    para = _module_paragraph(mdoc or "", p.id)
    return rst_to_md(para) if para else ""


# --------------------------------------------------------------------------------------
# Coverage matrix
# --------------------------------------------------------------------------------------
def module_columns(modules: Sequence[str]) -> List[str]:
    """Matrix columns of a VP's module list (non-module entries go to the UI/tools column)."""
    cols = []
    for m in modules:
        col = m.upper() if m.upper() in MODULE_COLUMNS else OTHER_COLUMN
        if col not in cols:
            cols.append(col)
    return cols


def coverage(runs: Sequence[VPRun]) -> Dict[str, Dict[str, int]]:
    """Per matrix column: number of VPs and how many passed / failed / errored / were skipped."""
    out = {c: {"vps": 0, STATUS_PASS: 0, STATUS_FAIL: 0, STATUS_ERROR: 0, STATUS_SKIPPED: 0}
           for c in MODULE_COLUMNS + (OTHER_COLUMN,)}
    for r in runs:
        for c in module_columns(r.modules):
            out[c]["vps"] += 1
            out[c][r.status] += 1
    return out


# --------------------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------------------
#: colours of the summary figure (reference palette of the dataviz method: one series in blue,
#: text in the primary/secondary inks, the tolerance line neutral)
_FIG = {"surface": "#fcfcfb", "series": "#2a78d6", "text": "#0b0b0b", "muted": "#52514e",
        "grid": "#e4e3df", "rule": "#52514e"}


def margin_figure(runs: Sequence[VPRun], path: Path) -> Optional[Path]:
    """Dot plot of the largest error/tolerance ratio of every VP (log scale).

    A dot left of the tolerance line passes.  Ratios below ``1e-12`` (including exact agreement)
    are drawn as open circles on the left edge.  VPs without numerical criteria, skipped VPs and VPs
    that raised are not drawn (they are listed in the summary table).  Returns the PNG path or None
    when matplotlib is not available or nothing can be drawn."""
    data = [(r.id, r.margin(), r.status) for r in runs if r.result is not None and r.margin() is not None]
    if not data:
        return None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:       # pragma: no cover - matplotlib is a dependency of the package
        return None
    lo, hi = 1e-12, 1e3
    ids = [d[0] for d in data][::-1]
    raw = [d[1] for d in data][::-1]
    stat = [d[2] for d in data][::-1]
    vals = [hi if not math.isfinite(v) else min(max(v, lo), hi) for v in raw]
    y = list(range(len(ids)))
    h = max(3.0, 0.2 * len(ids) + 1.3)
    fig, ax = plt.subplots(figsize=(7.0, h), dpi=100)
    fig.patch.set_facecolor(_FIG["surface"])
    ax.set_facecolor(_FIG["surface"])
    for yi in y:                                   # thin guide line per problem
        ax.plot([lo, hi], [yi, yi], color=_FIG["grid"], linewidth=0.5, zorder=1)
    filled = [i for i, v in enumerate(raw) if not (math.isfinite(v) and v < lo)]
    hollow = [i for i, v in enumerate(raw) if math.isfinite(v) and v < lo]
    ax.scatter([vals[i] for i in filled], [y[i] for i in filled], s=36, color=_FIG["series"], zorder=3,
               edgecolors=_FIG["surface"], linewidths=1.0)
    ax.scatter([vals[i] for i in hollow], [y[i] for i in hollow], s=36, facecolors=_FIG["surface"],
               edgecolors=_FIG["series"], linewidths=1.2, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(lo / 3, hi)
    ax.set_ylim(-1, len(ids))
    ax.axvline(1.0, color=_FIG["rule"], linewidth=1.0, linestyle="--", zorder=2)
    ax.text(1.0, len(ids) - 0.4, " tolerance", color=_FIG["muted"], fontsize=8, va="bottom", ha="left")
    for yi, v, s_ in zip(y, vals, stat):
        if s_ != STATUS_PASS:
            ax.text(v * 1.6, yi, s_, va="center", ha="left", fontsize=7, color=_FIG["text"], fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(ids, fontsize=7, color=_FIG["text"])
    ax.tick_params(axis="x", labelsize=8, colors=_FIG["muted"])
    ax.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(_FIG["grid"])
    ax.set_xlabel("largest error / tolerance over the numerical criteria (log scale; < 1 passes; "
                  "open circle: below 1e-12)", fontsize=8, color=_FIG["muted"])
    ax.set_title("Verification margins per problem", fontsize=10, color=_FIG["text"], loc="left")
    fig.tight_layout()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, facecolor=_FIG["surface"])
    plt.close(fig)
    return path


# --------------------------------------------------------------------------------------
# The verification manual
# --------------------------------------------------------------------------------------
def _environment() -> List[Tuple[str, str]]:
    from .. import PRODUCT, __version__
    env = [("Product", f"{PRODUCT} {__version__}"), ("Python", platform.python_version()),
           ("Platform", f"{platform.system()} {platform.release()} ({platform.machine()})")]
    for mod in ("numpy", "scipy", "matplotlib"):
        try:
            m = __import__(mod)
            env.append((mod, getattr(m, "__version__", "?")))
        except Exception:   # pragma: no cover
            env.append((mod, "not installed"))
    return env


def _rel(path: Path, base: Path) -> str:
    """Relative POSIX path from directory ``base`` to ``path`` (for Markdown links)."""
    try:
        import os
        return Path(os.path.relpath(Path(path).resolve(), Path(base).resolve())).as_posix()
    except ValueError:      # pragma: no cover - different drives on Windows
        return Path(path).as_posix()


def _status_line(r: VPRun) -> str:
    if r.status == STATUS_SKIPPED:
        return f"SKIPPED ({r.skipped})"
    if r.status == STATUS_ERROR:
        return "ERROR (the problem raised an exception; see below)"
    crit = r.criteria
    npass = sum(1 for c in crit if c.passed)
    word = "PASS" if r.status == STATUS_PASS else "FAIL"
    return f"{word}: {npass} of {len(crit)} criteria met, {len(r.informative)} informative, {r.elapsed:.1f} s"


def _decisions_for(vp_id: str, decisions: Sequence[Decision]) -> List[Decision]:
    pat = re.compile(r"(?<![\w-])" + re.escape(vp_id) + r"(?![\w])")
    return [d for d in decisions if pat.search(d.evidence) or pat.search(d.subject)]


def manual_markdown(runs: Sequence[VPRun], *, title: str = "SASSI-EDU Verification Manual",
                    command: str = "python -m sassi.verify.report", figure: Optional[Path] = None,
                    out_dir: Optional[Path] = None, requirements: Path = REQUIREMENTS_PATH,
                    sources: Optional[SourceIndex] = None, load_errors: Optional[Dict[str, str]] = None,
                    generated: Optional[str] = None, partial: bool = False) -> str:
    """Compose the verification manual for ``runs`` (see the module docstring for the content)."""
    out_dir = Path(out_dir) if out_dir else DEFAULT_MANUAL.parent
    req_rows = requirement_rows(requirements)
    sources = sources or SourceIndex.load(req=requirements)
    lead = lead_decision_sections(requirements)
    lead_rows = [d for _, _, rows in lead for d in rows]
    lead_txt = lead_sections_text([h.split()[0] for h, _, _ in lead] or lead_section_numbers(""))
    load_errors = dict(LOAD_ERRORS if load_errors is None else load_errors)
    generated = generated or datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    n = len(runs)
    by_status = {s: [r for r in runs if r.status == s] for s in (STATUS_PASS, STATUS_FAIL, STATUS_ERROR,
                                                                 STATUS_SKIPPED)}
    n_crit = sum(len(r.criteria) for r in runs)
    n_crit_pass = sum(1 for r in runs for c in r.criteria if c.passed)
    n_info = sum(len(r.informative) for r in runs)
    total_time = sum(r.elapsed for r in runs)
    L: List[str] = []
    add = L.append

    add(f"# {title}")
    add("")
    add(f"> Generated by `{command}` on {generated}. This file is produced by "
        "`sassi/verify/report.py` from the verification problems registered in `sassi.verify`; "
        "do not edit it by hand, run the generator again.")
    if partial:
        add(">")
        add("> **Partial report:** only a selection of the verification problems was run.")
    add("")
    add("SASSI-EDU is an educational re-implementation of the SASSI flexible-volume method with the "
        "module and command nomenclature of the ACS SASSI V3 manual. It is not affiliated with, endorsed "
        "by, or a substitute for ACS SASSI, and it is not qualified for licensing or design work. This "
        "manual is the equivalent of the vendor's verification and validation (V&V) manual: it shows, "
        "problem by problem, what the program computes and how close that is to an independent "
        "reference.")
    add("")
    add("## Contents")
    add("")
    for k, t in enumerate(["How the verification works", "Summary", "Module coverage matrix",
                           "Lead decisions and superseded references", "Verification problems",
                           "Reference sources", "Environment and run log"], start=1):
        add(f"{k}. [{t}](#{slug(f'{k}. {t}')})")
    add("")

    # 1 ----------------------------------------------------------------------------------
    add("## 1. How the verification works")
    add("")
    add("Each **verification problem (VP)** is a Python function in `sassi/verify/problems/` registered "
        "with `@problem`. It builds a model (through the command language, the module decks or the "
        "element library), runs the modules it verifies and compares computed quantities with "
        "reference values. The reference values and tolerances come from "
        "[R2 benchmarks](../spec/R2_benchmarks.md) (published **P**, exact **E**, approximate **A** "
        "or derived **D** solutions), from the theory checks of "
        "[R1](../spec/R1_sassi_theory.md) section 9, and from the verification plan of "
        "[requirements section 6](../spec/00_requirements.md).")
    add("")
    add("Every comparison is one row of the VP's table:")
    add("")
    add("* **rel**: relative error `|computed - reference| / |reference|` against a relative tolerance;")
    add("* **abs**: absolute error `|computed - reference|` against an absolute tolerance;")
    add("* **condition**: a logical requirement (for example *symmetric*, *monotone*, *file written*);")
    add("* **informative** (`info`): a comparison with a published *approximate* reference that the "
        f"lead decisions of requirements sections {lead_txt} have superseded by a rigorous one "
        "(section 4). It is printed with its deviation but is not a pass/fail criterion; the "
        "superseding comparison is always an ordinary criterion of the same VP.")
    add("")
    add("A VP passes when every criterion passes. Tolerances are never relaxed to make a VP pass "
        "(requirements section 6.1, ARCHITECTURE section 10): a VP that cannot meet its tolerance is "
        "reported as FAIL with the observed numbers.")
    add("")
    add("**Reproduce this report.** From the project root:")
    add("")
    add("```bash")
    add(".venv/bin/python -m sassi.verify.report               # all VPs -> docs/verification/VERIFICATION_MANUAL.md")
    add(".venv/bin/python -m sassi.verify.report --skip-slow   # leave out the slow VPs")
    add(".venv/bin/python -m pytest -q tests/verification      # the same VPs as pytest cases")
    add("```")
    add("")
    add("In the console or the GUI command line, `VERIFY,ALL` runs the VPs and prints the same rows, and "
        "`VERIFYREPORT` writes this report into the working directory.")
    add("")

    # 2 ----------------------------------------------------------------------------------
    add("## 2. Summary")
    add("")
    add(md_table(["Item", "Count"], [
        ["verification problems in the report", n],
        ["passed", len(by_status[STATUS_PASS])],
        ["failed", len(by_status[STATUS_FAIL])],
        ["errors (the problem raised an exception)", len(by_status[STATUS_ERROR])],
        ["skipped", len(by_status[STATUS_SKIPPED])],
        ["criteria checked (met)", f"{n_crit} ({n_crit_pass})"],
        ["informative comparisons", n_info],
        ["total run time", f"{total_time:.0f} s"],
    ], align="lr"))
    add("")
    if by_status[STATUS_FAIL] or by_status[STATUS_ERROR]:
        bad = ", ".join(f"[{r.id}](#{slug(r.id)})" for r in by_status[STATUS_FAIL] + by_status[STATUS_ERROR])
        add(f"**Not passing in this run:** {bad}. Their sections give the observed values and the reasons.")
        add("")
    if load_errors:
        add("**Problem modules that could not be imported** (their VPs are missing from this report):")
        add("")
        for mod, err in sorted(load_errors.items()):
            add(f"* `{mod}`: {md_escape(err.strip().splitlines()[-1] if err.strip() else err)}")
        add("")
    rows = []
    for r in runs:
        crit = r.criteria
        npass = sum(1 for c in crit if c.passed)
        margin = r.margin()
        rows.append([f"[{r.id}](#{slug(r.id)})", md_escape(r.title), r.tier, ", ".join(r.modules),
                     f"{npass}/{len(crit)}" if r.result is not None else "-",
                     str(len(r.informative)) if r.result is not None else "-",
                     fmt_small(margin) if margin is not None else "-",
                     f"{r.elapsed:.1f}" if not r.skipped else "-", r.status])
    add(md_table(["VP", "Title", "Tier", "Modules", "Criteria met", "Info", "Max err/tol", "Time (s)",
                  "Status"], rows, align="lllllrrrl", escape=False))
    add("")
    add("*Max err/tol* is the largest ratio of error to tolerance over the numerical criteria of the "
        "VP: below 1 the VP passes with that margin; `-` means the VP has only logical conditions.")
    add("")
    if figure is not None:
        add(f"![Largest error/tolerance ratio of each verification problem]({_rel(figure, out_dir)})")
        add("")
        add("*Figure 1. Largest error/tolerance ratio of every VP with numerical criteria (log scale). "
            "Dots left of the dashed tolerance line pass; open circles mark agreement better than 1e-12 "
            "of the tolerance; problems that do not pass are labelled.*")
        add("")

    # 3 ----------------------------------------------------------------------------------
    add("## 3. Module coverage matrix")
    add("")
    add("Which analysis modules each VP exercises (requirements section 6.5). The letter is the VP's "
        f"status: **P** pass, **F** fail, **E** error, **S** skipped. *{OTHER_COLUMN}* collects the "
        "interpreter, CHECK/AFWRITE, generation, section-cut, line-maths and converter tools.")
    add("")
    cols = list(MODULE_COLUMNS) + [OTHER_COLUMN]
    mrows = []
    for r in runs:
        have = set(module_columns(r.modules))
        mrows.append([f"[{r.id}](#{slug(r.id)})"] + [STATUS_LETTER[r.status] if c in have else "" for c in cols])
    add(md_table(["VP"] + cols, mrows, align="l" + "c" * len(cols), escape=False))
    add("")
    cov = coverage(runs)
    add(md_table(["Module", "VPs", "Passed", "Failed", "Errors", "Skipped"],
                 [[c, cov[c]["vps"], cov[c][STATUS_PASS], cov[c][STATUS_FAIL], cov[c][STATUS_ERROR],
                   cov[c][STATUS_SKIPPED]] for c in cols], align="lrrrrr"))
    add("")

    # 4 ----------------------------------------------------------------------------------
    add("## 4. Lead decisions and superseded references")
    add("")
    add("Some published references turned out to be approximate, internally inconsistent or affected "
        "by a documented defect when the implementation was verified. The project lead then decided, "
        "with evidence, how each affected VP is judged. These decisions are binding (requirements "
        f"sections {lead_txt}) and are reproduced here so that every informative comparison of "
        "this manual can be traced to its reason.")
    add("")
    if not lead:
        add("*The requirements document was not found; the decisions are listed in "
            "`docs/spec/00_requirements.md` section 7 (Lead decisions).*")
        add("")
    for heading, pre, drows in lead:
        num, _, rest = heading.partition(" ")
        add(f"### Requirements {num}: {rest}" if re.match(r"^\d+\.\d+$", num) else f"### {heading}")
        add("")
        if pre:
            add(pre)
            add("")
        add(md_table(["ID", "Subject", "Decision", "Evidence"],
                     [[d.id, d.subject, d.decision, d.evidence] for d in drows], escape=False))
        add("")
    info_rows = []
    for r in runs:
        for c in r.informative:
            ds = ", ".join(d.id for d in _decisions_for(r.id, lead_rows)) or "-"
            info_rows.append([f"[{r.id}](#{slug(r.id)})", md_escape(c.quantity), fmt_num(c.computed),
                              fmt_num(c.reference), fmt_small(c.error), ds])
    add("### Informative comparisons in this run")
    add("")
    if info_rows:
        add(md_table(["VP", "Quantity", "Computed", "Superseded reference", "Rel. difference", "Decision"],
                     info_rows, align="llrrrl", escape=False))
    else:
        add("No informative comparisons were recorded in this run.")
    add("")

    # 5 ----------------------------------------------------------------------------------
    add("## 5. Verification problems")
    add("")
    add("One section per problem, in natural order. *Purpose* and *Setup* are quoted from the "
        "verification plan (requirements sections 6.3 and 6.4) when the VP is listed there; "
        "*Model and method* is the description written with the VP.")
    add("")
    for r in runs:
        add(f"### {r.id}")
        add("")
        add(f"**{md_escape(r.title)}**")
        add("")
        row, relation = requirement_row_for(r.id, req_rows)
        info = [["Status", md_escape(_status_line(r))], ["Tier", r.tier],
                ["Modules verified", md_escape(", ".join(r.modules) or "-")]]
        if r.source:
            info.append(["Reference source", md_escape(r.source)])
        if row is not None:
            sec = row.get("section", "")
            label = f"requirements {sec}" + (f" (extends {row.get('ID', '').strip()})" if relation else "")
            info.append(["Verification plan", label])
            if row.get("Type"):
                info.append(["Reference type", row["Type"]])
            tol = row.get("Tolerance")
            if tol:
                info.append(["Planned tolerance", tol])
        decs = _decisions_for(r.id, lead_rows)
        if decs:
            info.append(["Lead decisions", ", ".join(f"[{d.id}](#{slug('4. Lead decisions and superseded references')})"
                                                      for d in decs)])
        add(md_table(["Item", "Value"], info, escape=False))
        add("")
        if row is not None:
            purpose = row.get("Problem") or row.get("Test") or ""
            expected = row.get("Expected reference") or row.get("Expected") or ""
            setup = row.get("Setup", "")
            if purpose:
                add(f"**Purpose.** {purpose}")
                add("")
            if setup:
                add(f"**Setup (plan).** {setup}")
                add("")
            if expected:
                add(f"**Expected (plan).** {expected}")
                add("")
        if r.description:
            add("**Model and method.**")
            add("")
            add(r.description)
            add("")
        cites = sources.cite(r.source) if r.source else []
        if cites:
            add("**References cited.**")
            add("")
            for lab, desc in cites:
                add(f"* {md_escape(lab)}: {md_escape(desc)}" if desc else f"* {md_escape(lab)}")
            add("")
        if r.status == STATUS_SKIPPED:
            add(f"*Not run: {md_escape(r.skipped)}.*")
            add("")
            continue
        if r.status == STATUS_ERROR:
            tail = "\n".join(r.error.strip().splitlines()[-12:])
            add("**The problem raised an exception** (last lines of the traceback):")
            add("")
            add("```text")
            add(tail)
            add("```")
            add("")
            continue
        add("**Results.**")
        add("")
        checks = r.checks
        with_note = any(c.note for c in checks)
        headers = ["#", "Quantity", "Computed", "Reference", "Error", "Tolerance", "Result"]
        if with_note:
            headers.append("Note")
        add(md_table(headers, [check_row(k, c, with_note) for k, c in enumerate(checks, start=1)],
                     align="rlrrrrl" + ("l" if with_note else "")))
        add("")
        if r.notes:
            add("**Notes.**")
            add("")
            for note in r.notes:
                add(f"* {md_escape(' '.join(note.split()))}")
            add("")
        for fig in r.figures:
            add(f"![{md_escape(r.id)} figure]({_rel(Path(fig), out_dir)})")
            add("")

    # 6 ----------------------------------------------------------------------------------
    add("## 6. Reference sources")
    add("")
    add("The VPs cite their references by section: **R2** is the benchmark collection "
        "[`docs/spec/R2_benchmarks.md`](../spec/R2_benchmarks.md) (its section L lists the publications "
        "with links), **R1** the SASSI theory note [`docs/spec/R1_sassi_theory.md`](../spec/R1_sassi_theory.md) "
        "(section 11 lists its sources; V1 to V11 are its numerical checks), **D-xxx-nn** the decisions "
        "of [requirements section 7](../spec/00_requirements.md), and \"spec NN\" the manual-derived "
        "specifications in `docs/spec/`. The independent boundary-element references of the impedance "
        "problems are described in [impedance_study.md](impedance_study.md).")
    add("")
    cited: Dict[str, Tuple[str, List[str]]] = {}
    for r in runs:
        for lab, desc in (sources.cite(r.source) if r.source else []):
            if lab.startswith("D-"):
                continue
            entry = cited.setdefault(lab, (desc, []))
            entry[1].append(r.id)
    if cited:
        def _ckey(lab: str):
            return (lab[:2], [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", lab[3:])])
        add(md_table(["Section", "Title", "Cited by"],
                     [[md_escape(lab), md_escape(desc), ", ".join(f"[{v}](#{slug(v)})" for v in vps)]
                      for lab, (desc, vps) in sorted(cited.items(), key=lambda kv: _ckey(kv[0]))],
                     escape=False))
        add("")

    # 7 ----------------------------------------------------------------------------------
    add("## 7. Environment and run log")
    add("")
    env = _environment() + [("Generated", generated), ("Command", f"`{command}`"),
                            ("Total VP run time", f"{total_time:.1f} s")]
    add(md_table(["Item", "Value"], env, escape=False))
    add("")
    slowest = sorted((r for r in runs if not r.skipped), key=lambda r: -r.elapsed)[:10]
    if slowest:
        add("Slowest problems:")
        add("")
        add(md_table(["VP", "Time (s)", "Registered as slow"],
                     [[f"[{r.id}](#{slug(r.id)})", f"{r.elapsed:.1f}", "yes" if r.slow else "no"]
                      for r in slowest], align="lrl", escape=False))
        add("")
    return "\n".join(L).rstrip() + "\n"


def generate_manual(out: Path = DEFAULT_MANUAL, *, only: Optional[Sequence[str]] = None,
                    tiers: Optional[Sequence[str]] = None, skip_slow: bool = False, figures: bool = False,
                    figure_dir: Optional[Path] = None, keep: bool = False,
                    progress: Optional[Callable[[int, int, VPRun], None]] = _default_progress,
                    command: Optional[str] = None, title: str = "SASSI-EDU Verification Manual",
                    registry: Optional[Dict[str, Problem]] = None) -> Tuple[Path, List[VPRun]]:
    """Run the selected VPs and write the manual to ``out``.  Returns (path, runs)."""
    reg = registry if registry is not None else load_all()
    selection = select_problems(reg, only=only, tiers=tiers, skip_slow=skip_slow)
    out = Path(out)
    figdir = Path(figure_dir) if figure_dir else out.parent / "figures"
    runs = run_vps(selection, keep=keep, figdir=figdir if figures else None, progress=progress)
    fig_path = None
    if figures:
        fig_path = margin_figure(runs, figdir / "vp_margins.png")
    text = manual_markdown(runs, title=title, command=command or "python -m sassi.verify.report",
                           figure=fig_path, out_dir=out.parent, partial=bool(only or tiers))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return out, runs


# --------------------------------------------------------------------------------------
# Command reference
# --------------------------------------------------------------------------------------
_CATEGORY_TITLES = {
    "3.4.A": "Session, files, models and global options", "3.4.B": "Frequency sets",
    "3.4.C": "Nodes and coordinate systems", "3.4.D": "Groups, elements and properties",
    "3.4.E": "Loads and masses", "3.4.F": "Module analysis options", "3.4.G": "Module run commands",
    "3.4.H": "Model checking", "3.4.I": "Model conditioning and generation",
    "3.4.J": "Cuts, submodels and section calculations", "3.4.K": "File conversion",
    "3.4.L": "Plotting and line mathematics", "3.4.M": "Programming (variables, loops, macros)",
    "3.4.N": "Water modelling", "3.4.O": "Option NON and nonlinear soil", "3.4.P": "Binary databases",
    "3.4.Q": "Thick shell", "3.4.R": "Extension commands of SASSI-EDU",
}
OTHER_CATEGORY = "Other commands registered by SASSI-EDU"


def catalogue_categories() -> Dict[str, str]:
    """Command name -> requirements category ('3.4.A' ...) from the registry catalogue."""
    from ..prep import registry
    text = getattr(registry, "_CATALOGUE_TEXT", "")
    out: Dict[str, str] = {}
    cat = ""
    for line in text.splitlines():
        m = re.match(r"^#\s*-+\s*(3\.4\.[A-Z])\b", line.strip())
        if m:
            cat = m.group(1)
            continue
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        out[s.split()[0].upper()] = cat
    return out


def requirement_syntax(path: Path = REQUIREMENTS_PATH) -> Dict[str, Tuple[str, str]]:
    """Command name -> (syntax, meaning) from the tables of requirements section 3.4."""
    text = _read(path)
    sec = _section(text, r"^3\.4\b")
    out: Dict[str, Tuple[str, str]] = {}
    for header, rows in parse_tables(sec):
        hl = [h.strip().lower() for h in header]
        if not hl or hl[0] != "command":
            continue
        i_syn = next((k for k, h in enumerate(hl) if h.startswith("syntax")), None)
        i_mean = next((k for k, h in enumerate(hl) if h in ("meaning", "module", "reports", "carries")), None)
        for r in rows:
            names = [n.strip().strip("*†").strip() for n in re.split(r"[/,]", r[0])]
            names = [re.sub(r"\s*\(.*\)$", "", n).strip().upper() for n in names if n.strip()]
            syn_cell = r[i_syn] if i_syn is not None and i_syn < len(r) else ""
            spans = re.findall(r"`([^`]+)`", syn_cell)
            meaning = r[i_mean] if i_mean is not None and i_mean < len(r) else ""
            for name in names:
                if not re.fullmatch(r"[A-Z0-9_]+", name):
                    continue
                syn = ""
                for s in spans:
                    head = re.split(r"[,\s]", s.strip(), maxsplit=1)[0].upper()
                    if head == name:
                        syn = s.strip()
                        break
                if not syn and spans:
                    s = spans[0].strip()
                    head = re.split(r"[,\s]", s, maxsplit=1)[0]
                    if head.upper() == "RUN<MOD>" and name.startswith("RUN"):
                        syn = name + s[len(head):]
                    elif len(names) > 1 and head.upper() in names:
                        syn = name + s[len(head):]
                out.setdefault(name, (syn, meaning))
    return out


def _doc_syntax(name: str, doc: str) -> Tuple[str, str]:
    """(syntax, meaning) from a handler docstring 'NAME,<a>,...: meaning' (first paragraph)."""
    if not doc:
        return "", ""
    para = " ".join(doc.strip().split("\n\n")[0].split())
    head = re.split(r"[,:\s(]", para, maxsplit=1)[0].upper()
    if head != name:
        return "", para
    m = re.search(r":\s", para)
    if m:
        syn, meaning = para[:m.start()], para[m.end():]
    else:
        syn, meaning = para.rstrip("."), ""
    syn = re.sub(r"\s*\(alias(?:es)?\b[^)]*\)", "", syn)
    return syn.strip(), meaning.strip()


_ABBREVIATIONS = {"e.g", "i.e", "etc", "vs", "eq", "eqs", "fig", "figs", "sec", "no", "approx", "cf", "incl",
                  "resp", "ref", "al", "max", "min", "rel", "abs", "deg"}


def _first_sentence(text: str, limit: int = 220) -> str:
    """First sentence of ``text`` (abbreviations such as 'e.g.' do not end it), capitalised and cut
    at ``limit`` characters."""
    text = " ".join(text.split())
    end = len(text)
    for m in re.finditer(r"(?<=[a-z0-9)\]])\.\s+(?=[A-Z(])", text):
        word = re.split(r"[\s(]", text[:m.start()])[-1].lower()
        if word in _ABBREVIATIONS or word.endswith("e.g") or word.endswith("i.e"):
            continue
        end = m.start() + 1
        break
    s = text[:end]
    if len(s) > limit:
        s = s[:limit].rsplit(" ", 1)[0] + " ..."
    return s[:1].upper() + s[1:]


@dataclass
class CommandInfo:
    name: str
    abbrev: str
    tier: str
    cls: str
    status: str            # 'yes', 'stored only', 'no'
    syntax: str
    meaning: str
    category: str
    module: str = ""


def command_infos() -> List[CommandInfo]:
    """Every registered command with its syntax, abbreviation, tier, implementation status and
    meaning (handler docstring first, requirements section 3.4 as the fallback)."""
    from ..prep import registry
    registry.load_commands()
    cats = catalogue_categories()
    req = requirement_syntax()
    out: List[CommandInfo] = []
    for spec in registry.all_commands():
        name = spec.name
        doc = "" if spec.placeholder else (inspect.getdoc(spec.handler) or "")
        dsyn, dmean = _doc_syntax(name, doc)
        rsyn, rmean = req.get(name, ("", ""))
        syntax = re.sub(r",\s+", ",", dsyn or rsyn or name)
        meaning = dmean or _plain(rmean) or spec.note or ""
        if spec.placeholder:
            if spec.cls in ("record", "record_keep", "indexed", "string"):
                status = "stored only"
                tail = f"parsed and stored for WRITE; no computation in this build (tier {spec.tier})"
            else:
                status = "no"
                tail = spec.note or f"prints '{name} is not available in this build (tier {spec.tier})'"
            meaning = (meaning + " -- " if meaning else "") + tail
        else:
            status = "yes"
        out.append(CommandInfo(name=name, abbrev=", ".join(spec.abbrev), tier=spec.tier, cls=spec.cls,
                               status=status, syntax=syntax, meaning=_first_sentence(meaning),
                               category=cats.get(name, ""), module=spec.module))
    return out


def command_reference_markdown(infos: Optional[Sequence[CommandInfo]] = None,
                               generated: Optional[str] = None) -> str:
    """The command reference: one table per command category, plus a summary by tier."""
    infos = list(infos) if infos is not None else command_infos()
    generated = generated or datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d")
    L: List[str] = []
    add = L.append
    add("# SASSI-EDU Command Reference")
    add("")
    add(f"> Generated on {generated} by `python -m sassi.verify.report --commands-only` from the "
        "command catalogue of the interpreter (`sassi/prep/registry.py`) and the docstrings of the "
        "command handlers (`sassi/prep/commands/`). Do not edit by hand.")
    add("")
    add("Every command of the ACS SASSI V3 manual and of the SASSI-EDU dialect is listed with its "
        "syntax, its documented abbreviation or alias, its priority tier and whether this build "
        "implements it. The rules of the command language (comma-separated fields, blank fields take "
        "defaults, `*` comment lines, abbreviations, variables and loops) are explained in the "
        "[User Guide](../user/USER_GUIDE.md#4-the-command-language); the normative detail of every "
        "command is in [requirements section 3](../spec/00_requirements.md) and the manual-derived "
        "specifications `docs/spec/07` to `11`.")
    add("")
    add("**Syntax notation.** `<arg>` is a required argument, `[arg]` an optional one; a blank field "
        "takes the documented default. Names are case-insensitive. Only the full name or the listed "
        "abbreviation is accepted (no prefix matching).")
    add("")
    add("**Tiers.** P0 = core linear SSI chain; P1 = important engineering features; P2 = advanced or "
        "cosmetic (requirements section 0.1).")
    add("")
    add("**Implemented.** *yes*: the command works in this build. *stored only*: the command is parsed "
        "and stored (WRITE writes it back, so `.pre` files of the original program load without loss) "
        "but the feature behind it is not available. *no*: the command prints "
        "`<CMD> is not available in this build (tier Pn)` or the manual's own message.")
    add("")
    order = list(_CATEGORY_TITLES)
    groups: Dict[str, List[CommandInfo]] = {}
    for ci in infos:
        groups.setdefault(ci.category if ci.category in _CATEGORY_TITLES else "", []).append(ci)
    add("## Summary")
    add("")
    srows = []
    cat_list = [c for c in order if c in groups] + ([""] if "" in groups else [])
    for c in cat_list:
        lst = groups[c]
        title = f"{c} {_CATEGORY_TITLES[c]}" if c else OTHER_CATEGORY
        srows.append([f"[{title}](#{slug(title)})", len(lst), sum(1 for x in lst if x.status == "yes"),
                      sum(1 for x in lst if x.status == "stored only"), sum(1 for x in lst if x.status == "no")])
    srows.append(["**Total**", len(infos), sum(1 for x in infos if x.status == "yes"),
                  sum(1 for x in infos if x.status == "stored only"), sum(1 for x in infos if x.status == "no")])
    add(md_table(["Category", "Commands", "Implemented", "Stored only", "Not available"], srows,
                 align="lrrrr", escape=False))
    add("")
    tiers = sorted({x.tier for x in infos})
    add(md_table(["Tier", "Commands", "Implemented"],
                 [[t, sum(1 for x in infos if x.tier == t), sum(1 for x in infos if x.tier == t and x.status == "yes")]
                  for t in tiers], align="lrr"))
    add("")
    for c in cat_list:
        lst = groups[c]
        title = f"{c} {_CATEGORY_TITLES[c]}" if c else OTHER_CATEGORY
        add(f"## {title}")
        add("")
        if not c:
            add("Commands registered by handler modules that are not in the requirements catalogue "
                "(SASSI-EDU extensions; full names only).")
            add("")
        rows = [[f"**{x.name}**", f"`{x.syntax}`" if x.syntax else "", x.abbrev or "-", x.tier, x.status,
                 md_escape(x.meaning)] for x in lst]
        add(md_table(["Command", "Syntax", "Abbreviation / alias", "Tier", "Implemented", "Meaning"], rows,
                     escape=False))
        add("")
    return "\n".join(L).rstrip() + "\n"


def write_command_reference(out: Path = DEFAULT_COMMAND_REFERENCE) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(command_reference_markdown(), encoding="utf-8")
    return out


# --------------------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="python -m sassi.verify.report",
                                 description="Run the SASSI-EDU verification problems and write the "
                                             "verification manual (and the command reference).")
    ap.add_argument("--out", default=None,
                    help="output Markdown file (default: docs/verification/VERIFICATION_MANUAL.md; with --only or "
                         "--tiers ./VERIFICATION_REPORT_partial.md, so a partial run never replaces the manual)")
    ap.add_argument("--only", default="", help="comma-separated VP ids (e.g. VP-01,VP-30); default all")
    ap.add_argument("--tiers", default="", help="comma-separated tiers to run (P0,P1,P2); default all")
    ap.add_argument("--skip-slow", action="store_true", help="do not run VPs registered with slow=True")
    ap.add_argument("--figures", action="store_true", help="write the summary figure (PNG, matplotlib Agg)")
    ap.add_argument("--figure-dir", default="", help="figure directory (default: <out dir>/figures)")
    ap.add_argument("--keep", action="store_true", help="keep the VP work directories (temporary dir)")
    ap.add_argument("--commands", default=str(DEFAULT_COMMAND_REFERENCE),
                    help="command reference output (Markdown); '' to skip")
    ap.add_argument("--commands-only", action="store_true", help="only write the command reference")
    ap.add_argument("--no-commands", action="store_true", help="do not write the command reference")
    ap.add_argument("--list", action="store_true", help="list the registered VPs and exit")
    ap.add_argument("--quiet", action="store_true", help="no per-VP progress lines")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.list:
        reg = load_all()
        for k in sorted(reg, key=vp_sort_key):
            p = reg[k]
            print(f"{k:8s} {p.tier}  {'slow ' if p.slow else '     '}{p.title}")
        return 0
    if (args.commands_only or not args.no_commands) and args.commands:
        path = write_command_reference(Path(args.commands))
        print(f"command reference written to {path}", file=sys.stderr)
    if args.commands_only:
        return 0
    only = [t for t in args.only.split(",") if t.strip()] or None
    tiers = [t for t in args.tiers.split(",") if t.strip()] or None
    shown = [str(a).replace(str(ROOT) + "/", "") for a in (argv if argv is not None else sys.argv[1:])]
    cmd = "python -m sassi.verify.report" + "".join(f" {a}" for a in shown)
    if args.out:
        out = Path(args.out)
    else:
        out = Path.cwd() / PARTIAL_REPORT if (only or tiers) else DEFAULT_MANUAL
    try:
        path, runs = generate_manual(out, only=only, tiers=tiers, skip_slow=args.skip_slow,
                                     figures=args.figures,
                                     figure_dir=Path(args.figure_dir) if args.figure_dir else None,
                                     keep=args.keep, progress=None if args.quiet else _default_progress,
                                     command=cmd)
    except KeyError as exc:
        print(f"error: {exc.args[0] if exc.args else exc}", file=sys.stderr)
        return 2
    bad = [r.id for r in runs if r.status in (STATUS_FAIL, STATUS_ERROR)]
    print(f"verification manual written to {path}: {len(runs)} problems, "
          f"{sum(1 for r in runs if r.status == STATUS_PASS)} passed, {len(bad)} not passing"
          + (f" ({', '.join(bad)})" if bad else ""), file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":     # pragma: no cover
    raise SystemExit(main())
