"""The command explainer of the Learn features (``GET /api/explain?line=<command line>[&cursor=k]``).

A learner clicks a command line (in a lesson, in the Command History, or while typing in Command
Entry) and sees what it does: the full command name (abbreviations resolved with the interpreter's own
:func:`sassi.prep.registry.lookup`, rule L10), the syntax, the one-line meaning, the tier and the
implementation status (:func:`sassi.verify.report.command_infos`, the source of the Command
Reference), what the command configures, and one row per argument: position, name, the value given in
the line and its meaning.

Argument names and meanings come from, in this order:

1. the typed option records of the module option commands (:data:`sassi.prep.options.OPTION_SPECS`:
   ``record.FIELDS`` name, default and dialog label) and the command records of
   :mod:`sassi.ui.cmdrecords` (EQL, NLSOIL, LOADGEN ...: ``RECORD_TYPES``) -- the names are the ones the
   dialogs and WRITE use;
2. the curated table ``command_args.json`` (next to this file), written from the argument definitions of
   ``docs/spec/07`` to ``11`` and checked against the command handlers: the plain-language meaning of
   every argument, defaults, the codes of coded arguments (``values``) and a summary;
3. otherwise the argument names of the syntax string.

The explanation is descriptive only: nothing is executed.
"""
from __future__ import annotations

import difflib
import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ARGS_FILE = Path(__file__).resolve().parent / "command_args.json"

TIER_TEXT = {"P0": "P0 - core linear SSI chain", "P1": "P1 - important engineering feature",
             "P2": "P2 - advanced or cosmetic feature", "OOS": "outside the scope of this build"}
STATUS_TEXT = {"yes": "implemented in this build",
               "stored only": "parsed and stored for WRITE (no computation in this build)",
               "no": "not available in this build (the command prints a message)"}
#: requirements section 3.4 categories (the headings of the Command Reference)
CATEGORY_TITLES = {
    "3.4.A": "Session, files, models and global options", "3.4.B": "Frequency sets",
    "3.4.C": "Nodes and coordinate systems", "3.4.D": "Groups, elements and properties",
    "3.4.E": "Loads and masses", "3.4.F": "Module analysis options", "3.4.G": "Module run commands",
    "3.4.H": "Model checking", "3.4.I": "Model conditioning and generation",
    "3.4.J": "Cuts, submodels and section calculations", "3.4.K": "File conversion",
    "3.4.L": "Plotting and line mathematics", "3.4.M": "Programming (variables, loops, macros)",
    "3.4.N": "Water modelling", "3.4.O": "Option NON and nonlinear soil", "3.4.P": "Binary databases",
    "3.4.Q": "Thick shell", "3.4.R": "Extension commands of SASSI-EDU",
}
REFERENCE_DOC = "docs/reference/COMMAND_REFERENCE.md"
#: command-record families of sassi.ui.cmdrecords -> the dialog that edits them
RECORD_DIALOGS = {"EQL": ("ANALYSIS", "NONLINEAR"), "P": ("ANALYSIS", "NONLINEAR"), "S": ("ANALYSIS", "NONLINEAR"),
                  "BBCI": ("ANALYSIS", "NONLINEAR"), "BBCX": ("ANALYSIS", "NONLINEAR"),
                  "BBCY": ("ANALYSIS", "NONLINEAR"), "NLSOIL": ("ANALYSIS", "SOIL"),
                  "NLSLAYER": ("ANALYSIS", "SOIL"), "LOADGEN": ("LOADGEN", ""), "LOADGENDYN": ("LOADGENDYN", ""),
                  "LGFILE": ("LOADGEN", ""), "LGNODE": ("LOADGEN", ""), "LGTIME": ("LOADGEN", ""),
                  "LGMAP": ("LOADGEN", ""), "LGOPT": ("LOADGEN", ""), "MOPT": ("MODEL", "")}
MAX_LINE = 4000


# ======================================================================================
# sources (cached: the catalogue does not change while the server runs)
# ======================================================================================
@lru_cache(maxsize=1)
def curated() -> Dict[str, Dict[str, Any]]:
    """The curated argument table (``command_args.json``); empty when the file is missing."""
    try:
        data = json.loads(ARGS_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return {str(k).upper(): v for k, v in data.items() if isinstance(v, dict) and not str(k).startswith("_")}


@lru_cache(maxsize=1)
def infos() -> Dict[str, Any]:
    """``command_infos()`` of the Command Reference, by command name."""
    from ..verify.report import command_infos
    return {i.name: i for i in command_infos()}


@lru_cache(maxsize=1)
def _names() -> Tuple[str, ...]:
    from ..prep.registry import all_commands
    out: List[str] = []
    for spec in all_commands():
        out.append(spec.name)
        out.extend(spec.abbrev)
    return tuple(sorted(set(out)))


def record_fields(name: str) -> Optional[List[Any]]:
    """The fields of the typed record of command ``name`` (option record or command record), or None."""
    from ..prep.options import OPTION_SPECS
    spec = OPTION_SPECS.get(name)
    if spec is not None and spec.record is not None:
        return list(spec.record.FIELDS)
    from ..model.options import RECORD_TYPES
    from . import cmdrecords
    cmdrecords._cmds()                       # the command modules register their record types
    cls = RECORD_TYPES.get(name)
    fields = list(getattr(cls, "FIELDS", ()) or ())
    return fields or None


def dialog_of(name: str) -> Optional[Dict[str, str]]:
    """The GUI dialog that edits command ``name`` (Options > Analysis tab, Model, LOADGEN), or None."""
    from ..prep.options import OPTION_SPECS
    spec = OPTION_SPECS.get(name)
    if spec is not None and spec.tab:
        if spec.tab == "AFWRITE":
            return {"name": "ANALYSIS", "tab": "AFWRITE", "label": "Options > Analysis > AFWRITE tab"}
        return {"name": "ANALYSIS", "tab": spec.tab, "label": f"Options > Analysis > {spec.tab} tab"}
    if name in RECORD_DIALOGS:
        d, tab = RECORD_DIALOGS[name]
        label = {"MODEL": "Options > Model", "LOADGEN": "Modules > ANSYS Eq. Static Load",
                 "LOADGENDYN": "Modules > ANSYS Dynamic Load"}.get(d, f"Options > Analysis > {tab} tab")
        return {"name": d, "tab": tab, "label": label}
    return None


# ======================================================================================
# line parsing
# ======================================================================================
def _split_raw(rest: str) -> List[str]:
    """Raw comma-separated fields of ``rest`` (double quotes keep commas; the quotes are kept)."""
    out, cur, quoted = [], [], False
    for ch in rest:
        if ch == '"':
            quoted = not quoted
        if ch == "," and not quoted:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return out


def _unquote(tok: str) -> str:
    t = tok.strip()
    if len(t) >= 2 and t[0] == t[-1] == '"':
        t = t[1:-1]
    return t


def split_line(line: str) -> Tuple[str, List[str], bool]:
    """``(head, argument tokens, has_comma)`` of a command line (tokens as written, unquoted)."""
    s = line.strip()
    m = re.match(r"^([^,\s]+)\s*(,?)(.*)$", s, re.S)
    if not m:
        return "", [], False
    head, comma, rest = m.group(1), m.group(2), m.group(3)
    if not comma:
        # 'NAME rest' (blank after the name) or just 'NAME'
        rest = rest.strip()
        return head, ([_unquote(t) for t in _split_raw(rest)] if rest else []), False
    return head, [_unquote(t) for t in _split_raw(rest)], True


def cursor_argument(line: str, cursor: int, text_from: Optional[int] = None) -> int:
    """Index of the argument the cursor is in: 0 = the command name, k = argument k."""
    prefix = line[:max(0, min(cursor, len(line)))]
    s = prefix.lstrip()
    m = re.match(r"^[^,\s]*", s)
    after = s[m.end():] if m else s
    if not after:
        return 0
    k, quoted = 0, False
    for ch in after:
        if ch == '"':
            quoted = not quoted
        elif ch == "," and not quoted:
            k += 1
    if k == 0:
        k = 1 if after.strip() else 0
    if text_from:
        k = min(k, text_from)
    return k


def _norm_code(tok: str) -> List[str]:
    """Spellings of a token to look up in a ``values`` map ('1.0' -> '1', words upper case)."""
    t = tok.strip()
    out = [t, t.upper()]
    try:
        f = float(t)
        if f == int(f):
            out.append(str(int(f)))
    except ValueError:
        pass
    return out


def _fmt_default(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if v == int(v) and abs(v) < 1e15:
            return str(int(v))
        return repr(v)
    return str(v)


def _syntax_args(syntax: str) -> List[Dict[str, Any]]:
    """Argument names of a syntax string ``NAME,<a>,[b],...``; a ``...`` repeats the previous one."""
    parts = [p.strip() for p in syntax.split(",")[1:]]
    out: List[Dict[str, Any]] = []
    for p in parts:
        if p in ("...", "…") or p.startswith("..."):
            if out:
                out[-1]["repeat"] = True
            break
        name = p.strip("[]<> ")
        if not name:
            continue
        out.append({"pos": len(out) + 1, "name": name, "meaning": "", "default": ""})
    return out


# ======================================================================================
# the explanation
# ======================================================================================
def _arg_table(name: str, cur: Dict[str, Any], syntax: str) -> Tuple[List[Dict[str, Any]], str]:
    """The argument definitions of command ``name`` and where they come from."""
    cargs = [dict(a) for a in (cur.get("args") or []) if isinstance(a, dict)]
    by_pos = {}
    for a in cargs:
        try:
            by_pos[int(a.get("pos"))] = a
        except (TypeError, ValueError):
            continue
    fields = record_fields(name)
    if fields:
        out = []
        for i, f in enumerate(fields, start=1):
            c = by_pos.get(i, {})
            label = getattr(f, "doc", "") or ""
            out.append({"pos": i, "name": f.name, "label": label,
                        "meaning": c.get("meaning") or label,
                        "default": c.get("default") if c.get("default") not in (None, "") else _fmt_default(f.default),
                        "values": c.get("values") or {}, "repeat": bool(c.get("repeat"))})
        # a curated variable tail beyond the record (rare) is kept
        for p in sorted(k for k in by_pos if k > len(fields)):
            c = by_pos[p]
            out.append({"pos": p, "name": c.get("name", f"arg{p}"), "label": "", "meaning": c.get("meaning", ""),
                        "default": c.get("default", ""), "values": c.get("values") or {}, "repeat": bool(c.get("repeat"))})
        return out, "record"
    if cargs:
        out = [{"pos": int(a.get("pos", i + 1)), "name": str(a.get("name", "")), "label": "",
                "meaning": str(a.get("meaning", "")), "default": str(a.get("default", "") or ""),
                "values": a.get("values") or {}, "repeat": bool(a.get("repeat"))} for i, a in enumerate(cargs)]
        out.sort(key=lambda a: a["pos"])
        return out, "curated"
    out = _syntax_args(syntax)
    for a in out:
        a.update(label="", values={}, repeat=bool(a.get("repeat")))
    return out, "syntax"


def explain(line: str, cursor: Optional[int] = None) -> Dict[str, Any]:
    """Explanation of one command line (see the module docstring); never raises for user input."""
    from ..prep.registry import lookup
    line = (line or "")[:MAX_LINE].replace("\r", " ").replace("\n", " ")
    text = line.strip()
    base: Dict[str, Any] = {"line": line, "known": False}
    if not text:
        base.update(kind="blank", message="Type a command, e.g. N,1,0,0,0")
        return base
    if text.startswith("*") or text.startswith("!"):
        base.update(kind="comment", known=True, name="*", summary=(
            "A comment line: it is echoed in the Command History (green) and in a .pre file documents the "
            "input; it is not executed."))
        return base
    head, tokens, _ = split_line(text)
    spec = lookup(head)
    if spec is None:
        cands = difflib.get_close_matches(head.upper(), _names(), n=5, cutoff=0.6)
        base.update(kind="unknown", typed=head, suggestions=cands,
                    message=(f"{head} is not a SASSI-EDU command. Names are matched exactly or by their "
                             f"documented abbreviation (no prefix matching)."))
        if cursor is not None:
            base["cursor_arg"] = 0
        return base
    name = spec.name
    info = infos().get(name)
    cur = curated().get(name, {})
    syntax = (info.syntax if info else "") or name
    meaning = (info.meaning if info else "") or spec.summary or ""
    tier = getattr(spec, "tier", "") or (info.tier if info else "")
    status = info.status if info else ("no" if spec.placeholder else "yes")
    category = info.category if info else ""
    args, origin = _arg_table(name, cur, syntax)
    # the text argument of TIT / GTIT ... takes the rest of the line (commas included, L7)
    if spec.text_from and len(tokens) >= spec.text_from:
        rest = text.split(",", 1)[1] if "," in text else ""
        fields = _split_raw(rest)
        k = spec.text_from - 1
        tokens = [_unquote(t) for t in fields[:k]] + [",".join(fields[k:]).strip()]
    rows: List[Dict[str, Any]] = []
    rep = next((a for a in args if a.get("repeat")), None)
    n = max(len(tokens), len(args))
    for i in range(1, n + 1):
        a = next((x for x in args if x["pos"] == i), None)
        repeated = False
        if a is None and rep is not None and i > rep["pos"]:
            a, repeated = rep, True
        tok = tokens[i - 1] if i <= len(tokens) else None
        row: Dict[str, Any] = {"pos": i, "name": a["name"] if a else "", "label": a.get("label", "") if a else "",
                               "meaning": a.get("meaning", "") if a else "", "default": a.get("default", "") if a else "",
                               "value": tok, "given": tok is not None and tok.strip() != "", "repeat": repeated}
        if a is None:
            row["extra"] = True
            row["meaning"] = "not an argument of this command (ignored or reported by the command)"
        vals = (a or {}).get("values") or {}
        if vals:
            row["values"] = {str(k): str(v) for k, v in vals.items()}
            if row["given"]:
                for cand in _norm_code(tok):
                    if cand in row["values"]:
                        row["value_meaning"] = row["values"][cand]
                        break
        if i > len(tokens) and repeated:
            continue
        rows.append(row)
    out = dict(base)
    out.update({
        "kind": "command", "known": True, "typed": head, "name": name,
        "abbreviated": head.upper() != name, "abbrev": list(spec.abbrev),
        "syntax": syntax, "meaning": meaning, "summary": cur.get("summary") or meaning,
        "notes": cur.get("notes", ""), "source": cur.get("source", ""),
        "tier": tier, "tier_text": TIER_TEXT.get(tier, tier), "status": status,
        "status_text": STATUS_TEXT.get(status, status), "class": spec.cls, "category": category,
        "category_title": CATEGORY_TITLES.get(category, ""),
        "configures": cur.get("configures") or _configures(name),
        "dialog": dialog_of(name), "args": rows, "args_from": origin,
        "help": {"doc": REFERENCE_DOC, "anchor": _category_slug(category)},
    })
    if cursor is not None:
        out["cursor_arg"] = cursor_argument(line, int(cursor), spec.text_from)
    return out


def _configures(name: str) -> str:
    """A fallback for commands without a curated entry: the module tab or the module run."""
    d = dialog_of(name)
    if d and d.get("tab") and d["name"] == "ANALYSIS":
        return f"{d['tab']} module options" if d["tab"] not in ("AFWRITE", "NONLINEAR") else d["label"]
    if name.startswith("RUN") and len(name) > 3:
        return f"Module run: {name[3:]}"
    return ""


def _category_slug(category: str) -> str:
    """The Command Reference anchor of a category heading ('3.4.A Session, ...' -> '34a-session-...')."""
    title = CATEGORY_TITLES.get(category)
    if not title:
        return ""
    from .markdown import slugify
    return slugify(f"{category} {title}")


def coverage(names: List[str]) -> List[str]:
    """The names of ``names`` without a curated entry or record fields (for the tests)."""
    out = []
    for n in names:
        if n in curated():
            continue
        if record_fields(n):
            continue
        out.append(n)
    return out
