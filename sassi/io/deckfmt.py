"""Generic keyword-structured text format of SASSI-EDU module input decks (D-GEN-04).

A deck is written by AFWRITE (``<model>.<ext>``) and read by exactly one module.  Layout::

    SASSI-EDU SITE DECK v1
    * comment lines start with '*' (or '#')
    [params]
    title = "Demo surface mat"
    opmode = 0
    delt = 0.005
    [table layers] no thick weight vp vs dp ds
    1 10.0 0.12 1500.0 800.0 0.05 0.05
    [end]

* Scalars are ``key = value``; strings are JSON-quoted when they contain blanks, quotes,
  ``=`` or are empty.  Lists are written as ``key = [v1, v2, ...]`` (JSON).
* ``[table <name>] <col1> <col2> ...`` starts a table; each following line is a row of
  whitespace-separated tokens (strings JSON-quoted when needed); ``[end]`` closes it.
* Types are imposed by the schema in :mod:`sassi.io.decks`.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np

PathLike = Union[str, os.PathLike]
HEADER_RE = re.compile(r"^SASSI-EDU\s+(\w+)\s+DECK\s+v(\d+)\s*$")
_NEEDS_QUOTE = re.compile(r'[\s"=\[\]\'*#,]')


@dataclass
class Table:
    columns: List[str]
    rows: List[List[Any]] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.rows)

    def col(self, name: str) -> np.ndarray:
        j = self.columns.index(name)
        return np.asarray([r[j] for r in self.rows])

    def as_dicts(self) -> List[Dict[str, Any]]:
        return [dict(zip(self.columns, r)) for r in self.rows]

    def append(self, row: Union[Sequence[Any], Dict[str, Any]]) -> None:
        if isinstance(row, dict):
            row = [row[c] for c in self.columns]
        if len(row) != len(self.columns):
            raise ValueError(f"row has {len(row)} values, table has {len(self.columns)} columns")
        self.rows.append(list(row))


@dataclass
class RawDeck:
    module: str
    version: int
    params: Dict[str, str] = field(default_factory=dict)       # raw text values
    tables: Dict[str, Table] = field(default_factory=dict)      # raw string tokens


def fmt_token(v: Any) -> str:
    """Format one scalar token."""
    if isinstance(v, (bool, np.bool_)):
        return "1" if v else "0"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    if isinstance(v, (list, tuple, np.ndarray)):
        return json.dumps([_py(x) for x in v])
    s = "" if v is None else str(v)
    if s == "" or _NEEDS_QUOTE.search(s):
        return json.dumps(s)
    return s


def _py(x):
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    return x


def split_tokens(line: str) -> List[str]:
    """Split a table row on whitespace, honouring JSON double-quoted strings."""
    out: List[str] = []
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c.isspace():
            i += 1
            continue
        if c == '"':
            j = i + 1
            while j < n:
                if line[j] == "\\":
                    j += 2
                    continue
                if line[j] == '"':
                    break
                j += 1
            out.append(json.loads(line[i:j + 1]))
            i = j + 1
        else:
            j = i
            while j < n and not line[j].isspace():
                j += 1
            out.append(line[i:j])
            i = j
    return out


def parse_value(text: str) -> Any:
    """Parse a raw ``key = value`` right-hand side into str/list (type coercion is the schema's job)."""
    t = text.strip()
    if t.startswith('"') or t.startswith("["):
        return json.loads(t)
    return t


def write_raw(path: PathLike, module: str, params: Dict[str, Any], tables: Dict[str, Table],
              version: int = 1, comments: Sequence[str] = ()) -> Path:
    path = Path(path)
    lines = [f"SASSI-EDU {module.upper()} DECK v{version}"]
    for c in comments:
        lines.append(f"* {c}")
    lines.append("[params]")
    for k, v in params.items():
        lines.append(f"{k} = {fmt_token(v)}")
    for name, tab in tables.items():
        lines.append(f"[table {name}] " + " ".join(tab.columns))
        for r in tab.rows:
            lines.append(" ".join(fmt_token(v) for v in r))
        lines.append("[end]")
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(tmp, path)
    return path


def read_raw(path: PathLike) -> RawDeck:
    path = Path(path)
    text = path.read_text(encoding="utf-8").splitlines()
    if not text:
        raise ValueError(f"{path.name}: empty deck")
    m = HEADER_RE.match(text[0].strip())
    if not m:
        raise ValueError(f"{path.name}: missing 'SASSI-EDU <MODULE> DECK v<n>' header")
    deck = RawDeck(module=m.group(1).upper(), version=int(m.group(2)))
    section: Optional[str] = None
    table: Optional[Table] = None
    for lineno, raw in enumerate(text[1:], start=2):
        line = raw.strip()
        if not line or line[0] in "*#":
            continue
        if line.lower() == "[params]":
            section, table = "params", None
            continue
        if line.lower().startswith("[table "):
            close = line.index("]")
            name = line[7:close].strip()
            cols = line[close + 1:].split()
            table = Table(columns=cols)
            deck.tables[name] = table
            section = "table"
            continue
        if line.lower() == "[end]":
            section, table = None, None
            continue
        if section == "params":
            if "=" not in line:
                raise ValueError(f"{path.name}:{lineno}: expected 'key = value'")
            k, v = line.split("=", 1)
            deck.params[k.strip()] = v.strip()
        elif section == "table" and table is not None:
            toks = split_tokens(line)
            if len(toks) != len(table.columns):
                raise ValueError(f"{path.name}:{lineno}: {len(toks)} values for {len(table.columns)} columns")
            table.rows.append(toks)
        else:
            raise ValueError(f"{path.name}:{lineno}: data outside a section")
    return deck
