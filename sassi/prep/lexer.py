"""Lexical rules of the ``.pre`` command language (requirements section 3.1, rules L1-L9).

A command line is ``Keyword, p1, p2, ..., pn`` (L1).  This module splits a line *after* macro
and variable substitution (L13) into the command name and its argument tokens:

* L1 / D-PAR-02 -- a single run of blanks or tabs may replace the **first** comma only
  (``EDGE 1,0,0,1``);
* L2 -- whitespace around every token is trimmed;
* L4 / D-PAR-03 -- a line whose first non-blank character is ``*`` is a comment; there are no
  inline comments;
* L5 -- an empty field (nothing or blanks between commas) is kept as an empty token ``''``, which
  handlers read as "use the documented default";
* L7 / D-PAR-06 -- a text argument in the **last** position takes the rest of the line verbatim,
  commas and unbalanced quotes included (``text_from``: only the fields before it are scanned);
  a double-quoted token ``"..."`` may be used anywhere to embed commas;
* L8 / D-PAR-08 -- node/element lists are integers or inclusive ranges ``a-b`` separated by
  commas, blanks, tabs or ``;`` (:func:`parse_id_list`);
* L9 / D-SOL-10 -- a final parenthesised token ``(...)`` of a legacy SOIL line is one token
  (``paren_last``).

Number parsing (L6) is in :mod:`sassi.model.values` and re-exported here.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

from ..model.values import NumberError, fmt_num, is_number, parse_float, parse_int  # noqa: F401 (re-export)

MAX_LINE = 3000        # L14 / D-PAR-11: maximum line length after substitution


class LexError(ValueError):
    """Malformed command text (unterminated quote, bad list ...)."""


def is_comment(line: str) -> bool:
    """L4: the first non-blank character is ``*``."""
    s = line.lstrip()
    return s.startswith("*")


def is_blank(line: str) -> bool:
    return not line.strip()


def split_head(line: str) -> Tuple[str, str]:
    """Split a command line into ``(name, rest)``.

    ``name`` is the keyword as typed (case kept); ``rest`` is the raw text after the first
    separator: a comma, or a run of blanks/tabs (optionally followed by the first comma), rule L1.
    """
    s = line.strip()
    i = 0
    n = len(s)
    while i < n and s[i] not in ", \t":
        i += 1
    name = s[:i]
    if i >= n:
        return name, ""
    if s[i] == ",":
        return name, s[i + 1:]
    # a run of blanks replaces the first comma; a comma right after the blanks is that comma
    j = i
    while j < n and s[j] in " \t":
        j += 1
    if j < n and s[j] == ",":
        return name, s[j + 1:]
    return name, s[j:]


@dataclass
class Field:
    start: int          # offset of the field in ``rest``
    end: int            # offset of the terminating comma (or len(rest))
    text: str           # token text (trimmed, quotes removed)
    quoted: bool = False


def _scan(rest: str, limit: Optional[int] = None) -> Tuple[List[Field], Optional[int]]:
    """Quote-aware comma split of ``rest``; stops after ``limit`` fields when ``limit`` is given.

    Returns ``(fields, tail)``: ``tail`` is the offset of the text after the ``limit``-th
    separating comma (``None`` when the line has no more fields).  Text after that offset is
    *not* scanned, so a rest-of-line text argument may hold anything, an unbalanced ``"``
    included (L7 / D-PAR-06).
    """
    if not rest.strip():
        return [], None
    if limit is not None and limit <= 0:
        return [], 0
    out: List[Field] = []
    n = len(rest)
    i = 0
    while True:
        j = i
        while j < n and rest[j] in " \t":
            j += 1
        if j < n and rest[j] == '"':
            k = rest.find('"', j + 1)
            if k < 0:
                raise LexError("unterminated quoted string")
            c = rest.find(",", k + 1)
            end = n if c < 0 else c
            text = rest[j + 1:k] + rest[k + 1:end].strip()
            out.append(Field(i, end, text, True))
        else:
            c = rest.find(",", i)
            end = n if c < 0 else c
            out.append(Field(i, end, rest[i:end].strip(), False))
        if c < 0:
            return out, None
        i = c + 1
        if limit is not None and len(out) >= limit:
            return out, i


def scan_fields(rest: str) -> List[Field]:
    """Quote-aware comma split of ``rest`` (no tokens for an empty/blank ``rest``).

    A field whose first non-blank character is ``"`` extends to the matching ``"``; commas inside
    the quotes do not split (D-PAR-06).
    """
    return _scan(rest)[0]


def _unquote_text(t: str) -> str:
    t = t.strip()
    if len(t) >= 2 and t[0] == '"' and t[-1] == '"' and t.count('"') == 2:
        return t[1:-1]
    return t


@dataclass
class Lexed:
    """Result of :func:`split_args`."""
    tokens: List[str] = field(default_factory=list)
    quoted: List[bool] = field(default_factory=list)
    legacy_paren: str = ""      # the parenthesised legacy format token, if one was found


def split_args(rest: str, text_from: Optional[int] = None, paren_last: bool = False) -> Lexed:
    """Split the argument text of a command into tokens.

    ``text_from`` (1-based): from that argument on, the rest of the line is one text token
    (L7).  ``paren_last``: a final ``(...)`` token is kept whole (L9, legacy SOIL).
    """
    res = Lexed()
    if text_from is not None:
        # only the fields *before* the text argument are scanned (an unbalanced quote inside
        # the text is data, L7); the text itself is the rest of the line verbatim
        head, tail = _scan(rest, text_from - 1)
        if tail is not None:
            res.tokens = [f.text for f in head] + [_unquote_text(rest[tail:])]
            res.quoted = [f.quoted for f in head] + [False]
            return res
        fields = head                     # fewer fields than text_from: ordinary tokens
    else:
        fields = scan_fields(rest)
    if paren_last and rest.rstrip().endswith(")"):
        for idx, f in enumerate(fields):
            if f.text.startswith("(") and not f.quoted:
                tail = rest[f.start:].strip()
                if tail.count("(") == tail.count(")"):
                    res.tokens = [g.text for g in fields[:idx]] + [tail]
                    res.quoted = [g.quoted for g in fields[:idx]] + [False]
                    res.legacy_paren = tail
                    return res
    res.tokens = [f.text for f in fields]
    res.quoted = [f.quoted for f in fields]
    return res


# --------------------------------------------------------------------------------------
# Node / element lists (L8, D-PAR-08)
# --------------------------------------------------------------------------------------
_RANGE_RE = re.compile(r"\s*-\s*")


def parse_id_list(items: Sequence[str]) -> List[int]:
    """Expand list tokens such as ``['1', '3-6 10;12']`` to ``[1, 3, 4, 5, 6, 10, 12]``.

    Items are integers or inclusive ranges ``a-b`` with ``a <= b``; separators are commas
    (already split), blanks, tabs and ``;``.  A descending range is an error (D-PAR-08).  The
    order is kept and duplicates are not removed (CHECK reports them).
    """
    text = " ".join(str(s) for s in items).replace(";", " ").replace("\t", " ")
    text = _RANGE_RE.sub("-", text.strip())
    out: List[int] = []
    for tok in text.split():
        k = tok.find("-", 1)
        if k > 0:
            a_s, b_s = tok[:k], tok[k + 1:]
            try:
                a, _ = parse_int(a_s)
                b, _ = parse_int(b_s)
            except NumberError:
                raise LexError(f"'{tok}' is not a range of integers") from None
            if a > b:
                raise LexError(f"descending range '{tok}' (D-PAR-08: a-b needs a <= b)")
            out.extend(range(a, b + 1))
        else:
            try:
                v, exact = parse_int(tok)
            except NumberError:
                raise LexError(f"'{tok}' is not an integer") from None
            if not exact:
                raise LexError(f"'{tok}' is not an integer")
            out.append(v)
    return out


def compress_ids(ids: Sequence[int]) -> List[str]:
    """Inverse of :func:`parse_id_list` for WRITE: ``[1,2,3,7]`` -> ``['1-3', '7']`` (order kept)."""
    out: List[str] = []
    ids = list(ids)
    i = 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        if j - i >= 2:
            out.append(f"{ids[i]}-{ids[j]}")
        else:
            out.extend(str(v) for v in ids[i:j + 1])
        i = j + 1
    return out


# --------------------------------------------------------------------------------------
# Formatting (WRITE)
# --------------------------------------------------------------------------------------
def format_token(value) -> str:
    """Token text for WRITE: numbers exactly (D-PAR-17); strings quoted when they hold commas."""
    if value is None:
        return ""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return fmt_num(value)
    s = str(value)
    if s and ("," in s or s != s.strip() or s.startswith('"')):
        return '"' + s + '"'
    return s


def join_command(name: str, tokens: Sequence, text_last: bool = False) -> str:
    """``NAME,t1,t2,...``; with ``text_last`` the last token is written verbatim (rest-of-line text)."""
    toks = list(tokens)
    if text_last and toks:
        last = "" if toks[-1] is None else str(toks[-1])
        parts = [format_token(t) for t in toks[:-1]] + [last]
    else:
        parts = [format_token(t) for t in toks]
    while parts and parts[-1] == "":
        parts.pop()
    return ",".join([name] + parts)
