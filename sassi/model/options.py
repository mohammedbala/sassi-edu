"""Module option records: the storage of the analysis-option commands (requirements section 3.3).

Most analysis options are set by commands that replace a whole *record* (SITE, SOIL, HOUSE,
ANALYS, MOTION, STRESS, EQUAKE, POINT, FORCE, AOPT, INCOH, WPASS, RELD, BINOUT, MOPT, the
extension X-commands ...) or that create/overwrite one *entry* of an indexed family (WAVE by wave
type, SPRO/SACC/SRS/SSTR/SSAF/SFOU by sublayer, DYNP by (label, point), RSIN/RSOUT/ACCIN/ACCOUT/
TPSD by spectrum number, SYMM by plane, ME by motion, EDUOPT by key ...).

Design (extension point for the option-command work package)
-----------------------------------------------------------
:class:`Record` is a *positional* record that keeps every argument **as the token text the user
typed** (``None`` = blank field = documented default).  Keeping text makes WRITE -> INP exact
(UT-03) and lets a record exist before its typed definition does.  A typed record is a subclass
that declares ``COMMAND`` and ``FIELDS`` (name, type, default per argument position) and is
registered with :func:`register_record_type`; it inherits storage, serialisation, equality and
WRITE support, and adds named, typed access::

    @register_record_type
    class SiteRecord(Record):
        COMMAND = "SITE"
        FIELDS = (Field("opmode", int, 0), Field("mode1", int, 1), Field("fstep", float, 0.0), ...)

    rec = model.options.record("SITE")      # a SiteRecord when the type is registered
    rec.fstep                                 # -> float, default when blank
    rec.set("nft", 4096)

Positional access is **1-based** (argument 1 is the first value after the command name), matching
the argument numbering of the manual and of docs/spec.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, ClassVar, Dict, Hashable, List, Optional, Sequence, Tuple, Type

from .values import NumberError, fmt_num, parse_float, parse_int

Token = Optional[str]


@dataclass(frozen=True)
class Field:
    """One positional argument of a typed option record."""
    name: str
    type: type = float          # int, float or str
    default: Any = None
    doc: str = ""


def _norm_token(v: Any) -> Token:
    """Normalise a value to stored token text (None for blank)."""
    if v is None:
        return None
    if isinstance(v, str):
        t = v.strip()
        return t if t else None
    if isinstance(v, (int, float)):
        return fmt_num(v)
    return str(v).strip() or None


class Record:
    """Positional option record (generic when ``FIELDS`` is empty)."""

    COMMAND: ClassVar[str] = ""
    FIELDS: ClassVar[Tuple[Field, ...]] = ()

    def __init__(self, command: Optional[str] = None, values: Sequence[Any] = (), legacy: str = ""):
        self.command: str = (command or self.COMMAND).upper()
        self.values: List[Token] = [_norm_token(v) for v in values]
        self.legacy: str = legacy or ""          # legacy format string (SOIL, rule L9)
        self._trim()

    # ------------------------------------------------------------------ positional access
    def _trim(self) -> None:
        while self.values and self.values[-1] is None:
            self.values.pop()

    def __len__(self) -> int:
        return len(self.values)

    def arg(self, k: int) -> Token:
        """Token text of argument ``k`` (1-based), or None when blank / not given."""
        if k < 1:
            raise IndexError("record arguments are 1-based")
        return self.values[k - 1] if k <= len(self.values) else None

    def given(self, k: int) -> bool:
        return self.arg(k) is not None

    def number(self, k: int, default: Optional[float] = None) -> Optional[float]:
        """Argument ``k`` as a float (default when blank or not numeric)."""
        t = self.arg(k)
        if t is None:
            return default
        try:
            return parse_float(t)
        except NumberError:
            return default

    def integer(self, k: int, default: Optional[int] = None) -> Optional[int]:
        """Argument ``k`` as an int (rounded, D-PAR-07); default when blank or not numeric."""
        t = self.arg(k)
        if t is None:
            return default
        try:
            return parse_int(t)[0]
        except NumberError:
            return default

    def set_arg(self, k: int, value: Any) -> None:
        """Set argument ``k`` (1-based); numbers are stored as exact text (D-PAR-17)."""
        if k < 1:
            raise IndexError("record arguments are 1-based")
        while len(self.values) < k:
            self.values.append(None)
        self.values[k - 1] = _norm_token(value)
        self._trim()

    # ------------------------------------------------------------------ named access (typed records)
    @classmethod
    def field_index(cls, name: str) -> int:
        """1-based position of field ``name``."""
        for i, f in enumerate(cls.FIELDS):
            if f.name == name:
                return i + 1
        raise KeyError(f"{cls.COMMAND or cls.__name__}: no field '{name}'")

    def get(self, name: str, default: Any = None) -> Any:
        """Typed value of field ``name`` (its default when the argument is blank)."""
        k = self.field_index(name)
        f = self.FIELDS[k - 1]
        dflt = f.default if default is None else default
        t = self.arg(k)
        if t is None:
            return dflt
        if f.type is int:
            v = self.integer(k)
            return dflt if v is None else v
        if f.type is float:
            v = self.number(k)
            return dflt if v is None else v
        return t

    def set(self, name: str, value: Any) -> None:
        self.set_arg(self.field_index(name), value)

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        for f in type(self).FIELDS:
            if f.name == name:
                return self.get(name)
        raise AttributeError(f"{type(self).__name__} has no attribute '{name}'")

    def is_default(self) -> bool:
        """True when every given argument equals its field default (typed records only).

        Used by WRITE for the extension X-commands, which are written only when non-default
        (requirements section 3.4.R).
        """
        if not self.FIELDS:
            return not self.values
        for i, f in enumerate(self.FIELDS):
            if self.arg(i + 1) is not None and self.get(f.name) != f.default:
                return False
        return len(self.values) <= len(self.FIELDS)

    # ------------------------------------------------------------------ serialisation
    def to_tokens(self) -> List[str]:
        """Argument tokens for WRITE (blank fields as empty strings)."""
        return ["" if v is None else v for v in self.values]

    def to_json(self) -> dict:
        d = {"command": self.command, "values": list(self.values)}
        if self.legacy:
            d["legacy"] = self.legacy
        return d

    @staticmethod
    def from_json(d: dict) -> "Record":
        return make_record(d["command"], d.get("values", []), d.get("legacy", ""))

    def copy(self) -> "Record":
        return make_record(self.command, list(self.values), self.legacy)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Record):
            return NotImplemented
        return (self.command, self.values, self.legacy) == (other.command, other.values, other.legacy)

    def __hash__(self):  # records are mutable; hash by identity
        return id(self)

    def __repr__(self) -> str:
        kind = type(self).__name__
        return f"{kind}({self.command}," + ",".join(self.to_tokens()) + (f", legacy={self.legacy!r}" if self.legacy else "") + ")"


# --------------------------------------------------------------------------------------
# Record type registry
# --------------------------------------------------------------------------------------
RECORD_TYPES: Dict[str, Type[Record]] = {}


def register_record_type(cls: Type[Record]) -> Type[Record]:
    """Class decorator: make ``cls`` the record type of ``cls.COMMAND``."""
    if not cls.COMMAND:
        raise ValueError("typed record needs COMMAND")
    RECORD_TYPES[cls.COMMAND.upper()] = cls
    return cls


def make_record(command: str, values: Sequence[Any] = (), legacy: str = "") -> Record:
    """Instantiate the registered record type of ``command`` (generic :class:`Record` otherwise)."""
    cls = RECORD_TYPES.get(command.upper(), Record)
    return cls(command.upper(), values, legacy)


@register_record_type
class MoptRecord(Record):
    """``MOPT,<incomp>,<matrix>,<mass>,<force>`` (spec 07 section 9.2.21, D-MDL-08).

    incomp 0 include / 1 suppress incompatible modes; matrix 0 mass / 1 weight units for MXM;
    mass and force 0 add / 1 set (overwrite) repeated definitions.  Defaults are the Options >
    Model dialog values (Suppress, Mass Matrix, Overwrite Mass, Overwrite Force).
    """
    COMMAND = "MOPT"
    FIELDS = (Field("incomp", int, 1, "0 include, 1 suppress incompatible modes"),
              Field("matrix", int, 0, "GENERAL mass units 0 mass, 1 weight"),
              Field("mass", int, 1, "repeated nodal masses 0 add, 1 overwrite"),
              Field("force", int, 1, "repeated nodal forces 0 add, 1 overwrite"))


# --------------------------------------------------------------------------------------
# Option store
# --------------------------------------------------------------------------------------
def _key_sort(k: Hashable):
    if isinstance(k, tuple):
        return (2, tuple(_key_sort(x) for x in k))
    if isinstance(k, (int, float)):
        return (0, k, "")
    return (1, 0, str(k))


def _key_to_json(k: Hashable):
    return {"t": list(k)} if isinstance(k, tuple) else k


def _key_from_json(k):
    if isinstance(k, dict) and "t" in k:
        return tuple(k["t"])
    return k


class OptionStore:
    """Record setters, indexed setters and string setters of one model (requirements section 3.3)."""

    def __init__(self) -> None:
        self.records: Dict[str, Record] = {}
        self.indexed: Dict[str, Dict[Hashable, Record]] = {}
        self.strings: Dict[str, str] = {}

    # -- record setters
    def record(self, name: str) -> Optional[Record]:
        return self.records.get(name.upper())

    def ensure_record(self, name: str) -> Record:
        """The record of ``name``, created empty (all defaults) when absent."""
        name = name.upper()
        if name not in self.records:
            self.records[name] = make_record(name)
        return self.records[name]

    def set_record(self, rec: Record) -> None:
        self.records[rec.command.upper()] = rec

    def delete_record(self, name: str) -> None:
        self.records.pop(name.upper(), None)

    # -- indexed setters
    def entries(self, name: str) -> List[Tuple[Hashable, Record]]:
        d = self.indexed.get(name.upper(), {})
        return [(k, d[k]) for k in sorted(d, key=_key_sort)]

    def entry(self, name: str, key: Hashable) -> Optional[Record]:
        return self.indexed.get(name.upper(), {}).get(key)

    def set_entry(self, name: str, key: Hashable, rec: Record) -> None:
        self.indexed.setdefault(name.upper(), {})[key] = rec

    def delete_entry(self, name: str, key: Hashable) -> bool:
        d = self.indexed.get(name.upper())
        if d is None or key not in d:
            return False
        del d[key]
        if not d:
            del self.indexed[name.upper()]
        return True

    # -- string setters
    def string(self, name: str, default: str = "") -> str:
        return self.strings.get(name.upper(), default)

    def set_string(self, name: str, value: str) -> None:
        if value:
            self.strings[name.upper()] = value
        else:
            self.strings.pop(name.upper(), None)

    # -- whole-store helpers
    def names(self) -> List[str]:
        return sorted(set(self.records) | set(self.indexed) | set(self.strings))

    def to_json(self) -> dict:
        return {
            "records": {k: self.records[k].to_json() for k in sorted(self.records)},
            "indexed": {n: [{"key": _key_to_json(k), "record": r.to_json()} for k, r in self.entries(n)]
                        for n in sorted(self.indexed)},
            "strings": {k: self.strings[k] for k in sorted(self.strings)},
        }

    @classmethod
    def from_json(cls, d: dict) -> "OptionStore":
        s = cls()
        for k, r in d.get("records", {}).items():
            s.records[k] = Record.from_json(r)
        for n, items in d.get("indexed", {}).items():
            for it in items:
                s.set_entry(n, _key_from_json(it["key"]), Record.from_json(it["record"]))
        s.strings.update(d.get("strings", {}))
        return s

    def copy(self) -> "OptionStore":
        return OptionStore.from_json(json.loads(json.dumps(self.to_json())))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, OptionStore):
            return NotImplemented
        return self.to_json() == other.to_json()
