"""Dialog bindings of the *command records* outside :data:`sassi.prep.options.OPTION_SPECS`.

Option NON (EQL, P, S, BBC*), SOIL-NON (NLSOIL, NLSLAYER) and Option A (LOADGEN, LOADGENDYN, LGFILE,
LGNODE, LGTIME, LGMAP, LGOPT) keep their data in records registered with
:func:`sassi.model.options.register_record_type` by their command modules -- not in the typed option
records of :mod:`sassi.prep.options` on which the rest of :mod:`sassi.ui.dialogs` is built.  This module
gives each of them a *family* adapter:

* ``read(model)`` -- what the dialog shows, i.e. what the analysis uses (a blank EQL <NonLinOpts> shows the
  element types of the P and S records; an absent LGTIME shows the default ``V,1``);
* ``lines(model, posted, errors)`` -- the command text of a dialog commit (rule L17): nothing when the
  posted values equal the shown ones, otherwise the commands a user would type, deletions first (PDEL,
  DELSPR, DELNLS, DELBBC).

Dialog paths (:mod:`sassi.ui.dialogs` lays out the fields, ``static/dialogs.js`` edits them):

======================  ==========================================================================
``%NAME.field``         field of a record family: ``%EQL.disp``, ``%NLSOIL.bedint``, ``%LOADGEN.data``,
                        ``%LGFILE.HOUSE`` (the file box of key HOUSE), ``%LGNODE.D`` (a node list),
                        ``%LGTIME.crit`` ...  (payload ``xrecords``)
``%NAME``               table family, one row per entry: ``%P``, ``%S``, ``%NLSLAYER`` (``xtables``)
``%NAME[sel].field``    entry of a keyed family chosen by the dialog selector ``sel``:
                        ``%BBC[bbc].points`` (``xindexed``; a ``null`` entry deletes the curve)
======================  ==========================================================================

Validation (UI-06): the commands of these records check their arguments themselves, so a commit is first
executed in a scratch interpreter holding a copy of the model (:func:`dry_run`); an error refuses the
whole commit with the command's own message.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

from ..model.options import RECORD_TYPES, Record, make_record
from ..model.values import NumberError, fmt_num, parse_float, parse_int
from ..prep.lexer import MAX_LINE, compress_ids, join_command


def _cmds():
    """The command modules of the families (importing them registers their record types)."""
    from ..prep.registry import load_commands
    load_commands()
    from ..prep.commands import loadgen_cmds, nonlinear_cmds, soilnon_cmds
    return nonlinear_cmds, soilnon_cmds, loadgen_cmds


def _typed(rec: Optional[Record]) -> Optional[Record]:
    """``rec`` as its registered record type (a record stored before the type was registered is generic)."""
    if rec is None:
        return None
    cls = RECORD_TYPES.get(rec.command.upper())
    if cls is not None and not isinstance(rec, cls):
        return make_record(rec.command, list(rec.values), rec.legacy)
    return rec


def record(model, name: str) -> Optional[Record]:
    return _typed(model.options.record(name))


def entries(model, name: str) -> List[Tuple[Any, Record]]:
    return [(k, _typed(r)) for k, r in model.options.entries(name)]


def _tok(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, (int, float)):
        return fmt_num(v)
    return str(v).strip()


def _int0(t: Any, default: int = 0) -> int:
    try:
        return parse_int(str(t))[0]
    except (NumberError, TypeError, ValueError):
        return default


class _Parser:
    """Reads dialog values (numbers, check boxes 0/1, text) and collects the messages of bad values."""

    def __init__(self, family: str, errors: List[str]):
        self.family = family
        self.errors = errors

    def bad(self, text: str) -> None:
        self.errors.append(f"{self.family}: {text}")

    def num(self, v: Any, what: str, integer: bool = False, default: Any = None) -> Any:
        t = _tok(v)
        if t == "":
            return default
        try:
            if integer:
                val, exact = parse_int(t)
                if not exact:
                    self.bad(f"{what}: {t} is not an integer")
                    return None
                return val
            return parse_float(t)
        except NumberError:
            self.bad(f"{what}: '{t}' is not a number")
            return None

    def flag(self, v: Any, what: str, default: int = 0) -> Optional[int]:
        x = self.num(v, what, integer=True, default=default)
        if x is not None and x not in (0, 1):
            self.bad(f"{what} must be 0 or 1")
            return None
        return x

    def word(self, v: Any, what: str, allowed: Dict[str, Any], default: Any) -> Any:
        t = _tok(v).upper()
        if not t:
            return default
        if t not in allowed:
            self.bad(f"{what} '{t}' is not one of {', '.join(sorted(set(str(k) for k in allowed)))}")
            return None
        return allowed[t]


# ======================================================================================
# Record families: %NAME.field (payload "xrecords")
# ======================================================================================
class RecordFamily:
    """A record shown field by field; on a change the complete command line is written (manual-identical)."""
    name = ""
    kind = "record"

    def read(self, model) -> Dict[str, Any]:            # pragma: no cover - abstract
        raise NotImplementedError

    def parse(self, vals: Dict[str, Any], P: _Parser) -> Dict[str, Any]:     # pragma: no cover - abstract
        raise NotImplementedError

    def emit(self, model, new: Dict[str, Any], old: Dict[str, Any]) -> List[str]:   # pragma: no cover
        raise NotImplementedError

    def lines(self, model, posted: Dict[str, Any], errors: List[str]) -> List[str]:
        if not isinstance(posted, dict):
            errors.append(f"{self.name}: the dialog values must be an object")
            return []
        shown = self.read(model)
        unknown = sorted(k for k in posted if k not in shown)
        if unknown:
            errors.append(f"{self.name}: unknown field(s) {', '.join(unknown)}")
            return []
        n0 = len(errors)
        new = self.parse(dict(shown, **posted), _Parser(self.name, errors))
        old = self.parse(shown, _Parser(self.name, []))
        if len(errors) > n0 or new == old:
            return []
        return self.emit(model, new, old)


class _TypedFamily(RecordFamily):
    """A family whose dialog fields are the fields of its typed record (numbers only)."""

    def cls(self):
        _cmds()                                          # the command module registers the record type
        return RECORD_TYPES[self.name]

    def read(self, model) -> Dict[str, Any]:
        cls = self.cls()                                 # first: record() types the stored record with it
        rec = record(model, self.name)
        return {f.name: (rec.get(f.name) if rec is not None else f.default) for f in cls.FIELDS}

    def parse(self, vals: Dict[str, Any], P: _Parser) -> Dict[str, Any]:
        return {f.name: P.num(vals.get(f.name), f"<{f.name}>", integer=f.type is int, default=f.default)
                for f in self.cls().FIELDS}

    def emit(self, model, new: Dict[str, Any], old: Dict[str, Any]) -> List[str]:
        return [join_command(self.name, [_tok(new[f.name]) for f in self.cls().FIELDS])]


class EqlFamily(_TypedFamily):
    """``EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>`` (NONLINEAR tab, Global Modeling Options).

    A blank <NonLinOpts> means "the element types of the P and S records" (build_eql): the dialog shows
    those bits and keeps the argument blank unless the user changes them."""
    name = "EQL"

    @staticmethod
    def derived(model) -> int:
        return (1 if model.options.entries("P") else 0) | (2 if model.options.entries("S") else 0)

    @staticmethod
    def opts_blank(model) -> bool:
        rec = record(model, "EQL")
        return rec is None or not rec.given(2)

    def read(self, model) -> Dict[str, Any]:
        vals = super().read(model)
        if self.opts_blank(model):
            vals["nonlinopts"] = self.derived(model)
        return vals

    def emit(self, model, new: Dict[str, Any], old: Dict[str, Any]) -> List[str]:
        toks = [_tok(new[f.name]) for f in self.cls().FIELDS]
        if self.opts_blank(model) and new["nonlinopts"] == old["nonlinopts"]:
            toks[1] = ""                                  # still "the types of the P/S records"
        return [join_command("EQL", toks)]


class NlsoilFamily(_TypedFamily):
    """``NLSOIL,<Opt>,<NSTimeSunInc>,<DispConv>,<ForceConv>,<EqualIt>,<BedInt>,<NLDampType>,<MMmult>,<SMmult>``
    (SOIL tab, Nonlinear Soil Behavior)."""
    name = "NLSOIL"


class LgoptFamily(_TypedFamily):
    """``LGOPT,<digits>,<opmode>,<rest>`` (LOADGEN dialogs)."""
    name = "LGOPT"

    def parse(self, vals: Dict[str, Any], P: _Parser) -> Dict[str, Any]:
        d = P.num(vals.get("digits"), "APDL significant digits", integer=True, default=12)
        if d is not None and not 6 <= d <= 17:
            P.bad("APDL significant digits must be 6..17")
        return {"digits": d, "opmode": P.flag(vals.get("opmode"), "<opmode>"),
                "rest": P.flag(vals.get("rest"), "<rest>", default=1)}

    def emit(self, model, new: Dict[str, Any], old: Dict[str, Any]) -> List[str]:
        return [join_command("LGOPT", [_tok(new["digits"]), _tok(new["opmode"]), _tok(new["rest"])])]


class LoadgenFamily(RecordFamily):
    """``LOADGEN,<data>,<multi>,<rotdisp>,<rotacc>,<masstype>,<genmass>,<source>`` (ANSYS Static Load Converter).

    <data> and <masstype> are shown as their codes (the words DISP, ACC ... LUMPED, MASTER are the same
    codes); a commit writes the codes."""
    name = "LOADGEN"
    FLAGS = ("multi", "rotdisp", "rotacc", "genmass")

    def read(self, model) -> Dict[str, Any]:
        _, _, LC = _cmds()
        rec = record(model, "LOADGEN")

        def t(k: int, d: str) -> str:
            v = rec.arg(k) if rec is not None else None
            return d if v is None else v
        return {"data": LC.DATA_WORDS.get(t(1, "2").upper(), 2), "multi": _int0(t(2, "0")),
                "rotdisp": _int0(t(3, "0")), "rotacc": _int0(t(4, "0")),
                "masstype": LC.MASS_WORDS.get(t(5, "1").upper(), 1), "genmass": _int0(t(6, "0")),
                "source": t(7, "RESULTS").upper()}

    def parse(self, vals: Dict[str, Any], P: _Parser) -> Dict[str, Any]:
        _, _, LC = _cmds()
        out = {"data": P.word(vals.get("data"), "<data>", LC.DATA_WORDS, 2),
               "masstype": P.word(vals.get("masstype"), "<masstype>", LC.MASS_WORDS, 1),
               "source": P.word(vals.get("source"), "<source>", {s: s for s in LC.SOURCES}, "RESULTS")}
        for k in self.FLAGS:
            out[k] = P.flag(vals.get(k), f"<{k}>")
        return out

    def emit(self, model, new: Dict[str, Any], old: Dict[str, Any]) -> List[str]:
        return [join_command("LOADGEN", [_tok(new[k]) for k in ("data", "multi", "rotdisp", "rotacc", "masstype",
                                                                "genmass", "source")])]


class LoadgenDynFamily(RecordFamily):
    """``LOADGENDYN,<alpha>,<beta>,<method>,<refnode>,<source>,<rotdisp>,<rotacc>,<zeta>,<f1>,<f2>,<gfopt>,<gmult>``
    (ANSYS Dynamic Load Converter).  With a damping ratio zeta > 0 the command computes alpha and beta, so
    they are written blank (given, they would only be reported as ignored)."""
    name = "LOADGENDYN"
    ORDER = ("alpha", "beta", "method", "refnode", "source", "rotdisp", "rotacc", "zeta", "f1", "f2", "gfopt", "gmult")

    def read(self, model) -> Dict[str, Any]:
        _cmds()
        rec = record(model, "LOADGENDYN")
        cls = RECORD_TYPES["LOADGENDYN"]
        vals = {f.name: (rec.get(f.name) if rec is not None else f.default) for f in cls.FIELDS}
        vals["method"] = str(vals["method"] or "REL").upper()
        vals["source"] = str(vals["source"] or "RESULTS").upper()
        return vals

    def parse(self, vals: Dict[str, Any], P: _Parser) -> Dict[str, Any]:
        _, _, LC = _cmds()
        out: Dict[str, Any] = {}
        for k in ("alpha", "beta", "zeta", "f1", "f2"):
            out[k] = P.num(vals.get(k), f"<{k}>", default=0.0)
        out["gmult"] = P.num(vals.get("gmult"), "<gmult>", default=1.0)
        out["refnode"] = P.num(vals.get("refnode"), "<refnode>", integer=True, default=0)
        out["gfopt"] = P.flag(vals.get("gfopt"), "<gfopt>")
        out["rotdisp"] = P.flag(vals.get("rotdisp"), "<rotdisp>")
        out["rotacc"] = P.flag(vals.get("rotacc"), "<rotacc>")
        out["method"] = P.word(vals.get("method"), "<method>", {m: m for m in LC.METHODS}, "REL")
        out["source"] = P.word(vals.get("source"), "<source>", {s: s for s in LC.SOURCES}, "RESULTS")
        return out

    def emit(self, model, new: Dict[str, Any], old: Dict[str, Any]) -> List[str]:
        toks = [_tok(new[k]) for k in self.ORDER]
        if (new["zeta"] or 0) > 0:
            toks[0] = toks[1] = ""
        return [join_command("LOADGENDYN", toks)]


class LgfileFamily(RecordFamily):
    """``LGFILE,<key>,<name>``: the file and path boxes of the LOADGEN dialogs, one field per key (blank =
    the default; a commit of a blank box removes the entry)."""
    name = "LGFILE"

    def keys(self) -> List[str]:
        _, _, LC = _cmds()
        return list(LC.FILE_KEYS)

    def read(self, model) -> Dict[str, Any]:
        out = {}
        for k in self.keys():
            e = model.options.entry("LGFILE", k)
            out[k] = (e.arg(2) or "") if e is not None else ""
        return out

    def parse(self, vals: Dict[str, Any], P: _Parser) -> Dict[str, Any]:
        out = {}
        for k in self.keys():
            v = "" if vals.get(k) is None else str(vals.get(k)).strip()
            if "\n" in v or "\r" in v:
                P.bad(f"{k}: one line expected")
            out[k] = v
        return out

    def emit(self, model, new: Dict[str, Any], old: Dict[str, Any]) -> List[str]:
        return [join_command("LGFILE", [k, new[k]], text_last=True) for k in self.keys() if new[k] != old[k]]


class LgnodeFamily(RecordFamily):
    """``LGNODE,<kind>,<nodes>`` (kinds D, M, A): one node-list field per kind (``1, 3-6 10``).  A changed list
    is written as WRITE writes it: ``LGNODE,<kind>,0`` (clear), then the nodes, 20 tokens per line."""
    name = "LGNODE"

    def read(self, model) -> Dict[str, Any]:
        _, _, LC = _cmds()
        from .dialogs import _ids_text
        lists = LC.lgnode_lists(model)
        return {k: _ids_text(lists.get(k, [])) for k in LC.NODE_KINDS}

    def parse(self, vals: Dict[str, Any], P: _Parser) -> Dict[str, Any]:
        _, _, LC = _cmds()
        from .dialogs import parse_ids
        out = {}
        for k in LC.NODE_KINDS:
            try:
                ids = sorted(set(parse_ids(vals.get(k) or "")))
            except ValueError as exc:
                P.bad(f"{k} nodes: {exc}")
                ids = []
            if any(n < 1 for n in ids):
                P.bad(f"{k} nodes: node numbers must be >= 1")
            out[k] = ids
        return out

    def emit(self, model, new: Dict[str, Any], old: Dict[str, Any]) -> List[str]:
        out: List[str] = []
        for k in new:
            if new[k] == old[k]:
                continue
            if old[k]:
                out.append(f"LGNODE,{k},0")
            toks = compress_ids(new[k])
            for i in range(0, len(toks), 20):
                out.append(",".join(["LGNODE", k] + toks[i:i + 20]))
        return out


class LgtimeFamily(RecordFamily):
    """``LGTIME,<crit>,...``: the critical times of the static loads as dialog fields (criterion, number of
    peaks, separation, moment point, node and DOF, or the list of times / steps).  Default ``V,1``."""
    name = "LGTIME"
    FIELDS = ("crit", "n", "tsep", "x0", "y0", "z0", "node", "dof", "times")

    def read(self, model) -> Dict[str, Any]:
        rec = model.options.record("LGTIME")
        toks = [t or "" for t in rec.to_tokens()] if rec is not None else ["V", "1"]
        crit = (toks[0] if toks else "V").upper() or "V"
        vals: Dict[str, Any] = {"crit": crit, "n": 1, "tsep": 0.0, "x0": "", "y0": "", "z0": "", "node": "", "dof": 1,
                                "times": ""}

        def num(i: int, default: Any) -> Any:
            if i >= len(toks) or not toks[i]:
                return default
            try:
                return parse_float(toks[i])
            except NumberError:
                return toks[i]
        if crit in ("TIME", "STEP"):
            vals["times"] = " ".join(t for t in toks[1:] if t)
        else:
            n = num(1, 1)
            vals["n"] = int(n) if isinstance(n, float) and n == int(n) else n
            vals["tsep"] = num(2, 0.0)
            if crit.startswith("M"):
                vals["x0"], vals["y0"], vals["z0"] = num(3, ""), num(4, ""), num(5, "")
            elif crit in ("ACC", "DISP"):
                node, dof = num(3, ""), num(4, 1)
                vals["node"] = int(node) if isinstance(node, float) and node == int(node) else node
                vals["dof"] = int(dof) if isinstance(dof, float) and dof == int(dof) else dof
        return vals

    def parse(self, vals: Dict[str, Any], P: _Parser) -> Dict[str, Any]:
        """The canonical token list of the values (or an empty list when they are invalid)."""
        _, _, LC = _cmds()
        from .dialogs import parse_numbers
        crit = _tok(vals.get("crit")).upper() or "V"
        if crit not in LC.CRITERIA:
            P.bad(f"criterion '{crit}' is not one of {' '.join(LC.CRITERIA)}")
            return {"tokens": []}
        toks = [crit]
        if crit in ("TIME", "STEP"):
            try:
                ts = parse_numbers(vals.get("times") or "", integer=(crit == "STEP"))
            except ValueError as exc:
                P.bad(f"{crit.lower()}s: {exc}")
                return {"tokens": []}
            if not ts:
                P.bad(f"give at least one {'time' if crit == 'TIME' else 'time step'}")
            if crit == "TIME" and any(t < 0 for t in ts) or crit == "STEP" and any(t < 1 for t in ts):
                P.bad("times must be >= 0" if crit == "TIME" else "time steps are >= 1 (1-based)")
            return {"tokens": toks + [_tok(t) for t in ts]}
        n = P.num(vals.get("n"), "number of peaks", integer=True, default=1)
        tsep = P.num(vals.get("tsep"), "separation", default=0.0)
        if n is not None and n < 1:
            P.bad("the number of peaks must be >= 1")
        if tsep is not None and tsep < 0:
            P.bad("the separation must be >= 0")
        toks += [_tok(n), _tok(tsep)]
        if crit.startswith("M"):
            xyz = [P.num(vals.get(k), f"moment point {k}") for k in ("x0", "y0", "z0")]
            if any(v is not None for v in xyz):
                if any(v is None for v in xyz):
                    P.bad("give all three coordinates of the moment point (or none: centroid of the interface)")
                toks += [_tok(v) for v in xyz]
        elif crit in ("ACC", "DISP"):
            node = P.num(vals.get("node"), "node", integer=True)
            dof = P.num(vals.get("dof"), "DOF", integer=True)
            if node is None or node < 1 or dof is None or not 1 <= dof <= 6:
                P.bad(f"{crit} needs a node >= 1 and a DOF 1..6")
            toks += [_tok(node), _tok(dof)]
        return {"tokens": toks}

    def emit(self, model, new: Dict[str, Any], old: Dict[str, Any]) -> List[str]:
        return [join_command("LGTIME", new["tokens"])]


class LgmapFamily(RecordFamily):
    """``LGMAP,<mode>,<tol>,<file>`` (node numbering of the ANSYS model).  IDENTITY is written as ``LGMAP``
    alone (the command removes the record); a tolerance 0 or blank is the default (1e-6 of the model size)."""
    name = "LGMAP"

    def read(self, model) -> Dict[str, Any]:
        rec = model.options.record("LGMAP")
        if rec is None:
            return {"mode": "IDENTITY", "tol": "", "file": ""}
        tol = rec.number(2)
        return {"mode": (rec.arg(1) or "IDENTITY").upper(), "tol": "" if tol is None else tol, "file": rec.arg(3) or ""}

    def parse(self, vals: Dict[str, Any], P: _Parser) -> Dict[str, Any]:
        _, _, LC = _cmds()
        mode = P.word(vals.get("mode"), "numbering mode", {m: m for m in LC.MAPMODES}, "IDENTITY")
        if mode in (None, "IDENTITY"):
            return {"tokens": ["IDENTITY"] if mode else []}
        f = "" if vals.get("file") is None else str(vals.get("file")).strip()
        if not f:
            P.bad(f"{mode} needs the numbering file")
        tol = P.num(vals.get("tol"), "COORD tolerance", default=0.0) if mode == "COORD" else 0.0
        if tol is not None and tol < 0:
            P.bad("the COORD tolerance must be >= 0")
        return {"tokens": [mode, _tok(tol) if tol else "", f]}

    def emit(self, model, new: Dict[str, Any], old: Dict[str, Any]) -> List[str]:
        if new["tokens"] == ["IDENTITY"]:
            return ["LGMAP"]
        return [join_command("LGMAP", new["tokens"], text_last=True)]


# ======================================================================================
# Table families: %NAME (payload "xtables")
# ======================================================================================
class TableFamily:
    """Indexed entries edited as table rows keyed by their number: removed rows are deleted with
    ``delete``, new and changed rows written as complete command lines."""
    name = ""
    kind = "table"
    what = "entry"
    #: (column = record field, integer, default when blank); the first column is the key
    columns: Tuple[Tuple[str, bool, Any], ...] = ()

    def delete(self, num: int) -> str:                   # pragma: no cover - abstract
        raise NotImplementedError

    def read(self, model) -> List[Dict[str, Any]]:
        _cmds()
        out = []
        for k, rec in entries(model, self.name):
            if not isinstance(k, int):
                continue
            out.append({c: rec.get(c) for c, _, _ in self.columns})
        return out

    def _rows(self, rows: Any, P: _Parser) -> Dict[int, Tuple[Any, ...]]:
        out: Dict[int, Tuple[Any, ...]] = {}
        if not isinstance(rows, list):
            P.bad("the table must be a list of rows")
            return out
        for i, r in enumerate(rows, start=1):
            if not isinstance(r, dict):
                P.bad(f"row {i} is not an object")
                continue
            vals = tuple(P.num(r.get(c), f"row {i} <{c}>", integer=integer, default=d) for c, integer, d in self.columns)
            key = vals[0]
            if key is None:
                if _tok(r.get(self.columns[0][0])) == "":
                    P.bad(f"row {i}: the {self.what} number is required")
                continue
            if key < 1:
                P.bad(f"row {i}: {self.what} numbers start at 1 (got {key})")
                continue
            if key in out:
                P.bad(f"{self.what} {key} is given twice")
                continue
            out[key] = vals
        return out

    def lines(self, model, rows: Any, errors: List[str]) -> List[str]:
        n0 = len(errors)
        new = self._rows(rows, _Parser(self.name, errors))
        if len(errors) > n0:
            return []
        old = self._rows(self.read(model), _Parser(self.name, []))
        out = [self.delete(k) for k in sorted(old) if k not in new]
        out += [join_command(self.name, [_tok(v) for v in new[k]]) for k in sorted(new) if old.get(k) != new[k]]
        return out


class PanelTable(TableFamily):
    """``P,<num>,<group>,<bbc>,<disp>,<force>`` wall panels; ``PDEL,<num>`` deletes one."""
    name = "P"
    what = "panel"
    columns = (("num", True, None), ("group", True, 0), ("bbc", True, 0), ("disp", True, 1), ("force", True, 1))

    def delete(self, num: int) -> str:
        return f"PDEL,{num}"


class SpringTable(TableFamily):
    """``S,<num>,<group>,<elem>,<bbc>,<disp>,<force>`` nonlinear springs; ``DELSPR,<num>`` deletes one."""
    name = "S"
    what = "spring"
    columns = (("num", True, None), ("group", True, 0), ("elem", True, 0), ("bbc", True, 0), ("disp", True, 1),
               ("force", True, 4))

    def delete(self, num: int) -> str:
        return f"DELSPR,{num}"


class NlslayerTable(TableFamily):
    """``NLSLAYER,<Num>,[curvefit],[B],[S],[refStrain],[Vis]`` (SOIL-NON layer data); ``DELNLS,<n>,<n>,1``."""
    name = "NLSLAYER"
    what = "layer"
    columns = (("num", True, None), ("curvefit", True, 0), ("b", False, 0.0), ("s", False, 0.0),
               ("refstrain", False, 0.0), ("vis", False, 0.0))

    def delete(self, num: int) -> str:
        return f"DELNLS,{num},{num},1"


# ======================================================================================
# Keyed family: %BBC[bbc].field (payload "xindexed")
# ======================================================================================
class BbcFamily:
    """Backbone curves (BBCI, BBCX, BBCY; the BBC grid of the NONLINEAR tab).  Entry ``{type, yield, points:
    [{x, y}, ...]}`` per curve number; ``null`` deletes the curve.  A changed curve is redefined completely
    (``DELBBC`` + ``BBCI`` + ``BBCX`` + ``BBCY``), so no yield-number synchronisation warning is printed and the
    stored curve is exactly the one of the dialog; a new curve without points is ``BBCI`` only."""
    name = "BBC"
    kind = "keyed"
    DEFAULT = {"type": 0, "yield": "", "points": []}

    def read(self, model) -> Dict[str, Dict[str, Any]]:
        NC, _, _ = _cmds()
        out: Dict[str, Dict[str, Any]] = {}
        for num in NC.curve_numbers(model):
            cv = NC.curve(model, num)
            xs, ys = list(cv["x"]), list(cv["y"])
            pts = [{"x": xs[i] if i < len(xs) else "", "y": ys[i] if i < len(ys) else ""}
                   for i in range(max(len(xs), len(ys)))]
            out[str(num)] = {"type": int(cv["type"] or 0), "yield": int(cv["yield_"]) if cv["yield_"] else "",
                             "points": pts}
        return out

    def _entry(self, num: int, ent: Any, P: _Parser):
        if not isinstance(ent, dict):
            P.bad(f"curve {num}: the dialog values must be an object")
            return None
        n0 = len(P.errors)
        typ = P.num(ent.get("type"), f"curve {num} Type", integer=True, default=0)
        if typ is not None and typ not in (0, 1, 2, 3, 4):
            P.bad(f"curve {num}: Type must be 1 CMS, 2 CMB, 3 TAK or 4 GMR (0 = not set)")
        yld = P.num(ent.get("yield"), f"curve {num} Yield Num.", integer=True)
        pts: List[Tuple[float, float]] = []
        rows = ent.get("points") or []
        if not isinstance(rows, list):
            P.bad(f"curve {num}: the points must be a list")
            rows = []
        for i, r in enumerate(rows, start=1):
            r = r if isinstance(r, dict) else {}
            x = P.num(r.get("x"), f"curve {num} point {i} X")
            y = P.num(r.get("y"), f"curve {num} point {i} Y")
            if (x is None or y is None) and len(P.errors) == n0:
                P.bad(f"curve {num} point {i}: X and Y are required (Delete the row to remove the point)")
            pts.append((x, y))
        if yld is not None and (yld < 1 or (pts and yld > len(pts))):
            P.bad(f"curve {num}: Yield Num. {yld} must lie in 1..{len(pts) or 'number of points'}")
        if len(P.errors) > n0:
            return None
        return typ or 0, yld, tuple(pts)

    def lines(self, model, posted: Any, errors: List[str]) -> List[str]:
        P = _Parser("BBC", errors)
        if not isinstance(posted, dict):
            P.bad("the curves must be an object {number: curve}")
            return []
        shown = self.read(model)
        out: List[str] = []
        n0 = len(errors)

        def key_of(k: Any) -> Tuple[int, str]:
            return (_int0(k, 0), str(k))
        for key in sorted(posted, key=key_of):
            try:
                num, exact = parse_int(str(key))
            except NumberError:
                num, exact = 0, False
            if not exact or num < 1:
                P.bad(f"curve number '{key}' must be an integer >= 1")
                continue
            old = shown.get(str(num))
            ent = posted[key]
            if ent is None:
                if old is not None:
                    out.append(f"DELBBC,{num}")
                continue
            new = self._entry(num, ent, P)
            if new is None:
                continue
            if old is not None and self._entry(num, old, _Parser("BBC", [])) == new:
                continue
            typ, yld, pts = new
            if old is not None:
                out.append(f"DELBBC,{num}")
            yt = "" if yld is None else str(yld)
            out.append(join_command("BBCI", [str(num), yt, str(typ) if typ else ""]))
            if pts:
                out.append(join_command("BBCX", [str(num), str(len(pts)), yt] + [_tok(x) for x, _ in pts]))
                out.append(join_command("BBCY", [str(num), str(len(pts)), yt] + [_tok(y) for _, y in pts]))
        if any(len(ln) > MAX_LINE for ln in out):
            P.bad(f"a curve has too many points for one BBCX/BBCY line ({MAX_LINE} characters)")
        return [] if len(errors) > n0 else out


# ======================================================================================
# Registry, values, commit
# ======================================================================================
FAMILIES: Dict[str, Any] = {f.name: f for f in (
    EqlFamily(), NlsoilFamily(), LoadgenFamily(), LoadgenDynFamily(), LgfileFamily(), LgnodeFamily(),
    LgtimeFamily(), LgmapFamily(), LgoptFamily(), PanelTable(), SpringTable(), NlslayerTable(), BbcFamily())}
PAYLOAD_KEY = {"record": "xrecords", "table": "xtables", "keyed": "xindexed"}


def values(model, names: Sequence[str]) -> Dict[str, Any]:
    """The dialog values of the families ``names`` (``GET /api/options/<NAME>``)."""
    out: Dict[str, Any] = {"xrecords": {}, "xtables": {}, "xindexed": {}, "xindexed_defaults": {}}
    for n in names:
        fam = FAMILIES[n]
        out[PAYLOAD_KEY[fam.kind]][n] = fam.read(model)
        if fam.kind == "keyed":
            out["xindexed_defaults"][n] = dict(fam.DEFAULT)
    return out


def context(model) -> Dict[str, Any]:
    """Facts the NONLINEAR and SOIL tabs show next to the families."""
    return {"eql_opts_blank": EqlFamily.opts_blank(model), "eql_opts_derived": EqlFamily.derived(model),
            "beams": len(model.options.entries("B")),
            "spro_layers": sorted(int(k) for k, _ in model.options.entries("SPRO") if isinstance(k, int))}


def commit(interp, payload: Dict[str, Any], names: Sequence[str], errors: List[str]) -> List[str]:
    """Command lines of the family parts of a dialog commit, in the order of ``names``.

    Unknown family names are refused; when every value is readable, the lines are dry-run
    (:func:`dry_run`) and the error messages of their commands refuse the commit (UI-06)."""
    model = interp.model
    lines: List[str] = []
    n0 = len(errors)
    for key in PAYLOAD_KEY.values():
        part = payload.get(key) or {}
        if not isinstance(part, dict):
            errors.append(f"{key} must be an object")
            continue
        for n in part:
            fam = FAMILIES.get(str(n).upper())
            if fam is None or fam.name not in names or PAYLOAD_KEY[fam.kind] != key:
                errors.append(f"{n} is not a{'n indexed' if key == 'xindexed' else ''} "
                              f"{'table' if key == 'xtables' else 'record'} of this dialog")
    if len(errors) > n0:
        return []
    for n in names:
        fam = FAMILIES[n]
        part = {str(k).upper(): v for k, v in (payload.get(PAYLOAD_KEY[fam.kind]) or {}).items()}
        if n in part:
            lines += fam.lines(model, part[n], errors)
    if lines and not errors:                   # every value readable (also the other parts of the commit)
        errors += dry_run(interp, lines)
    return lines


def dry_run(interp, lines: Sequence[str]) -> List[str]:
    """Error messages of ``lines`` executed in a scratch interpreter holding a copy of the active model (the
    commands validate their own arguments; nothing of the session changes)."""
    from ..prep import Interpreter
    from ..prep.messages import Kind, MessageSink
    sink = MessageSink(keep=2000)
    current: List[str] = []
    sink.subscribe(lambda m: current.append(m.text) if m.kind == Kind.ERROR else None)
    sb = Interpreter(cwd=interp.cwd, sink=sink)
    sb.models = {interp.active_model: interp.model.copy()}
    sb.active_model = interp.active_model
    out: List[str] = []
    for ln in lines:
        current.clear()
        if not sb.execute(ln):
            out += list(current) or [f"{ln}: refused"]
    return out
