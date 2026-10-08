"""The SASSI-EDU command interpreter (ARCHITECTURE sections 1.3 and 8; requirements section 3.1).

One :class:`Interpreter` serves the console, ``.pre`` files (INP), macros, FOREACH loops and the
GUI (UI-01, L17).  It holds the numbered models (ACTM), the working directory (CD/MDL), the
variables and macros of the programming commands, and a :class:`~sassi.prep.messages.MessageSink`.

Processing of one line (rule L13, D-PAR-10)::

    1. macro $k$ placeholders            (done by MACRO before the line reaches execute())
    2. variables  @NAME, @NAME[i], @NAME++ ...  and  # inside FOREACH bodies, left to right,
       each evaluated once (D-PAR-15/16); skipped for FOREACH, whose body is substituted afresh
       at every iteration
    3. tokenisation (sassi.prep.lexer)
    4. dispatch to the registered handler (sassi.prep.registry)

Errors never stop a ``.pre`` file (L12, D-PAR-09): an unknown name prints ``<X> Command not
found``, a failing command prints its error, and processing continues with the next line.  At the
end of every INP file a summary (commands, warnings, errors) is printed, and when the outermost
file ends ``INPUT FILE REACHED EOF, INPUT SWITCHED TO KEYBOARD`` (L16).  While INP or macro input
runs, command echo and confirmations are suppressed (L16).
"""
from __future__ import annotations

import os
import re
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Union

import numpy as np

from ..io import library as LIB
from ..model import SSIModel
from ..model.values import NumberError, parse_float, parse_int
from .lexer import MAX_LINE, LexError, Lexed, is_comment, parse_id_list, split_args, split_head
from .messages import Kind, MessageSink
from .registry import LOAD_ERRORS, CommandError, CommandReported, CommandSpec, load_commands, lookup

INP_DEPTH = 32        # L14 / D-PAR-11
MACRO_DEPTH = 64
EOF_MESSAGE = "INPUT FILE REACHED EOF, INPUT SWITCHED TO KEYBOARD"
BOM = "\ufeff"


def read_command_lines(path: Union[str, Path]) -> List[str]:
    """Lines of a command, macro or variable file (INP, LOADMACRO, LOADVAR).

    Spec 04 section 4.1: existing ``.pre`` files load unchanged.  Files saved by Windows editors
    often start with a UTF-8 byte-order mark; it is removed (``utf-8-sig``) so that it does not
    become part of the first command name.  Undecodable bytes are replaced, never fatal (L12).
    Raises :class:`OSError` when the file cannot be read.
    """
    text = Path(path).read_text(encoding="utf-8-sig", errors="replace")
    return text.splitlines()


# --------------------------------------------------------------------------------------
# Session objects
# --------------------------------------------------------------------------------------
@dataclass
class Variable:
    """A variable of the programming commands: a list of strings and an integer counter (spec 04 6.1)."""
    name: str
    items: List[str] = field(default_factory=list)
    counter: int = 0


@dataclass
class Macro:
    """A loaded macro (LOADMACRO): the cached lines of its file (spec 04 section 5)."""
    name: str
    file: str
    lines: List[str] = field(default_factory=list)


@dataclass
class LoopFrame:
    var: str          # lower-case variable name of the FOREACH
    index: int = 0    # 1-based iteration index ('#')


@dataclass
class InputFrame:
    path: Optional[Path]
    line: int = 0
    kind: str = "file"     # 'file', 'text' or 'macro'
    name: str = ""


@dataclass
class RunSummary:
    """Counts of one INP file (or text) run, nested files included."""
    source: str
    lines: int = 0
    commands: int = 0
    warnings: int = 0
    errors: int = 0
    ok: bool = True        # False when the file could not be read

    def text(self) -> str:
        return (f"INP {self.source}: {self.lines} lines, {self.commands} commands, "
                f"{self.warnings} warnings, {self.errors} errors")


# --------------------------------------------------------------------------------------
# Handler context
# --------------------------------------------------------------------------------------
class Call:
    """What a command handler receives: the parsed arguments plus access to the session.

    Arguments are addressed **1-based** as in the manual (``c.int(1)`` is the first value after
    the command name).  Blank or omitted fields return the ``default`` (rule L5); a *required*
    numeric field that is blank is read as 0 with a warning (D-PAR-05); a non-numeric value in
    a numeric field raises :class:`CommandError` (the command is skipped).
    """

    def __init__(self, interp: "Interpreter", spec: CommandSpec, typed: str, lexed: Lexed, rest: str, line: str):
        self.interp = interp
        self.spec = spec
        self.name = spec.name
        self.typed = typed
        self.lexed = lexed
        self.tokens: List[str] = lexed.tokens
        self.rest = rest
        self.line = line

    # ------------------------------------------------------------------ session access
    @property
    def model(self) -> SSIModel:
        return self.interp.model

    @property
    def nargs(self) -> int:
        """Number of argument tokens up to the last non-blank one."""
        n = len(self.tokens)
        while n and not self.tokens[n - 1].strip():
            n -= 1
        return n

    # ------------------------------------------------------------------ raw tokens
    def raw(self, k: int) -> str:
        return self.tokens[k - 1].strip() if 1 <= k <= len(self.tokens) else ""

    def given(self, k: int) -> bool:
        return bool(self.raw(k))

    def _label(self, k: int, what: Optional[str]) -> str:
        return f"argument {k}" + (f" <{what}>" if what else "")

    def float(self, k: int, default: Optional[float] = None, required: bool = False,
              what: Optional[str] = None) -> Optional[float]:
        t = self.raw(k)
        if not t:
            if required:
                self.warn(f"{self._label(k, what)} missing, 0 used")
                return 0.0
            return default
        try:
            return parse_float(t)
        except NumberError:
            raise CommandError(f"{self._label(k, what)}: '{t}' is not a number") from None

    def int(self, k: int, default: Optional[int] = None, required: bool = False,
            what: Optional[str] = None) -> Optional[int]:
        t = self.raw(k)
        if not t:
            if required:
                self.warn(f"{self._label(k, what)} missing, 0 used")
                return 0
            return default
        try:
            v, exact = parse_int(t)
        except NumberError:
            raise CommandError(f"{self._label(k, what)}: '{t}' is not an integer") from None
        if not exact:
            self.warn(f"{self._label(k, what)}: {t} rounded to {v}")
        return v

    def str(self, k: int, default: str = "") -> str:
        t = self.raw(k)
        return t if t else default

    def word(self, k: int, default: str = "") -> str:
        """Keyword argument, upper-cased (keywords are case-insensitive)."""
        return self.str(k, default).upper()

    def floats(self, k_from: int, k_to: Optional[int] = None, default: float = 0.0) -> List[float]:
        k_to = self.nargs if k_to is None else k_to
        return [self.float(k, default) for k in range(k_from, k_to + 1)]

    def ints(self, k_from: int, k_to: Optional[int] = None, default: int = 0) -> List[int]:
        k_to = self.nargs if k_to is None else k_to
        return [self.int(k, default) for k in range(k_from, k_to + 1)]

    def id_list(self, k_from: int) -> List[int]:
        """Node/element list from argument ``k_from`` to the end (rule L8)."""
        try:
            return parse_id_list([t for t in self.tokens[k_from - 1:] if t.strip()])
        except LexError as exc:
            raise CommandError(str(exc)) from None

    def text(self, k: int) -> str:
        """Rest-of-line text from argument ``k`` (for commands registered with ``text_from``)."""
        return self.raw(k)

    # ------------------------------------------------------------------ messages
    def _prefixed(self, text: str) -> str:
        """Messages start with the command name (added unless the text already starts with it)."""
        n = len(self.name)
        if text[:n].upper() == self.name and (len(text) == n or not text[n].isalnum()):
            return text
        return f"{self.name}: {text}"

    def warn(self, text: str) -> None:
        self.interp.emit(Kind.WARNING, self._prefixed(text), command=self.name)

    def error(self, text: str) -> None:
        """Report a non-fatal error (the command continues)."""
        self.interp.emit(Kind.ERROR, self._prefixed(text), command=self.name)

    def info(self, text: str) -> None:
        self.interp.emit(Kind.INFO, text, command=self.name)

    def confirm(self, text: str) -> None:
        self.interp.emit(Kind.CONFIRM, text, command=self.name)

    def fail(self, text: str) -> None:
        raise CommandError(text)

    # ------------------------------------------------------------------ paths
    def input_path(self, p: str) -> Path:
        return self.interp.resolve_path(p, must_exist=True)

    def output_path(self, p: str) -> Path:
        return self.interp.output_path(p)


# --------------------------------------------------------------------------------------
# Interpreter
# --------------------------------------------------------------------------------------
_OP_RE = re.compile(r"([+-])(\d+)")
_SET_RE = re.compile(r"=([+-]?\d+)")
_HASH_RE = re.compile(r"#(\+\+|--|[+-]\d+)?")


class Interpreter:
    """Command interpreter holding the models, variables, macros and the message sink.

    Parameters
    ----------
    cwd:
        Initial working directory (default: the process working directory).
    sink:
        Message sink (default: a collecting :class:`MessageSink`).
    seed:
        Seed of the RND/ADDRND generator (numpy PCG64, D-GEN-08); ``None`` = OS entropy until
        RNDSEED is given.
    """

    def __init__(self, cwd: Optional[Union[str, Path]] = None, sink: Optional[MessageSink] = None,
                 seed: Optional[int] = None):
        load_commands()
        self.sink = sink if sink is not None else MessageSink()
        self.models: Dict[int, SSIModel] = {0: SSIModel()}
        self.active_model = 0
        self.cwd = Path(cwd).resolve() if cwd is not None else Path(os.getcwd()).resolve()
        self.variables: Dict[str, Variable] = {}
        self.macros: Dict[str, Macro] = {}
        self.loops: List[LoopFrame] = []
        self.inputs: List[InputFrame] = []
        self.rng = np.random.Generator(np.random.PCG64(seed))
        #: session-global state of other packages (cuts D-MDL-13, plot lines, ...)
        self.session: Dict[str, Any] = {}
        #: commands accepted from the keyboard (GUI replay, L17)
        self.history: List[str] = []
        #: GUI progress hook: progress(current_line, total_lines, source, command_line) -- called before the
        #: line of an INP / macro file is executed
        self.progress: Optional[Callable[[int, int, str, str], None]] = None
        #: GUI hook for module runs in this interpreter (RUN<MODULE> inside an INP file):
        #: module_progress(module, fraction, text), the progress the module reports (ANALYS: frequency k/n ...)
        self.module_progress: Optional[Callable[[str, float, str], None]] = None
        #: module_step(module, key, data): the computation step a module enters, with its sizes (GUI)
        self.module_step: Optional[Callable[[str, str, Dict[str, Any]], None]] = None
        #: GUI hooks of the Run view: input_frame(event, frame, lines) when an INP / macro file starts ("start",
        #: before its first line) and ends ("end"); line_done(line, total, source, command_line, ok) after each
        #: of its lines.  Display aids: they never change or stop the run.
        self.input_frame: Optional[Callable[[str, "InputFrame", Sequence[str]], None]] = None
        self.line_done: Optional[Callable[[int, int, str, str, bool], None]] = None
        self.write_options: Dict[str, Any] = {"mdl": False, "afwr": False}
        self._summaries: List[RunSummary] = []
        self._direct = False
        self.last_ok = True
        for mod, tb in LOAD_ERRORS:
            self.emit(Kind.WARNING, f"command module {mod} failed to load:\n{tb.strip().splitlines()[-1]}")

    # ------------------------------------------------------------------ models
    @property
    def model(self) -> SSIModel:
        return self.models[self.active_model]

    def activate(self, number: int) -> SSIModel:
        """ACTM: make model ``number`` active, creating an empty model when absent."""
        if number not in self.models:
            self.models[number] = SSIModel()
        self.active_model = number
        return self.models[number]

    def new_model(self) -> int:
        """Model > New: create and activate an empty model with the lowest unused number (D-MDL-01).

        The GUI submits the equivalent ``ACTM,<n>`` command text (L17); this helper only finds n.
        """
        n = 0
        while n in self.models:
            n += 1
        self.activate(n)
        return n

    # ------------------------------------------------------------------ state
    @property
    def batch(self) -> bool:
        """True while INP, macro or FOREACH input runs (echo and confirmations suppressed, L16)."""
        return bool(self.inputs) or bool(self.loops)

    def _source(self) -> str:
        if not self.inputs:
            return "keyboard"
        f = self.inputs[-1]
        return f"{f.name}:{f.line}"

    # ------------------------------------------------------------------ messages
    def emit(self, kind: Kind, text: str, command: str = "") -> None:
        if kind in (Kind.ECHO, Kind.CONFIRM) and self.batch:
            return
        for s in self._summaries:
            if kind == Kind.WARNING:
                s.warnings += 1
            elif kind == Kind.ERROR:
                s.errors += 1
        self.sink.emit(kind, text, source=self._source(), command=command)

    # ------------------------------------------------------------------ substitution (L13, D-PAR-15)
    def _match_variable(self, text: str, pos: int) -> Optional[str]:
        """Longest defined variable name starting at ``text[pos]`` (case-insensitive)."""
        best = None
        for name in self.variables:
            if (best is None or len(name) > len(best)) and text[pos:pos + len(name)].lower() == name:
                best = name
        return best

    def _loop_index(self, var: Optional[str] = None) -> int:
        if not self.loops:
            raise CommandError("'#' used outside a FOREACH loop")
        if var is not None:
            for fr in reversed(self.loops):
                if fr.var == var:
                    return fr.index
        return self.loops[-1].index

    def _index_value(self, expr: str, var: str) -> int:
        e = expr.strip()
        if e.startswith("#"):
            m = _HASH_RE.fullmatch(e)
            if not m:
                raise CommandError(f"bad index expression [{expr}]")
            idx = self._loop_index(var)
            op = m.group(1)
            if op == "++":
                idx += 1
            elif op == "--":
                idx -= 1
            elif op:
                idx += int(op)
            return idx
        if e.startswith("@"):
            e = self.substitute(e)
        try:
            v, exact = parse_int(e)
        except NumberError:
            raise CommandError(f"bad index expression [{expr}]") from None
        if not exact:
            raise CommandError(f"index [{expr}] is not an integer")
        return v

    def substitute(self, line: str) -> str:
        """Replace ``@`` variable expressions and (inside FOREACH) ``#`` -- rule L13 step 2.

        Grammar (spec 04 section 6.2, D-PAR-15/16), with ``X`` the longest defined variable name:
        ``@X`` counter; ``@X+k``, ``@X-k``, ``@X++``, ``@X--`` change the counter first and give
        the new value; ``@X=k`` sets it; ``@X[i]`` is the i-th item (1-based; ``i`` may be ``#``,
        ``#+k`` or ``@Y``); ``#`` is the innermost loop index and ``#+k`` adds k without mutation.
        An unknown name is left literally with a warning; an index out of range is an error.  ``@`` followed
        by the name of a built-in library file (``@rg160h_030g.acc``) is a file name, never a variable, and is
        left as it is (D-W5-02).
        """
        if "@" not in line and not ("#" in line and self.loops):
            return line
        out: List[str] = []
        i, n = 0, len(line)
        while i < n:
            ch = line[i]
            if ch == "@":
                nlib = LIB.starts_with_name(line[i + 1:])
                if nlib:                        # a built-in input file name (@rg160h_030g.acc, D-W5-02)
                    out.append(line[i:i + 1 + nlib])
                    i += 1 + nlib
                    continue
                name = self._match_variable(line, i + 1)
                if name is None:
                    m = re.match(r"[A-Za-z_][A-Za-z0-9_]*", line[i + 1:])
                    if m:   # '@' not followed by a name is plain data (e.g. a title)
                        self.emit(Kind.WARNING, f"variable '{m.group(0)}' is not defined; "
                                                f"'@{m.group(0)}' left as is")
                    out.append("@")
                    i += 1
                    continue
                var = self.variables[name]
                j = i + 1 + len(name)
                tail = line[j:]
                if tail.startswith("["):
                    depth, k = 0, 0
                    for k, c2 in enumerate(tail):
                        if c2 == "[":
                            depth += 1
                        elif c2 == "]":
                            depth -= 1
                            if depth == 0:
                                break
                    if depth != 0:
                        raise CommandError(f"unterminated index in '@{var.name}['")
                    idx = self._index_value(tail[1:k], name)
                    if not 1 <= idx <= len(var.items):
                        raise CommandError(f"@{var.name}[{idx}]: index out of range 1..{len(var.items)} (D-PAR-16)")
                    out.append(var.items[idx - 1])
                    i = j + k + 1
                    continue
                if tail.startswith("++"):
                    var.counter += 1
                    j += 2
                elif tail.startswith("--"):
                    var.counter -= 1
                    j += 2
                else:
                    m = _OP_RE.match(tail)
                    if m:
                        var.counter += int(m.group(2)) * (1 if m.group(1) == "+" else -1)
                        j += m.end()
                    else:
                        m = _SET_RE.match(tail)
                        if m:
                            var.counter = int(m.group(1))
                            j += m.end()
                out.append(str(var.counter))
                i = j
                continue
            if ch == "#" and self.loops:
                m = _HASH_RE.match(line, i)
                idx = self._loop_index()
                op = m.group(1) if m else None
                if op == "++":
                    idx += 1
                elif op == "--":
                    idx -= 1
                elif op:
                    idx += int(op)
                out.append(str(idx))
                i = m.end() if m else i + 1
                continue
            out.append(ch)
            i += 1
        return "".join(out)

    # ------------------------------------------------------------------ execution
    def execute(self, line: str) -> bool:
        """Execute one command line; returns False when it failed (the session continues)."""
        text = line.rstrip("\r\n")
        if not text.strip():
            return True
        if is_comment(text):
            self.emit(Kind.COMMENT, text.strip())
            return True
        interactive = not self.batch
        if self._direct:          # a line of an INP file (loop iterations and macro lines not counted)
            self._direct = False
            for s in self._summaries:
                s.commands += 1
        if len(text) > MAX_LINE:
            self.emit(Kind.ERROR, f"command line longer than {MAX_LINE} characters (L14)")
            return False
        head, _ = split_head(text)
        spec = lookup(head) if ("@" not in head and "#" not in head) else None
        if spec is not None and spec.raw:
            sub = text
        else:
            try:
                sub = self.substitute(text)
            except CommandError as exc:
                if interactive:
                    self.emit(Kind.ECHO, text.strip())
                self.emit(Kind.ERROR, f"{head.upper()}: {exc}")
                return False
            if len(sub) > MAX_LINE:
                self.emit(Kind.ERROR, f"command line longer than {MAX_LINE} characters after substitution (L14)")
                return False
        name, rest = split_head(sub)
        spec = lookup(name)
        if interactive:
            self.emit(Kind.ECHO, text.strip())
        if spec is None:
            self.emit(Kind.ERROR, f"{name.upper()} Command not found")
            return False
        try:
            lexed = split_args(rest, spec.text_from, spec.paren) if not spec.raw else Lexed(tokens=[rest])
        except LexError as exc:
            self.emit(Kind.ERROR, f"{spec.name}: {exc}", command=spec.name)
            return False
        call = Call(self, spec, name, lexed, rest, sub)
        if spec.max_args is not None and call.nargs > spec.max_args:
            # the documented argument count is a hard limit: the extra arguments are really
            # dropped, so the handler cannot use what the warning says is ignored
            call.warn(f"arguments after argument {spec.max_args} ignored")
            lexed.tokens = lexed.tokens[:spec.max_args]
            lexed.quoted = lexed.quoted[:spec.max_args]
            call.tokens = lexed.tokens
        ok = True
        try:
            spec.handler(call)
        except CommandReported:           # the handler has printed its own error
            ok = False
        except CommandError as exc:
            call.error(str(exc))
            ok = False
        except RecursionError:
            call.error("recursion limit reached")
            ok = False
        except Exception as exc:  # a handler bug must not stop the session (L12)
            call.error(f"internal error: {exc!r}")
            self.sink.emit(Kind.INFO, traceback.format_exc(), source=self._source(), command=spec.name)
            ok = False
        if interactive and ok:
            self.history.append(text.strip())
        self.last_ok = ok
        return ok

    def execute_lines(self, lines: Sequence[str]) -> int:
        """Execute several lines as typed at the keyboard; returns the number of failures."""
        return sum(0 if self.execute(ln) else 1 for ln in lines)

    def _run_lines(self, lines: Sequence[str], frame: InputFrame, summary: RunSummary) -> RunSummary:
        if sum(1 for f in self.inputs if f.kind != "macro") >= INP_DEPTH:
            raise CommandError(f"INP nesting deeper than {INP_DEPTH} (L14)")
        self._summaries.append(summary)
        self.inputs.append(frame)
        self._hook(self.input_frame, "start", frame, lines)
        try:
            total = len(lines)
            for i, line in enumerate(lines, start=1):
                frame.line = i
                if self.progress is not None:
                    self.progress(i, total, summary.source, line)
                self._direct = True
                ok = self.execute(line)
                self._direct = False
                self._hook(self.line_done, i, total, summary.source, line, ok)
        finally:
            self.inputs.pop()
            self._summaries.pop()
            self._hook(self.input_frame, "end", frame, lines)
        summary.lines = len(lines)
        self.emit(Kind.INFO, summary.text())
        if not self.inputs:
            self.emit(Kind.INFO, EOF_MESSAGE)
        return summary

    @staticmethod
    def _hook(fn: Optional[Callable[..., None]], *args: Any) -> None:
        """Call a GUI display hook; a failing hook never changes or stops the run."""
        if fn is None:
            return
        try:
            fn(*args)
        except Exception:                 # noqa: BLE001 -- a display aid
            pass

    def run_file(self, path: Union[str, Path], resolve: bool = True) -> RunSummary:
        """INP: execute the commands of a ``.pre`` file (nested INP allowed, depth 32).

        A relative ``path`` is resolved by :meth:`resolve_path` (L15).  An unreadable file is the
        only condition that aborts a file (L12); it is reported as an error and the returned
        summary has ``ok = False`` (the INP command then fails).  A UTF-8 byte-order mark is
        removed (:func:`read_command_lines`).
        """
        try:
            p = self.resolve_path(str(path), must_exist=True) if resolve else Path(path)
        except CommandError as exc:
            self.emit(Kind.ERROR, f"INP: {exc}")
            return RunSummary(str(path), errors=1, ok=False)
        try:
            lines = read_command_lines(p)
        except OSError as exc:
            self.emit(Kind.ERROR, f"INP: cannot read {p}: {exc.strerror or exc}")
            return RunSummary(str(p), errors=1, ok=False)
        frame = InputFrame(p, kind="file", name=p.name)
        return self._run_lines(lines, frame, RunSummary(p.name))

    def run_text(self, text: str, name: str = "<text>") -> RunSummary:
        """Execute a block of command text with INP semantics (batch echo policy, summary)."""
        if text.startswith(BOM):          # text pasted from a file with a byte-order mark
            text = text[1:]
        frame = InputFrame(None, kind="text", name=name)
        return self._run_lines(text.splitlines(), frame, RunSummary(name))

    def run_macro(self, name: str, args: Sequence[str]) -> None:
        """MACRO: execute a loaded macro, replacing ``$k$`` by argument k as plain text (spec 04 5.3)."""
        m = self.macros.get(name.upper())
        if m is None:
            raise CommandError(f"macro {name.upper()} is not loaded (LOADMACRO)")
        depth = sum(1 for f in self.inputs if f.kind == "macro")
        if depth >= MACRO_DEPTH:
            raise CommandError(f"macro nesting deeper than {MACRO_DEPTH} (L14)")
        frame = InputFrame(Path(m.file), kind="macro", name=f"macro {m.name}")
        missing: List[int] = []

        def _sub(mo: "re.Match") -> str:
            k = int(mo.group(1))
            if 1 <= k <= len(args):
                return args[k - 1]
            if k not in missing:
                missing.append(k)
            return ""

        self.inputs.append(frame)
        try:
            for i, line in enumerate(m.lines, start=1):
                frame.line = i
                missing.clear()
                text = re.sub(r"\$([1-9][0-9]*)\$", _sub, line)
                for k in missing:
                    self.emit(Kind.WARNING, f"MACRO {m.name}: argument {k} missing, empty string used (D-PAR-13)")
                if len(text) > MAX_LINE:
                    self.emit(Kind.ERROR, f"MACRO {m.name} line {i}: longer than {MAX_LINE} characters")
                    continue
                self.execute(text)
        finally:
            self.inputs.pop()

    # ------------------------------------------------------------------ paths (L15, D-PAR-12)
    def _candidates(self, p: Path) -> List[Path]:
        cands: List[Path] = []
        if self.model.path:
            cands.append(Path(self.model.path) / p)
        cands.append(self.cwd / p)
        for f in reversed(self.inputs):
            if f.path is not None and f.kind in ("file", "macro"):
                cands.append(f.path.parent / p)
                break
        return cands

    def resolve_path(self, name: str, must_exist: bool = True) -> Path:
        """Resolve an input path: absolute; active model path; working directory; calling file's folder.

        Windows separators (``.\\Node-Macro.pre`` in the manual examples) are accepted.  ``@name`` is a
        file of the built-in input library (:mod:`sassi.io.library`, D-W5-01).
        """
        if LIB.is_library_name(name):
            p = LIB.library_path(name)
            if p is not None:
                return p
            if must_exist:
                raise CommandError(f"built-in file {name.strip()} not found (LIBRARY lists the built-in inputs: "
                                   f"{', '.join(LIB.names())})")
            return LIB.LIBRARY_DIR / LIB.bare(name)
        raw = os.path.expanduser(name.strip())
        variants = [raw]
        if "\\" in raw and os.sep == "/":
            variants.append(raw.replace("\\", "/"))
        tried: List[Path] = []
        for v in variants:
            p = Path(v)
            cands = [p] if p.is_absolute() else self._candidates(p)
            for c in cands:
                tried.append(c)
                if c.exists():
                    return c
        if must_exist:
            raise CommandError(f"file {name.strip()} not found (searched {', '.join(str(t) for t in tried)})")
        p = Path(variants[-1])
        return p if p.is_absolute() else self._candidates(p)[0]

    def output_path(self, name: str) -> Path:
        """Path of a file to write: absolute, else relative to the model path (if MDL set), else the CWD.
        A built-in ``@`` name is refused: the library is read-only (D-W5-01)."""
        if LIB.is_library_name(name):
            raise CommandError(f"{name.strip()}: built-in library files are read-only; give another file name")
        raw = os.path.expanduser(name.strip())
        if "\\" in raw and os.sep == "/":
            raw = raw.replace("\\", "/")
        p = Path(raw)
        if p.is_absolute():
            return p
        base = Path(self.model.path) if self.model.path else self.cwd
        return base / p
