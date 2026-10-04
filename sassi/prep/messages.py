"""Message sink of the command interpreter (ARCHITECTURE section 8, requirements section 5.6).

Every message the interpreter or a command handler produces goes through a :class:`MessageSink`
with one of the six classes of the ACS SASSI Command History (manual 6.6.3):

=============  =========  =============================================================
class          colour     typical content
=============  =========  =============================================================
ECHO           black      the command line as submitted
CONFIRM        blue       "output confirmation": what a command did
COMMENT        green      ``*`` comment lines (shown also while a ``.pre`` file runs)
INFO           purple     model information, listings (cannot be hidden)
WARNING        amber      warnings (processing continues)
ERROR          red        errors (the failing command is skipped, processing continues)
=============  =========  =============================================================

View > Command Display filters echo, confirmation, comments and warnings/errors; information
cannot be hidden.  While INP or a macro runs, echo and confirmation are suppressed by the
interpreter (rule L16), so they never reach the sink.

The console prints messages through :meth:`MessageSink.attach_stream`; the GUI subscribes with
:meth:`MessageSink.subscribe`; tests read :attr:`MessageSink.messages`.
"""
from __future__ import annotations

import enum
import sys
from collections import deque
from dataclasses import dataclass
from typing import Callable, Deque, Dict, List, Optional, TextIO


class Kind(enum.IntEnum):
    ECHO = 1
    CONFIRM = 2
    COMMENT = 3
    INFO = 4
    WARNING = 5
    ERROR = 6


#: default colours of the Command History (requirements section 5.6)
COLOURS = {Kind.ECHO: "black", Kind.CONFIRM: "blue", Kind.ERROR: "red", Kind.COMMENT: "green",
           Kind.INFO: "purple", Kind.WARNING: "amber"}

PREFIX = {Kind.ECHO: "> ", Kind.CONFIRM: "  ", Kind.COMMENT: "", Kind.INFO: "",
          Kind.WARNING: "*** WARNING: ", Kind.ERROR: "*** ERROR: "}


@dataclass
class Message:
    kind: Kind
    text: str
    source: str = ""      # 'keyboard', '<file>:<line>', 'macro NAME:<line>'
    command: str = ""     # canonical command name, when produced by a command

    def format(self) -> str:
        return PREFIX[self.kind] + self.text


Listener = Callable[[Message], None]


class MessageSink:
    """Collects, filters and dispatches interpreter messages.

    Parameters
    ----------
    keep:
        Maximum number of messages kept in :attr:`messages` (oldest dropped); ``None`` keeps all.
    stream:
        Optional text stream that receives every *shown* message (console front end).
    """

    #: View > Command Display toggles (INFO is always shown)
    FILTER_GROUPS = {"echo": (Kind.ECHO,), "confirm": (Kind.CONFIRM,), "comments": (Kind.COMMENT,),
                     "warnerr": (Kind.WARNING, Kind.ERROR)}

    def __init__(self, keep: Optional[int] = 100000, stream: Optional[TextIO] = None,
                 err_stream: Optional[TextIO] = None):
        self.messages: Deque[Message] = deque(maxlen=keep)
        self.show: Dict[Kind, bool] = {k: True for k in Kind}
        self.counts: Dict[Kind, int] = {k: 0 for k in Kind}
        self._listeners: List[Listener] = []
        self._stream = stream
        self._err_stream = err_stream

    # -- filters
    def set_filter(self, group: str, visible: bool) -> None:
        """Toggle one View > Command Display group: 'echo', 'confirm', 'comments', 'warnerr'."""
        for k in self.FILTER_GROUPS[group]:
            self.show[k] = bool(visible)

    def set_visible(self, kind: Kind, visible: bool) -> None:
        if kind == Kind.INFO:
            return           # information cannot be hidden (requirements section 5.6)
        self.show[kind] = bool(visible)

    def visible(self, kind: Kind) -> bool:
        return kind == Kind.INFO or self.show.get(kind, True)

    # -- dispatch
    def subscribe(self, fn: Listener) -> None:
        self._listeners.append(fn)

    def unsubscribe(self, fn: Listener) -> None:
        if fn in self._listeners:
            self._listeners.remove(fn)

    def attach_stream(self, stream: Optional[TextIO], err_stream: Optional[TextIO] = None) -> None:
        self._stream = stream
        self._err_stream = err_stream

    def emit(self, kind: Kind, text: str, source: str = "", command: str = "") -> Message:
        msg = Message(kind, text, source, command)
        self.counts[kind] += 1
        self.messages.append(msg)
        if self.visible(kind):
            if self._stream is not None:
                out = self._err_stream if (kind >= Kind.WARNING and self._err_stream is not None) else self._stream
                for line in msg.format().splitlines() or [""]:
                    out.write(line + "\n")
                out.flush()
            for fn in list(self._listeners):
                fn(msg)
        return msg

    # -- inspection helpers (tests, GUI)
    def texts(self, kind: Optional[Kind] = None) -> List[str]:
        return [m.text for m in self.messages if kind is None or m.kind == kind]

    def clear(self) -> None:
        self.messages.clear()
        self.counts = {k: 0 for k in Kind}

    def last(self, kind: Optional[Kind] = None) -> Optional[Message]:
        for m in reversed(self.messages):
            if kind is None or m.kind == kind:
                return m
        return None


def console_sink(quiet: bool = False) -> MessageSink:
    """Sink printing to stdout (warnings and errors to stderr).

    ``quiet`` prints only warnings and errors (batch runs, ``sassi run --quiet``); the messages
    are still collected, so nothing is lost for later inspection.
    """
    if not quiet:
        return MessageSink(stream=sys.stdout, err_stream=sys.stderr)
    s = MessageSink()

    def _print(msg: Message) -> None:
        if msg.kind >= Kind.WARNING:
            sys.stderr.write(msg.format() + "\n")
            sys.stderr.flush()
    s.subscribe(_print)
    return s
