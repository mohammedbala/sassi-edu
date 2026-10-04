"""Event log of a GUI session: what the browser polls (``GET /api/events?since=N``).

Every message of the interpreter (the six classes of the Command History, requirements 5.6), every
plot event of :class:`sassi.plotting.state.PlotState` (open, update, close, activate, dialog ...),
progress reports (INP current line / total lines, module runs) and job state changes are appended
with a sequence number.  The browser long-polls for the events after the last number it has seen,
so the Command History fills while a long command (INP, RUNxxx inside a ``.pre``) is still running
(UI-04: the GUI stays responsive).

The log keeps the last ``keep`` events; a client that asks for older ones receives ``reset`` and
re-reads the session state.

Listeners (:meth:`EventLog.subscribe`) are called after every new event: the browser version of the
GUI (:mod:`sassi.web.bridge`) cannot long-poll, so it pushes the new events to the page while a long
command or a module run keeps its Python engine busy.
"""
from __future__ import annotations

import threading
import time
from collections import deque
from typing import Any, Callable, Deque, Dict, List, Tuple


class EventLog:
    """Thread-safe, bounded, sequence-numbered event list with long polling."""

    def __init__(self, keep: int = 20000):
        self._events: Deque[Dict[str, Any]] = deque(maxlen=keep)
        self._seq = 0
        self._cond = threading.Condition()
        self._listeners: List[Callable[[Dict[str, Any]], None]] = []

    @property
    def last(self) -> int:
        """Sequence number of the newest event (0 when empty)."""
        return self._seq

    def subscribe(self, fn: Callable[[Dict[str, Any]], None]) -> None:
        """Call ``fn(event)`` after every new event (outside the lock)."""
        self._listeners.append(fn)

    def push(self, type_: str, **data: Any) -> int:
        """Append an event ``{'seq', 'type', 'time', ...data}``; returns its sequence number."""
        with self._cond:
            self._seq += 1
            ev = {"seq": self._seq, "type": type_, "time": time.time()}
            ev.update(data)
            self._events.append(ev)
            self._cond.notify_all()
            seq = self._seq
        for fn in self._listeners:
            fn(ev)
        return seq

    def since(self, seq: int, timeout: float = 0.0, limit: int = 5000) -> Tuple[List[Dict[str, Any]], int, bool]:
        """Events with a sequence number above ``seq``.

        Waits up to ``timeout`` seconds when there is none (long poll).  Returns
        ``(events, last_seq, reset)``; ``reset`` is True when events after ``seq`` were already
        dropped from the bounded log (the client then reloads the state).
        """
        deadline = time.time() + max(0.0, timeout)
        with self._cond:
            while self._seq <= seq:
                left = deadline - time.time()
                if left <= 0:
                    break
                self._cond.wait(left)
            if not self._events:
                return [], self._seq, False
            first = self._events[0]["seq"]
            reset = seq + 1 < first
            out = [e for e in self._events if e["seq"] > seq]
            if len(out) > limit:
                out = out[:limit]
            last = out[-1]["seq"] if out else self._seq
            return out, last, reset
