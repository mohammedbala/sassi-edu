"""Command handler modules of the interpreter.

Every module of this package is imported by :func:`sassi.prep.registry.load_commands`; a module
registers its handlers with :func:`sassi.prep.registry.command`.  See ``sassi/prep/README.md``.

This file also provides small helpers shared by the handler modules (range selection of
spec 08 section 1.4).
"""
from __future__ import annotations

from collections.abc import Mapping, Set
from typing import Iterable, Iterator, List, Optional


class UnionIds(Set):
    """Read-only union of the ids (keys) of several mappings or sets, without copying them.

    Membership is O(1) (one dict/set lookup per part), so a range command on "nodes or nodes with
    masses" (MUNITS) does not rebuild a set of every id at each call.  Iteration and ``len`` build
    the union on demand (O(N)); :func:`select_ids` avoids them for short ranges by using
    :meth:`size_bound`.
    """

    def __init__(self, *parts):
        self._parts = parts

    def __contains__(self, i) -> bool:
        return any(i in p for p in self._parts)

    def _union(self) -> set:
        out: set = set()
        for p in self._parts:
            out.update(p)
        return out

    def __iter__(self) -> Iterator[int]:
        return iter(sorted(self._union()))

    def __len__(self) -> int:
        return len(self._union())

    def size_bound(self) -> int:
        """Upper bound of ``len(self)`` in O(1) (the sum of the part sizes)."""
        return sum(len(p) for p in self._parts)


def _size(ids) -> int:
    return ids.size_bound() if isinstance(ids, UnionIds) else len(ids)


def select_ids(ids: Iterable[int], a: Optional[int], b: Optional[int] = None, inc: Optional[int] = None) -> List[int]:
    """Existing ids of the range ``a, a+inc, ... <= b`` in ascending order (spec 08 section 1.4).

    ``b`` defaults to ``a``; ``inc <= 0`` (or None) means 1; ``b < a`` swaps them; ids of the range
    that are not defined are skipped silently.  ``a = None`` selects all ids.
    """
    if a is None:
        return sorted(ids)
    if b is None:
        b = a
    if inc is None or inc <= 0:
        inc = 1
    if b < a:
        a, b = b, a
    n_range = (b - a) // inc + 1
    if isinstance(ids, (Mapping, Set)) and n_range <= _size(ids):
        # short range: O(1) membership tests (dicts, sets) instead of sorting every id
        return [i for i in range(a, b + 1, inc) if i in ids]
    return [i for i in sorted(ids) if a <= i <= b and (i - a) % inc == 0]


def range_args(c, k: int, ids: Iterable[int], all_default: bool = False):
    """Read ``<first>,[<last>],[<inc>]`` at arguments k, k+1, k+2 and select existing ids.

    With ``all_default`` an omitted first argument selects every id (listing commands);
    otherwise the first argument is required.
    """
    if not isinstance(ids, (Mapping, Set)):
        ids = set(ids)
    if not c.given(k):
        if all_default:
            if not ids:
                return []
            last = c.int(k + 1, default=max(ids))
            return select_ids(ids, min(ids), last, c.int(k + 2, default=1))
        c.int(k, required=True)
        a = 0
    else:
        a = c.int(k)
    b = c.int(k + 1, default=a)
    return select_ids(ids, a, b, c.int(k + 2, default=1))
