"""Load and mass commands (manual section 9.5; spec 08 section 7; requirements 3.4.E).

* F / MM give, per global direction, a load **factor** and an **arrival time** applied to the one
  reference load history used by FORCE: ``P(w) = a F_ref(w) exp(-i w t0)`` (e^{+i w t}).
* MT / MR are nodal masses; MUNITS sets per node whether they are given as mass (0) or weight (1,
  divided by g when the decks are written).  The default for a node never given MUNITS is 1
  (weight units), D-MDL-08.
* MOPT ``<mass>`` / ``<force>`` act when MT/MR/F/MM execute (D-MDL-08): 1 overwrites an existing
  definition (default), 0 adds to it.  In add mode the factors of F/MM are added and, for the
  components that receive a non-zero new factor, the new arrival time replaces the old one with a
  warning when they differ (a sum of two delayed loads cannot be stored as one factor/time pair).
* "Last two defined" defaults of FSCALE/MSCALE/MTSCALE/MRSCALE follow the definition order
  (D-MDL-05); scaling factors of 0 (or blank) mean 1.
"""
from __future__ import annotations

from typing import Dict, List

from ...model import fmt_num
from ...model.entities import NodalLoad
from ...model.ssimodel import touch_history
from ..registry import command
from . import UnionIds, range_args, select_ids


def _node(c) -> int:
    n = c.int(1, required=True, what="n")
    if n < 1:
        c.fail("node numbers must be positive")
    return n


# ======================================================================================
# Forces and moments
# ======================================================================================
def _load(c, table: Dict[int, NodalLoad], hist: List[int], what: str) -> None:
    m = c.model
    n = _node(c)
    fac = [c.float(k, 0.0) for k in (2, 3, 4)]
    arr = [c.float(k, 0.0) for k in (5, 6, 7)]
    old = table.get(n)
    if old is None or m.mopt.get("force") == 1:
        table[n] = NodalLoad(n, fac, arr)
        verb = "redefined" if old is not None else "defined"
    else:
        for k in range(3):
            if fac[k] != 0.0:
                if old.factor[k] != 0.0 and old.arrival[k] != arr[k]:
                    c.warn(f"node {n} component {k + 1}: arrival time {fmt_num(old.arrival[k])} replaced by "
                           f"{fmt_num(arr[k])} (add mode)")
                old.arrival[k] = arr[k]
            old.factor[k] += fac[k]
        verb = "added"
    touch_history(hist, n)
    c.confirm(f"{what} at node {n} {verb}")


@command("F", max_args=7)
def cmd_f(c):
    """F,<n>,<fx>,<fy>,<fz>,<tx>,<ty>,<tz>: force factors and arrival times (global X, Y, Z)."""
    _load(c, c.model.forces, c.model.force_history, "force")


@command("MM", max_args=7)
def cmd_mm(c):
    """MM,<n>,<fxx>,<fyy>,<fzz>,<txx>,<tyy>,<tzz>: moment factors and arrival times (about X, Y, Z)."""
    _load(c, c.model.moments, c.model.moment_history, "moment")


def _delete(c, table: dict, hists, what: str) -> List[int]:
    ids = range_args(c, 1, table)
    for i in ids:
        del table[i]
    for h in hists:
        h.discard(ids)
    c.confirm(f"{len(ids)} {what} deleted")
    return ids


@command("FDEL", max_args=3)
def cmd_fdel(c):
    """FDEL,<n1>,[<n2>],[<inc>]: delete forces."""
    _delete(c, c.model.forces, [c.model.force_history], "forces")


@command("MMDEL", max_args=3)
def cmd_mmdel(c):
    """MMDEL,<n1>,[<n2>],[<inc>]: delete moments."""
    _delete(c, c.model.moments, [c.model.moment_history], "moments")


def _list_loads(c, table: Dict[int, NodalLoad], what: str) -> None:
    ids = range_args(c, 1, table, all_default=True)
    c.info(f"{'node':>8} {what + ' factors X, Y, Z':>42}   arrival times X, Y, Z")
    for i in ids:
        ld = table[i]
        c.info(f"{i:>8} " + " ".join(f"{v:>13.6g}" for v in ld.factor) + "   "
               + " ".join(f"{v:>10.4g}" for v in ld.arrival))
    c.info(f"{len(ids)} {what}s listed")


@command("FLIST", max_args=3)
def cmd_flist(c):
    """FLIST,[<n1>],[<n2>],[<inc>]: list forces."""
    _list_loads(c, c.model.forces, "force")


@command("MMLIST", max_args=3)
def cmd_mmlist(c):
    """MMLIST,[<n1>],[<n2>],[<inc>]: list moments."""
    _list_loads(c, c.model.moments, "moment")


def _scale_range(c, table: dict, hist: List[int]) -> List[int]:
    """``[<n1>],[<n2>],[<inc>]`` with the "last two defined" default (spec 08 section 7.5)."""
    if not c.given(1) and not c.given(2):
        if not hist:
            return []
        a, b = (hist[-2], hist[-1]) if len(hist) >= 2 else (hist[-1], hist[-1])
        a, b = min(a, b), max(a, b)
    else:
        a = c.int(1, default=c.int(2))
        b = c.int(2, default=a)
    return select_ids(table, a, b, c.int(3, default=1))


def _factors(c):
    return [c.float(k, 0.0) or 1.0 for k in (4, 5, 6)]


def _scale_loads(c, table: Dict[int, NodalLoad], hist: List[int], what: str) -> None:
    ids = _scale_range(c, table, hist)
    f = _factors(c)
    for i in ids:
        table[i].factor = [v * s for v, s in zip(table[i].factor, f)]
    c.confirm(f"{len(ids)} {what} scaled by ({', '.join(fmt_num(v) for v in f)})")


@command("FSCALE", max_args=6)
def cmd_fscale(c):
    """FSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]: scale force factors (0 -> 1; arrival times kept)."""
    _scale_loads(c, c.model.forces, c.model.force_history, "forces")


@command("MSCALE", max_args=6)
def cmd_mscale(c):
    """MSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]: scale moment factors (0 -> 1)."""
    _scale_loads(c, c.model.moments, c.model.moment_history, "moments")


# ======================================================================================
# Masses
# ======================================================================================
def _set_mass(m, table: Dict[int, List[float]], n: int, vals: List[float]) -> bool:
    """MOPT <mass>: 1 overwrite (default), 0 add; returns True when added to an existing mass."""
    if n in table and m.mopt.get("mass") == 0:
        table[n] = [a + b for a, b in zip(table[n], vals)]
        return True
    table[n] = list(vals)
    return False


def _mass(c, table: Dict[int, List[float]], hist: List[int], what: str) -> None:
    m = c.model
    n = _node(c)
    vals = [c.float(k, 0.0) for k in (2, 3, 4)]
    existed = n in table
    added = _set_mass(m, table, n, vals)
    touch_history(hist, n)
    units = "weight" if m.mass_unit(n) == 1 else "mass"
    c.confirm(f"{what} at node {n} {'added' if added else ('redefined' if existed else 'defined')} ({units} units)")


@command("MT", max_args=4)
def cmd_mt(c):
    """MT,<n>,<mx>,<my>,<mz>: translational masses in global X, Y, Z (units by MUNITS)."""
    _mass(c, c.model.tmass, c.model.tmass_history, "translational mass")


@command("MR", max_args=4)
def cmd_mr(c):
    """MR,<n>,<mxx>,<myy>,<mzz>: rotational masses about global X, Y, Z (units by MUNITS)."""
    _mass(c, c.model.rmass, c.model.rmass_history, "rotational mass")


def _massgen(c, table: Dict[int, List[float]], hist: List[int], what: str) -> None:
    """MTGEN/MRGEN: copy pattern masses to n + k*ninc with increments k*(mx, my, mz)."""
    m = c.model
    itim = c.int(1, default=1)
    n1 = c.int(3, required=True, what="n1")
    n2 = c.int(4, default=n1)
    if n2 < n1:
        n1, n2 = n2, n1
    ninc = c.int(2, default=0) or (n2 - n1 + 1)
    inc = c.int(5, default=1)
    incr = [c.float(k, 0.0) for k in (6, 7, 8)]
    pattern = select_ids(table, n1, n2, inc)
    if not pattern:
        c.fail(f"no {what}es defined on nodes {n1}..{n2}")
    if itim < 1:
        c.warn("itim < 1: nothing generated")
        return
    low = pattern[0] + (itim if ninc < 0 else 1) * ninc
    if low < 1:
        c.fail(f"generated node numbers must be positive (node {low} would receive a mass)")
    base = {n: list(table[n]) for n in pattern}
    units = {n: m.mass_unit(n) for n in pattern}
    for k in range(1, itim + 1):
        for n in pattern:
            tgt = n + k * ninc
            _set_mass(m, table, tgt, [b + k * d for b, d in zip(base[n], incr)])
            _set_units(m, tgt, units[n])
            touch_history(hist, tgt)
    c.confirm(f"{itim * len(pattern)} {what}es generated")


@command("MTGEN", max_args=8)
def cmd_mtgen(c):
    """MTGEN,[<itim>],[<ninc>],<n1>,<n2>,[<inc>],[<mx>],[<my>],[<mz>]: generate translational masses."""
    _massgen(c, c.model.tmass, c.model.tmass_history, "translational mass")


@command("MRGEN", max_args=8)
def cmd_mrgen(c):
    """MRGEN,[<itim>],[<ninc>],<n1>,<n2>,[<inc>],[<mxx>],[<myy>],[<mzz>]: generate rotational masses."""
    _massgen(c, c.model.rmass, c.model.rmass_history, "rotational mass")


def _drop_orphan_units(m) -> None:
    for i in [i for i in m.mass_units if i not in m.nodes and i not in m.tmass and i not in m.rmass]:
        del m.mass_units[i]


@command("MTDEL", max_args=3)
def cmd_mtdel(c):
    """MTDEL,<n1>,[<n2>],[<inc>]: delete translational masses."""
    _delete(c, c.model.tmass, [c.model.tmass_history], "translational masses")
    _drop_orphan_units(c.model)


@command("MRDEL", max_args=3)
def cmd_mrdel(c):
    """MRDEL,<n1>,[<n2>],[<inc>]: delete rotational masses."""
    _delete(c, c.model.rmass, [c.model.rmass_history], "rotational masses")
    _drop_orphan_units(c.model)


def _scale_mass(c, table: Dict[int, List[float]], hist: List[int], what: str) -> None:
    ids = _scale_range(c, table, hist)
    f = _factors(c)
    for i in ids:
        table[i] = [v * s for v, s in zip(table[i], f)]
    c.confirm(f"{len(ids)} {what} scaled by ({', '.join(fmt_num(v) for v in f)})")


@command("MTSCALE", max_args=6)
def cmd_mtscale(c):
    """MTSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]: scale translational masses (0 -> 1)."""
    _scale_mass(c, c.model.tmass, c.model.tmass_history, "translational masses")


@command("MRSCALE", max_args=6)
def cmd_mrscale(c):
    """MRSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]: scale rotational masses (0 -> 1)."""
    _scale_mass(c, c.model.rmass, c.model.rmass_history, "rotational masses")


@command("MTLIST", max_args=3)
def cmd_mtlist(c):
    """MTLIST,[<n1>],[<n2>],[<inc>]: list translational and rotational masses with units and mass values."""
    m = c.model
    g = m.gravity
    ids = range_args(c, 1, UnionIds(m.tmass, m.rmass), all_default=True)
    try:
        masses = m.nodal_masses()
    except ValueError as exc:          # weight units with gravity <= 0 (CHECK Error 1)
        masses = None
        c.warn(f"masses in mass units not available ({exc}); entered values listed only")
    c.info(f"{'node':>8} {'units':>6} {'mx':>12} {'my':>12} {'mz':>12} {'mxx':>12} {'myy':>12} {'mzz':>12}")
    for i in ids:
        u = m.mass_unit(i)
        raw = list(m.tmass.get(i, [0.0] * 3)) + list(m.rmass.get(i, [0.0] * 3))
        c.info(f"{i:>8} {('weight' if u else 'mass'):>6} " + " ".join(f"{v:>12.5g}" for v in raw))
        if u == 1 and masses is not None:
            c.info(f"{'':>8} {'mass':>6} " + " ".join(f"{v:>12.5g}" for v in masses[i]) + f"   (g = {fmt_num(g)})")
    c.info(f"{len(ids)} nodes with masses listed")


def _set_units(m, n: int, u: int) -> None:
    if u == m.MUNITS_DEFAULT:
        m.mass_units.pop(n, None)
    else:
        m.mass_units[n] = u


@command("MUNITS", max_args=4)
def cmd_munits(c):
    """MUNITS,<n1>,[<n2>],[<step>],<units>: nodal mass units, 0 mass / 1 weight (default 1, D-MDL-08)."""
    m = c.model
    # nodes and nodes with masses are eligible; a view with O(1) membership instead of a new set
    # per command keeps INP of a WRITE output (one MUNITS line per run of units) linear
    ids = range_args(c, 1, UnionIds(m.nodes, m.tmass, m.rmass))
    u = c.int(4, required=True, what="units")
    if u not in (0, 1):
        c.fail("<units> must be 0 (mass) or 1 (weight)")
    for i in ids:
        _set_units(m, i, u)
    c.confirm(f"{len(ids)} nodes: masses in {'weight' if u else 'mass'} units")
