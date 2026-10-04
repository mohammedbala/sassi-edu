"""Group, element and element-attribute commands (manual section 9.4; spec 08 sections 4-5;
requirements 3.4.D).

* One element type per group; element numbers start at 1 in each group (gaps are compressed with
  ECOMPR before AFWRITE, CHECK Error 6).
* ``E`` gives a new element the active material ``MACT`` and property ``RACT`` (D-MDL-17: both
  start at 1 and are global), ETYPE 0, EINT 0, THICK 0 and no releases.
* MSET stores the raw index; whether it points to the M or the L table (excavated SOLID/PLANE) is
  decided at CHECK/AFWRITE so that a later ETYPE change reinterprets it (D-MDL-10).
* EGEN increments every node number (the BEAMS/GENERAL K node too) and copies the attributes of
  the pattern element (D-MDL-09).
* GDEL drops the EOUT requests of the deleted groups; ECOMPR remaps them (D-MDL-11).
"""
from __future__ import annotations

from ...conventions import ELEMENT_TYPE_CODES, ELEMENT_TYPE_NAMES
from ...model import fmt_num
from ...model.entities import MAX_ELEMENT_NODES, MIN_ELEMENT_NODES, Element, Group, type_name
from ...model.values import NumberError, parse_int
from ..registry import command
from . import range_args, select_ids

SOLID, BEAMS, SHELL, PLANE, TSHELL, SPRING, GENERAL = 1, 2, 3, 4, 5, 7, 9


def parse_group_type(c, token: str) -> int:
    """Group type given as a number (1 2 3 4 5 7 9) or a name (SOLID, BEAMS/BEAM, ...; D-MDL-16)."""
    t = token.strip().upper()
    try:
        v = parse_int(t)[0]
    except NumberError:
        v = ELEMENT_TYPE_CODES.get(t)
        if v is None:
            c.fail(f"unknown group type {token.strip()} (SOLID BEAMS SHELL PLANE TSHELL SPRING GENERAL)")
    if v not in ELEMENT_TYPE_NAMES:
        c.fail(f"group type {v} is not 1, 2, 3, 4, 5, 7 or 9")
    return v


def _group(c) -> Group:
    g = c.model.active_group
    if g is None:
        c.fail("no active group (use GROUP)")
    return g


# ======================================================================================
# Groups
# ======================================================================================
@command("GROUP", max_args=2)
def cmd_group(c):
    """GROUP,<ng>,<type>: create a group of the given type, or activate an existing group."""
    m = c.model
    ng = c.int(1, required=True, what="ng")
    if ng < 1:
        c.fail("group numbers are >= 1")
    g = m.groups.get(ng)
    if not c.given(2):
        if g is None:
            c.fail(f"group {ng} does not exist: the group type is required")
        m.group_active = ng
        c.confirm(f"group {ng} ({g.type_name}) active")
        return
    t = parse_group_type(c, c.raw(2))
    if g is not None and g.type != t:
        c.fail(f"group {ng} is of type {g.type_name}; use MTYPE to change it")
    if g is None:
        m.groups[ng] = Group(ng, t)
        c.confirm(f"group {ng} ({type_name(t)}) created and active")
    else:
        c.confirm(f"group {ng} ({type_name(t)}) active")
    m.group_active = ng


@command("MTYPE", max_args=2)
def cmd_mtype(c):
    """MTYPE,[<gr>],<type>: change the type of a group (default the active group)."""
    m = c.model
    if c.nargs <= 1:
        gid, tok = m.group_active, c.raw(1)
    else:
        gid = c.int(1, default=m.group_active)
        tok = c.raw(2)
    if gid is None or gid not in m.groups:
        c.fail("group not defined" if gid is not None else "no active group")
    if not tok:
        c.fail("type required")
    t = parse_group_type(c, tok)
    g = m.groups[gid]
    bad = [e.id for e in g.elements.values()
           if len(e.nodes) > MAX_ELEMENT_NODES[t] or len(e.nodes) < MIN_ELEMENT_NODES[t]]
    g.type = t
    if bad:
        c.warn(f"{len(bad)} elements do not have the node count of {type_name(t)} (CHECK Error 7)")
    c.confirm(f"group {gid} is now {type_name(t)}")


@command("GTIT", text_from=2)
def cmd_gtit(c):
    """GTIT,[gr],<title>: group title (rest of the line); without a group number the active group."""
    m = c.model
    first = c.raw(1)
    if not first:                         # GTIT,,<title>
        gid, title = m.group_active, c.text(2)
    else:
        try:
            v, exact = parse_int(first)
        except NumberError:
            v, exact = None, False
        if exact:                         # GTIT,<gr>,<title>
            gid, title = v, c.text(2)
        else:                             # GTIT,<title>  (active group; the whole rest is the title)
            gid, title = m.group_active, c.rest.strip()
    if gid is None or gid not in m.groups:
        c.fail(f"group {gid} is not defined" if gid is not None else "no active group")
    if len(title) >= 2 and title[0] == '"' and title[-1] == '"':
        title = title[1:-1]
    m.groups[gid].title = title
    c.confirm(f"group {gid} title: {title}")


@command("GDEL", max_args=3)
def cmd_gdel(c):
    """GDEL,<g1>,[<g2>],[<inc>]: delete groups and their elements (their EOUT requests are dropped)."""
    m = c.model
    ids = range_args(c, 1, m.groups)
    for gid in ids:
        del m.groups[gid]
        if m.group_active == gid:
            m.group_active = None
    s = set(ids)
    n_req = sum(1 for r in m.eout if r.group in s)
    if n_req:
        m.eout = [r for r in m.eout if r.group not in s]
        c.warn(f"{n_req} EOUT requests of the deleted groups removed")
    c.confirm(f"{len(ids)} groups deleted")


@command("GLIST", max_args=3)
def cmd_glist(c):
    """GLIST,[<g1>],[<g2>],[<inc>]: list groups."""
    m = c.model
    ids = range_args(c, 1, m.groups, all_default=True)
    c.info(f"{'group':>6} {'type':>3} {'name':<8} {'elements':>9}  title")
    for gid in ids:
        g = m.groups[gid]
        mark = "*" if gid == m.group_active else " "
        c.info(f"{mark}{gid:>5} {g.type:>3} {g.type_name:<8} {len(g.elements):>9}  {g.title}")
    c.info(f"{len(ids)} groups listed")


# ======================================================================================
# Elements
# ======================================================================================
@command("E", max_args=9)
def cmd_e(c):
    """E,<ne>,<n1>,...,<n8>: element of the active group (node count by group type)."""
    m = c.model
    g = _group(c)
    ne = c.int(1, required=True, what="ne")
    if ne < 1:
        c.fail("element numbers start at 1")
    nodes = [c.int(k, default=0) for k in range(2, c.nargs + 1)]
    while nodes and nodes[-1] == 0:
        nodes.pop()
    if not nodes:
        c.fail("node numbers required")
    if len(nodes) > MAX_ELEMENT_NODES[g.type]:
        c.fail(f"{g.type_name} elements have at most {MAX_ELEMENT_NODES[g.type]} nodes")
    if any(n < 0 for n in nodes):
        c.fail("node numbers must be positive")
    over = ne in g.elements
    g.elements[ne] = Element(ne, nodes, mat=m.mact, prop=m.ract)
    if over:
        c.warn(f"element {ne} of group {g.id} redefined")
    c.confirm(f"element {ne} of group {g.id} ({g.type_name})")


@command("EGEN", max_args=6)
def cmd_egen(c):
    """EGEN,[<itim>],<ninc1>,<e1>,[<e2>],[<inc>],[<ee>]: copy an element pattern itim times.

    Generated element k of a pattern element has node numbers ``n + k*ninc1`` (non-zero nodes,
    K node included) and the pattern's attributes; numbering continues from ``ee`` (default
    largest element number + 1).
    """
    g = _group(c)
    itim = c.int(1, default=1)
    ninc = c.int(2, required=True, what="ninc1")
    e1 = c.int(3, required=True, what="e1")
    e2 = c.int(4, default=e1)
    inc = c.int(5, default=1)
    pattern = select_ids(g.elements, e1, e2, inc)
    if not pattern:
        c.fail(f"no elements {e1}..{e2} in group {g.id}")
    nxt = c.int(6, default=max(g.elements) + 1)
    if itim < 1:
        c.warn("itim < 1: nothing generated")
        return
    if nxt < 1:
        c.fail("element numbers start at 1 (<ee> < 1)")
    src = [g.elements[i] for i in pattern]
    # a non-zero pattern node must stay a valid node number (>= 1) in every copy; 0 would turn it
    # into a blank slot and a negative number is not a node.  Checked before anything is created.
    given = [n for p in src for n in p.nodes if n != 0]
    if given:
        low = min(given) + (itim if ninc < 0 else 1) * ninc
        if low < 1:
            c.fail(f"generated node numbers must be positive (node {low} would be referenced)")
    over = 0
    first = nxt
    for k in range(1, itim + 1):
        for p in src:
            nodes = [n + k * ninc if n != 0 else 0 for n in p.nodes]
            over += nxt in g.elements
            g.elements[nxt] = p.copy(new_id=nxt, nodes=nodes)
            nxt += 1
    if over:
        c.warn(f"{over} existing elements overwritten")
    c.confirm(f"elements {first}..{nxt - 1} generated in group {g.id}")


@command("EDEL", max_args=3)
def cmd_edel(c):
    """EDEL,<e1>,[<e2>],[<inc>]: delete elements of the active group (compress with ECOMPR)."""
    g = _group(c)
    ids = range_args(c, 1, g.elements)
    for i in ids:
        del g.elements[i]
    c.confirm(f"{len(ids)} elements deleted from group {g.id}")


@command("ECOMPR", max_args=0)
def cmd_ecompr(c):
    """ECOMPR: renumber the elements of the active group 1..N (order kept; EOUT remapped)."""
    m = c.model
    g = _group(c)
    old = sorted(g.elements)
    mapping = {o: k for k, o in enumerate(old, start=1)}
    g.elements = {mapping[o]: g.elements[o] for o in old}
    for k, e in g.elements.items():
        e.id = k
    dropped = 0
    for r in m.eout:
        if r.group == g.id:
            keep = [mapping[e] for e in r.elements if e in mapping]
            dropped += len(r.elements) - len(keep)
            r.elements = keep
    m.eout = [r for r in m.eout if r.elements]
    if dropped:
        c.warn(f"{dropped} EOUT references to undefined elements removed")
    c.confirm(f"group {g.id}: {len(old)} elements numbered 1..{len(old)}")


@command("ELIST", max_args=3)
def cmd_elist(c):
    """ELIST,[<e1>],[<e2>],[<inc>]: list the elements of the active group."""
    g = _group(c)
    ids = range_args(c, 1, g.elements, all_default=True)
    c.info(f"group {g.id} ({g.type_name}) {g.title}")
    c.info(f"{'elem':>7}  {'nodes':<48} {'mat':>4} {'prop':>4} {'etype':>5} {'eint':>4} {'thick':>10}  KI     KJ")
    for i in ids:
        e = g.elements[i]
        nodes = " ".join(str(n) for n in e.nodes)
        c.info(f"{i:>7}  {nodes:<48} {e.mat:>4} {e.prop:>4} {e.etype:>5} {e.eint:>4} {e.thick:>10.4g}  "
               f"{''.join(str(v) for v in e.ki)} {''.join(str(v) for v in e.kj)}")
    c.info(f"{len(ids)} elements listed")


# ======================================================================================
# Element attributes
# ======================================================================================
def _elements(c):
    g = _group(c)
    return g, [g.elements[i] for i in range_args(c, 1, g.elements)]


@command("ETYPE", max_args=4)
def cmd_etype(c):
    """ETYPE,<e1>,<e2>,[<inc>],<type>: 0 implicit (by ground elevation), 1 structure, 2 excavated soil
    (SOLID/PLANE) or embedded/buried shell (SHELL/TSHELL)."""
    g, els = _elements(c)
    t = c.int(4, required=True, what="type")
    if t not in (0, 1, 2):
        c.fail("<type> must be 0, 1 or 2")
    if g.type in (BEAMS, SPRING, GENERAL) and t:
        c.warn(f"ETYPE has no effect on {g.type_name} elements")
    for e in els:
        e.etype = t
    c.confirm(f"{len(els)} elements of group {g.id}: ETYPE {t}")


@command("EINT", max_args=4)
def cmd_eint(c):
    """EINT,<e1>,<e2>,[<inc>],<order>: SOLID 0/1/2 (2, 3, 4 Gauss points); TSHELL 0 reduced / 1 selective."""
    g, els = _elements(c)
    order = c.int(4, required=True, what="order")
    if g.type == SOLID:
        if order not in (0, 1, 2):
            c.fail("SOLID integration order must be 0, 1 or 2")
    elif g.type == TSHELL:
        if order not in (0, 1):
            c.fail("TSHELL integration must be 0 (reduced) or 1 (selective)")
    else:
        c.warn(f"EINT applies to SOLID and TSHELL groups only; ignored for {g.type_name}")
        return
    for e in els:
        e.eint = order
    c.confirm(f"{len(els)} elements of group {g.id}: EINT {order}")


@command("THICK", max_args=4)
def cmd_thick(c):
    """THICK,<e1>,<e2>,[<inc>],<thick>: shell thickness (SHELL and TSHELL groups)."""
    g, els = _elements(c)
    if g.type not in (SHELL, TSHELL):
        c.fail(f"THICK needs a SHELL or TSHELL group (group {g.id} is {g.type_name})")
    t = c.float(4, required=True, what="thick")
    if t <= 0:
        c.fail("thickness must be > 0")
    for e in els:
        e.thick = t
    c.confirm(f"{len(els)} elements of group {g.id}: thickness {fmt_num(t)}")


def _releases(c, end: str):
    g, els = _elements(c)
    if g.type != BEAMS:
        c.fail(f"{c.name} needs a BEAMS group (group {g.id} is {g.type_name})")
    k = [c.int(i, default=0) for i in range(4, 10)]
    if any(v not in (0, 1) for v in k):
        c.fail("release codes are 0 (connected) or 1 (released)")
    for e in els:
        setattr(e, end, list(k))
    c.confirm(f"{len(els)} beams of group {g.id}: {c.name} {''.join(str(v) for v in k)} (P1 P2 P3 M1 M2 M3)")


@command("KI", max_args=9)
def cmd_ki(c):
    """KI,<e1>,[<e2>],[<inc>],<k1>,...,<k6>: beam end releases at node I (1 = released)."""
    _releases(c, "ki")


@command("KJ", max_args=9)
def cmd_kj(c):
    """KJ,<e1>,[<e2>],[<inc>],<k1>,...,<k6>: beam end releases at node J (1 = released)."""
    _releases(c, "kj")


@command("MSET", max_args=4)
def cmd_mset(c):
    """MSET,<e1>,[<e2>],[<inc>],<index>: material (M) index, or soil-layer (L) index for excavated
    SOLID/PLANE elements (interpreted at CHECK/AFWRITE, D-MDL-10)."""
    g, els = _elements(c)
    idx = c.int(4, required=True, what="index")
    if g.type in (SPRING, GENERAL):
        c.warn(f"MSET is not used by {g.type_name} elements; ignored")
        return
    for e in els:
        e.mat = idx
    c.confirm(f"{len(els)} elements of group {g.id}: MSET {idx}")


@command("RSET", max_args=4)
def cmd_rset(c):
    """RSET,<e1>,[<e2>],[<inc>],<index>: property index (R for BEAMS, SC for SPRING, MX for GENERAL)."""
    g, els = _elements(c)
    idx = c.int(4, required=True, what="index")
    if g.type not in (BEAMS, SPRING, GENERAL):
        c.warn(f"RSET is not used by {g.type_name} elements; ignored")
        return
    for e in els:
        e.prop = idx
    c.confirm(f"{len(els)} elements of group {g.id}: RSET {idx}")


@command("MACT", max_args=1)
def cmd_mact(c):
    """MACT,<index>: active material / soil-layer index for new elements (all groups)."""
    c.model.mact = c.int(1, required=True, what="index")
    c.confirm(f"active material index {c.model.mact}")


@command("RACT", max_args=1)
def cmd_ract(c):
    """RACT,<index>: active real / spring / matrix property index for new elements (all groups)."""
    c.model.ract = c.int(1, required=True, what="index")
    c.confirm(f"active property index {c.model.ract}")
