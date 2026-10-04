"""Node and coordinate-system commands (manual section 9.3; spec 08 sections 2-3; requirements 3.4.C).

Nodes are stored in the coordinate system that is active when they are defined; GLOBAL and
AFWRITE convert them to global coordinates (D-MDL-04).  The generation commands NGEN, FILL,
LMOVE, NMOVE, NSCALE work in the **active** system (ANSYS behaviour, spec 08 Q3): source
coordinates are expressed in the active system, the increments or factors are applied along the
active axes and the result is stored in the active system.  FILL and NMED give the same
geometric result in any Cartesian system.

"Last defined" defaults (FILL, NGEN) follow the definition order of nodes (D-MDL-05).  Generated
nodes do not inherit fixities or INT flags (spec 08 section 3.2).
"""
from __future__ import annotations

import numpy as np

from ...model import fmt_num
from ...model.entities import DOF_LABELS, DOF_NAMES, INT_CODES, INT_LETTERS, CoordSys
from ...model.geometry import loc_matrix, local_from_points
from ..registry import command
from . import range_args, select_ids


def _hist_defaults(m, k_n1: int, k_n2: int, c):
    """n1, n2 defaults: second-to-last and last defined node (spec 08 section 1.5)."""
    h = m.node_history
    d2 = h[-1] if h else None
    d1 = h[-2] if len(h) >= 2 else d2
    n1 = c.int(k_n1, default=d1)
    n2 = c.int(k_n2, default=d2)
    if n1 is None or n2 is None:
        c.fail("no nodes defined")
    return n1, n2


def _active(c) -> int:
    s = c.model.csys_active
    if s != 0 and s not in c.model.csys:
        c.fail(f"active coordinate system {s} is not defined")
    return s


def _check_new_ids(c, ids, what: str = "node") -> None:
    """Generated numbers must be valid ids (>= 1), as N and E require; checked before anything is
    created so that a failing generator leaves the model unchanged."""
    bad = sorted({i for i in ids if i < 1})
    if bad:
        c.fail(f"generated {what} numbers must be positive ({what} {bad[0]} would be created)")


def _set_coords(m, nid: int, xyz, s: int) -> None:
    """Replace the stored coordinates of an existing node without changing its definition order."""
    n = m.nodes[nid]
    n.x, n.y, n.z = (float(xyz[0]), float(xyz[1]), float(xyz[2]))
    n.csys = s


# ======================================================================================
# Definition and deletion
# ======================================================================================
@command("N", max_args=6)
def cmd_n(c):
    """N,<nd>,[<x>],[<y>],[<z>]: define a node in the active coordinate system (legacy fields 5-6 ignored)."""
    nd = c.int(1, required=True, what="nd")
    if nd <= 0:
        c.fail("node numbers must be positive")
    s = _active(c)
    xyz = (c.float(2, 0.0), c.float(3, 0.0), c.float(4, 0.0))
    existed = c.model.define_node(nd, xyz, s)
    c.confirm(f"node {nd} {'redefined' if existed else 'defined'} at ({fmt_num(xyz[0])}, {fmt_num(xyz[1])}, "
              f"{fmt_num(xyz[2])})" + (f" in system {s}" if s else ""))


@command("NDEL", max_args=3)
def cmd_ndel(c):
    """NDEL,<n1>,[<n2>],[<inc>]: delete nodes with their loads, masses, fixities and flags (D-MDL-11)."""
    m = c.model
    ids = range_args(c, 1, m.nodes)
    dangling = m.delete_nodes(ids)
    if dangling:
        c.warn(f"{dangling} element node references now point to deleted nodes (CHECK Error 41)")
    c.confirm(f"{len(ids)} nodes deleted")


# ======================================================================================
# Generation
# ======================================================================================
@command("NGEN", max_args=8)
def cmd_ngen(c):
    """NGEN,[itim],[step],[n1],[n2],[inc],[dx],[dy],[dz]: copy a node pattern itim times.

    New node ``n + k*step`` at ``x(n) + k*(dx, dy, dz)`` for k = 1..itim (active system).
    """
    m = c.model
    n1, n2 = _hist_defaults(m, 3, 4, c)
    if n2 < n1:
        n1, n2 = n2, n1
    itim = c.int(1, default=1)
    step = c.int(2, default=0) or (n2 - n1 + 1)
    inc = c.int(5, default=1)
    d = np.array([c.float(6, 0.0), c.float(7, 0.0), c.float(8, 0.0)])
    if itim < 1:
        c.warn("itim < 1: nothing generated")
        return
    s = _active(c)
    pattern = select_ids(m.nodes, n1, n2, inc)
    if not pattern:
        c.fail(f"no nodes defined in {n1}..{n2}")
    _check_new_ids(c, (n + k * step for k in (1, itim) for n in (pattern[0], pattern[-1])))
    base = {n: m.node_in_system(n, s) for n in pattern}
    over = 0
    for k in range(1, itim + 1):
        for n in pattern:
            new = n + k * step
            over += new in m.nodes
            m.define_node(new, base[n] + k * d, s)
    if over:
        c.warn(f"{over} existing nodes overwritten")
    c.confirm(f"{itim * len(pattern)} nodes generated from pattern {n1}..{n2}")


@command("FILL", max_args=3)
def cmd_fill(c):
    """FILL,[<n1>],[<n2>],[<nr>]: nr equally spaced nodes on the line n1-n2 (numbers n1 + k*d)."""
    m = c.model
    n1, n2 = _hist_defaults(m, 1, 2, c)
    for n in (n1, n2):
        if n not in m.nodes:
            c.fail(f"node {n} is not defined")
    if n2 < n1:
        n1, n2 = n2, n1
    nr = c.int(3, default=n2 - n1 - 1)
    if nr <= 0:
        c.info("FILL: no nodes to generate")
        return
    if (n2 - n1) % (nr + 1) == 0:
        d = (n2 - n1) // (nr + 1)
        ids = [n1 + k * d for k in range(1, nr + 1)]
    else:
        ids = [n1 + k for k in range(1, nr + 1)]
        c.warn(f"node increment (n2-n1)/(nr+1) is not an integer: nodes {ids[0]}..{ids[-1]} used")
    if ids[-1] >= n2:
        c.fail(f"generated numbers {ids[0]}..{ids[-1]} reach node {n2}")
    clash = [i for i in ids if i in m.nodes]
    if clash:
        c.fail(f"nodes {clash[:10]} already exist")
    s = _active(c)
    p1 = m.node_in_system(n1, s)
    p2 = m.node_in_system(n2, s)
    for k, nid in enumerate(ids, start=1):
        m.define_node(nid, p1 + (p2 - p1) * k / (nr + 1), s)
    c.confirm(f"{nr} nodes filled between {n1} and {n2}")


@command("NMED", max_args=9)
def cmd_nmed(c):
    """NMED,<nd>,<n1>,[<n2>,...,<n8>]: node at the mean of the global coordinates of 1 to 8 nodes.

    A 9th and later node are dropped with a warning (``max_args``) and are not used in the mean.
    """
    m = c.model
    nd = c.int(1, required=True, what="nd")
    if nd <= 0:
        c.fail("node numbers must be positive")
    src = [c.int(k) for k in range(2, c.nargs + 1) if c.given(k)]
    if not src:
        c.fail("at least one node is required")
    missing = [n for n in src if n not in m.nodes]
    if missing:
        c.fail(f"nodes {missing} are not defined")
    s = _active(c)
    avg = np.mean([m.node_global(n) for n in src], axis=0)
    existed = m.define_node(nd, m.to_system(avg, s), s)
    if existed:
        c.warn(f"node {nd} redefined")
    c.confirm(f"node {nd} at the mean of {len(src)} nodes")


def _copy_list(c, mode: str):
    """LMOVE (translation) and NMOVE (scaling) of an explicit node list (D-MDL-09)."""
    m = c.model
    if mode == "LMOVE":
        f = np.array([c.float(k, 0.0) for k in (1, 2, 3)])
    else:
        vals = []
        for k in (1, 2, 3):
            v = c.float(k, 1.0)
            if c.given(k) and v == 0.0:
                c.warn(f"factor {k} is 0: that coordinate collapses to 0 (spec 08 Q8)")
            vals.append(v)
        f = np.array(vals)
    nd = c.int(4, required=True, what="nd")
    if nd < 1:
        c.fail("node numbers must be positive (<nd> < 1)")
    # at most 15 list nodes (spec 08 Q6): later arguments are dropped by max_args = 19
    src = [c.int(k) for k in range(5, c.nargs + 1) if c.given(k)]
    if not src:
        c.fail("the node list must contain at least one node")
    missing = [n for n in src if n not in m.nodes]
    if missing:
        c.fail(f"source nodes {missing} are not defined")
    s = _active(c)
    base = [m.node_in_system(n, s) for n in src]
    over = 0
    for i, p in enumerate(base):
        new = nd + i
        over += new in m.nodes
        m.define_node(new, p + f if mode == "LMOVE" else p * f, s)
    if over:
        c.warn(f"{over} existing nodes overwritten")
    c.confirm(f"nodes {nd}..{nd + len(src) - 1} generated")


@command("LMOVE", max_args=19)
def cmd_lmove(c):
    """LMOVE,[<dx>],[<dy>],[<dz>],<nd>,<l1>,...,<l15>: node nd+i-1 = node l_i + (dx, dy, dz)."""
    _copy_list(c, "LMOVE")


@command("NMOVE", max_args=19)
def cmd_nmove(c):
    """NMOVE,[<dx>],[<dy>],[<dz>],<nd>,<l1>,...,<l15>: node nd+i-1 = (x dx, y dy, z dz) of node l_i
    (defaults 1)."""
    _copy_list(c, "NMOVE")


@command("NSCALE", max_args=6)
def cmd_nscale(c):
    """NSCALE,<n1>,<n2>,[<inc>],[<sfx>],[<sfy>],[<sfz>]: scale coordinates in place (factor 0 -> 1)."""
    m = c.model
    ids = range_args(c, 1, m.nodes)
    f = np.array([c.float(k, 0.0) or 1.0 for k in (4, 5, 6)])
    s = _active(c)
    for n in ids:
        _set_coords(m, n, m.node_in_system(n, s) * f, s)
    c.confirm(f"{len(ids)} nodes scaled by ({fmt_num(f[0])}, {fmt_num(f[1])}, {fmt_num(f[2])})")


# ======================================================================================
# Listing
# ======================================================================================
def _flag_text(n) -> str:
    return "".join(INT_LETTERS[k] for k in sorted(n.flags)) or "-"


@command("NLIST", max_args=3)
def cmd_nlist(c):
    """NLIST,[<n1>],[<n2>],[<inc>]: list nodes (stored coordinates, system, fixities, INT flags I/M/F/N)."""
    m = c.model
    ids = range_args(c, 1, m.nodes, all_default=True)
    c.info(f"{'node':>8} {'x':>14} {'y':>14} {'z':>14} {'csys':>5}  UX UY UZ RX RY RZ  flags")
    for i in ids:
        n = m.nodes[i]
        fx = "  ".join(str(v) for v in n.fix)
        c.info(f"{i:>8} {n.x:>14.6g} {n.y:>14.6g} {n.z:>14.6g} {n.csys:>5}  {fx}   {_flag_text(n)}")
    c.info(f"{len(ids)} nodes listed")


# ======================================================================================
# Fixities and node classification
# ======================================================================================
@command("D", max_args=10)
def cmd_d(c):
    """D,<n1>,<n2>,[<inc>],[<val>],<label1>,...,<label6>: set fixity codes (0 free, default; 1 fixed).

    At most six labels (spec 08 section 3.2); later ones are dropped with a warning (``max_args``).
    """
    m = c.model
    ids = range_args(c, 1, m.nodes)
    val = c.int(4, default=0)
    if val not in (0, 1):
        c.fail("<val> must be 0 (free) or 1 (fixed)")
    labels = [c.word(k) for k in range(5, c.nargs + 1) if c.given(k)]
    if not labels:
        c.fail("at least one DOF label is required (UX UY UZ ROTX ROTY ROTZ DISP ROT ALL)")
    dofs = set()
    for lab in labels:
        if lab not in DOF_LABELS:
            c.fail(f"unknown DOF label {lab}")
        dofs.update(DOF_LABELS[lab])
    for i in ids:
        fix = m.nodes[i].fix
        for k in dofs:
            fix[k] = val
    if not ids:
        c.warn("no defined node in the range")
    c.confirm(f"{len(ids)} nodes: {' '.join(DOF_NAMES[k] for k in sorted(dofs))} "
              f"{'fixed' if val else 'free'}")


@command("INT", max_args=5)
def cmd_int(c):
    """INT,<n1>,<n2>,[<inc>],<set>,[<code>]: set (1) / reset (0) node flag code 0 interaction (default),
    1 intermediate, 2 interface, 3 internal (D-MDL-07)."""
    m = c.model
    ids = range_args(c, 1, m.nodes)
    st = c.int(4, required=True, what="set")
    code = c.int(5, default=0)
    if st not in (0, 1):
        c.fail("<set> must be 1 (set) or 0 (reset)")
    if code not in INT_CODES:
        c.fail("<code> must be 0 (interaction), 1 (intermediate), 2 (interface) or 3 (internal)")
    for i in ids:
        if st:
            m.nodes[i].flags.add(code)
        else:
            m.nodes[i].flags.discard(code)
    if not ids:
        c.warn("no defined node in the range")
    c.confirm(f"{len(ids)} nodes: {INT_CODES[code]} flag {'set' if st else 'reset'}")


@command("INTLIST", max_args=7)
def cmd_intlist(c):
    """INTLIST,[<n1>],[<n2>],[<step>],[<c1>],[<c2>],[<c3>],[<c4>]: list classified nodes."""
    m = c.model
    ids = range_args(c, 1, m.nodes, all_default=True)
    want = {k for k in range(4) if c.int(4 + k, default=0) == 1}
    if not want:
        want = set(INT_CODES)
    rows = [i for i in ids if m.nodes[i].flags & want]
    c.info(f"{'node':>8}  flags  ({', '.join(INT_CODES[k] for k in sorted(want))})")
    for i in rows:
        c.info(f"{i:>8}  {_flag_text(m.nodes[i])}")
    c.info(f"{len(rows)} nodes listed")


# ======================================================================================
# Coordinate systems
# ======================================================================================
@command("CSYS", max_args=1)
def cmd_csys(c):
    """CSYS,<ns>: activate a coordinate system (0 = global)."""
    ns = c.int(1, required=True, what="ns")
    if ns != 0 and ns not in c.model.csys:
        c.fail(f"coordinate system {ns} is not defined")
    c.model.csys_active = ns
    c.confirm(f"coordinate system {ns} active")


def _check_type(c, k: int) -> int:
    t = c.int(k, default=0)
    if t == 1:
        c.fail("cylindrical systems (type 1) are not included in this version")
    if t != 0:
        c.fail("<type> must be 0 (Cartesian)")
    return t


def _store_system(c, cs: CoordSys) -> None:
    m = c.model
    if cs.id in m.csys:
        stored = m.nodes_in_system(cs.id)
        if stored:
            c.warn(f"system {cs.id} redefined: {len(stored)} nodes keep their local coordinates and move "
                   f"with the system (D-MDL-04)")
    m.csys[cs.id] = cs


@command("LOC", max_args=8)
def cmd_loc(c):
    """LOC,<ns>,<type>,<x0>,<y0>,<z0>,<txy>,<tyz>,<txz>: Cartesian system, R = Rz(txy) Rx(tyz) Ry(txz)."""
    ns = c.int(1, required=True, what="ns")
    if ns < 1:
        c.fail("system numbers are >= 1 (0 is the global system)")
    t = _check_type(c, 2)
    p = [c.float(k, 0.0) for k in range(3, 9)]
    cs = CoordSys(ns, origin=np.array(p[:3]), R=loc_matrix(p[3], p[4], p[5]), kind="LOC", type=t, params=p)
    _store_system(c, cs)
    c.confirm(f"system {ns} defined (LOC)")


@command("LOCAL", max_args=5)
def cmd_local(c):
    """LOCAL,<ns>,<type>,<n1>,<n2>,<n3>: system with origin n1, x toward n2, n3 in the +y half plane."""
    m = c.model
    ns = c.int(1, required=True, what="ns")
    if ns < 1:
        c.fail("system numbers are >= 1 (0 is the global system)")
    t = _check_type(c, 2)
    nn = [c.int(k, required=True, what=f"n{k - 2}") for k in (3, 4, 5)]
    missing = [n for n in nn if n not in m.nodes]
    if missing:
        c.fail(f"nodes {missing} are not defined")
    pts = [m.node_global(n) for n in nn]
    try:
        R = local_from_points(*pts)
    except ValueError as exc:
        c.fail(str(exc))
    cs = CoordSys(ns, origin=np.array(pts[0]), R=R, kind="LOCAL", type=t, nodes=tuple(nn),
                  points=[[float(v) for v in p] for p in pts])
    _store_system(c, cs)
    c.confirm(f"system {ns} defined from nodes {nn[0]}, {nn[1]}, {nn[2]}")


@command("GLOBAL", max_args=3)
def cmd_global(c):
    """GLOBAL,<n1>,<n2>,<inc>: convert stored coordinates to global (the active system is unchanged)."""
    m = c.model
    ids = range_args(c, 1, m.nodes, all_default=True)
    k = 0
    for i in ids:
        n = m.nodes[i]
        if n.csys != 0:
            _set_coords(m, i, m.node_global(i), 0)
            k += 1
    c.confirm(f"{k} nodes converted to global coordinates")


@command("SDEL", max_args=3)
def cmd_sdel(c):
    """SDEL,<s1>,[<s2>],[<inc>]: delete systems; their nodes are converted to global first."""
    m = c.model
    if c.int(1, default=None) == 0 and c.int(2, default=0) in (0, None):
        c.fail("the global system 0 cannot be deleted")
    ids = range_args(c, 1, m.csys)
    for s in ids:
        for i in m.nodes_in_system(s):
            _set_coords(m, i, m.node_global(i), 0)
        del m.csys[s]
        if m.csys_active == s:
            m.csys_active = 0
            c.warn(f"active system {s} deleted: system 0 is active")
    c.confirm(f"{len(ids)} coordinate systems deleted")


@command("SLIST", max_args=3)
def cmd_slist(c):
    """SLIST,[<s1>],[<s2>],[<inc>]: list coordinate systems."""
    m = c.model
    ids = range_args(c, 1, m.csys, all_default=True)
    for s in ids:
        cs = m.csys[s]
        o = ", ".join(f"{v:.6g}" for v in cs.origin)
        if cs.kind == "LOC":
            how = "angles txy, tyz, txz = " + ", ".join(f"{v:.6g}" for v in cs.params[3:6])
        else:
            how = "from nodes " + ", ".join(str(n) for n in cs.nodes)
        mark = " (active)" if s == m.csys_active else ""
        c.info(f"system {s}{mark}: {cs.kind}, Cartesian, origin ({o}), {how}")
        for lab, col in zip(("x", "y", "z"), np.asarray(cs.R).T):
            c.info(f"    local {lab} = ({col[0]:.9g}, {col[1]:.9g}, {col[2]:.9g})")
    c.info(f"{len(ids)} systems listed")
