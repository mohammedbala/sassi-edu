"""Node renumbering for the HOUSE "Optimize Model" option: reverse Cuthill-McKee with the interaction
nodes kept in their bottom-up order (requirements 4.4 item 5, D-HOU-03, D-FIL-09; spec 05b 3.7, OQ21 and
section 12 test 8).

Why renumber (for the structural engineer)
------------------------------------------
HOUSE numbers the equations node by node in the order of the node numbers.  A profile (skyline)
solver -- the solver of the original SASSI -- stores every row of the stiffness matrix from its first
non-zero to the diagonal, so its memory and time grow with the *profile*
``sum_i (i - first_i)`` and the *bandwidth* ``max |i - j|`` over coupled equations i, j.  Models
generated group by group (walls, slabs, the excavated soil...) often couple nodes whose numbers are far
apart; renumbering the nodes so that coupled nodes get close numbers can shrink the profile by orders of
magnitude.  SASSI-EDU's own sparse solver chooses its fill-reducing ordering internally, so the
optimizer does not change the results (VP-O1: identical transfer functions through the map, 1e-10);
it is kept for fidelity with the manual (the ``.hounew`` / ``.map`` files and the new numbering of
FILE4 and FILE8).

Algorithm
---------
1. Node graph: two nodes are adjacent when they carry DOFs of the same element (structure and
   excavated soil; orientation K nodes carry no DOF and are not connected).
2. The parts of the model that contain interaction nodes (normally all of it): one Cuthill-McKee
   breadth-first search (Cuthill & McKee 1969; George & Liu 1981) started from **all interaction nodes
   at once** -- level 0 is the interaction set in reverse interaction-table order -- visiting the
   neighbours by increasing degree; the visiting order is reversed (RCM).  The interaction nodes thus get
   the last numbers of their part, ascending in their bottom-up table order (D-HOU-03 by construction,
   no re-sorting), and the level structure grows from the interaction surface into the structure, so
   coupled nodes stay close even for wide, shallow meshes.  The other nodes are numbered from the far end
   of the structure towards the interaction surface (top-down above a foundation): the manual's bottom-up
   numbering (spec 05b R4) is a recommendation for hand numbering; only the interaction nodes must
   ascend bottom-up (EDU-21).
3. Components without interaction nodes: classical RCM from a pseudo-peripheral node (George-Liu
   iteration), started at the higher end (largest z) of the pseudo-diameter so the reversed order runs
   bottom-up (spec 05b R4); they are numbered first.
4. Nodes without element DOFs (gap, unused, orientation-only nodes) that are not interaction nodes are
   numbered last, in their original order.
5. Acceptance (spec 05b section 12 test 8, "bandwidth/profile does not increase"): the RCM order is
   used only if neither its bandwidth nor its profile exceeds that of the original order and it improves
   one of them (or makes non-ascending interaction numbers ascending); otherwise the original order is
   kept.  HOUSE measures its assembled equation matrix ``|K_s| + |K_e|`` (the equations of a node in
   ascending DOF order, nodes in the candidate order, :func:`equation_measures`); without it the node
   graph is measured.

Why the interaction nodes are the roots (and not re-sorted afterwards): ANALYS adds the free-field
impedance X_ff, which couples *every* interaction DOF with every other, to ``K_s - K_e``.  Numbered last and
contiguously, the interaction nodes confine that dense block to the last rows, where it adds little or
nothing to the profile (each interaction row already reaches back to its neighbours of the previous level,
numbered before the block).  A classical RCM swept
across a wide, shallow mesh scatters them through the order; permuting them back into table order then
couples rows far apart (bandwidth 159 -> 490 on a scrambled 12 x 12 x 4 hex grid) and the dense X_ff block
spans the whole profile (system profile 68722 with the rooted search, 74114 for the re-sorted RCM, although
the latter's ``K`` profile alone is smaller, 49821).

New node numbers are 1..N in the optimised order; :func:`write_map` writes the ``old new`` pairs
(one per line, D-FIL-09).
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

__all__ = ["node_adjacency", "pseudo_peripheral_pair", "cuthill_mckee", "reverse_cuthill_mckee", "bandwidth",
           "profile", "equation_measures", "keep_relative_order", "RenumberResult", "optimize_numbering",
           "write_map", "read_map"]

PathLike = Union[str, Path]


# ======================================================================================
# graph
# ======================================================================================
def node_adjacency(n: int, cliques: Iterable[Sequence[int]]) -> sp.csr_matrix:
    """Symmetric boolean adjacency (n, n) of nodes 0..n-1: every pair of nodes of one clique (the DOF
    nodes of one element) is connected; the diagonal is empty."""
    rows: List[np.ndarray] = []
    cols: List[np.ndarray] = []
    for c in cliques:
        idx = np.unique(np.asarray(list(c), dtype=np.int64))
        if idx.size < 2:
            continue
        r, s = np.meshgrid(idx, idx, indexing="ij")
        m = r != s
        rows.append(r[m])
        cols.append(s[m])
    if not rows:
        return sp.csr_matrix((n, n), dtype=bool)
    r = np.concatenate(rows)
    c = np.concatenate(cols)
    A = sp.csr_matrix((np.ones(r.size, dtype=bool), (r, c)), shape=(n, n))
    A.sum_duplicates()
    A.data[:] = True
    return A


def _neighbours(adj: sp.csr_matrix) -> List[np.ndarray]:
    return [adj.indices[adj.indptr[i]:adj.indptr[i + 1]] for i in range(adj.shape[0])]


def _levels(nbr: List[np.ndarray], start: int) -> List[List[int]]:
    """Breadth-first level structure rooted at ``start`` (its connected component only)."""
    seen = {start}
    levels = [[start]]
    while True:
        nxt: List[int] = []
        for v in levels[-1]:
            for u in nbr[v]:
                u = int(u)
                if u not in seen:
                    seen.add(u)
                    nxt.append(u)
        if not nxt:
            return levels
        levels.append(nxt)


def pseudo_peripheral_pair(adj: sp.csr_matrix, start: int, nbr: Optional[List[np.ndarray]] = None) -> Tuple[int, int]:
    """Two nodes at (nearly) maximal distance in the component of ``start`` (George & Liu 1979): from
    the last level of a rooted level structure take a node of minimum degree; repeat while the
    eccentricity grows.  Returns ``(root, far)``."""
    nbr = _neighbours(adj) if nbr is None else nbr
    deg = np.diff(adj.indptr)
    root = int(start)
    lv = _levels(nbr, root)
    while True:
        last = sorted(lv[-1], key=lambda v: (deg[v], v))
        far = last[0]
        lf = _levels(nbr, far)
        if len(lf) > len(lv):
            root, lv = far, lf
            continue
        return root, far


def cuthill_mckee(adj: sp.csr_matrix, start, nbr: Optional[List[np.ndarray]] = None) -> List[int]:
    """Cuthill-McKee order of the nodes reachable from ``start`` (one node, or a sequence of nodes that
    forms level 0 in the given order): breadth first, unvisited neighbours by increasing degree (ties by
    index)."""
    nbr = _neighbours(adj) if nbr is None else nbr
    deg = np.diff(adj.indptr)
    order = [int(start)] if np.ndim(start) == 0 else [int(v) for v in dict.fromkeys(int(v) for v in start)]
    seen = np.zeros(adj.shape[0], dtype=bool)
    seen[order] = True
    q = deque(order)
    while q:
        v = q.popleft()
        cand = [int(u) for u in nbr[v] if not seen[u]]
        cand.sort(key=lambda u: (deg[u], u))
        for u in cand:
            seen[u] = True
            order.append(u)
            q.append(u)
    return order


def reverse_cuthill_mckee(adj: sp.csr_matrix, z: Optional[Sequence[float]] = None,
                          isolated: Optional[Sequence[bool]] = None, roots: Sequence[int] = ()) -> np.ndarray:
    """RCM order (array of node indices) of all nodes.

    ``roots`` (the interaction nodes in their required order): the components that contain them are
    numbered by one CM search started from all of them at once (level 0 = ``roots`` reversed) and
    reversed, so the roots come last, in their given order.  The other components (in the order of their
    smallest node index) are numbered first, each by RCM from the higher end (largest ``z``, if given) of
    its pseudo-diameter.  Nodes flagged ``isolated`` (no element DOFs, not roots) are appended at the end
    in index order."""
    n = adj.shape[0]
    nbr = _neighbours(adj)
    deg = np.diff(adj.indptr)
    iso = np.zeros(n, dtype=bool) if isolated is None else np.asarray(isolated, dtype=bool)
    zz = np.zeros(n) if z is None else np.asarray(z, dtype=float)
    roots = [int(v) for v in dict.fromkeys(int(v) for v in roots)]
    done = np.zeros(n, dtype=bool)
    rooted: List[int] = []
    if roots:
        rooted = cuthill_mckee(adj, list(reversed(roots)), nbr)[::-1]
        done[rooted] = True
    out: List[int] = []
    tail: List[int] = []
    for i in range(n):
        if done[i]:
            continue
        if iso[i] and deg[i] == 0:
            done[i] = True
            tail.append(i)
            continue
        if deg[i] == 0:
            done[i] = True
            out.append(i)
            continue
        a, b = pseudo_peripheral_pair(adj, i, nbr)
        start = a if (zz[a], -a) >= (zz[b], -b) else b           # start CM at the top: RCM runs bottom-up
        comp = cuthill_mckee(adj, start, nbr)
        done[comp] = True
        out.extend(reversed(comp))
    return np.asarray(out + rooted + tail, dtype=np.int64)


# ======================================================================================
# measures
# ======================================================================================
def _positions(order: Sequence[int], n: int) -> np.ndarray:
    pos = np.empty(n, dtype=np.int64)
    pos[np.asarray(order, dtype=np.int64)] = np.arange(n)
    return pos


def bandwidth(adj: sp.csr_matrix, order: Sequence[int]) -> int:
    """max |p_i - p_j| over adjacent nodes, p = position in ``order``."""
    A = sp.coo_matrix(adj)
    if A.nnz == 0:
        return 0
    p = _positions(order, adj.shape[0])
    return int(np.max(np.abs(p[A.row] - p[A.col])))


def profile(adj: sp.csr_matrix, order: Sequence[int]) -> int:
    """Envelope size ``sum_i (p_i - min(p_i, min_{j adj i} p_j))`` of the node-level matrix in ``order``
    (the number of below-diagonal entries a profile solver stores, per DOF block)."""
    n = adj.shape[0]
    A = sp.coo_matrix(adj)
    p = _positions(order, n)
    first = p.copy()
    if A.nnz:
        np.minimum.at(first, A.row, p[A.col])
    return int(np.sum(p - first))


def equation_measures(order: Sequence[int], eq_node: np.ndarray, eq_dof: np.ndarray, pattern) -> Tuple[int, int]:
    """(bandwidth, profile) of an equation matrix after renumbering the nodes in ``order``.

    ``eq_node`` (neq,): the node *index* (0..n-1) of every equation of ``pattern`` (sparse (neq, neq), its
    stored non-zeros are the couplings), ``eq_dof`` its DOF label; the equations are renumbered node by
    node in ``order``, the DOFs of a node in ascending order (the HOUSE DOF map).  Bandwidth
    ``max |i - j|`` over the non-zeros; profile ``sum_i (i - first non-zero column of row i)``."""
    pos = _positions(order, len(order))
    key = pos[np.asarray(eq_node, dtype=np.int64)] * 8 + np.asarray(eq_dof, dtype=np.int64)
    new = np.empty(key.size, dtype=np.int64)
    new[np.argsort(key, kind="stable")] = np.arange(key.size)
    A = sp.coo_matrix(pattern)
    if A.nnz == 0:
        return 0, 0
    r, c = new[A.row], new[A.col]
    first = np.arange(key.size, dtype=np.int64)
    np.minimum.at(first, r, c)
    np.minimum.at(first, c, r)
    return int(np.max(np.abs(r - c))), int(np.sum(np.arange(key.size) - first))


def keep_relative_order(order: Sequence[int], sequence: Sequence[int]) -> np.ndarray:
    """``order`` with the nodes of ``sequence`` permuted among the positions they occupy so that they
    appear in the order of ``sequence`` (all other nodes keep their positions)."""
    order = np.array(order, dtype=np.int64, copy=True)
    seq = np.asarray(list(sequence), dtype=np.int64)
    if seq.size == 0:
        return order
    pos = _positions(order, order.size)
    slots = np.sort(pos[seq])
    order[slots] = seq
    return order


# ======================================================================================
# HOUSE optimizer
# ======================================================================================
@dataclass
class RenumberResult:
    """Outcome of :func:`optimize_numbering` (indices refer to the input node list)."""

    order: np.ndarray                   # node indices in the new order (new number k + 1 = order[k])
    new_number: np.ndarray              # (n,) new 1-based number of every input node
    bandwidth_before: int
    bandwidth_after: int
    profile_before: int
    profile_after: int
    method: str                         # 'RCM' or 'original order kept'
    n_components: int = 0
    n_isolated: int = 0
    notes: List[str] = field(default_factory=list)
    level: str = "node"                 # measures of the 'node' graph or of the 'equation' matrix
    candidates: List[str] = field(default_factory=list)   # measures of every candidate order

    @property
    def reduced(self) -> bool:
        return self.profile_after < self.profile_before


def _ascending(order: np.ndarray, keep: Sequence[int]) -> bool:
    if len(keep) < 2:
        return True
    pos = _positions(order, order.size)
    return bool(np.all(np.diff(pos[np.asarray(keep, dtype=np.int64)]) > 0))


def optimize_numbering(n: int, cliques: Iterable[Sequence[int]], keep: Sequence[int] = (),
                       z: Optional[Sequence[float]] = None, equations: Optional[Tuple[np.ndarray, np.ndarray]] = None,
                       pattern=None) -> RenumberResult:
    """New numbering of ``n`` nodes (indices 0..n-1 in their original order).

    ``cliques``: the DOF-node index lists of the elements; ``keep``: node indices whose relative order
    must be kept (the interaction nodes in interaction order; roots of the RCM search, they get ascending
    numbers); ``z``: elevations (components without interaction nodes run bottom-up); ``equations`` =
    ``(eq_node, eq_dof)`` and ``pattern``: the assembled equation matrix to measure (else the node graph,
    :func:`equation_measures`).  Nodes in no clique and not in ``keep`` are numbered last.  The RCM order is
    used if admissible (module docstring, step 5), else the original order -- spec 05b test 8."""
    cliques = [list(c) for c in cliques]
    adj = node_adjacency(n, cliques)
    deg = np.diff(adj.indptr)
    keep = [int(k) for k in keep]
    in_keep = np.zeros(n, dtype=bool)
    in_keep[keep] = True
    in_clique = np.zeros(n, dtype=bool)
    for c in cliques:
        in_clique[np.asarray(c, dtype=np.int64)] = True
    isolated = ~in_clique & ~in_keep
    original = np.arange(n, dtype=np.int64)
    candidates = {"RCM": reverse_cuthill_mckee(adj, z=z, isolated=isolated, roots=keep)}
    if equations is not None and pattern is not None:
        measure = lambda o: equation_measures(o, equations[0], equations[1], pattern)
        level = "equation"
    else:
        measure = lambda o: (bandwidth(adj, o), profile(adj, o))
        level = "node"
    bw0, pf0 = measure(original)
    asc0 = _ascending(original, keep)
    ncomp = connected_components(adj, directed=False)[0] - int(np.sum(deg == 0))
    best = None
    tried = []
    for name, o in candidates.items():
        bw1, pf1 = measure(o)
        tried.append(f"{name}: bandwidth {bw1}, profile {pf1}")
        asc1 = _ascending(o, keep)
        admissible = bw1 <= bw0 and pf1 <= pf0 and asc1 and (bw1 < bw0 or pf1 < pf0 or not asc0)
        if admissible and (best is None or (pf1, bw1) < (best[2], best[1])):
            best = (name, bw1, pf1, o)
    notes = []
    if best is not None:
        method, bwa, pfa, order = best
    else:
        method, bwa, pfa, order = "original order kept", bw0, pf0, original
        notes.append(f"reverse Cuthill-McKee does not reduce the {level} bandwidth/profile of the original order "
                     f"(bandwidth {bw0}, profile {pf0}) without increasing the other ({'; '.join(tried)}): the "
                     "original order is kept (spec 05b test 8)")
        if not asc0:
            notes.append("the interaction nodes keep their original, not ascending, numbers (EDU-21)")
    new = np.empty(n, dtype=np.int64)
    new[order] = np.arange(1, n + 1)
    return RenumberResult(order=np.asarray(order, dtype=np.int64), new_number=new, bandwidth_before=bw0,
                          bandwidth_after=bwa, profile_before=pf0, profile_after=pfa, method=method,
                          n_components=int(ncomp), n_isolated=int(np.sum(isolated)), notes=notes, level=level,
                          candidates=tried)


# ======================================================================================
# .map files (D-FIL-09)
# ======================================================================================
def write_map(path: PathLike, old: Sequence[int], new: Sequence[int]) -> Path:
    """``<model>.map``: one ``old new`` pair per line, ascending old number (D-FIL-09)."""
    pairs = sorted(zip((int(o) for o in old), (int(v) for v in new)))
    Path(path).write_text("".join(f"{o} {v}\n" for o, v in pairs), encoding="utf-8")
    return Path(path)


def read_map(path: PathLike) -> Dict[int, int]:
    """``{old: new}`` from a ``.map`` file (blank lines and ``#`` comments ignored)."""
    out: Dict[int, int] = {}
    for k, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        s = line.split("#", 1)[0].strip()
        if not s:
            continue
        parts = s.replace(",", " ").split()
        if len(parts) != 2:
            raise ValueError(f"{Path(path).name} line {k}: expected 'old new', got {line!r}")
        old, new = int(parts[0]), int(parts[1])
        if old in out:
            raise ValueError(f"{Path(path).name} line {k}: node {old} mapped twice")
        out[old] = new
    return out
