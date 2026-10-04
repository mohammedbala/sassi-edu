"""Helpers of the HOUSE module: geometry tolerance, soil-layer interfaces, element classification,
excavation sets, restart hashes and mass totals.

Normative sources
-----------------
* requirements 4.0.4 and D-GEN-06 (geometric tolerance), 4.4 items 1-4 (classification of
  SOLID/PLANE elements, interaction-set rules, outputs), D-HOU-01 (ETYPE 0 "below ground" test),
  D-W1-13 (FILE90 hashes), ARCHITECTURE section 4 (depth = gelev - z; user interface i is the top
  of TOPL layer i, the last one the top of the half-space).

The functions here are pure (no file access) so that ANALYS and STRESS can reuse them, e.g. to
recompute an interaction-node interface or to compare restart hashes.

Background for the student
--------------------------
In the flexible-volume method the soil impedance is computed by POINT/ANALYS from point-load
solutions of the *layered* site, evaluated only at the user interfaces of the SITE profile.  An
interaction node therefore has to lie exactly on one of those horizontal planes (EDU-01), and the
excavated soil elements between two planes must carry the properties of the layer they replace
(D-ELM-12, EDU-08) -- otherwise the "excavated soil" subtracted by the method is not the soil that
the free field contains, and the zero-SSI identity (VP-16) is lost.
"""
from __future__ import annotations

import hashlib
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import numpy as np
import scipy.sparse as sp

__all__ = ["SOLID_FACES", "PLANE_EDGES", "geometric_tolerance", "interface_depths", "interface_tolerance",
           "match_interfaces", "layer_index_at_depth", "resolve_etype", "excavation_sets", "stable_hash",
           "sparse_items", "translational_mass", "general_matrix_tables", "first_items", "DOF_NAMES",
           "SOIL_TYPES", "SHELL_TYPES", "excavation_layering", "excstrchk_kind"]

#: faces of the 8-node hexahedron, 1-based local node numbers (spec 08 section 4.4)
SOLID_FACES = ((1, 2, 3, 4), (5, 6, 7, 8), (1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (4, 1, 5, 8))
#: edges of the 4-node PLANE element
PLANE_EDGES = ((1, 2), (2, 3), (3, 4), (4, 1))
#: load / DOF component names (FORCE listing, FILE4 listing)
DOF_NAMES = {1: "FX", 2: "FY", 3: "FZ", 4: "MX", 5: "MY", 6: "MZ"}


# ======================================================================================
# Tolerances and soil-layer interfaces
# ======================================================================================
def geometric_tolerance(xyz: np.ndarray, override: float = 0.0) -> float:
    """D-GEN-06: ``tol = max(1e-6 L_ref, 1e-9)`` with L_ref the bounding-box diagonal of ``xyz``.

    ``override`` > 0 replaces the computed value (``EDUOPT,GEOMTOL``)."""
    if override and override > 0:
        return float(override)
    X = np.asarray(xyz, dtype=float).reshape(-1, 3)
    X = X[np.all(np.isfinite(X), axis=1)]
    lref = float(np.linalg.norm(X.max(axis=0) - X.min(axis=0))) if len(X) else 0.0
    return max(1e-6 * lref, 1e-9)


def interface_depths(topl_thick: Sequence[float]) -> np.ndarray:
    """Depths of the user interfaces 1 .. nTOPL+1 (ARCHITECTURE section 4).

    Interface i is the top of TOPL layer i (interface 1 = ground surface, depth 0); the last
    interface ``nTOPL + 1`` is the top of the half-space (or the rigid base)."""
    h = np.asarray(list(topl_thick), dtype=float)
    return np.concatenate([[0.0], np.cumsum(h)])


def interface_tolerance(depths: np.ndarray, tol: float) -> float:
    """Interface match tolerance ``max(tol, 1e-4 h_min)`` (requirements 4.0.4, D-GEN-06)."""
    d = np.asarray(depths, dtype=float)
    hmin = float(np.min(np.diff(d))) if d.size > 1 else 0.0
    return max(float(tol), 1e-4 * max(hmin, 0.0))


def match_interfaces(z: np.ndarray, gelev: float, depths: np.ndarray, itol: float):
    """Nearest user interface of nodes at elevations ``z``.

    Returns ``(iface, misfit, above)``: the 1-based nearest interface, ``|depth - d_iface|`` and a
    flag for nodes above the ground surface (``z > gelev + itol``).  A node is *on* an interface
    when ``misfit <= itol`` (EDU-01 otherwise)."""
    z = np.asarray(z, dtype=float).ravel()
    d = float(gelev) - z
    dep = np.asarray(depths, dtype=float)
    if dep.size == 0:
        return np.zeros(z.size, int), np.full(z.size, np.inf), z > gelev + itol
    k = np.argmin(np.abs(d[:, None] - dep[None, :]), axis=1)
    mis = np.abs(dep[k] - d)
    return (k + 1).astype(np.int64), mis, z > float(gelev) + float(itol)


def layer_index_at_depth(depths: np.ndarray, depth: float) -> int:
    """1-based row of the SITE layer table (``sitelayers``) containing ``depth``: 1..nTOPL for the
    TOPL layers, nTOPL+1 for the half-space (a depth on an interface belongs to the layer below)."""
    return int(np.searchsorted(np.asarray(depths, float), float(depth), side="right"))


# ======================================================================================
# Classification (D-HOU-01)
# ======================================================================================
SOIL_TYPES = (1, 4)      # SOLID, PLANE
SHELL_TYPES = (3, 5)     # SHELL, TSHELL


def resolve_etype(code: int, etype: int, z_nodes: Sequence[float], gelev: float, tol: float) -> int:
    """Resolved ETYPE (requirements 4.4 item 1, D-HOU-01).

    * SOLID/PLANE: ETYPE 1 structure (M table), 2 excavated soil (L table via MSET).  ETYPE 0 is
      excavated when *every* node has ``z <= gelev + tol`` and the centroid of the distinct nodes is
      strictly below grade (``< gelev - tol``); otherwise structure.
    * SHELL/TSHELL: ETYPE 2 = buried shell (still structure, its nodes are interaction nodes);
      anything else structure.
    * BEAMS, SPRING, GENERAL: always structure (1); ETYPE has no effect.
    ``z_nodes`` are the elevations of the element's distinct nodes."""
    code, etype = int(code), int(etype)
    if code in SOIL_TYPES:
        if etype in (1, 2):
            return etype
        z = np.asarray(list(z_nodes), dtype=float)
        if z.size and np.all(z <= gelev + tol) and float(z.mean()) < gelev - tol:
            return 2
        return 1
    if code in SHELL_TYPES:
        return 2 if etype == 2 else 1
    return 1


# ======================================================================================
# Excavation sets (spec 09 section 2.19; EXCSTRCHK, EDU-22)
# ======================================================================================
def excavation_sets(elements: Iterable[Tuple[int, Sequence[int]]], z_of: Mapping[int, float], gelev: float,
                    tol: float) -> Dict[str, Set[int]]:
    """Node sets of the excavated soil volume.

    ``elements`` are the excavated SOLID/PLANE elements as ``(code, nodes)`` (nodes in element
    order, 0 = unused slot).  The boundary of the excavation consists of the faces (edges in 2D)
    used by exactly one excavated element; faces with all nodes at grade form the top surface,
    the others the lateral and bottom surface (the FI-FSIN set).

    Returns ``nodes`` (all excavation nodes), ``boundary``, ``fsin`` (lateral/bottom boundary),
    ``top`` (nodes of the ground-surface faces) and ``interior`` (nodes - boundary)."""
    faces: Dict[Tuple[int, ...], int] = {}
    nodes: Set[int] = set()
    for code, ns in elements:
        ns = [int(n) for n in ns]
        nodes.update(n for n in ns if n)
        if int(code) == 1:
            pad = (ns + [0] * 8)[:8]
            for f in SOLID_FACES:
                key = tuple(sorted({pad[k - 1] for k in f} - {0}))
                if len(key) >= 3:
                    faces[key] = faces.get(key, 0) + 1
        else:
            pad = (ns + [0] * 4)[:4]
            tri = pad[3] == 0 or pad[3] == pad[2]
            edges = ((1, 2), (2, 3), (3, 1)) if tri else PLANE_EDGES
            for e in edges:
                key = tuple(sorted({pad[k - 1] for k in e} - {0}))
                if len(key) == 2:
                    faces[key] = faces.get(key, 0) + 1
    boundary_faces = [k for k, c in faces.items() if c == 1]
    boundary: Set[int] = set()
    top: Set[int] = set()
    fsin: Set[int] = set()
    for k in boundary_faces:
        boundary.update(k)
        if all(abs(float(z_of[n]) - gelev) <= tol for n in k):
            top.update(k)
        else:
            fsin.update(k)
    return {"nodes": nodes, "boundary": boundary, "fsin": fsin, "top": top, "interior": nodes - boundary}


def excstrchk_kind(owners: Iterable[Tuple[int, Sequence[int]]], excavation_nodes: Set[int],
                   interaction: bool) -> str:
    """Severity ('Error' or 'Warning') of an excavation-*interior* node shared with the structure
    (EXCSTRCHK, G-15).  One rule for CHECK and HOUSE (final audit, docs/verification/AUDIT_REPORT.md).

    ``owners`` are ``(code, nodes)`` of the structural elements at the node -- resolved ETYPE != 2
    SOLID/PLANE/SHELL/TSHELL and every BEAMS/SPRING/GENERAL element; a beam/GENERAL orientation (K)
    node is passed with code 0 -- ``excavation_nodes`` the nodes of the excavated soil elements and
    ``interaction`` whether the node is an interaction node.

    The manual forbids connecting interior excavation nodes to the structure (1.5.1 rules 4 and 11,
    9.8.1 "incorrect from a SASSI modeling point of view", application guideline 13a "none of the
    internal excavation volume nodes shall be connected to a structural node"):

    * **Error** when the node is not an interaction node (FI methods): the structure is tied to the
      subtracted excavated soil ``-(K*_e - w^2 M_e)`` without the soil impedance -- "very poor results,
      especially if FI methods are used" (rule 11);
    * **Error** when any owner is a BEAMS, SHELL, TSHELL, SPRING or GENERAL element (or a K node), or a
      SOLID/PLANE element not built on the excavation mesh: in FV the node is tied to the soil continuum
      through the discretisation residual of ``X_ff - (K*_e - w^2 M_e)``, which can suppress the dynamics
      of flexible subsystems (rule 11: "can affect significantly the accuracy of their SSI responses");
    * **Warning** when the node is an interaction node and *every* owner is a SOLID/PLANE element whose
      nodes all belong to the excavation (near-field / backfill soil or a solid basement modelled on the
      excavation mesh, manual rule 4 and 2.4): the manual's "FV method provides often close results for
      SSI models with ... unique mesh in the basement"; exact when the elements reproduce the excavated
      soil (zero-SSI identity, VP-16 / VP-T3)."""
    own = [(int(c), [int(n) for n in ns if n]) for c, ns in owners]
    if not own or not interaction:
        return "Error"
    volumetric = all(c in SOIL_TYPES and ns and set(ns) <= excavation_nodes for c, ns in own)
    return "Warning" if volumetric else "Error"


# ======================================================================================
# Excavation numbering and grouping rules R5/R6 (spec 05b section 1.5, D-HOU-04)
# ======================================================================================
def excavation_layering(elems: Iterable[Tuple[int, int, int, int]]) -> Dict[str, List[str]]:
    """Check the manual's layering rules for the excavated soil elements (D-HOU-04).

    ``elems`` are ``(group, element id, TOPL row index, layer number)`` of every excavated
    element, where the row index (1 = top TOPL layer, nTOPL + 1 = half-space) is the layer that
    contains the element's mid-depth (:func:`layer_index_at_depth`) and the layer number is the
    L number of that row (used in the messages only).

    * **R5** -- with embedment the excavation elements are numbered from the surface down, so the
      layer row must not decrease with the element number inside an excavation group.
    * **R6** -- one excavation group per embedment layer (a horizontal slab matching one far-field
      layer: a group must not span several layers and a layer must not be split over several
      groups) and the groups numbered from the surface to the foundation.

    Why it matters: STRESS recovers the soil pressures and the nodal contours of the excavation
    group by group, layer by layer; the manual makes R6 a strict rule when those outputs are
    wanted.  For the impedance and the SSI response the numbering is irrelevant, so HOUSE only
    warns (D-HOU-04).  Returns ``{"R5": [messages], "R6": [messages]}`` (empty lists = rules met)."""
    by_group: Dict[int, List[Tuple[int, int, int]]] = {}
    for g, e, idx, lay in elems:
        by_group.setdefault(int(g), []).append((int(e), int(idx), int(lay)))
    out: Dict[str, List[str]] = {"R5": [], "R6": []}
    if not by_group:
        return out
    name: Dict[int, int] = {}                                   # row index -> layer number
    for rows in by_group.values():
        for _, idx, lay in rows:
            name.setdefault(idx, lay)
    # ---- R5: layer row non-decreasing with the element number inside each group --------------
    for g in sorted(by_group):
        rows = sorted(by_group[g])
        idx = np.array([r[1] for r in rows])
        bad = np.flatnonzero(np.diff(idx) < 0)
        if bad.size:
            k = int(bad[0])
            out["R5"].append(f"excavation group {g}: element {rows[k + 1][0]} (layer {rows[k + 1][2]}) follows "
                             f"element {rows[k][0]} (layer {rows[k][2]}) -- number the excavation elements "
                             f"from the surface down ({bad.size} inversion(s))")
    # ---- R6: one group per layer, groups numbered surface -> foundation ---------------------------
    spans = {g: sorted({r[1] for r in rows}) for g, rows in by_group.items()}
    for g in sorted(spans):
        if len(spans[g]) > 1:
            out["R6"].append(f"excavation group {g} spans {len(spans[g])} layers "
                             f"({', '.join(str(name[i]) for i in spans[g])}): use one group per embedment layer")
    owners: Dict[int, List[int]] = {}
    for g in sorted(spans):
        for i in spans[g]:
            owners.setdefault(i, []).append(g)
    for i in sorted(owners):
        if len(owners[i]) > 1:
            out["R6"].append(f"embedment layer {name[i]} is split over excavation groups "
                             f"{', '.join(str(g) for g in owners[i])}: use one group per embedment layer")
    order = sorted(spans)
    top = [spans[g][0] for g in order]                          # shallowest layer of each group
    for k in range(1, len(order)):
        if top[k] < top[k - 1]:
            out["R6"].append(f"excavation groups are not numbered from the surface down: group {order[k]} "
                             f"(layer {name[top[k]]}) follows group {order[k - 1]} (layer {name[top[k - 1]]})")
            break
    return out


# ======================================================================================
# Hashes (FILE90, D-W1-13) and totals
# ======================================================================================
def stable_hash(items: Iterable[Tuple[str, object]]) -> str:
    """SHA-256 of named arrays/scalars (name, dtype, shape and raw bytes, in the given order).

    Bit-identical inputs give identical hashes on one platform (D-GEN-09).  Used for the FILE90
    restart keys ``int_hash``, ``layer_hash`` and ``file4_hash``."""
    h = hashlib.sha256()
    for name, value in items:
        a = np.ascontiguousarray(np.asarray(value))
        if a.dtype.kind in "US":
            a = np.ascontiguousarray(a.astype("U"))
        h.update(str(name).encode())
        h.update(a.dtype.str.encode())
        h.update(repr(a.shape).encode())
        h.update(a.tobytes())
    return h.hexdigest()


def sparse_items(name: str, mat) -> List[Tuple[str, object]]:
    """(name, array) items of a sparse matrix in canonical CSR form, for :func:`stable_hash`."""
    m = sp.csr_matrix(mat)
    m.sort_indices()
    return [(name + "__data", m.data), (name + "__indices", m.indices.astype(np.int64)),
            (name + "__indptr", m.indptr.astype(np.int64)), (name + "__shape", np.asarray(m.shape, np.int64))]


def translational_mass(M, eq_dof: np.ndarray) -> np.ndarray:
    """Rigid-body mass ``r_d^T M r_d`` for the three global translations d = X, Y, Z.

    ``r_d`` has ones on the equations of DOF d (``eq_dof == d``).  Both the consistent and the
    lumped mass preserve the element mass, so for a set of elements whose DOFs are all active this
    is exactly the sum of rho * volume (VP-H2) plus the translational nodal masses."""
    M = sp.csr_matrix(M)
    eq_dof = np.asarray(eq_dof)
    out = np.zeros(3)
    for d in (1, 2, 3):
        r = (eq_dof == d).astype(float)
        # element-wise product + sum (not a BLAS dot: the macOS Accelerate matmul raises spurious
        # floating-point flags, see sassi.elements.base.blas_quiet)
        out[d - 1] = float(np.sum(r * (M @ r))) if r.any() else 0.0
    return out


def general_matrix_tables(rows: Iterable[Mapping[str, object]]) -> Dict[int, Dict[str, np.ndarray]]:
    """GENERAL matrix properties from the HOUSE deck ``matrices`` table.

    Each row ``(prop, kind, row, t1..t12)`` holds row ``row`` (1..12) of the upper triangle of the
    real stiffness (kind R), imaginary stiffness (I) or mass/weight (M) of property ``prop``:
    ``t1 = A[r, r] ... t(13-r) = A[r, 12]`` (spec 08 sections 5.26-5.29).  The result maps
    ``prop -> {'R', 'I', 'M'} -> (12, 12)`` tables in that row layout (row r in its first 13 - r
    columns), ready for :func:`sassi.elements.general.from_upper_rows`.  A row given twice keeps the
    last definition."""
    out: Dict[int, Dict[str, np.ndarray]] = {}
    for r in rows:
        prop = int(r["prop"])
        kind = str(r["kind"]).strip().upper()
        if kind not in ("R", "I", "M"):
            raise ValueError(f"matrix property {prop}: kind '{r['kind']}' must be R, I or M")
        row = int(r["row"])
        if not 1 <= row <= 12:
            raise ValueError(f"matrix property {prop}: row {row} outside 1..12")
        terms = np.array([float(r[f"t{k}"]) for k in range(1, 13)])
        if np.any(terms[13 - row:] != 0.0):
            raise ValueError(f"matrix property {prop} {kind} row {row}: at most {13 - row} terms (upper triangle)")
        tab = out.setdefault(prop, {k: np.zeros((12, 12)) for k in ("R", "I", "M")})[kind]
        tab[row - 1, :] = terms
    return out


def first_items(values: Sequence[int], n: int = 12) -> str:
    """'1 2 3 ... (+k more)' for listings."""
    v = [str(int(x)) for x in values]
    txt = " ".join(v[:n])
    return txt + (f" ... (+{len(v) - n} more)" if len(v) > n else "")
