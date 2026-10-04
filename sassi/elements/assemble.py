"""DOF management and sparse assembly of the HOUSE model matrices.

Normative sources: ARCHITECTURE 6.2 (binding API), requirements 4.1 "DOF management" and 4.4
item 2 (frequency-independent complex sparse ``K*_s, M_s`` of the structure incl. near field and
``K*_e, M_e`` of the excavated soil on one global DOF map), D-ELM-09 (DOFs no attached element
defines are removed; defined but unstiffened DOFs are kept), D-ELM-02/03 (no incompatible modes
in excavated soil), D-MDL-08 and requirements 4.0.3 (nodal masses, divided by g when MUNITS = 1;
ignored on fixed DOFs, W5/W6).

DOF numbering
-------------
Equations are numbered node-major in the order of ``node_ids`` and, within a node, by DOF label
1 UX, 2 UY, 3 UZ, 4 ROTX, 5 ROTY, 6 ROTZ.  The active DOFs of a node are the union of the DOFs of
the elements attached to it (and of ``extra_dofs``), minus the D-fixed DOFs.  Orientation nodes
(BEAMS / GENERAL K) carry no DOFs.

Recovery operators
------------------
For every output-capable element type T the assembly returns
``recovery[T] = dict(S=(nE, nc, nd) complex, eq=(nE, nd) int, group=(nE,), id=(nE,), excavated=(nE,))``
so that a component transfer function is ``S[e] @ u[eq[e]]`` (``eq = -1``: fixed DOF, contributes 0).
SOLID and PLANE also carry ``B`` (strain operator).  SHELL and TSHELL triangles are padded to the
4-node size with ``eq = -1`` and zero columns.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp

from ..conventions import ELEMENT_TYPE_NAMES
from . import plane as _plane
from . import shell as _shell
from .base import ELEMENTS, ElementError, MaterialProps, blas_quiet, nodal_mass_to_mass_units

__all__ = ["ElemRecord", "DofMap", "AssembledModel", "build_dofmap", "assemble", "element_dof_nodes",
           "is_excavated_soil", "SOIL_CODES", "nodal_masses_from_table", "unstiffened_rotations",
           "natural_frequencies"]

_BATCH_KEYS = {1: ("eint", "incompatible", "mass"), 4: ("incompatible", "mass")}

#: element codes that can represent excavated soil (requirements 4.4 item 1: SOLID 3D, PLANE 2D)
SOIL_CODES = (1, 4)

#: at most this many element ids are listed in an assembly message
_MAX_LISTED = 10


@dataclass
class ElemRecord:
    """One element of the analysis model (ARCHITECTURE 6.2).

    ``nodes`` as in the E command (unused trailing slots may be 0): SOLID 8 (repeated for prisms
    and pyramids), PLANE 4 (L = 0 or L = K for a triangle), BEAMS (I, J, K), SHELL and TSHELL
    (I, J, K, L) (L = 0, or one pair of adjacent repeated nodes such as L = K, for a triangle; other
    repeated nodes are Error 9), SPRING (I, J), GENERAL (I, J) global or (I, J, K) local.
    ``props`` are the element keyword properties of :mod:`sassi.elements`:

    * SOLID: ``eint`` (0/1/2), ``incompatible`` (True = include the modes, i.e. MOPT <incomp> = 0;
      forced off for excavated elements), ``mass`` ('mixed' default);
    * PLANE: ``incompatible``, ``mass``;
    * BEAMS: ``section`` = dict(A, As2, As3, J, I2, I3), ``ki``, ``kj`` (6-digit release codes);
    * SHELL: ``thick``, ``membrane_output`` ('stress' default or 'force');
    * TSHELL: ``thick``, ``eint`` (0 reduced, default; 1 selective; D-ELM-08);
    * SPRING: ``k`` = (kx, ky, kz, kxx, kyy, kzz), ``damp``, ``form`` (CMODFORM);
    * GENERAL: ``KR``, ``KI``, ``MM``, ``upper_rows``, ``munits`` (MOPT <matrix>), ``gravity``.

    ``excavated`` = the element is **excavated soil**, assembled into ``Ke, Me`` (the matrices the
    flexible-volume method subtracts).  Only SOLID and PLANE elements can be excavated soil
    (requirements 4.4 item 1, :data:`SOIL_CODES`).  The deck ETYPE value 2 means "excavated soil"
    for SOLID/PLANE but "buried shell" (structure whose nodes are interaction nodes) for SHELL/TSHELL,
    and ETYPE has no effect on BEAMS (spec 08 sections 4.6 and 5.11).  :func:`assemble` therefore
    treats ``excavated=True`` on any other element type as structure (``Ks, Ms``, recovery flag 0)
    and reports it in ``AssembledModel.messages``; callers may map ETYPE 2 to ``excavated`` for
    every type, but should prefer ``excavated = (etype == 2 and code in SOIL_CODES)``.
    """

    group: int
    id: int
    code: int
    nodes: tuple
    excavated: bool = False
    mat: Optional[MaterialProps] = None
    props: dict = field(default_factory=dict)


def is_excavated_soil(rec: ElemRecord) -> bool:
    """True when ``rec`` is assembled as excavated soil: ``excavated`` set on a SOLID or PLANE.

    An ETYPE 2 SHELL is a buried shell, i.e. structure (spec 08 4.6), and ETYPE has no effect on
    BEAMS (spec 08 5.11); SPRING and GENERAL elements are never soil (requirements 4.4 item 1).
    """
    return bool(rec.excavated) and int(rec.code) in SOIL_CODES


def _shell_corners(n: List[int], where: str) -> List[int]:
    """Distinct corner nodes of a SHELL record (CHECK Error 9 otherwise).

    A 4-node record with exactly one pair of *cyclically adjacent* repeated nodes, e.g.
    (I, J, K, K) (the usual "L = K" triangle), (I, J, K, I) or (I, J, J, L), is the triangle of
    its three distinct nodes, kept in the order of their first appearance (so the circulation, i.e.
    the normal z', is unchanged; the triangle local axes follow spec 08 4.6 with these three
    nodes as I, J, K).  Any other repetition, e.g. (I, J, I, L), folds the facet onto itself and
    is Error 9.  The same rule is applied to coordinates by :func:`sassi.elements.shell.corner_rows`.
    """
    if len(n) == 3:
        if len(set(n)) != 3:
            raise ElementError(f"Error 9: {where} has coincident nodes {tuple(n)}")
        return n
    if len(set(n)) == 4:
        return n
    n_adjacent = sum(n[a] == n[(a + 1) % 4] for a in range(4))
    if len(set(n)) == 3 and n_adjacent == 1:
        return list(OrderedDict.fromkeys(n))
    raise ElementError(f"Error 9: {where} has coincident nodes {tuple(n)} (only one pair of adjacent "
                       "corners may coincide, giving a triangle)")


def _strip(nodes: Sequence[int]) -> List[int]:
    n = [int(v) for v in nodes]
    while n and n[-1] == 0:
        n.pop()
    return n


def element_dof_nodes(rec: ElemRecord) -> Tuple[List[int], List[int]]:
    """(dof_nodes, geometry_nodes) of an element record.

    dof_nodes carry the element DOFs (in element order, repeated nodes kept for degenerate
    SOLID/PLANE); geometry_nodes are passed (as coordinates) to the element functions.
    """
    code = int(rec.code)
    n = _strip(rec.nodes)
    where = f"{ELEMENT_TYPE_NAMES.get(code, code)} element {rec.id} of group {rec.group}"
    if code not in ELEMENTS:
        raise ElementError(f"unknown element type {code} ({where})")
    if any(v <= 0 for v in n):
        raise ElementError(f"{where}: node numbers must be > 0 (got {tuple(rec.nodes)})")
    if code == 1:
        if len(n) != 8:
            raise ElementError(f"Error 7: {where} needs 8 nodes")
        return n, n
    if code == 4:
        if len(n) == 3:
            n = n + [n[2]]
        if len(n) != 4:
            raise ElementError(f"Error 7: {where} needs 4 nodes (3 for a triangle)")
        return n, n
    if code == 2:
        if len(n) != 3:
            raise ElementError(f"Error 7: {where} needs nodes I, J and K")
        return n[:2], n
    if code in (3, 5):                                   # SHELL, TSHELL: triangles by a repeated corner
        if len(n) not in (3, 4):
            raise ElementError(f"Error 7: {where} needs 3 or 4 nodes")
        n = _shell_corners(n, where)
        return n, n
    if code == 7:
        if len(n) != 2:
            raise ElementError(f"Error 7: {where} needs nodes I and J")
        return n, n
    if code == 9:
        if len(n) not in (2, 3):
            raise ElementError(f"Error 7: {where} needs nodes I, J (global) or I, J, K (local)")
        return n[:2], n
    raise ElementError(f"unsupported element type {code}")


# ======================================================================================
# DOF map
# ======================================================================================
@dataclass
class DofMap:
    """Global equation numbering (ARCHITECTURE 6.2)."""

    node_ids: np.ndarray          # (nN,) node ids in numbering order
    eq_table: np.ndarray          # (nN, 6) equation number or -1
    defined: np.ndarray           # (nN, 6) DOF defined by an element / extra_dofs
    fixed: np.ndarray             # (nN, 6) D-fixed
    eq_node: np.ndarray           # (nEq,)
    eq_dof: np.ndarray            # (nEq,) 1..6
    _index: Dict[int, int] = field(default_factory=dict, repr=False)

    @property
    def neq(self) -> int:
        return int(self.eq_node.size)

    def node_index(self, node: int) -> int:
        try:
            return self._index[int(node)]
        except KeyError:
            raise ElementError(f"node {node} is not defined in the analysis model (Error 41)") from None

    def eq(self, node: int, dof: int) -> int:
        """Equation number of (node, dof), or -1 if the DOF is fixed, undefined or the node unknown."""
        i = self._index.get(int(node))
        if i is None or not 1 <= int(dof) <= 6:
            return -1
        return int(self.eq_table[i, int(dof) - 1])

    def eqs(self, nodes: Sequence[int], dofs: Sequence[int]) -> np.ndarray:
        """Node-major equation numbers of an element: [n1 dofs..., n2 dofs...] (-1 = absent)."""
        idx = np.array([self.node_index(n) for n in nodes], dtype=int)
        d = np.asarray(dofs, dtype=int) - 1
        return self.eq_table[idx[:, None], d[None, :]].ravel()

    @property
    def n_eliminated(self) -> int:
        """DOFs removed because no element defines them (D-ELM-09 information message)."""
        return int((~self.defined).sum())

    def summary(self) -> str:
        nN = self.node_ids.size
        return (f"{nN} nodes, {self.neq} equations; {int((self.defined & self.fixed).sum())} DOFs fixed by D, "
                f"{self.n_eliminated} DOFs eliminated (not defined by any attached element)")


def build_dofmap(node_ids: Sequence[int], node_fix, records: Iterable[ElemRecord],
                 extra_dofs: Optional[Mapping[int, Iterable[int]]] = None) -> DofMap:
    """Active DOFs = union of the DOFs of the attached elements (+ ``extra_dofs``) minus fixed.

    ``node_fix`` is (nN, 6) with 1 = fixed (D command) in the order of ``node_ids``.
    """
    ids = np.asarray(list(node_ids), dtype=np.int64)
    nN = ids.size
    index = {int(n): i for i, n in enumerate(ids)}
    if len(index) != nN:
        raise ElementError("duplicate node ids in the analysis model")
    fix = np.zeros((nN, 6), dtype=bool) if node_fix is None else np.asarray(node_fix).astype(bool)
    if fix.shape != (nN, 6):
        raise ElementError("node_fix must have shape (number of nodes, 6)")
    defined = np.zeros((nN, 6), dtype=bool)
    for rec in records:
        dn, gn = element_dof_nodes(rec)
        dofs = np.asarray(ELEMENTS[int(rec.code)].dofs, dtype=int) - 1
        for n in gn:
            if int(n) not in index:
                raise ElementError(f"element {rec.id} of group {rec.group} references undefined node {n} "
                                   "(Error 41)")
        for n in dn:
            defined[index[int(n)], dofs] = True
    if extra_dofs:
        for n, ds in extra_dofs.items():
            if int(n) not in index:
                raise ElementError(f"extra DOFs given for undefined node {n}")
            for d in ds:
                defined[index[int(n)], int(d) - 1] = True
    active = defined & ~fix
    eq_table = np.full((nN, 6), -1, dtype=np.int64)
    flat = active.ravel()
    eq_table.ravel()[flat] = np.arange(int(flat.sum()))
    rows, cols = np.nonzero(active)
    return DofMap(node_ids=ids, eq_table=eq_table, defined=defined, fixed=fix,
                  eq_node=ids[rows], eq_dof=(cols + 1).astype(np.int64), _index=index)


# ======================================================================================
# Assembly
# ======================================================================================
@dataclass
class AssembledModel:
    """Model matrices on the DofMap numbering (ARCHITECTURE 6.2)."""

    Ks: sp.csr_matrix                 # complex stiffness: structure incl. near-field soil
    Ms: sp.csr_matrix                 # real mass of the structure (+ nodal masses)
    Ke: sp.csr_matrix                 # complex stiffness of the excavated soil
    Me: sp.csr_matrix                 # real mass of the excavated soil
    recovery: Dict[str, Dict[str, np.ndarray]]
    ignored_masses: List[Tuple[int, int, float]] = field(default_factory=list)   # (node, dof, mass)
    messages: List[str] = field(default_factory=list)


class _Triplets:
    def __init__(self, dtype):
        self.rows: List[np.ndarray] = []
        self.cols: List[np.ndarray] = []
        self.vals: List[np.ndarray] = []
        self.dtype = dtype

    def add(self, eq: np.ndarray, A: np.ndarray) -> None:
        """Scatter element matrices A (nE, nd, nd) with equations eq (nE, nd)."""
        eq = np.atleast_2d(eq)
        A = A.reshape((eq.shape[0], eq.shape[1], eq.shape[1]))
        r = np.broadcast_to(eq[:, :, None], A.shape)
        c = np.broadcast_to(eq[:, None, :], A.shape)
        m = (r >= 0) & (c >= 0) & (A != 0)
        self.rows.append(r[m])
        self.cols.append(c[m])
        self.vals.append(A[m].astype(self.dtype))

    def matrix(self, n: int) -> sp.csr_matrix:
        if not self.rows:
            return sp.csr_matrix((n, n), dtype=self.dtype)
        M = sp.coo_matrix((np.concatenate(self.vals), (np.concatenate(self.rows), np.concatenate(self.cols))),
                          shape=(n, n), dtype=self.dtype).tocsr()
        M.sum_duplicates()
        M.eliminate_zeros()
        return M


def _list_elements(recs: Sequence[ElemRecord]) -> str:
    """'group g element e, ...' for the first few records of a message."""
    txt = ", ".join(f"group {r.group} element {r.id}" for r in recs[:_MAX_LISTED])
    return txt + (f" and {len(recs) - _MAX_LISTED} more" if len(recs) > _MAX_LISTED else "")


def _xyz_of(node_xyz, nodes: Sequence[int]) -> np.ndarray:
    try:
        return np.array([np.asarray(node_xyz[int(n)], dtype=float)[:3] for n in nodes])
    except KeyError as exc:
        raise ElementError(f"no coordinates for node {exc}") from None


def _undamped_record(rec: ElemRecord) -> ElemRecord:
    """Copy of a record with zero damping: real moduli, SPRING damp = 0, GENERAL K_I = 0 (the
    GENERAL K_R is kept as entered)."""
    props = dict(rec.props or {})
    if int(rec.code) == 7:
        props["damp"] = 0.0
    if int(rec.code) == 9:
        props["KI"] = None
    mat = rec.mat.undamped() if rec.mat is not None else None
    return ElemRecord(rec.group, rec.id, rec.code, rec.nodes, rec.excavated, mat, props)


def assemble(node_xyz: Mapping[int, Sequence[float]], records: Sequence[ElemRecord], dofmap: DofMap,
             masses: Optional[Mapping[int, Sequence[float]]] = None, undamped: bool = False) -> AssembledModel:
    """Assemble ``Ks, Ms`` (structure) and ``Ke, Me`` (excavated soil) and the recovery operators.

    ``masses``: nodal masses in **mass units**, ``{node: (mx, my, mz, mxx, myy, mzz)}`` (use
    :func:`nodal_masses_from_table` to convert MT/MR rows with their MUNITS flag).  Masses on
    fixed or undefined DOFs are ignored and reported in ``ignored_masses`` (W5/W6).

    ``undamped=True`` (SASSI-EDU addition) assembles the undamped stiffness K0 (all damping ratios
    set to 0, GENERAL K_I dropped), e.g. for fixed-base modal checks.  Note that Re(K*) is *not* K0
    for the default SASSI complex modulus: Re c(beta) = 1 - 2 beta^2.
    """
    if undamped:
        records = [_undamped_record(r) for r in records]
    neq = dofmap.neq
    Ks, Ke = _Triplets(complex), _Triplets(complex)
    Ms, Me = _Triplets(float), _Triplets(float)
    rec_data: Dict[str, Dict[str, list]] = OrderedDict()
    messages: List[str] = []
    n_incomp_off = 0
    not_soil: Dict[str, List[ElemRecord]] = OrderedDict()   # excavated=True on a non-soil type
    plane_flipped: List[ElemRecord] = []                    # PLANE records given clockwise

    # ---- group the records: batchable SOLID/PLANE by kernel options, others one by one --------
    groups: Dict[tuple, List[ElemRecord]] = OrderedDict()
    position = {id(rec): k for k, rec in enumerate(records)}      # recovery rows follow input order
    for rec in records:
        code = int(rec.code)
        element_dof_nodes(rec)                                   # validates node lists early
        props = dict(rec.props or {})
        soil = is_excavated_soil(rec)
        if rec.excavated and not soil:
            not_soil.setdefault(ELEMENTS[code].name, []).append(rec)
        if code in SOIL_CODES:
            if soil and props.get("incompatible"):
                n_incomp_off += 1
            if soil:
                props["incompatible"] = False                    # D-ELM-02 / D-ELM-03
            key = (code, soil) + tuple(props.get(k, None) for k in _BATCH_KEYS[code])
        else:
            key = ("single", id(rec))
        groups.setdefault(key, []).append(rec)
    if n_incomp_off:
        messages.append(f"incompatible modes suppressed on {n_incomp_off} excavated-soil elements (D-ELM-02)")
    for name, recs in not_soil.items():
        what = (f"an ETYPE 2 {name} is a buried shell, i.e. structure whose nodes are interaction nodes "
                "(spec 08 4.6)" if name in ("SHELL", "TSHELL") else
                "ETYPE classifies only SOLID/PLANE (soil) and SHELL (buried) elements (spec 08 5.11)")
        messages.append(f"{len(recs)} {name} elements flagged excavated assembled as structure (Ks, Ms): "
                        f"only SOLID/PLANE elements can be excavated soil (requirements 4.4 item 1); {what}: "
                        + _list_elements(recs))

    def _collect(name: str, recs: List[ElemRecord], S: np.ndarray, eq: np.ndarray, B=None):
        d = rec_data.setdefault(name, {"S": [], "eq": [], "group": [], "id": [], "excavated": [], "B": [],
                                       "pos": []})
        d["pos"].extend(position[id(r)] for r in recs)
        d["S"].append(S)
        d["eq"].append(eq)
        d["group"].extend(int(r.group) for r in recs)
        d["id"].extend(int(r.id) for r in recs)
        d["excavated"].extend(int(is_excavated_soil(r)) for r in recs)
        if B is not None:
            d["B"].append(B)

    for key, recs in groups.items():
        code = int(recs[0].code)
        spec = ELEMENTS[code]
        excav = is_excavated_soil(recs[0])          # same for every record of a SOLID/PLANE batch
        if code in SOIL_CODES:
            for r in recs:
                if r.mat is None:
                    raise ElementError(f"{spec.name} element {r.id} of group {r.group} has no material")
            nodes = [element_dof_nodes(r)[0] for r in recs]
            xyz = np.stack([_xyz_of(node_xyz, n) for n in nodes])
            if code == 4:                            # spec 08 4.7 / Q12: reorder clockwise input, warn
                plane_flipped.extend(r for r, det in zip(recs, _plane.centroid_jacobian(xyz)) if det < 0)
            props = dict(recs[0].props or {})
            if excav:
                props["incompatible"] = False
            out = spec.batch(xyz, [r.mat for r in recs], **props)
            eq = np.stack([dofmap.eqs(n, spec.dofs) for n in nodes])
            (Ke if excav else Ks).add(eq, out["K"])
            (Me if excav else Ms).add(eq, out["M"])
            _collect(spec.name, recs, out["S"], eq, out.get("B"))
            continue
        rec = recs[0]
        dn, gn = element_dof_nodes(rec)
        xyz = _xyz_of(node_xyz, gn)
        if code in (3, 5) and len(_shell.corner_rows(xyz)) != len(gn):
            raise ElementError(f"Error 9: {spec.name} element {rec.id} of group {rec.group}: distinct nodes "
                               f"{tuple(gn)} have coincident coordinates (merge the nodes or give a "
                               "triangle with a repeated node number)")
        props = dict(rec.props or {})
        if code in (2, 3, 5) and rec.mat is None:
            raise ElementError(f"{spec.name} element {rec.id} of group {rec.group} has no material")
        K, M = spec.matrices(xyz, rec.mat, **props)
        eq = dofmap.eqs(dn, spec.dofs)
        (Ke if excav else Ks).add(eq[None], K[None])
        (Me if excav else Ms).add(eq[None], M[None])
        if spec.components:
            S = spec.recovery(xyz, rec.mat, **props)
            nd_full = spec.nnodes * len(spec.dofs)
            if S.shape[1] < nd_full:                              # SHELL/TSHELL triangle -> 4-node layout
                S = np.pad(S, ((0, 0), (0, nd_full - S.shape[1])))
                eq = np.concatenate([eq, np.full(nd_full - eq.size, -1, dtype=eq.dtype)])
            _collect(spec.name, [rec], S[None], eq[None])

    if plane_flipped:
        messages.append(f"{len(plane_flipped)} PLANE elements are numbered clockwise in the X-Z view (x right, "
                        "z up); their node order was reversed internally (spec 08 4.7, Q12): "
                        + _list_elements(plane_flipped))

    # ---- nodal masses (MT/MR) -------------------------------------------------------------
    ignored: List[Tuple[int, int, float]] = []
    if masses:
        rows, vals = [], []
        for node, m6 in masses.items():
            m6 = np.asarray(m6, dtype=float).ravel()
            if m6.size != 6:
                raise ElementError(f"nodal mass of node {node} needs 6 values (mx my mz mxx myy mzz)")
            for d in range(6):
                if m6[d] == 0.0:
                    continue
                e = dofmap.eq(node, d + 1)
                if e < 0:
                    ignored.append((int(node), d + 1, float(m6[d])))
                else:
                    rows.append(e)
                    vals.append(m6[d])
        if rows:
            r = np.asarray(rows, dtype=np.int64)
            Ms.rows.append(r)
            Ms.cols.append(r)
            Ms.vals.append(np.asarray(vals, dtype=float))
        if ignored:
            messages.append(f"{len(ignored)} nodal mass terms on fixed or undefined DOFs ignored (W5/W6)")

    recovery: Dict[str, Dict[str, np.ndarray]] = OrderedDict()
    for name, d in rec_data.items():
        order = np.argsort(np.asarray(d["pos"]), kind="stable")
        out = {"S": np.concatenate(d["S"], axis=0).astype(complex)[order],
               "eq": np.concatenate(d["eq"], axis=0).astype(np.int64)[order],
               "group": np.asarray(d["group"], dtype=np.int64)[order],
               "id": np.asarray(d["id"], dtype=np.int64)[order],
               "excavated": np.asarray(d["excavated"], dtype=np.int64)[order]}
        if d["B"]:
            out["B"] = np.concatenate(d["B"], axis=0).astype(complex)[order]
        recovery[name] = out
    return AssembledModel(Ks=Ks.matrix(neq), Ms=Ms.matrix(neq), Ke=Ke.matrix(neq), Me=Me.matrix(neq),
                          recovery=recovery, ignored_masses=ignored, messages=messages)


# ======================================================================================
# Utilities
# ======================================================================================
def nodal_masses_from_table(rows: Iterable[Sequence[float]], gravity: float) -> Dict[int, np.ndarray]:
    """Convert MT/MR rows ``(node, mx, my, mz, mxx, myy, mzz, units)`` to mass units.

    ``units`` = MUNITS flag: 0 mass, 1 weight (divided by ``gravity``), D-MDL-08; repeated rows of
    one node are added.  UT-15.
    """
    out: Dict[int, np.ndarray] = {}
    for row in rows:
        row = list(row)
        node = int(row[0])
        units = int(row[7]) if len(row) > 7 else 1
        m = nodal_mass_to_mass_units(row[1:7], units, gravity)
        out[node] = out.get(node, np.zeros(6)) + m
    return out


def unstiffened_rotations(K, dofmap: DofMap, rtol: float = 1e-8) -> List[Tuple[int, np.ndarray]]:
    """Nodes whose active rotational DOFs have (nearly) no stiffness in some direction.

    For every node with active rotations, the eigenvalues of the rotational block of Re(K) are
    compared with the node's largest rotational diagonal over the whole model; an eigenvalue below
    ``rtol`` times that scale flags a rotation axis without stiffness (e.g. the drilling rotation
    of a Kirchhoff SHELL node, EDU-06).  Returns [(node, unit axis in global XYZ)].
    """
    K = sp.csr_matrix(K)
    diag = np.abs(K.diagonal().real)
    rot = dofmap.eq_table[:, 3:6]
    has = (rot >= 0).any(axis=1)
    if not has.any():
        return []
    scale = float(diag[rot[rot >= 0]].max()) if (rot >= 0).any() else 0.0
    out: List[Tuple[int, np.ndarray]] = []
    for i in np.nonzero(has)[0]:
        e = rot[i]
        act = np.nonzero(e >= 0)[0]
        blk = K[e[act]][:, e[act]].toarray().real
        w, v = np.linalg.eigh(0.5 * (blk + blk.T))
        for k in range(w.size):
            if w[k] <= rtol * max(scale, 1e-300):
                axis = np.zeros(3)
                axis[act] = v[:, k]
                out.append((int(dofmap.node_ids[i]), axis / np.linalg.norm(axis)))
    return out


def _eigh_singular_mass(K: np.ndarray, M: np.ndarray, rtol: float = 1e-12):
    """``K phi = w^2 M phi`` for a symmetric positive *semi*-definite M (dense): with ``M = V diag(lam) V^T``
    the directions ``V_0`` with ``lam <= rtol max(lam)`` carry no inertia and are condensed statically,
    ``K_c = K_pp - K_p0 K_00^-1 K_0p`` (``K_ab = V_a^T K V_b``), exact for massless directions; then
    ``K_c q = w^2 diag(lam_p) q`` and ``phi = V_p q - V_0 K_00^-1 K_0p q`` (M-orthonormal)."""
    import scipy.linalg as sla
    lam, V = np.linalg.eigh(0.5 * (M + M.T))
    keep = lam > rtol * max(float(lam.max()), 1e-300)
    Vp, V0 = V[:, keep], V[:, ~keep]
    Kp = K @ Vp
    Kpp, K0p, K00 = Vp.T @ Kp, V0.T @ Kp, V0.T @ K @ V0
    X = np.linalg.solve(0.5 * (K00 + K00.T), K0p) if V0.shape[1] else np.zeros((0, Vp.shape[1]))
    Kc = Kpp - K0p.T @ X
    w2, q = sla.eigh(0.5 * (Kc + Kc.T), np.diag(lam[keep]))
    return w2, Vp @ q - V0 @ (X @ q)


def natural_frequencies(K, M, nmodes: Optional[int] = None, return_modes: bool = False,
                        sparse_above: int = 1500, shift: Optional[float] = None):
    """Natural frequencies (Hz) of ``Re(K) phi = w^2 M phi``.

    Pass the undamped stiffness K0 (``assemble(..., undamped=True)``) for undamped frequencies:
    with hysteretic damping Re(K*) = K0 (1 - 2 beta^2) for the SASSI complex-modulus form.

    Small models (or ``nmodes`` None) are solved densely: DOFs without mass (e.g. shell rotations
    with lumped translational mass) are condensed statically first, which is exact for massless
    DOFs.  A mass matrix that is still singular -- massless *directions* that are not single DOFs, e.g. the
    TSHELL rotary inertia ``rho t^3/12 A/n (I - n n^T)`` of a facet in an oblique plane, which has no inertia
    about its normal n -- is handled by :func:`_eigh_singular_mass` (the massless directions of M are
    condensed statically, equally exact).  Larger models with ``nmodes`` given use sparse shift-invert
    Lanczos (scipy ``eigsh``) about a small negative shift, which tolerates rigid-body modes and massless
    DOFs.  Rigid-body modes give ~0 Hz.  Mode shapes are M-orthonormal on the massive DOFs.
    """
    import scipy.linalg as sla
    import scipy.sparse.linalg as spla

    n = K.shape[0]
    if nmodes is not None and n > sparse_above:
        Ks = sp.csr_matrix(K)
        if np.iscomplexobj(Ks.data):
            Ks = Ks.real
        Ks = 0.5 * (Ks + Ks.T)
        Msp = sp.csc_matrix(M)
        Msp = 0.5 * (Msp + Msp.T)
        if shift is None:
            dk = np.abs(Ks.diagonal())
            dm = np.abs(Msp.diagonal())
            act = dm > 0
            ratio = np.median(dk[act] / dm[act]) if act.any() else 1.0
            shift = -1e-6 * ratio
        with blas_quiet():
            w2, phi = spla.eigsh(Ks.tocsc(), k=int(nmodes), M=Msp.tocsc(), sigma=shift, which="LM")
        order = np.argsort(w2)
        w2, phi = w2[order], phi[:, order]
        f = np.sqrt(np.clip(w2, 0.0, None)) / (2.0 * np.pi)
        return (f, phi) if return_modes else f

    Kd = (K.toarray() if sp.issparse(K) else np.asarray(K))
    Kd = np.real(Kd).astype(float)
    Md = (M.toarray() if sp.issparse(M) else np.asarray(M)).astype(float)
    Kd = 0.5 * (Kd + Kd.T)
    Md = 0.5 * (Md + Md.T)
    mrow = np.abs(Md).sum(axis=1)
    m = np.nonzero(mrow > 1e-14 * max(mrow.max(), 1e-300))[0]
    z = np.setdiff1d(np.arange(Kd.shape[0]), m)
    with blas_quiet():
        if z.size:
            X = np.linalg.solve(Kd[np.ix_(z, z)], Kd[np.ix_(z, m)])
            Kc = Kd[np.ix_(m, m)] - Kd[np.ix_(m, z)] @ X
        else:
            X = None
            Kc = Kd
        Kc = 0.5 * (Kc + Kc.T)
        try:
            w2, phi = sla.eigh(Kc, Md[np.ix_(m, m)])
        except np.linalg.LinAlgError:              # M singular in directions that are not single DOFs
            w2, phi = _eigh_singular_mass(Kc, Md[np.ix_(m, m)])
    f = np.sqrt(np.clip(w2, 0.0, None)) / (2.0 * np.pi)
    if nmodes is not None:
        f, phi = f[:nmodes], phi[:, :nmodes]
    if not return_modes:
        return f
    full = np.zeros((Kd.shape[0], phi.shape[1]))
    full[m] = phi
    if X is not None:
        with blas_quiet():
            full[z] = -X @ phi
    return f, full
