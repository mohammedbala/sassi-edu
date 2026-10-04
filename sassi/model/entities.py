"""Entities of the model database: nodes, coordinate systems, groups, elements, property tables,
loads and output requests (spec 08 section 1.7; ARCHITECTURE section 8).

Every entity stores the **raw values the user entered** (spec 08 section 1.6); unit conversions
and physical interpretation (ETYPE 0 classification, MSET -> M or L table, weight -> mass) happen
when the analysis files are written (D-MDL-10, D-AFW-06).

Each class has ``to_json()`` / ``from_json()`` for SAVE/RESUME (D-MDL-02) and for the canonical
state used by the round-trip test UT-03 and the model hash.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from ..conventions import ELEMENT_TYPE_NAMES

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------
DOF_NAMES = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ")
#: D labels and the DOF indices (0-based) they expand to (spec 08 section 3.2)
DOF_LABELS: Dict[str, Tuple[int, ...]] = {
    "UX": (0,), "UY": (1,), "UZ": (2,), "ROTX": (3,), "ROTY": (4,), "ROTZ": (5,),
    "DISP": (0, 1, 2), "ROT": (3, 4, 5), "ALL": (0, 1, 2, 3, 4, 5),
}
#: INT codes (D-MDL-07): only code 0 affects the solution
INT_CODES = {0: "interaction", 1: "intermediate", 2: "interface", 3: "internal"}
INT_LETTERS = {0: "I", 1: "M", 2: "F", 3: "N"}

#: maximum node count of ``E`` per group type (spec 08 section 4.1)
MAX_ELEMENT_NODES = {1: 8, 2: 3, 3: 4, 4: 4, 5: 4, 7: 2, 9: 3}
#: minimum number of distinct-slot nodes per type (CHECK Error 7 threshold)
MIN_ELEMENT_NODES = {1: 8, 2: 3, 3: 3, 4: 3, 5: 3, 7: 2, 9: 2}
#: DOFs (1-based labels) defined by each element type at its nodes (spec 08 section 3.2 table)
ELEMENT_DOFS = {1: (1, 2, 3), 2: (1, 2, 3, 4, 5, 6), 3: (1, 2, 3, 4, 5, 6), 4: (1, 3),
                5: (1, 2, 3, 4, 5, 6), 7: (1, 2, 3, 4, 5, 6), 9: (1, 2, 3, 4, 5, 6)}


def type_name(code: int) -> str:
    return ELEMENT_TYPE_NAMES.get(code, f"TYPE{code}")


# --------------------------------------------------------------------------------------
# Nodes and coordinate systems
# --------------------------------------------------------------------------------------
@dataclass
class Node:
    """A node, stored in the coordinate system that was active when it was defined (D-MDL-04).

    ``fix`` holds the six fixity codes of D (0 free, 1 fixed) in the order UX UY UZ ROTX ROTY
    ROTZ; ``flags`` holds the INT codes that are set (0 interaction, 1 intermediate,
    2 interface, 3 internal).
    """
    id: int
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    csys: int = 0
    fix: List[int] = field(default_factory=lambda: [0] * 6)
    flags: Set[int] = field(default_factory=set)

    @property
    def xyz(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_interaction(self) -> bool:
        return 0 in self.flags

    def to_json(self) -> list:
        return [self.id, self.x, self.y, self.z, self.csys, list(self.fix), sorted(self.flags)]

    @classmethod
    def from_json(cls, d: Sequence) -> "Node":
        return cls(int(d[0]), float(d[1]), float(d[2]), float(d[3]), int(d[4]), [int(v) for v in d[5]],
                   set(int(v) for v in d[6]))


@dataclass
class CoordSys:
    """Local Cartesian coordinate system (LOC or LOCAL), D-MDL-04.

    ``origin`` is in global coordinates; ``R`` (3x3) has the local axes as **columns**.
    ``kind`` records how it was defined so WRITE can reproduce it exactly:

    * ``'LOC'``: ``params = [x0, y0, z0, txy, tyz, txz]``;
    * ``'LOCAL'``: ``nodes = (n1, n2, n3)`` and ``points`` = their global coordinates when the
      system was defined (the system does not follow later node moves).
    """
    id: int
    origin: np.ndarray
    R: np.ndarray
    kind: str = "LOC"
    type: int = 0
    params: List[float] = field(default_factory=list)
    nodes: Tuple[int, ...] = ()
    points: List[List[float]] = field(default_factory=list)

    def to_global(self, xyz) -> np.ndarray:
        return np.asarray(self.origin, float) + np.asarray(xyz, float) @ np.asarray(self.R, float).T

    def to_local(self, XYZ) -> np.ndarray:
        return (np.asarray(XYZ, float) - np.asarray(self.origin, float)) @ np.asarray(self.R, float)

    def to_json(self) -> dict:
        return dict(id=self.id, kind=self.kind, type=self.type, origin=[float(v) for v in self.origin],
                    R=[[float(v) for v in row] for row in np.asarray(self.R)],
                    params=[float(v) for v in self.params], nodes=list(self.nodes),
                    points=[[float(v) for v in p] for p in self.points])

    @classmethod
    def from_json(cls, d: dict) -> "CoordSys":
        return cls(id=int(d["id"]), origin=np.asarray(d["origin"], float), R=np.asarray(d["R"], float),
                   kind=d.get("kind", "LOC"), type=int(d.get("type", 0)), params=list(d.get("params", [])),
                   nodes=tuple(int(v) for v in d.get("nodes", [])), points=[list(p) for p in d.get("points", [])])


# --------------------------------------------------------------------------------------
# Groups and elements
# --------------------------------------------------------------------------------------
@dataclass
class Element:
    """An element of a group (spec 08 sections 5.5-5.36).

    ``nodes`` are the node numbers as given in ``E`` (trailing blanks removed; a 0 means "not
    given").  Attributes: ``mat`` (MSET/MACT: M index, or L index for excavated SOLID/PLANE --
    interpreted at CHECK/AFWRITE, D-MDL-10), ``prop`` (RSET/RACT: R, SC or MX index by group
    type), ``etype`` (ETYPE 0 implicit, 1 structure, 2 excavated soil / buried shell), ``eint``
    (EINT), ``thick`` (THICK), ``ki``/``kj`` (beam end releases P1 P2 P3 M1 M2 M3, 1 = released).
    """
    id: int
    nodes: List[int]
    mat: int = 1
    prop: int = 1
    etype: int = 0
    eint: int = 0
    thick: float = 0.0
    ki: List[int] = field(default_factory=lambda: [0] * 6)
    kj: List[int] = field(default_factory=lambda: [0] * 6)

    def copy(self, new_id: Optional[int] = None, nodes: Optional[List[int]] = None) -> "Element":
        return Element(self.id if new_id is None else new_id, list(self.nodes if nodes is None else nodes),
                       self.mat, self.prop, self.etype, self.eint, self.thick, list(self.ki), list(self.kj))

    def to_json(self) -> list:
        return [self.id, list(self.nodes), self.mat, self.prop, self.etype, self.eint, self.thick,
                list(self.ki), list(self.kj)]

    @classmethod
    def from_json(cls, d: Sequence) -> "Element":
        return cls(int(d[0]), [int(v) for v in d[1]], int(d[2]), int(d[3]), int(d[4]), int(d[5]), float(d[6]),
                   [int(v) for v in d[7]], [int(v) for v in d[8]])


@dataclass
class Group:
    """An element group: one element type per group, element numbers from 1 (spec 08 section 4.3)."""
    id: int
    type: int
    title: str = ""
    elements: Dict[int, Element] = field(default_factory=dict)

    @property
    def type_name(self) -> str:
        return type_name(self.type)

    def sorted_elements(self) -> List[Element]:
        return [self.elements[k] for k in sorted(self.elements)]

    def to_json(self) -> dict:
        return dict(id=self.id, type=self.type, title=self.title,
                    elements=[self.elements[k].to_json() for k in sorted(self.elements)])

    @classmethod
    def from_json(cls, d: dict) -> "Group":
        g = cls(int(d["id"]), int(d["type"]), d.get("title", ""))
        for e in d.get("elements", []):
            el = Element.from_json(e)
            g.elements[el.id] = el
        return g


# --------------------------------------------------------------------------------------
# Property tables
# --------------------------------------------------------------------------------------
@dataclass
class Material:
    """M table entry: ``M,<nm>,<val1>,<val2>,<weight>,<pdamp>,<sdamp>,<type>`` (raw values)."""
    id: int
    val1: float
    val2: float
    weight: float
    pdamp: float
    sdamp: float
    mtype: int = 1

    FIELDS = ("val1", "val2", "weight", "pdamp", "sdamp", "mtype")

    def constants(self, gravity: float):
        """Derived real elastic constants (:func:`sassi.model.materials.elastic_constants`)."""
        from .materials import elastic_constants
        return elastic_constants(self.mtype, self.val1, self.val2, self.weight, gravity)

    def to_json(self) -> list:
        return [self.id, self.val1, self.val2, self.weight, self.pdamp, self.sdamp, self.mtype]

    @classmethod
    def from_json(cls, d: Sequence) -> "Material":
        return cls(int(d[0]), float(d[1]), float(d[2]), float(d[3]), float(d[4]), float(d[5]), int(d[6]))


@dataclass
class SoilLayer:
    """L table entry: ``L,<nm>,<thick>,<weight>,<pveloc>,<sveloc>,<pdamp>,<sdamp>``."""
    id: int
    thick: float
    weight: float
    vp: float
    vs: float
    pdamp: float
    sdamp: float

    FIELDS = ("thick", "weight", "vp", "vs", "pdamp", "sdamp")

    def constants(self, gravity: float):
        from .materials import layer_constants
        return layer_constants(self.vp, self.vs, self.weight, gravity)

    def to_json(self) -> list:
        return [self.id, self.thick, self.weight, self.vp, self.vs, self.pdamp, self.sdamp]

    @classmethod
    def from_json(cls, d: Sequence) -> "SoilLayer":
        return cls(int(d[0]), *[float(v) for v in d[1:7]])


@dataclass
class BeamSection:
    """R table entry: ``R,<nm>,<axial>,<shear2>,<shear3>,<tors>,<flex2>,<flex3>`` (A, As2, As3, J, I2, I3)."""
    id: int
    axial: float
    shear2: float
    shear3: float
    tors: float
    flex2: float
    flex3: float

    FIELDS = ("axial", "shear2", "shear3", "tors", "flex2", "flex3")

    def to_json(self) -> list:
        return [self.id, self.axial, self.shear2, self.shear3, self.tors, self.flex2, self.flex3]

    @classmethod
    def from_json(cls, d: Sequence) -> "BeamSection":
        return cls(int(d[0]), *[float(v) for v in d[1:7]])


@dataclass
class SpringProp:
    """SC table entry: six uncoupled global spring constants and a damping ratio."""
    id: int
    scx: float
    scy: float
    scz: float
    scxx: float
    scyy: float
    sczz: float
    damp: float

    FIELDS = ("scx", "scy", "scz", "scxx", "scyy", "sczz", "damp")

    @property
    def k(self) -> Tuple[float, ...]:
        return (self.scx, self.scy, self.scz, self.scxx, self.scyy, self.sczz)

    def to_json(self) -> list:
        return [self.id, self.scx, self.scy, self.scz, self.scxx, self.scyy, self.sczz, self.damp]

    @classmethod
    def from_json(cls, d: Sequence) -> "SpringProp":
        return cls(int(d[0]), *[float(v) for v in d[1:8]])


MATRIX_KINDS = ("R", "I", "M")   # MXR real stiffness, MXI imaginary stiffness, MXM mass/weight


@dataclass
class MatrixProp:
    """GENERAL element matrix property (MXR / MXI / MXM), spec 08 sections 5.26-5.29.

    ``rows[kind][r]`` holds the ``13 - r`` upper-triangle terms of row ``r`` (1..12) as entered,
    ``A[r, r], A[r, r+1], ..., A[r, 12]``; rows never entered are zero.
    """
    id: int
    rows: Dict[str, Dict[int, List[float]]] = field(default_factory=lambda: {k: {} for k in MATRIX_KINDS})

    def set_row(self, kind: str, row: int, terms: Sequence[float]) -> None:
        n = 13 - row
        t = [float(v) for v in terms][:n]
        t += [0.0] * (n - len(t))
        self.rows.setdefault(kind, {})[row] = t

    def full(self, kind: str) -> np.ndarray:
        """Symmetric 12x12 matrix of ``kind`` ('R', 'I' or 'M'), lower triangle by symmetry."""
        A = np.zeros((12, 12))
        for r, terms in self.rows.get(kind, {}).items():
            i = r - 1
            A[i, i:] = terms
        return A + np.triu(A, 1).T

    def stiffness(self) -> np.ndarray:
        """Complex stiffness ``K* = K_R + i K_I`` as entered (D-CNV-04: GENERAL damping as entered)."""
        return self.full("R") + 1j * self.full("I")

    def mass(self, weight_units: int, gravity: float) -> np.ndarray:
        """Mass matrix in mass units; MXM terms are divided by g when MOPT ``<matrix>`` = 1."""
        from .materials import to_mass
        return to_mass(self.full("M"), weight_units, gravity)

    def to_json(self) -> dict:
        return dict(id=self.id, rows={k: {str(r): list(v) for r, v in sorted(self.rows.get(k, {}).items())}
                                      for k in MATRIX_KINDS})

    @classmethod
    def from_json(cls, d: dict) -> "MatrixProp":
        p = cls(int(d["id"]))
        for k in MATRIX_KINDS:
            for r, v in d.get("rows", {}).get(k, {}).items():
                p.rows[k][int(r)] = [float(x) for x in v]
        return p


# --------------------------------------------------------------------------------------
# Loads and output requests
# --------------------------------------------------------------------------------------
@dataclass
class NodalLoad:
    """F or MM entry: three factors (global X, Y, Z) and their arrival times (spec 08 section 7.1).

    With the e^{+i w t} convention FORCE builds ``P(w) = a F_ref(w) exp(-i w t0)`` per DOF.
    """
    node: int
    factor: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    arrival: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])

    def to_json(self) -> list:
        return [self.node, list(self.factor), list(self.arrival)]

    @classmethod
    def from_json(cls, d: Sequence) -> "NodalLoad":
        return cls(int(d[0]), [float(v) for v in d[1]], [float(v) for v in d[2]])


@dataclass
class NodalRequest:
    """NOUT request: ``NOUT,<dir>,<code1..6>,<node list>`` (MOTION output, spec 07 section 9.2.23)."""
    dir: int
    codes: List[int]
    nodes: List[int]

    def to_json(self) -> list:
        return [self.dir, list(self.codes), list(self.nodes)]

    @classmethod
    def from_json(cls, d: Sequence) -> "NodalRequest":
        return cls(int(d[0]), [int(v) for v in d[1]], [int(v) for v in d[2]])


@dataclass
class ElementRequest:
    """EOUT request: ``EOUT,<code1..12>,<group>,<element list>`` (STRESS output)."""
    codes: List[int]
    group: int
    elements: List[int]

    def to_json(self) -> list:
        return [list(self.codes), self.group, list(self.elements)]

    @classmethod
    def from_json(cls, d: Sequence) -> "ElementRequest":
        return cls([int(v) for v in d[0]], int(d[1]), [int(v) for v in d[2]])


@dataclass
class RelDispRequest:
    """RDND request: ``RDND,<node>,<X>,<Y>,<Z>,<XX>,<YY>,<ZZ>`` (flag >= 1 = on)."""
    node: int
    flags: List[int]

    def to_json(self) -> list:
        return [self.node, list(self.flags)]

    @classmethod
    def from_json(cls, d: Sequence) -> "RelDispRequest":
        return cls(int(d[0]), [int(v) for v in d[1]])
