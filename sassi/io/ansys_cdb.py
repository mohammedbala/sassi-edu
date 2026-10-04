"""ANSYS ``.cdb`` import: the CONVERT,ANSYS converter (requirements 3.4.K, D-ANS-01..05, D-ANS-07;
spec 04 sections 11.1-11.6; spec 09 section 6.4).

Two steps, kept separate so that each can be tested on its own:

1. :func:`read_cdb` parses the file into a :class:`CdbModel` (plain ANSYS data: nodes, element
   types and KEYOPTs, elements, MPDATA materials, real constant sets, sections, constraints and
   components).  It reads the *blocked* CDWRITE format (NBLOCK, EBLOCK, RLBLOCK, CMBLOCK,
   SECBLOCK) honouring the Fortran format line of every block, and it also reads the plain
   APDL commands of the same data (N, E/EN/EMORE with TYPE/MAT/REAL/SECNUM, MP, R/RMORE, D),
   so the APDL file written by the ``ANSYS`` command (:mod:`sassi.io.apdl`) can be read back
   without ANSYS (round-trip tests).
2. :func:`convert_cdb` maps the ANSYS data onto a SASSI-EDU :class:`~sassi.model.SSIModel`.

What is converted (D-ANS-01, spec 04 11.3) -- everything else is reported:

==============================  ================================================================
ANSYS                           SASSI-EDU
==============================  ================================================================
SOLID45/65/185 (186/95, 187:    SOLID (corner nodes; 20/10-node midside nodes dropped with a
corner nodes)                   warning; degenerate prisms/pyramids/tetrahedra keep their repeats)
SHELL63 (R: TK) / SHELL181      SHELL (thickness = THICK attribute; Kirchhoff facet; SHELL181's
(section) / SHELL281 (corners)  Mindlin behaviour is not carried over)
BEAM4/44 (R), BEAM188/189,      BEAMS (section -> R table by matching the ANSYS and SASSI local
PIPE288 (sections)              frames, D-ANS-03; an orientation node K is created when needed)
LINK8/180                       BEAMS with axial stiffness only (I2 = I3 = 0, tiny J); rotations
                                fixed only at nodes that no other element stiffens in rotation
COMBIN14                        SPRING (global axis) or 2-node GENERAL k e e^T (oblique, D-ANS-04)
MASS21                          MT / MR with MUNITS 0 (mass units, D-ANS-05)
MATRIX27 (stiffness / mass)     GENERAL (MXR / MXM, global axes)
PLANE42/182 (X-Y plane)         PLANE (the model is rotated into the SASSI X-Z plane, HOUSE dim 1)
MPDATA EX, PRXY/NUXY, DENS,     M,<mat>,EX,nu,DENS*g,beta,beta,1 with beta = DMPR/2 (D-ANS-01)
DMPR (DMPS)                     (DMPS = constant structural damping coefficient g: beta = g/2)
D (value 0)                     fixities (nodal coordinate rotations are reported, not applied)
CMBLOCK node component SSI_INT  interaction nodes (INT code 0) -- the component the export writes
==============================  ================================================================

Beam axes (D-ANS-03).  ANSYS (BEAM4/44/188/189): element x from I to J; the orientation node K
lies in the element x-z plane with z pointing towards K; without K the element y axis is
parallel to the global X-Y plane (global Y for a vertical member) and BEAM4's THETA turns y
towards z.  SASSI (spec 08 4.5): axis 1 from I to J, axis 2 towards its K node, axis 3 = 1 x 2.
The converter computes the ANSYS frame of every element and places the SASSI K node so that
SASSI axis 2 is ANSYS z (axis 3 = -y); then ``I2`` (about axis 2) = Izz, ``I3`` = Iyy, the shear
area along axis 2 is the ANSYS z shear area, and the releases map as P1 P2 P3 M1 M2 M3 =
UX UZ UY ROTX ROTZ ROTY.  An ASEC section with a product of inertia Iyz != 0 is rotated to its
principal axes (axis 2 along the principal direction nearest to z, assuming ANSYS
``Iyz = int(y z dA)``).

BEAM44 with the same component released at both ends (a pin-ended brace, KEYOPT(7) = KEYOPT(8)
= 11) is CHECK Error 10 in SASSI; it is mapped onto an equivalent member (a moment released at
both ends -> zero inertia of that bending plane; a force or torque -> released at I only).

Reporting (spec 04 11.1): every non-zero KEYOPT, real constant, SECCONTROL value or MPDATA label of
a supported element that changes the ANSYS stiffness or mass and is not converted gives a warning;
added masses (BEAM4/LINK180/BEAM188 ADDMAS, SHELL63 ADMSUA, shell SECCONTROL) are lumped like the
SASSI masses, with a warning.  The reader does not evaluate APDL: commands with parameters or
expressions are counted and quoted, ``D`` on a node component is expanded.

Units: ANSYS has none; ``weight = DENS * g`` with the CONVERT ``<gravity>`` (32.2 ft/s^2,
386.4 in/s^2, 9.81 m/s^2), which also becomes the model GRAVITY (D-ANS-01).
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..model import SSIModel
from ..model.entities import (BeamSection, Element, Group, Material, MatrixProp, SpringProp)
from ..model.geometry import cosd, sind
from ..model.materials import section_rectangle

PathLike = Union[str, Path]

DOF_LABELS = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ")

# ======================================================================================
# Element-type catalogue
# ======================================================================================
#: ANSYS element number -> (family, display name).  Families: SOLID, SHELL, BEAM, LINK, SPRING,
#: MASS, MATRIX, PLANE.  Anything else is reported as unsupported (spec 04 11.1).
ANSYS_ELEMENTS: Dict[int, Tuple[str, str]] = {
    45: ("SOLID", "SOLID45"), 65: ("SOLID", "SOLID65"), 185: ("SOLID", "SOLID185"),
    186: ("SOLID", "SOLID186"), 95: ("SOLID", "SOLID95"), 187: ("SOLID", "SOLID187"),
    63: ("SHELL", "SHELL63"), 181: ("SHELL", "SHELL181"), 281: ("SHELL", "SHELL281"),
    4: ("BEAM", "BEAM4"), 44: ("BEAM", "BEAM44"), 188: ("BEAM", "BEAM188"), 189: ("BEAM", "BEAM189"),
    288: ("BEAM", "PIPE288"),
    8: ("LINK", "LINK8"), 180: ("LINK", "LINK180"),
    14: ("SPRING", "COMBIN14"),
    21: ("MASS", "MASS21"),
    27: ("MATRIX", "MATRIX27"),
    42: ("PLANE", "PLANE42"), 182: ("PLANE", "PLANE182"),
}
#: element names known to be unsupported (better messages); everything else is "not supported"
_KNOWN_UNSUPPORTED = {
    50: "MATRIX50 (super-element; Option AA, P2)",
    184: "MPC184 (display-only spring in Option AA, P2)", 16: "PIPE16 (display-only beam in Option AA, P2)",
    18: "PIPE18 (display-only beam in Option AA, P2)", 80: "FLUID80 (display-only solid in Option AA, P2)",
    154: "SURF154", 200: "MESH200", 170: "TARGE170", 173: "CONTA173", 174: "CONTA174",
    3: "BEAM3", 1: "LINK1", 41: "SHELL41", 43: "SHELL43", 93: "SHELL93", 92: "SOLID92", 285: "SOLID285",
}

#: KEYOPTs per supported element: (interpreted by the converter, with their own messages;
#: without effect on the converted linear-elastic model -- output, stress-only, layer-storage or
#: large-deflection controls, or an ANSYS integration choice of an element whose formulation is
#: replaced anyway).  Any other non-zero KEYOPT is reported (spec 04 11.1: "a warning for every
#: unsupported ... KEYOPT").  When in doubt a KEYOPT is *not* listed as without effect.
_KEYOPTS: Dict[int, Tuple[frozenset, frozenset]] = {
    45: (frozenset({1, 2}), frozenset({4, 5, 6})),
    65: (frozenset({1}), frozenset({5, 6, 7})),
    185: (frozenset({2, 3, 6}), frozenset({8})),
    186: (frozenset({3, 6}), frozenset({2, 8})),
    187: (frozenset({6}), frozenset()),
    95: (frozenset(), frozenset({5, 6})),
    63: (frozenset({1}), frozenset({5, 6})),
    181: (frozenset({1}), frozenset({3, 8})),
    281: (frozenset({1}), frozenset({8})),
    4: (frozenset(), frozenset({6, 9})),
    44: (frozenset({7, 8}), frozenset({6, 9})),
    188: (frozenset({1, 3}), frozenset({2, 4, 6, 7, 9, 11, 12, 15})),
    189: (frozenset({1}), frozenset({2, 4, 6, 7, 9, 11, 12, 15})),
    288: (frozenset({1, 3}), frozenset({6, 7, 15})),
    8: (frozenset(), frozenset()),
    180: (frozenset({3}), frozenset({2})),
    14: (frozenset({2, 3}), frozenset({1})),
    21: (frozenset({1, 2, 3}), frozenset()),
    27: (frozenset({2, 3}), frozenset({4})),
    42: (frozenset({2, 3}), frozenset({1, 5, 6})),
    182: (frozenset({1, 3, 6}), frozenset()),
}

#: CDWRITE boilerplate and set-up commands that carry nothing to convert (silently skipped)
_SILENT = {
    "/COM", "/PREP7", "/NOPR", "/GOPR", "/GO", "FINISH", "*IF", "*ELSE", "*ENDIF", "*SET", "*DIM",
    "NUMOFF", "DOF", "EXTOPT", "TREF", "IRLF", "BFUNIF", "OMEGA", "DOMEGA", "CGLOC", "CGOMEGA",
    "DCGOMG", "KUSE", "TIME", "CRPLIM", "NCNV", "NEQIT", "ERESX", "MPTEMP", "/SOLU", "/POST1", "/CLEAR",
    "/BATCH", "/CONFIG", "/OUTPUT", "/FILNAME", "/SHOW", "/VIEW", "/PNUM", "/NUMBER", "/REPLOT", "/ESHAPE",
    "/AUTO", "/USER", "/DIST", "/FOCUS", "/ANG", "SECCONTROL", "MPTRES", "MPCOPY", "ALLSEL", "NSEL", "ESEL",
    "TBTEMP", "/EOF", "SHPP", "ANTYPE", "LUMPM", "/INPUT", "CMEDIT", "CMSEL", "CSYS", "LOCAL", "CLOCAL",
    "LSYM", "CS", "/TITLE", "/UNITS", "ACEL", "ALPHAD", "BETAD", "DMPRAT", "DMPSTR", "ENDRELEASE",
}
#: commands whose data is reported as not converted (with the reason)
_NOT_CONVERTED = {
    "CE": "constraint equations (CE) are not converted",
    "CP": "coupled DOF sets (CP) are not converted",
    "CERIG": "rigid regions (CERIG) are not converted",
    "RBE3": "RBE3 constraints are not converted",
    "CEINTF": "constraint equations (CEINTF) are not converted",
    "F": "nodal forces (F) are not converted (define F/MM loads in SASSI-EDU for vibration analyses)",
    "SF": "surface loads (SF) are not converted",
    "SFE": "surface loads (SFE) are not converted",
    "SFBEAM": "beam surface loads (SFBEAM) are not converted",
    "BF": "body loads (BF) are not converted",
    "BFE": "body loads (BFE) are not converted",
    "TB": "nonlinear / tabular material data (TB) is not converted",
    "TBDATA": "nonlinear / tabular material data (TBDATA) is not converted",
    "IC": "initial conditions (IC) are not converted",
}


def element_name(ename: int) -> str:
    """ANSYS element display name (``SOLID185``) of an element number."""
    if ename in ANSYS_ELEMENTS:
        return ANSYS_ELEMENTS[ename][1]
    if _KNOWN_UNSUPPORTED.get(ename):
        return _KNOWN_UNSUPPORTED[ename].split()[0]
    return f"element {ename}"


# ======================================================================================
# Data classes of the parsed file
# ======================================================================================
@dataclass
class CdbElement:
    """One ANSYS element: number, attributes (TYPE, MAT, REAL, SECNUM, ESYS) and nodes."""
    num: int
    type: int
    mat: int
    real: int
    secnum: int
    esys: int
    nodes: List[int]


@dataclass
class CdbSection:
    """SECTYPE + SECDATA / SECBLOCK / SECOFFSET / SECCONTROL of one section id."""
    id: int
    type: str = ""            # BEAM, SHELL, PIPE, LINK, ...
    subtype: str = ""         # RECT, ASEC, CSOLID, CTUBE, ...
    name: str = ""
    data: List[float] = field(default_factory=list)            # beam/link/pipe SECDATA values
    layers: List[Tuple[float, int, float, int]] = field(default_factory=list)   # shell (TK, MAT, THETA, NUMPT)
    offset: List[str] = field(default_factory=list)            # SECOFFSET fields
    control: List[float] = field(default_factory=list)         # SECCONTROL values
    block: bool = False                                        # beam SECBLOCK (mesh section)


@dataclass
class CdbModel:
    """Everything read from a ``.cdb`` (or APDL) file, in ANSYS terms."""
    nodes: Dict[int, np.ndarray] = field(default_factory=dict)          # global X, Y, Z
    node_rot: Dict[int, Tuple[float, float, float]] = field(default_factory=dict)   # non-zero THXY, THYZ, THZX
    etypes: Dict[int, int] = field(default_factory=dict)                # itype -> element number (185)
    keyopts: Dict[int, Dict[int, int]] = field(default_factory=dict)    # itype -> {k: value}
    elements: List[CdbElement] = field(default_factory=list)            # file order
    materials: Dict[int, Dict[str, float]] = field(default_factory=dict)  # mat -> {label: value at T1}
    temp_dependent: Dict[int, List[str]] = field(default_factory=dict)  # labels given at several temperatures
    reals: Dict[int, List[float]] = field(default_factory=dict)         # real set -> values
    sections: Dict[int, CdbSection] = field(default_factory=dict)
    constraints: List[Tuple[int, str, float, float]] = field(default_factory=list)  # node, label, value, value2
    components: Dict[str, Tuple[str, List[int]]] = field(default_factory=dict)       # name -> (NODE|ELEM, ids)
    title: str = ""
    units: str = ""
    skipped: Dict[str, int] = field(default_factory=dict)               # not-converted commands -> count
    unknown: Dict[str, int] = field(default_factory=dict)               # unknown commands -> count
    notes: List[str] = field(default_factory=list)                      # reader warnings
    unevaluated: List[str] = field(default_factory=list)                # first lines with APDL parameters
    comp_constraints: List[Tuple[str, str, float, float]] = field(default_factory=list)  # D on a component
    global_damping: Dict[str, float] = field(default_factory=dict)      # ALPHAD BETAD DMPRAT DMPSTR != 0
    source: str = ""

    def keyopt(self, itype: int, k: int) -> int:
        return int(self.keyopts.get(itype, {}).get(k, 0))

    def ename(self, itype: int) -> int:
        return int(self.etypes.get(itype, 0))


# ======================================================================================
# Fortran formats of the blocked commands
# ======================================================================================
_ITEM_RE = re.compile(r"^(\d*)([IEGFDX])(\d*)(?:\.(\d+))?(?:E(\d+))?$")


def parse_fortran_format(fmt: str) -> List[Tuple[str, int]]:
    """Field list of a simple Fortran format line such as ``(3i9,6e21.13e3)`` or ``(2i8,6g16.9)``.

    Returns ``[(kind, width), ...]`` with the repeat counts expanded; kind ``'i'`` (integer),
    ``'r'`` (real: E, G, F, D) or ``'x'`` (skip).  Raises ``ValueError`` for anything else.
    """
    s = fmt.strip()
    if not (s.startswith("(") and ")" in s):
        raise ValueError(f"not a format line: {fmt!r}")
    s = s[1:s.rindex(")")].replace(" ", "").upper()
    out: List[Tuple[str, int]] = []
    for item in s.split(","):
        if not item:
            continue
        m = _ITEM_RE.match(item)
        if not m:
            raise ValueError(f"unsupported format item {item!r} in {fmt!r}")
        rep = int(m.group(1)) if m.group(1) else 1
        kind = m.group(2)
        if kind == "X":                       # nX: skip n characters
            out.append(("x", rep))
            continue
        width = int(m.group(3)) if m.group(3) else 0
        if width <= 0:
            raise ValueError(f"format item {item!r} has no width")
        out.extend([("i" if kind == "I" else "r", width)] * rep)
    return out


#: a Fortran double-precision literal such as ``1.5D+03`` (the only tokens whose D becomes E)
_FORTRAN_D = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)[dD][+-]?\d+$")


def _num(tok: str) -> float:
    """Number of a field; ``1.5D+03`` (Fortran double) is accepted.  Anything else that is not a
    number (an APDL parameter name or expression) raises ``ValueError`` with the token unchanged."""
    t = tok.strip()
    if not t:
        return 0.0
    if _FORTRAN_D.match(t):
        t = t.replace("D", "E").replace("d", "e")
    return float(t)


def split_fixed(line: str, fields: Sequence[Tuple[str, int]]) -> List[str]:
    """Cut ``line`` into the fixed-width fields of a parsed format (short lines give fewer fields)."""
    out: List[str] = []
    pos = 0
    text = line.rstrip("\r\n")
    for kind, w in fields:
        if pos >= len(text):
            break
        chunk = text[pos:pos + w]
        pos += w
        if kind == "x":
            continue
        out.append(chunk.strip())
    while out and out[-1] == "":
        out.pop()
    return out


def read_fixed_numbers(line: str, fields: Optional[Sequence[Tuple[str, int]]]) -> List[float]:
    """Numbers of one data line: fixed-width per the block's format, whitespace as a fallback.

    Fixed widths matter because adjacent ``e21.13e3`` fields may touch
    (``-1.0000000000000E+000-2.0...``).  Lines of other writers that ignore the widths are read by
    whitespace when the fixed split does not give numbers.
    """
    if fields:
        toks = split_fixed(line, fields)
        try:
            return [_num(t) for t in toks]
        except ValueError:
            pass
    return [_num(t) for t in line.replace(",", " ").split()]


# ======================================================================================
# Reader
# ======================================================================================
def _strip_comment(line: str) -> str:
    i = line.find("!")
    return line if i < 0 else line[:i]


def _fields(line: str) -> List[str]:
    return [f.strip() for f in _strip_comment(line).split(",")]


def _fnum(fs: Sequence[str], k: int, default: float = 0.0) -> float:
    if k >= len(fs) or not fs[k].strip():
        return default
    return _num(fs[k])


def _fint(fs: Sequence[str], k: int, default: int = 0) -> int:
    if k >= len(fs) or not fs[k].strip():
        return default
    return int(round(_num(fs[k])))


#: ``skipped`` key of commands whose numeric fields hold APDL parameters or expressions
PARAM_KEY = "APDL parameters/expressions"

#: commands the reader handles, reachable by their 4-character APDL abbreviation (``KEYO``, ``SECD``)
_HANDLED = ("NBLOCK", "EBLOCK", "RLBLOCK", "CMBLOCK", "SECBLOCK", "KEYOPT", "MPDATA", "MPTEMP", "SECTYPE",
            "SECDATA", "SECOFFSET", "SECCONTROL", "SECNUM", "RMORE", "EMORE", "ENDRELEASE", "ALPHAD", "BETAD",
            "DMPRAT", "DMPSTR")
_ABBREV = {name[:4]: name for name in _HANDLED}


def _ename_from_token(tok: str) -> int:
    """``185`` or ``SOLID185`` -> 185."""
    t = tok.strip().upper()
    m = re.search(r"(\d+)$", t)
    if not m:
        raise ValueError(f"element name {tok!r} has no number")
    return int(m.group(1))


class _Reader:
    def __init__(self, lines: List[str], source: str = ""):
        self.lines = lines
        self.i = 0
        self.m = CdbModel(source=source)
        # state of the plain APDL element commands
        self.type = 1
        self.mat = 1
        self.real = 1
        self.secnum = 0
        self.esys = 0
        self.csys = 0
        self.last_real: Optional[int] = None
        self.last_elem: Optional[CdbElement] = None
        self.last_section: Optional[CdbSection] = None
        self.max_elem = 0
        self.mptemp_count = 1

    # ------------------------------------------------------------------ helpers
    def _next_line(self) -> Optional[str]:
        if self.i >= len(self.lines):
            return None
        line = self.lines[self.i]
        self.i += 1
        return line

    def _peek(self) -> Optional[str]:
        return self.lines[self.i] if self.i < len(self.lines) else None

    def _format(self) -> Optional[List[Tuple[str, int]]]:
        """Read the format line that follows a block header (None when absent or unreadable)."""
        nxt = self._peek()
        if nxt is None or not nxt.strip().startswith("("):
            return None
        self.i += 1
        try:
            return parse_fortran_format(nxt)
        except ValueError as exc:
            self.m.notes.append(f"{exc}; whitespace-separated fields assumed")
            return None

    def _count(self, table: Dict[str, int], key: str) -> None:
        table[key] = table.get(key, 0) + 1

    def _data(self, line: str, fmt) -> Optional[List[float]]:
        """Numbers of a block data line, or None (with a note) when the line is not numeric."""
        try:
            return read_fixed_numbers(line, fmt)
        except ValueError:
            self.m.notes.append(f"malformed block data line {line.strip()[:50]!r} skipped")
            return None

    def _add_element(self, el: CdbElement) -> None:
        self.m.elements.append(el)
        self.max_elem = max(self.max_elem, el.num)
        self.last_elem = el

    # ------------------------------------------------------------------ blocks
    def nblock(self, fs: List[str]) -> None:
        """NBLOCK,NUMFIELD,Solkey,NDMAX,NDSEL / format / node lines / ``N,R5.3,LOC,-1,``."""
        fmt = self._format()
        nints = sum(1 for k, _ in fmt if k == "i") if fmt else 3
        while True:
            line = self._peek()
            if line is None:
                break
            s = line.strip()
            if not s:
                self.i += 1
                continue
            if s.upper().startswith("N,") or s == "-1":
                self.i += 1                      # terminator 'N,R5.3,LOC,     -1,'
                break
            if s[0].isalpha() or s[0] in "/*":
                break                            # next command (block without a terminator)
            self.i += 1
            vals = self._data(line, fmt)
            if not vals:
                continue
            nid = int(round(vals[0]))
            coords = list(vals[nints:nints + 3]) + [0.0] * 3
            rot = list(vals[nints + 3:nints + 6]) + [0.0] * 3
            self.m.nodes[nid] = np.array(coords[:3], float)
            if any(abs(r) > 0.0 for r in rot[:3]):
                self.m.node_rot[nid] = (rot[0], rot[1], rot[2])

    def eblock(self, fs: List[str]) -> None:
        """EBLOCK,NUM_NODES,Solkey,NDMAX,NDSEL / format / element records / ``-1``.

        Solkey = SOLID (CDWRITE): mat, type, real, secnum, esys, birth/death, solid ref, shape,
        number of nodes, (unused), element number, nodes 1-8, continuation lines for more nodes.
        Solkey blank: element number, type, real, mat, esys, nodes (NUM_NODES on the first line).
        """
        solid = len(fs) > 2 and fs[2].upper().startswith("SOLID")
        num_first = _fint(fs, 1, 8)
        fmt = self._format()
        while True:
            line = self._next_line()
            if line is None:
                break
            s = line.strip()
            if not s:
                continue
            vals = self._data(line, fmt)
            if vals is None:
                continue
            if not vals or (len(vals) == 1 and int(vals[0]) == -1):
                break
            ints = [int(round(v)) for v in vals]
            if solid:
                if len(ints) < 11:
                    self.m.notes.append(f"EBLOCK record {s[:40]!r} too short; skipped")
                    continue
                mat, typ, real, sec, esys = ints[0:5]
                nn = ints[8]
                num = ints[10]
                nodes = ints[11:11 + nn]
                while len(nodes) < nn:
                    more = self._next_line()
                    if more is None:
                        break
                    nodes += [int(round(v)) for v in (self._data(more, fmt) or [])]
                nodes = nodes[:nn]
            else:
                if len(ints) < 6:
                    self.m.notes.append(f"EBLOCK record {s[:40]!r} too short; skipped")
                    continue
                num, typ, real, mat, esys = ints[0:5]
                sec = 0
                nodes = ints[5:5 + num_first]
                while nodes and nodes[-1] == 0:
                    nodes.pop()
            self._add_element(CdbElement(num, typ, mat, real, sec, esys, nodes))

    def rlblock(self, fs: List[str]) -> None:
        """RLBLOCK,NUMSETS,MAXSET,MAXITEMS,NPERLINE / (2i8,6g16.9) / (7g16.9) / sets."""
        nsets = _fint(fs, 1, 0)
        fmt1 = self._format()
        fmt2 = self._format()
        for _ in range(nsets):
            line = self._next_line()
            if line is None:
                break
            vals = self._data(line, fmt1) or []
            if len(vals) < 2:
                continue
            setno, n = int(round(vals[0])), int(round(vals[1]))
            data = list(vals[2:])
            while len(data) < n:
                more = self._next_line()
                if more is None:
                    break
                data += self._data(more, fmt2 or fmt1) or []
            self.m.reals[setno] = data[:n] + [0.0] * max(0, n - len(data))

    def cmblock(self, fs: List[str]) -> None:
        """CMBLOCK,Cname,Entity,NUMITEMS / (8i10) / items; a negative item closes a range."""
        name = fs[1].strip().upper() if len(fs) > 1 else ""
        kind = fs[2].strip().upper() if len(fs) > 2 else "NODE"
        n = _fint(fs, 3, 0)
        fmt = self._format()
        items: List[int] = []
        while len(items) < n:
            line = self._next_line()
            if line is None:
                break
            items += [int(round(v)) for v in (self._data(line, fmt) or [])]
        ids: List[int] = []
        for v in items[:n]:
            if v < 0 and ids:
                ids.extend(range(ids[-1] + 1, -v + 1))
            elif v > 0:
                ids.append(v)
        self.m.components[name] = (kind, ids)

    def secblock(self, fs: List[str]) -> None:
        """SECBLOCK,NLAYERS: shell layers ``TK MAT THETA NUMPT`` (one per line) follow."""
        sec = self.last_section
        n = _fint(fs, 1, 0)
        if sec is None:
            self.m.notes.append("SECBLOCK without a preceding SECTYPE; skipped")
        if sec is not None and sec.type != "SHELL":
            sec.block = True          # mesh-based beam section: data not usable
        fmt = self._format()
        for _ in range(n):
            line = self._next_line()
            if line is None:
                break
            vals = self._data(line, fmt) or []
            if sec is not None and sec.type == "SHELL" and vals:
                vals = vals + [0.0] * 4
                sec.layers.append((vals[0], int(round(vals[1])), vals[2], int(round(vals[3]))))

    # ------------------------------------------------------------------ commands
    def command(self, line: str) -> None:
        """Read one command.  APDL input may use parameters (``*SET,tid,4`` then ``ET,tid,170``),
        expressions (``N,2,L/2``) or names where a number is expected; the reader does not evaluate
        APDL, so such a command is counted, its line kept for the report, and reading goes on
        (spec 04 11.1: report every command that is skipped)."""
        try:
            self._dispatch(line)
        except ValueError:
            self._count(self.m.skipped, PARAM_KEY)
            if len(self.m.unevaluated) < 5:
                self.m.unevaluated.append(_strip_comment(line).strip()[:60])

    def _dispatch(self, line: str) -> None:
        fs = _fields(line)
        head = fs[0].upper()
        if not head:
            return
        if "=" in head:                       # parameter assignment (CDWRITE: _CDRDOFF=)
            return
        h = head
        if len(h) >= 4 and h not in _HANDLED and h[:4] in _ABBREV and _ABBREV[h[:4]].startswith(h):
            h = _ABBREV[h[:4]]                # APDL: commands may be abbreviated to 4 characters
        if h == "NBLOCK":
            return self.nblock(fs)
        if h == "EBLOCK":
            return self.eblock(fs)
        if h in ("RLBLOCK", "RBLOCK"):
            return self.rlblock(fs)
        if h == "CMBLOCK":
            return self.cmblock(fs)
        if h == "SECBLOCK":
            return self.secblock(fs)
        if h == "/TITLE":
            self.m.title = line.split(",", 1)[1].strip() if "," in line else ""
            return
        if h == "/UNITS":
            self.m.units = fs[1].upper() if len(fs) > 1 else ""
            return
        if h in ("ALPHAD", "BETAD", "DMPRAT", "DMPSTR"):
            v = _fnum(fs, 1)
            if v != 0.0:
                self.m.global_damping[h] = v
            return
        if h == "ACEL":
            if any(_fnum(fs, k) != 0.0 for k in (1, 2, 3)):
                self._count(self.m.skipped, "ACEL")
            return
        if h == "CSYS":
            self.csys = _fint(fs, 1, 0)
            return
        if h == "ENDRELEASE":
            self._count(self.m.skipped, "ENDRELEASE")
            return
        handler = getattr(self, "c_" + h.replace("/", "_").replace("*", "_"), None)
        if handler is not None:
            return handler(fs)
        if h in _NOT_CONVERTED:
            self._count(self.m.skipped, h)
            return
        if h in _SILENT:
            return
        self._count(self.m.unknown, h)

    def c_ET(self, fs):
        itype = _fint(fs, 1)
        try:
            self.m.etypes[itype] = _ename_from_token(fs[2]) if len(fs) > 2 else 0
        except ValueError as exc:
            self.m.notes.append(f"ET,{itype}: {exc}")
            return
        kop = self.m.keyopts.setdefault(itype, {})
        for k in range(1, 7):                  # ET,ITYPE,Ename,KOP1..KOP6
            if len(fs) > 2 + k and fs[2 + k].strip():
                kop[k] = _fint(fs, 2 + k)

    def c_KEYOPT(self, fs):
        self.m.keyopts.setdefault(_fint(fs, 1), {})[_fint(fs, 2)] = _fint(fs, 3)

    def c_N(self, fs):
        if len(fs) > 1 and fs[1].upper().startswith("R5"):
            return                               # NBLOCK terminator outside a block
        nid = _fint(fs, 1)
        if nid <= 0:
            return
        if self.csys != 0:
            self._count(self.m.skipped, "N in a local CSYS")
        self.m.nodes[nid] = np.array([_fnum(fs, 2), _fnum(fs, 3), _fnum(fs, 4)])
        rot = (_fnum(fs, 5), _fnum(fs, 6), _fnum(fs, 7))
        if any(r != 0.0 for r in rot):
            self.m.node_rot[nid] = rot

    def c_TYPE(self, fs):
        self.type = _fint(fs, 1, 1)

    def c_MAT(self, fs):
        self.mat = _fint(fs, 1, 1)

    def c_REAL(self, fs):
        self.real = _fint(fs, 1, 1)

    def c_SECNUM(self, fs):
        self.secnum = _fint(fs, 1, 0)

    c_SECN = c_SECNUM

    def c_ESYS(self, fs):
        self.esys = _fint(fs, 1, 0)

    def _elem(self, num: int, nodes: List[int]) -> None:
        while nodes and nodes[-1] == 0:
            nodes.pop()
        self._add_element(CdbElement(num, self.type, self.mat, self.real, self.secnum, self.esys, nodes))

    def c_E(self, fs):
        self._elem(self.max_elem + 1, [_fint(fs, k) for k in range(1, min(len(fs), 9))])

    def c_EN(self, fs):
        self._elem(_fint(fs, 1), [_fint(fs, k) for k in range(2, min(len(fs), 10))])

    def c_EMORE(self, fs):
        if self.last_elem is not None:
            more = [_fint(fs, k) for k in range(1, min(len(fs), 9))]
            while more and more[-1] == 0:
                more.pop()
            nodes = list(self.last_elem.nodes) + [0] * max(0, 8 - len(self.last_elem.nodes))
            self.last_elem.nodes = nodes + more

    def _mat_value(self, mat: int, lab: str, values: List[float], ntemp: int) -> None:
        lab = lab.upper()
        self.m.materials.setdefault(mat, {})[lab] = values[0] if values else 0.0
        if ntemp > 1 or len(values) > 1 and any(v != 0.0 for v in values[1:]):
            lst = self.m.temp_dependent.setdefault(mat, [])
            if lab not in lst:
                lst.append(lab)

    def c_MPTEMP(self, fs):
        if len(fs) > 1 and fs[1].upper().startswith("R5"):
            self.mptemp_count = _fint(fs, 2, 1)

    def c_MPDATA(self, fs):
        """``MPDATA,R5.0,LENGTH,Lab,MAT,STLOC,C1,...`` (CDWRITE) or ``MPDATA,Lab,MAT,STLOC,C1,...``."""
        if len(fs) > 1 and fs[1].upper().startswith("R5"):
            ntemp = _fint(fs, 2, 1)
            lab, mat, stloc, vals = fs[3], _fint(fs, 4), _fint(fs, 5, 1), fs[6:]
        else:
            ntemp = 1
            lab, mat, stloc, vals = fs[1], _fint(fs, 2), _fint(fs, 3, 1), fs[4:]
        values = [_num(v) for v in vals if v.strip()]
        if stloc == 1:
            self._mat_value(mat, lab, values, max(ntemp, self.mptemp_count))

    def c_MP(self, fs):
        """``MP,Lab,MAT,C0,C1,...`` (C0 = the value; temperature coefficients are reported)."""
        if len(fs) < 3:
            return
        values = [_num(v) for v in fs[3:] if v.strip()]
        self._mat_value(_fint(fs, 2), fs[1], values or [0.0], 1)

    def c_R(self, fs):
        setno = _fint(fs, 1)
        vals = [_fnum(fs, k) for k in range(2, 8)]
        self.m.reals[setno] = vals
        self.last_real = setno

    def c_RMORE(self, fs):
        if self.last_real is None:
            return
        cur = self.m.reals[self.last_real]
        pad = (-len(cur)) % 6
        cur.extend([0.0] * pad)
        cur.extend(_fnum(fs, k) for k in range(1, 7))

    def c_SECTYPE(self, fs):
        sid = _fint(fs, 1)
        sec = CdbSection(sid, (fs[2] if len(fs) > 2 else "").upper(), (fs[3] if len(fs) > 3 else "").upper(),
                         fs[4] if len(fs) > 4 else "")
        self.m.sections[sid] = sec
        self.last_section = sec

    def c_SECDATA(self, fs):
        sec = self.last_section
        if sec is None:
            return
        vals = [_num(v) if v.strip() else 0.0 for v in fs[1:]]
        if sec.type == "SHELL":
            v = vals + [0.0] * 4
            sec.layers.append((v[0], int(round(v[1])), v[2], int(round(v[3]))))
        else:
            sec.data.extend(vals)

    def c_SECOFFSET(self, fs):
        if self.last_section is not None:
            self.last_section.offset = [f.upper() for f in fs[1:] if f.strip()]

    def c_SECCONTROL(self, fs):
        if self.last_section is not None:
            self.last_section.control = [_num(v) if v.strip() else 0.0 for v in fs[1:]]

    def c_SECCONTROLS(self, fs):
        self.c_SECCONTROL(fs)

    def c_D(self, fs):
        """``D,NODE,Lab,VALUE,VALUE2,NEND,NINC,Lab2,...,Lab6``."""
        node = fs[1].strip().upper() if len(fs) > 1 else ""
        if node in ("ALL", "P", ""):
            self._count(self.m.skipped, "D on ALL/selected nodes")
            return
        val, val2 = _fnum(fs, 3), _fnum(fs, 4)
        labs = [fs[2] if len(fs) > 2 else ""] + [fs[k] for k in range(7, min(len(fs), 12))]
        labs = [l.strip().upper() for l in labs if l.strip()]
        if not _isnum(node):
            # D on a node component (D,SSI_BASE,ALL,0): resolved when the whole file is read,
            # because CDWRITE may write the CMBLOCK after the constraints
            for lab in labs:
                self.m.comp_constraints.append((node, lab, val, val2))
            return
        n1 = _fint(fs, 1)
        nend, ninc = _fint(fs, 5, n1), _fint(fs, 6, 1)
        nodes = range(n1, max(nend, n1) + 1, max(ninc, 1))
        for n in nodes:
            for lab in labs:
                self.m.constraints.append((n, lab, val, val2))

    def run(self) -> CdbModel:
        while True:
            line = self._next_line()
            if line is None:
                break
            s = line.strip()
            if not s or s.startswith("!"):
                continue
            self.command(s)
        self._resolve_component_constraints()
        return self.m

    def _resolve_component_constraints(self) -> None:
        """``D,<component>,...``: expand node components (an element component or an undefined
        name is reported)."""
        bad: Dict[str, int] = {}
        for name, lab, val, val2 in self.m.comp_constraints:
            kind, ids = self.m.components.get(name, ("", []))
            if not kind.startswith("NODE"):
                bad[name] = bad.get(name, 0) + 1
                continue
            for n in ids:
                self.m.constraints.append((n, lab, val, val2))
        for name, n in bad.items():
            self.m.notes.append(f"{n} D constraints on '{name}' not converted ({name} is not a node component "
                                "defined by CMBLOCK)")


def read_cdb_text(text: str, source: str = "") -> CdbModel:
    """Parse the text of a ``.cdb`` (or APDL) file; see :func:`read_cdb`."""
    return _Reader(text.splitlines(), source).run()


def read_cdb(path: PathLike) -> CdbModel:
    """Parse an ANSYS CDWRITE file (``CDWRITE,DB``; the manual's "CBD file").

    Undecodable bytes are replaced (Latin-1 titles in old files must not stop the import).
    """
    p = Path(path)
    return read_cdb_text(p.read_text(encoding="utf-8", errors="replace"), source=str(p))


# ======================================================================================
# Section properties
# ======================================================================================
@dataclass
class BeamProps:
    """Beam section in ANSYS element axes: A, Iyy (about y), Izz (about z), Iyz, J and the shear
    areas along y and z (0 = no shear deformation), plus added mass per length."""
    A: float
    Iyy: float
    Izz: float
    J: float
    Asy: float = 0.0
    Asz: float = 0.0
    Iyz: float = 0.0
    addmas: float = 0.0
    axial_only: bool = False


def rect_props(B: float, H: float) -> BeamProps:
    """ANSYS RECT section: width B along y, height H along z (Iyy = B H^3/12, Izz = H B^3/12).

    Torsion constant and shear areas as the SASSI rectangle (``sassi.model.section_rectangle``,
    Roark J, shear form factor 6/5), so a converted RECT gives exactly the R entry of the same
    rectangle in SASSI-EDU (b = B along axis 3, h = H along axis 2).
    """
    s = section_rectangle(B, H)
    return BeamProps(A=s["axial"], Iyy=s["flex3"], Izz=s["flex2"], J=s["tors"], Asy=s["shear3"], Asz=s["shear2"])


def tube_props(Do: float, Di: float) -> BeamProps:
    """Circular tube (PIPE: outer diameter Do, inner Di; CTUBE; CSOLID with Di = 0).

    A = pi (Do^2 - Di^2)/4, I = pi (Do^4 - Di^4)/64, J = 2 I; shear area 0.9 A for a solid circle
    (as ``section_circle``) and A/2 for a thin-walled tube (linear blend in Di/Do).
    """
    A = math.pi * (Do * Do - Di * Di) / 4.0
    I = math.pi * (Do ** 4 - Di ** 4) / 64.0
    r = Di / Do if Do > 0 else 0.0
    k = 0.9 * (1.0 - r) + 0.5 * r
    return BeamProps(A=A, Iyy=I, Izz=I, J=2.0 * I, Asy=k * A, Asz=k * A)


# ======================================================================================
# Beam axes (D-ANS-03)
# ======================================================================================
def ansys_beam_axes(xi, xj, xk=None, theta_deg: float = 0.0) -> Tuple[np.ndarray, float]:
    """ANSYS element axes of BEAM4/44/188/189 and LINK elements.

    Returns ``(A, L)`` with the rows of ``A`` the unit vectors ex, ey, ez (global components):

    * ex = (J - I)/L;
    * with an orientation node K: ez = the component of (K - I) normal to ex (K lies in the x-z
      plane, z towards K), ey = ez x ex;
    * without K: ey parallel to the global X-Y plane, ey = Z x ex normalised (global Y when the
      member is within a 0.01 % slope of the global Z axis), ez = ex x ey; then THETA (BEAM4)
      turns ey towards ez about ex.
    """
    xi = np.asarray(xi, float)
    d = np.asarray(xj, float) - xi
    L = float(np.linalg.norm(d))
    if L <= 0.0:
        raise ValueError("zero-length beam (I and J coincide)")
    ex = d / L
    if xk is not None:
        v = np.asarray(xk, float) - xi
        vz = v - float(v @ ex) * ex
        nz = float(np.linalg.norm(vz))
        if nz <= 1e-9 * max(L, float(np.linalg.norm(v))):
            raise ValueError("orientation node K is collinear with I-J")
        ez = vz / nz
        ey = np.cross(ez, ex)
        return np.vstack([ex, ey, ez]), L
    horiz = math.hypot(ex[0], ex[1])
    if horiz < 1e-4 * abs(ex[2]):
        y = np.array([0.0, 1.0, 0.0])
        ey = y - float(y @ ex) * ex
        ey /= np.linalg.norm(ey)
    else:
        ey = np.cross([0.0, 0.0, 1.0], ex)
        ey /= np.linalg.norm(ey)
    ez = np.cross(ex, ey)
    if theta_deg:
        c, s = cosd(theta_deg), sind(theta_deg)          # exact at multiples of 90 degrees
        ey, ez = c * ey + s * ez, -s * ey + c * ez
    return np.vstack([ex, ey, ez]), L


def sassi_section_from_ansys(p: BeamProps) -> Tuple[Dict[str, float], np.ndarray]:
    """SASSI R values from ANSYS section properties, and the SASSI axis-2 direction in ANSYS (y, z).

    Without a product of inertia SASSI axis 2 = ANSYS z and axis 3 = -y, so
    ``I2 = Izz``, ``I3 = Iyy``, ``As2 = Asz``, ``As3 = Asy``.  With ``Iyz != 0`` the principal
    axes of ``C = [[Izz, Iyz], [Iyz, Iyy]]`` (``C_ab = int a b dA`` in (y, z)) are used: axis 2 is
    the principal direction ``p2`` nearest to z, ``I3 = p2.C.p2``, ``I2`` the other eigenvalue.
    Returns ``(dict(axial, shear2, shear3, tors, flex2, flex3), p2)``.
    """
    scale = max(abs(p.Iyy), abs(p.Izz), 1e-300)
    if abs(p.Iyz) <= 1e-12 * scale:
        sec = dict(axial=p.A, shear2=p.Asz, shear3=p.Asy, tors=p.J, flex2=p.Izz, flex3=p.Iyy)
        return sec, np.array([0.0, 1.0])
    C = np.array([[p.Izz, p.Iyz], [p.Iyz, p.Iyy]])
    w, V = np.linalg.eigh(C)
    k = int(np.argmax(np.abs(V[1, :])))          # eigenvector with the largest z component
    p2 = V[:, k] * (1.0 if V[1, k] > 0 else -1.0)
    I3 = float(w[k])
    I2 = float(w[1 - k])
    sec = dict(axial=p.A, shear2=p.Asz, shear3=p.Asy, tors=p.J, flex2=I2, flex3=I3)
    return sec, p2


# ======================================================================================
# Conversion
# ======================================================================================
@dataclass
class CdbConversion:
    """Result of :func:`convert_cdb`."""
    model: SSIModel
    messages: List[Tuple[str, str]] = field(default_factory=list)     # (kind 'warning'|'info', text)
    element_map: List[Tuple[int, int, int]] = field(default_factory=list)   # (ANSYS elem, group, elem)
    created_nodes: List[int] = field(default_factory=list)            # K nodes created for beams
    counts: Dict[str, int] = field(default_factory=dict)

    @property
    def warnings(self) -> List[str]:
        return [t for k, t in self.messages if k == "warning"]

    @property
    def infos(self) -> List[str]:
        return [t for k, t in self.messages if k == "info"]

    def map_text(self, source: str = "") -> str:
        """Text of the ``<model>_cdb.map`` file (D-ANS-02): ``ansys_elem group elem`` per line."""
        lines = [f"! SASSI-EDU CONVERT,ANSYS element map{(': ' + source) if source else ''}",
                 "! ansys_elem  group  elem   (BEAM189 elements give two SASSI elements)"]
        lines += [f"{a:>12d} {g:>6d} {e:>6d}" for a, g, e in self.element_map]
        return "\n".join(lines) + "\n"


class _Converter:
    """Maps a :class:`CdbModel` onto an :class:`SSIModel` (D-ANS-01..05)."""

    def __init__(self, cdb: CdbModel, gravity: float, damping: Optional[float]):
        self.c = cdb
        self.g = float(gravity)
        self.beta0 = damping
        self.m = SSIModel()
        self.res = CdbConversion(self.m)
        self.warned: Dict[str, int] = {}
        self.next_node = (max(cdb.nodes) if cdb.nodes else 0) + 1
        self.sec_index: Dict[Tuple[float, ...], int] = {}
        self.sc_index: Dict[Tuple[float, ...], int] = {}
        self.mx_index: Dict[Tuple[Any, ...], int] = {}
        self.used_mats: Dict[int, List[str]] = {}      # mat -> element families using it
        self.groups: Dict[int, Group] = {}
        self.group_of_type: Dict[Tuple[int, int], int] = {}   # (itype, sassi code) -> group id
        self.next_group = (max(cdb.etypes) if cdb.etypes else 0) + 1
        self.link_nodes: set = set()                    # nodes of LINK-derived beams
        self.link_groups: set = set()                   # SASSI groups of LINK-derived beams
        self.k_cache: Optional[int] = None
        self.masses: Dict[int, np.ndarray] = {}
        self.solid_incomp: List[int] = []               # per SOLID/PLANE type: 0 include, 1 suppress
        # 2-D models: ANSYS X-Y plane -> SASSI X-Z plane, x_S = Q x_A (proper rotation about X)
        self.two_d = any(cdb.etypes.get(e.type) in (42, 182) for e in cdb.elements)
        self.Q = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]) if self.two_d else np.eye(3)
        self.Qa = np.abs(self.Q)
        self.dofmap = [int(np.argmax(self.Qa[:, a])) for a in range(3)]   # ANSYS axis a -> SASSI axis
        self.xyz: Dict[int, np.ndarray] = {}

    # ------------------------------------------------------------------ messages
    def warn(self, text: str) -> None:
        self.res.messages.append(("warning", text))

    def info(self, text: str) -> None:
        self.res.messages.append(("info", text))

    def warn_once(self, key: str, text: str) -> None:
        if key not in self.warned:
            self.warn(text)
        self.warned[key] = self.warned.get(key, 0) + 1

    # ------------------------------------------------------------------ helpers
    def _check_keyopts(self, itype: int, en: int) -> None:
        """Report the non-zero KEYOPTs of a supported type that the converter neither interprets nor
        knows to be without effect (spec 04 11.1); the interpreted ones get their own messages."""
        handled, no_effect = _KEYOPTS.get(en, (frozenset(), frozenset()))
        other = [(k, v) for k, v in sorted(self.c.keyopts.get(itype, {}).items())
                 if v != 0 and k not in handled and k not in no_effect]
        if other:
            self.warn(f"ET {itype} {element_name(en)}: " + ", ".join(f"KEYOPT({k}) = {v}" for k, v in other)
                      + " not converted (ignored; check that the option does not change the stiffness or mass)")

    def sassi_dof(self, a: int) -> int:
        """SASSI DOF index (0..5) of ANSYS DOF index a (UX..ROTZ) under the 2-D rotation."""
        return self.dofmap[a % 3] + (3 if a >= 3 else 0)

    def P(self, n: int) -> np.ndarray:
        return self.xyz[n]

    def group(self, itype: int, code: int, suffix: str = "") -> Group:
        key = (itype, code)
        gid = self.group_of_type.get(key)
        if gid is None:
            if not any(k[0] == itype for k in self.group_of_type):
                gid = itype
            else:
                gid = self.next_group
                self.next_group += 1
            self.group_of_type[key] = gid
            name = element_name(self.c.ename(itype))
            g = Group(gid, code, f"ANSYS ET {itype} {name}{suffix}")
            self.groups[gid] = g
        return self.groups[gid]

    def add(self, el: CdbElement, gkey: Tuple[int, int, str], nodes: List[int], **attrs) -> None:
        """Add a SASSI element to the group of ``gkey = (itype, code, title suffix)`` (created on
        first use, so element types whose elements are all rejected leave no empty group)."""
        g = self.group(*gkey)
        eid = len(g.elements) + 1
        g.elements[eid] = Element(eid, list(nodes), **attrs)
        self.res.element_map.append((el.num, g.id, eid))

    def new_node(self, xyz: np.ndarray) -> int:
        nid = self.next_node
        self.next_node += 1
        self.xyz[nid] = np.asarray(xyz, float)
        self.res.created_nodes.append(nid)
        return nid

    def section_id(self, sec: Dict[str, float]) -> int:
        key = tuple(float(sec[k]) for k in BeamSection.FIELDS)
        rid = self.sec_index.get(key)
        if rid is None:
            rid = len(self.sec_index) + 1
            self.sec_index[key] = rid
            self.m.sections[rid] = BeamSection(rid, *key)
        return rid

    def spring_id(self, k6: Sequence[float], damp: float = 0.0) -> int:
        key = tuple(float(v) for v in k6) + (float(damp),)
        sid = self.sc_index.get(key)
        if sid is None:
            sid = len(self.sc_index) + 1
            self.sc_index[key] = sid
            self.m.springs[sid] = SpringProp(sid, *key)
        return sid

    def matrix_id(self, key: Tuple[Any, ...], kind: str, A: np.ndarray) -> int:
        pid = self.mx_index.get(key)
        if pid is None:
            pid = len(self.mx_index) + 1
            self.mx_index[key] = pid
            p = MatrixProp(pid)
            for r in range(1, 13):
                terms = [float(v) + 0.0 for v in A[r - 1, r - 1:]]       # + 0.0: no -0.0 (exact WRITE/INP)
                if any(t != 0.0 for t in terms):
                    p.set_row(kind, r, terms)
            self.m.matrices[pid] = p
        return pid

    def use_mat(self, mat: int, family: str) -> None:
        self.used_mats.setdefault(mat, [])
        if family not in self.used_mats[mat]:
            self.used_mats[mat].append(family)

    def real(self, el: CdbElement) -> List[float]:
        return list(self.c.reals.get(el.real, []))

    def check_nodes(self, el: CdbElement, nodes: Sequence[int], name: str) -> bool:
        miss = [n for n in nodes if n and n not in self.xyz]
        if miss:
            self.warn_once(f"missing-{el.num}", f"{name} element {el.num}: nodes {miss} are not defined; "
                                                f"element not converted")
            return False
        return True

    # ------------------------------------------------------------------ nodes
    def nodes(self) -> None:
        for nid in sorted(self.c.nodes):
            self.xyz[nid] = self.Q @ self.c.nodes[nid] if self.two_d else np.array(self.c.nodes[nid], float)
        if self.two_d:
            zs = [abs(float(v[2])) for v in self.c.nodes.values()]
            if zs and max(zs) > 0.0:
                self.warn("PLANE elements found but some nodes have Z != 0: ANSYS 2-D models lie in the "
                          "X-Y plane; the model was rotated into the SASSI X-Z plane anyway")
            self.info("2-D model (PLANE elements): ANSYS (X, Y) -> SASSI (X, Z); y = -Z_ANSYS; DOF labels "
                      "UY <-> UZ and ROTY <-> ROTZ; HOUSE <dim> set to 1 (2-D)")
        if self.c.node_rot:
            self.warn(f"{len(self.c.node_rot)} nodes have nodal coordinate rotations (THXY, THYZ, THZX); "
                      "SASSI-EDU DOFs are in global axes: constraints, COMBIN14 KEYOPT(2) springs, MASS21 "
                      "and MATRIX27 data at these nodes are interpreted in global axes "
                      f"(nodes {_short(sorted(self.c.node_rot))})")

    # ------------------------------------------------------------------ element families
    def solid(self, itype: int, els: List[CdbElement]) -> None:
        en = self.c.ename(itype)
        name = element_name(en)
        kop = self.c.keyopts.get(itype, {})
        if en in (185, 186) and kop.get(3, 0) != 0:
            self.warn(f"ET {itype} {name}: KEYOPT(3) = {kop[3]} (layered solid) is not supported; "
                      "converted as a homogeneous solid (D-ANS-07)")
        if en in (185, 186, 187) and kop.get(6, 0) != 0:
            self.warn(f"ET {itype} {name}: KEYOPT(6) = {kop[6]} (mixed u-P formulation) is not available; the "
                      "SASSI SOLID is a pure displacement element (stiffer for nearly incompressible material)")
        if en == 185:
            k2 = kop.get(2, 0)
            self.solid_incomp.append(0 if k2 in (2, 3) else 1)
            if k2 == 1:
                self.warn(f"ET {itype} SOLID185: uniform reduced integration (KEYOPT(2) = 1) is not available; "
                          "SASSI SOLID uses full integration")
        elif en in (45, 65):
            self.solid_incomp.append(1 if kop.get(1, 0) == 1 else 0)
            if en == 45 and kop.get(2, 0) != 0:
                self.warn(f"ET {itype} SOLID45: uniform reduced integration (KEYOPT(2) = {kop[2]}) is not "
                          "available; SASSI SOLID uses full integration")
        if en == 65:
            self.warn(f"ET {itype} SOLID65: concrete cracking/crushing and smeared reinforcement are not "
                      "converted (linear SOLID)")
        if en in (186, 95, 187):
            self.warn(f"ET {itype} {name}: quadratic element converted to its corner nodes "
                      "(midside nodes are left unused)")
        g = (itype, 1, "")
        for el in els:
            n = list(el.nodes)
            if en == 187:
                if len(n) < 4:
                    self.warn_once(f"short-{itype}", f"{name} elements with fewer than 4 nodes skipped")
                    continue
                corners = [n[0], n[1], n[2], n[2], n[3], n[3], n[3], n[3]]
            else:
                if len(n) < 8:
                    self.warn_once(f"short-{itype}", f"{name} elements with fewer than 8 nodes skipped")
                    continue
                corners = n[:8]
            if not self.check_nodes(el, corners, name):
                continue
            self.add(el, g, corners, mat=el.mat)
            self.use_mat(el.mat, "SOLID")
        self._count_group("SOLID", g)

    def plane(self, itype: int, els: List[CdbElement]) -> None:
        en = self.c.ename(itype)
        name = element_name(en)
        kop = self.c.keyopts.get(itype, {})
        k3 = kop.get(3, 0)
        if k3 != 2:
            what = {0: "plane stress", 1: "axisymmetric", 3: "plane stress with thickness",
                    5: "generalized plane strain"}.get(k3, f"KEYOPT(3) = {k3}")
            self.warn(f"ET {itype} {name}: {what} is not available; SASSI PLANE is plane strain with unit "
                      "thickness")
        if en == 42:
            self.solid_incomp.append(1 if kop.get(2, 0) == 1 else 0)
        else:
            self.solid_incomp.append(0 if kop.get(1, 0) in (2, 3) else 1)
            if kop.get(1, 0) == 1:
                self.warn(f"ET {itype} PLANE182: uniform reduced integration (KEYOPT(1) = 1) is not available; "
                          "SASSI PLANE uses full integration")
            if kop.get(6, 0) != 0:
                self.warn(f"ET {itype} PLANE182: KEYOPT(6) = {kop[6]} (mixed u-P formulation) is not available; "
                          "the SASSI PLANE is a pure displacement element")
        g = (itype, 4, "")
        for el in els:
            n = list(el.nodes[:4])
            if len(n) < 3:
                continue
            if len(n) == 4 and n[2] == n[3]:
                n = n[:3]
            if not self.check_nodes(el, n, name):
                continue
            self.add(el, g, n, mat=el.mat)
            self.use_mat(el.mat, "PLANE")
        self._count_group("PLANE", g)

    def shell(self, itype: int, els: List[CdbElement]) -> None:
        en = self.c.ename(itype)
        name = element_name(en)
        kop = self.c.keyopts.get(itype, {})
        if kop.get(1, 0) != 0:
            what = "membrane-only or bending-only stiffness" if en == 63 else "membrane-only stiffness"
            self.warn(f"ET {itype} {name}: KEYOPT(1) = {kop[1]} ({what}) is not available; the SASSI SHELL "
                      "has membrane and bending stiffness (full shell stiffness used)")
        if en in (181, 281):
            self.info(f"ET {itype} {name}: converted to the thin (Kirchhoff) SASSI SHELL; the Mindlin "
                      f"transverse-shear flexibility of {name} is not carried over")
        if en == 281:
            self.warn(f"ET {itype} SHELL281: 8-node shell converted to its 4 corner nodes "
                      "(midside nodes are left unused)")
        g = (itype, 3, "")
        for el in els:
            n = list(el.nodes[:4])
            if len(n) < 3:
                continue
            thick, mat, added = self._shell_thickness(itype, en, el)
            if thick is None:
                continue
            if len(n) == 4 and n[2] == n[3]:
                n = n[:3]                         # triangle: SASSI wants 3 nodes (spec 04 11.3)
            if not self.check_nodes(el, n, name):
                continue
            self.add(el, g, n, mat=mat, thick=thick)
            self.use_mat(mat, "SHELL")
            if added > 0.0:
                # added mass per unit area, lumped like the SASSI shell mass (rho t A / n per corner,
                # D-ELM-07): ADMSUA A / n on the translations of every corner
                xyz = np.array([self.P(x) for x in n])
                if len(n) == 4:
                    area = 0.5 * float(np.linalg.norm(np.cross(xyz[2] - xyz[0], xyz[3] - xyz[1])))
                else:
                    area = 0.5 * float(np.linalg.norm(np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0])))
                for x in n:
                    self._add_mass(x, [added * area / len(n)] * 3 + [0.0] * 3, sassi=True)
        self._count_group("SHELL", g)

    def _shell_thickness(self, itype: int, en: int, el: CdbElement) -> Tuple[Optional[float], int, float]:
        """``(thickness, material, added mass per unit area)`` of a shell element (None = skip).

        SHELL63 real constants: TK(I..L) R1-R4, EFS R5 (elastic foundation), THETA R6 (element axes;
        no effect for isotropic material), RMI R7 (bending inertia ratio 12 I/TK^3, 0 = 1), CTOP/CBOT
        R8-R9 (stress output), ADMSUA R18 (added mass per unit area).  SHELL181/281 sections:
        SECDATA layers and SECCONTROL E11, E22, E12 (transverse shear stiffness), ADDMAS (added mass
        per unit area), hourglass and drilling factors.
        """
        name = element_name(en)
        if en == 63:
            r = self.real(el)
            if not r or r[0] <= 0.0:
                self.warn_once(f"thk-{itype}-{el.real}", f"{name} real set {el.real}: thickness TK(I) missing; "
                                                         "elements not converted")
                return None, el.mat, 0.0
            rr = r + [0.0] * 18
            key = f"{itype}-{el.real}"
            if rr[4] != 0.0:
                self.warn_once(f"efs-{key}", f"{name} real set {el.real}: elastic foundation stiffness EFS (R5) = "
                                             f"{_g(rr[4])} is not converted (model the foundation with springs)")
            if rr[6] not in (0.0, 1.0):
                self.warn_once(f"rmi-{key}", f"{name} real set {el.real}: bending moment of inertia ratio RMI "
                                             f"(R7) = {_g(rr[6])} is not converted (SASSI uses I = TK^3/12)")
            extra = [k + 1 for k in list(range(9, 17)) + list(range(18, len(r))) if rr[k] != 0.0]
            if extra:
                self.warn_once(f"s63x-{key}", f"{name} real set {el.real}: real constants {_rlist(extra)} are not "
                                              "converted")
            added = rr[17]
            if added != 0.0:
                self.warn_once(f"admsua-{key}", f"{name} real set {el.real}: added mass per unit area ADMSUA (R18) "
                                                f"= {_g(added)} lumped as ADMSUA*A/n at the corner nodes (MT in mass "
                                                "units, like the lumped SASSI shell mass)")
            tk = [v for v in r[:4] if v > 0.0]
            if len(set(tk)) > 1:
                self.warn_once(f"taper-{itype}-{el.real}", f"{name} real set {el.real}: tapered thickness "
                                                           f"{tk}; the average is used")
                return float(sum(tk) / len(tk)), el.mat, added
            return float(r[0]), el.mat, added
        sec = self.c.sections.get(el.secnum)
        if sec is None or sec.type != "SHELL" or not sec.layers:
            r = self.real(el)
            if r and r[0] > 0.0:
                self.warn_once(f"thk-{itype}-{el.real}", f"{name}: no SHELL section; thickness taken from "
                                                         f"real set {el.real}")
                extra = [k + 1 for k in range(1, len(r)) if r[k] != 0.0]
                if extra:
                    self.warn_once(f"shr-{itype}-{el.real}", f"{name} real set {el.real}: real constants "
                                                             f"{_rlist(extra)} are not converted")
                return float(r[0]), el.mat, 0.0
            self.warn_once(f"sec-{itype}-{el.secnum}", f"{name} element {el.num}: SHELL section {el.secnum} "
                                                       "missing; elements not converted")
            return None, el.mat, 0.0
        t = float(sum(l[0] for l in sec.layers))
        mats = [l[1] for l in sec.layers]
        if len(sec.layers) > 1:
            self.warn_once(f"layers-{sec.id}", f"SHELL section {sec.id}: {len(sec.layers)} layers "
                                               f"(materials {mats}) converted to one homogeneous layer of "
                                               f"thickness {_g(t)} and material {mats[0]}")
        off = sec.offset[0] if sec.offset else "MID"
        if off not in ("MID", ""):
            self.warn_once(f"shoff-{sec.id}", f"SHELL section {sec.id}: offset {off} is not represented "
                                              "(SASSI shells lie on their mid-surface)")
        ctl = list(sec.control) + [0.0] * 4
        if any(v != 0.0 for v in ctl[:3]):
            self.warn_once(f"shts-{sec.id}", f"SHELL section {sec.id}: SECCONTROL transverse shear stiffness "
                                             f"E11, E22, E12 = {_g(ctl[0])}, {_g(ctl[1])}, {_g(ctl[2])} not used (the "
                                             "Kirchhoff SASSI SHELL has no transverse shear deformation)")
        added = ctl[3]
        if added != 0.0:
            self.warn_once(f"shadd-{sec.id}", f"SHELL section {sec.id}: SECCONTROL added mass per unit area "
                                              f"ADDMAS = {_g(added)} lumped as ADDMAS*A/n at the corner nodes (MT "
                                              "in mass units, like the lumped SASSI shell mass)")
        mat = mats[0] if mats[0] > 0 else el.mat
        return t, mat, added

    # ---------------------------------------------------------------- beams
    def _beam_props(self, itype: int, en: int, el: CdbElement) -> Optional[BeamProps]:
        name = element_name(en)
        if en in (4, 44):
            r = self.real(el) + [0.0] * 30
            if el.real not in self.c.reals or r[0] <= 0.0:
                self.warn_once(f"r-{itype}-{el.real}", f"{name} real set {el.real} missing or AREA <= 0; "
                                                       "elements not converted")
                return None
            if en == 4:
                A, Izz, Iyy, ixx, shz, shy, addm = r[0], r[1], r[2], r[7], r[8], r[9], r[11]
                for k, what in ((6, "ISTRN (initial strain, a load)"), (10, "SPIN (rotational frequency, "
                                                                            "gyroscopic effects)")):
                    if r[k] != 0.0:
                        self.warn_once(f"b4r{k}-{itype}-{el.real}", f"{name} real set {el.real}: R{k + 1} {what} "
                                                                    f"= {_g(r[k])} is not converted")
                extra = [k + 1 for k in range(12, len(r)) if r[k] != 0.0]
                if extra:
                    self.warn_once(f"b4x-{itype}-{el.real}", f"{name} real set {el.real}: real constants "
                                                             f"{_rlist(extra)} beyond ADDMAS (R12) are not "
                                                             "converted")
            else:
                A1, Iz1, Iy1, Ix1 = r[0], r[1], r[2], r[5]
                A2, Iz2, Iy2, Ix2 = (r[6] or A1), (r[7] or Iz1), (r[8] or Iy1), (r[11] or Ix1)
                if (A2, Iz2, Iy2, Ix2) != (A1, Iz1, Iy1, Ix1):
                    self.warn_once(f"taper-{itype}-{el.real}", f"{name} real set {el.real}: tapered section; "
                                                               "SASSI beams are prismatic, the average of the "
                                                               "node I and J properties is used")
                A, Izz, Iyy, ixx = (A1 + A2) / 2, (Iz1 + Iz2) / 2, (Iy1 + Iy2) / 2, (Ix1 + Ix2) / 2
                shz, shy, addm = r[18], r[19], 0.0
                if any(v != 0.0 for v in r[12:18]):
                    self.warn_once(f"off-{itype}-{el.real}", f"{name} real set {el.real}: end offsets DX..DZ "
                                                             "are not represented")
                extra = [k + 1 for k in range(24, len(r)) if r[k] != 0.0]
                if extra:
                    # spec 04 11.3: the BEAM44 real-constant block is read up to R24 (TKYT2); later fields
                    # (shear-centre offsets, elastic foundation, THETA, ISTRN, ADDMAS, ...) are not
                    self.warn_once(f"b44x-{itype}-{el.real}", f"{name} real set {el.real}: real constants "
                                                              f"{_rlist(extra)} (beyond R24) are not converted")
            J = ixx if ixx > 0.0 else Iyy + Izz
            return BeamProps(A=A, Iyy=Iyy, Izz=Izz, J=J, Asy=A / shy if shy > 0 else 0.0,
                             Asz=A / shz if shz > 0 else 0.0, addmas=addm)
        sec = self.c.sections.get(el.secnum)
        if sec is None:
            self.warn_once(f"sec-{itype}-{el.secnum}", f"{name} element {el.num}: section {el.secnum} is not "
                                                       "defined; elements not converted")
            return None
        st = sec.subtype
        d = list(sec.data) + [0.0] * 12
        if sec.block:
            self.warn_once(f"secb-{sec.id}", f"section {sec.id}: mesh-based (SECBLOCK) beam sections are not "
                                             "supported; elements not converted")
            return None
        if sec.type == "BEAM" and st == "RECT":
            if d[0] <= 0 or d[1] <= 0:
                self.warn_once(f"rect-{sec.id}", f"RECT section {sec.id}: B and H must be > 0")
                return None
            p = rect_props(d[0], d[1])
        elif sec.type == "BEAM" and st == "ASEC":
            A, Iyy, Iyz, Izz, Iw, J = d[0], d[1], d[2], d[3], d[4], d[5]
            if A <= 0.0 or J <= 0.0:
                self.warn_once(f"asec-{sec.id}", f"ASEC section {sec.id}: A and J must be > 0; elements not "
                                                 "converted")
                return None
            p = BeamProps(A=A, Iyy=Iyy, Izz=Izz, J=J, Iyz=Iyz)
            if any(v != 0.0 for v in d[6:10]):
                self.warn_once(f"asecoff-{sec.id}", f"ASEC section {sec.id}: centroid / shear-centre offsets "
                                                    "CGy CGz SHy SHz are not represented")
            if Iw != 0.0:
                self.warn_once(f"asecw-{sec.id}", f"ASEC section {sec.id}: warping constant ignored")
        elif sec.type == "BEAM" and st == "CSOLID":
            p = tube_props(2.0 * d[0], 0.0)
        elif sec.type == "BEAM" and st == "CTUBE":
            p = tube_props(2.0 * d[1], 2.0 * d[0])
        elif sec.type == "PIPE":
            Do, tw = d[0], d[1]
            if Do <= 0 or tw <= 0 or 2 * tw > Do:
                self.warn_once(f"pipe-{sec.id}", f"PIPE section {sec.id}: invalid Do/tw")
                return None
            p = tube_props(Do, Do - 2.0 * tw)
            if f"pipe{sec.id}" not in self.warned:
                self.warned[f"pipe{sec.id}"] = 1
                self.info(f"PIPE section {sec.id}: pipe converted to an equivalent straight beam (internal "
                          "fluid, insulation and pressure effects are not converted)")
        else:
            self.warn_once(f"sect-{sec.id}", f"{name}: section {sec.id} type {sec.type} {st} is not supported "
                                             "(BEAM RECT, ASEC, CSOLID, CTUBE and PIPE are); elements not converted")
            return None
        if sec.type in ("BEAM", "PIPE"):
            off = sec.offset[0] if sec.offset else "CENT"
            if off not in ("CENT", "") and not (off == "ORIGIN" and st in ("RECT", "CSOLID", "CTUBE")) and \
                    not (off == "USER" and all(_num(v) == 0.0 for v in sec.offset[1:3] if _isnum(v))):
                self.warn_once(f"boff-{sec.id}", f"section {sec.id}: offset {off} is not represented (SASSI "
                                                 "beams are centroidal)")
            if sec.control and len(sec.control) >= 2 and (sec.control[0] > 0 or sec.control[1] > 0):
                G = self._shear_modulus(el.mat)
                if G > 0:
                    if sec.control[0] > 0:
                        p.Asz = sec.control[0] / G     # TXZ = k G A (shear along z)
                    if sec.control[1] > 0:
                        p.Asy = sec.control[1] / G
            if len(sec.control) >= 3 and sec.control[2] > 0:
                p.addmas = sec.control[2]
        return p

    def _shear_modulus(self, mat: int) -> float:
        mp = self.c.materials.get(mat, {})
        E = mp.get("EX", 0.0)
        nu = mp.get("PRXY", mp.get("NUXY", 0.3))
        return E / (2.0 * (1.0 + nu)) if E > 0 else 0.0

    def _k_node(self, I: int, J: int, e2: np.ndarray, L: float) -> int:
        """A K node for SASSI axis 2 along ``e2`` (reuses the last created K node when it fits)."""
        from ..elements.base import frame_from_three_points
        xi, xj = self.P(I), self.P(J)
        if self.k_cache is not None:
            try:
                lam, _ = frame_from_three_points(xi, xj, self.P(self.k_cache))
                if float(lam[1] @ e2) > 1.0 - 1e-12:
                    return self.k_cache
            except Exception:
                pass
        k = self.new_node(xi + L * e2)
        self.k_cache = k
        return k

    def beam(self, itype: int, els: List[CdbElement]) -> None:
        en = self.c.ename(itype)
        name = element_name(en)
        kop = self.c.keyopts.get(itype, {})
        link = en in (8, 180)
        if en in (188, 189, 288) and kop.get(1, 0) != 0:
            self.warn(f"ET {itype} {name}: warping DOF (KEYOPT(1) = {kop[1]}) is not available")
        if en in (188, 288) and kop.get(3, 0) != 3:
            self.info(f"ET {itype} {name}: KEYOPT(3) = {kop.get(3, 0)} ({'quadratic' if kop.get(3, 0) == 2 else 'linear'}"
                      " shape functions) -- each element becomes one SASSI Timoshenko beam, which is exact for end "
                      "loads, so coarse ANSYS meshes of this type give different results (use KEYOPT(3) = 3 or a "
                      "finer ANSYS mesh for cross-checks)")
        if en == 180 and kop.get(3, 0) != 0:
            self.warn(f"ET {itype} LINK180: KEYOPT(3) = {kop[3]} (tension-only or compression-only) is nonlinear; "
                      "converted as tension and compression")
        rel_i = rel_j = None
        if en == 44:
            rel_i, rel_j = _release_digits(kop.get(7, 0)), _release_digits(kop.get(8, 0))
        g = (itype, 2, " (axial only)" if link else "")
        for el in els:
            n = list(el.nodes)
            if en == 189:
                if len(n) < 3 or n[2] == 0:
                    self.warn_once(f"b189-{itype}", f"BEAM189 elements without a midside node skipped")
                    continue
                parts = [(n[0], n[2]), (n[2], n[1])]
                kn = n[3] if len(n) > 3 and n[3] else 0
            else:
                if len(n) < 2:
                    continue
                parts = [(n[0], n[1])]
                kn = n[2] if len(n) > 2 and n[2] and not link else 0
            if not self.check_nodes(el, [x for p in parts for x in p] + ([kn] if kn else []), name):
                continue
            if link:
                p = self._link_props(itype, en, el)
            else:
                p = self._beam_props(itype, en, el)
            if p is None:
                continue
            I0, J0 = (n[0], n[1])
            theta = self.real(el)[5] if en == 4 and len(self.real(el)) > 5 and not kn else 0.0
            # ANSYS frame (whole member for BEAM189: I..J), mapped to SASSI axes
            try:
                ax, _ = ansys_beam_axes(self.P(I0), self.P(J0), self.P(kn) if kn else None, theta)
            except ValueError as exc:
                self.warn(f"{name} element {el.num}: {exc}; element not converted")
                continue
            if link:
                sec = dict(axial=p.A, shear2=0.0, shear3=0.0, tors=1e-8 * p.A * p.A, flex2=0.0, flex3=0.0)
                p2 = np.array([0.0, 1.0])
            else:
                sec, p2 = sassi_section_from_ansys(p)
            e2 = p2[0] * ax[1] + p2[1] * ax[2]
            rotated = abs(p2[1]) < 1.0 - 1e-12
            if rotated:
                self.warn_once(f"princ-{itype}", f"{name}: sections with Iyz != 0 are converted in their "
                                                 "principal axes (SASSI axis 2 along the principal direction "
                                                 "nearest to ANSYS z; Iyz = int(y z dA) assumed)")
            ki = kj = [0] * 6
            if rel_i is not None:
                ki, kj = _sassi_release(rel_i), _sassi_release(rel_j)
                if rotated and (ki != [0] * 6 or kj != [0] * 6):
                    self.warn_once(f"relrot-{itype}", f"{name}: releases in rotated principal axes are "
                                                      "approximate")
                ki, kj, sec = self._same_release_both_ends(itype, kop, ki, kj, sec)
            rid = self.section_id(sec)
            for (a, b) in parts:
                L = float(np.linalg.norm(self.P(b) - self.P(a)))
                if L <= 0.0:
                    self.warn(f"{name} element {el.num}: zero length; not converted")
                    continue
                if kn and not rotated and not link:
                    knode = kn                       # ANSYS K: z towards K = SASSI axis 2
                else:
                    knode = self._k_node(a, b, e2, L)
                self.add(el, g, [a, b, knode], mat=el.mat, prop=rid, ki=list(ki), kj=list(kj))
                if link:
                    self.link_nodes.update((a, b))
                    self.link_groups.add(self.group_of_type[(itype, 2)])
                if p.addmas > 0.0:
                    for x in (a, b):
                        self._add_mass(x, [p.addmas * L / 2.0] * 3 + [0.0] * 3, sassi=True)
            if p.addmas > 0.0:
                src = "SECCONTROL" if (en in (180, 188, 189, 288) and self.c.sections.get(el.secnum)) else \
                    f"real set {el.real}"
                self.warn_once(f"addmas-{itype}-{src}",
                               f"ET {itype} {name} ({src}): added mass per unit length ADDMAS = {_g(p.addmas)} lumped "
                               "as ADDMAS*L/2 at both end nodes (MT in mass units; spec 04 11.4 optional "
                               "extension: ANSYS distributes it along the element)")
            self.use_mat(el.mat, "LINK" if link else "BEAM")
        grp = self.groups.get(self.group_of_type.get((itype, 2), -1))
        if grp is not None and any(any(e.ki) or any(e.kj) for e in grp.elements.values()):
            self.info(f"ET {itype} {name}: end releases KEYOPT(7) = {kop.get(7, 0)}, KEYOPT(8) = {kop.get(8, 0)} "
                      "mapped to KI/KJ (P1 P2 P3 M1 M2 M3 = UX UZ UY ROTX ROTZ ROTY)")
        self._count_group("BEAMS", g)

    def _same_release_both_ends(self, itype: int, kop: Dict[int, int], ki: List[int], kj: List[int],
                                sec: Dict[str, float]) -> Tuple[List[int], List[int], Dict[str, float]]:
        """BEAM44 with the same component released at I and J (spec 08 5.30: CHECK Error 10 in SASSI).

        ANSYS accepts it -- KEYOPT(7) = KEYOPT(8) = 11 is the usual pin-ended brace -- so the
        converter maps it onto an equivalent SASSI member:

        * M2 or M3 released at both ends: no end moment in that bending plane, hence no shear
          (V = (M_I + M_J)/L) and no bending stiffness at all -- the member only follows the rigid
          rotation of the plane.  SASSI represents exactly that with a zero inertia of the plane
          (M3 -> I3 = 0, M2 -> I2 = 0; beam.py: I = 0 is accepted and makes the plane's releases
          redundant), so the releases are dropped.
        * P1, P2, P3 or M1 released at both ends: one release already removes that force from the
          whole member (zero axial force, zero shear -> constant moment, zero torque); a second one
          only adds a rigid-body mechanism (singular in ANSYS too).  The release is kept at I only:
          same stiffness.
        """
        both = [c for c in range(6) if ki[c] and kj[c]]
        if not both:
            return ki, kj, sec
        ki, kj, sec = list(ki), list(kj), dict(sec)
        names = ["P1", "P2", "P3", "M1", "M2", "M3"]
        planes = []
        for c in both:
            if c in (4, 5):
                sec["flex2" if c == 4 else "flex3"] = 0.0
                ki[c] = kj[c] = 0
                planes.append(f"{names[c]} -> I{'2' if c == 4 else '3'} = 0")
            else:
                kj[c] = 0
                planes.append(f"{names[c]} at I only")
        why = []
        if any(c in (4, 5) for c in both):
            why.append("A moment released at both ends leaves its bending plane without stiffness (I = 0 "
                       "represents it exactly); the consistent mass of that plane then follows the SASSI beam "
                       "shape and the torsional mass rho (I2 + I3) loses the zeroed inertia.  As in ANSYS, a node "
                       "joined only by such members has no stiffness in those rotations: restrain them (FIXROT "
                       "or D).")
        if any(c < 4 for c in both):
            why.append("A force or torque released at one end is already zero along the whole member; the "
                       "second release only adds a mechanism, so the member keeps the same stiffness.")
        self.warn_once(f"relboth-{itype}-{''.join(map(str, both))}",
                       f"ET {itype} BEAM44: KEYOPT(7) = {kop.get(7, 0)}, KEYOPT(8) = {kop.get(8, 0)} release "
                       f"{', '.join(names[c] for c in both)} at both ends, which SASSI does not allow (CHECK "
                       f"Error 10); mapped to an equivalent member: {'; '.join(planes)}.  " + "  ".join(why))
        return ki, kj, sec

    def _link_props(self, itype: int, en: int, el: CdbElement) -> Optional[BeamProps]:
        """Area and added mass of LINK8 (R: AREA, ISTRN) and LINK180 (SECTYPE,LINK + SECDATA AREA and
        SECCONTROL ADDMAS, or the older R: AREA, ADDMAS, TENSKEY).  Tension- or compression-only
        behaviour (LINK180 KEYOPT(3) or TENSKEY) is nonlinear and is reported."""
        name = element_name(en)
        A = addm = 0.0
        sec = self.c.sections.get(el.secnum)
        r = self.real(el) + [0.0] * 3
        if en == 180 and sec is not None and sec.type == "LINK" and sec.data:
            A = sec.data[0]
            addm = sec.control[0] if sec.control else 0.0
        else:
            A = r[0]
            if en == 180:
                addm = r[1]
                if r[2] != 0.0:
                    self.warn_once(f"tens-{itype}-{el.real}", f"{name} real set {el.real}: TENSKEY = {_g(r[2])} "
                                                              "(tension or compression only) is nonlinear; "
                                                              "converted as tension and compression")
            elif r[1] != 0.0:
                self.warn_once(f"istrn-{itype}-{el.real}", f"{name} real set {el.real}: ISTRN (initial strain, "
                                                           f"a load) = {_g(r[1])} is not converted")
        if A <= 0.0:
            self.warn_once(f"link-{itype}", f"{name}: cross-sectional area missing; elements not converted")
            return None
        return BeamProps(A=A, Iyy=0.0, Izz=0.0, J=0.0, addmas=max(addm, 0.0), axial_only=True)

    # ---------------------------------------------------------------- springs, masses, matrices
    def spring(self, itype: int, els: List[CdbElement]) -> None:
        kop = self.c.keyopts.get(itype, {})
        k2, k3 = kop.get(2, 0), kop.get(3, 0)
        if k2 in (7, 8) or k2 > 8 or k3 not in (0, 1, 2):
            self.warn(f"ET {itype} COMBIN14: KEYOPT(2) = {k2}, KEYOPT(3) = {k3} is not supported (KEYOPT(2) = "
                      "1..6 or KEYOPT(3) = 0..2); elements not converted")
            return
        if k2 == 0 and k3 == 0:
            self.info(f"ET {itype} COMBIN14: 3-D longitudinal springs (KEYOPT(2) = KEYOPT(3) = 0) converted as "
                      "axial springs along I-J (SASSI-EDU extension of the manual's KEYOPT list)")
        g_sp = (itype, 7, "")
        g_gm = (itype, 9, " (oblique springs)")
        for el in els:
            n = list(el.nodes[:2])
            if len(n) < 2 or not self.check_nodes(el, n, "COMBIN14"):
                continue
            r = self.real(el) + [0.0] * 3
            k = r[0]
            if r[1] != 0.0 or r[2] != 0.0:
                self.warn_once(f"cv-{itype}", "COMBIN14: viscous damping CV1/CV2 has no SASSI counterpart "
                                              "(SC damping is a hysteretic ratio); ignored")
            if k == 0.0:
                self.warn_once(f"k0-{itype}", f"COMBIN14 real set {el.real}: K = 0; elements not converted")
                continue
            k6 = [0.0] * 6
            if k2 in range(1, 7):
                k6[self.sassi_dof(k2 - 1)] = k
                self.add(el, g_sp, n, prop=self.spring_id(k6))
                continue
            d = self.P(n[1]) - self.P(n[0])
            L = float(np.linalg.norm(d))
            if L <= 0.0:
                self.warn_once(f"zl-{itype}", "COMBIN14: zero-length axial springs need KEYOPT(2) (the direction "
                                              "is undefined); elements not converted")
                continue
            e = d / L
            off = 3 if k3 == 1 else 0                   # torsional spring: rotations
            a = int(np.argmax(np.abs(e)))
            if abs(e[a]) >= 1.0 - 1e-9:
                k6[off + a] = k
                self.add(el, g_sp, n, prop=self.spring_id(k6))
            else:
                # D-ANS-04: oblique axial spring -> 2-node GENERAL element K = k [ee^T -ee^T; -ee^T ee^T]
                kb = k * np.outer(e, e)
                K = np.zeros((12, 12))
                for (i0, j0, s) in ((0, 0, 1.0), (6, 6, 1.0), (0, 6, -1.0), (6, 0, -1.0)):
                    K[i0 + off:i0 + off + 3, j0 + off:j0 + off + 3] = s * kb
                pid = self.matrix_id(("C14", el.real, itype, tuple(np.round(e, 15))), "R", K)
                self.add(el, g_gm, n, prop=pid)
                self.warn_once(f"obl-{itype}", "COMBIN14: axial springs not parallel to a global axis converted "
                                               "to 2-node GENERAL elements k e e^T (D-ANS-04)")
        self._count_group("SPRING", g_sp)
        self._count_group("GENERAL", g_gm)

    def _add_mass(self, node: int, m6: Sequence[float], sassi: bool = False) -> None:
        """Accumulate a nodal mass (6 values in ANSYS DOF order unless ``sassi``)."""
        v = np.zeros(6)
        if sassi:
            v[:] = m6
        else:
            for a in range(6):
                v[self.sassi_dof(a)] += m6[a]
        self.masses[node] = self.masses.get(node, np.zeros(6)) + v

    def mass(self, itype: int, els: List[CdbElement]) -> None:
        kop = self.c.keyopts.get(itype, {})
        k1, k2, k3 = kop.get(1, 0), kop.get(2, 0), kop.get(3, 0)
        if k2 != 0:
            self.warn(f"ET {itype} MASS21: KEYOPT(2) = {k2} (nodal coordinate system); masses applied in "
                      "global axes")
        if k3 not in (0, 2, 3, 4):
            self.warn(f"ET {itype} MASS21: KEYOPT(3) = {k3} is not supported; elements not converted")
            return
        n_done = 0
        for el in els:
            if not el.nodes or not self.check_nodes(el, el.nodes[:1], "MASS21"):
                continue
            r = self.real(el) + [0.0] * 6
            if k3 == 0:
                m6 = r[:6]
            elif k3 == 2:
                m6 = [r[0]] * 3 + [0.0] * 3
            elif k3 == 3:
                m6 = [r[0], r[0], 0.0, 0.0, 0.0, r[1]]
            else:
                m6 = [r[0], r[0], 0.0, 0.0, 0.0, 0.0]
            if k1 == 1:
                dens = self.c.materials.get(el.mat, {}).get("DENS", 0.0)
                if dens <= 0.0:
                    self.warn_once(f"m21d-{itype}", f"MASS21 KEYOPT(1) = 1 needs DENS of material {el.mat}; "
                                                    "elements not converted")
                    continue
                m6 = [v * dens for v in m6]
            self._add_mass(el.nodes[0], m6)
            n_done += 1
        if k1 == 1 and n_done:
            self.info(f"ET {itype} MASS21: KEYOPT(1) = 1 (volume x density) converted to masses with DENS")
        self.res.counts["MASS21"] = self.res.counts.get("MASS21", 0) + n_done

    def matrix(self, itype: int, els: List[CdbElement]) -> None:
        kop = self.c.keyopts.get(itype, {})
        k2, k3 = kop.get(2, 0), kop.get(3, 0)
        kind = {4: "R", 2: "M"}.get(k3)
        if kind is None:
            what = "damping matrices" if k3 == 5 else f"KEYOPT(3) = {k3}"
            self.warn(f"ET {itype} MATRIX27: {what} are not supported (stiffness KEYOPT(3) = 4 and mass "
                      "KEYOPT(3) = 2 are); elements not converted")
            return
        if k2 != 0:
            self.warn(f"ET {itype} MATRIX27: unsymmetric matrices (KEYOPT(2) = 1) are not supported; elements "
                      "not converted")
            return
        g = (itype, 9, "")
        T = np.zeros((12, 12))
        for b in range(4):
            T[3 * b:3 * b + 3, 3 * b:3 * b + 3] = self.Q
        for el in els:
            n = list(el.nodes[:2])
            if len(n) < 2 or n[0] == n[1] or not self.check_nodes(el, n, "MATRIX27"):
                self.warn_once(f"m27n-{itype}", "MATRIX27 elements need two distinct nodes I and J; elements "
                                                "not converted")
                continue
            r = self.real(el) + [0.0] * 78
            A = np.zeros((12, 12))
            pos = 0
            for i in range(12):                    # C1..C78: upper triangle row by row
                A[i, i:] = r[pos:pos + 12 - i]
                pos += 12 - i
            A = np.triu(A) + np.triu(A, 1).T
            if self.two_d:
                A = T @ A @ T.T
            pid = self.matrix_id(("M27", itype, el.real), kind, A)
            self.add(el, g, n, prop=pid)
        if self.group_of_type.get((itype, 9)) is not None:
            self.info(f"ET {itype} MATRIX27 ({'stiffness' if kind == 'R' else 'mass'}): 2-node GENERAL elements "
                      "in global axes" + (" (mass in mass units, MOPT <matrix> = 0)" if kind == "M" else ""))
        self._count_group("GENERAL", g)

    def _count(self, key: str, n: int) -> None:
        self.res.counts[key] = self.res.counts.get(key, 0) + n

    def _count_group(self, key: str, gkey: Tuple[int, int, str]) -> None:
        gid = self.group_of_type.get(gkey[:2])
        if gid is not None:
            self._count(key, len(self.groups[gid].elements))

    # ------------------------------------------------------------------ materials
    def materials(self) -> None:
        missing_damp: List[int] = []
        for mat in sorted(self.used_mats):
            fam = self.used_mats[mat]
            mp = self.c.materials.get(mat)
            if mp is None:
                self.warn(f"material {mat} (used by {', '.join(fam)} elements) has no MPDATA; M,{mat} not defined")
                continue
            E = mp.get("EX", 0.0)
            if E <= 0.0:
                self.warn(f"material {mat}: EX missing or <= 0; M,{mat} not defined")
                continue
            if "PRXY" in mp:
                nu = mp["PRXY"]
            elif "NUXY" in mp:
                nu = mp["NUXY"]
            elif "GXY" in mp and mp["GXY"] > 0:
                nu = E / (2.0 * mp["GXY"]) - 1.0
                self.info(f"material {mat}: Poisson's ratio {_g(nu)} from EX and GXY")
            else:
                nu = 0.3
                self.warn(f"material {mat}: no PRXY/NUXY; ANSYS default 0.3 used")
            self._check_material_labels(mat, mp, E, nu)
            dens = mp.get("DENS", 0.0)
            if dens <= 0.0:
                self.warn(f"material {mat}: DENS missing or 0; weight 0 (massless elements)")
            if "DMPR" in mp:
                beta = mp["DMPR"] / 2.0
                if mp.get("DMPS", 0.0) != 0.0:
                    self.warn(f"material {mat}: both DMPR and DMPS given; beta = DMPR/2 used, DMPS = "
                              f"{_g(mp['DMPS'])} ignored")
            elif "DMPS" in mp:
                beta = mp["DMPS"] / 2.0
                self.info(f"material {mat}: beta = DMPS/2 = {_g(beta)} (constant structural damping coefficient)")
            else:
                beta = self.beta0 if self.beta0 is not None else 0.0
                missing_damp.append(mat)
            for lab in ("DAMP", "BETD", "ALPD"):
                if mp.get(lab, 0.0) != 0.0:
                    self.warn(f"material {mat}: Rayleigh damping {lab} = {_g(mp[lab])} has no frequency-independent "
                              "SASSI counterpart; ignored")
            if beta >= 0.5:
                self.warn(f"material {mat}: damping ratio {_g(beta)} >= 0.5 (EDU-04); check DMPR")
            if mat in self.c.temp_dependent:
                self.warn(f"material {mat}: temperature-dependent data ({', '.join(self.c.temp_dependent[mat])}); "
                          "the first temperature is used")
            self.m.materials[mat] = Material(mat, float(E), float(nu), float(dens * self.g), float(beta),
                                             float(beta), 1)
        if missing_damp:
            b = self.beta0 if self.beta0 is not None else 0.0
            self.warn(f"materials {_short(missing_damp)}: no DMPR; damping ratio {_g(b)} used "
                      "(CONVERT argument <damp>; D-ANS-01 maps beta = DMPR/2)")

    def _check_material_labels(self, mat: int, mp: Dict[str, float], E: float, nu: float) -> None:
        """D-ANS-01: MPDATA other than EX, PRXY/NUXY, DENS, DMPR/DMPS is ignored *with a warning*.

        Orthotropic labels that only repeat the isotropic values (EY = EZ = EX, PRYZ = PRXZ = nu,
        G = EX/(2(1 + nu)), relative tolerance 1e-6) are accepted silently; any difference is
        reported, because SASSI-EDU uses the isotropic EX and nu.  Rayleigh damping labels have
        their own message (in :meth:`materials`)."""
        G = E / (2.0 * (1.0 + nu))
        iso = {"EY": E, "EZ": E, "PRXY": nu, "NUXY": nu, "PRYZ": nu, "PRXZ": nu, "NUYZ": nu, "NUXZ": nu,
               "GXY": G, "GYZ": G, "GXZ": G}
        used = {"EX", "DENS", "DMPR", "DMPS", "DAMP", "BETD", "ALPD"}
        aniso, other = [], []
        for lab in sorted(mp):
            v = mp[lab]
            if lab in used:
                continue
            if lab in iso:
                ref = iso[lab]
                if abs(v - ref) > 1e-6 * max(abs(ref), 1e-300):
                    aniso.append(f"{lab} = {_g(v)} (isotropic value {_g(ref)})")
            elif v != 0.0:
                other.append(lab)
        if aniso:
            self.warn(f"material {mat}: orthotropic or inconsistent elastic data not converted: {', '.join(aniso)}; "
                      f"the isotropic EX = {_g(E)} and nu = {_g(nu)} are used")
        if other:
            self.warn(f"material {mat}: MPDATA {', '.join(other)} not converted (SASSI-EDU materials carry EX, nu, "
                      "density and a damping ratio only; D-ANS-01)")

    # ------------------------------------------------------------------ boundary data
    def constraints(self) -> None:
        nfix = 0
        bad: Dict[str, int] = {}
        for node, lab, val, val2 in self.c.constraints:
            if node not in self.xyz:
                bad["undefined node"] = bad.get("undefined node", 0) + 1
                continue
            if val != 0.0 or val2 != 0.0:
                bad["non-zero value"] = bad.get("non-zero value", 0) + 1
                continue
            if lab == "ALL":
                idx = list(range(6))
            elif lab in DOF_LABELS:
                idx = [self.sassi_dof(DOF_LABELS.index(lab))]
            else:
                bad[f"label {lab}"] = bad.get(f"label {lab}", 0) + 1
                continue
            fx = self.fix.setdefault(node, [0] * 6)
            for i in idx:
                fx[i] = 1
            nfix += 1
        for why, n in bad.items():
            self.warn(f"{n} D constraints not converted ({why}; only zero displacements become SASSI fixities)")
        if nfix:
            self.res.counts["D"] = nfix

    def components(self) -> None:
        names = []
        for name, (kind, ids) in sorted(self.c.components.items()):
            if name == "SSI_INT" and kind.startswith("NODE"):
                ok = [n for n in ids if n in self.xyz]
                for n in ok:
                    self.m.nodes[n].flags.add(0)
                self.info(f"node component SSI_INT: {len(ok)} interaction nodes (INT code 0)")
            else:
                names.append(f"{name} ({kind.lower()}, {len(ids)})")
        if names:
            self.info("components not converted (only the node component SSI_INT becomes interaction nodes): "
                      + ", ".join(names))

    # ------------------------------------------------------------------ driver
    def run(self) -> CdbConversion:
        c = self.c
        self.fix: Dict[int, List[int]] = {}
        self.nodes()
        by_type: Dict[int, List[CdbElement]] = {}
        for el in c.elements:
            by_type.setdefault(el.type, []).append(el)
        family_fn = {"SOLID": self.solid, "PLANE": self.plane, "SHELL": self.shell, "BEAM": self.beam,
                     "LINK": self.beam, "SPRING": self.spring, "MASS": self.mass, "MATRIX": self.matrix}
        for itype in sorted(by_type):
            en = c.ename(itype)
            els = by_type[itype]
            if itype not in c.etypes:
                self.warn(f"{len(els)} elements of type {itype} without an ET command; not converted")
                continue
            fam = ANSYS_ELEMENTS.get(en)
            if fam is None:
                why = _KNOWN_UNSUPPORTED.get(en) or f"element {en}"
                self.warn(f"ET {itype} {why}: unsupported element type; {len(els)} elements not converted")
                self._count("unsupported elements", len(els))
                continue
            self._check_keyopts(itype, en)
            family_fn[fam[0]](itype, els)
        # nodes (created K nodes included) and fixities
        for nid in sorted(self.xyz):
            p = self.xyz[nid]
            self.m.define_node(nid, (float(p[0]) + 0.0, float(p[1]) + 0.0, float(p[2]) + 0.0), 0)
        self.constraints()
        for n, fx in self.fix.items():
            self.m.nodes[n].fix = list(fx)
        self._link_rotations()
        for n in sorted(self.masses):
            v = self.masses[n]
            self.m.tmass[n] = [float(x) for x in v[:3]]
            self.m.tmass_history.touch(n)
            if np.any(v[3:] != 0.0):
                self.m.rmass[n] = [float(x) for x in v[3:]]
                self.m.rmass_history.touch(n)
            self.m.mass_units[n] = 0
        for gid in sorted(self.groups):
            self.m.groups[gid] = self.groups[gid]
        if self.groups:
            self.m.group_active = max(self.groups)
        self.materials()
        self.components()
        self._options()
        self._report()
        return self.res

    def _below_ground(self) -> int:
        """SOLID/PLANE elements that the CHECK rule D-HOU-01 resolves to excavated soil (ETYPE 0)."""
        if not any(g.type in (1, 4) for g in self.m.groups.values()):
            return 0
        from ..prep.check import ModelView          # authoritative ETYPE rule (lazy: io must not need prep)
        return sum(1 for r in ModelView(self.m).elems if r.type in (1, 4) and r.etype == 2)

    def _rotational_nodes(self) -> set:
        """Nodes whose rotations get stiffness or inertia from a converted element other than an
        axial-only LINK beam, or from a rotary mass.

        Built from the converted SASSI model itself (not per ANSYS family), so every element that
        carries rotational DOFs counts: BEAMS (nodes I and J; K only orients), SHELL/TSHELL (the
        plate-bending rotations of every corner), SPRING with a rotational constant, GENERAL with a
        non-zero rotational row, MR masses.  SOLID and PLANE nodes have translations only."""
        out: set = set()
        for gid, g in self.groups.items():
            if gid in self.link_groups:
                continue
            for e in g.elements.values():
                if g.type in (3, 5):
                    out.update(n for n in e.nodes if n)
                elif g.type == 2:
                    out.update(e.nodes[:2])
                elif g.type == 7:
                    if any(v != 0.0 for v in self.m.springs[e.prop].k[3:]):
                        out.update(n for n in e.nodes[:2] if n)
                elif g.type == 9:
                    p = self.m.matrices[e.prop]
                    for b, n in enumerate(e.nodes[:2]):
                        rows = slice(6 * b + 3, 6 * b + 6)
                        if any(np.any(p.full(kind)[rows, :] != 0.0) for kind in ("R", "M")):
                            out.add(n)
        out.update(n for n, v in self.masses.items() if np.any(v[3:] != 0.0))
        return out

    def _link_rotations(self) -> None:
        """Rotations of nodes that only LINK-derived beams (and possibly SOLID/PLANE or translational
        springs and masses) connect are fixed: a truss node has no rotational stiffness, but SASSI
        BEAMS carry 6 DOFs per node.  A node that any other element stiffens in rotation (a shell
        corner, a beam end, ...) keeps its rotations free -- fixing them would clamp that element."""
        if not self.link_nodes:
            return
        fixed = []
        for n in sorted(self.link_nodes - self._rotational_nodes()):
            nd = self.m.nodes[n]
            if nd.fix[3:] != [1, 1, 1]:
                nd.fix[3:] = [1, 1, 1]
                fixed.append(n)
        self.info("LINK elements converted to BEAMS with axial stiffness only (I2 = I3 = 0, J = 1e-8 A^2)"
                  + (f"; ROTX ROTY ROTZ fixed at {len(fixed)} nodes that no other element stiffens in rotation "
                     f"({_short(fixed)})" if fixed else ""))

    def _options(self) -> None:
        from ..model.options import make_record
        from ..model.values import fmt_num
        house = self.m.options.ensure_record("HOUSE")
        house.set_arg(1, fmt_num(self.g))
        if self.two_d:
            house.set_arg(4, "1")
        if self.solid_incomp:
            vals = set(self.solid_incomp)
            if vals == {0}:
                self.m.options.set_record(make_record("MOPT", ["0"]))
                self.info("solid/plane element types use extra (incompatible) displacement shapes: MOPT "
                          "<incomp> = 0 (include)")
            elif len(vals) > 1:
                self.warn("solid/plane element types differ in their extra-shape (incompatible mode) option; "
                          "SASSI-EDU sets it for the whole model (MOPT <incomp>), default 1 (suppress) kept")
        if self.c.title:
            self.m.title = self.c.title

    def _report(self) -> None:
        c = self.c
        for cmd, n in sorted(c.skipped.items()):
            why = _NOT_CONVERTED.get(cmd)
            if cmd == "ENDRELEASE":
                why = "beam end releases by ENDRELEASE (couplings) cannot be represented"
            elif cmd == "ACEL":
                why = "acceleration loads (ACEL) are not converted"
            elif cmd == "N in a local CSYS":
                why = "nodes defined while a local CSYS was active were read as global coordinates"
            elif cmd == "D on ALL/selected nodes":
                why = "D on ALL/selected nodes cannot be resolved without selection logic"
            elif cmd == PARAM_KEY:
                why = ("commands with APDL parameters, expressions or names in numeric fields are not "
                       "converted (the reader does not evaluate APDL; write the model with CDWRITE), e.g. "
                       + "; ".join(repr(t) for t in c.unevaluated[:3]))
            self.warn(f"{why or cmd + ' not converted'} ({n} commands)")
        if c.unknown:
            self.warn("commands not converted: " + ", ".join(f"{k} ({v})" for k, v in sorted(c.unknown.items())))
        for t in c.notes:
            self.warn(t)
        if c.global_damping:
            self.warn("global damping " + ", ".join(f"{k} = {_g(v)}" for k, v in sorted(c.global_damping.items()))
                      + " not converted (SASSI uses material damping ratios)")
        if c.units:
            self.info(f"/UNITS,{c.units}: check that <gravity> = {_g(self.g)} is in the same units")
        used = set(self.m.used_nodes())
        massn = set(self.masses)
        unused = [n for n in self.m.nodes if n not in used and n not in massn]
        if unused:
            self.info(f"{len(unused)} nodes are not used by converted elements (AFWRITE fixes them, Warning 4): "
                      f"{_short(sorted(unused))}")
        if self.res.created_nodes:
            self.info(f"{len(self.res.created_nodes)} beam orientation (K) nodes created from the ANSYS beam axes "
                      f"(nodes {_short(self.res.created_nodes)}; D-ANS-03)")
        n_soil = self._below_ground()
        if n_soil:
            self.warn(f"{n_soil} SOLID/PLANE elements lie below the ground elevation {_g(self.m.ground_elevation)} "
                      "and would be treated as excavated soil (ETYPE 0, D-HOU-01): set ETYPE,..,1 for structure "
                      "or GROUNDELEV")
        types = {g.type for g in self.m.groups.values()}
        if types & {1, 3, 7} or self.masses:
            self.info("run FIXROT to restrain the rotations that no element stiffens (solid-only nodes, shell "
                      "drilling rotations, spring-only nodes)")
        self.info("CONVERT,ANSYS: {} nodes, {} elements in {} groups ({}), {} materials, {} R, {} SC, {} MX, "
                  "masses at {} nodes, {} nodes with fixities; converted elements have ETYPE 0 (use ETYPEGEN "
                  "or ETYPE for embedded models)".format(
                      len(self.m.nodes), self.m.n_elements(), len(self.m.groups),
                      ", ".join(f"{k} {v}" for k, v in sorted(self.res.counts.items())
                                if k not in ("D",)) or "none",
                      len(self.m.materials), len(self.m.sections), len(self.m.springs), len(self.m.matrices),
                      len(self.masses), sum(1 for n in self.m.nodes.values() if any(n.fix))))


def _isnum(t: str) -> bool:
    try:
        _num(t)
        return True
    except ValueError:
        return False


def _g(x: float) -> str:
    return f"{x:.6g}"


def _rlist(positions: Sequence[int], n: int = 8) -> str:
    """``R13, R25, ...`` (1-based real-constant positions)."""
    s = ", ".join(f"R{k}" for k in positions[:n])
    return s + (", ..." if len(positions) > n else "")


def _short(ids: Sequence[int], n: int = 10) -> str:
    ids = list(ids)
    s = ", ".join(str(i) for i in ids[:n])
    return s + (f", ... ({len(ids)} in total)" if len(ids) > n else "")


def _release_digits(code: int) -> List[int]:
    """BEAM44 KEYOPT(7)/(8): six digits for UX UY UZ ROTX ROTY ROTZ (left to right, leading zeros
    implied; ANSYS lists ROTZ in column 1 = the units digit), 1 = released."""
    s = f"{int(code):06d}"[-6:]
    return [1 if ch != "0" else 0 for ch in s]


def _sassi_release(ans: Sequence[int]) -> List[int]:
    """ANSYS local release flags (UX UY UZ ROTX ROTY ROTZ) -> SASSI KI/KJ (P1 P2 P3 M1 M2 M3).

    SASSI axis 2 = ANSYS z, axis 3 = -y (D-ANS-03): P1 = UX, P2 = UZ, P3 = UY, M1 = ROTX,
    M2 = ROTZ, M3 = ROTY.
    """
    return [ans[0], ans[2], ans[1], ans[3], ans[5], ans[4]]


def convert_cdb(cdb: CdbModel, gravity: float, damping: Optional[float] = None) -> CdbConversion:
    """Convert parsed ANSYS data to a new :class:`SSIModel` (CONVERT,ANSYS; D-ANS-01..05).

    ``gravity`` turns DENS into specific weight (``weight = DENS * g``) and becomes the model
    GRAVITY (HOUSE ``<gravity>``).  ``damping`` is the damping ratio of materials without DMPR
    (or DMPS); ``None`` means 0 (reported).  Groups: one per ANSYS element type number (extra
    groups numbered after the largest type number for oblique COMBIN14 springs), elements
    renumbered 1..n in file order, node numbers kept, K nodes added after the largest node.
    """
    if gravity is None or gravity <= 0:
        raise ValueError("gravity must be > 0 (it converts DENS to specific weight)")
    return _Converter(cdb, gravity, damping).run()


def convert_cdb_file(path: PathLike, gravity: float, damping: Optional[float] = None) -> CdbConversion:
    """:func:`read_cdb` + :func:`convert_cdb`."""
    return convert_cdb(read_cdb(path), gravity, damping)
