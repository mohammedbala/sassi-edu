"""CHECK: the model and analysis-option checker (manual Chapter 10; requirements sections 1.9, 5.5,
D-CHK-01..09; spec 11 section 5; spec 07 section 6).

CHECK "simulates writing the analysis files": it checks the model and the options of every
module enabled by AOPT and lists the numbered messages of the manual's catalogue,

* ``Error 1`` ... ``Error 128`` and ``Warning 1`` ... ``Warning 11`` with the **exact titles** of the
  manual (they are the nomenclature, spec 11 section 5.2/5.3);
* the additional checks of this implementation ``EDU-nn`` (D-CHK-09), numbered separately so the
  manual numbering stays faithful.  EDU-01 ... EDU-22 are the decision-log checks (EDU-09, BBC
  slopes of Option NON, is P2 and EDU-18, the post-run ATF test, belongs to the module listings);
  this package adds the consistency checks EDU-23 duplicate frequency numbers (error, G-01),
  EDU-24 frequency numbers above NFFT/2 (warning), EDU-25 SOIL sublayer not defined by SPRO
  (error), EDU-26 option value or combination not allowed (error: D-PNT-01, G-17, coherency models
  2-7 / ME without wave passage, simultaneous cases with an angle, coherency models whose coefficients
  are not available (D-INC-04), missing user coherency files, incompatible stochastic / superposition
  options and simultaneous-case counts of incoherent runs), EDU-27 shell material with
  different P- and S-damping (warning, requirements 4.0.2) and EDU-28 recommended practice not
  followed (warning: D-SIT-03, D-SOL-02, SRP 3.7.1 duration, an L number repeated among the embedment
  layers of TOPL).  The EXCSTRCHK condition (G-15) is
  reported as ``Error EXCSTRCHK`` -- or ``Warning EXCSTRCHK`` when the shared interior node is an
  interaction node used only by SOLID/PLANE elements built on the excavation mesh (near-field soil, a
  solid basement, the zero-SSI identity): the rule of HOUSE, :func:`sassi.core.house_lib.excstrchk_kind`
  (final audit).

Attribution (D-CHK-01): model-level messages are grouped under a ``MODEL`` header; their errors
gate HOUSE and every module that reads FILE4 (ANALYS, STRESS).  Errors on shared variables
(gravity, time step, NFFT, frequency set, control-motion data ...) are reported under every enabled
module that uses them; the free-field layer profile (TOPL layers, SITE half-space) is checked
under SITE and under HOUSE, whose deck copies it.  Errors 42 (no nodes) and 43 (no groups) are
fatal: CHECK stops.

AFWRITE (:mod:`sassi.prep.afwrite`) runs CHECK and writes the deck of each enabled module that has
no error (D-AFW-01); warnings never block.

Geometry tolerance (D-GEN-06): ``tol = max(1e-6 L_ref, 1e-9)`` with ``L_ref`` the bounding-box
diagonal of the model (``EDUOPT,GEOMTOL`` overrides); an interaction node lies on a soil-layer
interface when ``|depth - d_i| <= max(tol, 1e-4 h_min)``.

The model-checking commands EXCSTRCHK, FIXEDINT, FREESPRING, HINGED, INTCOUNT, KINT and USED, the
FIXROT family and AFWRITE reuse the topology helpers of :class:`ModelView`.
"""
from __future__ import annotations

import re

import math
import os
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

import numpy as np

from ..conventions import ELEMENT_TYPE_NAMES, is_power_of_two, nearest_power_of_two
from ..core.house_lib import excstrchk_kind
from ..model.entities import MIN_ELEMENT_NODES
from ..model.materials import elastic_constants
from ..model.values import fmt_num
from .options import (AOPT_MODULES, eduopt, eduopt_float, get_entries, get_record, history_problems,
                      problem_fmt, time_grid_problems)

SOLID, BEAMS, SHELL, PLANE, TSHELL, SPRING, GENERAL = 1, 2, 3, 4, 5, 7, 9

#: modules in run order (AOPT order; requirements section 5.5 "module headers in run order")
RUN_ORDER = tuple(m for m in AOPT_MODULES if m is not None)
#: modules gated by model-level errors (HOUSE and the modules reading FILE4, D-CHK-01)
MODEL_GATED = ("HOUSE", "ANALYS", "STRESS")
MODEL = "MODEL"

# ======================================================================================
# Message catalogue (manual Chapter 10; titles exact, spec 11 section 5.2 / 5.3)
# ======================================================================================
ERRORS: Dict[int, str] = {
    1: "Illegal Acceleration of Gravity",
    2: "Illegal Node for Symmetry Plane/Line {i}",
    3: "Illegal Number of Nodes for Symmetry Plane/Line {i}",
    4: "Nodes from Symmetry Line {i} Are Equal",
    5: "Nodes from Symmetry Plane {i} Are Collinear",
    6: "Undefined Element {e}, Group {g}",
    7: "Element {e} from Group {g} Has Too Few Nodes",
    8: "Element {e}, Group {g} Has 0 Length",
    9: "Nodes from Element {e}, Group {g} Are Collinear",
    10: "Element {e}, Group {g} Has Improper Release Code",
    11: "Zero Area for Element {e}, Group {g}",
    12: "Warped Element {e}, Group {g}",
    13: "Material {m} Is not Defined",
    14: "Elasticity Modulus from Material {m} Is Illegal",
    15: "Poisson Coefficient from Material {m} Is Illegal",
    16: "Specific Weight from Material {m} Is Illegal",
    17: "P-Wave Damping Ratio from Material {m} Is Illegal",
    18: "S-Wave Damping Ratio from Material {m} Is Illegal",
    19: "Soil Layer {l} Is not Defined",
    20: "Thickness from Soil Layer {l} Is Illegal",
    21: "Specific Weight from Soil Layer {l} Is Illegal",
    22: "P-Wave Velocity from Soil Layer {l} Is Illegal",
    23: "S-Wave Velocity from Soil Layer {l} Is Illegal",
    24: "P-Wave Damping Ratio from Soil Layer {l} Is Illegal",
    25: "S-Wave Damping Ratio from Soil Layer {l} Is Illegal",
    26: "Property {p} Is not Defined",
    27: "Axial Area from Property {p} Is Illegal",
    28: "Shear Area 2 from Property {p} Is Illegal",
    29: "Shear Area 3 from Property {p} Is Illegal",
    30: "Torsion Inertia Moment from Property {p} Is Illegal",
    31: "Flexural Inertia Moment 2 from Property {p} Is Illegal",
    32: "Flexural Inertia Moment 3 from Property {p} Is Illegal",
    33: "Spring Property {p} Is not Defined",
    34: "Spring Constant X from Spring Property {p} Is Illegal",
    35: "Spring Constant Y from Spring Property {p} Is Illegal",
    36: "Spring Constant Z from Spring Property {p} Is Illegal",
    37: "Spring Constant XX from Spring Property {p} Is Illegal",
    38: "Spring Constant YY from Spring Property {p} Is Illegal",
    39: "Spring Constant ZZ from Spring Property {p} Is Illegal",
    40: "Mass out of Defined Nodes Range",
    41: "Node {n} Is Not Defined - Group {g} Element {e}",
    42: "No Nodes Defined",
    43: "No Groups Defined",
    44: "Frequency Set {s} Is Not Defined",
    45: "Mode 1 and Mode 2 Are Both Deselected",
    46: "No Top Layers",
    47: "Illegal Number of Layers for Halfspace Simulation",
    48: "Illegal Frequency Step",
    49: "Illegal Time Step of Control Motion",
    50: "Illegal Number of Values for Fourier Transform",
    51: "All Wave Fields Are Deselected",
    52: "Illegal Incident Angle of Wave {w}",
    53: "Illegal Value for Frequency {i}",
    54: "Illegal Value for Wave {w} Ratio at Frequency {i}",
    55: "Illegal Sum of Wave Ratios at Frequency {i}",
    56: "Illegal Last Layer Number in Near Field Zone",
    57: "Illegal Radius of Central Zone",
    58: "Illegal Coherence Parameter",
    59: "Illegal Mean Soil Shear Wave Velocity",
    60: "Illegal Number of Mesh Points / Embedment Level",
    61: "No Forces Defined",
    62: "Node for Force / Moment {i} is not defined",
    63: "Illegal Coordinate Transformation Angle",
    64: "No Nodal Output Request",
    65: "Illegal Nodal Output Request: {n}",
    66: "Nodal Output Request Defined More Than Once: {n}",
    67: "Illegal Output Time History Step",
    68: "Illegal Total Duration To Be Plotted",
    69: "Illegal First Frequency for RS Analysis",
    70: "Illegal Last Frequency for RS Analysis",
    71: "Illegal Number of Frequency Steps For RS Analysis",
    72: "Illegal Damping Ratio For RS Analysis",
    73: "Acceleration Time History File Does Not Exist",
    74: "Illegal First Record Number",
    75: "Illegal Last Record Number",
    76: "First Record Number Larger Than Last Record Number",
    77: "Multiplication Factor and Maximum Value Of Time History Are Both Zero",
    78: "Multiplication Factor and Maximum Value Of Time History Are Both Non-Zero",
    79: "No Element Output Request",
    80: "Illegal Group For Output Request: {g}",
    81: "Illegal Element Output Request: {e}, Group {g}",
    82: "Element Output Request Defined More Than Once: {e}, Group {g}",
    83: "Matrix Property {p} Is not Defined",
    84: "No RS Input Files Specified",
    85: "RS Input File {i} Does Not Exist",
    86: "Invalid RS Output File {i}",
    87: "Invalid Acceleration Output File {i}",
    88: "Invalid Acceleration Input File {i}",
    89: "Number of Frequencies Does Not Match RS Input File {i}",
    90: "Illegal Initial Random Number",
    91: "Illegal Number of Frequencies",
    92: "Illegal Duration",
    93: "No Correlation Factors Defined",
    94: "Illegal Correlation Factor",
    95: "No Dynamic Soil Properties Assigned",
    96: "Too Many Dynamic Soil Properties",
    97: "Dynamic property {p} has no shear modulus curve",
    98: "Dynamic property {p} has illegal shear modulus values",
    99: "Dynamic property {p} has no damping curve",
    100: "Number of Acceleration Values Is Illegal",
    101: "Cut-Off Frequency Is Illegal",
    102: "Illegal Reading Format",
    103: "Illegal Number of Header Lines",
    104: "Illegal Control Layer Number",
    105: "Illegal Number of Iterations",
    106: "Illegal Strain Ratio",
    107: "No Damping Ratios Defined",
    108: "Illegal Multiplier for Acceleration of Gravity",
    109: "Illegal Second Layer Number for Layer {i}",
    110: "Illegal Frequency Step for Layer {i}",
    111: "Illegal Number of Smoothings for Layer {i}",
    112: "Illegal Number of Values to Be Saved for Layer {i}",
    113: "Illegal Apparent Velocity for Line D",
    114: "Illegal Directional Coherence Factor",
    115: "No Multiple Excitation Data Defined",
    116: "Illegal First Node Number for Motion {i}",
    117: "Illegal Last Node Number for Motion {i}",
    118: "Illegal Spectral Amplification Ratio for Motion {i}",
    119: "Spectral Amplification Ratios for Motion {i} Do Not Match Frequencies",
    120: "Frequency set {i} is Empty. Change the active frequency set in the SITE tab of the Analysis Options",
    121: "No Panels Specified",
    122: "Material {i} referenced in panel {k} does not exist",
    123: "Group {i} referenced in panel {k} does not exist",
    124: "Node {i} is a fixed interaction node",
    125: "Not enough nonlinear soil properties have been defined for the Nonlinear soil model.",
    126: "Force Option {i} not supported for panel {k} in this version",
    127: "Force Option {i} not supported for spring {k}",
    128: "Displacement Option {i} not supported for panel {k}",
}

WARNINGS: Dict[int, str] = {
    1: "Gap Found at Node {n}",
    2: "Distorted Element {e}, Group {g}, Face {f}",
    3: "Warped Element {e}, Group {g}",
    4: "Unused Node {n}",
    5: "Translational Mass in Node {n} Is on Fixed DOF",
    6: "Rotational Mass in Node {n} Is on Fixed DOF",
    7: "Group {g} Has no Elements",
    8: "Too Many Top Layers",
    9: "Number of Values for Fourier Transform Is Not Power of 2",
    10: "Force in Node {n} Is on Fixed DOF",
    11: "Moment in Node {n} Is on Fixed DOF",
}

#: SASSI-EDU checks (D-CHK-09 and package additions): label -> (kind, title)
EDU_TEXT: Dict[str, Tuple[str, str]] = {
    "EDU-01": ("Error", "Interaction Node {n} Is Not on a Soil Layer Interface"),
    "EDU-02": ("Warning", "Size Limit of ACS SASSI V3 Exceeded"),
    "EDU-03": ("Error", "Rounded Number of Fourier Components Is Shorter Than the Records Used"),
    "EDU-04": ("Error", "Damping Ratio >= 0.5"),
    "EDU-05": ("Error", "Non-Positive Jacobian in Element {e}, Group {g}"),
    "EDU-06": ("Warning", "Unrestrained Shell Drilling Rotation at Node {n}"),
    "EDU-07": ("Error", "FILE88 Layers Do Not Match the Top Layers"),
    "EDU-08": ("Warning", "Excavated Element {e}, Group {g} Layer Differs from the Top Layer at Its Depth"),
    "EDU-10": ("Warning", "Passing Frequency Below the Cut-Off Frequency"),
    "EDU-11": ("Warning", "Poisson Ratio > 0.47 in Soil Layer {l}"),
    "EDU-12": ("Warning", "Non-FV Method Selected: Validate Against FV (ASCE 4-16, SRP 3.7.2)"),
    "EDU-13": ("Warning", "Deterministic AS/SRSS Incoherency Is Valid for Rigid Foundations Only"),
    "EDU-14": ("Warning", "Smoothing Parameter Not Zero"),
    "EDU-15": ("Warning", "Fewer Than 301 Response-Spectrum Frequencies"),
    "EDU-16": ("Warning", "Nodal Output Request Defined More Than Once: {n} (merged)"),
    "EDU-17": ("Warning", "Near-Coincident Horizontal Projections of Interaction Nodes"),
    "EDU-19": ("Warning", "Quiet Zone Too Short"),
    "EDU-20": ("Warning", "Memory Estimate of the Impedance Matrices"),
    "EDU-21": ("Error", "Interaction Nodes Not Numbered Bottom-Up"),
    "EDU-22": ("Warning", "Excavation Boundary Node {n} Is Not an Interaction Node"),
    "EDU-23": ("Error", "Duplicate Frequency Numbers in Frequency Set {s}"),
    "EDU-24": ("Warning", "Frequency Number Above NFFT/2 in Frequency Set {s}"),
    "EDU-25": ("Error", "SOIL Sublayer {i} Is Not Defined (SPRO)"),
    "EDU-26": ("Error", "Option Combination or Value Not Allowed"),
    "EDU-27": ("Warning", "Shell Material {m} Has Different P- and S-Wave Damping"),
    "EDU-28": ("Warning", "Recommended Practice Not Followed"),
    # non-linear soil SSI (requirements 4.4 item 7, 4.10 item 5)
    "EDU-41": ("Error", "Non-Linear Soil Group {g} Is Not Valid"),
    "EDU-42": ("Warning", "Non-Linear Soil SSI Input (.pin) Not Defined"),
    "EDU-43": ("Warning", "Non-Linear Soil SSI and STRESS Strain Computation Do Not Match"),
    "EDU-44": ("Error", "Option NON Input Not Valid"),
    "EDU-45": ("Warning", "Option NON Input Note"),
    "EXCSTRCHK": ("Error", "Excavation Interior Node {n} Is Shared with the Structure"),
}
#: severity of EDU-21 when the analysis is coherent (requirements section 4.4 item 3)
EDU21_COHERENT_KIND = "Warning"


# ======================================================================================
# Report
# ======================================================================================
@dataclass
class CheckOptions:
    """Options > Check (requirements section 5.4): session settings, not saved (UI-07)."""
    show_warnings: bool = True
    show_errors: bool = True
    suppress_window: bool = False
    break_at: int = 100              # per message type, per module (D-CHK-02)


@dataclass
class CheckMessage:
    module: str
    kind: str                        # 'Error' | 'Warning'
    number: Union[int, str]
    text: str
    detail: str = ""

    def line(self) -> str:
        s = f"{self.kind} {self.number} : {self.text}"
        return s + (f"  [{self.detail}]" if self.detail else "")


@dataclass
class CheckReport:
    """Result of CHECK: messages per module plus the gating information used by AFWRITE."""
    model: str = ""
    modules: List[str] = field(default_factory=list)      # enabled modules, run order
    sections: List[str] = field(default_factory=list)     # printed headers (MODEL first when checked)
    messages: List[CheckMessage] = field(default_factory=list)
    fatal: bool = False
    notes: List[str] = field(default_factory=list)

    # ---------------------------------------------------------------- adding
    def add(self, module: str, kind: str, number: Union[int, str], detail: str = "", **fmt) -> CheckMessage:
        if isinstance(number, int):
            title = (ERRORS if kind == "Error" else WARNINGS)[number]
        else:
            title = EDU_TEXT[number][1]
        text = format_title(title, fmt)
        msg = CheckMessage(module, kind, number, text, detail)
        self.messages.append(msg)
        if module not in self.sections:
            self.sections.append(module)
        return msg

    def error(self, module: str, number: Union[int, str], detail: str = "", **fmt) -> CheckMessage:
        return self.add(module, "Error", number, detail, **fmt)

    def warning(self, module: str, number: Union[int, str], detail: str = "", **fmt) -> CheckMessage:
        return self.add(module, "Warning", number, detail, **fmt)

    def edu(self, module: str, label: str, detail: str = "", kind: Optional[str] = None, **fmt) -> CheckMessage:
        return self.add(module, kind or EDU_TEXT[label][0], label, detail, **fmt)

    # ---------------------------------------------------------------- queries
    def of(self, module: Optional[str] = None, kind: Optional[str] = None) -> List[CheckMessage]:
        return [m for m in self.messages if (module is None or m.module == module) and
                (kind is None or m.kind == kind)]

    def errors(self, module: Optional[str] = None) -> List[CheckMessage]:
        return self.of(module, "Error")

    def warnings(self, module: Optional[str] = None) -> List[CheckMessage]:
        return self.of(module, "Warning")

    def numbers(self, kind: Optional[str] = None, module: Optional[str] = None) -> Set[Union[int, str]]:
        return {m.number for m in self.of(module, kind)}

    def has(self, kind: str, number: Union[int, str], module: Optional[str] = None) -> bool:
        return any(m.kind == kind and m.number == number for m in self.of(module))

    def blocked(self, module: str) -> bool:
        """True when the module's deck must not be written (D-AFW-01, D-CHK-01)."""
        if self.fatal:
            return True
        if self.errors(module):
            return True
        return module in MODEL_GATED and bool(self.errors(MODEL))

    def blocking_reason(self, module: str) -> str:
        if self.fatal:
            return "CHECK stopped on a fatal error"
        n = len(self.errors(module))
        if n:
            return f"{n} CHECK error(s)"
        if module in MODEL_GATED and self.errors(MODEL):
            return f"{len(self.errors(MODEL))} model error(s) (D-CHK-01)"
        return ""

    # ---------------------------------------------------------------- formatting
    def format(self, options: Optional[CheckOptions] = None, totals: bool = True) -> str:
        """The Check Errors window / ``.err`` text (requirements section 5.5).

        Options > Check "Break Check at N" (spec 05a section 3, OQ19; D-CHK-02): at most N messages
        of each *type* -- N errors and N warnings -- are printed per module (section); the hidden
        messages are counted per type and the totals still count every message.
        """
        opt = options or CheckOptions()
        lines = [f"CHECK: Errors and Warning for - {self.model or '(unnamed model)'}"]
        order = [s for s in [MODEL] + list(RUN_ORDER) if s in self.sections]
        for sec in order:
            lines.append(f"Errors and Warnings for {sec}")
            counts: Dict[str, int] = defaultdict(int)
            hidden: Dict[str, int] = defaultdict(int)
            for m in self.of(sec):
                if (m.kind == "Error" and not opt.show_errors) or (m.kind == "Warning" and not opt.show_warnings):
                    continue
                counts[m.kind] += 1
                if opt.break_at and counts[m.kind] > opt.break_at:
                    hidden[m.kind] += 1
                    continue
                lines.append(m.line())
            for kind in ("Error", "Warning"):
                if hidden.get(kind):
                    lines.append(f"  ... {hidden[kind]} more {kind.lower()}s not shown (Break Check at {opt.break_at})")
            if totals:
                lines.append(f"{sec}: {len(self.errors(sec))} errors, {len(self.warnings(sec))} warnings")
        if self.fatal:
            lines.append("CHECK stopped: fatal error (Error 42 / 43)")
        for n in self.notes:
            lines.append(n)
        return "\n".join(lines) + "\n"

    def summary(self) -> str:
        ne, nw = len(self.errors()), len(self.warnings())
        blocked = [m for m in self.modules if self.blocked(m)]
        s = f"CHECK: {ne} errors, {nw} warnings"
        if blocked:
            s += f"; modules with errors: {', '.join(blocked)}"
        return s


def write_err(report: CheckReport, path: Union[str, Path], options: Optional[CheckOptions] = None) -> Path:
    """Write the ``<model>.err`` file (requirements section 5.5)."""
    p = Path(path)
    p.write_text(report.format(options), encoding="utf-8")
    return p


# ======================================================================================
# Geometry helpers
# ======================================================================================
def _norm(v) -> float:
    return float(np.sqrt(np.dot(v, v)))


def newell_normal(P: Sequence[np.ndarray]) -> np.ndarray:
    """Polygon normal by Newell's method (robust for slightly warped quads)."""
    n = np.zeros(3)
    k = len(P)
    for i in range(k):
        a, b = P[i], P[(i + 1) % k]
        n += np.cross(a, b)
    return n


def polygon_angles(P: Sequence[np.ndarray]) -> List[float]:
    """Interior angles (degrees) of a polygon; reflex angles > 180 are detected with the normal."""
    k = len(P)
    nrm = newell_normal(P)
    ln = _norm(nrm)
    nh = nrm / ln if ln > 0 else nrm
    out = []
    for i in range(k):
        a = P[i - 1] - P[i]
        b = P[(i + 1) % k] - P[i]
        s = float(np.dot(nh, np.cross(b, a)))
        c = float(np.dot(b, a))
        ang = math.degrees(math.atan2(s, c)) if ln > 0 else math.degrees(math.atan2(_norm(np.cross(b, a)), c))
        out.append(ang % 360.0)
    return out


def quad_angles(Q: np.ndarray) -> np.ndarray:
    """Interior angles (degrees) of quadrilaterals ``Q`` (m, 4, 3), vectorised (reflex angles > 180)."""
    nrm = sum(np.cross(Q[:, i], Q[:, (i + 1) % 4]) for i in range(4))
    ln = np.linalg.norm(nrm, axis=1)
    nh = nrm / np.where(ln > 0, ln, 1.0)[:, None]
    a = np.roll(Q, 1, axis=1) - Q            # to the previous vertex
    b = np.roll(Q, -1, axis=1) - Q           # to the next vertex
    s = np.einsum("mk,mik->mi", nh, np.cross(b, a))
    c = np.einsum("mik,mik->mi", b, a)
    return np.degrees(np.arctan2(s, c)) % 360.0


def dedupe_polygon(ids: Sequence[int], P: Dict[int, np.ndarray], tol: float) -> Tuple[List[int], List[np.ndarray]]:
    """Remove consecutive (cyclic) repeated or coincident nodes (triangles entered as quads)."""
    oi: List[int] = []
    op: List[np.ndarray] = []
    for n in ids:
        p = P[n]
        if oi and (n == oi[-1] or _norm(p - op[-1]) <= tol):
            continue
        oi.append(n)
        op.append(p)
    while len(oi) > 1 and (oi[0] == oi[-1] or _norm(op[0] - op[-1]) <= tol):
        oi.pop()
        op.pop()
    return oi, op


#: SOLID faces (1-based slots): bottom, top, then the four sides (spec 08 Figure 9.1 node order)
SOLID_FACES = ((1, 2, 3, 4), (5, 6, 7, 8), (1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (4, 1, 5, 8))
PLANE_EDGES = ((1, 2), (2, 3), (3, 4), (4, 1))

_GP = 1.0 / math.sqrt(3.0)


def _hex_dN() -> np.ndarray:
    """Derivatives of the trilinear shape functions at the 2x2x2 Gauss points: (8 gp, 3, 8 nodes)."""
    nat = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    gps = nat * _GP
    out = np.zeros((8, 3, 8))
    for g, (x, y, z) in enumerate(gps):
        for a, (xa, ya, za) in enumerate(nat):
            out[g, 0, a] = xa * (1 + y * ya) * (1 + z * za) / 8.0
            out[g, 1, a] = ya * (1 + x * xa) * (1 + z * za) / 8.0
            out[g, 2, a] = za * (1 + x * xa) * (1 + y * ya) / 8.0
    return out


def _quad_dN() -> np.ndarray:
    nat = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)
    gps = nat * _GP
    out = np.zeros((4, 2, 4))
    for g, (x, y) in enumerate(gps):
        for a, (xa, ya) in enumerate(nat):
            out[g, 0, a] = xa * (1 + y * ya) / 4.0
            out[g, 1, a] = ya * (1 + x * xa) / 4.0
    return out


_HEX_DN = _hex_dN()
_QUAD_DN = _quad_dN()


def hex_jacobians(X: np.ndarray) -> np.ndarray:
    """det J at the 2x2x2 Gauss points of SOLID elements; ``X`` (n, 8, 3) -> (n, 8)."""
    J = np.einsum("gia,naj->ngij", _HEX_DN, X)
    return np.linalg.det(J)


def quad_jacobians_xz(X: np.ndarray) -> np.ndarray:
    """det J at the 2x2 Gauss points of PLANE elements in the X-Z plane; ``X`` (n, 4, 3) -> (n, 4)."""
    XZ = X[:, :, [0, 2]]
    J = np.einsum("gia,naj->ngij", _QUAD_DN, XZ)
    return np.linalg.det(J)


# ======================================================================================
# Model view: coordinates, topology, resolved ETYPE, excavation sets
# ======================================================================================
@dataclass
class ElemRef:
    """One element with its group data (resolved ETYPE, the nodes carrying DOFs, the K node)."""
    group: int
    type: int
    elem: Any                         # sassi.model.Element
    etype: int = 1                    # resolved: 1 structure, 2 excavated soil / buried shell

    @property
    def id(self) -> int:
        return self.elem.id

    @property
    def nodes(self) -> List[int]:
        return [n for n in self.elem.nodes]

    @property
    def dof_nodes(self) -> List[int]:
        """Nodes carrying DOFs: BEAMS and 3-node GENERAL elements use I and J (K is orientation)."""
        ns = [n for n in self.elem.nodes if n]
        if self.type in (BEAMS, GENERAL):
            return [n for n in self.elem.nodes[:2] if n]
        return ns

    @property
    def k_node(self) -> int:
        if self.type in (BEAMS, GENERAL) and len(self.elem.nodes) >= 3:
            return self.elem.nodes[2]
        return 0

    @property
    def excavated(self) -> bool:
        return self.type in (SOLID, PLANE) and self.etype == 2


class ModelView:
    """Derived views of a model used by CHECK, AFWRITE and the model-checking commands.

    * ``P[n]``: global coordinates of node n (local systems resolved);
    * ``tol``: geometric tolerance (D-GEN-06), ``lref`` the bounding-box diagonal;
    * ``elems``: :class:`ElemRef` of every element, ETYPE 0 resolved by D-HOU-01;
    * ``node_elems[n]``: elements carrying DOFs at node n; ``k_nodes``: beam / GENERAL K nodes.
    """

    def __init__(self, model):
        self.m = model
        self.P: Dict[int, np.ndarray] = {}
        try:
            ids, xyz = model.global_coordinates()            # vectorised per coordinate system
            self.P = {int(i): xyz[k] for k, i in enumerate(ids)}
        except KeyError:                                      # a node in an undefined system
            for nid, n in model.nodes.items():
                try:
                    self.P[nid] = model.to_global(n.xyz, n.csys)
                except KeyError:
                    self.P[nid] = np.array([np.nan, np.nan, np.nan])
        if self.P:
            A = np.array(list(self.P.values()))
            A = A[np.all(np.isfinite(A), axis=1)]
            self.lref = float(np.linalg.norm(A.max(axis=0) - A.min(axis=0))) if len(A) else 0.0
        else:
            self.lref = 0.0
        gt = eduopt_float(model, "GEOMTOL", 0.0)
        self.tol = gt if gt > 0 else max(1e-6 * self.lref, 1e-9)
        self.gelev = model.ground_elevation
        self.elems: List[ElemRef] = []
        self.node_elems: Dict[int, List[ElemRef]] = defaultdict(list)
        self.k_nodes: Set[int] = set()
        for g, e in model.iter_elements():
            r = ElemRef(g.id, g.type, e)
            r.etype = self.resolve_etype(r)
            self.elems.append(r)
            for n in set(r.dof_nodes):
                self.node_elems[n].append(r)
            if r.k_node:
                self.k_nodes.add(r.k_node)
        self._exc: Optional[Dict[str, Any]] = None

    # ---------------------------------------------------------------- ETYPE (D-HOU-01)
    def defined(self, nodes: Iterable[int]) -> bool:
        return all(n in self.P for n in nodes if n)

    def resolve_etype(self, r: ElemRef) -> int:
        """ETYPE 0 resolved (D-HOU-01, D-AFW-06): SOLID/PLANE excavated if every node has
        ``z <= gelev + tol`` and the centroid is strictly below ``gelev``; otherwise structure.
        SHELL/TSHELL ETYPE 2 = buried shell; BEAMS/SPRING/GENERAL are structure."""
        e = r.elem
        if r.type in (SOLID, PLANE):
            if e.etype in (1, 2):
                return e.etype
            ns = [n for n in e.nodes if n]
            if not ns or not self.defined(ns):
                return 1
            z = np.array([self.P[n][2] for n in ns])
            uniq = np.array([self.P[n][2] for n in dict.fromkeys(ns)])
            if np.all(z <= self.gelev + self.tol) and float(uniq.mean()) < self.gelev - self.tol:
                return 2
            return 1
        if r.type in (SHELL, TSHELL):
            return 2 if e.etype == 2 else 1
        return 1

    # ---------------------------------------------------------------- sets
    def interaction_nodes(self) -> List[int]:
        return sorted(n for n, nd in self.m.nodes.items() if 0 in nd.flags)

    def used_dof_nodes(self) -> Set[int]:
        return {n for n, lst in self.node_elems.items() if lst}

    def node_types(self, n: int) -> Set[int]:
        return {r.type for r in self.node_elems.get(n, [])}

    def excavation(self) -> Dict[str, Any]:
        """Excavation sets of spec 09 section 2.19 (resolved ETYPE): elements X, nodes N_X, boundary
        faces B (faces used by exactly one element of X), boundary nodes N_B, FSIN nodes (lateral and
        bottom boundary), interior nodes N_X minus N_B."""
        if self._exc is not None:
            return self._exc
        X = [r for r in self.elems if r.excavated and self.defined(r.nodes)]
        NX: Set[int] = set()
        faces: Dict[Tuple[int, ...], List[Tuple[ElemRef, Tuple[int, ...]]]] = defaultdict(list)
        for r in X:
            ns = r.elem.nodes
            NX.update(n for n in ns if n)
            if r.type == SOLID:
                pad = list(ns) + [0] * (8 - len(ns))
                for f in SOLID_FACES:
                    fn = tuple(pad[k - 1] for k in f)
                    u = tuple(sorted(set(x for x in fn if x)))
                    if len(u) >= 3:
                        faces[u].append((r, fn))
            else:
                pad = list(ns) + [0] * (4 - len(ns))
                if len(ns) == 3 or not pad[3]:
                    edges = ((1, 2), (2, 3), (3, 1))
                else:
                    edges = PLANE_EDGES
                for f in edges:
                    fn = tuple(pad[k - 1] for k in f)
                    u = tuple(sorted(set(x for x in fn if x)))
                    if len(u) == 2:
                        faces[u].append((r, fn))
        boundary = {k for k, v in faces.items() if len(v) == 1}
        NB: Set[int] = set()
        top: Set[Tuple[int, ...]] = set()
        for k in boundary:
            NB.update(k)
            if all(abs(self.P[n][2] - self.gelev) <= self.tol for n in k):
                top.add(k)
        fsin: Set[int] = set()
        for k in boundary - top:
            fsin.update(k)
        self._exc = dict(elements=X, nodes=NX, boundary_faces=boundary, boundary=NB, top_faces=top,
                         fsin=fsin, interior=NX - NB)
        return self._exc

    def buried_shell_nodes(self) -> Set[int]:
        return {n for r in self.elems if r.type in (SHELL, TSHELL) and r.etype == 2 for n in r.nodes if n}

    # ---------------------------------------------------------------- shells
    def shell_normal(self, r: ElemRef) -> Optional[np.ndarray]:
        ns = [n for n in r.elem.nodes if n]
        if not self.defined(ns):
            return None
        ids, pts = dedupe_polygon(ns, self.P, self.tol)
        if len(pts) < 3:
            return None
        nv = newell_normal(pts)
        ln = _norm(nv)
        return nv / ln if ln > 0 else None

    def coplanar_normal(self, n: int, eps: float = 1e-4) -> Optional[np.ndarray]:
        """Common unit normal of the SHELL elements at node n (None if not coplanar or no shells)."""
        shells = [r for r in self.node_elems.get(n, []) if r.type == SHELL]
        if not shells:
            return None
        ref = None
        for r in shells:
            nv = self.shell_normal(r)
            if nv is None:
                return None
            if ref is None:
                ref = nv
            elif abs(float(np.dot(nv, ref))) < 1.0 - eps:
                return None
        return ref

    # ---------------------------------------------------------------- soil profile
    def interfaces(self) -> Optional[np.ndarray]:
        """Depths of the user interfaces 1..nTOPL+1 (from the TOPL layers), or None if undefined."""
        m = self.m
        if not m.topl or any(l not in m.layers for l in m.topl):
            return None
        h = np.array([m.layers[l].thick for l in m.topl], float)
        if np.any(h <= 0):
            return None
        return np.concatenate([[0.0], np.cumsum(h)])

    def interface_of(self, z: float, depths: np.ndarray) -> Tuple[int, float]:
        """(1-based interface, |misfit|) nearest to elevation z."""
        d = self.gelev - z
        k = int(np.argmin(np.abs(depths - d)))
        return k + 1, float(abs(depths[k] - d))

    def interface_tol(self, depths: np.ndarray) -> float:
        hmin = float(np.min(np.diff(depths))) if len(depths) > 1 else 0.0
        return max(self.tol, 1e-4 * hmin)


# ======================================================================================
# Small utilities
# ======================================================================================
def physical_ram() -> Optional[int]:
    """Physical memory in bytes (None when unknown)."""
    try:
        return int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))
    except (ValueError, OSError, AttributeError):
        return None


def dense_matrix_bytes(n_int: int) -> int:
    """Bytes of one dense complex impedance / flexibility matrix of ``n_int`` interaction nodes:
    ``(3 N)^2 x 16`` (UT-18: 10,000 nodes -> 14.4 GB)."""
    return (3 * int(n_int)) ** 2 * 16


def analys_memory_bytes(n_int: int) -> int:
    """ANALYS estimate of D-ANL-10: three dense matrices ``3 (3 N)^2 x 16`` bytes."""
    return 3 * dense_matrix_bytes(n_int)


def resolve_file(name: str, dirs: Sequence[Union[str, Path]]) -> Optional[Path]:
    """Resolve a file name: absolute, else the first directory of ``dirs`` where it exists."""
    if not name:
        return None
    raw = os.path.expanduser(name.strip())
    variants = [raw] + ([raw.replace("\\", "/")] if "\\" in raw and os.sep == "/" else [])
    for v in variants:
        p = Path(v)
        if p.is_absolute():
            if p.exists():
                return p
            continue
        for d in dirs:
            if d and (Path(d) / p).exists():
                return Path(d) / p
    return None


def count_records(path: Path, fopt: int = 0, header: int = 0) -> Optional[int]:
    """Number of acceleration records of a history file (None if unreadable)."""
    from ..io import thfile
    try:
        if header > 0:
            return len(thfile.read_soil_history(path, 0, header))
        acc, _ = thfile.read_history(path, fopt=fopt, rec1=1, rec2=0)
        return len(acc)
    except (OSError, ValueError):
        return None


def count_xy_rows(path: Path) -> Optional[int]:
    """Number of data rows (lines with at least two numbers) of a 2-column spectrum file."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    n = 0
    for ln in text.splitlines():
        vals = []
        for tok in ln.replace(",", " ").split():
            try:
                vals.append(float(tok.replace("D", "E").replace("d", "e")))
            except ValueError:
                break
        if len(vals) >= 2:
            n += 1
    return n


# ======================================================================================
# Resolved analysis data shared by CHECK and AFWRITE
# ======================================================================================
NFFT_MAX = {"EQUAKE": 32768, "SOIL": 32768, "SITE": 32768, "HOUSE": 32768, "FORCE": 32768, "POINT": 32768,
            "ANALYS": 32768, "MOTION": 65536, "RELDISP": 65536, "STRESS": 65536}


class Resolved:
    """The analysis data of a model, resolved once for CHECK and AFWRITE.

    Holds the typed option records, the shared variables in their single storage location
    (spec 07 section 3), the NFFT written to the decks (nearest power of 2, D-CNV-10), the resolved
    frequency step and the frequency numbers of the selected set.
    """

    def __init__(self, model, dirs: Sequence[Union[str, Path]] = ()):
        self.m = model
        self.dirs = [Path(d) for d in dirs if d]
        R = self.rec = {name: get_record(model, name) for name in (
            "EQUAKE", "SOIL", "SOILX", "SITE", "SITEX", "POINT", "HOUSE", "HOUSEX", "INCOH", "WPASS", "FORCE",
            "ANALYS", "ANALYSX", "MOTION", "MOTIONX", "STRESS", "STRESSX", "RELD", "RELDX", "CMODFORM", "AOPT")}
        self.site, self.house, self.motion, self.analys = R["SITE"], R["HOUSE"], R["MOTION"], R["ANALYS"]
        self.limits = eduopt(model, "LIMITS")
        self.gravity = float(self.house.gravity)
        self.gelev = float(self.house.gelev)
        self.delt = float(self.site.delt)
        self.nft_raw = int(self.site.nft)
        self.fstep = float(self.site.fstep)
        rnd = eduopt(model, "NFFTROUND")
        if self.nft_raw > 0 and not is_power_of_two(self.nft_raw):
            if rnd == "UP":
                self.nft = 1 << (int(self.nft_raw - 1).bit_length())
            else:
                self.nft = nearest_power_of_two(self.nft_raw)
        else:
            self.nft = self.nft_raw
        self.df: Optional[float]
        if self.fstep > 0:
            self.df = self.fstep
        elif self.delt > 0 and self.nft > 0:
            self.df = 1.0 / (self.delt * self.nft)
        else:
            self.df = None
        self.freq_set = int(self.site.freq)
        raw = model.freq_sets.get(self.freq_set)
        self.fnums: Optional[List[int]] = None if raw is None else sorted(int(v) for v in raw)
        self.cmodform = int(R["CMODFORM"].form)
        self.type = int(self.analys.type)
        self.coh, self.wpass, self.me = int(self.house.coh), int(self.house.wpass), int(self.house.me)
        self.thfile = model.options.string("THFILE")
        self.thtit = model.options.string("THTIT")

    # ---------------------------------------------------------------- helpers
    def file(self, name: str) -> Optional[Path]:
        return resolve_file(name, self.dirs)

    def deck_file_name(self, name: str) -> str:
        """File name as written into a deck: relative names found outside the model directory
        are made absolute so the module (run in the model directory) finds them."""
        if not name:
            return ""
        p = self.file(name)
        if p is None:
            return name
        if Path(name).is_absolute() or (self.dirs and p.parent.resolve() == self.dirs[0].resolve()):
            return name
        return str(p.resolve())

    @property
    def fcut(self) -> Optional[float]:
        """Highest SSI frequency (Hz) of the selected set."""
        if not self.fnums or self.df is None:
            return None
        return max(self.fnums) * self.df

    def history_records(self) -> Optional[int]:
        """Number of records of the THFILE history (read once, cached)."""
        if not hasattr(self, "_nrec"):
            p = self.file(self.thfile)
            self._nrec = None if p is None else count_records(p, int(self.motion.fopt))
        return self._nrec

    def records_used(self) -> Optional[int]:
        """Number of control-motion records used by MOTION/STRESS/RELDISP (rec1..rec2 of THFILE)."""
        n = self.history_records()
        if n is None:
            return None
        r1 = max(1, int(self.motion.rec1))
        r2 = int(self.motion.rec2) or n
        return max(0, min(r2, n) - r1 + 1)


# ======================================================================================
# The checker
# ======================================================================================
class Checker:
    """Runs the catalogue for the enabled modules and fills a :class:`CheckReport`."""

    def __init__(self, model, modules: Optional[Sequence[str]] = None, dirs: Sequence[Union[str, Path]] = ()):
        self.m = model
        self.r = Resolved(model, dirs)
        self.enabled = [m for m in RUN_ORDER if m in (modules if modules is not None else self.r.rec["AOPT"].enabled())]
        self.rep = CheckReport(model=model.name, modules=list(self.enabled))
        self._view: Optional[ModelView] = None

    @property
    def view(self) -> ModelView:
        if self._view is None:
            self._view = ModelView(self.m)
        return self._view

    def on(self, module: str) -> bool:
        return module in self.enabled

    # ---------------------------------------------------------------- driver
    def run(self) -> CheckReport:
        rep = self.rep
        model_level = any(self.on(x) for x in MODEL_GATED)
        if model_level:
            rep.sections.append(MODEL)
            self.check_model()
            if rep.fatal:
                return rep
        for mod in self.enabled:
            if mod not in rep.sections:
                rep.sections.append(mod)
            fn = getattr(self, f"check_{mod.lower()}", None)
            if fn is not None:
                fn(mod)
        return rep

    # ================================================================== shared variables
    def _record_problems(self, module: str, rec, shared: bool = True) -> None:
        """Field-local rules of an option record; the title placeholders (``<i>`` of Error 53 ...)
        come from the rule that detected the problem (:class:`sassi.prep.options.Problem`)."""
        for p in rec.problems(shared=shared):
            kind, num, det = p
            self.rep.add(module, kind, num, det, **problem_fmt(p))

    def _gravity(self, module: str) -> None:
        if self.r.gravity <= 0:
            self.rep.error(module, 1, f"HOUSE <gravity> = {fmt_num(self.r.gravity)}")

    def _time_grid(self, module: str, harmonic_ok: bool) -> None:
        r = self.r
        for kind, num, det in time_grid_problems(r.fstep, r.delt, r.nft_raw, harmonic_ok):
            if num == 9:
                det = f"NFFT = {r.nft_raw}; {r.nft} is written (nearest power of 2, D-CNV-10)"
            self.rep.add(module, kind, num, det)
        if r.nft_raw > 0 and r.nft > NFFT_MAX.get(module, 65536) and r.limits != "UNLIMITED":
            self.rep.edu(module, "EDU-02", f"NFFT {r.nft} > {NFFT_MAX.get(module)} (D-CNV-11)")

    def _freq_set(self, module: str, full: bool = True) -> None:
        """Errors 44 / 120 and (``full``) EDU-23 duplicates, EDU-24 above Nyquist, EDU-02 > 500."""
        r = self.r
        s = r.freq_set
        if r.fnums is None:
            self.rep.error(module, 44, f"SITE <freq> = {s}", s=s)
            return
        if not r.fnums:
            self.rep.error(module, 120, i=s)
            return
        if not full:
            return
        raw = list(self.m.freq_sets.get(s, []))
        dups = sorted({v for v in raw if raw.count(v) > 1})
        if dups:
            self.rep.edu(module, "EDU-23", f"{dups[:10]} (SITE stops on duplicates, G-01)", s=s)
        if r.fstep <= 0 and r.nft > 0:
            over = [v for v in r.fnums if v > r.nft // 2]
            if over:
                self.rep.edu(module, "EDU-24", f"{len(over)} numbers > {r.nft // 2}, e.g. {over[0]}", s=s)
        if len(r.fnums) > 500 and r.limits != "UNLIMITED":
            self.rep.edu(module, "EDU-02", f"{len(r.fnums)} SSI frequencies (limit 500)")

    def _history(self, module: str) -> None:
        """Control-motion data shared by MOTION, STRESS and RELDISP (Errors 73-78, EDU-03)."""
        r = self.r
        mo = r.motion
        p = r.file(r.thfile)
        if p is None:
            self.rep.error(module, 73, f"THFILE '{r.thfile}'" if r.thfile else "THFILE not given")
        for kind, num, det in history_problems(mo.mult, mo.get("max"), mo.rec1, mo.rec2):
            self.rep.add(module, kind, num, det)
        if p is not None:
            n = r.history_records()
            if n is not None:
                if mo.rec1 > n:
                    self.rep.error(module, 74, f"<rec1> = {mo.rec1} > {n} records in {p.name}")
                used = r.records_used()
                if used is not None and r.nft > 0 and r.nft < used:
                    self._edu03(module, used)

    def _edu03(self, module: str, used: int) -> None:
        r = self.r
        nxt = 1 << (int(used - 1).bit_length())
        if r.nft_raw != r.nft:
            self.rep.edu(module, "EDU-03", f"NFFT {r.nft_raw} rounded to {r.nft} < {used} records; use {nxt}")
        else:
            self.rep.edu(module, "EDU-03", f"NFFT {r.nft} < {used} records; use {nxt} (D-CNV-10)")

    def _layer(self, module: str, l: int, halfspace: bool = False, nu_check: bool = True) -> bool:
        """Errors 19-25 (and EDU-04, EDU-11) of soil layer ``l``; False when it is undefined."""
        L = self.m.layers.get(l)
        if L is None:
            self.rep.error(module, 19, l=l)
            return False
        if not halfspace and L.thick <= 0:
            self.rep.error(module, 20, f"thickness {fmt_num(L.thick)}", l=l)
        if L.weight < 0:
            self.rep.error(module, 21, f"weight {fmt_num(L.weight)}", l=l)
        if L.vp < 0:
            self.rep.error(module, 22, f"Vp {fmt_num(L.vp)}", l=l)
        if L.vs < 0:
            self.rep.error(module, 23, f"Vs {fmt_num(L.vs)}", l=l)
        if L.pdamp < 0:
            self.rep.error(module, 24, f"Dp {fmt_num(L.pdamp)}", l=l)
        if L.sdamp < 0:
            self.rep.error(module, 25, f"Ds {fmt_num(L.sdamp)}", l=l)
        if L.pdamp >= 0.5 or L.sdamp >= 0.5:
            self.rep.edu(module, "EDU-04", f"soil layer {l}: Dp {fmt_num(L.pdamp)}, Ds {fmt_num(L.sdamp)}")
        if nu_check and L.vp > 0 and L.vs > 0 and L.vp ** 2 != L.vs ** 2:
            nu = (L.vp ** 2 - 2 * L.vs ** 2) / (2 * (L.vp ** 2 - L.vs ** 2))
            if nu > 0.47:
                self.rep.edu(module, "EDU-11", f"nu = {nu:.3f} (G-09)", l=l)
        return True

    # ================================================================== MODEL
    def check_model(self) -> None:
        m, rep = self.m, self.rep
        M = MODEL
        if not m.nodes:
            rep.error(M, 42)
            rep.fatal = True
            return
        if not m.groups:
            rep.error(M, 43)
            rep.fatal = True
            return
        v = self.view
        tol = v.tol
        # ---- nodes: gaps, unused, size
        ids = sorted(m.nodes)
        maxn = ids[-1]
        idset = set(ids)
        for n in range(1, maxn):
            if n not in idset:
                rep.warning(M, 1, n=n)
        used = v.used_dof_nodes()
        for n in ids:
            if n not in used and n not in v.k_nodes:
                rep.warning(M, 4, n=n)
        if maxn > 99999 and self.r.limits != "UNLIMITED":
            rep.edu(M, "EDU-02", f"largest node number {maxn} > 99,999 (IKTR9; 6-digit file names)")
        # ---- groups and elements
        by_group: Dict[int, List[ElemRef]] = defaultdict(list)
        for r in v.elems:
            by_group[r.group].append(r)
        for gid in sorted(m.groups):
            g = m.groups[gid]
            if not g.elements:
                rep.warning(M, 7, g=gid)
                continue
            mx = max(g.elements)
            for e in range(1, mx):
                if e not in g.elements:
                    rep.error(M, 6, e=e, g=gid)
        self._element_checks(by_group)
        self._tshell_checks(by_group)
        self._property_checks()
        self._mass_checks()
        self._interaction_checks()
        self._dimension_check()

    # ---------------------------------------------------------------- TSHELL (P2, D-ELM-08)
    def _tshell_checks(self, by_group: Dict[int, List[ElemRef]]) -> None:
        """TSHELL element data that HOUSE would reject (EDU-26, value not allowed): EINT must be 0
        (reduced, default) or 1 (selective) (manual 9.4.9) and the thickness THICK must be > 0 (spec 08
        4.6).  Listed per group with the first offending elements."""
        rep = self.rep
        for gid in sorted(by_group):
            refs = [r for r in by_group[gid] if r.type == TSHELL]
            if not refs:
                continue
            bad_eint = [r.id for r in refs if int(r.elem.eint) not in (0, 1)]
            bad_thick = [r.id for r in refs if not float(r.elem.thick) > 0.0]
            if bad_eint:
                rep.edu(MODEL, "EDU-26", f"TSHELL group {gid}: EINT must be 0 (reduced) or 1 (selective) for "
                                         f"{len(bad_eint)} elements, e.g. {bad_eint[:8]}")
            if bad_thick:
                rep.edu(MODEL, "EDU-26", f"TSHELL group {gid}: thickness (THICK) must be > 0 for {len(bad_thick)} "
                                         f"elements, e.g. {bad_thick[:8]}")

    # ---------------------------------------------------------------- elements
    def _element_checks(self, by_group: Dict[int, List[ElemRef]]) -> None:
        """Errors 7-12 and 41, Warnings 2-3 and EDU-05 (D-CHK-05 thresholds).

        Regular SOLIDs (8 distinct nodes) and quadrilaterals (4 distinct nodes) are checked in
        vectorised batches; degenerate shapes (repeated nodes: prisms, triangles entered as quads)
        go through the scalar path with consecutive duplicates removed.
        """
        rep, v = self.rep, self.view
        M = MODEL
        tol = v.tol
        P = v.P
        solids: List[ElemRef] = []
        planes: List[ElemRef] = []
        quads: List[Tuple[ElemRef, List[int]]] = []
        for gid in sorted(by_group):
            for r in by_group[gid]:
                e = r.elem
                ns = [n for n in e.nodes]
                bad = [n for n in ns if n and n not in P]
                for n in dict.fromkeys(bad):
                    rep.error(M, 41, n=n, g=gid, e=e.id)
                given = sum(1 for n in ns if n)
                if given < MIN_ELEMENT_NODES.get(r.type, 2):
                    rep.error(M, 7, f"{given} nodes, {ELEMENT_TYPE_NAMES.get(r.type)} needs "
                                    f"{MIN_ELEMENT_NODES.get(r.type)}", e=e.id, g=gid)
                    continue
                if bad:
                    continue
                t = r.type
                if t == BEAMS:
                    if any(e.ki[k] and e.kj[k] for k in range(6)):
                        rep.error(M, 10, e=e.id, g=gid)
                if t in (BEAMS, GENERAL):
                    I, J = ns[0], ns[1]
                    if I == J or _norm(P[I] - P[J]) <= tol:
                        rep.error(M, 8, e=e.id, g=gid)
                        continue
                    if len(ns) >= 3 and ns[2]:
                        a, b = P[J] - P[I], P[ns[2]] - P[I]
                        if _norm(np.cross(a, b)) <= 1e-9 * _norm(a) * max(_norm(b), 1e-300) or _norm(b) <= tol:
                            rep.error(M, 9, "K node collinear with I-J", e=e.id, g=gid)
                elif t in (SHELL, TSHELL, PLANE):
                    ids = [n for n in ns if n]
                    if len(ids) == 4 and len(set(ids)) == 4:
                        quads.append((r, ids))
                    else:
                        self._polygon_scalar(r, ids)
                elif t == SOLID:
                    if len(ns) == 8 and len(set(ns)) == 8:
                        solids.append(r)
                    else:
                        self._solid_faces_scalar(r, list(ns) + [0] * (8 - len(ns)))
        if quads:
            self._quad_batch(quads)
        # ---- solids: face angles (W2) and the Jacobian at the Gauss points (EDU-05)
        if solids:
            X = np.array([[P[n] for n in r.elem.nodes] for r in solids])
            fidx = np.array(SOLID_FACES) - 1
            F = X[:, fidx]                                            # (n, 6, 4, 3)
            edges = np.linalg.norm(F - np.roll(F, -1, axis=2), axis=3)
            ang = quad_angles(F.reshape(-1, 4, 3)).reshape(len(solids), 6, 4)
            for k, fi in zip(*np.nonzero(np.any(edges <= tol, axis=2))):
                self._face_scalar(solids[k], [solids[k].elem.nodes[j] for j in fidx[fi]], fi + 1)
            okf = ~np.any(edges <= tol, axis=2)
            badf = okf & ((ang.min(axis=2) < 15.0) | (ang.max(axis=2) > 165.0))
            for k, fi in zip(*np.nonzero(badf)):
                rep.warning(M, 2, f"angles {ang[k, fi].min():.1f}..{ang[k, fi].max():.1f} deg", e=solids[k].id,
                            g=solids[k].group, f=int(fi) + 1)
            dets = hex_jacobians(X)
            scale = np.maximum(np.max(np.abs(dets), axis=1), 1e-300)
            for k in np.flatnonzero(np.any(dets <= 1e-10 * scale[:, None], axis=1)):
                rep.edu(M, "EDU-05", f"min det J {dets[k].min():.3g}", e=solids[k].id, g=solids[k].group)
        planes = [r for r, _ in quads if r.type == PLANE and not getattr(r, "_degenerate", False)]
        if planes:
            X = np.array([[P[n] for n in r.elem.nodes[:4]] for r in planes])
            dets = quad_jacobians_xz(X)
            sg = np.sign(dets)
            for k in np.flatnonzero(~(np.all(sg > 0, axis=1) | np.all(sg < 0, axis=1))):
                rep.edu(M, "EDU-05", "det J changes sign or vanishes (X-Z plane)", e=planes[k].id,
                        g=planes[k].group)

    def _pts(self, r: ElemRef, ids: Sequence[int]) -> Dict[int, np.ndarray]:
        """Node coordinates of a polygon element (PLANE elements projected on the X-Z plane)."""
        P = self.view.P
        if r.type == PLANE:
            return {n: P[n] * np.array([1.0, 0.0, 1.0]) for n in ids}
        return {n: P[n] for n in ids}

    def _polygon_scalar(self, r: ElemRef, ids: List[int]) -> None:
        """Errors 9, 11, 12 and Warnings 2, 3 of one SHELL/TSHELL/PLANE element (any node pattern)."""
        rep, tol = self.rep, self.view.tol
        M, e, g = MODEL, r.id, r.group
        area_tol = tol * tol
        uid, upts = dedupe_polygon(ids, self._pts(r, ids), tol)
        r._degenerate = len(uid) != 4                                          # type: ignore[attr-defined]
        if len(upts) < 3:
            rep.error(M, 9, e=e, g=g)
            return
        if len(upts) == 3:
            if _norm(np.cross(upts[1] - upts[0], upts[2] - upts[0])) / 2 <= area_tol:
                rep.error(M, 9, e=e, g=g)
                return
        else:
            tri = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
            amax = max(_norm(np.cross(upts[b] - upts[a], upts[c] - upts[a])) / 2 for a, b, c in tri)
            if amax <= area_tol:
                rep.error(M, 9, e=e, g=g)
                return
            cv = np.cross(upts[2] - upts[0], upts[3] - upts[1])
            A = _norm(cv) / 2
            if A <= area_tol:
                rep.error(M, 11, e=e, g=g)
                return
            w = abs(float(np.dot(cv / _norm(cv), upts[3] - upts[0]))) / math.sqrt(A)
            if w > 5e-2:
                rep.error(M, 12, f"warp {w:.3g} > 5e-2 (D-CHK-05)", e=e, g=g)
                return
            if w > 1e-3:
                rep.warning(M, 3, f"warp {w:.3g} > 1e-3", e=e, g=g)
        ang = polygon_angles(upts)
        if min(ang) < 15.0 or max(ang) > 165.0:
            rep.warning(M, 2, f"angles {min(ang):.1f}..{max(ang):.1f} deg (D-CHK-05)", e=e, g=g, f=1)

    def _quad_batch(self, quads: List[Tuple[ElemRef, List[int]]]) -> None:
        """Vectorised Errors 9, 11, 12 and Warnings 2, 3 of quadrilaterals with 4 distinct nodes."""
        rep, tol = self.rep, self.view.tol
        P = self.view.P
        area_tol = tol * tol
        proj = np.array([[1.0, 0.0, 1.0] if r.type == PLANE else [1.0, 1.0, 1.0] for r, _ in quads])
        Q = np.array([[P[n] for n in ids] for _, ids in quads]) * proj[:, None, :]
        edges = np.linalg.norm(Q - np.roll(Q, -1, axis=1), axis=2)
        tri = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
        amax = np.max([np.linalg.norm(np.cross(Q[:, b] - Q[:, a], Q[:, c] - Q[:, a]), axis=1) / 2
                       for a, b, c in tri], axis=0)
        cv = np.cross(Q[:, 2] - Q[:, 0], Q[:, 3] - Q[:, 1])
        lc = np.linalg.norm(cv, axis=1)
        A = lc / 2
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.abs(np.einsum("ij,ij->i", cv / lc[:, None], Q[:, 3] - Q[:, 0])) / np.sqrt(A)
        ang = quad_angles(Q)
        for k, (r, ids) in enumerate(quads):
            if np.any(edges[k] <= tol):                       # coincident nodes: scalar path
                self._polygon_scalar(r, ids)
                continue
            e, g = r.id, r.group
            if amax[k] <= area_tol:
                rep.error(MODEL, 9, e=e, g=g)
            elif A[k] <= area_tol:
                rep.error(MODEL, 11, e=e, g=g)
            elif w[k] > 5e-2:
                rep.error(MODEL, 12, f"warp {w[k]:.3g} > 5e-2 (D-CHK-05)", e=e, g=g)
            else:
                if w[k] > 1e-3:
                    rep.warning(MODEL, 3, f"warp {w[k]:.3g} > 1e-3", e=e, g=g)
                if ang[k].min() < 15.0 or ang[k].max() > 165.0:
                    rep.warning(MODEL, 2, f"angles {ang[k].min():.1f}..{ang[k].max():.1f} deg (D-CHK-05)",
                                e=e, g=g, f=1)

    def _face_scalar(self, r: ElemRef, fn: Sequence[int], fi: int) -> None:
        uid, upts = dedupe_polygon(list(fn), self.view.P, self.view.tol)
        if len(upts) < 3:
            return
        ang = polygon_angles(upts)
        if min(ang) < 15.0 or max(ang) > 165.0:
            self.rep.warning(MODEL, 2, f"angles {min(ang):.1f}..{max(ang):.1f} deg", e=r.id, g=r.group, f=fi)

    def _solid_faces_scalar(self, r: ElemRef, pad: List[int]) -> None:
        """Warning 2 on every face of a degenerate SOLID (prism, pyramid: repeated nodes)."""
        for fi, f in enumerate(SOLID_FACES, start=1):
            self._face_scalar(r, [pad[k - 1] for k in f], fi)

    # ---------------------------------------------------------------- properties
    def _property_checks(self) -> None:
        m, rep, v = self.m, self.rep, self.view
        M = MODEL
        g = self.r.gravity if self.r.gravity > 0 else 1.0
        used_m: Dict[int, Set[int]] = defaultdict(set)      # material -> element types using it
        used_l: Set[int] = set()
        used_r: Set[int] = set()
        used_sc: Set[int] = set()
        for r in v.elems:
            e = r.elem
            if r.type in (SOLID, PLANE):
                if r.etype == 2:
                    used_l.add(e.mat)
                else:
                    used_m[e.mat].add(r.type)
            elif r.type in (SHELL, TSHELL, BEAMS):
                used_m[e.mat].add(r.type)
                if r.type == BEAMS:
                    used_r.add(e.prop)
            elif r.type == SPRING:
                used_sc.add(e.prop)
            elif r.type == GENERAL:
                if e.prop not in m.matrices:
                    rep.error(M, 83, p=e.prop)
        for mid in sorted(used_m):
            mat = m.materials.get(mid)
            if mat is None:
                rep.error(M, 13, m=mid)
                continue
            try:
                c = elastic_constants(mat.mtype, mat.val1, mat.val2, mat.weight, g)
                E, nu = c.E, c.nu
            except ValueError as exc:
                if mat.mtype not in (1, 2, 3):
                    rep.error(M, 14, f"material type {mat.mtype}", m=mid)
                else:
                    rep.error(M, 15, str(exc), m=mid)
                continue
            if not E > 0:
                rep.error(M, 14, f"E = {fmt_num(E)}", m=mid)
            if not nu > 0 or nu >= 0.5:
                rep.error(M, 15, f"nu = {nu:.6g}", m=mid)
            if mat.weight < 0:
                rep.error(M, 16, f"weight {fmt_num(mat.weight)}", m=mid)
            if mat.pdamp < 0:
                rep.error(M, 17, f"Dp {fmt_num(mat.pdamp)}", m=mid)
            if mat.sdamp < 0:
                rep.error(M, 18, f"Ds {fmt_num(mat.sdamp)}", m=mid)
            if mat.pdamp >= 0.5 or mat.sdamp >= 0.5:
                rep.edu(M, "EDU-04", f"material {mid}: Dp {fmt_num(mat.pdamp)}, Ds {fmt_num(mat.sdamp)}")
            if used_m[mid] & {SHELL, TSHELL} and mat.pdamp != mat.sdamp:
                rep.edu(M, "EDU-27", f"Dp {fmt_num(mat.pdamp)} != Ds {fmt_num(mat.sdamp)}: shells use "
                                     f"E* = E c(Ds) (requirements 4.0.2)", m=mid)
        for lid in sorted(used_l):
            self._layer(M, lid)
        for pid in sorted(used_r):
            s = m.sections.get(pid)
            if s is None:
                rep.error(M, 26, p=pid)
                continue
            for num, val, bad in ((27, s.axial, s.axial <= 0), (28, s.shear2, s.shear2 < 0),
                                  (29, s.shear3, s.shear3 < 0), (30, s.tors, s.tors <= 0),
                                  (31, s.flex2, s.flex2 < 0), (32, s.flex3, s.flex3 < 0)):
                if bad:
                    rep.error(M, num, f"value {fmt_num(val)}", p=pid)
        for pid in sorted(used_sc):
            s = m.springs.get(pid)
            if s is None:
                rep.error(M, 33, p=pid)
                continue
            for num, val in zip(range(34, 40), s.k):
                if val < 0:
                    rep.error(M, num, f"value {fmt_num(val)}", p=pid)
            if s.damp >= 0.5:
                rep.edu(M, "EDU-04", f"spring property {pid}: damping {fmt_num(s.damp)}")

    # ---------------------------------------------------------------- masses
    def _mass_checks(self) -> None:
        m, rep = self.m, self.rep
        M = MODEL
        maxn = max(m.nodes)
        out = sorted(n for n in set(m.tmass) | set(m.rmass) if n < 1 or n > maxn)
        if out:
            rep.error(M, 40, f"nodes {out[:10]} outside 1..{maxn}")
        for n, vals in sorted(m.tmass.items()):
            nd = m.nodes.get(n)
            if nd is not None and any(vals[k] != 0 and nd.fix[k] for k in range(3)):
                rep.warning(M, 5, n=n)
        for n, vals in sorted(m.rmass.items()):
            nd = m.nodes.get(n)
            if nd is not None and any(vals[k] != 0 and nd.fix[3 + k] for k in range(3)):
                rep.warning(M, 6, n=n)

    # ---------------------------------------------------------------- interaction nodes
    def _interaction_checks(self) -> None:
        m, rep, v = self.m, self.rep, self.view
        M = MODEL
        inter = v.interaction_nodes()
        for n in inter:
            if all(m.nodes[n].fix[:3]):
                rep.error(M, 124, "all translations fixed (D-CHK-06)", i=n)
        exc = v.excavation()
        # EXCSTRCHK (G-15): interior excavation nodes shared with the structure.  Exception (non-linear
        # soil SSI, requirements 4.4 item 7): the near-field soil SOLID/PLANE elements of a non-linear soil
        # group (PINGRP, HOUSEX <nlssi> = 1) replace the excavated soil at its own nodes (manual rule 4); at
        # interaction nodes (FV) such a node is not reported.  Other nodes: severity of excstrchk_severity.
        iset = set(inter)
        nl_groups = ({int(k) for k, _ in m.options.entries("PINGRP")}
                     if int(self.r.rec["HOUSEX"].nlssi) == 1 else set())
        exc_nodes = set(exc["nodes"])
        for n, owners in sorted(excstrchk(v).items()):
            if nl_groups and n in iset and n not in v.k_nodes and all(
                    r.group in nl_groups and r.type in (SOLID, PLANE) for r in v.node_elems.get(n, [])
                    if r.type in (BEAMS, SPRING, GENERAL) or (r.etype != 2 and r.type in (SOLID, PLANE, SHELL, TSHELL))):
                continue
            # severity: the rule HOUSE applies (sassi.core.house_lib.excstrchk_kind; final audit)
            rep.edu(M, "EXCSTRCHK", "; ".join(owners[:3]), kind=excstrchk_severity(v, n, iset, exc_nodes), n=n)
        # EDU-22: excavation boundary nodes that are not interaction nodes
        if exc["elements"]:
            imp = int(self.r.house.imp)
            need = set(exc["fsin"])
            if imp == 0:
                need |= exc["nodes"]
            missing = sorted(need - iset)
            for n in missing:
                rep.edu(M, "EDU-22", "FV needs every excavation node" if imp == 0 and n not in exc["fsin"]
                        else "lateral/bottom excavation surface (FSIN)", n=n)
        # EDU-01: interaction nodes on soil-layer interfaces 1..POINT<layer>+1, at or below grade
        depths = v.interfaces()
        if inter and depths is not None:
            itol = v.interface_tol(depths)
            last = max(int(get_record(m, "POINT").layer), 0) + 1
            for n in inter:
                p = v.P[n]
                if not np.all(np.isfinite(p)):
                    continue
                if p[2] > v.gelev + itol:
                    rep.edu(M, "EDU-01", f"z = {fmt_num(float(p[2]))} above the ground elevation "
                                         f"{fmt_num(v.gelev)} (G-12)", n=n)
                    continue
                k, mis = v.interface_of(float(p[2]), depths)
                if mis > itol:
                    rep.edu(M, "EDU-01", f"depth {fmt_num(float(v.gelev - p[2]))}, nearest interface {k} "
                                         f"at {fmt_num(float(depths[k - 1]))}", n=n)
                elif k > last:
                    rep.edu(M, "EDU-01", f"on interface {k}, below the POINT near-field zone "
                                         f"(interfaces 1..{last}; increase POINT <layer>)", n=n)
        # EDU-21: bottom-up ascending numbering (error when incoherent)
        if len(inter) > 1:
            z = np.array([v.P[n][2] for n in inter])
            bad = np.flatnonzero(np.diff(z) < -v.tol)
            if bad.size:
                k = int(bad[0])
                kind = "Error" if self.r.coh == 1 else EDU21_COHERENT_KIND
                rep.edu(M, "EDU-21", f"node {inter[k + 1]} (z = {fmt_num(float(z[k + 1]))}) follows node "
                                     f"{inter[k]} (z = {fmt_num(float(z[k]))}); {bad.size} inversions (G-13)",
                        kind=kind)
        # EDU-06: unrestrained shell drilling rotation
        for n, note in drilling_unrestrained(v).items():
            rep.edu(M, "EDU-06", note, n=n)
        # EDU-08 / EDU-10: excavated elements vs the TOPL layer at their depth
        if exc["elements"] and depths is not None:
            fcut = self.r.fcut
            for r in exc["elements"]:
                zs = np.array([v.P[n][2] for n in dict.fromkeys(x for x in r.elem.nodes if x)])
                dmid = v.gelev - float(zs.mean())
                k = int(np.searchsorted(depths, dmid, side="right"))      # layer index 1..nTOPL (+1 = half-space)
                expect = m.topl[k - 1] if 1 <= k <= len(m.topl) else int(self.r.site.hs)
                if r.elem.mat != expect:
                    rep.edu(M, "EDU-08", f"MSET {r.elem.mat}, TOPL layer {expect} at depth {dmid:.4g}",
                            e=r.id, g=r.group)
                L = m.layers.get(r.elem.mat)
                if fcut and L is not None and L.vs > 0:
                    hz = float(zs.max() - zs.min())
                    if hz > L.vs / (5.0 * fcut) * (1 + 1e-9):
                        rep.edu(M, "EDU-10", f"element {r.id} group {r.group}: height {hz:.4g} > Vs/(5 f_cut) = "
                                             f"{L.vs / (5 * fcut):.4g} (G-06)")
        # EDU-28: an L number repeated among the embedment layers of TOPL.  The manual allows a repeated
        # L number for identical layers, "but not for embedment layers" (spec 07 L, spec 05a TOPL
        # warning): each embedment layer needs its own L number, used as MSET by the excavated
        # elements of that layer (R5/R6), so that layer-dependent data (SITEX,1 strain-compatible
        # properties, STRESS soil pressures layer by layer) stay attached to the right layer.
        nemb = max(int(get_record(m, "POINT").layer), 0)
        if exc["elements"] and nemb > 1:
            emb = [int(l) for l in m.topl[:nemb]]
            dup = sorted({l for l in emb if emb.count(l) > 1})
            if dup:
                rep.edu(M, "EDU-28", f"L {', '.join(str(l) for l in dup)} repeated in the {len(emb)} embedment "
                                     f"layers of TOPL (POINT <layer> = {nemb}); give each embedment layer its own "
                                     "L number (spec 07 L, spec 05a TOPL)")

    def _dimension_check(self) -> None:
        """EDU-26 (D-PNT-01): HOUSE <dim> consistent with the element types."""
        types = {r.type for r in self.view.elems}
        dim = int(self.r.house.dim)
        if dim == 2 and PLANE in types:
            self.rep.edu(MODEL, "EDU-26", "PLANE elements in a 3D model (HOUSE <dim> = 2)")
        if dim == 1 and types & {SOLID, SHELL, TSHELL}:
            self.rep.edu(MODEL, "EDU-26", "3D elements in a 2D model (HOUSE <dim> = 1)")

    # ================================================================== EQUAKE
    def check_equake(self, mod: str) -> None:
        rep, m, r = self.rep, self.m, self.r
        eq = r.rec["EQUAKE"]
        self._gravity(mod)
        if r.delt <= 0:
            rep.error(mod, 49, f"<delt> = {fmt_num(r.delt)}")
        self._record_problems(mod, eq)
        files = {name: dict(get_entries(m, name)) for name in ("RSIN", "RSOUT", "ACCIN", "ACCOUT", "TPSD")}
        accopt = int(eq.accopt)
        comps = sorted({k for d in files.values() for k in d if isinstance(k, int)})
        if accopt != 2:
            rsin = {k: rec for k, rec in files["RSIN"].items() if rec.file}
            if not rsin:
                rep.error(mod, 84)
            nrfreq = int(eq.nrfreq)
            for i, rec in sorted(rsin.items()):
                p = r.file(rec.file)
                if p is None:
                    rep.error(mod, 85, rec.file, i=i)
                    continue
                nrec = count_xy_rows(p)
                if not eq.given(2):
                    if i == min(rsin):
                        nrfreq = nrec or 0
                        if nrfreq <= 0:
                            rep.error(mod, 91, f"RSIN {i} has no records")
                if nrec is not None and nrfreq > 0 and nrec != nrfreq:
                    rep.error(mod, 89, f"{nrec} records, <nrfreq> = {nrfreq}", i=i)
            for i in sorted(rsin):
                if not _file_of(files["RSOUT"], i):
                    rep.error(mod, 86, i=i)
                if not _file_of(files["ACCOUT"], i):
                    rep.error(mod, 87, i=i)
                if accopt == 1 and not _file_of(files["ACCIN"], i):
                    rep.error(mod, 88, i=i)
            if accopt == 1 and not rsin:
                for i in comps:
                    if not _file_of(files["ACCIN"], i):
                        rep.error(mod, 88, i=i)
        else:
            active = sorted(set(files["ACCIN"]) | set(files["RSOUT"]))
            if not active:
                rep.error(mod, 88, "External Accel selected: no ACCIN file", i=1)
            for i in active:
                if not _file_of(files["ACCIN"], i):
                    rep.error(mod, 88, i=i)
                elif r.file(_file_of(files["ACCIN"], i)) is None:
                    rep.error(mod, 88, f"{_file_of(files['ACCIN'], i)} not found", i=i)
                if not _file_of(files["RSOUT"], i):
                    rep.error(mod, 86, i=i)
        if int(eq.corr):
            pairs = get_entries(m, "CORR")
            if not pairs:
                rep.error(mod, 93)
            for _, c in pairs:
                self._record_problems(mod, c)
        if r.delt > 0 and eq.dur / r.delt > 32768 and r.limits != "UNLIMITED":
            rep.edu(mod, "EDU-02", f"{int(eq.dur / r.delt)} time steps > 32,768")
        if eq.dur < 20.0:
            rep.edu(mod, "EDU-28", f"total duration {fmt_num(eq.dur)} s < 20 s (SRP 3.7.1)")
        if int(eq.tpsd):
            for i in comps:
                f = _file_of(files["TPSD"], i)
                if f and r.file(f) is None:
                    rep.edu(mod, "EDU-28", f"target PSD file {f} not found (spectrum {i})")

    # ================================================================== SOIL
    def soil_profile(self) -> List[Tuple[int, Any]]:
        return [(k, rec) for k, rec in get_entries(self.m, "SPRO") if isinstance(k, int)]

    def check_soil(self, mod: str) -> None:
        rep, m, r = self.rep, self.m, self.r
        soil = r.rec["SOIL"]
        sx = r.rec["SOILX"]
        self._record_problems(mod, soil)
        if soil.legacy:
            return
        self._record_problems(mod, sx)
        for kind, num, det in time_grid_problems(0.0, r.delt, r.nft_raw, harmonic_ok=False):
            rep.add(mod, kind, num, det if num != 9 else f"NFFT = {r.nft_raw}; {r.nft} is written")
        prof = self.soil_profile()
        nlay = len(prof)
        layers_no = [k for k, _ in prof]
        expect = list(range(1, nlay + 1))
        if layers_no != expect:
            missing = sorted(set(range(1, max(layers_no or [0]) + 1)) - set(layers_no))
            for i in missing:
                rep.edu(mod, "EDU-25", "SPRO sublayers must be numbered 1..N (last = half-space, D-SOL-02)", i=i)
        if nlay < 2:
            rep.error(mod, 95, "the SOIL profile needs SPRO sublayers and the half-space (last SPRO)")
        for idx, (k, rec) in enumerate(prof):
            self._layer(mod, int(rec.prop), halfspace=(idx == nlay - 1))
        labels_used = [rec.dynprop for _, rec in prof if rec.dynprop]
        if prof and not labels_used:
            rep.error(mod, 95)
        dyn = {label for (label, _), _rec in get_entries(m, "DYNP")}
        distinct = list(dict.fromkeys(labels_used))
        if len(distinct) > 15 and r.limits == "PREP":
            rep.error(mod, 96, f"{len(distinct)} > 15 (LIMITS,PREP, D-SOL-05)")
        if len(distinct) > 100 and r.limits != "UNLIMITED":
            rep.edu(mod, "EDU-02", f"{len(distinct)} dynamic properties > 100")
        for lab in distinct:
            if lab not in dyn:
                rep.error(mod, 97, p=lab)
                rep.error(mod, 99, p=lab)
                continue
            pts = [rec for (l2, _), rec in get_entries(m, "DYNP") if l2 == lab]
            gpts = [p for p in pts if p.given(2) and p.given(3)]
            dpts = [p for p in pts if p.given(4) and p.given(5)]
            if not gpts:
                rep.error(mod, 97, p=lab)
            if any(not 0.0 <= p.g <= 1.0 for p in gpts):
                rep.error(mod, 98, p=lab)
            if not dpts:
                rep.error(mod, 99, p=lab)
            if len(pts) > 11 and r.limits != "UNLIMITED":
                rep.edu(mod, "EDU-02", f"dynamic property {lab}: {len(pts)} points > 11 (D-SOL-05)")
        cl = int(sx.cl) if int(sx.cl) > 0 else int(r.site.cl)
        if nlay and not 1 <= cl <= nlay:
            rep.error(mod, 104, f"control layer {cl} (1..{nlay})")
        # output requests
        for name in ("SACC", "SRS", "SSTR", "SSAF", "SFOU"):
            for k, rec in get_entries(m, name):
                if isinstance(k, int) and nlay and not 1 <= k <= nlay:
                    rep.edu(mod, "EDU-25", f"{name} request for sublayer {k} (1..{nlay})", i=k)
                if name in ("SSAF", "SFOU"):
                    self._record_problems(mod, rec)
                if name == "SSAF" and int(rec.save) and not 1 <= int(rec.layer2) <= max(nlay, 1):
                    rep.error(mod, 109, f"<layer2> = {rec.layer2} (1..{nlay})", i=k)
        rs_req = any(int(rec.save) for _, rec in get_entries(m, "SRS")) or \
            any(int(rec.save) for _, rec in get_entries(m, "SSAF"))
        if rs_req and not m.damp:
            rep.error(mod, 107)
        # input history
        fname = sx.file or r.thfile
        p = r.file(fname)
        if p is None:
            rep.error(mod, 73, f"'{fname}'" if fname else "THFILE not given")
        elif soil.nrval > 0:
            from ..io import thfile
            try:
                nvals: Optional[int] = len(thfile.read_soil_history(p, 0, max(int(soil.header), 0)))
            except (OSError, ValueError):
                nvals = None
            if nvals is not None and nvals < soil.nrval:
                rep.error(mod, 100, f"{soil.nrval} values requested, {nvals} in {p.name}")
        if soil.nrval > 0 and r.nft > 0 and r.nft < soil.nrval:
            self._edu03(mod, int(soil.nrval))
        # nonlinear soil (P2) and the half-space consistency (D-SOL-02)
        nls = m.options.record("NLSOIL")
        if nls is not None and nls.integer(1, 0) == 1:
            # one NLSLAYER set per soil sublayer; the last SPRO entry is the (elastic) half-space (D-SOL-02)
            have = {k for k, _ in m.options.entries("NLSLAYER") if isinstance(k, int)}
            missing = [k for k in range(1, nlay) if k not in have]
            if missing:
                rep.error(mod, 125, f"no NLSLAYER set for soil sublayer(s) {missing}")
        if self.on("SITE") and prof and int(r.site.hs) and int(prof[-1][1].prop) != int(r.site.hs):
            rep.edu(mod, "EDU-28", f"SOIL half-space (last SPRO) is layer {prof[-1][1].prop}, SITE <hs> is "
                                   f"{r.site.hs} (D-SOL-02)")

    # ================================================================== SITE
    def effective_waves(self) -> List[Any]:
        """WAVE entries, or the new-model default field (vertical SV; SH for <wopt> = 1)."""
        ents = [rec for _, rec in get_entries(self.m, "WAVE")]
        if ents:
            return ents
        from .options import WaveRecord
        return [WaveRecord("WAVE", ["2" if int(self.r.site.wopt) == 0 else "4", "1", "1", "1", "0"])]

    def check_site(self, mod: str) -> None:
        rep, m, r = self.rep, self.m, self.r
        s = r.site
        self._gravity(mod)
        self._record_problems(mod, s, shared=False)
        self._time_grid(mod, harmonic_ok=True)
        self._freq_set(mod)
        if s.cm not in (0, 1, 2):
            rep.edu(mod, "EDU-26", f"SITE <cm> = {s.cm} (0 x', 1 y', 2 z')")
        if s.wopt not in (0, 1):
            rep.edu(mod, "EDU-26", f"SITE <wopt> = {s.wopt} (0 R/SV/P, 1 SH/L)")
        if 0 < s.nl < 10:
            rep.edu(mod, "EDU-28", f"{s.nl} generated half-space layers (10-20 recommended, D-SIT-03)")
        if not m.topl:
            rep.error(mod, 46)
        else:
            if len(m.topl) > 100 and r.limits == "PREP":
                rep.warning(mod, 8, f"{len(m.topl)} top layers: only the first 100 are written (LIMITS,PREP)")
            if len(m.topl) > 200 and r.limits != "UNLIMITED":
                rep.edu(mod, "EDU-02", f"{len(m.topl)} top layers > 200")
            for l in dict.fromkeys(m.topl):
                self._layer(mod, l)
            nI = len(m.topl) + 1
            if not 1 <= int(s.cl) <= nI:
                rep.error(mod, 104, f"SITE <cl> = {s.cl} (1..{nI})")
        hs = int(s.hs)
        if hs or int(s.nl) > 0:
            self._layer(mod, hs, halfspace=True)
        # passing frequency of the layers (G-05)
        fcut = r.fcut
        if fcut:
            low = [l for l in m.topl if l in m.layers and m.layers[l].thick > 0 and m.layers[l].vs > 0
                   and m.layers[l].vs / (5 * m.layers[l].thick) < fcut]
            if low:
                L = m.layers[low[0]]
                rep.edu(mod, "EDU-10", f"layer {low[0]}: Vs/(5h) = {L.vs / (5 * L.thick):.4g} Hz < f_cut "
                                       f"{fcut:.4g} Hz ({len(low)} layers, G-05)")
        # waves (Mode 2)
        if int(s.mode2):
            fam = (1, 2, 3) if int(s.wopt) == 0 else (4, 5)
            waves = [w for w in self.effective_waves() if int(w.type) in fam and int(w.opt) != 0]
            if not waves:
                rep.error(mod, 51)
            for w in waves:
                for p in w.problems():
                    kind, num, det = p
                    rep.add(mod, kind, num, det, **{"w": w.type, **problem_fmt(p)})
            if waves:
                for i, attr in ((1, "ratio1"), (2, "ratio2")):
                    tot = sum(float(w.get(attr)) for w in waves)
                    if abs(tot - 1.0) > 1e-6:
                        rep.error(mod, 55, f"sum = {tot:.6g}", i=i)
        # non-linear soil properties (EDU-07, D-SOL-12)
        if int(r.rec["SITEX"].soilmode) == 1:
            n_topl = len(m.topl)
            p88 = r.file("FILE88")
            if self.on("SOIL"):
                n_soil = len(self.soil_profile()) - 1
                if n_soil != n_topl:
                    rep.edu(mod, "EDU-07", f"SOIL has {n_soil} sublayers (without the half-space), TOPL "
                                           f"{n_topl} (D-SOL-12)")
            elif p88 is None:
                rep.edu(mod, "EDU-07", "Non-Linear Soil selected (SITEX) but FILE88 does not exist and SOIL "
                                       "is not enabled")
            else:
                n88 = _file88_layers(p88)
                if n88 is not None and n88 != n_topl:
                    rep.edu(mod, "EDU-07", f"FILE88 has {n88} layers, TOPL {n_topl} (D-SOL-12)")

    # ================================================================== POINT
    def check_point(self, mod: str) -> None:
        rep, m, r = self.rep, self.m, self.r
        pt = r.rec["POINT"]
        self._record_problems(mod, pt)
        if m.topl and int(pt.layer) > len(m.topl):
            rep.error(mod, 56, f"<layer> = {pt.layer} > {len(m.topl)} top layers")
        if int(pt.layer) > 50 and r.limits != "UNLIMITED":
            rep.edu(mod, "EDU-02", f"{pt.layer} embedment layers > 50")
        if int(r.house.dim) not in (1, 2):
            rep.edu(mod, "EDU-26", f"HOUSE <dim> = {r.house.dim}: POINT needs 1 (POINT2) or 2 (POINT3)")
        self._freq_set(mod)

    # ================================================================== HOUSE
    def check_house(self, mod: str) -> None:
        rep, m, r = self.rep, self.m, self.r
        h = r.house
        self._gravity(mod)
        v = self.view
        # symmetry planes (Errors 2-5)
        symm = get_entries(m, "SYMM")
        for no, rec in symm:
            nodes = [int(rec.get(f"node{k}")) for k in (1, 2, 3)]
            given = [n for n in nodes if n]
            if any(n < 0 or (n and n not in m.nodes) for n in nodes):
                rep.error(mod, 2, f"nodes {nodes}", i=no)
                continue
            if len(given) < 2:
                rep.error(mod, 3, i=no)
                continue
            if len(given) == 2:
                if given[0] == given[1] or _norm(v.P[given[0]] - v.P[given[1]]) <= v.tol:
                    rep.error(mod, 4, i=no)
            else:
                a, b = v.P[given[1]] - v.P[given[0]], v.P[given[2]] - v.P[given[0]]
                if len(set(given)) < 3 or _norm(np.cross(a, b)) <= 1e-9 * _norm(a) * max(_norm(b), 1e-300):
                    rep.error(mod, 5, i=no)
        if len(symm) > 2:
            rep.edu(mod, "EDU-26", f"{len(symm)} symmetry planes (at most 2)")
        self._symmetry_geometry(mod, symm)
        # options and combinations
        if int(h.dim) not in (1, 2):
            rep.edu(mod, "EDU-26", f"HOUSE <dim> = {h.dim} (1 2D, 2 3D; 1D not available)")
        if int(h.imp) not in (0, 1, 2):
            rep.edu(mod, "EDU-26", f"HOUSE <imp> = {h.imp} (0 FV, 1 FFV, 2 FI)")
        elif int(h.imp) != 0 and v.excavation()["elements"]:
            rep.edu(mod, "EDU-12", f"<imp> = {h.imp} (G-21)")
        inter = v.interaction_nodes()
        if r.coh == 1:
            inc = r.rec["INCOH"]
            wp = r.rec["WPASS"]
            if int(h.dim) != 2:
                rep.edu(mod, "EDU-26", "incoherent analysis needs a 3D model (G-17)")
            if symm:
                rep.edu(mod, "EDU-26", "incoherent analysis is not allowed with SYMM (G-17)")
            self._record_problems(mod, inc)
            if int(wp.cohf) == 1 and inc.alpha <= 0:
                rep.error(mod, 59, f"<alpha> (mean Vs) = {fmt_num(inc.alpha)}")
            embedded = any(v.P[n][2] < v.gelev - v.tol for n in inter)
            if embedded and int(inc.ngp) <= 0:
                rep.error(mod, 60, f"<ngp> = {inc.ngp} (D-CHK-04)")
            if 2 <= int(wp.cohf) <= 7 and not r.wpass:
                rep.edu(mod, "EDU-26", f"coherency model {wp.cohf} needs wave passage (HOUSE <wpass> = 1)")
            if int(inc.hseed) == 0 and int(inc.vseed) == 0 and float(inc.randphz) == 0 or int(inc.nmodes) < 0:
                rep.edu(mod, "EDU-13")
            self._incoherency_model(mod, inc, wp)
            self._projections(mod, inter)
            self._time_grid(mod, harmonic_ok=True)
        if r.wpass:
            self._record_problems(mod, r.rec["WPASS"])
        if (r.wpass or r.me) and (symm or int(h.dim) == 1):        # full 3D models only (manual 2.7, 6.5.4)
            what = " and ".join(t for t, on in (("wave passage", r.wpass), ("multiple excitation", r.me)) if on)
            rep.edu(mod, "EDU-26", f"{what} not allowed " + ("with SYMM (G-17)" if symm else "in a 2D model (G-17)"))
        if r.me:
            if not r.wpass:
                rep.edu(mod, "EDU-26", "multiple excitation needs wave passage (HOUSE <wpass> = 1)")
            zones = get_entries(m, "ME")
            if not zones:
                rep.error(mod, 115)
            iset = set(inter)
            for no, z in zones:
                n1, n2 = int(z.nfirst), int(z.nlast)
                ok1 = n1 > 0 and n1 in m.nodes and n1 in iset
                if not ok1:
                    rep.error(mod, 116, f"node {n1}", i=no)
                if n2 <= 0 or n2 not in m.nodes or n2 not in iset or (ok1 and n2 < n1):
                    rep.error(mod, 117, f"node {n2}", i=no)
            nf = len(r.fnums or [])
            cplx = int(h.cmplxspec) == 1
            for no, vals in sorted(m.amp.items()):
                if cplx:
                    # complex SAR: [0, 10] is checked on the modulus (D-INC-10)
                    pairs = [complex(vals[k], vals[k + 1] if k + 1 < len(vals) else 0.0)
                             for k in range(0, len(vals), 2)]
                    bad = [abs(c) for c in pairs if not 0.0 <= abs(c) <= 10.0]
                    count = len(vals) / 2.0
                else:
                    bad = [a for a in vals if not 0.0 <= a <= 10.0]
                    count = len(vals)
                if bad:
                    rep.error(mod, 118, f"{bad[0]:.4g}", i=no)
                if r.fnums is not None and count != nf:
                    rep.error(mod, 119, f"{fmt_num(count)} ratios, {nf} frequencies", i=no)
            for no, z in zones:                      # a zone without AMP ratios has no excitation (D-INC-10)
                if int(no) not in m.amp and r.fnums is not None:
                    rep.error(mod, 119, f"0 ratios, {nf} frequencies", i=no)
        if r.coh == 1 or r.me or r.wpass:
            self._freq_set(mod, full=False)          # the frequencies of incoherency / WPASS / ME (HOUSE deck)
        self._house_site_layers(mod)
        if int(r.rec["HOUSEX"].nlssi) == 1:
            self._nonlinear_soil(mod)
        # memory (G-22, D-ANL-10)
        if inter:
            need = analys_memory_bytes(len(inter))
            ram = physical_ram()
            frac = eduopt_float(m, "MEMLIMIT", 0.8)
            if ram and need > frac * ram:
                rep.edu(mod, "EDU-20", f"{len(inter)} interaction nodes: about {need / 1e9:.3g} GB needed, "
                                       f"{ram / 1e9:.3g} GB physical memory (D-ANL-10)")
            if len(inter) > 20000 and r.limits != "UNLIMITED":
                rep.edu(mod, "EDU-02", f"{len(inter)} interaction nodes > 20,000 (practical limit)")

    def _nonlinear_soil(self, mod: str) -> None:
        """Non-linear soil SSI (HOUSEX <nlssi> = 1, requirements 4.4 item 7): the .pin records PIN / PINGRP /
        PINMAT (EDU-41 for an invalid group or material line: no SOLID/PLANE group of the model dimension,
        excavated elements -- ETYPE 2, or ETYPE 0 resolved as excavated soil (D-HOU-01) --, unknown
        material, unknown or unusable curve), EDU-42 when there is neither a PINGRP record nor a
        ``<model>.pin`` file, EDU-43 when STRESS does not compute the strains (<iter> = 0)."""
        from .commands.nlsoil_cmds import has_pin_data, pin_from_model
        rep, m, r = self.rep, self.m, self.r
        if has_pin_data(m):
            pin, probs, _notes, _mn = pin_from_model(m, view=self.view)
            for g, text in probs:
                rep.edu(mod, "EDU-41", text, g=g)
            dim = int(r.house.dim)
            want = PLANE if dim == 1 else SOLID
            for grp in (pin.groups if pin is not None else []):
                if int(m.groups[grp.igrp].type) != want:
                    rep.edu(mod, "EDU-41", f"group {grp.igrp}: a {'2D' if dim == 1 else '3D'} model needs "
                                           f"{'PLANE' if dim == 1 else 'SOLID'} nonlinear soil elements", g=grp.igrp)
        elif self.r.file(f"{m.name}.pin") is None:
            rep.edu(mod, "EDU-42", "HOUSEX <nlssi> = 1 but no PINGRP is defined and the model directory has no "
                                   f"{m.name}.pin: HOUSE analyses the model with its linear properties (no FILE78; "
                                   "STRESS <iter> = 1 then computes no strains)")
        if int(r.rec["STRESS"].iter) != 1:
            rep.edu(mod, "EDU-43", "HOUSEX <nlssi> = 1 but STRESS <iter> = 0: STRESS does not compute the effective "
                                   "strains (FILE74) of the next iteration")

    def _house_site_layers(self, mod: str) -> None:
        """Errors 19-25 / EDU-04 of the free-field layers that AFWRITE copies into the HOUSE deck.

        The HOUSE deck carries the TOPL layers with their L properties and the SITE half-space row
        (``sitelayers``: interface depths of the interaction nodes, layer of the excavated elements).
        They are shared data, so their errors are reported under HOUSE as well as under SITE and
        block both decks (D-CHK-01, D-AFW-01).  A model without top layers writes no site layers
        (Error 46 belongs to SITE).  The half-space row is needed unless the base is rigid without
        a half-space layer (``<hs>`` = 0 and ``<nl>`` = 0).  The Poisson-ratio warning (EDU-11) is
        reported under SITE only.
        """
        m, s = self.m, self.r.site
        if not m.topl:
            return
        for l in dict.fromkeys(m.topl):
            self._layer(mod, l, nu_check=False)
        hs = int(s.hs)
        if hs or int(s.nl) > 0:
            self._layer(mod, hs, halfspace=True, nu_check=False)

    def _incoherency_model(self, mod: str, inc, wp) -> None:
        """Incoherent analysis (HOUSE <coh> = 1, requirements 4.4 item 6): the coherency model must be
        available (D-INC-04: models 2 and 4 refuse to run, their coefficients could not be confirmed),
        the directionality factor of models 2-7 lies in [0, 1] (Error 114), model 7 needs its five
        user files, and the superposition / stochastic options must be compatible (spec 05b section 8)."""
        from ..core import coherency as COH
        from ..core.incoherency import MAX_SIM, stochastic
        rep, r = self.rep, self.r
        cohf = int(wp.cohf)
        if 2 <= cohf <= 6 and not COH.model_available(cohf):
            rep.edu(mod, "EDU-26", COH.unavailable_reason(cohf))
        if 2 <= cohf <= 7 and not 0.0 <= float(inc.alpha) <= 1.0:
            rep.error(mod, 114, f"<alpha> = {fmt_num(inc.alpha)}: directionality factor of coherency model {cohf} "
                                "must lie in [0, 1]")
        if cohf == 7:
            miss = [n for n in (COH.USER_FREQ, COH.USER_DIST) + COH.USER_TABLES if r.file(n) is None]
            if miss:
                rep.edu(mod, "EDU-26", "user coherency model 7 needs the files " + ", ".join(miss)
                        + " in the model directory (spec 05b 3.16)")
        hx = r.rec["HOUSEX"]
        nsim, supmode, nmodes = int(hx.nsim), int(hx.supmode), int(inc.nmodes)
        if stochastic(int(inc.hseed), int(inc.vseed), float(inc.randphz)):
            if not 1 <= nsim <= MAX_SIM:
                rep.edu(mod, "EDU-26", f"HOUSEX <nsim> = {nsim} stochastic simulations (1..{MAX_SIM})")
            if nmodes < 0:
                rep.edu(mod, "EDU-26", f"<nmodes> = {nmodes} (one mode) with stochastic simulation")
            if supmode == 1:
                rep.edu(mod, "EDU-26", "Quadratic (SRSS) superposition with stochastic simulation (HSeed/VSeed, "
                                       "RandPhz): SRSS is deterministic")
        elif supmode == 1 and nmodes >= 0:
            rep.edu(mod, "EDU-26", f"Quadratic (SRSS TF) superposition needs one run per mode: <nmodes> = -k, "
                                   f"not {nmodes}")

    def _projections(self, mod: str, inter: List[int]) -> None:
        """EDU-17 (G-18, D-INC-09): horizontal projections of interaction nodes of different levels
        closer than 0.1 x the mesh size."""
        if len(inter) < 2:
            return
        from scipy.spatial import cKDTree
        v = self.view
        X = np.array([v.P[n] for n in inter])
        t3 = cKDTree(X)
        d, _ = t3.query(X, k=2)
        h = float(np.median(d[:, 1])) if len(X) > 1 else 0.0
        if h <= 0:
            return
        t2 = cKDTree(X[:, :2])
        pairs = [(i, j) for i, j in t2.query_pairs(0.1 * h) if abs(X[i, 2] - X[j, 2]) > v.tol]
        if pairs:
            i, j = pairs[0]
            self.rep.edu(mod, "EDU-17", f"{len(pairs)} pairs, e.g. nodes {inter[i]} and {inter[j]}; "
                                        f"consider per-level incoherency (D-INC-09)")

    # ================================================================== FORCE
    def check_force(self, mod: str) -> None:
        rep, m, r = self.rep, self.m, self.r
        self._gravity(mod)
        self._time_grid(mod, harmonic_ok=True)
        self._freq_set(mod)
        if not m.forces and not m.moments:
            rep.error(mod, 61)
        for table, wno in ((m.forces, 10), (m.moments, 11)):
            for n, ld in sorted(table.items()):
                nd = m.nodes.get(n)
                if nd is None:
                    rep.error(mod, 62, i=n)
                    continue
                off = 0 if wno == 10 else 3
                if any(ld.factor[k] != 0 and nd.fix[off + k] for k in range(3)):
                    rep.warning(mod, wno, n=n)

    # ================================================================== ANALYS
    def _symmetry_geometry(self, mod: str, symm) -> None:
        """SYMM geometry rules of manual 9.2.39 (EDU-26; Errors 2-5 are reported by check_house): 3D planes
        parallel to XZ or YZ and orthogonal to each other, a 2D line parallel to Z (one line at most), and
        the model on one side of every plane (a half / quarter model; D-ANL-12, sassi.core.symmetry)."""
        if not symm:
            return
        from ..core import symmetry as SYM
        v = self.view
        rows = [dict(no=int(no), type=int(rec.get("type")), n1=int(rec.get("node1")), n2=int(rec.get("node2")),
                     n3=int(rec.get("node3"))) for no, rec in symm]
        planes, errs, _ = SYM.parse_planes(rows, v.P, int(self.r.house.dim), v.tol)
        for e in errs:                       # Errors 2-5 and the plane count are reported by check_house
            if not e.startswith("Error ") and "SYMM planes (at most" not in e:
                self.rep.edu(mod, "EDU-26", e.replace("EDU-26: ", ""))
        if errs or not planes:
            return
        nodes = sorted(set(v.node_elems) | set(v.interaction_nodes()))
        P = np.array([v.P[n] for n in nodes if n in v.P]).reshape(-1, 3)
        for no, neg, pos in SYM.side_violations(planes, P, v.tol):
            self.rep.edu(mod, "EDU-26", f"SYMM plane {no}: nodes on both sides of the plane ({neg} / {pos}); a half / "
                                        "quarter model lies on one side of every symmetry plane")

    def check_analys(self, mod: str) -> None:
        rep, m, r = self.rep, self.m, self.r
        a = r.analys
        self._gravity(mod)
        self._record_problems(mod, a)
        self._time_grid(mod, harmonic_ok=True)
        if int(a.fopt) == 0:
            self._freq_set(mod)
        if int(a.type) not in (0, 1):
            rep.edu(mod, "EDU-26", f"ANALYS <type> = {a.type}")
        if int(a.mode) not in (0, 1, 2, 3):
            rep.edu(mod, "EDU-26", f"ANALYS <mode> = {a.mode} (0..3)")
        if int(a.impe) not in (0, 1, 2):
            rep.edu(mod, "EDU-26", f"ANALYS <impe> = {a.impe} (0..2)")
        elif int(a.impe) > 0 and get_entries(m, "SYMM"):
            rep.edu(mod, "EDU-26", "global impedance is not allowed with SYMM (G-17)")
        elif int(a.impe) > 0 and int(r.house.dim) == 1:
            rep.edu(mod, "EDU-26", "global impedance is available for 3D models only (requirements 4.6 item 9)")
        if int(r.house.dim) == 1 and int(a.type) == 0:              # 2D: in-plane input (x' or z')
            ang = float(a.ang) % 360.0
            if min(abs(ang), abs(ang - 180.0)) > 1e-9:
                rep.edu(mod, "EDU-26", f"2D model: the coordinate transformation angle must be 0 or 180 deg "
                                       f"(<ang> = {fmt_num(a.ang)})")
            if int(r.site.cm) == 1 and int(a.simul) == 0:
                rep.edu(mod, "EDU-26", "2D model: SITE control direction y' (SH / Love) is anti-plane; the in-plane "
                                       "2D model needs x' or z' input (D-W2-07)")
        if int(a.simul) >= 1 and int(a.type) == 0 and r.coh == 0 and float(a.ang) != 0.0:
            rep.edu(mod, "EDU-26", "simultaneous X/Y/Z seismic cases need <ang> = 0 (spec 05a section 7.2)")
        if r.coh == 1 and int(a.simul) > 50:
            rep.edu(mod, "EDU-02", f"{a.simul} incoherent simulations > 50")
        self._analys_incoherency(mod)

    def _analys_incoherency(self, mod: str) -> None:
        """ANALYS with incoherent / wave-passage / ME input (requirements 4.6 item 5, D-ANL-06, D-INC-12):
        seismic <simul> > 1 needs incoherent stochastic simulation with the HOUSE number of samples; the
        stochastic X/Y/Z cases need the angle 0; FFM is recommended for surface foundations only."""
        from ..core.incoherency import stochastic
        rep, r = self.rep, self.r
        a = r.analys
        if int(a.type) != 0:
            return
        simul = int(a.simul)
        inc = r.rec["INCOH"]
        sto = r.coh == 1 and stochastic(int(inc.hseed), int(inc.vseed), float(inc.randphz))
        if simul > 1 and r.coh == 0:
            rep.edu(mod, "EDU-26", f"<simul> = {simul}: seismic simultaneous cases above 1 are incoherent stochastic "
                                   "simulations (HOUSE <coh> = 1)")
        if sto:
            nsim = int(r.rec["HOUSEX"].nsim)
            if simul != nsim:
                rep.edu(mod, "EDU-26", f"stochastic simulation: ANALYS <simul> = {simul} must equal the HOUSE number of "
                                       f"simulations HOUSEX <nsim> = {nsim} (spec 05c A.5.4)")
            if float(a.ang) != 0.0:
                rep.edu(mod, "EDU-26", "stochastic X/Y/Z cases (FILE1X/Y/Z) need <ang> = 0")
        elif r.coh == 1 and simul > 1:
            rep.edu(mod, "EDU-26", f"<simul> = {simul} with deterministic incoherency (HSeed = VSeed = 0 or RandPhz = 0): "
                                   "use 0 (FILE1) or 1 (FILE1X/Y/Z)")
        ffm = int(r.rec["ANALYSX"].ffm)
        if ffm and not (r.coh or r.wpass or r.me):
            rep.edu(mod, "EDU-28", "ANALYSX <ffm> = 1 (free-field motion) without incoherency, wave passage or "
                                   "multiple excitation: ignored")
        elif ffm:
            v = self.view
            if any(v.P[n][2] < v.gelev - v.tol for n in v.interaction_nodes()):
                rep.edu(mod, "EDU-28", "free-field motion (FFM) with embedded interaction nodes: FFM is recommended for "
                                       "surface foundations only, use FFL (D-INC-12, spec 05c A.5.9)")

    # ================================================================== MOTION
    def check_motion(self, mod: str) -> None:
        rep, m, r = self.rep, self.m, self.r
        mo, mx = r.motion, r.rec["MOTIONX"]
        self._gravity(mod)
        for kind, num, det in time_grid_problems(0.0, r.delt, r.nft_raw, harmonic_ok=False):
            rep.add(mod, kind, num, det if num != 9 else f"NFFT = {r.nft_raw}; {r.nft} is written")
        if r.nft > 65536 and r.limits != "UNLIMITED":
            rep.edu(mod, "EDU-02", f"NFFT {r.nft} > 65,536 (D-CNV-11)")
        self._record_problems(mod, mo, shared=False)
        allpts = any(int(mx.get(k)) for k in ("savetf", "saveacc", "savers", "saverot", "rsttf", "rstacc", "rstrs"))
        if not m.nout and not allpts and not int(mo.cnvrt):
            rep.error(mod, 64)
        seen: Dict[Tuple[int, int], int] = {}
        rs_req = bool(int(mx.savers) or int(mx.rstrs))
        for req in m.nout:
            if len(req.codes) > 4 and int(req.codes[4]):       # flag 5: Save Acceleration and Velocity R.S.
                rs_req = True
            for n in req.nodes:
                if n <= 0 or n not in m.nodes:
                    rep.error(mod, 65, n=n)
                key = (n, req.dir)
                if key in seen:
                    rep.edu(mod, "EDU-16", "Error 66 downgraded: flags OR'ed at AFWRITE (D-MOT-10)", n=n)
                seen[key] = 1
        rs_req = rs_req and int(mo.out) == 0
        for d in m.damp:
            if not 0.0 < d < 1.0:
                rep.error(mod, 72, f"damping {fmt_num(d)}")
        if len(m.damp) > 5 and r.limits != "UNLIMITED":
            rep.edu(mod, "EDU-02", f"{len(m.damp)} damping values > 5 (MOTION)")
        if rs_req:
            # MOTION needs a usable RS grid when spectra are requested (spec 05c B.8)
            if float(mo.freq1) <= 0 and float(mo.freq1) >= 0:
                rep.error(mod, 69, f"<freq1> = {fmt_num(mo.freq1)}: must be > 0 when RS are requested")
            if float(mo.freq2) <= float(mo.freq1) and float(mo.freq2) >= 0:
                rep.error(mod, 70, f"<freq2> = {fmt_num(mo.freq2)} <= <freq1>")
            if 0 <= int(mo.fstep) < 2:
                rep.error(mod, 71, f"<fstep> = {mo.fstep}: at least 2 RS frequencies")
            if not m.damp:
                rep.edu(mod, "EDU-28", "no RS damping ratios (DAMP): no response spectra are computed")
        if rs_req and 0 < int(mo.fstep) < 301:
            rep.edu(mod, "EDU-15", f"<fstep> = {mo.fstep} (SRP 3.7.1)")
        if float(mo.smo) != 0 and (r.coh == 0 or int(mo.interp) == 6):
            why = "interpolation option 6" if int(mo.interp) == 6 else "coherent input"
            rep.edu(mod, "EDU-14", f"<smo> = {fmt_num(mo.smo)} with {why}")
        if int(mo.out) == 0 or int(mo.cnvrt):
            self._history(mod)
            self._quiet_zone(mod)
        else:
            for kind, num, det in history_problems(mo.mult, mo.get("max"), mo.rec1, mo.rec2):
                rep.add(mod, kind, num, det)

    def _quiet_zone(self, mod: str) -> None:
        """EDU-19 (spec 01 section 5): the trailing zeros must let the free vibration decay to 1 %:
        ``NFFT dt - record length >= ln(100)/(2 pi f_min beta_min)``.

        ``f_min`` is "the lowest significant system frequency"; CHECK has no eigen-analysis, so it
        uses the fundamental frequency of the soil column ``1/(4 sum h/Vs)`` of the TOPL layers (not
        below the first SSI frequency).  ``beta_min`` is the smallest positive damping ratio of the
        layers, materials and springs.
        """
        r, m = self.r, self.m
        used = r.records_used()
        if used is None or not r.fnums or r.df is None or r.delt <= 0:
            return
        lay = [m.layers[l] for l in m.topl if l in m.layers]
        if not lay or any(L.vs <= 0 or L.thick <= 0 for L in lay):
            return
        f_soil = 1.0 / (4.0 * sum(L.thick / L.vs for L in lay))
        betas = [L.sdamp for L in lay if L.sdamp > 0] + [mt.sdamp for mt in m.materials.values() if mt.sdamp > 0] + \
            [s.damp for s in m.springs.values() if s.damp > 0]
        if not betas:
            return
        fmin = max(f_soil, min(r.fnums) * r.df)
        need = math.log(100.0) / (2 * math.pi * fmin * min(betas))
        have = (r.nft - used) * r.delt
        if have < need:
            self.rep.edu(mod, "EDU-19", f"{have:.4g} s < ln(100)/(2 pi f_min beta_min) = {need:.4g} s "
                                        f"(f_min = {fmin:.4g} Hz, beta_min = {min(betas):.4g})")

    # ================================================================== STRESS
    def check_stress(self, mod: str) -> None:
        rep, m, r = self.rep, self.m, self.r
        self._gravity(mod)
        for kind, num, det in time_grid_problems(0.0, r.delt, r.nft_raw, harmonic_ok=False):
            rep.add(mod, kind, num, det if num != 9 else f"NFFT = {r.nft_raw}; {r.nft} is written")
        if int(r.motion.step) < 0:
            rep.error(mod, 67, f"MOTION <step> = {r.motion.step}")
        nl_strains = int(r.rec["STRESS"].iter) == 1      # non-linear soil strains are an output (4.10 item 5) ...
        if not m.eout and not (nl_strains and int(r.rec["HOUSEX"].nlssi) == 1):   # ... of a non-linear HOUSE run only
            rep.error(mod, 79)
        if nl_strains and int(r.rec["HOUSEX"].nlssi) != 1:
            rep.edu(mod, "EDU-43", "STRESS <iter> = 1 needs FILE78 of a non-linear SSI HOUSE run (HOUSEX <nlssi> = 1)")
        seen: Set[Tuple[int, int]] = set()
        for req in m.eout:
            g = m.groups.get(req.group)
            if g is None or req.group <= 0:
                rep.error(mod, 80, g=req.group)
                continue
            for e in req.elements:
                if e <= 0 or e not in g.elements:
                    rep.error(mod, 81, e=e, g=req.group)
                if (req.group, e) in seen:
                    rep.error(mod, 82, e=e, g=req.group)
                seen.add((req.group, e))
        self._history(mod)
        self._thshlstr_file(mod)

    def _thshlstr_file(self, mod: str) -> None:
        """TSHELL face stresses (THSHLSTR, D-TSH-01): STRESS reads the flag from ``THSHLSTR.opt`` in the model
        directory (the STRESS deck has no ``thshlstr`` parameter); a missing or stale file is rewritten from
        the model record (:func:`sassi.prep.commands.thickshell.sync_thshlstr_file`), so AFWRITE always
        leaves the flag of the record for STRESS.  Without a model directory a THSHLSTR,1 cannot reach
        STRESS (EDU-28)."""
        from .commands.thickshell import THSHLSTR_FILE, model_thshlstr, sync_thshlstr_file
        mdir = Path(self.m.path) if getattr(self.m, "path", None) else None
        if mdir is None or not mdir.is_dir():
            if model_thshlstr(self.m) == 1:
                self.rep.edu(mod, "EDU-28", f"THSHLSTR,1 but the model directory is not defined (MDL): STRESS cannot "
                                            f"read {THSHLSTR_FILE} and computes no TSHELL face stresses")
            return
        try:
            sync_thshlstr_file(self.m, mdir)
        except OSError as exc:
            self.rep.edu(mod, "EDU-28", f"{THSHLSTR_FILE} could not be written in {mdir} ({exc}): STRESS uses the "
                                        "flag of the existing file or 0")

    # ================================================================== RELDISP
    def check_reldisp(self, mod: str) -> None:
        rep, m, r = self.rep, self.m, self.r
        self._gravity(mod)
        for kind, num, det in time_grid_problems(0.0, r.delt, r.nft_raw, harmonic_ok=False):
            rep.add(mod, kind, num, det if num != 9 else f"NFFT = {r.nft_raw}; {r.nft} is written")
        self._freq_set(mod, full=False)
        rd = r.rec["RELD"]
        if not m.rdnd and not int(rd.reldispsall):
            rep.error(mod, 64)
        seen: Set[int] = set()
        for req in m.rdnd:
            if req.node <= 0 or req.node not in m.nodes:
                rep.error(mod, 65, n=req.node)
            if req.node in seen:
                rep.error(mod, 66, n=req.node)
            seen.add(req.node)
        self._history(mod)

    # ================================================================== NONLINEAR (P2)
    def check_nonlinear(self, mod: str) -> None:
        """Option NON checks.  The manual's numbered errors (121-128) are applied here as catalogued
        (Error 121 when no panel and no spring is given; Errors 126/128 are lifted by EDUOPT,NONEXT,1 for
        force 3 / disp 2, D-NON-05); every further finding of
        :func:`sassi.prep.commands.nonlinear_cmds.build_eql` (the rules AFWRITE applies: D-NON-12,
        EDU-09, groups, materials, backbones) is reported as EDU-44 (error) or EDU-45 (warning)."""
        from .commands.nonlinear_cmds import build_eql
        from .options import eduopt
        rep, m = self.rep, self.m
        nonext = eduopt(m, "NONEXT") == "1"
        panels = m.options.entries("P")
        springs = m.options.entries("S")
        if not panels and not springs:
            rep.error(mod, 121)
        for key, rec in panels:
            k = rec.integer(1, 0)
            gid = rec.integer(2, 0)
            g = m.groups.get(gid)
            if g is None:
                rep.error(mod, 123, i=gid, k=k)
            else:
                for mid in sorted({e.mat for e in g.elements.values()}):
                    if mid not in m.materials:
                        rep.error(mod, 122, i=mid, k=k)
            disp, force = rec.integer(4, 1), rec.integer(5, 1)
            if force != 1 and not (nonext and force == 3):
                rep.error(mod, 126, i=force, k=k)
            if disp != 1 and not (nonext and disp == 2):
                rep.error(mod, 128, i=disp, k=k)
        for key, rec in springs:
            force = rec.integer(6, 4)
            if force != 4:
                rep.error(mod, 127, i=force, k=rec.integer(1, 0))
        numbered = rep.numbers("Error", mod)
        _, errs, warns = build_eql(m)
        for e in errs:
            mm = re.match(r"\s*Error\s+(\d+)", e)
            if mm and int(mm.group(1)) in numbered:
                continue                                   # already reported with its number
            if e.startswith("no nonlinear element type") and 121 in numbered:
                continue
            if "Force Option" in e and 126 in numbered:
                continue
            rep.edu(mod, "EDU-44", e)
        for w in warns:
            rep.edu(mod, "EDU-45", w)


# ======================================================================================
# Model-check algorithms shared with the commands (spec 09 section 3)
# ======================================================================================
def excstrchk(v: ModelView) -> Dict[int, List[str]]:
    """EXCSTRCHK: excavation interior nodes that also belong to a structure element (resolved
    ETYPE != 2) or to a BEAMS, SPRING or GENERAL element; node -> descriptions of those elements."""
    exc = v.excavation()
    out: Dict[int, List[str]] = {}
    for n in sorted(exc["interior"]):
        owners = []
        for r in v.node_elems.get(n, []):
            if r.type in (BEAMS, SPRING, GENERAL) or (r.etype != 2 and r.type in (SOLID, PLANE, SHELL, TSHELL)):
                owners.append(f"{ELEMENT_TYPE_NAMES.get(r.type)} element {r.id} group {r.group}")
        if n in v.k_nodes:
            owners.append("beam K node")
        if owners:
            out[n] = owners
    return out


def excstrchk_severity(v: ModelView, n: int, interaction: Optional[Set[int]] = None,
                       excavation_nodes: Optional[Set[int]] = None) -> str:
    """'Error' or 'Warning' for an EXCSTRCHK node ``n`` (the rule of HOUSE,
    :func:`sassi.core.house_lib.excstrchk_kind`): an error unless ``n`` is an interaction node and every
    structural element at it is a SOLID/PLANE element built on the excavation mesh.  Pass the interaction
    and excavation node sets when classifying many nodes."""
    iset = set(v.interaction_nodes()) if interaction is None else interaction
    exc = set(v.excavation()["nodes"]) if excavation_nodes is None else excavation_nodes
    owners = [(r.type, r.dof_nodes) for r in v.node_elems.get(n, [])
              if r.type in (BEAMS, SPRING, GENERAL) or (r.etype != 2 and r.type in (SOLID, PLANE, SHELL, TSHELL))]
    if n in v.k_nodes:
        owners.append((0, [n]))
    return excstrchk_kind(owners, exc, n in iset)


def fixed_interaction(v: ModelView) -> List[Tuple[int, List[str]]]:
    """FIXEDINT: interaction nodes with any fixed translation -> [(node, fixed labels)]."""
    out = []
    for n in v.interaction_nodes():
        fx = v.m.nodes[n].fix
        labs = [lab for lab, k in (("UX", 0), ("UY", 1), ("UZ", 2)) if fx[k]]
        if labs:
            out.append((n, labs))
    return out


def free_springs(v: ModelView) -> List[Tuple[int, bool]]:
    """FREESPRING: nodes without any fixed DOF connected only to SPRING elements -> [(node, has_mass)]."""
    out = []
    m = v.m
    for n, lst in sorted(v.node_elems.items()):
        if not lst or any(r.type != SPRING for r in lst):
            continue
        nd = m.nodes.get(n)
        if nd is None or any(nd.fix):
            continue
        mass = any(x != 0 for x in m.tmass.get(n, [])) or any(x != 0 for x in m.rmass.get(n, []))
        out.append((n, mass))
    return out


def hinges(v: ModelView, eps: float = 1e-3) -> List[Tuple[int, str]]:
    """HINGED: single-point connections of 6-DOF elements (BEAMS, SHELL, TSHELL) to SOLID/PLANE
    elements, and beam-shell drilling joints (one out-of-plane beam at coplanar shells)."""
    out: List[Tuple[int, str]] = []
    solid_nodes = {n for n, lst in v.node_elems.items() if any(r.type in (SOLID, PLANE) for r in lst)}

    def penetrates(r: ElemRef) -> bool:
        """A 6-DOF element with two or more nodes on the solids transmits moments as force couples."""
        return sum(1 for n in dict.fromkeys(r.dof_nodes) if n in solid_nodes) >= 2

    for r in v.elems:
        if r.type not in (BEAMS, SHELL, TSHELL):
            continue
        touch = [n for n in dict.fromkeys(r.dof_nodes) if n in solid_nodes]
        if len(touch) == 1:
            n = touch[0]
            if any(o.type in (BEAMS, SHELL, TSHELL) and penetrates(o) for o in v.node_elems.get(n, [])):
                continue          # an embedded (penetrating) beam or shell carries the moment
            out.append((n, f"{ELEMENT_TYPE_NAMES[r.type]} element {r.id} group {r.group} meets "
                           f"SOLID/PLANE elements at one node only: no moment is transmitted"))
    for n, lst in sorted(v.node_elems.items()):
        if not any(r.type == SHELL for r in lst):
            continue
        beams = [r for r in lst if r.type == BEAMS]
        if not beams:
            continue
        nrm = v.coplanar_normal(n)
        if nrm is None:
            continue
        outp = []
        for b in beams:
            I, J = b.dof_nodes[:2] if len(b.dof_nodes) >= 2 else (0, 0)
            if not (I in v.P and J in v.P):
                continue
            ax = v.P[J] - v.P[I]
            la = _norm(ax)
            if la > 0 and abs(float(np.dot(ax / la, nrm))) > eps:
                outp.append(b)
        if len(outp) == 1:
            b = outp[0]
            out.append((n, f"beam {b.id} group {b.group} meets the coplanar shells out of their plane: the "
                           f"rotation about the shell normal (drilling) is not transmitted"))
    return out


def kint(v: ModelView) -> List[Tuple[int, bool]]:
    """KINT: beam K nodes that are interaction nodes -> [(node, K-only)]."""
    inter = set(v.interaction_nodes())
    return [(n, n not in v.node_elems) for n in sorted(v.k_nodes & inter)]


def unused_nodes(v: ModelView) -> List[int]:
    """Nodes not referenced by any element (beam K nodes count as used, D-CHK-08)."""
    used = v.used_dof_nodes() | v.k_nodes
    return sorted(n for n in v.m.nodes if n not in used)


def drilling_unrestrained(v: ModelView, eps: float = 1e-4) -> Dict[int, str]:
    """EDU-06: nodes of coplanar Kirchhoff SHELL elements (shells with springs, SOLID or PLANE elements,
    which carry no rotations) whose drilling rotation (about the shell normal) is restrained neither by
    D nor by a rotational spring.  Nodes shared with SOLID/PLANE elements -- a basemat or wall on the
    excavated soil -- are included (final audit: HOUSE reported them as EDU-06 and ANALYS stopped on the
    singular system while CHECK said nothing); FIXROT/FIXSHLROT treat shell-only nodes (manual 9.7), so
    the message asks for D there."""
    out: Dict[int, str] = {}
    m = v.m
    for n, lst in v.node_elems.items():
        types = {r.type for r in lst}
        if SHELL not in types or not types <= {SHELL, SPRING, SOLID, PLANE}:
            continue
        with_solids = bool(types & {SOLID, PLANE})
        nrm = v.coplanar_normal(n, eps)
        if nrm is None:
            continue
        fix = m.nodes[n].fix
        if all(fix[3:6]):
            continue
        axis = [k for k in range(3) if abs(nrm[k]) >= 1.0 - eps]
        if axis and fix[3 + axis[0]]:
            continue
        kd = 0.0
        for r in lst:
            if r.type == SPRING:
                sc = m.springs.get(r.elem.prop)
                if sc is not None:
                    kd += sc.scxx * nrm[0] ** 2 + sc.scyy * nrm[1] ** 2 + sc.sczz * nrm[2] ** 2
        if kd > 0:
            continue
        if axis and with_solids:
            out[n] = (f"normal along global {'XYZ'[axis[0]]}, node shared with SOLID/PLANE elements: fix ROT"
                      f"{'XYZ'[axis[0]]} with D (FIXROT treats shell-only nodes)")
        elif axis:
            out[n] = f"normal along global {'XYZ'[axis[0]]}: fix ROT{'XYZ'[axis[0]]} (FIXROT)"
        elif with_solids:
            out[n] = (f"oblique shell normal ({nrm[0]:.3f}, {nrm[1]:.3f}, {nrm[2]:.3f}), node shared with "
                      f"SOLID/PLANE elements: add a rotational spring about the normal (FIXROT/FIXSHLROT treat "
                      f"shell-only nodes)")
        else:
            out[n] = (f"oblique shell normal ({nrm[0]:.3f}, {nrm[1]:.3f}, {nrm[2]:.3f}): D cannot restrain it; "
                      f"use FIXROT or FIXSHLROT")
    return out


# ======================================================================================
# helpers
# ======================================================================================
def _file_of(entries: Dict[Any, Any], i: int) -> str:
    rec = entries.get(i)
    return rec.file if rec is not None else ""


class _Placeholders(dict):
    """Leaves a catalogue placeholder whose value is not known as ``<i>``, ``<w>`` ... (the
    manual's notation); the message detail then gives the values."""

    def __missing__(self, key):
        return f"<{key}>"


def format_title(title: str, fmt: Dict[str, Any]) -> str:
    """Catalogue title with its placeholders filled from ``fmt`` (unknown ones kept as ``<key>``)."""
    try:
        return title.format_map(_Placeholders(fmt))
    except (IndexError, ValueError):
        return title


def _file88_layers(path: Path) -> Optional[int]:
    """Number of layer rows of a FILE88 text file (lines starting with an integer layer number)."""
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return None
    n = 0
    for ln in lines:
        t = ln.strip().split()
        if t and not t[0].startswith(("#", "*")) and t[0].isdigit() and len(t) >= 3:
            n += 1
    return n or None


def run_check(model, modules: Optional[Sequence[str]] = None, dirs: Sequence[Union[str, Path]] = ()) -> CheckReport:
    """Run CHECK for the AOPT-enabled modules (or ``modules``) of ``model``.

    ``dirs`` are the directories in which relative file names (THFILE, RSIN ...) are resolved:
    the model directory first, then the working directory.
    """
    return Checker(model, modules, dirs).run()


__all__ = ["ERRORS", "WARNINGS", "EDU_TEXT", "RUN_ORDER", "MODEL", "MODEL_GATED", "CheckOptions", "CheckMessage",
           "CheckReport", "ModelView", "ElemRef", "Resolved", "Checker", "run_check", "write_err", "excstrchk",
           "fixed_interaction", "free_springs", "hinges", "kint", "unused_nodes", "drilling_unrestrained",
           "dense_matrix_bytes", "analys_memory_bytes", "physical_ram", "resolve_file", "format_title"]
