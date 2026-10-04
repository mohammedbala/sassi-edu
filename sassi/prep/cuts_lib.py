"""Cuts, submodels and section-cut calculations (requirements sections 3.4.J and 4.14; spec 09
section 5; spec 10 sections 1.4 and 3; spec 04 section 8; decisions D-MDL-13, D-SEC-01 ... D-SEC-06,
D-STR-12).

Pure functions on :class:`sassi.model.SSIModel`.  The command handlers of
:mod:`sassi.prep.commands.cuts` parse the arguments, call them and print the collected notes; scripts
and the GUI can call them directly.

What a section cut computes (for the structural engineer)
---------------------------------------------------------
STRESS gives element-centre stresses: SOLID ``Sxx Syy Szz Sxy Sxz Syz`` in global axes, SHELL
membrane stresses ``Sx'x' Sy'y' Sx'y'`` and moments per unit length ``Mx'x' My'y' Mx'y'`` in the
element axes x' y' (D-STR-10).  To design a wall or to check the base shear of a building one needs
the *stress resultants* over a section: the forces and moments that the part of the structure on one
side of a plane exerts on the other part.  A cut is a set of elements; the plane through ``P`` with
normal ``n`` intersects them in *section pieces*:

* a SOLID gives a polygon of area ``A_i`` (the plane-hexahedron intersection);
* a SHELL gives a strip: the segment where the mid-surface meets the plane (length ``L``) times the
  width ``w = t / |m x ez|`` in the plane (``m`` the shell normal: an oblique wall is cut over a width
  larger than its thickness);
* a PLANE element (plane strain, unit thickness out of the X-Z plane) gives a strip of width 1
  (SASSI-EDU extension; the plane normal must lie in the X-Z plane).

With the element-centre stress tensor ``S_i`` taken as uniform over the piece, the traction on the
section is ``S_i . ez`` and (requirements 4.14, spec 10 section 3.2, normative)::

    f_i = A_i (S_i . ez)                                  force of piece i
    mb_i = L m x (Mb . nu_s)                              SHELL plate-bending couple (D-SEC-03)
    F = sum f_i ;  M = sum (c_i - C) x f_i + sum mb_i     about the area centroid C

``F`` and ``M`` are the actions of the **+ez side on the -ez side** (traction with outward normal +ez
of the -ez part), reported in the local cut axes of D-SEC-01::

    ez = n/|n| ;  ex = (r - (r.ez) ez)/|.| ;  ey = ez x ex ;  origin C

as ``Fx Fy`` (in-plane shears), ``Fz`` (normal force, + tension), ``Mx My`` (bending) and ``Mz``
(torsion).  Because every term is linear in the stresses, the resultants are ``R = sum_i T_i s_i``
with one 6x6 operator ``T_i`` per piece (:func:`resultant_operators`), built once and applied to every
stress frame of a history (CALCSECTHIST).

Section properties in the local axes about ``C``: ``A``, ``Ixx = int y'^2 dA``, ``Iyy = int x'^2 dA``,
``Ixy = int x'y' dA``, ``Izz = Ixx + Iyy`` (exact polygon formulas; strips with the rectangle tensor
``(A/12)(L^2 d d^T + w^2 nu nu^T)`` plus the parallel-axis terms).

Decisions applied
-----------------
* D-MDL-13: cuts are session-global sets of (group, element) kept in ``interp.session['cuts']`` as
  ``{cut: set of (group, element)}`` (the layout GCOM / MERGEGROUP and CUTPLOT read); CUTVOL / TRANVOL
  select elements whose DOF nodes all lie in the closed box (tolerance); CUT2SUB copies the elements,
  their nodes (fixities, interaction flags), the referenced M/L/R/SC/MX entries and the masses.
* D-SEC-01 axes, D-SEC-02 output order, D-SEC-03 shell bending on (``EDUOPT,SECTBEND,0`` excludes it),
  D-SEC-04 a face (SOLID) or edge (SHELL/PLANE) lying on the plane is counted once, the element on the
  -n side preferred, D-SEC-05 time ``(k-1) ts`` (step k when ts <= 0), D-SEC-06 warning for the
  ``ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`` maxima.
* D-STR-12: ``.ess`` frames use the ELEMENT_CENTER layout (spec 10 section 3.1) with signed values::

      <number of groups>
      <etype> <group#> <ordered group#> <#elements>       (one line per group)
      <etype> <group#> <ordered group#>                   (block header, per group)
      <element#> <c1> ... <c6>                            (one line per element)

  :func:`read_ess` is tolerant (comment lines ``#``/``*``/``!``, Fortran ``D`` exponents, a block
  whose element count differs from the group table, fewer than six columns).
* Geometric tolerance: D-GEN-06 (``tol = max(1e-6 L_ref, 1e-9)``, ``EDUOPT,GEOMTOL`` overrides), the
  ``1e-6 x model size`` of spec 10 section 3.2.

Where loaded results live (SASSI-EDU design, reported to the lead)
------------------------------------------------------------------
READSTR attaches element stresses to a model; they are *results*, not model definition, so they are
kept in ``model.ui_state['element_stress']`` (saved by SAVE, not written by WRITE, not part of the
model hash).  CSECT stores the exact section pieces of the cross-section model in
``model.ui_state['csect']`` (the pieces keep the original element geometry, so shell stresses stay in
their own x' y' axes), together with the BEAMS that cross the plane and the cut elements that give no
piece, so that CALCPAR / CALCMOI repeat the spec 10 section 3.2 warning listing them.

Which section CALCPAR / CALCMOI integrate (:func:`active_section`, spec 10 sections 3.7-3.8)
--------------------------------------------------------------------------------------------
1. a CSECT model (``ui_state['csect']``): its stored pieces; the normal must be the CSECT normal;
2. a model without CSECT data whose elements form a **unit-thickness extrusion along n** (every DOF
   node on one of two planes one unit apart, e.g. a CSECT model read back by WRITE -> INP): the
   extruded elements are intersected with their mid-plane (spec 10 section 3.7, "alternatively");
3. any other model (Q-11, "base forces and moments of a whole building"): the **whole structure**
   (excavated-soil elements, ETYPE 2, left out) is cut by the plane through its end on the -n side,
   i.e. its base when n points up.  The faces on that plane belong to the bottom row of elements
   (D-SEC-04); the resultants are the actions of the structure on its support.  When that plane only
   touches the structure (n oblique to the base) the command fails: another plane needs a cut and
   CSECT.  (SASSI-EDU decision reported to the lead.)

BEAMS are drawn in a cross-section model but are never section pieces (they have no ``.ess`` record),
so CALCC on a cross-section model (cases 1 and 2) leaves them out: its volume centroid is then the area
centroid of CALCPAR / CALCMOI (spec 10 section 3.3).
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

import numpy as np

from ..conventions import ELEMENT_TYPE_CODES, ELEMENT_TYPE_NAMES
from ..model import SSIModel
from ..model.entities import (BeamSection, CoordSys, Element, ElementRequest, Group, Material, MatrixProp, Node,
                              SoilLayer, SpringProp)
from ..model.geometry import euler_from_matrix, loc_matrix
from .check import BEAMS, GENERAL, PLANE, SHELL, SOLID, SPRING, TSHELL, ModelView
from .generation_lib import Notes, unique_points
from .options import eduopt

PathLike = Union[str, Path]
Key = Tuple[int, int]

#: names of the six section resultants and of the CALCPAR parameters (D-SEC-02)
RESULTANT_NAMES = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")
CALCPAR_NAMES = ("Area", "Xc", "Yc", "Zc", "Ixx", "Iyy", "Ixy", "Izz") + RESULTANT_NAMES
#: element types that give section pieces (SOLID polygons, SHELL / PLANE strips)
PIECE_TYPES = (SOLID, SHELL, PLANE)
#: the 12 edges of the 8-node SOLID (node slots; 1-4 bottom, 5-8 top, spec 09 section 1.3)
HEX_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7))
#: file names of the absolute-maxima tables (not simultaneous values, D-SEC-06)
MAXIMA_NAME = re.compile(r"(ELEMENT_CENTER_ABS_MAX|_ABS_MAX)", re.IGNORECASE)
#: relative tolerance of "parallel" normals (CALCPAR on a CSECT model)
PARALLEL_TOL = 1e-6


class SectionError(ValueError):
    """Invalid input of a cut / section command; nothing was changed."""


def _notes(notes: Optional[Notes]) -> Notes:
    return notes if notes is not None else Notes()


def _short(items: Sequence[Any], k: int = 10) -> str:
    items = list(items)
    s = ", ".join(str(i) for i in items[:k])
    return s + (f", ... ({len(items)} in total)" if len(items) > k else "")


def _ge(keys: Iterable[Key]) -> List[str]:
    return [f"{g}/{e}" for g, e in keys]


# ======================================================================================
# Session-global cuts (D-MDL-13)
# ======================================================================================
def cut_store(interp) -> Dict[int, Set[Key]]:
    """The cut table ``{cut: set of (group, element)}`` of an interpreter session (created when absent)."""
    store = interp.session.get("cuts")
    if not isinstance(store, dict):
        store = {}
        interp.session["cuts"] = store
    return store


def get_cut(interp, cut: int) -> Set[Key]:
    """Elements of cut ``cut``; :class:`SectionError` when the cut is not defined (spec 09 section 5.1)."""
    s = cut_store(interp).get(int(cut))
    if s is None:
        raise SectionError(f"cut {cut} is not defined (CUTADD / CUTVOL / SLICE)")
    return s


def model_elements(m: SSIModel, keys: Iterable[Key]) -> Tuple[List[Key], List[Key]]:
    """Split ``keys`` into (present in the model, missing), both sorted."""
    have, miss = [], []
    for g, e in sorted(set((int(a), int(b)) for a, b in keys)):
        grp = m.groups.get(g)
        (have if grp is not None and e in grp.elements else miss).append((g, e))
    return have, miss


# ======================================================================================
# Geometry helpers
# ======================================================================================
def model_view(m: SSIModel) -> ModelView:
    """:class:`sassi.prep.check.ModelView` (global node coordinates, resolved ETYPE, tolerance D-GEN-06)."""
    return ModelView(m)


def unit(v: Sequence[float], what: str = "vector") -> np.ndarray:
    a = np.asarray(v, dtype=float).reshape(3)
    n = float(np.linalg.norm(a))
    if not n > 0.0 or not math.isfinite(n):
        raise SectionError(f"the {what} must not be zero")
    return a / n


def section_axes(normal: Sequence[float], right: Sequence[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Local cut axes of D-SEC-01 / spec 10 section 1.4: ``ez = n/|n|``, ``ex`` = r projected on the
    plane (Gram-Schmidt), ``ey = ez x ex`` (right-handed, ``ex x ey = ez``)."""
    ez = unit(normal, "plane normal")
    r = np.asarray(right, dtype=float).reshape(3)
    rn = float(np.linalg.norm(r))
    if not rn > 0.0:
        raise SectionError("the right vector must not be zero")
    rp = r - float(r @ ez) * ez
    if float(np.linalg.norm(rp)) < 1e-8 * rn:
        raise SectionError("right vector is parallel to the plane normal")
    ex = rp / float(np.linalg.norm(rp))
    ey = np.cross(ez, ex)
    return ex, ey, ez


def in_plane_basis(ez: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Any orthonormal (u, v) of the plane with ``u x v = ez`` (used to order polygon vertices)."""
    a = np.array([1.0, 0.0, 0.0]) if abs(ez[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = a - float(a @ ez) * ez
    u /= float(np.linalg.norm(u))
    return u, np.cross(ez, u)


def element_dof_nodes(gtype: int, nodes: Sequence[int]) -> List[int]:
    """Distinct nodes carrying DOFs (BEAMS / 3-node GENERAL: I and J; the K node is orientation)."""
    ns = list(nodes[:2]) if gtype in (BEAMS, GENERAL) else list(nodes)
    return list(dict.fromkeys(n for n in ns if n))


def model_extents(v: ModelView) -> Tuple[np.ndarray, np.ndarray]:
    """(min, max) global coordinates of the model's nodes."""
    if not v.P:
        raise SectionError("the model has no nodes")
    A = np.array(list(v.P.values()))
    A = A[np.all(np.isfinite(A), axis=1)]
    if not len(A):
        raise SectionError("the model has no node with defined coordinates")
    return A.min(axis=0), A.max(axis=0)


def box_bounds(v: ModelView, given: Sequence[Optional[float]]) -> np.ndarray:
    """``[Xmin, Xmax, Ymin, Ymax, Zmin, Zmax]`` with blank (None) bounds replaced by the model extents
    (spec 09 section 5.7)."""
    lo, hi = model_extents(v)
    out = np.empty(6)
    for k in range(6):
        val = given[k] if k < len(given) else None
        out[k] = (lo if k % 2 == 0 else hi)[k // 2] if val is None else float(val)
    for a in range(3):
        if out[2 * a] > out[2 * a + 1]:
            raise SectionError(f"box: {'XYZ'[a]}min {out[2 * a]:g} > {'XYZ'[a]}max {out[2 * a + 1]:g}")
    return out


def elements_in_box(m: SSIModel, bounds: Sequence[float], v: Optional[ModelView] = None) -> List[Key]:
    """Elements whose DOF nodes **all** lie inside the closed box (tolerance), D-MDL-13."""
    v = v or model_view(m)
    b = np.asarray(bounds, float)
    lo = b[0::2] - v.tol
    hi = b[1::2] + v.tol
    out = []
    for g, e in m.iter_elements():
        ns = element_dof_nodes(g.type, e.nodes)
        if not ns or not v.defined(ns):
            continue
        X = np.array([v.P[n] for n in ns])
        if np.all(X >= lo) and np.all(X <= hi):
            out.append((g.id, e.id))
    return out


def elements_on_plane(m: SSIModel, point: Sequence[float], normal: Sequence[float],
                      v: Optional[ModelView] = None) -> List[Key]:
    """SLICE selection (spec 09 section 5.9): the signed distances of the DOF nodes have both signs, or
    one node lies within the tolerance of the plane."""
    v = v or model_view(m)
    ez = unit(normal, "plane normal")
    P = np.asarray(point, float)
    out = []
    for g, e in m.iter_elements():
        ns = element_dof_nodes(g.type, e.nodes)
        if not ns or not v.defined(ns):
            continue
        d = (np.array([v.P[n] for n in ns]) - P) @ ez
        if (d.min() < -v.tol and d.max() > v.tol) or np.any(np.abs(d) <= v.tol):
            out.append((g.id, e.id))
    return out


# ======================================================================================
# Element stresses: .ess frames (D-STR-12) and READSTR storage
# ======================================================================================
@dataclass
class EssData:
    """Content of one ELEMENT_CENTER-layout file: group types and per-element six values."""
    groups: Dict[int, str] = field(default_factory=dict)          # group -> element type name
    values: Dict[Key, np.ndarray] = field(default_factory=dict)   # (group, element) -> (6,)
    source: str = ""

    def __len__(self) -> int:
        return len(self.values)


def _is_word(tok: str) -> bool:
    return bool(tok) and tok[0].isalpha()


def _num(tok: str) -> float:
    return float(tok.replace("D", "E").replace("d", "e"))


def read_ess(path: PathLike) -> EssData:
    """Read an ``.ess`` frame / ``ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`` (layout of spec 10 section 3.1).

    Tolerant reader: blank and comment lines (``#``, ``*``, ``!``) are skipped, a group table that
    does not list every block is accepted, the rows of a block are read up to the next block header,
    Fortran ``D`` exponents are accepted and rows with fewer than six values are padded with zeros.
    A block whose element type differs from the group table is an error (:class:`SectionError`).
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise SectionError(f"{path}: cannot read ({exc.strerror or exc})") from None
    rows = []
    for ln in text.splitlines():
        s = ln.strip()
        if not s or s[0] in "#*!":
            continue
        rows.append(s.replace(",", " ").split())
    if not rows:
        raise SectionError(f"{path.name}: empty file")
    out = EssData(source=str(path))
    try:
        ng = int(_num(rows[0][0]))
    except (ValueError, IndexError):
        raise SectionError(f"{path.name}: the first line must be the number of groups") from None
    i = 1
    table: Dict[int, str] = {}
    while i < len(rows) and len(table) < ng and _is_word(rows[i][0]) and len(rows[i]) >= 4:
        try:
            table[int(_num(rows[i][1]))] = rows[i][0].upper()
        except ValueError:
            raise SectionError(f"{path.name}: malformed group table line {' '.join(rows[i])!r}") from None
        i += 1
    # blocks: a header line (element type first) followed by its element rows
    blocks: List[Tuple[int, List[List[str]]]] = []
    for r in rows[i:]:
        if _is_word(r[0]):
            try:
                g = int(_num(r[1]))
            except (ValueError, IndexError):
                raise SectionError(f"{path.name}: malformed block header {' '.join(r)!r}") from None
            t = r[0].upper()
            if g in table and table[g] != t:
                raise SectionError(f"{path.name}: block of group {g} is {t}, the group table says {table[g]}")
            out.groups[g] = t
            blocks.append((g, []))
            continue
        if not blocks:
            raise SectionError(f"{path.name}: element row {' '.join(r)!r} before the first block header")
        blocks[-1][1].append(r)
    for g, brows in blocks:
        if not brows:
            continue
        try:
            if all(len(r) == 7 for r in brows):          # fast path: one conversion per block
                arr = np.array(" ".join(" ".join(r) for r in brows).replace("D", "E").replace("d", "e").split(),
                               dtype=float).reshape(-1, 7)
            else:
                arr = np.array([[_num(x) for x in r[:7]] + [0.0] * (7 - len(r[:7])) for r in brows])
        except ValueError:
            raise SectionError(f"{path.name}: non-numeric element row in the block of group {g}") from None
        ids = arr[:, 0]
        if not np.all(ids == np.round(ids)):
            raise SectionError(f"{path.name}: non-integer element number in the block of group {g}")
        out.values.update(zip(((g, int(e)) for e in ids), arr[:, 1:7]))
    for g, t in table.items():
        out.groups.setdefault(g, t)
    return out


def write_ess(path: PathLike, groups: Dict[int, str], values: Dict[Key, Sequence[float]]) -> Path:
    """Write the ``.ess`` layout (used by the tests and the verification problems; STRESS writes the
    same layout with :func:`sassi.core.stress_lib.write_element_center`)."""
    path = Path(path)
    gids = sorted(groups)
    order: Dict[str, int] = {}
    ordered = {}
    for g in gids:
        order[groups[g]] = order.get(groups[g], 0) + 1
        ordered[g] = order[groups[g]]
    by_g: Dict[int, List[int]] = {g: [] for g in gids}
    for (g, e) in values:
        by_g.setdefault(g, []).append(e)
    lines = [str(len(gids))]
    for g in gids:
        lines.append(f"{groups[g]:<8s}{g:6d}{ordered[g]:6d}{len(by_g[g]):8d}")
    for g in gids:
        lines.append(f"{groups[g]:<8s}{g:6d}{ordered[g]:6d}")
        for e in sorted(by_g[g]):
            v = list(values[(g, e)]) + [0.0] * 6
            lines.append(f"{e:10d} " + " ".join(f"{float(x):.10e}" for x in v[:6]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def check_ess_against_model(m: SSIModel, ess: EssData, notes: Optional[Notes] = None) -> List[Key]:
    """READSTR rules (spec 10 section 3.9): the groups of the file must exist in the model with the same
    element type (error); records of elements not in the model are ignored (warning); SOLID / SHELL /
    PLANE elements without a record are listed (warning).  Returns the matched (group, element) keys."""
    notes = _notes(notes)
    bad = []
    for g, t in sorted(ess.groups.items()):
        grp = m.groups.get(g)
        if grp is None:
            bad.append(f"group {g} ({t}) is not in the model")
        elif ELEMENT_TYPE_NAMES.get(grp.type, "").upper() != t.upper() and \
                ELEMENT_TYPE_CODES.get(t.upper()) != grp.type:
            bad.append(f"group {g} is {t} in the file but {grp.type_name} in the model")
    if bad:
        raise SectionError("the stress file does not match the active model: " + "; ".join(bad[:10]))
    matched, extra = [], []
    for k in sorted(ess.values):
        grp = m.groups.get(k[0])
        (matched if grp is not None and k[1] in grp.elements else extra).append(k)
    if extra:
        notes.warn(f"{len(extra)} records of elements that are not in the model ignored: {_short(_ge(extra))}")
    have = set(matched)
    missing = [(g.id, e.id) for g, e in m.iter_elements() if g.type in PIECE_TYPES and (g.id, e.id) not in have]
    if missing:
        notes.warn(f"{len(missing)} SOLID/SHELL/PLANE elements have no stress record: {_short(_ge(missing))}")
    return matched


def attach_stresses(m: SSIModel, ess: EssData, keys: Optional[Iterable[Key]] = None, maxima: bool = False) -> int:
    """Store element stresses on a model (``ui_state['element_stress']``, replacing earlier data)."""
    keys = sorted(ess.values) if keys is None else sorted(keys)
    rows = [[int(g), int(e)] + [float(x) for x in ess.values[(g, e)]] for g, e in keys]
    groups = sorted({g for g, _ in keys})
    m.ui_state["element_stress"] = {
        "file": ess.source, "maxima": bool(maxima),
        "groups": [[g, ess.groups.get(g, ELEMENT_TYPE_NAMES.get(m.groups[g].type, "") if g in m.groups else "")]
                   for g in groups],
        "rows": rows}
    return len(rows)


def stress_table(m: SSIModel) -> Dict[Key, np.ndarray]:
    """Element stresses attached to a model by READSTR / CSECT: ``(group, element) -> (6,)``."""
    st = m.ui_state.get("element_stress")
    if not isinstance(st, dict):
        return {}
    out = {}
    for row in st.get("rows", []):
        out[(int(row[0]), int(row[1]))] = np.array([float(x) for x in row[2:8]] + [0.0] * max(0, 8 - len(row)))
    return out


def stress_info(m: SSIModel) -> Dict[str, Any]:
    st = m.ui_state.get("element_stress")
    return st if isinstance(st, dict) else {}


def stress_mismatch(table: Dict[Key, Any], keys: Iterable[Key], where: str = "are loaded",
                    remedy: str = "a renumbering after READSTR (GCOM, MERGEGROUP, ECOMPR) leaves the records under "
                                  "the old numbers; READSTR the file again") -> Optional[str]:
    """Warning text when element stresses are available but **none** belongs to the section elements
    ``keys`` (None otherwise).

    Records are matched by (group, element) (spec 10 section 3.9).  The usual cause is a renumbering
    after READSTR or after STRESS wrote the frames (GCOM, MERGEGROUP, ECOMPR): the records stay under the
    old numbers, so without this message the resultants would silently be zero."""
    keys = set(keys)
    if not table or not keys or any(k in table for k in keys):
        return None
    sg = sorted({g for g, _ in table})
    kg = sorted({g for g, _ in keys})
    return (f"{len(table)} element stress records {where} but none matches an element of the section "
            f"(groups in the stress table: {_short(sg)}; groups of the section: {_short(kg)}): {remedy}")


# ======================================================================================
# Section pieces (spec 10 section 3.2)
# ======================================================================================
@dataclass
class Piece:
    """One section piece: the intersection of an element with the section plane.

    ``kind`` SOLID: ``pts`` (k, 3) polygon vertices on the plane (counter-clockwise about the plane
    normal used to build it; the area formulas take the absolute value); SHELL / PLANE: ``pts`` (2, 3)
    the segment where the mid-surface meets the plane, ``thick`` the shell thickness (PLANE 1, the
    unit thickness of plane strain), ``frame`` (3, 3) rows = SHELL axes x', y', z' (z' = m, the normal)
    in which its stresses are given; PLANE ``frame`` rows X, Z, Y.
    """
    group: int
    element: int
    kind: str
    pts: np.ndarray
    thick: float = 0.0
    frame: Optional[np.ndarray] = None

    @property
    def key(self) -> Key:
        return (self.group, self.element)

    def to_json(self) -> dict:
        d = dict(g=self.group, e=self.element, kind=self.kind, pts=[[float(x) for x in p] for p in self.pts],
                 t=float(self.thick))
        if self.frame is not None:
            d["frame"] = [[float(x) for x in r] for r in self.frame]
        return d

    @classmethod
    def from_json(cls, d: dict) -> "Piece":
        fr = d.get("frame")
        return cls(int(d["g"]), int(d["e"]), str(d["kind"]), np.asarray(d["pts"], float), float(d.get("t", 0.0)),
                   None if fr is None else np.asarray(fr, float))


@dataclass
class BeamCrossing:
    """A BEAMS element crossing the plane at ``q`` (CSECT plots it as a unit-length beam along n)."""
    group: int
    element: int
    q: np.ndarray
    k: Optional[np.ndarray] = None      # K (orientation) node position


@dataclass
class SectionPieces:
    pieces: List[Piece] = field(default_factory=list)
    beams: List[BeamCrossing] = field(default_factory=list)
    skipped: Dict[str, List[Key]] = field(default_factory=dict)   # reason -> elements
    coincident: int = 0                                           # faces / edges on the plane


def _dedupe_points(pts: List[np.ndarray], tol: float) -> np.ndarray:
    out: List[np.ndarray] = []
    for p in pts:
        if all(float(np.linalg.norm(p - q)) > tol for q in out):
            out.append(p)
    return np.array(out).reshape(-1, 3)


def order_polygon(pts: np.ndarray, ez: np.ndarray) -> np.ndarray:
    """Order points of a convex polygon counter-clockwise about ``ez`` (angles about their mean)."""
    if len(pts) < 3:
        return pts
    u, v = in_plane_basis(ez)
    c = pts.mean(axis=0)
    ang = np.arctan2((pts - c) @ v, (pts - c) @ u)
    return pts[np.argsort(ang, kind="stable")]


def polygon_area_3d(pts: np.ndarray, ez: np.ndarray) -> float:
    if len(pts) < 3:
        return 0.0
    u, v = in_plane_basis(ez)
    x, y = pts @ u, pts @ v
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y)))


def _project(pts: np.ndarray, P: np.ndarray, ez: np.ndarray) -> np.ndarray:
    pts = np.asarray(pts, float).reshape(-1, 3)
    return pts - np.outer((pts - P) @ ez, ez)


def _segment_from_polygon(X: np.ndarray, d: np.ndarray, tol: float) -> Optional[np.ndarray]:
    """Segment where a (closed) polygon with signed corner distances ``d`` meets the plane."""
    n = len(X)
    pts = [X[k] for k in range(n) if abs(d[k]) <= tol]
    for k in range(n):
        a, b = k, (k + 1) % n
        if (d[a] < -tol and d[b] > tol) or (d[a] > tol and d[b] < -tol):
            t = d[a] / (d[a] - d[b])
            pts.append(X[a] + t * (X[b] - X[a]))
    P = _dedupe_points(pts, tol)
    if len(P) < 2:
        return None
    if len(P) > 2:                       # warped facet: keep the two farthest points
        D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
        i, j = np.unravel_index(int(np.argmax(D)), D.shape)
        P = P[[i, j]]
    return P


def _shell_corners(nodes: Sequence[int]) -> List[int]:
    """Corner node ids of a SHELL / PLANE facet (3 or 4 distinct; L = K etc. give triangles)."""
    return list(dict.fromkeys(n for n in nodes[:4] if n))


def section_pieces(m: SSIModel, keys: Iterable[Key], point: Sequence[float], normal: Sequence[float],
                   v: Optional[ModelView] = None, notes: Optional[Notes] = None) -> SectionPieces:
    """Section pieces of the elements ``keys`` of ``m`` cut by the plane (point, normal), spec 10
    section 3.2.

    * SOLID: crossed when ``min d < -tol`` and ``max d > tol`` (signed node distances ``d``); polygon =
      on-plane nodes + edge intersections, ordered counter-clockwise about ``ez``;
    * SHELL / PLANE: segment of the mid-surface; a shell lying in the plane is skipped (warning); PLANE
      needs a plane normal in the X-Z plane;
    * a SOLID face or a SHELL / PLANE edge lying on the plane (D-SEC-04) is counted once, keyed by its
      node set, preferring the element on the -ez side;
    * BEAMS crossings are returned separately (CSECT plots them; they carry no ``.ess`` record);
      SPRING, GENERAL and TSHELL elements are skipped (listed).
    """
    from ..elements.base import ElementError
    from ..elements.shell import local_frame

    notes = _notes(notes)
    v = v or model_view(m)
    tol = v.tol
    ez = unit(normal, "plane normal")
    P = np.asarray(point, float).reshape(3)
    out = SectionPieces()
    skip = out.skipped
    # coincident face / edge candidates: key -> list of (side, piece)
    faces: Dict[frozenset, List[Tuple[int, Piece]]] = {}
    beam_faces: Dict[frozenset, List[Tuple[int, BeamCrossing]]] = {}
    for g, e in sorted(set((int(a), int(b)) for a, b in keys)):
        grp = m.groups.get(g)
        el = grp.elements.get(e) if grp is not None else None
        if el is None:
            skip.setdefault("not in the model", []).append((g, e))
            continue
        t = grp.type
        if t in (SPRING, GENERAL, TSHELL):
            skip.setdefault(f"{grp.type_name} (no section piece)", []).append((g, e))
            continue
        if t == SOLID:
            ns = list(el.nodes) + [0] * (8 - len(el.nodes))
            if any(n == 0 for n in ns[:8]) or not v.defined(ns[:8]):
                skip.setdefault("undefined nodes", []).append((g, e))
                continue
            X = np.array([v.P[n] for n in ns[:8]])
            d = (X - P) @ ez
            on = [k for k in range(8) if abs(d[k]) <= tol]
            if d.min() < -tol and d.max() > tol:
                pts = [X[k] for k in on]
                seen = set()
                for a, b in HEX_EDGES:
                    pair = tuple(sorted((ns[a], ns[b])))
                    if ns[a] == ns[b] or pair in seen:
                        continue
                    seen.add(pair)
                    if (d[a] < -tol and d[b] > tol) or (d[a] > tol and d[b] < -tol):
                        tt = d[a] / (d[a] - d[b])
                        pts.append(X[a] + tt * (X[b] - X[a]))
                poly = order_polygon(_dedupe_points([p for p in _project(np.array(pts), P, ez)], tol), ez)
                if len(poly) >= 3 and polygon_area_3d(poly, ez) > tol * tol:
                    out.pieces.append(Piece(g, e, "SOLID", poly))
                continue
            ids = frozenset(ns[k] for k in on)
            if len(ids) >= 3:
                poly = order_polygon(_dedupe_points(list(_project(X[on], P, ez)), tol), ez)
                if len(poly) >= 3 and polygon_area_3d(poly, ez) > tol * tol:
                    side = -1 if d.max() <= tol else 1
                    faces.setdefault(ids, []).append((side, Piece(g, e, "SOLID", poly)))
            continue
        if t in (SHELL, PLANE):
            corners = _shell_corners(el.nodes)
            if len(corners) < 3 or not v.defined(corners):
                skip.setdefault("undefined nodes", []).append((g, e))
                continue
            X = np.array([v.P[n] for n in corners])
            if t == SHELL:
                try:
                    Lam, _, _ = local_frame(X)
                except ElementError:
                    skip.setdefault("degenerate SHELL", []).append((g, e))
                    continue
                thick = float(el.thick)
                if not thick > 0:
                    skip.setdefault("SHELL without thickness (THICK)", []).append((g, e))
                    continue
                frame = Lam
            else:
                if abs(ez[1]) > PARALLEL_TOL:
                    skip.setdefault("PLANE (the plane normal must lie in the X-Z plane)", []).append((g, e))
                    continue
                thick = 1.0
                frame = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
            d = (X - P) @ ez
            on = [k for k in range(len(corners)) if abs(d[k]) <= tol]
            if len(on) >= 3:
                skip.setdefault(f"{grp.type_name} lying in the section plane", []).append((g, e))
                continue
            if d.min() < -tol and d.max() > tol:
                seg = _segment_from_polygon(X, d, tol)
                if seg is not None:
                    out.pieces.append(Piece(g, e, grp.type_name, _project(seg, P, ez), thick, frame))
                continue
            if len(on) == 2 and (abs(on[0] - on[1]) == 1 or abs(on[0] - on[1]) == len(corners) - 1):
                side = -1 if d.max() <= tol else 1
                pc = Piece(g, e, grp.type_name, _project(X[on], P, ez), thick, frame)
                faces.setdefault(frozenset(corners[k] for k in on), []).append((side, pc))
            continue
        if t == BEAMS:
            ns = list(el.nodes) + [0, 0, 0]
            I, J, K = ns[0], ns[1], ns[2]
            if not I or not J or not v.defined([I, J]):
                skip.setdefault("undefined nodes", []).append((g, e))
                continue
            XI, XJ = v.P[I], v.P[J]
            dI, dJ = float((XI - P) @ ez), float((XJ - P) @ ez)
            kpos = v.P[K] if K and K in v.P else None
            if abs(dI) <= tol and abs(dJ) <= tol:
                skip.setdefault("BEAMS lying in the section plane", []).append((g, e))
            elif (dI < -tol and dJ > tol) or (dI > tol and dJ < -tol):
                q = XI + dI / (dI - dJ) * (XJ - XI)
                out.beams.append(BeamCrossing(g, e, _project(q, P, ez)[0], kpos))
            elif abs(dI) <= tol or abs(dJ) <= tol:
                node, other = (I, dJ) if abs(dI) <= tol else (J, dI)
                bc = BeamCrossing(g, e, _project(v.P[node], P, ez)[0], kpos)
                beam_faces.setdefault(frozenset([node]), []).append((-1 if other < 0 else 1, bc))
            continue
        skip.setdefault(f"type {t}", []).append((g, e))
    # D-SEC-04: a face / edge on the plane counts once, the element on the -ez side preferred
    for key in sorted(faces, key=lambda s: sorted(s)):
        cands = sorted(faces[key], key=lambda sp: (sp[0], sp[1].group, sp[1].element))
        out.pieces.append(cands[0][1])
        out.coincident += 1
    for key in sorted(beam_faces, key=lambda s: sorted(s)):
        cands = sorted(beam_faces[key], key=lambda sp: (sp[0], sp[1].group, sp[1].element))
        out.beams.append(cands[0][1])
    out.pieces.sort(key=lambda p: (p.group, p.element))
    out.beams.sort(key=lambda b: (b.group, b.element))
    if out.coincident:
        notes.info(f"{out.coincident} element faces/edges lie on the section plane: each is counted once (the "
                   f"element on the -n side preferred, D-SEC-04); move the plane to the mid-height of an element "
                   f"row so that the element-centre stresses represent the section")
    return out


def report_skipped(sp: SectionPieces, notes: Notes, beams_are_skipped: bool = True) -> None:
    """Warnings listing the elements of the cut that give no section piece."""
    for why, keys in sorted(sp.skipped.items()):
        notes.warn(f"{len(keys)} elements skipped, {why}: {_short(_ge(keys))}")
    if beams_are_skipped and sp.beams:
        notes.warn(f"{len(sp.beams)} BEAMS elements cross the plane; they carry end forces (no .ess record) and "
                   f"are not section pieces: {_short(_ge((b.group, b.element) for b in sp.beams))}")


# ======================================================================================
# Section properties and resultants (requirements 4.14, spec 10 section 3.2)
# ======================================================================================
@dataclass
class SectionProps:
    """Section properties in the local cut axes about the area centroid ``C`` (D-SEC-01)."""
    area: float
    C: np.ndarray                     # (3,) global
    Ixx: float
    Iyy: float
    Ixy: float
    ex: np.ndarray
    ey: np.ndarray
    ez: np.ndarray
    A: np.ndarray                     # (nP,) piece areas
    c: np.ndarray                     # (nP, 3) piece centroids (global, on the plane)

    @property
    def Izz(self) -> float:
        return self.Ixx + self.Iyy

    def principal(self) -> Tuple[float, float, float]:
        """Principal second moments (I1 >= I2) and the angle (deg) of the I1 axis from local x'...
        measured about ez (I1 is the moment about that axis)."""
        J = np.array([[self.Ixx, -self.Ixy], [-self.Ixy, self.Iyy]])
        w, V = np.linalg.eigh(J)
        k = int(np.argmax(w))
        ang = math.degrees(math.atan2(V[1, k], V[0, k]))
        return float(w.max()), float(w.min()), ang

    def values(self) -> List[float]:
        return [self.area, float(self.C[0]), float(self.C[1]), float(self.C[2]), self.Ixx, self.Iyy, self.Ixy,
                self.Izz]


def _strip_geometry(p: Piece, ez: np.ndarray) -> Tuple[float, float, np.ndarray, np.ndarray]:
    """(L, w, d unit along the segment, c midpoint) of a SHELL / PLANE strip; ``w = t/|m x ez|``."""
    a, b = p.pts[0], p.pts[1]
    L = float(np.linalg.norm(b - a))
    d = (b - a) / L if L > 0 else np.zeros(3)
    m = p.frame[2] if p.frame is not None else np.array([0.0, 1.0, 0.0])
    s = float(np.linalg.norm(np.cross(m, ez)))
    w = p.thick / s if s > 1e-12 else 0.0
    return L, w, d, 0.5 * (a + b)


def section_properties(pieces: Sequence[Piece], point: Sequence[float], ex: np.ndarray, ey: np.ndarray,
                       ez: np.ndarray) -> SectionProps:
    """Area, centroid and second moments of the section (exact polygon formulas, rectangle strips)."""
    if not pieces:
        raise SectionError("the plane does not cut any SOLID, SHELL or PLANE element of the section")
    O = np.asarray(point, float).reshape(3)
    nP = len(pieces)
    A = np.zeros(nP)
    c2 = np.zeros((nP, 2))
    I0 = np.zeros((nP, 3))            # (int x^2, int y^2, int x y) about O, local axes
    for i, p in enumerate(pieces):
        if p.kind == "SOLID":
            xy = np.column_stack([(p.pts - O) @ ex, (p.pts - O) @ ey])
            x, y = xy[:, 0], xy[:, 1]
            x1, y1 = np.roll(x, -1), np.roll(y, -1)
            a = x * y1 - x1 * y
            Asg = 0.5 * float(a.sum())
            if Asg == 0.0:
                continue
            s = 1.0 if Asg > 0 else -1.0
            A[i] = abs(Asg)
            c2[i] = [float(((x + x1) * a).sum()) / (6 * Asg), float(((y + y1) * a).sum()) / (6 * Asg)]
            I0[i] = [s * float(((x * x + x * x1 + x1 * x1) * a).sum()) / 12.0,
                     s * float(((y * y + y * y1 + y1 * y1) * a).sum()) / 12.0,
                     s * float(((x * y1 + 2 * x * y + 2 * x1 * y1 + x1 * y) * a).sum()) / 24.0]
        else:
            L, w, d, cm = _strip_geometry(p, ez)
            A[i] = L * w
            c2[i] = [float((cm - O) @ ex), float((cm - O) @ ey)]
            d2 = np.array([float(d @ ex), float(d @ ey)])
            nu2 = np.array([-d2[1], d2[0]])
            J = A[i] / 12.0 * (L * L * np.outer(d2, d2) + w * w * np.outer(nu2, nu2))
            I0[i] = [J[0, 0] + A[i] * c2[i, 0] ** 2, J[1, 1] + A[i] * c2[i, 1] ** 2,
                     J[0, 1] + A[i] * c2[i, 0] * c2[i, 1]]
    At = float(A.sum())
    if not At > 0:
        raise SectionError("the section has zero area")
    Cx, Cy = float(A @ c2[:, 0]) / At, float(A @ c2[:, 1]) / At
    Iyy = float(I0[:, 0].sum()) - At * Cx * Cx
    Ixx = float(I0[:, 1].sum()) - At * Cy * Cy
    Ixy = float(I0[:, 2].sum()) - At * Cx * Cy
    C = O + Cx * ex + Cy * ey
    cg = O[None, :] + c2[:, 0:1] * ex[None, :] + c2[:, 1:2] * ey[None, :]
    return SectionProps(At, C, Ixx, Iyy, Ixy, ex, ey, ez, A, cg)


def _skew(a: np.ndarray) -> np.ndarray:
    return np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])


def traction_operator(p: Piece, ez: np.ndarray) -> np.ndarray:
    """``G (3, 6)`` with ``S . ez = G s`` for the six stored stress components ``s`` of the piece.

    SOLID ``s = (Sxx, Syy, Szz, Sxy, Sxz, Syz)`` global; SHELL ``s = (Sx'x', Sy'y', Sx'y', ...)`` with the
    plane-stress tensor ``S = Sx'x' e1e1 + Sy'y' e2e2 + Sx'y' (e1e2 + e2e1)``; PLANE
    ``s = (Sxx, Szz, Txz, ...)`` in the X-Z plane.  The SHELL membrane columns of the ``.ess`` frames
    are *stresses* (D-STR-10, D-W1-14), so a strip of area ``A = L t/|m x ez|`` carries
    ``A (S . ez) = L t (S . nu_s)``: exact for oblique cuts.
    """
    G = np.zeros((3, 6))
    if p.kind == "SOLID":
        G[0] = [ez[0], 0, 0, ez[1], ez[2], 0]
        G[1] = [0, ez[1], 0, ez[0], 0, ez[2]]
        G[2] = [0, 0, ez[2], 0, ez[0], ez[1]]
    elif p.kind == "SHELL":
        e1, e2 = p.frame[0], p.frame[1]
        a1, a2 = float(e1 @ ez), float(e2 @ ez)
        G[:, 0] = e1 * a1
        G[:, 1] = e2 * a2
        G[:, 2] = e1 * a2 + e2 * a1
    else:                                         # PLANE: X-Z plane strain
        G[:, 0] = [ez[0], 0, 0]
        G[:, 1] = [0, 0, ez[2]]
        G[:, 2] = [ez[2], 0, ez[0]]
    return G


def bending_operator(p: Piece, ez: np.ndarray) -> np.ndarray:
    """``B (3, 6)``: the SHELL plate-bending couple ``mb = L m x (Mb . nu_s)`` (D-SEC-03; columns 3-5 =
    ``Mx'x' My'y' Mx'y'``; ``nu_s`` = in-shell unit normal of the cut line with ``nu_s . ez > 0``).

    With ``Mx'x' = int sigma_x'x' z' dz'`` the bending stresses over the thickness of the cut face exert
    the couple ``int (z' m) x (sigma . nu_s) dz' = m x (Mb . nu_s)`` per unit length of the cut line.
    """
    B = np.zeros((3, 6))
    if p.kind != "SHELL":
        return B
    e1, e2, m = p.frame[0], p.frame[1], p.frame[2]
    nu = ez - float(ez @ m) * m
    nn = float(np.linalg.norm(nu))
    if nn <= 1e-12:
        return B
    nu /= nn
    L = float(np.linalg.norm(p.pts[1] - p.pts[0]))
    b1, b2 = float(e1 @ nu), float(e2 @ nu)
    m1, m2 = np.cross(m, e1), np.cross(m, e2)
    B[:, 3] = L * b1 * m1
    B[:, 4] = L * b2 * m2
    B[:, 5] = L * (b2 * m1 + b1 * m2)
    return B


def resultant_operators(pieces: Sequence[Piece], props: SectionProps, bending: bool = True) -> np.ndarray:
    """``T (nP, 6, 6)`` with local ``(Fx, Fy, Fz, Mx, My, Mz) = sum_i T_i s_i`` (requirements 4.14)."""
    E = np.vstack([props.ex, props.ey, props.ez])
    T = np.zeros((len(pieces), 6, 6))
    for i, p in enumerate(pieces):
        AG = props.A[i] * traction_operator(p, props.ez)
        Mg = _skew(props.c[i] - props.C) @ AG
        if bending:
            Mg = Mg + bending_operator(p, props.ez)
        T[i, :3] = E @ AG
        T[i, 3:] = E @ Mg
    return T


def gather_stresses(pieces: Sequence[Piece], table: Dict[Key, np.ndarray],
                    parents: Optional[Dict[Key, Key]] = None) -> Tuple[np.ndarray, np.ndarray]:
    """(S (nP, 6), available (nP,) bool) of the pieces from a stress table."""
    S = np.zeros((len(pieces), 6))
    ok = np.zeros(len(pieces), dtype=bool)
    for i, p in enumerate(pieces):
        k = p.key
        val = table.get(k)
        if val is None and parents:
            val = table.get(parents.get(k, k))
        if val is not None:
            S[i] = val[:6]
            ok[i] = True
    return S, ok


def resultants(T: np.ndarray, S: np.ndarray, ok: Optional[np.ndarray] = None) -> np.ndarray:
    """Local resultants ``(6,)`` = ``sum_i T_i s_i`` over the pieces with data."""
    if ok is not None:
        S = np.where(ok[:, None], S, 0.0)
    return np.einsum("pij,pj->i", T, S)


def bending_on(m: SSIModel) -> bool:
    """D-SEC-03: SHELL plate bending enters the resultants unless ``EDUOPT,SECTBEND,0``."""
    return eduopt(m, "SECTBEND") != "0"


# ======================================================================================
# Local cut system (spec 10 section 1.4)
# ======================================================================================
def store_cut_system(m: SSIModel, sysno: int, origin: np.ndarray, ex: np.ndarray, ey: np.ndarray,
                     ez: np.ndarray, notes: Optional[Notes] = None) -> bool:
    """Store the local cut axes as Cartesian system ``sysno`` (not activated).  ``sysno <= 0`` stores
    nothing (information).  Two systems are never redefined (warning, nothing stored): one that
    already holds nodes (they would move) and the **active** system of CSYS, even without nodes yet,
    because later N commands are placed in it (spec 10 section 1.4: storing must not change the active
    coordinate system for later N commands)."""
    notes = _notes(notes)
    if sysno <= 0:
        notes.info("<sysno> 0 or blank: the local cut system is not stored")
        return False
    if int(sysno) == int(m.csys_active):
        notes.warn(f"system {sysno} is the active coordinate system (CSYS): not redefined (later N commands would "
                   f"be placed in the cut system); choose another <sysno>")
        return False
    held = m.nodes_in_system(int(sysno))
    if held:
        notes.warn(f"system {sysno} holds {len(held)} nodes: not redefined (its nodes would move); choose another "
                   f"<sysno>")
        return False
    R = np.column_stack([ex, ey, ez])
    txy, tyz, txz = euler_from_matrix(R)
    params = [float(x) + 0.0 for x in (origin[0], origin[1], origin[2], txy, tyz, txz)]   # + 0.0: no -0
    if int(sysno) in m.csys:
        notes.info(f"system {sysno} redefined as the local cut system")
    m.csys[int(sysno)] = CoordSys(int(sysno), origin=np.array(params[:3]), R=loc_matrix(txy, tyz, txz),
                                  kind="LOC", type=0, params=params)
    return True


# ======================================================================================
# Section definition of the active model for CALCPAR / CALCMOI
# ======================================================================================
@dataclass
class ActiveSection:
    """The section CALCPAR / CALCMOI integrate.  ``how``: ``"csect"`` (stored CSECT pieces),
    ``"extruded"`` (mid-plane of a unit-thickness extrusion) or ``"base"`` (base of the whole
    structure, Q-11); ``skipped`` lists the BEAMS and the elements without a piece (spec 10 section 3.2
    warning)."""
    pieces: List[Piece]
    point: np.ndarray
    from_csect: bool
    parents: Dict[Key, Key] = field(default_factory=dict)
    skipped: Optional[SectionPieces] = None
    how: str = "csect"


def _fmt_vec(a: Sequence[float]) -> str:
    return ", ".join(f"{float(x):g}" for x in a)


def is_csect_model(m: SSIModel) -> bool:
    """True when the model carries the section data of CSECT (``ui_state['csect']``)."""
    sec = m.ui_state.get("csect")
    return isinstance(sec, dict) and sec.get("pieces") is not None


def dof_offsets(v: ModelView, refs: Iterable[Any], ez: np.ndarray) -> np.ndarray:
    """Signed offsets ``n.x`` of the distinct nodes carrying DOFs of the elements ``refs``
    (:class:`sassi.prep.check.ElemRef`).  BEAMS / GENERAL K nodes only orient the element and nodes
    that belong to no element (e.g. lumped masses) say nothing about the section, so both are left out."""
    nodes = sorted({n for r in refs for n in r.dof_nodes if n in v.P})
    if not nodes:
        return np.zeros(0)
    X = np.array([v.P[n] for n in nodes])
    return X[np.all(np.isfinite(X), axis=1)] @ np.asarray(ez, float)


def unit_extrusion_midplane(d: np.ndarray, tol: float) -> Optional[float]:
    """Offset of the mid-plane when every offset ``d`` lies on one of two planes one unit apart (the
    +-0.5 extrusion of a CSECT model, spec 10 section 3.3); None otherwise."""
    if not len(d):
        return None
    lo, hi = float(d.min()), float(d.max())
    if abs(hi - lo - 1.0) > tol:
        return None
    if not np.all((np.abs(d - lo) <= tol) | (np.abs(d - hi) <= tol)):
        return None
    return 0.5 * (lo + hi) + 0.0


def extrusion_normal(m: SSIModel, v: Optional[ModelView] = None) -> Optional[np.ndarray]:
    """Unit normal ``n`` of a cross-section model without CSECT data (e.g. read back by WRITE -> INP),
    or None.

    CSECT extrudes every piece by the same vector ``n`` (SOLID: node k -> node k+4; SHELL / PLANE:
    node 1 -> node 4; BEAMS: I -> J), so the candidate is taken from the first element and accepted
    when it has unit length and every DOF node of the model lies on one of the two planes +-0.5 about
    the section (:func:`unit_extrusion_midplane`)."""
    v = v or model_view(m)
    for r in v.elems:
        ns = list(r.elem.nodes)
        if r.type == SOLID and len(ns) >= 8:
            a, b = ns[0], ns[4]
        elif r.type in (SHELL, PLANE) and len(ns) >= 4:
            a, b = ns[0], ns[3]
        elif r.type == BEAMS and len(ns) >= 2:
            a, b = ns[0], ns[1]
        else:
            return None
        if not a or not b or not v.defined([a, b]):
            return None
        w = v.P[b] - v.P[a]
        if not np.all(np.isfinite(w)) or abs(float(np.linalg.norm(w)) - 1.0) > v.tol:
            return None
        n = w / float(np.linalg.norm(w))
        return n if unit_extrusion_midplane(dof_offsets(v, v.elems, n), v.tol) is not None else None
    return None


def _stored_skips(m: SSIModel, sec: dict) -> SectionPieces:
    """BEAMS crossings and skipped elements recorded by CSECT (BEAMS deleted since are dropped)."""
    out = SectionPieces()
    for d in sec.get("beams", []):
        g, e = int(d["g"]), int(d["e"])
        if g in m.groups and e in m.groups[g].elements:
            out.beams.append(BeamCrossing(g, e, np.asarray(d.get("q", [np.nan] * 3), float)))
    for why, keys in sorted(sec.get("skipped", {}).items()):
        out.skipped[str(why)] = [(int(a), int(b)) for a, b in keys]
    return out


def active_section(m: SSIModel, normal: Sequence[float], notes: Optional[Notes] = None) -> ActiveSection:
    """Pieces of the active model for CALCPAR / CALCMOI (spec 10 sections 3.7-3.8; see the module
    docstring for the three cases).

    1. a CSECT model carries its exact pieces (``ui_state['csect']``); the normal must be parallel to
       the CSECT normal (either sense: the other sense views the section from the other side);
    2. a unit-thickness extrusion along ``n`` (a CSECT model after WRITE -> INP): every element is
       intersected with the mid-plane of the extrusion;
    3. otherwise the whole structure (excavated soil left out) is cut by the plane through its end on
       the -n side (Q-11: base forces and moments of a building); an error when that plane carries no
       section piece (n oblique to the base, or only BEAMS / SPRING elements reach it).
    """
    notes = _notes(notes)
    ez = unit(normal, "plane normal")
    sec = m.ui_state.get("csect")
    if is_csect_model(m):
        n0 = unit(sec["normal"], "stored CSECT normal")
        if float(np.linalg.norm(np.cross(n0, ez))) > PARALLEL_TOL:
            raise SectionError(f"the normal must be the CSECT plane normal ({_fmt_vec(n0)}) of this "
                               f"cross-section model")
        pieces = [Piece.from_json(d) for d in sec["pieces"]]
        present = [p for p in pieces if p.group in m.groups and p.element in m.groups[p.group].elements]
        if len(present) < len(pieces):
            notes.warn(f"{len(pieces) - len(present)} section pieces of deleted elements ignored")
        parents = {(int(a), int(b)): (int(a), int(c)) for a, b, c in sec.get("parents", [])}
        return ActiveSection(present, np.asarray(sec["point"], float), True, parents, _stored_skips(m, sec), "csect")
    v = model_view(m)
    model_extents(v)                                  # raises when the model has no defined node
    mid = unit_extrusion_midplane(dof_offsets(v, v.elems, ez), v.tol)
    if mid is not None:
        P = mid * ez
        notes.info(f"section plane: normal ({_fmt_vec(ez)}) through the mid-plane n.x = {mid:.6g} of the "
                   f"unit-thickness extrusion along it (cross-section model without CSECT data, spec 10 section 3.7)")
        sp = section_pieces(m, [(r.group, r.id) for r in v.elems], P, ez, v=v, notes=notes)
        return ActiveSection(sp.pieces, P, False, {}, sp, "extruded")
    # (A model extruded along another direction is not rejected: one layer of unit-thick elements is
    # geometrically indistinguishable from a CSECT extrusion, so it is treated as a whole structure.)
    # Q-11: the whole structure, cut at its end on the -n side (the base when n points up)
    struct = [r for r in v.elems if not r.excavated]
    use = struct or v.elems
    d = dof_offsets(v, use, ez)
    if not len(d):
        raise SectionError("the active model has no element with defined nodes")
    base = float(d.min()) + 0.0                     # + 0.0: no -0 in the message
    P = base * ez
    notes.info(f"section plane: normal ({_fmt_vec(ez)}) through the end of the structure on the -n side (its base "
               f"when n points up), n.x = {base:.6g}: the active model is not a cross-section, so the resultants "
               f"are the actions of the whole structure on its support (spec 10 section 3.8, Q-11); for any other "
               f"plane build a cut and CSECT")
    if struct and len(struct) < len(v.elems):
        notes.info(f"{len(v.elems) - len(struct)} excavated-soil elements (ETYPE 2) are not part of the structure: "
                   f"left out of the section")
    sp = section_pieces(m, [(r.group, r.id) for r in use], P, ez, v=v, notes=notes)
    if not sp.pieces:
        report_skipped(sp, notes)
        raise SectionError(f"no SOLID face or SHELL / PLANE edge of the structure lies on its base plane "
                           f"n.x = {base:.6g} (n is oblique to the base, or only BEAMS / SPRING elements reach "
                           f"it): CALCPAR / CALCMOI on a model that is not a cross-section give the base section "
                           f"(Q-11); for another plane build a cut and CSECT, then ACTM to the cross-section model")
    return ActiveSection(sp.pieces, P, False, {}, sp, "base")


# ======================================================================================
# CSECT (spec 09 section 5.2, spec 04 section 8.3)
# ======================================================================================
@dataclass
class CsectResult:
    model: SSIModel
    pieces: SectionPieces
    elements: int
    sub_elements: int
    nodes: int


def _copy_tables(dst: SSIModel, src: SSIModel) -> None:
    for k, s in src.materials.items():
        dst.materials[k] = Material(k, s.val1, s.val2, s.weight, s.pdamp, s.sdamp, s.mtype)
    for k, s in src.layers.items():
        dst.layers[k] = SoilLayer(k, s.thick, s.weight, s.vp, s.vs, s.pdamp, s.sdamp)
    for k, s in src.sections.items():
        dst.sections[k] = BeamSection(k, s.axial, s.shear2, s.shear3, s.tors, s.flex2, s.flex3)
    for k, s in src.springs.items():
        dst.springs[k] = SpringProp(k, s.scx, s.scy, s.scz, s.scxx, s.scyy, s.sczz, s.damp)
    for k, s in src.matrices.items():
        p = MatrixProp.from_json(s.to_json())
        dst.matrices[k] = p


def copy_settings(dst: SSIModel, src: SSIModel) -> None:
    """Model-wide settings a submodel needs to be interpreted like its source: HOUSE gravity and ground
    elevation (shared variables; weight -> mass, ETYPE 0 classification), MOPT, TOPL and the EDUOPT
    switches (e.g. GEOMTOL, SECTBEND)."""
    rec = dst.options.ensure_record("HOUSE")
    rec.set_arg(1, src.gravity)
    rec.set_arg(2, src.ground_elevation)
    mo = src.options.record("MOPT")
    if mo is not None:
        dst.options.set_record(mo.copy())
    for k, r in src.options.entries("EDUOPT"):
        dst.options.set_entry("EDUOPT", k, r.copy())
    dst.topl = list(src.topl)


def csect_model(src: SSIModel, keys: Iterable[Key], point: Sequence[float], normal: Sequence[float],
                notes: Optional[Notes] = None, source_model: int = -1, cut: int = 0) -> CsectResult:
    """Cross-section model of the elements ``keys`` (CSECT; spec 09 section 5.2, spec 04 section 8.3).

    Every section piece is extruded **+-0.5 along n** (unit thickness, centred on the plane, spec 10
    section 3.3) and keeps the group and element number of its element: a SOLID polygon of 3 / 4
    vertices becomes one prism / hexahedron, larger polygons a fan of hexahedra and prisms whose first
    sub-element keeps the element number and the others take numbers after the group's largest element
    number; a SHELL strip becomes a 4-node SHELL (A - n/2, B - n/2, B + n/2, A + n/2) of thickness
    ``t_eff = t/|m x n|``; a PLANE strip a 4-node PLANE (unit plane-strain thickness); a BEAMS crossing a
    unit-length beam along n (a node at the K position keeps the orientation).  Coincident nodes are
    merged; node numbers are new.  Materials, properties, settings, the resolved ETYPE and the element
    stresses (READSTR before CSECT) are copied; the exact pieces, the crossing BEAMS and the cut
    elements without a piece are kept in ``ui_state['csect']`` for CALCPAR / CALCMOI (which repeat the
    spec 10 section 3.2 warning listing them) and CALCC (which leaves the BEAMS out).

    BEAMS are drawn but are not section pieces: they carry end forces, not ``.ess`` stresses, so the
    section properties and resultants never include them (an information line lists them).
    """
    notes = _notes(notes)
    v = model_view(src)
    ez = unit(normal, "plane normal")
    P = np.asarray(point, float).reshape(3)
    keys = sorted(set(keys))
    sp = section_pieces(src, keys, P, ez, v=v, notes=notes)
    report_skipped(sp, notes, beams_are_skipped=False)
    if not sp.pieces and not sp.beams:
        raise SectionError("the plane does not cut any element of the cut")
    if sp.beams:
        notes.info(f"{len(sp.beams)} BEAMS elements cross the plane: drawn as unit-length beams along n, but they "
                   f"carry end forces (no .ess record) and are not part of the section properties or resultants "
                   f"(CALCPAR, CALCMOI) nor of CALCC: {_short(_ge((b.group, b.element) for b in sp.beams))}")
    etypes = {(r.group, r.id): r.etype for r in v.elems}
    out = SSIModel()
    copy_settings(out, src)
    _copy_tables(out, src)
    pts: List[np.ndarray] = []                  # new node coordinates (merged afterwards)

    def node(x: np.ndarray) -> int:
        pts.append(np.asarray(x, float))
        return len(pts) - 1                     # provisional index

    elems: List[Tuple[int, int, int, List[int], Element]] = []      # (group, number, type, provisional nodes, src)
    next_id: Dict[int, int] = {}
    parents: List[List[int]] = []
    h = 0.5 * ez
    nsub = 0
    for p in sp.pieces:
        g, e = p.group, p.element
        grp = src.groups[g]
        el = grp.elements[e]
        if p.kind == "SOLID":
            poly = p.pts
            k = len(poly)
            parts: List[List[int]] = []
            i = 1
            while i + 2 <= k - 1:
                parts.append([0, i, i + 1, i + 2])
                i += 2
            if i + 1 == k - 1:
                parts.append([0, i, i + 1])
            for j, part in enumerate(parts):
                bot = [node(poly[q] - h) for q in part]
                top = [node(poly[q] + h) for q in part]
                if len(part) == 3:
                    bot.append(bot[2])
                    top.append(top[2])
                if j == 0:
                    num = e
                else:
                    num = next_id.get(g, max(grp.elements) + 1)
                    next_id[g] = num + 1
                    parents.append([g, num, e])
                    nsub += 1
                elems.append((g, num, SOLID, bot + top, el))
        else:
            a, b = p.pts[0], p.pts[1]
            elems.append((g, e, grp.type, [node(a - h), node(b - h), node(b + h), node(a + h)], el))
    for bc in sp.beams:
        grp = src.groups[bc.group]
        el = grp.elements[bc.element]
        ns = [node(bc.q - h), node(bc.q + h)]
        if bc.k is not None:
            ns.append(node(bc.k))
        elems.append((bc.group, bc.element, BEAMS, ns, el))
    X = np.array(pts)
    lab = unique_points(X, v.tol)
    first: Dict[int, int] = {}
    for i, l in enumerate(lab):
        first.setdefault(int(l), i)
    for l in sorted(first):
        out.define_node(l + 1, tuple(float(c) for c in X[first[l]]), 0)
    shell_frames = {p.key: p.frame for p in sp.pieces if p.kind == "SHELL"}
    for g, num, t, prov, el in elems:
        grp = out.groups.get(g)
        if grp is None:
            sg = src.groups[g]
            grp = Group(g, t, sg.title)
            out.groups[g] = grp
        ne = el.copy(new_id=num, nodes=[int(lab[i]) + 1 for i in prov])
        if t in (SOLID, PLANE, SHELL):
            ne.etype = etypes.get((g, el.id), ne.etype)
        if t == SHELL:
            s = float(np.linalg.norm(np.cross(shell_frames[(g, el.id)][2], ez)))
            ne.thick = el.thick / s if s > 0 else el.thick
        grp.elements[num] = ne
    out.group_active = min(out.groups) if out.groups else None
    # element stresses (READSTR before CSECT) and the exact pieces
    table = stress_table(src)
    if table:
        rows = {}
        for g, num, t, prov, el in elems:
            val = table.get((g, el.id))
            if val is not None:
                rows[(g, num)] = val
        info = stress_info(src)
        out.ui_state["element_stress"] = {
            "file": info.get("file", ""), "maxima": bool(info.get("maxima", False)),
            "groups": [[g, ELEMENT_TYPE_NAMES.get(out.groups[g].type, "")] for g in sorted({k[0] for k in rows})],
            "rows": [[g, e] + [float(x) for x in rows[(g, e)]] for g, e in sorted(rows)]}
        mismatch = stress_mismatch(table, [p.key for p in sp.pieces])
        if mismatch:
            mismatch += ", then CSECT again"
            notes.warn(mismatch)
            out.ui_state["element_stress"]["unmatched"] = mismatch     # CALCPAR repeats it
    else:
        notes.info("no element stresses are loaded on the source model (READSTR before CSECT): CALCPAR will report "
                   "the section properties with zero resultants")
    out.ui_state["csect"] = {"source_model": int(source_model), "cut": int(cut), "point": [float(x) for x in P],
                             "normal": [float(x) for x in ez], "pieces": [p.to_json() for p in sp.pieces],
                             "parents": parents,
                             "beams": [{"g": b.group, "e": b.element, "q": [float(x) for x in b.q]} for b in sp.beams],
                             "skipped": {why: [[g, e] for g, e in ks] for why, ks in sorted(sp.skipped.items())}}
    if nsub:
        notes.info(f"{nsub} additional sub-elements (polygons with more than 4 vertices) numbered after the "
                   f"largest element number of their group: " + _short([f"{g}/{n} of {g}/{e}" for g, n, e in parents]))
    return CsectResult(out, sp, len(elems), nsub, len(first))


# ======================================================================================
# Submodels: CUT2SUB, EXTRACTEXCAV, TRANELEM / TRANVOL, SPLITGROUP
# ======================================================================================
def _copy_node(n: Node) -> Node:
    return Node(n.id, n.x, n.y, n.z, n.csys, list(n.fix), set(n.flags))


def submodel(src: SSIModel, keys: Iterable[Key], solid_shells: bool = False,
             notes: Optional[Notes] = None) -> SSIModel:
    """CUT2SUB (spec 09 section 5.3, D-MDL-13): the elements ``keys`` with their group and element
    numbers, every node they reference (fixities, interaction flags, masses and mass units), the
    coordinate systems, the referenced M / L / R / SC / MX entries, TOPL and its layers, and the model
    settings (:func:`copy_settings`).  ``solid_shells``: every SHELL becomes an 8-node SOLID extruded
    +-t/2 along its normal (thick-shell solids, **for plotting only**: there is no check that the
    solids do not intersect)."""
    from ..elements.base import ElementError
    from ..elements.shell import local_frame

    notes = _notes(notes)
    v = model_view(src)
    etypes = {(r.group, r.id): r.etype for r in v.elems}
    keys, miss = model_elements(src, keys)
    if miss:
        notes.warn(f"{len(miss)} elements of the cut are not in the active model: {_short(_ge(miss))}")
    if not keys:
        raise SectionError("no element of the cut is in the active model")
    out = SSIModel(title=src.title)
    copy_settings(out, src)
    for sid in sorted(src.csys):
        out.csys[sid] = CoordSys.from_json(src.csys[sid].to_json())
    next_node = max(src.nodes) + 1 if src.nodes else 1
    new_nodes: Dict[int, Tuple[float, float, float]] = {}
    converted = 0
    for g, e in keys:
        sg = src.groups[g]
        el = sg.elements[e]
        t = sg.type
        if solid_shells and t == SHELL:
            corners = _shell_corners(el.nodes)
            ok = len(corners) >= 3 and v.defined(corners) and el.thick > 0
            if ok:
                X = np.array([v.P[n] for n in corners])
                try:
                    Lam, _, _ = local_frame(X)
                except ElementError:
                    ok = False
            if not ok:
                notes.warn(f"SHELL {g}/{e} not converted (undefined nodes, no thickness or degenerate): not copied")
                continue
            else:
                mvec = 0.5 * el.thick * Lam[2]
                bot, top = [], []
                for x in X:
                    bot.append(next_node)
                    new_nodes[next_node] = tuple(float(c) for c in x - mvec)
                    top.append(next_node + 1)
                    new_nodes[next_node + 1] = tuple(float(c) for c in x + mvec)
                    next_node += 2
                if len(bot) == 3:
                    bot.append(bot[2])
                    top.append(top[2])
                grp = out.groups.get(g)
                if grp is None:
                    grp = Group(g, SOLID, sg.title)
                    out.groups[g] = grp
                grp.elements[e] = Element(e, bot + top, mat=el.mat, prop=el.prop, etype=1)
                converted += 1
                continue
        grp = out.groups.get(g)
        if grp is None:
            grp = Group(g, t, sg.title)
            out.groups[g] = grp
        grp.elements[e] = el.copy()
    # nodes referenced by the final element set
    used = sorted({n for gr in out.groups.values() for el in gr.elements.values() for n in el.nodes if n})
    for n in used:
        if n in src.nodes:
            out.nodes[n] = _copy_node(src.nodes[n])
            out.node_history.touch(n)
        elif n in new_nodes:
            out.define_node(n, new_nodes[n], 0)
    dangling = [n for n in used if n not in out.nodes]
    if dangling:
        notes.warn(f"{len(dangling)} referenced nodes are not defined in the source model: {_short(dangling)}")
    for n in used:
        if n in src.tmass:
            out.tmass[n] = list(src.tmass[n])
        if n in src.rmass:
            out.rmass[n] = list(src.rmass[n])
        if n in src.mass_units:
            out.mass_units[n] = src.mass_units[n]
    for name, hist in (("tmass", "tmass_history"), ("rmass", "rmass_history")):
        for n in getattr(src, hist):
            if n in getattr(out, name):
                getattr(out, hist).touch(n)
    _copy_referenced_tables(out, src, etypes)
    out.group_active = min(out.groups) if out.groups else None
    if converted:
        notes.info(f"{converted} SHELL elements converted to thick-shell SOLIDs (plotting only, no intersection "
                   f"check)")
    return out


def _copy_referenced_tables(out: SSIModel, src: SSIModel, etypes: Dict[Key, int]) -> None:
    """Copy the M / L / R / SC / MX entries referenced by the elements of ``out`` (same numbers) and the
    soil layers of TOPL."""
    mats: Set[int] = set()
    lays: Set[int] = set(int(x) for x in out.topl if x)
    secs: Set[int] = set()
    sprs: Set[int] = set()
    mxs: Set[int] = set()
    for g, el in out.iter_elements():
        t = g.type
        if t in (SOLID, PLANE):
            et = etypes.get((g.id, el.id), el.etype)
            if et == 2 or el.etype == 2:
                lays.add(el.mat)
            else:
                mats.add(el.mat)
        elif t in (SHELL, TSHELL, BEAMS):
            mats.add(el.mat)
        if t == BEAMS:
            secs.add(el.prop)
        elif t == SPRING:
            sprs.add(el.prop)
        elif t == GENERAL:
            mxs.add(el.prop)
    for k in sorted(mats):
        if k in src.materials:
            s = src.materials[k]
            out.materials[k] = Material(k, s.val1, s.val2, s.weight, s.pdamp, s.sdamp, s.mtype)
    for k in sorted(lays):
        if k in src.layers:
            s = src.layers[k]
            out.layers[k] = SoilLayer(k, s.thick, s.weight, s.vp, s.vs, s.pdamp, s.sdamp)
    for k in sorted(secs):
        if k in src.sections:
            s = src.sections[k]
            out.sections[k] = BeamSection(k, s.axial, s.shear2, s.shear3, s.tors, s.flex2, s.flex3)
    for k in sorted(sprs):
        if k in src.springs:
            s = src.springs[k]
            out.springs[k] = SpringProp(k, s.scx, s.scy, s.scz, s.scxx, s.scyy, s.sczz, s.damp)
    for k in sorted(mxs):
        if k in src.matrices:
            out.matrices[k] = MatrixProp.from_json(src.matrices[k].to_json())


def extract_excavation(src: SSIModel, notes: Optional[Notes] = None) -> Tuple[SSIModel, int]:
    """EXTRACTEXCAV (spec 09 section 5.8): every excavation element (SOLID / PLANE with ETYPE 2) with its
    nodes (interaction flags, fixities), the soil layers, TOPL and the settings, same group and element
    numbers.  The ETYPE must be explicit (spec 09 section 10): SOLID / PLANE elements with the implicit
    ETYPE 0 below grade (which D-HOU-01 would classify as excavated) are an error -- run ETYPEGEN."""
    notes = _notes(notes)
    v = model_view(src)
    implicit = [(r.group, r.id) for r in v.elems if r.type in (SOLID, PLANE) and r.elem.etype == 0 and r.etype == 2]
    if implicit:
        raise SectionError(f"{len(implicit)} SOLID/PLANE elements below grade have the implicit ETYPE 0 "
                           f"({_short(_ge(implicit))}): set the ETYPE explicitly (ETYPEGEN or ETYPE)")
    keys = [(r.group, r.id) for r in v.elems if r.type in (SOLID, PLANE) and r.elem.etype == 2]
    if not keys:
        raise SectionError("the active model has no excavation element (SOLID/PLANE with ETYPE 2)")
    out = submodel(src, keys, notes=notes)
    for k in sorted(src.layers):                     # every far-field layer (the L table is model-wide)
        s = src.layers[k]
        out.layers[k] = SoilLayer(k, s.thick, s.weight, s.vp, s.vs, s.pdamp, s.sdamp)
    return out, len(keys)


@dataclass
class TransferResult:
    elements: int
    nodes: int
    groups_changed: List[int]
    tables_added: int


def transfer_elements(src: SSIModel, dest: SSIModel, keys: Iterable[Key],
                      notes: Optional[Notes] = None) -> TransferResult:
    """TRANELEM / TRANVOL (spec 09 sections 5.11-5.12): copy elements of ``src`` into ``dest``.

    Overwrite semantics of the manual: destination elements with the same numbers are replaced, the
    destination group's type becomes the source type (warning), every referenced node is overwritten
    (global coordinates, stored in system 0; fixities and INT flags copied).  Referenced M / L / R / SC /
    MX entries missing in ``dest`` are added; existing entries are kept (warning when they differ:
    the destination's table governs, spec 09 open question).
    """
    notes = _notes(notes)
    v = model_view(src)
    etypes = {(r.group, r.id): r.etype for r in v.elems}
    keys, miss = model_elements(src, keys)
    if miss:
        notes.warn(f"{len(miss)} elements are not in the active model: {_short(_ge(miss))}")
    if not keys:
        raise SectionError("no element to transfer")
    nodes = sorted({n for g, e in keys for n in src.groups[g].elements[e].nodes if n})
    undefined = [n for n in nodes if n not in v.P or not np.all(np.isfinite(v.P[n]))]
    if undefined:
        raise SectionError(f"elements reference undefined nodes: {_short(undefined)}")
    changed = []
    for g in sorted({g for g, _ in keys}):
        sg = src.groups[g]
        dg = dest.groups.get(g)
        if dg is None:
            dest.groups[g] = Group(g, sg.type, sg.title)
        elif dg.type != sg.type:
            notes.warn(f"group {g} of the destination changed from {dg.type_name} to {sg.type_name} "
                       f"({len(dg.elements)} elements already in it keep their data)")
            dg.type = sg.type
            changed.append(g)
    for g, e in keys:
        dest.groups[g].elements[e] = src.groups[g].elements[e].copy()
    for n in nodes:
        sn = src.nodes[n]
        dest.define_node(n, tuple(float(c) for c in v.P[n]), 0)
        dn = dest.nodes[n]
        dn.fix = list(sn.fix)
        dn.flags = set(sn.flags)
    if dest.group_active is None:
        dest.group_active = min(g for g, _ in keys)
    tmp = SSIModel()
    for g, e in keys:
        grp = tmp.groups.setdefault(g, Group(g, src.groups[g].type))
        grp.elements[e] = src.groups[g].elements[e]
    tmp.topl = []
    _copy_referenced_tables(tmp, src, etypes)
    added = 0
    differ = []
    for name in ("materials", "layers", "sections", "springs", "matrices"):
        dt, st = getattr(dest, name), getattr(tmp, name)
        for k, obj in st.items():
            if k not in dt:
                dt[k] = obj
                added += 1
            elif dt[k].to_json() != obj.to_json():
                differ.append(f"{name} {k}")
    if differ:
        notes.warn(f"table entries kept as defined in the destination (they differ from the source): "
                   f"{_short(differ)}")
    return TransferResult(len(keys), len(nodes), changed, added)


def element_centroids(m: SSIModel, gid: int, v: Optional[ModelView] = None) -> Dict[int, np.ndarray]:
    """Centroids (mean of the distinct DOF nodes) of the elements of group ``gid``."""
    v = v or model_view(m)
    g = m.groups[gid]
    out = {}
    for e in g.sorted_elements():
        ns = element_dof_nodes(g.type, e.nodes)
        if ns and v.defined(ns):
            out[e.id] = np.mean([v.P[n] for n in ns], axis=0)
    return out


@dataclass
class SplitResult:
    new_group: Optional[int]
    moved: Dict[int, int]            # old element number -> new element number in the new group
    kept: int


def split_group(m: SSIModel, gid: int, point: Sequence[float], normal: Sequence[float],
                notes: Optional[Notes] = None) -> SplitResult:
    """SPLITGROUP (spec 09 section 5.10, inferred algorithm): elements whose centroid lies on the +n side
    of the plane (signed distance > tol) move to a new group ``max group + 1`` of the same type, title
    suffixed ``_B``, renumbered 1..k in their original order; the others (and those on the plane) stay.
    EOUT requests and element stresses of the moved elements follow them."""
    notes = _notes(notes)
    if gid not in m.groups:
        raise SectionError(f"group {gid} is not defined")
    v = model_view(m)
    ez = unit(normal, "plane normal")
    P = np.asarray(point, float)
    grp = m.groups[gid]
    cents = element_centroids(m, gid, v)
    undefined = [e for e in grp.elements if e not in cents]
    if undefined:
        notes.warn(f"{len(undefined)} elements with undefined nodes stay in group {gid}: {_short(undefined)}")
    pos = [e for e in sorted(cents) if float((cents[e] - P) @ ez) > v.tol]
    kept = len(grp.elements) - len(pos)
    if not pos:
        notes.warn(f"no element of group {gid} lies on the positive side of the plane: nothing split")
        return SplitResult(None, {}, kept)
    if kept == 0:
        notes.warn(f"every element of group {gid} lies on the positive side of the plane")
    new = max(m.groups) + 1
    ng = Group(new, grp.type, (grp.title + "_B") if grp.title else f"group {gid}_B")
    moved: Dict[int, int] = {}
    for k, e in enumerate(pos, start=1):
        ng.elements[k] = grp.elements.pop(e).copy(new_id=k)
        moved[e] = k
    m.groups[new] = ng
    # EOUT requests of moved elements follow them
    reqs = []
    for r in m.eout:
        if r.group != gid:
            reqs.append(r)
            continue
        stay = [e for e in r.elements if e not in moved]
        go = [moved[e] for e in r.elements if e in moved]
        if stay:
            r.elements = stay
            reqs.append(r)
        if go:
            reqs.append(ElementRequest(list(r.codes), new, go))
    m.eout = reqs
    st = m.ui_state.get("element_stress")
    if isinstance(st, dict):
        for row in st.get("rows", []):
            if int(row[0]) == gid and int(row[1]) in moved:
                row[0], row[1] = new, moved[int(row[1])]
        if any(int(r[0]) == new for r in st.get("rows", [])):
            st.setdefault("groups", []).append([new, ELEMENT_TYPE_NAMES.get(grp.type, "")])
    if kept and pos:
        notes.info(f"group {gid} keeps {kept} elements with their numbers (gaps: ECOMPR compresses them)")
    return SplitResult(new, moved, kept)


# ======================================================================================
# Volumes and masses: CALCC, CALCM (spec 10 sections 3.3, 3.6)
# ======================================================================================
_HEX_NAT = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                     [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=float)


def _hex_gauss():
    g = 1.0 / math.sqrt(3.0)
    pts = np.array([[a, b, c] for a in (-g, g) for b in (-g, g) for c in (-g, g)])
    N = np.prod(1.0 + pts[:, None, :] * _HEX_NAT[None, :, :], axis=2) / 8.0        # (8gp, 8)
    dN = np.zeros((8, 3, 8))
    for k in range(3):
        f = np.ones((8, 8))
        for j in range(3):
            f = f * (_HEX_NAT[None, :, j] if j == k else (1.0 + pts[:, None, j] * _HEX_NAT[None, :, j]))
        dN[:, k, :] = f / 8.0
    return N, dN


_HEX_N, _HEX_DN = _hex_gauss()


def solid_volume_centroid(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Volume and volume centroid of trilinear 8-node SOLIDs ``X (nE, 8, 3)`` (repeated nodes for prisms
    and pyramids).  2 x 2 x 2 Gauss integration of ``det J`` and ``x det J`` is exact for the trilinear
    geometry (``det J`` is at most quadratic in each natural coordinate)."""
    X = np.asarray(X, float)
    J = np.einsum("gkn,enj->egkj", _HEX_DN, X)          # (nE, 8, 3, 3)
    det = np.linalg.det(J)                               # (nE, 8)
    V = det.sum(axis=1)
    xg = np.einsum("gn,enj->egj", _HEX_N, X)             # (nE, 8, 3)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.einsum("eg,egj->ej", det, xg) / V[:, None]
    return np.abs(V), c


def facet_area_centroid(X: np.ndarray) -> Tuple[float, np.ndarray]:
    """Area and centroid of a flat facet (3 / 4 corners; a warped quad is projected on its mean plane)."""
    from ..elements.shell import local_frame
    Lam, c, xy = local_frame(X)
    x, y = xy[:, 0], xy[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    a = x * y1 - x1 * y
    A = 0.5 * float(a.sum())
    if A == 0.0:
        return 0.0, c
    cx, cy = float(((x + x1) * a).sum()) / (6 * A), float(((y + y1) * a).sum()) / (6 * A)
    return abs(A), c + cx * Lam[0] + cy * Lam[1]


@dataclass
class ElementVolume:
    group: int
    element: int
    type: int
    volume: float
    centroid: np.ndarray
    excavated: bool
    mat: int


def element_volumes(m: SSIModel, notes: Optional[Notes] = None,
                    v: Optional[ModelView] = None) -> List[ElementVolume]:
    """Volumes of the elements (spec 10 section 3.3): SOLID exact (trilinear), SHELL / TSHELL area x
    thickness, PLANE area x 1 (unit plane-strain thickness), BEAMS length x R ``<axial>`` area;
    SPRING, GENERAL and lumped masses have no volume."""
    from ..elements.base import ElementError
    notes = _notes(notes)
    v = v or model_view(m)
    out: List[ElementVolume] = []
    bad: Dict[str, List[Key]] = {}
    solids: List[Tuple[Any, np.ndarray]] = []
    for r in v.elems:
        el = r.elem
        key = (r.group, r.id)
        if r.type == SOLID:
            ns = list(el.nodes) + [0] * (8 - len(el.nodes))
            if any(n == 0 for n in ns[:8]) or not v.defined(ns[:8]):
                bad.setdefault("undefined nodes", []).append(key)
                continue
            solids.append((r, np.array([v.P[n] for n in ns[:8]])))
        elif r.type in (SHELL, TSHELL, PLANE):
            corners = _shell_corners(el.nodes)
            if len(corners) < 3 or not v.defined(corners):
                bad.setdefault("undefined nodes", []).append(key)
                continue
            try:
                A, c = facet_area_centroid(np.array([v.P[n] for n in corners]))
            except ElementError:
                bad.setdefault("degenerate facet", []).append(key)
                continue
            t = 1.0 if r.type == PLANE else float(el.thick)
            if r.type != PLANE and not t > 0:
                bad.setdefault("no thickness (THICK)", []).append(key)
                continue
            out.append(ElementVolume(r.group, r.id, r.type, A * t, c, r.excavated, el.mat))
        elif r.type == BEAMS:
            ns = el.nodes
            if len(ns) < 2 or not ns[0] or not ns[1] or not v.defined(ns[:2]):
                bad.setdefault("undefined nodes", []).append(key)
                continue
            sec = m.sections.get(el.prop)
            if sec is None:
                bad.setdefault("BEAMS without an R property", []).append(key)
                continue
            L = float(np.linalg.norm(v.P[ns[1]] - v.P[ns[0]]))
            out.append(ElementVolume(r.group, r.id, BEAMS, L * sec.axial, 0.5 * (v.P[ns[0]] + v.P[ns[1]]), False,
                                     el.mat))
    if solids:
        V, C = solid_volume_centroid(np.array([x for _, x in solids]))
        for (r, _), vol, c in zip(solids, V, C):
            out.append(ElementVolume(r.group, r.id, SOLID, float(vol), c, r.excavated, r.elem.mat))
    out.sort(key=lambda ev: (ev.group, ev.element))
    for why, keys in sorted(bad.items()):
        notes.warn(f"{len(keys)} elements without volume ({why}): {_short(_ge(keys))}")
    return out


def centroid_volumes(m: SSIModel, notes: Optional[Notes] = None) -> List[ElementVolume]:
    """Element volumes for CALCC (spec 10 section 3.3).

    On a cross-section model (CSECT data, or a unit-thickness extrusion read back by WRITE -> INP,
    :func:`extrusion_normal`) the volume centroid must equal the area centroid of CALCPAR / CALCMOI.
    The extruded pieces have volume ``A_i x 1`` centred on the plane, so that holds for them; the
    BEAMS crossing the plane are drawn but are not section pieces (no ``.ess`` record, spec 10 section
    3.2), so they are left out here too (information)."""
    notes = _notes(notes)
    v = model_view(m)
    vols = element_volumes(m, notes, v)
    if is_csect_model(m) or extrusion_normal(m, v) is not None:
        beams = [(ev.group, ev.element) for ev in vols if ev.type == BEAMS]
        if beams:
            notes.info(f"cross-section model: {len(beams)} BEAMS elements are drawn but are not section pieces, left "
                       f"out so that the volume centroid is the area centroid of CALCPAR / CALCMOI: "
                       f"{_short(_ge(beams))}")
            vols = [ev for ev in vols if ev.type != BEAMS]
    return vols


def volume_centroid(vols: Sequence[ElementVolume]) -> Tuple[float, np.ndarray]:
    """CALCC: total volume and ``Xc = sum V_e c_e / sum V_e``."""
    Vt = float(sum(ev.volume for ev in vols))
    if not Vt > 0:
        raise SectionError("the model has no element with a volume (SOLID, SHELL, PLANE, BEAMS)")
    c = sum(ev.volume * ev.centroid for ev in vols) / Vt
    return Vt, np.asarray(c, float)


@dataclass
class MassSummary:
    weight: float                    # structural elements
    mass: float
    soil_weight: float               # excavated soil elements (reported separately)
    soil_mass: float
    lumped: np.ndarray               # (6,) MT / MR in mass units
    gravity: float
    counts: Dict[str, int] = field(default_factory=dict)

    @property
    def total(self) -> np.ndarray:
        return self.mass + self.lumped[:3]


def model_masses(m: SSIModel, notes: Optional[Notes] = None) -> MassSummary:
    """CALCM (spec 10 section 3.6): element mass ``(w/g) V_e`` with ``w`` the M ``<weight>`` (structure) or
    L ``<weight>`` (excavated soil, reported separately), and the MT / MR lumped masses converted by
    MUNITS (weights divided by g)."""
    notes = _notes(notes)
    g = float(m.gravity)
    if not g > 0:
        raise SectionError(f"gravity {g:g} must be > 0 (GRAVITY)")
    v = model_view(m)
    vols = element_volumes(m, notes, v)
    W = Ws = 0.0
    nomat: List[Key] = []
    for ev in vols:
        if ev.excavated:
            lay = m.layers.get(ev.mat)
            if lay is None:
                nomat.append((ev.group, ev.element))
                continue
            Ws += lay.weight * ev.volume
        else:
            mat = m.materials.get(ev.mat)
            if mat is None:
                nomat.append((ev.group, ev.element))
                continue
            W += mat.weight * ev.volume
    if nomat:
        notes.warn(f"{len(nomat)} elements without a material / soil layer: {_short(_ge(nomat))}")
    ngen = sum(1 for gr in m.groups.values() if gr.type == GENERAL for _ in gr.elements)
    if ngen:
        notes.info(f"{ngen} GENERAL elements: their MXM mass matrices are not included")
    lumped = np.zeros(6)
    for _, mv in m.nodal_masses(g).items():
        lumped += mv
    return MassSummary(W, W / g, Ws, Ws / g, lumped, g, {"elements": len(vols)})


# ======================================================================================
# CALCSECTHIST (spec 10 section 3.4)
# ======================================================================================
def read_frame_list(path: PathLike, cwd: Optional[PathLike] = None) -> List[Path]:
    """Stress-history list file: one header line (ignored), then one ``.ess`` path per line; relative
    paths resolve against the list file's folder, then ``cwd``.  A first line that names an existing
    file is taken as a frame (lists written without a header); ``#`` comment lines are skipped."""
    path = Path(path)
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise SectionError(f"{path}: cannot read ({exc.strerror or exc})") from None
    entries = [ln.strip() for ln in lines]
    out: List[Path] = []

    def resolve(s: str) -> Optional[Path]:
        cands = [s]
        tok = s.split()[0] if s.split() else s
        if tok != s:
            cands.append(tok)
        for c in cands:
            p = Path(c)
            if p.is_absolute():
                if p.exists():
                    return p
                continue
            for base in (path.parent, Path(cwd) if cwd is not None else None):
                if base is not None and (base / p).exists():
                    return base / p
        return None

    first = True
    for s in entries:
        if not s or s[0] in "#*!":
            if first and s:
                first = False
            continue
        if first:
            first = False
            r = resolve(s)
            if r is None or not r.is_file():
                continue                          # the header line
            out.append(r)
            continue
        r = resolve(s)
        if r is None:
            raise SectionError(f"{path.name}: frame file {s!r} not found (relative to {path.parent} or the CWD)")
        out.append(r)
    if not out:
        raise SectionError(f"{path.name}: the list holds no frame file")
    return out


def history_rows(T: np.ndarray, pieces: Sequence[Piece], frames: Sequence[Path],
                 notes: Optional[Notes] = None) -> np.ndarray:
    """Resultant histories ``(nframes, 6)`` from ``.ess`` frames; pieces without a record in a frame
    contribute nothing (one summary warning)."""
    notes = _notes(notes)
    out = np.zeros((len(frames), 6))
    missing: Dict[Key, int] = {}
    kinds = {p.group: p.kind for p in pieces}
    matched = False
    mismatch = None
    for k, f in enumerate(frames):
        ess = read_ess(f)
        bad = [f"group {g} is {ess.groups[g]} in the frame but {t} in the model" for g, t in sorted(kinds.items())
               if g in ess.groups and ess.groups[g] != t]
        if bad:
            raise SectionError(f"{Path(f).name} does not belong to the active model: " + "; ".join(bad))
        S, ok = gather_stresses(pieces, ess.values)
        for i in np.nonzero(~ok)[0]:
            missing[pieces[i].key] = missing.get(pieces[i].key, 0) + 1
        matched = matched or bool(ok.any())
        if not ok.any() and mismatch is None:
            mismatch = stress_mismatch(ess.values, [p.key for p in pieces], f"are in frame {Path(f).name}",
                                       "the frames were written before a renumbering of the model (GCOM, "
                                       "MERGEGROUP, ECOMPR); use the model they were computed for")
        out[k] = resultants(T, S, ok)
    if not matched and mismatch:
        notes.warn(mismatch)
    if missing:
        notes.warn(f"{len(missing)} section elements have no record in some frames (they contribute nothing "
                   f"there): " + _short([f"{g}/{e} ({n} frames)" for (g, e), n in sorted(missing.items())]))
    return out


def signed_absmax(R: np.ndarray) -> np.ndarray:
    """Per column, the value of largest magnitude with its sign (ties: first occurrence), D-SEC-02."""
    R = np.asarray(R, float)
    idx = np.argmax(np.abs(R), axis=0)
    return R[idx, np.arange(R.shape[1])]


def write_section_csv(path: PathLike, R: np.ndarray, ts: float) -> Path:
    """7-column CSV ``Time,Fx,Fy,Fz,Mx,My,Mz`` (``Step,...`` when ts <= 0) with a final ``MAX`` row of the
    signed absolute maxima (D-SEC-02, D-SEC-05: time ``(k-1) ts``, step ``k``)."""
    path = Path(path)
    lines = [("Time" if ts > 0 else "Step") + "," + ",".join(RESULTANT_NAMES)]
    for k, row in enumerate(R, start=1):
        t = f"{(k - 1) * ts:.6E}" if ts > 0 else str(k)
        lines.append(t + "," + ",".join(f"{x:.6E}" for x in row))
    lines.append("MAX," + ",".join(f"{x:.6E}" for x in signed_absmax(R)))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
