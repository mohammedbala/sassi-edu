"""Numerical kernels and text formats of the STRESS module (requirements section 4.10).

What STRESS computes (for the structural engineer)
--------------------------------------------------
HOUSE stored, for every output-capable element, a *recovery operator* ``S`` (FILE4 arrays
``rec_<T>_S``): a complex matrix that turns the element's nodal displacement vector ``u_e`` into
its output components (centroid stresses of a SOLID, end forces of a BEAMS element, ...)::

    r_c = S_c . u_e                               (requirements 4.10 item 1, ARCHITECTURE 6.2)

ANALYS stored the nodal transfer functions ``H(f)`` (FILE8) at the SSI frequencies.  Because the
recovery is linear, the *stress transfer function* (STF) of a component at an SSI frequency is
simply ``STF_c(f_j) = S_c . U_e(f_j)`` with ``U_e`` the element's columns of FILE8.  Rigid-body
motion produces no strain, so the *total* motion TFs of a seismic FILE8 can be used directly
(spec 05d section 1.2 item 3): ``S`` annihilates rigid translations (and rigid rotations for the
elements that carry rotations consistently: SOLID, PLANE, SHELL, BEAMS).

This module holds the pieces that do not depend on the module's file handling:

* :func:`element_dof_tf` / :func:`element_stf` -- gather the element DOF TFs from FILE8 columns and
  apply the recovery operators (vectorised over elements and frequencies);
* :func:`octahedral_shear_stress`, :func:`octahedral_shear_strain`, :func:`max_shear_strain_2d` --
  the derived quantities that are non-linear in the components and therefore evaluated in the
  time domain (requirements 4.10 items 3 and 5);
* the ``ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`` / ``ESTRESS_n.ess`` / ``STATIC_SOIL_PRESSURES.TXT``
  layout (D-FIL-06, D-STR-12) with its reader, used later by READSTR and the section cuts;
* ``Frames.txt`` (D-FIL-07, D-STR-08), nodal averaging of element-centre values (plotting only,
  requirements 4.10 item 4) and the soil-pressure geometry of D-STR-09.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
import scipy.sparse as sp

__all__ = [
    "element_dof_tf", "element_stf", "rigid_body_dofs", "octahedral_shear_stress", "octahedral_shear_strain",
    "max_shear_strain_2d", "OCT_STRESS_PURE_SHEAR", "OCT_STRESS_UNIAXIAL",
    "CENTER_COLUMNS", "CENTER_TYPES", "center_columns", "ordered_group_numbers",
    "CenterBlock", "write_element_center", "read_element_center",
    "read_frames_file", "output_steps", "averaging_matrix", "average_histories", "membrane_to_global",
    "HEX_FACES", "solid_faces", "shared_faces", "face_normal", "normal_stress",
    "strain_component_names",
]

PathLike = Union[str, Path]

#: tau_oct / tau for pure shear and tau_oct / sigma for uniaxial stress (VP-35, spec 05d section 8)
OCT_STRESS_PURE_SHEAR = float(np.sqrt(6.0) / 3.0)      # 0.816497
OCT_STRESS_UNIAXIAL = float(np.sqrt(2.0) / 3.0)        # 0.471405


# ======================================================================================
# Stress transfer functions (requirements 4.10 item 1)
# ======================================================================================
def element_dof_tf(H: np.ndarray, eq: np.ndarray) -> np.ndarray:
    """Element DOF transfer functions ``U (nF, nE, nd)`` from FILE8 ``H (nF, nEq)``.

    ``eq (nE, nd)`` are the FILE4 equation numbers of the element DOFs (``rec_<T>_eq``); ``-1``
    marks a fixed DOF (zero total motion) or a padding column (SHELL triangles) and contributes 0.
    """
    H = np.asarray(H)
    eq = np.asarray(eq, dtype=np.int64)
    U = H[:, np.clip(eq, 0, None)]                            # (nF, nE, nd)
    if np.any(eq < 0):
        U = np.where((eq < 0)[None, :, :], 0.0, U)
    return U


def element_stf(S: np.ndarray, eq: np.ndarray, H: np.ndarray, rigid: Optional[np.ndarray] = None) -> np.ndarray:
    """Stress transfer functions ``STF (nF, nE, nc) = S_e . U_e(f)`` (requirements 4.10 item 1).

    ``S (nE, nc, nd)`` recovery operators (FILE4 ``rec_<T>_S``, or the strain operator
    ``x_rec_<T>_B``), ``eq (nE, nd)`` their equation numbers and ``H (nF, nEq)`` the FILE8 TFs.

    ``rigid`` (``(nd,)`` or ``(nE, nd)``, optional): a rigid-body translation of the element DOFs
    (seismic: the control direction projected on each DOF, :func:`rigid_body_dofs`) that is
    subtracted from *every* element DOF, fixed ones included, before the recovery:
    ``STF = S_e . (U_e - r_e)``.  Every recovery operator annihilates a uniform translation, so the
    result is ``S_e . U_e`` in exact arithmetic; numerically it avoids the cancellation of the
    large rigid-body terms of total-motion TFs (``H ~ 1``) at low frequency, where the deformation
    is ``O((f/f_n)^2)`` of the motion and the plain product loses that many digits.
    """
    U = element_dof_tf(H, eq)
    if rigid is not None:
        r = np.asarray(rigid)
        U = U - (r[None, None, :] if r.ndim == 1 else r[None, :, :])
    # numpy 2.0 + macOS Accelerate raises spurious floating-point flags in complex matmul
    # (sassi.elements.base.blas_quiet); the inputs are checked finite by the caller
    with np.errstate(divide="ignore", over="ignore", invalid="ignore", under="ignore"):
        return np.einsum("ecd,fed->fec", np.asarray(S), U, optimize=True)


def rigid_body_dofs(dofs: Sequence[int], nnodes: int, cm: int = 0, ang_deg: float = 0.0) -> np.ndarray:
    """Element-DOF vector ``(nnodes * len(dofs),)`` of the seismic rigid-body motion: the unit control
    motion (SITE direction ``cm`` rotated by the ANALYS angle) projected on each translational DOF,
    0 on rotations (the zero-frequency anchor of D-MOT-03, :func:`sassi.core.interp.rigid_body_anchor`).
    Element DOFs are node-major in ``dofs`` order (ARCHITECTURE 6.2)."""
    from .interp import rigid_body_anchor
    per_node = [rigid_body_anchor(int(k), cm, ang_deg) for k in dofs]
    return np.tile(np.asarray(per_node, dtype=float), int(nnodes))


# ======================================================================================
# Derived quantities, evaluated in the time domain (requirements 4.10 items 3 and 5)
# ======================================================================================
def octahedral_shear_stress(sxx, syy, szz, sxy, sxz, syz):
    """``tau_oct = 1/3 sqrt[(sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 + 6 (sxy^2 + syz^2 + sxz^2)]``.

    The shear stress on the octahedral planes (normals equally inclined to the principal axes);
    ``sqrt(3/2)`` times it is the von Mises stress.  Pure shear tau gives ``(sqrt 6/3) tau``,
    uniaxial sigma gives ``(sqrt 2/3) sigma`` (spec 05d section 8 item 1).  Being non-linear in the
    components it is computed from the component *histories* (requirements 4.10 item 3).
    """
    sxx, syy, szz = np.asarray(sxx), np.asarray(syy), np.asarray(szz)
    sxy, sxz, syz = np.asarray(sxy), np.asarray(sxz), np.asarray(syz)
    return np.sqrt((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2
                   + 6.0 * (sxy ** 2 + syz ** 2 + sxz ** 2)) / 3.0


def octahedral_shear_strain(exx, eyy, ezz, gxy, gxz, gyz):
    """Engineering octahedral shear strain (requirements 4.10 item 5, ISTR 1, 3D)::

        gamma_oct = 2/3 sqrt[(exx-eyy)^2 + (eyy-ezz)^2 + (ezz-exx)^2 + 1.5 (gxy^2 + gyz^2 + gxz^2)]

    with engineering shear strains ``g = 2 eps``.  Simple shear gamma gives ``(sqrt 6/3) gamma``.
    """
    exx, eyy, ezz = np.asarray(exx), np.asarray(eyy), np.asarray(ezz)
    gxy, gxz, gyz = np.asarray(gxy), np.asarray(gxz), np.asarray(gyz)
    return (2.0 / 3.0) * np.sqrt((exx - eyy) ** 2 + (eyy - ezz) ** 2 + (ezz - exx) ** 2
                                 + 1.5 * (gxy ** 2 + gyz ** 2 + gxz ** 2))


def max_shear_strain_2d(exx, ezz, gxz):
    """Maximum (engineering) shear strain of a plane-strain X-Z state, ``sqrt[(exx-ezz)^2 + gxz^2]``
    (requirements 4.10 item 5, 2D)."""
    exx, ezz, gxz = np.asarray(exx), np.asarray(ezz), np.asarray(gxz)
    return np.sqrt((exx - ezz) ** 2 + gxz ** 2)


def strain_component_names(etype: str, components: Sequence[str]) -> List[str]:
    """Strain output names of D-STR-04 (``EDUOPT,STRAINOUT,1``): the stress code with its first
    letter replaced by ``E`` (SOLID ``SXX -> EXX`` ... ``SYZ -> EYZ``, ``SOCT -> EOCT``; PLANE
    ``SXX -> EXX``, ``SZZ -> EZZ``, ``TXZ -> EXZ``).  Shear strains are engineering strains."""
    return ["E" + c[1:] for c in components]


# ======================================================================================
# ELEMENT_CENTER layout (D-FIL-06, D-STR-12, spec 02 section 3.9)
# ======================================================================================
#: element types written to ELEMENT_CENTER_ABS_MAX_STRESSES.TXT / .ess (BEAMS need end forces and
#: are excluded; TSHELL is P2) and the component of each of the six columns (None = column of zeros)
CENTER_COLUMNS: Dict[str, Tuple[Optional[str], ...]] = {
    "SOLID": ("SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ"),
    "SHELL": ("FXX", "FYY", "FXY", "MXX", "MYY", "MXY"),
    "SPRING": ("FX", "FY", "FZ", "MXX", "MYY", "MZZ"),
    "PLANE": ("SXX", "SZZ", "TXZ", None, None, None),
}
CENTER_TYPES = tuple(CENTER_COLUMNS)


def center_columns(etype: str, components: Sequence[str]) -> List[int]:
    """Indices into ``components`` of the six ELEMENT_CENTER columns (-1 = zero column)."""
    pos = {c: i for i, c in enumerate(components)}
    return [pos[c] if c is not None else -1 for c in CENTER_COLUMNS[etype]]


def ordered_group_numbers(groups: Iterable[Tuple[int, str]]) -> Dict[int, int]:
    """'ordered group #' of D-FIL-06: the 1-based order of a group among the groups of its element
    type, groups taken in ascending group number."""
    seen: Dict[str, int] = {}
    out: Dict[int, int] = {}
    for g, t in sorted(set((int(g), str(t)) for g, t in groups)):
        seen[t] = seen.get(t, 0) + 1
        out[g] = seen[t]
    return out


@dataclass
class CenterBlock:
    """One group block of the ELEMENT_CENTER layout."""
    etype: str
    group: int
    ordered: int
    elements: np.ndarray                      # (n,) element numbers within the group
    values: np.ndarray                        # (n, ncol) values (6 columns, 1 for soil pressures)

    def __post_init__(self):
        self.elements = np.asarray(self.elements, dtype=np.int64).reshape(-1)
        self.values = np.asarray(self.values, dtype=float).reshape(len(self.elements), -1)


def write_element_center(path: PathLike, blocks: Sequence[CenterBlock], fmt: str = "{:.9e}") -> Path:
    """Write the ELEMENT_CENTER layout (spec 02 section 3.9, D-FIL-06)::

        [# of groups]
        [elem type] [group #] [ordered group #] [# elements in group]     (one line per group)
        [elem type] [group #] [ordered group #]                           (block header per group)
        [element #] [comp1] ... [comp6]                                   (one line per element)

    Used for ``ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`` (absolute maxima), ``ESTRESS_n.ess`` (signed
    values at one time step, D-STR-12), the ``.sig/.tau/.bdsig/.bdtau`` maxima and, with one value
    column, ``STATIC_SOIL_PRESSURES.TXT`` / ``pres_max_ele``.
    """
    path = Path(path)
    lines = [str(len(blocks))]
    for b in blocks:
        lines.append(f"{b.etype:<8s}{b.group:6d}{b.ordered:6d}{len(b.elements):8d}")
    for b in blocks:
        lines.append(f"{b.etype:<8s}{b.group:6d}{b.ordered:6d}")
        for e, row in zip(b.elements, b.values):
            lines.append(f"{int(e):10d} " + " ".join(fmt.format(float(v)) for v in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def read_element_center(path: PathLike) -> List[CenterBlock]:
    """Read the ELEMENT_CENTER layout (any number of value columns); blank lines are skipped."""
    rows = [ln.split() for ln in Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")]
    if not rows:
        raise ValueError(f"{Path(path).name}: empty file")
    try:
        ng = int(rows[0][0])
        heads = [(r[0].upper(), int(r[1]), int(r[2]), int(r[3])) for r in rows[1:1 + ng]]
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{Path(path).name}: malformed group table ({exc})") from None
    blocks: List[CenterBlock] = []
    i = 1 + ng
    for etype, g, og, ne in heads:
        if i >= len(rows) or rows[i][0].upper() != etype or int(rows[i][1]) != g:
            raise ValueError(f"{Path(path).name}: block header of group {g} ({etype}) not found")
        i += 1
        el, vals = [], []
        for _ in range(ne):
            if i >= len(rows):
                raise ValueError(f"{Path(path).name}: group {g} has fewer than {ne} element rows")
            el.append(int(rows[i][0]))
            vals.append([float(v.replace("D", "E")) for v in rows[i][1:]])
            i += 1
        ncol = min((len(v) for v in vals), default=0)
        blocks.append(CenterBlock(etype, g, og, np.array(el), np.array([v[:ncol] for v in vals]).reshape(ne, ncol)))
    return blocks


# ======================================================================================
# Frames.txt and output steps (D-FIL-07, D-STR-07, D-STR-08)
# ======================================================================================
def read_frames_file(path: PathLike) -> Tuple[List[int], List[int]]:
    """Parse ``Frames.txt``: ``[# frames]``, the frame numbers, ``[# soil-pressure groups]`` and the
    groups (on one line in the manual's example; any whitespace is accepted).

    Frame numbers are 1-based time-step indices (D-STR-08: frame ``00001`` is t = 0).
    """
    toks: List[int] = []
    for ln in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.split("!")[0].strip()
        if not s or s[0] in "#*":
            continue
        for t in s.replace(",", " ").split():
            try:
                toks.append(int(float(t)))
            except ValueError:
                raise ValueError(f"{Path(path).name}: non-numeric entry {t!r}") from None
    if not toks:
        return [], []
    n = toks[0]
    if n < 0 or len(toks) < 1 + n:
        raise ValueError(f"{Path(path).name}: {n} frames announced, {len(toks) - 1} values found")
    frames = toks[1:1 + n]
    rest = toks[1 + n:]
    groups: List[int] = []
    if rest:
        m = rest[0]
        if m < 0 or len(rest) < 1 + m:
            raise ValueError(f"{Path(path).name}: {m} soil-pressure groups announced, {len(rest) - 1} found")
        groups = rest[1:1 + m]
    return frames, groups


def output_steps(nout: int, skip: int) -> np.ndarray:
    """0-based indices of the output samples: every ``skip``-th sample, 0/1 = all (D-STR-07)."""
    return np.arange(0, int(nout), max(1, int(skip)))


# ======================================================================================
# Nodal averages (plotting only, requirements 4.10 item 4) and soil pressures (D-STR-09)
# ======================================================================================
def averaging_matrix(elem_nodes: Sequence[Sequence[int]]) -> Tuple[np.ndarray, sp.csr_matrix]:
    """Plain-mean nodal averaging of element-centre values.

    ``elem_nodes``: node ids of each element (0 = unused slot; repeated ids of degenerate elements
    count once).  Returns ``(nodes, A)`` with ``nodes`` sorted ascending and ``A (nN, nE)`` such
    that ``A @ v_elem`` is the mean over the elements attached to each node.  No shape-function
    extrapolation is used: the values are for plotting only (spec 05d section 1.9).
    """
    rows, cols = [], []
    for e, ns in enumerate(elem_nodes):
        for n in sorted(set(int(v) for v in ns if int(v) > 0)):
            rows.append(n)
            cols.append(e)
    nodes = np.array(sorted(set(rows)), dtype=np.int64)
    pos = {int(n): i for i, n in enumerate(nodes)}
    r = np.array([pos[n] for n in rows], dtype=np.int64)
    c = np.array(cols, dtype=np.int64)
    A = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(nodes), len(elem_nodes)))
    cnt = np.asarray(A.sum(axis=1)).ravel()
    A = sp.diags(1.0 / np.where(cnt > 0, cnt, 1.0)) @ A
    return nodes, sp.csr_matrix(A)


def average_histories(A: sp.spmatrix, E: np.ndarray) -> np.ndarray:
    """Nodal averages ``N[t] = A @ E[t]`` of element-centre histories ``E (nt, nE, nc)`` -> ``(nt, nN, nc)``.

    One sparse product for all time steps and columns: ``A`` has only a few non-zeros per column
    (8 per SOLID, 4 per SHELL), so the cost is O(nnz(A) nt nc) instead of the O(nN nE nt nc) of a dense
    product (D-GEN-10 performance target; review finding on the nodal stress frames)."""
    E = np.asarray(E)
    nt, nE, nc = E.shape
    if A.shape[1] != nE:
        raise ValueError(f"averaging matrix has {A.shape[1]} columns for {nE} elements")
    N = sp.csr_matrix(A) @ E.transpose(1, 0, 2).reshape(nE, nt * nc)
    return np.asarray(N).reshape(A.shape[0], nt, nc).transpose(1, 0, 2)


def membrane_to_global(F: np.ndarray, Lam: np.ndarray) -> np.ndarray:
    """SHELL membrane stresses in the element local axes -> 3D stress tensor components in global axes.

    ``F (..., nE, 3)``: local ``FXX FYY FXY`` (stresses, D-STR-10 / D-W1-14); ``Lam (nE, 3, 3)``: rows
    ``x' y' z'`` of each element (spec 08 4.6, :func:`sassi.elements.shell.local_frame`).  A membrane
    state is plane stress in the shell plane, ``sigma_l = [[FXX, FXY, 0], [FXY, FYY, 0], [0, 0, 0]]``,
    and the tensor transforms as ``sigma_g = Lam^T sigma_l Lam``, i.e. with ``e1 = x'``, ``e2 = y'``::

        sigma_g = FXX e1 e1^T + FYY e2 e2^T + FXY (e1 e2^T + e2 e1^T)

    Returns ``(..., nE, 6)`` in the SOLID order ``SXX SYY SZZ SXY SXZ SYZ``, so that SHELL membrane stresses
    can be averaged with SOLID stresses at shared nodes (spec 05d 1.9: a model with SOLID and SHELL elements
    gets membrane-stress frames only).  The transformation is linear, so it applies to histories and
    transfer functions alike."""
    F = np.asarray(F)
    Lam = np.asarray(Lam, dtype=float)
    e1, e2 = Lam[:, 0, :], Lam[:, 1, :]
    pairs = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
    T = np.empty((Lam.shape[0], 6, 3))
    for k, (i, j) in enumerate(pairs):
        T[:, k, 0] = e1[:, i] * e1[:, j]
        T[:, k, 1] = e2[:, i] * e2[:, j]
        T[:, k, 2] = e1[:, i] * e2[:, j] + e2[:, i] * e1[:, j]
    return np.einsum("ekc,...ec->...ek", T, F)


#: the six quadrilateral faces of the 8-node SOLID (node slots 1-4 bottom, 5-8 top; spec 08 4.4)
HEX_FACES = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))


def solid_faces(nodes8: Sequence[int]) -> List[Tuple[int, ...]]:
    """Faces of a SOLID element as tuples of distinct node ids (>= 3); degenerate faces of prisms
    and pyramids (fewer than 3 distinct nodes) are dropped."""
    ns = [int(v) for v in nodes8]
    out = []
    for f in HEX_FACES:
        ids = []
        for k in f:
            if ns[k] > 0 and ns[k] not in ids:
                ids.append(ns[k])
        if len(ids) >= 3:
            out.append(tuple(ids))
    return out


def shared_faces(nodes8: Sequence[int], structure_nodes: Set[int]) -> List[Tuple[int, ...]]:
    """Faces of a soil SOLID whose nodes all belong to the structure (D-STR-09)."""
    return [f for f in solid_faces(nodes8) if all(n in structure_nodes for n in f)]


def face_normal(xyz: np.ndarray) -> np.ndarray:
    """Unit normal of a planar or slightly warped polygon (Newell's method)."""
    p = np.asarray(xyz, dtype=float)
    q = np.roll(p, -1, axis=0)
    n = np.array([np.sum((p[:, 1] - q[:, 1]) * (p[:, 2] + q[:, 2])),
                  np.sum((p[:, 2] - q[:, 2]) * (p[:, 0] + q[:, 0])),
                  np.sum((p[:, 0] - q[:, 0]) * (p[:, 1] + q[:, 1]))])
    norm = np.linalg.norm(n)
    if norm == 0.0:
        raise ValueError("degenerate face (zero area)")
    return n / norm


def normal_stress(sig6: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Normal stress ``n^T sigma n`` for stress vectors ``(..., 6)`` in the SOLID order
    ``SXX SYY SZZ SXY SXZ SYZ`` (the sign of ``n`` does not matter)."""
    s = np.asarray(sig6)
    nx, ny, nz = (float(v) for v in n)
    return (s[..., 0] * nx * nx + s[..., 1] * ny * ny + s[..., 2] * nz * nz
            + 2.0 * (s[..., 3] * nx * ny + s[..., 4] * nx * nz + s[..., 5] * ny * nz))
