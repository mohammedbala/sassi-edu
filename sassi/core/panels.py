"""Wall panels of Option NON: geometry, kinematics, shear capacities (SHEAR) and backbone curves (BBCGEN).

What a panel is (manual section 1.5.4, spec 05d section 3.4)
-----------------------------------------------------------
A *panel* is a group of coplanar SHELL elements in a vertical plane, assumed to deform in uniform
shear (or bending).  Only its **four corner nodes** are used: bottom-left BL, bottom-right BR,
top-right TR and top-left TL in the panel's local axes ``e_h`` (horizontal, in the panel plane) and
``e_v`` = global Z.  With the in-plane displacements ``u = d . e_h`` and ``w = d . e_v`` of the corners,
the width L and the height H, the deformations used by NONLINEAR are (requirements 4.15; inferred
rigid-body-invariant formulas of spec 05d section 3.4, OQ-N2)::

    gamma = 1/2 [(u_TL - u_BL) + (u_TR - u_BR)]/H + 1/2 [(w_BR - w_BL) + (w_TR - w_TL)]/L   shear strain
    kappa = (theta_t - theta_b)/H,  theta_b = (w_BR - w_BL)/L,  theta_t = (w_TR - w_TL)/L    curvature
    eps_v = 1/2 [(w_TL - w_BL) + (w_TR - w_BR)]/H                                            axial strain

A rigid translation or an in-plane rigid rotation of the panel gives zero for all three.

Shear capacities (SHEAR, BBCGEN; spec 10 section 3.11, decisions D-NON-09, D-NON-10)
------------------------------------------------------------------------------------
The empirical equations (as given by Gulec & Whittaker 2009) use sqrt(f'c) in **psi**, areas in in^2
and forces in lb.  The panel geometry comes from the model (feet or metres, detected from the gravity
constant: British if g > 20); inputs are ksi / kips (British) or kN/m^2 / kN (SI) and are converted
internally with the constants of spec 10 section 3.11; results are kips or kN::

    ACI 318-08:            V = (alpha_c sqrt(f'c) + rho_H f_y) A_W  <=  10 sqrt(f'c) A_W
    Wood 1990:             6 sqrt(f'c) A_W  <=  V = rho_V A_W f_y / 4  <=  10 sqrt(f'c) A_W
    Barda et al. 1977:     V = (8 sqrt(f'c) - 2.5 sqrt(f'c) h_W/l_W + N_U/(4 l_W t_W) + rho_V f_y) t_W (0.6 l_W)
    Gulec-Whittaker 2009:  V = (1.5 sqrt(f'c) A_W + 0.25 F_VW + 0.20 F_BE + 0.40 N_U) / sqrt(h_W/l_W)

with ``alpha_c`` = 3.0 (h_W/l_W <= 1.5), 2.0 (>= 2.0), linear in between (ACI 318-08 21.9.4.1),
``F_VW = rho_V A_W f_y`` and ``F_BE = A_BE f_y,BE`` (D-NON-10).

BBCGEN backbone (D-NON-09): 22 points -- point 1 the cracking point ``(gamma_cr, V_cr)`` with
``V_cr = 3 sqrt(f'c) A_W`` (ASCE 4-17 C.3.3.2) or ``CFL V_u``, ``gamma_cr = V_cr/(G A_W)``; points 2-21
equally spaced in strain up to the yield point ``(gamma_y, V_u)`` (gamma_y = 0.004) on
``V = V_cr + (V_u - V_cr)[1 - (1 - xi)^2]``, ``xi = (gamma - gamma_cr)/(gamma_y - gamma_cr)``; point 22
the failure point ``(0.02, 1.02 V_u)``; yield index 21.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

import numpy as np

#: corner labels in the order they are stored (deck, NONLINMOTDISP)
CORNERS = ("BL", "BR", "TR", "TL")
#: |n_z| above which a shell plane is not vertical (sin of 1 degree)
VERTICAL_TOL = math.sin(math.radians(1.0))
#: coplanarity tolerance relative to the panel size
COPLANAR_TOL = 1.0e-3
#: SHELL and TSHELL group types
PANEL_TYPES = (3, 5)

# unit conversions of spec 10 section 3.11 (to in, in^2, psi, lb) and back
US = dict(L=12.0, A=144.0, f=1000.0, N=1000.0, out=1.0e-3, force="kips", stress="ksi", length="ft")
SI = dict(L=39.3700787, A=1550.0031, f=0.1450377, N=224.808943, out=0.00444822162, force="kN",
          stress="kN/m^2", length="m")
SHEAR_MODELS = {1: "ACI 318-08", 2: "Wood 1990", 3: "Barda 1977", 4: "Gulec-Whittaker 2009"}
#: BBCGEN defaults (D-NON-09)
GAMMA_Y = 0.004
GAMMA_FAIL = 0.02
FAIL_FACTOR = 1.02
BBCGEN_STEPS = 20


class PanelError(ValueError):
    """Invalid panel geometry or capacity input."""


# ======================================================================================
# Units
# ======================================================================================
def unit_system(gravity: float) -> Tuple[str, Dict[str, float], str]:
    """``('british' | 'si', conversion table, warning text)`` from the gravity constant (spec 10 section 1.5)."""
    g = float(gravity)
    if g > 20.0:
        warn = "" if abs(g / 32.174 - 1.0) <= 0.05 else f"gravity {g:g} is not within 5 % of 32.174 ft/s2"
        return "british", US, warn
    warn = "" if abs(g / 9.80665 - 1.0) <= 0.05 else f"gravity {g:g} is not within 5 % of 9.80665 m/s2"
    return "si", SI, warn


# ======================================================================================
# Shear capacities (SHEAR) and BBCGEN
# ======================================================================================
def alpha_c(ratio: float) -> float:
    """ACI 318-08 21.9.4.1: 3.0 for h_W/l_W <= 1.5, 2.0 for >= 2.0, linear in between."""
    if ratio <= 1.5:
        return 3.0
    if ratio >= 2.0:
        return 2.0
    return 3.0 - 2.0 * (ratio - 1.5)


@dataclass
class ShearCapacities:
    """Capacities of one panel in the model's force unit (kips or kN)."""
    aci_raw: float
    aci: float               # min(V_ACI, upper bound)
    upper: float             # 10 sqrt(f'c) A_W
    wood_raw: float
    wood_lower: float        # 6 sqrt(f'c) A_W
    wood: float              # clamped to [lower, upper]
    barda: float
    gw: float
    cracking: float          # 3 sqrt(f'c) A_W (ASCE 4-17 C.3.3.2)
    alpha_c: float
    units: str
    force_unit: str

    def ultimate(self, model: int) -> float:
        """V_u of BBCGEN ShearModel 1 ACI (capped), 2 Wood (clamped), 3 Barda, 4 Gulec-Whittaker."""
        m = int(model)
        if m == 1:
            return self.aci
        if m == 2:
            return self.wood
        if m == 3:
            return self.barda
        if m == 4:
            return self.gw
        raise PanelError(f"ShearModel {model} must be 1 ACI 318-08, 2 Wood 1990, 3 Barda 1977 or 4 Gulec-Whittaker 2009")


def shear_capacities(h_w: float, l_w: float, t_w: float, fc: float, fy: float, rho: float, n_u: float = 0.0,
                     a_be: float = 0.0, fy_be: float = 0.0, gravity: float = 32.2, force_args: bool = False) -> ShearCapacities:
    """The four shear-strength equations for a panel h_w x l_w x t_w (model length unit).

    ``fc``, ``fy``, ``fy_be`` in ksi (British) or kN/m^2 (SI); ``n_u`` in kips or kN; ``a_be`` in ft^2 or
    m^2.  With ``force_args`` (EDUOPT,SHEARFORCEARGS,1) ``a_be`` and ``fy_be`` are the forces F_VW and F_BE.
    """
    if min(h_w, l_w, t_w) <= 0:
        raise PanelError(f"panel dimensions must be > 0 (h_W {h_w:g}, l_W {l_w:g}, t_W {t_w:g})")
    if fc <= 0:
        raise PanelError("f'c must be > 0")
    units, u, _ = unit_system(gravity)
    hw, lw, tw = h_w * u["L"], l_w * u["L"], t_w * u["L"]
    aw = lw * tw
    fc_psi, fy_psi = fc * u["f"], fy * u["f"]
    nu_lb = n_u * u["N"]
    sq = math.sqrt(fc_psi)
    r = hw / lw
    ac = alpha_c(r)
    upper = 10.0 * sq * aw
    lower = 6.0 * sq * aw
    aci_raw = (ac * sq + rho * fy_psi) * aw
    wood_raw = rho * aw * fy_psi / 4.0
    barda = (8.0 * sq - 2.5 * sq * r + nu_lb / (4.0 * lw * tw) + rho * fy_psi) * tw * (0.6 * lw)
    if force_args:
        f_vw, f_be = a_be * u["N"], fy_be * u["N"]
    else:
        f_vw = rho * aw * fy_psi
        f_be = a_be * u["A"] * fy_be * u["f"]
    gw = (1.5 * sq * aw + 0.25 * f_vw + 0.20 * f_be + 0.40 * nu_lb) / math.sqrt(r)
    o = u["out"]
    return ShearCapacities(aci_raw * o, min(aci_raw, upper) * o, upper * o, wood_raw * o, lower * o,
                           min(max(wood_raw, lower), upper) * o, barda * o, gw * o, 3.0 * sq * aw * o, ac,
                           units, u["force"])


def bbcgen_curve(v_u: float, v_cr: float, ga: float, gamma_y: float = GAMMA_Y, gamma_f: float = GAMMA_FAIL,
                 steps: int = BBCGEN_STEPS) -> Tuple[np.ndarray, np.ndarray, int]:
    """The 22-point BBC of BBCGEN (D-NON-09): ``(gamma, V, yield index)``; ``ga`` = G A_W."""
    if not v_u > 0:
        raise PanelError(f"ultimate shear {v_u:g} must be > 0")
    if not 0 < v_cr < v_u:
        raise PanelError(f"cracking force {v_cr:g} must lie in (0, V_u = {v_u:g})")
    if ga <= 0:
        raise PanelError("G A_W must be > 0")
    g_cr = v_cr / ga
    if not g_cr < gamma_y < gamma_f:
        raise PanelError(f"strains must increase: gamma_cr {g_cr:.4g} < gamma_y {gamma_y:g} < failure {gamma_f:g}")
    k = np.arange(1, steps + 1)
    xi = k / steps
    g = np.concatenate(([g_cr], g_cr + xi * (gamma_y - g_cr), [gamma_f]))
    v = np.concatenate(([v_cr], v_cr + (v_u - v_cr) * (1.0 - (1.0 - xi) ** 2), [FAIL_FACTOR * v_u]))
    return g, v, steps + 1


# ======================================================================================
# Panel geometry
# ======================================================================================
@dataclass
class PanelGeometry:
    """Geometry of one panel (SHELL group) in global coordinates."""
    nodes: List[int]
    corners: Tuple[int, int, int, int]       # BL, BR, TR, TL
    e_h: np.ndarray                          # unit horizontal axis in the panel plane
    normal: np.ndarray
    length: float                            # L (corner-based width)
    height: float                            # H (corner-based height)
    width_extent: float                      # l_W: horizontal extent of the group nodes
    height_extent: float                     # h_W: vertical extent of the group nodes
    residual: float                          # max distance of the nodes from the plane
    vertical: bool
    corner_method: str                       # 'count' (node-connection counting) or 'geometry'
    warnings: List[str] = field(default_factory=list)


def element_nodes(nodes: Sequence[int]) -> List[int]:
    """Distinct non-zero nodes of an element (repeated nodes of degenerate elements removed)."""
    out: List[int] = []
    for n in nodes:
        n = int(n)
        if n and n not in out:
            out.append(n)
    return out


def plane_fit(xyz: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Least-squares plane of points: (centroid, unit normal, max |distance|)."""
    p = np.asarray(xyz, dtype=float)
    c = p.mean(axis=0)
    if len(p) < 3:
        raise PanelError("a plane needs at least 3 points")
    _, s, vt = np.linalg.svd(p - c)
    n = vt[-1]
    n = n / np.linalg.norm(n)
    # canonical sign: first significant component positive
    for v in n:
        if abs(v) > 1e-12:
            if v < 0:
                n = -n
            break
    return c, n, float(np.max(np.abs((p - c) @ n)))


def horizontal_axis(normal: np.ndarray) -> np.ndarray:
    """``e_h``: the unit horizontal vector of a vertical plane, oriented toward +X (or +Y for a wall
    normal to X), so that 'left' means the smaller global coordinate."""
    n = np.asarray(normal, dtype=float)
    e = np.array([-n[1], n[0], 0.0])
    norm = np.linalg.norm(e)
    if norm < 1e-12:
        raise PanelError("the panel plane is horizontal (a floor), not a wall")
    e /= norm
    if abs(e[0]) > 1e-9:
        if e[0] < 0:
            e = -e
    elif e[1] < 0:
        e = -e
    return e


def count_corners(elements: Sequence[Sequence[int]]) -> List[int]:
    """Node-connection counting of NONLINMOTDISP (spec 11 section 2.4): nodes used by exactly one element."""
    count: Dict[int, int] = {}
    for el in elements:
        for n in element_nodes(el):
            count[n] = count.get(n, 0) + 1
    return sorted(n for n, k in count.items() if k == 1)


def panel_geometry(elements: Sequence[Sequence[int]], xyz: Dict[int, Sequence[float]]) -> PanelGeometry:
    """Corners, axes and dimensions of a panel made of ``elements`` (node lists) with node coordinates ``xyz``."""
    nodes = sorted({n for el in elements for n in element_nodes(el)})
    if len(nodes) < 3:
        raise PanelError("the panel has fewer than 3 nodes")
    missing = [n for n in nodes if n not in xyz]
    if missing:
        raise PanelError(f"nodes {missing[:5]} of the panel are not defined")
    P = np.array([xyz[n] for n in nodes], dtype=float)
    c, n, res = plane_fit(P)
    size = float(np.max(np.ptp(P, axis=0)))
    warnings: List[str] = []
    vertical = abs(n[2]) <= VERTICAL_TOL
    if not vertical:
        raise PanelError(f"the panel plane is not vertical (normal ({n[0]:.3f}, {n[1]:.3f}, {n[2]:.3f}))")
    if res > COPLANAR_TOL * max(size, 1e-30):
        warnings.append(f"the shells are not coplanar: max distance from the plane {res:.3g} ({res / size:.2e} of "
                        "the panel size)")
    e_h = horizontal_axis(n)
    h = (P - c) @ e_h
    v = P[:, 2]
    pos = {nd: k for k, nd in enumerate(nodes)}
    cand = count_corners(elements)
    method = "count"
    if len(cand) != 4:
        method = "geometry"
        warnings.append(f"node-connection counting found {len(cand)} corner nodes instead of 4 (triangular, warped or "
                        "non-rectangular panel): the nodes nearest to the corners of the panel's bounding rectangle "
                        "are used")
        box = [(h.min(), v.min()), (h.max(), v.min()), (h.max(), v.max()), (h.min(), v.max())]
        cand = []
        for (bh, bv) in box:
            d = (h - bh) ** 2 + (v - bv) ** 2
            for k in np.argsort(d, kind="stable"):
                if nodes[int(k)] not in cand:
                    cand.append(nodes[int(k)])
                    break
    pts = sorted(cand, key=lambda nd: (v[pos[nd]], h[pos[nd]]))
    bottom = sorted(pts[:2], key=lambda nd: h[pos[nd]])
    top = sorted(pts[2:], key=lambda nd: h[pos[nd]])
    bl, br = bottom
    tl, tr = top
    hv = {nd: (h[pos[nd]], v[pos[nd]]) for nd in (bl, br, tr, tl)}
    L = 0.5 * ((hv[br][0] - hv[bl][0]) + (hv[tr][0] - hv[tl][0]))
    H = 0.5 * ((hv[tl][1] - hv[bl][1]) + (hv[tr][1] - hv[br][1]))
    if L <= 0 or H <= 0:
        raise PanelError(f"degenerate panel corners (L = {L:.4g}, H = {H:.4g})")
    return PanelGeometry(nodes, (bl, br, tr, tl), e_h, n, float(L), float(H), float(np.ptp(h)), float(np.ptp(v)),
                         res, vertical, method, warnings)


def panel_deformations(corner_disp: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]], e_h: Sequence[float],
                       length: float, height: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Shear strain, curvature and axial strain histories of a panel (module docstring formulas).

    ``corner_disp[label] = (ux, uy, uz)`` global displacement histories of corner ``label`` (BL BR TR TL)."""
    eh = np.asarray(e_h, dtype=float)
    u, w = {}, {}
    for lab in CORNERS:
        ux, uy, uz = (np.asarray(a, dtype=float) for a in corner_disp[lab])
        u[lab] = eh[0] * ux + eh[1] * uy + eh[2] * uz
        w[lab] = uz
    L, H = float(length), float(height)
    gamma = 0.5 * ((u["TL"] - u["BL"]) + (u["TR"] - u["BR"])) / H + 0.5 * ((w["BR"] - w["BL"]) + (w["TR"] - w["TL"])) / L
    theta_b = (w["BR"] - w["BL"]) / L
    theta_t = (w["TR"] - w["TL"]) / L
    kappa = (theta_t - theta_b) / H
    eps_v = 0.5 * ((w["TL"] - w["BL"]) + (w["TR"] - w["BR"])) / H
    return gamma, kappa, eps_v


# ======================================================================================
# Geometry helpers of the panel-model commands (WALLFLR, PANELIZE, EDGE)
# ======================================================================================
def shell_plane(xyz: np.ndarray) -> Tuple[np.ndarray, float]:
    """Unit normal (canonical sign) and offset ``n . x`` of the plane of one shell element."""
    p = np.asarray(xyz, dtype=float)
    if len(p) >= 4:
        n = np.cross(p[2] - p[0], p[3] - p[1])
    else:
        n = np.cross(p[1] - p[0], p[2] - p[0])
    nn = np.linalg.norm(n)
    if nn == 0:
        raise PanelError("degenerate shell element (zero area)")
    n = n / nn
    for v in n:
        if abs(v) > 1e-9:
            if v < 0:
                n = -n
            break
    return n, float(n @ p.mean(axis=0))


def plane_clusters(normals: np.ndarray, offsets: np.ndarray, angle_tol: float = 1e-3,
                   offset_tol: float = 1e-6) -> np.ndarray:
    """Cluster element planes: same plane when the normals are parallel within ``angle_tol`` (radians,
    as |n1 x n2|) and the offsets differ by at most ``offset_tol``.  Returns a cluster label per element
    (labels in order of first appearance)."""
    N = np.asarray(normals, dtype=float)
    d = np.asarray(offsets, dtype=float)
    labels = -np.ones(len(N), dtype=int)
    reps: List[Tuple[np.ndarray, float]] = []
    for i in range(len(N)):
        for k, (rn, rd) in enumerate(reps):
            if np.linalg.norm(np.cross(N[i], rn)) <= angle_tol and abs(d[i] * float(N[i] @ rn) - rd) <= offset_tol:
                labels[i] = k
                break
        else:
            labels[i] = len(reps)
            reps.append((N[i], d[i]))
    return labels


def connected_components(n: int, pairs: Sequence[Tuple[int, int]]) -> np.ndarray:
    """Component label (smallest member index) of n items linked by ``pairs``."""
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for a, b in pairs:
        ra, rb = find(int(a)), find(int(b))
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    return np.array([find(i) for i in range(n)], dtype=int)


def element_edges(nodes: Sequence[int]) -> List[Tuple[int, int]]:
    """Edges (sorted node pairs) of a shell element (3 or 4 distinct nodes)."""
    nd = element_nodes(nodes)
    return [tuple(sorted((nd[k], nd[(k + 1) % len(nd)]))) for k in range(len(nd))]


def edge_cells(elem_xyz: Sequence[np.ndarray], e_h: np.ndarray, boundary_edges: Sequence[Tuple[np.ndarray, np.ndarray]],
               skip: Sequence[bool], line_tol: float) -> Tuple[np.ndarray, List[Tuple[float, float]]]:
    """EDGE cell assignment (spec 11 section 2.5): each kept boundary edge defines an infinite line in
    the panel plane (local h, v coordinates); collinear lines are merged; every element is labelled by
    the side of each line its centroid lies on.  ``skip[k]`` drops edge k (flags X/Y/Z = 1).

    Returns ``(cell key per element, cell centroids (h, v) are not needed)`` -- the key is an integer
    label of the sign pattern (labels in order of first appearance) and the list of distinct lines."""
    eh = np.asarray(e_h, dtype=float)
    lines: List[Tuple[np.ndarray, float]] = []          # (unit normal in (h, v), offset)
    for (a, b), sk in zip(boundary_edges, skip):
        if sk:
            continue
        pa = np.array([float(a @ eh), float(a[2])])
        pb = np.array([float(b @ eh), float(b[2])])
        t = pb - pa
        ln = np.linalg.norm(t)
        if ln == 0:
            continue
        nrm = np.array([-t[1], t[0]]) / ln
        for v in nrm:
            if abs(v) > 1e-9:
                if v < 0:
                    nrm = -nrm
                break
        off = float(nrm @ pa)
        if any(np.linalg.norm(nrm - q) <= 1e-6 and abs(off - o) <= line_tol for q, o in lines):
            continue
        lines.append((nrm, off))
    cent = np.array([[float(np.mean(p, axis=0) @ eh), float(np.mean(p, axis=0)[2])] for p in elem_xyz])
    keys: Dict[Tuple[int, ...], int] = {}
    labels = np.zeros(len(cent), dtype=int)
    for i, c in enumerate(cent):
        sig = tuple(int(np.sign(round((float(q @ c) - o) / max(line_tol, 1e-300), 6))) for q, o in lines)
        labels[i] = keys.setdefault(sig, len(keys))
    return labels, [(float(q[0]), float(o)) for q, o in lines]
