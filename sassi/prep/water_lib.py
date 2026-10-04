"""Water modelling algorithms: FILLPOOL, REFINEMODEL, LISTPOOLINTER and MERGEPOOL (manual section 9.16;
requirements section 3.4.N; spec 11 section 1; decisions D-WAT-01, D-CHK-07; VP-54W, VP-W1).

Pure functions on :class:`sassi.model.SSIModel`; the command handlers of
:mod:`sassi.prep.commands.water` call them and print the collected :class:`Notes`.  Every function
validates its input before it changes anything and raises :class:`GenerationError` (the model is
then unchanged).

Background for the structural engineer: water as "soft solids"
---------------------------------------------------------------
SASSI has no fluid element.  ACS SASSI models the water of a pool (spent-fuel pool, tank) with
ordinary 8-node SOLID elements that are given the properties of water:

* the **bulk modulus** of water, K = 2.2 GPa (P-wave velocity ~1,480 m/s ~ 4,870 ft/s), and its
  mass density (unit weight 9.81 kN/m3 = 62.4 pcf);
* a **near-zero shear modulus** G.  An inviscid fluid carries no shear stress; a solid with G -> 0
  carries only the pressure ``p = -K div u``.  Its equation of motion then reduces to the linear
  acoustic (wave) equation in displacement form, ``rho u_tt = K grad(div u)``.

For a pool whose walls move at frequencies well below the acoustic (compression) frequencies of
the water, ``f << c/(4H)`` (about 74 Hz for H = 5 m), the water is practically incompressible and
the solution of that equation is the classical **potential flow** of Westergaard (1933) and Housner
(1963): the part of the water that is forced to move with the walls -- the **impulsive** (added)
mass m_i -- pushes on the walls with the hydrodynamic pressure, and the rest of the water lags
behind because the free surface (p = 0) lets it slosh up and down.  The FE "water solids"
reproduce this impulsive behaviour (VP-W1: the base shear of a rigid tank converges to the exact
potential-flow value m_i a; Housner's closed form tanh(sqrt(3) L/H)/(sqrt(3) L/H) is a slightly
conservative engineering approximation of it).

Limits of the approximation (stated in docs/user/WATER.md):

1. **No sloshing (convective) response.**  Surface gravity waves need the restoring force rho g of
   the free surface, which the elastic solid does not have; the convective mass and the sloshing
   frequencies (a fraction of a hertz for pools) are missing.  Add them separately (Housner's
   convective spring-mass) when they matter (they usually do not for the structure, but they do for
   freeboard and the pool walls at low frequencies).
2. **Spurious low-frequency modes.**  Because G is not exactly zero the water has many "shear"
   modes at very low frequencies (with G = 1e-8 K, below about 0.05 Hz for metre-sized pools); far
   above them the water behaves as an inviscid fluid.  Do not interpret the water response below
   ~0.5 Hz.
3. **Volumetric locking.**  A trilinear SOLID with a nearly incompressible material "locks" (it
   cannot deform at constant volume and becomes far too stiff) unless it is enriched: the water
   needs the Wilson-Taylor **incompatible modes**, ``MOPT,0`` (D-ELM-02).  FILLPOOL and MERGEPOOL
   therefore set MOPT ``<incomp>`` = 0 and say so.  Without them the water moves rigidly with the
   walls (the whole mass instead of the impulsive mass, and spurious resonances).
4. **Interface springs.**  The water nodes are separate from the wall nodes and are attached by
   springs: stiff along the wall normal (``Stiff``) so that the water cannot penetrate or separate,
   free along the wall (``stiff2`` = 0) so that the inviscid water can slip.  The springs must be
   stiff compared with the water itself; FILLPOOL prints the frequency of the water mass on the
   springs, which must be well above the frequencies of interest.

Decisions (requirements section 7; spec 11 open questions OQ-1 ... OQ-6)
----------------------------------------------------------------------
* D-WAT-01 water material: K = 2.2 GPa and unit weight 9.81 kN/m3 (SI, ``g < 20``) or 62.4 pcf =
  0.0624 kcf (British, kip-ft); damping 0.5 %.  The shear modulus is **G = 1e-8 K** (Vs = 1e-4 Vp)
  instead of the nu = 0.49 of D-WAT-01: with nu = 0.49 (G = 0.02 K) the "water" is an elastic solid
  with shear modes in the seismic band (about 10 Hz for a 5 m deep pool) and moves rigidly with the
  walls below them; the lead instruction for this package (near-zero shear modulus) and the VP-W1
  study (docs/user/WATER.md, section "Choice of the shear modulus") require the smaller value.  The
  material is a type-3 (Vp, Vs) M-table entry that the user can redefine with ``M``.
* Oblique walls (OQ-2): SASSI spring constants (SC) are uncoupled global constants.  For an
  interface whose normals are the global axes the SPRING is exact (``Stiff`` along each normal,
  ``stiff2`` in the other directions).  For oblique walls the exact coupled 3x3 matrix
  ``Stiff P_N + stiff2 (I - P_N)`` (``P_N`` = projector on the wall normals) is replaced by its
  diagonal (OQ-2 (b)) -- the zero-length 2-node GENERAL elements of D-WAT-01 are rejected by
  CHECK Error 8 (D-CHK-07), see the work-package report.  The diagonal keeps the normal restraint
  but adds a tangential one (the water sticks to oblique walls).
* Corner and edge nodes (OQ-2): the normals of all wetted faces at the node are combined
  (``P_N`` spans them); floor nodes get Z springs.  Normals of adjacent boundary edges that differ
  by less than 30 degrees are averaged (a polygonal approximation of a curved wall is smooth).
* Spring damping 0 (OQ-3); interface-area shells (OQ-5): SHELL elements on the wall-side nodes of
  every wetted water face, ETYPE 1, a dummy material with E = 1e-6 K_water, nu = 0.3, zero weight
  and zero damping, thickness 1/100 of the mean wetted-face edge (negligible stiffness, no mass;
  delete the group with GDEL when it is not needed).
* Z-levels (OQ-1): single-linkage clustering of the node elevations with gap ``Sensitivity`` (the
  geometric tolerance when 0), level = cluster mean, as EXCAV (D-MDL-19).
* The floor template (the plan mesh of the water) is the lowest band of horizontal SHELL/TSHELL
  elements and of free, upward-facing SOLID faces; the water fills from that level to the level
  ``EmptyLevels`` below the highest level.
* Rotations: SPRING elements carry six DOFs per node.  The rotations of the water nodes (SOLID
  nodes) and of SOLID-wall interface nodes are restrained by nothing else, so FILLPOOL fixes them
  (``D ... ROTX ROTY ROTZ``); otherwise the system is singular (EDU-06).
* Every water SOLID gets ETYPE 1 explicitly (an implicit ETYPE 0 SOLID below grade would be
  classified as excavated soil, D-HOU-01).
* REFINEMODEL (OQ-6): SHELL, TSHELL and PLANE quadrilaterals (4 distinct nodes) are split into 4 and
  SOLID hexahedra (8 distinct nodes) into 8 at the edge midpoints, face centres and body centre
  (shared between neighbours); triangles, prisms/pyramids and the other types are unchanged.  New
  nodes are numbered after the last node; they inherit an INT code and a fixity when every parent
  corner node has it; the elements of a refined group are renumbered 1..n (children of a parent
  consecutive); EOUT requests and the session cuts follow.  Loads and masses are not redistributed.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from ..conventions import G_BRITISH, G_SI, unit_system
from ..model import SSIModel
from ..model.entities import Element, Group, Material, SpringProp
from ..model.options import Field, Record, register_record_type
from ..model.values import fmt_num
from .check import BEAMS, GENERAL, PLANE, SHELL, SOLID, SOLID_FACES, SPRING, TSHELL, ModelView, newell_normal
from .generation_lib import (HORIZONTAL_COS, GenerationError, Notes, _ordered_unique, _short, check_node_limit,
                             cluster_values, lowest_faces, signed_area, template_grid)

# ======================================================================================
# Water properties (D-WAT-01, amended: G = 1e-8 K)
# ======================================================================================
#: bulk modulus of water (Pa)
WATER_BULK_PA = 2.2e9
#: 1 ksf in Pa (1000 lbf / ft^2)
PA_PER_KSF = 1000.0 * 4.4482216152605 / 0.09290304
#: unit weight of water: kN/m^3 (SI) and kip/ft^3 (British, 62.4 pcf)
WATER_WEIGHT = {"SI": 9.81, "BS": 0.0624}
#: shear modulus of the water solids as a fraction of the bulk modulus (Vs = 1e-4 Vp)
WATER_G_OVER_K = 1.0e-8
#: hysteretic damping of the water solids (D-WAT-01)
WATER_DAMPING = 0.005
#: dummy interface-area shells: E as a fraction of the water bulk modulus, Poisson's ratio
SHELL_AREA_E_OVER_K = 1.0e-6
SHELL_AREA_NU = 0.3
#: boundary-edge normals closer than this angle are averaged (smooth curved wall)
SMOOTH_ANGLE_DEG = 30.0
#: default interface stiffness of FILLPOOL (manual 9.16.1)
DEFAULT_STIFF = 1.0e6
#: group titles written by FILLPOOL
TITLES = {"water": "FILLPOOL water", "springs": "FILLPOOL water-wall interface springs",
          "shells": "FILLPOOL interface areas"}

WALL_TYPES = (SHELL, TSHELL, SOLID)


@dataclass
class WaterProperties:
    """Water material in the model's units (D-WAT-01)."""
    units: str          # 'SI' (kN, m) or 'BS' (kip, ft)
    K: float            # bulk modulus
    G: float            # shear modulus (G = WATER_G_OVER_K K)
    weight: float       # unit weight (force / volume)
    rho: float          # mass density = weight / g
    vp: float           # sqrt((K + 4G/3) / rho)
    vs: float           # sqrt(G / rho)
    damping: float
    recognised: bool    # gravity close to 9.80665 or 32.174

    @property
    def label(self) -> str:
        return "kN-m" if self.units == "SI" else "kip-ft"


def water_properties(gravity: float) -> WaterProperties:
    """Water of D-WAT-01 in the unit system detected from the HOUSE gravity (D-CNV-08: British, kip-ft,
    when g > 20, else SI, kN-m): K = 2.2 GPa, unit weight 9.81 kN/m3 or 0.0624 kcf, G = 1e-8 K,
    damping 0.5 %.  ``recognised`` is False when g differs by more than 5 % from 9.80665 / 32.174."""
    g = float(gravity)
    if not g > 0:
        raise GenerationError("the gravity (HOUSE <gravity>, GRAVITY) must be > 0 to define the water mass")
    units = unit_system(g)
    K = WATER_BULK_PA / 1000.0 if units == "SI" else WATER_BULK_PA / PA_PER_KSF
    G = WATER_G_OVER_K * K
    w = WATER_WEIGHT[units]
    rho = w / g
    ref = G_SI if units == "SI" else G_BRITISH
    return WaterProperties(units, K, G, w, rho, math.sqrt((K + 4.0 * G / 3.0) / rho), math.sqrt(G / rho),
                           WATER_DAMPING, abs(g - ref) <= 0.05 * ref)


# ======================================================================================
# Pool data record (the "stored interface metadata" LISTPOOLINTER needs)
# ======================================================================================
@register_record_type
class PoolDataRecord(Record):
    """``POOLDATA`` -- what FILLPOOL created in a pool model (SASSI-EDU storage record; WRITE writes it after
    the groups, so a FILLPOOL pool survives WRITE -> INP and LISTPOOLINTER / MERGEPOOL still work)."""
    COMMAND = "POOLDATA"
    FIELDS = (Field("water", int, 0, "group of the water SOLIDs"),
              Field("springs", int, 0, "group of the interface SPRINGs"),
              Field("shells", int, 0, "group of the interface-area SHELLs (0 none)"),
              Field("stiff", float, DEFAULT_STIFF, "FILLPOOL <Stiff>"),
              Field("sensitivity", float, 0.0, "FILLPOOL <Sensitivity>"),
              Field("empty", int, 0, "FILLPOOL <EmptyLevels>"),
              Field("shellarea", int, -1, "FILLPOOL <ShellArea>"),
              Field("offset", int, 0, "resolved FILLPOOL <offset> (water nodes start at offset + 1)"),
              Field("stiff2", float, 0.0, "FILLPOOL <stiff2>"),
              Field("watermat", int, 0, "M-table entry of the water"),
              Field("zfloor", float, 0.0, "elevation of the pool floor (bottom of the water)"),
              Field("zsurface", float, 0.0, "elevation of the water surface"))


def pool_record(m: SSIModel) -> Optional[PoolDataRecord]:
    rec = m.options.record("POOLDATA")
    return rec if isinstance(rec, Record) else None


# ======================================================================================
# Analytical references of the impulsive (hydrodynamic) response of a rigid rectangular tank
# ======================================================================================
ZETA3 = 1.2020569031595942853997   # Apery's constant zeta(3)


def _lambdas(H: float, nterms: int) -> np.ndarray:
    n = np.arange(nterms, dtype=float)
    return (2.0 * n + 1.0) * math.pi / (2.0 * H)


def impulsive_ratio_exact(L: float, H: float, nterms: int = 400) -> float:
    """Impulsive mass ratio m_i/m of a **rigid rectangular tank** (length 2L in the direction of the motion,
    water depth H): incompressible inviscid water, free surface p = 0 (no gravity waves), slip walls.

    Potential flow ``phi = sum A_n sinh(l_n x) cos(l_n z)`` (z from the floor, ``l_n = (2n+1) pi/(2H)``)
    with ``d phi/dx = v`` on the walls x = +-L gives the wall pressure
    ``p(z) = rho a sum 2 (-1)^n tanh(l_n L) cos(l_n z) / (H l_n^2)`` and the base shear
    ``F = 2 b int p dz = rho a b sum 4 tanh(l_n L)/(H l_n^3)``, so::

        m_i / m = 2/(L H^2) sum_n tanh(l_n L) / l_n^3

    The slowly converging part ``sum 1/l_n^3 = (2H/pi)^3 (7/8) zeta(3)`` is summed in closed form and only
    the exponentially small remainder ``(1 - tanh(l_n L))/l_n^3`` numerically.  Limits: 1 for L/H -> 0,
    0.5428 H/L for L/H -> infinity (Westergaard's dam, exact series).  The series satisfies
    ``r(L/H) + r(H/L) = 1`` (checked numerically to round-off; hence exactly 1/2 for L = H)."""
    if not (L > 0 and H > 0):
        raise ValueError("L and H must be > 0")
    # enough terms for l_n L >= 20 (remainder terms below e^-40), at most 5e6
    nterms = int(min(max(nterms, math.ceil(20.0 * H / (math.pi * L)) + 1), 5_000_000))
    lam = _lambdas(H, nterms)
    s_inf = (2.0 * H / math.pi) ** 3 * 0.875 * ZETA3
    x = np.minimum(lam * L, 350.0)
    rem = np.sum((2.0 / (np.exp(2.0 * x) + 1.0)) / lam ** 3)     # 1 - tanh(x) = 2/(e^{2x} + 1)
    return float(2.0 / (L * H * H) * (s_inf - rem))


def impulsive_height_exact(L: float, H: float, nterms: int = 20000) -> float:
    """Height h_i/H of the resultant of the impulsive **wall** pressures above the floor (EBP, excluding the
    base pressure) for the tank of :func:`impulsive_ratio_exact`: ``M/F`` with
    ``int_0^H cos(l z) z dz = H (-1)^n / l - 1/l^2``."""
    lam = _lambdas(H, nterms)
    sgn = np.where(np.arange(nterms) % 2 == 0, 1.0, -1.0)
    t = np.tanh(lam * L)
    F = np.sum(2.0 * t / (H * lam ** 3))
    M = np.sum(2.0 * sgn * t / (H * lam ** 2) * (H * sgn / lam - 1.0 / lam ** 2))
    return float(M / F / H)


def housner_impulsive_ratio(L: float, H: float) -> float:
    """Housner (1963) rectangular tank: ``m_i/m = tanh(sqrt(3) L/H) / (sqrt(3) L/H)`` (L = half length)."""
    x = math.sqrt(3.0) * L / H
    return math.tanh(x) / x


#: Housner (1963) height of the impulsive force, wall pressures only (EBP)
HOUSNER_HEIGHT_EBP = 3.0 / 8.0


def housner_convective(L: float, H: float, g: float) -> Tuple[float, float, float]:
    """Housner (1963) first sloshing (convective) mode of a rigid rectangular tank -- the part FILLPOOL does
    **not** model: ``(m_c/m, f_c in Hz, h_c/H)`` with ``k = sqrt(5/2)``::

        m_c/m = (1/3) k (L/H) tanh(k H/L)        w_c^2 = (k g / L) tanh(k H/L)
        h_c/H = 1 - (cosh(k H/L) - 1) / ((k H/L) sinh(k H/L))      (EBP)

    (L = half length in the direction of the motion; the exact first sloshing frequency has pi/2 = 1.571
    in place of k = 1.581)."""
    k = math.sqrt(2.5)
    x = k * H / L
    mc = k / 3.0 * (L / H) * math.tanh(x)
    w = math.sqrt(k * g / L * math.tanh(x))
    hc = 1.0 - (math.cosh(x) - 1.0) / (x * math.sinh(x)) if x < 300 else 1.0 - 1.0 / x
    return mc, w / (2.0 * math.pi), hc


def sloshing_frequency_exact(L: float, H: float, g: float) -> float:
    """First sloshing frequency (Hz) of a rectangular tank of length 2L and depth H (linear potential theory):
    ``w^2 = (pi g / 2L) tanh(pi H / 2L)``."""
    return math.sqrt(math.pi * g / (2.0 * L) * math.tanh(math.pi * H / (2.0 * L))) / (2.0 * math.pi)


def westergaard_ratio(L: float, H: float) -> float:
    """Westergaard (1933) dam pressure ``p = 7/8 rho a sqrt(H y)`` on both end walls of a tank of length 2L:
    ``F = 2 (7/12) rho a H^2 b`` per tank width b, i.e. ``m_i/m = (7/12) H/L`` -- valid for long tanks
    (L >> H) only; it ignores the opposite wall."""
    return 7.0 / 12.0 * H / L


# ======================================================================================
# Small helpers
# ======================================================================================
def _wall_refs(v: ModelView):
    return [r for r in v.elems if r.type in WALL_TYPES]


def _solid_face_polys(ns: Sequence[int]) -> List[List[int]]:
    pad = list(ns) + [0] * (8 - len(ns))
    return [_ordered_unique(pad[k - 1] for k in f) for f in SOLID_FACES]


def _free_upward_solid_faces(v: ModelView, refs) -> List[List[int]]:
    """Free (unshared) SOLID faces whose outward normal points up (``n_z >= HORIZONTAL_COS |n|``)."""
    count: Dict[Tuple[int, ...], int] = {}
    owner: Dict[Tuple[int, ...], Tuple[List[int], np.ndarray]] = {}
    for r in refs:
        if r.type != SOLID or not v.defined(r.elem.nodes):
            continue
        ns = [n for n in r.elem.nodes if n]
        cen = np.mean([v.P[n] for n in dict.fromkeys(ns)], axis=0)
        for p in _solid_face_polys(r.elem.nodes):
            if len(p) < 3:
                continue
            key = tuple(sorted(p))
            count[key] = count.get(key, 0) + 1
            owner[key] = (p, cen)
    out = []
    for key, k in count.items():
        if k != 1:
            continue
        p, cen = owner[key]
        X = [v.P[n] for n in p]
        nv = newell_normal(X)
        if float(np.dot(nv, np.mean(X, axis=0) - cen)) < 0:
            nv = -nv
        ln = float(np.linalg.norm(nv))
        if ln > 0 and nv[2] >= HORIZONTAL_COS * ln:
            out.append(p)
    return out


def _horizontal_shells(v: ModelView, refs) -> List[List[int]]:
    out = []
    for r in refs:
        if r.type not in (SHELL, TSHELL) or not v.defined(r.elem.nodes):
            continue
        p = _ordered_unique(r.elem.nodes)
        if len(p) < 3:
            continue
        nv = newell_normal([v.P[n] for n in p])
        ln = float(np.linalg.norm(nv))
        if ln > 0 and abs(nv[2]) >= HORIZONTAL_COS * ln:
            out.append(p)
    return out


def _boundary_normals(xy: np.ndarray, cells: Sequence[Tuple[int, ...]]) -> Tuple[Dict[int, List[np.ndarray]],
                                                                                List[Tuple[int, int]]]:
    """Outward plan normals of the template boundary at every boundary point.

    Boundary edges are cell edges used by one cell; for a counter-clockwise cell the edge a -> b has the
    outward normal ``(dy, -dx)/|ab|``.  Normals of the edges meeting at a point that differ by less than
    :data:`SMOOTH_ANGLE_DEG` are averaged (smooth wall), the others are kept (corner).  Returns
    ``({point: [unit normals (3,)]}, boundary edges (a, b) in cell orientation)``."""
    use: Dict[Tuple[int, int], int] = {}
    oriented: Dict[Tuple[int, int], Tuple[int, int]] = {}
    for c in cells:
        k = len(c)
        for i in range(k):
            a, b = c[i], c[(i + 1) % k]
            key = (min(a, b), max(a, b))
            use[key] = use.get(key, 0) + 1
            oriented[key] = (a, b)
    edges = [oriented[k] for k in sorted(use) if use[k] == 1]
    raw: Dict[int, List[np.ndarray]] = {}
    for a, b in edges:
        d = xy[b] - xy[a]
        ln = float(np.hypot(d[0], d[1]))
        if ln == 0:
            continue
        nrm = np.array([d[1] / ln, -d[0] / ln, 0.0])
        raw.setdefault(a, []).append(nrm)
        raw.setdefault(b, []).append(nrm)
    cos_s = math.cos(math.radians(SMOOTH_ANGLE_DEG))
    out: Dict[int, List[np.ndarray]] = {}
    for p, lst in raw.items():
        groups: List[List[np.ndarray]] = []
        for nv in lst:
            for g in groups:
                mean = np.sum(g, axis=0)
                mean = mean / np.linalg.norm(mean)
                if float(np.dot(mean, nv)) >= cos_s:
                    g.append(nv)
                    break
            else:
                groups.append([nv])
        out[p] = [np.sum(g, axis=0) / np.linalg.norm(np.sum(g, axis=0)) for g in groups]
    return out, edges


def normal_projector(normals: Sequence[np.ndarray], rtol: float = 1e-6) -> np.ndarray:
    """Orthogonal projector P_N (3x3) on the span of the unit normals (SVD rank with ``rtol``)."""
    N = np.asarray(normals, float).reshape(-1, 3)
    if N.size == 0:
        return np.zeros((3, 3))
    _, s, vt = np.linalg.svd(N)
    r = int(np.sum(s > rtol * s[0]))
    B = vt[:r]
    return B.T @ B


def _round_sig(x: float, digits: int = 10) -> float:
    """``x`` rounded to ``digits`` significant digits (readable M-table values in WRITE)."""
    return float(f"{float(x):.{digits}g}")


#: element types that give a node rotational DOFs (and normally rotational stiffness)
ROTATION_TYPES = (BEAMS, SHELL, TSHELL, SPRING, GENERAL)


def fix_interface_rotations(m: SSIModel, v: ModelView, walls: Iterable[int]) -> Tuple[List[int], List[int]]:
    """Restrain the rotations that the interface SPRINGs (six DOFs per node, no rotational stiffness) add at
    wall nodes where nothing else restrains them (EDU-06; the rules of FIXROT, which no longer sees these
    nodes as solid-only / shell-only once the springs exist):

    * wall nodes of SOLID elements only: ROTX, ROTY, ROTZ fixed (FIXSLDROT rule);
    * wall nodes of coplanar SHELL elements only whose normal is a global axis: the drilling rotation about
      that axis fixed (FIXROT rule; oblique shells need FIXSHLROT springs, nodes where shells meet at an angle
      restrain each other).

    ``v`` is the view of ``m`` *before* the springs are added.  Returns (nodes with all rotations fixed now,
    nodes with the drilling rotation fixed now)."""
    eps = 1e-4
    solid_fixed: List[int] = []
    drill_fixed: List[int] = []
    for w in sorted(set(walls)):
        types = {r.type for r in v.node_elems.get(w, [])}
        fx = m.nodes[w].fix
        if not types & set(ROTATION_TYPES):
            if fx[3:6] != [1, 1, 1]:
                solid_fixed.append(w)
            fx[3:6] = [1, 1, 1]
        elif types == {SHELL}:
            nrm = v.coplanar_normal(w, eps)
            if nrm is None:
                continue
            axis = [k for k in range(3) if abs(nrm[k]) >= 1.0 - eps]
            if axis and not fx[3 + axis[0]]:
                fx[3 + axis[0]] = 1
                drill_fixed.append(w)
    return solid_fixed, drill_fixed


def interface_constants(P: np.ndarray, stiff: float, stiff2: float, atol: float = 1e-9) -> Tuple[Tuple[float, float, float], bool]:
    """Uncoupled spring constants (kx, ky, kz) of an interface with normal projector P and whether they are
    exact (P diagonal: every normal is a global axis).  Exact: ``Stiff`` along the normals, ``stiff2``
    across; oblique: the diagonal of ``Stiff P + stiff2 (I - P)`` (OQ-2 (b))."""
    off = P - np.diag(np.diag(P))
    exact = bool(np.max(np.abs(off)) <= atol)
    d = np.clip(np.diag(P), 0.0, 1.0)
    if exact:
        d = np.where(d > 0.5, 1.0, 0.0)
    k = stiff * d + stiff2 * (1.0 - d)
    return (float(k[0]), float(k[1]), float(k[2])), exact


# ======================================================================================
# FILLPOOL (manual 9.16.1; spec 11 section 1.1)
# ======================================================================================
@dataclass
class FillPoolResult:
    water_group: int
    spring_group: int
    shell_group: int
    water_material: int
    shell_material: int
    offset: int
    levels: List[float]                 # L0 (floor) ... Ln (top of the walls)
    filled: int                         # number of filled level intervals
    z_floor: float
    z_surface: float
    water_nodes: List[int]
    elements: int
    template_points: int
    cells: int
    pairs: List[Tuple[int, int]]        # (wall node, water node) of every spring
    spring_props: Dict[int, Tuple[float, float, float]]
    oblique: int
    unconnected: List[int]              # water nodes on the wetted boundary without a wall node
    multiple: List[int]                 # water nodes with several coincident wall nodes
    internal: List[int]                 # water nodes inside the water that coincide with a wall node
    surface: List[int]                  # interior nodes of the water surface that coincide with a wall node
    shells: int
    wall_rotations_fixed: List[int]     # SOLID-wall interface nodes: ROTX/ROTY/ROTZ fixed
    wall_drilling_fixed: List[int]      # coplanar SHELL-wall interface nodes: drilling rotation fixed
    volume: float
    mass: float
    wetted_area: float
    water: WaterProperties
    interface_freq: Dict[str, float]
    mopt_changed: bool
    #: sloshing NOT modelled -- rectangular-tank estimate from the plan extents per horizontal direction:
    #: {'x': (half length, first sloshing frequency Hz, Housner convective mass ratio m_c/m), 'y': ...}
    sloshing: Dict[str, Tuple[float, float, float]] = field(default_factory=dict)


def resolve_offset(m: SSIModel, offset: int) -> int:
    """FILLPOOL ``<offset>``: <= 0 uses the largest node number of the pool model; a positive value smaller
    than it is an error (manual 9.16.1).  Water nodes are numbered from offset + 1."""
    nmax = max(m.nodes) if m.nodes else 0
    if offset <= 0:
        return nmax
    if offset < nmax:
        raise GenerationError(f"<offset> {offset} is smaller than the largest node number {nmax} of the pool model "
                              f"(use 0 or a value >= {nmax}; to import the water into the original model use a "
                              f"value >= the original model's last node number)")
    return int(offset)


def _validate_pool(m: SSIModel, v: ModelView, notes: Notes):
    if not m.nodes:
        raise GenerationError("the active model has no nodes: activate the pool sub-model (walls and floor)")
    refs = _wall_refs(v)
    other = sorted({r.type for r in v.elems if r.type not in WALL_TYPES})
    if other:
        from ..conventions import ELEMENT_TYPE_NAMES
        names = ", ".join(ELEMENT_TYPE_NAMES.get(t, str(t)) for t in other)
        raise GenerationError(f"the pool model must contain only the walls and floor of one pool (SHELL/TSHELL or "
                              f"SOLID elements); it also has {names} elements"
                              + (" (already filled? delete the FILLPOOL groups first)" if SPRING in other else ""))
    if not refs:
        raise GenerationError("the pool model has no SHELL, TSHELL or SOLID elements (walls and floor)")
    rec = pool_record(m)
    if rec is not None and int(rec.water) in m.groups:
        raise GenerationError(f"the pool is already filled (POOLDATA: water group {rec.water}); delete the "
                              f"FILLPOOL groups (GDEL) before filling it again")
    bad = [r for r in refs if not v.defined(r.elem.nodes)]
    if bad:
        raise GenerationError(f"{len(bad)} wall/floor elements reference undefined nodes (CHECK Error 41), e.g. "
                              f"group {bad[0].group} element {bad[0].id}")
    exc = [r for r in refs if r.type == SOLID and r.elem.etype == 2]
    if exc:
        raise GenerationError(f"{len(exc)} SOLID elements have ETYPE 2 (excavated soil): pool walls are structure")
    implicit = [r for r in refs if r.type == SOLID and r.elem.etype == 0 and r.etype == 2]
    if implicit:
        notes.warn(f"{len(implicit)} SOLID wall/floor elements with the implicit ETYPE 0 lie below the ground "
                   f"elevation and are classified as excavated soil (D-HOU-01): give them ETYPE 1")
    return refs


def fill_pool(m: SSIModel, stiff: float = DEFAULT_STIFF, sensitivity: float = 0.0, empty_levels: int = 0,
              shell_area: int = -1, offset: int = -1, stiff2: float = 0.0,
              notes: Optional[Notes] = None) -> FillPoolResult:
    """FILLPOOL,<Stiff>,<Sensitivity>,<EmptyLevels>,<ShellArea>,<offset>,<stiff2> on pool model ``m``.

    Steps (spec 11 section 1.1 with the decisions of the module docstring):

    1. validation: only SHELL/TSHELL/SOLID elements; ``offset`` resolved (:func:`resolve_offset`);
    2. Z-levels of every wall/floor node, clustered with ``Sensitivity`` (level = cluster mean);
    3. floor template: lowest band of horizontal shells and free upward SOLID faces (its level is L0);
    4. water SOLIDs: every template cell extruded between consecutive levels L0 ... L(n - EmptyLevels)
       (hexahedron for a quad, prism for a triangle); nodes numbered from offset + 1, bottom-up, level by
       level; one SOLID group, ETYPE 1, a new water material (:func:`water_properties`);
    5. interface: every water node of the wetted boundary (floor level and the lateral template boundary;
       not the free surface) that coincides with a wall node (distance <= max(Sensitivity, tol)) gets a
       SPRING (wall node I, water node J) with the constants of :func:`interface_constants`; rotations of
       the water nodes and of SOLID-wall interface nodes are fixed;
    6. ``ShellArea`` = 1: dummy SHELLs on the wall-side nodes of every wetted water face;
    7. the POOLDATA record and MOPT ``<incomp>`` = 0 (incompatible modes, needed by the water solids).
    """
    notes = notes if notes is not None else Notes()
    stiff = float(stiff)
    stiff2 = float(stiff2)
    if not stiff > 0:
        raise GenerationError("<Stiff> must be > 0")
    if stiff2 < 0:
        raise GenerationError("<stiff2> must be >= 0")
    sens = max(float(sensitivity), 0.0)
    m_empty = int(empty_levels)
    if m_empty < 0:
        raise GenerationError("<EmptyLevels> must be >= 0")
    v = ModelView(m)
    refs = _validate_pool(m, v, notes)
    off = resolve_offset(m, int(offset))
    wp = water_properties(m.gravity)
    tol = v.tol
    thr = sens if sens > 0 else tol

    # ---- 2. Z-levels --------------------------------------------------------------------------------
    wall_nodes = sorted({n for r in refs for n in r.dof_nodes})
    z = np.array([v.P[n][2] for n in wall_nodes])
    means, lab = cluster_values(z, thr)
    level_of = {n: int(l) for n, l in zip(wall_nodes, lab)}

    # ---- 3. floor template ----------------------------------------------------------------------------
    faces = _horizontal_shells(v, refs) + _free_upward_solid_faces(v, refs)
    if not faces:
        raise GenerationError("no horizontal floor: the pool floor must be made of horizontal SHELL/TSHELL "
                              "elements or of SOLIDs with a free upward face")
    base = lowest_faces(v, faces, thr)
    labs = sorted({level_of[n] for f in base for n in f})
    if len(labs) > 1:
        raise GenerationError(f"the pool floor is not level: its nodes lie on {len(labs)} Z-levels "
                              f"({', '.join(f'{means[k]:g}' for k in labs)}); give <Sensitivity> >= the "
                              f"variation of the floor elevation")
    l0 = labs[0]
    levels = [float(x) for x in means[l0:]]
    n_int = len(levels) - 1
    if n_int < 1:
        raise GenerationError("no wall node above the pool floor: there is no volume to fill")
    if m_empty >= n_int:
        raise GenerationError(f"<EmptyLevels> = {m_empty} leaves no water: the walls have {n_int} level "
                              f"intervals above the floor (levels {', '.join(fmt_num(round(x, 10)) for x in levels)})")
    n_fill = n_int - m_empty
    tg = template_grid(v, base, notes)
    nT = len(tg.xy)
    if not tg.cells:
        raise GenerationError("the floor template has no valid cell")
    bnormals, bedges = _boundary_normals(tg.xy, tg.cells)

    # ---- 4. numbering --------------------------------------------------------------------------------
    gw = (max(m.groups) if m.groups else 0) + 1
    mat_w = (max(m.materials) if m.materials else 0) + 1

    def wnode(k: int, i: int) -> int:
        return off + 1 + k * nT + i

    water_nodes = [wnode(k, i) for k in range(n_fill + 1) for i in range(nT)]
    clash = [n for n in water_nodes if n in m.nodes]
    if clash:                      # cannot happen with offset >= max node; kept as a guard
        raise GenerationError(f"water node numbers {clash[:5]} are already used")

    # ---- 5. interface pairs ---------------------------------------------------------------------------
    from scipy.spatial import cKDTree
    wall_xyz = np.array([v.P[n] for n in wall_nodes])
    tree = cKDTree(wall_xyz)
    rmatch = max(sens, tol) * (1.0 + 1e-9)
    pairs: List[Tuple[int, int]] = []
    springs: List[Tuple[int, int, Tuple[float, float, float]]] = []
    unconnected: List[int] = []
    multiple: List[int] = []
    internal: List[int] = []
    surface: List[int] = []
    oblique = 0
    partner: Dict[int, int] = {}
    for k in range(n_fill + 1):
        for i in range(nT):
            wn = wnode(k, i)
            xyz = np.array([tg.xy[i, 0], tg.xy[i, 1], levels[k]])
            normals = list(bnormals.get(i, []))
            if k == 0:
                normals.append(np.array([0.0, 0.0, -1.0]))
            hits = tree.query_ball_point(xyz, rmatch)
            if not hits:
                if normals:
                    unconnected.append(wn)
                continue
            if not normals:
                (surface if k == n_fill else internal).append(wn)
                continue
            if len(hits) > 1:
                multiple.append(wn)
            d = np.linalg.norm(wall_xyz[hits] - xyz, axis=1)
            j = min(range(len(hits)), key=lambda q: (round(float(d[q]) / max(tol, 1e-300)), wall_nodes[hits[q]]))
            w = wall_nodes[hits[j]]
            kc, exact = interface_constants(normal_projector(normals), stiff, stiff2)
            oblique += int(not exact)
            pairs.append((w, wn))
            partner[wn] = w
            springs.append((w, wn, kc))

    # ---- 6. wetted faces (interface areas) ----------------------------------------------------------
    wet_faces: List[List[int]] = []                    # water-node polygons of the wetted surface
    for c in tg.cells:
        wet_faces.append([wnode(0, i) for i in c])
    for a, b in bedges:
        for k in range(n_fill):
            wet_faces.append([wnode(k, a), wnode(k, b), wnode(k + 1, b), wnode(k + 1, a)])
    shell_polys: List[List[int]] = []
    if int(shell_area) == 1:
        for f in wet_faces:
            if all(n in partner for n in f):
                p = _ordered_unique(partner[n] for n in f)
                if len(p) >= 3:
                    shell_polys.append(p)

    # ---- geometry summaries ---------------------------------------------------------------------------
    depth = levels[n_fill] - levels[0]
    plan_area = float(sum(abs(signed_area(tg.xy[list(c)])) for c in tg.cells))
    perim = float(sum(np.hypot(*(tg.xy[b] - tg.xy[a])) for a, b in bedges))
    volume = plan_area * depth
    mass = volume * wp.rho
    wetted = plan_area + perim * depth

    # ================================= modify the model =================================================
    for k in range(n_fill + 1):
        for i in range(nT):
            m.define_node(wnode(k, i), (float(tg.xy[i, 0]), float(tg.xy[i, 1]), float(levels[k])), 0)
    # water material: type 3 (Vp, Vs) -- D-WAT-01 (G = 1e-8 K)
    m.materials[mat_w] = Material(mat_w, _round_sig(wp.vp), _round_sig(wp.vs), float(wp.weight), wp.damping,
                                  wp.damping, 3)
    grp = Group(gw, SOLID, TITLES["water"])
    eid = 0
    for k in range(n_fill):
        for c in tg.cells:
            eid += 1
            bot = [wnode(k, i) for i in c]
            top = [wnode(k + 1, i) for i in c]
            ns = bot + top if len(c) == 4 else [bot[0], bot[1], bot[2], bot[2], top[0], top[1], top[2], top[2]]
            grp.elements[eid] = Element(eid, ns, mat=mat_w, etype=1)
    m.groups[gw] = grp
    n_elem = eid
    # springs
    gs = gw + 1
    sc_ids: Dict[Tuple[float, float, float], int] = {}
    sc_next = (max(m.springs) if m.springs else 0) + 1
    sgrp = Group(gs, SPRING, TITLES["springs"])
    for e, (w, wn, kc) in enumerate(springs, start=1):
        key = tuple(float(x) for x in kc)
        if key not in sc_ids:
            sc_ids[key] = sc_next
            m.springs[sc_next] = SpringProp(sc_next, key[0], key[1], key[2], 0.0, 0.0, 0.0, 0.0)
            sc_next += 1
        sgrp.elements[e] = Element(e, [w, wn], prop=sc_ids[key])
    if springs:
        m.groups[gs] = sgrp
    else:
        gs = 0
    # fix the rotations the springs introduce at nodes nothing else restrains (EDU-06)
    for _, wn in pairs:
        m.nodes[wn].fix[3:6] = [1, 1, 1]
    wall_rot, wall_drill = fix_interface_rotations(m, v, [w for w, _ in pairs])
    # interface-area shells
    gsh, mat_s, nsh = 0, 0, 0
    if shell_polys:
        gsh = (gs or gw) + 1
        mat_s = mat_w + 1
        edges = []
        for p in shell_polys:
            X = [v.P[n] for n in p]
            edges += [float(np.linalg.norm(X[(i + 1) % len(X)] - X[i])) for i in range(len(X))]
        thick = _round_sig(0.01 * float(np.mean(edges)), 6) if edges else 0.01
        m.materials[mat_s] = Material(mat_s, _round_sig(SHELL_AREA_E_OVER_K * wp.K), SHELL_AREA_NU, 0.0, 0.0, 0.0, 1)
        shg = Group(gsh, SHELL, TITLES["shells"])
        for e, p in enumerate(shell_polys, start=1):
            shg.elements[e] = Element(e, list(p), mat=mat_s, etype=1, thick=thick)
        m.groups[gsh] = shg
        nsh = len(shell_polys)
    elif int(shell_area) == 1:
        notes.warn("<ShellArea> = 1 but no wetted face has a wall node at every corner: no interface-area shells")
    # incompatible modes for the nearly incompressible water (module docstring, limit 3)
    mo = m.options.ensure_record("MOPT")
    changed = int(mo.get("incomp")) != 0
    if changed:
        mo.set("incomp", 0)
    # pool data
    rec = PoolDataRecord("POOLDATA", [gw, gs, gsh, fmt_num(stiff), fmt_num(sens), m_empty, int(shell_area), off,
                                      fmt_num(stiff2), mat_w, fmt_num(levels[0]), fmt_num(levels[n_fill])])
    m.options.set_record(rec)
    check_node_limit(m, notes)

    # ---- diagnostics --------------------------------------------------------------------------------
    kx = sum(kc[0] for _, _, kc in springs)
    ky = sum(kc[1] for _, _, kc in springs)
    kz = sum(kc[2] for _, _, kc in springs)
    freq = {}
    for lab_, kk in (("x", kx), ("y", ky), ("z", kz)):
        freq[lab_] = math.sqrt(kk / mass) / (2.0 * math.pi) if mass > 0 and kk > 0 else 0.0
    if unconnected:
        notes.warn(f"{len(unconnected)} water nodes on the wetted boundary have no coincident wall/floor node and are "
                   f"not attached (non-conforming wall mesh?): {_short(unconnected)}")
    if multiple:
        notes.warn(f"{len(multiple)} water nodes coincide with more than one wall/floor node (unwelded junction?); "
                   f"each is attached to the nearest (lowest-numbered) one -- WELD the pool model first: "
                   f"{_short(multiple)}")
    if internal:
        notes.warn(f"{len(internal)} water nodes inside the water coincide with wall nodes (internal walls or "
                   f"columns are not modelled: fill each compartment separately); not attached: {_short(internal)}")
    if surface:
        notes.warn(f"{len(surface)} nodes of the water surface coincide with wall/roof nodes and are not attached: "
                   f"FILLPOOL leaves the surface free (a roof touching the water is not modelled; use "
                   f"<EmptyLevels> >= 1 for an open pool): {_short(surface)}")
    if oblique:
        notes.warn(f"{oblique} interface springs on oblique walls: uncoupled SC constants = diagonal of "
                   f"Stiff P_N + stiff2 (I - P_N) (OQ-2 (b)); the water cannot slip freely along those walls")
    if not wp.recognised:
        notes.warn(f"the gravity {fmt_num(m.gravity)} is neither ~9.81 (kN-m) nor ~32.2 (kip-ft): the water material "
                   f"{mat_w} assumes {wp.label} units -- redefine it with M,{mat_w},<Vp>,<Vs>,<weight>,...,3")
    slosh: Dict[str, Tuple[float, float, float]] = {}
    for lab_, col in (("x", 0), ("y", 1)):
        half = 0.5 * float(np.ptp(tg.xy[:, col]))
        if half > 0 and depth > 0:
            slosh[lab_] = (half, sloshing_frequency_exact(half, depth, m.gravity),
                           housner_convective(half, depth, m.gravity)[0])
    res = FillPoolResult(gw, gs, gsh, mat_w, mat_s, off, levels, n_fill, levels[0], levels[n_fill],
                         water_nodes, n_elem, nT, len(tg.cells), pairs, {i: k for k, i in sc_ids.items()}, oblique,
                         unconnected, multiple, internal, surface, nsh, wall_rot, wall_drill, volume, mass, wetted,
                         wp, freq, changed, slosh)
    return res


# ======================================================================================
# Pool interface (LISTPOOLINTER, MERGEPOOL)
# ======================================================================================
def pool_interface(pool: SSIModel) -> List[Tuple[int, int]]:
    """(wall node, water node) of every FILLPOOL interface spring of a pool model (from its POOLDATA record);
    :class:`GenerationError` when the model was not filled by FILLPOOL or its groups were deleted."""
    rec = pool_record(pool)
    if rec is None:
        raise GenerationError("the pool model was not filled with FILLPOOL (no pool data)")
    gs = int(rec.springs)
    if gs == 0:
        return []
    g = pool.groups.get(gs)
    if g is None or g.type != SPRING:
        raise GenerationError(f"the FILLPOOL interface springs (group {gs}) are no longer in the pool model")
    return [(int(e.nodes[0]), int(e.nodes[1])) for e in g.sorted_elements() if len(e.nodes) >= 2]


@dataclass
class InterfaceMatch:
    matched: List[Tuple[int, np.ndarray]]          # (node, global xyz)
    missing: List[int]                              # interface node numbers not in the original model
    moved: List[Tuple[int, float]]                  # (node, distance) same number, other position


def match_interface(original: SSIModel, pool: SSIModel) -> InterfaceMatch:
    """LISTPOOLINTER: compare the pool's interface (wall) nodes with the nodes of the original model: the
    ones with the same number **and** the same position (distance <= the larger geometric tolerance of the
    two models) are listed; the others are returned as missing / moved."""
    pairs = pool_interface(pool)
    walls = sorted({w for w, _ in pairs})
    vo, vp = ModelView(original), ModelView(pool)
    tol = max(vo.tol, vp.tol)
    res = InterfaceMatch([], [], [])
    for n in walls:
        if n not in vp.P:
            continue
        if n not in vo.P:
            res.missing.append(n)
            continue
        d = float(np.linalg.norm(vo.P[n] - vp.P[n]))
        if d <= tol:
            res.matched.append((n, vo.P[n]))
        else:
            res.moved.append((n, d))
    return res


@dataclass
class MergePoolResult:
    groups: Dict[int, int]            # pool group -> new group
    materials: Dict[int, int]
    springs: Dict[int, int]
    nodes: List[int]                  # water nodes added
    wall_rotations_fixed: List[int]
    wall_drilling_fixed: List[int]
    mopt_changed: bool


def merge_pool(dest: SSIModel, pool: SSIModel, notes: Optional[Notes] = None) -> MergePoolResult:
    """MERGEPOOL (SASSI-EDU): import the FILLPOOL groups (water SOLIDs, interface SPRINGs, interface-area
    SHELLs) of ``pool`` into ``dest`` **keeping the node numbers** -- the wall nodes are the original model's
    nodes and the water nodes were numbered from the FILLPOOL offset + 1.

    Requires every interface wall node in ``dest`` with the same position (LISTPOOLINTER) and no water node
    number in ``dest``.  Groups are appended after the last group of ``dest``, the water / shell materials
    and the SC properties get numbers after the maxima of ``dest``; element numbers are kept.  The water
    nodes keep their fixities (rotations fixed); interface wall nodes without a rotational DOF from the
    elements of ``dest`` get their rotations fixed; MOPT ``<incomp>`` is set to 0."""
    notes = notes if notes is not None else Notes()
    rec = pool_record(pool)
    if rec is None:
        raise GenerationError("the pool model was not filled with FILLPOOL (no pool data)")
    gids = [int(g) for g in (rec.water, rec.springs, rec.shells) if int(g)]
    missing_g = [g for g in gids if g not in pool.groups]
    if missing_g:
        raise GenerationError(f"FILLPOOL groups {missing_g} are no longer in the pool model")
    match = match_interface(dest, pool)
    if match.missing or match.moved:
        bad = [n for n in match.missing] + [n for n, _ in match.moved]
        raise GenerationError(f"{len(bad)} pool interface nodes are not in the active model with the same number and "
                              f"position (LISTPOOLINTER): {_short(sorted(bad))}")
    walls = {w for w, _ in pool_interface(pool)}
    water_nodes = sorted({n for g in gids for e in pool.groups[g].elements.values() for n in e.nodes if n}
                         - walls)
    clash = [n for n in water_nodes if n in dest.nodes]
    if clash:
        raise GenerationError(f"{len(clash)} water node numbers are already used in the active model "
                              f"({_short(clash)}): fill the pool with <offset> >= {max(dest.nodes)} (the last node "
                              f"of the original model), or the water was already imported")
    vd = ModelView(dest)
    vpool = ModelView(pool)
    # nodes (global coordinates of the pool model)
    for n in water_nodes:
        p = vpool.P[n]
        dest.define_node(n, (float(p[0]), float(p[1]), float(p[2])), 0)
        dest.nodes[n].fix = list(pool.nodes[n].fix)
        dest.nodes[n].flags = set(pool.nodes[n].flags)
    # tables
    mmap: Dict[int, int] = {}
    smap: Dict[int, int] = {}
    mnext = (max(dest.materials) if dest.materials else 0) + 1
    snext = (max(dest.springs) if dest.springs else 0) + 1
    gnext = (max(dest.groups) if dest.groups else 0) + 1
    gmap: Dict[int, int] = {}
    for g in gids:
        src = pool.groups[g]
        ng = Group(gnext, src.type, src.title)
        gmap[g] = gnext
        gnext += 1
        for e in src.sorted_elements():
            ne = e.copy()
            if src.type in (SOLID, SHELL, TSHELL):
                if e.mat not in mmap:
                    pm = pool.materials.get(e.mat)
                    if pm is None:
                        raise GenerationError(f"material {e.mat} of pool group {g} is not defined")
                    mmap[e.mat] = mnext
                    mnext += 1
                ne.mat = mmap[e.mat]
            elif src.type == SPRING:
                if e.prop not in smap:
                    ps = pool.springs.get(e.prop)
                    if ps is None:
                        raise GenerationError(f"spring property {e.prop} of pool group {g} is not defined")
                    smap[e.prop] = snext
                    snext += 1
                ne.prop = smap[e.prop]
            ng.elements[ne.id] = ne
        dest.groups[ng.id] = ng
    for old, new in mmap.items():
        pm = pool.materials[old]
        dest.materials[new] = Material(new, pm.val1, pm.val2, pm.weight, pm.pdamp, pm.sdamp, pm.mtype)
    for old, new in smap.items():
        ps = pool.springs[old]
        dest.springs[new] = SpringProp(new, ps.scx, ps.scy, ps.scz, ps.scxx, ps.scyy, ps.sczz, ps.damp)
    # rotations of interface wall nodes that nothing in dest restrains (SOLID walls, flat SHELL walls)
    fixed, drill = fix_interface_rotations(dest, vd, walls)
    mo = dest.options.ensure_record("MOPT")
    changed = int(mo.get("incomp")) != 0
    if changed:
        mo.set("incomp", 0)
    check_node_limit(dest, notes)
    return MergePoolResult(gmap, mmap, smap, water_nodes, fixed, drill, changed)


# ======================================================================================
# REFINEMODEL (manual 9.16.2; spec 11 section 1.2)
# ======================================================================================
#: natural coordinates (0/2 lattice indices) of the 8 SOLID nodes (spec 08 section 4.4)
_HEX_LATTICE = ((0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0), (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2))
_QUAD_TYPES = (SHELL, TSHELL, PLANE)


@dataclass
class RefineResult:
    quads: int                              # parents split into 4
    hexes: int                              # parents split into 8
    new_nodes: List[int]
    groups: Dict[int, Dict[int, List[int]]]  # refined group -> {old element: [new elements]}
    untouched: int                          # triangles / prisms / other elements left unchanged
    hanging: int                            # refined edges shared with unrefined 2D/3D elements
    interaction_added: int


def _is_quad(t: int, ns: Sequence[int]) -> bool:
    return t in _QUAD_TYPES and len(ns) == 4 and all(ns) and len(set(ns)) == 4


def _is_hex(t: int, ns: Sequence[int]) -> bool:
    return t == SOLID and len(ns) == 8 and all(ns) and len(set(ns)) == 8


def refine_model(m: SSIModel, notes: Optional[Notes] = None) -> RefineResult:
    """REFINEMODEL: split every quadrilateral SHELL/TSHELL/PLANE into 4 and every hexahedral SOLID into 8.

    New nodes (global coordinates, numbered after the last node in order of creation): edge midpoints
    (shared through the sorted node pair), face centres = average of the 4 face corners (shared through the
    sorted node quadruple; also the centre of a split quad) and hexahedron body centres = average of the 8
    corners -- the bilinear / trilinear centres, so every child is an exact sub-patch of its parent and the
    area / volume is preserved.  Quad children ``(n1,m12,c,m41) (m12,n2,m23,c) (c,m23,n3,m34) (m41,c,m34,n4)``
    and the 2x2x2 hexahedron children keep the parent's node order (orientation, positive Jacobian).
    Children inherit material, property, ETYPE, EINT, THICK and the releases."""
    notes = notes if notes is not None else Notes()
    v = ModelView(m)
    todo = []
    untouched = 0
    for g, e in m.iter_elements():
        ns = list(e.nodes)
        if _is_quad(g.type, ns) or _is_hex(g.type, ns):
            if not v.defined(ns):
                raise GenerationError(f"element {e.id} of group {g.id} references undefined nodes (CHECK Error 41)")
            todo.append((g.id, e.id))
        else:
            untouched += 1
    if not todo:
        raise GenerationError("no quadrilateral SHELL/TSHELL/PLANE or hexahedral SOLID element to refine")
    next_id = (max(m.nodes) if m.nodes else 0) + 1
    new_nodes: List[int] = []
    parents: Dict[int, Tuple[int, ...]] = {}
    xyz_new: Dict[int, np.ndarray] = {}
    keyed: Dict[Tuple[int, ...], int] = {}

    def node_for(corners: Sequence[int], share: bool = True) -> int:
        nonlocal next_id
        key = tuple(sorted(corners))
        if share and key in keyed:
            return keyed[key]
        nid = next_id
        next_id += 1
        xyz_new[nid] = np.mean([v.P[c] for c in corners], axis=0)
        parents[nid] = key
        new_nodes.append(nid)
        if share:
            keyed[key] = nid
        return nid

    children: Dict[Tuple[int, int], List[List[int]]] = {}
    nq = nh = 0
    for gid, eid in todo:
        g = m.groups[gid]
        ns = list(g.elements[eid].nodes)
        if g.type in _QUAD_TYPES:
            n1, n2, n3, n4 = ns
            m12, m23, m34, m41 = node_for((n1, n2)), node_for((n2, n3)), node_for((n3, n4)), node_for((n4, n1))
            c = node_for((n1, n2, n3, n4))
            children[(gid, eid)] = [[n1, m12, c, m41], [m12, n2, m23, c], [c, m23, n3, m34], [m41, c, m34, n4]]
            nq += 1
        else:
            corner = {lat: ns[k] for k, lat in enumerate(_HEX_LATTICE)}
            lattice: Dict[Tuple[int, int, int], int] = dict(corner)
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        lat = (i, j, k)
                        if lat in lattice:
                            continue
                        cs = [corner[(a, b, c)] for a in ((i,) if i != 1 else (0, 2)) for b in ((j,) if j != 1 else (0, 2))
                              for c in ((k,) if k != 1 else (0, 2))]
                        lattice[lat] = node_for(cs, share=len(cs) < 8)
            kids = []
            for kk in (0, 1):
                for jj in (0, 1):
                    for ii in (0, 1):
                        kids.append([lattice[(ii + a // 2, jj + b // 2, kk + c // 2)] for a, b, c in _HEX_LATTICE])
            children[(gid, eid)] = kids
            nh += 1
    # hanging nodes: refined edges used by unrefined triangles/prisms/quads of other types
    refined_edges = {k for k in keyed if len(k) == 2}
    hanging = 0
    for g, e in m.iter_elements():
        if (g.id, e.id) in children or g.type in (BEAMS, SPRING, GENERAL):
            continue
        ns = [n for n in e.nodes if n]
        if any(tuple(sorted((a, b))) in refined_edges for a in ns for b in ns if a < b):
            hanging += 1
    # ================================= modify the model =================================================
    added_int = 0
    for nid in new_nodes:
        p = xyz_new[nid]
        m.define_node(nid, (float(p[0]), float(p[1]), float(p[2])), 0)
        par = [m.nodes[c] for c in parents[nid]]
        flags = set.intersection(*(set(nd.flags) for nd in par))
        m.nodes[nid].flags = flags
        added_int += int(0 in flags)
        m.nodes[nid].fix = [int(all(nd.fix[d] for nd in par)) for d in range(6)]
    maps: Dict[int, Dict[int, List[int]]] = {}
    for gid in sorted({g for g, _ in todo}):
        g = m.groups[gid]
        old = g.sorted_elements()
        new_elems: Dict[int, Element] = {}
        mp: Dict[int, List[int]] = {}
        k = 0
        for e in old:
            kids = children.get((gid, e.id))
            if kids is None:
                k += 1
                new_elems[k] = e.copy(new_id=k)
                mp[e.id] = [k]
                continue
            mp[e.id] = []
            for nodes in kids:
                k += 1
                new_elems[k] = e.copy(new_id=k, nodes=nodes)
                mp[e.id].append(k)
        g.elements = new_elems
        maps[gid] = mp
    # output requests follow the new element numbers
    for r in m.eout:
        mp = maps.get(r.group)
        if mp:
            r.elements = sorted({x for e in r.elements for x in mp.get(e, [e])})
    check_node_limit(m, notes)
    if hanging:
        notes.warn(f"{hanging} unrefined elements (triangles, prisms, ...) share a refined edge: the mesh has hanging "
                   f"nodes there (non-conforming)")
    loaded = sorted(n for n in set(m.forces) | set(m.moments) | set(m.tmass) | set(m.rmass)
                    if any(n in nodes for kids in children.values() for nodes in kids))
    if loaded:
        notes.warn(f"{len(loaded)} nodes of refined elements carry loads or masses; they are not redistributed to "
                   f"the new nodes")
    if added_int:
        depths = v.interfaces()
        if depths is not None:
            t = v.interface_tol(depths)
            off = sorted({round(float(xyz_new[n][2]), 9) for n in new_nodes if 0 in m.nodes[n].flags
                          and v.interface_of(float(xyz_new[n][2]), depths)[1] > t})
            if off:
                notes.warn(f"new interaction nodes lie at elevations {', '.join(f'{z:g}' for z in off[:6])} that are "
                           f"not soil-layer interfaces (EDU-01): refine the TOPL layers too")
    if pool_record(m) is not None:
        notes.warn("the model holds FILLPOOL groups: the interface springs are not regenerated for the new nodes -- "
                   "refine the pool walls before FILLPOOL")
    return RefineResult(nq, nh, new_nodes, maps, untouched, hanging, added_int)
