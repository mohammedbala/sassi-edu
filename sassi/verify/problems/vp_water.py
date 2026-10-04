"""Verification problems of the water modelling commands (ACS SASSI manual section 9.16; requirements
section 3.4.N; spec 11 sections 1 and 7; decisions D-WAT-01, D-CHK-07): VP-54W and VP-W1.

* VP-54W  water and refinement part of VP-54 (requirements 6.3; spec 11 section 7 items 1-2):
          REFINEMODEL -- one quadrilateral -> 4 with 5 new nodes, a 2 x 2 mesh -> 16 quadrilaterals and
          25 nodes (shared midpoints), a (distorted) hexahedron -> 8 with 19 new nodes, area and volume
          preserved to round-off, positive Jacobians, triangles and prisms untouched; FILLPOOL -- a
          rectangular SHELL pool a x b with walls H high, levels every h and ``EmptyLevels`` = m holds
          water SOLIDs of volume a b (H - m h) (computed independently from the element geometry), one
          spring per coincident wall/water node pair (counted independently), ``Stiff`` along the wall
          normals and ``stiff2`` along the walls, an error for 0 < offset < the largest node number (model
          unchanged), the water mass rho V with rho from the water material, the same for a pool with
          SOLID floor and walls; LISTPOOLINTER finds every interface node of a CUT2SUB pool in the
          original model; MERGEPOOL imports the water with the node numbers kept; WRITE -> INP of a
          filled pool gives the same model (UT-03).
* VP-W1   hydrodynamic (impulsive) mass of a rigid rectangular tank, length 2L = 10 m in the direction
          of the motion, width 2 m, walls 6 m high, filled by FILLPOOL to H = 5 m (EmptyLevels = 2, 0.5 m
          levels), on a practically rigid site (Vs = 10 km/s), run SITE -> POINT -> HOUSE -> ANALYS from
          the command language, in two variants: SHELL walls and floor at grade, and SOLID walls on a SOLID
          slab with the interface-area shells (``ShellArea`` = 1; the water on the slab top, the template
          from the free upward slab faces); harmonic horizontal (SV, X) input at 2, 4 and 8 Hz -- far below the
          compression frequency of the water column (c/4H = 74 Hz) and far above the spurious shear modes
          of the water solids (< 0.05 Hz).  Exact properties (``check``): the water mass in the HOUSE
          mass matrix equals rho a b H (1e-10); the base shear (sum of the X forces of the interface
          springs, FILE4 recovery operators) equals the inertia of the water (dynamic equilibrium, 1e-6);
          the tank follows the control motion (|H - 1| <= 1e-3).  Engineering criterion (``check``):
          F/(m a) equals the impulsive mass ratio of the **exact potential-flow solution** for a rigid
          rectangular tank, m_i/m = 2/(L H^2) sum tanh(l_n L)/l_n^3 = 0.5 for L = H, within 5 % (the FE
          discretisation error of the 0.5 m mesh is about +1.5 %, see docs/user/WATER.md).  Informative
          (``inform``): Housner's (1963) closed form tanh(sqrt 3 L/H)/(sqrt 3 L/H) = 0.5423, Westergaard's
          (1933) long-reservoir estimate (7/12) H/L = 0.5833, and the height of the resultant of the wall
          pressures (EBP) against the exact series (0.4047 H) and Housner (3/8 H).

Water model (D-WAT-01 as amended in :mod:`sassi.prep.water_lib`): SOLIDs with K = 2.2 GPa, G = 1e-8 K
(Vs = 1e-4 Vp), unit weight 9.81 kN/m3, 0.5 % damping, incompatible modes (MOPT,0); interface springs
along the wall normals only (frictionless).  Sloshing (the convective mass) is not modelled, so the
reference is the impulsive mass.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .. import VPResult, problem

# ======================================================================================
# Model builders (command text; also used by tests/unit/test_water*.py)
# ======================================================================================
WALL_E = 3.0e10            # kPa: 1000 x concrete -> practically rigid walls


class _Nodes:
    """Node numbering by coordinates (one N line per new point)."""

    def __init__(self, first: int = 1):
        self.ids: Dict[Tuple[float, float, float], int] = {}
        self.lines: List[str] = []
        self.next = first

    def __call__(self, x: float, y: float, z: float) -> int:
        key = (round(float(x), 9), round(float(y), 9), round(float(z), 9))
        if key not in self.ids:
            self.ids[key] = self.next
            self.lines.append(f"N,{self.next},{key[0]:.12g},{key[1]:.12g},{key[2]:.12g}")
            self.next += 1
        return self.ids[key]


def shell_pool_lines(a: float, b: float, H: float, nx: int, ny: int, nz: int, x0: float = 0.0, y0: float = 0.0,
                     z0: float = 0.0, group: int = 1, mat: int = 1, thick: float = 0.5, E: float = WALL_E,
                     weight: float = 0.0, first_node: int = 1, define_material: bool = True) -> List[str]:
    """Command lines of an open-top rectangular SHELL pool: floor ``[x0, x0+a] x [y0, y0+b]`` at ``z0`` and four
    vertical walls of height ``H``; ``nx x ny`` floor elements, ``nz`` wall elements over the height (levels
    every ``H/nz``).  The floor nodes are numbered first (row by row, x fastest), then the wall nodes level by
    level.  Material ``mat`` = (E, nu 0.2, ``weight``, no damping), thickness ``thick``."""
    xs = x0 + np.linspace(0.0, a, nx + 1)
    ys = y0 + np.linspace(0.0, b, ny + 1)
    zs = z0 + np.linspace(0.0, H, nz + 1)
    nid = _Nodes(first_node)
    quads: List[Tuple[int, int, int, int]] = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            nid(xs[i], ys[j], z0)
    for j in range(ny):
        for i in range(nx):
            quads.append((nid(xs[i], ys[j], z0), nid(xs[i + 1], ys[j], z0), nid(xs[i + 1], ys[j + 1], z0),
                          nid(xs[i], ys[j + 1], z0)))
    for k in range(nz):
        for i in range(nx):
            for y in (ys[0], ys[-1]):
                quads.append((nid(xs[i], y, zs[k]), nid(xs[i + 1], y, zs[k]), nid(xs[i + 1], y, zs[k + 1]),
                              nid(xs[i], y, zs[k + 1])))
        for j in range(ny):
            for x in (xs[0], xs[-1]):
                quads.append((nid(x, ys[j], zs[k]), nid(x, ys[j + 1], zs[k]), nid(x, ys[j + 1], zs[k + 1]),
                              nid(x, ys[j], zs[k + 1])))
    lines = list(nid.lines)
    if define_material:
        lines.append(f"M,{mat},{E:.12g},0.2,{weight:.12g},0,0,1")
    lines += [f"GROUP,{group},SHELL", f"GTIT,{group},pool walls and floor", f"MACT,{mat}"]
    lines += [f"E,{e},{q[0]},{q[1]},{q[2]},{q[3]}" for e, q in enumerate(quads, start=1)]
    lines.append(f"THICK,1,{len(quads)},1,{thick:.12g}")
    return lines


def solid_pool_lines(a: float, b: float, H: float, nx: int, ny: int, nz: int, tw: float = 0.5, tf: float = 0.5,
                     group: int = 1, mat: int = 1, E: float = WALL_E, weight: float = 0.0) -> List[str]:
    """Command lines of a rectangular pool of SOLID elements: inner plan ``[0, a] x [0, b]`` (``nx x ny``
    cells), a floor slab ``tf`` thick under the inner plan and the walls, walls ``tw`` thick and ``H`` high
    (``nz`` element layers) sitting on the slab -- the free upward slab faces inside the walls are the floor
    template of FILLPOOL.  The slab-bottom nodes are numbered first, row by row (x fastest)."""
    xs = np.concatenate([[-tw], np.linspace(0.0, a, nx + 1), [a + tw]])
    ys = np.concatenate([[-tw], np.linspace(0.0, b, ny + 1), [b + tw]])
    zs = np.concatenate([[-tf], np.linspace(0.0, H, nz + 1)])
    nid = _Nodes(1)
    for y in ys:                                 # the slab bottom first: nodes 1 .. (nx + 3)(ny + 3)
        for x in xs:
            nid(x, y, zs[0])
    hexes = []
    nX, nY = len(xs) - 1, len(ys) - 1
    for k in range(len(zs) - 1):
        for j in range(nY):
            for i in range(nX):
                inner = 0 < i < nX - 1 and 0 < j < nY - 1
                if k > 0 and inner:
                    continue                     # the pool
                c = [(xs[i], ys[j]), (xs[i + 1], ys[j]), (xs[i + 1], ys[j + 1]), (xs[i], ys[j + 1])]
                hexes.append([nid(x, y, zs[k]) for x, y in c] + [nid(x, y, zs[k + 1]) for x, y in c])
    lines = list(nid.lines)
    lines.append(f"M,{mat},{E:.12g},0.2,{weight:.12g},0,0,1")
    lines += [f"GROUP,{group},SOLID", f"GTIT,{group},pool walls and floor slab", f"MACT,{mat}"]
    lines += ["E," + ",".join(str(x) for x in [e] + h) for e, h in enumerate(hexes, start=1)]
    lines.append(f"ETYPE,1,{len(hexes)},1,1")
    return lines


# ======================================================================================
# Independent geometry helpers
# ======================================================================================
def hex_volume(X: np.ndarray) -> float:
    """Volume of a trilinear hexahedron (8, 3) by 2x2x2 Gauss quadrature of det J (exact for trilinear maps)."""
    from sassi.elements.solid import gauss_3d, hex8_shape
    pts, w = gauss_3d(2)
    _, dN = hex8_shape(pts)
    J = np.einsum("gai,aj->gij", dN, np.asarray(X, float))
    return float(np.sum(w * np.linalg.det(J)))


def min_hex_jacobian(X: np.ndarray) -> float:
    from sassi.elements.solid import gauss_3d, hex8_shape
    pts, _ = gauss_3d(2)
    _, dN = hex8_shape(pts)
    return float(np.min(np.linalg.det(np.einsum("gai,aj->gij", dN, np.asarray(X, float)))))


def quad_area(X: np.ndarray, n: int = 2) -> Tuple[float, np.ndarray]:
    """Area and mean unit normal of a bilinear quadrilateral (4, 3) by n x n Gauss quadrature of
    |dx/dxi x dx/deta| (2 points are exact for a planar quadrilateral; a warped one needs more)."""
    g, w = np.polynomial.legendre.leggauss(n)
    X = np.asarray(X, float)
    A = 0.0
    nsum = np.zeros(3)
    for xi, wi in zip(g, w):
        for eta, wj in zip(g, w):
            dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
            deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
            c = np.cross(dxi @ X, deta @ X)
            A += wi * wj * float(np.linalg.norm(c))
            nsum += wi * wj * c
    return A, nsum / np.linalg.norm(nsum)


def element_xyz(m, group: int, eid: int) -> np.ndarray:
    _, xyz = m.global_coordinates(m.groups[group].elements[eid].nodes)
    return xyz


# ======================================================================================
# VP-54W: REFINEMODEL and FILLPOOL properties
# ======================================================================================
def _ui(workdir):
    from sassi.prep import Interpreter
    return Interpreter(cwd=workdir)


def _errors(ui):
    from sassi.prep import Kind
    return ui.sink.texts(Kind.ERROR)


def _refine_checks(r: VPResult, workdir: Path) -> None:
    # (1) one quadrilateral (warped, rotated) -> 4, 5 new nodes, area preserved, same orientation
    ui = _ui(workdir)
    ui.run_text("\n".join(["N,1,0,0,0", "N,2,2.0,0.3,0.1", "N,3,2.4,1.9,0.3", "N,4,-0.2,1.5,0.05",
                           "M,1,3e7,0.2,24,0.05,0.05", "GROUP,1,SHELL", "E,1,1,2,3,4", "THICK,1,1,1,0.3"]))
    m = ui.model
    A0, n0 = quad_area(element_xyz(m, 1, 1), n=24)
    ui.execute("REFINEMODEL")
    r.require("REFINEMODEL on one quadrilateral runs without error", not _errors(ui), "; ".join(_errors(ui)))
    g = m.groups[1]
    r.check("one quadrilateral: number of SHELLs after REFINEMODEL", len(g.elements), 4, atol=0)
    r.check("one quadrilateral: new nodes (4 edge midpoints + 1 centre)", len(m.nodes) - 4, 5, atol=0)
    areas = [quad_area(element_xyz(m, 1, e), n=24) for e in sorted(g.elements)]
    r.check("one warped quadrilateral: sum of the child areas = parent area (bilinear surface, 24x24 Gauss)",
            sum(a for a, _ in areas), A0, rtol=1e-12)
    r.require("one quadrilateral: every child keeps the parent orientation (normal . parent normal > 0)",
              all(float(nn @ n0) > 0 for _, nn in areas))
    r.require("one quadrilateral: the children keep the material, thickness and ETYPE",
              all(e.mat == 1 and e.thick == 0.3 and e.etype == 0 for e in g.elements.values()))
    # (2) 2 x 2 mesh -> 16 quadrilaterals, 25 nodes (shared midpoints)
    ui = _ui(workdir)
    lines = [f"N,{1 + i + 3 * j},{i:.1f},{j:.1f},0" for j in range(3) for i in range(3)]
    lines += ["M,1,3e7,0.2,24,0.05,0.05", "GROUP,1,SHELL"]
    lines += [f"E,{1 + i + 2 * j},{1 + i + 3 * j},{2 + i + 3 * j},{5 + i + 3 * j},{4 + i + 3 * j}"
              for j in range(2) for i in range(2)]
    lines += ["THICK,1,4,1,0.2"]
    ui.run_text("\n".join(lines))
    ui.execute("REFINEMODEL")
    m = ui.model
    r.check("2 x 2 quadrilateral mesh: SHELLs after REFINEMODEL", len(m.groups[1].elements), 16, atol=0)
    r.check("2 x 2 quadrilateral mesh: nodes after REFINEMODEL (shared midpoints)", len(m.nodes), 25, atol=0)
    r.check("2 x 2 quadrilateral mesh: total area", sum(quad_area(element_xyz(m, 1, e))[0]
                                                         for e in m.groups[1].elements), 4.0, rtol=1e-12)
    r.require("2 x 2 quadrilateral mesh: element numbers 1..16 without gaps",
              sorted(m.groups[1].elements) == list(range(1, 17)))
    # (3) one distorted hexahedron -> 8, 19 new nodes, volume preserved, positive Jacobians
    ui = _ui(workdir)
    P = [(0, 0, 0), (2.0, 0.1, -0.1), (2.2, 1.8, 0.0), (-0.1, 1.6, 0.1), (0.1, -0.1, 1.5), (1.9, 0.0, 1.7),
         (2.3, 2.0, 1.6), (0.0, 1.7, 1.4)]
    ui.run_text("\n".join([f"N,{k + 1},{x},{y},{z}" for k, (x, y, z) in enumerate(P)]
                          + ["M,1,3e7,0.2,24,0.05,0.05", "GROUP,1,SOLID", "E,1,1,2,3,4,5,6,7,8", "ETYPE,1,1,1,1"]))
    m = ui.model
    V0 = hex_volume(element_xyz(m, 1, 1))
    ui.execute("REFINEMODEL")
    g = m.groups[1]
    r.check("one hexahedron: SOLIDs after REFINEMODEL", len(g.elements), 8, atol=0)
    r.check("one hexahedron: new nodes (12 edge + 6 face + 1 body)", len(m.nodes) - 8, 19, atol=0)
    vols = [hex_volume(element_xyz(m, 1, e)) for e in sorted(g.elements)]
    r.check("one hexahedron: sum of the child volumes = parent volume (trilinear map)", sum(vols), V0, rtol=1e-12)
    r.require("one hexahedron: every child has a positive Jacobian (node order kept)",
              min(min_hex_jacobian(element_xyz(m, 1, e)) for e in g.elements) > 0)
    r.require("one hexahedron: the children keep ETYPE 1", all(e.etype == 1 for e in g.elements.values()))
    # (4) triangles and prisms are untouched; a hexahedron beside a prism is split (hanging nodes reported)
    ui = _ui(workdir)
    pts = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
           (2, 0, 0), (2, 0, 1), (5, 0, 0), (6, 0, 0), (5, 1, 0)]
    ui.run_text("\n".join([f"N,{k + 1},{x},{y},{z}" for k, (x, y, z) in enumerate(pts)]
                          + ["M,1,3e7,0.2,24,0.05,0.05", "GROUP,1,SOLID", "E,1,1,2,3,4,5,6,7,8",
                             "E,2,2,9,3,3,6,10,7,7", "ETYPE,1,2,1,1", "GROUP,2,SHELL", "E,1,11,12,13,13",
                             "THICK,1,1,1,0.1"]))
    m = ui.model
    prism = list(m.groups[1].elements[2].nodes)
    tri = list(m.groups[2].elements[1].nodes)
    ui.sink.clear()
    ui.execute("REFINEMODEL")
    from sassi.prep import Kind
    warns = ui.sink.texts(Kind.WARNING)
    g1, g2 = m.groups[1], m.groups[2]
    r.check("hexahedron + prism: SOLIDs after REFINEMODEL (8 children + the prism)", len(g1.elements), 9, atol=0)
    r.require("the prism is unchanged (renumbered after the 8 children)", list(g1.elements[9].nodes) == prism)
    r.require("the triangular SHELL (repeated node) is unchanged", list(g2.elements[1].nodes) == tri)
    r.require("hanging nodes next to the unrefined prism are reported", any("hanging" in w for w in warns),
              "; ".join(warns))


def _fill_checks(r: VPResult, workdir: Path) -> None:
    from sassi.prep import Kind
    from sassi.prep import water_lib as wl
    from sassi.prep.writer import write_pre
    a, b, H, h, m_empty = 4.0, 3.0, 3.0, 0.5, 2
    stiff, stiff2 = 2.0e8, 3.0e3
    ui = _ui(workdir)
    ui.run_text("\n".join(["GRAVITY,9.81"] + shell_pool_lines(a, b, H, 4, 3, int(round(H / h)))))
    m = ui.model
    nmax = max(m.nodes)
    before = m.canonical()
    # offset 0 < offset < Nmax is an error and leaves the model unchanged
    ui.sink.clear()
    ui.execute(f"FILLPOOL,{stiff},0,{m_empty},1,{nmax - 1},{stiff2}")
    r.require("FILLPOOL with 0 < offset < the largest node number is an error",
              any("offset" in e for e in ui.sink.texts(Kind.ERROR)), "; ".join(ui.sink.texts(Kind.ERROR)))
    r.require("... and the model is unchanged", m.canonical() == before)
    ui.sink.clear()
    ui.execute(f"FILLPOOL,{stiff},0,{m_empty},1,{nmax + 100},{stiff2}")
    r.require("FILLPOOL on the SHELL pool runs without error", not _errors(ui), "; ".join(_errors(ui)))
    rec = wl.pool_record(m)
    r.require("FILLPOOL stores the pool data record", rec is not None)
    if rec is None:
        return
    gw, gs, gsh = int(rec.water), int(rec.springs), int(rec.shells)
    water = m.groups[gw]
    r.require("water nodes are numbered from offset + 1", min(n for e in water.elements.values() for n in e.nodes)
              == nmax + 101)
    vol = sum(hex_volume(element_xyz(m, gw, e)) for e in water.elements)
    r.check(f"water volume = a b (H - m h) = {a}*{b}*({H} - {m_empty}*{h})", vol, a * b * (H - m_empty * h),
            rtol=1e-12)
    # springs = coincident wall/water pairs, counted independently
    wall_nodes = sorted({n for e in m.groups[1].elements.values() for n in e.nodes})
    water_nodes = sorted({n for e in water.elements.values() for n in e.nodes})
    _, Pw = m.global_coordinates(wall_nodes)
    _, Pq = m.global_coordinates(water_nodes)
    d = np.linalg.norm(Pw[:, None, :] - Pq[None, :, :], axis=2)
    n_pairs = int(np.sum(d < 1e-9))
    springs = m.groups[gs]
    r.check("number of interface springs = number of coincident wall/water node pairs", len(springs.elements),
            n_pairs, atol=0)
    r.require("every spring joins a wall node (I) and a coincident water node (J)", all(
        e.nodes[0] in set(wall_nodes) and e.nodes[1] in set(water_nodes)
        and np.linalg.norm(m.node_global(e.nodes[0]) - m.node_global(e.nodes[1])) < 1e-9
        for e in springs.elements.values()))

    def k_at(x, y, z):
        for e in springs.elements.values():
            if np.linalg.norm(m.node_global(e.nodes[1]) - np.array([x, y, z])) < 1e-9:
                return m.springs[e.prop].k
        return None
    zmid = 1.0
    cases = [("x-wall node (normal X)", (0.0, 1.0, zmid), (stiff, stiff2, stiff2)),
             ("y-wall node (normal Y)", (2.0, b, zmid), (stiff2, stiff, stiff2)),
             ("floor node (normal Z)", (2.0, 1.0, 0.0), (stiff2, stiff2, stiff)),
             ("vertical corner (normals X, Y)", (a, b, zmid), (stiff, stiff, stiff2)),
             ("x-wall bottom edge (normals X, Z)", (a, 2.0, 0.0), (stiff, stiff2, stiff)),
             ("bottom corner (normals X, Y, Z)", (0.0, 0.0, 0.0), (stiff, stiff, stiff)),
             ("x-wall node at the water surface (no Z restraint)", (0.0, 1.0, H - m_empty * h), (stiff, stiff2, stiff2))]
    for label, p, ref in cases:
        k = k_at(*p)
        ok = k is not None and np.allclose(k[:3], ref, rtol=0, atol=0) and not any(k[3:])
        r.require(f"spring constants at the {label} = {ref} (Stiff along normals, stiff2 along the wall)", ok,
                  f"got {k}")
    mat = m.materials[water.elements[1].mat]
    c = mat.constants(m.gravity)
    r.check("water material: bulk modulus K = M - 4G/3 (kPa) = 2.2e6", c.M - 4.0 * c.G / 3.0, 2.2e6, rtol=1e-8)
    r.check("water material: G/K = 1e-8", c.G / (c.M - 4.0 * c.G / 3.0), 1.0e-8, rtol=1e-6)
    r.check("water mass = rho V (t) with rho = 9.81/9.81", c.rho * vol, a * b * (H - m_empty * h), rtol=1e-12)
    shells = m.groups[gsh]
    area = sum(quad_area(element_xyz(m, gsh, e))[0] for e in shells.elements)
    depth = H - m_empty * h
    r.check("interface-area shells: total area = floor + perimeter x depth", area, a * b + 2 * (a + b) * depth,
            rtol=1e-12)
    r.require("interface-area shells use wall nodes only",
              all(n in set(wall_nodes) for e in shells.elements.values() for n in e.nodes))
    r.require("water nodes with springs have their rotations fixed",
              all(m.nodes[e.nodes[1]].fix[3:6] == [1, 1, 1] for e in springs.elements.values()))
    r.require("MOPT <incomp> = 0 (incompatible modes for the water SOLIDs)", int(m.mopt.get("incomp")) == 0)
    # WRITE -> INP round trip (UT-03) of the filled pool
    text, _ = write_pre(m)
    ui2 = _ui(workdir)
    ui2.run_text(text)
    r.require("WRITE -> INP of the filled pool gives the same model (UT-03)", m.same_state(ui2.model) and
              not _errors(ui2), "; ".join(_errors(ui2)))
    r.require("... and LISTPOOLINTER data survives (POOLDATA)", wl.pool_record(ui2.model) == rec)

    # SOLID pool: free upward slab faces are the floor template
    ui = _ui(workdir)
    ui.run_text("\n".join(["GRAVITY,9.81"] + solid_pool_lines(3.0, 2.0, 2.0, 3, 2, 4)))
    m = ui.model
    ui.execute("FILLPOOL,1e9,0,1")
    r.require("FILLPOOL on the SOLID pool runs without error", not _errors(ui), "; ".join(_errors(ui)))
    rec = wl.pool_record(m)
    if rec is not None:
        gw, gs = int(rec.water), int(rec.springs)
        vol = sum(hex_volume(element_xyz(m, gw, e)) for e in m.groups[gw].elements)
        r.check("SOLID pool: water volume = 3 x 2 x (2 - 1 x 0.5)", vol, 3.0 * 2.0 * 1.5, rtol=1e-12)
        r.check("SOLID pool: water from the slab top (z = 0)", float(rec.zfloor), 0.0, atol=1e-12)
        # springs: floor 4 x 3 = 12 plan points at z = 0 and 10 perimeter points at 3 more levels
        r.check("SOLID pool: number of interface springs (12 floor + 10 x 3 wall pairs)", len(m.groups[gs].elements),
                42, atol=0)
        walls = {e.nodes[0] for e in m.groups[gs].elements.values()}
        r.require("SOLID pool: the rotations of the SOLID-wall interface nodes are fixed (EDU-06)",
                  all(m.nodes[w].fix[3:6] == [1, 1, 1] for w in walls))


def _workflow_checks(r: VPResult, workdir: Path) -> None:
    """CUT2SUB pool -> FILLPOOL -> LISTPOOLINTER -> MERGEPOOL into the original model."""
    from sassi.prep import Kind
    lines = ["GRAVITY,9.81"] + shell_pool_lines(4.0, 2.0, 2.0, 4, 2, 2)
    # a building around the pool: a column (BEAMS) with a K node, so the original has more nodes and groups
    lines += ["N,500,6,0,0", "N,501,6,0,3", "N,502,7,0,1.5", "M,2,3e7,0.2,24,0.05,0.05", "R,1,0.25,0.2,0.2,0.01,0.005,0.005",
              "GROUP,2,BEAMS", "MACT,2", "RACT,1", "E,1,500,501,502"]
    ui = _ui(workdir)
    ui.run_text("\n".join(lines))
    last = max(ui.model.nodes)
    ui.run_text("\n".join(["CUTADD,1,1,RANGE,1,32", "CUT2SUB,1,2", "ACTM,2", f"FILLPOOL,1e9,0,0,-1,{last},0",
                           "ACTM,0"]))
    r.require("CUT2SUB + FILLPOOL of the pool sub-model run without error", not _errors(ui), "; ".join(_errors(ui)))
    pool = ui.models.get(2)
    ui.sink.clear()
    ui.execute("LISTPOOLINTER,2")
    conf = " ".join(ui.sink.texts(Kind.CONFIRM))
    from sassi.prep import water_lib as wl
    n_if = len({w for w, _ in wl.pool_interface(pool)}) if pool is not None else -1
    r.require(f"LISTPOOLINTER lists all {n_if} interface nodes (same number and position in the original model)",
              f"{n_if} of {n_if} interface nodes" in conf and not ui.sink.texts(Kind.WARNING), conf)
    ui.sink.clear()
    ui.execute("MERGEPOOL,2")
    r.require("MERGEPOOL imports the pool without error", not _errors(ui), "; ".join(_errors(ui)))
    m = ui.model
    titles = [g.title for g in m.groups.values()]
    r.require("MERGEPOOL appends the water SOLIDs and the interface springs after the original groups",
              titles[:2] == ["pool walls and floor", ""] and wl.TITLES["water"] in titles
              and wl.TITLES["springs"] in titles, str(titles))
    water_nodes = sorted({n for g in m.groups.values() if g.title == wl.TITLES["water"]
                          for e in g.elements.values() for n in e.nodes})
    r.require("the water nodes keep the FILLPOOL numbers (offset = last node of the original + 1 ...)",
              bool(water_nodes) and water_nodes[0] == last + 1)
    ui.sink.clear()
    ui.execute("MERGEPOOL,2")
    r.require("a second MERGEPOOL is refused (water node numbers already used)", bool(_errors(ui)))


@problem("VP-54W", "Water modelling and refinement: REFINEMODEL counts/area/volume, FILLPOOL volume a*b*(H - m*h), "
                   "springs, offset, LISTPOOLINTER", tier="P2", modules=["UI"],
         source="requirements 6.3 VP-54; spec 11 sections 1 and 7 items 1-2; D-WAT-01")
def vp54w(workdir) -> VPResult:
    r = VPResult()
    wd = Path(workdir)
    _refine_checks(r, wd)
    _fill_checks(r, wd)
    _workflow_checks(r, wd)
    return r


# ======================================================================================
# VP-W1: hydrodynamic (impulsive) mass of a rigid rectangular tank
# ======================================================================================
#: tank of VP-W1: half length L along X (the motion), half width B, wall height, mesh, EmptyLevels
TANK = dict(L=5.0, B=1.0, wall=6.0, nx=20, ny=2, nz=12, empty=2, stiff=1.0e9)
#: frequencies of the harmonic input (Hz) and the SITE frequency step (exact frequency numbers)
W1_FREQS = (2.0, 4.0, 8.0)
W1_FSTEP = 0.25
#: engineering tolerance on F/(m a) against the exact potential-flow impulsive mass (docs/user/WATER.md)
W1_RTOL = 0.03          # D-W3-07 (audit): aligned with the pytest wrapper; observed +1.5 to +1.9 %


def tank_lines(model: str = "tank", tank: Dict = TANK, freqs: Sequence[float] = W1_FREQS, mopt: int = -1,
               walls: str = "shell") -> List[str]:
    """The VP-W1 model as command lines: a rigid, massless tank (E = 3e10 kPa) on a practically rigid site
    (Vs = 10 km/s), filled by FILLPOOL and analysed by SITE -> POINT -> HOUSE -> ANALYS for a vertically incident
    SV wave (control motion X at the surface).

    ``walls`` = 'shell': SHELL walls and floor (t = 0.5 m), floor at grade, interaction nodes = the floor nodes;
    'solid': SOLID walls 0.5 m thick on a 0.5 m SOLID slab whose bottom is at grade (GROUNDELEV -0.5), the
    water on the slab top (z = 0), interaction nodes = the slab-bottom nodes, interface-area shells made
    (``ShellArea`` = 1).  ``mopt`` >= 0 adds ``MOPT,<mopt>`` after FILLPOOL (1 = incompatible modes
    suppressed, the volumetric-locking demonstration of docs/user/WATER.md)."""
    L, B, Hw = tank["L"], tank["B"], tank["wall"]
    nx, ny, nz = tank["nx"], tank["ny"], tank["nz"]
    fnums = ",".join(str(int(round(f / W1_FSTEP))) for f in freqs)
    lines = [f"MDL,{model},.", f"TIT,VP-W1 rigid rectangular tank ({walls} walls) with FILLPOOL water on a stiff site"]
    if walls == "shell":
        lines += ["HOUSE,9.81,0,0,2,0,0,0,0,0"]
        lines += shell_pool_lines(2 * L, 2 * B, Hw, nx, ny, nz, x0=-L, y0=-B, z0=0.0)
        lines += [f"FILLPOOL,{tank['stiff']:.12g},0,{tank['empty']},-1,0,0"]
        nbase = (nx + 1) * (ny + 1)
    elif walls == "solid":
        tw = tf = 0.5
        lines += [f"HOUSE,9.81,{-tf:g},0,2,0,0,0,0,0"]
        lines += solid_pool_lines(2 * L, 2 * B, Hw, nx, ny, nz, tw=tw, tf=tf)
        lines += [f"FILLPOOL,{tank['stiff']:.12g},0,{tank['empty']},1,0,0"]
        nbase = (nx + 3) * (ny + 3)
    else:
        raise ValueError("walls must be 'shell' or 'solid'")
    if mopt >= 0:
        lines.append(f"MOPT,{mopt}")
    lines += ["FIXROT", f"INT,1,{nbase},1,1",
              "L,1,1.0,20.0,17320.508,10000.0,0.01,0.01", "TOPL,1,1,1,1,1",
              f"SITE,0,1,{W1_FSTEP},20,1,1,0,1,4096,1,0,0.01,1024,1", "WAVE,2,1,1,1,0",
              "POINT,0,0,0.6", "ANALYS,0,0,0,0,1,0,0,0,0,0,0", f"FREQ,1,{fnums}",
              "AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0", "CHECK", "AFWRITE", "RUNSITE", "RUNPOINT", "RUNHOUSE",
              "RUNANALYS"]
    return lines


def tank_response(workdir: Path, water_group: int, model: str = "tank", zfloor: float = 0.0) -> Dict[str, np.ndarray]:
    """Hydrodynamic response of the VP-W1 tank from the module files: per frequency the base shear (sum of the
    X forces of the FILLPOOL springs from the FILE4 recovery operators and FILE8), the inertia of the water
    (the SOLIDs of ``water_group``) ``w^2 sum (M_s H)_x`` over the water X equations, the height of the
    resultant of the X spring forces, the largest |H_x - 1| of the tank nodes and the water mass of the HOUSE
    mass matrix (COOSM)."""
    from sassi.io.files import read_container
    wd = Path(workdir)
    f4 = read_container(wd / f"{model}.N4", "FILE4")
    f8 = read_container(wd / "FILE8", "FILE8")
    Ms = read_container(wd / "COOSM", "COOSM").sparse("Ms").tocsr()
    H = np.asarray(f8["H"])
    freq = np.asarray(f8["freq"], float)
    eq_node, eq_dof = np.asarray(f4["eq_node"]), np.asarray(f4["eq_dof"])
    gid = np.asarray(f4["elem_group"])
    types = np.asarray(f4["elem_type"])
    nodes = np.asarray(f4["elem_nodes"])
    solid_nodes = {int(n) for k in np.flatnonzero((types == 1) & (gid == int(water_group))) for n in nodes[k] if n}
    water_x = np.flatnonzero(np.isin(eq_node, sorted(solid_nodes)) & (eq_dof == 1))
    tank_x = np.flatnonzero(~np.isin(eq_node, sorted(solid_nodes)) & (eq_dof == 1))
    mass = float(np.asarray(Ms[water_x].sum(axis=1)).sum())          # rows of the water X equations
    S = np.asarray(f4["rec_SPRING_S"])
    E = np.asarray(f4["rec_SPRING_eq"])
    idx = np.asarray(f4["rec_SPRING_idx"])
    xyz = {int(n): np.asarray(p, float) for n, p in zip(f4["node_id"], f4["node_xyz"])}
    zI = np.array([xyz[int(nodes[k][0])][2] for k in idx]) - float(zfloor)
    w2 = (2.0 * np.pi * freq) ** 2
    shear, inertia, height, dev = [], [], [], []
    for q in range(len(freq)):
        U = np.where(E >= 0, H[q][np.clip(E, 0, None)], 0.0)        # (nS, 12)
        F = np.einsum("scd,sd->sc", S, U)                           # (nS, 6) FX..MZZ
        Fx = F[:, 0]
        shear.append(Fx.sum())
        inertia.append(w2[q] * (Ms[water_x] @ H[q]).sum())
        height.append((Fx * zI).sum() / Fx.sum())
        dev.append(np.max(np.abs(H[q][tank_x] - 1.0)))
    return dict(freq=freq, shear=np.array(shear), inertia=np.array(inertia), height=np.array(height),
                dev=np.array(dev), mass=mass, w2=w2)


def run_tank(workdir: Path, model: str = "tank", mopt: int = -1, walls: str = "shell"):
    """Build and run the VP-W1 tank through the interpreter; returns (interpreter, error messages)."""
    from sassi.prep import Interpreter, Kind
    Path(workdir).mkdir(parents=True, exist_ok=True)
    ui = Interpreter(cwd=workdir)
    ui.run_text("\n".join(tank_lines(model, mopt=mopt, walls=walls)), name=f"{model}.pre")
    return ui, ui.sink.texts(Kind.ERROR)


#: binary module files deleted after a VP-W1 variant has been evaluated (disk space; listings are kept)
_W1_BINARIES = ("FILE1", "FILE2", "FILE3", "FILE8", "FILE90", "FILE91", "COOSK", "COOSM", "DOFSMAP", "COOX", "COOTK")


def _purge(wd: Path) -> None:
    for p in wd.iterdir():
        if p.is_file() and (p.name.startswith(_W1_BINARIES) or p.suffix.upper() == ".N4"):
            p.unlink()


def _w1_variant(r: VPResult, wd: Path, walls: str) -> None:
    from sassi.prep import water_lib as wl
    L, B, Hw, dz = TANK["L"], TANK["B"], TANK["wall"], TANK["wall"] / TANK["nz"]
    Hwater = Hw - TANK["empty"] * dz
    tag = f"{walls} walls"
    ui, errs = run_tank(wd, model="tank", walls=walls)
    r.require(f"{tag}: tank model, FILLPOOL, CHECK, AFWRITE and SITE/POINT/HOUSE/ANALYS run without error", not errs,
              "; ".join(errs)[-800:])
    if errs:
        return
    m = ui.model
    rec = wl.pool_record(m)
    c = m.materials[int(rec.watermat)].constants(m.gravity)
    m_ref = c.rho * (2 * L) * (2 * B) * Hwater
    res = tank_response(wd, int(rec.water), zfloor=float(rec.zfloor))
    _purge(wd)
    if walls == "solid":
        walls_if = {e.nodes[0] for e in m.groups[int(rec.springs)].elements.values()}
        r.require(f"{tag}: rotations of the SOLID-wall interface nodes fixed, interface-area shells made",
                  all(m.nodes[w].fix[3:6] == [1, 1, 1] for w in walls_if) and int(rec.shells) > 0)
    r.check(f"{tag}: water mass in the HOUSE mass matrix = rho 2L 2B H = {m_ref:g} t", res["mass"], m_ref, rtol=1e-10)
    ratio_exact = wl.impulsive_ratio_exact(L, Hwater)
    ratio_housner = wl.housner_impulsive_ratio(L, Hwater)
    ratio_west = wl.westergaard_ratio(L, Hwater)
    h_exact = wl.impulsive_height_exact(L, Hwater)
    for q, f in enumerate(res["freq"]):
        Fx, In = res["shear"][q], res["inertia"][q]
        r.check(f"{tag}, {f:g} Hz: base shear = inertia of the water (dynamic equilibrium, relative)",
                abs(Fx - In) / abs(Fx), 0.0, atol=1e-6)
        r.check(f"{tag}, {f:g} Hz: rigid tank follows the control motion, max |H_x - 1| of the tank nodes",
                res["dev"][q], 0.0, atol=1e-3)
        ratio = abs(Fx) / (res["w2"][q] * m_ref)
        r.check(f"{tag}, {f:g} Hz: F/(m a) = exact potential-flow impulsive mass ratio m_i/m (L/H = {L / Hwater:g})",
                ratio, ratio_exact, rtol=W1_RTOL)
        if walls == "shell":
            r.inform(f"{tag}, {f:g} Hz: F/(m a) vs Housner (1963) tanh(sqrt3 L/H)/(sqrt3 L/H)", ratio, ratio_housner)
            r.inform(f"{tag}, {f:g} Hz: F/(m a) vs Westergaard (1933) (7/12) H/L (long-reservoir estimate)", ratio,
                     ratio_west)
            r.inform(f"{tag}, {f:g} Hz: height of the wall-pressure resultant h_i/H vs exact series (EBP)",
                     float(np.real(res["height"][q])) / Hwater, h_exact)
            r.inform(f"{tag}, {f:g} Hz: height of the wall-pressure resultant h_i/H vs Housner 3/8 (EBP)",
                     float(np.real(res["height"][q])) / Hwater, wl.HOUSNER_HEIGHT_EBP)
    r.notes.append(f"{tag}: tank 2L x 2B x H = {2 * L:g} x {2 * B:g} x {Hwater:g} m of water (walls {Hw:g} m, "
                   f"EmptyLevels {TANK['empty']}), {TANK['nx']} x {TANK['ny']} x {int(round(Hwater / dz))} water SOLIDs, "
                   f"interface Stiff = {TANK['stiff']:g} kN/m; m = {m_ref:g} t; exact m_i/m = {ratio_exact:.6f}, Housner "
                   f"{ratio_housner:.6f}, Westergaard {ratio_west:.6f}; F/(m a) = "
                   + ", ".join(f"{abs(F) / (w2 * m_ref):.5f} ({f:g} Hz)" for F, w2, f in
                               zip(res["shear"], res["w2"], res["freq"])))


@problem("VP-W1", "Hydrodynamic (impulsive) mass of a rigid rectangular tank filled by FILLPOOL (stiff site, harmonic "
                  "input)", tier="P2", modules=["UI", "SITE", "POINT", "HOUSE", "ANALYS"],
         source="Westergaard (1933); Housner (1963); exact potential-flow series (sassi.prep.water_lib); D-WAT-01")
def vpw1(workdir) -> VPResult:
    r = VPResult()
    wd = Path(workdir)
    for walls in ("shell", "solid"):
        _w1_variant(r, wd / walls, walls)
    return r
