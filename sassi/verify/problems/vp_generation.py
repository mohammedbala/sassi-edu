"""VP-50 (INTGEN counts) and VP-51 (generation tools): requirements section 6.3, spec 09 section 8.

VP-50 -- INTGEN on a regular box excavation of n_x x n_y x n_z hexahedra (ETYPE 2, ground surface at
the top face).  Exact closed-form counts (spec 09 section 2.19):

* FV    (n_x+1)(n_y+1)(n_z+1)                 every excavation node
* EVBN  FV - (n_x-1)(n_y-1)(n_z-1)            boundary of the excavated volume (top face included)
* FSIN  EVBN - (n_x-1)(n_y-1)                 lateral + bottom faces (interior top-face nodes out)
* FFV   EVBN + (n_x-1)(n_y-1) * L_sel         L_sel = floor((n_z-1)/(skip+1)) internal levels
* surface (n_x+1)(n_y+1)                      foundation nodes at grade

options 1-5 are additive (union with the existing interaction nodes) and option 0 clears them.  The
FSIN and EVBN sets are also compared node by node with the analytic face sets.

VP-51 -- generation tools (exact counts, coordinates and number maps): MERGE offsets/translation;
MERGESOIL modes 1, 2, 3 (coincident interface nodes collapse to the lower number, the map lists every
excavation node; n springs of 1e7; split at SepLevel 1e7 / 10); EXCAV of a 4 x 4 basemat with 3
levels and z jitter +-0.001 (delta 0 -> extra levels; delta 0.01 -> 32 hexes and 75 nodes, D-W2-05; with
delta the excavation still joins all n interface nodes in MERGESOIL); WELD of
two adjacent 8-node cubes -> 12 connected nodes, 4 unused; NCOM after RMVUNUSED remaps masses,
fixities and interaction flags; ROTATE (1,0,0) by rxy = 90 -> (0,1,0); FIXROT on a flat XY plate ->
ROTZ fixed only, the plate rotated 30 deg about X -> drilling springs; FIXROT then ROTATE -> warning
that the global-axis ROTZ fixities are not rotated.

Note on the EXCAV reference: a stacked mesh on a 4 x 4-cell template (25 plan points) with L levels
has 16 (L-1) hexahedra and 25 L nodes, so "48 hexes and 75 nodes" cannot both hold (75 nodes -> L = 3
-> 32 hexes; 48 hexes -> L = 4 -> 100 nodes).  The model below has 3 levels as stated ("3 levels");
lead decision D-W2-05 fixed the reference at 32 hexahedra and 75 nodes, which is what is checked.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from .. import VPResult, problem


# ======================================================================================
# Model builders (also used by the unit tests)
# ======================================================================================
def box_node(nx: int, ny: int, i: int, j: int, k: int) -> int:
    """Node number of grid point (i, j, k) of a box, numbered bottom-up."""
    return 1 + i + (nx + 1) * j + (nx + 1) * (ny + 1) * k


def box_lines(nx: int, ny: int, nz: int, h: float = 5.0, etype: int = 2, group: int = 1,
              gelev: float = 0.0) -> List[str]:
    """n_x x n_y x n_z hexahedra of size h, top face at ``gelev``; nodes numbered bottom-up."""
    out = []
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                out.append(f"N,{box_node(nx, ny, i, j, k)},{h * i},{h * j},{gelev - h * nz + h * k}")
    out.append(f"GROUP,{group},SOLID")
    e = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                e += 1
                b = lambda a, c, d: box_node(nx, ny, a, c, d)   # noqa: E731
                ids = [b(i, j, k), b(i + 1, j, k), b(i + 1, j + 1, k), b(i, j + 1, k),
                       b(i, j, k + 1), b(i + 1, j, k + 1), b(i + 1, j + 1, k + 1), b(i, j + 1, k + 1)]
                out.append("E," + ",".join(str(x) for x in [e] + ids))
    if etype:
        out.append(f"ETYPE,1,{e},1,{etype}")
    out.append(f"GROUNDELEV,{gelev}")
    return out


def plate_lines(nx: int, ny: int, h: float = 5.0, z: float = 0.0, first: int = 1, group: int = 1) -> List[str]:
    """Flat n_x x n_y SHELL plate in the XY plane at elevation z (nodes from ``first``)."""
    out = []
    nid = lambda i, j: first + i + (nx + 1) * j       # noqa: E731
    for j in range(ny + 1):
        for i in range(nx + 1):
            out.append(f"N,{nid(i, j)},{h * i},{h * j},{z}")
    out.append(f"GROUP,{group},SHELL")
    e = 0
    for j in range(ny):
        for i in range(nx):
            e += 1
            out.append(f"E,{e},{nid(i, j)},{nid(i + 1, j)},{nid(i + 1, j + 1)},{nid(i, j + 1)}")
    out += ["M,1,4.32e6,0.2,0.15,0.05,0.05", f"THICK,1,{e},1,1.0"]
    return out


def basement_lines(n: int = 2, h: float = 5.0, levels: Tuple[float, ...] = (-10.0, -5.0, 0.0),
                   floors: bool = False, jitter: float = 0.0) -> List[str]:
    """Shell basement: n x n mat at levels[0], perimeter walls up to levels[-1] and (optionally) floor
    slabs at the other levels.  Nodes: level k, grid (i, j) -> 1 + i + (n+1) j + (n+1)^2 k (perimeter
    nodes only on wall levels without floors).  ``jitter`` adds +-jitter to z, alternating in a
    checkerboard pattern."""
    out = ["M,1,4.32e6,0.2,0.15,0.05,0.05"]
    m1 = n + 1
    nid = lambda i, j, k: 1 + i + m1 * j + m1 * m1 * k   # noqa: E731

    def on_perimeter(i, j):
        return i in (0, n) or j in (0, n)

    for k, z in enumerate(levels):
        for j in range(m1):
            for i in range(m1):
                if k == 0 or floors or on_perimeter(i, j):
                    dz = jitter * (1 if (i + j + k) % 2 == 0 else -1)
                    out.append(f"N,{nid(i, j, k)},{h * i},{h * j},{z + dz}")
    out.append("GROUP,1,SHELL")
    e = 0
    for k in range(len(levels)):
        if k > 0 and not floors:
            continue
        for j in range(n):
            for i in range(n):
                e += 1
                out.append(f"E,{e},{nid(i, j, k)},{nid(i + 1, j, k)},{nid(i + 1, j + 1, k)},{nid(i, j + 1, k)}")
    # perimeter walls: CCW loop of perimeter points
    loop = [(i, 0) for i in range(n)] + [(n, j) for j in range(n)] + [(i, n) for i in range(n, 0, -1)] + \
           [(0, j) for j in range(n, 0, -1)]
    for k in range(len(levels) - 1):
        for a in range(len(loop)):
            (i1, j1), (i2, j2) = loop[a], loop[(a + 1) % len(loop)]
            e += 1
            out.append(f"E,{e},{nid(i1, j1, k)},{nid(i2, j2, k)},{nid(i2, j2, k + 1)},{nid(i1, j1, k + 1)}")
    out += [f"THICK,1,{e},1,0.5", "ETYPE,1,{},1,1".format(e), f"GROUNDELEV,{levels[-1]}"]
    return out


def _ui(lines: List[str], workdir=None):
    from ...prep import Interpreter
    ui = Interpreter(cwd=workdir)
    ui.run_text("\n".join(lines))
    return ui


def _n_int(m) -> int:
    return sum(1 for nd in m.nodes.values() if 0 in nd.flags)


def _int_set(m) -> set:
    return {i for i, nd in m.nodes.items() if 0 in nd.flags}


# ======================================================================================
# VP-50
# ======================================================================================
def box_reference(nx: int, ny: int, nz: int) -> Dict[str, int]:
    fv = (nx + 1) * (ny + 1) * (nz + 1)
    evbn = fv - (nx - 1) * (ny - 1) * (nz - 1)
    fsin = evbn - (nx - 1) * (ny - 1)
    ref = dict(FV=fv, EVBN=evbn, FSIN=fsin, SURFACE=(nx + 1) * (ny + 1))
    for skip in (0, 1, 2):
        ref[f"FFV{skip}"] = evbn + (nx - 1) * (ny - 1) * ((nz - 1) // (skip + 1))
    return ref


def box_face_sets(nx: int, ny: int, nz: int) -> Tuple[set, set]:
    """Analytic EVBN and FSIN node sets of the box (FSIN: lateral or bottom face)."""
    evbn, fsin = set(), set()
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                lat = i in (0, nx) or j in (0, ny)
                n = box_node(nx, ny, i, j, k)
                if lat or k in (0, nz):
                    evbn.add(n)
                if lat or k == 0:
                    fsin.add(n)
    return evbn, fsin


@problem("VP-50", "INTGEN interaction-node counts on an nx x ny x nz box (FV, EVBN, FSIN, FFV, surface)",
         tier="P1", modules=["UI"], source="requirements 6.3 VP-50; spec 09 sections 2.19 and 8 item 1")
def vp50(workdir):
    r = VPResult()
    for nx, ny, nz in ((3, 3, 3), (4, 3, 5), (5, 4, 6)):
        tag = f"{nx}x{ny}x{nz}"
        ref = box_reference(nx, ny, nz)
        ui = _ui(box_lines(nx, ny, nz))
        m = ui.model
        for line, key in (("INTGEN,1", "FV"), ("INTGEN,2", "EVBN"), ("INTGEN,3", "FSIN"), ("INTGEN,5,0", "FFV0"),
                          ("INTGEN,5", "FFV1"), ("INTGEN,5,2", "FFV2"), ("INTGEN,4", "SURFACE")):
            ui.execute("INTGEN,0")
            ui.execute(line)
            r.check(f"{tag} {line} ({key}) count", _n_int(m), ref[key], atol=0)
        evbn, fsin = box_face_sets(nx, ny, nz)
        ui.execute("INTGEN,0")
        ui.execute("INTGEN,2")
        r.require(f"{tag} EVBN set = analytic boundary-node set", _int_set(m) == evbn)
        ui.execute("INTGEN,0")
        ui.execute("INTGEN,3")
        r.require(f"{tag} FSIN set = analytic lateral + bottom set", _int_set(m) == fsin)
        # additivity (union) and clearing
        ui.execute("INTGEN,4")
        r.check(f"{tag} INTGEN,3 then 4: FSIN + top interior = EVBN", _n_int(m), ref["EVBN"], atol=0)
        ui.execute("INTGEN,5")
        r.check(f"{tag} then INTGEN,5 adds the internal level(s)", _n_int(m), ref["FFV1"], atol=0)
        ui.execute("INTGEN,1")
        r.check(f"{tag} then INTGEN,1 gives FV", _n_int(m), ref["FV"], atol=0)
        ui.execute("INTGEN,3")
        r.check(f"{tag} INTGEN,3 after FV removes nothing", _n_int(m), ref["FV"], atol=0)
        ui.execute("INTGEN,0")
        r.check(f"{tag} INTGEN,0 clears", _n_int(m), 0, atol=0)
    # surface foundation: shell mat at grade, no excavation (option 4 needs no ETYPE)
    nx, ny = 4, 3
    ui = _ui(plate_lines(nx, ny) + ["GROUNDELEV,0"])
    ui.execute("INTGEN,4")
    r.check(f"surface mat {nx}x{ny}: INTGEN,4 count", _n_int(ui.model), (nx + 1) * (ny + 1), atol=0)
    # an existing interaction node outside the excavation is kept (options 1-5 never remove)
    nx, ny, nz = 3, 3, 3
    ui = _ui(box_lines(nx, ny, nz) + ["N,500,0,0,5", "INT,500,500,1,1"])
    ui.execute("INTGEN,2")
    r.check("pre-existing interaction node kept: EVBN + 1", _n_int(ui.model), box_reference(nx, ny, nz)["EVBN"] + 1,
            atol=0)
    # implicit ETYPE 0 below grade is classified by D-HOU-01 (with a warning)
    ui = _ui(box_lines(nx, ny, nz, etype=0))
    ui.execute("INTGEN,2")
    r.check("implicit ETYPE 0 box: EVBN count", _n_int(ui.model), box_reference(nx, ny, nz)["EVBN"], atol=0)
    # buried shell (ETYPE 2) outside the excavation: its nodes are added by every option
    ui = _ui(box_lines(nx, ny, nz) + plate_lines(2, 1, h=5.0, z=-5.0, first=600, group=2)
             + ["ETYPE,1,2,1,2"])
    ui.execute("INTGEN,3")
    r.check("FSIN + buried-shell nodes (6 nodes, 2 shared with the box face)", _n_int(ui.model),
            box_reference(nx, ny, nz)["FSIN"] + 6, atol=0)
    # no excavation: options 1, 2, 3, 5 refuse, nothing is flagged
    ui = _ui(box_lines(nx, ny, nz, etype=1))
    ok = ui.execute("INTGEN,1")
    r.require("structure-only model: INTGEN,1 fails and flags nothing", (not ok) and _n_int(ui.model) == 0)
    return r


# ======================================================================================
# VP-51
# ======================================================================================
def _cube(first: int, x0: float, group: int, mat: int = 1) -> List[str]:
    """Unit cube nodes first..first+7 at x0..x0+1 and one SOLID element in ``group``."""
    pts = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    out = [f"N,{first + k},{x0 + p[0]},{p[1]},{p[2]}" for k, p in enumerate(pts)]
    out += [f"GROUP,{group},SOLID", f"MACT,{mat}", "E,1," + ",".join(str(first + k) for k in range(8))]
    return out


@problem("VP-51", "Generation tools (MERGE, MERGESOIL modes, EXCAV, WELD, NCOM after RMVUNUSED, ROTATE, FIXROT)",
         tier="P1", modules=["UI"], source="requirements 6.3 VP-51; spec 09 section 8 items 2-9")
def vp51(workdir):
    from ...prep import Interpreter, Kind
    r = VPResult()
    workdir = Path(workdir)

    # ---- 1. MERGE: offsets by the Mdl1 maxima, translation of Mdl2, L table shared
    ui = Interpreter(cwd=workdir)
    ui.run_text("\n".join(["ACTM,1"] + _cube(1, 0.0, 1) + ["M,1,4e6,0.2,0.15,0.05,0.05", "L,1,10,0.12,1500,800,0.05,0.05",
                                                           "N,9,3,0,0", "N,10,4,0,0", "N,11,3,1,0", "GROUP,4,BEAMS",
                                                           "MACT,1", "RACT,1", "E,1,9,10,11", "R,1,1,1,1,1,1,1"]
                          + ["ACTM,2"] + _cube(1, 0.0, 1, mat=1) + ["M,1,3e6,0.25,0.14,0.04,0.04",
                                                                   "L,1,10,0.12,1500,800,0.05,0.05", "MT,8,7,7,7",
                                                                   "INT,1,1,1,1"]
                          + ["ACTM,3", "MERGE,1,2,100,0,-5"]))
    m = ui.model
    r.check("MERGE: node count", len(m.nodes), 11 + 8, atol=0)
    r.check("MERGE: Mdl2 node 1 -> 12 (offset = max node 11)", min(set(m.nodes) - set(range(1, 12))), 12, atol=0)
    r.require("MERGE: Mdl2 node 1 translated to (100, 0, -5)", np.allclose(m.node_global(12), [100.0, 0.0, -5.0],
                                                                           rtol=0, atol=1e-12))
    r.require("MERGE: Mdl2 group 1 -> 5 (offset = max group 4)", sorted(m.groups) == [1, 4, 5])
    r.check("MERGE: Mdl2 material 1 -> 2 (offset = max material 1)", m.groups[5].elements[1].mat, 2, atol=0)
    r.require("MERGE: Mdl2 element connectivity remapped", m.groups[5].elements[1].nodes == list(range(12, 20)))
    r.require("MERGE: L table not offset", sorted(m.layers) == [1])
    r.require("MERGE: mass and interaction flag remapped", 19 in m.tmass and 0 in m.nodes[12].flags
              and 0 not in m.nodes[1].flags)

    # ---- 2.-4. MERGESOIL: shell basement 2 x 2 cells, levels -10 / -5 / 0, excavation from EXCAV
    base = basement_lines(2, 5.0, (-10.0, -5.0, 0.0))
    soil = ["L,1,5,0.12,1500,800,0.05,0.05", "L,2,5,0.12,1500,800,0.05,0.05", "L,3,20,0.13,2000,1000,0.02,0.02",
            "TOPL,1,2,3"]
    n_pairs = 9 + 8 + 8            # mat grid + perimeter at -5 + perimeter at grade
    for mode, extra in ((1, ""), (2, ""), (3, ",,,-5")):
        ui = Interpreter(cwd=workdir)
        ui.run_text("\n".join(["ACTM,1"] + base + soil + ["EXCAV,2", "ACTM,2", "ETYPEGEN,2", "ACTM,3",
                                                          f"MERGESOIL,1,2,{mode}{extra}" +
                                                          (",excv.map" if mode == 1 else "")]))
        m = ui.model
        soil_g = max(g for g, grp in m.groups.items() if grp.type == 1)
        if mode == 1:
            used = set(m.used_nodes())
            soil_nodes = {n for e in m.groups[soil_g].elements.values() for n in e.nodes}
            struct_ids = set(ui.models[1].nodes)
            r.check("MERGESOIL 1: soil-element nodes that are structure nodes", len(soil_nodes & struct_ids), n_pairs,
                    atol=0)
            r.check("MERGESOIL 1: unused higher-numbered soil nodes", len(set(m.nodes) - used), n_pairs, atol=0)
            P = {i: tuple(np.round(m.node_global(i), 9)) for i in m.nodes}
            unused = set(m.nodes) - used
            r.require("MERGESOIL 1: every unused node is a soil node coincident with a lower-numbered "
                      "structure node used by the soil elements",
                      all(any(P[t] == P[u] and t < u and t in soil_nodes for t in struct_ids) for u in unused)
                      and min(unused) > max(struct_ids))
            lines = [ln for ln in (workdir / "excv.map").read_text().splitlines() if ln and not ln.startswith("#")]
            r.check("MERGESOIL 1: map lines = excavation nodes", len(lines), len(ui.models[2].nodes), atol=0)
            mp = dict(tuple(int(x) for x in ln.split()) for ln in lines)
            r.require("MERGESOIL 1: map sends interface nodes to structure nodes",
                      sum(1 for v in mp.values() if v in struct_ids) == n_pairs)
            r.require("MERGESOIL 1: soil elements ETYPE 2, MSET = soil layers",
                      all(e.etype == 2 and e.mat in (1, 2) for e in m.groups[soil_g].elements.values()))
            ui.execute("RMVUNUSED")
            r.check("MERGESOIL 1 + RMVUNUSED: nodes left", len(m.nodes), len(ui.models[1].nodes) + 27 - n_pairs,
                    atol=0)
        else:
            sg = [g for g, grp in m.groups.items() if grp.type == 7]
            r.require(f"MERGESOIL {mode}: one SPRING group", len(sg) == 1)
            els = list(m.groups[sg[0]].elements.values()) if sg else []
            r.check(f"MERGESOIL {mode}: spring count = interface pairs", len(els), n_pairs, atol=0)
            k = [m.springs[e.prop].scx for e in els]
            if mode == 2:
                r.require("MERGESOIL 2: all springs 1e7 (translational, rotation 0)",
                          all(m.springs[e.prop].k == (1e7, 1e7, 1e7, 0.0, 0.0, 0.0) for e in els))
            else:
                r.check("MERGESOIL 3: stiff springs (z <= SepLevel -5)", sum(1 for x in k if x == 1e7), 9 + 8, atol=0)
                r.check("MERGESOIL 3: soft springs (above SepLevel)", sum(1 for x in k if x == 10.0), 8, atol=0)

    # ---- 5. EXCAV: 4 x 4 basemat, 3 levels (-20, -10, 0 = grade), z jitter +-0.001
    jit = basement_lines(4, 5.0, (-20.0, -10.0, 0.0), floors=True, jitter=0.001) + \
        ["L,1,10,0.12,1500,800,0.05,0.05", "L,2,10,0.13,2000,1000,0.02,0.02", "L,3,30,0.14,3000,1500,0.02,0.02",
         "TOPL,1,2,3"]
    ui = Interpreter(cwd=workdir)
    ui.run_text("\n".join(["ACTM,1"] + jit + ["EXCAV,2,0"]))
    lev0 = len({round(n.z, 9) for n in ui.models[2].nodes.values()})
    r.require(f"EXCAV delta 0: extra levels ({lev0} > 3)", lev0 > 3)
    ui.execute("EXCAV,2,0.01")
    ex = ui.models[2]
    hexes = sum(len(g.elements) for g in ex.groups.values())
    r.check("EXCAV delta 0.01: hexes", hexes, 32, atol=0)   # D-W2-05: 16 (L-1) with L = 3 node levels
    r.check("EXCAV delta 0.01: nodes", len(ex.nodes), 75, atol=0)
    r.check("EXCAV delta 0.01: levels", len({round(n.z, 9) for n in ex.nodes.values()}), 3, atol=0)
    tpl = {(5.0 * i, 5.0 * j) for i in range(5) for j in range(5)}
    r.require("EXCAV: every node lies on the template (x, y) set",
              all((round(n.x, 9), round(n.y, 9)) in tpl for n in ex.nodes.values()))
    r.require("EXCAV: ETYPE 2 and MSET from the TOPL layer at mid-height (2 below -10, 1 above)",
              all(e.etype == 2 for e in ex.groups[1].elements.values())
              and sorted({e.mat for e in ex.groups[1].elements.values()}) == [1, 2])
    # EXCAV delta + MERGESOIL: the structure nodes scatter by +-0.001 around the snapped levels; the
    # excavation carries EDUOPT,GEOMTOL = delta and MERGESOIL still joins every interface node (spec 09
    # section 4.1 step 5 "MERGESOIL later welds coincident nodes to the structure")
    ui = Interpreter(cwd=workdir)
    ui.run_text("\n".join(["ACTM,1"] + basement_lines(2, 5.0, (-10.0, -5.0, 0.0), jitter=0.001) + soil
                           + ["EXCAV,2,0.01", "ACTM,2", "ETYPEGEN,2", "ACTM,3", "MERGESOIL,1,2,1"]))
    m = ui.model
    soil_g = max(g for g, grp in m.groups.items() if grp.type == 1)
    shared = {n for e in m.groups[soil_g].elements.values() for n in e.nodes} & set(ui.models[1].nodes)
    r.check("EXCAV delta 0.01 + MERGESOIL 1: interface nodes joined (z scatter +-0.001)", len(shared), n_pairs,
            atol=0)
    if hexes != 32:
        r.notes.append(f"EXCAV: observed {hexes} hexahedra and {len(ex.nodes)} nodes on 3 levels. The reference "
                       f"'48 hexes and 75 nodes' is internally inconsistent for a 16-cell / 25-point template "
                       f"(16 (L-1) hexes, 25 L nodes): 75 nodes -> L = 3 -> 32 hexes. Not weakened; reported.")

    # ---- 6. WELD: two adjacent unit cubes with separate nodes on the common face
    ui = Interpreter(cwd=workdir)
    ui.run_text("\n".join(_cube(1, 0.0, 1) + _cube(9, 1.0, 2) + ["WELD"]))
    m = ui.model
    r.check("WELD: connected nodes", len(m.used_nodes()), 12, atol=0)
    r.check("WELD: unused nodes", len(set(m.nodes) - set(m.used_nodes())), 4, atol=0)
    r.require("WELD: lowest numbers kept (9, 12, 13, 16 -> 2, 3, 6, 7)",
              m.groups[2].elements[1].nodes == [2, 10, 11, 3, 6, 14, 15, 7])

    # ---- 7. NCOM after RMVUNUSED: masses, fixities and interaction flags follow their nodes
    lines = []
    ids = [10, 20, 30, 40, 50, 57, 70, 80]
    pts = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    lines += [f"N,{n},{p[0]},{p[1]},{p[2]}" for n, p in zip(ids, pts)]
    lines += ["N,5,9,9,9", "N,33,8,8,8", "N,90,7,7,7", "GROUP,1,SOLID", "E,1," + ",".join(str(n) for n in ids),
              "MT,57,12.5,13.5,14.5", "D,20,20,1,1,UX", "INT,80,80,1,1", "INT,10,10,1,1", "RMVUNUSED", "NCOM"]
    ui = Interpreter(cwd=workdir)
    ui.run_text("\n".join(lines))
    m = ui.model
    r.require("RMVUNUSED + NCOM: nodes 1..8", sorted(m.nodes) == list(range(1, 9)))
    r.require("NCOM: connectivity 1..8 in order", m.groups[1].elements[1].nodes == list(range(1, 9)))
    r.require("NCOM: mass of old node 57 at new node 6", list(m.tmass) == [6] and m.tmass[6] == [12.5, 13.5, 14.5])
    r.require("NCOM: coordinates of new node 6 = old node 57", np.allclose(m.node_global(6), [1, 0, 1], atol=0))
    r.require("NCOM: fixity of old node 20 at new node 2", m.nodes[2].fix[0] == 1 and
              sum(sum(n.fix) for n in m.nodes.values()) == 1)
    r.require("NCOM: interaction flags of old 10 / 80 at new 1 / 8", _int_set(m) == {1, 8})

    # ---- 8. ROTATE: (1,0,0) by rxy = 90 about the origin -> (0,1,0); +30 then -30 restores
    ui = Interpreter(cwd=workdir)
    ui.run_text("N,1,1,0,0\nN,2,3.7,-2.1,5.3\nROTATE,0,0,0,90,0,0")
    p = ui.model.node_global(1)
    r.check("ROTATE rxy=90: x of (1,0,0)", p[0], 0.0, atol=1e-15)
    r.check("ROTATE rxy=90: y of (1,0,0)", p[1], 1.0, atol=1e-15)
    r.check("ROTATE rxy=90: z of (1,0,0)", p[2], 0.0, atol=1e-15)
    ui.run_text("ROTATE,1,2,3,0,30,0\nROTATE,1,2,3,0,-30,0")
    r.check("ROTATE +30 / -30 about X: round trip", float(np.max(np.abs(ui.model.node_global(2) - [2.1, 3.7, 5.3]))),
            0.0, atol=1e-12)

    # ---- 9. FIXROT interplay: flat XY plate -> ROTZ only; rotated 30 deg about X -> springs
    ui = Interpreter(cwd=workdir)
    ui.run_text("\n".join(plate_lines(2, 2) + ["FIXROT"]))
    m = ui.model
    r.require("FIXROT flat XY plate: every node ROTZ fixed, nothing else",
              all(n.fix == [0, 0, 0, 0, 0, 1] for n in m.nodes.values()))
    r.check("FIXROT flat XY plate: no spring group", sum(1 for g in m.groups.values() if g.type == 7), 0, atol=0)
    ui = Interpreter(cwd=workdir)
    ui.run_text("\n".join(plate_lines(2, 2) + ["ROTATE,0,0,0,0,30,0", "FIXROT"]))
    m = ui.model
    sg = [g for g in m.groups.values() if g.type == 7]
    r.check("FIXROT plate rotated 30 deg about X: drilling springs", sum(len(g.elements) for g in sg), 9, atol=0)
    r.require("FIXROT rotated plate: no rotation fixed with D on the shell nodes",
              all(m.nodes[n].fix[3:] == [0, 0, 0] for n in range(1, 10)))
    # the other order: FIXROT fixes ROTZ (drilling) on the flat plate, then ROTATE about X by 30 deg turns
    # that global-axis fixity into a bending restraint -> ROTATE warns (fixities are not rotated)
    ui = Interpreter(cwd=workdir)
    ui.run_text("\n".join(plate_lines(2, 2) + ["FIXROT", "ROTATE,0,0,0,0,30,0"]))
    fx = [w for w in ui.sink.texts(Kind.WARNING) if "D fixities" in w]
    r.require("FIXROT then ROTATE 30 deg about X: warning on the 9 global-axis ROTZ fixities",
              len(fx) == 1 and "9 nodes" in fx[0])
    return r
