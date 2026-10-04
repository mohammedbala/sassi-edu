"""Adversarial review of the model generation and combination package (requirements 3.4.I; spec 09
sections 2 and 4; spec 04 section 7; spec 01 rules 4-11; decisions D-HOU-01, D-MDL-12..15, D-MDL-19,
D-PNT-05).

The checks are derived independently of the implementation:

* INTGEN sets of irregular (L-shaped, stepped) excavations from a voxel argument: a node is an
  excavation-boundary node (EVBN) iff it touches at least one excavated voxel and at least one
  non-excavated voxel (soil below grade or air above it); it is a foundation-soil-interface node
  (FSIN) iff one of those non-excavated voxels is soil (below grade).  The 8 voxels around a node are
  face-connected, so this is equivalent to "node of a boundary face" / "node of a non-top boundary
  face" (spec 09 section 2.19).
* Splitting every voxel into two prisms (repeated nodes) adds no node, so every INTGEN set must be
  unchanged; storing the nodes in a rotated LOC system or numbering them at random must not change
  the sets either.
* ROTATE(a, b, c) about a point must equal the three single-axis rotations applied in the D-MDL-14
  order, and an explicit Ry Rx Rz product; NCOM / WELD / MERGESOIL must not move any element node
  (geometry invariance); MERGE offsets are the table maxima, not the counts; EXCAV fills
  footprint x depth exactly, numbers nodes bottom-up and takes MSET from the TOPL position.

The tests of the last section documented defects found by the review (strict xfail); the defects are
fixed and the markers removed.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.prep import Interpreter, Kind
from sassi.prep import generation_lib as gl
from sassi.prep.writer import write_pre
from sassi.verify.problems.vp_generation import basement_lines


# ======================================================================================
# helpers
# ======================================================================================
def run(lines, ui=None):
    ui = ui or Interpreter()
    ui.run_text("\n".join(lines))
    return ui


def int_set(m):
    return {i for i, nd in m.nodes.items() if 0 in nd.flags}


def gcoords(m):
    return {i: m.node_global(i) for i in m.nodes}


def rot(axis, deg):
    """Independent right-hand rotation matrix about a global axis."""
    c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    if axis == "z":
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
    if axis == "x":
        return np.array([[1.0, 0, 0], [0, c, -s], [0, s, c]])
    return np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])


SOIL3 = ["L,1,5,0.12,1500,800,0.05,0.05", "L,2,5,0.12,1500,800,0.05,0.05", "L,3,20,0.13,2000,1000,0.02,0.02",
         "TOPL,1,2,3"]


# ======================================================================================
# INTGEN: voxel excavations (independent derivation)
# ======================================================================================
NZ = 4


def lshape_stepped_voxels():
    """Shallow L-shaped part (top 2 layers) + a deep pit (all 4 layers) under one corner."""
    vox = set()
    for i in range(5):
        for j in range(4):
            if i >= 3 and j == 3:
                continue                          # notch: L-shaped plan
            for k in (2, 3):
                vox.add((i, j, k))
    for i in range(3):
        for j in range(3):
            for k in (0, 1):
                vox.add((i, j, k))
    return vox


def voxel_reference(vox, nz=NZ, h=2.0, gelev=0.0):
    """Independent FV / EVBN / FSIN / FFV / surface sets keyed by grid point (i, j, k)."""
    corners = {}
    for (i, j, k) in vox:
        for di in (0, 1):
            for dj in (0, 1):
                for dk in (0, 1):
                    corners.setdefault((i + di, j + dj, k + dk), None)
    fv, evbn, fsin = set(corners), set(), set()
    for (i, j, k) in fv:
        around = [(i - a, j - b, k - c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
        non = [v for v in around if v not in vox]
        if non:
            evbn.add((i, j, k))
        if any(v[2] < nz for v in non):           # a soil voxel (below grade) touches the node
            fsin.add((i, j, k))
    levels = sorted({k for (_, _, k) in fv})
    ffv = {}
    for skip in (0, 1, 2):
        sel = {levels[jj] for jj in range(1, len(levels) - 1) if jj % (skip + 1) == 0}
        ffv[skip] = evbn | {p for p in fv if p[2] in sel}
    surface = {p for p in fv if p[2] == nz}
    return dict(FV=fv, EVBN=evbn, FSIN=fsin, FFV0=ffv[0], FFV1=ffv[1], FFV2=ffv[2], SURFACE=surface)


def voxel_lines(vox, nz=NZ, h=2.0, gelev=0.0, prisms=False, seed=None, loc=None, etype=2):
    """Excavation SOLIDs of the voxels; nodes numbered at random when ``seed`` is given, stored in LOC
    system 1 (``loc`` = (x0, y0, z0, txy)) when given.  Returns (lines, point -> node id)."""
    pts = sorted({(i + a, j + b, k + c) for (i, j, k) in vox for a in (0, 1) for b in (0, 1) for c in (0, 1)})
    ids = list(range(1, len(pts) + 1))
    if seed is not None:
        ids = [int(x) + 1 for x in np.random.RandomState(seed).permutation(len(pts)) * 3 + 10]
    nid = dict(zip(pts, ids))
    out = []
    if loc is not None:
        x0, y0, z0, txy = loc
        out += [f"LOC,1,0,{x0},{y0},{z0},{txy},0,0", "CSYS,1"]
    zb = gelev - h * nz - (loc[2] if loc is not None else 0.0)
    for p in pts:
        out.append(f"N,{nid[p]},{h * p[0]!r},{h * p[1]!r},{zb + h * p[2]!r}")
    out += ["CSYS,0", "GROUP,1,SOLID"]
    e = 0
    for (i, j, k) in sorted(vox):
        b = [nid[(i, j, k)], nid[(i + 1, j, k)], nid[(i + 1, j + 1, k)], nid[(i, j + 1, k)]]
        t = [nid[(i, j, k + 1)], nid[(i + 1, j, k + 1)], nid[(i + 1, j + 1, k + 1)], nid[(i, j + 1, k + 1)]]
        if prisms:
            for tri in ((0, 1, 2), (0, 2, 3)):
                e += 1
                bb, tt = [b[q] for q in tri], [t[q] for q in tri]
                out.append("E," + ",".join(str(x) for x in [e] + bb + [bb[2]] + tt + [tt[2]]))
        else:
            e += 1
            out.append("E," + ",".join(str(x) for x in [e] + b + t))
    out += [f"ETYPE,1,{e},1,{etype}", f"GROUNDELEV,{gelev}"]
    return out, nid


INTGEN_CASES = (("INTGEN,1", "FV"), ("INTGEN,2", "EVBN"), ("INTGEN,3", "FSIN"), ("INTGEN,5,0", "FFV0"),
                ("INTGEN,5", "FFV1"), ("INTGEN,5,2", "FFV2"), ("INTGEN,4", "SURFACE"))


@pytest.mark.parametrize("variant", ["hex", "prism", "scrambled", "loc", "gelev"])
def test_intgen_irregular_excavation_matches_voxel_derivation(variant):
    vox = lshape_stepped_voxels()
    ref = voxel_reference(vox)
    kw = dict(hex={}, prism=dict(prisms=True), scrambled=dict(seed=7), loc=dict(loc=(3.0, -2.0, 0.0, 37.0), seed=3),
              gelev=dict(gelev=123.5, prisms=True, seed=11))[variant]
    lines, nid = voxel_lines(vox, **kw)
    ui = run(lines)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    m = ui.model
    for cmd, key in INTGEN_CASES:
        ui.execute("INTGEN,0")
        assert ui.execute(cmd), (cmd, ui.sink.texts(Kind.ERROR)[-1:])
        got = int_set(m)
        want = {nid[p] for p in ref[key]}
        assert got == want, (variant, cmd, len(got), len(want))
    # sanity of the reference itself: FSIN is strictly smaller than EVBN for this shape
    assert ref["FSIN"] < ref["EVBN"] < ref["FV"]


def test_intgen_options_are_a_union_and_clear_touches_code_0_only():
    vox = lshape_stepped_voxels()
    ref = voxel_reference(vox)
    lines, nid = voxel_lines(vox, seed=5)
    ui = run(lines + ["N,9999,50,50,-3", "INT,9999,9999,1,1,0", "INT,9999,9999,1,1,2"])
    m = ui.model
    acc = {9999}
    for cmd, key in (("INTGEN,3", "FSIN"), ("INTGEN,4", "SURFACE"), ("INTGEN,5,2", "FFV2"), ("INTGEN,2", "EVBN")):
        ui.execute(cmd)
        acc |= {nid[p] for p in ref[key]}
        assert int_set(m) == acc, cmd
    ui.execute("INTGEN,0")
    assert int_set(m) == set() and m.nodes[9999].flags == {2}


def test_intgen_2d_plane_lshape_matches_pixel_derivation():
    """PLANE (X-Z) excavation with a step: EVBN / FSIN from the pixel argument."""
    pix = {(i, k) for i in range(6) for k in (2, 3)} | {(i, k) for i in range(3) for k in (0, 1)}
    nz, h = 4, 1.5
    pts = sorted({(i + a, k + c) for (i, k) in pix for a in (0, 1) for c in (0, 1)})
    nid = {p: n for n, p in enumerate(pts, start=1)}
    lines = [f"N,{nid[p]},{h * p[0]},0,{-h * nz + h * p[1]}" for p in pts] + ["GROUP,1,PLANE"]
    e = 0
    for (i, k) in sorted(pix):
        e += 1
        lines.append(f"E,{e},{nid[(i, k)]},{nid[(i + 1, k)]},{nid[(i + 1, k + 1)]},{nid[(i, k + 1)]}")
    lines += [f"ETYPE,1,{e},1,2", "GROUNDELEV,0"]
    evbn, fsin = set(), set()
    for (i, k) in pts:
        non = [(i - a, k - c) for a in (0, 1) for c in (0, 1) if (i - a, k - c) not in pix]
        if non:
            evbn.add(nid[(i, k)])
        if any(q[1] < nz for q in non):
            fsin.add(nid[(i, k)])
    ui = run(lines)
    ui.execute("INTGEN,2")
    assert int_set(ui.model) == evbn
    ui.execute("INTGEN,0")
    ui.execute("INTGEN,3")
    assert int_set(ui.model) == fsin
    ui.execute("INTGEN,0")
    ui.execute("INTGEN,1")
    assert int_set(ui.model) == set(nid.values())


# ======================================================================================
# RADIUS (D-PNT-05)
# ======================================================================================
def test_radius_prism_mesh_uses_half_square_plan_area(tmp_path):
    vox = {(i, j, k) for i in range(3) for j in range(2) for k in range(2)}
    h = 2.0
    lines, _ = voxel_lines(vox, nz=2, h=h, prisms=True)
    ui = run(lines, Interpreter(cwd=tmp_path))
    radii = gl.excavation_radii(ui.model, 0.85)
    assert len(radii) == 2 * len(vox)
    assert all(math.isclose(r, 0.85 * h / math.sqrt(2.0), rel_tol=1e-12) for _, _, r in radii)
    assert ui.execute("RADIUS,0.9,r.txt")
    assert math.isclose(ui.session["radius"]["average"], 0.9 * h / math.sqrt(2.0), rel_tol=1e-12)
    rows = [ln.split() for ln in (tmp_path / "r.txt").read_text().splitlines() if ln and ln[0].isdigit()]
    assert len(rows) == 2 * len(vox) and all(math.isclose(float(r[2]), 0.9 * h / math.sqrt(2.0), rel_tol=1e-7) for r in rows)


# ======================================================================================
# ROTATE / TRANSLATE / GLB2LOC
# ======================================================================================
def _points():
    rs = np.random.RandomState(1)
    return rs.uniform(-10, 10, size=(6, 3))


def test_rotate_equals_sequential_single_axis_rotations_and_explicit_matrix():
    P = _points()
    base = [f"N,{k + 1},{float(p[0])!r},{float(p[1])!r},{float(p[2])!r}" for k, p in enumerate(P)]
    c = (1.5, -2.0, 0.75)
    a, b, g = 23.0, -41.0, 67.0
    ui1 = run(base + [f"ROTATE,{c[0]},{c[1]},{c[2]},{a},{b},{g}"])
    ui2 = run(base + [f"ROTATE,{c[0]},{c[1]},{c[2]},{a},0,0", f"ROTATE,{c[0]},{c[1]},{c[2]},0,{b},0",
                      f"ROTATE,{c[0]},{c[1]},{c[2]},0,0,{g}"])
    T = rot("y", g) @ rot("x", b) @ rot("z", a)
    cc = np.array(c)
    for k, p in enumerate(P):
        want = T @ (p - cc) + cc
        assert np.allclose(ui1.model.node_global(k + 1), want, rtol=0, atol=1e-12)
        assert np.allclose(ui2.model.node_global(k + 1), want, rtol=0, atol=1e-12)


def test_rotate_signs_of_ryz_and_rzx():
    ui = run(["N,1,0,1,0", "N,2,0,0,1", "ROTATE,0,0,0,0,90,0"])
    assert np.allclose(ui.model.node_global(1), [0, 0, 1], atol=1e-15)     # y -> z about +X
    ui = run(["N,2,0,0,1", "ROTATE,0,0,0,0,0,90"])
    assert np.allclose(ui.model.node_global(2), [1, 0, 0], atol=1e-15)     # z -> x about +Y


def test_rotate_moves_local_nodes_rigidly_and_round_trips():
    lines = ["LOC,1,0,4,5,6,30,20,10", "CSYS,1", "N,1,1,2,3", "N,2,-1,0.5,2", "CSYS,0", "N,3,7,8,9",
             "N,4,7,9,9", "N,5,8,8,9", "LOCAL,2,0,3,4,5", "CSYS,2", "N,6,1,1,1", "CSYS,0"]
    ui = run(lines)
    m = ui.model
    before = gcoords(m)
    c = np.array([2.0, -1.0, 3.0])
    ui.execute(f"ROTATE,{c[0]},{c[1]},{c[2]},15,25,35")
    T = rot("y", 35) @ rot("x", 25) @ rot("z", 15)
    after = gcoords(m)
    for i in before:
        assert np.allclose(after[i], T @ (before[i] - c) + c, rtol=0, atol=1e-11), i
    text, _ = write_pre(m)
    b = run([text])
    assert not b.sink.texts(Kind.ERROR)
    for i in before:
        assert np.allclose(b.model.node_global(i), after[i], rtol=0, atol=1e-11)


def test_translate_and_merge_shift_agree_for_local_nodes():
    m2 = ["LOC,1,0,10,0,0,90,0,0", "CSYS,1", "N,1,1,2,3", "CSYS,0", "N,2,4,5,6", "GROUP,1,SPRING", "E,1,1,2"]
    ui = run(["ACTM,1", "N,1,0,0,0", "ACTM,2"] + m2 + ["ACTM,3", "MERGE,1,2,3,-4,5"])
    merged = ui.model
    ui2 = run(["ACTM,2"] + m2 + ["TRANSLATE,3,-4,5"])
    t = ui2.model
    for i in (1, 2):
        assert np.allclose(merged.node_global(i + 1), t.node_global(i), rtol=0, atol=1e-12)


def test_glb2loc_two_systems_against_independent_transformation():
    o1, a1 = np.array([1.0, 2.0, 3.0]), (30.0, 0.0, 0.0)
    o2, a2 = np.array([-4.0, 5.0, 0.5]), (10.0, 20.0, 30.0)
    R1 = rot("z", a1[0]) @ rot("x", a1[1]) @ rot("y", a1[2])              # LOC = Rz Rx Ry (D-MDL-04)
    R2 = rot("z", a2[0]) @ rot("x", a2[1]) @ rot("y", a2[2])
    loc1 = np.array([2.0, -1.0, 0.5])
    ui = run([f"LOC,1,0,{o1[0]},{o1[1]},{o1[2]},{a1[0]},{a1[1]},{a1[2]}",
              f"LOC,2,0,{o2[0]},{o2[1]},{o2[2]},{a2[0]},{a2[1]},{a2[2]}",
              "CSYS,1", f"N,1,{loc1[0]},{loc1[1]},{loc1[2]}", "CSYS,0", "N,2,7,-3,2", "GLB2LOC,1,2,1,2"])
    m = ui.model
    xg1 = o1 + R1 @ loc1
    assert m.nodes[1].csys == 2 and m.nodes[2].csys == 2
    assert np.allclose(m.nodes[1].xyz, R2.T @ (xg1 - o2), rtol=0, atol=1e-12)
    assert np.allclose(m.nodes[2].xyz, R2.T @ (np.array([7.0, -3.0, 2.0]) - o2), rtol=0, atol=1e-12)
    assert np.allclose(m.node_global(1), xg1, atol=1e-12)


# ======================================================================================
# NCOM, RMVUNUSED, WELD, MERGE, MERGESOIL: invariants
# ======================================================================================
def _element_geometry(m):
    return {(g.id, e.id): [tuple(np.round(m.node_global(n), 9)) if n else None for n in e.nodes]
            for g, e in m.iter_elements()}


def test_ncom_preserves_geometry_and_attaches_every_attribute_to_the_same_point():
    lines = ["N,3,0,0,0", "N,17,1,0,0", "N,18,1,1,0", "N,40,0,1,0", "N,41,0,0,1", "N,77,1,0,1", "N,90,1,1,1",
             "N,91,0,1,1", "N,200,5,0,0", "N,201,6,0,0", "N,202,5,1,0",
             "GROUP,1,SOLID", "E,1,3,17,18,40,41,77,90,91", "GROUP,2,BEAMS", "E,1,200,201,202",
             "MT,90,1,2,3", "MR,201,4,5,6", "D,17,17,1,1,UX", "INT,41,41,1,1", "F,77,0,1,0", "NOUT,1,1,0,0,0,0,0,18",
             "RDND,40,1,0,0,0,0,0", "ME,1,17,41"]
    ui = run(lines)
    m = ui.model
    P0 = gcoords(m)
    key = {i: tuple(np.round(p, 9)) for i, p in P0.items()}
    geo0 = _element_geometry(m)
    me_set0 = {key[i] for i in m.nodes if 17 <= i <= 41}
    ui.execute("NCOM")
    assert sorted(m.nodes) == list(range(1, 12))
    assert _element_geometry(m) == geo0
    P1 = {i: tuple(np.round(p, 9)) for i, p in gcoords(m).items()}
    inv = {v: k for k, v in P1.items()}
    assert set(m.tmass) == {inv[key[90]]} and set(m.rmass) == {inv[key[201]]}
    assert m.nodes[inv[key[17]]].fix[0] == 1 and int_set(m) == {inv[key[41]]}
    assert set(m.forces) == {inv[key[77]]}
    assert m.nout[0].nodes == [inv[key[18]]] and m.rdnd[0].node == inv[key[40]]
    rec = m.options.entries("ME")[0][1]
    lo, hi = rec.integer(2), rec.integer(3)
    assert {P1[i] for i in m.nodes if lo <= i <= hi} == me_set0


def test_weld_grid_of_separate_cubes_is_connected_and_geometry_invariant():
    lines, first, e = [], 1, 0
    pts = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    lines.append("GROUP,1,SOLID")
    nodes = []
    for i in range(2):
        for j in range(2):
            ids = list(range(first, first + 8))
            nodes += [f"N,{n},{i + p[0]},{j + p[1]},{p[2]}" for n, p in zip(ids, pts)]
            e += 1
            lines.append("E," + ",".join(str(x) for x in [e] + ids))
            first += 8
    ui = run(nodes + lines)
    m = ui.model
    geo0 = _element_geometry(m)
    ui.execute("WELD")
    assert _element_geometry(m) == geo0
    used = m.used_nodes()
    assert len(used) == 18                                     # 3 x 3 x 2 grid
    P = {tuple(np.round(m.node_global(n), 9)) for n in used}
    assert len(P) == len(used)                                 # no two used nodes coincide
    # the kept number is the lowest of every coincident cluster
    allp = {}
    for n in m.nodes:
        allp.setdefault(tuple(np.round(m.node_global(n), 9)), []).append(n)
    assert set(used) == {min(v) for v in allp.values()}


def test_merge_offsets_are_table_maxima_not_counts():
    m1 = ["N,2,0,0,0", "N,50,1,0,0", "N,7,0,1,0", "M,3,4e6,0.2,0.15,0.05,0.05", "GROUP,9,BEAMS", "MACT,3",
          "E,4,2,50,7", "R,6,1,1,1,1,1,1"]
    m2 = ["N,1,0,0,0", "N,2,1,0,0", "N,3,0,1,0", "M,1,1e6,0.2,0.15,0.05,0.05", "R,1,2,2,2,2,2,2",
          "GROUP,1,BEAMS", "MACT,1", "RACT,1", "E,1,1,2,3"]
    ui = run(["ACTM,1"] + m1 + ["ACTM,2"] + m2 + ["ACTM,3", "MERGE,1,2,0,0,0"])
    m = ui.model
    assert sorted(m.nodes) == [2, 7, 50, 51, 52, 53]
    assert sorted(m.groups) == [9, 10]
    e = m.groups[10].elements[1]
    assert e.nodes == [51, 52, 53] and e.mat == 4 and e.prop == 7
    assert sorted(m.materials) == [3, 4] and sorted(m.sections) == [6, 7]


def _mergesoil(mode_args, floors=False):
    ui = run(["ACTM,1"] + basement_lines(2, 5.0, (-10.0, -5.0, 0.0), floors=floors) + SOIL3
             + ["EXCAV,2", "ACTM,2", "ETYPEGEN,2", "ACTM,3", f"MERGESOIL,1,2,{mode_args}"])
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    return ui


def test_mergesoil_mode1_keeps_soil_geometry_and_only_shares_interface_nodes():
    ui = _mergesoil("1", floors=True)
    m, soil = ui.model, ui.models[2]
    struct = ui.models[1]
    sg = [g for g in m.groups.values() if g.type == 1][0]
    # every soil element keeps its geometry (soil numbers offset by the structure maximum)
    off = max(struct.nodes)
    for e in sg.elements.values():
        e0 = soil.groups[1].elements[e.id]
        assert np.allclose([m.node_global(n) for n in e.nodes], [soil.node_global(n) for n in e0.nodes],
                           rtol=0, atol=1e-9)
    shared = {n for e in sg.elements.values() for n in e.nodes} & set(struct.nodes)
    # shared nodes are exactly the nodes of the lateral and bottom faces (25 of 27)
    v = gl.ModelView(m)
    fsin = v.excavation()["fsin"]
    assert shared == fsin and len(shared) == 25
    assert all(n > off for n in v.excavation()["interior"])


def test_excav_volume_numbering_and_mset_from_topl_position():
    xs, ys = [0.0, 3.0, 7.0, 12.0], [0.0, 4.0, 9.0]
    lines, nid = ["M,1,4.32e6,0.2,0.15,0.05,0.05"], 0
    ids = {}
    for z in (-12.0, -5.0, 0.0):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                if z == -12.0 or i in (0, 3) or j in (0, 2):
                    nid += 1
                    ids[(i, j, z)] = nid
                    lines.append(f"N,{nid},{x},{y},{z}")
    lines.append("GROUP,1,SHELL")
    e = 0
    for j in range(2):
        for i in range(3):
            e += 1
            lines.append(f"E,{e},{ids[(i, j, -12.0)]},{ids[(i + 1, j, -12.0)]},{ids[(i + 1, j + 1, -12.0)]},"
                         f"{ids[(i, j + 1, -12.0)]}")
    loop = [(i, 0) for i in range(3)] + [(3, j) for j in range(2)] + [(i, 2) for i in range(3, 0, -1)] + \
           [(0, j) for j in range(2, 0, -1)]
    for z0, z1 in ((-12.0, -5.0), (-5.0, 0.0)):
        for a in range(len(loop)):
            p, q = loop[a], loop[(a + 1) % len(loop)]
            e += 1
            lines.append(f"E,{e},{ids[(p[0], p[1], z0)]},{ids[(q[0], q[1], z0)]},{ids[(q[0], q[1], z1)]},"
                         f"{ids[(p[0], p[1], z1)]}")
    lines += [f"THICK,1,{e},1,0.5", f"ETYPE,1,{e},1,1", "GROUNDELEV,0",
              "L,7,5,0.12,1500,800,0.05,0.05", "L,4,7,0.12,1600,850,0.05,0.05", "L,9,20,0.13,2000,1000,0.02,0.02",
              "TOPL,7,4,9", "EXCAV,5"]
    ui = run(lines)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    ex = ui.models[5]
    els = list(ex.groups[1].elements.values())
    assert len(els) == 12 and len(ex.nodes) == 36
    vol = 0.0
    for el in els:
        P = np.array([ex.node_global(n) for n in el.nodes])
        vol += gl.hull_area(P[:4, :2]) * (P[4:, 2].mean() - P[:4, 2].mean())
        zmid = 0.5 * (P[:4, 2].mean() + P[4:, 2].mean())
        assert el.mat == (7 if zmid > -5.0 else 4) and el.etype == 2
        # positive orientation: bottom face counter-clockwise seen from above
        assert gl.signed_area(P[:4, :2]) > 0
    assert math.isclose(vol, 12.0 * 9.0 * 12.0, rel_tol=1e-12)
    z = [ex.node_global(n)[2] for n in sorted(ex.nodes)]
    assert all(b >= a for a, b in zip(z, z[1:]))                    # bottom-up numbering


# ======================================================================================
# Defects found by the review (fixed; formerly strict xfail)
# ======================================================================================
def test_mergesoil_converts_materials_of_an_etypegen2_soil_model():
    cube = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    soil = [f"N,{k + 1},{5 * p[0]},{5 * p[1]},{5 * p[2] - 5}" for k, p in enumerate(cube)]
    soil += ["M,1,9e5,0.3,0.11,0.04,0.03", "GROUP,1,SOLID", "MACT,1", "E,1,1,2,3,4,5,6,7,8", "ETYPEGEN,2",
             "GROUNDELEV,0"]
    ui = run(["ACTM,1"] + basement_lines(1, 5.0, (-5.0, 0.0)) + ["ETYPEGEN,1", "ACTM,2"] + soil
             + ["ACTM,3", "MERGESOIL,1,2,1"])
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    m = ui.model
    assert 1 in m.layers and math.isclose(m.layers[1].weight, 0.11)


def test_excav_with_delta_then_mergesoil_joins_the_interface():
    ui = run(["ACTM,1"] + basement_lines(2, 5.0, (-10.0, -5.0, 0.0), jitter=0.001) + SOIL3
             + ["EXCAV,2,0.01", "ACTM,2", "ETYPEGEN,2", "ACTM,3", "MERGESOIL,1,2,1"])
    m, struct = ui.model, ui.models[1]
    sg = [g for g in m.groups.values() if g.type == 1][0]
    shared = {n for e in sg.elements.values() for n in e.nodes} & set(struct.nodes)
    assert len(shared) == 25


def test_mergesoil_warns_when_interface_nodes_find_no_partner():
    ui = run(["ACTM,1"] + basement_lines(2, 5.0, (-10.0, -5.0, 0.0)) + SOIL3
             + ["EXCAV,2", "TRANSLATE,0.01,0,0", "ACTM,2", "ETYPEGEN,2", "ACTM,3", "MERGESOIL,1,2,1"])
    assert any("interface" in w.lower() for w in ui.sink.texts(Kind.WARNING)), ui.sink.texts(Kind.WARNING)


def test_weld_still_welds_other_nodes_coincident_with_a_protected_spring():
    ui = run(["N,1,0,0,0", "N,2,1,0,0", "N,3,1,1,0", "N,4,0,1,0",
              "N,5,1,0,0", "N,6,2,0,0", "N,7,2,1,0", "N,8,1,1,0", "N,9,1,0,0",
              "M,1,4e6,0.2,0.15,0.05,0.05", "GROUP,1,SHELL", "E,1,1,2,3,4", "E,2,5,6,7,8", "THICK,1,2,1,0.5",
              "SC,1,0,0,0,10,10,10,0", "GROUP,2,SPRING", "RACT,1", "E,1,2,9", "D,9,9,1,1,ALL", "WELD"])
    m = ui.model
    assert m.groups[2].elements[1].nodes == [2, 9]             # the spring is kept
    assert m.groups[1].elements[2].nodes == [2, 6, 7, 3]       # plate B welded to plate A at both nodes


def test_mergesoil_mode3_seplevel_uses_the_geometric_tolerance():
    ui = run(["ACTM,1"] + basement_lines(2, 5.0, (-10.0, -5.0 + 1e-12, 0.0)) + SOIL3
             + ["EXCAV,2", "ACTM,2", "ETYPEGEN,2", "ACTM,3", "MERGESOIL,1,2,3,,,-5"])
    m = ui.model
    sg = [g for g in m.groups.values() if g.type == 7][0]
    k = [m.springs[e.prop].scx for e in sg.elements.values()]
    assert sum(1 for x in k if x == 1e7) == 9 + 8 and sum(1 for x in k if x == 10.0) == 8


def test_mergesoil_mapping_file_in_the_mode_slot(tmp_path):
    ui = run(["ACTM,1"] + basement_lines(2, 5.0, (-10.0, -5.0, 0.0)) + SOIL3
             + ["EXCAV,2", "ACTM,2", "ETYPEGEN,2", "ACTM,3", "MERGESOIL,1,2,excv.map"], Interpreter(cwd=tmp_path))
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    assert (tmp_path / "excv.map").exists()


def test_rotate_warns_about_global_axis_fixities():
    from sassi.verify.problems.vp_generation import plate_lines
    ui = run(plate_lines(2, 2) + ["FIXROT", "ROTATE,0,0,0,0,30,0"])
    assert any("fix" in w.lower() for w in ui.sink.texts(Kind.WARNING)), ui.sink.texts(Kind.WARNING)
