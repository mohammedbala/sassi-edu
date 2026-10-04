"""Adversarial review of the cuts / submodel / section-cut package (requirements 3.4.J and 4.14; spec 09
section 5; spec 10 sections 1.4 and 3; decisions D-MDL-13, D-SEC-01 ... D-SEC-06).

The references are derived independently of :mod:`sassi.prep.cuts_lib`:

* SOLID sections of a non-uniform box mesh are rebuilt by **half-space clipping** of a large square in
  the cut plane against the six faces of every box (Sutherland-Hodgman) and integrated with the exact
  triangle formulas ``int x x^T dA = A/12 (sum v v^T + s s^T)`` (s = sum of the vertices), instead of the
  implementation's edge crossings + angular sort + shoelace formulas;
* SHELL stresses are produced by the element recovery itself (:func:`sassi.elements.shell.recovery`,
  the operator STRESS applies) from imposed membrane-strain and curvature fields, and the expected
  resultant is written with frame-free tensors (plane-stress law on the projected strain tensor,
  ``M = -D[(1-nu) K + nu tr(K) P]`` for the Hessian K of w), so any mismatch of the local x' y' axes
  or the component order between STRESS and the section cut shows up;
* invariances: Newton's third law (reversed normal), translation of everything, moment transfer between
  sub-sections, the closed-form discretisation error ``M/M_exact = 1 - 1/n^2`` of a linearly varying
  stress sampled at n strip centres, and CSECT element validity (det J > 0) on random planes.

Tests marked ``xfail(strict=True)`` document defects found in the review (the reason names the
defect); they turn into XPASS failures once the defect is fixed, so the marker must then be removed.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pytest

from sassi.elements import shell as SH
from sassi.elements.base import material_from_M
from sassi.model.geometry import loc_matrix
from sassi.prep import Interpreter, Kind
from sassi.prep import cuts_lib as cl
from sassi.prep.writer import write_pre

NAMES = ("Area", "Xc", "Yc", "Zc", "Ixx", "Iyy", "Ixy", "Fx", "Fy", "Fz", "Mx", "My", "Mz")


# ======================================================================================
# helpers (independent of the implementation)
# ======================================================================================
def num(x) -> str:
    return repr(float(x))


def vec(v) -> str:
    return ",".join(num(x) for x in v)


def run(ui: Interpreter, lines: Sequence[str]) -> None:
    for ln in lines:
        assert ui.execute(ln), (ln, ui.sink.texts(Kind.ERROR)[-1:])


def axes(n, r) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ez = np.asarray(n, float) / np.linalg.norm(n)
    ex = np.asarray(r, float) - (np.asarray(r, float) @ ez) * ez
    ex /= np.linalg.norm(ex)
    return ex, np.cross(ez, ex), ez


def box_mesh(xs, ys, zs, group: int = 1):
    """Hexahedra on the tensor grid xs x ys x zs; returns (.pre lines, {element: (lo, hi)})."""
    lines, nid, boxes = [], {}, {}
    k = 0
    for iz, z in enumerate(zs):
        for iy, y in enumerate(ys):
            for ix, x in enumerate(xs):
                k += 1
                nid[(ix, iy, iz)] = k
                lines.append(f"N,{k},{num(x)},{num(y)},{num(z)}")
    lines.append(f"GROUP,{group},SOLID")
    e = 0
    for iz in range(len(zs) - 1):
        for iy in range(len(ys) - 1):
            for ix in range(len(xs) - 1):
                e += 1
                c = [(ix, iy, iz), (ix + 1, iy, iz), (ix + 1, iy + 1, iz), (ix, iy + 1, iz)]
                ids = [nid[q] for q in c] + [nid[(a, b, cc + 1)] for a, b, cc in c]
                lines.append(f"E,{e}," + ",".join(map(str, ids)))
                boxes[e] = (np.array([xs[ix], ys[iy], zs[iz]], float), np.array([xs[ix + 1], ys[iy + 1], zs[iz + 1]],
                                                                                  float))
    return lines, boxes


def _clip(poly: List[np.ndarray], a: np.ndarray, b: float) -> List[np.ndarray]:
    """Sutherland-Hodgman: keep a.x <= b."""
    out = []
    for i in range(len(poly)):
        p, q = poly[i], poly[(i + 1) % len(poly)]
        fp, fq = float(a @ p) - b, float(a @ q) - b
        if fp <= 0:
            out.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            out.append(p + fp / (fp - fq) * (q - p))
    return out


def clipping_reference(boxes, stresses, P, n, r, keep=None) -> Dict[str, float]:
    """CALCPAR quantities of the plane (P, n) through the box mesh by half-space clipping (independent)."""
    ex, ey, ez = axes(n, r)
    P = np.asarray(P, float)
    sq = [P + 1e3 * (sx * ex + sy * ey) for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    A, S1, S2, data = 0.0, np.zeros(2), np.zeros((2, 2)), []
    for e, (lo, hi) in boxes.items():
        if keep is not None and not keep(e, lo, hi):
            continue
        poly = sq
        for k in range(3):
            a = np.zeros(3)
            a[k] = 1.0
            poly = _clip(poly, a, hi[k])
            if poly:
                poly = _clip(poly, -a, -lo[k])
            if not poly:
                break
        if len(poly) < 3:
            continue
        xy = [np.array([(p - P) @ ex, (p - P) @ ey]) for p in poly]
        Ae, s1, s2 = 0.0, np.zeros(2), np.zeros((2, 2))
        for i in range(1, len(xy) - 1):
            a, b, c = xy[0], xy[i], xy[i + 1]
            At = 0.5 * abs((b - a)[0] * (c - a)[1] - (b - a)[1] * (c - a)[0])
            s = a + b + c
            Ae += At
            s1 += At * s / 3.0
            s2 += At / 12.0 * (np.outer(a, a) + np.outer(b, b) + np.outer(c, c) + np.outer(s, s))
        if Ae < 1e-12:
            continue
        A += Ae
        S1 += s1
        S2 += s2
        data.append((e, Ae, s1 / Ae))
    c2 = S1 / A
    J = S2 - A * np.outer(c2, c2)
    C = P + c2[0] * ex + c2[1] * ey
    F, M = np.zeros(3), np.zeros(3)
    for e, Ae, ce in data:
        s = stresses[e]
        T = np.array([[s[0], s[3], s[4]], [s[3], s[1], s[5]], [s[4], s[5], s[2]]])
        f = Ae * T @ ez
        F += f
        M += np.cross(P + ce[0] * ex + ce[1] * ey - C, f)
    E = np.vstack([ex, ey, ez])
    Fl, Ml = E @ F, E @ M
    return dict(Area=A, Xc=C[0], Yc=C[1], Zc=C[2], Ixx=J[1, 1], Iyy=J[0, 0], Ixy=J[0, 1],
                Fx=Fl[0], Fy=Fl[1], Fz=Fl[2], Mx=Ml[0], My=Ml[1], Mz=Ml[2])


def read_csv(path: Path) -> Tuple[List[str], List[List[float]]]:
    rows = [ln.split(",") for ln in path.read_text().splitlines()]
    return rows[0], [[float(x) for x in r[1:]] for r in rows[1:]]


def assert_close(got: Dict[str, float], ref: Dict[str, float], rtol: float, keys=NAMES) -> None:
    scale = max(1.0, max(abs(ref[k]) for k in keys))
    bad = {k: (got[k], ref[k]) for k in keys if abs(got[k] - ref[k]) > rtol * scale}
    assert not bad, bad


XS, YS, ZS = [0.0, 0.7, 1.5, 2.2], [0.0, 1.3, 2.1], [0.0, 0.9, 1.6, 2.6, 3.1]


@pytest.fixture
def mesh(tmp_path):
    lines, boxes = box_mesh(XS, YS, ZS)
    rng = np.random.default_rng(20261001)
    st = {e: rng.normal(size=6) * 10.0 for e in boxes}
    cl.write_ess(tmp_path / "s.ess", {1: "SOLID"}, {(1, e): v for e, v in st.items()})
    (tmp_path / "s.lst").write_text("frames\ns.ess\n")
    return lines, boxes, st


# ======================================================================================
# 1. SOLID sections against an independent clipping integration
# ======================================================================================
OBLIQUE = [([1.1, 1.0, 1.7], [0.3, -0.5, 0.81], [1.0, 0.2, 0.0]),
           ([0.7, 1.3, 1.6], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]),          # through mesh nodes
           ([1.1, 1.05, 1.55], [0.05, -0.02, 1.0], [0.0, 1.0, 0.0])]


@pytest.mark.parametrize("P,n,r", OBLIQUE)
def test_oblique_solid_section_matches_halfspace_clipping(tmp_path, mesh, P, n, r):
    """READSTR -> CUTVOL -> CSECT -> ACTM -> CALCPAR and CALCSECTHIST on the original model both equal the
    clipping reference (area, centroid, inertias and the six resultants)."""
    lines, boxes, st = mesh
    ref = clipping_reference(boxes, st, P, n, r)
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["READSTR,s.ess", "CUTVOL,1", f"CSECT,5,1,{vec(P)},{vec(n)}", "ACTM,5", f"CALCPAR,{vec(n)},{vec(r)}"])
    assert_close(ui.session["calcpar"], ref, 1e-10)
    run(ui, ["ACTM,0", f"CALCSECTHIST,s.lst,1,{vec(P)},{vec(n)},{vec(r)},0,0,h.csv"])
    head, rows = read_csv(tmp_path / "h.csv")
    got = dict(zip(head[1:], rows[0]))
    assert_close({**ref, **got}, ref, 2e-6, keys=("Fx", "Fy", "Fz", "Mx", "My", "Mz"))     # %.6E output


def test_plane_on_a_mesh_level_counts_the_minus_side_element_once(tmp_path, mesh):
    """D-SEC-04: the plane z = 1.6 lies on a row of shared faces; the result equals the clipping reference
    of the elements just below the plane (the -n side) and, with the normal reversed, of those above."""
    lines, boxes, st = mesh
    P, r = [1.0, 1.0, 1.6], [1.0, 0.0, 0.0]
    below = clipping_reference(boxes, st, [1.0, 1.0, 1.6 - 1e-9], [0, 0, 1], r,
                               keep=lambda e, lo, hi: abs(hi[2] - 1.6) < 1e-12)
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["CUTVOL,1"])
    run(ui, [f"CALCSECTHIST,s.lst,1,{vec(P)},0,0,1,{vec(r)},0,0,h.csv"])
    head, rows = read_csv(tmp_path / "h.csv")
    got = dict(zip(head[1:], rows[0]))
    for k in ("Fx", "Fy", "Fz"):
        assert got[k] == pytest.approx(below[k], rel=2e-6, abs=1e-5), k
    # moments about the same centroid (the section geometry is the same from both sides)
    for k in ("Mx", "My", "Mz"):
        assert got[k] == pytest.approx(below[k], rel=2e-6, abs=1e-5), k
    above = clipping_reference(boxes, st, [1.0, 1.0, 1.6 + 1e-9], [0, 0, -1], r,
                               keep=lambda e, lo, hi: abs(lo[2] - 1.6) < 1e-12)
    run(ui, [f"CALCSECTHIST,s.lst,1,{vec(P)},0,0,-1,{vec(r)},0,0,g.csv"])
    head, rows = read_csv(tmp_path / "g.csv")
    got = dict(zip(head[1:], rows[0]))
    for k in ("Fx", "Fy", "Fz", "Mx", "My", "Mz"):
        assert got[k] == pytest.approx(above[k], rel=2e-6, abs=1e-5), k


# ======================================================================================
# 2. SHELL: contract with the element recovery (what STRESS writes) and exact wall cuts
# ======================================================================================
E_, NU_, T_ = 1000.0, 0.25, 0.3
Q_PLATE = loc_matrix(25.0, -35.0, 50.0)
O_PLATE = np.array([1.0, -2.0, 0.5])
G_STRAIN = np.array([[0.3, -0.2, 0.1], [0.4, -0.1, 0.25], [-0.15, 0.2, 0.05]]) * 1e-2
K_HESS = np.array([[0.2, 0.05, -0.1], [0.05, -0.3, 0.07], [-0.1, 0.07, 0.15]])


def _plate(tmp_path):
    """3 x 2 plate of 6 shells in a rotated plane; element 2 has a rotated node order (other local x'),
    element 5 a reversed one (normal -m).  Element stresses come from the shell recovery operator applied
    to the nodal displacements of a uniform membrane strain plus a uniform curvature field."""
    mat = material_from_M(1, E_, NU_, 0.0, 0.0, 0.0, 9.81)
    pts, lines = {}, []
    na, nb = 4, 3
    for j in range(nb):
        for i in range(na):
            k = j * na + i + 1
            pts[k] = O_PLATE + Q_PLATE @ np.array([float(i), float(j), 0.0])
            lines.append(f"N,{k},{vec(pts[k])}")
    lines.append("GROUP,1,SHELL")
    m = Q_PLATE[:, 2]
    Pp = np.eye(3) - np.outer(m, m)
    Kh = Pp @ K_HESS @ Pp

    def disp(p):
        x = p - O_PLATE
        return np.concatenate([G_STRAIN @ x + 0.5 * (x @ Kh @ x) * m, np.cross(Kh @ x, m)])

    st, e = {}, 0
    for j in range(nb - 1):
        for i in range(na - 1):
            e += 1
            n1 = j * na + i + 1
            ids = [n1, n1 + 1, n1 + 1 + na, n1 + na]
            if e == 2:
                ids = ids[1:] + ids[:1]
            if e == 5:
                ids = ids[::-1]
            lines.append(f"E,{e}," + ",".join(map(str, ids)))
            X = np.array([pts[q] for q in ids])
            u = np.concatenate([disp(pts[q]) for q in ids])
            st[(1, e)] = np.real(SH.recovery(X, mat, thick=T_) @ u)
    lines.append(f"THICK,1,{e},1,{num(T_)}")
    cl.write_ess(tmp_path / "p.ess", {1: "SHELL"}, st)
    (tmp_path / "p.lst").write_text("# frames\np.ess\n")
    eps = Pp @ (0.5 * (G_STRAIN + G_STRAIN.T)) @ Pp
    sig = E_ / (1 - NU_ ** 2) * ((1 - NU_) * eps + NU_ * np.trace(eps) * Pp)
    D = E_ * T_ ** 3 / (12 * (1 - NU_ ** 2))
    Mt = -D * ((1 - NU_) * Kh + NU_ * np.trace(Kh) * Pp)
    return lines, m, sig, Mt


@pytest.mark.parametrize("p_ab,n_loc", [((1.5, 1.0), (1.0, 0.2, 0.3)), ((1.3, 1.1), (0.4, 1.0, -0.6)),
                                        ((2.2, 0.7), (1.0, -1.0, 0.0))])
def test_shell_section_matches_the_stress_recovery_contract(tmp_path, p_ab, n_loc):
    """The membrane traction ``L t sigma.nu_s`` and the couple ``L m x (M.nu_s)`` (spec 10 section 3.2)
    computed with frame-free tensors equal CALCPAR on the CSECT model and CALCSECTHIST on the original
    model, for shells whose local axes differ element by element (rotated / reversed node order)."""
    lines, m, sig, Mt = _plate(tmp_path)
    P = O_PLATE + Q_PLATE @ np.array([p_ab[0], p_ab[1], 0.0])
    n = Q_PLATE @ np.array(n_loc, float)
    r = [0.0, 0.0, 1.0]
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["READSTR,p.ess", f"SLICE,1,{vec(P)},{vec(n)}", f"CSECT,4,1,{vec(P)},{vec(n)}", "ACTM,4",
                     f"CALCPAR,{vec(n)},{vec(r)}"])
    res = ui.session["calcpar"]
    ex, ey, ez = axes(n, r)
    E = np.vstack([ex, ey, ez])
    nu_s = ez - (ez @ m) * m
    nu_s /= np.linalg.norm(nu_s)
    d = np.cross(m, nu_s)
    # chord of the cut line across the 3 x 2 plate rectangle (plate coordinates a, b)
    p0, dl = np.array(p_ab, float), (Q_PLATE.T @ d)[:2]
    ts = sorted((bd - p0[k]) / dl[k] for k, lim in enumerate((3.0, 2.0)) for bd in (0.0, lim) if abs(dl[k]) > 1e-12)
    t0, t1 = max(x for x in ts if x <= 0), min(x for x in ts if x >= 0)
    L = t1 - t0
    F = L * T_ * sig @ nu_s
    C = np.array([res["Xc"], res["Yc"], res["Zc"]])
    assert np.allclose(C, P + 0.5 * (t0 + t1) * d, atol=1e-12)          # the strip midpoint
    assert res["Area"] == pytest.approx(L * T_ / np.linalg.norm(np.cross(m, ez)), rel=1e-12)
    M = L * np.cross(m, Mt @ nu_s)
    assert np.allclose([res["Fx"], res["Fy"], res["Fz"]], E @ F, atol=1e-10)
    assert np.allclose([res["Mx"], res["My"], res["Mz"]], E @ M, atol=1e-10)
    run(ui, ["ACTM,0", f"CALCSECTHIST,p.lst,1,{vec(P)},{vec(n)},{vec(r)},0,0,h.csv"])
    head, rows = read_csv(tmp_path / "h.csv")
    got = np.array(rows[0])
    assert np.allclose(got, np.concatenate([E @ F, E @ M]), rtol=2e-6, atol=1e-7)


@pytest.mark.parametrize("n", [(0, 0, 1), (0, 0.6, 0.8), (0.3, 0.4, 1.0), (-0.25, -0.5, 1.0)])
def test_doubly_oblique_wall_cut_transmits_the_wall_force_and_moment(tmp_path, n):
    """Wall y = 0 (x in [0, 4], z in [0, 2], t = 0.5, free vertical edges) with uniform Sy'y' = s and
    My'y' = M (y' = Z, z' = -Y): the field is equilibrated and traction-free on the free edges, so every
    plane crossing the whole wall transmits F = 4 t s Z and M = 4 M (z' x Z) = -4 M X about the centroid,
    whatever its obliquity."""
    t, s, M = 0.5, 37.0, 11.0
    lines = ["N,1,0,0,0", "N,2,2,0,0", "N,3,4,0,0", "N,4,0,0,2", "N,5,2,0,2", "N,6,4,0,2",
             "GROUP,1,SHELL", "E,1,1,2,5,4", "E,2,2,3,6,5", f"THICK,1,2,1,{t}"]
    cl.write_ess(tmp_path / "w.ess", {1: "SHELL"}, {(1, 1): [0, s, 0, 0, M, 0], (1, 2): [0, s, 0, 0, M, 0]})
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["READSTR,w.ess", f"SLICE,1,2,0,1,{vec(n)}", f"CSECT,4,1,2,0,1,{vec(n)}", "ACTM,4",
                     f"CALCPAR,{vec(n)},1,0,0"])
    r = ui.session["calcpar"]
    ex, ey, ez = axes(n, [1, 0, 0])
    E = np.vstack([ex, ey, ez])
    assert np.allclose(E.T @ [r["Fx"], r["Fy"], r["Fz"]], [0, 0, 4 * t * s], atol=1e-10)
    assert np.allclose(E.T @ [r["Mx"], r["My"], r["Mz"]], [-4 * M, 0, 0], atol=1e-10)
    assert [r["Xc"], r["Yc"], r["Zc"]] == pytest.approx([2.0, 0.0, 1.0], abs=1e-12)


# ======================================================================================
# 3. Invariances and closed-form discretisation error
# ======================================================================================
def test_reversed_normal_is_newtons_third_law(tmp_path, mesh):
    """With -n and the same r: ez' = -ez, ex' = ex, ey' = -ey and the force is the action of the other side
    (F' = -F, M' = -M about the same centroid), so (Fx, Fy, Fz, Mx, My, Mz)' = (-Fx, Fy, Fz, -Mx, My, Mz)."""
    lines, boxes, st = mesh
    ui = Interpreter(cwd=tmp_path)
    P, n, r = [1.1, 1.0, 1.7], np.array([0.3, -0.5, 0.81]), [1.0, 0.2, 0.0]
    run(ui, lines + ["CUTVOL,1", f"CALCSECTHIST,s.lst,1,{vec(P)},{vec(n)},{vec(r)},0,0,a.csv",
                     f"CALCSECTHIST,s.lst,1,{vec(P)},{vec(-n)},{vec(r)},0,0,b.csv"])
    a = np.array(read_csv(tmp_path / "a.csv")[1][0])
    b = np.array(read_csv(tmp_path / "b.csv")[1][0])
    assert np.allclose(b, a * [-1, 1, 1, -1, 1, 1], rtol=1e-6, atol=1e-6)


def test_moment_transfer_between_sub_sections(tmp_path, mesh):
    """Two disjoint sub-cuts of one section: the forces add, and the moments add once transferred to the
    common centroid, M = M1 + (C1 - C) x F1 + M2 + (C2 - C) x F2 (global)."""
    lines, boxes, st = mesh
    P, n, r = np.array([1.1, 1.0, 1.7]), np.array([0.3, -0.5, 0.81]), np.array([1.0, 0.2, 0.0])
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["READSTR,s.ess"])
    left = [e for e, (lo, hi) in boxes.items() if hi[0] <= 1.5 + 1e-12]
    right = [e for e in boxes if e not in left]
    run(ui, ["CUTADD,1,1," + ",".join(map(str, left)), "CUTADD,2,1," + ",".join(map(str, right)), "CUTVOL,3"])
    E = np.vstack(axes(n, r))
    out = {}
    for cut in (1, 2, 3):
        run(ui, ["ACTM,0", f"CSECT,{10 + cut},{cut},{vec(P)},{vec(n)}", f"ACTM,{10 + cut}", f"CALCPAR,{vec(n)},{vec(r)}"])
        q = ui.session["calcpar"]
        out[cut] = (np.array([q["Xc"], q["Yc"], q["Zc"]]), E.T @ [q["Fx"], q["Fy"], q["Fz"]],
                    E.T @ [q["Mx"], q["My"], q["Mz"]], q["Area"])
    C, F, M, A = out[3]
    (C1, F1, M1, A1), (C2, F2, M2, A2) = out[1], out[2]
    assert A1 + A2 == pytest.approx(A, rel=1e-12)
    assert np.allclose((A1 * C1 + A2 * C2) / A, C, atol=1e-12)
    assert np.allclose(F1 + F2, F, atol=1e-10)
    assert np.allclose(M1 + np.cross(C1 - C, F1) + M2 + np.cross(C2 - C, F2), M, atol=1e-10)


@pytest.mark.parametrize("nstrip", [1, 2, 4, 8])
def test_linear_stress_sampled_at_strip_centres(tmp_path, nstrip):
    """Spec 04 section 12 item 6: a bending stress Szz = k x' across a b x h rectangle of n strips, sampled
    at the strip centres, gives M = k sum A_i x_i^2 = (k I)(1 - 1/n^2) exactly and Fz = 0; the inertia is
    the exact b^3 h / 12."""
    b, h, k = 3.0, 0.8, 50.0
    xs = list(np.linspace(0.0, b, nstrip + 1))
    lines, boxes = box_mesh(xs, [0.0, h], [0.0, 1.0])
    st = {e: np.array([0, 0, k * (0.5 * (lo[0] + hi[0]) - b / 2), 0, 0, 0]) for e, (lo, hi) in boxes.items()}
    cl.write_ess(tmp_path / "b.ess", {1: "SOLID"}, {(1, e): v for e, v in st.items()})
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["READSTR,b.ess", "SLICE,1,0,0,0.5,0,0,1", "CSECT,2,1,0,0,0.5,0,0,1", "ACTM,2",
                     "CALCPAR,0,0,1,1,0,0"])
    r = ui.session["calcpar"]
    I = h * b ** 3 / 12.0
    assert r["Iyy"] == pytest.approx(I, rel=1e-12)
    assert r["Fz"] == pytest.approx(0.0, abs=1e-9)
    # +x fibres in tension -> moment about +Y is -(int x sigma dA): (x X) x (sigma Z) = -x sigma Y
    assert r["My"] == pytest.approx(-k * I * (1.0 - 1.0 / nstrip ** 2), rel=1e-12, abs=1e-12)


def test_translating_everything_moves_only_the_centroid(tmp_path, mesh):
    lines, boxes, st = mesh
    P, n, r = np.array([1.1, 1.0, 1.7]), np.array([0.3, -0.5, 0.81]), [1.0, 0.2, 0.0]
    shift = np.array([123.0, -45.5, 6.25])
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["READSTR,s.ess", "CUTVOL,1", f"CSECT,5,1,{vec(P)},{vec(n)}", "ACTM,5", f"CALCPAR,{vec(n)},{vec(r)}"])
    a = dict(ui.session["calcpar"])
    ui2 = Interpreter(cwd=tmp_path)
    run(ui2, lines + ["READSTR,s.ess", f"TRANSLATE,{vec(shift)}", "CUTVOL,1", f"CSECT,5,1,{vec(P + shift)},{vec(n)}",
                      "ACTM,5", f"CALCPAR,{vec(n)},{vec(r)}"])
    b = ui2.session["calcpar"]
    for k in NAMES:
        ref = a[k] + (shift["XYZ".index(k[0])] if k in ("Xc", "Yc", "Zc") else 0.0)
        assert b[k] == pytest.approx(ref, rel=1e-9, abs=1e-9), k


# ======================================================================================
# 4. Stored cut system, CSECT element validity, SAVE / RESUME
# ======================================================================================
@pytest.mark.parametrize("n,r", [((1, 0, 0), (0, 1, 0)), ((0, 1, 0), (0, 0, 1)), ((1, 1, 1), (1, -1, 0)),
                                 ((0, 0, -1), (1, 0, 0)), ((0.3, -0.4, 0.8), (0.2, 0.9, -0.1)),
                                 ((0, -1, 0), (1, 0, 0))])
def test_stored_cut_system_maps_local_coordinates(tmp_path, n, r):
    """CALCPAR stores LOC <sysno> with origin C and axes (ex, ey, ez) of D-SEC-01 (incl. gimbal-lock
    orientations): a node given in that system at local (1, 2, 3) lies at C + ex + 2 ey + 3 ez, also after
    WRITE -> INP; the active system is unchanged.

    (Adapted by the implementer: CALCMOI now runs on a CSECT model through the box centre.  It used to run on
    the original box, whose section was the mid-extent plane; after the review's defect 6 / Q-11 a model that
    is not a cross-section gives its base section, which an oblique n only touches at a corner.  The subject
    of the test, the stored system, is unchanged.)"""
    lines, boxes = box_mesh([0.0, 1.0, 2.0], [0.0, 1.0], [0.0, 1.0])
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["CUTVOL,1", f"CSECT,2,1,1,0.5,0.5,{vec(n)}", "ACTM,2", f"CALCMOI,{vec(n)},{vec(r)},9"])
    q = ui.session["calcmoi"]
    C = np.array([q["Xc"], q["Yc"], q["Zc"]])
    ex, ey, ez = axes(n, r)
    assert ui.model.csys_active == 0
    run(ui, ["CSYS,9", "N,500,1,2,3", "CSYS,0"])
    want = C + ex + 2 * ey + 3 * ez
    assert np.allclose(ui.model.node_global(500), want, atol=1e-12)
    b = Interpreter()
    b.run_text(write_pre(ui.model)[0])
    assert np.allclose(b.model.node_global(500), want, atol=1e-9)


def test_csect_random_planes_give_valid_unit_thickness_solids():
    """CSECT of one hexahedron by 120 random planes (3- to 6-vertex polygons): every extruded sub-element has
    det J > 0 at the 2x2x2 Gauss points, and CALCC volume / centroid of the cross-section model equal the
    CALCMOI area / centroid (spec 10 section 3.3: the extrusion is +-0.5 along n)."""
    rng = np.random.default_rng(5)
    cube = ["N,1,0,0,0", "N,2,1.2,0,0", "N,3,1.3,0.9,0", "N,4,-0.1,1.1,0", "N,5,0,0,1", "N,6,1,0.1,1.1",
            "N,7,1.1,1,1", "N,8,0,1,0.9", "GROUP,1,SOLID", "E,1,1,2,3,4,5,6,7,8", "CUTADD,1,1,1"]
    g = 1.0 / math.sqrt(3.0)
    nat = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])
    gps = [np.array([a, b, c]) * g for a in (-1, 1) for b in (-1, 1) for c in (-1, 1)]
    seen = set()
    for _ in range(120):
        n = rng.normal(size=3)
        P = rng.uniform(0.15, 0.85, size=3)
        ui = Interpreter()
        run(ui, cube + [f"CSECT,2,1,{vec(P)},{vec(n)}", "ACTM,2", f"CALCMOI,{vec(n)},{vec(np.cross(n, [0.3, 0.5, 0.7]))}",
                        "CALCC"])
        s = ui.model
        seen.add(len(s.ui_state["csect"]["pieces"][0]["pts"]))
        mo, cc = ui.session["calcmoi"], ui.session["calcc"]
        assert cc["Volume"] == pytest.approx(mo["Area"], rel=1e-12)
        assert [cc["Xc"], cc["Yc"], cc["Zc"]] == pytest.approx([mo["Xc"], mo["Yc"], mo["Zc"]], abs=1e-12)
        for grp, el in s.iter_elements():
            X = np.array([s.node_global(k) for k in el.nodes])
            for xi in gps:
                dN = np.array([[nat[a, k] * np.prod([1 + nat[a, j] * xi[j] for j in range(3) if j != k]) / 8
                                for a in range(8)] for k in range(3)])
                assert np.linalg.det(dN @ X) > 0, (grp.id, el.id)
    assert {4, 5, 6} <= seen


def test_save_resume_keeps_the_section_data(tmp_path):
    """The CSECT pieces and the copied element stresses live in ui_state: SAVE -> RESUME reproduces CALCPAR."""
    lines = ["N,1,0,0,0", "N,2,2,0,0", "N,3,4,0,0", "N,4,0,0,2", "N,5,2,0,2", "N,6,4,0,2",
             "GROUP,1,SHELL", "E,1,1,2,5,4", "E,2,2,3,6,5", "THICK,1,2,1,0.5"]
    cl.write_ess(tmp_path / "w.ess", {1: "SHELL"}, {(1, 1): [1, 100, 2, 3, 4, 5], (1, 2): [-1, -100, 2, 6, -4, 1]})
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["READSTR,w.ess", "SLICE,1,0,0,1,0,0.3,0.9", "CSECT,4,1,2,0,1,0.2,0.3,0.9", "ACTM,4",
                     f"MDL,cs,{tmp_path}", "CALCPAR,0.2,0.3,0.9,1,0,0", "SAVE"])
    ref = dict(ui.session["calcpar"])
    assert abs(ref["Mx"]) > 1 and abs(ref["Fx"]) > 1
    b = Interpreter(cwd=tmp_path)
    run(b, [f"MDL,cs,{tmp_path}", "RESUME", "CALCPAR,0.2,0.3,0.9,1,0,0"])
    for k, v in ref.items():
        assert b.session["calcpar"][k] == pytest.approx(v, rel=1e-12, abs=1e-12), k


def test_stress_module_frames_with_spring_blocks(tmp_path):
    """STRESS writes SPRING blocks into the .ess frames (sassi.core.stress_lib.CENTER_COLUMNS): READSTR
    accepts them and CALCSECTHIST lists the crossing spring as skipped without changing the solid result."""
    SL = pytest.importorskip("sassi.core.stress_lib")
    lines, boxes = box_mesh([0.0, 1.0, 2.0], [0.0, 1.0], [0.0, 1.0])
    lines += ["N,50,5,5,0", "N,51,5,5,1", "GROUP,2,SPRING", "E,1,50,51", "SC,1,1,1,1,0,0,0,0.05", "RSET,1,1,1,1"]
    vals = np.array([[0, 0, 10, 0, 3, 0], [0, 0, -10, 0, 3, 0]], float)
    blocks = [SL.CenterBlock("SOLID", 1, 1, [1, 2], vals), SL.CenterBlock("SPRING", 2, 1, [1], [[1, 2, 3, 4, 5, 6]])]
    d = tmp_path / "NSTRESS"
    d.mkdir()
    SL.write_element_center(d / "ESTRESS_00001.ess", blocks, fmt="{:.10e}")
    (d / "ESTRESS.lst").write_text("# ESTRESS frames of model m: dt = 0.01 s, frame n is time (n-1) dt\n"
                                   "ESTRESS_00001.ess\n")
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["READSTR,NSTRESS/ESTRESS_00001.ess", "SLICE,1,0,0,0.5,0,0,1",
                     "CALCSECTHIST,NSTRESS/ESTRESS.lst,1,0,0,0.5,0,0,1,1,0,0,0,0.01,h.csv"])
    assert len(cl.stress_table(ui.model)) == 3
    assert any("SPRING" in w and "2/1" in w for w in ui.sink.texts(Kind.WARNING))
    head, rows = read_csv(tmp_path / "h.csv")
    got = dict(zip(head[1:], rows[0]))
    assert got["Fx"] == pytest.approx(6.0) and got["My"] == pytest.approx(10.0) and got["Fz"] == pytest.approx(0.0)


# ======================================================================================
# 5. Defects found by the review (were strict xfail; fixed, markers removed)
# ======================================================================================
#: a column crossing z = 0.5; its K (orientation) node lies off the axis, 3.5 above the plane
FRAME_WITH_BEAM = ["N,20,5,5,0", "N,21,5,5,1", "N,22,6,5,4", "GROUP,2,BEAMS", "E,1,20,21,22", "R,1,0.1,0,0,0,0,0"]


def test_defect_calcpar_on_csect_model_warns_that_crossing_beams_are_excluded(tmp_path):
    lines, boxes = box_mesh([0.0, 1.0, 2.0], [0.0, 1.0], [0.0, 1.0])
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + FRAME_WITH_BEAM + ["CUTADD,1,1,1-2", "CUTADD,1,2,1", "CSECT,3,1,0,0,0.5,0,0,1", "ACTM,3"])
    n0 = len(ui.sink.texts(Kind.WARNING))
    run(ui, ["CALCPAR,0,0,1,1,0,0"])
    assert any("BEAMS" in w and "2/1" in w for w in ui.sink.texts(Kind.WARNING)[n0:])


def test_defect_calcmoi_after_write_inp_of_a_csect_model_with_a_beam(tmp_path):
    lines, boxes = box_mesh([0.0, 1.0, 2.0], [0.0, 1.0], [0.0, 1.0])
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + FRAME_WITH_BEAM + ["CUTADD,1,1,1-2", "CUTADD,1,2,1", "CSECT,3,1,0,0,0.5,0,0,1", "ACTM,3",
                                       "CALCMOI,0,0,1,1,0,0"])
    ref = dict(ui.session["calcmoi"])
    b = Interpreter()
    b.run_text(write_pre(ui.model)[0])
    assert b.execute("CALCMOI,0,0,1,1,0,0"), b.sink.texts(Kind.ERROR)[-1:]
    for k, v in ref.items():
        assert b.session["calcmoi"][k] == pytest.approx(v, abs=1e-12), k


def test_defect_cut_system_must_not_redefine_the_active_system(tmp_path):
    lines, boxes = box_mesh([0.0, 1.0, 2.0], [0.0, 1.0], [0.0, 1.0])
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["LOC,5,0,0,0,0,0,0", "CSYS,5", "CALCMOI,0,0,1,0,1,0,5", "N,100,1,0,0"])
    assert np.allclose(ui.model.node_global(100), [1.0, 0.0, 0.0])


def test_defect_calcc_on_csect_model_with_a_beam_equals_the_area_centroid(tmp_path):
    lines, boxes = box_mesh([0.0, 1.0, 2.0], [0.0, 1.0], [0.0, 1.0])
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + FRAME_WITH_BEAM + ["CUTADD,1,1,1-2", "CUTADD,1,2,1", "CSECT,3,1,0,0,0.5,0,0,1", "ACTM,3",
                                       "CALCMOI,0,0,1,1,0,0", "CALCC"])
    mo, cc = ui.session["calcmoi"], ui.session["calcc"]
    assert [cc["Xc"], cc["Yc"], cc["Zc"]] == pytest.approx([mo["Xc"], mo["Yc"], mo["Zc"]], abs=1e-12)


def test_defect_loaded_stresses_survive_gcom_or_the_mismatch_is_reported(tmp_path):
    lines, boxes = box_mesh([0.0, 1.0, 2.0], [0.0, 1.0], [0.0, 1.0], group=4)
    cl.write_ess(tmp_path / "g.ess", {4: "SOLID"}, {(4, 1): [0, 0, 10, 0, 0, 0], (4, 2): [0, 0, 10, 0, 0, 0]})
    ui = Interpreter(cwd=tmp_path)
    run(ui, lines + ["READSTR,g.ess", "GCOM", "CUTADD,1,1,1-2", "CSECT,3,1,0,0,0.5,0,0,1", "ACTM,3",
                     "CALCPAR,0,0,1,1,0,0"])
    fz = ui.session["calcpar"]["Fz"]
    reported = any("group" in w.lower() and "stress" in w.lower() and "match" in w.lower()
                   for w in ui.sink.texts(Kind.WARNING))
    assert fz == pytest.approx(20.0) or reported
