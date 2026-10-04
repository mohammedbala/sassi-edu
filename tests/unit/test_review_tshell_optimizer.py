"""Adversarial review tests of work package 'tshell_optimizer' (TSHELL thick shell, HOUSE node optimizer,
THSHLSTR / THSHLSMH).

These tests probe the package from angles independent of the implementer's own tests and VPs:

* TSHELL statics against closed forms the implementation does not use: the Navier series of the
  simply supported Mindlin plate under uniform load (thin and thick, spec 08 section 12 test 15), the
  statically determinate resultants of a cantilever plate in an oblique plane (signs of MXX and QXZ),
  the membrane response with the free drilling rotation against SHELL with the drilling rotation fixed,
  the rigid-body modes of a folded (wall + slab) TSHELL assembly and the rigid-body STF (spec 05d
  test 3);
* the NAFEMS Test 21 reference re-derived as a 3 x 3 generalised eigenproblem per (m, n) mode (not the
  implementer's quadratic), and TSHELL convergence to the Mindlin closed form with the element's own
  shear factor 5/6;
* the node optimizer on a 2D PLANE model (dim = 1) with excavated soil and a scrambled numbering: the
  four matrices are exact permutations through the ``.map``, the interaction nodes keep their
  bottom-up table order (which removes the EDU-21 warning), FILE4 ``x_node_old_id`` is consistent with
  the coordinates, the equation profile does not grow, and the recovery operators (STRESS) are the
  same element by element.

Tests marked ``xfail(strict=True)`` document defects found in the review (the reason names the defect);
remove the marker when the defect is fixed.

Requirements exercised: 1.5 / 4.1 TSHELL, 4.4 item 5, 2.6 (optimizer check), D-ELM-08, D-HOU-03,
D-FIL-09, D-TSH-01, spec 05d test 3, spec 08 test 15, spec 11 section 4.1, R2 I.2.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from sassi.core import renumber as RN
from sassi.elements import (ElemRecord, assemble, build_dofmap, material_from_E_nu, natural_frequencies, shell,
                            tshell)
from sassi.io import decks
from sassi.io.files import read_container

pytestmark = [pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning"),
              pytest.mark.filterwarnings("ignore:.*encountered in.*:RuntimeWarning")]

E_ST, NU_ST, RHO_ST = 200e9, 0.3, 8000.0


# ======================================================================================
# helpers (independent of sassi.verify.problems.vp_tshell)
# ======================================================================================
def _rotation(seed: int) -> np.ndarray:
    Q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    return Q if np.linalg.det(Q) > 0 else -Q


def _grid(n: int, a: float, b: float, mat, t: float, eint: int, code: int = 5, triangles: str = ""):
    """n x n mesh of the rectangle [0, a] x [0, b] in the X-Y plane: quads, or two triangles per cell
    (``triangles`` 'right': all diagonals (i, j)-(i+1, j+1); 'cross': alternating diagonals)."""
    nid = lambda i, j: 1 + i + j * (n + 1)
    xyz = {nid(i, j): (a * i / n, b * j / n, 0.0) for j in range(n + 1) for i in range(n + 1)}
    props = dict(thick=t) if code == 3 else dict(thick=t, eint=eint)
    recs = []
    for j in range(n):
        for i in range(n):
            c = (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1))
            if not triangles:
                tris = [c]
            elif triangles == "right" or (i + j) % 2 == 0:
                tris = [(c[0], c[1], c[2]), (c[0], c[2], c[3])]
            else:
                tris = [(c[0], c[1], c[3]), (c[1], c[2], c[3])]
            for nodes in tris:
                recs.append(ElemRecord(1, len(recs) + 1, code, nodes, False, mat, dict(props)))
    return xyz, recs


def _assemble(xyz, recs, fixfun, undamped=True):
    ids = sorted(xyz)
    fix = np.array([fixfun(*xyz[k]) for k in ids], dtype=int)
    dm = build_dofmap(ids, fix, recs)
    return dm, assemble(xyz, recs, dm, undamped=undamped)


def _navier_center(a: float, t: float, E: float, nu: float, q: float = 1.0, kappa: float = 5.0 / 6.0,
                   nterms: int = 199) -> float:
    """Centre deflection of a simply supported square Mindlin plate under uniform load q (hard simple
    supports): w_mn = q_mn (1/(D k^4) + 1/(kappa G t k^2)), q_mn = 16 q/(pi^2 m n), m, n odd."""
    D = E * t ** 3 / (12.0 * (1.0 - nu * nu))
    G = E / (2.0 * (1.0 + nu))
    w = 0.0
    for m in range(1, nterms + 1, 2):
        for n in range(1, nterms + 1, 2):
            k2 = (m * math.pi / a) ** 2 + (n * math.pi / a) ** 2
            w += (16.0 * q / (math.pi ** 2 * m * n) * (1.0 / (D * k2 * k2) + 1.0 / (kappa * G * t * k2))
                  * math.sin(m * math.pi / 2) * math.sin(n * math.pi / 2))
    return w


def _ss_plate_center(n: int, t: float, eint: int, a: float = 10.0, code: int = 5, triangles: str = "") -> float:
    """Centre deflection of the hard simply supported plate (w = 0 and tangential rotation 0 on the
    edges, in-plane restrained) under a unit uniform load lumped to the nodes."""
    mat = material_from_E_nu(E_ST, NU_ST, RHO_ST)
    xyz, recs = _grid(n, a, a, mat, t, eint, code, triangles)

    def fixfun(x, y, z):
        ex = abs(x) < 1e-9 or abs(x - a) < 1e-9
        ey = abs(y) < 1e-9 or abs(y - a) < 1e-9
        return [1, 1, int(ex or ey), int(ex), int(ey), int(code == 3)]
    dm, am = _assemble(xyz, recs, fixfun)
    h = a / n
    F = np.zeros(dm.neq)
    for k, (x, y, _) in xyz.items():
        e = dm.eq(k, 3)
        if e >= 0:
            F[e] += (h if 0 < x < a else h / 2) * (h if 0 < y < a else h / 2)
    u = spla.spsolve(sp.csc_matrix(am.Ks.real), F)
    c = next(k for k, (x, y, _) in xyz.items() if abs(x - a / 2) < 1e-9 and abs(y - a / 2) < 1e-9)
    return float(u[dm.eq(c, 3)])


def _mindlin_3x3(a: float, t: float, E: float, nu: float, rho: float, kappa: float, nmodes: int = 8) -> np.ndarray:
    """Flexural frequencies (Hz) of the hard simply supported square Mindlin plate with rotary inertia from
    the 3 x 3 generalised eigenproblem in (W, Psi_x, Psi_y) of each (m, n) Navier mode (Mindlin 1951):
    w = W sin sin, psi_x = X cos sin, psi_y = Y sin cos; the lowest root of each (m, n) is flexural."""
    G = E / (2 * (1 + nu))
    D = E * t ** 3 / (12 * (1 - nu * nu))
    kGt, I = kappa * G * t, rho * t ** 3 / 12
    out = []
    for m in range(1, 7):
        for n in range(1, 7):
            al, be = m * math.pi / a, n * math.pi / a
            k2 = al * al + be * be
            K = np.array([[kGt * k2, kGt * al, kGt * be],
                          [kGt * al, D / 2 * ((1 - nu) * k2 + (1 + nu) * al * al) + kGt, D / 2 * (1 + nu) * al * be],
                          [kGt * be, D / 2 * (1 + nu) * al * be, D / 2 * ((1 - nu) * k2 + (1 + nu) * be * be) + kGt]])
            w2 = sla.eigh(K, np.diag([rho * t, I, I]), eigvals_only=True)
            out.append(math.sqrt(w2[0]) / (2 * math.pi))
    return np.sort(out)[:nmodes]


# ======================================================================================
# 1. TSHELL statics and dynamics against independent closed forms
# ======================================================================================
@pytest.mark.parametrize("ta", [1e-3, 0.1])
@pytest.mark.parametrize("eint", [0, 1])
def test_ss_plate_uniform_load_converges_to_navier(ta, eint):
    """Spec 08 section 12 test 15: square plate, simply supported, uniform load -> Navier (Kirchhoff in the
    thin limit, Mindlin with kappa = 5/6 for t/a = 1/10); both EINT converge at O(h^2), EINT 0 (one-point
    bending) from above, EINT 1 (2x2 bending) from below."""
    a = 10.0
    ref = _navier_center(a, ta * a, E_ST, NU_ST)
    e8, e16 = (_ss_plate_center(n, ta * a, eint) / ref - 1.0 for n in (8, 16))
    assert abs(e16) < 2e-3, (e8, e16)
    assert 3.0 < e8 / e16 < 5.0                              # O(h^2)
    assert (e16 > 0) if eint == 0 else (e16 < 0)


@pytest.mark.parametrize("t", [0.01, 1.0])
@pytest.mark.parametrize("eint", [0, 1])
def test_cantilever_resultants_equal_statics_in_an_oblique_plane(t, eint):
    """Cantilever plate L x b (nu = 0) clamped at x' = 0, unit total transverse load along the plate normal
    at the tip, rotated into an oblique plane: the recovered local resultants are the statically
    determinate values MXX = -P (L - x_c)/b (tension on the -z' face, D-STR-10 sign of M = int sigma z dz)
    and QXZ = +P/b at every element centre; membrane forces, MYY and QYZ vanish."""
    L, b, nx, ny = 10.0, 2.0, 8, 2
    R = _rotation(21)
    mat = material_from_E_nu(E_ST, 0.0, RHO_ST)
    nid = lambda i, j: 1 + i + j * (nx + 1)
    xyz0 = {nid(i, j): (L * i / nx, b * j / ny, 0.0) for j in range(ny + 1) for i in range(nx + 1)}
    recs = [ElemRecord(1, k + 1, 5, (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)), False, mat,
                       dict(thick=t, eint=eint)) for k, (i, j) in enumerate((i, j) for j in range(ny) for i in range(nx))]
    xyz = {k: tuple(R @ np.asarray(p)) for k, p in xyz0.items()}
    ids = sorted(xyz)
    fix = np.array([[1] * 6 if abs(xyz0[k][0]) < 1e-12 else [0] * 6 for k in ids])
    dm = build_dofmap(ids, fix, recs)
    am = assemble(xyz, recs, dm, undamped=True)
    F = np.zeros(dm.neq)
    for k in ids:
        x, y, _ = xyz0[k]
        if abs(x - L) < 1e-12:
            w = (b / ny) * (0.5 if (abs(y) < 1e-12 or abs(y - b) < 1e-12) else 1.0) / b
            for d in range(3):
                F[dm.eq(k, d + 1)] += w * R[d, 2]
    u = spla.spsolve(sp.csc_matrix(am.Ks.real), F)
    rec = am.recovery["TSHELL"]
    s = np.einsum("ecd,ed->ec", rec["S"].real, np.where(rec["eq"] >= 0, u[np.maximum(rec["eq"], 0)], 0.0))
    xc = np.array([np.mean([xyz0[n][0] for n in r.nodes]) for r in recs])
    M_ref, Q_ref = -(L - xc) / b, 1.0 / b
    assert np.allclose(s[:, 5], M_ref, rtol=0, atol=1e-6 * np.abs(M_ref).max())
    assert np.allclose(s[:, 3], Q_ref, rtol=1e-6)
    assert np.abs(s[:, [0, 1, 2, 4, 6]]).max() < 1e-6 * np.abs(M_ref).max()


def test_drilling_penalty_does_not_stiffen_the_membrane():
    """In-plane bending of a cantilever strip: TSHELL with the drilling rotation free (penalty tied to the
    membrane rotation, D-ELM-08) against SHELL (same Q4 + incompatible-mode membrane) with the drilling
    rotation fixed: the penalty may only stiffen marginally (< 1e-3)."""
    L, b, n_x, n_y, t = 10.0, 2.0, 16, 4, 0.1
    mat = material_from_E_nu(E_ST, 0.25, RHO_ST)
    out = {}
    for code in (3, 5):
        nid = lambda i, j: 1 + i + j * (n_x + 1)
        xyz = {nid(i, j): (L * i / n_x, b * j / n_y, 0.0) for j in range(n_y + 1) for i in range(n_x + 1)}
        props = dict(thick=t) if code == 3 else dict(thick=t, eint=0)
        recs = [ElemRecord(1, k + 1, code, (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)), False, mat,
                           dict(props)) for k, (i, j) in enumerate((i, j) for j in range(n_y) for i in range(n_x))]
        dm, am = _assemble(xyz, recs, lambda x, y, z: [1] * 6 if x < 1e-9 else [0, 0, 1, 1, 1, int(code == 3)])
        F = np.zeros(dm.neq)
        tip = [k for k, p in xyz.items() if abs(p[0] - L) < 1e-9]
        for k in tip:
            F[dm.eq(k, 2)] += (b / n_y) * (0.5 if xyz[k][1] in (0.0, b) else 1.0) / b
        u = spla.spsolve(sp.csc_matrix(am.Ks.real), F)
        out[code] = np.mean([u[dm.eq(k, 2)] for k in tip])
    rel = 1.0 - out[5] / out[3]
    assert 0.0 <= rel < 1e-3, rel


@pytest.mark.parametrize("t", [0.05, 0.5])
def test_folded_plate_has_exactly_six_rigid_body_modes(t):
    """A free wall + slab TSHELL assembly meeting at 90 degrees: the drilling penalty of each facet couples
    with the bending rotation of the other, no spurious mechanism, exactly six zero-energy modes."""
    mat = material_from_E_nu(E_ST, NU_ST, RHO_ST)
    n = 3
    xyz, key = {}, {}
    for side in (0, 1):
        for j in range(n + 1):
            for i in range(n + 1):
                p = (float(i), float(j), 0.0) if side == 0 else (float(i), 0.0, float(j))
                if p not in key:
                    key[p] = len(key) + 1
                    xyz[key[p]] = p
    recs = []
    for side in (0, 1):
        pt = (lambda i, j: (float(i), float(j), 0.0)) if side == 0 else (lambda i, j: (float(i), 0.0, float(j)))
        for j in range(n):
            for i in range(n):
                c = tuple(key[pt(*ij)] for ij in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)))
                recs.append(ElemRecord(1, len(recs) + 1, 5, c, False, mat, dict(thick=t, eint=0)))
    _, am = _assemble(xyz, recs, lambda *p: [0] * 6)
    K = am.Ks.toarray().real
    w = np.linalg.eigvalsh(0.5 * (K + K.T))
    assert int(np.sum(np.abs(w) < 1e-9 * np.abs(w).max())) == 6


@pytest.mark.parametrize("n_rows", [3, 4])
def test_rigid_body_motion_gives_zero_resultants(n_rows):
    """Spec 05d test 3: the six rigid-body motions (incl. the drilling rotation theta = omega) of a distorted
    TSHELL quad / triangle in an oblique plane give zero recovered resultants (STF = 0)."""
    R = _rotation(5)
    X = (np.array([[0.0, 0.0, 0.0], [1.3, 0.2, 0.0], [1.1, 0.95, 0.0], [-0.15, 1.05, 0.0]]) @ R.T)[:n_rows]
    S = tshell.recovery(X, material_from_E_nu(3.0e7, 0.2, 2.4), thick=0.3, eint=0)
    c = X.mean(axis=0)
    for k in range(6):
        u = np.zeros(6 * n_rows)
        for a in range(n_rows):
            if k < 3:
                u[6 * a + k] = 1.0
            else:
                ax = np.eye(3)[k - 3]
                u[6 * a:6 * a + 3] = np.cross(ax, X[a] - c)
                u[6 * a + 3:6 * a + 6] = ax
        assert np.abs(S @ u).max() < 1e-12 * np.abs(S).max() * np.abs(u).max()


def test_test21_reference_by_an_independent_mindlin_formulation():
    """NAFEMS Test 21 (R2 I.2): the Mindlin closed form re-derived as a 3x3 generalised eigenproblem per mode
    (kappa = pi^2/12, rotary inertia) reproduces the targets (2e-4) and the implementer's quadratic."""
    from sassi.verify.problems.vp_tshell import TEST21_REF, mindlin_ss_frequencies
    f = _mindlin_3x3(10.0, 1.0, E_ST, NU_ST, RHO_ST, math.pi ** 2 / 12)
    assert np.allclose(f, TEST21_REF, rtol=2e-4)
    assert np.allclose(f, mindlin_ss_frequencies(10.0, 10.0, 1.0, E_ST, NU_ST, RHO_ST), rtol=1e-12)


@pytest.mark.parametrize("eint", [0, 1])
def test_test21_tshell_converges_to_its_own_mindlin_theory(eint):
    """Test 21 with TSHELL (kappa = 5/6): the first four frequencies converge at O(h^2) from below to the
    kappa = 5/6 Mindlin closed form (16 x 16 within 1.5 %, error ratio 16/32 about 4)."""
    from sassi.verify.problems.vp_tshell import test21
    ref = _mindlin_3x3(10.0, 1.0, E_ST, NU_ST, RHO_ST, 5.0 / 6.0)
    e16 = test21(16, eint)[:4] / ref[:4] - 1.0
    e32 = test21(32, eint)[:4] / ref[:4] - 1.0
    assert np.all(e16 < 0) and np.all(np.abs(e16) < 0.015)
    assert np.all((3.5 < e16 / e32) & (e16 / e32 < 4.5))


# ======================================================================================
# 2. THSHLSTR face stresses: a hand-computed case
# ======================================================================================
def test_face_stresses_hand_computed():
    """D-TSH-01 with numbers worked by hand: t = 0.5, |N| = (100, 40, 10), |M| = (2, 1, 0.5), E = 1e4,
    nu = 0.25.  '++': SXX = 200 + 48 = 248, SYY = 80 + 24 = 104, TXY = 20 + 12 = 32; '+-': SXX = 152,
    SYY = 56; the four rows contain both actual faces of a state with N > 0 and M < 0."""
    t, E, nu = 0.5, 1.0e4, 0.25
    out = tshell.face_stresses([100.0, -40.0, 10.0, 3.0, -6.0, -2.0, 1.0, 0.5], t, E, nu)
    pp, pm = out["rows"]["++"], out["rows"]["+-"]
    assert np.allclose(pp[:3], [248.0, 104.0, 32.0]) and np.allclose(pm[:3], [152.0, 56.0, 32.0])
    c, r = (248.0 + 104.0) / 2, math.hypot((248.0 - 104.0) / 2, 32.0)
    assert np.allclose(pp[3:5], [c + r, c - r])
    assert np.allclose(pp[5:8], [(248.0 - 0.25 * 104.0) / E, (104.0 - 0.25 * 248.0) / E, 32.0 * 2.5 / E])
    assert out["TXZ"] == pytest.approx(9.0) and out["TYZ"] == pytest.approx(18.0)
    assert np.allclose(out["rows"]["--"][:2], -pp[:2]) and np.allclose(out["rows"]["-+"][:2], -pm[:2])


# ======================================================================================
# 3. Node optimizer on a 2D PLANE model (dim = 1) with excavated soil
# ======================================================================================
def _plane2d_deck(seed: int = 7) -> decks.Deck:
    """2D (dim = 1) model: 4 excavated PLANE elements (layer 1, z = -2..0, all nodes interaction nodes on
    interfaces 1 and 2, listed bottom-up) under 4 structural PLANE elements (z = 0..2); node numbers
    scrambled at random, so the interaction numbers are not ascending (EDU-21 warning)."""
    from sassi.verify.problems.vp_house import add_element, add_node, house_deck
    d = house_deck(dim=1)
    ids = {}
    for iz, z in enumerate((-2.0, 0.0, 2.0)):
        for ix in range(5):
            ids[(ix, iz)] = len(ids) + 1
    perm = np.random.default_rng(seed).permutation(len(ids)) + 1
    num = {k: int(perm[v - 1]) for k, v in ids.items()}
    for (ix, iz), n in sorted(num.items(), key=lambda kv: kv[1]):
        add_node(d, n, float(ix), 0.0, (-2.0, 0.0, 2.0)[iz])
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("groups").append([1, 4, "excavated soil"])
    d.table("groups").append([2, 4, "structure"])
    for ix in range(4):
        add_element(d, 1, ix + 1, [num[(ix, 0)], num[(ix + 1, 0)], num[(ix + 1, 1)], num[(ix, 1)]], etype=2, mat=1)
        add_element(d, 2, ix + 1, [num[(ix, 1)], num[(ix + 1, 1)], num[(ix + 1, 2)], num[(ix, 2)]], etype=1, mat=1)
    for iz in (0, 1):
        for ix in range(5):
            d.table("interaction").append([num[(ix, iz)]])
    return d


def _profile(A) -> int:
    A = sp.coo_matrix(A)
    first = np.arange(A.shape[0])
    np.minimum.at(first, A.row, A.col)
    return int(np.sum(np.arange(A.shape[0]) - first))


def test_optimizer_2d_plane_model_is_an_exact_permutation(tmp_path):
    from sassi.verify.problems.vp_house import run_house
    a, b = tmp_path / "a", tmp_path / "b"
    d = _plane2d_deck()
    rc, out_a = run_house(a, d)
    assert rc == 0, out_a[-2000:]
    d["optimize"] = 1
    rc, out_b = run_house(b, d)
    assert rc == 0, out_b[-2000:]
    mp = RN.read_map(b / "m.map")
    lines = (b / "m.map").read_text().split("\n")
    assert all(len(ln.split()) == 2 for ln in lines if ln.strip()) and len(mp) == 15
    assert sorted(mp.values()) == list(range(1, 16))
    fa, fb = read_container(a / "m.N4", "FILE4"), read_container(b / "m.N4", "FILE4")
    ka, kb = read_container(a / "COOSK", "COOSK"), read_container(b / "COOSK", "COOSK")
    ma, mb = read_container(a / "COOSM", "COOSM"), read_container(b / "COOSM", "COOSM")
    assert set(int(v) for v in fb["eq_dof"]) == {1, 3}             # 2D: UX, UZ only
    eqb = {(int(n), int(k)): i for i, (n, k) in enumerate(zip(fb["eq_node"], fb["eq_dof"]))}
    p = np.array([eqb[(mp[int(n)], int(k))] for n, k in zip(fa["eq_node"], fa["eq_dof"])])
    for A, B in ((ka.sparse("Ks"), kb.sparse("Ks")), (ka.sparse("Ke"), kb.sparse("Ke")),
                 (ma.sparse("Ms"), mb.sparse("Ms")), (ma.sparse("Me"), mb.sparse("Me"))):
        assert abs(B.tocsr()[p][:, p] - A).max() <= 1e-14 * max(abs(A).max(), 1e-300)
    # interaction nodes: same physical nodes in the bottom-up table order, now with ascending new numbers
    ia, ib = [int(v) for v in fa["int_node"]], [int(v) for v in fb["int_node"]]
    assert not all(np.diff(ia) > 0)                            # scrambled originally
    assert ib == [mp[n] for n in ia] and all(np.diff(ib) > 0)
    assert "EDU-21" not in out_b
    # x_node_old_id: the new node k sits where the old node x_node_old_id[k] was
    xa = {int(n): np.asarray(x) for n, x in zip(fa["node_id"], fa["node_xyz"])}
    for n, o, x in zip(fb["node_id"], fb["x_node_old_id"], fb["node_xyz"]):
        assert mp[int(o)] == int(n) and np.allclose(x, xa[int(o)])
    # profile of the structure + excavated-soil pattern does not grow
    pat = lambda c: abs(c.sparse("Ks")) + abs(c.sparse("Ke"))
    assert _profile(pat(kb)) <= _profile(pat(ka))
    # STRESS operators: same S per element, equations of the same (node, dof) through the map
    for name in ("PLANE",):
        Sa, Sb = fa[f"rec_{name}_S"], fb[f"rec_{name}_S"]
        assert np.allclose(Sa, Sb, rtol=0, atol=1e-14 * np.abs(Sa).max())
        qa, qb = fa[f"rec_{name}_eq"], fb[f"rec_{name}_eq"]
        node_a = lambda e: (-1, -1) if e < 0 else (mp[int(fa["eq_node"][e])], int(fa["eq_dof"][e]))
        node_b = lambda e: (-1, -1) if e < 0 else (int(fb["eq_node"][e]), int(fb["eq_dof"][e]))
        assert all(node_a(x) == node_b(y) for x, y in zip(qa.ravel(), qb.ravel()))


def test_optimizer_is_idempotent_on_its_own_output(tmp_path):
    """Optimising the .hounew deck again changes nothing (identity map)."""
    from sassi.verify.problems.vp_house import run_house
    d = _plane2d_deck(seed=11)
    d["optimize"] = 1
    rc, out = run_house(tmp_path / "a", d)
    assert rc == 0, out[-1500:]
    dn = decks.read(tmp_path / "a" / "m.hounew", "HOUSE")
    assert int(dn["optimize"]) == 0
    dn["optimize"] = 1
    rc, out = run_house(tmp_path / "b", dn)
    assert rc == 0, out[-1500:]
    assert all(k == v for k, v in RN.read_map(tmp_path / "b" / "m.map").items())


def test_thshlstr_file_recreated_by_inp_after_mdl(tmp_path):
    """WRITE -> MDL -> INP: the .pre restores the THSHLSTR record and re-creates the THSHLSTR.opt that
    STRESS reads (WRITE does not write MDL, so MDL must come first; see the strict xfail below)."""
    from sassi.modules.stress import thshlstr_flag
    from sassi.prep import Interpreter
    ui = Interpreter(cwd=tmp_path)
    assert ui.execute(f"MDL,m,{tmp_path}") and ui.execute("THSHLSTR,1")
    assert ui.execute(f"WRITE,{tmp_path / 'w.pre'}")
    (tmp_path / "THSHLSTR.opt").unlink()
    ui2 = Interpreter(cwd=tmp_path)
    ui2.execute(f"MDL,m,{tmp_path}")
    ui2.execute(f"INP,{tmp_path / 'w.pre'}")
    assert thshlstr_flag(None, tmp_path)[0] == 1


# ======================================================================================
# 4. Defects found by the review (strict xfail)
# ======================================================================================
def test_tshell_triangles_do_not_lock_in_thin_plate_bending():
    a, t = 10.0, 0.01
    ref = _navier_center(a, t, E_ST, NU_ST)
    assert _ss_plate_center(16, t, 0, triangles="right") / ref > 0.95


def test_natural_frequencies_of_an_oblique_tshell_plate():
    mat = material_from_E_nu(E_ST, NU_ST, RHO_ST)
    xyz, recs = _grid(3, 10.0, 10.0, mat, 1.0, 1)
    R = _rotation(2)
    _, flat = _assemble(xyz, recs, lambda *p: [0] * 6)
    f0 = natural_frequencies(flat.Ks, flat.Ms, 10)
    _, obl = _assemble({k: tuple(R @ np.asarray(p)) for k, p in xyz.items()}, recs, lambda *p: [0] * 6)
    f1 = natural_frequencies(obl.Ks, obl.Ms, 10)
    assert np.allclose(f1[6:], f0[6:], rtol=1e-8)


def test_failed_optimised_run_leaves_no_stale_map(tmp_path):
    from sassi.verify.problems.vp_house import mixed_house_deck, run_house
    rc, _ = run_house(tmp_path, mixed_house_deck())
    assert rc == 0
    d = mixed_house_deck()
    d["optimize"] = 1
    first = int(d.table("interaction").rows[0][0])
    for r in d.table("nodes").rows:
        if int(r[0]) == first:
            r[3] = float(r[3]) + 0.37                             # EDU-01: off the interface -> HOUSE fails
    rc, out = run_house(tmp_path, d)
    assert rc == 1 and "EDU-01" in out
    f4 = read_container(tmp_path / "m.N4", "FILE4")
    assert not (tmp_path / "m.map").exists() or int(f4.meta.get("x_optimized", 0)) == 1


def test_vibration_load_follows_the_physical_node_with_the_optimizer(tmp_path):
    from sassi.verify import builders as B
    from sassi.verify.problems.vp_tshell import optimizer_model
    site, mdl, d0 = optimizer_model()
    fs = B.FrequencySet.harmonic(1.0, [3, 9])
    B.run_soil(tmp_path, "m", site, fs, layer=mdl.layer, rad=mdl.rad)
    nodes = {int(r[0]): float(r[3]) for r in d0.table("nodes").rows}
    roof = sorted(n for n, z in nodes.items() if abs(z - max(nodes.values())) < 1e-9)
    load = roof[len(roof) // 2]
    d0["optimize"] = 0
    B.write_deck(tmp_path, "m", d0)
    B.run("HOUSE", tmp_path, "m")
    B.run_force(tmp_path, "m", [(load, 1, 1.0, 0.0)], fs)
    B.run_analys(tmp_path, "m", fs, type=1)
    Ha = B.tf(B.read_file8(tmp_path), load, 1)
    d1 = decks.read(tmp_path / "m.hou", "HOUSE")
    d1["optimize"] = 1
    B.write_deck(tmp_path, "m", d1)
    B.run("HOUSE", tmp_path, "m")
    mp = RN.read_map(tmp_path / "m.map")
    assert mp[load] != load                                        # the loaded node is renumbered
    B.run_force(tmp_path, "m", [(load, 1, 1.0, 0.0)], fs)          # FORCE deck from the interpreter: old number
    B.run_analys(tmp_path, "m", fs, type=1)
    Hb = B.tf(B.read_file8(tmp_path), mp[load], 1)
    out = B.listing(tmp_path, "m", "ANALYS").lower()
    warned = any(k in out for k in ("optimi", ".map", "renumber"))
    assert np.allclose(Hb, Ha, rtol=1e-8, atol=0) or warned


def test_optimizer_does_not_increase_the_bandwidth(tmp_path):
    from sassi.verify.problems.vp_house import mixed_house_deck, run_house
    d = mixed_house_deck()
    rc, _ = run_house(tmp_path / "a", d)
    assert rc == 0
    d["optimize"] = 1
    rc, out = run_house(tmp_path / "b", d)
    assert rc == 0
    ka, kb = (read_container(tmp_path / w / "COOSK", "COOSK") for w in ("a", "b"))
    bw = lambda c: int(np.max(np.abs(np.diff(sp.find(abs(c.sparse("Ks")) + abs(c.sparse("Ke")))[:2], axis=0))))
    assert _profile(abs(kb.sparse("Ks")) + abs(kb.sparse("Ke"))) <= _profile(abs(ka.sparse("Ks")) + abs(ka.sparse("Ke")))
    assert bw(kb) <= bw(ka)


def test_rcm_keeps_a_narrow_band_with_a_bottom_interaction_layer():
    nx, ny, nz = 12, 12, 4
    idx = lambda i, j, k: i + nx * (j + ny * k)
    cl = [[idx(i, j, k), idx(i + 1, j, k), idx(i + 1, j + 1, k), idx(i, j + 1, k), idx(i, j, k + 1), idx(i + 1, j, k + 1),
           idx(i + 1, j + 1, k + 1), idx(i, j + 1, k + 1)] for k in range(nz - 1) for j in range(ny - 1) for i in range(nx - 1)]
    z = np.array([k for k in range(nz) for j in range(ny) for i in range(nx)], dtype=float)
    keep = [idx(i, j, 0) for j in range(ny) for i in range(nx)]
    r = RN.optimize_numbering(nx * ny * nz, cl, keep=keep, z=z)
    assert all(np.diff(r.new_number[keep]) > 0)                    # the constraint itself holds
    assert r.profile_after <= r.profile_before
    assert r.bandwidth_after <= r.bandwidth_before


def test_tshell_outputs_have_no_nan_rows(tmp_path):
    from sassi.prep.commands import thickshell as TH
    from sassi.verify.problems.vp_house import run_house
    from sassi.verify.problems.vp_stress import add_eout, field_file8, harmonic_file, run_stress, stress_deck
    from sassi.verify.problems.vp_tshell import DF54, DT54, FNUM54, NFFT54, tshell_plate_deck
    d, (o, e1, e2) = tshell_plate_deck()
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out[-1500:]
    harmonic_file(tmp_path / "eq.acc", NFFT54, DT54, 7, 0.1)
    field_file8(tmp_path, FNUM54, DF54, lambda node, x, f: np.concatenate(
        [2e-4 * ((np.asarray(x) - o) @ e1) * e1, np.zeros(3)]).astype(complex), NFFT54, DT54)
    TH.write_thshlstr_file(tmp_path, 1)
    dk = stress_deck(NFFT54, DT54)
    add_eout(dk, 1, [1], [1] * 8)
    add_eout(dk, 1, [2], [0] * 8)
    rc, out = run_stress(tmp_path, dk)
    assert rc == 0, out[-1500:]
    for name in ("TSHELL_ELEMENT_MAX.TXT", "TSHELL_FACE_STRESSES.TXT"):
        assert "nan" not in (tmp_path / name).read_text().lower()


def test_housex_optimize_is_not_reported_unavailable(tmp_path):
    from sassi.prep import Interpreter, Kind
    ui = Interpreter(cwd=tmp_path)
    ui.execute(f"MDL,m,{tmp_path}")
    ui.execute("HOUSEX,1")
    warns = [m.text for m in ui.sink.messages if m.kind == Kind.WARNING]
    assert not any("optimizer" in w and "not available" in w for w in warns)


@pytest.mark.parametrize("scenario", ["mdl_change", "inp_then_mdl"])
def test_thshlstr_flag_reaches_stress(tmp_path, scenario):
    from sassi.modules.stress import thshlstr_flag
    from sassi.prep import Interpreter
    ui = Interpreter(cwd=tmp_path)
    assert ui.execute(f"MDL,m,{tmp_path / 'A'}") and ui.execute("THSHLSTR,1")
    if scenario == "mdl_change":
        assert ui.execute(f"MDL,m,{tmp_path / 'B'}")
    else:
        assert ui.execute(f"WRITE,{tmp_path / 'w.pre'}")
        ui = Interpreter(cwd=tmp_path)
        ui.execute(f"INP,{tmp_path / 'w.pre'}")
        assert ui.execute(f"MDL,m,{tmp_path / 'B'}")
    assert ui.model.options.record("THSHLSTR").values == ["1"]
    # lead fix: the flag travels in the STRESS deck written by AFWRITE (decks.STRESS 'thshlstr')
    from sassi.io import decks
    from sassi.prep.afwrite import DeckBuilder
    (tmp_path / "B").mkdir(exist_ok=True)
    deck = DeckBuilder(ui.model, [tmp_path / "B"]).stress()
    path = decks.write(tmp_path / "B" / "m.str", deck)
    assert thshlstr_flag(path, tmp_path / "B")[0] == 1
