"""TSHELL thick-shell element (sassi.elements.tshell; requirements 1.5 / 4.1, D-ELM-08, D-TSH-01)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.conventions import ELEMENT_COMPONENTS, cfactor
from sassi.elements import ELEMENTS, ElemRecord, ElementError, assemble, build_dofmap, material_from_E_nu, material_from_M
from sassi.elements import shell, tshell
from sassi.elements.assemble import element_dof_nodes

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")

MAT = material_from_E_nu(200e9, 0.3, 8000.0)
QUAD = np.array([[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [1.0, 0.9, 0.0], [-0.1, 1.1, 0.0]])
SQUARE = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])


def _rot(seed=3):
    Q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    return Q if np.linalg.det(Q) > 0 else -Q


def _nzero(K, rtol=1e-9):
    w = np.linalg.eigvalsh(0.5 * (K + K.T))
    return int(np.sum(np.abs(w) < rtol * np.abs(w).max()))


# ---------------------------------------------------------------------------------------- registry
def test_registered_with_eight_components():
    spec = ELEMENTS[5]
    assert spec.name == "TSHELL" and spec.nnodes == 4 and spec.dofs == (1, 2, 3, 4, 5, 6)
    assert spec.components == ELEMENT_COMPONENTS["TSHELL"] == ["NXX", "NYY", "NXY", "QXZ", "QYZ", "MXX", "MYY", "MXY"]


def test_dof_nodes_quads_and_triangles():
    rec = lambda n: ElemRecord(1, 1, 5, n, False, MAT, dict(thick=0.2))
    assert element_dof_nodes(rec((1, 2, 3, 4)))[0] == [1, 2, 3, 4]
    assert element_dof_nodes(rec((1, 2, 3, 0)))[0] == [1, 2, 3]
    assert element_dof_nodes(rec((1, 2, 3, 3)))[0] == [1, 2, 3]
    with pytest.raises(ElementError, match="Error 9"):
        element_dof_nodes(rec((1, 2, 1, 4)))
    with pytest.raises(ElementError, match="Error 7"):
        element_dof_nodes(rec((1, 2)))


# ---------------------------------------------------------------------------------------- stiffness
@pytest.mark.parametrize("eint", [0, 1])
@pytest.mark.parametrize("t", [0.01, 0.5])
def test_six_rigid_body_modes_and_symmetry(eint, t):
    R = _rot()
    for xyz in (QUAD @ R.T, (QUAD @ R.T)[:3]):
        K, M = tshell.matrices(xyz, MAT, thick=t, eint=eint)
        assert np.allclose(K, K.T, rtol=0, atol=1e-9 * np.abs(K).max())
        assert _nzero(K.real) == 6
        assert np.allclose(M, M.T)


def test_no_hourglass_stabilisation_leaves_spurious_modes():
    """EINT 0 without stabilisation: the two rotation hourglass modes (beta = xi eta) have no stiffness."""
    K, _ = tshell.matrices(SQUARE, MAT, thick=0.1, eint=0, hourglass=0.0)
    assert _nzero(K.real) == 8
    K, _ = tshell.matrices(SQUARE, MAT, thick=0.1, eint=0)
    assert _nzero(K.real) == 6


def test_eint0_equals_eint1_with_full_stabilisation_on_parallelograms():
    para = np.array([[0.0, 0, 0], [2.0, 0, 0], [2.6, 1.5, 0], [0.6, 1.5, 0]])
    K1, _ = tshell.matrices(para, MAT, thick=0.2, eint=1)
    K0, _ = tshell.matrices(para, MAT, thick=0.2, eint=0, hourglass=1.0)
    assert np.allclose(K0, K1, rtol=1e-12, atol=1e-12 * np.abs(K1).max())
    K0d, _ = tshell.matrices(para, MAT, thick=0.2, eint=0)
    assert not np.allclose(K0d, K1, rtol=1e-6)


def test_rotation_invariance():
    R = _rot(5)
    K, M = tshell.matrices(QUAD, MAT, thick=0.3, eint=0)
    Kr, Mr = tshell.matrices(QUAD @ R.T, MAT, thick=0.3, eint=0)
    T = np.kron(np.eye(8), R)
    assert np.allclose(T @ K @ T.T, Kr, atol=1e-9 * np.abs(K).max())
    assert np.allclose(T @ M @ T.T, Mr, atol=1e-12 * np.abs(M).max())


def test_complex_modulus():
    m = material_from_M(1, 3.0e7, 0.2, 24.0, 0.05, 0.05, 9.81)
    Kc, _ = tshell.matrices(QUAD, m, thick=0.3)
    K0, _ = tshell.matrices(QUAD, m.undamped(), thick=0.3)
    assert np.allclose(Kc, cfactor(0.05) * K0.real, rtol=1e-12, atol=1e-12 * np.abs(K0).max())


def test_drilling_penalty_free_for_rigid_rotation_and_scaled():
    """theta_z tied to omega = (v,x - u,y)/2: a rigid in-plane rotation has no energy; the drilling
    diagonal is 1e-4 x min membrane diagonal x area (D-ELM-08)."""
    K, _ = tshell.matrices(SQUARE, MAT, thick=0.2)
    u = np.zeros(24)
    for a, (x, y, _) in enumerate(SQUARE):
        u[6 * a:6 * a + 3] = [-y, x, 0.0]
        u[6 * a + 5] = 1.0
    assert abs(u @ K.real @ u) < 1e-9 * np.abs(K).max()
    D = tshell.plate_rigidities(MAT.E0, MAT.nu0, 0.2)[0]
    Km, _ = shell.membrane_q4(SQUARE[:, :2], D)
    kd = 1e-4 * np.min(np.diag(Km)) * 1.0
    xi = np.zeros(24)
    xi[5] = 1.0                                         # theta_z of node 1 alone (omega = 0)
    assert xi @ K.real @ xi == pytest.approx(kd, rel=1e-12)


def test_triangle_three_rows_equals_repeated_node():
    K3, M3 = tshell.matrices(QUAD[:3], MAT, thick=0.2)
    K4, M4 = tshell.matrices(QUAD[[0, 1, 2, 2]], MAT, thick=0.2)
    assert np.allclose(K4[:18, :18], K3) and np.all(K4[18:] == 0) and np.all(K4[:, 18:] == 0)
    assert np.allclose(M4[:18, :18], M3)
    S3 = tshell.recovery(QUAD[:3], MAT, thick=0.2)
    S4 = tshell.recovery(QUAD[[0, 1, 2, 2]], MAT, thick=0.2)
    assert S4.shape == (8, 24) and np.allclose(S4[:, :18], S3)


def test_invalid_input():
    with pytest.raises(ElementError, match="thickness"):
        tshell.matrices(QUAD, MAT, thick=0.0)
    with pytest.raises(ElementError, match="EINT"):
        tshell.matrices(QUAD, MAT, thick=0.1, eint=2)
    with pytest.raises(ElementError, match="EDU-05"):                       # re-entrant (dart) quadrilateral
        tshell.matrices(np.array([[0, 0, 0], [2, 0, 0], [0.6, 0.6, 0], [0, 2, 0.0]]), MAT, thick=0.1)
    with pytest.raises(ElementError, match="TSHELL"):                       # bow-tie: no x' axis
        tshell.matrices(np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0.0]]), MAT, thick=0.1)


# ---------------------------------------------------------------------------------------- mass
def test_lumped_mass_with_rotary_inertia():
    R = _rot(7)
    t = 0.4
    _, M = tshell.matrices(QUAD @ R.T, MAT, thick=t)
    A = shell.polygon_area(QUAD[:, :2])
    for a in range(4):
        assert np.allclose(M[6 * a:6 * a + 3, 6 * a:6 * a + 3], MAT.rho * t * A / 4 * np.eye(3))
        n = R[:, 2]
        Irot = MAT.rho * t ** 3 / 12 * A / 4
        assert np.allclose(M[6 * a + 3:6 * a + 6, 6 * a + 3:6 * a + 6], Irot * (np.eye(3) - np.outer(n, n)))
    assert np.count_nonzero(M[:6, 6:]) == 0                     # lumped: no coupling between nodes


# ---------------------------------------------------------------------------------------- recovery
def _local_field(xy, Lam, f):
    """Global element DOF vector of a field given in the element's local axes: f(x', y') -> [u' v' w' tx' ty' tz']."""
    u = np.zeros(6 * len(xy))
    for a, (x, y) in enumerate(xy):
        v = np.asarray(f(x, y), dtype=float)
        u[6 * a:6 * a + 3] = Lam.T @ v[:3]
        u[6 * a + 3:6 * a + 6] = Lam.T @ v[3:]
    return u


@pytest.mark.parametrize("eint", [0, 1])
@pytest.mark.parametrize("n", [3, 4])
def test_recovery_of_constant_states(eint, n):
    """Constant membrane, bending and transverse-shear states (local axes) are recovered exactly at the
    centre, in an oblique plane."""
    X = QUAD[:n] @ _rot(11).T
    S = tshell.recovery(X, MAT, thick=0.2, eint=eint).real
    Dm, Db, Ds = tshell.plate_rigidities(MAT.E0, MAT.nu0, 0.2)
    Lam, _, xyp = shell.local_frame(X)
    eps = np.array([1e-4, -2e-4, 3e-4])
    u = _local_field(xyp, Lam, lambda x, y: [eps[0] * x + 0.5 * eps[2] * y, eps[1] * y + 0.5 * eps[2] * x, 0, 0, 0, 0])
    assert np.allclose(S @ u, np.concatenate([Dm @ eps, [0, 0, 0, 0, 0]]), atol=1e-9 * np.abs(Dm @ eps).max())
    k1, k2, k12 = 1e-3, 2e-3, 0.5e-3
    u = _local_field(xyp, Lam, lambda x, y: [0, 0, 0.5 * k1 * x * x + 0.5 * k2 * y * y + k12 * x * y,
                                            k2 * y + k12 * x, -(k1 * x + k12 * y), 0])
    m = Db @ np.array([-k1, -k2, -2 * k12])
    out = S @ u
    assert np.allclose(out[5:], m, atol=1e-9 * np.abs(m).max()) and np.abs(out[3:5]).max() < 1e-9 * np.abs(m).max()
    assert np.abs(out[:3]).max() < 1e-9 * np.abs(m).max()
    g = np.array([2e-4, -1e-4])
    shear = lambda x, y: [0, 0, g[0] * x + g[1] * y, 0, 0, 0]
    if n == 4:
        u = _local_field(xyp, Lam, shear)
        assert np.allclose((S @ u)[3:5], Ds @ g, atol=1e-10 * np.abs(Ds @ g).max())
    else:
        # triangle: the strain operators reproduce constant shear with the exact rotation field (bubble 0);
        # the condensed recovery relaxes this unloaded, non-equilibrium state through the bubble (the
        # equilibrium state is test_triangle_constant_shear_equilibrium_patch)
        op = tshell.tri_plate_operators(xyp, Db, Ds, 0.2)
        ub = np.concatenate([np.asarray(shear(x, y), float)[2:5] for x, y in xyp] + [np.zeros(2)])
        assert np.allclose(Ds @ op["Bs"] @ ub, Ds @ g, rtol=1e-12, atol=0)


def test_assembly_through_records_and_recovery_padding():
    node_xyz = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0), 5: (2, 0, 0)}
    recs = [ElemRecord(1, 1, 5, (1, 2, 3, 4), False, MAT, dict(thick=0.2, eint=1)),
            ElemRecord(1, 2, 5, (2, 5, 3, 0), False, MAT, dict(thick=0.2))]
    dm = build_dofmap(sorted(node_xyz), None, recs)
    am = assemble(node_xyz, recs, dm)
    rec = am.recovery["TSHELL"]
    assert rec["S"].shape == (2, 8, 24) and list(rec["eq"][1][-6:]) == [-1] * 6
    K1, _ = tshell.matrices(np.array([node_xyz[n] for n in (1, 2, 3, 4)], float), MAT, thick=0.2, eint=1)
    K2, _ = tshell.matrices(np.array([node_xyz[n] for n in (2, 5, 3)], float), MAT, thick=0.2)
    ref = np.zeros((dm.neq, dm.neq), complex)
    for K, nodes in ((K1, (1, 2, 3, 4)), (K2, (2, 5, 3))):
        e = dm.eqs(nodes, range(1, 7))
        ref[np.ix_(e, e)] += K
    assert np.allclose(am.Ks.toarray(), ref)


# ---------------------------------------------------------------------------------------- locking
def _cantilever(t, nx=8, ny=2, eint=0, L=10.0, b=2.0):
    mat = material_from_E_nu(200e9, 0.0, 8000.0)
    nid = lambda i, j: 1 + i + j * (nx + 1)
    node_xyz = {nid(i, j): (L * i / nx, b * j / ny, 0.0) for j in range(ny + 1) for i in range(nx + 1)}
    recs = [ElemRecord(1, k + 1, 5, (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)), False, mat,
                       dict(thick=t, eint=eint)) for k, (i, j) in enumerate((i, j) for j in range(ny) for i in range(nx))]
    ids = sorted(node_xyz)
    fix = np.array([[1] * 6 if node_xyz[n][0] == 0 else [0] * 6 for n in ids])
    dm = build_dofmap(ids, fix, recs)
    am = assemble(node_xyz, recs, dm)
    F = np.zeros(dm.neq)
    tip = [n for n in ids if node_xyz[n][0] == L]
    for n in tip:
        y = node_xyz[n][1]
        F[dm.eq(n, 3)] = (b / ny) * (0.5 if y in (0, b) else 1.0) / b
    u = np.linalg.solve(am.Ks.toarray().real, F)
    return np.mean([u[dm.eq(n, 3)] for n in tip]) / (L ** 3 / (3 * 200e9 * b * t ** 3 / 12))


@pytest.mark.parametrize("eint", [0, 1])
def test_no_shear_locking_thin_cantilever(eint):
    assert _cantilever(0.01, eint=eint) == pytest.approx(1.0 - 1.0 / 256.0, rel=1e-5)    # t/L = 1/1000


# ---------------------------------------------------------------------------------------- triangle
TRI = np.array([[0.0, 0.0], [1.3, 0.2], [0.4, 0.9]])


def _old_mitc3(xy, Db, Ds):
    """Independent plain MITC3 plate (linear rotations, RT0 shear integrated at the side mid-points)."""
    A = shell.polygon_area(xy)
    x, y = xy[:, 0], xy[:, 1]
    dN = np.column_stack([[y[1] - y[2], y[2] - y[0], y[0] - y[1]], [x[2] - x[1], x[0] - x[2], x[1] - x[0]]]) / (2 * A)
    Bb = np.zeros((3, 9))
    Bb[0, 2::3], Bb[1, 1::3], Bb[2, 1::3], Bb[2, 2::3] = dN[:, 0], -dN[:, 1], -dN[:, 0], dN[:, 1]
    K = A * Bb.T @ Db @ Bb
    # field a + c (-(y - yc), x - xc): its tangential integral along each side d = x_j - x_i equals that of
    # the displacement-based strain, (w_j - w_i) + (beta_i + beta_j) . d / 2  (beta_x = theta_y, beta_y = -theta_x)
    xc = xy.mean(axis=0)
    rows, rhs = [], []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        d, mid = xy[j] - xy[i], 0.5 * (xy[i] + xy[j]) - xc
        rows.append([d[0], d[1], -mid[1] * d[0] + mid[0] * d[1]])
        r = np.zeros(9)
        r[3 * j], r[3 * i] = 1.0, -1.0
        for k in (i, j):
            r[3 * k + 2] += d[0] / 2
            r[3 * k + 1] -= d[1] / 2
        rhs.append(r)
    T = np.linalg.solve(np.array(rows), np.array(rhs))
    for i, j in ((0, 1), (1, 2), (2, 0)):                  # mid-side rule: exact for the quadratic energy
        p = 0.5 * (xy[i] + xy[j]) - xc
        P = np.array([[1.0, 0.0, -p[1]], [0.0, 1.0, p[0]]])
        K += A / 3 * (P @ T).T @ Ds @ (P @ T)
    return K


def test_triangle_without_bubble_and_weight_is_mitc3():
    """With the bubble rotation fixed at 0 (linear rotations) and the full weight on the linear shear part
    (alpha = 0) the triangle is the MITC3 plate (Lee & Bathe 2004)."""
    for t in (0.01, 0.5):
        _, Db, Ds = tshell.plate_rigidities(1e6, 0.25, t)
        op = tshell.tri_plate_operators(TRI, Db, Ds, t, alpha=0.0)
        K3 = _old_mitc3(TRI, Db, Ds)
        assert np.allclose(op["K"][:9, :9], K3, rtol=0, atol=1e-13 * np.abs(K3).max())
        assert float(op["weight"]) == 1.0
        op = tshell.tri_plate_operators(TRI, Db, Ds, t)
        h = max(np.linalg.norm(TRI[i] - TRI[j]) for i, j in ((0, 1), (1, 2), (2, 0)))
        assert float(op["weight"]) == pytest.approx(t * t / (t * t + tshell.TRI_STAB_ALPHA * h * h), rel=1e-14)


def test_triangle_bubble_decouples_from_constant_curvature():
    """int grad f4 dA = 0: no bending coupling between the bubble and the linear rotations, and the bubble
    enters the shear with its element mean 27/60."""
    _, Db, Ds = tshell.plate_rigidities(1e6, 0.25, 0.1)
    Ds0 = 0.0 * Ds
    op = tshell.tri_plate_operators(TRI, Db, Ds0, 0.1)
    assert np.abs(op["K"][:9, 9:]).max() < 1e-12 * np.abs(op["K"]).max()
    assert np.allclose(tshell.tri_plate_operators(TRI, Db, Ds, 0.1)["Bs"][:, 9:], 27 / 60 * np.eye(2))


@pytest.mark.parametrize("t", [0.01, 0.5, 2.0])
def test_triangle_has_no_near_zero_energy_mode(t):
    """Six rigid-body modes; the next eigenvalue is the drilling penalty level, not a near-mechanism (the
    MITC3+ tying d = 1e-4 would leave one at ~3e-9 of the maximum)."""
    K, _ = tshell.matrices(QUAD[:3], MAT, thick=t)
    w = np.linalg.eigvalsh(0.5 * (K.real + K.real.T))
    assert np.all(np.abs(w[:6]) < 1e-12 * w.max()) and w[6] > 1e-7 * w.max()


def _mh_triangle_patch():
    out = [(0.0, 0.0), (0.24, 0.0), (0.24, 0.12), (0.0, 0.12)]
    inn = [(0.04, 0.02), (0.18, 0.03), (0.16, 0.08), (0.08, 0.08)]
    xy = np.array(out + inn)
    tris = []
    for c in ((0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7)):
        tris += [(c[0], c[1], c[2]), (c[0], c[2], c[3])]
    return xy, tris


@pytest.mark.parametrize("t", [0.001, 0.01, 0.2])
def test_triangle_constant_shear_equilibrium_patch(t):
    """Constant transverse shear Q = D_s gamma balanced by the uniform distributed couple m = Q (moment
    equilibrium Q = div M + m with M = 0): distorted MacNeal-Harder patch of triangles, boundary nodes
    prescribed, interior free, work-consistent loads (corner rotations: A/3 m, bubble: 27/60 A m,
    condensed): the solution and Q are exact."""
    xy, tris = _mh_triangle_patch()
    _, Db, Ds = tshell.plate_rigidities(1e6, 0.25, t)
    gam = np.array([1e-3, -2e-3])
    m = Ds @ gam
    K = np.zeros((24, 24))
    F = np.zeros(24)
    keep = []
    for tri in tris:
        X = xy[list(tri)]
        op = tshell.tri_plate_operators(X, Db, Ds, t)
        A = shell.polygon_area(X)
        f = np.zeros(11)
        f[2:9:3], f[1:9:3] = A / 3 * m[0], -A / 3 * m[1]          # theta_y = beta_x, theta_x = -beta_y
        f[9:] = 27 / 60 * A * m
        Kf = op["K"]
        Xc = np.linalg.solve(Kf[9:, 9:], np.column_stack([Kf[9:, :9], f[9:]]))
        dofs = np.array([3 * k + d for k in tri for d in range(3)])
        K[np.ix_(dofs, dofs)] += Kf[:9, :9] - Kf[:9, 9:] @ Xc[:, :9]
        F[dofs] += f[:9] - Kf[:9, 9:] @ Xc[:, 9]
        keep.append((dofs, op, Xc))
    u_ex = np.zeros(24)
    u_ex[0::3] = xy @ gam
    free = np.repeat(np.arange(8) >= 4, 3)
    u = u_ex.copy()
    u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free] - K[np.ix_(free, ~free)] @ u_ex[~free])
    assert np.abs(u - u_ex).max() < 1e-10 * np.abs(u_ex).max()
    for dofs, op, Xc in keep:
        bubble = Xc[:, 9] - Xc[:, :9] @ u[dofs]
        Q = Ds @ (op["Bs"][:, :9] @ u[dofs] + op["Bs"][:, 9:] @ bubble)
        assert np.allclose(Q, m, rtol=1e-10, atol=0)


@pytest.mark.parametrize("t", [0.001, 0.01, 0.2])
def test_triangle_bending_patch(t):
    """Constant curvature (w quadratic, beta = -grad w) on the distorted triangle patch, interior nodes free
    (condensed elements): exact nodal values and moments, zero shear."""
    xy, tris = _mh_triangle_patch()
    _, Db, Ds = tshell.plate_rigidities(1e6, 0.25, t)
    k1, k2, k12 = 1e-3, 2e-3, 0.5e-3
    w = lambda x, y: 0.5 * k1 * x * x + 0.5 * k2 * y * y + k12 * x * y
    u_ex = np.concatenate([[w(x, y), k2 * y + k12 * x, -(k1 * x + k12 * y)] for x, y in xy])
    K = np.zeros((24, 24))
    ops = []
    for tri in tris:
        Kc, Bk, Bs = tshell._tri_plate(xy[list(tri)], Db, Ds, t)
        dofs = np.array([3 * k + d for k in tri for d in range(3)])
        K[np.ix_(dofs, dofs)] += Kc
        ops.append((dofs, Bk, Bs))
    free = np.repeat(np.arange(8) >= 4, 3)
    u = u_ex.copy()
    u[free] = np.linalg.solve(K[np.ix_(free, free)], -K[np.ix_(free, ~free)] @ u_ex[~free])
    assert np.abs(u - u_ex).max() < 1e-9 * np.abs(u_ex).max()
    m = Db @ np.array([-k1, -k2, -2 * k12])
    for dofs, Bk, Bs in ops:
        assert np.allclose(Db @ Bk @ u[dofs], m, rtol=1e-9, atol=0)
        assert np.abs(Ds @ Bs @ u[dofs]).max() < 1e-8 * np.abs(m).max() / 0.24


def _ss_triangle_plate(n, t, diag="right", a=10.0, E=200e9, nu=0.3):
    """Hard simply supported square plate on an n x n mesh split into triangles, unit uniform load: centre w."""
    mat = material_from_E_nu(E, nu, 8000.0)
    nid = lambda i, j: 1 + i + j * (n + 1)
    node_xyz = {nid(i, j): (a * i / n, a * j / n, 0.0) for j in range(n + 1) for i in range(n + 1)}
    recs = []
    for j in range(n):
        for i in range(n):
            c = (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1))
            tri = ([(c[0], c[1], c[2]), (c[0], c[2], c[3])] if diag == "right" or (i + j) % 2 == 0
                   else [(c[0], c[1], c[3]), (c[1], c[2], c[3])])
            recs += [ElemRecord(1, len(recs) + k + 1, 5, nodes, False, mat, dict(thick=t)) for k, nodes in enumerate(tri)]
    ids = sorted(node_xyz)
    edge = lambda v: abs(v) < 1e-9 or abs(v - a) < 1e-9
    fix = np.array([[1, 1, int(edge(node_xyz[k][0]) or edge(node_xyz[k][1])), int(edge(node_xyz[k][0])),
                     int(edge(node_xyz[k][1])), 0] for k in ids])
    dm = build_dofmap(ids, fix, recs)
    am = assemble(node_xyz, recs, dm, undamped=True)
    h = a / n
    F = np.zeros(dm.neq)
    for k, (x, y, _) in node_xyz.items():
        e = dm.eq(k, 3)
        if e >= 0:
            F[e] += (h if 0 < x < a else h / 2) * (h if 0 < y < a else h / 2)
    import scipy.sparse.linalg as spla
    u = spla.spsolve(am.Ks.real.tocsc(), F)
    centre = nid(n // 2, n // 2)
    D = E * t ** 3 / (12 * (1 - nu * nu))
    G = E / (2 * (1 + nu))
    ref = sum(16 / (math.pi ** 2 * p * q) * (1 / (D * k2 * k2) + 1 / (5 / 6 * G * t * k2)) * math.sin(p * math.pi / 2)
              * math.sin(q * math.pi / 2) for p in range(1, 120, 2) for q in range(1, 120, 2)
              for k2 in [((p * math.pi / a) ** 2 + (q * math.pi / a) ** 2)])
    return float(u[dm.eq(centre, 3)]) / ref


@pytest.mark.parametrize("t", [0.01, 1.0])
def test_triangles_do_not_lock(t):
    """Simply supported plate, t/a = 1/1000 and 1/10, 16 x 16 cells with every diagonal the same way (the plain
    MITC3 triangle gives 0.457 for t/a = 1/1000): within 1.5 % of Navier/Mindlin, and the 'cross' pattern too."""
    assert abs(_ss_triangle_plate(16, t) - 1.0) < 0.015
    assert abs(_ss_triangle_plate(16, t, "cross") - 1.0) < 0.015


def test_mitc3_reproduces_constant_shear_and_linear_edge_tying():
    xy = np.array([[0.0, 0.0], [1.3, 0.2], [0.4, 0.9]])
    T, _ = tshell._mitc3_coefficients(xy)
    g = np.array([3.0, -2.0])
    u = np.zeros(9)
    u[0::3] = xy @ g                                    # w linear, theta = 0
    assert np.allclose(T @ u, [g[0], g[1], 0.0])
    # rigid rotation of the plate (beta = -grad w, gamma = 0) gives no shear
    u = np.zeros(9)
    u[0::3] = xy @ g
    u[1::3] = g[1]                                      # theta_x = w,y
    u[2::3] = -g[0]                                     # theta_y = -w,x
    assert np.allclose(T @ u, 0.0, atol=1e-12)


# ---------------------------------------------------------------------------------------- faces
def test_face_stresses_membrane_and_bending():
    t, E, nu = 0.5, 3.0e7, 0.2
    out = tshell.face_stresses([100.0, 0, 0, 0, 0, 0, 0, 0], t, E, nu)
    assert out["rows"]["++"][0] == pytest.approx(200.0) and out["rows"]["+-"][0] == pytest.approx(200.0)
    assert out["rows"]["--"][0] == pytest.approx(-200.0)
    out = tshell.face_stresses([0, 0, 0, 0, 0, 10.0, 0, 0], t, E, nu)
    assert out["rows"]["++"][0] == pytest.approx(6 * 10.0 / t ** 2)
    assert out["rows"]["+-"][0] == pytest.approx(-6 * 10.0 / t ** 2)
    assert out["rows"]["-+"][0] == pytest.approx(6 * 10.0 / t ** 2)
    v = tshell.face_stresses([10, 20, 5, 6, 8, 1, 2, 0.5], t, E, nu)
    for name, sn, sm in tshell.FACE_PERMUTATIONS:
        sxx, syy, txy, s1, s2, exx, eyy, gxy, e1, e2 = v["rows"][name]
        assert sxx == pytest.approx(sn * 10 / t + sm * 6 / t ** 2)
        assert txy == pytest.approx(5 / t + 6 * 0.5 / t ** 2)
        assert s1 + s2 == pytest.approx(sxx + syy) and s1 * s2 == pytest.approx(sxx * syy - txy ** 2)
        assert exx == pytest.approx((sxx - nu * syy) / E) and gxy == pytest.approx(txy * 2 * (1 + nu) / E)
        assert e1 + e2 == pytest.approx(exx + eyy)
    assert v["TXZ"] == pytest.approx(1.5 * 6 / t) and v["TYZ"] == pytest.approx(1.5 * 8 / t)
    # maxima are magnitudes: the sign of the input does not matter
    w = tshell.face_stresses([-10, -20, -5, -6, -8, -1, -2, -0.5], t, E, nu)
    assert np.allclose(w["rows"]["++"], v["rows"]["++"])
    with pytest.raises(ElementError):
        tshell.face_stresses([1, 2, 3], t, E, nu)


def test_mindlin_reference_reproduces_test21():
    from sassi.verify.problems.vp_tshell import TEST21_REF, mindlin_ss_frequencies
    f = mindlin_ss_frequencies(10.0, 10.0, 1.0, 200e9, 0.3, 8000.0)
    assert np.allclose(f, TEST21_REF, rtol=2e-4)
    # thin limit = Kirchhoff  f = pi/2 (2/a^2) sqrt(D / (rho t))
    t = 0.01
    D = 200e9 * t ** 3 / (12 * 0.91)
    kir = math.pi / 2 * 2 / 100.0 * math.sqrt(D / (8000.0 * t))
    assert mindlin_ss_frequencies(10.0, 10.0, t, 200e9, 0.3, 8000.0, nmodes=1)[0] == pytest.approx(kir, rel=1e-5)


# ---------------------------------------------------------------------------------------- eigenvalues
def test_natural_frequencies_with_the_singular_rotary_inertia_of_an_oblique_plate():
    """The rotary inertia rho t^3/12 A/n (I - n n^T) of an oblique facet is singular without a zero row:
    the dense eigen-solution condenses its massless directions (the normal rotations) and gives the
    frequencies of the same plate in the X-Y plane; the mode shapes are M-orthonormal eigenvectors."""
    from sassi.elements import natural_frequencies
    mat = material_from_E_nu(200e9, 0.3, 8000.0)
    n = 3
    nid = lambda i, j: 1 + i + j * (n + 1)
    xyz = {nid(i, j): (2.0 * i, 2.0 * j, 0.0) for j in range(n + 1) for i in range(n + 1)}
    recs = [ElemRecord(1, k + 1, 5, (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)), False, mat,
                       dict(thick=0.8, eint=1)) for k, (i, j) in enumerate((i, j) for j in range(n) for i in range(n))]
    R = _rot(13)
    out = []
    for X in (xyz, {k: tuple(R @ np.asarray(p)) for k, p in xyz.items()}):
        ids = sorted(X)
        dm = build_dofmap(ids, None, recs)
        am = assemble(X, recs, dm, undamped=True)
        out.append((am, natural_frequencies(am.Ks, am.Ms, 12, return_modes=True)))
    (f0, _), (am, (f1, phi)) = out[0][1], out[1]
    assert np.allclose(f1[6:], f0[6:], rtol=1e-8) and np.all(f1[:6] < 1e-4 * f1[6])
    K, M = am.Ks.real.toarray(), am.Ms.toarray()
    assert np.allclose(phi.T @ M @ phi, np.eye(12), atol=1e-9)
    w2 = (2 * np.pi * f1[6:]) ** 2
    res = K @ phi[:, 6:] - M @ phi[:, 6:] * w2
    assert np.abs(res).max() < 1e-7 * np.abs(K @ phi[:, 6:]).max()
