"""Unit tests of the SYMM helpers (sassi/core/symmetry.py; spec 07 9.2.39, D-ANL-12): plane parsing and the
manual's Errors 2-5 / EDU-26 geometry rules, the symmetry boundary conditions, image sums of the reduced
flexibility, the removal of constrained DOFs before the inversion and the free-field symmetry test."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.core import symmetry as SYM

TOL = 1e-6


def _rows(*specs):
    out = []
    for s in specs:
        s = tuple(s) + (0,) * (5 - len(s))
        out.append(dict(no=s[0], type=s[1], n1=s[2], n2=s[3], n3=s[4]))
    return out


NODES = {1: (0.0, 0.0, 0.0), 2: (0.0, 5.0, 0.0), 3: (0.0, 0.0, 9.0), 4: (5.0, 0.0, 0.0), 5: (5.0, 5.0, 0.0),
         6: (0.0, 0.0, -3.0), 7: (3.0, 4.0, 0.0), 8: (0.0, 10.0, 0.0), 9: (2.0, 0.0, 0.0), 10: (5.0, 0.0, -3.0)}


# ---------------------------------------------------------------------------------------------- parsing
def test_planes_parallel_to_yz_and_xz():
    planes, errs, warns = SYM.parse_planes(_rows((1, 1, 1, 2, 3), (2, 0, 1, 4, 3)), NODES, 2, TOL)
    assert errs == [] and warns == []
    p1, p2 = planes
    assert (p1.no, p1.type, p1.axis, p1.coord) == (1, 1, 0, 0.0)          # x = 0 (normal X), antisymmetric
    assert (p2.no, p2.type, p2.axis, p2.coord) == (2, 0, 1, 0.0)          # y = 0 (normal Y), symmetric
    assert "x = 0" in p1.describe() and "antisymmetry" in p1.describe() and "line" in p1.describe(dim=1)


def test_two_nodes_define_a_vertical_plane_in_3d():
    planes, errs, _ = SYM.parse_planes(_rows((1, 0, 1, 2)), NODES, 2, TOL)     # (0,0,0)-(0,5,0): x = 0
    assert errs == [] and planes[0].axis == 0
    planes, errs, _ = SYM.parse_planes(_rows((1, 0, 1, 3)), NODES, 2, TOL)     # vertical line: ambiguous
    assert planes == [] and any("one vertical line" in e for e in errs)


@pytest.mark.parametrize("row,code", [((1, 0, 1, 2, 99), 2), ((1, 0, -1, 2, 3), 2), ((1, 0, 1), 3), ((1, 0, 1, 1), 4),
                                      ((1, 0, 1, 2, 8), 5)])
def test_manual_errors_2_to_5(row, code):
    planes, errs, _ = SYM.parse_planes(_rows(row), NODES, 2, TOL)
    assert planes == [] and errs and errs[0].startswith(f"Error {code}:")


def test_geometry_rules_edu26():
    # oblique vertical plane (through (0,0,0), (3,4,0), (0,0,9)): not parallel to XZ or YZ
    _, errs, _ = SYM.parse_planes(_rows((1, 0, 1, 7, 3)), NODES, 2, TOL)
    assert any("not parallel to the XZ or YZ plane" in e for e in errs)
    # horizontal plane
    _, errs, _ = SYM.parse_planes(_rows((1, 0, 1, 4, 5)), NODES, 2, TOL)
    assert any("not parallel" in e for e in errs)
    # two parallel planes
    _, errs, _ = SYM.parse_planes(_rows((1, 0, 1, 2, 3), (2, 1, 1, 8, 6)), NODES, 2, TOL)
    assert any("parallel" in e and "orthogonal" in e for e in errs)
    # plane number, type, duplicates
    _, errs, _ = SYM.parse_planes(_rows((3, 0, 1, 2, 3)), NODES, 2, TOL)
    assert any("must be 1 or 2" in e for e in errs)
    _, errs, _ = SYM.parse_planes(_rows((1, 2, 1, 2, 3)), NODES, 2, TOL)
    assert any("<type> = 2" in e for e in errs)
    _, errs, _ = SYM.parse_planes(_rows((1, 0, 1, 2, 3), (1, 0, 1, 4, 3)), NODES, 2, TOL)
    assert any("defined twice" in e for e in errs)


def test_2d_line_parallel_to_z_and_reset_rows():
    planes, errs, _ = SYM.parse_planes(_rows((1, 1, 1, 6)), NODES, 1, TOL)     # (0,0,0)-(0,0,-3)
    assert errs == [] and planes[0].axis == 0 and planes[0].coord == 0.0
    _, errs, _ = SYM.parse_planes(_rows((1, 1, 1, 4)), NODES, 1, TOL)          # horizontal line
    assert any("not parallel to the Z axis" in e for e in errs)
    _, errs, _ = SYM.parse_planes(_rows((1, 1, 1, 6), (2, 0, 4, 10)), NODES, 1, TOL)
    assert any("at most one SYMM line" in e for e in errs)
    planes, errs, warns = SYM.parse_planes(_rows((1, 0, 0, 0)), NODES, 2, TOL)
    assert planes == [] and errs == [] and any("resets" in w for w in warns)


# ------------------------------------------------------------------------------------ boundary conditions
@pytest.mark.parametrize("stype,axis,dofs", [(0, 0, (1, 5, 6)), (0, 1, (2, 4, 6)), (1, 0, (2, 3, 4)), (1, 1, (1, 3, 5))])
def test_fixed_dofs_d_anl_12(stype, axis, dofs):
    """Symmetric: normal translation + in-plane rotations; antisymmetric: in-plane translations + normal
    rotation (D-ANL-12)."""
    assert SYM.SymmetryPlane(1, stype, axis, 0.0).fixed_dofs == dofs


def test_symmetry_fixity_on_two_planes():
    planes = [SYM.SymmetryPlane(1, 1, 0, 0.0), SYM.SymmetryPlane(2, 0, 1, 0.0)]       # X input: x anti, y sym
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 2.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0]])
    fix = SYM.symmetry_fixity(planes, xyz, TOL)
    assert fix[0].tolist() == [False, True, True, True, False, True]    # both planes: only UX, ROTY free
    assert fix[1].tolist() == [False, True, True, True, False, False]   # x = 0: UY, UZ, ROTX
    assert fix[2].tolist() == [False, True, False, True, False, True]   # y = 0: UY, ROTX, ROTZ
    assert not fix[3].any()
    assert SYM.constrained_translations(planes, xyz, TOL)[0].tolist() == [False, True, True]


def test_side_violations_and_mirror():
    p = SYM.SymmetryPlane(1, 0, 0, 1.0)
    assert SYM.side_violations([p], np.array([[1.0, 0, 0], [2.0, 0, 0]]), TOL) == []
    assert SYM.side_violations([p], np.array([[0.5, 0, 0], [2.0, 0, 0]]), TOL) == [(1, 1, 1)]
    np.testing.assert_allclose(p.mirror(np.array([[3.0, 2.0, -1.0]])), [[-1.0, 2.0, -1.0]])
    terms = SYM.image_terms([p, SYM.SymmetryPlane(2, 1, 1, 0.0)])
    assert [len(S) for S, _, _ in terms] == [0, 1, 1, 2]
    assert [s for _, s, _ in terms] == [1.0, 1.0, -1.0, -1.0]
    np.testing.assert_array_equal(terms[3][2], [-1.0, -1.0, 1.0])


def test_required_type_for_vertical_input():
    assert SYM.required_type(0, 0) == SYM.ANTISYMMETRIC and SYM.required_type(1, 0) == SYM.SYMMETRIC
    assert SYM.required_type(0, 2) == SYM.SYMMETRIC and SYM.required_type(1, 1) == SYM.ANTISYMMETRIC


def test_planes_array_round_trip():
    planes = [SYM.SymmetryPlane(1, 1, 0, 2.5), SYM.SymmetryPlane(2, 0, 1, -1.0)]
    back = SYM.planes_from_array(SYM.planes_array(planes))
    assert [(p.no, p.type, p.axis, p.coord) for p in back] == [(1, 1, 0, 2.5), (2, 0, 1, -1.0)]
    assert SYM.planes_array([]).shape == (0, 4)


# ---------------------------------------------------------------------------------- image sums
def kelvin_block(obs: np.ndarray, load: np.ndarray, a: float = 0.7, b: float = 0.3) -> np.ndarray:
    """A reciprocal, reflection-invariant model flexibility (Kelvin-like, regularised): F_ij = f(r) I +
    g(r) d d^T / r^2 with f = a / (1 + r), g = b / (1 + r)."""
    no, nl = obs.shape[0], load.shape[0]
    F = np.zeros((no, 3, nl, 3), complex)
    for i in range(no):
        for j in range(nl):
            d = obs[i] - load[j]
            r = float(np.linalg.norm(d))
            dd = np.outer(d, d) / r ** 2 if r > 0 else np.zeros((3, 3))
            F[i, :, j, :] = (a / (1 + r)) * np.eye(3) + (b / (1 + r)) * dd + 0.05j * np.eye(3)
    return F.reshape(3 * no, 3 * nl)


def test_reduced_flexibility_zero_rows_and_symmetry():
    """F_red(i, j) = sum_S s_S F(i, j_S) P_S is symmetric and its rows and columns vanish at the
    translations the symmetry constrains on the plane nodes (zero force, zero displacement)."""
    planes = [SYM.SymmetryPlane(1, 1, 0, 0.0), SYM.SymmetryPlane(2, 0, 1, 0.0)]
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 1.5, 0.0], [1.0, 0.0, -1.0], [2.0, 1.0, 0.0], [1.0, 2.5, -2.0]])
    Fr = SYM.reduced_flexibility(lambda L: kelvin_block(xyz, L), xyz, planes)
    np.testing.assert_allclose(Fr, Fr.T, atol=1e-13)
    cons = SYM.constrained_translations(planes, xyz, TOL).reshape(-1)
    assert cons.sum() == 2 + 2 + 1                          # node 0 (UY, UZ), node 1 (UY, UZ), node 2 (UY)
    assert np.abs(Fr[cons]).max() < 1e-13 and np.abs(Fr[:, cons]).max() < 1e-13


def test_reduced_flexibility_reproduces_the_full_set():
    """Image sum vs the full point set: for a symmetric load pattern f(j') = s P f(j) the displacements of
    the reduced nodes from the full flexibility equal F_red g (g = full force, halved on the plane)."""
    p = SYM.SymmetryPlane(1, 0, 0, 0.0)                     # symmetric about x = 0
    red = np.array([[0.0, 0.0, 0.0], [1.0, 0.5, 0.0], [2.0, -0.5, -1.0]])
    full = np.vstack([red, p.mirror(red[1:])])              # node 0 lies on the plane (no image)
    rng = np.random.default_rng(3)
    g = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    g[0, 0] = 0.0                                           # no normal force on the plane node
    f_full = np.vstack([2.0 * g[0], g[1], g[2], g[1:] * p.reflection[None, :]])
    u_full = (kelvin_block(full, full) @ f_full.reshape(-1)).reshape(-1, 3)[:3]
    Fr = SYM.reduced_flexibility(lambda L: kelvin_block(red, L), red, [p])
    np.testing.assert_allclose(Fr @ g.reshape(-1), u_full.reshape(-1), rtol=1e-13, atol=1e-13)


def test_reduced_impedance_removes_dofs_before_inversion():
    xyz = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.5, 0.5, 0.0]])
    F = kelvin_block(xyz, xyz)
    keep = SYM.keep_mask(1, None, 3)                         # 2D: UX, UZ only
    assert keep.tolist() == [True, False, True] * 3
    X, rc = SYM.reduced_impedance(F.copy(), keep)
    k = np.flatnonzero(keep)
    np.testing.assert_allclose(X[np.ix_(k, k)], np.linalg.inv(0.5 * (F + F.T)[np.ix_(k, k)]), rtol=1e-10)
    assert np.all(X[~keep] == 0) and np.all(X[:, ~keep] == 0) and rc > 0
    cons = np.zeros((3, 3), bool)
    cons[0, 1:] = True
    assert SYM.keep_mask(2, cons, 3).tolist() == [True, False, False] + [True] * 6
    X0, _ = SYM.reduced_impedance(F.copy(), np.ones(9, bool))
    np.testing.assert_allclose(X0, np.linalg.inv(0.5 * (F + F.T)), rtol=1e-10)


def test_field_symmetry_error():
    """A uniform horizontal X field is antisymmetric about x = const and symmetric about y = const."""
    U = np.tile([1.0 + 0.2j, 0.0, 0.0], (4, 1))
    planes = [SYM.SymmetryPlane(1, 1, 0, 0.0), SYM.SymmetryPlane(2, 0, 1, 0.0)]
    (s1, a1), (s2, a2) = SYM.field_symmetry_error(planes, U, [U, U])
    assert a1 == 0.0 and s1 == pytest.approx(2.0) and s2 == 0.0 and a2 == pytest.approx(2.0)
