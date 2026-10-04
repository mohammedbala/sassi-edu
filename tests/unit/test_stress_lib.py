"""Unit tests of sassi.core.stress_lib (requirements 4.10, D-FIL-06, D-STR-07/08/09/12)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.core import stress_lib as SL


# ---------------------------------------------------------------------------------------------
# STF kernels
# ---------------------------------------------------------------------------------------------
def test_element_dof_tf_masks_fixed_dofs():
    H = np.arange(12, dtype=complex).reshape(2, 6) + 1j
    eq = np.array([[0, -1, 5], [2, 2, -1]])
    U = SL.element_dof_tf(H, eq)
    assert U.shape == (2, 2, 3)
    assert np.array_equal(U[:, 0, 0], H[:, 0]) and np.all(U[:, 0, 1] == 0) and np.all(U[:, 1, 2] == 0)
    assert np.array_equal(U[:, 1, 1], H[:, 2])


def test_element_stf_matches_loop_and_rigid_subtraction_is_exact():
    rng = np.random.default_rng(1)
    nE, nc, nd, nF, neq = 3, 4, 6, 5, 9
    S = rng.standard_normal((nE, nc, nd)) + 1j * rng.standard_normal((nE, nc, nd))
    # make S annihilate a uniform translation along DOF label 1 of a 2-node, 3-DOF element
    rig = SL.rigid_body_dofs((1, 2, 3), 2, cm=0)
    S -= np.einsum("ecd,d->ec", S, rig)[:, :, None] * rig[None, None, :] / rig.dot(rig)
    eq = rng.integers(-1, neq, size=(nE, nd))
    H = rng.standard_normal((nF, neq)) + 1j * rng.standard_normal((nF, neq))
    ref = np.zeros((nF, nE, nc), complex)
    for f in range(nF):
        for e in range(nE):
            u = np.array([H[f, k] if k >= 0 else 0.0 for k in eq[e]])
            ref[f, e] = S[e] @ u
    assert np.allclose(SL.element_stf(S, eq, H), ref, rtol=1e-13, atol=1e-13)
    assert np.allclose(SL.element_stf(S, eq, H, rig), ref, rtol=1e-12, atol=1e-12)


def test_rigid_body_dofs():
    assert np.array_equal(SL.rigid_body_dofs((1, 2, 3), 2, cm=0), [1, 0, 0, 1, 0, 0])
    assert np.array_equal(SL.rigid_body_dofs((1, 2, 3, 4, 5, 6), 1, cm=2), [0, 0, 1, 0, 0, 0])
    r = SL.rigid_body_dofs((1, 3), 2, cm=0, ang_deg=30.0)            # PLANE: UX, UZ
    assert np.allclose(r, [np.cos(np.pi / 6), 0, np.cos(np.pi / 6), 0])
    r = SL.rigid_body_dofs((1, 2, 3), 1, cm=1, ang_deg=90.0)          # y' = -x
    assert np.allclose(r, [-1, 0, 0])


# ---------------------------------------------------------------------------------------------
# derived quantities
# ---------------------------------------------------------------------------------------------
def test_octahedral_stress_reference_states():
    tau, sig = 3.7, -2.2
    assert SL.octahedral_shear_stress(0, 0, 0, tau, 0, 0) == pytest.approx(0.816497 * tau, rel=1e-6)
    assert SL.octahedral_shear_stress(0, 0, 0, 0, tau, 0) == pytest.approx(np.sqrt(6) / 3 * tau, rel=1e-14)
    assert SL.octahedral_shear_stress(sig, 0, 0, 0, 0, 0) == pytest.approx(np.sqrt(2) / 3 * abs(sig), rel=1e-14)
    assert SL.octahedral_shear_stress(5.0, 5.0, 5.0, 0, 0, 0) == 0.0          # hydrostatic
    assert SL.OCT_STRESS_PURE_SHEAR == pytest.approx(0.8165, abs=5e-5)
    # invariance: rotate a stress tensor, tau_oct unchanged (it is sqrt(2/3 J2))
    s = np.array([[1.0, 0.4, -0.3], [0.4, -2.0, 0.7], [-0.3, 0.7, 0.5]])
    th = 0.7
    R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    s2 = R @ s @ R.T
    a = SL.octahedral_shear_stress(s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[0, 2], s[1, 2])
    b = SL.octahedral_shear_stress(s2[0, 0], s2[1, 1], s2[2, 2], s2[0, 1], s2[0, 2], s2[1, 2])
    dev = s - np.trace(s) / 3 * np.eye(3)
    assert a == pytest.approx(b, rel=1e-13) and a == pytest.approx(np.sqrt(2 / 3 * 0.5 * np.sum(dev * dev)), rel=1e-13)


def test_octahedral_strain_and_2d_max_shear():
    g = 1e-3
    assert SL.octahedral_shear_strain(0, 0, 0, g, 0, 0) == pytest.approx(np.sqrt(6) / 3 * g, rel=1e-14)
    e = 2e-4
    assert SL.octahedral_shear_strain(e, -0.25 * e, -0.25 * e, 0, 0, 0) == pytest.approx(
        2 / 3 * np.sqrt(2 * 1.25 ** 2) * e, rel=1e-14)
    assert SL.max_shear_strain_2d(3e-4, -1e-4, 3e-4) == pytest.approx(5e-4, rel=1e-14)
    assert SL.strain_component_names("SOLID", ["SXX", "SXY", "SOCT"]) == ["EXX", "EXY", "EOCT"]
    assert SL.strain_component_names("PLANE", ["SXX", "SZZ", "TXZ"]) == ["EXX", "EZZ", "EXZ"]


# ---------------------------------------------------------------------------------------------
# ELEMENT_CENTER layout
# ---------------------------------------------------------------------------------------------
MANUAL_EXAMPLE = """3
SOLID           1    1       3
SHELL           2    1       1
SOLID           3    1       2
SOLID           1    1
              1 0.69504 0.61290 0.93326 0.21454 1.36011 0.45008
              2 0.82394 0.70086 0.68225 0.20217 0.65360 0.34301
              3 1.63296 1.09535 1.41437 0.49395 1.40079 0.39915
              SHELL    2    1
              1 6.98477 12.51727 9.28106 0.87311 0.58065 0.44423
              SOLID    3    1
              1 1.07909 1.12969 1.91468 0.10359 1.03872 0.25992
              2 1.0 2.0 3.0 4.0 5.0 6.0
"""


def test_read_manual_example_layout(tmp_path):
    p = tmp_path / "ELEMENT_CENTER_ABS_MAX_STRESSES.TXT"
    p.write_text(MANUAL_EXAMPLE)
    blocks = SL.read_element_center(p)
    assert [(b.etype, b.group, b.ordered, len(b.elements)) for b in blocks] == [
        ("SOLID", 1, 1, 3), ("SHELL", 2, 1, 1), ("SOLID", 3, 1, 2)]
    assert blocks[1].values[0, 1] == pytest.approx(12.51727)
    assert blocks[0].values.shape == (3, 6)


def test_element_center_round_trip_and_ordered_groups(tmp_path):
    og = SL.ordered_group_numbers([(3, "SOLID"), (1, "SOLID"), (2, "SHELL"), (7, "SPRING"), (5, "SHELL")])
    assert og == {1: 1, 2: 1, 3: 2, 5: 2, 7: 1}                     # D-FIL-06: order among the same type
    blocks = [SL.CenterBlock("SOLID", 1, 1, [1, 2], np.array([[1.5, -2, 3, 4, 5, 6], [0, 0, 0, 0, 0, 1e-12]])),
              SL.CenterBlock("SPRING", 7, 1, [10], np.arange(6.0)[None, :]),
              SL.CenterBlock("SOLID", 3, 2, [4], [[7.0]])]
    p = SL.write_element_center(tmp_path / "x.ess", blocks)
    lines = p.read_text().splitlines()
    assert lines[0] == "3" and lines[1].split() == ["SOLID", "1", "1", "2"] and lines[4].split() == ["SOLID", "1", "1"]
    back = SL.read_element_center(p)
    for a, b in zip(blocks, back):
        assert (a.etype, a.group, a.ordered) == (b.etype, b.group, b.ordered)
        assert np.array_equal(a.elements, b.elements)
        assert np.allclose(a.values, b.values, rtol=1e-9, atol=0)
    with pytest.raises(ValueError):
        (tmp_path / "bad").write_text("2\nSOLID 1 1 1\n")
        SL.read_element_center(tmp_path / "bad")


def test_center_columns():
    assert SL.center_columns("PLANE", ["SXX", "SZZ", "TXZ"]) == [0, 1, 2, -1, -1, -1]
    assert SL.center_columns("SOLID", ["SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ", "SOCT"]) == [0, 1, 2, 3, 4, 5]
    assert "BEAMS" not in SL.CENTER_COLUMNS


# ---------------------------------------------------------------------------------------------
# Frames.txt, output steps
# ---------------------------------------------------------------------------------------------
def test_frames_file_manual_example(tmp_path):
    p = tmp_path / "Frames.txt"
    p.write_text("10\n" + "\n".join(str(i) for i in range(1, 11)) + "\n3\n8 9 10\n")
    frames, groups = SL.read_frames_file(p)
    assert frames == list(range(1, 11)) and groups == [8, 9, 10]
    p.write_text("2\n5\n7\n")
    assert SL.read_frames_file(p) == ([5, 7], [])
    p.write_text("4\n1\n2\n")
    with pytest.raises(ValueError):
        SL.read_frames_file(p)


def test_output_steps():
    assert np.array_equal(SL.output_steps(7, 0), np.arange(7))
    assert np.array_equal(SL.output_steps(7, 1), np.arange(7))
    assert np.array_equal(SL.output_steps(7, 3), [0, 3, 6])


# ---------------------------------------------------------------------------------------------
# nodal averages and soil-pressure geometry
# ---------------------------------------------------------------------------------------------
def test_averaging_matrix_plain_mean():
    nodes, A = SL.averaging_matrix([[1, 2, 3, 3, 0], [3, 4, 2]])
    assert np.array_equal(nodes, [1, 2, 3, 4])
    v = A @ np.array([2.0, 6.0])
    assert np.allclose(v, [2.0, 4.0, 4.0, 6.0])


def test_solid_faces_shared_faces_normals():
    hexa = [1, 2, 3, 4, 5, 6, 7, 8]
    faces = SL.solid_faces(hexa)
    assert len(faces) == 6 and all(len(f) == 4 for f in faces)
    prism = [1, 2, 3, 3, 5, 6, 7, 7]                     # nodes 3 = 4 and 7 = 8
    pf = SL.solid_faces(prism)
    assert len(pf) == 5 and sorted(len(f) for f in pf) == [3, 3, 4, 4, 4]
    assert SL.shared_faces(hexa, {2, 3, 6, 7, 99}) == [(2, 3, 7, 6)]
    xyz = np.array([[1.0, 0, 0], [1.0, 1, 0], [1.0, 1, 1], [1.0, 0, 1]])
    n = SL.face_normal(xyz)
    assert np.allclose(np.abs(n), [1, 0, 0])
    sig = np.array([[-3.0, 1.0, 2.0, 0.5, 0.25, 0.1]])
    assert SL.normal_stress(sig, n)[0] == pytest.approx(-3.0)
    m = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)
    assert SL.normal_stress(sig, m)[0] == pytest.approx(0.5 * (-3 + 1) + 0.5)
    with pytest.raises(ValueError):
        SL.face_normal(np.zeros((3, 3)))


def test_average_histories_equals_per_step_product_and_stays_sparse():
    """Review finding (D-GEN-10): the nodal averaging is one sparse product, never a dense nN x nE matrix."""
    rng = np.random.default_rng(3)
    elems = [rng.choice(np.arange(1, 40), size=8, replace=False) for _ in range(25)]
    nodes, A = SL.averaging_matrix(elems)
    E = rng.standard_normal((7, 25, 3))
    N = SL.average_histories(A, E)
    assert N.shape == (7, len(nodes), 3)
    assert np.allclose(N, np.stack([A.toarray() @ E[t] for t in range(7)]), rtol=1e-13, atol=1e-13)
    # 60,000 single-node 'elements': a dense A would hold 3.6e9 entries; the sparse product is immediate
    n = 60000
    nodes, A = SL.averaging_matrix([[i + 1, 0] for i in range(n)])
    E = np.arange(4 * n * 2, dtype=float).reshape(4, n, 2)
    assert np.array_equal(SL.average_histories(A, E), E)
    with pytest.raises(ValueError):
        SL.average_histories(A, E[:, :10, :])


def test_membrane_to_global_is_the_tensor_rotation():
    """sigma_g = Lam^T sigma_l Lam with sigma_l = [[FXX, FXY, 0], [FXY, FYY, 0], [0, 0, 0]] (plane stress in the
    shell plane), returned in the SOLID order SXX SYY SZZ SXY SXZ SYZ; invariants are preserved."""
    rng = np.random.default_rng(5)
    nE = 6
    Lam = np.stack([np.linalg.qr(rng.standard_normal((3, 3)))[0].T for _ in range(nE)])
    F = rng.standard_normal((4, nE, 3))                         # 4 time steps
    G = SL.membrane_to_global(F, Lam)
    assert G.shape == (4, nE, 6)
    for t in range(4):
        for e in range(nE):
            fxx, fyy, fxy = F[t, e]
            sl = np.array([[fxx, fxy, 0.0], [fxy, fyy, 0.0], [0.0, 0.0, 0.0]])
            sg = Lam[e].T @ sl @ Lam[e]
            ref = [sg[0, 0], sg[1, 1], sg[2, 2], sg[0, 1], sg[0, 2], sg[1, 2]]
            assert np.allclose(G[t, e], ref, rtol=1e-13, atol=1e-13)
            assert np.trace(sg) == pytest.approx(fxx + fyy)
    # identity axes: FXX FYY 0 / FXY 0 0
    G = SL.membrane_to_global(np.array([[1.0, 2.0, 3.0]]), np.eye(3)[None])
    assert np.allclose(G, [[1.0, 2.0, 0.0, 3.0, 0.0, 0.0]])
