"""Unit tests of sassi.elements.assemble: DOF map (requirements 4.1 DOF management, D-ELM-09),
sparse assembly of structure / excavated soil (requirements 4.4), nodal masses (UT-15, D-MDL-08)
and the recovery operators (ARCHITECTURE 6.2)."""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from sassi.conventions import cfactor
from sassi.elements import (ElemRecord, ElementError, assemble, beam, build_dofmap, material_from_E_nu,
                            material_from_layer, natural_frequencies, nodal_masses_from_table, solid,
                            unstiffened_rotations)

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")

MAT = material_from_E_nu(1000.0, 0.25, 2.0)
SEC = dict(A=0.02, As2=0.0, As3=0.0, J=3e-5, I2=4e-5, I3=1e-4)


def _brick_nodes(offset=0, z0=0.0):
    xyz = (solid.NODE_NAT + 1.0) / 2.0 + np.array([0, 0, z0])
    return {offset + i + 1: xyz[i] for i in range(8)}


def test_dofmap_union_fix_and_k_node():
    node_xyz = _brick_nodes()
    node_xyz[20] = (0.0, 0.0, 3.0)            # beam end node
    node_xyz[30] = (5.0, 5.0, 5.0)            # K node only
    node_xyz[40] = (9.0, 9.0, 9.0)            # unused node
    recs = [ElemRecord(1, 1, 1, tuple(range(1, 9)), False, MAT, {}),
            ElemRecord(2, 1, 2, (5, 20, 30), False, MAT, dict(section=SEC))]
    ids = sorted(node_xyz)
    fix = np.zeros((len(ids), 6), dtype=int)
    fix[ids.index(1)] = 1
    dm = build_dofmap(ids, fix, recs)
    assert dm.eq(1, 1) == -1                                      # fixed
    assert dm.eq(2, 1) == 0                                       # first free equation, node-major
    assert dm.eq(2, 4) == -1                                      # SOLID-only node: no rotations
    assert all(dm.eq(5, d) >= 0 for d in range(1, 7))             # beam + solid node: 6 DOFs
    assert all(dm.eq(30, d) == -1 for d in range(1, 7))           # K node carries no DOF
    assert all(dm.eq(40, d) == -1 for d in range(1, 7))           # unused node
    assert dm.neq == 6 * 3 + 6 + 6                                # nodes 2,3,4,6,7,8 (x3), 5 and 20 (x6)
    assert np.array_equal(np.sort(dm.eq_table[dm.eq_table >= 0]), np.arange(dm.neq))
    for e in range(dm.neq):
        assert dm.eq(dm.eq_node[e], dm.eq_dof[e]) == e
    assert dm.n_eliminated > 0 and "equations" in dm.summary()
    dm2 = build_dofmap(ids, fix, recs, extra_dofs={1: [4], 2: [4]})
    assert dm2.eq(2, 4) >= 0 and dm2.eq(1, 4) == -1               # extra DOF, but fixed stays fixed


def test_dofmap_count():
    node_xyz = _brick_nodes()
    recs = [ElemRecord(1, 1, 1, tuple(range(1, 9)), False, MAT, {})]
    dm = build_dofmap(sorted(node_xyz), None, recs)
    assert dm.neq == 24 and list(dm.eq_dof[:3]) == [1, 2, 3]


def test_dofmap_errors():
    recs = [ElemRecord(1, 1, 1, tuple(range(1, 9)), False, MAT, {})]
    with pytest.raises(ElementError, match="Error 41"):
        build_dofmap([1, 2, 3], None, recs)
    # TSHELL (type 5) is available since wave 3: it gets six DOFs per node like SHELL
    dm = build_dofmap([1, 2, 3, 4], None, [ElemRecord(1, 1, 5, (1, 2, 3, 4), False, MAT, {})])
    assert dm.neq == 24
    with pytest.raises(ElementError, match="Error 7"):
        build_dofmap(list(range(1, 9)), None, [ElemRecord(1, 1, 1, (1, 2, 3), False, MAT, {})])


def test_structure_and_excavated_matrices_are_separate():
    soil = material_from_layer(400.0, 200.0, 19.62, 0.02, 0.05, 9.81)
    node_xyz = {**_brick_nodes(0, 0.0), **_brick_nodes(8, 1.0)}
    # brick 2 shares the top face of brick 1 (nodes 5..8)
    n2 = (5, 6, 7, 8, 13, 14, 15, 16)
    recs = [ElemRecord(1, 1, 1, tuple(range(1, 9)), True, soil, dict(incompatible=True)),
            ElemRecord(2, 1, 1, n2, False, MAT, dict(incompatible=True))]
    ids = sorted(node_xyz)
    dm = build_dofmap(ids, None, recs)
    am = assemble(node_xyz, recs, dm)
    assert any("D-ELM-02" in m for m in am.messages)               # no incompatible modes in soil
    K1, M1 = solid.matrices(np.array([node_xyz[n] for n in range(1, 9)]), soil, incompatible=False)
    eq1 = dm.eqs(range(1, 9), (1, 2, 3))
    assert np.allclose(am.Ke.toarray()[np.ix_(eq1, eq1)], K1)
    assert np.allclose(am.Me.toarray()[np.ix_(eq1, eq1)], M1)
    K2, _ = solid.matrices(np.array([node_xyz[n] for n in n2]), MAT, incompatible=True)
    eq2 = dm.eqs(n2, (1, 2, 3))
    assert np.allclose(am.Ks.toarray()[np.ix_(eq2, eq2)], K2)
    assert am.Ks.dtype == complex and am.Ms.dtype == float and sp.isspmatrix_csr(am.Ks)
    rec = am.recovery["SOLID"]
    assert rec["S"].shape == (2, 7, 24) and rec["eq"].shape == (2, 24) and rec["B"].shape == (2, 6, 24)
    assert list(rec["group"]) == [1, 2] and list(rec["excavated"]) == [1, 0]


def test_ut15_nodal_mass_units_and_ignored_masses():
    rows_w = [(2, 32.2, 32.2, 32.2, 0, 0, 0, 1)]
    rows_m = [(2, 1.0, 1.0, 1.0, 0, 0, 0, 0)]
    mw = nodal_masses_from_table(rows_w, 32.2)
    mm = nodal_masses_from_table(rows_m, 32.2)
    assert np.array_equal(mw[2], mm[2])
    node_xyz = {1: (0, 0, 0), 2: (1, 0, 0)}
    recs = [ElemRecord(1, 1, 7, (1, 2), False, None, dict(k=(1000, 0, 0, 0, 0, 0), damp=0.0))]
    fix = np.array([[1] * 6, [0, 1, 1, 1, 1, 1]])
    dm = build_dofmap([1, 2], fix, recs)
    am_w = assemble(node_xyz, recs, dm, masses=mw)
    am_m = assemble(node_xyz, recs, dm, masses=mm)
    assert (am_w.Ms != am_m.Ms).nnz == 0
    assert am_w.Ms.toarray()[0, 0] == 1.0
    assert sorted(d for _, d, _ in am_w.ignored_masses) == [2, 3]   # MT on fixed UY, UZ ignored (W5/W6)
    f = natural_frequencies(am_w.Ks, am_w.Ms)
    assert f[0] == pytest.approx(np.sqrt(1000.0) / (2 * np.pi), rel=1e-12)


def test_undamped_assembly_and_complex_spring():
    node_xyz = {1: (0, 0, 0), 2: (1, 0, 0)}
    recs = [ElemRecord(1, 1, 7, (1, 2), False, None, dict(k=(1000, 0, 0, 0, 0, 0), damp=0.1))]
    fix = np.array([[1] * 6, [0, 1, 1, 1, 1, 1]])
    dm = build_dofmap([1, 2], fix, recs)
    am = assemble(node_xyz, recs, dm)
    am0 = assemble(node_xyz, recs, dm, undamped=True)
    assert am.Ks.toarray()[0, 0] == pytest.approx(1000.0 * cfactor(0.1))
    assert am0.Ks.toarray()[0, 0] == pytest.approx(1000.0)
    rec = am.recovery["SPRING"]
    assert rec["S"].shape == (1, 6, 12) and rec["eq"][0, 6] == 0 and rec["eq"][0, 0] == -1


def test_shell_triangles_padded_and_drilling_detection():
    node_xyz = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0), 5: (2, 0, 0)}
    recs = [ElemRecord(3, 1, 3, (1, 2, 3, 4), False, MAT, dict(thick=0.1)),
            ElemRecord(3, 2, 3, (2, 5, 3, 3), False, MAT, dict(thick=0.1)),   # triangle: L = K
            ElemRecord(3, 3, 3, (2, 5, 3, 0), False, MAT, dict(thick=0.1))]   # triangle: L = 0
    ids = sorted(node_xyz)
    dm = build_dofmap(ids, None, recs)
    am = assemble(node_xyz, recs, dm)
    rec = am.recovery["SHELL"]
    assert rec["S"].shape == (3, 6, 24)
    assert np.all(rec["eq"][1:, 18:] == -1) and np.allclose(rec["S"][1:, :, 18:], 0)
    flagged = unstiffened_rotations(am.Ks, dm)
    assert sorted(n for n, _ in flagged) == ids                   # every flat-shell node: drilling about Z
    assert all(abs(abs(ax[2]) - 1.0) < 1e-9 for _, ax in flagged)
    # a beam along Z at node 1 restrains its drilling rotation (torsion)
    node_xyz[6] = (0, 0, 1)
    node_xyz[7] = (1, 0, 1)
    recs.append(ElemRecord(4, 1, 2, (1, 6, 7), False, MAT, dict(section=SEC)))
    ids = sorted(node_xyz)
    dm = build_dofmap(ids, None, recs)
    am = assemble(node_xyz, recs, dm)
    assert 1 not in [n for n, _ in unstiffened_rotations(am.Ks, dm)]


def test_beam_recovery_through_assembly_matches_element():
    node_xyz = {1: (0, 0, 0), 2: (3, 4, 0), 3: (0, 0, 1)}
    recs = [ElemRecord(2, 7, 2, (1, 2, 3), False, MAT, dict(section=SEC, ki="000011"))]
    dm = build_dofmap([1, 2, 3], None, recs)
    am = assemble(node_xyz, recs, dm)
    S = beam.recovery(np.array([node_xyz[n] for n in (1, 2, 3)]), MAT, section=SEC, ki="000011")
    rec = am.recovery["BEAMS"]
    assert np.allclose(rec["S"][0], S) and rec["id"][0] == 7
    assert list(rec["eq"][0]) == list(range(12))


def test_plane_triangle_and_general_assembly():
    node_xyz = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 0, 1), 4: (5, 0, 5)}
    KR = np.zeros((12, 12))
    KR[2, 2] = KR[8, 8] = 10.0
    KR[2, 8] = KR[8, 2] = -10.0
    recs = [ElemRecord(1, 1, 4, (1, 2, 3, 0), False, MAT, {}),
            ElemRecord(9, 1, 9, (3, 4), False, None, dict(KR=KR))]
    dm = build_dofmap([1, 2, 3, 4], None, recs)
    assert dm.eq(1, 2) == -1 and dm.eq(1, 1) >= 0 and dm.eq(1, 3) >= 0     # PLANE: UX, UZ only
    am = assemble(node_xyz, recs, dm)
    assert "PLANE" in am.recovery and "GENERAL" not in am.recovery
    assert am.recovery["PLANE"]["S"].shape == (1, 3, 8)
    assert am.Ks.toarray()[dm.eq(4, 3), dm.eq(3, 3)] == pytest.approx(-10.0)


def test_natural_frequencies_dense_and_sparse_agree():
    # chain of springs with masses
    n = 30
    node_xyz = {i: (float(i), 0.0, 0.0) for i in range(1, n + 2)}
    recs = [ElemRecord(1, i, 7, (i, i + 1), False, None, dict(k=(100.0, 0, 0, 0, 0, 0))) for i in range(1, n + 1)]
    fix = np.ones((n + 1, 6), dtype=int)
    fix[1:, 0] = 0
    dm = build_dofmap(range(1, n + 2), fix, recs)
    am = assemble(node_xyz, recs, dm, masses={i: (1.0, 0, 0, 0, 0, 0) for i in range(2, n + 2)})
    fd = natural_frequencies(am.Ks, am.Ms, 5)
    fs = natural_frequencies(am.Ks, am.Ms, 5, sparse_above=10)
    assert np.allclose(fd, fs, rtol=1e-9)
    f, modes = natural_frequencies(am.Ks, am.Ms, 3, return_modes=True)
    M = am.Ms.toarray()
    assert np.allclose(modes.T @ M @ modes, np.eye(3), atol=1e-10)


def test_recovery_rows_follow_input_order_across_batches():
    node_xyz = {**_brick_nodes(0, 0.0), **_brick_nodes(8, 1.0)}
    n2 = (5, 6, 7, 8, 13, 14, 15, 16)
    recs = [ElemRecord(5, 1, 1, tuple(range(1, 9)), False, MAT, dict(eint=1)),
            ElemRecord(5, 2, 1, n2, False, MAT, dict(eint=0)),
            ElemRecord(5, 3, 1, tuple(range(1, 9)), False, MAT, dict(eint=1))]
    dm = build_dofmap(sorted(node_xyz), None, recs)
    rec = assemble(node_xyz, recs, dm).recovery["SOLID"]
    assert list(rec["id"]) == [1, 2, 3]
    S2 = solid.recovery(np.array([node_xyz[n] for n in n2]), MAT, eint=0)
    assert np.allclose(rec["S"][1], S2)


# ---- review fixes: excavated flag scope, PLANE orientation warning, SHELL repeated corners ----
def test_excavated_flag_counts_only_for_solid_and_plane():
    """Only SOLID/PLANE can be excavated soil (requirements 4.4 item 1).  ETYPE 2 on a SHELL is a
    buried shell, i.e. structure (spec 08 4.6), and ETYPE has no effect on BEAMS (spec 08 5.11):
    such records go to Ks/Ms with a message, never into the subtracted Ke/Me."""
    node_xyz = {**_brick_nodes(), 20: (0.0, 0.0, 3.0), 30: (5.0, 5.0, 5.0)}
    KR = np.zeros((12, 12))
    KR[0, 0] = KR[6, 6] = 10.0
    KR[0, 6] = KR[6, 0] = -10.0

    def model(flag):
        return [ElemRecord(1, 1, 1, tuple(range(1, 9)), True, MAT, {}),                    # soil
                ElemRecord(3, 1, 3, (5, 6, 7, 8), flag, MAT, dict(thick=0.1)),            # buried shell
                ElemRecord(2, 1, 2, (5, 20, 30), flag, MAT, dict(section=SEC)),
                ElemRecord(7, 1, 7, (6, 20), flag, None, dict(k=(1, 2, 3, 0, 0, 0))),
                ElemRecord(9, 1, 9, (7, 20), flag, None, dict(KR=KR))]

    dm = build_dofmap(sorted(node_xyz), None, model(False))
    ref = assemble(node_xyz, model(False), dm)
    am = assemble(node_xyz, model(True), dm)
    assert (am.Ks != ref.Ks).nnz == 0 and (am.Ms != ref.Ms).nnz == 0
    assert (am.Ke != ref.Ke).nnz == 0 and (am.Me != ref.Me).nnz == 0
    Ksoil, _ = solid.matrices(np.array([node_xyz[n] for n in range(1, 9)]), MAT)
    eq = dm.eqs(range(1, 9), (1, 2, 3))
    assert np.allclose(am.Ke.toarray()[np.ix_(eq, eq)], Ksoil) and am.Ke.nnz == ref.Ke.nnz
    for name in ("SHELL", "BEAMS", "SPRING", "GENERAL"):
        assert any(m.startswith(f"1 {name} elements flagged excavated assembled as structure") for m in am.messages)
    assert any("buried shell" in m for m in am.messages)
    assert not ref.messages
    assert list(am.recovery["SOLID"]["excavated"]) == [1]
    assert list(am.recovery["SHELL"]["excavated"]) == [0] and list(am.recovery["BEAMS"]["excavated"]) == [0]
    from sassi.elements.assemble import is_excavated_soil
    assert [is_excavated_soil(r) for r in model(True)] == [True, False, False, False, False]


def test_clockwise_plane_elements_are_reported():
    """Requirements 4.1 / spec 08 Q12: clockwise PLANE input is reordered internally *with a warning*."""
    node_xyz = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 0, 1), 4: (0, 0, 1), 5: (2, 0, 0), 6: (2, 0, 1)}
    ccw = [ElemRecord(1, 1, 4, (1, 2, 3, 4), False, MAT, {}), ElemRecord(1, 2, 4, (2, 5, 6, 3), False, MAT, {})]
    dm = build_dofmap(sorted(node_xyz), None, ccw)
    assert not any("clockwise" in m for m in assemble(node_xyz, ccw, dm).messages)
    mixed = [ccw[0], ElemRecord(1, 2, 4, (2, 3, 6, 5), False, MAT, {}),             # clockwise
             ElemRecord(1, 3, 4, (1, 4, 3), False, MAT, {})]                       # clockwise triangle
    am = assemble(node_xyz, mixed, dm)
    msg = [m for m in am.messages if "clockwise" in m]
    assert len(msg) == 1 and msg[0].startswith("2 PLANE elements")
    assert "group 1 element 2" in msg[0] and "group 1 element 3" in msg[0] and "element 1," not in msg[0]
    # the physics does not depend on the orientation: same Ks for the reordered element 2
    am_ccw = assemble(node_xyz, ccw, dm)
    am_cw = assemble(node_xyz, [ccw[0], mixed[1]], dm)
    assert np.abs(am_ccw.Ks - am_cw.Ks).max() <= 1e-12 * np.abs(am_ccw.Ks).max()


def test_shell_records_with_repeated_corner_nodes():
    """A 4-node SHELL record with one pair of adjacent repeated nodes is the triangle of its
    distinct nodes (generalising L = K); other repetitions and distinct nodes at the same point
    are CHECK Error 9 (spec 08 5.5), never a ZeroDivisionError."""
    from sassi.elements.assemble import element_dof_nodes
    node_xyz = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0.2, 1, 0), 4: (0, 0, 0)}
    mk = lambda nodes: ElemRecord(3, 1, 3, nodes, False, MAT, dict(thick=0.1))  # noqa: E731
    ref = [mk((1, 2, 3, 0))]
    dm = build_dofmap([1, 2, 3, 4], None, ref)
    am_ref = assemble(node_xyz, ref, dm)
    for nodes in ((1, 2, 3, 3), (1, 2, 3, 1), (1, 1, 2, 3), (1, 2, 2, 3)):
        assert element_dof_nodes(mk(nodes))[0] == [1, 2, 3]
        am = assemble(node_xyz, [mk(nodes)], dm)
        assert np.abs(am.Ks - am_ref.Ks).max() <= 1e-12 * np.abs(am_ref.Ks).max()
        assert (am.Ms != am_ref.Ms).nnz == 0
        assert np.allclose(am.recovery["SHELL"]["S"], am_ref.recovery["SHELL"]["S"])
    for nodes in ((1, 2, 1, 3), (1, 1, 2, 2), (1, 2, 2), (1, 1, 2, 0)):
        with pytest.raises(ElementError, match="Error 9"):
            build_dofmap([1, 2, 3, 4], None, [mk(nodes)])
    rec = mk((1, 2, 3, 4))                       # node 4 is distinct from node 1 but coincides with it
    dm4 = build_dofmap([1, 2, 3, 4], None, [rec])
    with pytest.raises(ElementError, match="Error 9"):
        assemble(node_xyz, [rec], dm4)
