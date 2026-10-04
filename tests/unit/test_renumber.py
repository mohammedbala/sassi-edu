"""Reverse Cuthill-McKee renumbering and the HOUSE "Optimize Model" option (sassi.core.renumber,
sassi.modules.house; requirements 4.4 item 5, D-HOU-03, D-FIL-09)."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.core import renumber as RN
from sassi.io import decks
from sassi.io.files import read_container
from sassi.modules.house import renumbered_deck
from sassi.verify.problems.vp_house import add_element, add_node, house_deck, mixed_house_deck, run_house

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")


def _grid(nx, ny, perm=None):
    perm = np.arange(nx * ny) if perm is None else perm
    idx = lambda i, j: int(perm[i + j * nx])
    cl = [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)] for j in range(ny - 1) for i in range(nx - 1)]
    z = np.zeros(nx * ny)
    for j in range(ny):
        for i in range(nx):
            z[idx(i, j)] = j
    return cl, z, idx


# ---------------------------------------------------------------------------------------- graph tools
def test_adjacency_symmetric_without_diagonal():
    A = RN.node_adjacency(4, [[0, 1, 2], [2, 3], [3]])
    D = A.toarray()
    assert (D == D.T).all() and not D.diagonal().any()
    assert D[0, 1] and D[0, 2] and D[2, 3] and not D[0, 3]


def test_pseudo_peripheral_pair_of_a_path():
    A = RN.node_adjacency(5, [[k, k + 1] for k in range(4)])
    assert sorted(RN.pseudo_peripheral_pair(A, 2)) == [0, 4]


def test_measures():
    A = RN.node_adjacency(4, [[0, 3], [1, 2]])
    assert RN.bandwidth(A, [0, 1, 2, 3]) == 3 and RN.profile(A, [0, 1, 2, 3]) == 3 + 1
    assert RN.bandwidth(A, [0, 3, 1, 2]) == 1 and RN.profile(A, [0, 3, 1, 2]) == 2


def test_keep_relative_order():
    out = RN.keep_relative_order([3, 0, 4, 1, 2], [1, 3])
    assert list(out) == [1, 0, 4, 3, 2]
    assert list(RN.keep_relative_order([2, 1, 0], [])) == [2, 1, 0]


def test_rcm_reduces_a_scrambled_grid_and_runs_bottom_up():
    nx, ny = 12, 9
    perm = np.random.default_rng(4).permutation(nx * ny)
    cl, z, idx = _grid(nx, ny, perm)
    r = RN.optimize_numbering(nx * ny, cl, z=z)
    assert r.method == "RCM" and r.reduced
    assert 3 * r.bandwidth_after < r.bandwidth_before and 3 * r.profile_after < r.profile_before
    assert sorted(r.new_number) == list(range(1, nx * ny + 1))
    # bottom-up: the mean new number grows with the row (z)
    means = [np.mean([r.new_number[idx(i, j)] for i in range(nx)]) for j in range(ny)]
    assert all(np.diff(means) > 0)


def test_original_order_kept_when_not_better():
    cl = [[k, k + 1] for k in range(9)]                       # a path numbered along itself
    r = RN.optimize_numbering(10, cl)
    assert r.method == "original order kept" and list(r.new_number) == list(range(1, 11))
    assert r.profile_after == r.profile_before and r.notes


def test_interaction_order_and_isolated_nodes():
    nx, ny = 8, 6
    perm = np.random.default_rng(9).permutation(nx * ny)
    cl, z, idx = _grid(nx, ny, perm)
    n = nx * ny + 3                                          # 3 isolated nodes (no element)
    z = np.concatenate([z, [0.0, 0.0, 0.0]])
    keep = sorted((idx(i, j) for j in (0, 1) for i in range(nx)), key=lambda k: (z[k], k))
    r = RN.optimize_numbering(n, cl, keep=keep, z=z)
    assert all(np.diff(r.new_number[keep]) > 0)              # relative order kept
    assert sorted(r.new_number[-3:]) == [n - 2, n - 1, n]    # isolated nodes last, in order
    assert list(r.new_number[-3:]) == [n - 2, n - 1, n] and r.n_isolated == 3
    # an interaction node without elements is not 'isolated': it keeps its place in the sequence
    r2 = RN.optimize_numbering(n, cl, keep=keep + [n - 1], z=z)
    assert r2.n_isolated == 2 and r2.new_number[n - 1] > max(r2.new_number[keep])


def test_map_file_round_trip(tmp_path):
    p = RN.write_map(tmp_path / "m.map", [30, 10, 20], [2, 3, 1])
    assert p.read_text() == "10 3\n20 1\n30 2\n"
    assert RN.read_map(p) == {10: 3, 20: 1, 30: 2}
    (tmp_path / "bad.map").write_text("1 2\n1 3\n")
    with pytest.raises(ValueError, match="mapped twice"):
        RN.read_map(tmp_path / "bad.map")
    (tmp_path / "bad2.map").write_text("1 2 3\n")
    with pytest.raises(ValueError, match="expected 'old new'"):
        RN.read_map(tmp_path / "bad2.map")


# ---------------------------------------------------------------------------------------- HOUSE
def test_renumbered_deck_maps_every_node_column():
    d = house_deck()
    for k, z in ((1, 0.0), (2, 0.0), (3, 1.0), (4, 1.0), (5, 2.0)):
        add_node(d, k, 0.0, float(k), z)
    d.table("interaction").append([2])
    d.table("interaction").append([1])
    d.table("groups").append([1, 2, "beams"])
    add_element(d, 1, 1, [1, 3, 5])
    d.table("masses").append([4, 1, 1, 1, 0, 0, 0, 0])
    d.table("symm").append([1, 1, 3, 4, 0])
    d["optimize"] = 1
    mp = {1: 5, 2: 4, 3: 3, 4: 2, 5: 1}
    out = renumbered_deck(d, mp)
    assert [r[0] for r in out.table("nodes").rows] == [1, 2, 3, 4, 5]
    assert [r[3] for r in out.table("nodes").rows] == [2.0, 1.0, 1.0, 0.0, 0.0]     # sorted by the new id
    assert [r[0] for r in out.table("interaction").rows] == [4, 5]                    # order kept
    assert out.rows("elements")[0]["n1"] == 5 and out.rows("elements")[0]["n3"] == 1
    assert out.rows("elements")[0]["n4"] == 0 and out.rows("masses")[0]["node"] == 2
    assert out.rows("symm")[0]["n2"] == 2 and out.rows("symm")[0]["n3"] == 0
    assert out["optimize"] == 0 and d["optimize"] == 1


def _coosk(wd):
    return read_container(wd / "COOSK", "COOSK").sparse("Ks"), read_container(wd / "m.N4", "FILE4")


def test_house_optimizer_matrices_are_a_permutation(tmp_path):
    """mixed_house_deck (SOLID, BEAMS, SHELL, SPRING, GENERAL with K nodes, masses, excavation and
    interaction nodes): the optimised Ks is P Ks P^T of the original one through the .map."""
    a, b = tmp_path / "a", tmp_path / "b"
    d = mixed_house_deck()
    rc, out = run_house(a, d)
    assert rc == 0, out[-1500:]
    d["optimize"] = 1
    rc, out = run_house(b, d)
    assert rc == 0, out[-1500:]
    assert "node-numbering optimizer" in out and "use the NEW numbers" in out
    assert (b / "m.hounew").exists() and (b / "m.map").exists() and not (a / "m.map").exists()
    mp = RN.read_map(b / "m.map")
    Ka, fa = _coosk(a)
    Kb, fb = _coosk(b)
    assert sorted(mp.values()) == list(range(1, len(mp) + 1)) == sorted(int(v) for v in fb["node_id"])
    eqb = {(int(n), int(k)): i for i, (n, k) in enumerate(zip(fb["eq_node"], fb["eq_dof"]))}
    perm = [eqb[(mp[int(n)], int(k))] for n, k in zip(fa["eq_node"], fa["eq_dof"])]
    assert abs(Kb[perm][:, perm] - Ka).max() <= 1e-14 * abs(Ka).max()
    assert [int(v) for v in fb["int_node"]] == [mp[int(v)] for v in fa["int_node"]]
    assert all(np.diff([int(v) for v in fb["int_node"]]) > 0)
    assert fb.meta["x_optimized"] == 1 and fb.meta["x_node_map"] == "m.map"
    assert [mp[int(o)] for o in fb["x_node_old_id"]] == [int(v) for v in fb["node_id"]]
    hn = decks.read(b / "m.hounew", "HOUSE")
    assert int(hn["optimize"]) == 0 and len(hn.table("nodes")) == len(mp)


def test_house_optimizer_data_check_and_errors(tmp_path):
    d = mixed_house_deck()
    d["optimize"], d["opmode"] = 1, 1
    rc, out = run_house(tmp_path / "c", d)
    assert rc == 0 and "Data check only: .hounew and .map are not written" in out
    assert not (tmp_path / "c" / "m.map").exists()
    d = mixed_house_deck()
    d["optimize"] = 1
    add_element(d, 1, 99, [1, 2, 3, 4, 5, 6, 7, 999])          # undefined node: input error
    rc, out = run_house(tmp_path / "e", d)
    assert rc == 1 and "optimizer (HOUSEX <optimize> = 1) skipped" in out
    assert not (tmp_path / "e" / "m.map").exists()


def test_stale_optimizer_files_are_retired(tmp_path):
    d = mixed_house_deck()
    d["optimize"] = 1
    rc, _ = run_house(tmp_path, d)
    assert rc == 0 and (tmp_path / "m.map").exists()
    d["optimize"] = 0
    rc, out = run_house(tmp_path, d)
    assert rc == 0 and not (tmp_path / "m.map").exists() and (tmp_path / "m.map.prev").exists()
    assert (tmp_path / "m.hounew.prev").exists() and "m.map.prev (renamed" in out
    f4 = read_container(tmp_path / "m.N4", "FILE4")
    assert "x_optimized" not in f4.meta and "x_node_old_id" not in f4.arrays


# ---------------------------------------------------------------------------------------- interaction-rooted RCM
def _hex_grid(nx, ny, nz, seed=None):
    """nx x ny x nz node grid of hexahedra, bottom layer = interaction nodes (in x, y order), optionally
    numbered at random (``seed``)."""
    idx = lambda i, j, k: i + nx * (j + ny * k)
    cl = [[idx(i, j, k), idx(i + 1, j, k), idx(i + 1, j + 1, k), idx(i, j + 1, k), idx(i, j, k + 1), idx(i + 1, j, k + 1),
           idx(i + 1, j + 1, k + 1), idx(i, j + 1, k + 1)] for k in range(nz - 1) for j in range(ny - 1) for i in range(nx - 1)]
    z = np.array([k for k in range(nz) for j in range(ny) for i in range(nx)], dtype=float)
    keep = [idx(i, j, 0) for j in range(ny) for i in range(nx)]
    if seed is not None:
        perm = np.random.default_rng(seed).permutation(nx * ny * nz)
        cl = [[int(perm[v]) for v in c] for c in cl]
        zz = np.empty_like(z)
        zz[perm] = z
        z, keep = zz, [int(perm[v]) for v in keep]
    return nx * ny * nz, cl, keep, z


def test_rcm_from_the_interaction_nodes_numbers_them_last_in_order():
    """Roots of the CM search: the interaction nodes get the last numbers of their component, ascending in
    their given order, and every level structure follows (no re-sorting needed)."""
    n, cl, keep, z = _hex_grid(6, 5, 4, seed=3)
    adj = RN.node_adjacency(n, cl)
    order = RN.reverse_cuthill_mckee(adj, z=z, roots=keep)
    assert sorted(order) == list(range(n)) and list(order[-len(keep):]) == keep


def test_scrambled_wide_shallow_mesh_gets_a_narrow_band():
    """Wide, shallow hex mesh numbered at random, bottom layer = interaction nodes: the band shrinks by an order
    of magnitude and the interaction nodes become ascending.  The classical RCM with the interaction nodes
    re-sorted into their slots is much wider, and with the dense free-field impedance (all interaction nodes
    coupled, as in the ANALYS system) also has the larger profile."""
    n, cl, keep, z = _hex_grid(20, 20, 8, seed=1)
    r = RN.optimize_numbering(n, cl, keep=keep, z=z)
    assert r.method == "RCM" and all(np.diff(r.new_number[keep]) > 0)
    assert 7 * r.bandwidth_after < r.bandwidth_before and 3 * r.profile_after < r.profile_before
    adj, sysadj = RN.node_adjacency(n, cl), RN.node_adjacency(n, cl + [keep])
    resorted = RN.keep_relative_order(RN.reverse_cuthill_mckee(adj, z=z), keep)
    assert RN.bandwidth(adj, resorted) > 5 * r.bandwidth_after
    assert RN.profile(sysadj, r.order) == r.profile_after < RN.profile(sysadj, resorted)


def test_a_candidate_that_widens_the_band_is_rejected():
    """Spec 05b test 8: hex grid numbered bottom-up layer by layer (the manual's recommendation, already
    optimal): no candidate beats it in both measures, the original order is kept and the note says why."""
    n, cl, keep, z = _hex_grid(12, 12, 4)
    r = RN.optimize_numbering(n, cl, keep=keep, z=z)
    assert r.method == "original order kept" and list(r.new_number) == list(range(1, n + 1))
    assert r.bandwidth_after == r.bandwidth_before and r.profile_after == r.profile_before
    assert "spec 05b test 8" in r.notes[0]


def test_equation_measures_of_an_assembled_pattern():
    """Two nodes with 3 and 2 DOFs and a coupling: equation numbering node by node, DOFs ascending."""
    import scipy.sparse as sp
    eq_node = np.array([0, 0, 0, 1, 1])
    eq_dof = np.array([1, 2, 3, 1, 3])
    A = sp.csr_matrix(np.array([[1, 1, 0, 0, 1], [1, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 1]]))
    assert RN.equation_measures([0, 1], eq_node, eq_dof, A) == (4, 1 + 4)
    # node 1 first: its equations become 0, 1 and node 0's become 2, 3, 4: coupling (0 dof1, 1 dof3) -> (2, 1)
    assert RN.equation_measures([1, 0], eq_node, eq_dof, A) == (1, 1 + 1)


def test_house_lists_the_measures_of_its_assembled_matrices(tmp_path):
    """The optimizer section quotes the bandwidth and profile of the COOSK |Ks| + |Ke| of the run without
    and with the optimizer (VP-O1 model, which the optimizer renumbers)."""
    import re
    import scipy.sparse as sp
    from sassi.verify import builders as B
    from sassi.verify.problems.vp_tshell import equation_profile, optimizer_model
    site, mdl, d0 = optimizer_model()
    rc_a, out_a = run_house(tmp_path / "a", d0)
    d0["optimize"] = 1
    rc_b, out_b = run_house(tmp_path / "b", d0)
    assert rc_a == 0 and rc_b == 0, out_b[-1500:]
    pat = lambda w: (abs(read_container(tmp_path / w / "COOSK", "COOSK").sparse("Ks"))
                     + abs(read_container(tmp_path / w / "COOSK", "COOSK").sparse("Ke"))).tocsr()
    sym = lambda A: (A + A.T).tocsr()          # |Ks| + |Ke| can lose a round-off cancellation on one side only
    (bwa, pfa), (bwb, pfb) = equation_profile(sym(pat("a"))), equation_profile(sym(pat("b")))
    bw = re.search(r"Bandwidth before/after\s*: (\d+) / (\d+)  \(equation matrix", out_b).groups()
    pf = re.search(r"Profile before/after\s*: (\d+) / (\d+)", out_b).groups()
    assert (int(bw[0]), int(bw[1]), int(pf[0]), int(pf[1])) == (bwa, bwb, pfa, pfb)
    assert bwb < bwa and pfb < pfa


def test_failed_optimised_run_quotes_the_model_numbers(tmp_path):
    """An error found by the HOUSE checks quotes the node number of the model (.hou), and neither .map nor
    .hounew is written."""
    d = mixed_house_deck()
    d["optimize"] = 1
    first = int(d.table("interaction").rows[0][0])
    for r in d.table("nodes").rows:
        if int(r[0]) == first:
            r[3] = float(r[3]) + 0.37
    rc, out = run_house(tmp_path, d)
    assert rc == 1 and f"EDU-01: interaction node {first} " in out
    assert not (tmp_path / "m.map").exists() and not (tmp_path / "m.hounew").exists()
