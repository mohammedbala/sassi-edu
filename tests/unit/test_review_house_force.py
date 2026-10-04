"""Adversarial review tests of the HOUSE and FORCE modules (work package house_force).

These tests are written independently of the implementer's builders: the decks are filled
with ``sassi.io.decks.new`` and run through ``sassi.modules.base.run_module``; the references
are physical invariants (rigid-body null space, translation invariance, the zero-SSI identity
K_s = K_e for a structure made of the excavated soil), hand-built matrices or time-domain
identities.  Tests marked ``xfail(strict=True)`` document defects found in the review (the
reason names the defect); they flip to XPASS -> FAIL once the defect is fixed, so the marker
must then be removed.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from sassi.conventions import cfactor
from sassi.elements import material_from_layer, solid
from sassi.io import decks
from sassi.io.files import read_container
from sassi.modules.base import run_module

G = 9.81
#: (L no, thick, weight, vp, vs, dp, ds); the last row is the half-space
LAYERS = [(1, 2.0, 18.0, 400.0, 200.0, 0.05, 0.05),
          (2, 3.0, 19.0, 600.0, 300.0, 0.04, 0.04),
          (3, 0.0, 20.0, 1000.0, 500.0, 0.02, 0.02)]


# ======================================================================================
# small deck helpers (independent of sassi/verify/problems/vp_house.py)
# ======================================================================================
def _site(wd: Path, gravity=G, nl=20, cmodform=0):
    s = decks.new("SITE")
    s["model"] = "r"
    s["gravity"] = gravity
    s["nl"] = nl
    s["cmodform"] = cmodform
    s["hs"] = LAYERS[-1][0]
    for row in LAYERS[:-1]:
        s.table("layers").append(list(row))
    s.table("halfspace").append(list(LAYERS[-1]))
    decks.write(wd / "r.sit", s)


def _house(gelev=0.0, gravity=G, **params) -> decks.Deck:
    d = decks.new("HOUSE")
    d["model"] = "r"
    d["gravity"] = gravity
    d["gelev"] = gelev
    for k, v in params.items():
        d[k] = v
    for row in LAYERS:
        d.table("layers").append(list(row))
    for row in LAYERS[:-1]:
        d.table("sitelayers").append(list(row))
    d.table("sitelayers").append([LAYERS[-1][0], 0.0] + list(LAYERS[-1][2:]))
    return d


def _node(d, n, x, y, z, fix=(0,) * 6):
    d.table("nodes").append([n, float(x), float(y), float(z)] + list(fix))


def _elem(d, g, e, nodes, etype=1, mat=1, prop=1, eint=0, thick=0.0, ki="000000", kj="000000"):
    d.table("elements").append([g, e, etype, mat, prop, eint, thick] + (list(nodes) + [0] * 8)[:8] + [ki, kj])


def _run(wd: Path, d: decks.Deck, **site):
    wd.mkdir(parents=True, exist_ok=True)
    _site(wd, gravity=float(d["gravity"]), cmodform=int(d["cmodform"]), **site)
    decks.write(wd / "r.hou", d)
    rc = run_module("HOUSE", "r", wd)
    return rc, (wd / "r_HOUSE.out").read_text()


def _block(d, x0=0.0, y0=0.0, zs=(-5.0, -2.0, 0.0), nx=2, ny=2, h=2.0, etype=2, mats=(2, 1), group=1,
           start=1, interaction=True, skew=0.0):
    """nx x ny x len(zs)-1 hexahedra; layer k (bottom first) uses mats[k]; returns {(i,j,k): id}."""
    ids = {}
    n = start
    for k, z in enumerate(zs):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ids[(i, j, k)] = n
                x = x0 + i * h + (skew * z if 0 < i < nx else 0.0)
                _node(d, n, x, y0 + j * h, z)
                if interaction:
                    d.table("interaction").append([n])
                n += 1
    d.table("groups").append([group, 1, "block"])
    e = 0
    for k in range(len(zs) - 1):
        for j in range(ny):
            for i in range(nx):
                e += 1
                c = [ids[(i, j, k)], ids[(i + 1, j, k)], ids[(i + 1, j + 1, k)], ids[(i, j + 1, k)]]
                t = [ids[(i, j, k + 1)], ids[(i + 1, j, k + 1)], ids[(i + 1, j + 1, k + 1)], ids[(i, j + 1, k + 1)]]
                _elem(d, group, e, c + t, etype=etype, mat=mats[k])
    return ids


def _files(wd: Path):
    return (read_container(wd / "r.N4", "FILE4"), read_container(wd / "COOSK", "COOSK"),
            read_container(wd / "COOSM", "COOSM"))


def _eq(f4, node, dof):
    k = np.flatnonzero((f4["eq_node"] == node) & (f4["eq_dof"] == dof))
    return int(k[0]) if k.size else -1


def _rigid_modes(f4) -> np.ndarray:
    """6 rigid-body modes (3 translations, 3 small rotations about the origin) on the FILE4 numbering."""
    xyz = {int(n): p for n, p in zip(f4["node_id"], f4["node_xyz"])}
    R = np.zeros((f4["eq_node"].size, 6))
    for e, (n, c) in enumerate(zip(f4["eq_node"], f4["eq_dof"])):
        x, y, z = xyz[int(n)]
        R[e, c - 1] = 1.0
        if c == 1:
            R[e, 4], R[e, 5] = z, -y
        elif c == 2:
            R[e, 3], R[e, 5] = -z, x
        elif c == 3:
            R[e, 3], R[e, 4] = y, -x
    return R


# ======================================================================================
# HOUSE: physical invariants of the assembled matrices
# ======================================================================================
def test_excavated_block_rigid_body_null_space_and_mass(tmp_path):
    """K*_e of a free excavated block annihilates the 6 rigid-body modes (any damping, any
    distortion, 3x3x3 Gauss), is complex symmetric (not Hermitian), and r^T M_e r = sum rho V."""
    d = _house()
    _block(d, skew=0.15)
    for r in d.tables["elements"].rows:
        r[5] = 1                                                     # EINT 1: 3x3x3 Gauss
    rc, out = _run(tmp_path, d)
    assert rc == 0, out
    f4, ck, cm = _files(tmp_path)
    Ke, Me = ck.sparse("Ke"), cm.sparse("Me")
    R = _rigid_modes(f4)
    assert np.abs(Ke @ R).max() <= 1e-12 * abs(Ke).max() * 10.0
    assert abs(Ke - Ke.T).max() <= 1e-14 * abs(Ke).max()
    assert abs(Ke.imag).max() > 0 and abs(Ke - Ke.conj().T).max() > 1e-3 * abs(Ke).max()
    assert abs(Me - Me.T).max() <= 1e-14 * abs(Me).max() and np.all(Me.diagonal() > 0)
    # skew is a shear of the interior columns: volumes unchanged = 16 x (2 + 3)
    m_ref = 16.0 * (2.0 * 18.0 + 3.0 * 19.0) / G
    for c in range(3):
        np.testing.assert_allclose(R[:, c] @ (Me @ R[:, c]), m_ref, rtol=1e-12)
    m91 = read_container(tmp_path / "FILE91").meta
    np.testing.assert_allclose(m91["mass_excavated"], [m_ref] * 3, rtol=1e-12)


def test_structure_of_soil_equals_excavated_soil(tmp_path):
    """Zero-SSI identity at matrix level (D-ELM-12, D-CNV-04, VP-16 precondition): a structure
    (ETYPE 1) with M-table type-3 materials equal to the L layers (Vp, Vs, weight, damping) has
    K*_s = K*_e and M_s = M_e exactly when incompatible modes are suppressed; with MOPT <incomp> = 0
    only the structure gets them (D-ELM-02)."""
    de = _house(incomp=0)
    _block(de, etype=2)
    rc, out = _run(tmp_path / "exc", de)
    assert rc == 0, out
    Ke = read_container(tmp_path / "exc" / "COOSK").sparse("Ke")
    Me = read_container(tmp_path / "exc" / "COOSM").sparse("Me")
    for incomp, equal in ((1, True), (0, False)):
        ds = _house(incomp=incomp)
        _block(ds, etype=1)
        for row in LAYERS:
            ds.table("materials").append([row[0], 3, row[3], row[4], row[2], row[5], row[6]])
        wd = tmp_path / f"str{incomp}"
        rc, out = _run(wd, ds)
        assert rc == 0, out
        Ks = read_container(wd / "COOSK").sparse("Ks")
        Ms = read_container(wd / "COOSM").sparse("Ms")
        assert abs(Ms - Me).max() == 0.0
        if equal:
            assert abs(Ks - Ke).max() == 0.0
        else:
            assert abs(Ks - Ke).max() > 1e-6 * abs(Ke).max()


def test_ground_elevation_shift_invariance(tmp_path):
    """Moving the whole model and the ground elevation by the same amount changes nothing but the
    coordinates: same interfaces, same K*_e, ETYPE 0 still resolved to excavated (D-HOU-01 with
    gelev != 0), same layer hash, different interaction-node hash (coordinates)."""
    d0 = _house()
    _block(d0, etype=2)
    rc, out = _run(tmp_path / "a", d0)
    assert rc == 0, out
    d1 = _house(gelev=12.5)
    _block(d1, zs=(7.5, 10.5, 12.5), etype=0)
    rc, out = _run(tmp_path / "b", d1)
    assert rc == 0, out
    a, b = _files(tmp_path / "a"), _files(tmp_path / "b")
    np.testing.assert_array_equal(a[0]["int_iface"], b[0]["int_iface"])
    np.testing.assert_array_equal(b[0]["elem_excav"], np.ones(8, int))
    Ka, Kb = a[1].sparse("Ke"), b[1].sparse("Ke")
    assert abs(Ka - Kb).max() <= 1e-13 * abs(Ka).max()
    ha, hb = read_container(tmp_path / "a" / "FILE90").meta, read_container(tmp_path / "b" / "FILE90").meta
    assert ha["layer_hash"] == hb["layer_hash"] and ha["int_hash"] != hb["int_hash"]


def test_etype0_boundary_cases(tmp_path):
    """D-HOU-01: an element whose top face is at grade is excavated; an element with one node
    above grade by much more than the tolerance is structure (M table)."""
    d = _house()
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("groups").append([1, 1, "etype 0"])
    xy = ((0, 0), (1, 0), (1, 1), (0, 1))
    n = 0
    for e, top in enumerate((0.0, 0.0)):
        ns = []
        for z in (-2.0, top):
            for x, y in xy:
                n += 1
                zz = z + (0.5 if (e == 1 and z == top and (x, y) == (1, 1)) else 0.0)
                _node(d, n, x + 3 * e, y, zz)
                ns.append(n)
        _elem(d, 1, e + 1, ns, etype=0, mat=1)
    for k in range(1, 9):
        d.table("interaction").append([k])
    rc, out = _run(tmp_path, d)
    assert rc == 0, out
    f4 = read_container(tmp_path / "r.N4")
    np.testing.assert_array_equal(f4["elem_excav"], [1, 0])


# ======================================================================================
# HOUSE: interaction-node interfaces (D-GEN-06, EDU-01)
# ======================================================================================
def _mat_on_soil(d, zb, n0=1):
    """Structural slab SOLID 2x2x1 above grade; its 4 bottom nodes (elevations zb) are interaction nodes."""
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("groups").append([1, 1, "slab"])
    xy = ((0, 0), (2, 0), (2, 2), (0, 2))
    for k, ((x, y), z) in enumerate(zip(xy, zb)):
        _node(d, n0 + k, x, y, z)
        d.table("interaction").append([n0 + k])
    for k, (x, y) in enumerate(xy):
        _node(d, n0 + 4 + k, x, y, 1.0)
    _elem(d, 1, 1, list(range(n0, n0 + 8)), etype=1, mat=1)


@pytest.mark.parametrize("factor,ok", [(0.9, True), (1.1, False)])
def test_interface_tolerance_boundary(tmp_path, factor, ok):
    """itol = max(1e-6 L_ref, 1e-9, 1e-4 h_min) = 2e-4 here (h_min = 2): a node 0.9 itol below the
    surface is on interface 1, 1.1 itol is EDU-01."""
    d = _house()
    _mat_on_soil(d, (0.0, -factor * 2e-4, 0.0, 0.0))
    rc, out = _run(tmp_path, d)
    if ok:
        assert rc == 0, out
        np.testing.assert_array_equal(read_container(tmp_path / "r.N4")["int_iface"], [1, 1, 1, 1])
    else:
        assert rc == 1 and "EDU-01: interaction node 2" in out


def test_interaction_nodes_on_deep_interfaces_and_in_halfspace(tmp_path):
    """A mat on the half-space top (depth 5 = interface 3) is accepted; a mat inside the half-space
    (depth 6, no user interface) is EDU-01."""
    d = _house(gelev=0.0)
    _mat_on_soil(d, (-5.0,) * 4)
    for r in d.tables["nodes"].rows[4:]:
        r[3] = -4.0
    rc, out = _run(tmp_path / "a", d)
    assert rc == 0, out
    np.testing.assert_array_equal(read_container(tmp_path / "a" / "r.N4")["int_iface"], [3, 3, 3, 3])
    d = _house(gelev=0.0)
    _mat_on_soil(d, (-6.0,) * 4)
    for r in d.tables["nodes"].rows[4:]:
        r[3] = -5.0
    rc, out = _run(tmp_path / "b", d)
    assert rc == 1 and "EDU-01: interaction node 1 at depth 6" in out


def test_interaction_order_is_the_deck_order(tmp_path):
    """int_node keeps the interaction order of the deck (not sorted); a top-down order is EDU-21
    (warning when coherent, D-HOU-05)."""
    d = _house()
    _block(d, etype=2, interaction=False)
    order = [n for n in range(27, 0, -1)]                         # top-down
    for n in order:
        d.table("interaction").append([n])
    rc, out = _run(tmp_path, d)
    assert rc == 0, out
    f4 = read_container(tmp_path / "r.N4")
    np.testing.assert_array_equal(f4["int_node"], order)
    np.testing.assert_array_equal(f4["int_iface"], [1] * 9 + [2] * 9 + [3] * 9)
    assert "EDU-21" in out


# ======================================================================================
# HOUSE: element data through the deck
# ======================================================================================
def test_cmodform1_reaches_layers_and_springs(tmp_path):
    """CMODFORM 1 (1 + 2 i beta): the excavated element is solid.matrices with material_from_layer(form=1)
    and a SPRING is k (1 + 2 i beta) (D-CNV-03)."""
    d = _house(cmodform=1)
    ids = _block(d, nx=1, ny=1, zs=(-2.0, 0.0), etype=2, mats=(1,))
    d.table("springprops").append([1, 100.0, 200.0, 300.0, 0.0, 0.0, 0.0, 0.1])
    d.table("groups").append([2, 7, "spring"])
    _node(d, 50, 5.0, 0.0, 0.0, fix=(1,) * 6)
    _elem(d, 2, 1, [ids[(1, 0, 1)], 50], prop=1)
    rc, out = _run(tmp_path, d)
    assert rc == 0, out
    f4, ck, _ = _files(tmp_path)
    lay = material_from_layer(400.0, 200.0, 18.0, 0.05, 0.05, G, form=1)
    xyz = np.array([f4["node_xyz"][list(f4["node_id"]).index(ids[k])] for k in
                    ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))])
    Kel, _ = solid.matrices(xyz, lay, incompatible=False, eint=0)
    eqs = [_eq(f4, ids[k], c) for k in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1),
                                         (1, 1, 1), (0, 1, 1)) for c in (1, 2, 3)]
    Ke = ck.sparse("Ke").toarray()
    np.testing.assert_allclose(Ke[np.ix_(eqs, eqs)], Kel, rtol=0, atol=1e-12 * np.abs(Kel).max())
    Ks = ck.sparse("Ks").toarray()
    n = ids[(1, 0, 1)]
    for c, k in ((1, 100.0), (2, 200.0), (3, 300.0)):
        e = _eq(f4, n, c)
        assert Ks[e, e] == pytest.approx(k * (1 + 0.2j), rel=1e-14)
        assert cfactor(0.1, 1) == pytest.approx(1 + 0.2j)


def test_beam_end_release_through_the_deck(tmp_path):
    """KJ '000011' (M2, M3 released at J): the recovered end moments MYJ and MZJ vanish for every
    nodal displacement, the other end forces do not (spec 08 5.16-5.17, D-ELM-04)."""
    d = _house()
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("beamprops").append([1, 0.25, 0.2, 0.2, 0.01, 0.004, 0.006])
    d.table("groups").append([1, 2, "beam"])
    _node(d, 1, 0.0, 0.0, 0.0)
    _node(d, 2, 3.0, 1.0, 2.0)
    _node(d, 3, 0.0, 5.0, 0.0, fix=(1,) * 6)
    _elem(d, 1, 1, [1, 2, 3], kj="000011")
    rc, out = _run(tmp_path, d)
    assert rc == 0, out
    f4 = read_container(tmp_path / "r.N4")
    S = f4["rec_BEAMS_S"][0]                                         # (12 comps, 12 dofs)
    scale = np.abs(S).max()
    comps = f4.meta["components"]["BEAMS"]
    assert np.abs(S[comps.index("MYJ")]).max() <= 1e-12 * scale
    assert np.abs(S[comps.index("MZJ")]).max() <= 1e-12 * scale
    for c in ("FXI", "MXJ", "MYI", "MZI"):
        assert np.abs(S[comps.index(c)]).max() > 1e-6 * scale


def test_general_upper_rows_and_weight_units(tmp_path):
    """GENERAL (2 nodes, global): MXR/MXI/MXM upper-triangle rows give K = K_R + i K_I with the lower
    triangle filled by symmetry, and the mass in weight units is divided by g (MOPT <matrix> = 1)."""
    rng = np.random.default_rng(7)
    A = rng.standard_normal((12, 12))
    KR = A @ A.T * 1e3
    KI = 0.05 * KR
    MM = np.diag(rng.uniform(1.0, 3.0, 12)) * G
    d = _house(gmunits=1)
    _node(d, 1, 0.0, 0.0, 0.0)
    _node(d, 2, 1.0, 2.0, 3.0)
    d.table("groups").append([1, 9, "gm"])
    for kind, Mx in (("R", KR), ("I", KI), ("M", MM)):
        for r in range(12):
            d.table("matrices").append([1, kind, r + 1] + list(np.triu(Mx)[r, r:]) + [0.0] * r)
    _elem(d, 1, 1, [1, 2], prop=1)
    rc, out = _run(tmp_path, d)
    assert rc == 0, out
    f4, ck, cm = _files(tmp_path)
    eqs = [_eq(f4, n, c) for n in (1, 2) for c in range(1, 7)]
    K = ck.sparse("Ks").toarray()[np.ix_(eqs, eqs)]
    M = cm.sparse("Ms").toarray()[np.ix_(eqs, eqs)]
    np.testing.assert_allclose(K, KR + 1j * KI, rtol=0, atol=1e-12 * np.abs(KR).max())
    np.testing.assert_allclose(M, MM / G, rtol=1e-14)


def test_rotary_nodal_mass_weight_units(tmp_path):
    """MR in weight units (MUNITS 1) is divided by g on ROTX..ROTZ; MT in mass units is not."""
    d = _house()
    d.table("materials").append([1, 1, 3.0e7, 0.2, 0.0, 0.0, 0.0])          # massless beam
    d.table("beamprops").append([1, 0.25, 0.0, 0.0, 0.01, 0.004, 0.006])
    d.table("groups").append([1, 2, "beam"])
    _node(d, 1, 0.0, 0.0, 0.0, fix=(1,) * 6)
    _node(d, 2, 0.0, 0.0, 4.0)
    _node(d, 3, 1.0, 0.0, 0.0, fix=(1,) * 6)
    _elem(d, 1, 1, [1, 2, 3])
    d.table("masses").append([2, 0.0, 0.0, 0.0, 3.0 * G, 4.0 * G, 5.0 * G, 1])
    d.table("masses").append([2, 7.0, 7.0, 7.0, 0.0, 0.0, 0.0, 0])
    rc, out = _run(tmp_path, d)
    assert rc == 0, out
    f4, _, cm = _files(tmp_path)
    Ms = cm.sparse("Ms")
    for dof, ref in ((1, 7.0), (2, 7.0), (3, 7.0), (4, 3.0), (5, 4.0), (6, 5.0)):
        e = _eq(f4, 2, dof)
        assert Ms[e, e] == pytest.approx(ref, rel=1e-14)


def test_mass_on_fixed_dof_in_total_only(tmp_path):
    """A nodal mass on a fixed DOF is ignored in M_s (W5) but counted in the FILE91 total mass."""
    d = _house()
    d.table("materials").append([1, 1, 3.0e7, 0.2, 0.0, 0.0, 0.0])
    d.table("beamprops").append([1, 0.25, 0.0, 0.0, 0.01, 0.004, 0.006])
    d.table("groups").append([1, 2, "beam"])
    _node(d, 1, 0.0, 0.0, 0.0, fix=(1, 1, 1, 1, 1, 1))
    _node(d, 2, 0.0, 0.0, 4.0, fix=(0, 0, 1, 0, 0, 0))
    _node(d, 3, 1.0, 0.0, 0.0, fix=(1,) * 6)
    _elem(d, 1, 1, [1, 2, 3])
    d.table("masses").append([2, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0])
    rc, out = _run(tmp_path, d)
    assert rc == 0, out
    f4, _, cm = _files(tmp_path)
    assert _eq(f4, 2, 3) == -1 and "W5/W6" in out
    m91 = read_container(tmp_path / "FILE91").meta
    np.testing.assert_allclose(m91["mass_structure"], [2.0, 2.0, 2.0], rtol=1e-14)
    assert abs(cm.sparse("Ms").sum() - 4.0) < 1e-12                        # UX + UY only


# ======================================================================================
# HOUSE: outputs, hashes, determinism (D-GEN-09, D-W1-13)
# ======================================================================================
def _array_map(c):
    return {k: c.arrays[k] for k in c.arrays}


def test_outputs_are_bit_identical_between_runs(tmp_path):
    for k in ("a", "b"):
        d = _house(incomp=0)
        _block(d, skew=0.1)
        rc, out = _run(tmp_path / k, d)
        assert rc == 0, out
    for name in ("r.N4", "DOFSMAP"):
        A, B = read_container(tmp_path / "a" / name), read_container(tmp_path / "b" / name)
        assert set(A.arrays) == set(B.arrays)
        for k in A.arrays:
            if sp.issparse(A.arrays[k]):
                continue
            assert np.array_equal(np.asarray(A[k]), np.asarray(B[k])), (name, k)
    for name, keys in (("COOSK", ("Ks", "Ke")), ("COOSM", ("Ms", "Me"))):
        A, B = read_container(tmp_path / "a" / name), read_container(tmp_path / "b" / name)
        for k in keys:
            assert (A.sparse(k) != B.sparse(k)).nnz == 0
    assert read_container(tmp_path / "a" / "FILE90").meta == read_container(tmp_path / "b" / "FILE90").meta


def test_restart_hash_selectivity(tmp_path):
    """requirements 2.5: a structural change (nodal mass) changes only file4_hash; gravity changes
    the layer hash; moving a surface interaction node in plan changes the interaction hash."""
    def model(mass=1.0, gravity=G, dx=0.0):
        d = _house(gravity=gravity)
        _mat_on_soil(d, (0.0,) * 4)
        d.tables["nodes"].rows[1][1] += dx
        d.table("masses").append([5, mass, mass, mass, 0.0, 0.0, 0.0, 0])
        return d
    hs = {}
    for k, kw in (("base", {}), ("mass", dict(mass=2.0)), ("grav", dict(gravity=9.80665)), ("move", dict(dx=0.1))):
        rc, out = _run(tmp_path / k, model(**kw))
        assert rc == 0, out
        hs[k] = read_container(tmp_path / k / "FILE90").meta
    b = hs["base"]
    assert hs["mass"]["int_hash"] == b["int_hash"] and hs["mass"]["layer_hash"] == b["layer_hash"]
    assert hs["mass"]["file4_hash"] != b["file4_hash"]
    assert hs["grav"]["layer_hash"] != b["layer_hash"] and hs["grav"]["int_hash"] == b["int_hash"]
    assert hs["move"]["int_hash"] != b["int_hash"] and hs["move"]["layer_hash"] == b["layer_hash"]


def test_2d_interaction_equations_are_in_plane_only(tmp_path):
    d = _house(dim=1)
    d.table("groups").append([1, 4, "soil"])
    ids = {}
    n = 0
    for iz, z in enumerate((-2.0, 0.0)):
        for ix, x in enumerate((0.0, 1.0, 2.0)):
            n += 1
            ids[(ix, iz)] = n
            _node(d, n, x, 0.0, z)
            d.table("interaction").append([n])
    for ix in range(2):
        _elem(d, 1, ix + 1, [ids[(ix, 0)], ids[(ix + 1, 0)], ids[(ix + 1, 1)], ids[(ix, 1)]], etype=2, mat=1)
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("beamprops").append([1, 0.25, 0.0, 0.0, 0.01, 0.004, 0.006])
    d.table("groups").append([2, 2, "column"])
    _node(d, 50, 1.0, 0.0, 3.0)
    _node(d, 99, 5.0, 0.0, 0.0, fix=(1,) * 6)
    _elem(d, 2, 1, [ids[(1, 1)], 50, 99])
    rc, out = _run(tmp_path, d)
    assert rc == 0, out
    f4 = read_container(tmp_path / "r.N4")
    assert np.all(f4["int_eq"][:, 1] == -1)


def test_edu10_uses_delt_nfft_when_df_is_zero(tmp_path):
    d = _house(df=0.0, delt=0.005, nft=4096)
    _block(d, etype=2)
    d.table("freqs").append([800])                                    # 39 Hz > Vs/(5h) = 20 Hz
    rc, out = _run(tmp_path, d)
    assert rc == 0, out
    assert "EDU-10" in out


# ======================================================================================
# FORCE
# ======================================================================================
def _force(wd: Path, loads, fnum, mforce=1, **params):
    wd.mkdir(parents=True, exist_ok=True)
    d = decks.new("FORCE")
    d["model"] = "r"
    d["mforce"] = mforce
    for k, v in params.items():
        d[k] = v
    for row in loads:
        d.table("loads").append(list(row))
    for n in fnum:
        d.table("freqs").append([n])
    decks.write(wd / "r.frc", d)
    rc = run_module("FORCE", "r", wd)
    return rc, (wd / "r_FORCE.out").read_text()


def test_file9_reproduces_delayed_scaled_histories(tmp_path):
    """D-FRC-01 end to end: with every frequency 1..N/2-1 in FILE9, irfft(P F_ref) for each load is the
    reference history scaled by the factor and delayed by the arrival time (circular shift), for
    loads in any order, including moments and a negative factor."""
    n, dt = 512, 0.01
    t = np.arange(n) * dt
    f_ref = np.exp(-((t - 0.6) / 0.08) ** 2) * np.cos(2 * np.pi * 5 * t)
    f_ref -= f_ref.mean()                                              # no f = 0 component
    loads = [(30, 5, -1.5, 0.37), (4, 1, 2.0, 0.0), (4, 3, 0.5, 1.23), (17, 2, 1.0, 2.0)]
    rc, out = _force(tmp_path, loads, fnum=list(range(n // 2 - 1, 0, -1)), delt=dt, nft=n)
    assert rc == 0, out
    f9 = read_container(tmp_path / "FILE9", "FILE9")
    np.testing.assert_allclose(f9["freq"], np.arange(1, n // 2) / (n * dt), rtol=1e-14)
    F = np.fft.rfft(f_ref)
    for node, dof, a, t0 in loads:
        k = int(np.flatnonzero((f9["load_node"] == node) & (f9["load_dof"] == dof))[0])
        spec = np.zeros(n // 2 + 1, complex)
        spec[1:n // 2] = f9["P"][:, k] * F[1:n // 2]
        y = np.fft.irfft(spec, n)
        expect = np.fft.irfft(np.r_[0.0, F[1:n // 2], 0.0], n)
        np.testing.assert_allclose(y, a * np.roll(expect, int(round(t0 / dt))), atol=1e-12)
    assert "MY" in out and "FY" in out


def test_moment_dof_labels_and_amplitude(tmp_path):
    rc, out = _force(tmp_path, [(9, 4, 2.0, 0.0), (9, 6, -3.0, 0.5)], fnum=(1, 2), delt=0.01, nft=256)
    assert rc == 0, out
    f9 = read_container(tmp_path / "FILE9", "FILE9")
    np.testing.assert_array_equal(f9["load_dof"], [4, 6])
    np.testing.assert_allclose(np.abs(f9["P"]), [[2.0, 3.0], [2.0, 3.0]], rtol=1e-15)
    lines = [ln for ln in out.splitlines() if ln.strip().startswith("9 ")]
    assert any(" MX " in ln for ln in lines) and any(" MZ " in ln for ln in lines)


@pytest.mark.parametrize("params,fnum,msg", [
    (dict(delt=0.01, nft=256), (1, 2, 2), "duplicate"),
    (dict(delt=0.01, nft=256), (0, 1), "positive"),
    (dict(delt=0.001, nft=65536), (1,), "exceeds the FORCE maximum"),
    (dict(delt=0.0, nft=0, fstep=-1.0), (1,), "Error 48"),
])
def test_force_input_errors(tmp_path, params, fnum, msg):
    rc, out = _force(tmp_path, [(1, 1, 1.0, 0.0)], fnum=fnum, **params)
    assert rc == 1 and msg in out and not (tmp_path / "FILE9").exists()


def test_force_nfft_not_power_of_two_warns(tmp_path):
    rc, out = _force(tmp_path, [(1, 1, 1.0, 0.0)], fnum=(1,), delt=0.01, nft=3000)
    assert rc == 0 and "Warning 9" in out
    assert read_container(tmp_path / "FILE9").meta["df"] == pytest.approx(1.0 / 30.0, rel=1e-15)


def test_add_mode_zero_factor_does_not_change_the_load(tmp_path):
    rc, out = _force(tmp_path, [(3, 1, 1.0, 0.1), (3, 1, 0.0, 0.5)], fnum=(1, 2, 3), mforce=0, delt=0.01, nft=256)
    assert rc == 0, out
    f9 = read_container(tmp_path / "FILE9")
    np.testing.assert_allclose(f9["P"][:, 0], np.exp(-2j * np.pi * f9["freq"] * 0.1), rtol=1e-14)


def test_add_mode_arrival_matches_the_interpreter(tmp_path):
    from sassi.prep import Interpreter
    ui = Interpreter(cwd=tmp_path)
    ui.run_text("MDL,r,{d}\nN,3,0,0,0\nMOPT,1,0,1,0\nF,3,1.0,0,0,0.3,0,0\nF,3,4.0,0,0,0.2,0,0\n".format(d=tmp_path))
    t_ui = ui.model.forces[3].arrival[0]
    rc, out = _force(tmp_path / "f", [(3, 1, 1.0, 0.3), (3, 1, 4.0, 0.2)], fnum=(1,), mforce=0, delt=0.01, nft=256)
    assert rc == 0, out
    assert read_container(tmp_path / "f" / "FILE9")["x_arrival"][0] == pytest.approx(t_ui)
