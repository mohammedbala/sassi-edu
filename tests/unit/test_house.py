"""Unit tests of the HOUSE module (requirements 4.1, 4.4; D-HOU-01..05, D-W1-13).

The decks are built with the public builders of sassi/verify/problems/vp_house.py exactly as AFWRITE
writes them (``decks.new('HOUSE')``) and run through ``run_module``.
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from sassi.conventions import ELEMENT_COMPONENTS
from sassi.core import house_lib as hl
from sassi.elements import material_from_layer, material_from_M
from sassi.io import decks
from sassi.io.files import SCHEMAS, read_container, validate
from sassi.modules import house
from sassi.verify.problems.vp_house import (G_SI, SITE_LAYERS, add_element, add_matrix_property, add_node,
                                            house_deck, mixed_house_deck, mixed_reference, run_house,
                                            soil_block_deck, write_site_deck)


def _dense(A):
    return A.toarray() if sp.issparse(A) else np.asarray(A)


def _close(A, B, rtol=1e-12):
    A, B = _dense(A), _dense(B)
    return np.abs(A - B).max() <= rtol * max(np.abs(B).max(), 1e-300)


def surface_block(z_bottom=(0.0, 0.0, 0.0, 0.0), etype=1, gelev=0.0, z_top=1.0, **params):
    """One structural SOLID with its four bottom nodes as interaction nodes (surface foundation)."""
    d = house_deck(gelev=gelev, **params)
    xy = [(0, 0), (2, 0), (2, 2), (0, 2)]
    for k, ((x, y), zb) in enumerate(zip(xy, z_bottom)):
        add_node(d, k + 1, x, y, zb)
        d.table("interaction").append([k + 1])
    for k, (x, y) in enumerate(xy):
        add_node(d, k + 5, x, y, z_top)
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("groups").append([1, 1, "mat"])
    add_element(d, 1, 1, list(range(1, 9)), etype=etype, mat=1)
    return d


# ======================================================================================
# deck -> FILE4 / COOSK / COOSM of the mixed model
# ======================================================================================
@pytest.fixture(scope="module")
def mixed(tmp_path_factory):
    wd = tmp_path_factory.mktemp("mixed")
    d = mixed_house_deck()
    rc, out = run_house(wd, d)
    assert rc == 0, out
    return wd, d, out


def test_mixed_files_and_schemas(mixed):
    wd, d, out = mixed
    for name, kind in (("m.N4", "FILE4"), ("COOSK", "COOSK"), ("COOSM", "COOSM"), ("DOFSMAP", "DOFSMAP"),
                       ("FILE90", "FILE90"), ("FILE91", "FILE91")):
        c = read_container(wd / name, kind)
        assert validate(c) == [], (name, validate(c))
        arrays, meta = SCHEMAS[kind]
        assert set(arrays) <= set(c.arrays) and set(meta) <= set(c.meta)
    f4 = read_container(wd / "m.N4")
    assert f4.meta["nEq"] == f4["eq_node"].size and f4.meta["gravity"] == G_SI
    assert f4.meta["model_hash"] == read_container(wd / "FILE90").meta["file4_hash"]
    ds = read_container(wd / "DOFSMAP")
    np.testing.assert_array_equal(ds["eq_node"], f4["eq_node"])
    np.testing.assert_array_equal(ds["int_node"], f4["int_node"])
    assert "Excavated soil elements: layer assignment" in out and "Interaction nodes per user interface" in out


def test_mixed_matrices_equal_direct_assembly(mixed):
    wd, d, _ = mixed
    dm, am, recs = mixed_reference(d)
    f4 = read_container(wd / "m.N4")
    ck, cm = read_container(wd / "COOSK"), read_container(wd / "COOSM")
    np.testing.assert_array_equal(f4["eq_node"], dm.eq_node)
    np.testing.assert_array_equal(f4["eq_dof"], dm.eq_dof)
    assert _close(ck.sparse("Ks"), am.Ks) and _close(ck.sparse("Ke"), am.Ke)
    assert _close(cm.sparse("Ms"), am.Ms) and _close(cm.sparse("Me"), am.Me)
    assert np.iscomplexobj(ck.sparse("Ks").data) and not np.iscomplexobj(cm.sparse("Ms").data)
    # Ks and Ke are complex symmetric (never Hermitian with damping)
    Ks = ck.sparse("Ks")
    assert abs(Ks - Ks.T).max() <= 1e-12 * abs(Ks).max()
    assert abs(Ks.imag).max() > 0


def test_mixed_recovery_operators(mixed):
    wd, d, _ = mixed
    dm, am, _ = mixed_reference(d)
    f4 = read_container(wd / "m.N4")
    comps = f4.meta["components"]
    assert set(comps) == {"SOLID", "SHELL", "BEAMS", "SPRING"}            # GENERAL has no components
    for t, rd in am.recovery.items():
        np.testing.assert_allclose(f4[f"rec_{t}_S"], rd["S"], rtol=0, atol=1e-12 * np.abs(rd["S"]).max())
        np.testing.assert_array_equal(f4[f"rec_{t}_eq"], rd["eq"])
        idx = f4[f"rec_{t}_idx"]
        np.testing.assert_array_equal(f4["elem_group"][idx], rd["group"])
        np.testing.assert_array_equal(f4["elem_id"][idx], rd["id"])
        assert comps[t] == ELEMENT_COMPONENTS[t]
        assert f4[f"rec_{t}_S"].shape[1] == len(ELEMENT_COMPONENTS[t])
    # elem arrays: excavated flag and node slots
    exc = f4["elem_excav"].astype(bool)
    assert exc.sum() == 4 and np.all(f4["elem_type"][exc] == 1) and np.all(f4["elem_group"][exc] == 6)
    assert f4["elem_nodes"].shape == (f4["elem_id"].size, 8)
    # SOLID strain operator for STRESS / nonlinear soil
    assert f4["x_rec_SOLID_B"].shape == (5, 6, 24)


def test_mixed_spring_recovery_is_force(mixed):
    """SPRING recovery F = k* (u_J - u_I) applied to a unit relative displacement (D-STR-05)."""
    wd, _, _ = mixed
    f4 = read_container(wd / "m.N4")
    S, eq = f4["rec_SPRING_S"][0], f4["rec_SPRING_eq"][0]
    u = np.zeros(f4.meta["nEq"], complex)
    u[f4["int_eq"][list(f4["int_node"]).index(15), 0]] = -1.0          # node I = 15 moves -1 in X
    F = S @ np.where(eq >= 0, u[np.clip(eq, 0, None)], 0)
    from sassi.conventions import cfactor
    np.testing.assert_allclose(F[0], 1.0e5 * cfactor(0.03), rtol=1e-14)
    assert np.all(F[1:] == 0)


def test_mixed_interaction_nodes(mixed):
    wd, _, _ = mixed
    f4 = read_container(wd / "m.N4")
    np.testing.assert_array_equal(f4["int_node"], np.arange(1, 19))
    np.testing.assert_array_equal(f4["int_iface"], np.repeat([3, 2, 1], 6))   # z = -5, -2, 0
    assert np.all(f4["int_eq"] >= 0)
    for k, n in enumerate(f4["int_node"]):
        for c in range(3):
            e = f4["int_eq"][k, c]
            assert f4["eq_node"][e] == n and f4["eq_dof"][e] == c + 1
    np.testing.assert_allclose(f4["x_iface_depth"], [0.0, 2.0, 5.0])


def test_mixed_mass_totals_and_ignored_mass(mixed):
    wd, _, out = mixed
    m91 = read_container(wd / "FILE91").meta
    soil = (2 * 1 * 2) * 18.0 / G_SI + (2 * 1 * 3) * 19.0 / G_SI
    np.testing.assert_allclose(m91["mass_excavated"], [soil] * 3, rtol=1e-13)
    assert m91["n_interaction"] == 18 and m91["n_excavated"] == 4
    assert "node 26 dof 1" in out and "W5/W6" in out                   # mass on the fixed K-only node
    assert "EDU-06" in out                                              # shell drilling at nodes 20, 21


def test_incompatible_modes_only_structural(tmp_path):
    """MOPT <incomp> = 0 changes the structural SOLID, never the excavated soil (D-ELM-02)."""
    rc0, _ = run_house(tmp_path / "a", mixed_house_deck(incomp=0))
    rc1, _ = run_house(tmp_path / "b", mixed_house_deck(incomp=1))
    assert rc0 == rc1 == 0
    ka, kb = read_container(tmp_path / "a" / "COOSK"), read_container(tmp_path / "b" / "COOSK")
    assert _close(ka.sparse("Ke"), kb.sparse("Ke"), rtol=0)
    assert not _close(ka.sparse("Ks"), kb.sparse("Ks"), rtol=1e-6)


# ======================================================================================
# interaction interfaces: tolerance and errors (EDU-01, G-12)
# ======================================================================================
def test_interface_tolerance_accepts_small_misfit(tmp_path):
    # itol = max(tol, 1e-4 h_min) = 2e-4 (h_min = 2): misfits of 1e-4 are on the surface
    rc, out = run_house(tmp_path, surface_block(z_bottom=(0.0, -1e-4, 1e-4, 0.0)))
    assert rc == 0, out
    f4 = read_container(tmp_path / "m.N4")
    np.testing.assert_array_equal(f4["int_iface"], [1, 1, 1, 1])


def test_interface_error_off_interface(tmp_path):
    rc, out = run_house(tmp_path, surface_block(z_bottom=(0.0, -1e-3, 0.0, 0.0), etype=1))
    assert rc == 1
    assert "EDU-01: interaction node 2 at depth 0.001 is not on a soil-layer interface" in out
    assert not (tmp_path / "m.N4").exists()


def test_interface_error_above_grade(tmp_path):
    rc, out = run_house(tmp_path, surface_block(gelev=-0.5))
    assert rc == 1 and "above the ground elevation" in out


def test_interfaces_with_ground_elevation(tmp_path):
    """depth = gelev - z: a mat at z = 2 with gelev = 4 sits on interface 2 (depth 2)."""
    d = surface_block(z_bottom=(2.0,) * 4, z_top=3.0, gelev=4.0)
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out
    np.testing.assert_array_equal(read_container(tmp_path / "m.N4")["int_iface"], [2, 2, 2, 2])


def test_match_interfaces_helper():
    dep = hl.interface_depths([2.0, 3.0])
    itol = hl.interface_tolerance(dep, 1e-6)
    assert itol == pytest.approx(2e-4)
    k, mis, above = hl.match_interfaces(np.array([0.0, -2.0001, -5.0, -3.0, 0.001]), 0.0, dep, itol)
    np.testing.assert_array_equal(k[:3], [1, 2, 3])
    assert np.all(mis[:3] <= itol) and mis[3] == pytest.approx(1.0) and above[4] and not above[0]
    assert hl.layer_index_at_depth(dep, 1.0) == 1 and hl.layer_index_at_depth(dep, 4.0) == 2
    assert hl.layer_index_at_depth(dep, 7.0) == 3


def test_rigid_base_interface_error(tmp_path):
    d, _ = soil_block_deck()
    write_site_deck(tmp_path, nl=0)
    rc, out = run_house(tmp_path, d, site=False)
    assert rc == 1 and "rigid base" in out


# ======================================================================================
# classification (D-HOU-01), materials, units
# ======================================================================================
def test_etype0_resolution(tmp_path):
    d = house_deck()
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("groups").append([1, 1, "etype 0"])
    zlev = {"below": (-2.0, 0.0), "above": (0.0, 1.0), "straddle": (-1.0, 1.0)}
    nid = 0
    for e, (key, (z0, z1)) in enumerate(zlev.items()):
        ns = []
        for z in (z0, z1):
            for x, y in ((0, 0), (1, 0), (1, 1), (0, 1)):
                nid += 1
                add_node(d, nid, x + 3 * e, y, z)
                ns.append(nid)
        add_element(d, 1, e + 1, ns, etype=0, mat=1)
    for n in (1, 2, 3, 4, 5, 6, 7, 8):              # excavated block: FV interaction set
        d.table("interaction").append([n])
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out
    f4 = read_container(tmp_path / "m.N4")
    np.testing.assert_array_equal(f4["elem_excav"], [1, 0, 0])
    np.testing.assert_array_equal(f4["x_elem_etype"], [2, 1, 1])
    assert "ETYPE 0 resolved for 3 SOLID/PLANE elements: 1 excavated, 2 structure (D-HOU-01)" in out
    # the excavated element uses L layer 1 (MSET 1), the others M material 1
    lay = material_from_layer(400.0, 200.0, 18.0, 0.05, 0.05, G_SI)
    ck = read_container(tmp_path / "COOSK")
    from sassi.elements import solid
    xyz = np.array([[x, y, z] for z in (-2.0, 0.0) for x, y in ((0, 0), (1, 0), (1, 1), (0, 1))], float)
    Ke_el, _ = solid.matrices(xyz, lay)
    eq = [int(np.flatnonzero((f4["eq_node"] == n) & (f4["eq_dof"] == c))[0]) for n in range(1, 9) for c in (1, 2, 3)]
    assert _close(ck.sparse("Ke").toarray()[np.ix_(eq, eq)], Ke_el)


def test_resolve_etype_helper():
    tol = 1e-9
    assert hl.resolve_etype(1, 0, [-1, 0], 0.0, tol) == 2
    assert hl.resolve_etype(1, 0, [0, 1], 0.0, tol) == 1
    assert hl.resolve_etype(1, 0, [0, 0], 0.0, tol) == 1          # centroid not strictly below
    assert hl.resolve_etype(4, 2, [5, 6], 0.0, tol) == 2           # explicit ETYPE wins
    assert hl.resolve_etype(3, 2, [0], 0.0, tol) == 2              # buried shell
    assert hl.resolve_etype(2, 2, [0], 0.0, tol) == 1              # BEAMS: no effect


def _two_node_model(gmunits=0, mm=1.0, munits=1, mt=G_SI):
    """GENERAL element (diagonal spring pattern) between nodes 1 (fixed) and 2, plus an MT mass at 2."""
    d = house_deck(gmunits=gmunits)
    add_node(d, 1, 0, 0, 0, fix=(1,) * 6)
    add_node(d, 2, 1, 0, 0)
    d.table("groups").append([1, 9, "gm"])
    K = np.zeros((12, 12))
    for i in range(6):
        K[i, i] = K[i + 6, i + 6] = 1000.0 * (i + 1)
        K[i, i + 6] = K[i + 6, i] = -1000.0 * (i + 1)
    add_matrix_property(d, 1, "R", K)
    add_matrix_property(d, 1, "M", np.eye(12) * mm)
    add_element(d, 1, 1, [1, 2], prop=1)
    d.table("masses").append([2, mt, mt, mt, 0.0, 0.0, 0.0, munits])
    return d


def test_gmunits_and_munits(tmp_path):
    """UT-15: GENERAL mass in weight units (MOPT <matrix> = 1) and MT in weight units (MUNITS = 1)
    equal the same masses in mass units."""
    runs = {"weight": _two_node_model(gmunits=1, mm=G_SI, munits=1, mt=G_SI),
            "mass": _two_node_model(gmunits=0, mm=1.0, munits=0, mt=1.0)}
    M = {}
    for k, d in runs.items():
        rc, out = run_house(tmp_path / k, d)
        assert rc == 0, out
        M[k] = read_container(tmp_path / k / "COOSM").sparse("Ms").toarray()
    np.testing.assert_allclose(M["weight"], M["mass"], rtol=1e-14)
    np.testing.assert_allclose(np.diag(M["mass"]), [2, 2, 2, 1, 1, 1])     # GM 1 + MT 1 on translations


def test_material_types_and_layer_material():
    """Structural elements use material_from_M with the deck gravity; excavated ones the L layer."""
    d = mixed_house_deck()
    hm = house.build_house_model(d)
    assert not hm.errors
    m3 = material_from_M(2, 4.0e7, 1.5e7, 23.0, 0.03, 0.02, G_SI)
    assert hm.materials[3].G == m3.G and hm.materials[3].M == m3.M and hm.materials[3].rho == m3.rho
    l1 = material_from_layer(400.0, 200.0, 18.0, 0.05, 0.05, G_SI)
    assert hm.layers[1].G == l1.G and hm.layers[1].M == l1.M
    rec = [r for r in hm.records if r.code == 2][1]
    assert rec.props["ki"] == "000011" and rec.props["section"]["I3"] == 0.006
    gm = [r for r in hm.records if r.code == 9][0]
    assert gm.props["munits"] == 1 and gm.props["upper_rows"]


def test_input_errors_are_listed(tmp_path):
    d = mixed_house_deck()
    d.tables["materials"].rows = [r for r in d.tables["materials"].rows if r[0] != 2]   # SHELL material
    d.table("groups").append([7, 5, "thick shell"])
    add_element(d, 7, 1, [20, 23, 24, 21], mat=1, thick=0.3)
    rc, out = run_house(tmp_path, d)
    assert rc == 1
    assert "Error 13: material 2 of element 1 of group 2 is not defined" in out
    assert "TSHELL elements are tier P2" not in out          # TSHELL is supported since wave 3
    assert "error(s) in the HOUSE input" in out


# ======================================================================================
# options, prerequisites, data check
# ======================================================================================
def test_requires_sit(tmp_path):
    rc, out = run_house(tmp_path, mixed_house_deck(), site=False)
    assert rc == 1 and "HOUSE needs m.sit (run AFWRITE with SITE enabled)" in out


@pytest.mark.parametrize("param,msg", [("me", "Error 115"), ("ansys", "Option AA")])
def test_unavailable_options_stop(tmp_path, param, msg):
    # incoherency, wave passage and multiple excitation are implemented since wave 3 (package incoherency):
    # <me> = 1 without ME zones stops with Error 115 (requirements 4.4 item 6); Option AA is still P2
    d = mixed_house_deck()
    d[param] = 1
    rc, out = run_house(tmp_path, d)
    assert rc == 1 and msg in out and not (tmp_path / "m.N4").exists()


@pytest.mark.parametrize("param,msg", [("optimize", "node-numbering optimizer"), ("nlssi", "non-linear soil"),
                                       ("wpass", "wave passage"),
                                       ("coh", "Incoherency, wave passage and multiple excitation")])
def test_p1_options_warn_and_continue(tmp_path, param, msg):
    d = mixed_house_deck()
    d[param] = 1
    rc, out = run_house(tmp_path, d)
    assert rc == 0 and msg in out and (tmp_path / "m.N4").exists()


def test_data_check_mode_writes_no_files(tmp_path):
    d = mixed_house_deck()
    d["opmode"] = 1
    rc, out = run_house(tmp_path, d)
    assert rc == 0 and "Data check only" in out
    assert not (tmp_path / "m.N4").exists() and not (tmp_path / "COOSK").exists()


def test_site_table_mismatch_warning(tmp_path):
    layers = [list(r) for r in SITE_LAYERS]
    layers[0][4] = 210.0
    write_site_deck(tmp_path, layers=[tuple(r) for r in layers])
    rc, out = run_house(tmp_path, mixed_house_deck(), site=False)
    assert rc == 0 and "differs between the HOUSE and .sit decks (vs)" in out


# ======================================================================================
# interaction-set rules (requirements 4.4 item 3)
# ======================================================================================
def test_excstrchk_and_edu22(tmp_path):
    d, nid = soil_block_deck()
    centre = nid[(1, 1, 1)]                            # interior node of the 2x2x2 block
    top = nid[(1, 1, 2)]
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("beamprops").append([1, 0.25, 0.0, 0.0, 0.01, 0.004, 0.006])
    d.table("groups").append([3, 2, "pile"])
    add_node(d, 100, 10.0, 0.0, 0.0, fix=(1,) * 6)
    add_element(d, 3, 1, [centre, top, 100], mat=1)
    rc, out = run_house(tmp_path / "fv", d)
    # final audit (tests/unit/test_audit_excstrchk.py): a beam on an interior excavation node is an error
    # also when the node interacts (FV), as in CHECK; only SOLID/PLANE elements on the excavation mesh warn
    assert rc == 1, out
    assert f"EXCSTRCHK: excavation interior node {centre} is shared with the structure" in out
    assert "ties that structure to the soil continuum" in out
    # the same interior node without interaction flag: EXCSTRCHK error and EDU-22 warning
    d.tables["interaction"].rows = [r for r in d.tables["interaction"].rows if r[0] != centre]
    rc, out = run_house(tmp_path / "fi", d)
    assert rc == 1
    assert f"EXCSTRCHK: excavation interior node {centre} is shared with the structure" in out
    assert "EDU-22: 1 excavation nodes are not interaction nodes (FV needs every excavation node)" in out


def test_g14_bottom_up_and_fixed_interaction(tmp_path):
    """A pile (BEAMS) with interaction nodes and no excavation: G-14 warning; nodes numbered top-down:
    EDU-21 warning (coherent); a fixed translation: FIXEDINT; all fixed: Error 124."""
    d = house_deck()
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("beamprops").append([1, 0.25, 0.0, 0.0, 0.01, 0.004, 0.006])
    d.table("groups").append([1, 2, "pile"])
    add_node(d, 1, 0, 0, 0.0)
    add_node(d, 2, 0, 0, -2.0, fix=(0, 0, 1, 0, 0, 0))
    add_node(d, 3, 0, 0, -5.0)
    add_node(d, 9, 1, 0, 0.0, fix=(1,) * 6)
    add_element(d, 1, 1, [1, 2, 9])
    add_element(d, 1, 2, [2, 3, 9])
    for n in (1, 2, 3):
        d.table("interaction").append([n])
    rc, out = run_house(tmp_path / "a", d)
    assert rc == 0, out
    assert "G-14" in out and "EDU-21" in out and "FIXEDINT: interaction node 2 has fixed translations ['UZ']" in out
    f4 = read_container(tmp_path / "a" / "m.N4")
    np.testing.assert_array_equal(f4["int_iface"], [1, 2, 3])
    assert f4["int_eq"][1, 2] == -1 and f4["int_eq"][1, 0] >= 0
    d.tables["nodes"].rows[1][4:7] = [1, 1, 1]
    rc, out = run_house(tmp_path / "b", d)
    assert rc == 1 and "Error 124: interaction node 2 has all translations fixed" in out


def test_edu08_layer_assignment_warning(tmp_path):
    d, _ = soil_block_deck()
    for r in d.tables["elements"].rows:
        r[3] = 1                                       # every excavated element MSET 1
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out
    assert "EDU-08: excavated element 1 of group 2 uses layer 1, the TOPL layer at its mid-depth 3.5 is 2" in out
    assert out.count("EDU-08") >= 4


def test_plane_2d_model(tmp_path):
    """2D (dim = 1): PLANE excavated soil in the X-Z plane; interaction DOFs UX, UZ only."""
    d = house_deck(dim=1)
    d.table("groups").append([1, 4, "soil"])
    ids = {}
    n = 0
    for iz, z in enumerate((-5.0, -2.0, 0.0)):
        for ix, x in enumerate((0.0, 1.0, 2.0)):
            n += 1
            ids[(ix, iz)] = n
            add_node(d, n, x, 0.0, z)
            d.table("interaction").append([n])
    e = 0
    for iz in range(2):
        for ix in range(2):
            e += 1
            add_element(d, 1, e, [ids[(ix, iz)], ids[(ix + 1, iz)], ids[(ix + 1, iz + 1)], ids[(ix, iz + 1)]],
                        etype=2, mat=2 if iz == 0 else 1)
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out
    f4 = read_container(tmp_path / "m.N4")
    assert set(f4["eq_dof"]) == {1, 3} and np.all(f4["int_eq"][:, 1] == -1)
    np.testing.assert_array_equal(f4["int_iface"], np.repeat([3, 2, 1], 3))
    Me = read_container(tmp_path / "COOSM").sparse("Me")
    m = hl.translational_mass(Me, f4["eq_dof"])
    np.testing.assert_allclose(m[[0, 2]], 2 * 2 * 18.0 / G_SI + 2 * 3 * 19.0 / G_SI, rtol=1e-13)
    assert m[1] == 0.0


# ======================================================================================
# FILE90 hashes (D-W1-13)
# ======================================================================================
def _hashes(wd):
    return read_container(wd / "FILE90").meta


def test_hashes_deterministic_and_selective(tmp_path):
    base = mixed_house_deck()
    rc, _ = run_house(tmp_path / "a", base)
    rc2, _ = run_house(tmp_path / "b", mixed_house_deck())
    assert rc == rc2 == 0
    ha, hb = _hashes(tmp_path / "a"), _hashes(tmp_path / "b")
    assert all(ha[k] == hb[k] for k in ("int_hash", "layer_hash", "file4_hash"))
    # new structure (material change): only file4_hash changes (requirements 2.5 "New Structure")
    d = mixed_house_deck()
    d.tables["materials"].rows[0][2] = 3.1e7
    run_house(tmp_path / "c", d)
    hc = _hashes(tmp_path / "c")
    assert hc["int_hash"] == ha["int_hash"] and hc["layer_hash"] == ha["layer_hash"]
    assert hc["file4_hash"] != ha["file4_hash"]
    # damping form changes the layer hash
    d = mixed_house_deck()
    d["cmodform"] = 1
    run_house(tmp_path / "d", d)
    assert _hashes(tmp_path / "d")["layer_hash"] != ha["layer_hash"]
    # a different interaction set changes the interaction hash
    d = mixed_house_deck()
    d.tables["interaction"].rows = d.tables["interaction"].rows[6:]
    rc, out = run_house(tmp_path / "e", d)
    assert rc == 0, out
    assert _hashes(tmp_path / "e")["int_hash"] != ha["int_hash"]


def test_stable_hash_sensitivity():
    a = np.arange(5)
    h = hl.stable_hash([("a", a)])
    assert h == hl.stable_hash([("a", a.copy())])
    assert h != hl.stable_hash([("a", a.astype(float))]) and h != hl.stable_hash([("b", a)])


def test_excavation_sets_helper():
    z = {}
    els = []
    n = 0
    nid = {}
    for iz, zz in enumerate((-2.0, -1.0, 0.0)):
        for iy in range(3):
            for ix in range(3):
                n += 1
                nid[(ix, iy, iz)] = n
                z[n] = zz
    for iz in range(2):
        for iy in range(2):
            for ix in range(2):
                els.append((1, [nid[(ix, iy, iz)], nid[(ix + 1, iy, iz)], nid[(ix + 1, iy + 1, iz)],
                                nid[(ix, iy + 1, iz)], nid[(ix, iy, iz + 1)], nid[(ix + 1, iy, iz + 1)],
                                nid[(ix + 1, iy + 1, iz + 1)], nid[(ix, iy + 1, iz + 1)]]))
    ex = hl.excavation_sets(els, z, 0.0, 1e-9)
    assert ex["interior"] == {nid[(1, 1, 1)]} and len(ex["nodes"]) == 27
    assert ex["top"] == {nid[(ix, iy, 2)] for ix in range(3) for iy in range(3)}
    assert nid[(1, 1, 2)] not in ex["fsin"] and nid[(0, 0, 2)] in ex["fsin"]


# ======================================================================================
# end to end: interpreter -> AFWRITE -> RUNHOUSE / RUNFORCE
# ======================================================================================
E2E_MODEL = """
MDL,demo,{dir}
TIT,HOUSE end-to-end
N,1,0,0,-5
N,2,10,0,-5
N,3,10,10,-5
N,4,0,10,-5
N,5,0,0,0
N,6,10,0,0
N,7,10,10,0
N,8,0,10,0
N,10,5,5,0
N,11,5,5,10
N,12,15,5,0
INT,1,8,1,1
M,1,4.32e5,0.25,0.15,0.05,0.05
R,1,1,0.8,0.8,0.2,0.1,0.1
L,1,5,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
TOPL,1,1
GROUP,1,SOLID
E,1,1,2,3,4,5,6,7,8
GROUP,2,BEAMS
E,1,10,11,12
GROUP,3,SPRING
SC,1,1e4,1e4,1e4,0,0,0,0.02
RACT,1
E,1,7,10
FREQ,1,2,4,8
SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,4096,1
POINT,0,1,4.5
HOUSE,32.2,0,0,2,0,0,0,0,0
FORCE,0
F,11,1,2,0,0.1,0,0
MT,11,32.2,32.2,32.2
AOPT,0,0,0,1,1,1,0,1,0,0,0,0,0,0
AFWRITE
RUNHOUSE
RUNFORCE
"""


def test_end_to_end_afwrite_runhouse_runforce(tmp_path):
    from sassi.prep import Interpreter, Kind
    mdir = tmp_path / "demo"
    mdir.mkdir()
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(E2E_MODEL.format(dir=mdir))
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    out = (mdir / "demo_HOUSE.out").read_text()
    assert "HOUSE finished with status OK" in out, out[-2000:]
    f4 = read_container(mdir / "demo.N4", "FILE4")
    assert f4.meta["afwrite_model_hash"] == ui.model.model_hash()
    np.testing.assert_array_equal(f4["int_node"], np.arange(1, 9))
    np.testing.assert_array_equal(f4["int_iface"], [2, 2, 2, 2, 1, 1, 1, 1])      # bottom-up numbering
    np.testing.assert_array_equal(f4["elem_excav"], [1, 0, 0])                    # ETYPE 0 below grade
    assert "EDU-21" not in out
    # beam: I = 10, J = 11 with all 6 DOFs; K node 12 has none (fixed by AFWRITE, D-CHK-08)
    assert set(f4["eq_dof"][f4["eq_node"] == 11]) == {1, 2, 3, 4, 5, 6}
    assert not np.any(f4["eq_node"] == 12)
    Ms = read_container(mdir / "COOSM").sparse("Ms")
    e = int(np.flatnonzero((f4["eq_node"] == 11) & (f4["eq_dof"] == 3))[0])   # axial DOF of the vertical beam
    beam_end = 0.15 / 32.2 * 1.0 * 10.0 / 3.0                                     # consistent axial mass rho A L / 3
    assert Ms[e, e] == pytest.approx(1.0 + beam_end, rel=1e-12)                    # MT 32.2 (weight) -> 1
    f9 = read_container(mdir / "FILE9", "FILE9")
    np.testing.assert_array_equal(f9["load_node"], [11, 11])
    np.testing.assert_array_equal(f9["load_dof"], [1, 2])
    np.testing.assert_allclose(f9["P"][:, 0], np.exp(-2j * np.pi * f9["freq"] * 0.1), rtol=1e-14)
    np.testing.assert_allclose(f9["P"][:, 1], 2.0)


def _edit_gravity(d):
    d["gravity"] = 0.0


def _edit_dim(d):
    d["dim"] = 0


def _edit_duplicate_element(d):
    d.tables["elements"].rows.append(list(d.tables["elements"].rows[0]))


def _edit_missing_layer(d):
    d.tables["elements"].rows[-1][3] = 7                    # excavated element MSET 7: no L layer 7


def _edit_shell_thickness(d):
    for r in d.tables["elements"].rows:
        if r[0] == 2:
            r[6] = 0.0


def _edit_mx_row(d):
    for r in d.tables["matrices"].rows:
        if r[1] == "R" and r[2] == 12:
            r[4] = 1.0                                       # row 12 may hold one term only


@pytest.mark.parametrize("edit,msg", [(_edit_gravity, "Error 1: acceleration of gravity"),
                                      (_edit_dim, "1D analysis is not available"),
                                      (_edit_duplicate_element, "element 1 of group 1 is defined twice"),
                                      (_edit_missing_layer, "Error 19: soil layer 7"),
                                      (_edit_shell_thickness, "shell thickness 0.0 must be > 0"),
                                      (_edit_mx_row, "R row 12: at most 1 terms")])
def test_input_error_paths(tmp_path, edit, msg):
    d = mixed_house_deck()
    edit(d)
    rc, out = run_house(tmp_path, d)
    assert rc == 1 and msg in out, out[-1500:]
    assert "unexpected failure" not in out and not (tmp_path / "m.N4").exists()


def test_soilmode_and_gravity_warnings(tmp_path):
    d, _ = soil_block_deck()
    sit = decks.read(write_site_deck(tmp_path), "SITE")
    sit["soilmode"] = 1
    decks.write(tmp_path / "m.sit", sit)
    rc, out = run_house(tmp_path, d, site=False)
    assert rc == 0 and "strain-compatible FILE88 properties" in out
    d, _ = soil_block_deck()
    d["gravity"] = 12.0
    rc, out = run_house(tmp_path / "g", d)
    assert rc == 0 and "differs by more than 5 %" in out
