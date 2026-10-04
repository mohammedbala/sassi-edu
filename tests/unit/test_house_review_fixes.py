"""Unit tests of the HOUSE fixes made after the adversarial review of work package house_force.

* 2D interaction set: only UX, UZ are interaction equations (requirements 4.1; FILE4 meta int_dofs);
* EDU-10 frequency step falls back to 1/(delt NFFT) when the deck df is 0 (G-06);
* EDU-26 is an error in HOUSE, as in CHECK (D-PNT-01);
* EXCSTRCHK owners exclude buried shells (ETYPE 2), as sassi/prep/check.py:excstrchk and spec 09 3.1;
* excavation numbering/grouping rules R5/R6 (spec 05b section 1.5, D-HOU-04): warnings + FILE4 meta.
"""
from __future__ import annotations

import numpy as np
import pytest

from sassi.core import house_lib as hl
from sassi.io import decks
from sassi.io.files import read_container
from sassi.modules import house
from sassi.verify.problems.vp_house import add_element, add_node, house_deck, run_house, soil_block_deck


# ======================================================================================
# 2D interaction equations
# ======================================================================================
def _plane_with_column(dim: int = 1):
    """Two excavated PLANE elements (one layer, 2 m) and a BEAMS column on the top-middle node."""
    d = house_deck(dim=dim)
    d.table("groups").append([1, 4, "soil"])
    ids = {}
    n = 0
    for iz, z in enumerate((-2.0, 0.0)):
        for ix, x in enumerate((0.0, 1.0, 2.0)):
            n += 1
            ids[(ix, iz)] = n
            add_node(d, n, x, 0.0, z)
            d.table("interaction").append([n])
    for ix in range(2):
        add_element(d, 1, ix + 1, [ids[(ix, 0)], ids[(ix + 1, 0)], ids[(ix + 1, 1)], ids[(ix, 1)]], etype=2, mat=1)
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("beamprops").append([1, 0.25, 0.0, 0.0, 0.01, 0.004, 0.006])
    d.table("groups").append([2, 2, "column"])
    add_node(d, 50, 1.0, 0.0, 3.0)
    add_node(d, 99, 5.0, 0.0, 0.0, fix=(1,) * 6)
    add_element(d, 2, 1, [ids[(1, 1)], 50, 99])
    return d, ids


def test_2d_interaction_equations_are_ux_uz_only(tmp_path):
    d, ids = _plane_with_column()
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out
    f4 = read_container(tmp_path / "m.N4")
    top_mid = list(f4["int_node"]).index(ids[(1, 1)])
    # the column node has a UY equation (6-DOF BEAMS element) ...
    k = np.flatnonzero((f4["eq_node"] == ids[(1, 1)]) & (f4["eq_dof"] == 2))
    assert k.size == 1
    # ... but it is structural only: the interaction set is UX, UZ of every node
    assert np.all(f4["int_eq"][:, 1] == -1)
    assert np.all(f4["int_eq"][:, [0, 2]] >= 0) and f4["int_eq"][top_mid, 0] >= 0
    assert f4.meta["int_dofs"] == [1, 3]
    assert "6 interaction nodes, 12 interaction equations" in out
    assert house.interaction_dofs(1) == (1, 3) and house.interaction_dofs(2) == (1, 2, 3)


def test_3d_interaction_meta(tmp_path):
    d, _ = soil_block_deck()
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out
    f4 = read_container(tmp_path / "m.N4")
    assert f4.meta["int_dofs"] == [1, 2, 3] and np.all(f4["int_eq"] >= 0)


# ======================================================================================
# EDU-10 frequency step
# ======================================================================================
def test_house_frequency_step_rule():
    d = decks.new("HOUSE")
    d["df"], d["delt"], d["nft"] = 0.25, 0.01, 1024
    assert house.house_frequency_step(d) == 0.25                     # resolved df wins
    d["df"] = 0.0
    assert house.house_frequency_step(d) == pytest.approx(1.0 / 10.24, rel=1e-15)
    d["delt"] = 0.0
    assert house.house_frequency_step(d) == 0.0                      # nothing given: check skipped


@pytest.mark.parametrize("df,delt,nft", [(0.0, 0.005, 4096), (1.0 / (0.005 * 4096), 0.005, 4096)])
def test_edu10_with_and_without_resolved_df(tmp_path, df, delt, nft):
    """Element heights 3 m (layer 2, Vs 300) and 2 m (layer 1, Vs 200): f_pass = 20 Hz; frequency
    number 800 at df = 1/20.48 Hz is 39 Hz -> EDU-10 for all 8 excavated elements in both cases."""
    d, _ = soil_block_deck()
    d["df"], d["delt"], d["nft"] = df, delt, nft
    d.table("freqs").append([800])
    rc, out = run_house(tmp_path, d)
    assert rc == 0, out
    assert "EDU-10: 8 excavated elements are taller than Vs/(5 f_cut) at f_cut = 39.06 Hz" in out


# ======================================================================================
# EDU-26 severity (D-PNT-01)
# ======================================================================================
def test_edu26_plane_in_3d_is_an_error(tmp_path):
    d, _ = _plane_with_column(dim=2)
    rc, out = run_house(tmp_path, d)
    assert rc == 1 and "EDU-26: PLANE elements in a 3D model" in out and "*** ERROR" in out
    assert not (tmp_path / "m.N4").exists() and not (tmp_path / "COOSK").exists()


def test_edu26_solid_in_2d_is_an_error(tmp_path):
    d, _ = soil_block_deck()
    d["dim"] = 1
    rc, out = run_house(tmp_path, d)
    assert rc == 1 and "EDU-26: SOLID/SHELL/TSHELL elements in a 2D model" in out
    assert not (tmp_path / "m.N4").exists()


# ======================================================================================
# EXCSTRCHK and buried shells
# ======================================================================================
def _buried_shell_deck(etype: int):
    """FI block (only the boundary nodes interact) with a SHELL on the plane z = -2 that uses the
    interior excavation node 14 (the centre of the 2 x 2 x 2 block)."""
    d, nid = soil_block_deck()
    d["imp"] = 2
    centre = nid[(1, 1, 1)]
    d.tables["interaction"].rows = [r for r in d.tables["interaction"].rows if r[0] != centre]
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("groups").append([3, 3, "buried wall" if etype == 2 else "slab"])
    add_element(d, 3, 1, [centre, nid[(2, 1, 1)], nid[(2, 2, 1)], nid[(1, 2, 1)]], etype=etype, mat=1, thick=0.2)
    return d, centre


def test_excstrchk_ignores_buried_shells(tmp_path):
    d, centre = _buried_shell_deck(etype=2)
    rc, out = run_house(tmp_path / "buried", d)
    assert rc == 0, out
    assert "EXCSTRCHK: excavation interior node" not in out
    assert f"buried shell (ETYPE 2) nodes that are not interaction nodes: {centre}" in out
    # the same shell as structure (ETYPE 1) is a basement element: EXCSTRCHK error (spec 09 3.1)
    d, centre = _buried_shell_deck(etype=1)
    rc, out = run_house(tmp_path / "struct", d)
    assert rc == 1
    assert f"EXCSTRCHK: excavation interior node {centre} is shared with the structure (SHELL element 1 group 3)" in out


# ======================================================================================
# R5 / R6 (D-HOU-04)
# ======================================================================================
def test_excavation_layering_helper():
    ok = [(1, 1, 1, 1), (1, 2, 1, 1), (2, 1, 2, 2), (2, 2, 2, 2)]
    assert hl.excavation_layering(ok) == {"R5": [], "R6": []}
    assert hl.excavation_layering([]) == {"R5": [], "R6": []}
    # one group, numbered bottom-up across two layers: R5 and R6 (span)
    bad = [(1, 1, 2, 2), (1, 2, 2, 2), (1, 3, 1, 1), (1, 4, 1, 1)]
    res = hl.excavation_layering(bad)
    assert len(res["R5"]) == 1 and "element 3 (layer 1) follows element 2 (layer 2)" in res["R5"][0]
    assert res["R6"] == ["excavation group 1 spans 2 layers (1, 2): use one group per embedment layer"]
    # one group, numbered top-down across two layers: R5 met, R6 (span) only
    res = hl.excavation_layering([(1, 1, 1, 1), (1, 2, 2, 2)])
    assert res["R5"] == [] and len(res["R6"]) == 1
    # groups numbered from the foundation up
    res = hl.excavation_layering([(1, 1, 2, 2), (2, 1, 1, 1)])
    assert res["R5"] == [] and res["R6"] == ["excavation groups are not numbered from the surface down: "
                                             "group 2 (layer 1) follows group 1 (layer 2)"]
    # one layer split over two groups
    res = hl.excavation_layering([(1, 1, 1, 1), (2, 1, 1, 1), (3, 1, 2, 2)])
    assert res["R6"] == ["embedment layer 1 is split over excavation groups 1, 2: use one group per "
                         "embedment layer"]


def test_r5_r6_in_house_listing_and_file4(tmp_path):
    # the VP block follows the rules: no R5/R6 warnings, empty lists in FILE4 meta
    d, _ = soil_block_deck()
    rc, out = run_house(tmp_path / "ok", d)
    assert rc == 0, out
    assert "R5/R6 satisfied" in out and "WARNING: R5" not in out and "WARNING: R6" not in out
    assert read_container(tmp_path / "ok" / "m.N4").meta["excav_layering"] == {"R5": [], "R6": []}
    # the same block as one group numbered bottom-up: R5 and R6 warnings, the run continues
    d, _ = soil_block_deck()
    rows = d.tables["elements"].rows
    lower = [r for r in rows if r[0] == 2]                  # group 2 = layer 2 (z -5 .. -2)
    upper = [r for r in rows if r[0] == 1]
    for k, r in enumerate(lower + upper):
        r[0], r[1] = 1, k + 1
    d.tables["groups"].rows = [r for r in d.tables["groups"].rows if r[0] == 1]
    rc, out = run_house(tmp_path / "bad", d)
    assert rc == 0, out
    assert "WARNING: R5: excavation group 1: element 5 (layer 1) follows element 4 (layer 2)" in out
    assert "WARNING: R6: excavation group 1 spans 2 layers (1, 2)" in out and "D-HOU-04" in out
    lay = read_container(tmp_path / "bad" / "m.N4").meta["excav_layering"]
    assert len(lay["R5"]) == 1 and len(lay["R6"]) == 1
    # no excavated soil: nothing to check
    d = house_deck()
    d.table("materials").append([1, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05])
    d.table("beamprops").append([1, 0.25, 0.0, 0.0, 0.01, 0.004, 0.006])
    d.table("groups").append([1, 2, "column"])
    add_node(d, 1, 0, 0, 0.0, fix=(1,) * 6)
    add_node(d, 2, 0, 0, 3.0)
    add_node(d, 3, 1, 0, 0.0, fix=(1,) * 6)
    add_element(d, 1, 1, [1, 2, 3])
    rc, out = run_house(tmp_path / "none", d)
    assert rc == 0, out
    assert "Excavation groups" not in out
    assert read_container(tmp_path / "none" / "m.N4").meta["excav_layering"] == {"R5": [], "R6": []}
