"""Unit tests of the shared model/deck builders (sassi/verify/builders.py)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.io import decks
from sassi.io.container import read_container
from sassi.verify import builders as B

FS = B.FrequencySet.fourier(0.01, 1024, [4, 20])


def _site():
    return B.layered_site([(2.0, 150.0, 300.0, 1.9, 0.05), (3.0, 250.0, 500.0, 2.0, 0.04, 0.03),
                           (5.0, 350.0, 700.0, 2.1, 0.03)], (600.0, 1200.0, 2.2, 0.02))


def test_site_geometry_and_rows():
    s = _site()
    assert s.nI == 4
    np.testing.assert_allclose(s.depths(), [0, 2, 5, 10])
    s.gelev = 3.0
    np.testing.assert_allclose(s.elevations(), [3, 1, -2, -7])
    rows = s.layer_rows()
    assert rows[1] == [2, 3.0, 2.0 * s.gravity, 500.0, 250.0, 0.03, 0.04]       # no thick weight vp vs dp ds
    assert rows[-1][0] == 4 and rows[-1][1] == 0.0
    sub = s.subdivided(1.0)
    assert len(sub.layers) == 2 + 3 + 5
    np.testing.assert_allclose(sub.depths()[-1], 10.0)
    lay = B.SoilLayer.from_nu(1.0, 100.0, 2.0, nu=1.0 / 3.0)
    assert math.isclose(lay.vp, 200.0) and lay.dp == lay.beta_s
    u = B.uniform_site(20.0, 4, 200.0, 2.0)
    assert len(u.layers) == 4 and u.halfspace.vs == 200.0


def test_frequency_sets():
    fs = B.FrequencySet.fourier(0.005, 4096, [8, 2, 4])
    assert fs.fnum == [2, 4, 8] and math.isclose(fs.df, 1 / 20.48)
    h = B.FrequencySet.harmonic(1e-3, range(4995, 5006))
    assert h.nft == 16384 and math.isclose(1.0 / (h.delt * h.nft), 1e-3) and h.fstep == 1e-3
    np.testing.assert_allclose(h.subset([5000]).freq, [5.0])
    d = h.fill(decks.new("SITE"))
    assert d["fstep"] == 1e-3 and [r["number"] for r in d.rows("freqs")] == h.fnum
    d2 = h.fill(decks.new("ANALYS"))
    assert d2["df"] == 1e-3


@pytest.mark.parametrize("wave,cm,wopt,wtype", [("SV", 0, 0, 2), ("SH", 1, 1, 4), ("P", 2, 0, 3)])
def test_site_deck_waves(wave, cm, wopt, wtype):
    d = B.site_deck(_site(), FS, wave=wave, cl=2)
    assert (d["cm"], d["wopt"], d["cl"]) == (cm, wopt, 2)
    assert d.rows("waves")[0]["type"] == wtype and len(d.rows("layers")) == 3 and len(d.rows("halfspace")) == 1


def test_decks_reject_unknown_parameters():
    with pytest.raises(KeyError):
        B.analys_deck(FS, bogus=1)
    with pytest.raises(KeyError):
        B.motion_deck(FS, bogus=1)
    d = B.motion_deck(FS, "a.th", [(1, 1, 1, 0, 0, 0, 0, 0)], damp=[0.05], type=1)
    assert d["type"] == 1 and math.isclose(d["df"], FS.df) and len(d.rows("damp")) == 1
    f = B.force_deck([(3, 1, 0.5, 0.1), (4, 2, 1.0)], FS)
    assert f.rows("loads")[1]["arrival"] == 0.0


def _house_ok(tmp_path, mdl, name="m"):
    wd = tmp_path / name
    B.write_deck(wd, "m", B.site_deck(_site(), FS, model="m"))
    f4 = B.run_house(wd, "m", mdl)
    text = B.listing(wd, "m", "HOUSE")
    return f4, text


@pytest.mark.parametrize("method,count", [("FV", 75), ("FSIN", 75 - 9 - 9), ("EVBN", 75 - 9), ("FFV", 75 - 9)])
def test_embedded_box_interaction_sets(tmp_path, method, count):
    """R1 4.5: FV every node; FSIN lateral + bottom; EVBN + top face; FFV + every skip-th plane
    (n_emb = 2, skip 2: the interior plane is level 1, so FFV = EVBN here)."""
    structure = "soil" if method == "FV" else "shell"
    mdl = B.embedded_box(_site(), 3.0, 2, ndiv=4, method=method, structure=structure)
    assert len(mdl.house.interaction) == count
    f4, text = _house_ok(tmp_path, mdl)
    assert len(f4["int_node"]) == count
    z = np.asarray(f4["x_int_xyz"])[:, 2]
    assert np.all(np.diff(z) >= 0)                          # bottom-up interaction order (EDU-21)
    assert "EDU-21" not in text and "EDU-06" not in text and "EXCSTRCHK: excavation interior node" not in text
    assert mdl.layer == 2 and math.isclose(mdl.rad, 0.9 * 1.5)


def test_embedded_box_ffv_skip_and_errors():
    mdl = B.embedded_box(_site(), 3.0, 3, ndiv=2, method="FFV", structure="none", skip=2)
    lev = mdl["levels"]                                     # top first: level 2 is an interior plane
    assert set(lev[2]) <= set(mdl.house.interaction)
    assert not set(lev[1]) <= set(mdl.house.interaction)
    with pytest.raises(ValueError, match="only method 'FV'"):
        B.embedded_box(_site(), 3.0, 2, method="FSIN", structure="soil")
    with pytest.raises(ValueError):
        B.embedded_box(_site(), 3.0, 2, method="XYZ")
    with pytest.raises(ValueError):
        B.embedded_box(_site(), 3.0, 4)
    with pytest.raises(ValueError, match="top centre"):
        B.embedded_box(_site(), 3.0, 1, ndiv=3, structure="none", stick=B.Stick([1.0], [1.0], 1e7, 1.0, 1.0))


def test_zero_ssi_box_has_identical_matrices(tmp_path):
    """structure='soil': K*_s = K*_e and M_s = M_e (VP-16 precondition, D-ELM-12)."""
    mdl = B.embedded_box(_site(), 3.0, 2, ndiv=2, method="FV", structure="soil")
    _house_ok(tmp_path, mdl)
    ck, cm = read_container(tmp_path / "m" / "COOSK"), read_container(tmp_path / "m" / "COOSM")
    assert abs(ck.sparse("Ks") - ck.sparse("Ke")).max() == 0.0
    assert abs(cm.sparse("Ms") - cm.sparse("Me")).max() == 0.0


@pytest.mark.parametrize("shape,rigid,n", [("square", "beams", 25), ("square", "shell", 25), ("disk", "beams", 1 + 3 * 2 * 3)])
def test_surface_rigid_mat(tmp_path, shape, rigid, n):
    mdl = B.surface_rigid_mat(_site(), half_width=4.0, ndiv=4 if shape == "square" else 2, shape=shape, rigid=rigid,
                              mass=10.0)
    assert len(mdl["mat"]) == n and mdl["centre"] in mdl["mat"]
    f4, text = _house_ok(tmp_path, mdl)
    assert len(f4["int_node"]) == n and "EDU-06" not in text
    m91 = read_container(tmp_path / "m" / "FILE91").meta
    np.testing.assert_allclose(m91["mass_structure"], [10.0] * 3)          # rigid links are massless
    assert math.isclose(mdl.rad, (0.9 if shape == "square" else 0.85) * mdl["h"])


def test_surface_mat_errors():
    with pytest.raises(ValueError, match="centre node"):
        B.surface_rigid_mat(_site(), ndiv=3)
    with pytest.raises(ValueError):
        B.surface_rigid_mat(_site(), shape="disk", rigid="shell")
    with pytest.raises(ValueError):
        B.surface_rigid_mat(_site(), shape="hexagon")
    m = B.surface_rigid_mat(_site(), ndiv=3, centre_height=1.0)
    assert m["centre"] not in m["mat"]


def test_stick_on_mat_and_orientation_nodes(tmp_path):
    stick = B.Stick(heights=[3.0, 6.0], masses=[40.0, 30.0], E=3e7, A=1.0, I=0.5, rotary=[2.0, 1.0])
    mdl = B.stick_on_mat(_site(), stick, half_width=2.0, ndiv=2)
    hb = mdl.house
    assert len(mdl["stick"]) == 2 and mdl["top"] == mdl["stick"][-1]
    # orientation (K) nodes are fixed and never re-used as structural nodes
    knodes = set(hb._knodes.values())
    assert knodes and all(all(hb.fixity[k]) for k in knodes) and not knodes & set(mdl["stick"])
    f4, text = _house_ok(tmp_path, mdl)
    m91 = read_container(tmp_path / "m" / "FILE91").meta
    np.testing.assert_allclose(m91["mass_structure"], [70.0] * 3)


@pytest.mark.parametrize("direction", [1, 2])
def test_sdof_elements_equivalent(tmp_path, direction):
    """SPRING, 2-node GENERAL (global) and 3-node GENERAL (local) give the same K* (VP-39)."""
    K = {}
    for element in ("spring", "general", "general3"):
        mdl = B.sdof_on_node(_site(), 1000.0, 1.0, 0.05, element=element, direction=direction)
        _house_ok(tmp_path, mdl, element)
        K[element] = read_container(tmp_path / element / "COOSK").sparse("Ks").toarray()
        eqd = read_container(tmp_path / element / "m.N4")["eq_dof"]
        assert sorted(set(eqd.tolist())) == [1, 2, 3]
    for element in ("general", "general3"):
        np.testing.assert_allclose(K[element], K["spring"], rtol=0, atol=1e-9)
    with pytest.raises(ValueError):
        B.sdof_on_node(_site(), 1.0, 1.0, 0.05, element="general3", direction=3)
    with pytest.raises(ValueError):
        B.sdof_on_node(_site(), 1.0, 1.0, 0.05, element="beam")


def test_house_builder_basics():
    hb = B.HouseBuilder(_site())
    a = hb.node(0, 0, 0)
    assert hb.node(0, 0, 1e-9) == a and hb.find(0, 0, 0) == a
    with pytest.raises(ValueError):
        hb.add_node(1, 1, 1, nid=a)
    m1 = hb.elastic(3e7, 0.2, 2.4, 0.05)
    assert hb.elastic(3e7, 0.2, 2.4, 0.05) == m1
    assert hb.materials[m1 - 1][4] == pytest.approx(2.4 * hb.site.gravity)
    assert hb.soil_material(2) != m1
    g1, g2 = hb.group(2, "x"), hb.group(2, "x")
    assert g1 == g2
    d = hb.deck("q")
    assert d["model"] == "q" and len(d.rows("sitelayers")) == 4 and d.rows("sitelayers")[-1]["thick"] == 0.0


def test_runner_reports_failures(tmp_path):
    B.write_deck(tmp_path, "m", B.analys_deck(FS, model="m"))
    with pytest.raises(B.ChainError, match="ANALYS failed"):
        B.run("ANALYS", tmp_path, "m")
    assert B.run("ANALYS", tmp_path, "m", check=False) == 1


def test_tf_reader(tmp_path):
    mdl = B.surface_rigid_mat(_site(), half_width=2.0, ndiv=2)
    wd = tmp_path / "w"
    B.run_soil(wd, "m", _site(), FS, layer=0, rad=mdl.rad)
    B.run_house(wd, "m", mdl)
    B.run_analys(wd, "m", FS)
    f8 = B.read_file8(wd)
    knode = next(iter(mdl.house._knodes.values()))
    assert np.array_equal(B.tf(f8, knode, 1), np.zeros(2))
    np.testing.assert_allclose(B.tf(f8, mdl["centre"], 1), 1.0, atol=1e-8)
    names = B.run_site_xyz(wd, "m", _site(), FS)
    assert names == ["FILE1X", "FILE1Y", "FILE1Z"]
    assert read_container(wd / "FILE1").meta["cm"] == 0 and read_container(wd / "FILE1Y").meta["cm"] == 1
