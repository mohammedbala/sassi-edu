"""SYMM half / quarter models through HOUSE, ANALYS and CHECK (D-ANL-12, spec 07 9.2.39, G-17): the symmetry
boundary conditions and FILE4 data of HOUSE, the image-sum impedance and the input-symmetry test of
ANALYS, the option combinations refused with SYMM (incoherency, wave passage, multiple excitation,
global impedance), restarts, the node optimizer and the CHECK geometry rules."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from sassi.core import symmetry as SYM
from sassi.io import decks
from sassi.io.container import read_container
from sassi.modules import analys
from sassi.modules.base import ModuleError
from sassi.prep import Interpreter
from sassi.prep.check import run_check
from sassi.verify import builders as B
from sassi.verify.problems.vp_twod_symm import (T2_SITE_HS, T2_SITE_LAYERS, compare_reduced, symmetric_model)

FS = B.FrequencySet.fourier(0.01, 1024, [10, 41])


def _site() -> B.Site:
    return B.layered_site(T2_SITE_LAYERS, T2_SITE_HS)


def _run(wd, part, types, wave="SV", embedded=False, analys_params=None, check=True):
    site = _site()
    mdl = symmetric_model(site, part, types, embedded)
    B.run_soil(wd, "m", site, FS, layer=mdl.layer, rad=mdl.rad, wave=wave)
    f4 = B.run_house(wd, "m", mdl)
    rc = B.run_analys(wd, "m", FS, check=check, **(analys_params or {}))
    return mdl, f4, rc


def _nodes_xyz(f4):
    return {int(n): np.asarray(p, float) for n, p in zip(np.asarray(f4["node_id"]), np.asarray(f4["node_xyz"]))}


@pytest.fixture(scope="module")
def full_x(tmp_path_factory):
    wd = tmp_path_factory.mktemp("symm_full")
    mdl, f4, _ = _run(wd, "full", (1, 0))
    return wd, mdl, f4, B.read_file8(wd)


# ---------------------------------------------------------------------------------------------- HOUSE
def test_house_applies_symmetry_conditions_and_writes_file4_data(tmp_path, full_x):
    mdl, f4, _ = _run(tmp_path, "half", (SYM.ANTISYMMETRIC, SYM.SYMMETRIC))
    out = B.listing(tmp_path, "m", "HOUSE")
    assert "Symmetry planes (SYMM, D-ANL-12): half model" in out and "antisymmetry" in out
    assert "FIXEDINT" not in out                      # symmetry conditions are not user fixities
    planes = SYM.planes_from_array(f4["x_symm"])
    assert [(p.no, p.type, p.axis, p.coord) for p in planes] == [(1, 1, 0, 0.0)]
    xyz = np.asarray(f4["x_int_xyz"])
    cons = np.asarray(f4["x_int_symm"]) != 0
    on = np.abs(xyz[:, 0]) < 1e-9
    assert np.array_equal(cons[on], np.tile([False, True, True], (on.sum(), 1))) and not cons[~on].any()
    assert np.all(np.asarray(f4["int_eq"])[cons] == -1)
    eq = {(int(n), int(d)) for n, d in zip(f4["eq_node"], f4["eq_dof"])}
    c = mdl["centre"]
    assert (c, 1) in eq and (c, 5) in eq and not {(c, 2), (c, 3), (c, 4)} & eq     # antisymmetry: UY UZ ROTX
    # the planes enter the restart hash of the interaction set
    h_half = read_container(tmp_path / "FILE90").meta["int_hash"]
    d = decks.read(tmp_path / "m.hou", "HOUSE")
    d.tables["symm"].rows = []
    B.write_deck(tmp_path, "m", d)
    B.run("HOUSE", tmp_path, "m")
    assert read_container(tmp_path / "FILE90").meta["int_hash"] != h_half
    assert "x_symm" not in read_container(tmp_path / "m.N4")


def test_house_symmetry_errors(tmp_path):
    site = _site()
    mdl = symmetric_model(site, "full", (1, 0), False)
    hb = mdl.house
    hb.symmetry(1, 1, [hb.find(0.0, 0.0, 0.0), hb.find(0.0, 4.0, 0.0), mdl["top"]])   # full model + plane
    B.run_soil(tmp_path, "m", site, FS, layer=0, rad=mdl.rad)
    B.write_deck(tmp_path, "m", mdl.deck("m"))
    assert B.run("HOUSE", tmp_path, "m", check=False) == 1
    assert "nodes on both sides" in B.listing(tmp_path, "m", "HOUSE")
    half = symmetric_model(site, "half", (1, 0), False)
    hb = half.house
    hb.symmetry(1, 1, [hb.find(0.0, 0.0, 0.0), hb.find(2.0, 2.0, 0.0), half["top"]])  # oblique plane
    B.write_deck(tmp_path, "m", half.deck("m"))
    assert B.run("HOUSE", tmp_path, "m", check=False) == 1
    assert "not parallel to the XZ or YZ plane" in B.listing(tmp_path, "m", "HOUSE")


@pytest.mark.parametrize("params,msg", [
    (dict(coh=1, wpass=1), "incoherent motion and wave passage (HOUSE <coh>/<wpass>) not allowed with SYMM planes"),
    (dict(wpass=1), "wave passage (HOUSE <wpass>) not allowed with SYMM"),
    (dict(wpass=1, me=1), "wave passage and multiple excitation (HOUSE <wpass>/<me>) not allowed with SYMM")])
def test_house_refuses_incoherency_with_symm(tmp_path, params, msg):
    site = _site()
    half = symmetric_model(site, "half", (1, 0), False)
    B.run_soil(tmp_path, "m", site, FS, layer=0, rad=half.rad)
    d = FS.fill(half.deck("m"))
    for k, v in params.items():
        d[k] = v
    B.write_deck(tmp_path, "m", d)
    assert B.run("HOUSE", tmp_path, "m", check=False) == 1
    assert msg in B.listing(tmp_path, "m", "HOUSE")
    assert not (tmp_path / "FILE77").exists()


# --------------------------------------------------------------------------------------------- ANALYS
def test_half_and_quarter_equal_full_model(tmp_path, full_x):
    _, _, f4f, f8f = full_x
    for part in ("half", "quarter"):
        wd = tmp_path / part
        wd.mkdir()
        _, f4, _ = _run(wd, part, (SYM.ANTISYMMETRIC, SYM.SYMMETRIC))
        err, ndof = compare_reduced((f8f, _nodes_xyz(f4f)), (B.read_file8(wd), _nodes_xyz(f4)))
        assert err < 1e-8 and ndof > 20
        out = B.listing(wd, "m", "ANALYS")
        assert f"{'half' if part == 'half' else 'quarter'} model" in out and "image terms" in out
        assert B.read_file8(wd).meta["x_symm"]


def test_wrong_symmetry_type_is_refused(tmp_path):
    """X input with a *symmetric* plane normal to X: the free field is antisymmetric (D-ANL-12)."""
    _, _, rc = _run(tmp_path, "half", (SYM.SYMMETRIC, SYM.SYMMETRIC), check=False)
    out = B.listing(tmp_path, "m", "ANALYS")
    assert rc == 1 and "does not have the symmetry of SYMM plane 1: x = 0" in out
    assert "it is antisymmetry: use SYMM type 1" in out and not (tmp_path / "FILE8").exists()


def test_inclined_wave_across_the_plane_is_refused(tmp_path):
    site = _site()
    half = symmetric_model(site, "half", (1, 0), False)
    B.run_soil(tmp_path, "m", site, FS, layer=0, rad=half.rad, waves=[(2, 1, 1.0, 1.0, 30.0)], cm=0)
    B.run_house(tmp_path, "m", half)
    assert B.run_analys(tmp_path, "m", FS, check=False) == 1
    assert "neither symmetric nor antisymmetric" in B.listing(tmp_path, "m", "ANALYS")


def test_global_impedance_refused_with_symm(tmp_path):
    _, _, rc = _run(tmp_path, "half", (1, 0), analys_params=dict(impe=1), check=False)
    assert rc == 1 and "global impedance (<impe> = 1) is not allowed with SYMM" in B.listing(tmp_path, "m", "ANALYS")


def test_model_rules_refuse_factored_input_with_symm():
    md = SimpleNamespace(symm=[SYM.SymmetryPlane(1, 1, 0, 0.0)], dim=2)
    with pytest.raises(ModuleError, match="not allowed with SYMM planes"):
        analys.model_rules(SimpleNamespace(factored=True, impe=0), md)
    md2 = SimpleNamespace(symm=[], dim=1)
    with pytest.raises(ModuleError, match="need a 3D model"):
        analys.model_rules(SimpleNamespace(factored=True, impe=0), md2)
    analys.model_rules(SimpleNamespace(factored=False, impe=0), md)          # coherent: allowed


def test_restart_with_symm_equals_initiation(tmp_path):
    """New Structure restart (COOX reuse) of a half model equals its initiation (VP-22 for SYMM)."""
    _, _, _ = _run(tmp_path, "half", (1, 0), analys_params=dict(save=1))
    H0 = np.asarray(B.read_file8(tmp_path)["H"]).copy()
    assert B.run_analys(tmp_path, "m", FS, mode=1) == 0
    np.testing.assert_allclose(np.asarray(B.read_file8(tmp_path)["H"]), H0, rtol=0, atol=1e-12 * np.abs(H0).max())


def test_optimizer_keeps_the_symmetry_planes(tmp_path):
    site = _site()
    half = symmetric_model(site, "half", (1, 0), False)
    B.run_soil(tmp_path, "m", site, FS, layer=0, rad=half.rad)
    B.run_house(tmp_path, "m", half)
    B.run_analys(tmp_path, "m", FS)
    f4a, f8a = read_container(tmp_path / "m.N4"), B.read_file8(tmp_path)
    d = decks.read(tmp_path / "m.hou", "HOUSE")
    d["optimize"] = 1
    B.write_deck(tmp_path, "m", d)
    B.run("HOUSE", tmp_path, "m")
    B.run_analys(tmp_path, "m", FS)
    f4b = read_container(tmp_path / "m.N4")
    assert "x_symm" in f4b and int(f4b.meta.get("x_optimized", 0)) == 1
    err, _ = compare_reduced((f8a, _nodes_xyz(f4a)), (B.read_file8(tmp_path), _nodes_xyz(f4b)))
    assert err < 1e-10


# ---------------------------------------------------------------------------------------------- CHECK
BASE = """
N,1,0,0,0
N,2,0,4,0
N,3,4,0,0
N,4,4,4,0
N,5,0,0,6
N,9,-4,0,0
N,6,4,0,6
INT,1,4,1,1
L,1,5,1.9,300,150,0.05,0.05
L,2,10,2.0,500,250,0.04,0.04
TOPL,1,1
M,1,3e7,0.2,24,0.05,0.05
GROUP,1,SHELL
E,1,1,3,4,2
THICK,1,1,1,0.5
HOUSE,9.81,0,0,2,0,0,0,0,0
"""


def _check(extra: str):
    ui = Interpreter()
    ui.run_text(BASE + extra)
    rep = run_check(ui.model, modules=["HOUSE"])
    return [m.line() for m in rep.errors("HOUSE")]


def test_check_symmetry_geometry():
    assert not [e for e in _check("SYMM,1,1,1,2,5\n") if "SYMM" in e or "EDU-26" in e]
    errs = _check("SYMM,1,1,1,4,5\n")                          # oblique plane through (0,0,0), (4,4,0)
    assert any("not parallel to the XZ or YZ plane" in e for e in errs)
    errs = _check("GROUP,1\nE,2,9,1,2,2\nSYMM,1,1,1,2,5\n")    # a shell on the negative side of x = 0
    assert any("nodes on both sides of the plane" in e for e in errs)
    errs = _check("SYMM,1,1,1,2,5\nSYMM,2,0,3,4,6\n")          # x = 0 and x = 4: parallel planes
    assert any("parallel" in e and "orthogonal" in e for e in errs)
    errs = _check("SYMM,1,1,1,2,5\nWPASS,1e9,0,1\nHOUSE,9.81,0,0,2,0,0,1,0,0\n")
    assert any("wave passage not allowed with SYMM" in e for e in errs)
