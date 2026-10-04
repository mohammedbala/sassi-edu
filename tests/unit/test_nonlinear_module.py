"""Unit tests of the NONLINEAR module (sassi/modules/nonlinear.py; requirements 4.15, spec 05d section 3,
D-NON-01, D-NON-04, D-NON-06 ... D-NON-08) with synthetic relative-displacement histories: the .eql deck,
the state machine (PANEL.NON / SPRING.NON, *_EQL_Matl_Prop.txt), panel kinematics from the corner .THD
files, the new HOUSE deck, the output files, COMB_XYZ_THD and the batch entry point."""
from __future__ import annotations

import io
import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from sassi.conventions import nodal_result_name
from sassi.core import hysteresis as HY
from sassi.io import decks
from sassi.modules import nonlinear as NL
from sassi.verify import builders as B

DT = 0.01
T = np.arange(1000) * DT


def site():
    return B.layered_site([(10.0, 300.0, 600.0, 2.0, 0.05)], (800.0, 1600.0, 2.2, 0.02))


def write_thd(wd: Path, node: int, dof: int, values):
    NL.write_history(wd / nodal_result_name(node, dof, "THD"), np.asarray(values, float), DT)


# ======================================================================================
# A spring model
# ======================================================================================
K_EL, XI_EL, XY, FY = 1000.0, 0.02, 0.01, 10.0


def spring_case(wd: Path, amp: float = 0.03, edf: float = 1.0, elasticd: int = 1) -> NL.EqlDeck:
    hb = B.HouseBuilder(site())
    a = hb.add_node(0.0, 0.0, 0.0, fix=(0, 0, 0, 1, 1, 1))
    b = hb.add_node(0.0, 0.0, 1.0, fix=(0, 1, 1, 1, 1, 1))
    hb.set_interaction([a])
    hb.spring(a, b, (K_EL, 0, 0, 0, 0, 0), XI_EL)
    hb.mass(b, 1.0)
    B.write_deck(wd, "m", hb.deck("m"))
    e = NL.EqlDeck()
    e.params.update(model="m", edf=edf, nonlinopts=NL.OPT_SPRINGS, elasticd=elasticd, gravity=9.81)
    for p, (x, y) in enumerate(((XY, FY), (20 * XY, FY + 0.1 * K_EL * 19 * XY)), start=1):
        e.add("bbc", **{"bbc": 1, "type": 4, "yield": 1, "point": p, "x": x, "y": y})
    e.add("springs", num=1, group=1, elem=1, bbc=1, disp=1, force=4, prop=1, node_i=a, node_j=b, scx=K_EL, scy=0.0,
          scz=0.0, scxx=0.0, scyy=0.0, sczz=0.0, damp=XI_EL, k_el=K_EL)
    NL.write_eql(wd / "m.eql", e)
    write_thd(wd, a, 1, 0.002 * np.sin(2 * np.pi * T))                # base: any common motion cancels
    write_thd(wd, b, 1, 0.002 * np.sin(2 * np.pi * T) + amp * np.sin(2 * np.pi * 2 * T))
    return e


def test_eql_deck_round_trip(tmp_path):
    e = spring_case(tmp_path)
    d = NL.read_eql(tmp_path / "m.eql")
    assert d.params["edf"] == 1.0 and d.params["nonlinopts"] == 2 and d.params["maxit"] == 10
    assert d.rows("springs") == e.rows("springs") and d.rows("bbc") == e.rows("bbc")
    (typ, yi, x, y), = d.backbones().values()
    assert (typ, yi) == (4, 1) and list(x) == [XY, 20 * XY]


def test_spring_state_machine_and_outputs(tmp_path):
    spring_case(tmp_path)
    assert NL.run_nonlinear("m", tmp_path) == 0
    lst = (tmp_path / "m_NONLINEAR.out").read_text()
    assert "ELASTIC run" in lst and "status OK" in lst
    assert (tmp_path / "SPRING.NON").read_text().strip() == "1"
    it, rows = NL.read_props(tmp_path / "SPRING_EQL_Matl_Prop.txt")
    r = rows[1]
    bb = HY.Backbone([XY, 20 * XY], [FY, FY + 0.1 * K_EL * 19 * XY], 1)
    xm = float(np.max(np.abs(0.03 * np.sin(2 * np.pi * 2 * T))))         # sampled peak
    assert it == 1 and r.x_max == pytest.approx(xm, rel=1e-7)
    assert r.ratio == pytest.approx(bb.secant(xm) / K_EL, rel=1e-7)
    assert r.xi_h == pytest.approx(HY.masing_loop_energy(bb, xm) / (2 * math.pi * xm * bb.force(xm)), rel=1e-7)
    assert r.xi_new == pytest.approx(r.xi_h + XI_EL, rel=1e-12)
    # the new HOUSE deck holds k_el * ratio and the new damping; everything else unchanged
    new = decks.read(tmp_path / "m_new.hou", "HOUSE")
    sc = new.rows("springprops")[0]
    assert sc["scx"] == pytest.approx(K_EL * r.ratio) and sc["damp"] == pytest.approx(r.xi_new)
    old = decks.read(tmp_path / "m.hou", "HOUSE")
    assert new.rows("nodes") == old.rows("nodes") and new.rows("elements") == old.rows("elements")
    for f in ("SPRING0001.thd", "SPRING0001.ths", "SPRING0001.crv", "Spring.fmu", NL.CONV_FILE):
        assert (tmp_path / f).exists(), f
    x, dt = NL.read_history(tmp_path / "SPRING0001.thd")
    assert dt == DT and np.allclose(x, 0.03 * np.sin(2 * np.pi * 2 * T), atol=1e-15)
    rows_c = NL.read_convergence(tmp_path)
    assert rows_c[-1]["iteration"] == 0 and rows_c[-1]["converged"] == 0
    # ---- second run: SPRING.NON = 1 -> iteration; the properties of the props file are 'used'
    p = tmp_path / "SPRING_EQL_Matl_Prop.txt"
    text = p.read_text().replace(f"{r.xi_new:>23.16e}", f"{0.123:>23.16e}")
    p.write_text(text)
    shutil.copyfile(tmp_path / "m_new.hou", tmp_path / "m.hou")
    assert NL.run_nonlinear("m", tmp_path) == 0
    lst = (tmp_path / "m_NONLINEAR.out").read_text()
    assert "nonlinear iteration" in lst and "SSI analysis processed: 1" in lst
    assert "0.12300" in lst                                            # the damping 'used' comes from the file
    assert "does not hold the properties" in lst                        # the HOUSE deck has another damping
    it2, rows2 = NL.read_props(p)
    assert it2 == 2
    c = NL.read_convergence(tmp_path)[-1]
    assert c["iteration"] == 1 and c["max_de_pct"] == pytest.approx(0.0, abs=1e-9)
    assert c["max_dxi_pct"] == pytest.approx(100 * abs(0.123 - rows2[1].xi_new), rel=1e-5)


def test_elastic_response_converges_immediately(tmp_path):
    spring_case(tmp_path, amp=0.005, elasticd=1)
    assert NL.run_nonlinear("m", tmp_path) == 0
    c = NL.read_convergence(tmp_path)[-1]
    assert c["converged"] == 1 and c["max_de_pct"] == 0.0
    _, rows = NL.read_props(tmp_path / "SPRING_EQL_Matl_Prop.txt")
    assert rows[1].ratio == 1.0 and rows[1].xi_new == XI_EL


def test_missing_thd_and_inconsistent_state(tmp_path):
    spring_case(tmp_path)
    (tmp_path / nodal_result_name(2, 1, "THD")).unlink()
    assert NL.run_nonlinear("m", tmp_path) == 1
    assert "not found" in (tmp_path / "m_NONLINEAR.out").read_text()
    spring_case(tmp_path)
    (tmp_path / "SPRING.NON").write_text("0\n")
    assert NL.run_nonlinear("m", tmp_path) == 1
    assert "does not contain 1" in (tmp_path / "m_NONLINEAR.out").read_text()


# ======================================================================================
# A wall panel
# ======================================================================================
def panel_case(wd: Path, gamma, rigid=True):
    """4 x 2 shell wall 12 m x 4 m in the XZ plane; corner histories u = gamma(t) z (+ rigid-body motion)."""
    hb = B.HouseBuilder(site())
    nid = {}
    for k in range(3):
        for i in range(5):
            nid[(i, k)] = hb.add_node(3.0 * i, 0.0, 2.0 * k)
    hb.set_interaction([nid[(i, 0)] for i in range(5)])
    mat = hb.material(1, 3.0e7, 0.2, 24.0, 0.04, 0.04)
    g = hb.group(3, "wall")
    for k in range(2):
        for i in range(4):
            hb.shell((nid[(i, k)], nid[(i + 1, k)], nid[(i + 1, k + 1)], nid[(i, k + 1)]), mat, 0.3, group=g)
    B.write_deck(wd, "m", hb.deck("m"))
    corners = (nid[(0, 0)], nid[(4, 0)], nid[(4, 2)], nid[(0, 2)])
    G = 3.0e7 / 2.4
    ga = G * 0.3 * 12.0
    e = NL.EqlDeck()
    e.params.update(model="m", edf=0.8, nonlinopts=NL.OPT_PANELS, elasticd=1, dampcutoff=7.0, gravity=9.81)
    gcr, vcr, vu = 4911.84 / ga, 4911.84, 7527.35
    from sassi.core.panels import bbcgen_curve
    gg, vv, yi = bbcgen_curve(vu, vcr, ga)
    for p, (x, y) in enumerate(zip(gg, vv), start=1):
        e.add("bbc", **{"bbc": 1, "type": 1, "yield": yi, "point": p, "x": x, "y": y})
    e.add("panels", num=1, group=g, bbc=1, disp=1, force=1, mat=mat, mtype=1, val1=3.0e7, val2=0.2, weight=24.0,
          pdamp=0.04, sdamp=0.04, e_el=3.0e7, nu_el=0.2, g_el=G, thick=0.3, length=12.0, height=4.0, area=3.6,
          inertia=0.3 * 12.0 ** 3 / 12.0, bl=corners[0], br=corners[1], tr=corners[2], tl=corners[3], ex=1.0, ey=0.0,
          nelem=8)
    NL.write_eql(wd / "m.eql", e)
    xyz = {n: np.array(hb.nodes[n]) for n in corners}
    tr = np.array([0.01, -0.02, 0.005])
    theta = 1e-3
    for n in corners:
        x, _, z = xyz[n]
        ux = gamma * z + (tr[0] - theta * (z - 2.0) if rigid else 0.0)       # rigid rotation about Y at (6, 2)
        uz = (tr[2] + theta * (x - 6.0)) if rigid else 0.0 * gamma
        uy = (tr[1] + 0.0 * gamma) if rigid else 0.0 * gamma
        write_thd(wd, n, 1, ux * np.ones_like(T) if np.ndim(ux) == 0 else ux)
        write_thd(wd, n, 2, uy * np.ones_like(T) if np.ndim(uy) == 0 else uy)
        write_thd(wd, n, 3, uz * np.ones_like(T) if np.ndim(uz) == 0 else uz)
    return gcr, mat


def test_panel_strain_and_cms_update(tmp_path):
    gamma = 3e-4 * np.sin(2 * np.pi * T) * np.minimum(1.0, T)
    gcr, mat = panel_case(tmp_path, gamma)
    assert NL.run_nonlinear("m", tmp_path) == 0, (tmp_path / "m_NONLINEAR.out").read_text()[-2000:]
    g, _ = NL.read_history(tmp_path / "Panel0001.thd")
    np.testing.assert_allclose(g, gamma, rtol=0, atol=1e-15)          # rigid-body motion removed exactly
    ax, _ = NL.read_history(tmp_path / "Panel0001_AXIAL.thd")
    assert np.max(np.abs(ax)) < 1e-15
    f, _ = NL.read_history(tmp_path / "Panel0001.ths")
    assert np.max(np.abs(f)) > 4911.0                                  # cracked: the force passed V_cr
    _, rows = NL.read_props(tmp_path / "Panel_EQL_Matl_Prop.txt")
    r = rows[1]
    assert r.x_max == pytest.approx(3e-4, rel=1e-3) and r.x_eq == pytest.approx(0.8 * r.x_max)
    assert r.k_new == pytest.approx(3.0e7 * r.ratio) and 0.0 < r.ratio < 1.0
    assert r.xi_new == pytest.approx(min(0.07, 0.04 + r.xi_h))       # elastic + hysteretic, cut off at 7 %
    new = decks.read(tmp_path / "m_new.hou", "HOUSE")
    m = [row for row in new.rows("materials") if row["id"] == mat][0]
    assert m["val1"] == pytest.approx(3.0e7 * r.ratio) and m["val2"] == 0.2
    assert m["pdamp"] == m["sdamp"] == pytest.approx(r.xi_new)
    crv = np.loadtxt(tmp_path / "Panel0001.crv")
    assert crv.shape[1] == 5 and np.all(crv[crv[:, 0] <= gcr, 1] == 1.0)
    fmu = np.loadtxt(tmp_path / "Panel.fmu")
    assert fmu[3] == pytest.approx(r.x_max / gcr)                      # ductility = max|x| / x_cr


def test_panel_options_need_nonext(tmp_path):
    panel_case(tmp_path, 1e-5 * np.sin(2 * np.pi * T))
    e = NL.read_eql(tmp_path / "m.eql")
    e.tables["panels"][0]["force"] = 3
    NL.write_eql(tmp_path / "m.eql", e)
    assert NL.run_nonlinear("m", tmp_path) == 1
    assert "126" in (tmp_path / "m_NONLINEAR.out").read_text()
    e.params["nonext"] = 1
    NL.write_eql(tmp_path / "m.eql", e)
    assert NL.run_nonlinear("m", tmp_path) == 0


# ======================================================================================
# COMB_XYZ_THD and the batch entry
# ======================================================================================
def test_comb_xyz_thd(tmp_path):
    for k, s in enumerate("XYZ", start=1):
        NL.write_history(tmp_path / f"{s}_a.thd", np.full(10, float(k)), 0.01)
    (tmp_path / "c.inp").write_text("* comment\n2\nout1.thd X_a.thd Y_a.thd Z_a.thd\nout2.thd X_a.thd - Z_a.thd\n")
    written, warns = NL.comb_xyz_thd(tmp_path, "c.inp")
    assert written == ["out1.thd", "out2.thd"] and not warns
    assert np.all(NL.read_history(tmp_path / "out1.thd")[0] == 6.0)
    assert np.all(NL.read_history(tmp_path / "out2.thd")[0] == 4.0)
    NL.write_history(tmp_path / "Y_a.thd", np.full(10, 1.0), 0.02)
    with pytest.raises(NL.ModuleError, match="dt"):
        NL.comb_xyz_thd(tmp_path, "c.inp")
    (tmp_path / "bad.inp").write_text("3\nout1.thd X_a.thd\n")
    with pytest.raises(NL.ModuleError, match="announced"):
        NL.comb_xyz_thd(tmp_path, "bad.inp")


def test_batch_main(tmp_path, monkeypatch):
    spring_case(tmp_path)
    monkeypatch.chdir(tmp_path)
    rc = NL.main([], stdin=io.StringIO("m\nm.eql\nm_NONLINEAR.out\n"))
    assert rc == 0 and (tmp_path / "SPRING.NON").exists()
    for k, s in enumerate("XYZ", start=1):
        NL.write_history(tmp_path / f"{s}_b.thd", np.full(4, float(k)), 0.01)
    (tmp_path / NL.COMB_INP).write_text("1\nb.thd X_b.thd Y_b.thd Z_b.thd\n")
    assert NL.main(["COMB_XYZ_THD"]) == 0 and (tmp_path / "b.thd").exists()
