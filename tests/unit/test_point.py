"""Unit tests of the POINT module (POINT3 axisymmetric core, POINT2 plane-strain strip)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sassi.core import axisym, greens_tlm, strip2d, tlm
from sassi.core.flexibility import flexibility_matrix
from sassi.io import decks
from sassi.io.files import SCHEMAS, read_container, validate
from sassi.modules.base import run_module

G = 9.81
LAYERS = [(0.5, 150, 300, 1.9, 0.05)] * 8 + [(1.0, 250, 500, 2.0, 0.04)] * 6


def write_site(tmp_path: Path, nl=0, fnum=(1, 2), df=2.5, layers=LAYERS):
    d = decks.new("SITE")
    d["gravity"] = G
    d["nl"] = nl
    d["df"] = df
    d["mode2"] = 0
    for i, (h, vs, vp, rho, b) in enumerate(layers):
        d.table("layers").append([i + 1, h, rho * G, vp, vs, b, b])
    d.table("halfspace").append([99, 0.0, 2.1 * G, 800.0, 400.0, 0.03, 0.03])
    for n in fnum:
        d.table("freqs").append([n])
    decks.write(tmp_path / "m.sit", d)
    assert run_module("SITE", "m", tmp_path) == 0


def run_point(tmp_path: Path, layer=2, rad=0.9, dim=2, fnum=(1, 2), df=2.5, **params):
    d = decks.new("POINT")
    d["layer"] = layer
    d["rad"] = rad
    d["dim"] = dim
    d["df"] = df
    for k, v in params.items():
        d[k] = v
    for n in fnum:
        d.table("freqs").append([n])
    decks.write(tmp_path / "m.poi", d)
    rc = run_module("POINT", "m", tmp_path)
    return rc, (tmp_path / "m_POINT.out").read_text()


def test_point3_file3_schema_and_shapes(tmp_path):
    write_site(tmp_path, nl=20)
    rc, out = run_point(tmp_path, layer=2)
    assert rc == 0, out
    f3 = read_container(tmp_path / "FILE3", "FILE3")
    assert validate(f3) == []
    arrays, meta = SCHEMAS["FILE3"]
    assert set(arrays) <= set(f3.arrays) and set(meta) <= set(f3.meta)
    Nf = len(LAYERS) + 1 + 20
    nF, nL = 2, 3
    assert f3["alpha0"].shape == (nF, nL, 3 * Nf) and f3["alpha1"].shape == (nF, nL, 3 * Nf)
    assert f3["axis0"].shape == (nF, nL, nL, 3) and f3["phix_obs"].shape == (nF, nL, 2 * Nf)
    assert f3["phiy_obs"].shape == (nF, nL, Nf) and f3["kL"].shape == (nF, Nf)
    np.testing.assert_array_equal(f3["load_iface"], [1, 2, 3])
    assert f3.meta["R0"] == 0.9 and f3.meta["dim"] == 2
    # axis constraints: mu = 0 -> u_rho = u_theta = 0; mu = 1 -> u_rho = u_theta, u_z = 0
    assert np.all(f3["axis0"][..., :2] == 0) and np.all(f3["axis1"][..., 2] == 0)
    np.testing.assert_array_equal(f3["axis1"][..., 0], f3["axis1"][..., 1])
    assert "Self-flexibility" in out


@pytest.mark.parametrize("base_nl", [0, 20])
def test_point3_far_field_matches_kausel_series(tmp_path, base_nl):
    """FE core + transmitting boundary vs the exact point-load series of the same discrete medium
    (R1 V8 setup is VP-09; here a shallow 10 m site near resonance, r = 32 R0): <= 0.5 % on all
    components."""
    write_site(tmp_path, nl=base_nl, fnum=(2,), df=2.5)
    rc, out = run_point(tmp_path, layer=2, rad=0.5, fnum=(2,))
    assert rc == 0, out
    f2 = read_container(tmp_path / "FILE2")
    f3 = read_container(tmp_path / "FILE3")
    md = tlm.modes_from_file2(f2, 0)
    xy = np.array([[0, 0], [0, 0], [16.0, 0], [16.0, 0]])
    F = flexibility_matrix(f3, 0, xy, [1, 3, 1, 3], symmetrize=False)
    for jl, n in ((0, 1), (1, 3)):
        for io, m in ((2, 1), (3, 3)):
            g = greens_tlm.point_load(md, m - 1, n - 1, 16.0)
            blk = F[3 * io:3 * io + 3, 3 * jl:3 * jl + 3]
            got = {"u": blk[0, 0], "v": blk[1, 1], "w": blk[2, 0], "p": blk[0, 2], "q": blk[2, 2]}
            for c in "uvwpq":
                assert abs(got[c] - g[c]) <= 5e-3 * abs(g[c]), (n, m, c)


def test_point3_reciprocity_between_interfaces(tmp_path):
    write_site(tmp_path, nl=20, fnum=(2,))
    rc, out = run_point(tmp_path, layer=3, rad=0.9, fnum=(2,))
    assert rc == 0, out
    f3 = read_container(tmp_path / "FILE3")
    xy = np.array([[0, 0], [2.7, 1.8], [-3.1, 4.4], [5.0, -1.0]])
    F = flexibility_matrix(f3, 0, xy, [1, 4, 3, 2], symmetrize=False)
    assert np.linalg.norm(F - F.T) / np.linalg.norm(F) < 1e-3


def test_point3_core_solution_api():
    """Library-level check of solve_point3: axis values and amplitudes are consistent with the rim."""
    layers = LAYERS
    h = [l[0] for l in layers]
    rho = [l[3] for l in layers]
    Gs, Ms = tlm.layer_moduli(rho, [l[1] for l in layers], [l[2] for l in layers], [l[4] for l in layers],
                              [l[4] for l in layers])
    col = tlm.Column(h=h, rho=rho, G=Gs, M=Ms)
    om = 2 * np.pi * 4.0
    md = tlm.column_modes(col, om)
    sol = axisym.solve_point3(col, md, om, 0.9, [0, 4])
    for mu in (0, 1):
        P, _ = axisym.psi_matrices(md, mu, 0.9, 0.9)
        np.testing.assert_allclose(P @ sol.alpha[mu].T, sol.rim[mu].T, rtol=1e-8, atol=1e-14)
    # vertical load: the response is axisymmetric (no torsional Love part)
    nR = md.kR.size
    assert np.abs(sol.alpha[0][:, nR:]).max() < 1e-8 * np.abs(sol.alpha[0][:, :nR]).max()


@pytest.mark.parametrize("change,message", [
    (dict(layer=-1), "Error 56"),
    (dict(rad=0.0), "Error 57"),
    (dict(dim=3), "<dim>"),
    (dict(fnum=(1, 7)), "not in FILE2"),
    (dict(layer=40), "needs interface"),
])
def test_point_input_errors(tmp_path, change, message):
    write_site(tmp_path)
    rc, out = run_point(tmp_path, **change)
    assert rc == 1 and message in out


def test_point_rigid_base_interface_cannot_be_loaded(tmp_path):
    write_site(tmp_path, nl=0)
    rc, out = run_point(tmp_path, layer=len(LAYERS))
    assert rc == 1 and "rigid base" in out


def test_point_requires_file2(tmp_path):
    rc, out = run_point(tmp_path)
    assert rc == 1 and "FILE2 missing" in out


def test_point_data_check_only(tmp_path):
    write_site(tmp_path)
    rc, out = run_point(tmp_path, opmode=1)
    assert rc == 0 and not (tmp_path / "FILE3").exists()


def test_point_uses_all_file2_frequencies_when_deck_has_none(tmp_path):
    write_site(tmp_path, fnum=(1, 2, 3))
    rc, out = run_point(tmp_path, fnum=())
    assert rc == 0, out
    np.testing.assert_array_equal(read_container(tmp_path / "FILE3")["fnum"], [1, 2, 3])


# ------------------------------------------------------------------ POINT2 (P1)
def test_point2_far_field_matches_line_load(tmp_path):
    """POINT2 strip + Waas boundaries vs Kausel's line-load functions (R1 V7), |x| = 6 R0: <= 1 %."""
    write_site(tmp_path, nl=0, fnum=(2,), df=2.0)
    rc, out = run_point(tmp_path, layer=2, rad=1.0, dim=1, fnum=(2,), df=2.0)
    assert rc == 0, out
    f2 = read_container(tmp_path / "FILE2")
    f3 = read_container(tmp_path / "FILE3")
    assert validate(f3) == [] and f3.meta["dim"] == 1
    md = tlm.modes_from_file2(f2, 0)
    xy = np.array([[0.0, 0.0], [6.0, 0.0], [-6.0, 0.0]])
    F = flexibility_matrix(f3, 0, xy, [3, 1, 1], symmetrize=False)
    for io, x in ((1, 6.0), (2, -6.0)):
        g = greens_tlm.line_load(md, 0, 2, x)
        blk = F[3 * io:3 * io + 3, 0:3]
        got = {"xx": blk[0, 0], "zx": blk[2, 0], "xz": blk[0, 2], "zz": blk[2, 2], "yy": blk[1, 1]}
        for c, val in got.items():
            assert abs(val - g[c]) <= 1e-2 * abs(g[c]), (x, c)
    # odd coupling terms change sign with x, even terms do not
    assert F[5, 0] == pytest.approx(-F[8, 0]) and F[3, 0] == pytest.approx(F[6, 0])


def test_point2_core_matches_prototype_accuracy():
    layers = [(1.0, 150, 300, 1.9, 0.05)] * 6 + [(2.0, 300, 600, 2.0, 0.05)] * 6
    h = [l[0] for l in layers]
    rho = [l[3] for l in layers]
    Gs, Ms = tlm.layer_moduli(rho, [l[1] for l in layers], [l[2] for l in layers], [l[4] for l in layers],
                              [l[4] for l in layers])
    col = tlm.Column(h=h, rho=rho, G=Gs, M=Ms)
    om = 2 * np.pi * 4.0
    md = tlm.column_modes(col, om, tlm.MASS_CONSISTENT)
    sol = strip2d.solve_point2(col, md, om, 1.0, [3], tlm.MASS_CONSISTENT)
    e = strip2d.exterior_2d(md.kR, md.kL, md.phix[[0]], md.phiz[[0]], md.phiy[[0]], sol.alpha_x, sol.alpha_z,
                            sol.alpha_y, [5.0, 12.0], 1.0)
    for i, x in enumerate((5.0, 12.0)):
        g = greens_tlm.line_load(md, 0, 3, x)
        for c in ("xx", "zx", "xz", "zz", "yy"):
            assert abs(e[c][i, 0, 0] - g[c]) <= 3e-3 * abs(g[c])


def test_vp09_tolerance_plan_follows_requirements(tmp_path):
    """VP-09 (requirements §6.3, R1 V8, D-W1-04): all components are checked at 3.3 R0 (2.35 %) and beyond 13 R0
    (0.25 %); only at 6.7 R0 are the minor components (< 10 % of the dominant one) exempt (0.8 %)."""
    import collections

    from sassi.verify import run_problem
    res = run_problem("VP-09", tmp_path)
    by_dist = collections.defaultdict(list)
    for c in res.checks:
        by_dist[c.quantity.split("r = ")[1].split(" R0")[0]].append(c)
    n_all = 2 * 2 * 5                      # 2 load x 2 observation interfaces x (u, v, w, p, q)
    assert [len(by_dist[d]) for d in ("3.3", "13.3", "27.8")] == [n_all] * 3
    assert 0 < len(by_dist["6.7"]) < n_all and not any("minor" in c.quantity for c in by_dist["6.7"])
    assert {c.tolerance for c in by_dist["3.3"]} == {0.0235} and {c.tolerance for c in by_dist["6.7"]} == {0.008}
    assert all(c.passed for d in ("3.3", "6.7", "13.3", "27.8") for c in by_dist[d])
