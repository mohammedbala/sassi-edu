"""Adversarial review tests of the 'nonlinear_soil' work package (near-field soil SSI iterations:
sassi.core.nlsoil, HOUSE <nlssi> = 1, STRESS <iter> = 1, COMBXYZSTRAIN, NLSSIRESET, NLSSIITER, CHECK).

The implementer's tests compare the modules with references built from the same helpers
(``near_field_column``, ``strain_compatible``, the FILE4 recovery operators).  These tests check the
package from other angles:

* the FILE74 effective strains of STRESS against an independent recomputation from the FILE8 nodal
  transfer functions alone: a centroid strain operator of the box element written out by hand, the
  convolution with the ground displacement ``-g A / w^2`` and ESF x max|gamma(t)| in percent, for
  ISTR 0 and ISTR 1, with a near-field column softer than the excavated soil (true SSI);
* the HOUSE property update against the plain linear path: the iteration-1 matrices COOSK/COOSM must
  equal those of a linear run whose near-field elements carry M-table materials with
  ``Vs sqrt(G/Gmax)``, ``Vp sqrt(G/Gmax)`` (constant Poisson's ratio) and ``beta_p = beta_s = D``,
  with G/Gmax and D interpolated by hand (linear in log10 strain); the .pin material lines map to the
  group's materials in ascending number (spec 05b OQ22), not in element order;
* requirements 2.5 for the iterations: a New Structure restart (ANALYS <mode> 1) with *changed*
  near-field properties equals a fresh initiation run with the same FILE4;
* the 2D branch (PLANE groups, HOUSE <dim> = 1) through HOUSE and STRESS with a synthetic FILE8 of a
  uniform strain field (ISTR 1 / ISTR 0 = sqrt((ex - ez)^2 + gxz^2) / |gxz|), since 2D ANALYS is not
  available;
* UT-20 by hand and the COMB_XYZ_STRAIN update with per-element G_max.

Tests marked ``xfail(strict=True)`` document defects found by the review (the reason names the
defect); they turn into XPASS failures once the defect is fixed, so the marker must then be removed.
"""
from __future__ import annotations

import math
from itertools import permutations
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pytest

from sassi.core import nlsoil as NL
from sassi.core import shake as SH
from sassi.io import decks
from sassi.io.container import read_container
from sassi.io.thfile import write_history
from sassi.prep import Interpreter, registry
from sassi.prep.messages import Kind
from sassi.verify import builders as B
from sassi.verify.problems.vp_motion import write_file8

GRAV = 9.81
RHO, NU, BETA0 = 2.0, 1.0 / 3.0, 0.02
DT, NFT = 0.02, 256
DF = 1.0 / (DT * NFT)                                   # 0.1953125 Hz

#: a hand-made curve (strain %, G/Gmax, damping %): log-linear interpolation is easy to do by hand
CURVE = dict(strain=[1e-4, 1e-2, 1e-1, 1.0], g=[1.0, 0.8, 0.5, 0.2], d=[0.5, 3.0, 8.0, 15.0])


def _vp(vs: float) -> float:
    return vs * math.sqrt(2.0 * (1.0 - NU) / (1.0 - 2.0 * NU))


def _curve() -> SH.DynamicProperty:
    s = np.array(CURVE["strain"])
    return SH.DynamicProperty("Hand", s, np.array(CURVE["g"]), s, np.array(CURVE["d"]))


def _by_hand(strain_pct: float, ys: Sequence[float]) -> float:
    """Linear interpolation in log10(strain) between the bracketing points of CURVE (no numpy.interp)."""
    xs = CURVE["strain"]
    if strain_pct <= xs[0]:
        return ys[0]
    for k in range(len(xs) - 1):
        if xs[k] <= strain_pct <= xs[k + 1]:
            t = (math.log10(strain_pct) - math.log10(xs[k])) / (math.log10(xs[k + 1]) - math.log10(xs[k]))
            return ys[k] + t * (ys[k + 1] - ys[k])
    return ys[-1]


def _site(vs_layers=(200.0, 200.0)) -> B.Site:
    return B.layered_site([(1.0, vs, _vp(vs), RHO, BETA0) for vs in vs_layers], (1000.0, 2000.0, 2.2, 0.01),
                          gravity=GRAV, nl=10)


def _column(site: B.Site, nl_mats: Sequence[Tuple[float, float]]) -> Tuple[B.HouseBuilder, int, Dict[int, int]]:
    """FV excavation of the top len(nl_mats) TOPL layers, 1 m x 1 m in plan: excavated SOLIDs (ETYPE 2, MSET =
    layer) and, at the same nodes, near-field SOLIDs (ETYPE 1) of type-3 materials ``(vs, beta)`` (level 0 =
    top).  The materials are created bottom-up, so material numbers do not follow the element numbers.
    Returns (builder, near-field group, {element: material id})."""
    hb = B.HouseBuilder(site, imp=0, incomp=1, title="review column")
    z = site.elevations()
    n = len(nl_mats)
    ids = {}
    for lev in range(n, -1, -1):                                      # bottom-up node numbering
        for j, y in enumerate((-0.5, 0.5)):
            for i, x in enumerate((-0.5, 0.5)):
                ids[(i, j, lev)] = hb.node(x, y, z[lev])
    mids = {}
    for lev in range(n - 1, -1, -1):                                  # materials bottom-up
        vs, beta = nl_mats[lev]
        mids[lev] = hb.material(3, _vp(vs), vs, RHO * GRAV, beta, beta)
    gexc = hb.group(1, "excavated soil")
    gnl = hb.group(1, "near-field soil")
    emat = {}
    for lev in range(n):
        bot = [ids[(0, 0, lev + 1)], ids[(1, 0, lev + 1)], ids[(1, 1, lev + 1)], ids[(0, 1, lev + 1)]]
        top = [ids[(0, 0, lev)], ids[(1, 0, lev)], ids[(1, 1, lev)], ids[(0, 1, lev)]]
        hb.solid(bot + top, mat=lev + 1, etype=2, group=gexc)
        e = hb.solid(bot + top, mat=mids[lev], etype=1, group=gnl)
        emat[e] = mids[lev]
    hb.set_interaction(sorted(ids.values()))
    return hb, gnl, emat


def _stress_deck(thfile: str) -> decks.Deck:
    d = decks.new("STRESS")
    d.params.update(dict(model="m", title="review: non-linear soil strains", iter=1, save=0, thfile=thfile,
                         mult=1.0, max=0.0, gravity=GRAV, delt=DT, nft=NFT, df=DF, interopt=1))
    return d


def _motion(n: int = 180, pga: float = 0.3, seed: int = 7) -> np.ndarray:
    """Deterministic band-limited record (g), not the implementer's generator."""
    rng = np.random.default_rng(seed)
    f = np.fft.rfftfreq(n, DT)
    a = np.fft.irfft(np.fft.rfft(rng.standard_normal(n)) * ((f > 0.4) & (f < 10.0)), n)
    t = np.arange(n) * DT
    a *= np.minimum(1.0, t / 0.5) * np.exp(-np.maximum(0.0, t - 2.0))
    return a * (pga / np.max(np.abs(a)))


# =============================================================================================
# 1. UT-20 by hand; COMB_XYZ_STRAIN with per-element G_max
# =============================================================================================
def test_ut20_by_hand_and_order_independence(tmp_path):
    assert float(NL.srss(0.10, 0.08, 0.02)) == pytest.approx(math.sqrt(0.01 + 0.0064 + 0.0004), rel=1e-15)
    assert round(float(NL.srss(0.10, 0.08, 0.02)), 4) == 0.1296                          # UT-20
    gmax = {1: 8.0e4, 2: 1.25e5}
    files = []
    for lab, gam in (("X", (0.10, 0.004)), ("Y", (0.08, 0.02)), ("Z", (0.02, 0.0))):
        nl = NL.NLFile("FILE74", 3, 0.65, 1, direction=lab, groups=[NL.NLGroupInfo(5, "SOLID", 2, 0)])
        nl.rows = [NL.NLRow(5, e, gam[e - 1], 0.0, 0.0, 1, gmax[e], gam[e - 1] / 0.65) for e in (1, 2)]
        files.append(nl)
    curves = {1: _curve()}
    results = []
    for order in permutations(range(3)):
        comb, warns = NL.combine_files([files[k] for k in order], curves)
        assert not warns
        results.append([(r.element, r.gamma_eff, r.G, r.beta) for r in comb.rows])
    assert all(r == results[0] for r in results)                                       # SRSS is symmetric
    comb, _ = NL.combine_files(files, curves)
    for r in comb.rows:
        g = math.sqrt(sum(f.rows[r.element - 1].gamma_eff ** 2 for f in files))
        assert r.gamma_eff == pytest.approx(g, rel=1e-14)
        assert r.G == pytest.approx(gmax[r.element] * _by_hand(g, CURVE["g"]), rel=1e-12)     # own G_max
        assert r.beta == pytest.approx(_by_hand(g, CURVE["d"]) / 100.0, rel=1e-12)


# =============================================================================================
# 2. HOUSE property update == the linear path with equivalent M-table materials
# =============================================================================================
def test_house_update_equals_linear_model_with_equivalent_materials(tmp_path):
    """Iteration 0: the two .pin material lines go to the materials in ascending number (here the lower
    element gets the higher material number).  Iteration 1: G = G_max (G/Gmax)(gamma), beta_s = beta_p = D,
    Poisson's ratio and density kept -- COOSK and COOSM must equal a linear HOUSE run of the same model with
    type-3 materials (Vs sqrt(r), Vp sqrt(r), D, D)."""
    wd = tmp_path
    site = _site()
    nl_mats = [(250.0, 0.03), (200.0, 0.02)]                         # level 0 (element 1), level 1 (element 2)
    hb, gnl, emat = _column(site, nl_mats)
    assert emat[1] > emat[2]                                         # element 1 has the higher material number
    B.write_deck(wd, "m", B.site_deck(site, B.FrequencySet.fourier(DT, NFT, [1, 5]), model="m"))
    d = hb.deck("m")
    d["nlssi"] = 1
    B.write_deck(wd, "m", d)
    SH.write_file73(wd / "FILE73", {"Hand": _curve()})
    lines = [NL.PinMaterial(0.5, 2.0, 1), NL.PinMaterial(0.8, 1.5, 1)]      # ascending material numbers
    NL.write_pin(wd / "m.pin", NL.PinData(0.6, 1, [NL.PinGroup(gnl, 2, 2, 0, lines)]))
    assert B.run("HOUSE", wd, "m") == 0
    f78 = {r.element: r for r in NL.read_nlfile(wd / "FILE78", "FILE78").rows}
    order = sorted(emat.values())
    for e, (vs, beta) in zip((1, 2), nl_mats):
        line = lines[order.index(emat[e])]
        # requirements 4.4 item 7 (normative): iteration 0 scales the free-field TOPL layer at the element
        # (both layers of _site(): Vs 200 m/s, BETA0), not the element material; G_max is the material.
        # (Amended by the implementer: the original lines checked GFAC G_mat / DFAC beta_mat, the rule the
        # review itself reported as a spec deviation -- test_iteration0_scales_the_free_field_layer_properties.)
        assert f78[e].G == pytest.approx(line.gfac * RHO * 200.0 ** 2, rel=1e-12)
        assert f78[e].beta == pytest.approx(line.dfac * BETA0, rel=1e-12)
        assert f78[e].gmax == pytest.approx(RHO * vs ** 2, rel=1e-12)
    # iteration 1 from strains between curve points
    strains = {1: 0.003, 2: 0.05}
    NL.write_nlfile(wd / "FILE74", NL.NLFile("FILE74", 0, 0.6, 1, direction="X", rows=[
        NL.NLRow(gnl, e, g, 0.0, 0.0, 1, RHO * nl_mats[e - 1][0] ** 2, g / 0.6) for e, g in strains.items()]))
    assert B.run("HOUSE", wd, "m") == 0
    Ks1 = read_container(wd / "COOSK").sparse("Ks").toarray()
    Ms1 = read_container(wd / "COOSM").sparse("Ms").toarray()
    f78 = {r.element: r for r in NL.read_nlfile(wd / "FILE78", "FILE78").rows}
    # the same model, linear, with equivalent materials for the near-field elements
    hb2 = B.HouseBuilder(site, imp=0, incomp=1)
    hb2.nodes, hb2.fixity = dict(hb.nodes), dict(hb.fixity)
    hb2.groups, hb2.interaction = [list(g) for g in hb.groups], list(hb.interaction)
    hb2.materials, hb2._mat_cache = [list(m) for m in hb.materials], dict(hb._mat_cache)
    for row in hb.elements:
        row = list(row)
        if row[0] == gnl:
            e = row[1]
            vs, _ = nl_mats[e - 1]
            r = _by_hand(strains[e], CURVE["g"])
            dmp = _by_hand(strains[e], CURVE["d"]) / 100.0
            assert f78[e].G == pytest.approx(r * RHO * vs ** 2, rel=1e-12)
            assert f78[e].beta == pytest.approx(dmp, rel=1e-12)
            row[3] = hb2.material(3, _vp(vs) * math.sqrt(r), vs * math.sqrt(r), RHO * GRAV, dmp, dmp)
        hb2.elements.append(row)
    d2 = hb2.deck("m")
    d2["nlssi"] = 0
    B.write_deck(wd, "m", d2)
    assert B.run("HOUSE", wd, "m") == 0
    Ks0 = read_container(wd / "COOSK").sparse("Ks").toarray()
    Ms0 = read_container(wd / "COOSM").sparse("Ms").toarray()
    assert np.max(np.abs(Ks1 - Ks0)) <= 1e-12 * np.max(np.abs(Ks0))
    assert np.max(np.abs(Ms1 - Ms0)) <= 1e-12 * np.max(np.abs(Ms0))


# =============================================================================================
# 3. FILE74 strains of STRESS == an independent recomputation from FILE8
# =============================================================================================
def _chain(wd: Path, istr: int, gfac: float):
    """SITE + POINT on all Fourier bins up to 11.9 Hz (interpolation exact on the grid), HOUSE (iteration 0,
    near-field column GFAC softer than the excavated soil), ANALYS, STRESS <iter> = 1."""
    site = _site()
    fs = B.FrequencySet.fourier(DT, NFT, range(1, int(12.0 / DF) + 1))
    hb, gnl, _ = _column(site, [(200.0, BETA0), (200.0, BETA0)])
    B.run_soil(wd, "m", site, fs, layer=2, rad=0.9)
    d = hb.deck("m")
    d["nlssi"] = 1
    B.write_deck(wd, "m", d)
    SH.write_file73(wd / "FILE73", {"Hand": _curve()})
    NL.write_pin(wd / "m.pin", NL.PinData(0.6, 1, [NL.PinGroup(gnl, 1, 2, istr, [NL.PinMaterial(gfac, 1.0, 1)])]))
    acc = _motion()
    write_history(wd / "ctrl.acc", acc, DT)
    B.write_deck(wd, "m", _stress_deck("ctrl.acc"))
    B.run("HOUSE", wd, "m")
    B.run_analys(wd, "m", fs, gravity=GRAV, save=0)
    B.run("STRESS", wd, "m")
    return hb, gnl, acc


@pytest.mark.parametrize("istr", [0, 1])
def test_file74_strains_match_an_independent_recomputation_from_file8(tmp_path, istr):
    hb, gnl, acc = _chain(tmp_path, istr, gfac=0.5)
    f8 = read_container(tmp_path / "FILE8", "FILE8")
    H, fnum = np.asarray(f8["H"]), np.asarray(f8["fnum"])
    eqn, eqd = np.asarray(f8["eq_node"]), np.asarray(f8["eq_dof"])

    def tf(n, k):
        hit = np.flatnonzero((eqn == n) & (eqd == k))
        return H[:, hit[0]] if hit.size else np.zeros(len(fnum), complex)

    padded = np.zeros(NFT)
    padded[:len(acc)] = acc
    w = 2.0 * np.pi * np.arange(NFT // 2 + 1) * DF
    Ug = np.zeros(NFT // 2 + 1, complex)
    Ug[1:] = -GRAV * np.fft.rfft(padded)[1:] / w[1:] ** 2                # ground displacement spectrum
    f74 = {r.element: r for r in NL.read_nlfile(tmp_path / "FILE74", "FILE74").rows}
    assert f74 and all(r.group == gnl for r in f74.values())
    for row in (e for e in hb.elements if e[0] == gnl):
        nodes = row[7:15]
        X = np.array([hb.nodes[n] for n in nodes])
        c, half = X.mean(axis=0), (X.max(axis=0) - X.min(axis=0)) / 2.0
        s = np.sign(X - c)                                              # corner natural coordinates (+-1)
        U = np.stack([np.stack([tf(n, k) for k in (1, 2, 3)], axis=1) for n in nodes], axis=1)   # (nF, 8, 3)
        grad = np.einsum("aj,fai->fij", s, U) / (8.0 * half[None, None, :])     # du_i/dx_j at the centroid
        e6 = np.stack([grad[:, 0, 0], grad[:, 1, 1], grad[:, 2, 2], grad[:, 0, 1] + grad[:, 1, 0],
                       grad[:, 0, 2] + grad[:, 2, 0], grad[:, 1, 2] + grad[:, 2, 1]], axis=1)
        grid = np.zeros((NFT // 2 + 1, 6), complex)
        grid[fnum] = e6
        h = np.fft.irfft(grid * Ug[:, None], n=NFT, axis=0)
        if istr == 0:
            gam = np.max(np.abs(h[:, 3:]), axis=1)
        else:
            ex, ey, ez, gxy, gxz, gyz = (h[:, k] for k in range(6))
            gam = (2.0 / 3.0) * np.sqrt((ex - ey) ** 2 + (ey - ez) ** 2 + (ez - ex) ** 2
                                        + 1.5 * (gxy ** 2 + gxz ** 2 + gyz ** 2))
        gmax_pct = 100.0 * float(np.max(gam))
        r = f74[row[1]]
        assert r.extra == pytest.approx(gmax_pct, rel=1e-8)
        assert r.gamma_eff == pytest.approx(0.6 * gmax_pct, rel=1e-8)
        assert r.G == pytest.approx(r.gmax * _by_hand(r.gamma_eff, CURVE["g"]), rel=1e-9)
        assert r.beta == pytest.approx(_by_hand(r.gamma_eff, CURVE["d"]) / 100.0, rel=1e-9)


# =============================================================================================
# 4. New Structure restart with changed near-field properties == fresh initiation (requirements 2.5)
# =============================================================================================
def test_new_structure_restart_with_changed_properties_equals_a_fresh_initiation(tmp_path):
    wd = tmp_path
    site = _site()
    fs = B.FrequencySet.fourier(DT, NFT, [1, 4, 9, 17, 30, 45])
    hb, gnl, _ = _column(site, [(200.0, BETA0), (200.0, BETA0)])
    B.run_soil(wd, "m", site, fs, layer=2, rad=0.9)
    d = hb.deck("m")
    d["nlssi"] = 1
    B.write_deck(wd, "m", d)
    SH.write_file73(wd / "FILE73", {"Hand": _curve()})
    NL.write_pin(wd / "m.pin", NL.PinData(0.6, 1, [NL.PinGroup(gnl, 1, 2, 0, [NL.PinMaterial(1.0, 1.0, 1)])]))
    B.run("HOUSE", wd, "m")
    B.run_analys(wd, "m", fs, gravity=GRAV, mode=0, save=1)
    H0 = np.asarray(read_container(wd / "FILE8", "FILE8")["H"])
    gm = RHO * 200.0 ** 2
    NL.write_nlfile(wd / "FILE74", NL.NLFile("FILE74", 0, 0.6, 1, direction="X", rows=[
        NL.NLRow(gnl, 1, 0.03, 0, 0, 1, gm, 0.05), NL.NLRow(gnl, 2, 0.3, 0, 0, 1, gm, 0.5)]))
    B.run("HOUSE", wd, "m")                                            # iteration 1: properties change
    B.run_analys(wd, "m", fs, gravity=GRAV, mode=1, save=0)            # New Structure restart (COOX reused)
    c1 = read_container(wd / "FILE8", "FILE8")
    H1 = np.asarray(c1["H"])
    assert int(c1.meta.get("mode", -1)) == 1
    for p in wd.glob("COO[XT]*"):
        p.unlink()
    B.run_analys(wd, "m", fs, gravity=GRAV, mode=0, save=0)            # fresh initiation, same FILE4
    H2 = np.asarray(read_container(wd / "FILE8", "FILE8")["H"])
    assert np.max(np.abs(H1 - H0)) > 1e-3 * np.max(np.abs(H0))         # the iteration did change FILE8
    assert np.max(np.abs(H1 - H2)) <= 1e-10 * np.max(np.abs(H2))


# =============================================================================================
# 5. 2D: a PLANE nonlinear group through HOUSE (<dim> = 1) and STRESS (synthetic FILE8)
# =============================================================================================
def test_plane_group_strain_measures_through_house_and_stress(tmp_path):
    wd = tmp_path
    site = _site()
    fs = B.FrequencySet.fourier(DT, NFT, range(1, 40))
    B.write_deck(wd, "m", B.site_deck(site, fs, model="m"))
    hb = B.HouseBuilder(site, dim=1, imp=0, incomp=1)
    z = site.elevations()
    ids = {(i, lev): hb.node(x, 0.0, z[lev]) for lev in (2, 1, 0) for i, x in enumerate((-0.5, 0.5))}
    gexc, gnl = hb.group(4, "excavated"), hb.group(4, "near field")
    mid = hb.material(3, _vp(200.0), 200.0, RHO * GRAV, BETA0, BETA0)
    for lev in range(2):
        quad = [ids[(0, lev + 1)], ids[(1, lev + 1)], ids[(1, lev)], ids[(0, lev)]]
        hb._element(gexc, quad, etype=2, mat=lev + 1)
        hb._element(gnl, quad, etype=1, mat=mid)
    hb.set_interaction(sorted(ids.values()))
    d = hb.deck("m")
    d["nlssi"] = 1
    B.write_deck(wd, "m", d)
    SH.write_file73(wd / "FILE73", {"Hand": _curve()})
    NL.write_pin(wd / "m.pin", NL.PinData(0.6, 1, [NL.PinGroup(gnl, 1, 2, 0, [NL.PinMaterial(1.0, 1.0, 1)])]))
    assert B.run("HOUSE", wd, "m") == 0
    f78 = NL.read_nlfile(wd / "FILE78", "FILE78")
    assert [(g.group, g.etype) for g in f78.groups] == [(gnl, "PLANE")]
    f4 = read_container(wd / "m.N4", "FILE4")
    eqn, eqd = np.asarray(f4["eq_node"]), np.asarray(f4["eq_dof"])
    assert set(eqd.tolist()) == {1, 3}
    # uniform strain field: ux = e1 x + a z, uz = b x + e2 z  ->  exx = e1, ezz = e2, gxz = a + b
    e1, e2, a, b = 0.3, -0.1, 0.7, 0.2
    vals = np.array([(e1 * hb.nodes[n][0] + a * hb.nodes[n][2]) if k == 1 else (b * hb.nodes[n][0] + e2 * hb.nodes[n][2])
                     for n, k in zip(eqn, eqd)], dtype=complex)
    write_file8(wd / "FILE8", fs.fnum, DF, eqn, eqd, np.tile(vals, (len(fs.fnum), 1)), nfft=NFT, delt=DT,
                model_hash=str(f4.meta["model_hash"]))
    write_history(wd / "ctrl.acc", _motion(150), DT)
    B.write_deck(wd, "m", _stress_deck("ctrl.acc"))
    out = {}
    for istr in (0, 1):
        f78.groups[0].istr = istr
        NL.write_nlfile(wd / "FILE78", f78)
        assert B.run("STRESS", wd, "m") == 0
        out[istr] = NL.read_nlfile(wd / "FILE74", "FILE74")
    for r0, r1 in zip(out[0].rows, out[1].rows):
        assert r1.extra / r0.extra == pytest.approx(math.hypot(e1 - e2, a + b) / abs(a + b), rel=1e-9)
    assert out[0].rows[0].extra == pytest.approx(out[0].rows[1].extra, rel=1e-9)        # uniform field


# =============================================================================================
# 6. Defects found by the review (fixed in the package -- their strict xfail markers were removed; the
#    HOUSEX warning lives in sassi/prep/commands/extensions.py, outside the package, and stays a strict xfail)
# =============================================================================================
@pytest.fixture
def fake_iteration():
    """REVIEWNLSTEP,<k>: writes the FILE78/FILE74 of iteration k - 1 of a converging sequence (flat curve
    G/Gmax 0.5): the change in G is 30 % x 0.5^(k-1)."""
    calls: List[int] = []

    @registry.command("REVIEWNLSTEP", tier="P0")
    def _step(c):
        k = c.int(1)
        calls.append(k)
        d = Path(c.model.path)
        it = k - 1
        G = 100.0 * (0.5 + 0.3 * 0.5 ** it)
        NL.write_nlfile(d / "FILE78", NL.NLFile("FILE78", it, 0.65, 1, rows=[NL.NLRow(3, 1, 0.01, G, 0.05, 1, 100.0, 1)]))
        NL.write_nlfile(d / "FILE74", NL.NLFile("FILE74", it, 0.65, 1, rows=[NL.NLRow(3, 1, 0.02, 0, 0, 1, 100.0, 0.03)]))

    yield calls
    registry.REGISTRY.pop("REVIEWNLSTEP", None)


def test_nlssiiter_after_nlssireset_runs_the_new_analysis(tmp_path, fake_iteration):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text("MDL,m,m\n")
    d = tmp_path / "m"
    d.mkdir(exist_ok=True)
    SH.write_file73(d / "FILE73", {"Flat": SH.DynamicProperty("Flat", np.array([1e-4, 10.0]), np.array([0.5, 0.5]),
                                                             np.array([1e-4, 10.0]), np.array([5.0, 5.0]))})
    # end state of a previous, converged analysis (iteration 6)
    NL.write_nlfile(d / "FILE78", NL.NLFile("FILE78", 6, 0.65, 1, rows=[NL.NLRow(3, 1, 0.01, 50.1, 0.05, 1, 100.0, 1)]))
    NL.write_nlfile(d / "FILE74", NL.NLFile("FILE74", 6, 0.65, 1, rows=[NL.NLRow(3, 1, 0.02, 0, 0, 1, 100.0, 0.03)]))
    NL.write_liq(d / "m.liq", 1)
    assert ui.execute("NLSSIRESET") and NL.read_liq(d / "m.liq") == 0
    assert ui.execute('VAR,BODY,"REVIEWNLSTEP,#"')
    assert ui.execute("NLSSIITER,BODY,8")
    assert fake_iteration and fake_iteration[0] == 1                   # the new analysis must run


MODEL_ETYPE0 = """MDL,m,m
L,1,1.0,19.62,400,200,0.02,0.02
L,2,1.0,19.62,400,200,0.02,0.02
L,3,1.0,21.58,2000,1000,0.01,0.01
TOPL,1,2
FREQ,1,4,20
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.01,1024,1
WAVE,2,1,1,1,0
HOUSE,9.81,0,0,2,0,0,0,0,0
DYNP,1,0.0001,1.0,0.0001,0.5,Soft
DYNP,2,1.0,0.2,1.0,15.0,Soft
N,1,0,0,-2
N,2,1,0,-2
N,3,1,1,-2
N,4,0,1,-2
NGEN,2,4,1,4,1,0,0,1
GROUP,1,SOLID
MACT,1
E,1,5,6,7,8,9,10,11,12
ETYPE,1,1,1,2
GROUP,2,SOLID
MACT,2
E,1,1,2,3,4,5,6,7,8
ETYPE,1,1,1,2
M,1,400,200,19.62,0.02,0.02,3
GROUP,3,SOLID
MACT,1
E,1,1,2,3,4,5,6,7,8
E,2,5,6,7,8,9,10,11,12
INT,1,12,1,1
POINT,0,2,0.9
PIN,0.65
PINGRP,3,0,1.0,1.0,Soft
HOUSEX,0,0,1,1
STRESS,0,1,0,0,1
AOPT,0,0,0,1,0,1,0,0,0,0,0,0,0,0
"""


def test_check_reports_etype0_elements_of_a_nonlinear_group(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(MODEL_ETYPE0)
    assert not ui.sink.texts(Kind.ERROR)
    ui.execute("CHECK")
    err = (tmp_path / "m" / "m.err").read_text()
    assert "EDU-41" in err


def test_iteration0_scales_the_free_field_layer_properties(tmp_path):
    """The free field (TOPL layer 1, L table) carries strain-compatible soil, Vs 150 m/s, 4 % damping; the
    near-field material holds the low-strain soil (Vs 200 m/s, 2 %).  GFAC 0.9, DFAC 1.5."""
    wd = tmp_path
    site = _site((150.0, 150.0))
    site.layers = [B.SoilLayer(l.thick, l.vs, l.vp, l.rho, 0.04) for l in site.layers]
    hb, gnl, _ = _column(site, [(200.0, BETA0)])
    B.write_deck(wd, "m", B.site_deck(site, B.FrequencySet.fourier(DT, NFT, [1, 5]), model="m"))
    d = hb.deck("m")
    d["nlssi"] = 1
    B.write_deck(wd, "m", d)
    NL.write_pin(wd / "m.pin", NL.PinData(0.65, 1, [NL.PinGroup(gnl, 1, 1, 0, [NL.PinMaterial(0.9, 1.5, 1)])]))
    assert B.run("HOUSE", wd, "m") == 0
    r = NL.read_nlfile(wd / "FILE78", "FILE78").rows[0]
    assert r.G == pytest.approx(0.9 * RHO * 150.0 ** 2, rel=1e-9)
    assert r.beta == pytest.approx(1.5 * 0.04, rel=1e-9)


def _model_etype1(**replace) -> str:
    text = MODEL_ETYPE0.replace("E,2,5,6,7,8,9,10,11,12\n", "E,2,5,6,7,8,9,10,11,12\nETYPE,1,2,1,1\n")
    for old, new in replace.items():
        text = text.replace(old, new)
    return text


def test_check_keeps_error_79_when_no_nonlinear_house_run_feeds_stress(tmp_path):
    (tmp_path / "m").mkdir()
    (tmp_path / "m" / "acc.th").write_text("0.01\n0.0\n0.1\n0.0\n")
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(_model_etype1(**{"HOUSEX,0,0,1,1": "HOUSEX,0,0,1,0",
                                 "AOPT,0,0,0,1,0,1,0,0,0,0,0,0,0,0": "AOPT,0,0,0,1,0,1,0,0,0,0,0,1,0,0\nTHFILE,acc.th"}))
    ui.execute("CHECK")
    assert "Error 79" in (tmp_path / "m" / "m.err").read_text()


def test_check_validates_the_curves_referenced_by_the_pin(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(_model_etype1(**{"PINGRP,3,0,1.0,1.0,Soft": "DYNP,1,0.0001,1.0,,,Back\nDYNP,2,1.0,0.3,,,Back\n"
                                                           "PINGRP,3,0,1.0,1.0,Back"}))
    assert not ui.sink.texts(Kind.ERROR)
    ui.execute("CHECK")
    err = (tmp_path / "m" / "m.err").read_text()
    assert any("Back" in ln for ln in err.splitlines())


def test_stress_refuses_a_file8_of_another_iteration(tmp_path):
    _chain(tmp_path, istr=0, gfac=1.0)                                 # iteration 0: HOUSE, ANALYS, STRESS
    assert NL.read_nlfile(tmp_path / "FILE74", "FILE74").iteration == 0
    B.run("HOUSE", tmp_path, "m")                                      # iteration 1 -- ANALYS not run
    assert NL.read_nlfile(tmp_path / "FILE78", "FILE78").iteration == 1
    rc = B.run("STRESS", tmp_path, "m", check=False)
    assert rc != 0 or NL.read_nlfile(tmp_path / "FILE74", "FILE74").iteration != 1


def test_housex_nlssi_is_not_reported_as_unavailable(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text("MDL,m,m\n")
    assert ui.execute("HOUSEX,0,0,1,1")
    assert not any("not available" in t.lower() and "non-linear ssi" in t.lower() for t in ui.sink.texts(Kind.WARNING))
