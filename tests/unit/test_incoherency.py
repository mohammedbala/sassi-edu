"""Unit tests of the incoherency core (sassi.core.incoherency) and of the incoherent free-field load
(sassi.core.ssi_solver.incoherent_seismic_load): requirements 4.4 item 6, 4.6 item 5, D-INC-02/03/06/
07/09/10; spec 05b section 11 tests 1-7."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.core import coherency as C
from sassi.core import incoherency as I
from sassi.core import ssi_solver as ss
from sassi.io.container import Container

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")


def _grid(n: int, h: float) -> np.ndarray:
    g = np.arange(n) * h
    X, Y = np.meshgrid(g, g, indexing="xy")
    return np.column_stack([X.ravel(), Y.ravel()])


def _lw(xy, f, gamma=0.3, vs=300.0):
    return C.CoherencySpec(model=1, gamma=(gamma,) * 3, alpha=vs).matrix_builder(xy, 0)(f)


# --------------------------------------------------------------------------------------
# decomposition (Eqs. 6.1-6.4, D-INC-03)
# --------------------------------------------------------------------------------------
def test_decomposition_trace_order_and_reconstruction():
    xy = _grid(4, 3.0)
    S = _lw(xy, 8.0)
    dec = I.decompose(S)
    n = xy.shape[0]
    assert dec.trace_error < 1e-12 and np.all(np.diff(dec.lam) <= 0) and np.all(dec.lam >= 0)
    assert np.allclose(dec.phi @ np.diag(dec.lam) @ dec.phi.T, S, atol=1e-12)
    assert np.allclose(dec.phi.T @ dec.phi, np.eye(n), atol=1e-12)
    ups = dec.contributions()
    assert ups.sum() == pytest.approx(100.0) and ups[0] == pytest.approx(100 * dec.lam[0] / n)
    m90 = dec.modes_for(0.9)
    cum = np.cumsum(dec.lam) / n
    assert cum[m90 - 1] >= 0.9 - 1e-12 and (m90 == 1 or cum[m90 - 2] < 0.9)


def test_sign_convention_and_ties():
    xy = _grid(3, 5.0)                                  # symmetric layout: zero-sum (antisymmetric) modes
    dec = I.decompose(_lw(xy, 6.0))
    sums = dec.phi.sum(axis=0)
    scale = np.abs(dec.phi).sum(axis=0)
    zero = np.abs(sums) <= 1e-9 * scale
    assert np.all(sums[~zero] > 0)
    assert zero.any()
    for k in np.flatnonzero(zero):                      # first significant component positive
        col = dec.phi[:, k]
        first = np.flatnonzero(np.abs(col) > 1e-6 * np.abs(col).max())[0]
        assert col[first] > 0
    raw = I.decompose(_lw(xy, 6.0), sign="RAW")
    assert np.allclose(np.abs(raw.lam), np.abs(dec.lam))


def test_canonical_basis_is_independent_of_the_eigensolver_basis():
    rng = np.random.default_rng(3)
    Q, _ = np.linalg.qr(rng.standard_normal((12, 3)))
    R, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    A = I.canonical_basis(Q)
    B = I.canonical_basis(Q @ R)
    assert np.allclose(A, B, atol=1e-12)
    assert np.allclose(A.T @ A, np.eye(3), atol=1e-12)
    assert np.allclose(A @ A.T, Q @ Q.T, atol=1e-12)                 # same subspace
    one = np.ones(12)
    assert np.allclose(A[:, 0], Q @ Q.T @ one / np.linalg.norm(Q.T @ one))   # first: projection of 1


def test_as_factors_do_not_depend_on_round_off_in_degenerate_eigenspaces():
    """A symmetric grid has equal eigenvalues; a perturbation of the order of round-off must not change
    the AS factors (canonical basis), as it would with the raw eigensolver basis."""
    xy = _grid(5, 3.0)
    S = _lw(xy, 10.0)
    assert any(j - i > 1 for i, j in I.eigen_clusters(I.decompose(S, "RAW").lam))
    E = np.random.default_rng(5).standard_normal(S.shape) * 1e-15
    S2 = S + 0.5 * (E + E.T)
    a = I.synthesize(I.decompose(S), np.arange(25))
    b = I.synthesize(I.decompose(S2), np.arange(25))
    assert np.abs(a - b).max() < 1e-9


def test_zero_frequency_limit_gives_coherent_motion():
    xy = _grid(4, 4.0)
    dec = I.decompose(np.ones((16, 16)))
    assert dec.lam[0] == pytest.approx(16.0) and np.allclose(dec.lam[1:], 0.0, atol=1e-12)
    assert np.allclose(I.synthesize(dec, np.arange(16)), 1.0, atol=1e-12)
    assert np.allclose(I.synthesize(dec, np.array([0])), 1.0, atol=1e-12)        # mode 1 alone
    assert np.allclose(I.synthesize(dec, np.array([3])), 0.0, atol=1e-6)         # a higher mode vanishes
    s = I.synthesize(I.decompose(_lw(xy, 1e-3)), np.arange(16))
    assert np.abs(s - 1.0).max() < 1e-4


def test_select_modes():
    assert I.select_modes(0, 5).tolist() == [0, 1, 2, 3, 4]
    assert I.select_modes(2, 5).tolist() == [0, 1]
    assert I.select_modes(9, 5).tolist() == [0, 1, 2, 3, 4]
    assert I.select_modes(-3, 5).tolist() == [2]
    with pytest.raises(I.IncoherencyError, match="mode 6 does not exist"):
        I.select_modes(-6, 5)


def test_synthesis_forms():
    dec = I.decompose(_lw(_grid(3, 4.0), 12.0))
    modes = np.arange(9)
    a = dec.phi * np.sqrt(dec.lam)[None, :]
    assert np.allclose(I.synthesize(dec, modes), a.sum(axis=1))                  # AS
    assert np.allclose(I.synthesize(dec, np.array([1])), a[:, 1])                # single mode
    th = np.linspace(-1, 1, 9)
    assert np.allclose(I.synthesize(dec, modes, th), a @ np.exp(1j * th))        # random phases


# --------------------------------------------------------------------------------------
# random phases (D-INC-07)
# --------------------------------------------------------------------------------------
def test_phase_generators_are_reproducible_and_independent():
    g1 = I.phase_generators(11, 17, 3)
    g2 = I.phase_generators(11, 17, 3)
    a = [[g.uniform(size=4) for g in row] for row in g1]
    b = [[g.uniform(size=4) for g in row] for row in g2]
    assert all(np.array_equal(x, y) for ra, rb in zip(a, b) for x, y in zip(ra, rb))
    flat = np.array([x for row in a for x in row])
    assert len({tuple(v) for v in flat}) == 9                                    # all streams differ
    c = I.phase_generators(11, 99, 3)
    assert np.array_equal(c[0][0].uniform(size=4), a[0][0]) and not np.array_equal(c[2][0].uniform(size=4), a[2][0])
    assert not np.array_equal(I.phase_generators(-11, 17, 1)[0][0].uniform(size=4), a[0][0])


def _run(xy, freq, **kw):
    st = I.IncoherencySettings(**kw)
    spec = C.CoherencySpec(model=st.cohf, gamma=st.gamma, alpha=st.alpha) if st.coh else None
    return I.FactorRun(st, spec, np.column_stack([xy, np.zeros(len(xy))]), np.arange(1, len(freq) + 1),
                       np.asarray(freq, float))


def test_stochastic_samples_do_not_depend_on_chunking():
    run = _run(_grid(3, 5.0), [2.0, 6.0, 9.0], coh=1, cohf=1, alpha=250.0, gamma=(0.3, 0.3, 0.3), hseed=4,
               vseed=8, randphz=180.0, nsim=5)
    all5, _, _ = I.compute_factors(run, list(range(5)))
    tail, _, _ = I.compute_factors(run, [3, 4])
    assert np.array_equal(all5[3:], tail)
    assert not np.allclose(all5[0], all5[1])
    # deterministic input (zero seeds or zero phase): the AS factors
    det, _, _ = I.compute_factors(_run(_grid(3, 5.0), [2.0, 6.0, 9.0], coh=1, cohf=1, alpha=250.0,
                                       gamma=(0.3, 0.3, 0.3), hseed=4, vseed=8, randphz=0.0), [0])
    dec = I.decompose(_lw(_grid(3, 5.0), 6.0, 0.3, 250.0))
    assert np.allclose(det[0, 1, 0], I.synthesize(dec, np.arange(9)))


def test_ensemble_mean_converges_to_the_coherency_matrix():
    xy = _grid(3, 5.0)
    run = _run(xy, [8.0], coh=1, cohf=1, alpha=250.0, gamma=(0.3, 0.3, 0.3), hseed=1, vseed=2, randphz=180.0,
               nsim=50)
    S, _, _ = I.compute_factors(run, list(range(50)))
    Sig = _lw(xy, 8.0, 0.3, 250.0)
    est = np.einsum("ri,rj->ij", S[:, 0, 0], np.conj(S[:, 0, 0])) / 50
    assert np.abs(est - Sig).max() < 0.5 and np.allclose(np.diag(est).real.mean(), 1.0, atol=0.2)


# --------------------------------------------------------------------------------------
# wave passage, multiple excitation, merging of plan positions
# --------------------------------------------------------------------------------------
def test_wave_passage_delays_and_phase_difference():
    xy = np.array([[2.0, 1.0], [2.0 + 10 * math.cos(0.5), 1.0 + 10 * math.sin(0.5)]])
    tau = I.arrival_delays(xy, 500.0, math.degrees(0.5), xc=2.0, yc=1.0)
    assert tau[0] == pytest.approx(0.0, abs=1e-15) and tau[1] == pytest.approx(10.0 / 500.0)
    w = I.wave_passage_factors([4.0], tau)
    assert np.angle(w[0, 1] / w[0, 0]) == pytest.approx(-2 * math.pi * 4.0 * 10.0 / 500.0)
    with pytest.raises(I.IncoherencyError, match="Error 113"):
        I.arrival_delays(xy, 0.0, 0.0)
    run = _run(xy, [4.0], coh=0, wpass=1, appv=500.0, wang=math.degrees(0.5), xc=2.0, yc=1.0)
    S, stats, _ = I.compute_factors(run, [0])
    assert stats[0][0] is None
    assert np.allclose(S[0, 0, 0], w[0]) and np.allclose(S[0, 0, 2], w[0])       # every direction


def test_multiple_excitation_zones_and_errors():
    ids = [1, 2, 3, 4, 5, 6]
    rows = [{"no": 1, "nfirst": 1, "nlast": 3}, {"no": 2, "nfirst": 5, "nlast": 6}]
    amp = [{"no": 1, "idx": k, "re": 1.5, "im": 0.0} for k in (1, 2)] + \
          [{"no": 2, "idx": k, "re": 0.0, "im": 2.0} for k in (1, 2)]
    zones, errs, warns = I.me_zones(rows, amp, ids, 2, cmplxspec=1)
    assert not errs and [z.nodes.tolist() for z in zones] == [[0, 1, 2], [4, 5]]
    sar = I.sar_factors(zones, 2, 6)
    assert np.allclose(sar[:, :3], 1.5) and np.allclose(sar[:, 3], 1.0) and np.allclose(sar[:, 4:], 2.0j)
    # real ratios ignore the imaginary parts (warning)
    zones, errs, warns = I.me_zones(rows, amp, ids, 2, cmplxspec=0)
    assert not errs and any("imaginary" in w for w in warns)
    assert I.me_zones([], [], ids, 2, 0)[1][0].startswith("Error 115")
    assert I.me_zones([{"no": 1, "nfirst": 9, "nlast": 3}], amp, ids, 2, 0)[1][0].startswith("Error 116")
    assert I.me_zones([{"no": 1, "nfirst": 1, "nlast": 9}], amp, ids, 2, 0)[1][0].startswith("Error 117")
    bad = [{"no": 1, "idx": k, "re": 11.0, "im": 0.0} for k in (1, 2)]
    assert I.me_zones(rows[:1], bad, ids, 2, 0)[1][0].startswith("Error 118")
    assert I.me_zones(rows[:1], amp[:1], ids, 2, 0)[1][0].startswith("Error 119")
    over = [{"no": 1, "nfirst": 1, "nlast": 4}, {"no": 2, "nfirst": 3, "nlast": 6}]
    assert any("overlaps" in e for e in I.me_zones(over, amp, ids, 2, 1)[1])


def test_merged_plan_positions_share_one_factor():
    xy = np.array([[0.0, 0.0], [5.0, 0.0], [0.0, 0.0], [5.0, 5.0]])     # nodes 1 and 3 above each other
    run = _run(xy, [10.0], coh=1, cohf=1, alpha=200.0, gamma=(0.3, 0.3, 0.3), merge=1)
    S, stats, info = I.compute_factors(run, [0])
    assert info["n_positions"] == 3 and info["merged"] and stats[0][0].n == 3
    assert S[0, 0, 0, 0] == S[0, 0, 0, 2]
    pos, idx = I.unique_positions(xy, 1e-9)
    assert pos.shape == (3, 2) and idx.tolist() == [0, 1, 0, 2]


def test_settings_validation():
    ok = I.IncoherencySettings(coh=1, cohf=1, alpha=300.0)
    assert ok.validate(9, 2, False, 0)[0] == [] and ok.method == "AS" and not ok.stochastic
    errs = I.IncoherencySettings(coh=1, cohf=5).validate(9, 1, True, 0)[0]
    assert any("3D model" in e for e in errs) and any("SYMM" in e for e in errs) \
        and any("wave passage" in e for e in errs)
    sto = I.IncoherencySettings(coh=1, cohf=1, alpha=300.0, hseed=1, vseed=2, randphz=180.0, nsim=60, nmodes=-2,
                                supmode=1)
    assert sto.stochastic and sto.method == "SS"
    e = " | ".join(sto.validate(9, 2, False, 0)[0])
    assert "1..50" in e and "one mode" in e and "Quadratic" in e
    quad = I.IncoherencySettings(coh=1, cohf=1, alpha=300.0, supmode=1, nmodes=0)
    assert any("one HOUSE + ANALYS run per mode" in x for x in quad.validate(9, 2, False, 0)[0])
    assert I.IncoherencySettings(coh=1, cohf=1, alpha=300.0, nmodes=-3).method == "SINGLE"
    assert any("Error 60" in x for x in I.IncoherencySettings(coh=1, cohf=1, alpha=300.0, ngp=0)
               .validate(9, 2, False, 2)[0])
    me = I.IncoherencySettings(me=1, wpass=0).validate(9, 2, False, 0)[0]
    assert any("multiple excitation needs wave passage" in x for x in me)
    assert any("Error 113" in x for x in I.IncoherencySettings(wpass=1, appv=0.0).validate(9, 2, False, 0)[0])
    assert I.IncoherencySettings(wpass=1).method == "COHERENT"


def test_option_file_roundtrip(tmp_path):
    assert I.read_options(tmp_path) == {"INCOHSIGN": "ADJUST", "INCOHMERGE": "0"}
    (tmp_path / I.OPTION_FILE).write_text(I.options_text("raw", "1"))
    assert I.read_options(tmp_path) == {"INCOHSIGN": "RAW", "INCOHMERGE": "1"}


# --------------------------------------------------------------------------------------
# FILE77 and Build_FILE77
# --------------------------------------------------------------------------------------
def _write77(path, ids, xyz, fnum=(3, 5), value=1.0, **meta):
    m = dict(df=0.5, method="AS", stochastic=0)
    m.update(meta)
    s = np.full((len(fnum), 3, len(ids)), value, complex) * (1 + np.arange(len(ids)))[None, None, :]
    return I.write_file77(path, list(fnum), np.asarray(fnum) * 0.5, ids, xyz, s, m)


def test_file77_roundtrip_and_node_permutation(tmp_path):
    xyz = np.array([[0, 0, 0], [1, 0, 0], [2, 0, -1.0]])
    _write77(tmp_path / "FILE77", [10, 20, 30], xyz)
    c = I.read_file77(tmp_path / "FILE77")
    assert np.asarray(c["s"]).shape == (2, 3, 3)
    perm = I.node_permutation(c, [30, 10], xyz[[2, 0]])
    assert perm.tolist() == [2, 0]
    with pytest.raises(I.IncoherencyError, match="no factors for interaction nodes \\[40\\]"):
        I.node_permutation(c, [10, 40])
    with pytest.raises(I.IncoherencyError, match="another position"):
        I.node_permutation(c, [10], np.array([[0.5, 0.0, 0.0]]))


def test_build_file77(tmp_path):
    _write77(tmp_path / "L1", [1, 2], np.array([[0, 0, -5.0], [1, 0, -5.0]]), value=1.0)
    _write77(tmp_path / "L2", [3, 4], np.array([[0, 0, 0.0], [1, 0, 0.0]]), value=10.0)
    a, b = I.read_file77(tmp_path / "L1"), I.read_file77(tmp_path / "L2")
    arrays, meta, notes = I.build_file77([a, b], ["L1", "L2"])
    assert arrays["int_node"].tolist() == [1, 2, 3, 4] and meta["method"] == "BUILT"
    assert np.allclose(arrays["s"][0, 0], [1, 2, 10, 20])
    arrays, _, notes = I.build_file77([a, b], ["L1", "L2"], order=[3, 1, 4])
    assert arrays["int_node"].tolist() == [3, 1, 4] and np.allclose(arrays["s"][0, 0], [10, 1, 20])
    assert any("dropped" in n for n in notes)
    with pytest.raises(I.IncoherencyError, match="already in an earlier input"):
        I.build_file77([a, a], ["L1", "L1b"])
    _write77(tmp_path / "L3", [7], np.zeros((1, 3)), fnum=(3, 6))
    with pytest.raises(I.IncoherencyError, match="frequency numbers differ"):
        I.build_file77([a, I.read_file77(tmp_path / "L3")], ["L1", "L3"])
    with pytest.raises(I.IncoherencyError, match="none of the input"):
        I.build_file77([a], ["L1"], order=[1, 9])


def test_inco_table_has_the_search_key():
    lines = []
    dec = I.decompose(_lw(_grid(3, 4.0), 12.0))
    I.inco_table(lines.append, 0, 7, 3.5, 2, dec, max_rows=4)
    text = "\n".join(lines)
    assert "I N C O" in text and "direction Z" in text and "5 further modes" in text


# --------------------------------------------------------------------------------------
# the incoherent free-field load (ANALYS, requirements 4.6 item 5)
# --------------------------------------------------------------------------------------
def test_incoherent_seismic_load_ffl_ffm():
    rng = np.random.default_rng(7)
    n = 4
    A = rng.standard_normal((3 * n, 3 * n)) + 1j * rng.standard_normal((3 * n, 3 * n))
    X = A + A.T
    Up = rng.standard_normal((3 * n, 2)) + 1j * rng.standard_normal((3 * n, 2))
    int_eq = np.arange(3 * n).reshape(n, 3)
    int_eq[1, 2] = -1                                     # a fixed interaction translation
    part = ss.Partition.from_int_eq(int_eq, 3 * n + 2)
    ones = np.ones((n, 2))
    base = ss.seismic_load(X, Up, part)
    assert np.array_equal(ss.incoherent_seismic_load(X, Up, ones, part), base)
    assert np.array_equal(ss.incoherent_seismic_load(X, Up, ones, part, ffm=True), base)
    s = rng.standard_normal((n, 2)) + 1j * rng.standard_normal((n, 2))
    s3 = np.repeat(s, 3, axis=0)
    ffl = ss.incoherent_seismic_load(X, Up, s, part)
    ffm = ss.incoherent_seismic_load(X, Up, s, part, ffm=True)
    assert np.allclose(ffl, (s3 * (X @ Up))[part.x_pos])
    assert np.allclose(ffm, (X @ (s3 * Up))[part.x_pos])
    # one factor set for every case
    one = ss.incoherent_seismic_load(X, Up, s[:, 0], part)
    assert np.allclose(one[:, 1], (np.repeat(s[:, 0], 3) * (X @ Up[:, 1]))[part.x_pos])
