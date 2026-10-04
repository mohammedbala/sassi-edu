"""Unit tests of the unlagged coherency models (sassi.core.coherency; requirements 4.4 item 6,
D-INC-01, D-INC-04, D-INC-05, D-INC-11).

The Abrahamson coefficients are data files (sassi/data/coherency/*.json) with their citations; the
"plot tests" below evaluate the model curves as the published figures show them (5 separations,
0-50 Hz) and check that they equal 1 at zero separation and zero frequency, decay with frequency and
separation, stay within [0, 1] and reproduce the values digitised from the published figures."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from sassi.core import coherency as C

pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")


# --------------------------------------------------------------------------------------
# coefficient files
# --------------------------------------------------------------------------------------
def test_every_model_has_a_cited_coefficient_file():
    for n, name in C.MODEL_FILES.items():
        info = json.loads((C.DATA_DIR / name).read_text())
        assert info["model"] == n
        assert "primary" in info["source"] and len(info["source"]["primary"]) > 40
        if info["verified"] and n != 1:
            for comp in ("horizontal", "vertical"):
                assert set(info["components"][comp]) == {"a1", "a2", "a3", "n1", "n2", "fc"}
            src = info["source"]
            assert "Equation" in src["equation"] and "Table" in src["tables"]
            assert info["published_points"]


@pytest.mark.parametrize("n", [2, 4])
def test_unconfirmed_models_refuse_to_run(n):
    assert not C.model_available(n)
    with pytest.raises(C.CoherencyError, match="coefficients not available"):
        C.load_model(n)
    spec = C.CoherencySpec(model=n, alpha=0.5)
    with pytest.raises(C.CoherencyError, match="coefficients not available"):
        spec.validate()
    assert "use instead" in C.unavailable_reason(n)


@pytest.mark.parametrize("n", [1, 3, 5, 6, 7])
def test_available_models(n):
    assert C.model_available(n)


# --------------------------------------------------------------------------------------
# published example values (digitised from the figures of EPRI 1012968 / 1015110)
# --------------------------------------------------------------------------------------
def _points():
    for n in (3, 5, 6):
        for p in C.model_info(n)["published_points"]:
            yield n, p


@pytest.mark.parametrize("n,p", list(_points()), ids=lambda v: str(v) if isinstance(v, int) else
                         f"{v['component'][0]}{v['xi_m']:g}m{v['f_hz']:g}Hz")
def test_published_values(n, p):
    g = float(C.load_model(n)(p["f_hz"], p["xi_m"], p["component"]))
    assert abs(g - p["gamma"]) <= p["tol"], (n, p, g)


def test_hand_computed_values_of_the_equations():
    """Independent evaluation of the published equations (EPRI 1015110 Eqs. 5-2, 6-1, 7-1) at one point."""
    t = math.tanh(0.4 * 10.0)
    L = math.log(11.0)
    # 2007 hard rock, horizontal, 10 m, 20 Hz
    n1 = 3.80 - 0.040 * L + 0.0105 * (L - 3.6) ** 2
    fc = 27.9 - 4.82 * L + 1.24 * (L - 3.6) ** 2
    g = (1 + (20 * t / fc) ** n1) ** -0.5 * (1 + (20 * t / 40.0) ** 16.4) ** -0.5
    assert C.load_model(5)(20.0, 10.0, "horizontal") == pytest.approx(g, rel=1e-12)
    # 2005 all sites, vertical, 10 m, 20 Hz
    fc = math.exp(2.43 - 0.025 * L - 0.048 * L ** 2)
    g = (1 + (20 * t / (3.15 * fc)) ** 4.95) ** -0.5 * (1 + (20 * t / (1.0 * fc)) ** 1.685) ** -0.5
    assert C.load_model(3)(20.0, 10.0, "vertical") == pytest.approx(g, rel=1e-12)
    # 2007 soil, horizontal, 10 m, 10 Hz
    fc = 14.3 - 2.35 * L
    a2 = 15.8 - 0.044 * 10.0
    g = (1 + (10 * t / fc) ** 3) ** -0.5 * (1 + (10 * t / a2) ** 15) ** -0.5
    assert C.load_model(6)(10.0, 10.0, "horizontal") == pytest.approx(g, rel=1e-12)


# --------------------------------------------------------------------------------------
# "plot tests": the curve families of the published figures
# --------------------------------------------------------------------------------------
FREQ = np.linspace(0.0, 50.0, 501)
SEPS = (10.0, 25.0, 50.0, 100.0, 150.0)


@pytest.mark.parametrize("n", [3, 5, 6])
@pytest.mark.parametrize("comp", ["horizontal", "vertical"])
def test_curves_decay_with_frequency_and_separation(n, comp):
    m = C.load_model(n)
    curves = np.array([m(FREQ, xi, comp) for xi in SEPS])          # (5, nf)
    assert np.all((curves >= 0.0) & (curves <= 1.0))
    assert np.allclose(curves[:, 0], 1.0)                          # f = 0
    assert np.all(np.diff(curves, axis=1) <= 1e-15)                # non-increasing in f
    assert np.all(curves[:, -1] < curves[:, 0])                    # and decaying
    # larger separation, smaller coherency.  The published hard-rock vertical model (EPRI 1015110 Table
    # 6-2) is the one exception, by at most 3e-5 below 1 Hz: its exponent n1(xi) grows with xi, so at
    # f << fc a larger separation is very slightly MORE coherent -- a property of the coefficients
    tiny = 1e-4 if (n, comp) == (5, "vertical") else 1e-12
    assert np.all(np.diff(curves, axis=0) <= tiny)
    assert np.all(np.diff(curves[:, FREQ >= 2.0], axis=0) <= 1e-12)
    assert np.allclose(m(FREQ, 0.0, comp), 1.0)                    # zero separation: 1 at every frequency
    xi = np.linspace(0.0, 150.0, 301)
    for f in (5.0, 10.0, 20.0, 40.0):
        g = m(f, xi, comp)
        assert np.all(np.diff(g) <= 1e-12) and g[0] == pytest.approx(1.0)


@pytest.mark.parametrize("n", [5, 6])
def test_separations_beyond_the_published_range_are_held(n):
    m = C.load_model(n)
    assert m.formula_max_m == 150.0
    for comp in ("horizontal", "vertical"):
        assert np.allclose(m(10.0, [150.0, 300.0, 1000.0], comp), m(10.0, 150.0, comp))


def test_figure_like_plot(tmp_path):
    """The curve families render (matplotlib Agg) -- the figures a student compares with the reports."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(9, 3))
    for ax, n in zip(axes, (3, 5, 6)):
        m = C.load_model(n)
        for xi in SEPS:
            ax.plot(FREQ, m(FREQ, xi, "horizontal"), label=f"{xi:g} m")
        ax.set_title(m.data["short"], fontsize=8)
    axes[0].legend(fontsize=6)
    p = tmp_path / "coherency.png"
    fig.savefig(p, dpi=40)
    plt.close(fig)
    assert p.stat().st_size > 1000


# --------------------------------------------------------------------------------------
# Luco-Wong
# --------------------------------------------------------------------------------------
def test_luco_wong():
    g = C.luco_wong(20.0, 10.0, 0.5, 1920.0)
    assert g == pytest.approx(math.exp(-(0.5 * 2 * math.pi * 20.0 * 10.0 / 1920.0) ** 2), rel=1e-14)
    assert C.luco_wong(30.0, 0.0, 0.3, 300.0) == 1.0
    assert C.luco_wong(0.0, 50.0, 0.3, 300.0) == 1.0
    f = np.linspace(0, 20, 41)
    assert np.all(np.diff(C.luco_wong(f, 10.0, 0.3, 300.0)) < 0)
    assert np.all(np.diff(C.luco_wong(5.0, np.linspace(0, 50, 21), 0.3, 300.0)) < 0)
    with pytest.raises(C.CoherencyError, match="Error 59"):
        C.luco_wong(1.0, 1.0, 0.3, 0.0)


# --------------------------------------------------------------------------------------
# distances
# --------------------------------------------------------------------------------------
def test_directional_distances_unit_separations():
    for ang in (0.0, 30.0, 135.0):
        a = math.radians(ang)
        ex = np.array([[0.0, 0.0], [math.cos(a), math.sin(a)]])        # along Line D
        ey = np.array([[0.0, 0.0], [-math.sin(a), math.cos(a)]])       # across Line D
        for alpha in (0.1, 0.5, 0.9):
            dx = C.directional_distances(ex, alpha, ang)[0, 1]
            dy = C.directional_distances(ey, alpha, ang)[0, 1]
            assert dx == pytest.approx(math.sqrt(2 * alpha), rel=1e-12)
            assert dy == pytest.approx(math.sqrt(2 * (1 - alpha)), rel=1e-12)
    assert (C.directional_distances(np.array([[0, 0], [1, 0.0]]), 0.1)[0, 1]
            / C.directional_distances(np.array([[0, 0], [0, 1.0]]), 0.1)[0, 1]) == pytest.approx(1 / 3)
    xy = np.random.default_rng(1).uniform(-10, 10, (6, 2))
    assert np.allclose(C.directional_distances(xy, 0.5, 37.0), C.horizontal_distances(xy))


# --------------------------------------------------------------------------------------
# user tables (model 7)
# --------------------------------------------------------------------------------------
def _write_tables(d, freq, dist, fun, comment=True):
    (d / "FREQCOH").write_text(("# frequencies\n" if comment else "") + " ".join(map(str, freq)))
    (d / "DISTCOH").write_text(", ".join(map(str, dist)))
    for k, nm in enumerate(C.USER_TABLES):
        (d / nm).write_text("\n".join(" ".join(f"{fun(f, x, k):.12g}" for x in dist) for f in freq))


def test_user_tables_bilinear_clamped(tmp_path):
    freq = [1.0, 5.0, 10.0]
    dist = [2.0, 10.0, 20.0]
    _write_tables(tmp_path, freq, dist, lambda f, x, k: 1.0 - 0.01 * f - 0.02 * x - 0.001 * k)
    u = C.UserCoherency.from_dir(tmp_path)
    assert u.freq.tolist() == freq and u.dist.tolist() == dist
    lin = lambda f, x, k: 1.0 - 0.01 * f - 0.02 * x - 0.001 * k
    assert u.evaluate(3.0, np.array([6.0]), 1)[0] == pytest.approx(lin(3.0, 6.0, 1), rel=1e-12)     # bilinear
    assert u.evaluate(50.0, np.array([15.0]), 0)[0] == pytest.approx(lin(10.0, 15.0, 0), rel=1e-12)  # clamped f
    assert u.evaluate(3.0, np.array([99.0]), 2)[0] == pytest.approx(lin(3.0, 20.0, 2), rel=1e-12)    # clamped D
    assert u.evaluate(3.0, np.array([0.0]), 0)[0] == 1.0                                             # D = 0
    # between D = 0 (anchor 1) and the first tabulated distance: linear
    g2 = lin(3.0, 2.0, 0)
    assert u.evaluate(3.0, np.array([1.0]), 0)[0] == pytest.approx(0.5 * (1.0 + g2), rel=1e-12)
    assert any("100 x 100" in n for n in u.notes)


def test_user_tables_errors(tmp_path):
    with pytest.raises(C.CoherencyError, match="needs the files"):
        C.UserCoherency.from_dir(tmp_path)
    _write_tables(tmp_path, [1.0, 2.0], [0.0, 5.0], lambda f, x, k: 1.5)
    with pytest.raises(C.CoherencyError, match=r"\[-1, 1\]"):
        C.UserCoherency.from_dir(tmp_path)
    _write_tables(tmp_path, [2.0, 1.0], [0.0, 5.0], lambda f, x, k: 0.5)
    with pytest.raises(C.CoherencyError, match="strictly increasing"):
        C.UserCoherency.from_dir(tmp_path)
    _write_tables(tmp_path, [1.0, 2.0], [0.0, 5.0], lambda f, x, k: 0.5)
    (tmp_path / "COHYUSER").write_text("0.5 0.5 0.5")
    with pytest.raises(C.CoherencyError, match="COHYUSER: 3 values"):
        C.UserCoherency.from_dir(tmp_path)


def test_user_tables_force_unit_coherency_at_zero_distance(tmp_path):
    _write_tables(tmp_path, [1.0, 2.0], [0.0, 5.0], lambda f, x, k: 0.8)
    u = C.UserCoherency.from_dir(tmp_path)
    assert np.all(u.tables[0][:, 0] == 1.0) and any("set to 1" in n for n in u.notes)


# --------------------------------------------------------------------------------------
# the specification of one HOUSE run
# --------------------------------------------------------------------------------------
def test_spec_validation_errors():
    with pytest.raises(C.CoherencyError, match="Error 58"):
        C.CoherencySpec(model=1, gamma=(0.05, 0.1, 0.2), alpha=300.0).validate()
    with pytest.raises(C.CoherencyError, match="Error 59"):
        C.CoherencySpec(model=1, alpha=0.0).validate()
    with pytest.raises(C.CoherencyError, match="Error 114"):
        C.CoherencySpec(model=5, alpha=1.5).validate()
    with pytest.raises(C.CoherencyError, match="Error 114"):
        C.CoherencySpec(model=8).validate()
    C.CoherencySpec(model=5, alpha=0.1).validate()


def test_matrix_builder_unit_diagonal_symmetric_and_units():
    xy = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 20.0], [15.0, 15.0]])
    spec = C.CoherencySpec(model=5, alpha=0.5, length_to_m=C.FT_TO_M)
    S = spec.matrix_builder(xy, 0)(20.0)
    assert np.allclose(np.diag(S), 1.0) and np.allclose(S, S.T)
    # distances in feet are converted to metres before the empirical model is evaluated
    assert S[0, 1] == pytest.approx(float(C.load_model(5)(20.0, 10.0 * C.FT_TO_M, "horizontal")), rel=1e-12)
    Z = spec.matrix_builder(xy, 2)(20.0)
    assert Z[0, 1] == pytest.approx(float(C.load_model(5)(20.0, 10.0 * C.FT_TO_M, "vertical")), rel=1e-12)
    assert spec.same_matrix(0, 1) and not spec.same_matrix(0, 2)
    lw = C.CoherencySpec(model=1, gamma=(0.2, 0.3, 0.3), alpha=300.0)
    assert not lw.same_matrix(0, 1) and lw.same_matrix(1, 2)
    S1 = lw.matrix_builder(xy, 1)(5.0)
    assert S1[0, 1] == pytest.approx(C.luco_wong(5.0, 10.0, 0.3, 300.0), rel=1e-14)
    assert C.length_to_metres(32.2) == C.FT_TO_M and C.length_to_metres(9.81) == 1.0
