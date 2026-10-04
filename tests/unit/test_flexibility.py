"""Unit tests of the F_ff assembly from FILE3 (sassi.core.flexibility, D-ANL-02)."""
from __future__ import annotations

import time

import numpy as np
import pytest

from sassi.core import flexibility as fl
from sassi.io import decks
from sassi.io.files import read_container
from sassi.modules.base import run_module

G = 9.81


@pytest.fixture(scope="module")
def file3(tmp_path_factory):
    wd = tmp_path_factory.mktemp("flex")
    d = decks.new("SITE")
    d["gravity"] = G
    d["nl"] = 20
    d["df"] = 2.0
    d["mode2"] = 0
    for i, (h, vs) in enumerate([(0.5, 150.0)] * 6 + [(1.0, 250.0)] * 6):
        d.table("layers").append([i + 1, h, 1.9 * G, 2 * vs, vs, 0.05, 0.05])
    d.table("halfspace").append([13, 0.0, 2.1 * G, 800.0, 400.0, 0.02, 0.02])
    for n in (1, 3):
        d.table("freqs").append([n])
    decks.write(wd / "m.sit", d)
    assert run_module("SITE", "m", wd) == 0
    p = decks.new("POINT")
    p["layer"] = 3
    p["rad"] = 0.9
    p["df"] = 2.0
    decks.write(wd / "m.poi", p)
    assert run_module("POINT", "m", wd) == 0
    return read_container(wd / "FILE3", "FILE3")


def _rot(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def test_rotation_invariance(file3):
    """Rotating the whole node set by an angle rotates every 3x3 block: F' = R F R^T."""
    xy = np.array([[0, 0], [1.1, 0.3], [-2.0, 3.5], [4.0, -1.0], [0.4, 0.2]])
    ifc = [1, 2, 4, 1, 3]
    F = fl.flexibility_matrix(file3, 1, xy, ifc, symmetrize=False)
    a = np.radians(37.0)
    R = _rot(a)
    F2 = fl.flexibility_matrix(file3, 1, xy @ R[:2, :2].T, ifc, symmetrize=False)
    n = len(xy)
    Rb = np.kron(np.eye(n), R)
    with np.errstate(all="ignore"):                 # spurious Accelerate FPE flags in matmul
        ref = Rb @ F @ Rb.T
    np.testing.assert_allclose(F2, ref, atol=1e-12 * np.abs(F).max())


def test_same_vertical_line_uses_axis_values(file3):
    pd = fl.point_data(file3, 0)
    F = fl.flexibility_matrix(file3, 0, np.array([[2.0, 3.0], [2.0, 3.0]]), [1, 3], symmetrize=False)
    blk = F[3:6, 0:3]                      # node on interface 3 due to loads at the node on interface 1
    assert blk[0, 0] == pd.axis1[0, 2, 0] and blk[1, 1] == pd.axis1[0, 2, 0]
    assert blk[2, 2] == pd.axis0[0, 2, 2]
    assert blk[0, 1] == 0 and blk[0, 2] == 0 and blk[2, 0] == 0
    self_blk = F[0:3, 0:3]
    assert self_blk[0, 0] == pd.axis1[0, 0, 0] and self_blk[2, 2] == pd.axis0[0, 0, 2]


def test_inner_zone_is_linear_between_axis_and_rim(file3):
    pd = fl.point_data(file3, 0)
    r = np.array([0.0, 0.3, 0.6, 0.9])
    c = fl.cylindrical_components(pd, r)
    for key in fl.COMPONENTS:
        lin = c[key][0] + (r / 0.9)[:, None, None] * (c[key][3] - c[key][0])
        np.testing.assert_allclose(c[key], lin, rtol=1e-12, atol=1e-20)


def test_tabulated_far_field_matches_exact(file3):
    pd = fl.point_data(file3, 1)
    r = np.linspace(0.9, 60.0, 400)
    exact = fl.cylindrical_components(pd, r)
    table = fl.cylindrical_components(pd, r, max_exact=10)
    for key in fl.COMPONENTS:
        scale = np.abs(exact[key]).max(axis=0, keepdims=True)
        assert np.abs(table[key] - exact[key]).max() / scale.max() < 1e-5


def test_symmetrisation_and_blocking(file3):
    rng = np.random.default_rng(3)
    xy = rng.uniform(-6, 6, size=(12, 2))
    ifc = rng.integers(1, 5, size=12)
    F = fl.flexibility_matrix(file3, 0, xy, ifc, symmetrize=False)
    Fs = fl.flexibility_matrix(file3, 0, xy, ifc)
    np.testing.assert_allclose(Fs, 0.5 * (F + F.T), rtol=0, atol=1e-15 * np.abs(F).max())
    np.testing.assert_array_equal(Fs, Fs.T)
    Fb = fl.flexibility_matrix(file3, 0, xy, ifc, block=1)
    np.testing.assert_allclose(Fb, Fs, rtol=1e-13)
    # random pairs include r < R0 and r ~ R0 on coarse layers, where the FE-core reciprocity is only
    # approximate (VP-07 checks the 1e-3 criterion of requirements §6.3 on its model)
    assert np.linalg.norm(F - F.T) / np.linalg.norm(F) < 1e-2


def test_errors_and_frequency_row(file3):
    with pytest.raises(ValueError, match="POINT <layer>"):
        fl.flexibility_matrix(file3, 0, np.zeros((2, 2)) + [[0, 0], [1, 0]], [1, 6])
    with pytest.raises(IndexError):
        fl.flexibility_matrix(file3, 5, np.zeros((1, 2)), [1])
    assert fl.frequency_row(file3, 3) == 1
    with pytest.raises(KeyError):
        fl.frequency_row(file3, 2)


def test_large_grid_performance(file3):
    """D-GEN-10 scale check: a 30 x 30 surface grid (900 nodes) assembles in a few seconds."""
    g = np.arange(30) * 1.0
    X, Y = np.meshgrid(g, g)
    xy = np.c_[X.ravel(), Y.ravel()]
    t = time.time()
    F = fl.flexibility_matrix(file3, 0, xy, np.ones(len(xy), int))
    assert F.shape == (2700, 2700)
    assert time.time() - t < 30.0


def test_file3_carries_column_loss(file3):
    """POINT stores the column's material loss (x_loss) so that ANALYS classifies propagating modes
    with the same damping-widened sector as SITE (D-SIT-09)."""
    beta = 0.05                      # largest damping ratio of the fixture site (layers 5 %, half-space 2 %)
    np.testing.assert_allclose(file3["x_loss"], beta / np.sqrt(1 - beta ** 2), rtol=1e-12)
    assert fl.point_data(file3, 0).loss == pytest.approx(beta / np.sqrt(1 - beta ** 2))


def test_table_grid_resolves_heavily_damped_soil():
    """beta = 0.46 (EDU-04 admits beta < 0.5): every wave decays by |Im k|/Re k = 0.52 from material
    damping alone.  The tabulation grid must still resolve the shortest propagating wavelength
    (lambda_min/48) instead of falling back to the coarse rmax/48 spacing."""
    from sassi.core import tlm
    n, beta = 40, 0.46
    G, M = tlm.layer_moduli(2.0, 100.0, 200.0, beta, beta)
    col = tlm.Column(h=[0.25] * n, rho=[2.0] * n, G=[complex(G)] * n, M=[complex(M)] * n)
    md = tlm.column_modes(col, 2 * np.pi * 20.0)
    lam_min = 2 * np.pi / md.kR[tlm.select_mode(md.kR, 1, tlm.column_loss(col))].real
    for loss in (tlm.column_loss(col), None):
        pd = fl.PointData(R0=0.9, load_iface=np.array([1]), kR=md.kR, kL=md.kL, phix=None, phiz=None, phiy=None,
                          alpha0=None, alpha1=None, axis0=None, axis1=None, loss=loss)
        grid = fl._table_grid(pd, 60.0)
        assert grid[0] == 0.9 and grid[-1] == 60.0
        assert np.diff(grid).max() <= lam_min / 48 * (1 + 1e-12)
