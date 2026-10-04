"""Unit tests of the 2D SSI helpers (sassi/core/ssi2d.py) and of the flexibility blocks used by the SYMM
image sums (sassi/core/flexibility.flexibility_block, sassi/core/strip2d.flexibility_block_2d;
requirements 4.6 item 2, D-PNT-01, D-W2-07)."""
from __future__ import annotations


import numpy as np
import pytest

from sassi.core import ssi2d as S2
from sassi.core import strip2d
from sassi.core.flexibility import flexibility_block, flexibility_matrix
from sassi.io.container import read_container
from sassi.verify import builders as B

SITE = [(2.0, 150.0, 300.0, 1.9, 0.05), (3.0, 250.0, 500.0, 2.0, 0.04)]
HS = (400.0, 800.0, 2.1, 0.03)


def test_dimension_and_input_rules():
    assert S2.point_dimension_error(1, 1) == "" and S2.point_dimension_error(2, 2) == ""
    assert "re-run POINT with <dim> = 1" in S2.point_dimension_error(1, 2)
    assert "POINT2 (2D line-load)" in S2.point_dimension_error(2, 1)
    assert S2.inplane_input_error(0, 0.0) == "" and S2.inplane_input_error(2, 180.0) == ""
    assert "anti-plane" in S2.inplane_input_error(1, 0.0)
    assert "0 deg" in S2.inplane_input_error(0, 30.0) and S2.inplane_input_error(0, 540.0) == ""


def test_plane_coordinates_and_duplicates():
    xyz = np.array([[1.0, 5.0, -2.0], [1.0, -3.0, -2.0], [2.0, 0.0, 0.0]])
    np.testing.assert_array_equal(S2.plane_coordinates(xyz, 0.5)[:, 1], [0.5, 0.5, 0.5])
    assert np.unique(S2.duplicate_key(xyz, 1), axis=0).shape[0] == 2            # y ignored in 2D
    assert np.unique(S2.duplicate_key(xyz, 2), axis=0).shape[0] == 3


def test_rigid_transform_2d_is_a_rigid_body_motion():
    xyz = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 0.0, -1.0]])
    T = S2.rigid_transform_2d(xyz, (1.0, 0.0))
    u = (T @ np.array([0.0, 0.0, 0.01])).reshape(-1, 3)          # rotation theta_y = 0.01 about (1, 0)
    np.testing.assert_allclose(u[:, 0], 0.01 * (xyz[:, 2] - 0.0))
    np.testing.assert_allclose(u[:, 2], -0.01 * (xyz[:, 0] - 1.0))
    assert np.all(u[:, 1] == 0.0)


def test_strip_reference_values():
    """R2 B.4 (Jakub & Roesset 1977) rows for nu = 0.30 and 0.45."""
    kx, kp = S2.strip_stiffness_reference(1.0, 4.0, 2.0, 0.30)
    assert kx == pytest.approx(2.0 * 1.175 * (1 + 2.15 / 4)) and kp == pytest.approx(2.0 * 2.394 * (1 + 0.17 / 4))
    kx, kp = S2.strip_stiffness_reference(2.0, 4.0, 1.0, 0.45)
    assert kx == pytest.approx(1.419 * 1.975) and kp == pytest.approx(4.0 * 2.907 * 1.12)
    with pytest.raises(ValueError):
        S2.strip_stiffness_reference(1.0, 4.0, 1.0, 0.25)


@pytest.fixture(scope="module")
def point_files(tmp_path_factory):
    out = {}
    fs = B.FrequencySet.harmonic(2.0, [2, 5])
    site = B.layered_site(SITE, HS)
    for dim in (2, 1):
        wd = tmp_path_factory.mktemp(f"flexblock{dim}")
        B.run_soil(wd, "m", site, fs, layer=1, rad=0.9, mode2=False, dim=dim)
        out[dim] = read_container(wd / "FILE3", "FILE3")
    return out


@pytest.mark.parametrize("dim", [2, 1])
def test_flexibility_block_equals_the_matrix(point_files, dim):
    """With observation = load set, the block is the unsymmetrised F_ff bit for bit; blocks are
    consistent (obs, load) <-> the full matrix of the union of both sets."""
    f3 = point_files[dim]
    xy = np.array([[0.0, 0.0], [1.3, 0.4], [2.6, -0.7], [0.0, 1.1], [1.3, 0.4]])
    iface = [1, 1, 2, 2, 2]
    if dim == 1:
        xy[:, 1] = 0.0
        xy[3, 0] = 4.1
    for q in (0, 1):
        F = flexibility_matrix(f3, q, xy, iface, symmetrize=False)
        Fb = flexibility_block(f3, q, xy, iface, xy, iface)
        assert np.array_equal(F, Fb)
        Fo = flexibility_block(f3, q, xy[:2], iface[:2], xy[2:], iface[2:])
        np.testing.assert_allclose(Fo, F[:6, 6:], rtol=1e-9, atol=1e-12 * np.abs(F).max())


def test_flexibility_block_2d_parity(point_files):
    """P-SV couplings are odd in x - x_load: mirrored load points flip their sign (2D image terms)."""
    f3 = point_files[1]
    xo = np.array([1.0, 2.5])
    A = strip2d.flexibility_block_2d(f3, 0, xo, [0, 1], np.array([0.0]), [1]).reshape(2, 3, 1, 3)
    Bm = strip2d.flexibility_block_2d(f3, 0, -xo, [0, 1], np.array([0.0]), [1]).reshape(2, 3, 1, 3)
    np.testing.assert_allclose(Bm[:, 2, 0, 0], -A[:, 2, 0, 0], rtol=1e-14)
    np.testing.assert_allclose(Bm[:, 0, 0, 0], A[:, 0, 0, 0], rtol=1e-14)
    np.testing.assert_allclose(Bm[:, 0, 0, 2], -A[:, 0, 0, 2], rtol=1e-14)
