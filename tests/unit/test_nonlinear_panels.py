"""Unit tests of the Option NON panel helpers (sassi/core/panels.py): corner detection, panel kinematics,
SHEAR equations and the BBCGEN curve (requirements 4.15, spec 10 section 3.11, D-NON-09, D-NON-10)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.core import panels as P


def wall_mesh(nx=4, nz=2, L=12.0, H=4.0, angle=0.0, x0=(0.0, 0.0, 0.0)):
    """Quad mesh of a vertical wall rotated by ``angle`` (degrees) about Z: (elements, xyz)."""
    e = np.array([math.cos(math.radians(angle)), math.sin(math.radians(angle)), 0.0])
    xyz = {}
    nid = lambda i, k: k * (nx + 1) + i + 1
    for k in range(nz + 1):
        for i in range(nx + 1):
            xyz[nid(i, k)] = np.asarray(x0) + e * (L * i / nx) + np.array([0, 0, H * k / nz])
    els = [(nid(i, k), nid(i + 1, k), nid(i + 1, k + 1), nid(i, k + 1)) for k in range(nz) for i in range(nx)]
    return els, xyz


def test_corner_counting_and_ordering():
    els, xyz = wall_mesh()
    assert P.count_corners(els) == [1, 5, 11, 15]
    g = P.panel_geometry(els, xyz)
    assert g.corners == (1, 5, 15, 11) and g.corner_method == "count"
    assert g.length == pytest.approx(12.0) and g.height == pytest.approx(4.0)
    assert (g.width_extent, g.height_extent) == (pytest.approx(12.0), pytest.approx(4.0))
    np.testing.assert_allclose(g.e_h, [1, 0, 0], atol=1e-12)
    assert g.vertical and not g.warnings


def test_oblique_wall_axes():
    els, xyz = wall_mesh(angle=150.0, x0=(5.0, -3.0, 1.0))
    g = P.panel_geometry(els, xyz)
    # e_h oriented toward +X: 'left' is the smaller x
    assert g.e_h[0] > 0 and g.length == pytest.approx(12.0) and g.height == pytest.approx(4.0)
    assert xyz[g.corners[0]][0] < xyz[g.corners[1]][0]


def test_horizontal_plane_rejected_and_triangles_fall_back():
    els, xyz = wall_mesh()
    flat = {k: np.array([v[0], v[2], 0.0]) for k, v in xyz.items()}
    with pytest.raises(P.PanelError, match="not vertical"):
        P.panel_geometry(els, flat)
    tri = []
    for a, b, c, d in els:
        tri += [(a, b, c), (a, c, d)]
    g = P.panel_geometry(tri, xyz)
    assert g.corner_method == "geometry" and g.corners == (1, 5, 15, 11)
    assert any("counting" in w for w in g.warnings)


def _corner_disp(g, xyz, field):
    out = {}
    for lab, n in zip(P.CORNERS, g.corners):
        u = field(np.asarray(xyz[n]))
        out[lab] = (np.array([u[0]]), np.array([u[1]]), np.array([u[2]]))
    return out


def test_panel_deformations_rigid_body_invariant():
    els, xyz = wall_mesh(angle=30.0)
    g = P.panel_geometry(els, xyz)
    th = 1e-3
    c = np.mean(list(xyz.values()), axis=0)
    n = g.normal
    for field in (lambda p: np.array([0.3, -0.2, 0.7]),                         # translation
                  lambda p: np.cross(th * n, p - c),                           # in-plane rigid rotation
                  lambda p: np.cross(np.array([0.0, 0.0, th]), p - c)):        # rotation about Z
        gam, kap, eps = P.panel_deformations(_corner_disp(g, xyz, field), g.e_h, g.length, g.height)
        assert abs(gam[0]) < 1e-14 and abs(kap[0]) < 1e-14 and abs(eps[0]) < 1e-14


def test_panel_deformations_pure_states():
    els, xyz = wall_mesh()
    g = P.panel_geometry(els, xyz)
    z0 = 0.0
    shear = lambda p: np.array([0.002 * (p[2] - z0), 0.0, 0.0])              # u = gamma z
    gam, kap, eps = P.panel_deformations(_corner_disp(g, xyz, shear), g.e_h, g.length, g.height)
    assert gam[0] == pytest.approx(0.002) and abs(kap[0]) < 1e-15 and abs(eps[0]) < 1e-15
    axial = lambda p: np.array([0.0, 0.0, -0.001 * p[2]])
    gam, kap, eps = P.panel_deformations(_corner_disp(g, xyz, axial), g.e_h, g.length, g.height)
    assert eps[0] == pytest.approx(-0.001) and abs(gam[0]) < 1e-15
    # bending: w = kappa x z (the top edge rotates relative to the bottom edge)
    bend = lambda p: np.array([-0.5e-4 * p[2] ** 2, 0.0, 1e-4 * (p[0] - 6.0) * p[2]])
    gam, kap, eps = P.panel_deformations(_corner_disp(g, xyz, bend), g.e_h, g.length, g.height)
    assert kap[0] == pytest.approx(1e-4)


def test_unit_system_detection():
    assert P.unit_system(32.2)[0] == "british" and P.unit_system(9.81)[0] == "si"
    assert P.unit_system(30.0)[2] and not P.unit_system(32.174)[2]


def test_shear_equations_british_reference():
    """spec 10 section 3.11 worked example (VP-46 inputs): h_W 20 ft, l_W 30 ft, t_W 2 ft, f'c 5 ksi,
    fy 60 ksi, rho 0.005, N_U 1000 kips."""
    c = P.shear_capacities(20, 30, 2, 5, 60, 0.005, 1000, gravity=32.2)
    assert c.units == "british" and c.force_unit == "kips" and c.alpha_c == 3.0
    for got, ref in ((c.aci, 4424.82), (c.upper, 6109.40), (c.wood_raw, 648.00), (c.wood_lower, 3665.64),
                     (c.barda, 4026.77), (c.gw, 2405.90)):
        assert got == pytest.approx(ref, rel=1e-6)
    assert c.wood == c.wood_lower                                     # Wood clamped to its lower bound
    c2 = P.shear_capacities(20, 30, 2, 5, 60, 0.005, 1000, a_be=0.1, fy_be=60, gravity=32.2)
    assert c2.gw == pytest.approx(2617.54, rel=1e-6)
    # forces given directly (EDUOPT,SHEARFORCEARGS,1): F_VW = 2592 kips, F_BE = 864 kips
    c3 = P.shear_capacities(20, 30, 2, 5, 60, 0.005, 1000, a_be=2592.0, fy_be=864.0, gravity=32.2, force_args=True)
    assert c3.gw == pytest.approx(c2.gw, rel=1e-12)


def test_alpha_c_interpolation_and_aci_cap():
    assert P.alpha_c(1.0) == 3.0 and P.alpha_c(2.5) == 2.0 and P.alpha_c(1.75) == pytest.approx(2.5)
    c = P.shear_capacities(20, 30, 2, 5, 60, 0.05, 0, gravity=32.2)        # heavy reinforcement: ACI capped
    assert c.aci_raw > c.upper and c.aci == c.upper
    assert c.wood == c.upper                                               # Wood clamped from above


def test_ultimate_selection():
    c = P.shear_capacities(20, 30, 2, 5, 60, 0.005, 1000, gravity=32.2)
    assert [c.ultimate(k) for k in (1, 2, 3, 4)] == [c.aci, c.wood, c.barda, c.gw]
    with pytest.raises(P.PanelError):
        c.ultimate(5)


def test_bbcgen_curve():
    vu, vcr, ga = 2405.9, 1832.82, 1.2e7
    g, v, yi = P.bbcgen_curve(vu, vcr, ga)
    assert len(g) == 22 and yi == 21
    assert g[0] == pytest.approx(vcr / ga) and v[0] == vcr
    assert (g[20], v[20]) == (pytest.approx(0.004), pytest.approx(vu))
    assert (g[21], v[21]) == (0.02, pytest.approx(1.02 * vu))
    np.testing.assert_allclose(np.diff(g[:21]), (0.004 - vcr / ga) / 20)   # 20 equal strain steps
    assert np.all(np.diff(v) > 0)
    xi = (g[1:21] - g[0]) / (0.004 - g[0])
    np.testing.assert_allclose(v[1:21], vcr + (vu - vcr) * (1 - (1 - xi) ** 2))
    with pytest.raises(P.PanelError):
        P.bbcgen_curve(vu, vu * 1.1, ga)
    with pytest.raises(P.PanelError):
        P.bbcgen_curve(vu, vcr, 1.0)                                       # gamma_cr above gamma_y


def test_plane_clusters_and_components():
    n = np.array([[0, 0, 1.0], [0, 0, 1.0], [0, 1.0, 0], [0, 0, 1.0]])
    d = np.array([0.0, 0.0, 2.0, 3.0])
    lab = P.plane_clusters(n, d, offset_tol=1e-9)
    assert list(lab) == [0, 0, 1, 2]
    comp = P.connected_components(5, [(0, 1), (3, 4)])
    assert list(comp) == [0, 0, 2, 3, 3]
