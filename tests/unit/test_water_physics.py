"""Physics of the FILLPOOL water SOLIDs (docs/user/WATER.md "Why it works" and "Choice of the shear modulus").

A slice of water ``[-L, L] x [0, 1] x [0, H]`` of SOLID elements in a rigid tank with slip walls (the normal
displacement prescribed: ``u_x = U`` on x = +-L, ``u_y = 0`` on y = 0, 1, ``u_z = 0`` on the floor, the
surface free) is driven harmonically; the effective mass ``sum (M u)_x / U`` divided by the water mass must
equal the impulsive mass ratio of the exact potential-flow solution (:func:`sassi.prep.water_lib.
impulsive_ratio_exact`) when the water has the FILLPOOL properties (G = 1e-8 K) and incompatible modes.
The same slice shows why the alternatives fail: without incompatible modes the nearly incompressible
elements lock, and with nu = 0.49 (G = 0.02 K, the original D-WAT-01 proposal) the "water" is an elastic
solid that moves with the walls at low frequency.
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from sassi.conventions import cfactor
from sassi.elements.solid import solid_kernel
from sassi.prep import water_lib as wl


def effective_mass_ratio(L, H, nx, nz, G_over_K, freqs, incompatible=True, beta=0.005):
    wp = wl.water_properties(9.81)                       # kN, m, t: K = 2.2e6 kPa, rho = 1 t/m3
    K, rho = wp.K, wp.rho
    G = G_over_K * K
    xs, ys, zs = np.linspace(-L, L, nx + 1), np.array([0.0, 1.0]), np.linspace(0.0, H, nz + 1)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    P = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    nid = lambda i, j, k: (i * 2 + j) * (nz + 1) + k
    E = np.array([[nid(i, 0, k), nid(i + 1, 0, k), nid(i + 1, 1, k), nid(i, 1, k),
                   nid(i, 0, k + 1), nid(i + 1, 0, k + 1), nid(i + 1, 1, k + 1), nid(i, 1, k + 1)]
                  for i in range(nx) for k in range(nz)])
    Gs, Ms = G * cfactor(beta), (K + 4.0 * G / 3.0) * cfactor(beta)
    r = solid_kernel(P[E], Ms - 2.0 * Gs, Gs, rho, eint=0, incompatible=incompatible)
    dofs = (3 * E[:, :, None] + np.arange(3)[None, None, :]).reshape(len(E), 24)
    I, J = np.repeat(dofs, 24, axis=1).ravel(), np.tile(dofs, (1, 24)).ravel()
    n = 3 * len(P)
    Kg = sp.csr_matrix((r["K"].ravel(), (I, J)), shape=(n, n))
    Mg = sp.csr_matrix((r["M"].ravel(), (I, J)), shape=(n, n))
    pres = {}
    for a, (x, y, z) in enumerate(P):
        if abs(abs(x) - L) < 1e-9:
            pres[3 * a] = 1.0
        if y in (0.0, 1.0):
            pres[3 * a + 1] = 0.0
        if z == 0.0:
            pres[3 * a + 2] = 0.0
    p = np.array(sorted(pres))
    up = np.array([pres[k] for k in p])
    f = np.setdiff1d(np.arange(n), p)
    out = []
    for fr in freqs:
        D = (Kg - (2 * np.pi * fr) ** 2 * Mg).tocsc()
        u = np.zeros(n, complex)
        u[p] = up
        u[f] = spla.spsolve(D[f][:, f], -D[f][:, p] @ up)
        out.append(abs((Mg @ u)[0::3].sum()) / (rho * 2 * L * H))
    return np.array(out)


@pytest.mark.parametrize("L", [2.5, 5.0, 10.0])
def test_fillpool_water_gives_the_impulsive_mass(L):
    """G = 1e-8 K with incompatible modes: within 3 % of the exact impulsive mass from 1 to 10 Hz (H = 5 m,
    0.5 m elements) -- far below the compression frequency c/4H = 74 Hz."""
    H = 5.0
    ref = wl.impulsive_ratio_exact(L, H)
    nx = int(round(2 * L / 0.5))
    got = effective_mass_ratio(L, H, nx, 10, wl.WATER_G_OVER_K, [1.0, 2.0, 5.0, 10.0])
    assert np.all(np.abs(got / ref - 1.0) < 0.03), (got, ref)
    assert np.all(got > ref)                     # the FE solution converges from above (stiffer mesh)


def test_mesh_refinement_converges_to_the_exact_series():
    ref = wl.impulsive_ratio_exact(5.0, 5.0)
    e8 = effective_mass_ratio(5.0, 5.0, 8, 8, 1e-8, [2.0])[0] / ref - 1.0
    e16 = effective_mass_ratio(5.0, 5.0, 16, 16, 1e-8, [2.0])[0] / ref - 1.0
    assert 0 < e16 < e8 < 0.06 and e8 / e16 > 2.5          # about O(h^2)
    assert e16 < 0.02


def test_without_incompatible_modes_the_water_locks():
    """MOPT,1 (incompatible modes suppressed): the volumetric constraint of the 2x2x2-integrated trilinear
    SOLID locks the nearly incompressible water: wrong (here resonant) effective masses in the seismic band."""
    ref = wl.impulsive_ratio_exact(5.0, 5.0)
    got = effective_mass_ratio(5.0, 5.0, 16, 16, 1e-8, [1.0, 2.0, 5.0], incompatible=False)
    assert np.max(np.abs(got / ref - 1.0)) > 0.5, got


def test_nu_049_water_moves_with_the_walls():
    """nu = 0.49 (G = 0.0201 K): an elastic solid with shear modes near 10 Hz for 5 m of water; below them it
    moves with the tank (effective mass = the whole water), which is not the hydrodynamic behaviour."""
    G_over_K = 3 * (1 - 2 * 0.49) / (2 * (1 + 0.49))
    got = effective_mass_ratio(5.0, 5.0, 16, 16, G_over_K, [0.5, 1.0, 2.0])
    assert np.all(np.abs(got - 1.0) < 0.01), got
    ref = wl.impulsive_ratio_exact(5.0, 5.0)
    assert got[0] / ref > 1.9


def test_compressibility_raises_the_effective_mass_towards_the_acoustic_frequency():
    """Above ~20 Hz the water column's compression (first frequency c/4H = 74 Hz for H = 5 m) adds to the
    impulsive mass (the incompressible reference holds only well below it)."""
    got = effective_mass_ratio(5.0, 5.0, 16, 16, 1e-8, [2.0, 30.0])
    assert got[1] > got[0] * 1.04
