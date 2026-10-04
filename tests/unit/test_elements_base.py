"""Unit tests of sassi.elements.base: materials (UT-14), sections (UT-14), complex moduli
(requirements 4.0.2, D-CNV-04), mass schemes (D-CNV-05) and helpers."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.conventions import CMOD_SIMPLE, cfactor
from sassi.elements import ELEMENTS
from sassi.elements.base import (ElementError, as_bool6, frame_from_three_points,
                                 material_from_layer, material_from_M, mixed_mass, nodal_mass_to_mass_units,
                                 rectangle_torsion_constant, section_circle, section_rectangle)

# numpy 2.0 + macOS Accelerate raises spurious FP flags in matmul (results are finite and exact)
pytestmark = pytest.mark.filterwarnings("ignore:.*encountered in matmul:RuntimeWarning")


# ---------------------------------------------------------------- UT-14 materials
def test_ut14_material_types_agree():
    E, nu, gam, g = 30000.0, 0.25, 0.15, 32.2
    rho = gam / g
    m1 = material_from_M(1, E, nu, gam, 0.0, 0.0, g)
    G = E / (2 * (1 + nu))
    M = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    m2 = material_from_M(2, M, G, gam, 0.0, 0.0, g)
    m3 = material_from_M(3, math.sqrt(M / rho), math.sqrt(G / rho), gam, 0.0, 0.0, g)
    for m in (m1, m2, m3):
        assert m.G0 == pytest.approx(12000.0, rel=1e-10)
        assert m.M0 == pytest.approx(36000.0, rel=1e-10)
        assert m.vs == pytest.approx(math.sqrt(G / rho), rel=1e-10)
        assert m.vp == pytest.approx(math.sqrt(M / rho), rel=1e-10)
        assert m.E0 == pytest.approx(E, rel=1e-10)
        assert m.nu0 == pytest.approx(nu, rel=1e-10)
        assert m.rho == pytest.approx(rho, rel=1e-14)


def test_ut14_sections():
    c = section_circle(1.0)
    assert c["I2"] == pytest.approx(0.785398, abs=5e-7)
    assert c["I3"] == pytest.approx(0.785398, abs=5e-7)
    assert c["J"] == pytest.approx(1.570796, abs=5e-7)
    assert c["As2"] == pytest.approx(math.pi * 0.9)
    r = section_rectangle(0.5, 1.0)
    assert r["J"] == pytest.approx(0.028610, abs=5e-7)
    # manual example R,1,0.5,0.41667,0.41667,0.02861,0.010417,0.041667
    assert r["A"] == 0.5
    assert r["As2"] == pytest.approx(0.41667, abs=5e-6)
    assert r["I2"] == pytest.approx(0.010417, abs=5e-7)
    assert r["I3"] == pytest.approx(0.041667, abs=5e-7)
    assert rectangle_torsion_constant(1.0, 0.5) == rectangle_torsion_constant(0.5, 1.0)   # b > h swapped


# ---------------------------------------------------------------- complex moduli
def test_complex_moduli_sassi_form():
    m = material_from_M(1, 1000.0, 0.3, 2.0, 0.04, 0.02, 1.0)
    assert m.G == pytest.approx(m.G0 * cfactor(0.02))
    assert m.M == pytest.approx(m.M0 * cfactor(0.04))
    assert m.lam == pytest.approx(m.M - 2 * m.G)
    assert m.nu == pytest.approx((m.M - 2 * m.G) / (2 * (m.M - m.G)))
    assert m.E == pytest.approx(2 * m.G * (1 + m.nu))


def test_beam_and_shell_moduli_equal_E_c_when_dampings_equal():
    b = 0.05
    m = material_from_M(1, 1000.0, 0.3, 2.0, b, b, 1.0)
    E, G = m.beam_moduli()
    assert E == pytest.approx(1000.0 * cfactor(b), rel=1e-13)
    assert G == pytest.approx(m.G0 * cfactor(b), rel=1e-13)
    Es, nus = m.shell_moduli()
    assert Es == pytest.approx(1000.0 * cfactor(b), rel=1e-13)
    assert isinstance(nus, float) and nus == pytest.approx(0.3)


def test_cmodform_simple_and_undamped_copy():
    m = material_from_M(1, 1000.0, 0.25, 1.0, 0.1, 0.1, 1.0, form=CMOD_SIMPLE)
    assert m.G == pytest.approx(m.G0 * (1 + 0.2j))
    u = m.undamped()
    assert u.G == m.G0 and u.M == m.M0 and u.beta_s == 0.0


def test_material_from_layer():
    m = material_from_layer(400.0, 200.0, 19.62, 0.02, 0.05, 9.81)
    assert m.rho == pytest.approx(2.0)
    assert m.G0 == pytest.approx(2.0 * 200.0 ** 2)
    assert m.M0 == pytest.approx(2.0 * 400.0 ** 2)
    assert m.G == pytest.approx(m.G0 * cfactor(0.05))
    assert m.M == pytest.approx(m.M0 * cfactor(0.02))


@pytest.mark.parametrize("args", [(1, -1.0, 0.3), (1, 1.0, 0.5), (2, 1.0, 1.0), (3, 1.0, 1.0), (4, 1.0, 0.2)])
def test_material_errors(args):
    with pytest.raises(ElementError):
        material_from_M(args[0], args[1], args[2], 1.0, 0.0, 0.0, 1.0)


def test_damping_limit_edu04():
    with pytest.raises(ElementError, match="EDU-04"):
        material_from_M(1, 1.0, 0.2, 1.0, 0.5, 0.0, 1.0)


# ---------------------------------------------------------------- helpers
def test_mixed_mass_definition():
    Mc = np.array([[2.0, 1.0], [1.0, 2.0]]) / 6.0
    Ml = mixed_mass(Mc, "lumped")
    assert np.allclose(Ml, np.diag([0.5, 0.5]))
    assert np.allclose(mixed_mass(Mc), 0.5 * Ml + 0.5 * Mc)
    assert np.allclose(mixed_mass(Mc, "consistent"), Mc)
    with pytest.raises(ElementError):
        mixed_mass(Mc, "diagonal")


def test_frame_from_three_points_and_error9():
    Lam, L = frame_from_three_points([0, 0, 0], [0, 2, 0], [0, 0, 5])
    assert L == 2.0
    assert np.allclose(Lam, [[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    with pytest.raises(ElementError, match="Error 9"):
        frame_from_three_points([0, 0, 0], [1, 0, 0], [3, 0, 0])
    with pytest.raises(ElementError, match="Error 8"):
        frame_from_three_points([1, 1, 1], [1, 1, 1], [0, 0, 0])


def test_release_codes():
    assert as_bool6("000011") == (False, False, False, False, True, True)
    assert as_bool6("11") == (False, False, False, False, True, True)
    assert as_bool6([0, 1, 0, 0, 0, 0])[1] is True
    assert as_bool6(None) == (False,) * 6
    with pytest.raises(ElementError):
        as_bool6("000021")


def test_nodal_mass_units():
    assert np.allclose(nodal_mass_to_mass_units([32.2, 32.2, 0, 0, 0, 0], 1, 32.2), [1, 1, 0, 0, 0, 0])
    assert np.allclose(nodal_mass_to_mass_units([1, 1, 0, 0, 0, 0], 0, 32.2), [1, 1, 0, 0, 0, 0])


def test_registry():
    assert sorted(ELEMENTS) == [1, 2, 3, 4, 5, 7, 9]     # TSHELL (5) registered in wave 3
    assert ELEMENTS[1].name == "SOLID" and ELEMENTS[1].dofs == (1, 2, 3) and ELEMENTS[1].nnodes == 8
    assert ELEMENTS[2].nnodes == 2 and ELEMENTS[2].dofs == (1, 2, 3, 4, 5, 6)
    assert ELEMENTS[4].dofs == (1, 3)
    assert ELEMENTS[3].components == ["FXX", "FYY", "FXY", "MXX", "MYY", "MXY"]
    assert ELEMENTS[7].components == ["FX", "FY", "FZ", "MXX", "MYY", "MZZ"]
    assert ELEMENTS[9].components == []
