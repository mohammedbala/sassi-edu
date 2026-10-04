"""Property-table physics: material conversions, section formulas, weight -> mass (spec 08 section 6).

The interpreter stores the *raw* user values of the M, L, R, SC and MX tables (spec 08 section 1.6:
conversions happen when the analysis files are written, because the gravity may be set after
the materials).  This module provides the conversions so that MLIST, AFWRITE and the tests use
one implementation.

Material types of ``M,<nm>,<val1>,<val2>,<weight>,<pdamp>,<sdamp>,<type>`` (spec 08 section 6.1),
with mass density ``rho = weight / g`` (``weight`` is a *specific weight*, force/volume):

======  ==================  ==========================================================
type    given               derived
======  ==================  ==========================================================
1       E, nu               G = E/(2(1+nu)); M = E(1-nu)/((1+nu)(1-2nu))
2       M (= lam+2G), G     nu = (M-2G)/(2(M-G)); E = G(3M-4G)/(M-G)
3       Vp, Vs              G = rho Vs^2; M = rho Vp^2; nu = (Vp^2-2Vs^2)/(2(Vp^2-Vs^2))
======  ==================  ==========================================================

and ``Vs = sqrt(G/rho)``, ``Vp = sqrt(M/rho)``, ``lam = M - 2G``.  Damping is applied
separately with the complex modulus of :func:`sassi.conventions.complex_lame` (D-CNV-03/04).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict

G_DEFAULT = 32.2   # HOUSE <gravity> default (sassi.io.decks HOUSE schema)


@dataclass(frozen=True)
class ElasticConstants:
    """Real (undamped) elastic constants of a material or soil layer."""
    E: float
    nu: float
    G: float
    M: float        # constrained (P-wave) modulus lam + 2G
    lam: float
    Vp: float
    Vs: float
    rho: float      # mass density = weight / g

    def as_dict(self) -> Dict[str, float]:
        return dict(E=self.E, nu=self.nu, G=self.G, M=self.M, lam=self.lam, Vp=self.Vp, Vs=self.Vs,
                    rho=self.rho)


def _sqrt(x: float) -> float:
    return math.sqrt(x) if x > 0 else 0.0


def elastic_constants(mtype: int, val1: float, val2: float, weight: float, gravity: float) -> ElasticConstants:
    """Convert the raw M-table values of material type 1, 2 or 3 (spec 08 section 6.1).

    Raises ``ValueError`` for an unknown type, a non-positive gravity, or values that make a
    conversion singular (e.g. nu = 0.5 for type 1, M = G for type 2, Vp = Vs for type 3).
    """
    if gravity <= 0:
        raise ValueError("gravity must be > 0 to convert weight to mass (CHECK Error 1)")
    rho = weight / gravity
    if mtype == 1:
        E, nu = float(val1), float(val2)
        if (1.0 + nu) == 0 or (1.0 - 2.0 * nu) == 0:
            raise ValueError("Poisson's ratio gives a singular conversion")
        G = E / (2.0 * (1.0 + nu))
        M = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
    elif mtype == 2:
        M, G = float(val1), float(val2)
        if M == G:
            raise ValueError("M = G gives a singular conversion")
        nu = (M - 2.0 * G) / (2.0 * (M - G))
        E = G * (3.0 * M - 4.0 * G) / (M - G)
    elif mtype == 3:
        vp, vs = float(val1), float(val2)
        if vp * vp == vs * vs:
            raise ValueError("Vp = Vs gives a singular conversion")
        G = rho * vs * vs
        M = rho * vp * vp
        nu = (vp * vp - 2.0 * vs * vs) / (2.0 * (vp * vp - vs * vs))
        E = 2.0 * G * (1.0 + nu)
    else:
        raise ValueError(f"material type {mtype} is not 1, 2 or 3")
    lam = M - 2.0 * G
    if rho > 0:
        Vs, Vp = _sqrt(G / rho), _sqrt(M / rho)
    else:
        Vs = Vp = 0.0
    return ElasticConstants(E=E, nu=nu, G=G, M=M, lam=lam, Vp=Vp, Vs=Vs, rho=rho)


def layer_constants(vp: float, vs: float, weight: float, gravity: float) -> ElasticConstants:
    """Elastic constants of an L-table soil layer (spec 08 section 6.2): G = rho Vs^2, M = rho Vp^2."""
    return elastic_constants(3, vp, vs, weight, gravity)


# --------------------------------------------------------------------------------------
# Beam sections (R table, spec 08 section 6.3; manual figure of section 9.4.30)
# --------------------------------------------------------------------------------------
def section_circle(r: float) -> Dict[str, float]:
    """Solid circle of radius r: A = pi r^2, I2 = I3 = pi r^4/4, J = pi r^4/2, As = A/(10/9)."""
    A = math.pi * r * r
    I = math.pi * r ** 4 / 4.0
    As = A / (10.0 / 9.0)
    return dict(axial=A, shear2=As, shear3=As, tors=2.0 * I, flex2=I, flex3=I)


def section_rectangle(b: float, h: float) -> Dict[str, float]:
    """Solid rectangle, ``b`` along local axis 3 and ``h`` along local axis 2 (spec 08 section 6.3).

    I2 = h b^3/12 (about axis 2), I3 = b h^3/12 (about axis 3), shear form factor 6/5, and the
    Roark torsion constant ``J = [1/3 - 0.21 (b/h)(1 - b^4/(12 h^4))] h b^3`` with b <= h (the
    manual figure prints ``0.21(bh)``, a typo for b/h -- spec 08 Q16; b and h are swapped in the
    J formula only when b > h).
    """
    A = b * h
    I2 = h * b ** 3 / 12.0
    I3 = b * h ** 3 / 12.0
    s, l = (b, h) if b <= h else (h, b)
    J = (1.0 / 3.0 - 0.21 * (s / l) * (1.0 - s ** 4 / (12.0 * l ** 4))) * l * s ** 3
    As = A / 1.2
    return dict(axial=A, shear2=As, shear3=As, tors=J, flex2=I2, flex3=I3)


# --------------------------------------------------------------------------------------
# Mass units (MUNITS, MOPT <matrix>)
# --------------------------------------------------------------------------------------
def to_mass(value: float, units: int, gravity: float) -> float:
    """Nodal (or GENERAL matrix) mass in mass units: ``value / g`` when given as weight (units 1).

    MUNITS 0 = mass units, 1 = weight units (spec 08 section 7.19, D-MDL-08); MOPT ``<matrix>``
    uses the same codes for the GENERAL element mass matrix (MXM).  Rotational "weights" are
    weight x length^2 and are divided by g in the same way.
    """
    if units == 1:
        if gravity <= 0:
            raise ValueError("gravity must be > 0 to convert weight to mass")
        return value / gravity
    return value
