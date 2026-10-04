"""Element library foundation: materials, section properties, the element registry and
small numerical helpers shared by the HOUSE finite elements.

Normative sources
-----------------
* requirements section 1.5 (element table), 4.0.2 (complex moduli per element type),
  4.1 (formulations), decisions D-CNV-03/04/05, D-ELM-01..09, D-MDL-08;
* spec 08 section 5.20 (M material types 1/2/3), 6.1 (conversions), 6.3 (section formulas);
* R1 section 1 (complex-modulus damping convention).

Complex moduli (requirements 4.0.2, D-CNV-04)
---------------------------------------------
Hysteretic damping is represented by complex moduli.  With the SASSI/SHAKE91 factor
``c(b) = 1 - 2b^2 + 2ib sqrt(1-b^2)`` (|c| = 1) the shear modulus carries the S-wave damping and
the *constrained* modulus ``M = lam + 2G = rho Vp^2`` carries the P-wave damping::

    G* = G c(beta_s)        M* = M c(beta_p)        lam* = M* - 2 G*

Every element type derives its complex stiffness from these two numbers:

* SOLID / PLANE : complex Lame constants (lam*, G*);
* BEAMS         : E* and G* from (M*, G*) with the type-2 formulas, so E* = E c(beta) if beta_p = beta_s;
* SHELL         : E* = E c(beta), nu real (CHECK warns when beta_p != beta_s; beta_s is used);
* SPRING        : k* = k c(damp);  GENERAL: K_R + i K_I as entered.

The mass is always real (rho = weight / gravity).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, Optional, Sequence, Tuple

import numpy as np

from ..conventions import CMOD_SASSI, ELEMENT_COMPONENTS, cfactor

__all__ = [
    "ElementError", "MaterialProps", "material_from_M", "material_from_layer", "material_from_E_nu",
    "ElementSpec", "ELEMENTS", "register", "section_circle", "section_rectangle",
    "rectangle_torsion_constant", "nodal_mass_to_mass_units", "gauss_legendre", "mixed_mass",
    "lumped_from_consistent", "frame_from_three_points", "cfactor",
]


class ElementError(ValueError):
    """Invalid element data or geometry (the message cites the CHECK/EDU code when one applies)."""


# ======================================================================================
# Materials
# ======================================================================================
@dataclass
class MaterialProps:
    """Mass density and complex moduli of one material (ARCHITECTURE section 6.1).

    Attributes
    ----------
    rho : mass density (weight / gravity)
    G   : complex shear modulus G* = G0 c(beta_s)
    M   : complex constrained modulus M* = M0 c(beta_p) = lam* + 2 G*
    beta_s, beta_p : hysteretic damping ratios (fractions)
    G0, M0 : real (undamped) moduli
    form : complex-modulus form (0 SASSI default, 1 = 1 + 2 i beta; CMODFORM)
    """

    rho: float
    G: complex
    M: complex
    beta_s: float
    beta_p: float
    G0: float
    M0: float
    form: int = CMOD_SASSI

    # ---- complex derived constants (requirements 4.0.2) --------------------------------
    @property
    def lam(self) -> complex:
        """Complex Lame constant lam* = M* - 2 G*."""
        return self.M - 2.0 * self.G

    @property
    def nu(self) -> complex:
        """Complex Poisson ratio nu* = (M* - 2G*) / (2 (M* - G*)) (type-2 formula)."""
        return (self.M - 2.0 * self.G) / (2.0 * (self.M - self.G))

    @property
    def E(self) -> complex:
        """Complex Young's modulus E* = 2 G* (1 + nu*) (equals E0 c(beta) when beta_p = beta_s)."""
        return 2.0 * self.G * (1.0 + self.nu)

    # ---- real (undamped) constants ---------------------------------------------------------
    @property
    def lam0(self) -> float:
        return self.M0 - 2.0 * self.G0

    @property
    def nu0(self) -> float:
        return (self.M0 - 2.0 * self.G0) / (2.0 * (self.M0 - self.G0))

    @property
    def E0(self) -> float:
        return 2.0 * self.G0 * (1.0 + self.nu0)

    @property
    def vs(self) -> float:
        """S-wave velocity sqrt(G0 / rho)."""
        return math.sqrt(self.G0 / self.rho) if self.rho > 0 else float("inf")

    @property
    def vp(self) -> float:
        """P-wave velocity sqrt(M0 / rho)."""
        return math.sqrt(self.M0 / self.rho) if self.rho > 0 else float("inf")

    # ---- per element-type moduli -------------------------------------------------------
    def shell_moduli(self) -> Tuple[complex, float]:
        """(E*, nu) for SHELL elements: E* = E0 c(beta_s), nu = nu0 real (requirements 4.0.2).

        The manual requires beta_p = beta_s for shell materials; CHECK warns otherwise and the
        shear damping beta_s is used here.
        """
        return self.E0 * cfactor(self.beta_s, self.form), self.nu0

    def beam_moduli(self) -> Tuple[complex, complex]:
        """(E*, G*) for BEAMS elements, derived from (M*, G*) (D-CNV-04)."""
        return self.E, self.G

    def undamped(self) -> "MaterialProps":
        """Copy with zero damping (real moduli): used for modal checks and mass transformations."""
        return MaterialProps(self.rho, complex(self.G0), complex(self.M0), 0.0, 0.0, self.G0, self.M0, self.form)


def _check_damping(beta: float, what: str) -> None:
    if beta < 0:
        raise ElementError(f"negative {what} damping ratio {beta}")
    if beta >= 0.5:
        raise ElementError(f"{what} damping ratio {beta} >= 0.5 (EDU-04)")


def _make(rho: float, G0: float, M0: float, pdamp: float, sdamp: float, form: int) -> MaterialProps:
    _check_damping(pdamp, "P-wave")
    _check_damping(sdamp, "S-wave")
    G = complex(G0 * cfactor(sdamp, form))
    M = complex(M0 * cfactor(pdamp, form))
    return MaterialProps(rho=float(rho), G=G, M=M, beta_s=float(sdamp), beta_p=float(pdamp),
                         G0=float(G0), M0=float(M0), form=int(form))


def material_from_M(mtype: int, val1: float, val2: float, weight: float, pdamp: float, sdamp: float,
                    gravity: float, form: int = CMOD_SASSI) -> MaterialProps:
    """Material of the M table (spec 08 section 5.20 and 6.1).

    ``mtype`` 1: (E, nu); 2: (constrained modulus M = lam + 2G, shear modulus G); 3: (Vp, Vs).
    ``weight`` is the specific weight; rho = weight / gravity (requirements 4.0.3).
    """
    if gravity <= 0:
        raise ElementError("gravity must be > 0")
    rho = float(weight) / float(gravity)
    mtype = int(mtype) if mtype else 1
    if mtype == 1:
        E, nu = float(val1), float(val2)
        if E <= 0:
            raise ElementError("Error 14: Young's modulus must be > 0")
        if not (-1.0 < nu < 0.5):
            raise ElementError(f"Poisson's ratio {nu} outside (-1, 0.5)")
        G0 = E / (2.0 * (1.0 + nu))
        M0 = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
    elif mtype == 2:
        M0, G0 = float(val1), float(val2)
        if G0 <= 0 or M0 <= 4.0 * G0 / 3.0:
            raise ElementError("type-2 material needs G > 0 and M > 4G/3 (positive bulk modulus)")
    elif mtype == 3:
        vp, vs = float(val1), float(val2)
        if vs <= 0 or vp <= 0:
            raise ElementError("type-3 material needs Vp > 0 and Vs > 0")
        G0 = rho * vs * vs
        M0 = rho * vp * vp
        if M0 <= 4.0 * G0 / 3.0:
            raise ElementError("type-3 material needs Vp > (2/sqrt(3)) Vs (positive bulk modulus)")
    else:
        raise ElementError(f"material type {mtype} must be 1, 2 or 3")
    return _make(rho, G0, M0, pdamp, sdamp, form)


def material_from_E_nu(E: float, nu: float, rho: float, beta: float = 0.0,
                       form: int = CMOD_SASSI) -> MaterialProps:
    """Convenience constructor in mass-density units (tests, verification problems)."""
    return material_from_M(1, E, nu, rho, beta, beta, 1.0, form)


def material_from_layer(vp: float, vs: float, weight: float, dp: float, ds: float, gravity: float,
                        form: int = CMOD_SASSI) -> MaterialProps:
    """Material of an L soil layer for excavated SOLID/PLANE elements (spec 08 section 6.2,
    D-ELM-12): rho = weight/g, G = rho Vs^2, M = rho Vp^2, damping dp on M and ds on G."""
    if gravity <= 0:
        raise ElementError("gravity must be > 0")
    rho = float(weight) / float(gravity)
    return _make(rho, rho * vs * vs, rho * vp * vp, dp, ds, form)


def nodal_mass_to_mass_units(values: Sequence[float], units: int, gravity: float) -> np.ndarray:
    """Nodal masses MT/MR in mass units: divided by g when MUNITS = 1 (weight units).

    requirements 4.0.3 and D-MDL-08 (default MUNITS = 1).  UT-15.
    """
    v = np.asarray(values, dtype=float)
    if int(units) == 1:
        if gravity <= 0:
            raise ElementError("gravity must be > 0 to convert weights to masses")
        return v / float(gravity)
    return v.copy()


# ======================================================================================
# Beam section formulas (spec 08 section 6.3; UT-14)
# ======================================================================================
def section_circle(r: float) -> Dict[str, float]:
    """Solid circle of radius r: A = pi r^2, I2 = I3 = pi r^4/4, J = pi r^4/2, As = A/f, f = 10/9."""
    A = math.pi * r * r
    I = math.pi * r ** 4 / 4.0
    return dict(A=A, As2=A * 0.9, As3=A * 0.9, J=2.0 * I, I2=I, I3=I)


def rectangle_torsion_constant(b: float, h: float) -> float:
    """Roark: J = [1/3 - 0.21 (b/h)(1 - b^4/(12 h^4))] h b^3 with b <= h (spec 08 Q16)."""
    b, h = (b, h) if b <= h else (h, b)
    return (1.0 / 3.0 - 0.21 * (b / h) * (1.0 - b ** 4 / (12.0 * h ** 4))) * h * b ** 3


def section_rectangle(b: float, h: float) -> Dict[str, float]:
    """Solid rectangle, b along local axis 3 and h along local axis 2 (spec 08 section 6.3).

    I2 = h b^3/12 (about axis 2), I3 = b h^3/12 (about axis 3), As = A/f with f = 6/5.
    """
    A = b * h
    return dict(A=A, As2=A / 1.2, As3=A / 1.2, J=rectangle_torsion_constant(b, h),
                I2=h * b ** 3 / 12.0, I3=b * h ** 3 / 12.0)


# ======================================================================================
# Element registry (ARCHITECTURE section 6.2)
# ======================================================================================
@dataclass
class ElementSpec:
    """Description of one HOUSE element type.

    ``matrices(xyz, mat, **props) -> (K complex (nd, nd), M real (nd, nd))`` and
    ``recovery(xyz, mat, **props) -> S complex (ncomp, nd)`` with the element DOF vector in
    node-major order ``[node1 dofs..., node2 dofs..., ...]`` (``dofs`` order).

    ``batch`` (optional, SASSI-EDU addition) evaluates many elements of the type at once:
    ``batch(xyz (nE, nn, 3), mats: list, **props) -> dict(K=(nE,nd,nd), M=(nE,nd,nd), S=(nE,nc,nd)
    [, B=(nE,nstrain,nd) strain operator])``.
    """

    code: int
    name: str
    nnodes: int
    dofs: tuple
    components: list
    matrices: Callable
    recovery: Callable
    batch: Optional[Callable] = None
    description: str = ""

    @property
    def ndof(self) -> int:
        return self.nnodes * len(self.dofs)


ELEMENTS: Dict[int, ElementSpec] = {}


def register(spec: ElementSpec) -> ElementSpec:
    """Add an element type to :data:`ELEMENTS` (called by each element module)."""
    if spec.components is None:
        spec.components = list(ELEMENT_COMPONENTS.get(spec.name, []))
    ELEMENTS[spec.code] = spec
    return spec


# ======================================================================================
# Numerical helpers
# ======================================================================================
def gauss_legendre(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """n-point Gauss-Legendre rule on [-1, 1] (points, weights)."""
    x, w = np.polynomial.legendre.leggauss(int(n))
    return x, w


def lumped_from_consistent(Mc: np.ndarray) -> np.ndarray:
    """Row-sum lumped mass: diag(sum_j Mc_ij) (works on stacks (..., n, n))."""
    Mc = np.asarray(Mc)
    rs = Mc.sum(axis=-1)
    out = np.zeros_like(Mc)
    idx = np.arange(Mc.shape[-1])
    out[..., idx, idx] = rs
    return out


def mixed_mass(Mc: np.ndarray, scheme: str = "mixed") -> np.ndarray:
    """Element mass by scheme (D-CNV-05): 'mixed' = 1/2 lumped (row sum) + 1/2 consistent,
    'consistent', or 'lumped'.

    Why the mix: for linear elements the consistent mass over-estimates and the lumped mass
    under-estimates the frequencies; their average cancels the leading (kh)^2 dispersion error
    and leaves a 4th-order error ``1 - (kh)^4/480`` (R2 A.3).  This is why the SASSI lambda/5
    mesh rule is sufficient.
    """
    scheme = (scheme or "mixed").lower()
    if scheme == "consistent":
        return np.array(Mc, dtype=float, copy=True)
    ML = lumped_from_consistent(Mc)
    if scheme == "lumped":
        return ML
    if scheme == "mixed":
        return 0.5 * ML + 0.5 * np.asarray(Mc, dtype=float)
    raise ElementError(f"unknown mass scheme '{scheme}' (mixed, consistent, lumped)")


def frame_from_three_points(xi: np.ndarray, xj: np.ndarray, xk: np.ndarray) -> Tuple[np.ndarray, float]:
    """Local axes of BEAMS and 3-node GENERAL elements (spec 08 section 4.5).

    e1 = (X_J - X_I)/L;  e2 = component of (X_K - X_I) normal to e1, normalised;  e3 = e1 x e2.
    Returns (Lam, L) with the rows of ``Lam`` equal to e1, e2, e3 (local = Lam @ global).
    Raises ElementError (CHECK Error 9) when K lies on the I-J line, and Error 8 for zero length.
    """
    xi = np.asarray(xi, float)
    d = np.asarray(xj, float) - xi
    L = float(np.linalg.norm(d))
    if L <= 0.0:
        raise ElementError("Error 8: zero-length element (I and J coincide)")
    e1 = d / L
    v = np.asarray(xk, float) - xi
    v2 = v - np.dot(v, e1) * e1
    n2 = float(np.linalg.norm(v2))
    scale = max(L, float(np.linalg.norm(v)), 1e-300)
    if n2 <= 1e-9 * scale:
        raise ElementError("Error 9: orientation node K is collinear with I-J")
    e2 = v2 / n2
    e3 = np.cross(e1, e2)
    return np.vstack([e1, e2, e3]), L


def as_bool6(code) -> Tuple[bool, ...]:
    """Normalise an end-release code: '000011', (0,0,0,0,1,1), None -> 6 bools (1 = released)."""
    if code is None:
        return (False,) * 6
    if isinstance(code, str):
        s = code.strip()
        if s == "":
            return (False,) * 6
        s = s.zfill(6) if len(s) < 6 else s
        if len(s) != 6 or any(c not in "01" for c in s):
            raise ElementError(f"release code '{code}' must be six 0/1 digits")
        return tuple(c == "1" for c in s)
    vals = list(code)
    if len(vals) != 6:
        raise ElementError(f"release code {code!r} must have six entries")
    out = []
    for v in vals:
        iv = int(v)
        if iv not in (0, 1):
            raise ElementError(f"release code {code!r}: values must be 0 or 1")
        out.append(bool(iv))
    return tuple(out)


def blas_quiet():
    """Context manager silencing spurious floating-point flags of matmul.

    numpy 2.0 linked against the macOS Accelerate BLAS raises "divide by zero / overflow /
    invalid value encountered in matmul" warnings for some perfectly finite complex products.
    The element kernels run their linear algebra inside this context and then verify the
    results with :func:`check_finite`, so genuine problems are still reported.
    """
    return np.errstate(divide="ignore", over="ignore", invalid="ignore", under="ignore")


def check_finite(a: np.ndarray, what: str) -> np.ndarray:
    """Raise ElementError if ``a`` contains NaN or Inf (singular/degenerate element data)."""
    if not np.all(np.isfinite(a)):
        raise ElementError(f"non-finite values in {what} (degenerate geometry or material?)")
    return a
