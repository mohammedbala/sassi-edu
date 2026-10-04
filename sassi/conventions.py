"""Global technical conventions shared by every SASSI-EDU module.

Normative source: docs/spec/00_requirements.md section 0.3 and section 4.0.

* Harmonic factor ``exp(+i w t)``; horizontally outgoing waves ``exp(-i k x)``;
  cylindrical outgoing waves use Hankel functions of the second kind ``H^(2)``.
* FFT: ``numpy.fft.rfft`` (forward, no scaling) and ``numpy.fft.irfft`` (inverse, 1/N).
* Complex modulus (SASSI / SHAKE91 form, default)::

      c(beta) = 1 - 2 beta^2 + 2 i beta sqrt(1 - beta^2)       |c| = 1

  ``CMODFORM,1`` selects the simpler ``1 + 2 i beta`` form (benchmarks only).
* Geometry: right-handed global Cartesian, Z up.  2D models lie in the X-Z plane.
* DOF labels per node: 1 UX, 2 UY, 3 UZ, 4 ROTX, 5 ROTY, 6 ROTZ.
"""
from __future__ import annotations

import math
from typing import Union

import numpy as np

ArrayLike = Union[float, np.ndarray]

# --------------------------------------------------------------------------------------
# Complex modulus
# --------------------------------------------------------------------------------------
CMOD_SASSI = 0   # 1 - 2b^2 + 2ib sqrt(1-b^2)  (default)
CMOD_SIMPLE = 1  # 1 + 2ib


def cfactor(beta: ArrayLike, form: int = CMOD_SASSI) -> ArrayLike:
    """Return the complex-modulus factor c(beta) for hysteretic damping ratio ``beta``.

    ``G* = G * cfactor(beta)``.  ``form`` is ``CMOD_SASSI`` (default) or ``CMOD_SIMPLE``.
    """
    b = np.asarray(beta, dtype=float)
    if np.any(b < 0) or np.any(b >= 0.5):
        raise ValueError("damping ratio must satisfy 0 <= beta < 0.5 (EDU-04)")
    if form == CMOD_SIMPLE:
        out = 1.0 + 2.0j * b
    else:
        out = 1.0 - 2.0 * b * b + 2.0j * b * np.sqrt(1.0 - b * b)
    if np.ndim(out) == 0:
        return complex(out)
    return out


def complex_velocity(v: ArrayLike, beta: ArrayLike, form: int = CMOD_SASSI) -> ArrayLike:
    """Complex wave velocity ``V* = V sqrt(c(beta))`` (principal root)."""
    return np.asarray(v) * np.sqrt(cfactor(beta, form))


# --------------------------------------------------------------------------------------
# Elastic constant conversions (real parts; damping applied separately)
# --------------------------------------------------------------------------------------
def moduli_from_velocities(vp: float, vs: float, rho: float):
    """Return (G, M, lam, nu, E) from P/S velocities and mass density."""
    G = rho * vs * vs
    M = rho * vp * vp
    lam = M - 2.0 * G
    nu = lam / (2.0 * (lam + G)) if (lam + G) != 0 else 0.0
    E = 2.0 * G * (1.0 + nu)
    return G, M, lam, nu, E


def moduli_from_E_nu(E: float, nu: float):
    """Return (G, M, lam) from Young's modulus and Poisson's ratio."""
    G = E / (2.0 * (1.0 + nu))
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    M = lam + 2.0 * G
    return G, M, lam


def complex_lame(G: float, M: float, beta_s: float, beta_p: float, form: int = CMOD_SASSI):
    """Complex (G*, M*, lam*) with beta_s applied to G and beta_p to M = lam + 2G (D-CNV-04)."""
    Gs = G * cfactor(beta_s, form)
    Ms = M * cfactor(beta_p, form)
    return Gs, Ms, Ms - 2.0 * Gs


# --------------------------------------------------------------------------------------
# Units
# --------------------------------------------------------------------------------------
G_BRITISH = 32.174
G_SI = 9.80665


def unit_system(gravity: float) -> str:
    """'BS' (British, ft) when g > 20, else 'SI' (m) (D-CNV-08)."""
    return "BS" if gravity > 20.0 else "SI"


# --------------------------------------------------------------------------------------
# Frequencies
# --------------------------------------------------------------------------------------
def frequency_step(delt: float, nfft: int, fstep: float = 0.0) -> float:
    """Frequency step: ``fstep`` if > 0 else ``1/(delt*nfft)``."""
    if fstep and fstep > 0:
        return float(fstep)
    if delt <= 0 or nfft <= 0:
        raise ValueError("need delt > 0 and nfft > 0 (or fstep > 0)")
    return 1.0 / (delt * nfft)


def nearest_power_of_two(n: int) -> int:
    """Nearest power of two (ties go up) -- AFWRITE rule for Warning 9 (D-CNV-10)."""
    if n <= 1:
        return 1
    lo = 1 << (int(n).bit_length() - 1)
    hi = lo << 1
    return lo if (n - lo) < (hi - n) else hi


def is_power_of_two(n: int) -> bool:
    return n > 0 and (n & (n - 1)) == 0


# --------------------------------------------------------------------------------------
# Naming formulas (requirements section 1.8, normative)
# --------------------------------------------------------------------------------------
DOF_TAGS = {1: "TR_X", 2: "TR_Y", 3: "TR_Z", 4: "R_XX", 5: "R_YY", 6: "R_ZZ"}
DOF_LABELS = {1: "UX", 2: "UY", 3: "UZ", 4: "ROTX", 5: "ROTY", 6: "ROTZ"}
ELEMENT_TYPE_NAMES = {1: "SOLID", 2: "BEAMS", 3: "SHELL", 4: "PLANE", 5: "TSHELL", 7: "SPRING", 9: "GENERAL"}
ELEMENT_TYPE_CODES = {v: k for k, v in ELEMENT_TYPE_NAMES.items()}
ELEMENT_TYPE_CODES.update({"BEAM": 2, "GM": 9, "GENERALMATRIX": 9})

# Element output component codes (D-FIL-04), in EOUT code order.
ELEMENT_COMPONENTS = {
    "SOLID": ["SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ", "SOCT"],
    "BEAMS": ["FXI", "FYI", "FZI", "MXI", "MYI", "MZI", "FXJ", "FYJ", "FZJ", "MXJ", "MYJ", "MZJ"],
    "SHELL": ["FXX", "FYY", "FXY", "MXX", "MYY", "MXY"],
    "TSHELL": ["NXX", "NYY", "NXY", "QXZ", "QYZ", "MXX", "MYY", "MXY"],
    "PLANE": ["SXX", "SZZ", "TXZ"],
    "SPRING": ["FX", "FY", "FZ", "MXX", "MYY", "MZZ"],
}


def node_width(max_node: int) -> int:
    """Digits used in node-based file names: 5, or 6 when any node id exceeds 99,999."""
    return 6 if max_node > 99999 else 5


def nodal_result_name(node: int, dof: int, ext: str, max_node: int = 0) -> str:
    """e.g. ``nodal_result_name(12, 1, 'TFU') -> '00012TR_X.TFU'``."""
    w = node_width(max(max_node, node))
    return f"{node:0{w}d}{DOF_TAGS[dof]}.{ext}"


def nodal_rs_name(node: int, dof: int, idamp: int, max_node: int = 0) -> str:
    """e.g. ``nodal_rs_name(12, 1, 1) -> '00012TR_X01.RS'`` (idamp is 1-based)."""
    w = node_width(max(max_node, node))
    return f"{node:0{w}d}{DOF_TAGS[dof]}{idamp:02d}.RS"


def element_result_name(etype: str, group: int, elem: int, comp: str, ext: str) -> str:
    """e.g. ``element_result_name('BEAMS', 3, 45, 'MXJ', 'THS') -> 'BEAMS_003_00045_MXJ.THS'``."""
    return f"{etype}_{group:03d}_{elem:05d}_{comp}.{ext}"


def layer_th_name(prefix: str, layer: int) -> str:
    """SOIL layer histories: ``ACC001.TH``, ``SN001.TH``, ``SS001.TH``."""
    return f"{prefix}{layer:03d}.TH"


def incoherent_file8_name(sim: int, direction: int) -> str:
    """FILE8 name for incoherent simulation ``sim`` (1..50) and direction 1 X, 2 Y, 3 Z."""
    return f"FILE8{3 * (sim - 1) + direction:03d}"


def loadcase_name(base: str, case: int) -> str:
    """``FILE8001``, ``FILE9012``, ``FILE77003`` ..."""
    return f"{base}{case:03d}"


def restart_names(order: int):
    """(COOXqqq, COOTKqqq) for 1-based frequency order number ``order``."""
    return f"COOX{order:03d}", f"COOTK{order:03d}"


def deg2rad(a: float) -> float:
    return a * math.pi / 180.0
