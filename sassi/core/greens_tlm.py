"""Exact thin-layer Green functions of a layered medium (Kausel 1981): the reference oracle.

For a *discrete* layered medium (the same thin-layer column as SITE) the response to a point or
line load is an exact, closed-form sum over the Rayleigh and Love modes (R1 §3.2, §3.5;
[S4] Eq. 47-48, 88-89).  POINT approximates the point load by a finite central zone (an FE core of
radius R0, R1 §3.3); comparing POINT's far field with these series isolates the core
approximation from the layer discretisation (VP-09, R1 V8).  They are also used by the unit tests
and by VP-07 (static limits) as an independent check.

Conventions: ``exp(+i w t)``, outgoing ``H^(2)``, z up, interfaces numbered from the surface.
Cylindrical Fourier decomposition (R1 §3.2):

* horizontal unit load ``P_x = 1`` (harmonic mu = 1): ``u_rho = u~ cos(theta)``,
  ``u_theta = -v~ sin(theta)``, ``u_z = w~ cos(theta)``, theta measured from the load direction;
* vertical unit load ``P_z = +1`` (up, mu = 0): ``u_rho = p~``, ``u_z = q~`` (z up).

With the Kausel normalisation of the modes (sassi.core.tlm) these are true physical displacements.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
from scipy.special import hankel2

from .tlm import Modes, quiet_fpe


def _dh1(x):
    """d/dx H_1^(2)(x) = H_0 - H_1/x."""
    return hankel2(0, x) - hankel2(1, x) / x


@quiet_fpe
def point_load(modes: Modes, m, n, r) -> Dict[str, np.ndarray]:
    """Kausel point-load Green functions ([S4] Eq. 88-89, R1 §3.2).

    ``m`` = observation free-interface index (or array), ``n`` = load free-interface index (or
    array), ``r`` = horizontal distance(s) > 0.  Arrays broadcast against each other.
    Returns a dict with the cylindrical components

    * ``u``, ``v``, ``w``: horizontal unit load (u~_rho, u~_theta, u~_z)
    * ``p``, ``q``: vertical (upward) unit load (u_rho, u_z)
    """
    kR, phx, phz, kL, phy = modes.kR, modes.phix, modes.phiz, modes.kL, modes.phiy
    m = np.asarray(m)
    n = np.asarray(n)
    r = np.asarray(r, dtype=float)
    shape = np.broadcast(m, n, r).shape
    m_, n_, r_ = (np.broadcast_to(a, shape).ravel() for a in (m, n, r))
    xR = kR[None, :] * r_[:, None]          # (N, 2Nf)
    xL = kL[None, :] * r_[:, None]          # (N, Nf)
    H0R, H1R = hankel2(0, xR), hankel2(1, xR)
    H1L = hankel2(1, xL)
    dH1R = H0R - H1R / xR
    dH1L = hankel2(0, xL) - H1L / xL
    cxx = phx[m_, :] * phx[n_, :]            # Rayleigh x-x products
    czx = phz[m_, :] * phx[n_, :]
    cxz = phx[m_, :] * phz[n_, :]
    czz = phz[m_, :] * phz[n_, :]
    cyy = phy[m_, :] * phy[n_, :]
    a = 1.0 / 4j
    rr = r_
    u = a * (np.sum(cxx * dH1R, 1) + np.sum(cyy * H1L / kL[None, :], 1) / rr)
    v = a * (np.sum(cxx * H1R / kR[None, :], 1) / rr + np.sum(cyy * dH1L, 1))
    w = -a * np.sum(czx * H1R, 1)
    p = a * np.sum(cxz * H1R, 1)
    q = a * np.sum(czz * H0R, 1)
    return {key: val.reshape(shape) for key, val in (("u", u), ("v", v), ("w", w), ("p", p), ("q", q))}


@quiet_fpe
def line_load(modes: Modes, m, n, x) -> Dict[str, np.ndarray]:
    """Kausel line-load Green functions of the plane-strain problem ([S4] Eq. 47-48, R1 §3.5).

    Physical displacements at interface ``m`` and abscissa ``x`` due to unit line loads at
    interface ``n``, ``x = 0``:

    * ``xx``: u_x due to P_x; ``zx``: u_z due to P_x (odd in x); ``xz``: u_x due to P_z (odd);
      ``zz``: u_z due to P_z; ``yy``: u_y due to P_y (SH).
    """
    kR, phx, phz, kL, phy = modes.kR, modes.phix, modes.phiz, modes.kL, modes.phiy
    m = np.asarray(m)
    n = np.asarray(n)
    x = np.asarray(x, dtype=float)
    shape = np.broadcast(m, n, x).shape
    m_, n_, x_ = (np.broadcast_to(a, shape).ravel() for a in (m, n, x))
    sgn = np.where(x_ >= 0, 1.0, -1.0)
    ER = np.exp(-1j * kR[None, :] * np.abs(x_)[:, None]) / kR[None, :]
    EL = np.exp(-1j * kL[None, :] * np.abs(x_)[:, None]) / kL[None, :]
    a = 1.0 / 2j
    xx = a * np.sum(phx[m_] * phx[n_] * ER, 1)
    # i-scaled quantities: u~_z = i u_z, p~_z = i p_z (R1 §3.5)
    zx = -1j * sgn * a * np.sum(phz[m_] * phx[n_] * ER, 1)
    xz = 1j * sgn * a * np.sum(phx[m_] * phz[n_] * ER, 1)
    zz = a * np.sum(phz[m_] * phz[n_] * ER, 1)
    yy = a * np.sum(phy[m_] * phy[n_] * EL, 1)
    return {key: val.reshape(shape) for key, val in (("xx", xx), ("zx", zx), ("xz", xz), ("zz", zz), ("yy", yy))}


def _I3_at_R(k: np.ndarray, R: float) -> np.ndarray:
    """``I_3 = int_0^inf J_1(k' R)^2 / (k' (k'^2 - k^2)) dk'`` (R1 §3.2 rho = R branch):
    ``pi/(2 i k^2) J_1(kR) H_1(kR) - 1/(2 k^2)``."""
    from scipy.special import jv
    return np.pi / (2j * k ** 2) * jv(1, k * R) * hankel2(1, k * R) - 1.0 / (2.0 * k ** 2)


@quiet_fpe
def disk_load_average(modes: Modes, m: int, n: int, R: float) -> Dict[str, complex]:
    """Average displacement under a uniform disk load of radius R and unit resultant at interface n,
    observed on the disk area at interface m ([S4] Eq. 65, 72; R1 §3.2):

    ``u_z,avg = 2 q sum phi_z phi_z I3R``, ``u_x,avg = q [sum phi_x phi_x I3R + sum phi_y phi_y I3L]``
    with ``q = 1/(pi R^2)``.  Returns ``{'vertical': ..., 'horizontal': ...}``.
    """
    q = 1.0 / (np.pi * R * R)
    I3R = _I3_at_R(modes.kR, R)
    I3L = _I3_at_R(modes.kL, R)
    uz = 2.0 * q * np.sum(modes.phiz[m] * modes.phiz[n] * I3R)
    ux = q * (np.sum(modes.phix[m] * modes.phix[n] * I3R) + np.sum(modes.phiy[m] * modes.phiy[n] * I3L))
    return {"vertical": complex(uz), "horizontal": complex(ux)}
