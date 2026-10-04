"""POINT3: point-load solutions of the layered soil by an axisymmetric central zone (R1 §3.3).

How POINT obtains the soil flexibility
--------------------------------------
Each column of the free-field flexibility matrix is the response of the layered soil to a unit
harmonic point load at one interaction node ([S3], R1 §3.1).  A true point load has an infinite
displacement under the load, so SASSI spreads it over a small *central zone* of radius R0:

* **Core** (R1 §3.3 i, D-PNT-03): one radial axisymmetric finite element (nodes on the axis and
  at rho = R0) per sublayer of the SITE column, bilinear in (rho, z), mixed mass in z, dashpots at
  the base when the column has a viscous base.  The load is expanded in Fourier harmonics around
  the vertical axis: a vertical load excites mu = 0, a horizontal load mu = 1.  Axis constraints:
  mu = 0 -> u_rho = u_theta = 0; mu = 1 -> u_rho = u_theta (one unknown), u_z = 0.
* **Exterior** (R1 §3.3 ii): for rho >= R0 the field is a superposition of outgoing Rayleigh and
  Love modes (from SITE Mode 1), ``u(rho) = Psi_mu(rho) alpha`` with Hankel functions H^(2).
* **Transmitting boundary** (R1 §3.3 iii): the consistent boundary stiffness of the exterior,
  ``R_mu = -T_mu Psi_mu(R0)^-1``, is added to the core at rho = R0.
* **Solve** (R1 §3.3 iv): ``[K_c - w^2 M_c + R_mu] U = F`` with ``F = P/pi`` (mu = 1, on the tied
  axis unknown) or ``P/(2 pi)`` (mu = 0); then ``alpha = Psi_mu(R0)^-1 U_boundary``.

Numerical scaling: the exterior basis is stored *scaled to R0*, i.e. column j of Psi is
``H_mu(k_j rho) exp(i k_j R0)`` (and the matching derivatives), computed with the exponentially
scaled Hankel function ``hankel2e``.  For rho >= R0 this never overflows (|exp(-i k (rho-R0))| <= 1
because Im k <= 0) and evanescent modes do not underflow at rho = R0.  The boundary stiffness R
is independent of this column scaling; the stored amplitudes ``alpha`` refer to the scaled basis.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
from scipy import linalg as sla
from scipy.special import hankel2e

from .tlm import (BASE_DASHPOT, MASS_MIXED, Column, Modes, boundary_integrals, mass_pattern,
                  quiet_fpe)


# ---------------------------------------------------------------------------------------
# Hankel functions scaled to R0
# ---------------------------------------------------------------------------------------
def scaled_hankel(mu: int, k: np.ndarray, rho: np.ndarray, R0: float):
    """``H_mu(k rho) e^{i k R0}``, its first and second derivatives with respect to the argument.

    ``k`` (nk,), ``rho`` (nr,) -> arrays (nr, nk).  Uses ``H_mu' = (H_{mu-1} - H_{mu+1})/2`` and
    ``H_mu'' = -H_mu'/x - (1 - mu^2/x^2) H_mu`` (R1 §3.3 ii).
    """
    k = np.asarray(k, complex)[None, :]
    rho = np.asarray(rho, float).reshape(-1)[:, None]
    x = k * rho
    E = np.exp(-1j * k * (rho - R0))
    Hm1 = hankel2e(mu - 1, x) * E
    H = hankel2e(mu, x) * E
    Hp1 = hankel2e(mu + 1, x) * E
    dH = 0.5 * (Hm1 - Hp1)
    ddH = -dH / x - (1.0 - mu * mu / (x * x)) * H
    return x, H, dH, ddH


@quiet_fpe
def psi_matrices(modes: Modes, mu: int, rho: float, R0: float, rows=None) -> Tuple[np.ndarray, np.ndarray]:
    """Exterior mode matrix ``Psi_mu(rho)`` and its radial derivative (R1 §3.3 ii), scaled to R0.

    Rows: ``[u~_rho at the interfaces; u~_theta; u~_z]``; columns: 2Nf Rayleigh then Nf Love modes.
    ``rows`` optionally selects interfaces (default all free interfaces).
    """
    phx, phz, phy = modes.phix, modes.phiz, modes.phiy
    if rows is not None:
        phx, phz, phy = phx[rows], phz[rows], phy[rows]
    kR, kL = modes.kR, modes.kL
    xR, HR, dHR, ddHR = (a[0] for a in scaled_hankel(mu, kR, [rho], R0))
    xL, HL, dHL, ddHL = (a[0] for a in scaled_hankel(mu, kL, [rho], R0))
    zL = np.zeros_like(phy)
    Psi = np.block([[phx * dHR, phy * (mu * HL / xL)],
                    [phx * (mu * HR / xR), phy * dHL],
                    [-phz * HR, zL]])
    dPsi = np.block([[phx * (kR * ddHR), phy * (mu * kL * (dHL / xL - HL / xL ** 2))],
                     [phx * (mu * kR * (dHR / xR - HR / xR ** 2)), phy * (kL * ddHL)],
                     [-phz * (kR * dHR), zL]])
    return Psi, dPsi


@quiet_fpe
def boundary_stiffness(col: Column, modes: Modes, mu: int, R0: float) -> Tuple[np.ndarray, np.ndarray]:
    """Consistent cylindrical transmitting boundary ``R_mu = -T_mu Psi_mu(R0)^-1`` (R1 §3.3 iii).

    ``T_mu`` holds, per mode column, the consistent generalised forces (per radian) that the
    exterior exerts on the core at rho = R0:

        f_rho   = R0 [E_lp2 d(u_rho)/drho + E_lam (u_rho - mu u_theta)/R0 + Q_lam u_z]
        f_theta = R0 [E_G (mu u_rho/R0 + d(u_theta)/drho - u_theta/R0)]
        f_z     = R0 [Q_G u_rho + E_G d(u_z)/drho]

    Returns ``(R, Psi(R0))`` (3Nf x 3Nf each).
    """
    n = modes.n_free
    Elp, _ = boundary_integrals(col, "lp2")
    El, Ql = boundary_integrals(col, "lam")
    EG, QG = boundary_integrals(col, "G")
    P, dP = psi_matrices(modes, mu, R0, R0)
    ur, ut, uz = P[:n], P[n:2 * n], P[2 * n:]
    dur, dut, duz = dP[:n], dP[n:2 * n], dP[2 * n:]
    fr = R0 * (Elp @ dur + El @ (ur - mu * ut) / R0 + Ql @ uz)
    ft = R0 * (EG @ (mu * ur / R0 + dut - ut / R0))
    fz = R0 * (QG @ ur + EG @ duz)
    T = np.vstack([fr, ft, fz])
    R = -sla.solve(P.T, T.T, check_finite=False).T          # -T P^-1
    return R, P


# ---------------------------------------------------------------------------------------
# Axisymmetric core element (one radial element per sublayer)
# ---------------------------------------------------------------------------------------
_GP, _GW = np.polynomial.legendre.leggauss(3)


def _core_templates(mu: int, R0: float):
    """Layer-independent parts of the core element stiffness (R1 §3.3 i).

    The 12 local unknowns are ordered ``(a, b, c)`` = (radial node 0 axis / 1 rim, vertical node
    0 top / 1 bottom, component rho/theta/z), index ``(2a + b)*3 + c``.  The strain-displacement
    matrix splits into ``B = B0 + B1/h`` (B1 = the d/dz terms), and ``D = lam D_lam + G D_G``, so
    ``K_e = sum_c c (h K0_c + K1_c + K2_c/h)`` for c in (lam, G).  3x3 Gauss points: rho = 0 is
    never evaluated, and after the axis constraints the integrand is a polynomial (exact).
    """
    Dlam = np.zeros((6, 6))
    Dlam[:3, :3] = 1.0
    DG = np.diag([2.0, 2.0, 2.0, 1.0, 1.0, 1.0])
    out = {key: np.zeros((12, 12)) for key in ("l0", "l1", "l2", "g0", "g1", "g2")}
    for gr, wr in zip(_GP, _GW):
        rho = 0.5 * R0 * (gr + 1.0)
        L = np.array([1.0 - rho / R0, rho / R0])
        dL = np.array([-1.0 / R0, 1.0 / R0])
        for gz, wz in zip(_GP, _GW):
            s = 0.5 * (gz + 1.0)                  # 0 at the top node, 1 at the bottom node
            Nz = np.array([1.0 - s, s])
            dNz = np.array([1.0, -1.0])           # times 1/h: d/dz with z up
            B0 = np.zeros((6, 12))
            B1 = np.zeros((6, 12))
            for a in range(2):
                for b in range(2):
                    base = (2 * a + b) * 3
                    ir, it, iz = base, base + 1, base + 2
                    N = L[a] * Nz[b]
                    Nr = dL[a] * Nz[b]
                    Nzz = L[a] * dNz[b]
                    # strains [e_rr, e_tt, e_zz, g_rz, g_rt, g_tz] (R1 §3.3 i)
                    B0[0, ir] += Nr
                    B0[1, ir] += N / rho
                    B0[1, it] += -mu * N / rho
                    B1[2, iz] += Nzz
                    B1[3, ir] += Nzz
                    B0[3, iz] += Nr
                    B0[4, ir] += mu * N / rho
                    B0[4, it] += Nr - N / rho
                    B1[5, it] += Nzz
                    B0[5, iz] += mu * N / rho
            w = wr * wz * 0.5 * R0 * 0.5 * rho     # dz = (h/2) ds' -> the h factor is applied later
            for D, tag in ((Dlam, "l"), (DG, "g")):
                out[tag + "0"] += w * (B0.T @ D @ B0)
                out[tag + "1"] += w * (B0.T @ D @ B1 + B1.T @ D @ B0)
                out[tag + "2"] += w * (B1.T @ D @ B1)
    return out


def _radial_mass(R0: float) -> np.ndarray:
    """``int_0^R0 L_a L_b rho drho`` (consistent in rho)."""
    return R0 * R0 * np.array([[1.0 / 12.0, 1.0 / 12.0], [1.0 / 12.0, 1.0 / 4.0]])


@quiet_fpe
def core_dynamic_stiffness(col: Column, mu: int, R0: float, omega: float, mass: str = MASS_MIXED) -> np.ndarray:
    """Per-radian dynamic stiffness ``K_c - w^2 M_c (+ i w C_base)`` of the core on all free
    interfaces, full unknown order ``(a*3 + c)*Nf + i`` (a radial node, c component, i interface).

    Mass: mixed (1/2 lumped + 1/2 consistent) in z exactly as the SITE layers, consistent in rho
    (R1 §3.3 i "same mass lumping in z as in SITE").  A viscous base adds
    ``i w c int L_a L_b rho drho`` on the bottom-interface unknowns (c_s on rho, theta; c_p on z).
    """
    t = _core_templates(mu, R0)
    h = col.h[:, None, None]
    lam = col.lam[:, None, None]
    G = col.G[:, None, None]
    Ke = (lam * (h * t["l0"] + t["l1"] + t["l2"] / h) + G * (h * t["g0"] + t["g1"] + t["g2"] / h)).astype(complex)
    Me = np.kron(np.kron(_radial_mass(R0), mass_pattern(mass)), np.eye(3))
    Ke = Ke - omega ** 2 * (col.rho[:, None, None] * h) * Me[None, :, :]
    ni = col.n_iface
    # global index of each local unknown, per sublayer j
    a = np.repeat(np.arange(2), 6)
    b = np.tile(np.repeat(np.arange(2), 3), 2)
    c = np.tile(np.arange(3), 4)
    j = np.arange(col.n_sub)[:, None]
    gidx = (a * 3 + c)[None, :] * ni + (j + b[None, :])          # (nsub, 12)
    K = np.zeros((6 * ni, 6 * ni), complex)
    rows = np.broadcast_to(gidx[:, :, None], Ke.shape)
    cols = np.broadcast_to(gidx[:, None, :], Ke.shape)
    np.add.at(K, (rows.ravel(), cols.ravel()), Ke.ravel())
    if col.base == BASE_DASHPOT:
        Rm = _radial_mass(R0)
        i = ni - 1
        for comp, cc in ((0, col.cs), (1, col.cs), (2, col.cp)):
            for aa in range(2):
                for bb in range(2):
                    K[(aa * 3 + comp) * ni + i, (bb * 3 + comp) * ni + i] += 1j * omega * cc * Rm[aa, bb]
    # keep the free interfaces only
    nf = col.n_free
    keep = (np.arange(6)[:, None] * ni + np.arange(nf)[None, :]).ravel()
    return K[np.ix_(keep, keep)]


def _reduction(mu: int, nf: int) -> np.ndarray:
    """Boolean/tie matrix T (6Nf x 4Nf) imposing the axis constraints (R1 §3.3 i).

    Reduced unknowns: ``[axis (Nf); rim u_rho (Nf); rim u_theta (Nf); rim u_z (Nf)]`` where the axis
    unknown is u_z for mu = 0 and the tied u_rho = u_theta for mu = 1.
    """
    T = np.zeros((6 * nf, 4 * nf))
    i = np.arange(nf)
    if mu == 0:
        T[2 * nf + i, i] = 1.0
    else:
        T[0 * nf + i, i] = 1.0
        T[1 * nf + i, i] = 1.0
    for comp in range(3):
        T[(3 + comp) * nf + i, (1 + comp) * nf + i] = 1.0
    return T


@dataclass
class Point3Solution:
    """POINT3 results at one frequency for loads at the free interfaces ``load_idx``.

    ``alpha[mu]`` (nL, 3Nf): exterior amplitudes (basis scaled to R0) for a unit load at each
    load interface; ``axis[mu]`` (nL, nL, 3): core-axis displacements ``(u~_rho, u~_theta, u~_z)``
    at the load interfaces (index ``[load, obs, comp]``); ``rim[mu]`` (nL, 3Nf): rim unknowns.
    """

    alpha: Dict[int, np.ndarray]
    axis: Dict[int, np.ndarray]
    rim: Dict[int, np.ndarray]


@quiet_fpe
def solve_point3(col: Column, modes: Modes, omega: float, R0: float, load_idx, mass: str = MASS_MIXED
                 ) -> Point3Solution:
    """Unit point loads at the axis of the core, at every free interface in ``load_idx``
    (R1 §3.3 iv): vertical (upward) load -> mu = 0, horizontal (x) load -> mu = 1."""
    load_idx = np.asarray(load_idx, int)
    nf = col.n_free
    if modes.n_free != nf:
        raise ValueError("modes and column have different numbers of free interfaces")
    nL = len(load_idx)
    alpha, axis, rim = {}, {}, {}
    for mu in (0, 1):
        T = _reduction(mu, nf)
        K = T.T @ core_dynamic_stiffness(col, mu, R0, omega, mass) @ T
        R, P = boundary_stiffness(col, modes, mu, R0)
        K[nf:, nf:] += R
        F = np.zeros((4 * nf, nL), complex)
        F[load_idx, np.arange(nL)] = 1.0 / (np.pi if mu == 1 else 2.0 * np.pi)
        U = sla.lu_solve(sla.lu_factor(K, check_finite=False), F, check_finite=False)
        if not np.all(np.isfinite(U)):
            raise FloatingPointError("POINT3 core solution is not finite (singular core system)")
        Ub = U[nf:, :]
        a = np.linalg.solve(P, Ub)                      # (3Nf, nL)
        ax = np.zeros((nL, nL, 3), complex)
        uax = U[:nf, :][load_idx, :].T                  # [load, obs]
        if mu == 0:
            ax[:, :, 2] = uax
        else:
            ax[:, :, 0] = uax
            ax[:, :, 1] = uax
        alpha[mu], axis[mu], rim[mu] = a.T, ax, Ub.T
    return Point3Solution(alpha, axis, rim)


# ---------------------------------------------------------------------------------------
# Exterior field evaluation (used by flexibility.py and the verification problems)
# ---------------------------------------------------------------------------------------
@quiet_fpe
def exterior_field(kR, kL, phix_obs, phiz_obs, phiy_obs, alpha, mu: int, r, R0: float) -> np.ndarray:
    """Cylindrical displacement amplitudes ``(u~_rho, u~_theta, u~_z)`` for r >= R0 (R1 §3.6).

    ``phi*_obs`` hold the mode rows at the observation interfaces (nO, ...), ``alpha`` (nL, 3Nf)
    the amplitudes of each load (basis scaled to R0).  Returns ``(nr, nL, nO, 3)``.
    """
    kR = np.asarray(kR, complex)
    kL = np.asarray(kL, complex)
    r = np.asarray(r, float).reshape(-1)
    alpha = np.asarray(alpha, complex)
    nO = phix_obs.shape[0]
    nL = alpha.shape[0]
    nR = kR.size
    a = alpha[:, :nR]                   # Rayleigh amplitudes (nL, 2Nf)
    b = alpha[:, nR:]                   # Love amplitudes (nL, Nf)
    xR, HR, dHR, _ = scaled_hankel(mu, kR, r, R0)
    xL, HL, dHL, _ = scaled_hankel(mu, kL, r, R0)
    # coefficient tensors c[n, m, j] = amplitude(load n, mode j) * mode value(obs m, mode j)
    cx = (a[:, None, :] * phix_obs[None, :, :]).reshape(nL * nO, -1).T     # (2Nf, nL*nO)
    cz = (a[:, None, :] * phiz_obs[None, :, :]).reshape(nL * nO, -1).T
    cy = (b[:, None, :] * phiy_obs[None, :, :]).reshape(nL * nO, -1).T     # (Nf, nL*nO)
    out = np.empty((r.size, nL * nO, 3), complex)
    if mu == 0:
        out[:, :, 0] = dHR @ cx
        out[:, :, 1] = dHL @ cy
    else:
        out[:, :, 0] = dHR @ cx + (mu * HL / xL) @ cy
        out[:, :, 1] = (mu * HR / xR) @ cx + dHL @ cy
    out[:, :, 2] = -(HR @ cz)
    return out.reshape(r.size, nL, nO, 3)
