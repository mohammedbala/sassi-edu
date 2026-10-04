"""POINT2 (tier P1): line-load solutions of the layered soil for plane-strain SSI (R1 §3.5, §2.6).

The 2D counterpart of POINT3.  The central zone is a strip |x| <= R0 of two plane-strain
elements per SITE sublayer (nodes at x = -R0, 0, +R0), linear in z, mixed mass in z.  At
x = +-R0 the semi-infinite layered regions are represented by Waas-Lysmer consistent transmitting
boundaries built from the Rayleigh (P-SV) and Love (SH) modes of FILE2::

    outgoing mode j at the right (s = +1) / left (s = -1) boundary:
        u_x = phi_xj exp(-i s k_j (x - x_b)),  u_z = -i s phi_zj exp(...)    (P-SV)
        u_y = phi_yj exp(-i s k_j (x - x_b))                                (SH)
    R = -T Psi^-1 (consistent boundary forces T on the core, R1 §2.6)

Unit line loads (per unit length out of plane) at the centre node of interface n give the far
field ``U(x) = sum_j alpha_j psi_j exp(-i k_j (|x| - R0))`` for |x| >= R0.  By symmetry the
horizontal displacement is even and the vertical odd in x for a horizontal load, and vice versa
for a vertical load (requirements §4.3).  The exact reference is Kausel's line-load Green
function (:func:`sassi.core.greens_tlm.line_load`, R1 V7).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
from scipy import linalg as sla

from .tlm import BASE_DASHPOT, MASS_MIXED, Column, Modes, boundary_integrals, mass_pattern, quiet_fpe

_GP, _GW = np.polynomial.legendre.leggauss(2)


def _plane_templates(R0: float):
    """Layer-independent parts of the 4-node plane-strain element of width R0 (P-SV).

    Local unknowns ``(a, b, c)``: a = left/right node in x, b = top/bottom node in z, c = x/z
    component, index ``(2a + b)*2 + c``.  ``K_e = sum_c c (h K0_c + K1_c + K2_c/h)``.
    """
    Dlam = np.zeros((3, 3))
    Dlam[:2, :2] = 1.0
    DG = np.diag([2.0, 2.0, 1.0])
    out = {key: np.zeros((8, 8)) for key in ("l0", "l1", "l2", "g0", "g1", "g2")}
    dL = np.array([-1.0, 1.0]) / R0
    for gx, wx in zip(_GP, _GW):
        xi = 0.5 * (gx + 1.0)
        L = np.array([1.0 - xi, xi])
        for gz, wz in zip(_GP, _GW):
            s = 0.5 * (gz + 1.0)
            Nz = np.array([1.0 - s, s])
            dNz = np.array([1.0, -1.0])           # times 1/h (z up)
            B0 = np.zeros((3, 8))
            B1 = np.zeros((3, 8))
            for a in range(2):
                for b in range(2):
                    iu = (2 * a + b) * 2
                    iw = iu + 1
                    B0[0, iu] += dL[a] * Nz[b]                 # e_xx
                    B1[1, iw] += L[a] * dNz[b]                 # e_zz
                    B1[2, iu] += L[a] * dNz[b]                 # g_xz = du/dz + dw/dx
                    B0[2, iw] += dL[a] * Nz[b]
            w = wx * wz * 0.5 * R0 * 0.5
            for D, tag in ((Dlam, "l"), (DG, "g")):
                out[tag + "0"] += w * (B0.T @ D @ B0)
                out[tag + "1"] += w * (B0.T @ D @ B1 + B1.T @ D @ B0)
                out[tag + "2"] += w * (B1.T @ D @ B1)
    return out


def _antiplane_templates(R0: float):
    """Anti-plane (SH) 4-node element of width R0: ``K_e = G (h K0 + K2/h)``."""
    K0 = np.zeros((4, 4))
    K2 = np.zeros((4, 4))
    dL = np.array([-1.0, 1.0]) / R0
    for gx, wx in zip(_GP, _GW):
        xi = 0.5 * (gx + 1.0)
        L = np.array([1.0 - xi, xi])
        for gz, wz in zip(_GP, _GW):
            s = 0.5 * (gz + 1.0)
            Nz = np.array([1.0 - s, s])
            dNz = np.array([1.0, -1.0])
            bx = np.array([dL[a] * Nz[b] for a in range(2) for b in range(2)])
            bz = np.array([L[a] * dNz[b] for a in range(2) for b in range(2)])
            w = wx * wz * 0.5 * R0 * 0.5
            K0 += w * np.outer(bx, bx)
            K2 += w * np.outer(bz, bz)
    return K0, K2


def _x_mass(R0: float) -> np.ndarray:
    return R0 / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])


def _assemble_strip(col: Column, Ke: np.ndarray, ncomp: int) -> np.ndarray:
    """Assemble element matrices of both strip elements (x in [-R0, 0] and [0, R0]).

    Global unknown ``c*(3*ni) + ix*ni + i`` (component c, x node ix = 0 left / 1 centre / 2 right,
    interface i); returns the matrix restricted to the free interfaces.
    """
    ni = col.n_iface
    nloc = 4 * ncomp
    a = np.repeat(np.arange(2), 2 * ncomp)
    b = np.tile(np.repeat(np.arange(2), ncomp), 2)
    c = np.tile(np.arange(ncomp), 4)
    K = np.zeros((ncomp * 3 * ni, ncomp * 3 * ni), complex)
    j = np.arange(col.n_sub)[:, None]
    for e in range(2):
        g = c[None, :] * (3 * ni) + (e + a)[None, :] * ni + (j + b[None, :])
        rows = np.broadcast_to(g[:, :, None], (col.n_sub, nloc, nloc))
        cols = np.broadcast_to(g[:, None, :], (col.n_sub, nloc, nloc))
        np.add.at(K, (rows.ravel(), cols.ravel()), Ke.ravel())
    nf = col.n_free
    keep = (np.arange(ncomp * 3)[:, None] * ni + np.arange(nf)[None, :]).ravel()
    return K[np.ix_(keep, keep)]


def _base_dashpot(col: Column, K: np.ndarray, omega: float, coefs, R0: float) -> None:
    """Add ``i w c int L_a L_b dx`` of the base dashpots on the bottom free interface."""
    if col.base != BASE_DASHPOT:
        return
    nf = col.n_free
    Xm = _x_mass(R0)
    i = nf - 1
    for comp, cc in enumerate(coefs):
        for e in range(2):
            for aa in range(2):
                for bb in range(2):
                    K[comp * 3 * nf + (e + aa) * nf + i, comp * 3 * nf + (e + bb) * nf + i] += 1j * omega * cc * Xm[aa, bb]


@quiet_fpe
def core_psv(col: Column, R0: float, omega: float, mass: str = MASS_MIXED) -> np.ndarray:
    """Dynamic stiffness of the P-SV strip core on the free interfaces."""
    t = _plane_templates(R0)
    h = col.h[:, None, None]
    lam, G = col.lam[:, None, None], col.G[:, None, None]
    Ke = (lam * (h * t["l0"] + t["l1"] + t["l2"] / h) + G * (h * t["g0"] + t["g1"] + t["g2"] / h)).astype(complex)
    Me = np.kron(np.kron(_x_mass(R0), mass_pattern(mass)), np.eye(2))
    Ke = Ke - omega ** 2 * (col.rho[:, None, None] * h) * Me[None]
    K = _assemble_strip(col, Ke, 2)
    _base_dashpot(col, K, omega, (col.cs, col.cp), R0)
    return K


@quiet_fpe
def core_sh(col: Column, R0: float, omega: float, mass: str = MASS_MIXED) -> np.ndarray:
    """Dynamic stiffness of the anti-plane (SH) strip core on the free interfaces."""
    K0, K2 = _antiplane_templates(R0)
    h = col.h[:, None, None]
    G = col.G[:, None, None]
    Ke = (G * (h * K0 + K2 / h)).astype(complex)
    Me = np.kron(_x_mass(R0), mass_pattern(mass))
    Ke = Ke - omega ** 2 * (col.rho[:, None, None] * h) * Me[None]
    K = _assemble_strip(col, Ke, 1)
    _base_dashpot(col, K, omega, (col.cs,), R0)
    return K


@quiet_fpe
def boundary_psv(col: Column, modes: Modes, side: int):
    """Waas-Lysmer P-SV boundary ``R = -T Psi^-1`` on ``[u_x(ifaces); u_z(ifaces)]`` (R1 §2.6)."""
    Elp, _ = boundary_integrals(col, "lp2")
    _, Ql = boundary_integrals(col, "lam")
    EG, QG = boundary_integrals(col, "G")
    s = float(side)
    k = modes.kR
    ux = modes.phix.astype(complex)
    uz = -1j * s * modes.phiz
    dux = -1j * s * k * ux
    duz = -1j * s * k * uz
    T = s * np.vstack([Elp @ dux + Ql @ uz, QG @ ux + EG @ duz])
    Psi = np.vstack([ux, uz])
    R = -sla.solve(Psi.T, T.T, check_finite=False).T
    return R, Psi


@quiet_fpe
def boundary_sh(col: Column, modes: Modes):
    """Waas-Lysmer SH boundary ``R = i E_G Phi_y K_L Phi_y^-1`` (same on both sides)."""
    EG, _ = boundary_integrals(col, "G")
    Psi = modes.phiy.astype(complex)
    T = -1j * EG @ (Psi * modes.kL[None, :])
    R = -sla.solve(Psi.T, T.T, check_finite=False).T
    return R, Psi


@dataclass
class Point2Solution:
    """POINT2 results: right-boundary amplitudes ``alpha_x``, ``alpha_z`` (nL, 2Nf) for x/z line
    loads, ``alpha_y`` (nL, Nf) for y line loads, and centre-node displacements ``centre[...]``
    (nL load, nL obs) for the components 'xx', 'zz', 'yy'."""

    alpha_x: np.ndarray
    alpha_z: np.ndarray
    alpha_y: np.ndarray
    centre: Dict[str, np.ndarray]


@quiet_fpe
def solve_point2(col: Column, modes: Modes, omega: float, R0: float, load_idx, mass: str = MASS_MIXED
                 ) -> Point2Solution:
    """Unit line loads at the centre node (x = 0) of every free interface in ``load_idx``."""
    load_idx = np.asarray(load_idx, int)
    nf = col.n_free
    nL = load_idx.size
    # P-SV
    K = core_psv(col, R0, omega, mass)
    def bidx(ix):
        return np.concatenate([0 * 3 * nf + ix * nf + np.arange(nf), 1 * 3 * nf + ix * nf + np.arange(nf)])
    Rr, Pr = boundary_psv(col, modes, +1)
    Rl, _ = boundary_psv(col, modes, -1)
    K[np.ix_(bidx(2), bidx(2))] += Rr
    K[np.ix_(bidx(0), bidx(0))] += Rl
    F = np.zeros((K.shape[0], 2 * nL), complex)
    F[0 * 3 * nf + nf + load_idx, np.arange(nL)] = 1.0           # x loads at the centre node
    F[1 * 3 * nf + nf + load_idx, nL + np.arange(nL)] = 1.0      # z loads
    U = sla.lu_solve(sla.lu_factor(K, check_finite=False), F, check_finite=False)
    a = np.linalg.solve(Pr, U[bidx(2)])                          # (2Nf, 2nL)
    cx = U[0 * 3 * nf + nf + load_idx]                           # u_x at centre (obs, load)
    cz = U[1 * 3 * nf + nf + load_idx]
    # SH
    Ks = core_sh(col, R0, omega, mass)
    Rs, Ps = boundary_sh(col, modes)
    Ks[np.ix_(2 * nf + np.arange(nf), 2 * nf + np.arange(nf))] += Rs
    Ks[np.ix_(np.arange(nf), np.arange(nf))] += Rs
    Fs = np.zeros((Ks.shape[0], nL), complex)
    Fs[nf + load_idx, np.arange(nL)] = 1.0
    Us = sla.lu_solve(sla.lu_factor(Ks, check_finite=False), Fs, check_finite=False)
    ay = np.linalg.solve(Ps, Us[2 * nf + np.arange(nf)])
    cy = Us[nf + load_idx]
    if not (np.all(np.isfinite(U)) and np.all(np.isfinite(Us))):
        raise FloatingPointError("POINT2 core solution is not finite")
    centre = {"xx": cx[:, :nL].T, "zz": cz[:, nL:].T, "yy": cy.T,
              "zx": cz[:, :nL].T, "xz": cx[:, nL:].T}            # [load, obs]; zx, xz ~ 0 by symmetry
    return Point2Solution(a[:, :nL].T, a[:, nL:].T, ay.T, centre)


@quiet_fpe
def exterior_2d(kR, kL, phix_obs, phiz_obs, phiy_obs, alpha_x, alpha_z, alpha_y, a, R0: float
                ) -> Dict[str, np.ndarray]:
    """Far field for x = a >= R0 (right side): arrays (na, nL load, nO obs) for 'xx' (u_x due to
    P_x), 'zx' (u_z due to P_x), 'xz' (u_x due to P_z), 'zz', 'yy'."""
    a = np.asarray(a, float).reshape(-1)
    ER = np.exp(-1j * np.asarray(kR)[None, :] * (a[:, None] - R0))    # (na, 2Nf)
    EL = np.exp(-1j * np.asarray(kL)[None, :] * (a[:, None] - R0))
    def comb(amp, phi, E):
        c = (amp[:, None, :] * phi[None, :, :]).reshape(amp.shape[0] * phi.shape[0], -1).T
        return (E @ c).reshape(a.size, amp.shape[0], phi.shape[0])
    return {"xx": comb(alpha_x, phix_obs, ER), "zx": comb(alpha_x, -1j * phiz_obs, ER),
            "xz": comb(alpha_z, phix_obs, ER), "zz": comb(alpha_z, -1j * phiz_obs, ER),
            "yy": comb(alpha_y, phiy_obs, EL)}


def components_2d(file3, q: int, a) -> Dict[str, np.ndarray]:
    """POINT2 components at distances ``a = |x_i - x_j| >= 0`` (right-side sign convention):
    a = 0 centre values, 0 < a < R0 linear interpolation to the boundary, a >= R0 exterior."""
    A = file3.arrays
    R0 = float(file3.meta["R0"])
    nR = A["kR"].shape[1]
    al1, al0 = A["alpha1"][q], A["alpha0"][q]
    a = np.asarray(a, float).reshape(-1)
    ext = a >= R0 * (1 - 1e-12)
    args = (A["kR"][q], A["kL"][q], A["phix_obs"][q], A["phiz_obs"][q], A["phiy_obs"][q],
            al1[:, :nR], al0[:, :nR], al1[:, nR:])
    nL = al1.shape[0]
    c1, c0 = A["axis1"][q], A["axis0"][q]
    centre = {"xx": c1[..., 0], "yy": c1[..., 1], "zx": c1[..., 2], "xz": c0[..., 0], "zz": c0[..., 2]}
    out = {key: np.zeros((a.size, nL, nL), complex) for key in centre}
    if np.any(ext):
        e = exterior_2d(*args, np.maximum(a[ext], R0), R0)
        for key in out:
            out[key][ext] = e[key]
    inner = ~ext
    if np.any(inner):
        rim = exterior_2d(*args, np.array([R0]), R0)
        t = (a[inner] / R0)[:, None, None]
        for key in out:
            out[key][inner] = centre[key][None] + t * (rim[key] - centre[key][None])
    return out


@quiet_fpe
def flexibility_matrix_2d(file3, q: int, xy, pos, symmetrize: bool = True) -> np.ndarray:
    """Plane-strain flexibility (3n x 3n, [ux, uy, uz] per node) from a POINT2 FILE3.

    P-SV block ``[[F_xx, sgn F_xz], [sgn F_zx, F_zz]]`` with ``sgn = sgn(x_i - x_j)`` and the SH
    term ``F_yy`` (R1 §4.1, 2D)."""
    x = np.asarray(xy, float)[:, 0]
    F = flexibility_block_2d(file3, q, x, pos, x, pos)
    if symmetrize:
        F += F.T
        F *= 0.5
    return F


@quiet_fpe
def flexibility_block_2d(file3, q: int, x_obs, pos_obs, x_load, pos_load) -> np.ndarray:
    """Plane-strain flexibility of observation points ``x_obs`` (interface rows ``pos_obs`` of FILE3)
    under unit line loads at ``x_load`` (rows ``pos_load``): (3 n_obs, 3 n_load), node-major
    ``[ux, uy, uz]``, not symmetrised.  P-SV terms with the odd couplings times ``sgn(x_i - x_j)``,
    SH term ``F_yy`` (R1 §4.1, 2D)."""
    xo = np.asarray(x_obs, float).reshape(-1)
    xl = np.asarray(x_load, float).reshape(-1)
    no, nl = xo.size, xl.size
    d = xo[:, None] - xl[None, :]
    a = np.abs(d)
    sgn = np.sign(d)
    scale = max(float(a.max()) if a.size else 0.0, float(file3.meta["R0"]))
    key = np.round(a / scale * 1e10).astype(np.int64)
    ukey, inv = np.unique(key, return_inverse=True)
    inv = inv.reshape(no, nl)
    comp = components_2d(file3, q, ukey * (scale * 1e-10))
    jl = np.asarray(pos_load)[None, :]
    io = np.asarray(pos_obs)[:, None]
    F = np.zeros((no, 3, nl, 3), complex)
    F[:, 0, :, 0] = comp["xx"][inv, jl, io]
    F[:, 2, :, 0] = sgn * comp["zx"][inv, jl, io]
    F[:, 0, :, 2] = sgn * comp["xz"][inv, jl, io]
    F[:, 2, :, 2] = comp["zz"][inv, jl, io]
    F[:, 1, :, 1] = comp["yy"][inv, jl, io]
    F = F.reshape(3 * no, 3 * nl)
    if not np.all(np.isfinite(F)):
        raise FloatingPointError("plane-strain flexibility contains non-finite values")
    return F
