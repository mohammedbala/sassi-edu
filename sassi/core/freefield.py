"""Free-field motion of the layered site (SITE Mode 2) and its evaluation at interaction nodes.

Binding API (docs/ARCHITECTURE.md §6.3)::

    free_field_at_nodes(file1, q, xyz, depth_iface, ang_deg, xc, yc) -> (n, 3) complex

What SITE Mode 2 computes (R1 §2.7, requirements §4.2)
------------------------------------------------------
The seismic environment is a superposition of plane *wave fields* ``w`` (body waves P, SV, SH,
vertical or inclined, and surface waves Rayleigh/Love).  For every wave and frequency SITE builds
the field at the user interfaces in the SITE axes x'y'z' (z' up, x' in the vertical plane of
propagation, y' = z' x x'), **normalised so that the motion at the control point** (top of TOPL
layer ``cl``, direction ``cm``) **is 1** (D-SIT-07, within motion).  The total motion is

    U(z, x') = sum_w r_w(f) U_w(z) exp(-i k_w (x' - x'_cp))

with participation ratios ``r_w`` interpolated linearly between Frequency 1 and 2 (frequency
numbers, D-SIT-05) and summing to 1.

* Vertical body waves (P0): k = 0, the column decouples into 1-D shear (SV x', SH y') and
  compression (P z') columns on the variable-depth half-space with base dashpots (R1 §2.7a).
* Surface waves (P1): one Rayleigh or Love mode of FILE2, chosen by ``opt`` (D-SIT-09);
  physical vertical component ``u_z = -i phi_z`` (R1 §2.7c).
* Inclined body waves (P1): apparent wavenumber ``k = w sin(theta)/V_hs``; the half-space is the
  *exact* elastic half-space stiffness at k (Kausel-Roesset form, derived here from potentials),
  the generated buffer sublayers are not used (D-SIT-08).

ANALYS evaluates FILE1 at its interaction nodes with :func:`free_field_at_nodes`: the phase uses
``x' = (x - xc) cos(ang) + (y - yc) sin(ang)`` (the control point (xc, yc) has x' = 0) and the
x'y' components are rotated to the global axes by the angle ``ang`` (R1 §2.7d).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from ..conventions import deg2rad
from ..io.container import Container
from .tlm import (MASS_MIXED, Column, ModeChoice, Modes, choose_mode, column_loss, column_matrices,
                  outcrop_response, quiet_fpe, vertical_response)

WAVE_R, WAVE_SV, WAVE_P, WAVE_SH, WAVE_L = 1, 2, 3, 4, 5
WAVE_NAMES = {WAVE_R: "R", WAVE_SV: "SV", WAVE_P: "P", WAVE_SH: "SH", WAVE_L: "L"}
WAVE_FAMILY = {WAVE_R: 0, WAVE_SV: 0, WAVE_P: 0, WAVE_SH: 1, WAVE_L: 1}   # SITE <wopt>
DIRECTION_NAMES = ("x'", "y'", "z'")


@dataclass
class WaveField:
    """Raw (un-normalised) field of one wave at the user interfaces at one frequency.

    ``U`` (nI, 3): physical components x', y', z' (z' up); ``k``: horizontal wavenumber along +x';
    ``outcrop`` (3,): discrete outcrop motion at the top of the half-space for the same incident
    wave (body waves; zeros for surface waves).
    """

    U: np.ndarray
    k: complex
    outcrop: np.ndarray


# ---------------------------------------------------------------------------------------
# Wave ratios (requirements §4.2, D-SIT-05)
# ---------------------------------------------------------------------------------------
def wave_ratio(fnum, freq1: float, freq2: float, ratio1: float, ratio2: float) -> np.ndarray:
    """Participation ratio at frequency numbers ``fnum``: linear between (freq1, ratio1) and
    (freq2, ratio2), constant outside that range."""
    n = np.asarray(fnum, float)
    if freq2 == freq1:
        return np.full(n.shape, float(ratio1))
    t = np.clip((n - freq1) / (freq2 - freq1), 0.0, 1.0)
    return ratio1 + t * (ratio2 - ratio1)


# ---------------------------------------------------------------------------------------
# Vertical body waves (P0)
# ---------------------------------------------------------------------------------------
def vertical_body_wave(col: Column, omega: float, wave_type: int, mass: str = MASS_MIXED) -> WaveField:
    """Vertically incident SV (x'), SH (y') or P (z') wave in the 1-D column (R1 §2.7a)."""
    comp = {WAVE_SV: 0, WAVE_SH: 1, WAVE_P: 2}[wave_type]
    kind = "p" if wave_type == WAVE_P else "s"
    u = vertical_response(col, omega, kind, mass)
    U = np.zeros((col.iface_user.size, 3), complex)
    U[:, comp] = u[col.iface_user]
    out = np.zeros(3, complex)
    out[comp] = outcrop_response(col, omega, kind, mass)
    return WaveField(U=U, k=0j, outcrop=out)


# ---------------------------------------------------------------------------------------
# Surface waves (P1)
# ---------------------------------------------------------------------------------------
def surface_mode_choice(col: Column, modes: Modes, wave_type: int, opt: int) -> ModeChoice:
    """Mode selection of a Rayleigh (``wave_type`` 1) or Love (5) wave field (D-SIT-09).

    Rayleigh: ``opt`` 1 shortest wavelength (largest Re k of the propagating modes), 2 least decay
    (smallest |Im k|, ties broken by the largest Re k); Love always uses the shortest-wavelength
    rule.  The propagating sector is widened by the column's material loss angle
    (:func:`sassi.core.tlm.propagating_ratio`), so heavily damped soil keeps its surface waves.
    """
    loss = column_loss(col)
    if wave_type == WAVE_R:
        return choose_mode(modes.kR, 2 if opt == 2 else 1, loss)
    if wave_type == WAVE_L:
        return choose_mode(modes.kL, 1, loss)
    raise ValueError("surface_wave: wave type must be 1 (R) or 5 (L)")


def surface_wave(col: Column, modes: Modes, wave_type: int, opt: int) -> Tuple[WaveField, int]:
    """Rayleigh (``wave_type`` 1) or Love (5) mode field (R1 §2.7c, D-SIT-09); the mode is chosen
    by :func:`surface_mode_choice`.  Returns the field and the selected mode index."""
    j = surface_mode_choice(col, modes, wave_type, opt).index
    idx = col.user_free_index()
    ok = idx >= 0
    U = np.zeros((idx.size, 3), complex)
    if wave_type == WAVE_R:
        U[ok, 0] = modes.phix[idx[ok], j]
        U[ok, 2] = -1j * modes.phiz[idx[ok], j]          # physical u_z = -i phi_z
        k = complex(modes.kR[j])
    else:
        U[ok, 1] = modes.phiy[idx[ok], j]
        k = complex(modes.kL[j])
    return WaveField(U=U, k=k, outcrop=np.zeros(3, complex)), j


# ---------------------------------------------------------------------------------------
# Inclined body waves (P1): exact half-space stiffness
# ---------------------------------------------------------------------------------------
def _vertical_wavenumber(kb2: complex, k: float) -> complex:
    """``nu = sqrt(kb^2 - k^2)`` with Im nu < 0 (Re nu > 0 when real): down-going waves
    ``exp(+i nu z)`` decay with depth (z up, z < 0 in the half-space)."""
    nu = np.sqrt(complex(kb2 - k * k))
    if abs(nu.imag) <= 1e-14 * abs(nu):
        return complex(abs(nu.real), 0.0) if nu.real < 0 else nu
    return -nu if nu.imag > 0 else nu


def halfspace_stiffness(rho: float, G: complex, M: complex, omega: float, k: float, inplane: bool = True
                        ) -> Tuple[np.ndarray, np.ndarray]:
    """Exact stiffness of a visco-elastic half-space at horizontal wavenumber ``k`` (D-SIT-08).

    Half-space below z = 0 (z up), fields ``exp(i w t - i k x)``.  With P and SV potentials,
    ``u = grad(phi) + curl(psi e_y)``, a down-going field ``phi = A e^{i nu_p z}``,
    ``psi = B e^{i nu_s z}`` has surface displacement ``[u_x, u_z] = U_d [A, B]`` and surface
    traction ``[s_xz, s_zz] = S_d [A, B]`` with

        U_d = [[-i k, -i nu_s], [i nu_p, -i k]]
        S_d = G [[2 k nu_p, ks^2 - 2k^2], [-(ks^2 - 2k^2), 2 k nu_s]],   ks^2 = rho w^2 / G*

    so ``K_out = S_d U_d^-1`` (force on the half-space surface per unit displacement of a
    down-going, radiating field).  ``K_in`` is the same for an up-going field (nu -> -nu).  Both
    are returned in Kausel variables ``{u_x, i u_z}`` (``T K T^-1``, T = diag(1, i)), the form of
    the thin-layer matrices.  Vertical incidence gives the dashpots ``i w rho V*`` (R1 §2.4).
    For ``inplane=False`` the SH values ``K_out = i nu_s G*``, ``K_in = -i nu_s G*`` are returned.
    """
    ks2 = rho * omega ** 2 / G
    nus = _vertical_wavenumber(ks2, k)
    if not inplane:
        return np.array([[1j * nus * G]]), np.array([[-1j * nus * G]])
    kp2 = rho * omega ** 2 / M
    nup = _vertical_wavenumber(kp2, k)
    a = ks2 - 2.0 * k * k
    T = np.diag([1.0, 1j])
    Ti = np.diag([1.0, -1j])
    out = []
    for sgn in (1.0, -1.0):         # +1 down-going (radiating), -1 up-going (incident)
        Ud = np.array([[-1j * k, -1j * sgn * nus], [1j * sgn * nup, -1j * k]])
        Sd = G * np.array([[2.0 * k * sgn * nup, a], [-a, 2.0 * k * sgn * nus]])
        Kp = Sd @ np.linalg.inv(Ud)
        out.append(T @ Kp @ Ti)
    return out[0], out[1]


def incident_polarisation(wave_type: int, theta: float, k: float = None, rho: float = None, G: complex = None,
                          M: complex = None, omega: float = None) -> np.ndarray:
    """Unit displacement of the incident plane wave in x'y'z' (z' up), ``theta`` from the vertical.

    Geometric (undamped) form: propagation ``(sin theta, 0, cos theta)``; P along it, SV
    ``(cos theta, 0, -sin theta)`` (= +x' at vertical incidence), SH +y'.  When the half-space
    properties and the real horizontal wavenumber ``k`` are given, the *pure* up-going P or SV wave
    of the damped medium is returned instead: ``(k, 0, nu_p)/k_p*`` and ``(nu_s, 0, -k)/k_s*``
    (complex direction cosines; with damping the geometric vector would add a small parasitic wave
    of the other type).  The two forms coincide for an elastic half-space.
    """
    if wave_type == WAVE_SH:
        return np.array([0.0, 1.0, 0.0], complex)
    if k is None:
        if wave_type == WAVE_P:
            return np.array([np.sin(theta), 0.0, np.cos(theta)], complex)
        return np.array([np.cos(theta), 0.0, -np.sin(theta)], complex)
    kb2 = rho * omega ** 2 / (M if wave_type == WAVE_P else G)
    nu = _vertical_wavenumber(kb2, k)
    kb = np.sqrt(complex(kb2))
    if wave_type == WAVE_P:
        return np.array([k / kb, 0.0, nu / kb], complex)
    return np.array([nu / kb, 0.0, -k / kb], complex)


@quiet_fpe
def inclined_body_wave(col: Column, omega: float, wave_type: int, angle_deg: float, hs_rho: float,
                       hs_G: complex, hs_M: complex, vs_hs: float, vp_hs: float,
                       rigid_base: bool = False, mass: str = MASS_MIXED) -> WaveField:
    """Inclined P, SV or SH plane wave incident from the half-space (R1 §2.7b, D-SIT-08).

    ``col`` must be the column of the *user layers only* (its bottom interface is the top of the
    half-space; the base flag of ``col`` is ignored, ``rigid_base`` selects the base model).
    ``k = w sin(theta)/V_hs`` (V = Vs for SV/SH, Vp for P, undamped velocity).  With
    ``K_out``/``K_in`` from :func:`halfspace_stiffness` and the displacement ``u_I`` of the pure
    incident wave at the top of the half-space (:func:`incident_polarisation`), the column
    equation is

        (A k^2 + B k + G - w^2 M + K_out e_N e_N^T) U = (K_out - K_in) u_I  (at the base)

    A rigid base (``nl = 0``) instead moves with the incident polarisation ``u_I``.
    """
    theta = deg2rad(angle_deg)
    v = vp_hs if wave_type == WAVE_P else vs_hs
    k = omega * np.sin(theta) / v
    pol = incident_polarisation(wave_type, theta, k, hs_rho, hs_G, hs_M, omega)
    mats = column_matrices(col, mass)
    n = col.n_iface
    rigid = bool(rigid_base)
    U = np.zeros((n, 3), complex)
    out = np.zeros(3, complex)
    if wave_type == WAVE_SH:
        K = mats.Az * k * k + mats.Gs - omega ** 2 * mats.Mm
        uI = np.array([pol[1]], complex)
        Kout, Kin = halfspace_stiffness(hs_rho, hs_G, hs_M, omega, k, inplane=False)
        dofs = [n - 1]
        comps = [1]
    else:
        K = np.block([[mats.Ax * k * k + mats.Gs - omega ** 2 * mats.Mm, mats.Bxz * k],
                      [mats.Bxz.T * k, mats.Az * k * k + mats.Gp - omega ** 2 * mats.Mm]])
        uI = np.array([pol[0], 1j * pol[2]], complex)             # Kausel variables {u_x, i u_z}
        Kout, Kin = halfspace_stiffness(hs_rho, hs_G, hs_M, omega, k, inplane=True)
        dofs = [n - 1, 2 * n - 1]
        comps = [0, 2]
    dofs = np.asarray(dofs)
    if rigid:
        free = np.setdiff1d(np.arange(K.shape[0]), dofs)
        x = np.zeros(K.shape[0], complex)
        x[dofs] = uI
        x[free] = np.linalg.solve(K[np.ix_(free, free)], -K[np.ix_(free, dofs)] @ uI)
        out[comps] = pol[comps]
    else:
        K = K.copy()
        K[np.ix_(dofs, dofs)] += Kout
        f = np.zeros(K.shape[0], complex)
        f[dofs] = (Kout - Kin) @ uI
        x = np.linalg.solve(K, f)
        # outcrop: free surface of the half-space alone, K_out u = (K_out - K_in) u_I
        uo = np.linalg.solve(Kout, (Kout - Kin) @ uI)
        out[comps] = uo
    if wave_type == WAVE_SH:
        U[:, 1] = x
    else:
        U[:, 0] = x[:n]
        U[:, 2] = -1j * x[n:]                                      # back to physical u_z
        if not rigid:
            out[2] = -1j * out[2]
    return WaveField(U=U[col.iface_user], k=complex(k), outcrop=out)


# ---------------------------------------------------------------------------------------
# Normalisation to the control point (D-SIT-07)
# ---------------------------------------------------------------------------------------
def normalise(field: WaveField, cl: int, cm: int, rtol: float = 1e-10) -> Tuple[np.ndarray, complex]:
    """Divide the field by its motion at the top of control layer ``cl`` (1-based user interface)
    in direction ``cm`` (0 x', 1 y', 2 z').  Returns ``(U_normalised, u_cp)``; raises ValueError
    when the wave has no motion in the control direction at the control point."""
    ucp = complex(field.U[cl - 1, cm])
    scale = float(np.abs(field.U).max())
    if scale == 0 or abs(ucp) <= rtol * scale:
        raise ValueError(f"no motion in the control direction {DIRECTION_NAMES[cm]} at the control point")
    return field.U / ucp, ucp


# ---------------------------------------------------------------------------------------
# Binding API: free field at the interaction nodes (ANALYS)
# ---------------------------------------------------------------------------------------
@quiet_fpe
def free_field_at_nodes(file1: Container, q: int, xyz, depth_iface, ang_deg: float, xc: float, yc: float
                        ) -> np.ndarray:
    """Free-field motion (global x, y, z components) at nodes for unit control motion.

    ``file1``: FILE1 container; ``q``: 0-based frequency row of FILE1; ``xyz`` (n, 3) global node
    coordinates (only x, y are used: the depth is given by the interface); ``depth_iface`` (n,):
    1-based user interface of each node; ``ang_deg``: ANALYS coordinate transformation angle from
    x' to the global x axis; ``(xc, yc)``: control point (phase reference, x' = 0).

        U'(node) = R_z(ang) sum_w r_w U_w(iface) exp(-i k_w x'),  x' = (x-xc) cos a + (y-yc) sin a
    """
    U = np.asarray(file1["U"])
    nW, nF, nI, _ = U.shape
    if not 0 <= q < nF:
        raise IndexError(f"frequency row {q} outside FILE1 (nF = {nF})")
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    it = np.asarray(depth_iface, int).reshape(-1) - 1
    if it.size != xyz.shape[0]:
        raise ValueError("xyz and depth_iface have different lengths")
    if np.any(it < 0) or np.any(it >= nI):
        raise ValueError(f"interface numbers must lie in 1..{nI}")
    a = deg2rad(ang_deg)
    ca, sa = np.cos(a), np.sin(a)
    xp = (xyz[:, 0] - xc) * ca + (xyz[:, 1] - yc) * sa
    k = np.asarray(file1["k"])[:, q]
    r = np.asarray(file1["ratio"])[:, q]
    phase = np.exp(-1j * k[:, None] * xp[None, :])                  # (nW, n)
    Up = np.einsum("w,wn,wnc->nc", r, phase, U[:, q, it, :])          # x'y'z' components
    out = np.empty_like(Up)
    out[:, 0] = ca * Up[:, 0] - sa * Up[:, 1]
    out[:, 1] = sa * Up[:, 0] + ca * Up[:, 1]
    out[:, 2] = Up[:, 2]
    return out
