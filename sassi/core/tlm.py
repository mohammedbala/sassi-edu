"""Thin-layer method (TLM) for horizontally layered soil: the numerical core of SITE and POINT.

Physics in one paragraph
------------------------
A horizontally layered site is cut into thin sublayers.  Inside a sublayer the displacement
varies *linearly* with depth (finite elements in z) while the horizontal dependence is kept
analytic, ``exp(i w t - i k x)``.  For every sublayer this gives four small matrices A, B, G, M
(Kausel 1981, Table 1; requirements R1 §2.2) and the interface load-displacement relation

    P = (A k^2 + B k + G - w^2 M) U                                       (R1 §2.2)

Setting P = 0 gives eigenproblems whose eigenvalues are the horizontal wavenumbers k of the
Rayleigh (in-plane, P-SV) and Love (anti-plane, SH) modes of the *discrete* layered medium.
These modes are the building blocks of the free field (SITE Mode 2), of the transmitting
boundaries and of the point-load solutions (POINT).

Conventions (R1 §0, D-CNV-01): time factor ``exp(+i w t)``; horizontal propagation
``exp(-i k x)``, outgoing waves have Im k < 0 (Re k > 0 when k is real); **z up**; interfaces are
numbered from the ground surface downwards; in-plane vector ``{u_x, i u_z}`` ("i-scaled vertical",
Kausel form) so that all matrices are complex *symmetric*.  Physical vertical displacement is
``u_z = -i * phi_z``.

Half-space (R1 §2.4, D-SIT-02): ``nl`` extra sublayers with the half-space properties and total
thickness ``1.5 * Vs_hs / f`` are added below the user layers ("variable-depth method"), and
Lysmer-Kuhlemeyer dashpots ``c_s = rho V*_s`` (x, y) and ``c_p = rho V*_p`` (z) are attached at the
base ("viscous boundary").  ``nl = 0`` means a rigid base (the base interface is fixed).

Mass (D-CNV-05, R1 §2.1): every sublayer uses the mixed mass ``M = 1/2 M_consistent + 1/2
M_lumped`` (4th-order accurate shear-wave dispersion, R2 A.3).
"""
from __future__ import annotations

import functools
from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np
from scipy import linalg as sla

from ..conventions import CMOD_SASSI, cfactor

MASS_MIXED = "mixed"
MASS_CONSISTENT = "consistent"
MASS_LUMPED = "lumped"

BASE_RIGID = "rigid"
BASE_DASHPOT = "dashpot"

HS_GEOMETRIC = "geometric"
HS_UNIFORM = "uniform"
HS_LINEAR = "linear"
HS_LAWS = (HS_GEOMETRIC, HS_UNIFORM, HS_LINEAR)


def quiet_fpe(func):
    """Silence floating-point *flags* raised inside ``func``.

    numpy 2.0 linked to Apple Accelerate sets spurious divide/overflow/invalid flags in
    ``matmul`` on perfectly finite data (sizes >= ~50), which numpy reports as RuntimeWarnings.
    The SASSI-EDU numerics check finiteness explicitly where it matters, so these flags are
    ignored inside the layered-soil routines.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with np.errstate(divide="ignore", over="ignore", invalid="ignore", under="ignore"):
            return func(*args, **kwargs)
    return wrapper


# 2x2 shape-function integrals of a linear element of unit thickness (top node first)
_NN = np.array([[2.0, 1.0], [1.0, 2.0]]) / 6.0          # int N^T N dz / h
_DD = np.array([[1.0, -1.0], [-1.0, 1.0]])              # h * int N'^T N' dz


def mass_pattern(mass: str = MASS_MIXED) -> np.ndarray:
    """2x2 sublayer mass pattern (to be multiplied by rho*h).  R1 §2.1/§2.2, D-CNV-05."""
    if mass == MASS_CONSISTENT:
        return _NN.copy()
    if mass == MASS_LUMPED:
        return 0.5 * np.eye(2)
    if mass == MASS_MIXED:
        return 0.5 * _NN + 0.25 * np.eye(2)                # [[5/12, 1/12], [1/12, 5/12]]
    raise ValueError(f"unknown mass option {mass!r}")


# ---------------------------------------------------------------------------------------
# Material helpers
# ---------------------------------------------------------------------------------------
def layer_moduli(rho, vs, vp, ds, dp, form: int = CMOD_SASSI):
    """Complex (G*, M*) of layers from density, velocities and damping ratios (D-CNV-03/04).

    ``G* = rho Vs^2 c(ds)`` and ``M* = lam* + 2G* = rho Vp^2 c(dp)``; ``c`` is
    :func:`sassi.conventions.cfactor` (SASSI form by default).
    """
    rho = np.asarray(rho, float)
    G = rho * np.asarray(vs, float) ** 2 * np.asarray(cfactor(np.asarray(ds, float), form))
    M = rho * np.asarray(vp, float) ** 2 * np.asarray(cfactor(np.asarray(dp, float), form))
    return np.asarray(G, complex), np.asarray(M, complex)


def complex_speeds(rho, G, M):
    """Complex wave speeds ``V*_s = sqrt(G*/rho)``, ``V*_p = sqrt(M*/rho)`` (principal roots)."""
    rho = np.asarray(rho, float)
    return np.sqrt(np.asarray(G, complex) / rho), np.sqrt(np.asarray(M, complex) / rho)


# ---------------------------------------------------------------------------------------
# Variable-depth half-space (D-SIT-02, UT-19)
# ---------------------------------------------------------------------------------------
def halfspace_depth(vs_hs: float, f: float) -> float:
    """Total thickness of the generated half-space sublayers: ``1.5 * lambda_s = 1.5 Vs_hs/f``
    (R1 §2.4, manual "Rigid Base Rock vs. Halfspace Condition")."""
    if f <= 0:
        raise ValueError("frequency must be positive for the variable-depth half-space")
    return 1.5 * vs_hs / f


def halfspace_sublayers(vs_hs: float, f: float, nl: int, h_last: float = 0.0, vs_last: float = 0.0,
                        law: str = HS_GEOMETRIC) -> Tuple[np.ndarray, str]:
    """Thicknesses of the ``nl`` generated half-space sublayers at frequency ``f`` (D-SIT-02).

    * total thickness ``H = 1.5 lambda_s`` with ``lambda_s = Vs_hs / f``;
    * ``law='geometric'`` (default): ``h_i = h1 q^(i-1)``, ``h1 = min(h_last Vs_hs/Vs_last, H/nl)``,
      ``q >= 1`` found by bisection so that ``sum h_i = H``; every ``h_i <= lambda_s/8``, otherwise
      (or when ``h_last``/``vs_last`` are not given) the uniform law is used ("uniform fallback");
    * ``law='linear'``: ``h_i = h1 + (i-1) d`` with the same ``h1``, cap and fallback;
    * ``law='uniform'``: ``h_i = H/nl``.

    The thicknesses are non-decreasing with depth.  Returns ``(h, law_used)``.
    """
    if nl <= 0:
        return np.zeros(0), HS_UNIFORM
    law = (law or HS_GEOMETRIC).lower()
    if law not in HS_LAWS:
        raise ValueError(f"unknown half-space sublayer law {law!r} (geometric, uniform, linear)")
    H = halfspace_depth(vs_hs, f)
    lam = vs_hs / f
    uniform = np.full(nl, H / nl)
    if law == HS_UNIFORM or nl == 1 or h_last <= 0 or vs_last <= 0:
        return uniform, HS_UNIFORM
    h1 = min(h_last * vs_hs / vs_last, H / nl)
    if h1 >= H / nl * (1.0 - 1e-12):
        return uniform, HS_UNIFORM
    i = np.arange(nl)
    if law == HS_LINEAR:
        d = (H - nl * h1) / (nl * (nl - 1) / 2.0)
        h = h1 + i * d
    else:
        # sum_{i<n} h1 q^i = H  -> bisection on q in (1, qmax); the sum is increasing in q
        def total(q):
            return h1 * np.sum(q ** i)
        lo, hi = 1.0, 2.0
        while total(hi) < H:
            hi *= 2.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if total(mid) < H:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-15 * hi:
                break
        q = 0.5 * (lo + hi)
        h = h1 * q ** i
    h *= H / h.sum()                       # remove the residual of the bisection exactly
    if h.max() > lam / 8.0 * (1.0 + 1e-12):
        return uniform, HS_UNIFORM          # cap violated -> uniform fallback
    return h, law


# ---------------------------------------------------------------------------------------
# Soil column
# ---------------------------------------------------------------------------------------
@dataclass
class Column:
    """Discretised soil column, sublayers listed from the ground surface downwards.

    ``h, rho, G, M`` are per sublayer (complex moduli ``G* = rho Vs^2 c(bs)`` and
    ``M* = lam* + 2G*``).  ``base`` is ``'rigid'`` (bottom interface fixed) or ``'dashpot'``
    (Lysmer-Kuhlemeyer dashpots ``cs`` on x, y and ``cp`` on z at the bottom interface).
    ``iface_user[i]`` is the column interface index (0-based) of user interface ``i+1``.
    ``n_gen`` is the number of generated half-space sublayers at the bottom of the column.
    """

    h: np.ndarray
    rho: np.ndarray
    G: np.ndarray
    M: np.ndarray
    base: str = BASE_RIGID
    cs: complex = 0j
    cp: complex = 0j
    iface_user: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    n_gen: int = 0

    def __post_init__(self):
        self.h = np.asarray(self.h, float)
        self.rho = np.asarray(self.rho, float)
        self.G = np.asarray(self.G, complex)
        self.M = np.asarray(self.M, complex)
        if self.iface_user is None or len(self.iface_user) == 0:
            self.iface_user = np.arange(len(self.h) + 1)
        self.iface_user = np.asarray(self.iface_user, int)
        if np.any(self.h <= 0):
            raise ValueError("sublayer thicknesses must be positive")

    @property
    def lam(self) -> np.ndarray:
        return self.M - 2.0 * self.G

    @property
    def n_sub(self) -> int:
        return len(self.h)

    @property
    def n_iface(self) -> int:
        return len(self.h) + 1

    @property
    def n_free(self) -> int:
        """Number of free interfaces (the rigid base is removed)."""
        return self.n_iface if self.base == BASE_DASHPOT else self.n_iface - 1

    @property
    def depth(self) -> np.ndarray:
        """Depth of every column interface (positive down from the surface)."""
        return np.concatenate([[0.0], np.cumsum(self.h)])

    def user_free_index(self) -> np.ndarray:
        """Free-interface index of each user interface (-1 for a removed rigid base)."""
        idx = self.iface_user.copy()
        idx[idx >= self.n_free] = -1
        return idx


def build_column(thick, rho, G, M, hs_rho: float, hs_G: complex, hs_M: complex,
                 h_gen=None) -> Column:
    """Column of the user layers (as entered, D-SIT-01) plus generated half-space sublayers.

    ``thick, rho, G, M`` describe the TOPL layers (top first).  ``h_gen`` are the generated
    sublayer thicknesses (empty or None -> rigid base at the bottom of the last user layer).
    With generated sublayers the base carries dashpots ``c = rho_hs V*`` (R1 §2.4 item 1).
    """
    thick = np.asarray(thick, float)
    n_user = len(thick)
    h_gen = np.zeros(0) if h_gen is None else np.asarray(h_gen, float)
    n_gen = len(h_gen)
    h = np.concatenate([thick, h_gen])
    r = np.concatenate([np.asarray(rho, float), np.full(n_gen, hs_rho)])
    g = np.concatenate([np.asarray(G, complex), np.full(n_gen, hs_G, dtype=complex)])
    m = np.concatenate([np.asarray(M, complex), np.full(n_gen, hs_M, dtype=complex)])
    if n_gen > 0:
        cs = np.sqrt(hs_rho * complex(hs_G))     # rho V*_s = sqrt(rho G*)
        cp = np.sqrt(hs_rho * complex(hs_M))     # rho V*_p = sqrt(rho M*)
        base = BASE_DASHPOT
    else:
        cs = cp = 0j
        base = BASE_RIGID
    return Column(h=h, rho=r, G=g, M=m, base=base, cs=cs, cp=cp,
                  iface_user=np.arange(n_user + 1), n_gen=n_gen)


def sub_column(col: Column, first: int) -> Column:
    """The part of ``col`` below column interface ``first`` (same base), e.g. the half-space
    sublayers only, used for the discrete outcrop motion (R1 §2.7a)."""
    return Column(h=col.h[first:], rho=col.rho[first:], G=col.G[first:], M=col.M[first:], base=col.base,
                  cs=col.cs, cp=col.cp, iface_user=np.array([0]), n_gen=min(col.n_gen, col.n_sub - first))


# ---------------------------------------------------------------------------------------
# Assembled layer matrices
# ---------------------------------------------------------------------------------------
def _assemble(blocks: np.ndarray, n: int) -> np.ndarray:
    """Overlap per-sublayer 2x2 blocks at shared interfaces (FE assembly, R1 §2.2)."""
    out = np.zeros((n, n), dtype=complex)
    j = np.arange(blocks.shape[0])
    out[j, j] += blocks[:, 0, 0]
    out[j, j + 1] += blocks[:, 0, 1]
    out[j + 1, j] += blocks[:, 1, 0]
    out[j + 1, j + 1] += blocks[:, 1, 1]
    return out


@dataclass
class ColumnMatrices:
    """Assembled TLM matrices of a column on *all* interfaces (rigid base included).

    Direction sub-blocks (R1 §2.2): ``Ax`` uses lam+2G, ``Az = Ay`` use G; ``Gs = Gx = Gy`` use
    G, ``Gp = Gz`` uses lam+2G; ``Bxz`` is the x-row/z-column coupling block; ``Mm`` the mass;
    ``ds``/``dp`` the dashpot coefficients on the bottom interface (0 for a rigid base).
    Use :meth:`free` to obtain the blocks on the free interfaces.
    """

    Ax: np.ndarray
    Az: np.ndarray
    Gs: np.ndarray
    Gp: np.ndarray
    Bxz: np.ndarray
    Mm: np.ndarray
    ds: complex
    dp: complex
    n_free: int

    def free(self) -> "ColumnMatrices":
        s = slice(0, self.n_free)
        return ColumnMatrices(self.Ax[s, s], self.Az[s, s], self.Gs[s, s], self.Gp[s, s], self.Bxz[s, s],
                              self.Mm[s, s], self.ds, self.dp, self.n_free)

    @property
    def Ay(self) -> np.ndarray:
        return self.Az

    def dashpot(self, which: str, omega: float) -> np.ndarray:
        """``i w c`` on the bottom interface (diagonal matrix), ``which`` = 's' or 'p'."""
        n = self.Ax.shape[0]
        D = np.zeros((n, n), dtype=complex)
        c = self.ds if which == "s" else self.dp
        if c != 0 and n == self.n_free:
            D[-1, -1] = 1j * omega * c
        return D

    def C(self, omega: float):
        """``(Cx, Cz, Cy) = (G_dir - w^2 M + i w D_dir)`` (dashpots included, R1 §2.5)."""
        w2M = omega ** 2 * self.Mm
        Cx = self.Gs - w2M + self.dashpot("s", omega)
        Cz = self.Gp - w2M + self.dashpot("p", omega)
        Cy = Cx
        return Cx, Cz, Cy

    def K_inplane(self, omega: float, k: complex) -> np.ndarray:
        """In-plane dynamic stiffness ``[[Ax k^2 + Cx, Bxz k], [Bxz^T k, Az k^2 + Cz]]`` acting on
        ``{u_x; i u_z}`` (Kausel form, R1 §2.2)."""
        Cx, Cz, _ = self.C(omega)
        return np.block([[self.Ax * k * k + Cx, self.Bxz * k], [self.Bxz.T * k, self.Az * k * k + Cz]])

    def K_antiplane(self, omega: float, k: complex) -> np.ndarray:
        """Anti-plane (SH) dynamic stiffness ``Ay k^2 + Cy``."""
        _, _, Cy = self.C(omega)
        return self.Ay * k * k + Cy


def column_matrices(col: Column, mass: str = MASS_MIXED) -> ColumnMatrices:
    """Assemble the Kausel-form TLM matrices (R1 §2.2, Table 1 of Kausel 1981) of ``col``."""
    h = col.h[:, None, None]
    G = col.G[:, None, None]
    lam = col.lam[:, None, None]
    lp2 = col.M[:, None, None]
    n = col.n_iface
    Ax = _assemble(h * lp2 * _NN, n)
    Az = _assemble(h * G * _NN, n)
    Gs = _assemble(G / h * _DD, n)
    Gp = _assemble(lp2 / h * _DD, n)
    bxz = 0.5 * np.stack([np.stack([lam - G, -(lam + G)], -1), np.stack([lam + G, -(lam - G)], -1)], -2)
    Bxz = _assemble(bxz.reshape(-1, 2, 2), n)
    Mm = _assemble(col.rho[:, None, None] * h * mass_pattern(mass), n)
    return ColumnMatrices(Ax, Az, Gs, Gp, Bxz, Mm, complex(col.cs), complex(col.cp), col.n_free)


def boundary_integrals(col: Column, coef: str) -> Tuple[np.ndarray, np.ndarray]:
    """Consistent boundary integrals on the free interfaces (R1 §2.6):

    ``E_c = int c N N^T dz`` (per sublayer ``c h/6 [2 1; 1 2]``) and
    ``Q_c = int c N N'^T dz`` (z up, top node first, per sublayer ``c/2 [1 -1; 1 -1]``),
    with ``c`` = ``'lp2'`` (lam+2G), ``'lam'`` or ``'G'``.  Used by the transmitting boundaries
    of POINT3 (R1 §3.3 iii) and POINT2 (R1 §2.6).
    """
    c = {"lp2": col.M, "lam": col.lam, "G": col.G}[coef][:, None, None]
    E = _assemble(c * col.h[:, None, None] * _NN, col.n_iface)
    Q = _assemble(c * 0.5 * np.array([[1.0, -1.0], [1.0, -1.0]]), col.n_iface)
    s = slice(0, col.n_free)
    return E[s, s], Q[s, s]


# ---------------------------------------------------------------------------------------
# Eigenproblems (Mode 1)
# ---------------------------------------------------------------------------------------
def outgoing_root(k2: np.ndarray) -> np.ndarray:
    """``k = sqrt(k^2)`` on the branch Im k < 0, or Re k > 0 when k is real (R1 §2.5, [S4] p.15)."""
    k = np.sqrt(np.asarray(k2, dtype=complex))
    scale = np.maximum(np.abs(k), 1e-300)
    real_like = np.abs(k.imag) <= 1e-12 * scale
    flip = np.where(real_like, k.real < 0, k.imag > 0)
    return np.where(flip, -k, k)


@dataclass
class Modes:
    """Rayleigh and Love eigen-solutions of a column at one frequency (R1 §2.5).

    ``kR`` (2Nf,), ``phix`` and ``phiz`` (Nf, 2Nf): Rayleigh wavenumbers and Kausel-normalised mode
    shapes (``phiz`` is the i-scaled vertical component, physical ``u_z = -i phiz``);
    ``kL`` (Nf,), ``phiy`` (Nf, Nf): Love modes normalised by ``phiy^T Ay phiy = 1``.
    """

    kR: np.ndarray
    phix: np.ndarray
    phiz: np.ndarray
    kL: np.ndarray
    phiy: np.ndarray

    @property
    def n_free(self) -> int:
        return self.phix.shape[0]


def _sort_modes(k: np.ndarray) -> np.ndarray:
    """Deterministic order: decreasing Re k (slowest = shortest-wavelength modes first)."""
    return np.lexsort((np.round(k.imag, 12), -np.round(k.real, 12)))


@quiet_fpe
def rayleigh_modes(mats: ColumnMatrices, omega: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generalised Rayleigh (P-SV) modes, linearised without doubling the size (R1 §2.5, [S4] Eq. 16):

        (k^2 [Ax 0; Bxz^T Az] + [Cx Bxz; 0 Cz]) {phi_x; k phi_z} = 0

    solved as the standard eigenproblem of ``-Abar^-1 Cbar`` (Abar is block lower triangular).
    Each mode is normalised by ``Y^T Abar Z = k`` with ``Y = {k phi_x; phi_z}`` ([S4] Eq. 22a),
    i.e. ``phi_x^T Ax phi_x + phi_z^T Az phi_z + phi_z^T Bxz^T phi_x / k = 1``.
    """
    m = mats.free()
    n = m.n_free
    Cx, Cz, _ = m.C(omega)
    Ax, Az, Bxz = m.Ax, m.Az, m.Bxz
    # Y = Abar^-1 Cbar by block forward substitution
    top = sla.solve(Ax, np.hstack([Cx, Bxz]), assume_a="sym", check_finite=False)
    rhs = np.hstack([np.zeros((n, n), complex), Cz]) - Bxz.T @ top
    bot = sla.solve(Az, rhs, assume_a="sym", check_finite=False)
    X = -np.vstack([top, bot])
    k2, Z = sla.eig(X, check_finite=False)
    k = outgoing_root(k2)
    order = _sort_modes(k)
    k, Z = k[order], Z[:, order]
    phx = Z[:n, :].astype(complex)
    phz = Z[n:, :] / k[None, :]
    # normalisation (Kausel Eq. 22a)
    nrm = (np.einsum("ij,ik,kj->j", phx, Ax, phx) + np.einsum("ij,ik,kj->j", phz, Az, phz)
           + np.einsum("ij,ki,kj->j", phz, Bxz, phx) / k)
    s = 1.0 / np.sqrt(nrm)
    return k, phx * s[None, :], phz * s[None, :]


@quiet_fpe
def love_modes(mats: ColumnMatrices, omega: float) -> Tuple[np.ndarray, np.ndarray]:
    """Love (SH) modes ``(Ay k^2 + Cy) phi_y = 0`` normalised by ``phi_y^T Ay phi_y = 1`` (R1 §2.5)."""
    m = mats.free()
    _, _, Cy = m.C(omega)
    X = -sla.solve(m.Ay, Cy, assume_a="sym", check_finite=False)
    k2, V = sla.eig(X, check_finite=False)
    k = outgoing_root(k2)
    order = _sort_modes(k)
    k, V = k[order], V[:, order].astype(complex)
    nrm = np.einsum("ij,ik,kj->j", V, m.Ay, V)
    return k, V / np.sqrt(nrm)[None, :]


@quiet_fpe
def column_modes(col: Column, omega: float, mass: str = MASS_MIXED) -> Modes:
    """Rayleigh + Love modes of ``col`` at circular frequency ``omega`` (SITE Mode 1)."""
    mats = column_matrices(col, mass)
    kR, phx, phz = rayleigh_modes(mats, omega)
    kL, phy = love_modes(mats, omega)
    return Modes(kR, phx, phz, kL, phy)


# ---------------------------------------------------------------------------------------
# Wavenumber-domain flexibility (verification of the modal identities, R1 §2.5, [S4] Eq. 49-52)
# ---------------------------------------------------------------------------------------
@quiet_fpe
def modal_flexibility(modes: Modes, k: complex) -> Tuple[np.ndarray, np.ndarray]:
    """In-plane ``[[Fxx, Fxz], [Fzx, Fzz]]`` and anti-plane ``Fyy`` from the modes:

    ``Fxx = Phi_x D Phi_x^T``, ``Fxz = k Phi_x K^-1 D Phi_z^T``, ``Fzz = Phi_z D Phi_z^T``,
    ``Fyy = Phi_y D_L Phi_y^T`` with ``D = diag(1/(k^2 - k_j^2))``.
    """
    kR, phx, phz = modes.kR, modes.phix, modes.phiz
    D = 1.0 / (k * k - kR * kR)
    Fxx = (phx * D) @ phx.T
    Fxz = k * (phx * (D / kR)) @ phz.T
    Fzx = k * (phz * (D / kR)) @ phx.T
    Fzz = (phz * D) @ phz.T
    DL = 1.0 / (k * k - modes.kL ** 2)
    Fyy = (modes.phiy * DL) @ modes.phiy.T
    return np.block([[Fxx, Fxz], [Fzx, Fzz]]), Fyy


@quiet_fpe
def direct_flexibility(col: Column, omega: float, k: complex, mass: str = MASS_MIXED):
    """Direct inverse of the in-plane and anti-plane dynamic stiffness at wavenumber ``k``."""
    m = column_matrices(col, mass).free()
    return np.linalg.inv(m.K_inplane(omega, k)), np.linalg.inv(m.K_antiplane(omega, k))


# ---------------------------------------------------------------------------------------
# Mode selection for surface-wave fields (D-SIT-09)
# ---------------------------------------------------------------------------------------
# Which modes are "propagating"?
#
# In *elastic* soil a propagating mode has a real wavenumber; evanescent (near-field) modes are
# imaginary or complex.  SASSI-EDU calls a mode propagating when it decays by less than exp(-pi)
# per wavelength, |Im k| <= 0.5 Re k, i.e. |arg k| <= atan(0.5).
#
# Material damping rotates the wavenumbers: a plane body wave in a sublayer with complex modulus
# G* = |G*| exp(i delta) has k = w sqrt(rho/G*), so arg k = -delta/2 and
# |Im k|/Re k = tan(delta/2) (= beta/sqrt(1 - beta^2) for the SASSI form, 0.577 at beta = 0.5).
# For a uniform column *every* eigenvalue is rotated by exactly this angle, so the elastic
# classification is preserved if the propagating sector is widened by the material loss angle:
#
#     |arg k| <= atan(0.5) + atan(loss),   loss = max tan(delta/2) over the column.
#
# A fixed sector would reject every mode of a heavily damped column (beta > 0.447, which EDU-04
# admits) although the modes are the physical surface waves.
PROPAGATING_RATIO = 0.5   # elastic soil: |Im k| <= 0.5 Re k, decays less than exp(-pi) per wavelength
TIE_RTOL = 1e-8           # |Im k| values within TIE_RTOL * max|k| are equal ("least decay" ties)
MAX_MATERIAL_LOSS = float(np.tan(np.pi / 6))   # tan(delta/2) at beta -> 0.5 (SASSI form), EDU-04


def material_loss(G, M=None) -> float:
    """Largest material loss ratio ``tan(delta/2)`` of the complex moduli G* (and M*).

    ``delta = arg G*``; ``tan(delta/2) = |Im k| / Re k`` of a plane body wave in that material
    (D-CNV-03/04: ``beta / sqrt(1 - beta^2)`` for the SASSI form, ``tan(atan(2 beta)/2)`` for
    ``1 + 2 i beta``).  This is the decay that material damping alone causes (D-SIT-09 note).
    """
    parts = [np.ravel(np.asarray(G, complex))]
    if M is not None:
        parts.append(np.ravel(np.asarray(M, complex)))
    z = np.concatenate(parts)
    if z.size == 0:
        return 0.0
    return float(np.max(np.abs(np.tan(0.5 * np.angle(z)))))


def column_loss(col: Column) -> float:
    """:func:`material_loss` over all sublayers of ``col`` (user + generated half-space)."""
    return material_loss(col.G, col.M)


def spectrum_loss(k: np.ndarray) -> float:
    """Estimate of the material loss when the column is not at hand: the decay ratio
    ``|Im k| / Re k`` of the least-attenuated mode with Re k > 0.

    Exact for a uniform column (all modes are rotated by the same loss angle), a lower bound for
    a layered one -- callers that know the column pass :func:`column_loss` instead.  The estimate
    is capped at the largest material loss EDU-04 admits (beta < 0.5: ``tan 30 deg``), so a
    spectrum without any propagating mode is not mistaken for heavy damping.
    """
    k = np.asarray(k)
    pos = k.real > 0
    if not np.any(pos):
        return 0.0
    return float(min(np.min(np.abs(k.imag[pos]) / k.real[pos]), MAX_MATERIAL_LOSS))


def propagating_ratio(loss: float = 0.0) -> float:
    """Threshold of ``|Im k| / Re k`` for a propagating mode in a column of material loss ``loss``:
    the elastic sector ``|arg k| <= atan(0.5)`` widened by the loss angle ``atan(loss)``
    (``0.5`` for elastic soil, ``0.564`` at beta = 5 %, ``1.47`` at beta = 0.49)."""
    # tan(a + b) = (tan a + tan b) / (1 - tan a tan b), exact 0.5 for elastic soil
    t, L = PROPAGATING_RATIO, max(float(loss), 0.0)
    return (t + L) / (1.0 - t * L) if t * L < 1.0 else np.inf


def propagating(k: np.ndarray, ratio: float = PROPAGATING_RATIO) -> np.ndarray:
    """Boolean mask of propagating modes: Re k > 0 and |Im k| <= ratio * Re k
    (``ratio`` from :func:`propagating_ratio` for damped soil)."""
    k = np.asarray(k)
    return (k.real > 0) & (np.abs(k.imag) <= ratio * k.real)


@dataclass(frozen=True)
class ModeChoice:
    """Result of the surface-wave mode selection (D-SIT-09), for the SITE listing.

    ``index`` selected mode; ``rule`` 1 shortest wavelength / 2 least decay; ``ratio`` the
    propagating threshold |Im k|/Re k used by rule 1; ``n_candidates`` modes eligible under the
    rule; ``n_tied`` modes whose |Im k| equals the least decay within the tolerance (rule 2; 1
    when the choice is unique).
    """

    index: int
    rule: int
    ratio: float
    n_candidates: int
    n_tied: int = 1


def choose_mode(k: np.ndarray, rule: int, loss: Optional[float] = None) -> ModeChoice:
    """Surface-wave mode used by SITE Mode 2 (D-SIT-09), with selection diagnostics.

    ``rule = 1`` "Shortest Wavelength": the propagating mode (|Im k| <= ratio Re k, ratio from
    :func:`propagating_ratio` of the column loss) with the largest Re k.

    ``rule = 2`` "Least Decay": the mode with Re k > 0 and the smallest |Im k|.  Values of |Im k|
    within ``TIE_RTOL * max|k|`` of the smallest one are *equal* -- e.g. every propagating mode of
    an undamped column on a rigid base has Im k = 0 in exact arithmetic, and only rounding noise
    would separate them.  Ties are broken deterministically by the largest Re k (the fundamental,
    shortest-wavelength mode), so the field does not depend on floating-point noise (D-GEN-09).

    ``loss`` is the column's material loss (:func:`column_loss`); ``None`` estimates it from the
    spectrum (:func:`spectrum_loss`).
    """
    k = np.asarray(k)
    if loss is None:
        loss = spectrum_loss(k)
    ratio = propagating_ratio(loss)
    if rule == 2:
        cand = np.flatnonzero(k.real > 0)
        if cand.size == 0:
            raise ValueError("no mode with Re k > 0")
        decay = np.abs(k.imag[cand])
        tol = TIE_RTOL * float(np.max(np.abs(k)))
        tied = cand[decay <= decay.min() + tol]
        j = tied[np.argmax(k.real[tied])]
        return ModeChoice(int(j), 2, ratio, int(cand.size), int(tied.size))
    cand = np.flatnonzero(propagating(k, ratio))
    if cand.size == 0:
        raise ValueError(f"no propagating mode (|Im k| <= {ratio:.3g} Re k)")
    return ModeChoice(int(cand[np.argmax(k.real[cand])]), 1, ratio, int(cand.size), 1)


def select_mode(k: np.ndarray, rule: int, loss: Optional[float] = None) -> int:
    """Index of the surface-wave mode used by SITE Mode 2 (D-SIT-09); see :func:`choose_mode`."""
    return choose_mode(k, rule, loss).index


# ---------------------------------------------------------------------------------------
# 1-D column response (vertically propagating waves, R1 §2.7a)
# ---------------------------------------------------------------------------------------
@quiet_fpe
def vertical_response(col: Column, omega: float, kind: str, mass: str = MASS_MIXED,
                      mats: Optional[ColumnMatrices] = None) -> np.ndarray:
    """Displacement at every column interface for a vertically incident wave (k = 0).

    ``kind = 's'`` (SV or SH: shear column with G*) or ``'p'`` (P: column with M* = lam*+2G*).
    Dashpot base: incident (up-going) wave of unit displacement amplitude at the dashpot level,

        (G - w^2 M + i w c e_N e_N^T) u = 2 i w c E_b e_N,   E_b = 1        (R1 §2.7a)

    Rigid base: the base interface moves with unit amplitude.  Returns ``u`` (n_iface,) complex.
    """
    m = mats if mats is not None else column_matrices(col, mass)
    Gd = m.Gs if kind == "s" else m.Gp
    K = Gd - omega ** 2 * m.Mm
    n = col.n_iface
    if col.base == BASE_DASHPOT:
        c = col.cs if kind == "s" else col.cp
        K = K.copy()
        K[-1, -1] += 1j * omega * c
        f = np.zeros(n, complex)
        f[-1] = 2j * omega * c
        return np.linalg.solve(K, f)
    u = np.zeros(n, complex)
    u[-1] = 1.0
    nf = n - 1
    u[:nf] = np.linalg.solve(K[:nf, :nf], -K[:nf, nf])
    return u


def outcrop_response(col: Column, omega: float, kind: str, mass: str = MASS_MIXED) -> complex:
    """Discrete outcrop motion at the top of the half-space (R1 §2.7a, [V] test4b):
    the surface motion of the half-space-only column (same generated sublayers, dashpot and
    incident wave).  For a rigid base the outcrop motion equals the base motion (1)."""
    if col.base != BASE_DASHPOT or col.n_gen == 0:
        return 1.0 + 0j
    sub = sub_column(col, col.n_sub - col.n_gen)
    return complex(vertical_response(sub, omega, kind, mass)[0])


# ---------------------------------------------------------------------------------------
# FILE2 helpers (SITE Mode 1 output, read by SITE Mode 2 and POINT)
# ---------------------------------------------------------------------------------------
def column_from_file2(file2, q: int) -> Column:
    """Rebuild the discretised column of frequency row ``q`` (0-based) of a FILE2 container:
    user layers (``layer_*``), generated sublayers ``h_gen[q]`` with the half-space properties
    (last entry of ``layer_*``) and the base dashpots (requirements §4.3: POINT uses the SITE
    sublayers, generated sublayers and base dashpots)."""
    a = file2.arrays
    thick = np.asarray(a["layer_thick"], float)
    rho = np.asarray(a["layer_rho"], float)
    G = np.asarray(a["layer_G"], complex)
    M = np.asarray(a["layer_M"], complex)
    h_gen = np.asarray(a["h_gen"], float)
    hg = h_gen[q] if h_gen.size else np.zeros(0)
    return build_column(thick, rho[:-1], G[:-1], M[:-1], rho[-1], G[-1], M[-1], hg)


def modes_from_file2(file2, q: int) -> Modes:
    """Rayleigh/Love modes of frequency row ``q`` of a FILE2 container."""
    a = file2.arrays
    return Modes(kR=a["kR"][q], phix=a["phix"][q], phiz=a["phiz"][q], kL=a["kL"][q], phiy=a["phiy"][q])
