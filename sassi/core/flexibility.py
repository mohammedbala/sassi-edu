"""Free-field flexibility matrix F_ff at the interaction nodes, assembled from FILE3 (ANALYS step 2).

Binding API (docs/ARCHITECTURE.md §6.3)::

    flexibility_matrix(file3, q, xy, iface) -> (3n, 3n) complex      # node-major [ux1,uy1,uz1,ux2,...]

How it works (R1 §4.1, requirements §4.6 item 2, D-ANL-02)
--------------------------------------------------------
POINT stored, for every load interface n and frequency, the response of the layered soil to a
unit point load on the axis of its central zone, split into Fourier harmonics: mu = 1
(horizontal load) and mu = 0 (vertical load).  For an interaction node i on interface m and a
node j on interface n, at horizontal distance r and angle ``theta = atan2(y_i - y_j, x_i - x_j)``
(from the load point j to the observation point i), the POINT values

* ``u~, v~, w~`` = (u_rho, u_theta, u_z) amplitudes for the horizontal load (mu = 1),
* ``p~, q~``     = (u_rho, u_z) for the vertical load (mu = 0)

are rotated into the Cartesian 3x3 block (rows u_x, u_y, u_z at i; columns P_x, P_y, P_z at j)::

    F_ij = [ u c^2 + v s^2    (u - v) s c     p c ]
           [ (u - v) s c      u s^2 + v c^2   p s ]        c = cos(theta), s = sin(theta)
           [ w c              w s             q   ]

* ``r = 0`` (same vertical line): ``diag(u_axis, u_axis, q_axis)`` from the core-axis values;
* ``0 < r < R0`` (mesh finer than the central zone): linear in r between the axis value and the
  value at R0, consistent with the core shape functions;
* ``r >= R0``: exact evaluation of the outgoing-mode expansion ``Psi_mu(r) alpha`` (D-PNT-02).

Reciprocity ``F_ij = F_ji^T`` holds only approximately with the FE core (~1e-4, R1 V10); the
matrix is symmetrised ``F <- (F + F^T)/2`` (D-ANL-02) unless ``symmetrize=False``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
from scipy.interpolate import CubicSpline

from ..io.container import Container
from .axisym import exterior_field
from .tlm import propagating, propagating_ratio, quiet_fpe, spectrum_loss

COMPONENTS = ("u", "v", "w", "p", "q")
MAX_EXACT = 6000          # unique distances evaluated exactly; above this a dense table is used


@dataclass
class PointData:
    """One frequency row of FILE3 (POINT3)."""

    R0: float
    load_iface: np.ndarray        # (nL,) user interfaces 1..L+1
    kR: np.ndarray
    kL: np.ndarray
    phix: np.ndarray              # (nL, 2Nf) mode rows at the observation (= load) interfaces
    phiz: np.ndarray
    phiy: np.ndarray              # (nL, Nf)
    alpha0: np.ndarray            # (nL, 3Nf)
    alpha1: np.ndarray
    axis0: np.ndarray             # (nL load, nL obs, 3)
    axis1: np.ndarray
    dim: int = 2
    loss: Optional[float] = None  # material loss tan(delta/2) of the column (FILE3 x_loss), None: unknown


def point_data(file3: Container, q: int) -> PointData:
    """Extract frequency row ``q`` (0-based index into ``file3['fnum']``) of FILE3."""
    nF = len(file3["fnum"])
    if not 0 <= q < nF:
        raise IndexError(f"frequency row {q} outside FILE3 (nF = {nF})")
    a = file3.arrays
    return PointData(R0=float(file3.meta["R0"]), load_iface=np.asarray(a["load_iface"], int),
                     kR=a["kR"][q], kL=a["kL"][q], phix=a["phix_obs"][q], phiz=a["phiz_obs"][q],
                     phiy=a["phiy_obs"][q], alpha0=a["alpha0"][q], alpha1=a["alpha1"][q],
                     axis0=a["axis0"][q], axis1=a["axis1"][q], dim=int(file3.meta.get("dim", 2)),
                     loss=float(a["x_loss"][q]) if "x_loss" in a else None)


def frequency_row(container: Container, fnum: int) -> int:
    """0-based row of frequency number ``fnum`` in a FILE1/FILE2/FILE3 container (D-ANL-08)."""
    idx = np.flatnonzero(np.asarray(container["fnum"]) == int(fnum))
    if idx.size == 0:
        raise KeyError(f"frequency number {fnum} not in {container.kind}")
    return int(idx[0])


# ---------------------------------------------------------------------------------------
# Cylindrical components as functions of r
# ---------------------------------------------------------------------------------------
def _exterior(pd: PointData, r: np.ndarray) -> Dict[str, np.ndarray]:
    """Exact exterior components at r >= R0, arrays (nr, nL load, nL obs)."""
    e1 = exterior_field(pd.kR, pd.kL, pd.phix, pd.phiz, pd.phiy, pd.alpha1, 1, r, pd.R0)
    e0 = exterior_field(pd.kR, pd.kL, pd.phix, pd.phiz, pd.phiy, pd.alpha0, 0, r, pd.R0)
    return {"u": e1[..., 0], "v": e1[..., 1], "w": e1[..., 2], "p": e0[..., 0], "q": e0[..., 2]}


def _table_grid(pd: PointData, rmax: float) -> np.ndarray:
    """Interpolation grid on [R0, rmax]: geometric (0.25 %) near R0, at most lambda_min/48 apart.

    ``lambda_min`` is the shortest wavelength of the propagating modes; the propagating sector is
    widened by the column's material loss (D-SIT-09, :func:`sassi.core.tlm.propagating_ratio`).
    Without any propagating mode every mode with Re k > 0 is resolved (conservative).
    """
    k = np.concatenate([pd.kR, pd.kL])
    loss = pd.loss if pd.loss is not None else spectrum_loss(k)
    kp = k[propagating(k, propagating_ratio(loss))]
    if kp.size == 0:
        kp = k[k.real > 0]
    lam_min = 2 * np.pi / np.max(kp.real) if kp.size else rmax
    dmax = lam_min / 48.0
    pts = [pd.R0]
    r = pd.R0
    while r < rmax:
        r = r + min(0.0025 * r, dmax)
        pts.append(min(r, rmax))
    return np.asarray(pts)


@quiet_fpe
def cylindrical_components(pd: PointData, r, max_exact: int = MAX_EXACT) -> Dict[str, np.ndarray]:
    """POINT values u~, v~, w~ (mu = 1) and p~, q~ (mu = 0) at distances ``r`` (any r >= 0).

    Returns arrays (nr, nL load, nL obs).  r = 0 gives the core-axis values, 0 < r < R0 the linear
    interpolation between the axis and R0 (R1 §4.1).  When more than ``max_exact`` distinct
    distances >= R0 are requested, the exact expansion is tabulated on a dense grid and
    interpolated by cubic splines of ``r * u(r)`` (R1 §3.6 tabulation option).
    """
    r = np.asarray(r, float).reshape(-1)
    nL = pd.load_iface.size
    out = {c: np.zeros((r.size, nL, nL), complex) for c in COMPONENTS}
    tol = 1e-9 * max(pd.R0, 1.0)
    ext = r >= pd.R0 * (1 - 1e-12)
    inner = (~ext) & (r > tol)
    axis = r <= tol
    ax = {"u": pd.axis1[..., 0], "v": pd.axis1[..., 1], "w": pd.axis1[..., 2],
          "p": pd.axis0[..., 0], "q": pd.axis0[..., 2]}
    if np.any(ext):
        re = np.maximum(r[ext], pd.R0)
        if re.size <= max_exact:
            vals = _exterior(pd, re)
        else:
            grid = _table_grid(pd, re.max())
            g = _exterior(pd, grid)
            vals = {}
            for c in COMPONENTS:
                sp = CubicSpline(grid, g[c] * grid[:, None, None], axis=0)
                vals[c] = sp(re) / re[:, None, None]
        for c in COMPONENTS:
            out[c][ext] = vals[c]
    if np.any(inner):
        rim = _exterior(pd, np.array([pd.R0]))
        t = (r[inner] / pd.R0)[:, None, None]
        for c in COMPONENTS:
            out[c][inner] = ax[c][None] + t * (rim[c] - ax[c][None])
    if np.any(axis):
        for c in COMPONENTS:
            out[c][axis] = ax[c][None]
    return out


# ---------------------------------------------------------------------------------------
# Binding API
# ---------------------------------------------------------------------------------------
def _iface_positions(pd: PointData, iface) -> np.ndarray:
    iface = np.asarray(iface, int).reshape(-1)
    lookup = {int(v): i for i, v in enumerate(pd.load_iface)}
    try:
        return np.array([lookup[int(v)] for v in iface], int)
    except KeyError as exc:
        raise ValueError(f"interaction node on interface {exc.args[0]} but FILE3 holds load interfaces "
                         f"{pd.load_iface.min()}..{pd.load_iface.max()} -- increase POINT <layer>") from None


@quiet_fpe
def flexibility_matrix(file3: Container, q: int, xy, iface, symmetrize: bool = True,
                       max_exact: int = MAX_EXACT, block: Optional[int] = None) -> np.ndarray:
    """Free-field flexibility matrix ``F_ff`` (3n x 3n, node-major ``[ux1,uy1,uz1,ux2,...]``).

    ``file3``: FILE3 container; ``q``: 0-based frequency row of FILE3 (see :func:`frequency_row`);
    ``xy`` (n, 2): horizontal coordinates of the interaction nodes; ``iface`` (n,): 1-based user
    interface of each node.  Entry ``F[3i+a, 3j+b]`` is the displacement of node i in direction a
    due to a unit force on node j in direction b.  Symmetrised ``(F + F^T)/2`` (D-ANL-02) unless
    ``symmetrize=False``.  POINT2 files (dim = 1) give the plane-strain blocks (requirements §4.6
    item 2): u_x, u_z from the P-SV solution (odd coupling terms times sgn(x_i - x_j)) and u_y from
    the SH solution; the y coordinate is ignored.
    """
    pd = point_data(file3, q)
    xy = np.asarray(xy, float).reshape(-1, 2)
    n = xy.shape[0]
    if n == 0:
        return np.zeros((0, 0), complex)
    pos = _iface_positions(pd, iface)
    if pos.size != n:
        raise ValueError("xy and iface have different lengths")
    if pd.dim == 1:
        from .strip2d import flexibility_matrix_2d
        return flexibility_matrix_2d(file3, q, xy, pos, symmetrize=symmetrize)
    F = _block_3d(pd, xy, pos, xy, pos, max_exact, block)
    if symmetrize:
        F += F.T
        F *= 0.5
    return F


@quiet_fpe
def flexibility_block(file3: Container, q: int, xy_obs, iface_obs, xy_load, iface_load,
                      max_exact: int = MAX_EXACT, block: Optional[int] = None) -> np.ndarray:
    """Flexibility between two point sets (3 n_obs x 3 n_load, node-major ``[ux, uy, uz]``): entry
    ``[3i+a, 3j+b]`` is the free-field displacement of observation point i (on interface
    ``iface_obs[i]``) in direction a due to a unit force at load point j (on ``iface_load[j]``) in
    direction b.  Not symmetrised.  With ``obs == load`` it equals ``flexibility_matrix(...,
    symmetrize=False)``; the SYMM image sums of ANALYS use it with the load points mirrored in the
    symmetry planes (:func:`sassi.core.symmetry.reduced_flexibility`).  POINT2 files give the
    plane-strain blocks (the y coordinates are ignored)."""
    pd = point_data(file3, q)
    xo = np.asarray(xy_obs, float).reshape(-1, 2)
    xl = np.asarray(xy_load, float).reshape(-1, 2)
    if xo.shape[0] == 0 or xl.shape[0] == 0:
        return np.zeros((3 * xo.shape[0], 3 * xl.shape[0]), complex)
    po, pl = _iface_positions(pd, iface_obs), _iface_positions(pd, iface_load)
    if po.size != xo.shape[0] or pl.size != xl.shape[0]:
        raise ValueError("coordinates and interfaces have different lengths")
    if pd.dim == 1:
        from .strip2d import flexibility_block_2d
        return flexibility_block_2d(file3, q, xo[:, 0], po, xl[:, 0], pl)
    return _block_3d(pd, xo, po, xl, pl, max_exact, block)


def _block_3d(pd: PointData, xo: np.ndarray, po: np.ndarray, xl: np.ndarray, pl: np.ndarray,
              max_exact: int, block: Optional[int]) -> np.ndarray:
    """3D flexibility of observation points ``xo`` (interface rows ``po``) under unit loads at ``xl``
    (``pl``): the rotated POINT3 blocks of the module docstring (not symmetrised)."""
    no, nl = xo.shape[0], xl.shape[0]
    dx = xo[:, 0][:, None] - xl[:, 0][None, :]          # x_i - x_j (obs - load)
    dy = xo[:, 1][:, None] - xl[:, 1][None, :]
    r = np.hypot(dx, dy)
    scale = max(float(r.max()), pd.R0)
    key = np.round(r / scale * 1e10).astype(np.int64)    # distances equal to 1e-10 relative
    ukey, inv = np.unique(key, return_inverse=True)
    inv = inv.reshape(no, nl)
    ru = ukey * (scale * 1e-10)
    comp = cylindrical_components(pd, ru, max_exact=max_exact)
    F = np.zeros((no, 3, nl, 3), complex)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(r > 0, dx / np.where(r > 0, r, 1.0), 1.0)
        s = np.where(r > 0, dy / np.where(r > 0, r, 1.0), 0.0)
    step = block or max(1, int(2_000_000 // max(nl, 1)))
    jload = pl[None, :]
    for i0 in range(0, no, step):
        i1 = min(no, i0 + step)
        rid = inv[i0:i1]
        iobs = po[i0:i1, None]
        u = comp["u"][rid, jload, iobs]
        v = comp["v"][rid, jload, iobs]
        w = comp["w"][rid, jload, iobs]
        p = comp["p"][rid, jload, iobs]
        qq = comp["q"][rid, jload, iobs]
        cc, ss = c[i0:i1], s[i0:i1]
        F[i0:i1, 0, :, 0] = u * cc * cc + v * ss * ss
        F[i0:i1, 0, :, 1] = (u - v) * ss * cc
        F[i0:i1, 0, :, 2] = p * cc
        F[i0:i1, 1, :, 0] = (u - v) * ss * cc
        F[i0:i1, 1, :, 1] = u * ss * ss + v * cc * cc
        F[i0:i1, 1, :, 2] = p * ss
        F[i0:i1, 2, :, 0] = w * cc
        F[i0:i1, 2, :, 1] = w * ss
        F[i0:i1, 2, :, 2] = qq
    F = F.reshape(3 * no, 3 * nl)
    if not np.all(np.isfinite(F)):
        raise FloatingPointError("flexibility matrix contains non-finite values")
    return F
