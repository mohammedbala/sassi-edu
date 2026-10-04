"""FFT helpers, frequency-domain convolution and baseline correction (MOTION, RELDISP, STRESS).

Conventions (requirements section 4.0.1, D-CNV-01/02):

* harmonic factor ``exp(+i w t)``; forward transform ``A_k = sum_n a_n exp(-2 pi i k n / N)``
  (``numpy.fft.rfft``, no scaling) and inverse with ``1/N`` (``numpy.fft.irfft``);
* Fourier frequencies ``f_k = k df``, ``k = 0 .. N/2``, ``df = 1/(N dt)``; the last bin is the
  Nyquist frequency ``1/(2 dt)``.

A transfer function computed by ANALYS under ``exp(+i w t)`` therefore applies directly:
the response is ``r = irfft(H * rfft(a))``.  The record is periodic with period ``N dt``; the
zero-padded tail (the *quiet zone*) must be long enough for the free vibration to decay,
otherwise it wraps around to the start of the record (requirements section 4.0.1, spec 05c
B.5.10).

Normalisation of seismic TFs (D-CNV-06, requirements section 4.0.5): FILE8 holds the
dimensionless ratio ``H = U/U_cp``, the same for total displacement and total acceleration.
Accelerations are convolved with the control-acceleration spectrum ``A`` (in g, giving
results in g), quantities linear in displacement with the control *displacement* spectrum
``U_g = -g A / w^2`` (the f = 0 term set to zero).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

__all__ = [
    "fourier_grid", "omega_grid", "pad_record", "forward", "convolve", "hermitian_bins",
    "displacement_spectrum", "derivative_factor", "output_length", "integrate_trapz",
    "BaselineResult", "baseline_correction", "hudson_baseline",
]


# ---------------------------------------------------------------------------------------
# Grids and transforms
# ---------------------------------------------------------------------------------------
def fourier_grid(nfft: int, delt: float) -> np.ndarray:
    """Fourier frequencies ``f_k = k/(nfft*delt)``, k = 0 .. nfft//2 (Hz), the rfft bins."""
    nfft = int(nfft)
    if nfft < 2 or delt <= 0:
        raise ValueError("need nfft >= 2 and delt > 0")
    return np.arange(nfft // 2 + 1) / (nfft * float(delt))


def omega_grid(nfft: int, delt: float) -> np.ndarray:
    """Circular frequencies ``w_k = 2 pi f_k`` of the rfft bins (rad/s)."""
    return 2.0 * np.pi * fourier_grid(nfft, delt)


def pad_record(a: np.ndarray, nfft: int) -> np.ndarray:
    """Zero-pad a record to ``nfft`` samples (requirements section 4.8 item 1).

    Raises ValueError when the record is longer than the Fourier period (the error MOTION
    reports as "time history longer than Fourier period").
    """
    a = np.asarray(a, dtype=float)
    if a.ndim != 1:
        raise ValueError("record must be one-dimensional")
    if len(a) > nfft:
        raise ValueError(f"time history has {len(a)} values, longer than the Fourier period NFFT = {nfft}")
    out = np.zeros(int(nfft))
    out[:len(a)] = a
    return out


def forward(a: np.ndarray, nfft: int) -> np.ndarray:
    """``rfft`` of the zero-padded record: the Fourier coefficients ``A_k`` (k = 0 .. nfft/2)."""
    return np.fft.rfft(pad_record(a, nfft))


def hermitian_bins(H_grid: np.ndarray, nfft: int) -> np.ndarray:
    """Make the Nyquist bin of a TF on the rfft grid real (D-MOT-03).

    A real time signal has a real Nyquist coefficient; ``irfft`` would silently drop the
    imaginary part, so we make that explicit.  The f = 0 bin of a seismic TF is the real
    rigid-body anchor; a vibration TF keeps H_1 there (its imaginary part is discarded by the
    inverse transform, which is again made explicit by taking the real part).
    """
    H = np.array(H_grid, dtype=complex, copy=True)
    if int(nfft) % 2 == 0 and H.shape[0] == int(nfft) // 2 + 1:
        H[-1] = H[-1].real
    H[0] = H[0].real
    return H


def convolve(H_grid: np.ndarray, A: np.ndarray, nfft: Optional[int] = None) -> np.ndarray:
    """Frequency-domain convolution ``r = irfft(H_k A_k)`` (requirements section 4.8 item 7).

    ``H_grid`` is (nK,) or (nK, m) on the rfft grid, ``A`` (nK,) the Fourier coefficients of the
    input; returns (nfft,) or (nfft, m).  ``nfft`` defaults to ``2 (nK - 1)``.  The DC and
    Nyquist bins of the product are taken real (the inverse of a real record).
    """
    H_grid = np.asarray(H_grid, dtype=complex)
    A = np.asarray(A, dtype=complex)
    nK = A.shape[0]
    if H_grid.shape[0] != nK:
        raise ValueError(f"H has {H_grid.shape[0]} bins, A has {nK}")
    n = int(nfft) if nfft is not None else 2 * (nK - 1)
    R = H_grid * (A if H_grid.ndim == 1 else A[:, None])
    R = hermitian_bins(R, n)
    return np.fft.irfft(R, n=n, axis=0)


def displacement_spectrum(A: np.ndarray, f: np.ndarray, gravity: float) -> np.ndarray:
    """Control displacement spectrum ``U_g = -g A / w^2`` with the f = 0 term zero (D-CNV-06).

    ``A`` is the Fourier transform of the control acceleration in g; the result is the Fourier
    transform of the ground displacement in model length units (``g`` = HOUSE gravity).  The
    zero-frequency term is undetermined (a constant displacement offset) and set to 0.
    """
    A = np.asarray(A, dtype=complex)
    w = 2.0 * np.pi * np.asarray(f, dtype=float)
    out = np.zeros_like(A)
    nz = w > 0
    out[nz] = -float(gravity) * A[nz] / (w[nz] ** 2)
    return out


def derivative_factor(f: np.ndarray, order: int) -> np.ndarray:
    """Factor ``(i w)^order`` that turns a displacement spectrum into velocity (1) or acceleration (2)."""
    w = 2.0 * np.pi * np.asarray(f, dtype=float)
    return (1j * w) ** int(order)


def output_length(nfft: int, delt: float, dur: float) -> int:
    """Number of output samples: ``min(nfft, round(1.2 dur/delt))``, all ``nfft`` when dur = 0.

    The 20 % extension shows part of the free vibration after the excitation (D-MOT-09).
    """
    if dur is None or dur <= 0:
        return int(nfft)
    return int(min(int(nfft), max(1, int(round(1.2 * float(dur) / float(delt))))))


# ---------------------------------------------------------------------------------------
# Integration and baseline correction
# ---------------------------------------------------------------------------------------
def integrate_trapz(acc: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    """Trapezoidal velocity and displacement starting from rest (v_0 = d_0 = 0)."""
    acc = np.asarray(acc, dtype=float)
    vel = np.concatenate([[0.0], np.cumsum(0.5 * (acc[1:] + acc[:-1]) * dt)])
    dis = np.concatenate([[0.0], np.cumsum(0.5 * (vel[1:] + vel[:-1]) * dt)])
    return vel, dis


@dataclass
class BaselineResult:
    """Output of the Hudson-Housner correction (D-MOT-08, spec 05c B.6.13).

    ``acc`` corrected acceleration (input units), ``vel`` and ``dis`` the re-integrated
    corrected velocity and displacement (input units times ``scale``), ``c1``/``c2`` the removed
    acceleration baseline ``c1 + c2 t`` (input units), ``final_disp`` = ``dis[-1]``.
    """
    acc: np.ndarray
    vel: np.ndarray
    dis: np.ndarray
    c1: float
    c2: float
    final_disp: float


def baseline_correction(acc: np.ndarray, dt: float, scale: float = 1.0) -> BaselineResult:
    """Classical Hudson-Housner time-domain baseline correction (D-MOT-08, spec 05c B.6.13).

    1. integrate the acceleration (trapezoidal rule) to the velocity ``v`` (v_0 = 0);
    2. least-squares fit ``v(t) ~ c1 t + c2 t^2/2`` -- the velocity produced by an acceleration
       baseline error ``c1 + c2 t`` (minimum mean-square velocity criterion);
    3. subtract the baseline from the acceleration and re-integrate to velocity and
       displacement.

    ``scale`` converts acceleration units before integration (e.g. g -> length/s^2 with the
    SSI gravity).  The correction only approximates absolute displacements: the user should
    check that the final displacement is about zero (RELDISP gives accurate relative
    displacements).
    """
    a = np.asarray(acc, dtype=float)
    n = len(a)
    t = np.arange(n) * float(dt)
    vel, _ = integrate_trapz(a, dt)
    if n < 3 or t[-1] <= 0:
        c1 = c2 = 0.0
    else:
        B = np.stack([t, 0.5 * t * t], axis=1)
        # scale the columns for conditioning, then solve the 2x2 least-squares problem
        sc = np.array([t[-1], 0.5 * t[-1] ** 2])
        coef, *_ = np.linalg.lstsq(B / sc, vel, rcond=None)
        c1, c2 = (coef / sc).tolist()
    ac = a - (c1 + c2 * t)
    vc, dc = integrate_trapz(ac, dt)
    s = float(scale)
    return BaselineResult(acc=ac, vel=vc * s, dis=dc * s, c1=float(c1), c2=float(c2),
                          final_disp=float(dc[-1] * s) if n else 0.0)


def hudson_baseline(acc: np.ndarray, dt: float) -> np.ndarray:
    """Baseline-corrected acceleration (binding API, ARCHITECTURE section 6.4; D-MOT-08)."""
    return baseline_correction(acc, dt).acc
