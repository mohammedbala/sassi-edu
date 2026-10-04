"""Control-motion (acceleration history) file reading shared by SOIL, MOTION, STRESS, RELDISP.

Formats (MOTION <fopt>, requirements section 4.8 item 1):

* ``fopt = 0``: the first value is the time step, then one acceleration value per
  record (free format: several values per line are accepted, Fortran ``D``
  exponents allowed);
* ``fopt = 1``: (time, acceleration) pairs, one pair per line; dt is taken from the
  first two times and must be constant (relative tolerance 1e-6).

SOIL uses ``header`` lines to skip and ``nrval`` values (no dt in the file).
Accelerations are in g.  Scaling: if ``mult != 0`` multiply by it, else scale so that
max|a| = ``maxval`` (exactly one must be non-zero; Errors 77/78 are CHECK's job).
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np


def _numbers(text: str) -> List[float]:
    vals: List[float] = []
    for tok in text.replace(",", " ").split():
        try:
            vals.append(float(tok.replace("D", "E").replace("d", "e")))
        except ValueError:
            continue
    return vals


def read_history(path: Union[str, Path], fopt: int = 0, rec1: int = 1, rec2: int = 0,
                 header: int = 0) -> Tuple[np.ndarray, Optional[float]]:
    """Return (acc, dt).  ``dt`` is None when the format carries no time step.

    ``rec1``/``rec2`` are 1-based record numbers after the format header; ``rec2 = 0``
    means the last record.
    """
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()[header:]
    dt: Optional[float] = None
    if fopt == 1:
        t, a = [], []
        for ln in lines:
            v = _numbers(ln)
            if len(v) >= 2:
                t.append(v[0])
                a.append(v[1])
        if len(t) >= 2:
            dts = np.diff(np.asarray(t))
            dt = float(dts[0])
            if np.any(np.abs(dts - dt) > 1e-6 * abs(dt)):
                raise ValueError(f"{Path(path).name}: time step is not constant")
        acc = np.asarray(a, dtype=float)
    else:
        vals = _numbers("\n".join(lines))
        if header == 0 and vals:
            dt, acc = float(vals[0]), np.asarray(vals[1:], dtype=float)
        else:
            acc = np.asarray(vals, dtype=float)
    n = len(acc)
    r1 = max(1, int(rec1 or 1))
    r2 = n if not rec2 or rec2 <= 0 else min(int(rec2), n)
    return acc[r1 - 1:r2].copy(), dt


def scale_history(acc: np.ndarray, mult: float = 0.0, maxval: float = 0.0) -> np.ndarray:
    """MOTION/SOIL scaling rule (``mult`` has priority when non-zero)."""
    if mult:
        return acc * mult
    if maxval:
        peak = float(np.max(np.abs(acc))) if acc.size else 0.0
        return acc * (maxval / peak) if peak > 0 else acc.copy()
    return acc.copy()


def read_soil_history(path: Union[str, Path], nrval: int, header: int = 0) -> np.ndarray:
    """SOIL input: skip ``header`` lines then read ``nrval`` free-format values."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()[header:]
    vals = _numbers("\n".join(lines))
    if nrval and len(vals) < nrval:
        raise ValueError(f"{Path(path).name}: {len(vals)} values found, {nrval} requested")
    return np.asarray(vals[:nrval] if nrval else vals, dtype=float)


def write_history(path: Union[str, Path], acc: np.ndarray, dt: float, fopt: int = 0) -> Path:
    """Write a history in MOTION format (fopt 0: dt then values; 1: pairs)."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        if fopt == 1:
            for i, a in enumerate(acc):
                fh.write(f"{i * dt:.6f} {a:.8e}\n")
        else:
            fh.write(f"{dt:.8g}\n")
            for a in acc:
                fh.write(f"{a:.8e}\n")
    return path
