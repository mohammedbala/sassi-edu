"""Text result files (decision D-FIL-02) shared by MOTION, RELDISP, STRESS, SOIL, EQUAKE and the UI.

* ``.TFU / .TFI / .TFD``: optional ``#`` header line(s), then rows ``f amp [phase_rad]``
  (phase present when the complex option is on).  Phase is ``angle(H)`` in radians.
* ``.ACC / .THD / .THS / .TH / .acc / .vel / .dis``: first line the time step ``dt``,
  then one value per line (re-usable as a THFILE with fopt 0).
* ``.RS / .RSO / .rso / .rsi``: rows ``f SA`` (Hz, g).
* ``.psd``: rows ``f PSD``; ``.fft``: rows ``f Re Im``.

Readers skip non-numeric lines.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

import numpy as np

PathLike = Union[str, Path]


def _rows(path: PathLike):
    out = []
    for ln in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s[0] in "#*":
            continue
        try:
            out.append([float(t.replace("D", "E")) for t in s.replace(",", " ").split()])
        except ValueError:
            continue
    return out


def write_tf(path: PathLike, f: np.ndarray, H: np.ndarray, complex_: bool = True, header: str = "") -> Path:
    """Write a transfer function (.TFU/.TFI/.TFD)."""
    path = Path(path)
    H = np.asarray(H)
    with open(path, "w", encoding="utf-8") as fh:
        if header:
            fh.write(f"# {header}\n")
        if complex_:
            for fi, h in zip(f, H):
                fh.write(f"{fi:.6f} {abs(h):.8e} {np.angle(h):.8e}\n")
        else:
            for fi, h in zip(f, H):
                fh.write(f"{fi:.6f} {abs(h):.8e}\n")
    return path


def read_tf(path: PathLike) -> Tuple[np.ndarray, np.ndarray, bool]:
    """Return (f, H, is_complex).  For amplitude-only files H is real (the amplitude)."""
    rows = _rows(path)
    if not rows:
        return np.zeros(0), np.zeros(0), False
    cplx = all(len(r) >= 3 for r in rows)
    f = np.asarray([r[0] for r in rows])
    amp = np.asarray([r[1] for r in rows])
    if cplx:
        ph = np.asarray([r[2] for r in rows])
        return f, amp * np.exp(1j * ph), True
    return f, amp, False


def write_history(path: PathLike, values: np.ndarray, dt: float) -> Path:
    """Write a time history: first line dt, then one value per line."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"{dt:.8g}\n")
        for v in np.asarray(values, dtype=float):
            fh.write(f"{v:.8e}\n")
    return path


def read_history(path: PathLike) -> Tuple[np.ndarray, float]:
    rows = _rows(path)
    if not rows:
        return np.zeros(0), 0.0
    flat = [v for r in rows for v in r]
    return np.asarray(flat[1:], dtype=float), float(flat[0])


def write_xy(path: PathLike, x: np.ndarray, *ys: np.ndarray, header: str = "", fmt: str = "{:.8e}") -> Path:
    """Write columns (x, y1, y2, ...), e.g. .RS (f, SA), .psd, .fft."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        if header:
            fh.write(f"# {header}\n")
        for i in range(len(x)):
            fh.write(" ".join([f"{x[i]:.6f}"] + [fmt.format(float(y[i])) for y in ys]) + "\n")
    return path


def read_xy(path: PathLike) -> np.ndarray:
    """Return a 2D array of the numeric rows (ragged rows truncated to the shortest)."""
    rows = _rows(path)
    if not rows:
        return np.zeros((0, 0))
    n = min(len(r) for r in rows)
    return np.asarray([r[:n] for r in rows], dtype=float)
