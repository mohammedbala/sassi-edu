"""Deformed shapes made from the analysis results (Plot > Deformed Shape, Load Frame Data).

An animation plays frames.  The modules write frames only when a restart option asks for them, and no
example does, so after running an example the Load Frame Data dialog had nothing to offer: a deformed
shape needed ``HARMFRAME`` and ``PROCFRAME`` by hand (the lessons' Animate buttons type them).  The dialog
now also offers the steady-state motion at one computed frequency from the transfer functions ANALYS
writes (``FILE8``), with the frequency of the largest deformation preselected.  This module lists those files
and turns a choice into command text (rule L17):

* :func:`sources` -- the FILE8-type files in the active model's folder, their computed frequencies and the
  frequency of the largest deformation (``GET /api/harmonic``);
* :func:`plan` -- ``HARMFRAME,<file>,<f>,<folder>,,[0]`` into a frame folder of its own, ``HARM_<f>``
  (``HARM_<f>R`` relative to the free field) (``POST /api/harmonic/plan``);
* :func:`show` -- after HARMFRAME has run: ``PROCFRAME`` of that folder and ``DEFORMPLOT`` with the
  automatic scale, as a lesson's Animate button (:func:`sassi.ui.learn.animate_lines`)
  (``POST /api/harmonic/show``).
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from ..prep.lexer import join_command
from . import learn

#: names of the FILE8-type files ANALYS and COMBIN write (FILE8, FILE8X, FILE81 ...)
FILE8_NAME = re.compile(r"FILE8[A-Z0-9_]*", re.IGNORECASE)
#: the frame folder of one choice: HARM_<f with p for the point>[R]
FOLDER_PREFIX = "HARM_"


class HarmonicError(ValueError):
    """A request that cannot be served (status: 400 bad request, 404 not found)."""

    def __init__(self, message: str, status: int = 400):
        super().__init__(message)
        self.status = status


def _folder(interp) -> Optional[Path]:
    m = interp.model
    if not m.path:
        return None
    p = Path(m.path)
    return p if p.is_dir() else None


def _summary(path: Path) -> Optional[Dict[str, Any]]:
    """Frequencies and the largest deformation of a FILE8-type file (None when it is not one).

    The preselected frequency is where the motion differs most from node to node: the largest
    ``|H_n - mean_n H|`` over the nodes, in the direction of the control motion (seismic) or in each
    translation (vibration).  The largest ``|H|`` would not do: an embedded box never moves more than the
    free field, and its largest total motion is the rigid translation at the lowest frequency; the spread
    picks the SSI resonance of a structure (3.49 Hz for example 1) and the rocking and scattering of a box."""
    from ..modules.base import ModuleError
    from ..modules.motion import read_file8
    if not zipfile.is_zipfile(path):
        return None
    try:
        f8 = read_file8(path, path.name)
    except ModuleError:
        return None
    meta = dict(f8.meta)
    fnum = np.asarray(f8["fnum"], dtype=np.int64)
    freqs = fnum * float(meta["df"])
    H = np.asarray(f8["H"])
    eq_node = np.asarray(f8["eq_node"], dtype=np.int64)
    eq_dof = np.asarray(f8["eq_dof"], dtype=np.int64)
    seismic = int(meta.get("type", 0)) == 0
    cm = int(meta.get("cm", 0))
    # the spread of the motion over the nodes, per direction: the control direction (seismic), else X, Y, Z
    dirs = [cm + 1] if seismic else [1, 2, 3]
    cols_all: List[int] = []
    spread: List[np.ndarray] = []
    for d in dirs:
        cols = np.nonzero(eq_dof == d)[0]
        if len(cols):
            Hd = H[:, cols]
            spread.append(np.abs(Hd - Hd.mean(axis=1, keepdims=True)))
            cols_all.extend(int(c) for c in cols)
    if not cols_all or not len(freqs):
        return None
    A = np.hstack(spread)
    q, c = np.unravel_index(int(np.argmax(A)), A.shape)
    col = cols_all[c]
    return {"file": path.name, "freqs": [float(f"{f:.6g}") for f in freqs], "peak": float(f"{freqs[q]:.6g}"),
            "peak_node": int(eq_node[col]), "peak_dir": "XYZ"[int(eq_dof[col]) - 1],
            "peak_amp": float(A[q, c]), "seismic": seismic, "direction": "XYZ"[cm] if seismic else "",
            "title": str(meta.get("title", "")), "content": str(meta.get("content", ""))}


def sources(interp) -> Dict[str, Any]:
    """The FILE8-type files of the active model's folder (``GET /api/harmonic``)."""
    folder = _folder(interp)
    out: List[Dict[str, Any]] = []
    if folder is not None:
        for p in sorted(folder.iterdir(), key=lambda x: (x.name.upper() != "FILE8", x.name.upper())):
            if p.is_file() and FILE8_NAME.fullmatch(p.name):
                s = _summary(p)
                if s is not None:
                    out.append(s)
    return {"model": interp.model.name or "", "folder": str(folder) if folder else "", "sources": out}


def folder_name(freq: float, relative: bool) -> str:
    """``HARM_3p491`` / ``HARM_3p491R``: the frame folder of a frequency (no point in the name)."""
    return FOLDER_PREFIX + f"{float(freq):.4g}".replace(".", "p").replace("-", "m") + ("R" if relative else "")


def plan(interp, file: str, freq: float, relative: bool = False) -> Dict[str, Any]:
    """``HARMFRAME`` of one choice (``POST /api/harmonic/plan``): the line and the frame folder."""
    folder = _folder(interp)
    name = (file or "").strip()
    if folder is None:
        raise HarmonicError("the active model has no folder yet: run its analysis first", 404)
    if not FILE8_NAME.fullmatch(name) or not (folder / name).is_file():
        raise HarmonicError(f"{name or '(none)'}: no FILE8-type file of that name in {folder}", 404)
    try:
        f = float(freq)
    except (TypeError, ValueError):
        raise HarmonicError(f"frequency {freq!r} is not a number") from None
    if not f > 0:
        raise HarmonicError("the frequency must be > 0 Hz")
    out = folder_name(f, relative)
    args: List[Any] = [name, float(f"{f:.6g}"), out]
    if relative:
        args += ["", 0]
    return {"lines": [join_command("HARMFRAME", args)], "folder": out}


def show(interp, folder: str, title: str = "") -> Dict[str, Any]:
    """``PROCFRAME`` + ``DEFORMPLOT`` of a frame folder HARMFRAME wrote (``POST /api/harmonic/show``)."""
    base = _folder(interp)
    if base is None:
        raise HarmonicError("the active model has no folder", 404)
    name = (folder or "").strip()
    if not re.fullmatch(FOLDER_PREFIX + r"[0-9A-Za-z]+", name):
        raise HarmonicError(f"{name!r} is not a frame folder of the Deformed Shape dialog")
    spec = f"{name} | deformed" + (f" | {title.replace('|', '/').strip()}" if title and title.strip() else "")
    try:
        lines = learn.animate_lines(spec, base, interp)
    except learn.LearnError as exc:
        raise HarmonicError(str(exc), exc.status) from None
    return {"lines": lines}
