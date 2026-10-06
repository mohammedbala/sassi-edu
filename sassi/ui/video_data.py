"""Curves for the explainer videos (``sassi/ui/static/videos``): run a lesson headlessly and print the
result files a video plots, downsampled, as JSON to paste into its script.

Usage (from the project root)::

    python -m sassi.ui.video_data 01-why-ssi ex01/00085TR_X.TFI ex01/00085TR_X02.RS ex01/00041TR_X.ACC
    python -m sassi.ui.video_data 01-why-ssi --list ex01          # the files a run writes in ex01/
    python -m sassi.ui.video_data 01-why-ssi --keep DIR ...        # keep the workspace in DIR
    python -m sassi.ui.video_data 01-why-ssi --from DIR ...        # read a kept workspace, no run
    python -m sassi.ui.video_data 01-why-ssi --js 01 FILES > sassi/ui/static/videos/d01.js   # the data file

The lesson runs like a learner pressing every "Run step" (:func:`sassi.ui.lessons.run_lesson_headless`) in a
temporary folder that is deleted afterwards (the runs write 20-200 MB).  Paths are relative to the lesson
workspace.  What is printed per file (numbers rounded to 4 significant digits):

* ``.TFI`` / ``.TFU`` / ``.TFD`` (transfer functions): ``{"f": [...], "amp": [...], "ph": [...]}``; a ``.TFI``
  is thinned to at most ``--points`` frequencies (every k-th, keeping the peak);
* ``.RS`` (response spectra): ``{"f": [...], "sa": [...]}`` thinned the same way;
* ``.ACC`` / ``.THD`` / ``.THS`` and other one-column histories (dt first): ``{"dt": ..., "t": [...], "v": [...]}``
  with a min/max envelope per bin, so the peaks survive the thinning;
* anything else with numeric rows: ``{"cols": [[...], ...]}`` (columns, thinned; the peak of every column kept).
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Dict, List

import numpy as np

from . import lessons as L


def _r(v: float) -> float:
    return float(f"{float(v):.4g}")


def _thin_xy(x: np.ndarray, ys: List[np.ndarray], n: int) -> List[np.ndarray]:
    """Every k-th point, the last one and the largest |value| of every y column, so that ~n remain."""
    m = len(x)
    if m <= n:
        return [x] + ys
    k = int(np.ceil(m / n))
    idx = set(range(0, m, k)) | {m - 1} | {int(np.argmax(np.abs(y))) for y in ys}
    # local peaks and troughs of the first column, with their neighbours, so that peak shapes survive
    y0 = ys[0]
    ext = [i for i in range(1, m - 1) if (y0[i] - y0[i - 1]) * (y0[i + 1] - y0[i]) < 0]
    ext.sort(key=lambda i: -abs(y0[i]))
    for i in ext[: max(4, n // 4)]:
        idx |= {i - 1, i, i + 1}
    sel = np.array(sorted(idx))
    return [x[sel]] + [y[sel] for y in ys]


def _rows(path: Path) -> List[List[float]]:
    rows = []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        try:
            rows.append([float(t) for t in s.replace(",", " ").split()])
        except ValueError:
            continue
    return rows


def extract(path: Path, points: int = 240) -> Dict:
    """The JSON-ready content of one result file (see the module docstring)."""
    suf = path.suffix.upper()
    rows = _rows(path)
    if not rows:
        return {"error": "no numeric rows"}
    if suf in (".TFI", ".TFU", ".TFD"):
        A = np.asarray([r[:3] for r in rows if len(r) >= 3])
        f, amp, ph = _thin_xy(A[:, 0], [A[:, 1], A[:, 2]], points)
        return {"f": [_r(v) for v in f], "amp": [_r(v) for v in amp], "ph": [_r(v) for v in ph]}
    if suf == ".RS":
        A = np.asarray([r[:2] for r in rows if len(r) >= 2])
        f, sa = _thin_xy(A[:, 0], [A[:, 1]], points)
        return {"f": [_r(v) for v in f], "sa": [_r(v) for v in sa]}
    if all(len(r) == 1 for r in rows):
        flat = [r[0] for r in rows]
        dt, v = flat[0], np.asarray(flat[1:])
        bins = max(1, points // 2)
        k = max(1, int(np.ceil(len(v) / bins)))
        t_out, v_out = [], []
        for i0 in range(0, len(v), k):
            seg = v[i0:i0 + k]
            a, b = int(np.argmin(seg)), int(np.argmax(seg))
            for j in sorted((a, b)) if a != b else (a,):
                t_out.append(_r((i0 + j) * dt))
                v_out.append(_r(seg[j]))
        return {"dt": dt, "n": len(v), "peak": _r(np.max(np.abs(v))), "t": t_out, "v": v_out}
    n = min(len(r) for r in rows)
    A = np.asarray([r[:n] for r in rows])
    cols = _thin_xy(A[:, 0], [A[:, j] for j in range(1, n)], points)
    return {"cols": [[_r(v) for v in c] for c in cols]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="run a lesson headlessly and print result curves for an explainer video")
    ap.add_argument("lesson", help="lesson id, e.g. 01-why-ssi")
    ap.add_argument("files", nargs="*", help="result files, relative to the lesson workspace (e.g. ex01/00085TR_X.TFI)")
    ap.add_argument("--points", type=int, default=240, help="at most about this many points per curve (default 240)")
    ap.add_argument("--list", metavar="DIR", help="list the files of this workspace folder instead")
    ap.add_argument("--keep", metavar="DIR", help="run in DIR and keep it (default: a temporary folder, deleted)")
    ap.add_argument("--from", dest="ws", metavar="DIR", help="read the files of a workspace kept with --keep (no run)")
    ap.add_argument("--js", metavar="KEY", help="print a data script for the videos: SV.DATA[KEY] = {...}")
    a = ap.parse_intermixed_args(argv)
    les = L.lesson_by_id(a.lesson)
    if a.ws:
        root = Path(a.ws)
        a.keep = a.ws
    else:
        root = Path(a.keep) if a.keep else Path(tempfile.mkdtemp(prefix="sassi-video-"))
    try:
        if a.ws:
            ws = root / les.id if (root / les.id).is_dir() else root
        else:
            rep = L.run_lesson_headless(les, root)
            ws = Path(rep.workspace)
            if rep.errors:
                print(f"warning: the lesson run reported {len(rep.errors)} error(s): {rep.errors[:3]}", file=sys.stderr)
        if a.list:
            for p in sorted((ws / a.list).iterdir()):
                print(p.name)
            return 0
        out = {}
        for rel in a.files:
            p = ws / rel
            out[rel] = extract(p, a.points) if p.is_file() else {"error": "not found"}
        text = json.dumps(out, separators=(",", ":"))
        if a.js:
            text = text.replace(',"', ',\n"')
            av = list(argv if argv is not None else sys.argv[1:])
            for flag in ("--from", "--keep"):          # the printed command reruns the lesson instead
                while flag in av:
                    i = av.index(flag)
                    del av[i:i + 2]
            cmd = " ".join(["python -m sassi.ui.video_data"] + av)
            print(f"/* Curves of lesson {a.lesson} for explainer video {a.js} (generated; do not edit):\n *   {cmd}\n */")
            print(f'"use strict";\n(window.SV = window.SV || {{}}).DATA = window.SV.DATA || {{}};\nSV.DATA[{json.dumps(a.js)}] = {text};')
        else:
            print(text)
        return 0
    finally:
        if not a.keep:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
