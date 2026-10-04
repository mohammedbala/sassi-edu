"""Thick shell commands THSHLSTR and THSHLSMH (manual section 9.19; requirements 3.4.Q; spec 11 section 4,
decisions D-TSH-01; TSHELL element :mod:`sassi.elements.tshell`).

``THSHLSTR,<flag>`` -- output of the TSHELL face stresses and strains by STRESS:

* ``0`` (default): STRESS writes only the 8 basic element components NXX NYY NXY QXZ QYZ MXX MYY MXY
  (local element axes; N and Q per unit length, M per unit length);
* ``1``: STRESS also computes the top/bottom **face** stresses and strains from the **maximum values**
  of the basic components (not from the time histories): ``sigma = N/t +- 6 M/t^2`` with the four sign
  permutations ``++ -- +- -+`` of (NXX, MXX) and (NYY, MYY), the in-plane shear with ``++`` only, the
  principal stresses and the plane-stress strains of each permutation, and the mid-surface transverse
  shear stresses ``1.5 Q/t`` (D-TSH-01; listing and ``TSHELL_FACE_STRESSES.TXT``).

The flag is a record of the model (WRITE writes it back, INP restores it).  The STRESS deck has no
``thshlstr`` parameter yet (the deck schema is lead-owned; requested in the work-package report), so the
command also writes the one-line file ``THSHLSTR.opt`` (``THSHLSTR,<flag>``) in the model directory,
which STRESS reads (:func:`sassi.modules.stress.thshlstr_flag`).  CHECK -- and therefore AFWRITE, which
writes the STRESS deck -- rewrites that file from the model record when it is missing or differs
(:func:`sync_thshlstr_file`), so a THSHLSTR given before MDL, or before MDL to another directory, still
reaches STRESS once AFWRITE has run in the model directory.

``THSHLSMH,<passes>,[type],[workdir]`` -- smoothing ("smear") of the TSHELL transverse shear forces
(manual 9.19.1, **NON VALIDATED in this version**: the command runs and prints that warning):

* input: ``TSHELL_ELEMENT_MAX.TXT``, the maxima of the basic components of the requested TSHELL
  elements written by STRESS in the model directory;
* neighbours: elements of the **same group** that share a node are ring 1, their neighbours ring 2, ...
  up to ``<passes>`` rings (spec 11 OQ-24: within a group, so that walls and slabs are not mixed);
* ``type`` 0 (default) Parzen-weighted average with ``M = passes + 1``:
  ``w(k) = 1 - 6 (k/M)^2 + 6 (k/M)^3`` for ``k <= M/2``, ``2 (1 - k/M)^3`` for ``M/2 < k <= M``;
  ``type`` 1 plain average (``w = 1``); ``Q_e' = sum w(k) Q_k / sum w(k)`` over the rings (QXZ and QYZ
  separately; they are in each element's local axes, so a group should share one orientation);
* output: the frame ``TSHELL_SMOOTHED_SHEAR.txt`` in ``workdir`` (default: the model directory), D-FIL-03
  header ``nrows ncols`` then rows ``group element QXZ QYZ`` (smoothed maxima).
"""
from __future__ import annotations

from collections import deque
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ...model import make_record
from ..check import TSHELL
from ..registry import command

#: sidecar read by STRESS (the deck schema has no thshlstr parameter)
THSHLSTR_FILE = "THSHLSTR.opt"
#: STRESS output read by THSHLSMH, and the THSHLSMH output frame
TSHELL_MAX_FILE = "TSHELL_ELEMENT_MAX.TXT"
SMOOTHED_FILE = "TSHELL_SMOOTHED_SHEAR.txt"
TYPES = {0: "Parzen-weighted average", 1: "plain average"}


# ======================================================================================
# THSHLSTR
# ======================================================================================
def write_thshlstr_file(model_dir, flag: int) -> Path:
    """``<model dir>/THSHLSTR.opt`` with the line ``THSHLSTR,<flag>`` (read by STRESS)."""
    p = Path(model_dir) / THSHLSTR_FILE
    p.write_text(f"THSHLSTR,{int(flag)}\n", encoding="utf-8")
    return p


def read_thshlstr_file(model_dir) -> Optional[int]:
    """Flag stored in ``<model dir>/THSHLSTR.opt`` (as :func:`sassi.modules.stress.thshlstr_flag` reads it), or
    None when the file is absent."""
    from ...modules.stress import thshlstr_flag
    d = Path(model_dir)
    if not d.is_dir() or not any(p.name.upper() == THSHLSTR_FILE.upper() for p in d.iterdir()):
        return None
    return thshlstr_flag(None, d)[0]


def model_thshlstr(model) -> int:
    """THSHLSTR flag of the model record (0, the manual default, without a record)."""
    rec = model.options.record("THSHLSTR")
    try:
        return 1 if rec is not None and int(float(rec.values[0])) == 1 else 0
    except (IndexError, TypeError, ValueError):
        return 0


def sync_thshlstr_file(model, model_dir=None) -> Optional[str]:
    """Bring ``THSHLSTR.opt`` in the model directory in line with the THSHLSTR record of the model.

    STRESS reads the flag from that file (the STRESS deck has no ``thshlstr`` parameter yet), which the
    THSHLSTR command writes when it runs.  After MDL to another directory, or INP of a ``.pre`` (WRITE does
    not write MDL) before MDL, the file is missing or stale; CHECK -- and so AFWRITE -- calls this function
    so that the STRESS run uses the flag of the model record.  The file is written when the record says 1
    and the file does not, or when an existing file differs from the record (no record = 0).  Returns a
    description of what was done, or None."""
    d = Path(model_dir) if model_dir else (Path(model.path) if getattr(model, "path", None) else None)
    flag = model_thshlstr(model)
    if d is None or not d.is_dir():
        return None
    current = read_thshlstr_file(d)
    if current == flag or (current is None and flag == 0):
        return None
    write_thshlstr_file(d, flag)
    was = "missing" if current is None else f"THSHLSTR,{current}"
    return f"{THSHLSTR_FILE} written in {d} from the model record (THSHLSTR,{flag}; the file was {was})"


@command("THSHLSTR", max_args=1)
def cmd_thshlstr(c):
    """THSHLSTR,<flag>: TSHELL face stresses/strains in STRESS (0 the 8 basic components only, 1 also the
    top/bottom face values from their maxima, D-TSH-01)."""
    flag = c.int(1, default=0, what="flag")
    if flag not in (0, 1):
        c.fail("<flag> must be 0 (only the 8 basic components) or 1 (also the face stresses and strains)")
    m = c.model
    m.options.set_record(make_record("THSHLSTR", [str(flag)]))
    what = ("STRESS also writes the TSHELL face stresses and strains (from the maxima of the basic components)"
            if flag else "STRESS writes only the 8 basic TSHELL components")
    if m.path and Path(m.path).is_dir():
        p = write_thshlstr_file(m.path, flag)
        c.confirm(f"THSHLSTR {flag}: {what} ({p.name} written in the model directory)")
    else:
        c.confirm(f"THSHLSTR {flag}: {what}")
        c.info("THSHLSTR: the model path is not defined yet (MDL): repeat THSHLSTR after MDL so that STRESS "
               f"finds {THSHLSTR_FILE} in the model directory")


# ======================================================================================
# THSHLSMH
# ======================================================================================
def parzen_weights(passes: int) -> np.ndarray:
    """Weights w(k), k = 0..passes, of the Parzen window with M = passes + 1 (spec 11 section 4.2)."""
    M = float(int(passes) + 1)
    k = np.arange(int(passes) + 1, dtype=float)
    r = k / M
    return np.where(r <= 0.5, 1.0 - 6.0 * r ** 2 + 6.0 * r ** 3, 2.0 * (1.0 - r) ** 3)


def ring_distances(neighbours: Dict[int, Sequence[int]], start: int, passes: int) -> Dict[int, int]:
    """Ring number (0 = itself) of every element reachable from ``start`` in at most ``passes`` rings."""
    dist = {start: 0}
    q = deque([start])
    while q:
        e = q.popleft()
        if dist[e] >= passes:
            continue
        for f in neighbours.get(e, ()):
            if f not in dist:
                dist[f] = dist[e] + 1
                q.append(f)
    return dist


def smooth_shear(values: Dict[Tuple[int, int], np.ndarray], nodes: Dict[Tuple[int, int], Sequence[int]],
                 passes: int, typ: int = 0) -> Dict[Tuple[int, int], np.ndarray]:
    """Smoothed QXZ, QYZ of every element of ``values`` ({(group, element): (QXZ, QYZ)}): weighted mean
    over the rings of neighbours of the same group (elements sharing a node), ``typ`` 0 Parzen, 1 plain."""
    keys = sorted(values)
    by_node: Dict[Tuple[int, int], List[Tuple[int, int]]] = {}
    for k in keys:
        for n in set(int(v) for v in nodes.get(k, ()) if int(v) > 0):
            by_node.setdefault((k[0], n), []).append(k)          # same group only
    neighbours: Dict[Tuple[int, int], set] = {k: set() for k in keys}
    for lst in by_node.values():
        for a in lst:
            neighbours[a].update(b for b in lst if b != a)
    w = parzen_weights(passes) if int(typ) == 0 else np.ones(int(passes) + 1)
    out = {}
    for k in keys:
        dist = ring_distances(neighbours, k, int(passes))
        ww = np.array([w[d] for d in dist.values()])
        vv = np.array([values[e] for e in dist])
        out[k] = (ww[:, None] * vv).sum(axis=0) / ww.sum()
    return out


def read_tshell_max(path: Path) -> Tuple[List[str], Dict[Tuple[int, int], np.ndarray]]:
    """Components and ``{(group, element): values}`` of ``TSHELL_ELEMENT_MAX.TXT`` (written by STRESS)."""
    comps: List[str] = []
    rows: Dict[Tuple[int, int], np.ndarray] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            t = s[1:].split()
            if len(t) > 2 and t[0] == "group" and t[1] == "element":
                comps = t[2:]
            continue
        p = s.split()
        rows[(int(p[0]), int(p[1]))] = np.array([float(x) for x in p[2:]])
    return comps, rows


@command("THSHLSMH", max_args=3)
def cmd_thshlsmh(c):
    """THSHLSMH,<passes>,[type],[workdir]: smooth the TSHELL transverse shear forces of STRESS over rings of
    neighbouring elements (type 0 Parzen weights, 1 plain average; NON VALIDATED in this version)."""
    passes = c.int(1, required=True, what="passes")
    typ = c.int(2, default=0, what="type")
    if passes < 0:
        c.fail("<passes> must be >= 0 (number of adjacent element layers)")
    if typ not in TYPES:
        c.fail("<type> must be 0 (Parzen-weighted average) or 1 (average)")
    c.warn("THSHLSMH is NON VALIDATED in this version (manual 9.19.1)")
    m = c.model
    mdir = Path(m.path) if m.path else Path(c.interp.cwd)
    src = next((p for p in (mdir.iterdir() if mdir.is_dir() else ())
                if p.name.upper() == TSHELL_MAX_FILE.upper()), None)
    if src is None:
        c.fail(f"{TSHELL_MAX_FILE} not found in {mdir}: run STRESS with TSHELL element output requests first")
    comps, rows = read_tshell_max(src)
    try:
        jx, jy = comps.index("QXZ"), comps.index("QYZ")
    except ValueError:
        c.fail(f"{src.name} has no QXZ/QYZ columns")
    nodes: Dict[Tuple[int, int], Sequence[int]] = {}
    missing = []
    for (g, e) in rows:
        grp = m.groups.get(g)
        if grp is None or grp.type != TSHELL or e not in grp.elements:
            missing.append((g, e))
            continue
        nodes[(g, e)] = list(grp.elements[e].nodes)
    if missing:
        c.warn(f"{len(missing)} elements of {src.name} are not TSHELL elements of the active model and are not "
               f"smoothed: " + ", ".join(f"group {g} element {e}" for g, e in missing[:8]))
    vals = {k: rows[k][[jx, jy]] for k in nodes}
    if not vals:
        c.fail("no TSHELL element of the active model in " + src.name)
    sm = smooth_shear(vals, nodes, passes, typ)
    out_dir = Path(c.str(3, "")) if c.given(3) else mdir
    if not out_dir.is_absolute():
        out_dir = mdir / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    keys = sorted(sm)
    lines = [f"{len(keys)} 4"] + [f"{g} {e} {sm[(g, e)][0]:.10e} {sm[(g, e)][1]:.10e}" for g, e in keys]
    out = out_dir / SMOOTHED_FILE
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    change = max(float(np.max(np.abs(sm[k] - vals[k]))) for k in keys)
    c.confirm(f"THSHLSMH: QXZ, QYZ of {len(keys)} TSHELL elements smoothed over {passes} ring(s) "
              f"({TYPES[typ]}); largest change {change:.6g}; frame {out}")
