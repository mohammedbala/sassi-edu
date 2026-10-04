"""Model-directory files for the GUI: listing with result types, text read/write, path safety.

The Results browser, the Line Selection dialog (Plot > Spectrum TFU-TFI / Time History), the File
Editor (File > Open) and the Check Errors window read files through this module.

**Path safety.**  The server only reads and writes files inside the *allowed roots*: the model
directories (MDL paths) of the models in memory and the interpreter's working directory.  A name
is resolved (symbolic links followed) before the check, so ``..`` components or links cannot leave
the roots (:func:`safe_path`).  Nothing here executes a file.

**Result types** follow requirements section 1.7 (files) and D-FIL-02 (text formats):

=====================  ==========  =================================================================
pattern                plot        content
=====================  ==========  =================================================================
``.TFU .TFI .TFD``     spectrum    ``f amp [phase]`` transfer functions (MOTION, STRESS, RELDISP)
``.RS .RSO .rsi .psd`` spectrum    ``f value`` response spectra, target spectra, PSD
``.fft``               spectrum    ``f Re Im``
``.ACC .THD .THS .TH`` history     first line dt, then one value per line (READTH Pair 0)
``.acc .vel .dis``     history     EQUAKE histories (same format)
=====================  ==========  =================================================================
"""
from __future__ import annotations

import os
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

MAX_TEXT_BYTES = 8 * 1024 * 1024          # File Editor / viewer limit
SPEC_EXT = {".tfu": "Transfer function (computed)", ".tfi": "Transfer function (interpolated)",
            ".tfd": "Relative-displacement TF", ".rs": "Response spectrum", ".rso": "Response spectrum (output)",
            ".rsi": "Target response spectrum", ".psd": "Power spectral density", ".fft": "Fourier spectrum"}
TH_EXT = {".acc": "Acceleration history", ".thd": "Relative-displacement history", ".ths": "Stress history",
          ".th": "Soil layer history", ".vel": "Velocity history", ".dis": "Displacement history"}
DECK_EXT = {".equ": "EQUAKE deck", ".soi": "SOIL deck", ".sit": "SITE deck", ".poi": "POINT deck",
            ".hou": "HOUSE deck", ".frc": "FORCE deck", ".anl": "ANALYS deck", ".mot": "MOTION deck",
            ".str": "STRESS deck", ".rdi": "RELDISP deck", ".eql": "NONLINEAR deck",
            ".lgn": "LOADGEN deck"}
OTHER_EXT = {".pre": "Command file (.pre)", ".err": "CHECK errors", ".out": "Module listing", ".sdb": "Model database",
             ".png": "Image", ".bmp": "Image", ".csv": "Table", ".txt": "Text", ".inp": "Input file",
             ".pin": "Nonlinear soil input", ".map": "Node map", ".bak": "Backup"}
BINARY_FILES = re.compile(r"^(FILE\d+[XYZ]?\d*|COOS[KM]|DOFSMAP|COOX\w*|COOTK\w*|.*\.N4|.*\.n4)$")
TEXT_FILES = re.compile(r"^(FILE1[2-5]|FILE7[34]|FILE78|FILE88|FOUN(STIF|DASH|DAMP|IMPD)|SRSSTF\.txt|CONTTRS\.txt|"
                        r"Frames\.txt|COMBIN\.opt)$", re.I)


class PathError(ValueError):
    """A file name outside the allowed roots or otherwise not acceptable."""


def allowed_roots(interp, extra: Sequence[Path] = ()) -> List[Path]:
    """Model directories of the models in memory (active model first), the working directory and
    ``extra`` directories (the GUI's start directory, ``sassi-gui --model-dir``)."""
    roots: List[Path] = []
    order = [interp.active_model] + [n for n in sorted(interp.models) if n != interp.active_model]
    for n in order:
        p = interp.models[n].path
        if p:
            rp = Path(p).expanduser().resolve()
            if rp not in roots:
                roots.append(rp)
    for d in [Path(interp.cwd)] + [Path(x) for x in extra]:
        rd = d.expanduser().resolve()
        if rd not in roots:
            roots.append(rd)
    return roots


def _inside(p: Path, roots: Sequence[Path]) -> bool:
    for r in roots:
        try:
            p.relative_to(r)
            return True
        except ValueError:
            continue
    return False


def safe_path(name: str, roots: Sequence[Path], must_exist: bool = False) -> Path:
    """Resolve ``name`` (absolute, or relative to the first root containing it) inside ``roots``.

    Relative names are tried against every root in order; for a file to be created the first root
    is used.  The resolved path (links followed) must lie inside one of the roots, otherwise
    :class:`PathError` -- this rejects ``../`` traversal and links pointing elsewhere.
    """
    if name is None or not str(name).strip():
        raise PathError("file name required")
    raw = os.path.expanduser(str(name).strip())
    if "\x00" in raw:
        raise PathError("invalid file name")
    if "\\" in raw and os.sep == "/":
        raw = raw.replace("\\", "/")
    p = Path(raw)
    refused = PathError(f"{name}: outside the model directory and the working directory (refused)")
    if p.is_absolute():
        rc = p.resolve()
        if not _inside(rc, roots):
            raise refused
        if must_exist and not rc.exists():
            raise PathError(f"{name}: file not found")
        return rc
    outside = False
    for r in roots:
        rc = (r / p).resolve()
        if not _inside(rc, roots):
            outside = True
            continue
        if rc.exists():
            return rc
    if outside:
        raise refused
    if must_exist:
        raise PathError(f"{name}: file not found")
    return (roots[0] / p).resolve()             # a new file: in the first root (model directory)


def classify(path: Path) -> Dict[str, Any]:
    """Result type of a file: ``kind`` (spectrum / history / deck / binary / text / image / other),
    ``label`` and ``plot`` ('spec', 'th' or None)."""
    name = path.name
    ext = path.suffix.lower()
    if ext in SPEC_EXT:
        return {"kind": "spectrum", "label": SPEC_EXT[ext], "plot": "spec", "type": ext[1:].upper()}
    if ext in TH_EXT:
        return {"kind": "history", "label": TH_EXT[ext], "plot": "th", "type": ext[1:].upper()}
    if ext in DECK_EXT:
        return {"kind": "deck", "label": DECK_EXT[ext], "plot": None, "type": ext[1:].upper()}
    if BINARY_FILES.match(name):
        return {"kind": "binary", "label": "Binary inter-module file", "plot": None, "type": "BIN"}
    if TEXT_FILES.match(name):
        return {"kind": "text", "label": "Module text file", "plot": None, "type": "TXT"}
    if ext in OTHER_EXT:
        kind = "image" if ext in (".png", ".bmp") else ("binary" if ext == ".sdb" else "text")
        return {"kind": kind, "label": OTHER_EXT[ext], "plot": None, "type": ext[1:].upper()}
    return {"kind": "other", "label": "", "plot": None, "type": ext[1:].upper() if ext else ""}


def listing(directory: Path, roots: Sequence[Path]) -> Dict[str, Any]:
    """Directory listing for the Results browser: files with size, time and result type, sub-dirs."""
    d = safe_path(str(directory), roots, must_exist=True)
    if not d.is_dir():
        raise PathError(f"{directory} is not a directory")
    files, dirs = [], []
    for p in sorted(d.iterdir(), key=lambda q: q.name.lower()):
        try:
            st = p.stat()
        except OSError:
            continue
        if p.is_dir():
            dirs.append({"name": p.name, "path": str(p)})
            continue
        info = classify(p)
        info.update(name=p.name, path=str(p), size=st.st_size,
                    mtime=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime)))
        files.append(info)
    parent = d.parent if _inside(d.parent, roots) and d.parent != d else None
    return {"dir": str(d), "parent": str(parent) if parent else None, "files": files, "dirs": dirs,
            "roots": [str(r) for r in roots]}


def _numeric_rows(path: Path, max_rows: Optional[int] = None) -> Tuple[List[List[float]], int]:
    rows: List[List[float]] = []
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s[0] in "#*":
                continue
            try:
                vals = [float(t.replace("D", "E")) for t in s.replace(",", " ").split()]
            except ValueError:
                continue
            n += 1
            if max_rows is None or len(rows) < max_rows:
                rows.append(vals)
    return rows, n


def file_info(path: Path) -> Dict[str, Any]:
    """What the Line Selection dialog needs: data columns (abscissa excluded), pair option, points.

    Spectrum files: columns of the first numeric row minus the frequency column (READSPEC
    ``<numLines>``).  Histories: ``pair`` 0 when the first numeric row holds one value (dt, then
    values: ACS SASSI output format) and 1 for time/value pairs (READTH ``<Pair>``).
    """
    info = classify(path)
    info.update(name=path.name, path=str(path))
    if not path.is_file():
        raise PathError(f"{path} is not a file")
    if path.stat().st_size > 64 * 1024 * 1024:
        info.update(columns=0, points=0, note="file too large to inspect")
        return info
    rows, n = _numeric_rows(path, max_rows=4)
    if not rows:
        info.update(columns=0, points=0, pair=0, note="no numeric data")
        return info
    if info["plot"] == "th" or (info["plot"] is None and len(rows[0]) == 1):
        pair = 0 if len(rows[0]) == 1 else 1
        info.update(plot="th", pair=pair, columns=1, points=(n - 1) if pair == 0 else n)
    else:
        info.update(plot=info["plot"] or "spec", columns=max(len(rows[0]) - 1, 1), points=n)
    return info


def read_text(path: Path) -> str:
    if not path.is_file():
        raise PathError(f"{path} is not a file")
    if path.stat().st_size > MAX_TEXT_BYTES:
        raise PathError(f"{path.name} is larger than {MAX_TEXT_BYTES // (1024 * 1024)} MB (open it with another editor)")
    with open(path, "rb") as fh:
        head = fh.read(4096)
    if b"\x00" in head:
        raise PathError(f"{path.name} is a binary file")
    return path.read_text(encoding="utf-8-sig", errors="replace")


def write_text(path: Path, text: str) -> Path:
    if len(text.encode("utf-8")) > MAX_TEXT_BYTES:
        raise PathError("text too large")
    if not path.parent.is_dir():
        raise PathError(f"directory {path.parent} does not exist")
    path.write_text(text, encoding="utf-8")
    return path

