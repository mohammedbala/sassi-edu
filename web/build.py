"""Build the static browser version of SASSI-EDU (GitHub Pages) into ``web/site/``.

Usage (from the project root; needs only the standard library and the ``plotly`` package, for Plotly.js)::

    python web/build.py                 # writes web/site/ and prints its size
    python web/build.py --out DIR       # another output directory (emptied first)

The site (see ``web/README.md`` and ``docs/internal/web_version.md``)::

    index.html                  sassi/ui/static/index.html without the session token, with web/boot.js
    .nojekyll                   GitHub Pages: serve the files as they are
    static/                     the GUI front end (app.js, dialogs.js, plots.js, learn.js, styles.css, katex/,
                                examples/ -- the gallery pictures); every file of sassi/ui/static, subfolders too
    static/plotly.min.js        Plotly.js of the installed plotly package
    web/boot.js, web/worker.js  the loading overlay and transport; the Web Worker running Pyodide
    web/sassi-edu.zip           what Python needs: the sassi package, README.md, examples/ and the public
                                documentation (docs/index.md, user, theory, verification, reference)
    docs/...                    the figures of the Help pages (PNG / JPEG / GIF)

Never bundled: ``reference/`` (the ACS SASSI manual text), ``docs/spec/``, ``docs/internal/``, ``tests/``,
``examples/sassi-course/`` (the course workspaces of local GUI runs), virtual environments and caches --
:data:`FORBIDDEN` is checked for every file written.  The build is reproducible (sorted files, fixed zip
time stamps); its id (a hash of the content) is put into the URLs of the page so that a browser never mixes
files of two builds.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import re
import shutil
import sys
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

ROOT = Path(__file__).resolve().parents[1]           # the project root (sassi/, docs/, examples/)
WEB = ROOT / "web"
STATIC = ROOT / "sassi" / "ui" / "static"
DEFAULT_OUT = WEB / "site"
#: the public documentation Python reads for Help and the lessons (relative to the project root)
DOCS = ("docs/index.md", "docs/user", "docs/theory", "docs/verification", "docs/reference")
#: figures of the Help pages, also published as files of the site (the page shows them as docs/<path>)
IMAGE_SUFFIXES = (".png", ".jpg", ".jpeg", ".gif")
#: never part of the site (copyright, internal specifications, tests, local runs, environments)
FORBIDDEN = ("reference/", "docs/spec/", "docs/internal/", "tests/", "examples/sassi-course/", ".venv/", "venv/",
             "web/site/", "web/node_modules/")
#: path parts that are never copied (caches, editor and OS files)
SKIP_PARTS = ("__pycache__", ".pytest_cache", ".DS_Store", ".git", ".ipynb_checkpoints")
CSP = ("default-src 'self'; script-src 'self' 'unsafe-eval'; style-src 'self' 'unsafe-inline'; "
       "img-src 'self' data: blob:; connect-src 'self'; worker-src 'self'; object-src 'none'; base-uri 'self'; "
       "form-action 'none'")
ZIP_TIME = (1980, 1, 1, 0, 0, 0)


class BuildError(RuntimeError):
    """The site cannot be built (a missing input, a forbidden file)."""


def _skip(rel: str) -> bool:
    parts = rel.split("/")
    return any(p in SKIP_PARTS or p.endswith(".pyc") for p in parts) or parts[-1].startswith(".")


def _check(rel: str) -> str:
    if any(rel == f.rstrip("/") or rel.startswith(f) for f in FORBIDDEN):
        raise BuildError(f"{rel} must never be published")
    return rel


def _files(rel_dir: str, exclude: Iterable[str] = ()) -> List[str]:
    """The files under ``rel_dir`` (project-relative POSIX paths, sorted), without caches and ``exclude``."""
    base = ROOT / rel_dir
    if base.is_file():
        return [rel_dir]
    if not base.is_dir():
        raise BuildError(f"{rel_dir} not found")
    out = []
    for p in sorted(base.rglob("*")):
        rel = p.relative_to(ROOT).as_posix()
        if p.is_file() and not _skip(rel) and not any(rel.startswith(e) for e in exclude):
            out.append(rel)
    return out


def bundle_files() -> List[str]:
    """What ``web/sassi-edu.zip`` holds (project-relative paths)."""
    rels = _files("sassi", exclude=("sassi/ui/static/",))         # the front end is a part of the site itself
    rels += ["README.md"] + _files("examples", exclude=("examples/sassi-course/",))
    for d in DOCS:
        rels += _files(d)
    return [_check(r) for r in rels]


def make_zip(rels: List[str]) -> bytes:
    """A reproducible zip of the project files ``rels``."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for rel in rels:
            info = zipfile.ZipInfo(rel, date_time=ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, (ROOT / rel).read_bytes())
    return buf.getvalue()


def plotly_js() -> Path:
    """``plotly.min.js`` of the installed plotly package (found without importing it)."""
    spec = importlib.util.find_spec("plotly")
    if spec is None or not spec.origin:
        raise BuildError("the plotly package is not installed (pip install plotly): Plotly.js comes from it")
    p = Path(spec.origin).parent / "package_data" / "plotly.min.js"
    if not p.is_file():
        raise BuildError(f"{p} not found")
    return p


def page(build_id: str) -> str:
    """``index.html`` of the site: the GUI page without the session token, with web/boot.js first and the build
    id on every local URL."""
    html = (STATIC / "index.html").read_text(encoding="utf-8")
    html, n = re.subn(r'\s*<meta name="sassi-token" content="\{\{TOKEN\}\}">', "", html)
    if n != 1:
        raise BuildError("index.html: the session-token meta tag was not found")
    head = (f'\n<meta http-equiv="Content-Security-Policy" content="{CSP}">'
            '\n<meta name="description" content="SASSI-EDU: learn soil-structure interaction (the SASSI flexible-volume '
            'method) with a guided course, examples and plots; Python runs in your browser.">')
    html = html.replace('<meta name="viewport" content="width=device-width, initial-scale=1">',
                        '<meta name="viewport" content="width=device-width, initial-scale=1">' + head, 1)
    html = re.sub(r'((?:src|href)="static/[^"?]+)"', lambda m: f'{m.group(1)}?v={build_id}"', html)
    first = html.index('<script src="static/')
    html = html[:first] + f'<script src="web/boot.js?v={build_id}"></script>\n' + html[first:]
    if "{{" in html or 'src="/' in html or 'href="/' in html:
        raise BuildError("index.html: a placeholder or a root-absolute URL is left")
    return html


def build(out: Path = DEFAULT_OUT) -> Dict[str, object]:
    """Write the site into ``out`` (emptied first); returns a summary (files, sizes, build id)."""
    out = Path(out).resolve()
    # only an empty folder or an earlier build is emptied (never the project or a folder of other files)
    if out == ROOT or out in ROOT.parents or (out.is_dir() and any(out.iterdir())
                                              and not (out / "web" / "sassi-edu.zip").is_file()):
        raise BuildError(f"refusing to empty {out}: not an earlier build of the site")
    rels = bundle_files()
    zbytes = make_zip(rels)
    statics: List[Tuple[str, Path]] = [("static/" + p.relative_to(STATIC).as_posix(), p)
                                       for p in sorted(STATIC.rglob("*"))
                                       if p.is_file() and p.name != "index.html"
                                       and not _skip(p.relative_to(STATIC).as_posix())]
    statics.append(("static/plotly.min.js", plotly_js()))
    images = [(r, ROOT / r) for r in rels if r.startswith("docs/") and r.lower().endswith(IMAGE_SUFFIXES)]
    h = hashlib.sha256(zbytes)
    for name, p in statics + [("web/boot.js", WEB / "boot.js"), ("web/worker.js", WEB / "worker.js")]:
        h.update(name.encode() + b"\0" + p.read_bytes())
    build_id = h.hexdigest()[:12]

    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    written: Dict[str, int] = {}

    def put(name: str, data: bytes) -> None:
        _check(name)
        dest = out / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
        written[name] = len(data)

    put("index.html", page(build_id).encode("utf-8"))
    put(".nojekyll", b"")
    for name, p in statics + images:
        put(name, p.read_bytes())
    for name in ("boot.js", "worker.js"):
        text = (WEB / name).read_text(encoding="utf-8").replace("__SASSI_BUILD__", build_id)
        put("web/" + name, text.encode("utf-8"))
    put("web/sassi-edu.zip", zbytes)
    total = sum(written.values())
    return {"out": str(out), "build": build_id, "files": sorted(written), "bytes": total, "zip_bytes": len(zbytes),
            "zip_files": len(rels), "sizes": written}


def _mb(n: int) -> str:
    return f"{n / 1e6:.2f} MB"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="build the static browser version of SASSI-EDU (GitHub Pages)")
    ap.add_argument("--out", default=str(DEFAULT_OUT), help="output directory, emptied first (default web/site)")
    args = ap.parse_args(argv)
    try:
        s = build(Path(args.out))
    except BuildError as exc:
        print(f"build failed: {exc}", file=sys.stderr)
        return 1
    sizes = s["sizes"]
    print(f"SASSI-EDU web version {s['build']} written to {s['out']}")
    print(f"  {len(s['files'])} files, {_mb(s['bytes'])} in all")
    print(f"  web/sassi-edu.zip    {_mb(s['zip_bytes'])} ({s['zip_files']} files: sassi, examples, docs)")
    print(f"  static/plotly.min.js {_mb(sizes['static/plotly.min.js'])}")
    katex = sum(v for k, v in sizes.items() if k.startswith("static/katex/"))
    print(f"  static/katex/        {_mb(katex)}")
    print("  + Pyodide 314.0.7, NumPy and SciPy from cdn.jsdelivr.net at the first visit (cached by the browser)")
    print(f"test locally: python -m http.server 8811 -d {s['out']}   then open http://127.0.0.1:8811/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
