"""Help > Help: the project documentation rendered to HTML (requirements 5.3 Help (F1), UI-08).

The Help tab lists and shows the Markdown documents of the project *as they appear on disk* (no
copy is made, so a regenerated Verification Manual is shown at once):

==============  ==============================================================================
group           documents
==============  ==============================================================================
Start here      ``docs/index.md`` (the documentation map)
User guides     ``docs/user/*.md`` (User Guide, GUI guide, ANSYS interface, and any later guide)
Theory          ``docs/theory/*.md`` (Theory Manual)
Verification    ``docs/verification/*.md`` (Verification Manual first)
Reference       ``docs/reference/*.md`` (Command Reference)
Examples        ``examples/README.md`` (the tutorial examples)
==============  ==============================================================================

Links between documents (``[Theory Manual §3](../theory/THEORY_MANUAL.md#3-flexible-volume-...)``)
are rendered as in-GUI links (``data-doc`` / ``data-anchor``): the browser loads the target
document and scrolls to the heading with that GitHub slug.  Every Markdown file under ``docs/``
(except ``docs/internal/``) and ``examples/README.md``, ``sassi/prep/README.md``, ``README.md`` can be
opened that way (the specification cross-links of the documentation index).  Images under
``docs/`` (PNG, JPEG, GIF; e.g. ``docs/verification/figures/*.png``) are served at ``/docs/<path>``
(:func:`image_file`); the page refers to them with the *relative* URL ``docs/<path>``, which also
works for the static browser version under ``https://<user>.github.io/<repo>/``.  Anything else -- a directory, a source file, a path leaving the project --
is shown as plain text.  Rendering: :mod:`sassi.ui.markdown` (no raw HTML is passed through).
"""
from __future__ import annotations

import posixpath
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import quote, unquote

from .markdown import default_link, render

PROJECT_ROOT = Path(__file__).resolve().parents[2]
HOME = "docs/index.md"
#: the groups of the Help document list: (title, glob patterns relative to the project root; the
#: first pattern of a group is listed first)
GROUPS: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("Start here", ("docs/index.md",)),
    ("User guides", ("docs/user/USER_GUIDE.md", "docs/user/GUI.md", "docs/user/*.md")),
    ("Theory", ("docs/theory/THEORY_MANUAL.md", "docs/theory/*.md")),
    ("Verification", ("docs/verification/VERIFICATION_MANUAL.md", "docs/verification/*.md")),
    ("Reference", ("docs/reference/COMMAND_REFERENCE.md", "docs/reference/*.md")),
    ("Examples", ("examples/README.md",)),
)
#: Markdown files outside docs/ that the documentation links to
EXTRA_READABLE = ("examples/README.md", "sassi/prep/README.md", "README.md")
#: image types served at /docs/<path> (no SVG: an SVG opened directly could run script)
IMAGE_TYPES = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg", ".gif": "image/gif"}
MAX_DOC_BYTES = 8 * 1024 * 1024
_SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.\-]*:")
_H1 = re.compile(r"^#\s+(.+?)\s*#*\s*$")
_SECTION = re.compile(r'<h([1-6]) data-slug="([^"]*)">(.*?)</h\1>', re.S)
_TAGS = re.compile(r"<[^>]+>")


class HelpError(ValueError):
    """A document that cannot be shown (``status``: 404 not found / not a help document)."""

    def __init__(self, message: str, status: int = 404):
        super().__init__(message)
        self.status = status


def _inside(p: Path, root: Path) -> bool:
    try:
        p.relative_to(root)
        return True
    except ValueError:
        return False


class HelpDocs:
    """The documentation of one project tree (``root``: the directory holding ``docs/``)."""

    def __init__(self, root: Optional[Path] = None):
        self.root = Path(root or PROJECT_ROOT).resolve()
        self._cache: Dict[str, Tuple[Tuple[int, int], Dict[str, Any]]] = {}

    # ------------------------------------------------------------------ catalogue
    @property
    def available(self) -> bool:
        return (self.root / "docs").is_dir()

    def normalize(self, name: str) -> Optional[str]:
        """``name`` as a normalised POSIX path relative to the root, or None when it leaves the root."""
        s = unquote(str(name or "")).strip().replace("\\", "/")
        if not s or "\x00" in s or s.startswith("/") or _SCHEME.match(s):
            return None
        s = posixpath.normpath(s)
        if s == "." or s == ".." or s.startswith("../"):
            return None
        return s

    def is_readable(self, rel: Optional[str]) -> bool:
        """A Markdown document of the help set (docs/ except docs/internal/, and EXTRA_READABLE)."""
        if not rel or not rel.lower().endswith(".md"):
            return False
        if not ((rel.startswith("docs/") and not rel.startswith("docs/internal/")) or rel in EXTRA_READABLE):
            return False
        p = (self.root / rel).resolve()
        return p.is_file() and _inside(p, self.root)

    def title_of(self, rel: str) -> str:
        """The first level-1 heading of a document (its file name when it has none)."""
        try:
            with open(self.root / rel, "r", encoding="utf-8", errors="replace") as fh:
                in_fence = False
                for k, ln in enumerate(fh):
                    if k > 80:
                        break
                    if ln.lstrip().startswith(("```", "~~~")):
                        in_fence = not in_fence
                    m = None if in_fence else _H1.match(ln)
                    if m:
                        return re.sub(r"[`*]", "", m.group(1)).strip()
        except OSError:
            pass
        return posixpath.basename(rel)

    def catalogue(self) -> List[Dict[str, Any]]:
        """The documents listed by Help, in group order: ``{"id", "title", "group", "size"}``."""
        out: List[Dict[str, Any]] = []
        seen = set()
        for group, patterns in GROUPS:
            for pat in patterns:
                if "*" in pat:
                    rels = sorted(p.relative_to(self.root).as_posix() for p in self.root.glob(pat) if p.is_file())
                else:
                    rels = [pat]
                for rel in rels:
                    if rel in seen or not self.is_readable(rel):
                        continue
                    seen.add(rel)
                    out.append({"id": rel, "title": self.title_of(rel), "group": group,
                                "size": (self.root / rel).stat().st_size})
        return out

    # ------------------------------------------------------------------ links and images
    def _resolver(self, doc: str):
        base = posixpath.dirname(doc)

        def resolve(path: str) -> Optional[str]:
            p = posixpath.normpath(posixpath.join(base, unquote(path)))
            return None if p == ".." or p.startswith("../") else p

        def link(href: str) -> Optional[Dict[str, str]]:
            h = href.strip()
            if not h:
                return None
            if h.startswith("#") or _SCHEME.match(h):
                return default_link(h)
            path, _, anchor = h.partition("#")
            target = resolve(path) if path else doc
            if target is None:
                return None
            if self.is_readable(target):
                d = {"href": "#" + anchor, "data-doc": target}
                if anchor:
                    d["data-anchor"] = anchor
                return d
            if image_file(self.root, target) is not None:
                return {"href": "docs/" + quote(target[len("docs/"):]), "target": "_blank",
                        "rel": "noopener noreferrer"}
            return None

        def image(src: str) -> Optional[str]:
            s = src.strip()
            if not s or _SCHEME.match(s) or s.startswith(("/", "#")):
                return None                 # no external images: the GUI works offline (CSP img-src 'self')
            target = resolve(s.partition("#")[0])
            if target is None or image_file(self.root, target) is None:
                return None
            return "docs/" + quote(target[len("docs/"):])

        return link, image

    # ------------------------------------------------------------------ documents
    def document(self, name: str) -> Dict[str, Any]:
        """The rendered document ``name`` (project-relative path): ``{"id", "title", "group", "html",
        "toc", "size"}``.  Raises :class:`HelpError` for anything that is not a help document."""
        rel = self.normalize(name)
        if rel is None or not self.is_readable(rel):
            raise HelpError(f"{name}: not a help document (Markdown files under docs/ and examples/README.md)")
        p = self.root / rel
        st = p.stat()
        if st.st_size > MAX_DOC_BYTES:
            raise HelpError(f"{rel}: larger than {MAX_DOC_BYTES // (1024 * 1024)} MB", 413)
        key = (st.st_mtime_ns, st.st_size)
        hit = self._cache.get(rel)
        if hit is not None and hit[0] == key:
            return hit[1]
        link, image = self._resolver(rel)
        r = render(p.read_text(encoding="utf-8", errors="replace"), link=link, image=image)
        group = next((d["group"] for d in self.catalogue() if d["id"] == rel), "")
        doc = {"id": rel, "title": r.title or posixpath.basename(rel), "group": group, "html": r.html,
               "toc": r.toc, "size": st.st_size}
        if len(self._cache) > 24:
            self._cache.pop(next(iter(self._cache)))
        self._cache[rel] = (key, doc)
        return doc

    def search(self, query: str, limit: int = 200) -> Dict[str, Any]:
        """Sections of the listed documents whose text contains ``query`` (case-insensitive):
        ``{"query", "results": [{"doc", "title", "slug", "heading", "snippet"}], "truncated"}``."""
        q = (query or "").strip()
        if len(q) < 2:
            return {"query": q, "results": [], "truncated": False}
        ql = q.lower()
        results: List[Dict[str, Any]] = []
        for d in self.catalogue():
            try:
                doc = self.document(d["id"])
            except (HelpError, OSError):
                continue
            html = doc["html"]
            heads = list(_SECTION.finditer(html))
            spans = [(None, "", 0, heads[0].start() if heads else len(html))] if not heads or heads[0].start() > 0 else []
            for k, m in enumerate(heads):
                end = heads[k + 1].start() if k + 1 < len(heads) else len(html)
                spans.append((m.group(2), _plain(m.group(3)), m.start(), end))
            for slug, heading, a, b in spans:
                text = _plain(html[a:b])
                pos = text.lower().find(ql)
                if pos < 0:
                    continue
                lo = max(0, pos - 70)
                snippet = ("..." if lo else "") + text[lo:pos + len(q) + 90].strip() + "..."
                results.append({"doc": d["id"], "title": d["title"], "slug": slug or "", "heading": heading,
                                "snippet": snippet})
                if len(results) >= limit:
                    return {"query": q, "results": results, "truncated": True}
        return {"query": q, "results": results, "truncated": False}


def _plain(html: str) -> str:
    import html as _html
    text = re.sub(r"\s+", " ", _html.unescape(_TAGS.sub(" ", html))).strip()
    return re.sub(r"\s+([,.;:)])(?=\s|$)", r"\1", text).replace("( ", "(")


def image_file(root: Path, rel: str) -> Optional[Path]:
    """The image ``rel`` (project-relative, under ``docs/``) when it may be served, else None."""
    if not rel or not rel.startswith("docs/") or rel.startswith("docs/internal/"):
        return None
    root = Path(root).resolve()
    p = (root / rel).resolve()
    if p.suffix.lower() not in IMAGE_TYPES or not p.is_file() or not _inside(p, root / "docs"):
        return None
    return p
