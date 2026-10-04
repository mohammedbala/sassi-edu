"""A small, safe Markdown -> HTML renderer for the Help pages (requirements 5.3 Help, UI-08).

Help > Help shows the documentation of the project (``docs/``, ``examples/README.md``) inside the
GUI.  The documents are GitHub-flavoured Markdown; this module renders the subset they use, with no
dependency outside the standard library:

* blocks: ATX headings (``#`` ... ``######``), paragraphs (hard breaks: two trailing spaces or a
  trailing backslash), fenced code (``````` and ``~~~``), block quotes, bullet and ordered lists
  (nested by indentation, lazy continuation lines, tight and loose lists), GFM tables (alignment,
  ``\\|`` escapes), thematic breaks (``---``);
* inline: code spans, backslash escapes, ``**strong**`` / ``__strong__``, ``*em*`` / ``_em_``
  (never inside a word, so ``K*s`` and ``FILE_8`` stay literal), ``~~del~~``, links
  ``[text](href "title")``, images ``![alt](src)``, autolinks ``<https://...>`` and bare
  ``http(s)://`` URLs, HTML entities (``&lt;``);
* **LaTeX math**: ``$...$`` inline, ``$$...$$`` (on its own lines) and fenced ``math`` blocks as display
  equations.  The TeX source is kept verbatim (no emphasis or escape processing inside) in a
  ``data-tex`` attribute of ``<span class="math">`` / ``<div class="math-display">``, with the source as
  the element's text; the browser typesets it with KaTeX (``static/katex``).  ``\\$`` is a literal
  dollar; a ``$`` inside a code span is literal.  :func:`math_spans` lists the formulas of a text (used
  to check them, ``tests/unit/test_lessons.py``).

**Safety.**  Raw HTML is *never* passed through: every ``<``, ``>`` and ``&`` of the text is escaped
(only well-formed character references such as ``&lt;`` are kept), attribute values are quoted and
escaped, and link / image targets go through resolver callbacks -- a target the resolver does not
accept is rendered as plain text, so ``javascript:`` and other schemes cannot appear.  The page's
Content-Security-Policy additionally blocks inline scripts.

**Anchors.**  Headings get the GitHub slug of their text (lower case, punctuation removed, spaces to
``-``, ``-1``, ``-2`` ... for repeats) in a ``data-slug`` attribute -- not as an ``id``, so a heading
can never collide with an element id of the GUI page.  ``[text](#slug)`` and links to other
documents (``../user/GUI.md#3-menus``) are resolved by the browser with these slugs.
"""
from __future__ import annotations

import html as _html
import re
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

#: link resolver: Markdown href -> attributes of the <a> element (None: render the text only)
LinkResolver = Callable[[str], Optional[Dict[str, str]]]
#: image resolver: Markdown src -> URL of the <img> (None: show the alt text)
ImageResolver = Callable[[str], Optional[str]]

#: attributes a resolver may set on a link
LINK_ATTRS = ("href", "title", "target", "rel", "class", "data-doc", "data-anchor")
#: href prefixes accepted on a rendered link (after the resolver); ``docs/`` is the relative URL of the
#: documentation figures (the page also works under a sub-path, e.g. GitHub Pages)
SAFE_HREF = re.compile(r"^(#|/?docs/|https?://|mailto:)", re.I)

_ENTITY = re.compile(r"&(?:#[0-9]{1,7}|#[xX][0-9a-fA-F]{1,6}|[A-Za-z][A-Za-z0-9]{1,31});")
_FENCE = re.compile(r"^( {0,3})(`{3,}|~{3,})[ \t]*([\w+#.\-]*)[^`]*$")
_ATX = re.compile(r"^ {0,3}(#{1,6})(?=[ \t]|$)(.*)$")
_HR = re.compile(r"^ {0,3}([-*_])(?:[ \t]*\1){2,}[ \t]*$")
_QUOTE = re.compile(r"^ {0,3}> ?(.*)$")
_ITEM = re.compile(r"^( {0,3})([-*+]|\d{1,9}[.)])([ \t]+.*|[ \t]*)$")
_DELIM_CELL = re.compile(r"^:?-+:?$")

_CODE_SPAN = re.compile(r"(`+)(.+?)(?<!`)\1(?!`)", re.S)
_ESCAPABLE = re.compile(r"\\([!\"#$%&'()*+,\-./:;<=>?@\[\\\]^_`{|}~])")
_HARD_BREAK = re.compile(r"(?: {2,}|\\)\n")
_AUTOLINK = re.compile(r"<((?:https?|mailto):[^<>\s]+)>", re.I)
_TARGET = r"(?:<([^<>\n]*)>|((?:[^\s()<>]|\([^\s()<>]*\))+))(?:\s+(?:\"([^\"]*)\"|'([^']*)'))?"
_IMAGE = re.compile(r"!\[((?:[^\[\]]|\[[^\[\]]*\])*)\]\(\s*" + _TARGET + r"\s*\)")
_LINK = re.compile(r"\[((?:[^\[\]]|\[[^\[\]]*\])*)\]\(\s*" + _TARGET + r"\s*\)")
_BARE_URL = re.compile(r"(?<![\w/\"'=@])(https?://[^\s<>\"'\x00]+)", re.I)
_STRONG_STAR = re.compile(r"\*\*(?=[^\s*])(.+?)(?<=[^\s*])\*\*", re.S)
_STRONG_UNDER = re.compile(r"(?<![\w_])__(?=[^\s_])(.+?)(?<=[^\s_])__(?![\w_])", re.S)
_EM_STAR = re.compile(r"(?<![\w*])\*(?=[^\s*])(.+?)(?<=[^\s*])\*(?![\w*])", re.S)
_EM_UNDER = re.compile(r"(?<![\w_])_(?=[^\s_])(.+?)(?<=[^\s_])_(?![\w_])", re.S)
_DEL = re.compile(r"~~(?=\S)(.+?)(?<=\S)~~", re.S)
_PLACEHOLDER = re.compile("\x00(\\d+)\x00")
#: inline math: ``$`` not preceded by a word character, ``$`` or ``\``, not followed by a space; the
#: closing ``$`` not preceded by a space or ``\`` and not followed by a word character or a digit
_MATH_INLINE = re.compile(r"(?<![\w$\\])\$(?=[^\s$])((?:\\.|[^$\\])+?)(?<![\s\\])\$(?![\w$])", re.S)
_MATH_FENCE_LANGS = ("math", "latex", "tex")
_TAG = re.compile(r"<[^>]*>")


@dataclass
class Rendered:
    """Result of :func:`render`: the HTML, the title (first level-1 heading) and the headings."""
    html: str
    title: str = ""
    toc: List[Dict[str, Any]] = field(default_factory=list)   # {"level", "text", "slug"}


def escape_text(s: str) -> str:
    """HTML-escape text (``&``, ``<``, ``>``) keeping well-formed character references (``&lt;``)."""
    out: List[str] = []
    pos = 0
    for m in _ENTITY.finditer(s):
        out.append(_html.escape(s[pos:m.start()], quote=False))
        out.append(m.group(0))
        pos = m.end()
    out.append(_html.escape(s[pos:], quote=False))
    return "".join(out)


def attr(s: str) -> str:
    """An attribute value (quotes escaped as well)."""
    return _html.escape(str(s), quote=True)


def slugify(text: str) -> str:
    """GitHub heading slug: lower case, characters other than letters, digits, ``_``, ``-`` and
    spaces removed, spaces replaced by ``-`` (``1.3 A first analysis`` -> ``13-a-first-analysis``)."""
    s = re.sub(r"[^\w\- ]", "", text.strip().lower())
    return s.replace(" ", "-")


def default_link(href: str) -> Optional[Dict[str, str]]:
    """Links accepted without a resolver: anchors of the same page and external web / mail links."""
    h = href.strip()
    if h.startswith("#") and len(h) > 1:
        return {"href": h, "data-anchor": h[1:]}
    if re.match(r"^(https?://|mailto:)", h, re.I):
        return {"href": h, "target": "_blank", "rel": "noopener noreferrer"}
    return None


def _indent(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def _blank(line: str) -> bool:
    return not line.strip()


class _Renderer:
    def __init__(self, link: Optional[LinkResolver], image: Optional[ImageResolver]):
        self.link = link or default_link
        self.image = image
        self.toc: List[Dict[str, Any]] = []
        self.title = ""
        self._slugs: Dict[str, int] = {}

    # ================================================================== blocks
    def blocks(self, lines: List[str], tight: bool = False) -> str:
        out: List[str] = []
        i, n = 0, len(lines)
        while i < n:
            line = lines[i]
            if _blank(line):
                i += 1
                continue
            m = _FENCE.match(line)
            if m:
                i = self._fence(lines, i, m, out)
                continue
            if line.lstrip().startswith("$$"):
                i = self._display_math(lines, i, out)
                continue
            m = _ATX.match(line)
            if m:
                out.append(self._heading(m))
                i += 1
                continue
            if _HR.match(line):
                out.append("<hr>")
                i += 1
                continue
            if _QUOTE.match(line):
                i = self._quote(lines, i, out)
                continue
            if self._table_start(lines, i):
                i = self._table(lines, i, out)
                continue
            if _ITEM.match(line):
                i = self._list(lines, i, out)
                continue
            i = self._paragraph(lines, i, out, tight)
        return "\n".join(out)

    def _interrupts(self, line: str) -> bool:
        """True when ``line`` ends a paragraph (CommonMark: a bullet item with content, or an
        ordered item numbered 1, can interrupt a paragraph; other numbers cannot -- a wrapped line
        starting with ``1994. The`` stays text)."""
        if (_blank(line) or _FENCE.match(line) or _ATX.match(line) or _HR.match(line) or _QUOTE.match(line)
                or line.lstrip().startswith("$$")):
            return True
        m = _ITEM.match(line)
        if m and m.group(3).strip():
            mk = m.group(2)
            return mk[:-1] == "1" if mk[0].isdigit() else True
        return False

    def _paragraph(self, lines: List[str], i: int, out: List[str], tight: bool) -> int:
        buf = [lines[i].lstrip()]
        i += 1
        while i < len(lines) and not self._interrupts(lines[i]) and not self._table_start(lines, i):
            buf.append(lines[i].lstrip())
            i += 1
        inner = self.inline("\n".join(buf).rstrip())
        out.append(inner if tight else f"<p>{inner}</p>")
        return i

    def _fence(self, lines: List[str], i: int, m, out: List[str]) -> int:
        indent, fence, info = len(m.group(1)), m.group(2), m.group(3)
        close = re.compile(r"^ {0,3}" + re.escape(fence[0]) + "{" + str(len(fence)) + r",}[ \t]*$")
        body: List[str] = []
        i += 1
        while i < len(lines) and not close.match(lines[i]):
            ln = lines[i]
            k = 0
            while k < indent and k < len(ln) and ln[k] == " ":
                k += 1
            body.append(ln[k:])
            i += 1
        if info.lower() in _MATH_FENCE_LANGS:
            out.append(math_html("\n".join(body).strip(), display=True))
            return i + 1
        cls = f' class="language-{attr(info)}"' if info else ""
        out.append(f"<pre><code{cls}>{_html.escape(chr(10).join(body), quote=False)}</code></pre>")
        return i + 1

    def _display_math(self, lines: List[str], i: int, out: List[str]) -> int:
        """``$$ ... $$``: from a line starting with ``$$`` to the line ending with ``$$`` (the same line
        for a one-line equation).  An unclosed ``$$`` runs to the next blank line."""
        first = lines[i].strip()[2:]
        body: List[str] = []
        if first.rstrip().endswith("$$"):
            out.append(math_html(first.rstrip()[:-2].strip(), display=True))
            return i + 1
        body.append(first)
        i += 1
        while i < len(lines) and not _blank(lines[i]):
            ln = lines[i].rstrip()
            if ln.endswith("$$"):
                body.append(ln[:-2])
                i += 1
                break
            body.append(ln)
            i += 1
        out.append(math_html("\n".join(body).strip(), display=True))
        return i

    def _heading(self, m) -> str:
        level = len(m.group(1))
        text = m.group(2).strip()
        text = re.sub(r"(?:^|[ \t]+)#+[ \t]*$", "", text)          # optional closing sequence
        inner = self.inline(text)
        plain = _html.unescape(_TAG.sub("", inner)).strip()
        slug = self._unique(slugify(plain))
        self.toc.append({"level": level, "text": plain, "slug": slug})
        if level == 1 and not self.title:
            self.title = plain
        return f'<h{level} data-slug="{attr(slug)}">{inner}</h{level}>'

    def _unique(self, slug: str) -> str:
        k = self._slugs.get(slug, 0)
        self._slugs[slug] = k + 1
        return slug if k == 0 else f"{slug}-{k}"

    def _quote(self, lines: List[str], i: int, out: List[str]) -> int:
        buf: List[str] = []
        while i < len(lines):
            m = _QUOTE.match(lines[i])
            if m:
                buf.append(m.group(1))
                i += 1
                continue
            if buf and not _blank(buf[-1]) and not self._interrupts(lines[i]):     # lazy continuation
                buf.append(lines[i])
                i += 1
                continue
            break
        out.append("<blockquote>\n" + self.blocks(buf) + "\n</blockquote>")
        return i

    # ------------------------------------------------------------------ tables (GFM)
    @staticmethod
    def _cells(line: str) -> List[str]:
        s = line.strip().replace("\\|", "\x02")
        if s.startswith("|"):
            s = s[1:]
        if s.endswith("|"):
            s = s[:-1]
        return [c.strip().replace("\x02", "|") for c in s.split("|")]

    def _table_start(self, lines: List[str], i: int) -> bool:
        if i + 1 >= len(lines) or "|" not in lines[i] or "-" not in lines[i + 1]:
            return False
        delim = self._cells(lines[i + 1])
        if not delim or not all(_DELIM_CELL.match(c.replace(" ", "")) for c in delim):
            return False
        if "|" not in lines[i + 1] and len(delim) == 1:
            return False                                   # a bare '---' under text is not a table
        return len(self._cells(lines[i])) == len(delim)

    def _table(self, lines: List[str], i: int, out: List[str]) -> int:
        head = self._cells(lines[i])
        aligns = []
        for c in self._cells(lines[i + 1]):
            c = c.replace(" ", "")
            aligns.append("center" if c.startswith(":") and c.endswith(":") else
                          "right" if c.endswith(":") else "left" if c.startswith(":") else "")
        i += 2
        rows: List[List[str]] = []
        while i < len(lines) and not _blank(lines[i]) and "|" in lines[i] and not (
                _FENCE.match(lines[i]) or _ATX.match(lines[i]) or _QUOTE.match(lines[i])):
            rows.append(self._cells(lines[i]))
            i += 1
        ncol = len(head)

        def cell(tag: str, text: str, k: int) -> str:
            cls = f' class="al-{aligns[k]}"' if aligns[k] else ""
            return f"<{tag}{cls}>{self.inline(text)}</{tag}>"

        parts = ['<div class="md-table"><table>', "<thead><tr>" + "".join(cell("th", h, k) for k, h in enumerate(head))
                 + "</tr></thead>", "<tbody>"]
        for r in rows:
            r = (r + [""] * ncol)[:ncol]
            parts.append("<tr>" + "".join(cell("td", c, k) for k, c in enumerate(r)) + "</tr>")
        parts.append("</tbody></table></div>")
        out.append("".join(parts))
        return i

    # ------------------------------------------------------------------ lists
    def _list(self, lines: List[str], i: int, out: List[str]) -> int:
        m = _ITEM.match(lines[i])
        mk0 = m.group(2)
        ordered = mk0[0].isdigit()
        delim = mk0[-1]
        start = int(mk0[:-1]) if ordered else 1
        base = len(m.group(1))
        n = len(lines)
        items: List[List[str]] = []
        loose = False

        def same_list(mm) -> bool:
            mk = mm.group(2)
            return mk[0].isdigit() == ordered and mk[-1] == delim and len(mm.group(1)) <= base + 3

        while i < n:
            m = _ITEM.match(lines[i])
            if not m or _HR.match(lines[i]) or not same_list(m):
                break
            ind, mk, rest = len(m.group(1)), m.group(2), m.group(3)
            content = rest.lstrip(" \t")
            sp = len(rest) - len(content)
            if not content:
                width = ind + len(mk) + 1
            elif sp > 4:                                   # indented code inside the item
                width = ind + len(mk) + 1
                content = rest[1:]
            else:
                width = ind + len(mk) + sp
            item = [content]
            i += 1
            while i < n:
                ln = lines[i]
                if _blank(ln):
                    j = i
                    while j < n and _blank(lines[j]):
                        j += 1
                    if j < n and _indent(lines[j]) >= width:
                        item.extend([""] * (j - i))
                        loose = True
                        i = j
                        continue
                    break
                if _indent(ln) >= width:
                    item.append(ln[width:])
                    i += 1
                    continue
                if (item[-1].strip() and not self._interrupts(ln) and not _ITEM.match(ln)
                        and not self._table_start(lines, i)):
                    item.append(ln.strip())                # lazy paragraph continuation
                    i += 1
                    continue
                break
            items.append(item)
            j = i
            while j < n and _blank(lines[j]):
                j += 1
            if j > i and j < n:
                m2 = _ITEM.match(lines[j])
                if m2 and not _HR.match(lines[j]) and same_list(m2):
                    loose = True
                    i = j
                    continue
                break
        tag = "ol" if ordered else "ul"
        attrs = f' start="{start}"' if ordered and start != 1 else ""
        body = "".join(f"<li>{self.blocks(it, tight=not loose)}</li>" for it in items)
        out.append(f"<{tag}{attrs}>{body}</{tag}>")
        return i

    # ================================================================== inline
    def inline(self, text: str) -> str:
        stash: List[str] = []

        def put(h: str) -> str:
            stash.append(h)
            return f"\x00{len(stash) - 1}\x00"

        text = text.replace("\x00", "")

        def code(m) -> str:
            c = m.group(2).replace("\n", " ")
            if len(c) > 2 and c[0] == " " and c[-1] == " " and c.strip():
                c = c[1:-1]
            return put("<code>" + _html.escape(c, quote=False) + "</code>")

        text = _CODE_SPAN.sub(code, text)
        text = text.replace("\\$", put("$"))                    # \$ is a literal dollar
        text = _MATH_INLINE.sub(lambda m: put(math_html(m.group(1).strip(), display=False)), text)
        text = _HARD_BREAK.sub(lambda m: put("<br>") + "\n", text)
        text = _ESCAPABLE.sub(lambda m: put(_html.escape(m.group(1), quote=False)), text)
        text = _AUTOLINK.sub(lambda m: put(self._anchor(m.group(1), escape_text(m.group(1)))), text)
        text = _IMAGE.sub(lambda m: put(self._image(m, stash)), text)
        text = _LINK.sub(lambda m: put(self._link(m)), text)
        text = _BARE_URL.sub(lambda m: self._bare_url(m.group(1), put), text)
        return self._restore(self._spans(text), stash)

    @staticmethod
    def _target(m) -> str:
        return (m.group(2) if m.group(2) is not None else m.group(3) or "").strip()

    @staticmethod
    def _title(m) -> Optional[str]:
        return m.group(4) if m.group(4) is not None else m.group(5)

    @staticmethod
    def _spans(text: str) -> str:
        """Escape the text and apply the emphasis rules (placeholders are left untouched)."""
        t = escape_text(text)
        t = _STRONG_STAR.sub(r"<strong>\1</strong>", t)
        t = _STRONG_UNDER.sub(r"<strong>\1</strong>", t)
        t = _EM_STAR.sub(r"<em>\1</em>", t)
        t = _EM_UNDER.sub(r"<em>\1</em>", t)
        return _DEL.sub(r"<del>\1</del>", t)

    @staticmethod
    def _restore(text: str, stash: List[str]) -> str:
        for _ in range(8):                                 # stashed HTML may hold placeholders itself
            if "\x00" not in text:
                break
            text = _PLACEHOLDER.sub(lambda m: stash[int(m.group(1))], text)
        return text.replace("\x00", "")

    def _link(self, m) -> str:
        """``[text](href)``; a target that is not a link at all (no ``/``, ``.``, ``#`` or ``:``, e.g. the
        ``d[a+a0](T)`` of a formula) keeps its source text instead of losing ``(T)``."""
        href = self._target(m)
        if self.link(href) is None and not re.search(r"[/.#:]", href):
            return self._spans(m.group(0))
        return self._anchor(href, self._spans(m.group(1)), self._title(m))

    def _bare_url(self, url: str, put) -> str:
        tail = ""
        while url and url[-1] in ".,:;!?)":
            if url[-1] == ")" and url.count("(") >= url.count(")"):
                break
            tail = url[-1] + tail
            url = url[:-1]
        return put(self._anchor(url, escape_text(url))) + tail

    def _anchor(self, href: str, inner_html: str, title: Optional[str] = None) -> str:
        target = self.link(href)
        if not target or not SAFE_HREF.match(str(target.get("href", ""))):
            return f'<span class="md-nolink" title="{attr(href)}">{inner_html}</span>'
        attrs = dict(target)
        if title and "title" not in attrs:
            attrs["title"] = title
        return "<a" + "".join(f' {k}="{attr(v)}"' for k, v in attrs.items() if k in LINK_ATTRS) + f">{inner_html}</a>"

    def _image(self, m, stash: List[str]) -> str:
        alt = _html.unescape(_TAG.sub("", self._restore(escape_text(m.group(1)), stash)))
        src = self._target(m)
        url = self.image(src) if self.image is not None else None
        if not url or not SAFE_HREF.match(url) or url.startswith(("#", "mailto:")):
            return f'<span class="md-noimage" title="{attr(src)}">[image: {escape_text(alt)}]</span>'
        title = self._title(m)
        return (f'<img src="{attr(url)}" alt="{attr(alt)}"' + (f' title="{attr(title)}"' if title else "")
                + ' loading="lazy">')


def math_html(tex: str, display: bool) -> str:
    """The element of a formula: its TeX in ``data-tex`` (typeset by KaTeX in the browser) and as text."""
    if display:
        return f'<div class="math-display" data-tex="{attr(tex)}">{_html.escape(tex, quote=False)}</div>'
    return f'<span class="math" data-tex="{attr(tex)}">{_html.escape(tex, quote=False)}</span>'


def math_spans(text: str) -> List[Dict[str, Any]]:
    """The formulas of Markdown ``text`` as the renderer finds them: ``{"tex", "display", "line"}`` in
    document order (``line``: 1-based line of the formula's start)."""
    out: List[Dict[str, Any]] = []
    for m in re.finditer(r'<(span|div) class="math(-display)?" data-tex="([^"]*)"', render(text).html):
        out.append({"tex": _html.unescape(m.group(3)), "display": bool(m.group(2)), "line": 0})
    # source lines: search each formula's TeX in the text, in order
    pos = 0
    for d in out:
        k = text.find(d["tex"].split("\n")[0], pos)
        if k >= 0:
            d["line"] = text.count("\n", 0, k) + 1
            pos = k
    return out


def render(text: str, link: Optional[LinkResolver] = None, image: Optional[ImageResolver] = None) -> Rendered:
    """Render Markdown ``text`` to safe HTML.

    ``link(href)`` returns the attributes of the ``<a>`` element for a link target (``href`` must
    start with ``#``, ``docs/`` (or ``/docs/``), ``http(s)://`` or ``mailto:``; ``data-doc`` / ``data-anchor`` let
    the browser load another document); None renders the link text only.  ``image(src)`` returns
    the URL of an image (None shows ``[image: alt]``).  Default: same-page anchors and external links,
    no images.
    """
    r = _Renderer(link, image)
    lines = text.replace("\r\n", "\n").replace("\r", "\n").expandtabs(4).split("\n")
    html = r.blocks(lines)
    return Rendered(html=html, title=r.title, toc=r.toc)
