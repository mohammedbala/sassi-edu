"""Help pages: the safe Markdown renderer (sassi.ui.markdown) and the documentation set
(sassi.ui.helpdocs) shown by Help > Help (requirements 5.3 Help (F1), UI-08)."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from sassi.ui.helpdocs import HelpDocs, HelpError, image_file
from sassi.ui.markdown import escape_text, render, slugify

ROOT = Path(__file__).resolve().parents[2]


# ------------------------------------------------------------------ blocks
def test_headings_paragraphs_and_slugs():
    r = render("# Title\n\nSome *em* and **strong** text\nwrapped.\n\n## 1.3 A first analysis in five minutes\n\n"
               "## Repeat\n\n## Repeat\n")
    assert r.title == "Title"
    assert [t["slug"] for t in r.toc] == ["title", "13-a-first-analysis-in-five-minutes", "repeat", "repeat-1"]
    assert '<h2 data-slug="13-a-first-analysis-in-five-minutes">1.3 A first analysis in five minutes</h2>' in r.html
    assert "<p>Some <em>em</em> and <strong>strong</strong> text\nwrapped.</p>" in r.html
    assert " id=" not in r.html               # headings never get ids that could clash with the GUI page


def test_github_slugs():
    assert slugify("2.2 The module chain (manual Fig. 1.1)") == "22-the-module-chain-manual-fig-11"
    assert slugify("13. Common errors and how to fix them") == "13-common-errors-and-how-to-fix-them"
    assert slugify("VP-02a") == "vp-02a" and slugify("my_heading") == "my_heading"


def test_fenced_code_is_escaped_and_kept_verbatim():
    r = render("```python\nif a < b and c & d:\n    print('<b>')\n```\n\n~~~\n*not em* `x`\n~~~\n")
    assert '<pre><code class="language-python">if a &lt; b and c &amp; d:\n    print(\'&lt;b&gt;\')</code></pre>' in r.html
    assert "<pre><code>*not em* `x`</code></pre>" in r.html


def test_lists_nested_lazy_and_ordered():
    md = ("* **VP-A1** (`vp_ansys.py`) is converted\n  from a file. Checks:\n  * first check;\n  * second check,\n"
          "    continued;\n* Unit tests.\n\n3. third\n4. fourth\n   continued\n")
    h = render(md).html
    assert h.count("<ul>") == 2 and h.count("<li>") == 6
    assert "<li><strong>VP-A1</strong> (<code>vp_ansys.py</code>) is converted\nfrom a file. Checks:\n<ul>" in h
    assert "<li>second check,\ncontinued;</li>" in h
    assert '<ol start="3"><li>third</li><li>fourth\ncontinued</li></ol>' in h


def test_loose_list_and_paragraph_interruption():
    h = render("- a\n\n- b\n").html
    assert h == "<ul><li><p>a</p></li><li><p>b</p></li></ul>"
    # only an ordered item numbered 1 interrupts a paragraph: a wrapped '1994. The' stays text
    h = render("Published in\n1994. The paper\n").html
    assert h == "<p>Published in\n1994. The paper</p>"
    h = render("Steps:\n1. one\n2. two\n").html
    assert h == "<p>Steps:</p>\n<ol><li>one</li><li>two</li></ol>"


def test_blockquote_with_list_and_hr():
    h = render("> **Note.** text\n> more\n>\n> * a\n> * b\n\n---\n").html
    assert h.startswith("<blockquote>\n<p><strong>Note.</strong> text\nmore</p>\n<ul><li>a</li><li>b</li></ul>\n</blockquote>")
    assert h.endswith("<hr>")


def test_tables_alignment_escaped_pipes_and_inline():
    md = ("| VP | Quantity | Result |\n|---|:---:|---:|\n| [VP-01](#vp-01) | `VERIFY,[id\\|P0]` | **PASS** |\n"
          "| short |\n")
    h = render(md).html
    assert h.startswith('<div class="md-table"><table><thead><tr><th>VP</th><th class="al-center">Quantity</th>'
                        '<th class="al-right">Result</th></tr></thead>')
    assert '<td><a href="#vp-01" data-anchor="vp-01">VP-01</a></td>' in h
    assert '<td class="al-center"><code>VERIFY,[id|P0]</code></td>' in h
    assert "<tr><td>short</td><td class=\"al-center\"></td><td class=\"al-right\"></td></tr>" in h


# ------------------------------------------------------------------ inline and safety
def test_inline_code_escapes_entities_and_breaks():
    h = render("`a<b` and \\*literal\\* and G\\* and &lt;maxit> and a & b  \nnext line\\\nthird").html
    assert "<code>a&lt;b</code>" in h and "*literal*" in h and "G*" in h
    assert "&lt;maxit&gt;" in h and "a &amp; b" in h
    assert h.count("<br>") == 2


def test_no_emphasis_inside_words():
    h = render("K*s - w Ms, FILE_8 and snake_case_name, 2*3*4").html
    assert "<em>" not in h and "K*s" in h and "FILE_8" in h and "2*3*4" in h


def test_raw_html_and_unsafe_links_are_neutralised():
    md = ('<script>alert(1)</script> <img src=x onerror=alert(1)> [x](javascript:alert(1)) '
          '[y](data:text/html,hi) ![z](javascript:1) [ok](https://example.org "T") <https://a.b/c?d=1&e=2>')
    h = render(md).html
    assert "<script" not in h and "<img" not in h and "javascript:" not in h.replace('title="javascript:alert(1)"', "") \
        .replace('title="javascript:1"', "")
    assert "&lt;script&gt;" in h and "onerror=alert(1)&gt;" in h
    assert '<a href="https://example.org" target="_blank" rel="noopener noreferrer" title="T">ok</a>' in h
    assert '<a href="https://a.b/c?d=1&amp;e=2" target="_blank" rel="noopener noreferrer">https://a.b/c?d=1&amp;e=2</a>' in h
    assert '[image: z]' in h
    # a resolver cannot smuggle an unsafe href either
    h = render("[a](x.md)", link=lambda href: {"href": "javascript:alert(1)"}).html
    assert "<a" not in h
    # attribute values are escaped
    h = render('[a](https://e.org/"onmouseover="x)').html
    assert 'onmouseover="x' not in h


def test_formula_brackets_keep_their_text():
    h = render("|d_c[a+a0](T) - d_c[a](T)|").html
    assert "d_c[a+a0](T) - d_c[a](T)" in h


def test_bare_urls_and_trailing_punctuation():
    h = render("See https://doi.org/10.1/x. Done (https://e.org/a_(b)).").html
    assert '<a href="https://doi.org/10.1/x" target="_blank" rel="noopener noreferrer">' in h
    assert '>https://e.org/a_(b)</a>).' in h


def test_escape_text_keeps_entities():
    assert escape_text("&lt; & &amp; &#60; &#x3C; <") == "&lt; &amp; &amp; &#60; &#x3C; &lt;"


# ------------------------------------------------------------------ the documentation set
@pytest.fixture(scope="module")
def docs():
    return HelpDocs(ROOT)


def test_catalogue_lists_the_manuals(docs):
    ids = [d["id"] for d in docs.catalogue()]
    for want in ("docs/index.md", "docs/user/USER_GUIDE.md", "docs/user/GUI.md", "docs/user/ANSYS.md",
                 "docs/theory/THEORY_MANUAL.md", "docs/verification/VERIFICATION_MANUAL.md",
                 "docs/reference/COMMAND_REFERENCE.md", "examples/README.md"):
        assert want in ids, want
    assert ids[0] == "docs/index.md"
    assert not any(i.startswith("docs/internal/") for i in ids)
    groups = [d["group"] for d in docs.catalogue()]
    assert groups.index("User guides") < groups.index("Theory") < groups.index("Verification") < groups.index("Examples")


def test_every_listed_document_renders_and_its_links_resolve(docs):
    """Every link between documents of the help set reaches an existing document and heading."""
    for d in docs.catalogue():
        doc = docs.document(d["id"])
        assert doc["html"] and doc["toc"] and doc["title"], d["id"]
        for m in re.finditer(r"<a ([^>]*)>", doc["html"]):
            a = dict(re.findall(r'([\w-]+)="([^"]*)"', m.group(1)))
            assert re.match(r"^(#|docs/|https?://|mailto:)", a["href"]), a      # docs/: relative (web version)
            target = a.get("data-doc", d["id"])
            if "data-doc" in a:
                assert docs.is_readable(target), (d["id"], target)
            if "data-anchor" in a:
                slugs = {t["slug"] for t in docs.document(target)["toc"]}
                assert a["data-anchor"] in slugs, (d["id"], target, a["data-anchor"])


def test_index_links_to_the_manuals_and_spec(docs):
    h = docs.document("docs/index.md")["html"]
    assert 'data-doc="docs/user/USER_GUIDE.md"' in h and 'data-doc="docs/theory/THEORY_MANUAL.md"' in h
    assert 'data-doc="examples/README.md"' in h and 'data-doc="docs/spec/00_requirements.md"' in h
    assert 'href="#13-a-first-analysis-in-five-minutes" data-doc="docs/user/USER_GUIDE.md" ' \
           'data-anchor="13-a-first-analysis-in-five-minutes"' in h


def test_verification_figure_is_served_from_docs(docs, tmp_path):
    h = docs.document("docs/verification/VERIFICATION_MANUAL.md")["html"]
    m = re.search(r'<img src="(docs/verification/figures/[^"]+\.png)"', h)       # relative URL (web version)
    assert m, "Figure 1 of the Verification Manual"
    rel = m.group(1)
    assert image_file(ROOT, rel) is not None
    assert image_file(ROOT, "docs/../README.md") is None and image_file(ROOT, "docs/index.md") is None
    assert image_file(ROOT, "sassi/ui/static/app.js") is None


@pytest.mark.parametrize("name", ["../etc/passwd", "/etc/passwd", "docs/../../x.md", "docs/internal/wave3_packages.md",
                                  "sassi/ui/api.py", "docs/user", "docs/missing.md", "javascript:x"])
def test_documents_outside_the_help_set_are_refused(docs, name):
    with pytest.raises(HelpError):
        docs.document(name)


def test_search_finds_sections(docs):
    r = docs.search("Line Selection")
    assert r["results"] and all("line selection" in x["snippet"].lower() for x in r["results"])
    gui = [x for x in r["results"] if x["doc"] == "docs/user/GUI.md"]
    assert gui and all(x["slug"] for x in gui)
    assert docs.search("x")["results"] == []


def test_a_private_tree(tmp_path):
    (tmp_path / "docs" / "user").mkdir(parents=True)
    (tmp_path / "docs" / "img").mkdir()
    (tmp_path / "docs" / "index.md").write_text("# Home\n\n[Guide](user/G.md#two) ![f](img/a.png) ![g](img/b.svg)\n")
    (tmp_path / "docs" / "user" / "G.md").write_text("# Guide\n\n## Two\n\n[back](../index.md) [dir](../img/)\n")
    (tmp_path / "docs" / "img" / "a.png").write_bytes(b"\x89PNG\r\n")
    (tmp_path / "docs" / "img" / "b.svg").write_text("<svg/>")
    h = HelpDocs(tmp_path)
    assert [d["id"] for d in h.catalogue()] == ["docs/index.md", "docs/user/G.md"]
    home = h.document("docs/index.md")["html"]
    assert 'href="#two" data-doc="docs/user/G.md" data-anchor="two"' in home
    assert '<img src="docs/img/a.png" alt="f" loading="lazy">' in home and "[image: g]" in home
    g = h.document("docs/user/G.md")["html"]
    assert 'data-doc="docs/index.md"' in g and '<span class="md-nolink" title="../img/">dir</span>' in g
    (tmp_path / "docs" / "user" / "G.md").write_text("# Guide v2\n")
    assert h.document("docs/user/G.md")["title"] == "Guide v2"            # re-rendered when the file changes


def test_latex_math_is_kept_verbatim_and_escaped():
    from sassi.ui.markdown import render
    h = render('A $\\frac{a_1}{b_2} < c^*$ and $x<script>alert(1)</script>$ "q".\n\n$$\nK_{ff}^{-1} & X\n$$\n').html
    # emphasis and escapes are not applied inside a formula; < > & " are escaped in the attribute and the text
    assert 'data-tex="\\frac{a_1}{b_2} &lt; c^*"' in h and "<em>" not in h
    assert "<script>" not in h and "&lt;script&gt;" in h
    assert '<div class="math-display" data-tex="K_{ff}^{-1} &amp; X">' in h
    # a $ that does not open a formula stays text; tables keep their cells
    assert render("a $ 5 fee and 10$, US$5").html == "<p>a $ 5 fee and 10$, US$5</p>"
    t = render("| a | b |\n|---|---|\n| $x_1$ | $\\lvert y \\rvert$ |\n").html
    assert t.count('class="math"') == 2 and t.count("<td>") == 2
