"""Unit tests of the verification-manual and command-reference generator (sassi.verify.report).

The run/report machinery is tested with small synthetic problems (not registered in the global
registry), so the tests are fast and independent of the physics VPs.  One end-to-end test runs a
real, fast VP (VP-47, line mathematics) through the module command line and the VERIFYREPORT
command.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import pytest

from sassi.verify import Problem, VPResult
from sassi.verify import report as rep


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def _p(pid, func, tier="P0", modules=("SITE",), source="R2 C.1, C.2; D-ANL-07", slow=False, title=None):
    return Problem(pid, title or f"synthetic problem {pid}", func, tier, list(modules), source, slow)


def _vp_pass(workdir):
    """Synthetic passing problem: one rel, one abs, one condition and one informative check."""
    r = VPResult()
    r.check("peak |H|, beta = 5 %", 10.0125, 10.0125, rtol=1e-6)
    r.check("H_u - 1", 1e-4, 0.0, atol=1e-3, note="abs check")
    r.require("F symmetric", True)
    r.inform("Veletsos-Verbic c_r (superseded)", 0.30, 0.20)
    r.notes.append("a note with <angle> brackets and a | pipe")
    return r


def _vp_fail(workdir):
    r = VPResult()
    r.check("resonance A0", 0.55, 0.50, rtol=0.02)
    return r


def _vp_raise(workdir):
    raise RuntimeError("deliberate failure inside the VP")


def _vp_figure(workdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(1, 1), dpi=20)
    ax.plot([0, 1], [0, 1])
    p = Path(workdir) / "tiny.png"
    fig.savefig(p)
    plt.close(fig)
    r = VPResult()
    r.check("x", 1.0, 1.0, rtol=1e-9)
    r.figures.append(str(p))
    return r


def _table_rows_consistent(md: str) -> None:
    """Every Markdown table of ``md`` has the same number of cells in each row as in its header."""
    for header, rows in rep.parse_tables(md):
        for row in rows:
            assert len(row) == len(header), (header, row)


# --------------------------------------------------------------------------------------
# ids and selection
# --------------------------------------------------------------------------------------
def test_vp_sort_key_natural_order():
    ids = ["VP-E1", "VP-38T", "VP-02b", "VP-10", "VP-02a", "VP-A1", "VP-38", "VP-01", "VP-39b", "VP-39"]
    assert sorted(ids, key=rep.vp_sort_key) == ["VP-01", "VP-02a", "VP-02b", "VP-10", "VP-38", "VP-38T",
                                                 "VP-39", "VP-39b", "VP-A1", "VP-E1"]


def test_normalise_and_resolve_ids():
    assert rep.normalise_id("1") == "VP-01"
    assert rep.normalise_id("vp-2a") == "VP-02a"
    assert rep.normalise_id("a1") == "VP-A1"
    reg = {k: _p(k, _vp_pass) for k in ("VP-01", "VP-02a", "VP-A1")}
    assert rep.resolve_ids(["1", "vp-02A", "VP-a1", "VP-01"], reg) == ["VP-01", "VP-02a", "VP-A1"]
    with pytest.raises(KeyError):
        rep.resolve_ids(["VP-99"], reg)


def test_select_problems_tiers_only_and_slow():
    reg = {"VP-01": _p("VP-01", _vp_pass), "VP-02": _p("VP-02", _vp_pass, tier="P1"),
           "VP-03": _p("VP-03", _vp_pass, slow=True)}
    sel = rep.select_problems(reg, skip_slow=True)
    assert [p.id for p, _ in sel] == ["VP-01", "VP-02", "VP-03"]
    assert [bool(s) for _, s in sel] == [False, False, True]
    assert [p.id for p, _ in rep.select_problems(reg, tiers=["p1"])] == ["VP-02"]
    assert [p.id for p, _ in rep.select_problems(reg, only=["3", "VP-01"])] == ["VP-01", "VP-03"]


# --------------------------------------------------------------------------------------
# Markdown helpers
# --------------------------------------------------------------------------------------
def test_md_escape_keeps_code_spans_and_escapes_specials():
    s = rep.md_escape("k* = k c(beta), <nl> and `G* <x>` and beta_s, _x_")
    assert "`G* <x>`" in s                     # code span untouched
    assert "k\\* = k" in s and "&lt;nl>" in s     # only a tag-like '<' is escaped
    assert "beta_s" in s and "\\_x\\_" in s    # intraword underscore kept, boundary ones escaped
    assert rep.md_escape("written *by hand* as decks; k*/(k* - m)") == "written *by hand* as decks; k\\*/(k\\* - m)"


def test_md_table_escapes_pipes_and_keeps_column_count():
    t = rep.md_table(["a", "b"], [["|H_u - 1|", "x"], ["`a|b`", "y"]])
    _table_rows_consistent(t)
    assert "\\|H_u - 1\\|" in t and "`a\\|b`" in t


def test_split_row_keeps_escaped_pipes():
    assert rep.split_row("| VP-15 | \\|H_u − 1\\| < 1 % | x |") == ["VP-15", "\\|H_u − 1\\| < 1 %", "x"]


def test_rst_to_md_roles_bullets_and_literal_blocks():
    doc = ("Model built with :func:`sassi.elements.assemble` and ``k*``.\n\n"
           "* VP-01   first item\n  continued\n* second item\n\nExample::\n\n    N,1,0,0,0\n    N,2,1,0,0\n")
    md = rep.rst_to_md(doc)
    assert "`sassi.elements.assemble`" in md and "`k*`" in md
    assert "- VP-01 first item continued" in md and "- second item" in md
    assert "```text\nN,1,0,0,0\nN,2,1,0,0\n```" in md
    assert "Example:" in md and "Example::" not in md


def test_fmt_num_and_small():
    assert rep.fmt_num(10.0125) == "10.0125"
    assert rep.fmt_num(1.5e-7) == "1.5000e-07"
    assert rep.fmt_num(0.0) == "0"
    assert rep.fmt_num(float("nan")) == "nan"
    assert rep.fmt_small(1.234567e-3) == "0.00123"


def test_slug_matches_heading_anchors():
    assert rep.slug("VP-02a") == "vp-02a"
    assert rep.slug("4. Lead decisions and superseded references") == "4-lead-decisions-and-superseded-references"


# --------------------------------------------------------------------------------------
# specification lookups (the project's own requirements, R1 and R2)
# --------------------------------------------------------------------------------------
def test_requirement_rows_and_base_ids():
    rows = rep.requirement_rows()
    assert "VP-11" in rows and "±10 %" in rows["VP-11"]["Tolerance"]
    assert rows["VP-47"]["section"] == "6.4"
    row, rel = rep.requirement_row_for("VP-39b", rows)        # no own row: the base id's row
    assert row is rows["VP-39"] and rel == "extends"
    row, rel = rep.requirement_row_for("VP-38T", rows)        # own row in section 6.3 since wave 3
    assert row is rows["VP-38T"] and rel == ""
    assert rep.requirement_row_for("VP-02b", rows)[0] is rows["VP-02b"]
    assert rep.requirement_row_for("VP-A1", rows) == (None, "none")


def test_decisions_and_lead_sections():
    dec = rep.parse_decisions()
    assert dec["D-W2-01"].subject.startswith("VP-11") and "VP-11" in dec["D-W2-01"].evidence
    assert "Half-space sublayer law" in dec["D-SIT-02"].subject      # id cell with a bold suffix
    lead = rep.lead_decision_sections()
    assert [h.split()[0] for h, _, _ in lead] == ["7.16", "7.17", "7.18", "7.19", "7.20"]
    assert any(d.id == "D-W3-05" for d in lead[2][2])
    assert any(d.id == "D-W5-12" for d in lead[3][2])          # built-in inputs and defaults (wave 5)
    assert [d.id for d in lead[4][2]] == ["D-W6-01", "D-W6-02", "D-W6-03", "D-W6-04"]   # the soil island picture
    assert rep.lead_sections_text(["7.16", "7.17"]) == "7.16 and 7.17"
    assert rep.lead_sections_text(["7.16", "7.17", "7.18"]) == "7.16 to 7.18"
    assert all(len(rows) >= 10 for _, _, rows in lead[:4])
    assert "informative" in lead[1][1]


def test_source_citations_resolve_r1_r2_and_decisions():
    si = rep.SourceIndex.load()
    labels = [lab for lab, _ in si.cite("R2 C.1, C.2; D-ANL-07")]
    assert labels == ["R2 C.1", "R2 C.2", "D-ANL-07"]
    cites = dict(si.cite("R2 A.2, R1 §2.4/V4, D-SIT-02"))
    assert "Halfspace" in cites["R1 §2.4"] and "R1 V4" in cites and "R2 A.2" in cites
    assert dict(si.cite("R2 section 0 (Ostadan et al. 2004)"))["R2 §0"].startswith("0 Global conventions")
    assert "Kinematic interaction" in dict(si.cite("R2 E6"))["R2 E6"]
    assert si.cite("spec 05c C3") == []


# --------------------------------------------------------------------------------------
# running and the manual
# --------------------------------------------------------------------------------------
def _fake_runs(tmp_path, figdir=None):
    sel = [(_p("VP-11", _vp_pass, modules=("SITE", "POINT", "ANALYS")), ""),
           (_p("VP-13", _vp_fail, tier="P1", modules=("ANALYS", "UI")), ""),
           (_p("VP-90", _vp_raise, modules=("MOTION",), source=""), ""),
           (_p("VP-91", _vp_pass, slow=True), "slow (skipped with --skip-slow)"),
           (_p("VP-92", _vp_figure, modules=("CALCPAR",)), "")]
    seen = []
    runs = rep.run_vps(sel, workroot=tmp_path / "work", figdir=figdir,
                       progress=lambda k, n, r: seen.append((k, n, r.id)))
    assert seen == [(1, 5, "VP-11"), (2, 5, "VP-13"), (3, 5, "VP-90"), (4, 5, "VP-91"), (5, 5, "VP-92")]
    return runs


def test_run_vps_statuses_and_cleanup(tmp_path):
    runs = _fake_runs(tmp_path)
    assert [r.status for r in runs] == ["PASS", "FAIL", "ERROR", "SKIPPED", "PASS"]
    assert "deliberate failure" in runs[2].error
    assert runs[3].result is None and runs[3].skipped
    # work directories are removed after each problem (D-W2-11)
    assert not any((tmp_path / "work").iterdir())
    r = runs[0]
    assert len(r.criteria) == 3 and len(r.informative) == 1
    assert r.margin() == pytest.approx(0.1)          # abs check 1e-4 / 1e-3
    assert runs[1].margin() == pytest.approx(0.1 / 0.02)


def test_run_vps_keep_and_figures(tmp_path):
    figdir = tmp_path / "figs"
    sel = [(_p("VP-92", _vp_figure), "")]
    runs = rep.run_vps(sel, workroot=tmp_path / "work", keep=True, figdir=figdir, progress=None)
    assert (tmp_path / "work" / "VP-92" / "tiny.png").exists()            # kept
    assert len(runs[0].figures) == 1 and Path(runs[0].figures[0]).parent == figdir


def test_check_ratio_edge_cases():
    r = VPResult()
    r.check("zero tol exact", 1.0, 1.0, atol=0.0)
    r.check("zero tol miss", 1.1, 1.0, atol=0.0)
    r.check("nan", float("nan"), 1.0, rtol=0.1)
    ratios = [rep.check_ratio(c) for c in r.checks]
    assert ratios[0] == 0.0 and math.isinf(ratios[1]) and math.isinf(ratios[2])


def test_manual_markdown_structure(tmp_path):
    runs = _fake_runs(tmp_path, figdir=tmp_path / "out" / "figures")
    fig = rep.margin_figure(runs, tmp_path / "out" / "figures" / "vp_margins.png")
    assert fig is not None and fig.exists() and fig.stat().st_size < 300_000
    md = rep.manual_markdown(runs, figure=fig, out_dir=tmp_path / "out", generated="2026-01-01 00:00 UTC",
                             load_errors={"sassi.verify.problems.vp_broken": "ImportError('x')"})
    _table_rows_consistent(md)
    for heading in ("## 1. How the verification works", "## 2. Summary", "## 3. Module coverage matrix",
                    "## 4. Lead decisions and superseded references", "## 5. Verification problems",
                    "## 6. Reference sources", "## 7. Environment and run log"):
        assert heading in md
    for vp in ("VP-11", "VP-13", "VP-90", "VP-91", "VP-92"):
        assert f"\n### {vp}\n" in md
    assert "**Not passing in this run:** [VP-13](#vp-13), [VP-90](#vp-90)" in md
    assert "deliberate failure inside the VP" in md                   # traceback of the ERROR
    assert "*Not run: slow (skipped with --skip-slow).*" in md
    assert "vp_broken" in md                                          # import errors are listed
    assert "![Largest error/tolerance ratio of each verification problem](figures/vp_margins.png)" in md
    # coverage matrix: VP-11 passes and covers SITE, POINT, ANALYS; VP-13 fails and covers UI/tools
    row11 = next(l for l in md.splitlines() if l.startswith("| [VP-11](#vp-11) | ") and "| P |" in l)
    assert row11.count("| P ") == 3
    row13 = next(l for l in md.splitlines() if l.startswith("| [VP-13](#vp-13) | ") and "| F |" in l)
    assert row13.count("| F ") == 2
    # requirements plan text, lead decision link and the informative comparison table
    assert "**Purpose.** Rigid disk dynamic impedance" in md
    assert "D-W2-01" in md and "| [VP-11](#vp-11) | Veletsos-Verbic c\\_r (superseded) |" not in md
    info_line = next(l for l in md.splitlines() if "Veletsos-Verbic c_r (superseded)" in l and "D-W2-01" in l)
    assert "0.3" in info_line and "0.2" in info_line
    # checks table: condition and informative rows, notes escaped
    assert "| 3 | F symmetric | true | true |  | condition | pass |" in md
    assert "informative | info |" in md
    assert "&lt;angle> brackets and a | pipe" in md
    assert "* R2 C.1:" in md and "* D-ANL-07:" in md


def test_coverage_counts():
    runs = [rep.VPRun("VP-01", "a", modules=["SITE", "ANALYS"], result=VPResult()),
            rep.VPRun("VP-02", "b", modules=["READSTR", "UI"], skipped="slow")]
    runs[0].result.check("x", 1, 1, rtol=1e-9)
    cov = rep.coverage(runs)
    assert cov["SITE"]["vps"] == 1 and cov["SITE"]["PASS"] == 1
    assert cov[rep.OTHER_COLUMN]["vps"] == 1 and cov[rep.OTHER_COLUMN]["SKIPPED"] == 1
    assert rep.module_columns(["READSTR", "CSECT", "MOTION"]) == [rep.OTHER_COLUMN, "MOTION"]


def test_describe_problem_from_docstrings():
    def documented(workdir):
        """Disk R = 5 m, ``R0 = 0.85 h``; see :func:`sassi.core.axisym.solve`."""
    assert rep.describe_problem(_p("VP-77", documented)) == \
        "Disk R = 5 m, `R0 = 0.85 h`; see `sassi.core.axisym.solve`."
    doc = "Problems.\n\n* VP-01   fixed-base SDOF peak\n          second line\n* VP-15   surface mat\n"
    assert rep._module_paragraph(doc, "VP-01") == "fixed-base SDOF peak second line"
    assert rep._module_paragraph(doc, "VP-1") == ""
    assert rep._module_paragraph("VP-A1  A cantilever read from a .cdb file.\n\nChecks:", "VP-A1") == \
        "A cantilever read from a .cdb file."
    # a bullet that names the problem in bold, as a label
    doc2 = "Intro.\n\n* stages things;\n* defines **VP-E1**: example 1 through the\n  interpreter.\n\nWhy."
    assert rep._module_paragraph(doc2, "VP-E1") == "example 1 through the interpreter."
    # a line listing several problems is not a description; an id list with commas is skipped too
    assert rep._module_paragraph("VP-25 COMBIN merge, VP-28 interpolation, VP-29 identity.\n", "VP-25") == ""
    assert rep._module_paragraph("Problems: VP-38T, VP-54T.\n\n* VP-38T  thick plate.\n", "VP-38T") == "thick plate."
    assert rep._module_paragraph("* VP-47 (T-L1): A x = [0, 1]", "VP-47") == "(T-L1): A x = [0, 1]"


def test_generate_manual_with_injected_registry(tmp_path):
    reg = {"VP-01": _p("VP-01", _vp_pass), "VP-02": _p("VP-02", _vp_pass, slow=True)}
    out = tmp_path / "docs" / "VM.md"
    path, runs = rep.generate_manual(out, skip_slow=True, figures=True, progress=None, registry=reg)
    assert path == out and out.exists()
    assert (tmp_path / "docs" / "figures" / "vp_margins.png").exists()
    text = out.read_text()
    assert [r.status for r in runs] == ["PASS", "SKIPPED"]
    assert "Partial report" not in text                       # skip-slow keeps every problem listed
    _table_rows_consistent(text)
    path2, _ = rep.generate_manual(tmp_path / "VM2.md", only=["VP-02"], progress=None, registry=reg)
    assert "Partial report" in path2.read_text()


# --------------------------------------------------------------------------------------
# command reference
# --------------------------------------------------------------------------------------
def test_command_reference_covers_every_command(tmp_path):
    from sassi.prep import registry
    infos = rep.command_infos()
    names = [i.name for i in infos]
    assert len(names) == len(set(names)) == len(registry.all_commands())
    assert set(registry.CATALOGUE) <= set(names)
    by = {i.name: i for i in infos}
    assert by["N"].syntax == "N,<nd>,[<x>],[<y>],[<z>]" and by["N"].status == "yes"
    assert by["AFWRITE"].abbrev == "AFWR" and by["AFWRITE"].category == "3.4.A"
    assert by["RUNSOIL"].syntax.startswith("RUNSOIL")
    assert by["DELL"].syntax.startswith("DELL,")
    for i in infos:
        spec = registry.REGISTRY[i.name]
        assert (i.status == "yes") == (not spec.placeholder)
        assert i.syntax and i.tier in ("P0", "P1", "P2", "OOS")
    md = rep.command_reference_markdown(infos, generated="2026-01-01")
    _table_rows_consistent(md)
    rows = [l for l in md.splitlines() if l.startswith("| **") and not l.startswith("| **Total")]
    assert len(rows) == len(infos)
    out = rep.write_command_reference(tmp_path / "ref" / "COMMANDS.md")
    assert out.exists() and out.read_text().startswith("# SASSI-EDU Command Reference")


def test_doc_syntax_parsing():
    assert rep._doc_syntax("FCOPY", "FCOPY,<src>,<dst>: copy a file.") == ("FCOPY,<src>,<dst>", "copy a file.")
    assert rep._doc_syntax("ANALYS", "ANALYS,<opmode>,<type>.") == ("ANALYS,<opmode>,<type>", "")
    assert rep._doc_syntax("FIXSPRROT", "FIXSPRROT (alias FIXSPROT): fix DOFs.")[0] == "FIXSPRROT"
    assert rep._doc_syntax("X", "Something else.") == ("", "Something else.")
    assert rep._first_sentence("copy a file (e.g. FILE1 -> FILE1X). More.") == "Copy a file (e.g. FILE1 -> FILE1X)."


# --------------------------------------------------------------------------------------
# end to end with a real (fast) problem: module command line and VERIFYREPORT
# --------------------------------------------------------------------------------------
def test_main_runs_a_real_vp_and_lists(tmp_path, capsys):
    out = tmp_path / "VM.md"
    rc = rep.main(["--only", "VP-47", "--out", str(out), "--no-commands", "--quiet"])
    assert rc == 0
    text = out.read_text()
    assert "\n### VP-47\n" in text and "Partial report" in text and "| PASS |" in text
    _table_rows_consistent(text)
    assert rep.main(["--list"]) == 0
    assert "VP-47" in capsys.readouterr().out
    assert rep.main(["--only", "VP-999", "--out", str(tmp_path / "x.md"), "--no-commands", "--quiet"]) == 2


def test_partial_run_never_replaces_the_manual(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    manual = rep.DEFAULT_MANUAL
    before = manual.stat().st_mtime if manual.exists() else None
    assert rep.main(["--only", "VP-47", "--no-commands", "--quiet"]) == 0
    assert (tmp_path / rep.PARTIAL_REPORT).exists()
    assert (manual.stat().st_mtime if manual.exists() else None) == before


def test_verifyreport_command(tmp_path):
    from sassi.prep import Interpreter
    from sassi.prep.messages import Kind
    ui = Interpreter(cwd=tmp_path)
    ui.execute("N,1,0,0,0")
    before = ui.model.canonical() if hasattr(ui.model, "canonical") else None
    assert ui.execute("VERIFYREPORT,VP-47,rep.md")
    text = (tmp_path / "rep.md").read_text()
    assert text.startswith("# SASSI-EDU Verification Report") and "\n### VP-47\n" in text
    assert any("report written to" in t for t in ui.sink.texts(Kind.CONFIRM) + ui.sink.texts(Kind.INFO))
    assert not ui.execute("VERIFYREPORT,VP-999")                   # unknown problem: error, no file
    assert ui.sink.texts(Kind.ERROR)
    if before is not None:
        assert ui.model.canonical() == before                     # an action: the model is unchanged
    # WRITE -> INP round trip is unaffected (the command stores nothing)
    assert ui.execute("WRITE,m.pre")
    ui2 = Interpreter(cwd=tmp_path)
    ui2.execute("INP,m.pre")
    assert ui2.model.same_state(ui.model)
