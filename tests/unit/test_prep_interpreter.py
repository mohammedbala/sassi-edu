"""Interpreter core: command lookup (UT-01), lexical rules through the interpreter (UT-02), error
handling and EOF summaries (L12, L16), INP nesting and path resolution (L14, L15), placeholders of
lower-tier commands (section 0.1), the registry."""
from __future__ import annotations

import pytest

from sassi.prep import Interpreter, Kind
from sassi.prep.interpreter import EOF_MESSAGE
from sassi.prep.registry import ABBREVIATIONS, CATALOGUE, REGISTRY, command, lookup, register


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


# ------------------------------------------------------------------ UT-01
def test_ut01_every_full_name_and_abbreviation_resolves():
    for full, abbr in ABBREVIATIONS.items():
        assert lookup(full).name == full
        assert lookup(abbr).name == full, abbr
        assert lookup(abbr.lower()).name == full
        assert lookup(full.title()).name == full


def test_ut01_name_aliases_and_catalogue_complete():
    aliases = {"FIXSPROT": "FIXSPRROT", "PANELGEN": "PNLGEN", "NONLINMODISP": "NONLINMOTDISP",
               "ACCANIDB": "ACCDBANI", "LOADACCDBANI": "LOADACCDB", "COMDISPDB": "COMBDISPDB"}
    for a, full in aliases.items():
        assert lookup(a).name == full
    for name in CATALOGUE:
        assert lookup(name) is not None, name
    # no prefix matching; RELFILE has no abbreviation; UI commands need their full name
    for bad in ("ACC", "ANALY", "RELF", "WRI", "ANALYSI", "NGE", "INTGE", "CUTA"):
        assert lookup(bad) is None, bad
    # full names win over abbreviations
    assert lookup("GROUPMAT").name == "GROUPMAT" and lookup("GROU").name == "GROUP"
    assert lookup("ETYPEGEN").name == "ETYPEGEN" and lookup("ETYP").name == "ETYPE"


def test_ut01_not_found_and_processing_continues(ui):
    s = ui.run_text("ACC,1\nANALY,0\nRELF,x.tfi\nN,1,1,2,3\n")
    assert errors(ui) == ["ACC Command not found", "ANALY Command not found", "RELF Command not found"]
    assert 1 in ui.model.nodes                       # later lines still run
    assert s.errors == 3 and s.commands == 4


def test_ut01_names_case_insensitive_paths_case_kept(ui):
    ui.execute("thfile,/Data/Records/H1_Upper.ACC")
    ui.execute("tHtIt,Mixed Case Title")
    ui.execute("n,7,1,2,3")
    assert ui.model.options.string("THFILE") == "/Data/Records/H1_Upper.ACC"
    assert ui.model.options.string("THTIT") == "Mixed Case Title"
    assert 7 in ui.model.nodes


# ------------------------------------------------------------------ UT-02 through the interpreter
def test_ut02_blank_fields_titles_and_legacy_soil(ui):
    ui.execute("NGEN,,,,,,,,")                       # all blank, but no nodes: error only
    ui.execute("N,1,,,5")                            # blanks take defaults (0)
    assert ui.model.nodes[1].xyz == (0.0, 0.0, 5.0)
    ui.execute("TIT,Reactor building, best-estimate soil, run 3")
    assert ui.model.title == "Reactor building, best-estimate soil, run 3"
    ui.execute("SOIL,4000,32.2,0,1,1,8,0.65,1.0,(5X,F10.4)")
    rec = ui.model.options.record("SOIL")
    assert rec.legacy == "(5X,F10.4)"
    assert rec.values[-1] == "(5X,F10.4)"            # one token
    assert any("legacy" in w and "Error 102" in w for w in warnings(ui))


def test_required_numeric_blank_is_zero_with_warning(ui):
    ui.execute("N,")                                  # nd missing -> 0 with a warning, then refused
    assert any("missing, 0 used" in w for w in warnings(ui))
    assert errors(ui) == ["N: node numbers must be positive"]


def test_bad_number_is_an_error_and_continues(ui):
    s = ui.run_text("N,1,abc,0,0\nN,2,1,0,0\n")
    assert len(errors(ui)) == 1 and "not a number" in errors(ui)[0]
    assert 1 not in ui.model.nodes and 2 in ui.model.nodes
    assert s.errors == 1


def test_integer_field_given_as_real_rounds_with_warning(ui):
    ui.execute("N,3.0,1,1,1")
    assert 3 in ui.model.nodes and not warnings(ui)
    ui.execute("N,4.4,1,1,1")
    assert 4 in ui.model.nodes
    assert any("rounded to 4" in w for w in warnings(ui))


def test_extra_arguments_warn(ui):
    ui.execute("CSYS,0,5")
    assert any("ignored" in w for w in warnings(ui))
    ui.sink.clear()
    ui.execute("N,1,0,0,0,7,8")                       # legacy pile/control fields: silently accepted
    assert not warnings(ui)


def test_line_length_guard(ui):
    ui.execute("TIT," + "x" * 3100)
    assert any("3000" in e for e in errors(ui))


def test_comment_lines_and_echo_policy(ui, tmp_path):
    ui.execute("N,1,0,0,0")
    assert ui.sink.texts(Kind.ECHO) == ["N,1,0,0,0"]
    assert ui.sink.texts(Kind.CONFIRM)
    ui.sink.clear()
    (tmp_path / "a.pre").write_text("* a comment\nN,2,0,0,0\nBOGUS\n")
    ui.execute("INP,a.pre")
    # L16: echo and confirmations suppressed inside INP; comments, errors and info shown
    assert ui.sink.texts(Kind.ECHO) == ["INP,a.pre"]
    assert ui.sink.texts(Kind.CONFIRM) == []
    assert ui.sink.texts(Kind.COMMENT) == ["* a comment"]
    assert errors(ui) == ["BOGUS Command not found"]
    info = ui.sink.texts(Kind.INFO)
    assert info[-2].startswith("INP a.pre: 3 lines, 2 commands, 0 warnings, 1 errors")
    assert info[-1] == EOF_MESSAGE


def test_internal_errors_do_not_stop_the_session(ui):
    @command("ZZTESTBOOM", replaces=True)
    def boom(c):
        raise ZeroDivisionError("boom")
    try:
        assert ui.execute("ZZTESTBOOM") is False
        assert any("internal error" in e for e in errors(ui))
        assert ui.execute("N,1,0,0,0")
    finally:
        REGISTRY.pop("ZZTESTBOOM", None)


# ------------------------------------------------------------------ INP nesting and paths (L14, L15)
def test_nested_inp_and_relative_paths(ui, tmp_path):
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "inner.pre").write_text("N,2,0,0,0\nBAD\n")
    # inner.pre is found relative to the calling file's folder (last resort of L15)
    (sub / "outer.pre").write_text("N,1,0,0,0\nINP,inner.pre\nN,3,0,0,0\n")
    s = ui.run_file("sub/outer.pre")
    assert set(ui.model.nodes) == {1, 2, 3}
    assert s.errors == 1 and s.commands == 5
    info = ui.sink.texts(Kind.INFO)
    assert info.count(EOF_MESSAGE) == 1 and info[-1] == EOF_MESSAGE
    assert any(t.startswith("INP inner.pre: 2 lines") for t in info)


def test_inp_model_path_has_priority(ui, tmp_path):
    (tmp_path / "m").mkdir()
    (tmp_path / "x.pre").write_text("N,1,0,0,0\n")
    (tmp_path / "m" / "x.pre").write_text("N,2,0,0,0\n")
    ui.execute("MDL,demo,m")
    ui.execute("INP,x.pre")
    assert set(ui.model.nodes) == {2}


def test_inp_missing_file_and_recursion_guard(ui, tmp_path):
    ui.execute("INP,nothere.pre")
    assert any("not found" in e for e in errors(ui))
    (tmp_path / "loop.pre").write_text("INP,loop.pre\n")
    ui.sink.clear()
    ui.execute("INP,loop.pre")
    assert any("nesting deeper than 32" in e for e in errors(ui))


def test_windows_separators_accepted(ui, tmp_path):
    (tmp_path / "d").mkdir()
    (tmp_path / "d" / "w.pre").write_text("N,9,0,0,0\n")
    ui.execute(r"INP,.\d\w.pre")
    assert 9 in ui.model.nodes


# ------------------------------------------------------------------ placeholders (tier messages)
def test_lower_tier_action_prints_not_available(ui):
    ui.execute("SETENV,X,1")                          # still a P2 placeholder (FILLPOOL is implemented)
    assert warnings(ui) == ["SETENV is not available in this build (tier P2)"]


def test_former_storage_placeholder_now_has_a_real_handler(ui):
    # NLSLAYER was a P2 storage placeholder; SOIL-NON implements it (wave 3): no tier message any more
    ui.execute("NLSLAYER,3,1,0.5")
    assert not any("not available in this build" in w for w in warnings(ui))


def test_p0_storage_placeholders_are_silent(ui):
    ui.execute("SITE,0,1,0,20,3,1,0,1,2048,1,0,0.005,4096,1")
    ui.execute("WAVE,2,1,1.0,1.0,0.0")
    ui.execute("EDUOPT,hslaw,UNIFORM")
    ui.execute("BINOUT,1,1,0")
    ui.execute("BINOUT,,0")                           # BINOUT: blank = unchanged
    ui.execute("SYMM,1,0,1,2,5")
    ui.execute("SYMM,2,1,3,4,6")
    ui.execute("SYMM,1,0,0")                          # node1 = 0 resets plane 1
    assert not warnings(ui) and not errors(ui)
    o = ui.model.options
    assert o.record("SITE").integer(13) == 4096
    assert o.entry("WAVE", 2).number(3) == 1.0
    assert o.entry("EDUOPT", "HSLAW").values == ["hslaw", "UNIFORM"]
    assert o.record("BINOUT").values == ["1", "0", "0"]
    assert [k for k, _ in o.entries("SYMM")] == [2]
    ui.execute("DYNP,1,0.0001,1.0,0.0001,0.4,Rock")
    ui.execute("DYNP,1,0.0001,0.9,0.0001,0.4,Sand")
    assert [k for k, _ in o.entries("DYNP")] == [("Rock", 1), ("Sand", 1)]


# ------------------------------------------------------------------ registry
def test_registry_duplicate_and_override():
    def h1(c):
        pass

    def h2(c):
        pass
    register("ZZTESTDUP", h1)
    try:
        with pytest.raises(ValueError):
            register("ZZTESTDUP", h2)
        spec = register("ZZTESTDUP", h2, replaces=True)
        assert spec.handler is h2
        # a placeholder never replaces a real handler
        assert register("ZZTESTDUP", h1, placeholder=True).handler is h2
        with pytest.raises(ValueError):
            register("ZZTESTAB", h1, abbrev=("NGEN",))        # abbreviation = an existing name
        with pytest.raises(ValueError):
            register("ZZTESTAB", h1, abbrev=("NSCA",))        # abbreviation used by NSCALE
    finally:
        REGISTRY.pop("ZZTESTDUP", None)
        REGISTRY.pop("ZZTESTAB", None)


def test_catalogue_defaults_apply_to_handlers():
    assert lookup("TIT").text_from == 1 and lookup("TIT").placeholder is False
    assert lookup("FOREACH").raw is True
    assert lookup("SOIL").paren is True
    assert lookup("NGEN").abbrev == () and lookup("NSCALE").abbrev == ("NSCA",)
    assert lookup("NGEN").module.endswith("commands.nodes")


# ------------------------------------------------------------------ max_args is a hard limit
def test_extra_arguments_are_really_ignored(ui):
    """The 'ignored' warning and the handler agree: arguments after max_args are dropped."""
    ui.run_text("\n".join([f"N,{i},{4 * (i - 1)},0,0" for i in range(1, 9)] + ["N,9,1000,0,0"]))
    ui.execute("NMED,500,1,2,3,4,5,6,7,8,9")           # spec 08 3.14: 1 to 8 nodes
    assert ui.model.nodes[500].x == pytest.approx(14.0)
    assert sum("ignored" in w for w in warnings(ui)) == 1
    ui.sink.clear()
    ui.execute("D,1,1,1,1,UX,UY,UZ,ROTX,ROTY,ROTZ,ALL")    # six labels at most: ALL is dropped
    assert ui.model.nodes[1].fix == [1] * 6 and any("ignored" in w for w in warnings(ui))
    ui.execute("D,2,2,1,1,UX,UX,UX,UX,UX,UX,ROTZ")
    assert ui.model.nodes[2].fix == [1, 0, 0, 0, 0, 0]


def test_lmove_takes_at_most_15_list_nodes(ui):
    """Spec 08 Q6: LMOVE/NMOVE accept up to 15 list nodes; a 16th is dropped with one warning."""
    ui.run_text("\n".join(f"N,{i},{i},0,0" for i in range(1, 17)))
    ui.sink.clear()
    ui.execute("LMOVE,0,0,1,101," + ",".join(str(i) for i in range(1, 17)))
    new = sorted(i for i in ui.model.nodes if i > 100)
    assert new == list(range(101, 116))
    w = warnings(ui)
    assert len(w) == 1 and "ignored" in w[0]


# ------------------------------------------------------------------ files with a byte-order mark
def test_inp_and_run_text_ignore_utf8_bom(ui, tmp_path):
    """Spec 04 4.1: .pre files saved with a UTF-8 BOM (Windows editors) load unchanged."""
    (tmp_path / "bom.pre").write_bytes("* first line comment\r\nN,1,1,2,3\r\n".encode("utf-8-sig"))
    (tmp_path / "bom2.pre").write_bytes("N,2,4,5,6\n".encode("utf-8-sig"))
    ui.execute("INP,bom.pre")
    ui.execute("INP,bom2.pre")
    assert sorted(ui.model.nodes) == [1, 2] and not errors(ui)
    assert ui.sink.texts(Kind.COMMENT) == ["* first line comment"]
    ui.run_text("﻿N,3")
    assert 3 in ui.model.nodes and not errors(ui)


# ------------------------------------------------------------------ INP of an unreadable file fails
def test_inp_missing_file_fails_with_one_message(ui):
    assert ui.execute("N,1") is True
    assert ui.execute("INP,missing.pre") is False
    assert ui.last_ok is False
    assert len(errors(ui)) == 1 and "missing.pre not found" in errors(ui)[0]
    assert ui.history == ["N,1"]                      # the failed INP is not replayed (L17)


def test_inp_with_errors_inside_still_succeeds(ui, tmp_path):
    """L12: errors inside a file do not stop it, and do not make the INP command fail."""
    (tmp_path / "e.pre").write_text("NOSUCH\nN,1\n")
    assert ui.execute("INP,e.pre") is True and 1 in ui.model.nodes
