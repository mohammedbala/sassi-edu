"""Lexer of the .pre language: rules L1-L9 (requirements section 3.1), UT-02 (lexer part)."""
from __future__ import annotations

import pytest

from sassi.model.values import NumberError, fmt_num, parse_float, parse_int, round_half_away
from sassi.prep.lexer import (LexError, compress_ids, format_token, is_comment, join_command, parse_id_list,
                              split_args, split_head)


def test_split_head_comma_and_blank_separator():
    assert split_head("N,1,2,3") == ("N", "1,2,3")
    # L1 / D-PAR-02: a run of blanks may replace the first comma only
    assert split_head("EDGE 1,0,0,1") == ("EDGE", "1,0,0,1")
    assert split_head("  edge \t 1,0") == ("edge", "1,0")
    # blanks followed by the first comma: that comma is the separator (L2 trimming)
    assert split_head("N , 1, 2") == ("N", " 1, 2")
    assert split_head("NLIST") == ("NLIST", "")
    assert split_head("NLIST   ") == ("NLIST", "")


def test_blank_fields_are_kept_ut02():
    # UT-02: cutvol,3,,,,,2.53 -> blank fields take defaults
    name, rest = split_head("cutvol,3,,,,,2.53")
    lx = split_args(rest)
    assert name == "cutvol"
    assert lx.tokens == ["3", "", "", "", "", "2.53"]
    # whitespace-only fields are blank too, and tokens are trimmed (L2, L5)
    assert split_args(" 1 ,  , 3 ").tokens == ["1", "", "3"]
    assert split_args("").tokens == []
    assert split_args("1,").tokens == ["1", ""]


def test_text_from_takes_rest_of_line_ut02():
    # L7: last-position text keeps its commas
    lx = split_args("Demo surface mat, coherent, BS units", text_from=1)
    assert lx.tokens == ["Demo surface mat, coherent, BS units"]
    lx = split_args("3,Reactor building walls, east", text_from=2)
    assert lx.tokens == ["3", "Reactor building walls, east"]
    lx = split_args('1,0.0001,1.0,0.0001,0.4,Rock, hard', text_from=6)
    assert lx.tokens[-1] == "Rock, hard"
    # a fully quoted text loses its quotes
    assert split_args('"A, B"', text_from=1).tokens == ["A, B"]
    # fewer fields than text_from: ordinary split
    assert split_args("3", text_from=2).tokens == ["3"]


def test_quoted_tokens_embed_commas():
    lx = split_args('1,"a,b",c')
    assert lx.tokens == ["1", "a,b", "c"]
    assert lx.quoted == [False, True, False]
    with pytest.raises(LexError):
        split_args('1,"abc')
    # a quote inside a token is data
    assert split_args('ab"c,d').tokens == ['ab"c', "d"]


def test_legacy_parenthesised_token_ut02():
    # L9: a final (...) token is one token, inner commas do not split
    lx = split_args("4000,32.2,0,(5X,F10.4)", paren_last=True)
    assert lx.tokens == ["4000", "32.2", "0", "(5X,F10.4)"]
    assert lx.legacy_paren == "(5X,F10.4)"
    # without the option the commas split
    assert split_args("4000,(5X,F10.4)").tokens == ["4000", "(5X", "F10.4)"]


def test_comments():
    assert is_comment("* comment")
    assert is_comment("   *indented comment")
    assert not is_comment("N,1,2*3")      # no inline comments (L4)


def test_numbers_free_format_and_fortran_exponents():
    assert parse_float("1.0D9") == 1.0e9
    assert parse_float("1e+008") == 1.0e8
    assert parse_float(" -2.5d-3 ") == -2.5e-3
    assert parse_float(".5") == 0.5
    with pytest.raises(NumberError):
        parse_float("abc")
    with pytest.raises(NumberError):
        parse_float("")
    # D-PAR-07: integer fields given as reals are rounded; inexact flagged
    assert parse_int("5.0") == (5, True)
    assert parse_int("5.5") == (6, False)
    assert parse_int("2.5") == (3, False)          # half away from zero, not banker's rounding
    assert round_half_away(-2.5) == -3


def test_fmt_num_is_exact_and_compact():
    for x in (0.1, 1 / 3, 13.535533905932738, 1e-17, -4.0e9, 2.0 ** 60):
        assert float(fmt_num(x)) == x
    assert fmt_num(10.0) == "10"
    assert fmt_num(-3) == "-3"
    assert fmt_num(0.15) == "0.15"


def test_id_lists_and_ranges():
    assert parse_id_list(["1", "3-6 10;12"]) == [1, 3, 4, 5, 6, 10, 12]
    assert parse_id_list(["1 - 3"]) == [1, 2, 3]
    assert parse_id_list(["5\t7"]) == [5, 7]
    with pytest.raises(LexError):
        parse_id_list(["6-3"])          # D-PAR-08: descending range is an error
    with pytest.raises(LexError):
        parse_id_list(["a-b"])
    with pytest.raises(LexError):
        parse_id_list(["2.5"])
    assert compress_ids([1, 2, 3, 7, 9, 10, 11, 12]) == ["1-3", "7", "9-12"]
    assert compress_ids([4, 5]) == ["4", "5"]


def test_format_and_join():
    assert format_token("a,b") == '"a,b"'
    assert format_token(0.5) == "0.5"
    assert format_token(None) == ""
    assert join_command("TIT", ["x, y"], text_last=True) == "TIT,x, y"
    assert join_command("SITE", ["0", "1", "", ""]) == "SITE,0,1"


# ------------------------------------------------------------------ L7 / D-PAR-06: verbatim text
def test_text_argument_is_not_scanned_for_quotes():
    """Only the fields before a rest-of-line text argument are scanned (L7)."""
    lx = split_args('1,Walls, "north', text_from=2)
    assert lx.tokens == ["1", 'Walls, "north']
    lx = split_args('acc X, "raw', text_from=1)
    assert lx.tokens == ['acc X, "raw']
    lx = split_args('"Reactor building, 12 in. wall', text_from=1)
    assert lx.tokens == ['"Reactor building, 12 in. wall']
    # a text that is one quoted token is unquoted; an embedded pair of quotes is kept
    assert split_args('"A, B"', text_from=1).tokens == ["A, B"]
    assert split_args('say "hi", then', text_from=1).tokens == ['say "hi", then']


def test_text_argument_head_fields_still_quote_aware():
    lx = split_args('"a,b",x, "y', text_from=3)
    assert lx.tokens == ["a,b", "x", '"y'] and lx.quoted == [True, False, False]
    with pytest.raises(LexError):
        split_args('"a,b,x', text_from=3)            # the unbalanced quote is in a head field
    assert split_args("1", text_from=2).tokens == ["1"]          # fewer fields than text_from
    assert split_args("1,", text_from=2).tokens == ["1", ""]
    assert split_args("", text_from=1).tokens == []
