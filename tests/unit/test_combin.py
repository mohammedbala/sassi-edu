"""Unit tests of the COMBIN module (requirements section 4.7, D-CMB-01/02)."""
from __future__ import annotations

import numpy as np

from sassi.io.container import read_container
from sassi.modules.base import run_module
from sassi.modules.combin import read_options
from sassi.verify.problems.vp_motion import listing_text, write_file8

DF = 0.048828125


def _H(fnum, neq=2):
    return np.outer(fnum, np.arange(1, neq + 1)) * (1 + 0.5j)


def test_merge_sorted_by_frequency_number(tmp_path):
    a, b = np.array([1, 4, 9, 20]), np.array([2, 3, 10, 30, 31])
    write_file8(tmp_path / "FILE81", a, DF, [5, 5], [1, 3], _H(a))
    write_file8(tmp_path / "FILE82", b, DF, [5, 5], [1, 3], _H(b))
    assert run_module("COMBIN", "c", tmp_path) == 0
    c = read_container(tmp_path / "FILE8", "FILE8")
    assert c["fnum"].tolist() == sorted(a.tolist() + b.tolist())
    assert np.array_equal(c["H"], _H(c["fnum"])) and np.allclose(c["freq"], c["fnum"] * DF)
    assert c.meta["module"] == "COMBIN" and c.meta["df"] == DF
    assert "FILE8 written: 9 frequencies" in listing_text(tmp_path, "c", "COMBIN")


def test_incompatible_files(tmp_path):
    a, b = np.array([1, 3]), np.array([2, 4])
    write_file8(tmp_path / "FILE81", a, DF, [5, 5], [1, 3], _H(a))
    for kw, msg in ((dict(eq_dof=[1, 2]), "DOF maps"), (dict(df=2 * DF), "frequency steps"), (dict(type=1), "'type'"),
                    (dict(cm=1), "'cm'"), (dict(ang=30.0), "angle"),
                    (dict(case=2), "different simultaneous / load cases")):      # D-ANL-06: FILE8001 vs FILE8002
        args = dict(eq_node=[5, 5], eq_dof=[1, 3], df=DF)
        args.update(kw)
        ex = {k: v for k, v in kw.items() if k in ("type", "cm", "ang", "case")}
        write_file8(tmp_path / "FILE82", b, args["df"], args["eq_node"], args["eq_dof"], _H(b), **ex)
        assert run_module("COMBIN", "c", tmp_path) == 1
        assert msg in listing_text(tmp_path, "c", "COMBIN")


def test_missing_input(tmp_path):
    write_file8(tmp_path / "FILE81", [1], DF, [1], [1], np.ones((1, 1)))
    assert run_module("COMBIN", "c", tmp_path) == 1
    assert "FILE82 missing -- run ANALYS" in listing_text(tmp_path, "c", "COMBIN")


def test_options_file(tmp_path):
    assert read_options(tmp_path)["COMBINDUP"] == "ERROR"
    (tmp_path / "COMBIN.opt").write_text("* comment\ncombindup prefer82\n")
    assert read_options(tmp_path)["COMBINDUP"] == "PREFER82"
    (tmp_path / "COMBIN.opt").write_text("EDUOPT,COMBINDUP,ERROR\n")
    assert read_options(tmp_path)["COMBINDUP"] == "ERROR"


def test_same_case_number_merges(tmp_path):
    """Two partial solutions of the same load case (case 3) merge; the case number is kept."""
    a, b = np.array([1, 3]), np.array([2, 4])
    write_file8(tmp_path / "FILE81", a, DF, [5, 5], [1, 3], _H(a), type=1, case=3)
    write_file8(tmp_path / "FILE82", b, DF, [5, 5], [1, 3], _H(b), type=1, case=3)
    assert run_module("COMBIN", "c", tmp_path) == 0
    assert int(read_container(tmp_path / "FILE8", "FILE8").meta["case"]) == 3
