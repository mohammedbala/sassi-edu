"""``sassi`` console entry point (UI-02 b): --version, run, -c, the REPL loop."""
from __future__ import annotations

import io

import pytest

from sassi import __version__
from sassi.cli import main, repl
from sassi.prep import Interpreter


def test_version(capsys):
    with pytest.raises(SystemExit) as e:
        main(["--version"])
    assert e.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_run_file_exit_codes(tmp_path, capsys):
    good = tmp_path / "good.pre"
    good.write_text("* demo\nN,1,0,0,0\nN,2,1,0,0\nFILL\n")
    assert main(["--cwd", str(tmp_path), "run", "good.pre"]) == 0
    out = capsys.readouterr().out
    assert "* demo" in out and "INPUT FILE REACHED EOF" in out
    bad = tmp_path / "bad.pre"
    bad.write_text("N,1,0,0,0\nNOSUCH\n")
    assert main(["--cwd", str(tmp_path), "run", "bad.pre", "--quiet"]) == 1
    cap = capsys.readouterr()
    assert "NOSUCH Command not found" in cap.err and "1 errors" in cap.out
    assert main(["--cwd", str(tmp_path), "run", "missing.pre"]) == 2


def test_command_lines(tmp_path, capsys):
    assert main(["--cwd", str(tmp_path), "-c", "N,1,0,0,0", "-c", "NLIST"]) == 0
    assert "1 nodes listed" in capsys.readouterr().out
    assert main(["-c", "BOGUS"]) == 1


def test_repl_reads_until_exit(tmp_path, capsys):
    import sys
    ui = Interpreter(cwd=tmp_path)
    ui.sink.attach_stream(sys.stdout)
    assert repl(ui, stdin=io.StringIO("N,1,0,0,0\nNLIST\nexit\nN,2,0,0,0\n")) == 0
    assert set(ui.model.nodes) == {1}
    out = capsys.readouterr().out
    assert "1 nodes listed" in out and "> N,1,0,0,0" not in out      # no echo of typed lines


def test_command_line_inp_of_missing_file_fails(tmp_path, capsys):
    assert main(["--cwd", str(tmp_path), "-c", "INP,does_not_exist.pre"]) == 1
    assert "not found" in capsys.readouterr().err
