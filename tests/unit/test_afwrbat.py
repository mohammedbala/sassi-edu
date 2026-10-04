"""AFWRBAT (frequency-split runs): unit tests and an end-to-end equivalence check."""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from sassi.io.container import read_container
from sassi.prep import Interpreter
from sassi.prep.commands.afwrbat import split_contiguous

ROOT = Path(__file__).resolve().parents[2]


def test_split_contiguous_sizes_differ_by_at_most_one():
    nums = list(range(1, 12))
    for k in range(1, 12):
        blocks = split_contiguous(nums, k)
        assert sum(blocks, []) == nums
        sizes = [len(b) for b in blocks]
        assert max(sizes) - min(sizes) <= 1 and len(blocks) == k


def test_afwrbat_requires_mdl(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    assert not ui.execute("AFWRBAT,2")


@pytest.mark.slow
def test_split_run_equals_full_run(tmp_path):
    """ex01 run in full and as 3 frequency blocks combined with COMBIN: identical FILE8 (and MOTION output)."""
    ex = tmp_path / "examples"
    shutil.copytree(ROOT / "examples", ex)
    ui = Interpreter(cwd=ex)
    ui.run_file(ex / "ex01_surface_stick.pre")
    base = Path(ui.model.path)
    ref = read_container(base / "FILE8", "FILE8")
    ref_acc = (base / "00085TR_X.ACC").read_text()
    for f in base.glob("FILE8*"):
        f.unlink()
    (base / "00085TR_X.ACC").unlink()
    assert ui.execute("AFWRBAT,3")
    folders = sorted(p for p in base.iterdir() if p.is_dir() and p.name.startswith(ui.model.name + "_"))
    assert [p.name for p in folders] == [f"{ui.model.name}_{k}" for k in (1, 2, 3)]
    nblock = []
    for k, folder in enumerate(folders, start=1):
        assert (folder / f"run_{k}.sh").is_file() and (folder / f"run_{k}.bat").is_file()
        sub = Interpreter(cwd=folder)
        sub.run_file(folder / f"run_{k}.pre")
        part = read_container(folder / "FILE8", "FILE8")
        nblock.append(len(part["fnum"]))
    assert sum(nblock) == len(ref["fnum"]) and max(nblock) - min(nblock) <= 1
    top = Interpreter(cwd=base)
    top.run_file(base / "combine.pre")
    comb = read_container(base / "FILE8", "FILE8")
    assert np.array_equal(comb["fnum"], ref["fnum"])
    assert np.array_equal(comb["eq_node"], ref["eq_node"]) and np.array_equal(comb["eq_dof"], ref["eq_dof"])
    assert np.max(np.abs(comb["H"] - ref["H"])) <= 1e-12 * np.max(np.abs(ref["H"]))
    acc = np.array([float(x) for x in (base / "00085TR_X.ACC").read_text().split()])
    acc0 = np.array([float(x) for x in ref_acc.split()])
    assert np.allclose(acc, acc0, rtol=0, atol=1e-10 * np.max(np.abs(acc0)))
