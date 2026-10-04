"""AFWRBAT: split an SSI run by frequency subsets (manual section 9.7.2; spec 09 section 2.2; Q-G1).

``AFWRBAT,<splits>`` prepares a frequency-parallel run of the active model:

1. The SSI frequency set used by the analysis (SITE ``<freq>``) is split into ``<splits>`` contiguous
   blocks whose sizes differ by at most one.
2. For block k a folder ``<Path>/<Model>_<k>/`` receives the input decks of the frequency-dependent
   modules enabled by AOPT (SITE, POINT, HOUSE, FORCE, ANALYS) with only block k in the frequency
   set, a command file ``run_<k>.pre`` (``MDL,<Model>,.`` then ``RUN<MODULE>`` in run order) and two
   launchers ``run_<k>.sh`` / ``run_<k>.bat`` that change to the folder and run it.  The folders are
   self-contained and can be moved to other computers (paths inside are relative).
3. In the model folder AFWRBAT writes the decks of the post-processors enabled by AOPT (MOTION,
   STRESS, RELDISP) for the full frequency set, and ``combine.pre`` / ``combine.sh`` / ``combine.bat``.
   ``combine`` copies the FILE8 of every folder and folds them pairwise with COMBIN (which reads
   FILE81 + FILE82 and writes FILE8, requirements 4.7), once per solution file (FILE8; FILE8X/Y/Z for
   simultaneous X/Y/Z input; FILE8001... for load cases or incoherent simulations), copies the
   frequency-independent HOUSE files (FILE4 = <Model>.N4, COOSK, COOSM, DOFSMAP, FILE90, FILE91) from
   folder 1 when STRESS needs them, and then runs MOTION, STRESS and RELDISP.

EQUAKE and SOIL are frequency-independent free-field steps: run them before AFWRBAT; a FILE88 present
in the model folder (SITEX non-linear soil) is copied into every folder.

The launchers call the Python interpreter that is running SASSI-EDU now (``sys.executable``); edit
them when the folders are moved to a computer with another installation (as the manual notes for
the ACS SASSI installation path).
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path
from typing import List

from ...conventions import incoherent_file8_name, loadcase_name
from ..afwrite import AfwriteError, afwrite
from ..options import aopt_modules, get_record
from ..registry import command

#: frequency-dependent modules run in every folder (run order)
FOLDER_MODULES = ("SITE", "POINT", "HOUSE", "FORCE", "ANALYS")
#: post-processors run after the FILE8 files are combined (run order)
POST_MODULES = ("MOTION", "STRESS", "RELDISP")
#: frequency-independent HOUSE outputs copied from folder 1 for STRESS
HOUSE_FILES = ("COOSK", "COOSM", "DOFSMAP", "FILE90", "FILE91")


def split_contiguous(nums: List[int], k: int) -> List[List[int]]:
    """``k`` contiguous blocks of the sorted list ``nums`` whose sizes differ by at most one."""
    n = len(nums)
    base, extra = divmod(n, k)
    out, i = [], 0
    for b in range(k):
        size = base + (1 if b < extra else 0)
        out.append(nums[i:i + size])
        i += size
    return out


def solution_files(model) -> List[str]:
    """FILE8 names written by ANALYS for the model's options (single, X/Y/Z, load cases, simulations)."""
    an = get_record(model, "ANALYS")
    simul = int(an.get("simul") or 0)
    atype = int(an.get("type") or 0)
    coh = int(get_record(model, "HOUSE").get("coh") or 0)
    if atype == 1 and simul > 1:
        return [loadcase_name("FILE8", i) for i in range(1, simul + 1)]
    if atype == 0 and coh == 1 and simul > 1:
        return [incoherent_file8_name(s, d) for s in range(1, simul + 1) for d in (1, 2, 3)]
    if atype == 0 and simul >= 1:
        return ["FILE8X", "FILE8Y", "FILE8Z"]
    return ["FILE8"]


def _launchers(folder: Path, stem: str) -> None:
    py = sys.executable
    sh = folder / f"{stem}.sh"
    sh.write_text("#!/bin/sh\n# SASSI-EDU AFWRBAT launcher: run from anywhere\n"
                  f'cd "$(dirname "$0")" && "{py}" -m sassi.cli run {stem}.pre\n', encoding="utf-8")
    sh.chmod(0o755)
    (folder / f"{stem}.bat").write_text("@echo off\r\nrem SASSI-EDU AFWRBAT launcher\r\n"
                                        f'cd /d "%~dp0"\r\n"{py}" -m sassi.cli run {stem}.pre\r\n',
                                        encoding="utf-8")


@command("AFWRBAT", max_args=1)
def cmd_afwrbat(c):
    """AFWRBAT,<splits>: split the SSI frequency set into <splits> folders with run and combine scripts."""
    splits = c.int(1, required=True, what="splits")
    m = c.model
    if not m.name or not m.path:
        c.fail("model name and path are not defined -- use MDL,<Model>,<Path>")
    base = Path(m.path)
    if not base.is_dir():
        c.fail(f"model directory {base} does not exist")
    set_no = int(get_record(m, "SITE").get("freq") or 1)
    nums = sorted(set(int(n) for n in m.freq_sets.get(set_no, [])))
    if not nums:
        c.fail(f"frequency set {set_no} (SITE <freq>) is not defined")
    if not 1 <= splits <= len(nums):
        c.fail(f"<splits> must be 1..{len(nums)} (frequency set {set_no} has {len(nums)} frequencies)")
    enabled = aopt_modules(m)
    folder_mods = [x for x in FOLDER_MODULES if x in enabled]
    post_mods = [x for x in POST_MODULES if x in enabled]
    if "ANALYS" not in folder_mods:
        c.fail("AOPT does not enable ANALYS: nothing to split")
    blocks = split_contiguous(nums, splits)
    names = solution_files(m)
    dirs = [c.interp.cwd]
    for k, block in enumerate(blocks, start=1):
        folder = base / f"{m.name}_{k}"
        folder.mkdir(parents=True, exist_ok=True)
        mk = m.copy()
        mk.path = str(folder)
        mk.freq_sets[set_no] = list(block)
        try:
            res = afwrite(mk, dirs=dirs + [base], modules=folder_mods)
        except AfwriteError as exc:
            c.fail(f"folder {k}: {exc}")
        if res.blocked:
            c.fail(f"folder {k}: AFWRITE blocked {', '.join(sorted(res.blocked))} -- run CHECK and fix the model")
        if (base / "FILE88").is_file():
            shutil.copy2(base / "FILE88", folder / "FILE88")
        lines = [f"* AFWRBAT block {k} of {splits}: frequency numbers {block[0]}..{block[-1]} ({len(block)})",
                 f"MDL,{m.name},."] + [f"RUN{mod}" for mod in folder_mods]
        (folder / f"run_{k}.pre").write_text("\n".join(lines) + "\n", encoding="utf-8")
        _launchers(folder, f"run_{k}")
        c.info(f"AFWRBAT: {folder.name}: {len(block)} frequencies ({block[0]}..{block[-1]}), decks "
               f"{', '.join(sorted(res.written))}")
    if post_mods:
        try:
            res = afwrite(m, dirs=dirs, modules=post_mods)
        except AfwriteError as exc:
            c.fail(str(exc))
        if res.blocked:
            c.warn(f"AFWRITE blocked {', '.join(sorted(res.blocked))} in {base} (post-processing decks)")
    comb = [f"* AFWRBAT combine: fold the FILE8 files of {splits} folders with COMBIN (FILE81 + FILE82 -> FILE8)",
            f"MDL,{m.name},."]
    for name in names:
        first = f"{m.name}_1/{name}"
        if splits == 1:
            comb.append(f"FCOPY,{first},{name}")
            continue
        comb.append(f"FCOPY,{first},FILE81")
        for k in range(2, splits + 1):
            comb.append(f"FCOPY,{m.name}_{k}/{name},FILE82")
            comb.append("RUNCOMBIN")
            if k < splits:
                comb.append("FMOVE,FILE8,FILE81")
        if name != "FILE8":
            comb.append(f"FMOVE,FILE8,{name}")
    if "STRESS" in post_mods:
        comb.append(f"FCOPY,{m.name}_1/{m.name}.N4,{m.name}.N4")
        comb += [f"FCOPY,{m.name}_1/{f},{f}" for f in HOUSE_FILES]
    comb += [f"RUN{mod}" for mod in post_mods]
    (base / "combine.pre").write_text("\n".join(comb) + "\n", encoding="utf-8")
    _launchers(base, "combine")
    c.confirm(f"AFWRBAT: {splits} folders written in {base}; run run_<k>.sh (or .bat) in each folder, then "
              f"combine.sh in {base.name}")
