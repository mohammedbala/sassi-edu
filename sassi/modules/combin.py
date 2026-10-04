"""COMBIN: merge two SSI solutions computed for different frequency sets into one FILE8.

Why (requirements section 4.7, spec 02 section 3.11)
---------------------------------------------------
Interpolation quality depends on the SSI frequencies.  After MOTION, the user often finds
that more frequencies are needed (TFU vs TFI checks, CRITFREQ).  Instead of repeating the
whole ANALYS run, ANALYS is run on the *additional* frequencies only; the old and new FILE8
are renamed FILE81 and FILE82 and COMBIN merges them by frequency number into a new FILE8
that MOTION and STRESS then read.

Rules (D-CMB-01, D-CMB-02):

* the DOF maps (``eq_node``, ``eq_dof``), the frequency step, the analysis type, the
  simultaneous / load case number ``case`` (D-ANL-06: FILE8001, FILE8002 ... of different
  incoherent or vibration load cases must not be mixed) and the seismic environment (control
  direction ``cm``, angle ``ang``) must be identical (error);
* frequencies are merged and sorted ascending; a frequency number present in both files is
  an error unless ``EDUOPT,COMBINDUP,PREFER82`` (then the FILE82 value is used);
* the output may hold up to 1,500 frequencies (more -> EDU-02 warning).

COMBIN has no input deck (D-CMB-02).  Because ``run_module`` passes no options to deck-less
modules, the duplicate policy is read from an optional text file ``COMBIN.opt`` in the model
directory holding the EDUOPT line (``EDUOPT,COMBINDUP,PREFER82`` or ``COMBINDUP PREFER82``).

Batch use (``python -m sassi.modules.combin < COMBIN.inp``): the shared three-line protocol
skips blank lines, so give an existing placeholder such as ``.`` as the (unused) deck line::

    modelname
    .
    modelname_COMBIN.out
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Tuple

import numpy as np

from ..io.container import Container, read_container, write_container
from ..io.files import validate
from .base import ModuleContext, ModuleError, batch_main

NAME = "COMBIN"
MAX_FREQ = 1500
OPTION_FILE = "COMBIN.opt"


def read_options(workdir: Path) -> dict:
    """Optional ``COMBIN.opt``: EDUOPT key/value lines (``EDUOPT,COMBINDUP,PREFER82``)."""
    opts = {"COMBINDUP": "ERROR"}
    p = Path(workdir) / OPTION_FILE
    if not p.exists():
        return opts
    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s[0] in "*#!":
            continue
        toks = [t for t in re.split(r"[,\s]+", s) if t]
        if toks and toks[0].upper() == "EDUOPT":
            toks = toks[1:]
        if len(toks) >= 2:
            opts[toks[0].upper()] = toks[1].upper()
    return opts


def _load(path: Path, name: str) -> Container:
    try:
        c = read_container(path, kind="FILE8")
    except (ValueError, FileNotFoundError) as exc:
        raise ModuleError(f"{name}: {exc}") from exc
    probs = validate(c)
    if probs:
        raise ModuleError(f"{name} does not satisfy the FILE8 schema: " + "; ".join(probs))
    if np.any(np.diff(np.asarray(c["fnum"])) <= 0):
        raise ModuleError(f"{name}: frequency numbers are not strictly increasing")
    return c


def _same_value(a, b) -> bool:
    """Equality of two meta values, numerically when both are numbers (``1`` == ``1.0``)."""
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return a == b


def combine(c81: Container, c82: Container, prefer82: bool = False) -> Tuple[dict, dict, dict]:
    """Merge two FILE8 containers (D-CMB-01).

    Returns ``(arrays, meta, info)`` for the merged FILE8; ``info`` reports the counts and the
    duplicate frequency numbers.  Raises ModuleError on incompatible files or duplicates
    (unless ``prefer82``).
    """
    for key in ("eq_node", "eq_dof"):
        if not np.array_equal(np.asarray(c81[key]), np.asarray(c82[key])):
            raise ModuleError(f"FILE81 and FILE82 have different DOF maps ({key}); they come from different models")
    m1, m2 = c81.meta, c82.meta
    df1, df2 = float(m1["df"]), float(m2["df"])
    if abs(df1 - df2) > 1e-9 * max(abs(df1), abs(df2)):
        raise ModuleError(f"FILE81 and FILE82 have different frequency steps ({df1:.9g} / {df2:.9g} Hz)")
    for key in ("type", "cm"):
        if int(m1[key]) != int(m2[key]):
            raise ModuleError(f"FILE81 and FILE82 differ in '{key}' ({m1[key]} / {m2[key]}): not the same analysis")
    if not _same_value(m1.get("case"), m2.get("case")):
        raise ModuleError(f"FILE81 and FILE82 belong to different simultaneous / load cases (case {m1.get('case')} / "
                          f"{m2.get('case')}, D-ANL-06): merging them would mix the TFs of different excitations")
    if abs(float(m1["ang"]) - float(m2["ang"])) > 1e-9:
        raise ModuleError(f"FILE81 and FILE82 differ in the transformation angle ({m1['ang']} / {m2['ang']})")
    f1 = np.asarray(c81["fnum"], dtype=np.int64)
    f2 = np.asarray(c82["fnum"], dtype=np.int64)
    dup = np.intersect1d(f1, f2)
    if dup.size and not prefer82:
        raise ModuleError("frequency numbers present in both FILE81 and FILE82: "
                          + ", ".join(str(int(v)) for v in dup[:30]) + (" ..." if dup.size > 30 else "")
                          + " (EDUOPT,COMBINDUP,PREFER82 keeps the FILE82 values)")
    keep1 = ~np.isin(f1, dup) if prefer82 else np.ones(len(f1), bool)
    fnum = np.concatenate([f1[keep1], f2])
    src = np.concatenate([np.zeros(int(keep1.sum()), np.int64), np.ones(len(f2), np.int64)])
    row = np.concatenate([np.where(keep1)[0], np.arange(len(f2))])
    order = np.argsort(fnum, kind="stable")
    fnum, src, row = fnum[order], src[order], row[order]
    H1 = np.asarray(c81["H"])
    H2 = np.asarray(c82["H"])
    H = np.empty((len(fnum), H1.shape[1]), dtype=np.result_type(H1, H2))
    freq = np.empty(len(fnum))
    a, b = src == 0, src == 1
    H[a], H[b] = H1[row[a]], H2[row[b]]
    freq[a] = np.asarray(c81["freq"], dtype=float)[row[a]]
    freq[b] = np.asarray(c82["freq"], dtype=float)[row[b]]
    arrays = {"fnum": fnum, "freq": freq, "eq_node": np.asarray(c81["eq_node"]),
              "eq_dof": np.asarray(c81["eq_dof"]), "H": H}
    meta = {k: v for k, v in m1.items() if k not in ("format", "kind", "version", "module")}
    meta["combined_from"] = ["FILE81", "FILE82"]
    info = {"n81": len(f1), "n82": len(f2), "n": len(fnum), "dup": dup.tolist(),
            "hash_mismatch": m1.get("model_hash") != m2.get("model_hash")}
    return arrays, meta, info


def run(ctx: ModuleContext) -> int:
    """COMBIN module entry point: FILE81 + FILE82 -> FILE8 (requirements section 4.7)."""
    lst = ctx.listing
    p81 = ctx.require("FILE81", "ANALYS (rename its FILE8 to FILE81)")
    p82 = ctx.require("FILE82", "ANALYS (rename its FILE8 to FILE82)")
    opts = read_options(ctx.workdir)
    prefer82 = opts.get("COMBINDUP", "ERROR") == "PREFER82"
    c81 = _load(p81, "FILE81")
    c82 = _load(p82, "FILE82")
    lst.section("COMBIN: merge FILE81 and FILE82 by frequency number")
    for nm, c in (("FILE81", c81), ("FILE82", c82)):
        fn = np.asarray(c["fnum"])
        lst.write(f"   {nm}: {len(fn)} frequencies, numbers {int(fn[0]) if len(fn) else 0}..{int(fn[-1]) if len(fn) else 0},"
                  f" {len(c['eq_node'])} equations, df {float(c.meta['df']):.9g} Hz")
    arrays, meta, info = combine(c81, c82, prefer82=prefer82)
    if info["hash_mismatch"]:
        lst.warning("FILE81 and FILE82 carry different model hashes: check that both come from the same model")
    if info["dup"]:
        lst.warning(f"EDUOPT,COMBINDUP,PREFER82: {len(info['dup'])} duplicate frequencies taken from FILE82: "
                    + ", ".join(str(v) for v in info["dup"][:30]))
    if info["n"] > MAX_FREQ:
        lst.warning(f"EDU-02: {info['n']} frequencies exceed the manual limit of {MAX_FREQ}")
    out = write_container(ctx.path("FILE8"), "FILE8", arrays, meta=meta, module=NAME)
    lst.write(f"   FILE8 written: {info['n']} frequencies ({info['n81']} + {info['n82']}"
              + (f" - {len(info['dup'])} duplicates" if info["dup"] else "") + ")")
    lst.write(f"   output file {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
