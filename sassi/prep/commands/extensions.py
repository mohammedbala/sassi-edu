"""Extension commands of the SASSI-EDU dialect (requirements section 3.4.R; full names only, L11).

* X-commands (SITEX, SOILX, HOUSEX, ANALYSX, MOTIONX, STRESSX, RELDX) carry dialog fields that have
  no argument in the manual's commands.  They are record setters; a record equal to its defaults is
  **not stored**, so that WRITE (which writes an X-command only when it is non-default) gives a
  canonical round trip (UT-03) and default models produce manual-identical ``.pre`` files.
* CMODFORM selects the complex-modulus form (D-CNV-03); EDUOPT sets the algorithm switches of the
  decision log (requirements section 7).
* FCOPY / FMOVE copy or rename files in the model directory (FILE1 -> FILE1X, FILE8 -> FILE81 ...).
* VERIFY runs the verification problems of requirements section 6 (``sassi.verify``).
* REMOVEFREQ (P1) is ``Remove_Frequencies_from_FILE8``.
"""
from __future__ import annotations

import os
import shutil
import tempfile
import traceback
from pathlib import Path

import numpy as np

from ...model.values import NumberError, parse_float
from ..options import EDUOPT_KEYS
from ..registry import command
from .options import build_record, report_problems, store_entry


# ======================================================================================
# X-commands
# ======================================================================================
def store_xrecord(c, name: str):
    """Store an X-command record; a default record is removed instead (canonical state)."""
    rec = build_record(c, name, c.tokens)
    if rec.is_default():
        existed = c.model.options.record(name) is not None
        c.model.options.delete_record(name)
        c.confirm(f"{name}: default values" + (" (record removed)" if existed else ""))
        return rec
    c.model.options.set_record(rec)
    report_problems(c, rec)
    c.confirm(f"{name} options set")
    return rec


def _tier(c, what: str, tier: str) -> None:
    c.warn(f"{what} is not available in this build (tier {tier}): stored only")


@command("SITEX")
def cmd_sitex(c):
    """SITEX,<soilmode>: 0 Linear Soil (L table), 1 Non-Linear Soil (FILE88 from SOIL), D-SIT-06."""
    rec = store_xrecord(c, "SITEX")
    if rec.soilmode not in (0, 1):
        c.warn("<soilmode> must be 0 (linear) or 1 (non-linear, FILE88)")


@command("SOILX", text_from=5)
def cmd_soilx(c):
    """SOILX,<indir>,<mult>,<max>,[<cl>],[<file>]: SOIL input direction, SOIL scaling (exactly one of
    mult/max non-zero) and optional SOIL-only control layer and history file (D-SOL-11, D-SOL-13)."""
    rec = store_xrecord(c, "SOILX")
    if rec.indir not in (0, 1):
        c.warn("<indir> must be 0 (horizontal, Vs) or 1 (vertical, Vp)")
    if rec.indir == 1:
        c.info("SOILX: vertical input uses Vp and P-wave damping; set SOIL <iter> = 0 (no iterations)")


@command("HOUSEX")
def cmd_housex(c):
    """HOUSEX,<optimize>,<supmode>,<nsim>,<nlssi>,<ansys>: HOUSE dialog fields without an argument."""
    rec = store_xrecord(c, "HOUSEX")
    if not 1 <= rec.nsim <= 50:
        c.warn("<nsim> must be 1..50 (G-27)")
    # the node optimizer and non-linear soil SSI are implemented in HOUSE (wave 3)
    if rec.supmode or rec.nsim > 1:
        _tier(c, "HOUSEX quadratic superposition / stochastic simulations", "P1")
    if rec.ansys:
        _tier(c, "HOUSEX ANSYS model input (Option AA)", "P2")


@command("ANALYSX")
def cmd_analysx(c):
    """ANALYSX,<ffm>,<delrst>: Free-Field Load (0) / Free-Field Motion (1); delete restart files."""
    rec = store_xrecord(c, "ANALYSX")
    if rec.ffm:
        _tier(c, "ANALYSX free-field motion (FFM)", "P1")


@command("MOTIONX")
def cmd_motionx(c):
    """MOTIONX,<f1213>,<resp>,<srss>,<savetf>,<saveacc>,<savers>,<saverot>,<rsttf>,<rstacc>,<rstrs>."""
    rec = store_xrecord(c, "MOTIONX")
    if rec.f1213 not in (0, 1, 2):
        c.warn("<f1213> must be 0 (none), 1 (FILE13) or 2 (FILE12)")
    if rec.resp not in (0, 1, 2):
        c.warn("<resp> must be 0 (displacement), 1 (velocity) or 2 (acceleration)")


@command("STRESSX")
def cmd_stressx(c):
    """STRESSX,<pzadj>,<smo>,<skip>,<savemax>,<saveth>,<rstns>,<rstsp>: STRESS dialog fields."""
    rec = store_xrecord(c, "STRESSX")
    if rec.smo < 0:
        c.warn("<smo> must be >= 0")


@command("RELDX")
def cmd_reldx(c):
    """RELDX,<saverot>,<rstframes>: Save Rotations for ANSYS; Restart for Frame Generation (P1)."""
    store_xrecord(c, "RELDX")             # both fields are read by RELDISP (THDR files, THD/THDR frames)


@command("CMODFORM")
def cmd_cmodform(c):
    """CMODFORM,<form>: complex modulus 0 = 1-2b^2+2ib sqrt(1-b^2) (default), 1 = 1+2ib (benchmarks)."""
    rec = build_record(c, "CMODFORM", c.tokens)
    if rec.form not in (0, 1):
        c.fail("<form> must be 0 (SASSI/SHAKE91 form) or 1 (1 + 2 i beta)")
    c.model.options.set_record(rec)
    c.confirm("complex modulus: " + ("1 + 2 i beta (CMODFORM 1, benchmarks only)" if rec.form
                                     else "1 - 2 beta^2 + 2 i beta sqrt(1 - beta^2)"))


@command("EDUOPT")
def cmd_eduopt(c):
    """EDUOPT,<key>,<value>: algorithm switch of requirements section 7; a blank value resets the key."""
    key = c.word(1)
    if not key:
        c.fail("key required (" + ", ".join(sorted(EDUOPT_KEYS)) + ")")
    if not c.given(2):
        existed = c.model.options.delete_entry("EDUOPT", key)
        c.confirm(f"EDUOPT {key} reset to its default" if existed else f"EDUOPT {key} was not set")
        return
    val = c.str(2)
    spec = EDUOPT_KEYS.get(key)
    if spec is None:
        c.warn(f"unknown key {key} (stored; known: {', '.join(sorted(EDUOPT_KEYS))})")
    else:
        default, kind, dec = spec
        if isinstance(kind, tuple) and val.upper() not in kind:
            c.warn(f"{key} = {val}: expected one of {', '.join(kind)} ({dec})")
        elif kind is float:
            try:
                parse_float(val)
            except NumberError:
                c.fail(f"{key} needs a number ({dec})")
    store_entry(c, "EDUOPT")


# ======================================================================================
# Files of the model directory
# ======================================================================================
def model_dir(c) -> Path:
    """The model directory (MDL path) or, without MDL, the working directory."""
    m = c.model
    return Path(m.path) if m.path else c.interp.cwd


def _in_dir(c, name: str) -> Path:
    raw = os.path.expanduser(name.strip())
    if "\\" in raw and os.sep == "/":
        raw = raw.replace("\\", "/")
    p = Path(raw)
    return p if p.is_absolute() else model_dir(c) / p


def _copy_or_move(c, move: bool) -> None:
    src_n, dst_n = c.str(1), c.str(2)
    if not src_n or not dst_n:
        c.fail("source and destination file names are required")
    src, dst = _in_dir(c, src_n), _in_dir(c, dst_n)
    if not src.is_file():
        c.fail(f"{src} not found")
    if dst.is_dir():
        dst = dst / src.name
    if dst.resolve() == src.resolve():
        c.fail("source and destination are the same file")
    if dst.exists():
        c.warn(f"{dst.name} overwritten")
    try:
        if move:
            os.replace(src, dst)
        else:
            shutil.copy2(src, dst)
    except OSError as exc:
        c.fail(f"cannot {'move' if move else 'copy'} {src.name}: {exc.strerror or exc}")
    c.confirm(f"{src.name} {'renamed' if move else 'copied'} to {dst}")


@command("FCOPY", max_args=2)
def cmd_fcopy(c):
    """FCOPY,<src>,<dst>: copy a file of the model directory (e.g. FILE1 -> FILE1X, FILE8 -> FILE81)."""
    _copy_or_move(c, move=False)


@command("FMOVE", max_args=2)
def cmd_fmove(c):
    """FMOVE,<src>,<dst>: rename a file of the model directory (e.g. FILE8 -> FILE82)."""
    _copy_or_move(c, move=True)


@command("REMOVEFREQ")
def cmd_removefreq(c):
    """REMOVEFREQ,<infile>,<outfile>,<n1>,...: Remove_Frequencies_from_FILE8 (P1).

    The frequency numbers ``n1 ...`` are removed from the FILE8-type container ``<infile>`` and the
    result is written to ``<outfile>``; when both are the same file the original is kept as
    ``<infile>.bak``.
    """
    from ...io.container import read_container, write_container
    src_n, dst_n = c.str(1), c.str(2)
    if not src_n or not dst_n:
        c.fail("REMOVEFREQ,<infile>,<outfile>,<n1>,... (input and output file names required)")
    nums = set(c.id_list(3))
    if not nums:
        c.fail("no frequency numbers to remove")
    src, dst = _in_dir(c, src_n), _in_dir(c, dst_n)
    try:
        cont = read_container(src)
    except (OSError, ValueError) as exc:
        c.fail(f"cannot read {src}: {exc}")
    if "fnum" not in cont.arrays or "H" not in cont.arrays:
        c.fail(f"{src.name} is not a FILE8-type file (no fnum/H)")
    fnum = np.asarray(cont["fnum"])
    keep = ~np.isin(fnum, list(nums))
    missing = sorted(nums - set(int(v) for v in fnum))
    if missing:
        c.warn(f"frequency numbers not in {src.name}: {missing[:10]}")
    nf = len(fnum)
    arrays = {}
    for k, a in cont.arrays.items():
        arrays[k] = a[keep] if (k in ("fnum", "freq", "H") and a.shape[:1] == (nf,)) else a
    if dst.resolve() == src.resolve():
        bak = src.with_name(src.name + ".bak")
        shutil.copy2(src, bak)
        c.info(f"REMOVEFREQ: backup {bak.name}")
    meta = {k: v for k, v in cont.meta.items() if k not in ("format", "kind", "version", "module")}
    write_container(dst, cont.kind, arrays, meta=meta, module=str(cont.meta.get("module", "")))
    c.confirm(f"{int((~keep).sum())} frequencies removed: {dst.name} has {int(keep.sum())} frequencies")


# ======================================================================================
# VERIFY
# ======================================================================================
def _load_problems(c):
    """``sassi.verify.load_all`` with a fallback that reports problem modules failing to import."""
    import importlib
    import pkgutil
    from ... import verify
    try:
        return verify.load_all()
    except Exception:            # another package's VP module may be incomplete
        from ...verify import problems
        for mm in pkgutil.iter_modules(problems.__path__):
            try:
                importlib.import_module(f"{problems.__name__}.{mm.name}")
            except Exception as exc:
                c.warn(f"verification module {mm.name} failed to load: {exc!r}")
        return verify.REGISTRY


def _run_registered(p, workdir: Path):
    """Run one registered verification problem in ``workdir``.

    This is :func:`sassi.verify.run_problem` without its ``load_all()`` call: the problems are
    already registered by :func:`_load_problems`, and a verification module of another package
    that fails to import must not stop the problems that did register (concurrent development).
    """
    import time
    workdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    res = p.func(workdir)
    res.elapsed = time.time() - t0
    return res


@command("VERIFY", max_args=1)
def cmd_verify(c):
    """VERIFY,[id|P0|P1|P2|ALL|LIST]: run verification problems and print computed, reference, error,
    tolerance and pass/fail (requirements section 6.1); blank = P0."""
    reg = _load_problems(c)
    sel = c.str(1, "P0").strip().upper()
    if sel == "LIST":
        for pid in sorted(reg):
            p = reg[pid]
            c.info(f"{pid:8s} {p.tier}  {p.title}  [{', '.join(p.modules)}]")
        return
    if sel in ("ALL",):
        ids = sorted(reg)
    elif sel in ("P0", "P1", "P2"):
        ids = sorted(k for k, p in reg.items() if p.tier == sel)
    else:
        pid = sel if sel.startswith("VP-") else f"VP-{sel.zfill(2)}" if sel.isdigit() else sel
        match = [k for k in reg if k.upper() == pid]
        if not match:
            c.fail(f"unknown verification problem {c.str(1)} (VERIFY,LIST lists them)")
        ids = match
    if not ids:
        c.info(f"VERIFY: no verification problem of tier {sel}")
        return
    base = Path(tempfile.mkdtemp(prefix="sassi_verify_"))
    n_pass = n_fail = 0
    for pid in ids:
        p = reg[pid]
        c.info(f"{pid}: {p.title}")
        try:
            res = _run_registered(p, base / pid)
        except Exception as exc:
            n_fail += 1
            c.error(f"{pid} failed to run: {exc!r}")
            c.info(traceback.format_exc().strip().splitlines()[-1])
            continue
        for ch in res.checks:
            kind = "rel" if ch.kind == "rel" else ("abs" if ch.kind == "abs" else "")
            status = "pass" if ch.passed else "FAIL"
            if ch.kind == "bool":
                c.info(f"   {status}  {ch.quantity}")
            elif ch.kind == "info":
                c.info(f"   info  {ch.quantity}: computed {ch.computed:.6g}, published {ch.reference:.6g} "
                       f"(difference {ch.error:.3g}, informative only)")
            else:
                c.info(f"   {status}  {ch.quantity}: computed {ch.computed:.6g}, reference {ch.reference:.6g}, "
                       f"error {ch.error:.3g} ({kind}), tolerance {ch.tolerance:.3g}")
        for n in res.notes:
            c.info(f"   note: {n}")
        if res.passed:
            n_pass += 1
        else:
            n_fail += 1
        c.info(f"   {pid} {'PASSED' if res.passed else 'FAILED'} ({res.elapsed:.2f} s)")
    if n_fail:
        # keep the work files of a failed run for inspection; a clean run leaves nothing behind
        msg = f"VERIFY {sel}: {n_pass} passed, {n_fail} failed (work files kept in {base})"
        c.warn(msg)
    else:
        shutil.rmtree(base, ignore_errors=True)
        msg = f"VERIFY {sel}: {n_pass} passed, {n_fail} failed"
        c.confirm(msg)
        c.info(msg)
