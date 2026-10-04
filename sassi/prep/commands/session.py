"""Session, file and model commands; frequency sets; model-level lists and output requests.

Requirements section 3.4.A (INP, WRITE, SAVE, RESUME, MDL, MDLNAME, TIT, STATUS, ACTM, CPMODEL,
DMODEL, MODELLIST, CD, MKDIR, MOPT, GRAVITY, GROUNDELEV), section 3.4.B (FREQ, LFREQ), and the
storage of the list and request commands of section 3.4.F (DAMP, TOPL, AMP, NOUT, EOUT, RDND),
whose semantics are fixed by section 3.3 and D-PAR-14:

* list append/delete -- a non-zero first value appends, a first value of 0 deletes the list
  (FREQ, DAMP, TOPL drop later zeros; AMP keeps them, see :func:`cmd_amp`);
* request append -- each command appends one request; ``NOUT,0`` / ``EOUT,0`` / ``RDND,0`` clear.

Decisions: D-MDL-01 (model numbers from 0), D-MDL-02 (SAVE/RESUME ``<model>.sdb`` container,
versioned), D-MDL-03 (SAVE/RESUME need MDL), D-PAR-12 (path resolution), D-PAR-17 (WRITE order).
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Optional

import numpy as np

from ...model import SSIModel, fmt_num
from ...model.entities import ElementRequest, NodalRequest, RelDispRequest
from ...model.options import MoptRecord
from ...model.ssimodel import HOUSE_GELEV, HOUSE_GRAVITY, SCHEMA_VERSION
from ..registry import CommandReported, command

SDB_KIND = "SDB"


# ======================================================================================
# Input / output of command files
# ======================================================================================
@command("INP")
def cmd_inp(c):
    """INP,<filename>: execute the commands of a .pre file (spec 07 section 9.2.18).

    Errors *inside* the file do not make INP fail (L12: processing continues).  A file that
    cannot be found or read does: ``run_file`` has reported it, and INP is then marked as failed
    (``execute()`` returns False, ``sassi -c`` exits with status 1) without a second message.
    """
    name = c.str(1)
    if not name:
        c.fail("file name required")
    summary = c.interp.run_file(name)
    if not summary.ok:
        raise CommandReported(f"INP: {name} could not be read")


@command("WRITE")
def cmd_write(c):
    """WRITE,[<file>],[<path>]: write the model as commands (default <model>.pre in the model path)."""
    from ..writer import write_pre
    m = c.model
    fname = c.str(1)
    if not fname:
        if not m.name:
            c.fail("no file name given and the model has no name (use MDL or WRITE,<file>)")
        fname = m.name + ".pre"
    dirname = c.str(2)
    if dirname:
        d = Path(dirname).expanduser()
        if not d.is_absolute():
            d = c.interp.cwd / d
        target = d / fname
    else:
        target = c.output_path(fname)
    if not target.parent.exists():
        c.fail(f"directory {target.parent} does not exist")
    opts = c.interp.write_options
    text, notes = write_pre(m, mdl=bool(opts.get("mdl")), afwr=bool(opts.get("afwr")), filename=target.name)
    try:
        target.write_text(text, encoding="utf-8")
    except OSError as exc:
        c.fail(f"cannot write {target}: {exc.strerror or exc}")
    for n in notes:
        c.warn(n)
    c.confirm(f"{text.count(chr(10))} lines written to {target}")


def _sdb_path(c) -> Path:
    m = c.model
    if not m.name or not m.path:
        c.fail("the model has no name and path (use MDL first, D-MDL-03)")
    return Path(m.path) / f"{m.name}.sdb"


@command("SAVE")
def cmd_save(c):
    """SAVE: save the active model to <path>/<model>.sdb (D-MDL-02)."""
    from ...io.container import write_container
    p = _sdb_path(c)
    m = c.model
    text = json.dumps(m.to_json(), sort_keys=True)
    try:
        write_container(p, SDB_KIND, {"model": np.asarray(text)},
                        meta={"sdb_version": SCHEMA_VERSION, "model_name": m.name, "model_hash": m.model_hash()},
                        module="UI")
    except OSError as exc:
        c.fail(f"cannot write {p}: {exc.strerror or exc}")
    c.confirm(f"model saved to {p}")


@command("RESUME")
def cmd_resume(c):
    """RESUME: reload the last SAVE of the active model (D-MDL-02; unknown major versions refused)."""
    from ...io.container import read_container
    p = _sdb_path(c)
    if not p.exists():
        c.fail(f"{p} not found (SAVE first)")
    try:
        cont = read_container(p, kind=SDB_KIND)
    except (OSError, ValueError) as exc:
        c.fail(f"cannot read {p}: {exc}")
    ver = str(cont.meta.get("sdb_version", "0"))
    if ver.split(".")[0] != SCHEMA_VERSION.split(".")[0]:
        c.fail(f"{p.name} has database version {ver}; this build reads {SCHEMA_VERSION.split('.')[0]}.x")
    try:
        new = SSIModel.from_json(json.loads(str(cont["model"])))
    except (ValueError, KeyError) as exc:
        c.fail(f"{p.name}: {exc}")
    new.name, new.path = c.model.name, c.model.path
    c.interp.models[c.interp.active_model] = new
    c.confirm(f"model resumed from {p}")


# ======================================================================================
# Model identity and the models in memory
# ======================================================================================
@command("MDL", max_args=2)
def cmd_mdl(c):
    """MDL,<Model>,<Path>: set the model name and directory; also the working directory."""
    name = c.str(1)
    if not name:
        c.fail("model name required")
    pth = c.str(2)
    m = c.model
    if pth:
        d = Path(pth).expanduser()
        if "\\" in pth and not d.exists():
            d = Path(pth.replace("\\", "/")).expanduser()
        if not d.is_absolute():
            if re.match(r"^[A-Za-z]:[\\/]", pth) and os.sep == "/":
                c.warn(f"Windows path {pth} is interpreted relative to the working directory")
            d = c.interp.cwd / d
    elif m.path:
        d = Path(m.path)
    else:
        d = c.interp.cwd
    d = d.resolve()
    if not d.exists():
        try:
            d.mkdir(parents=True)
        except OSError as exc:
            c.fail(f"cannot create directory {d}: {exc.strerror or exc}")
        c.info(f"MDL: directory {d} created")
    m.name, m.path = name, str(d)
    c.interp.cwd = d
    c.confirm(f"model {c.interp.active_model}: name {name}, path {d}")


@command("MDLNAME", max_args=1)
def cmd_mdlname(c):
    """MDLNAME,<name>: change the model name only (path and title unchanged)."""
    name = c.str(1)
    if not name:
        c.fail("model name required")
    c.model.name = name
    c.confirm(f"model name {name}")


@command("TIT", text_from=1)
def cmd_tit(c):
    """TIT,<title>: model title (the rest of the line, commas included)."""
    c.model.title = c.text(1)
    c.confirm(f"title: {c.model.title}")


@command("ACTM", max_args=1)
def cmd_actm(c):
    """ACTM,<Model>: activate a model (an empty one is created when absent)."""
    n = c.int(1, required=True, what="Model")
    if n < 0:
        c.fail("model numbers are >= 0")
    existed = n in c.interp.models
    c.interp.activate(n)
    c.confirm(f"model {n} active" + ("" if existed else " (new empty model)"))


@command("CPMODEL", max_args=1)
def cmd_cpmodel(c):
    """CPMODEL,<Mdl>: copy the active model (name and path included) to model <Mdl>."""
    n = c.int(1, required=True, what="Mdl")
    if n < 0:
        c.fail("model numbers are >= 0")
    if n == c.interp.active_model:
        c.fail("destination is the active model")
    if n in c.interp.models:
        c.warn(f"model {n} overwritten")
    cp = c.model.copy()
    c.interp.models[n] = cp
    if cp.name:
        c.warn(f"model {n} has the same name and path as model {c.interp.active_model}; AFWRITE from it "
               f"overwrites the files of the original (use MDL)")
    c.confirm(f"model {c.interp.active_model} copied to model {n}")


@command("DMODEL", max_args=1)
def cmd_dmodel(c):
    """DMODEL,<Mdl>: remove a model from memory (files untouched)."""
    n = c.int(1, required=True, what="Mdl")
    if n not in c.interp.models:
        c.fail(f"model {n} is not in memory")
    if n == c.interp.active_model:
        c.interp.models[n] = SSIModel()
        c.confirm(f"model {n} removed; the active model {n} is now empty")
    else:
        del c.interp.models[n]
        c.confirm(f"model {n} removed from memory")


@command("MODELLIST", max_args=0)
def cmd_modellist(c):
    """MODELLIST: list the models in memory."""
    c.info(f"{'model':>6}  {'nodes':>8} {'elements':>9}  name / path")
    for n in sorted(c.interp.models):
        m = c.interp.models[n]
        mark = "*" if n == c.interp.active_model else " "
        ident = f"{m.name} / {m.path}" if m.name else "(no name)"
        c.info(f"{mark}{n:>5}  {len(m.nodes):>8} {m.n_elements():>9}  {ident}")


# ======================================================================================
# Working directory
# ======================================================================================
@command("CD", max_args=1)
def cmd_cd(c):
    """CD,<dir>: change the working directory (it must exist), spec 10 section 5.2."""
    d = c.str(1)
    if not d:
        c.info(f"working directory: {c.interp.cwd}")
        return
    p = Path(d).expanduser()
    if not p.is_absolute():
        p = c.interp.cwd / p
    if not p.is_dir():
        c.fail(f"directory {p} does not exist")
    c.interp.cwd = p.resolve()
    c.confirm(f"working directory {c.interp.cwd}")


@command("MKDIR", max_args=1)
def cmd_mkdir(c):
    """MKDIR,<dir>: create a directory (intermediate ones too); warning if it exists."""
    d = c.str(1)
    if not d:
        c.fail("directory name required")
    p = Path(d).expanduser()
    if not p.is_absolute():
        p = c.interp.cwd / p
    if p.exists():
        c.warn(f"{p} already exists")
        return
    try:
        p.mkdir(parents=True)
    except OSError as exc:
        c.fail(f"directory {p} was not created: {exc.strerror or exc}")
    c.confirm(f"directory {p} created")


# ======================================================================================
# Global model options
# ======================================================================================
@command("MOPT", max_args=4)
def cmd_mopt(c):
    """MOPT,<incomp>,<matrix>,<mass>,<force>: model options (record setter; blank = default)."""
    for k, name in enumerate(("incomp", "matrix", "mass", "force"), start=1):
        v = c.int(k, what=name)
        if v is not None and v not in (0, 1):
            c.fail(f"<{name}> must be 0 or 1")
    rec = MoptRecord("MOPT", [c.raw(k) for k in range(1, 5)])
    c.model.options.set_record(rec)
    c.confirm("MOPT: incompatible modes {}, GENERAL mass {}, masses {}, forces {}".format(
        "suppressed" if rec.incomp else "included", "weight units" if rec.matrix else "mass units",
        "overwrite" if rec.mass else "add", "overwrite" if rec.force else "add"))


def _set_house_arg(c, k: int, value: float, label: str) -> None:
    rec = c.model.options.ensure_record("HOUSE")
    rec.set_arg(k, value)
    c.confirm(f"HOUSE <{label}> = {fmt_num(value)}")


@command("GRAVITY", max_args=1)
def cmd_gravity(c):
    """GRAVITY,<grav>: set HOUSE <gravity> only (other HOUSE options untouched)."""
    g = c.float(1, required=True, what="grav")
    if g <= 0:
        c.warn("gravity <= 0 (CHECK Error 1)")
    _set_house_arg(c, HOUSE_GRAVITY, g, "gravity")


@command("GROUNDELEV", max_args=1)
def cmd_groundelev(c):
    """GROUNDELEV,<elev>: set HOUSE <gelev> only."""
    _set_house_arg(c, HOUSE_GELEV, c.float(1, required=True, what="elev"), "gelev")


# ======================================================================================
# STATUS
# ======================================================================================
@command("STATUS", max_args=0)
def cmd_status(c):
    """STATUS: global model information (spec 07 section 9.2.37)."""
    from ...conventions import unit_system
    m = c.model
    k = m.counts()
    c.info(f"Model {c.interp.active_model}: name '{m.name or '-'}', path '{m.path or '-'}'")
    c.info(f"Title: {m.title or '-'}")
    c.info(f"Working directory: {c.interp.cwd}")
    c.info(f"Nodes {k['nodes']} (interaction {k['interaction']}), coordinate systems {k['csys']} "
           f"(active {m.csys_active})")
    c.info(f"Groups {k['groups']} (active {m.group_active if m.group_active is not None else '-'}), "
           f"elements {k['elements']}; MACT {m.mact}, RACT {m.ract}")
    c.info(f"Tables: M {k['materials']}, L {k['layers']}, R {k['sections']}, SC {k['springs']}, "
           f"MX {k['matrices']}")
    c.info(f"Loads: F {k['forces']}, MM {k['moments']}; masses: MT {k['tmass']}, MR {k['rmass']}")
    g = m.gravity
    c.info(f"Gravity {fmt_num(g)} ({'British' if unit_system(g) == 'BS' else 'SI'} units), "
           f"ground elevation {fmt_num(m.ground_elevation)}")
    mo = m.mopt
    c.info(f"MOPT: incomp {mo.incomp}, matrix {mo.matrix}, mass {mo.mass}, force {mo.force}")
    if m.freq_sets:
        sets = ", ".join(f"{s} ({len(v)})" for s, v in sorted(m.freq_sets.items()))
        c.info(f"Frequency sets (count): {sets}; df = {_df_text(m)}")
    else:
        c.info("Frequency sets: none")
    c.info(f"Lists: DAMP {len(m.damp)}, TOPL {len(m.topl)}, AMP {len(m.amp)}; requests: NOUT {len(m.nout)}, "
           f"EOUT {len(m.eout)}, RDND {len(m.rdnd)}")
    aopt = m.options.record("AOPT")
    c.info(f"AOPT: {','.join(aopt.to_tokens()) if aopt else '(default)'}")
    symm = m.options.entries("SYMM")
    if symm:
        for key, rec in symm:
            c.info(f"SYMM {key}: {','.join(rec.to_tokens())}")
    names = [n for n in m.options.names() if n not in ("MOPT", "AOPT", "SYMM")]
    c.info("Stored options: " + (", ".join(names) if names else "none"))


# ======================================================================================
# Frequency sets (3.4.B)
# ======================================================================================
def _df_text(m: SSIModel) -> str:
    """The resolved frequency step for listings, or why it is not available (SITE fstep/delt/nft)."""
    try:
        return f"{m.frequency_step():.9g} Hz"
    except ValueError as exc:
        return f"unavailable ({exc})"


@command("FREQ")
def cmd_freq(c):
    """FREQ,<ndx>,<f1>,...,<f10>: append frequency numbers to set ndx; f1 = 0 deletes the set."""
    ndx = c.int(1, required=True, what="ndx")
    f1 = c.int(2, required=True, what="f1")
    m = c.model
    if f1 == 0:
        existed = m.freq_sets.pop(ndx, None) is not None
        c.confirm(f"frequency set {ndx} deleted" if existed else f"frequency set {ndx} cleared (was not defined)")
        return
    vals = [v for v in c.ints(2, default=0) if v != 0]
    if any(v < 0 for v in vals):
        c.fail("frequency numbers must be positive integers")
    if c.nargs > 11:
        c.warn("more than 10 frequency numbers on one line (manual limit); all are used")
    lst = m.freq_sets.setdefault(ndx, [])
    dup = sorted(set(v for v in vals if v in lst) | set(v for v in vals if vals.count(v) > 1))
    lst.extend(vals)
    if dup:
        c.warn(f"duplicate frequency numbers in set {ndx}: {dup} (SITE stops on duplicates)")
    c.confirm(f"frequency set {ndx}: {len(lst)} numbers")


@command("LFREQ", max_args=3)
def cmd_lfreq(c):
    """LFREQ,[start],[end],[step]: list frequency sets with numbers and Hz values (current df)."""
    m = c.model
    if not m.freq_sets:
        c.info("no frequency sets defined")
        return
    start = c.int(1, default=min(1, min(m.freq_sets)))
    end = c.int(2, default=max(m.freq_sets))
    step = c.int(3, default=1)
    if step <= 0:
        step = 1
    try:
        df: Optional[float] = m.frequency_step()
    except ValueError as exc:
        df = None
        c.warn(f"frequency step unknown ({exc}): frequency numbers listed without Hz values")
    for s in sorted(m.freq_sets):
        if s < start or s > end or (s - start) % step:
            continue
        nums = sorted(m.freq_sets[s])
        if df is None:
            c.info(f"Frequency set {s}: {len(nums)} frequencies, df unavailable")
            cells = [f"{n:>6d}" for n in nums]
            per_row = 10
        else:
            c.info(f"Frequency set {s}: {len(nums)} frequencies, df = {df:.9g} Hz")
            cells = [f"{n:>6d} {n * df:>10.4f}" for n in nums]
            per_row = 6
        for k in range(0, len(cells), per_row):
            c.info("  " + "  ".join(cells[k:k + per_row]))


# ======================================================================================
# List commands of section 3.4.F (list append / delete)
# ======================================================================================
@command("DAMP")
def cmd_damp(c):
    """DAMP,<d1>,...,<d10>: append RS damping ratios (fractions); d1 = 0 clears the list."""
    d1 = c.float(1, required=True, what="d1")
    m = c.model
    if d1 == 0:
        m.damp.clear()
        c.confirm("damping list cleared")
        return
    vals = [v for v in c.floats(1) if v != 0]
    if c.nargs > 10:
        c.warn("more than 10 values on one line (manual limit); all are used")
    m.damp.extend(vals)
    c.confirm(f"damping list: {', '.join(fmt_num(v) for v in m.damp)}")


@command("TOPL")
def cmd_topl(c):
    """TOPL,<l1>,...: append L-layer numbers (top first) to the free-field list; l1 = 0 clears it."""
    l1 = c.int(1, required=True, what="l1")
    m = c.model
    if l1 == 0:
        m.topl.clear()
        c.confirm("top-layer list cleared")
        return
    vals = [v for v in c.ints(1) if v != 0]
    if any(v < 0 for v in vals):
        c.fail("layer numbers must be positive")
    m.topl.extend(vals)
    if len(m.topl) > 200:
        c.warn(f"{len(m.topl)} top layers (manual limit 200)")
    c.confirm(f"top-layer list: {len(m.topl)} layers")


@command("AMP")
def cmd_amp(c):
    """AMP,<no>,<a1>,...,<a100>: append spectral amplification ratios of motion <no>; a1 = 0 clears.

    Zero values after the first are **kept** (implementation decision of this package, proposed
    for the requirements decision log).  The generic list rule of spec 07 section 1.5 ("a non-zero first value
    appends non-zero values") would drop them, but for AMP a 0 is data:

    * D-INC-10: with HOUSE ``<cmplxspec>`` = 1 the ratios are consecutive (Re, Im) pairs, and a
      real-valued complex SAR has Im = 0; dropping it would shift every later pair;
    * spec 07 section 9.2.4: a SAR of 0 is inside the admissible range [0, 10] (Error 118), and
      the number of ratios must equal the number of frequencies of the set (Error 119), so a
      dropped value would change the count.

    Only the *first* value of a line keeps its list-class meaning (0 clears the list); WRITE never
    starts a continuation line with 0 (writer ``_w_amp``).  Blank fields are skipped as usual.
    """
    no = c.int(1, required=True, what="no")
    a1 = c.float(2, required=True, what="a1")
    m = c.model
    if a1 == 0:
        m.amp.pop(no, None)
        c.confirm(f"amplification ratios of motion {no} cleared")
        return
    vals = [c.float(k) for k in range(2, c.nargs + 1) if c.given(k)]
    m.amp.setdefault(no, []).extend(vals)
    c.confirm(f"motion {no}: {len(m.amp[no])} amplification ratios")


# ======================================================================================
# Output requests (request append, D-PAR-14)
# ======================================================================================
def _flags(c, k_from: int, n: int):
    return [c.int(k, default=0) for k in range(k_from, k_from + n)]


@command("NOUT")
def cmd_nout(c):
    """NOUT,<dir>,<code1>,...,<code6>,<node list>: nodal output request; NOUT,0 clears the list."""
    d = c.int(1, required=True, what="dir")
    m = c.model
    if d == 0:
        m.nout.clear()
        c.confirm("nodal output requests cleared")
        return
    if not 1 <= d <= 6:
        c.fail("<dir> must be 1..6 (x, y, z, xx, yy, zz)")
    codes = _flags(c, 2, 6)
    nodes = c.id_list(8)
    if not nodes:
        c.fail("node list is empty")
    m.nout.append(NodalRequest(d, codes, nodes))
    c.confirm(f"NOUT request {len(m.nout)}: direction {d}, {len(nodes)} nodes")


@command("EOUT")
def cmd_eout(c):
    """EOUT,<code1>,...,<code12>,<group>,<element list>: element output request; EOUT,0 clears."""
    m = c.model
    if c.nargs <= 1 and c.int(1, required=True, what="code1") == 0:
        m.eout.clear()
        c.confirm("element output requests cleared")
        return
    codes = _flags(c, 1, 12)
    if any(v not in (0, 1, 2) for v in codes):
        c.fail("output codes are 0 (none), 1 (maximum) or 2 (maximum and time history)")
    group = c.int(13, required=True, what="group")
    elems = c.id_list(14)
    if not elems:
        c.fail("element list is empty")
    m.eout.append(ElementRequest(codes, group, elems))
    c.confirm(f"EOUT request {len(m.eout)}: group {group}, {len(elems)} elements")


@command("RDND", max_args=7)
def cmd_rdnd(c):
    """RDND,<NodeNum>,<X>,<Y>,<Z>,<XX>,<YY>,<ZZ>: RELDISP node request (flag >= 1 = on); RDND,0 clears."""
    node = c.int(1, required=True, what="NodeNum")
    m = c.model
    if node == 0:
        m.rdnd.clear()
        c.confirm("RELDISP node requests cleared")
        return
    m.rdnd.append(RelDispRequest(node, _flags(c, 2, 6)))
    c.confirm(f"RDND request {len(m.rdnd)}: node {node}")
