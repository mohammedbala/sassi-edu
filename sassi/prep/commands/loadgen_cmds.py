"""Option A commands: the LOADGEN module (ACS SASSI-ANSYS two-step approach, manual 6.4.15, spec 04
section 15.6, requirements 1.1 and 5.3; :mod:`sassi.modules.loadgen`).

The manual runs LOADGEN from two Modules-menu dialogs only ("ANSYS Eq. Static Load", "ANSYS Dynamic
Load"): on Ok the UI writes the LOADGEN input file and starts LOADGEN.  SASSI-EDU keeps the dialog
fields as option records, so that a session can be replayed from a ``.pre`` file (L17) and WRITE -> INP
keeps them (UT-03); ``RUNLOADGEN`` is the dialog's Ok button (D-LGN-01).  All the records are written by
WRITE under "Other stored options".

Dialog fields
-------------
``LOADGEN,<data>,<multi>,<rotdisp>,<rotacc>,<masstype>,<genmass>,<source>`` -- "ANSYS Static Load
Converter":

* ``<data>`` Data to Add From ACS SASSI: 1 (or DISP) Displacement, 2 (ACC, default) Acceleration,
  3 (DISPACC) Disp. and Accel., 4 (SOILDISP) Disp. for Soil Module;
* ``<multi>`` Use Multiple File List Inputs: 1 = one APDL file per critical time (0: one file);
* ``<rotdisp>``, ``<rotacc>`` the Rotational Disp. / Rotational Accel. check boxes (0/1);
* ``<masstype>`` Mass Type: 1 (LUMPED, default) Lumped Mass, 2 (MASTER) Master Node Mass;
* ``<genmass>`` Generate Mass Data (0/1): the mass file is written (overwritten) from the HOUSE mass matrix;
* ``<source>`` (SASSI-EDU) RESULTS (default: MOTION ``.ACC`` and RELDISP ``.THD`` files, or their frames)
  or FILE8 (the same histories computed from FILE8 and the control motion of the MOTION deck).

``LOADGENDYN,<alpha>,<beta>,<method>,<refnode>,<source>,<rotdisp>,<rotacc>,<zeta>,<f1>,<f2>,<gfopt>,<gmult>``
-- "ANSYS Dynamic Load Converter": Rayleigh damping ``<alpha>``, ``<beta>`` (C = alpha M + beta K); and
the SASSI-EDU fields ``<method>`` REL (default, the manual's: ACEL = ground motion + D = interface
displacements relative to it) or ACC (fixed interface, ACEL = absolute acceleration of ``<refnode>``),
``<refnode>`` the node whose motion is the reference (0 = the control motion, the free field),
``<source>``, ``<rotdisp>``, ``<rotacc>`` as above, ``<zeta>``, ``<f1>``, ``<f2>``: alpha and beta from a
damping ratio zeta at two frequencies (``alpha = 2 zeta w1 w2/(w1 + w2)``, ``beta = 2 zeta/(w1 + w2)``,
used when zeta > 0), ``<gfopt>``, ``<gmult>`` format (0 dt then values, 1 (t, a) pairs) and scale factor
of the Ground Acceleration File.

``LGFILE,<key>,<name>`` -- the file and path boxes (the name is the rest of the line; blank = default):

====================  =============================================================  =======================
key                   dialog box                                                     default
====================  =============================================================  =======================
SSIPATH               SASSI Model and Results Input: Path                            model directory
HOUSE                 HOUSE Module Input                                             ``<model>.hou``
DISP, DISPROT         Displacement Results (translational / rotational frames)       ``THD``, ``THDR``
ACC, ACCROT           Acceleration Results (translational / rotational frames)       ``ACC``, ``ACCR``
ANSYSPATH             ANSYS Model and Data Input: Path                               model directory
LUMPED, MASTER        Lumped node / Master Node Mass file                            ``<model>.masl/.masm``
APDL, APDLDYN         ANSYS Output File: APDL file (static / dynamic)                ``<model>_LGS.inp`` /
                                                                                     ``<model>_LGD.inp``
GROUND                Ground Acceleration File (dynamic, in g)                       the control motion
====================  =============================================================  =======================

``LGNODE,<kind>,<nodes>`` -- node lists (model numbering; ``a-b`` ranges; each command adds,
``LGNODE,<kind>,0`` clears): D = nodes that receive D (default: the interaction nodes; "Disp. for Soil
Module": also the nodes of the excavated soil), M = master nodes (Master Node Mass generation), A = nodes
whose SSI absolute accelerations are exported as tables for comparison (dynamic).

``LGTIME,<crit>,...`` -- critical times of the static loads (D-LGN-07):
``LGTIME,V|VX|VY|VZ,<n>,<tsep>`` the n largest peaks (at least tsep s apart) of the base shear of the
inertia forces (V: along the control direction); ``LGTIME,MX|MY|MZ,<n>,<tsep>,<x0>,<y0>,<z0>`` of the
overturning moment about (x0, y0, z0) (blank: centroid of the interface nodes);
``LGTIME,ACC|DISP,<n>,<tsep>,<node>,<dof>`` of a node's absolute acceleration / relative displacement;
``LGTIME,TIME,<t1>,<t2>,...`` given times (s); ``LGTIME,STEP,<k1>,...`` given time steps (1-based, the
frame numbers).  Default: ``LGTIME,V,1``.

``LGMAP,<mode>,<tol>,<file>`` -- node numbering of the ANSYS model (D-LGN-05): IDENTITY (default; the
``ANSYS`` command and ``CONVERT,ANSYS`` keep the node numbers), PAIRS (``<file>`` lines ``sassi_node
ansys_node``) or COORD (the ANSYS node at the same position in the ``.cdb`` / APDL ``<file>``, within
``<tol>``, default 1e-6 of the model size).  The HOUSE optimizer numbering is translated automatically.

``LGOPT,<digits>,<opmode>,<rest>`` -- significant digits of the APDL values (default 12), the data-check
mode (1 = no file written) and <rest> (1, default: the relative displacements start at rest, D-LGN-06).

``LGLIST`` lists the settings; ``RUNLOADGEN,[STATIC|DYNAMIC],[model]`` (or, in the RUN<MODULE> order,
``RUNLOADGEN,<model>,[STATIC|DYNAMIC]``) writes ``<model>.lgn`` and runs LOADGEN in the model directory
(STATIC by default), like the dialog's Ok.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional

from ...model.options import Field, Record, make_record, register_record_type
from ...model.values import NumberError, parse_float, parse_int
from ..lexer import compress_ids
from ..registry import CommandReported, command
from .. import writer

TIER = "P2"

DATA_WORDS = {"1": 1, "DISP": 1, "DISPLACEMENT": 1, "2": 2, "ACC": 2, "ACCEL": 2, "ACCELERATION": 2,
              "3": 3, "DISPACC": 3, "BOTH": 3, "4": 4, "SOILDISP": 4, "SOIL": 4}
MASS_WORDS = {"1": 1, "LUMPED": 1, "LUMP": 1, "2": 2, "MASTER": 2}
SOURCES = ("RESULTS", "FILE8")
METHODS = ("REL", "ACC")
FILE_KEYS = {"SSIPATH": "ssipath", "HOUSE": "housefile", "DISP": "dispdir", "DISPROT": "disprotdir",
             "ACC": "accdir", "ACCROT": "accrotdir", "ANSYSPATH": "ansyspath", "LUMPED": "lumpfile",
             "MASTER": "masterfile", "APDL": "apdlfile", "APDLDYN": "apdlfile", "GROUND": "groundfile"}
NODE_KINDS = {"D": "dnodes", "M": "masters", "A": "checknodes"}
NODE_KIND_NAMES = {"D": "interface nodes that receive D", "M": "master nodes", "A": "acceleration check nodes"}
CRITERIA = ("V", "VX", "VY", "VZ", "MX", "MY", "MZ", "ACC", "DISP", "TIME", "STEP")
MAPMODES = ("IDENTITY", "PAIRS", "COORD")


# ======================================================================================
# Records
# ======================================================================================
@register_record_type
class LoadgenRecord(Record):
    """``LOADGEN,<data>,<multi>,<rotdisp>,<rotacc>,<masstype>,<genmass>,<source>``: ANSYS Static Load Converter."""
    COMMAND = "LOADGEN"
    FIELDS = (Field("data", str, "2", "1 Displacement, 2 Acceleration, 3 Disp. and Accel., 4 Disp. for Soil Module"),
              Field("multi", int, 0, "Use Multiple File List Inputs"),
              Field("rotdisp", int, 0, "Rotational Disp."),
              Field("rotacc", int, 0, "Rotational Accel."),
              Field("masstype", str, "1", "1 Lumped Mass, 2 Master Node Mass"),
              Field("genmass", int, 0, "Generate Mass Data"),
              Field("source", str, "RESULTS", "RESULTS or FILE8"))


@register_record_type
class LoadgenDynRecord(Record):
    """``LOADGENDYN,<alpha>,<beta>,<method>,<refnode>,<source>,<rotdisp>,<rotacc>,<zeta>,<f1>,<f2>,<gfopt>,<gmult>``:
    ANSYS Dynamic Load Converter."""
    COMMAND = "LOADGENDYN"
    FIELDS = (Field("alpha", float, 0.0, "Rayleigh alpha"), Field("beta", float, 0.0, "Rayleigh beta"),
              Field("method", str, "REL", "REL or ACC"), Field("refnode", int, 0, "reference node (0 = free field)"),
              Field("source", str, "RESULTS", "RESULTS or FILE8"), Field("rotdisp", int, 0, "rotations prescribed"),
              Field("rotacc", int, 0, "rotational accelerations in the check tables"),
              Field("zeta", float, 0.0, "damping ratio for alpha/beta"), Field("f1", float, 0.0, "frequency 1, Hz"),
              Field("f2", float, 0.0, "frequency 2, Hz"), Field("gfopt", int, 0, "ground file format"),
              Field("gmult", float, 1.0, "ground file scale factor"))


@register_record_type
class LgmapRecord(Record):
    """``LGMAP,<mode>,<tol>,<file>``: SASSI-EDU -> ANSYS node numbering."""
    COMMAND = "LGMAP"
    FIELDS = (Field("mode", str, "IDENTITY", "IDENTITY, PAIRS or COORD"), Field("tol", float, 0.0, "COORD tolerance"),
              Field("file", str, "", "pairs file or ANSYS .cdb/APDL file"))


@register_record_type
class LgoptRecord(Record):
    """``LGOPT,<digits>,<opmode>,<rest>``: APDL digits, data-check mode, relative displacements started at rest."""
    COMMAND = "LGOPT"
    FIELDS = (Field("digits", int, 12, "significant digits (6..17)"), Field("opmode", int, 0, "1 data check"),
              Field("rest", int, 1, "1 subtract the initial value of the relative displacements"))


# ======================================================================================
# helpers
# ======================================================================================
def _num(c, k: int, what: str, integer: bool = False, default=None):
    t = c.raw(k)
    if not t:
        return default
    try:
        if integer:
            v, exact = parse_int(t)
            if not exact:
                c.fail(f"argument {k} <{what}>: {t} must be an integer")
            return v
        return parse_float(t)
    except NumberError:
        c.fail(f"argument {k} <{what}>: '{t}' is not a number")


def _flag(c, k: int, what: str) -> int:
    v = _num(c, k, what, integer=True, default=0)
    if v not in (0, 1):
        c.fail(f"argument {k} <{what}> must be 0 or 1")
    return v


def _tokens(c, n: int) -> List[str]:
    return [c.raw(k) for k in range(1, n + 1)]


def rayleigh(zeta: float, f1: float, f2: float):
    """Rayleigh coefficients giving the damping ratio ``zeta`` at ``f1`` and ``f2`` (Hz):
    ``alpha = 2 zeta w1 w2/(w1 + w2)``, ``beta = 2 zeta/(w1 + w2)`` (zeta(w) = alpha/(2 w) + beta w/2)."""
    w1, w2 = 2.0 * math.pi * f1, 2.0 * math.pi * f2
    return 2.0 * zeta * w1 * w2 / (w1 + w2), 2.0 * zeta / (w1 + w2)


def _rec(m, name: str) -> Optional[Record]:
    return m.options.record(name)


def _get(rec: Optional[Record], k: int, default=None):
    return default if rec is None or rec.arg(k) is None else rec.arg(k)


# ======================================================================================
# LOADGEN, LOADGENDYN
# ======================================================================================
@command("LOADGEN", tier=TIER, max_args=7)
def cmd_loadgen(c):
    """LOADGEN,<data>,<multi>,<rotdisp>,<rotacc>,<masstype>,<genmass>,<source>: ANSYS Eq. Static Load options
    (Option A): data 1 DISP / 2 ACC / 3 DISPACC / 4 SOILDISP, Use Multiple File List Inputs, Rotational Disp.,
    Rotational Accel., Mass Type 1 LUMPED / 2 MASTER, Generate Mass Data, history source RESULTS / FILE8."""
    toks = _tokens(c, 7)
    data = (toks[0] or "2").upper()
    if data not in DATA_WORDS:
        c.fail(f"<data> must be 1 DISP, 2 ACC, 3 DISPACC or 4 SOILDISP (got '{toks[0]}')")
    for k, what in ((2, "multi"), (3, "rotdisp"), (4, "rotacc"), (6, "genmass")):
        _flag(c, k, what)
    mt = (toks[4] or "1").upper()
    if mt not in MASS_WORDS:
        c.fail(f"<masstype> must be 1 LUMPED or 2 MASTER (got '{toks[4]}')")
    src = (toks[6] or "RESULTS").upper()
    if src not in SOURCES:
        c.fail(f"<source> must be RESULTS or FILE8 (got '{toks[6]}')")
    toks[0] = toks[0].upper()
    toks[4] = toks[4].upper()
    toks[6] = toks[6].upper()
    c.model.options.set_record(LoadgenRecord("LOADGEN", toks))
    from ...modules.loadgen import DATA_NAMES, MASS_NAMES
    c.confirm(f"LOADGEN static options set: {DATA_NAMES[DATA_WORDS[data]]}, {MASS_NAMES[MASS_WORDS[mt]]}, source {src}")


@command("LOADGENDYN", tier=TIER, max_args=12)
def cmd_loadgendyn(c):
    """LOADGENDYN,<alpha>,<beta>,<method>,<refnode>,<source>,<rotdisp>,<rotacc>,<zeta>,<f1>,<f2>,<gfopt>,<gmult>:
    ANSYS Dynamic Load options (Option A): Rayleigh alpha/beta (or from zeta at f1, f2), method REL (ground
    ACEL + relative D, manual) / ACC (fixed base driven by <refnode>), reference node (0 = control motion),
    source RESULTS / FILE8, rotations, Ground Acceleration File format and factor."""
    toks = _tokens(c, 12)
    alpha = _num(c, 1, "alpha", default=0.0)
    beta = _num(c, 2, "beta", default=0.0)
    if alpha < 0 or beta < 0:
        c.fail("Rayleigh coefficients <alpha> and <beta> must be >= 0")
    method = (toks[2] or "REL").upper()
    if method not in METHODS:
        c.fail(f"<method> must be REL or ACC (got '{toks[2]}')")
    ref = _num(c, 4, "refnode", integer=True, default=0)
    if ref < 0:
        c.fail("<refnode> must be >= 0")
    if method == "ACC" and ref == 0:
        c.fail("method ACC needs the reference node <refnode> whose absolute acceleration drives the fixed base")
    src = (toks[4] or "RESULTS").upper()
    if src not in SOURCES:
        c.fail(f"<source> must be RESULTS or FILE8 (got '{toks[4]}')")
    _flag(c, 6, "rotdisp")
    _flag(c, 7, "rotacc")
    zeta = _num(c, 8, "zeta", default=0.0)
    f1 = _num(c, 9, "f1", default=0.0)
    f2 = _num(c, 10, "f2", default=0.0)
    if zeta:
        if not 0 < zeta < 1:
            c.fail("<zeta> must be a damping ratio in (0, 1)")
        if not 0 < f1 < f2:
            c.fail("Rayleigh damping from <zeta> needs 0 < <f1> < <f2> (Hz)")
        if c.given(1) or c.given(2):
            c.warn("<alpha>/<beta> ignored: they are computed from <zeta>, <f1>, <f2>")
    gfopt = _num(c, 11, "gfopt", integer=True, default=0)
    if gfopt not in (0, 1):
        c.fail("<gfopt> must be 0 (dt then values) or 1 ((t, a) pairs)")
    _num(c, 12, "gmult", default=1.0)
    toks[2] = toks[2].upper()
    toks[4] = toks[4].upper()
    c.model.options.set_record(LoadgenDynRecord("LOADGENDYN", toks))
    if zeta:
        a, b = rayleigh(zeta, f1, f2)
        c.info(f"LOADGENDYN: Rayleigh alpha = {a:.6g}, beta = {b:.6g} (zeta = {zeta:g} at {f1:g} and {f2:g} Hz)")
        alpha, beta = a, b
    if alpha == 0 and beta == 0:
        c.warn("no Rayleigh damping (alpha = beta = 0)")
    c.confirm(f"LOADGEN dynamic options set: method {method}, alpha {alpha:.6g}, beta {beta:.6g}, source {src}")


# ======================================================================================
# LGFILE, LGNODE, LGTIME, LGMAP, LGOPT
# ======================================================================================
@command("LGFILE", tier=TIER, text_from=2)
def cmd_lgfile(c):
    """LGFILE,<key>,<name>: a file or path box of the LOADGEN dialogs (keys SSIPATH HOUSE DISP DISPROT ACC ACCROT
    ANSYSPATH LUMPED MASTER APDL APDLDYN GROUND); <name> is the rest of the line, blank = default."""
    key = c.word(1)
    if key not in FILE_KEYS:
        c.fail(f"<key> must be one of {' '.join(FILE_KEYS)} (got '{c.raw(1)}')")
    name = c.text(2).strip()
    if name:
        c.model.options.set_entry("LGFILE", key, make_record("LGFILE", [key, name]))
        c.confirm(f"LOADGEN file {key} = {name}")
    else:
        existed = c.model.options.delete_entry("LGFILE", key)
        c.confirm(f"LOADGEN file {key}: default" + (" (entry removed)" if existed else ""))


def lgnode_lists(m) -> Dict[str, List[int]]:
    """{kind: sorted node list} of the LGNODE records."""
    out: Dict[str, List[int]] = {}
    from ..lexer import parse_id_list
    for key, rec in m.options.entries("LGNODE"):
        toks = [t for t in rec.to_tokens()[1:] if t]
        out[str(key)] = sorted(set(parse_id_list(toks))) if toks else []
    return out


@command("LGNODE", tier=TIER)
def cmd_lgnode(c):
    """LGNODE,<kind>,<nodes>: add nodes to a LOADGEN node list (kind D interface nodes that receive D, M master
    nodes, A acceleration check nodes); LGNODE,<kind>,0 clears the list.  Model node numbers, ranges a-b."""
    kind = c.word(1)
    if kind not in NODE_KINDS:
        c.fail(f"<kind> must be D, M or A (got '{c.raw(1)}')")
    if c.nargs < 2:
        c.fail("node list missing (LGNODE,<kind>,0 clears the list)")
    if c.nargs == 2 and c.raw(2) in ("0", "-0"):
        existed = c.model.options.delete_entry("LGNODE", kind)
        c.confirm(f"LOADGEN {NODE_KIND_NAMES[kind]}: list cleared" + ("" if existed else " (was empty)"))
        return
    nodes = c.id_list(2)
    if any(n <= 0 for n in nodes):
        c.fail("node numbers must be >= 1 (LGNODE,<kind>,0 alone clears the list)")
    known = set(c.model.nodes)
    unknown = [n for n in nodes if n not in known]
    if unknown and known:
        c.warn(f"nodes not in the model: {', '.join(map(str, unknown[:10]))}" + (" ..." if len(unknown) > 10 else ""))
    cur = set(lgnode_lists(c.model).get(kind, []))
    allnodes = sorted(cur | set(nodes))
    c.model.options.set_entry("LGNODE", kind, make_record("LGNODE", [kind] + compress_ids(allnodes)))
    c.confirm(f"LOADGEN {NODE_KIND_NAMES[kind]}: {len(allnodes)} nodes")


def _w_lgnode(m) -> List[str]:
    """WRITE of the LGNODE lists: a clear, then the nodes in lines of at most 20 tokens (LGNODE adds)."""
    out: List[str] = []
    for kind, nodes in sorted(lgnode_lists(m).items()):
        if not nodes:
            continue
        toks = compress_ids(nodes)
        out.append(f"LGNODE,{kind},0")
        for i in range(0, len(toks), 20):
            out.append(",".join(["LGNODE", kind] + toks[i:i + 20]))
    return out


writer.register_writer("LGNODE", _w_lgnode)


@command("LGTIME", tier=TIER)
def cmd_lgtime(c):
    """LGTIME,<crit>,...: critical times of the static loads -- V|VX|VY|VZ,<n>,<tsep>; MX|MY|MZ,<n>,<tsep>,<x0>,
    <y0>,<z0>; ACC|DISP,<n>,<tsep>,<node>,<dof>; TIME,<t1>,<t2>,...; STEP,<k1>,<k2>,... (default V,1)."""
    crit = c.word(1)
    if crit not in CRITERIA:
        c.fail(f"<crit> must be one of {' '.join(CRITERIA)} (got '{c.raw(1)}')")
    toks = [c.raw(k) for k in range(1, c.nargs + 1)]
    toks[0] = crit
    if crit in ("TIME", "STEP"):
        if c.nargs < 2:
            c.fail(f"LGTIME,{crit} needs at least one {'time' if crit == 'TIME' else 'time step'}")
        for k in range(2, c.nargs + 1):
            v = _num(c, k, "time" if crit == "TIME" else "step", integer=(crit == "STEP"))
            if v is None:
                c.fail(f"argument {k} is blank")
            if (crit == "TIME" and v < 0) or (crit == "STEP" and v < 1):
                c.fail(f"argument {k}: {'times must be >= 0' if crit == 'TIME' else 'steps are >= 1'}")
    else:
        n = _num(c, 2, "n", integer=True, default=1)
        tsep = _num(c, 3, "tsep", default=0.0)
        if n < 1:
            c.fail("<n> must be >= 1")
        if tsep < 0:
            c.fail("<tsep> must be >= 0")
        if crit in ("ACC", "DISP"):
            node = _num(c, 4, "node", integer=True, default=0)
            dof = _num(c, 5, "dof", integer=True, default=0)
            if node <= 0 or not 1 <= dof <= 6:
                c.fail(f"LGTIME,{crit} needs <node> >= 1 and <dof> 1..6")
            if c.nargs > 5:
                c.fail(f"LGTIME,{crit} takes 5 arguments")
        elif crit.startswith("M"):
            given = [c.given(k) for k in (4, 5, 6)]
            if any(given) and not all(given):
                c.fail("give all three coordinates <x0>,<y0>,<z0> of the moment point (or none)")
            for k in (4, 5, 6):
                _num(c, k, "x0")
            if c.nargs > 6:
                c.fail(f"LGTIME,{crit} takes 6 arguments")
        elif c.nargs > 3:
            c.fail(f"LGTIME,{crit} takes 3 arguments")
    c.model.options.set_record(make_record("LGTIME", toks))
    c.confirm(f"LOADGEN critical times: {','.join(t for t in toks if t)}")


@command("LGMAP", tier=TIER, text_from=3)
def cmd_lgmap(c):
    """LGMAP,<mode>,<tol>,<file>: ANSYS node numbering -- IDENTITY (default), PAIRS (file of 'sassi_node
    ansys_node' lines) or COORD (ANSYS node at the same position in a .cdb / APDL file, tolerance <tol>)."""
    mode = c.word(1, "IDENTITY")
    if mode not in MAPMODES:
        c.fail(f"<mode> must be IDENTITY, PAIRS or COORD (got '{c.raw(1)}')")
    tol = _num(c, 2, "tol", default=0.0)
    if tol < 0:
        c.fail("<tol> must be >= 0")
    name = c.text(3).strip()
    if mode != "IDENTITY" and not name:
        c.fail(f"LGMAP,{mode} needs the file name (argument 3)")
    if mode == "IDENTITY" and not c.given(2) and not name:
        existed = c.model.options.record("LGMAP") is not None
        c.model.options.delete_record("LGMAP")
        c.confirm("LOADGEN node map: identity" + (" (record removed)" if existed else ""))
        return
    c.model.options.set_record(LgmapRecord("LGMAP", [mode, c.raw(2), name]))
    c.confirm(f"LOADGEN node map: {mode}" + (f", {name}" if name else ""))


@command("LGOPT", tier=TIER, max_args=3)
def cmd_lgopt(c):
    """LGOPT,<digits>,<opmode>,<rest>: significant digits of the APDL values (6..17, default 12), data-check mode
    (1 = check the inputs, no file written) and <rest> (1, default: the relative displacements start at rest --
    their initial value, the zero-mean offset of the periodic SSI solution, is subtracted; 0: as computed)."""
    digits = _num(c, 1, "digits", integer=True, default=12)
    if not 6 <= digits <= 17:
        c.fail("<digits> must be 6..17")
    _flag(c, 2, "opmode")
    rest = _num(c, 3, "rest", integer=True, default=1)
    if rest not in (0, 1):
        c.fail("argument 3 <rest> must be 0 or 1")
    rec = LgoptRecord("LGOPT", _tokens(c, 3))
    if rec.is_default():
        c.model.options.delete_record("LGOPT")
    else:
        c.model.options.set_record(rec)
    c.confirm(f"LOADGEN output options: {digits} digits, opmode {rec.opmode}, start at rest {rest}")


# ======================================================================================
# deck from the records, LGLIST, RUNLOADGEN
# ======================================================================================
def build_deck(m, analysis: str):
    """The LOADGEN deck of model ``m`` for ``analysis`` STATIC / DYNAMIC from its records (D-LGN-01)."""
    from ...modules import loadgen as LG
    d = LG.new_deck()
    d["title"] = m.title or ""
    d["model"] = m.name or ""
    d["analysis"] = analysis
    st = _rec(m, "LOADGEN")
    dyn = _rec(m, "LOADGENDYN")
    d["data"] = DATA_WORDS[str(_get(st, 1, "2")).upper()]
    d["multi"] = int(parse_int(_get(st, 2, "0"))[0])
    d["masstype"] = MASS_WORDS[str(_get(st, 5, "1")).upper()]
    d["genmass"] = int(parse_int(_get(st, 6, "0"))[0])
    if analysis == "STATIC":
        d["rotdisp"] = int(parse_int(_get(st, 3, "0"))[0])
        d["rotacc"] = int(parse_int(_get(st, 4, "0"))[0])
        d["source"] = str(_get(st, 7, "RESULTS")).upper()
    else:
        d["rotdisp"] = int(parse_int(_get(dyn, 6, "0"))[0])
        d["rotacc"] = int(parse_int(_get(dyn, 7, "0"))[0])
        d["source"] = str(_get(dyn, 5, "RESULTS")).upper()
    d["method"] = str(_get(dyn, 3, "REL")).upper()
    d["refnode"] = int(parse_int(_get(dyn, 4, "0"))[0])
    zeta = parse_float(_get(dyn, 8, "0"))
    if zeta > 0:
        d["alpha"], d["beta"] = rayleigh(zeta, parse_float(_get(dyn, 9, "0")), parse_float(_get(dyn, 10, "0")))
    else:
        d["alpha"] = parse_float(_get(dyn, 1, "0"))
        d["beta"] = parse_float(_get(dyn, 2, "0"))
    d["groundfopt"] = int(parse_int(_get(dyn, 11, "0"))[0])
    d["groundmult"] = parse_float(_get(dyn, 12, "1"))
    for key, rec in m.options.entries("LGFILE"):
        key = str(key)
        name = rec.arg(2) or ""
        if key == "APDL" and analysis != "STATIC" or key == "APDLDYN" and analysis != "DYNAMIC":
            continue
        d[FILE_KEYS[key]] = name
    lists = lgnode_lists(m)
    for kind, tab in NODE_KINDS.items():
        d.tables[tab] = list(lists.get(kind, []))
    tr = _rec(m, "LGTIME")
    if tr is not None:
        toks = [t for t in tr.to_tokens()]
        crit = toks[0].upper()
        d["crit"] = crit
        if crit in ("TIME", "STEP"):
            d["times"] = [parse_float(t) for t in toks[1:] if t]
        else:
            d["ncrit"] = int(parse_int(tr.arg(2) or "1")[0])
            d["tsep"] = parse_float(tr.arg(3) or "0")
            if crit in ("ACC", "DISP"):
                d["critnode"] = int(parse_int(tr.arg(4))[0])
                d["critdof"] = int(parse_int(tr.arg(5))[0])
            elif crit.startswith("M") and tr.arg(4) is not None:
                d["refpoint"] = [parse_float(tr.arg(k)) for k in (4, 5, 6)]
    mp = _rec(m, "LGMAP")
    if mp is not None:
        d["mapmode"] = str(mp.arg(1) or "IDENTITY").upper()
        d["maptol"] = parse_float(mp.arg(2) or "0")
        d["mapfile"] = mp.arg(3) or ""
    op = _rec(m, "LGOPT")
    if op is not None:
        d["digits"] = int(parse_int(op.arg(1) or "12")[0])
        d["opmode"] = int(parse_int(op.arg(2) or "0")[0])
        d["rest"] = int(parse_int(op.arg(3) or "1")[0])
    return d


@command("LGLIST", tier=TIER, max_args=1)
def cmd_lglist(c):
    """LGLIST,[STATIC|DYNAMIC]: list the LOADGEN (Option A) settings and the deck RUNLOADGEN would write."""
    from ...modules import loadgen as LG
    analysis = c.word(1, "STATIC")
    if analysis not in ("STATIC", "DYNAMIC"):
        c.fail("argument 1 must be STATIC or DYNAMIC")
    m = c.model
    try:
        d = build_deck(m, analysis)
    except (NumberError, KeyError, ValueError) as exc:
        c.fail(f"invalid LOADGEN records: {exc}")
    lines = [f"LOADGEN (Option A) settings of model {m.name or '(no MDL)'}, {analysis} "
             f"(RUNLOADGEN,{analysis} writes {m.name or '<model>'}{LG.DECK_EXT} and runs LOADGEN):"]
    if analysis == "STATIC":
        lines.append(f"  data {d['data']} {LG.DATA_NAMES[d['data']]}; multiple files {d['multi']}; mass "
                     f"{LG.MASS_NAMES[d['masstype']]}, generate {d['genmass']}; source {d['source']}; rotations disp "
                     f"{d['rotdisp']} acc {d['rotacc']}")
        crit = d["crit"]
        if crit in ("TIME", "STEP"):
            lines.append(f"  critical {'times' if crit == 'TIME' else 'steps'}: " + ", ".join(f"{t:g}" for t in d["times"]))
        else:
            what = {"V": "base shear along the control direction"}.get(crit, crit)
            extra = (f" of node {d['critnode']} dof {d['critdof']}" if crit in ("ACC", "DISP") else
                     (f" about {d['refpoint']}" if d["refpoint"] else ""))
            lines.append(f"  critical times: {d['ncrit']} largest peak(s) of {what}{extra}, separation {d['tsep']:g} s")
    else:
        lines.append(f"  method {d['method']}; Rayleigh alpha {d['alpha']:.6g}, beta {d['beta']:.6g}; reference node "
                     f"{d['refnode'] or '0 (control motion)'}; source {d['source']}; rotations disp {d['rotdisp']} "
                     f"acc {d['rotacc']}")
    for p in LG.PARAMS:
        if p.name in ("ssipath", "housefile", "dispdir", "disprotdir", "accdir", "accrotdir", "ansyspath", "lumpfile",
                      "masterfile", "apdlfile", "groundfile") and d[p.name] != p.default:
            lines.append(f"  {p.name:<10s} = {d[p.name]}")
    lines.append(f"  node map: {d['mapmode']}" + (f", {d['mapfile']}" if d["mapfile"] else "")
                 + (f", tolerance {d['maptol']:g}" if d["maptol"] else ""))
    for kind, tab in NODE_KINDS.items():
        nodes = d.tables[tab]
        lines.append(f"  {NODE_KIND_NAMES[kind]}: " + (" ".join(compress_ids(nodes)) if nodes else
                                                       ("default (interaction nodes)" if kind == "D" else "none")))
    lines.append(f"  APDL digits {d['digits']}; opmode {d['opmode']}; relative displacements start at rest {d['rest']}")
    c.info("\n".join(lines))


@command("RUNLOADGEN", tier=TIER, max_args=2)
def cmd_runloadgen(c):
    """RUNLOADGEN,[STATIC|DYNAMIC],[model]: write <model>.lgn from the LOADGEN records and run LOADGEN in the
    model directory (the Ok of the ANSYS Eq. Static Load / ANSYS Dynamic Load dialogs; STATIC by default).
    The RUN<MODULE> order RUNLOADGEN,<model>,[STATIC|DYNAMIC] is accepted too."""
    from ...modules import loadgen as LG
    interp = c.interp
    try:
        model_first = bool(c.raw(1)) and parse_int(c.raw(1)) is not None
    except NumberError:
        model_first = False
    ka, km = (2, 1) if model_first else (1, 2)
    analysis = c.word(ka, "STATIC")
    if analysis not in ("STATIC", "DYNAMIC"):
        c.fail(f"argument {ka} must be STATIC or DYNAMIC (got '{c.raw(ka)}')")
    num = c.int(km, default=-1, what="model")
    if num is None or num < 0:
        num = interp.active_model
    if num not in interp.models:
        c.fail(f"model {num} is not in memory")
    m = interp.models[num]
    if not m.name or not m.path:
        c.fail("Model name/path not defined -- use MDL")
    mdir = Path(m.path)
    try:
        d = build_deck(m, analysis)
    except (NumberError, KeyError, ValueError) as exc:
        c.fail(f"invalid LOADGEN records: {exc}")
    ssidir = Path(d["ssipath"]).expanduser() if d["ssipath"] else mdir
    if not ssidir.is_absolute():
        ssidir = mdir / ssidir
    missing = []
    if not (ssidir / f"{m.name}.N4").exists():
        missing.append(f"{m.name}.N4 (FILE4) missing -- run HOUSE")
    if d["source"] == "FILE8" and not (ssidir / f"{m.name}.mot").exists():
        missing.append(f"{m.name}.mot missing -- source FILE8 needs the MOTION deck (AFWRITE with MOTION enabled)")
    if missing:
        for t in missing:
            c.error(t)
        raise CommandReported("LOADGEN: prerequisites missing")
    mdir.mkdir(parents=True, exist_ok=True)
    dp = LG.deck_path(mdir, m.name)
    LG.write_deck(dp, d, comments=[f"LOADGEN input written by RUNLOADGEN,{analysis} (SASSI-EDU Option A)"])
    c.info(f"RUNLOADGEN: {dp.name} written; LOADGEN {analysis} for model {m.name} in {mdir}")
    rc = LG.run_loadgen(m.name, mdir, deck_path=dp, echo=lambda line: c.info(line))
    listing = mdir / f"{m.name}_{LG.NAME}.out"
    if rc != 0:
        c.fail(f"LOADGEN finished with status FAILED ({rc}); see {listing.name}")
    c.confirm(f"LOADGEN finished (listing {listing.name})")
