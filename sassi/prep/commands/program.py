"""Programming commands: variables, loops, macros and random lists (manual sections 5.6, 5.9, 9.15;
spec 04 sections 5-6; spec 10 section 5; requirements 3.4.M).

A variable is a name (case-insensitive), a list of strings (1-based) and an integer counter that
VAR resets to 0.  The ``@`` / ``#`` substitution itself is done by the interpreter
(:meth:`sassi.prep.interpreter.Interpreter.substitute`, rule L13).

* FOREACH,<Var>,<command>: runs the command once per item of Var; the command text is
  substituted afresh at each iteration; a nested FOREACH over the same variable is an error and
  does not run (manual 5.9).
* LOADMACRO caches a macro file under an upper-case name (<= 50 characters); MACRO replaces
  ``$k$`` with argument k as plain text, also inside words, before variable processing; a missing
  argument gives an empty string and a warning (D-PAR-13).
* RND / ADDRND use numpy PCG64 seeded by RNDSEED (D-GEN-08); values are stored as strings
  (floats via repr, integers via str).  LOGNORM takes the mean and standard deviation of the
  lognormal variable itself: ``s^2 = ln(1 + (sd/mean)^2)``, ``mu = ln(mean) - s^2/2`` (spec 10 5.10).
"""
from __future__ import annotations

import math
import re

import numpy as np

from ...model.values import NumberError, parse_float, parse_int
from ..interpreter import LoopFrame, Macro, Variable, read_command_lines
from ..registry import command

_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
MACRO_NAME_MAX = 50
MACRO_ARG_MAX = 3000


def _varname(c, k: int = 1) -> str:
    name = c.raw(k)
    if name.startswith("@"):
        name = name[1:]
    if not name:
        c.fail("variable name required")
    if not _NAME_RE.match(name):
        c.fail(f"'{name}' is not a valid variable name (letters, digits, underscore)")
    return name


def _get(c, name: str) -> Variable:
    v = c.interp.variables.get(name.lower())
    if v is None:
        c.fail(f"variable {name} not defined")
    return v


def _set(c, name: str, items, reset: bool = True) -> Variable:
    vs = c.interp.variables
    old = vs.get(name.lower())
    v = Variable(name, [str(x) for x in items], 0 if (reset or old is None) else old.counter)
    vs[name.lower()] = v
    return v


# ======================================================================================
# Variables
# ======================================================================================
@command("VAR")
def cmd_var(c):
    """VAR,<Name>,<X1>,...,<Xn>: define or redefine a variable (the counter is reset to 0)."""
    name = _varname(c)
    items = [t.strip() for t in c.tokens[1:c.nargs]]
    v = _set(c, name, items)
    c.confirm(f"variable {v.name}: {len(items)} items, counter 0")


@command("SETVAR")
def cmd_setvar(c):
    """SETVAR,...: no effect on the model; its arguments are substituted, so counters change."""
    c.info(c.line.strip())


@command("SHOWVAR", max_args=1)
def cmd_showvar(c):
    """SHOWVAR,<varname>: show the items (1-based) and the counter of a variable."""
    v = _get(c, _varname(c))
    c.info(f"variable {v.name}: {len(v.items)} items, counter {v.counter}")
    for i, it in enumerate(v.items, start=1):
        c.info(f"  [{i}] {it}")


@command("VARLIST", max_args=0)
def cmd_varlist(c):
    """VARLIST: list the variables (value of single-item variables, item count otherwise)."""
    vs = c.interp.variables
    if not vs:
        c.info("no variables defined")
    for key in sorted(vs):
        v = vs[key]
        if len(v.items) == 1:
            c.info(f"{v.name} = {v.items[0]}   (counter {v.counter})")
        else:
            c.info(f"{v.name} : {len(v.items)} elements   (counter {v.counter})")


@command("LOADVAR", max_args=2)
def cmd_loadvar(c):
    """LOADVAR,<filename>,[name]: variable from a text file, one item per line (name = file stem)."""
    fname = c.str(1)
    if not fname:
        c.fail("file name required")
    p = c.input_path(fname)
    try:
        lines = read_command_lines(p)
    except OSError as exc:
        c.fail(f"cannot read {p}: {exc.strerror or exc}")
    name = _varname(c, 2) if c.given(2) else p.stem.upper()
    if not _NAME_RE.match(name):
        c.fail(f"file stem '{name}' is not a valid variable name; give the name as argument 2")
    items = [ln.strip() for ln in lines if ln.strip()]
    v = _set(c, name, items)
    c.confirm(f"variable {v.name}: {len(items)} items loaded from {p.name}")


@command("REDUCESET", max_args=2)
def cmd_reduceset(c):
    """REDUCESET,<var>,[sorttype]: sort ascending and remove duplicates (STRING default, INT, FLOAT)."""
    v = _get(c, _varname(c))
    kind = c.word(2, "STRING")
    if kind == "STRING":
        items = sorted(set(v.items))
    elif kind in ("INT", "FLOAT"):
        try:
            if kind == "INT":
                vals = []
                for s in v.items:
                    n, exact = parse_int(s)
                    if not exact:
                        raise NumberError(s)
                    vals.append(n)
                items = [str(x) for x in sorted(set(vals))]
            else:
                items = [repr(x) for x in sorted(set(parse_float(s) for s in v.items))]
        except NumberError as exc:
            c.fail(f"item {exc} is not {'an integer' if kind == 'INT' else 'a number'}; variable unchanged")
    else:
        c.fail("sort type must be STRING, INT or FLOAT")
    v.items = items
    v.counter = 0
    c.confirm(f"variable {v.name}: {len(items)} unique items")


# ======================================================================================
# Random lists
# ======================================================================================
_DIST_ALIASES = {"UNI": "UNI", "UNIINT": "UNIINT", "UNINT": "UNIINT", "NORM": "NORM", "LOGNORM": "LOGNORM",
                 "POISSON": "POISSON", "POSSION": "POISSON"}


def _samples(c, n: int):
    dist = _DIST_ALIASES.get(c.word(3))
    if dist is None:
        c.fail("distribution must be UNI, UNIINT, NORM, LOGNORM or POISSON")
    rng = c.interp.rng
    p1 = c.float(4, required=True, what="param1")
    if dist == "UNI":
        p2 = c.float(5, required=True, what="param2")
        if p2 < p1:
            c.fail("UNI: maximum < minimum")
        return [repr(float(x)) for x in rng.uniform(p1, p2, n)]
    if dist == "UNIINT":
        p2 = c.float(5, required=True, what="param2")
        lo, hi = int(math.ceil(p1)), int(math.floor(p2))
        if hi < lo:
            c.fail("UNIINT: empty integer range")
        return [str(int(x)) for x in rng.integers(lo, hi + 1, n)]
    if dist == "NORM":
        sd = c.float(5, required=True, what="param2")
        if sd < 0:
            c.fail("NORM: standard deviation < 0")
        return [repr(float(x)) for x in rng.normal(p1, sd, n)]
    if dist == "LOGNORM":
        sd = c.float(5, required=True, what="param2")
        if p1 <= 0 or sd < 0:
            c.fail("LOGNORM: mean must be > 0 and standard deviation >= 0")
        s2 = math.log(1.0 + (sd / p1) ** 2)
        mu = math.log(p1) - 0.5 * s2
        return [repr(float(x)) for x in rng.lognormal(mu, math.sqrt(s2), n)]
    if p1 < 0:
        c.fail("POISSON: mean must be >= 0")
    return [str(int(x)) for x in rng.poisson(p1, n)]


@command("RND")
def cmd_rnd(c):
    """RND,<var>,<numsamples>,<dist>,<param1>,...: fill a variable with random values (overwrites it)."""
    name = _varname(c)
    n = c.int(2, required=True, what="numsamples")
    if n < 0:
        c.fail("numsamples must be >= 0")
    v = _set(c, name, _samples(c, n))
    c.confirm(f"variable {v.name}: {n} random values ({c.word(3)})")


@command("ADDRND")
def cmd_addrnd(c):
    """ADDRND,<var>,<numsamples>,<dist>,<param1>,...: append random values until the list has numsamples items."""
    name = _varname(c)
    n = c.int(2, required=True, what="numsamples")
    old = c.interp.variables.get(name.lower())
    have = len(old.items) if old else 0
    add = max(0, n - have)
    new = _samples(c, add)
    if add == 0:
        c.warn(f"{name} already has {have} items: nothing appended (numsamples is the final count)")
    items = (old.items if old else []) + new
    v = _set(c, old.name if old else name, items, reset=False)
    c.confirm(f"variable {v.name}: {add} random values appended ({len(items)} items)")


@command("RNDSEED", max_args=1)
def cmd_rndseed(c):
    """RNDSEED,<seed>: seed the random generator of RND/ADDRND (positive integer)."""
    seed = c.int(1, required=True, what="seed")
    if seed <= 0:
        c.fail("the seed must be a positive integer")
    c.interp.rng = np.random.Generator(np.random.PCG64(seed))
    c.confirm(f"random seed {seed}")


# ======================================================================================
# Loops
# ======================================================================================
@command("FOREACH", raw=True)
def cmd_foreach(c):
    """FOREACH,<Var>,<Command>: run Command once per item of Var ('#' = 1-based loop index)."""
    rest = c.rest
    k = rest.find(",")
    if k < 0:
        c.fail("syntax FOREACH,<Var>,<Command>")
    name = rest[:k].strip()
    if name.startswith("@"):
        name = name[1:]
    body = rest[k + 1:].strip()
    if not name or not body:
        c.fail("syntax FOREACH,<Var>,<Command>")
    var = c.interp.variables.get(name.lower())
    if var is None:
        c.fail(f"variable {name} not defined")
    loops = c.interp.loops
    if any(fr.var == name.lower() for fr in loops):
        c.fail(f"nested loops cannot reuse the variable {var.name}; the loop is not run")
    frame = LoopFrame(name.lower(), 0)
    loops.append(frame)
    try:
        for i in range(1, len(var.items) + 1):
            frame.index = i
            c.interp.execute(body)
    finally:
        loops.pop()
    c.confirm(f"FOREACH {var.name}: {len(var.items)} iterations")


# ======================================================================================
# Macros
# ======================================================================================
@command("LOADMACRO", max_args=2)
def cmd_loadmacro(c):
    """LOADMACRO,<Name>,<MACROFILE>: load and cache a macro file under an upper-case name."""
    name = c.str(1).upper()
    if not name:
        c.fail("macro name required")
    if len(name) > MACRO_NAME_MAX:
        c.fail(f"macro names have at most {MACRO_NAME_MAX} characters")
    fname = c.str(2)
    if not fname:
        c.fail("macro file required")
    p = c.input_path(fname)
    try:
        lines = read_command_lines(p)
    except OSError as exc:
        c.fail(f"cannot read {p}: {exc.strerror or exc}")
    if name in c.interp.macros:
        c.warn(f"macro {name} redefined")
    c.interp.macros[name] = Macro(name, str(p), lines)
    c.confirm(f"macro {name} loaded from {p} ({len(lines)} lines)")


@command("MACRO")
def cmd_macro(c):
    """MACRO,<Name>,<1>,...,<N>: run a loaded macro, replacing `$k$` with argument k."""
    name = c.str(1).upper()
    if not name:
        c.fail("macro name required")
    args = [t.strip() for t in c.tokens[1:c.nargs]]
    for k, a in enumerate(args, start=1):
        if len(a) > MACRO_ARG_MAX:
            c.fail(f"argument {k} longer than {MACRO_ARG_MAX} characters")
    c.interp.run_macro(name, args)


@command("MACROLIST", max_args=0)
def cmd_macrolist(c):
    """MACROLIST: list the loaded macros and their files."""
    ms = c.interp.macros
    if not ms:
        c.info("no macros loaded")
    for name in sorted(ms):
        c.info(f"{name:<20} {ms[name].file}")
