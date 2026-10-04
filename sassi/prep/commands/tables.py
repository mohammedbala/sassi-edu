"""Property-table commands: M, L, R, SC, MXR/MXI/MXM, their deletion and listing commands
(manual section 9.4; spec 08 sections 5.1-5.35 and 6; requirements 3.4.D).

The tables store the raw values as entered; the conversions (material types 1/2/3, weight to
mass, complex moduli) are done when the analysis files are written (spec 08 section 1.6) and are
implemented once in :mod:`sassi.model.materials`.  The listing commands show the derived values
with the current HOUSE gravity so the user can check them.

"Last defined" in the listing defaults means the highest index (D-MDL-05).
"""
from __future__ import annotations

from ...model import fmt_num
from ...model.entities import BeamSection, MatrixProp, Material, SoilLayer, SpringProp
from ..registry import command
from . import range_args, select_ids


def _index(c) -> int:
    nm = c.int(1, required=True, what="nm")
    if nm < 1:
        c.fail("table numbers start at 1")
    return nm


def _damping(c, *vals: float) -> None:
    for v in vals:
        if v > 0.5:
            c.warn(f"damping ratio {fmt_num(v)} > 0.5: damping ratios are fractions (0.05 = 5 %)")
            return


# ======================================================================================
# Definitions
# ======================================================================================
@command("M", max_args=7)
def cmd_m(c):
    """M,<nm>,<val1>,<val2>,<weight>,<pdamp>,<sdamp>,<type>: material; type 1 (E, nu), 2 (M, G), 3 (Vp, Vs)."""
    nm = _index(c)
    v1 = c.float(2, required=True, what="val1")
    v2 = c.float(3, required=True, what="val2")
    w = c.float(4, required=True, what="weight")
    pd, sd = c.float(5, 0.0), c.float(6, 0.0)
    t = c.int(7, default=1)
    if t not in (1, 2, 3):
        c.fail("material <type> must be 1 (E, nu), 2 (M, G) or 3 (Vp, Vs)")
    _damping(c, pd, sd)
    existed = nm in c.model.materials
    c.model.materials[nm] = Material(nm, v1, v2, w, pd, sd, t)
    c.confirm(f"material {nm} {'redefined' if existed else 'defined'} (type {t})")


@command("L", max_args=7)
def cmd_l(c):
    """L,<nm>,<thick>,<weight>,<pveloc>,<sveloc>,<pdamp>,<sdamp>: soil layer (SITE and excavated soil)."""
    nm = _index(c)
    vals = [c.float(k, required=True, what=w) for k, w in zip(range(2, 6), ("thick", "weight", "pveloc", "sveloc"))]
    pd, sd = c.float(6, 0.0), c.float(7, 0.0)
    _damping(c, pd, sd)
    existed = nm in c.model.layers
    c.model.layers[nm] = SoilLayer(nm, vals[0], vals[1], vals[2], vals[3], pd, sd)
    c.confirm(f"soil layer {nm} {'redefined' if existed else 'defined'}")


@command("R", max_args=7)
def cmd_r(c):
    """R,<nm>,<axial>,<shear2>,<shear3>,<tors>,<flex2>,<flex3>: beam section (A, As2, As3, J, I2, I3)."""
    nm = _index(c)
    a = c.float(2, required=True, what="axial")
    rest = [c.float(k, 0.0) for k in range(3, 8)]
    existed = nm in c.model.sections
    c.model.sections[nm] = BeamSection(nm, a, *rest)
    c.confirm(f"beam section {nm} {'redefined' if existed else 'defined'}")


@command("SC", max_args=8)
def cmd_sc(c):
    """SC,<nm>,<scx>,<scy>,<scz>,<scxx>,<scyy>,<sczz>,<damp>: spring constants (global axes) and damping."""
    nm = _index(c)
    vals = [c.float(k, 0.0) for k in range(2, 9)]
    _damping(c, vals[6])
    existed = nm in c.model.springs
    c.model.springs[nm] = SpringProp(nm, *vals)
    c.confirm(f"spring property {nm} {'redefined' if existed else 'defined'}")


def _matrix_row(c, kind: str) -> None:
    p = _index(c)
    row = c.int(2, required=True, what="row")
    if not 1 <= row <= 12:
        c.fail("<row> must be 1..12")
    n = 13 - row
    if c.nargs - 2 > n:
        c.fail(f"row {row} of the upper triangle has {n} terms; {c.nargs - 2} given")
    terms = [c.float(k, 0.0) for k in range(3, 3 + n)]
    m = c.model
    prop = m.matrices.get(p)
    if prop is None:
        prop = m.matrices[p] = MatrixProp(p)
    prop.set_row(kind, row, terms)
    c.confirm(f"matrix property {p}: {('real stiffness', 'imaginary stiffness', 'mass')['RIM'.index(kind)]} row {row}")


@command("MXR", max_args=14)
def cmd_mxr(c):
    """MXR,<p>,<row>,<t1>,...,<t12>: real stiffness row of GENERAL matrix property p (upper triangle)."""
    _matrix_row(c, "R")


@command("MXI", max_args=14)
def cmd_mxi(c):
    """MXI,<p>,<row>,<t1>,...,<t12>: imaginary stiffness row of GENERAL matrix property p."""
    _matrix_row(c, "I")


@command("MXM", max_args=14)
def cmd_mxm(c):
    """MXM,<p>,<row>,<t1>,...,<t12>: mass (or weight, MOPT <matrix> = 1) row of GENERAL matrix property p."""
    _matrix_row(c, "M")


# ======================================================================================
# Deletion
# ======================================================================================
def _delete(c, table: dict, what: str, check: str) -> None:
    ids = range_args(c, 1, table)
    for i in ids:
        del table[i]
    c.confirm(f"{len(ids)} {what} deleted" + (f" (remaining references: CHECK {check})" if ids else ""))


@command("DELM", max_args=3)
def cmd_delm(c):
    """DELM,<m1>,[<m2>],[<inc>]: delete materials."""
    _delete(c, c.model.materials, "materials", "Error 13")


@command("DELL", max_args=3)
def cmd_dell(c):
    """DELL,<m1>,[<m2>],[<inc>]: delete soil layers."""
    _delete(c, c.model.layers, "soil layers", "Error 19")


@command("DELR", max_args=3)
def cmd_delr(c):
    """DELR,<r1>,[<r2>],[<inc>]: delete beam sections."""
    _delete(c, c.model.sections, "beam sections", "Error 26")


@command("DELSC", max_args=3)
def cmd_delsc(c):
    """DELSC,<r1>,[<r2>],[<inc>]: delete spring properties."""
    _delete(c, c.model.springs, "spring properties", "Error 33")


@command("MXDEL", max_args=3)
def cmd_mxdel(c):
    """MXDEL,<p1>,[<p2>],[<step>]: delete matrix properties."""
    _delete(c, c.model.matrices, "matrix properties", "Error 83")


# ======================================================================================
# Listing
# ======================================================================================
def _list_ids(c, table: dict):
    """``<m1>,[<m2>],[<step>]``: m2 defaults to the last defined (highest) index; no argument = all."""
    if not table:
        return []
    if not c.given(1):
        return sorted(table)
    a = c.int(1)
    b = c.int(2, default=max(table))
    return select_ids(table, a, b, c.int(3, default=1))


def _g(x: float) -> str:
    return f"{x:>12.5g}"


@command("MLIST", max_args=3)
def cmd_mlist(c):
    """MLIST,<m1>,[<m2>],[<step>]: list materials with derived E, nu, G, M, Vp, Vs, rho."""
    m = c.model
    g = m.gravity
    c.info(f"{'mat':>5} {'type':>4} {'val1':>12} {'val2':>12} {'weight':>12} {'pdamp':>8} {'sdamp':>8}")
    for i in _list_ids(c, m.materials):
        t = m.materials[i]
        c.info(f"{i:>5} {t.mtype:>4} {_g(t.val1)} {_g(t.val2)} {_g(t.weight)} {t.pdamp:>8.4g} {t.sdamp:>8.4g}")
        try:
            k = t.constants(g)
            c.info(f"      derived (g = {fmt_num(g)}): E {k.E:.6g}  nu {k.nu:.6g}  G {k.G:.6g}  M {k.M:.6g}  "
                   f"Vp {k.Vp:.6g}  Vs {k.Vs:.6g}  rho {k.rho:.6g}")
        except ValueError as exc:
            c.info(f"      derived values not available: {exc}")


@command("LLIST", max_args=3)
def cmd_llist(c):
    """LLIST,<m1>,[<m2>],[<step>]: list soil layers with derived G, M, rho."""
    m = c.model
    g = m.gravity
    c.info(f"{'layer':>5} {'thick':>12} {'weight':>12} {'Vp':>12} {'Vs':>12} {'pdamp':>8} {'sdamp':>8}")
    for i in _list_ids(c, m.layers):
        t = m.layers[i]
        c.info(f"{i:>5} {_g(t.thick)} {_g(t.weight)} {_g(t.vp)} {_g(t.vs)} {t.pdamp:>8.4g} {t.sdamp:>8.4g}")
        try:
            k = t.constants(g)
            c.info(f"      derived (g = {fmt_num(g)}): G {k.G:.6g}  M {k.M:.6g}  nu {k.nu:.6g}  rho {k.rho:.6g}")
        except ValueError as exc:
            c.info(f"      derived values not available: {exc}")


@command("RLIST", max_args=3)
def cmd_rlist(c):
    """RLIST,<r1>,[<r2>],[<step>]: list beam sections."""
    m = c.model
    c.info(f"{'sect':>5} {'A':>12} {'As2':>12} {'As3':>12} {'J':>12} {'I2':>12} {'I3':>12}")
    for i in _list_ids(c, m.sections):
        t = m.sections[i]
        c.info(f"{i:>5} {_g(t.axial)} {_g(t.shear2)} {_g(t.shear3)} {_g(t.tors)} {_g(t.flex2)} {_g(t.flex3)}")


@command("SCLIST", max_args=3)
def cmd_sclist(c):
    """SCLIST,<r1>,[<r2>],[<step>]: list spring properties."""
    m = c.model
    c.info(f"{'sc':>5} {'kx':>12} {'ky':>12} {'kz':>12} {'kxx':>12} {'kyy':>12} {'kzz':>12} {'damp':>8}")
    for i in _list_ids(c, m.springs):
        t = m.springs[i]
        c.info(f"{i:>5} " + " ".join(_g(v) for v in t.k) + f" {t.damp:>8.4g}")


@command("MXLIST", max_args=1)
def cmd_mxlist(c):
    """MXLIST,<p>: list the real stiffness, imaginary stiffness and mass/weight matrices of property p."""
    m = c.model
    p = c.int(1, required=True, what="p")
    prop = m.matrices.get(p)
    if prop is None:
        c.fail(f"matrix property {p} is not defined")
    units = "weight" if m.mopt.get("matrix") == 1 else "mass"
    for kind, title in (("R", "real stiffness K_R"), ("I", "imaginary stiffness K_I"), ("M", f"{units} matrix")):
        c.info(f"matrix property {p}: {title} (12 x 12, symmetric)")
        A = prop.full(kind)
        for r in range(12):
            c.info("  " + " ".join(f"{v:>10.4g}" for v in A[r]))
