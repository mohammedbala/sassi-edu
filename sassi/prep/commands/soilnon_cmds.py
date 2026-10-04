"""SOIL-NON commands: the nonlinear time-domain option of the SOIL module (requirements 1.2 and 3.4.O,
tier P2; spec 05a section 6.6, spec 11 sections 2.9-2.10; manual 9.17.12, 9.17.19, 9.17.20).

Commands
--------
``NLSOIL,<Opt>,<NSTimeSunInc>,<DispConv>,<ForceConv>,<EqualIt>,<BedInt>,<NLDampType>,<MMmult>,<SMmult>``
    Global SOIL-NON options (record setter; every argument defaults to 0; it does not change the
    linear / equivalent-linear soil data).  ``Opt`` 1 runs SOIL as SOIL-NON instead of SOIL-EQL;
    ``NSTimeSunInc`` sub-increments per input time step (0 = flexible sub-stepping, decision SN-9 of
    :mod:`sassi.core.soilnon`); ``DispConv`` / ``ForceConv`` displacement / force convergence errors
    (relative 2-norms, 0 = 1e-6, SN-10); ``EqualIt`` maximum Newton-Raphson equilibrium iterations per
    sub-step (0 = 25); ``BedInt`` bedrock interface 0 rigid (the base moves with the input motion, a
    within motion) / 1 viscoelastic (Joyner-Chen dashpot with the half-space properties, outcrop
    motion -- "currently disabled" in ACS SASSI V3, implemented here, SN-7); ``NLDampType`` small-strain
    damping 1 frequency independent, 2 visco-elastic (needs Vis), 3 Rayleigh (needs the multipliers;
    both 0 = set at f1 and 5 f1); ``MMmult`` / ``SMmult`` the Rayleigh mass and stiffness multipliers.
``NLSLAYER,<Num>,[curvefit],[B],[S],[refStrain],[Vis]``
    Nonlinear layer data set of SOIL sublayer ``Num`` (indexed setter; SPRO numbering, 1 = top, SN-1)
    for the hyperbolic law ``tau = G0 g / (1 + B (g / g_r)^S) + Vis dg/dt``: ``curvefit`` 1 fits the
    parameters to the sublayer's DYNP G/Gmax curve (B, S, refStrain and Vis are then ignored, manual);
    0 uses B (typically 1-1.4), S (typically 0.7-0.9) and refStrain (in percent, typically 0.03 %,
    SN-2).  B = 0 makes the sublayer linear elastic.  Every soil sublayer needs a set (Error 125); the
    half-space (last SPRO entry) needs none -- it is the elastic base (D-W3-10).
``DELNLS,<start>,[end],[stride]``
    Delete the NLSLAYER sets ``start, start+stride, ... <= end`` (end defaults to start; a blank stride
    or a stride < 1 uses 1 with a warning).  The linear soil data (L, SOIL, SPRO ...) are not touched
    (manual 9.17.12, whose syntax line is misprinted, spec 11 section 2.10).

AFWRITE writes the NLSOIL record and the NLSLAYER sets to ``<model>.nls`` next to the SOIL deck
(:func:`write_nls_file`; the SOIL deck schema has no SOIL-NON fields) and the SOIL module reads it.
All records survive WRITE -> INP (they are written in the NONLINEAR options section).
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional

from ...core import soilnon as SN
from ...model.options import Field, Record, register_record_type
from ...model.values import NumberError, parse_float
from ..options import get_entries, get_record
from ..registry import command


# ======================================================================================
# Records
# ======================================================================================
@register_record_type
class NlsoilRecord(Record):
    """``NLSOIL,<Opt>,<NSTimeSunInc>,<DispConv>,<ForceConv>,<EqualIt>,<BedInt>,<NLDampType>,<MMmult>,<SMmult>``."""
    COMMAND = "NLSOIL"
    FIELDS = (Field("opt", int, 0, "Nonlinear Time Domain (0 no, 1 yes)"),
              Field("nsub", int, 0, "Subincrements per Timestep (0 flexible)"),
              Field("dispconv", float, 0.0, "Displacement Convergence Error"),
              Field("forceconv", float, 0.0, "Force Convergence Error"),
              Field("equalit", int, 0, "Equilibrium Iterations"),
              Field("bedint", int, 0, "Bedrock Interface (0 rigid, 1 viscoelastic)"),
              Field("damptype", int, 0, "Damping Type (1 frequency independent, 2 visco-elastic, 3 Rayleigh)"),
              Field("mmmult", float, 0.0, "Mass Matrix Mult."),
              Field("smmult", float, 0.0, "Stiff Matrix Mult."))


@register_record_type
class NlslayerRecord(Record):
    """``NLSLAYER,<Num>,[curvefit],[B],[S],[refStrain],[Vis]``."""
    COMMAND = "NLSLAYER"
    FIELDS = (Field("num", int, 0, "Layer Number"),
              Field("curvefit", int, 0, "Curve Fit Hyperbolic Parameters (0 user, 1 fit)"),
              Field("b", float, 0.0, "Beta"),
              Field("s", float, 0.0, "S exponent"),
              Field("refstrain", float, 0.0, "Reference Strain (%)"),
              Field("vis", float, 0.0, "Viscosity"))


def _numbers(c, toks: List[str], fields) -> List[Optional[float]]:
    """Parsed values of the given tokens (None = blank); non-numbers and non-integers fail the command."""
    out: List[Optional[float]] = []
    for i, (t, f) in enumerate(zip(toks, fields), start=1):
        t = t.strip()
        if not t:
            out.append(None)
            continue
        try:
            v = parse_float(t)
        except NumberError:
            c.fail(f"argument {i} <{f.name}>: '{t}' is not a number")
        if f.type is int and v != round(v):
            c.fail(f"argument {i} <{f.name}>: {t} must be an integer")
        out.append(v)
    return out


def options_from_record(rec: Optional[Record]) -> SN.NLSoilOptions:
    """The NLSOIL record as :class:`sassi.core.soilnon.NLSoilOptions` (defaults when absent)."""
    if rec is None:
        return SN.NLSoilOptions()
    return SN.NLSoilOptions.from_values([rec.number(k) for k in range(1, 10)])


def layer_from_record(rec: Record) -> SN.NLLayerData:
    return SN.NLLayerData.from_values([rec.number(k) for k in range(1, 7)])


# ======================================================================================
# Handlers
# ======================================================================================
@command("NLSOIL", max_args=9)
def cmd_nlsoil(c):
    """NLSOIL,<Opt>,<NSTimeSunInc>,<DispConv>,<ForceConv>,<EqualIt>,<BedInt>,<NLDampType>,<MMmult>,<SMmult>: global
    SOIL-NON options (Opt 1 = nonlinear time-domain SOIL; sub-increments, convergence errors, equilibrium
    iterations, bedrock 0 rigid / 1 viscoelastic, damping type 1 frequency independent / 2 visco-elastic /
    3 Rayleigh, Rayleigh mass and stiffness multipliers; all default 0)."""
    toks = list(c.tokens[:9])
    vals = _numbers(c, toks, NlsoilRecord.FIELDS)
    o = SN.NLSoilOptions.from_values(vals)
    probs = o.problems()
    if probs:
        c.fail("; ".join(probs))
    rec = NlsoilRecord("NLSOIL", toks)
    c.model.options.set_record(rec)
    if o.opt == 1:
        if o.bedint == 1:
            c.warn("BedInt 1 (viscoelastic bedrock) is 'currently disabled' in ACS SASSI V3; SASSI-EDU implements it "
                   "with the Joyner & Chen (1975) dashpot (the input is then the outcrop motion at bedrock)")
        if o.damptype == 0:
            c.info("NLSOIL: NLDampType 0 is not a documented choice; SOIL uses 1 (frequency independent)")
        if o.damptype == 3 and o.mmmult == 0 and o.smmult == 0:
            c.info("NLSOIL: Rayleigh multipliers both 0: SOIL sets them for the small-strain damping at f1 and 5 f1")
        nspro = len([k for k, _ in get_entries(c.model, "SPRO") if isinstance(k, int)])
        if nspro:
            sx = get_record(c.model, "SOILX")
            cl = int(sx.cl) if int(sx.cl) > 0 else int(get_record(c.model, "SITE").cl)
            if cl != nspro:
                c.info(f"NLSOIL: SOIL-NON needs the control motion at bedrock -- the SOIL control layer is {cl}, "
                       f"the half-space is SPRO sublayer {nspro} (SOILX <cl>)")
    c.confirm(f"NLSOIL: {'nonlinear time domain (SOIL-NON)' if o.opt == 1 else 'equivalent linear (SOIL-EQL)'}"
              + (f", bedrock {SN.BEDROCK_TYPES[o.bedint].split(' (')[0]}, damping type {o.damptype or 1}"
                 if o.opt == 1 else ""))


@command("NLSLAYER", max_args=6)
def cmd_nlslayer(c):
    """NLSLAYER,<Num>,[curvefit],[B],[S],[refStrain],[Vis]: SOIL-NON data of soil sublayer <Num> -- curvefit 1 fits
    the MKZ parameters to the sublayer's DYNP G/Gmax curve (B, S, refStrain, Vis ignored), 0 uses Beta, S exponent,
    Reference Strain (in %) and Viscosity of tau = G0 g / (1 + B (g/g_r)^S) + Vis dg/dt (B = 0: linear sublayer)."""
    toks = list(c.tokens[:6])
    vals = _numbers(c, toks, NlslayerRecord.FIELDS)
    if not vals or vals[0] is None:
        c.fail("the soil layer number <Num> is required")
    num = int(round(vals[0]))
    if num < 1:
        c.fail(f"<Num> = {num}: soil sublayer numbers start at 1")
    lay = SN.NLLayerData.from_values(vals)
    if lay.curvefit not in (0, 1):
        c.fail(f"<curvefit> = {lay.curvefit} (0 parameters given, 1 curve fit)")
    if lay.beta < 0 or lay.s < 0 or lay.refstrain < 0 or lay.vis < 0:
        c.fail("Beta, S exponent, Reference Strain and Viscosity must be >= 0")
    rec = NlslayerRecord("NLSLAYER", toks)
    c.model.options.set_entry("NLSLAYER", num, rec)
    if lay.curvefit == 1:
        if any(v for v in (lay.beta, lay.s, lay.refstrain, lay.vis)):
            c.info(f"NLSLAYER {num}: curve fit -- Beta, S exponent, Reference Strain and Viscosity are ignored (manual)")
        c.confirm(f"NLSLAYER {num}: parameters fitted to the DYNP G/Gmax curve of sublayer {num}")
        return
    if lay.beta == 0:
        c.warn(f"NLSLAYER {num}: Beta = 0 -- the sublayer is linear elastic in SOIL-NON")
    elif lay.s == 0 or lay.refstrain == 0:
        c.warn(f"NLSLAYER {num}: curvefit 0 needs Beta, S exponent and Reference Strain > 0 (SOIL stops otherwise)")
    if lay.s > 1.0:
        c.warn(f"NLSLAYER {num}: S = {lay.s:g} > 1 -- the backbone stress decreases at large strains (typical 0.7-0.9)")
    c.confirm(f"NLSLAYER {num}: Beta {lay.beta:g}, S {lay.s:g}, Reference Strain {lay.refstrain:g} %, "
              f"Viscosity {lay.vis:g}")


@command("DELNLS", max_args=3)
def cmd_delnls(c):
    """DELNLS,<start>,[end],[stride]: delete NLSLAYER sets start, start+stride, ... <= end (end = start by default;
    a blank stride or a stride < 1 uses 1 with a warning); linear soil data are not changed."""
    start = c.int(1)
    if start is None:
        c.fail("the first soil layer set <start> is required")
    end = c.int(2, default=start)
    stride = c.int(3, default=None)
    if stride is None or stride < 1:
        c.warn("default stride 1 used")
        stride = 1
    opts = c.model.options
    deleted = [k for k in range(start, end + 1, stride) if opts.delete_entry("NLSLAYER", k)]
    if deleted:
        c.confirm(f"DELNLS: {len(deleted)} nonlinear soil layer set(s) deleted ({', '.join(map(str, deleted))})")
    else:
        c.info(f"DELNLS: no NLSLAYER set in {start}..{end} (stride {stride})")


# ======================================================================================
# Model -> <model>.nls (AFWRITE)
# ======================================================================================
def nls_from_model(model) -> Optional[SN.NlsData]:
    """The SOIL-NON data of the model, or None without an NLSOIL record."""
    rec = model.options.record("NLSOIL")
    if rec is None:
        return None
    layers = {}
    for k, r in model.options.entries("NLSLAYER"):
        try:
            lay = layer_from_record(r)
        except (TypeError, ValueError):
            continue
        layers[lay.num] = lay
    return SN.NlsData(options_from_record(rec), layers, model.name)


def write_nls_file(model, mdir: Path, res=None) -> Optional[Path]:
    """Write ``<model>.nls`` next to the SOIL deck (AFWRITE, requirements 3.4.O).

    Without an NLSOIL record an existing ``.nls`` is stale and is renamed ``.nls.bak`` (D-AFW-01), so SOIL
    runs SOIL-EQL.  ``res`` is the :class:`sassi.prep.afwrite.AfwriteResult` (notes, warnings, model hash).
    """
    path = Path(mdir) / f"{model.name}.nls"
    data = nls_from_model(model)
    if data is None:
        if path.exists():
            bak = path.with_name(path.name + ".bak")
            os.replace(path, bak)
            if res is not None:
                res.warnings.append(f"SOIL: stale {path.name} renamed {bak.name} (no NLSOIL record: SOIL-EQL)")
        return None
    data.model_hash = getattr(res, "model_hash", "") if res is not None else ""
    SN.write_nls(path, data)
    if res is not None:
        o = data.options
        res.notes.append(f"{path.name} written (SOIL-NON options: NLSOIL Opt {o.opt}, {len(data.layers)} NLSLAYER "
                         f"set(s)); SOIL runs {'SOIL-NON' if o.opt == 1 else 'SOIL-EQL'}")
    return path
