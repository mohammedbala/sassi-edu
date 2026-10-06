"""Options dialogs of the GUI: Options > Model, Write, Check and Analysis (requirements 5.4).

The dialogs are *driven by the typed option records* of :mod:`sassi.prep.options` (field names =
command argument names, defaults = new-model defaults D-UI-03, dialog labels).  This module

1. describes the layout of every dialog as JSON (:func:`form`) -- tabs, group boxes and fields,
   each field bound to **one** storage location (a *path*, below), so a shared variable such as the
   time step is one value shown on several tabs (spec 07 section 3);
2. reads the current values of the active model (:func:`values`);
3. turns the values of a dialog commit into the equivalent command text (:func:`commit_commands`):
   only the records / entries / lists that changed produce a command, and the commands are those a
   user would type (rule L17, so the Command History replays the dialog);
4. validates the commit with the same field-local CHECK rules as the option commands (UI-06): an
   invalid value refuses the whole commit with the Chapter 10 message text (``Error 47 : ...``).

Value paths
-----------
==================  ===========================================================================
``REC.field``       field of a record setter (SITE, HOUSE, ..., X-commands SITEX, HOUSEX ...)
``IDX[sel].field``  field of the entry of an indexed setter whose key is the value of the dialog
                    selector ``sel`` (WAVE[wave], SPRO[layer], RSIN[spec], ME[motion] ...)
``$NAME``           string setter (THFILE, THTIT, EQTIT, RELFILE)
``#NAME``           list command as text (DAMP, TOPL; ``#AMP[motion]`` per motion)
``@NAME``           request table (NOUT, EOUT, RDND) or indexed table (CORR)
``%NAME...``        command records outside OPTION_SPECS (Option NON EQL / P / S / BBC, SOIL-NON
                    NLSOIL / NLSLAYER, Option A LOADGEN ...): ``%NAME.field``, the table ``%NAME``
                    and ``%NAME[sel].field``, see :mod:`sassi.ui.cmdrecords`
==================  ===========================================================================

Modules > ANSYS Eq. Static Load / ANSYS Dynamic Load are the one-page dialogs LOADGEN and LOADGENDYN
(:data:`MODULE_DIALOGS`): their OK submits the LOADGEN record commands, their Run also RUNLOADGEN.

Shown = used
------------
A value displayed for something that is *not stored* is what AFWRITE writes for it, so the dialog
never shows an analysis different from the one that runs: the implicit vertical SV (SH) wave field
while no WAVE entry exists (and it is stored with the first other wave page); SOIL layer output
requests "not requested" for a layer without entry; SITE Frequency 2 blank = NFFT/2; the
stochastic incoherency input read and validated with the D-INC-02 rule.  A blank input file with a built-in
default (requirements section 7.19, :mod:`sassi.prep.defaults`) shows that default as its placeholder
("built-in: RG 1.60, 0.30 g (@rg160h_030g.acc)"), computed with the policy AFWRITE applies
(``context.input_defaults`` and the browser's ``defaultPlaceholder``), and the SOIL profile of a model without
SPRO entries shows the default profile on every layer page; it is stored with the first SPRO page the user
changes (like the implicit wave field).  File fields get a Library button listing the built-in files that
fit them (``library``: record, load, spectrum, psd).

Enable/disable rules (D-UI-04) are carried by the fields as ``enable`` / ``visible`` conditions
``[path, op, value]`` (ops ``==``, ``!=``, ``in``, ``notin``) and evaluated by the browser; the
rules that *force* a value (2D -> coherent; unlagged coherency model 2-7 or multiple excitation ->
wave passage on) are also applied here on commit, with a message.
"""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ..conventions import ELEMENT_COMPONENTS, ELEMENT_TYPE_NAMES
from ..model.values import NumberError, fmt_num, parse_float, parse_int
from ..io import library as LIB
from ..prep import defaults as DEF
from ..prep.check import EDU_TEXT, ERRORS, WARNINGS, CheckOptions, format_title
from ..prep.lexer import join_command
from ..prep.options import (OPTION_SPECS, TABS, get_entries, get_record, problem_fmt, record_class,
                            to_command, typed)
from . import cmdrecords as CR

# ======================================================================================
# Field constructors
# ======================================================================================
SOL_DC = [[0, "Solution"], [1, "Data Check"]]
TYPE_CH = [[0, "Seismic"], [1, "Foundation Vibration"]]
INTERP_CH = [[0, "0"], [1, "1"], [2, "2"], [3, "3"], [4, "4"], [5, "5"], [6, "6 (complex bicubic spline)"]]
DIR_CH = [[1, "X"], [2, "Y"], [3, "Z"], [4, "XX"], [5, "YY"], [6, "ZZ"]]
NOUT_FLAGS = ["Printed Plot of Transfer Function", "Save Time History of Requested Response",
              "Plot Time History of Requested Response", "Plot Acceleration and Velocity R.S.",
              "Save Acceleration and Velocity R.S.", "Print Maximum Requested Response"]
WAVE_NAMES = {1: "R", 2: "SV", 3: "P", 4: "SH", 5: "L"}
#: wave pages per <wopt> (spec 05a section 7.2)
WAVE_PAGES = {0: [1, 2, 3], 1: [4, 5]}
WAVE_OPT = {1: [[0, "No Wave Field"], [1, "Shortest Wavelength"], [2, "Least Decay"]]}
for _t in (2, 3, 4, 5):
    WAVE_OPT[_t] = [[0, f"No {WAVE_NAMES[_t]}-Wave Field"], [1, f"{WAVE_NAMES[_t]}-Wave Field"]]

#: record fields whose blank token means "computed" (the dialog shows them empty, not their default):
#: EQUAKE <nrfreq> = records of RSIN 1; SITE <freq2> = NFFT/2 (requirements 5.4 default, AFWRITE writes
#: nft // 2 for a blank <freq2>), so an unrelated SITE edit must not freeze it at the record default
BLANK_DEFAULT = {("EQUAKE", "nrfreq"): "records of RSIN 1", ("SITE", "freq2"): "NFFT/2"}

#: SOIL layer output requests (SACC, SRS, SSTR, SSAF, SFOU): AFWRITE writes a row only for a stored
#: entry and SOIL computes nothing for the other layers, so a layer without an entry is shown as
#: "not requested"; the record defaults (requirements 5.4: Compute Maximum, Save RS, compute
#: stresses + strains) are what the browser preselects when the user first requests output there
LAYER_REQUESTS = ("SACC", "SRS", "SSTR", "SSAF", "SFOU")

#: 5.4 default of the Random Phase Angle once "Stochastically Simulated Incoherency Input" is chosen
STOCHASTIC_PHASE = 180.0


def f(path: str, label: Optional[str] = None, w: Optional[str] = None, **kw) -> Dict[str, Any]:
    """A dialog field bound to ``path``; ``w`` the widget (int real text path check radio select ...)."""
    d: Dict[str, Any] = {"path": path}
    if label is not None:
        d["label"] = label
    if w is not None:
        d["w"] = w
    d.update(kw)
    return d


def radio(path: str, label: str, choices, **kw) -> Dict[str, Any]:
    return f(path, label, "radio", choices=choices, **kw)


def check(path: str, label: Optional[str] = None, **kw) -> Dict[str, Any]:
    return f(path, label, "check", **kw)


def sel(name: str, label: str, lo: int = 1, hi: Optional[int] = None, **kw) -> Dict[str, Any]:
    """A selector (spin box) choosing the key of the indexed entries shown below it."""
    d = {"w": "selector", "sel": name, "label": label, "min": lo, "max": hi}
    d.update(kw)
    return d


def button(label: str, action: str, **kw) -> Dict[str, Any]:
    d = {"w": "button", "label": label, "action": action}
    d.update(kw)
    return d


def greyed(label: str, value: str = "", note: str = "") -> Dict[str, Any]:
    """A control shown disabled (not available in this version / tier P2)."""
    return {"w": "disabled", "label": label, "value": value, "note": note}


def info(label: str, context: str = "", keys_of: str = "", **kw) -> Dict[str, Any]:
    """Read-only text: ``context[<context>]`` of the values, or the keys of a keyed family (``%BBC``)."""
    d = {"w": "info", "label": label, "context": context, "keys_of": keys_of}
    d.update(kw)
    return d


def group(title: str, items: List[Dict[str, Any]], col: int = 0, note: str = "", **kw) -> Dict[str, Any]:
    """A group box in column ``col`` (0-2); ``full=True`` puts it below the columns, across the page."""
    d = {"title": title, "items": items, "col": col, "note": note}
    d.update(kw)
    return d


def _acc_history(freq: bool = False) -> List[Dict[str, Any]]:
    """The Acceleration Time History Data group shared by MOTION, STRESS and RELDISP (spec 07 sec. 3)."""
    items = [f("SITE.nft", "Nr. of Fourier Components"), f("SITE.delt", "Time Step of Control Motion")]
    if freq:
        items.append(f("SITE.freq", "Frequency Set Number"))
    items += [f("MOTION.mult", "Multiplication Factor"), f("MOTION.max", "Max Value for Time History"),
              f("MOTION.rec1", "First Record"), f("MOTION.rec2", "Last Record (0 = last)"),
              f("$THTIT", "Title", "text"),
              f("$THFILE", "File", "path", edit=True, library=["record", "load"], default_ph="thfile"),
              check("MOTION.fopt", "File Contains Pairs Time Step - Accel.")]
    return items


INCOH_ON = ["HOUSE.coh", "==", 1]


# ======================================================================================
# Options > Analysis: tabs EQUAKE ... AFWRITE (requirements 5.4)
# ======================================================================================
def analysis_tabs() -> List[Dict[str, Any]]:
    """The Analysis Options window: one entry per tab in the manual's order (TABS)."""
    tabs: Dict[str, List[Dict[str, Any]]] = {}
    tabs["EQUAKE"] = [
        group("Spectrum Files", [
            sel("spec", "Spectrum Number", 1, 3),
            f("RSIN[spec].file", "Spectrum Input File", "path", edit=True, library=["spectrum"], default_ph="rsin"),
            f("RSOUT[spec].file", "Spectrum Output File", "path", default_ph="rsout"),
            f("ACCOUT[spec].file", "Acceleration Output File", "path", default_ph="accout")]),
        group("Optional Spectrum Files", [
            f("EQUAKE.accopt", "", "excl", choices=[[1, "Accel. Record"], [2, "External Accel"]]),
            f("ACCIN[spec].file", "Acceleration Input File", "path", enable=["EQUAKE.accopt", "!=", 0],
              library=["record"])]),
        group("Target PSD", [
            check("EQUAKE.tpsd", "Use Target PSD"),
            f("TPSD[spec].file", "PSD File", "path", enable=["EQUAKE.tpsd", "==", 1], library=["psd"])]),
        group("Simulation Parameters", [
            f("EQUAKE.nrfreq", "Number of Frequencies", placeholder="records of RSIN 1"),
            f("EQUAKE.rand", "Initial Random SEED"), f("EQUAKE.damp", "Damping Value"),
            f("SITE.delt", "Time Step"), f("EQUAKE.dur", "Total Duration"), f("EQUAKE.seeds", "Number of SEEDs"),
            f("HOUSE.gravity", "Acceleration of Gravity (units: 32.2 BS, 9.81 IS)")], col=1),
        group("Correlation", [
            check("EQUAKE.corr", "Correlated"),
            f("@CORR", "Correlation", "table", key="no", enable=["EQUAKE.corr", "==", 1],
              columns=[{"name": "no", "label": "Row", "w": "int"}, {"name": "time", "label": "Time", "w": "real"},
                       {"name": "val", "label": "Corr.", "w": "real"}])], col=1, note="CORR is tier P1"),
        group("", [f("$EQTIT", "Spectra Title", "text")], col=1),
    ]
    tabs["SOIL"] = [
        group("Input Motion", [
            f("SITE.nft", "Number of Fourier Components"), f("SITE.delt", "Time Step of Input Motion"),
            f("SOIL.nrval", "Number of Values"),
            f("SOILX.mult", "Multiplication Factor"), f("SOILX.max", "Max Value for Time History (g)"),
            f("SOIL.grav", "Gravity Accel. (ft/s^2 or m/s^2) (used for free-field analysis)"),
            f("SOIL.header", "Number of Header Lines"),
            radio("SOILX.indir", "Input Direction", [[0, "Horizontal (Vs)"], [1, "Vertical (Vp)"]]),
            f("SITE.cl", "Control Point Layer"),
            f("SOILX.cl", "SOIL-only Control Point Layer (0 = SITE value)"),
            f("$THFILE", "File", "path", edit=True, library=["record"], default_ph="soil_thfile"),
            f("SOILX.file", "SOIL-only File (blank = File above)", "path", library=["record"]),
            check("SOIL.outcrop", "Assign as Outcrop Motion")]),
        group("Iteration Parameters (Equivalent-Linear Soil Behavior)", [
            check("SOIL.save", "Save Strain-Compatible Soil Properties"),
            f("SOIL.iter", "Number of Iterations"), f("SOIL.ratio", "Equiv. Uniform / Max Strain")]),
        group("Soil Profile", [
            sel("layer", "Layer Number", 1, 200),
            f("SPRO[layer].prop", "Property Number"),
            f("SPRO[layer].dynprop", "Dynamic Soil Property", "dynp"),
            info("Without SPRO", "spro_default",
                 tip="no SPRO entry stored: the default profile AFWRITE writes (EDU-29); it is stored with the first "
                     "layer you change")], col=1),
        group("Accelerations", [
            radio("SACC[layer].opt", "", [[0, "No Computation"], [1, "Compute Maximum"],
                                          [2, "Compute Maximum + Time History"]]),
            check("SACC[layer].outcrop", "Outcropping")], col=1),
        group("Response Spectrum", [
            check("SRS[layer].save", "Save Response Spectrum"), check("SRS[layer].outcrop", "Outcropping"),
            f("SOIL.gravmult", "Multiplier for Acceleration of Gravity"),
            f("#DAMP", "Damping Ratios", "list", tip="fractions, e.g. 0.02 0.05 (at most 10)")], col=1),
        group("Stresses / Strains", [
            check("SSTR[layer].opt1", "Compute Stresses"), check("SSTR[layer].opt2", "Save Stress Time History"),
            check("SSTR[layer].opt3", "Compute Strains"), check("SSTR[layer].opt4", "Save Strain Time History")],
            col=2),
        group("Spectral Amplification Factor", [
            check("SSAF[layer].save", "Save Spectral Amplification Factor"),
            check("SSAF[layer].outcrop1", "Outcropping of First Layer"),
            check("SSAF[layer].outcrop2", "Outcropping of Second Layer"),
            f("SSAF[layer].layer2", "Second Layer Number"), f("SSAF[layer].freqstep", "Frequency Step"),
            f("SSAF[layer].title", "Title", "text")], col=2),
        group("Fourier Spectrum", [greyed("Compute Fourier Spectrum"), greyed("Save to File"), greyed("Outcropping"),
                                   greyed("Nr. of Smoothings", "0"), greyed("Nr. of Values to be Saved", "0")],
              col=2, note="Not usable in this version (compute Fourier spectra with EQUAKE)"),
        group("Nonlinear Soil Behavior", [
            check("%NLSOIL.opt", "Nonlinear Time Domain", tip="SOIL-NON (hyperbolic model, time domain) instead of "
                                                              "SOIL-EQL; the control motion must be at bedrock"),
            f("%NLSOIL.nsub", "Subincrements per Timestep", "int", tip="0 = flexible sub-stepping"),
            f("%NLSOIL.dispconv", "Displacement Convergence Error", "real", tip="relative 2-norm; 0 = 1e-6"),
            f("%NLSOIL.forceconv", "Force Convergence Error", "real", tip="relative 2-norm; 0 = 1e-6"),
            f("%NLSOIL.equalit", "Equilibrium Iterations", "int",
              tip="Newton-Raphson iterations per sub-step before SOIL stops; 0 = 25"),
            radio("%NLSOIL.bedint", "Bedrock Interface", [[0, "Rigid"], [1, "Viscoelastic"]],
                  tip="Viscoelastic: Joyner-Chen dashpot with the half-space properties (outcrop input motion); "
                      "'not applicable' in ACS SASSI V3, implemented by SASSI-EDU"),
            f("%NLSOIL.damptype", "Damping Type", "select",
              choices=[[0, "0 (not set: 1 is used)"], [1, "1 Frequency independent"],
                       [2, "2 Visco-elastic (Viscosity)"], [3, "3 Rayleigh (multipliers)"]]),
            f("%NLSOIL.mmmult", "Mass Matrix Mult.", "real", enable=["%NLSOIL.damptype", "==", 3],
              tip="Rayleigh alpha (C = alpha M + beta K); both 0: set at f1 and 5 f1"),
            f("%NLSOIL.smmult", "Stiff Matrix Mult.", "real", enable=["%NLSOIL.damptype", "==", 3])],
            col=2, note="NLSOIL (SOIL-NON, tier P2); AFWRITE writes it with the NLSLAYER sets to <model>.nls"),
        group("Nonlinear Soil Layers (NLSLAYER)", [
            f("%NLSLAYER", "", "table", key="num",
              columns=[{"name": "num", "label": "Layer", "w": "int"},
                       {"name": "curvefit", "label": "Curve Fit (DYNP G/Gmax)", "w": "check"},
                       {"name": "b", "label": "Beta", "w": "real"},
                       {"name": "s", "label": "S exponent", "w": "real"},
                       {"name": "refstrain", "label": "Reference Strain (%)", "w": "real"},
                       {"name": "vis", "label": "Viscosity", "w": "real"}],
              fill={"context": "spro_layers", "label": "Add the SOIL sublayers without a row",
                    "row": {"curvefit": 1, "b": 0, "s": 0, "refstrain": 0, "vis": 0}},
              note="OK: DELNLS for a removed row, NLSLAYER for a new or changed row")],
            full=True,
            note="One row per SOIL sublayer (SPRO numbering, 1 = top; the half-space too, Error 125). "
                 "tau = G0 g / (1 + Beta (g/g_r)^S) + Viscosity dg/dt; Curve Fit: the parameters are fitted to "
                 "the sublayer's G/Gmax curve and the values given are ignored; Beta = 0: linear sublayer."),
    ]
    tabs["SITE"] = [
        group("Operation Mode", [
            radio("SITEX.soilmode", "", [[0, "Linear Soil"], [1, "Non-Linear Soil"]]),
            check("SITE.mode1", "Mode 1"), check("SITE.mode2", "Mode 2"),
            radio("SITE.opmode", "Run", SOL_DC)]),
        group("Mode 1", [
            f("HOUSE.gravity", "Gravity Accel. (ft/s^2 or m/s^2)"), f("SITE.fstep", "Frequency Step"),
            f("SITE.delt", "Time Step Control Motion"), f("SITE.nft", "Nr. of Fourier Component"),
            f("SITE.freq", "Frequency Set Number"), f("SITE.nl", "Number of Generated Layers"),
            f("SITE.hs", "Halfspace Layer"),
            f("#TOPL", "Top Layers", "list", rows=4, tip="L layer numbers, topmost first (at most 200)")]),
        group("Mode 2", [
            radio("SITE.wopt", "", [[0, "R-, SV- and P-Waves"], [1, "SH- and L-Waves"]]),
            sel("wave", "Wave", choices_by={"path": "SITE.wopt",
                                             "map": {str(k): [[t, f"{WAVE_NAMES[t]}-Wave"] for t in v]
                                                     for k, v in WAVE_PAGES.items()}}),
            f("WAVE[wave].opt", "", "radio", choices_by_key=WAVE_OPT),
            f("WAVE[wave].ratio1", "Wave Ratio 1"), f("WAVE[wave].ratio2", "Wave Ratio 2"),
            f("WAVE[wave].angle", "Incident Angle", visible=["sel.wave", "in", [2, 3, 4]]),
            f("SITE.freq1", "Frequency 1", hz=True),
            f("SITE.freq2", "Frequency 2", hz=True, placeholder="NFFT/2", blank_half="SITE.nft",
              tip="blank = NFFT/2 (follows Nr. of Fourier Components)"),
            f("SITE.cl", "Control Point Layer"),
            radio("SITE.cm", "Direction", [[0, "X"], [1, "Y"], [2, "Z"]])], col=1),
    ]
    tabs["POINT"] = [
        group("POINT Module Options", [
            radio("POINT.opmode", "Operation Mode", SOL_DC),
            f("POINT.layer", "Number of Embedment Soil Layers"),
            f("POINT.rad", "Point Load Central Zone Radius"),
            button("From mesh (RADIUS)", "radius",
                   tip="RADIUS: average r_e = 0.9 sqrt(A_plan) of the excavation elements (D-PNT-05)")]),
    ]
    tabs["HOUSE"] = [
        group("Model", [
            radio("HOUSE.opmode", "Operation Mode", SOL_DC),
            radio("HOUSE.dim", "Dimension of Analysis", [[0, "1D"], [1, "2D"], [2, "3D"]], disabled_choices=[0]),
            radio("HOUSE.imp", "Flexible Volume Method", [[0, "Flexible Volume (FV)"], [1, "Fast Flexible Volume (FFV)"],
                                                         [2, "Flexible Interface (FI)"]]),
            f("HOUSE.gravity", "Acceleration of Gravity"), f("HOUSE.gelev", "Ground Elevation"),
            check("HOUSEX.nlssi", "Non-Linear SSI"), button("Input Data (.pin)", "edit:{model}.pin"),
            check("HOUSEX.optimize", "Optimize Model")]),
        group("Soil Motion", [
            radio("HOUSE.coh", "", [[0, "Coherent"], [1, "Incoherent"]],
                  disable_choice_if={"1": ["HOUSE.dim", "==", 1]}),
            f("INCOH.gammax", "Coherence Parameter X Dir", enable=INCOH_ON),
            f("INCOH.gammay", "Coherence Parameter Y Dir", enable=INCOH_ON),
            f("INCOH.gammaz", "Coherence Parameter Z Dir", enable=INCOH_ON),
            f("INCOH.alpha", "Alpha Directionality Factor (Vs for model 1)", enable=INCOH_ON),
            f("INCOH.ngp", "Number of Embedded Layers", enable=INCOH_ON),
            f("SITE.delt", "Time Step of Seismic Motion"), f("SITE.nft", "Nr. of Fourier Components"),
            f("SITE.freq", "Frequency Set Number"),
            f("INCOH.nmodes", "Number of Incoh. Modes", enable=INCOH_ON),
            check("INCOH.ipr", "Print Coherency Matrix", enable=INCOH_ON)], col=1,
            note="incoherency (INCOH) is tier P1"),
        group("Multiple Excitation", [
            check("HOUSE.me", "Use Multiple Excitation"),
            sel("motion", "Input Motion Number", 1, 10),
            f("ME[motion].nfirst", "First Foundation Node", enable=["HOUSE.me", "==", 1]),
            f("ME[motion].nlast", "Last Foundation Node", enable=["HOUSE.me", "==", 1]),
            f("ME[motion].xc", "X Coord. of Control Point (unused)", enable=["HOUSE.me", "==", 1]),
            f("ME[motion].yc", "Y Coord. of Control Point (unused)", enable=["HOUSE.me", "==", 1]),
            f("ME[motion].zc", "Z Coord. of Control Point (unused)", enable=["HOUSE.me", "==", 1]),
            f("#AMP[motion]", "Spectral Amplification", "list", rows=3,
              tip="ratios of motion <no>, one per frequency of the set (AMP)"),
            check("HOUSE.cmplxspec", "Use Complex Spectral Amp."),
            greyed("Non-Uniform Motion"), greyed("Non-Uniform Soil")], col=2, note="ME and AMP are tier P1"),
        group("Wave Passage", [
            check("HOUSE.wpass", "Use Wave Passage",
                  force_on=[["WPASS.cohf", "in", [2, 3, 4, 5, 6, 7]], ["HOUSE.me", "==", 1]]),
            f("WPASS.appv", "Apparent Velocity for Line D"), f("WPASS.ang", "Angle Line D with X Axis"),
            f("WPASS.cohf", "Unlagged Coherency Model (1-7)")]),
        group("Motion Incoherency Simulation", [
            radio("INCOH.@stoch", "", [[0, "Deterministic (Median) Incoherency Input"],
                                       [1, "Stochastically Simulated Incoherency Input"]], enable=INCOH_ON),
            f("INCOH.hseed", "Horizontal SEED", enable=[INCOH_ON, ["INCOH.@stoch", "==", 1]]),
            f("INCOH.vseed", "Vertical SEED", enable=[INCOH_ON, ["INCOH.@stoch", "==", 1]]),
            f("INCOH.randphz", "Random Phase Angle", enable=[INCOH_ON, ["INCOH.@stoch", "==", 1]]),
            f("HOUSEX.nsim", "Number of Simulations", enable=[INCOH_ON, ["INCOH.@stoch", "==", 1]])], col=1),
        group("Superposition Mode", [
            radio("HOUSEX.supmode", "", [[0, "Linear"], [1, "Quadratic"]], enable=INCOH_ON),
            check("HOUSEX.ansys", "Ansys Model Input"),
            greyed("Ansys Model Type: Embedded / Surface", "Embedded",
                   "ANSYSMODELTYPE is tier P2 (enabled only with Ansys Model Input)")], col=2),
    ]
    tabs["FORCE"] = [
        group("FORCE Module Options", [
            radio("FORCE.opmode", "Operation Mode", SOL_DC),
            f("HOUSE.gravity", "Acceleration of Gravity"), f("SITE.fstep", "Frequency Step"),
            f("SITE.delt", "Time Step of Motion Control"), f("SITE.nft", "Nr. of Fourier Components"),
            f("SITE.freq", "Frequency Set Number")]),
    ]
    tabs["ANALYS"] = [
        group("Analysis", [
            radio("ANALYS.opmode", "Operation Mode", SOL_DC),
            radio("ANALYS.type", "Type of Analysis", TYPE_CH),
            radio("ANALYS.mode", "Mode of Analysis", [[0, "Initiation"], [1, "New Structure"],
                                                      [2, "New Seismic Environment"], [3, "New Dynamic Loading"]]),
            f("ANALYS.simul", "Simultaneous Cases"), check("ANALYS.save", "Save Restart Files"),
            check("ANALYSX.delrst", "Delete Restart Files"), check("ANALYS.prnt", "Print Amplitude Only")]),
        group("Frequency Numbers", [
            check("ANALYS.fopt", "Take Frequency Numbers from File1/File9"),
            f("SITE.freq", "Frequency Set Number")]),
        group("Control Motion Foundation Reference Point", [
            f("ANALYS.xc", "X-Coordinate of Control Point"), f("ANALYS.yc", "Y-Coordinate of Control Point"),
            f("ANALYS.zc", "Z-Coordinate of Control Point"), f("ANALYS.ang", "Coordinate Transformation Angle")],
            col=1),
        group("Incoherency", [
            radio("HOUSE.coh", "", [[0, "Coherent"], [1, "Incoherent"]],
                  disable_choice_if={"1": ["HOUSE.dim", "==", 1]}),
            check("HOUSE.wpass", "Wave Passage Effects Included",
                  force_on=[["WPASS.cohf", "in", [2, 3, 4, 5, 6, 7]], ["HOUSE.me", "==", 1]]),
            radio("ANALYSX.ffm", "", [[0, "Free-Field Load"], [1, "Free-Field Motion"]], enable=INCOH_ON)], col=1),
        group("Multiple Excitation", [
            check("HOUSE.me", "Use Multiple Excitation"), sel("motion", "Input Motion Number", 1, 10),
            f("ME[motion].nfirst", "First Foundation Node", enable=["HOUSE.me", "==", 1]),
            f("ME[motion].nlast", "Last Foundation Node", enable=["HOUSE.me", "==", 1]),
            f("ME[motion].xc", "X Coord. of Control Point", enable=["HOUSE.me", "==", 1]),
            f("ME[motion].yc", "Y Coord. of Control Point", enable=["HOUSE.me", "==", 1]),
            f("ME[motion].zc", "Z Coord. of Control Point", enable=["HOUSE.me", "==", 1])], col=2),
        group("Global Impedance Calculations", [
            radio("ANALYS.impe", "", [[0, "No Impedance Calculations"], [1, "Only Decoupled (Diagonal) Impedances"],
                                      [2, "Full Rigid Body Impedance Matrix 6X6"]])], col=2),
    ]
    tabs["MOTION"] = [
        group("Analysis", [
            radio("MOTION.opmode", "Operation Mode", SOL_DC), radio("ANALYS.type", "Type of Analysis", TYPE_CH),
            radio("MOTION.bl", "Baseline Correction", [[0, "No Correction"], [1, "With Correction"]])]),
        group("Response Spectrum Data", [
            f("MOTION.freq1", "First Frequency"), f("MOTION.freq2", "Last Frequency"),
            f("MOTION.fstep", "Total Number of Freq. Steps"),
            f("#DAMP", "Damping Ratios", "list", tip="fractions, e.g. 0.02 0.05 (at most 10)")]),
        group("Output Control", [
            check("MOTION.out", "Output Only Transfer Functions"), check("MOTION.cplx", "Save Complex Transfer Functions"),
            radio("MOTIONX.f1213", "Save FILE 12 or FILE 13", [[0, "None"], [1, "FILE13"], [2, "FILE12"]]),
            f("MOTION.dur", "Total Duration to be Plotted"),
            radio("MOTIONX.resp", "Response type for Foundation Vibration",
                  [[0, "Displacement"], [1, "Velocity"], [2, "Acceleration"]], visible=["ANALYS.type", "==", 1]),
            check("MOTIONX.srss", "Incoherent SSI: use SRSS"), button("Input (SRSSTF.txt)", "edit:SRSSTF.txt"),
            f("MOTION.interp", "Interpolation Option", "select", choices=INTERP_CH),
            radio("MOTION.pzadj", "Phase Adjustment", [[0, "0"], [1, "1"]]),
            f("MOTION.smo", "Smoothing Parameter", enable=["MOTION.interp", "!=", 6])], col=1),
        group("Nodal Output", [
            f("@NOUT", "Node List", "table",
              columns=[{"name": "dir", "label": "Direction", "w": "select", "choices": DIR_CH},
                       *[{"name": f"c{k}", "label": lab, "w": "check"} for k, lab in enumerate(NOUT_FLAGS, start=1)],
                       {"name": "nodes", "label": "Node list (1, 3-6 10)", "w": "text"}])], col=1, wide=True),
        group("Acceleration Time History Data", _acc_history(), col=2),
        group("Convert Time History to Response Spectrum", [
            check("MOTION.cnvrt", "Select External Files"),
            button("Input Time History Files (CONTTRS.txt)", "edit:CONTTRS.txt")], col=2),
        group("Post Processing Options", [
            check("MOTIONX.savetf", "Save TF in All Points"), check("MOTIONX.saveacc", "Save ACC in All Points"),
            check("MOTIONX.savers", "Save RS in All Points"), check("MOTIONX.saverot", "Save Rotation for Ansys"),
            check("MOTIONX.rsttf", "Restart for TF"), check("MOTIONX.rstacc", "Restart for ACC"),
            check("MOTIONX.rstrs", "Restart for RS")], col=2, note="post-processing flags are tier P1"),
        group("Binary Output Option", [greyed("Save Binary Database", "off", "BINOUT databases are tier P2")], col=2),
    ]
    tabs["STRESS"] = [
        group("Analysis", [
            radio("STRESS.opmode", "Operation Mode", SOL_DC), radio("ANALYS.type", "Type of Analysis", TYPE_CH)]),
        group("Output Control", [
            check("STRESS.iter", "Auto Computation of Strains in Soil El."),
            check("STRESS.save", "Save Stress Time Histories"), check("STRESS.itran", "Output Transfer Function"),
            radio("STRESSX.pzadj", "Phase Adjustment", [[0, "0"], [1, "1"]]),
            f("STRESS.interopt", "Interpolation Option", "select", choices=INTERP_CH),
            f("STRESSX.smo", "Smoothing Option", enable=["STRESS.interopt", "!=", 6]),
            f("STRESSX.skip", "Skip Time History Steps")]),
        group("Acceleration Time History Data", _acc_history(freq=True), col=1),
        group("Element Output Data", [
            f("@EOUT", "Element Output Data", "table", eout=True,
              columns=[{"name": "group", "label": "Group", "w": "int"},
                       {"name": "elements", "label": "Element List", "w": "text"},
                       {"name": "codes", "label": "Output Code (0 none, 1 max, 2 max + TH per component)",
                        "w": "codes"}])], col=2, wide=True),
        group("Post Processing Options", [
            check("STRESSX.savemax", "Save Max Value"), check("STRESSX.saveth", "Save Time History"),
            check("STRESSX.rstns", "Restart for Nodal Stress Contours"),
            check("STRESSX.rstsp", "Restart for Soil Pressure Contours"),
            button("Frame Selection (Frames.txt)", "edit:Frames.txt")], col=2, note="tier P1"),
        group("Section Cut Options", [greyed("Save Time History", "off", "SECDATAOPT is tier P1")], col=2),
        group("Binary Processing Option", [greyed("Save Binary Database", "off", "BINOUT databases are tier P2")],
              col=2),
    ]
    tabs["RELDISP"] = [
        group("Reference Location and Direction", [f("$RELFILE", "Complex TF File Name", "path")]),
        group("Output Control", [f("RELD.reldisoutput", "Save Rel Disp Complex TF", "bit", bit=1)]),
        group("Acceleration Time History Data", _acc_history(freq=True)),
        group("Nodal Output Data", [
            f("@RDND", "Nodal Output", "table",
              columns=[{"name": "node", "label": "Node Num.", "w": "int"},
                       *[{"name": lab.lower(), "label": lab, "w": "check"} for lab in ("X", "Y", "Z", "XX", "YY", "ZZ")]])],
            col=1, wide=True),
        group("Post Processing Options", [
            check("RELD.reldispsall", "Save Relative Displacement in All Nodes"),
            check("RELDX.saverot", "Save Rotations for ANSYS V11.0"),
            check("RELDX.rstframes", "Restart For Frame Generation")], col=1),
        group("Binary Disp. Option", [greyed("No Binary / TFD Binary / THD Binary", "No Binary",
                                             "BINOUT databases are tier P2")], col=1),
    ]
    tabs["NONLINEAR"] = [
        group("Global Modeling Options", [
            f("%EQL.disp", "Disp. Factor", "real", tip="Equivalent-Linear Displacement Factor EDF: x_eq = EDF x_max "
                                                       "(0.7-0.9, best about 0.8)"),
            f("%EQL.dampcutoff", "Damping Cutoff %", "real", tip="upper limit of the damping in percent (0 = none), "
                                                                 "e.g. 7 for the ASCE 4 level of cracked concrete"),
            f("%EQL.dampscale", "Damping Scale Factor", "real", tip="scales the damping (0 = 1)"),
            greyed("Material Parameter", "", "disabled in this version (manual)"),
            f("%EQL.nonlinopts", "Use Non-linear Panels", "bit", bit=1),
            f("%EQL.nonlinopts", "Use Non-linear Springs", "bit", bit=2),
            f("%EQL.nonlinopts", "Use Non-linear Beams", "bit", bit=4, disabled=True,
              tip="nonlinear beams are not available in this version"),
            check("%EQL.elasticd", "Include Elastic Damping",
                  tip="damping = hysteretic + elastic; unchecked: hysteretic only")],
            note="EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>. While <NonLinOpts> was never given, "
                 "the element types are those of the P and S records (shown); a click on a type stores it."),
        group("Backbone Curve Data", [
            sel("bbc", "Backbone Curve", 1, 99999),
            info("Defined curves", keys_of="%BBC"),
            f("%BBC[bbc].type", "Type", "select", choices=[[0, "0 (not set)"], [1, "1 CMS"], [2, "2 CMB"], [3, "3 TAK"],
                                                           [4, "4 GMR"]],
              tip="hysteretic model of the curve: used for the titles of the printout only"),
            f("%BBC[bbc].yield", "Yield Num.", "int", tip="1-based number of the yield point among the points"),
            button("New Curve", "newkey:bbc:%BBC", tip="the next free curve number (stored when a field is set)"),
            button("Delete Curve", "delentry:%BBC[bbc]", tip="DELBBC on OK")],
            col=1, note="BBC (file of X Y pairs), BBCGEN (22 points from the shear capacity of a panel) and BBCP also "
                        "define curves (Command Entry)."),
        group("Beam Data", [greyed("Beam"), greyed("Group Num."), greyed("Spring Gr."), greyed("BBC Num"),
                            greyed("Force Opt"), greyed("Beam End 1"), greyed("Beam End 2"),
                            info("Stored B records", context="beams")],
              col=1, note="NOT AVAILABLE IN THIS VERSION (B records are stored, not used)"),
        group("Backbone Curve Points", [
            f("%BBC[bbc].points", "", "table",
              columns=[{"name": "#", "label": "Point", "w": "index"}, {"name": "x", "label": "X", "w": "real"},
                       {"name": "y", "label": "Y", "w": "real"}],
              note="OK: DELBBC + BBCI, BBCX, BBCY of a new or changed curve")],
            col=2, note="Origin (0, 0) omitted; point 1 = cracking point. X: shear strain (Disp. Opt 1) or "
                        "displacement, Y: force (consistent units, never converted)."),
        group("Panel Data", [
            f("%P", "", "table", key="num",
              columns=[{"name": "num", "label": "Panel", "w": "int"},
                       {"name": "group", "label": "Group Num.", "w": "int"},
                       {"name": "bbc", "label": "BBC Num.", "w": "int"},
                       {"name": "disp", "label": "Disp. Opt.", "w": "select",
                        "choices": [[1, "1 shear strain"], [2, "2 bending (experimental)"]]},
                       {"name": "force", "label": "Force Opt.", "w": "select",
                        "choices": [[1, "1 CMS"], [2, "2 CMB (not included)"], [3, "3 TAK (experimental)"]]}],
              note="OK: PDEL for a removed row, P for a new or changed row")],
            full=True, note="P,<num>,<group>,<bbc>,<disp>,<force>: a SHELL group in a vertical plane with a material of "
                            "its own; PNLGEN makes one panel per vertical shell group, PLIST checks them."),
        group("Spring Data", [
            f("%S", "", "table", key="num",
              columns=[{"name": "num", "label": "Spring", "w": "int"},
                       {"name": "group", "label": "Group Num.", "w": "int"},
                       {"name": "elem", "label": "Elem Num.", "w": "int"},
                       {"name": "bbc", "label": "BBC Num.", "w": "int"},
                       {"name": "disp", "label": "Dof.", "w": "select", "choices": [[1, "1 X"], [2, "2 Y"], [3, "3 Z"]]},
                       {"name": "force", "label": "Force Opt.", "w": "select", "default": 4,
                        "choices": [[4, "4 GMR"], [1, "1 (not applicable)"], [2, "2 (not applicable)"],
                                    [3, "3 (not applicable)"]]}],
              note="OK: DELSPR for a removed row, S for a new or changed row")],
            full=True, note="S,<num>,<group>,<elem>,<bbc>,<disp>,<force>: an existing SPRING element (E command) with "
                            "an SC property of its own, acting on one translation."),
    ]
    aopt = record_class("AOPT")
    tabs["AFWRITE"] = [
        group("Modules processed by CHECK and AFWRITE (AOPT)", [
            check(f"AOPT.{fl.name}", fl.doc.split(" (")[0], disabled=fl.name in ("dep1", "dep2"))
            for fl in aopt.FIELDS]),
    ]
    notes = {"NONLINEAR": "Option NON (tier P2, Help > Option NON): the panel model is built with PNLGEN, WALLFLR, "
                          "PANELIZE ..., the curves with SHEAR / BBCGEN. After RELDISP, run Modules > NONLINEAR "
                          "(RUNNONLINEAR writes <model>.eql from this tab).",
             "EQUAKE": "", "SOIL": ""}
    return [{"name": t, "title": t, "groups": tabs[t], "note": notes.get(t, "")} for t in TABS]


# ======================================================================================
# Modules > ANSYS Eq. Static Load / ANSYS Dynamic Load (Option A, LOADGEN; spec 04 section 15.6)
# ======================================================================================
#: dialog -> (RUNLOADGEN analysis, the families it shows in commit order)
MODULE_DIALOGS = {"LOADGEN": ("STATIC", ("LOADGEN", "LGFILE", "LGNODE", "LGTIME", "LGMAP", "LGOPT")),
                  "LOADGENDYN": ("DYNAMIC", ("LOADGENDYN", "LGFILE", "LGNODE", "LGMAP", "LGOPT"))}
#: families of the Analysis Options window, in commit order
ANALYSIS_FAMILIES = ("NLSOIL", "NLSLAYER", "EQL", "BBC", "P", "S")


def _lg_common() -> List[Dict[str, Any]]:
    """Node lists, ANSYS numbering and output options (SASSI-EDU additions of both LOADGEN dialogs)."""
    return [
        f("%LGMAP.mode", "ANSYS Numbering (LGMAP)", "select",
          choices=[["IDENTITY", "IDENTITY (same numbers)"], ["PAIRS", "PAIRS (file)"], ["COORD", "COORD (positions)"]],
          tip="IDENTITY: the SASSI node numbers; PAIRS: a file of 'sassi_node ansys_node' lines; COORD: the ANSYS "
              "node at the same position in a .cdb / APDL file"),
        f("%LGMAP.file", "Pairs / ANSYS Model File", "path", enable=["%LGMAP.mode", "!=", "IDENTITY"]),
        f("%LGMAP.tol", "COORD Tolerance", "real", enable=["%LGMAP.mode", "==", "COORD"],
          placeholder="1e-6 of the model size"),
        f("%LGOPT.digits", "APDL Significant Digits", "int", tip="6..17 (default 12)"),
        check("%LGOPT.opmode", "Data Check Only (no file written)"),
        check("%LGOPT.rest", "Relative Displacements Start at Rest",
              tip="subtract the initial value of the relative displacements (default)")]


def loadgen_form(name: str) -> Dict[str, Any]:
    """"ANSYS Static Load Converter" (LOADGEN) / "ANSYS Dynamic Load Converter" (LOADGENDYN) of spec 04
    section 15.6, with the SASSI-EDU additions of docs/user/OPTION_A.md (history source, critical times,
    node lists, numbering, output options)."""
    if name == "LOADGEN":
        disp_on = ["%LOADGEN.data", "in", [1, 3, 4]]
        acc_on = ["%LOADGEN.data", "in", [2, 3]]
        groups = [
            group("Data to Add From ACS SASSI to the ANSYS model", [
                radio("%LOADGEN.data", "", [[1, "Displacement"], [2, "Acceleration"], [3, "Disp. and Accel."],
                                            [4, "Disp. for Soil Module"]]),
                check("%LOADGEN.multi", "Use Multiple File List Inputs", tip="one APDL file per critical time")]),
            group("SSI Model and Results Input", [
                f("%LGFILE.SSIPATH", "Path", "path", folder=True, placeholder="the model directory"),
                f("%LGFILE.HOUSE", "HOUSE Module Input", "path", placeholder="<model>.hou"),
                f("%LGFILE.DISP", "Displacement Results", "path", folder=True, placeholder="THD (frame folder)",
                  enable=disp_on),
                check("%LOADGEN.rotdisp", "Rotational Disp.", enable=disp_on),
                f("%LGFILE.DISPROT", "Rotational Displacement Results", "path", folder=True, placeholder="THDR",
                  enable=[disp_on, ["%LOADGEN.rotdisp", "==", 1]]),
                f("%LGFILE.ACC", "Trans. Acceleration Results", "path", folder=True, placeholder="ACC (frame folder)",
                  enable=acc_on),
                check("%LOADGEN.rotacc", "Rotational Accel.", enable=acc_on),
                f("%LGFILE.ACCROT", "Rotational Acceleration Results", "path", folder=True, placeholder="ACCR",
                  enable=[acc_on, ["%LOADGEN.rotacc", "==", 1]]),
                f("%LOADGEN.source", "History Source", "select",
                  choices=[["RESULTS", "RESULTS (MOTION / RELDISP files)"], ["FILE8", "FILE8 (computed)"]],
                  tip="RESULTS: the MOTION .ACC and RELDISP .THD files (or frames); FILE8: the same histories "
                      "computed from FILE8 and the control motion of the MOTION deck")]),
            group("Ansys Model and Data Input", [
                f("%LGFILE.ANSYSPATH", "Path", "path", folder=True, placeholder="the model directory")]),
            group("Mass Data for Internal Load (Ignore for Displacement)", [
                radio("%LOADGEN.masstype", "Mass Type", [[1, "Lumped Mass"], [2, "Master Node Mass"]], enable=acc_on),
                check("%LOADGEN.genmass", "Generate Mass Data", enable=acc_on,
                      tip="the mass file below is written (overwritten) from the HOUSE mass matrix"),
                f("%LGFILE.LUMPED", "For Lumped Mass: Lumped node", "path", placeholder="<model>.masl",
                  enable=[acc_on, ["%LOADGEN.masstype", "==", 1]]),
                f("%LGFILE.MASTER", "For Master Mass: Master Node Mass", "path", placeholder="<model>.masm",
                  enable=[acc_on, ["%LOADGEN.masstype", "==", 2]]),
                f("%LGNODE.M", "Master Nodes (LGNODE M)", "text", tip="1, 3-6 10",
                  enable=[acc_on, ["%LOADGEN.masstype", "==", 2]])], col=1),
            group("ANSYS Output File", [
                f("%LGFILE.APDL", "ADPL File", "path", placeholder="<model>_LGS.inp")], col=1),
            group("Critical Times (LGTIME)", [
                f("%LGTIME.crit", "Criterion", "select",
                  choices=[["V", "V: base shear (control dir.)"], ["VX", "VX: base shear X"],
                           ["VY", "VY: base shear Y"], ["VZ", "VZ: base shear Z"], ["MX", "MX: moment about X"],
                           ["MY", "MY: moment about Y"], ["MZ", "MZ: moment about Z"],
                           ["ACC", "ACC: node acceleration"], ["DISP", "DISP: node displacement"],
                           ["TIME", "TIME: given times"], ["STEP", "STEP: given time steps"]],
                  tip="the n largest peaks of the base shear along the control direction (V) or of a component, of "
                      "an overturning moment, of the absolute acceleration / relative displacement of a node; or "
                      "given times / time steps"),
                f("%LGTIME.n", "Number of Peaks", "int", enable=["%LGTIME.crit", "notin", ["TIME", "STEP"]]),
                f("%LGTIME.tsep", "Minimum Separation (s)", "real", enable=["%LGTIME.crit", "notin", ["TIME", "STEP"]]),
                f("%LGTIME.x0", "Moment Point X", "real", placeholder="centroid",
                  enable=["%LGTIME.crit", "in", ["MX", "MY", "MZ"]]),
                f("%LGTIME.y0", "Moment Point Y", "real", placeholder="centroid",
                  enable=["%LGTIME.crit", "in", ["MX", "MY", "MZ"]]),
                f("%LGTIME.z0", "Moment Point Z", "real", placeholder="centroid",
                  enable=["%LGTIME.crit", "in", ["MX", "MY", "MZ"]]),
                f("%LGTIME.node", "Node", "int", enable=["%LGTIME.crit", "in", ["ACC", "DISP"]]),
                f("%LGTIME.dof", "DOF (1-6)", "int", enable=["%LGTIME.crit", "in", ["ACC", "DISP"]]),
                f("%LGTIME.times", "Times (s) / Steps", "text", tip="t1 t2 ... or k1 k2 ... (1-based)",
                  enable=["%LGTIME.crit", "in", ["TIME", "STEP"]])], col=2),
            group("Nodes, Numbering and Output", [
                f("%LGNODE.D", "D Nodes (LGNODE D)", "text", placeholder="interaction nodes",
                  tip="the nodes that receive D: 1, 3-6 10")] + _lg_common(), col=2),
        ]
        return {"name": "LOADGEN", "title": "ANSYS Static Load Converter", "columns": 3, "width": "min(1120px, 96vw)",
                "run": {"label": "Run", "line": "RUNLOADGEN,STATIC"},
                "tabs": [{"name": "LOADGEN", "title": "ANSYS Static Load Converter", "groups": groups,
                          "note": "Option A, step 2 = equivalent static analysis in ANSYS (Help > Option A). Ok "
                                  "stores the settings (LOADGEN, LGFILE, LGNODE, LGTIME, LGMAP, LGOPT); Run also "
                                  "runs RUNLOADGEN,STATIC (needs HOUSE, MOTION and RELDISP results or FILE8)."}]}
    if name == "LOADGENDYN":
        groups = [
            group("SASSI Model and Results Input", [
                f("%LGFILE.SSIPATH", "Path", "path", folder=True, placeholder="the model directory"),
                f("%LGFILE.HOUSE", "HOUSE Module Input", "path", placeholder="<model>.hou"),
                f("%LGFILE.GROUND", "Ground Acceleration File", "path", placeholder="the control motion",
                  tip="ground or kinematic SSI acceleration, in g"),
                radio("%LOADGENDYN.gfopt", "Ground File Format", [[0, "dt, then one value per line"],
                                                                  [1, "(t, a) pairs"]]),
                f("%LOADGENDYN.gmult", "Ground File Scale Factor", "real")]),
            group("Ansys Model and Data Input", [
                f("%LGFILE.ANSYSPATH", "Path", "path", folder=True, placeholder="the model directory")]),
            group("Raleigh Damping Coeff.", [
                f("%LOADGENDYN.alpha", "Alpha", "real", enable=["%LOADGENDYN.zeta", "==", 0]),
                f("%LOADGENDYN.beta", "Beta", "real", enable=["%LOADGENDYN.zeta", "==", 0]),
                f("%LOADGENDYN.zeta", "Damping Ratio zeta (0 = Alpha, Beta)", "real",
                  tip="alpha and beta giving the damping ratio zeta at f1 and f2"),
                f("%LOADGENDYN.f1", "Frequency f1 (Hz)", "real", enable=["%LOADGENDYN.zeta", "!=", 0]),
                f("%LOADGENDYN.f2", "Frequency f2 (Hz)", "real", enable=["%LOADGENDYN.zeta", "!=", 0])], col=1,
                note="C = alpha M + beta K. WARNING (manual): Rayleigh damping is a significant limitation of the "
                     "dynamic stress analysis."),
            group("ANSYS Output File", [
                f("%LGFILE.APDLDYN", "ADPL File", "path", placeholder="<model>_LGD.inp")], col=1),
            group("Method", [
                radio("%LOADGENDYN.method", "", [["REL", "REL: ACEL = ground motion, D = interface displacements "
                                                         "relative to it (manual)"],
                                                 ["ACC", "ACC: fixed interface, ACEL = absolute acceleration of the "
                                                         "reference node"]]),
                f("%LOADGENDYN.refnode", "Reference Node (0 = control motion)", "int"),
                f("%LOADGENDYN.source", "History Source", "select",
                  choices=[["RESULTS", "RESULTS (RELDISP / MOTION files)"], ["FILE8", "FILE8 (computed)"]],
                  tip="RESULTS: the RELDISP .THD and MOTION files (or frames); FILE8: the same histories computed "
                      "from FILE8 and the control motion of the MOTION deck"),
                check("%LOADGENDYN.rotdisp", "Rotations Prescribed (D)"),
                check("%LOADGENDYN.rotacc", "Rotational Accelerations in the Check Tables")], col=2),
            group("Nodes, Numbering and Output", [
                f("%LGNODE.D", "D Nodes (LGNODE D)", "text", placeholder="interaction nodes",
                  tip="the nodes that receive D: 1, 3-6 10"),
                f("%LGNODE.A", "Check Nodes (LGNODE A)", "text",
                  tip="nodes whose SSI absolute accelerations are exported as tables: 1, 3-6 10")] + _lg_common(), col=2),
        ]
        return {"name": "LOADGENDYN", "title": "ANSYS Dynamic Load Converter", "columns": 3,
                "width": "min(1120px, 96vw)", "run": {"label": "Run", "line": "RUNLOADGEN,DYNAMIC"},
                "tabs": [{"name": "LOADGENDYN", "title": "ANSYS Dynamic Load Converter", "groups": groups,
                          "note": "Option A, step 2 = transient analysis in ANSYS (Help > Option A). Ok stores the "
                                  "settings (LOADGENDYN, LGFILE, LGNODE, LGMAP, LGOPT); Run also runs "
                                  "RUNLOADGEN,DYNAMIC (needs the RELDISP relative displacements or FILE8)."}]}
    raise KeyError(f"unknown dialog {name}")


def form(name: str) -> Dict[str, Any]:
    """Layout of dialog ``name``: ANALYSIS, MODEL, WRITE, CHECK, LOADGEN or LOADGENDYN."""
    name = name.upper()
    if name == "ANALYSIS":
        return {"name": "ANALYSIS", "title": "Analysis Options", "tabs": analysis_tabs()}
    if name in MODULE_DIALOGS:
        return loadgen_form(name)
    if name == "MODEL":
        return {"name": "MODEL", "title": "Model Options", "tabs": [{"name": "MODEL", "title": "Model Options", "groups": [
            group("Incompatible Modes", [radio("MOPT.incomp", "", [[0, "Include Incompatible Modes"],
                                                                   [1, "Suppress Incompatible Modes"]])]),
            group("General Matrix", [radio("MOPT.matrix", "", [[0, "Mass Matrix"], [1, "Weight Matrix"]])]),
            group("", [check("MOPT.mass", "Overwrite Mass"), check("MOPT.force", "Overwrite Force")])]}]}
    if name == "WRITE":
        return {"name": "WRITE", "title": "Extended Write Options", "tabs": [{"name": "WRITE", "title": "Extended Write Options", "groups": [
            group("", [check("WRITE.mdl", "MDL command in *.Pre"),
                       check("WRITE.ext", "Extend integer fields (for Models with more than 100000 nodes)",
                             tip="accepted and recorded; no effect on SASSI-EDU decks (D-AFW-03)"),
                       check("WRITE.afwr", "AFWR Command in *.Pre"),
                       check("WRITE.sim", "Simulation Commands")]),
            group("Simulation Command Location", [
                radio("WRITE.sim_location", "", [["sim", "Write to *-Sim.Pre"], ["pre", "Write to *.Pre"]],
                      enable=["WRITE.sim", "==", 1])])]}]}
    if name == "CHECK":
        return {"name": "CHECK", "title": "Check Options", "tabs": [{"name": "CHECK", "title": "Check Options", "groups": [
            group("", [check("CHECK.show_warnings", "Show Warnings"), check("CHECK.show_errors", "Show Errors"),
                       check("CHECK.suppress_window", "Suppress Error Window"),
                       f("CHECK.break_at", "Break Check at [N] Messages", "int")],
                  note="session settings, not saved (UI-07)")]}]}
    raise KeyError(f"unknown dialog {name}")


# ======================================================================================
# Current values
# ======================================================================================
def _typed_value(rec, fl) -> Any:
    """Dialog value of one record field: typed value, or '' for a blank field with a computed default."""
    if (rec.command, fl.name) in BLANK_DEFAULT and rec.arg(rec.field_index(fl.name)) is None:
        return ""
    return rec.get(fl.name)


def record_values(rec) -> Dict[str, Any]:
    return {fl.name: _typed_value(rec, fl) for fl in rec.FIELDS}


def _record_names() -> List[str]:
    return [s.command for s in OPTION_SPECS.values() if s.kind == "record" and s.record is not None]


def _indexed_names() -> List[str]:
    return [s.command for s in OPTION_SPECS.values() if s.kind == "indexed" and s.record is not None]


def _key_str(k) -> str:
    return "|".join(str(x) for x in k) if isinstance(k, tuple) else str(k)


def _numbers_text(vals: Sequence[float]) -> str:
    return " ".join(fmt_num(v) for v in vals)


def main_wave(wopt: int) -> int:
    """The implicit wave field of a model without WAVE entries: vertical SV (SH for ``<wopt>`` = 1)."""
    return 4 if int(wopt) == 1 else 2


def wave_effective_default(model, key: int, wopt: Optional[int] = None) -> Dict[str, Any]:
    """WAVE page of a model without that WAVE entry -- what AFWRITE writes for it
    (:meth:`sassi.prep.check.Checker.effective_waves`):

    * no WAVE entry stored at all: the new-model field, vertical SV (SH for ``<wopt>`` = 1), so that
      page shows "field on" and the others "no field";
    * some WAVE entry stored: only the stored entries are written, so every other page is "no field".

    ``wopt`` overrides the stored SITE ``<wopt>`` (the dialog may be switching it).
    """
    if wopt is None:
        wopt = get_record(model, "SITE").wopt
    on = not get_entries(model, "WAVE") and int(key) == main_wave(wopt)
    return {"type": int(key), "opt": 1 if on else 0, "ratio1": 1.0, "ratio2": 1.0, "angle": 0.0}


def request_default(name: str) -> Dict[str, Any]:
    """Record defaults of an indexed setter (for a SOIL layer request: the 5.4 values preselected
    when output is first requested at a layer)."""
    return record_values(record_class(name)(name))


def absent_entry(model, name: str, key, wopt: Optional[int] = None) -> Dict[str, Any]:
    """Values shown for an indexed entry that is not stored: what the analysis uses for it.

    WAVE pages follow :func:`wave_effective_default`; SOIL layer requests (LAYER_REQUESTS) are "not
    requested" (every flag 0: AFWRITE writes no row, SOIL computes nothing at that layer); other
    entries show the record defaults.  The key field is set.
    """
    if name == "WAVE":
        return wave_effective_default(model, int(key), wopt)
    if name == "SPRO" and key is not None:
        dflt = spro_defaults(model).get(str(key))
        if dflt is not None:
            return dict(dflt)
    cls = record_class(name)
    vals = request_default(name)
    if name in LAYER_REQUESTS:
        for fl in cls.FIELDS:
            if fl.name not in cls.KEY:
                vals[fl.name] = "" if fl.type is str else fl.type(0)
    if cls.KEY:
        vals[cls.KEY[0]] = key
    return vals


def spro_defaults(model) -> Dict[str, Dict[str, Any]]:
    """The default SOIL profile of a model without SPRO entries, by layer key (what AFWRITE writes, D-W5-10);
    empty when SPRO entries are stored, no default applies or EDUOPT,DEFAULTS,OFF."""
    prof, use, _ = DEF.soil_profile(model)
    if use is None:
        return {}
    return {str(k): {"layer": int(k), "prop": int(rec.prop), "dynprop": rec.dynprop or ""} for k, rec in prof}


def _spro_short(dflt: Dict[str, Dict[str, Any]]) -> str:
    """One line for the SOIL tab: "default: Sand 1-4, Rock 5; 6 = half-space L 3 (EDU-29)"."""
    keys = sorted(int(k) for k in dflt)
    if not keys:
        return ""
    last = keys[-1]
    groups = []
    for lab in ("Sand", "Rock"):
        ks = [k for k in keys[:-1] if dflt[str(k)]["dynprop"] == lab]
        if ks:
            groups.append(f"{lab} {DEF._ranges(ks)}")
    return f"default: {', '.join(groups)}; {last} = half-space L {dflt[str(last)]['prop']} (EDU-29)"


def input_defaults_context(model) -> Optional[Dict[str, Any]]:
    """Placeholder texts of the blank input files with a built-in default (requirements 7.19); the browser
    applies the same conditions as :mod:`sassi.prep.defaults` to the values being edited
    (``OptionsForm.defaultPlaceholder``).  None with EDUOPT,DEFAULTS,OFF (no defaults, no placeholders)."""
    if not DEF.enabled(model):
        return None
    rec = LIB.entry(DEF.SEISMIC_RECORD)
    return {
        "dt": rec.dt if rec is not None else 0.005,
        "thfile": {"0": f"built-in: RG 1.60, 0.30 g ({DEF.SEISMIC_RECORD})",
                   "1": f"built-in: 5 Hz Ricker pulse ({DEF.VIBRATION_LOAD})"},
        "soil_thfile": f"built-in: RG 1.60, 0.30 g ({DEF.SEISMIC_RECORD})",
        "no_default": "no built-in default: it needs pairs off (fopt 0) and time step 0.005 s",
        "no_default_soil": "no built-in default: it needs time step 0.005 s",
        "rsin": f"built-in: RG 1.60 H, 0.30 g ({DEF.TARGET_SPECTRUM})",
        "rsout": DEF.default_output_name(model, "RSOUT", 0).replace("_eq0.", "_eq{i}."),
        "accout": DEF.default_output_name(model, "ACCOUT", 0).replace("_eq0.", "_eq{i}."),
    }


def stochastic_incoherency(hseed, vseed, randphz) -> bool:
    """D-INC-02: the incoherency input is stochastically simulated iff (HSeed != 0 or VSeed != 0) and
    RandPhz > 0 (random phases uniform in [-RandPhz, +RandPhz]); otherwise deterministic (median)."""
    return (int(hseed) != 0 or int(vseed) != 0) and float(randphz) > 0


def values(interp, name: str = "ANALYSIS") -> Dict[str, Any]:
    """Current values of dialog ``name`` for the active model (``GET /api/options/<name>``)."""
    model = interp.model
    name = name.upper()
    if name == "MODEL":
        rec = model.mopt
        return {"records": {"MOPT": {fl.name: rec.get(fl.name) for fl in rec.FIELDS}}}
    if name == "WRITE":
        wo = interp.write_options
        return {"records": {"WRITE": {"mdl": int(bool(wo.get("mdl"))), "ext": int(bool(wo.get("ext"))),
                                      "afwr": int(bool(wo.get("afwr"))), "sim": int(bool(wo.get("sim"))),
                                      "sim_location": str(wo.get("sim_location", "sim"))}}}
    if name == "CHECK":
        from ..prep.commands.checks import check_options
        co = check_options(interp)
        return {"records": {"CHECK": {"show_warnings": int(co.show_warnings), "show_errors": int(co.show_errors),
                                      "suppress_window": int(co.suppress_window), "break_at": co.break_at}}}
    if name in MODULE_DIALOGS:
        out = CR.values(model, MODULE_DIALOGS[name][1])
        out["context"] = {"model": model.name, "path": model.path, "analysis": MODULE_DIALOGS[name][0]}
        return out
    records = {}
    defaults = {}
    for cmd in _record_names():
        records[cmd] = record_values(get_record(model, cmd))
        defaults[cmd] = record_values(record_class(cmd)(cmd))
    inc = records.get("INCOH")
    if inc is not None:
        inc["@stoch"] = int(stochastic_incoherency(inc["hseed"], inc["vseed"], inc["randphz"]))
    indexed: Dict[str, Dict[str, Any]] = {}
    for cmd in _indexed_names():
        if cmd in ("EDUOPT", "SYMM", "DYNP"):
            continue
        indexed[cmd] = {_key_str(k): record_values(r) for k, r in get_entries(model, cmd)}
    strings = {n: model.options.string(n) for n in ("THFILE", "THTIT", "EQTIT", "RELFILE")}
    lists = {"DAMP": _numbers_text(model.damp), "TOPL": " ".join(str(v) for v in model.topl)}
    for no, vals in sorted(model.amp.items()):
        lists[f"AMP[{no}]"] = _numbers_text(vals)
    requests = {
        "NOUT": [dict({"dir": r.dir, "nodes": _ids_text(r.nodes)},
                      **{f"c{k}": int(r.codes[k - 1]) if k <= len(r.codes) else 0 for k in range(1, 7)})
                 for r in model.nout],
        "EOUT": [{"group": r.group, "elements": _ids_text(r.elements), "codes": list(r.codes)} for r in model.eout],
        "RDND": [dict({"node": r.node}, **{lab: int(r.flags[k] >= 1) if k < len(r.flags) else 0
                                            for k, lab in enumerate(("x", "y", "z", "xx", "yy", "zz"))})
                 for r in model.rdnd],
        "CORR": [record_values(r) for _, r in get_entries(model, "CORR")],
    }
    dynp = sorted({k[0] for k, _ in model.options.entries("DYNP") if isinstance(k, tuple)})
    groups = {str(g.id): {"type": ELEMENT_TYPE_NAMES.get(g.type, str(g.type)),
                          "components": ELEMENT_COMPONENTS.get(ELEMENT_TYPE_NAMES.get(g.type, ""), []),
                          "n": len(g.elements)} for g in model.groups.values()}
    wave_defaults = {str(t): wave_effective_default(model, t) for t in WAVE_NAMES}
    # the dialog may switch <wopt> before OK: the absent pages for either value
    wave_defaults_by_wopt = {str(w): {str(t): wave_effective_default(model, t, w) for t in WAVE_NAMES}
                             for w in WAVE_PAGES}
    # absent entries as the analysis sees them (shown), and the values preselected on a first edit
    # (WAVE pages are shown from wave_defaults_by_wopt; its entry here only lists the fields)
    indexed_defaults = {cmd: absent_entry(model, cmd, None) if cmd != "WAVE" else request_default(cmd)
                        for cmd in indexed}
    for cmd in indexed_defaults:
        cls = record_class(cmd)
        if cls.KEY:
            indexed_defaults[cmd].pop(cls.KEY[0], None)
    request_defaults = {cmd: request_default(cmd) for cmd in LAYER_REQUESTS if cmd in indexed}
    # entries shown per key: the default SOIL profile while no SPRO is stored (D-W5-10)
    spro_dflt = spro_defaults(model)
    indexed_defaults_by_key = {"SPRO": spro_dflt} if spro_dflt else {}
    prof_use = DEF.soil_profile(model)[1]
    lib_dynp = [lab for lab in LIB.dynp_labels() if lab not in dynp] if DEF.enabled(model) else []
    out = {"records": records, "defaults": defaults, "indexed": indexed, "strings": strings, "lists": lists,
           "requests": requests, "wave_defaults": wave_defaults, "wave_defaults_by_wopt": wave_defaults_by_wopt,
           "indexed_defaults": indexed_defaults, "request_defaults": request_defaults,
           "indexed_defaults_by_key": indexed_defaults_by_key,
           "context": {"dynp": dynp, "groups": groups, "model": model.name, "path": model.path,
                       "layers": sorted(model.layers), "topl": list(model.topl),
                       "blank_default": {f"{a}.{b}": v for (a, b), v in BLANK_DEFAULT.items()},
                       # built-in inputs (requirements 7.19): the Library picker, the library soil curves the
                       # model does not define, the placeholders of blank files and the default profile
                       "library": LIB.catalogue_json(), "dynp_library": lib_dynp,
                       "input_defaults": input_defaults_context(model),
                       "spro_default": (_spro_short(spro_dflt) if prof_use is not None else
                                        "none (SPRO entries are stored)" if get_entries(model, "SPRO") else
                                        "none (needs TOPL, SITE <hs>, EDUOPT,DEFAULTS,ON)")}}
    # command records outside OPTION_SPECS: SOIL-NON (SOIL tab) and Option NON (NONLINEAR tab)
    out.update(CR.values(model, ANALYSIS_FAMILIES))
    out["context"].update(CR.context(model))
    return out


def _ids_text(ids: Sequence[int]) -> str:
    """Compact id list ``1-5 7 9-12`` (the dialog's list syntax, L8)."""
    ids = sorted(set(int(v) for v in ids))
    out: List[str] = []
    i = 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        out.append(str(ids[i]) if j == i else f"{ids[i]}-{ids[j]}")
        i = j + 1
    return " ".join(out)


# ======================================================================================
# Commit: values -> validated command text (L17, UI-06)
# ======================================================================================
class DialogError(ValueError):
    """A dialog commit refused by validation; ``problems`` holds the Chapter 10 messages."""

    def __init__(self, problems: List[str]):
        super().__init__("; ".join(problems))
        self.problems = problems


def problem_message(p) -> str:
    """``Error 47 : <manual title>  [detail]`` for a record problem (CHECK catalogue text)."""
    kind, num, detail = p[0], p[1], p[2]
    if isinstance(num, int):
        title = (ERRORS if kind == "Error" else WARNINGS).get(num, "")
    else:
        title = EDU_TEXT.get(num, ("", ""))[1]
    return f"{kind} {num} : {format_title(title, problem_fmt(p))}" + (f"  [{detail}]" if detail else "")


def _token(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, (int, float)):
        return fmt_num(v)
    return str(v).strip()


def _same_typed(a, b, fields) -> bool:
    for fl in fields:
        va, vb = a.get(fl.name), b.get(fl.name)
        if isinstance(va, float) or isinstance(vb, float):
            try:
                if float(va) != float(vb):
                    return False
            except (TypeError, ValueError):
                return False
        elif va != vb:
            return False
    return True


def _same_record(a, b) -> bool:
    """Same typed values *and* the same blank/given state of the 'computed default' fields
    (BLANK_DEFAULT): blanking an explicit SITE <freq2> = 2048 changes the deck (NFFT/2) although
    both read 2048 through the record default."""
    if not _same_typed(a, b, a.FIELDS):
        return False
    for (cmd, fname) in BLANK_DEFAULT:
        if a.command == cmd:
            k = a.field_index(fname)
            if a.given(k) != b.given(k):
                return False
    return True


def _new_record(name: str, base, changes: Dict[str, Any], errors: List[str]):
    """Record ``name`` = ``base`` tokens with the dialog ``changes`` applied (blank -> default)."""
    cls = record_class(name)
    toks = base.to_tokens() if base is not None else []
    toks = toks + [""] * (len(cls.FIELDS) - len(toks))
    for fname, v in changes.items():
        if fname.startswith("@"):
            continue                          # dialog-only pseudo fields (INCOH.@stoch)
        try:
            k = cls.field_index(fname)
        except KeyError:
            errors.append(f"{name}: unknown field {fname}")
            continue
        toks[k - 1] = _token(v)
    # a dialog writes every field (manual-identical lines): blank fields get their default text,
    # except the fields whose blank means "computed" (BLANK_DEFAULT) and text fields
    for i, fl in enumerate(cls.FIELDS):
        if not toks[i].strip() and fl.type is not str and (name, fl.name) not in BLANK_DEFAULT:
            toks[i] = _token(fl.default)
    rec = cls(name, toks)
    for msg in rec.token_errors():
        errors.append(f"{name} {msg}")
    return rec


def _record_line(rec) -> str:
    return to_command(rec)


_SPLIT = re.compile(r"[\s,;]+")


def parse_numbers(text: str, integer: bool = False) -> List[float]:
    """List-field text (separators blank, tab, ',', ';', Enter; spec 05a section 7.2)."""
    out = []
    for t in _SPLIT.split(str(text or "").strip()):
        if not t:
            continue
        try:
            out.append(parse_int(t)[0] if integer else parse_float(t))
        except NumberError:
            raise ValueError(f"'{t}' is not a number") from None
    return out


def parse_ids(text: str) -> List[int]:
    """Node / element list ``1, 3-6 10`` -> ids (L8; descending ranges rejected)."""
    out: List[int] = []
    s = re.sub(r"\s*-\s*", "-", str(text or "").strip())
    for t in _SPLIT.split(s):
        if not t:
            continue
        m = re.fullmatch(r"(\d+)-(\d+)", t)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            if b < a:
                raise ValueError(f"descending range {t}")
            out.extend(range(a, b + 1))
        elif t.isdigit():
            out.append(int(t))
        else:
            raise ValueError(f"'{t}' is not a node/element number or range")
    return out


def _chunks(name: str, vals: Sequence[Any], n: int, prefix: Sequence[Any] = ()) -> List[str]:
    lines = []
    for i in range(0, len(vals), n):
        lines.append(join_command(name, list(prefix) + [fmt_num(v) if isinstance(v, (int, float)) else v
                                                         for v in vals[i:i + n]]))
    return lines


def commit_commands(interp, payload: Dict[str, Any], name: str = "ANALYSIS") -> Tuple[List[str], List[str]]:
    """Command lines of a dialog commit and the informational notes (forced values ...).

    ``payload`` has the structure returned by :func:`values` (only the parts the dialog changed are
    needed).  Raises :class:`DialogError` with the Chapter 10 messages when a value is invalid.
    """
    name = name.upper()
    model = interp.model
    errors: List[str] = []
    notes: List[str] = []
    lines: List[str] = []
    recs_in: Dict[str, Dict[str, Any]] = {k.upper(): dict(v) for k, v in (payload.get("records") or {}).items()}
    if name == "MODEL":
        from ..model.options import MoptRecord
        cur = model.mopt
        vals = {fl.name: cur.get(fl.name) for fl in MoptRecord.FIELDS}
        for k, v in recs_in.get("MOPT", {}).items():
            if k not in vals:
                errors.append(f"MOPT: unknown field {k}")
                continue
            try:
                vals[k] = parse_int(_token(v))[0]
            except NumberError:
                errors.append(f"MOPT {k}: '{v}' is not an integer")
                continue
            if vals[k] not in (0, 1):
                errors.append(f"MOPT {k}: {vals[k]} (0 or 1)")
        if errors:
            raise DialogError(errors)
        if any(vals[k] != cur.get(k) for k in vals):
            lines.append(join_command("MOPT", [vals[fl.name] for fl in MoptRecord.FIELDS]))
        return lines, notes
    if name in MODULE_DIALOGS:                      # ANSYS Eq. Static Load / Dynamic Load (Option A)
        for key in ("records", "indexed", "strings", "lists", "requests"):
            if payload.get(key):
                errors.append(f"the {name} dialog has no {key} (its fields are LOADGEN records)")
        if not errors:
            lines = CR.commit(interp, payload, MODULE_DIALOGS[name][1], errors)
        if errors:
            raise DialogError(errors)
        return lines, notes

    # ---- forced values (D-UI-04)
    house = recs_in.setdefault("HOUSE", {})
    cur_house = get_record(model, "HOUSE")
    cur_wpass = get_record(model, "WPASS")

    def eff(rec_name, fname, cur):
        v = recs_in.get(rec_name, {}).get(fname, cur.get(fname))
        try:
            return parse_int(_token(v))[0]
        except NumberError:
            return v
    if eff("HOUSE", "dim", cur_house) == 1 and eff("HOUSE", "coh", cur_house) == 1:
        house["coh"] = 0
        notes.append("2D analysis: incoherent motion is not available -- Coherent used (D-UI-04)")
    if eff("HOUSE", "wpass", cur_house) == 0:
        why = []
        if eff("WPASS", "cohf", cur_wpass) in (2, 3, 4, 5, 6, 7):
            why.append(f"unlagged coherency model {eff('WPASS', 'cohf', cur_wpass)}")
        if eff("HOUSE", "me", cur_house) == 1:
            why.append("multiple excitation")
        if why:
            house["wpass"] = 1
            notes.append("Use Wave Passage switched on: required by " + " and ".join(why) + " (D-UI-04)")
    if not house:
        recs_in.pop("HOUSE")
    inc = recs_in.get("INCOH")
    if inc is not None and "@stoch" in inc:
        _incoherency_input(get_record(model, "INCOH"), inc, errors, notes)

    # ---- record setters, in the order of the option specifications
    order = [c for c in _record_names() if c in recs_in]
    for cmd in recs_in:
        if cmd not in order:
            errors.append(f"{cmd} is not an option record of this dialog")
    for cmd in order:
        stored = typed(model.options.record(cmd))
        cur = get_record(model, cmd)
        rec = _new_record(cmd, stored, recs_in[cmd], errors)
        if errors:
            continue
        if _same_record(rec, cur):          # unchanged (a never-stored record: still all defaults)
            continue
        for p in rec.problems():
            (errors if p[0] == "Error" else notes).append(f"{cmd}: {problem_message(p)}")
        lines.append(_record_line(rec))

    # ---- indexed setters
    wopt = eff("SITE", "wopt", get_record(model, "SITE"))     # the <wopt> this commit leaves
    if wopt not in (0, 1):
        wopt = int(get_record(model, "SITE").wopt)
    for cmd, entries in (payload.get("indexed") or {}).items():
        cmd = cmd.upper()
        if cmd not in _indexed_names():
            errors.append(f"{cmd} is not an indexed option of this dialog")
            continue
        cls = record_class(cmd)
        stored_map = {_key_str(k): r for k, r in get_entries(model, cmd)}
        emitted: Dict[str, Any] = {}

        def entry_record(key, vals):
            vals = dict(vals)
            kname = cls.KEY[0] if cls.KEY else None
            if kname:
                vals[kname] = key
            stored = stored_map.get(str(key))
            return _new_record(cmd, stored if stored is not None else cls(cmd, [str(key)]), vals, errors)

        for key in sorted(entries, key=lambda s: (len(str(s)), str(s))):
            stored = stored_map.get(str(key))
            rec = entry_record(key, entries[key])
            if errors:
                continue
            if stored is not None:
                if _same_typed(rec, stored, rec.FIELDS):
                    continue
            else:
                ref = absent_entry(model, cmd, parse_int(str(key))[0] if str(key).lstrip("-").isdigit() else key,
                                   wopt)
                if all(_eq(rec.get(fl.name), ref.get(fl.name)) for fl in rec.FIELDS):
                    continue
            if cmd in ("RSIN", "RSOUT", "ACCIN", "ACCOUT", "TPSD") and not rec.file and stored is None:
                continue
            emitted[str(key)] = rec
        if cmd == "WAVE" and emitted and not stored_map:
            # No WAVE stored yet: AFWRITE writes the implicit vertical SV (SH) field only while no WAVE
            # entry exists (Checker.effective_waves).  The first stored entry would silently drop it,
            # so the field the dialog shows on that page is stored with it.
            main = str(main_wave(wopt))
            if main not in emitted:
                emitted[main] = entry_record(main, entries.get(main) or wave_effective_default(model, int(main), wopt))
                notes.append(f"WAVE: the {WAVE_NAMES[int(main)]}-wave field shown on its page (new-model default) is "
                             f"stored with the other wave pages -- once a WAVE entry exists AFWRITE writes only "
                             f"the stored ones")
        if cmd == "SPRO" and emitted and not stored_map:
            # No SPRO stored yet: AFWRITE writes the default profile (D-W5-10) only while no SPRO entry exists;
            # the first stored entry would silently drop it, so the layers the dialog shows are stored with it.
            dflt = spro_defaults(model)
            added = [k for k in dflt if k not in emitted]
            for k in added:
                emitted[k] = entry_record(k, entries.get(k) or dflt[k])
            if added:
                notes.append("SPRO: the default SOIL profile shown on the layer pages (built-in curves, EDU-29) is "
                             "stored with the layer you changed -- once a SPRO entry exists AFWRITE writes only the "
                             "stored ones")
        if errors:
            continue
        for key in sorted(emitted, key=lambda s: (len(s), s)):
            rec = emitted[key]
            for p in rec.problems():
                (errors if p[0] == "Error" else notes).append(f"{cmd} {key}: {problem_message(p)}")
            lines.append(_record_line(rec))

    # ---- string setters
    for sname, v in (payload.get("strings") or {}).items():
        sname = sname.upper()
        if sname not in ("THFILE", "THTIT", "EQTIT", "RELFILE"):
            errors.append(f"{sname} is not a string option of this dialog")
            continue
        v = "" if v is None else str(v).strip()
        if v != model.options.string(sname):
            lines.append(join_command(sname, [v], text_last=True) if v else f"{sname},")

    # ---- list commands
    for lname, text in (payload.get("lists") or {}).items():
        lname = lname.upper()
        try:
            if lname == "DAMP":
                vals = parse_numbers(text)
                if any(not 0.0 < v < 1.0 for v in vals):
                    errors.append(problem_message(("Error", 72, "damping ratios must lie in (0, 1)")))
                    continue
                if len(vals) > 10:
                    errors.append("DAMP: at most 10 damping ratios")
                    continue
                if [float(v) for v in vals] != [float(v) for v in model.damp]:
                    lines.append("DAMP,0")
                    lines += _chunks("DAMP", vals, 10)
            elif lname == "TOPL":
                vals = [int(v) for v in parse_numbers(text, integer=True)]
                if any(v <= 0 for v in vals):
                    errors.append("TOPL: layer numbers must be positive")
                    continue
                if len(vals) > 200:
                    errors.append(problem_message(("Warning", 8, f"{len(vals)} top layers (at most 200)")))
                    continue
                if vals != list(model.topl):
                    lines.append("TOPL,0")
                    lines += _chunks("TOPL", vals, 20)
            elif lname.startswith("AMP[") and lname.endswith("]"):
                no = int(lname[4:-1])
                vals = parse_numbers(text)
                if [float(v) for v in vals] != [float(v) for v in model.amp.get(no, [])]:
                    lines.append(f"AMP,{no},0")
                    if vals:
                        if vals[0] == 0:
                            errors.append("AMP: the first ratio of a motion cannot be 0 (0 clears the list)")
                            continue
                        lines += _chunks("AMP", vals, 20, prefix=[no])
            else:
                errors.append(f"{lname} is not a list option of this dialog")
        except ValueError as exc:
            errors.append(f"{lname}: {exc}")

    # ---- request tables
    for rname, rows in (payload.get("requests") or {}).items():
        rname = rname.upper()
        try:
            new_lines = _request_lines(model, rname, rows or [])
        except ValueError as exc:
            errors.append(f"{rname}: {exc}")
            continue
        lines += new_lines

    # ---- command records outside OPTION_SPECS (NLSOIL, NLSLAYER, EQL, BBC, P, S): validated by their commands
    lines += CR.commit(interp, payload, ANALYSIS_FAMILIES, errors)
    if errors:
        raise DialogError(errors)
    return lines, notes


def _eq(a, b) -> bool:
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return (a or "") == (b or "")


def _incoherency_input(cur, inc: Dict[str, Any], errors: List[str], notes: List[str]) -> None:
    """HOUSE tab "Motion Incoherency Simulation" radio (``INCOH.@stoch``, a dialog-only field) applied
    to the INCOH values of a commit (requirements 5.4; 05b section 3; D-INC-02).

    * Deterministic (Median) chosen while the stored input is stochastic: HSeed = VSeed = RandPhz = 0.
    * Stochastically Simulated: a Random Phase Angle of 0 becomes the 5.4 default 180 degrees, and
      the commit is refused unless HSeed or VSeed is non-zero and RandPhz > 0 -- otherwise the input
      would be stored (and analysed, D-INC-02) as deterministic while the dialog showed "stochastic".
    """
    try:
        st = _token(inc.pop("@stoch"))
        stoch = bool(parse_int(st)[0]) if st else False
    except NumberError:
        errors.append("INCOH: the incoherency input choice must be 0 (deterministic) or 1 (stochastic)")
        return
    was = stochastic_incoherency(cur.hseed, cur.vseed, cur.randphz)
    if not stoch:
        # deterministic (median) input: zero phases -- unless the stored values are left untouched
        # (a mixed state typed as a command is kept: posting the dialog back changes nothing, L17)
        touched = any(name in inc and not _eq(inc[name], cur.get(name)) for name in ("hseed", "vseed", "randphz"))
        if was or touched:
            inc.update(hseed=0, vseed=0, randphz=0)
        return

    fresh = type(cur)(cur.command)

    def val(name, conv):
        t = _token(inc.get(name, cur.get(name)))
        try:
            return conv(t) if t else fresh.get(name)         # a blank field takes its default
        except NumberError:
            return None                                      # reported by the INCOH record check
    hseed = val("hseed", lambda t: parse_int(t)[0])
    vseed = val("vseed", lambda t: parse_int(t)[0])
    phase = val("randphz", parse_float)
    if phase is not None and phase == 0:
        inc["randphz"] = STOCHASTIC_PHASE
        phase = STOCHASTIC_PHASE
        notes.append(f"Stochastically Simulated Incoherency Input: Random Phase Angle {fmt_num(STOCHASTIC_PHASE)} "
                     f"degrees (requirements 5.4 default)")
    if hseed is not None and vseed is not None and hseed == 0 and vseed == 0:
        errors.append("Stochastically Simulated Incoherency Input needs a non-zero Horizontal or Vertical SEED "
                      "(D-INC-02: with both SEEDs 0 the input is deterministic) -- enter the SEED integers or "
                      "choose Deterministic (Median) Incoherency Input")
    if phase is not None and phase < 0:
        errors.append(f"Stochastically Simulated Incoherency Input: Random Phase Angle {fmt_num(phase)} must be > 0 "
                      f"(phases uniform in [-RandPhz, +RandPhz], D-INC-02)")


def _request_lines(model, name: str, rows: List[Dict[str, Any]]) -> List[str]:
    """Clear + one command per row when a request table changed (request append class, D-PAR-14)."""
    if name == "NOUT":
        new = []
        for r in rows:
            d = int(parse_int(_token(r.get("dir", 1)))[0])
            if not 1 <= d <= 6:
                raise ValueError(f"direction {d} (1..6)")
            codes = [1 if _truthy(r.get(f"c{k}", 0)) else 0 for k in range(1, 7)]
            nodes = parse_ids(r.get("nodes", ""))
            if not nodes:
                raise ValueError("a node list is empty")
            new.append((d, codes, nodes))
        old = [(q.dir, [int(v) for v in q.codes], list(q.nodes)) for q in model.nout]
        if new == old:
            return []
        return ["NOUT,0"] + [join_command("NOUT", [d] + codes + _ids_text(nodes).split())
                             for d, codes, nodes in new]
    if name == "EOUT":
        new = []
        for r in rows:
            g = int(parse_int(_token(r.get("group", 0)))[0])
            elems = parse_ids(r.get("elements", ""))
            codes = [int(parse_int(_token(c))[0]) if _token(c) else 0 for c in (list(r.get("codes") or []) + [0] * 12)[:12]]
            if any(c not in (0, 1, 2) for c in codes):
                raise ValueError("output codes are 0, 1 or 2")
            if not elems:
                raise ValueError("an element list is empty")
            new.append((codes, g, elems))
        old = [([int(v) for v in q.codes], q.group, list(q.elements)) for q in model.eout]
        if new == old:
            return []
        return ["EOUT,0"] + [join_command("EOUT", codes + [g] + _ids_text(elems).split())
                             for codes, g, elems in new]
    if name == "RDND":
        new = []
        for r in rows:
            node = int(parse_int(_token(r.get("node", 0)))[0])
            if node <= 0:
                raise ValueError("node numbers must be positive")
            new.append((node, [1 if _truthy(r.get(k, 0)) else 0 for k in ("x", "y", "z", "xx", "yy", "zz")]))
        old = [(q.node, [1 if v >= 1 else 0 for v in q.flags]) for q in model.rdnd]
        if new == old:
            return []
        return ["RDND,0"] + [join_command("RDND", [n] + fl) for n, fl in new]
    if name == "CORR":
        out = []
        stored = {k: r for k, r in get_entries(model, "CORR")}
        cls = record_class("CORR")
        for r in rows:
            rec = cls("CORR", [_token(r.get("no")), _token(r.get("time")), _token(r.get("val"))])
            if rec.token_errors():
                raise ValueError(rec.token_errors()[0])
            for p in rec.problems():
                if p[0] == "Error":
                    raise ValueError(problem_message(p))
            old = stored.get(rec.no)
            if old is None or not _same_typed(rec, typed(old), rec.FIELDS):
                out.append(to_command(rec))
        return out
    raise ValueError("not a request table of this dialog")


def _truthy(v) -> bool:
    if isinstance(v, str):
        return v.strip() not in ("", "0", "false", "False", "off")
    return bool(v)


# ======================================================================================
# Options > Write / Check (session settings without a command)
# ======================================================================================
def apply_session_dialog(interp, name: str, payload: Dict[str, Any]) -> List[str]:
    """OK of Options > Write / Check: set the session settings; returns confirmation notes."""
    rec = dict((payload.get("records") or {}).get(name.upper(), {}))
    if name.upper() == "WRITE":
        wo = interp.write_options
        for k in ("mdl", "ext", "afwr", "sim"):
            if k in rec:
                wo[k] = _truthy(rec[k])
        if "sim_location" in rec:
            loc = str(rec["sim_location"]).lower()
            if loc not in ("sim", "pre"):
                raise DialogError([f"Simulation Command Location '{loc}' (sim or pre)"])
            wo["sim_location"] = loc
        return [f"Extended Write Options: {', '.join(k for k in ('mdl', 'ext', 'afwr', 'sim') if wo.get(k)) or 'none'}"
                + (f"; simulation commands to {'*-Sim.Pre' if wo.get('sim_location', 'sim') == 'sim' else '*.Pre'}"
                   if wo.get("sim") else "")]
    if name.upper() == "CHECK":
        from ..prep.commands.checks import check_options
        co = check_options(interp)
        new = CheckOptions(co.show_warnings, co.show_errors, co.suppress_window, co.break_at)
        for k in ("show_warnings", "show_errors", "suppress_window"):
            if k in rec:
                setattr(new, k, _truthy(rec[k]))
        if "break_at" in rec:
            try:
                b = parse_int(_token(rec["break_at"]))[0]
            except NumberError:
                raise DialogError([f"Break Check at: '{rec['break_at']}' is not an integer"]) from None
            if b < 0:
                raise DialogError(["Break Check at: must be >= 0"])
            new.break_at = b
        interp.session["check_options"] = new
        return [f"Check Options: warnings {'shown' if new.show_warnings else 'hidden'}, errors "
                f"{'shown' if new.show_errors else 'hidden'}, error window "
                f"{'suppressed' if new.suppress_window else 'shown'}, break at {new.break_at}"]
    raise KeyError(name)
