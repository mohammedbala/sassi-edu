"""Module input decks: the normative contract between AFWRITE (sassi.prep) and the modules.

Each SSI module reads exactly one deck ``<model>.<ext>`` written by AFWRITE.  The
deck holds every value the module needs, already *resolved* from the model database
(e.g. the SITE deck carries the TOPL layers with their L-table properties, the
frequency numbers of the selected FREQ set and the resolved frequency step ``df``).
Modules never read the ``.pre`` file or the UI model.

Usage::

    from sassi.io import decks
    d = decks.new("SITE")                 # Deck with all params at their defaults
    d.params["nl"] = 20
    d.table("layers").append([1, 10.0, 0.12, 1500.0, 800.0, 0.05, 0.05])
    decks.write(path, d)
    d2 = decks.read(path)                 # types coerced, missing params defaulted

The schemas below list, for every module, the parameter names (named after the
command arguments in docs/spec/00_requirements.md section 3.4.F and the extension
commands of section 3.4.R), their types and defaults, and the table layouts.
Units are the model's consistent unit system; ``weight`` is specific weight
(force/length^3) and ``gravity`` converts weight to mass.  Damping values are ratios
(fractions), never percent, except DYNP curve points which keep the manual's percent
convention (strain % and damping %).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from .deckfmt import RawDeck, Table, parse_value, read_raw, write_raw

DECK_VERSION = 1


@dataclass(frozen=True)
class P:
    """Scalar or list parameter."""
    name: str
    type: type          # int, float, str, list
    default: Any
    doc: str = ""
    item: type = float  # element type when type is list


@dataclass(frozen=True)
class T:
    """Table specification: ordered (column, type) pairs."""
    name: str
    columns: Tuple[Tuple[str, type], ...]
    doc: str = ""


@dataclass(frozen=True)
class DeckSchema:
    module: str
    ext: str
    params: Tuple[P, ...]
    tables: Tuple[T, ...] = ()
    doc: str = ""

    def param(self, name: str) -> P:
        for p in self.params:
            if p.name == name:
                return p
        raise KeyError(name)

    def table_spec(self, name: str) -> T:
        for t in self.tables:
            if t.name == name:
                return t
        raise KeyError(name)


def _cols(spec: str, default: type = float) -> Tuple[Tuple[str, type], ...]:
    """'node:int x y z label:str' -> ((node,int),(x,float),...)."""
    out = []
    for tok in spec.split():
        if ":" in tok:
            n, t = tok.split(":")
            out.append((n, {"int": int, "float": float, "str": str}[t]))
        else:
            out.append((tok, default))
    return tuple(out)


# Common parameter groups ---------------------------------------------------------------
_COMMON = (
    P("title", str, "", "model title (TIT)"),
    P("model", str, "", "model name (MDL)"),
    P("opmode", int, 0, "0 complete solution, 1 data check only"),
)
_FFT = (
    P("delt", float, 0.005, "time step of the control motion, s (SITE <delt>)"),
    P("nft", int, 4096, "number of Fourier components NFFT, power of 2 (SITE <nft>)"),
    P("df", float, 0.0, "resolved frequency step, Hz (fstep if > 0 else 1/(delt*nft))"),
)
_HIST = (  # control-motion history data shared by MOTION, STRESS, RELDISP
    P("thfile", str, "", "control-motion acceleration history file (THFILE); accelerations in g"),
    P("thtit", str, "", "history title (THTIT)"),
    P("mult", float, 1.0, "multiplication factor (exactly one of mult/max non-zero)"),
    P("max", float, 0.0, "target peak value in g (history rescaled)"),
    P("rec1", int, 1, "first record (1-based)"),
    P("rec2", int, 0, "last record (0 = last in file)"),
    P("fopt", int, 0, "0: first line dt then one value per line; 1: (t, a) pairs per line"),
    P("dur", float, 0.0, "duration of output histories, s (0 = NFFT*delt); output is 1.2*dur"),
    P("type", int, 0, "analysis type (ANALYS <type>): 0 seismic, 1 foundation vibration"),
    P("gravity", float, 32.2, "SSI acceleration of gravity (HOUSE <gravity>)"),
    P("file8", str, "FILE8", "name of the solution file to read"),
    P("maxnode", int, 0, "largest node id of the model (file-name width 5 or 6)"),
)
_FREQS = T("freqs", _cols("number:int"), "SSI frequency numbers n (f = n*df), sorted ascending")
_DAMP = T("damp", _cols("value"), "response-spectrum damping ratios (DAMP), fractions")
_LAYER_COLS = _cols("no:int thick weight vp vs dp ds")

# ---------------------------------------------------------------------------------------
SCHEMAS: Dict[str, DeckSchema] = {}


def _reg(s: DeckSchema) -> DeckSchema:
    SCHEMAS[s.module] = s
    return s


EQUAKE = _reg(DeckSchema("EQUAKE", ".equ", _COMMON + (
    P("accopt", int, 0, "0 random phases, 1 seed record (keep Fourier phases of ACCIN), 2 external history (RS/PSD/FFT only)"),
    P("nrfreq", int, 0, "number of records (frequency points) in each RSIN file"),
    P("rand", int, 11975, "initial random SEED"),
    P("damp", float, 0.05, "damping ratio of the target spectrum"),
    P("dur", float, 20.0, "total motion duration, s"),
    P("corr", int, 0, "1 = correlated X-Y components (CORR pairs)"),
    P("seeds", int, 1, "number of random-seed trials (best fit kept)"),
    P("tpsd", int, 0, "1 = check against the target PSD files"),
    P("delt", float, 0.005, "time step, s (SITE <delt>)"),
    P("gravity", float, 32.2, "SSI gravity (sets British/SI output units)"),
    P("eqtit", str, "", "spectra title (EQTIT)"),
), (
    T("spectra", _cols("no:int rsin:str rsout:str accin:str accout:str tpsd_file:str"),
      "per spectrum 1..3: target RS file, output RS file, seed/external input history, output history, target PSD"),
    T("corr", _cols("no:int time val"), "time-varying X-Y correlation pairs (P1)"),
)))

SOIL = _reg(DeckSchema("SOIL", ".soi", _COMMON + (
    P("nrval", int, 0, "number of acceleration values read from the history file"),
    P("grav", float, 32.2, "gravity for the free-field SOIL analysis"),
    P("header", int, 0, "number of header lines in the history file"),
    P("outcrop", int, 1, "1 input is an outcrop motion at the control layer top, 0 within motion"),
    P("save", int, 1, "1 save strain-compatible properties to FILE88"),
    P("iter", int, 8, "number of equivalent-linear iterations"),
    P("ratio", float, 0.65, "effective / maximum strain ratio"),
    P("gravmult", float, 1.0, "multiplier for g in the RS output"),
    P("cof", float, 0.0, "cut-off frequency (always Nyquist); see soilcutoff"),
    P("soilcutoff", float, 0.0, "EDUOPT,SOILCUTOFF: >0 removes Fourier components above this frequency (SHAKE91 runs)"),
    P("delt", float, 0.005, "time step of the input history, s"),
    P("nft", int, 4096, "NFFT"),
    P("cl", int, 1, "control point layer (input applied at top of this sublayer)"),
    P("thfile", str, "", "input acceleration history (g), free format after the header lines"),
    P("thtit", str, "", "history title"),
    P("mult", float, 0.0, "scale factor (exactly one of mult/max non-zero)"),
    P("max", float, 0.1, "target PGA in g"),
    P("indir", int, 0, "0 horizontal (Vs, iterations), 1 vertical (Vp, no iterations)"),
    P("vppolicy", str, "NU", "EDUOPT,VPPOLICY: 'NU' Vp of iterated sublayers at constant Poisson's ratio and beta_p = beta_s; 'VP' keep input Vp and beta_p (D-SOL-06)"),
    P("cmodform", int, 0, "complex-modulus form (CMODFORM)"),
), (
    T("profile", _cols("layer:int thick weight vp vs dp ds dynprop:str"),
      "SPRO sublayers top first, properties resolved from L; the LAST row is the half-space"),
    T("dynp", _cols("label:str kind:str strain value"),
      "curve points: kind 'G' (strain %, G/Gmax) or 'D' (strain %, damping %)"),
    T("sacc", _cols("layer:int opt:int outcrop:int"), "acceleration output: opt 0 none, 1 max, 2 max + history"),
    T("srs", _cols("layer:int save:int outcrop:int"), "response-spectrum output at layer tops"),
    T("sstr", _cols("layer:int opt1:int opt2:int opt3:int opt4:int"),
      "stress/strain output: opt1 compute stress, opt2 save SSxxx.TH, opt3 compute strain, opt4 save SNxxx.TH (spec 07 9.2.36)"),
    T("ssaf", _cols("layer:int save:int outcrop1:int outcrop2:int layer2:int freqstep title:str"),
      "spectral amplification factor RS(layer)/RS(layer2)"),
    _DAMP,
)))

SITE = _reg(DeckSchema("SITE", ".sit", _COMMON + (
    P("mode1", int, 1, "1 compute Rayleigh/Love eigen-solutions -> FILE2"),
    P("fstep", float, 0.0, "frequency step as entered (0 -> 1/(delt*nft))"),
    P("nl", int, 20, "number of generated half-space sublayers (0 rigid base, else 4..20)"),
    P("hs", int, 0, "L number of the half-space layer"),
    P("mode2", int, 1, "1 compute free-field motion -> FILE1"),
    P("wopt", int, 0, "0 R/SV/P (in-plane x'z'), 1 SH/L (anti-plane y')"),
    P("freq1", int, 1, "frequency number 1 of the wave-ratio curve"),
    P("freq2", int, 2048, "frequency number 2 of the wave-ratio curve"),
    P("cl", int, 1, "control point layer (top of this TOPL layer; 1 = surface)"),
    P("cm", int, 0, "control motion direction 0 x', 1 y', 2 z'"),
    P("freq", int, 1, "frequency set number"),
    P("soilmode", int, 0, "SITEX: 0 linear soil (L table), 1 non-linear soil (FILE88 properties)"),
    P("gravity", float, 32.2, "SSI gravity (weight -> mass)"),
    P("cmodform", int, 0, "complex-modulus form"),
    P("hslaw", str, "uniform", "EDUOPT,HSLAW: half-space sublayer law 'uniform' (default, gated by VP-02b), 'geometric' or 'linear' (D-SIT-02)"),
) + _FFT, (
    T("layers", _LAYER_COLS, "TOPL layers, top first (L-table properties); interface i is the top of layer i"),
    T("halfspace", _LAYER_COLS, "one row: the half-space layer (thick ignored)"),
    T("waves", _cols("type:int opt:int ratio1 ratio2 angle"),
      "WAVE records: type 1 R, 2 SV, 3 P, 4 SH, 5 L; opt 0 off/1 field (R: 1 shortest wavelength, 2 least decay)"),
    _FREQS,
)))

POINT = _reg(DeckSchema("POINT", ".poi", _COMMON + (
    P("layer", int, 0, "last layer number of the near-field zone (0 surface): loads at interfaces 1..layer+1"),
    P("rad", float, 1.0, "radius of the central zone R0"),
    P("dim", int, 2, "HOUSE <dim>: 1 2D (POINT2), 2 3D (POINT3)"),
    P("freq", int, 1, "frequency set number (informative)"),
    P("df", float, 0.0, "resolved frequency step, Hz"),
), (_FREQS,)))

HOUSE = _reg(DeckSchema("HOUSE", ".hou", _COMMON + (
    P("gravity", float, 32.2, "SSI acceleration of gravity"),
    P("gelev", float, 0.0, "ground surface elevation"),
    P("dim", int, 2, "1 2D, 2 3D"),
    P("imp", int, 0, "0 FV, 1 FFV, 2 FI (informative; the interaction set defines the method)"),
    P("coh", int, 0, "0 coherent, 1 incoherent"),
    P("wpass", int, 0, "1 wave passage"),
    P("me", int, 0, "1 multiple excitation"),
    P("cmplxspec", int, 0, "1 complex spectral amplification ratios"),
    P("optimize", int, 0, "HOUSEX: node renumbering optimizer (.hounew/.map)"),
    P("supmode", int, 0, "HOUSEX: incoherent superposition 0 linear (AS), 1 quadratic"),
    P("nsim", int, 1, "HOUSEX: number of stochastic simulations"),
    P("nlssi", int, 0, "HOUSEX: non-linear soil SSI (.pin)"),
    P("ansys", int, 0, "HOUSEX: ANSYS model input (Option AA, P2)"),
    P("incomp", int, 1, "MOPT <incomp>: 0 include incompatible modes, 1 suppress"),
    P("gmunits", int, 0, "MOPT <matrix>: GENERAL mass given as 0 mass, 1 weight"),
    P("cmodform", int, 0, "complex-modulus form"),
    # incoherency (INCOH) and wave passage (WPASS) -- P1
    P("gammax", float, 0.1), P("gammay", float, 0.1), P("gammaz", float, 0.2),
    P("alpha", float, 0.5, "directionality factor (or Vs for Luco-Wong)"),
    P("ngp", int, 1), P("ipr", int, 0), P("nmodes", int, 0), P("met", int, 0),
    P("hseed", int, 0), P("vseed", int, 0), P("randphz", float, 0.0),
    P("appv", float, 1.0e9, "apparent velocity for Line D"),
    P("wang", float, 0.0, "angle of Line D with the X axis, deg"),
    P("cohf", int, 1, "unlagged coherency model 1..7"),
    P("xc", float, 0.0, "control point (ANALYS <xc>) for wave passage"),
    P("yc", float, 0.0), P("zc", float, 0.0),
) + _FFT, (
    T("nodes", _cols("id:int x y z fx:int fy:int fz:int frx:int fry:int frz:int"),
      "global coordinates; fixity flags 1 = fixed (D command)"),
    T("interaction", _cols("id:int"), "interaction nodes (INT code 0) in interaction order"),
    T("groups", _cols("id:int type:int title:str"), "type 1 SOLID, 2 BEAMS, 3 SHELL, 4 PLANE, 5 TSHELL, 7 SPRING, 9 GENERAL"),
    T("elements", _cols("group:int id:int etype:int mat:int prop:int eint:int thick n1:int n2:int n3:int n4:int "
                        "n5:int n6:int n7:int n8:int ki:str kj:str"),
      "unused node slots = 0; ki/kj are 6-digit release strings (beams) e.g. '000011'; etype 0 implicit, 1 structure, 2 excavated/buried"),
    T("materials", _cols("id:int type:int val1 val2 weight pdamp sdamp"),
      "M table: type 1 (E, nu), 2 (constrained M, G), 3 (Vp, Vs)"),
    T("layers", _LAYER_COLS, "L table (all defined soil layers)"),
    T("sitelayers", _LAYER_COLS, "TOPL layers top first (from the SITE options) incl. a final half-space row with thick 0"),
    T("beamprops", _cols("id:int axial shear2 shear3 tors flex2 flex3"), "R table"),
    T("springprops", _cols("id:int scx scy scz scxx scyy sczz damp"), "SC table (global axes)"),
    T("matrices", _cols("prop:int kind:str row:int t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12"),
      "GENERAL element matrix rows, kind R/I/M, upper triangle (row r has 13-r terms, rest 0)"),
    T("masses", _cols("node:int mx my mz mxx myy mzz units:int"), "nodal masses (MT/MR); units 0 mass, 1 weight"),
    T("me", _cols("no:int nfirst:int nlast:int"), "multiple-excitation zones (P1)"),
    T("amp", _cols("no:int idx:int re im"), "spectral amplification ratios per zone and frequency index (P1)"),
    T("symm", _cols("no:int type:int n1:int n2:int n3:int"), "symmetry planes (P1)"),
    _FREQS,
)))

FORCE = _reg(DeckSchema("FORCE", ".frc", _COMMON + (
    P("gravity", float, 32.2),
    P("mforce", int, 1, "MOPT <force>: 0 add repeated definitions, 1 overwrite"),
    P("freq", int, 1, "frequency set number"),
    P("fstep", float, 0.0),
) + _FFT, (
    T("loads", _cols("node:int dof:int factor arrival"), "F (dof 1-3) and MM (dof 4-6) factors and arrival times, s"),
    _FREQS,
)))

ANALYS = _reg(DeckSchema("ANALYS", ".anl", _COMMON + (
    P("type", int, 0, "0 seismic (FILE1), 1 foundation vibration (FILE9)"),
    P("mode", int, 0, "0 initiation, 1 new structure, 2 new seismic environment, 3 new dynamic loading"),
    P("save", int, 0, "1 save restart files COOXqqq/COOTKqqq"),
    P("prnt", int, 1, "1 print amplitudes only, 0 complex"),
    P("fopt", int, 0, "0 frequency set, 1 all frequencies of FILE1/FILE9"),
    P("ang", float, 0.0, "coordinate transformation angle x' -> x, deg"),
    P("xc", float, 0.0), P("yc", float, 0.0), P("zc", float, 0.0),
    P("impe", int, 0, "global impedance 0 none, 1 diagonal, 2 full 6x6"),
    P("simul", int, 0, "simultaneous cases: 0 single; seismic coherent 1 = X/Y/Z (FILE1X/Y/Z); incoherent Ns; vibration Nl"),
    P("ffm", int, 0, "ANALYSX: 0 free-field load (FFL), 1 free-field motion (FFM)"),
    P("delrst", int, 0, "ANALYSX: delete restart files after success"),
    P("coh", int, 0), P("wpass", int, 0), P("me", int, 0),
    P("freq", int, 1, "frequency set number"),
    P("gravity", float, 32.2),
) + _FFT, (_FREQS,)))

MOTION = _reg(DeckSchema("MOTION", ".mot", _COMMON + (
    P("out", int, 0, "0 full output, 1 transfer functions only"),
    P("step", int, 0, "print step of histories in the listing (0 table only)"),
    P("res", int, 0, "not used"),
    P("freq1", float, 0.1, "first RS frequency, Hz"),
    P("freq2", float, 100.0, "last RS frequency, Hz"),
    P("fstep", int, 301, "number of RS frequencies (log spaced)"),
    P("bl", int, 0, "0 no baseline correction, 1 Hudson-Housner correction"),
    P("smo", float, 0.0, "smoothing parameter S (options 0-5)"),
    P("cplx", int, 1, "1 save complex TFU/TFI (needed by RELDISP)"),
    P("cnvrt", int, 0, "1 compute RS of external histories listed in CONTTRS.txt"),
    P("pzadj", int, 0, "1 phase adjustment"),
    P("interp", int, 1, "interpolation option 0..6"),
    P("f1213", int, 0, "MOTIONX: 0 none, 1 FILE13, 2 FILE12"),
    P("resp", int, 2, "MOTIONX: vibration response 0 displacement, 1 velocity, 2 acceleration"),
    P("srss", int, 0, "MOTIONX: SRSS TF from SRSSTF.txt"),
    P("savetf", int, 0), P("saveacc", int, 0), P("savers", int, 0), P("saverot", int, 0),
    P("rsttf", int, 0), P("rstacc", int, 0), P("rstrs", int, 0),
    P("ang", float, 0.0, "ANALYS <ang> (zero-frequency anchor of rotated inputs)"),
    P("cm", int, 0, "SITE <cm> control direction (zero-frequency anchor)"),
) + _FFT + _HIST, (
    _DAMP,
    T("nout", _cols("node:int dir:int c1:int c2:int c3:int c4:int c5:int c6:int"),
      "nodal requests; dir 1 x .. 6 zz; flags: TF print, save TH, plot TH, plot RS (deprecated), save RS, print max"),
)))

STRESS = _reg(DeckSchema("STRESS", ".str", _COMMON + (
    P("iter", int, 0, "1 automatic strains in soil elements (nonlinear soil)"),
    P("save", int, 1, "1 save stress time histories (.THS)"),
    P("itran", int, 0, "1 output stress TFs (.TFU/.TFI)"),
    P("interopt", int, 1, "interpolation option 0..6"),
    P("pzadj", int, 0, "STRESSX phase adjustment"),
    P("smo", float, 0.0, "STRESSX smoothing parameter"),
    P("skip", int, 0, "STRESSX: skip time-history steps in output"),
    P("savemax", int, 0, "STRESSX: all-element maxima (ELEMENT_CENTER_ABS_MAX_STRESSES.TXT)"),
    P("saveth", int, 0, "STRESSX: all-element histories"),
    P("rstns", int, 0, "STRESSX: restart for nodal stress contours (frames)"),
    P("rstsp", int, 0, "STRESSX: restart for soil pressure contours (frames)"),
    P("secdataopt", int, 0, "SECDATAOPT: write ESTRESS_n.ess frames"),
    P("thshlstr", int, -1, "THSHLSTR: 0 the 8 basic TSHELL components, 1 also face stresses/strains; "
      "-1 = not set (decks not written by AFWRITE: STRESS falls back to THSHLSTR.opt)"),
    P("cm", int, 0), P("ang", float, 0.0),
) + _FFT + _HIST, (
    T("eout", _cols("group:int element:int c1:int c2:int c3:int c4:int c5:int c6:int c7:int c8:int c9:int c10:int c11:int c12:int"),
      "element requests; code per component 0 none, 1 max, 2 max + history (component order of ELEMENT_COMPONENTS)"),
)))

RELDISP = _reg(DeckSchema("RELDISP", ".rdi", _COMMON + (
    P("relfile", str, "", "reference-node complex .TFI file (RELFILE); empty = free-field unit reference"),
    P("reldisoutput", int, 1, "1 save relative-displacement complex TF (.TFD)"),
    P("reldispsall", int, 0, "1 relative displacements at all nodes"),
    P("numfiles", int, 1, "number of output files"),
    P("saverot", int, 0, "RELDX"), P("rstframes", int, 0, "RELDX"),
    P("cm", int, 0), P("ang", float, 0.0),
) + _FFT + _HIST, (
    T("rdnd", _cols("node:int x:int y:int z:int xx:int yy:int zz:int"), "nodes/DOFs requested (>= 1 = on)"),
)))


# ---------------------------------------------------------------------------------------
@dataclass
class Deck:
    """Typed deck: ``params`` (name -> value) and ``tables`` (name -> Table)."""

    module: str
    params: Dict[str, Any] = field(default_factory=dict)
    tables: Dict[str, Table] = field(default_factory=dict)

    @property
    def schema(self) -> DeckSchema:
        return SCHEMAS[self.module]

    def __getitem__(self, name: str) -> Any:
        return self.params[name]

    def __setitem__(self, name: str, value: Any) -> None:
        self.params[name] = value

    def get(self, name: str, default=None):
        return self.params.get(name, default)

    def table(self, name: str) -> Table:
        if name not in self.tables:
            spec = self.schema.table_spec(name)
            self.tables[name] = Table(columns=[c for c, _ in spec.columns])
        return self.tables[name]

    def rows(self, name: str) -> List[Dict[str, Any]]:
        return self.table(name).as_dicts()


def new(module: str) -> Deck:
    """A deck with every parameter at its default and empty tables."""
    s = SCHEMAS[module.upper()]
    d = Deck(module=s.module, params={p.name: (list(p.default) if isinstance(p.default, list) else p.default)
                                       for p in s.params})
    for t in s.tables:
        d.tables[t.name] = Table(columns=[c for c, _ in t.columns])
    return d


def deck_path(model_dir: Union[str, Path], model: str, module: str) -> Path:
    return Path(model_dir) / f"{model}{SCHEMAS[module.upper()].ext}"


def _coerce(v: Any, typ: type, item: type = float):
    if typ is list:
        if isinstance(v, str):
            v = parse_value(v)
            if isinstance(v, str):
                v = [x for x in v.replace(",", " ").split() if x]
        return [_coerce(x, item) for x in v]
    if typ is str:
        if isinstance(v, str):
            t = v.strip()
            if t.startswith('"'):
                return parse_value(t)
            return v
        return str(v)
    if typ is int:
        if isinstance(v, str):
            t = v.strip().replace("D", "E").replace("d", "e")
            f = float(t) if t else 0.0
        else:
            f = float(v)
        if f != int(round(f)):
            raise ValueError(f"expected an integer, got {v!r}")
        return int(round(f))
    if typ is float:
        if isinstance(v, str):
            t = v.strip().replace("D", "E").replace("d", "e")
            return float(t) if t else 0.0
        return float(v)
    raise TypeError(typ)


def write(path: Union[str, Path], deck: Deck, comments: Sequence[str] = ()) -> Path:
    """Validate ``deck`` against its schema and write it."""
    s = deck.schema
    params: Dict[str, Any] = {}
    known = {p.name for p in s.params}
    for p in s.params:
        params[p.name] = _coerce(deck.params.get(p.name, p.default), p.type, p.item)
    unknown = set(deck.params) - known
    if unknown:
        raise KeyError(f"{s.module} deck: unknown parameters {sorted(unknown)}")
    tables: Dict[str, Table] = {}
    for t in s.tables:
        tab = deck.tables.get(t.name)
        cols = [c for c, _ in t.columns]
        out = Table(columns=cols)
        if tab is not None:
            if list(tab.columns) != cols:
                raise ValueError(f"{s.module} deck table {t.name}: columns {tab.columns} != {cols}")
            for r in tab.rows:
                out.rows.append([_coerce(v, typ) for v, (_, typ) in zip(r, t.columns)])
        tables[t.name] = out
    extra = set(deck.tables) - {t.name for t in s.tables}
    if extra:
        raise KeyError(f"{s.module} deck: unknown tables {sorted(extra)}")
    return write_raw(path, s.module, params, tables, version=DECK_VERSION, comments=comments)


def read(path: Union[str, Path], module: Optional[str] = None) -> Deck:
    """Read and type a deck; missing parameters take their defaults (with a warning)."""
    raw: RawDeck = read_raw(path)
    if module is not None and raw.module != module.upper():
        raise ValueError(f"{Path(path).name}: expected a {module.upper()} deck, found {raw.module}")
    s = SCHEMAS[raw.module]
    d = new(raw.module)
    for p in s.params:
        if p.name in raw.params:
            d.params[p.name] = _coerce(raw.params[p.name], p.type, p.item)
        else:
            warnings.warn(f"{raw.module} deck: parameter '{p.name}' missing, default {p.default!r} used")
    for name in raw.params:
        if name not in d.params:
            warnings.warn(f"{raw.module} deck: unknown parameter '{name}' ignored")
    for t in s.tables:
        rt = raw.tables.get(t.name)
        if rt is None:
            continue
        cols = [c for c, _ in t.columns]
        if rt.columns != cols:
            raise ValueError(f"{raw.module} deck table {t.name}: columns {rt.columns} != {cols}")
        tab = d.tables[t.name]
        for r in rt.rows:
            tab.rows.append([_coerce(v, typ) for v, (_, typ) in zip(r, t.columns)])
    return d
