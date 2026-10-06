"""Typed analysis-option records of the module option commands (requirements sections 3.3, 3.4.F,
3.4.R and 5.4; spec 07 sections 3-4; specs 05a-05d).

Every module option command of ACS SASSI stores its arguments in a *record* of the model
(:class:`sassi.model.OptionStore`).  This module declares, for each command, the argument names
(the names used in the manual's syntax line, which are also the parameter names of the module
decks in :mod:`sassi.io.decks`), their types, their **new-model defaults** (requirements section
5.4 "Default" column, decision D-UI-03) and their dialog labels.  The records keep the token text
the user typed, so WRITE -> INP is exact (UT-03), and give typed, named access::

    rec = get_record(model, "SITE")       # typed SiteRecord (defaults when SITE was never given)
    rec.nft, rec.get("delt")              # 4096, 0.005
    rec.problems()                        # field-local CHECK rules: [('Error', 47, '...'), ...]

The GUI (Options > Analysis dialog) works with plain dataclasses whose field names are the
command argument names, whose defaults are the new-model defaults and whose field metadata carry
the dialog label::

    opts = from_args("SITE,0,1,0,20,3")   # SITEOptions(opmode=0, mode1=1, fstep=0.0, nl=20, hs=3, ...)
    opts.nl = 10
    to_command(opts)                      # 'SITE,0,1,0,10,3,1,0,1,2048,1,0,0.005,4096,1'

Rule L17: a dialog commit submits the text of :func:`to_command` to the interpreter, so a GUI
session can be replayed as a ``.pre`` file.

Shared analysis variables (spec 07 section 3) have **one** storage location (their owning
command argument).  :data:`SHARED_VARIABLES` lists them with the tabs that show them, and
:func:`shared_value` / :func:`shared_command` read them and build the command that changes them.
"""
from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass
from typing import Any, ClassVar, Dict, List, Optional, Sequence, Tuple, Type, Union

from ..conventions import is_power_of_two
from ..model.options import Field, Record, make_record, register_record_type
from ..model.values import NumberError, fmt_num, parse_float, parse_int

class Problem(tuple):
    """A violated CHECK rule: the tuple ``(kind, number, detail)`` plus the values ``fmt`` of the
    catalogue placeholders of its title (manual Chapter 10, spec 11 section 5).

    ``kind`` is ``'Error'`` or ``'Warning'``, ``number`` the catalogue number (or an ``EDU-nn``
    label), ``detail`` the offending values.  Titles such as Error 53 "Illegal Value for Frequency
    <i>" name an index that cannot be recovered reliably from the detail text, so the rule that
    detects the problem states it explicitly: ``Problem('Error', 53, '<freq1> = 0', i=1)``.
    A ``Problem`` compares and unpacks like a plain 3-tuple.
    """

    def __new__(cls, kind: str, number: Union[int, str], detail: str = "", **fmt: Any) -> "Problem":
        self = tuple.__new__(cls, (kind, number, detail))
        self.fmt = dict(fmt)
        return self

    @property
    def kind(self) -> str:
        return self[0]

    @property
    def number(self) -> Union[int, str]:
        return self[1]

    @property
    def detail(self) -> str:
        return self[2]

    def __getnewargs_ex__(self):                 # copy / pickle rebuild with the same arguments
        return (self[0], self[1], self[2]), dict(self.fmt)

    def __repr__(self) -> str:
        extra = "".join(f", {k}={v!r}" for k, v in self.fmt.items())
        return f"Problem({self[0]!r}, {self[1]!r}, {self[2]!r}{extra})"


def problem_fmt(p: Sequence[Any]) -> Dict[str, Any]:
    """Placeholder values of a problem (empty for a plain tuple)."""
    return dict(getattr(p, "fmt", {}) or {})


# ======================================================================================
# Field and record base classes
# ======================================================================================
@dataclass(frozen=True)
class OptField(Field):
    """A dialog field / command argument: :class:`sassi.model.Field` plus GUI metadata.

    ``doc`` is the dialog label (spec 05a-05d), ``choices`` the admissible codes of radio buttons
    and check boxes (informative for the GUI; CHECK enforces the rules), ``shared`` the name of
    the shared variable this argument owns (:data:`SHARED_VARIABLES`).
    """
    choices: Tuple[Any, ...] = ()
    shared: str = ""


def F(name: str, type_: type, default: Any, label: str = "", choices: Tuple[Any, ...] = (),
      shared: str = "") -> OptField:
    """Shorthand constructor of :class:`OptField`."""
    return OptField(name, type_, default, label, choices, shared)


class OptionRecord(Record):
    """Base class of the typed option records.

    Class attributes:

    * ``KIND`` -- ``'record'`` (record setter: the command replaces the record) or ``'indexed'``
      (indexed setter: the command creates/overwrites the entry at its key);
    * ``KEY`` -- field names forming the key of an indexed setter;
    * ``TAB`` -- the Options > Analysis tab that shows the record;
    * ``OPTIONAL_FROM`` -- first trailing *optional* argument (1-based), omitted by
      :func:`to_command` while it and the following ones are at their defaults;
    * ``XCOMMAND`` -- extension X-command (requirements section 3.4.R): written by WRITE only when
      non-default, and :func:`to_command` omits trailing defaults.
    """
    KIND: ClassVar[str] = "record"
    KEY: ClassVar[Tuple[str, ...]] = ()
    TAB: ClassVar[str] = ""
    TITLE: ClassVar[str] = ""
    OPTIONAL_FROM: ClassVar[Optional[int]] = None
    XCOMMAND: ClassVar[bool] = False

    def problems(self, shared: bool = True) -> List[Problem]:
        """Field-local CHECK rules of this record (manual Chapter 10 numbers).

        ``shared=False`` leaves out the rules on shared variables (CHECK reports those under
        every module that uses them, D-CHK-01).  Rules that need other records or the model are
        in :mod:`sassi.prep.check`.
        """
        return []

    # ---------------------------------------------------------------- token validation
    def token_errors(self) -> List[str]:
        """Messages for given tokens that cannot be read with their field type."""
        out = []
        for i, f in enumerate(self.FIELDS, start=1):
            t = self.arg(i)
            if t is None or f.type is str:
                continue
            try:
                parse_float(t)
            except NumberError:
                out.append(f"argument {i} <{f.name}>: '{t}' is not a number")
        return out

    def rounded_ints(self) -> List[str]:
        """Messages for integer fields given as reals with a fractional part (D-PAR-07)."""
        out = []
        for i, f in enumerate(self.FIELDS, start=1):
            t = self.arg(i)
            if t is None or f.type is not int:
                continue
            try:
                v, exact = parse_int(t)
            except NumberError:
                continue
            if not exact:
                out.append(f"argument {i} <{f.name}>: {t} rounded to {v}")
        return out


def _flag(v) -> bool:
    return bool(v) and v != 0


# ======================================================================================
# EQUAKE tab (spec 07 sections 9.2.1-9.2.13, 05a section 5)
# ======================================================================================
@register_record_type
class EquakeRecord(OptionRecord):
    """``EQUAKE,<accopt>,<nrfreq>,<rand>,<damp>,<dur>,<corr>,<seeds>,[tpsd]`` (D-EQK-01, D-UI-03).

    ``<nrfreq>`` blank means "the number of records of RSIN 1" (requirements section 5.4 default);
    ``<seeds>`` 0 is treated as 1 (D-EQK-05); the duration default is 20 s (D-UI-03).
    """
    COMMAND = "EQUAKE"
    TAB = "EQUAKE"
    OPTIONAL_FROM = 8
    FIELDS = (F("accopt", int, 0, "Accel. Record (1) / External Accel (2)", (0, 1, 2)),
              F("nrfreq", int, 0, "Number of Frequencies"),
              F("rand", int, 11975, "Initial Random SEED"),
              F("damp", float, 0.05, "Damping Value"),
              F("dur", float, 20.0, "Total Duration"),
              F("corr", int, 0, "Correlated", (0, 1)),
              F("seeds", int, 1, "Number of SEEDs"),
              F("tpsd", int, 0, "Use Target PSD", (0, 1)))

    def problems(self, shared: bool = True) -> List[Problem]:
        out: List[Problem] = []
        if self.rand <= 0:
            out.append(("Error", 90, f"<rand> = {self.rand}"))
        if self.given(2) and self.nrfreq <= 0 and self.accopt != 2:
            out.append(("Error", 91, f"<nrfreq> = {self.nrfreq}"))
        if self.dur <= 0:
            out.append(("Error", 92, f"<dur> = {fmt_num(self.dur)}"))
        return out


class _SpectrumFile(OptionRecord):
    """``<CMD>,<no>,<file>``: one file of spectrum ``<no>`` (1-3, one per direction)."""
    KIND = "indexed"
    KEY = ("no",)
    TAB = "EQUAKE"


@register_record_type
class RsinRecord(_SpectrumFile):
    """``RSIN,<no>,<file>``: target response-spectrum file (2 columns, ``<nrfreq>`` records)."""
    COMMAND = "RSIN"
    FIELDS = (F("no", int, 1, "Spectrum Number", (1, 2, 3)), F("file", str, "", "Spectrum Input File"))


@register_record_type
class RsoutRecord(_SpectrumFile):
    """``RSOUT,<no>,<file>``: response spectrum of the generated motion (``.rso``)."""
    COMMAND = "RSOUT"
    FIELDS = (F("no", int, 1, "Spectrum Number", (1, 2, 3)), F("file", str, "", "Spectrum Output File"))


@register_record_type
class AccinRecord(_SpectrumFile):
    """``ACCIN,<no>,<file>``: seed record (accopt 1) or external history (accopt 2)."""
    COMMAND = "ACCIN"
    FIELDS = (F("no", int, 1, "Spectrum Number", (1, 2, 3)), F("file", str, "", "Acceleration Input File"))


@register_record_type
class AccoutRecord(_SpectrumFile):
    """``ACCOUT,<no>,<file>``: generated acceleration history (``.acc``; ``.vel``/``.dis`` alike)."""
    COMMAND = "ACCOUT"
    FIELDS = (F("no", int, 1, "Spectrum Number", (1, 2, 3)), F("file", str, "", "Acceleration Output File"))


@register_record_type
class TpsdRecord(_SpectrumFile):
    """``TPSD,<num>,<file>``: target PSD of spectrum ``<num>`` (EQUAKE ``[tpsd]`` = 1)."""
    COMMAND = "TPSD"
    FIELDS = (F("no", int, 1, "Spectrum Number", (1, 2, 3)), F("file", str, "", "PSD File"))


@register_record_type
class CorrRecord(OptionRecord):
    """``CORR,<no>,<time>,<val>``: X-Y correlation pair (P1; EQUAKE ``<corr>`` = 1)."""
    COMMAND = "CORR"
    KIND = "indexed"
    KEY = ("no",)
    TAB = "EQUAKE"
    FIELDS = (F("no", int, 1, "Row"), F("time", float, 0.0, "Time"), F("val", float, 0.0, "Corr."))

    def problems(self, shared: bool = True) -> List[Problem]:
        return [("Error", 94, f"CORR {self.no}: {fmt_num(self.val)}")] if abs(self.val) > 1.0 else []


# ======================================================================================
# SOIL tab (spec 07 sections 9.2.10, 9.2.28-9.2.36; 05a section 6)
# ======================================================================================
@register_record_type
class SoilRecord(OptionRecord):
    """``SOIL,<nrval>,<grav>,<header>,<outcrop>,<save>,<iter>,<ratio>,<gravmult>,<cof>``.

    ``<cof>`` is always replaced by the Nyquist frequency (manual); the legacy PREP form with a
    final ``(format)`` token is stored but not mapped (D-SOL-10, Error 102).
    """
    COMMAND = "SOIL"
    TAB = "SOIL"
    FIELDS = (F("nrval", int, 0, "Number of Values"),
              F("grav", float, 32.2, "Gravity Accel. (used for free-field analysis)", shared="soil_gravity"),
              F("header", int, 0, "Number of Header Lines"),
              F("outcrop", int, 1, "Assign as Outcrop Motion", (0, 1)),
              F("save", int, 1, "Save Strain-Compatible Soil Properties", (0, 1)),
              F("iter", int, 8, "Number of Iterations"),
              F("ratio", float, 0.65, "Equiv. Uniform / Max Strain"),
              F("gravmult", float, 1.0, "Multiplier for Acceleration of Gravity"),
              F("cof", float, 0.0, "Cut-Off Frequency (always Nyquist)"))

    def problems(self, shared: bool = True) -> List[Problem]:
        out: List[Problem] = []
        if self.legacy:
            out.append(("Error", 102, f"legacy format {self.legacy} (arguments not mapped, D-SOL-10)"))
            return out
        if self.grav <= 0:
            out.append(("Error", 1, f"SOIL <grav> = {fmt_num(self.grav)}"))
        if self.nrval <= 0:
            out.append(("Error", 100, f"<nrval> = {self.nrval}"))
        if self.cof < 0:
            out.append(("Error", 101, f"<cof> = {fmt_num(self.cof)}"))
        if self.header < 0:
            out.append(("Error", 103, f"<header> = {self.header}"))
        if self.iter < 0:
            out.append(("Error", 105, f"<iter> = {self.iter}"))
        if not 0.0 < self.ratio < 1.0:
            out.append(("Error", 106, f"<ratio> = {fmt_num(self.ratio)}"))
        if self.gravmult <= 0:
            out.append(("Error", 108, f"<gravmult> = {fmt_num(self.gravmult)}"))
        return out


class _SoilLayerRequest(OptionRecord):
    KIND = "indexed"
    KEY = ("layer",)
    TAB = "SOIL"


@register_record_type
class SproRecord(_SoilLayerRequest):
    """``SPRO,<layer>,<prop>,<dynprop>``: SOIL sublayer -> L property + DYNP label (last = half-space)."""
    COMMAND = "SPRO"
    FIELDS = (F("layer", int, 1, "Layer Number"), F("prop", int, 0, "Property Number"),
              F("dynprop", str, "", "Dynamic Soil Property"))


@register_record_type
class DynpRecord(OptionRecord):
    """``DYNP,<no>,<sg>,<g>,<sd>,<d>,<label>``: curve point (strain %, G/Gmax; strain %, damping %)."""
    COMMAND = "DYNP"
    KIND = "indexed"
    KEY = ("label", "no")
    TAB = "SOIL"
    FIELDS = (F("no", int, 1, "Point"), F("sg", float, 0.0, "Strain (%)"), F("g", float, 1.0, "Mod. Red."),
              F("sd", float, 0.0, "Strain (%)"), F("d", float, 0.0, "Damp (%)"), F("label", str, "", "Title"))

    def problems(self, shared: bool = True) -> List[Problem]:
        if self.given(3) and not 0.0 <= self.g <= 1.0:
            return [Problem("Error", 98, f"{self.label} point {self.no}: G/Gmax = {fmt_num(self.g)}", p=self.label)]
        return []


@register_record_type
class SaccRecord(_SoilLayerRequest):
    """``SACC,<layer>,<opt>,<outcrop>``: acceleration at the top of a sublayer (0 none, 1 max, 2 + history)."""
    COMMAND = "SACC"
    FIELDS = (F("layer", int, 1, "Layer Number"), F("opt", int, 1, "Accelerations", (0, 1, 2)),
              F("outcrop", int, 0, "Outcropping", (0, 1)))


@register_record_type
class SrsRecord(_SoilLayerRequest):
    """``SRS,<layer>,<save>,<outcrop>``: response spectra at the top of a sublayer (DAMP ratios)."""
    COMMAND = "SRS"
    FIELDS = (F("layer", int, 1, "Layer Number"), F("save", int, 1, "Save Response Spectrum", (0, 1)),
              F("outcrop", int, 0, "Outcropping", (0, 1)))


@register_record_type
class SstrRecord(_SoilLayerRequest):
    """``SSTR,<layer>,<opt1..opt4>``: compute stress, save stress, compute strain, save strain."""
    COMMAND = "SSTR"
    FIELDS = (F("layer", int, 1, "Layer Number"), F("opt1", int, 1, "Compute Stresses", (0, 1)),
              F("opt2", int, 0, "Save Stress Time History", (0, 1)), F("opt3", int, 1, "Compute Strains", (0, 1)),
              F("opt4", int, 0, "Save Strain Time History", (0, 1)))


@register_record_type
class SsafRecord(_SoilLayerRequest):
    """``SSAF,<layer>,<save>,<outcrop1>,<outcrop2>,<layer2>,<freqstep>,<title>``: RS(layer)/RS(layer2)."""
    COMMAND = "SSAF"
    FIELDS = (F("layer", int, 1, "Layer Number"), F("save", int, 0, "Save Spectral Amplification Factor", (0, 1)),
              F("outcrop1", int, 0, "Outcropping of First Layer", (0, 1)),
              F("outcrop2", int, 0, "Outcropping of Second Layer", (0, 1)),
              F("layer2", int, 0, "Second Layer Number"), F("freqstep", float, 0.0, "Frequency Step"),
              F("title", str, "", "Title"))

    def problems(self, shared: bool = True) -> List[Problem]:
        if self.save and self.freqstep <= 0:
            return [Problem("Error", 110, f"layer {self.layer}: <freqstep> = {fmt_num(self.freqstep)}", i=self.layer)]
        return []


@register_record_type
class SfouRecord(_SoilLayerRequest):
    """``SFOU,<layer>,<out>,<save>,<outcrop>,<smooth>,<nrval>``: parsed only ("not usable in this version")."""
    COMMAND = "SFOU"
    FIELDS = (F("layer", int, 1, "Layer Number"), F("out", int, 0, "Compute Fourier Spectrum", (0, 1)),
              F("save", int, 0, "Save to File", (0, 1)), F("outcrop", int, 0, "Outcropping", (0, 1)),
              F("smooth", int, 0, "Nr. of Smoothings"), F("nrval", int, 0, "Nr. of Values to be Saved"))

    def problems(self, shared: bool = True) -> List[Problem]:
        out: List[Problem] = []
        if self.smooth < 0:
            out.append(Problem("Error", 111, f"layer {self.layer}: <smooth> = {self.smooth}", i=self.layer))
        if self.nrval < 0:
            out.append(Problem("Error", 112, f"layer {self.layer}: <nrval> = {self.nrval}", i=self.layer))
        return out


# ======================================================================================
# SITE tab (spec 07 sections 9.2.31, 9.2.44, 9.2.45; 05a section 7)
# ======================================================================================
@register_record_type
class SiteRecord(OptionRecord):
    """``SITE,<opmode>,<mode1>,<fstep>,<nl>,<hs>,<mode2>,<wopt>,<freq1>,<freq2>,<cl>,<cm>,<delt>,<nft>,<freq>``.

    SITE owns the shared time step, NFFT, frequency step, frequency set and control layer
    (spec 07 section 3).  ``<freq1>``/``<freq2>`` are frequency numbers (D-SIT-05); the dialog
    default of ``<freq2>`` is NFFT/2 (2048 for the default NFFT 4096).
    """
    COMMAND = "SITE"
    TAB = "SITE"
    FIELDS = (F("opmode", int, 0, "Operation Mode (0 solution, 1 data check)", (0, 1)),
              F("mode1", int, 1, "Mode 1", (0, 1)),
              F("fstep", float, 0.0, "Frequency Step", shared="fstep"),
              F("nl", int, 20, "Number of Generated Layers"),
              F("hs", int, 0, "Halfspace Layer"),
              F("mode2", int, 1, "Mode 2", (0, 1)),
              F("wopt", int, 0, "R-, SV- and P-Waves (0) / SH- and L-Waves (1)", (0, 1)),
              F("freq1", int, 1, "Frequency 1"),
              F("freq2", int, 2048, "Frequency 2"),
              F("cl", int, 1, "Control Point Layer", shared="cl"),
              F("cm", int, 0, "Direction (X, Y, Z)", (0, 1, 2)),
              F("delt", float, 0.005, "Time Step Control Motion", shared="delt"),
              F("nft", int, 4096, "Nr. of Fourier Component", shared="nft"),
              F("freq", int, 1, "Frequency Set Number", shared="freq"))

    def problems(self, shared: bool = True) -> List[Problem]:
        out: List[Problem] = []
        if not self.mode1 and not self.mode2:
            out.append(("Error", 45, ""))
        if self.nl != 0 and not 4 <= self.nl <= 20:
            out.append(("Error", 47, f"<nl> = {self.nl} (0 or 4..20)"))
        if self.freq1 <= 0:
            out.append(Problem("Error", 53, f"<freq1> = {self.freq1}", i=1))
        if self.freq2 <= 0:
            out.append(Problem("Error", 53, f"<freq2> = {self.freq2}", i=2))
        if shared:
            out += time_grid_problems(self.fstep, self.delt, self.nft, harmonic_ok=True)
        return out


def time_grid_problems(fstep: float, delt: float, nft: int, harmonic_ok: bool) -> List[Problem]:
    """Rules on the shared frequency grid (Errors 48-50, Warning 9; spec 07 section 9.2.31).

    A time-history analysis needs ``delt > 0`` and ``nft > 0``; a single-harmonic analysis
    (``harmonic_ok`` and ``fstep > 0``) needs only the frequency step.
    """
    out: List[Problem] = []
    if fstep < 0:
        out.append(("Error", 48, f"<fstep> = {fmt_num(fstep)}"))
    need_th = not (harmonic_ok and fstep > 0)
    if delt < 0 or (need_th and delt == 0):
        out.append(("Error", 49, f"<delt> = {fmt_num(delt)}"))
    if nft < 0 or (need_th and nft == 0):
        out.append(("Error", 50, f"<nft> = {nft}"))
    elif nft > 0 and not is_power_of_two(nft):
        out.append(("Warning", 9, f"NFFT = {nft}"))
    return out


@register_record_type
class WaveRecord(OptionRecord):
    """``WAVE,<type>,<opt>,<ratio1>,<ratio2>,<angle>``: wave field (1 R, 2 SV, 3 P, 4 SH, 5 L).

    The dialog default of a wave page is "field on, ratios 1 and 1, angle 0" (requirements
    section 5.4); a model without WAVE commands gets the vertical SV field (SH for ``<wopt>`` = 1)
    at AFWRITE.
    """
    COMMAND = "WAVE"
    KIND = "indexed"
    KEY = ("type",)
    TAB = "SITE"
    FIELDS = (F("type", int, 2, "Wave", (1, 2, 3, 4, 5)),
              F("opt", int, 1, "Wave Field (R: 1 shortest wavelength, 2 least decay)", (0, 1, 2)),
              F("ratio1", float, 1.0, "Wave Ratio 1"), F("ratio2", float, 1.0, "Wave Ratio 2"),
              F("angle", float, 0.0, "Incident Angle"))

    def problems(self, shared: bool = True) -> List[Problem]:
        out: List[Problem] = []
        if self.opt == 0:
            return out
        if not 0.0 <= self.angle < 360.0:
            out.append(Problem("Error", 52, f"wave {self.type}: angle {fmt_num(self.angle)} (D-CHK-03: [0, 360))",
                               w=self.type))
        for i, r in ((1, self.ratio1), (2, self.ratio2)):
            if not 0.0 < r <= 1.0:
                out.append(Problem("Error", 54, f"wave {self.type}, frequency {i}: ratio {fmt_num(r)}", w=self.type, i=i))
        return out


@register_record_type
class SitexRecord(OptionRecord):
    """``SITEX,<soilmode>``: 0 Linear Soil (L table) / 1 Non-Linear Soil (FILE88), D-SIT-06."""
    COMMAND = "SITEX"
    TAB = "SITE"
    XCOMMAND = True
    FIELDS = (F("soilmode", int, 0, "Linear Soil (0) / Non-Linear Soil (1)", (0, 1)),)


# ======================================================================================
# POINT, HOUSE, FORCE, ANALYS tabs (spec 07 sections 9.2.5, 9.2.14-9.2.20, 9.2.24, 9.2.39, 9.2.46)
# ======================================================================================
@register_record_type
class PointRecord(OptionRecord):
    """``POINT,<opmode>,<layer>,<rad>``: embedment layers (0 surface) and central-zone radius.

    The radius has no default (requirements section 5.4: "—"); a model without it gets Error 57.
    """
    COMMAND = "POINT"
    TAB = "POINT"
    FIELDS = (F("opmode", int, 0, "Operation Mode (0 solution, 1 data check)", (0, 1)),
              F("layer", int, 0, "Number of Embedment Soil Layers"),
              F("rad", float, 0.0, "Point Load Central Zone Radius"))

    def problems(self, shared: bool = True) -> List[Problem]:
        out: List[Problem] = []
        if self.layer < 0:
            out.append(("Error", 56, f"<layer> = {self.layer}"))
        if self.rad <= 0:
            out.append(("Error", 57, f"<rad> = {fmt_num(self.rad)}"))
        return out


@register_record_type
class HouseRecord(OptionRecord):
    """``HOUSE,<gravity>,<gelev>,<opmode>,<dim>,<imp>,<coh>,<wpass>,<me>,<cmplxspec>``.

    HOUSE owns the SSI gravity and the ground elevation (GRAVITY and GROUNDELEV set arguments 1
    and 2 only) and the coherency / wave passage / multiple excitation switches shown on the ANALYS
    tab (spec 07 section 3).
    """
    COMMAND = "HOUSE"
    TAB = "HOUSE"
    FIELDS = (F("gravity", float, 32.2, "Acceleration of Gravity", shared="gravity"),
              F("gelev", float, 0.0, "Ground Elevation", shared="gelev"),
              F("opmode", int, 0, "Operation Mode (0 solution, 1 data check)", (0, 1)),
              F("dim", int, 2, "Dimension (1 2D, 2 3D)", (1, 2)),
              F("imp", int, 0, "Flexible Volume Method (0 FV, 1 FFV, 2 FI)", (0, 1, 2)),
              F("coh", int, 0, "Soil Motion (0 Coherent, 1 Incoherent)", (0, 1), shared="coh"),
              F("wpass", int, 0, "Use Wave Passage", (0, 1), shared="wpass"),
              F("me", int, 0, "Use Multiple Excitation", (0, 1), shared="me"),
              F("cmplxspec", int, 0, "Use Complex Spectral Amp.", (0, 1)))

    def problems(self, shared: bool = True) -> List[Problem]:
        if shared and self.gravity <= 0:
            return [("Error", 1, f"HOUSE <gravity> = {fmt_num(self.gravity)}")]
        return []


@register_record_type
class HousexRecord(OptionRecord):
    """``HOUSEX,<optimize>,<supmode>,<nsim>,<nlssi>,<ansys>``: HOUSE dialog fields without an argument."""
    COMMAND = "HOUSEX"
    TAB = "HOUSE"
    XCOMMAND = True
    FIELDS = (F("optimize", int, 0, "Optimize Model", (0, 1)),
              F("supmode", int, 0, "Superposition Mode (0 Linear, 1 Quadratic)", (0, 1)),
              F("nsim", int, 1, "Number of Simulations"),
              F("nlssi", int, 0, "Non-Linear SSI", (0, 1)),
              F("ansys", int, 0, "Ansys Model Input", (0, 1)))


@register_record_type
class IncohRecord(OptionRecord):
    """``INCOH,<gammax>,<gammay>,<gammaz>,<alpha>,<ngp>,<ipr>,<nmodes>,<met>,<HSeed>,<VSeed>,<RandPhz>`` (P1)."""
    COMMAND = "INCOH"
    TAB = "HOUSE"
    FIELDS = (F("gammax", float, 0.1, "Coherence Parameter X Dir"),
              F("gammay", float, 0.1, "Coherence Parameter Y Dir"),
              F("gammaz", float, 0.2, "Coherence Parameter Z Dir"),
              F("alpha", float, 0.5, "Alpha Directionality Factor (Vs for model 1)"),
              F("ngp", int, 1, "Number of Embedded Layers"),
              F("ipr", int, 0, "Print Coherency Matrix", (0, 1)),
              F("nmodes", int, 0, "Number of Incoh. Modes"),
              F("met", int, 0, "Units (0 British, 1 SI)", (0, 1)),
              F("hseed", int, 0, "Horizontal SEED"),
              F("vseed", int, 0, "Vertical SEED"),
              F("randphz", float, 0.0, "Random Phase Angle"))

    def problems(self, shared: bool = True) -> List[Problem]:
        out: List[Problem] = []
        for name in ("gammax", "gammay", "gammaz"):
            if self.get(name) < 0.1:
                out.append(("Error", 58, f"<{name}> = {fmt_num(self.get(name))} (>= 0.1)"))
        return out


@register_record_type
class WpassRecord(OptionRecord):
    """``WPASS,<appv>,<ang>,<cohf>``: apparent velocity and angle of Line D; unlagged coherency model 1-7."""
    COMMAND = "WPASS"
    TAB = "HOUSE"
    FIELDS = (F("appv", float, 1.0e9, "Apparent Velocity for Line D"),
              F("ang", float, 0.0, "Angle Line D with X Axis"),
              F("cohf", int, 1, "Unlagged Coherency Model", (1, 2, 3, 4, 5, 6, 7)))

    def problems(self, shared: bool = True) -> List[Problem]:
        out: List[Problem] = []
        if self.appv <= 0:
            out.append(("Error", 113, f"<appv> = {fmt_num(self.appv)}"))
        if not 1 <= self.cohf <= 7:
            out.append(("Error", 114, f"<cohf> = {self.cohf} (unlagged coherency model 1..7, D-INC-01)"))
        return out


@register_record_type
class MeRecord(OptionRecord):
    """``ME,<no>,<nfirst>,<nlast>,<xc>,<yc>,<zc>``: multiple-excitation zone (P1; control point unused)."""
    COMMAND = "ME"
    KIND = "indexed"
    KEY = ("no",)
    TAB = "HOUSE"
    FIELDS = (F("no", int, 1, "Input Motion Number"), F("nfirst", int, 0, "First Foundation Node"),
              F("nlast", int, 0, "Last Foundation Node"), F("xc", float, 0.0, "X Coord."),
              F("yc", float, 0.0, "Y Coord."), F("zc", float, 0.0, "Z Coord."))


@register_record_type
class SymmRecord(OptionRecord):
    """``SYMM,<no>,[<type>],[<node1>],[<node2>],[<node3>]``: symmetry (0) / antisymmetry (1) plane (P1)."""
    COMMAND = "SYMM"
    KIND = "indexed"
    KEY = ("no",)
    TAB = ""
    FIELDS = (F("no", int, 1, "Plane / line number", (1, 2)), F("type", int, 0, "Symmetry (0) / Anti-symmetry (1)", (0, 1)),
              F("node1", int, 0, "Node 1"), F("node2", int, 0, "Node 2"), F("node3", int, 0, "Node 3"))


@register_record_type
class ForceRecord(OptionRecord):
    """``FORCE,<opmode>``: the other FORCE fields are shared (gravity, df, delt, NFFT, set)."""
    COMMAND = "FORCE"
    TAB = "FORCE"
    FIELDS = (F("opmode", int, 0, "Operation Mode (0 solution, 1 data check)", (0, 1)),)


@register_record_type
class AnalysRecord(OptionRecord):
    """``ANALYS,<opmode>,<type>,<mode>,<save>,<prnt>,<fopt>,<ang>,<xc>,<yc>,<zc>,<impe>,[simul]``."""
    COMMAND = "ANALYS"
    TAB = "ANALYS"
    OPTIONAL_FROM = 12
    FIELDS = (F("opmode", int, 0, "Operation Mode (0 solution, 1 data check)", (0, 1)),
              F("type", int, 0, "Type of Analysis (0 Seismic, 1 Foundation Vibration)", (0, 1), shared="type"),
              F("mode", int, 0, "Mode of Analysis (0 Initiation, 1 New Structure, 2 New Seismic Env., "
                                "3 New Dynamic Loading)", (0, 1, 2, 3)),
              F("save", int, 0, "Save Restart Files", (0, 1)),
              F("prnt", int, 1, "Print Amplitude Only", (0, 1)),
              F("fopt", int, 0, "Take Frequency Numbers from File1/File9", (0, 1)),
              F("ang", float, 0.0, "Coordinate Transformation Angle"),
              F("xc", float, 0.0, "X-Coordinate of Control Point"),
              F("yc", float, 0.0, "Y-Coordinate of Control Point"),
              F("zc", float, 0.0, "Z-Coordinate of Control Point"),
              F("impe", int, 0, "Global Impedance (0 None, 1 Diagonal, 2 Full 6x6)", (0, 1, 2)),
              F("simul", int, 0, "Simultaneous Cases"))

    def problems(self, shared: bool = True) -> List[Problem]:
        if not 0.0 <= self.ang < 360.0:
            return [("Error", 63, f"<ang> = {fmt_num(self.ang)} (D-CHK-03: [0, 360))")]
        return []


@register_record_type
class AnalysxRecord(OptionRecord):
    """``ANALYSX,<ffm>,<delrst>``: Free-Field Load (0) / Free-Field Motion (1); delete restart files."""
    COMMAND = "ANALYSX"
    TAB = "ANALYS"
    XCOMMAND = True
    FIELDS = (F("ffm", int, 0, "Free-Field Load (0) / Free-Field Motion (1)", (0, 1)),
              F("delrst", int, 0, "Delete Restart Files", (0, 1)))


# ======================================================================================
# MOTION, STRESS, RELDISP tabs (spec 07 sections 9.2.22, 9.2.38, 9.2.48; 05c B; 05d)
# ======================================================================================
def history_problems(mult: float, maxval: float, rec1: int, rec2: int) -> List[Problem]:
    """Errors 74-78 of the control-motion scaling and record selection (MOTION, STRESS, RELDISP, SOIL)."""
    out: List[Problem] = []
    if rec1 < 0:
        out.append(("Error", 74, f"<rec1> = {rec1}"))
    if rec2 < 0:
        out.append(("Error", 75, f"<rec2> = {rec2}"))
    if rec2 > 0 and rec1 > rec2:
        out.append(("Error", 76, f"<rec1> = {rec1} > <rec2> = {rec2}"))
    if mult == 0 and maxval == 0:
        out.append(("Error", 77, ""))
    elif mult != 0 and maxval != 0:
        out.append(("Error", 78, f"mult = {fmt_num(mult)}, max = {fmt_num(maxval)}"))
    return out


@register_record_type
class MotionRecord(OptionRecord):
    """``MOTION,<opmode>,<out>,<step>,<dur>,<res>,<freq1>,<freq2>,<fstep>,<mult>,<max>,<rec1>,<rec2>,
    <fopt>,<bl>,<smo>,<cplx>,<cnvrt>,<pzadj>,<interp>`` (19 arguments).

    MOTION owns the control-motion scaling and record selection shared by STRESS and RELDISP
    (spec 07 section 3).  New-model defaults: RS 0.1-100 Hz with 301 points, factor 1, records
    1..last, Save Complex TF on (D-UI-03, RELDISP needs it), interpolation option 1.
    """
    COMMAND = "MOTION"
    TAB = "MOTION"
    FIELDS = (F("opmode", int, 0, "Operation Mode (0 solution, 1 data check)", (0, 1)),
              F("out", int, 0, "Output Only Transfer Functions", (0, 1)),
              F("step", int, 0, "Output time-history print step"),
              F("dur", float, 0.0, "Total Duration to be Plotted"),
              F("res", int, 0, "(not used)"),
              F("freq1", float, 0.1, "First Frequency"),
              F("freq2", float, 100.0, "Last Frequency"),
              F("fstep", int, 301, "Total Number of Freq. Steps"),
              F("mult", float, 1.0, "Multiplication Factor", shared="mult"),
              F("max", float, 0.0, "Max Value for Time History", shared="max"),
              F("rec1", int, 1, "First Record", shared="rec1"),
              F("rec2", int, 0, "Last Record (0 = last)", shared="rec2"),
              F("fopt", int, 0, "File Contains Pairs Time Step - Accel.", (0, 1), shared="fopt"),
              F("bl", int, 0, "Baseline Correction (0 No, 1 With)", (0, 1)),
              F("smo", float, 0.0, "Smoothing Parameter"),
              F("cplx", int, 1, "Save Complex Transfer Functions", (0, 1)),
              F("cnvrt", int, 0, "Convert TH to RS: Select External Files", (0, 1)),
              F("pzadj", int, 0, "Phase Adjustment", (0, 1)),
              F("interp", int, 1, "Interpolation Option", (0, 1, 2, 3, 4, 5, 6)))

    def problems(self, shared: bool = True) -> List[Problem]:
        out: List[Problem] = []
        if self.step < 0:
            out.append(("Error", 67, f"<step> = {self.step}"))
        if self.dur < 0:
            out.append(("Error", 68, f"<dur> = {fmt_num(self.dur)}"))
        if self.freq1 < 0:
            out.append(("Error", 69, f"<freq1> = {fmt_num(self.freq1)}"))
        if self.freq2 < 0:
            out.append(("Error", 70, f"<freq2> = {fmt_num(self.freq2)}"))
        if self.fstep < 0:
            out.append(("Error", 71, f"<fstep> = {self.fstep}"))
        if shared:
            out += history_problems(self.mult, self.get("max"), self.rec1, self.rec2)
        return out


@register_record_type
class MotionxRecord(OptionRecord):
    """``MOTIONX,<f1213>,<resp>,<srss>,<savetf>,<saveacc>,<savers>,<saverot>,<rsttf>,<rstacc>,<rstrs>``."""
    COMMAND = "MOTIONX"
    TAB = "MOTION"
    XCOMMAND = True
    FIELDS = (F("f1213", int, 0, "Save FILE 12 or FILE 13 (0 none, 1 FILE13, 2 FILE12)", (0, 1, 2)),
              F("resp", int, 2, "Response type for Foundation Vibration (0 Disp., 1 Vel., 2 Acc.)", (0, 1, 2)),
              F("srss", int, 0, "Incoherent SSI: use SRSS (SRSSTF.txt)", (0, 1)),
              F("savetf", int, 0, "Save TF in All Points", (0, 1)),
              F("saveacc", int, 0, "Save ACC in All Points", (0, 1)),
              F("savers", int, 0, "Save RS in All Points", (0, 1)),
              F("saverot", int, 0, "Save Rotation for Ansys", (0, 1)),
              F("rsttf", int, 0, "Restart for TF", (0, 1)),
              F("rstacc", int, 0, "Restart for ACC", (0, 1)),
              F("rstrs", int, 0, "Restart for RS", (0, 1)))


@register_record_type
class StressRecord(OptionRecord):
    """``STRESS,<opmode>,<iter>,<save>,<itran>,<interopt>``."""
    COMMAND = "STRESS"
    TAB = "STRESS"
    FIELDS = (F("opmode", int, 0, "Operation Mode (0 solution, 1 data check)", (0, 1)),
              F("iter", int, 0, "Auto Computation of Strains in Soil El.", (0, 1)),
              F("save", int, 1, "Save Stress Time Histories", (0, 1)),
              F("itran", int, 0, "Output Transfer Function", (0, 1)),
              F("interopt", int, 1, "Interpolation Option", (0, 1, 2, 3, 4, 5, 6)))


@register_record_type
class StressxRecord(OptionRecord):
    """``STRESSX,<pzadj>,<smo>,<skip>,<savemax>,<saveth>,<rstns>,<rstsp>``."""
    COMMAND = "STRESSX"
    TAB = "STRESS"
    XCOMMAND = True
    FIELDS = (F("pzadj", int, 0, "Phase Adjustment", (0, 1)),
              F("smo", float, 0.0, "Smoothing Option"),
              F("skip", int, 0, "Skip Time History Steps"),
              F("savemax", int, 0, "Save Max Value", (0, 1)),
              F("saveth", int, 0, "Save Time History", (0, 1)),
              F("rstns", int, 0, "Restart for Nodal Stress Contours", (0, 1)),
              F("rstsp", int, 0, "Restart for Soil Pressure Contours", (0, 1)))


@register_record_type
class ReldRecord(OptionRecord):
    """``RELD,<RelDisOutput>,<RelDispSAll>,<RelDispNumFiles>`` (D-RDP-03: bit 0 complex TFD, bit 1 use TFU)."""
    COMMAND = "RELD"
    TAB = "RELDISP"
    FIELDS = (F("reldisoutput", int, 1, "Save Rel Disp Complex TF"),
              F("reldispsall", int, 0, "Save Relative Displacement in All Nodes", (0, 1)),
              F("reldispnumfiles", int, 0, "Number of output nodes (written = RDND count)"))


@register_record_type
class ReldxRecord(OptionRecord):
    """``RELDX,<saverot>,<rstframes>`` (P1): Save Rotations for ANSYS; Restart for Frame Generation."""
    COMMAND = "RELDX"
    TAB = "RELDISP"
    XCOMMAND = True
    FIELDS = (F("saverot", int, 0, "Save Rotations for ANSYS V11.0", (0, 1)),
              F("rstframes", int, 0, "Restart For Frame Generation", (0, 1)))


# ======================================================================================
# AFWRITE tab and SASSI-EDU extensions
# ======================================================================================
#: modules of the AOPT flags, in argument order (spec 07 section 9.2.6); None = DEP slot
AOPT_MODULES: Tuple[Optional[str], ...] = ("EQUAKE", "SOIL", None, "SITE", "POINT", "HOUSE", None, "FORCE",
                                           "ANALYS", "COMBIN", "MOTION", "STRESS", "RELDISP", "NONLINEAR")


@register_record_type
class AoptRecord(OptionRecord):
    """``AOPT,<EQUAKE>,<SOIL>,<DEP>,<SITE>,<POINT>,<HOUSE>,<DEP>,<FORCE>,<ANALYS>,<COMBIN>,<MOTION>,
    <STRESS>,<RELDISP>,<PANEL>``: modules processed by CHECK and AFWRITE.

    Default (D-AFW-05): SITE, POINT, HOUSE, ANALYS and MOTION enabled.  The DEP slots (LIQUEF and
    PINT in the dialog, not in V3) must be 0.  ``<PANEL>`` is the NONLINEAR module (P2).
    """
    COMMAND = "AOPT"
    TAB = "AFWRITE"
    FIELDS = (F("equake", int, 0, "EQUAKE", (0, 1)), F("soil", int, 0, "SOIL", (0, 1)),
              F("dep1", int, 0, "LIQUEF (not in this version)", (0,)), F("site", int, 1, "SITE", (0, 1)),
              F("point", int, 1, "POINT", (0, 1)), F("house", int, 1, "HOUSE", (0, 1)),
              F("dep2", int, 0, "PINT (not in this version)", (0,)), F("force", int, 0, "FORCE", (0, 1)),
              F("analys", int, 1, "ANALYS", (0, 1)), F("combin", int, 0, "COMBIN", (0, 1)),
              F("motion", int, 1, "MOTION", (0, 1)), F("stress", int, 0, "STRESS", (0, 1)),
              F("reldisp", int, 0, "RELDISP", (0, 1)), F("panel", int, 0, "NONLINEAR", (0, 1)))

    def enabled(self) -> List[str]:
        """Enabled modules in AOPT (run) order; a non-zero flag includes the module."""
        out = []
        for i, mod in enumerate(AOPT_MODULES, start=1):
            if mod is not None and self.integer(i, self.FIELDS[i - 1].default) not in (0, None):
                out.append(mod)
        return out


@register_record_type
class SoilxRecord(OptionRecord):
    """``SOILX,<indir>,<mult>,<max>,[<cl>],[<file>]``: SOIL input direction, SOIL's own scaling
    (exactly one of mult/max non-zero, D-SOL-11) and optional SOIL-only control layer and input file
    (0 / blank = SITE ``<cl>`` / THFILE, D-SOL-13)."""
    COMMAND = "SOILX"
    TAB = "SOIL"
    XCOMMAND = True
    FIELDS = (F("indir", int, 0, "Input Direction (0 horizontal, 1 vertical)", (0, 1)),
              F("mult", float, 0.0, "Multiplication Factor"),
              F("max", float, 0.1, "Max Value for Time History (g)"),
              F("cl", int, 0, "Control Point Layer (0 = SITE <cl>)"),
              F("file", str, "", "File (blank = THFILE)"))

    def problems(self, shared: bool = True) -> List[Problem]:
        return [p for p in history_problems(self.mult, self.get("max"), 1, 0) if p[1] in (77, 78)]


@register_record_type
class CmodformRecord(OptionRecord):
    """``CMODFORM,<form>``: complex modulus 0 ``1-2b^2+2ib sqrt(1-b^2)`` (default) / 1 ``1+2ib`` (D-CNV-03)."""
    COMMAND = "CMODFORM"
    TAB = ""
    FIELDS = (F("form", int, 0, "Complex-modulus form", (0, 1)),)


@register_record_type
class EduoptRecord(OptionRecord):
    """``EDUOPT,<key>,<value>``: algorithm switches and tolerances of requirements section 7."""
    COMMAND = "EDUOPT"
    KIND = "indexed"
    KEY = ("key",)
    TAB = ""
    FIELDS = (F("key", str, "", "Key"), F("value", str, "", "Value"))


#: EDUOPT keys of requirements section 7: key -> (default, admissible values or a type, decision)
EDUOPT_KEYS: Dict[str, Tuple[str, Any, str]] = {
    "HSLAW": ("UNIFORM", ("GEOMETRIC", "UNIFORM", "LINEAR"), "D-SIT-02"),
    "HSMAXDEPTH": ("0", float, "D-SIT-02"),
    "NFFTROUND": ("NEAREST", ("NEAREST", "UP"), "D-CNV-10"),
    "LIMITS": ("ACS", ("ACS", "PREP", "UNLIMITED"), "D-GEN-07"),
    "GEOMTOL": ("0", float, "D-GEN-06"),
    "SOILCUTOFF": ("0", float, "D-SOL-08"),
    "SOILTOL": ("0", float, "D-SOL-04"),
    "VPPOLICY": ("NU", ("NU", "VP"), "D-SOL-06"),
    "INCOHSIGN": ("ADJUST", ("ADJUST", "RAW"), "D-INC-03"),
    "INCOHMERGE": ("0", ("0", "1"), "D-INC-09"),
    "MEMLIMIT": ("0.8", float, "D-ANL-10"),
    "COMBINDUP": ("ERROR", ("ERROR", "PREFER82"), "D-CMB-01"),
    "SRSSORDER": ("INTERPFIRST", ("INTERPFIRST", "COMBINEFIRST"), "D-MOT-06"),
    "TFFILE": ("FILE8", str, "D-MOT-11"),
    "STRAINOUT": ("0", ("0", "1"), "D-STR-04"),
    "BROADENGRID": ("AUGMENTED", ("AUGMENTED", "ACS"), "D-LIN-01"),
    "BRIDGE": ("VALLEY", ("VALLEY", "AMPLITUDE"), "D-LIN-01"),
    "SECTBEND": ("1", ("0", "1"), "D-SEC-03"),
    "ANSYSMODERN": ("0", ("0", "1"), "D-ANS-06"),
    "NONEXT": ("0", ("0", "1"), "D-NON-05"),
    "SHEARFORCEARGS": ("0", ("0", "1"), "D-NON-10"),
    "DEFAULTS": ("ON", ("ON", "OFF"), "D-W5-12"),
}


def eduopt(model, key: str, default: Optional[str] = None) -> str:
    """Value of ``EDUOPT,<key>`` (upper-cased keyword values; the documented default when not set)."""
    key = key.upper()
    rec = model.options.entry("EDUOPT", key)
    if rec is not None and rec.arg(2) is not None:
        val = rec.arg(2)
        spec = EDUOPT_KEYS.get(key)
        return val.upper() if spec is not None and isinstance(spec[1], tuple) else val
    if default is not None:
        return default
    spec = EDUOPT_KEYS.get(key)
    return spec[0] if spec else ""


def eduopt_float(model, key: str, default: float = 0.0) -> float:
    try:
        return float(parse_float(eduopt(model, key, fmt_num(default))))
    except (NumberError, ValueError):
        return default


# ======================================================================================
# Option specifications (GUI registry)
# ======================================================================================
@dataclass(frozen=True)
class OptionSpec:
    """How the GUI stores one option command: kind, record class, tab, key fields."""
    command: str
    kind: str                                   # 'record' | 'indexed' | 'string' | 'list' | 'request'
    record: Optional[Type[OptionRecord]]
    tab: str
    label: str = ""
    key: Tuple[str, ...] = ()


def _specs() -> Dict[str, OptionSpec]:
    out: Dict[str, OptionSpec] = {}
    for cls in (EquakeRecord, RsinRecord, RsoutRecord, AccinRecord, AccoutRecord, TpsdRecord, CorrRecord,
                SoilRecord, SoilxRecord, SproRecord, DynpRecord, SaccRecord, SrsRecord, SstrRecord, SsafRecord,
                SfouRecord, SiteRecord, SitexRecord, WaveRecord, PointRecord, HouseRecord, HousexRecord,
                IncohRecord, WpassRecord, MeRecord, SymmRecord, ForceRecord, AnalysRecord, AnalysxRecord,
                MotionRecord, MotionxRecord, StressRecord, StressxRecord, ReldRecord, ReldxRecord, AoptRecord,
                CmodformRecord, EduoptRecord):
        out[cls.COMMAND] = OptionSpec(cls.COMMAND, cls.KIND, cls, cls.TAB, (cls.__doc__ or "").split("\n")[0],
                                      cls.KEY)
    for name, tab, label in (("EQTIT", "EQUAKE", "Spectra Title"), ("THFILE", "MOTION", "File"),
                             ("THTIT", "MOTION", "Title"), ("RELFILE", "RELDISP", "Complex TF File Name")):
        out[name] = OptionSpec(name, "string", None, tab, label)
    for name, tab, label in (("DAMP", "MOTION", "Damping Ratios"), ("TOPL", "SITE", "Top Layers"),
                             ("AMP", "HOUSE", "Spectral Amplification"), ("FREQ", "SITE", "Frequency set")):
        out[name] = OptionSpec(name, "list", None, tab, label)
    for name, tab, label in (("NOUT", "MOTION", "Node List"), ("EOUT", "STRESS", "Element Output Data"),
                             ("RDND", "RELDISP", "Nodal Output")):
        out[name] = OptionSpec(name, "request", None, tab, label)
    return out


OPTION_SPECS: Dict[str, OptionSpec] = _specs()

#: Options > Analysis tab order (requirements section 5.4)
TABS = ("EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "MOTION", "STRESS", "RELDISP",
        "NONLINEAR", "AFWRITE")


def tab_specs(tab: str) -> List[OptionSpec]:
    """Option commands stored by one Options > Analysis tab (shared variables excluded)."""
    return [s for s in OPTION_SPECS.values() if s.tab == tab.upper()]


def tab_snapshot(model, tab: str) -> Dict[str, Any]:
    """Current values of a dialog tab for the GUI (``GET /api/options/<MODULE>``).

    Record setters give their dataclass view (defaults when never set), indexed setters a list of
    dataclass views in key order, string setters their text, list commands their values; the
    shared variables shown on the tab are added under ``'shared'``.
    """
    out: Dict[str, Any] = {}
    for s in tab_specs(tab):
        if s.kind == "record":
            out[s.command] = options_from_record(get_record(model, s.command))
        elif s.kind == "indexed":
            out[s.command] = [options_from_record(r) for _, r in get_entries(model, s.command)]
        elif s.kind == "string":
            out[s.command] = model.options.string(s.command)
    out["shared"] = {name: shared_value(model, name) for name, v in SHARED_VARIABLES.items()
                     if tab.upper() in v.tabs}
    return out


def record_class(command: str) -> Type[OptionRecord]:
    spec = OPTION_SPECS.get(command.upper())
    if spec is None or spec.record is None:
        raise KeyError(f"{command.upper()} has no typed option record")
    return spec.record


def typed(rec: Optional[Record]) -> Optional[Record]:
    """The typed record of ``rec`` (generic records stored before the types were registered)."""
    if rec is None:
        return None
    cls = OPTION_SPECS.get(rec.command)
    if cls is not None and cls.record is not None and not isinstance(rec, cls.record):
        return make_record(rec.command, list(rec.values), rec.legacy)
    return rec


def get_record(model, command: str) -> OptionRecord:
    """Stored record of a record setter (typed), or a fresh record with the new-model defaults."""
    command = command.upper()
    rec = typed(model.options.record(command))
    if rec is None:
        rec = record_class(command)(command)
    return rec  # type: ignore[return-value]


def get_entries(model, command: str) -> List[Tuple[Any, OptionRecord]]:
    """Entries of an indexed setter as typed records, in key order."""
    return [(k, typed(r)) for k, r in model.options.entries(command.upper())]  # type: ignore[misc]


def aopt_modules(model) -> List[str]:
    """Modules enabled by AOPT (default D-AFW-05), in run order."""
    return get_record(model, "AOPT").enabled()  # type: ignore[attr-defined]


# ======================================================================================
# Shared analysis variables (spec 07 section 3): one storage location each
# ======================================================================================
@dataclass(frozen=True)
class SharedVariable:
    """A dialog value shown on several tabs and stored in exactly one command argument."""
    name: str
    owner: str                # owning command
    field: Optional[str]      # field of the owner record (None for string / list commands)
    label: str
    tabs: Tuple[str, ...]


SHARED_VARIABLES: Dict[str, SharedVariable] = {v.name: v for v in (
    SharedVariable("gravity", "HOUSE", "gravity", "Acceleration of Gravity (SSI)",
                   ("EQUAKE", "SITE", "HOUSE", "FORCE")),
    SharedVariable("soil_gravity", "SOIL", "grav", "Gravity Accel. (free field)", ("SOIL",)),
    SharedVariable("gelev", "HOUSE", "gelev", "Ground Elevation", ("HOUSE",)),
    SharedVariable("delt", "SITE", "delt", "Time Step", ("EQUAKE", "SOIL", "SITE", "HOUSE", "FORCE", "MOTION",
                                                         "STRESS", "RELDISP")),
    SharedVariable("nft", "SITE", "nft", "Number of Fourier Components", ("SOIL", "SITE", "HOUSE", "FORCE",
                                                                          "MOTION", "STRESS", "RELDISP")),
    SharedVariable("fstep", "SITE", "fstep", "Frequency Step", ("SITE", "FORCE")),
    SharedVariable("freq", "SITE", "freq", "Frequency Set Number", ("SITE", "HOUSE", "FORCE", "ANALYS", "STRESS",
                                                                    "RELDISP")),
    SharedVariable("cl", "SITE", "cl", "Control Point Layer", ("SITE", "SOIL")),
    SharedVariable("type", "ANALYS", "type", "Type of Analysis", ("ANALYS", "MOTION", "STRESS")),
    SharedVariable("coh", "HOUSE", "coh", "Coherent / Incoherent", ("HOUSE", "ANALYS")),
    SharedVariable("wpass", "HOUSE", "wpass", "Wave Passage Effects Included", ("HOUSE", "ANALYS")),
    SharedVariable("me", "HOUSE", "me", "Multiple Excitation", ("HOUSE", "ANALYS")),
    SharedVariable("mult", "MOTION", "mult", "Multiplication Factor", ("MOTION", "STRESS", "RELDISP")),
    SharedVariable("max", "MOTION", "max", "Max Value for Time History", ("MOTION", "STRESS", "RELDISP")),
    SharedVariable("rec1", "MOTION", "rec1", "First Record", ("MOTION", "STRESS", "RELDISP")),
    SharedVariable("rec2", "MOTION", "rec2", "Last Record", ("MOTION", "STRESS", "RELDISP")),
    SharedVariable("fopt", "MOTION", "fopt", "File Contains Pairs", ("MOTION", "STRESS", "RELDISP")),
    SharedVariable("thfile", "THFILE", None, "File", ("SOIL", "MOTION", "STRESS", "RELDISP")),
    SharedVariable("thtit", "THTIT", None, "Title", ("SOIL", "MOTION", "STRESS", "RELDISP")),
    SharedVariable("damp", "DAMP", None, "Damping Ratios", ("SOIL", "MOTION")),
)}


def shared_value(model, name: str) -> Any:
    """Current value of a shared variable (its owner's default when the owner was never given)."""
    v = SHARED_VARIABLES[name]
    if v.owner in ("THFILE", "THTIT"):
        return model.options.string(v.owner)
    if v.owner == "DAMP":
        return list(model.damp)
    return get_record(model, v.owner).get(v.field)


def shared_command(model, name: str, value: Any) -> List[str]:
    """Command lines that set shared variable ``name`` to ``value`` in its one storage location (L17).

    The owner record keeps its other arguments; gravity and ground elevation use GRAVITY and
    GROUNDELEV (which change HOUSE arguments 1 and 2 only).
    """
    v = SHARED_VARIABLES[name]
    if name == "gravity":
        return [f"GRAVITY,{_tok(value)}"]
    if name == "gelev":
        return [f"GROUNDELEV,{_tok(value)}"]
    if v.owner in ("THFILE", "THTIT"):
        return [f"{v.owner},{value}"]
    if v.owner == "DAMP":
        vals = list(value) if isinstance(value, (list, tuple)) else [value]
        out = ["DAMP,0"]
        for i in range(0, len(vals), 10):
            out.append("DAMP," + ",".join(_tok(x) for x in vals[i:i + 10]))
        return out
    rec = get_record(model, v.owner).copy()
    rec.set(v.field, value)
    return [to_command(rec)]


# ======================================================================================
# Dataclass views for the GUI: from_args / to_command
# ======================================================================================
_DATACLASSES: Dict[str, type] = {}


def options_class(command: str) -> type:
    """Dataclass of command ``command``: field names = argument names, defaults = new-model defaults,
    ``metadata['label']`` = dialog label, ``metadata['position']`` = argument number."""
    command = command.upper()
    if command in _DATACLASSES:
        return _DATACLASSES[command]
    cls = record_class(command)
    fields = []
    for i, f in enumerate(cls.FIELDS, start=1):
        md = {"label": f.doc, "position": i, "choices": getattr(f, "choices", ()), "shared": getattr(f, "shared", "")}
        fields.append((f.name, f.type, dataclasses.field(default=f.default, metadata=md)))
    dc = dataclasses.make_dataclass(f"{command}Options", fields,
                                    namespace={"COMMAND": command, "KIND": cls.KIND, "KEY": cls.KEY,
                                               "__doc__": cls.__doc__})
    _DATACLASSES[command] = dc
    return dc


def dialog_fields(command: str) -> List[Dict[str, Any]]:
    """The dialog fields of a command: name, type, default, label, choices, position, shared."""
    cls = record_class(command)
    return [dict(name=f.name, type=f.type.__name__, default=f.default, label=f.doc,
                 choices=list(getattr(f, "choices", ())), position=i, shared=getattr(f, "shared", ""))
            for i, f in enumerate(cls.FIELDS, start=1)]


def _convert(f: Field, token: Optional[str]) -> Any:
    if token is None or token == "":
        return f.default
    if f.type is str:
        return token
    if f.type is int:
        return parse_int(token)[0]
    return parse_float(token)


def options_from_record(rec: Record):
    """Typed dataclass view of a stored record (blank fields -> defaults)."""
    dc = options_class(rec.command)
    cls = record_class(rec.command)
    vals = {}
    for i, f in enumerate(cls.FIELDS, start=1):
        try:
            vals[f.name] = _convert(f, rec.arg(i))
        except NumberError:
            vals[f.name] = f.default
    return dc(**vals)


def _command_spec(name: str):
    """The command specification the interpreter uses for ``name`` (rule L17: the GUI parses and
    emits command text exactly like the console).

    :func:`sassi.prep.registry.lookup` loads the command handlers and returns the *registered*
    specification, which carries the options declared at registration (e.g. the rest-of-line text
    argument ``SOILX <file>``) and resolves abbreviations; the static catalogue is the fallback.
    """
    from .registry import CATALOGUE, lookup
    return lookup(name) or CATALOGUE.get(name.strip().upper())


def _split_command(text: str) -> Tuple[str, List[str]]:
    from .lexer import split_args, split_head
    name, rest = split_head(text)
    spec = _command_spec(name)
    lexed = split_args(rest, spec.text_from if spec else None, bool(spec and spec.paren))
    return (spec.name if spec else name.upper()), list(lexed.tokens)


def from_args(args: Union[str, Sequence[str]], command: Optional[str] = None):
    """Dataclass of a command from its text (``'SITE,0,1,...'``) or from its argument tokens.

    ``args`` is either a command line (the command name first) or a sequence of tokens; with a
    sequence, ``command`` names the command unless the first token is a command name.  Blank or
    omitted arguments take the new-model defaults; non-numeric numbers raise ``ValueError``.
    """
    if isinstance(args, str):
        name, toks = _split_command(args)
    else:
        toks = [("" if t is None else str(t)) for t in args]
        if command is None:
            if not toks:
                raise ValueError("command name required")
            name, toks = toks[0].upper(), toks[1:]
        else:
            name = command.upper()
    cls = record_class(name)
    dc = options_class(name)
    vals = {}
    for i, f in enumerate(cls.FIELDS):
        t = toks[i].strip() if i < len(toks) else ""
        try:
            vals[f.name] = _convert(f, t)
        except NumberError as exc:
            raise ValueError(f"{name} argument {i + 1} <{f.name}>: {exc}") from None
    return dc(**vals)


def _tok(v: Any) -> str:
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, (int, float)):
        if isinstance(v, float) and not math.isfinite(v):
            return repr(v)
        return fmt_num(v)
    return "" if v is None else str(v)


def record_from_options(obj) -> OptionRecord:
    """Record (token text) of a dataclass view."""
    name = obj.COMMAND
    cls = record_class(name)
    return cls(name, [_tok(getattr(obj, f.name)) for f in cls.FIELDS])  # type: ignore[return-value]


def to_command(obj) -> str:
    """Command text of a dataclass view or of a record (rule L17: the GUI submits this text).

    Every argument is written (manual-identical lines) except trailing *optional* arguments at
    their defaults (EQUAKE ``[tpsd]``, ANALYS ``[simul]``) and, for X-commands, trailing defaults.
    A record keeps the token text the user typed.
    """
    from .lexer import join_command
    if isinstance(obj, Record):
        name = obj.command
        cls = OPTION_SPECS[name].record if name in OPTION_SPECS else None
        values = obj.to_tokens()                    # the user's token text; blank = default
        defaults = [_tok(f.default) for f in cls.FIELDS] if cls is not None else []
    else:
        name = obj.COMMAND
        cls = record_class(name)
        values = [_tok(getattr(obj, f.name)) for f in cls.FIELDS]
        defaults = [_tok(f.default) for f in cls.FIELDS]
    trim_from = None
    if cls is not None:
        if cls.XCOMMAND:
            trim_from = 1
        elif cls.OPTIONAL_FROM:
            trim_from = cls.OPTIONAL_FROM
    if trim_from is not None:
        while len(values) >= trim_from and len(values) <= len(defaults) and \
                _same(values[-1], defaults[len(values) - 1]):
            values.pop()
    spec = _command_spec(name)
    text_last = bool(spec and spec.text_from and len(values) == spec.text_from)
    return join_command(name, values, text_last=text_last)


def _same(a: Optional[str], b: Optional[str]) -> bool:
    if (a or "") == (b or "") or not a:             # a blank token is the default
        return True
    try:
        return parse_float(a or "0") == parse_float(b or "0")
    except NumberError:
        return False


__all__ = ["OptField", "OptionRecord", "Problem", "problem_fmt", "OPTION_SPECS", "OptionSpec", "TABS", "SHARED_VARIABLES",
           "SharedVariable", "AOPT_MODULES", "EDUOPT_KEYS", "record_class", "typed", "get_record", "get_entries",
           "aopt_modules", "shared_value", "shared_command", "options_class", "dialog_fields", "from_args",
           "to_command", "options_from_record", "record_from_options", "eduopt", "eduopt_float",
           "tab_specs", "tab_snapshot",
           "time_grid_problems", "history_problems"]
