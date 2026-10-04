"""Nonlinear (near-field) soil SSI: the equivalent-linear SSI iterations of ACS SASSI.

Normative sources: requirements section 4.4 item 7 (HOUSE), 4.10 item 5 (STRESS), 2.4 (run sequence
"NONLINEAR SOIL SSI ITERATIONS"), 3.4.R (COMBXYZSTRAIN); decisions D-NLS-01 ... D-NLS-04 and
D-FIL-08; spec 01 section 8.1, spec 03 item 18 and section 4.5, spec 05b section 6 (the ``.pin``
file), spec 05d section 1.8; theory R1 section 6.3 (the SHAKE property update, which is the same
rule applied element by element).

What the iterations do (for the structural engineer)
----------------------------------------------------
The SASSI solution is linear (complex frequency domain).  Soil nonlinearity is approximated with
the *equivalent-linear* method of SHAKE: each soil element gets a shear modulus and a damping ratio
that are compatible with its own strain level.  There are two kinds of soil nonlinearity
(spec 03 item 18):

* **primary** (free field): the vertically propagating waves in the horizontally layered site;
  SOIL (SHAKE) iterates the layer properties, SITE uses them (SITEX,1 and FILE88);
* **secondary** (local, near field): the extra straining caused by the structure, e.g. a soft
  backfill zone next to an embedded wall, or several heavy neighbouring buildings.  The near-field
  soil is modelled with SOLID (3D) or PLANE (2D) elements that are part of the *structure* model
  (HOUSE ETYPE 1, next to the excavated free-field soil of the flexible-volume method) and their
  properties are iterated element by element.

One iteration (requirements 2.4, spec 03 section 4.5)::

    HOUSE   reads the .pin file (groups, ESF, ISTR, GFAC, DFAC, ICURVE); iteration 0 uses
            G = GFAC G_layer, beta = DFAC beta_layer of the free field; later iterations (.liq = 1)
            read FILE74 and set G = G_max (G/G_max)(gamma_eff), beta = D(gamma_eff) from the FILE73
            curve ICURVE, one internal material per element; writes FILE78 (the properties used)
    ANALYS  "New Structure" restart (deck <mode> 1, manual Mode 2): the soil impedance X_ff of the
            initiation run (COOXqqq) is reused, only the structure matrices change
    STRESS  (<iter> = 1, once per input direction) strain histories of every nonlinear element in
            the time domain, gamma(t) per ISTR, gamma_eff = ESF max|gamma(t)| -> FILE74
    COMBXYZSTRAIN  SRSS of the directional effective strains -> FILE74 (D-NLS-04)
    converged?     max |dG/G| < 2 % and max |d beta| < 0.5 % (absolute), at most 8 iterations
                   (D-NLS-02)

Strain measures (D-NLS-03; requirements 4.10 item 5)::

    ISTR 0   3D: gamma(t) = max(|g_xy|, |g_xz|, |g_yz|)    2D: |g_xz|
    ISTR 1   3D: octahedral  2/3 sqrt[(e_x-e_y)^2 + (e_y-e_z)^2 + (e_z-e_x)^2 + 1.5 (g_xy^2 + g_yz^2 + g_xz^2)]
             2D: maximum shear strain  sqrt[(e_x-e_z)^2 + g_xz^2]

(engineering shear strains ``g = 2 eps``).  Simple shear gamma gives ``gamma`` with ISTR 0 and
``sqrt(2/3) gamma`` with ISTR 1, so ISTR 0 is the measure of SHAKE (R1 section 6.2).

Curve interpolation is that of SOIL: linear in log10(strain), constant outside the tabulated range
(D-SOL-03); curve strains and damping are in percent.  When G changes, Poisson's ratio is kept
(the constrained modulus follows G) and the P-wave damping equals the shear damping -- the default
policy of SOIL for iterated sublayers (D-SOL-06).

What the near-field material and the .pin factors mean (requirements 4.4 item 7; manual: GFAC
"1.00 indicates same shear modulus as in free-field")::

    M-table material of the element   the LOW-STRAIN soil of the near field: density, Poisson's
                                      ratio and G_max = rho Vs^2, the reference of the G/G_max curve
    G_layer, beta_layer               the free field at the element: the TOPL layer (L properties,
                                      as for the excavated soil, D-ELM-12) containing the element's
                                      mid-depth (centroid of its distinct nodes; above grade: TOPL
                                      layer 1; below the last TOPL layer: the half-space row)
    iteration 0                       G = GFAC G_layer, beta_s = DFAC ds_layer, beta_p = DFAC dp_layer,
                                      Poisson's ratio and density of the element material
    iteration k >= 1                  G = G_max (G/G_max)(gamma_eff), beta_s = beta_p = D(gamma_eff)

GFAC = DFAC = 1 starts the near field from the free-field state; a soft backfill typically starts
lower (GFAC < 1).  Without a free-field layer table the element material is the reference (a
warning).  GFAC/DFAC choose the starting point only: as in SHAKE, the converged (strain-compatible)
state is set by G_max, the curves and the motion (VP-N1 starts from both ends).

Iteration state (which files belong to the analysis in progress)
---------------------------------------------------------------
FILE78 carries the FILE4 hash HOUSE wrote with its properties, a digest of the HOUSE deck and of the
.pin; FILE74 carries the FILE4 hash of the properties whose response it holds.  STRESS (<iter> = 1)
accepts only a FILE8 computed from that FILE4; COMBXYZSTRAIN only directional FILE74s of one
iteration and FILE4; :func:`status` only a FILE74 of the FILE78 run; NLSSIITER reports "already
converged" only for the analysis in progress (:func:`current_status`).  NLSSIRESET and a linear
HOUSE run retire FILE78 (and NLSSIRESET FILE74) as ``FILE78.prev`` / ``FILE74.prev``.

File formats (D-FIL-08; text, written and read here)
----------------------------------------------------
``<model>.pin`` (user input, spec 05b section 6; free format, commas or blanks, ``*``/``#`` comment
lines and trailing ``!`` comments allowed)::

    NGRP, ESF, NCURV
    IGRP, NMAT, NELEM, ISTR          (per group)
    GFAC, DFAC, ICURVE               (NMAT lines: the materials of the group in ascending number)

``<model>.liq``: ``1`` -> the next HOUSE run reads FILE74 (anything else: initial run from the .pin).

``FILE78`` (HOUSE -> STRESS) and ``FILE74`` (STRESS / COMBXYZSTRAIN -> HOUSE)::

    # comment lines
    ITERATION <k>                     iteration number (0 = initial properties of the .pin)
    ESF <esf>
    NCURV <ncurv>
    DIRECTION <X|Y|Z|XYZ ...>         FILE74: input direction(s) of the strains
    SOURCE <PIN | FILE74 ...>         FILE78: where the properties come from
    FILE4 <hash>                      FILE78: model hash of the FILE4 written with these properties;
                                      FILE74: the same hash (the strains are the response to them)
    DECK <digest>                     FILE78: digest of the HOUSE deck data (comments excluded)
    PIN <digest>                      FILE78: digest of the .pin content
    GROUP <igrp> <SOLID|PLANE> <nelem> <istr>
    <group> <element> <gamma_eff %> <G> <beta> <curve> <Gmax> <extra>

``extra`` is the material number (FILE78) or the maximum strain max|gamma(t)| in percent
(FILE74).  The first six columns are those of D-FIL-08.

``NLSOIL_CONVERGENCE.TXT``: one row per iteration k -- the change from the properties *used* in
iteration k (FILE78) to the strain-compatible properties of its response (FILE74).
"""
from __future__ import annotations

import hashlib
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..conventions import cfactor
from ..elements.base import MaterialProps
from .shake import DynamicProperty, read_file73

PathLike = Union[str, Path]

# ======================================================================================
# Constants (D-NLS-02, D-NLS-03)
# ======================================================================================
#: convergence: max |dG/G| below this many percent ...
TOL_G_PCT = 2.0
#: ... and max |d beta| below this many percentage points (absolute damping change)
TOL_BETA_PCT = 0.5
#: default iteration limit
MAX_ITERATIONS = 8
#: ISTR codes of the .pin file
STRAIN_FLAGS = {0: "maximum shear-strain component (3D: max |g_xy|, |g_xz|, |g_yz|; 2D: |g_xz|)",
                1: "octahedral shear strain (3D) / maximum shear strain (2D)"}
#: element types of nonlinear soil groups: SOLID in 3D models (HOUSE <dim> 2), PLANE in 2D (<dim> 1)
NL_TYPES = {1: "SOLID", 4: "PLANE"}
#: strain components of the FILE4 strain operators x_rec_<T>_B (engineering shear strains)
STRAIN_COMPONENTS = {"SOLID": ("EXX", "EYY", "EZZ", "GXY", "GXZ", "GYZ"), "PLANE": ("EXX", "EZZ", "GXZ")}
#: control direction <cm> -> FILE74 direction label
DIRECTIONS = ("X", "Y", "Z")
CONVERGENCE_FILE = "NLSOIL_CONVERGENCE.TXT"
FILE_KINDS = ("FILE74", "FILE78")


class NLSoilError(ValueError):
    """Invalid nonlinear-soil input (.pin, FILE74, FILE78, curves) -- the message says what to fix."""


# ======================================================================================
# The .pin file (spec 05b section 6)
# ======================================================================================
@dataclass
class PinMaterial:
    """One material line of a nonlinear group: initial factors and the FILE73 curve number.

    The factors scale the *free field* at the element (requirements 4.4 item 7): iteration 0 has
    ``G = GFAC G_layer`` and ``beta = DFAC beta_layer``; the element's M-table material holds the
    low-strain soil (G_max of the curve)."""
    gfac: float = 1.0          # initial shear-modulus factor (1.0 = the free-field G of the TOPL layer)
    dfac: float = 1.0          # initial damping factor (1.0 = the free-field damping of the TOPL layer)
    icurve: int = 1            # pair of G/G_max and D curves in FILE73 (1-based)


@dataclass
class PinGroup:
    """One nonlinear soil element group of the .pin file."""
    igrp: int
    nmat: int
    nelem: int
    istr: int
    materials: List[PinMaterial] = field(default_factory=list)


@dataclass
class PinData:
    """Content of ``<model>.pin``: effective strain factor ESF, number of curves NCURV, groups."""
    esf: float
    ncurv: int
    groups: List[PinGroup] = field(default_factory=list)

    @property
    def ngrp(self) -> int:
        return len(self.groups)

    def group(self, igrp: int) -> Optional[PinGroup]:
        for g in self.groups:
            if g.igrp == int(igrp):
                return g
        return None


_SPLIT = re.compile(r"[,\s;]+")


def _pin_records(text: str) -> List[Tuple[int, List[str]]]:
    """(line number, tokens) of the non-comment lines (``*``, ``#`` comment lines; ``!`` comments)."""
    out = []
    for k, ln in enumerate(text.splitlines(), start=1):
        s = ln.split("!", 1)[0].strip()
        if not s or s[0] in "*#":
            continue
        toks = [t for t in _SPLIT.split(s) if t]
        if toks:
            out.append((k, toks))
    return out


def _num(tok: str, what: str, line: int, integer: bool = False) -> float:
    try:
        v = float(tok.replace("D", "E").replace("d", "e"))
    except ValueError:
        raise NLSoilError(f".pin line {line}: {what} '{tok}' is not a number") from None
    if integer:
        if v != round(v):
            raise NLSoilError(f".pin line {line}: {what} = {tok} must be an integer")
        return int(round(v))
    return v


def parse_pin(text: str) -> PinData:
    """Parse the .pin text (spec 05b section 6) and validate its counts and codes.

    Raises :class:`NLSoilError` naming the line of the first problem."""
    recs = _pin_records(text)
    if not recs:
        raise NLSoilError(".pin is empty: line 1 must hold NGRP, ESF, NCURV")
    ln, t = recs[0]
    if len(t) < 3:
        raise NLSoilError(f".pin line {ln}: line 1 needs 3 items NGRP, ESF, NCURV (found {len(t)})")
    ngrp = _num(t[0], "NGRP", ln, True)
    esf = _num(t[1], "ESF", ln)
    ncurv = _num(t[2], "NCURV", ln, True)
    if ngrp < 1:
        raise NLSoilError(f".pin line {ln}: NGRP = {ngrp} must be >= 1")
    if not 0.0 < esf <= 1.0:
        raise NLSoilError(f".pin line {ln}: the effective strain factor ESF = {esf:g} must lie in (0, 1] "
                          "(typically 0.6-0.7)")
    if ncurv < 1:
        raise NLSoilError(f".pin line {ln}: NCURV = {ncurv} must be >= 1 (curves of FILE73)")
    pin = PinData(esf=float(esf), ncurv=int(ncurv))
    i = 1
    for _ in range(int(ngrp)):
        if i >= len(recs):
            raise NLSoilError(f".pin: {ngrp} groups announced on line 1, only {pin.ngrp} found")
        ln, t = recs[i]
        i += 1
        if len(t) < 4:
            raise NLSoilError(f".pin line {ln}: a group line needs 4 items IGRP, NMAT, NELEM, ISTR")
        g = PinGroup(igrp=_num(t[0], "IGRP", ln, True), nmat=_num(t[1], "NMAT", ln, True),
                     nelem=_num(t[2], "NELEM", ln, True), istr=_num(t[3], "ISTR", ln, True))
        if g.igrp < 1:
            raise NLSoilError(f".pin line {ln}: group number {g.igrp} must be >= 1")
        if pin.group(g.igrp) is not None:
            raise NLSoilError(f".pin line {ln}: group {g.igrp} is listed twice")
        if g.nmat < 1:
            raise NLSoilError(f".pin line {ln}: NMAT = {g.nmat} must be >= 1")
        if g.istr not in STRAIN_FLAGS:
            raise NLSoilError(f".pin line {ln}: ISTR = {g.istr} must be 0 (maximum shear-strain component) or "
                              "1 (octahedral / 2D maximum shear strain)")
        for _m in range(g.nmat):
            if i >= len(recs):
                raise NLSoilError(f".pin: group {g.igrp} announces {g.nmat} materials, only "
                                  f"{len(g.materials)} lines found")
            ln, t = recs[i]
            i += 1
            if len(t) < 3:
                raise NLSoilError(f".pin line {ln}: a material line needs 3 items GFAC, DFAC, ICURVE")
            m = PinMaterial(gfac=_num(t[0], "GFAC", ln), dfac=_num(t[1], "DFAC", ln),
                            icurve=_num(t[2], "ICURVE", ln, True))
            if m.gfac <= 0:
                raise NLSoilError(f".pin line {ln}: GFAC = {m.gfac:g} must be > 0")
            if m.dfac < 0:
                raise NLSoilError(f".pin line {ln}: DFAC = {m.dfac:g} must be >= 0")
            if not 1 <= m.icurve <= pin.ncurv:
                raise NLSoilError(f".pin line {ln}: ICURVE = {m.icurve} must lie in 1..NCURV = {pin.ncurv}")
            g.materials.append(m)
        pin.groups.append(g)
    if i < len(recs):
        raise NLSoilError(f".pin line {recs[i][0]}: unexpected data after the {pin.ngrp} groups announced on line 1")
    return pin


def read_pin(path: PathLike) -> PinData:
    return parse_pin(Path(path).read_text(encoding="utf-8", errors="replace"))


def format_pin(pin: PinData, comments: Sequence[str] = (), material_notes: Optional[Dict[int, List[str]]] = None) -> str:
    """The .pin text in the manual's layout (comma separated), with ``*`` comment lines first and
    optional trailing ``!`` notes per material line (``material_notes[igrp][k]``)."""
    lines = [f"* {c}" if not c.startswith("*") else c for c in comments]
    lines.append(f"{pin.ngrp}, {float(pin.esf)!r}, {pin.ncurv}")
    for g in pin.groups:
        lines.append(f"{g.igrp}, {g.nmat}, {g.nelem}, {g.istr}")
        notes = (material_notes or {}).get(g.igrp, [])
        for k, m in enumerate(g.materials):
            note = f"   ! {notes[k]}" if k < len(notes) and notes[k] else ""
            lines.append(f"{float(m.gfac)!r}, {float(m.dfac)!r}, {m.icurve}{note}")
    return "\n".join(lines) + "\n"


def write_pin(path: PathLike, pin: PinData, comments: Sequence[str] = (),
              material_notes: Optional[Dict[int, List[str]]] = None) -> Path:
    path = Path(path)
    path.write_text(format_pin(pin, comments, material_notes), encoding="utf-8")
    return path


def pin_digest(pin: PinData) -> str:
    """Digest of the .pin *content* (comments and layout ignored): recorded in FILE78 so that a later
    change of the .pin is recognised (:func:`current_status`)."""
    return hashlib.sha256(format_pin(pin).encode("utf-8")).hexdigest()[:32]


def deck_digest(path: PathLike) -> str:
    """Digest of the data lines of a module deck (D-GEN-04; blank and ``*``/``#`` comment lines are
    ignored, so the model-hash comment AFWRITE rewrites at every run does not count).  '' when the deck
    cannot be read."""
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""
    h = hashlib.sha256()
    for ln in text.splitlines():
        s = ln.strip()
        if s and s[0] not in "*#":
            h.update(s.encode("utf-8"))
            h.update(b"\n")
    return h.hexdigest()[:32]


def _fmt(v: float) -> str:
    """Shortest exact text of a float (repr), integers without a trailing '.0'."""
    v = float(v)
    return str(int(v)) if v == int(v) and abs(v) < 1e15 else repr(v)


# ======================================================================================
# .liq (iteration flag)
# ======================================================================================
def read_liq(path: PathLike) -> int:
    """1 when ``<model>.liq`` exists and its first item is 1 (HOUSE then reads FILE74), else 0."""
    p = Path(path)
    if not p.is_file():
        return 0
    toks = p.read_text(encoding="utf-8", errors="replace").split()
    try:
        return 1 if toks and int(float(toks[0])) == 1 else 0
    except ValueError:
        return 0


def write_liq(path: PathLike, value: int = 1) -> Path:
    p = Path(path)
    p.write_text(f"{int(value)}\n", encoding="utf-8")
    return p


RETIRED_SUFFIX = ".prev"


def retire(workdir: PathLike, names: Sequence[str] = ("FILE78", "FILE74")) -> List[str]:
    """Rename iteration files that no longer belong to the analysis in progress to ``<name>.prev``
    (an older ``.prev`` is replaced), so that no module or command picks them up; returns the new
    names.  Used by NLSSIRESET (FILE78 and FILE74) and by a linear HOUSE run (FILE78)."""
    out = []
    for nm in names:
        p = Path(workdir) / nm
        if p.is_file():
            q = p.with_name(p.name + RETIRED_SUFFIX)
            os.replace(p, q)
            out.append(q.name)
    return out


def same_file4(a: str, b: str) -> bool:
    """True unless both FILE4 hashes are known and differ (hand-made files carry none)."""
    return not (a and b and a != b)


# ======================================================================================
# FILE74 / FILE78 (D-FIL-08)
# ======================================================================================
@dataclass
class NLRow:
    """One element row of FILE74/FILE78 (strains in percent, damping as a fraction)."""
    group: int
    element: int
    gamma_eff: float           # effective strain % (FILE78: the strain the properties come from; 0 initially)
    G: float                   # shear modulus G (real, undamped part of G*)
    beta: float                # shear damping ratio
    curve: int                 # FILE73 curve number ICURVE
    gmax: float                # low-strain modulus G_max of the element material
    extra: float = 0.0         # FILE78: material number; FILE74: max|gamma(t)| %

    @property
    def key(self) -> Tuple[int, int]:
        return (self.group, self.element)

    @property
    def ratio(self) -> float:
        """G / G_max."""
        return self.G / self.gmax if self.gmax > 0 else float("nan")


@dataclass
class NLGroupInfo:
    group: int
    etype: str
    nelem: int
    istr: int


@dataclass
class NLFile:
    """Content of FILE74 or FILE78."""
    kind: str                  # 'FILE74' or 'FILE78'
    iteration: int
    esf: float
    ncurv: int = 0
    direction: str = ""
    source: str = ""
    file4: str = ""            # FILE4 hash of the properties (FILE78) / of the response (FILE74)
    model: str = ""
    groups: List[NLGroupInfo] = field(default_factory=list)
    rows: List[NLRow] = field(default_factory=list)
    deck: str = ""             # FILE78: digest of the HOUSE deck data (deck_digest)
    pin: str = ""              # FILE78: digest of the .pin content (pin_digest)

    def by_key(self) -> Dict[Tuple[int, int], NLRow]:
        return {r.key: r for r in self.rows}

    def group_info(self, g: int) -> Optional[NLGroupInfo]:
        for gi in self.groups:
            if gi.group == int(g):
                return gi
        return None

    def arrays(self) -> Dict[str, np.ndarray]:
        """Columns as arrays (group, element, gamma_eff, G, beta, curve, gmax, extra)."""
        cols = ("group", "element", "gamma_eff", "G", "beta", "curve", "gmax", "extra")
        out = {c: np.array([getattr(r, c) for r in self.rows], dtype=float) for c in cols}
        for c in ("group", "element", "curve"):
            out[c] = out[c].astype(np.int64)
        return out


_HEADERS = {
    "FILE78": "nonlinear-soil properties used by HOUSE (D-FIL-08, requirements 4.4 item 7)",
    "FILE74": "effective strains and strain-compatible properties (D-FIL-08, requirements 4.10 item 5)",
}
_EXTRA = {"FILE78": "material", "FILE74": "gamma_max_pct"}


def format_nlfile(nl: NLFile, title: str = "") -> str:
    if nl.kind not in FILE_KINDS:
        raise ValueError(f"unknown nonlinear-soil file kind {nl.kind!r}")
    lines = [f"# SASSI-EDU {nl.kind} v1: {_HEADERS[nl.kind]}"]
    if nl.model or title:
        lines.append(f"# model {nl.model}" + (f" -- {title}" if title else ""))
    lines.append("# gamma_eff in percent; G and Gmax in model stress units (real part of G*); beta = shear "
                 "damping ratio (fraction)")
    lines.append(f"ITERATION {int(nl.iteration)}")
    lines.append(f"ESF {_fmt(nl.esf)}")
    lines.append(f"NCURV {int(nl.ncurv)}")
    if nl.direction:
        lines.append(f"DIRECTION {nl.direction}")
    if nl.source:
        lines.append(f"SOURCE {nl.source}")
    if nl.file4:
        lines.append(f"FILE4 {nl.file4}")
    if nl.deck:
        lines.append(f"DECK {nl.deck}")
    if nl.pin:
        lines.append(f"PIN {nl.pin}")
    if nl.model:
        lines.append(f"MODEL {nl.model}")
    for gi in nl.groups:
        lines.append(f"GROUP {gi.group} {gi.etype} {gi.nelem} {gi.istr}")
    lines.append(f"# group element gamma_eff_pct G beta curve Gmax {_EXTRA[nl.kind]}")
    for r in nl.rows:
        ex = f"{int(round(r.extra)):d}" if nl.kind == "FILE78" else f"{r.extra:.10e}"
        lines.append(f"{r.group:6d} {r.element:7d} {r.gamma_eff:.10e} {r.G:.12e} {r.beta:.10e} {r.curve:4d} "
                     f"{r.gmax:.12e} {ex}")
    return "\n".join(lines) + "\n"


def write_nlfile(path: PathLike, nl: NLFile, title: str = "") -> Path:
    p = Path(path)
    p.write_text(format_nlfile(nl, title), encoding="utf-8")
    return p


def read_nlfile(path: PathLike, kind: Optional[str] = None) -> NLFile:
    """Read FILE74 or FILE78 (``kind`` checks the header; ``None`` accepts both)."""
    p = Path(path)
    if not p.is_file():
        raise NLSoilError(f"{p.name} not found")
    text = p.read_text(encoding="utf-8", errors="replace")
    found = ""
    m = re.search(r"SASSI-EDU (FILE7[48])", text[:400])
    if m:
        found = m.group(1)
    if kind is not None and found and found != kind:
        raise NLSoilError(f"{p.name} is a {found} file, expected {kind}")
    nl = NLFile(kind=found or (kind or "FILE74"), iteration=-1, esf=0.0)
    for k, ln in enumerate(text.splitlines(), start=1):
        s = ln.strip()
        if not s or s[0] in "#*":
            continue
        t = s.split()
        head = t[0].upper()
        try:
            if head == "ITERATION":
                nl.iteration = int(t[1])
            elif head == "ESF":
                nl.esf = float(t[1])
            elif head == "NCURV":
                nl.ncurv = int(t[1])
            elif head == "DIRECTION":
                nl.direction = " ".join(t[1:])
            elif head == "SOURCE":
                nl.source = " ".join(t[1:])
            elif head == "FILE4":
                nl.file4 = t[1] if len(t) > 1 else ""
            elif head == "DECK":
                nl.deck = t[1] if len(t) > 1 else ""
            elif head == "PIN":
                nl.pin = t[1] if len(t) > 1 else ""
            elif head == "MODEL":
                nl.model = " ".join(t[1:])
            elif head == "GROUP":
                nl.groups.append(NLGroupInfo(int(t[1]), t[2].upper(), int(t[3]), int(t[4])))
            else:
                if len(t) < 7:
                    raise ValueError("a row needs group element gamma_eff G beta curve Gmax")
                nl.rows.append(NLRow(int(t[0]), int(t[1]), float(t[2]), float(t[3]), float(t[4]), int(t[5]),
                                     float(t[6]), float(t[7]) if len(t) > 7 else 0.0))
        except (ValueError, IndexError) as exc:
            raise NLSoilError(f"{p.name} line {k}: cannot read '{s[:60]}' ({exc})") from None
    if nl.iteration < 0:
        raise NLSoilError(f"{p.name}: no ITERATION line (not a SASSI-EDU FILE74/FILE78 file)")
    keys = [r.key for r in nl.rows]
    if len(set(keys)) != len(keys):
        raise NLSoilError(f"{p.name}: an element is listed twice")
    return nl


# ======================================================================================
# Curves (D-SOL-03) and the property update (requirements 4.4 item 7)
# ======================================================================================
def load_curves(path: PathLike) -> Dict[int, DynamicProperty]:
    """FILE73 curves {ICURVE: DynamicProperty} (written by SOIL, D-FIL-10)."""
    p = Path(path)
    if not p.is_file():
        raise NLSoilError(f"{p.name} missing -- run SOIL (it writes the G/Gmax and damping curves of the DYNP "
                          "labels to FILE73)")
    curves = read_file73(p)
    if not curves:
        raise NLSoilError(f"{p.name} holds no curves")
    return curves


def curve_label(curves: Optional[Dict[int, DynamicProperty]], icurve: int) -> str:
    if not curves or int(icurve) not in curves:
        return "?"
    return curves[int(icurve)].label


def curve_problems(curve: DynamicProperty) -> List[str]:
    """What makes a FILE73 curve unusable for the property update: no G/G_max point or no damping
    point with a positive strain (log-strain interpolation, D-SOL-03), G/G_max outside [0, 1], negative
    damping (the SOIL checks Errors 97-99 for curves that SPRO does not use)."""
    out = []
    gs, gr = np.asarray(curve.g_strain, float), np.asarray(curve.g_ratio, float)
    ds, dp = np.asarray(curve.d_strain, float), np.asarray(curve.d_pct, float)
    if not np.any(gs > 0):
        out.append("no G/G_max point with a positive strain (Error 97)")
    elif np.any(gr[gs > 0] < 0) or np.any(gr[gs > 0] > 1):
        out.append("G/G_max values outside [0, 1] (Error 98)")
    if not np.any(ds > 0):
        out.append("no damping point with a positive strain (Error 99)")
    elif np.any(dp[ds > 0] < 0):
        out.append("negative damping values")
    return out


def check_curves(curves: Dict[int, DynamicProperty], numbers: Iterable[int]) -> None:
    """Raise :class:`NLSoilError` when one of the FILE73 curves ``numbers`` is missing or unusable."""
    errs = []
    for k in sorted(set(int(n) for n in numbers)):
        cv = curves.get(k)
        if cv is None:
            errs.append(f"curve {k} is not in FILE73 (curves {sorted(curves)})")
            continue
        errs += [f"FILE73 curve {k} ({cv.label}): {t}" for t in curve_problems(cv)]
    if errs:
        raise NLSoilError("; ".join(errs) + " -- correct the DYNP points and re-run SOIL")


def strain_compatible(curve: DynamicProperty, gamma_eff_pct, gmax) -> Tuple[np.ndarray, np.ndarray]:
    """``G = G_max (G/G_max)(gamma_eff)`` and ``beta = D(gamma_eff)`` (fraction) for effective strains
    in percent (requirements 4.4 item 7; linear in log10 strain, constant outside, D-SOL-03).

    Raises :class:`NLSoilError` (not a bare ValueError) for an unusable curve."""
    probs = curve_problems(curve)
    if probs:
        raise NLSoilError(f"curve {curve.label}: " + "; ".join(probs) + " -- correct the DYNP points and re-run SOIL")
    g = np.atleast_1d(np.asarray(gamma_eff_pct, dtype=float))
    try:
        ratio = np.atleast_1d(curve.g_over_gmax(g))
        beta = np.atleast_1d(curve.damping(g))
    except ValueError as exc:                       # pragma: no cover - curve_problems catches these
        raise NLSoilError(f"curve {curve.label}: {exc}") from None
    return np.asarray(gmax, dtype=float) * ratio, beta


def update_rows(rows: Sequence[NLRow], curves: Dict[int, DynamicProperty]) -> List[NLRow]:
    """Rows with G and beta recomputed from their effective strain and curve (FILE73)."""
    out = []
    for r in rows:
        cv = curves.get(int(r.curve))
        if cv is None:
            raise NLSoilError(f"curve {r.curve} of group {r.group} element {r.element} is not in FILE73 "
                              f"(curves {sorted(curves)})")
        G, b = strain_compatible(cv, r.gamma_eff, r.gmax)
        out.append(NLRow(r.group, r.element, r.gamma_eff, float(G[0]), float(b[0]), r.curve, r.gmax, r.extra))
    return out


# ======================================================================================
# Strain measures (requirements 4.10 item 5, D-NLS-03)
# ======================================================================================
def shear_strain_history(strains: np.ndarray, istr: int, etype: str = "SOLID") -> np.ndarray:
    """Shear strain measure gamma(t) from strain-component histories (engineering shear strains).

    ``strains`` (..., nc): SOLID [EXX EYY EZZ GXY GXZ GYZ] (nc = 6), PLANE [EXX EZZ GXZ] (nc = 3), the
    rows of the FILE4 strain operators.  ISTR 0: maximum shear-strain component; ISTR 1: octahedral
    (3D) or maximum shear strain (2D) -- see the module docstring."""
    e = np.asarray(strains, dtype=float)
    et = etype.upper()
    if et == "SOLID":
        if e.shape[-1] != 6:
            raise ValueError("SOLID strains need 6 components (EXX EYY EZZ GXY GXZ GYZ)")
        if int(istr) == 0:
            return np.max(np.abs(e[..., 3:6]), axis=-1)
        exx, eyy, ezz, gxy, gxz, gyz = (e[..., k] for k in range(6))
        return (2.0 / 3.0) * np.sqrt((exx - eyy) ** 2 + (eyy - ezz) ** 2 + (ezz - exx) ** 2
                                     + 1.5 * (gxy ** 2 + gyz ** 2 + gxz ** 2))
    if et == "PLANE":
        if e.shape[-1] != 3:
            raise ValueError("PLANE strains need 3 components (EXX EZZ GXZ)")
        if int(istr) == 0:
            return np.abs(e[..., 2])
        return np.sqrt((e[..., 0] - e[..., 1]) ** 2 + e[..., 2] ** 2)
    raise ValueError(f"nonlinear soil elements must be SOLID or PLANE, not {etype}")


def effective_strain(gamma_t: np.ndarray, esf: float) -> Tuple[np.ndarray, np.ndarray]:
    """``(gamma_max, gamma_eff)`` per element (axis 0 = time): ``gamma_eff = ESF max|gamma(t)|``."""
    gmax = np.max(np.abs(np.asarray(gamma_t, dtype=float)), axis=0)
    return gmax, float(esf) * gmax


def srss(*gammas) -> np.ndarray:
    """SRSS of directional effective strains, ``sqrt(gx^2 + gy^2 + gz^2)`` (D-NLS-04, spec 03 item 18).

    UT-20: 0.10, 0.08, 0.02 % -> 0.1296 %."""
    arrs = [np.asarray(g, dtype=float) for g in gammas if g is not None]
    if not arrs:
        raise ValueError("SRSS of nothing")
    return np.sqrt(np.sum([a * a for a in arrs], axis=0))


def combine_files(files: Sequence[NLFile], curves: Optional[Dict[int, DynamicProperty]] = None,
                  labels: Sequence[str] = ()) -> Tuple[NLFile, List[str]]:
    """COMB_XYZ_STRAIN: SRSS of the effective strains of directional FILE74s, element by element
    (D-NLS-04).  The directional files must be the responses of one HOUSE run: the same elements, the
    same iteration, ESF and FILE4 hash (an error otherwise -- a stale directional file would mix two
    iterations).  G and beta are recomputed from FILE73 at the combined strain when ``curves`` is
    given (else they are left 0 and HOUSE recomputes them anyway).  Returns ``(combined FILE74,
    warnings)``."""
    if not files:
        raise NLSoilError("COMBXYZSTRAIN needs at least one FILE74")
    warn: List[str] = []
    ref = files[0]
    keys = [r.key for r in ref.rows]
    kset = set(keys)
    file4 = ref.file4
    for k, f in enumerate(files[1:], start=1):
        ks = set(r.key for r in f.rows)
        if ks != kset:
            miss = sorted(kset ^ ks)[:5]
            raise NLSoilError(f"the FILE74 files hold different elements (e.g. group/element {miss}); they must "
                              "come from STRESS runs of the same HOUSE model")
        if f.iteration != ref.iteration:
            raise NLSoilError(f"the directional FILE74 files are of different iterations ({ref.iteration} and "
                              f"{f.iteration}): run STRESS for every direction of the current iteration")
        if abs(f.esf - ref.esf) > 1e-12:
            raise NLSoilError(f"the directional FILE74 files have different ESF ({ref.esf:g}, {f.esf:g}): they are "
                              "not of the same HOUSE run")
        if not same_file4(file4, f.file4):
            raise NLSoilError("the directional FILE74 files are responses of different HOUSE runs (FILE4 hash): run "
                              "STRESS for every direction after the same HOUSE and ANALYS runs")
        file4 = file4 or f.file4
    maps = [f.by_key() for f in files]
    g_dir = [np.array([m[k].gamma_eff for k in keys]) for m in maps]
    gmx_dir = [np.array([m[k].extra for k in keys]) for m in maps]
    g = srss(*g_dir)
    gm = srss(*gmx_dir)
    dirs = [lab for lab in labels if lab] or [f.direction or f"#{i + 1}" for i, f in enumerate(files)]
    out = NLFile("FILE74", ref.iteration, ref.esf, ref.ncurv, direction="SRSS " + "+".join(dirs),
                 source=f"COMBXYZSTRAIN of {len(files)} directional FILE74", file4=file4, model=ref.model,
                 groups=list(ref.groups))
    for j, k in enumerate(keys):
        r0 = maps[0][k]
        out.rows.append(NLRow(r0.group, r0.element, float(g[j]), 0.0, 0.0, r0.curve, r0.gmax, float(gm[j])))
    if curves:
        out.rows = update_rows(out.rows, curves)
    else:
        warn.append("FILE73 not found: the combined FILE74 holds the strains only (G, beta = 0); HOUSE recomputes "
                    "the properties from the strains and FILE73")
    return out, warn


# ======================================================================================
# Convergence (D-NLS-02)
# ======================================================================================
@dataclass
class Convergence:
    """Change from the properties used in iteration ``iteration`` to the strain-compatible ones."""
    iteration: int
    n: int
    dG_pct: float              # signed max-magnitude (G_new - G_used)/G_new x 100
    dG_elem: Tuple[int, int]
    dbeta_pct: float           # signed max-magnitude (beta_new - beta_used) x 100 (percentage points)
    dbeta_elem: Tuple[int, int]
    mean_ratio: float          # mean G_new/G_max
    mean_beta: float           # mean beta_new
    tol_g: float = TOL_G_PCT
    tol_beta: float = TOL_BETA_PCT

    @property
    def converged(self) -> bool:
        return abs(self.dG_pct) < self.tol_g and abs(self.dbeta_pct) < self.tol_beta

    def text(self) -> str:
        return (f"iteration {self.iteration}: max |dG/G| = {abs(self.dG_pct):.3f} % (group {self.dG_elem[0]} element "
                f"{self.dG_elem[1]}), max |d beta| = {abs(self.dbeta_pct):.3f} % (group {self.dbeta_elem[0]} element "
                f"{self.dbeta_elem[1]}) -> " + ("CONVERGED" if self.converged else "not converged")
                + f" (D-NLS-02: < {self.tol_g:g} % and < {self.tol_beta:g} %)")


def compare(used: Sequence[NLRow], new: Sequence[NLRow], iteration: int, tol_g: float = TOL_G_PCT,
            tol_beta: float = TOL_BETA_PCT) -> Convergence:
    """Convergence measures between the properties *used* (FILE78) and the *new* strain-compatible
    properties (FILE74, or HOUSE's update) of the same elements.  ``dG`` is relative to the new
    value as in SHAKE (R1 section 6.3); ``d beta`` is absolute (percentage points)."""
    nm = {r.key: r for r in new}
    keys = [r.key for r in used if r.key in nm]
    if not keys:
        raise NLSoilError("no common elements between the properties used (FILE78) and the new ones (FILE74)")
    um = {r.key: r for r in used}
    Gu = np.array([um[k].G for k in keys])
    Gn = np.array([nm[k].G for k in keys])
    bu = np.array([um[k].beta for k in keys])
    bn = np.array([nm[k].beta for k in keys])
    gm = np.array([nm[k].gmax for k in keys])
    dG = np.where(Gn != 0, (Gn - Gu) / np.where(Gn != 0, Gn, 1.0), 0.0) * 100.0
    db = (bn - bu) * 100.0
    iG = int(np.argmax(np.abs(dG)))
    ib = int(np.argmax(np.abs(db)))
    ratio = np.where(gm > 0, Gn / np.where(gm > 0, gm, 1.0), np.nan)
    return Convergence(int(iteration), len(keys), float(dG[iG]), keys[iG], float(db[ib]), keys[ib],
                       float(np.nanmean(ratio)), float(np.mean(bn)), float(tol_g), float(tol_beta))


_CONV_COLS = ("iteration", "n_elem", "max_dG_pct", "group", "element", "max_dbeta_pct", "group", "element",
              "mean_G_Gmax", "mean_beta", "converged")


def record_convergence(workdir: PathLike, conv: Convergence) -> Path:
    """Add (or replace) the row of ``conv.iteration`` in ``NLSOIL_CONVERGENCE.TXT``."""
    p = Path(workdir) / CONVERGENCE_FILE
    rows = {int(r["iteration"]): r for r in read_convergence(workdir)}
    rows[int(conv.iteration)] = {"iteration": conv.iteration, "n_elem": conv.n, "max_dG_pct": conv.dG_pct,
                                 "g1": conv.dG_elem[0], "e1": conv.dG_elem[1], "max_dbeta_pct": conv.dbeta_pct,
                                 "g2": conv.dbeta_elem[0], "e2": conv.dbeta_elem[1], "mean_G_Gmax": conv.mean_ratio,
                                 "mean_beta": conv.mean_beta, "converged": int(conv.converged)}
    lines = ["# SASSI-EDU nonlinear soil SSI convergence history (D-NLS-02: converged when max |dG/G| < "
             f"{conv.tol_g:g} % and max |d beta| < {conv.tol_beta:g} %)",
             "# row k: properties used in iteration k (FILE78) -> strain-compatible properties of its response "
             "(FILE74); dG relative to the new G, d beta in percentage points",
             "# " + " ".join(_CONV_COLS)]
    for k in sorted(rows):
        r = rows[k]
        lines.append(f"{int(r['iteration']):4d} {int(r['n_elem']):7d} {float(r['max_dG_pct']):12.5f} {int(r['g1']):5d} "
                     f"{int(r['e1']):7d} {float(r['max_dbeta_pct']):12.5f} {int(r['g2']):5d} {int(r['e2']):7d} "
                     f"{float(r['mean_G_Gmax']):11.6f} {float(r['mean_beta']):10.6f} {int(r['converged']):3d}")
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def read_convergence(workdir: PathLike) -> List[Dict[str, float]]:
    """Rows of ``NLSOIL_CONVERGENCE.TXT`` (empty when absent)."""
    p = Path(workdir) / CONVERGENCE_FILE
    if not p.is_file():
        return []
    out = []
    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s[0] == "#":
            continue
        t = s.split()
        if len(t) < 11:
            continue
        try:
            out.append({"iteration": int(t[0]), "n_elem": int(t[1]), "max_dG_pct": float(t[2]), "g1": int(t[3]),
                        "e1": int(t[4]), "max_dbeta_pct": float(t[5]), "g2": int(t[6]), "e2": int(t[7]),
                        "mean_G_Gmax": float(t[8]), "mean_beta": float(t[9]), "converged": int(t[10])})
        except ValueError:
            continue
    return out


def reset_convergence(workdir: PathLike) -> None:
    """Start a new history (NLSSIRESET)."""
    p = Path(workdir) / CONVERGENCE_FILE
    if p.exists():
        p.unlink()


def status(workdir: PathLike, tol_g: float = TOL_G_PCT, tol_beta: float = TOL_BETA_PCT) -> Optional[Convergence]:
    """Convergence of the last complete iteration in ``workdir``: FILE78 (properties used) against
    FILE74 (strain-compatible properties of the response), when both belong to the same iteration --
    the same iteration number and, when both files carry it, the same FILE4 hash (the strains are the
    response to these properties) -- and FILE73 is present.  ``None`` when the files do not describe
    a complete iteration."""
    wd = Path(workdir)
    p74, p78, p73 = wd / "FILE74", wd / "FILE78", wd / "FILE73"
    if not (p74.is_file() and p78.is_file() and p73.is_file()):
        return None
    try:
        f78 = read_nlfile(p78, "FILE78")
        f74 = read_nlfile(p74, "FILE74")
        if f74.iteration != f78.iteration or not same_file4(f78.file4, f74.file4):
            return None
        new = update_rows(f74.rows, load_curves(p73))
        return compare(f78.rows, new, f78.iteration, tol_g, tol_beta)
    except NLSoilError:
        return None


def current_status(workdir: PathLike, model: str, tol_g: float = TOL_G_PCT,
                   tol_beta: float = TOL_BETA_PCT) -> Tuple[Optional[Convergence], str]:
    """:func:`status` of the analysis *in progress*, or ``(None, reason)``.

    The last complete iteration (FILE78 + FILE74) is that of the analysis in progress only when

    * ``<model>.liq`` = 1 (NLSSIRESET sets 0: a new analysis has been requested);
    * FILE78 was written with the FILE4 now in the directory (its FILE4 hash = FILE90 ``file4_hash``);
    * the HOUSE deck and the .pin are those HOUSE read for FILE78 (deck and .pin digests: the model,
      the near-field materials, ESF, ISTR or the curves have not changed since).

    NLSSIITER reports "already converged" only then (requirements 2.4, D-NLS-02)."""
    wd = Path(workdir)
    if read_liq(wd / f"{model}.liq") != 1:
        return None, (f"{model}.liq is not 1 (NLSSIRESET requested a new analysis, or no non-linear HOUSE run "
                      "has been made)")
    st = status(wd, tol_g, tol_beta)
    if st is None:
        return None, "FILE78 and FILE74 do not describe one complete iteration"
    try:
        f78 = read_nlfile(wd / "FILE78", "FILE78")
    except NLSoilError as exc:                       # pragma: no cover - status() read it already
        return None, str(exc)
    try:
        from ..io.container import read_container
        h90 = str(read_container(wd / "FILE90").meta.get("file4_hash", ""))
    except Exception:                                # noqa: BLE001 -- missing or unreadable FILE90
        h90 = ""
    if not f78.file4 or f78.file4 != h90:
        return None, "FILE78 does not belong to the FILE4 in the model directory (FILE90 hash)"
    from ..io.decks import deck_path
    hou = deck_path(wd, model, "HOUSE")
    if not f78.deck or deck_digest(hou) != f78.deck:
        return None, f"{hou.name} has changed since the HOUSE run of FILE78 (iteration {f78.iteration})"
    try:
        pin = read_pin(wd / f"{model}.pin")
    except (NLSoilError, OSError):
        return None, f"{model}.pin is missing or invalid"
    if not f78.pin or pin_digest(pin) != f78.pin:
        return None, f"{model}.pin has changed since the HOUSE run of FILE78 (iteration {f78.iteration})"
    return st, ""


# ======================================================================================
# Element materials (requirements 4.4 item 7; constant Poisson's ratio, D-SOL-06)
# ======================================================================================
def nonlinear_material(base: MaterialProps, G0: float, beta_s: float, beta_p: Optional[float] = None) -> MaterialProps:
    """Internal material of one nonlinear element: shear modulus ``G0`` (real part), the constrained
    modulus scaled with it (``M0 = M0_base G0 / G0_base``: Poisson's ratio and the density are those of
    the base material), damping ``beta_s`` on G and ``beta_p`` (default ``beta_s``) on M (D-SOL-06)."""
    bp = float(beta_s if beta_p is None else beta_p)
    bs = float(beta_s)
    for b, what in ((bs, "shear"), (bp, "P-wave")):
        if not 0.0 <= b < 0.5:
            raise NLSoilError(f"{what} damping ratio {b:g} outside [0, 0.5) (EDU-04): check the damping curve")
    if base.G0 <= 0 or G0 <= 0:
        raise NLSoilError("nonlinear soil elements need a positive shear modulus")
    M0 = base.M0 * (float(G0) / base.G0)
    return MaterialProps(rho=base.rho, G=complex(G0 * cfactor(bs, base.form)), M=complex(M0 * cfactor(bp, base.form)),
                         beta_s=bs, beta_p=bp, G0=float(G0), M0=float(M0), form=base.form)


def initial_material(base: MaterialProps, ref: MaterialProps, gfac: float, dfac: float) -> MaterialProps:
    """Iteration-0 material of a nonlinear element (requirements 4.4 item 7): ``G = GFAC G_ref``,
    ``beta_s = DFAC beta_s,ref``, ``beta_p = DFAC beta_p,ref`` with ``ref`` the free-field layer at the
    element (``G_layer``, ``beta_layer``); density and Poisson's ratio of the element material ``base``
    (its constrained modulus follows G, as in the later iterations).  Returns ``base`` itself when the
    result is identical to it (the zero-SSI identity stays bit-exact, VP-16)."""
    G0 = float(gfac) * ref.G0
    bs, bp = float(dfac) * ref.beta_s, float(dfac) * ref.beta_p
    if G0 == base.G0 and bs == base.beta_s and bp == base.beta_p:
        return base
    return nonlinear_material(base, G0, bs, bp)


# ======================================================================================
# HOUSE: the properties of one iteration (requirements 4.4 item 7)
# ======================================================================================
@dataclass
class NLElementIn:
    """A HOUSE element as seen by the nonlinear plan (built by sassi.modules.house).

    ``layer`` is the free-field material at the element (the TOPL layer containing its mid-depth,
    L properties, built like the excavated soil, D-ELM-12), ``layer_no`` its L number and ``depth``
    the mid-depth; ``layer`` None = no free-field layer table (the element material is used)."""
    index: int                 # position in the HOUSE record list
    group: int
    element: int
    code: int                  # 1 SOLID, 4 PLANE ...
    excavated: bool
    mat: int                   # material number (M table; L layer when excavated)
    base: Optional[MaterialProps]
    layer: Optional[MaterialProps] = None
    layer_no: int = 0
    depth: float = float("nan")


@dataclass
class NLElement:
    """One nonlinear element and the properties HOUSE assigns to it (``base`` = the element's M-table
    material, the low-strain soil: G_max = base.G0)."""
    index: int
    group: int
    element: int
    mat: int
    base: MaterialProps
    material: MaterialProps
    gamma_eff: float           # % (0 in the initial iteration)
    curve: int
    istr: int
    gfac: float
    dfac: float
    ref_layer: int = 0         # iteration 0: L number of the free-field layer scaled by GFAC/DFAC (0 = material)
    depth: float = float("nan")

    @property
    def ratio(self) -> float:
        return self.material.G0 / self.base.G0


@dataclass
class NLPlan:
    iteration: int
    pin: PinData
    source: str                            # 'PIN' or 'FILE74 ...'
    elements: List[NLElement] = field(default_factory=list)
    groups: List[NLGroupInfo] = field(default_factory=list)
    curves: Optional[Dict[int, DynamicProperty]] = None
    convergence: Optional[Convergence] = None
    previous: Optional[NLFile] = None      # FILE78 of the previous iteration (when it matches FILE74)
    file74: Optional[NLFile] = None
    warnings: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def file78(self, model: str = "", file4: str = "", deck: str = "") -> NLFile:
        """FILE78 of this plan (``file4`` = FILE90 hash of the FILE4 written with these properties,
        ``deck`` = :func:`deck_digest` of the HOUSE deck; the .pin digest is added here)."""
        nl = NLFile("FILE78", self.iteration, self.pin.esf, self.pin.ncurv, source=self.source, file4=file4,
                    model=model, groups=list(self.groups), deck=deck, pin=pin_digest(self.pin))
        for e in self.elements:
            nl.rows.append(NLRow(e.group, e.element, e.gamma_eff, e.material.G0, e.material.beta_s, e.curve,
                                 e.base.G0, float(e.mat)))
        return nl

    def rows(self) -> List[NLRow]:
        return self.file78().rows


def plan_iteration(pin: PinData, elements: Sequence[NLElementIn], dim: int, liq: int,
                   file74: Optional[NLFile] = None, curves: Optional[Dict[int, DynamicProperty]] = None,
                   previous78: Optional[NLFile] = None) -> NLPlan:
    """Properties of the nonlinear elements for one HOUSE run (requirements 4.4 item 7).

    * ``liq`` != 1 (initial run, iteration 0): ``G = GFAC G_layer``, ``beta = DFAC beta_layer`` of the
      free-field layer at the element (:attr:`NLElementIn.layer`; :func:`initial_material`); an element
      without a layer falls back to its own material (a warning);
    * ``liq`` == 1: ``G = G_max (G/G_max)(gamma_eff)``, ``beta_s = beta_p = D(gamma_eff)`` with
      ``G_max`` the low-strain modulus of the element's M-table material, gamma_eff from ``file74`` and
      the curve ICURVE of ``curves`` (FILE73); iteration = FILE74 iteration + 1.  When ``previous78``
      (the FILE78 of the run that produced the FILE74 strains: same iteration and FILE4 hash) is
      given, the change of the properties is the convergence measure of that iteration.

    ``elements`` are all HOUSE elements (only those of the .pin groups are used).  The NMAT
    material lines of a group map to its distinct material numbers in ascending order (spec 05b
    OQ 22).  The curves referenced by the .pin must be usable (:func:`curve_problems`) whenever FILE73
    is given.  Raises :class:`NLSoilError` with every problem found."""
    errs: List[str] = []
    plan = NLPlan(iteration=0, pin=pin, source="PIN")
    want = 4 if int(dim) == 1 else 1
    by_group: Dict[int, List[NLElementIn]] = {}
    for e in elements:
        by_group.setdefault(int(e.group), []).append(e)
    chosen: List[Tuple[NLElementIn, PinGroup, PinMaterial]] = []
    for g in pin.groups:
        els = by_group.get(g.igrp, [])
        if not els:
            errs.append(f"nonlinear soil group {g.igrp} of the .pin has no elements in the HOUSE model")
            continue
        codes = sorted({e.code for e in els})
        if codes != [want]:
            errs.append(f"nonlinear soil group {g.igrp} must hold {NL_TYPES[want]} elements in a "
                        f"{'2D' if want == 4 else '3D'} model (found type codes {codes})")
            continue
        exc = [e.element for e in els if e.excavated]
        if exc:
            errs.append(f"nonlinear soil group {g.igrp} holds excavated soil elements ({_ids(exc)}): the near-field "
                        "soil is part of the structure model (ETYPE 1, an M-table material with the soil properties), "
                        "placed in the excavated volume next to the excavated free-field soil (its own ETYPE 2 group)")
            continue
        mats = sorted({e.mat for e in els})
        if len(mats) != g.nmat:
            errs.append(f"nonlinear soil group {g.igrp}: the .pin gives NMAT = {g.nmat} material lines but the "
                        f"group uses {len(mats)} materials {mats} (lines map to them in ascending order)")
            continue
        if g.nelem != len(els):
            plan.warnings.append(f"nonlinear soil group {g.igrp}: NELEM = {g.nelem} in the .pin, the group has "
                                 f"{len(els)} elements (all of them are nonlinear)")
        if g.istr == 1 and want == 4:
            plan.notes.append(f"group {g.igrp}: ISTR 1 = maximum shear strain of the 2D (X-Z) strain state")
        mline = {m: g.materials[k] for k, m in enumerate(mats)}
        for e in sorted(els, key=lambda x: x.element):
            if e.base is None:
                errs.append(f"element {e.element} of group {g.igrp} has no material")
                continue
            chosen.append((e, g, mline[e.mat]))
        plan.groups.append(NLGroupInfo(g.igrp, NL_TYPES[want], len(els), g.istr))
    if curves is not None:
        for g in pin.groups:
            for m in g.materials:
                if m.icurve not in curves:
                    errs.append(f"group {g.igrp}: ICURVE = {m.icurve} is not a curve of FILE73 (curves {sorted(curves)})")
                else:
                    errs += [f"group {g.igrp}: FILE73 curve {m.icurve} ({curves[m.icurve].label}) is not usable: {t}"
                             for t in curve_problems(curves[m.icurve])]
        if len(curves) != pin.ncurv:
            plan.warnings.append(f"the .pin says NCURV = {pin.ncurv}, FILE73 holds {len(curves)} curves")
    if errs:
        raise NLSoilError("; ".join(dict.fromkeys(errs)))
    plan.curves = curves
    if int(liq) != 1:
        no_layer, above, stiff = [], [], []
        for e, g, m in chosen:
            ref = e.layer if e.layer is not None else e.base
            if e.layer is None:
                no_layer.append(e.element)
            elif e.depth < 0:
                above.append(e.element)
            mat = initial_material(e.base, ref, m.gfac, m.dfac)
            if mat.G0 > e.base.G0 * (1.0 + 1e-6):
                stiff.append(e.element)
            plan.elements.append(NLElement(e.index, e.group, e.element, e.mat, e.base, mat, 0.0, m.icurve, g.istr,
                                           m.gfac, m.dfac, ref_layer=int(e.layer_no) if e.layer is not None else 0,
                                           depth=float(e.depth)))
        if no_layer:
            plan.warnings.append(f"no free-field layer table: iteration 0 scales the element material instead of the "
                                 f"free-field layer (G = GFAC G_mat) for {len(no_layer)} elements ({_ids(no_layer)})")
        if above:
            plan.notes.append(f"{len(above)} nonlinear elements lie above the ground surface ({_ids(above)}): their "
                              "iteration-0 properties scale TOPL layer 1")
        if stiff:
            plan.warnings.append(f"iteration 0: {len(stiff)} nonlinear elements ({_ids(stiff)}) start stiffer than the "
                                 "low-strain modulus G_max of their material (GFAC G_layer > G_mat): GFAC scales the "
                                 "free-field layer (requirements 4.4 item 7) and the M-table material must hold the "
                                 "low-strain soil -- check GFAC and the material")
        return plan
    # ---- later iterations: strain-compatible properties from FILE74 and FILE73
    if file74 is None:
        raise NLSoilError("FILE74 missing -- run STRESS with <iter> = 1 (Auto Computation of Strains in Soil El.) "
                          "and, for X/Y/Z input, COMBXYZSTRAIN; or start a new analysis (NLSSIRESET / delete the .liq)")
    if curves is None:
        raise NLSoilError("FILE73 (soil curves) is needed to update the near-field properties -- run SOIL")
    plan.file74 = file74
    plan.iteration = int(file74.iteration) + 1
    plan.source = f"FILE74 iteration {file74.iteration}" + (f" ({file74.direction})" if file74.direction else "")
    if abs(file74.esf - pin.esf) > 1e-12:
        plan.warnings.append(f"FILE74 was computed with ESF = {file74.esf:g}, the .pin now says {pin.esf:g}: the "
                             "strains are used as they are (the new ESF applies from the next STRESS run)")
    rows = file74.by_key()
    missing = [(e.group, e.element) for e, _, _ in chosen if (e.group, e.element) not in rows]
    if missing:
        raise NLSoilError(f"FILE74 has no strain for {len(missing)} nonlinear elements (e.g. group/element "
                          f"{missing[:5]}): it belongs to another model or nonlinear group definition -- start a "
                          "new analysis (NLSSIRESET, or delete the .liq file)")
    for e, g, m in chosen:
        r = rows[(e.group, e.element)]
        if r.curve != m.icurve:
            plan.warnings.append(f"group {e.group} element {e.element}: FILE74 refers to curve {r.curve}, the .pin to "
                                 f"{m.icurve}; the .pin curve is used")
        if r.gmax > 0 and abs(r.gmax - e.base.G0) > 1e-9 * e.base.G0:
            plan.warnings.append(f"group {e.group} element {e.element}: G_max {r.gmax:.6g} in FILE74 differs from "
                                 f"the material value {e.base.G0:.6g} (material changed?)")
        G, b = strain_compatible(curves[m.icurve], r.gamma_eff, e.base.G0)
        mat = nonlinear_material(e.base, float(G[0]), float(b[0]))
        plan.elements.append(NLElement(e.index, e.group, e.element, e.mat, e.base, mat, float(r.gamma_eff),
                                       m.icurve, g.istr, m.gfac, m.dfac, depth=float(e.depth)))
    if previous78 is not None and previous78.iteration == file74.iteration and same_file4(previous78.file4,
                                                                                            file74.file4):
        plan.previous = previous78
        plan.convergence = compare(previous78.rows, plan.rows(), file74.iteration)
    elif previous78 is not None and previous78.iteration == file74.iteration:
        plan.notes.append(f"the FILE74 strains of iteration {file74.iteration} are not the response to the FILE78 "
                          "properties in the directory (FILE4 hash): no convergence measure for this update")
    elif previous78 is not None:
        plan.notes.append(f"FILE78 in the directory is of iteration {previous78.iteration}, the FILE74 strains of "
                          f"iteration {file74.iteration}: no convergence measure for this update")
    return plan


def _ids(ids: Iterable[int], n: int = 8) -> str:
    ids = sorted(set(int(i) for i in ids))
    return ", ".join(str(i) for i in ids[:n]) + (" ..." if len(ids) > n else "")


# ======================================================================================
# STRESS: FILE74 of one directional run (requirements 4.10 item 5)
# ======================================================================================
def strain_file(file78: NLFile, gamma_max_pct: Dict[Tuple[int, int], float], curves: Dict[int, DynamicProperty],
                direction: str = "") -> NLFile:
    """FILE74 of a STRESS run: ``gamma_eff = ESF gamma_max`` (percent) and the strain-compatible G and
    beta from FILE73 for every element of FILE78 (``gamma_max_pct`` holds max|gamma(t)| in percent).
    The FILE4 hash of FILE78 is carried over: these strains are the response to its properties."""
    out = NLFile("FILE74", file78.iteration, file78.esf, file78.ncurv, direction=direction,
                 source="STRESS", file4=file78.file4, model=file78.model, groups=list(file78.groups))
    rows = []
    for r in file78.rows:
        gm = float(gamma_max_pct[r.key])
        rows.append(NLRow(r.group, r.element, file78.esf * gm, 0.0, 0.0, r.curve, r.gmax, gm))
    out.rows = update_rows(rows, curves)
    return out


def group_summary(rows: Sequence[NLRow]) -> List[Tuple[int, int, float, float, float, float, float, float]]:
    """Per group: (group, n, min/mean/max G/G_max, min/mean/max beta) for the listings."""
    out = []
    for g in sorted({r.group for r in rows}):
        rr = [r for r in rows if r.group == g]
        ratio = np.array([r.ratio for r in rr])
        beta = np.array([r.beta for r in rr])
        out.append((g, len(rr), float(np.min(ratio)), float(np.mean(ratio)), float(np.max(ratio)),
                    float(np.min(beta)), float(np.mean(beta)), float(np.max(beta))))
    return out


__all__ = [
    "TOL_G_PCT", "TOL_BETA_PCT", "MAX_ITERATIONS", "STRAIN_FLAGS", "STRAIN_COMPONENTS", "DIRECTIONS",
    "CONVERGENCE_FILE", "NLSoilError", "PinMaterial", "PinGroup", "PinData", "parse_pin", "read_pin", "format_pin",
    "write_pin", "pin_digest", "deck_digest", "read_liq", "write_liq", "RETIRED_SUFFIX", "retire", "same_file4",
    "NLRow", "NLGroupInfo", "NLFile", "format_nlfile", "write_nlfile", "read_nlfile", "load_curves", "curve_label",
    "curve_problems", "check_curves", "strain_compatible", "update_rows", "shear_strain_history",
    "effective_strain", "srss", "combine_files", "Convergence", "compare", "record_convergence", "read_convergence",
    "reset_convergence", "status", "current_status", "nonlinear_material", "initial_material", "NLElementIn",
    "NLElement", "NLPlan", "plan_iteration", "strain_file", "group_summary",
]
