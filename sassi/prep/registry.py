"""Command registry and the catalogue of every ``.pre`` command (requirements sections 3.2-3.4).

Handlers are plain functions registered with the :func:`command` decorator in a module of
:mod:`sassi.prep.commands` (every module of that package is imported automatically)::

    from sassi.prep.registry import command, CommandError

    @command("NGEN")                       # abbreviation, tier, class come from the catalogue
    def ngen(c):                            # c is a sassi.prep.interpreter.Call
        itim = c.int(1, default=1)
        ...

The **catalogue** below lists every command of the manual and of the SASSI-EDU dialect with its
documented abbreviation (section 3.2 and D-PAR-04), priority tier (section 0.1), state class
(section 3.3) and lexical options (rest-of-line text position L7, raw FOREACH body L13, legacy
SOIL parenthesis L9).  After the handler modules are imported, every catalogued command without
a handler gets a *placeholder*:

* storage classes (record, indexed, string setters) store the command generically in
  ``model.options`` so WRITE writes it back (lower tiers must be "parsed and stored", section 0.1);
  P1/P2 ones also print ``<CMD> is not available in this build (tier Pn)``;
* action commands print ``<CMD> is not available in this build (tier Pn)``.

A real handler always replaces a placeholder.  Registering a second real handler for the same
name is an error unless ``replaces=True`` is given.

Command lookup (rule L10, D-PAR-04): exact full name first, then the exact documented
abbreviation or alias; no prefix matching.
"""
from __future__ import annotations

import importlib
import pkgutil
import traceback
from dataclasses import dataclass, field, replace
from typing import Callable, Dict, List, Optional, Tuple


class CommandError(Exception):
    """Raised by a handler: the message is reported as an error and the command is skipped."""


class CommandReported(CommandError):
    """Raised by a handler that has **already reported** its error (e.g. INP of an unreadable file).

    The interpreter marks the command as failed (``execute()`` returns False, ``last_ok`` is
    False, the line is not added to the replay history) without printing a second message.
    """


Handler = Callable[["object"], None]

CLASSES = ("action", "record", "record_keep", "indexed", "list", "request", "string", "ui")
TIERS = ("P0", "P1", "P2", "OOS")


@dataclass
class CommandSpec:
    """Everything the interpreter needs to dispatch one command."""
    name: str
    handler: Optional[Handler] = None
    abbrev: Tuple[str, ...] = ()        # accepted abbreviations and aliases (exact match)
    tier: str = "P0"
    cls: str = "action"
    text_from: Optional[int] = None     # 1-based position of a rest-of-line text argument (L7)
    raw: bool = False                   # handler gets the unsubstituted rest of line (FOREACH, L13)
    paren: bool = False                 # final '(...)' token kept whole (legacy SOIL, L9)
    key: str = "arg1"                   # key rule of indexed setters
    max_args: Optional[int] = None      # arguments after max_args are dropped with a warning
    silent: bool = False                # placeholder stores without the tier message (P0 parsing)
    note: str = ""                      # message printed by the placeholder
    summary: str = ""
    placeholder: bool = False
    module: str = ""

    @property
    def stores(self) -> bool:
        return self.cls in ("record", "record_keep", "indexed", "string")


# --------------------------------------------------------------------------------------
# Catalogue: NAME  ABBR(,ALIAS...)  TIER  CLASS  [options]
#   options: text=k  raw  paren  key=label_no|arg12|word|symm  silent
# --------------------------------------------------------------------------------------
_CATALOGUE_TEXT = r"""
# ---- 3.4.A session, files, models, global options
INP          -     P0 action
WRITE        WRIT  P0 action
SAVE         -     P0 action
RESUME       RESU  P0 action
MDL          -     P0 action
MDLNAME      -     P0 action
TIT          -     P0 string   text=1
STATUS       STAT  P0 action
ACTM         -     P0 action
CPMODEL      -     P0 action
DMODEL       -     P0 action
MODELLIST    -     P0 action
CD           -     P0 action
MKDIR        -     P0 action
MOPT         -     P0 record
GRAVITY      -     P0 action
GROUNDELEV   -     P0 action
AOPT         -     P0 record
CHECK        CHEC  P0 action
AFWRITE      AFWR  P0 action
AFWRBAT      -     P1 action
SETENV       -     P2 action
GETENV       -     P2 action
# ---- 3.4.B frequency sets
FREQ         -     P0 list
LFREQ        LFRE  P0 action
# ---- 3.4.C nodes and coordinate systems
N            -     P0 action
NDEL         -     P0 action
NGEN         -     P0 action
FILL         -     P0 action
NMED         -     P0 action
LMOVE        LMOV  P0 action
NMOVE        NMOV  P0 action
NSCALE       NSCA  P0 action
NLIST        NLIS  P0 action
D            -     P0 action
INT          -     P0 action
INTLIST      INTL  P0 action
CSYS         -     P0 action
LOC          -     P0 action
LOCAL        LOCA  P0 action
GLOBAL       GLOB  P0 action
SDEL         -     P0 action
SLIST        SLIS  P0 action
# ---- 3.4.D groups, elements, properties
GROUP        GROU  P0 action
MTYPE        MTYP  P0 action
GTIT         -     P0 action   text=2
GDEL         -     P0 action
GLIST        GLIS  P0 action
E            -     P0 action
EGEN         -     P0 action
EDEL         -     P0 action
ECOMPR       ECOM  P0 action
ELIST        ELIS  P0 action
ETYPE        ETYP  P0 action
EINT         -     P0 action
THICK        THIC  P0 action
KI           -     P0 action
KJ           -     P0 action
MSET         -     P0 action
MACT         -     P0 action
RSET         -     P0 action
RACT         -     P0 action
M            -     P0 action
L            -     P0 action
R            -     P0 action
SC           -     P0 action
MXR          -     P0 action
MXI          -     P0 action
MXM          -     P0 action
DELM         -     P0 action
DELL         -     P0 action
DELR         -     P0 action
DELSC        DELS  P0 action
MXDEL        MXDE  P0 action
MLIST        MLIS  P0 action
LLIST        LLIS  P0 action
RLIST        RLIS  P0 action
SCLIST       SCLI  P0 action
MXLIST       MXLI  P0 action
# ---- 3.4.E loads and masses
F            -     P0 action
MM           -     P0 action
FDEL         -     P0 action
MMDEL        MMDE  P0 action
FLIST        FLIS  P0 action
MMLIST       MMLI  P0 action
FSCALE       FSCA  P0 action
MSCALE       MSCA  P0 action
MT           -     P0 action
MR           -     P0 action
MTGEN        MTGE  P0 action
MRGEN        MRGE  P0 action
MTDEL        MTDE  P0 action
MRDEL        MRDE  P0 action
MTSCALE      MTSC  P0 action
MRSCALE      MRSC  P0 action
MTLIST       MTLI  P0 action
MUNITS       MUNI  P0 action
FREAD        -     P2 action
MREAD        -     P2 action
# ---- 3.4.F module analysis options
EQUAKE       EQUA  P0 record
RSIN         -     P0 indexed  text=2
RSOUT        RSOU  P0 indexed  text=2
ACCIN        ACCI  P0 indexed  text=2
ACCOUT       ACCO  P0 indexed  text=2
TPSD         -     P0 indexed  text=2
CORR         -     P1 indexed
EQTIT        EQTI  P0 string   text=1
SOIL         -     P0 record   paren
SPRO         -     P0 indexed  text=3
DYNP         -     P0 indexed  text=6 key=label_no
SACC         -     P0 indexed
SRS          -     P0 indexed
SSTR         -     P0 indexed
SSAF         -     P0 indexed  text=7
SFOU         -     P0 indexed
DAMP         -     P0 list
THFILE       THFI  P0 string   text=1
THTIT        THTI  P0 string   text=1
SITE         -     P0 record
TOPL         -     P0 list
WAVE         -     P0 indexed
POINT        POIN  P0 record
HOUSE        HOUS  P0 record
INCOH        INCO  P1 record
WPASS        WPAS  P1 record
ME           -     P1 indexed
AMP          -     P1 list
SYMM         -     P1 indexed  key=symm silent
FORCE        FORC  P0 record
ANALYS       ANAL  P0 record
MOTION       MOTI  P0 record
NOUT         -     P0 request
STRESS       STRE  P0 record
EOUT         -     P0 request
RELD         -     P0 record
RELFILE      -     P0 string   text=1
RDND         -     P0 request
# ---- 3.4.G module runs
RUNEQUAKE    -     P0 action
RUNSOIL      -     P0 action
RUNSITE      -     P0 action
RUNPOINT     -     P0 action
RUNHOUSE     -     P0 action
RUNFORCE     -     P0 action
RUNANALYS    -     P0 action
RUNCOMBIN    -     P0 action
RUNMOTION    -     P0 action
RUNRELDISP   -     P0 action
RUNSTRESS    -     P0 action
RUNNONLINEAR -     P2 action
# ---- 3.4.H model checking
EXCSTRCHK    -     P0 action
FIXEDINT     -     P0 action
FREESPRING   -     P0 action
HINGED       -     P0 action
INTCOUNT     -     P0 action
KINT         -     P0 action
USED         -     P0 action
# ---- 3.4.I conditioning and generation
FIXSLDROT    -        P0 action
FIXSHLROT    -        P0 action
FIXSPRROT    FIXSPROT P0 action
FIXROT       -        P0 action
ETYPEGEN     -        P1 action
INTGEN       -        P1 action
RADIUS       -        P1 action
GLB2LOC      -        P1 action
GROUPMAT     -        P2 action
RMVUNUSED    -        P1 action
NCOM         -        P1 action
GCOM         -        P1 action
WELD         -        P1 action
MERGE        -        P1 action
MERGESOIL    -        P1 action
MERGEGROUP   -        P1 action
ROTATE       -        P1 action
TRANSLATE    -        P1 action
EXCAV        -        P1 action
SOILMESH     -        P1 action
CRITFREQ     -        P1 action
FRAMECOMBIN  -        P1 action
FRAMESEL     -        P1 action
MODFRAMES    -        P1 action
# ---- 3.4.J cuts, submodels, section calculations
CUTADD         -   P1 action
CUTRMV         -   P1 action
CUTVOL         -   P1 action
SLICE          -   P1 action
CUTCLR         -   P1 action
CUT2SUB        -   P1 action
CSECT          -   P1 action
EXTRACTEXCAV   -   P1 action
SPLITGROUP     -   P1 action
TRANELEM       -   P1 action
TRANVOL        -   P1 action
READSTR        -   P1 action
SECDATAOPT     -   P1 record
CALCPAR        -   P1 action
CALCMOI        -   P1 action
CALCC          -   P1 action
CALCM          -   P1 action
CALCSECTHIST   -   P1 action
CALCSECTHISTDB -   P2 action
SHEAR          -   P2 action
# ---- 3.4.K file conversion
CONVERT        -   P1 action
ANSYS          -   P1 action
ANSYSREFORMAT  -   P1 action
ANSYSMODELTYPE -   P2 record
GENMATRIXDAMP  -   P2 action
# ---- 3.4.L plotting and line mathematics
READSPEC       -   P1 action
READTH         -   P1 action
WRITESPEC      -   P1 action
WRITETH        -   P1 action
ADDITION       -   P1 action
SUBTRACTION    -   P1 action
LINECOMBIN     -   P1 action
AVERAGE        -   P1 action
SRSS           -   P1 action
BROADEN        -   P1 action
LBINCORS       -   P2 action
SPECPLOT       -   P1 ui
THPLOT         -   P1 ui
LAYERPLOT      -   P1 ui
SOILPROPPLOT   -   P1 ui
MODELPLOT      -   P1 ui
NODEPLOT       -   P1 ui
CUTPLOT        -   P1 ui
PROCFRAME      -   P1 ui
BUBBLEPLOT     -   P1 ui
CONTOURPLOT    -   P1 ui
VECTORPLOT     -   P1 ui
DEFORMPLOT     -   P1 ui
AXES           -   P1 ui
PLOTRANGE      -   P1 ui
PLOTTITLE      -   P1 ui  text=1
XTITLE         -   P1 ui  text=1
YTITLE         -   P1 ui  text=1
YTITLE2        -   P1 ui  text=1
LINENAME       -   P1 ui  text=2
MARKERS        -   P1 ui
CAPTUREPLOT    -   P1 ui
CLOSEPLOT      -   P1 ui
WINDOWSETTINGS -   P1 ui
CNGVIEW        -   P1 ui
RSTVIEW        -   P1 ui
CNGCENTER      -   P1 ui
RSTCENTER      -   P1 ui
ELECOLOR       -   P1 ui
ELENUM         -   P1 ui
GROUPNUM       -   P1 ui
NODENUM        -   P1 ui
NODESEL        -   P1 ui
ELEMSEL        -   P1 ui
SELCLR         -   P1 ui
SHOWDOF        -   P1 ui
SHOWMASS       -   P1 ui
SHRINK         -   P1 ui
WIREFRAME      -   P1 ui
PAUSE          -   P1 ui
COLOR          -   P2 ui
SHADEROPTIONS  -   P2 ui
STIPPLE        -   P2 ui
DEBUG          -   P2 ui
# ---- 3.4.M programming
VAR          -   P1 action
SETVAR       -   P1 action
SHOWVAR      -   P1 action
VARLIST      -   P1 action
FOREACH      -   P1 action  raw
LOADMACRO    -   P1 action
MACRO        -   P1 action
MACROLIST    -   P1 action
LOADVAR      -   P1 action
RND          -   P1 action
ADDRND       -   P1 action
RNDSEED      -   P1 action
REDUCESET    -   P1 action
# ---- 3.4.N water modelling
FILLPOOL       -   P2 action
REFINEMODEL    -   P2 action
LISTPOOLINTER  -   P2 action
MERGEPOOL      -   P2 action
POOLDATA       -   P2 record
# ---- 3.4.O option NON and nonlinear soil
EQL            -   P2 record
P              -   P2 indexed
S              -   P2 indexed
B              -   P2 indexed
BBC            -   P2 indexed
BBCI           -   P2 indexed
BBCP           -   P2 indexed  key=arg12
BBCX           -   P2 indexed
BBCY           -   P2 indexed
BBCGEN         -   P2 action
DELBBC         -   P2 action
DELBM          -   P2 action
DELNLS         -   P2 action
DELSPR         -   P2 action
PDEL           -   P2 action
PLIST          -   P2 action
PNLGEN         PANELGEN P2 action
PANELIZE       -   P2 action
WALLFLR        -   P2 action
EDGE           -   P2 action
EDGEMODEL      -   P2 action
UNIPNL         -   P2 action
DGRDFLR        -   P2 action
MERGEPANEL     -   P2 action
NONLINMOTDISP  NONLINMODISP P2 action
NONLINBAT      -   P2 action
SOLIDPILE      -   P2 action
BEAMPILE       -   P2 action
DCOUPLEBEAM    -   P2 action
SOILREDEF      -   P2 action
NLSOIL         -   P2 record
NLSLAYER       -   P2 indexed
# ---- 3.4.P binary databases
BINOUT         -   P0 record_keep
LOADACCDB      LOADACCDBANI P2 action
LOADDISPDB     -   P2 action
LOADTHSDB      -   P2 action
DELDB          -   P2 action
COMBACCDB      -   P2 action
COMBDISPDB     COMDISPDB P2 action
COMBDISPDIR    -   P2 action
COMBTHSDB      -   P2 action
BINFRAMEOUT    -   P2 action
BINSTRTBL      -   P2 action
MAXDBFRAME     -   P2 action
ACCDBANI       ACCANIDB P2 action
DISPDBANI      -   P2 action
THSDBANI       -   P2 action
# ---- 3.4.Q thick shell
THSHLSTR       -   P2 record
THSHLSMH       -   P2 action
# ---- 3.4.R extension commands of SASSI-EDU
NONLINITER     -   P2 action
NONLINSAVE     -   P2 action
NONLINRESET    -   P2 action
NONLINTHD      -   P2 action
SITEX          -   P0 record  silent
SOILX          -   P0 record  silent
HOUSEX         -   P0 record  silent
ANALYSX        -   P0 record  silent
MOTIONX        -   P0 record  silent
STRESSX        -   P0 record  silent
RELDX          -   P1 record
CMODFORM       -   P0 record
EDUOPT         -   P0 indexed key=word
FCOPY          -   P0 action
FMOVE          -   P0 action
REMOVEFREQ     -   P1 action
HARMFRAME      -   P1 action
COMBXYZSTRAIN  -   P1 action
BUILDFILE77    -   P1 action
COMBXYZTHD     -   P2 action
VERIFY         -   P0 action
VERIFYREPORT   -   P0 action
LIBRARY        -   P0 action
ACTIVATEPLOT   -   P1 ui
SHOWSOIL       -   P1 ui
PIN            -   P1 action
PINDEL         -   P1 action
PINGRP         -   P1 action
PINLIST        -   P1 action
PINMAT         -   P1 action
NLSSIITER      -   P1 action
NLSSIRESET     -   P1 action
LOADGEN        -   P2 action
LOADGENDYN     -   P2 action
LGFILE         -   P2 action
LGNODE         -   P2 action
LGTIME         -   P2 action
LGMAP          -   P2 action
LGOPT          -   P2 action
LGLIST         -   P2 action
RUNLOADGEN     -   P2 action
"""

#: messages printed by placeholders (manual statements and decisions)
NOTES = {
    "SFOU": "SFOU is not usable in this version (stored only; compute Fourier spectra with EQUAKE)",
    "GENMATRIXDAMP": "GENMATRIXDAMP is not usable in this version",
    "B": "B: nonlinear beams are not usable in this version (stored only)",
    "BEAMPILE": "BEAMPILE is not applicable in this version",
    "DCOUPLEBEAM": "DCOUPLEBEAM is not usable in this version",
    "LBINCORS": "LBINCORS: no algorithm is specified by the manual",
    "FREAD": "FREAD is not available (syntax not documented, D-FRC-04)",
    "MREAD": "MREAD is not available (syntax not documented, D-FRC-04)",
}


def _parse_catalogue(text: str) -> Dict[str, CommandSpec]:
    out: Dict[str, CommandSpec] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        name, abbr, tier, cls = parts[0].upper(), parts[1], parts[2], parts[3]
        spec = CommandSpec(name=name, tier=tier, cls=cls)
        if abbr != "-":
            spec.abbrev = tuple(a.upper() for a in abbr.split(","))
        for opt in parts[4:]:
            if opt.startswith("text="):
                spec.text_from = int(opt[5:])
            elif opt == "raw":
                spec.raw = True
            elif opt == "paren":
                spec.paren = True
            elif opt.startswith("key="):
                spec.key = opt[4:]
            elif opt == "silent":
                spec.silent = True
            else:
                raise ValueError(f"catalogue: bad option {opt!r} for {name}")
        if tier not in TIERS or cls not in CLASSES:
            raise ValueError(f"catalogue: bad tier/class for {name}")
        spec.note = NOTES.get(name, "")
        out[name] = spec
    return out


CATALOGUE: Dict[str, CommandSpec] = _parse_catalogue(_CATALOGUE_TEXT)

#: abbreviations of requirements section 3.2 (used by UT-01; a subset of CATALOGUE abbreviations)
ABBREVIATIONS: Dict[str, str] = {
    "ACCIN": "ACCI", "ACCOUT": "ACCO", "AFWRITE": "AFWR", "ANALYS": "ANAL", "CHECK": "CHEC", "EQTIT": "EQTI",
    "EQUAKE": "EQUA", "FORCE": "FORC", "ECOMPR": "ECOM", "GROUP": "GROU", "MRGEN": "MRGE", "MTDEL": "MTDE",
    "MUNITS": "MUNI", "HOUSE": "HOUS", "INCOH": "INCO", "LFREQ": "LFRE", "MOTION": "MOTI", "POINT": "POIN",
    "RESUME": "RESU", "RSOUT": "RSOU", "STATUS": "STAT", "ELIST": "ELIS", "LLIST": "LLIS", "MRDEL": "MRDE",
    "MTGEN": "MTGE", "THFILE": "THFI", "THTIT": "THTI", "WPASS": "WPAS", "WRITE": "WRIT", "GLOBAL": "GLOB",
    "LOCAL": "LOCA", "NLIST": "NLIS", "NMOVE": "NMOV", "NSCALE": "NSCA", "GLIST": "GLIS", "MMDEL": "MMDE",
    "MRSCALE": "MRSC", "MTLIST": "MTLI", "MLIST": "MLIS", "MTYPE": "MTYP", "MXDEL": "MXDE", "MXLIST": "MXLI",
    "RLIST": "RLIS", "SCLIST": "SCLI", "THICK": "THIC", "FLIST": "FLIS", "FSCALE": "FSCA", "MMLIST": "MMLI",
    "MSCALE": "MSCA", "MTSCALE": "MTSC",
    # legacy aliases of D-PAR-04
    "INTLIST": "INTL", "LMOVE": "LMOV", "SLIST": "SLIS", "DELSC": "DELS", "ETYPE": "ETYP", "STRESS": "STRE",
}

# --------------------------------------------------------------------------------------
# Registry
# --------------------------------------------------------------------------------------
REGISTRY: Dict[str, CommandSpec] = {}
_ALIASES: Dict[str, str] = {}
LOAD_ERRORS: List[Tuple[str, str]] = []     # (module, traceback) of handler modules that failed
_loaded = False


def _merge(name: str, opts: dict) -> CommandSpec:
    base = CATALOGUE.get(name)
    spec = replace(base) if base is not None else CommandSpec(name=name)
    for k, v in opts.items():
        if v is not None:
            setattr(spec, k, v)
    return spec


def _index_aliases(spec: CommandSpec) -> None:
    """Index the abbreviations of ``spec``; they must not collide with a name or another abbreviation."""
    for a in spec.abbrev:
        if (a in CATALOGUE or a in REGISTRY) and a != spec.name:
            raise ValueError(f"abbreviation {a} of {spec.name} collides with a command name")
        other = _ALIASES.get(a)
        if other is not None and other != spec.name:
            raise ValueError(f"abbreviation {a} of {spec.name} is already used by {other}")
        _ALIASES[a] = spec.name


def register(name: str, handler: Handler, *, abbrev: Optional[Tuple[str, ...]] = None, tier: Optional[str] = None,
             cls: Optional[str] = None, text_from: Optional[int] = None, raw: Optional[bool] = None,
             paren: Optional[bool] = None, key: Optional[str] = None, max_args: Optional[int] = None,
             summary: str = "", replaces: bool = False, placeholder: bool = False, module: str = "") -> CommandSpec:
    """Register ``handler`` for command ``name`` (see :func:`command`)."""
    name = name.upper()
    old = REGISTRY.get(name)
    reload = (old is not None and old.handler is not None
              and getattr(old.handler, "__qualname__", None) == getattr(handler, "__qualname__", "")
              and old.module == (module or getattr(handler, "__module__", "")))
    if old is not None and not old.placeholder and not replaces and not placeholder and not reload:
        raise ValueError(f"command {name} is already registered by {old.module or '?'} "
                         f"(use replaces=True to override)")
    if old is not None and not old.placeholder and placeholder:
        return old
    spec = _merge(name, dict(abbrev=tuple(a.upper() for a in abbrev) if abbrev else None, tier=tier, cls=cls,
                             text_from=text_from, raw=raw, paren=paren, key=key, max_args=max_args))
    spec.handler = handler
    spec.summary = summary or (handler.__doc__ or "").strip().split("\n")[0]
    spec.placeholder = placeholder
    spec.module = module or getattr(handler, "__module__", "")
    _index_aliases(spec)
    REGISTRY[name] = spec
    return spec


def command(name: str, **opts):
    """Decorator registering a command handler.

    Options (all optional; defaults come from the catalogue): ``abbrev`` (tuple of accepted
    abbreviations), ``tier``, ``cls``, ``text_from`` (1-based rest-of-line text argument),
    ``raw`` (receive the unsubstituted rest of the line), ``paren`` (legacy parenthesised last
    token), ``max_args`` (arguments after it are dropped before the handler runs, with a warning),
    ``summary`` and ``replaces`` (override an existing real handler).
    """
    def deco(fn: Handler) -> Handler:
        register(name, fn, **opts)
        return fn
    return deco


def lookup(token: str) -> Optional[CommandSpec]:
    """Resolve a typed command name (rule L10): full name, else documented abbreviation/alias."""
    load_commands()
    t = token.strip().upper()
    spec = REGISTRY.get(t)
    if spec is not None:
        return spec
    full = _ALIASES.get(t)
    return REGISTRY.get(full) if full else None


def all_commands() -> List[CommandSpec]:
    load_commands()
    return [REGISTRY[k] for k in sorted(REGISTRY)]


# --------------------------------------------------------------------------------------
# Placeholders
# --------------------------------------------------------------------------------------
def _tier_message(spec: CommandSpec) -> str:
    return f"{spec.name} is not available in this build (tier {spec.tier})"


def _announce(c, spec: CommandSpec) -> None:
    if spec.note:
        c.warn(spec.note)
    elif not spec.silent and spec.tier != "P0":
        c.warn(_tier_message(spec))


def _ph_action(c) -> None:
    """Placeholder of an action command that is not implemented in this build."""
    spec = c.spec
    c.warn(spec.note or _tier_message(spec))


def _ph_record(c) -> None:
    """Generic record setter: the record is replaced; every positional argument is stored."""
    from ..model.options import make_record
    spec = c.spec
    legacy = c.lexed.legacy_paren if spec.paren else ""
    rec = make_record(spec.name, c.tokens, legacy=legacy)
    if legacy:
        c.warn(f"{spec.name}: legacy (PREP) format {legacy} detected; the arguments are stored but not "
               f"mapped to the current argument order (Error 102)")
    c.model.options.set_record(rec)
    _announce(c, spec)


def _ph_record_keep(c) -> None:
    """Generic record setter where a blank field leaves the stored value unchanged (BINOUT)."""
    from ..model.options import make_record
    spec = c.spec
    old = c.model.options.record(spec.name)
    rec = old.copy() if old is not None else make_record(spec.name)
    for k, tok in enumerate(c.tokens, start=1):
        if tok.strip():
            rec.set_arg(k, tok)
    c.model.options.set_record(rec)
    _announce(c, spec)


def indexed_key(spec: CommandSpec, tokens: List[str]):
    """Key of an indexed setter entry from its tokens (see the catalogue ``key`` option)."""
    from ..model.values import NumberError, parse_int

    def _int(i: int, what: str) -> int:
        t = tokens[i].strip() if i < len(tokens) else ""
        if not t:
            raise CommandError(f"{spec.name}: {what} is required")
        try:
            return parse_int(t)[0]
        except NumberError:
            raise CommandError(f"{spec.name}: {what} '{t}' is not an integer") from None

    if spec.key == "label_no":
        label = tokens[5].strip() if len(tokens) > 5 else ""
        if not label:
            raise CommandError(f"{spec.name}: label (argument 6) is required")
        return (label, _int(0, "point number"))
    if spec.key == "arg12":
        return (_int(0, "argument 1"), _int(1, "argument 2"))
    if spec.key == "word":
        t = tokens[0].strip() if tokens else ""
        if not t:
            raise CommandError(f"{spec.name}: key is required")
        return t.upper()
    t = tokens[0].strip() if tokens else ""
    if not t:
        raise CommandError(f"{spec.name}: index (argument 1) is required")
    try:
        return parse_int(t)[0]
    except NumberError:
        return t.upper()


def _ph_indexed(c) -> None:
    """Generic indexed setter: creates or overwrites the entry at its key."""
    from ..model.options import make_record
    from ..model.values import NumberError, parse_int
    spec = c.spec
    toks = list(c.tokens)
    key = indexed_key(spec, toks)
    if spec.key == "symm":
        # SYMM,<no>,[<type>],[<node1>]...: node1 = 0 (or blank) resets plane <no>
        n1 = toks[2].strip() if len(toks) > 2 else ""
        try:
            reset = (not n1) or parse_int(n1)[0] == 0
        except NumberError:
            reset = False
        if reset:
            existed = c.model.options.delete_entry(spec.name, key)
            c.confirm(f"{spec.name} {key} reset" if existed else f"{spec.name} {key} not defined (nothing reset)")
            return
    c.model.options.set_entry(spec.name, key, make_record(spec.name, toks))
    _announce(c, spec)


def _ph_string(c) -> None:
    """Generic string setter (THFILE, THTIT, EQTIT, RELFILE ...)."""
    spec = c.spec
    value = c.tokens[0] if c.tokens else ""
    c.model.options.set_string(spec.name, value)
    _announce(c, spec)


_PLACEHOLDERS = {"record": _ph_record, "record_keep": _ph_record_keep, "indexed": _ph_indexed,
                 "string": _ph_string}


def install_placeholders() -> None:
    for name, base in CATALOGUE.items():
        if name in REGISTRY:
            continue
        fn = _PLACEHOLDERS.get(base.cls, _ph_action)
        register(name, fn, placeholder=True, module=__name__)


def load_commands(force: bool = False) -> None:
    """Import every handler module of :mod:`sassi.prep.commands`, then install placeholders.

    A module that fails to import is recorded in :data:`LOAD_ERRORS` (the interpreter reports it)
    and does not stop the others.
    """
    global _loaded
    if _loaded and not force:
        return
    _loaded = True
    from . import commands as pkg
    for m in sorted(pkgutil.iter_modules(pkg.__path__), key=lambda m: m.name):
        modname = f"{pkg.__name__}.{m.name}"
        try:
            importlib.import_module(modname)
        except Exception:  # pragma: no cover - depends on other packages' modules
            LOAD_ERRORS.append((modname, traceback.format_exc()))
    install_placeholders()
