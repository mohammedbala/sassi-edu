"""The built-in input library (``sassi/data/library``; requirements section 7.19, D-W5-01 / D-W5-02).

Standard inputs that ship with the package -- the RG 1.60 0.30 g record and spectrum, the RG 1.60 1 g
spectra, the SRP 3.7.1 Appendix A target PSDs, a 5 Hz Ricker load pulse and the SHAKE91 soil curves --
so that a model can be analysed before its project inputs exist.  A file of the library is named with
the prefix ``@`` wherever a file name is read::

    THFILE,@rg160h_030g.acc          RSIN,1,@rg160h_030g.rsi          INP,@dynp_library.pre

The name resolves to the package folder in every installation: a source checkout, a pip install
(package data ``data/**/*``, ``pyproject.toml``) and the browser version (Pyodide, where the package
is unpacked under ``/home/pyodide/sassi-edu/sassi``) -- the folder is found from this module's own
location.  Library files are read-only: an ``@`` name is never an output file.

Every reader of a file name uses :func:`resolve` (or :func:`library_path`): the interpreter
(:meth:`sassi.prep.interpreter.Interpreter.resolve_path`), CHECK (:func:`sassi.prep.check.resolve_file`),
AFWRITE (decks keep the portable ``@`` name) and the modules (:func:`module_path`).  :data:`CATALOGUE`
describes each file (kind, units, source); the defaults of blank inputs that use these files are in
:mod:`sassi.prep.defaults`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

#: the library folder inside the package
LIBRARY_DIR = Path(__file__).resolve().parents[1] / "data" / "library"
#: prefix of a library file name
PREFIX = "@"
#: kinds of library files (and the input fields they fit)
KINDS = {
    "record": "acceleration record (THFILE, ACCIN, SOILX file)",
    "load": "load history (THFILE of a forced-vibration analysis)",
    "spectrum": "target response spectrum (RSIN)",
    "psd": "target power spectral density (TPSD)",
    "dynp": "strain-dependent soil curves (DYNP; INP)",
}
_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]*$")


@dataclass(frozen=True)
class LibraryFile:
    """One file of the library."""
    name: str                     # file name in LIBRARY_DIR (the @ name without the prefix)
    kind: str                     # key of KINDS
    title: str                    # one line, e.g. "RG 1.60 horizontal record, 0.30 g"
    units: str
    source: str
    dt: Optional[float] = None    # time step of a history (MOTION <fopt> 0 layout: dt on line 1)
    values: int = 0               # values of a history / rows of a two-column file
    short: str = ""               # the short form used in the default-report lines

    @property
    def ref(self) -> str:
        """The name to type: ``@<name>``."""
        return PREFIX + self.name

    @property
    def path(self) -> Path:
        return LIBRARY_DIR / self.name

    @property
    def fields(self) -> Tuple[str, ...]:
        return {"record": ("THFILE", "ACCIN", "SOILX"), "load": ("THFILE",), "spectrum": ("RSIN",),
                "psd": ("TPSD",), "dynp": ("INP",)}[self.kind]

    def describe(self) -> str:
        """``@name: title (units)``."""
        return f"{self.ref}: {self.title} ({self.units})"

    def to_dict(self) -> Dict[str, object]:
        return {"name": self.ref, "file": self.name, "kind": self.kind, "title": self.title, "units": self.units,
                "source": self.source, "dt": self.dt, "values": self.values, "fields": list(self.fields)}


CATALOGUE: Tuple[LibraryFile, ...] = (
    LibraryFile("rg160h_030g.acc", "record", "RG 1.60 horizontal record matched by EQUAKE, PGA 0.324 g, 20 s",
                "g; dt 0.005 s, 4001 values, MOTION <fopt> 0 layout",
                "SASSI-EDU EQUAKE (seed 11975) matched to rg160h_030g.rsi; examples/data/rg160h_030g.acc",
                dt=0.005, values=4001, short="RG 1.60 record (0.30 g, 20 s)"),
    LibraryFile("rg160h_030g.rsi", "spectrum", "RG 1.60 horizontal design spectrum anchored to 0.30 g, 5 % damping",
                "Hz, SA in g; 27 rows, 0.1-100 Hz",
                "US NRC RG 1.60 (sassi.core.equake_lib.rg160_spectrum); examples/data/rg160h_030g.rsi",
                values=27, short="RG 1.60 H spectrum (0.30 g, 5 %, 27 frequencies)"),
    LibraryFile("rg160h_1g.rsi", "spectrum", "RG 1.60 horizontal design spectrum anchored to 1.0 g, 5 % damping",
                "Hz, SA in g; 27 rows, 0.1-100 Hz; scale linearly for another PGA",
                "US NRC RG 1.60; sassi/data/rg160/RG160H_5pct_1g.rsi", values=27,
                short="RG 1.60 H spectrum (1 g, 5 %)"),
    LibraryFile("rg160v_1g.rsi", "spectrum", "RG 1.60 vertical design spectrum anchored to 1.0 g, 5 % damping",
                "Hz, SA in g; 27 rows, 0.1-100 Hz; scale linearly for another PGA",
                "US NRC RG 1.60; sassi/data/rg160/RG160V_5pct_1g.rsi", values=27,
                short="RG 1.60 V spectrum (1 g, 5 %)"),
    LibraryFile("rg160h_030g_cm2s3.tpsd", "psd", "SRP 3.7.1 App. A minimum PSD, RG 1.60 horizontal, 0.30 g, SI",
                "Hz, one-sided PSD in cm^2/s^3; 124 rows, 0.05-100 Hz",
                "SRP 3.7.1 Rev. 4 Appendix A (sassi.core.equake_lib.rg160_target_psd, pga 0.30)", values=124,
                short="SRP App. A PSD (0.30 g, SI)"),
    LibraryFile("rg160h_030g_in2s3.tpsd", "psd", "SRP 3.7.1 App. A minimum PSD, RG 1.60 horizontal, 0.30 g, British",
                "Hz, one-sided PSD in in^2/s^3; 124 rows, 0.05-100 Hz",
                "SRP 3.7.1 Rev. 4 Appendix A (sassi.core.equake_lib.rg160_target_psd, pga 0.30)", values=124,
                short="SRP App. A PSD (0.30 g, British)"),
    LibraryFile("rg160h_appa_1g_cm2s3.tpsd", "psd", "SRP 3.7.1 App. A minimum PSD, RG 1.60 horizontal, 1.0 g, SI",
                "Hz, one-sided PSD in cm^2/s^3; 124 rows; scale by PGA^2",
                "SRP 3.7.1 Rev. 4 Appendix A; sassi/data/rg160/RG160H_AppA_1g_cm2s3.tpsd", values=124,
                short="SRP App. A PSD (1 g, SI)"),
    LibraryFile("rg160h_appa_1g_in2s3.tpsd", "psd", "SRP 3.7.1 App. A minimum PSD, RG 1.60 horizontal, 1.0 g, British",
                "Hz, one-sided PSD in in^2/s^3; 124 rows; scale by PGA^2",
                "SRP 3.7.1 Rev. 4 Appendix A; sassi/data/rg160/RG160H_AppA_1g_in2s3.tpsd", values=124,
                short="SRP App. A PSD (1 g, British)"),
    LibraryFile("ricker_5hz.th", "load", "5 Hz Ricker wavelet, peak 1 at t = 0.5 s, 2 s",
                "load factor (dimensionless); dt 0.005 s, 401 values, MOTION <fopt> 0 layout",
                "SASSI-EDU; examples/data/ricker_5hz.th", dt=0.005, values=401,
                short="5 Hz Ricker pulse (peak 1, 2 s)"),
    LibraryFile("dynp_library.pre", "dynp", "SHAKE91 soil curves Clay, Sand and Rock (DYNP commands)",
                "strain %, G/Gmax, damping %",
                "SHAKE91 INP.DAT (Idriss & Sun 1992, public domain): Seed & Sun 1989 / Seed & Idriss 1970 / "
                "Schnabel 1973 with Idriss 1990 damping; sassi/data/dynp_library.pre",
                short="SHAKE91 soil curves"),
)

#: DYNP labels of dynp_library.pre and their sources (the curves themselves: :func:`dynp_curves`)
DYNP_LABELS = {"Clay": "Seed & Sun (1989) upper range, Idriss (1990) damping",
               "Sand": "Seed & Idriss (1970) upper range, Idriss (1990) damping",
               "Rock": "Schnabel (1973) average"}
DYNP_FILE = "dynp_library.pre"


# ======================================================================================
# names and resolution
# ======================================================================================
def is_library_name(name: Union[str, Path, None]) -> bool:
    """True for a name of the form ``@<file>`` (surrounding blanks and quotes ignored)."""
    if name is None:
        return False
    s = str(name).strip().strip('"').strip()
    return s.startswith(PREFIX) and len(s) > 1


def bare(name: str) -> str:
    """The file name of an ``@`` name (prefix, blanks and quotes removed)."""
    s = str(name).strip().strip('"').strip()
    return s[len(PREFIX):] if s.startswith(PREFIX) else s


def entry(name: Union[str, Path, None]) -> Optional[LibraryFile]:
    """The catalogue entry of a library name (``@x`` or ``x``; case-insensitive), or None."""
    if name is None:
        return None
    b = bare(str(name)).lower()
    for e in CATALOGUE:
        if e.name.lower() == b:
            return e
    return None


def library_path(name: Union[str, Path, None]) -> Optional[Path]:
    """Path of the library file ``@name`` (a :data:`CATALOGUE` name, case-insensitive), or None when ``name``
    is not one.  Names with a path separator or ``..`` are never library files."""
    if not is_library_name(name):
        return None
    b = bare(str(name))
    if not _NAME_RE.match(b) or ".." in b:
        return None
    e = entry(b)
    return e.path if e is not None and e.path.is_file() else None


def names() -> List[str]:
    """The ``@`` names of the catalogue."""
    return [e.ref for e in CATALOGUE]


def starts_with_name(text: str) -> int:
    """Length of the library file name at the start of ``text`` (the text after an ``@``), 0 when none.

    The interpreter uses it to leave ``@rg160h_030g.acc`` alone instead of reading ``@rg160h_030g`` as a
    variable (D-W5-02): the longest catalogue name that matches case-insensitively and is followed by
    the end of the text or a delimiter."""
    low = text.lower()
    best = 0
    for e in CATALOGUE:
        n = len(e.name)
        if low[:n] == e.name.lower() and (len(text) == n or not (text[n].isalnum() or text[n] in "_.-")):
            best = max(best, n)
    return best


def resolve(name: str, dirs: Sequence[Union[str, Path]]) -> Optional[Path]:
    """Resolve an input file name like CHECK does: ``@name`` from the library, an absolute path, else
    the first directory of ``dirs`` that holds it.  None when not found."""
    if is_library_name(name):
        return library_path(name)
    raw = str(name).strip()
    if not raw:
        return None
    p = Path(raw)
    if p.is_absolute():
        return p if p.exists() else None
    for d in dirs:
        if d and (Path(d) / p).exists():
            return Path(d) / p
    return None


def module_path(name: str, workdir: Union[str, Path]) -> Path:
    """Input path of a deck file name for a module run in ``workdir``: an ``@`` name is the library
    file (or, when unknown, a path that does not exist, so the module's "not found" message names it),
    an absolute path is kept, a relative one is taken in the model directory."""
    raw = str(name).strip().strip('"')
    if is_library_name(raw):
        p = library_path(raw)
        return p if p is not None else Path(workdir) / raw
    p = Path(raw)
    return p if p.is_absolute() else Path(workdir) / p


def note(name: Union[str, Path, None]) -> str:
    """Listing note for a library input, ``"built-in input @x: title (units)"``; '' for other names."""
    if not is_library_name(name):
        return ""
    e = entry(name)
    if e is None:
        return f"built-in input {str(name).strip()}: not in the library (LIBRARY lists the built-in inputs)"
    return f"built-in input {e.describe()}; source {e.source}"


def history_dt(name: Union[str, Path, None]) -> Optional[float]:
    """Time step of a library history (dt on its first line), None for other names."""
    e = entry(name) if is_library_name(name) else None
    return e.dt if e is not None else None


def entries(kind: Optional[str] = None) -> List[LibraryFile]:
    return [e for e in CATALOGUE if kind is None or e.kind == kind]


def catalogue_json() -> List[Dict[str, object]]:
    """The catalogue for the GUI (dialog Library picker)."""
    return [e.to_dict() for e in CATALOGUE if e.path.is_file()]


# ======================================================================================
# DYNP curves
# ======================================================================================
@dataclass
class DynpPoint:
    """One DYNP command of the library: ``DYNP,<no>,<sg>,<g>,<sd>,<d>,<label>`` (values as text)."""
    no: int
    sg: str
    g: str
    sd: str
    d: str
    label: str

    def tokens(self) -> List[str]:
        return [str(self.no), self.sg, self.g, self.sd, self.d, self.label]


@lru_cache(maxsize=1)
def _dynp_points() -> Tuple[DynpPoint, ...]:
    p = LIBRARY_DIR / DYNP_FILE
    out: List[DynpPoint] = []
    try:
        text = p.read_text(encoding="utf-8")
    except OSError:
        return ()
    for ln in text.splitlines():
        s = ln.strip()
        if not s.upper().startswith("DYNP,"):
            continue
        t = [x.strip() for x in s.split(",")]
        if len(t) < 7:
            continue
        out.append(DynpPoint(int(t[1]), t[2], t[3], t[4], t[5], ",".join(t[6:])))
    return tuple(out)


def dynp_labels() -> List[str]:
    """Labels of the library curves in file order (Clay, Sand, Rock)."""
    return list(dict.fromkeys(pt.label for pt in _dynp_points()))


def dynp_points(label: str) -> List[DynpPoint]:
    """The DYNP points of library curve ``label`` (case-sensitive, like DYNP labels), [] when unknown."""
    return sorted((pt for pt in _dynp_points() if pt.label == label), key=lambda pt: pt.no)


def dynp_curves() -> Dict[str, List[Dict[str, float]]]:
    """``{label: [{no, sg, g, sd, d}, ...]}`` of the library curves (GUI, plots)."""
    out: Dict[str, List[Dict[str, float]]] = {}
    for lab in dynp_labels():
        out[lab] = [{"no": pt.no, "sg": float(pt.sg), "g": float(pt.g), "sd": float(pt.sd), "d": float(pt.d)}
                    for pt in dynp_points(lab)]
    return out


__all__ = ["LIBRARY_DIR", "PREFIX", "KINDS", "LibraryFile", "CATALOGUE", "DYNP_LABELS", "DYNP_FILE",
           "is_library_name", "bare", "entry", "library_path", "names", "starts_with_name", "resolve",
           "module_path", "note", "history_dt", "entries", "catalogue_json", "DynpPoint", "dynp_labels",
           "dynp_points", "dynp_curves"]
