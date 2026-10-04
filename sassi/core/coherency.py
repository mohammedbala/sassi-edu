"""Unlagged coherency models of incoherent seismic input (HOUSE, requirements 4.4 item 6; decisions
D-INC-01, D-INC-04, D-INC-05, D-INC-11; spec 05b sections 3.9-3.16; ACS SASSI manual section 6.5.4).

What a coherency function is (for the structural engineer)
---------------------------------------------------------
Two points of the free-field ground surface a horizontal distance ``D`` apart do not move exactly
alike: scattered waves make the motion partly *incoherent*.  The (plane-wave) coherency
``gamma(f, D)`` is the correlation coefficient of the Fourier amplitudes of the two motions at
frequency ``f`` after the wave-passage time shift has been removed.  It is 1 for ``D = 0`` or
``f = 0`` (identical motions) and decays with frequency and separation.  HOUSE builds, for every SSI
frequency and motion component, the coherency matrix ``Sigma_ij = gamma(f, D_ij)`` of the
interaction nodes (:mod:`sassi.core.incoherency`); wave passage (the time lag along Line D) is
applied separately (WPASS), so the functions below are *unlagged* (requirements glossary).

Models (``WPASS <cohf>``, D-INC-01)
-----------------------------------
====  ===================================================  ==========================================
cohf  model                                                status in SASSI-EDU
====  ===================================================  ==========================================
1     Luco and Wong (1986) ``exp[-(gamma_c w D / Vs)^2]``    available (D-INC-05; w in rad/s)
2     Abrahamson (1993), all soil types                    coefficients not available (refused)
3     Abrahamson (2005), all sites, surface foundations    available (EPRI 1012968)
4     Abrahamson (2006), all sites, embedded foundations   coefficients not available (refused)
5     Abrahamson (2007), hard rock                         available (EPRI 1015110, chapter 6)
6     Abrahamson (2007), soil sites, surface foundations   available (EPRI 1015110, chapter 7)
7     user tables COHXUSER/COHYUSER/COHZUSER, FREQCOH,     available (D-INC-11)
      DISTCOH
====  ===================================================  ==========================================

The Abrahamson coefficients are **data**: ``sassi/data/coherency/*.json``, each file with the
citation (report, equation, table and page) and example values read from the published figures,
which the unit tests reproduce (D-INC-04: never from memory).  A model whose file is marked
``"verified": false`` refuses to run with a "coefficients not available" error.  All Abrahamson
models share the plane-wave form (EPRI 1012968 Eq. 3-1; EPRI 1015110 Eqs. 5-2, 6-1, 7-1)::

    gamma_pw(f, xi) = [1 + (f tanh(a3 xi) / (a1 fc(xi)))^n1]^(-1/2) [1 + (f tanh(a3 xi) / (a2 [fc(xi)]))^n2]^(-1/2)

with ``f`` in Hz and the separation ``xi`` in metres (the 2005 model has ``a2 fc(xi)`` in the second
factor, the 2007 models ``a2`` only), separate horizontal (X, Y) and vertical (Z) coefficient sets.
Model distances are converted to metres with the unit system of the model (D-CNV-08: g > 20 ->
feet).  Separations beyond the published range of a model (``formula_max_m``: 150 m for the 2007
models, whose coefficient expressions turn non-physical beyond it -- e.g. the soil-site ``a2`` is
negative above 359 m) are evaluated at ``formula_max_m``: the coherency is held at its value at the
end of the range, the conservative (more coherent) choice; HOUSE lists a warning.

Distances (requirements 4.4 item 6, spec 05b 3.10): model 1 uses the plain horizontal distance;
models 2-7 the directional distance ``D = sqrt(2 (alpha dX'^2 + (1 - alpha) dY'^2))`` with ``dX'``,
``dY'`` measured along and across Line D (``WPASS <ang>``); ``alpha = 0.5`` is the isotropic case,
``alpha = 0.1`` weights Y' separations 3 times more than X' separations (the formula is
authoritative, conflict C-22).
"""
from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

#: directory of the coefficient files (requirements D-INC-04)
DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "coherency"
#: coefficient / citation file of each model number
MODEL_FILES = {1: "luco_wong_1986.json", 2: "abrahamson_1993.json", 3: "abrahamson_2005_all_sites.json",
               4: "abrahamson_2006_embedded.json", 5: "abrahamson_2007_hard_rock.json",
               6: "abrahamson_2007_soil.json"}
#: dialog names of the "Unlagged Coherency Model" types (ACS SASSI manual section 6.5.4)
MODEL_NAMES = {1: "1986 Luco-Wong", 2: "1993 Abrahamson, all soil types",
               3: "2005 Abrahamson, all sites, surface foundations",
               4: "2006 Abrahamson, all sites, embedded foundations",
               5: "2007 Abrahamson, hard-rock sites, surface/embedded foundations",
               6: "2007 Abrahamson, soil sites, surface foundations",
               7: "user-defined coherency functions (COHXUSER, COHYUSER, COHZUSER, FREQCOH, DISTCOH)"}
#: user coherency files (model 7, D-INC-11): one table per motion component X, Y, Z
USER_TABLES = ("COHXUSER", "COHYUSER", "COHZUSER")
USER_FREQ, USER_DIST = "FREQCOH", "DISTCOH"
DIRECTIONS = ("X", "Y", "Z")
FT_TO_M = 0.3048


class CoherencyError(ValueError):
    """Invalid coherency input, or a model whose coefficients are not available."""


# ======================================================================================
# Distances
# ======================================================================================
def horizontal_distances(xy) -> np.ndarray:
    """Plain horizontal distances ``|r_i - r_j|`` of plan points ``xy`` (n, 2) -> (n, n)."""
    p = np.asarray(xy, float).reshape(-1, 2)
    d = p[:, None, :] - p[None, :, :]
    return np.sqrt(np.einsum("ijk,ijk->ij", d, d))


def directional_distances(xy, alpha: float, angle_deg: float = 0.0) -> np.ndarray:
    """Directional distance of models 2-7 (spec 05b 3.10, manual "homotopic relationship")::

        D_ij = sqrt(2 (alpha dX'^2 + (1 - alpha) dY'^2))

    ``dX'``, ``dY'`` are the separations along and across Line D, the horizontal line at
    ``angle_deg`` from the global X axis (WPASS <ang>).  ``alpha = 0.5`` gives the Euclidean
    distance; for unit separations along X' and Y' the distances are ``sqrt(2 alpha)`` and
    ``sqrt(2 (1 - alpha))`` (VP-27)."""
    p = np.asarray(xy, float).reshape(-1, 2)
    a = math.radians(float(angle_deg))
    ca, sa = math.cos(a), math.sin(a)
    xp = p[:, 0] * ca + p[:, 1] * sa
    yp = -p[:, 0] * sa + p[:, 1] * ca
    dx = xp[:, None] - xp[None, :]
    dy = yp[:, None] - yp[None, :]
    return np.sqrt(2.0 * (alpha * dx * dx + (1.0 - alpha) * dy * dy))


# ======================================================================================
# Model 1: Luco-Wong
# ======================================================================================
def luco_wong(f, D, gamma_c: float, vs: float) -> np.ndarray:
    """Luco and Wong (1986) coherency ``exp[-(gamma_c w D / Vs)^2]``, ``w = 2 pi f`` in rad/s
    (D-INC-05; EPRI/NRC 2007 CLASSI incoherency session ML072620217, slide 6)."""
    if not vs > 0:
        raise CoherencyError(f"Luco-Wong model: mean shear-wave velocity <alpha> = {vs} must be > 0 (Error 59)")
    w = 2.0 * math.pi * np.asarray(f, float)
    z = gamma_c * w * np.asarray(D, float) / vs
    return np.exp(-z * z)


# ======================================================================================
# Models 2-6: Abrahamson plane-wave coherency (coefficient files)
# ======================================================================================
def _expr(spec: dict, xi: np.ndarray) -> np.ndarray:
    """Evaluate a coefficient expression of a coefficient file at separations ``xi`` (m)."""
    if len(spec) != 1:
        raise CoherencyError(f"coefficient expression {spec!r} must have exactly one kind")
    kind, c = next(iter(spec.items()))
    xi = np.asarray(xi, float)
    if kind == "const":
        return np.full(xi.shape, float(c))
    c = [float(v) for v in c]
    if kind == "lin":                                   # c0 + c1 xi
        return c[0] + c[1] * xi
    if kind == "lnq":                                   # c0 + c1 ln(xi+1) + c2 [ln(xi+1) - c3]^2
        L = np.log(xi + 1.0)
        return c[0] + c[1] * L + c[2] * (L - c[3]) ** 2
    if kind == "exp_lnq":                               # exp(c0 + c1 ln(xi+1) + c2 [ln(xi+1) - c3]^2)
        L = np.log(xi + 1.0)
        return np.exp(c[0] + c[1] * L + c[2] * (L - c[3]) ** 2)
    if kind == "ln_shift":                              # c0 + c1 ln(xi + c2)
        return c[0] + c[1] * np.log(xi + c[2])
    if kind == "ln_inv":                                # c0 + c1 ln(c2/(xi+1) + c3)
        return c[0] + c[1] * np.log(c[2] / (xi + 1.0) + c[3])
    if kind == "exp_lin":                               # exp(c0 + c1 xi)
        return np.exp(c[0] + c[1] * xi)
    raise CoherencyError(f"unknown coefficient expression kind {kind!r}")


@dataclass
class PreparedPW:
    """Separation-dependent coefficients of one component at a set of separations (evaluated once
    per coherency matrix pattern, then reused at every frequency)."""

    t: np.ndarray          # tanh(a3 xi)
    d1: np.ndarray         # a1 fc
    d2: np.ndarray         # a2 [fc]
    n1: np.ndarray
    n2: np.ndarray


@dataclass
class AbrahamsonModel:
    """A plane-wave coherency model read from its coefficient file (D-INC-04)."""

    number: int
    name: str
    data: dict
    path: Path

    @property
    def second_uses_fc(self) -> bool:
        return bool(self.data.get("second_factor_uses_fc", False))

    @property
    def formula_max_m(self) -> float:
        return float(self.data.get("formula_max_m", 150.0))

    @property
    def validated_range_m(self) -> Tuple[float, float]:
        lo, hi = self.data.get("validated_range_m", [0.0, 150.0])
        return float(lo), float(hi)

    @property
    def citation(self) -> str:
        s = self.data.get("source", {})
        return "; ".join(str(s[k]) for k in ("primary", "equation", "tables") if k in s)

    def prepare(self, xi_m, component: str) -> PreparedPW:
        """Coefficients at separations ``xi_m`` (metres; clamped to ``formula_max_m``)."""
        comp = self.data["components"][component]
        xi = np.minimum(np.asarray(xi_m, float), self.formula_max_m)
        a1, a2, a3 = (_expr(comp[k], xi) for k in ("a1", "a2", "a3"))
        fc = _expr(comp["fc"], xi)
        return PreparedPW(t=np.tanh(a3 * xi), d1=a1 * fc, d2=a2 * fc if self.second_uses_fc else a2,
                          n1=_expr(comp["n1"], xi), n2=_expr(comp["n2"], xi))

    @staticmethod
    def evaluate_prepared(f: float, p: PreparedPW) -> np.ndarray:
        """``gamma_pw`` at frequency ``f`` (Hz) from prepared coefficients."""
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            x1 = (float(f) * p.t) / p.d1
            x2 = (float(f) * p.t) / p.d2
            g1 = np.power(1.0 + np.power(x1, p.n1), -0.5)
            g2 = np.power(1.0 + np.power(x2, p.n2), -0.5)
        out = g1 * g2
        out[~np.isfinite(out)] = 0.0
        return out

    def __call__(self, f, xi_m, component: str = "horizontal") -> np.ndarray:
        """Plane-wave coherency at frequencies ``f`` (Hz) and separations ``xi_m`` (m), broadcast."""
        f_b, xi_b = np.broadcast_arrays(np.asarray(f, float), np.asarray(xi_m, float))
        p = self.prepare(xi_b, component)
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            x1 = f_b * p.t / p.d1
            x2 = f_b * p.t / p.d2
            out = np.power(1.0 + np.power(x1, p.n1), -0.5) * np.power(1.0 + np.power(x2, p.n2), -0.5)
        out = np.where(np.isfinite(out), out, 0.0)
        return out if out.ndim else float(out)


def model_file(number: int) -> Path:
    if number not in MODEL_FILES:
        raise CoherencyError(f"coherency model {number} has no coefficient file")
    return DATA_DIR / MODEL_FILES[number]


def model_info(number: int) -> dict:
    """The coefficient / citation file of model ``number`` (1-6) as a dict."""
    p = model_file(number)
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise CoherencyError(f"coherency model {number}: cannot read {p.name}: {exc}") from None


def model_available(number: int) -> bool:
    """True when model ``number`` can be evaluated (models 1 and 7 always; 2-6 when their
    coefficient file is marked verified)."""
    if number in (1, 7):
        return True
    try:
        return bool(model_info(number).get("verified")) and model_info(number).get("components") is not None
    except CoherencyError:
        return False


def unavailable_reason(number: int) -> str:
    """The 'coefficients not available' message of a refused model (D-INC-04)."""
    try:
        info = model_info(number)
    except CoherencyError as exc:
        return str(exc)
    msg = (f"coherency model {number} ({MODEL_NAMES.get(number, '?')}): "
           + str(info.get("reason", "coefficients not available")))
    alt = info.get("alternative")
    return msg + (f" -- use instead: {alt}" if alt else "")


@lru_cache(maxsize=None)
def load_model(number: int) -> AbrahamsonModel:
    """The Abrahamson plane-wave model ``number`` (2-6) from its coefficient file.

    Raises :class:`CoherencyError` with the 'coefficients not available' explanation when the file
    is marked unverified (D-INC-04).  The coefficient expressions are checked to be positive and
    finite over the separation range ``[0, formula_max_m]`` (a guard against transcription errors)."""
    if number not in (2, 3, 4, 5, 6):
        raise CoherencyError(f"coherency model {number} is not an Abrahamson model (2-6)")
    info = model_info(number)
    if not info.get("verified") or info.get("components") is None:
        raise CoherencyError(unavailable_reason(number))
    m = AbrahamsonModel(number, str(info.get("name", MODEL_NAMES[number])), info, model_file(number))
    xi = np.linspace(0.0, m.formula_max_m, 301)
    for comp in ("horizontal", "vertical"):
        p = m.prepare(xi, comp)
        for nm, a in (("a1 fc", p.d1), ("a2 [fc]", p.d2), ("n1", p.n1), ("n2", p.n2)):
            if not np.all(np.isfinite(a)) or np.any(a <= 0):
                raise CoherencyError(f"coherency model {number} ({m.path.name}): {comp} {nm} is not positive over "
                                     f"0..{m.formula_max_m:g} m -- check the coefficient file")
    return m


# ======================================================================================
# Model 7: user tables
# ======================================================================================
_NUM = re.compile(r"[,\s;]+")


def _read_numbers(path: Path) -> List[float]:
    """Numbers of a free-format text file (comma/blank separated; lines starting with '#', '*' or
    '!' are comments; Fortran 'D' exponents accepted)."""
    out: List[float] = []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s[0] in "#*!":
            continue
        for tok in _NUM.split(s):
            if tok:
                try:
                    out.append(float(tok.replace("D", "E").replace("d", "e")))
                except ValueError:
                    raise CoherencyError(f"{path.name}: '{tok}' is not a number") from None
    return out


@dataclass
class UserCoherency:
    """User-defined coherency tables of model 7 (D-INC-11).

    ``freq`` (NF,) Hz ascending, ``dist`` (ND,) separations in model length units ascending, and one
    table (NF, ND) per component X, Y, Z (rows = frequencies, columns = distances).  Evaluation is
    bilinear in (f, D), clamped at the table ends; ``gamma(f, 0) = 1``: when the first distance is
    larger than 0 the point (0, 1) is added so that the coherency rises linearly to 1 at D = 0."""

    freq: np.ndarray
    dist: np.ndarray
    tables: Tuple[np.ndarray, np.ndarray, np.ndarray]
    notes: List[str] = field(default_factory=list)

    @classmethod
    def from_dir(cls, workdir: Union[str, Path]) -> "UserCoherency":
        wd = Path(workdir)
        missing = [n for n in (USER_FREQ, USER_DIST) + USER_TABLES if not (wd / n).exists()]
        if missing:
            raise CoherencyError("user coherency model 7 needs the files " + ", ".join(missing)
                                 + " in the model directory (no extension; spec 05b 3.16)")
        freq = np.asarray(_read_numbers(wd / USER_FREQ), float)
        dist = np.asarray(_read_numbers(wd / USER_DIST), float)
        if freq.size < 1 or dist.size < 1:
            raise CoherencyError(f"{USER_FREQ} and {USER_DIST} must hold at least one value each")
        if freq.size > 1 and np.any(np.diff(freq) <= 0):
            raise CoherencyError(f"{USER_FREQ}: frequencies must be strictly increasing")
        if dist.size > 1 and np.any(np.diff(dist) <= 0):
            raise CoherencyError(f"{USER_DIST}: distances must be strictly increasing")
        if np.any(freq < 0) or np.any(dist < 0):
            raise CoherencyError(f"{USER_FREQ}/{USER_DIST}: values must be >= 0")
        notes: List[str] = []
        tabs = []
        for name in USER_TABLES:
            vals = np.asarray(_read_numbers(wd / name), float)
            if vals.size != freq.size * dist.size:
                raise CoherencyError(f"{name}: {vals.size} values, {freq.size} x {dist.size} = "
                                     f"{freq.size * dist.size} expected (rows = {USER_FREQ} frequencies, "
                                     f"columns = {USER_DIST} distances)")
            t = vals.reshape(freq.size, dist.size)
            if np.any(np.abs(t) > 1.0 + 1e-9) or not np.all(np.isfinite(t)):
                raise CoherencyError(f"{name}: coherency values must lie in [-1, 1]")
            if dist[0] == 0.0 and np.any(np.abs(t[:, 0] - 1.0) > 1e-12):
                notes.append(f"{name}: coherency at distance 0 set to 1 (D-INC-11)")
                t = t.copy()
                t[:, 0] = 1.0
            tabs.append(t)
        if freq.size * dist.size != 100 * 100:
            notes.append(f"user coherency tables of {freq.size} frequencies x {dist.size} distances "
                         "(the manual's default size is 100 x 100; sizes are taken from the files)")
        return cls(freq, dist, (tabs[0], tabs[1], tabs[2]), notes)

    def _row(self, f: float, direction: int) -> Tuple[np.ndarray, np.ndarray]:
        """Distances and coherencies of the table interpolated at frequency ``f`` (clamped),
        with the (0, 1) anchor."""
        t = self.tables[direction]
        if self.freq.size == 1:
            row = t[0]
        else:
            fc = min(max(float(f), float(self.freq[0])), float(self.freq[-1]))
            k = int(np.searchsorted(self.freq, fc, side="right")) - 1
            k = min(max(k, 0), self.freq.size - 2)
            w = (fc - self.freq[k]) / (self.freq[k + 1] - self.freq[k])
            row = (1.0 - w) * t[k] + w * t[k + 1]
        d = self.dist
        if d[0] > 0.0:
            d = np.concatenate([[0.0], d])
            row = np.concatenate([[1.0], row])
        return d, row

    def evaluate(self, f: float, D, direction: int) -> np.ndarray:
        """Coherency of component ``direction`` (0 X, 1 Y, 2 Z) at frequency ``f`` and distances ``D``."""
        d, row = self._row(f, direction)
        D = np.asarray(D, float)
        out = np.interp(D.reshape(-1), d, row).reshape(D.shape)
        out[D == 0.0] = 1.0
        return out


# ======================================================================================
# The coherency field of one HOUSE run
# ======================================================================================
@dataclass
class CoherencySpec:
    """Everything needed to build the coherency matrices of one HOUSE run.

    ``model`` = WPASS <cohf>; ``gamma`` = INCOH (<gammax>, <gammay>, <gammaz>) (model 1);
    ``alpha`` = INCOH <alpha> (directionality factor, models 2-7; mean Vs, model 1); ``angle_deg``
    = WPASS <ang> (Line D); ``length_to_m`` converts model lengths to metres (models 2-6);
    ``user`` = the model-7 tables."""

    model: int
    gamma: Tuple[float, float, float] = (0.1, 0.1, 0.2)
    alpha: float = 0.5
    angle_deg: float = 0.0
    length_to_m: float = 1.0
    user: Optional[UserCoherency] = None

    def validate(self) -> None:
        """Errors 58, 59, 114 and the availability of the model (D-INC-01, D-INC-04)."""
        if self.model not in MODEL_NAMES:
            raise CoherencyError(f"Error 114: unlagged coherency model <cohf> = {self.model} must be 1..7 (D-INC-01)")
        if self.model == 1:
            bad = [f"{n} = {v:g}" for n, v in zip(("gammax", "gammay", "gammaz"), self.gamma) if not v >= 0.1]
            if bad:
                raise CoherencyError("Error 58: coherence parameter " + ", ".join(bad) + " must be >= 0.1")
            if not self.alpha > 0:
                raise CoherencyError(f"Error 59: Luco-Wong mean shear-wave velocity <alpha> = {self.alpha:g} "
                                     "must be > 0")
        elif not 0.0 <= self.alpha <= 1.0:
            raise CoherencyError(f"Error 114: directionality factor <alpha> = {self.alpha:g} must lie in [0, 1] for "
                                 f"coherency model {self.model} (0.5 isotropic)")
        if self.model in (2, 3, 4, 5, 6):
            load_model(self.model)
        if self.model == 7 and self.user is None:
            raise CoherencyError("user coherency model 7: tables not loaded")

    def distances(self, xy) -> np.ndarray:
        """Separations in model length units: plain horizontal (model 1), directional (2-7)."""
        if self.model == 1:
            return horizontal_distances(xy)
        return directional_distances(xy, self.alpha, self.angle_deg)

    def same_matrix(self, d1: int, d2: int) -> bool:
        """True when components ``d1`` and ``d2`` have identical coherency matrices (lets HOUSE reuse
        one eigen-decomposition for X and Y)."""
        if self.model == 1:
            return self.gamma[d1] == self.gamma[d2]
        if self.model == 7:
            return np.array_equal(self.user.tables[d1], self.user.tables[d2])
        return (d1 == 2) == (d2 == 2)

    def matrix_builder(self, xy, direction: int) -> Callable[[float], np.ndarray]:
        """A function ``f -> Sigma(f)`` (n, n) for component ``direction`` (0 X, 1 Y, 2 Z) of the plan
        points ``xy``; the distance-dependent parts are evaluated once."""
        D = self.distances(xy)
        if self.model == 1:
            g = float(self.gamma[direction])
            vs = float(self.alpha)

            def build(f: float) -> np.ndarray:
                S = luco_wong(f, D, g, vs)
                np.fill_diagonal(S, 1.0)
                return S
            return build
        if self.model == 7:
            user = self.user

            def build(f: float) -> np.ndarray:
                S = user.evaluate(f, D, direction)
                np.fill_diagonal(S, 1.0)
                return S
            return build
        m = load_model(self.model)
        prep = m.prepare(D * self.length_to_m, "vertical" if direction == 2 else "horizontal")

        def build(f: float) -> np.ndarray:
            S = m.evaluate_prepared(f, prep)
            np.fill_diagonal(S, 1.0)
            return S
        return build

    def gamma_at(self, f, D, direction: int) -> np.ndarray:
        """The coherency function itself at frequency ``f`` and model-unit distances ``D``."""
        D = np.asarray(D, float)
        if self.model == 1:
            return luco_wong(f, D, float(self.gamma[direction]), float(self.alpha))
        if self.model == 7:
            return self.user.evaluate(float(f), D, direction)
        m = load_model(self.model)
        return m(f, D * self.length_to_m, "vertical" if direction == 2 else "horizontal")

    def describe(self) -> str:
        name = MODEL_NAMES.get(self.model, "?")
        if self.model == 1:
            return (f"model 1 ({name}): gamma_x, gamma_y, gamma_z = {self.gamma[0]:g}, {self.gamma[1]:g}, "
                    f"{self.gamma[2]:g}; mean Vs = {self.alpha:g}; plain horizontal distances")
        txt = f"model {self.model} ({name}); directionality alpha = {self.alpha:g}, Line D at {self.angle_deg:g} deg"
        if self.model in (2, 3, 4, 5, 6):
            txt += f"; distances x {self.length_to_m:g} to metres"
        return txt


def length_to_metres(gravity: float) -> float:
    """Model length unit in metres from the unit system of the gravity (D-CNV-08: g > 20 -> feet)."""
    return FT_TO_M if float(gravity) > 20.0 else 1.0
