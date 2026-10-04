"""SOIL-NON: nonlinear time-domain 1D site response (numerical core of the SOIL module's
nonlinear option, ``NLSOIL,1``).

Normative sources: requirements sections 1.2 (SOIL-NON) and 3.4.O (NLSOIL / NLSLAYER / DELNLS);
spec 05a section 6.6 (SOIL dialog "Nonlinear Soil" group), spec 02 section 3.2 (SOIL-NON warnings),
spec 11 section 2.9 and open questions OQ-18 / OQ-19; ACS SASSI V3 manual sections 9.17.19-9.17.20
(the hyperbolic stress-strain law, transcribed from the equation image).  The manual cites DEEPSOIL
(Hashash & Park 2001; Hashash et al. 2009, DEEPSOIL V3.7 user manual) as the theory of SOIL-NON;
the details the manual leaves open are taken from that literature and listed below as SN-1 ... SN-12.

Physics in one paragraph
------------------------
The soil column above the half-space is a chain of shear elements with lumped masses (a
"multi-degree-of-freedom lumped-mass shear beam", Hashash & Park 2001).  Each element obeys the
modified hyperbolic (MKZ) backbone of the manual (Matasovic & Vucetic 1993)::

    tau = G0 gamma / (1 + beta (|gamma| / gamma_r)^s)  [+ eta d(gamma)/dt]          manual 9.17.19

with the extended Masing rules for unloading and reloading (Masing 1926; Pyke 1979; Vucetic 1990;
Kramer 1996 section 6.4.3): after a reversal at (gamma_R, tau_R) the curve is
``tau = tau_R + 2 F((gamma - gamma_R)/2)`` (F = backbone); a branch that reaches the previous
reversal continues on the curve of the previous cycle; a branch that exceeds the largest past strain
continues on the backbone.  The small-strain (viscous) damping, which the hyperbolic model does not
give (zero hysteretic damping at small strains), is added as a viscous damping matrix (NLDampType).
The base is either rigid (BedInt 0: the base moves with the input motion, a *within* motion) or an
elastic half-space represented by the Joyner & Chen (1975) dashpot ``c = rho_r V_r`` per unit area
(BedInt 1, the input is the *outcrop* motion).  The equation of motion, written relative to the
input frame (the base for BedInt 0, the outcrop motion for BedInt 1)::

    M u'' + C u' + R(u) + c_b u'_base = - M 1 a_in(t)

is integrated with the implicit Newmark average-acceleration method (gamma = 1/2, beta = 1/4) and
Newton-Raphson equilibrium iterations with the consistent tangent of the Masing branches.

Implementation decisions (SASSI-EDU interpretation of OQ-18 / OQ-19, with sources)
---------------------------------------------------------------------------------
SN-1  ``NLSLAYER <Num>`` is the SOIL sublayer number (SPRO numbering, 1 = top); every soil sublayer
      needs a set (Error 125); a set for the half-space (last SPRO entry) is not used.
SN-2  ``refStrain`` is in percent, like the DYNP strains and DEEPSOIL's "reference strain (%)"
      (manual: "typically around 0.03 %").  ``B`` = 0 with ``curvefit`` = 0 makes the sublayer
      linear elastic (the backbone degenerates to tau = G0 gamma).
SN-3  Curve fit (``curvefit`` = 1): least squares on the sublayer's DYNP G/Gmax points (equal weights).
      beta and gamma_r enter the backbone only through beta / gamma_r^s, so they cannot both be fitted:
      the fit sets beta = 1 and returns gamma_r (G/Gmax = 1/2 at gamma_r) and s, with s restricted to
      [0.2, 1] so that the backbone stress never decreases (s > 1 gives a peak stress and softening).
      The manual: with curve fit B, S, refStrain and Vis are ignored.
SN-4  Small-strain damping ratio xi_min of a sublayer: the DYNP damping at the smallest strain of its
      curve (D_min, DEEPSOIL's "small strain damping") when the sublayer has a DYNP label, otherwise
      the S-wave damping of its L property.
SN-5  NLDampType 1 "frequency independent" (DEEPSOIL builds it from the eigenvalues and eigenvectors of
      the system without user-selected modes or frequencies, Phillips & Hashash 2009; DEEPSOIL V3.7 manual
      2.10.1).  SASSI-EDU uses the classical modal damping matrix
      ``C = M Phi diag(2 xi_n w_n) Phi^T M`` of the small-strain column on a fixed base (all modes,
      mass-normalised; Clough & Penzien 1993; Chopra 2012 section 11.5), with xi_n the strain-energy
      weighted average of the element ratios xi_min (composite modal damping, Roesset, Whitman &
      Dobry 1973).  Every fixed-base mode then has its damping ratio at its own frequency, which is
      the response of SHAKE's frequency-independent complex modulus at resonance.  For BedInt 1 the
      matrix acts on the deformation relative to the base node (C 1 = 0).
      NLDampType 2 "visco-elastic": element dashpots eta / h (the manual's eta d(gamma)/dt term);
      eta = NLSLAYER Vis, or, when Vis is ignored (curve fit), eta = 2 xi_min G0 / w_1 (xi_min matched at
      the fundamental frequency w_1 of the small-strain column).
      NLDampType 3 "Rayleigh": C = MMmult M + SMmult K0 (K0 small-strain stiffness, DEEPSOIL default
      "update K matrix: no"); with both multipliers 0 they are set to give the thickness-weighted
      average xi_min at f_1 and 5 f_1 (DEEPSOIL's recommended frequencies).
      NLDampType 0 (dialog default, not a documented choice) is treated as 1.
SN-6  The stiffness-proportional, visco-elastic and modal damping forces depend only on the deformation
      of the column (C 1 = 0).  The mass-proportional Rayleigh term acts on the velocity relative to the
      input frame (rigid base: the usual relative formulation; elastic base: the outcrop motion), a
      documented approximation for BedInt 1.
SN-7  BedInt 1 (viscoelastic bedrock, "currently disabled" in ACS SASSI V3) is implemented with the
      Joyner & Chen (1975) dashpot: in the outcrop frame the base dashpot acts on the base node's
      relative velocity and the half-space material damping is not represented.
SN-8  Input: band-limited (FFT) interpolation of the control motion for the sub-increments (the
      Fourier representation SOIL-EQL uses), after the optional EDUOPT,SOILCUTOFF low-pass; the
      response is computed over NFFT samples (zero padding), like SOIL-EQL.
SN-9  Sub-increments: ``NSTimeSunInc`` = n >= 1 gives n fixed sub-increments per input time step
      (DEEPSOIL "fixed step"); 0 gives the flexible step: base sub-steps of at most 1/400 s (16 per
      period at 25 Hz, so the Newmark period elongation stays below 1.3 % up to 25 Hz), and a sub-step
      whose largest strain increment exceeds 0.005 % is repeated with a smaller step (DEEPSOIL's
      flexible step, Hashash & Park 2001); the step grows again where the strain rate falls.
SN-10 Convergence (per sub-step, both required): force ``||r||_2 <= ForceConv * F_ref`` with F_ref the
      largest of the norms of the inertia, internal, damping and input forces; displacement
      ``||du||_2 <= DispConv * max(||u_n+1 - u_n||_2, 1e-6 ||u_n+1||_2)`` with du the correction the
      current residual still asks for and the displacements taken relative to the column base (the
      rigid-body motion of BedInt 1 does not enter).  Defaults (0 = default): DispConv = ForceConv =
      1e-6, EqualIt = 25.  The Newton unknown is the step increment.  A sub-step that does not
      converge within EqualIt iterations is bisected (up to 2^12 sub-steps per input step) before the
      run stops (manual: "before the program terminates").
SN-11 Element size: each SPRO sublayer is divided into an odd number of equal elements of thickness
      h <= Vs / (10 f_max), f_max = min(25 Hz, Nyquist) (10 elements per wavelength at DEEPSOIL's
      minimum recommended maximum frequency); the middle element gives the mid-height strain and
      stress.  The DEEPSOIL frequency check f = Vs / (4 h) is printed.
SN-12 Outputs: ACCxxx.TH absolute acceleration at sublayer tops (within; the outcrop motion exists
      only at the base for BedInt 1, where it is the input), SNxxx.TH strain (%) and SSxxx.TH total
      shear stress (hysteretic + viscous, from equilibrium: the inertia of the soil above) at mid-height;
      FILE88 equivalent-linear properties at the effective strain ratio x max|gamma|: secant G of the
      backbone and damping xi_min + Masing damping (the properties of the model that was integrated).

References
----------
Hashash, Y.M.A., Park, D. (2001) Non-linear one-dimensional seismic ground motion propagation in the
Mississippi embayment. Engineering Geology 62, 185-206.
Hashash, Y.M.A., Groholski, D.R., Phillips, C.A., Park, D. (2009) DEEPSOIL V3.7beta, User Manual and
Tutorial (sections 2.5.2, 2.6, 2.7, 2.8.1, 2.10).
Phillips, C., Hashash, Y.M.A. (2009) Damping formulation for nonlinear 1D site response analyses.
Soil Dynamics and Earthquake Engineering 29(7), 1143-1158.
Matasovic, N., Vucetic, M. (1993) Cyclic characterization of liquefiable sands. J. Geotech. Eng.
119(11), 1805-1822.
Joyner, W.B., Chen, A.T.F. (1975) Calculation of nonlinear ground response in earthquakes. BSSA 65(5),
1315-1336.
Kramer, S.L. (1996) Geotechnical Earthquake Engineering, section 6.4.3 (Masing / extended Masing rules).
Clough, R.W., Penzien, J. (1993) Dynamics of Structures; Chopra, A.K. (2012) Dynamics of Structures,
section 11.5 (classical damping matrices).
"""
from __future__ import annotations

import math
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy import linalg as _sla
from scipy import optimize as _opt
from scipy import special as _sp

PathLike = Union[str, Path]

# Newmark average acceleration (unconditionally stable, no numerical damping)
NEWMARK_GAMMA = 0.5
NEWMARK_BETA = 0.25
#: SN-10 defaults (0 in NLSOIL means "default")
DEFAULT_DISPCONV = 1.0e-6
DEFAULT_FORCECONV = 1.0e-6
DEFAULT_EQUALIT = 25
#: SN-9 flexible sub-incrementation: maximum strain increment (fraction; DEEPSOIL 0.005 %)
FLEX_STRAIN_INCREMENT = 0.005e-2
#: SN-9 flexible scheme: largest base sub-step (16 steps per period at 25 Hz: Newmark period
#: elongation (w dt)^2 / 12 < 1.3 % up to 25 Hz)
DT_MAX_FLEX = 1.0 / (16.0 * 25.0)
#: bisection depth of a sub-step (2**KMAX sub-steps per base increment at most)
KMAX = 12
#: SN-11 discretisation
F_MAX_DISCRETISATION = 25.0
ELEMENTS_PER_WAVELENGTH = 10
#: SN-8 band-limited input: samples of the fine input grid per time step
INPUT_UPSAMPLING = 16
#: SN-3 curve-fit bounds of the exponent s
S_BOUNDS = (0.2, 1.0)

DAMPING_TYPES = {1: "frequency independent (modal, all fixed-base modes)",
                 2: "visco-elastic (element dashpots eta/h)",
                 3: "Rayleigh C = alpha M + beta K0"}
BEDROCK_TYPES = {0: "rigid (within motion at the base)",
                 1: "viscoelastic half-space, Joyner-Chen dashpot (outcrop motion)"}


class SoilNonError(RuntimeError):
    """Input or convergence error of a SOIL-NON analysis (the SOIL module stops with the message)."""


# =======================================================================================
# NLSOIL / NLSLAYER data and the <model>.nls side file
# =======================================================================================
@dataclass
class NLSoilOptions:
    """``NLSOIL,<Opt>,<NSTimeSunInc>,<DispConv>,<ForceConv>,<EqualIt>,<BedInt>,<NLDampType>,<MMmult>,<SMmult>``
    (manual 9.17.20; every argument defaults to 0)."""

    opt: int = 0
    nsub: int = 0
    dispconv: float = 0.0
    forceconv: float = 0.0
    equalit: int = 0
    bedint: int = 0
    damptype: int = 0
    mmmult: float = 0.0
    smmult: float = 0.0

    FIELDS = ("opt", "nsub", "dispconv", "forceconv", "equalit", "bedint", "damptype", "mmmult", "smmult")
    INTS = ("opt", "nsub", "equalit", "bedint", "damptype")

    @classmethod
    def from_values(cls, values: Sequence[Optional[float]]) -> "NLSoilOptions":
        kw = {}
        for name, v in zip(cls.FIELDS, values):
            if v is None:
                continue
            kw[name] = int(round(float(v))) if name in cls.INTS else float(v)
        return cls(**kw)

    def values(self) -> List[float]:
        return [getattr(self, n) for n in self.FIELDS]

    def problems(self) -> List[str]:
        """Values the analysis cannot use (the SOIL module stops with these messages)."""
        out = []
        if self.opt not in (0, 1):
            out.append(f"NLSOIL <Opt> = {self.opt} (0 or 1)")
        if self.nsub < 0:
            out.append(f"NLSOIL <NSTimeSunInc> = {self.nsub} (>= 0; 0 = flexible sub-incrementation)")
        if self.dispconv < 0 or self.forceconv < 0:
            out.append("NLSOIL <DispConv> and <ForceConv> must be >= 0 (0 = default 1e-6)")
        if self.equalit < 0:
            out.append(f"NLSOIL <EqualIt> = {self.equalit} (>= 0; 0 = default {DEFAULT_EQUALIT})")
        if self.bedint not in (0, 1):
            out.append(f"NLSOIL <BedInt> = {self.bedint} (0 rigid, 1 viscoelastic)")
        if self.damptype not in (0, 1, 2, 3):
            out.append(f"NLSOIL <NLDampType> = {self.damptype} (1, 2 or 3)")
        if self.mmmult < 0 or self.smmult < 0:
            out.append("NLSOIL Rayleigh multipliers <MMmult>, <SMmult> must be >= 0")
        return out


@dataclass
class NLLayerData:
    """``NLSLAYER,<Num>,[curvefit],[B],[S],[refStrain],[Vis]`` (manual 9.17.19; refStrain in %, SN-2)."""

    num: int
    curvefit: int = 0
    beta: float = 0.0
    s: float = 0.0
    refstrain: float = 0.0
    vis: float = 0.0

    FIELDS = ("num", "curvefit", "beta", "s", "refstrain", "vis")

    @classmethod
    def from_values(cls, values: Sequence[Optional[float]]) -> "NLLayerData":
        vals = list(values) + [None] * 6
        if vals[0] is None:
            raise ValueError("NLSLAYER needs the layer number <Num>")
        kw = {"num": int(round(float(vals[0])))}
        for name, v in zip(cls.FIELDS[1:], vals[1:6]):
            if v is not None:
                kw[name] = int(round(float(v))) if name == "curvefit" else float(v)
        return cls(**kw)

    def values(self) -> List[float]:
        return [getattr(self, n) for n in self.FIELDS]

    def problems(self) -> List[str]:
        out = []
        if self.num < 1:
            out.append(f"NLSLAYER <Num> = {self.num} (sublayer numbers start at 1)")
        if self.curvefit not in (0, 1):
            out.append(f"NLSLAYER {self.num}: <curvefit> = {self.curvefit} (0 or 1)")
        if self.curvefit == 0 and self.beta < 0:
            out.append(f"NLSLAYER {self.num}: Beta = {self.beta:g} < 0")
        if self.curvefit == 0 and self.beta > 0 and (self.s <= 0 or self.refstrain <= 0):
            out.append(f"NLSLAYER {self.num}: S exponent and Reference Strain must be > 0 when Beta > 0 "
                       f"(curvefit 0: the user must define all curve parameters)")
        if self.vis < 0:
            out.append(f"NLSLAYER {self.num}: Viscosity = {self.vis:g} < 0")
        return out


@dataclass
class NlsData:
    """Content of the ``<model>.nls`` side file (NLSOIL record + NLSLAYER sets)."""

    options: NLSoilOptions
    layers: Dict[int, NLLayerData] = field(default_factory=dict)
    model: str = ""
    model_hash: str = ""


NLS_HEADER = "# SASSI-EDU NLS v1  SOIL-NON options (NLSOIL / NLSLAYER) for the SOIL module"


def format_value(v: float) -> str:
    """Number as written in the .nls file and the listings (integers without a decimal point)."""
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return f"{float(v):.10g}"


def nls_text(options: NLSoilOptions, layers: Sequence[NLLayerData], model: str = "", model_hash: str = "") -> str:
    """Text of the ``<model>.nls`` side file written by AFWRITE next to ``<model>.soi``.

    The SOIL deck schema has no SOIL-NON fields, so AFWRITE carries the NLSOIL record and the NLSLAYER
    sets in this file (one command line each, the manual's syntax, every argument written)::

        # SASSI-EDU NLS v1 ...
        # model: <name>
        # model_hash = <hash of the model that wrote the deck>
        NLSOIL,1,0,0,0,0,1,1,0,0
        NLSLAYER,1,1,0,0,0,0
    """
    lines = [NLS_HEADER]
    if model:
        lines.append(f"# model: {model}")
    if model_hash:
        lines.append(f"# model_hash = {model_hash}")
    lines.append("NLSOIL," + ",".join(format_value(v) for v in options.values()))
    for lay in sorted(layers, key=lambda x: x.num):
        lines.append("NLSLAYER," + ",".join(format_value(v) for v in lay.values()))
    return "\n".join(lines) + "\n"


def _num(tok: str) -> Optional[float]:
    t = tok.strip()
    if not t:
        return None
    return float(t.replace("D", "E").replace("d", "e"))


def parse_nls(text: str) -> NlsData:
    """Parse the text of a ``.nls`` file (see :func:`nls_text`)."""
    opts: Optional[NLSoilOptions] = None
    layers: Dict[int, NLLayerData] = {}
    model = mhash = ""
    for k, raw in enumerate(text.splitlines(), start=1):
        s = raw.strip()
        if not s:
            continue
        if s.startswith("#") or s.startswith("*"):
            m = re.match(r"[#*]\s*model:\s*(.*)$", s)
            if m:
                model = m.group(1).strip()
            m = re.match(r"[#*]\s*model_hash\s*=\s*(\S+)", s)
            if m:
                mhash = m.group(1)
            continue
        toks = s.split(",")
        name = toks[0].strip().upper()
        try:
            vals = [_num(t) for t in toks[1:]]
        except ValueError:
            raise ValueError(f".nls line {k}: non-numeric argument in '{s}'") from None
        if name == "NLSOIL":
            opts = NLSoilOptions.from_values(vals)
        elif name == "NLSLAYER":
            lay = NLLayerData.from_values(vals)
            layers[lay.num] = lay
        else:
            raise ValueError(f".nls line {k}: unknown record '{toks[0].strip()}' (NLSOIL or NLSLAYER expected)")
    if opts is None:
        raise ValueError(".nls file has no NLSOIL record")
    return NlsData(opts, layers, model, mhash)


def read_nls(path: PathLike) -> NlsData:
    return parse_nls(Path(path).read_text(encoding="utf-8"))


def write_nls(path: PathLike, data: NlsData) -> Path:
    path = Path(path)
    path.write_text(nls_text(data.options, list(data.layers.values()), data.model, data.model_hash), encoding="utf-8")
    return path


# =======================================================================================
# MKZ backbone, Masing damping, curve fit
# =======================================================================================
def mkz_stress(g, G0, gr, beta, s):
    """Backbone ``F(gamma) = G0 gamma / (1 + beta (|gamma|/gamma_r)^s)`` (strains as fractions)."""
    g = np.asarray(g, dtype=float)
    return G0 * g / (1.0 + beta * (np.abs(g) / gr) ** s)


def mkz_tangent(g, G0, gr, beta, s):
    """``dF/dgamma = G0 (1 + beta (1 - s) x^s) / (1 + beta x^s)^2``, x = |gamma| / gamma_r."""
    z = (np.abs(np.asarray(g, dtype=float)) / gr) ** s
    return G0 * (1.0 + beta * (1.0 - s) * z) / (1.0 + beta * z) ** 2


def mkz_modulus_ratio(strain_pct, gr_pct, beta, s):
    """Secant ratio ``G/G0 = 1 / (1 + beta (gamma/gamma_r)^s)`` (strain and gamma_r in the same units)."""
    x = np.abs(np.asarray(strain_pct, dtype=float))
    if beta == 0:
        return np.ones_like(x) if np.ndim(x) else 1.0
    out = 1.0 / (1.0 + beta * (x / gr_pct) ** s)
    return out if np.ndim(out) else float(out)


def masing_damping(strain, gr, beta, s):
    """Hysteretic damping ratio (fraction) of a Masing loop of amplitude ``strain`` on the MKZ backbone.

    The Masing loop of amplitude gamma_a has the secant modulus of the backbone and the area
    ``A = 8 int_0^ga F dg - 4 ga F(ga)``, so ``D = A / (4 pi W)`` with ``W = ga F(ga) / 2``::

        D = (2/pi) [ 2 int_0^ga F(g) dg / (ga F(ga)) - 1 ]

    For the MKZ backbone ``int_0^ga F dg = G0 ga^2/2 * 2F1(1, 2/s; 1 + 2/s; -z)`` with
    ``z = beta (ga/gamma_r)^s`` (Gauss hypergeometric function).  For s = 1 this is the closed form of
    the hyperbolic model used by Darendeli (2001) (R2 H.2, "D_Masing,a=1").  ``strain`` and ``gr`` in
    the same units.
    """
    g = np.abs(np.asarray(strain, dtype=float))
    if beta == 0:
        return np.zeros_like(g) if np.ndim(g) else 0.0
    z = beta * (g / gr) ** s
    with np.errstate(divide="ignore", invalid="ignore"):
        if abs(s - 1.0) < 1e-14:
            # int/(g F) = (1+z)(z - ln(1+z))/z^2 ; series for small z
            lz = np.log1p(z)
            ratio = np.where(z > 1e-4, (1.0 + z) * (z - lz) / np.where(z > 0, z * z, 1.0),
                             0.5 + z / 6.0 - z * z / 12.0)
        else:
            h = _sp.hyp2f1(1.0, 2.0 / s, 1.0 + 2.0 / s, -z)
            ratio = 0.5 * h * (1.0 + z)
    d = (2.0 / np.pi) * (2.0 * ratio - 1.0)
    d = np.where(g > 0, d, 0.0)
    return d if np.ndim(d) else float(d)


@dataclass
class MKZFit:
    """Result of the curve fit of the MKZ backbone to a G/Gmax curve (SN-3)."""

    beta: float
    s: float
    gr_pct: float
    rms: float
    max_err: float
    npts: int
    note: str = ""


def fit_mkz(strain_pct: Sequence[float], g_ratio: Sequence[float],
            s_bounds: Tuple[float, float] = S_BOUNDS) -> MKZFit:
    """Fit ``G/Gmax = 1/(1 + (gamma/gamma_r)^s)`` (beta = 1, SN-3) to curve points by least squares.

    ``strain_pct`` in percent (DYNP units); points with zero strain are ignored.  Returns beta = 1,
    s and gamma_r (percent) with the rms and largest absolute error of G/Gmax at the points.  A curve
    without reduction (all G/Gmax >= 1) gives the linear model beta = 0.
    """
    x = np.asarray(strain_pct, dtype=float)
    y = np.asarray(g_ratio, dtype=float)
    keep = x > 0
    x, y = x[keep], y[keep]
    if x.size == 0:
        raise ValueError("curve fit: no G/Gmax point with a positive strain")
    order = np.argsort(x)
    x, y = x[order], y[order]
    if np.all(y >= 1.0 - 1e-12):
        return MKZFit(0.0, 1.0, float("inf"), 0.0, 0.0, int(x.size), "no modulus reduction: linear")
    if x.size < 2:
        raise ValueError("curve fit: at least two G/Gmax points with a positive strain are needed")
    s0 = 0.9
    below = np.nonzero(y < 0.5)[0]
    if below.size and below[0] > 0:
        j = below[0]
        lx = np.log10(x[j - 1]) + (0.5 - y[j - 1]) * (np.log10(x[j]) - np.log10(x[j - 1])) / (y[j] - y[j - 1])
        gr0 = 10.0 ** lx
    else:
        j = int(np.argmin(y))
        yy = min(max(y[j], 1e-3), 1.0 - 1e-6)
        gr0 = x[j] * (1.0 / yy - 1.0) ** (-1.0 / s0)
    lo = (math.log10(x[0]) - 4.0, s_bounds[0])
    hi = (math.log10(x[-1]) + 4.0, s_bounds[1])
    p0 = np.array([min(max(math.log10(gr0), lo[0] + 1e-6), hi[0] - 1e-6), min(max(s0, lo[1]), hi[1])])
    lx = np.log10(x)

    def resid(p):
        return 1.0 / (1.0 + 10.0 ** (p[1] * (lx - p[0]))) - y

    sol = _opt.least_squares(resid, p0, bounds=(lo, hi), method="trf", xtol=1e-14, ftol=1e-14, gtol=1e-14,
                             max_nfev=2000)
    e = resid(sol.x)
    return MKZFit(1.0, float(sol.x[1]), float(10.0 ** sol.x[0]), float(np.sqrt(np.mean(e * e))),
                  float(np.max(np.abs(e))), int(x.size))


# =======================================================================================
# Extended Masing hysteresis (vectorised over elements)
# =======================================================================================
@dataclass
class MasingTrial:
    """Trial state of :class:`MasingMKZ` at strains ``g`` (not yet committed)."""

    g: np.ndarray
    tau: np.ndarray
    kt: np.ndarray
    top: np.ndarray
    rev: np.ndarray
    d: np.ndarray
    og: Optional[np.ndarray] = None      # origin (reversal) of the current branch
    ot: Optional[np.ndarray] = None
    tg: Optional[np.ndarray] = None      # strain at which the current branch closes (where top >= 1)


class MasingMKZ:
    """MKZ backbone with the extended Masing rules, one material point per element.

    State per element: committed strain/stress, last loading direction and a stack of reversal points
    ``(gamma_R, tau_R)``.  ``top`` = 0 means "on the backbone"; otherwise the current branch starts at
    the last reversal R_top-1 with ``tau = tau_R + 2 F((gamma - gamma_R)/2)``.  Rules (Kramer 1996,
    section 6.4.3):

    1. first loading follows the backbone;
    2. a strain reversal starts a Masing branch (the backbone scaled by 2 about the reversal point);
    3. a branch that passes the previous reversal R_top-2 closes its loop: both reversals are removed
       and the curve of the previous cycle is followed (rule 4 of the extended rules);
    4. the first branch from the backbone (one reversal R_0) rejoins the backbone at -gamma_R0 (the
       largest past strain on the other side), after which the backbone is followed (rule 3).

    A step from the committed strain to a trial strain is a straight path, so at most one reversal is
    added per step (at the committed point) and any number of loops can close.  Linear elements
    (beta = 0) keep no reversal stack.
    """

    def __init__(self, G0, gr, beta, s, cap: int = 64):
        self.G0 = np.array(G0, dtype=float, ndmin=1)
        n = self.G0.size
        gr = np.broadcast_to(np.asarray(gr, dtype=float), (n,)).astype(float)
        beta = np.broadcast_to(np.asarray(beta, dtype=float), (n,)).astype(float)
        s = np.broadcast_to(np.asarray(s, dtype=float), (n,)).astype(float)
        self.linear = (beta <= 0) | ~np.isfinite(gr) | ~np.isfinite(beta)
        if np.any(~self.linear & ((gr <= 0) | (s <= 0))):
            raise ValueError("MKZ: gamma_r and s must be > 0 for nonlinear elements")
        self.gr = np.where(self.linear, 1.0, gr)
        self.beta = np.where(self.linear, 0.0, beta)
        self.s = np.where(self.linear, 1.0, s)
        self.n = n
        self._ar = np.arange(n)
        self.g = np.zeros(n)
        self.t = np.zeros(n)
        self.d = np.zeros(n, dtype=np.int8)
        self.top = np.zeros(n, dtype=np.int64)
        self.rg = np.zeros((n, cap))
        self.rt = np.zeros((n, cap))
        # current branch (cached from the stack): origin R_top-1 and closing strain R_top-2 (or -R_0)
        self.og = np.zeros(n)
        self.ot = np.zeros(n)
        self.tg = np.zeros(n)              # meaningful only where top >= 1
        self.any_nonlinear = bool(np.any(~self.linear))

    # -- backbone ----------------------------------------------------------------------
    def F(self, x: np.ndarray) -> np.ndarray:
        return self.G0 * x / (1.0 + self.beta * (np.abs(x) / self.gr) ** self.s)

    def dF(self, x: np.ndarray) -> np.ndarray:
        z = (np.abs(x) / self.gr) ** self.s
        return self.G0 * (1.0 + self.beta * (1.0 - self.s) * z) / (1.0 + self.beta * z) ** 2

    # -- trial / commit -------------------------------------------------------------------
    def trial(self, g: np.ndarray) -> MasingTrial:
        """Stress and tangent at strains ``g`` reached from the committed state along a straight path."""
        g = np.asarray(g, dtype=float)
        dg = g - self.g
        d = np.sign(dg).astype(np.int8)
        moving = d != 0
        if not self.any_nonlinear:
            tau = self.G0 * g
            return MasingTrial(g.copy(), tau, self.G0.copy(), self.top.copy(), np.zeros(self.n, bool), d)
        rev = moving & (self.d != 0) & (d != self.d) & ~self.linear
        # a reversal at the committed point starts a new branch: origin = committed point, closing strain =
        # the origin of the previous branch (or the mirror point -gamma_R0 of the first branch)
        top = self.top + rev.astype(np.int64)
        og = np.where(rev, self.g, self.og)
        ot = np.where(rev, self.t, self.ot)
        tg = np.where(rev, np.where(self.top >= 1, self.og, -self.g), self.tg)
        passed = moving & (top >= 1) & (d * (g - tg) > 0)
        if passed.any():
            # loop closures (rules 3-4): pop reversals until the branch reaches beyond the trial strain
            top, og, ot, tg = top.copy(), og.copy(), ot.copy(), tg.copy()
            while passed.any():
                i = np.nonzero(passed)[0]
                ti = top[i]
                ti = np.where(ti >= 2, ti - 2, ti - 1)
                top[i] = ti
                kk = np.maximum(ti - 1, 0)
                vg, vt = self._rev_at_idx(i, kk, rev)
                og[i] = vg
                ot[i] = vt
                k2 = np.maximum(ti - 2, 0)
                tgi, _ = self._rev_at_idx(i, np.where(ti >= 2, k2, 0), rev)
                tg[i] = np.where(ti >= 2, tgi, np.where(ti == 1, -tgi, 0.0))
                passed = np.zeros(self.n, bool)
                passed[i] = (ti >= 1) & (d[i] * (g[i] - tg[i]) > 0)
        bb = top == 0
        x = np.where(bb, g, 0.5 * (g - og))
        Fx = self.F(x)
        tau = np.where(bb, Fx, ot + 2.0 * Fx)
        tau = np.where(moving, tau, self.t)
        kt = self.dF(x)
        return MasingTrial(g.copy(), tau, kt, top, rev, d, og, ot, tg)

    def _rev_at_idx(self, i: np.ndarray, k: np.ndarray, rev: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Reversal ``k`` of elements ``i`` (the point pushed by a reversal of this trial included)."""
        cap = self.rg.shape[1]
        kk = np.clip(k, 0, cap - 1)
        vg = self.rg[i, kk]
        vt = self.rt[i, kk]
        ov = rev[i] & (k == self.top[i])
        if ov.any():
            vg = np.where(ov, self.g[i], vg)
            vt = np.where(ov, self.t[i], vt)
        return vg, vt

    def commit(self, tr: MasingTrial) -> None:
        if tr.rev.any():
            idx = np.nonzero(tr.rev)[0]
            pos = self.top[idx]
            cap = self.rg.shape[1]
            if int(pos.max()) >= cap:
                new = max(2 * cap, int(pos.max()) + 1)
                if new > 1 << 20:
                    raise SoilNonError("Masing reversal memory exceeded")
                self.rg = np.concatenate([self.rg, np.zeros((self.n, new - cap))], axis=1)
                self.rt = np.concatenate([self.rt, np.zeros((self.n, new - cap))], axis=1)
            self.rg[idx, pos] = self.g[idx]
            self.rt[idx, pos] = self.t[idx]
        self.top = tr.top.copy()
        if tr.og is not None:
            self.og, self.ot, self.tg = tr.og, tr.ot, tr.tg
        self.g = tr.g.copy()
        self.t = tr.tau.copy()
        self.d = np.where(tr.d != 0, tr.d, self.d).astype(np.int8)

    def drive(self, strains: Sequence[float]) -> np.ndarray:
        """Stress history of a strain path (each value committed; used by tests and VP-SN2)."""
        out = []
        for gv in strains:
            tr = self.trial(np.full(self.n, float(gv)) if np.ndim(gv) == 0 else np.asarray(gv, float))
            self.commit(tr)
            out.append(tr.tau.copy())
        return np.asarray(out)


# =======================================================================================
# Column: discretisation, mass, stiffness
# =======================================================================================
@dataclass
class SublayerModel:
    """Nonlinear model of one soil sublayer (resolved from NLSLAYER and DYNP, SN-2 ... SN-4)."""

    beta: float          # 0 = linear
    s: float
    gr_pct: float        # reference strain, percent
    xi_min: float        # small-strain viscous damping ratio
    eta: float = 0.0     # viscosity (NLDampType 2), stress x time
    source: str = ""     # 'user', 'fit', 'linear'
    fit: Optional[MKZFit] = None
    xi_source: str = ""

    @property
    def linear(self) -> bool:
        return self.beta <= 0 or not math.isfinite(self.gr_pct)


@dataclass
class SoilNonColumn:
    """Lumped-mass shear column of SOIL-NON (soil sublayers top first; the half-space is the base).

    Sublayer ``j`` (0-based) is divided into ``nel[j]`` equal elements (odd, SN-11).  Node 0 is the
    ground surface, node ``ne`` the top of the half-space.  Masses: half of each adjacent element
    (``rho h / 2``); the half-space adds no mass (it is the Joyner-Chen dashpot for BedInt 1).
    """

    thick: np.ndarray            # (ns,) sublayer thickness
    rho: np.ndarray              # (ns,) mass density
    vs: np.ndarray               # (ns,) small-strain shear-wave velocity
    models: List[SublayerModel]
    base_rho: float
    base_vs: float
    nel: np.ndarray              # (ns,) elements per sublayer

    def __post_init__(self):
        self.thick = np.asarray(self.thick, dtype=float)
        self.rho = np.asarray(self.rho, dtype=float)
        self.vs = np.asarray(self.vs, dtype=float)
        self.nel = np.asarray(self.nel, dtype=int)
        ns = len(self.thick)
        lay = np.repeat(np.arange(ns), self.nel)
        self.elem_layer = lay
        self.h = (self.thick / self.nel)[lay]
        self.elem_rho = self.rho[lay]
        self.G0 = (self.rho * self.vs ** 2)[lay]
        self.ne = int(lay.size)
        self.nn = self.ne + 1
        m = np.zeros(self.nn)
        half = 0.5 * self.elem_rho * self.h
        m[:-1] += half
        m[1:] += half
        self.mass = m
        first = np.concatenate([[0], np.cumsum(self.nel)])
        self.top_node = first                         # node at the top of sublayer j (j = ns: base)
        self.mid_elem = first[:-1] + self.nel // 2    # element at mid-height of sublayer j
        self.elem_beta = np.array([self.models[j].beta for j in lay])
        self.elem_s = np.array([self.models[j].s for j in lay])
        self.elem_gr = np.array([self.models[j].gr_pct for j in lay]) / 100.0
        self.elem_xi = np.array([self.models[j].xi_min for j in lay])
        self.elem_eta = np.array([self.models[j].eta for j in lay])

    @property
    def ns(self) -> int:
        return len(self.thick)

    @property
    def k0(self) -> np.ndarray:
        """Small-strain element stiffnesses G0 / h (per unit area)."""
        return self.G0 / self.h

    @property
    def depth_top(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.thick)])

    def material(self) -> MasingMKZ:
        return MasingMKZ(self.G0, self.elem_gr, self.elem_beta, self.elem_s)

    def deepsoil_fmax(self) -> np.ndarray:
        """DEEPSOIL's maximum frequency of a sublayer's elements ``Vs / (4 h)`` (Hz)."""
        return self.vs / (4.0 * self.thick / self.nel)


def elements_per_sublayer(thick: Sequence[float], vs: Sequence[float], fmax: float,
                          per_wavelength: int = ELEMENTS_PER_WAVELENGTH) -> np.ndarray:
    """SN-11: smallest odd number of equal elements with ``h <= Vs / (per_wavelength * fmax)``."""
    thick = np.asarray(thick, dtype=float)
    vs = np.asarray(vs, dtype=float)
    hmax = vs / (per_wavelength * fmax)
    n = np.maximum(1, np.ceil(thick / hmax - 1e-9).astype(int))
    n = np.where(n % 2 == 0, n + 1, n)
    return n


def fixed_base_modes(col: SoilNonColumn, stiffness: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Natural circular frequencies and mass-normalised modes of the column on a fixed base.

    Returns ``(w, Phi)`` with ``Phi`` of shape (ne, ne) (nodes 0..ne-1, base fixed), ``Phi^T M Phi = I``.
    The generalised problem is reduced to the symmetric tridiagonal ``M^-1/2 K M^-1/2``.
    """
    k = col.k0 if stiffness is None else np.asarray(stiffness, dtype=float)
    ne = col.ne
    m = col.mass[:ne]
    kd = np.zeros(ne)
    kd += k                      # element i below node i (node ne is the fixed base)
    kd[1:] += k[:-1]             # element i-1 above node i
    sm = np.sqrt(m)
    a = kd / m
    b = -k[:-1] / (sm[:-1] * sm[1:])
    lam, y = _sla.eigh_tridiagonal(a, b)
    lam = np.maximum(lam, 0.0)
    phi = y / sm[:, None]
    return np.sqrt(lam), phi


# =======================================================================================
# Viscous damping (SN-5, SN-6)
# =======================================================================================
@dataclass
class DampingModel:
    """Viscous damping matrix on the free DOFs: tridiagonal (``diag``, ``off``) or dense (``full``)."""

    kind: str                        # 'tri' or 'dense'
    diag: Optional[np.ndarray] = None
    off: Optional[np.ndarray] = None
    full: Optional[np.ndarray] = None
    type: int = 1
    alpha: float = 0.0
    beta: float = 0.0
    modal_w: Optional[np.ndarray] = None
    modal_xi: Optional[np.ndarray] = None
    eta: Optional[np.ndarray] = None
    notes: List[str] = field(default_factory=list)

    def matvec(self, v: np.ndarray) -> np.ndarray:
        if self.kind == "dense":
            with np.errstate(all="ignore"):
                return self.full @ v
        out = self.diag * v
        out[:-1] += self.off * v[1:]
        out[1:] += self.off * v[:-1]
        return out


def _tri_from_elements(c_el: np.ndarray, nfree: int) -> Tuple[np.ndarray, np.ndarray]:
    """Assemble element 'springs' c_el (between nodes e, e+1) into a tridiagonal matrix on nodes 0..nfree-1."""
    ne = c_el.size
    d = np.zeros(ne + 1)
    d[:-1] += c_el
    d[1:] += c_el
    off = -c_el.copy()
    return d[:nfree].copy(), off[:nfree - 1].copy()


def build_damping(col: SoilNonColumn, damptype: int, bedint: int, mmmult: float = 0.0,
                  smmult: float = 0.0) -> DampingModel:
    """Viscous damping matrix of NLDampType 1 / 2 / 3 on the free DOFs (SN-5, SN-6).

    Free DOFs: nodes 0..ne-1 (rigid base, BedInt 0) or 0..ne (BedInt 1).  The Joyner-Chen base dashpot
    is *not* included (it is added by the integrator).
    """
    nfree = col.ne + (1 if bedint == 1 else 0)
    t = 1 if damptype == 0 else int(damptype)
    notes: List[str] = []
    if damptype == 0:
        notes.append("NLDampType 0 is not a documented choice: type 1 (frequency independent) used")
    if t == 2:
        eta = col.elem_eta.copy()
        d, o = _tri_from_elements(eta / col.h, nfree)
        return DampingModel("tri", d, o, type=2, eta=eta, notes=notes)
    if t == 3:
        a, b = float(mmmult), float(smmult)
        if a == 0.0 and b == 0.0:
            w, _ = fixed_base_modes(col)
            w1 = float(w[0])
            w2 = 5.0 * w1
            xi = float(np.sum(col.elem_xi * col.h) / np.sum(col.h))
            a = 2.0 * xi * w1 * w2 / (w1 + w2)
            b = 2.0 * xi / (w1 + w2)
            notes.append(f"Rayleigh multipliers 0: alpha = {a:.6g}, beta = {b:.6g} give xi = {xi:.4f} at "
                         f"f1 = {w1 / (2 * np.pi):.3f} Hz and 5 f1 (DEEPSOIL recommended frequencies)")
        d, o = _tri_from_elements(b * col.k0, nfree)
        d = d + a * col.mass[:nfree]
        return DampingModel("tri", d, o, type=3, alpha=a, beta=b, notes=notes)
    # type 1: modal damping of the fixed-base column, strain-energy weighted ratios
    w, phi = fixed_base_modes(col)
    ne = col.ne
    dphi = np.vstack([phi[1:] - phi[:-1], -phi[-1:]])     # element strains x h of each mode (base fixed)
    ek = col.k0[:, None] * dphi ** 2                      # element strain energies (x2)
    tot = ek.sum(axis=0)
    xi_n = (col.elem_xi[:, None] * ek).sum(axis=0) / np.where(tot > 0, tot, 1.0)
    mphi = col.mass[:ne, None] * phi
    with np.errstate(all="ignore"):              # spurious FP flags of Accelerate's matmul (numpy 2.0)
        cf = (mphi * (2.0 * xi_n * w)) @ mphi.T
    cf = 0.5 * (cf + cf.T)
    if bedint == 1:
        full = np.zeros((ne + 1, ne + 1))
        full[:ne, :ne] = cf
        c1 = cf.sum(axis=1)
        full[:ne, ne] = -c1
        full[ne, :ne] = -c1
        full[ne, ne] = c1.sum()
    else:
        full = cf
    return DampingModel("dense", full=full, type=1, modal_w=w, modal_xi=xi_n, notes=notes)


# =======================================================================================
# Linear transfer function of the discrete model (checks and listing)
# =======================================================================================
def linear_transfer(col: SoilNonColumn, damping: DampingModel, bedint: int, freqs: Sequence[float],
                    nodes: Sequence[int] = (0,)) -> np.ndarray:
    """Steady-state absolute-acceleration transfer functions of the small-strain discrete model.

    ``H[f, j]`` = absolute acceleration of node ``nodes[j]`` per unit input acceleration (base motion
    for BedInt 0, outcrop motion for BedInt 1): ``(K0 + i w (C + c_b) - w^2 M) U = -M 1``,
    ``H = 1 - w^2 U`` (e^{+i w t}, D-CNV-01).  Exact for the discrete system, so the time integration of
    the linear model can be checked against it.
    """
    nfree = col.ne + (1 if bedint == 1 else 0)
    kd, ko = _tri_from_elements(col.k0, nfree)
    K = np.diag(kd) + np.diag(ko, 1) + np.diag(ko, -1)
    C = damping.full if damping.kind == "dense" else (np.diag(damping.diag) + np.diag(damping.off, 1)
                                                      + np.diag(damping.off, -1))
    C = C.copy()
    if bedint == 1:
        C[-1, -1] += col.base_rho * col.base_vs
    m = col.mass[:nfree]
    out = np.zeros((len(freqs), len(nodes)), dtype=complex)
    with np.errstate(all="ignore"):         # spurious FP flags of Accelerate's BLAS (numpy 2.0)
        for i, f in enumerate(np.asarray(freqs, dtype=float)):
            w = 2.0 * np.pi * f
            if w == 0.0:
                out[i, :] = 1.0             # static limit: the column moves with the input
                continue
            A = K + 1j * w * C - (w * w) * np.diag(m)
            U = np.linalg.solve(A, -m.astype(complex))
            for j, nd in enumerate(nodes):
                out[i, j] = 1.0 - w * w * U[nd] if nd < nfree else 1.0
    return out


# =======================================================================================
# Input motion (SN-8)
# =======================================================================================
def band_limited_input(acc: np.ndarray, dt: float, nfft: int, cutoff: float = 0.0,
                       up: int = INPUT_UPSAMPLING) -> Tuple[np.ndarray, np.ndarray]:
    """``(a_dt, a_fine)``: the control motion on the NFFT grid and its band-limited interpolation at
    ``dt / up`` (FFT zero padding; the Nyquist bin is split so ``a_fine[::up] == a_dt``).  Fourier
    components above ``cutoff`` (> 0) are removed first (EDUOPT,SOILCUTOFF, as SOIL-EQL)."""
    acc = np.asarray(acc, dtype=float)
    if len(acc) > nfft:
        raise ValueError(f"{len(acc)} values do not fit in NFFT = {nfft}")
    freqs = np.fft.rfftfreq(nfft, dt)
    A = np.fft.rfft(acc, n=nfft)
    if cutoff and cutoff > 0:
        A = np.where(freqs > cutoff + 1e-12, 0.0, A)
    a_dt = np.fft.irfft(A, n=nfft)
    B = A.copy()
    if nfft % 2 == 0:
        B[-1] *= 0.5
    a_fine = np.fft.irfft(B, n=up * nfft) * up
    return a_dt, a_fine


# =======================================================================================
# Time integration
# =======================================================================================
@dataclass
class Controls:
    """Resolved integration controls (SN-9, SN-10)."""

    nsub: int                 # base sub-increments per time step (1 for the flexible scheme)
    flexible: bool
    dg_max: float             # flexible scheme: maximum strain increment (fraction)
    tol_d: float
    tol_f: float
    maxit: int

    @classmethod
    def from_options(cls, o: NLSoilOptions, dt: float) -> "Controls":
        """SN-9: ``NSTimeSunInc`` = n >= 1 fixed sub-increments; 0 = flexible scheme with base sub-steps
        ``dt / ceil(dt / DT_MAX_FLEX)`` (16 sub-steps per period at 25 Hz) refined on strain increments."""
        if int(o.nsub) > 0:
            nsub, flex = int(o.nsub), False
        else:
            nsub, flex = max(1, int(math.ceil(dt / DT_MAX_FLEX - 1e-9))), True
        return cls(nsub=nsub, flexible=flex, dg_max=FLEX_STRAIN_INCREMENT,
                   tol_d=o.dispconv if o.dispconv > 0 else DEFAULT_DISPCONV,
                   tol_f=o.forceconv if o.forceconv > 0 else DEFAULT_FORCECONV,
                   maxit=o.equalit if o.equalit > 0 else DEFAULT_EQUALIT)


@dataclass
class SoilNonResult:
    """Histories and statistics of a SOIL-NON run (accelerations in length units / s^2)."""

    dt: float
    nfft: int
    input_motion: np.ndarray        # (nfft,) control motion actually applied (band-limited)
    acc_top: np.ndarray             # (nfft, ns+1) absolute acceleration at sublayer tops (+ base)
    strain_mid: np.ndarray          # (nfft, ns) strain (fraction) at mid-height
    stress_mid: np.ndarray          # (nfft, ns) total shear stress at mid-height (equilibrium)
    tau_h_mid: np.ndarray           # (nfft, ns) hysteretic (backbone/Masing) stress at mid-height
    gmax_elem: np.ndarray           # (ne,) max |strain| of every element over all sub-steps
    gmax_mid: np.ndarray            # (ns,) max |strain| at mid-height over all sub-steps
    disp_top: np.ndarray            # (nfft, ns+1) displacement relative to the input frame at tops
    stats: Dict[str, float] = field(default_factory=dict)

    def strain_history(self, layer: int) -> np.ndarray:
        return self.strain_mid[:, layer - 1]

    def stress_history(self, layer: int) -> np.ndarray:
        return self.stress_mid[:, layer - 1]

    def acc_history(self, layer: int) -> np.ndarray:
        return self.acc_top[:, layer - 1]


class _Step:
    """Converged (or failed) trial of one sub-step."""

    __slots__ = ("ok", "u", "v", "a", "mat", "iters", "dgmax", "rnorm")

    def __init__(self, ok, u=None, v=None, a=None, mat=None, iters=0, dgmax=0.0, rnorm=0.0):
        self.ok, self.u, self.v, self.a, self.mat = ok, u, v, a, mat
        self.iters, self.dgmax, self.rnorm = iters, dgmax, rnorm


class _Integrator:
    """Newmark average acceleration + Newton-Raphson for the SOIL-NON column (relative coordinates)."""

    def __init__(self, col: SoilNonColumn, damping: DampingModel, bedint: int, controls: Controls,
                 force_scale: float, disp_scale: float):
        self.col = col
        self.damp = damping
        self.bedint = int(bedint)
        self.ctl = controls
        self.nfree = col.ne + (1 if self.bedint == 1 else 0)
        self.m = col.mass[: self.nfree].copy()
        self.cb = col.base_rho * col.base_vs if self.bedint == 1 else 0.0
        self.mat = col.material()
        self.h = col.h
        self.u = np.zeros(self.nfree)
        self.v = np.zeros(self.nfree)
        self.a = np.zeros(self.nfree)
        self.ffloor = force_scale
        self.dfloor = disp_scale
        self.mmin = float(np.min(self.m))
        # |deformation(x)| <= dfac |x| (subtracting the base value, BedInt 1)
        self.dfac = 1.0 + math.sqrt(self.nfree) if self.bedint == 1 else 1.0
        if damping.kind == "tri":
            self.cdiag = damping.diag.copy()
            self.coff = damping.off.copy()
            if self.bedint == 1:
                self.cdiag[-1] += self.cb
        else:
            self.cfull = damping.full.copy()
            if self.bedint == 1:
                self.cfull[-1, -1] += self.cb

    # -- helpers -------------------------------------------------------------------------
    def strains(self, u: np.ndarray) -> np.ndarray:
        if self.bedint == 1:
            return (u[1:] - u[:-1]) / self.h
        full = np.append(u, 0.0)
        return (full[1:] - full[:-1]) / self.h

    def internal(self, tau: np.ndarray) -> np.ndarray:
        R = np.zeros(self.col.nn)
        R[:-1] -= tau
        R[1:] += tau
        return R[: self.nfree]

    def cmul(self, v: np.ndarray) -> np.ndarray:
        if self.damp.kind == "dense":
            return self.cfull @ v
        out = self.cdiag * v
        out[:-1] += self.coff * v[1:]
        out[1:] += self.coff * v[:-1]
        return out

    def solve(self, a0: float, a1: float, kt: np.ndarray, r: np.ndarray) -> np.ndarray:
        ke = kt / self.h
        kd = np.zeros(self.col.nn)
        kd[:-1] += ke
        kd[1:] += ke
        kd = kd[: self.nfree]
        ko = -ke[: self.nfree - 1]
        if self.damp.kind == "dense":
            J = a1 * self.cfull
            J[np.diag_indices(self.nfree)] += a0 * self.m + kd
            idx = np.arange(self.nfree - 1)
            J[idx, idx + 1] += ko
            J[idx + 1, idx] += ko
            return np.linalg.solve(J, r)
        ab = np.zeros((3, self.nfree))
        off = a1 * self.coff + ko
        ab[0, 1:] = off
        ab[1, :] = a0 * self.m + a1 * self.cdiag + kd
        ab[2, :-1] = off
        return _sla.solve_banded((1, 1), ab, r, overwrite_ab=True, overwrite_b=False, check_finite=False)

    def deformation(self, x: np.ndarray) -> np.ndarray:
        """Nodal displacements relative to the column base node (rigid-body part removed, SN-10)."""
        return x - x[-1] if self.bedint == 1 else x

    # -- one sub-step --------------------------------------------------------------------
    def try_step(self, dts: float, ag: float) -> _Step:
        gN, bN = NEWMARK_GAMMA, NEWMARK_BETA
        a0 = 1.0 / (bN * dts * dts)
        a1 = gN / (bN * dts)
        a2 = 1.0 / (bN * dts)
        a3 = 1.0 / (2.0 * bN) - 1.0
        un, vn, an = self.u, self.v, self.a
        # The Newton unknown is the step increment (not u_n+1): the accelerations and strains are then
        # computed without cancellation against a large rigid-body offset (BedInt 1, outcrop frame).
        inc = dts * vn + 0.5 * dts * dts * an
        g_committed = self.mat.g
        undef = self.deformation(un)
        fin = abs(ag) * float(np.sqrt(self.m @ self.m))
        for it in range(self.ctl.maxit + 1):
            g = g_committed + self.strains(inc)
            tr = self.mat.trial(g)
            acc = a0 * inc - a2 * vn - a3 * an
            vel = vn + dts * ((1.0 - gN) * an + gN * acc)
            R = self.internal(tr.tau)
            inert = self.m * (acc + ag)
            cv = self.cmul(vel)
            r = -inert - cv - R
            rn = float(np.sqrt(r @ r))
            fref = max(float(np.sqrt(inert @ inert)), float(np.sqrt(R @ R)), float(np.sqrt(cv @ cv)), fin,
                       self.ffloor)
            dinc = self.deformation(inc)
            uu = undef + dinc
            dref = max(float(np.sqrt(dinc @ dinc)), 1e-6 * float(np.sqrt(uu @ uu)), self.dfloor)
            if rn <= self.ctl.tol_f * fref:
                # The remaining correction du = J^-1 r is bounded by |r| / (a0 m_min) because
                # J = a0 M + a1 C + K_t with C and K_t positive semi-definite (tangents >= 0, s <= 1):
                # when the bound already meets the displacement criterion no solve is needed.
                if float(np.min(tr.kt)) >= 0.0 and rn * self.dfac / (a0 * self.mmin) <= self.ctl.tol_d * dref:
                    dg = float(np.max(np.abs(g - g_committed))) if g.size else 0.0
                    return _Step(True, un + inc, vel, acc, tr, it, dg, rn)
                du = self.solve(a0, a1, tr.kt, r)
                ddu = self.deformation(du)
                if float(np.sqrt(ddu @ ddu)) <= self.ctl.tol_d * dref:
                    dg = float(np.max(np.abs(g - g_committed))) if g.size else 0.0
                    return _Step(True, un + inc, vel, acc, tr, it, dg, rn)
            else:
                du = self.solve(a0, a1, tr.kt, r)
            if it == self.ctl.maxit or not np.all(np.isfinite(du)):
                break
            inc = inc + du
        return _Step(False, iters=self.ctl.maxit)

    def accept(self, st: _Step) -> None:
        self.u, self.v, self.a = st.u, st.v, st.a
        self.mat.commit(st.mat)


def run_soilnon(acc: np.ndarray, dt: float, nfft: int, col: SoilNonColumn, damping: DampingModel,
                bedint: int, controls: Controls, cutoff: float = 0.0, up: int = INPUT_UPSAMPLING,
                progress: Optional[Callable[[float], None]] = None,
                cancelled: Optional[Callable[[], bool]] = None) -> SoilNonResult:
    """Integrate the SOIL-NON column for the control motion ``acc`` (length units / s^2, ``len <= nfft``).

    ``acc`` is the base (within) motion for BedInt 0 and the outcrop motion for BedInt 1.  The response
    is computed at the NFFT time samples (the record is zero-padded as in SOIL-EQL); sub-steps follow
    SN-9 (fixed ``controls.nsub`` or flexible) with bisection on non-convergence (SN-10).
    """
    t0 = time.time()
    a_dt, a_fine = band_limited_input(acc, dt, nfft, cutoff, up)
    nf = len(a_fine)
    dtf = dt / up

    def ag_at(t: float) -> float:
        x = t / dtf
        i = int(math.floor(x + 1e-9))
        if i >= nf - 1:
            return float(a_fine[min(i, nf - 1)])
        w = x - i
        if w < 1e-9:
            return float(a_fine[i])
        return float(a_fine[i] * (1.0 - w) + a_fine[i + 1] * w)

    amax = float(np.max(np.abs(a_dt))) if a_dt.size else 0.0
    mtot = float(np.sum(col.mass))
    force_scale = 1e-12 * max(mtot * amax, 1e-300)
    disp_scale = 1e-14 * max(float(np.sum(col.thick)), 1e-300)
    integ = _Integrator(col, damping, bedint, controls, force_scale, disp_scale)
    ns = col.ns
    tops = col.top_node
    mids = col.mid_elem
    nfree = integ.nfree
    acc_top = np.zeros((nfft, ns + 1))
    disp_top = np.zeros((nfft, ns + 1))
    strain_mid = np.zeros((nfft, ns))
    stress_mid = np.zeros((nfft, ns))
    tau_mid = np.zeros((nfft, ns))
    gmax_elem = np.zeros(col.ne)

    ag0 = float(a_dt[0])
    integ.a = -ag0 * np.ones(nfree)          # equilibrium at rest: M a = -M 1 a_g(0)

    def record(k: int, ag: float) -> None:
        full_a = np.zeros(col.nn)
        full_a[:nfree] = integ.a
        abs_a = full_a + ag
        if bedint == 0:
            abs_a[-1] = ag
        acc_top[k] = abs_a[tops]
        full_u = np.zeros(col.nn)
        full_u[:nfree] = integ.u
        disp_top[k] = full_u[tops]
        g = integ.mat.g
        strain_mid[k] = g[mids]
        tau_mid[k] = integ.mat.t[mids]
        cum = np.cumsum(col.mass * abs_a)        # total stress in element e = inertia of nodes 0..e
        stress_mid[k] = cum[mids]

    record(0, ag0)
    units = controls.nsub * (1 << KMAX)       # sub-step units per time step
    k_level = 0
    stats = dict(substeps=0, iterations=0, max_iterations=0, bisections=0, strain_rejections=0,
                 min_substep=dt / controls.nsub, max_level=0)
    # Accelerate's BLAS raises spurious floating-point flags in matmul (numpy 2.0): the loop runs in one
    # errstate context; non-finite corrections are detected explicitly (a failed step is bisected).
    with np.errstate(over="ignore", invalid="ignore", divide="ignore", under="ignore"):
        for n in range(nfft - 1):
            if cancelled is not None and n % 64 == 0 and cancelled():
                raise SoilNonError("SOIL-NON run cancelled")
            t_n = n * dt
            p = 0
            while p < units:
                size = 1 << (KMAX - k_level)
                dts = dt * size / units
                ag = ag_at(t_n + dt * (p + size) / units)
                st = integ.try_step(dts, ag)
                stats["iterations"] += st.iters + 1
                if not st.ok:
                    if k_level >= KMAX:
                        raise SoilNonError(
                            f"equilibrium not reached within {controls.maxit} iterations at t = {t_n + dt * p / units:.5f} s "
                            f"even with {controls.nsub << KMAX} sub-increments per time step (increase EqualIt or "
                            f"the convergence tolerances)")
                    k_level += 1
                    stats["bisections"] += 1
                    continue
                if controls.flexible and st.dgmax > controls.dg_max and k_level < KMAX:
                    k_level = min(KMAX, k_level + max(1, int(math.ceil(math.log2(st.dgmax / controls.dg_max)))))
                    stats["strain_rejections"] += 1
                    continue
                integ.accept(st)
                np.maximum(gmax_elem, np.abs(st.mat.g), out=gmax_elem)
                stats["substeps"] += 1
                stats["max_iterations"] = max(stats["max_iterations"], st.iters + 1)
                stats["max_level"] = max(stats["max_level"], k_level)
                stats["min_substep"] = min(stats["min_substep"], dts)
                p += size
                # coarsen again when allowed (aligned with the coarser grid)
                if k_level > 0 and p % (2 * size) == 0:
                    if not controls.flexible or st.dgmax < 0.25 * controls.dg_max:
                        k_level -= 1
            record(n + 1, float(a_dt[n + 1]) if n + 1 < len(a_dt) else 0.0)
            if progress is not None and n % 256 == 0:
                progress((n + 1) / nfft)
    stats["cpu_s"] = time.time() - t0
    gmax_mid = gmax_elem[mids]
    return SoilNonResult(dt=dt, nfft=nfft, input_motion=a_dt, acc_top=acc_top, strain_mid=strain_mid,
                         stress_mid=stress_mid, tau_h_mid=tau_mid, gmax_elem=gmax_elem, gmax_mid=gmax_mid,
                         disp_top=disp_top, stats=stats)


# =======================================================================================
# Model resolution (NLSLAYER + DYNP -> SublayerModel) and equivalent-linear properties
# =======================================================================================
def resolve_sublayer(layer: int, data: Optional[NLLayerData], curve=None, ds: float = 0.0,
                     ) -> SublayerModel:
    """Sublayer model from its NLSLAYER set and DYNP curve (SN-1 ... SN-4).

    ``curve`` is a :class:`sassi.core.shake.DynamicProperty` (or None for a sublayer without DYNP label);
    ``ds`` the S-wave damping of its L property.
    """
    if data is None:
        raise SoilNonError(f"Error 125: no NLSLAYER set for soil sublayer {layer}")
    probs = data.problems()
    if probs:
        raise SoilNonError("; ".join(probs))
    if curve is not None and len(curve.d_strain):
        j = int(np.argmin(curve.d_strain))
        xi = float(curve.d_pct[j]) / 100.0
        xi_src = f"DYNP {curve.label} damping at {curve.d_strain[j]:g} %"
    else:
        xi = float(ds)
        xi_src = "L property S-wave damping"
    if data.curvefit == 1:
        if curve is None or len(curve.g_strain) == 0:
            raise SoilNonError(f"NLSLAYER {layer}: curve fit requested but sublayer {layer} has no DYNP "
                               f"G/Gmax curve (SPRO dynamic property)")
        fit = fit_mkz(curve.g_strain, curve.g_ratio)
        if fit.beta == 0:
            return SublayerModel(0.0, 1.0, float("inf"), xi, 0.0, "linear", fit, xi_src)
        return SublayerModel(fit.beta, fit.s, fit.gr_pct, xi, 0.0, "fit", fit, xi_src)
    if data.beta == 0:
        return SublayerModel(0.0, 1.0, float("inf"), xi, data.vis, "linear", None, xi_src)
    return SublayerModel(data.beta, data.s, data.refstrain, xi, data.vis, "user", None, xi_src)


def equivalent_properties(model: SublayerModel, gamma_eff_pct: float) -> Tuple[float, float]:
    """``(G/G0, damping)`` of the sublayer model at an effective strain (SN-12): secant modulus of the
    backbone and ``xi_min + D_Masing``."""
    if model.linear:
        return 1.0, model.xi_min
    gg = float(mkz_modulus_ratio(gamma_eff_pct, model.gr_pct, model.beta, model.s))
    dm = float(masing_damping(gamma_eff_pct, model.gr_pct, model.beta, model.s))
    return gg, model.xi_min + dm
