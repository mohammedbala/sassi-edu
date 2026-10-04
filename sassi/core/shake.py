"""SHAKE equivalent-linear free-field site response (numerical core of the SOIL module).

Normative sources: requirements section 4.11, decisions D-SOL-01 ... D-SOL-12, D-CNV-03;
theory R1 section 6 (SHAKE 1972 report and the public-domain SHAKE91 code).

Physics in one paragraph
------------------------
The soil column is a stack of horizontal, homogeneous, visco-elastic sublayers resting on a
half-space (the *last* SPRO entry, D-SOL-02).  A vertically propagating shear wave (or a
P wave for vertical input) is a superposition, in every sublayer m, of an up-going wave
``E_m`` and a down-going wave ``F_m``::

    u_m(z, t) = [E_m exp(+i k*_m z) + F_m exp(-i k*_m z)] exp(i w t)       z down from the layer top

with the complex wave number ``k*_m = w / V*_m`` and the complex velocity
``V*_m = sqrt(G*_m / rho_m)``.  Hysteretic damping enters through the complex modulus
``G* = G c(beta)`` (``c = 1 - 2 beta^2 + 2 i beta sqrt(1 - beta^2)`` by default, ``1 + 2 i beta``
with CMODFORM,1; D-CNV-03).  Continuity of displacement and shear stress at each
interface gives the SHAKE recursion (R1 section 6.1)::

    E_{m+1} = 1/2 E_m (1 + a_m) exp(+i k_m h_m) + 1/2 F_m (1 - a_m) exp(-i k_m h_m)
    F_{m+1} = 1/2 E_m (1 - a_m) exp(+i k_m h_m) + 1/2 F_m (1 + a_m) exp(-i k_m h_m)
    a_m = rho_m V*_m / (rho_{m+1} V*_{m+1})          (complex impedance ratio)

started from the free-surface condition ``E_1 = F_1`` (= 1).  The *within* motion at the
top of sublayer m is ``E_m + F_m``; the *outcrop* motion is ``2 E_m`` (the motion the
up-going wave alone would produce at a free surface).  Transfer functions are ratios of
these quantities, so any common scaling cancels.

The equivalent-linear method (Seed & Idriss) iterates: compute the response with the
current (G, beta), take the peak shear strain at mid-height of every sublayer, form the
effective strain ``gamma_eff = R_gamma * gamma_max`` and read new (G/Gmax, D) from the
strain-dependent curves, interpolated linearly in log10(strain) (D-SOL-03).  The
half-space is not iterated (SHAKE91 convention).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..conventions import cfactor

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------------------
# Strain-dependent curves
# ---------------------------------------------------------------------------------------
def interp_log_strain(strain_pct: Union[float, np.ndarray], xs: Sequence[float],
                      ys: Sequence[float]) -> Union[float, np.ndarray]:
    """Interpolate a curve y(strain) linearly in log10(strain), constant outside (D-SOL-03).

    ``xs`` are strains in percent (any order; duplicates are merged), ``ys`` the curve values
    (G/Gmax, or damping in percent).  This is SHAKE91's ``GN = AS*log10(SS) + BS`` rule.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    keep = xs > 0
    xs, ys = xs[keep], ys[keep]
    if xs.size == 0:
        raise ValueError("curve has no point with positive strain")
    order = np.argsort(xs, kind="stable")
    xs, ys = xs[order], ys[order]
    # merge exactly repeated abscissae (e.g. a curve padded with its last point)
    uniq, idx = np.unique(xs, return_index=True)
    xs, ys = uniq, ys[idx]
    x = np.log10(np.maximum(np.asarray(strain_pct, dtype=float), 1e-300))
    out = np.interp(x, np.log10(xs), ys, left=ys[0], right=ys[-1])
    if np.ndim(out) == 0:
        return float(out)
    return out


@dataclass
class DynamicProperty:
    """One DYNP label: modulus-reduction and damping curves (strain and damping in percent)."""

    label: str
    g_strain: np.ndarray        # strain %, G curve
    g_ratio: np.ndarray         # G/Gmax
    d_strain: np.ndarray        # strain %, damping curve
    d_pct: np.ndarray           # damping %

    def g_over_gmax(self, strain_pct):
        return interp_log_strain(strain_pct, self.g_strain, self.g_ratio)

    def damping(self, strain_pct):
        """Damping ratio (fraction) at the given strain (percent)."""
        return interp_log_strain(strain_pct, self.d_strain, self.d_pct) / 100.0


def curves_from_dynp_rows(rows: Sequence[Dict]) -> Dict[str, DynamicProperty]:
    """Build curves from the SOIL deck ``dynp`` table rows (label, kind 'G'|'D', strain, value).

    Labels keep the order of their first appearance (this is the FILE73 curve number).
    """
    tmp: Dict[str, Dict[str, List[Tuple[float, float]]]] = {}
    for r in rows:
        lab = str(r["label"])
        kind = str(r["kind"]).strip().upper()[:1]
        d = tmp.setdefault(lab, {"G": [], "D": []})
        if kind not in ("G", "D"):
            raise ValueError(f"DYNP {lab}: unknown curve kind {r['kind']!r} (expected G or D)")
        d[kind].append((float(r["strain"]), float(r["value"])))
    out: Dict[str, DynamicProperty] = {}
    for lab, d in tmp.items():
        g = np.asarray(d["G"], dtype=float).reshape(-1, 2)
        dd = np.asarray(d["D"], dtype=float).reshape(-1, 2)
        out[lab] = DynamicProperty(lab, g[:, 0], g[:, 1], dd[:, 0], dd[:, 1])
    return out


# ---------------------------------------------------------------------------------------
# Wave solution
# ---------------------------------------------------------------------------------------
def complex_velocity(G: np.ndarray, rho: np.ndarray, beta: np.ndarray, form: int = 0) -> np.ndarray:
    """``V* = sqrt(G c(beta) / rho)`` (principal root, Re V* > 0, Im V* >= 0)."""
    return np.sqrt(np.asarray(G, float) * cfactor(np.asarray(beta, float), form) / np.asarray(rho, float))


def wave_amplitudes(freqs: np.ndarray, thick: np.ndarray, rho: np.ndarray,
                    vstar: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Up-going/down-going amplitudes ``E, F`` (shape (nf, n)) with ``E_1 = F_1 = 1``.

    ``thick``, ``rho``, ``vstar`` have one entry per sublayer, the last being the half-space
    (its thickness is not used).  Recursion of R1 section 6.1 (SHAKE91 ``CXSOIL``/``AMP``),
    vectorised over the frequencies.
    """
    freqs = np.asarray(freqs, dtype=float)
    n = len(rho)
    w = 2.0 * np.pi * freqs
    E = np.empty((len(freqs), n), dtype=complex)
    F = np.empty((len(freqs), n), dtype=complex)
    E[:, 0] = 1.0
    F[:, 0] = 1.0
    for m in range(n - 1):
        alpha = rho[m] * vstar[m] / (rho[m + 1] * vstar[m + 1])
        ikh = 1j * w * thick[m] / vstar[m]          # i k*_m h_m
        ep = np.exp(ikh)
        em = np.exp(-ikh)
        E[:, m + 1] = 0.5 * E[:, m] * (1 + alpha) * ep + 0.5 * F[:, m] * (1 - alpha) * em
        F[:, m + 1] = 0.5 * E[:, m] * (1 - alpha) * ep + 0.5 * F[:, m] * (1 + alpha) * em
    if not (np.all(np.isfinite(E)) and np.all(np.isfinite(F))):
        raise FloatingPointError("SHAKE recursion overflow: reduce the cut-off frequency or the column depth")
    return E, F


def motion_at_top(E: np.ndarray, F: np.ndarray, layer: int, outcrop: bool) -> np.ndarray:
    """Motion amplitude at the top of sublayer ``layer`` (1-based): ``2E`` (outcrop) or ``E+F`` (within)."""
    j = layer - 1
    return 2.0 * E[:, j] if outcrop else E[:, j] + F[:, j]


def mid_layer_strain_factor(freqs: np.ndarray, E: np.ndarray, F: np.ndarray, thick: np.ndarray,
                            vstar: np.ndarray) -> np.ndarray:
    """Shear strain at mid-height of sublayers 1..n-1 per unit acceleration amplitude.

    With acceleration-wave amplitudes (R1 section 6.2, SHAKE91 ``STRT``)::

        gamma(w) = (E e^{i k h/2} - F e^{-i k h/2}) / (i w V*)

    The f = 0 term is set to zero (no static strain).  Shape (nf, n-1).
    """
    w = 2.0 * np.pi * np.asarray(freqs, dtype=float)
    ns = E.shape[1] - 1
    out = np.zeros((len(w), ns), dtype=complex)
    pos = w > 0
    for m in range(ns):
        ikh2 = 1j * w[pos] * 0.5 * thick[m] / vstar[m]
        out[pos, m] = (E[pos, m] * np.exp(ikh2) - F[pos, m] * np.exp(-ikh2)) / (1j * w[pos] * vstar[m])
    return out


# ---------------------------------------------------------------------------------------
# Column description and analysis
# ---------------------------------------------------------------------------------------
@dataclass
class SoilColumn:
    """Soil column for SOIL: sublayers top first, the last entry is the half-space (D-SOL-02).

    ``vel0`` and ``beta0`` are the low-strain wave velocity and damping used by the analysis
    (Vs, beta_s for horizontal input; Vp, beta_p for vertical input, indir = 1).
    ``labels`` are the DYNP labels ('' = linear sublayer, never iterated).
    """

    thick: np.ndarray
    weight: np.ndarray          # specific weight (force / length^3)
    vel0: np.ndarray
    beta0: np.ndarray
    gravity: float
    labels: List[str] = field(default_factory=list)

    def __post_init__(self):
        self.thick = np.asarray(self.thick, dtype=float)
        self.weight = np.asarray(self.weight, dtype=float)
        self.vel0 = np.asarray(self.vel0, dtype=float)
        self.beta0 = np.asarray(self.beta0, dtype=float)
        if not self.labels:
            self.labels = [""] * len(self.thick)

    @property
    def n(self) -> int:
        return len(self.thick)

    @property
    def rho(self) -> np.ndarray:
        return self.weight / self.gravity

    @property
    def gmax(self) -> np.ndarray:
        """Low-strain modulus rho V^2 (G_max for S waves, constrained modulus for P waves)."""
        return self.rho * self.vel0 ** 2

    @property
    def depth_top(self) -> np.ndarray:
        """Depth of the top of each sublayer (the last value is the top of the half-space)."""
        return np.concatenate([[0.0], np.cumsum(self.thick[:-1])])

    @property
    def depth_mid(self) -> np.ndarray:
        """Depth of the mid-height of the soil sublayers 1..n-1."""
        return self.depth_top[:-1] + 0.5 * self.thick[:-1]

    def average_velocity(self, vel: Optional[np.ndarray] = None) -> float:
        """Thickness-weighted average velocity of the soil sublayers (SHAKE91 'average shear velocity')."""
        v = self.vel0 if vel is None else np.asarray(vel)
        h = self.thick[:-1]
        return float(np.sum(h * v[:-1]) / np.sum(h)) if np.sum(h) > 0 else float(v[0])


def column_transfer(freqs: np.ndarray, col: SoilColumn, G: np.ndarray, beta: np.ndarray,
                    layer_out: int, outcrop_out: bool, layer_in: int, outcrop_in: bool,
                    form: int = 0) -> np.ndarray:
    """Transfer function  X(layer_out) / X(layer_in)  for given properties (G per sublayer).

    Each X is the outcrop (2E) or within (E+F) motion at the top of the sublayer.  Exact at
    any frequency (no FFT grid involved), e.g. SSAF output and VP-02a.
    """
    vstar = complex_velocity(G, col.rho, beta, form)
    E, F = wave_amplitudes(freqs, col.thick, col.rho, vstar)
    return motion_at_top(E, F, layer_out, outcrop_out) / motion_at_top(E, F, layer_in, outcrop_in)


@dataclass
class IterationRecord:
    """Properties and strains of one equivalent-linear iteration (all arrays over soil sublayers)."""

    number: int
    gamma_max: np.ndarray       # % (max |gamma(t)| at mid-height)
    gamma_eff: np.ndarray       # %
    G_used: np.ndarray
    G_new: np.ndarray
    beta_used: np.ndarray
    beta_new: np.ndarray

    @property
    def dG_pct(self) -> np.ndarray:
        """(new - used)/new * 100 (SHAKE91 error definition, R1 section 6.3)."""
        return np.where(self.G_new != 0, (self.G_new - self.G_used) / np.where(self.G_new != 0, self.G_new, 1) * 100, 0.0)

    @property
    def dbeta_pct(self) -> np.ndarray:
        return np.where(self.beta_new != 0,
                        (self.beta_new - self.beta_used) / np.where(self.beta_new != 0, self.beta_new, 1) * 100, 0.0)


@dataclass
class ShakeResult:
    """Final state of a SOIL analysis (everything needed for the outputs)."""

    freqs: np.ndarray           # FFT frequencies (nfft//2 + 1)
    nfft: int
    dt: float
    input_spectrum: np.ndarray  # rfft of the (scaled, cut-off) control motion, length units / s^2
    G: np.ndarray               # final properties used (all n sublayers incl. half-space)
    beta: np.ndarray
    E: np.ndarray               # wave amplitudes scaled to the input motion (nf, n), acceleration units
    F: np.ndarray
    vstar: np.ndarray
    strain_tf: np.ndarray       # mid-layer strain spectrum (nf, n-1), fraction
    gamma_max: np.ndarray       # % max |gamma(t)| of the final response
    gamma_eff: np.ndarray       # % ratio x gamma_max of the final response
    iterations: List[IterationRecord] = field(default_factory=list)

    @property
    def gamma_eff_compatible(self) -> np.ndarray:
        """Effective strain (%) with which the final ``G`` and ``beta`` are strain-compatible.

        The final properties were read from the curves at the effective strain of the *last
        iteration* (the response computed with the previous properties), so that strain is the
        one to pair with them, as SHAKE91 Table B-2 does (FILE88 and the final-properties
        table, R1 section 6.3).  :attr:`gamma_eff` is the strain of one more response computed
        with the final properties; the two agree once the iterations have converged.  Without
        any iteration (``iter`` = 0, vertical input) the properties are the low-strain ones and
        the strain of the final (linear) response is returned.
        """
        return (self.iterations[-1].gamma_eff if self.iterations else self.gamma_eff).copy()

    # -- frequency-domain motions --------------------------------------------------------
    def acc_spectrum(self, layer: int, outcrop: bool) -> np.ndarray:
        return motion_at_top(self.E, self.F, layer, outcrop)

    def acc_history(self, layer: int, outcrop: bool) -> np.ndarray:
        """Acceleration history at the top of ``layer`` (length units / s^2), NFFT samples."""
        return np.fft.irfft(self.acc_spectrum(layer, outcrop), n=self.nfft)

    def strain_history(self, layer: int) -> np.ndarray:
        """Shear (or axial, indir = 1) strain at mid-height of soil sublayer ``layer`` (fraction)."""
        return np.fft.irfft(self.strain_tf[:, layer - 1], n=self.nfft)

    def stress_history(self, layer: int, form: int = 0) -> np.ndarray:
        """Stress at mid-height: ``tau(w) = G* gamma(w)`` (model stress units)."""
        gs = self.G[layer - 1] * cfactor(self.beta[layer - 1], form)
        return np.fft.irfft(gs * self.strain_tf[:, layer - 1], n=self.nfft)


def _response(freqs, col: SoilColumn, G, beta, A_in, cl, outcrop, form):
    """Wave amplitudes scaled to the input spectrum ``A_in`` given at the top of sublayer ``cl``."""
    vstar = complex_velocity(G, col.rho, beta, form)
    E, F = wave_amplitudes(freqs, col.thick, col.rho, vstar)
    ctrl = motion_at_top(E, F, cl, outcrop)
    scale = np.zeros_like(A_in)
    nz = np.abs(ctrl) > 0
    scale[nz] = A_in[nz] / ctrl[nz]
    E = E * scale[:, None]
    F = F * scale[:, None]
    strain = mid_layer_strain_factor(freqs, E, F, col.thick, vstar)
    return E, F, vstar, strain


def run_shake(acc: np.ndarray, dt: float, nfft: int, col: SoilColumn,
              curves: Dict[str, DynamicProperty], cl: int, outcrop: bool, ratio: float,
              iterations: int, form: int = 0, cutoff: float = 0.0, iterate: bool = True,
              tol_pct: float = 0.0) -> ShakeResult:
    """Equivalent-linear SHAKE analysis (R1 sections 6.1-6.3, D-SOL-01/03/04).

    Parameters
    ----------
    acc : control motion in length units / s^2 (already scaled; ``len(acc) <= nfft``).
    cl, outcrop : the motion is applied at the top of sublayer ``cl`` as outcrop (2E) or
        within (E+F) motion; all transfer functions are scaled so the control motion equals
        the input exactly.
    ratio : effective/maximum strain ratio R_gamma.
    iterations : number of property updates; exactly this many are made (D-SOL-04) unless
        ``tol_pct > 0`` stops earlier (``EDUOPT,SOILTOL``).  The final response is always
        computed with the last updated properties ("adjusted last iteration", SHAKE91), so the
        reported motions, strains and FILE88 are mutually consistent.
    cutoff : > 0 removes Fourier components above this frequency (``EDUOPT,SOILCUTOFF``).
    iterate : False for vertical input (indir 1, R1 section 6.3): no property update.
    """
    acc = np.asarray(acc, dtype=float)
    if len(acc) > nfft:
        raise ValueError(f"{len(acc)} values do not fit in NFFT = {nfft}")
    freqs = np.fft.rfftfreq(nfft, dt)
    A_in = np.fft.rfft(acc, n=nfft)
    if cutoff and cutoff > 0:
        A_in = np.where(freqs > cutoff + 1e-12, 0.0, A_in)
    n = col.n
    if not (1 <= cl <= n):
        raise ValueError(f"control layer {cl} is not in 1..{n}")
    G = col.gmax.copy()
    beta = col.beta0.copy()
    records: List[IterationRecord] = []
    soil = np.arange(n - 1)
    iterable = np.array([bool(col.labels[i]) and col.labels[i] in curves for i in soil], dtype=bool)
    nit = int(iterations) if iterate else 0
    for it in range(1, nit + 1):
        E, F, vstar, strain = _response(freqs, col, G, beta, A_in, cl, outcrop, form)
        gmax_t = np.max(np.abs(np.fft.irfft(strain, n=nfft, axis=0)), axis=0) * 100.0
        geff = ratio * gmax_t
        G_new = G.copy()
        b_new = beta.copy()
        for i in soil[iterable]:
            cv = curves[col.labels[i]]
            G_new[i] = col.gmax[i] * cv.g_over_gmax(geff[i])
            b_new[i] = cv.damping(geff[i])
        records.append(IterationRecord(it, gmax_t, geff, G[:-1].copy(), G_new[:-1].copy(),
                                       beta[:-1].copy(), b_new[:-1].copy()))
        G, beta = G_new, b_new
        if tol_pct and tol_pct > 0:
            rec = records[-1]
            if max(np.max(np.abs(rec.dG_pct)), np.max(np.abs(rec.dbeta_pct))) < tol_pct:
                break
    E, F, vstar, strain = _response(freqs, col, G, beta, A_in, cl, outcrop, form)
    gmax_t = np.max(np.abs(np.fft.irfft(strain, n=nfft, axis=0)), axis=0) * 100.0
    return ShakeResult(freqs=freqs, nfft=nfft, dt=dt, input_spectrum=A_in, G=G, beta=beta, E=E, F=F,
                       vstar=vstar, strain_tf=strain, gamma_max=gmax_t, gamma_eff=ratio * gmax_t,
                       iterations=records)


def max_amplification(col: SoilColumn, G: np.ndarray, beta: np.ndarray, layer_out: int, outcrop_out: bool,
                      layer_in: int, outcrop_in: bool, fmax: float, form: int = 0,
                      df: float = 0.01) -> Tuple[float, float]:
    """Peak of |X(layer_out)/X(layer_in)| on (0, fmax]: grid search at ``df`` then golden refinement."""
    f = np.arange(df, fmax + 0.5 * df, df)
    H = np.abs(column_transfer(f, col, G, beta, layer_out, outcrop_out, layer_in, outcrop_in, form))
    k = int(np.argmax(H))
    a = f[max(k - 1, 0)]
    b = f[min(k + 1, len(f) - 1)]
    g = (np.sqrt(5.0) - 1) / 2

    def amp(x):
        return float(np.abs(column_transfer(np.array([x]), col, G, beta, layer_out, outcrop_out,
                                            layer_in, outcrop_in, form))[0])
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = amp(c), amp(d)
    for _ in range(60):
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = amp(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = amp(d)
    x = 0.5 * (a + b)
    return amp(x), x


# ---------------------------------------------------------------------------------------
# FILE73 / FILE88 (text, D-FIL-10)
# ---------------------------------------------------------------------------------------
FILE88_COLUMNS = ("layer", "thick", "gamma_eff_pct", "G", "Vs", "beta_s", "Vp", "beta_p")


def write_file88(path: PathLike, thick, gamma_eff, G, vs, beta_s, vp, beta_p, units: str = "",
                 halfspace: Optional[Sequence[float]] = None, title: str = "") -> Path:
    """Write FILE88: strain-compatible properties of the SOIL sublayers (half-space excluded).

    Format (D-FIL-10, read by SITE when SITEX soil mode = 1, matched by position, D-SOL-12)::

        # SASSI-EDU FILE88 v1 ...           comment lines start with '#'
        # columns: layer thick gamma_eff_pct G Vs beta_s Vp beta_p
        1  5.0  7.7e-04  3851.5  996.2  0.0071  1863.8  0.0071
        ...
        # halfspace thick weight Vp Vs beta_p beta_s   (informative comment)

    Damping values are ratios (fractions); gamma_eff in percent; G, Vs, Vp in model units.
    """
    path = Path(path)
    lines = ["# SASSI-EDU FILE88 v1  strain-compatible soil properties written by SOIL"]
    if title:
        lines.append(f"# title: {title}")
    if units:
        lines.append(f"# units: {units}")
    lines.append(f"# nlayers = {len(thick)}")
    lines.append("# columns: " + " ".join(FILE88_COLUMNS))
    for i in range(len(thick)):
        lines.append(f"{i + 1:5d} {thick[i]:.8e} {gamma_eff[i]:.8e} {G[i]:.8e} {vs[i]:.8e} "
                     f"{beta_s[i]:.8e} {vp[i]:.8e} {beta_p[i]:.8e}")
    if halfspace is not None:
        lines.append("# halfspace (not iterated): " + " ".join(f"{v:.8e}" for v in halfspace))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def read_file88(path: PathLike) -> Dict[str, np.ndarray]:
    """Read FILE88 into a dict of arrays keyed by :data:`FILE88_COLUMNS` (layer as int)."""
    rows = []
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        s = ln.strip()
        if not s or s[0] in "#*":
            continue
        rows.append([float(t) for t in s.split()])
    arr = np.asarray(rows, dtype=float).reshape(-1, len(FILE88_COLUMNS))
    out = {c: arr[:, j] for j, c in enumerate(FILE88_COLUMNS)}
    out["layer"] = out["layer"].astype(int)
    return out


def write_file73(path: PathLike, curves: Dict[str, DynamicProperty], title: str = "") -> Path:
    """Write FILE73: the soil material curves (text, D-FIL-10), read by STRESS (nonlinear SSI).

    Curve number ``ICURVE`` = 1-based order of the label in the DYNP table.  Layout::

        # SASSI-EDU FILE73 v1 ... strain in %, G/Gmax, damping in %
        CURVE <icurve> <label> <nG> <nD>
        G <strain %> <G/Gmax>        (nG lines)
        D <strain %> <damping %>     (nD lines)
    """
    path = Path(path)
    lines = ["# SASSI-EDU FILE73 v1  soil material curves written by SOIL",
             "# strain in percent; G/Gmax dimensionless; damping in percent"]
    if title:
        lines.append(f"# title: {title}")
    lines.append(f"# ncurves = {len(curves)}")
    for k, (lab, cv) in enumerate(curves.items(), start=1):
        lines.append(f"CURVE {k} {lab} {len(cv.g_strain)} {len(cv.d_strain)}")
        for x, y in zip(cv.g_strain, cv.g_ratio):
            lines.append(f"G {x:.8e} {y:.8e}")
        for x, y in zip(cv.d_strain, cv.d_pct):
            lines.append(f"D {x:.8e} {y:.8e}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def read_file73(path: PathLike) -> Dict[int, DynamicProperty]:
    """Read FILE73; returns {icurve: DynamicProperty} (insertion order = curve number)."""
    out: Dict[int, DynamicProperty] = {}
    cur = None
    g: List[Tuple[float, float]] = []
    d: List[Tuple[float, float]] = []

    def flush():
        if cur is not None:
            ga = np.asarray(g, float).reshape(-1, 2)
            da = np.asarray(d, float).reshape(-1, 2)
            out[cur[0]] = DynamicProperty(cur[1], ga[:, 0], ga[:, 1], da[:, 0], da[:, 1])

    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        s = ln.strip()
        if not s or s[0] in "#*":
            continue
        tok = s.split()
        if tok[0].upper() == "CURVE":
            flush()
            cur = (int(tok[1]), " ".join(tok[2:-2]))
            g, d = [], []
        elif tok[0].upper() == "G":
            g.append((float(tok[1]), float(tok[2])))
        elif tok[0].upper() == "D":
            d.append((float(tok[1]), float(tok[2])))
    flush()
    return out
