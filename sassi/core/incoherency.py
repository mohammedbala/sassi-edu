"""Incoherent seismic input, wave passage and multiple excitation: the incoherency factors of the
interaction nodes and the file FILE77 (requirements 4.4 item 6 and 4.6 item 5; decisions D-INC-02,
D-INC-03, D-INC-06 ... D-INC-12; spec 05b section 4, spec 05c A.5.9; ACS SASSI manual 6.5.4).

The idea (for the structural engineer)
--------------------------------------
For every SSI frequency and motion component ``c`` (X, Y, Z) HOUSE builds the coherency matrix of
the N interaction nodes, ``Sigma_ij = gamma_c(f, D_ij)`` (:mod:`sassi.core.coherency`; horizontal
projections of the nodes), and factorises it spectrally (Eq. 6.1 of the manual)::

    Sigma = Phi Lambda Phi^T,   lambda_1 >= lambda_2 >= ... >= 0,   sum_j lambda_j = trace(Sigma) = N

Column ``phi_k`` is the k-th *incoherent spatial mode*, ``lambda_k`` its variance.  The incoherent
free field is the coherent one (SITE, FILE1) times a complex factor per interaction node,
``u_i,inc(w) = s_i(w) u_i,coh(w)``, with ``s`` synthesised from the modes (requirements 4.4 item 6):

* stochastic simulation (SS, "Simulation Mean"):  ``s = Phi_m Lambda_m^(1/2) exp(i theta)`` with
  independent random phases ``theta_k ~ U[-RandPhz, +RandPhz]`` per frequency, mode, direction and
  sample, so that the ensemble mean of ``s s^H`` is ``Sigma`` (D-INC-07; VP-I2);
* deterministic algebraic sum (AS, Linear):  ``s = sum_k sqrt(lambda_k) phi_k`` (zero phases);
* single mode (``nmodes = -k``, SRSS approaches):  ``s = sqrt(lambda_k) phi_k``.

Eigenvectors are defined up to their sign; D-INC-03 fixes it so that ``sum_i phi_ik >= 0``
("with ATF phase adjustment", ``EDUOPT,INCOHSIGN,ADJUST``, the default; ``RAW`` keeps the eigensolver
signs).  Then for ``f -> 0`` (``Sigma -> 1 1^T``, ``lambda_1 = N``, ``phi_1 = 1/sqrt(N)``) the AS
factor is ``s = 1``: coherent motion and a zero-frequency ATF of 1.00, the manual's check.  A
mode whose sum vanishes (an antisymmetric mode of a symmetric layout) is oriented by its first
significant component (positive).  Equal eigenvalues (mirror-image modes of a symmetric layout of
interaction nodes) leave the eigenvectors undetermined within their eigenspace, and the AS sum
depends on that choice: with ADJUST the eigenspace gets a canonical basis that depends only on the
subspace -- first the projection of the uniform field, then pivoted orthogonalisation in node order
(:func:`canonical_basis`; an implementation decision completing D-INC-03), so the factors do not
depend on round-off or on the eigensolver.  With wave passage or near-equal eigenvalues the AS
factors remain sensitive to small changes of the coherency matrix, as the manual warns for the
deterministic approaches.

Wave passage multiplies ``s_i`` by ``exp(-i w tau_i)``, ``tau_i = ((x_i - x_c) cos(a) + (y_i - y_c)
sin(a)) / V_app``, the arrival delay along Line D (angle ``a`` = WPASS <ang>) measured from the
ANALYS control point (D-INC-06); multiple excitation multiplies the factors of the interaction nodes
of zone ``k`` by its spectral amplification ratio ``SAR_k(w)`` (``SAR = 1`` outside every zone),
after incoherency and wave passage (D-INC-10).

FILE77 (``FILE77``, or ``FILE77001`` ... ``FILE77050`` for the stochastic samples) stores ``s``
per frequency, direction X/Y/Z and interaction node (:data:`sassi.io.files.SCHEMAS`); ANALYS
multiplies the free-field load (FFL) or motion (FFM) by it (requirements 4.6 item 5).

Random numbers (D-INC-07): the X and Y streams are the two children of ``SeedSequence(HSeed)``, the
Z stream ``SeedSequence(VSeed)``; every stream has one PCG64 child generator per sample, which
draws N phases (one per mode, descending lambda) per frequency in ascending frequency order.  The
phases of sample ``r`` are therefore independent of the number of samples computed together and
the results are reproducible.  A negative seed ``-n`` uses the entropy ``(n, 1)``, ``n >= 0`` the
entropy ``(n, 0)``.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..io.container import Container, read_container, write_container
from ..io.files import validate
from .coherency import CoherencyError, CoherencySpec, DIRECTIONS
from .tlm import quiet_fpe

#: requirements 4.4 item 6 / spec 05b 3.17: at most 50 stochastic simulations (FILE77001..FILE77050)
MAX_SIM = 50
#: trace check of Eq. 6.1 (requirements 4.4 item 6, VP-27)
TRACE_TOL = 1e-8
#: cumulative variance of the mode-count criterion (manual 6.5.4, Eq. 6.4)
MODE_FRACTION = 0.90
SIGN_ADJUST, SIGN_RAW = "ADJUST", "RAW"
#: optional EDUOPT file of HOUSE written by AFWRITE (INCOHSIGN, INCOHMERGE; D-INC-03, D-INC-09)
OPTION_FILE = "HOUSE.opt"
#: relative tolerance of a vanishing mode sum (sign convention tie)
SIGN_TIE = 1e-9


class IncoherencyError(ValueError):
    """Invalid incoherency / wave-passage / multiple-excitation input."""


# ======================================================================================
# Names
# ======================================================================================
def file77_name(sim: Optional[int] = None) -> str:
    """``FILE77`` (deterministic, wave passage, multiple excitation) or ``FILE77sss`` (stochastic
    sample ``sim`` = 1..50; requirements 1.8)."""
    return "FILE77" if sim is None else f"FILE77{int(sim):03d}"


def stochastic(hseed: int, vseed: int, randphz: float) -> bool:
    """D-INC-02: stochastic simulation iff (HSeed != 0 or VSeed != 0) and RandPhz > 0."""
    return (int(hseed) != 0 or int(vseed) != 0) and float(randphz) > 0.0


# ======================================================================================
# Spectral factorisation of the coherency matrix
# ======================================================================================
@dataclass
class Decomposition:
    """``Sigma = Phi Lambda Phi^T`` of one frequency and direction.

    ``lam`` (N,) descending, negative round-off clipped to 0; ``phi`` (N, N) columns = modes with the
    sign convention applied; ``trace_error = |sum(lam) - N| / N`` (Eq. 6.1); ``min_raw`` the most
    negative raw eigenvalue divided by N (0 when Sigma is positive semi-definite)."""

    lam: np.ndarray
    phi: np.ndarray
    trace_error: float
    min_raw: float
    n_clipped: int

    @property
    def n(self) -> int:
        return int(self.lam.size)

    def contributions(self) -> np.ndarray:
        """Percent contribution ``upsilon_j = 100 lambda_j / N`` of every mode (Eq. 6.3)."""
        return 100.0 * self.lam / max(self.n, 1)

    def modes_for(self, fraction: float = MODE_FRACTION) -> int:
        """Smallest m with a cumulative contribution >= ``fraction`` (Eq. 6.4, 90 % criterion)."""
        c = np.cumsum(self.lam) / max(self.n, 1)
        hit = np.flatnonzero(c >= fraction - 1e-12)
        return int(hit[0]) + 1 if hit.size else self.n


def normalise_signs(phi: np.ndarray) -> np.ndarray:
    """Orient every column so that its sum is >= 0 (D-INC-03), in place.  A column whose sum vanishes
    (``|sum| <= 1e-9 sum |phi|``) is oriented so that its first significant component
    (``|phi_ik| > 1e-6 max_i |phi_ik|``) is positive, which makes the orientation deterministic."""
    if phi.size == 0:
        return phi
    sums = phi.sum(axis=0)
    scale = np.abs(phi).sum(axis=0)
    sgn = np.where(sums < 0.0, -1.0, 1.0)
    ties = np.flatnonzero(np.abs(sums) <= SIGN_TIE * np.maximum(scale, 1e-300))
    for k in ties:
        col = phi[:, k]
        big = np.flatnonzero(np.abs(col) > 1e-6 * np.abs(col).max()) if np.any(col) else np.zeros(0, int)
        sgn[k] = -1.0 if big.size and col[big[0]] < 0.0 else 1.0
    phi *= sgn[None, :]
    return phi


@quiet_fpe
def canonical_basis(Q: np.ndarray, tie: float = 1e-9) -> np.ndarray:
    """A basis of the subspace spanned by the orthonormal columns ``Q`` (N, m) that depends on the
    subspace only, not on the basis the eigensolver happened to return.

    Degenerate eigenvalues are common (a symmetric layout of interaction nodes gives equal lambda for
    mirror-image modes) and any rotation of their eigenvectors is an equally valid eigen-solution,
    but the AS sum ``sum_k sqrt(lambda_k) phi_k`` and the single modes of ``<nmodes> = -k`` depend on
    it.  The canonical basis is: first the projection of the uniform field ``1`` (the coherent
    motion) on the subspace, when it is not zero; then, for the rest of the subspace, pivoted
    orthogonalisation in node order -- the node with the largest projector diagonal (the first
    within a relative tie ``tie``) gives the next vector ``P e_j / |P e_j|``, which is removed before
    the next pivot.  Every step uses only the projector ``P = Q Q^T``, so the result is invariant
    under ``Q -> Q R`` (R orthogonal)."""
    Q = np.asarray(Q, float)
    n, m = Q.shape
    if m <= 1:
        return Q.copy()
    out = []
    p = Q.T @ np.ones(n)
    pn = float(np.linalg.norm(p))
    if pn > 1e-8 * math.sqrt(n):
        c = p / pn
        out.append(Q @ c)
        # orthonormal complement of c in R^m (Householder reflector mapping c to e_1)
        e1 = np.zeros(m)
        e1[0] = 1.0
        v = c - e1 if c[0] <= 0 else c + e1
        H = np.eye(m) - 2.0 * np.outer(v, v) / float(v @ v)
        Q = Q @ H[:, 1:]
    A = Q.T.copy()                                   # (k, N): coordinates of the projector columns
    for _ in range(A.shape[0]):
        norms = np.einsum("ij,ij->j", A, A)
        top = float(norms.max())
        if not top > 0:
            break
        j = int(np.flatnonzero(norms >= top * (1.0 - tie))[0])
        a = A[:, j] / math.sqrt(norms[j])
        out.append(Q @ a)
        A -= np.outer(a, a @ A)
    return np.column_stack(out)


def eigen_clusters(lam: np.ndarray, rtol: float = 1e-9, floor: float = 1e-12) -> List[Tuple[int, int]]:
    """Index ranges [i, j) of (numerically) equal eigenvalues among ``lam`` (descending) larger than
    ``floor`` x lam[0]: consecutive values within ``rtol`` x lam[0]."""
    out = []
    n = lam.size
    if n == 0:
        return out
    scale = max(float(lam[0]), 1e-300)
    i = 0
    while i < n:
        j = i + 1
        while j < n and abs(lam[j - 1] - lam[j]) <= rtol * scale:
            j += 1
        if j - i > 1 and lam[i] > floor * scale:
            out.append((i, j))
        i = j
    return out


@quiet_fpe
def decompose(S: np.ndarray, sign: str = SIGN_ADJUST) -> Decomposition:
    """Spectral factorisation of a real symmetric coherency matrix (LAPACK syevd through
    ``numpy.linalg.eigh``; requirements 4.4 item 6): lambda descending, negative round-off clipped
    to 0, and with ``sign`` = ``ADJUST`` (D-INC-03, the default) the deterministic modes: every group
    of equal eigenvalues gets its canonical basis (:func:`canonical_basis`) and every mode the sign
    convention ``sum_i phi_ik >= 0`` (:func:`normalise_signs`); ``RAW`` keeps the eigensolver's
    vectors and signs."""
    S = np.asarray(S, float)
    n = S.shape[0]
    if n == 0:
        return Decomposition(np.zeros(0), np.zeros((0, 0)), 0.0, 0.0, 0)
    w, V = np.linalg.eigh(S)
    lam = w[::-1].copy()
    phi = V[:, ::-1].copy()
    neg = lam < 0.0
    min_raw = float(min(0.0, lam.min())) / n
    lam[neg] = 0.0
    if str(sign).upper() != SIGN_RAW:
        for i, j in eigen_clusters(lam):
            phi[:, i:j] = canonical_basis(phi[:, i:j])
        normalise_signs(phi)
    return Decomposition(lam, phi, abs(float(lam.sum()) - n) / n, min_raw, int(neg.sum()))


def select_modes(nmodes: int, n: int) -> np.ndarray:
    """0-based mode indices of INCOH <nmodes> (requirements 4.4 item 6, spec 05b 3.13): 0 = all N
    modes; k > 0 = the k largest (all when k >= N); -k = mode k only (k <= N)."""
    nmodes = int(nmodes)
    if nmodes == 0:
        return np.arange(n)
    if nmodes > 0:
        return np.arange(min(nmodes, n))
    k = -nmodes
    if k > n:
        raise IncoherencyError(f"<nmodes> = {nmodes}: mode {k} does not exist ({n} incoherent modes = number of "
                               "interaction nodes)")
    return np.array([k - 1])


@quiet_fpe
def synthesize(dec: Decomposition, modes: np.ndarray, phases: Optional[np.ndarray] = None) -> np.ndarray:
    """Incoherency factors ``s`` (N,) complex from the selected modes: ``sum_k sqrt(lambda_k) phi_k
    exp(i theta_k)`` (``phases`` = theta of every mode, radians; None = zero phases, the AS / single
    mode cases)."""
    amp = np.sqrt(dec.lam[modes])
    if phases is None:
        return (dec.phi[:, modes] @ amp).astype(complex)
    return dec.phi[:, modes] @ (amp * np.exp(1j * np.asarray(phases, float)[modes]))


# ======================================================================================
# Random phases (D-INC-07)
# ======================================================================================
def _entropy(seed: int) -> List[int]:
    s = int(seed)
    return [abs(s), 1 if s < 0 else 0]


def phase_generators(hseed: int, vseed: int, nsim: int) -> List[List[np.random.Generator]]:
    """``gen[d][r]``: the PCG64 generator of direction d (0 X, 1 Y, 2 Z) and sample r (0-based)."""
    hx, hy = np.random.SeedSequence(_entropy(hseed)).spawn(2)
    hz = np.random.SeedSequence(_entropy(vseed))
    return [[np.random.Generator(np.random.PCG64(c)) for c in ss.spawn(int(nsim))] for ss in (hx, hy, hz)]


# ======================================================================================
# Wave passage and multiple excitation
# ======================================================================================
def arrival_delays(xy, appv: float, angle_deg: float, xc: float = 0.0, yc: float = 0.0) -> np.ndarray:
    """Arrival delay along Line D, ``tau_i = ((x_i - x_c) cos a + (y_i - y_c) sin a) / V_app``
    (requirements 4.4 item 6, D-INC-06: tau = 0 at the ANALYS control point)."""
    if not appv > 0:
        raise IncoherencyError(f"Error 113: apparent velocity for Line D <appv> = {appv:g} must be > 0")
    p = np.asarray(xy, float).reshape(-1, 2)
    a = math.radians(float(angle_deg))
    return ((p[:, 0] - xc) * math.cos(a) + (p[:, 1] - yc) * math.sin(a)) / float(appv)


def wave_passage_factors(freq: Sequence[float], tau: np.ndarray) -> np.ndarray:
    """``exp(-i w tau_i)`` (nF, N) with ``w = 2 pi f`` (e^{+i w t} convention: a later arrival is a
    phase lag)."""
    w = 2.0 * math.pi * np.asarray(freq, float)
    return np.exp(-1j * w[:, None] * np.asarray(tau, float)[None, :])


@dataclass
class MEZone:
    """One multiple-excitation zone: ME <no>, the interaction nodes ``nfirst..nlast`` and its SAR."""

    no: int
    nfirst: int
    nlast: int
    nodes: np.ndarray            # positions in the interaction-node list
    sar: np.ndarray              # (nF,) complex spectral amplification ratios


def me_zones(zone_rows: Sequence[dict], amp_rows: Sequence[dict], int_node: Sequence[int], nF: int,
             cmplxspec: int) -> Tuple[List[MEZone], List[str], List[str]]:
    """Multiple-excitation zones from the HOUSE deck tables ``me`` (no, nfirst, nlast) and ``amp``
    (no, idx, re, im) (requirements 4.4 item 6, D-INC-10, spec 05b 3.19).

    A zone holds the interaction nodes with ``nfirst <= id <= nlast``; ``nfirst`` and ``nlast`` must
    be interaction nodes (Errors 116/117); zones may not overlap; every zone needs one ratio per SSI
    frequency (Error 119) with modulus in [0, 10] (Error 118); real ratios unless ``cmplxspec`` = 1.
    Returns (zones, errors, warnings)."""
    errors: List[str] = []
    warnings: List[str] = []
    ids = [int(n) for n in int_node]
    pos = {n: k for k, n in enumerate(ids)}
    if not zone_rows:
        errors.append("Error 115: multiple excitation is on (HOUSE <me> = 1) but no ME zone is defined")
        return [], errors, warnings
    amp: Dict[int, Dict[int, complex]] = {}
    for r in amp_rows:
        amp.setdefault(int(r["no"]), {})[int(r["idx"])] = complex(float(r["re"]), float(r["im"]))
    zones: List[MEZone] = []
    owner: Dict[int, int] = {}
    for r in sorted(zone_rows, key=lambda r: int(r["no"])):
        no, n1, n2 = int(r["no"]), int(r["nfirst"]), int(r["nlast"])
        if n1 not in pos:
            errors.append(f"Error 116: Illegal First Node Number for Motion {no} (node {n1} is not an interaction node)")
            continue
        if n2 not in pos or n2 < n1:
            errors.append(f"Error 117: Illegal Last Node Number for Motion {no} (node {n2})")
            continue
        members = np.array(sorted(pos[n] for n in ids if n1 <= n <= n2), dtype=np.int64)
        clash = sorted({owner[k] for k in members if k in owner})
        if clash:
            errors.append(f"multiple-excitation motion {no} (nodes {n1}..{n2}) overlaps motion(s) {clash}")
            continue
        for k in members:
            owner[int(k)] = no
        vals = amp.get(no, {})
        if len(vals) != nF or sorted(vals) != list(range(1, nF + 1)):
            errors.append(f"Error 119: Spectral Amplification Ratios for Motion {no} Do Not Match Frequencies "
                          f"({len(vals)} ratios, {nF} SSI frequencies)")
            continue
        sar = np.array([vals[k] for k in range(1, nF + 1)], complex)
        if not int(cmplxspec):
            if np.any(sar.imag != 0.0):
                warnings.append(f"motion {no}: imaginary parts of the SAR ignored (HOUSE <cmplxspec> = 0)")
            sar = sar.real.astype(complex)
        bad = np.flatnonzero((np.abs(sar) > 10.0) | (sar.real < 0.0) & (not int(cmplxspec)))
        if bad.size:
            errors.append(f"Error 118: Illegal Spectral Amplification Ratio for Motion {no} "
                          f"({sar[bad[0]]:.4g} at frequency index {int(bad[0]) + 1}; [0, 10] on the modulus)")
            continue
        if no > 10:
            warnings.append(f"ME motion number {no} > 10 (manual command range 1-10; accepted, D-INC-10)")
        zones.append(MEZone(no, n1, n2, members, sar))
    if len(zone_rows) > 5000:
        warnings.append(f"{len(zone_rows)} multiple-excitation zones exceed the manual limit of 5000")
    return zones, errors, warnings


def sar_factors(zones: Sequence[MEZone], nF: int, n: int) -> np.ndarray:
    """``SAR_zone(i)(w)`` (nF, N): the ratio of the zone of each interaction node, 1 outside zones."""
    out = np.ones((nF, n), complex)
    for z in zones:
        out[:, z.nodes] = z.sar[:, None]
    return out


# ======================================================================================
# Settings of one HOUSE run
# ======================================================================================
@dataclass
class IncoherencySettings:
    """HOUSE <coh>, <wpass>, <me>, <cmplxspec>, INCOH, WPASS, HOUSEX <supmode>/<nsim> and the
    EDUOPT switches INCOHSIGN / INCOHMERGE of one HOUSE run."""

    coh: int = 0
    wpass: int = 0
    me: int = 0
    cmplxspec: int = 0
    cohf: int = 1
    gamma: Tuple[float, float, float] = (0.1, 0.1, 0.2)
    alpha: float = 0.5
    ngp: int = 1
    ipr: int = 0
    nmodes: int = 0
    met: int = 0
    hseed: int = 0
    vseed: int = 0
    randphz: float = 0.0
    supmode: int = 0
    nsim: int = 1
    appv: float = 1.0e9
    wang: float = 0.0
    xc: float = 0.0
    yc: float = 0.0
    sign: str = SIGN_ADJUST
    merge: int = 0

    @classmethod
    def from_deck(cls, d, options: Optional[Dict[str, str]] = None) -> "IncoherencySettings":
        o = options or {}
        return cls(coh=int(d["coh"]), wpass=int(d["wpass"]), me=int(d["me"]), cmplxspec=int(d["cmplxspec"]),
                   cohf=int(d["cohf"]), gamma=(float(d["gammax"]), float(d["gammay"]), float(d["gammaz"])),
                   alpha=float(d["alpha"]), ngp=int(d["ngp"]), ipr=int(d["ipr"]), nmodes=int(d["nmodes"]),
                   met=int(d["met"]), hseed=int(d["hseed"]), vseed=int(d["vseed"]), randphz=float(d["randphz"]),
                   supmode=int(d["supmode"]), nsim=int(d["nsim"]), appv=float(d["appv"]), wang=float(d["wang"]),
                   xc=float(d["xc"]), yc=float(d["yc"]),
                   sign=str(o.get("INCOHSIGN", SIGN_ADJUST)).upper(), merge=int(str(o.get("INCOHMERGE", "0")) == "1"))

    @property
    def active(self) -> bool:
        """FILE77 is written when the motion is incoherent, with wave passage or multiple excitation."""
        return bool(self.coh or self.wpass or self.me)

    @property
    def stochastic(self) -> bool:
        return bool(self.coh) and stochastic(self.hseed, self.vseed, self.randphz)

    @property
    def method(self) -> str:
        """'SS' stochastic simulation, 'AS' algebraic sum, 'SINGLE' one mode, 'COHERENT' (wave
        passage / multiple excitation of coherent motion)."""
        if not self.coh:
            return "COHERENT"
        if self.stochastic:
            return "SS"
        return "SINGLE" if self.nmodes < 0 else "AS"

    @property
    def n_files(self) -> int:
        return int(self.nsim) if self.stochastic else 1

    def describe_method(self) -> str:
        m = self.method
        if m == "COHERENT":
            return "coherent motion (no incoherency): FILE77 carries the wave-passage / multiple-excitation factors"
        if m == "SS":
            return (f"stochastic simulation (Simulation Mean): {self.nsim} sample(s), random modal phases in "
                    f"[-{self.randphz:g}, +{self.randphz:g}] deg, HSeed {self.hseed} (X, Y), VSeed {self.vseed} (Z)")
        if m == "SINGLE":
            return (f"single incoherent mode {-self.nmodes} (<nmodes> = {self.nmodes}; "
                    f"{'Quadratic: SRSS TF' if self.supmode else 'Linear: SRSS FRS'} approach, one run per mode)")
        return "deterministic algebraic sum of the scaled modes (AS, Linear; zero modal phases)"

    def validate(self, n_int: int, dim: int, symm: bool, embedded_levels: int) -> Tuple[List[str], List[str]]:
        """Errors and warnings of the option combination (requirements 4.4 item 6, spec 05b section 8,
        D-INC-02, D-INC-08, D-INC-10)."""
        errors: List[str] = []
        warnings: List[str] = []
        if self.coh:
            if dim != 2:
                errors.append("EDU-26: incoherent analysis needs a 3D model (HOUSE <dim> = 2; G-17)")
            if symm:
                errors.append("EDU-26: incoherent analysis is not allowed with SYMM planes (full models only; G-17)")
            if n_int == 0:
                errors.append("incoherent analysis needs interaction nodes")
            if not 1 <= self.cohf <= 7:
                errors.append(f"Error 114: unlagged coherency model <cohf> = {self.cohf} must be 1..7 (D-INC-01)")
            elif 2 <= self.cohf <= 7 and not self.wpass:
                errors.append(f"EDU-26: coherency model {self.cohf} is applied only with wave passage (HOUSE <wpass> = 1; "
                              "manual 6.5.4) -- set <wpass> = 1 (V_app = 1e9 suppresses the delays)")
            if self.supmode not in (0, 1):
                errors.append(f"HOUSEX <supmode> = {self.supmode} must be 0 (Linear) or 1 (Quadratic)")
            if self.stochastic:
                if not 1 <= self.nsim <= MAX_SIM:
                    errors.append(f"EDU-26: {self.nsim} stochastic simulations: HOUSEX <nsim> must be 1..{MAX_SIM}")
                if self.nmodes < 0:
                    errors.append(f"EDU-26: <nmodes> = {self.nmodes} (one mode) cannot be combined with stochastic "
                                  "simulation (spec 05b section 8)")
                if self.supmode == 1:
                    errors.append("EDU-26: Quadratic (SRSS) superposition is deterministic: set HSeed = VSeed = 0 or "
                                  "RandPhz = 0, or use Linear superposition for stochastic simulation")
                if self.nmodes > 0:
                    warnings.append(f"stochastic simulation with the first {self.nmodes} modes only (all modes are "
                                    "recommended, spec 05b 3.13)")
                if self.randphz > 180.0:
                    warnings.append(f"random phase angle {self.randphz:g} deg > 180 (phases wrap around)")
                if int(self.hseed) == 0 or int(self.vseed) == 0:
                    warnings.append("stochastic simulation with a zero seed: the manual asks for non-zero HSeed and VSeed")
            else:
                if self.supmode == 1 and self.nmodes >= 0:
                    errors.append(f"EDU-26: Quadratic (SRSS TF) superposition needs one HOUSE + ANALYS run per mode: "
                                  f"set <nmodes> = -k (k = 1, 2, ...) instead of {self.nmodes}")
                if self.nsim > 1:
                    warnings.append(f"HOUSEX <nsim> = {self.nsim} is used only by stochastic simulation "
                                    "(HSeed/VSeed and RandPhz are zero): one deterministic FILE77 is written")
                warnings.append("EDU-13: Deterministic AS/SRSS Incoherency Is Valid for Rigid Foundations Only")
            if self.nmodes < 0 and -self.nmodes > max(n_int, 1) and n_int:
                errors.append(f"<nmodes> = {self.nmodes}: mode {-self.nmodes} does not exist ({n_int} interaction nodes)")
            if embedded_levels > 0:
                if self.ngp <= 0:
                    errors.append(f"Error 60: Illegal Number of Mesh Points / Embedment Level (<ngp> = {self.ngp}; "
                                  "embedded interaction nodes, D-INC-08)")
                elif self.ngp != embedded_levels:
                    warnings.append(f"INCOH <ngp> = {self.ngp} but the interaction nodes lie on {embedded_levels} "
                                    "embedded level(s) below grade (D-INC-08: informative)")
        if self.wpass and not self.appv > 0:
            errors.append(f"Error 113: Illegal Apparent Velocity for Line D (<appv> = {self.appv:g})")
        if self.me:
            if not self.wpass:
                errors.append("EDU-26: multiple excitation needs wave passage (HOUSE <wpass> = 1; manual 6.5.4)")
        if self.sign not in (SIGN_ADJUST, SIGN_RAW):
            errors.append(f"EDUOPT,INCOHSIGN,{self.sign}: must be ADJUST or RAW (D-INC-03)")
        return errors, warnings


def read_options(workdir: Union[str, Path]) -> Dict[str, str]:
    """Optional ``HOUSE.opt`` of the model directory (written by AFWRITE): EDUOPT key/value lines,
    e.g. ``EDUOPT,INCOHSIGN,RAW`` and ``EDUOPT,INCOHMERGE,1`` (D-INC-03, D-INC-09)."""
    opts = {"INCOHSIGN": SIGN_ADJUST, "INCOHMERGE": "0"}
    p = Path(workdir) / OPTION_FILE
    if not p.exists():
        return opts
    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s[0] in "*#!":
            continue
        toks = [t for t in re.split(r"[,\s]+", s) if t]
        if toks and toks[0].upper() == "EDUOPT":
            toks = toks[1:]
        if len(toks) >= 2:
            opts[toks[0].upper()] = toks[1].upper()
    return opts


def options_text(sign: str, merge: str) -> str:
    """Content of ``HOUSE.opt`` (AFWRITE)."""
    return ("* SASSI-EDU HOUSE options without a deck parameter (written by AFWRITE)\n"
            f"EDUOPT,INCOHSIGN,{str(sign).upper()}\nEDUOPT,INCOHMERGE,{merge}\n")


# ======================================================================================
# Plan positions (D-INC-09)
# ======================================================================================
def unique_positions(xy: np.ndarray, tol: float) -> Tuple[np.ndarray, np.ndarray]:
    """Distinct plan positions of ``xy`` (n, 2) within ``tol`` (``EDUOPT,INCOHMERGE,1``): returns
    (positions (m, 2), index of the position of every point (n,)).  The first point of every group
    gives its position (deterministic)."""
    p = np.asarray(xy, float).reshape(-1, 2)
    if p.shape[0] == 0:
        return p, np.zeros(0, np.int64)
    key = np.round(p / max(tol, 1e-300)).astype(np.int64)
    uniq, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    order = np.argsort(first)                       # keep the first-appearance order
    rank = np.empty_like(order)
    rank[order] = np.arange(order.size)
    return p[first[order]], rank[np.asarray(inv).reshape(-1)]


def close_projections(xyz: np.ndarray, tol: float) -> List[Tuple[int, int]]:
    """Pairs (i, j) of points on different levels whose horizontal projections are closer than
    0.1 x the median mesh size but not coincident (EDU-17, D-INC-09)."""
    p = np.asarray(xyz, float).reshape(-1, 3)
    if p.shape[0] < 2:
        return []
    from scipy.spatial import cKDTree
    d, _ = cKDTree(p).query(p, k=2)
    h = float(np.median(d[:, 1]))
    if not h > 0:
        return []
    pairs = cKDTree(p[:, :2]).query_pairs(0.1 * h)
    out = []
    for i, j in sorted(pairs):
        if abs(p[i, 2] - p[j, 2]) > tol and np.linalg.norm(p[i, :2] - p[j, :2]) > tol:
            out.append((int(i), int(j)))
    return out


# ======================================================================================
# The incoherency factors of one HOUSE run
# ======================================================================================
@dataclass
class DirectionStats:
    """Listing summary of one frequency and direction."""

    n: int
    lam1_pct: float
    n90: int
    trace_error: float
    min_raw: float
    n_clipped: int


@dataclass
class FactorRun:
    """Inputs of :func:`compute_factors`."""

    settings: IncoherencySettings
    spec: Optional[CoherencySpec]
    int_xyz: np.ndarray                  # (N, 3) interaction-node coordinates
    fnum: np.ndarray                     # (nF,) frequency numbers
    freq: np.ndarray                     # (nF,) Hz
    sar: Optional[np.ndarray] = None     # (nF, N) multiple-excitation ratios
    merge_tol: float = 1e-9


@quiet_fpe
def compute_factors(run: FactorRun, samples: Sequence[int] = (0,),
                    on_decomposition: Optional[Callable[[int, int, Decomposition], None]] = None
                    ) -> Tuple[np.ndarray, List[List[Optional[DirectionStats]]], Dict[str, object]]:
    """Incoherency factors of the samples ``samples`` (0-based) -> ``s`` (len(samples), nF, 3, N).

    For incoherent motion every frequency q and direction d gets its coherency matrix and its
    decomposition (``on_decomposition(q, d, dec)`` is called once per pair, e.g. for the I N C O
    table); X and Y share one decomposition when their matrices are identical.  Then
    ``s <- s exp(-i w tau)`` (wave passage) and ``s <- s SAR`` (multiple excitation).  Returns
    (s, stats[q][d], info)."""
    st = run.settings
    xy = np.asarray(run.int_xyz, float)[:, :2]
    N = xy.shape[0]
    nF = len(run.fnum)
    S_out = np.ones((len(samples), nF, 3, N), complex)
    stats: List[List[Optional[DirectionStats]]] = [[None, None, None] for _ in range(nF)]
    info: Dict[str, object] = {"n_positions": N, "merged": False, "max_trace_error": 0.0, "min_raw": 0.0}
    if st.coh:
        if run.spec is None:
            raise IncoherencyError("incoherent motion needs a coherency model")
        if st.merge:
            pos, idx = unique_positions(xy, run.merge_tol)
            info.update(n_positions=int(pos.shape[0]), merged=pos.shape[0] < N)
        else:
            pos, idx = xy, np.arange(N)
        builders = [run.spec.matrix_builder(pos, d) for d in range(3)]
        same_xy = run.spec.same_matrix(0, 1)
        modes = select_modes(st.nmodes, pos.shape[0])
        gens = phase_generators(st.hseed, st.vseed, st.nsim) if st.stochastic else None
        R = math.radians(st.randphz)
        for q in range(nF):
            f = float(run.freq[q])
            dec_x = None
            for d in range(3):
                if d == 1 and same_xy and dec_x is not None:
                    dec = dec_x
                else:
                    dec = decompose(builders[d](f), st.sign)
                if d == 0:
                    dec_x = dec
                if on_decomposition is not None:
                    on_decomposition(q, d, dec)
                stats[q][d] = DirectionStats(dec.n, float(100.0 * dec.lam[0] / dec.n) if dec.n else 0.0,
                                             dec.modes_for(), dec.trace_error, dec.min_raw, dec.n_clipped)
                info["max_trace_error"] = max(float(info["max_trace_error"]), dec.trace_error)
                info["min_raw"] = min(float(info["min_raw"]), dec.min_raw)
                for j, r in enumerate(samples):
                    if gens is not None:
                        # every generator draws one phase per mode at every frequency, in frequency order:
                        # sample r's phases do not depend on which samples are computed together
                        theta = gens[d][r].uniform(-R, R, size=dec.n)
                        s = synthesize(dec, modes, theta)
                    else:
                        s = synthesize(dec, modes)
                    S_out[j, q, d] = s[idx]
    if st.wpass:
        tau = arrival_delays(xy, st.appv, st.wang, st.xc, st.yc)
        S_out *= wave_passage_factors(run.freq, tau)[None, :, None, :]
        info["tau_range"] = (float(tau.min()), float(tau.max())) if tau.size else (0.0, 0.0)
    if st.me and run.sar is not None:
        S_out *= np.asarray(run.sar, complex)[None, :, None, :]
    return S_out, stats, info


# ======================================================================================
# FILE77
# ======================================================================================
def write_file77(path: Union[str, Path], fnum, freq, int_node, int_xyz, s: np.ndarray, meta: Dict[str, object],
                 extra: Optional[Dict[str, np.ndarray]] = None, module: str = "HOUSE") -> Path:
    """Write one FILE77 (``s`` (nF, 3, nInt); files.SCHEMAS['FILE77'] plus ``x_int_xyz``)."""
    arrays: Dict[str, np.ndarray] = {"fnum": np.asarray(fnum, np.int64), "freq": np.asarray(freq, float),
                                     "int_node": np.asarray(int_node, np.int64),
                                     "s": np.asarray(s, complex).reshape(len(fnum), 3, -1),
                                     "x_int_xyz": np.asarray(int_xyz, float).reshape(-1, 3)}
    if extra:
        arrays.update(extra)
    probs = validate(Container("FILE77", dict(meta), arrays))
    if probs:                                            # pragma: no cover - programming error
        raise IncoherencyError("FILE77: " + "; ".join(probs))
    return write_container(path, "FILE77", arrays, meta, module=module)


def read_file77(path: Union[str, Path]) -> Container:
    """Read and check a FILE77 (schema, shapes, ascending frequency numbers)."""
    p = Path(path)
    try:
        c = read_container(p, kind="FILE77")
    except (ValueError, OSError) as exc:
        raise IncoherencyError(f"{p.name}: {exc}") from None
    probs = validate(c)
    if probs:
        raise IncoherencyError(f"{p.name} does not satisfy the FILE77 schema: " + "; ".join(probs))
    fn = np.asarray(c["fnum"])
    s = np.asarray(c["s"])
    n = np.asarray(c["int_node"]).size
    if s.shape != (fn.size, 3, n):
        raise IncoherencyError(f"{p.name}: s has shape {s.shape}, expected ({fn.size}, 3, {n})")
    if fn.size > 1 and np.any(np.diff(fn) <= 0):
        raise IncoherencyError(f"{p.name}: frequency numbers are not strictly increasing")
    return c


def node_permutation(c77: Container, int_node: Sequence[int], int_xyz: Optional[np.ndarray] = None,
                     tol: float = 1e-6, name: str = "FILE77") -> np.ndarray:
    """Positions in ``c77`` of the interaction nodes ``int_node`` (ANALYS order): every node must be
    in the file, and when the file stores coordinates (``x_int_xyz``) they must agree within ``tol``
    (a FILE77 of another model or numbering is refused)."""
    ids = np.asarray(c77["int_node"], np.int64)
    pos = {int(n): k for k, n in enumerate(ids)}
    miss = [int(n) for n in int_node if int(n) not in pos]
    if miss:
        raise IncoherencyError(f"{name} has no factors for interaction nodes {miss[:10]}"
                               + (" ..." if len(miss) > 10 else "") + " -- re-run HOUSE (or BUILDFILE77)")
    perm = np.array([pos[int(n)] for n in int_node], np.int64)
    if int_xyz is not None and "x_int_xyz" in c77:
        a = np.asarray(c77["x_int_xyz"], float)[perm]
        b = np.asarray(int_xyz, float).reshape(-1, 3)
        if a.shape == b.shape and a.size:
            dev = float(np.abs(a - b).max())
            if dev > tol:
                k = int(np.argmax(np.abs(a - b).max(axis=1)))
                raise IncoherencyError(f"{name}: interaction node {int(int_node[k])} is at another position than in the "
                                       f"HOUSE model (difference {dev:.3g}) -- re-run HOUSE")
    return perm


def build_file77(inputs: Sequence[Container], names: Sequence[str], order: Optional[Sequence[int]] = None
                 ) -> Tuple[Dict[str, np.ndarray], Dict[str, object], List[str]]:
    """Build_FILE77 (requirements 1.1, 3.4.R BUILDFILE77; spec 05b section 7): one FILE77 covering
    the interaction nodes of several per-level FILE77s.

    The inputs must have the same frequency numbers and step; their node sets must not overlap.  The
    nodes are concatenated in input order, or arranged in ``order`` (the interaction nodes of the
    final FILE4) when given -- nodes of ``order`` missing from every input are an error, input nodes
    not in ``order`` are dropped with a note.  Returns (arrays, meta, notes)."""
    if not inputs:
        raise IncoherencyError("BUILDFILE77 needs at least one input FILE77")
    ref = inputs[0]
    fnum = np.asarray(ref["fnum"], np.int64)
    df = float(ref.meta.get("df", 0.0))
    notes: List[str] = []
    ids: List[int] = []
    s_parts, xyz_parts = [], []
    for c, nm in zip(inputs, names):
        if not np.array_equal(np.asarray(c["fnum"], np.int64), fnum):
            raise IncoherencyError(f"{nm}: frequency numbers differ from {names[0]}")
        dfc = float(c.meta.get("df", df))
        if df > 0 and abs(dfc - df) > 1e-6 * df:
            raise IncoherencyError(f"{nm}: frequency step {dfc:.9g} Hz differs from {names[0]} ({df:.9g} Hz)")
        if c.meta.get("stochastic") != ref.meta.get("stochastic"):
            notes.append(f"{nm}: {'stochastic' if c.meta.get('stochastic') else 'deterministic'} factors combined with "
                         f"{'stochastic' if ref.meta.get('stochastic') else 'deterministic'} ones of {names[0]}")
        if c.meta.get("sim") != ref.meta.get("sim"):
            notes.append(f"{nm}: simulation {c.meta.get('sim')} combined with simulation {ref.meta.get('sim')} of {names[0]}")
        nodes = [int(n) for n in np.asarray(c["int_node"])]
        dup = sorted(set(nodes) & set(ids))
        if dup:
            raise IncoherencyError(f"{nm}: interaction nodes {dup[:10]} are already in an earlier input")
        ids += nodes
        s_parts.append(np.asarray(c["s"], complex))
        xyz_parts.append(np.asarray(c["x_int_xyz"], float) if "x_int_xyz" in c else np.full((len(nodes), 3), np.nan))
    s = np.concatenate(s_parts, axis=2)
    xyz = np.concatenate(xyz_parts, axis=0)
    ids_a = np.asarray(ids, np.int64)
    if order is not None:
        pos = {n: k for k, n in enumerate(ids)}
        order = [int(n) for n in order]
        miss = [n for n in order if n not in pos]
        if miss:
            raise IncoherencyError(f"interaction nodes {miss[:10]} of FILE4 are in none of the input FILE77s")
        extra = sorted(set(ids) - set(order))
        if extra:
            notes.append(f"{len(extra)} input nodes are not interaction nodes of FILE4 and are dropped: {extra[:10]}")
        perm = np.array([pos[n] for n in order], np.int64)
        s, xyz, ids_a = s[:, :, perm], xyz[perm], np.asarray(order, np.int64)
    arrays = {"fnum": fnum, "freq": np.asarray(ref["freq"], float), "int_node": ids_a, "s": s, "x_int_xyz": xyz}
    meta = {k: v for k, v in ref.meta.items() if k not in ("format", "kind", "version", "module")}
    meta.update(method="BUILT", built_from=list(names), df=df,
                note="Build_FILE77: per-level factors combined (spec 05b section 7); the modes of different levels "
                     "were computed separately (sign convention D-INC-03 at every level)")
    return arrays, meta, notes


# ======================================================================================
# Listing helpers
# ======================================================================================
def inco_table(write: Callable[[str], None], q: int, fnum: int, f: float, direction: int, dec: Decomposition,
               max_rows: int = 100) -> None:
    """The ``I N C O`` table of one frequency and direction (manual Eqs. 6.1-6.4; the header contains
    the documented search key ``I N C O``)."""
    write("")
    write(f" I N C O H E R E N T   M O D E   C O N T R I B U T I O N S  (I N C O)   frequency {q + 1} "
          f"(number {fnum}, f = {f:.6g} Hz), direction {DIRECTIONS[direction]}")
    write(f"   N = {dec.n} interaction nodes (plan positions); trace check |sum(lambda) - N| / N = "
          f"{dec.trace_error:.2e} (Eq. 6.1)" + (f"; {dec.n_clipped} negative eigenvalue(s) set to 0 (most negative "
                                                 f"{dec.min_raw * dec.n:.3g})" if dec.n_clipped else ""))
    write(f"{'mode':>9s}{'lambda':>16s}{'upsilon (%)':>14s}{'cumulative (%)':>16s}")
    ups = dec.contributions()
    cum = np.cumsum(ups)
    m = min(dec.n, max_rows)
    for j in range(m):
        write(f"{j + 1:>9d}{dec.lam[j]:>16.6e}{ups[j]:>14.6f}{cum[j]:>16.6f}")
    if dec.n > m:
        write(f"   ... {dec.n - m} further modes: {100.0 - cum[m - 1]:.6f} % of the variance")
