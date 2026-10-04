"""Hysteresis models and equivalent-linear properties of Option NON (NONLINEAR module).

What this module does (for the structural engineer)
---------------------------------------------------
Option NON linearises a nonlinear wall panel or spring for the frequency-domain SSI solution
(requirements 4.15, spec 05d section 3.5).  For each nonlinear element and SSI iteration:

1. the deformation history ``x(t)`` (panel shear strain, or spring relative displacement) is run
   through a **hysteresis model** in the time domain, which gives the force history ``F(t)``;
2. the equivalent-linear amplitude is ``x_eq = EDF max|x(t)|`` (EDF ~ 0.8, manual section 6.5.4);
3. the **secant stiffness** of the backbone curve (BBC) at ``x_eq``, normalised by the BBC initial
   slope ``K_el = Y1/X1`` (D-NON-11), scales the elastic modulus: ``E_new = E_el K_sec/K_el``;
4. the **hysteretic damping** ``xi_h = E_D/(4 pi E_S)`` is the area of the stabilised cyclic loop of
   the model at ``x_eq`` (E_D) over 4 pi times the secant strain energy ``E_S = x_eq F_eq/2``;
5. the damping used by the next iteration is ``xi = min(cutoff, scale xi_h + [xi_el])`` with
   scale 0 -> 1 and cutoff 0 -> none (D-NON-04);
6. convergence: max relative change of E < 2 % and max change of xi < 0.5 % (absolute) between
   iterations, at most 10 iterations (D-NON-06).

Backbone curves (BBC, spec 05d section 3.7)
-------------------------------------------
The user gives the points ``(X_i, Y_i)``, i = 1..n, **without the origin**; point 1 is the cracking
point and ``Yield Num.`` is the index of the yield point.  :class:`Backbone` is the odd, piecewise
linear curve through (0, 0), (X_1, Y_1), ..., (X_n, Y_n); beyond X_n the force is held constant
(SASSI-EDU: no failure branch is modelled; the NONLINEAR listing flags elements that pass X_n).

Hysteresis models (BBC Type / Force Opt codes of the manual; D-NON-03)
---------------------------------------------------------------------
* **4 GMR -- General Masing Rule** (:class:`MasingGMR`): virgin loading follows the BBC; an
  unloading or reloading branch from the reversal point (x_r, F_r) is ``F = F_r + 2 F_bb((x - x_r)/2)``
  (Masing 1926, factor 2); extended Masing memory rules (Pyke 1979; Kramer 1996, section 6.4.2): a
  branch that reaches the previous reversal point continues on the branch it left, and a branch that
  reaches the largest past excursion joins the backbone.  Closed symmetric loops, so for an elastic-
  perfectly-plastic BBC ``xi_h = 2 (x - x_y)/(pi x)`` and ``K_sec = F_y/x`` (VP-45).
* **1 CMS -- Cheng-Mertz Shear** (:class:`ChengMertzShear`): the S1 shear hysteresis model of
  Cheng & Mertz (1989) for low-rise RC shear walls, transcribed rule by rule from subroutine HYST04 of
  the INRESB-3D-SUP program (Cheng, F.Y. and Mertz, G.E., "A computer program for inelastic analysis of
  3-dimensional reinforced-concrete and steel seismic buildings", Civil Engineering Study 89-31,
  University of Missouri-Rolla / NSF, 1989, Appendix C, listing pp. 151-153; NTIS PB90-123225).  The
  model uses only the multi-linear backbone (cracking point = BBC point 1, ``SI = PC/DC``), exactly the
  data the manual's BBC holds.  Its rules (numbers as in HYST04 and manual Fig. 1.3):

  - rule 1, loading on the backbone; elastic (``K = SI``) before the first cracking;
  - rules 2/3/4, unloading in three force bands measured from the peak load PM of the current sign
    (above ``PA = PM - PC``, between PA and ``PB = PC/2``, below PB) with the degrading stiffnesses
    ``S1 = SI min(1.4675 (DC/DMAX)^0.345, 1)``, ``S2 = SI min(0.7761 (DC/DMAX)^0.5195, 1)``,
    ``S3 = SI min(0.0707 (DC/DMAX)^1.369, 1)`` (DMAX = largest past displacement in either direction),
    steepened when needed so that the unloading does not pass ``DO``, the zero-force intercept of the
    peak-to-peak line;
  - rule 5, unloading inside small loops toward the stored reversal point of the loop;
  - rules 6/7, reloading toward ``(D2, P2)``: the point at ``0.95 PMAX`` on the first unloading branch
    from the largest past peak (cracked wall) or the previous peak (uncracked);
  - rules 8/9, **pinched** reloading after a load reversal below 0.75 PC with the slip stiffness
    ``SR = SI min((DC/DMAX)^1.02, 1)`` up to PC/4 and the harmonic mean of the slip and the reloading
    stiffness up to 0.75 PC (no pinching when the direct reloading line is steeper: a high cracking
    point gives unpinched loops, manual Fig. 1.3);
  - rule 10, loading from (D2, P2) toward ``(1.04 DMAX, PMAX)`` (degrading factor ALPHA = 1.04 when
    the current direction holds the largest displacement) until the backbone is met;
  - rule 11, reloading inside small loops toward their stored reversal points (up to NI = 10 loops).

  The original is load-incremental (it returns the tangent stiffness and the load limit of the
  current rule); here every rule is a straight branch followed exactly in displacement control, the
  rule changing at its load limit or at a reversal (same rule logic, exact event detection).
* **3 TAK -- Takeda** (:class:`Takeda`): Takeda, Sozen and Nielsen (1970), "Reinforced concrete
  response to simulated earthquakes", J. Struct. Div. ASCE 96(ST12), 2557-2573, with the rules of
  Otani (1974, UILU-ENG-74-2029) as implemented in HYST06 of the same program: trilinear backbone
  through the cracking point (BBC point 1) and the yield point (BBC yield point), flat beyond yield
  (the manual: "the yield point force is used in the Takeda model for the wall panel peak
  capacity"); unloading before yield toward the cracking point of the opposite direction; unloading
  after yield with ``K_u = K_cy (D_y/D_max)^0.4``, ``K_cy = (P_y + P_c)/(D_y + D_c)`` (exponent 0.4,
  D-NON-03); reloading after a zero crossing toward the peak of the opposite direction (or its yield
  point when that line is steeper and the direction has not yielded); inner loops reload toward the
  most recent unloading point first (Otani's points U0, U1 ... as a stack).  The manual: TAK is usable
  only with the elastic damping (its hysteretic damping is not used, see :func:`combine_damping`).
* **2 CMB -- Cheng-Mertz Bending**: "not included in this version" (manual); asking for it is an error.

Every model is driven by displacements and records the exact vertices of its force-deformation path
(``model.vertices``), so loop areas are exact polygon areas.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

#: BBC Type / Force Opt codes (manual section 6.5.4; D-NON-03)
MODEL_CODES: Dict[int, str] = {1: "CMS", 2: "CMB", 3: "TAK", 4: "GMR"}
MODEL_NAMES: Dict[int, str] = {1: "Cheng-Mertz Shear (CMS)", 2: "Cheng-Mertz Bending (CMB)",
                               3: "Takeda (TAK)", 4: "General Masing Rule (GMR)"}
#: convergence of the iterations (D-NON-06): relative change of E, absolute change of the damping ratio
TOL_E = 0.02
TOL_XI = 0.005
MAX_ITERATIONS = 10
#: number of cycles of the cyclic test that gives the stabilised loop (the last cycle is used)
LOOP_CYCLES = 3


class HysteresisError(ValueError):
    """Invalid backbone curve or model request."""


def fsign(a: float, b: float) -> float:
    """Fortran ``SIGN(A, B)``: |A| with the sign of B (B = 0 counts as positive)."""
    return abs(a) if b >= 0.0 else -abs(a)


def kmin(dd: float, dp: float, si: float) -> float:
    """``KMIN`` of INRESB-3D-SUP: the slope ``dp/dd`` when positive, limited to the initial stiffness SI;
    SI when the two differences have opposite signs (or one is zero)."""
    if dd * dp > 0.0:
        return min(dp / dd, si)
    return si


# ======================================================================================
# Backbone curve
# ======================================================================================
@dataclass
class Backbone:
    """Backbone curve: points (x_i, y_i) > 0 without the origin; ``yield_index`` is 1-based.

    The curve is odd (``F(-x) = -F(x)``), linear between the points and from the origin to point 1,
    and constant beyond the last point.
    """

    x: np.ndarray
    y: np.ndarray
    yield_index: int = 1

    def __post_init__(self) -> None:
        self.x = np.asarray(self.x, dtype=float).ravel()
        self.y = np.asarray(self.y, dtype=float).ravel()
        self.yield_index = int(self.yield_index)
        problems = backbone_problems(self.x, self.y, self.yield_index)
        if problems:
            raise HysteresisError("; ".join(problems))
        self._xp = np.concatenate(([0.0], self.x))
        self._fp = np.concatenate(([0.0], self.y))
        # cumulative area under the curve at the breakpoints (exact integral of the polyline)
        self._area = np.concatenate(([0.0], np.cumsum(np.diff(self._xp) * (self._fp[1:] + self._fp[:-1]) / 2.0)))

    # ---- named points --------------------------------------------------------------------
    @property
    def n(self) -> int:
        return len(self.x)

    @property
    def k_el(self) -> float:
        """Elastic stiffness: the slope of the first segment ``Y1/X1`` (D-NON-11)."""
        return float(self.y[0] / self.x[0])

    @property
    def x_cr(self) -> float:
        return float(self.x[0])

    @property
    def f_cr(self) -> float:
        return float(self.y[0])

    @property
    def x_y(self) -> float:
        return float(self.x[self.yield_index - 1])

    @property
    def f_y(self) -> float:
        return float(self.y[self.yield_index - 1])

    # ---- evaluation ----------------------------------------------------------------------
    def force(self, d):
        """F(d), odd, piecewise linear, constant beyond the last point (scalar or array)."""
        a = np.abs(np.asarray(d, dtype=float))
        f = np.interp(a, self._xp, self._fp)
        out = np.sign(np.asarray(d, dtype=float)) * f
        return float(out) if np.ndim(out) == 0 else out

    def integral(self, a: float) -> float:
        """``int_0^a F(x) dx`` for a >= 0 (exact for the polyline)."""
        a = float(abs(a))
        k = int(np.searchsorted(self._xp, a, side="right")) - 1
        if k >= len(self._xp) - 1:
            return float(self._area[-1] + self._fp[-1] * (a - self._xp[-1]))
        f_a = float(np.interp(a, self._xp, self._fp))
        return float(self._area[k] + (a - self._xp[k]) * (self._fp[k] + f_a) / 2.0)

    def secant(self, a: float) -> float:
        """Secant stiffness ``F(a)/a`` (``K_el`` for a -> 0)."""
        a = abs(float(a))
        if a <= self.x[0]:
            return self.k_el
        return float(self.force(a)) / a

    def breakpoints(self) -> np.ndarray:
        return self.x.copy()

    def scaled(self, fx: float = 1.0, fy: float = 1.0) -> "Backbone":
        return Backbone(self.x * fx, self.y * fy, self.yield_index)

    def max_slope_ratio(self) -> float:
        """Largest segment slope after the first one divided by the first slope (should be <= 1)."""
        if self.n < 2:
            return 0.0
        s = np.diff(self.y) / np.diff(self.x)
        return float(np.max(s) / self.k_el)


def backbone_problems(x: Sequence[float], y: Sequence[float], yield_index: int, strict: bool = False) -> List[str]:
    """Why points (x, y) and the yield index do not form a usable BBC (spec 11 section 2.2 checks).

    ``strict`` (CMS and TAK) also requires strictly increasing forces: their rules need a positive
    tangent on every segment; GMR accepts flat segments (elastic-perfectly-plastic curves)."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    out: List[str] = []
    if x.size != y.size:
        out.append(f"{x.size} X values but {y.size} Y values")
        return out
    if x.size < 2:
        out.append(f"{x.size} point(s): at least 2 points are needed (point 1 = cracking point)")
        return out
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        out.append("non-finite values")
        return out
    if np.any(x <= 0):
        out.append("X values must be > 0 (the origin is implied and not given)")
    if np.any(np.diff(x) <= 0):
        out.append("X values must be strictly increasing")
    if np.any(y <= 0):
        out.append("Y values must be > 0")
    if np.any(np.diff(y) < 0):
        out.append("Y values must not decrease (softening branches are not supported by the hysteresis rules)")
    elif strict and np.any(np.diff(y) <= 0):
        out.append("Y values must increase strictly for the CMS and TAK models (use GMR for a flat plateau)")
    if not 1 <= int(yield_index) <= x.size:
        out.append(f"yield point number {yield_index} outside 1..{x.size}")
    return out


# ======================================================================================
# Base class
# ======================================================================================
class HysteresisModel:
    """Displacement-driven hysteresis model; ``step(x)`` moves to displacement x and returns F."""

    code = 0
    #: maximum number of rule changes inside one displacement step (safety net)
    MAX_EVENTS = 400

    def __init__(self, backbone: Backbone):
        self.bb = backbone
        self.reset()

    def reset(self) -> None:
        self.x = 0.0
        self.f = 0.0
        self.vertices: List[Tuple[float, float]] = [(0.0, 0.0)]
        self.warnings: List[str] = []
        self.beyond_bbc = False

    def _vertex(self, x: float, f: float) -> None:
        if self.vertices and self.vertices[-1][0] == x and self.vertices[-1][1] == f:
            return
        self.vertices.append((float(x), float(f)))
        if abs(x) > self.bb.x[-1] * (1.0 + 1e-12):
            self.beyond_bbc = True

    def step(self, x_new: float) -> float:   # pragma: no cover - abstract
        raise NotImplementedError

    def run(self, xs: Sequence[float]) -> np.ndarray:
        """Force history for the displacement history ``xs`` (continuing from the current state)."""
        out = np.empty(len(xs))
        for i, v in enumerate(np.asarray(xs, dtype=float)):
            out[i] = self.step(float(v))
        return out

    def backbone_force(self, d: float) -> float:
        """The model's own envelope: the BBC (GMR, CMS) or the trilinear Takeda curve (TAK)."""
        return float(self.bb.force(d))

    def backbone_secant(self, a: float) -> float:
        a = abs(float(a))
        if a <= self.bb.x_cr:
            return self.bb.k_el
        return self.backbone_force(a) / a


# ======================================================================================
# General Masing Rule (GMR)
# ======================================================================================
class MasingGMR(HysteresisModel):
    """General Masing Rule with the extended (memory) rules; see the module docstring.

    Branch through the base point (x_b, F_b) with factor c: ``F = F_b + c F_bb((x - x_b)/c)``;
    c = 1 for the backbone (base = origin), c = 2 for every unloading/reloading branch."""

    code = 4

    def reset(self) -> None:
        super().reset()
        self.stack: List[Tuple[float, float]] = []      # reversal points, oldest first
        self.base = (0.0, 0.0, 1.0)
        self.d = 0

    def _branch(self, x: float) -> float:
        xb, fb, c = self.base
        return fb + c * float(self.bb.force((x - xb) / c))

    def _closure(self) -> Optional[Tuple[float, float]]:
        if len(self.stack) >= 2:
            return self.stack[-2]
        if len(self.stack) == 1:
            x1, f1 = self.stack[0]
            return (-x1, -f1)          # the branch from the largest excursion rejoins the backbone
        return None

    def _close(self) -> None:
        if len(self.stack) >= 2:
            self.stack.pop()
            self.stack.pop()
        else:
            self.stack.pop()
        if self.stack:
            xb, fb = self.stack[-1]
            self.base = (xb, fb, 2.0)
        else:
            self.base = (0.0, 0.0, 1.0)

    def _breaks(self, x0: float, x1: float) -> None:
        """Record the kinks of the current branch strictly between x0 and x1 (exact loop polygons)."""
        xb, fb, c = self.base
        pts = np.concatenate((xb + c * self.bb.x, xb - c * self.bb.x))
        lo, hi = (x0, x1) if x0 < x1 else (x1, x0)
        sel = np.sort(pts[(pts > lo) & (pts < hi)])
        if x1 < x0:
            sel = sel[::-1]
        for p in sel:
            self._vertex(float(p), self._branch(float(p)))

    def step(self, x_new: float) -> float:
        x_new = float(x_new)
        x = self.x
        if x_new == x:
            return self.f
        d = 1 if x_new > x else -1
        if self.d != 0 and d != self.d:
            self.stack.append((x, self.f))
            self.base = (x, self.f, 2.0)
        self.d = d
        for _ in range(self.MAX_EVENTS):
            cl = self._closure()
            if cl is not None and (x_new - cl[0]) * d >= 0.0:
                self._breaks(x, cl[0])
                x, self.f = cl
                self._close()
                self._vertex(x, self.f)
                if x == x_new:
                    break
                continue
            self._breaks(x, x_new)
            self.f = self._branch(x_new)
            x = x_new
            self._vertex(x, self.f)
            break
        self.x = x
        return self.f


# ======================================================================================
# Cheng-Mertz Shear (CMS): HYST04 of INRESB-3D-SUP
# ======================================================================================
class ChengMertzShear(HysteresisModel):
    """Cheng-Mertz S1 shear hysteresis model (rules 1-11 of HYST04, see the module docstring)."""

    code = 1
    ALPHA = 1.04          # degrading stiffness factor (DATA ALPHA/1.04/)
    NI = 10               # small-amplitude loop reversal points stored (manual input NI, "usually 10")

    def __init__(self, backbone: Backbone, ni: int = NI):
        probs = backbone_problems(backbone.x, backbone.y, backbone.yield_index, strict=True)
        if probs:
            raise HysteresisError("CMS: " + "; ".join(probs))
        self.ni = int(ni)
        self.DS = np.asarray(backbone.x, dtype=float)
        self.PS = np.asarray(backbone.y, dtype=float)
        # beyond the last BBC point HYST04 lets the wall fail (K = 0, rule < 0); SASSI-EDU holds the force
        # (a practically flat extension, slope 1e-6 SI) so that the equivalent linearisation continues
        self.PC, self.DC = float(self.PS[0]), float(self.DS[0])
        self.SI = self.PC / self.DC
        xe = self.DS[-1] * 1.0e3
        self.DS = np.append(self.DS, xe)
        self.PS = np.append(self.PS, self.PS[-1] + 1.0e-6 * self.SI * (xe - self.DS[-2]))
        super().__init__(backbone)

    def reset(self) -> None:
        super().reset()
        self.rule = 0.0
        self.K, self.A = self.SI, self.PC
        self.DIR = 0
        self.active = False
        self.PMG, self.PMH, self.DMG, self.DMH = self.PC, -self.PC, self.DC, -self.DC
        self.BACKB = False
        self.PR: List[float] = []
        self.DR: List[float] = []
        self.FR: List[str] = []
        self.SR1 = self.SR2 = self.SR3 = self.SRM = 0.0
        self.failures = 0
        # direction cracked by the motion: HYST04 tests CRACKD = |DM| > DC, which the load-incremental
        # original reaches by the small overshoot of the step that crosses the cracking point; in exact
        # displacement control the flag is set when the motion reaches |D| = DC
        self.cracked = {1: False, -1: False}

    # ---- statement functions of HYST04 ------------------------------------------------------
    def FS1(self, X: float) -> float:
        return self.SI * min(1.4675 * (self.DC / X) ** 0.345, 1.0)

    def FS2(self, X: float) -> float:
        return self.SI * min(0.7761 * (self.DC / X) ** 0.5195, 1.0)

    def FS3(self, X: float) -> float:
        return self.SI * min(0.0707 * (self.DC / X) ** 1.369, 1.0)

    def FSR(self, X: float) -> float:
        return self.SI * min(abs(self.DC / X) ** 1.02, 1.0)

    # ---- helpers -----------------------------------------------------------------------------
    @staticmethod
    def _direction(P: float, s: int) -> int:
        """DIR of HYST04: 1 positive loading, 2 positive unloading, 3 negative loading, 4 negative unloading."""
        if P > 0.0:
            return 1 if s > 0 else 2
        if P < 0.0:
            return 3 if s < 0 else 4
        return 1 if s > 0 else 3

    @staticmethod
    def _passed(F: float, A: float, DIR: int) -> bool:
        if DIR in (1, 4):
            return F > A
        return F < A

    def _maxima(self, P: float, D: float) -> None:
        self.PMG = max(self.PMG, self.PC, P)
        self.DMG = max(self.DMG, self.DC, D)
        self.PMH = min(self.PMH, -self.PC, P)
        self.DMH = min(self.DMH, -self.DC, D)
        if abs(D) >= self.DC * (1.0 - 1e-12):
            self.cracked[1 if D > 0 else -1] = True

    def step(self, x_new: float) -> float:
        x_new = float(x_new)
        x, P = self.x, self.f
        if x_new == x:
            return P
        s = 1 if x_new > x else -1
        for _ in range(self.MAX_EVENTS):
            DIR = self._direction(P, s)
            if (not self.active) or DIR != self.DIR:
                DIRL = self.DIR
                self._evaluate(P, x, DIR, DIRL, s)
                self.DIR = DIR
                self.active = True
            K, A = self.K, self.A
            Ft = P + K * (x_new - x)
            if self._passed(Ft, A, DIR) and K > 0.0:
                xa = x + (A - P) / K
                if (xa - x) * s < 0.0:            # limit behind the current point (round-off): no move
                    xa = x
                if (x_new - xa) * s <= 0.0:       # the limit is reached exactly at the end of the step
                    x, P = x_new, A
                    self._maxima(P, x)
                    self._vertex(x, P)
                    self.active = False
                    break
                x, P = xa, A
                self._maxima(P, x)
                self._vertex(x, P)
                self.active = False
                continue
            x, P = x_new, Ft
            self._maxima(P, x)
            self._vertex(x, P)
            break
        else:      # pragma: no cover - safety net
            self.warnings.append("CMS: too many rule changes in one step; the step was completed on the last rule")
            P = P + self.K * (x_new - x)
            x = x_new
            self._maxima(P, x)
            self._vertex(x, P)
        self.x, self.f = x, P
        return P

    # ---- rule logic (HYST04 from "DETERMINE IF ELASTIC" to label 9999) -------------------------
    def _evaluate(self, P: float, D: float, DIR: int, DIRL: int, s: int) -> None:
        PC, DC, SI = self.PC, self.DC, self.SI
        AP = abs(P)
        # elastic before the first cracking (RULE 0): K = SI in both directions up to +-PC
        if self.rule == 0.0 and AP < PC * (1.0 - 1e-12):
            self.K = SI
            self.A = PC if s > 0 else -PC
            return
        LRULE = int(self.rule)
        DMAX = max(self.DMG, -self.DMH)
        PMAX = max(self.PMG, -self.PMH)
        if DIR in (1, 2):
            PM, DM = self.PMG, self.DMG
        else:
            PM, DM = self.PMH, self.DMH
        # a signed reference for SIGN(., P) when P is exactly 0 (zero crossing): the side being loaded
        Pref = P if P != 0.0 else (1e-300 if DIR in (1, 4) else -1e-300)
        # ---- erase small loop flags (loops that the response has passed)
        for i in range(len(self.PR)):
            if ((self.FR[i] == "L" and (P <= self.PR[i] or D <= self.DR[i]))
                    or (self.FR[i] == "7" and (P >= self.PR[i] or D >= self.DR[i]))):
                del self.PR[i:], self.DR[i:], self.FR[i:]
                break
        IR = len(self.PR)
        if DIR in (1, 3):
            self._loading(P, D, DIR, DIRL, LRULE, PM, DM, PMAX, DMAX, AP, Pref, IR)
        else:
            self._unloading(P, D, DIR, PM, DM, AP, IR)
        # ---- label 100: save reversal points for small loops
        IRULE = int(self.rule)
        REVRSL = IRULE in (1, 10) or LRULE in (1, 10)
        if (DIRL == 1 and DIR == 2 and not REVRSL) or (DIRL == 4 and DIR == 3):
            if len(self.PR) < self.ni:
                self.PR.append(P), self.DR.append(D), self.FR.append("7")
        elif (DIRL == 3 and DIR == 4 and not REVRSL) or (DIRL == 2 and DIR == 1):
            if len(self.PR) < self.ni:
                self.PR.append(P), self.DR.append(D), self.FR.append("L")
        if not (1.05 * SI > self.K > 0.0):
            # HYST04: K = 0 and RULE = -RULE (wall inactive).  SASSI-EDU keeps a tiny stiffness so the
            # history can continue; the listing reports the count.
            self.failures += 1
            self.K = 1.0e-9 * SI
            self.A = fsign(1.0e30, s)

    def _rule1(self, P: float, D: float, AP: float, Pref: float) -> None:
        """Rule 1: loading on the backbone (steepest line <= SI to a backbone point ahead)."""
        K, A = 0.0, None
        for J in range(1, len(self.DS)):
            DD = self.DS[J] - abs(D)
            DP = self.PS[J] - AP
            if AP < self.PS[J] and DD * DP > 0.0:
                S = DP / DD
                if K < S <= self.SI:
                    K = S
                    A = fsign(self.PS[J], Pref)
        self.K = K
        self.A = A if A is not None else fsign(self.PS[-1], Pref)
        self.rule = 1.0
        self.BACKB = True
        del self.PR[:], self.DR[:], self.FR[:]

    def _loading(self, P, D, DIR, DIRL, LRULE, PM, DM, PMAX, DMAX, AP, Pref, IR) -> None:
        PC, DC, SI = self.PC, self.DC, self.SI
        PC2 = fsign(PC / 2.0, PM)
        PR, DR, FR = self.PR, self.DR, self.FR
        label = 30
        I = IR + 1
        CRACKD = False
        REVRSL = False
        P2 = D2 = S = 0.0
        for _ in range(100):
            if label == 30:
                if self.BACKB or self.rule == 0.0:
                    self._rule1(P, D, AP, Pref)
                    return
                CRACKD = abs(DM) > DC or self.cracked[1 if DIR == 1 else -1]
                REVRSL = ((DIR == 1 and DIRL == 4) or LRULE == 8 or (DIR == 3 and DIRL == 2) or LRULE == 9
                          or LRULE == 5 or LRULE == 11)
                if CRACKD:
                    P2 = 0.95 * fsign(PMAX, PM)
                    D2 = fsign(DMAX, PM) - 0.05 * fsign(PMAX, PM) / self.FS1(DMAX)      # FD2(DMAX)
                else:
                    P2, D2 = PM, DM
                S = kmin(D2 - D, P2 - P, SI)
                label = 20
                continue
            # ---- label 20: rule 11 (reloading inside small loops) and the rule chain 6 / 7 / 8-9 / 10
            I -= 1
            if I > 0:
                i = I - 1
                if not ((FR[i] == "7" and DIR == 1) or (FR[i] == "L" and DIR == 3)):
                    continue
                if self.SR1 <= 0 or self.SR2 <= 0 or self.SR3 <= 0 or self.SRM <= 0:
                    I = 1
                    continue
                if AP < abs(PC2) and S > self.SRM:
                    REVRSL = False
                    I = 1
                    continue
                if abs(PR[i]) <= PC / 4.0:
                    self.K = kmin(D - DR[i], P - PR[i], SI)
                    self.A = PR[i]
                    self.rule = 11.1 + I / 1000.0
                elif abs(PR[i]) <= 3.0 * PC / 4.0:
                    X2 = D2 - (P2 - 1.5 * PC2) / self.SR3
                    X1 = X2 - PC2 / self.SR2
                    X = D + (PC2 / 2.0 - P) / self.SR2
                    if AP < PC / 4.0 and ((DIR == 1 and X < X1) or (DIR == 3 and X > X1)):
                        self.K = kmin(X1 - D, PC2 / 2.0 - P, SI)
                        self.A = fsign(PC / 4.0, PM)
                        self.rule = 11.2 + I / 1000.0
                    else:
                        self.K = kmin(DR[i] - D, PR[i] - P, SI)
                        self.A = PR[i]
                        self.rule = 11.3 + I / 1000.0
                else:
                    X2 = DR[i] - (PR[i] - 1.5 * PC2) / self.SRM
                    X1 = X2 - PC2 / self.SR2
                    X = D + (PC2 / 2.0 - P) / self.SR2
                    if AP < PC / 4.0 and ((DIR == 1 and X < X1) or (DIR == 3 and X > X1)):
                        self.K = kmin(X1 - D, PC2 / 2.0 - P, SI)
                        self.A = fsign(PC / 4.0, PM)
                        self.rule = 11.4 + I / 1000.0
                    else:
                        X = D + (1.5 * PC2 - P) / self.SR3
                        if AP < 3.0 * PC / 4.0 and ((DIR == 1 and X < X2) or (DIR == 3 and X > X2)):
                            self.K = kmin(X2 - D, 1.5 * PC2 - P, SI)
                            self.A = fsign(3.0 * PC / 4.0, PM)
                            self.rule = 11.5 + I / 1000.0
                        else:
                            self.SR3 = self.SRM
                            self.K = kmin(DR[i] - D, PR[i] - P, SI)
                            self.A = PR[i]
                            self.rule = 11.6 + I / 1000.0
                return
            if AP < abs(PC2) and not REVRSL:                                   # rule 6
                self.K = self.FS1(DMAX)
                self.A = PC2
                self.rule = 6.0
                return
            if AP < abs(P2) and ((not REVRSL) or AP >= PC * 0.75):            # rule 7
                self.K = S
                self.A = P2
                self.rule = 7.0
                return
            if AP < PC * 0.75:                                                 # rules 8 and 9
                SR = self.FSR(DMAX)
                S1, S2, S3 = self.FS1(DMAX), self.FS2(DMAX), self.FS3(DMAX)
                PA = fsign(PMAX, PM) - fsign(PC, PM)
                DA = fsign(DMAX, PM) - fsign(PC, PM) / S1
                PB = fsign(min(abs(PA), PC / 2.0), PM)
                DB = DA - (PA - PB) / S2
                DOP = DB - PB / S3
                DO = DM - PM * (self.DMG - self.DMH) / (self.PMG - self.PMH)
                if DIR == 1:
                    DOP = max(DOP, DO)
                if DIR == 3:
                    DOP = min(DOP, DO)
                DC2 = DOP + PC2 / S1
                self.SRM = kmin(DC2 - D2, PC2 - P2, S1)
                if S > self.SRM:
                    REVRSL = False
                    label = 20
                    continue
                SRP = kmin(DC2 - D, PC2 - P, S1)
                self.SR1 = max(SRP, min(SR, S))
                DC2P = D + (PC2 - P) / self.SR1
                self.SR3 = kmin(D2 - DC2P, P2 - PC2, SI)
                self.SR2 = 2.0 / (1.0 / self.SR1 + 1.0 / self.SR3)
                if AP < PC / 4.0:
                    self.K = self.SR1
                    self.A = fsign(PC / 4.0, PM)
                    self.rule = 8.0
                else:
                    self.K = self.SR2
                    self.A = fsign(PC, PM) * 0.75
                    self.rule = 9.0
                return
            # rule 10: loading towards the backbone curve
            self.BACKB = True
            if not CRACKD:
                label = 30
                continue
            del self.PR[:], self.DR[:], self.FR[:]
            ALPHAP = self.ALPHA if abs(DM) == abs(DMAX) else 1.0
            K = kmin(ALPHAP * fsign(DMAX, D2) - D2, fsign(PMAX, P2) - P2, SI)
            for J in range(1, len(self.DS)):
                if self.PS[J] <= abs(PMAX):
                    continue
                J0 = J - 1
                DD = self.DS[J] - self.DS[J0]
                DP = self.PS[J] - self.PS[J0]
                if DP * DD <= 0.0:
                    continue
                DX, PX = _intersect(D, P, K, fsign(self.DS[J], Pref), fsign(self.PS[J], Pref), DP / DD)
                # HYST04 tests |PX| > |PMAX| strictly; when the rule-10 target (DMAX, PMAX) lies on the
                # backbone (ALPHAP = 1) the intersection is that point and the test is decided by the
                # overshoot of the original's load steps: a relative tolerance makes it deterministic
                tol = 1e-9
                if (abs(DX) >= self.DS[J0] * (1 - tol) and abs(DX) <= self.DS[J] * (1 + tol)
                        and abs(PX) > abs(PMAX) * (1 - tol) and PM * PX > 0.0):
                    if AP >= abs(PX):
                        continue
                    self.K = kmin(DX - D, PX - P, SI)
                    self.A = PX
                    self.rule = 10.0
                    return
            label = 30
        # safety net (never reached in practice): load on the backbone
        self._rule1(P, D, AP, Pref)   # pragma: no cover

    def _unloading(self, P, D, DIR, PM, DM, AP, IR) -> None:
        PC, SI = self.PC, self.SI
        DMAX = max(self.DMG, -self.DMH)
        self.BACKB = False
        S1, S2, S3 = self.FS1(DMAX), self.FS2(DMAX), self.FS3(DMAX)
        DO = DM - PM * (self.DMG - self.DMH) / (self.PMG - self.PMH)
        PA = PM - fsign(PC, PM)
        DA = DM - fsign(PC, PM) / S1
        PB = fsign(min(abs(PA), PC / 2.0), PM)
        SO2 = kmin(DA - DO, PA, SI)
        DB = DA - (PA - PB) / max(S2, SO2)
        SO3 = kmin(DB - DO, PB, SI)
        SU = 0.0
        # rule 5: unloading inside small loops
        I = IR + 1
        while True:
            I -= 1
            if I <= 0:
                break
            i = I - 1
            if not ((DIR == 2 and self.FR[i] == "L") or (DIR == 4 and self.FR[i] == "7")):
                continue
            pr, dr = self.PR[i], self.DR[i]
            same_band = ((DIR == 2 and ((P > PA and pr > PA) or (PA >= P > PB and pr > PB) or (PB >= P and pr > 0.0)))
                         or (DIR == 4 and ((P < PA and pr < PA) or (PA <= P < PB and pr < PB)
                                           or (PB <= P and pr < 0.0))))
            if same_band:
                if dr == D:
                    continue
                K = (pr - P) / (dr - D)
                if K <= 0.0 or K > SI:
                    continue
                self.K, self.A, self.rule = K, pr, 5.0 + I / 1000.0
                return
            SU = kmin(dr - D, pr - P, SI)
            break
        if AP > abs(PA):                                     # rule 2: top segment
            self.K, self.A, self.rule = max(S1, SU), PA, 2.0
        elif AP > abs(PB):                                   # rule 3: middle segment
            self.K, self.A, self.rule = max(S2, SO2, SU), PB, 3.0
        else:                                                # rule 4: bottom segment
            self.K, self.A, self.rule = max(S3, SO3, SU), 0.0, 4.0


def _intersect(x1: float, y1: float, s1: float, x2: float, y2: float, s2: float) -> Tuple[float, float]:
    """``INTSCT``: intersection of the lines through (x1, y1) with slope s1 and (x2, y2) with slope s2."""
    if s1 == s2:
        return 0.0, 0.0
    b1 = y1 - s1 * x1
    b2 = y2 - s2 * x2
    x = (b1 - b2) / (s2 - s1)
    return x, (b1 * s2 - b2 * s1) / (s2 - s1)


# ======================================================================================
# Takeda (TAK)
# ======================================================================================
class Takeda(HysteresisModel):
    """Takeda model (Takeda et al. 1970; Otani 1974 rules of HYST06), see the module docstring.

    States: ``elastic`` (before the first cracking, ``K0 = Pc/Dc`` both ways), ``primary`` (on the
    trilinear envelope, moving outward), ``unload`` (straight line with the unloading stiffness toward
    zero force) and ``reload`` (straight line toward the next target point: the most recent unloading
    point of an interrupted reloading first, then the peak ``UM`` of the direction)."""

    code = 3
    EXPONENT = 0.4

    def __init__(self, backbone: Backbone):
        probs = backbone_problems(backbone.x, backbone.y, backbone.yield_index, strict=True)
        if probs:
            raise HysteresisError("TAK: " + "; ".join(probs))
        self.Dc, self.Pc = float(backbone.x[0]), float(backbone.y[0])
        self.Dy, self.Py = float(backbone.x_y), float(backbone.f_y)
        if self.Dy <= self.Dc:            # yield at the cracking point: bilinear Takeda
            self.Dy, self.Py = self.Dc, self.Pc
        self.K0 = self.Pc / self.Dc
        self.K1 = (self.Py - self.Pc) / (self.Dy - self.Dc) if self.Dy > self.Dc else 0.0
        self.Kcy = (self.Py + self.Pc) / (self.Dy + self.Dc)
        super().__init__(backbone)

    def reset(self) -> None:
        super().reset()
        self.mode = "elastic"
        self.um = {1: (self.Dc, self.Pc), -1: (-self.Dc, -self.Pc)}     # peaks UM (cracking points at first)
        self.cracked = {1: False, -1: False}
        self.yielded = {1: False, -1: False}
        self.ku = self.K0
        self.side = 1
        self.line = (0.0, 0.0, self.K0)        # (x_s, f_s, slope) of the current straight branch
        self.limit = 0.0                        # force at the end of the current straight branch
        self.targets: List[Tuple[float, float]] = []
        self.d = 0

    # ---- envelope ------------------------------------------------------------------------------
    def backbone_force(self, d: float) -> float:
        a = abs(float(d))
        if a <= self.Dc:
            f = self.K0 * a
        elif a <= self.Dy:
            f = self.Pc + self.K1 * (a - self.Dc)
        else:
            f = self.Py
        return f if d >= 0 else -f

    def _corner_ahead(self, x: float, s: int) -> Optional[float]:
        """Next kink of the envelope beyond x in direction s (cracking or yield displacement)."""
        a = abs(x)
        if a < self.Dc * (1.0 - 1e-12):
            return s * self.Dc
        if a < self.Dy * (1.0 - 1e-12):
            return s * self.Dy
        return None

    def _unloading_stiffness(self, x: float, f: float, s: int) -> float:
        if not self.yielded[s]:
            return (f + s * self.Pc) / (x + s * self.Dc)          # toward the opposite cracking point
        dmax = max(abs(x), abs(self.um[1][0]), abs(self.um[-1][0]), self.Dy)
        return self.Kcy * (self.Dy / dmax) ** self.EXPONENT

    def _mark(self, x: float, f: float) -> None:
        """Bookkeeping of a point on the envelope: cracking, yielding and the peak UM of its side."""
        s = 1 if x > 0 else -1
        a = abs(x)
        if a >= self.Dc * (1.0 - 1e-12):
            self.cracked[s] = True
            if a >= abs(self.um[s][0]):
                self.um[s] = (x, f)
        if a >= self.Dy * (1.0 - 1e-12):
            self.yielded[s] = True

    # ---- straight branches -------------------------------------------------------------------------
    def _line_to(self, x: float, f: float, xt: float, ft: float) -> None:
        k = (ft - f) / (xt - x) if xt != x else self.K0
        if k <= 0.0:
            k = self.ku
        self.line = (x, f, k)
        self.limit = ft

    def _start_reload(self, x: float, f: float) -> None:
        while self.targets:
            xt, ft = self.targets[-1]
            if (xt - x) * self.side > 0.0 and (ft - f) * self.side > 0.0:
                self.mode = "reload"
                self._line_to(x, f, xt, ft)
                return
            self.targets.pop()                 # already passed: drop it
        self.mode = "primary"

    def _zero_crossing(self, x: float, s_new: int) -> None:
        """Unloading reached zero force at x; the response continues into direction s_new."""
        self.side = s_new
        self.targets = []
        if not self.cracked[s_new]:
            # the opposite direction never cracked: continue the unloading line to its cracking force,
            # then toward its yield point and on the envelope
            self.mode = "reload"
            self.line = (x, 0.0, self.ku)
            self.limit = s_new * self.Pc
            self.targets = [(s_new * self.Dy, s_new * self.Py), (float("nan"), s_new * self.Pc)]
            return
        xm, fm = self.um[s_new]
        if not self.yielded[s_new]:
            k_um = fm / (xm - x)
            k_y = s_new * self.Py / (s_new * self.Dy - x)
            if abs(fm) < self.Py and k_y > k_um:          # aim at the yield point (HYST06 rule 4.33)
                self.um[s_new] = (s_new * self.Dy, s_new * self.Py)
        self.targets = [self.um[s_new]]
        self._start_reload(x, 0.0)

    def _target_reached(self, x: float, f: float) -> None:
        if self.targets:
            self.targets.pop()
        s = self.side
        if abs(f) >= self.Pc * (1.0 - 1e-12):
            self.cracked[s] = True
        fb = self.backbone_force(x)
        if abs(f - fb) <= 1e-9 * max(abs(fb), self.Pc) and abs(x) >= self.Dc * (1.0 - 1e-12):
            # on the envelope (the peak UM, or the cracking point of a virgin direction): primary curve
            self.targets = []
            self.mode = "primary"
            self._mark(x, f)
            return
        if self.targets:
            self._start_reload(x, f)
            return
        if abs(x) < self.Dy:
            # a reload that ended off the envelope (cracking force reached on a degraded line): go on to
            # the yield point first
            self.targets = [(s * self.Dy, s * self.Py)]
            self._start_reload(x, f)
            return
        self.mode = "primary"
        self._mark(x, f)

    # ---- driver ----------------------------------------------------------------------------------------
    def _reverse(self, x: float, f: float, d: int) -> None:
        """Displacement reversal at (x, f); the new motion direction is d."""
        if self.mode == "elastic":
            return
        if self.mode == "unload":
            if f == 0.0:
                self._zero_crossing(x, d)
                return
            # partial unloading reversed: reload toward the start of the unloading branch first
            xs, fs, _ = self.line
            self.side = 1 if f > 0 else -1
            if not self.targets or self.targets[-1] != (xs, fs):
                self.targets.append((xs, fs))
            self._start_reload(x, f)
            return
        s = 1 if f > 0 else (-1 if f < 0 else -d)
        if self.mode == "reload":
            if f == 0.0:
                self._zero_crossing(x, d)
                return
            # interrupted reloading: the point becomes the first target of the next reloading of this side
            self.targets.append((x, f))
        else:                                   # primary
            self._mark(x, f)
            self.targets = []
            self.ku = self._unloading_stiffness(x, f, s)
        self.side = s
        self.mode = "unload"
        self.line = (x, f, self.ku)
        self.limit = 0.0

    def step(self, x_new: float) -> float:
        x_new = float(x_new)
        x, f = self.x, self.f
        if x_new == x:
            return f
        d = 1 if x_new > x else -1
        if self.d != 0 and d != self.d:
            self._reverse(x, f, d)
        self.d = d
        for _ in range(self.MAX_EVENTS):
            if self.mode == "elastic":
                xc = d * self.Dc
                if (x_new - xc) * d > 0.0:
                    x, f = xc, d * self.Pc
                    self._vertex(x, f)
                    self.mode, self.side = "primary", d
                    self._mark(x, f)
                    continue
                x, f = x_new, self.K0 * x_new
                self._vertex(x, f)
                break
            if self.mode == "primary":
                xc = self._corner_ahead(x, d)
                if xc is not None and (x_new - xc) * d > 0.0:
                    x, f = xc, self.backbone_force(xc)
                    self._vertex(x, f)
                    self._mark(x, f)
                    continue
                x, f = x_new, self.backbone_force(x_new)
                self._vertex(x, f)
                self._mark(x, f)
                break
            xs, fs, k = self.line
            ft = fs + k * (x_new - xs)
            lim = self.limit
            if k > 0.0 and (ft - lim) * d >= 0.0:
                xl = xs + (lim - fs) / k
                if (xl - x) * d < 0.0:
                    xl = x
                x, f = xl, lim
                self._vertex(x, f)
                if self.mode == "unload":
                    self._zero_crossing(x, d)
                else:
                    self._target_reached(x, f)
                if x == x_new:
                    break
                continue
            x, f = x_new, ft
            self._vertex(x, f)
            break
        self.x, self.f = x, f
        return f


# ======================================================================================
# Factory, loops and equivalent-linear properties
# ======================================================================================
def make_model(code: int, backbone: Backbone) -> HysteresisModel:
    """Hysteresis model of BBC Type / Force Opt ``code`` (1 CMS, 3 TAK, 4 GMR)."""
    code = int(code)
    if code == 4:
        return MasingGMR(backbone)
    if code == 1:
        return ChengMertzShear(backbone)
    if code == 3:
        return Takeda(backbone)
    if code == 2:
        raise HysteresisError("the Cheng-Mertz Bending (CMB) model is not included in this version (manual)")
    raise HysteresisError(f"unknown hysteresis model code {code} (1 CMS, 3 TAK, 4 GMR)")


@dataclass
class LoopProperties:
    """Stabilised symmetric cyclic loop of a model at amplitude ``amplitude``."""
    amplitude: float
    f_pos: float            # force at +amplitude at the end of the last cycle
    f_neg: float            # force at -amplitude in the last cycle
    energy: float           # dissipated energy E_D of the last cycle (loop area)
    closed: float           # |F(end) - F(start)| of the last cycle (0 for a closed loop)
    vertices: np.ndarray    # (m, 2) path of the last cycle

    @property
    def f_amp(self) -> float:
        return 0.5 * (self.f_pos - self.f_neg)

    @property
    def secant(self) -> float:
        return self.f_amp / self.amplitude if self.amplitude > 0 else 0.0

    @property
    def damping(self) -> float:
        """``xi_h = E_D/(4 pi E_S)``, E_S = amplitude F_amp / 2."""
        es = 0.5 * self.amplitude * self.f_amp
        return self.energy / (4.0 * math.pi * es) if es > 0 else 0.0


def cyclic_loop(model: HysteresisModel, amplitude: float, cycles: int = LOOP_CYCLES) -> LoopProperties:
    """Cycle ``model`` (reset first) between +-amplitude ``cycles`` times; the last cycle is measured.

    The path is 0 -> +a -> -a -> +a ...; each leg is one displacement step (the models detect every
    rule change exactly), so the vertices are the exact polygon of the loop."""
    a = float(abs(amplitude))
    model.reset()
    model.step(a)
    for _ in range(max(1, int(cycles)) - 1):
        model.step(-a)
        model.step(a)
    i0 = len(model.vertices) - 1
    f_start = model.f
    f_neg = model.step(-a)
    f_pos = model.step(a)
    path = np.asarray(model.vertices[i0:], dtype=float)
    xs, fs = path[:, 0], path[:, 1]
    energy = float(np.sum(np.diff(xs) * (fs[1:] + fs[:-1]) / 2.0))
    return LoopProperties(a, f_pos, f_neg, -energy if energy < 0 else energy, abs(f_pos - f_start), path)


def masing_loop_energy(backbone: Backbone, amplitude: float) -> float:
    """Closed-form Masing loop area ``E_D = 8 int_0^a F dx - 4 a F(a)`` (independent check of GMR)."""
    a = abs(float(amplitude))
    return 8.0 * backbone.integral(a) - 4.0 * a * float(backbone.force(a))


def hysteretic_damping(code: int, backbone: Backbone, amplitude: float) -> float:
    """``xi_h(amplitude)`` from the stabilised cyclic loop of model ``code`` (0 below cracking)."""
    a = abs(float(amplitude))
    if a <= backbone.x_cr * (1.0 + 1e-12):
        return 0.0
    loop = cyclic_loop(make_model(code, backbone), a)
    return max(0.0, loop.damping)


def secant_ratio(code: int, backbone: Backbone, amplitude: float, model: Optional[HysteresisModel] = None) -> float:
    """``K_sec/K_el`` at ``amplitude`` on the model's backbone (the BBC; trilinear for TAK)."""
    m = model if model is not None else make_model(code, backbone)
    return m.backbone_secant(amplitude) / backbone.k_el


def eql_curve(code: int, backbone: Backbone, amplitudes: Sequence[float]) -> np.ndarray:
    """The ``.crv`` table: rows ``(amplitude, E/E_el, xi_h, K_sec, F)`` (manual section 6.5.4)."""
    model = make_model(code, backbone)
    rows = []
    for a in amplitudes:
        a = abs(float(a))
        r = secant_ratio(code, backbone, a, model)
        xi = hysteretic_damping(code, backbone, a)
        rows.append((a, r, xi, r * backbone.k_el, model.backbone_force(a)))
    return np.asarray(rows, dtype=float)


def crv_amplitudes(backbone: Backbone, n: int = 60) -> np.ndarray:
    """Amplitude grid of the ``.crv`` curves: log spaced from X1/10 to 2 X_n, plus the BBC points."""
    lo, hi = backbone.x_cr / 10.0, backbone.x[-1] * 2.0
    a = np.geomspace(lo, hi, max(int(n), 2))
    return np.unique(np.concatenate((a, backbone.x)))


def combine_damping(xi_h: float, xi_el: float, cutoff_pct: float = 0.0, scale: float = 0.0,
                    include_elastic: bool = False) -> float:
    """``xi = min(cutoff, scale xi_h + [xi_el])`` with scale 0 -> 1 and cutoff 0 -> none (D-NON-04).

    ``cutoff_pct`` is the dialog's *Damping Cutoff %* (percent); damping ratios are fractions."""
    s = float(scale) if float(scale) != 0.0 else 1.0
    xi = s * float(xi_h) + (float(xi_el) if include_elastic else 0.0)
    if float(cutoff_pct) > 0.0:
        xi = min(xi, float(cutoff_pct) / 100.0)
    return xi


@dataclass
class EquivalentLinear:
    """Equivalent-linear properties of one nonlinear element for the next SSI iteration."""
    x_max: float
    x_eq: float
    ratio: float            # E_new / E_el = K_sec(x_eq)/K_el
    xi_h: float             # hysteretic damping at x_eq
    xi: float               # damping of the next iteration (combine_damping)
    ductility: float        # max|x| / x_cr
    f_mu: float             # K_el max|x| / |F(t*)|, t* = time of max|x|
    force: np.ndarray       # nonlinear force history F(t)
    beyond_bbc: bool        # the history passed the last BBC point
    notes: List[str]


def equivalent_linear(code: int, backbone: Backbone, history: Sequence[float], edf: float, xi_el: float,
                      cutoff_pct: float = 0.0, scale: float = 0.0, include_elastic: bool = False) -> EquivalentLinear:
    """Run the hysteresis model on ``history`` and linearise it at ``x_eq = edf max|x|`` (requirements 4.15).

    TAK: the manual allows the Takeda model only with the constant elastic damping, so ``xi = xi_el``."""
    x = np.asarray(history, dtype=float)
    model = make_model(code, backbone)
    F = model.run(x) if x.size else np.zeros(0)
    notes: List[str] = list(model.warnings)
    if getattr(model, "failures", 0):
        notes.append(f"{model.failures} rule evaluation(s) without a valid stiffness (HYST04 failure state)")
    x_max = float(np.max(np.abs(x))) if x.size else 0.0
    x_eq = float(edf) * x_max
    ratio = secant_ratio(code, backbone, x_eq, model)
    xi_h = hysteretic_damping(code, backbone, x_eq)
    if int(code) == 3:
        xi = float(xi_el)
        notes.append("TAK: equivalent damping = the constant elastic damping (manual warning)")
    else:
        xi = combine_damping(xi_h, xi_el, cutoff_pct, scale, include_elastic)
    k = int(np.argmax(np.abs(x))) if x.size else 0
    f_star = abs(float(F[k])) if x.size else 0.0
    f_mu = (backbone.k_el * x_max / f_star) if f_star > 0 else 1.0
    return EquivalentLinear(x_max, x_eq, ratio, xi_h, xi, x_max / backbone.x_cr, f_mu, F, model.beyond_bbc, notes)


@dataclass
class Convergence:
    """Change of the properties between two iterations (D-NON-06)."""
    max_de: float          # max |E_new - E_old| / E_old
    max_dxi: float         # max |xi_new - xi_old| (damping ratio)
    converged: bool

    def text(self) -> str:
        return (f"max |dE/E| = {100 * self.max_de:.3f} % (tolerance {100 * TOL_E:g} %), max |d xi| = "
                f"{100 * self.max_dxi:.3f} % (tolerance {100 * TOL_XI:g} %): "
                + ("converged" if self.converged else "not converged"))


def convergence(e_old: Sequence[float], e_new: Sequence[float], xi_old: Sequence[float], xi_new: Sequence[float],
                tol_e: float = TOL_E, tol_xi: float = TOL_XI) -> Convergence:
    """Max relative change of E and max absolute change of the damping ratio (D-NON-06)."""
    eo, en = np.asarray(e_old, float), np.asarray(e_new, float)
    xo, xn = np.asarray(xi_old, float), np.asarray(xi_new, float)
    de = float(np.max(np.abs(en - eo) / np.abs(eo))) if eo.size else 0.0
    dx = float(np.max(np.abs(xn - xo))) if xo.size else 0.0
    return Convergence(de, dx, de < tol_e and dx < tol_xi)
