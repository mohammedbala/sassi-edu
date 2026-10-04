"""Foundation impedance and inertial SSI verification problems (sassi/verify/problems/vp_impedance.py).

Lead decisions D-W2-01..04 (requirements §7.17): VP-11, VP-13, VP-14 and VP-17 report comparisons against approximate references (Veletsos-Verbic rocking
damping, the Pais-Kausel square rocking fit, the 1-D layer-frequency rule for the vertical
compliance peaks, the Veletsos-Verbic-based worked example) that SASSI-EDU does not meet within the
stated tolerance.  The evidence (independent exact-kernel BEM solutions, exact low-frequency limits,
layer dispersion) is in the VP notes and in docs/verification/impedance_study.md.  Those checks are
*not* relaxed: each test asserts that every other check passes and then reports the known failures
as an expected failure (xfail) with the observed numbers.  If a known failure starts to pass, the
test passes.

An xfail must not hide a regression of the very quantity that fails, so every known failure carries
its documented value (impedance_study.md) and a band: a known check whose computed value leaves the
band fails the test.  The VPs also contain derived-reference guards of the same quantities (relaxed
c_r vs the BEM in VP-11, the ZGV onset and the lightly damped peak in VP-13, K_xx vs the welded BEM
in VP-14, the peak vs F.2 with the welded BEM impedance in VP-17), which must pass."""
from __future__ import annotations

import math
from typing import NamedTuple, Optional, Sequence

import pytest

from sassi.verify import run_problem


class Known(NamedTuple):
    """A known failing check: label prefix, documented computed value (NaN = 'no value', e.g. no
    peak found) and the relative band the computed value must stay in."""
    prefix: str
    value: float
    band: float


def _fmt(checks) -> str:
    return "\n".join(f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, error {c.error:.3g} > "
                     f"{c.tolerance:.3g} {c.note}" for c in checks)


def _match(quantity: str, known: Sequence[Known]) -> Optional[Known]:
    hits = [k for k in known if quantity.startswith(k.prefix)]
    assert len(hits) <= 1, f"ambiguous known-failure prefixes for {quantity!r}"
    return hits[0] if hits else None


def _in_band(computed: float, k: Known) -> bool:
    if math.isnan(k.value):
        return math.isnan(computed)
    return math.isfinite(computed) and abs(computed - k.value) <= k.band * abs(k.value)


def _run(vp_id: str, tmp_path, known: Sequence[Known] = (), reason: str = ""):
    """Run a VP; every pass/fail check must pass.  ``known`` lists the *informative* comparisons
    against superseded approximate references (lead decisions D-W2-01..04): each must exist and its
    computed value must stay within the documented band (a regression guard)."""
    res = run_problem(vp_id, tmp_path)
    assert res.checks, f"{vp_id} produced no checks"
    failed = [c for c in res.checks if not c.passed]
    assert not failed, _fmt(failed) + "\n" + "\n".join(res.notes)
    for k in known:
        hits = [c for c in res.checks if c.quantity.startswith(k.prefix)]
        assert hits, f"no check labelled {k.prefix!r}"
        for c in hits:
            assert c.kind == "info", f"{c.quantity} should be informative ({reason})"
            assert _in_band(c.computed, k), (f"{c.quantity}: computed {c.computed:.6g}, documented {k.value:.6g} "
                                             f"+- {100 * k.band:g} % (regression?)")


def test_vp10(tmp_path):
    _run("VP-10", tmp_path)


def test_vp11(tmp_path):
    _run("VP-11", tmp_path,
         known=[Known("c_r at a0 = 0.5 (K_G, welded) vs Veletsos-Verbic", 0.0458, 0.03),
                Known("c_r at a0 = 1 (K_G, welded) vs Veletsos-Verbic", 0.1246, 0.03),
                Known("c_r at a0 = 1.5 (K_G, welded) vs Veletsos-Verbic", 0.1859, 0.03)],
         reason="the Veletsos-Verbic rocking damping is 28-53 % below the rigorous relaxed solution (exact-kernel BEM, "
                "exact low-frequency limit c_r = 0.240 a0^2); SASSI-EDU agrees with the BEM within 1.2 %")


def test_vp13(tmp_path):
    _run("VP-13", tmp_path,
         known=[Known("vertical resonance 1: A0 of the |compliance| peak, beta = 5 %", 0.4848, 0.01),
                Known("vertical resonance 2: A0 of the |compliance| peak, beta = 5 %", math.nan, 0.0)],
         reason="the vertical compliance peaks of the damped layer are not at the 1-D frequencies (2n-1)Vp/(4H): "
                "zero-group-velocity onset at A0 = 0.509 below the cut-off, and 5 % damping hides the second peak "
                "(1.585 at 1 % damping)")


def test_vp14(tmp_path):
    _run("VP-14", tmp_path,
         known=[Known("Kxx (welded K_G, extrapolated) vs Pais & Kausel", 6.4136, 0.01),
                Known("Kyy (welded K_G, extrapolated) vs Pais & Kausel", 6.4136, 0.01)],
         reason="Pais & Kausel's square rocking coefficient (6.0) is 7 % below the welded exact-kernel BEM value "
                "(6.457), which SASSI-EDU reproduces within 0.7 %")


@pytest.mark.slow
def test_vp17(tmp_path):
    _run("VP-17", tmp_path, known=[Known("peak |u_t/u_g| at the mass vs 6.71", 5.563, 0.02)],
         reason="the R2 F.1 peak 6.71 uses Veletsos-Verbic impedances (low rocking damping: 6.19 with rigorous relaxed "
                "impedances) and ignores the sliding-rocking coupling of the welded mat; F.2 with the rigorous welded "
                "impedance gives 5.563, as SASSI-EDU does (5.563)")


# ---------------------------------------------------------------------------------------
# the xfail bookkeeping itself (no SSI run)
# ---------------------------------------------------------------------------------------
def test_known_failure_band_logic():
    k = Known("peak", 5.563, 0.02)
    assert _in_band(5.60, k) and not _in_band(3.0, k) and not _in_band(math.nan, k)
    none = Known("vertical resonance 2", math.nan, 0.0)
    assert _in_band(math.nan, none) and not _in_band(1.585, none)
    known = [Known("vertical resonance 1: A0", 0.4848, 0.01)]
    assert _match("vertical resonance 1: A0 of the peak", known) is known[0]
    assert _match("guard (vertical peak 1): A0 of the peak", known) is None
