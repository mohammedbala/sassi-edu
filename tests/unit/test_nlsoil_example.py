"""Tutorial example 6 (examples/ex06_nonlinear_soil.pre): a loose backfill block behind an embedded wall,
non-linear soil SSI iterations driven by NLSSIITER.  The example is staged in a temporary copy of the
examples directory and run as ``sassi --cwd examples run`` does (sassi.verify.problems.vp_examples)."""
from __future__ import annotations

import pytest

from sassi.core import nlsoil as NL
from sassi.prep import Interpreter
from sassi.prep.messages import Kind
from sassi.verify.problems.vp_examples import EXAMPLES_DIR, listing_status, purge_binaries, run_example

pytestmark = pytest.mark.skipif(not (EXAMPLES_DIR / "ex06_nonlinear_soil.pre").exists(),
                                reason="the tutorial examples ship with the source tree")


@pytest.fixture(scope="module")
def ex06(tmp_path_factory):
    root = tmp_path_factory.mktemp("ex06")
    ui, summary, pre = run_example("ex06_nonlinear_soil", root)
    md = pre.parent / "ex06"
    yield ui, summary, md
    purge_binaries(md)


def test_ex06_runs_cleanly(ex06):
    ui, summary, md = ex06
    errors = ui.sink.texts(Kind.ERROR)
    assert summary.ok and summary.errors == 0 and not errors, "\n".join(errors[:20])
    status = listing_status(md)
    for mod in ("SOIL", "SITE", "POINT", "HOUSE", "ANALYS", "STRESS", "MOTION"):
        assert "status OK" in status.get(f"ex06_{mod}.out", "missing"), mod
    pin = NL.read_pin(md / "ex06.pin")
    assert (pin.ngrp, pin.esf) == (1, 0.65) and pin.group(3).istr == 1 and pin.group(3).nelem == 8
    # GFAC scales the free field (requirements 4.4 item 7): 0.4 x the native soil ~ the backfill's G_max
    assert (pin.group(3).materials[0].gfac, pin.group(3).materials[0].dfac) == (0.4, 1.0)


def test_ex06_iterations_converge(ex06):
    _, _, md = ex06
    rows = NL.read_convergence(md)
    assert rows and rows[-1]["converged"] == 1 and len(rows) <= NL.MAX_ITERATIONS + 1
    assert all(r["converged"] == 0 for r in rows[:-1])
    f78 = NL.read_nlfile(md / "FILE78", "FILE78")
    f74 = NL.read_nlfile(md / "FILE74", "FILE74")
    assert f78.iteration == f74.iteration == rows[-1]["iteration"] >= 2
    assert f74.direction == "SRSS X+Y+Z"
    st = NL.status(md)
    assert st is not None and st.converged
    # the loose backfill degrades under the strong shaking, more next to the wall (x = 0..2: odd elements)
    ratio = {r.element: r.ratio for r in f78.rows}
    assert all(0.3 < v < 0.8 for v in ratio.values())
    assert all(ratio[e] < ratio[e + 1] for e in (1, 3, 5, 7))
    # the directional strains combine to the FILE74 strains (SRSS)
    fx, fy, fz = (NL.read_nlfile(md / f"FILE74{d}") for d in "XYZ")
    for r, a, b, c in zip(f74.rows, fx.rows, fy.rows, fz.rows):
        assert r.gamma_eff == pytest.approx(float(NL.srss(a.gamma_eff, b.gamma_eff, c.gamma_eff)), rel=1e-9)
        assert a.gamma_eff > c.gamma_eff                               # horizontal X input dominates the vertical


def test_ex06_wall_top_spectrum(ex06):
    from sassi.io import textfiles
    _, _, md = ex06
    zpa = {n: float(textfiles.read_xy(md / f"000{n}TR_X01.RS")[-1, 1]) for n in (22, 24)}
    assert 0.5 < zpa[24] < 0.7                                         # back of the backfill ~ free field 0.58 g
    assert zpa[22] > 1.1 * zpa[24]                                     # the wall top carries the deck mass


def test_ex06_write_inp_round_trip(ex06):
    ui, _, md = ex06
    b = Interpreter(cwd=md)
    b.run_file(str(md / "ex06.pre"), resolve=False)
    assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)[:10]
    assert b.model.same_state(ui.model)
    assert b.model.options.entry("PINGRP", 3) is not None
