"""Tutorial example 7 (examples/ex07_option_non.pre): Option NON on a two-storey RC shear-wall building --
eight wall panels with BBCGEN backbones, X/Y/Z input, COMB_XYZ_THD and the NONLINITER iterations.  The example
is staged in a temporary copy of the examples directory and run as ``sassi --cwd examples run`` does."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.io import decks
from sassi.modules import nonlinear as NL
from sassi.prep import Interpreter
from sassi.prep.messages import Kind
from sassi.verify.problems.vp_examples import EXAMPLES_DIR, listing_status, purge_binaries, run_example

pytestmark = [pytest.mark.skipif(not (EXAMPLES_DIR / "ex07_option_non.pre").exists(),
                                 reason="the tutorial examples ship with the source tree"), pytest.mark.slow]

STOREY1, STOREY2 = (1, 2, 3, 4), (5, 6, 7, 8)


@pytest.fixture(scope="module")
def ex07(tmp_path_factory):
    root = tmp_path_factory.mktemp("ex07")
    ui, summary, pre = run_example("ex07_option_non", root)
    md = pre.parent / "ex07"
    yield ui, summary, md
    purge_binaries(md)


def test_ex07_runs_cleanly(ex07):
    ui, summary, md = ex07
    errors = ui.sink.texts(Kind.ERROR)
    assert summary.ok and summary.errors == 0 and not errors, "\n".join(errors[:20])
    status = listing_status(md)
    for mod in ("SITE", "POINT", "HOUSE", "ANALYS", "MOTION", "RELDISP", "NONLINEAR"):
        assert "status OK" in status.get(f"ex07_{mod}.out", "missing"), mod
    # only the expected warnings: FCOPY replacing the HOUSE deck in every iteration
    others = [w for w in ui.sink.texts(Kind.WARNING) if "overwritten" not in w]
    assert not others, others


def test_ex07_panels_and_backbones(ex07):
    ui, _, md = ex07
    m = ui.model
    from sassi.prep.commands import nonlinear_cmds as NC
    assert [k for k, _ in NC.panels(m)] == list(range(1, 9))
    for k in range(1, 9):
        cv = NC.curve(m, k)
        assert len(cv["x"]) == 22 and cv["yield_"] == 21 and cv["type"] == 1
        assert cv["y"][0] == pytest.approx(0.3 * cv["y"][20])                  # CrackingForceLevel 0.3
        assert cv["y"][0] / cv["x"][0] == pytest.approx(3.0e7 / 2.4 * 0.3 * 12.0)   # G A_W
    inp = (md / "COMB_XYZ_THD.inp").read_text().split("\n")
    assert inp[0] == "36"                                                     # 12 corner nodes x 3 DOFs
    assert (md / "ex07_NONLINBAT.pre").exists()


def test_ex07_iterations_converge(ex07):
    _, _, md = ex07
    rows = NL.read_convergence(md)
    assert rows[0]["iteration"] == 0 and rows[-1]["converged"] == 1
    assert all(r["converged"] == 0 for r in rows[:-1])
    assert 3 <= len(rows) - 1 <= 10                                           # D-NON-06: at most 10 iterations
    assert rows[-1]["max_de_pct"] < 2.0 and rows[-1]["max_dxi_pct"] < 0.5
    # the changes decrease after the first iterations (no divergence)
    de = [r["max_de_pct"] for r in rows]
    assert de[-1] < de[2] < de[0]


def test_ex07_cracked_first_storey(ex07):
    _, _, md = ex07
    it, rows = NL.read_props(md / "Panel_EQL_Matl_Prop.txt")
    assert it == len(NL.read_convergence(md))
    for k in STOREY1:                                                         # cracked walls
        r = rows[k]
        assert 0.2 < r.ratio < 0.8 and r.xi_new > r.xi_el and r.x_max > 3.0e-4 * 0.5
    for k in STOREY2:                                                         # still elastic
        assert rows[k].ratio == 1.0 and rows[k].xi_new == pytest.approx(0.04)
    fmu = np.loadtxt(md / "Panel.fmu")
    assert np.all(fmu[np.isin(fmu[:, 0], STOREY1), 3] > 2.0)                  # ductility (cracking based)
    assert np.all(fmu[np.isin(fmu[:, 0], STOREY1), 6] > 1.5)                  # force-reduction factor F_mu
    # the next HOUSE deck holds exactly these properties
    new = decks.read(md / "ex07_new.hou", "HOUSE")
    mats = {int(r["id"]): r for r in new.rows("materials")}
    for k in range(1, 9):
        assert mats[10 + k]["val1"] == pytest.approx(rows[k].k_new, rel=1e-12)
        assert mats[10 + k]["sdamp"] == pytest.approx(rows[k].xi_new, rel=1e-12)
    # hysteresis loops: strain and force histories of panel 1, and the elastic results kept by NONLINSAVE
    g, dt = NL.read_history(md / "Panel0001.thd")
    f, _ = NL.read_history(md / "Panel0001.ths")
    assert dt == 0.005 and len(g) == len(f) and np.max(np.abs(f)) > 0
    assert (md / "Panel0001_elastic.thd").exists() and (md / "Panel_EQL_Matl_Prop_elastic.txt").exists()


def test_ex07_write_inp_round_trip(ex07):
    ui, _, md = ex07
    b = Interpreter(cwd=md)
    b.run_file(str(md / "ex07.pre"), resolve=False)
    assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)[:10]
    assert b.model.same_state(ui.model)
