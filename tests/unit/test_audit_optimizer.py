"""Final audit, finding F-02 (docs/verification/AUDIT_REPORT.md): HOUSE node optimizer and post-processing
node requests (requirements 2.6, "HOUSE optimizer used: MOTION/RELDISP node requests are in the new
numbering -- warning").

With Optimize Model, FILE4/FILE8 and the MOTION/RELDISP files use the new node numbers and the manual
has the user select the post-processing nodes from the ``.hounew``.  Before the audit nothing told the user
that a request in the model numbering names another node: example 1 with HOUSEX,1 reported the mat nodes
37 and 78-81 for the requested 41 and 82-85 (the stick) with "0 warning(s)".  MOTION and RELDISP now list
every requested node whose model number differs (``new = model``) and requests absent from the map.

Finding F-04: Option NON (NONLINEAR) reads the ``.THD`` histories of its panel / spring nodes by model number,
so after a renumbering it would silently use other nodes' histories; it now stops with a clear message.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from sassi.core.renumber import read_map
from sassi.io import decks
from sassi.verify import builders as B

FS = B.FrequencySet.fourier(0.01, 1024, [4, 20, 41, 82])


def _optimized_stick(wd):
    site = B.layered_site([(2.0, 200.0, 400.0, 2.0, 0.05)], (400.0, 800.0, 2.1, 0.02))
    stick = B.Stick(heights=[3.0, 6.0], masses=[50.0, 50.0], E=3.0e7, A=1.0, I=0.1, beta=0.05)
    mdl = B.stick_on_mat(site, stick, half_width=2.0, ndiv=2)
    B.run_soil(wd, "m", site, FS, layer=0, rad=mdl.rad)
    d = mdl.deck("m")
    d["optimize"] = 1
    B.write_deck(wd, "m", d)
    B.run("HOUSE", wd, "m")
    B.write_deck(wd, "m", B.analys_deck(FS))
    B.run("ANALYS", wd, "m")
    return mdl, read_map(wd / "m.map")


def test_motion_and_reldisp_list_the_model_numbers_of_renumbered_requests(tmp_path):
    mdl, old_new = _optimized_stick(tmp_path)
    top = mdl["top"]                                   # a request in the model numbering ...
    assert old_new[top] != top                         # ... names another node in the new numbering
    new_old = {n: o for o, n in old_new.items()}
    B.write_history(tmp_path / "acc.th", 0.1 * np.sin(np.linspace(0, 40, 400)), 0.01)
    md = B.motion_deck(FS, thfile="acc.th", nout=[(top, 1, 1, 1, 0, 0, 0, 1)], dur=4.0, cplx=1)
    B.write_deck(tmp_path, "m", md)
    B.run("MOTION", tmp_path, "m")
    out = B.listing(tmp_path, "m", "MOTION")
    assert "HOUSE node optimizer (m.map)" in out and f"{top} = {new_old[top]}" in out, out[-3000:]
    # RELDISP: the same note for the RDND node and the RELFILE reference
    rd = decks.new("RELDISP")
    for k, v in dict(model="m", title="m", delt=FS.delt, nft=FS.nft, df=FS.df, thfile="acc.th", mult=1.0,
                     gravity=9.81, relfile=f"{top:05d}TR_X.TFI").items():
        rd[k] = v
    rd.table("rdnd").append([top, 1, 0, 0, 0, 0, 0])
    B.write_deck(tmp_path, "m", rd)
    B.run("RELDISP", tmp_path, "m", check=False)          # (whether the reference exists depends on the map)
    out = B.listing(tmp_path, "m", "RELDISP")
    assert "HOUSE node optimizer (m.map)" in out and f"{top} = {new_old[top]}" in out, out[-3000:]


def test_no_note_without_the_optimizer(tmp_path):
    site = B.layered_site([(2.0, 200.0, 400.0, 2.0, 0.05)], (400.0, 800.0, 2.1, 0.02))
    mdl = B.surface_rigid_mat(site, half_width=2.0, ndiv=2)
    B.run_soil(tmp_path, "m", site, FS, layer=0, rad=mdl.rad)
    B.run_house(tmp_path, "m", mdl)
    B.write_deck(tmp_path, "m", B.analys_deck(FS))
    B.run("ANALYS", tmp_path, "m")
    B.write_history(tmp_path / "acc.th", 0.1 * np.sin(np.linspace(0, 40, 400)), 0.01)
    B.write_deck(tmp_path, "m", B.motion_deck(FS, thfile="acc.th", nout=[(mdl["centre"], 1, 1, 1, 0, 0, 0, 1)],
                                              dur=4.0))
    B.run("MOTION", tmp_path, "m")
    assert "HOUSE node optimizer" not in B.listing(tmp_path, "m", "MOTION")


def test_option_non_refuses_a_renumbered_model(tmp_path):
    """F-04: NONLINEAR reads the .THD histories of its panel / spring nodes by model number; after the HOUSE
    optimizer renumbered the model (non-identity map) those files hold other nodes' histories, so the run
    stops with a clear message instead of computing equivalent-linear properties from the wrong nodes.  An
    identity map (the optimizer found nothing to improve, as for example 7) is accepted."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_nonlinear_module import spring_case
    from sassi.modules import nonlinear as NL
    spring_case(tmp_path)
    (tmp_path / "m.map").write_text("1 1\n2 2\n", encoding="utf-8")            # identity: accepted
    assert NL.run_nonlinear("m", tmp_path) == 0
    for p in tmp_path.glob("SPRING*"):
        p.unlink()
    (tmp_path / "m.map").write_text("1 2\n2 1\n", encoding="utf-8")            # renumbered: refused
    assert NL.run_nonlinear("m", tmp_path) == 1
    out = (tmp_path / "m_NONLINEAR.out").read_text(encoding="utf-8")
    assert "HOUSE node optimizer renumbered 2 nodes" in out and "HOUSEX <optimize> = 0" in out
