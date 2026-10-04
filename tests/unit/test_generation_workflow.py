"""Structure + excavation workflow of spec 09 section 7 (EXCAV, MERGESOIL, RMVUNUSED, NCOM, INTGEN, FIXROT)
ending with AFWRITE: the HOUSE deck carries the generated interaction nodes (code 0 flags, D-MDL-07)."""
from __future__ import annotations

from sassi.io import decks
from sassi.prep import Interpreter, Kind
from sassi.verify.problems.vp_generation import basement_lines

SOIL = ["L,1,5,0.12,1500,800,0.05,0.05", "L,2,5,0.12,1600,850,0.05,0.05", "L,3,20,0.13,2000,1000,0.02,0.02",
        "TOPL,1,2,3"]
OPTIONS = ["FREQ,1,2,4,8", "SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,4096,1", "POINT,0,3,4.5",
           "HOUSE,32.2,0,0,2,0,0,0,0,0", "AOPT,0,0,0,0,0,1,0,0,0,0,0,0,0,0"]


def test_structure_excavation_workflow_to_afwrite(tmp_path):
    mdir = tmp_path / "ssi"
    mdir.mkdir()
    ui = Interpreter(cwd=tmp_path)
    lines = (["ACTM,1"] + basement_lines(2, 5.0, (-10.0, -5.0, 0.0)) + SOIL
             + ["EXCAV,2", "ACTM,2", "ETYPEGEN,2", "ACTM,3", f"MDL,ssi,{mdir}", "MERGESOIL,1,2,1,,,,ssi_Excv.map",
                "RMVUNUSED", "NCOM", "INTGEN,0", "INTGEN,2", "FIXROT", "EXCSTRCHK"] + OPTIONS + ["AFWRITE"])
    for ln in lines:
        assert ui.execute(ln), (ln, ui.sink.texts(Kind.ERROR))
    m = ui.model
    assert sorted(m.nodes) == list(range(1, 28))                    # 25 structure + 2 soil-only nodes
    flagged = sorted(i for i, n in m.nodes.items() if 0 in n.flags)
    assert len(flagged) == 26                                       # EVBN: all but the interior node at -5
    d = decks.read(decks.deck_path(mdir, "ssi", "HOUSE"), "HOUSE")
    assert sorted(int(r[0]) for r in d.table("interaction").rows) == flagged
    assert (mdir / "ssi_Excv.map").exists()
    assert not any("shared" in w for w in ui.sink.texts(Kind.WARNING))   # EXCSTRCHK clean
