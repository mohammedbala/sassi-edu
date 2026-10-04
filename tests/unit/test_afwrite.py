"""AFWRITE: gating (UT-05, D-AFW-01), shared variables (UT-06, spec 07 section 3), resolved deck contents
(ARCHITECTURE section 3, D-AFW-06, D-HOU-01, D-CNV-10, D-MOT-10) and the frequency grid (UT-16)."""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from sassi.io import decks
from sassi.prep import Interpreter, Kind
from sassi.prep.afwrite import AfwriteError, afwrite

MODEL = """
MDL,demo,{dir}
TIT,AFWRITE demo, embedded block
LOC,1,0,100,0,0,0,0,0
CSYS,1
N,1,0,0,0
N,2,10,0,0
N,3,10,10,0
N,4,0,10,0
N,5,0,0,-5
N,6,10,0,-5
N,7,10,10,-5
N,8,0,10,-5
CSYS,0
N,10,100,0,5
N,11,100,0,10
N,12,105,0,5
INT,1,8,1,1
M,1,4.32e5,0.25,0.15,0.05,0.05
R,1,1,0.8,0.8,0.2,0.1,0.1
L,1,5,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
TOPL,1,1
GROUP,1,SOLID
E,1,5,6,7,8,1,2,3,4
GROUP,2,BEAMS
E,1,10,11,12
GROUP,3,SHELL
GROUP,4,SOLID
ETYPE,1,1,1,1
FREQ,1,2,4,8,16
SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,4096,1
POINT,0,1,4.5
HOUSE,32.2,0,0,2,0,0,0,0,0
ANALYS,0,0,0,0,1,0,30,1,2,3,0
THFILE,H1.acc
DAMP,0.02,0.05
MOTION,0,0,0,0,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1
NOUT,1,1,1,0,0,1,1,1-4
NOUT,1,0,0,1,0,0,0,2-3
NOUT,3,1,1,0,0,1,1,1
STRESS,0,0,1,0,1
EOUT,1,1,1,1,1,1,0,0,0,0,0,0,1,1
RDND,1,1,0,0,0,0,0
FORCE,0
F,10,1,0,0,0.1
F,11,0,2,0
MT,11,32.2,32.2,32.2
"""


def _setup(tmp_path, extra="", aopt="AOPT,1,1,0,1,1,1,0,1,1,1,1,1,1,0"):
    mdir = tmp_path / "demo"
    mdir.mkdir(exist_ok=True)
    (mdir / "H1.acc").write_text("0.005\n" + "\n".join(f"{0.1 * math.sin(k / 5)}" for k in range(300)) + "\n")
    (mdir / "H1.rsi").write_text("1 0.2\n5 0.5\n10 0.4\n")
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(MODEL.format(dir=mdir) + "\n" + aopt + "\n" + extra)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    return ui, mdir


def _deck(mdir, module):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")          # every schema parameter must be present
        return decks.read(decks.deck_path(mdir, "demo", module), module)


def test_afwrite_needs_model_name_and_path(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.execute("N,1,0,0,0")
    with pytest.raises(AfwriteError):
        afwrite(ui.model)
    ui.execute("AFWRITE")
    assert any("MDL" in e for e in ui.sink.texts(Kind.ERROR))


def test_ut05_gating_and_stale_deck(tmp_path):
    ui, mdir = _setup(tmp_path, aopt="AOPT,0,0,0,1,1,1,0,0,1,0,1,0,0,0")
    ui.execute("AFWRITE")
    assert (mdir / "demo.sit").exists() and (mdir / "demo.err").exists()
    ui.execute("SITE,0,1,0,3,2,1,0,1,2048,1,0,0.005,4096,1")          # <nl> = 3: Error 47
    ui.sink.clear()
    ui.execute("AFWRITE")
    assert not (mdir / "demo.sit").exists()
    assert (mdir / "demo.sit.bak").exists()
    for ext in (".poi", ".hou", ".anl", ".mot"):
        assert (mdir / f"demo{ext}").exists()
    warns = ui.sink.texts(Kind.WARNING)
    assert any("demo.sit renamed demo.sit.bak" in w for w in warns)
    assert any(w.startswith("AFWRITE: SITE: input deck not written") for w in warns)
    assert "Error 47 : Illegal Number of Layers for Halfspace Simulation" in (mdir / "demo.err").read_text()


def test_ut06_shared_variables_one_storage_location(tmp_path):
    ui, mdir = _setup(tmp_path, extra="SOILX,0,0,0.1\nSPRO,1,1,S\nSPRO,2,2,S\nDYNP,1,0.001,1,0.001,1,S\n"
                                      "SOIL,100,32.2,0,1,1,8,0.65,1,0\nEQUAKE,0,3,11975,0.05,20,0,1\n"
                                      "RSIN,1,H1.rsi\nRSOUT,1,o.rso\nACCOUT,1,o.acc")
    ui.execute("SITE,0,1,0,20,2,1,0,1,2048,1,0,0.01,4096,1")
    ui.execute("GRAVITY,9.81")
    ui.execute("AFWRITE")
    for mod in ("SOIL", "MOTION", "STRESS", "RELDISP", "EQUAKE", "SITE", "FORCE", "ANALYS", "HOUSE"):
        assert _deck(mdir, mod)["delt"] == 0.01, mod
    house = ui.model.options.record("HOUSE")
    assert house.to_tokens() == ["9.81", "0", "0", "2", "0", "0", "0", "0", "0"]   # only argument 1 changed
    for mod in ("SITE", "HOUSE", "FORCE", "ANALYS", "MOTION", "STRESS", "RELDISP", "EQUAKE"):
        assert _deck(mdir, mod)["gravity"] == 9.81, mod
    assert _deck(mdir, "SOIL")["grav"] == 32.2                                # SOIL keeps its own gravity
    for mod in ("MOTION", "STRESS", "RELDISP"):
        d = _deck(mdir, mod)
        assert (d["thfile"], d["mult"], d["rec1"], d["type"], d["ang"]) == ("H1.acc", 1.0, 1, 0, 30.0)
        assert d["df"] == pytest.approx(1 / (0.01 * 4096), rel=1e-15)


def test_house_deck_resolved_data(tmp_path):
    ui, mdir = _setup(tmp_path)
    ui.execute("AFWRITE")
    d = _deck(mdir, "HOUSE")
    nodes = {r["id"]: r for r in d.rows("nodes")}
    assert nodes[1]["x"] == 100.0                                          # LOC system 1 resolved
    assert [nodes[9][k] for k in ("fx", "fy", "fz", "frx", "fry", "frz")] == [1] * 6   # gap node (W1)
    assert [nodes[12][k] for k in ("fx", "frz")] == [1, 1]                 # K-only node fixed (D-CHK-08)
    assert [nodes[10][k] for k in ("fx", "frz")] == [0, 0]
    assert [r["id"] for r in d.rows("interaction")] == list(range(1, 9))
    assert [r["id"] for r in d.rows("groups")] == [1, 2]                  # empty groups 3, 4 skipped (W7)
    el = {(r["group"], r["id"]): r for r in d.rows("elements")}
    assert el[(1, 1)]["etype"] == 2                                       # below grade: excavated (D-HOU-01)
    assert el[(2, 1)]["etype"] == 1 and el[(2, 1)]["n3"] == 12 and el[(2, 1)]["ki"] == "000000"
    assert [(r["no"], r["thick"]) for r in d.rows("sitelayers")] == [(1, 5.0), (1, 5.0), (2, 0.0)]
    assert [r["number"] for r in d.rows("freqs")] == [2, 4, 8, 16]
    assert d["df"] == 1 / (0.005 * 4096)
    m = d.rows("masses")[0]
    assert (m["node"], m["mx"], m["units"]) == (11, 32.2, 1)
    assert (d["xc"], d["yc"], d["zc"]) == (1.0, 2.0, 3.0)
    assert d["title"] == "AFWRITE demo, embedded block" and d["model"] == "demo"
    text = (mdir / "demo.hou").read_text()
    assert "model_hash = " + ui.model.model_hash() in text


def test_site_point_force_analys_decks(tmp_path):
    ui, mdir = _setup(tmp_path)
    ui.execute("AFWRITE")
    s = _deck(mdir, "SITE")
    assert [r["no"] for r in s.rows("layers")] == [1, 1] and s.rows("halfspace")[0]["no"] == 2
    assert s.rows("halfspace")[0]["vs"] == 1200.0
    assert s.rows("waves") == [dict(type=2, opt=1, ratio1=1.0, ratio2=1.0, angle=0.0)]   # new-model default
    assert s["hslaw"] == "uniform" and s["freq2"] == 2048
    p = _deck(mdir, "POINT")
    assert (p["layer"], p["rad"], p["dim"], p["df"]) == (1, 4.5, 2, 1 / (0.005 * 4096))
    f = _deck(mdir, "FORCE")
    assert sorted((r["node"], r["dof"], r["factor"], r["arrival"]) for r in f.rows("loads")) == \
        [(10, 1, 1.0, 0.1), (11, 2, 2.0, 0.0)]
    a = _deck(mdir, "ANALYS")
    assert (a["ang"], a["xc"], a["coh"], a["freq"]) == (30.0, 1.0, 0, 1)


def test_motion_duplicates_merged_and_stress_reldisp(tmp_path):
    ui, mdir = _setup(tmp_path)
    ui.execute("AFWRITE")
    infos = ui.sink.texts(Kind.INFO)
    assert any("duplicate MOTION output requests merged" in t for t in infos)
    d = _deck(mdir, "MOTION")
    rows = {(r["node"], r["dir"]): [r[f"c{k}"] for k in range(1, 7)] for r in d.rows("nout")}
    assert rows[(2, 1)] == [1, 1, 1, 0, 1, 1]                             # OR of the two X requests
    assert rows[(1, 3)] == [1, 1, 0, 0, 1, 1] and len(rows) == 5
    assert [r["value"] for r in d.rows("damp")] == [0.02, 0.05]
    st = _deck(mdir, "STRESS")
    assert st.rows("eout")[0]["group"] == 1 and st["interopt"] == 1
    rd = _deck(mdir, "RELDISP")
    assert rd["numfiles"] == 1 and rd.rows("rdnd")[0]["x"] == 1


def test_ut16_nfft_rounding_and_records(tmp_path):
    ui, mdir = _setup(tmp_path)
    ui.execute("SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,3000,1")
    ui.execute("AFWRITE")
    s = _deck(mdir, "SITE")
    assert s["nft"] == 2048 and s["df"] == 1 / (0.005 * 2048)
    assert any("NFFT 3000 written as 2048" in t for t in ui.sink.texts(Kind.INFO))
    rep = ui.session["check_report"]
    assert rep.has("Warning", 9, "SITE") and not rep.has("Error", "EDU-03")
    # 3000 records > 2048: EDU-03 error for the modules reading the history
    (mdir / "H1.acc").write_text("0.005\n" + "\n".join("0.01" for _ in range(3000)) + "\n")
    ui.execute("AFWRITE")
    rep = ui.session["check_report"]
    assert rep.has("Error", "EDU-03", "MOTION") and rep.has("Error", "EDU-03", "STRESS")
    assert rep.blocked("MOTION") and not rep.blocked("SITE")


def test_ut16_frequency_grid_values(tmp_path):
    ui, mdir = _setup(tmp_path)
    ui.execute("SITE,0,1,0.9765625,20,2,1,0,1,2048,1,0,0.005,4096,1")
    ui.execute("FREQ,1,0")
    ui.execute("FREQ,1,1,3,5,7,9,11,13,15,16,18")
    ui.execute("AFWRITE")
    s = _deck(mdir, "SITE")
    hz = [r["number"] * s["df"] for r in s.rows("freqs")]
    assert [round(f, 3) for f in hz] == [0.977, 2.930, 4.883, 6.836, 8.789, 10.742, 12.695, 14.648, 15.625, 17.578]
    ui.execute("FREQ,1,18")
    ui.execute("AFWRITE")
    rep = ui.session["check_report"]
    assert rep.has("Error", "EDU-23", "SITE") and rep.blocked("SITE")


def test_soil_and_equake_decks(tmp_path):
    ui, mdir = _setup(tmp_path, extra="SPRO,1,1,Sand\nSPRO,2,1,Sand\nSPRO,3,2,Rock\n"
                                      "DYNP,1,0.0001,1,0.0001,0.5,Sand\nDYNP,2,0.1,0.5,0.1,10,Sand\n"
                                      "DYNP,1,0.0001,1,0.0001,0.4,Rock\nSOIL,100,32.2,0,1,1,8,0.65,1,0\n"
                                      "SRS,1,1,0\nSSAF,1,1,0,1,3,0.125,surface / rock, 5%\nEDUOPT,SOILCUTOFF,25\n"
                                      "EQUAKE,0,,11975,0.05,20,0,1\nRSIN,1,H1.rsi\nRSOUT,1,o.rso\nACCOUT,1,o.acc")
    ui.execute("AFWRITE")
    assert (mdir / "demo.soi").exists(), ui.sink.texts(Kind.WARNING)
    from sassi.modules.soil import read_deck
    d = read_deck(mdir / "demo.soi")
    assert [r["layer"] for r in d.rows("profile")] == [1, 2, 3] and d.rows("profile")[2]["vs"] == 1200.0
    assert d.rows("ssaf")[0]["freqstep"] == 0.125 and d.rows("ssaf")[0]["title"] == "surface / rock, 5%"
    assert d["cof"] == 100.0 and d["soilcutoff"] == 25.0 and d["cl"] == 1 and d["max"] == 0.1
    kinds = [(r["label"], r["kind"]) for r in d.rows("dynp")]
    assert kinds[:4] == [("Sand", "G"), ("Sand", "G"), ("Sand", "D"), ("Sand", "D")]
    e = _deck(mdir, "EQUAKE")
    assert e["nrfreq"] == 3                                                # blank: records of RSIN 1
    sp = e.rows("spectra")[0]
    assert (sp["no"], sp["rsin"], sp["rsout"], sp["accout"]) == (1, "H1.rsi", "o.rso", "o.acc")


def test_simulation_commands_file(tmp_path):
    ui, mdir = _setup(tmp_path, aopt="AOPT,0,0,0,1,1,1,0,0,1,0,1,0,0,0")
    ui.write_options["sim"] = True
    ui.execute("AFWRITE")
    lines = (mdir / "demo-Sim.pre").read_text().splitlines()
    assert [ln for ln in lines if not ln.startswith("*")] == ["RUNSITE", "RUNPOINT", "RUNHOUSE", "RUNANALYS",
                                                              "RUNMOTION"]
    ui.write_options["sim_location"] = "pre"
    ui.execute("WRITE")
    ui.execute("AFWRITE")
    ui.execute("AFWRITE")                                   # the block is replaced, not repeated
    text = (mdir / "demo.pre").read_text()
    assert text.count("RUNSITE") == 1 and text.rstrip().endswith("* ---- end of simulation commands")
    assert text.index("AOPT,") < text.index("RUNSITE")


def test_every_written_deck_reads_back_with_all_parameters(tmp_path):
    ui, mdir = _setup(tmp_path, extra="SPRO,1,1,Sand\nSPRO,2,2,Sand\nDYNP,1,0.0001,1,0.0001,0.5,Sand\n"
                                      "SOIL,100,32.2,0,1,1,8,0.65,1,0\nEQUAKE,0,3,11975,0.05,20,0,1\n"
                                      "RSIN,1,H1.rsi\nRSOUT,1,o.rso\nACCOUT,1,o.acc")
    ui.execute("AFWRITE")
    st = ui.session["afwrite"][0]
    assert set(st["written"]) == {"EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "MOTION",
                                  "STRESS", "RELDISP"}, ui.session["check_report"].format()
    for mod in st["written"]:
        if mod != "SOIL":
            _deck(mdir, mod)


def test_builders_assign_every_schema_parameter(tmp_path, monkeypatch):
    """ARCHITECTURE section 3: AFWRITE fills every schema parameter with resolved data (none is left at
    the schema default by omission)."""
    import sassi.prep.afwrite as AF
    sentinel = object()
    real_new = decks.new

    def fake_new(module):
        d = real_new(module)
        for k in d.params:
            d.params[k] = sentinel
        return d

    ui, mdir = _setup(tmp_path, extra="SPRO,1,1,Sand\nSPRO,2,2,Sand\nDYNP,1,0.0001,1,0.0001,0.5,Sand\n"
                                      "SOIL,100,32.2,0,1,1,8,0.65,1,0\nEQUAKE,0,3,11975,0.05,20,0,1\n"
                                      "RSIN,1,H1.rsi\nRSOUT,1,o.rso\nACCOUT,1,o.acc")
    monkeypatch.setattr(AF.decks, "new", fake_new)
    b = AF.DeckBuilder(ui.model, [mdir])
    for mod in AF.DECK_MODULES:
        d = b.build(mod)
        left = sorted(k for k, v in d.params.items() if v is sentinel)
        assert not left, f"{mod}: {left}"


def test_tab_snapshot_for_the_gui(tmp_path):
    from sassi.prep import options as O
    ui, mdir = _setup(tmp_path)
    snap = O.tab_snapshot(ui.model, "SITE")
    assert snap["SITE"].nl == 20 and snap["SITE"].hs == 2 and snap["shared"]["delt"] == 0.005
    assert snap["SITEX"].soilmode == 0 and snap["WAVE"] == []
    mot = O.tab_snapshot(ui.model, "MOTION")
    assert mot["THFILE"] == "H1.acc" and mot["shared"]["damp"] == [0.02, 0.05]


def test_combin_option_file(tmp_path):
    ui, mdir = _setup(tmp_path, extra="EDUOPT,COMBINDUP,PREFER82", aopt="AOPT,0,0,0,0,0,0,0,0,0,1,0,0,0,0")
    ui.execute("AFWRITE")
    assert (mdir / "COMBIN.opt").read_text().strip() == "EDUOPT,COMBINDUP,PREFER82"
    from sassi.modules.combin import read_options
    assert read_options(mdir)["COMBINDUP"] == "PREFER82"


# ======================================================================================
# Shared free-field layers gate HOUSE as well as SITE (D-CHK-01); one failing deck never stops
# the others (D-AFW-01)
# ======================================================================================
def test_undefined_topl_layer_blocks_site_and_house_only(tmp_path):
    """TOPL names an undefined L layer: Error 19 under SITE *and* HOUSE (HOUSE copies the TOPL layers
    into ``sitelayers``); both decks are withheld (stale ones renamed), the others are written."""
    ui, mdir = _setup(tmp_path, aopt="AOPT,0,0,0,1,1,1,0,1,1,0,1,0,0,0")
    assert ui.execute("AFWRITE")
    assert (mdir / "demo.sit").exists() and (mdir / "demo.hou").exists()
    ui.execute("TOPL,0")
    ui.execute("TOPL,1,3")
    ui.sink.clear()
    assert ui.execute("AFWRITE")
    assert not any("internal error" in e for e in ui.sink.texts(Kind.ERROR) + ui.sink.texts(Kind.WARNING))
    err = (mdir / "demo.err").read_text()
    house = err.split("Errors and Warnings for HOUSE")[1].split("Errors and Warnings for")[0]
    assert "Error 19 : Soil Layer 3 Is not Defined" in house
    for ext in (".sit", ".hou"):
        assert not (mdir / f"demo{ext}").exists() and (mdir / f"demo{ext}.bak").exists(), ext
    for ext in (".poi", ".frc", ".anl", ".mot"):
        assert (mdir / f"demo{ext}").exists(), ext


def test_house_site_layer_values_and_halfspace_are_checked(tmp_path):
    from sassi.prep.check import run_check
    ui, mdir = _setup(tmp_path, aopt="AOPT,0,0,0,1,0,1,0,0,0,0,0,0,0,0")
    rep = run_check(ui.model)
    assert not rep.errors("HOUSE")
    ui.execute("L,1,5,0.12,1500,-800,0.05,0.05")                      # Vs < 0 in a TOPL layer
    rep = run_check(ui.model)
    assert rep.has("Error", 23, "SITE") and rep.has("Error", 23, "HOUSE")
    assert not any(m.number == "EDU-11" for m in rep.of("HOUSE"))      # nu warning only under SITE
    ui.execute("L,1,5,0.12,1500,800,0.05,0.05")
    ui.execute("SITE,0,1,0,20,7,1,0,1,2048,1,0,0.005,4096,1")          # <hs> = 7 undefined
    rep = run_check(ui.model)
    assert rep.has("Error", 19, "SITE") and rep.has("Error", 19, "HOUSE")
    assert any(m.text == "Soil Layer 7 Is not Defined" for m in rep.errors("HOUSE"))
    ui.execute("SITE,0,1,0,0,0,1,0,1,2048,1,0,0.005,4096,1")           # rigid base, no half-space layer
    rep = run_check(ui.model)
    assert not rep.errors("SITE") and not rep.errors("HOUSE")
    ui.sink.clear()
    assert ui.execute("AFWRITE")
    assert any("rigid base (SITE <nl> = 0)" in t for t in ui.sink.texts(Kind.INFO))
    rows = _deck(mdir, "HOUSE").rows("sitelayers")
    assert [r["no"] for r in rows] == [1, 1, 1] and rows[-1]["thick"] == 0.0
    assert [r["no"] for r in _deck(mdir, "SITE").rows("halfspace")] == [1]


def test_house_without_top_layers_writes_no_site_layers(tmp_path):
    ui, mdir = _setup(tmp_path, extra="TOPL,0", aopt="AOPT,0,0,0,0,0,1,0,0,0,0,0,0,0,0")
    assert ui.execute("AFWRITE")
    assert _deck(mdir, "HOUSE").rows("sitelayers") == []


def test_a_failing_deck_builder_blocks_only_its_module(tmp_path, monkeypatch):
    """An unexpected failure while building one deck (CHECK passed) blocks only that module: no deck,
    stale deck renamed, reason in the window and in the .err file; the other decks are written."""
    from sassi.prep.afwrite import DeckBuilder
    ui, mdir = _setup(tmp_path, aopt="AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0")
    assert ui.execute("AFWRITE")
    assert (mdir / "demo.poi").exists()

    def boom(self):
        raise RuntimeError("simulated builder failure")

    monkeypatch.setattr(DeckBuilder, "point", boom)
    ui.sink.clear()
    assert ui.execute("AFWRITE")
    assert not (mdir / "demo.poi").exists() and (mdir / "demo.poi.bak").exists()
    for ext in (".sit", ".hou", ".anl"):
        assert (mdir / f"demo{ext}").exists(), ext
    warns = ui.sink.texts(Kind.WARNING)
    assert any(w.startswith("AFWRITE: POINT: input deck not written (internal error while building the deck")
               for w in warns), warns
    assert "POINT input deck not written -- internal error" in (mdir / "demo.err").read_text()
    st = ui.session["afwrite"][ui.active_model]
    assert "POINT" in st["blocked"] and "POINT" not in st["written"]


def test_eduopt_soilcutoff_accepts_fortran_exponents(tmp_path):
    ui, mdir = _setup(tmp_path, extra="SOILX,0,0,0.1\nSPRO,1,1,S\nSPRO,2,2,S\nDYNP,1,0.001,1,0.001,1,S\n"
                                      "SOIL,100,32.2,0,1,1,8,0.65,1,0\nEDUOPT,SOILCUTOFF,1.25D1",
                      aopt="AOPT,0,1,0,0,0,0,0,0,0,0,0,0,0,0")
    assert ui.execute("AFWRITE"), ui.sink.texts(Kind.ERROR)
    from sassi.modules.soil import read_deck
    assert read_deck(mdir / "demo.soi")["soilcutoff"] == 12.5
