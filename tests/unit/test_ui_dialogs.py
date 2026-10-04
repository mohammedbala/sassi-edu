"""Options dialogs of the GUI (requirements 5.4, D-UI-03, D-UI-04, UI-06, rule L17).

* every Analysis tab exists in the manual's order and shows every field of the 5.4 table, bound to
  the command (record / entry / list / request) that stores it;
* the values of a new model are the new-model defaults of the 5.4 "Default" column (D-UI-03);
* a commit produces the command text a user would type; shared variables are stored once;
* the enable/disable rules that force values (D-UI-04) are applied on commit, with a message;
* invalid values are refused with the Chapter 10 text (UI-06).
"""
from __future__ import annotations

import pytest

from sassi.prep import Interpreter
from sassi.prep.options import TABS
from sassi.ui import dialogs as D


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def _paths(tab):
    out = []
    for g in tab["groups"]:
        for it in g["items"]:
            out.append(it.get("path") or (("sel." + it["sel"]) if it.get("sel") else it.get("label")))
    return out


#: requirements 5.4: the storage of every Analysis tab field (dialog label -> path)
REQUIRED = {
    "EQUAKE": ["sel.spec", "RSIN[spec].file", "RSOUT[spec].file", "ACCOUT[spec].file", "EQUAKE.accopt", "ACCIN[spec].file",
               "EQUAKE.tpsd", "TPSD[spec].file", "EQUAKE.nrfreq", "EQUAKE.rand", "EQUAKE.damp", "SITE.delt",
               "EQUAKE.dur", "EQUAKE.seeds", "EQUAKE.corr", "@CORR", "$EQTIT"],
    "SOIL": ["SITE.nft", "SITE.delt", "SOIL.nrval", "SOILX.mult", "SOILX.max", "SOIL.grav", "SOIL.header", "SOILX.indir",
             "SITE.cl", "SOILX.cl", "$THFILE", "SOILX.file", "SOIL.outcrop", "SOIL.save", "SOIL.iter", "SOIL.ratio",
             "sel.layer", "SPRO[layer].prop", "SPRO[layer].dynprop", "SACC[layer].opt", "SACC[layer].outcrop",
             "SRS[layer].save", "SRS[layer].outcrop", "SOIL.gravmult", "#DAMP", "SSTR[layer].opt1", "SSTR[layer].opt2",
             "SSTR[layer].opt3", "SSTR[layer].opt4", "SSAF[layer].save", "SSAF[layer].outcrop1", "SSAF[layer].outcrop2",
             "SSAF[layer].layer2", "SSAF[layer].freqstep", "SSAF[layer].title", "Compute Fourier Spectrum",
             # Nonlinear Soil Behavior (spec 05a section 6.6): NLSOIL and the NLSLAYER table
             "%NLSOIL.opt", "%NLSOIL.nsub", "%NLSOIL.dispconv", "%NLSOIL.forceconv", "%NLSOIL.equalit",
             "%NLSOIL.bedint", "%NLSOIL.damptype", "%NLSOIL.mmmult", "%NLSOIL.smmult", "%NLSLAYER"],
    "SITE": ["SITEX.soilmode", "SITE.mode1", "SITE.mode2", "HOUSE.gravity", "SITE.fstep", "SITE.delt", "SITE.nft",
             "SITE.freq", "SITE.nl", "SITE.hs", "#TOPL", "SITE.wopt", "sel.wave", "WAVE[wave].opt", "WAVE[wave].ratio1",
             "WAVE[wave].ratio2", "WAVE[wave].angle", "SITE.freq1", "SITE.freq2", "SITE.cl", "SITE.cm"],
    "POINT": ["POINT.opmode", "POINT.layer", "POINT.rad", "From mesh (RADIUS)"],
    "HOUSE": ["HOUSE.opmode", "HOUSE.dim", "HOUSE.imp", "HOUSE.gravity", "HOUSE.gelev", "HOUSEX.nlssi", "Input Data (.pin)",
              "HOUSEX.optimize", "HOUSE.coh", "INCOH.gammax", "INCOH.gammay", "INCOH.gammaz", "INCOH.alpha", "INCOH.ngp",
              "INCOH.nmodes", "INCOH.ipr", "HOUSE.me", "sel.motion", "ME[motion].nfirst", "ME[motion].nlast",
              "ME[motion].xc", "ME[motion].yc", "ME[motion].zc", "#AMP[motion]", "HOUSE.cmplxspec", "Non-Uniform Motion",
              "Non-Uniform Soil", "HOUSE.wpass", "WPASS.appv", "WPASS.ang", "WPASS.cohf", "INCOH.@stoch", "INCOH.hseed",
              "INCOH.vseed", "INCOH.randphz", "HOUSEX.nsim", "HOUSEX.supmode", "HOUSEX.ansys"],
    "FORCE": ["FORCE.opmode", "HOUSE.gravity", "SITE.fstep", "SITE.delt", "SITE.nft", "SITE.freq"],
    "ANALYS": ["ANALYS.opmode", "ANALYS.type", "ANALYS.mode", "ANALYS.simul", "ANALYS.save", "ANALYSX.delrst",
               "ANALYS.fopt", "SITE.freq", "ANALYS.xc", "ANALYS.yc", "ANALYS.zc", "ANALYS.ang", "HOUSE.coh", "HOUSE.wpass",
               "ANALYSX.ffm", "ANALYS.prnt", "HOUSE.me", "ME[motion].nfirst", "ANALYS.impe"],
    "MOTION": ["MOTION.opmode", "ANALYS.type", "MOTION.bl", "MOTION.freq1", "MOTION.freq2", "MOTION.fstep", "#DAMP",
               "MOTION.out", "MOTION.cplx", "MOTIONX.f1213", "MOTION.dur", "MOTIONX.resp", "MOTIONX.srss",
               "MOTION.interp", "MOTION.pzadj", "MOTION.smo", "@NOUT", "SITE.nft", "SITE.delt", "MOTION.mult",
               "MOTION.max", "MOTION.rec1", "MOTION.rec2", "$THTIT", "$THFILE", "MOTION.fopt", "MOTION.cnvrt",
               "MOTIONX.savetf", "MOTIONX.saveacc", "MOTIONX.savers", "MOTIONX.saverot", "MOTIONX.rsttf",
               "MOTIONX.rstacc", "MOTIONX.rstrs", "Save Binary Database"],
    "STRESS": ["STRESS.opmode", "ANALYS.type", "STRESS.iter", "STRESS.save", "STRESS.itran", "STRESSX.pzadj",
               "STRESS.interopt", "STRESSX.smo", "STRESSX.skip", "SITE.nft", "SITE.delt", "SITE.freq", "MOTION.mult",
               "MOTION.max", "MOTION.rec1", "MOTION.rec2", "$THTIT", "$THFILE", "MOTION.fopt", "@EOUT",
               "STRESSX.savemax", "STRESSX.saveth", "STRESSX.rstns", "STRESSX.rstsp", "Save Time History"],
    "RELDISP": ["$RELFILE", "RELD.reldisoutput", "SITE.nft", "SITE.delt", "MOTION.mult", "MOTION.max", "MOTION.rec1",
                "MOTION.rec2", "$THTIT", "$THFILE", "MOTION.fopt", "@RDND", "RELD.reldispsall", "RELDX.saverot",
                "RELDX.rstframes"],
    # spec 05d section 3.8: EQL, the BBC selector / type / yield / X-Y grid, panel and spring data, beams greyed
    "NONLINEAR": ["%EQL.disp", "%EQL.dampcutoff", "%EQL.dampscale", "Material Parameter", "%EQL.nonlinopts",
                  "%EQL.elasticd", "sel.bbc", "%BBC[bbc].type", "%BBC[bbc].yield", "%BBC[bbc].points", "%P", "%S",
                  "Beam", "Beam End 2"],
    "AFWRITE": ["AOPT.equake", "AOPT.soil", "AOPT.dep1", "AOPT.site", "AOPT.point", "AOPT.house", "AOPT.dep2",
                "AOPT.force", "AOPT.analys", "AOPT.combin", "AOPT.motion", "AOPT.stress", "AOPT.reldisp", "AOPT.panel"],
}


def test_analysis_tabs_in_manual_order_with_every_field():
    form = D.form("ANALYSIS")
    assert [t["name"] for t in form["tabs"]] == list(TABS) == ["EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE",
                                                             "ANALYS", "MOTION", "STRESS", "RELDISP", "NONLINEAR",
                                                             "AFWRITE"]
    for tab in form["tabs"]:
        have = _paths(tab)
        missing = [p for p in REQUIRED[tab["name"]] if p not in have]
        assert not missing, (tab["name"], missing)


def _check_paths(form, v):
    for tab in form["tabs"]:
        for g in tab["groups"]:
            for it in g["items"]:
                p = it.get("path")
                if not p or p[0] in "$#@" or p.startswith("sel."):
                    continue
                if p.startswith("%"):              # command records outside OPTION_SPECS (sassi/ui/cmdrecords.py)
                    if "." not in p:
                        assert p[1:] in v["xtables"], p
                        continue
                    rec, field = p[1:].split(".")
                    if "[" in rec:
                        assert field in v["xindexed_defaults"][rec.split("[")[0]], p
                    else:
                        assert field in v["xrecords"][rec], p
                    continue
                rec, field = p.split(".")
                if "[" in rec:
                    name = rec.split("[")[0]
                    assert field in v["indexed_defaults"][name], p
                else:
                    assert field.startswith("@") or field in v["records"][rec], p


def test_every_record_path_is_a_real_field(ui):
    """Each bound path names a field of the typed option record or of the command-record family (so OK stores it
    in its command); also for the Modules-menu dialogs ANSYS Eq. Static Load / ANSYS Dynamic Load."""
    _check_paths(D.form("ANALYSIS"), D.values(ui))
    for name in ("LOADGEN", "LOADGENDYN"):
        _check_paths(D.form(name), D.values(ui, name))


def test_nonlinear_tab_is_bound_and_only_unavailable_controls_are_disabled():
    """The NONLINEAR tab (Option NON) and the SOIL Nonlinear Soil group are editable now: only the controls the
    manual itself greys (Material Parameter, nonlinear beams) stay disabled; no "not editable" note is left."""
    tabs = {t["name"]: t for t in D.form("ANALYSIS")["tabs"]}
    items = [it for g in tabs["NONLINEAR"]["groups"] for it in g["items"]]
    greyed = [it["label"] for it in items if it.get("w") == "disabled"]
    assert greyed == ["Material Parameter", "Beam", "Group Num.", "Spring Gr.", "BBC Num", "Force Opt", "Beam End 1",
                      "Beam End 2"]
    bits = [(it["label"], it["bit"], bool(it.get("disabled"))) for it in items if it.get("w") == "bit"]
    assert bits == [("Use Non-linear Panels", 1, False), ("Use Non-linear Springs", 2, False),
                    ("Use Non-linear Beams", 4, True)]
    soil = [it for g in tabs["SOIL"]["groups"] for it in g["items"]]
    assert not [it for it in soil if it.get("w") == "disabled" and "Nonlinear" in it.get("label", "")]
    text = repr(D.form("ANALYSIS"))
    assert "not editable" not in text and "Not editable" not in text
    aopt = [it for g in tabs["AFWRITE"]["groups"] for it in g["items"]]
    assert [it["path"] for it in aopt if it.get("disabled")] == ["AOPT.dep1", "AOPT.dep2"]     # LIQUEF, PINT greyed


def test_new_model_defaults_D_UI_03(ui):
    r = D.values(ui)["records"]
    assert (r["EQUAKE"]["rand"], r["EQUAKE"]["damp"], r["EQUAKE"]["dur"], r["EQUAKE"]["seeds"]) == (11975, 0.05, 20.0, 1)
    assert r["EQUAKE"]["nrfreq"] == ""                   # blank = number of records of RSIN 1
    assert (r["SITE"]["nft"], r["SITE"]["delt"], r["SITE"]["nl"], r["SITE"]["freq"]) == (4096, 0.005, 20, 1)
    assert (r["SITE"]["mode1"], r["SITE"]["mode2"], r["SITE"]["freq1"], r["SITE"]["freq2"], r["SITE"]["cm"]) == (1, 1, 1, "", 0)
    assert D.values(ui)["context"]["blank_default"]["SITE.freq2"] == "NFFT/2"     # blank = NFFT/2 (5.4)
    assert (r["SOIL"]["grav"], r["SOIL"]["iter"], r["SOIL"]["ratio"], r["SOIL"]["outcrop"], r["SOIL"]["save"]) == (32.2, 8, 0.65, 1, 1)
    assert (r["SOILX"]["mult"], r["SOILX"]["max"], r["SOILX"]["indir"]) == (0.0, 0.1, 0)
    assert (r["HOUSE"]["gravity"], r["HOUSE"]["gelev"], r["HOUSE"]["dim"], r["HOUSE"]["imp"], r["HOUSE"]["coh"]) == (32.2, 0.0, 2, 0, 0)
    assert (r["INCOH"]["gammax"], r["INCOH"]["gammay"], r["INCOH"]["gammaz"], r["INCOH"]["alpha"], r["INCOH"]["ngp"]) == (0.1, 0.1, 0.2, 0.5, 1)
    assert (r["WPASS"]["appv"], r["WPASS"]["ang"], r["WPASS"]["cohf"]) == (1e9, 0.0, 1)
    assert (r["ANALYS"]["prnt"], r["ANALYS"]["impe"], r["ANALYS"]["simul"], r["ANALYSX"]["ffm"]) == (1, 0, 0, 0)
    m = r["MOTION"]
    assert (m["freq1"], m["freq2"], m["fstep"], m["mult"], m["max"], m["rec1"], m["interp"], m["pzadj"], m["smo"], m["cplx"]) == \
        (0.1, 100.0, 301, 1.0, 0.0, 1, 1, 0, 0.0, 1)
    assert r["MOTIONX"]["resp"] == 2 and r["STRESS"]["save"] == 1 and r["STRESS"]["interopt"] == 1
    assert r["RELD"]["reldisoutput"] == 1
    a = r["AOPT"]
    assert [k for k, v in a.items() if v] == ["site", "point", "house", "analys", "motion"]        # D-AFW-05
    w = D.values(ui)["wave_defaults"]
    assert w["2"]["opt"] == 1 and w["1"]["opt"] == 0 and w["3"]["opt"] == 0        # vertical SV only
    assert D.values(ui, "MODEL")["records"]["MOPT"] == {"incomp": 1, "matrix": 0, "mass": 1, "force": 1}
    assert D.values(ui, "CHECK")["records"]["CHECK"] == {"show_warnings": 1, "show_errors": 1, "suppress_window": 0,
                                                         "break_at": 100}


@pytest.mark.parametrize("payload, expected", [
    ({"records": {"EQUAKE": {"nrfreq": 24, "dur": 25}}}, ["EQUAKE,0,24,11975,0.05,25,0,1"]),
    ({"records": {"EQUAKE": {"tpsd": 1}}, "indexed": {"TPSD": {"1": {"file": "t.psd"}}}},
     ["EQUAKE,0,,11975,0.05,20,0,1,1", "TPSD,1,t.psd"]),
    ({"indexed": {"RSIN": {"2": {"file": "rs y.rsi"}}}}, ["RSIN,2,rs y.rsi"]),
    ({"strings": {"EQTIT": "Design, horizontal"}}, ["EQTIT,Design, horizontal"]),
    ({"records": {"SOIL": {"nrval": 3000, "iter": 6}}}, ["SOIL,3000,32.2,0,1,1,6,0.65,1,0"]),
    ({"records": {"SOILX": {"indir": 1}}}, ["SOILX,1"]),
    ({"indexed": {"SPRO": {"2": {"prop": 3, "dynprop": "Sand, dense"}}}}, ["SPRO,2,3,Sand, dense"]),
    ({"indexed": {"SACC": {"1": {"opt": 2}}}}, ["SACC,1,2,0"]),
    ({"indexed": {"WAVE": {"1": {"opt": 2, "ratio1": 0.5, "ratio2": 0.5}, "2": {"ratio1": 0.5, "ratio2": 0.5}}}},
     ["WAVE,1,2,0.5,0.5,0", "WAVE,2,1,0.5,0.5,0"]),
    ({"records": {"SITEX": {"soilmode": 1}}}, ["SITEX,1"]),
    ({"records": {"POINT": {"layer": 2, "rad": 4.5}}}, ["POINT,0,2,4.5"]),
    ({"records": {"HOUSE": {"gelev": -10, "imp": 1}}}, ["HOUSE,32.2,-10,0,2,1,0,0,0,0"]),
    ({"records": {"HOUSEX": {"supmode": 1}}}, ["HOUSEX,0,1"]),
    ({"records": {"FORCE": {"opmode": 1}}}, ["FORCE,1"]),
    ({"records": {"ANALYS": {"type": 1, "impe": 2}}}, ["ANALYS,0,1,0,0,1,0,0,0,0,0,2"]),
    ({"records": {"ANALYSX": {"delrst": 1}}}, ["ANALYSX,0,1"]),
    ({"records": {"MOTION": {"interp": 6, "bl": 1}}}, ["MOTION,0,0,0,0,0,0.1,100,301,1,0,1,0,0,1,0,1,0,0,6"]),
    ({"records": {"MOTIONX": {"rsttf": 1}}}, ["MOTIONX,0,2,0,0,0,0,0,1"]),
    ({"records": {"STRESS": {"itran": 1}}, "records_": None}, ["STRESS,0,0,1,1,1"]),
    ({"records": {"STRESSX": {"skip": 2}}}, ["STRESSX,0,0,2"]),
    ({"records": {"RELD": {"reldisoutput": 3}}}, ["RELD,3,0,0"]),
    ({"records": {"RELDX": {"rstframes": 1}}}, ["RELDX,0,1"]),
    ({"strings": {"RELFILE": "00415TR_X.TFI"}}, ["RELFILE,00415TR_X.TFI"]),
    ({"records": {"AOPT": {"stress": 1, "motion": 0}}}, ["AOPT,0,0,0,1,1,1,0,0,1,0,0,1,0,0"]),
    ({"lists": {"DAMP": "0.02 0.05;0.07"}}, ["DAMP,0", "DAMP,0.02,0.05,0.07"]),
    ({"lists": {"AMP[2]": "1 1.2 0 1.5"}}, ["AMP,2,0", "AMP,2,1,1.2,0,1.5"]),
    ({"requests": {"RDND": [{"node": 15, "x": 1, "z": 1, "yy": 1}]}}, ["RDND,0", "RDND,15,1,0,1,0,1,0"]),
    ({"requests": {"CORR": [{"no": 1, "time": 0, "val": 0.3}]}}, ["CORR,1,0,0.3"]),
])
def test_commit_command_text(ui, payload, expected):
    payload = {k: v for k, v in payload.items() if v is not None}
    lines, notes = D.commit_commands(ui, payload)
    assert lines == expected
    for ln in lines:                                     # the text is accepted by the interpreter
        assert ui.execute(ln), (ln, ui.sink.texts())


def test_commit_round_trip_values(ui):
    lines, _ = D.commit_commands(ui, {"records": {"SITE": {"nl": 10, "freq2": 1024, "delt": 0.01, "nft": 2048}},
                                      "lists": {"TOPL": "1 2 3"}})
    for ln in lines:
        assert ui.execute(ln)
    v = D.values(ui)
    assert (v["records"]["SITE"]["nl"], v["records"]["SITE"]["freq2"], v["records"]["SITE"]["nft"]) == (10, 1024, 2048)
    assert v["lists"]["TOPL"] == "1 2 3"
    # shared variable: one storage location, shown on several tabs
    assert D.commit_commands(ui, {"records": {"SITE": {"delt": 0.01}}})[0] == []


def test_requests_nout_eout_replace_the_list(ui):
    ui.execute("GROUP,3,BEAMS")
    lines, _ = D.commit_commands(ui, {"requests": {
        "NOUT": [{"dir": 1, "c1": 1, "c2": 1, "c5": 1, "nodes": "1, 3-6 10"}, {"dir": 3, "c6": 1, "nodes": "7"}],
        "EOUT": [{"group": 3, "elements": "1-28", "codes": [1, 1, 0, 0, 0, 2]}]}})
    assert lines == ["NOUT,0", "NOUT,1,1,1,0,0,1,0,1,3-6,10", "NOUT,3,0,0,0,0,0,1,7",
                     "EOUT,0", "EOUT,1,1,0,0,0,2,0,0,0,0,0,0,3,1-28"]
    for ln in lines:
        assert ui.execute(ln)
    m = ui.model
    assert m.nout[0].nodes == [1, 3, 4, 5, 6, 10] and m.nout[1].dir == 3 and m.eout[0].codes[5] == 2
    v = D.values(ui)["requests"]
    assert v["NOUT"][0]["nodes"] == "1 3-6 10" and v["NOUT"][0]["c5"] == 1
    assert v["EOUT"][0]["group"] == 3 and v["EOUT"][0]["elements"] == "1-28"
    with pytest.raises(D.DialogError):
        D.commit_commands(ui, {"requests": {"NOUT": [{"dir": 1, "nodes": "6-3"}]}})       # descending range


def test_forced_values_D_UI_04(ui):
    lines, notes = D.commit_commands(ui, {"records": {"HOUSE": {"dim": 1, "coh": 1}}})
    assert lines == ["HOUSE,32.2,0,0,1,0,0,0,0,0"] and "Coherent used" in notes[0]          # 2D -> coherent
    lines, notes = D.commit_commands(ui, {"records": {"WPASS": {"cohf": 3}}})
    assert lines == ["HOUSE,32.2,0,0,2,0,0,1,0,0", "WPASS,1000000000,0,3"] and "unlagged coherency model 3" in notes[0]
    lines, notes = D.commit_commands(ui, {"records": {"HOUSE": {"me": 1}}})
    assert lines == ["HOUSE,32.2,0,0,2,0,0,1,1,0"] and "multiple excitation" in notes[0]
    # deterministic incoherency input: seeds and random phase are zero
    lines, _ = D.commit_commands(ui, {"records": {"INCOH": {"@stoch": 0, "hseed": 5, "randphz": 180}}})
    assert lines == []
    lines, _ = D.commit_commands(ui, {"records": {"INCOH": {"@stoch": 1, "hseed": 5, "randphz": 180}}})
    assert lines == ["INCOH,0.1,0.1,0.2,0.5,1,0,0,0,5,0,180"]


@pytest.mark.parametrize("payload, text", [
    ({"records": {"SITE": {"nl": 25}}}, "Error 47 : Illegal Number of Layers for Halfspace Simulation"),
    ({"records": {"SITE": {"mode1": 0, "mode2": 0}}}, "Error 45 :"),
    ({"records": {"POINT": {"rad": -1}}}, "Error 57 :"),
    ({"records": {"EQUAKE": {"rand": 0}}}, "Error 90 :"),
    ({"records": {"SOIL": {"nrval": 100, "ratio": 1.5}}}, "Error 106 :"),
    ({"records": {"INCOH": {"gammax": 0.05}}}, "Error 58 :"),
    ({"records": {"WPASS": {"appv": 0}}}, "Error 113 :"),
    ({"records": {"ANALYS": {"ang": 400}}}, "Error 63 :"),
    ({"records": {"MOTION": {"mult": 1, "max": 0.3}}}, "Error 78 :"),
    ({"indexed": {"WAVE": {"2": {"ratio1": 1.5}}}}, "Error 54 :"),
    ({"lists": {"DAMP": "0.05 1.2"}}, "Error 72 :"),
    ({"records": {"SITE": {"nft": "4096x"}}}, "not a number"),
])
def test_invalid_values_refused_UI_06(ui, payload, text):
    with pytest.raises(D.DialogError) as ei:
        D.commit_commands(ui, payload)
    assert any(text in p for p in ei.value.problems), ei.value.problems


def test_session_dialogs_write_and_check(ui):
    notes = D.apply_session_dialog(ui, "WRITE", {"records": {"WRITE": {"afwr": 1, "ext": 1}}})
    assert ui.write_options["afwr"] and ui.write_options["ext"] and "afwr" in notes[0]
    with pytest.raises(D.DialogError):
        D.apply_session_dialog(ui, "WRITE", {"records": {"WRITE": {"sim_location": "elsewhere"}}})
    D.apply_session_dialog(ui, "CHECK", {"records": {"CHECK": {"suppress_window": 1, "break_at": "7"}}})
    co = ui.session["check_options"]
    assert co.suppress_window and co.break_at == 7
    form = D.form("WRITE")["tabs"][0]
    assert "WRITE.sim_location" in _paths(form)


def test_parse_helpers():
    assert D.parse_ids("1, 3-6 10;12 - 13") == [1, 3, 4, 5, 6, 10, 12, 13]
    assert D.parse_numbers("0.02\t0.05,\n0.1") == [0.02, 0.05, 0.1]
    with pytest.raises(ValueError):
        D.parse_ids("a-b")
    assert D._ids_text([5, 1, 2, 3, 9]) == "1-3 5 9"


# ------------------------------------------------------------------ the dialog shows what AFWRITE writes
def _run(ui, lines):
    for ln in lines:
        assert ui.execute(ln), (ln, ui.sink.texts()[-3:])


def _effective_waves(ui):
    from sassi.prep.check import Checker
    return sorted((int(w.type), int(w.opt)) for w in Checker(ui.model).effective_waves() if int(w.opt) != 0)


def test_first_wave_page_keeps_the_implicit_field(ui):
    """No WAVE stored: AFWRITE writes the vertical SV field; storing another page must store the SV
    field the dialog shows with it (Checker.effective_waves uses only stored entries once one exists)."""
    v = D.values(ui)
    assert v["wave_defaults"]["2"]["opt"] == 1 and v["wave_defaults"]["1"]["opt"] == 0
    lines, notes = D.commit_commands(ui, {"indexed": {"WAVE": {"1": dict(v["wave_defaults"]["1"], opt=1)}}})
    assert lines == ["WAVE,1,1,1,1,0", "WAVE,2,1,1,1,0"] and "SV-wave field" in notes[0]
    _run(ui, lines)
    assert _effective_waves(ui) == [(1, 1), (2, 1)]
    # the SV page edited in the same commit is stored with its new values (not the default)
    ui2 = Interpreter(cwd=ui.cwd)
    lines, _ = D.commit_commands(ui2, {"indexed": {"WAVE": {"3": {"opt": 1}, "2": {"opt": 1, "ratio1": 0.5, "ratio2": 0.5}}}})
    assert lines == ["WAVE,2,1,0.5,0.5,0", "WAVE,3,1,1,1,0"]
    # switching to SH/L waves in the same commit: the implicit field is then SH (<wopt> = 1)
    ui3 = Interpreter(cwd=ui.cwd)
    v3 = D.values(ui3)
    assert v3["wave_defaults_by_wopt"]["1"]["4"]["opt"] == 1 and v3["wave_defaults_by_wopt"]["1"]["2"]["opt"] == 0
    lines, _ = D.commit_commands(ui3, {"records": {"SITE": {"wopt": 1}}, "indexed": {"WAVE": {"5": {"opt": 1}}}})
    assert lines[1:] == ["WAVE,4,1,1,1,0", "WAVE,5,1,1,1,0"]
    # turning the SV page off is stored explicitly (no wave field at all is then what AFWRITE writes)
    ui4 = Interpreter(cwd=ui.cwd)
    lines, _ = D.commit_commands(ui4, {"indexed": {"WAVE": {"2": {"opt": 0}}}})
    assert lines == ["WAVE,2,0,1,1,0"]


def test_wave_pages_without_entry_are_off_once_a_wave_is_stored(ui):
    _run(ui, ["WAVE,1,2,1,1,0"])                         # R waves only, typed as a command
    v = D.values(ui)
    assert [v["wave_defaults"][str(t)]["opt"] for t in range(1, 6)] == [0, 0, 0, 0, 0]
    assert _effective_waves(ui) == [(1, 2)]
    lines, notes = D.commit_commands(ui, {"indexed": {"WAVE": {"3": {"opt": 1}}}})
    assert lines == ["WAVE,3,1,1,1,0"] and not notes   # no implicit SV added: entries exist already


def test_soil_layer_requests_absent_shown_not_requested(ui):
    """A layer without SACC/SRS/SSTR is 'not requested' (AFWRITE writes no row); the 5.4 values are
    preselected on a first edit, and selecting them emits the command."""
    v = D.values(ui)
    assert v["indexed_defaults"]["SACC"] == {"opt": 0, "outcrop": 0}
    assert v["indexed_defaults"]["SRS"] == {"save": 0, "outcrop": 0}
    assert v["indexed_defaults"]["SSTR"] == {"opt1": 0, "opt2": 0, "opt3": 0, "opt4": 0}
    assert v["indexed_defaults"]["SPRO"] == {"prop": 0, "dynprop": ""}       # other entries: record defaults
    assert (v["request_defaults"]["SACC"]["opt"], v["request_defaults"]["SRS"]["save"],
            v["request_defaults"]["SSTR"]["opt3"]) == (1, 1, 1)
    lines, _ = D.commit_commands(ui, {"indexed": {"SACC": {"2": {"opt": 1, "outcrop": 0}},
                                                  "SRS": {"2": {"save": 1, "outcrop": 0}},
                                                  "SSTR": {"2": {"opt1": 1, "opt2": 0, "opt3": 1, "opt4": 0}}}})
    assert lines == ["SACC,2,1,0", "SRS,2,1,0", "SSTR,2,1,0,1,0"]
    # "No Computation" on a layer without entry is what the deck already does: nothing emitted
    assert D.commit_commands(ui, {"indexed": {"SACC": {"4": {"opt": 0, "outcrop": 0}}}}) == ([], [])
    _run(ui, lines)
    assert D.values(ui)["indexed"]["SACC"]["2"]["opt"] == 1


def test_stochastic_incoherency_input(ui):
    """05b / D-INC-02: stochastic iff (HSeed or VSeed non-zero) and RandPhz > 0; the 5.4 default
    Random Phase Angle 180 is applied when stochastic is chosen; a seedless choice is refused."""
    v = D.values(ui)
    house = dict(v["records"]["HOUSE"], coh=1)
    inc = dict(v["records"]["INCOH"], **{"@stoch": 1})
    with pytest.raises(D.DialogError) as ei:
        D.commit_commands(ui, {"records": {"HOUSE": house, "INCOH": inc}})
    assert "non-zero Horizontal or Vertical SEED" in ei.value.problems[0]
    lines, notes = D.commit_commands(ui, {"records": {"HOUSE": house, "INCOH": dict(inc, hseed=11)}})
    assert lines[-1] == "INCOH,0.1,0.1,0.2,0.5,1,0,0,0,11,0,180" and "180" in notes[0]
    _run(ui, lines)
    after = D.values(ui)["records"]["INCOH"]
    assert after["@stoch"] == 1 and after["randphz"] == 180
    with pytest.raises(D.DialogError):                   # a negative phase is refused
        D.commit_commands(ui, {"records": {"INCOH": dict(after, randphz=-5)}})
    # back to Deterministic (Median): zero phases
    lines, _ = D.commit_commands(ui, {"records": {"INCOH": dict(after, **{"@stoch": 0})}})
    assert lines == ["INCOH,0.1,0.1,0.2,0.5,1,0,0,0,0,0,0"]


def test_incoherency_mixed_state_is_read_as_deterministic_and_kept(ui):
    _run(ui, ["INCOH,0.1,0.1,0.2,0.5,1,0,0,0,0,0,180"])     # phase without seed: deterministic (D-INC-02)
    v = D.values(ui)
    assert v["records"]["INCOH"]["@stoch"] == 0
    assert D.commit_commands(ui, {"records": {"INCOH": v["records"]["INCOH"]}}) == ([], [])     # L17 idempotent
    lines, _ = D.commit_commands(ui, {"records": {"INCOH": dict(v["records"]["INCOH"], randphz=90)}})
    assert lines == ["INCOH,0.1,0.1,0.2,0.5,1,0,0,0,0,0,0"]  # deterministic: a phase edit is not stored


def test_site_frequency2_blank_follows_nfft(ui):
    """<freq2> blank = NFFT/2 (AFWRITE); an unrelated SITE edit keeps it blank, blanking an explicit
    value is a change, and an explicit value is kept."""
    from sassi.prep.options import get_record
    _run(ui, ["SITE,0,1,0,20,0,1,0,1,,1,0,0.005,8192,1"])
    v = D.values(ui)
    assert v["records"]["SITE"]["freq2"] == ""
    lines, _ = D.commit_commands(ui, {"records": {"SITE": dict(v["records"]["SITE"], nl=10, nft=16384)}})
    assert lines == ["SITE,0,1,0,10,0,1,0,1,,1,0,0.005,16384,1"]
    _run(ui, ["SITE,0,1,0,20,0,1,0,1,2048,1,0,0.005,8192,1"])
    v = D.values(ui)
    assert v["records"]["SITE"]["freq2"] == 2048
    lines, _ = D.commit_commands(ui, {"records": {"SITE": dict(v["records"]["SITE"], freq2="")}})
    assert lines == ["SITE,0,1,0,20,0,1,0,1,,1,0,0.005,8192,1"]
    _run(ui, lines)
    assert not get_record(ui.model, "SITE").given(9)
    it = [i for t in D.form("ANALYSIS")["tabs"] for g in t["groups"] for i in g["items"] if i.get("path") == "SITE.freq2"][0]
    assert it["placeholder"] == "NFFT/2" and it["blank_half"] == "SITE.nft"


# ------------------------------------------------------------------ command records outside OPTION_SPECS
#: a 12 m x 4 m wall in the XZ plane (4 x 2 shells, its own material) and a spring (tests/unit/test_nonlinear_cmds.py)
WALL = """N,1,0,0,0
N,5,12,0,0
FILL,1,5
NGEN,2,5,1,5,1,0,0,2
M,1,3.0E7,0.2,24.0,0.04,0.04,1
GROUP,1,SHELL
MACT,1
E,1,1,2,7,6
EGEN,3,1,1
EGEN,1,5,1,4
THICK,1,8,1,0.3
N,20,0,0,-1
SC,1,1000,0,0,0,0,0,0.02
GROUP,2,SPRING
RACT,1
E,1,20,1
GRAVITY,9.81"""

X_KEYS = ("records", "indexed", "strings", "lists", "requests", "xrecords", "xtables", "xindexed")


def _posted_back(v):
    import copy
    return {k: copy.deepcopy(v[k]) for k in X_KEYS if k in v}


def _written(ui):
    from sassi.prep.writer import write_pre
    return write_pre(ui.model)[0].splitlines()


@pytest.fixture
def wall(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(WALL)
    return ui


def test_soil_nonlinear_group_round_trip(ui):
    """SOIL tab, Nonlinear Soil Behavior (NLSOIL) and the NLSLAYER table: open (defaults 0, no rows), change, OK
    (the command text), reopen shows the values, WRITE holds the commands, an unchanged post emits nothing."""
    _run(ui, ["SPRO,1,1,Sand", "SPRO,2,1,Sand", "SPRO,3,2"])
    v = D.values(ui)
    assert v["xrecords"]["NLSOIL"] == {"opt": 0, "nsub": 0, "dispconv": 0.0, "forceconv": 0.0, "equalit": 0,
                                       "bedint": 0, "damptype": 0, "mmmult": 0.0, "smmult": 0.0}
    assert v["xtables"]["NLSLAYER"] == [] and v["context"]["spro_layers"] == [1, 2, 3]
    assert D.commit_commands(ui, _posted_back(v)) == ([], [])
    payload = {"xrecords": {"NLSOIL": dict(v["xrecords"]["NLSOIL"], opt=1, bedint=1, damptype="3", mmmult="0.1")},
               "xtables": {"NLSLAYER": [{"num": 1, "curvefit": 1}, {"num": "2", "b": "1.2", "s": 0.8, "refstrain": 0.03},
                                        {"num": 3, "b": 0}]}}
    lines, _ = D.commit_commands(ui, payload)
    assert lines == ["NLSOIL,1,0,0,0,0,1,3,0.1,0", "NLSLAYER,1,1,0,0,0,0", "NLSLAYER,2,0,1.2,0.8,0.03,0",
                     "NLSLAYER,3,0,0,0,0,0"]
    _run(ui, lines)
    v = D.values(ui)
    assert (v["xrecords"]["NLSOIL"]["opt"], v["xrecords"]["NLSOIL"]["damptype"], v["xrecords"]["NLSOIL"]["mmmult"]) == (1, 3, 0.1)
    assert [r["num"] for r in v["xtables"]["NLSLAYER"]] == [1, 2, 3] and v["xtables"]["NLSLAYER"][1]["b"] == 1.2
    assert D.commit_commands(ui, _posted_back(v)) == ([], [])          # L17: a fixed point
    assert all(ln in _written(ui) for ln in lines)                      # WRITE writes what OK submitted
    # a removed row is deleted, a changed one rewritten (whole line)
    rows = [dict(r) for r in v["xtables"]["NLSLAYER"] if r["num"] != 2]
    rows[0]["curvefit"] = 0
    lines, _ = D.commit_commands(ui, {"xtables": {"NLSLAYER": rows}})
    assert lines == ["DELNLS,2,2,1", "NLSLAYER,1,0,0,0,0,0"]
    _run(ui, lines)
    assert sorted(k for k, _ in ui.model.options.entries("NLSLAYER")) == [1, 3]


def test_nonlinear_tab_round_trip(wall):
    """NONLINEAR tab: EQL (a blank <NonLinOpts> shows the types of the P/S records and stays blank), the BBC grid
    of a curve (complete redefinition, no yield-synchronisation warning), panel and spring tables."""
    from sassi.prep.messages import Kind
    v = D.values(wall)
    assert v["xrecords"]["EQL"] == {"disp": 0.8, "nonlinopts": 0, "dampcutoff": 0.0, "dampscale": 0.0, "elasticd": 0}
    assert (v["xtables"]["P"], v["xtables"]["S"], v["xindexed"]["BBC"]) == ([], [], {})
    assert v["xindexed_defaults"]["BBC"] == {"type": 0, "yield": "", "points": []}
    assert D.commit_commands(wall, _posted_back(v)) == ([], [])
    payload = {"xrecords": {"EQL": dict(v["xrecords"]["EQL"], disp="0.75", dampcutoff=7, elasticd=1)},
               "xtables": {"P": [{"num": 1, "group": 1, "bbc": 1, "disp": 1, "force": 1}],
                           "S": [{"num": 1, "group": 2, "elem": 1, "bbc": 2, "disp": 1, "force": 4}]},
               "xindexed": {"BBC": {"1": {"type": 1, "yield": 2, "points": [{"x": 0.0002, "y": 1500},
                                                                            {"x": "0.004", "y": "3600"},
                                                                            {"x": 0.02, "y": 3700}]},
                                    "2": {"type": 4, "yield": 1, "points": [{"x": 0.01, "y": 10}, {"x": 0.2, "y": 29}]}}}}
    lines, _ = D.commit_commands(wall, payload)
    assert lines == ["EQL,0.75,,7,0,1", "BBCI,1,2,1", "BBCX,1,3,2,0.0002,0.004,0.02", "BBCY,1,3,2,1500,3600,3700",
                     "BBCI,2,1,4", "BBCX,2,2,1,0.01,0.2", "BBCY,2,2,1,10,29", "P,1,1,1,1,1", "S,1,2,1,2,1,4"]
    wall.sink.clear()
    _run(wall, lines)
    assert not wall.sink.texts(Kind.WARNING), wall.sink.texts(Kind.WARNING)
    v = D.values(wall)
    assert v["xrecords"]["EQL"]["nonlinopts"] == 3 and v["context"]["eql_opts_blank"]     # shown = used (P and S)
    assert v["xindexed"]["BBC"]["1"] == {"type": 1, "yield": 2, "points": [{"x": 0.0002, "y": 1500.0},
                                                                        {"x": 0.004, "y": 3600.0},
                                                                        {"x": 0.02, "y": 3700.0}]}
    assert v["xtables"]["S"] == [{"num": 1, "group": 2, "elem": 1, "bbc": 2, "disp": 1, "force": 4}]
    assert D.commit_commands(wall, _posted_back(v)) == ([], [])
    assert all(ln in _written(wall) for ln in lines)
    # a changed curve is redefined; a deleted curve, panel and spring are removed; a type click stores <NonLinOpts>
    pts = v["xindexed"]["BBC"]["1"]["points"][:2]
    lines, _ = D.commit_commands(wall, {"xindexed": {"BBC": {"1": {"type": 1, "yield": 2, "points": pts}, "2": None}},
                                        "xtables": {"P": [], "S": []},
                                        "xrecords": {"EQL": dict(v["xrecords"]["EQL"], nonlinopts=1)}})
    assert lines == ["EQL,0.75,1,7,0,1", "DELBBC,1", "BBCI,1,2,1", "BBCX,1,2,2,0.0002,0.004", "BBCY,1,2,2,1500,3600",
                     "DELBBC,2", "PDEL,1", "DELSPR,1"]
    wall.sink.clear()
    _run(wall, lines)
    assert not wall.sink.texts(Kind.WARNING)
    v = D.values(wall)
    assert sorted(v["xindexed"]["BBC"]) == ["1"] and not v["context"]["eql_opts_blank"]
    assert v["xrecords"]["EQL"]["nonlinopts"] == 1 and v["xtables"]["P"] == []


@pytest.mark.parametrize("payload, text", [
    ({"xrecords": {"EQL": {"disp": 2}}}, "EQL: <disp> = 2: the equivalent-linear displacement factor must lie in (0, 1]"),
    ({"xrecords": {"EQL": {"disp": "abc"}}}, "EQL: <disp>: 'abc' is not a number"),
    ({"xtables": {"P": [{"num": 1, "group": 0}]}}, "P: panel number <num> and group <group> must be >= 1"),
    ({"xtables": {"P": [{"num": 1, "group": 1}, {"num": 1, "group": 1}]}}, "P: panel 1 is given twice"),
    ({"xtables": {"S": [{"num": 1, "group": 2, "elem": 1, "bbc": 1, "disp": 5}]}}, "S: <disp> must be 1 (X), 2 (Y) or 3 (Z)"),
    ({"xtables": {"S": [{"group": 2}]}}, "S: row 1: the spring number is required"),
    ({"xtables": {"NLSLAYER": [{"num": 1, "curvefit": 2}]}}, "NLSLAYER: <curvefit> = 2 (0 parameters given, 1 curve fit)"),
    ({"xtables": {"NLSLAYER": [{"num": 1, "b": -1}]}}, "Beta, S exponent, Reference Strain and Viscosity must be >= 0"),
    ({"xrecords": {"NLSOIL": {"damptype": 7}}}, "NLSOIL <NLDampType> = 7"),
    ({"xindexed": {"BBC": {"1": {"yield": 3, "points": [{"x": 0.01, "y": 10}]}}}}, "BBC: curve 1: Yield Num. 3 must lie in 1..1"),
    ({"xindexed": {"BBC": {"1": {"points": [{"x": 0.01, "y": ""}]}}}}, "BBC: curve 1 point 1: X and Y are required"),
    ({"xindexed": {"BBC": {"0": {"type": 1}}}}, "BBC: curve number '0' must be an integer >= 1"),
    ({"xrecords": {"LOADGEN": {"data": 2}}}, "LOADGEN is not a record of this dialog"),
    ({"xtables": {"EQL": []}}, "EQL is not a table of this dialog"),
    ({"xrecords": {"EQL": {"nosuch": 1}}}, "EQL: unknown field(s) nosuch"),
])
def test_command_record_values_refused_UI_06(wall, payload, text):
    """Bad values of the command records refuse the whole commit with the message of their command (dry run in a
    copy of the model) or of the dialog; nothing reaches the session."""
    before = wall.model.copy()
    with pytest.raises(D.DialogError) as ei:
        D.commit_commands(wall, payload)
    assert any(text in p for p in ei.value.problems), ei.value.problems
    assert wall.model.same_state(before)


def test_command_records_refuse_the_whole_commit(wall):
    """A family error refuses the option-record part of the same commit too (one validated commit, UI-06)."""
    with pytest.raises(D.DialogError) as ei:
        D.commit_commands(wall, {"records": {"SITE": {"nl": 12}}, "xrecords": {"EQL": {"disp": 5}}})
    assert len(ei.value.problems) == 1 and ei.value.problems[0].startswith("EQL:")
    with pytest.raises(D.DialogError) as ei:            # both parts reported when both are wrong
        D.commit_commands(wall, {"records": {"SITE": {"nl": 25}}, "xtables": {"P": [{"num": "x"}]}})
    assert any("Error 47" in p for p in ei.value.problems) and any(p.startswith("P:") for p in ei.value.problems)


def test_loadgen_dialogs_layout(ui):
    """Modules > ANSYS Eq. Static Load / ANSYS Dynamic Load (spec 04 section 15.6): the manual's groups and fields
    plus the SASSI-EDU settings, an Ok and a Run (RUNLOADGEN) button; the enable rules of the manual."""
    st = D.form("LOADGEN")
    assert st["title"] == "ANSYS Static Load Converter" and st["run"]["line"] == "RUNLOADGEN,STATIC"
    titles = [g["title"] for g in st["tabs"][0]["groups"]]
    assert titles[:5] == ["Data to Add From ACS SASSI to the ANSYS model", "SSI Model and Results Input",
                          "Ansys Model and Data Input", "Mass Data for Internal Load (Ignore for Displacement)",
                          "ANSYS Output File"]
    items = {it["path"]: it for g in st["tabs"][0]["groups"] for it in g["items"] if it.get("path")}
    assert [c[1] for c in items["%LOADGEN.data"]["choices"]] == ["Displacement", "Acceleration", "Disp. and Accel.",
                                                                 "Disp. for Soil Module"]
    assert items["%LGFILE.DISP"]["enable"] == ["%LOADGEN.data", "in", [1, 3, 4]]          # disabled for Acceleration
    assert items["%LGFILE.LUMPED"]["enable"][1] == ["%LOADGEN.masstype", "==", 1]
    assert {"%LGFILE.SSIPATH", "%LGFILE.HOUSE", "%LGFILE.ACC", "%LGFILE.ANSYSPATH", "%LGFILE.APDL", "%LOADGEN.genmass",
            "%LOADGEN.multi", "%LGTIME.crit", "%LGNODE.D", "%LGMAP.mode", "%LGOPT.digits"} <= set(items)
    dy = D.form("LOADGENDYN")
    assert dy["title"] == "ANSYS Dynamic Load Converter" and dy["run"]["line"] == "RUNLOADGEN,DYNAMIC"
    items = {it["path"]: it for g in dy["tabs"][0]["groups"] for it in g["items"] if it.get("path")}
    assert {"%LGFILE.SSIPATH", "%LGFILE.HOUSE", "%LGFILE.GROUND", "%LGFILE.ANSYSPATH", "%LOADGENDYN.alpha",
            "%LOADGENDYN.beta", "%LGFILE.APDLDYN", "%LOADGENDYN.method", "%LGNODE.A"} <= set(items)
    assert items["%LOADGENDYN.alpha"]["enable"] == ["%LOADGENDYN.zeta", "==", 0]          # zeta given: computed
    with pytest.raises(KeyError):
        D.form("NOSUCH")


def test_loadgen_dialogs_round_trip(ui):
    """Open (the defaults: Acceleration, Lumped Mass, critical time V,1, identity numbering), change, OK (the
    LOADGEN record commands), reopen, WRITE; a shared file box (HOUSE) is one LGFILE entry for both dialogs."""
    v = D.values(ui, "LOADGEN")
    x = v["xrecords"]
    assert x["LOADGEN"] == {"data": 2, "multi": 0, "rotdisp": 0, "rotacc": 0, "masstype": 1, "genmass": 0,
                            "source": "RESULTS"}
    assert x["LGTIME"]["crit"] == "V" and x["LGTIME"]["n"] == 1 and x["LGMAP"]["mode"] == "IDENTITY"
    assert x["LGOPT"] == {"digits": 12, "opmode": 0, "rest": 1} and set(x["LGNODE"]) == {"D", "M", "A"}
    assert v["context"]["analysis"] == "STATIC" and D.commit_commands(ui, _posted_back(v), "LOADGEN") == ([], [])
    payload = {"xrecords": {
        "LOADGEN": dict(x["LOADGEN"], data=3, genmass=1, masstype=2, source="file8"),
        "LGFILE": dict(x["LGFILE"], HOUSE="my model, X.hou", MASTER="m.masm"),
        "LGNODE": dict(x["LGNODE"], M="5 7-9", D="1-4"),
        "LGTIME": dict(x["LGTIME"], crit="ACC", n=2, tsep=0.5, node=12, dof=1),
        "LGMAP": {"mode": "COORD", "tol": "0.001", "file": "ansys.cdb"},
        "LGOPT": dict(x["LGOPT"], digits=10)}}
    lines, _ = D.commit_commands(ui, payload, "LOADGEN")
    assert lines == ["LOADGEN,3,0,0,0,2,1,FILE8", "LGFILE,HOUSE,my model, X.hou", "LGFILE,MASTER,m.masm",
                     "LGNODE,D,1-4", "LGNODE,M,5,7-9", "LGTIME,ACC,2,0.5,12,1", "LGMAP,COORD,0.001,ansys.cdb",
                     "LGOPT,10,0,1"]
    _run(ui, lines)
    v = D.values(ui, "LOADGEN")
    x = v["xrecords"]
    assert (x["LOADGEN"]["data"], x["LOADGEN"]["masstype"], x["LOADGEN"]["source"]) == (3, 2, "FILE8")
    assert x["LGFILE"]["HOUSE"] == "my model, X.hou" and x["LGNODE"]["M"] == "5 7-9"
    assert (x["LGTIME"]["crit"], x["LGTIME"]["node"], x["LGMAP"]["tol"]) == ("ACC", 12, 0.001)
    assert D.commit_commands(ui, _posted_back(v), "LOADGEN") == ([], [])
    assert all(ln in _written(ui) for ln in lines if not ln.startswith("LGNODE"))
    assert "LGNODE,M,5,7-9" in _written(ui)                          # WRITE: clear + list (LGNODE adds)
    dyn = D.values(ui, "LOADGENDYN")["xrecords"]
    assert dyn["LGFILE"]["HOUSE"] == "my model, X.hou"                 # the same file box in both dialogs
    lines, _ = D.commit_commands(ui, {"xrecords": {
        "LOADGENDYN": dict(dyn["LOADGENDYN"], zeta=0.05, f1=1, f2=20, method="acc", refnode=12),
        "LGNODE": dict(dyn["LGNODE"], M="", A="12 14 16"), "LGMAP": {"mode": "IDENTITY"}}}, "LOADGENDYN")
    assert lines == ["LOADGENDYN,,,ACC,12,RESULTS,0,0,0.05,1,20,0,1", "LGNODE,M,0", "LGNODE,A,12,14,16", "LGMAP"]
    _run(ui, lines)
    assert ui.model.options.record("LGMAP") is None and D.values(ui, "LOADGENDYN")["xrecords"]["LOADGENDYN"]["zeta"] == 0.05


@pytest.mark.parametrize("payload, text", [
    ({"xrecords": {"LOADGEN": {"data": 7}}}, "LOADGEN: <data> '7' is not one of"),
    ({"xrecords": {"LOADGEN": {"multi": 2}}}, "LOADGEN: <multi> must be 0 or 1"),
    ({"xrecords": {"LGTIME": {"crit": "ACC"}}}, "LGTIME: ACC needs a node >= 1 and a DOF 1..6"),
    ({"xrecords": {"LGTIME": {"crit": "TIME", "times": ""}}}, "LGTIME: give at least one time"),
    ({"xrecords": {"LGTIME": {"crit": "MX", "x0": 1}}}, "give all three coordinates of the moment point"),
    ({"xrecords": {"LGNODE": {"D": "5-3"}}}, "LGNODE: D nodes: descending range 5-3"),
    ({"xrecords": {"LGMAP": {"mode": "PAIRS"}}}, "LGMAP: PAIRS needs the numbering file"),
    ({"xrecords": {"LGOPT": {"digits": 30}}}, "LGOPT: APDL significant digits must be 6..17"),
    ({"xrecords": {"LOADGENDYN": {"method": "ACC", "refnode": 0}}}, "method ACC needs the reference node"),
    ({"xrecords": {"LOADGENDYN": {"zeta": 0.05, "f1": 5, "f2": 1}}}, "needs 0 < <f1> < <f2>"),
    ({"records": {"SITE": {"nl": 10}}}, "the LOADGEN dialog has no records"),
    ({"xrecords": {"EQL": {"disp": 0.7}}}, "EQL is not a record of this dialog"),
])
def test_loadgen_dialog_values_refused(ui, payload, text):
    with pytest.raises(D.DialogError) as ei:
        D.commit_commands(ui, payload, "LOADGEN" if "LOADGENDYN" not in payload.get("xrecords", {}) else "LOADGENDYN")
    assert any(text in p for p in ei.value.problems), ei.value.problems


def test_command_records_posted_back_emit_nothing(wall):
    """L17 fixed point over every part of the Analysis window and the LOADGEN dialogs, command records included."""
    _run(wall, ["EQL,0.8,1,7,0,1", "P,1,1,1,1,1", "BBCX,1,2,1,0.01,0.2", "BBCY,1,2,1,10,29", "BBCI,1,1,4",
                "S,1,2,1,1,1,4", "NLSOIL,1,4,0.0001,0,30,0,3,0.2,0.001", "NLSLAYER,1,0,1.2,0.8,0.03,0",
                "LOADGEN,DISPACC,1,0,1,MASTER,1,FILE8", "LGFILE,APDL,out dir/x, y.inp", "LGNODE,D,1-5,7",
                "LGTIME,TIME,1.5,3", "LGMAP,PAIRS,,pairs.txt", "LGOPT,8,1,0", "LOADGENDYN,0.4,0.002,ACC,7,RESULTS,1,1"])
    assert D.commit_commands(wall, _posted_back(D.values(wall))) == ([], [])
    for name in ("LOADGEN", "LOADGENDYN"):
        assert D.commit_commands(wall, _posted_back(D.values(wall, name)), name) == ([], [])
