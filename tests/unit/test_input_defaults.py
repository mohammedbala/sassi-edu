"""Built-in inputs and the defaults of blank inputs (requirements section 7.19, D-W5-01 ... D-W5-14).

* the library ``sassi/data/library`` and its ``@`` names: catalogue, copies identical to their originals, the
  resolution in the interpreter, CHECK, AFWRITE, the modules and the GUI file API, in the local layout, the
  bundled layout of the browser version and the pip package data;
* a new model with a blank THFILE passes CHECK and runs SITE ... MOTION with the RG 1.60 record, identical to the
  run with the record given explicitly; a forced-vibration model takes the Ricker pulse;
* EQUAKE runs with every input blank; SOIL runs with library curves (no INP) and with the default profile;
* every default used is reported: CHECK Warning EDU-29, AFWRITE notes, the RUN listing, LIBRARY,DEFAULTS;
* ``EDUOPT,DEFAULTS,OFF`` gives the ACS errors back; given but missing files stay errors;
* the GUI: placeholders of blank files (the browser logic run in node against the Python values), the Library
  picker data, the library soil curves and the default SOIL profile of the SPRO pages.
"""
from __future__ import annotations

import filecmp
import glob
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from sassi.io import decks, textfiles
from sassi.io import library as LIB
from sassi.prep import Interpreter, Kind
from sassi.prep import defaults as DEF
from sassi.prep.check import run_check
from sassi.ui import dialogs as D

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples" / "data"
DATA = ROOT / "sassi" / "data"

#: a two-node stick on a two-layer site: the smallest model that runs SITE ... MOTION (SI-like numbers, BS gravity)
TINY = """MDL,m,m
TIT,built-in defaults
L,1,5,0.12,600,300,0.05,0.05
L,2,5,0.13,2400,1200,0.02,0.02
TOPL,1,1
SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,4096,1
FREQ,1,4,20,41,82
N,1,0,0,0
N,2,0,0,10
D,1,2,1,1,ROTX,ROTY,ROTZ
GROUP,1,SPRING
SC,1,1e5,1e5,1e5,0,0,0,0.02
E,1,1,2
MT,2,10,10,10
INT,1,1,1,1
POINT,0,0,4.5
HOUSE,32.2,0,0,2,0,0,0,0,0
NOUT,1,1,1,0,0,1,1,2
DAMP,0.05
"""
SEISMIC_LINE = "THFILE blank: the built-in RG 1.60 record @rg160h_030g.acc (0.30 g, 20 s, dt 0.005 s) is used"


def _ui(tmp_path, text: str) -> Interpreter:
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(text)
    return ui


def _errors(ui):
    return ui.sink.texts(Kind.ERROR)


def _infos(ui):
    return ui.sink.texts(Kind.INFO)


def _edu29(rep, module=None):
    return [m.detail for m in rep.warnings(module) if m.number == "EDU-29"]


def _run(ui, *modules):
    for mod in modules:
        assert ui.execute(f"RUN{mod}"), (mod, _errors(ui)[-3:])


# ====================================================================== the library and its @ names
def test_catalogue_files_readme_and_copies():
    files = sorted(p.name for p in LIB.LIBRARY_DIR.iterdir() if p.is_file())
    assert set(files) - {"README.txt"} == {e.name for e in LIB.CATALOGUE}
    readme = (LIB.LIBRARY_DIR / "README.txt").read_text(encoding="utf-8")
    for e in LIB.CATALOGUE:                           # README names every file with its units and source
        assert e.name in readme, e.name
    copies = {"rg160h_030g.acc": EXAMPLES / "rg160h_030g.acc", "rg160h_030g.rsi": EXAMPLES / "rg160h_030g.rsi",
              "ricker_5hz.th": EXAMPLES / "ricker_5hz.th", "rg160h_1g.rsi": DATA / "rg160" / "RG160H_5pct_1g.rsi",
              "rg160v_1g.rsi": DATA / "rg160" / "RG160V_5pct_1g.rsi",
              "rg160h_appa_1g_cm2s3.tpsd": DATA / "rg160" / "RG160H_AppA_1g_cm2s3.tpsd",
              "rg160h_appa_1g_in2s3.tpsd": DATA / "rg160" / "RG160H_AppA_1g_in2s3.tpsd",
              "dynp_library.pre": DATA / "dynp_library.pre"}
    for name, orig in copies.items():
        assert filecmp.cmp(LIB.LIBRARY_DIR / name, orig, shallow=False), name
    from sassi.core import equake_lib as EL
    for tag, units in (("cm2s3", "SI"), ("in2s3", "BS")):            # 0.30 g PSD = 0.09 x the 1 g one
        f, p03 = EL.read_two_columns(LIB.LIBRARY_DIR / f"rg160h_030g_{tag}.tpsd")
        _, p1 = EL.read_two_columns(LIB.LIBRARY_DIR / f"rg160h_appa_1g_{tag}.tpsd")
        np.testing.assert_allclose(p03, 0.09 * p1, rtol=1e-7)
        np.testing.assert_allclose(p03, EL.rg160_target_psd(f, pga=0.3, units=units), rtol=1e-7)
    for e in LIB.CATALOGUE:                           # the catalogue numbers are those of the files
        if e.dt is not None:
            a, dt = textfiles.read_history(e.path)
            assert dt == e.dt and len(a) == e.values, e.name
        elif e.kind in ("spectrum", "psd"):
            assert len(EL.read_two_columns(e.path)[0]) == e.values, e.name
    a, _ = textfiles.read_history(LIB.LIBRARY_DIR / "rg160h_030g.acc")
    assert abs(np.max(np.abs(a)) - 0.324) < 5e-4
    assert LIB.dynp_labels() == ["Clay", "Sand", "Rock"] and len(LIB.dynp_points("Rock")) == 8


def test_library_names_resolve():
    p = LIB.library_path("@rg160h_030g.acc")
    assert p == LIB.LIBRARY_DIR / "rg160h_030g.acc" and p.is_file()
    assert LIB.library_path(' "@RG160H_030G.ACC" ') == p                 # case and quotes do not matter
    for bad in ("@missing.acc", "@../rg160/RG160H_5pct_1g.rsi", "@sub/x.acc", "rg160h_030g.acc", "@", "@README.txt"):
        assert LIB.library_path(bad) is None, bad
    assert LIB.resolve("@ricker_5hz.th", []) == LIB.LIBRARY_DIR / "ricker_5hz.th"
    assert LIB.module_path("@rg160h_030g.rsi", "/x") == LIB.LIBRARY_DIR / "rg160h_030g.rsi"
    assert LIB.module_path("@nope.rsi", "/x") == Path("/x/@nope.rsi")       # not found: the module names it
    assert LIB.starts_with_name("rg160h_030g.acc,1") == len("rg160h_030g.acc")
    assert LIB.starts_with_name("rg160h_030g") == 0 and LIB.starts_with_name("rg160h_030g.accx") == 0


def test_bundled_and_installed_layouts(tmp_path):
    """The browser bundle (web/build.py) carries the library and ``@`` names resolve from the unpacked package
    (/home/pyodide/sassi-edu/sassi in Pyodide; here a copy in tmp_path run by a fresh interpreter); the pip
    package data globs of pyproject.toml include the folder."""
    spec = importlib.util.spec_from_file_location("sassi_web_build", ROOT / "web" / "build.py")
    build = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build)
    rels = [r for r in build.bundle_files() if r.startswith("sassi/")]
    for e in LIB.CATALOGUE:
        assert f"sassi/data/library/{e.name}" in rels, e.name
    assert "sassi/data/library/README.txt" in rels
    bundle = tmp_path / "sassi-edu"
    for rel in rels:
        if rel.endswith(".py") or rel.startswith("sassi/data/library/") or rel == "sassi/ui/command_args.json":
            dst = bundle / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / rel, dst)
    code = ("import sassi, sassi.io.library as L; from sassi.prep import Interpreter; "
            "ui = Interpreter(); ui.execute('INP,@dynp_library.pre'); "
            "print(L.library_path('@rg160h_030g.acc')); print(sassi.__file__); "
            "print(len(ui.model.options.entries('DYNP')))")
    env = dict(os.environ, PYTHONPATH=str(bundle))
    r = subprocess.run([sys.executable, "-c", code], cwd=str(tmp_path), env=env, capture_output=True, text=True,
                       timeout=120)
    assert r.returncode == 0, r.stderr
    lib_path, pkg, ndynp = r.stdout.split()
    assert Path(lib_path) == bundle / "sassi" / "data" / "library" / "rg160h_030g.acc"
    assert Path(pkg).parent == bundle / "sassi" and int(ndynp) == 30
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert '"data/**/*"' in text
    matched = {Path(p).name for p in glob.glob(str(ROOT / "sassi" / "data" / "**" / "*"), recursive=True)}
    assert {e.name for e in LIB.CATALOGUE} <= matched


def test_interpreter_reads_at_names(tmp_path):
    ui = _ui(tmp_path, "THFILE,@rg160h_030g.acc\nRSIN,1,@rg160h_030g.rsi\nINP,@dynp_library.pre\n"
                       "READSPEC,@rg160h_030g.rsi,1,2\nREADTH,@ricker_5hz.th,0,1")
    assert not _errors(ui) and not ui.sink.texts(Kind.WARNING), (_errors(ui), ui.sink.texts(Kind.WARNING))
    assert ui.model.options.string("THFILE") == "@rg160h_030g.acc"
    assert len(ui.model.options.entries("DYNP")) == 30
    # a library name wins over a variable of the same stem; other @ text is still a variable
    ui = _ui(tmp_path, "VAR,rg160h_030g,7\nTHFILE,@rg160h_030g.acc\nTIT,run @rg160h_030g[1]")
    assert ui.model.options.string("THFILE") == "@rg160h_030g.acc" and ui.model.title == "run 7"
    assert not ui.execute("INP,@nothing.pre") and "LIBRARY lists the built-in inputs" in _errors(ui)[-1]
    assert not ui.execute("WRITE,@rg160h_030g.acc") and "read-only" in _errors(ui)[-1]
    assert filecmp.cmp(LIB.LIBRARY_DIR / "rg160h_030g.acc", EXAMPLES / "rg160h_030g.acc", shallow=False)


# ====================================================================== MOTION, STRESS, RELDISP: THFILE
def test_new_model_motion_runs_with_the_default_record(tmp_path):
    """A new model without THFILE: CHECK has no error and names the default (EDU-29), AFWRITE writes the portable
    @ name, SITE ... MOTION run, the listing names the default, and the result equals the explicit file."""
    ui = _ui(tmp_path, TINY + "ANALYS,0,0,0,0,1,0,0,0,0,0,0\nAOPT,0,0,0,1,1,1,0,0,1,0,1,0,0,0\nCHECK\nAFWRITE")
    assert not _errors(ui), _errors(ui)
    rep = ui.session["check_report"]
    assert rep.errors() == [], rep.format()
    assert _edu29(rep, "MOTION") == [SEISMIC_LINE]
    assert f"Warning EDU-29 : Built-In Default Input Used  [{SEISMIC_LINE}]" in _infos(ui)
    assert any(t.startswith(f"AFWRITE: {SEISMIC_LINE} (EDU-29)") for t in _infos(ui))
    mdir = tmp_path / "m"
    assert decks.read(mdir / "m.mot", "MOTION")["thfile"] == "@rg160h_030g.acc"
    _run(ui, "SITE", "POINT", "HOUSE", "ANALYS", "MOTION")
    out = (mdir / "m_MOTION.out").read_text(encoding="utf-8")
    assert "Built-in defaults of blank inputs" in out and SEISMIC_LINE in out
    assert "built-in input @rg160h_030g.acc: RG 1.60 horizontal record" in out
    assert SEISMIC_LINE in [t.strip() for t in _infos(ui)]             # the RUN command shows it too
    acc_default = (mdir / "00002TR_X.ACC").read_bytes()
    # the same run with the record given explicitly (a copy of the original in the model folder)
    shutil.copyfile(EXAMPLES / "rg160h_030g.acc", mdir / "motion.acc")
    assert ui.execute("THFILE,motion.acc") and ui.execute("AFWRITE")
    _run(ui, "MOTION")
    assert (mdir / "00002TR_X.ACC").read_bytes() == acc_default
    out = (mdir / "m_MOTION.out").read_text(encoding="utf-8")
    assert "Built-in defaults" not in out and "built-in input" not in out


def test_forced_vibration_takes_the_ricker_pulse(tmp_path):
    ui = _ui(tmp_path, TINY + "FORCE,0\nF,2,1,0,0\nANALYS,0,1,0,0,1,0,0,0,0,0,0\n"
                              "AOPT,0,0,0,1,1,1,0,1,1,0,1,0,0,0\nCHECK\nAFWRITE")
    rep = ui.session["check_report"]
    assert rep.errors() == [], rep.format()
    line = ("THFILE blank (forced vibration, ANALYS <type> 1): the built-in 5 Hz Ricker pulse @ricker_5hz.th "
            "(peak 1, 2 s, dt 0.005 s) is the load history")
    assert _edu29(rep, "MOTION") == [line]
    _run(ui, "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "MOTION")
    out = (tmp_path / "m" / "m_MOTION.out").read_text(encoding="utf-8")
    assert line in out and "built-in input @ricker_5hz.th: 5 Hz Ricker wavelet" in out
    assert "Reference load history" in out


@pytest.mark.parametrize("extra,why", [
    ("MOTION,0,0,0,0,0,0.1,100,301,1,0,1,0,1,0,0,1,0,0,1", "<fopt> 0 layout"),          # pairs: no default
    ("SITE,0,1,0,20,2,1,0,1,2048,1,0,0.01,4096,1", "dt = 0.005 s, the time step SITE <delt> is 0.01 s"),
    ("EDUOPT,DEFAULTS,OFF", "EDUOPT,DEFAULTS,OFF"),
])
def test_no_default_keeps_error_73_with_the_reason(tmp_path, extra, why):
    ui = _ui(tmp_path, TINY + "ANALYS,0,0,0,0,1,0,0,0,0,0,0\n" + extra)
    rep = run_check(ui.model, modules=["MOTION", "STRESS", "RELDISP"], dirs=[tmp_path / "m"])
    e73 = [m for m in rep.errors() if m.number == 73]
    assert {m.module for m in e73} == {"MOTION", "STRESS", "RELDISP"}
    assert all(why in m.detail and m.detail.startswith("THFILE not given") for m in e73), e73[0].detail
    assert not _edu29(rep)


def test_motion_without_history_reports_no_default(tmp_path):
    """MOTION <out> = 1 (transfer functions only) reads no history: nothing is written or reported."""
    ui = _ui(tmp_path, TINY + "ANALYS,0,0,0,0,1,0,0,0,0,0,0\nMOTION,0,1,0,0,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1\n"
                              "AOPT,0,0,0,1,1,1,0,0,1,0,1,0,0,0\nAFWRITE")
    rep = ui.session["check_report"]
    assert rep.errors() == [] and not _edu29(rep)
    assert decks.read(tmp_path / "m" / "m.mot", "MOTION")["thfile"] == ""
    assert not any("THFILE blank" in t for t in _infos(ui))


# ====================================================================== EQUAKE
def test_equake_runs_with_every_input_blank(tmp_path):
    ui = _ui(tmp_path, "MDL,e,e\nAOPT,1,0,0,0,0,0,0,0,0,0,0,0,0,0\nCHECK\nAFWRITE\nLIBRARY,DEFAULTS")
    rep = ui.session["check_report"]
    assert rep.errors() == [], rep.format()
    lines = _edu29(rep, "EQUAKE")
    assert lines == [
        "RSIN blank: the built-in RG 1.60 H spectrum @rg160h_030g.rsi (0.30 g, 5 %, 27 frequencies) is the target "
        "of spectrum 1",
        "RSOUT 1 blank: the output spectrum is written to e_eq1.rso (model folder)",
        "ACCOUT 1 blank: the generated record is written to e_eq1.acc (.vel, .dis, .psd, .fft next to it; model "
        "folder)"]
    assert [f"  EQUAKE: {x}" for x in lines] == [t for t in _infos(ui) if t.startswith("  EQUAKE: ")]
    d = decks.read(tmp_path / "e" / "e.equ", "EQUAKE")
    assert int(d["nrfreq"]) == 27
    assert [list(r.values()) for r in d.rows("spectra")] == [[1, "@rg160h_030g.rsi", "e_eq1.rso", "", "e_eq1.acc", ""]]
    _run(ui, "EQUAKE")
    mdir = tmp_path / "e"
    for ext in (".acc", ".rso", ".vel", ".dis", ".psd"):
        assert (mdir / f"e_eq1{ext}").is_file(), ext
    out = (mdir / "e_EQUAKE.out").read_text(encoding="utf-8")
    assert lines[0] in out and "built-in input @rg160h_030g.rsi: RG 1.60 horizontal design spectrum" in out
    a, dt = textfiles.read_history(mdir / "e_eq1.acc")
    assert dt == 0.005 and 0.2 < np.max(np.abs(a)) < 0.45                  # matched to the 0.30 g spectrum
    assert not (LIB.LIBRARY_DIR / "e_eq1.acc").exists()


def test_equake_acs_mode_and_read_only_library(tmp_path):
    ui = _ui(tmp_path, "MDL,e,e\nAOPT,1,0,0,0,0,0,0,0,0,0,0,0,0,0\nEDUOPT,DEFAULTS,OFF")
    rep = run_check(ui.model, dirs=[tmp_path / "e"])
    assert rep.numbers("Error", "EQUAKE") == {84}
    ui = _ui(tmp_path, "MDL,e,e\nAOPT,1,0,0,0,0,0,0,0,0,0,0,0,0,0\nRSOUT,1,@rg160h_030g.rsi\nACCOUT,1,@x.acc")
    rep = run_check(ui.model, dirs=[tmp_path / "e"])
    assert rep.numbers("Error", "EQUAKE") == {86, 87}
    assert all("read-only" in m.detail for m in rep.errors("EQUAKE"))
    ui = _ui(tmp_path, "MDL,e,e\nAOPT,1,0,0,0,0,0,0,0,0,0,0,0,0,0\nRSIN,2,@rg160v_1g.rsi")   # spectrum 2 only
    files, uses = DEF.equake_files(ui.model, 0)
    assert files["RSIN"] == {2: "@rg160v_1g.rsi"} and [u.item for u in uses] == ["RSOUT 2", "ACCOUT 2"]


# ====================================================================== SOIL
SOIL_SITE = """MDL,s,s
L,1,2,0.12,600,300,0.05,0.05
L,2,2,0.13,800,400,0.05,0.05
L,3,1,0.14,4000,2000,0.02,0.02
TOPL,1,1,2,2,2
SITE,0,1,0,20,3,1,0,1,2048,1,0,0.005,4096,1
SOIL,4001,32.2,0,1,1,3,0.65,1,0
SACC,1,2,0
AOPT,0,1,0,0,0,0,0,0,0,0,0,0,0,0
"""


def _soil_deck(path):
    from sassi.modules.soil import read_deck
    return read_deck(path)


def test_soil_runs_with_library_curves_without_inp(tmp_path):
    prof = "SPRO,1,1,Sand\nSPRO,2,1,Sand\nSPRO,3,2,Clay\nSPRO,4,2,Clay\nSPRO,5,2,Clay\nSPRO,6,3\n"
    ui = _ui(tmp_path / "a", SOIL_SITE + prof + "CHECK\nAFWRITE")
    rep = ui.session["check_report"]
    assert rep.errors() == [], rep.format()
    got = _edu29(rep, "SOIL")
    assert any(x.startswith("DYNP Sand not defined in the model: the built-in curve Sand") for x in got)
    assert any(x.startswith("DYNP Clay not defined in the model: the built-in curve Clay") for x in got)
    assert ("SOIL <header> = 0 with the built-in record @rg160h_030g.acc: its first line is the time step, SOIL "
            "skips it (read as <header> = 1)") in got
    _run(ui, "SOIL")
    sdir = tmp_path / "a" / "s"
    out = (sdir / "s_SOIL.out").read_text(encoding="utf-8")
    assert "line 1 is the time step; SOIL skips it" in out and "Values read        : 4001 after 1 header" in out
    # the same analysis with the curves read by INP and the record given explicitly (header 1): same decks
    # (but the file name) and the same strain-compatible properties
    shutil.copyfile(EXAMPLES / "rg160h_030g.acc", tmp_path / "rec.acc")
    ui2 = _ui(tmp_path / "b", SOIL_SITE.replace("SOIL,4001,32.2,0,", "SOIL,4001,32.2,1,") +
              f"INP,@dynp_library.pre\nTHFILE,{tmp_path / 'rec.acc'}\n" + prof + "AFWRITE")
    assert not _errors(ui2) and not _edu29(ui2.session["check_report"])
    d1, d2 = _soil_deck(sdir / "s.soi"), _soil_deck(tmp_path / "b" / "s" / "s.soi")
    used = [r for r in d2.rows("dynp") if r["label"] in ("Sand", "Clay")]
    assert d1.rows("dynp") == used and d1.rows("profile") == d2.rows("profile")
    _run(ui2, "SOIL")
    assert (sdir / "FILE88").read_text() == (tmp_path / "b" / "s" / "FILE88").read_text()
    np.testing.assert_array_equal(textfiles.read_history(sdir / "ACC001.TH")[0],
                                  textfiles.read_history(tmp_path / "b" / "s" / "ACC001.TH")[0])


def test_soil_default_profile(tmp_path):
    ui = _ui(tmp_path, SOIL_SITE + "AFWRITE\nLIBRARY,DEFAULTS")
    rep = ui.session["check_report"]
    assert rep.errors() == [], rep.format()
    line = ("SPRO blank: the default SOIL profile is used -- sublayers 1-5 = the TOPL layers 1, 1, 2, 2, 2 with the "
            "built-in curves Sand (Vs < 2493.4 ft/s) and Rock (Vs >= 2493.4 ft/s): Sand 1-5; sublayer 6 = the "
            "half-space L 3 (linear)")
    assert line in _edu29(rep, "SOIL") and f"  SOIL: {line}" in _infos(ui)
    d = _soil_deck(tmp_path / "s" / "s.soi")
    assert [(r["layer"], r["dynprop"]) for r in d.rows("profile")] == \
        [(1, "Sand"), (2, "Sand"), (3, "Sand"), (4, "Sand"), (5, "Sand"), (6, "")]
    _run(ui, "SOIL")
    assert line in (tmp_path / "s" / "s_SOIL.out").read_text(encoding="utf-8")
    # SI units: Vs >= 760 m/s is Rock; the FILE73 (PIN) curve order follows the deck
    ui = _ui(tmp_path / "si", SOIL_SITE.replace("32.2", "9.81"))
    prof, use, _ = DEF.soil_profile(ui.model)
    assert [rec.dynprop for _, rec in prof] == ["Sand"] * 5 + [""] and "760 m/s" in use.text
    ui.execute("L,2,2,0.13,1600,800,0.05,0.05")
    prof, use, _ = DEF.soil_profile(ui.model)
    assert [rec.dynprop for _, rec in prof] == ["Sand", "Sand", "Rock", "Rock", "Rock", ""]
    from sassi.prep.commands.nlsoil_cmds import file73_labels
    assert file73_labels(ui.model) == ["Sand", "Rock"]


def test_given_model_curve_wins_and_case_matters(tmp_path):
    ui = _ui(tmp_path, SOIL_SITE + "DYNP,1,0.0001,1,0.0001,1,Sand\nDYNP,2,1,0.5,1,10,Sand\nSPRO,1,1,Sand\nSPRO,2,3")
    rep = run_check(ui.model, dirs=[tmp_path / "s"])
    assert rep.errors() == [] and not [x for x in _edu29(rep) if x.startswith("DYNP")]
    rows = DEF.dynp_table(ui.model, ["Sand"])
    assert [r[1:4] for r in rows] == [(1, 0.0001, 1.0), (2, 1.0, 0.5)]
    ui = _ui(tmp_path, SOIL_SITE + "SPRO,1,1,sand\nSPRO,2,3")
    rep = run_check(ui.model, dirs=[tmp_path / "s"])
    assert {97, 99} <= rep.numbers("Error") and "built-in curves: Clay, Sand, Rock" in rep.errors("SOIL")[0].detail


# ====================================================================== LIBRARY
def test_library_command(tmp_path):
    ui = _ui(tmp_path, "LIBRARY\nLIBRARY,SPECTRUM\nLIBRARY,@ricker_5hz.th\nLIBRARY,DEFAULTS")
    info = _infos(ui)
    for e in LIB.CATALOGUE:
        assert any(t.strip().startswith(e.ref) for t in info), e.ref
    spectra = [t for t in info if t.startswith("  @") and " spectrum " in t]
    assert len(spectra) == 6                     # 3 in LIBRARY and 3 in LIBRARY,SPECTRUM
    assert "LIBRARY @ricker_5hz.th: 5 Hz Ricker wavelet, peak 1 at t = 0.5 s, 2 s" in info
    assert any(t.startswith("  file      ") and t.endswith("ricker_5hz.th") for t in info)
    assert not ui.execute("LIBRARY,NOPE") and not ui.execute("LIBRARY,@nope.acc")
    ui = _ui(tmp_path, "EDUOPT,DEFAULTS,OFF\nLIBRARY,DEFAULTS")
    assert any("EDUOPT,DEFAULTS,OFF -- no defaults" in t for t in _infos(ui))


# ====================================================================== GUI
def test_dialog_values_placeholders_and_spro_pages(tmp_path):
    ui = _ui(tmp_path, SOIL_SITE)
    v = D.values(ui)
    ctx = v["context"]
    assert {e["name"] for e in ctx["library"]} == set(LIB.names())
    assert ctx["dynp_library"] == ["Clay", "Sand", "Rock"]
    idf = ctx["input_defaults"]
    assert idf["thfile"]["0"] == "built-in: RG 1.60, 0.30 g (@rg160h_030g.acc)"
    assert idf["thfile"]["1"] == "built-in: 5 Hz Ricker pulse (@ricker_5hz.th)"
    assert idf["rsout"] == "s_eq{i}.rso" and idf["accout"] == "s_eq{i}.acc"
    assert v["indexed_defaults_by_key"]["SPRO"]["3"] == {"layer": 3, "prop": 2, "dynprop": "Sand"}
    assert ctx["spro_default"] == "default: Sand 1-5; 6 = half-space L 3 (EDU-29)"
    # the field specs: the defaults and the Library kinds
    fields = {}
    for tab in D.analysis_tabs():
        for g in tab["groups"]:
            for it in g["items"]:
                fields.setdefault(it.get("path"), []).append(it)
    assert {it.get("default_ph") for it in fields["$THFILE"]} == {"thfile", "soil_thfile"}
    assert fields["RSIN[spec].file"][0]["library"] == ["spectrum"] and fields["TPSD[spec].file"][0]["library"] == ["psd"]
    # changing one SPRO page stores the default profile with it (Shown = used)
    lines, notes = D.commit_commands(ui, {"indexed": {"SPRO": {"2": {"layer": 2, "prop": 1, "dynprop": "Clay"}}}})
    assert lines == ["SPRO,1,1,Sand", "SPRO,2,1,Clay", "SPRO,3,2,Sand", "SPRO,4,2,Sand", "SPRO,5,2,Sand", "SPRO,6,3"]
    assert any("default SOIL profile" in n for n in notes)
    ui.execute("EDUOPT,DEFAULTS,OFF")
    v = D.values(ui)
    assert v["context"]["input_defaults"] is None and not v["indexed_defaults_by_key"]


_NODE = r"""
const vm = require("vm"), fs = require("fs");
const ctx = {SASSI: {el: () => ({})}, console};
vm.createContext(ctx);
vm.runInContext(fs.readFileSync(process.argv[2], "utf8"), ctx);
const data = JSON.parse(fs.readFileSync(process.argv[3], "utf8"));
const F = new ctx.SASSI.D.OptionsForm({tabs: []}, data.values, "ANALYSIS");
const out = {};
for (const [name, it, spec, edits] of data.cases) {
  F.sel.spec = spec;
  for (const [path, val] of edits) F.set(path, val);
  out[name] = F.defaultPlaceholder(it);
}
out.spro3 = F.absentEntry("SPRO", "3");
console.log(JSON.stringify(out));
"""


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_browser_placeholders_follow_the_policy(tmp_path):
    """dialogs.js defaultPlaceholder (the browser) applied to the Python values: the placeholder of every blank file
    is the default sassi/prep/defaults.py applies for the values being edited."""
    ui = _ui(tmp_path, SOIL_SITE)
    th = {"default_ph": "thfile"}
    cases = [
        ["seismic", th, 1, []],
        ["vibration", th, 1, [["ANALYS.type", 1]]],
        ["pairs", th, 1, [["MOTION.fopt", 1]]],
        ["rsin1", {"default_ph": "rsin"}, 1, []],
        ["rsin2", {"default_ph": "rsin", "placeholder": "own"}, 2, []],
        ["rsout1", {"default_ph": "rsout"}, 1, []],
        ["accout2_none", {"default_ph": "accout"}, 2, []],
        ["accout2", {"default_ph": "accout"}, 2, [["RSIN[spec].file", "@rg160v_1g.rsi"]]],
        ["rsin1_after", {"default_ph": "rsin"}, 1, []],
        ["soil_dt", {"default_ph": "soil_thfile"}, 1, [["SITE.delt", 0.01]]],
    ]
    data = {"values": D.values(ui), "cases": cases}
    (tmp_path / "data.json").write_text(json.dumps(data), encoding="utf-8")
    (tmp_path / "t.js").write_text(_NODE, encoding="utf-8")
    r = subprocess.run(["node", str(tmp_path / "t.js"), str(ROOT / "sassi" / "ui" / "static" / "dialogs.js"),
                        str(tmp_path / "data.json")], capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr
    got = json.loads(r.stdout)
    assert got["seismic"] == "built-in: RG 1.60, 0.30 g (@rg160h_030g.acc)"
    assert got["vibration"] == "built-in: 5 Hz Ricker pulse (@ricker_5hz.th)"
    assert got["pairs"].startswith("no built-in default")
    assert got["rsin1"] == "built-in: RG 1.60 H, 0.30 g (@rg160h_030g.rsi)"
    assert got["rsin2"] == "own" and got["accout2_none"] is None
    assert got["rsout1"] == "default: s_eq1.rso" and got["accout2"] == "default: s_eq2.acc"
    assert got["rsin1_after"] is None             # RSIN 2 given in the dialog: no default for RSIN 1 (D-W5-07)
    assert got["soil_dt"] == "no built-in default: it needs time step 0.005 s"
    assert got["spro3"] == {"layer": 3, "prop": 2, "dynprop": "Sand"}


def test_gui_file_api_reads_library_files(tmp_path):
    from sassi.ui.api import GuiSession
    S = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    st, f = S.route("GET", "/api/file", {"name": ["@rg160h_030g.rsi"]}, {})
    assert st == 200 and f["readonly"] and f["text"].startswith("0.100000")
    st, info = S.route("GET", "/api/fileinfo", {"name": ["@rg160h_030g.acc"]}, {})
    assert st == 200 and info["kind"] == "history"
    st, _ = S.route("POST", "/api/file", {}, {"name": "@rg160h_030g.rsi", "text": "x"})
    assert st == 403                                                     # read-only
    st, op = S.route("POST", "/api/file", {}, {"name": "@rg160h_030g.rsi"})        # File > Open, the Edit button
    assert st == 200 and op["path"] == "@rg160h_030g.rsi" and op["readonly"] and not op["created"]
    assert S.route("GET", "/api/file", {"name": ["@nope.rsi"]}, {})[0] == 404
    st, dy = S.route("GET", "/api/dynp", {}, {})
    assert st == 200 and sorted(dy["library"]) == ["Clay", "Rock", "Sand"] and dy["properties"] == {}
    assert filecmp.cmp(LIB.LIBRARY_DIR / "rg160h_030g.rsi", EXAMPLES / "rg160h_030g.rsi", shallow=False)


def test_soil_property_plot_of_a_library_curve(tmp_path):
    from sassi.plotting.state import soil_property_curves
    ui = _ui(tmp_path, "")
    c = soil_property_curves(ui.model, "Sand")
    assert c["builtin"] and c["name"] == "Sand (built-in, @dynp_library.pre)" and c["g"][6] == 0.37
    assert ui.execute("SOILPROPPLOT,Clay")
