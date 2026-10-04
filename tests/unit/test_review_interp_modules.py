"""Adversarial review tests of work package A2 (interpreter part 2: option commands, CHECK, AFWRITE,
RUN commands, model checks and conditioning).

These tests are written independently of the implementer's tests: the expected values are derived
by hand from the requirements (sections 1.9, 2.6, 3.3, 3.4.F-R, 4.0.1, 5.4, 5.5, decisions D-CHK-*,
D-AFW-*, D-RUN-*, D-HOU-01, D-CNV-10, D-MOT-10) and from spec 09 / spec 11, not from the
implementation.  Tests that fail document defects found in the review.
"""
from __future__ import annotations

import math
import shutil
import warnings
from pathlib import Path

import numpy as np
import pytest

from sassi.io import decks
from sassi.prep import Interpreter, Kind
from sassi.prep import options as O
from sassi.prep.afwrite import afwrite
from sassi.prep.check import MODEL, CheckOptions, ModelView, run_check


# ======================================================================================
# helpers
# ======================================================================================
def _ui(text: str, cwd=None) -> Interpreter:
    ui = Interpreter(cwd=cwd)
    ui.run_text(text)
    return ui


def _deck(mdir: Path, model: str, module: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return decks.read(decks.deck_path(mdir, model, module), module)


def _history(path: Path, n: int, dt: float = 0.01) -> None:
    path.write_text(f"{dt}\n" + "\n".join(f"{0.1 * math.sin(0.3 * k):.6f}" for k in range(n)) + "\n")


#: an excavated 1x1x5 block (bottom nodes 1-4 at z = -5, top nodes 5-8 at grade), all
#: interaction nodes, on two soil layers; a frequency set and the SITE/POINT options
BLOCK = """
N,1,0,0,-5
N,2,1,0,-5
N,3,1,1,-5
N,4,0,1,-5
N,5,0,0,0
N,6,1,0,0
N,7,1,1,0
N,8,0,1,0
GROUP,1,SOLID
E,1,1,2,3,4,5,6,7,8
L,1,5,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
TOPL,1
INT,1,8,1,1
FREQ,1,1,2,3
SITE,0,1,0,20,2,1,0,1,,1,0,{delt},{nft},1
POINT,0,1,1
"""


def _block(tmp_path: Path, extra: str = "", aopt: str = "AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0",
           delt: float = 0.01, nft: int = 4096, name: str = "blk"):
    mdir = tmp_path / name
    mdir.mkdir(exist_ok=True)
    ui = _ui(f"MDL,{name},{mdir}\n" + BLOCK.format(delt=delt, nft=nft) + extra + "\n" + aopt + "\n",
             cwd=tmp_path)
    return ui, mdir


# ======================================================================================
# AFWRITE robustness: a CHECK error in one module must not crash the other decks
# ======================================================================================
def test_afwrite_undefined_topl_layer_blocks_site_without_crashing(tmp_path):
    """A typo in TOPL (layer 3 undefined) is Error 19 under SITE (spec 11 section 5.2).  SITE gets no
    deck (D-AFW-01); AFWRITE itself must not fail with an internal error, and the modules without
    errors (POINT, ANALYS ...) still get their decks."""
    ui, mdir = _block(tmp_path, extra="TOPL,0\nTOPL,1,3\n")
    ok = ui.execute("AFWRITE")
    errors = ui.sink.texts(Kind.ERROR)
    assert ok and not any("internal error" in e for e in errors), errors
    assert not (mdir / "blk.sit").exists()
    assert "Error 19 : Soil Layer 3 Is not Defined" in (mdir / "blk.err").read_text()
    assert (mdir / "blk.poi").exists()


def test_afwrite_eduopt_soilcutoff_with_fortran_exponent(tmp_path):
    """Rule L6 / D-PAR-07: numbers accept Fortran ``D`` exponents.  EDUOPT accepts ``2.5D1`` for the
    float key SOILCUTOFF, so AFWRITE must write soilcutoff = 25 instead of failing."""
    mdir = tmp_path / "s"
    mdir.mkdir()
    _history(mdir / "H1.acc", 300)
    ui = _ui(f"""MDL,s,{mdir}
L,1,5,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
THFILE,H1.acc
SPRO,1,1,S
SPRO,2,2,S
DYNP,1,0.001,1,0.001,1,S
SOIL,100,32.2,0,1,1,8,0.65,1,0
SITE,0,1,0,20,2,1,0,1,,1,0,0.01,4096,1
EDUOPT,SOILCUTOFF,2.5D1
AOPT,0,1,0,0,0,0,0,0,0,0,0,0,0,0
""", cwd=tmp_path)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    assert ui.execute("AFWRITE"), ui.sink.texts(Kind.ERROR)
    from sassi.modules.soil import read_deck
    assert read_deck(mdir / "s.soi")["soilcutoff"] == pytest.approx(25.0)


# ======================================================================================
# Frequency grid and NFFT rounding (section 4.0.1, D-CNV-10, UT-16)
# ======================================================================================
@pytest.mark.parametrize("nft, expect", [(3000, 2048), (3072, 4096), (6144, 8192), (5000, 4096),
                                         (4096, 4096), (12289, 16384)])
def test_nfft_nearest_power_of_two_and_df_in_every_deck(tmp_path, nft, expect):
    """Nearest power of 2, ties up; df = 1/(delt NFFT_written) copied into every deck (shared)."""
    ui, mdir = _block(tmp_path, delt=0.01, nft=nft)
    assert ui.execute("AFWRITE")
    df = 1.0 / (0.01 * expect)
    for mod in ("SITE", "POINT", "HOUSE", "ANALYS"):
        d = _deck(mdir, "blk", mod)
        assert d["df"] == pytest.approx(df, rel=1e-12), mod
        if "nft" in d.params:
            assert d["nft"] == expect, mod
    assert _deck(mdir, "blk", "SITE")["freq2"] == expect // 2           # blank <freq2> = NFFT/2
    err = (mdir / "blk.err").read_text()
    assert ("Warning 9 : Number of Values for Fourier Transform Is Not Power of 2" in err) == (nft != expect)


def test_nfftround_up(tmp_path):
    ui, mdir = _block(tmp_path, extra="EDUOPT,NFFTROUND,UP\n", nft=3000)
    assert ui.execute("AFWRITE")
    assert _deck(mdir, "blk", "SITE")["nft"] == 4096


def test_ut16_harmonic_grid_with_frequency_step(tmp_path):
    """UT-16: df = 0.9765625 Hz with numbers {1,3,...,18} -> f = n df; numbers sorted in the deck."""
    nums = [18, 1, 3, 5, 7, 9, 11, 13, 15, 16]
    ui, mdir = _block(tmp_path, extra="FREQ,1,0\nFREQ,1," + ",".join(map(str, nums)) + "\n"
                      "SITE,0,1,0.9765625,20,2,1,0,1,,1,0,0.01,4096,1\n")
    assert ui.execute("AFWRITE")
    d = _deck(mdir, "blk", "SITE")
    fn = [r["number"] for r in d.rows("freqs")]
    assert fn == sorted(nums)
    f = np.array(fn) * d["df"]
    ref = [0.977, 2.930, 4.883, 6.836, 8.789, 10.742, 12.695, 14.648, 15.625, 17.578]
    assert np.allclose(np.round(f, 3), ref, atol=1e-12)


def test_duplicate_frequency_numbers_are_an_error(tmp_path):
    """Section 4.0.1: duplicate SSI frequency numbers are an error (EDU-23, G-01)."""
    ui, mdir = _block(tmp_path, extra="FREQ,1,4,4\n")
    rep = run_check(ui.model)
    assert rep.has("Error", "EDU-23", "SITE")
    assert rep.blocked("SITE")


# ======================================================================================
# HOUSE deck: global coordinates, ETYPE resolution, gap / unused nodes (D-AFW-06, D-HOU-01)
# ======================================================================================
def test_house_deck_global_coordinates_of_local_nodes(tmp_path):
    """UT-12 values: LOC,1,0,10,0,0,45,0,0 then N,100,5,0,0 in CSYS 1 -> (13.5355, 3.5355, 0)."""
    ui, mdir = _block(tmp_path, extra="LOC,1,0,10,0,0,45,0,0\nCSYS,1\nN,100,5,0,0\nCSYS,0\n"
                                      "GROUP,2,BEAMS\nE,1,5,100,6\nM,1,4e5,0.25,0.15,0.05,0.05\n"
                                      "R,1,1,0.8,0.8,0.2,0.1,0.1\nD,100,100,1,1,ALL\n")
    assert ui.execute("AFWRITE")
    rows = {r["id"]: r for r in _deck(mdir, "blk", "HOUSE").rows("nodes")}
    r = rows[100]
    assert (r["x"], r["y"], r["z"]) == pytest.approx((10 + 5 / math.sqrt(2), 5 / math.sqrt(2), 0.0), abs=1e-9)


def test_house_deck_etype_resolution_with_ground_elevation(tmp_path):
    """D-HOU-01 with gelev = 2: an element touching grade from below is excavated (2), one that
    crosses grade is structure (1); explicit ETYPE is kept."""
    text = """
GROUNDELEV,2
L,1,5,0.12,1500,800,0.05,0.05
N,11,5,0,1
N,12,6,0,1
N,13,6,1,1
N,14,5,1,1
N,15,5,0,2
N,16,6,0,2
N,17,6,1,2
N,18,5,1,2
N,21,8,0,1.5
N,22,9,0,1.5
N,23,9,1,1.5
N,24,8,1,1.5
N,25,8,0,2.5
N,26,9,0,2.5
N,27,9,1,2.5
N,28,8,1,2.5
GROUP,2,SOLID
E,1,11,12,13,14,15,16,17,18
E,2,21,22,23,24,25,26,27,28
M,1,4e5,0.25,0.15,0.05,0.05
GROUP,3,SOLID
E,1,21,22,23,24,25,26,27,28
ETYPE,1,1,1,2
AOPT,0,0,0,0,0,1,0,0,0,0,0,0,0,0
"""
    mdir = tmp_path / "g"
    mdir.mkdir()
    ui = _ui(f"MDL,g,{mdir}\n" + text, cwd=tmp_path)
    v = ModelView(ui.model)
    by = {(r.group, r.id): r.etype for r in v.elems}
    assert by[(2, 1)] == 2 and by[(2, 2)] == 1 and by[(3, 1)] == 2
    assert ui.execute("AFWRITE"), ui.sink.texts(Kind.ERROR)
    el = {(r["group"], r["id"]): r["etype"] for r in _deck(mdir, "g", "HOUSE").rows("elements")}
    assert el == {(2, 1): 2, (2, 2): 1, (3, 1): 2}


def test_house_deck_fixes_gap_and_unused_nodes_but_not_unused_interaction_nodes(tmp_path):
    ui, mdir = _block(tmp_path, extra="N,12,3,3,0\nN,14,4,4,0\nINT,14,14,1,1\n")
    assert ui.execute("AFWRITE")
    rows = {r["id"]: r for r in _deck(mdir, "blk", "HOUSE").rows("nodes")}
    fix = lambda n: [rows[n][k] for k in ("fx", "fy", "fz", "frx", "fry", "frz")]   # noqa: E731
    assert sorted(rows) == list(range(1, 15))                    # gaps 9-11, 13 written (Warning 1)
    for n in (9, 10, 11, 13):
        assert fix(n) == [1] * 6
    assert fix(12) == [1] * 6                                    # unused node (Warning 4)
    assert fix(14) == [0] * 6                                    # unused interaction node: kept
    assert fix(1) == [0] * 6


# ======================================================================================
# MOTION deck: duplicate output requests merged with OR (D-MOT-10)
# ======================================================================================
def test_nout_duplicates_merged_with_or_and_reported_as_edu16(tmp_path):
    ui, mdir = _block(tmp_path, aopt="AOPT,0,0,0,0,0,0,0,0,0,0,1,0,0,0")
    _history(mdir / "H1.acc", 400)
    ui.run_text("THFILE,H1.acc\nDAMP,0.05\n"
                "NOUT,1,1,0,0,0,0,0,1-3\n"
                "NOUT,1,0,1,0,0,0,1,2\n"
                "NOUT,2,0,0,0,0,1,0,2\n")
    rep = run_check(ui.model, dirs=[mdir])
    assert not rep.has("Error", 66) and rep.has("Warning", "EDU-16", "MOTION")
    assert ui.execute("AFWRITE")
    rows = [(r["node"], r["dir"], [r[f"c{k}"] for k in range(1, 7)]) for r in _deck(mdir, "blk", "MOTION").rows("nout")]
    got = {(n, d): c for n, d, c in rows}
    assert len(rows) == len(got) == 4
    assert got[(2, 1)] == [1, 1, 0, 0, 0, 1]
    assert got[(1, 1)] == [1, 0, 0, 0, 0, 0] and got[(3, 1)] == [1, 0, 0, 0, 0, 0]
    assert got[(2, 2)] == [0, 0, 0, 0, 1, 0]


# ======================================================================================
# Shared variables: one storage location, copied into every deck (spec 07 section 3)
# ======================================================================================
def test_history_data_copied_into_motion_stress_reldisp(tmp_path):
    ui, mdir = _block(tmp_path, aopt="AOPT,0,0,0,0,0,0,0,0,0,0,1,1,1,0")
    _history(mdir / "H1.acc", 400)
    ui.run_text("THFILE,H1.acc\nTHTIT,my motion\nDAMP,0.05\nNOUT,1,1,0,0,0,0,0,1\n"
                "MOTION,0,0,0,0,0,0.1,100,301,0,0.35,2,300,0,0,0,1,0,0,1\n"
                "EOUT,1,1,1,1,1,1,0,0,0,0,0,0,1,1\nRDND,1,1,0,0,0,0,0\n")
    assert ui.execute("AFWRITE"), ui.sink.texts(Kind.ERROR)
    for mod in ("MOTION", "STRESS", "RELDISP"):
        d = _deck(mdir, "blk", mod)
        assert (d["thtit"], d["mult"], d["max"], d["rec1"], d["rec2"], d["fopt"]) == \
            ("my motion", 0.0, 0.35, 2, 300, 0), mod
        assert Path(d["thfile"]).name == "H1.acc" and d["delt"] == 0.01 and d["nft"] == 4096, mod


def test_soil_control_layer_shared_with_site_and_soilx_override(tmp_path):
    mdir = tmp_path / "s"
    mdir.mkdir()
    _history(mdir / "H1.acc", 300)
    base = f"""MDL,s,{mdir}
L,1,5,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
THFILE,H1.acc
SPRO,1,1,S
SPRO,2,1,S
SPRO,3,2,S
DYNP,1,0.001,1,0.001,1,S
SOIL,100,32.2,0,1,1,8,0.65,1,0
SITE,0,1,0,20,2,1,0,1,,2,0,0.01,4096,1
AOPT,0,1,0,0,0,0,0,0,0,0,0,0,0,0
"""
    from sassi.modules.soil import read_deck
    ui = _ui(base, cwd=tmp_path)
    assert ui.execute("AFWRITE")
    assert read_deck(mdir / "s.soi")["cl"] == 2                 # SITE <cl>
    ui.execute("SOILX,0,0,0.1,3")
    assert ui.execute("AFWRITE")
    assert read_deck(mdir / "s.soi")["cl"] == 3                 # SOIL-only override


def test_gravity_error_reported_under_every_module_that_uses_it(tmp_path):
    """D-CHK-01: an error on a shared variable is reported under every module using it."""
    ui, mdir = _block(tmp_path, extra="GRAVITY,-1\n", aopt="AOPT,0,0,0,1,1,1,0,1,1,0,1,1,1,0")
    rep = run_check(ui.model, dirs=[mdir])
    for mod in ("SITE", "HOUSE", "FORCE", "ANALYS", "MOTION", "STRESS", "RELDISP"):
        assert rep.has("Error", 1, mod), mod


def test_undefined_frequency_set_reported_under_its_users(tmp_path):
    ui, mdir = _block(tmp_path, extra="SITE,0,1,0,20,2,1,0,1,,1,0,0.01,4096,7\n",
                      aopt="AOPT,0,0,0,1,1,0,0,1,1,0,0,0,0,0")
    rep = run_check(ui.model, dirs=[mdir])
    for mod in ("SITE", "POINT", "FORCE", "ANALYS"):
        assert rep.has("Error", 44, mod), mod
        assert any("Frequency Set 7 Is Not Defined" in m.text for m in rep.errors(mod))


# ======================================================================================
# CHECK report: .err format and break count (section 5.5, Options > Check)
# ======================================================================================
def test_err_file_format_headers_in_run_order(tmp_path):
    ui, mdir = _block(tmp_path, extra="POINT,0,1,0\n", aopt="AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0")
    ui.execute("CHECK")
    lines = (mdir / "blk.err").read_text().splitlines()
    assert lines[0] == "CHECK: Errors and Warning for - blk"
    heads = [ln[len("Errors and Warnings for "):] for ln in lines if ln.startswith("Errors and Warnings for ")]
    assert heads == ["MODEL", "SITE", "POINT", "HOUSE", "ANALYS"]
    i = lines.index("Errors and Warnings for POINT")
    assert lines[i + 1].startswith("Error 57 : Illegal Radius of Central Zone")
    assert "POINT: 1 errors, 0 warnings" in lines


def test_break_check_limits_errors_and_warnings_per_module(tmp_path):
    """Options > Check "Break Check at N": at most N messages of each type (errors and warnings)
    printed per module (spec 05a section 3: "Most messages of each type (errors and warnings)"),
    totals still counted (D-CHK-02)."""
    text = """
N,1,0,0,0
N,2,1,0,0
N,3,1,1,0
N,4,0,1,0
N,9,5,5,5
N,10,6,5,5
N,11,7,5,5
GROUP,1,SHELL
E,1,1,2,3,40
E,2,1,2,3,41
E,3,1,2,3,42
MSET,1,1,1,5
MSET,2,2,1,6
MSET,3,3,1,7
"""
    ui = _ui(text)
    rep = run_check(ui.model, modules=["HOUSE"])
    n_err = len(rep.errors(MODEL))
    assert n_err >= 6                                    # three Error 41 and three Error 13
    out = rep.format(CheckOptions(break_at=2))
    sec = out.split("Errors and Warnings for MODEL")[1].split("Errors and Warnings for")[0]
    printed_errors = [ln for ln in sec.splitlines() if ln.startswith("Error ")]
    printed_warnings = [ln for ln in sec.splitlines() if ln.startswith("Warning ")]
    assert len(printed_errors) == 2
    assert len(printed_warnings) <= 2
    assert f"MODEL: {n_err} errors" in sec


# ======================================================================================
# Model checks (spec 09 section 3, D-CHK-06)
# ======================================================================================
def test_error_124_needs_all_translations_fixed_fixedint_lists_any(tmp_path):
    ui, mdir = _block(tmp_path, extra="D,1,1,1,1,UX\nD,2,2,1,1,DISP\n")
    rep = run_check(ui.model, modules=["HOUSE"])
    e124 = [m for m in rep.errors(MODEL) if m.number == 124]
    assert len(e124) == 1 and "Node 2 is a fixed interaction node" in e124[0].text
    ui.sink.clear()
    ui.execute("FIXEDINT")
    infos = "\n".join(ui.sink.texts(Kind.INFO))
    assert "2 fixed interaction nodes" in infos and "node 1: UX fixed" in infos


def test_excstrchk_interior_node_shared_with_beam(tmp_path):
    """A 2x2x2 excavated cube has one interior node (the centre); a structural beam attached to it is
    reported by EXCSTRCHK; a beam attached to a boundary node is not."""
    lines = []
    nid = {}
    k = 0
    for iz in range(3):
        for iy in range(3):
            for ix in range(3):
                k += 1
                nid[(ix, iy, iz)] = k
                lines.append(f"N,{k},{ix},{iy},{iz - 2}")
    lines.append("GROUP,1,SOLID")
    e = 0
    for iz in range(2):
        for iy in range(2):
            for ix in range(2):
                e += 1
                c = [nid[(ix + a, iy + b, iz + cz)] for cz in (0, 1) for (a, b) in ((0, 0), (1, 0), (1, 1), (0, 1))]
                lines.append("E,{},{}".format(e, ",".join(map(str, c))))
    lines.append("ETYPE,1,8,1,2")
    lines += ["M,1,4e5,0.25,0.15,0.05,0.05", "R,1,1,0.8,0.8,0.2,0.1,0.1", "N,100,1,1,5", "N,101,0,5,5",
              "GROUP,2,BEAMS", f"E,1,{nid[(1, 1, 1)]},100,101"]
    centre = nid[(1, 1, 1)]
    ui = _ui("\n".join(lines))
    ui.execute("EXCSTRCHK")
    infos = "\n".join(ui.sink.texts(Kind.INFO))
    assert "1 excavation interior nodes" in infos and f"node {centre}:" in infos
    ui2 = _ui("\n".join(lines[:-1] + [f"E,1,{nid[(1, 1, 2)]},100,101"]))
    ui2.execute("EXCSTRCHK")
    assert "0 excavation interior nodes" in "\n".join(ui2.sink.texts(Kind.INFO))


def _plate(theta_deg: float) -> str:
    """3x3-node shell plate through the X axis, rotated by theta about X (normal (0, -sin, cos))."""
    th = math.radians(theta_deg)
    lines, ids, k = [], {}, 0
    for j in range(3):
        for i in range(3):
            k += 1
            ids[(i, j)] = k
            lines.append(f"N,{k},{i},{j * math.cos(th):.15g},{j * math.sin(th):.15g}")
    lines.append("GROUP,1,SHELL")
    e = 0
    for j in range(2):
        for i in range(2):
            e += 1
            lines.append(f"E,{e},{ids[(i, j)]},{ids[(i + 1, j)]},{ids[(i + 1, j + 1)]},{ids[(i, j + 1)]}")
    lines += ["M,1,3e6,0.25,0.15,0.05,0.05", "D,1,9,1,1,DISP"]
    return "\n".join(lines) + "\n"


def test_fixrot_oblique_plate_springs_follow_the_normal():
    """FIXROT,[Stiff] (spec 09 section 2.8/2.10): oblique shells get a drilling spring to a new
    coincident fixed node with scxx = k nx^2, scyy = k ny^2, sczz = k nz^2; the command is idempotent
    and removes EDU-06."""
    ui = _ui(_plate(40.0))
    assert len([m for m in run_check(ui.model, modules=["HOUSE"]).warnings(MODEL) if m.number == "EDU-06"]) == 9
    ui.execute("FIXROT,8")
    m = ui.model
    th = math.radians(40.0)
    assert len(m.springs) == 1
    sc = next(iter(m.springs.values()))
    assert (sc.scx, sc.scy, sc.scz) == (0.0, 0.0, 0.0)
    assert sc.scxx == pytest.approx(0.0, abs=1e-12)
    assert sc.scyy == pytest.approx(8 * math.sin(th) ** 2, rel=1e-10)
    assert sc.sczz == pytest.approx(8 * math.cos(th) ** 2, rel=1e-10)
    spring_groups = [g for g in m.groups.values() if g.type == 7]
    assert len(spring_groups) == 1 and len(spring_groups[0].elements) == 9
    for el in spring_groups[0].elements.values():
        a, b = el.nodes[:2]
        pa, pb = m.nodes[a], m.nodes[b]
        assert (pa.x, pa.y, pa.z) == pytest.approx((pb.x, pb.y, pb.z), abs=1e-12)
        assert pb.fix == [1] * 6 and pa.fix[3:] == [0, 0, 0]
    assert not any(x.number == "EDU-06" for x in run_check(m, modules=["HOUSE"]).messages)
    n_before = sum(len(g.elements) for g in m.groups.values())
    ui.execute("FIXROT,8")
    assert sum(len(g.elements) for g in m.groups.values()) == n_before


def test_fixrot_plate_in_xz_plane_fixes_roty_only():
    ui = _ui(_plate(90.0))
    ui.execute("FIXROT")
    m = ui.model
    assert not m.springs
    for n in range(1, 10):
        assert m.nodes[n].fix == [1, 1, 1, 0, 1, 0], n


def test_fixsprrot_fixes_unstiffened_dofs_of_spring_only_nodes():
    """Spring-only node with kx only and no mass: UY, UZ and the rotations are fixed, UX stays free;
    a node carrying a translational mass keeps its translations."""
    ui = _ui("N,1,0,0,0\nN,2,0,0,0\nN,3,1,0,0\nN,4,1,0,0\nD,2,2,1,1,ALL\nD,4,4,1,1,ALL\n"
             "SC,1,100,0,0,0,0,0,0\nGROUP,1,SPRING\nE,1,1,2\nE,2,3,4\nRSET,1,2,1,1\nMT,3,0,5,0\n")
    ui.execute("FIXSPRROT")
    m = ui.model
    assert m.nodes[1].fix == [0, 1, 1, 1, 1, 1]
    assert m.nodes[3].fix == [0, 0, 1, 1, 1, 1]


# ======================================================================================
# SITE wave field rules (Errors 51, 54, 55)
# ======================================================================================
def test_wave_ratio_rules(tmp_path):
    ui, mdir = _block(tmp_path, extra="WAVE,1,1,0.3,0.5,0\nWAVE,2,1,0.3,0.5,0\n")
    rep = run_check(ui.model)
    e55 = [m for m in rep.errors("SITE") if m.number == 55]
    assert len(e55) == 1 and "Frequency 1" in e55[0].text          # 0.6 at f1; 1.0 at f2
    ui.execute("WAVE,3,1,0.4,0,0")
    rep = run_check(ui.model)
    assert rep.has("Error", 54, "SITE") and not rep.has("Error", 55, "SITE")
    ui.execute("SITE,0,1,0,20,2,1,1,1,,1,0,0.01,4096,1")              # SH/L family: none defined
    assert run_check(ui.model).has("Error", 51, "SITE")


def test_error_53_names_the_frequency_index():
    """Error 53 "Illegal Value for Frequency <i>": i is 1 for <freq1> and 2 for <freq2> (spec 11)."""
    ui = _ui("L,1,5,0.12,1500,800,0.05,0.05\nL,2,10,0.13,2400,1200,0.02,0.02\nTOPL,1\nFREQ,1,1,2\n"
             "SITE,0,1,0,20,2,1,0,0,-4,1,0,0.01,4096,1\n")
    texts = sorted(m.text for m in run_check(ui.model, modules=["SITE"]).errors("SITE") if m.number == 53)
    assert texts == ["Illegal Value for Frequency 1", "Illegal Value for Frequency 2"]


# ======================================================================================
# MOTION checks (Errors 72, EDU-03, EDU-19)
# ======================================================================================
def test_motion_damping_ratio_range(tmp_path):
    ui, mdir = _block(tmp_path, aopt="AOPT,0,0,0,0,0,0,0,0,0,0,1,0,0,0")
    _history(mdir / "H1.acc", 400)
    ui.run_text("THFILE,H1.acc\nNOUT,1,1,0,0,0,1,0,1\nDAMP,0.05,1.0\n")
    rep = run_check(ui.model, dirs=[mdir])
    assert rep.has("Error", 72, "MOTION")


def test_edu03_rounded_nfft_shorter_than_records(tmp_path):
    """D-CNV-10: NFFT 3000 is written as 2048; a 2500-record history is then truncated: EDU-03 error
    for the modules that read the history; a 2000-record history is fine."""
    ui, mdir = _block(tmp_path, nft=3000, aopt="AOPT,0,0,0,0,0,0,0,0,0,0,1,1,1,0")
    _history(mdir / "H1.acc", 2500)
    ui.run_text("THFILE,H1.acc\nDAMP,0.05\nNOUT,1,1,0,0,0,0,0,1\nEOUT,1,1,1,1,1,1,0,0,0,0,0,0,1,1\n"
                "RDND,1,1,0,0,0,0,0\n")
    rep = run_check(ui.model, dirs=[mdir])
    for mod in ("MOTION", "STRESS", "RELDISP"):
        assert rep.has("Error", "EDU-03", mod), mod
        assert rep.has("Warning", 9, mod), mod
    _history(mdir / "H1.acc", 2000)
    rep = run_check(ui.model, dirs=[mdir])
    assert not any(m.number == "EDU-03" for m in rep.messages)


def test_edu19_quiet_zone(tmp_path):
    """EDU-19: trailing zeros NFFT dt - T_record >= ln(100)/(2 pi f_min beta_min).

    One 5 m layer, Vs = 800: f_soil = 800/20 = 40 Hz; first SSI frequency 1 x df = 1/(0.01 4096)
    = 0.0244 Hz, so f_min = 40 Hz; beta_min = 0.05 -> need = ln(100)/(2 pi 40 0.05) = 0.3665 s.
    A 4000-record history leaves 0.96 s (fine); 4090 records leave 0.06 s (warning)."""
    ui, mdir = _block(tmp_path, aopt="AOPT,0,0,0,0,0,0,0,0,0,0,1,0,0,0")
    ui.run_text("THFILE,H1.acc\nDAMP,0.05\nNOUT,1,1,0,0,0,0,0,1\nM,1,4e5,0.25,0.15,0.05,0.05\n")
    _history(mdir / "H1.acc", 4000)
    assert not any(m.number == "EDU-19" for m in run_check(ui.model, dirs=[mdir]).messages)
    _history(mdir / "H1.acc", 4090)
    assert any(m.number == "EDU-19" for m in run_check(ui.model, dirs=[mdir]).messages)


# ======================================================================================
# Option records: GUI API round trips, X-commands, new-model defaults (section 5.4, L17, UT-03)
# ======================================================================================
def test_every_dialog_dataclass_round_trips_through_command_text():
    import sassi.prep.commands  # noqa: F401  (registers the handlers)
    for name, spec in O.OPTION_SPECS.items():
        if spec.record is None:
            continue
        dc = O.options_class(name)()
        assert O.from_args(O.to_command(dc)) == dc, name
        ui = Interpreter()
        ui.execute("SITE,0,1,0,20,2" if name == "WAVE" else "N,1,0,0,0")
        line = O.to_command(dc)
        if spec.kind == "record":
            assert ui.execute(line), (name, line, ui.sink.texts(Kind.ERROR))
            assert O.options_from_record(O.get_record(ui.model, name)) == dc, name


def test_default_x_commands_are_not_written(tmp_path):
    """Section 3.4.R: X-commands appear in WRITE only when a value differs from its default."""
    ui = _ui("N,1,0,0,0\nMOTIONX,0,2,0,0,0,0,0,0,0,0\nSITEX,0\nANALYSX,0,0\nRELDX\nSTRESSX,1\n")
    ui.execute(f"WRITE,x.pre,{tmp_path}")
    body = [ln for ln in (tmp_path / "x.pre").read_text().splitlines() if not ln.startswith("*")]
    assert not any(ln.startswith(("MOTIONX", "SITEX", "ANALYSX", "RELDX")) for ln in body)
    assert "STRESSX,1" in body
    assert O.to_command(O.options_class("MOTIONX")(resp=1)) == "MOTIONX,0,1"


def test_new_model_dialog_defaults_section_5_4():
    """A few new-model defaults of the section 5.4 table (D-UI-03)."""
    rec = O.get_record(Interpreter().model, "MOTION")
    assert (rec.freq1, rec.freq2, rec.fstep, rec.cplx, rec.interp, rec.mult, rec.rec1) == (0.1, 100.0, 301, 1, 1, 1.0, 1)
    assert O.get_record(Interpreter().model, "ANALYS").prnt == 1
    assert O.get_record(Interpreter().model, "SOIL").ratio == 0.65
    assert O.get_record(Interpreter().model, "EQUAKE").dur == 20.0
    assert O.aopt_modules(Interpreter().model) == ["SITE", "POINT", "HOUSE", "ANALYS", "MOTION"]
    sx = O.get_record(Interpreter().model, "SOILX")
    assert (sx.mult, sx.get("max")) == (0.0, 0.1)


def test_from_args_reads_soilx_file_like_the_interpreter():
    """L17: the GUI parses and emits the same text as the console.  SOILX <file> is a rest-of-line
    text argument (the interpreter keeps 'a, b.acc' as one token); from_args must agree."""
    ui = Interpreter()
    ui.execute("SOILX,0,0,0.1,0,a, b.acc")
    assert O.get_record(ui.model, "SOILX").file == "a, b.acc"
    assert O.from_args("SOILX,0,0,0.1,0,a, b.acc").file == "a, b.acc"


def test_shared_command_changes_only_its_storage_location():
    ui = Interpreter()
    ui.execute("SITE,1,0,0.5,10,2,1,1,2,100,2,1,0.01,8192,3")
    before = O.get_record(ui.model, "SITE").to_tokens()
    for line in O.shared_command(ui.model, "nft", 2048):
        assert ui.execute(line)
    after = O.get_record(ui.model, "SITE").to_tokens()
    assert after[12] == "2048"
    assert after[:12] == before[:12] and after[13:] == before[13:]
    assert O.shared_value(ui.model, "nft") == 2048


# ======================================================================================
# RUN commands (section 2.6, D-RUN-03)
# ======================================================================================
def test_run_refuses_module_blocked_in_last_afwrite_even_if_a_deck_file_exists(tmp_path):
    ui, mdir = _block(tmp_path)
    assert ui.execute("AFWRITE")
    ui.execute("SITE,0,1,0,3,2,1,0,1,,1,0,0.01,4096,1")              # <nl> = 3: Error 47
    ui.execute("AFWRITE")
    assert (mdir / "blk.sit.bak").exists()
    shutil.copy(mdir / "blk.sit.bak", mdir / "blk.sit")              # user restores the stale deck
    ui.sink.clear()
    assert not ui.execute("RUNSITE")
    assert any("CHECK errors in the last AFWRITE" in e for e in ui.sink.texts(Kind.ERROR))


def test_simulation_commands_file(tmp_path):
    ui, mdir = _block(tmp_path, extra="POINT,0,1,0\n")                # POINT blocked (Error 57)
    ui.write_options["sim"] = True
    assert ui.execute("AFWRITE")
    lines = [ln for ln in (mdir / "blk-Sim.pre").read_text().splitlines() if not ln.startswith("*")]
    assert lines == ["RUNSITE", "RUNHOUSE", "RUNANALYS"]


def test_verify_survives_a_broken_problem_module(tmp_path, monkeypatch):
    """VERIFY must still run the registered problems when one verification module (written by
    another package) fails to import."""
    import sassi.verify as V
    import sassi.verify.problems.vp_checks  # noqa: F401  (registers VP-52)

    def broken():
        raise ImportError("vp_broken: simulated import failure")

    monkeypatch.setattr(V, "load_all", broken)
    ui = Interpreter(cwd=tmp_path)
    ui.execute("VERIFY,VP-52")
    assert any("VP-52 PASSED" in t for t in ui.sink.texts(Kind.INFO)), ui.sink.texts(Kind.ERROR)
