"""Tutorial examples run through the command language, and VP-E1 (AFWRITE writes the intended decks).

The tutorial models of ``examples/`` are ordinary ``.pre`` files: a user runs them with
``sassi run <file>.pre`` (from the ``examples`` directory), the interpreter builds the model, CHECK and
AFWRITE write the module decks and the RUN<MODULE> commands run the modules (ARCHITECTURE sections
3 and 8, requirements 2.4 and 3.3).  This module

* stages an example in a work directory and runs it through :class:`sassi.prep.Interpreter`
  (:func:`stage_example`, :func:`run_example`; also used by ``tests/integration/test_examples.py``,
  which runs the examples listed in :data:`EXAMPLES`);
* defines **VP-E1**: example 1 (lumped-mass stick on a rigid surface mat, vertical SV input) run
  through the interpreter gives the same FILE8 (relative difference <= 1e-12) as the same model
  written *by hand* as module decks with the primitives of :mod:`sassi.verify.builders`
  (:class:`~sassi.verify.builders.HouseBuilder` with the node numbers of the ``.pre`` file,
  :func:`~sassi.verify.builders.site_deck`, :func:`~sassi.verify.builders.point_deck`,
  :func:`~sassi.verify.builders.analys_deck`).  The deck contents themselves are compared first, so a
  failure says *which* AFWRITE resolution differs (global coordinates, fixities, ETYPE, the L/M/R
  tables, masses with their units flag, interaction nodes, the TOPL layers with the half-space row,
  the frequency numbers and ``df``; requirements 3.3, D-AFW-01..06).

Why it matters: every example and every user model reaches the modules through AFWRITE.  The module
verification problems (VP-01 ... VP-41) write their decks directly; VP-E1 closes the loop between the
command language and those decks.
"""
from __future__ import annotations

import math

import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ...io import decks
from .. import VPResult, problem
from .. import builders as B

#: repository directory of the tutorial models (``<repo>/examples``)
EXAMPLES_DIR = Path(__file__).resolve().parents[3] / "examples"
#: DYNP curve library that example 4 reads with INP,../sassi/data/dynp_library.pre
DYNP_LIBRARY = Path(__file__).resolve().parents[2] / "data" / "dynp_library.pre"

#: example file stems and their model directories (MDL paths relative to ``examples``)
EXAMPLES: Dict[str, Tuple[str, ...]] = {
    "ex01_surface_stick": ("ex01",),
    "ex02_embedded_box": ("ex02", "ex02_fsin", "ex02_evbn"),
    "ex03_forced_vibration": ("ex03",),
    "ex04_site_response": ("ex04",),
    "ex05_xyz_simultaneous": ("ex05",),
    "ex08_embedded_building": ("ex08", "ex08_fsin", "ex08_evbn"),
    "ex09_steel_frame": ("ex09", "ex09_fixed"),
}

#: binary intermediate files deleted by :func:`purge_binaries` (decks, listings and text results kept)
BINARY_PREFIXES = ("FILE1", "FILE2", "FILE3", "FILE8", "FILE9", "FILE11", "COO", "DOFSMAP")


# ======================================================================================
# Staging and running an example through the interpreter
# ======================================================================================
def stage_example(name: str, root: Path) -> Path:
    """Copy ``examples/<name>.pre`` and ``examples/data`` to ``<root>/examples`` (and the DYNP library
    to ``<root>/sassi/data``, the relative location example 4 reads it from).  Returns the staged
    ``.pre`` path; the example's model directories are created next to it by MDL."""
    src = EXAMPLES_DIR / f"{name}.pre"
    if not src.exists():
        raise FileNotFoundError(f"tutorial example {src} not found (the examples ship with the source tree)")
    ex = Path(root) / "examples"
    ex.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, ex / src.name)
    if (EXAMPLES_DIR / "data").is_dir() and not (ex / "data").exists():
        shutil.copytree(EXAMPLES_DIR / "data", ex / "data")
    lib = Path(root) / "sassi" / "data"
    if DYNP_LIBRARY.exists() and not (lib / DYNP_LIBRARY.name).exists():
        lib.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(DYNP_LIBRARY, lib / DYNP_LIBRARY.name)
    return ex / src.name


def run_example(name: str, root: Path, stop_after: Optional[str] = None):
    """Stage example ``name`` under ``root`` and run it as ``sassi run`` does from the examples
    directory: an :class:`~sassi.prep.Interpreter` whose working directory is ``examples`` executes the
    file (INP semantics).  ``stop_after`` (a command line such as ``'RUNANALYS'``) runs the file only up
    to and including the first line equal to it (the staged file is shortened accordingly).
    Returns ``(interpreter, summary, staged .pre path)``."""
    from ...prep import Interpreter
    pre = stage_example(name, root)
    if stop_after:
        lines = pre.read_text(encoding="utf-8").splitlines()
        stops = [i for i, ln in enumerate(lines) if ln.strip().upper() == stop_after.strip().upper()]
        if not stops:
            raise ValueError(f"{pre.name} has no line {stop_after!r}")
        pre.write_text("\n".join(lines[:stops[0] + 1]) + "\n", encoding="utf-8")
    ui = Interpreter(cwd=pre.parent)
    summary = ui.run_file(str(pre), resolve=False)
    return ui, summary, pre


def purge_binaries(directory: Path) -> int:
    """Delete the binary inter-module files of a finished run (they are large and the results have been
    read); decks, listings and text results stay -- FILE88 of SOIL, for example, is a text file and is
    kept.  Returns the number of bytes freed."""
    freed = 0
    d = Path(directory)
    if not d.is_dir():
        return 0
    for p in d.iterdir():
        if p.is_file() and (p.name.startswith(BINARY_PREFIXES) or p.suffix == ".N4") and zipfile.is_zipfile(p):
            freed += p.stat().st_size
            p.unlink()
    return freed


def listing_status(directory: Path) -> Dict[str, str]:
    """``{listing name: final status line}`` of the module listings ``*_<MODULE>.out`` of a directory."""
    out = {}
    for p in sorted(Path(directory).glob("*_*.out")):
        lines = [ln for ln in p.read_text(encoding="utf-8", errors="replace").splitlines() if "finished with status" in ln]
        out[p.name] = lines[-1].strip() if lines else "no status line"
    return out


# ======================================================================================
# Example 1 written by hand (VP-E1)
# ======================================================================================
GRAVITY = 9.81
#: ex01 frequency numbers (FREQ set 1), dt 0.005 s, NFFT 8192
EX01_FNUM = [4, 20, 41, 61, 82, 102, 123, 143, 164, 184, 205, 225, 246, 287, 328, 369, 410, 492, 573, 655, 737, 819]
#: (thick, weight, Vp, Vs, damping) of L 1, L 2 and the half-space L 3
EX01_L = {1: (0.5, 19.0, 600.0, 300.0, 0.05), 2: (1.0, 20.0, 1000.0, 500.0, 0.04), 3: (1.0, 21.0, 2000.0, 1000.0, 0.02)}
EX01_TOPL = [1] * 10 + [2] * 12
#: mat nodes around the stick base 41 tied to it by the rigid spider
EX01_SPIDER = (31, 32, 33, 40, 42, 49, 50, 51)


def ex01_site() -> B.Site:
    """The TOPL layers and the half-space of example 1 (weights converted with g = 9.81)."""
    def layer(no: int, thick: Optional[float] = None) -> B.SoilLayer:
        t, w, vp, vs, beta = EX01_L[no]
        return B.SoilLayer(t if thick is None else thick, vs, vp, w / GRAVITY, beta, beta)
    return B.Site([layer(n) for n in EX01_TOPL], layer(3, 0.0), gravity=GRAVITY, nl=20)


def ex01_house(site: B.Site) -> B.HouseBuilder:
    """Example 1 as HOUSE deck data, numbered exactly as the ``.pre`` file numbers it.

    * nodes 1-81: 9 x 9 mat grid (2.5 m) at grade, ROTZ (shell drilling) fixed by FIXROT at the
      shell-only nodes, i.e. all but the stick base 41 and the spider ends; 82-85: the stick at
      z = 5, 10, 15, 20; 86-89: beam K nodes at (10, 0, 2.5 + 5 k); 90: K node of the rigid links at
      (0, 0, 2.5).  K-only nodes are written with every DOF fixed (D-CHK-08), as AFWRITE does;
    * group 1 SHELL basemat (M 1, 1.5 m), group 2 BEAMS stick (M 2, R 1), group 3 BEAMS rigid spider
      (M 3, R 2) from the stick base 41 to its 8 neighbours;
    * masses: 9810 kN at each floor as *weights* (MT with the MUNITS default 1).
    """
    hb = B.HouseBuilder(site, dim=2, imp=0, incomp=1, title="Ex01 - 4-mass stick on a rigid surface mat, "
                                                                "vertical SV (X) input")
    beam_ends = {41} | set(EX01_SPIDER)            # mat nodes with beams: drilling rotation stays free
    for j in range(9):
        for i in range(9):
            n = 9 * j + i + 1
            hb.add_node(-10.0 + 2.5 * i, -10.0 + 2.5 * j, 0.0, fix=(0, 0, 0, 0, 0, int(n not in beam_ends)), nid=n)
    for k in range(4):
        hb.add_node(0.0, 0.0, 5.0 * (k + 1), nid=82 + k)
    for k in range(4):
        hb.add_node(10.0, 0.0, 2.5 + 5.0 * k, fix=(1,) * 6, nid=86 + k, register=False)
    hb.add_node(0.0, 0.0, 2.5, fix=(1,) * 6, nid=90, register=False)
    m_mat = hb.material(1, 3.0e8, 0.2, 24.0, 0.05, 0.05)
    m_stick = hb.material(1, 3.0e7, 0.2, 0.0, 0.05, 0.05)
    m_rigid = hb.material(1, 3.0e10, 0.2, 0.0, 0.0, 0.0)
    r_core = hb.section(20.0, 10.0, 10.0, 400.0, 200.0, 200.0)
    r_link = hb.section(6.25, 5.2, 5.2, 5.5, 3.26, 3.26)
    g_mat, g_stick, g_link = hb.group(3, "basemat"), hb.group(2, "stick"), hb.group(2, "rigid links")
    for j in range(8):
        for i in range(8):
            n1 = 9 * j + i + 1
            hb._element(g_mat, (n1, n1 + 1, n1 + 10, n1 + 9), etype=1, mat=m_mat, prop=1, thick=1.5)
    for (i, j, k) in ((41, 82, 86), (82, 83, 87), (83, 84, 88), (84, 85, 89)):
        hb._element(g_stick, (i, j, k), etype=1, mat=m_stick, prop=r_core)
    for n in EX01_SPIDER:
        hb._element(g_link, (41, n, 90), etype=1, mat=m_rigid, prop=r_link)
    for n in range(82, 86):
        hb.mass(n, 9810.0, 9810.0, 9810.0, units=1)
    hb.set_interaction(range(1, 82))
    return hb


def write_ex01_by_hand(wd: Path, model: str = "ex01") -> Dict[str, decks.Deck]:
    """Write the SITE, POINT, HOUSE and ANALYS decks of example 1 by hand into ``wd``."""
    site = ex01_site()
    fs = B.FrequencySet.fourier(0.005, 8192, EX01_FNUM)
    title = "Ex01 - 4-mass stick on a rigid surface mat, vertical SV (X) input"
    out = {"SITE": B.site_deck(site, fs, wave="SV", cl=1, model=model, title=title),
           "POINT": B.point_deck(fs, layer=0, rad=2.25, model=model, title=title),
           "HOUSE": fs.fill(ex01_house(site).deck(model, title)),
           "ANALYS": B.analys_deck(fs, gravity=GRAVITY, model=model, title=title, type=0, mode=0, save=0,
                                   prnt=1, fopt=0, impe=0, simul=0)}
    for d in out.values():
        B.write_deck(wd, model, d)
    return out


# ======================================================================================
# Deck comparison helpers
# ======================================================================================
def _rows(d: decks.Deck, table: str, cols: Sequence[str]) -> List[list]:
    return [[r[c] for c in cols] for r in d.rows(table)]


def compare_rows(r: VPResult, label: str, a: List[list], b: List[list], rtol: float = 1e-12) -> None:
    """Row-by-row equality of two tables (numbers to ``rtol``, strings exactly)."""
    if not r.require(f"{label}: same number of rows ({len(a)} vs {len(b)})", len(a) == len(b)):
        return
    worst, where = 0.0, ""
    same_text = True
    for i, (ra, rb) in enumerate(zip(a, b)):
        for j, (x, y) in enumerate(zip(ra, rb)):
            if isinstance(x, str) or isinstance(y, str):
                if str(x) != str(y):
                    same_text = False
                    where = where or f"row {i + 1} column {j + 1}: {x!r} vs {y!r}"
                continue
            e = abs(float(x) - float(y)) / max(abs(float(y)), 1.0)
            if math.isnan(worst):
                continue                       # a NaN stays the result (final audit: it used to be skipped)
            if math.isnan(e) or e > worst:
                worst, where = e, f"row {i + 1} column {j + 1}: {x!r} vs {y!r}"
    r.require(f"{label}: text columns identical", same_text, where if not same_text else "")
    r.check(f"{label}: max relative difference", worst, 0.0, atol=rtol, note=where)


def compare_decks(r: VPResult, afw: Dict[str, decks.Deck], hand: Dict[str, decks.Deck]) -> None:
    """The quantities each module uses must be identical in the AFWRITE and the hand-written decks."""
    a, h = afw["SITE"], hand["SITE"]
    lay = ("thick", "weight", "vp", "vs", "dp", "ds")
    compare_rows(r, "SITE layers (TOPL resolved from L)", _rows(a, "layers", lay), _rows(h, "layers", lay))
    compare_rows(r, "SITE half-space row", _rows(a, "halfspace", lay[1:]), _rows(h, "halfspace", lay[1:]))
    compare_rows(r, "SITE wave records", _rows(a, "waves", ("type", "opt", "ratio1", "ratio2", "angle")),
                 _rows(h, "waves", ("type", "opt", "ratio1", "ratio2", "angle")))
    for mod in ("SITE", "POINT", "HOUSE", "ANALYS"):
        compare_rows(r, f"{mod} frequency numbers", _rows(afw[mod], "freqs", ("number",)),
                     _rows(hand[mod], "freqs", ("number",)))
    scal = {"SITE": ("mode1", "mode2", "nl", "cl", "cm", "wopt", "delt", "nft", "df", "gravity", "cmodform"),
            "POINT": ("layer", "rad", "dim", "df"),
            "HOUSE": ("gravity", "gelev", "dim", "coh", "wpass", "me", "cmodform", "delt", "nft", "df"),
            "ANALYS": ("type", "mode", "save", "fopt", "ang", "xc", "yc", "zc", "impe", "simul", "delt", "nft", "df",
                       "gravity")}
    for mod, names in scal.items():
        compare_rows(r, f"{mod} parameters {', '.join(names)}", [[afw[mod][n] for n in names]],
                     [[hand[mod][n] for n in names]])
    a, h = afw["HOUSE"], hand["HOUSE"]
    compare_rows(r, "HOUSE nodes (global coordinates, fixities)",
                 _rows(a, "nodes", ("id", "x", "y", "z", "fx", "fy", "fz", "frx", "fry", "frz")),
                 _rows(h, "nodes", ("id", "x", "y", "z", "fx", "fy", "fz", "frx", "fry", "frz")))
    ecols = ("group", "id", "etype", "mat", "prop", "eint", "thick", "n1", "n2", "n3", "n4", "n5", "n6", "n7",
             "n8", "ki", "kj")
    compare_rows(r, "HOUSE elements", _rows(a, "elements", ecols), _rows(h, "elements", ecols))
    compare_rows(r, "HOUSE groups", _rows(a, "groups", ("id", "type")), _rows(h, "groups", ("id", "type")))
    mcols = ("id", "type", "val1", "val2", "weight", "pdamp", "sdamp")
    compare_rows(r, "HOUSE materials (M table)", _rows(a, "materials", mcols), _rows(h, "materials", mcols))
    bcols = ("id", "axial", "shear2", "shear3", "tors", "flex2", "flex3")
    compare_rows(r, "HOUSE beam sections (R table)", _rows(a, "beamprops", bcols), _rows(h, "beamprops", bcols))
    scols = ("node", "mx", "my", "mz", "mxx", "myy", "mzz", "units")
    compare_rows(r, "HOUSE masses (MT, weight units)", _rows(a, "masses", scols), _rows(h, "masses", scols))
    compare_rows(r, "HOUSE interaction nodes", _rows(a, "interaction", ("id",)), _rows(h, "interaction", ("id",)))
    compare_rows(r, "HOUSE site layers (TOPL + half-space row)", _rows(a, "sitelayers", lay),
                 _rows(h, "sitelayers", lay))


# ======================================================================================
# VP-E1
# ======================================================================================
@problem("VP-E1", "Example 1 through the interpreter = the same decks written by hand (FILE8 1e-12)", tier="P0",
         modules=["AFWRITE", "SITE", "POINT", "HOUSE", "ANALYS"],
         source="examples/ex01_surface_stick.pre; requirements 3.3, D-AFW-01..06; ARCHITECTURE 3")
def vp_e1(workdir: Path) -> VPResult:
    r = VPResult()
    workdir = Path(workdir)
    src = EXAMPLES_DIR / "ex01_surface_stick.pre"
    if not r.require("tutorial example examples/ex01_surface_stick.pre available", src.exists(),
                     f"{src} not found: the tutorial examples ship with the source tree, not with an installed "
                     "package; run VP-E1 from a source checkout"):
        r.notes.append("VP-E1 not run: the examples directory of the source tree is not available")
        return r
    # the interpreter run: example 1 up to RUNANALYS (MOTION, STRESS and RELDISP only read FILE8)
    ui, summary, pre = run_example("ex01_surface_stick", workdir / "interpreter", stop_after="RUNANALYS")
    mdir = pre.parent / "ex01"
    r.require("ex01 runs through the interpreter without errors", summary.ok and summary.errors == 0,
              summary.text())
    status = listing_status(mdir)
    for mod in ("SITE", "POINT", "HOUSE", "ANALYS"):
        line = status.get(f"ex01_{mod}.out", "listing missing")
        r.require(f"RUN{mod} finished with status OK", "status OK" in line, line)
    if not r.require("the interpreter run wrote FILE8", (mdir / "FILE8").exists(), str(mdir / "FILE8")):
        r.notes.append("the interpreter run produced no FILE8; nothing to compare")
        return r
    afw = {mod: decks.read(decks.deck_path(mdir, "ex01", mod), mod) for mod in ("SITE", "POINT", "HOUSE", "ANALYS")}
    fa = B.read_file8(mdir)
    purge_binaries(mdir)                          # keep the disk footprint small (results are in memory)
    # the same model as decks written by hand, run module by module (what RUN<MODULE> does)
    hand_dir = workdir / "by_hand"
    hand = write_ex01_by_hand(hand_dir)
    compare_decks(r, afw, hand)
    for mod in ("SITE", "POINT"):
        B.run(mod, hand_dir, "ex01")
    (hand_dir / "FILE2").unlink()                 # FILE2 is read by POINT only
    for mod in ("HOUSE", "ANALYS"):
        B.run(mod, hand_dir, "ex01")
    fh = B.read_file8(hand_dir)
    purge_binaries(hand_dir)
    same_map = (np.array_equal(fa["eq_node"], fh["eq_node"]) and np.array_equal(fa["eq_dof"], fh["eq_dof"])
                and np.array_equal(fa["fnum"], fh["fnum"]))
    r.require("FILE8 equation numbering and frequency numbers identical", same_map)
    r.check("FILE8 frequency step df (Hz)", float(fa.meta["df"]), float(fh.meta["df"]), rtol=1e-14)
    if same_map:
        Ha, Hh = np.asarray(fa["H"]), np.asarray(fh["H"])
        err = float(np.max(np.abs(Ha - Hh)) / np.max(np.abs(Hh)))
        r.check("FILE8 max |H_interpreter - H_by_hand| / max |H|", err, 0.0, atol=1e-12)
        top = B.tf(fa, 85, 1)
        r.notes.append(f"FILE8: {Ha.shape[1]} equations x {Ha.shape[0]} frequencies; relative difference {err:.2e}; "
                       f"top-floor |ATF_x| peak {np.max(np.abs(top)):.3f} at "
                       f"{float(fa['freq'][int(np.argmax(np.abs(top)))]):.3f} Hz")
    return r
