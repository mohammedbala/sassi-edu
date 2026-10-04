"""Final audit: one EXCSTRCHK rule for CHECK and HOUSE (G-15; manual 1.5.1 rules 4 and 11, 9.8.1,
application guideline 13a).  See docs/verification/AUDIT_REPORT.md, finding F-01.

Before the audit CHECK reported every excavation-interior node shared with the structure as an error,
while HOUSE accepted every such node that was an interaction node with a warning claiming that "the
flexible-volume general assembly is exact there".  The audit measured both configurations (scratch runs
documented in the report): in FV a flexible floor slab on the interior excavation nodes loses its resonance
completely (peak |H| 1.0 instead of 9.7), while a stiff or soil-like SOLID block on the excavation mesh
changes the transfer functions by 0.5-2.5 %; in FI-EVBN sharing changes them by 5-480 %.  Decision (one
helper, :func:`sassi.core.house_lib.excstrchk_kind`, used by both programs):

* error: the node is not an interaction node, or any structural element at it is not a SOLID/PLANE
  element built on the excavation mesh (beams, shells, springs, GENERAL elements, K nodes);
* warning: an interaction node used only by SOLID/PLANE elements whose nodes all belong to the
  excavation (near-field soil, a solid basement, the zero-SSI identity of VP-16), which the FV method
  tolerates ("often close results", rule 11) and reproduces exactly when the elements equal the soil.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sassi.core.house_lib import excstrchk_kind
from sassi.io.container import read_container
from sassi.prep import Interpreter
from sassi.prep.afwrite import DeckBuilder, write_deck
from sassi.prep.check import run_check
from sassi.verify import builders as B

# 2 x 2 x 2 block of excavated soil (27 nodes, interior node 14) in two 1 m layers over a half-space
BLOCK = """MDL,b,b
L,1,1.0,19.62,400,200,0.02,0.02
L,2,1.0,19.62,400,200,0.02,0.02
L,3,1.0,21.58,2000,1000,0.01,0.01
TOPL,1,2
FREQ,1,4,20,41,61
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.01,1024,1
WAVE,2,1,1,1,0
HOUSE,9.81,0,0,2,0,0,0,0,0
POINT,0,2,0.9
ANALYS,0,0,0,0,1,0,0,0,0,0,0
N,1,0,0,-2
N,3,2,0,-2
FILL,1,3
NGEN,2,3,1,3,1,0,1,0
NGEN,2,9,1,9,1,0,0,1
GROUP,1,SOLID
MACT,1
E,1,10,11,14,13,19,20,23,22
EGEN,1,1,1
EGEN,1,3,1,2
ETYPE,1,4,1,2
GROUP,2,SOLID
MACT,2
E,1,1,2,5,4,10,11,14,13
EGEN,1,1,1
EGEN,1,3,1,2
ETYPE,1,4,1,2
"""
CENTRE = 14
#: structure SOLIDs on the excavation mesh with the soil properties (K*_s = K*_e: zero-SSI identity)
SOIL_STRUCTURE = """M,1,400,200,19.62,0.02,0.02,3
GROUP,3,SOLID
MACT,1
E,1,1,2,5,4,10,11,14,13
EGEN,1,1,1
EGEN,1,3,1,2
EGEN,1,9,1,4
ETYPE,1,8,1,1
"""
#: a structural pile (BEAMS) from the interior node 14 to a node above grade
PILE = """M,2,3.0e7,0.2,24.0,0.05,0.05,1
R,1,0.25,0.2,0.2,0.01,0.005,0.005
N,30,1,1,3
N,31,5,1,3
GROUP,4,BEAMS
MACT,2
RACT,1
E,1,14,30,31
D,30,ALL
"""
#: a floor SHELL on the interior level z = -1 (nodes 10..18), structure
FLOOR = """M,2,3.0e7,0.2,24.0,0.05,0.05,1
GROUP,5,SHELL
MACT,2
E,1,10,11,14,13
EGEN,1,1,1
EGEN,1,3,1,2
THICK,1,4,1,0.1
"""


def _ui(tmp_path: Path, extra: str, fv: bool) -> Interpreter:
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(BLOCK + extra + ("INT,1,27,1,1\n" if fv else "INT,1,27,1,1\nINT,14,14,1,0\n"))
    return ui


def _check_kind(ui: Interpreter) -> str:
    rep = run_check(ui.model, modules=["HOUSE"])
    kinds = {m.kind for m in rep.messages if m.number == "EXCSTRCHK"}
    assert len(kinds) == 1, rep.format()
    return kinds.pop()


def _house(ui: Interpreter, wd: Path):
    """Run HOUSE on the deck AFWRITE would write for this model (written without the CHECK gate, so
    a model that CHECK rejects can be given to HOUSE as well)."""
    wd.mkdir(parents=True, exist_ok=True)
    db = DeckBuilder(ui.model)
    write_deck(wd / "b.sit", db.build("SITE"))
    write_deck(wd / "b.hou", db.build("HOUSE"))
    rc = B.run("HOUSE", wd, "b", check=False)
    return rc, B.listing(wd, "b", "HOUSE")


def test_rule_of_the_shared_helper():
    exc = set(range(1, 28))
    hexa = (1, 2, 5, 4, 10, 11, 14, 13)
    assert excstrchk_kind([(1, hexa)], exc, interaction=True) == "Warning"       # solid on the excavation mesh
    assert excstrchk_kind([(1, hexa)], exc, interaction=False) == "Error"        # FI: tied to -K_e
    assert excstrchk_kind([(2, (14, 30))], exc, interaction=True) == "Error"     # beam (flexible subsystem)
    assert excstrchk_kind([(3, (10, 11, 14, 13))], exc, interaction=True) == "Error"   # shell floor
    assert excstrchk_kind([(1, hexa), (2, (14, 30))], exc, interaction=True) == "Error"
    assert excstrchk_kind([(1, hexa[:7] + (99,))], exc, interaction=True) == "Error"   # solid off the mesh
    assert excstrchk_kind([(0, (14,))], exc, interaction=True) == "Error"        # K node
    assert excstrchk_kind([(4, (1, 2, 5, 4))], exc, interaction=True) == "Warning"     # PLANE (2D)


@pytest.mark.parametrize("name, extra, fv, kind", [
    ("solid-FV", SOIL_STRUCTURE, True, "Warning"),
    ("solid-FI", SOIL_STRUCTURE, False, "Error"),
    ("pile-FV", PILE, True, "Error"),
    ("pile-FI", PILE, False, "Error"),
    ("floor-FV", FLOOR, True, "Error"),
])
def test_check_and_house_classify_alike(tmp_path, name, extra, fv, kind):
    """The same model through CHECK and through HOUSE: the same severity for node 14."""
    ui = _ui(tmp_path, extra, fv)
    assert _check_kind(ui) == kind
    rc, out = _house(ui, tmp_path / name)
    if kind == "Error":
        assert rc == 1 and f"EXCSTRCHK: excavation interior node {CENTRE} is shared with the structure" in out, out
    else:
        assert rc == 0, out
        assert "EXCSTRCHK: 1 excavation interior nodes are shared with SOLID/PLANE structure" in out
        assert "exact there" not in out                     # the pre-audit claim is gone


def test_zero_ssi_identity_runs_through_the_command_path(tmp_path):
    """VP-16 (structure = excavated soil on the same mesh, FV) built with commands: CHECK warns but does
    not block, AFWRITE writes the HOUSE deck, and ANALYS reproduces the free field U = U'_f at every
    node (compared with FILE1 directly: vertical SV, k = 0, so U'_f(node) = U(interface))."""
    ui = _ui(tmp_path, SOIL_STRUCTURE + "AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0\n", fv=True)
    for cmd in ("CHECK", "AFWRITE", "RUNSITE", "RUNPOINT", "RUNHOUSE", "RUNANALYS"):
        ui.execute(cmd)
    md = tmp_path / "b"
    err = (md / "b.err").read_text()
    assert "Warning EXCSTRCHK" in err and "Error EXCSTRCHK" not in err
    f8 = read_container(md / "FILE8", "FILE8")
    f1 = read_container(md / "FILE1", "FILE1")
    f4 = read_container(md / "b.N4")
    U1 = np.asarray(f1["U"])[0]                              # (nF, nI, 3), single wave
    rows = {int(n): q for q, n in enumerate(np.asarray(f1["fnum"]))}
    H = np.asarray(f8["H"])
    eq_node, eq_dof = np.asarray(f8["eq_node"]), np.asarray(f8["eq_dof"])
    nodes, iface = np.asarray(f4["int_node"]), np.asarray(f4["int_iface"])
    assert len(nodes) == 27
    worst = 0.0
    for q, n in enumerate(np.asarray(f8["fnum"])):
        for nd, it in zip(nodes, iface):
            for d in (1, 2, 3):
                k = np.flatnonzero((eq_node == nd) & (eq_dof == d))
                worst = max(worst, abs(H[q, k[0]] - U1[rows[int(n)], it - 1, d - 1]))
    assert worst < 1e-8
