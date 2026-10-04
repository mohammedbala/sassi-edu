"""Property tables (spec 08 sections 5-6): UT-14 material conversions and section formulas, M/L/R/SC/MX
commands, deletion and listing commands."""
from __future__ import annotations

import math

import numpy as np
import pytest

from sassi.model import elastic_constants, section_circle, section_rectangle
from sassi.prep import Interpreter, Kind


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def run(ui, text):
    ui.run_text(text)
    return ui.model


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


# ------------------------------------------------------------------ UT-14
def test_ut14_material_types_agree(ui):
    E, nu, w, g = 30000.0, 0.25, 0.15, 32.2
    rho = w / g
    G = E / (2 * (1 + nu))
    M = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    vp, vs = math.sqrt(M / rho), math.sqrt(G / rho)
    m = run(ui, f"GRAVITY,32.2\nM,1,{E!r},{nu!r},{w!r},0.05,0.05,1\nM,2,{M!r},{G!r},{w!r},0.05,0.05,2\n"
                f"M,3,{vp!r},{vs!r},{w!r},0.05,0.05,3\n")
    ks = [m.materials[i].constants(m.gravity) for i in (1, 2, 3)]
    for k in ks:
        for name, ref in (("G", G), ("M", M), ("Vs", vs), ("Vp", vp), ("E", E), ("nu", nu)):
            assert getattr(k, name) == pytest.approx(ref, rel=1e-10), name
    assert ks[0].rho == pytest.approx(rho, rel=1e-15)
    assert ks[0].lam == pytest.approx(M - 2 * G, rel=1e-12)


def test_ut14_section_formulas():
    c = section_circle(1.0)
    assert c["flex2"] == pytest.approx(0.785398, abs=5e-7)
    assert c["flex3"] == pytest.approx(0.785398, abs=5e-7)
    assert c["tors"] == pytest.approx(1.570796, abs=5e-7)
    assert c["shear2"] == pytest.approx(0.9 * math.pi)
    r = section_rectangle(0.5, 1.0)
    assert r["tors"] == pytest.approx(0.028610, abs=5e-7)
    assert r["flex2"] == pytest.approx(1.0 * 0.5 ** 3 / 12) and r["flex3"] == pytest.approx(0.5 / 12)
    assert r["shear2"] == pytest.approx(0.5 / 1.2)
    # b > h: b and h swapped in the torsion formula only
    assert section_rectangle(1.0, 0.5)["tors"] == pytest.approx(r["tors"])
    # spec 08 section 5.30 example R,1,0.5,0.41667,0.41667,0.02861,0.010417,0.041667
    assert (r["axial"], round(r["shear2"], 5), round(r["tors"], 5), round(r["flex2"], 6), round(r["flex3"], 6)) == \
        (0.5, 0.41667, 0.02861, 0.010417, 0.041667)


def test_conversion_errors():
    with pytest.raises(ValueError):
        elastic_constants(4, 1, 1, 1, 32.2)
    with pytest.raises(ValueError):
        elastic_constants(1, 1, 0.5, 1, 32.2)
    with pytest.raises(ValueError):
        elastic_constants(1, 1, 0.2, 1, 0.0)


# ------------------------------------------------------------------ table commands
def test_table_commands(ui):
    m = run(ui, "M,1,519120,0.17,0.150,0.04,0.04,1\nM,2,1,2,3\nL,1,10.0,0.120,2000.0,1000.0,0.02,0.02\n"
                "R,1,0.5,0.41667,0.41667,0.02861,0.010417,0.041667\nSC,1,1.0E5,1.0E5,2.0E5,0,0,0,0.05\n")
    assert m.materials[1].val1 == 519120.0 and m.materials[1].mtype == 1
    assert m.materials[2].mtype == 1 and m.materials[2].pdamp == 0.0      # blank type -> 1, damping -> 0
    assert m.layers[1].vs == 1000.0 and m.layers[1].sdamp == 0.02
    assert m.sections[1].tors == 0.02861
    assert m.springs[1].k == (1e5, 1e5, 2e5, 0, 0, 0) and m.springs[1].damp == 0.05
    run(ui, "M,3,1,2,3,0.01,0.01,4\nM,0,1\nL,2,1,1,1,1,5,0\n")
    assert any("type> must be 1" in e for e in errors(ui))
    assert any("start at 1" in e for e in errors(ui))
    assert any("fractions" in w for w in ui.sink.texts(Kind.WARNING))


def test_matrix_rows_upper_triangle(ui):
    m = run(ui, "MXR,1,1,1000,0,0,0,0,0,-1000\nMXR,1,7,1000\nMXI,1,1,40,0,0,0,0,0,-40\nMXI,1,7,40\n"
                "MXM,1,1,2\nMXM,1,12,3\nMXR,1,12,1,2\nMXR,1,13,1\n")
    p = m.matrices[1]
    K = p.full("R")
    assert K[0, 0] == 1000 and K[0, 6] == -1000 and K[6, 0] == -1000 and K[6, 6] == 1000
    assert np.allclose(K, K.T)
    Ks = p.stiffness()
    assert Ks[0, 0] == 1000 + 40j
    assert p.full("M")[11, 11] == 3
    e = errors(ui)
    assert any("has 1 terms; 2 given" in t for t in e)
    assert any("<row> must be 1..12" in t for t in e)


def test_deletions(ui):
    m = run(ui, "M,1,1,0.2,1\nM,2,1,0.2,1\nM,3,1,0.2,1\nL,1,1,1,1,1\nR,1,1\nSC,1,1\nSC,2,1\nMXR,1,1,1\nMXR,4,1,1\n"
                "DELM,1,3,2\nDELL,1\nDELR,1\nDELSC,1,2\nMXDEL,1,4,3\n")
    assert sorted(m.materials) == [2]
    assert not m.layers and not m.sections and not m.springs and not m.matrices


def test_list_commands(ui):
    run(ui, "GRAVITY,32.2\nM,1,30000,0.25,0.15,0,0,1\nM,4,4000,1500,0.12,0.05,0.05,3\nL,1,10,0.12,2000,1000,0.02,0.02\n"
            "R,1,25\nSC,1,1\nMXR,1,1,5\nMLIST,4\nLLIST\nRLIST\nSCLIST\nMXLIST,1\nMXLIST,2\n")
    info = ui.sink.texts(Kind.INFO)
    rows = [t for t in info if t.strip().split(" ")[0] in ("1", "4")]
    assert any(t.lstrip().startswith("4 ") for t in rows) and not any(t.lstrip().startswith("1    1") for t in rows)
    assert any("derived (g = 32.2)" in t and "Vs 1500" in t for t in info)
    assert any("real stiffness" in t for t in info)
    assert any("matrix property 2 is not defined" in e for e in errors(ui))
