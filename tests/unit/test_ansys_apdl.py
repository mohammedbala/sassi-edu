"""ANSYS command: APDL export (requirements 3.4.K, D-ANS-06, R1 section 1.2; spec 04 section 12.3).

The exported text is checked line by line, and read back with the CONVERT,ANSYS reader (which
reads the plain APDL commands the export writes): model -> APDL -> model must keep every
exported field, and the assembled structural stiffness and mass must be unchanged.
"""
from __future__ import annotations

import numpy as np
import pytest

from sassi.conventions import cfactor
from sassi.io.ansys_cdb import convert_cdb, read_cdb_text
from sassi.io.apdl import (ApdlOptions, _release_keyopt, ansys_damping, cmblock_lines, export_apdl)
from sassi.prep import Interpreter, Kind
from sassi.verify.problems.vp_ansys import assemble_model

MODEL = """
TIT,APDL export test
GRAVITY,32.2
GROUNDELEV,0
MOPT,0
N,1,0,0,0
N,2,1,0,0
N,3,1,1,0
N,4,0,1,0
N,5,0,0,1
N,6,1,0,1
N,7,1,1,1
N,8,0,1,1
N,21,0,0,-1
N,22,1,0,-1
N,23,1,1,-1
N,24,0,1,-1
N,9,2,0,1
N,10,2,1,1
N,11,1,0,4
N,12,3,0,1
N,13,3,0,4
N,15,2,0,0
N,16,2,2,1
N,17,1,2,1
M,1,4.32e5,0.25,0.15,0.05,0.05,1
M,2,3e6,0.2,0.16,0.03,0.07,1
L,1,10,0.12,1500,800,0.05,0.05
R,1,1,0.8,0.7,0.2,0.1,0.3
SC,1,100,0,200,0,0,5,0.02
MXR,1,1,1000,0,0,0,0,0,-1000
MXR,1,7,1000
MXR,1,2,500,0,0,0,0,0,-500
MXR,1,8,500
MXM,1,1,2
MXM,1,7,2
GROUP,1,SOLID
E,1,1,2,3,4,5,6,7,8
E,2,21,22,23,24,1,2,3,4
GROUP,2,SHELL
E,1,5,6,7,8
E,2,6,9,10,7
E,3,6,9,10
THICK,1,3,1,0.3
MSET,1,3,1,2
GROUP,3,BEAMS
E,1,6,11,12
E,2,11,13,12
KI,2,2,1,0,0,0,0,1,1
GROUP,4,SPRING
E,1,15,9
GROUP,5,GENERAL
E,1,10,16,17
MT,11,2,2,2
MR,11,0.1,0.2,0.3
MUNITS,11,11,1,0
MT,13,64.4,64.4,64.4
D,21,24,1,1,ALL
D,15,15,1,1,ALL
D,1,4,1,1,UZ
INT,1,4,1,1
INT,21,24,1,1
"""


def _model(text=MODEL):
    ui = Interpreter()
    ui.run_text(text)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    return ui.model


def _lines(text):
    return [ln.split("!")[0].strip() for ln in text.splitlines() if ln.split("!")[0].strip()]


# ======================================================================================
# damping mapping (R1 1.2) and small helpers
# ======================================================================================
@pytest.mark.parametrize("beta", [0.0, 0.02, 0.05, 0.1, 0.2])
def test_exact_damping_mapping_reproduces_the_sassi_complex_modulus(beta):
    E = 3.0e7
    Ea, g = ansys_damping(E, beta, "exact")
    assert Ea * (1 + 1j * g) == pytest.approx(E * cfactor(beta), rel=1e-14)
    Er, gr = ansys_damping(E, beta, "ratio")
    assert (Er, gr) == (E, 2 * beta)                                      # D-ANS-06: DMPR = 2 beta
    assert ansys_damping(E, beta, "none") == (E, 0.0)


def test_release_keyopt_digits():
    # SASSI KI P1 P2 P3 M1 M2 M3 -> BEAM44 digits UX UY UZ ROTX ROTY ROTZ (UY = P3, UZ = P2, ROTY = M3, ROTZ = M2)
    assert _release_keyopt([0, 0, 0, 0, 1, 1]) == 11
    assert _release_keyopt([0, 0, 0, 0, 0, 1]) == 10                       # M3 = about axis 3 = -y -> ROTY
    assert _release_keyopt([0, 0, 0, 0, 1, 0]) == 1                        # M2 -> ROTZ
    assert _release_keyopt([1, 1, 0, 1, 0, 0]) == 101100


def test_cmblock_ranges():
    assert cmblock_lines("SSI_INT", "NODE", [4, 1, 2, 3, 10, 21, 22]) == [
        "CMBLOCK,SSI_INT,NODE,       5", "(8i10)", "         1        -4        10        21       -22"]


# ======================================================================================
# legacy export
# ======================================================================================
def test_legacy_export_text():
    m = _model()
    res = export_apdl(m)
    L = _lines(res.text)
    assert L[0] == "/PREP7" and L[1] == "/TITLE,APDL export test" and L[-1] == "FINISH"
    # element types (one per group; beams by release pattern; springs per component; MATRIX27 K and M; MASS21)
    ets = [ln for ln in L if ln.startswith("ET,")]
    assert ets == ["ET,1,SOLID45", "ET,2,SHELL63", "ET,3,BEAM4", "ET,4,BEAM44", "ET,5,COMBIN14", "ET,6,COMBIN14",
                   "ET,7,COMBIN14", "ET,8,MATRIX27", "ET,9,MATRIX27", "ET,10,MASS21"]
    for k in ("KEYOPT,1,1,0", "KEYOPT,4,7,11", "KEYOPT,4,8,0", "KEYOPT,5,2,1", "KEYOPT,6,2,3", "KEYOPT,7,2,6",
              "KEYOPT,8,3,4", "KEYOPT,9,3,2", "KEYOPT,10,3,0"):
        assert k in L, k
    # materials: DENS = weight/g, DMPR = 2 beta_s
    assert "MP,EX,1,432000" in L and "MP,PRXY,1,0.25" in L and f"MP,DENS,1,{0.15 / 32.2!r}" in L
    assert "MP,DMPR,1,0.1" in L and "MP,DMPR,2,0.14" in L
    assert any("beta_p = 0.03 != beta_s = 0.07" in w for w in res.warnings)
    # real constants: SHELL63 thickness, BEAM4 (AREA IZZ=I2 IYY=I3 ... IXX=J SHEARZ=A/As2 SHEARY=A/As3)
    assert "R,1,0.3" in L
    i = L.index(next(ln for ln in L if ln.startswith("R,2,")))
    assert L[i].split(",")[2:5] == ["1", "0.1", "0.3"]
    assert [float(v) for v in L[i + 1].split(",")[1:]] == pytest.approx([0, 0.2, 1 / 0.8, 1 / 0.7])   # ISTRN IXX SHEAR
    # nodes in global coordinates, elements with EN
    assert "N,13,3,0,4" in L and "EN,1,1,2,3,4,5,6,7,8" in L
    assert "EN,4,6,9,10,10" in L                                            # shell triangle I J K K
    # excavated element 2 of group 1 not exported
    assert not any(ln.startswith("EN,") and ln.endswith("21,22,23,24,1,2,3,4") for ln in L)
    assert any("excavated-soil" in n for n in res.notes)
    # constraints only on DOFs the node has in ANSYS
    assert "D,1,UZ,0" in L and "D,15,ALL,0" in L and not any(ln.startswith("D,21,") for ln in L)
    # components: interaction nodes and groups
    k = L.index("CMBLOCK,SSI_INT,NODE,       4")
    assert L[k + 2].split() == ["1", "-4", "21", "-24"]
    assert any(ln.startswith("CMBLOCK,SSI_G3,ELEM") for ln in L)
    # spring damping and MXI warnings
    assert any("SC spring damping" in w for w in res.warnings)
    assert (1, 1, 1) in res.element_map


def test_export_is_deterministic():
    m = _model()
    assert export_apdl(m).text == export_apdl(m).text


def test_exact_damping_option_text():
    m = _model()
    L = _lines(export_apdl(m, ApdlOptions(damping="exact")).text)
    b = 0.05
    E = float(next(ln for ln in L if ln.startswith("MP,EX,1,")).split(",")[3])
    g = float(next(ln for ln in L if ln.startswith("MP,DMPR,1,")).split(",")[3])
    assert E == pytest.approx(4.32e5 * (1 - 2 * b * b), rel=1e-15)
    assert g == pytest.approx(2 * b * np.sqrt(1 - b * b) / (1 - 2 * b * b), rel=1e-15)
    L0 = _lines(export_apdl(m, ApdlOptions(damping="none")).text)
    assert not any(ln.startswith("MP,DMPR") for ln in L0)
    with pytest.raises(ValueError):
        export_apdl(m, ApdlOptions(damping="viscous"))


def _same_matrices(m1, m2):
    d1, a1 = assemble_model(m1)
    d2, a2 = assemble_model(m2)
    k1 = list(zip(d1.eq_node.tolist(), d1.eq_dof.tolist()))
    k2 = list(zip(d2.eq_node.tolist(), d2.eq_dof.tolist()))
    assert sorted(k1) == sorted(k2)
    p = [k2.index(k) for k in k1]
    for A1, A2 in ((a1.Ks.toarray().real, a2.Ks.toarray().real), (a1.Ms.toarray(), a2.Ms.toarray())):
        B = A2[np.ix_(p, p)]
        np.testing.assert_allclose(B, A1, rtol=1e-11, atol=1e-11 * np.abs(A1).max())


def test_legacy_round_trip_keeps_the_structure():
    """model -> APDL -> CONVERT,ANSYS reader -> model: fields and assembled K, M unchanged."""
    m = _model()
    res = export_apdl(m)
    back = convert_cdb(read_cdb_text(res.text), 32.2).model
    for n in m.nodes:
        np.testing.assert_allclose(back.node_global(n), m.node_global(n), atol=1e-14)
    assert back.nodes[1].fix == [0, 0, 1, 0, 0, 0] and back.nodes[15].fix == [1] * 6
    assert sorted(n for n, nd in back.nodes.items() if 0 in nd.flags) == [1, 2, 3, 4, 21, 22, 23, 24]
    mt1, mt2 = back.materials[1], back.materials[2]
    assert (mt1.val1, mt1.val2, mt1.sdamp, mt1.pdamp) == pytest.approx((4.32e5, 0.25, 0.05, 0.05))
    assert mt1.weight == pytest.approx(0.15, rel=1e-14)
    assert (mt2.sdamp, mt2.pdamp) == pytest.approx((0.07, 0.07))       # one ANSYS damping value
    # beams: same section; the released element keeps its releases (BEAM44)
    b = {e.id: e for e in back.groups[4].elements.values()}
    assert b[1].ki == [0, 0, 0, 0, 1, 1] and b[1].nodes == [11, 13, 12]
    s = back.sections[b[1].prop]
    assert (s.axial, s.shear2, s.shear3, s.tors, s.flex2, s.flex3) == pytest.approx((1, 0.8, 0.7, 0.2, 0.1, 0.3),
                                                                                     rel=1e-14)
    # shells: thickness and material
    sh = back.groups[2].sorted_elements()
    assert [e.thick for e in sh] == [0.3] * 3 and [e.mat for e in sh] == [2] * 3 and sh[2].nodes == [6, 9, 10]
    # masses in mass units
    assert back.tmass[11] == [2, 2, 2] and back.rmass[11] == [0.1, 0.2, 0.3]
    assert back.tmass[13] == pytest.approx([2, 2, 2]) and back.mass_units[13] == 0
    # GENERAL: the local 3-node element comes back in global axes (e1 = +Y, e2 = -X)
    gk = [g for g in back.groups.values() if g.type == 9]
    K = sum(back.matrices[g.elements[1].prop].full("R") for g in gk)
    assert K[1, 1] == pytest.approx(1000) and K[1, 7] == pytest.approx(-1000) and K[0, 0] == pytest.approx(500)
    # the structure: identical stiffness and mass (excavated soil excluded on both sides)
    _same_matrices(m, back)


def test_modern_export_and_round_trip():
    m = _model()
    res = export_apdl(m, ApdlOptions(modern=True))
    L = _lines(res.text)
    ets = [ln for ln in L if ln.startswith("ET,")]
    assert ets[:3] == ["ET,1,SOLID185", "ET,2,SHELL181", "ET,3,BEAM188"]
    assert "KEYOPT,1,2,3" in L                                       # incompatible modes (MOPT 0)
    assert "SECTYPE,1,SHELL" in L and "SECDATA,0.3,2,0,3" in L
    assert any(ln.startswith("SECTYPE,") and ln.endswith("BEAM,ASEC") for ln in L)
    assert any("releases are not written for BEAM188" in w for w in res.warnings)
    back = convert_cdb(read_cdb_text(res.text), 32.2).model
    s = [back.sections[e.prop] for _, e in back.iter_elements() if e.prop in back.sections][0]
    assert (s.axial, s.shear2, s.shear3, s.tors, s.flex2, s.flex3) == pytest.approx((1, 0.8, 0.7, 0.2, 0.1, 0.3),
                                                                                     rel=1e-12)
    m.groups[3].elements[2].ki = [0] * 6                             # modern export drops the release
    _same_matrices(m, back)


def test_two_d_export_rotates_back_to_xy():
    text = """GRAVITY,9.81
N,1,0,0,0
N,2,2,0,0
N,3,2,0,1
N,4,0,0,1
N,5,4,0,0
M,1,2e5,0.3,19.62,0.05,0.05,1
GROUP,1,PLANE
E,1,1,2,3,4
E,2,2,5,3
D,1,1,1,1,UX,UZ
D,5,5,1,1,UZ
"""
    m = _model(text)
    res = export_apdl(m)
    L = _lines(res.text)
    assert "ET,1,PLANE42" in L and "KEYOPT,1,3,2" in L
    assert "N,3,2,1,0" in L                                          # (x, y, z)_SASSI -> (x, z, -y)
    assert "D,1,UX,0" in L and "D,1,UY,0" in L and "D,5,UY,0" in L   # SASSI UZ -> ANSYS UY
    assert "EN,2,2,5,3,3" in L
    back = convert_cdb(read_cdb_text(res.text), 9.81).model
    for n in m.nodes:
        np.testing.assert_allclose(back.node_global(n), m.node_global(n), atol=1e-15)
    assert back.nodes[1].fix == [1, 0, 1, 0, 0, 0] and back.nodes[5].fix == [0, 0, 1, 0, 0, 0]
    _same_matrices(m, back)
