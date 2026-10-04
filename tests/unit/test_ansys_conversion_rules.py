"""CONVERT,ANSYS / ANSYS: conversion rules added after the adversarial review of the ANSYS interfaces
(requirements 3.4.K; D-ANS-01, D-ANS-03, D-ANS-06; spec 04 sections 11.1-11.4; R1 section 1.2).

* rotations are fixed only at nodes that no element other than an axial-only LINK stiffens
  (a LINK tie at a shell corner must not clamp the plate);
* BEAM44 with the same release at both ends (pin-ended brace) becomes an equivalent valid
  SASSI member, checked against a closed-form frame;
* every option that changes the ANSYS stiffness or mass and is not converted is reported
  (KEYOPTs, real constants, SECCONTROL, MPDATA labels); added masses are lumped and reported;
* APDL parameters, expressions and component names do not stop the reader;
* the exact damping mapping of the export follows the model's CMODFORM.
"""
from __future__ import annotations

import math
import re

import numpy as np
import pytest
import scipy.sparse.linalg as spla

from sassi.conventions import cfactor
from sassi.io.ansys_cdb import PARAM_KEY, _num, convert_cdb, read_cdb_text
from sassi.io.apdl import ApdlOptions, ansys_damping, export_apdl
from sassi.prep import Interpreter, Kind
from sassi.verify.problems.vp_ansys import assemble_model, tip_deflection


def _conv(text: str, g: float = 9.81, damp=0.05):
    return convert_cdb(read_cdb_text(text), g, damp)


def _has(res, fragment: str, kind: str = "warning") -> bool:
    return any(fragment in t for k, t in res.messages if k == kind)


def _native(text: str):
    ui = Interpreter()
    ui.run_text(text)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    return ui.model


def _flex(model, node: int, dirs) -> np.ndarray:
    """Static flexibility block of ``node`` for global DOF directions ``dirs`` (1..6), real K."""
    dm, am = assemble_model(model)
    K = am.Ks.real.tocsc()
    F = np.zeros((dm.neq, len(dirs)))
    for c, d in enumerate(dirs):
        F[dm.eq(node, d), c] = 1.0
    U = spla.spsolve(K, F).reshape(dm.neq, len(dirs))
    return np.array([[U[dm.eq(node, d), c] for c in range(len(dirs))] for d in dirs])


# ======================================================================================
# LINK nodes: which rotations are fixed
# ======================================================================================
_PLATE = """/PREP7
ET,1,63
R,1,0.2
MP,EX,1,3e10
MP,PRXY,1,0.2
MP,DENS,1,2400
MP,DMPR,1,0.1
N,1,0,0,0
N,2,4,0,0
N,3,4,4,0
N,4,0,4,0
N,5,2,2,-3
TYPE,1
REAL,1
MAT,1
EN,1,1,2,3,4
{support}
D,1,ALL,0
D,2,ALL,0
D,4,ALL,0
D,5,ALL,0
D,3,ROTZ,0
"""
_EA_TIE = 3.0e10 * 0.01                    # tie of the plate material, A = 0.01


def test_link_at_a_shell_corner_acts_as_an_axial_spring_only():
    """SHELL63 plate supported at one corner by a LINK8 tie: the plate keeps its corner rotations,
    and the tie adds exactly its axial stiffness -- the same plate with the tie replaced by an
    oblique COMBIN14 of stiffness EA/L along the same line (a route without rotational DOFs,
    converted to a GENERAL k e e^T element) has the same corner flexibility."""
    L = math.sqrt(2.0 ** 2 + 2.0 ** 2 + 3.0 ** 2)
    link = _conv(_PLATE.format(support="ET,2,8\nR,2,0.01\nTYPE,2\nREAL,2\nEN,2,3,5"))
    spring = _conv(_PLATE.format(support=f"ET,2,14\nR,2,{_EA_TIE / L!r}\nTYPE,2\nREAL,2\nEN,2,3,5"))
    assert link.model.nodes[3].fix == [0, 0, 0, 0, 0, 1]                     # only the D,3,ROTZ
    f_link, f_spring = _flex(link.model, 3, (1, 2, 3, 4, 5)), _flex(spring.model, 3, (1, 2, 3, 4, 5))
    np.testing.assert_allclose(f_link, f_spring, rtol=1e-8, atol=1e-12 * np.abs(f_spring).max())
    # with the corner rotations clamped (the defect) the plate would be about 16 % stiffer
    clamped = link.model.copy()
    clamped.nodes[3].fix = [0, 0, 0, 1, 1, 1]
    assert tip_deflection(clamped, 3, 3) < 0.9 * tip_deflection(link.model, 3, 3)


def test_rotations_fixed_only_where_no_element_stiffens_them():
    """A chain of LINK180 ties 1-2-...-9: the rotations of LINK-only, LINK + SOLID, LINK + translational
    COMBIN14 and LINK + translational MATRIX27 nodes are fixed; LINK + rotational COMBIN14, LINK +
    rotary MASS21, LINK + BEAM4 and LINK + rotational MATRIX27 nodes keep them free."""
    txt = """ET,1,180
SECTYPE,1,LINK
SECDATA,0.01
ET,2,45
ET,3,14
KEYOPT,3,2,1
ET,4,14
KEYOPT,4,2,4
ET,5,21
ET,6,27
KEYOPT,6,3,4
ET,7,4
R,3,1e5
R,4,1e5
R,5,1,1,1,0.1,0.1,0.1
R,6,1000
R,7,0.01,1e-5,1e-5,0,0,0
RMORE,0,2e-5
R,8,0,0,0,0,0,0
RMORE,0,0,0,0,0,0
RMORE,0,0,0,0,0,0
RMORE,0,0,0,0,0,0
RMORE,0,0,0,0,0,0
RMORE,0,0,0,1000,0,0
MP,EX,1,2e11
MP,PRXY,1,0.3
MP,DENS,1,7850
MP,DMPR,1,0.04
N,1,0,0,0
N,2,1,0,0
N,3,2,0,0
N,4,3,0,0
N,5,4,0,0
N,6,5,0,0
N,7,6,0,0
N,8,7,0,0
N,9,8,0,0
N,10,1,1,0
N,11,1,1,1
N,12,0,1,1
N,13,0,0,1
N,14,1,0,1
N,15,9,0,0
N,16,0,1,0
N,17,1,-1,0
TYPE,1
SECNUM,1
MAT,1
EN,1,1,2
EN,2,2,3
EN,3,3,4
EN,4,4,5
EN,5,5,6
EN,6,6,7
EN,7,7,8
EN,8,8,9
TYPE,2
EN,9,1,17,10,16,13,14,11,12
TYPE,3
REAL,3
EN,10,3,15
TYPE,4
REAL,4
EN,11,4,15
TYPE,5
REAL,5
EN,12,5
TYPE,6
REAL,6
EN,13,6,15
REAL,8
EN,14,9,15
TYPE,7
REAL,7
EN,15,7,8,11
"""
    r = _conv(txt)
    fix = {n: r.model.nodes[n].fix[3:] for n in range(1, 10)}
    assert fix[1] == [1, 1, 1]                       # LINK + SOLID45
    assert fix[2] == [1, 1, 1]                       # LINK only
    assert fix[3] == [1, 1, 1]                       # LINK + translational COMBIN14 (KEYOPT(2) = 1)
    assert fix[4] == [0, 0, 0]                       # LINK + rotational COMBIN14 (KEYOPT(2) = 4)
    assert fix[5] == [0, 0, 0]                       # LINK + MASS21 with rotary inertia
    assert fix[6] == [1, 1, 1]                       # LINK + MATRIX27 with C1 (UX) only
    assert fix[7] == [0, 0, 0] and fix[8] == [0, 0, 0]   # LINK + BEAM4
    assert fix[9] == [0, 0, 0]                       # LINK + MATRIX27 with C34 (ROTX)
    assert _has(r, "ROTX ROTY ROTZ fixed at 4 nodes that no other element stiffens in rotation (1, 2, 3, 6)",
                "info")


# ======================================================================================
# BEAM44 released at both ends
# ======================================================================================
E7, NU7 = 2.0e11, 0.3
H_COL, B_BR = 4.0, 3.0
A_COL, I_COL, J_COL = 0.02, 2.0e-4, 3.0e-4
A_BR = 5.0e-3


def _braced_cantilever(k78: int) -> str:
    """Vertical BEAM4 cantilever 1-2 (fixed at 1) with a BEAM44 brace 3-2 (fixed at 3) in the X-Z
    plane; the brace has large inertias and KEYOPT(7) = KEYOPT(8) = ``k78``."""
    return f"""ET,1,4
ET,2,44
KEYOPT,2,7,{k78}
KEYOPT,2,8,{k78}
R,1,{A_COL},{I_COL},{I_COL},0,0,0
RMORE,0,{J_COL}
R,2,{A_BR},1e-2,2e-2,0,0,1e-2
RMORE,{A_BR},1e-2,2e-2,0,0,1e-2
MP,EX,1,{E7}
MP,PRXY,1,{NU7}
MP,DENS,1,7850
MP,DMPR,1,0.04
N,1,0,0,0
N,2,0,0,{H_COL}
N,3,{B_BR},0,0
N,4,0,1,{H_COL}
TYPE,1
REAL,1
EN,1,1,2
TYPE,2
REAL,2
EN,2,3,2,4
D,1,ALL,0
D,3,ALL,0
"""


def _column_k() -> np.ndarray:
    """Euler-Bernoulli tip stiffness of the fixed-base column in (UX, UZ, ROTY); ROTY = +dUX/dZ."""
    EI = E7 * I_COL
    H = H_COL
    return np.array([[12 * EI / H ** 3, 0.0, -6 * EI / H ** 2],
                     [0.0, E7 * A_COL / H, 0.0],
                     [-6 * EI / H ** 2, 0.0, 4 * EI / H]])


def test_pin_ended_brace_is_an_axial_member_closed_form():
    """KEYOPT(7) = KEYOPT(8) = 11 (ROTY, ROTZ released at both ends): the brace carries axial force
    only, so the tip stiffness is the column's plus EA/L e e^T of the brace (closed form), whatever
    the brace inertias.  SASSI gets it with I2 = I3 = 0 and no release (CHECK Error 10 avoided)."""
    r = _conv(_braced_cantilever(11))
    m = r.model
    br = m.groups[2].elements[1]
    s = m.sections[br.prop]
    assert br.ki == [0] * 6 and br.kj == [0] * 6
    assert (s.flex2, s.flex3) == (0.0, 0.0) and s.tors == pytest.approx(1e-2) and s.axial == A_BR
    assert _has(r, "at both ends") and _has(r, "Error 10")
    L = math.hypot(B_BR, H_COL)
    e = np.array([-B_BR, H_COL]) / L                           # brace axis in (X, Z), node 3 -> node 2
    K = _column_k()
    K[np.ix_([0, 1], [0, 1])] += E7 * A_BR / L * np.outer(e, e)
    np.testing.assert_allclose(_flex(m, 2, (1, 3, 5)), np.linalg.inv(K), rtol=1e-9)


def test_axial_and_moment_releases_at_both_ends_leave_the_column_alone():
    """KEYOPT(7) = KEYOPT(8) = 100011: UX (P1) released at both ends is kept at I only (one release
    already removes the axial force), ROTY/ROTZ give I2 = I3 = 0: the brace carries no in-plane
    force and the frame is the bare column (closed form)."""
    r = _conv(_braced_cantilever(100011))
    br = r.model.groups[2].elements[1]
    assert br.ki == [1, 0, 0, 0, 0, 0] and br.kj == [0] * 6
    assert _has(r, "P1 at I only")
    np.testing.assert_allclose(_flex(r.model, 2, (1, 3, 5)), np.linalg.inv(_column_k()), rtol=1e-9)


def test_shear_released_at_both_ends_closed_form():
    """KEYOPT(7) = KEYOPT(8) = 10000 releases UY (the in-plane shear: the brace's K node puts ANSYS z
    along global Y, so y lies in the X-Z plane) at both ends.  Without shear the bending moment is
    constant, so in-plane the brace adds its axial EA/L e e^T and a pure rotational stiffness
    E Izz / L at node 2 (node 3 is fixed) -- the closed form; kept as P3 released at I only."""
    r = _conv(_braced_cantilever(10000))
    br = r.model.groups[2].elements[1]
    assert br.ki == [0, 0, 1, 0, 0, 0] and br.kj == [0] * 6 and _has(r, "P3 at I only")
    L = math.hypot(B_BR, H_COL)
    e = np.array([-B_BR, H_COL]) / L
    K = _column_k()
    K[np.ix_([0, 1], [0, 1])] += E7 * A_BR / L * np.outer(e, e)
    K[2, 2] += E7 * 1e-2 / L                                    # IZZ of the brace = R2 = 1e-2
    np.testing.assert_allclose(_flex(r.model, 2, (1, 3, 5)), np.linalg.inv(K), rtol=1e-9)


def test_pin_ended_brace_model_passes_the_model_check():
    """The converted pin-ended BEAM44 model has no MODEL error in CHECK (in particular no Error 10,
    which the same member with KI = KJ = 000011 gets)."""
    def model_lines(model):
        ui = Interpreter()
        ui.models[ui.active_model] = model
        ui.execute("CHECK")
        return [t for k in (Kind.ERROR, Kind.WARNING, Kind.INFO) for t in ui.sink.texts(k)]
    m = _conv(_braced_cantilever(11)).model
    lines = model_lines(m)
    assert "MODEL: 0 errors, 0 warnings" in lines and not any("Error 10" in t for t in lines)
    bad = m.copy()
    bad.groups[2].elements[1].ki = [0, 0, 0, 0, 1, 1]
    bad.groups[2].elements[1].kj = [0, 0, 0, 0, 1, 1]
    assert any("Error 10" in t for t in model_lines(bad))


# ======================================================================================
# Unsupported options are reported; added masses are lumped and reported
# ======================================================================================
_QUAD = """ET,1,{et}
{kop}
SECTYPE,1,SHELL
SECDATA,0.2,1,0,3
{secctl}
R,1,{real}
MP,EX,1,3e10
MP,PRXY,1,0.2
MP,DENS,1,2400
MP,DMPR,1,0.1
N,1,0,0,0
N,2,2,0,0
N,3,2,3,0
N,4,0,3,0
TYPE,1
SECNUM,1
EN,1,1,2,3,4
"""


@pytest.mark.parametrize("et,kop,real,secctl,fragment", [
    ("181", "KEYOPT,1,1,1", "0.2", "", "KEYOPT(1) = 1 (membrane-only stiffness)"),
    ("281", "KEYOPT,1,1,1", "0.2", "", "KEYOPT(1) = 1 (membrane-only stiffness)"),
    ("63", "KEYOPT,1,1,2", "0.2", "", "KEYOPT(1) = 2 (membrane-only or bending-only"),
    ("63", "", "0.2,0,0,0,1e6", "", "elastic foundation stiffness EFS (R5)"),
    ("63", "", "0.2,0,0,0,0,0\nRMORE,0.5", "", "RMI (R7) = 0.5"),
    ("63", "", "0.2,0,0,0,0,0\nRMORE,0,0,0,7", "", "real constants R10 are not converted"),
    ("181", "", "0.2", "SECCONTROL,1e9,1e9,0", "transverse shear stiffness"),
    ("181", "KEYOPT,1,9,1", "0.2", "", "KEYOPT(9) = 1 not converted"),
])
def test_shell_options_that_change_stiffness_are_reported(et, kop, real, secctl, fragment):
    r = _conv(_QUAD.format(et=et, kop=kop, real=real, secctl=secctl))
    assert _has(r, fragment), [t for k, t in r.messages]


def test_shell_options_without_effect_are_silent():
    """SHELL181 KEYOPT(3) = 2 (integration) and KEYOPT(8) = 2 (layer data storage), SHELL63 RMI = 1
    (the default ratio) and THETA (R6, isotropic material) give no warning."""
    r = _conv(_QUAD.format(et="181", kop="KEYOPT,1,3,2\nKEYOPT,1,8,2", real="0.2", secctl=""))
    assert not r.warnings, r.warnings
    r = _conv(_QUAD.format(et="63", kop="", real="0.2,0,0,0,0,30\nRMORE,1", secctl=""))
    assert not r.warnings, r.warnings


@pytest.mark.parametrize("et,real,secctl,label", [
    ("63", "0.2,0,0,0,0,0\nRMORE,0,0,0,0,0,0\nRMORE,0,0,0,0,0,150", "", "ADMSUA (R18)"),
    ("181", "0.2", "SECCONTROL,0,0,0,150", "SECCONTROL added mass per unit area"),
])
def test_shell_added_mass_per_area_is_lumped_like_the_shell_mass(et, real, secctl, label):
    """150 per unit area on a 2 x 3 plate: 150 * 6 / 4 = 225 at each corner (mass units)."""
    r = _conv(_QUAD.format(et=et, kop="", real=real, secctl=secctl))
    assert _has(r, label)
    assert r.model.tmass == {n: [225.0, 225.0, 225.0] for n in (1, 2, 3, 4)}
    assert all(r.model.mass_units[n] == 0 for n in (1, 2, 3, 4))


def _cant(et: str, extra: str, real: str, ne: int = 4, L: float = 8.0) -> str:
    lines = [f"ET,1,{et}", extra, real, "MP,EX,1,2e11", "MP,PRXY,1,0.3", "MP,DENS,1,7850", "MP,DMPR,1,0.04"]
    lines += [f"N,{i + 1},{L * i / ne!r},0,0" for i in range(ne + 1)]
    lines += [f"N,{ne + 2},0,0,5", "TYPE,1", "REAL,1", "MAT,1", "SECNUM,1"]
    lines += [f"EN,{i + 1},{i + 1},{i + 2}" + ("" if et in ("8", "180") else f",{ne + 2}") for i in range(ne)]
    return "\n".join(lines) + "\nD,1,ALL,0\n"


@pytest.mark.parametrize("et,extra,real", [
    ("4", "", "R,1,0.08,2e-4,1e-3,0,0,0\nRMORE,0,1e-4,0,0,0,50"),                 # BEAM4 ADDMAS = R12
    ("188", "SECTYPE,1,BEAM,RECT\nSECDATA,0.2,0.4\nSECCONTROL,0,0,50", ""),       # BEAM188 SECCONTROL ADDMAS
    ("180", "", "R,1,0.01,50"),                                                   # LINK180 R2 = ADDMAS
    ("180", "SECTYPE,1,LINK\nSECDATA,0.01\nSECCONTROL,50", ""),                   # LINK180 section ADDMAS
])
def test_added_mass_per_length_is_lumped_and_reported(et, extra, real):
    """ADDMAS = 50 per unit length on 4 elements of length 2: 50 at the ends, 100 inside."""
    r = _conv(_cant(et, extra, real))
    assert _has(r, "ADDMAS = 50 lumped as ADDMAS*L/2")
    tm = r.model.tmass
    assert [tm[n][0] for n in range(1, 6)] == pytest.approx([50.0, 100.0, 100.0, 100.0, 50.0])


@pytest.mark.parametrize("et,extra,real,fragment", [
    ("4", "", "R,1,0.08,2e-4,1e-3,0,0,0\nRMORE,1e-3,1e-4", "ISTRN (initial strain, a load) = 0.001"),
    ("4", "", "R,1,0.08,2e-4,1e-3,0,0,0\nRMORE,0,1e-4,0,0,10", "SPIN"),
    ("44", "", "R,1,0.08,2e-4,1e-3,0,0,1e-4\nRMORE,0.08,2e-4,1e-3,0,0,1e-4\nRMORE\nRMORE\nRMORE\nRMORE,0,0,0,0,0,7",
     "R36 (beyond R24) are not converted"),
    ("8", "", "R,1,0.01,1e-3", "ISTRN (initial strain, a load) = 0.001"),
    ("180", "", "R,1,0.01,0,1", "TENSKEY = 1 (tension or compression only)"),
    ("180", "KEYOPT,1,3,1", "R,1,0.01", "KEYOPT(3) = 1 (tension-only or compression-only)"),
    ("188", "KEYOPT,1,13,1\nSECTYPE,1,BEAM,RECT\nSECDATA,0.2,0.4", "", "KEYOPT(13) = 1 not converted"),
])
def test_beam_and_link_options_are_reported(et, extra, real, fragment):
    r = _conv(_cant(et, extra, real))
    assert _has(r, fragment), [t for k, t in r.messages]


@pytest.mark.parametrize("kop,fragment", [
    ("KEYOPT,1,6,1", "KEYOPT(6) = 1 (mixed u-P formulation)"),
    ("KEYOPT,1,15,1", "KEYOPT(15) = 1 not converted"),
])
def test_solid_keyopts_are_reported(kop, fragment):
    txt = f"""ET,1,185
{kop}
MP,EX,1,1e7
MP,PRXY,1,0.3
MP,DENS,1,2
MP,DMPR,1,0.1
N,1,0,0,1
N,2,1,0,1
N,3,1,1,1
N,4,0,1,1
N,5,0,0,2
N,6,1,0,2
N,7,1,1,2
N,8,0,1,2
EN,1,1,2,3,4,5,6,7,8
"""
    assert _has(_conv(txt), fragment)


_MAT = """ET,1,185
MP,EX,1,3e10
MP,PRXY,1,0.2
MP,DENS,1,2400
MP,DMPR,1,0.1
{extra}
N,1,0,0,1
N,2,1,0,1
N,3,1,1,1
N,4,0,1,1
N,5,0,0,2
N,6,1,0,2
N,7,1,1,2
N,8,0,1,2
EN,1,1,2,3,4,5,6,7,8
"""


@pytest.mark.parametrize("extra,fragment", [
    ("MP,GXY,1,5e9", "GXY = 5e+09 (isotropic value 1.25e+10)"),
    ("MP,EY,1,2e10", "EY = 2e+10"),
    ("MP,NUXY,1,0.25", "NUXY = 0.25"),
    ("MP,ALPX,1,1e-5\nMP,KXX,1,1.7", "MPDATA ALPX, KXX not converted"),
    ("MP,DMPS,1,0.08", "both DMPR and DMPS given"),
])
def test_material_labels_not_converted_are_reported(extra, fragment):
    r = _conv(_MAT.format(extra=extra))
    assert _has(r, fragment), r.warnings
    assert r.model.materials[1].val1 == 3e10 and r.model.materials[1].val2 == 0.2


def test_isotropic_data_written_in_full_is_silent():
    """EY = EZ = EX, PRYZ = PRXZ = PRXY and G = EX/(2(1 + nu)) (e.g. Workbench output) is accepted."""
    extra = "MP,EY,1,3e10\nMP,EZ,1,3e10\nMP,PRYZ,1,0.2\nMP,PRXZ,1,0.2\nMP,GXY,1,1.25e10\nMP,GYZ,1,1.25e10\nMP,GXZ,1,1.25e10"
    r = _conv(_MAT.format(extra=extra))
    assert not r.warnings, r.warnings


# ======================================================================================
# APDL parameters, expressions and component names
# ======================================================================================
def test_fortran_double_literals_and_parameter_names():
    assert _num("1.5D+03") == 1500.0 and _num("-2.d-1") == -0.2 and _num(" 7 ") == 7.0
    with pytest.raises(ValueError, match="'tid'"):                       # no D -> E mangling of names
        _num("tid")


def test_commands_with_parameters_are_counted_and_reported():
    txt = """/PREP7
*SET,tid,4
ET,1,185
ET,tid,170
MP,EX,1,1e7
MP,PRXY,1,0.3
MP,DENS,1,2
MP,DMPR,1,0.1
N,1,0,0,1
N,2,L/2,0,1
N,3,1,1,1
N,4,0,1,1
N,5,0,0,2
N,6,1,0,2
N,7,1,1,2
N,8,0,1,2
N,9,0.5,0,1
EN,1,1,9,3,4,5,6,7,8
"""
    c = read_cdb_text(txt)
    assert c.skipped[PARAM_KEY] == 2 and c.unevaluated == ["ET,tid,170", "N,2,L/2,0,1"]
    assert sorted(c.nodes) == [1, 3, 4, 5, 6, 7, 8, 9]                    # reading went on
    r = convert_cdb(c, 9.81)
    assert _has(r, "APDL parameters, expressions or names") and _has(r, "'ET,tid,170'")


def test_d_on_a_node_component_is_expanded():
    """``D,BASE,ALL,0`` before the CMBLOCK that defines BASE (CDWRITE order) fixes its nodes; a D on
    an element component or an undefined name is reported."""
    txt = """ET,1,185
MP,EX,1,1e7
MP,PRXY,1,0.3
MP,DENS,1,2
MP,DMPR,1,0.1
N,1,0,0,1
N,2,1,0,1
N,3,1,1,1
N,4,0,1,1
N,5,0,0,2
N,6,1,0,2
N,7,1,1,2
N,8,0,1,2
EN,1,1,2,3,4,5,6,7,8
D,BASE,UZ,0
D,BASE,UX,0,,,,UY
D,TOPEL,ALL,0
D,NOPE,ALL,0
CMBLOCK,BASE,NODE,       2
(8i10)
         1        -4
CMBLOCK,TOPEL,ELEM,       1
(8i10)
         1
"""
    r = _conv(txt)
    assert all(r.model.nodes[n].fix[:3] == [1, 1, 1] for n in (1, 2, 3, 4))
    assert not any(r.model.nodes[n].fix[:3] != [0, 0, 0] for n in (5, 6, 7, 8))
    assert _has(r, "on 'TOPEL' not converted") and _has(r, "on 'NOPE' not converted")


def test_malformed_block_line_is_reported_not_fatal():
    txt = """NBLOCK,6,SOLID,         3,         3
(3i9,6e21.13e3)
        1        0        0
        2        0        0 1.0000000000000E+00x
        3        0        0 2.0000000000000E+000
N,R5.3,LOC,       -1,
"""
    c = read_cdb_text(txt)
    assert sorted(c.nodes) == [1, 3] and any("malformed block data line" in t for t in c.notes)


def test_convert_command_with_workbench_style_input(tmp_path):
    """Through the interpreter: no 'internal error', the parameter lines are warned about."""
    p = tmp_path / "ds.cdb"
    p.write_text("/PREP7\n*SET,tid,4\nET,1,185\nET,tid,170\nMP,EX,1,1e7\nMP,PRXY,1,0.3\nMP,DENS,1,2\n"
                 "N,1,0,0,1\nN,2,1,0,1\nN,3,1,1,1\nN,4,0,1,1\nN,5,0,0,2\nN,6,1,0,2\nN,7,1,1,2\nN,8,0,1,2\n"
                 "EN,1,1,2,3,4,5,6,7,8\n")
    ui = Interpreter(cwd=tmp_path)
    ok = ui.execute(f"CONVERT,ANSYS,,{p},9.81")
    errs = ui.sink.texts(Kind.ERROR)
    assert ok and not errs, errs
    assert any("'ET,tid,170'" in t for t in ui.sink.texts(Kind.WARNING))
    assert ui.model.n_elements() == 1


# ======================================================================================
# Export: CMODFORM and round trips with axial-only members at shell nodes
# ======================================================================================
@pytest.mark.parametrize("beta", [0.02, 0.1, 0.3])
def test_exact_mapping_for_both_complex_modulus_forms(beta):
    E = 2.5e7
    for form in (0, 1):
        Ea, g = ansys_damping(E, beta, "exact", form)
        assert Ea * (1 + 1j * g) == pytest.approx(E * cfactor(beta, form), rel=1e-14)
    assert ansys_damping(E, beta, "exact", 1) == (E, 2 * beta)


def test_export_exact_damping_follows_cmodform():
    txt = """GRAVITY,32.2
CMODFORM,{form}
N,1,0,0,0
N,2,1,0,0
N,3,0,1,0
M,1,1e6,0.25,0.15,0.1,0.1,1
R,1,1,0,0,0.1,0.1,0.1
GROUP,1,BEAMS
E,1,1,2,3
"""
    for form in (0, 1):
        text = export_apdl(_native(txt.format(form=form)), ApdlOptions(damping="exact")).text
        ex = float(re.search(r"^MP,EX,1,([^\s!]+)", text, re.M).group(1))
        g = float(re.search(r"^MP,DMPR,1,([^\s!]+)", text, re.M).group(1))
        assert ex * (1 + 1j * g) == pytest.approx(1e6 * cfactor(0.1, form), rel=1e-12)
        assert ("CMODFORM,1" in text) == (form == 1)


_SHELL_TIE = """GRAVITY,9.81
N,1,0,0,0
N,2,4,0,0
N,3,4,4,0
N,4,0,4,0
N,5,2,2,-3
M,1,3e10,0.2,24000,0.05,0.05,1
M,2,2e11,0.3,77000,0.02,0.02,1
R,1,0.01,0,0,{J},0,0
GROUP,1,SHELL
E,1,1,2,3,4
THICK,1,1,1,0.2
GROUP,2,BEAMS
E,1,3,5,1
MSET,1,1,1,2
D,1,1,1,1,ALL
D,2,2,1,1,ALL
D,4,4,1,1,ALL
D,5,5,1,1,ALL
D,3,3,1,1,ROTZ
"""


@pytest.mark.parametrize("modern", [False, True])
def test_shell_with_axial_member_round_trip(modern):
    """Native shell + axial-only member -> APDL (LINK8 / LINK180) -> CONVERT: same fixities and the
    same corner flexibility (the LINK J of 1e-8 A^2 equals the native J here)."""
    m = _native(_SHELL_TIE.format(J=repr(1e-8 * 0.01 ** 2)))
    back = convert_cdb(read_cdb_text(export_apdl(m, ApdlOptions(modern=modern)).text), 9.81).model
    assert back.nodes[3].fix == m.nodes[3].fix == [0, 0, 0, 0, 0, 1]
    np.testing.assert_allclose(_flex(back, 3, (1, 2, 3, 4, 5)), _flex(m, 3, (1, 2, 3, 4, 5)), rtol=1e-9)
