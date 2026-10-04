"""Adversarial review of the ANSYS interfaces (CONVERT,ANSYS / ANSYS / ANSYSREFORMAT / CONVERT,SSI).

Independent checks of sassi/io/ansys_cdb.py, sassi/io/apdl.py and sassi/prep/commands/conversion.py
against requirements 3.4.K, D-ANS-01..08, R1 1.2 and spec 04 sections 11-12:

* closed-form Euler-Bernoulli / Timoshenko references derived here (not taken from the package);
* limiting cases and symmetry (swapped bending planes, oblique springs, 2-D rotation);
* round trips model -> APDL -> model on models the package's own tests do not cover
  (shell + axial-only members, pin-ended BEAM44, CMODFORM 1);
* robustness of the reader on APDL input that uses parameters.

Tests named ``test_defect_*`` reproduce defects found in the review; they fail until the defect is
fixed (they assert the correct behaviour).
"""
from __future__ import annotations

import math
import re

import numpy as np
import pytest
import scipy.sparse.linalg as spla

from sassi.conventions import cfactor
from sassi.io.ansys_cdb import ansys_beam_axes, convert_cdb, read_cdb_text
from sassi.io.apdl import ApdlOptions, cmblock_lines, export_apdl
from sassi.prep import Interpreter, Kind
from sassi.verify.problems.vp_ansys import (DATA, assemble_model, fixed_base_frequencies, tip_deflection)

E_ST, NU_ST = 2.0e11, 0.3
G_ST = E_ST / (2.0 * (1.0 + NU_ST))


def _conv(text: str, g: float = 9.81, damp=None):
    return convert_cdb(read_cdb_text(text), g, damp)


def _native(text: str):
    ui = Interpreter()
    ui.run_text(text)
    errs = ui.sink.texts(Kind.ERROR)
    assert not errs, errs
    return ui.model


def _flexibility(model, node: int, dirs=(2, 3)) -> np.ndarray:
    """Static flexibility block of ``node`` (rows/cols = global directions ``dirs``), real K."""
    dm, am = assemble_model(model)
    K = am.Ks.real.tocsc()
    F = np.zeros((dm.neq, len(dirs)))
    for c, d in enumerate(dirs):
        F[dm.eq(node, d), c] = 1.0
    U = spla.spsolve(K, F)
    U = U.reshape(dm.neq, len(dirs))
    return np.array([[U[dm.eq(node, d), c] for c in range(len(dirs))] for d in dirs])


def _cantilever_apdl(section_lines, ne=10, L=10.0, k_node=(0.0, 0.0, 5.0), et="188", real=None):
    lines = ["/PREP7", f"ET,1,{et}", "MP,EX,1,2e11", "MP,PRXY,1,0.3", "MP,DENS,1,7850"] + list(section_lines)
    if real is not None:
        lines.append(real)
    for i in range(ne + 1):
        lines.append(f"N,{i + 1},{L * i / ne!r},0,0")
    if k_node is not None:
        lines.append(f"N,{ne + 2},{k_node[0]!r},{k_node[1]!r},{k_node[2]!r}")
    lines += ["TYPE,1", "MAT,1", "REAL,1", "SECNUM,1"]
    for i in range(ne):
        lines.append(f"E,{i + 1},{i + 2}" + (f",{ne + 2}" if k_node is not None else ""))
    lines.append("D,1,ALL,0")
    return "\n".join(lines) + "\n"


# ======================================================================================
# Beam axes and sections: independent references
# ======================================================================================
def test_asec_product_of_inertia_matches_unsymmetric_bending_theory():
    """ASEC with Iyz != 0 (principal-axis conversion): the Euler-Bernoulli tip flexibility of a
    cantilever is L^3/(3E) C^-1 with C = [[Izz, Iyz], [Iyz, Iyy]] in the ANSYS (y, z) axes
    (unsymmetric bending, Iyz = int y z dA).  Derived from the strain energy, independent of the
    eigen-decomposition used by the converter."""
    Iyy, Izz, Iyz, L = 1.0e-3, 2.0e-4, 3.0e-4, 10.0
    txt = _cantilever_apdl(["SECTYPE,1,BEAM,ASEC", f"SECDATA,0.08,{Iyy!r},{Iyz!r},{Izz!r},0,1e-4"], L=L)
    m = _conv(txt).model
    flex = _flexibility(m, 11)
    C = np.array([[Izz, Iyz], [Iyz, Iyy]])           # beam along X, K above: y = global Y, z = global Z
    ref = L ** 3 / (3.0 * E_ST) * np.linalg.inv(C)
    np.testing.assert_allclose(flex, ref, rtol=1e-9, atol=1e-12 * np.abs(ref).max())


def test_asec_principal_axes_are_invariant_under_a_rotated_k_node():
    """The same ASEC member with its K node rotated by 90 degrees about the member axis (z -> -y)
    and the section data re-expressed in the new axes (Iyy <-> Izz, Iyz -> -Iyz) is the same
    structure: identical global flexibility (symmetry check of D-ANS-03)."""
    Iyy, Izz, Iyz = 1.0e-3, 2.0e-4, 3.0e-4
    a = _conv(_cantilever_apdl(["SECTYPE,1,BEAM,ASEC", f"SECDATA,0.08,{Iyy!r},{Iyz!r},{Izz!r},0,1e-4"]))
    # K along -Y: ANSYS z' = -Y, y' = z' x x = -Y x X = Z; so y' = z, z' = -y -> Iy'y' = Izz, Iz'z' = Iyy,
    # Iy'z' = int z (-y) dA = -Iyz
    b = _conv(_cantilever_apdl(["SECTYPE,1,BEAM,ASEC", f"SECDATA,0.08,{Izz!r},{-Iyz!r},{Iyy!r},0,1e-4"],
                               k_node=(0.0, -5.0, 0.0)))
    np.testing.assert_allclose(_flexibility(a.model, 11), _flexibility(b.model, 11), rtol=1e-9)


def test_default_orientation_tip_deflections_use_the_ansys_section_axes():
    """VP-A1 checks the default-orientation (no K node) file only through the frequency list, which
    cannot see a swap of the bending planes (see the next test).  Close the gap: tip deflections of
    the no-K-node cantilever along Y and Z must use Izz = H B^3/12 and Iyy = B H^3/12."""
    m = convert_cdb(read_cdb_text((DATA / "beam188_default_orientation.cdb").read_text()), 9.81).model
    L, B, H = 20.0, 0.2, 0.4
    A = B * H
    for d, I in ((2, H * B ** 3 / 12.0), (3, B * H ** 3 / 12.0)):
        ref = L ** 3 / (3.0 * E_ST * I) + L / (G_ST * A / 1.2)
        assert tip_deflection(m, 21, d) == pytest.approx(ref, rel=1e-9)


def test_frequency_comparison_is_blind_to_swapped_bending_planes():
    """Documents a limitation of VP-A1 checks (1) and (4): swapping I2 and I3 (the classical axis
    mapping error) leaves the sorted frequency list of the straight cantilever unchanged, while
    the static check sees it.  Frequencies alone do not verify D-ANS-03."""
    m = convert_cdb(read_cdb_text((DATA / "beam188_cantilever.cdb").read_text()), 9.81).model
    sw = m.copy()
    for s in sw.sections.values():
        s.flex2, s.flex3 = s.flex3, s.flex2
        s.shear2, s.shear3 = s.shear3, s.shear2
    f0, f1 = fixed_base_frequencies(m)[:8], fixed_base_frequencies(sw)[:8]
    np.testing.assert_allclose(f0, f1, rtol=1e-8)
    assert tip_deflection(m, 21, 2) != pytest.approx(tip_deflection(sw, 21, 2), rel=1e-3)


def test_beam4_real_constants_equal_beam188_rect():
    """BEAM4 with AREA, IZZ, IYY, IXX, SHEARZ, SHEARY of the RECT section (real-constant positions
    1, 2, 3, 8, 9, 10) gives the same stiffness as BEAM188 RECT B x H (independent route)."""
    from sassi.model.materials import section_rectangle
    B, H = 0.2, 0.4
    s = section_rectangle(B, H)
    A, Izz, Iyy, J = B * H, H * B ** 3 / 12, B * H ** 3 / 12, s["tors"]
    real = (f"R,1,{A!r},{Izz!r},{Iyy!r},0,0,0\n"
            f"RMORE,0,{J!r},1.2,1.2,0,0")                    # ISTRN, IXX, SHEARZ, SHEARY, SPIN, ADDMAS
    m4 = _conv(_cantilever_apdl([], et="4", real=real)).model
    m188 = _conv(_cantilever_apdl(["SECTYPE,1,BEAM,RECT", f"SECDATA,{B!r},{H!r}"])).model
    np.testing.assert_allclose(_flexibility(m4, 11), _flexibility(m188, 11), rtol=1e-9)
    # and Timoshenko closed form: P L^3/(3 E I) + P L/(G A/1.2)
    L = 10.0
    ref = np.diag([L ** 3 / (3 * E_ST * Izz) + L / (G_ST * A / 1.2), L ** 3 / (3 * E_ST * Iyy) + L / (G_ST * A / 1.2)])
    np.testing.assert_allclose(_flexibility(m4, 11), ref, rtol=1e-9, atol=1e-15)


def test_ansys_axes_are_right_handed_and_k_in_xz_plane():
    rng = np.random.default_rng(3)
    for _ in range(20):
        xi, xj, xk = rng.normal(size=(3, 3))
        A, L = ansys_beam_axes(xi, xj, xk)
        np.testing.assert_allclose(A @ A.T, np.eye(3), atol=1e-12)
        assert np.linalg.det(A) == pytest.approx(1.0)
        v = xk - xi
        assert abs(v @ A[1]) < 1e-10 and v @ A[2] > 0                   # K in x-z plane, z towards K
        A0, _ = ansys_beam_axes(xi, xj)
        assert abs(A0[1, 2]) < 1e-12                                       # default y parallel to X-Y plane


# ======================================================================================
# Springs, masses, matrices, 2-D rotation
# ======================================================================================
def test_oblique_combin14_is_k_e_eT_in_global_axes():
    k = 5.0e6
    txt = """ET,1,14
R,1,5e6
N,1,0,0,0
N,2,1,2,2
EN,1,1,2
"""
    m = _conv(txt).model
    p = m.matrices[m.groups[1].elements[1].prop]
    e = np.array([1.0, 2.0, 2.0]) / 3.0
    K = p.full("R")
    np.testing.assert_allclose(K[:3, :3], k * np.outer(e, e), rtol=1e-14)
    np.testing.assert_allclose(K[:3, 6:9], -k * np.outer(e, e), rtol=1e-14)
    assert not np.any(K[3:6, :]) and np.allclose(K, K.T)


def test_matrix27_upper_triangle_row_order():
    vals = list(range(1, 79))                                         # C1..C78 distinct
    real = "R,1," + ",".join(str(v) for v in vals[:6]) + "\n"
    for k in range(6, 78, 6):
        real += "RMORE," + ",".join(str(v) for v in vals[k:k + 6]) + "\n"
    txt = "ET,1,27\nKEYOPT,1,3,4\n" + real + "N,1,0,0,0\nN,2,1,0,0\nEN,1,1,2\n"
    m = _conv(txt).model
    K = m.matrices[m.groups[1].elements[1].prop].full("R")
    pos = 0
    for i in range(12):
        for j in range(i, 12):
            pos += 1
            assert K[i, j] == pos and K[j, i] == pos


def test_two_d_plane_model_dof_mapping():
    """ANSYS 2-D (X-Y) -> SASSI X-Z: a vertical (ANSYS UY) spring, mass and fixity become UZ."""
    txt = """ET,1,182
KEYOPT,1,3,2
ET,2,14
KEYOPT,2,2,2
ET,3,21
KEYOPT,3,3,4
R,2,1000
R,3,7
MP,EX,1,1e7
MP,PRXY,1,0.3
MP,DENS,1,2
N,1,0,0
N,2,1,0
N,3,1,1
N,4,0,1
N,5,0,2
TYPE,1
MAT,1
EN,1,1,2,3,4
TYPE,2
REAL,2
EN,2,4,5
TYPE,3
REAL,3
EN,3,3
D,1,UY,0
"""
    r = _conv(txt)
    m = r.model
    np.testing.assert_allclose(m.node_global(3), [1, 0, 1])
    sp = m.springs[m.groups[2].elements[1].prop]
    assert sp.k == (0.0, 0.0, 1000.0, 0.0, 0.0, 0.0)
    assert m.tmass[3] == [7.0, 0.0, 7.0] and m.mass_units[3] == 0
    assert m.nodes[1].fix == [0, 0, 1, 0, 0, 0]
    from sassi.elements.plane import centroid_jacobian
    xyz = np.array([m.node_global(n) for n in m.groups[1].elements[1].nodes])
    assert float(np.ravel(centroid_jacobian(xyz))[0]) > 0.0               # CCW kept: no reordering


def test_cmblock_writer_reader_round_trip():
    ids = [3, 4, 5, 6, 9, 11, 12, 100]
    text = "\n".join(cmblock_lines("SSI_INT", "NODE", ids)) + "\n"
    c = read_cdb_text(text)
    assert c.components["SSI_INT"] == ("NODE", ids)


def test_apdl_beam44_release_round_trip_and_stiffness():
    """A SASSI beam with M3 released at I -> BEAM44 KEYOPT(7) -> CONVERT: same KI, same K."""
    txt = """GRAVITY,9.81
N,1,0,0,0
N,2,3,0,0
N,3,0,1,0
M,1,1e7,0.25,9.81,0,0,1
R,1,0.02,0,0,1e-4,2e-4,5e-4
GROUP,1,BEAMS
E,1,1,2,3
KI,1,1,1,0,0,0,0,0,1
D,1,1,1,1,ALL
"""
    m = _native(txt)
    back = _conv(export_apdl(m).text).model
    assert back.groups[1].elements[1].ki == [0, 0, 0, 0, 0, 1]
    np.testing.assert_allclose(_flexibility(back, 2, (1, 2, 3)), _flexibility(m, 2, (1, 2, 3)), rtol=1e-10)


# ======================================================================================
# Defects found in the review (assert the correct behaviour)
# ======================================================================================
SHELL_LINK_APDL = """/PREP7
ET,1,63
ET,2,8
R,1,0.2
R,2,0.01
MP,EX,1,3e10
MP,PRXY,1,0.2
MP,DENS,1,2400
N,1,0,0,0
N,2,4,0,0
N,3,4,4,0
N,4,0,4,0
N,5,2,2,-3
TYPE,1
REAL,1
MAT,1
EN,1,1,2,3,4
TYPE,2
REAL,2
EN,2,3,5
D,1,ALL,0
D,2,ALL,0
D,4,ALL,0
D,5,ALL,0
"""


def test_defect_link_at_a_shell_node_must_not_clamp_the_shell_rotations():
    """LINK8 tie attached to a SHELL63 corner: the shell gives that node bending (rotational)
    stiffness, so its rotations are not 'truss-only' and must stay free.  The converter fixes
    ROTX ROTY ROTZ there (SHELL nodes are never added to ``rot_nodes``), clamping the plate."""
    r = _conv(SHELL_LINK_APDL)
    assert r.model.nodes[3].fix[3:5] == [0, 0], r.model.nodes[3].fix


def test_defect_shell_plus_axial_member_apdl_round_trip_changes_stiffness():
    """Native model (shell + I2 = I3 = 0 member) -> ANSYS (LINK8) -> CONVERT must keep the
    structure; the clamped shell rotation at the link node makes the plate ~16 % stiffer."""
    txt = """GRAVITY,9.81
N,1,0,0,0
N,2,4,0,0
N,3,4,4,0
N,4,0,4,0
N,5,2,2,-3
M,1,3e10,0.2,24000,0.05,0.05,1
R,1,0.01,0,0,1e-10,0,0
GROUP,1,SHELL
E,1,1,2,3,4
THICK,1,1,1,0.2
GROUP,2,BEAMS
E,1,3,5,1
D,1,1,1,1,ALL
D,2,2,1,1,ALL
D,4,4,1,1,ALL
D,5,5,1,1,ALL
D,3,3,1,1,ROTZ
"""
    m = _native(txt)
    back = _conv(export_apdl(m).text).model
    assert back.nodes[3].fix == m.nodes[3].fix
    assert tip_deflection(back, 3, 3) == pytest.approx(tip_deflection(m, 3, 3), rel=1e-8)


def test_defect_pin_ended_beam44_converts_to_a_valid_model_or_is_reported():
    """BEAM44 with the same moment releases at both ends (KEYOPT(7) = KEYOPT(8) = 11, a pin-ended
    brace) is ordinary ANSYS input.  SASSI rejects the same release at I and J (CHECK Error 10), so
    the converter must either map it (no bending stiffness in the released planes, e.g. I = 0) or
    report it; today it converts silently into a model that fails CHECK and cannot be assembled."""
    txt = """ET,1,BEAM44
KEYOPT,1,7,11
KEYOPT,1,8,11
R,1,0.02,2e-4,5e-4,0,0,1e-4
RMORE,0.02,2e-4,5e-4,0,0,1e-4
N,1,0,0,0
N,2,3,0,0
N,3,0,0,1
MP,EX,1,1e7
MP,PRXY,1,0.25
MP,DENS,1,1
EN,1,1,2,3
D,1,ALL,0
"""
    r = _conv(txt)
    reported = any(re.search(r"both|Error 10|I and J", t) for t in r.warnings)
    try:
        assemble_model(r.model)
        assembles = True
    except Exception:
        assembles = False
    assert reported or assembles


def test_defect_exact_damping_mapping_ignores_cmodform():
    """``<dmap>`` = 1 promises E_ANSYS (1 + i g) = the SASSI complex modulus.  With CMODFORM,1 the
    model's modulus is E (1 + 2 i beta), so the exact mapping is EX = E, DMPR = 2 beta; the export
    writes E (1 - 2 beta^2) and g = 2 beta sqrt(1 - beta^2)/(1 - 2 beta^2) regardless."""
    txt = """GRAVITY,32.2
CMODFORM,1
N,1,0,0,0
N,2,1,0,0
N,3,0,1,0
M,1,1e6,0.25,0.15,0.1,0.1,1
R,1,1,0,0,0.1,0.1,0.1
GROUP,1,BEAMS
E,1,1,2,3
"""
    m = _native(txt)
    text = export_apdl(m, ApdlOptions(damping="exact")).text
    ex = float(re.search(r"^MP,EX,1,([^\s!]+)", text, re.M).group(1))
    g = float(re.search(r"^MP,DMPR,1,([^\s!]+)", text, re.M).group(1))
    beta = 0.1
    assert ex * (1 + 1j * g) == pytest.approx(1e6 * cfactor(beta, 1), rel=1e-12)


@pytest.mark.parametrize("snippet", [
    "*SET,tid,4\nET,tid,170\n",                 # Workbench contact definitions
    "N,1,0,0,0\nD,SSI_BASE,ALL,0\n",            # D on a component name
    "N,1,0,0,0\nN,2,L/2,0,0\n",                 # expression
])
def test_defect_reader_reports_instead_of_crashing_on_apdl_parameters(snippet):
    """Unsupported APDL (parameters, component names, expressions) must be reported, not raise
    (the CONVERT command shows 'internal error: could not convert string to float: tie' -- the
    D->E exponent substitution even mangles the parameter name)."""
    c = read_cdb_text("/PREP7\nET,1,185\n" + snippet)
    assert c.notes or c.unknown or c.skipped


def test_defect_added_mass_is_reported():
    """Spec 04 11.4: ADDMAS lumped as ADDMAS L/2 is an optional extension *with a warning*; the
    converter adds the masses silently."""
    real = "R,1,0.08,2.6667e-4,1.0667e-3,0,0,0\nRMORE,0,1e-4,0,0,0,50"     # ADDMAS = 50 (R12)
    r = _conv(_cantilever_apdl([], et="4", real=real))
    assert r.model.tmass                                                   # masses were added ...
    assert any("ADDMAS" in t or "added mass" in t.lower() for t in r.warnings + r.infos)   # ... and reported


_SHELL_BASE = """ET,1,{et}
{kop}
SECTYPE,1,SHELL
SECDATA,0.2,1,0,3
R,1,{real}
MP,EX,1,3e10
MP,PRXY,1,0.2
MP,DENS,1,2400
{extra}
N,1,0,0,0
N,2,1,0,0
N,3,1,1,0
N,4,0,1,0
TYPE,1
SECNUM,1
EN,1,1,2,3,4
"""


@pytest.mark.parametrize("et,kop,real,extra,what", [
    ("181", "KEYOPT,1,1,1", "0.2", "", "SHELL181 KEYOPT(1) = 1 (membrane only)"),
    ("63", "", "0.2,0,0,0,1e6", "", "SHELL63 EFS (R5, elastic foundation stiffness)"),
    ("63", "", "0.2,0,0,0,0,0\nRMORE,0.5", "", "SHELL63 RMI (R7, bending inertia ratio)"),
    ("63", "", "0.2", "MP,GXY,1,5e9", "GXY inconsistent with EX/PRXY (D-ANS-01: others ignored with warning)"),
])
def test_defect_unsupported_options_that_change_stiffness_are_reported(et, kop, real, extra, what):
    """Spec 04 11.1: 'emit a warning for every unsupported element type, KEYOPT, real-constant
    layout or command it skips'.  These inputs change the ANSYS stiffness but convert silently."""
    r = _conv(_SHELL_BASE.format(et=et, kop=kop, real=real, extra=extra), damp=0.05)
    other = [t for t in r.warnings if "no DMPR" not in t]
    assert other, f"{what}: converted without a warning"
