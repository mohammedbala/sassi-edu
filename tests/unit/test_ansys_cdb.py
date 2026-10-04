"""CONVERT,ANSYS: the ``.cdb`` reader and the ANSYS -> SASSI-EDU mapping (requirements 3.4.K,
D-ANS-01..05, D-ANS-07; spec 04 sections 11.1-11.6).

The sample files in ``tests/data/ansys`` are hand-written in the exact CDWRITE layout (ANSYS
2021 R1 style: NBLOCK ``(3i9,6e21.13e3)`` with trailing zero coordinates omitted, EBLOCK ``(19i9)``
with continuation lines, RLBLOCK ``(2i8,6g16.9)`` / ``(7g16.9)``, MPTEMP/MPDATA R5.0 records,
CMBLOCK ``(8i10)`` ranges, SECTYPE/SECDATA/SECBLOCK/SECOFFSET).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from sassi.io.ansys_cdb import (BeamProps, ansys_beam_axes, convert_cdb, element_name, parse_fortran_format,
                                read_cdb, read_cdb_text, read_fixed_numbers, rect_props, sassi_section_from_ansys,
                                split_fixed)
from sassi.verify.problems.vp_ansys import assemble_model, fixed_base_frequencies, tip_deflection

DATA = Path(__file__).resolve().parents[1] / "data" / "ansys"


def _conv(name, g=9.81, damp=0.05):
    return convert_cdb(read_cdb(DATA / name), g, damp)


def _texts(res, kind="warning"):
    return [t for k, t in res.messages if k == kind]


def _has(res, fragment, kind="warning"):
    return any(fragment in t for t in _texts(res, kind))


# ======================================================================================
# Fortran formats and fixed-width fields
# ======================================================================================
def test_fortran_format_parsing():
    assert parse_fortran_format("(3i9,6e21.13e3)") == [("i", 9)] * 3 + [("r", 21)] * 6
    assert parse_fortran_format("(2i8,6g16.9)") == [("i", 8)] * 2 + [("r", 16)] * 6
    assert parse_fortran_format("(19i9)") == [("i", 9)] * 19
    assert parse_fortran_format("(8i10)") == [("i", 10)] * 8
    assert parse_fortran_format("(1x,2i8)") == [("x", 1), ("i", 8), ("i", 8)]
    with pytest.raises(ValueError):
        parse_fortran_format("(3a9)")


def test_touching_e21_fields_are_split_by_width():
    """Negative e21.13e3 fields fill all 21 columns: whitespace splitting would merge them."""
    fmt = parse_fortran_format("(3i9,6e21.13e3)")
    line = "       17        0        0-1.0000000000000E+000-2.5000000000000E-001 3.0000000000000E+002"
    assert read_fixed_numbers(line, fmt) == [17.0, 0.0, 0.0, -1.0, -0.25, 300.0]
    assert split_fixed("        5        0        0", fmt) == ["5", "0", "0"]


def test_whitespace_fallback_for_free_format_data():
    fmt = parse_fortran_format("(3i8,6e16.9)")
    # a writer that ignored the widths: the fixed split gives no numbers -> whitespace split
    assert read_fixed_numbers("12 0 0 1.5 2.5 3.5", fmt) == [12.0, 0.0, 0.0, 1.5, 2.5, 3.5]


# ======================================================================================
# Reader: blocks and APDL commands
# ======================================================================================
def test_reader_blocks_of_the_solid_sample():
    c = read_cdb(DATA / "solids.cdb")
    assert len(c.nodes) == 44 and len(c.elements) == 6
    np.testing.assert_array_equal(c.nodes[1], [0.0, 0.0, 0.0])          # all coordinates omitted
    np.testing.assert_array_equal(c.nodes[2], [1.0, 0.0, 0.0])          # trailing zeros omitted
    assert c.node_rot == {12: (30.0, 0.0, 0.0)}
    assert c.etypes == {1: 185, 2: 45, 3: 186, 4: 187}
    e5 = c.elements[4]
    assert (e5.num, e5.type, e5.mat) == (5, 3, 2) and len(e5.nodes) == 20     # continuation line read
    assert c.materials[1] == {"EX": 3.0e7, "PRXY": 0.2, "DENS": 2.4, "DMPR": 0.1}
    assert c.materials[2]["NUXY"] == 0.25
    assert c.components["SSI_INT"] == ("NODE", [1, 2, 3, 4, 15, 16])         # range 1..4 from '1 -4'
    assert c.skipped == {"CP": 1}
    assert c.title.startswith("SASSI-EDU sample solids")


def test_reader_rlblock_continuation_and_sections():
    c = read_cdb(DATA / "frame_beams_links.cdb")
    assert c.reals[1] == pytest.approx([0.12, 1.6e-3, 9.0e-4, 0.4, 0.3, 90.0, 0.0, 2.0e-3, 1.2, 1.2])
    assert len(c.reals[2]) == 20 and c.reals[2][18:] == pytest.approx([1.2, 1.2])
    assert c.keyopt(2, 7) == 11
    s5 = c.sections[5]
    assert (s5.type, s5.subtype) == ("BEAM", "ASEC") and s5.data[:6] == pytest.approx([0.15, 3e-3, 0, 1e-3, 0, 2e-3])
    assert c.sections[6].type == "PIPE" and c.sections[3].type == "LINK"
    sh = read_cdb(DATA / "shells.cdb").sections
    assert sh[1].layers == [(0.5, 1, 0.0, 3)]                                 # SECDATA layer
    assert sh[2].layers == [(0.15, 1, 0.0, 3), (0.25, 2, 0.0, 3)]             # SECBLOCK layers


def test_reader_plain_apdl_commands():
    text = """
/PREP7
ET,1,SOLID185
ET,2,BEAM4,,,,,,
KEYO,2,7,1
MP,EX,1,2e8
MP,PRXY,1,0.3
MPDATA,DENS,1,1,7.85
R,3,0.1,1e-3,2e-3,0,0,0
RMORE,0,5e-3
N,1,0,0,0
N,2,1,0,0
N,3,1,1,0,30
TYPE,2
MAT,1
REAL,3
E,1,2,3
EN,10,2,3,1
TYPE,1
EN,11,1,2,3,3,1,2,3,3
EMORE,9,9
E,3,2
D,1,ALL,0,,3,2
D,2,UX,0,,,,UY,ROTZ
CE,1,0,1,UX,1
"""
    c = read_cdb_text(text)
    assert c.etypes == {1: 185, 2: 4}
    assert c.materials[1] == {"EX": 2e8, "PRXY": 0.3, "DENS": 7.85}
    assert c.reals[3] == pytest.approx([0.1, 1e-3, 2e-3, 0, 0, 0, 0, 5e-3, 0, 0, 0, 0])
    assert c.node_rot == {3: (30.0, 0.0, 0.0)}
    nums = [(e.num, e.type, e.real, e.nodes) for e in c.elements]
    assert nums[0] == (1, 2, 3, [1, 2, 3])
    assert nums[1] == (10, 2, 3, [2, 3, 1])
    assert nums[2] == (11, 1, 3, [1, 2, 3, 3, 1, 2, 3, 3, 9, 9])       # EMORE appended
    assert nums[3][0] == 12                                            # E numbers after the largest
    assert (1, "ALL", 0.0, 0.0) in c.constraints and (3, "ALL", 0.0, 0.0) in c.constraints   # NEND/NINC
    assert {(2, l) for l in ("UX", "UY", "ROTZ")} <= {(n, l) for n, l, _, _ in c.constraints}
    assert c.skipped == {"CE": 1}
    assert c.keyopt(2, 7) == 1 and c.unknown == {}                      # KEYO = KEYOPT (4-character abbreviation)


def test_reader_workbench_style_lowercase_and_short_nblock():
    """Workbench input (ds.dat style): lower case, abbreviated commands, NBLOCK with one integer
    field and a ``-1`` terminator, EBLOCK solid key in lower case."""
    text = """/prep7
et,1,185
keyo,1,2,3
nblock,3,,4
(1i9,3e20.9e3)
        1     0.000000000E+00     0.000000000E+00     0.000000000E+00
        2     1.000000000E+00     0.000000000E+00     0.000000000E+00
        3     1.000000000E+00     1.000000000E+00     0.000000000E+00
        4     0.000000000E+00     1.000000000E+00    -2.500000000E-01
-1
eblock,19,solid,,1
(19i9)
        1        1        1        1        0        0        0        0        4        0        7        1        2        3        4
-1
mp,ex,1,2e11
secd,0.1
"""
    c = read_cdb_text(text)
    assert c.etypes == {1: 185} and c.keyopt(1, 2) == 3
    np.testing.assert_array_equal(c.nodes[4], [0.0, 1.0, -0.25])
    assert [(e.num, e.nodes) for e in c.elements] == [(7, [1, 2, 3, 4])]
    assert c.materials == {1: {"EX": 2e11}} and not c.unknown


def test_element_names():
    assert element_name(185) == "SOLID185" and element_name(14) == "COMBIN14" and element_name(154) == "SURF154"


# ======================================================================================
# Beam axes (D-ANS-03)
# ======================================================================================
def test_ansys_default_orientation():
    A, L = ansys_beam_axes([0, 0, 0], [2, 0, 0])
    np.testing.assert_allclose(A, np.eye(3), atol=1e-15)                   # x = X, y = Y, z = Z
    assert L == 2.0
    A, _ = ansys_beam_axes([0, 0, 0], [0, 3, 0])                           # along Y: y = -X, z = Z
    np.testing.assert_allclose(A, [[0, 1, 0], [-1, 0, 0], [0, 0, 1]], atol=1e-15)
    A, _ = ansys_beam_axes([0, 0, 0], [0, 0, 5])                           # vertical: y = global Y
    np.testing.assert_allclose(A, [[0, 0, 1], [0, 1, 0], [-1, 0, 0]], atol=1e-15)
    A, _ = ansys_beam_axes([0, 0, 0], [1e-5, 0, 1])                        # within 0.01 % slope: vertical rule
    np.testing.assert_allclose(A[1], [0, 1, 0], atol=1e-12)
    A, _ = ansys_beam_axes([0, 0, 0], [2, 0, 0], theta_deg=90.0)           # THETA turns y towards z
    np.testing.assert_allclose(A, [[1, 0, 0], [0, 0, 1], [0, -1, 0]], atol=1e-15)


def test_ansys_k_node_orientation():
    A, _ = ansys_beam_axes([0, 0, 0], [4, 0, 0], [-3, 2, 0])               # K in the x-z plane: z towards K
    np.testing.assert_allclose(A, [[1, 0, 0], [0, 0, -1], [0, 1, 0]], atol=1e-15)
    with pytest.raises(ValueError):
        ansys_beam_axes([0, 0, 0], [4, 0, 0], [8, 0, 0])
    with pytest.raises(ValueError):
        ansys_beam_axes([1, 1, 1], [1, 1, 1])


def test_rect_section_maps_to_sassi_rectangle():
    """RECT B (y) x H (z): SASSI axis 2 = z, so I2 = Izz = H B^3/12, I3 = Iyy = B H^3/12."""
    from sassi.model.materials import section_rectangle
    p = rect_props(0.2, 0.5)
    assert p.Iyy == pytest.approx(0.2 * 0.5 ** 3 / 12) and p.Izz == pytest.approx(0.5 * 0.2 ** 3 / 12)
    sec, p2 = sassi_section_from_ansys(p)
    assert sec == section_rectangle(0.2, 0.5) and tuple(p2) == (0.0, 1.0)


def test_principal_axes_of_an_asec_section():
    p = BeamProps(A=1.0, Iyy=3.0, Izz=2.0, J=1.0, Iyz=0.5)
    sec, p2 = sassi_section_from_ansys(p)
    C = np.array([[2.0, 0.5], [0.5, 3.0]])                                 # [[Izz, Iyz], [Iyz, Iyy]]
    w = np.linalg.eigvalsh(C)
    assert sorted([sec["flex2"], sec["flex3"]]) == pytest.approx(sorted(w))
    assert sec["flex3"] == pytest.approx(p2 @ C @ p2)                      # I3 = int (r.e2)^2
    assert abs(p2[1]) > abs(p2[0]) and p2[1] > 0                           # nearest to +z


# ======================================================================================
# Conversion of the samples
# ======================================================================================
def test_cantilever_k_node_kept_and_materials():
    r = _conv("beam188_cantilever.cdb")
    m = r.model
    assert r.created_nodes == [] and not r.warnings
    g = m.groups[1]
    assert g.type == 2 and len(g.elements) == 20 and g.elements[7].nodes == [7, 8, 22]
    mt = m.materials[1]
    assert (mt.val1, mt.val2, mt.mtype) == (2.0e11, 0.3, 1)
    assert mt.weight == pytest.approx(7850.0 * 9.81)                       # weight = DENS g
    assert mt.pdamp == mt.sdamp == pytest.approx(0.02)                     # beta = DMPR/2 (D-ANS-01)
    assert m.gravity == 9.81                                               # model GRAVITY set
    assert m.nodes[1].fix == [1] * 6
    assert r.element_map[:2] == [(1, 1, 1), (2, 1, 2)]


def test_default_orientation_creates_one_shared_k_node():
    r = _conv("beam188_default_orientation.cdb")
    assert r.created_nodes == [22]
    np.testing.assert_allclose(r.model.node_global(22), [0, 0, 1])          # I + L e2, e2 = ANSYS z = +Z
    assert all(e.nodes[2] == 22 for e in r.model.groups[1].elements.values())
    assert _has(r, "no DMPR; damping ratio 0.05")


def test_rect_axis_mapping_static_closed_form():
    """D-ANS-03 unit test with B != H: tip deflections along Y and Z use Izz and Iyy (ANSYS axes)."""
    m = _conv("beam188_cantilever.cdb").model
    E, nu, L, B, H = 2e11, 0.3, 20.0, 0.2, 0.4
    G, A = E / (2 * (1 + nu)), B * H
    for d, I in ((2, H * B ** 3 / 12), (3, B * H ** 3 / 12)):
        ref = L ** 3 / (3 * E * I) + L / (G * A / 1.2)
        assert tip_deflection(m, 21, d) == pytest.approx(ref, rel=1e-9)


def _beam44_apdl(k7: int, k_node=(0.0, 0.0, 1.0)):
    return f"""
ET,1,BEAM44
KEYOPT,1,7,{k7}
R,1,0.02,2e-4,5e-4,0,0,1e-4
RMORE,0.02,2e-4,5e-4,0,0,1e-4
N,1,0,0,0
N,2,3,0,0
N,3,{k_node[0]},{k_node[1]},{k_node[2]}
MP,EX,1,1e7
MP,PRXY,1,0.25
MP,DENS,1,1
EN,1,1,2,3
D,1,ALL,0
D,2,ROTX,0,,,,ROTY,ROTZ
"""


@pytest.mark.parametrize("k_node", [(0.0, 0.0, 1.0), (0.0, -2.0, 0.0)])
def test_beam44_release_follows_the_ansys_axes(k_node):
    """KEYOPT(7) = 10 releases ROTY at I: bending about the element y axis (governed by IY, deflection
    along element z) becomes pinned-guided 3EI/L^3; bending about z stays fixed-guided 12EI/L^3."""
    r = convert_cdb(read_cdb_text(_beam44_apdl(10, k_node)), 9.81, 0.0)
    m = r.model
    assert m.groups[1].elements[1].ki == [0, 0, 0, 0, 0, 1]             # ROTY -> M3 (axis 3 = -y)
    E, L, IZ, IY = 1e7, 3.0, 2e-4, 5e-4
    ax, _ = ansys_beam_axes([0, 0, 0], [3, 0, 0], k_node)
    ez = ax[2]
    dz = int(np.argmax(np.abs(ez))) + 1                                    # global direction of element z
    dy = int(np.argmax(np.abs(ax[1]))) + 1
    assert 1.0 / tip_deflection(m, 2, dz) == pytest.approx(3 * E * IY / L ** 3, rel=1e-9)
    assert 1.0 / tip_deflection(m, 2, dy) == pytest.approx(12 * E * IZ / L ** 3, rel=1e-9)


def test_frame_sample_beams_links_pipes():
    r = _conv("frame_beams_links.cdb")
    m = r.model
    # BEAM4 column with THETA = 90 and no K: K node created, R from real constants
    col = m.groups[1].elements[1]
    s1 = m.sections[col.prop]
    assert (s1.axial, s1.flex2, s1.flex3, s1.tors) == pytest.approx((0.12, 1.6e-3, 9e-4, 2e-3))
    assert s1.shear2 == pytest.approx(0.1) and s1.shear3 == pytest.approx(0.1)      # A / SHEARZ, A / SHEARY
    from sassi.elements.base import frame_from_three_points
    lam, _ = frame_from_three_points(m.node_global(1), m.node_global(2), m.node_global(col.nodes[2]))
    np.testing.assert_allclose(lam[1], [0, -1, 0], atol=1e-14)              # vertical default y=Y, THETA 90 -> z=-Y
    # BEAM44: K node kept, KEYOPT(7) = 11 releases ROTY ROTZ at I -> KI M2 M3
    gird = m.groups[2].elements[1]
    assert gird.nodes == [2, 5, 6] and gird.ki == [0, 0, 0, 0, 1, 1] and gird.kj == [0] * 6
    assert _has(r, "KEYOPT(7) = 11", "info")
    # LINK180 / LINK8: axial only, rotations of truss-only nodes are not stiffened
    for gid, A in ((3, 5e-3), (4, 1e-2)):
        s = m.sections[m.groups[gid].elements[1].prop]
        assert (s.axial, s.flex2, s.flex3, s.shear2) == pytest.approx((A, 0, 0, 0)) and 0 < s.tors < 1e-9
    # BEAM189 I-J with midside K and orientation L -> two BEAMS I-K, K-J with K node 9
    g5 = m.groups[5]
    assert [e.nodes for e in g5.sorted_elements()] == [[3, 8, 9], [8, 7, 9]]
    assert [a for a, g, e in r.element_map if g == 5] == [7, 7]
    s5 = m.sections[g5.elements[1].prop]
    assert (s5.flex2, s5.flex3, s5.tors, s5.shear2) == pytest.approx((1e-3, 3e-3, 2e-3, 0.0))   # ASEC: Izz, Iyy
    # PIPE288 Do = 0.3, tw = 0.02
    sp = m.sections[m.groups[6].elements[1].prop]
    Do, Di = 0.3, 0.26
    assert sp.axial == pytest.approx(math.pi * (Do ** 2 - Di ** 2) / 4)
    assert sp.flex2 == sp.flex3 == pytest.approx(math.pi * (Do ** 4 - Di ** 4) / 64)
    assert sp.tors == pytest.approx(2 * sp.flex2)
    # the model is a valid, stable structure for the HOUSE element library
    f = fixed_base_frequencies(m)
    assert np.all(np.isfinite(f)) and f[0] > 0


def test_solids_sample():
    r = _conv("solids.cdb")
    m = r.model
    assert m.groups[1].elements[3].nodes == [9, 12, 13, 13, 10, 11, 14, 14]          # prism keeps repeats
    assert m.groups[3].elements[1].nodes == [15, 19, 20, 16, 17, 21, 22, 18]          # SOLID186 corners
    assert m.groups[4].elements[1].nodes == [35, 37, 36, 36, 38, 38, 38, 38]          # SOLID187 tetra
    assert all(e.etype == 0 for _, e in m.iter_elements())                          # ETYPE 0 (implicit)
    assert sorted(n for n, nd in m.nodes.items() if 0 in nd.flags) == [1, 2, 3, 4, 15, 16]   # SSI_INT
    assert m.nodes[19].fix == [0, 0, 1, 0, 0, 0]
    assert m.materials[2].pdamp == 0.05                                              # no DMPR -> default
    for frag in ("nodal coordinate rotations", "SOLID186: quadratic", "SOLID187: quadratic", "CP", "no DMPR"):
        assert _has(r, frag), frag
    assert _has(r, "BASE (node, 4)", "info")
    # every converted solid has a positive Jacobian (assembly would raise otherwise)
    dm, am = assemble_model(m)
    assert am.Ks.shape[0] == dm.neq > 0


def test_shells_sample():
    r = _conv("shells.cdb")
    m = r.model
    g1 = m.groups[1]
    assert g1.elements[4].nodes == [5, 6, 9]                                         # triangle: 3 nodes
    assert all(e.thick == 0.5 and e.mat == 1 for e in g1.elements.values())
    wall = m.groups[2].elements[1]
    assert (wall.thick, wall.mat) == (0.25, 2)                                       # SHELL63 R1
    s281 = m.groups[3].elements[1]
    assert s281.nodes == [10, 14, 15, 11] and s281.thick == pytest.approx(0.4)       # 2 layers summed
    assert _has(r, "2 layers") and _has(r, "SHELL281: 8-node")
    assert m.materials[1].sdamp == pytest.approx(0.04)
    dm, am = assemble_model(m)
    assert dm.neq > 0


def test_springs_masses_matrices_sample():
    r = _conv("springs_masses_matrices.cdb")
    m = r.model
    sc = {g: m.springs[m.groups[g].elements[1].prop].k for g in (1, 2, 3)}
    assert sc[1] == (0, 0, 5e4, 0, 0, 0)                     # KEYOPT(2) = 3: UZ
    assert sc[2] == (2e3, 0, 0, 0, 0, 0)                     # longitudinal along X
    assert sc[3] == (0, 0, 0, 0, 0, 700.0)                   # torsional about Z (KEYOPT(3) = 1)
    # oblique spring (4 -> 5 along (3, 4, 0)/5) -> GENERAL k e e^T (D-ANS-04)
    gobl = [g for g in m.groups.values() if g.type == 9 and "oblique" in g.title][0]
    K = m.matrices[gobl.elements[1].prop].full("R")
    e = np.array([0.6, 0.8, 0.0])
    np.testing.assert_allclose(K[:3, :3], 2e3 * np.outer(e, e), atol=1e-12)
    np.testing.assert_allclose(K[:3, 6:9], -2e3 * np.outer(e, e), atol=1e-12)
    # MASS21: two elements at node 3 summed, point mass at node 4, MUNITS 0
    assert m.tmass == {3: [3.0, 3.0, 3.0], 4: [3.0, 3.0, 3.0]}
    assert m.rmass == {3: [0.5, 1.0, 1.5]} and m.mass_units == {3: 0, 4: 0}
    # MATRIX27: C1..C78 upper triangle row by row
    g6, g7 = m.groups[6], m.groups[7]
    KR = m.matrices[g6.elements[1].prop].full("R")
    assert KR[0, 0] == KR[6, 6] == 1000.0 and KR[0, 6] == KR[6, 0] == -1000.0
    assert KR[5, 5] == KR[11, 11] == 50.0 and KR[5, 11] == -50.0
    MM = m.matrices[g7.elements[1].prop].full("M")
    np.testing.assert_array_equal(np.diag(MM), [2, 2, 2, 0, 0, 0] * 2)
    assert g6.elements[1].nodes == [7, 8]


def test_plane_2d_sample_is_rotated_into_xz():
    r = _conv("plane_2d.cdb")
    m = r.model
    np.testing.assert_allclose(m.node_global(5), [0.0, 0.0, 0.0])            # ANSYS (0, 0) -> (0, 0, 0)
    np.testing.assert_allclose(m.node_global(1), [0.0, 0.0, -1.0])           # ANSYS (0, -1) -> z = -1
    assert m.nodes[1].fix == [1, 0, 1, 0, 0, 0]                               # UX, UY -> UX, UZ
    assert m.nodes[4].fix == [0, 0, 1, 0, 0, 0]
    assert m.groups[2].elements[1].nodes == [3, 4, 8]                         # PLANE42 triangle
    assert m.options.record("HOUSE").integer(4) == 1                          # HOUSE <dim> = 1
    dm, am = assemble_model(m)
    assert not any("clockwise" in t for t in am.messages)                     # counter-clockwise in X-Z


def test_unsupported_items_are_all_reported():
    r = _conv("unsupported.cdb")
    for frag in ("SURF154: unsupported element type", "section 7 type BEAM I is not supported",
                 "KEYOPT(2) = 8", "offset TOP", "non-zero value", "label TEMP", "material 3",
                 "constraint equations (CE)", "nodal forces (F)", "(TB)", "ALPHAD = 0.5"):
        assert _has(r, frag), frag
    assert list(r.model.groups) == [4]                                        # no empty groups
    assert r.counts["unsupported elements"] == 1


def test_material_rules():
    text = """ET,1,SOLID185
KEYOPT,1,3,1
MP,EX,1,100
MP,GXY,1,40
MP,DENS,1,2
MP,DAMP,1,0.01
MP,EX,2,100
MP,EY,2,50
MP,DMPS,2,0.06
MPTEMP,R5.0, 2, 1,  0.00000000    , 100.0
MPDATA,R5.0, 2,EX  ,       3, 1,  200.0, 180.0,
N,1,0,0,0
N,2,1,0,0
N,3,1,1,0
N,4,0,1,0
N,5,0,0,1
N,6,1,0,1
N,7,1,1,1
N,8,0,1,1
MAT,1
EN,1,1,2,3,4,5,6,7,8
MAT,2
EN,2,1,2,3,4,5,6,7,8
MAT,3
EN,3,1,2,3,4,5,6,7,8
"""
    r = convert_cdb(read_cdb_text(text), 32.2, None)
    m = r.model
    assert m.materials[1].val2 == pytest.approx(100 / 80 - 1)            # nu from EX and GXY
    assert m.materials[1].weight == pytest.approx(2 * 32.2)
    assert m.materials[2].sdamp == pytest.approx(0.03)                   # DMPS / 2
    assert m.materials[2].val2 == 0.3                                     # ANSYS default PRXY
    assert m.materials[3].val1 == 200.0                                   # first temperature
    assert m.materials[1].sdamp == 0.0                                    # no DMPR, no default -> 0
    for frag in ("DAMP = 0.01", "orthotropic", "temperature-dependent", "default 0.3", "KEYOPT(3) = 1 (layered",
                 "DENS missing"):
        assert _has(r, frag), frag


def test_gravity_must_be_positive():
    with pytest.raises(ValueError):
        convert_cdb(read_cdb(DATA / "beam188_cantilever.cdb"), 0.0)


def test_map_text():
    r = _conv("beam188_cantilever.cdb")
    lines = r.map_text("x.cdb").splitlines()
    assert lines[0].startswith("!") and lines[2].split() == ["1", "1", "1"] and len(lines) == 22
