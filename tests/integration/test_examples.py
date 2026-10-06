"""Integration tests: the tutorial examples of ``examples/`` run through the command language.

Each example is staged in a temporary copy (``sassi.verify.problems.vp_examples.stage_example``) and
executed by an :class:`~sassi.prep.Interpreter` whose working directory is the staged ``examples``
directory -- exactly what ``sassi --cwd examples run exNN.pre`` does (examples/README.md).  The
interpreter builds the model, CHECK and AFWRITE write the decks, the RUN<MODULE> commands run the
modules.  The tests assert that every module finishes with status OK and that the key results are
physically sensible (requirements 2.4, 3.3; ARCHITECTURE sections 3 and 8):

* low-frequency transfer functions -> 1 (everything moves with the ground), rocking antisymmetric
  on a symmetric model, the ISRS zero-period acceleration of the base close to the input PGA, the
  relative displacement of the base with respect to itself zero, storey shears growing downwards and
  the top-storey shear equal to top mass x top acceleration (ex01);
* FV and FI-EVBN alike, FI-FSIN with its spurious roof resonance (ex02);
* inverse compliance of six load cases = global impedance K_G, symmetric and close to Pais-Kausel (ex03);
* EQUAKE acceptance, strain-compatible soil softer than the low-strain soil, SITE with FILE88
  reproducing the SOIL amplification (ex04);
* simultaneous X/Y/Z cases with their control directions, the X-to-Y coupling of the eccentric mass
  kept under its own file names, RELDISP with the free-field reference in Y (ex05);
* an embedded shear-wall building: clean CHECK, INTGEN counts and masses, FI-EVBN close to FV, FI-FSIN
  with the spurious resonance of the enclosed excavated soil near 15.5 Hz, its ISRS and wall forces (ex08);
* WRITE -> INP gives back the same model (UT-03) for every model of every example.

The binary inter-module files are deleted as soon as each example has run and the results the tests
need are in memory (``purge_binaries``), to keep the disk footprint of the test session small.
"""
from __future__ import annotations

import numpy as np
import pytest

from sassi.io import decks, textfiles
from sassi.io.container import read_container
from sassi.prep import Interpreter, Kind
from sassi.prep.check import run_check
from sassi.verify import builders as B
from sassi.verify.problems.vp_examples import (EXAMPLES, EXAMPLES_DIR, listing_status, purge_binaries,
                                               run_example, stage_example)

pytestmark = pytest.mark.skipif(not (EXAMPLES_DIR / "ex01_surface_stick.pre").exists(),
                                reason="the tutorial examples ship with the source tree only")

G = 9.81


# ======================================================================================
# helpers
# ======================================================================================
def _run(name, tmp_path_factory):
    """Run example ``name`` in a fresh temporary root; returns (interpreter, summary, examples dir)."""
    root = tmp_path_factory.mktemp(name)
    ui, summary, pre = run_example(name, root)
    return ui, summary, pre.parent


def _assert_clean_run(ui, summary, ex_dir, name, modules):
    errors = ui.sink.texts(Kind.ERROR)
    assert summary.ok and summary.errors == 0 and not errors, "\n".join(errors[:20])
    for d, mods in modules.items():
        status = listing_status(ex_dir / d)
        model = {"ex02_fsin": "ex02fsin", "ex02_evbn": "ex02evbn", "ex08_fsin": "ex08fsin",
                 "ex08_evbn": "ex08evbn"}.get(d, d)
        for mod in mods:
            line = status.get(f"{model}_{mod}.out", "listing missing")
            assert "status OK" in line, f"{name} {d} {mod}: {line}"


def _assert_write_inp_roundtrip(ui, ex_dir, written):
    """``written`` = {model number: .pre written by the example's WRITE}; INP of each gives the same model."""
    for num, rel in written.items():
        path = ex_dir / rel
        assert path.exists(), path
        b = Interpreter(cwd=path.parent)
        b.run_file(str(path), resolve=False)
        assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)[:10]
        assert b.model.same_state(ui.models[num]), f"WRITE -> INP of {rel} differs from model {num}"


def _peak(f, h):
    k = int(np.argmax(np.abs(h)))
    return float(f[k]), float(np.abs(h[k]))


def _rs(path):
    """(frequencies, SA) of a .RS file (first damping column)."""
    t = textfiles.read_xy(path)
    return t[:, 0], t[:, 1]


def _pga(path):
    acc, _ = textfiles.read_history(path)
    return float(np.max(np.abs(acc)))


# ======================================================================================
# Example 1: stick on a rigid surface mat
# ======================================================================================
@pytest.fixture(scope="module")
def ex01(tmp_path_factory):
    ui, summary, ex = _run("ex01_surface_stick", tmp_path_factory)
    md = ex / "ex01"
    f8 = B.read_file8(md) if (md / "FILE8").exists() else None
    purge_binaries(md)                                       # results needed are in memory or text files
    return ui, summary, ex, md, f8


def test_ex01_runs_every_module(ex01):
    ui, summary, ex, md, f8 = ex01
    _assert_clean_run(ui, summary, ex, "ex01", {"ex01": ("SITE", "POINT", "HOUSE", "ANALYS", "MOTION", "STRESS",
                                                         "RELDISP")})
    assert f8 is not None
    _assert_write_inp_roundtrip(ui, ex, {0: "ex01/ex01.pre"})


def test_ex01_low_frequency_transfer_functions_are_rigid_body(ex01):
    *_, f8 = ex01
    assert float(f8["freq"][0]) < 0.11
    for node in (41, 82, 83, 84, 85):
        assert abs(abs(B.tf(f8, node, 1)[0]) - 1.0) < 0.01, node
    for node in (37, 45):                                   # no rocking at 0.1 Hz
        assert abs(B.tf(f8, node, 3)[0]) < 0.01


def test_ex01_rocking_is_antisymmetric(ex01):
    *_, f8 = ex01
    h37, h45 = B.tf(f8, 37, 3), B.tf(f8, 45, 3)
    assert np.max(np.abs(h37)) > 0.01                       # the mat rocks under the SV wave
    assert np.max(np.abs(h37 + h45)) <= 1e-8 * np.max(np.abs(h37))


def test_ex01_ssi_frequency_below_fixed_base(ex01):
    _, _, _, md, _ = ex01
    f, h, _ = textfiles.read_tf(md / "00085TR_X.TFI")
    fp, hp = _peak(f, h)
    assert 3.0 < fp < 4.0, fp                               # fixed base 5.05 Hz: SSI lowers it to ~3.5 Hz
    assert 3.0 < hp < 15.0, hp                              # amplified, with radiation damping


def test_ex01_isrs(ex01):
    _, _, ex, md, _ = ex01
    pga = _pga(ex / "data" / "rg160h_030g.acc")
    zpa = {n: _rs(md / f"0{n:04d}TR_X01.RS")[1][-1] for n in (41, 82, 83, 84, 85)}
    assert 0.8 < zpa[41] / pga < 1.3, (zpa[41], pga)       # base ISRS ZPA ~ input PGA (inertial SSI)
    z = [zpa[n] for n in (41, 82, 83, 84, 85)]
    assert all(b > a for a, b in zip(z, z[1:])), z          # amplification up the stick
    f, sa = _rs(md / "00085TR_X01.RS")                      # 2 % damping
    assert 3.0 < f[int(np.argmax(sa))] < 4.0
    for n in (41, 85):                                      # ZPA = peak of the saved floor history
        assert abs(zpa[n] / _pga(md / f"0{n:04d}TR_X.ACC") - 1.0) < 0.02


def test_ex01_relative_displacements(ex01):
    _, _, _, md, _ = ex01
    d = {n: np.max(np.abs(textfiles.read_history(md / f"0{n:04d}TR_X.THD")[0])) for n in (41, 82, 83, 84, 85)}
    assert d[41] <= 1e-12 * d[85]                           # the base relative to itself
    v = [d[n] for n in (82, 83, 84, 85)]
    assert all(b > a for a, b in zip(v, v[1:])) and 0.005 < d[85] < 0.1, d


def test_ex01_beam_forces(ex01):
    _, _, _, md, _ = ex01
    shear = [np.max(np.abs(textfiles.read_history(md / f"BEAMS_002_{e:05d}_FYI.THS")[0])) for e in (1, 2, 3, 4)]
    assert all(a > b for a, b in zip(shear, shear[1:])), shear      # storey shear grows downwards
    a_top = _pga(md / "00085TR_X.ACC")
    assert abs(shear[3] / (1000.0 * G * a_top) - 1.0) < 0.02        # top storey: V = m a (massless stick)
    moment = [np.max(np.abs(textfiles.read_history(md / f"BEAMS_002_{e:05d}_MZI.THS")[0])) for e in (1, 2, 3, 4)]
    assert all(a > b for a, b in zip(moment, moment[1:])), moment


# ======================================================================================
# Example 2: embedded box, FV / FI-FSIN / FI-EVBN
# ======================================================================================
@pytest.fixture(scope="module")
def ex02(tmp_path_factory):
    ui, summary, ex = _run("ex02_embedded_box", tmp_path_factory)
    f8 = {d: B.read_file8(ex / d) for d in EXAMPLES["ex02_embedded_box"] if (ex / d / "FILE8").exists()}
    for d in EXAMPLES["ex02_embedded_box"]:
        purge_binaries(ex / d)
    return ui, summary, ex, f8


def test_ex02_runs_every_model(ex02):
    ui, summary, ex, f8 = ex02
    full = ("HOUSE", "ANALYS", "MOTION", "STRESS")
    _assert_clean_run(ui, summary, ex, "ex02", {"ex02": ("SITE", "POINT") + full, "ex02_fsin": full,
                                                "ex02_evbn": full})
    assert set(f8) == {"ex02", "ex02_fsin", "ex02_evbn"}
    _assert_write_inp_roundtrip(ui, ex, {0: "ex02/ex02.pre", 2: "ex02_fsin/ex02fsin.pre",
                                         3: "ex02_evbn/ex02evbn.pre"})


def test_ex02_embedment_layers_and_excavated_soil_are_consistent(ex02):
    """Manual rules: embedment layers are not repeated in TOPL, and each excavated group uses the L
    number of its layer (CHECK EDU-28 / EDU-08 silent)."""
    ui = ex02[0]
    for num in (0, 2, 3):
        rep = run_check(ui.models[num], dirs=[ui.models[num].path])
        bad = [m.line() for m in rep.messages if m.number in ("EDU-08", "EDU-28")]
        assert not bad, bad


def test_ex02_low_frequency_and_antisymmetry(ex02):
    f8 = ex02[3]
    for d, f in f8.items():
        for node in (13, 138):
            assert abs(abs(B.tf(f, node, 1)[0]) - 1.0) < 0.01, (d, node)
        h136, h140 = B.tf(f, 136, 3), B.tf(f, 140, 3)
        assert np.max(np.abs(h136 + h140)) <= 1e-8 * np.max(np.abs(h136)), d


def test_ex02_interaction_sets(ex02):
    """FI-EVBN follows FV; FI-FSIN has the spurious roof resonance (roof slab on non-interaction nodes)."""
    f8 = ex02[3]
    fv = f8["ex02"]

    def dev(d, node, dof):
        a, b = B.tf(f8[d], node, dof), B.tf(fv, node, dof)
        return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))

    for node, dof in ((13, 1), (138, 1), (136, 3)):
        assert dev("ex02_evbn", node, dof) < 0.03, (node, dof)
    assert dev("ex02_fsin", 138, 1) > 0.2
    freq = np.asarray(f8["ex02_fsin"]["freq"])
    fp, hp = _peak(freq, B.tf(f8["ex02_fsin"], 138, 1))
    assert 5.0 < fp < 7.5 and hp > 1.2, (fp, hp)
    assert np.max(np.abs(B.tf(fv, 138, 1))) < 1.05          # FV: no resonance of the roof


# ======================================================================================
# Example 3: forced vibration, compliance and global impedance
# ======================================================================================
@pytest.fixture(scope="module")
def ex03(tmp_path_factory):
    ui, summary, ex = _run("ex03_forced_vibration", tmp_path_factory)
    md = ex / "ex03"
    C = None
    for k in range(6):
        f8 = B.read_file8(md, f"FILE800{k + 1}")
        col = np.column_stack([B.tf(f8, 25, d) for d in range(1, 7)])
        if C is None:
            C = np.zeros((col.shape[0], 6, 6), complex)
            freq = np.asarray(f8["freq"])
        C[:, :, k] = col
    KG = np.asarray(read_container(md / "FILE11")["KG"])
    purge_binaries(md)
    return ui, summary, ex, md, freq, C, KG


def test_ex03_runs_every_module(ex03):
    ui, summary, ex, *_ = ex03
    _assert_clean_run(ui, summary, ex, "ex03", {"ex03": ("SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "MOTION")})
    for f in ("FOUNSTIF", "FOUNDASH", "FOUNDAMP", "FOUNIMPD"):
        assert (ex / "ex03" / f).exists(), f
    _assert_write_inp_roundtrip(ui, ex, {0: "ex03/ex03.pre"})


def test_ex03_inverse_compliance_equals_global_impedance(ex03):
    *_, freq, C, KG = ex03
    K = np.linalg.inv(C)
    err = np.max(np.abs(K - KG), axis=(1, 2)) / np.max(np.abs(KG), axis=(1, 2))
    assert np.max(err) < 0.005, err                          # rigid links: stiff, not infinitely rigid
    assert np.max(np.abs(K - np.transpose(K, (0, 2, 1)))) <= 1e-8 * np.max(np.abs(K))
    assert np.max(np.abs(K[:, 0, 0] - K[:, 1, 1])) <= 1e-8 * np.max(np.abs(K[:, 0, 0]))     # square: Kx = Ky
    assert np.max(np.abs(K[:, 3, 3] - K[:, 4, 4])) <= 1e-8 * np.max(np.abs(K[:, 3, 3]))     # Kxx = Kyy


def test_ex03_static_stiffness_close_to_pais_kausel(ex03):
    """Coarse 2 m mesh: stiffer than Pais-Kausel (R2 B.2), about +6 % for the forces and +18..20 % for
    the moments, as the tutorial text states."""
    *_, freq, C, KG = ex03
    K0 = np.real(np.linalg.inv(C[0]))
    assert freq[0] < 0.25
    Gs, Bh, nu = 80000.0, 6.0, 1.0 / 3.0
    ref = {0: Gs * Bh / (2 - nu) * 9.2, 2: Gs * Bh / (1 - nu) * 4.7, 3: Gs * Bh ** 3 / (1 - nu) * 4.0,
           5: Gs * Bh ** 3 * 8.31}
    over = {i: K0[i, i] / v - 1.0 for i, v in ref.items()}
    assert 0.0 < over[0] < 0.10 and 0.0 < over[2] < 0.10, over
    assert 0.10 < over[3] < 0.25 and 0.10 < over[5] < 0.25, over


def test_ex03_ricker_pulse_response(ex03):
    _, _, _, md, freq, C, KG = ex03
    u, dt = textfiles.read_history(md / "00025TR_Z.ACC")             # MOTIONX <resp> = 0: displacement
    static = 1000.0 / np.real(np.linalg.inv(C[0]))[2, 2]
    k = int(np.argmax(np.abs(u)))
    assert 0.5 < np.max(np.abs(u)) / static < 1.0                    # |K(w)| grows with frequency
    assert 0.4 < k * dt < 0.7                                          # peak near the pulse peak (0.5 s)
    assert np.max(np.abs(u[int(1.0 / dt):])) < 0.01 * np.max(np.abs(u))   # radiation damping


# ======================================================================================
# Example 4: EQUAKE -> SOIL -> SITEX 1 -> SITE
# ======================================================================================
@pytest.fixture(scope="module")
def ex04(tmp_path_factory):
    ui, summary, ex = _run("ex04_site_response", tmp_path_factory)
    md = ex / "ex04"
    f1 = read_container(md / "FILE1") if (md / "FILE1").exists() else None
    purge_binaries(md)
    return ui, summary, ex, md, f1


@pytest.mark.slow                     # EQUAKE spectral matching: about 15 s
def test_ex04_runs_the_free_field_chain(ex04):
    ui, summary, ex, md, f1 = ex04
    assert f1 is not None
    _assert_clean_run(ui, summary, ex, "ex04", {"ex04": ("EQUAKE", "SOIL", "SITE")})
    _assert_write_inp_roundtrip(ui, ex, {0: "ex04/ex04.pre"})


@pytest.mark.slow
def test_ex04_equake_meets_the_acceptance_criteria(ex04):
    md = ex04[3]
    txt = (md / "ex04_EQUAKE.out").read_text(encoding="utf-8")
    assert "PASS" in txt and "FAIL" not in txt
    assert 0.27 < _pga(md / "rg160h_eq.acc") < 0.40           # anchored to 0.30 g


@pytest.mark.slow
def test_ex04_strain_compatible_soil_is_softer(ex04):
    md = ex04[3]
    t = np.loadtxt(md / "FILE88", comments="#")
    assert t.shape[0] == 22
    vs0 = np.where(t[:, 0] <= 10, 200.0, 300.0)
    assert np.all(t[:, 4] < vs0) and np.all(t[:, 5] > 0.01)   # G reduced, damping increased
    assert t[:10, 4].min() < 0.5 * 200.0                      # strongly softened sand near 10 m


@pytest.mark.slow
def test_ex04_site_reproduces_the_soil_amplification(ex04):
    """SITE (thin layers, FILE88 properties) vs the final SHAKE amplification of SOIL (SSAF)."""
    md, f1 = ex04[3], ex04[4]
    U = np.asarray(f1["U"])[0]
    freq = np.asarray(f1["freq"])
    amp = U[:, 0, 0] / U[:, -1, 0]                             # surface / top of the half-space (within)
    f, H, cplx = textfiles.read_tf(md / "SAF001W_023W.TFU")
    assert cplx
    for fi, a in zip(freq, amp):
        k = int(np.argmin(np.abs(f - fi)))
        if abs(f[k] - fi) > 1e-3:
            continue                                           # 0.1 Hz is not on the SSAF grid
        err = abs(a - H[k]) / abs(H[k])
        assert err < (0.01 if fi < 7.5 else 0.06), (fi, err)  # discretisation error grows with (k h)^2
    fp, hp = _peak(freq, amp)
    assert 1.5 < fp < 2.5 and hp > 3.0                         # softened site: ~1.8 Hz (low strain: 3.1 Hz)


# ======================================================================================
# Example 5: X, Y, Z simultaneous cases
# ======================================================================================
@pytest.fixture(scope="module")
def ex05(tmp_path_factory):
    ui, summary, ex = _run("ex05_xyz_simultaneous", tmp_path_factory)
    md = ex / "ex05"
    f8 = {c: B.read_file8(md, f"FILE8{c}") for c in "XYZ" if (md / f"FILE8{c}").exists()}
    purge_binaries(md)
    return ui, summary, ex, md, f8


def test_ex05_runs_every_module(ex05):
    ui, summary, ex, md, f8 = ex05
    _assert_clean_run(ui, summary, ex, "ex05", {"ex05": ("SITE", "POINT", "HOUSE", "ANALYS", "MOTION", "RELDISP")})
    assert set(f8) == set("XYZ")
    assert [int(f8[c].meta["cm"]) for c in "XYZ"] == [0, 1, 2]
    _assert_write_inp_roundtrip(ui, ex, {0: "ex05/ex05.pre"})


def test_ex05_each_case_moves_with_its_control_motion(ex05):
    f8 = ex05[4]
    for c, dof in (("X", 1), ("Y", 2), ("Z", 3)):
        for node in (13, 26, 27, 30):
            assert abs(abs(B.tf(f8[c], node, dof)[0]) - 1.0) < 0.01, (c, node)


def test_ex05_eccentric_mass_couples_x_to_y(ex05):
    _, _, _, md, f8 = ex05
    x2y = B.tf(f8["X"], 30, 2)
    assert abs(x2y[0]) < 1e-3 and np.max(np.abs(x2y)) > 0.5          # coupling only by dynamics
    assert np.max(np.abs(B.tf(f8["X"], 27, 6))) > 0.05                 # torsion of the stick top
    f, h, _ = textfiles.read_tf(md / "X2Y_00030TR_Y.TFU")              # the X-run file, kept by FCOPY
    assert np.allclose(h, x2y, rtol=0, atol=1e-6 * np.max(np.abs(x2y)))
    f, h, _ = textfiles.read_tf(md / "00030TR_Y.TFU")                  # the Y run's own (direct) response
    assert np.allclose(h, B.tf(f8["Y"], 30, 2), rtol=0, atol=1e-6 * np.max(np.abs(h)))
    assert np.max(np.abs(h)) > 2.0 * np.max(np.abs(x2y))


def test_ex05_decks_carry_the_direction_of_their_case(ex05):
    md = ex05[3]
    assert decks.read(decks.deck_path(md, "ex05", "RELDISP"), "RELDISP")["cm"] == 1     # FILE8Y
    mot = decks.read(decks.deck_path(md, "ex05", "MOTION"), "MOTION")                    # last: FILE8Z
    assert mot["cm"] == 2 and mot["file8"] == "FILE8Z"
    d, _ = textfiles.read_history(md / "00027TR_Y.THD")
    assert 1e-4 < np.max(np.abs(d)) < 0.05                             # roof drift relative to the ground


# ======================================================================================
# Example 8: embedded shear-wall building, FV / FI-FSIN (SM) / FI-EVBN (MSM)
# ======================================================================================
EX08_DIRS = {"fv": "ex08", "fsin": "ex08_fsin", "evbn": "ex08_evbn"}
#: the structure outputs: basemat, basement slab, grade slab, floor, main roof corner and tower roof (X);
#: basemat and grade edges, the grade slab under a 150 t item and the two roof corners (Z)
EX08_X = (41, 675, 757, 446, 567, 605)
EX08_Z = (45, 369, 739, 567, 617)
EX08_SOIL = 365                                  # excavated soil, centre of the excavation top face


@pytest.fixture(scope="module")
def ex08(tmp_path_factory):
    ui, summary, ex = _run("ex08_embedded_building", tmp_path_factory)
    f8 = {k: B.read_file8(ex / d) for k, d in EX08_DIRS.items() if (ex / d / "FILE8").exists()}
    for d in EX08_DIRS.values():
        purge_binaries(ex / d)
    return ui, summary, ex, f8


def _ex08_dev(f8, method, node, dof):
    """max |H_method - H_FV| / max |H_FV| per frequency (complex difference, relative to the FV peak)."""
    a, b = B.tf(f8[method], node, dof), B.tf(f8["fv"], node, dof)
    return np.abs(a - b) / np.max(np.abs(b))


@pytest.mark.slow                     # three ANALYS runs of 41 frequencies: about 1 minute
def test_ex08_runs_every_model(ex08):
    ui, summary, ex, f8 = ex08
    full = ("HOUSE", "ANALYS", "MOTION", "STRESS")
    _assert_clean_run(ui, summary, ex, "ex08", {"ex08": ("SITE", "POINT") + full, "ex08_fsin": full,
                                                "ex08_evbn": full})
    assert set(f8) == set(EX08_DIRS)
    _assert_write_inp_roundtrip(ui, ex, {0: "ex08/ex08.pre", 2: "ex08_fsin/ex08fsin.pre",
                                         3: "ex08_evbn/ex08evbn.pre"})


@pytest.mark.slow
def test_ex08_model_checks_counts_and_masses(ex08):
    """CHECK: no error; no warning for FV, only the EDU-12 reminder for the reduced sets.  The interior
    structure has its own nodes (no excavation-interior node shared), the numbering has no gap, the
    INTGEN sets have the counts of the documentation and HOUSE has the masses of CALCM."""
    ui, _, ex, _ = ex08
    counts = {0: 405, 2: 209, 3: 258}
    for num, n_int in counts.items():
        m = ui.models[num]
        rep = run_check(m, dirs=[m.path])
        assert not rep.errors(), [x.line() for x in rep.errors()][:5]
        warn = {x.number for x in rep.warnings()}
        assert warn == (set() if num == 0 else {"EDU-12"}), warn
        assert sum(1 for nd in m.nodes.values() if 0 in nd.flags) == n_int
        assert len(m.nodes) == max(m.nodes) == 781
        assert m.ui_state.get("hide_groups") == [1, 2, 3, 4]       # the soil is hidden in the 3D views
    for d in EX08_DIRS.values():
        txt = (ex / d / f"{d.replace('_', '')}_HOUSE.out").read_text(encoding="utf-8")
        line = next(ln for ln in txt.splitlines() if "Total structural mass" in ln)
        assert abs(float(line.split(":")[1].split()[0]) - 16555.7) < 1.0, line     # 15,746 t + 810 t equipment
        line = next(ln for ln in txt.splitlines() if "Total excavated soil mass" in ln)
        assert abs(float(line.split(":")[1].split()[0]) - 8924.8) < 1.0, line      # 24 x 24 x 8 m x 1.937 t/m3


@pytest.mark.slow
def test_ex08_low_frequency_transfer_functions_are_rigid_body(ex08):
    f8 = ex08[3]
    for k, f in f8.items():
        assert float(f["freq"][0]) < 0.11
        for node in EX08_X + (EX08_SOIL,):
            assert abs(abs(B.tf(f, node, 1)[0]) - 1.0) < 0.01, (k, node)
        for node in EX08_Z:
            assert abs(B.tf(f, node, 3)[0]) < 0.01, (k, node)


@pytest.mark.slow
def test_ex08_modified_subtraction_follows_fv(ex08):
    """FI-EVBN: within 3 % of the FV peak in X and 5 % in Z at every computed frequency (observed 1.3 % and
    3.7 %, the grade slab under a 150 t item at 20 Hz)."""
    f8 = ex08[3]
    for node in EX08_X:
        assert np.max(_ex08_dev(f8, "evbn", node, 1)) < 0.03, node
    for node in EX08_Z:
        assert np.max(_ex08_dev(f8, "evbn", node, 3)) < 0.05, node


@pytest.mark.slow
def test_ex08_subtraction_method_has_the_spurious_resonance(ex08):
    """FI-FSIN follows FV below 10 Hz, then resonates near 15.5-16 Hz: the soil enclosed by the walls and the
    basemat (the non-interaction excavated nodes) has its first natural frequency there (15.5 Hz)."""
    f8 = ex08[3]
    freq = np.asarray(f8["fv"]["freq"])
    low, band = freq < 10.0, (freq > 14.5) & (freq < 17.5)
    for node, dof in [(n, 1) for n in EX08_X] + [(n, 3) for n in EX08_Z]:
        dev = _ex08_dev(f8, "fsin", node, dof)
        assert np.max(dev[low]) < 0.03, (node, dof)
        k = int(np.argmax(dev))
        assert 14.5 < freq[k] < 17.5 and dev[k] > 0.15, (node, dof, freq[k], dev[k])
    res = (freq > 15.2) & (freq < 16.3)                     # the computed points at 15.5 and 16 Hz
    for node in (41, 757, 605):                             # basemat, grade, tower roof: about 2 times FV
        ratio = np.abs(B.tf(f8["fsin"], node, 1)) / np.abs(B.tf(f8["fv"], node, 1))
        assert np.max(ratio[res]) > 1.7, (node, ratio[res])
    soil = {k: np.abs(B.tf(f8[k], EX08_SOIL, 1)) for k in f8}
    assert np.max(soil["fsin"][res] / soil["fv"][res]) > 5.0  # the enclosed soil rings (12.6 against 1.6)
    assert np.max(np.abs(B.tf(f8["fv"], 41, 1))) < 1.01     # FV: the basemat never exceeds the free surface
    fp, _ = _peak(freq, B.tf(f8["fv"], 605, 1))
    assert 6.0 < fp < 7.0, fp                               # SSI peak of the tower roof (fixed base 12.4 Hz)


@pytest.mark.slow
def test_ex08_isrs(ex08):
    """5 % ISRS: amplification up the building; MSM within 2 % of FV; SM within 6 % below 10 Hz (observed
    4.9 %) and off by more than 5 % in the 13-17 Hz band of its spurious resonance (observed up to 31 %)."""
    _, _, ex, _ = ex08
    pga = _pga(ex / "data" / "rg160h_030g.acc")
    zpa = [_rs(ex / "ex08" / f"0{n:04d}TR_X01.RS")[1][-1] for n in (41, 757, 605)]
    assert 0.7 < zpa[0] / pga < 1.0, (zpa, pga)             # kinematic interaction: basemat below the PGA
    assert zpa[0] < zpa[1] < zpa[2], zpa
    files = [f"0{n:04d}TR_X01.RS" for n in EX08_X] + [f"0{n:04d}TR_Z01.RS" for n in EX08_Z]
    for fn in files:
        f, fv = _rs(ex / "ex08" / fn)
        r_evbn = _rs(ex / "ex08_evbn" / fn)[1] / fv
        r_fsin = _rs(ex / "ex08_fsin" / fn)[1] / fv
        assert np.max(np.abs(r_evbn - 1.0)) < 0.02, fn
        assert np.max(np.abs(r_fsin[f < 10.0] - 1.0)) < 0.06, fn
        band = (f > 13.0) & (f < 17.0)
        assert np.max(np.abs(r_fsin[band] - 1.0)) > 0.05, fn


@pytest.mark.slow
def test_ex08_basement_wall_forces(ex08):
    """Peak in-plane shear of the south wall and bending of the west wall: SM and MSM within 3 % of FV
    (they are governed by the 6.5 Hz SSI response, far below the spurious resonance)."""
    _, _, ex, _ = ex08
    for e, comp in ((4, "FXY"), (5, "FXY"), (68, "MYY"), (69, "MYY")):
        peak = {k: np.max(np.abs(textfiles.read_history(ex / d / f"SHELL_012_{e:05d}_{comp}.THS")[0]))
                for k, d in EX08_DIRS.items()}
        for k in ("fsin", "evbn"):
            assert abs(peak[k] / peak["fv"] - 1.0) < 0.03, (e, comp, peak)


# ======================================================================================
# WRITE -> INP of every model of every example (no module runs: fast)
# ======================================================================================
@pytest.mark.parametrize("name", sorted(EXAMPLES))
def test_write_inp_roundtrip_of_the_example_models(name, tmp_path):
    pre = stage_example(name, tmp_path)
    keep = [ln for ln in pre.read_text(encoding="utf-8").splitlines()
            if not ln.strip().upper().startswith(("RUN", "FCOPY"))]
    pre.write_text("\n".join(keep) + "\n", encoding="utf-8")
    ui = Interpreter(cwd=pre.parent)
    ui.run_file(str(pre), resolve=False)
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)[:10]
    from sassi.prep.writer import write_pre
    for num, m in sorted(ui.models.items()):
        if not m.nodes and not m.layers:
            continue
        text, _ = write_pre(m)
        b = Interpreter(cwd=tmp_path)
        b.run_text(text)
        assert not b.sink.texts(Kind.ERROR), b.sink.texts(Kind.ERROR)[:10]
        assert b.model.same_state(m), f"{name}: model {num}"


# ======================================================================================
# The documented command line (examples/README.md)
# ======================================================================================
def test_documented_command_line_runs_from_the_repository_root(tmp_path, monkeypatch):
    """``sassi --cwd examples run ex01_surface_stick.pre`` from the repository root: the model directory
    and the control-motion file resolve inside ``examples`` (the decks get the absolute data path)."""
    from sassi.cli import main
    pre = stage_example("ex01_surface_stick", tmp_path)
    lines = pre.read_text(encoding="utf-8").splitlines()
    stop = next(i for i, ln in enumerate(lines) if ln.strip().upper() == "AFWRITE")
    pre.write_text("\n".join(lines[:stop + 1]) + "\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)                                    # the "repository root"
    assert main(["--cwd", "examples", "run", pre.name, "--quiet"]) == 0
    md = tmp_path / "examples" / "ex01"
    mot = decks.read(decks.deck_path(md, "ex01", "MOTION"), "MOTION")
    assert (md / mot["thfile"]).resolve() == (tmp_path / "examples" / "data" / "rg160h_030g.acc").resolve()
    assert not (tmp_path / "ex01").exists()                        # no stray model directory
