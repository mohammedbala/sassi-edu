"""Unit tests of the ANALYS module (sassi/modules/analys.py; requirements 2.1, 2.5, 2.6, 4.6;
D-ANL-01..13; spec 05c A.5/A.7)."""
from __future__ import annotations

import math
import shutil
from types import SimpleNamespace

import numpy as np
import pytest

import sassi.modules.analys as analys
from sassi.core import ssi_solver as ss
from sassi.core.flexibility import flexibility_matrix, frequency_row
from sassi.core.freefield import free_field_at_nodes
from sassi.io.container import read_container, write_container
from sassi.io.files import validate
from sassi.modules.base import Listing
from sassi.verify import builders as B

FS = B.FrequencySet.fourier(0.01, 1024, [4, 20, 60])


def _site() -> B.Site:
    return B.layered_site([(2.0, 150.0, 300.0, 1.9, 0.05), (3.0, 250.0, 500.0, 2.0, 0.04)],
                          (400.0, 800.0, 2.1, 0.03))


def _model() -> B.Model:
    stick = B.Stick(heights=[3.0, 6.0], masses=[40.0, 30.0], E=3.0e7, A=1.0, I=0.5, beta=0.05)
    return B.stick_on_mat(_site(), stick, half_width=2.0, ndiv=2, mat_mass=20.0)


@pytest.fixture(scope="module")
def base_dir(tmp_path_factory):
    """SITE + POINT + HOUSE of a small stick-on-mat model, run once."""
    wd = tmp_path_factory.mktemp("analys_base")
    mdl = _model()
    B.run_soil(wd, "m", _site(), FS, layer=0, rad=mdl.rad)
    B.run_house(wd, "m", mdl)
    return wd, mdl


@pytest.fixture
def wd(base_dir, tmp_path):
    src, mdl = base_dir
    dst = tmp_path / "m"
    shutil.copytree(src, dst)
    return dst, mdl


def _rc(wd, **params):
    return B.run_analys(wd, "m", params.pop("fs", FS), check=False, **params)


def _listing(wd) -> str:
    return B.listing(wd, "m", "ANALYS")


# ---------------------------------------------------------------------------------------
# reference solutions assembled independently with dense matrices
# ---------------------------------------------------------------------------------------
def _dense_system(wd, q, fnum):
    f3 = read_container(wd / "FILE3")
    f4 = read_container(wd / "m.N4")
    ck, cm = read_container(wd / "COOSK"), read_container(wd / "COOSM")
    w = 2 * math.pi * fnum * FS.df
    xyz, iface = np.asarray(f4["x_int_xyz"]), np.asarray(f4["int_iface"])
    F = flexibility_matrix(f3, frequency_row(f3, fnum), xyz[:, :2], iface)
    X = np.linalg.inv(F)
    A = (ck.sparse("Ks") - ck.sparse("Ke")).toarray() - w * w * (cm.sparse("Ms") - cm.sparse("Me")).toarray()
    eq = np.asarray(f4["int_eq"]).reshape(-1)
    act = eq >= 0
    A = A.astype(complex)
    A[np.ix_(eq[act], eq[act])] += X[np.ix_(act, act)]
    return f4, A, X, xyz, iface, eq, act


def test_seismic_solution_equals_dense_flexible_volume_equation(wd):
    """R1 4.4 general assembly solved densely: (C_s - C_e + A^T X A) U = A^T X U'."""
    wd, mdl = wd
    assert _rc(wd) == 0
    f8 = B.read_file8(wd)
    f1 = read_container(wd / "FILE1")
    for q, n in enumerate(FS.fnum):
        f4, A, X, xyz, iface, eq, act = _dense_system(wd, q, n)
        Up = free_field_at_nodes(f1, frequency_row(f1, n), xyz, iface, 0.0, 0.0, 0.0).reshape(-1)
        b = np.zeros(A.shape[0], complex)
        with np.errstate(all="ignore"):                  # spurious Accelerate matmul flags (numpy 2.0)
            b[eq[act]] = (X @ Up)[act]
        U = np.linalg.solve(A, b)
        np.testing.assert_allclose(f8["H"][q], U, rtol=0, atol=1e-9 * np.abs(U).max())


def test_vibration_solution_equals_dense(wd):
    wd, mdl = wd
    B.run_force(wd, "m", [(mdl["top"], 1, 1.0, 0.0), (mdl["centre"], 3, 2.0, 0.05)], FS)
    assert _rc(wd, type=1) == 0
    f8 = B.read_file8(wd)
    assert f8.meta["type"] == 1 and f8.meta["case"] == 0
    f9 = read_container(wd / "FILE9")
    for q, n in enumerate(FS.fnum):
        f4, A, X, *_ = _dense_system(wd, q, n)
        col = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(f4["eq_node"], f4["eq_dof"]))}
        b = np.zeros(A.shape[0], complex)
        for k, (nd, d) in enumerate(zip(f9["load_node"], f9["load_dof"])):
            b[col[(int(nd), int(d))]] += f9["P"][q, k]
        U = np.linalg.solve(A, b)
        np.testing.assert_allclose(f8["H"][q], U, rtol=0, atol=1e-9 * np.abs(U).max())


def test_fixed_base_vibration_without_interaction_nodes(tmp_path):
    """No interaction nodes: plain (K* - w^2 M) U = P, no FILE3 needed (warning in the listing)."""
    site = _site()
    hb = B.HouseBuilder(site)
    b0 = hb.add_node(0.0, 0.0, 0.0, fix=(1,) * 6)
    stick = B.Stick(heights=[3.0, 6.0], masses=[40.0, 30.0], E=3.0e7, A=1.0, I=0.5, beta=0.05)
    nodes = B.add_stick(hb, b0, stick)
    B.write_deck(tmp_path, "m", B.site_deck(site, FS, model="m"))
    B.run_house(tmp_path, "m", hb)
    B.run_force(tmp_path, "m", [(nodes[-1], 2, 1.0, 0.0)], FS)
    assert _rc(tmp_path, type=1) == 0
    assert "no interaction nodes" in _listing(tmp_path)
    ck, cm = read_container(tmp_path / "COOSK"), read_container(tmp_path / "COOSM")
    f4 = read_container(tmp_path / "m.N4")
    f8 = B.read_file8(tmp_path)
    k = int(np.flatnonzero((f4["eq_node"] == nodes[-1]) & (f4["eq_dof"] == 2))[0])
    for q, f in enumerate(FS.freq):
        A = ck.sparse("Ks").toarray() - (2 * math.pi * f) ** 2 * cm.sparse("Ms").toarray()
        b = np.zeros(A.shape[0], complex)
        b[k] = 1.0
        np.testing.assert_allclose(f8["H"][q], np.linalg.solve(A, b), rtol=1e-10, atol=1e-16)
    assert _rc(tmp_path) != 0 and "no interaction nodes" in _listing(tmp_path)


# ---------------------------------------------------------------------------------------
# options, tiers and input checks
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("params,msg", [
    # incoherency, wave passage and multiple excitation are implemented since wave 3 (package incoherency):
    # without the HOUSE factors (FILE77) the run stops naming the producer (requirements 2.6)
    (dict(coh=1), "FILE77 missing -- run HOUSE"),
    (dict(me=1), "FILE77 missing -- run HOUSE"),
    (dict(wpass=1), "FILE77 missing -- run HOUSE"),
    (dict(simul=2), "it needs incoherent motion (HOUSE <coh> = 1)"),
    (dict(mode=6), "Mode 6 (New Load Vector) is not available"),
    (dict(simul=1, ang=30.0), "coordinate transformation angle 0"),
    (dict(type=1, simul=501), "exceed the maximum 500"),
    (dict(impe=3), "<impe> = 3"),
])
def test_option_errors(wd, params, msg):
    wd, _ = wd
    assert _rc(wd, **params) == 1
    assert msg in _listing(wd)


def test_2d_model_needs_inplane_interaction_dofs(wd):
    """2D analysis is implemented (wave 3, package twod_symm; tests/unit/test_twod_modules.py): a FILE4 marked
    2D must have the in-plane interaction DOFs UX, UZ (D-W2-07), otherwise the run stops."""
    wd, _ = wd
    f4 = read_container(wd / "m.N4")
    meta = dict(f4.meta, dim=1)
    write_container(wd / "m.N4", "FILE4", f4.arrays, meta, module="HOUSE")
    assert _rc(wd) == 1
    assert "a 2D model has the interaction DOFs UX, UZ" in _listing(wd)


def test_angle_normalised_with_warning(wd):
    wd, mdl = wd
    assert _rc(wd, ang=-30.0) == 0
    assert "normalised to 330" in _listing(wd)
    f8 = B.read_file8(wd)
    assert f8.meta["ang"] == 330.0
    H330 = B.tf(f8, mdl["top"], 1)
    assert _rc(wd, ang=0.0) == 0
    np.testing.assert_allclose(H330, math.cos(math.radians(330)) * B.tf(B.read_file8(wd), mdl["top"], 1),
                               rtol=1e-8, atol=1e-12)


def test_frequency_survey_stops_on_missing_frequency(wd):
    wd, _ = wd
    assert _rc(wd, fs=FS.subset([4, 20, 999])) == 1
    text = _listing(wd)
    assert "Frequency survey" in text and "MISSING" in text
    assert "Frequency 999 not in FILE1 -- re-run SITE" in text and "Frequency 999 not in FILE3 -- re-run POINT" in text
    assert not (wd / "FILE8").exists()


def test_fopt1_uses_every_file1_frequency(wd):
    wd, _ = wd
    assert _rc(wd, fs=FS.subset([]), fopt=1) == 0
    assert list(B.read_file8(wd)["fnum"]) == FS.fnum
    assert _rc(wd, fs=FS.subset([])) == 1 and "Error 120" in _listing(wd)


def test_too_many_frequencies_and_df_mismatch(wd, monkeypatch):
    wd, _ = wd
    monkeypatch.setattr(analys, "MAX_FREQ", 2)
    assert _rc(wd) == 0 and "EDU-02: 3 SSI frequencies exceed the manual limit of 2" in _listing(wd)
    monkeypatch.setattr(analys, "MAX_FREQ", 500)
    bad = B.FrequencySet(FS.fnum, FS.df * 1.01, FS.delt, FS.nft)
    assert _rc(wd, fs=bad) == 1 and "frequency step of FILE1" in _listing(wd)


def test_other_site_is_refused(wd):
    wd, _ = wd
    f1 = read_container(wd / "FILE1")
    arr = dict(f1.arrays, depth_user=np.asarray(f1["depth_user"]) * 1.5)
    write_container(wd / "FILE1", "FILE1", arr, f1.meta, module="SITE")
    assert _rc(wd) == 1 and "another site layering" in _listing(wd)


def test_data_check_writes_nothing(wd):
    wd, _ = wd
    assert _rc(wd, opmode=1) == 0
    assert "Data check only" in _listing(wd)
    assert not (wd / "FILE8").exists()


def test_resource_check_refuses_oversized_runs(wd, monkeypatch):
    wd, _ = wd
    monkeypatch.setattr(analys, "physical_memory", lambda: 1000)
    assert _rc(wd) == 1 and "EDU-20" in _listing(wd)


def test_file8_schema_and_meta(wd):
    wd, mdl = wd
    assert _rc(wd, save=1) == 0
    f8 = B.read_file8(wd)
    f4 = read_container(wd / "m.N4")
    assert validate(f8) == []
    assert f8.meta["model_hash"] == f4.meta["model_hash"]
    assert f8.meta["nfft"] == FS.nft and f8.meta["delt"] == FS.delt and f8.meta["cm"] == 0
    np.testing.assert_array_equal(f8["eq_node"], f4["eq_node"])
    np.testing.assert_allclose(f8["freq"], FS.freq)
    assert set(np.asarray(f8["x_paths"]).tolist()) <= {0, 1, 2}
    text = _listing(wd)
    for key in ("Frequency survey", "Resource check", "Per-frequency summary", "Low-frequency check",
                "Transfer functions FILE8", "COOXI", "COOTKI"):
        assert key in text


def test_print_options_and_node_limit(wd, monkeypatch):
    wd, _ = wd
    monkeypatch.setattr(analys, "MAX_PRINT_NODES", 2)
    assert _rc(wd, prnt=0) == 0
    text = _listing(wd)
    assert "real and imaginary parts" in text and "the first 2 of" in text
    blk = text[text.index("Transfer functions FILE8"):]
    assert "  Re" in blk and "  Im" in blk


def test_full_fallback_path_in_module(wd, monkeypatch):
    wd, _ = wd
    assert _rc(wd) == 0
    ref = B.read_file8(wd)["H"].copy()
    monkeypatch.setattr(ss, "RCOND_MIN", 2.0)
    assert _rc(wd) == 0
    assert "full sparse LU" in _listing(wd)
    np.testing.assert_allclose(B.read_file8(wd)["H"], ref, rtol=0, atol=1e-10 * np.abs(ref).max())


def test_mode_type_mismatch_warns(wd):
    wd, mdl = wd
    assert _rc(wd, save=1) == 0
    B.run_force(wd, "m", [(mdl["top"], 1, 1.0, 0.0)], FS)
    assert _rc(wd, type=1, mode=2) == 0
    assert "solved as New Dynamic Loading" in _listing(wd)


# ---------------------------------------------------------------------------------------
# restarts (requirements 2.5, D-ANL-05)
# ---------------------------------------------------------------------------------------
def test_restart_needs_index_files(wd):
    wd, _ = wd
    assert _rc(wd, mode=1) == 1 and "COOXI missing" in _listing(wd)
    assert _rc(wd, save=1, fs=FS.subset([4, 20])) == 0
    assert _rc(wd, mode=1) == 1 and "Frequency 60 not in COOXI" in _listing(wd)
    assert _rc(wd, mode=2, fs=FS.subset([20])) == 0


def test_restart_hash_validation(wd):
    wd, _ = wd
    assert _rc(wd, save=1) == 0
    c = read_container(wd / "COOXI")
    write_container(wd / "COOXI", "COOXI", c.arrays, dict(c.meta, int_hash="other"), module="ANALYS")
    assert _rc(wd, mode=1) == 1 and "different interaction nodes" in _listing(wd)
    write_container(wd / "COOXI", "COOXI", c.arrays, c.meta, module="ANALYS")
    x = read_container(wd / "COOX002")
    write_container(wd / "COOX002", "COOX", x.arrays, dict(x.meta, fnum=999), module="ANALYS")
    assert _rc(wd, mode=1) == 1 and "the restart files were overwritten" in _listing(wd)


def test_restart_index_must_have_the_run_frequency_step(wd):
    """D-ANL-05 / requirements 4.6 item 8: records are keyed by frequency value.  An index saved at
    another df (NFFT or time step changed since the initiation) is refused before any record is used,
    as is an index without df or with an entry at another frequency; the same df reached through
    another route (deck df = 1/(NFFT dt) given as a free step) is accepted."""
    wd, _ = wd
    assert _rc(wd, save=1) == 0
    ref = B.read_file8(wd)["H"].copy()
    c = read_container(wd / "COOXI")
    write_container(wd / "COOXI", "COOXI", c.arrays, dict(c.meta, df=2.0 * FS.df), module="ANALYS")
    assert _rc(wd, mode=1) == 1
    text = _listing(wd)
    assert "COOXI: the restart files were saved for the frequency step df = 0.1953125 Hz" in text
    assert "this run uses df = 0.09765625 Hz" in text
    meta = {k: v for k, v in c.meta.items() if k != "df"}
    write_container(wd / "COOXI", "COOXI", c.arrays, meta, module="ANALYS")
    assert _rc(wd, mode=1) == 1 and "COOXI does not record the frequency step" in _listing(wd)
    arr = dict(c.arrays, freq=np.asarray(c["freq"]) * (1.0 + 1e-3))
    write_container(wd / "COOXI", "COOXI", arr, c.meta, module="ANALYS")
    assert _rc(wd, mode=1) == 1 and "COOXI: frequency number 4 was saved at f =" in _listing(wd)
    write_container(wd / "COOXI", "COOXI", c.arrays, c.meta, module="ANALYS")
    t = read_container(wd / "COOTKI")
    write_container(wd / "COOTKI", "COOTKI", t.arrays, dict(t.meta, df=0.5 * FS.df), module="ANALYS")
    assert _rc(wd, mode=2) == 1 and "COOTKI: the restart files were saved for the frequency step" in _listing(wd)
    write_container(wd / "COOTKI", "COOTKI", t.arrays, t.meta, module="ANALYS")
    same = B.FrequencySet.harmonic(1.0 / (FS.delt * FS.nft), FS.fnum, nft=FS.nft)
    assert _rc(wd, mode=2, fs=same) == 0
    np.testing.assert_allclose(B.read_file8(wd)["H"], ref, rtol=0, atol=1e-10 * np.abs(ref).max())


@pytest.mark.parametrize("record", ["COOX002", "COOTK002"])
def test_restart_record_must_hold_the_frequency_of_its_number(wd, record):
    """Every record stores f (requirements 4.6 item 8); a record whose f is not fnum x df of the run is
    refused even when the index agrees (a record overwritten by another initiation)."""
    wd, _ = wd
    assert _rc(wd, save=1) == 0
    kind = "COOX" if record.startswith("COOX0") else "COOTK"
    x = read_container(wd / record)
    write_container(wd / record, kind, x.arrays, dict(x.meta, freq=2.0 * float(x.meta["freq"])), module="ANALYS")
    assert _rc(wd, mode=2) == 1
    text = _listing(wd)
    assert f"{record} holds frequency f = " in text and "frequency number 20 of this run is" in text
    write_container(wd / record, kind, x.arrays, {k: v for k, v in x.meta.items() if k != "freq"}, module="ANALYS")
    assert _rc(wd, mode=2) == 1 and f"{record} holds frequency f = None Hz" in _listing(wd)


def test_restart_modes_reproduce_initiation(wd):
    wd, _ = wd
    assert _rc(wd, save=1) == 0
    ref = B.read_file8(wd)["H"].copy()
    for mode in (1, 2):
        assert _rc(wd, mode=mode) == 0
        np.testing.assert_allclose(B.read_file8(wd)["H"], ref, rtol=0, atol=1e-10 * np.abs(ref).max())
    assert _rc(wd, mode=2, fs=FS.subset([60, 4])) == 0
    f8 = B.read_file8(wd)
    assert list(f8["fnum"]) == [4, 60]
    np.testing.assert_allclose(f8["H"], ref[[0, 2]], rtol=0, atol=1e-10 * np.abs(ref).max())
    assert "factorisation reused (COOTK)" in _listing(wd)


def test_restart_file_pattern_covers_four_digit_order_numbers(tmp_path):
    """More than 999 saved frequencies give COOX1000 / COOTK1000 (conventions.restart_names); Delete
    Restart Files must remove them too, and nothing of HOUSE (COOSK, COOSM) or other files."""
    from sassi.conventions import restart_names
    assert restart_names(1000) == ("COOX1000", "COOTK1000")
    keep = ["COOSK", "COOSM", "COOX", "COOX12", "COOXII", "COOTKX", "FILE8", "COOX001.tmp"]
    gone = ["COOX001", "COOTK001", "COOX999", "COOX1000", "COOTK1000", "COOX12345", "COOXI", "COOTKI"]
    for name in keep + gone:
        (tmp_path / name).write_bytes(b"x")
    deleted = analys._delete_restart_files(SimpleNamespace(workdir=tmp_path))
    assert sorted(deleted) == sorted(gone)
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(keep)


def test_delete_restart_files(wd):
    wd, _ = wd
    assert _rc(wd, save=1) == 0
    assert (wd / "COOX001").exists() and (wd / "COOTK003").exists()
    assert _rc(wd, mode=2, delrst=1, save=1) == 0
    assert (wd / "COOX001").exists() and "restart files are kept" in _listing(wd)
    assert _rc(wd, mode=2, delrst=1) == 0
    assert not any((wd / n).exists() for n in ("COOX001", "COOTK003", "COOXI", "COOTKI"))
    assert (wd / "FILE8").exists()


# ---------------------------------------------------------------------------------------
# loads, impedance, listing checks
# ---------------------------------------------------------------------------------------
def test_vibration_loads_on_fixed_dofs_and_unknown_nodes(wd):
    wd, mdl = wd
    hb = mdl.house
    knode = next(n for n, f in hb.fixity.items() if all(f) and n not in hb.interaction)
    B.run_force(wd, "m", [(mdl["top"], 1, 1.0, 0.0), (knode, 1, 1.0, 0.0)], FS)
    assert _rc(wd, type=1) == 0
    assert f"load on node {knode} UX ignored" in _listing(wd)
    B.run_force(wd, "m", [(9999, 1, 1.0, 0.0)], FS)
    assert _rc(wd, type=1) == 1 and "Error 62" in _listing(wd)


def test_simultaneous_vibration_cases_need_their_files(wd):
    wd, mdl = wd
    B.run_force(wd, "m", [(mdl["top"], 1, 1.0, 0.0)], FS, copy_to="FILE9001")
    assert _rc(wd, type=1, simul=2) == 1 and "FILE9002 missing" in _listing(wd)
    B.copy_file(wd, "FILE9001", "FILE9002")
    assert _rc(wd, type=1, simul=2) == 0
    a, b = B.read_file8(wd, "FILE8001"), B.read_file8(wd, "FILE8002")
    assert np.array_equal(a["H"], b["H"]) and a.meta["case"] == 1 and b.meta["case"] == 2


@pytest.mark.parametrize("impe", [1, 2])
def test_global_impedance_outputs(wd, impe):
    """D-ANL-07: K_G = T^T X_ff T; FOUNSTIF = Re K, FOUNDASH = Im K / w, FOUNDAMP = Im K/(2|Re K|),
    FOUNIMPD = |K|; FILE11 holds X_ff, T and K_G."""
    wd, _ = wd
    xc, yc, zc = 0.3, -0.2, 0.0
    assert _rc(wd, impe=impe, save=1, xc=xc, yc=yc, zc=zc) == 0
    c11 = read_container(wd / "FILE11", "FILE11")
    f4 = read_container(wd / "m.N4")
    T = ss.rigid_body_transform(f4["x_int_xyz"], (xc, yc, zc))
    np.testing.assert_allclose(c11["T"], T)
    for q in range(len(FS.fnum)):
        X = read_container(wd / f"COOX{q + 1:03d}")["X"]
        np.testing.assert_allclose(c11["X"][q], X)
        np.testing.assert_allclose(c11["KG"][q], T.T @ X @ T, rtol=1e-12, atol=1e-9 * np.abs(c11["KG"][q]).max())
    KG = c11["KG"]
    w = 2 * np.pi * FS.freq
    files = {}
    for name in ("FOUNSTIF", "FOUNDASH", "FOUNDAMP", "FOUNIMPD"):
        rows = [ln.split() for ln in (wd / name).read_text().splitlines() if not ln.startswith("#")]
        files[name] = np.array([[float(v) for v in r] for r in rows])
        assert files[name].shape == (len(FS.fnum), 7 if impe == 1 else 37)
        np.testing.assert_allclose(files[name][:, 0], FS.freq, rtol=1e-12)
    sel = (lambda A: np.diagonal(A, axis1=1, axis2=2)) if impe == 1 else (lambda A: A.reshape(len(FS.fnum), 36))
    np.testing.assert_allclose(files["FOUNSTIF"][:, 1:], sel(KG.real), rtol=1e-14)
    np.testing.assert_allclose(files["FOUNDASH"][:, 1:], sel(KG.imag / w[:, None, None]), rtol=1e-14)
    np.testing.assert_allclose(files["FOUNIMPD"][:, 1:], sel(np.abs(KG)), rtol=1e-14)
    d = np.diagonal(KG, axis1=1, axis2=2)
    damp = files["FOUNDAMP"][:, 1:]
    dd = damp if impe == 1 else np.diagonal(damp.reshape(len(FS.fnum), 6, 6), axis1=1, axis2=2)
    np.testing.assert_allclose(dd, d.imag / (2 * np.abs(d.real)), rtol=1e-14)
    assert np.all(d.real[:, :3] > 0)                       # translational stiffnesses positive
    assert "Global soil impedance" in _listing(wd)


def test_global_impedance_needs_interaction_nodes(tmp_path):
    site = _site()
    hb = B.HouseBuilder(site)
    b0 = hb.add_node(0.0, 0.0, 0.0, fix=(1,) * 6)
    B.add_stick(hb, b0, B.Stick(heights=[3.0], masses=[10.0], E=3e7, A=1.0, I=0.5))
    B.write_deck(tmp_path, "m", B.site_deck(site, FS, model="m"))
    B.run_house(tmp_path, "m", hb)
    B.run_force(tmp_path, "m", [(2, 1, 1.0, 0.0)], FS)
    assert _rc(tmp_path, type=1, impe=1) == 1 and "no interaction nodes" in _listing(tmp_path)


def test_low_frequency_check_warns():
    """G-19 / EDU-18: a translational ATF in the control direction far from 1 at f_1 is reported."""
    lst = Listing(None)
    opt = SimpleNamespace(seismic=True, df=0.1, ang=0.0)
    md = SimpleNamespace(eq_dof=np.array([1, 2, 1, 5]), eq_node=np.array([1, 1, 2, 2]))
    case = SimpleNamespace(outfile="FILE8", cm=0)
    H = {"FILE8": np.array([[1.01, 0.0, 1.2, 0.3]])}
    analys._lowfreq_check(lst, opt, md, [case], H, [1])
    assert any("EDU-18" in w and "2 UX" in w for w in lst.warnings)
    lst2 = Listing(None)
    analys._lowfreq_check(lst2, opt, md, [case], {"FILE8": np.array([[1.01, 0.0, 0.98, 0.3]])}, [1])
    assert not lst2.warnings


# ---------------------------------------------------------------------------------------
# end to end through the interpreter: AFWRITE writes the ANALYS deck, RUNANALYS runs it
# ---------------------------------------------------------------------------------------
E2E_MODEL = """
MDL,demo,{dir}
TIT,ANALYS end-to-end
N,1,0,0,-5
N,2,10,0,-5
N,3,10,10,-5
N,4,0,10,-5
N,5,0,0,0
N,6,10,0,0
N,7,10,10,0
N,8,0,10,0
N,10,5,5,0
N,11,5,5,10
N,12,15,5,0
INT,1,8,1,1
M,1,4.32e5,0.25,0.15,0.05,0.05
R,1,1,0.8,0.8,0.2,0.1,0.1
L,1,5,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
TOPL,1,1
GROUP,1,SOLID
E,1,1,2,3,4,5,6,7,8
GROUP,2,BEAMS
E,1,10,11,12
GROUP,3,SPRING
SC,1,1e4,1e4,1e4,1e8,1e8,1e8,0.02
RACT,1
E,1,7,10
FREQ,1,2,4,8
SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,4096,1
POINT,0,1,4.5
HOUSE,32.2,0,0,2,0,0,0,0,0
MT,11,32.2,32.2,32.2
D,7,7,1,1,ROT
ANALYS,0,0,0,1,1,0,0,0,0,0,1
AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
"""


def test_end_to_end_interpreter_runanalys(tmp_path):
    """.pre -> AFWRITE -> RUNSITE/RUNPOINT/RUNHOUSE/RUNANALYS: the AFWRITE ANALYS deck (frequency
    set, df, options) is what ANALYS needs; global impedance about the origin is written.  (The SPRING
    clamps the beam to the SOLID node 7, whose rotations exist only through the spring: they are fixed
    with D, EDU-06.)"""
    from sassi.prep import Interpreter, Kind
    mdir = tmp_path / "demo"
    mdir.mkdir()
    ui = Interpreter(cwd=tmp_path)
    ui.run_text(E2E_MODEL.format(dir=mdir))
    assert not ui.sink.texts(Kind.ERROR), ui.sink.texts(Kind.ERROR)
    out = (mdir / "demo_ANALYS.out").read_text()
    assert "ANALYS finished with status OK" in out, out[-3000:]
    f8 = read_container(mdir / "FILE8", "FILE8")
    assert list(f8["fnum"]) == [2, 4, 8] and np.isclose(f8.meta["df"], 1 / (0.005 * 4096))
    assert (mdir / "COOX003").exists() and (mdir / "FOUNSTIF").exists()
    H11 = B.tf(f8, 11, 1)
    assert np.all(np.isfinite(H11)) and abs(abs(H11[0]) - 1.0) < 0.05        # 0.1 Hz: rigid-body motion


def test_singular_system_names_the_unrestrained_dofs(tmp_path):
    """A SPRING without rotational constants adds rotations to a SOLID-only node: no stiffness, no
    mass -> singular system; the error names the equations (EDU-06 hint)."""
    from sassi.prep import Interpreter, Kind
    mdir = tmp_path / "demo"
    mdir.mkdir()
    ui = Interpreter(cwd=tmp_path)
    text = E2E_MODEL.replace("SC,1,1e4,1e4,1e4,1e8,1e8,1e8,0.02", "SC,1,1e4,1e4,1e4,0,0,0,0.02")
    text = text.replace("D,7,7,1,1,ROT\n", "")
    ui.run_text(text.format(dir=mdir))
    out = (mdir / "demo_ANALYS.out").read_text()
    assert "singular" in out and "node 7 ROTX" in out and "EDU-06" in out


def test_batch_protocol(wd, monkeypatch):
    """D-RUN-04: python -m sassi.modules.analys < ANALYS.inp (model, deck, listing)."""
    import io
    from sassi.modules.base import batch_main
    wd, _ = wd
    B.write_deck(wd, "m", B.analys_deck(FS, model="m"))
    monkeypatch.chdir(wd)
    rc = batch_main("ANALYS", stdin=io.StringIO("m\nm.anl\nbatch_ANALYS.out\n"))
    assert rc == 0 and (wd / "FILE8").exists()
    assert "ANALYS finished with status OK" in (wd / "batch_ANALYS.out").read_text()


def test_frequency_subsets_merge_with_combin(wd):
    """Requirements 2.4 frequency refinement: ANALYS on two subsets, FILE81 + FILE82 -> COMBIN ->
    the FILE8 of the full set (ANALYS -> COMBIN contract, VP-25 with real ANALYS files)."""
    wd, _ = wd
    assert _rc(wd) == 0
    full = B.read_file8(wd)
    assert _rc(wd, fs=FS.subset([4, 60])) == 0
    B.copy_file(wd, "FILE8", "FILE81")
    assert _rc(wd, fs=FS.subset([20])) == 0
    B.copy_file(wd, "FILE8", "FILE82")
    (wd / "FILE8").unlink()
    assert B.run("COMBIN", wd, "m") == 0
    merged = B.read_file8(wd)
    assert list(merged["fnum"]) == FS.fnum
    np.testing.assert_array_equal(merged["H"], full["H"])
