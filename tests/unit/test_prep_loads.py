"""Loads and masses (spec 08 section 7): UT-15 mass units, MOPT add/overwrite (D-MDL-08), generation,
scaling with "last two defined" defaults, deletion and listing."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.prep import Interpreter, Kind


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def run(ui, text):
    ui.run_text(text)
    return ui.model


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


# ------------------------------------------------------------------ UT-15
def test_ut15_munits_weight_equals_mass(ui, tmp_path):
    m = run(ui, "GRAVITY,32.2\nN,1,0,0,0\nN,2,0,0,0\nMT,1,32.2,32.2,32.2\nMR,1,64.4,64.4,64.4\nMUNITS,1,1,1,1\n"
                "MT,2,1,1,1\nMR,2,2,2,2\nMUNITS,2,2,1,0\n")
    mm = m.nodal_masses()
    assert np.array_equal(mm[1], mm[2])
    assert mm[2].tolist() == [1, 1, 1, 2, 2, 2]
    # the default flag of a node never given MUNITS is 1 (weight), D-MDL-08
    assert m.mass_unit(1) == 1 and 1 not in m.mass_units and m.mass_units[2] == 0


def test_ut15_general_matrix_weight_units(tmp_path):
    a = Interpreter(cwd=tmp_path)
    a.run_text("GRAVITY,32.2\nMOPT,1,1,1,1\nMXM,1,1,32.2\nMXM,1,7,64.4\nMXR,1,1,1000")
    b = Interpreter(cwd=tmp_path)
    b.run_text("GRAVITY,32.2\nMOPT,1,0,1,1\nMXM,1,1,1\nMXM,1,7,2\nMXR,1,1,1000")
    Ka, Ma = a.model.general_matrices(1)
    Kb, Mb = b.model.general_matrices(1)
    assert np.array_equal(Ka, Kb)
    assert np.array_equal(Ma, Mb)
    assert Mb[0, 0] == 1.0 and Mb[6, 6] == 2.0


# ------------------------------------------------------------------ forces and moments
def test_force_overwrite_and_add_modes(ui):
    m = run(ui, "F,50,1.0,0,0\nF,51,0.5,0,0,0.1\nF,50,2,0,0\n")
    assert m.forces[50].factor == [2.0, 0.0, 0.0]                 # MOPT default: overwrite
    assert m.forces[51].arrival == [0.1, 0.0, 0.0]
    run(ui, "MOPT,1,0,1,0\nF,50,1,3,0,0.5,0.2\nMM,13,0,0,1,0,0,0.2\nMM,13,0,0,1\n")
    assert m.forces[50].factor == [3.0, 3.0, 0.0]
    assert m.forces[50].arrival == [0.5, 0.2, 0.0]
    assert any("arrival time 0 replaced by 0.5" in w for w in warnings(ui))
    assert m.moments[13].factor == [0.0, 0.0, 2.0]
    assert m.force_history == [51, 50]


def test_fscale_last_two_defined_and_zero_rule(ui):
    m = run(ui, "F,10,1,1,1\nF,30,1,1,1\nF,20,1,1,1\nFSCALE,,,,2\n")
    # last two defined: 30 and 20 -> range 20..30
    assert m.forces[20].factor == [2, 1, 1] and m.forces[30].factor == [2, 1, 1]
    assert m.forces[10].factor == [1, 1, 1]
    run(ui, "FSCALE,10,30,10,0,3,0\nMM,5,1,1,1,0.3\nMSCALE,5,,,4,4,4\n")
    assert m.forces[10].factor == [1, 3, 1] and m.forces[30].factor == [2, 3, 1]
    assert m.moments[5].factor == [4, 4, 4] and m.moments[5].arrival == [0.3, 0, 0]


def test_fdel_mmdel_and_lists(ui):
    m = run(ui, "F,1,1\nF,2,1\nF,3,1\nMM,1,1\nFDEL,1,3,2\nMMDEL,1\nFLIST\nMMLIST\n")
    assert sorted(m.forces) == [2] and not m.moments
    assert m.force_history == [2]
    info = ui.sink.texts(Kind.INFO)
    assert "1 forces listed" in info and "0 moments listed" in info


# ------------------------------------------------------------------ masses
def test_mass_add_mode_and_generation(ui):
    m = run(ui, "MT,100,1,1,1\nMT,100,2,2,2\n")
    assert m.tmass[100] == [2, 2, 2]
    run(ui, "MOPT,1,0,0,1\nMT,100,1,1,1\n")
    assert m.tmass[100] == [3, 3, 3]
    run(ui, "MOPT\nN,100,0,0,0\nMT,101,1,1,1\nMUNITS,100,101,1,0\nMTGEN,3,10,100,101,1,0.5\n")
    for k in (1, 2, 3):
        assert m.tmass[100 + 10 * k] == [3 + 0.5 * k, 3, 3]
        assert m.tmass[101 + 10 * k] == [1 + 0.5 * k, 1, 1]
        assert m.mass_unit(100 + 10 * k) == 0                      # units copied from the pattern
    run(ui, "MR,5,1,2,3\nMRGEN,2,,5,5,,1,1,1\n")                   # ninc default n2-n1+1 = 1
    assert m.rmass[6] == [2, 3, 4] and m.rmass[7] == [3, 4, 5]


def test_mass_scale_delete_list(ui):
    m = run(ui, "MT,1,1,1,1\nMT,2,1,1,1\nMT,3,1,1,1\nMTSCALE,,,,0,0,5\nMR,4,1,1,1\nMRSCALE,4,4,1,2\n")
    assert m.tmass[2] == [1, 1, 5] and m.tmass[3] == [1, 1, 5] and m.tmass[1] == [1, 1, 1]
    assert m.rmass[4] == [2, 1, 1]
    run(ui, "MUNITS,3,3,1,0\nMTDEL,1,3\nMTLIST\n")
    assert not m.tmass and 3 not in m.mass_units                   # orphan flag dropped
    assert "1 nodes with masses listed" in ui.sink.texts(Kind.INFO)


def test_munits_validation_and_default_normalised(ui):
    m = run(ui, "N,1,0,0,0\nMUNITS,1,1,1,0\nMUNITS,1,1,1,1\nMUNITS,1,1,1,2\n")
    assert m.mass_units == {}                                       # setting the default removes the entry
    assert any("must be 0 (mass) or 1 (weight)" in e for e in ui.sink.texts(Kind.ERROR))


def test_mtlist_shows_converted_mass(ui):
    run(ui, "GRAVITY,32.2\nMT,1,32.2,0,0\nMTLIST\n")
    info = ui.sink.texts(Kind.INFO)
    assert any("weight" in t and "32.2" in t for t in info)
    assert any(t.strip().startswith("mass") and " 1 " in f" {t} " for t in info)


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


# ------------------------------------------------------------------ generated ids, robustness, speed
def test_mtgen_refuses_nonpositive_targets(ui):
    m = run(ui, "MT,3,1,1,1\nMR,3,1,1,1")
    ui.execute("MTGEN,1,-5,3,3")
    ui.execute("MRGEN,2,-2,3,3")                     # second copy would be node -1
    assert sorted(m.tmass) == [3] and sorted(m.rmass) == [3]
    assert sum("must be positive" in e for e in errors(ui)) == 2
    ui.execute("MTGEN,2,-1,3,3")                     # 2 and 1 are valid
    assert sorted(m.tmass) == [1, 2, 3]


def test_mtlist_with_zero_gravity_lists_raw_values(ui):
    run(ui, "GRAVITY,0\nMT,1,1,2,3\nMT,2,5,5,5\nMUNITS,2,2,1,0\nMTLIST")
    assert not any("internal error" in e for e in errors(ui))
    assert any("not available" in w for w in warnings(ui))
    info = ui.sink.texts(Kind.INFO)
    assert "2 nodes with masses listed" in info


def test_munits_selection_nodes_and_mass_only_ids(ui):
    """MUNITS selects defined nodes and ids that only carry a mass (range short and long)."""
    m = run(ui, "N,1\nN,5\nMT,7,1,1,1\nMR,9,1,1,1\nMUNITS,1,9,1,0")
    assert sorted(m.mass_units) == [1, 5, 7, 9]
    run(ui, "MUNITS,1,100000,4,1")                    # range longer than the id count: 1, 5, 9
    assert sorted(m.mass_units) == [7]
    run(ui, "MUNITS,6,6,1,0")                         # 6 is neither a node nor a mass: nothing set
    assert 6 not in m.mass_units


def test_munits_cost_does_not_grow_with_model_size(ui):
    """One MUNITS line costs O(range), not O(number of nodes): INP of a large WRITE stays linear."""
    import time
    from sassi.prep.commands import UnionIds
    m = run(ui, "\n".join(f"N,{i}" for i in range(1, 20001)))
    u = UnionIds(m.nodes, m.tmass, m.rmass)
    assert 7 in u and 20001 not in u and u.size_bound() == 20000 and len(u) == 20000
    t0 = time.perf_counter()
    ui.run_text("\n".join(f"MUNITS,{i},{i},1,0" for i in range(1, 20001, 2)))
    dt = time.perf_counter() - t0
    assert len(m.mass_units) == 10000
    assert dt < 2.0, f"10000 MUNITS lines took {dt:.2f} s"
