"""Model database (sassi.model): geometry helpers, SSIModel derived views, serialisation and hash,
option records (generic and typed), materials helpers."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from sassi.model import (Field, MoptRecord, OptionStore, Record, SSIModel, make_record, register_record_type)
from sassi.model.entities import Element, Group, MatrixProp
from sassi.model.geometry import cosd, euler_from_matrix, loc_matrix, local_from_points, sind
from sassi.model.options import RECORD_TYPES
from sassi.prep import Interpreter


# ------------------------------------------------------------------ geometry
def test_trig_exact_at_right_angles():
    assert (cosd(90), sind(90), cosd(180), sind(270), cosd(-90), sind(-90)) == (0.0, 1.0, -1.0, -1.0, 0.0, -1.0)
    assert cosd(45) == pytest.approx(math.sqrt(0.5))


@pytest.mark.parametrize("angles", [(30, 20, 10), (-120, 45, 170), (0, 0, 0), (10, 90, 0), (25, -90, 0)])
def test_euler_round_trip(angles):
    R = loc_matrix(*angles)
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-14) and np.linalg.det(R) == pytest.approx(1.0)
    back = loc_matrix(*euler_from_matrix(R))
    assert np.allclose(back, R, atol=1e-13)


def test_local_from_points():
    R = local_from_points((0, 0, 0), (2, 0, 0), (1, 5, 0))
    assert np.allclose(R, np.eye(3))
    with pytest.raises(ValueError):
        local_from_points((0, 0, 0), (0, 0, 0), (1, 0, 0))
    with pytest.raises(ValueError):
        local_from_points((0, 0, 0), (1, 0, 0), (2, 0, 0))


# ------------------------------------------------------------------ SSIModel helpers
@pytest.fixture
def model(tmp_path):
    ui = Interpreter(cwd=tmp_path)
    ui.run_text("LOC,1,0,10,0,0,90,0,0\nN,1,0,0,0\nN,2,1,0,0\nCSYS,1\nN,3,1,0,0\nN,4,0,1,5\nCSYS,0\n"
                "GROUP,1,SOLID\nE,1,1,2,3,4,1,2,3,4\nGROUP,2,BEAMS\nE,1,1,2,3\nE,2,2,3,4\nGROUP,3,SHELL\nE,1,1,2,3\n"
                "GRAVITY,32.2\nMT,4,32.2,0,0\nMUNITS,4,4,1,1\nMR,3,1,1,1\nMUNITS,3,3,1,0\n")
    return ui.model


def test_global_coordinates_vectorised(model):
    ids, xyz = model.global_coordinates()
    assert ids.tolist() == [1, 2, 3, 4]
    for i, n in enumerate(ids):
        assert np.allclose(xyz[i], model.node_global(int(n)))
    assert np.allclose(model.node_global(3), [10, 1, 0])
    assert np.allclose(model.node_global(4), [9, 0, 5])
    assert np.allclose(model.node_in_system(1, 1), [0, 10, 0])


def test_element_tables(model):
    by = model.elements_by_type()
    assert sorted(by) == [1, 2, 3] and len(by[2]) == 2
    g, e, nodes = model.element_node_table(2)
    assert g.tolist() == [2, 2] and e.tolist() == [1, 2] and nodes.tolist() == [[1, 2, 3], [2, 3, 4]]
    _, _, sh = model.element_node_table(3)
    assert sh.shape == (1, 3)
    assert model.element_node_table(7)[2].shape == (0, 0)
    assert model.used_nodes() == [1, 2, 3, 4] and model.n_elements() == 4


def test_nodal_masses_in_mass_units(model):
    mm = model.nodal_masses()
    assert mm[4].tolist() == [1, 0, 0, 0, 0, 0]
    assert mm[3].tolist() == [0, 0, 0, 1, 1, 1]
    assert model.nodal_masses(gravity=16.1)[4][0] == 2.0


def test_serialisation_copy_and_hash(model):
    d = model.to_json()
    back = SSIModel.from_json(json.loads(json.dumps(d)))
    assert back.to_json() == d
    cp = model.copy()
    cp.nodes[1].x = 7.0
    assert model.nodes[1].x == 0.0
    h = model.model_hash()
    assert len(h) == 64
    model.name, model.title, model.path = "other", "new title", "/x"
    assert model.model_hash() == h                       # identity and title do not change the content hash
    model.nodes[1].z = 1e-9
    assert model.model_hash() != h
    d["schema"] = "2.0"
    with pytest.raises(ValueError):
        SSIModel.from_json(d)


def test_definition_history():
    from sassi.model import History
    h = History([3, 1])
    h.touch(5)
    h.touch(3)                                           # redefinition moves to the end
    assert h == [1, 5, 3] and h[-1] == 3 and h[-2] == 5 and h[0] == 1 and len(h) == 3 and 5 in h
    h.discard([5, 99])
    assert list(h) == [1, 3] and h == History([1, 3])


def test_frequency_step_defaults():
    m = SSIModel()
    assert m.frequency_step() == pytest.approx(1 / (0.005 * 4096))
    m.options.set_record(make_record("SITE", ["0", "1", "0.5"]))
    assert m.frequency_step() == 0.5
    m.options.set_record(make_record("SITE", ["0", "1", "0", None, None, None, None, None, None, None, None,
                                              "0.01", "2048"]))
    assert m.frequency_step() == 0.048828125


def test_matrix_prop():
    p = MatrixProp(1)
    p.set_row("R", 1, [1, 2])
    p.set_row("R", 12, [5])
    A = p.full("R")
    assert A[0, 0] == 1 and A[0, 1] == 2 and A[1, 0] == 2 and A[11, 11] == 5
    assert len(p.rows["R"][1]) == 12 and len(p.rows["R"][12]) == 1
    assert np.array_equal(p.mass(1, 2.0), p.full("M") / 2.0)
    assert MatrixProp.from_json(p.to_json()).to_json() == p.to_json()


def test_group_element_json():
    g = Group(3, 2, "frame, east")
    g.elements[2] = Element(2, [1, 2, 3], 4, 5, 1, 0, 0.0, [0, 0, 0, 0, 1, 1], [0] * 6)
    g.elements[1] = Element(1, [2, 3, 4])
    back = Group.from_json(g.to_json())
    assert [e.id for e in back.sorted_elements()] == [1, 2] and back.elements[2].ki[-1] == 1
    assert back.type_name == "BEAMS"


# ------------------------------------------------------------------ option records
def test_generic_record_positional_access():
    r = make_record("SITE", ["0", "1", "", "20"])
    assert isinstance(r, Record) and r.values == ["0", "1", None, "20"]  # typed subclass once sassi.prep.options is imported
    assert r.arg(4) == "20" and r.arg(3) is None and r.arg(99) is None
    assert r.integer(4) == 20 and r.number(3, 7.5) == 7.5 and r.number(1) == 0.0
    r.set_arg(14, 3)
    assert len(r) == 14 and r.to_tokens()[-1] == "3"
    r.set_arg(14, None)
    assert len(r) == 4
    with pytest.raises(IndexError):
        r.arg(0)
    assert r == make_record("SITE", ["0", "1", None, "20"]) and r != make_record("SITE", ["0"])
    assert Record.from_json(r.to_json()) == r


def test_typed_record_registration():
    @register_record_type
    class ZZTestRecord(Record):
        COMMAND = "ZZTESTREC"
        FIELDS = (Field("a", int, 1), Field("b", float, 0.5), Field("name", str, "x"))
    try:
        r = make_record("ZZTESTREC", ["", "2.5"])
        assert isinstance(r, ZZTestRecord)
        assert (r.a, r.b, r.name) == (1, 2.5, "x")
        r.set("name", "abc")
        r.set("a", 3)
        assert r.to_tokens() == ["3", "2.5", "abc"]
        assert not r.is_default()
        assert make_record("ZZTESTREC", ["1"]).is_default()
        assert type(Record.from_json(r.to_json())) is ZZTestRecord
        with pytest.raises(AttributeError):
            r.nothere
        with pytest.raises(KeyError):
            r.get("nothere")
    finally:
        RECORD_TYPES.pop("ZZTESTREC", None)


def test_mopt_defaults():
    r = MoptRecord()
    assert (r.incomp, r.matrix, r.mass, r.force) == (1, 0, 1, 1)
    assert SSIModel().mopt.force == 1


def test_option_store_serialisation():
    s = OptionStore()
    s.set_record(make_record("HOUSE", ["32.2"]))
    s.set_entry("DYNP", ("Rock", 2), make_record("DYNP", ["2", "0.1", "1", "0.1", "1", "Rock"]))
    s.set_entry("DYNP", ("Rock", 1), make_record("DYNP", ["1", "0.1", "1", "0.1", "1", "Rock"]))
    s.set_entry("WAVE", 2, make_record("WAVE", ["2", "1"]))
    s.set_string("THFILE", "a.acc")
    back = OptionStore.from_json(json.loads(json.dumps(s.to_json())))
    assert back == s and back.copy() == s
    assert [k for k, _ in back.entries("DYNP")] == [("Rock", 1), ("Rock", 2)]
    assert back.names() == ["DYNP", "HOUSE", "THFILE", "WAVE"]
    assert s.delete_entry("WAVE", 2) and not s.delete_entry("WAVE", 2)
    s.set_string("THFILE", "")
    assert s.string("THFILE", "none") == "none"
