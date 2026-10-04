"""Frame files (D-FIL-03), PROCFRAME frame store and SASSIani.xml (requirements 5.8, D-UI-15), animation
commands (BUBBLEPLOT, VECTORPLOT, CONTOURPLOT, DEFORMPLOT, PAUSE) and the frame / frequency tools CRITFREQ,
FRAMECOMBIN, FRAMESEL, MODFRAMES (spec 09 sections 2.4, 2.11, 2.12, 2.24)."""
from __future__ import annotations

import json

import numpy as np
import pytest

from sassi.plotting import state as S
from sassi.plotting.state import plot_state, write_frame, read_frame
from sassi.prep import Interpreter, Kind


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


@pytest.fixture
def ui(tmp_path, monkeypatch):
    monkeypatch.setenv("SASSI_EDU_SETTINGS_DIR", str(tmp_path / "settings"))
    u = Interpreter(cwd=tmp_path)
    plot_state(u).auto_render = False
    return u


def box_model(ui, n=2):
    """(n+1)^3 nodes, n^3 hexes of unit size."""
    lines = []
    nid = lambda i, j, k: 1 + i + (n + 1) * j + (n + 1) ** 2 * k
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                lines.append(f"N,{nid(i, j, k)},{i},{j},{k}")
    lines.append("GROUP,1,SOLID")
    e = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                c = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k)]
                c += [v + (n + 1) ** 2 for v in c]
                lines.append(f"E,{e}," + ",".join(map(str, c)))
                e += 1
    ui.run_text("\n".join(lines))
    return (n + 1) ** 3


def make_frames(tmp_path, nnodes, nframes=6, sub="ACC", cols=3, start_node=1):
    d = tmp_path / sub
    d.mkdir(exist_ok=True)
    nodes = np.arange(start_node, start_node + nnodes)
    names = []
    for s in range(nframes):
        vals = np.column_stack([(s + 1) * 0.01 * nodes + c for c in range(cols)])
        nm = f"{sub}_{s * 0.005:06.3f}_{s + 1:05d}" if sub != "TFU" else f"TFU_{0.5 * (s + 1):06.2f}_{s + 1:05d}"
        write_frame(d / nm, nodes, vals)
        names.append(nm)
    lst = tmp_path / f"{sub.lower()}.dispani"
    lst.write_text("PREP frame list\n" + "\n".join(f"{sub}/{nm}" for nm in names) + "\n")
    return lst, nodes


# ------------------------------------------------------------------ frame files
def test_frame_file_layout_and_errors(tmp_path):
    p = write_frame(tmp_path / "f1", [5, 7], [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    assert p.read_text().splitlines()[0] == "2 4"            # header nrows ncols, node id = column 1
    fr = read_frame(p)
    assert fr.nodes.tolist() == [5, 7] and fr.values.tolist() == [[1, 2, 3], [4, 5, 6]] and fr.ncols == 4
    (tmp_path / "bad").write_text("2 3\n1 1.0 2.0 3.0\n2 1.0 2.0 3.0\n")
    with pytest.raises(S.PlotError, match="MODFRAMES"):
        read_frame(tmp_path / "bad")
    (tmp_path / "legacy").write_text("1.0 2.0\n")
    with pytest.raises(S.PlotError):
        read_frame(tmp_path / "legacy")
    text, ch = S.modframes_text("2 3\n1 1 2 3\n", 4)
    assert ch and text.startswith("2 4\n")
    text, ch = S.modframes_text("2\n1 1 2 3\n", 4)
    assert ch and text.startswith("2 4\n")
    assert S.frame_name_info("ACC_00.015_00004") == {"tag": "ACC", "value": 0.015, "frame": 4, "comp": None}
    assert S.frame_name_info("stress_00.000_00001_sig")["comp"] == "sig"
    assert S.frame_name_info("RS01_000.10_00001")["tag"] == "RS01"


def test_combine_frames_ops(tmp_path):
    a = S.Frame(np.array([1, 2]), np.array([[3.0, 1.0], [0.0, -2.0]]))
    b = S.Frame(np.array([1, 2]), np.array([[4.0, 1.0], [1.0, 2.0]]))
    assert S.combine_frames([a, b], 0).values.tolist() == [[5.0, np.sqrt(2)], [1.0, np.sqrt(8)]]
    assert S.combine_frames([a, b], 1).values.tolist() == [[7.0, 2.0], [1.0, 0.0]]
    assert S.combine_frames([a, b], 2).values.tolist() == [[3.5, 1.0], [0.5, 0.0]]
    with pytest.raises(S.PlotError):
        S.combine_frames([a, S.Frame(np.array([1, 3]), b.values)], 1)
    with pytest.raises(S.PlotError):
        S.combine_frames([a], 3)


def test_frame_list_relative_paths_and_directories(tmp_path):
    lst, _ = make_frames(tmp_path, 4, nframes=3)
    files = S.read_frame_list(lst)
    assert [f.name for f in files] == ["ACC_00.000_00001", "ACC_00.005_00002", "ACC_00.010_00003"]
    assert all(f.parent == tmp_path / "ACC" for f in files)
    assert [f.name for f in S.read_frame_list(tmp_path / "ACC")] == [f.name for f in files]
    (tmp_path / "empty.dispani").write_text("header only\n")
    with pytest.raises(S.PlotError):
        S.read_frame_list(tmp_path / "empty.dispani")


# ------------------------------------------------------------------ PROCFRAME and the frame store
def test_procframe_store_and_animation_database(ui, tmp_path):
    lst, nodes = make_frames(tmp_path, 27, nframes=6)
    events = []
    plot_state(ui).subscribe(events.append)
    assert ui.execute("PROCFRAME,acc.dispani,buf,Accelerations, X Y Z,0")
    idx = json.loads((tmp_path / "buf" / "index.json").read_text())
    assert idx["nframes"] == 6 and idx["data_columns"] == 3 and idx["layout"] == "xyz"
    assert idx["description"] == "Accelerations, X Y Z" and idx["shared_nodes"]
    assert (tmp_path / "buf" / "frame_00006.npy").exists() and (tmp_path / "buf" / "nodes.npy").exists()
    assert idx["frames"][2]["value"] == pytest.approx(0.010) and idx["frames"][2]["tag"] == "ACC"
    store = S.FrameStore(tmp_path / "buf")
    assert np.allclose(store.values(2)[:, 0], 0.02 * nodes, rtol=1e-6)          # float32 store
    assert store.column_range([1, 6], 1) == pytest.approx((0.01, 0.06 * 27), rel=1e-6)
    assert "ACC 0.005 s" in store.frame_label(2)
    db = S.AnimationDB(tmp_path / "settings" / "SASSIani.xml")
    ents = db.entries()
    assert len(ents) == 1 and ents[0]["frames"] == "6" and ents[0]["type"] == "0"
    assert ents[0]["description"] == "Accelerations, X Y Z"
    assert [e.kind for e in events if e.kind in ("progress", "animations")][-1] == "animations"
    assert any(e.kind == "progress" and e.data["fraction"] == 1.0 for e in events)
    assert ui.execute("PROCFRAME,acc.dispani,buf,again")        # same directory: entry replaced
    assert len(db.entries()) == 1 and db.entries()[0]["description"] == "again"
    assert db.remove(tmp_path / "buf", delete_files=True) and not (tmp_path / "buf").exists()
    assert not ui.execute("PROCFRAME")                           # Parse Frame Data dialog: batch error
    assert not ui.execute("PROCFRAME,acc.dispani")


def test_procframe_varying_node_lists(ui, tmp_path):
    d = tmp_path / "THD"
    d.mkdir()
    write_frame(d / "THD_00.000_00001", [1, 2], [[1, 0, 0], [2, 0, 0]])
    write_frame(d / "THD_00.010_00002", [1, 2, 3], [[1, 0, 0], [2, 0, 0], [3, 0, 0]])
    (tmp_path / "thd.dispani").write_text("x\nTHD/THD_00.000_00001\nTHD/THD_00.010_00002\n")
    assert ui.execute("PROCFRAME,thd.dispani,bufthd")
    st = S.FrameStore(tmp_path / "bufthd")
    assert not st.index["shared_nodes"] and st.nodes(1).tolist() == [1, 2] and st.nodes(2).tolist() == [1, 2, 3]


# ------------------------------------------------------------------ animations
def test_bubble_and_contour_requests(ui, tmp_path):
    nn = box_model(ui)
    lst, nodes = make_frames(tmp_path, nn, nframes=6)
    ui.execute("PROCFRAME,acc.dispani,buf")
    st = plot_state(ui)
    assert ui.execute("BUBBLEPLOT,buf,2,6,2,0,0.5,1")
    p = st.active_plot
    assert p.params["frames"] == [2, 4, 6] and p.params["current"] == 2
    assert (p.params["vmin"], p.params["vmax"], p.params["col"]) == (0.0, 0.5, 1)
    d = st.plot_data(p.id, ui)
    fr = d["frame"]
    v = np.asarray(fr["value"], float)
    assert np.allclose(v, 0.02 * np.asarray(d["scene"]["node_id"]), rtol=1e-6)
    t = np.asarray(fr["t"], float)
    assert np.allclose(t, np.clip(v / 0.5, 0, 1), rtol=1e-6)
    assert np.allclose(fr["size"], 10 * np.maximum(0.2, t))
    assert len(fr["colorbar"]) == 5
    assert ui.execute("CONTOURPLOT,buf,1,6,1")                  # partial arguments: defaults + warning
    assert any("partial arguments" in w for w in warnings(ui))
    q = st.active_plot
    assert q.params["vmin"] == pytest.approx(0.01) and q.params["vmax"] == pytest.approx(0.06 * 27, rel=1e-6)
    assert not ui.execute("BUBBLEPLOT,buf,1,6,1,0,1,4")          # Col 4 > 3 data columns
    assert not ui.execute("BUBBLEPLOT,buf,1,7,1,0,1,1")          # MxF > frames
    assert not ui.execute("BUBBLEPLOT,buf,1,6,1,2,1,1")          # MnR > MxR
    assert not ui.execute("BUBBLEPLOT,nodir,1,6,1,0,1,1")
    assert not ui.execute("BUBBLEPLOT")
    json.dumps(d)


def test_vector_deform_and_pause(ui, tmp_path):
    nn = box_model(ui)
    make_frames(tmp_path, nn, nframes=4, sub="TFU", cols=6)
    lst, nodes = make_frames(tmp_path, nn, nframes=4, sub="THD", cols=3)
    (tmp_path / "tfu.dispani").write_text("h\n" + "\n".join(f"TFU/{f.name}" for f in sorted((tmp_path / "TFU").iterdir())))
    ui.execute("PROCFRAME,tfu.dispani,buft,TF,1")
    ui.execute("PROCFRAME,thd.dispani,bufd,disp,0")
    st = plot_state(ui)
    assert ui.execute("VECTORPLOT,buft,1,4,1,2")
    p = st.active_plot
    assert p.params["layout"] == "complex_xyz" and p.params["scale"] == 2.0 and p.view.direction == "X"
    d = st.plot_data(p.id, ui, raw=True)
    k = 5                                                       # 6th scene node (node id 6)
    nid = int(d["scene"]["node_id"][k])
    # frame 1 columns: ReX = 0.01 n, ImX = 0.01 n + 1, ReY = 0.01 n + 2, ImY ..., ReZ ..., ImZ = 0.01 n + 5
    a = 0.01 * nid
    assert np.allclose(d["frame"]["vectors"][k, 0], 2 * np.array([a, a + 1, a + 1]), rtol=1e-6)
    assert np.allclose(d["frame"]["vectors"][k, 2], 2 * np.array([a + 5, a + 5, a + 4]), rtol=1e-6)
    assert ui.execute("WINDOWSETTINGS,DIRECTION,ALL") and p.view.direction == "ALL"
    ui.execute("PAUSE")
    assert p.view.paused
    ui.execute("PAUSE,0")
    assert not p.view.paused
    assert not ui.execute("PAUSE,3")
    assert ui.execute("DEFORMPLOT,bufd,1,4,1,10")
    q = st.active_plot
    d = st.plot_data(q.id, ui, raw=True)
    xyz0 = d["scene"]["xyz"]
    u = 0.01 * d["scene"]["node_id"]
    assert np.allclose(d["frame"]["xyz"][:, 0], xyz0[:, 0] + 10 * u, rtol=1e-6)
    assert ui.execute("WINDOWSETTINGS,FRAME,3") and q.params["current"] == 3
    d = st.plot_data(q.id, ui, raw=True)
    assert np.allclose(d["frame"]["xyz"][:, 0], xyz0[:, 0] + 30 * u, rtol=1e-6)
    assert not ui.execute("WINDOWSETTINGS,FRAME,9")
    assert ui.execute("WINDOWSETTINGS,SCALE,1") and q.params["scale"] == 1.0
    assert ui.execute("SHADEROPTIONS,,,,4") and q.params["scale"] == 4.0
    assert ui.execute("WINDOWSETTINGS,UNDEFORMED,1") and q.view.show_undeformed


def test_animation_frames_with_unknown_nodes_and_rendering(ui, tmp_path):
    from sassi.plotting.render_mpl import render_plot
    nn = box_model(ui)
    make_frames(tmp_path, nn + 3, nframes=2)                     # 3 nodes not in the model
    ui.execute("PROCFRAME,acc.dispani,buf")
    st = plot_state(ui)
    for cmd in ("BUBBLEPLOT,buf,1,2,1,0,1,2", "CONTOURPLOT,buf,1,2,1,0,1,3", "VECTORPLOT,buf,1,2,1,1",
                "DEFORMPLOT,buf,1,2,1,1"):
        assert ui.execute(cmd), errors(ui)
        out, notes = render_plot(st, st.active_plot, ui, tmp_path / (cmd.split(",")[0] + ".png"))
        assert out.stat().st_size > 2000
    assert any("3 of 30 frame nodes" in w for w in warnings(ui))
    ui.execute("ACTM,2")
    assert not ui.execute("BUBBLEPLOT,buf,1,2,1,0,1,1")          # model without elements


# ------------------------------------------------------------------ frame / frequency tools
def test_framecombin_command(ui, tmp_path):
    write_frame(tmp_path / "a", [1, 2], [[3.0, 0.0], [1.0, 1.0]])
    write_frame(tmp_path / "b", [1, 2], [[4.0, 2.0], [1.0, 1.0]])
    assert ui.execute("FRAMECOMBIN,0,2,a,b,srss")
    assert np.allclose(read_frame(tmp_path / "srss").values, [[5.0, 2.0], [np.sqrt(2), np.sqrt(2)]], rtol=1e-10)
    assert ui.execute("FRAMECOMBIN,2,2,a,b,avg")
    assert read_frame(tmp_path / "avg").values.tolist() == [[3.5, 1.0], [1.0, 1.0]]
    assert not ui.execute("FRAMECOMBIN,1,3,a,b,out")             # num does not match the file count
    assert not ui.execute("FRAMECOMBIN,5,2,a,b,out")
    write_frame(tmp_path / "c", [1, 3], [[3.0, 0.0], [1.0, 1.0]])
    assert not ui.execute("FRAMECOMBIN,1,2,a,c,out")


def test_framesel_command(ui, tmp_path):
    a = [0.0, 1.0, 0.0, -0.5, 0.0, 0.2, 0.0, -2.0, 0.0]
    (tmp_path / "n.acc").write_text("0.01\n" + "\n".join(str(v) for v in a))
    assert ui.execute("FRAMESEL,40,n.acc,CF")
    assert ui.variables["cf"].items == ["2", "8"] and ui.variables["cf"].counter == 0
    (tmp_path / "p.th").write_text("\n".join(f"{0.01 * k} {v}" for k, v in enumerate(a)))
    assert ui.execute("FRAMESEL,40,p.th,CF2") and ui.variables["cf2"].items == ["2", "8"]
    ui.execute("VAR,X,1")
    assert ui.execute("N,@CF[2],1,2,3") and 8 in ui.model.nodes   # usable with @Var[i]
    assert not ui.execute("FRAMESEL,120,n.acc,CF") and not ui.execute("FRAMESEL,10,n.acc,9bad")


def test_modframes_command(ui, tmp_path):
    d = tmp_path / "F"
    d.mkdir()
    (d / "f1").write_text("2 3\n1 1 2 3\n2 1 2 3\n")
    (d / "f2").write_text("2\n1 1 2 3\n2 1 2 3\n")
    (tmp_path / "f.zpani").write_text("legacy\nF/f1\nF/f2\n")
    assert ui.execute("MODFRAMES,4,f.zpani")
    assert read_frame(d / "f1").ncols == 4 and read_frame(d / "f2").ncols == 4
    assert not ui.execute("MODFRAMES,1,f.zpani")


def test_critfreq_command_writes_frequency_numbers(ui, tmp_path):
    F = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    AU = np.array([1.0, 1.2, 2.0, 2.5, 1.5, 1.0])
    fi = np.round(np.arange(1.0, 6.0 + 1e-9, 0.05), 10)
    ai = np.interp(fi, F, AU) + 1.25 * np.clip(1.0 - np.abs(fi - 3.5) / 0.5, 0.0, None)
    (tmp_path / "00010TR_X.TFU").write_text("\n".join(f"{a} {b} 0.1" for a, b in zip(F, AU)))
    (tmp_path / "00010TR_X.tfi").write_text("\n".join(f"{a} {b} 0.1" for a, b in zip(fi, ai)))
    assert ui.execute("CRITFREQ,20,50,00010TR_X,FR")
    assert ui.variables["fr"].items == ["70"]                    # 3.5 Hz / 0.05 Hz
    assert ui.execute("CRITFREQ,50,50,00010TR_X.TFU,FR")         # extension given: accepted
    assert ui.variables["fr"].items == []
    assert any("3.5 Hz" in t for t in ui.sink.texts(Kind.INFO))
    assert not ui.execute("CRITFREQ,20,50,nothere,FR")
    assert not ui.execute("CRITFREQ,20,150,00010TR_X,FR")
