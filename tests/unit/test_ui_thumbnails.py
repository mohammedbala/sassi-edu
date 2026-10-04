"""Preview pictures of the examples gallery (sassi.ui.thumbnails, sassi/ui/static/examples/*.png).

Every example has its picture and every picture an example; the pictures are small PNGs of the expected
size; the generator builds an example's model part without module runs in a temporary folder that it
deletes, draws it on the Agg canvas without pyplot, and writes the same bytes every time; the local
server serves the pictures (the web build ships them: tests/unit/test_web.py)."""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import threading
import urllib.error
import urllib.request
from pathlib import Path

import pytest

from sassi.prep import Interpreter
from sassi.plotting.state import plot_data, plot_state
from sassi.ui import learn
from sassi.ui import thumbnails as T
from sassi.ui.lessons import EXAMPLES_DIR

ROOT = Path(__file__).resolve().parents[2]
PICTURES = ROOT / "sassi" / "ui" / "static" / "examples"
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _examples():
    return sorted(p.stem for p in EXAMPLES_DIR.glob("*.pre"))


# ====================================================================== the committed pictures
def test_every_example_has_a_picture_and_every_picture_an_example():
    pictures = sorted(p.stem for p in PICTURES.glob("*.png"))
    assert _examples() and pictures == _examples()
    assert T.example_names() == _examples()
    assert not [p.name for p in PICTURES.iterdir() if p.suffix != ".png"]           # nothing else in the folder


def test_pictures_are_valid_small_pngs_of_the_expected_size():
    from PIL import Image
    total = 0
    for p in sorted(PICTURES.glob("*.png")):
        data = p.read_bytes()
        assert data[:8] == PNG_SIGNATURE, p.name
        with Image.open(p) as im:
            im.verify()                                     # a complete, uncorrupted PNG
        with Image.open(p) as im:
            assert im.format == "PNG" and im.size == (T.WIDTH, T.HEIGHT) == (480, 300), p.name
        assert len(data) < 60_000, (p.name, len(data))
        total += len(data)
    assert total < 1_000_000


# ====================================================================== the model part of an example
def test_model_lines_leave_out_checks_writes_and_module_runs():
    for name in _examples():
        text = (EXAMPLES_DIR / f"{name}.pre").read_text(encoding="utf-8")
        lines = learn.example_model_lines(text)
        names = [learn._command_name(ln) for ln in lines]
        assert lines and names[0] == "MDL", name
        assert not [n for n in names if n.startswith("RUN") or n in learn.NOT_MODEL or n in learn.MODEL_SWITCH], name
        assert not [ln for ln in lines if ln.upper().startswith("FOREACH") and learn._foreach_command(ln).startswith("@")]
        assert not [ln for ln in lines if ln.startswith("*")]
        assert learn.example_view_command(lines) == ("LAYERPLOT" if name == "ex04_site_response" else "MODELPLOT")
    ex2 = learn.example_model_lines((EXAMPLES_DIR / "ex02_embedded_box.pre").read_text(encoding="utf-8"))
    assert ex2[0] == "MDL,ex02,ex02" and not [ln for ln in ex2 if "fsin" in ln.lower()]    # model 1 only
    ex1 = learn.example_model_lines((EXAMPLES_DIR / "ex01_surface_stick.pre").read_text(encoding="utf-8"))
    assert "FOREACH,SPIDER,E,#,41,@SPIDER[#],90" in ex1                     # element generation is kept
    assert learn.example_model_lines("* c\nN,1\nrunsite\nafwr\nCHECK\nE,1,1,2\nACTM,2\nN,3\n") == ["N,1", "E,1,1,2"]


def test_every_example_model_builds_without_error_and_matches_its_view():
    """The view chosen from the command text (MODELPLOT / LAYERPLOT) fits the model the commands build."""
    for name in _examples():
        d = T.example_plot(name)
        if name == "ex04_site_response":
            assert d["kind"] == "LAYERPLOT" and len(d["table"]["layers"]) == 22 and d["table"]["halfspace"]
        else:
            assert d["kind"] == "MODELPLOT" and d["scene"]["n_elements"] > 0, name
            assert len(d["scene"]["interaction"]) > 0, name                  # every structure has interaction nodes


def test_shown_faces_keep_coincident_solids_and_hide_shared_faces(tmp_path):
    """Two solids on the same nodes (excavated soil and backfill) show the volume once (the later element);
    two solids sharing a face hide that face.  The scene's rule (GUI 3D view) and the pictures agree."""
    ui = Interpreter(cwd=tmp_path)
    plot_state(ui).auto_render = False
    for ln in ("N,1,0,0,0", "N,2,1,0,0", "N,3,1,1,0", "N,4,0,1,0", "NGEN,2,4,1,4,1,0,0,1",
               "GROUP,1,SOLID", "E,1,1,2,3,4,5,6,7,8", "GROUP,2,SOLID", "E,1,1,2,3,4,5,6,7,8", "MODELPLOT"):
        assert ui.execute(ln), ln
    sc = plot_data(plot_state(ui), plot_state(ui).active_plot, ui, raw=True)["scene"]
    assert int(sc["face_boundary"].sum()) == 6                                # the volume once, as the GUI
    shown = T.shown_faces(sc)
    assert len(shown) == 6 and {int(sc["elem_group"][sc["face_elem"][r]]) for r in shown} == {2}
    for ln in ("NGEN,1,4,5,8,1,0,0,1", "GROUP,3,SOLID", "E,1,5,6,7,8,9,10,11,12", "MODELPLOT"):
        assert ui.execute(ln), ln
    sc = plot_data(plot_state(ui), plot_state(ui).active_plot, ui, raw=True)["scene"]
    assert len(T.shown_faces(sc)) == 10                                      # 6 + 6 - the shared face twice


# ====================================================================== the generator
def _thumb_dirs():
    return sorted(p.name for p in Path(tempfile.gettempdir()).glob("sassi-thumb-*"))


def test_generator_writes_one_example_into_a_folder_and_cleans_up(tmp_path):
    from PIL import Image
    before = _thumb_dirs()
    out = tmp_path / "pics"
    sizes = T.generate(["ex03_forced_vibration"], out=out)
    assert sizes == {"ex03_forced_vibration.png": (out / "ex03_forced_vibration.png").stat().st_size}
    assert sorted(p.name for p in out.iterdir()) == ["ex03_forced_vibration.png"]
    with Image.open(out / "ex03_forced_vibration.png") as im:
        assert im.size == (480, 300)
    assert T.thumbnail("ex03_forced_vibration") == (out / "ex03_forced_vibration.png").read_bytes()   # deterministic
    assert _thumb_dirs() == before                                             # the temporary copy is deleted
    assert not [p for p in EXAMPLES_DIR.iterdir() if p.is_dir() and p.name.startswith("ex03")]       # nothing in examples/


def test_generator_command_line(tmp_path, capsys):
    assert T.main(["ex04_site_response", "--out", str(tmp_path)]) == 0
    assert "ex04_site_response.png" in capsys.readouterr().out and (tmp_path / "ex04_site_response.png").is_file()
    assert T.main(["nosuch", "--out", str(tmp_path)]) == 2


def test_rendering_uses_the_agg_canvas_without_pyplot():
    code = ("import sys\n"
            "from sassi.ui import thumbnails as T\n"
            "fig = T.figure(T.example_plot('ex01_surface_stick'))\n"
            "print(type(fig.canvas).__name__, 'matplotlib.pyplot' in sys.modules)\n")
    env = dict(os.environ, MPLBACKEND="macosx", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=ROOT, env=env, timeout=120)
    assert r.returncode == 0, r.stderr
    assert r.stdout.split() == ["FigureCanvasAgg", "False"]


# ====================================================================== the local server
@pytest.fixture
def server(tmp_path):
    from sassi.ui.server import make_server
    srv = make_server(0, model_dir=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    th = threading.Thread(target=srv.serve_forever, kwargs={"poll_interval": 0.05}, daemon=True)
    th.start()
    yield srv
    srv.shutdown()
    srv.server_close()


def _get(srv, path):
    try:
        with urllib.request.urlopen(srv.url.rstrip("/") + path, timeout=30) as r:
            return r.status, r.read(), r.headers
    except urllib.error.HTTPError as e:
        return e.code, e.read(), e.headers


def test_server_serves_the_pictures(server):
    for name in _examples():
        st, body, hd = _get(server, f"/{learn.THUMB_URL}/{name}.png")
        assert st == 200 and hd["Content-Type"] == "image/png" and body[:8] == PNG_SIGNATURE, name
    for bad in ("/static/examples/missing.png", "/static/examples/../app.js", "/static/examples/.x.png",
                "/static/examples/ex01_surface_stick.pre", "/static/examples/sub/x.png", "/static/examples/"):
        assert _get(server, bad)[0] == 404, bad
