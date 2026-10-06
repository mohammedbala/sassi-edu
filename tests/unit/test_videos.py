"""Explainer videos of the guided course (sassi/ui/static/videos, docs/internal/explainer_videos.md): one per
lesson, their pages, data files and recorded narration (ElevenLabs clips, web/voice_videos.mjs), the GUI routes
and the Learn buttons, the data tool (sassi.ui.video_data), the narration (short, no prose of the ACS SASSI
manual, no API key anywhere) and -- with Node.js and Chrome -- every beat of every video in headless Chrome
(web/test_videos.mjs)."""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import threading
import urllib.request
from pathlib import Path

import numpy as np
import pytest

from sassi.ui import lessons as L
from sassi.ui import video_data

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"
VIDEOS = STATIC / "videos"
NODE = shutil.which("node")
CHROME = os.environ.get("CHROME", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
LESSONS = L.load_lessons()


def _nn(i: int) -> str:
    return f"{i:02d}"


def _says(js: str):
    """The narration of a video script: the `say:` strings (double or single quotes, template literals)."""
    out = []
    for m in re.finditer(r"""say:\s*(?:"((?:[^"\\]|\\.)*)"|'((?:[^'\\]|\\.)*)'|`((?:[^`\\]|\\.)*)`)""", js):
        out.append(next(g for g in m.groups() if g is not None))
    return out


def test_every_lesson_has_a_video():
    assert len(LESSONS) == 11
    for i, les in enumerate(LESSONS, start=1):
        nn = _nn(i)
        page, script, data = VIDEOS / f"{nn}.html", VIDEOS / f"v{nn}.js", VIDEOS / f"d{nn}.js"
        assert page.is_file() and script.is_file() and data.is_file(), nn
        html = page.read_text(encoding="utf-8")
        refs = re.findall(r'(?:src|href)="([^"]+)"', html)
        assert refs == ["../katex/katex.min.css", "player.css", "../katex/katex.min.js", "../figures.js", "player.js",
                        "kit.js", f"d{nn}.js", f"a{nn}.js", f"v{nn}.js"], (nn, refs)
        for r in refs:
            assert (VIDEOS / r).resolve().is_file(), (nn, r)
        js = script.read_text(encoding="utf-8")
        assert f'lesson: "{les.id}"' in js and re.search(rf"\bn: {i}\b", js), nn
        assert f'SV.DATA["{nn}"]' in js, nn
        if i < len(LESSONS):
            assert f'href: "{_nn(i + 1)}.html"' in js, nn                     # the end card links the next one
        else:
            assert "next:" not in js, nn
        says = _says(js)
        assert len(says) >= 25, (nn, len(says))
        words = sum(len(re.sub(r"\[\w+\]|\{([^{}|]*)\|[^{}]*\}", r"\1", s).split()) for s in says)
        assert 400 <= words <= 850, (nn, words)                                # about 3-4 minutes at 1.15x
        assert max(len(re.sub(r"\[\w+\]", "", s).split()) for s in says) <= 32, nn  # short beats


def test_video_data_files_are_generated():
    for i, les in enumerate(LESSONS, start=1):
        nn = _nn(i)
        text = (VIDEOS / f"d{nn}.js").read_text(encoding="utf-8")
        assert text.startswith(f"/* Curves of lesson {les.id} for explainer video {nn} (generated"), nn
        assert f"python -m sassi.ui.video_data {les.id}" in text and "--from" not in text and "--keep" not in text
        assert f'SV.DATA["{nn}"] = {{' in text, nn
        assert len(text.encode()) < 90_000, (nn, len(text))


def test_video_pages_use_relative_urls_only():
    for f in sorted(VIDEOS.glob("*.html")) + sorted(VIDEOS.glob("*.js")) + [VIDEOS / "player.css"]:
        text = f.read_text(encoding="utf-8")
        assert not re.search(r"""(?:src|href)=["'](?:/|https?:)""", text), f.name
        assert not re.findall(r"https?://(?!www\.w3\.org/2000/svg)", text), f.name           # no external requests
        assert "setTimeout(" not in text or f.name in ("player.js", "index.html"), f.name   # scenes use the clock


def test_narration_is_recorded():
    """Every video has its manifest aNN.js; every clip it names exists, nothing else is in audio/NN; the clips
    are small MP3s.  (That every beat's current text has its clip is checked with Node below.)"""
    import json
    total = 0
    for i, les in enumerate(LESSONS, start=1):
        nn = _nn(i)
        text = (VIDEOS / f"a{nn}.js").read_text(encoding="utf-8")
        assert text.startswith(f"/* ElevenLabs narration of explainer video {nn} (generated"), nn
        man = json.loads(text[text.index(f'SV.AUDIO["{nn}"] = ') + len(f'SV.AUDIO["{nn}"] = '):text.rindex("};") + 1])
        assert man["voice"] and man["model"] and man["beats"], nn
        files = set()
        for h, clip in man["beats"].items():
            assert re.fullmatch(r"[0-9a-f]{8}", h) and clip["f"] == f"audio/{nn}/{h}.mp3" and 0.3 < clip["d"] < 30, (nn, h)
            f = VIDEOS / clip["f"]
            assert f.is_file() and f.read_bytes()[:3] in (b"ID3", b"\xff\xfb", b"\xff\xf3"), (nn, h)
            files.add(f.name)
            total += f.stat().st_size
        assert {f.name for f in (VIDEOS / "audio" / nn).iterdir()} == files, nn          # no orphans
    assert total < 45_000_000, total


@pytest.mark.skipif(NODE is None, reason="Node.js is not installed")
def test_every_beat_has_its_recording():
    r = subprocess.run([NODE, str(ROOT / "web" / "voice_videos.mjs"), "--dry"], cwd=str(ROOT), capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr
    lines = re.findall(r"^(\d\d): (\d+) beats, (\d+) to record", r.stdout, re.M)
    assert len(lines) == len(LESSONS) and all(todo == "0" for _, _, todo in lines), r.stdout


def test_no_api_key_in_the_project():
    """The ElevenLabs key is given to web/voice_videos.mjs through the environment, never stored."""
    for base in (VIDEOS, ROOT / "web", ROOT / "docs", ROOT / "sassi" / "ui"):
        for f in base.rglob("*"):
            if f.is_file() and f.suffix in (".js", ".mjs", ".py", ".md", ".html", ".json", ".css") and "node_modules" not in f.parts:
                assert not re.search(r"\bsk_[0-9a-f]{32,}", f.read_text(encoding="utf-8", errors="replace")), f


def test_index_lists_every_lesson():
    html = (VIDEOS / "index.html").read_text(encoding="utf-8")
    for i, les in enumerate(LESSONS, start=1):
        assert les.title in html, les.id
    assert 'href="01.html"' in html and "SV.catalog = []" in html


def test_learn_buttons_and_menu_open_the_videos():
    learn = (STATIC / "learn.js").read_text(encoding="utf-8")
    assert 'L.videoUrl = (id) => (id ? `static/videos/${String(L.lessonNumber(id)).padStart(2, "0")}.html` : "static/videos/index.html");' in learn
    assert 'btn("▶ Explainer", () => L.watchVideo(les.id)' in learn
    assert 'btn("▶ Watch the explainer video", () => L.watchVideo(les.id)' in learn
    assert '{label: "Explainer Videos", action: () => L.watchVideo(null)' in learn


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
        with urllib.request.urlopen(srv.url.rstrip("/") + path, timeout=10) as r:
            return r.status, r.read(), r.headers.get("Content-Type", "")
    except urllib.error.HTTPError as e:
        return e.code, e.read(), ""


def test_server_serves_the_videos(server):
    for name, ctype in (("01.html", "text/html"), ("player.js", "javascript"), ("player.css", "text/css"),
                        ("d01.js", "javascript"), ("index.html", "text/html")):
        st, body, ct = _get(server, "/static/videos/" + name)
        assert st == 200 and ctype in ct and body == (VIDEOS / name).read_bytes(), name
    for bad in ("/static/videos/missing.js", "/static/videos/../app.js", "/static/videos/.x.js",
                "/static/videos/sub/x.js", "/static/videos/", "/static/videos/x.py"):
        assert _get(server, bad)[0] == 404, bad
    assert _get(server, "/static/figures.js")[0] == 200 and _get(server, "/static/katex/katex.min.js")[0] == 200
    clip = next((VIDEOS / "audio" / "01").glob("*.mp3"))
    st, body, ct = _get(server, f"/static/videos/audio/01/{clip.name}")
    assert st == 200 and ct == "audio/mpeg" and body == clip.read_bytes()
    for bad in ("/static/videos/audio/01/x.mp3", "/static/videos/audio/1/0000000a.mp3", "/static/videos/audio/01/../a01.js"):
        assert _get(server, bad)[0] == 404, bad


def test_video_data_extract(tmp_path):
    f = np.linspace(0, 20, 821)
    H = 1 + 12 * np.exp(-((f - 3.49) / 0.2) ** 2)
    tfi = tmp_path / "00085TR_X.TFI"
    tfi.write_text("# node 85\n" + "".join(f"{a:.6f} {b:.8e} {-0.1:.8e}\n" for a, b in zip(f, H)))
    d = video_data.extract(tfi, points=100)
    assert len(d["f"]) <= 103 and d["f"][0] == 0 and d["f"][-1] == 20
    assert max(d["amp"]) == pytest.approx(H.max(), rel=1e-3)                # the peak survives the thinning
    rs = tmp_path / "x.RS"
    rs.write_text("# spectrum\n" + "".join(f"{a:.6f} {b:.6e}\n" for a, b in zip(f[1:], H[1:])))
    assert set(video_data.extract(rs, 50)) == {"f", "sa"}
    acc = tmp_path / "x.ACC"
    v = np.sin(np.arange(4800) * 0.05) * np.exp(-np.arange(4800) / 2000)
    v[1234] = 3.0
    acc.write_text("0.005\n" + "".join(f"{x:.6e}\n" for x in v))
    h = video_data.extract(acc, 200)
    assert h["dt"] == 0.005 and h["n"] == 4800 and h["peak"] == 3.0 and 3.0 in h["v"] and len(h["t"]) <= 200


@pytest.mark.skipif(NODE is None, reason="Node.js is not installed")
def test_video_js_syntax():
    files = [VIDEOS / "player.js", VIDEOS / "kit.js", ROOT / "web" / "test_videos.mjs"]
    files += sorted(VIDEOS.glob("v??.js")) + sorted(VIDEOS.glob("d??.js"))
    for f in files:
        r = subprocess.run([NODE, "--check", str(f)], capture_output=True, text=True, timeout=60)
        assert r.returncode == 0, (f.name, r.stderr)


def test_narration_has_no_prose_of_the_manual():
    """The narration paraphrases the lessons; it never copies the ACS SASSI manual (reference/acs-sassi.txt is
    the user's local, copyrighted copy: the check runs only where it exists)."""
    manual = ROOT / "reference" / "acs-sassi.txt"
    if not manual.is_file():
        pytest.skip("reference/acs-sassi.txt is not present")
    norm = lambda t: re.findall(r"[a-z0-9]+", t.lower())
    words = norm(manual.read_text(encoding="utf-8", errors="replace"))
    n = 10
    shingles = {" ".join(words[i:i + n]) for i in range(len(words) - n + 1)}
    for f in sorted(VIDEOS.glob("v??.js")):
        for s in _says(f.read_text(encoding="utf-8")):
            w = norm(re.sub(r"\[\w+\]|\{([^{}|]*)\|[^{}]*\}", r"\1", s))
            hits = [" ".join(w[i:i + n]) for i in range(len(w) - n + 1) if " ".join(w[i:i + n]) in shingles]
            assert not hits, (f.name, hits[0])


@pytest.mark.slow
@pytest.mark.skipif(NODE is None or not Path(CHROME).is_file(), reason="needs Node.js and Google Chrome")
def test_every_beat_of_every_video_in_chrome():
    r = subprocess.run([NODE, str(ROOT / "web" / "test_videos.mjs")], cwd=str(ROOT), capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
    assert len(re.findall(r"^\d\d: \d+ scenes, \d+ beats, ~[\d.]+ min  ok$", r.stdout, re.M)) == len(LESSONS), r.stdout
