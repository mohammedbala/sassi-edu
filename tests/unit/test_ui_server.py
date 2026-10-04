"""The ``sassi-gui`` HTTP server (ARCHITECTURE section 9): in-process server driven with urllib.

Static front end (page with the session token, JavaScript, CSS, Plotly from the installed
package), the JSON API over HTTP (command round trip, options, model, files, a module run), and the
safety rules: loopback only, session token, Host header check, no path traversal in static files or
the file API."""
from __future__ import annotations

import json
import os
import re
import signal
import socket
import stat
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

from sassi.ui import server as srvmod
from sassi.ui.server import make_server


@pytest.fixture
def server(tmp_path):
    srv = make_server(0, model_dir=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    th = threading.Thread(target=srv.serve_forever, kwargs={"poll_interval": 0.05}, daemon=True)
    th.start()
    yield srv
    srv.shutdown()
    srv.server_close()


def http(srv, method, path, body=None, token=True, headers=None, raw=False):
    h = {"Content-Type": "application/json"}
    if token:
        h["X-SASSI-Token"] = srv.session.token
    h.update(headers or {})
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(srv.url.rstrip("/") + path, data=data, method=method, headers=h)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            content = r.read()
            return r.status, (content if raw else json.loads(content)), r.headers
    except urllib.error.HTTPError as e:
        content = e.read()
        try:
            return e.code, json.loads(content), e.headers
        except ValueError:
            return e.code, content, e.headers


def test_binds_to_loopback_only(server):
    assert server.server_address[0] == "127.0.0.1"
    assert server.url.startswith("http://127.0.0.1:")


def test_page_and_static_files(server):
    st, html, hd = http(server, "GET", "/", token=False, raw=True)
    html = html.decode()
    assert st == 200 and hd["Content-Type"].startswith("text/html")
    assert f'content="{server.session.token}"' in html and "{{TOKEN}}" not in html
    assert "default-src 'self'" in hd["Content-Security-Policy"]
    # relative URLs (static/...): the same page also works under a sub-path (the web version on GitHub Pages)
    assert 'src="/' not in html and 'href="/' not in html
    for js in re.findall(r'src="(static/[^"]+)"', html) + re.findall(r'href="(static/[^"]+)"', html):
        st, body, hd = http(server, "GET", "/" + js, token=False, raw=True)
        assert st == 200 and len(body) > 100, js
    st, body, hd = http(server, "GET", "/static/plotly.min.js", token=False, raw=True)
    assert st == 200 and b"plotly" in body[:2000].lower() and "javascript" in hd["Content-Type"]
    st, body, hd = http(server, "GET", "/static/app.js", token=False, raw=True)
    # the menu bar of requirements 5.2 and the items of 5.3
    text = body.decode()
    for item in ('title: "Model"', 'title: "File"', 'title: "Plot"', 'title: "Modules"', 'title: "Options"',
                 'title: "View"', 'title: "Help"', '"Converters"', '"Export to ANSYS"', '"Export to STRUDL"',
                 '"Export Image"', '"Export Table"', '"Spectrum TFU-TFI"', '"Soil Layers"', '"Soil Properties"',
                 '"Non Uniform Soil Field"', '"Process Animation Frame List"', '"Deformed Shape"', '"Location"',
                 '"Extension"', '"NONLINEAR"', '"Windows Settings"', '"Shader Options"', '"Reset Plot"',
                 '"Check Errors"', '"Command Window"', '"Command Echo"', '"Output Confirmation"', '"Comments"',
                 '"Warnings & Errors"', '"Main Toolbar"', '"Plot Toolbar"', '"About"', '"Verification"'):
        assert item in text, item


def test_help_pages_and_documentation_figures(server):
    """Help > Help: documents rendered by the server; figures under docs/ served as images only."""
    st, h, _ = http(server, "GET", "/api/help")
    assert st == 200 and h["home"] == "docs/index.md" and any(d["id"] == "docs/user/GUI.md" for d in h["docs"])
    st, d, _ = http(server, "GET", "/api/help/doc?name=docs/verification/VERIFICATION_MANUAL.md")
    assert st == 200 and d["title"] == "SASSI-EDU Verification Manual"
    src = re.search(r'<img src="(docs/[^"]+\.png)"', d["html"]).group(1)          # relative to the page at /
    st, body, hd = http(server, "GET", "/" + src, token=False, raw=True)
    assert st == 200 and hd["Content-Type"] == "image/png" and body[:4] == b"\x89PNG"
    for p in ("/docs/index.md", "/docs/../README.md", "/docs/..%2F..%2Fpyproject.toml", "/docs/internal/wave3_packages.md",
              "/docs/verification/figures/missing.png", "/docs/"):
        st, body, hd = http(server, "GET", p, token=False, raw=True)
        assert st == 404, p
    st, r, _ = http(server, "GET", "/api/help/doc?name=../pyproject.toml")
    assert st == 404
    st, r, _ = http(server, "GET", "/api/help/doc?name=docs/index.md", token=False)
    assert st == 403                                                     # documents need the session token too


def test_disconnected_browser_is_not_reported(server, capsys):
    """A page reload drops the long poll: the server keeps quiet about the broken pipe."""
    for exc in (BrokenPipeError(32, "Broken pipe"), ConnectionResetError(54, "reset")):
        try:
            raise exc
        except OSError:
            server.handle_error(None, ("127.0.0.1", 1))
    assert "Traceback" not in capsys.readouterr().err
    try:
        raise RuntimeError("a real bug")
    except RuntimeError:
        server.handle_error(None, ("127.0.0.1", 1))
    assert "a real bug" in capsys.readouterr().err


def test_katex_files(server):
    """KaTeX typesets the LaTeX formulas of the lessons and Help pages: script, style sheet and the woff2
    fonts its style sheet names are served from sassi/ui/static/katex (offline, no CDN)."""
    st, body, hd = http(server, "GET", "/static/katex/katex.min.js", token=False, raw=True)
    assert st == 200 and b"katex" in body[:5000].lower() and "javascript" in hd["Content-Type"]
    st, css, hd = http(server, "GET", "/static/katex/katex.min.css", token=False, raw=True)
    assert st == 200 and hd["Content-Type"].startswith("text/css")
    fonts = sorted(set(re.findall(r"url\((fonts/[^)]+)\)", css.decode())))
    assert fonts and all(f.endswith(".woff2") for f in fonts)
    for f in fonts:
        st, body, hd = http(server, "GET", "/static/katex/" + f, token=False, raw=True)
        assert st == 200 and hd["Content-Type"] == "font/woff2" and body[:4] == b"wOF2", f
    for bad in ("/static/katex/../app.js", "/static/katex/fonts/../../server.py", "/static/katex/.x",
                "/static/katex/fonts/sub/x.woff2", "/static/katex/missing.js"):
        assert http(server, "GET", bad, token=False, raw=True)[0] == 404, bad


def test_static_traversal_refused(server):
    for p in ("/static/../server.py", "/static/..%2Fserver.py", "/static/.hidden", "/static/sub/x.js", "/etc/passwd"):
        st, body, hd = http(server, "GET", p, token=False, raw=True)
        assert st == 404, p


def test_token_and_host_checks(server):
    st, body, hd = http(server, "GET", "/api/state", token=False)
    assert st == 403 and "token" in body["error"]
    st, body, hd = http(server, "POST", "/api/command", {"line": "N,1,0,0,0"}, token=False)
    assert st == 403 and not server.session.interp.model.nodes        # CSRF: nothing executed
    st, body, hd = http(server, "GET", "/api/state", headers={"Host": "evil.example.com"})
    assert st == 403 and "loopback" in body["error"]
    st, body, hd = http(server, "GET", "/api/state", headers={"Host": f"localhost:{server.server_port}"})
    assert st == 200


def test_bad_bodies(server):
    req = urllib.request.Request(server.url + "api/command", data=b"not json", method="POST",
                                 headers={"X-SASSI-Token": server.session.token})
    with pytest.raises(urllib.error.HTTPError) as ei:
        urllib.request.urlopen(req, timeout=10)
    assert ei.value.code == 400


def test_command_options_model_round_trip_over_http(server, tmp_path):
    lines = ["N,1,0,0,0", "N,2,1,0,0", "N,3,1,1,0", "N,4,0,1,0", "GROUP,1,SHELL", "E,1,1,2,3,4", "INT,1,4,1,1"]
    st, r, _ = http(server, "POST", "/api/command", {"lines": lines})
    assert st == 200 and r["ok"] and len(r["results"]) == 7
    st, m, _ = http(server, "GET", "/api/model")
    assert m["counts"]["nodes"] == 4 and m["interaction"] == [1, 2, 3, 4]
    st, o, _ = http(server, "GET", "/api/options/HOUSE")
    assert o["values"]["records"]["HOUSE"]["gravity"] == 32.2
    st, r, _ = http(server, "POST", "/api/options/ANALYSIS", {"records": {"HOUSE": {"gravity": 9.81, "gelev": -2}}})
    assert st == 200 and r["commands"] == ["HOUSE,9.81,-2,0,2,0,0,0,0,0"]
    assert server.session.interp.model.gravity == 9.81
    st, r, _ = http(server, "POST", "/api/options/ANALYSIS", {"records": {"HOUSE": {"gravity": -1}}})
    assert st == 422 and "Error 1 :" in r["problems"][0]
    st, r, _ = http(server, "POST", "/api/command", {"line": "MODELPLOT"})
    st, d, _ = http(server, "GET", "/api/plot/1")
    assert d["family"] == "3d" and d["scene"]["faces"][0] == [0, 1, 2, 3]
    st, ev, _ = http(server, "GET", "/api/events?since=0&timeout=0")
    assert any(e["type"] == "plot" and e["event"] == "open" for e in ev["events"])


def test_file_api_over_http(server, tmp_path):
    (tmp_path / "a.RS").write_text("0.1 0.2\n1 0.5\n", encoding="utf-8")
    st, lst, _ = http(server, "GET", "/api/files")
    assert [f["name"] for f in lst["files"]] == ["a.RS"] and lst["files"][0]["plot"] == "spec"
    st, f, _ = http(server, "GET", "/api/file?name=a.RS")
    assert f["text"].startswith("0.1")
    st, f, _ = http(server, "GET", "/api/file?name=..%2F..%2F..%2Fetc%2Fpasswd")
    assert st == 403


def test_module_run_over_http(server, tmp_path):
    model = f"""MDL,run,{tmp_path / 'run'}
N,1,0,0,0
N,2,1,0,0
N,3,1,1,0
N,4,0,1,0
INT,1,4,1,1
M,1,3e6,0.25,0.15,0.05,0.05
GROUP,1,SHELL
E,1,1,2,3,4
THICK,1,1,1,0.5
L,1,10,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
TOPL,1
SITE,0,1,0,10,2,1,0,1,2048,1,0,0.005,4096,1
FREQ,1,4,8
AOPT,0,0,0,1,0,0,0,0,0,0,0,0,0,0
AFWRITE"""
    st, r, _ = http(server, "POST", "/api/command", {"lines": model.splitlines()})
    assert r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    st, job, _ = http(server, "POST", "/api/run/SITE")
    assert st == 200
    t0 = time.time()
    while True:
        st, d, _ = http(server, "GET", f"/api/jobs/{job['id']}?since=0")
        if d["state"] not in ("starting", "running") or time.time() - t0 > 120:
            break
        time.sleep(0.1)
    assert d["state"] == "done" and (tmp_path / "run" / "FILE1").exists()
    assert any("SITE finished with status OK" in m["text"] for m in d["messages"])


def _free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def test_main_entry_point_and_exit(tmp_path, capsys):
    """``sassi-gui --port P --no-browser --model-dir D``: serves until Model > Exit (POST /api/exit)."""
    port = _free_port()
    rc = {}
    th = threading.Thread(target=lambda: rc.setdefault("rc", srvmod.main(
        ["--port", str(port), "--no-browser", "--model-dir", str(tmp_path), "--settings-dir", str(tmp_path / "s")])),
        daemon=True)
    th.start()
    url = f"http://127.0.0.1:{port}/"
    html = None
    for _ in range(100):
        try:
            with urllib.request.urlopen(url, timeout=2) as r:
                html = r.read().decode()
            break
        except OSError:
            time.sleep(0.05)
    assert html is not None
    token = re.search(r'name="sassi-token" content="([0-9a-f]+)"', html).group(1)
    req = urllib.request.Request(url + "api/state", headers={"X-SASSI-Token": token})
    with urllib.request.urlopen(req, timeout=5) as r:
        assert json.loads(r.read())["cwd"] == str(tmp_path.resolve())
    req = urllib.request.Request(url + "api/exit", data=b"{}", method="POST",
                                 headers={"X-SASSI-Token": token, "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=5) as r:
        assert json.loads(r.read())["ok"]
    th.join(10)
    assert not th.is_alive() and rc["rc"] == 0
    assert (tmp_path / "s" / "SASSIini.xml").exists()
    assert "SASSI-EDU GUI running at" in capsys.readouterr().out


def test_main_rejects_missing_model_dir(tmp_path):
    with pytest.raises(SystemExit):
        srvmod.main(["--no-browser", "--model-dir", str(tmp_path / "missing")])


# ------------------------------------------------------------------ HTTP/1.1 keep-alive and request bodies
def _raw(srv, data: bytes, wait: float = 1.0) -> bytes:
    s = socket.create_connection(("127.0.0.1", srv.server_port))
    s.sendall(data)
    s.settimeout(wait)
    out = b""
    t0 = time.time()
    try:
        while time.time() - t0 < 3 * wait:
            chunk = s.recv(65536)
            if not chunk:
                break
            out += chunk
    except socket.timeout:
        pass
    s.close()
    return out


def _post(path: bytes, body: bytes, headers: bytes = b"") -> bytes:
    return (b"POST " + path + b" HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: application/json\r\n" + headers +
            b"Content-Length: %d\r\n\r\n" % len(body)) + body


def test_refused_requests_keep_the_connection_in_step(server):
    """Every refused request (405 static POST, 403 Host, 403 token) has its body consumed: the next
    request of a keep-alive connection is answered normally."""
    tok = server.session.token.encode()
    body = b'{"line": "N,1,0,0,0"}'
    data = (_post(b"/index.html", body) +
            _post(b"/api/command", body, b"Host: evil.example.com\r\nX-SASSI-Token: " + tok + b"\r\n").replace(
                b"Host: 127.0.0.1\r\n", b"", 1) +
            _post(b"/api/command", body, b"X-SASSI-Token: wrong\r\n") +
            b"GET /api/state HTTP/1.1\r\nHost: 127.0.0.1\r\nX-SASSI-Token: " + tok + b"\r\nConnection: close\r\n\r\n")
    out = _raw(server, data)
    assert re.findall(rb"HTTP/1\.[01] (\d{3})", out) == [b"405", b"403", b"403", b"200"], out[:400]
    assert not server.session.interp.model.nodes


@pytest.mark.parametrize("length", [b"-5", b"abc", b"1e3"])
def test_malformed_content_length_is_refused_without_blocking(server, length):
    tok = server.session.token.encode()
    req = (b"POST /api/command HTTP/1.1\r\nHost: 127.0.0.1\r\nX-SASSI-Token: " + tok +
           b"\r\nContent-Length: " + length + b"\r\n\r\n{}")
    t0 = time.time()
    out = _raw(server, req, wait=2.0)
    assert out.startswith(b"HTTP/1.1 400") and b"Content-Length" in out and b"Connection: close" in out
    assert time.time() - t0 < 5.0                     # the connection was closed by the server


def test_oversized_body_closes_the_connection_unread(server):
    tok = server.session.token.encode()
    req = (b"POST /api/command HTTP/1.1\r\nHost: 127.0.0.1\r\nX-SASSI-Token: " + tok +
           b"\r\nContent-Length: %d\r\n\r\n" % (srvmod.MAX_BODY + 1))
    t0 = time.time()
    out = _raw(server, req, wait=2.0)
    assert out.startswith(b"HTTP/1.1 400") and b"too large" in out and b"Connection: close" in out
    assert time.time() - t0 < 5.0


# ------------------------------------------------------------------ Ctrl-C in the terminal
SURFACE = """MDL,run,{dir}
N,1,0,0,0
N,2,1,0,0
N,3,1,1,0
N,4,0,1,0
INT,1,4,1,1
D,1,4,1,1,ROTZ
M,1,3e6,0.25,0.15,0.05,0.05
GROUP,1,SHELL
E,1,1,2,3,4
THICK,1,1,1,0.5
L,1,10,0.12,1500,800,0.05,0.05
L,2,10,0.13,2400,1200,0.02,0.02
TOPL,1
SITE,0,1,0,10,2,1,0,1,2048,1,0,0.005,4096,1
FREQ,1,4,8
AOPT,0,0,0,1,1,0,0,0,0,0,0,0,0,0
AFWRITE
"""


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


@pytest.mark.skipif(os.name != "posix", reason="SIGINT and process groups are POSIX")
def test_ctrl_c_stops_the_running_module(tmp_path):
    """Ctrl-C in the terminal of ``sassi-gui`` terminates a running module run: the worker leads its
    own process group (no terminal SIGINT), so the server must stop it before it exits."""
    port = _free_port()
    pidfile = tmp_path / "module.pid"
    exe = tmp_path / "slow_site"
    exe.write_text(f"#!{sys.executable}\nimport os, sys, time\nsys.stdin.read()\n"
                   f"open({str(pidfile)!r}, 'w').write(str(os.getpid()))\nprint('working', flush=True)\n"
                   f"time.sleep(120)\n", encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    proc = subprocess.Popen([sys.executable, "-m", "sassi.ui.server", "--port", str(port), "--no-browser",
                             "--model-dir", str(tmp_path), "--settings-dir", str(tmp_path / "s")],
                            cwd=str(Path(srvmod.__file__).resolve().parents[2]), stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    try:
        url = f"http://127.0.0.1:{port}/"
        html = None
        for _ in range(200):
            try:
                with urllib.request.urlopen(url, timeout=2) as r:
                    html = r.read().decode()
                break
            except OSError:
                time.sleep(0.05)
        token = re.search(r'name="sassi-token" content="([0-9a-f]+)"', html).group(1)

        def api(path, body):
            req = urllib.request.Request(url + path, data=json.dumps(body).encode(), method="POST",
                                         headers={"X-SASSI-Token": token, "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())
        assert api("api/settings/locations", {"locations": {"SITE": str(exe)}})["ok"]
        lines = [ln for ln in SURFACE.format(dir=tmp_path / "run").splitlines() if ln.strip()]
        assert api("api/command", {"lines": lines})["ok"]
        api("api/run/SITE", {})
        t0 = time.time()
        while not (pidfile.exists() and pidfile.read_text()) and time.time() - t0 < 30:
            time.sleep(0.05)
        child = int(pidfile.read_text())
        assert _alive(child)
        proc.send_signal(signal.SIGINT)                    # Ctrl-C
        assert proc.wait(30) == 0
        t0 = time.time()
        while _alive(child) and time.time() - t0 < 10:
            time.sleep(0.05)
        assert not _alive(child)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait(10)
