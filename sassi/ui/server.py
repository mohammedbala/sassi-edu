"""``sassi-gui``: the local web server of the SASSI-EDU GUI (ARCHITECTURE section 9, requirements 5).

Usage::

    sassi-gui [--port 8765] [--no-browser] [--model-dir DIR]
    python -m sassi.ui.server ...

The server is a standard-library :class:`http.server.ThreadingHTTPServer` bound to **127.0.0.1
only**.  It serves the front end (``sassi/ui/static``), Plotly.js from the installed ``plotly``
package (so the GUI works offline; the file is not copied into this package) and the JSON API of
:class:`sassi.ui.api.GuiSession`.

Safety:

* loopback only (127.0.0.1); requests whose ``Host`` header is not a loopback name are refused
  (protection against DNS rebinding);
* every ``/api/`` request must carry the session token (``X-SASSI-Token`` header), which only the
  page served by this server knows -- another web site cannot drive the interpreter (CSRF);
* files are read and written only inside the model directories and the working directory
  (:func:`sassi.ui.files.safe_path`); static files only from ``sassi/ui/static``; ``/docs/<path>``
  serves only the PNG / JPEG / GIF figures of the documentation under ``docs/`` (Help pages,
  :func:`sassi.ui.helpdocs.image_file`);
* no shell is ever executed: module runs use a Python worker process (or the executable chosen
  in Modules > Location) with an argument list.
"""
from __future__ import annotations

import argparse
import json
import mimetypes
import posixpath
import re
import sys
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Optional, Tuple
from urllib.parse import parse_qs, unquote, urlsplit

from .api import GuiSession
from .helpdocs import IMAGE_TYPES, image_file

STATIC_DIR = Path(__file__).resolve().parent / "static"
#: content types of the KaTeX files (static/katex)
KATEX_TYPES = {".js": "application/javascript; charset=utf-8", ".css": "text/css; charset=utf-8",
               ".woff2": "font/woff2", ".txt": "text/plain; charset=utf-8"}
DEFAULT_PORT = 8765
MAX_BODY = 32 * 1024 * 1024
#: a refused request's body up to this size is read and discarded (keep-alive stays in step);
#: a larger one closes the connection instead of being read
MAX_DRAIN = 1024 * 1024
LOOPBACK_HOSTS = ("127.0.0.1", "localhost", "[::1]")


def plotly_js_path() -> Optional[Path]:
    """``plotly.min.js`` of the installed plotly package (served at /static/plotly.min.js)."""
    try:
        import plotly
    except ImportError:
        return None
    p = Path(plotly.__file__).resolve().parent / "package_data" / "plotly.min.js"
    return p if p.is_file() else None


def _host_ok(host: str, port: int) -> bool:
    host = (host or "").strip().lower()
    if not host:
        return True                      # HTTP/1.0 clients without Host
    if host.startswith("["):
        name = host.split("]")[0] + "]"
    else:
        name = host.split(":")[0]
    return name in LOOPBACK_HOSTS


class Handler(BaseHTTPRequestHandler):
    """Request handler: static files, index page with the session token, JSON API."""

    server_version = "SASSI-EDU-GUI"
    protocol_version = "HTTP/1.1"

    @property
    def session(self) -> GuiSession:
        return self.server.session          # type: ignore[attr-defined]

    def log_message(self, fmt, *args):          # quiet console (errors only)
        if getattr(self.server, "verbose", False):
            sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    # ------------------------------------------------------------------ helpers
    def _send(self, status: int, body: bytes, ctype: str, extra: Optional[dict] = None) -> None:
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        if self.close_connection and self.request_version == "HTTP/1.1":
            self.send_header("Connection", "close")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, status: int, payload) -> None:
        body = json.dumps(payload, allow_nan=False, default=str).encode("utf-8")
        self._send(status, body, "application/json; charset=utf-8")

    def _content_length(self) -> Optional[int]:
        """Declared body length: 0 when absent, None when malformed or negative."""
        raw = (self.headers.get("Content-Length") or "").strip()
        if not raw:
            return 0
        if not raw.isdigit():
            return None
        return int(raw)

    def _discard_body(self) -> None:
        """A refused request: consume its body so that the next request of a keep-alive connection is
        parsed from its own first line; a body that is malformed or too large closes the connection."""
        if self._body_done:
            return
        self._body_done = True
        n = self._content_length()
        if n is None or n > MAX_DRAIN:
            self.close_connection = True
            return
        while n > 0:
            chunk = self.rfile.read(min(n, 65536))
            if not chunk:
                self.close_connection = True
                return
            n -= len(chunk)

    def _read_body(self) -> Tuple[Optional[dict], Optional[str]]:
        n = self._content_length()
        if n is None:
            self._body_done = True
            self.close_connection = True
            return None, "malformed Content-Length"
        if n > MAX_BODY:
            self._discard_body()                  # too large to read: the connection is closed
            return None, "request body too large"
        self._body_done = True
        raw = self.rfile.read(n) if n else b""
        if len(raw) < n:
            self.close_connection = True
            return None, "incomplete request body"
        if not raw:
            return {}, None
        try:
            data = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, ValueError):
            return None, "request body is not JSON"
        if not isinstance(data, dict):
            return None, "request body must be a JSON object"
        return data, None

    # ------------------------------------------------------------------ verbs
    def do_GET(self):
        self._dispatch("GET")

    def do_HEAD(self):
        self._dispatch("GET")

    def do_POST(self):
        self._dispatch("POST")

    def _dispatch(self, method: str) -> None:
        self._body_done = False
        try:
            self._route_request(method)
        finally:
            # every answer leaves the connection at the start of the next request (HTTP/1.1
            # keep-alive): a request refused before its body was read has it discarded here
            self._discard_body()

    def _route_request(self, method: str) -> None:
        if self._content_length() is None:
            self._body_done = True
            self.close_connection = True
            self._json(400, {"error": "malformed Content-Length"})
            return
        if not _host_ok(self.headers.get("Host", ""), self.server.server_port):
            self._json(403, {"error": "this server only answers loopback requests"})
            return
        url = urlsplit(self.path)
        path = url.path
        if path.startswith("/api/"):
            if self.headers.get("X-SASSI-Token", "") != self.session.token:
                self._json(403, {"error": "missing or wrong session token (open the GUI from its own page)"})
                return
            body: Optional[dict] = {}
            if method == "POST":
                body, err = self._read_body()
                if err:
                    self._json(400, {"error": err})
                    return
            status, payload = self.session.route(method, path, parse_qs(url.query), body or {})
            try:
                self._json(status, payload)
            except (ValueError, TypeError) as exc:      # NaN or a non-serialisable value
                self._json(500, {"error": f"response could not be encoded: {exc}"})
            return
        if method != "GET":
            self._json(405, {"error": "method not allowed"})
            return
        self._static(path)

    def _static(self, path: str) -> None:
        if path in ("/", "/index.html"):
            p = STATIC_DIR / "index.html"
            text = p.read_text(encoding="utf-8").replace("{{TOKEN}}", self.session.token)
            self._send(200, text.encode("utf-8"), "text/html; charset=utf-8",
                       {"Content-Security-Policy": "default-src 'self'; script-src 'self' 'unsafe-eval'; "
                                                   "style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; "
                                                   "connect-src 'self'; object-src 'none'; frame-ancestors 'none'"})
            return
        if path == "/static/plotly.min.js":
            p = plotly_js_path()
            if p is None:
                self._json(404, {"error": "plotly is not installed (pip install plotly)"})
                return
            self._send(200, p.read_bytes(), "application/javascript; charset=utf-8")
            return
        m = re.match(r"^/static/katex/((?:fonts/)?[A-Za-z0-9_.\-]+)$", path)
        if m and not m.group(1).split("/")[-1].startswith("."):
            # KaTeX (typesetting of the LaTeX formulas): katex.min.js / .css and the woff2 fonts
            p = STATIC_DIR / "katex" / m.group(1)
            if not p.is_file():
                self._json(404, {"error": "not found"})
                return
            ctype = KATEX_TYPES.get(p.suffix.lower(), "application/octet-stream")
            self._send(200, p.read_bytes(), ctype, {"Cache-Control": "max-age=86400"})
            return
        if path.startswith("/static/"):
            name = path[len("/static/"):]
            if not name or "/" in name or "\\" in name or name.startswith("."):
                self._json(404, {"error": "not found"})
                return
            p = STATIC_DIR / name
            if not p.is_file():
                self._json(404, {"error": "not found"})
                return
            ctype = mimetypes.guess_type(name)[0] or "application/octet-stream"
            if ctype.startswith("text/") or ctype.endswith("javascript"):
                ctype += "; charset=utf-8"
            self._send(200, p.read_bytes(), ctype)
            return
        if path.startswith("/docs/"):
            # figures of the Help pages (docs/verification/figures/*.png ...): images under docs/ only
            rel = posixpath.normpath("docs/" + unquote(path[len("/docs/"):]))
            p = image_file(self.session.helpdocs.root, rel)
            if p is None:
                self._json(404, {"error": "not found"})
                return
            self._send(200, p.read_bytes(), IMAGE_TYPES[p.suffix.lower()])
            return
        if path == "/favicon.ico":
            self._send(204, b"", "image/x-icon")
            return
        self._json(404, {"error": "not found"})


class GuiServer(ThreadingHTTPServer):
    """The HTTP server holding one :class:`GuiSession`."""

    daemon_threads = True
    allow_reuse_address = True

    def __init__(self, port: int = DEFAULT_PORT, session: Optional[GuiSession] = None, verbose: bool = False):
        super().__init__(("127.0.0.1", port), Handler)
        self.session = session or GuiSession()
        self.session.shutdown = self._stop
        self.verbose = verbose

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.server_port}/"

    def handle_error(self, request, client_address) -> None:
        """A browser that left (page reloaded or closed during a long poll) is not an error worth a
        traceback on the console; anything else is reported as usual."""
        if isinstance(sys.exc_info()[1], (BrokenPipeError, ConnectionResetError, ConnectionAbortedError)):
            return
        super().handle_error(request, client_address)

    def _stop(self) -> None:
        threading.Thread(target=self.shutdown, daemon=True).start()


def make_server(port: int = DEFAULT_PORT, model_dir: Optional[str] = None, settings_dir: Optional[str] = None,
                verbose: bool = False) -> GuiServer:
    """Create the server (port 0 = any free port); a busy port falls back to a free one."""
    session = GuiSession(cwd=model_dir, settings_dir=settings_dir)
    try:
        return GuiServer(port, session, verbose)
    except OSError:
        if port == 0:
            raise
        return GuiServer(0, session, verbose)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="sassi-gui", description="SASSI-EDU graphical user interface (local web app)")
    ap.add_argument("--port", type=int, default=DEFAULT_PORT, help=f"port on 127.0.0.1 (default {DEFAULT_PORT}; "
                                                                     "0 = any free port)")
    ap.add_argument("--no-browser", action="store_true", help="do not open the web browser")
    ap.add_argument("--model-dir", default=None, help="initial working directory (default: current directory)")
    ap.add_argument("--settings-dir", default=None, help="directory of SASSIini.xml / SASSIdb.xml (default: the "
                                                          "per-user settings directory)")
    ap.add_argument("--verbose", action="store_true", help="log every HTTP request")
    args = ap.parse_args(argv)
    if args.model_dir and not Path(args.model_dir).expanduser().is_dir():
        ap.error(f"--model-dir {args.model_dir} is not a directory")
    srv = make_server(args.port, str(Path(args.model_dir).expanduser()) if args.model_dir else None,
                      args.settings_dir, args.verbose)
    print(f"SASSI-EDU GUI running at {srv.url}  (Ctrl-C or Model > Exit to stop)", flush=True)
    if not args.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(srv.url)).start()
    try:
        srv.serve_forever(poll_interval=0.25)
    except KeyboardInterrupt:
        print("\nstopping", flush=True)
    finally:
        stop_session(srv)
    return 0


def stop_session(srv: GuiServer) -> None:
    """Server stop (Ctrl-C in the terminal, Model > Exit): a running module worker is terminated --
    it leads its own process group, so the terminal's SIGINT does not reach it and it would keep
    writing files after the GUI has gone (UI-04) -- then SASSIini.xml is written and the socket closed."""
    try:
        srv.session.jobs.shutdown()
    except Exception as exc:                       # pragma: no cover - stop must always complete
        print(f"could not stop the running job: {exc!r}", file=sys.stderr, flush=True)
    try:
        srv.session.save_settings()
    except OSError:
        pass
    srv.server_close()


if __name__ == "__main__":
    raise SystemExit(main())
