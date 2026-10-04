"""Adversarial review tests of the web GUI work package (ARCHITECTURE section 9; requirements 5.1-5.10,
UI-01..UI-08; D-UI-03, D-UI-04, D-UI-07; rule L17).

The checks are written from the requirements, not from the implementation:

* L17 idempotence: posting the current dialog values back must emit no command, and what a commit
  emits must read back unchanged (fixed point of ``values`` / ``commit_commands``);
* the dialog must show what the analysis will use: a value displayed for an entry that is *not*
  stored must equal what AFWRITE writes for it (SITE waves, SOIL layer requests, SITE Frequency 2);
* every widget choice of the Options dialogs is an admissible code of its option record;
* every command text the front end submits names a registered command;
* result-file classification follows the naming formulas of requirements 1.8;
* path safety, Host / token checks and the event log from first principles.

Tests marked ``xfail(strict=True)`` document defects found by the review (the reason names the
defect); they turn into XPASS failures once the defect is fixed, so the marker must then be removed.
"""
from __future__ import annotations

import copy
import re
import socket
import threading
import time
from pathlib import Path

import pytest

from sassi.conventions import element_result_name, layer_th_name, nodal_result_name, nodal_rs_name
from sassi.plotting.state import plot_state
from sassi.prep import Interpreter
from sassi.prep.check import Checker
from sassi.prep.options import get_entries, get_record, record_class
from sassi.prep.registry import load_commands, lookup
from sassi.ui import dialogs as D
from sassi.ui import files as F
from sassi.ui.api import GuiSession
from sassi.ui.events import EventLog
from sassi.ui.server import _host_ok, make_server

STATIC = Path(D.__file__).resolve().parent / "static"

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
WAVE,2,1,1,1,0
FREQ,1,4,8
POINT,0,0,0.45
AOPT,0,0,0,1,1,0,0,0,0,0,0,0,0,0
"""


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


@pytest.fixture
def S(tmp_path):
    s = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    yield s
    run = s.jobs.running()
    if run is not None:
        s.jobs.cancel(run.id)
        t0 = time.time()
        while s.jobs.get(run.id).active and time.time() - t0 < 20:
            time.sleep(0.05)


def _run_all(interp, lines):
    for ln in lines:
        assert interp.execute(ln), (ln, interp.sink.texts()[-3:])


def _commit_and_run(interp, payload):
    lines, notes = D.commit_commands(interp, payload)
    _run_all(interp, lines)
    return lines, notes


# ======================================================================================
# 1. L17 idempotence of the Options dialogs (fixed point of values -> commit)
# ======================================================================================
RICH = ["SITE,0,1,0,10,2,1,0,1,1024,1,0,0.01,2048,1", "WAVE,2,1,0.5,0.7,10", "WAVE,1,2,1,1,0",
        "HOUSE,9.81,-3,0,2,1,1,0,0,0", "INCOH,0.2,0.2,0.3,0.5,2,1,3,0,5,7,180", "WPASS,2000,30,1",
        "ME,1,1,10,0,0,0", "AMP,1,1,1.2", "DAMP,0.02,0.05", "TOPL,1,2", "NOUT,1,1,0,0,0,1,0,1,3-6",
        "RDND,15,1,0,1,0,0,0", "GROUP,3,BEAMS", "EOUT,1,1,0,0,0,0,0,0,0,0,0,0,3,1-4", "SACC,1,2,1", "SRS,1,1,0",
        "SPRO,1,2,Sand", "THFILE,acc.acc", "EQTIT,my title", "CORR,1,0,0.3", "RSIN,1,a.rsi",
        "MOTIONX,1,2,0,1", "AOPT,1,1,0,1,1,1,0,0,1,0,1,1,1,0", "EQUAKE,0,24,11975,0.05,25,0,1",
        "SOIL,3000,32.2,0,1,1,6,0.65,1,0", "RELD,3,1,0", "ANALYS,0,0,0,0,1,0,0,0,0,0,2,0"]


def test_posting_the_current_values_back_emits_nothing(ui):
    """Every part of the dialog payload, unchanged, must produce no command (OK sends only changes)."""
    _run_all(ui, RICH)
    v = D.values(ui)
    payload = {k: copy.deepcopy(v[k]) for k in ("records", "indexed", "strings", "lists", "requests")}
    assert D.commit_commands(ui, payload) == ([], [])


def test_commit_read_back_is_a_fixed_point(ui):
    """What a commit emits reads back as the committed values, and committing them again is a no-op."""
    _run_all(ui, ["TOPL,1,2,3"])
    payload = {"records": {"SITE": {"nl": 12, "nft": 8192, "delt": 0.0025, "freq2": 3000, "cm": 2},
                           "HOUSE": {"gravity": 9.81, "gelev": -1.5}, "MOTION": {"freq1": 0.2, "fstep": 200}},
               "indexed": {"SPRO": {"2": {"prop": 1, "dynprop": "Clay, soft"}}, "ME": {"2": {"nfirst": 3, "nlast": 9}}},
               "strings": {"THTIT": "Design motion, X"}, "lists": {"DAMP": "0.03, 0.07", "AMP[2]": "1 1.1 1.2"},
               "requests": {"RDND": [{"node": 4, "x": 1, "zz": 1}]}}
    lines, _ = _commit_and_run(ui, payload)
    assert lines
    v = D.values(ui)
    s = v["records"]["SITE"]
    assert (s["nl"], s["nft"], s["delt"], s["freq2"], s["cm"]) == (12, 8192, 0.0025, 3000, 2)
    assert v["indexed"]["SPRO"]["2"]["dynprop"] == "Clay, soft" and v["indexed"]["ME"]["2"]["nlast"] == 9
    assert v["strings"]["THTIT"] == "Design motion, X" and v["lists"]["DAMP"] == "0.03 0.07"
    assert v["lists"]["AMP[2]"] == "1 1.1 1.2"
    assert D.commit_commands(ui, payload) == ([], [])


def test_ids_text_and_parse_ids_are_inverse():
    for ids in ([1], [5, 1, 2, 3, 9], list(range(1, 40)) + [100, 102, 103], [7, 7, 7, 3], [10, 9, 8, 1]):
        assert D.parse_ids(D._ids_text(ids)) == sorted(set(ids))


def test_dialog_choices_are_admissible_codes():
    """Each radio / select / exclusive-check code of the Analysis dialog is an admissible value of the
    option field it is bound to (the typed records carry the admissible codes)."""
    bad = []
    for tab in D.form("ANALYSIS")["tabs"]:
        for g in tab["groups"]:
            for it in g["items"]:
                if it.get("w") not in ("radio", "select", "excl") or not it.get("choices"):
                    continue
                m = re.fullmatch(r"(\w+)(?:\[\w+\])?\.(\w+)", it["path"])
                if not m:
                    continue
                cls = record_class(m.group(1))
                fl = cls.FIELDS[cls.field_index(m.group(2)) - 1]
                if not getattr(fl, "choices", ()):
                    continue
                for val, _lab in it["choices"]:
                    if val in (it.get("disabled_choices") or []):
                        continue
                    if val not in fl.choices and not (it["w"] == "excl" and val == 0):
                        bad.append((it["path"], val, fl.choices))
    assert not bad, bad


# ======================================================================================
# 2. D-UI-04 forcing rules applied on the stored state, not only on the payload
# ======================================================================================
def test_forced_values_use_the_stored_state(ui):
    _run_all(ui, ["HOUSE,32.2,0,0,2,0,1,0,0,0"])                       # stored: 3D, incoherent
    lines, notes = D.commit_commands(ui, {"records": {"HOUSE": {"dim": 1}}})
    assert lines == ["HOUSE,32.2,0,0,1,0,0,0,0,0"] and notes                  # 2D -> coherent
    ui2 = Interpreter(cwd=ui.cwd)
    _run_all(ui2, ["WPASS,1000000000,0,4"])                                   # stored: unlagged model 4
    lines, notes = D.commit_commands(ui2, {"records": {"HOUSE": {"gelev": -2}}})
    assert lines == ["HOUSE,32.2,-2,0,2,0,0,1,0,0"] and "unlagged coherency model 4" in notes[0]


# ======================================================================================
# 3. The dialog shows what AFWRITE writes (defects)
# ======================================================================================
def _waves_on(interp):
    return sorted((int(w.type), int(w.opt)) for w in Checker(interp.model).effective_waves() if int(w.opt) != 0)


def _waves_shown(interp):
    v = D.values(interp)
    out = []
    for t in range(1, 6):
        e = v["indexed"].get("WAVE", {}).get(str(t)) or v["wave_defaults"][str(t)]
        if int(e["opt"]) != 0:
            out.append((t, int(e["opt"])))
    return sorted(out)


def test_wave_pages_shown_match_the_effective_waves_initially(ui):
    assert _waves_shown(ui) == _waves_on(ui) == [(2, 1)]


def test_adding_r_waves_keeps_the_displayed_sv_field(ui):
    v = D.values(ui)
    _commit_and_run(ui, {"indexed": {"WAVE": {"1": dict(v["wave_defaults"]["1"], opt=1)}}})
    assert _waves_shown(ui) == _waves_on(ui)
    assert (2, 1) in _waves_on(ui)


@pytest.mark.parametrize("cmd, entry", [
    ("SACC", {"opt": 1, "outcrop": 0}),                              # Compute Maximum
    ("SRS", {"save": 1, "outcrop": 0}),                              # Save Response Spectrum
    ("SSTR", {"opt1": 1, "opt2": 0, "opt3": 1, "opt4": 0}),          # Compute Stresses + Strains
])
def test_soil_layer_requests_at_dialog_defaults_reach_the_deck(ui, cmd, entry):
    _commit_and_run(ui, {"indexed": {cmd: {"3": dict(entry)}}})
    stored = {int(k): r for k, r in get_entries(ui.model, cmd)}
    assert 3 in stored
    assert all(int(stored[3].get(f)) == v for f, v in entry.items())


def test_stochastic_incoherency_selection_is_stored(ui):
    # Changed by the web_gui fix round: the original payload chose "Stochastic" with both SEEDs 0 and
    # expected it stored as stochastic.  D-INC-02 (normative) makes the input stochastic only when
    # HSeed or VSeed is non-zero and RandPhz > 0 (spec 05b: "arbitrary non-zero SEED integers"), so
    # SEEDs 0 + phase 180 is deterministic for the analysis; storing it would show "stochastic" for a
    # deterministic input.  The choice is now stored with a SEED, and a seedless choice is refused
    # with a message (never silently lost).  The 5.4 default phase 180 is still what is checked.
    v = D.values(ui)
    inc = dict(v["records"]["INCOH"], **{"@stoch": 1})
    with pytest.raises(D.DialogError):
        D.commit_commands(ui, {"records": {"HOUSE": dict(v["records"]["HOUSE"], coh=1), "INCOH": inc}})
    _commit_and_run(ui, {"records": {"HOUSE": dict(v["records"]["HOUSE"], coh=1), "INCOH": dict(inc, hseed=11, vseed=13)}})
    after = D.values(ui)["records"]["INCOH"]
    assert after["@stoch"] == 1
    assert after["randphz"] == 180


def _effective_freq2(interp):
    s = get_record(interp.model, "SITE")
    return int(s.freq2) if s.given(9) else max(int(s.nft) // 2, 1)          # AFWRITE rule (requirements 5.4)


def test_site_frequency2_follows_nfft(ui):
    _run_all(ui, ["SITE,0,1,0,20,0,1,0,1,,1,0,0.005,8192,1"])
    assert _effective_freq2(ui) == 4096
    v = D.values(ui)
    assert v["records"]["SITE"]["freq2"] in ("", 4096)                         # shown = used
    _commit_and_run(ui, {"records": {"SITE": dict(v["records"]["SITE"], nl=10)}})
    assert _effective_freq2(ui) == 4096                                       # unrelated edit keeps it


# ======================================================================================
# 4. Front-end command text (L17): every submitted command exists
# ======================================================================================
_CMD_CTX = re.compile(r"(?:S\.command|cmdOk|cmds\.push|lines\.push|calc)\(\s*(?:\(d, sl\) =>\s*)?\[?\s*"
                      r"[`\"]([A-Z][A-Z0-9]+)")


def test_front_end_submits_registered_commands_only():
    load_commands()
    names = set()
    for js in ("app.js", "dialogs.js", "plots.js"):
        text = (STATIC / js).read_text(encoding="utf-8")
        names.update(_CMD_CTX.findall(text))
        for arr in re.findall(r"S\.command\(\[([^\]]*)\]", text):
            names.update(re.findall(r"[`\"]([A-Z][A-Z0-9]+)", arr))
    assert {"ACTM", "SAVE", "INP", "MODELPLOT", "CUTPLOT", "CLOSEPLOT", "ACTIVATEPLOT", "CNGVIEW", "LINECOMBIN"} <= names
    missing = sorted(n for n in names if lookup(n) is None)
    assert not missing, missing


def test_browser_camera_uses_pan_and_centre():
    text = (STATIC / "plots.js").read_text(encoding="utf-8")
    body = re.search(r"function cameraOf\(d\) \{(.*?)\n  \}", text, re.S).group(1)
    assert re.search(r"\.px|extent|\.center", body), body


# ======================================================================================
# 5. Result files (requirements 1.8 naming formulas) and path safety
# ======================================================================================
@pytest.mark.parametrize("name, plot", [
    (nodal_result_name(12, 1, "TFU"), "spec"), (nodal_result_name(12, 6, "TFI"), "spec"),
    (nodal_rs_name(12, 3, 2), "spec"), (nodal_result_name(415, 2, "ACC"), "th"),
    (nodal_result_name(415, 2, "TFD"), "spec"), (nodal_result_name(415, 2, "THD"), "th"),
    (element_result_name("BEAMS", 3, 45, "MXJ", "THS"), "th"), (layer_th_name("ACC", 1), "th"),
    (layer_th_name("SN", 12), "th"), ("FILE8", None), ("run.sit", None),
])
def test_result_files_classified_by_naming_formula(name, plot):
    assert F.classify(Path(name))["plot"] == plot


def test_pair_history_file_info(tmp_path):
    p = tmp_path / "x.THD"
    p.write_text("0.0 1.0\n0.01 2.0\n0.02 3.0\n", encoding="utf-8")
    info = F.file_info(p)
    assert info["pair"] == 1 and info["points"] == 3


def test_safe_path_sibling_prefix_and_traversal(tmp_path):
    root = tmp_path / "model"
    sib = tmp_path / "model_other"
    root.mkdir()
    sib.mkdir()
    (sib / "x.txt").write_text("x", encoding="utf-8")
    (root / "a.txt").write_text("a", encoding="utf-8")
    roots = [root.resolve()]
    assert F.safe_path("a.txt", roots, must_exist=True) == (root / "a.txt").resolve()
    for bad in (str(sib / "x.txt"), "../model_other/x.txt", "sub/../../model_other/x.txt"):
        with pytest.raises(F.PathError):
            F.safe_path(bad, roots, must_exist=True)
    with pytest.raises(F.PathError):
        F.safe_path("a\x00.txt", roots)
    assert F.safe_path("new.pre", roots) == (root / "new.pre").resolve()


@pytest.mark.parametrize("host, ok", [
    ("127.0.0.1:8765", True), ("LOCALHOST:8765", True), ("[::1]:8765", True), ("localhost", True),
    ("localhost.evil.com", False), ("127.0.0.1.nip.io:8765", False), ("evil.com:127.0.0.1", False),
    ("[::2]:8765", False), ("0.0.0.0:8765", False),
])
def test_host_header_rule(host, ok):
    assert _host_ok(host, 8765) is ok


# ======================================================================================
# 6. Event log, token on every route, keep-alive
# ======================================================================================
def test_event_log_bounds_and_reset():
    log = EventLog(keep=5)
    for k in range(12):
        log.push("message", text=str(k))
    evs, last, reset = log.since(0)
    assert reset and [e["seq"] for e in evs] == [8, 9, 10, 11, 12] and last == 12
    evs, last, reset = log.since(10)
    assert not reset and [e["seq"] for e in evs] == [11, 12]
    t0 = time.time()
    evs, last, reset = log.since(12, timeout=0.2)
    assert evs == [] and last == 12 and time.time() - t0 >= 0.15

    def later():
        time.sleep(0.1)
        log.push("state")
    threading.Thread(target=later).start()
    evs, last, reset = log.since(12, timeout=5)
    assert [e["type"] for e in evs] == ["state"] and last == 13


@pytest.fixture
def server(tmp_path):
    srv = make_server(0, model_dir=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    th = threading.Thread(target=srv.serve_forever, kwargs={"poll_interval": 0.05}, daemon=True)
    th.start()
    yield srv
    srv.shutdown()
    srv.server_close()


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


def test_every_api_route_needs_the_token(server):
    import urllib.error
    import urllib.request
    for method, path in [("GET", "/api/state"), ("GET", "/api/model"), ("GET", "/api/files"), ("GET", "/api/file?name=x"),
                         ("POST", "/api/command"), ("POST", "/api/run/SITE"), ("POST", "/api/file"),
                         ("POST", "/api/settings/locations"), ("POST", "/api/db"), ("POST", "/api/exit"),
                         ("GET", "/api/events?since=0"), ("POST", "/api/verify")]:
        req = urllib.request.Request(server.url.rstrip("/") + path, method=method,
                                     data=b"{}" if method == "POST" else None,
                                     headers={"Content-Type": "application/json"})
        with pytest.raises(urllib.error.HTTPError) as ei:
            urllib.request.urlopen(req, timeout=10)
        assert ei.value.code == 403, path
    assert not server.session.interp.model.nodes and not server.session.jobs.jobs


def test_refused_post_does_not_desync_the_connection(server):
    body = b'{"line": "N,1,0,0,0"}'
    req1 = (b"POST /api/command HTTP/1.1\r\nHost: 127.0.0.1\r\nX-SASSI-Token: wrong\r\n"
            b"Content-Type: application/json\r\nContent-Length: %d\r\n\r\n" % len(body)) + body
    req2 = (b"GET /api/state HTTP/1.1\r\nHost: 127.0.0.1\r\nX-SASSI-Token: " + server.session.token.encode()
            + b"\r\nConnection: close\r\n\r\n")
    out = _raw(server, req1 + req2)
    status = re.findall(rb"HTTP/1\.[01] (\d{3})", out)
    assert status == [b"403", b"200"], out[:300]


# ======================================================================================
# 7. UI-04: module runs in a worker process (Cancel, progress) -- also when typed
# ======================================================================================
def _surface(S, tmp_path):
    r = S.execute([ln for ln in SURFACE.format(dir=tmp_path / "run").splitlines() if ln.strip()] + ["AFWRITE"])
    assert r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]


def test_typed_run_command_uses_a_worker_job(S, tmp_path):
    _surface(S, tmp_path)
    st, r = S.route("POST", "/api/command", {}, {"line": "RUNSITE"})
    assert st == 200
    assert any(j.line == "RUNSITE" for j in S.jobs.jobs.values())


def test_menu_run_job_snapshot_is_independent_of_later_edits(S, tmp_path):
    """The worker holds a copy of the models: a model edit after the start does not reach the run,
    and the run's model hash equals the AFWRITE hash (no spurious 'model changed' warning)."""
    _surface(S, tmp_path)
    st, job = S.route("POST", "/api/run/SITE", {}, {})
    assert st == 200
    S.execute(["N,99,5,5,5"])
    t0 = time.time()
    while S.jobs.get(job["id"]).active and time.time() - t0 < 120:
        time.sleep(0.05)
    d = S.jobs.get(job["id"]).to_dict()
    assert d["ok"], [m["text"] for m in d["messages"]]
    assert not any("model changed since the last AFWRITE" in m["text"] for m in d["messages"])


def test_load_model_open_with_comma_in_path(S, tmp_path):
    d = tmp_path / "Smith, J"
    d.mkdir()
    r = S.execute(S.db_open_lines("m", str(d)))
    assert r["ok"]
    assert Path(S.interp.model.path) == d


# ======================================================================================
# 8. Persistence (UI-07, D-UI-07, requirements 5.10)
# ======================================================================================
def test_check_options_and_toolbars_not_persisted(S, tmp_path):
    S.route("POST", "/api/options/CHECK", {}, {"records": {"CHECK": {"break_at": 7, "show_warnings": 0}}})
    S.exit()
    S2 = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    co = S2.route("GET", "/api/options/CHECK", {}, {})[1]["values"]["records"]["CHECK"]
    assert co["break_at"] == 100 and co["show_warnings"] == 1
    xml = (tmp_path / "settings" / "SASSIini.xml").read_text(encoding="utf-8")
    assert "break" not in xml.lower() and "toolbar" not in xml.lower()


def test_shader_options_persist(S, tmp_path):
    assert S.execute(["SHADEROPTIONS,12,0.03,0.1,2"])["ok"]
    S.exit()
    S2 = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    sh = plot_state(S2.interp).shader
    assert (sh.points, sh.linew, sh.shrink, sh.scale) == (12.0, 0.03, 0.1, 2.0)
