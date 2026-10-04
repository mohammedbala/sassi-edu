"""Adversarial review of work package 'gui_polish' (requirements 5.3-5.9, spec 06 section 5.2).

The implementer of this package changed no file (its report says the run was spent on a disk survey),
so these tests probe the defects the package was meant to fix, from a different angle than the
existing tests/unit/test_ui_*.py:

* backend contract tests (pass): the JSON API already carries everything the Line Selection dialog
  needs to pre-fill itself from the active plot, and READSPEC publishes a "lines" event while the
  command is still running (the ingredient of the start-number race in the dialog);
* front-end probes (strict xfail until the package is done): the dialog pre-fills from the active
  plot, reads the starting number before awaiting the command, its axis rows fit their grid column,
  and Help serves the docs/ Markdown rendered to HTML.

Static JS/CSS probes are used because no headless browser is part of the test environment; each one
states the numbers or source lines it checks.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from sassi.ui.api import GuiSession

STATIC = Path(__file__).resolve().parents[2] / "sassi" / "ui" / "static"


@pytest.fixture
def S(tmp_path):
    s = GuiSession(cwd=str(tmp_path), settings_dir=str(tmp_path / "settings"))
    yield s
    run = s.jobs.running()
    if run is not None:
        s.jobs.cancel(run.id)


def call(S, method, path, body=None, query=None):
    q = {k: [str(v)] for k, v in (query or {}).items()}
    return S.route(method, path, q, body or {})


def _two_spectra(tmp_path):
    a = tmp_path / "00001TR_X.TFU"
    a.write_text("# a\n0.5 1.0\n1.0 2.0\n2.0 1.5\n", encoding="utf-8")
    b = tmp_path / "00041TR_X01.RS"
    b.write_text("# b\n0.5 0.2\n1.0 0.4\n2.0 0.3\n", encoding="utf-8")
    return a, b


def _line_selection_source() -> str:
    js = (STATIC / "dialogs.js").read_text(encoding="utf-8")
    i = js.index("D.lineSelection = async function")
    j = js.index("D.onLinesChanged", i)
    return js[i:j]


# ------------------------------------------------------------------ backend contract (should pass)
def test_active_plot_record_carries_everything_needed_to_prefill(S, tmp_path):
    """Spec 06 section 5.2 + work-package item (1): reopening Line Selection must show the active
    plot's lines, titles, log axes and ranges. The API already exposes all of them, so the fix is
    front-end only."""
    a, b = _two_spectra(tmp_path)
    st, r = call(S, "POST", "/api/command", {"lines": [
        f"READSPEC,{a.name},1,1", f"READSPEC,{b.name},1,2", "SPECPLOT,1,2", "PLOTTITLE,My RS",
        "XTITLE,f", "YTITLE,Sa", "AXES,1,1,0,0,1,0", "PLOTRANGE,0.5,2,0,3"]})
    assert st == 200 and r["ok"], [m["text"] for m in r["messages"] if m["kind"] == "ERROR"]
    st, plots = call(S, "GET", "/api/plots")
    act = [p for p in plots["plots"] if p["id"] == plots["active"]][0]
    assert act["kind"] == "SPECPLOT" and act["params"]["lines"] == [1, 2]
    s = act["settings"]
    assert s["log_x"] is True and s["log_y"] is False
    assert (s["xtitle"], s["ytitle"]) == ("f", "Sa")
    assert (s["xmin"], s["xmax"], s["ymin"], s["ymax"]) == (0.5, 2.0, 0.0, 3.0)
    assert "My RS" in (act["title"] or "")


def test_readspec_publishes_lines_event_during_the_command(S, tmp_path):
    """The browser long-polls /api/events; READSPEC publishes a "lines" event while POST
    /api/command is still in flight. The dialog's onLinesChanged -> refresh() then rewrites the
    'Starting Number' box to max+1, which addLines() reads only after awaiting the command
    (dialogs.js addLines). This is the ingredient of the 'new line not ticked' race."""
    a, b = _two_spectra(tmp_path)
    call(S, "POST", "/api/command", {"line": f"READSPEC,{a.name},1,1"})
    st, ev0 = call(S, "GET", "/api/events", query={"since": 0})
    since = max([e.get("seq", e.get("id", 0)) for e in ev0["events"]] or [0])
    call(S, "POST", "/api/command", {"line": f"READSPEC,{b.name},1,2"})
    ev = call(S, "GET", "/api/events", query={"since": since})[1]["events"]
    lines_ev = [e for e in ev if e.get("type") == "plot" and e.get("event") == "lines"]
    assert lines_ev and lines_ev[-1]["data"]["numbers"] == [2]


# ------------------------------------------------------------------ front-end probes (package not done)
def test_line_selection_prefills_from_active_plot():
    src = _line_selection_source()
    assert "plots.active" in src or "activePlot" in src
    assert "log_x" in src and "log_y" in src
    assert "xtitle" in src and "params" in src


def test_add_lines_captures_start_number_before_awaiting():
    src = _line_selection_source()
    add = src[src.index("const addLines"):src.index("const body")]
    i_read = add.index("Number(start.value)")
    i_await = add.index("await S.command")
    assert i_read < i_await


def _css_px(css: str, selector: str, prop: str) -> float:
    m = re.search(re.escape(selector) + r"\s*\{[^}]*?" + re.escape(prop) + r"\s*:\s*([0-9.]+)px", css)
    assert m, (selector, prop)
    return float(m.group(1))


def test_line_selection_axis_row_fits_its_grid_column():
    """Layout arithmetic from styles.css and dialogs.js (no browser): the right grid column of a
    640 px body with an 8 px gap is 316 px; the row 'Min [num] Max [num]' needs
    label.lab min-width + 2 x input.num basis + 3 gaps + the 'Max' label + fieldset padding."""
    css = (STATIC / "styles.css").read_text(encoding="utf-8")
    src = _line_selection_source()
    m = re.search(r'const body = el\("div", \{style: \{width: "(\d+)px"\}\}', src)
    body_w = float(m.group(1)) if m else 0.0
    lab = _css_px(css, ".row > label.lab", "min-width")
    num = float(re.search(r"\.row input\.num\s*\{\s*flex:\s*0 0 (\d+)px", css).group(1))
    gap = _css_px(css, ".row", "gap")
    max_label = 28.0                     # 'Max' at 13 px UI font, conservative
    fieldset_pad = 2 * 8 + 2             # padding 4px 8px + 1 px borders
    need = lab + 2 * num + 3 * gap + max_label + fieldset_pad
    column = (body_w - 8) / 2 if body_w else float("inf")
    assert need <= column, f"row needs {need:.0f} px, column is {column:.0f} px"


def test_help_serves_docs_rendered_to_html(S):
    st, h = call(S, "GET", "/api/help")
    assert st == 200
    blob = repr(h)
    assert "<h1" in blob or "<h2" in blob
    assert "impedance_study" in blob or "verification" in str(h.get("docs", "")).lower()


def test_gui_guide_has_a_section_per_top_level_menu():
    """Work package: docs/user/GUI.md keeps a menu-by-menu description (requirements 5.3)."""
    guide = (Path(__file__).resolve().parents[2] / "docs" / "user" / "GUI.md").read_text(encoding="utf-8")
    for menu in ("Model", "File", "Plot", "Modules", "Options", "View", "Help"):
        assert re.search(r"\b" + menu + r"\b", guide), menu
