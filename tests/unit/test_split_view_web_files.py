"""The split view of the GUI and the files of the browser version (requirements 7.20, D-W6-07 to D-W6-10).

* Split view: plots and the Results browser open in a tab group of their own beside the work (Command
  History, Learn, the lesson, File Editors, module output); View > Split View switches it off.
* Browser version: the workspace and the settings are IndexedDB file systems kept between visits (File >
  Clear Saved Files deletes them); every file dialog copies files from the visitor's computer ("From your
  computer..."); Model > Output and Export Table download what they write.

The behaviour was checked in a browser (local GUI and a local build of the site); these tests pin the
contracts between the scripts."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"
WEB = ROOT / "web"


def _read(p: Path) -> str:
    return p.read_text(encoding="utf-8")


def test_two_tab_groups_in_the_page():
    page = _read(STATIC / "index.html")
    for ident in ('id="group-main"', 'id="tabbar"', 'id="tabcontent"', 'id="split-resizer"', 'id="group-view"',
                  'id="tabbar2"', 'id="tabcontent2"'):
        assert ident in page, ident
    assert page.index('id="group-main"') < page.index('id="split-resizer"') < page.index('id="group-view"')
    css = _read(STATIC / "styles.css")
    assert "#group-view[hidden], #split-resizer[hidden] { display: none; }" in css
    assert ".tabcontent { flex: 1; position: relative; min-height: 0; }" in css
    assert re.search(r"@media \(max-width: 760px\) \{\s+/\* a narrow window: the plots and results below", css)


def test_tab_groups_in_app_js():
    js = _read(STATIC / "app.js")
    # plots, the Results browser and the Run view go right while the view is split; everything else stays in the work group
    assert 'S.groupFor = (t) => (S.split && (t.kind === "plot" || t.kind === "results" || t.kind === "runview") ? "view" : "main");' in js
    assert "S.isShown = (t) =>" in js and "S.layoutSplit = function" in js and "S.setSplit = function" in js
    assert '{label: "Split View", check: () => S.split, action: () => S.setSplit(!S.split)' in js
    assert 'lsGet(SPLIT_KEY) !== "0"' in js                       # on by default, kept by the browser
    # selecting a tab only touches its own group; removing the last right tab closes the group
    sel = js[js.index("S.selectTab = function"):js.index("S.requestCloseTab = function")]
    assert "if (o.group !== t.group) continue;" in sel and "S.groups[t.group].active = t;" in sel
    rem = js[js.index("S.removeTab = function"):js.index("// ------------------------------------------------------------------ status bar")]
    assert "S.layoutSplit();" in rem and "const same = S.tabs.filter((o) => o.group === t.group);" in rem
    # the 3D plot keys act after a click in the plot group
    assert 'S.pointerGroup = g.id === "group-view" ? "view" : "main";' in js
    assert "setupSplitResizer();" in js


def test_plots_and_lessons_use_the_group_visibility():
    plots = _read(STATIC / "plots.js")
    assert "S.activeTab !== t" not in plots and "S.activeTab === t" not in plots
    assert plots.count("S.isShown(t)") >= 5
    assert 'if (!S.keepFocus || t.group === "view") S.selectTab(t.id);' in plots     # a new plot beside the lesson
    assert 'S.split ? (S.pointerGroup === "view" ? S.groups.view.active : null) : S.activeTab' in plots
    learn = _read(STATIC / "learn.js")
    assert "const front = S.groups ? S.groups.main.active : S.activeTab;" in learn
    assert 'S.$("#group-main") || S.$("#workarea")' in learn


def test_browser_workspace_is_kept_in_indexeddb():
    w = _read(WEB / "worker.js")
    assert 'const WORK = "/home/pyodide/work";' in w and 'const SETTINGS = "/home/pyodide/.sassi-edu";' in w
    assert "FS.mount(FS.filesystems.IDBFS, {}, d);" in w
    assert "FS.syncfs(true," in w and "pyFS.syncfs(false," in w
    # restored before the session (the bridge) starts; saved after a non-GET request, a job, and on flush
    assert w.index("await restoreKept(py)") < w.index('py.pyimport("sassi.web.bridge")')
    assert 'if (m.method !== "GET") saveSoon();' in w and 'if (m.type === "flush") { save(); return; }' in w
    b = _read(WEB / "boot.js")
    dbs = re.search(r"const KEPT_DBS = \[([^\]]+)\]", b).group(1)
    assert '"/home/pyodide/work"' in dbs and '"/home/pyodide/.sassi-edu"' in dbs    # the worker's mount points
    assert "indexedDB.deleteDatabase(name)" in b and "worker.terminate()" in b
    assert 'window.addEventListener("pagehide", flush);' in b
    from sassi.web import bridge
    assert bridge.WORK_DIR == "/home/pyodide/work" and bridge.SETTINGS_DIR == "/home/pyodide/.sassi-edu"


def test_files_of_the_visitors_computer_and_downloads():
    d = _read(STATIC / "dialogs.js")
    pick = d[d.index("D.pickFile = function"):d.index("D.pickFile = function") + 6000]
    assert "if (S.web && !opt.folder && !opt.noUpload)" in pick and "From your computer…" in pick
    assert "await uploadTexts(chosen, dir);" in pick
    assert "async function uploadTexts(files, folder)" in d and "async function downloadWorkspaceFile(path, what)" in d
    out = d[d.index("D.output = function"):d.index("D.exportAnsys = async function")]
    assert 'if (S.web && w) await downloadWorkspaceFile(w[1].trim(), "Output");' in out
    assert 'if (S.web) await downloadWorkspaceFile(r.path, "Export Table");' in d
    assert "D.clearSavedFiles = async function" in d and "web.clearSaved" in d
    app = _read(STATIC / "app.js")
    assert '{label: "Clear Saved Files...", action: () => D().clearSavedFiles(), webOnly: true' in app


def test_built_site_keeps_the_layout():
    """web/build.py writes the GUI's page with both tab groups and the boot script of the browser version."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("webbuild", WEB / "build.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    html = mod.page("testbuild")
    assert 'id="group-main"' in html and 'id="group-view"' in html and 'id="split-resizer"' in html
    assert "web/boot.js" in html
