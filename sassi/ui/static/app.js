/* SASSI-EDU GUI -- main window (requirements 5.2): menu bar, toolbars, tabbed document area,
 * Command Entry, status bar, event polling.  Rule L17: every action that changes the model or the
 * plots submits command text to the session interpreter (POST /api/command); the Command History
 * therefore replays as a .pre file.  Vanilla JavaScript, no build step.
 *
 * Browser version (GitHub Pages): web/boot.js sets window.SASSI_TRANSPORT before this file runs; the
 * requests then go to the Python session in a Web Worker (Pyodide) instead of the local server, the event
 * poll never long-polls, and the worker pushes the events and the job output while it is busy (S.push). */
"use strict";

const SASSI = (() => {
  const S = {};
  S.token = (document.querySelector('meta[name="sassi-token"]') || {}).content || "";
  // browser version (web/boot.js): requests go to Python in a Web Worker, there is no server
  S.transport = window.SASSI_TRANSPORT || null;
  S.web = !!S.transport;
  S.state = null;            // GET /api/state
  S.lastEvent = 0;
  S.tabs = [];               // {id, title, kind, pane, closable, plotId, group, ...}
  S.activeTab = null;        // the tab brought to the front last (either group)
  S.recall = [];             // Command Entry Up/Down recall
  S.recallPos = 0;
  S.queue = Promise.resolve();
  S.busy = 0;
  S.connectedEditor = null;  // File Editor connected to Command Entry (one at a time)
  S.toolbarsVisible = {main: true, plot: true};   // not persisted (UI-07)
  const $ = (sel, root) => (root || document).querySelector(sel);
  S.$ = $;

  // ------------------------------------------------------------------ helpers
  S.el = function (tag, attrs, ...children) {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v === undefined || v === null || v === false) continue;
      if (k === "class") e.className = v;
      else if (k === "text") e.textContent = v;
      else if (k === "html") e.innerHTML = v;
      else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
      else if (k === "style" && typeof v === "object") Object.assign(e.style, v);
      else if (v === true) e.setAttribute(k, "");
      else e.setAttribute(k, v);
    }
    for (const c of children.flat()) {
      if (c === null || c === undefined || c === false) continue;
      e.appendChild(typeof c === "string" || typeof c === "number" ? document.createTextNode(String(c)) : c);
    }
    return e;
  };
  const el = S.el;
  S.sleep = (ms) => new Promise((r) => setTimeout(r, ms));
  S.esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}[c]));
  S.fmt = function (v, d) {
    if (v === null || v === undefined || v === "") return "";
    if (typeof v !== "number") return String(v);
    if (Number.isInteger(v)) return String(v);
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e9)) return v.toExponential(Math.max((d || 6) - 1, 1)).replace(/\.?0+e/, "e");
    return Number(v.toPrecision(d || 6)).toString();
  };

  // ------------------------------------------------------------------ API
  S.api = async function (method, path, body) {
    if (S.transport) return webApi(method, path, body);
    const opt = {method, headers: {"X-SASSI-Token": S.token}};
    if (body !== undefined) {
      opt.headers["Content-Type"] = "application/json";
      opt.body = JSON.stringify(body);
    }
    const r = await fetch(path, opt);
    let data;
    try { data = await r.json(); } catch (e) { data = {error: r.statusText || "bad response"}; }
    if (!r.ok) {
      const err = new Error((data && data.error) || r.statusText);
      err.status = r.status;
      err.data = data;
      throw err;
    }
    return data;
  };
  /** The browser version: the transport (web/boot.js) answers {status, payload} from the Python session in the
   *  Web Worker; an error status is raised as for HTTP. */
  async function webApi(method, path, body) {
    const r = await S.transport(method, path, body);
    const data = r.payload;
    if (r.status >= 400) {
      const err = new Error((data && data.error) || `error ${r.status}`);
      err.status = r.status;
      err.data = data;
      throw err;
    }
    return data;
  }
  S.get = (path) => S.api("GET", path);
  S.post = (path, body) => S.api("POST", path, body || {});

  // ------------------------------------------------------------------ commands (L17)
  /** Submit command text (one line or several) in order; resolves with the API result. */
  S.command = function (lines, opts) {
    opts = opts || {};
    if (typeof lines === "string") lines = [lines];
    lines = lines.filter((l) => l !== null && l !== undefined && String(l).trim() !== "");
    if (!lines.length) return Promise.resolve({ok: true, results: [], messages: []});
    const p = S.queue.then(() => runCommand(lines, opts));
    S.queue = p.catch(() => {});
    return p;
  };
  async function runCommand(lines, opts) {
    setBusy(1);
    try {
      const res = await S.post("/api/command", {lines});
      // Input > Connect to Command Entry: every accepted command goes into the connected editor --
      // except the INP of its own Run button (the file would INP itself on the next run)
      if (S.connectedEditor && S.connectedEditor.tabId !== opts.skipEditor) {
        for (const r of res.results) if (r.ok) S.connectedEditor.append(r.line);
      }
      return res;
    } catch (e) {
      S.local("ERROR", `${lines.join(" | ")}: ${e.message}`);
      throw e;
    } finally {
      setBusy(-1);
    }
  }
  function setBusy(d) {
    S.busy += d;
    $("#cmd").classList.toggle("busy", S.busy > 0);
    if (S.busy === 0 && S.activity) S.activity.onIdle();
    if (S.busy > 0) S.status("Running ...");
    else if (!S.state || !S.state.job) {
      S.progress(null);
      if ($("#status-text").textContent === "Running ...") S.status("Ready");
    }
  }

  // ------------------------------------------------------------------ Command History
  const HIST_MAX = 30000;
  const PREFIX = {ECHO: "> ", WARNING: "*** WARNING: ", ERROR: "*** ERROR: "};
  S.historyEl = null;
  /* Lines are appended in batches: a module run streams thousands of listing lines, and measuring the scroll
   * position for every line forces a layout each time (seconds per batch on a long history, during which the
   * page cannot repaint -- the activity panel froze).  A batch is flushed once the current events have been
   * handled (a microtask): one measure, one append, one scroll. */
  let pendingLines = null;
  S.appendMessage = function (kind, text) {
    if (!S.historyEl) return;
    const div = document.createElement("div");
    div.className = "m-" + kind;
    div.textContent = (PREFIX[kind] || "") + text;
    if (kind === "ECHO") { div.dataset.line = text; div.title = "click: explain this command"; }
    if (!pendingLines) {
      pendingLines = document.createDocumentFragment();
      queueMicrotask(flushMessages);
    }
    pendingLines.appendChild(div);
  };
  function flushMessages() {
    const h = S.historyEl, frag = pendingLines;
    pendingLines = null;
    if (!h || !frag) return;
    const sc = h.parentElement;
    const atEnd = sc.clientHeight ? sc.scrollTop + sc.clientHeight >= sc.scrollHeight - 30 : (!S.historyAtEnd || S.historyAtEnd());
    h.appendChild(frag);
    const extra = h.childNodes.length - HIST_MAX;
    for (let i = 0; i < extra; i++) h.removeChild(h.firstChild);
    if (atEnd) sc.scrollTop = sc.scrollHeight;
  }
  /** A GUI-side message (not from the interpreter). */
  S.local = function (kind, text) { S.appendMessage(kind === "ERROR" || kind === "WARNING" ? kind : "LOCAL", text); };
  S.applyDisplayFilters = function () {
    const d = (S.state && S.state.settings && S.state.settings.display) || {};
    const h = S.historyEl;
    if (!h) return;
    h.classList.toggle("hide-ECHO", d.echo === false);
    h.classList.toggle("hide-CONFIRM", d.confirm === false);
    h.classList.toggle("hide-COMMENT", d.comments === false);
    h.classList.toggle("hide-WARNING", d.warnerr === false);
    h.classList.toggle("hide-ERROR", d.warnerr === false);
  };
  S.applyColours = function () {
    const c = (S.state && S.state.settings && S.state.settings.colours) || {};
    const map = {ECHO: "--c-echo", CONFIRM: "--c-confirm", COMMENT: "--c-comment", INFO: "--c-info",
                 WARNING: "--c-warning", ERROR: "--c-error"};
    for (const [k, v] of Object.entries(map)) if (c[k]) document.documentElement.style.setProperty(v, c[k]);
  };
  function makeHistoryTab() {
    const body = el("div", {class: "pane-body"});
    const hist = el("div", {class: "history", "aria-live": "polite"});
    body.appendChild(hist);
    const pane = el("div", {class: "pane"},
      el("div", {class: "pane-bar"},
        el("span", {class: "grow", text: "Messages of the command interpreter: echo, confirmation, comments, information, warnings, errors"}),
        el("button", {class: "btn small", text: "Copy", onclick: () => navigator.clipboard && navigator.clipboard.writeText(hist.innerText)}),
        el("button", {class: "btn small", text: "Clear", onclick: () => { hist.innerHTML = ""; }})),
      body);
    S.historyEl = hist;
    // the history follows new messages while it is scrolled to its end, also while hidden behind
    // another tab: bringing it back shows the latest messages
    let atEnd = true;
    body.addEventListener("scroll", () => { if (body.clientHeight) atEnd = body.scrollTop + body.clientHeight >= body.scrollHeight - 30; });
    const t = S.addTab({id: "history", title: "Command History", kind: "history", pane, closable: false});
    t.onActivate = () => { if (atEnd) body.scrollTop = body.scrollHeight; };
    S.historyAtEnd = () => atEnd;
    return t;
  }
  /** View > Command Window: show (re)open the Command History; reopening clears its output. */
  S.commandWindow = function () {
    if (S.tabs.find((t) => t.id === "history")) {
      if (S.historyEl) S.historyEl.innerHTML = "";
      S.selectTab("history");
      return;
    }
    makeHistoryTab();
    S.applyDisplayFilters();
    S.selectTab("history");
  };

  // ------------------------------------------------------------------ tabs (two groups: split view)
  /** Split view (View > Split View, on by default): plots and the Results browser open in a group of their
   *  own beside the work -- Command History, Learn, the lesson, File Editors, module output -- instead of a tab
   *  in front of it; a narrow window stacks the two groups.  Off: one group, as in ACS SASSI. */
  const SPLIT_KEY = "sassi-edu.split", SPLIT_W_KEY = "sassi-edu.splitWidth";
  const lsGet = (k) => { try { return window.localStorage.getItem(k); } catch (e) { return null; } };
  const lsSet = (k, v) => { try { window.localStorage.setItem(k, v); } catch (e) { /* private mode */ } };
  S.split = lsGet(SPLIT_KEY) !== "0";
  S.groups = {main: {id: "main", active: null}, view: {id: "view", active: null}};
  const groupEls = (g) => g === "view" ? {bar: $("#tabbar2"), content: $("#tabcontent2")} : {bar: $("#tabbar"), content: $("#tabcontent")};
  /** The group of a tab: plots and results right while the view is split, everything else left. */
  S.groupFor = (t) => (S.split && (t.kind === "plot" || t.kind === "results" || t.kind === "runview") ? "view" : "main");
  /** Is tab t on the screen (the front tab of its group)? */
  S.isShown = (t) => !!t && S.groups[t.group] && S.groups[t.group].active === t;
  function placeTab(t) {
    t.group = S.groupFor(t);
    const g = groupEls(t.group);
    g.bar.appendChild(t.head);
    g.content.appendChild(t.pane);
  }
  /** Show the right group while it holds a tab; its width (a fraction of the work area) is kept per browser. */
  S.layoutSplit = function () {
    const view = $("#group-view"), rz = $("#split-resizer");
    if (!view) return;
    const any = S.tabs.some((t) => t.group === "view");
    const was = !view.hidden;
    view.hidden = !any;
    if (rz) rz.hidden = !any;
    const f = Number(lsGet(SPLIT_W_KEY)) || 0.5;
    view.style.width = `${Math.round(Math.max(0.2, Math.min(0.8, f)) * 100)}%`;
    if (was !== any) window.dispatchEvent(new Event("resize"));      // the plots follow their new width
  };
  S.setSplit = function (on) {
    S.split = !!on;
    lsSet(SPLIT_KEY, on ? "1" : "0");
    const front = S.activeTab;
    for (const t of S.tabs) placeTab(t);
    // one front tab per group (selectTab marks it and unmarks the rest of its group)
    for (const g of [S.groups.view, S.groups.main]) {
      const cands = S.tabs.filter((t) => t.group === g.id);
      const keep = g.active && g.active.group === g.id ? g.active : (cands.includes(front) ? front : cands[cands.length - 1]);
      g.active = null;
      if (keep) S.selectTab(keep.id);
    }
    if (front && S.tab(front.id) && !S.isShown(front)) S.selectTab(front.id);
    S.layoutSplit();
    window.dispatchEvent(new Event("resize"));
    if (S.rebuildMenus) S.rebuildMenus();
  };
  function setupSplitResizer() {
    const rz = $("#split-resizer"), wa = $("#workarea"), view = $("#group-view");
    if (!rz || !wa || !view) return;
    rz.addEventListener("mousedown", (ev) => {
      ev.preventDefault();
      const box = wa.getBoundingClientRect();
      const move = (e) => {
        const w = Math.max(260, Math.min(box.right - e.clientX, box.width - 260));
        view.style.width = `${Math.round(100 * w / box.width)}%`;
      };
      const up = () => {
        document.removeEventListener("mousemove", move);
        document.removeEventListener("mouseup", up);
        document.body.classList.remove("resizing");
        lsSet(SPLIT_W_KEY, String(view.getBoundingClientRect().width / Math.max(1, box.width)));
        window.dispatchEvent(new Event("resize"));
      };
      document.body.classList.add("resizing");
      document.addEventListener("mousemove", move);
      document.addEventListener("mouseup", up);
    });
    rz.addEventListener("dblclick", () => { lsSet(SPLIT_W_KEY, "0.5"); S.layoutSplit(); window.dispatchEvent(new Event("resize")); });
  }
  S.addTab = function (t) {
    t.closable = t.closable !== false;
    const head = el("div", {class: "tab", role: "tab", title: t.title},
      el("span", {class: "t", text: t.title}),
      t.closable ? el("button", {class: "x", title: "Close", text: "×", onclick: (ev) => { ev.stopPropagation(); S.requestCloseTab(t.id); }}) : null);
    head.addEventListener("click", () => S.selectTab(t.id, {user: true}));
    head.addEventListener("auxclick", (ev) => { if (ev.button === 1 && t.closable) S.requestCloseTab(t.id); });
    t.head = head;
    placeTab(t);
    S.tabs.push(t);
    S.layoutSplit();
    return t;
  };
  S.tab = (id) => S.tabs.find((t) => t.id === id);
  S.setTabTitle = function (id, title) {
    const t = S.tab(id);
    if (!t) return;
    t.title = title;
    t.head.querySelector(".t").textContent = title;
    t.head.title = title;
  };
  /** Bring a tab to the front of its group (the other group keeps its front tab on the screen).  A plot tab
   *  chosen by the user submits ACTIVATEPLOT,<id> (L17). */
  S.selectTab = function (id, opts) {
    opts = opts || {};
    const t = S.tab(id);
    if (!t) return;
    for (const o of S.tabs) {
      if (o.group !== t.group) continue;
      o.head.classList.toggle("active", o === t);
      o.pane.classList.toggle("active", o === t);
    }
    S.groups[t.group].active = t;
    S.activeTab = t;
    if (t.head.scrollIntoView) t.head.scrollIntoView({block: "nearest", inline: "nearest"});   // a crowded tab bar
    if (opts.user && t.plotId && S.state && S.state.plots && S.state.plots.active !== t.plotId) {
      S.command(`ACTIVATEPLOT,${t.plotId}`);
    }
    if (t.onActivate) t.onActivate(opts);
    S.updateToolbars();
  };
  S.requestCloseTab = function (id) {
    const t = S.tab(id);
    if (!t) return;
    if (t.plotId) {             // plots are closed by command (CLOSEPLOT closes the active plot)
      const lines = [];
      if (S.state.plots.active !== t.plotId) lines.push(`ACTIVATEPLOT,${t.plotId}`);
      lines.push("CLOSEPLOT");
      S.command(lines);
      return;
    }
    if (t.beforeClose) {
      const ok = t.beforeClose();             // true / false, or a promise of it (a question dialog)
      if (ok && typeof ok.then === "function") { ok.then((yes) => { if (yes) S.removeTab(id); }); return; }
      if (ok === false) return;
    }
    S.removeTab(id);
  };
  S.removeTab = function (id) {
    const i = S.tabs.findIndex((t) => t.id === id);
    if (i < 0) return;
    const t = S.tabs[i];
    if (t.onClose) t.onClose();
    t.head.remove();
    t.pane.remove();
    S.tabs.splice(i, 1);
    if (S.connectedEditor && S.connectedEditor.tabId === id) S.connectedEditor = null;
    const g = S.groups[t.group];
    if (g && g.active === t) {
      // the next tab of the same group comes to the front; an emptied right group closes
      g.active = null;
      const same = S.tabs.filter((o) => o.group === t.group);
      const next = same.length ? same[Math.min(S.tabs.slice(0, i).filter((o) => o.group === t.group).length, same.length - 1)] : null;
      if (next) S.selectTab(next.id);
    }
    if (S.activeTab === t) S.activeTab = S.groups.main.active || S.groups.view.active || null;
    S.layoutSplit();
  };

  // ------------------------------------------------------------------ status bar
  S.status = function (text) { $("#status-text").textContent = text; };
  S.progress = function (fraction, text) {
    const p = $("#status-progress");
    if (fraction === null || fraction === undefined) { p.classList.remove("on"); return; }
    p.classList.add("on");
    p.value = Math.max(0, Math.min(1, fraction));
    if (text) S.status(text);
  };
  function updateTitle() {
    const st = S.state;
    if (!st) return;
    const m = st.models.find((x) => x.active) || {};
    const name = m.name ? `${m.name}` : "(no name)";
    $("#model-badge").textContent = `Model ${st.active}: ${name}${m.path ? "  —  " + m.path : ""}`;
    $("#status-model").textContent = `Model ${st.active} | ${m.nodes || 0} nodes, ${m.elements || 0} elements | cwd ${st.cwd}`;
    document.title = `SASSI-EDU — Model ${st.active} ${m.name || ""}`;
  }
  S.updateTitle = updateTitle;

  // ------------------------------------------------------------------ menus (requirements 5.3)
  S.lowestUnusedModel = function () {
    const used = new Set((S.state.models || []).map((m) => m.number));
    let n = 0;
    while (used.has(n)) n++;
    return n;
  };
  const D = () => SASSI.D;     // dialogs.js
  const P = () => SASSI.P;     // plots.js
  function menuDefs() {
    const runItems = ["EQUAKE", "SOIL", "LIQUEF", "SITE", "POINT", "HOUSE", "PINT", "FORCE", "ANALYS", "COMBIN", "MOTION", "STRESS", "RELDISP"]
      .map((m) => (m === "LIQUEF" || m === "PINT") ? {label: m, disabled: true, tag: "not in V3"}
        : {label: m, action: () => S.runModule(m), tip: `RUN${m} on the active model`});
    return [
      {title: "Model", items: [
        {label: "New", action: () => S.command(`ACTM,${S.lowestUnusedModel()}`), tip: "ACTM,<lowest unused number>"},
        {label: "Open", key: "Ctrl+O", action: () => D().loadModel()},
        {label: "Save", action: () => S.command("SAVE")},
        "-",
        {label: "Input", action: () => D().pickFile({title: "Input (.pre)", filter: [".pre", ".txt", ".mac"], onOk: (p) => S.command(`INP,${S.q(p)}`)})},
        {label: "Open Example...", action: () => SASSI.Learn.examplesDialog(), tag: "Learn", tip: "copy a tutorial example into a fresh workspace: open it in the File Editor or run it"},
        {label: "Converters", sub: [
          {label: "SASSI .hou", action: () => D().converter("SSI")},
          {label: "ANSYS .cdb", action: () => D().converter("ANSYS")},
          {label: "GT-STRUDL Database", disabled: true, tag: "not included"}]},
        "-",
        {label: "Output", action: () => D().output()},
        {label: "Export to ANSYS", action: () => D().exportAnsys()},
        {label: "Export to STRUDL", disabled: true, tag: "not included"},
        "-",
        {label: "Exit", action: () => D().exit(), web: false}]},
      {title: "File", items: [
        {label: "Open", action: () => D().pickFile({title: "File > Open (text editor)", create: true, onOk: (p) => S.openEditor(p)})},
        // browser version: the workspace is in the memory of this tab -- copy files in and out
        {label: "Upload to Workspace...", action: () => D().uploadFiles(), webOnly: true, tip: "copy text files from your computer into the workspace of this browser tab"},
        {label: "Download...", action: () => D().downloadFile(), webOnly: true, tip: "save a text file of the workspace (listing, deck, .pre, spectrum ...) on your computer"},
        {label: "Clear Saved Files...", action: () => D().clearSavedFiles(), webOnly: true, tip: "delete the workspace and settings this browser keeps between visits, and start afresh"},
        {label: "Export Image", action: () => D().exportImage(), enabled: () => !!(S.state && S.state.plots.active)},
        {label: "Export Table", action: () => D().exportTable(), enabled: () => P().activeKind() && ["SPECPLOT", "THPLOT", "SOILPROPPLOT"].includes(P().activeKind())},
        "-",
        {label: "Results Browser", action: () => S.openResults(), tag: "extension"}]},
      {title: "Plot", items: [
        {label: "Model", sub: [
          {label: "Elements", action: () => S.command("MODELPLOT")},
          {label: "Nodes", action: () => S.command("NODEPLOT")}]},
        {label: "Cuts", action: () => D().cutPlot()},
        {label: "Spectrum TFU-TFI", action: () => D().lineSelection("SPECPLOT")},
        {label: "Time History", action: () => D().lineSelection("THPLOT")},
        {label: "Soil Layers", action: () => S.command("LAYERPLOT")},
        {label: "Soil Properties", action: () => D().soilProperty()},
        {label: "Non Uniform Soil Field", disabled: true, tag: "not active"},
        "-",
        {label: "Process Animation Frame List", action: () => D().procFrame()},
        "-",
        {label: "Bubble", action: () => D().loadFrameData("BUBBLEPLOT")},
        {label: "Vector", action: () => D().loadFrameData("VECTORPLOT")},
        {label: "Contour", action: () => D().loadFrameData("CONTOURPLOT")},
        {label: "Deformed Shape", action: () => D().loadFrameData("DEFORMPLOT")}]},
      {title: "Modules", items: [
        {label: "Location", action: () => D().moduleLocation(), web: false},
        {label: "Extension", action: () => D().moduleExtension()},
        "-",
        ...runItems,
        // Option NON (tier P2): enabled when the server reports RUNNONLINEAR as implemented
        {label: "NONLINEAR", action: () => S.runModule("NONLINEAR"), tag: "P2",
         enabled: () => !!(S.state && (S.state.run_modules || []).includes("NONLINEAR")), tip: "RUNNONLINEAR on the active model"},
        "-",
        // Option A (LOADGEN, P2): the converter dialogs; Ok stores the settings, Run also RUNLOADGEN (a worker job)
        {label: "ANSYS Eq. Static Load", action: () => D().optionsDialog("LOADGEN"), tag: "P2",
         tip: "ANSYS Static Load Converter: LOADGEN, LGFILE, LGNODE, LGTIME ...; Run = RUNLOADGEN,STATIC"},
        {label: "ANSYS Dynamic Load", action: () => D().optionsDialog("LOADGENDYN"), tag: "P2",
         tip: "ANSYS Dynamic Load Converter: LOADGENDYN, LGFILE, LGNODE ...; Run = RUNLOADGEN,DYNAMIC"},
        {label: "ANSYS Super Element Utilities", disabled: true, tag: "P2"}]},
      {title: "Options", items: [
        {label: "Model", action: () => D().optionsDialog("MODEL")},
        {label: "Write", action: () => D().optionsDialog("WRITE")},
        {label: "Check", action: () => D().optionsDialog("CHECK")},
        {label: "Analysis", action: () => D().optionsDialog("ANALYSIS")},
        "-",
        {label: "Windows Settings", action: () => D().windowSettings(), enabled: () => !!(S.state && S.state.plots.active)},
        {label: "Colors", action: () => D().colours()},
        {label: "Font", disabled: true, tag: "P2"},
        "-",
        {label: "Shader Options", action: () => D().shaderOptions()},
        {label: "Reset Plot", action: () => S.command(["RSTVIEW", "RSTCENTER"]), enabled: () => P().activeIs3D()}]},
      {title: "View", items: [
        {label: "Check Errors", action: () => D().checkErrors()},
        {label: "Command Window", action: () => S.commandWindow()},
        {label: "Command Display", sub: [
          {label: "Command Echo", check: () => S.display("echo"), action: () => S.toggleDisplay("echo")},
          {label: "Output Confirmation", check: () => S.display("confirm"), action: () => S.toggleDisplay("confirm")},
          {label: "Comments", check: () => S.display("comments"), action: () => S.toggleDisplay("comments")},
          {label: "Warnings & Errors", check: () => S.display("warnerr"), action: () => S.toggleDisplay("warnerr")}]},
        {label: "Toolbars", sub: [
          {label: "Main Toolbar", check: () => S.toolbarsVisible.main, action: () => S.toggleToolbar("main")},
          {label: "Plot Toolbar", check: () => S.toolbarsVisible.plot, action: () => S.toggleToolbar("plot")}]},
        {label: "Split View", check: () => S.split, action: () => S.setSplit(!S.split), tip: "plots and results beside the work instead of a tab in front of it"},
        {label: "Run Summary", action: () => S.openRunSummary && S.openRunSummary(), tip: "key inputs, key outputs and graphs of the active model's run"},
        {label: "Run View", action: () => SASSI.RunView && SASSI.RunView.open(), tip: "an input file as it runs: its listing beside the model, the part being processed highlighted"},
        {label: "Properties", check: () => !!(SASSI.Props && SASSI.Props.isShown()), action: () => SASSI.Props && SASSI.Props.toggle(), tip: "the selected nodes and elements: their properties, edited in the model and in its input file"},
        "-",
        {label: "Results Browser", action: () => S.openResults()}]},
      // Learn: the guided course (static/learn.js); its lessons are read from the server
      {title: "Learn", items: (SASSI.Learn && SASSI.Learn.menuItems()) || []},
      {title: "Help", items: [
        {label: "Help", key: "F1", action: () => S.openHelp()},
        {label: "About", action: () => D().about()},
        {label: "Verification", action: () => S.openVerification(), tag: "VERIFY"}]},
    ];
  }

  function buildMenu(items) {
    const box = el("div", {class: "menu-items"});
    for (const it of items) {
      if (it === "-") { box.appendChild(el("div", {class: "menu-sep"})); continue; }
      // items that need the local server (processes, executables) or only make sense in the browser version
      if ((S.web && it.web === false) || (!S.web && it.webOnly)) continue;
      const cls = ["mi"];
      if (it.sub) cls.push("sub");
      const mi = el("div", {class: cls.join(" "), title: it.tip || null}, el("span", {text: it.label}),
        it.key ? el("span", {class: "key", text: it.key}) : null,
        it.tag ? el("span", {class: "tag", text: it.tag}) : null);
      mi._def = it;
      if (it.sub) {
        const subBox = buildMenu(it.sub);
        mi.appendChild(subBox);
        // a submenu that would leave the window on the right opens on the side with more room, and is
        // narrowed to it (its labels are ellipsized) when neither side is wide enough (narrow windows)
        mi.addEventListener("mouseenter", () => {
          subBox.classList.remove("flip");
          subBox.style.maxWidth = "";
          subBox.style.minWidth = "";
          const r = mi.getBoundingClientRect(), w = subBox.offsetWidth;
          const right = window.innerWidth - r.right - 4, left = r.left - 4;
          if (w <= right) return;
          if (left > right) subBox.classList.add("flip");
          if (w > Math.max(left, right)) { subBox.style.maxWidth = Math.max(left, right) + "px"; subBox.style.minWidth = "0"; }
        });
      }
      mi.addEventListener("click", (ev) => {
        ev.stopPropagation();
        if (it.sub || mi.classList.contains("disabled")) return;
        closeMenus();
        if (it.action) it.action();
      });
      box.appendChild(mi);
    }
    while (box.lastChild && box.lastChild.classList.contains("menu-sep")) box.lastChild.remove();   // after a hidden last item
    return box;
  }
  function refreshMenuStates(root) {
    root.querySelectorAll(".mi").forEach((mi) => {
      const d = mi._def;
      if (!d) return;
      const dis = d.disabled || (d.enabled && !d.enabled());
      mi.classList.toggle("disabled", !!dis);
      if (d.check) mi.classList.toggle("checked", !!d.check());
    });
  }
  function closeMenus() { document.querySelectorAll(".menu.open").forEach((m) => m.classList.remove("open")); }
  S.closeMenus = closeMenus;
  function buildMenubar() {
    const bar = $("#menubar");
    bar.innerHTML = "";
    for (const m of menuDefs()) {
      const menu = el("div", {class: "menu"}, el("div", {class: "menu-title", text: m.title}));
      menu.appendChild(buildMenu(m.items));
      const title = menu.firstChild;
      title.addEventListener("click", (ev) => {
        ev.stopPropagation();
        const open = menu.classList.contains("open");
        closeMenus();
        if (!open) { refreshMenuStates(menu); menu.classList.add("open"); }
      });
      title.addEventListener("mouseenter", () => {
        if (document.querySelector(".menu.open") && !menu.classList.contains("open")) {
          closeMenus();
          refreshMenuStates(menu);
          menu.classList.add("open");
        }
      });
      bar.appendChild(menu);
    }
    document.addEventListener("click", closeMenus);
  }
  /** Rebuild the menu bar (the Learn menu lists the lessons and the progress). */
  S.rebuildMenus = buildMenubar;

  // ------------------------------------------------------------------ View settings
  S.display = (g) => !(S.state && S.state.settings && S.state.settings.display && S.state.settings.display[g] === false);
  S.toggleDisplay = async function (g) {
    const d = Object.assign({}, S.state.settings.display);
    d[g] = !S.display(g);
    try {
      const r = await S.post("/api/settings/display", {display: d});
      S.state.settings = r.settings;
      S.applyDisplayFilters();
    } catch (e) { S.local("ERROR", e.message); }
  };
  S.toggleToolbar = function (which) {
    S.toolbarsVisible[which] = !S.toolbarsVisible[which];
    $("#toolbar-" + which).classList.toggle("hidden", !S.toolbarsVisible[which]);
  };

  // ------------------------------------------------------------------ toolbars (5.9)
  S.updateToolbars = function () { if (SASSI.Toolbars) SASSI.Toolbars.update(); };

  // ------------------------------------------------------------------ module runs (UI-04)
  S.jobs = {};
  S.runModule = async function (module) {
    try {
      const job = await S.post(`/api/run/${module}`);
      S.openJobTab(job);
    } catch (e) {
      // a refusal the server already wrote into the Command History is not repeated
      if (!(e.data && e.data.reported)) S.local("ERROR", `RUN${module}: ${e.message}`);
    }
  };
  /** Module output tab of a job (Modules menu, a RUN<MODULE> typed in Command Entry, Help >
   *  Verification); a job that already has its tab is brought to the front. */
  S.openJobTab = function (job) {
    const id = "job" + job.id;
    if (S.tab(id)) { S.selectTab(id); return S.tab(id); }
    const pre = el("pre", {class: "textview"});
    const info = el("span", {class: "grow"});
    // browser version: a module run executes in the Python engine of the page and cannot be interrupted
    const cancel = el("button", {class: "btn small", text: "Cancel", hidden: S.web, onclick: () => S.post(`/api/jobs/${job.id}/cancel`).catch((e) => S.local("ERROR", e.message))});
    const pane = el("div", {class: "pane"}, el("div", {class: "pane-bar"}, info, cancel), el("div", {class: "pane-body"}, pre));
    // a module run with an analysis word (RUNLOADGEN,STATIC / DYNAMIC) names it in the tab title
    const word = String(job.line || "").split(",")[1] || "";
    const title = job.kind === "verify" ? "Verification output" : `${job.module}${/^[A-Za-z]+$/.test(word.trim()) ? " " + word.trim().toUpperCase() : ""} output`;
    // a new run of a module replaces the finished output tab of its previous run (the listing file keeps
    // the text), so a lesson or a long .pre does not pile up tabs
    for (const o of S.tabs.slice()) {
      if (o.kind === "job" && o.title === title && !o.head.classList.contains("running")) S.removeTab(o.id);
    }
    const t = S.addTab({id, title, kind: "job", pane, jobId: job.id});
    t.head.classList.add("running");
    S.jobs[job.id] = {tab: t, pre, info, cancel, next: 0, job};
    if (!S.keepFocus) S.selectTab(id);          // a lesson step in the full-width layout keeps the lesson in front
    followJob(job.id);
    return t;
  };
  async function followJob(jid) {
    const j = S.jobs[jid];
    while (j) {
      let d;
      try { d = await S.get(`/api/jobs/${jid}?since=${j.next}`); } catch (e) { await S.sleep(1000); continue; }
      applyJob(j, d);
      if (!["starting", "running"].includes(j.job.state)) break;
      await S.sleep(400);
    }
  }
  /** One state of a job (GET /api/jobs/<id>?since=k, or pushed by the browser version's worker): new listing
   *  lines into its output tab, the info line, the end.  The messages are numbered (d.next - d.messages.length
   *  is the number of the first): those already shown -- pushed while this answer was waiting -- are skipped. */
  function applyJob(j, d) {
    const first = d.next - d.messages.length;
    const msgs = d.messages.slice(Math.max(0, j.next - first));
    if (msgs.length) {
      const atEnd = j.pre.parentElement.scrollTop + j.pre.parentElement.clientHeight >= j.pre.parentElement.scrollHeight - 30;
      j.pre.appendChild(document.createTextNode(msgs.map((m) => (m.kind === "ERROR" ? "*** ERROR: " : m.kind === "WARNING" ? "*** WARNING: " : "") + m.text).join("\n") + "\n"));
      if (atEnd) j.pre.parentElement.scrollTop = j.pre.parentElement.scrollHeight;
    }
    j.next = Math.max(j.next, d.next);
    if (j.job && !["starting", "running"].includes(j.job.state) && ["starting", "running"].includes(d.state)) return;   // an older state
    j.job = d;
    j.info.textContent = `${d.line}: ${d.state}${d.progress_text ? " — " + d.progress_text : ""} (${d.elapsed.toFixed(1)} s)` +
      (d.listing ? `  listing ${d.listing}` : "");
    if (j.onUpdate) j.onUpdate(d);
    if (!["starting", "running"].includes(d.state)) {
      j.cancel.disabled = true;
      j.tab.head.classList.remove("running");
    }
  }
  function onJobEvent(job) {
    if (S.activity) { try { S.activity.onJob(job); } catch (e) { console.error(e); } }
    if (SASSI.RunView) { try { SASSI.RunView.onJob(job); } catch (e) { console.error(e); } }
    S.state.job = ["starting", "running"].includes(job.state) ? job : null;
    // a job started by typed command text (RUNSITE in Command Entry, commands after it) gets its tab
    if (S.state.job && !S.jobs[job.id] && job.kind === "module") S.openJobTab(job);
    const cancel = $("#status-cancel");
    if (S.state.job) {
      if (!S.web) cancel.classList.add("on");           // browser version: a running module cannot be cancelled
      cancel.onclick = () => S.post(`/api/jobs/${job.id}/cancel`).catch((e) => S.local("ERROR", e.message));
      S.progress(job.progress || 0, `${job.line}: ${job.progress_text || job.state}`);
    } else {
      cancel.classList.remove("on");
      S.progress(null);
      S.status(`${job.line}: ${job.state}${job.ok === false && job.state !== "cancelled" ? " (see the output tab)" : ""}`);
      if (job.kind === "module") P().refreshResults();
    }
  }

  // ------------------------------------------------------------------ editor tabs (File > Open, 5.6)
  S.editors = {};             // open File Editors by path (the Run view marks the line being executed)
  /** Mark the section and the line of file ``path`` being executed in its File Editor (m = {line, start, end},
   *  0-based; null clears).  The path of the run (resolved by the interpreter) or, failing that, the only
   *  open editor of a file of that name. */
  S.editorHighlight = function (path, m) {
    let e = S.editors[path];
    if (!e && path) {
      const base = path.split(/[\\/]/).pop();
      const same = Object.values(S.editors).filter((x) => x.path.split(/[\\/]/).pop() === base);
      if (same.length === 1) e = same[0];
    }
    if (e) e.highlight(m);
  };
  S.openEditor = async function (path) {
    let d;
    try {
      d = await S.post("/api/file", {name: path});      // creates the file when missing
    } catch (e) { S.local("ERROR", `File > Open ${path}: ${e.message}`); return; }
    const id = "edit:" + d.path;
    if (S.tab(id)) { S.selectTab(id); return; }
    const ta = el("textarea", {class: "editor", spellcheck: "false", wrap: "off"});
    ta.value = d.text || "";
    // the Run view marks the section and the line of this file being executed (behind the transparent text)
    const hl = el("div", {class: "editor-hl", "aria-hidden": "true"}, el("div", {class: "hl-sec", hidden: true}), el("div", {class: "hl-cur", hidden: true}));
    const chg = el("div", {class: "hl-chgs"});      // lines the Properties panel changed (D-W6-16), for a while
    hl.appendChild(chg);
    let mark = null, changed = [], chgTimer = null;
    const metrics = () => { const cs = getComputedStyle(ta); return {lh: parseFloat(cs.lineHeight) || 17.4, pt: parseFloat(cs.paddingTop) || 8}; };
    const placeMark = () => {
      const [sec, cur] = hl.children;
      const {lh, pt} = metrics(), y = (i) => pt + i * lh - ta.scrollTop;
      chg.replaceChildren(...changed.map((i) => el("div", {class: "hl-chg", style: {top: `${y(i)}px`, height: `${lh}px`}})));
      sec.hidden = cur.hidden = !mark;
      if (!mark) return;
      sec.style.top = `${y(mark.start)}px`; sec.style.height = `${(mark.end - mark.start + 1) * lh}px`;
      cur.style.top = `${y(mark.line)}px`; cur.style.height = `${lh}px`;
    };
    ta.addEventListener("scroll", placeMark);
    let dirty = false;
    const status = el("span", {class: "grow", text: d.created ? "new file" : ""});
    const connect = el("button", {class: "btn small", text: "Input ▸ Connect to Command Entry"});
    const save = async () => {
      try {
        await S.post("/api/file", {name: d.path, text: ta.value});
        dirty = false;
        status.textContent = "saved " + new Date().toLocaleTimeString();
        S.setTabTitle(id, `File Editor - ${d.path}`);
      } catch (e) { status.textContent = "save failed: " + e.message; }
    };
    const run = async () => { await save(); S.command(`INP,${S.q(d.path)}`, {skipEditor: id}); };
    const pane = el("div", {class: "pane"},
      el("div", {class: "pane-bar"},
        el("button", {class: "btn small", text: "File ▸ Save", onclick: save}),
        connect,
        el("button", {class: "btn small", text: "Run (INP)", title: "save, then execute the buffer with INP (extension)", onclick: run}),
        status),
      el("div", {class: "editor-wrap"}, hl, ta));
    const t = S.addTab({id, title: `File Editor - ${d.path}`, kind: "editor", pane, onClose: () => { delete S.editors[d.path]; },
      beforeClose: () => !dirty || SASSI.D.confirm("File Editor", `${d.path} has unsaved changes. Close anyway?`, "Close", "Cancel")});
    ta.addEventListener("input", () => { if (!dirty) { dirty = true; S.setTabTitle(id, `File Editor - ${d.path} *`); } });
    ta.addEventListener("keydown", (ev) => { if ((ev.ctrlKey || ev.metaKey) && ev.key === "s") { ev.preventDefault(); save(); } });
    const editor = {tabId: id, path: d.path, append: (line) => {
      if (ta.value && !ta.value.endsWith("\n")) ta.value += "\n";
      ta.value += line + "\n";
      if (!dirty) { dirty = true; S.setTabTitle(id, `File Editor - ${d.path} *`); }
    }, highlight: (m) => {
      mark = m;
      if (m) {
        const {lh, pt} = metrics(), y = pt + m.line * lh;
        if (y < ta.scrollTop + lh || y > ta.scrollTop + ta.clientHeight - 2 * lh) ta.scrollTop = Math.max(0, y - ta.clientHeight / 3);
      }
      placeMark();
    }, text: () => ta.value, isDirty: () => dirty,
    /** Replace the buffer (the Properties panel changed the file): saved = the text is the file's. */
    setText: (text, saved) => {
      const top = ta.scrollTop;
      ta.value = text;
      ta.scrollTop = top;
      dirty = !saved;
      S.setTabTitle(id, `File Editor - ${d.path}${dirty ? " *" : ""}`);
      if (saved) status.textContent = "updated " + new Date().toLocaleTimeString();
    },
    /** Mark lines (0-based) for a few seconds and bring the first one into view. */
    flash: (lines) => {
      changed = lines.slice();
      if (changed.length) {
        const {lh, pt} = metrics(), y0 = pt + changed[0] * lh;
        if (y0 < ta.scrollTop || y0 > ta.scrollTop + ta.clientHeight - 2 * lh) ta.scrollTop = Math.max(0, y0 - ta.clientHeight / 3);
      }
      placeMark();
      clearTimeout(chgTimer);
      chgTimer = setTimeout(() => { changed = []; placeMark(); }, 6000);
    }};
    S.editors[d.path] = editor;
    connect.addEventListener("click", () => {
      if (S.connectedEditor && S.connectedEditor.tabId === id) {
        S.connectedEditor = null;
        connect.classList.remove("primary");
        connect.textContent = "Input ▸ Connect to Command Entry";
      } else {
        document.querySelectorAll(".connected-btn").forEach((b) => { b.classList.remove("primary", "connected-btn"); b.textContent = "Input ▸ Connect to Command Entry"; });
        S.connectedEditor = editor;                  // one window at a time
        connect.classList.add("primary", "connected-btn");
        connect.textContent = "Connected (click to disconnect)";
      }
    });
    S.selectTab(t.id);
  };

  // ------------------------------------------------------------------ results browser, help, verification tabs
  S.openResults = function () { P().openResults(); };
  S.openHelp = function () { D().helpTab(); };
  S.openVerification = function () { D().verificationTab(); };

  // ------------------------------------------------------------------ events
  async function refreshState() {
    S.state = await S.get("/api/state");
    updateTitle();
    S.applyDisplayFilters();
    S.applyColours();
  }
  S.refreshState = refreshState;
  function handleEvent(ev) {
    if (S.activity) { try { S.activity.onEvent(ev); } catch (e) { console.error(e); } }   // the activity panel
    if (SASSI.RunView && (ev.type === "run" || ev.type === "progress")) {                 // the Run view
      try { if (ev.type === "run") SASSI.RunView.onEvent(ev); else SASSI.RunView.onProgress(ev); } catch (e) { console.error(e); }
    }
    switch (ev.type) {
      case "message": S.appendMessage(ev.kind, ev.text); break;
      case "plot": P().onEvent(ev); break;
      case "progress":
        if (ev.line >= ev.total) { S.progress(null); S.status(`${ev.source}: ${ev.total} lines done`); }
        else S.progress(ev.fraction, ev.text);
        break;
      case "state":
        S.state.models = ev.models;
        S.state.active = ev.active;
        S.state.cwd = ev.cwd;
        S.state.revision = ev.revision;
        updateTitle();
        P().onModelChanged();
        if (SASSI.Props) SASSI.Props.refresh();          // the Properties panel (selection and values)
        break;
      case "job": onJobEvent(ev.job); break;
      case "check":         // CHECK / AFWRITE produced messages (requirements 5.5, spec 05a section 3)
        if (!ev.suppress) D().checkErrors({refresh: true});
        break;
      default: break;
    }
  }
  /** Events after the last one handled (those a push already delivered are skipped by number). */
  function handleEvents(events, last) {
    for (const ev of events) {
      if (ev.seq <= S.lastEvent) continue;
      try { handleEvent(ev); } catch (e) { console.error(e); }
      S.lastEvent = ev.seq;
    }
    S.lastEvent = Math.max(S.lastEvent, last);
  }
  async function pollLoop() {
    let failures = 0;
    for (;;) {
      try {
        // browser version: no long poll (it would block the Python engine); a short poll every 250 ms when
        // idle, and at once after a push (webWait)
        const r = await S.get(`/api/events?since=${S.lastEvent}&timeout=${S.web ? 0 : 20}`);
        failures = 0;
        if (r.reset) await refreshState();
        handleEvents(r.events, r.last);
        if (S.web) await webWait(250);
      } catch (e) {
        failures++;
        S.status(`connection to the SASSI-EDU ${S.web ? "Python engine" : "server"} lost (${e.message}); retrying ...`);
        await S.sleep(Math.min(5000, 500 * failures));
      }
    }
  }
  // ------------------------------------------------------------------ pushes (browser version)
  let wake = null;
  /** Wait ms milliseconds, or less when a push arrives. */
  function webWait(ms) {
    return new Promise((resolve) => {
      const t = setTimeout(() => { wake = null; resolve(); }, ms);
      wake = () => { clearTimeout(t); wake = null; resolve(); };
    });
  }
  /** A push of the browser version's worker (web/boot.js), sent while a long command or a module run keeps
   *  the Python engine busy (no request is answered then): {events, last, reset} -- the new events (Command
   *  History, status-bar progress, job states) -- and {job} -- the job with its new listing lines, the same
   *  object as GET /api/jobs/<id>?since=k. */
  S.push = function (p) {
    if (!S.state) return;                         // not started yet: the first poll reads everything
    if (p.events && !p.reset) handleEvents(p.events, p.last || 0);
    if (p.job && S.jobs[p.job.id]) applyJob(S.jobs[p.job.id], p.job);
    if (wake) wake();
  };

  // ------------------------------------------------------------------ Command Entry
  function setupCommandEntry() {
    const inp = $("#cmd");
    inp.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") {
        const line = inp.value;
        if (!line.trim()) return;
        S.recall.push(line);
        S.recallPos = S.recall.length;
        inp.value = "";
        S.command(line).catch(() => {});
      } else if (ev.key === "ArrowUp") {
        if (S.recallPos > 0) { S.recallPos--; inp.value = S.recall[S.recallPos]; }
        ev.preventDefault();
      } else if (ev.key === "ArrowDown") {
        if (S.recallPos < S.recall.length - 1) { S.recallPos++; inp.value = S.recall[S.recallPos]; }
        else { S.recallPos = S.recall.length; inp.value = ""; }
        ev.preventDefault();
      }
    });
  }
  function setupKeys() {
    // the group last clicked in: with the split view the plot keys (Insert, Home, PageUp ...) rotate the plot
    // only after a click in the plots group, so that they still scroll the Command History or a lesson
    document.addEventListener("pointerdown", (ev) => {
      const g = ev.target && ev.target.closest ? ev.target.closest(".tabgroup") : null;
      if (g) S.pointerGroup = g.id === "group-view" ? "view" : "main";
    }, true);
    document.addEventListener("keydown", (ev) => {
      if ((ev.ctrlKey || ev.metaKey) && ev.key.toLowerCase() === "o") { ev.preventDefault(); D().loadModel(); return; }
      if (ev.key === "F1") { ev.preventDefault(); S.openHelp(); return; }
      if (ev.key === "Escape") { closeMenus(); if (D().closeTop) D().closeTop(); return; }
      const tag = (ev.target && ev.target.tagName) || "";
      if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return;
      if (P().onKey) P().onKey(ev);
    });
  }

  // ------------------------------------------------------------------ LaTeX formulas (KaTeX)
  /** Typeset the formulas under root: the elements with data-tex that sassi/ui/markdown.py writes for
   *  $...$, $$...$$ and math blocks (lessons, Help pages).  Without KaTeX, or for a formula KaTeX rejects,
   *  the TeX source stays as text (an error is outlined, the message is its tooltip). */
  S.renderMath = function (root) {
    root = root || document;
    if (!window.katex || !root.querySelectorAll) return;
    for (const e of root.querySelectorAll("[data-tex]:not([data-math])")) {
      e.setAttribute("data-math", "1");
      try {
        window.katex.render(e.getAttribute("data-tex"), e,
          {displayMode: e.classList.contains("math-display"), throwOnError: true, strict: "ignore", trust: false});
      } catch (err) {
        e.classList.add("math-error");
        e.title = "LaTeX: " + String((err && err.message) || err);
      }
    }
  };
  /** Every page part that shows Markdown (lesson, Help, explainer, dialogs) is typeset when it appears. */
  function watchMath() {
    let pending = false;
    const run = () => { pending = false; S.renderMath(document); };
    const obs = new MutationObserver(() => { if (!pending) { pending = true; requestAnimationFrame(run); } });
    for (const id of ["app", "modal-root"]) { const n = document.getElementById(id); if (n) obs.observe(n, {childList: true, subtree: true}); }
    run();
  }

  // ------------------------------------------------------------------ start
  S.start = async function () {
    watchMath();
    buildMenubar();
    setupSplitResizer();
    makeHistoryTab();
    S.selectTab("history");
    setupCommandEntry();
    setupKeys();
    if (SASSI.Toolbars) SASSI.Toolbars.build();
    try {
      await refreshState();
    } catch (e) {
      S.status(`cannot reach the ${S.web ? "Python engine" : "server"}: ${e.message}`);
      return;
    }
    S.recall = (S.state.history || []).slice();
    S.recallPos = S.recall.length;
    // messages emitted before this page was loaded (start-up warnings, earlier commands)
    try {
      const r = await S.get(`/api/events?since=0&timeout=0`);
      for (const ev of r.events) if (ev.type === "message" && ev.seq <= S.state.events) S.appendMessage(ev.kind, ev.text);
    } catch (e) { /* not fatal */ }
    S.lastEvent = S.state.events;
    P().restorePlots(S.state.plots);
    if (S.state.job) { onJobEvent(S.state.job); S.openJobTab(S.state.job); }
    // browser version: how Python started (web/boot.js) -- from this computer, or what was downloaded
    const web = window.SASSI_WEB;
    if (S.web && web && web.readyTime) S.local("LOCAL", `Python ready in ${web.readyTime.toFixed(1)} s. ${web.startNote || ""}`.trim());
    S.updateToolbars();
    S.status("Ready");
    if (S.web) S.local("INFO", "SASSI-EDU in the browser: Python runs in this tab (Pyodide). The workspace " +
      `${S.state.cwd} is kept in memory: a reload starts afresh (course progress is kept). Use File > Upload to ` +
      "Workspace... and File > Download... to copy files in and out.");
    $("#cmd").focus();
    pollLoop();
    // Learn: the start page (opened automatically) and the lesson of this session after a page reload
    if (SASSI.Learn) SASSI.Learn.init({hadPlots: S.tabs.some((t) => t.plotId)}).catch((e) => S.local("WARNING", `Learn: ${e.message}`));
  };
  return S;
})();
