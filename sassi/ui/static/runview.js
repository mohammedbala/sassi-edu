/* SASSI-EDU GUI -- the Run view (requirements 7.20, D-W6-15): an input file as it runs.
 *
 * While an INP file runs, its File Editor -- opened by itself when it is not open yet -- marks the section
 * being processed and the line being executed, and the Run view (right-hand group of the split view) draws
 * the model as the file builds it (sassi/ui/runview.py: snapshots of the active model, with the free-field
 * soil around it), the part of the model the last line added or changed highlighted: nodes, elements, soil
 * layers, interaction nodes, fixities, masses, output requests.  While a module runs it highlights the part
 * the module works on: the soil layers (SITE, SOIL), the interaction nodes (POINT, the impedance and the
 * solution of ANALYS), the elements (HOUSE, the dynamic stiffness of ANALYS), the loaded nodes (FORCE), the
 * output nodes and elements (MOTION, STRESS, RELDISP).  Its bar names the file, the line and the section.
 *
 * Speed: "full speed" draws what the run reaches (the commands that build a model take milliseconds); "watch"
 * and "slow" pause after every line that changes the model (POST /api/runview {pace}).  The editor and the
 * view open by itself when a file starts (unless "open on run" is off, or while a lesson is open).  It shows;
 * it never changes the model (no command text).
 *
 * Events: run {file, end, line, scene} (the server), progress {line, source} and {kind: "step"} and job (the
 * activity panel's: the module steps).
 */
"use strict";

const SASSI_RUNVIEW = (function () {
  const S = SASSI, el = S.el;
  const R = {};
  const P = () => SASSI.P;
  const FOLLOW_KEY = "sassi-edu.runview.follow", PACE_KEY = "sassi-edu.runview.pace";
  const lsGet = (k) => { try { return localStorage.getItem(k); } catch (e) { return null; } };
  const lsSet = (k, v) => { try { localStorage.setItem(k, v); } catch (e) { /* private mode */ } };
  const HI = "#ff8c00", HI_DARK = "#8a4500";
  const PACES = [[0, "full speed"], [0.12, "watch"], [0.4, "slow"]];
  R.follow = () => lsGet(FOLLOW_KEY) !== "0";
  R.pace = () => { const v = Number(lsGet(PACE_KEY)); return PACES.some((p) => p[0] === v) ? v : 0; };

  // ---------------------------------------------------------------- the sections of an input file
  const isComment = (s) => /^\s*[*!]/.test(s);
  const isBreak = (s) => /^\s*[*!]?\s*$/.test(s);               // a blank line or a lone "*"
  const PART = /^\s*[*!]\s*[-=~#]{4,}\s*(.*?)\s*[-=~#]*\s*$/;    // "* ------- structure"
  const clean = (s) => s.replace(/^\s*[*!]+\s*/, "").replace(/^[-=~#\s]+|[-=~#\s]+$/g, "").trim();
  /** The sections of a listing: blocks of lines between blank lines and lone "*" lines.  A block's title is
   *  its first comment text; a "* ---- name" line starts a part, shown before the titles of its blocks.
   *  Returns [{start, end, title, part}] (0-based line indices, end included). */
  R.sections = function (lines) {
    const out = [];
    let part = "", cur = null;
    lines.forEach((ln, i) => {
      if (isBreak(ln)) { cur = null; return; }
      const m = PART.exec(ln);
      if (m) {
        part = clean(m[1] || "") || part;
        cur = null;
        if (!m[1]) return;
      }
      if (!cur) {
        cur = {start: i, end: i, title: "", part};
        out.push(cur);
      }
      cur.end = i;
      if (!cur.title && isComment(ln) && !m) cur.title = clean(ln);
    });
    for (const s of out) if (!s.title) s.title = s.part || `line ${s.start + 1}`;
    return out;
  };
  const sectionOf = (f, i) => f.sections.find((s) => i >= s.start && i <= s.end) || null;

  // ---------------------------------------------------------------- state
  const run = {files: [], last: null, scene: null, focus: null, text: "", step: null, active: false};
  let view = null;                    // the tab: {id, pane, list, plot (a plot-like object for P.drawModel), ...}
  let drawTimer = null;
  const top = () => run.files[run.files.length - 1] || run.last;     // the running file, else the last one run

  // ---------------------------------------------------------------- module steps -> the part of the model
  const STEP_FOCUS = {
    "SITE.modes": () => ({soil: true}),
    "SITE.field": () => ({soil: true}),
    "POINT.green": (r) => ({soil: true, nodes: r.interaction}),
    "HOUSE.elements": () => ({allElements: true}),
    "HOUSE.assembled": (r) => ({allElements: true, nodes: r.interaction}),
    "FORCE.loads": (r) => ({nodes: r.loads}),
    "ANALYS.impedance": (r) => ({nodes: r.interaction}),
    "ANALYS.dynamic": () => ({allElements: true}),
    "ANALYS.factor": (r) => ({allElements: true, nodes: r.interaction}),
    "ANALYS.solve": (r) => ({allElements: true, nodes: r.interaction}),
    "MOTION.interp": (r) => ({nodes: r.nout}),
    "MOTION.convolve": (r) => ({nodes: r.nout}),
    "MOTION.spectra": (r, d) => ({nodes: [Number(d.node)]}),
    "STRESS.elements": (r) => ({elements: r.eout}),
    "RELDISP.relative": (r, d) => ({nodes: [Number(d.node)]}),
    "SOIL.iterate": () => ({soil: true}),
  };
  const WHAT = {soil: "the soil layers", interaction: "the interaction nodes", elements: "the elements"};
  R.stepFocus = function (key, data, lists) {
    const f = STEP_FOCUS[key];
    return f ? f(lists || {}, data || {}) : {};
  };

  // ---------------------------------------------------------------- events
  /** The File Editor of a file (the run's resolved path, else the only open editor of a file of that name). */
  const editorFor = (path) => {
    const E = S.editors || {};
    if (!path) return null;
    if (E[path]) return E[path];
    const b = path.split(/[\\/]/).pop(), same = Object.values(E).filter((x) => x.path.split(/[\\/]/).pop() === b);
    return same.length === 1 ? same[0] : null;
  };
  const lessonOpen = () => !!(SASSI.Learn && SASSI.Learn.isLessonOpen && SASSI.Learn.isLessonOpen());
  R.onEvent = function (ev) {
    if (ev.type !== "run") return;
    if (ev.event === "file") {
      run.files.push({path: ev.path, name: ev.name, lines: ev.lines || [], sections: R.sections(ev.lines || []), cur: -1, depth: ev.depth});
      if (ev.depth === 1) {
        run.active = true; run.last = null; run.focus = null; run.step = null; run.text = "";
        if (R.follow() && !lessonOpen()) {
          // the file being run in its File Editor (opened from the lines the run sent: no request, which the
          // browser version could answer only after the run), the model being built beside it
          if (ev.kind === "file" && ev.path && !editorFor(ev.path) && S.openEditor) S.openEditor(ev.path, {text: (ev.lines || []).join("\n") + "\n"});
          if (!view) R.open({auto: true});
        }
      }
      paintWhere();
    } else if (ev.event === "end") {
      const f = run.files.pop();
      if (f && S.editorHighlight) S.editorHighlight(f.path, null);
      if (!run.files.length) { run.active = false; run.step = null; if (f) { f.done = true; run.last = f; } }
      paintWhere();
      paintHead();
      scheduleDraw();
    } else if (ev.event === "line") {
      const f = run.files.find((x) => x.path === ev.path) || top();
      if (f) f.cur = ev.line - 1;
      run.focus = ev.focus || null;
      run.text = ev.text || "";
      run.step = null;
      paintWhere();
      paintHead();
      scheduleDraw();
    } else if (ev.event === "scene") {
      if (!ev.data || ev.data.stale) return;
      run.scene = ev.data;
      scheduleDraw();
    }
  };
  /** progress events of the activity panel: the line being executed (all lines, every 0.15 s) and the steps
   *  of a module run inside the file. */
  R.onProgress = function (ev) {
    if (ev.kind === "step") { setStep(ev.module, ev.key, ev.data); return; }
    if (ev.kind === "module") return;
    if (ev.line === undefined) return;
    const f = top();
    if (f && (!ev.source || ev.source === f.name)) {
      f.cur = ev.line - 1;
      paintWhere();
    }
  };
  /** A module run in a worker (Modules menu, RUN<MODULE> typed): its steps; at its start the model is drawn
   *  again (nothing runs in the server's interpreter, so the snapshot can be read). */
  R.onJob = function (job) {
    if (!job || job.kind !== "module") return;
    if (job.state === "starting" || job.state === "running") {
      if (!run.active && view && !run.scene) loadScene();
      if (job.step) setStep(job.module, job.step.key, job.step.data);
    } else if (run.step && !run.active) { run.step = null; paintHead(); scheduleDraw(); }
  };
  function setStep(module, key, data) {
    run.step = {module, key, data: data || {}};
    paintHead();
    scheduleDraw();
  }

  // ---------------------------------------------------------------- the tab: the model being built
  /** Snapshots are made while a file runs if the view is open or opens by itself (the browser version answers
   *  a request only after the run: the setting must be on before it starts). */
  const wantScenes = () => R.follow() || !!view;
  const sendSettings = () => S.post("/api/runview", {enabled: wantScenes(), pace: R.pace()}).catch(() => {});
  R.open = function (opts) {
    opts = opts || {};
    if (view && S.tab(view.id)) { if (!opts.auto) S.selectTab(view.id); return; }
    const plotHost = el("div", {class: "rv-plot"});
    const title = el("span", {class: "rv-file"});
    const where = el("span", {class: "rv-where"});
    const now = el("div", {class: "rv-now"});
    const speed = el("select", {class: "rv-speed", title: "pause after every line that changes the model, to watch it being built"},
      ...PACES.map(([v, label]) => el("option", {value: String(v), text: label})));
    speed.value = String(R.pace());
    const follow = el("input", {type: "checkbox"});
    follow.checked = R.follow();
    const pane = el("div", {class: "pane runview"},
      el("div", {class: "pane-bar rv-bar"}, title, where, el("span", {class: "grow"}),
        el("label", {class: "rv-opt", title: "speed of the next lines of a running file"}, "Speed ", speed),
        el("label", {class: "rv-opt", title: "open the input file and the Run view when an input file starts"}, follow, " open on run")),
      now,
      el("div", {class: "rv-body"}, plotHost));
    speed.addEventListener("change", () => { lsSet(PACE_KEY, speed.value); sendSettings(); });
    follow.addEventListener("change", () => { lsSet(FOLLOW_KEY, follow.checked ? "1" : "0"); sendSettings(); });
    view = {id: "runview", title: "Run view", kind: "runview", pane, title_: title, where, now, plot: {host: plotHost}};
    const t = S.addTab({id: view.id, title: "Run view", kind: "runview", pane,
      onClose: () => { view = null; document.body.classList.remove("runview-shown"); sendSettings(); },
      onActivate: () => scheduleDraw()});
    view.tab = t;
    if (window.ResizeObserver) {
      const ro = new ResizeObserver(() => { if (view && view.plot.plotDiv) P().resizeHost(view.plot); });
      ro.observe(pane);
    }
    sendSettings();
    if (!opts.auto || !S.keepFocus || t.group === "view") S.selectTab(t.id);
    if (!run.active) loadScene();
    paintWhere();
    paintHead();
    scheduleDraw();
  };
  R.isOpen = () => !!view;
  // the activity panel moves to the left while the Run view is on the screen (styles.css: runview-shown)
  const shown = () => document.body.classList.toggle("runview-shown", !!view && !!view.tab && S.isShown(view.tab));
  const selectTab = S.selectTab;
  S.selectTab = function (id, opts) { selectTab(id, opts); shown(); };
  async function loadScene() {
    try { run.scene = await S.get("/api/runview"); } catch (e) { /* the server is busy: the next snapshot comes */ }
    scheduleDraw();
  }

  // ---------------------------------------------------------------- the line being run: in the File Editor
  /** Mark the section and the line being executed in the File Editor of the running file, and name them in
   *  the Run view's bar. */
  function paintWhere() {
    const f = top();
    const sec = f && f.cur >= 0 ? sectionOf(f, f.cur) : null;
    if (f && !f.done && S.editorHighlight) S.editorHighlight(f.path, f.cur >= 0 ? {line: f.cur, start: sec ? sec.start : f.cur, end: sec ? sec.end : f.cur} : null);
    if (!view) return;
    view.title_.textContent = f ? f.name + (run.files.length > 1 ? ` (in ${run.files[run.files.length - 2].name})` : "") : "";
    view.where.textContent = !f ? "" : f.done ? `finished · ${f.lines.length} lines`
      : f.cur >= 0 ? `line ${f.cur + 1} of ${f.lines.length}` + (sec ? ` · ${sec.part && sec.part !== sec.title ? sec.part + " › " : ""}${sec.title}` : "") : "";
  }

  // ---------------------------------------------------------------- the caption
  function paintHead() {
    if (!view) return;
    let text = "";
    if (run.step) {
      const st = S.activity && S.activity.STEPS && S.activity.STEPS[run.step.key];
      const what = stepWhat(R.stepFocus(run.step.key, run.step.data, run.scene && run.scene.run));
      text = `${run.step.module} · ${st ? st.title : run.step.key}${what ? ` — ${what}` : ""}`;
    } else if (run.text) {
      text = run.text;
    } else if (!run.active) {
      text = run.scene && !run.scene.empty ? `Model ${run.scene.model}${run.scene.name ? " (" + run.scene.name + ")" : ""}` : "";
    }
    view.now.textContent = text;
    view.now.classList.toggle("on", !!(run.step || run.text));
  }
  function stepWhat(f) {
    const out = [];
    if (f.soil) out.push(WHAT.soil);
    if (f.allElements) out.push(WHAT.elements);
    if (f.elements && f.elements.length) out.push(`${f.elements.length} output element${f.elements.length > 1 ? "s" : ""}`);
    if (f.nodes && f.nodes.length) {
      const r = run.scene && run.scene.run;
      out.push(r && f.nodes === r.interaction ? WHAT.interaction : f.nodes.length === 1 ? `node ${f.nodes[0]}` : `${f.nodes.length} nodes`);
    }
    return out.join(", ");
  }

  // ---------------------------------------------------------------- the model and its highlight
  function scheduleDraw() {
    if (!view || drawTimer) return;
    drawTimer = setTimeout(() => { drawTimer = null; draw(); }, 60);
  }
  function draw() {
    shown();
    if (!view || !S.isShown(view.tab)) return;
    const d = run.scene, t = view.plot;
    if (!d || d.empty) {
      if (t.kind !== "empty") { try { Plotly.purge(t.plotDiv); } catch (e) { /* none */ } t.host.innerHTML = ""; t.plotDiv = null; t.kind = "empty"; }
      t.host.replaceChildren(el("div", {class: "empty", text: run.active ? "Building the model: it appears with its first nodes or soil layers."
        : d && d.model !== null && d.model !== undefined ? "The model has no nodes or soil layers yet." : "Run an input file (File Editor > Run, Model > Input, an example's Run all): the model is drawn here as the file builds it."}));
      return;
    }
    if (t.kind !== d.family) { try { Plotly.purge(t.plotDiv); } catch (e) { /* none */ } t.host.innerHTML = ""; t.plotDiv = null; t.kind = d.family; }
    const focus = currentFocus(d);
    if (d.family === "layer") {
      P().drawLayers(t, d);
      layerHighlight(t, d, focus);
      return;
    }
    const any = focus.allElements || focus.elements.length || focus.nodes.length;
    P().drawModel(t, d, focusTraces(d, focus), {note: false, dim: !!any});   // the soil's hover text says what it is
  }
  /** What to highlight now: the module step's part, else the last line's. */
  function currentFocus(d) {
    if (run.step) {
      const f = R.stepFocus(run.step.key, run.step.data, d.run);
      return {nodes: f.nodes || [], elements: f.elements || [], allElements: !!f.allElements, soil: !!f.soil, layers: null};
    }
    const f = run.focus || {};
    if (f.all) return {nodes: [], elements: [], allElements: false, soil: false, layers: null};
    const req = f.requests || {};
    return {nodes: [].concat(f.nodes || [], f.interaction || [], f.fixed || [], f.mass || [], f.loads || [], req.nout || [], req.rdnd || []),
      elements: [].concat(f.elements || [], req.eout || []), allElements: false, soil: !!(f.layers && f.layers.length), layers: f.layers || null};
  }
  function indexes(d) {
    if (d._ix) return d._ix;
    const sc = d.scene, node = new Map(), elem = new Map();
    sc.node_id.forEach((id, i) => node.set(id, i));
    sc.elem_id.forEach((id, k) => elem.set(sc.elem_group[k] + ":" + id, k));
    d._ix = {node, elem};
    return d._ix;
  }
  function focusTraces(d, f) {
    const sc = d.scene, H = P().sceneHelpers, ix = indexes(d), out = [];
    // nodes no element uses yet (a mesh before its elements, orientation nodes): small grey dots
    if (sc.free && sc.free.length) out.push(H.markerTrace(sc, sc.free, {marker: {size: 3, color: "#7a7a7a"}, hover: " (no element yet)"}));
    const dim = f.allElements || f.elements.length || f.nodes.length;       // drawModel then fades the model
    if (dim && sc.interaction && sc.interaction.length) out.push(H.markerTrace(sc, sc.interaction, {marker: {size: 2.5, color: "#e00000"}, hover: " (interaction)"}));
    if (f.soil && d.soil) out.push(...soilHighlight(d, f.layers));
    const ks = new Set();
    if (f.allElements) sc.elem_id.forEach((_, k) => ks.add(k));
    for (const [g, e] of f.elements) { const k = ix.elem.get(g + ":" + e); if (k !== undefined) ks.add(k); }
    if (ks.size) {
      const faces = H.faceList(sc, {all: !f.allElements, only: (k) => ks.has(k)});
      if (faces.length) {
        out.push(Object.assign(H.meshTrace(sc, faces, {color: () => HI, opacity: f.allElements ? 0.35 : 0.6}), {hoverinfo: "skip"}));
        out.push(H.outlineTrace(sc, faces, {color: HI, width: f.allElements ? 2 : 4}));
      }
      const x = [], y = [], z = [];
      sc.edges.forEach((ed, k) => {
        if (!ks.has(sc.edge_elem[k])) return;
        for (const i of ed) { x.push(sc.xyz[i][0]); y.push(sc.xyz[i][1]); z.push(sc.xyz[i][2]); }
        x.push(null); y.push(null); z.push(null);
      });
      if (x.length) out.push({type: "scatter3d", mode: "lines", x, y, z, line: {color: HI, width: 10}, opacity: 0.85, hoverinfo: "skip"});
    }
    const ni = [];
    for (const id of new Set(f.nodes)) { const i = ix.node.get(Number(id)); if (i !== undefined) ni.push(i); }
    if (ni.length) out.push(H.markerTrace(sc, ni, {marker: {size: ni.length > 200 ? 4 : 6, color: HI, line: {color: HI_DARK, width: 1}}, hover: " (now)"}));
    return out;
  }
  /** Outline the soil bands of the highlighted layers (L numbers; all of them for a module step). */
  function soilHighlight(d, layers) {
    const so = d.soil, want = layers ? new Set(layers.map(Number)) : null;
    const x = [], y = [], z = [];
    so.quads.forEach((q, i) => {
      const b = so.bands[so.quad_band[i]];
      if (!b || (want && !want.has(Number(b.layer)))) return;
      for (const p of q.concat([q[0]])) { x.push(p[0]); y.push(p[1]); z.push(p[2]); }
      x.push(null); y.push(null); z.push(null);
    });
    return x.length ? [{type: "scatter3d", mode: "lines", x, y, z, line: {color: HI, width: 3}, hoverinfo: "skip"}] : [];
  }
  /** The soil column (no node yet): outline the highlighted layers. */
  function layerHighlight(t, d, f) {
    if (!f.soil || !t.plotDiv || !d.table) return;
    const rows = d.table.layers, want = f.layers ? new Set(f.layers.map(Number)) : null;
    const top0 = rows.length ? rows[0].top : 0;
    const shapes = ((t.plotDiv.layout && t.plotDiv.layout.shapes) || []).filter((s) => !(s.line && s.line.color === HI));
    rows.forEach((r) => {
      if (want && !want.has(Number(r.layer))) return;
      const y0 = r.top - top0;
      shapes.push({type: "rect", x0: 0, x1: 1, y0, y1: y0 + r.thick, fillcolor: "rgba(255,140,0,0.25)", line: {color: HI, width: 3}});
    });
    try { Plotly.relayout(t.plotDiv, {shapes}); } catch (e) { /* hidden */ }
  }

  setTimeout(sendSettings, 800);             // on before the first run (see wantScenes)
  return R;
})();
SASSI.RunView = SASSI_RUNVIEW;
