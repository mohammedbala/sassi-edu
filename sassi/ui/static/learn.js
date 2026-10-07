/* SASSI-EDU GUI -- Learn: the guided course for structural engineers who are new to SASSI.
 *
 *  * the Learn tab (start page): welcome, "How SASSI works in 2 minutes" (the module chain, clickable),
 *    the course outline by part with progress, the examples gallery;
 *  * the lesson player, docked to the right of the tab area so that the plots, listings and the Command
 *    History stay visible next to the lesson: introduction, steps with their commands (explain any line),
 *    Run step / Run steps 1-k, the labelled sections, action buttons, Previous / Next, Reset;
 *  * the command explainer (GET /api/explain): a popover from lesson lines, Command History lines and the
 *    "?" of Command Entry, which also gives a live hint of the argument under the cursor;
 *  * Model > Open Example... and the gallery buttons (copy an example into a fresh workspace, CD there,
 *    open it in the File Editor or run it with INP); the picture of a card builds the model there without
 *    module runs and opens the 3D model view.
 *
 * Rule L17: every step, setup and action submits command text through S.command (POST /api/command), so the
 * Command History of a lesson replays as a .pre file.  Progress (steps completed) is kept per browser in
 * localStorage (wrapped in try/catch: a private window simply forgets it).  Vanilla JavaScript. */
"use strict";

(function (S) {
  const el = S.el;
  const L = {};
  S.Learn = L;
  const D = () => S.D;
  const enc = encodeURIComponent;

  // ================================================================== small helpers
  const store = {
    get(key, dflt) {
      try { const v = window.localStorage.getItem("sassi-edu.learn." + key); return v === null ? dflt : JSON.parse(v); } catch (e) { return dflt; }
    },
    set(key, value) {
      try { window.localStorage.setItem("sassi-edu.learn." + key, JSON.stringify(value)); } catch (e) { /* storage unavailable */ }
    },
  };
  /** HTML rendered by the server (sassi.ui.markdown: escaped, no raw HTML passes) into a container. */
  function html(cls, markup) { const d = el("div", {class: cls}); d.innerHTML = markup || ""; return d; }
  function btn(text, onclick, attrs) { return el("button", Object.assign({class: "btn small", text, onclick}, attrs || {})); }
  /** Links inside lesson / gallery HTML: documents open in the Help tab. */
  function wireDocLinks(root) {
    root.addEventListener("click", (ev) => {
      const a = ev.target.closest ? ev.target.closest("a") : null;
      if (!a || !root.contains(a)) return;
      const doc = a.getAttribute("data-doc");
      if (doc) { ev.preventDefault(); D().helpTab({doc, anchor: a.getAttribute("data-anchor") || undefined}); return; }
      if ((a.getAttribute("href") || "").startsWith("#")) ev.preventDefault();
    });
  }
  const PART_ORDER = ["Fundamentals", "Design applications", "Advanced"];

  // ================================================================== course data
  L.outline = null;
  L.examples = null;
  L.load = async function (force) {
    if (L.outline && !force) return L.outline;
    try { L.outline = await S.get("/api/lessons"); }
    catch (e) { L.outline = {parts: [], count: 0, errors: [], root: ""}; S.local("WARNING", `Learn: ${e.message}`); }
    if (S.rebuildMenus) S.rebuildMenus();
    return L.outline;
  };
  L.loadExamples = async function (force) {
    if (L.examples && !force) return L.examples;
    try { L.examples = (await S.get("/api/examples")).examples; }
    catch (e) { L.examples = []; S.local("WARNING", `Learn: ${e.message}`); }
    return L.examples;
  };
  L.lessons = () => (L.outline ? L.outline.parts.flatMap((p) => p.lessons) : []);
  L.lessonNumber = (id) => L.lessons().findIndex((x) => x.id === id) + 1;
  /** The explainer video of a lesson (static/videos/NN.html, NN = its number in the course) or the list of
   *  them; it opens in a browser tab of its own (narrated by the browser's text-to-speech). */
  L.videoUrl = (id) => (id ? `static/videos/${String(L.lessonNumber(id)).padStart(2, "0")}.html` : "static/videos/index.html");
  L.watchVideo = function (id) {
    const w = window.open(L.videoUrl(id), "_blank");
    if (w) w.opener = null;
    else S.local("WARNING", `The browser blocked the new tab: open ${L.videoUrl(id)} yourself.`);
  };
  L.findLesson = (id) => L.lessons().find((x) => x.id === id) || null;

  // ------------------------------------------------------------------ progress (per browser)
  function progressAll() { const p = store.get("progress", {}); return p && typeof p === "object" ? p : {}; }
  L.progress = function (id) {
    const p = progressAll()[id] || {};
    return {done: Array.isArray(p.done) ? p.done : [], last: Number(p.last) || 0};
  };
  function saveProgress(id, fn) {
    const all = progressAll();
    const p = L.progress(id);
    fn(p);
    all[id] = {done: Array.from(new Set(p.done)).sort((a, b) => a - b), last: p.last, t: Date.now()};
    all._last = id;
    store.set("progress", all);
  }
  L.markDone = (id, k) => saveProgress(id, (p) => { if (!p.done.includes(k)) p.done.push(k); });
  L.setLast = (id, k) => saveProgress(id, (p) => { p.last = k; });
  L.resetProgress = (id) => { const all = progressAll(); if (id) delete all[id]; else for (const k of Object.keys(all)) delete all[k]; store.set("progress", all); };
  L.lastLesson = () => { const id = progressAll()._last; return id && L.findLesson(id) ? id : null; };

  // ================================================================== command runner (L17)
  /** Wait until a worker job has ended; resolves with its final state. */
  L.waitJob = async function (jid) {
    for (;;) {
      let d = null;
      try { d = await S.get(`/api/jobs/${jid}?since=1000000000`); } catch (e) { /* retry */ }
      if (d && !["starting", "running"].includes(d.state)) return d;
      if (L.stopRequested && d && ["starting", "running"].includes(d.state) && !L._cancelSent) {
        L._cancelSent = true;
        S.post(`/api/jobs/${jid}/cancel`).catch(() => {});
      }
      await S.sleep(500);
    }
  };
  /** Run command lines one after the other exactly as if typed in Command Entry: a RUN<MODULE> starts its
   *  worker job (listing streamed into its output tab and the Command History) and the next line waits for
   *  its end.  hooks.onLine(i, state, info) with state running | ok | job | error | skipped.  Stops at the
   *  first error or when L.stop() is called.  Resolves true when every line succeeded. */
  L.runLines = async function (lines, hooks) {
    hooks = hooks || {};
    const on = hooks.onLine || (() => {});
    L.stopRequested = false;
    L._cancelSent = false;
    for (let i = 0; i < lines.length; i++) {
      if (L.stopRequested) { for (let j = i; j < lines.length; j++) on(j, "skipped"); return false; }
      const line = lines[i];
      if (!String(line).trim()) { on(i, "ok"); continue; }
      on(i, "running");
      let res;
      try { res = await S.command(line); }
      catch (e) { on(i, "error", e.message); return false; }
      const r = (res.results || [])[0];
      if (res.job) {
        on(i, "job", res.job);
        const end = await L.waitJob(res.job.id);
        if (!end.ok) { on(i, "error", `${end.line}: ${end.state}${end.state === "cancelled" ? "" : " (see its output tab)"}`); return false; }
      } else if (!r || !r.ok) {
        const errs = (res.messages || []).filter((m) => m.kind === "ERROR").map((m) => m.text);
        on(i, "error", errs.join("; ") || "the command failed");
        return false;
      }
      on(i, "ok");
    }
    return true;
  };
  L.stop = function () { L.stopRequested = true; };

  /** POST a prepare request; when the server answers needs_confirm, ask the learner first. */
  async function prepareWithConfirm(path, body, what) {
    let r = await S.post(path, body || {});
    if (r.needs_confirm) {
      const list = (r.unsaved || []).map((u) => `  Model ${u.number}: ${u.name || "(no name)"}`).join("\n");
      const yes = await D().confirm(`${what}: start a fresh model?`,
        `${what} starts a fresh model in its own workspace folder.\n\nThese models in memory were changed since their last ` +
        `Model > Save:\n${list}\n\nThey will be removed from memory (their files on disk are kept). Continue?`, "Continue", "Cancel");
      if (!yes) return null;
      r = await S.post(path, Object.assign({}, body || {}, {confirm: true}));
    }
    return r;
  }

  // ================================================================== command explainer
  let pop = null, popAnchor = null;
  L.closeExplain = function () { if (pop) { pop.remove(); pop = null; popAnchor = null; } };
  function argValueCell(a) {
    if (a.extra) return el("td", {class: "l xp-extra", text: a.value === null ? "" : a.value});
    if (!a.given) return el("td", {class: "l xp-blank", text: a.value === "" ? "(blank)" : "—"});
    const td = el("td", {class: "l"}, el("code", {text: a.value}));
    if (a.value_meaning) td.appendChild(el("span", {class: "xp-code", text: ` = ${a.value_meaning}`}));
    return td;
  }
  /** The explanation as DOM nodes (opts.compact: no argument table). */
  L.renderExplain = function (r, opts) {
    opts = opts || {};
    const out = [];
    if (r.kind === "blank") { out.push(el("p", {class: "xp-msg", text: r.message})); return out; }
    if (r.kind === "comment") { out.push(el("p", {text: r.summary})); return out; }
    if (r.kind === "unknown") {
      out.push(el("p", {class: "xp-msg", text: r.message}));
      if (r.suggestions && r.suggestions.length) {
        out.push(el("div", {class: "xp-sugg"}, "Did you mean: ", ...r.suggestions.map((s) => el("button", {class: "linkish", text: s, onclick: () => L.explain(s, popAnchor)}))));
      }
      return out;
    }
    const head = el("div", {class: "xp-head"}, el("code", {class: "xp-name", text: r.name}));
    if (r.abbreviated) head.appendChild(el("span", {class: "xp-abbr", text: `typed as ${r.typed} (documented abbreviation)`}));
    out.push(head);
    out.push(el("p", {class: "xp-summary", text: r.summary || r.meaning || ""}));
    if (r.summary && r.meaning && r.summary !== r.meaning && !opts.compact) out.push(el("p", {class: "xp-meaning", text: r.meaning}));
    out.push(el("div", {class: "xp-syntax"}, el("span", {class: "k", text: "Syntax"}), el("code", {text: r.syntax})));
    const chips = el("div", {class: "xp-chips"});
    if (r.configures) chips.appendChild(el("span", {class: "chip xp-conf", text: r.configures}));
    if (r.tier) chips.appendChild(el("span", {class: "chip", title: r.tier_text, text: `tier ${r.tier}`}));
    if (r.status) chips.appendChild(el("span", {class: `chip xp-st-${String(r.status).replace(/\s+/g, "-")}`, title: r.status_text, text: r.status === "yes" ? "implemented" : r.status}));
    if (r.category_title) chips.appendChild(el("span", {class: "chip", title: r.category, text: r.category_title}));
    out.push(chips);
    if (opts.compact) return out;
    if (r.args && r.args.length) {
      const tb = el("tbody");
      for (const a of r.args) {
        const name = a.name ? el("code", {text: a.name}) : el("span", {class: "xp-blank", text: "—"});
        const meaning = el("td", {class: "l"}, a.meaning || (a.label ? a.label : el("span", {class: "xp-blank", text: "see the Command Reference"})));
        if (a.default !== undefined && a.default !== "" && !a.extra) meaning.appendChild(el("span", {class: "xp-def", text: ` (default ${a.default})`}));
        if (a.values && !a.value_meaning && Object.keys(a.values).length && Object.keys(a.values).length <= 8) {
          meaning.appendChild(el("div", {class: "xp-vals", text: Object.entries(a.values).map(([k, v]) => `${k} = ${v}`).join(";  ")}));
        }
        const tr = el("tr", {class: (a.extra ? "xp-extra-row " : "") + (r.cursor_arg === a.pos ? "xp-cur" : "")},
          el("td", {text: a.pos}), el("td", {class: "l"}, name), argValueCell(a), meaning);
        tb.appendChild(tr);
      }
      out.push(el("div", {class: "xp-tablewrap"}, el("table", {class: "grid xp-args"},
        el("thead", {}, el("tr", {}, el("th", {text: "#"}), el("th", {class: "l", text: "Argument"}), el("th", {class: "l", text: "Value given"}), el("th", {class: "l", text: "Meaning"}))), tb)));
    } else {
      out.push(el("p", {class: "xp-blank", text: "This command takes no arguments."}));
    }
    if (r.notes) out.push(el("p", {class: "xp-notes", text: r.notes}));
    const foot = el("div", {class: "xp-foot"});
    if (r.dialog) foot.appendChild(btn(`Open ${r.dialog.label}`, () => { L.closeExplain(); D().optionsDialog(r.dialog.name, r.dialog.tab || undefined); }));
    foot.appendChild(btn("Command Reference", () => { L.closeExplain(); D().helpTab({doc: r.help.doc, anchor: r.help.anchor || undefined}); }));
    if (r.source) foot.appendChild(el("span", {class: "xp-src", text: `source: ${r.source}`}));
    out.push(foot);
    return out;
  };
  /** Show the explainer popover for a command line next to anchor (an element, or null: centred). */
  L.explain = async function (line, anchor) {
    let r;
    try { r = await S.get(`/api/explain?line=${enc(line)}`); } catch (e) { S.local("ERROR", `Explain: ${e.message}`); return; }
    L.closeExplain();
    const close = el("button", {class: "x", title: "Close (Esc)", text: "×", onclick: () => L.closeExplain()});
    pop = el("div", {class: "explain-pop", role: "dialog", "aria-label": "Command explainer"},
      el("div", {class: "xp-bar"}, el("span", {class: "xp-title", text: "Explain"}), el("code", {class: "xp-line", text: line, title: line}), close),
      el("div", {class: "xp-body"}, ...L.renderExplain(r)));
    popAnchor = anchor || null;
    document.body.appendChild(pop);
    placePopover(pop, anchor);
    close.focus();
  };
  function placePopover(p, anchor) {
    const vw = window.innerWidth, vh = window.innerHeight, m = 8;
    const w = Math.min(600, vw - 2 * m);
    p.style.width = w + "px";
    const ph = Math.min(p.offsetHeight, vh - 2 * m);
    let x = (vw - w) / 2, y = Math.max(m, (vh - ph) / 3);
    if (anchor && anchor.getBoundingClientRect) {
      const a = anchor.getBoundingClientRect();
      x = Math.min(Math.max(m, a.left), vw - w - m);
      y = a.bottom + 6;
      if (y + ph > vh - m) y = a.top - ph - 6 >= m ? a.top - ph - 6 : Math.max(m, vh - ph - m);
    }
    p.style.left = x + "px";
    p.style.top = y + "px";
    p.style.maxHeight = (vh - 2 * m) + "px";
  }
  document.addEventListener("mousedown", (ev) => {
    if (pop && !pop.contains(ev.target) && !(popAnchor && popAnchor.contains && popAnchor.contains(ev.target))) L.closeExplain();
  });
  document.addEventListener("keydown", (ev) => { if (ev.key === "Escape" && pop) { L.closeExplain(); ev.stopPropagation(); } }, true);
  window.addEventListener("resize", () => { if (pop) placePopover(pop, popAnchor); });

  // ------------------------------------------------------------------ Command Entry: "?" and the live hint
  let hintOn = store.get("hint", true) !== false;
  let hintSeq = 0, hintTimer = null;
  function setupHint() {
    const inp = S.$("#cmd"), help = S.$("#cmd-help"), box = S.$("#cmd-hint");
    if (!inp || !help || !box) return;
    const paint = () => { help.classList.toggle("on", hintOn); help.setAttribute("aria-pressed", hintOn ? "true" : "false"); };
    paint();
    help.addEventListener("click", () => {
      if (inp.value.trim() && hintOn) { L.explain(inp.value, help); return; }       // a second click: full explanation
      hintOn = !hintOn;
      store.set("hint", hintOn);
      paint();
      update();
      inp.focus();
    });
    const update = () => {
      clearTimeout(hintTimer);
      hintTimer = setTimeout(async () => {
        const line = inp.value;
        if (!hintOn || !line.trim()) { box.hidden = true; return; }
        const seq = ++hintSeq;
        let r;
        try { r = await S.get(`/api/explain?line=${enc(line)}&cursor=${inp.selectionStart || line.length}`); } catch (e) { return; }
        if (seq !== hintSeq) return;
        renderHint(box, r, line);
      }, 120);
    };
    for (const ev of ["input", "keyup", "click", "focus"]) inp.addEventListener(ev, update);
    inp.addEventListener("keydown", (ev) => { if (ev.key === "Enter") setTimeout(() => { box.hidden = true; }, 0); });
  }
  function renderHint(box, r, line) {
    box.innerHTML = "";
    box.hidden = false;
    if (r.kind !== "command") {
      box.appendChild(el("span", {class: "hint-msg", text: r.kind === "unknown" ? `${r.typed}: not a command` + (r.suggestions && r.suggestions.length ? ` (did you mean ${r.suggestions.slice(0, 3).join(", ")}?)` : "") : (r.summary || r.message || "")}));
      return;
    }
    const cur = r.cursor_arg || 0;
    const sig = el("span", {class: "hint-sig"}, el("code", {class: cur === 0 ? "hint-arg on" : "hint-name", text: r.name}));
    const shown = r.args.filter((a) => !a.repeat || a.pos === cur);
    for (const a of shown.slice(0, 24)) {
      sig.appendChild(document.createTextNode(","));
      sig.appendChild(el("code", {class: "hint-arg" + (a.pos === cur ? " on" : "") + (a.extra ? " extra" : ""), text: a.name || `arg${a.pos}`}));
    }
    const a = r.args.find((x) => x.pos === cur);
    const desc = cur === 0 ? (r.summary || r.meaning) :
      a ? `${a.name ? a.name + " — " : ""}${a.meaning || a.label || "(see the Command Reference)"}${a.default ? ` (default ${a.default})` : ""}${a.value_meaning ? `  →  ${a.value} = ${a.value_meaning}` : ""}`
        : "no further argument";
    box.append(sig, el("span", {class: "hint-desc", text: desc}),
      el("button", {class: "linkish", text: "details", title: "full explanation (argument table)", onclick: () => L.explain(line, S.$("#cmd-help"))}));
  }

  // ------------------------------------------------------------------ Command History: click a command line
  function setupHistoryClicks() {
    document.addEventListener("click", (ev) => {
      const d = ev.target.closest ? ev.target.closest(".history .m-ECHO") : null;
      if (!d) return;
      const sel = window.getSelection ? String(window.getSelection()) : "";
      if (sel.trim()) return;                                  // selecting text to copy
      const line = d.dataset.line || d.textContent.replace(/^> /, "");
      L.explain(line, d);
    });
  }

  // ================================================================== How SASSI works (module chain)
  const THEORY = "docs/theory/THEORY_MANUAL.md";
  const MODULES = {
    EQUAKE: {stage: "Seismic input", sub: "design motion", lesson: "06-seismic-input", anchor: "21-spectrum-compatible-motions-equake",
      text: "Generates an acceleration history whose response spectrum matches a target design spectrum (for example an RG 1.60 spectrum anchored to 0.3 g), by frequency-domain matching with a wavelet refinement, and reports the acceptance checks of the manual and of SRP 3.7.1 (spectrum envelope, strong-motion duration, PSD). The history becomes the control motion of the SSI analysis.",
      io: "in: target spectrum (RSIN) - out: acceleration history (ACCOUT), its spectrum (RSOUT)"},
    SOIL: {stage: "Seismic input", sub: "site response", lesson: "06-seismic-input", anchor: "12-the-shake-equivalent-linear-method-soil",
      text: "One-dimensional equivalent-linear site response (the SHAKE method): vertically propagating shear waves in the layered column. The shear modulus and damping of every sublayer are iterated on the strain-dependent curves (DYNP: G/Gmax and damping against shear strain) until they agree with the effective strain. The strain-compatible properties (FILE88) feed SITE with SITEX,1, so the SSI free field is consistent with the earthquake level.",
      io: "in: layers, DYNP curves, rock motion - out: strain-compatible properties (FILE88), motions and spectra at the layers"},
    SITE: {stage: "Free field", sub: "layered site", lesson: "02-free-field", anchor: "4-the-thin-layer-method-site-mode-1",
      text: "Solves the horizontally layered site without the structure (thin-layer method): the layers are cut into sublayers thin compared with the shortest wavelength (h ≤ Vs/(5 f_max)), and the half-space is a stack of extra layers on a viscous boundary (variable-depth method). For the wave field you choose (vertically incident SV, SH or P) it gives the free-field motion at every sublayer interface per unit control motion at the control point.",
      io: "out: free-field motion (FILE1), layer eigen-solutions (FILE2)"},
    POINT: {stage: "Soil impedance", sub: "springs + dashpots", lesson: "03-impedance", anchor: "7-point-load-solutions-central-zone-and-transmitting-boundary-point",
      text: "Computes how the free field deforms under unit harmonic point loads at the depth of each interaction-node level (an axisymmetric model with a central zone of radius R0 and a transmitting boundary). ANALYS assembles from it the flexibility of all interaction nodes and inverts it into the soil impedance X_ff(ω): the frequency-dependent springs and dashpots of the foundation, radiation damping included.",
      io: "in: FILE2 of SITE - out: point-load solutions (FILE3)"},
    HOUSE: {stage: "Structure", sub: "your FE model", lesson: "04-surface-ssi", anchor: "13-element-formulations-house",
      text: "Builds the stiffness and mass matrices of the structure (beams, shells, solids, springs, lumped masses) and of the excavated soil from the finite-element model, with complex-modulus (hysteretic) damping G* = G(1 - 2β² + 2iβ√(1-β²)). This is the part you know from ANSYS: the fixed-base model, but without supports at the interaction nodes. HOUSE does not need POINT; the order is only the run order.",
      io: "in: nodes, elements, materials, masses - out: structure and excavated-soil matrices (FILE4)"},
    ANALYS: {stage: "SSI solution", sub: "per frequency", lesson: "04-surface-ssi", anchor: "9-solution-of-the-ssi-equation-analys",
      text: "For each SSI frequency it solves [(K_s - ω²M_s) - (K_e - ω²M_e) + X_ff] U = X_ff U'_f: structure, minus the excavated soil, plus the soil impedance at the interaction nodes, loaded by the free-field motion U'_f. The result is the transfer function of every degree of freedom per unit control motion. In ANSYS terms: a harmonic analysis with a frequency-dependent complex soil impedance and a seismic load from the free field.",
      io: "in: FILE1, FILE3, FILE4 - out: transfer functions (FILE8)"},
    MOTION: {stage: "Results", sub: "histories + ISRS", lesson: "04-surface-ssi", anchor: "11-convolution-response-spectra-relative-displacements-and-stresses",
      text: "Interpolates the transfer functions between the solved frequencies, multiplies them by the Fourier transform of the control motion and transforms back: acceleration time histories and in-structure response spectra (ISRS) at the output nodes (NOUT), with the damping ratios you request (DAMP).",
      io: "in: FILE8, control motion (THFILE) - out: .TFU/.TFI transfer functions, .ACC histories, .RS spectra"},
    RELDISP: {stage: "Results", sub: "relative displacement", lesson: "08-design-outputs", anchor: "114-relative-displacements-reldisp",
      text: "Relative displacement histories of a node with respect to a reference node or to the free field, from the difference of their transfer functions (no double integration drift): seismic gaps, piping and equipment anchors, building separation.",
      io: "in: .TFI of the nodes (RDND, RELFILE) - out: relative displacement histories"},
    STRESS: {stage: "Results", sub: "member forces", lesson: "08-design-outputs", anchor: "115-stresses-and-forces-stress",
      text: "Element forces, moments and stresses (beams, shells, solids, springs) from the transfer functions and the control motion: maxima and histories at the requested elements (EOUT), the input to member design.",
      io: "in: FILE8, element recovery data of HOUSE - out: element force and stress results"},
  };
  const FLOW = [["EQUAKE", "SOIL"], ["SITE"], ["POINT"], ["HOUSE"], ["ANALYS"], ["MOTION", "RELDISP", "STRESS"]];
  function howPanel() {
    const detail = el("div", {class: "how-detail", "aria-live": "polite"});
    const boxes = [];
    const show = (m) => {
      const d = MODULES[m];
      boxes.forEach((b) => b.classList.toggle("sel", b.dataset.m === m));
      const les = d.lesson && L.findLesson(d.lesson);
      detail.innerHTML = "";
      detail.append(el("h4", {}, el("code", {text: m}), ` — ${d.stage}: ${d.sub}`), el("p", {text: d.text}), el("p", {class: "how-io", text: d.io}),
        el("div", {class: "how-btns"},
          les ? btn(`Lesson ${L.lessonNumber(les.id)}: ${les.title}`, () => L.openLesson(les.id), {class: "btn small primary"}) : null,
          btn("Theory Manual", () => D().helpTab({doc: THEORY, anchor: d.anchor})),
          btn(`Explain RUN${m}`, (ev) => L.explain(`RUN${m}`, ev.currentTarget))));
    };
    const flow = el("div", {class: "how-flow", role: "group", "aria-label": "SASSI module chain"});
    FLOW.forEach((col, i) => {
      if (i) flow.appendChild(el("div", {class: "how-arrow", "aria-hidden": "true", text: "→"}));
      const c = el("div", {class: "how-col"}, el("div", {class: "how-stage", text: MODULES[col[0]].stage}));
      for (const m of col) {
        const b = el("button", {class: "how-box", "data-m": m, title: `${m}: ${MODULES[m].sub} (click for details)`}, el("strong", {text: m}), el("small", {text: MODULES[m].sub}));
        b.addEventListener("click", () => show(m));
        boxes.push(b);
        c.appendChild(b);
      }
      flow.appendChild(c);
    });
    show("ANALYS");
    return el("section", {class: "learn-sec", id: "learn-how"},
      el("h2", {text: "How SASSI works in 2 minutes"}),
      el("p", {text: "A fixed-base analysis shakes the base of your model with the design motion. SASSI instead splits the soil-structure system into the free field (the layered site without the structure), the soil impedance at the interaction nodes, and the structure, and solves their coupling in the frequency domain, one frequency at a time. Each step is a module; the model commands write one input deck per module (AFWRITE) and RUN<MODULE> runs it. Click a module:"}),
      flow, detail,
      el("p", {class: "how-eq"}, "The substructuring equation solved by ANALYS: ", el("code", {text: "[(K_s - ω²M_s) - (K_e - ω²M_e) + X_ff] U = X_ff U'_f"}), " ",
        el("button", {class: "linkish", text: "Theory Manual §3", onclick: () => D().helpTab({doc: THEORY, anchor: "3-flexible-volume-substructuring"})})));
  }

  // ================================================================== the Learn tab (start page)
  let startBody = null;
  L.openStart = async function (opts) {
    opts = opts || {};
    if (!S.tab("learn")) {
      startBody = el("div", {class: "pane-body learn-body"});
      wireDocLinks(startBody);
      const go = (id) => () => { const t = startBody.querySelector("#" + id); if (t) startBody.scrollTop = t.offsetTop - 8; };
      const pane = el("div", {class: "pane"},
        el("div", {class: "pane-bar"}, el("strong", {text: "Learn"}), el("span", {class: "grow learn-bar-note", text: "guided course, examples and the command explainer"}),
          btn("How SASSI works", go("learn-how")), btn("Course", go("learn-course")), btn("Examples", go("learn-examples")),
          btn("Explain a command…", () => L.askExplain())),
        startBody);
      S.addTab({id: "learn", title: "Learn", kind: "learn", pane});
    }
    if (opts.select !== false) S.selectTab("learn");
    await Promise.all([L.load(opts.force), L.loadExamples(opts.force)]);
    renderStart();
    if (opts.section) {
      const t = startBody.querySelector("#learn-" + opts.section);
      if (t) startBody.scrollTop = t.offsetTop - 8;
    }
  };
  function renderStart() {
    if (!startBody) return;
    const keep = startBody.scrollTop;
    startBody.innerHTML = "";
    const page = el("div", {class: "learn-page"});
    page.append(heroPanel(), howPanel(), coursePanel(), examplesPanel(), tipsPanel());
    startBody.appendChild(page);
    startBody.scrollTop = keep;
  }
  L.refreshStart = renderStart;
  function heroPanel() {
    const last = L.lastLesson();
    const all = L.lessons();
    const firstOpen = all.find((x) => L.progress(x.id).done.length < x.nsteps) || all[0];
    const cont = last ? L.findLesson(last) : null;
    const cta = el("div", {class: "learn-cta"});
    if (cont) {
      const p = L.progress(cont.id);
      cta.appendChild(btn(`Continue: lesson ${L.lessonNumber(cont.id)}, ${p.last ? "step " + p.last : "introduction"}`, () => L.openLesson(cont.id, {step: p.last}), {class: "btn primary"}));
    }
    if (firstOpen && (!cont || firstOpen.id !== cont.id)) cta.appendChild(btn(cont ? `Next lesson: ${firstOpen.title}` : "Start the course", () => L.openLesson(firstOpen.id), {class: cont ? "btn" : "btn primary"}));
    cta.appendChild(btn("Open an example", () => L.examplesDialog(), {class: "btn"}));
    const showAtStart = el("input", {type: "checkbox", id: "learn-autostart"});
    showAtStart.checked = store.get("autostart", true) !== false;
    showAtStart.addEventListener("change", () => store.set("autostart", showAtStart.checked));
    const root = (L.outline && L.outline.root) || (S.state && S.state.course_root) || "";
    return el("section", {class: "learn-hero"},
      el("h1", {text: "Learn soil-structure interaction with SASSI-EDU"}),
      el("p", {text: "You design safety-related structures with fixed-base finite-element models: a design spectrum or a time history at the base, modes, in-structure response spectra, member forces. Soil-structure interaction changes all of these. The soil is flexible and radiates energy away, the foundation does not follow the free-field motion, and the structure's frequencies and damping shift."}),
      el("p", {text: "This course takes you through SASSI one command at a time. Each step shows the commands and runs them in this window, so the plots and listings appear next to the lesson. It then explains what the commands do, why they matter for your design deliverables (ISRS, member forces, relative displacements), their technical basis, and how they map to what you would do in ANSYS. Click any command to see its arguments explained."}),
      cta,
      el("div", {class: "learn-meta"},
        el("span", {}, "Workspaces: ", el("code", {text: root || "(not set)"})), btn("Change…", () => L.workspaceDialog()),
        btn("Free disk space…", () => L.freeSpaceDialog(), {title: "show the size of the lesson and example workspaces and delete those you no longer need"}),
        el("label", {class: "learn-auto"}, showAtStart, " Show this page when the GUI starts")));
  }
  function coursePanel() {
    const sec = el("section", {class: "learn-sec", id: "learn-course"}, el("h2", {text: "The course"}),
      el("p", {}, "Every lesson has a short narrated explainer video (3-4 minutes) with its big idea in plain words: ",
        btn("▶ Explainer videos", () => L.watchVideo(null), {class: "btn small linkbtn"})));
    const o = L.outline || {parts: []};
    if (!o.parts.length) sec.appendChild(el("p", {class: "empty", text: "No lessons are installed yet (sassi/ui/lessons/*.md). The examples below can be loaded and run already."}));
    let n = 0;
    for (const part of o.parts) {
      const grid = el("div", {class: "lcards"});
      for (const les of part.lessons) grid.appendChild(lessonCard(les, ++n));
      sec.append(el("h3", {class: "learn-part", text: part.part}), grid);
    }
    for (const e of (o.errors || [])) sec.appendChild(el("p", {class: "learn-err", text: `${e.file}: ${e.error}`}));
    return sec;
  }
  function lessonCard(les, n) {
    const p = L.progress(les.id);
    const done = p.done.filter((k) => k >= 1 && k <= les.nsteps).length;
    const frac = les.nsteps ? done / les.nsteps : 0;
    const state = !p.done.length && !p.last ? "Start" : done >= les.nsteps ? "Review" : "Continue";
    const pre = (les.prerequisites || []).map((id) => { const x = L.findLesson(id); return x ? `${L.lessonNumber(id)}` : id; });
    const card = el("article", {class: "lcard" + (done >= les.nsteps && les.nsteps ? " complete" : "")},
      el("div", {class: "lcard-head"}, el("span", {class: "lnum", text: n}), el("h4", {text: les.title}), el("span", {class: "lmin", text: les.minutes ? `${les.minutes} min` : ""})),
      el("p", {class: "lsum", text: les.summary || ""}),
      les.objectives && les.objectives.length ? el("ul", {class: "lobj"}, ...les.objectives.map((o) => el("li", {text: o}))) : null,
      el("div", {class: "lmeta"}, les.example ? `Example ${les.example}` : "", pre.length ? `${les.example ? " · " : ""}after lesson ${pre.join(", ")}` : ""),
      el("div", {class: "lprog", title: `${done} of ${les.nsteps} steps completed`}, el("div", {class: "bar", style: {width: `${Math.round(100 * frac)}%`}})),
      el("div", {class: "lcard-foot"}, el("span", {class: "lsteps", text: `${done} / ${les.nsteps} steps`}),
        p.done.length || p.last ? btn("Reset progress", async () => {
          if (await D().confirm("Reset progress", `Forget the completed steps of "${les.title}" in this browser?`, "Reset", "Cancel")) { L.resetProgress(les.id); renderStart(); if (S.rebuildMenus) S.rebuildMenus(); }
        }, {class: "btn small linkbtn"}) : null,
        btn("▶ Explainer", () => L.watchVideo(les.id), {class: "btn small linkbtn", title: "the narrated explainer video of this lesson (new browser tab)"}),
        btn(state, () => L.openLesson(les.id, {step: state === "Continue" ? p.last : 0}), {class: "btn small primary"})));
    return card;
  }
  function examplesPanel() {
    const sec = el("section", {class: "learn-sec", id: "learn-examples"}, el("h2", {text: "Examples"}),
      el("p", {text: "Nine complete models, each commented line by line. Click a picture to build the model in a fresh workspace and turn it around in the 3D view. Load one into a fresh workspace to read it in the File Editor (and run it with Run (INP)), or run it all at once and explore the results (File > Results Browser, Plot > Spectrum). When a lesson is built on the example, the guided lesson is the best way in."}));
    const grid = el("div", {class: "xcards"});
    for (const x of (L.examples || [])) grid.appendChild(exampleCard(x));
    if (!(L.examples || []).length) grid.appendChild(el("p", {class: "empty", text: "No examples found."}));
    sec.appendChild(grid);
    return sec;
  }
  /** The preview picture of an example card (python -m sassi.ui.thumbnails); a click shows the model in the
   *  3D view (the soil layers for an example without a structure).  A missing picture removes the band. */
  function exampleThumb(x) {
    if (!x.thumb) return null;
    const layers = x.view === "LAYERPLOT";
    const label = x.number ? `example ${x.number}` : x.name;
    return el("button", {class: "xthumb", type: "button", "data-hint": layers ? "Open the layer plot" : "Open in 3D view",
      title: `Show ${layers ? "the soil layers" : "the model"} of ${label}: copy it into a fresh workspace folder, build the model (no module runs) and open ${x.view || "MODELPLOT"}`,
      onclick: () => { closeExamplesDialog(); L.showExampleModel(x.name); }},
    el("img", {loading: "lazy", decoding: "async", width: 480, height: 300, src: x.thumb,
      alt: `${layers ? "Soil layer column" : "3D model"} of ${label}: ${x.title}`,
      onerror: (ev) => { const b = ev.target.closest(".xthumb"); if (b) b.remove(); }}));
  }
  function exampleCard(x, compact) {
    const lessonsHere = (x.lessons || []).filter((l) => L.findLesson(l.id));
    return el("article", {class: "xcard"},
      exampleThumb(x),
      el("div", {class: "xhead"}, el("span", {class: "xnum", text: x.number ? `Example ${x.number}` : x.name}), el("h4", {text: x.title})),
      el("div", {class: "xmeta"}, el("code", {text: x.file}), x.runtime ? ` · run time ${x.runtime}` : "", x.results ? ` · results in ${x.results}/` : ""),
      x.topic ? el("p", {class: "xtopic", text: x.topic}) : null,
      x.modules ? el("p", {class: "xmods"}, el("span", {class: "k", text: "Modules: "}), x.modules) : null,
      !compact && x.learn && x.learn.length ? el("div", {class: "xlearn"}, el("div", {class: "k", text: "What you learn"}), el("ul", {}, ...x.learn.map((t) => el("li", {text: t})))) : null,
      !compact && x.physics ? el("details", {class: "xphys"}, el("summary", {text: "The physics / the model"}), el("p", {text: x.physics})) : null,
      el("div", {class: "xbtns"},
        ...lessonsHere.map((l) => btn(lessonsHere.length > 1 ? `Lesson ${L.lessonNumber(l.id)}` : "Open guided lesson", () => { closeExamplesDialog(); L.openLesson(l.id); }, {class: "btn small primary", title: `Guided lesson ${L.lessonNumber(l.id)}: ${l.title}`})),
        btn("Load into workspace", () => { closeExamplesDialog(); L.loadExample(x.name, false); }, {title: "copy it into a fresh workspace folder, CD there and open the .pre in the File Editor"}),
        btn("Run all", () => { closeExamplesDialog(); L.loadExample(x.name, true); }, {title: `copy it into a fresh workspace folder, CD there and run it with INP${x.runtime ? " (" + x.runtime + ")" : ""}`})));
  }
  function tipsPanel() {
    return el("section", {class: "learn-sec learn-tips", id: "learn-tips"}, el("h2", {text: "While you work"}),
      el("ul", {},
        el("li", {}, el("strong", {text: "Explain any command. "}), "Click a command line in a lesson or in the Command History. While you type in Command Entry, the bar above it names the argument under the cursor (the ", el("code", {text: "?"}), " button switches it on and off; click it again for the full table)."),
        el("li", {}, el("strong", {text: "Everything is a command. "}), "Lessons, menus and dialogs all send command text to one interpreter, so the Command History is a replayable .pre file of what you did."),
        el("li", {}, el("strong", {text: "Your own variations. "}), "The ✎ button of a lesson line copies it into Command Entry; edit it and press Enter. Reset lesson starts again from a fresh workspace."),
        el("li", {}, el("strong", {text: "Reference. "}), "Help (F1) has the User Guide, the Theory Manual, the Verification Manual and the Command Reference.")));
  }

  // ------------------------------------------------------------------ examples
  let exDlg = null;
  function closeExamplesDialog() { if (exDlg) { exDlg.close(); exDlg = null; } }
  /** Model > Open Example...: the gallery as a dialog. */
  L.examplesDialog = async function () {
    await Promise.all([L.load(), L.loadExamples()]);
    const list = el("div", {class: "xcards xcards-dlg"}, ...(L.examples || []).map((x) => exampleCard(x, true)));
    exDlg = D().modal({title: "Open Example", width: "860px", body: el("div", {class: "dlg xdlg"},
      el("p", {class: "opt-note", text: "Each example is copied into a fresh folder under the course workspace folder; the original in examples/ is not changed."}), list),
      buttons: [{label: "Close", action: () => true}], onClose: () => { exDlg = null; }});
  };
  /** Copy example name into a fresh workspace, CD there, fresh model; then open it in the File Editor or run it. */
  L.loadExample = async function (name, run) {
    let r;
    try { r = await prepareWithConfirm(`/api/examples/${enc(name)}/prepare`, {run: !!run}, `Example ${name}`); }
    catch (e) { S.local("ERROR", `Example ${name}: ${e.message}`); D().alert("Example", e.message); return; }
    if (!r) return;
    S.local("LOCAL", `Example ${name}: workspace ${r.workspace}`);
    if (run) S.selectTab("history");
    const ok = await L.runLines(r.commands);
    if (!ok) { S.local("WARNING", `Example ${name}: stopped at an error (see the Command History)`); return; }
    if (!run) await S.openEditor(r.pre);
    else S.status(`Example ${name} done -- File > Results Browser shows its results`);
  };
  /** The picture of an example card: copy it into a fresh workspace, CD there, fresh model, then its model
   *  commands (without CHECK, AFWRITE, RUN<MODULE> ...) and MODELPLOT -- LAYERPLOT for an example without a
   *  structure -- all as command text (L17), so the Command History rebuilds the model. */
  L.showExampleModel = async function (name) {
    let r;
    try { r = await prepareWithConfirm(`/api/examples/${enc(name)}/prepare`, {model: true}, `Example ${name}`); }
    catch (e) { S.local("ERROR", `Example ${name}: ${e.message}`); D().alert("Example", e.message); return; }
    if (!r) return;
    S.local("LOCAL", `Example ${name}: workspace ${r.workspace}; building the model (no module runs), then ${r.view}`);
    const ok = await L.runLines(r.commands);
    if (!ok) { S.local("WARNING", `Example ${name}: stopped at an error (see the Command History)`); return; }
    S.status(`Example ${name}: ${r.view === "LAYERPLOT" ? "soil layers" : "model"} shown -- Load into workspace or Run all to analyse it`);
  };
  L.askExplain = async function () {
    const line = await D().ask("Explain a command", "Command line (e.g. NGEN,3,5,1,5,1,0,0,10)", "");
    if (line && line.trim()) L.explain(line.trim(), null);
  };
  const fmtBytes = (b) => b >= 1e9 ? `${(b / 1e9).toFixed(2)} GB` : b >= 1e6 ? `${(b / 1e6).toFixed(1)} MB` : `${Math.max(0, Math.round(b / 1e3))} kB`;
  /** Learn > Free Disk Space...: the course workspaces on disk and their size; delete those not in use
   *  (only folders the course created are ever deleted; a lesson recreates its workspace when opened). */
  L.freeSpaceDialog = async function () {
    let d;
    try { d = await S.get("/api/course/workspaces"); }
    catch (e) { D().alert("Free Disk Space", e.message); return; }
    const rows = d.workspaces.map((w) => el("tr", {},
      el("td", {class: "l"}, el("code", {text: w.name})), el("td", {text: fmtBytes(w.bytes)}),
      el("td", {class: "l", text: w.in_use ? "in use (kept)" : ""})));
    const free = d.workspaces.filter((w) => !w.in_use);
    const freeBytes = free.reduce((a, w) => a + w.bytes, 0);
    const body = el("div", {class: "dlg"},
      el("p", {}, "Course workspaces in ", el("code", {text: d.root}), `: ${d.workspaces.length}, ${fmtBytes(d.bytes)} in total. Lessons and examples recreate their workspace when you open them again, so deleting them loses only the results of runs you made there.`),
      d.workspaces.length ? el("table", {class: "grid"}, el("thead", {}, el("tr", {}, el("th", {class: "l", text: "Workspace"}), el("th", {text: "Size"}), el("th", {text: ""}))), el("tbody", {}, ...rows))
        : el("p", {class: "empty", text: "There are no course workspaces on disk."}));
    D().modal({title: "Free Disk Space", width: "560px", body, buttons: [
      ...(free.length ? [{label: `Delete ${free.length} workspace${free.length > 1 ? "s" : ""} (${fmtBytes(freeBytes)})`, primary: true, action: async () => {
        try {
          const r = await S.post("/api/course/workspaces/delete", {names: free.map((w) => w.name)});
          S.local("LOCAL", `Free Disk Space: deleted ${r.deleted.length} course workspace(s), ${fmtBytes(r.freed)} freed` +
            (r.skipped.length ? `; kept ${r.skipped.map((s) => `${s.name} (${s.reason})`).join(", ")}` : ""));
          await L.load(true); renderStart();
        } catch (e) { D().alert("Free Disk Space", e.message); }
        return true;
      }}] : []),
      {label: "Close", action: () => true}]});
  };
  L.workspaceDialog = async function () {
    const cur = (L.outline && L.outline.root) || (S.state && S.state.course_root) || "";
    const v = await D().ask("Course Workspace Folder", "Folder of the lesson and example workspaces (blank = <start folder>/sassi-course)", (S.state && S.state.settings && S.state.settings.course && S.state.settings.course.root) || cur);
    if (v === null) return;
    try {
      const r = await S.post("/api/settings/course", {course: {root: v}});
      S.state.settings = r.settings;
      await S.refreshState();
      await L.load(true);
      S.local("LOCAL", `Course workspace folder: ${S.state.course_root}`);
      renderStart();
    } catch (e) { D().alert("Course Workspace Folder", e.message); }
  };

  // ================================================================== lesson player
  const P = {id: null, lesson: null, page: 0, ran: new Set(), lineState: {}, running: false, status: "", statusKind: ""};
  L.player = P;
  const dock = () => S.$("#lessondock");
  /** Layout of the lesson player: "full" -- a full-width tab of its own (the default; the results open in
   *  their tabs and "Back to lesson" returns) -- or "side" -- a panel docked beside the results. */
  L.layout = () => (store.get("layout", "full") === "side" ? "side" : "full");
  const isFull = () => L.layout() === "full";
  const LESSON_TAB = "lesson";
  function moveDockHome() {
    const d = dock(), main = S.$("main");
    if (!d) return;
    d.classList.remove("full");
    if (main && d.parentElement !== main) main.appendChild(d);
  }
  function ensureLessonTab() {
    let t = S.tab(LESSON_TAB);
    if (!t) {
      const pane = el("div", {class: "pane lesson-pane"});
      t = S.addTab({id: LESSON_TAB, title: "Lesson", kind: "lesson", pane,
        beforeClose: () => { L.closeLesson(); return false; },          // closing the tab closes the lesson
        onClose: () => moveDockHome()});
    }
    const d = dock();
    if (d.parentElement !== t.pane) t.pane.appendChild(d);
    d.classList.add("full");
    return t;
  }
  function lessonTabTitle() {
    if (!P.lesson || !S.tab(LESSON_TAB)) return;
    const n = P.lesson.steps.length;
    S.setTabTitle(LESSON_TAB, `Lesson ${L.lessonNumber(P.id) || ""} · ${P.page ? `step ${P.page}/${n}` : "introduction"}`);
  }
  function dockVisible(on) {
    const d = dock(), rz = S.$("#dock-resizer");
    if (!d) return;
    if (on && isFull()) {
      const t = ensureLessonTab();
      d.hidden = false;
      if (rz) rz.hidden = true;
      document.body.classList.remove("dock-open");
      lessonTabTitle();
      S.selectTab(t.id);
    } else {
      if (S.tab(LESSON_TAB)) S.removeTab(LESSON_TAB);                // its onClose moves the panel home
      moveDockHome();
      d.hidden = !on;
      if (rz) rz.hidden = !on;
      document.body.classList.toggle("dock-open", !!on);
    }
    updateBackChip();
  }
  /** Full-width layout: a "Back to lesson" button over the tab area while another tab is in front. */
  let backChip = null;
  function updateBackChip() {
    const d = dock();
    // the lesson tab is in the work group: the chip shows while another tab of that group hides it (with the
    // split view a plot opens beside the lesson and does not)
    const front = S.groups ? S.groups.main.active : S.activeTab;
    const show = !!(isFull() && P.lesson && d && !d.hidden && front && front.id !== LESSON_TAB);
    if (!backChip) {
      backChip = el("button", {class: "lesson-back", type: "button", hidden: true, title: "back to the lesson tab",
        onclick: () => S.selectTab(LESSON_TAB)});
      const wa = S.$("#group-main") || S.$("#workarea");
      if (wa) wa.appendChild(backChip);
    }
    backChip.hidden = !show;
    if (show) backChip.textContent = `◀ Back to lesson ${L.lessonNumber(P.id) || ""} · ${P.page ? "step " + P.page : "introduction"}`;
  }
  const selectTab0 = S.selectTab;
  S.selectTab = function (id, opts) { selectTab0(id, opts); updateBackChip(); };
  /** Switch the layout (Learn > Lesson Layout, the button in the lesson header). */
  L.setLayout = function (mode) {
    store.set("layout", mode === "side" ? "side" : "full");
    const d = dock();
    if (P.lesson && d && !d.hidden) { showDock(); renderPlayer(); }
    if (S.rebuildMenus) S.rebuildMenus();
  };
  /** While a step runs in the full-width layout the lesson stays in front: module output tabs and plots
   *  made by the step's commands open behind it (S.keepFocus); the step's buttons bring them forward. */
  async function keepingFocus(fn) {
    const before = S.keepFocus;
    S.keepFocus = isFull();
    try { return await fn(); } finally { S.keepFocus = before; }
  }
  function applyDockWidth() {
    const d = dock();
    if (!d) return;
    if (isFull()) { d.style.width = ""; d.classList.remove("collapsed"); return; }
    const w = Number(store.get("dockWidth", 0)) || 0;
    if (w) d.style.width = Math.max(300, Math.min(w, Math.round(window.innerWidth * 0.62))) + "px";
    else d.style.width = "";
    d.classList.toggle("collapsed", !!store.get("dockCollapsed", false));
  }
  function setupResizer() {
    const rz = S.$("#dock-resizer");
    if (!rz) return;
    rz.addEventListener("mousedown", (ev) => {
      ev.preventDefault();
      const move = (e) => {
        const w = Math.max(300, Math.min(window.innerWidth - e.clientX, Math.round(window.innerWidth * 0.62)));
        dock().style.width = w + "px";
      };
      const up = () => {
        document.removeEventListener("mousemove", move);
        document.removeEventListener("mouseup", up);
        document.body.classList.remove("resizing");
        store.set("dockWidth", dock().getBoundingClientRect().width);
      };
      document.body.classList.add("resizing");
      document.addEventListener("mousemove", move);
      document.addEventListener("mouseup", up);
    });
    rz.addEventListener("dblclick", () => { store.set("dockWidth", 0); applyDockWidth(); });
  }
  L.closeLesson = function () { dockVisible(false); if (S.rebuildMenus) S.rebuildMenus(); };
  /** Is a lesson open (its tab or its side panel)?  The run summary window does not open by itself then. */
  L.isLessonOpen = () => !!(P.lesson && dock() && !dock().hidden);
  L.toggleCollapse = function () { store.set("dockCollapsed", !store.get("dockCollapsed", false)); applyDockWidth(); };

  /** Open lesson id in the dock.  opts.step: the page to show (0 = introduction); opts.fresh: always start
   *  from a fresh workspace (Reset).  A lesson that is the open lesson of this server session (e.g. after a
   *  page reload) is shown as it is; any other lesson gets a fresh workspace, a fresh model and its setup. */
  L.openLesson = async function (id, opts) {
    opts = opts || {};
    if (P.running) { D().alert("Lesson", "A lesson step is still running; wait for it to finish or press Stop."); return; }
    L.closeExplain();
    let d;
    try { d = await S.get(`/api/lessons/${enc(id)}`); } catch (e) { D().alert("Lesson", e.message); return; }
    if (!L.outline) await L.load();
    const resume = d.session && !opts.fresh;
    if (!resume) {
      let r;
      try { r = await prepareWithConfirm(`/api/lessons/${enc(id)}/open`, {}, `Lesson "${d.title}"`); }
      catch (e) { D().alert("Lesson", e.message); return; }
      if (!r) return;
      await closeLessonTabs();
      Object.assign(P, {id, lesson: d, ran: new Set(), lineState: {}, page: 0, baseTabs: new Set(S.tabs.map((t) => t.id))});
      d.session = {id, workspace: r.workspace, ran: []};
      showDock();
      setStatus("Preparing the workspace: CD, fresh model" + (r.setup.length ? ", setup commands" : "") + " …", "run");
      P.running = true;
      renderPlayer();
      const first = r.commands.length - r.setup.length;          // CD and the fresh model come first
      P.lineState.setup = [];
      const ok = await keepingFocus(() => L.runLines(r.commands, {onLine: (i, st, info) => { if (i >= first) { P.lineState.setup[i - first] = {st, info}; paintLines(); } }}));
      P.running = false;
      setStatus(ok ? `Workspace ready: ${r.workspace}` : "The setup stopped at an error: see the Command History, then Reset the lesson.", ok ? "ok" : "err");
    } else {
      Object.assign(P, {id, lesson: d, ran: new Set(d.session.ran || []), lineState: {}, page: 0, baseTabs: P.baseTabs || new Set(S.tabs.map((t) => t.id))});
      setStatus(`Workspace: ${d.session.workspace}`, "");
      showDock();
    }
    const want = opts.step !== undefined && opts.step !== null ? Number(opts.step) : (resume ? L.progress(id).last : 0);
    goPage(Math.max(0, Math.min(want || 0, d.steps.length)));
    if (S.rebuildMenus) S.rebuildMenus();
  };
  /** Opening a lesson (or Reset) closes what the previous lesson opened, so the tab bar does not fill up:
   *  its plots (ACTIVATEPLOT + CLOSEPLOT, rule L17), its finished module output tabs and its File Editor
   *  tabs of course-workspace files (a file with unsaved edits asks first).  Tabs that were open before
   *  the lesson are kept. */
  async function closeLessonTabs() {
    const base = P.baseTabs;
    if (!base) return;
    const root = (L.outline && L.outline.root) || (S.state && S.state.course_root) || "";
    const plots = [];
    for (const t of S.tabs.slice()) {
      if (base.has(t.id)) continue;
      if (t.plotId) plots.push(t.plotId);
      else if (t.kind === "job" && !t.head.classList.contains("running")) S.removeTab(t.id);
      else if (t.kind === "editor" && root && t.id.startsWith("edit:" + root)) S.requestCloseTab(t.id);
    }
    if (plots.length) {
      try { await S.command(plots.flatMap((id) => [`ACTIVATEPLOT,${id}`, "CLOSEPLOT"])); } catch (e) { /* reported in the Command History */ }
    }
  }
  function showDock() {
    const d = dock();
    if (!d) return;
    store.set("dockCollapsed", false);
    applyDockWidth();
    dockVisible(true);
  }
  function setStatus(text, kind) {
    P.status = text;
    P.statusKind = kind || "";
    const s = dock() && dock().querySelector(".dock-status");
    if (s) { s.textContent = text; s.title = text; s.className = "dock-status " + (kind || ""); }
  }
  function goPage(k) {
    P.page = k;
    if (P.id) L.setLast(P.id, k);
    const st = k > 0 ? P.lesson.steps[k - 1] : null;
    if (st && !st.commands.length) L.markDone(P.id, k);         // a reading step is done once read
    renderPlayer();
    lessonTabTitle();
    updateBackChip();
    const body = dock().querySelector(".dock-body");
    if (body) body.scrollTop = 0;
    if (startBody && S.isShown(S.tab("learn"))) renderStart();
  }
  function renderPlayer() {
    const d = dock();
    if (!d || !P.lesson) return;
    const les = P.lesson, n = les.steps.length, k = P.page;
    const prog = L.progress(P.id);
    d.innerHTML = "";
    // ---- header
    const pills = el("div", {class: "dock-steps", role: "tablist", "aria-label": "Steps"});
    const pill = (i, label, title) => {
      const b = el("button", {class: "pill" + (i === k ? " cur" : "") + (i > 0 && prog.done.includes(i) ? " done" : "") + (P.ran.has(i) ? " ran" : ""),
        title, text: label, role: "tab", "aria-selected": i === k ? "true" : "false"});
      b.addEventListener("click", () => goPage(i));
      return b;
    };
    pills.appendChild(pill(0, "i", "Introduction"));
    les.steps.forEach((s, i) => pills.appendChild(pill(i + 1, String(i + 1), `${i + 1}. ${s.title}${prog.done.includes(i + 1) ? " (done)" : ""}${P.ran.has(i + 1) ? " - has run in this workspace" : ""}`)));
    const done = prog.done.filter((x) => x >= 1 && x <= n).length;
    const head = el("div", {class: "dock-head"},
      el("div", {class: "dock-top"},
        el("span", {class: "dock-kicker", text: `${les.part || "Lesson"} · lesson ${L.lessonNumber(P.id) || ""}${les.minutes ? ` · ${les.minutes} min` : ""}`}),
        el("span", {class: "grow"}),
        el("button", {class: "tb dock-tool", title: "Reset lesson: fresh workspace and model", "aria-label": "Reset lesson", text: "⟲", onclick: () => resetLesson()}),
        isFull() ? btn("Side panel", () => L.setLayout("side"), {class: "btn small dock-layout", title: "show the lesson as a panel beside the plots and listings"})
          : btn("Full width", () => L.setLayout("full"), {class: "btn small dock-layout", title: "show the lesson full width, in a tab of its own"}),
        isFull() ? null : el("button", {class: "tb dock-tool", title: "Collapse the lesson panel", "aria-label": "Collapse the lesson panel", text: "»", onclick: () => L.toggleCollapse()}),
        el("button", {class: "tb dock-tool", title: "Close the lesson (Learn > Continue reopens it)", "aria-label": "Close the lesson", text: "×", onclick: () => L.closeLesson()})),
      el("div", {class: "dock-title", text: les.title}),
      el("div", {class: "lprog", title: `${done} of ${n} steps completed`}, el("div", {class: "bar", style: {width: `${n ? Math.round(100 * done / n) : 0}%`}})),
      pills);
    const expand = el("button", {class: "dock-expand", title: "Show the lesson panel", "aria-label": "Show the lesson panel", onclick: () => L.toggleCollapse()},
      el("span", {text: "«"}), el("span", {class: "dock-expand-t", text: les.title}));
    const body = el("div", {class: "dock-body"});
    wireDocLinks(body);
    if (k === 0) body.appendChild(introPage(les));
    else body.appendChild(stepPage(les, k));
    d.append(expand, head, body, footer(les, k), el("div", {class: "dock-status " + P.statusKind, "aria-live": "polite", title: P.status, text: P.status}));
    paintLines();
  }
  function introPage(les) {
    const pre = (les.prerequisites || []).map((id) => {
      const x = L.findLesson(id);
      return x ? btn(`${L.lessonNumber(id)}. ${x.title}`, () => L.openLesson(id), {class: "btn small linkbtn"}) : el("code", {text: id});
    });
    const intro = html("doc lesson-md", les.intro_html);
    fillBlocks(intro, les.intro_blocks, "setup");
    const ws = (les.session && les.session.workspace) || les.workspace;
    return el("div", {class: "lesson-page"},
      el("div", {class: "step-kicker", text: "Introduction"}),
      el("h2", {text: les.title}),
      les.summary ? el("p", {class: "lesson-summary", text: les.summary}) : null,
      les.objectives && les.objectives.length ? el("section", {class: "lpanel lpanel-obj"}, el("h4", {text: "You will be able to"}), el("ul", {}, ...les.objectives.map((o) => el("li", {text: o})))) : null,
      el("div", {class: "lesson-meta"},
        les.minutes ? el("span", {text: `About ${les.minutes} min`}) : null,
        les.example ? el("span", {}, "Example ", el("code", {text: les.example})) : null,
        pre.length ? el("span", {class: "lesson-pre"}, "Builds on: ", ...pre) : null,
        btn("▶ Watch the explainer video", () => L.watchVideo(les.id), {class: "btn small", title: "a narrated motion-graphics explainer of this lesson, in a new browser tab"})),
      intro,
      el("p", {class: "lesson-ws"}, "Workspace: ", el("code", {text: ws})),
      el("section", {class: "lpanel lpanel-how"}, el("h4", {text: "How to use the lesson panel"}),
        el("ul", {},
          el("li", {text: "Read the step, then press Run step: its commands go to the interpreter exactly as if you typed them. Module runs stream their listing into an output tab."}),
          el("li", {text: "Click a command line (or its ?) to see every argument explained; ✎ copies it into Command Entry for your own variation."}),
          el("li", {text: "After a step, its buttons open the plots and listings it produced. The numbered pills jump between steps; a green pill is a completed step."}))));
  }
  function stepPage(les, k) {
    const st = les.steps[k - 1];
    const narrative = html("doc lesson-md", st.narrative_html);
    fillBlocks(narrative, st.blocks, k);
    const page = el("div", {class: "lesson-page"},
      el("div", {class: "step-kicker", text: `Step ${k} of ${les.steps.length}`}),
      el("h2", {text: st.title}),
      narrative);
    if (st.actions.length) page.appendChild(actionsBar(st, k));
    for (const sec of st.sections) {
      const body = html("doc lesson-md", sec.html);
      fillBlocks(body, st.blocks, k);
      const panel = el("section", {class: `lpanel lpanel-${sec.key}`}, el("h4", {text: sec.title}), body);
      if (sec.key === "check" && sec.answer_html) {
        const ans = html("doc lesson-md lesson-answer", sec.answer_html);
        ans.hidden = true;
        const b = btn("Show answer", () => { ans.hidden = !ans.hidden; b.textContent = ans.hidden ? "Show answer" : "Hide answer"; });
        panel.append(b, ans);
      }
      page.appendChild(panel);
    }
    return page;
  }
  function actionsBar(st, k) {
    const ran = P.ran.has(k) || !st.commands.length;
    const bar = el("div", {class: "lesson-actions" + (ran ? " ready" : "")},
      el("div", {class: "k", text: ran ? "See the results" : "After running this step"}));
    for (const a of st.actions) {
      const b = btn(a.label, (ev) => doAction(st, k, a, ev.currentTarget), {class: "btn small" + (ran ? " primary" : ""), title: `${a.verb}${a.args ? ": " + a.args : ""}`});
      if (!a.known) b.disabled = true;
      bar.appendChild(b);
    }
    return bar;
  }
  const NEEDS_RUN = ["plot-model", "plot-nodes", "plot-soil", "plot-layers", "plot-soilprops", "plot-spectrum", "plot-history", "open-file", "open-listing", "animate"];
  async function doAction(st, k, a, anchor) {
    if (NEEDS_RUN.includes(a.verb) && st.commands.length && !P.ran.has(k)) {
      const yes = await D().confirm("Run the step first?", `"${a.label}" shows what step ${k} produces, and the step has not run in this workspace yet.`, "Run the step, then show", "Show anyway");
      if (yes) { const ok = await runStep(k); if (!ok) return; }
    }
    let r;
    try { r = await S.post(`/api/lessons/${enc(P.id)}/action`, {verb: a.verb, args: a.args}); }
    catch (e) { setStatus(e.message, "err"); S.local("WARNING", `${a.label}: ${e.message}`); return; }
    switch (r.kind) {
      case "commands": {
        const res = await S.command(r.lines).catch((e) => ({ok: false, messages: [{kind: "ERROR", text: e.message}]}));
        if (!res.ok) setStatus(`${a.label}: ${(res.messages || []).filter((m) => m.kind === "ERROR").map((m) => m.text).join("; ") || "failed"}`, "err");
        break;
      }
      case "file": S.openEditor(r.path); break;
      case "dialog": D().optionsDialog(r.name, r.tab || undefined); break;
      case "doc": D().helpTab({doc: r.doc, anchor: r.anchor || undefined}); break;
      case "explain": L.explain(r.line, anchor); break;
      default: break;
    }
  }
  /** Fill the command-block placeholders of a rendered piece with the interactive blocks. */
  function fillBlocks(root, blocks, stepKey) {
    const byK = {};
    for (const b of blocks || []) byK[b.k] = b;
    // the runnable lines are numbered in document order, as the step's commands
    let idx = 0;
    const order = (blocks || []).filter((b) => b.kind === "sassi").map((b) => b.k);
    const start = {};
    for (const kk of order) { start[kk] = idx; idx += byK[kk].lines.length; }
    root.querySelectorAll(".lesson-block").forEach((ph) => {
      const b = byK[Number(ph.dataset.block)];
      if (!b) { ph.remove(); return; }
      const label = b.kind === "sassi" ? "Commands of this step" : b.kind === "sassi-setup" ? "Setup (ran when the lesson opened)" : "Syntax / variation (not run by Run step)";
      const box = el("div", {class: `cmdblock cb-${b.kind}`}, el("div", {class: "cb-label", text: label}));
      b.lines.forEach((ln, i) => {
        const isComment = /^\s*[*!]/.test(ln);
        const code = el("code", {class: "cb-text", text: ln});
        const row = el("div", {class: "cmdline" + (isComment ? " comment" : ""), title: isComment ? "" : "click: explain this command"},
          el("span", {class: "cb-st", "aria-hidden": "true"}), code);
        if (b.kind === "sassi") { row.dataset.cmd = String(start[b.k] + i); row.dataset.step = String(stepKey); }
        if (b.kind === "sassi-setup") { row.dataset.cmd = String(i); row.dataset.step = "setup"; }
        if (!isComment) {
          row.append(el("button", {class: "cb-btn", title: "Explain this command", "aria-label": `Explain ${ln}`, text: "?", onclick: (ev) => { ev.stopPropagation(); L.explain(ln, row); }}),
            el("button", {class: "cb-btn", title: "Copy into Command Entry (edit, then Enter)", "aria-label": `Copy ${ln} into Command Entry`, text: "✎", onclick: (ev) => {
              ev.stopPropagation();
              const inp = S.$("#cmd");
              inp.value = ln;
              inp.focus();
              inp.dispatchEvent(new Event("input"));
            }}));
          row.addEventListener("click", () => { if (!String(window.getSelection ? window.getSelection() : "").trim()) L.explain(ln, row); });
        }
        box.appendChild(row);
      });
      ph.replaceWith(box);
    });
    // concept figures (```figure blocks): static/figures.js draws them; the caption is server-rendered HTML
    if (S.Figures) S.Figures.fill(root, blocks, (b) => html("doc lesson-md lfig-caption", b.caption_html), `${P.id}/${stepKey}`);
  }
  const ST_ICON = {running: "▶", job: "⟳", ok: "✓", error: "✗", skipped: "–"};
  function paintLines() {
    const d = dock();
    if (!d) return;
    d.querySelectorAll(".cmdline[data-cmd]").forEach((row) => {
      const key = row.dataset.step;
      const s = (P.lineState[key] || [])[Number(row.dataset.cmd)];
      row.classList.remove("st-running", "st-job", "st-ok", "st-error", "st-skipped");
      const icon = row.querySelector(".cb-st");
      if (!s) { if (icon) icon.textContent = ""; row.removeAttribute("data-info"); return; }
      row.classList.add("st-" + s.st);
      if (icon) icon.textContent = ST_ICON[s.st] || "";
      if (s.st === "error" && s.info) row.setAttribute("data-info", s.info);
      else if (s.st === "job") row.setAttribute("data-info", S.web ? "module run: see its output tab" : "running in a worker process: see its output tab");
      else row.removeAttribute("data-info");
    });
  }
  function footer(les, k) {
    const n = les.steps.length;
    const st = k > 0 ? les.steps[k - 1] : null;
    const foot = el("div", {class: "dock-foot"});
    if (P.running) {
      foot.appendChild(el("div", {class: "run-row"}, el("span", {class: "run-busy", text: "Running …"}), btn("Stop", () => { L.stop(); setStatus("Stopping after the current command …", "run"); }, {class: "btn small"})));
    } else if (st && st.commands.length) {
      const pending = les.steps.slice(0, k).filter((s, i) => s.commands.length && !P.ran.has(i + 1)).length;
      foot.appendChild(el("div", {class: "run-row"},
        btn(P.ran.has(k) ? "▶ Run step again" : "▶ Run step", () => runStep(k), {class: "btn primary", title: "run this step's commands in order (as if typed in Command Entry)"}),
        k > 1 ? btn(`Run steps 1–${k}`, () => runUpTo(k), {class: "btn", disabled: pending === 0, title: pending ? `run the ${pending} step(s) up to here that have not run in this workspace yet, in order` : "every step up to here has run in this workspace"}) : null));
    }
    const prev = btn("◀ Previous", () => goPage(k - 1), {class: "btn", disabled: k === 0});
    const next = k < n ? btn(k === 0 ? "Start ▶" : "Next ▶", () => { if (k > 0 && st && !st.commands.length) L.markDone(P.id, k); goPage(k + 1); }, {class: k === 0 ? "btn primary" : "btn"})
      : btn("Finish ✓", () => { L.markDone(P.id, k); setStatus("Lesson complete. The Learn page lists the next lesson.", "ok"); L.openStart(); renderPlayer(); }, {class: "btn"});
    foot.appendChild(el("div", {class: "nav-row"}, prev, el("span", {class: "grow nav-pos", text: k === 0 ? "Introduction" : `${k} / ${n}`}), next));
    return foot;
  }
  async function runStep(k) {
    if (P.running) return false;
    const st = P.lesson.steps[k - 1];
    if (!st.commands.length) { L.markDone(P.id, k); return true; }
    P.running = true;
    P.lineState[k] = [];
    setStatus(`Step ${k}: running ${st.commands.length} command(s) …`, "run");
    renderPlayer();
    const ok = await keepingFocus(() => L.runLines(st.commands, {onLine: (i, s, info) => { P.lineState[k][i] = {st: s, info}; paintLines(); if (s === "job" && info) setStatus(`Step ${k}: ${info.line} ${S.web ? "running (output tab)" : "running in a worker process (output tab, Cancel in the status bar)"}`, "run"); }}));
    P.running = false;
    if (ok) {
      P.ran.add(k);
      L.markDone(P.id, k);
      S.post(`/api/lessons/${enc(P.id)}/progress`, {step: k}).catch(() => {});
      setStatus(`Step ${k} done.` + (st.actions.length ? " Its buttons show the results." : ""), "ok");
    } else {
      setStatus(L.stopRequested ? `Step ${k} stopped.` : `Step ${k} stopped at an error (marked ✗; the Command History has the message).`, "err");
    }
    if (P.page === k || !L._batch) renderPlayer();
    return ok;
  }
  async function runUpTo(k) {
    if (P.running) return;
    L._batch = true;
    try {
      for (let j = 1; j <= k; j++) {
        const st = P.lesson.steps[j - 1];
        if (!st.commands.length || P.ran.has(j)) { if (!st.commands.length) L.markDone(P.id, j); continue; }
        if (P.page !== j) goPage(j);
        const ok = await runStep(j);
        if (!ok) return;
      }
      if (P.page !== k) goPage(k);
      setStatus(`Steps 1–${k} have run.`, "ok");
    } finally { L._batch = false; }
  }
  async function resetLesson() {
    if (P.running) return;
    const yes = await D().confirm("Reset lesson", `Delete the lesson workspace and start again from a fresh model?\n\n${(P.lesson.session && P.lesson.session.workspace) || P.lesson.workspace}`, "Reset", "Cancel");
    if (!yes) return;
    await L.openLesson(P.id, {fresh: true, step: 0});
  }

  // ================================================================== menu, start-up
  /** Items of the Learn menu (app.js builds the menu bar with them). */
  L.menuItems = function () {
    const items = [{label: "Start Page", action: () => L.openStart(), tip: "welcome, how SASSI works, the course, the examples"}];
    const cont = P.lesson ? P.id : L.lastLesson();
    if (cont && L.findLesson(cont)) {
      const x = L.findLesson(cont);
      items.push({label: `Continue: ${L.lessonNumber(cont)}. ${x.title.length > 40 ? x.title.slice(0, 38).replace(/[\s,:;]+$/, "") + "…" : x.title}`, tip: x.title, action: () => {
        if (P.lesson && P.id === cont) { showDock(); renderPlayer(); } else L.openLesson(cont, {step: L.progress(cont).last});
      }});
    }
    items.push("-");
    const o = L.outline || {parts: []};
    let n = 0;
    const parts = o.parts.slice().sort((a, b) => PART_ORDER.indexOf(a.part) - PART_ORDER.indexOf(b.part));
    for (const part of parts) {
      items.push({label: part.part, sub: part.lessons.map((les) => {
        n += 1;
        const p = L.progress(les.id);
        const done = p.done.filter((k) => k >= 1 && k <= les.nsteps).length;
        const t = les.title.length > 46 ? les.title.slice(0, 44).replace(/[\s,:;]+$/, "") + "…" : les.title;
        return {label: `${n}. ${t}`, tip: `${n}. ${les.title}`, tag: les.nsteps ? `${done}/${les.nsteps}` : "",
          action: () => L.openLesson(les.id, {step: p.last && done < les.nsteps ? p.last : 0})};
      })});
    }
    if (!parts.length) items.push({label: "(no lessons installed)", disabled: true});
    items.push("-",
      {label: "How SASSI Works", action: () => L.openStart({section: "how"})},
      {label: "Explainer Videos", action: () => L.watchVideo(null), tip: "one narrated motion-graphics video per lesson (new browser tab)"},
      {label: "Examples Gallery", action: () => L.openStart({section: "examples"})},
      {label: "Open Example...", action: () => L.examplesDialog()},
      {label: "Explain a Command...", action: () => L.askExplain(), tip: "the arguments of a command line explained"},
      "-",
      {label: "Lesson Layout: " + (isFull() ? "Full Width ✓" : "Side Panel ✓"), tip: "switch between a full-width lesson tab and a panel beside the results",
        action: () => L.setLayout(isFull() ? "side" : "full")},
      {label: "Course Workspace Folder...", action: () => L.workspaceDialog()},
      {label: "Free Disk Space...", action: () => L.freeSpaceDialog()},
      {label: "Reset Course Progress", action: async () => {
        if (await D().confirm("Reset Course Progress", "Forget the completed steps of every lesson in this browser? (The workspaces on disk are not touched.)", "Reset", "Cancel")) { L.resetProgress(); if (S.rebuildMenus) S.rebuildMenus(); renderStart(); }
      }});
    return items;
  };

  /** Called by S.start once the state is loaded. */
  L.init = async function (opts) {
    opts = opts || {};
    setupHint();
    setupHistoryClicks();
    setupResizer();
    applyDockWidth();
    window.addEventListener("resize", applyDockWidth);
    await L.load();
    // the start page opens with the GUI (unless switched off on it); it comes to the front only when no
    // plot was restored (a page reload keeps the learner's plots in front)
    if (store.get("autostart", true) !== false) await L.openStart({select: !opts.hadPlots});
    const ls = S.state && S.state.lesson;
    if (ls && ls.id && L.findLesson(ls.id)) L.openLesson(ls.id, {step: L.progress(ls.id).last});
  };
})(SASSI);
