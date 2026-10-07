/* SASSI-EDU GUI -- the run summary window: after a run, the key inputs, the key outputs and their graphs of the
 * model (GET /api/run_summary, sassi/ui/runsummary.py: read from the decks the modules used and the result files
 * they wrote).
 *
 * It opens by itself when a run that ran modules ends (an input file, an example's Run all, a module from the
 * Modules menu or Command Entry; the activity panel reports it: S.onRunFinished) -- not while a lesson is open,
 * which has its own result buttons, not after a single intermediate module (SITE, POINT, HOUSE, FORCE, ANALYS),
 * and not when "Show after every run" is off (kept by the browser).  View > Run Summary opens it at any time.  "Open as plot" turns a graph into a plot tab
 * with command text (READSPEC / READTH and SPECPLOT / THPLOT, rule L17). */
"use strict";

const SASSI_RUNSUMMARY = (function () {
  const S = SASSI, el = S.el;
  const R = {};
  const KEY = "sassi-edu.runSummary";
  R.auto = () => { try { return localStorage.getItem(KEY) !== "0"; } catch (e) { return true; } };
  R.setAuto = (on) => { try { localStorage.setItem(KEY, on ? "1" : "0"); } catch (e) { /* private mode */ } };
  let current = null;                           // the open window

  function table(sec) {
    const body = el("tbody");
    for (const [k, v] of sec.rows) body.appendChild(el("tr", {}, el("th", {text: k}), el("td", {text: v})));
    return el("section", {class: "rs-sec"}, el("h4", {text: sec.title}), el("table", {class: "rs-table"}, body));
  }
  function moduleChips(mods) {
    const box = el("div", {class: "rs-chips"});
    for (const m of mods) {
      const bad = !m.ok || (m.errors || 0) > 0;
      const t = m.seconds === null || m.seconds === undefined ? "" : ` ${m.seconds < 10 ? m.seconds.toFixed(1) : m.seconds.toFixed(0)} s`;
      const w = m.warnings ? ` · ${m.warnings} warning${m.warnings > 1 ? "s" : ""}` : "";
      box.appendChild(el("span", {class: "rs-chip" + (bad ? " bad" : m.warnings ? " warn" : ""), title: `${m.listing}${w}`,
        text: `${bad ? "✗" : "✓"} ${m.module}${t}${w}`}));
    }
    return box;
  }
  function chartCard(ch, d, api) {
    const plot = el("div", {class: "rs-plot"});
    const head = el("div", {class: "rs-chart-head"}, el("span", {class: "grow", text: ch.title}));
    if (ch.kind && ch.files && ch.files.length) {
      head.appendChild(el("button", {class: "btn small", text: "Open as plot", title: `${ch.kind === "spec" ? "READSPEC + SPECPLOT" : "READTH + THPLOT"} of ${ch.files.join(", ")}`,
        onclick: async () => {
          try {
            const r = await S.post("/api/run_summary/plot", {model: d.model.name || d.model.number, chart: ch.id});
            await S.command(r.lines);
            api.close();
          } catch (e) { api.setMessage(e.message); }
        }}));
    }
    const card = el("div", {class: "rs-chart"}, head, plot);
    card._draw = () => {
      if (!window.Plotly) { plot.textContent = "(Plotly is not loaded)"; return; }
      const traces = ch.lines.map((ln) => ({type: "scatter", mode: "lines", name: ln.name, x: ln.x, y: ln.y,
        line: Object.assign({width: ln.dash ? 1.4 : 1.6, dash: ln.dash ? "dash" : "solid"}, ln.color ? {color: ln.color} : {}), hovertemplate: `${S.esc(ln.name)}: %{y:.4g} at %{x:.4g}<extra></extra>`}));
      const layout = {height: 250, margin: {l: 58, r: 10, t: 6, b: 40}, paper_bgcolor: "#fff", plot_bgcolor: "#fff",
        font: {size: 11}, showlegend: true, legend: {orientation: "h", y: -0.32, font: {size: 10}},
        xaxis: {title: {text: ch.x}, type: ch.logx ? "log" : "linear", gridcolor: "#eceff3", zeroline: false,
          range: ch.range && !ch.logx ? ch.range : undefined},
        yaxis: {title: {text: ch.y}, gridcolor: "#eceff3", zeroline: false, autorange: ch.reversey ? "reversed" : true}};
      window.Plotly.newPlot(plot, traces, layout, {displaylogo: false, responsive: true, modeBarButtonsToRemove: ["select2d", "lasso2d"]});
    };
    return card;
  }

  /** Open the window for a model (name or number; default: the active model). */
  R.open = async function (opts) {
    opts = opts || {};
    let d;
    try { d = await S.get("/api/run_summary" + (opts.model !== undefined && opts.model !== null ? `?model=${encodeURIComponent(opts.model)}` : "")); }
    catch (e) { S.local("WARNING", `Run Summary: ${e.message}`); return null; }
    if (current) current.close();
    const m = d.model || {};
    const head = el("div", {class: "rs-head"},
      el("div", {class: "rs-model"}, el("b", {text: `Model ${m.number}${m.name ? " " + m.name : ""}`}), m.title ? el("span", {text: ` · ${m.title}`}) : null),
      moduleChips(d.modules || []));
    const others = (d.others || []).filter((o) => o.name);
    if (others.length > 1 || (opts.models && opts.models.length > 1)) {
      const sel = el("select", {"aria-label": "model"}, ...others.map((o) => el("option", {value: o.name, text: `${o.name}${o.title ? " · " + o.title : ""}`})));
      sel.value = m.name;
      sel.addEventListener("change", () => R.open({model: sel.value, models: opts.models}));
      head.appendChild(el("label", {class: "rs-pick"}, "Model ", sel));
    }
    const ins = el("div", {class: "rs-col"}, el("h3", {text: "Key inputs"}), ...(d.inputs || []).map(table));
    const outs = el("div", {class: "rs-col"}, el("h3", {text: "Key outputs"}), ...(d.outputs || []).map(table));
    if (!(d.outputs || []).length) outs.appendChild(el("p", {class: "opt-note", text: "No result files in the model folder yet."}));
    const charts = el("div", {class: "rs-charts"});
    const autoBox = el("input", {type: "checkbox"});
    autoBox.checked = R.auto();
    autoBox.addEventListener("change", () => R.setAuto(autoBox.checked));
    const body = el("div", {class: "rs-body"}, head,
      d.note ? el("p", {class: "opt-note", text: d.note}) : null,
      el("div", {class: "rs-grid"}, ins, outs),
      (d.charts || []).length ? el("h3", {class: "rs-graphs", text: "Graphs"}) : null, charts,
      el("label", {class: "chk rs-auto"}, autoBox, " Show this summary after every run (View > Run Summary opens it at any time)"));
    const api = S.D.modal({title: "Run Summary", width: "min(1180px, calc(100vw - 24px))", body, buttons: [
      {label: "Results Browser", action: () => { S.openResults(); return true; }},
      {label: "Close", primary: true, action: () => true}],
      onClose: () => { if (current === api) current = null; }});
    current = api;
    const cards = (d.charts || []).map((ch) => chartCard(ch, d, api));
    for (const c of cards) charts.appendChild(c);
    // after the window is laid out (a timer: an animation frame waits while the page is hidden); the window opens
    // at its top (the dialog focuses its first input, the checkbox at the bottom)
    setTimeout(() => {
      if (document.activeElement === autoBox) autoBox.blur();
      api.body.scrollTop = 0;
      for (const c of cards) { try { c._draw(); } catch (e) { console.error(e); } }
    }, 30);
    return api;
  };

  /** Modules whose results the summary is about: a single module run from the menu or Command Entry opens the
   *  window after one of these, not after SITE, POINT, HOUSE, FORCE or ANALYS (the steps on the way to them). */
  const RESULT_MODULES = ["MOTION", "STRESS", "RELDISP", "SOIL", "EQUAKE", "COMBIN", "NONLINEAR"];
  /** The activity panel reports the end of a run (static/activity.js). */
  S.onRunFinished = function (run) {
    if (!R.auto() || S.keepFocus || !run || !run.modules.length) return;
    // while a lesson is open its steps run the modules, and the lesson has its own result buttons (the run's end
    // is reported by the event poll, after the step: its running flag is no guide)
    if (SASSI.Learn && SASSI.Learn.isLessonOpen && SASSI.Learn.isLessonOpen()) return;
    if (run.kind === "job" && !run.modules.some((m) => RESULT_MODULES.includes(m))) return;
    // the first model the run ran (an example runs its main model first, then its variants), else the active one
    const model = run.models && run.models.length ? run.models[0] : undefined;
    setTimeout(() => R.open({model, models: run.models}), 400);   // after the last messages and the panel's update
  };
  S.openRunSummary = (model) => R.open({model});
  return R;
})();
SASSI.RunSummary = SASSI_RUNSUMMARY;
