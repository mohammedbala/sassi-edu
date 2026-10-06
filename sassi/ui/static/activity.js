/* SASSI-EDU GUI -- the activity panel: animated feedback while an input file runs (INP: Model > Input,
 * Run (INP) in the File Editor, an example's Run all, a lesson step) or a module runs (RUN<MODULE> typed or
 * from the Modules menu).
 *
 * The panel only shows what the interpreter reports (it sends no command; rule L17 is not affected).  It
 * listens to the session's events (app.js handleEvent / onJobEvent / setBusy):
 *  - progress {line, total, source, cmd, name, what}: the INP / macro file and the command it is executing
 *    (name: the canonical command, what: its one-line meaning from the command explainer);
 *  - progress {kind: "module", module, fraction, text}: a module run inside the file (ANALYS: frequency k/n);
 *  - message INFO "RUN<MODULE>: model ..." (a module starts), "<MODULE> finished with status <S> in <t> s; ..."
 *    (it ends: the listing's last line), "INP <file>: n lines, m commands, w warnings, e errors" (the file
 *    ends); WARNING / ERROR messages are counted;
 *  - job: a module run in a worker process, with its progress;
 *  - the end of the request (setBusy back to 0): a file that stopped early is closed as well.
 * Under prefers-reduced-motion the spinner and the pulsing are still (styles.css). */
const SASSI_ACTIVITY = (function () {
  "use strict";
  const S = SASSI;
  const el = S.el;
  const A = {};

  // plain-language phrases for the commands of a typical input file (else the explainer's meaning)
  const FRIENDLY = {
    MDL: "Naming the model and its folder", TIT: "Setting the model title", MUNITS: "Setting the mass units",
    GRAVITY: "Setting gravity", L: "Defining a soil layer type", TOPL: "Stacking the soil layers of the site",
    FREQ: "Choosing the analysis frequencies", N: "Defining nodes", NGEN: "Generating nodes", FILL: "Filling nodes in a line",
    E: "Defining elements", EGEN: "Generating elements", GROUP: "Starting an element group", GTIT: "Naming an element group",
    M: "Defining a material", MACT: "Selecting the material", R: "Defining section properties", RACT: "Selecting the section",
    MT: "Adding a lumped mass", MR: "Adding a rotary mass", D: "Fixing degrees of freedom", INT: "Marking interaction nodes",
    INTGEN: "Generating the interaction nodes", MSET: "Linking excavated soil to its layer", ETYPE: "Setting the element type",
    SITE: "Setting the SITE options", WAVE: "Defining the seismic wave field", POINT: "Setting the POINT options",
    HOUSE: "Setting the HOUSE options", ANALYS: "Setting the ANALYS options", MOTION: "Setting the MOTION options",
    STRESS: "Setting the STRESS options", RELD: "Choosing relative displacements", NOUT: "Choosing the output nodes",
    EOUT: "Choosing the output elements", AOPT: "Choosing the modules to run", CHECK: "Checking the model",
    AFWRITE: "Writing the module input files", CPMODEL: "Copying the model", ACTM: "Switching the active model",
    FCOPY: "Copying a result file", VAR: "Defining a list", FOREACH: "Repeating a command over a list",
    INP: "Reading another input file", CALCM: "Totalling the masses", HARMFRAME: "Writing animation frames",
    WINDOWSETTINGS: "Setting what the plots show", MODELPLOT: "Plotting the model", EQUAKE: "Setting the EQUAKE options",
    SOIL: "Setting the SOIL options", DYNP: "Defining a soil degradation curve", SAVE: "Saving the model",
  };
  // what each module computes (the "Now" line while it runs)
  const MODULES = {
    EQUAKE: "spectrum-compatible motion", SOIL: "strain-compatible soil properties (SHAKE)",
    SITE: "free field of the layered site", POINT: "soil point-load solutions", HOUSE: "structure and excavated-soil matrices",
    FORCE: "load vectors", ANALYS: "soil impedance and the SSI solution, frequency by frequency",
    COMBIN: "combination of analyses", MOTION: "time histories and in-structure spectra", STRESS: "element forces",
    RELDISP: "relative displacements", NONLINEAR: "nonlinear iterations", LOADGEN: "ANSYS load files",
  };
  const PARENTS = ["NONLINEAR"];
  const RE_RUN = /^RUN([A-Z]+): model /;
  const RE_DONE = /^\s*([A-Z]+) finished with status (\w+) in ([\d.]+) s/;
  const RE_INP = /^INP (.+): (\d+) lines, (\d+) commands, (\d+) warnings, (\d+) errors/;

  let box = null, ui = null;
  let run = null;          // {kind: "inp" | "job", source, t0, line, total, warnings, errors, modules: Map, current}
  let tick = null, hideTimer = null;

  function build() {
    if (box) return;
    ui = {
      spin: el("span", {class: "act-spin", "aria-hidden": "true"}),
      title: el("span", {class: "act-title"}),
      time: el("span", {class: "act-time"}),
      close: el("button", {class: "act-x", title: "Hide (the run goes on; the status bar still shows it)", "aria-label": "Hide the activity panel", text: "×"}),
      fill: el("div"),
      line: el("div", {class: "act-line"}),
      now: el("div", {class: "act-now"}),
      cmd: el("code", {class: "act-cmd"}),
      modFill: el("div"),
      modText: el("span", {class: "act-mod-t"}),
      chain: el("div", {class: "act-chain"}),
      counts: el("span", {class: "act-counts"}),
      hist: el("button", {class: "btn small", text: "Command History"}),
    };
    ui.mod = el("div", {class: "act-mod", hidden: true}, el("div", {class: "act-mod-bar"}, ui.modFill), ui.modText);
    box = el("div", {id: "activity", class: "activity", role: "status", "aria-live": "polite", hidden: true},
      el("div", {class: "act-head"}, ui.spin, ui.title, ui.time, ui.close),
      el("div", {class: "act-bar"}, ui.fill), ui.line,
      el("div", {class: "act-nowrow"}, el("span", {class: "act-k", text: "Now"}), ui.now), ui.cmd, ui.mod,
      ui.chain,
      el("div", {class: "act-foot"}, ui.counts, ui.hist));
    ui.close.addEventListener("click", () => hide(true));
    ui.hist.addEventListener("click", () => {
      S.selectTab("history");
      const first = run && run.errors ? document.querySelector(".m-ERROR") : null;
      if (first && first.scrollIntoView) first.scrollIntoView({block: "center"});
    });
    document.body.appendChild(box);
  }

  function show() {
    build();
    clearTimeout(hideTimer);
    box.hidden = false;
    box.classList.remove("done", "failed", "leaving");
    requestAnimationFrame(() => box.classList.add("in"));
  }
  function hide(now) {
    if (!box) return;
    clearTimeout(hideTimer);
    box.classList.add("leaving");
    box.classList.remove("in");
    hideTimer = setTimeout(() => { box.hidden = true; box.classList.remove("leaving"); }, now ? 200 : 400);
  }

  function start(kind, source) {
    run = {kind, source, t0: Date.now(), line: 0, total: 0, warnings: 0, errors: 0, modules: new Map(), current: null, dismissed: false};
    show();
    ui.chain.innerHTML = "";
    ui.mod.hidden = true;
    ui.cmd.textContent = "";
    box.classList.toggle("is-job", kind === "job");
    ui.title.textContent = kind === "inp" ? `Running ${source}` : `Running ${source}`;
    setNow(kind === "inp" ? "Reading the input file" : `Running ${source}`);
    paint();
    clearInterval(tick);
    tick = setInterval(paintTime, 500);
    paintTime();
  }
  function paintTime() {
    if (!run || !ui) return;
    const s = (Date.now() - run.t0) / 1000;
    ui.time.textContent = s < 60 ? `${s.toFixed(0)} s` : `${Math.floor(s / 60)} min ${String(Math.round(s % 60)).padStart(2, "0")} s`;
  }
  function setNow(text, cmd) {
    if (ui.now.textContent !== text) {
      ui.now.textContent = text;
      ui.now.classList.remove("flash");
      void ui.now.offsetWidth;              // restart the fade-in
      ui.now.classList.add("flash");
    }
    ui.cmd.textContent = cmd || "";
    ui.cmd.hidden = !cmd;
  }
  function paint() {
    if (!run || !ui) return;
    if (run.kind === "inp") {
      const f = run.total ? Math.min(1, run.line / run.total) : 0;
      ui.fill.style.width = `${(100 * f).toFixed(1)}%`;
      ui.line.textContent = run.total ? `line ${run.line} of ${run.total}` : "";
    } else {
      const m = run.current && run.modules.get(run.current);
      ui.fill.style.width = `${(100 * ((m && m.fraction) || 0)).toFixed(1)}%`;
      ui.line.textContent = "";
    }
    const w = run.warnings, e = run.errors;
    ui.counts.textContent = `${e} error${e === 1 ? "" : "s"} · ${w} warning${w === 1 ? "" : "s"}`;
    ui.counts.classList.toggle("bad", e > 0);
    ui.counts.classList.toggle("warn", e === 0 && w > 0);
  }

  // ---------------------------------------------------------------- modules (chips)
  function moduleStart(name) {
    if (!run) return;
    let m = run.modules.get(name);
    if (m && m.state === "run") { run.current = name; return; }   // the same run, announced again
    // a module that runs others (NONLINEAR iterates HOUSE, ANALYS ...) stays running as their parent
    if (run.current && run.current !== name && !PARENTS.includes(run.current)) finishCurrentModule("ok");
    if (!m) {
      const chip = el("span", {class: "act-chip run", title: MODULES[name] ? `${name}: ${MODULES[name]}` : name},
        el("span", {class: "act-chip-dot", "aria-hidden": "true"}), el("span", {class: "act-chip-n", text: name}), el("span", {class: "act-chip-s"}));
      ui.chain.appendChild(chip);
      m = {name, chip, t0: Date.now(), fraction: 0, state: "run", count: 1};
      run.modules.set(name, m);
    } else {                                // run again (NONLINEAR iterations, a second analysis)
      m.t0 = Date.now(); m.state = "run"; m.fraction = 0; m.count += 1;
      m.chip.className = "act-chip run";
      m.chip.querySelector(".act-chip-s").textContent = `×${m.count}`;
    }
    run.current = name;
    setNow(`Running ${name}${MODULES[name] ? ": " + MODULES[name] : ""}`, `RUN${name}`);
    ui.mod.hidden = false;
    ui.mod.classList.add("indet");             // until the module reports its progress (some report none)
    ui.modFill.style.width = "0%";
    ui.modText.textContent = "working ...";
  }
  function moduleProgress(name, fraction, text) {
    if (!run) return;
    const known = run.modules.get(name);
    if (known && known.state !== "run" && fraction >= 0.999) return;   // the final report of a finished run
    if (!known || known.state !== "run") moduleStart(name);
    const m = run.modules.get(name);
    m.fraction = fraction;
    ui.mod.classList.remove("indet");
    ui.modFill.style.width = `${(100 * fraction).toFixed(1)}%`;
    ui.modText.textContent = text || `${Math.round(100 * fraction)} %`;
    m.chip.querySelector(".act-chip-s").textContent = `${Math.round(100 * fraction)}%`;
    if (run.kind === "job") paint();
  }
  function moduleDone(name, state, seconds) {
    if (!run) return;
    const m = run.modules.get(name);
    if (!m) return;
    m.state = state;
    m.chip.className = `act-chip ${state}`;
    const s = seconds !== undefined ? seconds : (Date.now() - m.t0) / 1000;
    const times = m.count > 1 ? ` ×${m.count}` : "";
    m.chip.querySelector(".act-chip-s").textContent = (state === "ok" ? `✓ ${s < 10 ? s.toFixed(1) : s.toFixed(0)} s` : "✗") + times;
    if (run.current === name) { run.current = null; ui.mod.hidden = true; }
  }
  function finishCurrentModule(state) {
    if (run && run.current) moduleDone(run.current, state);
  }

  // ---------------------------------------------------------------- end of a run
  function finish(summary) {
    if (!run) return;
    finishCurrentModule(run.errors ? "fail" : "ok");
    clearInterval(tick);
    paintTime();
    const bad = run.errors > 0;
    box.classList.add(bad ? "failed" : "done");
    if (run.kind === "inp") {
      ui.fill.style.width = "100%";
      ui.title.textContent = `${bad ? "Finished with errors" : "Finished"}: ${run.source}`;
      const t = (Date.now() - run.t0) / 1000;
      ui.line.textContent = summary ? `${summary.commands} commands in ${t.toFixed(1)} s` : `${t.toFixed(1)} s`;
      setNow(bad ? "Some commands failed: see the Command History" : "The input file has run", "");
    } else {
      ui.title.textContent = `${bad ? "Failed" : "Finished"}: ${run.source}`;
      setNow(bad ? "The run failed: see its output tab" : "The module has run", "");
    }
    ui.mod.hidden = true;
    paint();
    const ended = run;
    run = null;
    // a clean run fades out on its own; with errors the panel stays until it is closed
    if (!bad) hideTimer = setTimeout(() => { if (!run && ended) hide(); }, 4500);
  }

  // ---------------------------------------------------------------- the session's events
  A.onEvent = function (ev) {
    if (!ev) return;
    if (ev.type === "progress") {
      if (ev.kind === "module") {
        // only inside a run the panel already shows (a late report after the file's end opens nothing)
        if (run) moduleProgress(String(ev.module || ""), Number(ev.fraction) || 0, ev.text);
        return;
      }
      if (ev.source === undefined || ev.total === undefined) return;
      if (!run || run.kind !== "inp") start("inp", ev.source);
      if (ev.source !== run.source) {                       // a nested INP / macro
        if (ev.cmd) setNow(`${ev.source}: ${describe(ev)}`, ev.cmd);
        return;
      }
      run.line = ev.line;
      run.total = ev.total;
      if (ev.name && ev.name.startsWith("RUN") && MODULES[ev.name.slice(3)] !== undefined) moduleStart(ev.name.slice(3));
      else if (ev.cmd && !(run.current && ev.name === "")) setNow(describe(ev), ev.name ? ev.cmd : "");
      paint();
      return;
    }
    if (ev.type !== "message" || !run) return;
    const text = String(ev.text || "");
    if (ev.kind === "WARNING") run.warnings++;
    else if (ev.kind === "ERROR") run.errors++;
    else if (ev.kind === "INFO") {
      let m = RE_RUN.exec(text);
      if (m && MODULES[m[1]] !== undefined) { moduleStart(m[1]); paint(); return; }
      m = RE_DONE.exec(text);
      if (m && run.modules.has(m[1])) { moduleDone(m[1], /^OK$/i.test(m[2]) ? "ok" : "fail", Number(m[3])); paint(); return; }
      m = RE_INP.exec(text);
      if (m && run.kind === "inp" && m[1] === run.source) {
        run.warnings = Math.max(run.warnings, Number(m[4]));
        run.errors = Math.max(run.errors, Number(m[5]));
        finish({lines: Number(m[2]), commands: Number(m[3])});
        return;
      }
    }
    if (ev.kind === "WARNING" || ev.kind === "ERROR") paint();
  };
  function describe(ev) {
    if (FRIENDLY[ev.name]) return FRIENDLY[ev.name];
    if (ev.what) return ev.what.replace(/\s*\((requirements|D-|spec)[^)]*\)\s*\.?$/i, "").replace(/\.$/, "");
    if (!ev.name) return "Reading comments";
    return ev.name;
  }
  /** A module run in a worker process (app.js onJobEvent). */
  A.onJob = function (job) {
    if (!job || job.kind !== "module") return;
    const active = ["starting", "running"].includes(job.state);
    if (active) {
      if (!run) start("job", job.module || job.line);
      if (run.kind === "job") {
        if (!run.modules.has(job.module) || run.current !== job.module) moduleStart(job.module);
        if (job.progress) moduleProgress(job.module, Number(job.progress) || 0, job.progress_text || "");
      }
    } else if (run && run.kind === "job") {
      if (job.ok === false) run.errors = Math.max(run.errors, 1);
      moduleDone(job.module, job.ok === false ? "fail" : "ok");
      finish(null);
    }
  };
  /** The request has ended (app.js setBusy back to 0): close an INP that stopped without its summary line. */
  A.onIdle = function () {
    if (run && run.kind === "inp") setTimeout(() => { if (run && run.kind === "inp" && !(S.state && S.state.job)) finish(null); }, 300);
  };
  return A;
})();
SASSI.activity = SASSI_ACTIVITY;
