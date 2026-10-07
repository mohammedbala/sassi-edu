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
 *  - progress {kind: "step", module, key, data}: the computation step a module enters (ctx.announce in
 *    sassi/modules: ANALYS.impedance ...) with its sizes: the panel shows the step's equation (STEPS, the
 *    notation of the Theory Manual, with a link to its section) and the numbers -- the matrices being formed,
 *    inverted or factorised and their dimensions;
 *  - job: a module run in a worker process, with its progress (and job.step);
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

  // ---------------------------------------------------------------- the computation steps (the math being done)
  const n = (v) => (v === undefined || v === null ? "?" : Number(v).toLocaleString("en-US"));
  const hz = (f) => `${S.fmt(Number(f), 4)} Hz`;
  const freq = (d) => `frequency ${d.q}/${d.nF} · ${hz(d.f)}`;
  const TH = "docs/theory/THEORY_MANUAL.md";
  /** key -> {short (the step strip), title, tex (string or function of the data), detail (function), ref
   *  (Theory Manual section, anchor)}.  The modules report the keys and the sizes (sassi/modules: ctx.announce). */
  const STEPS = {
    "SITE.modes": {short: "Wave modes", title: "Wave modes of the layered soil (thin-layer method)", ref: ["4.3", "43-eigenproblems"],
      tex: String.raw`\left(k^2\mathbf{A} + k\,\mathbf{B} + \mathbf{C} - \omega^2\mathbf{M}\right)\boldsymbol{\phi} = \mathbf{0}`,
      detail: (d) => `${freq(d)} · ${n(d.n)} sublayers (${n(d.ngen)} generated for the half-space, ${d.base} base) → eigenproblems in k² for the Rayleigh and Love wavenumbers`},
    "SITE.field": {short: "Free field", title: "Free-field motion, normalised to the control point", ref: ["6.2", "62-normalisation-to-the-control-point"],
      tex: String.raw`\mathbf{U}'(z,\omega) = \frac{\mathbf{u}(z,\omega)}{u_{cp}(\omega)}`,
      detail: (d) => `${d.wave} wave · ${freq(d)} · motion at ${n(d.nI)} layer interfaces, unit motion at the top of layer ${d.cl} in ${d.cm}`},
    "POINT.green": {short: "Point loads", title: "Point-load solutions: central zone plus outgoing wave modes", ref: ["7.1", "71-what-point-computes"],
      tex: String.raw`\tilde{\mathbf{u}}(\rho) = \boldsymbol{\Psi}_\mu(\rho)\,\boldsymbol{\alpha},\qquad \rho \ge R_0`,
      detail: (d) => `${freq(d)} · ${n(d.nR)} Rayleigh + ${n(d.nLove)} Love modes · loads at ${n(d.nL)} interface(s), central zone R₀ = ${S.fmt(d.R0, 4)}`},
    "HOUSE.elements": {short: "Element matrices", title: "Element stiffness and mass matrices", ref: ["13.1", "131-the-element-library"],
      tex: String.raw`\mathbf{K}_e = \int_{V_e}\mathbf{B}^{\mathsf T}\mathbf{D}\,\mathbf{B}\,dV,\qquad \mathbf{M}_e = \int_{V_e}\rho\,\mathbf{N}^{\mathsf T}\mathbf{N}\,dV`,
      detail: (d) => `${n(d.ne)} elements (${d.dim === 1 ? "2D" : "3D"})`},
    "HOUSE.assembled": {short: "Assembly", title: "Assembly: complex (hysteretic) stiffness and mass", ref: ["2.1", "21-the-sassi-form"],
      tex: String.raw`\mathbf{K}^{*} = \sum_e \mathbf{K}_e\left(1 - 2\beta^{2} + 2i\beta\sqrt{1-\beta^{2}}\right),\qquad \mathbf{M} = \sum_e \mathbf{M}_e`,
      detail: (d) => `${n(d.neq)} equations · ${n(d.nnz)} non-zeros (structure)${d.nnze ? ` · ${n(d.nnze)} (excavated soil)` : ""} · ${n(d.nint)} interaction nodes`},
    "FORCE.loads": {short: "Loads", title: "Load vectors of the foundation vibration", ref: ["9.1", "91-per-frequency"],
      tex: String.raw`P_k(\omega) = a_k\,e^{-i\omega t_k}`,
      detail: (d) => `${n(d.nloads)} loaded DOF(s) on ${n(d.nnodes)} node(s) at ${n(d.nF)} frequencies`},
    "ANALYS.impedance": {short: "Impedance", title: "Soil impedance at the interaction nodes", ref: ["8.2", "82-impedance"],
      tex: String.raw`\mathbf{X}_{ff}(\omega) = \mathbf{F}_{ff}(\omega)^{-1}`,
      detail: (d) => `${freq(d)} · ${d.restart ? "read from the restart file" : `dense ${n(d.n)} × ${n(d.n)} complex inverse of the POINT flexibility`}`},
    "ANALYS.dynamic": {short: "Dynamic stiffness", title: "Dynamic stiffness: structure minus excavated soil", ref: ["9.1", "91-per-frequency"],
      tex: String.raw`\mathbf{C}(\omega) = \left(\mathbf{K}^{*}_{s} - \omega^{2}\mathbf{M}_{s}\right) - \left(\mathbf{K}^{*}_{e} - \omega^{2}\mathbf{M}_{e}\right)`,
      detail: (d) => `${freq(d)} · sparse ${n(d.neq)} × ${n(d.neq)} (${n(d.nnz)} non-zeros in K*ₛ)`},
    "ANALYS.factor": {short: "Condense + factorise", title: "Condense onto the interaction DOFs (Schur complement) and factorise", ref: ["9.2", "92-schur-complement"],
      tex: String.raw`\mathbf{S} = \mathbf{C}_{ff} - \mathbf{C}_{fn}\,\mathbf{C}_{nn}^{-1}\,\mathbf{C}_{nf},\qquad \mathrm{LU}\!\left(\mathbf{S} + \mathbf{X}_{ff}\right)`,
      detail: (d) => `${freq(d)} · sparse LU of C_nn (${n(d.nn)} DOFs) · dense LU of S + X (${n(d.nf)} × ${n(d.nf)})`},
    "ANALYS.solve": {short: "Solve", title: "Load and solution of the SSI equation", ref: ["9.2", "92-schur-complement"],
      tex: (d) => d.seismic
        ? String.raw`\mathbf{b}_f = \mathbf{X}_{ff}\,\mathbf{U}'_f,\qquad \mathbf{U}_f = \left(\mathbf{S} + \mathbf{X}_{ff}\right)^{-1}\!\left(\mathbf{b}_f - \mathbf{C}_{fn}\mathbf{C}_{nn}^{-1}\mathbf{b}_n\right)`
        : String.raw`\begin{bmatrix}\mathbf{C}_{nn} & \mathbf{C}_{nf}\\ \mathbf{C}_{fn} & \mathbf{C}_{ff} + \mathbf{X}_{ff}\end{bmatrix}\begin{bmatrix}\mathbf{U}_n\\ \mathbf{U}_f\end{bmatrix} = \mathbf{P}`,
      detail: (d) => `${freq(d)} · ${n(d.ncases)} load case(s) · ${n(d.neq)} transfer functions ${d.seismic ? "per unit control motion" : "per unit load"}${d.path && d.path !== "schur" ? ` (path: ${d.path})` : ""}`},
    "MOTION.fft": {short: "FFT", title: "Fourier transform of the input", ref: ["11.1", "111-convolution"],
      tex: String.raw`A(\omega_k) = \sum_{j=0}^{N-1} a(t_j)\,e^{-i\omega_k t_j}`,
      detail: (d) => `${d.seismic ? "control motion" : "load history"} padded to N = ${n(d.nfft)} points, Δt = ${S.fmt(d.dt, 4)} s`},
    "MOTION.interp": {short: "Interpolate", title: "Transfer functions at every Fourier frequency", ref: ["10", "10-transfer-function-interpolation"],
      tex: String.raw`H(\omega_k) = \mathcal{I}\left\{H(\omega_1), \ldots, H(\omega_{N_f})\right\}`,
      detail: (d) => `${n(d.nF)} SSI frequencies → ${n(d.nK)} Fourier frequencies · ${n(d.mb)} node-directions (option ${d.option})`},
    "MOTION.convolve": {short: "Convolve", title: "Response histories by inverse Fourier transform", ref: ["11.1", "111-convolution"],
      tex: (d) => d.seismic ? String.raw`a_j(t) = \mathcal{F}^{-1}\left\{H_j(\omega)\,A(\omega)\right\}`
        : String.raw`u_j(t) = \mathcal{F}^{-1}\left\{H_j(\omega)\,F(\omega)\right\},\quad a_j(t) = \mathcal{F}^{-1}\left\{-\omega^2 H_j F\right\}`,
      detail: (d) => `${n(d.mb)} histories of ${n(d.nfft)} points`},
    "MOTION.spectra": {short: "Spectra", title: "In-structure response spectrum (Nigam–Jennings)", ref: ["11.3", "113-response-spectra-nigam-jennings"],
      tex: String.raw`\ddot{x} + 2\zeta\omega\dot{x} + \omega^{2}x = -a(t),\qquad S_a = \max_t\left|2\zeta\omega\dot{x} + \omega^{2}x\right|`,
      detail: (d) => `node ${d.node} ${d.dof} · ${n(d.nfreq)} oscillator frequencies × ${n(d.ndamp)} damping value(s)`},
    "STRESS.elements": {short: "Element forces", title: "Element forces and stresses from the nodal motions", ref: ["11.5", "115-stresses-and-forces-stress"],
      tex: String.raw`\sigma_c(\omega) = \mathbf{S}_c\,\mathbf{u}_e(\omega),\qquad \sigma_c(t) = \mathcal{F}^{-1}\left\{\sigma_c(\omega)\,U_g(\omega)\right\}`,
      detail: (d) => `${d.type}: elements ${n(Number(d.done) + 1)}–${n(Number(d.done) + Number(d.mb))} of ${n(d.total)}`},
    "RELDISP.relative": {short: "Relative", title: "Relative displacement", ref: ["11.4", "114-relative-displacements-reldisp"],
      tex: String.raw`D(f) = \left(H_{node}(f) - H_{ref}(f)\right)U_g(f),\qquad d(t) = \mathcal{F}^{-1}\{D\}`,
      detail: (d) => `node ${d.node} ${d.dof} (${d.i}/${d.n}) relative to ${d.freefield ? "the free field" : "the reference node"}`},
    "SOIL.iterate": {short: "Equivalent linear", title: "Equivalent-linear site response (SHAKE)", ref: ["12.2", "122-strains-and-the-equivalent-linear-iteration"],
      tex: String.raw`\gamma_{\mathrm{eff}} = R_\gamma\,\gamma_{\max},\qquad G = G_{\max}\,\frac{G}{G_{\max}}(\gamma_{\mathrm{eff}}),\quad \beta = D(\gamma_{\mathrm{eff}})`,
      detail: (d) => d.vertical ? `${n(d.nsub)} sublayers · vertical input: no iterations` : `${n(d.nsub)} sublayers · ${n(d.niter)} iterations · R_γ = ${S.fmt(d.ratio, 3)} · N = ${n(d.nfft)}`},
    "EQUAKE.fit": {short: "Spectrum matching", title: "Spectrum-compatible motion: frequency-domain matching", ref: ["21", "21-spectrum-compatible-motions-equake"],
      tex: String.raw`|A_k|^{(i+1)} = |A_k|^{(i)}\,\frac{S_a^{\,\mathrm{target}}(f_k)}{S_a^{(i)}(f_k)}`,
      detail: (d) => `spectrum ${d.k}/${d.n} · ${S.fmt(d.dur, 4)} s at Δt = ${S.fmt(d.dt, 4)} s · damping ${S.fmt(100 * d.zeta, 3)} %`},
  };
  /** The steps of each module, in order (the strip). */
  const MODULE_STEPS = {
    SITE: ["SITE.modes", "SITE.field"], POINT: ["POINT.green"], HOUSE: ["HOUSE.elements", "HOUSE.assembled"],
    FORCE: ["FORCE.loads"], ANALYS: ["ANALYS.impedance", "ANALYS.dynamic", "ANALYS.factor", "ANALYS.solve"],
    MOTION: ["MOTION.fft", "MOTION.interp", "MOTION.convolve", "MOTION.spectra"], STRESS: ["STRESS.elements"],
    RELDISP: ["RELDISP.relative"], SOIL: ["SOIL.iterate"], EQUAKE: ["EQUAKE.fit"],
  };
  A.STEPS = STEPS;
  A.MODULE_STEPS = MODULE_STEPS;
  const MATH_KEY = "sassi-edu.activityMath";
  const mathOn = () => { try { return localStorage.getItem(MATH_KEY) !== "0"; } catch (e) { return true; } };
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
      mathStrip: el("div", {class: "act-math-strip"}),
      mathTitle: el("span", {class: "act-math-t"}),
      mathRef: el("button", {class: "act-math-ref", type: "button", title: "open this step in the Theory Manual"}),
      mathEq: el("div", {class: "act-math-eq"}),
      mathDetail: el("div", {class: "act-math-d"}),
      mathToggle: el("button", {class: "act-x act-math-toggle", type: "button", title: "show / hide the equations of the steps", text: "∑"}),
      counts: el("span", {class: "act-counts"}),
      hist: el("button", {class: "btn small", text: "Command History"}),
    };
    ui.mod = el("div", {class: "act-mod", hidden: true}, el("div", {class: "act-mod-bar"}, ui.modFill), ui.modText);
    ui.math = el("div", {class: "act-math", hidden: true}, ui.mathStrip,
      el("div", {class: "act-math-head"}, ui.mathTitle, ui.mathRef), ui.mathEq, ui.mathDetail);
    box = el("div", {id: "activity", class: "activity", role: "status", "aria-live": "polite", hidden: true},
      el("div", {class: "act-head"}, ui.spin, ui.title, ui.time, ui.mathToggle, ui.close),
      el("div", {class: "act-bar"}, ui.fill), ui.line,
      el("div", {class: "act-nowrow"}, el("span", {class: "act-k", text: "Now"}), ui.now), ui.cmd, ui.mod,
      ui.math, ui.chain,
      el("div", {class: "act-foot"}, ui.counts, ui.hist));
    ui.close.addEventListener("click", () => hide(true));
    ui.mathToggle.addEventListener("click", () => {
      const on = !mathOn();
      try { localStorage.setItem(MATH_KEY, on ? "1" : "0"); } catch (e) { /* private mode */ }
      ui.mathToggle.classList.toggle("on", on);
      if (!on) clearMath(); else if (lastStep) showStep(lastStep.module, lastStep.key, lastStep.data, true);
    });
    ui.mathRef.addEventListener("click", () => {
      const st = lastStep && STEPS[lastStep.key];
      if (st && S.D && S.D.helpTab) S.D.helpTab({doc: TH, anchor: st.ref[1]});
    });
    ui.mathToggle.classList.toggle("on", mathOn());
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
    run = {kind, source, t0: Date.now(), line: 0, total: 0, warnings: 0, errors: 0, modules: new Map(), models: [], current: null, dismissed: false};
    show();
    clearMath();
    lastStep = null;
    if (ui) delete ui.mathStrip.dataset.module;
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

  // ---------------------------------------------------------------- the computation step (the math)
  let lastStep = null, shownKey = "", shownTex = "", stepTimer = null;
  function clearMath() {
    if (!ui) return;
    ui.math.hidden = true;
    box.classList.remove("with-math");
    shownKey = shownTex = "";
  }
  /** A module reports the step it enters (at most one repaint per 120 ms: ANALYS reports four steps per
   *  frequency; the newest is shown). */
  function showStep(module, key, data, now) {
    lastStep = {module: String(module || "").toUpperCase(), key: String(key || ""), data: data || {}};
    if (!mathOn() || !STEPS[lastStep.key]) return;
    if (now) { paintStep(); return; }
    if (!stepTimer) stepTimer = setTimeout(() => { stepTimer = null; paintStep(); }, 120);
  }
  A.showStep = showStep;
  function paintStep() {
    if (!run || !lastStep || !ui) return;
    const {module, key, data} = lastStep;
    const st = STEPS[key];
    if (!st) return;
    ui.math.hidden = false;
    box.classList.add("with-math");
    // the strip: the module's steps, the current one lit
    const order = MODULE_STEPS[module] || [key];
    if (ui.mathStrip.dataset.module !== module) {
      ui.mathStrip.innerHTML = "";
      ui.mathStrip.dataset.module = module;
      for (const k of order) ui.mathStrip.appendChild(el("span", {class: "act-math-s", "data-key": k, text: (STEPS[k] || {}).short || k}));
    }
    const idx = order.indexOf(key);
    ui.mathStrip.querySelectorAll(".act-math-s").forEach((c, i) => {
      c.classList.toggle("on", i === idx);
      c.classList.toggle("past", idx >= 0 && i < idx);
    });
    if (shownKey !== key) {
      shownKey = key;
      ui.mathTitle.textContent = `${module}: ${st.title}`;
      ui.mathRef.textContent = `Theory §${st.ref[0]}`;
    }
    const tex = typeof st.tex === "function" ? st.tex(data) : st.tex;
    if (tex !== shownTex) {
      shownTex = tex;
      if (window.katex) {
        try { window.katex.render(tex, ui.mathEq, {displayMode: true, throwOnError: false, strict: "ignore"}); }
        catch (e) { ui.mathEq.textContent = tex; }
      } else ui.mathEq.textContent = tex;
    }
    try { ui.mathDetail.textContent = st.detail(data); } catch (e) { ui.mathDetail.textContent = ""; }
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
    clearMath();
    lastStep = null;
    paint();
    const ended = run;
    run = null;
    // a clean run fades out on its own; with errors the panel stays until it is closed
    if (!bad) hideTimer = setTimeout(() => { if (!run && ended) hide(); }, 4500);
    // the run summary window (static/summary.js) after a run that ran modules
    if (S.onRunFinished && ended.modules.size) {
      try {
        S.onRunFinished({kind: ended.kind, source: ended.source, modules: [...ended.modules.keys()], models: ended.models.slice(),
          warnings: ended.warnings, errors: ended.errors, seconds: (Date.now() - ended.t0) / 1000});
      } catch (e) { console.error(e); }
    }
  }

  // ---------------------------------------------------------------- the session's events
  A.onEvent = function (ev) {
    if (!ev) return;
    if (ev.type === "progress") {
      if (ev.kind === "step") {
        if (run) showStep(ev.module, ev.key, ev.data);
        return;
      }
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
      if (m && MODULES[m[1]] !== undefined) {
        const mm = /^RUN[A-Z]+: model (\S+)/.exec(text);          // the models the run ran, in order
        if (mm && !run.models.includes(mm[1])) run.models.push(mm[1]);
        moduleStart(m[1]);
        paint();
        return;
      }
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
        if (job.step && job.step.key) showStep(job.module, job.step.key, job.step.data);
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
