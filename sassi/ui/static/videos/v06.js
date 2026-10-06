/* Explainer video, lesson 6 (newcomer edition): Seismic input -- the two inputs every SSI run starts from: a
 * record that matches the design spectrum (EQUAKE) and the soil as that earthquake shakes it (SOIL), then the
 * soil cases and their envelope.
 * Numbers: sassi/ui/lessons/06_seismic_input.md; curves: d06.js (the lesson's run, python -m sassi.ui.video_data);
 * the library Sand curves are SV.K.PH.SAND (figures.js).
 * Soil columns move with the exact 1-D wave solution of the lesson's site (PH.columnWaves / columnU): 10 m sand,
 * 12 m clay, rock (Vs 1200 m/s, 1 %), with the strain-compatible sublayers of FILE88 (best estimate, lower and
 * upper bound) or the low-strain table of the lesson.  The drawn motion is the steady state at the frequency
 * named on screen, scaled so that the surface acceleration amplitude is that case's surface PGA (0.544 / 0.605 /
 * 0.690 g, the SOIL listings): u0 = PGA g / (2 pi f)^2; one displacement factor per scene, stated on screen. */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["06"], PH = K.PH;
  const TAU = 2 * Math.PI, G = 9.81;
  const CLAY = "#8d7a63";

  // ---------------------------------------------------------------- data
  const TGT = D["data/rg160h_030g.rsi"].cols;          // [f, SA (g)]: the RG 1.60 target, 27 points
  const RSO = D["ex04/rg160h_eq.rso"].cols;            // [f, SA (g)]: 5 % spectrum of the EQUAKE record (rock)
  const REC = D["ex04/rg160h_eq.acc"];                 // the EQUAKE record (g), 20 s, peak 0.324 g
  const SURF = D["ex04/ACC001.TH"];                    // SOIL surface motion (g), peak 0.605 g
  const RS1 = D["ex04/RS001_01.RS"];                   // 5 % surface spectrum, best estimate
  const LB = D["ex04lb/RS001_01.RS"], UB = D["ex04ub/RS001_01.RS"], ENV = D["ex04/surface_envelope.rs"];
  const F88 = D["ex04/FILE88"].cols;                   // FILE88: layer thick gamma_eff_pct G Vs beta_s Vp beta_p

  /** The largest y and where it is. */
  const argmax = (xs, ys) => { let i = 0; ys.forEach((y, j) => { if (y > ys[i]) i = j; }); return {x: xs[i], y: ys[i], i}; };
  const PK_TGT = argmax(TGT[0], TGT[1]);               // 2.5 Hz: the corner of the RG 1.60 spectrum
  const PK_ROCK = argmax(RSO[0], RSO[1]);              // 2.57 Hz, 1.005 g
  const PK_BE = argmax(RS1.f, RS1.sa);                 // 1.86 Hz, 2.24 g
  const PK_UB = argmax(UB.f, UB.sa);                   // 2.57 Hz, 3.14 g

  /** Log-log interpolation of (X, Y) at x. */
  function llInterp(X, Y, x) {
    if (x <= X[0]) return Y[0];
    if (x >= X[X.length - 1]) return Y[Y.length - 1];
    let i = 1;
    while (X[i] < x) i++;
    const t = (Math.log(x) - Math.log(X[i - 1])) / (Math.log(X[i]) - Math.log(X[i - 1]));
    return Math.exp(Math.log(Y[i - 1]) + t * (Math.log(Y[i]) - Math.log(Y[i - 1])));
  }

  /** How EQUAKE tunes noise (illustration, computed here): Gaussian noise under an envelope (dt 0.01 s, 20 s),
   *  then frequency-domain passes that multiply its Fourier amplitudes by target / its 5 % spectrum, phases
   *  kept (Levy-Wilkinson; EQUAKE then refines in the time domain).  Computed once. */
  let TUNE = null;
  function tuning() {
    if (TUNE) return TUNE;
    const dt = 0.01, n = 2000, nfft = 4096, df = 1 / (nfft * dt), NP = 4;
    let st = 11975;
    const rnd = () => { st = (Math.imul(1664525, st) + 1013904223) >>> 0; return st / 4294967296; };
    const env = (t) => (t < 2 ? (t / 2) ** 2 : t < 12 ? 1 : Math.exp(-0.35 * (t - 12)));
    let a = new Float64Array(n), pk = 0;
    for (let i = 0; i < n; i++) { a[i] = Math.sqrt(-2 * Math.log(1 - rnd())) * Math.cos(TAU * rnd()) * env(i * dt); pk = Math.max(pk, Math.abs(a[i])); }
    for (let i = 0; i < n; i++) a[i] *= 0.26 / pk;
    const fs = PH.logspace(0.1, 50, 48), passes = [];
    for (let p = 0; p <= NP; p++) {
      const sa = PH.spectrum(a, dt, fs, 0.05);
      passes.push({a: Array.from(a), sa});
      if (p === NP) break;
      const ratio = fs.map((f, i) => llInterp(TGT[0], TGT[1], f) / sa[i]);
      const F = PH.rfft(a, nfft);
      for (let k = 0; k <= nfft / 2; k++) {
        const f = k * df, r = f < 0.1 ? 0 : llInterp(fs, ratio, Math.min(50, f));
        F.re[k] *= r; F.im[k] *= r;
      }
      const b = PH.irfft(F.re, F.im, nfft);
      a = new Float64Array(n);
      for (let i = 0; i < n; i++) a[i] = b[i];
    }
    TUNE = {t: Array.from({length: n}, (_, i) => i * dt), fs, passes, NP};
    return TUNE;
  }

  // ---------------------------------------------------------------- the site and its 1-D wave solution
  const ROCK = {vs: 1200, w: 22, beta: 0.01};          // L 3, the half-space (linear)
  /** The soil report (low strain): L 1 sand and L 2 clay of the lesson, 1 % damping. */
  const LOWCOL = {g: G, layers: [{h: 10, vs: 200, w: 19, beta: 0.01}, {h: 12, vs: 300, w: 18.5, beta: 0.01}], hs: ROCK};
  /** A strain-compatible column from FILE88 (22 sublayers of 1 m; unit weight from G and Vs). */
  const f88col = (F) => ({g: G, hs: ROCK, layers: F[0].map((_, i) => ({h: F[1][i], vs: F[4][i], w: (F[3][i] * G) / (F[4][i] * F[4][i]), beta: F[5][i]}))});
  const COL = {be: f88col(F88), lb: f88col(D["ex04lb/FILE88"].cols), ub: f88col(D["ex04ub/FILE88"].cols)};
  /** Site frequencies (peaks of the SOIL surface / base amplification) and surface PGAs: the lesson's table. */
  const CASE = {lb: {f: 1.37, pga: 0.544}, be: {f: 1.83, pga: SURF.peak}, ub: {f: 2.56, pga: 0.690}};
  const F_LOW = 3.10;                                  // the low-strain amplification peak (lesson: 68 at 3.10 Hz)
  const SLOW = 3;                                      // every soil motion is slowed down 3 x
  /** Displacement (px) at depth d px of a column drawn at pxm px/m, at scene time t: the steady state at f,
   *  surface amplitude pga g / (2 pi f)^2 times the scene factor X, slowed down SLOW x. */
  function motion(col, f, pga, pxm, X, zmax) {
    const wv = PH.columnWaves(col, f), U = [];
    for (let i = 0; i <= Math.ceil(zmax * 10); i++) U.push(PH.columnU(wv, col, i / 10));
    const amp = X * pxm * (pga * G) / ((TAU * f) * (TAU * f));
    return (d, t) => {
      const u = U[Math.max(0, Math.min(U.length - 1, Math.round((d / pxm) * 10)))], ph = (TAU * f * t) / SLOW;
      return amp * (u[0] * Math.cos(ph) - u[1] * Math.sin(ph));       // Re[U e^{i w t}]
    };
  }

  /** The lesson's column: 10 m sand, 12 m clay, rock; pxm px per metre; names: labels on the right. */
  function column(g, o) {
    const n = o.names;
    return K.soil(g, {x: o.x, y: o.y, w: o.w, hidden: o.hidden,
      layers: [{h: 10 * o.pxm, kind: "sand", name: n ? "Sand" : "", vs: n ? 200 : null},
        {h: 12 * o.pxm, kind: "gravel", color: CLAY, name: n ? "Clay" : "", vs: n ? 300 : null}],
      hs: {h: o.hs || 100, kind: "rock", name: n ? "Rock" : "", vs: n ? 1200 : null},
      labels: n ? "right" : false, labelSize: 32});
  }
  /** A small muted note (the scale of a picture). */
  const note = (g, x, y, str, o) => K.text(g, x, y, str, Object.assign({cls: "t-label", size: 24, color: "var(--muted)", hidden: true}, o));
  /** A centred pill of width w. */
  const pillC = (root, html, cx, y, w, color, size) =>
    K.h("div", {cls: "sv-pill", html, x: cx - w / 2, y, w, align: "center", color, size: size || 34, in: "pop"}, root);

  const scenes = [];

  // ---------------------------------------------------------------- 0. hook
  scenes.push({
    id: "hook", title: "Two inputs first",
    build(s) {
      const g = K.g(s.svg);
      s.l1 = K.heading(s.root, "Design spectrum", {x: 230, y: 120, size: "h3"});
      s.p = K.plot(s.svg, {x: 250, y: 210, w: 600, h: 400, xr: [0.1, 100], xlog: true, yr: [0, 1.2], xticks: [0.1, 1, 10, 100], yticks: [0, 0.4, 0.8, 1.2],
        xlabel: "frequency (Hz)", ylabel: "SA (g)", ylabelOffset: 80});
      s.tgt = s.p.line(TGT[0], TGT[1], {color: C.ref, width: 6, draw: true});
      s.l2 = K.heading(s.root, "Soil report", {x: 1150, y: 120, size: "h3"});
      const pxm = 17, top = 200;
      s.soil = column(g, {x: 1150, y: top, w: 300, pxm, hs: 90, names: true, hidden: true});
      s.dims = K.g(g, {hidden: true});
      K.dim(s.dims, 1120, top, 1120, top + 10 * pxm, "10 m", {side: "left"});
      K.dim(s.dims, 1120, top + 10 * pxm, 1120, top + 22 * pxm, "12 m", {side: "left"});
      s.mv = motion(COL.be, CASE.be.f, CASE.be.pga, pxm, 25, 22 + 90 / pxm);
      s.note = note(g, 1150, 702, `this quake: 1-D waves at ${CASE.be.f} Hz\ndisplacements × 25 · slowed ${SLOW} ×`);
      s.q = K.h("div", {x: 0, y: 800, w: 1920, align: "center", html: "Straight into SASSI?", size: 76, in: "pop", style: {fontWeight: "800"}}, s.root);
      s.ar1 = K.arrow(g, 550, 712, 550, 805, {color: C.wave, width: 6, head: 20, hidden: true});
      s.ar2 = K.arrow(g, 1300, 760, 1300, 805, {color: C.ssi, width: 6, head: 20, hidden: true});
      s.t1 = pillC(s.root, "1 · turn it into a record", 550, 822, 520, "var(--wave)");
      s.t2 = pillC(s.root, "2 · soften it for that quake", 1300, 822, 560, "var(--ssi)");
      s.A = 0;
    },
    tick(s, t) {
      s.soil.shear((d) => s.A * s.mv(d, t));
    },
    beats: [
      {say: "[a]You have a design spectrum, [b]and a soil report with the stiffness of every layer.",
        a(k) { k.show([k.s.l1, k.s.p.g]); k.draw(k.s.tgt, 1000); },
        b(k) { k.show([k.s.l2, k.s.soil.g, k.s.dims]); }},
      {say: "[q]Can you feed them straight into SASSI?", q(k) { k.show(k.s.q); }},
      {say: "[n]Not quite: [a]the spectrum must become a record, [b]and the soil must match how hard that quake shakes it.",
        n(k) { k.hide(k.s.q); k.show([k.s.ar1, k.s.ar2]); }, a(k) { k.show(k.s.t1); },
        b(k) { k.show([k.s.t2, k.s.note]); k.tween(k.s, {A: 1}, 800); }},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 6",
    build(s) {
      s.t = K.titleCard(s.root, {n: 6, part: "Design applications", title: "Seismic input", sub: "A record that matches your design spectrum, and the soil as that earthquake shakes it."});
    },
    beats: [
      {say: "Lesson six: seismic input. The record, and the soil, that every SSI run starts from.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 100}); k.show(t.title, {delay: 200}); k.show(t.sub, {delay: 350}); k.show(t.rule, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. a matching record (EQUAKE)
  scenes.push({
    id: "record", title: "Input 1: a matching record",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Input 1: a matching record", {x: 120, y: 70, size: "h2"});
      s.tp = K.plot(s.svg, {x: 230, y: 240, w: 770, h: 220, xr: [0, 20], yr: [-0.4, 0.4], xticks: [0, 5, 10, 15, 20], yticks: [-0.3, 0, 0.3],
        xlabel: "time (s)", ylabel: "acc. (g)", ylabelOffset: 96});
      s.recLab = K.text(g, 230, 212, "the shaking, second by second", {cls: "t-label", size: 32, color: "var(--wave)", hidden: true});
      s.rec = s.tp.line(REC.t, REC.v, {color: C.wave, width: 3, draw: true});
      // the peak of the record, from the data
      let ip = 0;
      REC.v.forEach((v, i) => { if (Math.abs(v) > Math.abs(REC.v[ip])) ip = i; });
      s.pk = K.g(s.tp.g, {hidden: true});
      const px = s.tp.X(REC.t[ip]), py = s.tp.Y(REC.v[ip]);
      K.circle(s.pk, px, py, 10, {fill: C.wave, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.pk, px + 20, py + (REC.v[ip] > 0 ? 4 : 28), `${REC.peak.toFixed(3)} g`, {cls: "t-label", size: 34, color: "var(--wave)", weight: 800});
      // EQUAKE: noise, tuned pass by pass (illustration computed here), then EQUAKE's record
      const T = tuning();
      s.NP = T.NP;
      s.nTr = s.tp.line(T.t, T.passes[0].a, {color: C.muted, width: 2.5, hidden: true});
      s.dT = T.passes.map((p) => K.d(s.tp.pts(T.t, p.a)));
      s.eqk = K.h("div", {cls: "sv-pill", html: "EQUAKE", x: 230, y: 572, color: "var(--wave)", size: 38, in: "pop", style: {fontFamily: "var(--mono)"}}, s.root);
      s.passLab = K.text(g, 452, 616, "", {cls: "t-label", size: 32, color: "var(--ink2)", hidden: true});
      s.illus = note(g, 452, 656, "passes: an illustration computed here");
      s.card = K.card(s.root, {x: 230, y: 700, w: 770, kind: "check", title: "Spectrum-compatible", size: 34, body: "Its spectrum matches your design spectrum."});
      // the spectra
      s.p = K.plot(s.svg, {x: 1200, y: 210, w: 600, h: 430, xr: [0.1, 100], xlog: true, yr: [0, 1.2], xticks: [0.1, 1, 10, 100], yticks: [0, 0.4, 0.8, 1.2],
        xlabel: "frequency (Hz)", ylabel: "SA (g), 5 % damping", ylabelOffset: 80});
      s.tgt = s.p.line(TGT[0], TGT[1], {color: C.ref, width: 7, draw: true});
      s.nSp = s.p.line(T.fs, T.passes[0].sa, {color: C.muted, width: 4, hidden: true});
      s.dS = T.passes.map((p) => K.d(s.p.pts(T.fs, p.sa)));
      s.rs = s.p.line(RSO[0], RSO[1], {color: C.wave, width: 4, draw: true});
      s.leg0 = s.p.legend([{label: "design spectrum", color: C.ref}], {x: 1230, y: 770, size: 26, hidden: true});
      s.legN = s.p.legend([{label: "noise, tuned", color: C.muted}], {x: 1230, y: 810, size: 26, hidden: true});
      s.leg1 = s.p.legend([{label: "EQUAKE record", color: C.wave}], {x: 1230, y: 810, size: 26, hidden: true});
      s.P = 0;
    },
    tick(s) {
      if (s.final) return;
      const p = Math.min(s.NP, Math.floor(s.P + 1e-6));
      if (p !== s.pShown) {
        s.pShown = p;
        s.nTr.setAttribute("d", s.dT[p]);
        s.nSp.setAttribute("d", s.dS[p]);
        s.passLab.textContent = p === 0 ? "start: random noise" : `pass ${p} of ${s.NP}: tuned`;
      }
    },
    beats: [
      {say: "[a]SASSI can't shake a model with a spectrum. [b]It needs the shaking itself, second by second.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.tp.g); }, b(k) { k.show(k.s.recLab); k.draw(k.s.rec, 1000); }},
      {say: "[e]So where does such a record come from? [q]A SASSI module called EQUAKE makes one.",
        e(k) { k.hide([k.s.rec, k.s.recLab]); }, q(k) { k.show(k.s.eqk); }},
      {say: "[a]It starts from random noise, like static on a radio, [b]then tunes it, note by note.",
        a(k) { const s = k.s; k.show([s.nTr, s.passLab, s.illus, s.p.g, s.nSp, s.leg0, s.legN]); k.draw(s.tgt, 800); },
        b(k) { k.tween(k.s, {P: k.s.NP}, 2600, {ease: "linear"}); }},
      {say: "[s]It keeps tuning until the record's spectrum sits right on top of your design spectrum.",
        s(k) {
          const s = k.s;
          s.final = true;
          k.hide([s.nTr, s.nSp, s.illus, s.legN]); k.text(s.passLab, `EQUAKE's record, ${Math.round(REC.t[REC.t.length - 1])} s`);
          k.draw(s.rec, 1000, {from: 0}); k.draw(s.rs, 1000); k.show(s.leg1);
        }},
      {say: "[p]The result: a 20 second record, peaking at 0.324 g. [c]Engineers call it spectrum-compatible.",
        p(k) { k.show(k.s.pk); }, c(k) { k.show(k.s.card); }},
    ],
  });

  // ---------------------------------------------------------------- 3. shaking softens soil
  // the largest effective strain of the sand in this earthquake (FILE88: 0.51 % at 9.5 m)
  const SAND_ROWS = F88[0].map((_, i) => i).filter((i) => i < 10);
  const GAM_MAX = Math.max(...SAND_ROWS.map((i) => F88[2][i]));
  scenes.push({
    id: "soften", title: "Input 2: shaking softens soil",
    build(s) {
      K.defs();
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Input 2: shaking softens soil", {x: 120, y: 70, size: "h2"});
      // a piece of soil in cyclic shear: the shear angle is the strain on the curve, drawn x 30
      s.el = K.g(g, {hidden: true});
      s.elPoly = K.path(s.el, "", {fill: "var(--sand)", "fill-opacity": 0.85, stroke: "rgba(255,240,215,.6)", "stroke-width": 3});
      s.elTex = K.path(s.el, "", {fill: "url(#sv-pat-sand)"});
      s.elTop = K.arrow(s.el, 0, 0, 1, 0, {color: C.ssi, width: 6, head: 20});
      s.elBot = K.arrow(s.el, 0, 0, 1, 0, {color: C.ssi, width: 6, head: 20});
      K.text(s.el, 340, 700, "a piece of soil", {cls: "t-label", anchor: "middle", size: 32});
      s.gam = K.text(s.el, 340, 750, "", {cls: "t-label", anchor: "middle", size: 32, color: "var(--wave)", weight: 750});
      K.text(s.el, 340, 796, "strain drawn × 30", {cls: "t-label", anchor: "middle", size: 24, color: "var(--muted)"});
      // the library Sand curves (SHAKE91): stiffness and damping against strain
      const S = PH.SAND;
      const xfmt = (v) => (v >= 0.01 ? String(v) : v === 0.001 ? "0.001" : "0.0001");
      s.pg = K.plot(s.svg, {x: 860, y: 210, w: 800, h: 460, xr: [1e-4, 10], xlog: true, yr: [0, 1], xticks: [1e-4, 1e-3, 0.01, 0.1, 1, 10], yticks: [0, 0.5, 1],
        xfmt, xlabel: "shear strain (%)", ylabel: "stiffness (G / Gmax)", ylabelOffset: 80});
      s.ax2 = K.g(s.pg.g, {hidden: true});
      K.line(s.ax2, 1660, 210, 1660, 670, {stroke: C.dash, "stroke-width": 2});
      [0, 10, 20, 30].forEach((d) => K.text(s.ax2, 1674, s.pg.Y(d / 30) + 8, String(d), {cls: "tick", color: "var(--dash)"}));
      const al = K.text(s.ax2, 1750, 440, "damping (%)", {cls: "axlab", anchor: "middle", color: "var(--dash)"});
      al.setAttribute("transform", "rotate(90 1750 440)");
      s.gl = s.pg.line(S.g, S.G, {color: C.spring, width: 6, draw: true});
      s.dl = s.pg.line(S.g, S.D.map((d) => d / 30), {color: C.dash, width: 6, draw: true});
      s.rep = s.pg.text(S.g[0], S.G[0], "soil report", {cls: "t-label", size: 30, color: "var(--wave)", dx: 16, dy: 44, hidden: true});
      s.legG = K.g(s.pg.g, {hidden: true});
      K.line(s.legG, 900, 392, 946, 392, {stroke: C.spring, "stroke-width": 6, "stroke-linecap": "round"});
      K.text(s.legG, 962, 402, "stiffness falls", {cls: "t-label", size: 32, color: "var(--spring)"});
      s.legD = K.g(s.pg.g, {hidden: true});
      K.line(s.legD, 900, 446, 946, 446, {stroke: C.dash, "stroke-width": 6, "stroke-linecap": "round"});
      K.text(s.legD, 962, 456, "damping rises", {cls: "t-label", size: 32, color: "var(--dash)"});
      s.gdot = K.circle(s.pg.g, 0, 0, 12, {fill: C.ink, stroke: "#0a111d", "stroke-width": 3});
      s.ddot = K.circle(s.pg.g, 0, 0, 12, {fill: C.ink, stroke: "#0a111d", "stroke-width": 3});
      s.vl = K.line(s.pg.g, 0, 210, 0, 670, {stroke: "rgba(238,243,249,.35)", "stroke-width": 2, "stroke-dasharray": "6 6"});
      [s.gdot, s.ddot, s.vl].forEach((e) => e.classList.add("sv-in", "sv-fade"));
      s.lg = Math.log10(S.g[0]); s.A = 0;
    },
    tick(s, t) {
      const S = PH.SAND, lgs = S.g.map(Math.log10);
      const G0 = s.lg <= lgs[0] ? S.G[0] : PH.interpLin(lgs, S.G, s.lg), Dm = s.lg <= lgs[0] ? S.D[0] : PH.interpLin(lgs, S.D, s.lg);
      const xg = s.pg.X(Math.pow(10, s.lg));
      s.gdot.setAttribute("cx", xg); s.gdot.setAttribute("cy", s.pg.Y(G0));
      s.ddot.setAttribute("cx", xg); s.ddot.setAttribute("cy", s.pg.Y(Dm / 30));
      s.vl.setAttribute("x1", xg); s.vl.setAttribute("x2", xg);
      // the piece of soil: shear angle = the strain (rad) x 30, cycling
      const x0 = 210, y0 = 350, w = 260, h = 260;
      const ph = TAU * t / 1.2, off = s.A * (Math.pow(10, s.lg) / 100) * h * 30 * Math.sin(ph);
      const d = `M${x0 + off},${y0}L${x0 + w + off},${y0}L${x0 + w},${y0 + h}L${x0},${y0 + h}Z`;
      s.elPoly.setAttribute("d", d); s.elTex.setAttribute("d", d);
      const dir = Math.cos(ph) >= 0 ? 1 : -1;
      s.elTop.update(x0 + w / 2 + off - 50 * dir, y0 - 28, x0 + w / 2 + off + 50 * dir, y0 - 28);
      s.elBot.update(x0 + w / 2 + 50 * dir, y0 + h + 28, x0 + w / 2 - 50 * dir, y0 + h + 28);
      const gp = Math.pow(10, s.lg), txt = `strain ${gp < 0.01 ? gp.toFixed(4) : gp < 0.1 ? gp.toFixed(3) : gp.toFixed(2)} %`;
      if (s.gamTxt !== txt) {
        s.gamTxt = txt;
        s.gam.textContent = txt;
        s.gam.style.fill = gp > 0.01 ? "var(--ssi)" : "var(--wave)";
      }
    },
    beats: [
      {say: "[a]Now the soil. Its stiffness in the soil report is measured with very gentle shaking.",
        go(k) { k.show(k.s.head); },
        a(k) { k.show(k.s.el); k.tween(k.s, {A: 1}, 800); k.show(k.s.pg.g, {delay: 200}); k.draw(k.s.gl, 1000, {delay: 300}); k.show([k.s.gdot, k.s.vl, k.s.rep], {delay: 400}); }},
      {say: "[b]But soil is a pile of grains. [c]Shake it hard, and the grains slip: it turns softer.",
        b(k) { k.pulse(k.s.el, {amp: 0.05}); }, c(k) { k.hide(k.s.rep); k.tween(k.s, {lg: Math.log10(GAM_MAX)}, 2600); }},
      {say: "[d]The slipping also burns energy, so the soil damps the shaking more.",
        d(k) { k.show([k.s.ax2, k.s.ddot]); k.draw(k.s.dl, 1000); }},
      {say: "Engineers plot both against strain, how far the soil is pushed out of shape: [g]stiffness falls, [h]damping rises.",
        g(k) { k.show(k.s.legG); k.pulse(k.s.gl, {amp: 0.03}); }, h(k) { k.show(k.s.legD); k.pulse(k.s.dl, {amp: 0.03}); }},
    ],
  });

  // ---------------------------------------------------------------- 4. guess, check, repeat (SOIL)
  scenes.push({
    id: "soil", title: "Guess, check, repeat",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Guess, check, repeat: SOIL", {x: 120, y: 70, size: "h2"});
      const pxm = 24, top = 250, X = 40;
      s.soil = column(g, {x: 150, y: top, w: 280, pxm, hs: 100, hidden: true});
      s.labs = K.g(g, {hidden: true});
      [["sand", 5], ["clay", 16], ["rock", 22 + 50 / pxm]].forEach(([t, z]) => K.text(s.labs, 500, top + z * pxm + 10, t, {cls: "t-label", size: 32}));
      s.inp = K.g(g, {hidden: true});
      K.arrow(s.inp, 290, 990, 290, top + 22 * pxm + 100 + 12, {color: C.wave, width: 6, head: 20});
      K.text(s.inp, 322, 962, "the record, on rock", {cls: "t-label", size: 30, color: "var(--wave)"});
      // the column's motion: soil-report stiffness (round 1), then the converged soil
      s.mLow = motion(LOWCOL, F_LOW, CASE.be.pga, pxm, X, 22 + 100 / pxm);
      s.mSc = motion(COL.be, CASE.be.f, CASE.be.pga, pxm, X, 22 + 100 / pxm);
      s.nf = note(g, 150, 196, "", {color: "var(--ink2)", size: 26});
      s.nx = note(g, 150, 230, `1-D waves · displacements × ${X} · slowed ${SLOW} ×`);
      // chicken and egg
      s.ce = K.h("div", {x: 950, y: 500, w: 700, align: "center", html: "softness ⇄ shaking", size: 56, in: "pop", style: {fontWeight: "800", color: "var(--ink)"}}, s.root);
      // the loop
      const cx = 1300, cy = 560, R = 260;
      s.cx = cx; s.cy = cy; s.R = R;
      const P = (a) => [cx + R * Math.cos((a * Math.PI) / 180), cy + R * Math.sin((a * Math.PI) / 180)];
      s.ring = K.circle(g, cx, cy, R, {fill: "none", stroke: "#2a3b52", "stroke-width": 3, "stroke-dasharray": "10 10", hidden: true});
      const node = (a, html, color) => {
        const [x, y] = P(a);
        return K.h("div", {x: x - 180, y: y - 36, w: 360, align: "center", html, in: "pop", size: 32,
          style: {background: "rgba(15,26,43,.96)", border: `3px solid ${color}`, borderRadius: "16px", padding: "12px 14px", fontWeight: "650"}}, s.root);
      };
      s.n = [node(-90, "Shake the column", "var(--wave)"), node(30, "Measure the strain", "var(--ssi)"), node(150, "Soften each layer", "var(--spring)")];
      s.arrs = [[-50, -5, -40], [65, 115, -48], [185, 230, -40]].map(([a0, a1, bend]) => {
        const [x1, y1] = P(a0), [x2, y2] = P(a1);
        return K.curveArrow(g, x1, y1, x2, y2, {color: C.muted, width: 5, bend, head: 20, hidden: true});
      });
      s.tag = K.h("div", {x: cx - 150, y: cy - 34, w: 300, align: "center", html: "SOIL", size: 52, in: "pop", style: {fontFamily: "var(--mono)", fontWeight: "700", color: "var(--ssi)"}}, s.root);
      s.tok = K.circle(g, cx, cy - R, 14, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3, hidden: true});
      s.iter = K.h("div", {x: cx - 150, y: cy - 28, w: 300, align: "center", html: "round <b>1</b>", size: 40, in: "fade"}, s.root);
      s.eqlin = pillC(s.root, "the equivalent-linear method", cx, 880, 560, "var(--ink2)", 34);
      // the result: the strain-compatible profile, on the column's depth scale
      s.pv = K.plot(s.svg, {x: 760, y: top, w: 560, h: 22 * pxm, xr: [0, 350], yr: [22, 0], xticks: [0, 100, 200, 300], yticks: [0, 5, 10, 15, 20],
        xlabel: "shear-wave speed Vs (m/s)", ylabel: "depth (m)", ylabelOffset: 66});
      const PZ = [], PV = [];
      let z = 0;
      F88[0].forEach((_, i) => { PZ.push(z, z + F88[1][i]); PV.push(F88[4][i], F88[4][i]); z += F88[1][i]; });
      s.v0 = s.pv.line([200, 200, 300, 300], [0, 10, 10, 22], {color: C.ref, width: 6, draw: true});
      s.v1 = s.pv.line(PV, PZ, {color: C.ssi, width: 7, draw: true});
      s.l200 = s.pv.text(200, 2.5, "200 m/s", {cls: "t-label", size: 30, color: "var(--ref)", dx: 14, hidden: true});
      s.leg = s.pv.legend([{label: "soil report", color: C.ref}, {label: "this earthquake", color: C.ssi}], {x: 785, y: top + 22 * pxm - 68, size: 28, hidden: true});
      // the softest sand sublayer, from FILE88
      let im = 0;
      SAND_ROWS.forEach((i) => { if (F88[4][i] < F88[4][im]) im = i; });
      const vMin = F88[4][im], zMid = PZ[2 * im] + F88[1][im] / 2;
      s.vMin = Math.round(vMin);
      s.a74 = K.g(s.pv.g, {hidden: true});
      const ax = s.pv.X(vMin), ay = s.pv.Y(zMid);
      K.circle(s.a74, ax, ay, 11, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      K.line(s.a74, ax - 4, ay + 14, ax - 24, ay + 52, {stroke: C.ssi, "stroke-width": 2.5});
      K.text(s.a74, ax - 60, ay + 84, `${s.vMin} m/s at ${zMid} m`, {cls: "t-label", size: 32, color: "var(--ssi)", weight: 800});
      s.st = K.stat(s.root, {x: 1390, y: 290, w: 440, value: `200 → ${s.vMin}`, label: `metres per second,<br>sand at ${zMid} m deep`, color: "var(--ssi)", vsize: 84});
      s.card = K.card(s.root, {x: 1400, y: 560, w: 420, title: "Strain-compatible", size: 32, body: "the soil as this earthquake shakes it"});
      s.A = 0; s.ph = 0;
    },
    tick(s, t) {
      // the last round blends the soil-report solution into the converged one
      const W = Math.max(0, Math.min(1, s.ph - 7));
      s.soil.shear((d) => s.A * ((1 - W) * s.mLow(d, t) + W * s.mSc(d, t)));
      const txt = W < 0.5 ? `soil-report stiffness: ${F_LOW.toFixed(2)} Hz` : `softened by this quake: ${CASE.be.f} Hz`;
      if (s.nfTxt !== txt) { s.nfTxt = txt; s.nf.textContent = txt; s.nf.style.fill = W < 0.5 ? "var(--ink2)" : "var(--ssi)"; }
      const lap = Math.min(8, 1 + Math.floor(s.ph)), a = -Math.PI / 2 + TAU * (s.ph >= 8 ? 0 : s.ph % 1);
      s.tok.setAttribute("cx", s.cx + s.R * Math.cos(a)); s.tok.setAttribute("cy", s.cy + s.R * Math.sin(a));
      if (s.lap !== lap) { s.lap = lap; s.iter.innerHTML = `round <b>${lap}</b>`; }
      s.iter.style.color = s.ph >= 8 ? "var(--good)" : "";
    },
    beats: [
      {say: "[a]That's a chicken-and-egg problem: the softness sets the shaking, and the shaking sets the softness.",
        go(k) { k.show(k.s.head); },
        a(k) { k.show([k.s.soil.g, k.s.labs, k.s.inp, k.s.nf, k.s.nx]); k.tween(k.s, {A: 1}, 800); k.show(k.s.ce, {delay: 400}); }},
      {say: "[b]The module SOIL solves it the way you set a shower: try, check, adjust.",
        b(k) { k.hide(k.s.ce); k.show([k.s.ring, k.s.tag], {delay: 300}); }},
      {say: "[c]It shakes the column with the record, [d]measures each layer's strain, [e]and softens each layer from the curves.",
        c(k) { k.show(k.s.n[0]); k.pulse(k.s.inp, {amp: 0.06}); }, d(k) { k.show([k.s.arrs[0].g, k.s.n[1]]); }, e(k) { k.show([k.s.arrs[1].g, k.s.n[2]]); }},
      {say: "[f]Then it shakes again, with the softer soil. [g]After eight rounds here, the answers barely change.",
        f(k) { k.show([k.s.arrs[2].g, k.s.tok]); k.hide(k.s.tag); k.show(k.s.iter); k.tween(k.s, {ph: 8}, 2800, {ease: "linear"}); },
        g(k) { k.pulse(k.s.iter, {amp: 0.08}); k.show(k.s.eqlin, {delay: 400}); }},
      {say: "[p]The result? The soil report gave the sand a shear-wave speed of 200 metres per second.",
        p(k) { k.hide([k.s.ring, k.s.n, k.s.arrs.map((a) => a.g), k.s.tok, k.s.iter, k.s.eqlin]); k.show(k.s.pv.g, {delay: 200}); k.draw(k.s.v0, 1000, {delay: 300}); k.show(k.s.l200, {delay: 400}); }},
      {say: "[q]In this earthquake, at 9.5 metres deep, it drops to just 74.",
        q(k) { k.draw(k.s.v1, 1000); k.show(k.s.leg); k.show([k.s.a74, k.s.st.el], {delay: 400}); }},
      {say: "[c]That's the strain-compatible soil: the soil as your earthquake really shakes it.", c(k) { k.show(k.s.card); }},
    ],
  });

  // ---------------------------------------------------------------- 5. what the soft soil does to the motion
  scenes.push({
    id: "motion", title: "Soft soil reshapes the shaking",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Soft soil reshapes the shaking", {x: 120, y: 70, size: "h2"});
      const TX = 200, TW = 800, PX = 120 / SURF.peak, Y1 = 300, Y2 = 610;           // one scale for both: PX px per g
      s.TX = TX; s.TW = TW;
      /** A history from 0 to 20 s across TW px, PX px per g, centred on yc. */
      const trace = (gg, d, yc, color) => {
        const pts = [];
        d.t.forEach((t, i) => { if (t <= 20) pts.push([TX + (TW * t) / 20, yc - PX * d.v[i]]); });
        return K.path(gg, K.d(pts), {stroke: color, "stroke-width": 3, fill: "none", "stroke-linejoin": "round", draw: true});
      };
      s.hr = K.g(g, {hidden: true});
      K.line(s.hr, TX, Y1, TX + TW, Y1, {stroke: "#50657f", "stroke-width": 2});
      K.text(s.hr, TX, 215, "rock below", {cls: "t-label", size: 34, color: "var(--wave)", weight: 750});
      s.tr = trace(s.hr, REC, Y1, C.wave);
      s.vr = K.text(g, TX + TW, 215, `peak ${REC.peak.toFixed(3)} g`, {cls: "t-label", size: 34, anchor: "end", color: "var(--wave)", weight: 800, hidden: true});
      s.hs = K.g(g, {hidden: true});
      K.line(s.hs, TX, Y2, TX + TW, Y2, {stroke: "#50657f", "stroke-width": 2});
      K.text(s.hs, TX, 470, "ground surface", {cls: "t-label", size: 34, color: "var(--ssi)", weight: 750});
      s.ts = trace(s.hs, SURF, Y2, C.ssi);
      s.vs = K.text(g, TX + TW, 470, `peak ${SURF.peak.toFixed(3)} g`, {cls: "t-label", size: 34, anchor: "end", color: "var(--ssi)", weight: 800, hidden: true});
      s.x2 = K.h("div", {x: 560, y: 432, w: 260, align: "center", html: `× ${(SURF.peak / REC.peak).toFixed(2)}`, size: 40, in: "pop", style: {fontWeight: "800", color: "var(--ssi)"}}, s.root);
      // the scales: 0.5 g bar, time axis
      s.sc = K.g(g, {hidden: true});
      K.line(s.sc, TX - 30, Y2, TX - 30, Y2 - 0.5 * PX, {stroke: C.ink2, "stroke-width": 4});
      [Y2, Y2 - 0.5 * PX].forEach((y) => K.line(s.sc, TX - 40, y, TX - 20, y, {stroke: C.ink2, "stroke-width": 3}));
      K.text(s.sc, TX - 48, Y2 - 0.25 * PX + 9, "0.5 g", {cls: "tick", anchor: "end"});
      K.line(s.sc, TX, 756, TX + TW, 756, {stroke: "#50657f", "stroke-width": 2});
      [0, 5, 10, 15, 20].forEach((v) => {
        K.line(s.sc, TX + (TW * v) / 20, 756, TX + (TW * v) / 20, 768, {stroke: "#50657f", "stroke-width": 2});
        K.text(s.sc, TX + (TW * v) / 20, 796, String(v), {cls: "tick", anchor: "middle"});
      });
      K.text(s.sc, TX + TW / 2, 832, "time (s) · same scale for both", {cls: "axlab", anchor: "middle"});
      s.cur = K.line(g, TX, 200, TX, 756, {stroke: "rgba(238,243,249,.35)", "stroke-width": 2});
      s.C = 0;
      // the spectra
      s.p = K.plot(s.svg, {x: 1200, y: 200, w: 600, h: 460, xr: [0.1, 100], xlog: true, yr: [0, 2.5], xticks: [0.1, 1, 10, 100], yticks: [0, 1, 2],
        xlabel: "frequency (Hz)", ylabel: "SA (g), 5 %", ylabelOffset: 70});
      s.r0 = s.p.line(RSO[0], RSO[1], {color: C.wave, width: 5, draw: true});
      s.r1 = s.p.line(RS1.f, RS1.sa, {color: C.ssi, width: 6, draw: true});
      s.l0 = s.p.legend([{label: "rock", color: C.wave}], {x: 1580, y: 236, size: 28, hidden: true});
      s.l1 = s.p.legend([{label: "surface", color: C.ssi}], {x: 1580, y: 278, size: 28, hidden: true});
      s.m0 = K.g(s.p.g, {hidden: true});
      K.line(s.m0, s.p.X(PK_ROCK.x), s.p.Y(0), s.p.X(PK_ROCK.x), s.p.Y(PK_ROCK.y), {stroke: C.wave, "stroke-width": 3, "stroke-dasharray": "8 7"});
      K.circle(s.m0, s.p.X(PK_ROCK.x), s.p.Y(PK_ROCK.y), 10, {fill: C.wave, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.m0, s.p.X(PK_ROCK.x) + 14, s.p.Y(0.3), "≈ 2.5 Hz", {cls: "t-label", size: 32, color: "var(--wave)", weight: 800});
      s.m1 = K.g(s.p.g, {hidden: true});
      K.circle(s.m1, s.p.X(PK_BE.x), s.p.Y(PK_BE.y), 11, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.m1, s.p.X(PK_BE.x) - 20, s.p.Y(PK_BE.y) + 8, `${PK_BE.x.toFixed(2)} Hz`, {cls: "t-label", size: 32, anchor: "end", color: "var(--ssi)", weight: 800});
      s.warn = K.card(s.root, {x: 140, y: 858, w: 1660, kind: "warn", title: "Why it matters", size: 34, body: "Skip SOIL, and your ISRS miss both changes."});
    },
    tick(s, t) {
      const x = s.TX + s.TW * ((t * 0.09) % 1);
      s.cur.setAttribute("x1", x); s.cur.setAttribute("x2", x);
      s.cur.style.opacity = String(s.C);
    },
    beats: [
      {say: "[r]The rock peaks at 0.324 g. [s]The surface peaks at 0.605 g: nearly double.",
        go(k) { k.show(k.s.head); },
        r(k) { k.show(k.s.hr); k.draw(k.s.tr, 900); k.show(k.s.vr, {delay: 300}); },
        s(k) { k.show([k.s.hs, k.s.sc]); k.draw(k.s.ts, 900); k.show([k.s.vs, k.s.x2], {delay: 300}); k.tween(k.s, {C: 1}, 600, {delay: 400}); }},
      {say: "[p]Now the spectra. [t]On rock, the peak sits at about 2.5 hertz.",
        p(k) { k.show(k.s.p.g); k.draw(k.s.r0, 1000); k.show(k.s.l0); }, t(k) { k.show(k.s.m0); }},
      {say: "[u]At the surface, it moves down to 1.86 hertz, and it's much taller.",
        u(k) { k.draw(k.s.r1, 1000); k.show(k.s.l1); k.show(k.s.m1, {delay: 400}); }},
      {say: "[w]That's the shaking your building feels. Skip SOIL, and your floor spectra, the ISRS, miss both changes.",
        w(k) { k.show(k.s.warn); }},
    ],
  });

  // ---------------------------------------------------------------- 6. how SASSI does it
  scenes.push({
    id: "sassi", title: "How SASSI does it",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "How SASSI does it", {x: 120, y: 70, size: "h2"});
      const BX = [150, 580, 1010, 1440], BY = 200, BW = 330, BH = 130;
      const box = (i, name, sub, color) => K.h("div", {x: BX[i], y: BY, w: BW, h: BH, in: "up",
        html: `<div style="font:700 40px var(--mono);color:${color}">${name}</div><div style="font-size:28px;color:var(--ink2);margin-top:6px">${sub}</div>`,
        style: {border: `3px solid ${color}`, borderRadius: "16px", background: "rgba(20,32,51,.92)", display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "center"}}, s.root);
      s.b = [box(0, "EQUAKE", "makes the record", "var(--wave)"), box(1, "SOIL", "softens the soil", "var(--ssi)"),
        box(2, "SITE", "builds the site", "var(--spring)"), box(3, "SSI run", "your building", "var(--ink)")];
      s.a = [0, 1, 2].map((i) => K.arrow(g, BX[i] + BW + 16, BY + BH / 2, BX[i + 1] - 14, BY + BH / 2, {color: "#5f7da6", width: 5, head: 18, hidden: true}));
      s.code = K.code(s.root, ["SITEX,1", "RUNSITE"], {x: 1010, y: 410, w: 420, size: 40});
      s.expl = K.h("div", {x: 1010, y: 650, w: 820, html: "→ use SOIL's layers, not the report's", size: 34, in: "left", style: {color: "var(--ssi)", fontWeight: "650"}}, s.root);
      // the trap: the SSI motion is defined at the ground surface here
      const pxm = 10, X = 60;
      s.mc = K.g(g, {hidden: true});
      s.soil = column(s.mc, {x: 180, y: 470, w: 220, pxm, hs: 60});
      s.cp = K.circle(s.mc, 290, 470, 13, {fill: C.wave, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.mc, 450, 482, "SSI motion defined here", {cls: "t-label", size: 32, color: "var(--wave)", weight: 750});
      K.text(s.mc, 180, 790, `1-D waves, ${CASE.be.f} Hz\ndisplacements × ${X} · slowed ${SLOW} ×`, {cls: "t-label", size: 24, color: "var(--muted)"});
      s.mv = motion(COL.be, CASE.be.f, CASE.be.pga, pxm, X, 22 + 60 / pxm);
      s.ok = K.h("div", {cls: "sv-pill", html: `✓ SOIL's surface motion · ${SURF.peak.toFixed(3)} g`, x: 450, y: 560, size: 30, color: "var(--good)", in: "pop"}, s.root);
      s.bad = K.h("div", {cls: "sv-pill", html: `<s>the rock record · ${REC.peak.toFixed(3)} g</s>`, x: 450, y: 650, size: 30, color: "var(--bad)", in: "pop"}, s.root);
      s.A = 0;
    },
    tick(s, t) {
      const u = (d) => s.A * s.mv(d, t);
      s.soil.shear(u);
      s.cp.setAttribute("cx", 290 + u(0));
    },
    beats: [
      {say: "[a]In SASSI, it's a short chain of modules, run one after the other.",
        go(k) { k.show(k.s.head); },
        a(k) { const s = k.s; k.show([s.b[0], s.a[0], s.b[1], s.a[1], s.b[2], s.a[2], s.b[3]], {stagger: 120}); }},
      {say: "[c]SITE then builds the site for the SSI run, [d]and one short command hands it the softened soil.",
        c(k) { k.pulse(k.s.b[2], {amp: 0.06}); }, d(k) { k.show(k.s.code.el); k.type(k.s.code.lines, {cps: 18}); }},
      {say: "[k]{SITEX,1|Site X one} means: use SOIL's softened layers, not the soil report's.", k(k) { k.show(k.s.expl); }},
      {say: "[t]One trap: in this example, the SSI motion is defined at the ground surface.",
        t(k) { k.show(k.s.mc); k.tween(k.s, {A: 1}, 800); }},
      {say: "[x]So drive the SSI run with SOIL's surface motion, 0.605 g, [y]not the 0.324 g rock record.",
        x(k) { k.show(k.s.ok); }, y(k) { k.show(k.s.bad); }},
    ],
  });

  // ---------------------------------------------------------------- 7. soil cases and their envelope
  scenes.push({
    id: "cases", title: "Bracket the soil",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Nobody knows the soil exactly", {x: 120, y: 70, size: "h2"});
      // three soil cases: each column moves with its own strain-compatible solution at its own site frequency
      const pxm = 14, X = 30, hs = 60;
      s.cases = [["softer", C.violet, "lb"], ["best estimate", C.ssi, "be"], ["stiffer", C.wave, "ub"]].map(([lab, col, key], i) => {
        const x = 160 + i * 290, cs = CASE[key];
        const grp = K.g(g, {hidden: true});
        K.rect(grp, x - 40, 270, 250, 420, {rx: 16, fill: "none", stroke: col, "stroke-width": 3});
        const so = column(grp, {x, y: 300, w: 170, pxm, hs});
        K.text(grp, x + 85, 740, lab, {cls: "t-label", anchor: "middle", size: 32, color: col, weight: 750});
        const fl = K.text(grp, x + 85, 784, `${cs.f.toFixed(2)} Hz`, {cls: "t-label", anchor: "middle", size: 30, color: col});
        return {g: grp, so, fl, mv: motion(COL[key], cs.f, cs.pga, pxm, X, 22 + hs / pxm)};
      });
      s.note = note(g, 120, 850, `each at its site frequency: 1-D waves\ndisplacements × ${X} · slowed ${SLOW} ×`);
      // the surface spectra of the three SOIL runs, and their envelope (BROADEN, no broadening)
      s.p = K.plot(s.svg, {x: 1080, y: 200, w: 720, h: 470, xr: [0.1, 100], xlog: true, yr: [0, 3.5], xticks: [0.1, 1, 10, 100], yticks: [0, 1, 2, 3],
        xlabel: "frequency (Hz)", ylabel: "surface SA (g), 5 %", ylabelOffset: 70});
      s.env = s.p.line(ENV.f, ENV.sa, {color: "rgba(238,243,249,.3)", width: 20, draw: true});
      s.lb = s.p.line(LB.f, LB.sa, {color: C.violet, width: 5, draw: true});
      s.be = s.p.line(RS1.f, RS1.sa, {color: C.ssi, width: 5, draw: true});
      s.ub = s.p.line(UB.f, UB.sa, {color: C.wave, width: 5, draw: true});
      s.leg = s.p.legend([{label: "softer", color: C.violet}, {label: "best estimate", color: C.ssi}, {label: "stiffer", color: C.wave}], {x: 1560, y: 300, dy: 42, size: 28, hidden: true});
      s.envLeg = s.p.legend([{label: "envelope", color: "rgba(238,243,249,.45)"}], {x: 1560, y: 430, size: 28, hidden: true});
      // the input's 2.5 Hz peak (the RG 1.60 corner) and the two peaks it names
      s.v25 = K.g(s.p.g, {hidden: true});
      K.line(s.v25, s.p.X(PK_TGT.x), s.p.Y(0), s.p.X(PK_TGT.x), s.p.Y(3.5), {stroke: C.muted, "stroke-width": 3, "stroke-dasharray": "10 8"});
      K.text(s.v25, s.p.X(PK_TGT.x) + 12, s.p.Y(0.3), `input peak ${PK_TGT.x} Hz`, {cls: "t-label", size: 26, color: "var(--muted)"});
      s.pk = K.g(s.p.g, {hidden: true});
      K.circle(s.pk, s.p.X(PK_BE.x), s.p.Y(PK_BE.y), 10, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      K.circle(s.pk, s.p.X(PK_UB.x), s.p.Y(PK_UB.y), 11, {fill: C.wave, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.pk, s.p.X(PK_UB.x) + 20, s.p.Y(PK_UB.y) + 10, `+${Math.round((PK_UB.y / PK_BE.y - 1) * 100)} %`, {cls: "t-label", size: 38, color: "var(--wave)", weight: 800});
      s.A = 0;
    },
    tick(s, t) {
      s.cases.forEach((c) => c.so.shear((d) => s.A * c.mv(d, t)));
    },
    beats: [
      {say: "[a]Last idea: nobody knows the soil exactly.", go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.cases[1].g); k.tween(k.s, {A: 1}, 800); }},
      {say: "[b]So you bracket it, the way you pack for colder and warmer weather: [c]softer, best estimate, stiffer.",
        b(k) { k.show(k.s.cases[0].g); }, c(k) { k.show(k.s.cases[2].g); k.show(k.s.note, {delay: 400}); }},
      {say: "[c]Each case gets its own SOIL run, and its own surface spectrum.",
        c(k) { k.show(k.s.p.g); k.show(k.s.leg, {delay: 150}); k.draw(k.s.lb, 900); k.draw(k.s.be, 900, {delay: 200}); k.draw(k.s.ub, 900, {delay: 400}); }},
      {say: "[u]The stiffer soil's natural frequency lands on the input's 2.5 hertz peak, [p]so its spectrum peaks 40 percent higher.",
        u(k) { k.pulse([k.s.cases[2].g, k.s.cases[2].fl], {amp: 0.04}); k.show(k.s.v25); }, p(k) { k.show(k.s.pk); }},
      {say: "[e]So you keep the envelope: the highest of the three, at every frequency.",
        e(k) { k.draw(k.s.env, 1000); k.show(k.s.envLeg); }},
    ],
  });

  // ---------------------------------------------------------------- 8. recap
  scenes.push({
    id: "recap", title: "Recap",
    build(s) {
      s.head = K.heading(s.root, "Recap", {x: 120, y: 90, size: "h1"});
      s.list = K.bullets(s.root, [
        {t: "EQUAKE: a record that matches the spectrum", sub: "noise, tuned until its spectrum fits"},
        {t: "SOIL: how soft the soil gets", sub: "guess, check, repeat: equivalent-linear"},
        {t: "Bracket the soil, keep the envelope", sub: "softer, best estimate, stiffer"},
      ], {x: 120, y: 240, w: 1000, num: true, size: 40});
      s.q = K.card(s.root, {x: 1180, y: 240, w: 640, kind: "check", title: "Check yourself", size: 34,
        body: "The soil report says 200 m/s. Should SITE use it for this earthquake?"});
      s.a = K.card(s.root, {x: 1180, y: 560, w: 640, title: "Answer", size: 34,
        body: "No: that's for gentle shaking. Run SOIL first: here the sand drops to 74 m/s."});
      s.next = K.pill(s.root, "Next · Lesson 7: three earthquake directions and design ISRS →", {x: 120, y: 900, size: 30, color: "var(--wave)"});
    },
    beats: [
      {say: "[a]Let's recap. One: EQUAKE turns your design spectrum into a matching record.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[a]Two: SOIL finds how soft the soil gets in that earthquake. Guess, check, repeat.", a(k) { k.show(k.s.list.items[1]); }},
      {say: "[a]Three: run softer, best and stiffer soil, and keep the envelope.", a(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Check yourself. The soil report says 200 metres per second. Should SITE use that for this earthquake?",
        q(k) { k.show(k.s.q); }, gap: 1500},
      {say: "[a]No. That's the gentle-shaking value. Run SOIL first: here the sand drops to 74, and the motion changes.",
        a(k) { k.show(k.s.a); }},
      {say: "[n]Next: earthquakes shake in three directions at once, and you'll turn the results into design floor spectra.",
        n(k) { k.show(k.s.next); }},
    ],
  });

  SV.video({
    id: "06", n: 6, part: "Design applications", lesson: "06-seismic-input",
    title: "Seismic input: spectrum-compatible motion and strain-compatible soil",
    subtitle: "The two inputs of every SSI run: a record that matches the design spectrum, and the soil as that earthquake shakes it.",
    next: {href: "07.html", title: "Three earthquake directions and design ISRS"},
    scenes,
  });
})();
