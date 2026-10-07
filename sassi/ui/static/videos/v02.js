/* Explainer video, lesson 2 (newcomer edition): The site -- layers, half-space and free-field motion.
 * One big idea: the earthquake reaches the building through layers of soil that amplify some frequencies;
 * SITE computes this free field for every layer and frequency, and the half-space below lets waves leave.
 * Numbers: sassi/ui/lessons/02_free_field.md; curves: d02.js (the lesson's SOIL runs, python -m sassi.ui.video_data)
 * and the exact one-dimensional column of the lesson's soil-column figure (SV.K.PH.columnU, figures.js).
 * Geometry to scale: the example 1 site (16 ft sand, 39 ft gravel, rock) and its building (64 ft mat 5 ft thick,
 * four 16 ft storeys).  Every moving column follows the 1D wave solution at the frequency stated on screen,
 * with one displacement factor per scene (a surface or rock motion of 0.04 in, drawn larger); the building of
 * the hook moves with the one-mode SSI model of example 1 (PH.ssiResponse, lessons 1 and 4) at the same
 * frequency. */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["02"], PH = K.PH, COL = PH.COLUMN;
  const TAU = 2 * Math.PI;
  const OUT = D["ex01/SAF001W_023O.TFU"], WIT = D["ex01/SAF001W_023W.TFU"], ACC = D["data/rg160h_030g.acc"];
  const DF = 0.0244140625;                                  // SITE frequency step (Hz)
  const F01 = 4 * DF, F20 = 819 * DF;                        // lowest / highest SITE frequency
  const SLOW = 10;                                          // animations at 7.3-7.4 Hz run 10 x slower
  /** Mode 1 of the example 1 stick (floors at 16, 32, 48, 64 ft, 2,200 kips each; Timoshenko beams of the input,
   *  clamped at the base): Gamma1 phi1 at the floors.  4.99 Hz, M1 = 6,750 kips, h1 = 49.8 ft, 1.34 at the roof. */
  const GPHI = [0.2157, 0.5578, 0.9542, 1.3423];
  /** Pixels of a displayed displacement: `inch` inches drawn `times` larger, at `pxf` px per ft. */
  const disp = (inch, times, pxf) => inch / 12 * times * pxf;

  // ---------------------------------------------------------------- shared pieces
  /** Exact free field of the lesson's column (SHAKE recursion, figures.js) per unit surface motion at f
   *  (rounded to 0.1 Hz): a function of the depth z (ft) returning the complex amplitude [re, im]. */
  const PROF = {};
  function profile(f) {
    const key = Math.max(0.1, Math.round(f * 10) / 10).toFixed(1);
    if (PROF[key]) return PROF[key];
    const wv = PH.columnRatios(COL, +key).wv, dz = 0.8, us = [];
    for (let z = 0; z <= 131.2001; z += dz) us.push(PH.columnU(wv, COL, z));
    const fn = (z) => {
      const x = Math.max(0, Math.min(us.length - 1.001, z / dz)), i = Math.floor(x), a = x - i;
      return [us[i][0] * (1 - a) + us[i + 1][0] * a, us[i][1] * (1 - a) + us[i + 1][1] * a];
    };
    fn.wv = wv;
    PROF[key] = fn;
    return fn;
  }
  const re = (c, ph) => c[0] * Math.cos(ph) - c[1] * Math.sin(ph);      // Re(c e^{i ph})
  const mag = (c) => Math.hypot(c[0], c[1]);
  /** Linear interpolation of a computed curve. */
  function interp(X, Y, x) {
    if (x <= X[0]) return Y[0];
    for (let i = 1; i < X.length; i++) if (X[i] >= x) { const a = (x - X[i - 1]) / (X[i] - X[i - 1]); return Y[i - 1] * (1 - a) + Y[i] * a; }
    return Y[Y.length - 1];
  }
  /** The largest value of a computed curve between a and b: {f, a}. */
  function peak(X, Y, a, b) {
    let best = null;
    X.forEach((x, i) => { if (x >= a && x <= b && (!best || Y[i] > best.a)) best = {f: x, a: Y[i]}; });
    return best;
  }
  /** The site of example 1: 16 ft of sand on 39 ft of gravel on rock (pxf: px per ft). */
  function column(p, o) {
    const pxf = o.pxf;
    const col = K.soil(p, {x: o.x, y: o.y, w: o.w, layers: [{h: 16 * pxf, kind: "sand"}, {h: 39 * pxf, kind: "gravel"}],
      hs: o.hsH === 0 ? null : {h: o.hsH || 120, kind: "rock"}, labels: false, hidden: o.hidden});
    col.base = o.y + 55 * pxf;
    return col;
  }
  /** Shear a column with the free field U (per unit surface motion): surface amplitude A px, phase ph. Returns u(depth px). */
  function shake(col, U, pxf, A, ph) {
    const u = (d) => A * re(U(d / pxf), ph);
    col.shear(u);
    return u;
  }
  function halo(e, w) {
    e.style.paintOrder = "stroke"; e.style.stroke = "rgba(10,17,29,.85)"; e.style.strokeWidth = (w || 6) + "px"; e.style.strokeLinejoin = "round";
    return e;
  }
  /** A word on a layer (halo, large). */
  const onLayer = (p, x, y, str, o) => halo(K.text(p, x, y, str, Object.assign({cls: "t-label", size: 32, color: "var(--ink)"}, o || {})));
  /** The scale note of an animated scene (one displacement factor, the time factor). */
  const note = (p, x, y, str, o) => halo(K.text(p, x, y, str, Object.assign({cls: "t-small", size: 23}, o || {})), 5);
  /** An arrow built with {draw: true} whose head waits for the line: hide the head at build, then rise(). */
  function arrowD(p, x1, y1, x2, y2, o) {
    const a = K.arrow(p, x1, y1, x2, y2, Object.assign({draw: true}, o));
    a.head.classList.add("sv-in", "sv-fade");
    return a;
  }
  const rise = (k, a, ms) => k.draw(a.line, ms, {done: () => k.show(a.head)});
  const IFC = [].concat(K.linspace(0, 16, 11), K.linspace(19.25, 55, 12));   // the 23 interface depths of SITE (ft)

  // ================================================================== scenes
  const scenes = [];

  // ---------------------------------------------------------------- 0. hook
  // to scale, 5 px per ft: the example 1 building (64 ft mat, four 16 ft storeys) on its site
  const HP = 5, HA = disp(0.04, 2000, HP), HF = 7.4;          // px per ft; 0.04 in of surface motion x 2000 = 33 px
  scenes.push({
    id: "hook", title: "Where the shaking comes from",
    build(s) {
      const g = K.g(s.svg);
      s.U = profile(HF);
      s.ssi = PH.ssiResponse(PH.ssiModel(1), HF);            // the building's motion per unit free-field motion
      s.head = K.heading(s.root, "Where does the shaking come from?", {x: 120, y: 80});
      const stick = (p, x, y) => K.stick(p, {x, y, z: [16 * HP, 32 * HP, 48 * HP, 64 * HP], matW: 64 * HP, matH: 5 * HP, r: 16, slabW: 100});
      // left: the fixed-base model
      s.L = K.g(g, {hidden: true});
      K.text(s.L, 470, 210, "Fixed base", {cls: "t-label", anchor: "middle", size: 36, color: "var(--ref)"});
      K.ground(s.L, 290, 640, 360);
      s.fb = stick(s.L, 470, 640 - 5 * HP);
      s.acc = K.signal(s.L, {x: 290, y: 690, w: 360, h: 100, color: C.wave, width: 3, draw: true, n: 260,
        fn: (u) => { const i = Math.min(ACC.v.length - 1, Math.floor(u * ACC.v.length * 0.55)); return ACC.v[i] / ACC.peak; }});
      s.accLab = K.text(s.L, 470, 840, "earthquake in, at the base", {cls: "t-label", anchor: "middle", size: 28, hidden: true});
      // right: the same building on the soil site
      s.R = K.g(g, {hidden: true});
      K.text(s.R, 1340, 210, "On a soil site", {cls: "t-label", anchor: "middle", size: 36, color: "var(--wave)"});
      s.col = column(s.R, {x: 1040, y: 600, w: 600, pxf: HP, hsH: 115});
      s.glue = [K.path(s.R, "", {fill: C.sand, "fill-opacity": 0.78}), K.path(s.R, "", {fill: "url(#sv-pat-sand)"})];   // soil stuck to the mat
      s.st = stick(s.R, 1340, 600 - 5 * HP);
      onLayer(s.R, 1060, 650, "sand", {size: 30}); onLayer(s.R, 1060, 785, "gravel", {size: 30}); onLayer(s.R, 1060, 940, "rock", {size: 30});
      s.rise = arrowD(s.R, 1700, 985, 1700, 612, {color: C.wave, width: 6, head: 20});
      K.text(s.R, 1722, 780, "waves\nrise", {cls: "t-label", size: 28, color: "var(--wave)"});
      note(s.R, 1830, 290, `${HF} Hz · surface 0.04 in`, {anchor: "end"});
      note(s.R, 1830, 318, "drawn × 2000 · 10 × slower", {anchor: "end"});
      note(s.R, 1830, 346, "building: simplified model", {anchor: "end"});
      s.q = K.heading(s.root, "Weaker…<br>or <span style='color:var(--ssi)'>stronger</span>?", {x: 160, y: 430, size: "h1"});
      s.S = 0;
    },
    tick(s, t) {
      const ph = TAU * HF * t / SLOW, A = HA * s.S, r = s.ssi;
      shake(s.col, s.U, HP, A, ph);
      // th is the rotation per unit length of ground motion in the length unit of PH.SSI, whose mat width is B:
      // the ground motion drawn as A px is A * B / (64 * HP) of those units
      s.st.set({sway: A * re([1 + r.u0[0], r.u0[1]], ph), rock: A * PH.SSI.B / (64 * HP) * re(r.th, ph), bend: GPHI.map((gp) => A * gp * re(r.u, ph))});
      // the soil under the mat follows it (welded contact): fill between the ground surface and the mat's underside
      const st = s.st.state, c = Math.cos(st.rock), sn = Math.sin(st.rock), h = 2.5 * HP;
      const P = (lx) => [1340 + st.sway + lx * c - h * sn, 600 - h + lx * sn + h * c];
      const a = P(-32 * HP), b = P(32 * HP), d = `M${a}L${b}L${b[0]},600L${a[0]},600Z`;
      s.glue.forEach((e) => e.setAttribute("d", d));
    },
    beats: [
      {say: "In a fixed-base model, the earthquake goes [acc]straight into the base of your building.",
        go(k) { k.show([k.s.head, k.s.L]); }, acc(k) { k.draw(k.s.acc, 1000); k.show(k.s.accLab, {delay: 400}); }},
      {say: "On a real soil site, [r]the waves start deeper down, in the rock.",
        r(k) { k.show(k.s.R); rise(k, k.s.rise, 900); k.tween(k.s, {S: 1}, 1000, {delay: 300}); }},
      {say: "So does the soil [q]weaken the shaking, or make it stronger? And how does SASSI know?",
        q(k) { k.hide(k.s.L); k.show(k.s.q, {delay: 200}); }},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 2",
    build(s) {
      s.t = K.titleCard(s.root, {n: 2, part: "Fundamentals", title: "The site: layers, half-space and free-field motion",
        sub: "How layers of soil change the earthquake, and how SASSI computes it."});
    },
    beats: [
      {say: "Lesson two: the site. Layers, half-space, and the free-field motion.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 100}); k.show(t.title, {delay: 200}); k.show(t.rule, {delay: 300}); k.show(t.sub, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. the site and the wave speed
  scenes.push({
    id: "site", title: "Meet the site",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Meet the site", {x: 120, y: 70});
      const SF = 9.25;                                                          // px per ft
      s.col = column(g, {x: 250, y: 230, w: 420, pxf: SF, hsH: 211});         // sand 230-378, gravel 378-739, rock 739-950
      s.col.layers.forEach((L) => L.g.classList.add("sv-in", "sv-fade"));
      const [Ls, Lg, Lr] = s.col.layers;
      onLayer(Ls.g, 280, 318, "Sand");
      K.dim(Ls.g, 216, 230, 216, 230 + 16 * SF, "");
      K.text(Ls.g, 196, 316, "16 ft", {cls: "t-label", anchor: "end", size: 32, color: "var(--ink)"});
      onLayer(Lg.g, 280, 572, "Gravel");
      K.dim(Lg.g, 216, 230 + 16 * SF, 216, 230 + 55 * SF, "");
      K.text(Lg.g, 196, 572, "39 ft", {cls: "t-label", anchor: "end", size: 32, color: "var(--ink)"});
      onLayer(Lr.g, 280, 830, "Rock");
      // the race of the shear waves: speeds to scale
      s.lanes = [["Sand", 1000, "1,000 ft/s", 400], ["Gravel", 1650, "1,650 ft/s", 560], ["Rock", 3300, "3,300 ft/s", 720]].map(([name, vs, lab, y]) => {
        const lg = K.g(g, {hidden: true});
        K.text(lg, 900, y + 10, name, {cls: "t-label", size: 32, color: "var(--ink)"});
        K.line(lg, 1070, y, 1630, y, {stroke: "rgba(185,196,208,.35)", "stroke-width": 3, "stroke-linecap": "round"});
        const pulse = K.path(lg, "", {stroke: C.wave, "stroke-width": 5, fill: "none", "stroke-linecap": "round"});
        const v = K.text(lg, 1660, y + 10, lab, {cls: "t-label", size: 32, color: "var(--ink)", hidden: true});
        return {g: lg, pulse, v, vs, y};
      });
      s.lanesHead = K.text(g, 900, 300, "how fast shear waves run (speeds to scale)", {cls: "t-label", size: 30, hidden: true});
      s.term = K.card(s.root, {x: 900, y: 800, w: 900, title: "Shear-wave velocity, Vs",
        body: "slow means soft soil · fast means stiff rock"});
    },
    tick(s, t) {
      for (const L of s.lanes) {
        const x0 = 1070, len = 560, xp = x0 + ((t * L.vs * 0.096 + 120) % len), pts = [];
        for (let dx = -60; dx <= 60; dx += 4) {
          const x = xp + dx;
          if (x < x0 || x > x0 + len) continue;
          pts.push([x, L.y - 34 * Math.exp(-(dx / 22) * (dx / 22))]);
        }
        L.pulse.setAttribute("d", pts.length > 1 ? K.d(pts) : "");
      }
    },
    beats: [
      {say: "Here's the site of example one: [s]16 ft of sand, [g]on 39 ft of gravel, [r]on rock.",
        go(k) { k.show(k.s.head); }, s(k) { k.show(k.s.col.layers[0].g); }, g(k) { k.show(k.s.col.layers[1].g); }, r(k) { k.show(k.s.col.layers[2].g); }},
      {say: "How stiff is each layer? We measure it by [v]how fast shear waves run through it.",
        v(k) { k.show(k.s.lanesHead); k.show(k.s.lanes.map((L) => L.g), {stagger: 150}); k.show(k.s.lanes.map((L) => L.v), {stagger: 150, delay: 400}); }},
      {say: "[t]That speed is called the shear-wave velocity, Vs: slow means soft, fast means stiff.",
        t(k) { k.show(k.s.term); k.pulse([k.s.lanes[0].v, k.s.lanes[2].v], {amp: 0.2}); }},
    ],
  });

  // ---------------------------------------------------------------- 3. waves rise: the free field
  // to scale, 5.5 px per ft; 7.4 Hz, 0.04 in of surface motion drawn x 2000 = 37 px
  const FP = 5.5, FA = disp(0.04, 2000, FP);
  scenes.push({
    id: "free-field", title: "The free field",
    build(s) {
      const g = K.g(s.svg);
      s.U = profile(HF);
      s.head = K.heading(s.root, "Waves rise through the soil", {x: 120, y: 70});
      s.col = column(g, {x: 160, y: 580, w: 520, pxf: FP, hsH: 80});         // base 882, rock to 962
      s.env = K.path(g, "", {stroke: C.ink2, "stroke-width": 2.5, "stroke-dasharray": "7 7", fill: "none", hidden: true});
      // the vertically incident shear wave: travels up, moves the soil sideways
      s.rise = arrowD(g, 730, 962, 730, 600, {color: C.wave, width: 6, head: 20});
      s.riseLab = K.g(g, {hidden: true});
      K.text(s.riseLab, 760, 700, "shear wave travels up", {cls: "t-label", size: 30, color: "var(--wave)"});
      K.text(s.riseLab, 760, 734, "(vertical SV wave)", {cls: "t-small"});
      s.side = K.g(g, {hidden: true});
      K.arrow(s.side, 760, 820, 880, 820, {color: C.ink, width: 5, head: 16, both: true});
      K.text(s.side, 900, 830, "soil moves sideways", {cls: "t-label", size: 30, color: "var(--ink)"});
      s.scale = note(g, 160, 990, `${HF} Hz · surface 0.04 in, drawn × 2000 · 10 × slower`, {hidden: true});
      // the building of example 1 to scale, later (64 ft mat 5 ft thick, four 16 ft storeys)
      s.ghost = K.g(g, {hidden: true});
      const gs = {stroke: C.muted, "stroke-width": 3, "stroke-dasharray": "10 8", fill: "none"};
      const mt = 580 - 5 * FP;
      K.rect(s.ghost, 420 - 32 * FP, mt, 64 * FP, 5 * FP, gs);
      K.line(s.ghost, 420, mt, 420, mt - 64 * FP, gs);
      for (const z of [16, 32, 48, 64]) { K.line(s.ghost, 370, mt - z * FP, 470, mt - z * FP, gs); K.circle(s.ghost, 420, mt - z * FP, 14, gs); }
      K.text(s.ghost, 500, 300, "your building, later", {cls: "t-label", size: 28, color: "var(--ink2)"});
      s.card = K.card(s.root, {x: 1060, y: 380, w: 720, title: "The free field",
        body: "how the ground shakes with <b>no building</b> on it"});
      s.A = 0;
    },
    tick(s, t) {
      const A = FA * s.A;
      shake(s.col, s.U, FP, A, TAU * HF * t / SLOW);
      // the envelope of the standing wave: how far each depth swings
      const xr = 680, L = [], R = [];
      for (let z = 0; z <= 69.1; z += 1.5) { const a = A * mag(s.U(z)), y = 580 + z * FP; L.push([xr - a, y]); R.push([xr + a, y]); }
      s.env.setAttribute("d", K.d(L) + K.d(R));
    },
    beats: [
      {say: "When an earthquake strikes, [w]shear waves rise from the rock, up through the layers.",
        go(k) { k.show(k.s.head); }, w(k) { rise(k, k.s.rise, 1000); k.show(k.s.riseLab, {delay: 300}); }},
      {say: "[s]The layers sway from side to side, and the top swings more than the bottom.",
        s(k) { k.tween(k.s, {A: 1}, 1000); k.show([k.s.side, k.s.scale]); k.show(k.s.env, {delay: 400}); }},
      {say: "This shaking of the bare ground, [ff]with no building on it, is called the free field.",
        ff(k) { k.show(k.s.card); }},
      {say: "[b]It's the input of every SASSI analysis: the motion waiting for your building.",
        b(k) { k.show(k.s.ghost); k.pulse(k.s.card, {amp: 0.04}); }},
    ],
  });

  // ---------------------------------------------------------------- 4. amplification: the column's own note
  // to scale, 8 px per ft; rock outcrop 0.04 in at every frequency, drawn x 500 = 13 px
  const AP = 8, AA = disp(0.04, 500, AP);
  scenes.push({
    id: "amplification", title: "Stronger at some frequencies",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Stronger or weaker?", {x: 120, y: 70});
      s.col = column(g, {x: 150, y: 300, w: 360, pxf: AP, hsH: 150});         // base 740, rock to 890
      note(g, 150, 940, "rock outcrop 0.04 in at every frequency");
      note(g, 150, 970, "drawn × 500 · 10 × slower");
      const P = K.plot(s.svg, {x: 760, y: 190, w: 1020, h: 520, xr: [0, 25], yr: [0, 2.5], xticks: [0, 5, 10, 15, 20, 25], yticks: [0, 1, 2],
        xlabel: "frequency (Hz)", ylabel: "surface ÷ bare rock", ylabelOffset: 70});
      s.P = P;
      s.curve = P.line([0, 0.01], [1, 1], {color: C.wave, width: 6});
      s.cur = K.g(P.g, {hidden: true});
      s.curLine = K.line(s.cur, 0, 190, 0, 710, {stroke: C.ink2, "stroke-width": 2.5, "stroke-dasharray": "8 8"});
      s.curDot = K.circle(s.cur, 0, 0, 11, {fill: C.ink, stroke: "#0a111d", "stroke-width": 3});
      s.curTxt = halo(K.text(s.cur, 0, 690, "", {cls: "t-label", size: 28, color: "var(--ink)"}));
      // annotations at the data points they name (SOIL, SAF001W_023O.TFU)
      const a0 = interp(OUT.f, OUT.amp, 0.6);
      s.one = K.g(P.g, {hidden: true});
      halo(K.text(s.one, P.X(0.6) + 24, P.Y(a0) + 46, `${a0.toFixed(0)} × at low frequency`, {cls: "t-label", size: 30, color: "var(--ink)"}));
      s.pk1 = peak(OUT.f, OUT.amp, 5, 10); s.pk2 = peak(OUT.f, OUT.amp, 14, 20);
      s.pk = K.g(P.g, {hidden: true});
      halo(K.text(s.pk, P.X(s.pk1.f), P.Y(s.pk1.a) - 30, `${s.pk1.a.toFixed(1)} × at ${s.pk1.f.toFixed(1)} Hz`, {cls: "t-big", size: 44, anchor: "middle", color: "var(--wave)"}));
      s.m1 = halo(K.text(P.g, P.X(s.pk1.f), P.Y(s.pk1.a) - 84, "1st mode", {cls: "t-label", size: 30, anchor: "middle", color: "var(--ink)", hidden: true}));
      s.m2 = K.g(P.g, {hidden: true});
      K.circle(s.m2, P.X(s.pk2.f), P.Y(s.pk2.a), 10, {fill: C.wave, stroke: "#0a111d", "stroke-width": 3});
      halo(K.text(s.m2, P.X(s.pk2.f), P.Y(s.pk2.a) - 30, `2nd mode · ${s.pk2.f.toFixed(1)} Hz`, {cls: "t-label", size: 30, anchor: "middle", color: "var(--ink)"}));
      s.note = K.pill(s.root, "one frequency at a time", {x: 1460, y: 120, color: "var(--wave)"});
      s.why = K.pill(s.root, "Things that vibrate near 7 Hz get an extra push", {x: 760, y: 850, color: "var(--ssi)", size: 32});
      s.A = 0; s.f = 0.1; s.fmax = 0; s.ph = 0;
    },
    tick(s, t, dt) {
      const f = s.f, a = interp(OUT.f, OUT.amp, f), P = s.P;
      s.ph += TAU * f / SLOW * (dt || 0);
      shake(s.col, profile(f), AP, AA * a * s.A, s.ph);
      const n = OUT.f.findIndex((x) => x > s.fmax);
      const m = n < 0 ? OUT.f.length : Math.max(2, n);
      s.curve.setAttribute("d", K.d(P.pts(OUT.f.slice(0, m), OUT.amp.slice(0, m))));
      const x = P.X(f);
      s.curLine.setAttribute("x1", x); s.curLine.setAttribute("x2", x); s.curLine.setAttribute("y1", P.Y(a) + 14);
      s.curDot.setAttribute("cx", x); s.curDot.setAttribute("cy", P.Y(a));
      s.curTxt.setAttribute("x", f > 18 ? x - 12 : x + 12);
      s.curTxt.setAttribute("text-anchor", f > 18 ? "end" : "start");
      s.curTxt.textContent = f.toFixed(1) + " Hz";
    },
    beats: [
      {say: "So, weaker or stronger? [p]SASSI answers one frequency at a time, like the notes of a chord.",
        go(k) { k.show(k.s.head); k.tween(k.s, {A: 1}, 600); }, p(k) { k.show([k.s.P.g, k.s.note]); }},
      {say: "[c]This curve shows how much the surface shakes compared with bare rock, at each frequency.",
        c(k) { k.show(k.s.cur); k.show(k.s.one, {delay: 200}); k.tween(k.s, {fmax: 25, f: 24.9}, 2400, {ease: "linear"}); }},
      {say: "[pk]But at 7.4 hertz, the surface shakes 2.2 times harder than the rock.",
        pk(k) { k.tween(k.s, {f: k.s.pk1.f}, 1000, {done: () => k.show(k.s.pk)}); }},
      {say: "[n]That's the soil column's first natural frequency; [m]its second mode, at 17.8 hertz, makes the second hump.",
        n(k) { k.show(k.s.m1); }, m(k) { k.show(k.s.m2); k.pulse(k.s.m2, {amp: 0.1}); }},
      {say: "[w]So anything in your building that vibrates near 7 hertz gets an extra push from the soil.",
        w(k) { k.show(k.s.why); k.pulse(k.s.pk, {amp: 0.06}); }},
    ],
  });

  // ---------------------------------------------------------------- 5. the half-space
  // to scale, 8 px per ft; 7.3 Hz (the peak of the rigid-base curve, within 0.2 % of the outcrop peak), the same
  // 0.04 in of rock motion in both cases (the rigid base, or the rock outcrop), displacements x 150 = 4 px
  const HSP = 8, HSA = disp(0.04, 150, HSP), HSF = 7.3;
  scenes.push({
    id: "half-space", title: "Rock that goes on forever",
    build(s) {
      const g = K.g(s.svg);
      s.U = profile(HSF);
      const rr = PH.columnRatios(COL, HSF);
      s.rOut = mag(rr.outcrop); s.rWit = mag(rr.within);       // 2.16 and 17.06: surface per unit rock motion
      s.head = K.heading(s.root, "What's below the rock?", {x: 120, y: 70});
      s.col = column(g, {x: 150, y: 230, w: 380, pxf: HSP, hsH: 260});        // base 670, rock to 930
      s.rockG = s.col.layers[2].g;
      s.rockG.classList.add("sv-in", "sv-fade", "sv-on");
      s.rigid = K.g(g, {hidden: true});
      K.ground(s.rigid, 140, s.col.base, 400, {color: "var(--ink2)"});
      K.text(s.rigid, 340, 750, "rigid base", {cls: "t-label", anchor: "middle", size: 32, color: "var(--ink2)"});
      s.trap = K.g(g, {hidden: true});
      K.arrow(s.trap, 625, 640, 625, 290, {color: C.wave, width: 5, head: 18});
      K.arrow(s.trap, 660, 290, 660, 640, {color: C.ink2, width: 5, head: 18});
      K.text(s.trap, 682, 476, "trapped", {cls: "t-label", size: 30, color: "var(--ink)"});
      s.leave = K.g(g, {hidden: true});
      K.arrow(s.leave, 600, 560, 600, 900, {color: C.dash, width: 5, head: 18});
      K.text(s.leave, 622, 800, "leaves", {cls: "t-label", size: 30, color: "var(--dash)"});
      // energy flowing down into the rock (an illustration of the flow, not of a wave)
      s.chev = [0, 1, 2, 3].map(() => K.path(g, "", {stroke: C.dash, "stroke-width": 5, fill: "none", "stroke-linecap": "round", "stroke-linejoin": "round"}));
      s.chevLab = note(g, 150, 965, "energy flows down (schematic)", {hidden: true, color: "var(--dash)"});
      s.hsPill = K.pill(s.root, "half-space", {x: 190, y: 800, color: "var(--wave)", size: 34});
      note(g, 150, 190, `${HSF} Hz · the same rock motion in both, 0.04 in`);
      note(g, 150, 218, "drawn × 150 · 10 × slower");
      const P = K.plot(s.svg, {x: 900, y: 190, w: 880, h: 460, xr: [0, 25], yr: [0, 18], xticks: [0, 5, 10, 15, 20, 25], yticks: [0, 5, 10, 15],
        xlabel: "frequency (Hz)", ylabel: "amplification", ylabelOffset: 64});
      s.P = P;
      s.wit = P.line(WIT.f, WIT.amp, {color: C.ref, width: 5, draw: true});
      s.out = P.line(OUT.f, OUT.amp, {color: C.wave, width: 6, draw: true});
      s.legW = P.legend([{label: "rigid base: surface ÷ base", color: C.ref}], {x: 1290, y: 222, hidden: true});
      s.legO = P.legend([{label: "half-space: surface ÷ rock outcrop", color: C.wave}], {x: 1290, y: 264, hidden: true});
      const pW = peak(WIT.f, WIT.amp, 5, 10), pO = peak(OUT.f, OUT.amp, 5, 10);
      s.dW = P.dot(pW.f, pW.a, {color: C.ref, hidden: true});
      s.dO = P.dot(pO.f, pO.a, {color: C.wave, hidden: true});
      s.s17 = K.stat(s.root, {x: 900, y: 750, w: 380, value: pW.a.toFixed(0), unit: "×", label: "rigid base", color: "var(--ref)", vsize: 104});
      s.s22 = K.stat(s.root, {x: 1360, y: 750, w: 380, value: pO.a.toFixed(1), unit: "×", label: "real rock: a half-space", color: "var(--wave)", vsize: 104});
      s.A = 0; s.R = 0; s.Lv = 0;
    },
    tick(s, t) {
      const ph = TAU * HSF * t / SLOW, ratio = (1 - s.R) * s.rOut + s.R * s.rWit;
      const u = shake(s.col, s.U, HSP, HSA * ratio * s.A, ph);
      s.rigid.setAttribute("transform", `translate(${u(55 * HSP).toFixed(2)},0)`);   // the rigid base moves with the input
      s.chev.forEach((c, i) => {
        const p = (t * 0.55 + i / 4) % 1, y = 700 + p * 200, x = 450;
        c.setAttribute("d", `M${x - 16},${y}L${x},${y + 16}L${x + 16},${y}`);
        c.style.opacity = String(s.Lv * Math.min(1, 4 * p) * (1 - p) * 0.95);
      });
    },
    beats: [
      {say: "Why does the boost stop at 2.2? [r]The answer lies in the rock below.",
        go(k) { k.show([k.s.head, k.s.P.g, k.s.legO]); k.draw(k.s.out, 900); k.tween(k.s, {A: 1}, 800); },
        r(k) { k.pulse(k.s.rockG, {amp: 0.04}); }},
      {say: "[box]Imagine the rock were perfectly rigid, like the bottom of a closed box.",
        box(k) { k.hide(k.s.rockG); k.show(k.s.rigid); k.dim(k.s.out, true); }},
      {say: "[b]Waves would bounce up and down, trapped, with nowhere to go.",
        b(k) { k.show(k.s.trap); k.tween(k.s, {R: 1}, 1200); }},
      {say: "With only the soil's own damping to stop them, [n]the surface would shake 17 times harder than the base.",
        n(k) { k.draw(k.s.wit, 1000); k.show([k.s.legW, k.s.dW, k.s.s17.el], {delay: 400}); }},
      {say: "[h]But real rock goes down a long way: waves that head down keep going, and never come back.",
        h(k) { k.hide([k.s.rigid, k.s.trap]); k.show([k.s.rockG, k.s.leave]); k.tween(k.s, {R: 0}, 1200); }},
      {say: "[r]Like ripples from a stone in a pond, they carry energy away, [p]and the peak stays at 2.2.",
        r(k) { k.tween(k.s, {Lv: 1}, 600); k.show(k.s.chevLab); },
        p(k) { k.dim(k.s.out, false); k.dim(k.s.wit, true); k.show([k.s.dO, k.s.s22.el]); }},
      {say: "[t]Engineers call this bottomless rock a half-space.",
        t(k) { k.show(k.s.hsPill); }},
    ],
  });

  // ---------------------------------------------------------------- 6. how SASSI does it: SITE
  // to scale, 6.8 px per ft: the 22 sublayers, then the generated rock sublayers as at 20 Hz (20 of 12.38 ft,
  // 248 ft in all, SITE listing: the first two and the last are drawn); the column moves at 20 Hz,
  // 0.04 in x 1000 = 23 px
  const SP = 6.8, SA = disp(0.04, 1000, SP), ST = 230, SB = ST + 55 * SP;   // surface 230, top of the rock 604
  const HL = 12.378, HB = 55 + 20 * HL;                                     // rock sublayer at 20 Hz; their base (ft)
  const BUF = [[55, 55 + HL, SB, SB + HL * SP], [55 + HL, 55 + 2 * HL, SB + HL * SP, SB + 2 * HL * SP], [HB - HL, HB, 800, 800 + HL * SP]];
  const BOT = 800 + HL * SP;                                      // 884: the bottom of the 20th rock sublayer
  scenes.push({
    id: "site-module", title: "How SASSI does it: SITE",
    build(s) {
      const g = K.g(s.svg);
      s.U = profile(20);
      s.head = K.heading(s.root, "How SASSI does it: <span style='color:var(--wave)'>SITE</span>", {x: 120, y: 70});
      s.col = column(g, {x: 180, y: ST, w: 400, pxf: SP, hsH: 0, hidden: true});
      s.subs = K.path(g, IFC.slice(1, 22).map((z) => `M180,${ST + SP * z}L580,${ST + SP * z}`).join(""), {stroke: "rgba(255,244,225,.85)", "stroke-width": 2.5, fill: "none", draw: true});
      s.subLab = K.text(g, 640, 420, "22 thin slices", {cls: "t-label", size: 32, color: "var(--ink)", hidden: true});
      // the generated half-space sublayers (to scale at 20 Hz) and the viscous boundary at their base
      K.defs();
      s.buf = K.g(g, {hidden: true});
      s.bufP = BUF.map(() => [K.path(s.buf, "", {fill: C.rock, "fill-opacity": 0.6}), K.path(s.buf, "", {fill: "url(#sv-pat-rock)"}),
        K.path(s.buf, "", {stroke: "rgba(255,255,255,.55)", "stroke-width": 2, fill: "none", "stroke-dasharray": "8 6"})]);
      K.text(s.buf, 380, 792, "⋮", {cls: "t-label", anchor: "middle", size: 34, color: "var(--ink2)"});
      K.text(s.buf, 640, 690, "20 rock slices,", {cls: "t-label", size: 30, color: "var(--ink)"});
      K.text(s.buf, 640, 726, "248 ft deep at 20 Hz", {cls: "t-label", size: 30, color: "var(--ink)"});
      K.text(s.buf, 640, 760, "(17 not drawn)", {cls: "t-small"});
      s.dps = K.g(g, {hidden: true});
      s.dpV = K.dashpot(s.dps, 300, BOT, 300, 950, {w: 26, cyl: 36});
      K.ground(s.dps, 260, 950, 80, {color: "var(--muted)"});
      s.dpLink = K.line(s.dps, 470, BOT, 470, 918, {stroke: C.dash, "stroke-width": 4});
      s.dpH = K.dashpot(s.dps, 470, 918, 565, 918, {w: 26, cyl: 36});
      K.line(s.dps, 568, 890, 568, 946, {stroke: C.muted, "stroke-width": 5, "stroke-linecap": "round"});
      for (let y = 894; y < 946; y += 14) K.line(s.dps, 572, y, 588, y + 14, {cls: "hatch"});
      K.text(s.dps, 276, 930, "ρVp", {cls: "t-small", anchor: "end", color: "var(--dash)"});
      K.text(s.dps, 520, 965, "ρVs", {cls: "t-small", anchor: "middle", color: "var(--dash)"});
      K.text(s.dps, 640, 905, "shock absorbers", {cls: "t-label", size: 30, color: "var(--dash)"});
      K.text(s.dps, 640, 937, "dashpots, schematic", {cls: "t-small"});
      s.ifc = K.nodes(g, IFC.map((z) => [580, ST + SP * z]), {r: 6, hidden: true});
      s.cp = K.circle(g, 580, ST, 14, {fill: C.wave, stroke: C.ink, "stroke-width": 3, hidden: true, in: "pop"});
      s.cpPill = K.pill(s.root, "control point", {x: 625, y: 196, color: "var(--wave)", size: 30});
      s.scale = note(g, 640, 985, "20 Hz · surface 0.04 in, drawn × 1000 · 30 × slower", {hidden: true});
      s.code = K.code(s.root, ["L,1,1.6,0.120,2000,1000,0.05,0.05", "TOPL,1,1,1,1,1,1,1,1,1,1"], {x: 1000, y: 230, w: 800, size: 34});
      s.codeL = [K.label(s.root, "L: the sand · Vs 1,000 ft/s · 5 % damping", {x: 1020, y: 470, size: 30}),
        K.label(s.root, "TOPL: ten 1.6 ft slices of it, from the top", {x: 1020, y: 530, size: 30})];
      // the check against the exact solution (SOIL, the lesson's SHAKE run) at four SITE frequencies (listing)
      const P = K.plot(s.svg, {x: 1080, y: 230, w: 700, h: 380, xr: [0, 25], yr: [0, 2.5], xticks: [0, 5, 10, 15, 20, 25], yticks: [0, 1, 2],
        xlabel: "frequency (Hz)", ylabel: "surface ÷ bare rock", ylabelOffset: 70});
      s.P = P;
      s.exact = P.line(OUT.f, OUT.amp, {color: C.wave, width: 5});
      s.dots = [[F01, 1.000], [287 * DF, 2.126], [492 * DF, 1.480], [737 * DF, 2.181]].map(([f, a]) => P.dot(f, a, {color: C.ink, r: 11, hidden: true}));
      s.leg = P.legend([{label: "exact solution", color: C.wave}, {label: "SITE", color: C.ink, dash: "0.1 14"}], {x: 1480, y: 560});
      s.stat = K.stat(s.root, {x: 1080, y: 760, w: 700, value: "0.5", unit: "%", label: "SITE and the exact solution agree within this", color: "var(--good)", vsize: 104});
      s.A = 0;
    },
    tick(s, t) {
      const ph = TAU * 20 * t / 30, A = SA * s.A;
      const u = shake(s.col, s.U, SP, A, ph);
      s.subs.setAttribute("d", IFC.slice(1, 22).map((z) => { const du = u(z * SP), y = ST + SP * z; return `M${180 + du},${y}L${580 + du},${y}`; }).join(""));
      s.ifc.dots.forEach((dot, i) => dot.setAttribute("cx", 580 + u(IFC[i] * SP)));
      s.cp.setAttribute("cx", 580 + u(0));
      // the rock sublayers move with the exact solution at their real depth
      const ub = (z) => A * re(PH.columnU(s.U.wv, COL, z), ph);
      BUF.forEach(([z0, z1, y0, y1], i) => {
        const Lp = [], Rp = [];
        for (let j = 0; j <= 4; j++) { const z = z0 + (z1 - z0) * j / 4, y = y0 + (y1 - y0) * j / 4, du = ub(z); Lp.push([180 + du, y]); Rp.push([580 + du, y]); }
        const d = K.d(Lp.concat(Rp.reverse())) + "Z";
        s.bufP[i][0].setAttribute("d", d); s.bufP[i][1].setAttribute("d", d);
        const a = ub(z0), b = ub(z1);
        s.bufP[i][2].setAttribute("d", `M${180 + a},${y0}L${580 + a},${y0}M${180 + b},${y1}L${580 + b},${y1}`);
      });
      const uB = ub(HB);
      s.dpV.update(300 + uB, BOT, 300 + uB, 950);                       // on rollers: it moves sideways with the base
      s.dpLink.setAttribute("x1", 470 + uB); s.dpLink.setAttribute("x2", 470 + uB);
      s.dpH.update(470 + uB, 918, 565, 918);
    },
    beats: [
      {say: "Inside SASSI, one module computes this free field: [s]SITE.",
        go(k) { k.show(k.s.head); }, s(k) { k.pulse(k.s.head, {amp: 0.04}); k.show(k.s.col.g); }},
      {say: "[a]It cuts the soil into thin slices, like a finite-element mesh.",
        a(k) { k.draw(k.s.subs, 1000); k.show(k.s.subLab, {delay: 400}); }},
      {say: "[c]In the input, L describes the sand, and TOPL stacks ten thin slices of it.",
        c(k) { k.show(k.s.code.el); k.type(k.s.code.lines, {cps: 40}); k.show(k.s.codeL, {stagger: 150, delay: 400}); }},
      {say: "[h]Below the soil, SITE adds extra rock slices, [d]with shock absorbers at the bottom, so waves can leave.",
        h(k) { k.show(k.s.buf); }, d(k) { k.show(k.s.dps); }},
      {say: "[cp]Your design motion goes in at one point, the control point: usually the ground surface.",
        cp(k) { k.show([k.s.cp, k.s.cpPill]); k.pulse(k.s.cp, {amp: 0.5}); }},
      {say: "[f]SITE then works out how every slice moves, at every frequency, for one unit of motion there.",
        f(k) { k.show([k.s.ifc.g, k.s.scale]); k.tween(k.s, {A: 1}, 1000); }},
      {say: "Can you trust it? [ok]Against an exact solution of the same column, SITE agrees within half a percent.",
        go(k) { k.hide([k.s.code.el].concat(k.s.codeL)); },
        ok(k) { k.show(k.s.P.g); k.show(k.s.dots, {stagger: 150, delay: 300}); k.show(k.s.stat.el, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 7. what it means: the motion changes with depth
  // to scale, 11 px per ft; surface motion 0.04 in drawn x 1000 = 37 px
  const DP = 11, DA = disp(0.04, 1000, DP);
  scenes.push({
    id: "depth", title: "Every layer moves differently",
    build(s) {
      const g = K.g(s.svg);
      s.lo = profile(F01); s.hi = profile(F20);
      s.head = K.heading(s.root, "Every layer moves differently", {x: 120, y: 70});
      s.col = column(g, {x: 200, y: 220, w: 400, pxf: DP, hsH: 120});         // base 825, rock to 945
      note(g, 200, 990, "surface 0.04 in, drawn × 1000 · time not to scale");
      // a 32 ft x 16 ft basement (the size of example 2), outline only
      s.bsmt = K.g(g, {hidden: true});
      K.rect(s.bsmt, 400 - 16 * DP, 220, 32 * DP, 16 * DP, {rx: 4, fill: "rgba(10,17,29,.35)", stroke: C.ssi, "stroke-width": 4, "stroke-dasharray": "12 8"});
      halo(K.text(s.bsmt, 400, 296, "basement", {cls: "t-label", anchor: "middle", size: 30, color: "var(--ssi)"}));
      halo(K.text(s.bsmt, 400, 330, "16 ft deep · outline", {cls: "t-small", anchor: "middle", color: "var(--ink2)"}));
      const P = K.plot(s.svg, {x: 860, y: 220, w: 480, h: 55 * DP, xr: [0, 1.1], yr: [55, 0], xticks: [0, 0.5, 1], yticks: [0, 16, 55],
        yfmt: (v) => v + " ft", xlabel: "motion (surface = 1)", ylabel: "depth", ylabelOffset: 86});
      s.P = P;
      const zs = K.linspace(0, 55, 111);
      s.lLo = P.line(zs.map((z) => mag(s.lo(z))), zs, {color: C.ref, width: 5, draw: true});
      s.lHi = P.line(zs.map((z) => mag(s.hi(z))), zs, {color: C.wave, width: 6, draw: true});
      s.legLo = P.legend([{label: "0.098 Hz", color: C.ref}], {x: 1410, y: 520, size: 30, hidden: true});
      s.legHi = P.legend([{label: "20 Hz", color: C.wave}], {x: 1410, y: 570, size: 30, hidden: true});
      // annotations at the computed points they name (Z13: interface 9 of SITE, 12.8 ft)
      const Z13 = 12.8, u0 = mag(s.hi(0)), u4 = mag(s.hi(Z13));
      s.cS = K.g(P.g, {hidden: true});
      K.line(s.cS, P.X(u0) + 12, P.Y(0), 1400, P.Y(0) + 14, {stroke: C.ink2, "stroke-width": 2});
      K.text(s.cS, 1410, P.Y(0) + 24, "surface: moves fully", {cls: "t-label", size: 30, color: "var(--ink)"});
      s.c4 = K.g(P.g, {hidden: true});
      K.circle(s.c4, P.X(u4), P.Y(Z13), 10, {fill: C.wave, stroke: "#0a111d", "stroke-width": 3});
      K.line(s.c4, P.X(u4) + 14, P.Y(Z13), 1400, P.Y(Z13), {stroke: C.ink2, "stroke-width": 2});
      K.text(s.c4, 1410, P.Y(Z13) - 6, "13 ft down:", {cls: "t-label", size: 30, color: "var(--ink)"});
      K.text(s.c4, 1410, P.Y(Z13) + 32, `almost still (${u4.toFixed(2)})`, {cls: "t-label", size: 30, color: "var(--wave)"});
      s.kin = K.card(s.root, {x: 1400, y: 670, w: 430, title: "Kinematic interaction", size: 30,
        body: "a stiff basement averages these motions · lesson 5"});
      s.kin.querySelector(".sv-card-t").style.color = "var(--kin)";
      s.A = 0; s.mix = 0;
    },
    tick(s, t) {
      const ph = TAU * t / 1.6, A = DA * s.A, m = s.mix;
      s.col.shear((d) => { const z = d / DP; return A * ((1 - m) * re(s.lo(z), ph) + m * re(s.hi(z), ph)); });
    },
    beats: [
      {say: "[a]At very low frequency, the whole column moves together, as one block.",
        go(k) { k.show([k.s.head, k.s.P.g]); }, a(k) { k.tween(k.s, {A: 1}, 800); k.draw(k.s.lLo, 1000); k.show(k.s.legLo, {delay: 400}); }},
      {say: "[b]But at 20 hertz, the surface moves fully, [c]while 13 ft down, the soil is almost still.",
        b(k) { k.tween(k.s, {mix: 1}, 1000); k.dim(k.s.lLo, true); k.draw(k.s.lHi, 1000); k.show([k.s.legHi, k.s.cS], {delay: 400}); },
        c(k) { k.show(k.s.c4); }},
      {say: "[k]A stiff basement averages them: that's kinematic interaction, the topic of lesson five.",
        k(k) { k.show(k.s.bsmt); k.show(k.s.kin, {delay: 300}); }},
    ],
  });

  // ---------------------------------------------------------------- 8. recap
  scenes.push({
    id: "recap", title: "Recap",
    build(s) {
      s.head = K.heading(s.root, "Recap", {x: 120, y: 80, size: "h1"});
      s.list = K.bullets(s.root, [
        {t: "Soil layers amplify their own frequencies", sub: "2.2 × the rock at 7.4 Hz"},
        {t: "SITE computes the free field", sub: "every layer, every frequency"},
        {t: "The half-space lets waves leave", sub: "peak 2.2 ×, not 17 × on a rigid base"},
      ], {x: 120, y: 230, w: 900, num: true});
      s.q = K.card(s.root, {x: 1100, y: 230, w: 720, kind: "check", title: "Check yourself",
        body: "Same soil, on much stiffer rock. Higher or lower surface peak?"});
      s.a = K.card(s.root, {x: 1100, y: 540, w: 720, title: "Answer",
        body: "<b>Higher</b>: less energy escapes, so the peak climbs towards the rigid-base value."});
      s.next = K.pill(s.root, "Next · Lesson 3: soil springs and dashpots →", {x: 120, y: 880, color: "var(--wave)", size: 32});
    },
    beats: [
      {say: "Let's recap. [a]Layers of soil amplify the shaking at their own natural frequencies.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[a]SITE computes this free field, for every layer and every frequency.",
        a(k) { k.show(k.s.list.items[1]); }},
      {say: "[a]And the half-space below lets waves leave, which keeps the peaks down.",
        a(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Check yourself: same soil, but on much stiffer rock. Will the surface peak be higher or lower?",
        q(k) { k.show(k.s.q); }, gap: 1500},
      {say: "[a]Higher. Stiffer rock lets less energy escape, so the peak climbs towards the rigid-base value.",
        a(k) { k.show(k.s.a); }},
      {say: "[n]Next: the soil under a foundation acts like springs and shock absorbers. That's lesson three.",
        n(k) { k.show(k.s.next); }},
    ],
  });

  SV.video({
    id: "02", n: 2, part: "Fundamentals", lesson: "02-free-field",
    title: "The site: layers, half-space and free-field motion",
    subtitle: "How layers of soil amplify the earthquake, how SITE computes the free field, and why the rock below lets waves leave.",
    next: {href: "03.html", title: "Foundation impedance: soil springs and dashpots"},
    pron: [[/\bSITE\b/g, "site"], [/\bTOPL\b/g, "top L"]],
    scenes,
  });
})();
