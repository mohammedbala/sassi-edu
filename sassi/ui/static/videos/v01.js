/* Explainer video, lesson 1: From fixed base to SSI -- what changes and why (newcomer edition).
 * Numbers: sassi/ui/lessons/01_why_ssi.md; curves: d01.js (the lesson's runs, python -m sassi.ui.video_data). */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["01"];
  const TAU = 2 * Math.PI;
  const DEG = Math.PI / 180;

  /** Vibration of the example stick, total motion per unit ground motion (the lesson's runs, .TFU files):
   *  fixed base at 5.005 Hz: mat 1.00 in phase, roof 12.7 at -107 deg; SSI at 3.49 Hz: mat 1.42 at -54 deg,
   *  roof 13.0 at -94 deg, rocking +-1.58 at the mat edges (32 ft from the centre) at -98 deg.
   *  px: pixels per unit of ground motion; period: visual period (s). */
  function vibrate(st, t, o) {
    const ph = TAU * t / o.period;
    const ug = Math.cos(ph) * o.px;
    const mat = o.ssi ? 1.42 * Math.cos(ph - 54 * DEG) * o.px : ug;
    const roof = (o.ssi ? 13.0 * Math.cos(ph - 94 * DEG) : 12.7 * Math.cos(ph - 107 * DEG)) * o.px;
    const rock = o.ssi ? (1.58 * o.px / o.halfW) * Math.cos(ph - 98 * DEG) * (o.rockGain || 1) : 0;
    const H = o.zs[o.zs.length - 1];
    const bendRoof = roof - mat - Math.sin(rock) * H;
    const shape = (z) => Math.sin((Math.PI / 2) * z / H);
    st.set({sway: mat * o.amp, rock: rock * o.amp, bend: (z) => bendRoof * shape(z) * o.amp});
    return ug * o.amp;
  }
  const big = (root, text, o) => K.h("div", Object.assign({cls: "sv-h2", html: text}, o), root);
  const PH = K.PH, CX = PH.C;
  /** Free field of the lesson-1 site (16 ft sand, 39 ft gravel, rock) at f: u(z ft, phase) per unit surface motion. */
  function freeField(f) {
    const wv = PH.columnWaves(PH.COLUMN, f), cache = new Map();
    const U = (z) => { const k = Math.round(z * 4); if (!cache.has(k)) cache.set(k, PH.columnU(wv, PH.COLUMN, Math.max(0, Math.min(54.99, z)))); return cache.get(k); };
    return {U, at: (z, ph) => { const u = U(z); return u[0] * Math.cos(ph) - u[1] * Math.sin(ph); }};
  }
  /** Two decimals of a spectral value stored with four digits: a final 5 (8.035, 1.395) is rounded down, as the
   *  lesson reads the full-precision files (8.03 g, 1.39 g). */
  const g2 = (v) => (v - 5e-5).toFixed(2);
  /** Largest value of a curve between f1 and f2: {f, v}. */
  function peak(xs, ys, f1, f2) {
    let b = {f: xs[0], v: -Infinity};
    xs.forEach((x, i) => { if (x >= (f1 || -Infinity) && x <= (f2 || Infinity) && ys[i] > b.v) b = {f: x, v: ys[i]}; });
    return b;
  }

  const scenes = [];

  // ---------------------------------------------------------------- 0. the model you know
  scenes.push({
    id: "fixed-base", title: "The model you know",
    build(s) {
      const g = K.g(s.svg);
      // to scale, 5 px per ft: 64 ft mat 5 ft thick, 16 ft storeys, 16 ft sand on 39 ft gravel
      s.zs = [80, 160, 240, 320];
      s.ground = K.ground(g, 380, 644, 360, {hidden: true});
      s.soil = K.soil(g, {x: 300, y: 644, w: 520, layers: [{h: 80, kind: "sand"}, {h: 195, kind: "gravel"}], hs: {h: 77, kind: "rock"}, labels: false, hidden: true});
      s.stick = K.stick(g, {x: 560, y: 619, z: s.zs, matW: 320, matH: 25, r: 18, slabW: 110, hidden: true});
      s.ff = freeField(3.49);
      s.note = K.text(g, 560, 1005, "computed motion per unit ground motion · slowed down 10×", {cls: "t-small", anchor: "middle", hidden: true});
      s.acc = K.signal(g, {x: 330, y: 700, w: 460, h: 80, fn: (u) => {
        const a = D["data/rg160h_030g.acc"]; const i = Math.min(a.v.length - 1, Math.floor(u * a.v.length * 0.55));
        return a.v[i] / a.peak; }, color: C.wave, width: 3, draw: true, n: 260});
      s.accLab = K.text(g, 560, 820, "the design earthquake, at the base", {cls: "t-small", anchor: "middle", hidden: true});
      s.head = K.heading(s.root, "The fixed-base model", {x: 1020, y: 210, size: "h1"});
      s.sub = K.para(s.root, "A building clamped to rigid ground.", {x: 1020, y: 310, w: 800, size: 38});
      s.outs = ["Floor spectra", "Forces", "Displacements"].map((t, i) =>
        K.pill(s.root, t, {x: 1020 + [0, 250, 400][i], y: 420, color: "var(--ink2)"}));
      s.q = big(s.root, 'What changes when<br>the ground is <span style="color:var(--ssi)">soft</span>?', {x: 1020, y: 560, w: 820, size: 64});
      s.A = 0; s.ssi = 0;
    },
    tick(s, t) {
      const ssi = s.ssi > 0, period = ssi ? 10 / 3.49 : 10 / 5.0;           // slowed down 10 times
      vibrate(s.stick, t, {period, px: 3.4, halfW: 160, zs: s.zs, amp: s.A, ssi});
      if (ssi) { const ph = TAU * t / period; s.soil.shear((d) => 3.4 * s.A * s.ff.at(d / 5, ph)); }
    },
    beats: [
      {say: "Picture the seismic model you already know: a building clamped to perfectly rigid ground.",
        go(k) { const s = k.s; k.show([s.ground, s.stick.g, s.note]); k.show([s.head, s.sub], {stagger: 150}); k.tween(s, {A: 1}, 1000); }},
      {say: "You shake its base with the [acc]design earthquake, and read off [o]floor spectra, forces and displacements.",
        acc(k) { k.draw(k.s.acc, 1500); k.show(k.s.accLab, {delay: 400}); }, o(k) { k.show(k.s.outs, {stagger: 220}); }},
      {say: "But real buildings don't stand on rigid ground. [soil]They stand on soil, and soil gives.",
        soil(k) { const s = k.s; k.hide([s.ground, s.acc, s.accLab]); k.show(s.soil.g); k.tween(s, {ssi: 1}, 10); k.dim([s.sub, ...s.outs], true); }},
      {say: "[q]So what changes when the ground can move with the building?", q(k) { k.show(k.s.q); }, hold: 600},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 1",
    build(s) {
      s.t = K.titleCard(s.root, {n: 1, part: "Fundamentals", title: "From fixed base to SSI", sub: "What changes when a building stands on soil, and why."});
    },
    beats: [
      {say: "Lesson one: from fixed base to soil-structure interaction, or SSI.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 150}); k.show(t.title, {delay: 300}); k.show(t.sub, {delay: 600}); k.show(t.rule, {delay: 800}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. the soil is a suspension
  scenes.push({
    id: "suspension", title: "Soil acts like a suspension",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "The soil acts like a <span style='color:var(--spring)'>suspension</span>", {x: 120, y: 90, size: "h2"});
      s.zs = [95, 190, 285, 380];                      // 5.94 px per ft: 16 ft storeys, like the 64 ft mat (380 px)
      s.soil = K.soil(g, {x: 150, y: 726, w: 960, layers: [{h: 250, kind: "gravel"}], labels: false});
      s.stick = K.stick(g, {x: 630, y: 700, z: s.zs, matW: 380, matH: 28, r: 22, slabW: 140});
      s.sp = K.g(g, {hidden: true});
      s.springs = [500, 760].map((x) => K.spring(s.sp, x, 728, x, 880, {coils: 6, amp: 17}));
      s.hspring = K.spring(s.sp, 250, 700, 440, 700, {coils: 6, amp: 15});
      K.ground(s.sp, 420, 880, 420, {color: "var(--muted)"});
      K.line(s.sp, 250, 640, 250, 790, {stroke: "var(--muted)", "stroke-width": 6});
      s.dp = K.g(g, {hidden: true});
      s.dashes = [590, 670].map((x) => K.dashpot(s.dp, x, 728, x, 880, {w: 30, cyl: 80}));
      s.hdash = K.dashpot(s.dp, 440, 760, 250, 760, {w: 28, cyl: 80});
      s.ripples = K.ripples(g, {cx: 630, cy: 890, r0: 60, r1: 480, n: 4, color: C.dash, squash: 0.5, hidden: true});
      s.forces = K.g(g, {hidden: true});
      s.farrows = [0, 1, 2, 3].map(() => K.arrow(s.forces, 0, 0, 1, 0, {color: C.ssi, width: 6, head: 20}));
      s.sch = K.text(g, 630, 1000, "springs and dashpots: schematic · motion computed at 3.49 Hz, slowed down 10×", {cls: "t-small", anchor: "middle"});
      s.l1 = K.card(s.root, {x: 1200, y: 250, w: 620, title: "Springs", body: "The building slides and rocks: it sways <b>more slowly</b>.", size: 32});
      s.l1.querySelector(".sv-card-t").style.color = "var(--spring)";
      s.l2 = K.card(s.root, {x: 1200, y: 480, w: 620, title: "Shock absorbers", body: "Energy leaves as waves, <b>like ripples in a pond</b>.", size: 32});
      s.l2.querySelector(".sv-card-t").style.color = "var(--dash)";
      s.tag = K.pill(s.root, "Inertial interaction", {x: 1200, y: 740, color: "var(--ssi)"});
      s.tagSub = K.para(s.root, "the building's sway pushes back on the soil", {x: 1200, y: 810, w: 620, size: 30});
      s.A = 0; s.F = 0; s.Rp = 0;
    },
    tick(s, t) {
      vibrate(s.stick, t, {period: 10 / 3.49, px: 3.6, halfW: 190, zs: s.zs, amp: s.A, ssi: true});
      const st = s.stick.state;
      s.springs.forEach((sp, i) => { const x = [500, 760][i]; sp.update(x + st.sway, 728 + (x - 630) * Math.sin(st.rock), x, 880); });
      s.dashes.forEach((dp, i) => { const x = [590, 670][i]; dp.update(x + st.sway, 728 + (x - 630) * Math.sin(st.rock), x, 880); });
      s.hspring.update(250, 700, 440 + st.sway, 700);
      s.hdash.update(440 + st.sway, 760, 250, 760);
      // inertia force of a floor = m * omega^2 * u: along its total displacement u, in proportion to it (equal masses)
      const nodes = s.stick.nodes();
      s.farrows.forEach((a, i) => {
        const q = nodes[i + 1], f = (q[0] - 630) * 2.2 * s.F;
        const sg = Math.sign(f || 1);
        a.update(q[0] + sg * 26, q[1], q[0] + sg * 26 + f, q[1]);
      });
      s.ripples.amp = s.Rp;
      s.ripples.set(t, 2.0);
    },
    beats: [
      {say: "Under the foundation, the soil behaves like a car's suspension.",
        go(k) { k.show(k.s.head); k.tween(k.s, {A: 1}, 1200); }},
      {say: "It has [sp]springs, and [dp]shock absorbers. Engineers call them springs and dashpots.",
        sp(k) { k.show(k.s.sp); k.dim(k.s.soil.g, true); }, dp(k) { k.show(k.s.dp); }},
      {say: "[s]The springs let the foundation slide and rock, so the whole building sways more slowly.",
        s(k) { k.show(k.s.l1); }},
      {say: "[d]The shock absorbers drain energy away, as waves spreading into the ground, [r]like ripples in a pond.",
        d(k) { k.show(k.s.l2); }, r(k) { k.show(k.s.ripples.g); k.tween(k.s, {Rp: 1}, 1000); }},
      {say: "[t]This is called inertial interaction: the building's own swaying [f]pushes back on the soil.",
        t(k) { k.show([k.s.tag, k.s.tagSub], {stagger: 200}); }, f(k) { k.show(k.s.forces); k.tween(k.s, {F: 1}, 800); }},
    ],
  });

  // ---------------------------------------------------------------- 3. a stiff foundation averages
  const KS = 9.25;                                     // px per ft
  scenes.push({
    id: "kinematic", title: "A stiff foundation averages",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "A second effect: <span style='color:var(--kin)'>averaging</span>", {x: 120, y: 90, size: "h2"});
      // to scale, 9.25 px per ft: the lesson-1 site (16 ft sand on 39 ft gravel); its free field at 12 Hz
      s.ksoil = K.soil(g, {x: 150, y: 420, w: 960, layers: [{h: 16 * KS, kind: "sand"}, {h: 39 * KS, kind: "gravel"}], labels: false});
      s.block = K.rect(g, 640 - 16 * KS, 400, 32 * KS, 20 + 16 * KS, {fill: "rgba(95,225,208,.22)", stroke: C.kin, "stroke-width": 6, rx: 3, hidden: true});
      s.ff = K.path(g, "", {stroke: C.wave, "stroke-width": 6, fill: "none"});
      s.ffLab = K.text(g, 300, 400, "the ground", {cls: "t-label", anchor: "middle", color: "var(--wave)"});
      s.blLab = K.text(g, 640, 385, "stiff box, 32 ft wide, 16 ft deep", {cls: "t-label", anchor: "middle", color: "var(--kin)", hidden: true});
      s.kNote = K.text(g, 630, 990, "free field at 12 Hz, computed for this site · box: average over its depth (simplified)", {cls: "t-small", anchor: "middle"});
      s.kf = freeField(12);
      // the box takes the average of the free field over its 16 ft depth (complex mean)
      let re = 0, im = 0, n = 0;
      for (let z = 0; z <= 16.001; z += 0.8) { const u = s.kf.U(z); re += u[0]; im += u[1]; n++; }
      s.kavg = [re / n, im / n];
      s.c1 = K.para(s.root, "Waves make the ground move <b>differently at each depth</b>.", {x: 1200, y: 300, w: 620, size: 36});
      s.c2 = K.para(s.root, "A stiff foundation <b style='color:var(--kin)'>can't follow every wiggle</b>: it averages them.", {x: 1200, y: 480, w: 620, size: 36});
      s.tag = K.pill(s.root, "Kinematic interaction", {x: 1200, y: 700, color: "var(--kin)"});
      s.tagSub = K.para(s.root, "most important for buried basements: lesson 5", {x: 1200, y: 770, w: 620, size: 30});
      s.W = 0;
    },
    tick(s, t) {
      const ph = TAU * t / 2.4, A = 40 * s.W;                                  // 12 Hz slowed down about 30 times
      const pts = [];
      for (let d = 0; d <= 55 * KS; d += 10) pts.push([300 + A * s.kf.at(d / KS, ph), 420 + d]);
      s.ff.setAttribute("d", K.d(pts, true));
      s.ksoil.shear((d) => A * s.kf.at(d / KS, ph));
      const ub = A * (s.kavg[0] * Math.cos(ph) - s.kavg[1] * Math.sin(ph));
      s.block.setAttribute("transform", `translate(${ub} 0)`);
    },
    beats: [
      {say: "[w]Seismic waves make the ground move differently at different depths.", go(k) { k.show(k.s.head); }, w(k) { k.tween(k.s, {W: 1}, 1000); k.show(k.s.c1); }},
      {say: "[b]A stiff foundation can't follow every wiggle, so it averages them.", b(k) { k.show([k.s.block, k.s.blLab, k.s.c2]); }},
      {say: "[t]That's kinematic interaction. It matters most for buried basements, in lesson five.", t(k) { k.show([k.s.tag, k.s.tagSub], {stagger: 200}); }},
    ],
  });

  // ---------------------------------------------------------------- 4. how SASSI does it
  scenes.push({
    id: "split", title: "How SASSI does it",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "SASSI splits the problem in two", {x: 120, y: 90, size: "h2"});
      const P = (x, title, color) => {
        const pg = K.g(g, {hidden: true});
        K.rect(pg, x, 220, 460, 470, {rx: 20, fill: "rgba(20,32,51,.6)", stroke: "#2a3b52", "stroke-width": 2});
        K.text(pg, x + 230, 280, title, {cls: "t-label", anchor: "middle", color, size: 34});
        return pg;
      };
      s.p1 = P(120, "The soil layers", "var(--wave)");
      K.soil(s.p1, {x: 160, y: 470, w: 380, layers: [{h: 60, kind: "sand"}, {h: 90, kind: "gravel"}], hs: {h: 50, kind: "rock"}, labels: false});
      s.n1 = K.nodes(s.p1, K.linspace(250, 450, 7).map((x) => [x, 470]), {r: 8});
      s.p2 = P(730, "The building", "var(--steel)");
      K.stick(s.p2, {x: 960, y: 560, z: [55, 110, 165, 220], matW: 210, matH: 20, r: 14, slabW: 80});
      s.n2 = K.nodes(s.p2, K.linspace(860, 1060, 7).map((x) => [x, 580]), {r: 8});
      s.p3 = P(1340, "Joined", "var(--ssi)");
      K.soil(s.p3, {x: 1380, y: 520, w: 380, layers: [{h: 50, kind: "sand"}, {h: 70, kind: "gravel"}], hs: {h: 40, kind: "rock"}, labels: false});
      K.stick(s.p3, {x: 1570, y: 494, z: [50, 100, 150, 200], matW: 210, matH: 20, r: 13, slabW: 76});
      s.n3 = K.nodes(s.p3, K.linspace(1470, 1670, 7).map((x) => [x, 520]), {r: 8});
      s.ops = ["+", "="].map((t, i) => K.text(g, [655, 1265][i], 470, t, {cls: "t-big", anchor: "middle", size: 72, hidden: true}));
      s.c1 = K.label(s.root, "exact wave formulas", {x: 120, y: 710, w: 460, align: "center", color: "var(--muted)"});
      s.c2 = K.label(s.root, "finite elements, as in ANSYS", {x: 730, y: 710, w: 460, align: "center", color: "var(--muted)"});
      s.nodes = K.pill(s.root, "● interaction nodes: where they touch", {x: 1290, y: 780, color: "var(--wave)"});
    },
    beats: [
      {say: "How does SASSI compute all this? It splits the problem in two.", go(k) { k.show(k.s.head); }},
      {say: "[a]The soil layers are solved on their own, with exact wave formulas.", a(k) { k.show([k.s.p1, k.s.c1]); }},
      {say: "[b]The building is an ordinary finite-element model, much like your ANSYS model.", b(k) { k.show([k.s.ops[0], k.s.p2, k.s.c2]); }},
      {say: "[c]Then the two are joined at the points where they touch, [n]called interaction nodes.",
        c(k) { k.show([k.s.ops[1], k.s.p3]); }, n(k) { k.show(k.s.nodes); k.pulse([k.s.n1.g, k.s.n2.g, k.s.n3.g], {amp: 0.1}); }},
    ],
  });

  // ---------------------------------------------------------------- 5. one frequency at a time
  scenes.push({
    id: "notes", title: "One frequency at a time",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "One frequency at a time", {x: 120, y: 90, size: "h2"});
      const acc = D["data/rg160h_030g.acc"];
      const at = (u) => { const i = Math.min(acc.v.length - 1, Math.floor(u * acc.v.length * 0.6)); return acc.v[i] / acc.peak; };
      s.mix = K.signal(g, {x: 120, y: 260, w: 560, h: 220, fn: at, color: C.wave, width: 3, draw: true, n: 300});
      s.mixLab = K.text(g, 400, 530, "an earthquake: a mix of frequencies", {cls: "t-label", anchor: "middle", hidden: true});
      s.arrow = K.arrow(g, 720, 370, 820, 370, {color: C.muted, width: 5, hidden: true});
      const notes = [1, 2.5, 4, 7];
      s.notes = notes.map((f, i) => {
        const ng = K.g(g, {hidden: true});
        K.signal(ng, {x: 880, y: 200 + i * 95, w: 420, h: 70, fn: (u) => Math.sin(TAU * f * 1.4 * u), color: i === 2 ? C.ssi : C.wave, width: 4, n: 240});
        K.text(ng, 1330, 245 + i * 95, `${f} Hz`, {cls: "t-label", color: i === 2 ? "var(--ssi)" : "var(--ink2)"});
        return ng;
      });
      s.notesLab = K.text(g, 1090, 600, "its notes (illustration)", {cls: "t-small", anchor: "middle", hidden: true});
      const tf = D["ex01/00085TR_X.TFI"];
      s.p = K.plot(s.svg, {x: 1480, y: 230, w: 330, h: 300, xr: [0, 10], yr: [0, 14], xticks: [0, 5, 10], yticks: [0, 7, 14], grid: false});
      s.tf = s.p.line(tf.f, tf.amp, {color: C.ssi, width: 5, draw: true});
      s.tfLab = K.label(s.root, "roof transfer function", {x: 1430, y: 600, w: 430, align: "center", color: "var(--ssi)"});
      s.eqz = K.card(s.root, {x: 1100, y: 700, w: 700, title: "Transfer function", size: 32,
        body: "How much a point <b>amplifies</b> the ground motion, at each frequency. Like an equaliser."});
      s.back = K.card(s.root, {x: 120, y: 700, w: 900, title: "Then put the notes back together", size: 32,
        body: "with the real earthquake: floor motions and floor spectra."});
    },
    beats: [
      {say: "[g]An earthquake record is a messy mix of frequencies.", go(k) { k.show(k.s.head); }, g(k) { k.draw(k.s.mix, 1400); k.show(k.s.mixLab, {delay: 500}); }},
      {say: "[n]SASSI takes it apart, like a chord split into its notes, and solves one note at a time.",
        n(k) { k.show(k.s.arrow); k.show(k.s.notes, {stagger: 220}); k.show(k.s.notesLab, {delay: 900}); }},
      {say: "[t]For every point it finds a transfer function: how much that point amplifies the ground, at each frequency.",
        t(k) { k.show(k.s.p.g); k.draw(k.s.tf, 1600); k.show(k.s.tfLab); k.show(k.s.eqz, {delay: 600}); }},
      {say: "[b]Then it puts the notes back together with the real earthquake, to get floor motions and floor spectra.",
        b(k) { k.show(k.s.back); }},
    ],
  });

  // ---------------------------------------------------------------- 6. the example
  scenes.push({
    id: "model", title: "Our example",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Our example", {x: 120, y: 80, size: "h2"});
      // to scale, 6.25 px per ft: 64 ft mat 5 ft thick, 16 ft storeys, 16 ft sand on 39 ft gravel (rock cut with a fade)
      const MS = 6.25, gy = 590, top = gy - 5 * MS;
      s.soil = K.soil(g, {x: 120, y: gy, w: 760, layers: [
        {h: 16 * MS, kind: "sand", name: "Sand, 16 ft"}, {h: 39 * MS, kind: "gravel", name: "Gravel, 39 ft"}],
        hs: {h: 110, kind: "rock", name: "Rock"}, labels: "right", labelSize: 32});
      s.soil.layers.forEach((L) => { L.g.classList.add("sv-in", "sv-fade"); });
      const zs = [16, 32, 48, 64].map((z) => z * MS);
      s.stick = K.stick(g, {x: 500, y: top, z: zs, matW: 64 * MS, matH: 5 * MS, r: 24, slabW: 170, hidden: true});
      s.dim = K.g(g, {hidden: true});
      K.dim(s.dim, 500 - 32 * MS, top - 45, 500 + 32 * MS, top - 45, "");
      K.text(s.dim, 360, top - 55, "64 ft", {cls: "t-label", anchor: "middle", size: 28});
      s.masses = K.g(g, {hidden: true});
      zs.forEach((z) => K.text(s.masses, 610, top - z + 9, "2,200 kips", {cls: "t-label", size: 30, color: "var(--ink)"}));
      s.code = K.code(s.root, ["INP,ex01_surface_stick.pre"], {x: 1180, y: 520, w: 640, size: 32, title: "One command"});
      s.codeLab = K.para(s.root, "runs the whole SASSI chain in about <b>4 seconds</b>", {x: 1180, y: 680, w: 640, size: 32});
    },
    beats: [
      {say: "[m]Our example is a four-storey stick on a 64-foot square mat, with 2,200 kips per floor.",
        go(k) { k.show(k.s.head); }, m(k) { k.show([k.s.stick.g, k.s.dim]); k.show(k.s.masses, {delay: 500}); }},
      {say: "It sits on [s]16 ft of sand, [g]39 ft of gravel, [r]and weathered rock.",
        s(k) { k.show(k.s.soil.g); k.show(k.s.soil.layers[0].g); }, g(k) { k.show(k.s.soil.layers[1].g); }, r(k) { k.show(k.s.soil.layers[2].g); }},
      {say: "[c]One command runs the whole SASSI chain, in about four seconds.",
        c(k) { k.show(k.s.code.el); k.type(k.s.code.lines, {cps: 34}); k.show(k.s.codeLab, {delay: 800}); }},
    ],
  });

  // ---------------------------------------------------------------- 7. the roof's transfer function
  scenes.push({
    id: "roof", title: "The roof's transfer function",
    build(s) {
      s.p = K.plot(s.svg, {x: 240, y: 170, w: 1000, h: 660, xr: [0, 20], yr: [0, 14], xticks: [0, 5, 10, 15, 20], yticks: [0, 2, 4, 6, 8, 10, 12, 14],
        xlabel: "frequency (Hz)", ylabel: "roof motion ÷ ground motion", ylabelOffset: 80});
      const roof = D["ex01/00085TR_X.TFI"], pk = peak(roof.f, roof.amp);
      s.roof = s.p.line(roof.f, roof.amp, {color: C.ssi, width: 7, draw: true});
      s.one = K.g(s.p.g, {hidden: true});
      K.circle(s.one, s.p.X(0.098), s.p.Y(1), 12, {fill: C.ink});
      K.line(s.one, s.p.X(0.098) + 10, s.p.Y(1) + 10, s.p.X(0.9), s.p.Y(0.45), {stroke: C.ink, "stroke-width": 3});
      K.text(s.one, s.p.X(1.0), s.p.Y(0.3), "1 at low frequency", {cls: "t-label", color: "var(--ink)"});
      s.pk = K.g(s.p.g, {hidden: true});
      K.circle(s.pk, s.p.X(pk.f), s.p.Y(pk.v), 14, {fill: C.ssi});
      K.text(s.pk, s.p.X(pk.f) + 30, s.p.Y(pk.v) + 12, `${pk.v.toFixed(1)} × the ground at ${pk.f.toFixed(2)} Hz`, {cls: "t-label", color: "var(--ssi)", size: 34});
      s.s1 = K.stat(s.root, {x: 1340, y: 260, w: 460, value: "1", label: "everything moves together", color: "var(--ink)", vsize: 120});
      s.s2 = K.stat(s.root, {x: 1340, y: 500, w: 460, value: "13×", label: "the roof at 3.5 Hz:<br>the system's sway", color: "var(--ssi)", vsize: 120});
      s.ansys = K.card(s.root, {x: 1320, y: 760, w: 500, kind: "ansys", title: "In ANSYS terms", size: 30, body: "A harmonic analysis with a <b>unit base motion</b>."});
    },
    beats: [
      {say: "[p]Here is the roof's transfer function.", p(k) { k.show(k.s.p.g); k.draw(k.s.roof, 2200); }},
      {say: "[o]At very low frequency it's exactly 1: everything moves together with the ground.", o(k) { k.show([k.s.one, k.s.s1.el]); }},
      {say: "[k]Then it peaks at 3.5 hertz, where the roof moves 13 times as much as the ground.", k(k) { k.show([k.s.pk, k.s.s2.el]); }},
      {say: "[a]In ANSYS terms, it's what a harmonic analysis gives you for a unit shake of the base, frequency by frequency.",
        a(k) { k.show(k.s.ansys); }},
    ],
  });

  // ---------------------------------------------------------------- 8. what the soil changed
  scenes.push({
    id: "compare", title: "What the soil changed",
    build(s) {
      s.head = K.heading(s.root, "Fixed base against soil", {x: 120, y: 70, size: "h2"});
      const tfS = D["ex01/00085TR_X.TFI"], tfF = D["ex01_fixed/00085TR_X.TFI"];
      s.p1 = K.plot(s.svg, {x: 200, y: 200, w: 700, h: 470, xr: [0, 10], yr: [0, 14], xticks: [0, 2, 4, 6, 8, 10], yticks: [0, 7, 14], xlabel: "frequency (Hz)", ylabel: "roof ÷ ground", ylabelOffset: 70});
      s.lF = s.p1.line(tfF.f, tfF.amp, {color: C.ref, width: 6, draw: true});
      s.lS = s.p1.line(tfS.f, tfS.amp, {color: C.ssi, width: 7, draw: true});
      s.p1.legend([{label: "fixed base", color: C.ref}, {label: "on soil", color: C.ssi}], {x: 640, y: 240, size: 28});
      const pF = peak(tfF.f, tfF.amp), pS = peak(tfS.f, tfS.amp);
      s.shift = K.arrow(s.p1.g, s.p1.X(pF.f) - 6, s.p1.Y(14.6), s.p1.X(pS.f) + 6, s.p1.Y(14.6), {color: C.ink, width: 6, head: 22, hidden: true});
      s.pct = K.stat(s.root, {x: 250, y: 790, w: 600, value: "−30 %", label: "sway frequency: 5 → 3.5 Hz", color: "var(--ssi)", vsize: 100});
      const rsS = D["ex01/00085TR_X02.RS"], rsF = D["ex01_fixed/00085TR_X02.RS"];
      s.p2 = K.plot(s.svg, {x: 1120, y: 200, w: 700, h: 470, xr: [0.5, 30], xlog: true, yr: [0, 9], xticks: [0.5, 1, 2, 5, 10, 20], yticks: [0, 3, 6, 9], xlabel: "frequency (Hz)", ylabel: "roof spectrum (g)", ylabelOffset: 70});
      s.rF = s.p2.line(rsF.f, rsF.sa, {color: C.ref, width: 6, draw: true});
      s.rS = s.p2.line(rsS.f, rsS.sa, {color: C.ssi, width: 7, draw: true});
      s.rsLab = K.g(s.p2.g, {hidden: true});
      const qF = peak(rsF.f, rsF.sa, 1, 20), qS = peak(rsS.f, rsS.sa, 1, 20);
      K.text(s.rsLab, s.p2.X(qF.f) + 18, s.p2.Y(qF.v) - 18, `${g2(qF.v)} g`, {cls: "t-label", color: "var(--ref)"});
      K.text(s.rsLab, s.p2.X(qS.f) - 18, s.p2.Y(qS.v) - 18, `${g2(qS.v)} g`, {cls: "t-label", anchor: "end", color: "var(--ssi)"});
      s.warn = K.card(s.root, {x: 1120, y: 790, w: 700, kind: "warn", title: "Not automatically good news", size: 32, body: "The peak moved and <b>stayed as high</b>. Compute it."});
    },
    beats: [
      {say: "[f]Now run the same building on rigid ground: the fixed-base case.", go(k) { k.show(k.s.head); }, f(k) { k.show(k.s.p1.g); k.draw(k.s.lF, 1500); k.draw(k.s.lS, 1500, {delay: 600}); }},
      {say: "[s]The soil lowered the sway frequency by 30 percent, from 5 to 3.5 hertz.", s(k) { k.show(k.s.shift); k.draw(k.s.shift.line, 700); k.show(k.s.pct.el, {delay: 300}); }},
      {say: "[r]And the floor spectrum peak moved right along with it.", r(k) { k.show(k.s.p2.g); k.draw(k.s.rF, 1300); k.draw(k.s.rS, 1300, {delay: 500}); }},
      {say: "[h]Notice the peak did not drop: 8.03 g on soil, the same as on a fixed base.", h(k) { k.show(k.s.rsLab); }},
      {say: "[w]So SSI is not automatically good news. You have to compute it.", w(k) { k.show(k.s.warn); }},
    ],
  });

  // ---------------------------------------------------------------- 9. the basemat shakes too
  scenes.push({
    id: "basemat", title: "The foundation shakes too",
    build(s) {
      s.head = K.heading(s.root, "The foundation shakes too", {x: 120, y: 80, size: "h2"});
      const rsS = D["ex01/00041TR_X02.RS"], rsF = D["ex01_fixed/00041TR_X02.RS"];
      s.p = K.plot(s.svg, {x: 240, y: 210, w: 900, h: 560, xr: [0.5, 30], xlog: true, yr: [0, 1.6], xticks: [0.5, 1, 2, 5, 10, 20], yticks: [0, 0.4, 0.8, 1.2, 1.6], xlabel: "frequency (Hz)", ylabel: "basemat spectrum (g)", ylabelOffset: 80});
      s.lF = s.p.line(rsF.f, rsF.sa, {color: C.wave, width: 6, draw: true});
      s.lS = s.p.line(rsS.f, rsS.sa, {color: C.ssi, width: 7, draw: true});
      s.p.legend([{label: "the ground alone", color: C.wave}, {label: "basemat, on soil", color: C.ssi}], {x: 820, y: 250, size: 28});
      s.pk = K.g(s.p.g, {hidden: true});
      const q = peak(rsS.f, rsS.sa, 1, 10);
      K.circle(s.pk, s.p.X(q.f), s.p.Y(q.v), 12, {fill: C.ssi});
      K.text(s.pk, s.p.X(q.f) + 24, s.p.Y(q.v) - 16, `${g2(q.v)} g`, {cls: "t-label", color: "var(--ssi)", size: 34});
      s.why = K.card(s.root, {x: 1240, y: 300, w: 580, title: "Why?", size: 34, body: "The swaying building <b>rocks its own foundation</b>."});
      s.eq = K.card(s.root, {x: 1240, y: 560, w: 580, kind: "warn", title: "Equipment on the basemat", size: 32, body: "feels more than the ground: don't use the free-field spectrum there."});
    },
    beats: [
      {say: "[p]One more surprise: look at the basemat itself.", go(k) { k.show(k.s.head); }, p(k) { k.show(k.s.p.g); k.draw(k.s.lF, 1300); k.draw(k.s.lS, 1300, {delay: 500}); }},
      {say: "[k]Its spectrum peaks at 1.39 g near 3.2 hertz, where the ground alone gives about 1.0.", k(k) { k.show(k.s.pk); }},
      {say: "[w]Why? The swaying building rocks its own foundation.", w(k) { k.show(k.s.why); k.show(k.s.eq, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 10. recap
  scenes.push({
    id: "recap", title: "Recap",
    build(s) {
      s.head = K.heading(s.root, "Recap", {x: 120, y: 90, size: "h1"});
      s.list = K.bullets(s.root, [
        "On soil, the ground acts like a <b style='color:var(--spring)'>suspension</b>",
        "SASSI joins soil and building at the <b style='color:var(--wave)'>interaction nodes</b>, one frequency at a time",
        "Here the sway frequency fell from <b style='color:var(--ssi)'>5 to 3.5 Hz</b>, and the spectra moved with it",
      ], {x: 120, y: 240, w: 1000, num: true, size: 40});
      s.q = K.card(s.root, {x: 1200, y: 240, w: 620, kind: "check", title: "Quick check", size: 34, body: "Does SSI always make floor spectra smaller?"});
      s.a = K.card(s.root, {x: 1200, y: 520, w: 620, title: "Answer", size: 34, body: "<b>No.</b> Here the peak moved and stayed as high. That's why you compute it."});
      s.next = K.pill(s.root, "Next · Lesson 2: the soil layers and the free field →", {x: 120, y: 860, color: "var(--wave)"});
    },
    beats: [
      {say: "[a]On soil, the ground acts like a suspension, with springs and shock absorbers.", go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[a]SASSI joins the soil and the building at the interaction nodes, one frequency at a time.", a(k) { k.show(k.s.list.items[1]); }},
      {say: "[a]Here, the soil moved the sway frequency from 5 to 3.5 hertz, and the floor spectra moved with it.", a(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Quick check: does SSI always make floor spectra smaller?", q(k) { k.show(k.s.q); }, gap: 1400},
      {say: "[a]No. Here the peak just moved, and stayed as high. That's why you compute it.", a(k) { k.show(k.s.a); }},
      {say: "[n]Next, in lesson two: the soil layers, and how the earthquake travels up through them.", n(k) { k.show(k.s.next); }},
    ],
  });

  SV.video({
    id: "01", n: 1, part: "Fundamentals", lesson: "01-why-ssi",
    title: "From fixed base to SSI: what changes and why",
    subtitle: "Soil acts like a suspension: what that does to a building's sway and its floor spectra, and how SASSI computes it.",
    next: {href: "02.html", title: "The site: layers, half-space and free-field motion"},
    scenes,
  });
})();
