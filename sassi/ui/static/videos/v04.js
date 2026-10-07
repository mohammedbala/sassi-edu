/* Explainer video, lesson 4 (newcomer edition): Your first SSI analysis -- a stick on a surface mat.
 * The big idea: an SSI model is your ANSYS model plus a frequency list and the nodes that touch the soil; the
 * SASSI chain turns it into floor spectra, displacements and forces; the frequency drops from 5 to 3.5 Hz mostly
 * because the mat rocks on the soil, and the peak stays high because rocking radiates almost no energy here.
 * Numbers: sassi/ui/lessons/04_surface_ssi.md; curves: d04.js (the lesson's run, python -m sassi.ui.video_data).
 * Units: ft, kip, s.  Pictures to scale (64 ft mat, 16 ft storeys, 16 ft sand on 39 ft gravel); every moving
 * picture is the computed steady state for 0.04 in of ground motion, all displacements x 400, slowed down 7 times. */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["04"];
  const PH = K.PH, Cx = PH.C;
  const TAU = 2 * Math.PI;
  const F = (name) => D["ex01/" + name];
  const argmax = (a) => a.reduce((b, v, i) => (Math.abs(v) > Math.abs(a[b]) ? i : b), 0);

  // ---------------------------------------------------------------- the run's values at 3.49 Hz
  /** Complex transfer function of a .TFU file at the analysis frequency nearest f. */
  function tfAt(name, f) {
    const o = F(name);
    let i = 0;
    o.f.forEach((v, j) => { if (Math.abs(v - f) < Math.abs(o.f[i] - f)) i = j; });
    return [o.amp[i] * Math.cos(o.ph[i]), o.amp[i] * Math.sin(o.ph[i])];
  }
  /** Steady state at 3.49 Hz per unit control motion (MOTION output of the run): mat centre X (node 41), floors X
   *  (82..85), and the rocking theta = uz(37) / 32 ft (node 37 is the mat edge at x = -32 ft). */
  const HM = {
    mat: tfAt("00041TR_X.TFU", 3.49),
    fl: ["00082TR_X.TFU", "00083TR_X.TFU", "00084TR_X.TFU", "00085TR_X.TFU"].map((n) => tfAt(n, 3.49)),
    th: Cx.scale(tfAt("00037TR_Z.TFU", 3.49), 1 / 32),
  };
  /** The same stick clamped to rigid ground (lesson 4: E = 576,000 ksf, I = 23,000 ft4, shear area 108 ft2,
   *  16 ft storeys, 2,200 kips per floor, 5 % hysteretic damping; kip, ft, s): floor motion relative to the base
   *  per unit base motion at 4.97 Hz, its first mode (4.99 Hz elastic; roof 13.4, as lesson 1's fixed-base run). */
  const HFB = (function () {
    const E = 576000, I = 23000, GAs = (576000 / 2.4) * 108, z = [16, 32, 48, 64], m = 2200 / 32.2, w = TAU * 4.97;
    const flex = z.map((zi) => z.map((zj) => { const a = Math.min(zi, zj), b = Math.max(zi, zj); return a * a * (3 * b - a) / (6 * E * I) + a / GAs; }));
    const cols = z.map((_, j) => PH.csolve(flex.map((r) => r.map((v) => [v, 0])), z.map((_, i) => [i === j ? 1 : 0, 0])).x);
    const cf = PH.cfac(0.05);
    const A = z.map((_, i) => z.map((_, j) => Cx.sub(Cx.scale(cf, cols[j][i][0]), [i === j ? w * w * m : 0, 0])));
    return PH.csolve(A, z.map(() => [w * w * m, 0])).x;
  })();
  /** The free field of the site (lesson 2's column) at 3.49 Hz: displacement at depth z (ft) per unit surface motion. */
  const WV = PH.columnWaves(PH.COLUMN, 3.49);
  const colU = (z) => PH.columnU(WV, PH.COLUMN, z);
  // foundation impedance of the run (FOUNSTIF, FOUNDAMP): sliding X and rocking about Y
  const imp = {f: F("FOUNDAMP").cols[0], dx: F("FOUNDAMP").cols[1], dyy: F("FOUNDAMP").cols[29],
    kx: F("FOUNSTIF").cols[1], kyy: F("FOUNSTIF").cols[29]};
  const iF = (f) => { let i = 0; imp.f.forEach((v, j) => { if (Math.abs(v - f) < Math.abs(imp.f[i] - f)) i = j; }); return i; };
  // soil flexibility for a push at the roof (h = 64 ft), in per 1,000 kips: 1/Kx = 0.0119 and h^2/Kyy = 0.0416 (lesson 4)
  const i349 = iF(3.49);
  const FLEX = {slide: 12e3 / imp.kx[i349], rock: 12e3 * 64 * 64 / imp.kyy[i349]};

  // ---------------------------------------------------------------- the stick on its mat, to scale
  const XD = 400;                                   // displacement exaggeration of every moving picture
  const UG = 0.04 / 12;                             // ground motion of every moving picture: 0.04 in, in ft
  const TV = 2.0;                                   // visual period at 3.49 Hz (s): slowed down 7 times
  const CAP = "computed at 3.49 Hz: ground motion 0.04 in, displacements × 400, slowed down 7 ×";
  const CAPFB = "computed at 4.97 Hz: ground motion 0.04 in, displacements × 400, slowed down 7 ×";
  /** The stick of example 1 on its mat over the site, with sliding and rocking springs (r.springs) and dashpots
   *  (r.dashes), both hidden.  o: x (centre), y (top of the mat), s (px per foot), soilW, soilHidden, hsH (px of
   *  half-space shown).  Displacements: o.A = px for 0.04 in of ground motion (XD x s x UG). */
  function rig(g, o) {
    const r = {o};
    o.A = XD * o.s * UG;
    const half = 32 * o.s, matH = 5 * o.s, ys = o.y + matH, yb = ys + (o.springL || 118);
    r.half = half; r.ys = ys; r.yb = yb;
    r.soil = K.soil(g, {x: o.x - o.soilW / 2, y: ys, w: o.soilW, labels: "inside", labelSize: 24, hidden: o.soilHidden,
      layers: [{h: 16 * o.s, kind: "sand", name: "sand · 16 ft"}, {h: 39 * o.s, kind: "gravel", name: "gravel · 39 ft"}],
      hs: {h: o.hsH || 16 * o.s, kind: "rock", name: "rock"}});
    r.springs = K.g(g, {hidden: true});
    r.dashes = K.g(g, {hidden: true});
    r.vx = [-(half - 30), half - 30]; r.dx = [-(half - 82), half - 82];
    r.vsp = r.vx.map((d) => K.spring(r.springs, o.x + d, ys, o.x + d, yb, {coils: 5, amp: 12}));
    r.vdp = r.dx.map((d) => K.dashpot(r.dashes, o.x + d, ys, o.x + d, yb, {w: 24, cyl: 56}));
    r.gndS = K.g(r.springs, {}); r.gndD = K.g(r.dashes, {});
    r.vx.forEach((d) => K.ground(r.gndS, o.x + d - 26, yb, 52, {color: "var(--muted)"}));
    r.dx.forEach((d) => K.ground(r.gndD, o.x + d - 26, yb, 52, {color: "var(--muted)"}));
    r.wx = o.x - half - 130;
    r.hy = [o.y + 2, ys + 8];                       // the sliding spring and dashpot: at the mat, horizontal
    r.hsp = K.spring(r.springs, r.wx, r.hy[0], o.x - half, r.hy[0], {coils: 5, amp: 9});
    r.hdp = K.dashpot(r.dashes, o.x - half, r.hy[1], r.wx, r.hy[1], {w: 18, cyl: 50});
    r.wall = K.line(r.springs, r.wx, o.y - 30, r.wx, ys + 30, {stroke: "var(--muted)", "stroke-width": 5});
    r.stick = K.stick(g, {x: o.x, y: o.y, z: [1, 2, 3, 4].map((i) => i * 16 * o.s), matW: 2 * half, matH, r: o.r || 20, slabW: o.slabW || 22 * o.s});
    return r;
  }
  /** Pose the rig at time t: the computed steady state at 3.49 Hz split into its parts, weighted by w = {slide,
   *  rock, bend, ground} (1 = as computed), plus w.fixed x the clamped stick at 4.97 Hz (hook). */
  function pose(r, t, w) {
    const o = r.o, A = o.A, ph = TAU * t / TV, c = Math.cos(ph), sn = Math.sin(ph);
    const re = (z) => z[0] * c - z[1] * sn;
    const pf = TAU * t / (TV * 3.49 / 4.97), cf = Math.cos(pf), sf = Math.sin(pf), wf = w.fixed || 0;
    const ug = A * ((w.ground || 0) * c + wf * cf);                       // the ground surface (free field)
    const sway = A * (w.slide * re(HM.mat) + wf * cf), rock = A * w.rock * re(HM.th) / o.s;
    const bend = HM.fl.map((H, i) => A * w.bend * re(Cx.sub(Cx.sub(H, HM.mat), Cx.scale(HM.th, 16 * (i + 1))))
      + A * wf * (HFB[i][0] * cf - HFB[i][1] * sf));
    r.stick.set({sway, rock, bend});
    const sr = Math.sin(rock);
    r.vsp.forEach((sp, i) => { const d = r.vx[i]; sp.update(o.x + d + sway, r.ys + d * sr, o.x + d + ug, r.yb); });
    r.vdp.forEach((dp, i) => { const d = r.dx[i]; dp.update(o.x + d + sway, r.ys + d * sr, o.x + d + ug, r.yb); });
    r.hsp.update(r.wx + ug, r.hy[0], o.x - r.half + sway, r.hy[0] - r.half * sr);
    r.hdp.update(o.x - r.half + sway, r.hy[1] - r.half * sr, r.wx + ug, r.hy[1]);
    for (const e of [r.wall, r.gndS, r.gndD]) e.setAttribute("transform", `translate(${ug} 0)`);
    // the soil column: the one-dimensional wave solution of the site at 3.49 Hz (a rigid ground for the clamped case)
    r.soil.shear((d) => { const u = colU(d / o.s); return A * ((w.ground || 0) * re(u) + wf * cf); });
  }
  /** A module box (HTML): name in mono, a plain-words caption under it (hidden until named). */
  function modBox(root, x, y, w, h, name, cap, color) {
    const b = K.h("div", {x, y, w, h, in: "pop", style: {border: "3px solid #33507a", borderRadius: "16px", background: "rgba(20,32,51,.94)",
      display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "center", gap: "8px"}}, root);
    const n = document.createElement("div");
    n.style.cssText = `font: 700 42px var(--mono); color: ${color || "var(--ink)"}`;
    n.textContent = name;
    b.appendChild(n);
    const c = document.createElement("div");
    c.className = "sv-in sv-fade";
    c.style.cssText = "font: 500 29px var(--ui); color: var(--ink2)";
    c.textContent = cap;
    b.appendChild(c);
    b.cap = c;
    return b;
  }
  /** A small arrow (HTML) between two boxes of a row. */
  const rowArrow = (root, x, y, w) => K.h("div", {x, y, w, h: 36, in: "fade",
    html: `<svg viewBox="0 0 44 36" width="${w}" height="36"><path d="M2 18H36M26 8l12 10-12 10" stroke="#5f7da6" stroke-width="4" fill="none" stroke-linecap="round" stroke-linejoin="round"/></svg>`}, root);

  // ================================================================== scenes
  const scenes = [];

  // ---------------------------------------------------------------- 0. hook
  scenes.push({
    id: "hook", title: "A familiar building",
    build(s) {
      const g = K.g(s.svg);
      s.pic = K.g(g, {hidden: true});
      s.r = rig(s.pic, {x: 540, y: 520, s: 5.625, soilW: 760, soilHidden: true, hsH: 80});
      s.hatch = K.ground(s.pic, 540 - 260, s.r.ys, 520, {color: "var(--ref)"});
      K.dim(s.pic, 790, 520, 790, 520 - 64 * 5.625, "64 ft");
      s.capF = K.text(g, 540, 975, CAPFB, {cls: "t-small", anchor: "middle", hidden: true});
      s.capS = K.text(g, 540, 975, CAP, {cls: "t-small", anchor: "middle", hidden: true});
      s.fb = K.stat(s.root, {x: 1080, y: 210, w: 680, value: "4.97", unit: "Hz", label: "clamped to rigid ground", color: "var(--ref)", vsize: 120});
      s.ssi = K.stat(s.root, {x: 1080, y: 470, w: 680, value: "3.49", unit: "Hz", label: "on real soil", color: "var(--ssi)", vsize: 120});
      s.q = K.heading(s.root, '<span style="color:var(--ssi)">What slows it down?</span>', {x: 1080, y: 770, w: 680, align: "center", size: "h1"});
      s.w = {slide: 0, rock: 0, bend: 0, ground: 0, fixed: 0};
    },
    tick(s, t) {
      pose(s.r, t, s.w);
      s.hatch.setAttribute("transform", `translate(${s.r.o.A * s.w.fixed * Math.cos(TAU * t / (TV * 3.49 / 4.97))} 0)`);
    },
    beats: [
      {say: "You've built this model before: a four-storey stick on a concrete mat, with heavy floors.",
        go(k) { k.show(k.s.pic); }},
      {say: "Clamped to rigid ground, it sways [a]at 5 Hz. On real soil, it sways [b]at 3.5.",
        a(k) { k.show(k.s.fb.el); k.show(k.s.capF); k.tween(k.s.w, {fixed: 1}, 900); },
        b(k) {
          const s = k.s;
          k.show(s.ssi.el); k.hide([s.hatch, s.capF]); k.show([s.r.soil.g, s.capS]);
          k.tween(s.w, {fixed: 0, slide: 1, rock: 1, bend: 1, ground: 1}, 1000);
        }},
      {say: "Today you'll build that SSI model yourself, [q]and find out what slows it down.", hold: 600,
        q(k) { k.show(k.s.q); }},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 4",
    build(s) {
      s.t = K.titleCard(s.root, {n: 4, part: "Fundamentals", title: "Your first SSI analysis: a stick on a surface mat",
        sub: "Build it, run it, read it, and see why it sways at 3.5 Hz."});
      s.t.title.style.maxWidth = "1180px";
    },
    beats: [
      {say: "Lesson four. Your first SSI analysis: a stick on a surface mat.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 100}); k.show(t.title, {delay: 200}); k.show(t.sub, {delay: 350}); k.show(t.rule, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. the model (isometric view, to scale)
  scenes.push({
    id: "model", title: "Build it like ANSYS",
    build(s) {
      K.defs();
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Build it like an ANSYS model", {x: 120, y: 70, size: "h2"});
      const CX = 600, CY = 580, S = 5.9375;           // px per foot
      const P = (x, y, z) => [CX + (x - y) * S * 0.82, CY + (x + y) * S * 0.36 - (z || 0) * S];
      const quad = (gg, pts, a) => K.path(gg, K.d(pts) + "Z", a);
      // the site of lesson 2 to scale: 16 ft of sand on 39 ft of gravel (cut, fading out below 30 ft)
      s.soil = K.g(g, {hidden: true});
      const E = 42;
      quad(s.soil, [P(-E, -E, 0), P(E, -E, 0), P(E, E, 0), P(-E, E, 0)], {fill: C.sand, "fill-opacity": 0.32, stroke: "rgba(255,240,215,.35)", "stroke-width": 2});
      const strips = [[0, -16, C.sand, "sv-pat-sand", 1]].concat([0, 1, 2, 3, 4, 5, 6, 7, 8].map((i) => [-16 - 1.6 * i, -17.6 - 1.6 * i, C.gravel, "sv-pat-gravel", 1 - i / 9]));
      strips.forEach(([z1, z2, col, pat, op]) => {
        const fr = [P(E, -E, z1), P(E, E, z1), P(E, E, z2), P(E, -E, z2)], fl = [P(-E, E, z1), P(E, E, z1), P(E, E, z2), P(-E, E, z2)];
        quad(s.soil, fr, {fill: col, "fill-opacity": 0.8 * op}); quad(s.soil, fr, {fill: `url(#${pat})`, opacity: op});
        quad(s.soil, fl, {fill: col, "fill-opacity": 0.62 * op}); quad(s.soil, fl, {fill: `url(#${pat})`, opacity: op});
      });
      K.text(s.soil, P(0, E, -8)[0], P(0, E, -8)[1] + 10, "sand · 16 ft", {cls: "t-label", anchor: "middle", size: 26});
      K.text(s.soil, P(0, E, -22)[0], P(0, E, -22)[1] + 10, "gravel · 39 ft", {cls: "t-label", anchor: "middle", size: 26});
      // the mat: 64 ft x 64 ft, 5 ft thick
      const M = 32, T = 5;
      s.mat = K.g(g, {hidden: true});
      quad(s.mat, [P(M, -M, 0), P(M, M, 0), P(M, M, -T), P(M, -M, -T)], {fill: "#56677d", stroke: C.concrete, "stroke-width": 2});
      quad(s.mat, [P(-M, M, 0), P(M, M, 0), P(M, M, -T), P(-M, M, -T)], {fill: "#435267", stroke: C.concrete, "stroke-width": 2});
      quad(s.mat, [P(-M, -M, 0), P(M, -M, 0), P(M, M, 0), P(-M, M, 0)], {fill: "var(--concrete-fill)", stroke: C.concrete, "stroke-width": 3});
      s.grid = K.g(g, {hidden: true});
      for (let i = 1; i < 8; i++) {
        const v = -M + 8 * i;
        K.path(s.grid, K.d([P(v, -M), P(v, M)]), {stroke: "rgba(223,230,238,.42)", "stroke-width": 2});
        K.path(s.grid, K.d([P(-M, v), P(M, v)]), {stroke: "rgba(223,230,238,.42)", "stroke-width": 2});
      }
      const rowPts = (j) => K.linspace(-M, M, 9).map((x) => P(x, -M + 8 * j));
      s.nodes = K.nodes(g, [].concat(...K.linspace(0, 8, 9).map(rowPts)), {r: 6, color: "var(--ink2)", hidden: true});
      s.irows = K.linspace(0, 8, 9).map((j) => K.nodes(g, rowPts(j), {r: 9, color: C.wave, hidden: true}));
      // the stick (4 storeys of 16 ft) and the floor masses
      s.stick = K.g(g, {hidden: true});
      K.path(s.stick, K.d([P(0, 0, 0), P(0, 0, 64)]), {cls: "struct"});
      s.masses = [16, 32, 48, 64].map((z) => { const q = P(0, 0, z); return K.circle(g, q[0], q[1], 20, {cls: "mass", fill: "url(#sv-grad-mass)", hidden: true, in: "pop"}); });
      s.labs = K.g(g, {});
      s.labMat = K.text(s.labs, P(M, -M)[0] + 28, P(M, -M)[1] - 8, "shell mat", {cls: "t-label", hidden: true});
      s.labMat2 = K.text(s.labs, P(M, -M)[0] + 28, P(M, -M)[1] + 26, "64 ft × 64 ft", {cls: "t-small", hidden: true});
      s.labMat3 = K.text(s.labs, P(M, -M)[0] + 28, P(M, -M)[1] + 56, "nodes every 8 ft", {cls: "t-small", hidden: true});
      s.labStick = K.text(s.labs, CX - 40, P(0, 0, 40)[1], "beam stick", {cls: "t-label", anchor: "end", hidden: true});
      s.labStick2 = K.text(s.labs, CX - 40, P(0, 0, 40)[1] + 32, "4 storeys × 16 ft", {cls: "t-small", anchor: "end", hidden: true});
      s.labMass = K.text(s.labs, CX + 40, P(0, 0, 64)[1] + 10, "2,200 kips per floor", {cls: "t-label", hidden: true});
      // right column
      s.ansys = K.card(s.root, {x: 1150, y: 190, w: 650, kind: "ansys", title: "In ANSYS terms", size: 32,
        body: "The same model: shells, beams and point masses."});
      s.noSoil = K.card(s.root, {x: 1150, y: 190, w: 650, kind: "warn", title: "No soil mesh", size: 32,
        body: "Not a single soil element."});
      s.code = K.code(s.root, ["GROUP,1,SHELL", "GROUP,2,BEAMS", "INT,1,81,1,1"], {x: 1150, y: 410, w: 650, size: 36, title: "SASSI commands"});
      s.inPill = K.pill(s.root, "81 interaction nodes", {x: 1150, y: 700, color: "var(--wave)", style: {fontSize: "32px"}});
      s.hands = K.label(s.root, "where building and soil meet", {x: 1154, y: 790, size: 32, color: "var(--ink2)", in: "fade"});
    },
    beats: [
      {say: "[m]Start with the mat: a {64 ft|64-foot} square concrete slab, [s]made of shell elements.",
        go(k) { k.show(k.s.head); },
        m(k) { k.show(k.s.mat); },
        s(k) { k.show(k.s.grid); k.show(k.s.nodes.g, {delay: 300}); k.show([k.s.labMat, k.s.labMat2, k.s.labMat3], {delay: 200, stagger: 120}); }},
      {say: "[b]Then the building: a stick of beam elements, [m]with 2,200 kips at each floor.",
        b(k) { k.show([k.s.stick, k.s.labStick, k.s.labStick2]); },
        m(k) { k.show(k.s.masses, {stagger: 140}); k.show(k.s.labMass, {delay: 400}); }},
      {say: "[a]So far, it's exactly the model you'd build in ANSYS.",
        a(k) { k.show(k.s.ansys); k.show(k.s.code.el, {delay: 200}); k.type(k.s.code.lines.slice(0, 2), {cps: 40}); }},
      {say: "[i]Then one new line: INT marks the nodes that touch the soil, [n]all 81 nodes of the mat.",
        i(k) { k.type(k.s.code.lines.slice(2), {cps: 30}); k.pulse(k.s.code.lines[2], {amp: 0.05, ms: 900}); },
        n(k) { k.show(k.s.irows.map((r) => r.g), {stagger: 90}); k.show(k.s.inPill, {delay: 400}); }},
      {say: "Engineers call them interaction nodes: [h]the points where the building and the soil hold hands.",
        h(k) { k.show(k.s.hands); k.pulse(k.s.irows.map((r) => r.g), {amp: 0.05, ms: 800}); }},
      {say: "[x]And notice what's missing: there is no soil mesh at all.",
        x(k) { k.hide([k.s.ansys, k.s.labMat, k.s.labMat2, k.s.labMat3, k.s.labStick, k.s.labStick2, k.s.labMass]); k.show(k.s.noSoil, {delay: 250}); }},
      {say: "[s]SASSI brings in the site you built in lesson 2, and joins it only at those nodes.", hold: 600,
        s(k) { k.show(k.s.soil); k.pulse(k.s.irows.map((r) => r.g), {amp: 0.06, ms: 900}); }},
    ],
  });

  // ---------------------------------------------------------------- 3. the frequency list
  scenes.push({
    id: "frequencies", title: "Choose the frequencies",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Choose the frequencies", {x: 120, y: 70, size: "h2"});
      // a chord and its notes (right column, an illustration)
      s.ch = K.g(g, {hidden: true});
      K.text(s.ch, 1565, 214, "a chord", {cls: "t-label", anchor: "middle"});
      s.chord = K.path(s.ch, "", {stroke: C.wave, "stroke-width": 4, fill: "none", "stroke-linejoin": "round"});
      K.arrow(s.ch, 1565, 322, 1565, 372, {color: C.muted, width: 4, head: 16});
      s.notes = [0, 1, 2].map(() => K.path(s.ch, "", {stroke: C.ssi, "stroke-width": 3.5, fill: "none"}));
      K.text(s.ch, 1565, 548, "its notes (illustration)", {cls: "t-label", anchor: "middle"});
      s.n22 = K.stat(s.root, {x: 1330, y: 600, w: 470, value: "22", label: "frequencies, 0.1 to 20 Hz", color: "var(--ink)", vsize: 110});
      // the roof transfer function of the run and its 22 analysis frequencies
      const tfi = F("00085TR_X.TFI"), tfu = F("00085TR_X.TFU");
      const PX = 200, PY = 250, PW = 1000, PH = 520;
      s.p = K.plot(s.svg, {x: PX, y: PY, w: PW, h: PH, xr: [0, 25], yr: [0, 14], xticks: [0, 5, 10, 15, 20, 25], yticks: [0, 4, 8, 12],
        xlabel: "frequency (Hz)", ylabel: "roof ÷ ground motion", ylabelOffset: 70});
      // the dense band: where the analysis frequencies are about 0.5 Hz apart (from the run's list)
      const fd = tfu.f.filter((f, i) => i > 0 && f - tfu.f[i - 1] < 0.6);
      s.bandA = K.g(s.p.g, {hidden: true});
      K.rect(s.bandA, s.p.X(0), PY, s.p.X(fd[fd.length - 1]) - s.p.X(0), PH, {fill: "rgba(76,195,255,.09)"});
      K.text(s.bandA, s.p.X(fd[fd.length - 1] / 2), PY - 16, "every 0.5 Hz", {cls: "t-label", anchor: "middle", color: "var(--wave)"});
      s.ticks = tfu.f.map((f) => K.line(s.p.g, s.p.X(f), PY + PH, s.p.X(f), PY + PH - 28, {stroke: C.ink, "stroke-width": 4, hidden: true}));
      s.tf = s.p.line(tfi.f, tfi.amp, {color: C.ssi, width: 5, draw: true});
      s.dots = K.g(s.p.data, {hidden: true});
      tfu.f.forEach((f, i) => K.circle(s.dots, s.p.X(f), s.p.Y(tfu.amp[i]), 8, {fill: "#0a111d", stroke: C.ink, "stroke-width": 3}));
      const ip = argmax(tfi.amp), fp = tfi.f[ip], ap = tfi.amp[ip];
      const il = tfi.amp.findIndex((a, i) => i > ip && a < ap / 2);           // the curve's right flank, half height
      s.tfLab = s.p.text(tfi.f[il], tfi.amp[il], "roof transfer function", {cls: "t-label", dx: 26, dy: 6, color: "var(--ssi)", hidden: true});
      s.pk = K.g(s.p.g, {hidden: true});
      K.circle(s.pk, s.p.X(fp), s.p.Y(ap), 13, {fill: "none", stroke: C.ssi, "stroke-width": 4});
      K.text(s.pk, s.p.X(fp) + 28, s.p.Y(ap) + 10, `peak at ${fp.toFixed(1)} Hz`, {cls: "t-label", color: "var(--ssi)"});
      const fc = tfu.f[tfu.f.length - 1];
      s.cut = K.g(s.p.g, {hidden: true});
      K.rect(s.cut, s.p.X(fc), PY, s.p.X(25) - s.p.X(fc), PH, {fill: "rgba(255,107,107,.13)"});
      K.line(s.cut, s.p.X(fc), PY, s.p.X(fc), PY + PH, {stroke: C.bad, "stroke-width": 4});
      K.text(s.cut, s.p.X(22.5), PY + 60, "cut-off", {cls: "t-label", anchor: "middle", color: "var(--bad)", size: 32});
      K.text(s.cut, s.p.X(22.5), PY + 110, "nothing", {cls: "t-label", anchor: "middle", color: "var(--bad)"});
      K.text(s.cut, s.p.X(22.5), PY + 146, "gets through", {cls: "t-label", anchor: "middle", color: "var(--bad)"});
      s.chordOn = false;
    },
    tick(s, t) {
      if (!s.chordOn) return;
      const N = 120, x0 = 1345, w = 440, fr = [2, 3.3, 5.1], am = [1, 0.6, 0.42];
      const wave = (i, u) => am[i] * Math.sin(TAU * (fr[i] * u - 0.35 * t * (i + 1)));
      const c = [], ns = [[], [], []];
      for (let j = 0; j <= N; j++) {
        const u = j / N, x = x0 + w * u;
        c.push([x, 268 - 22 * (wave(0, u) + wave(1, u) + wave(2, u))]);
        for (let i = 0; i < 3; i++) ns[i].push([x, 402 + 46 * i - 17 * wave(i, u) / am[i]]);
      }
      s.chord.setAttribute("d", K.d(c));
      s.notes.forEach((p, i) => p.setAttribute("d", K.d(ns[i])));
    },
    beats: [
      {say: "[c]Next, the frequencies. Remember: SASSI works one frequency at a time, like taking a chord apart into notes.",
        go(k) { k.show(k.s.head); },
        c(k) { k.s.chordOn = true; k.show(k.s.ch); }},
      {say: "[t]So you choose the notes: here, [n]22 frequencies, from 0.1 to 20 Hz.",
        t(k) { k.show(k.s.p.g); k.show(k.s.ticks, {stagger: 45, delay: 200}); },
        n(k) { k.show(k.s.n22.el); }},
      {say: "[s]SASSI solves the problem at each one, [l]then fills in a smooth curve between them.",
        s(k) { k.show(k.s.dots); }, l(k) { k.draw(k.s.tf, 1000); }},
      {say: "[r]That's the roof's transfer function: how much the roof amplifies the ground motion at each frequency.",
        r(k) { k.show(k.s.tfLab); }},
      {say: "[p]Its peak sits near 3.5 Hz, right where the points are close together.",
        p(k) { k.show([k.s.pk, k.s.bandA]); }},
      {say: "[k]And the last frequency is a ceiling: nothing above 20 Hz gets through.",
        k(k) { k.show(k.s.cut); }},
    ],
  });

  // ---------------------------------------------------------------- 4. run the chain
  scenes.push({
    id: "run", title: "Run the chain",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Run the chain", {x: 120, y: 70, size: "h2"});
      const BW = 370, BH = 160, GAP = 50, X0 = 145, Y1 = 260;
      const row1 = [["SITE", "the ground shakes"], ["POINT", "the soil gives"], ["HOUSE", "your structure"], ["ANALYS", "joins the two"]];
      s.b1 = row1.map(([n, c], i) => modBox(s.root, X0 + i * (BW + GAP), Y1, BW, BH, n, c));
      s.a1 = [0, 1, 2].map((i) => rowArrow(s.root, X0 + i * (BW + GAP) + BW + 4, Y1 + BH / 2 - 18, GAP - 8));
      s.freq = K.label(s.root, "one frequency at a time: 22 solutions", {x: 960, y: Y1 + BH + 40, w: 600, align: "right", size: 28, color: "var(--ssi)", in: "fade"});
      // connector from ANALYS down to the post-processors
      const ax = X0 + 3 * (BW + GAP) + BW / 2, B2W = 420, B2G = 60, X2 = 960 - (3 * B2W + 2 * B2G) / 2, Y2 = 740, yc = 680;
      s.conn = K.path(g, `M${ax},${Y1 + BH}L${ax},${yc}L${X2 + B2W / 2},${yc}`, {stroke: "#5f7da6", "stroke-width": 4, fill: "none", draw: true, "stroke-linejoin": "round"});
      s.drops = K.g(g, {hidden: true});
      [0, 1, 2].forEach((i) => K.arrow(s.drops, X2 + B2W / 2 + i * (B2W + B2G), yc, X2 + B2W / 2 + i * (B2W + B2G), Y2 - 4, {color: "#5f7da6", width: 4, head: 16}));
      const row2 = [["MOTION", "floor spectra"], ["RELDISP", "displacements"], ["STRESS", "member forces"]];
      s.b2 = row2.map(([n, c], i) => modBox(s.root, X2 + i * (B2W + B2G), Y2, B2W, BH, n, c, "var(--ssi)"));
      s.run = 0; s.done = false;
    },
    tick(s, t) {
      const on = s.run ? Math.floor(t / 0.3) % 4 : -1;
      s.b1.forEach((b, i) => {
        b.style.borderColor = i === on ? "var(--ssi)" : s.done ? "var(--good)" : "#33507a";
        b.style.boxShadow = i === on ? "0 0 34px rgba(255,181,71,.5)" : "none";
      });
    },
    beats: [
      {say: "[c]Then press run. Four modules solve the problem, each with one job.",
        go(k) { k.show(k.s.head); },
        c(k) { const s = k.s; k.show([s.b1[0], s.a1[0], s.b1[1], s.a1[1], s.b1[2], s.a1[2], s.b1[3]], {stagger: 110}); }},
      {say: "[a]SITE works out how the ground shakes, [b]POINT how the soil gives under a push.",
        a(k) { k.show(k.s.b1[0].cap); k.pulse(k.s.b1[0], {amp: 0.05}); },
        b(k) { k.show(k.s.b1[1].cap); k.pulse(k.s.b1[1], {amp: 0.05}); }},
      {say: "[c]HOUSE builds your structure, [d]and ANALYS joins the two, one frequency at a time.",
        c(k) { k.show(k.s.b1[2].cap); k.pulse(k.s.b1[2], {amp: 0.05}); },
        d(k) { k.show(k.s.b1[3].cap); k.show(k.s.freq, {delay: 300}); k.s.run = 1; }},
      {say: "[p]Then three more modules turn the answer into design results.",
        go(k) { k.s.run = 0; k.s.done = true; },
        p(k) { k.draw(k.s.conn, 600); k.show(k.s.drops, {delay: 300}); k.show(k.s.b2, {stagger: 150, delay: 400}); k.show(k.s.b2.map((b) => b.cap), {stagger: 150, delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 5. the design results
  scenes.push({
    id: "results", title: "Read the results",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Read the results", {x: 120, y: 70, size: "h2"});
      // the stick to scale (7.5 px per foot: 64 ft mat, 16 ft storeys, 5 ft mat)
      const SX = 400, SY = 830, SC = 7.5;
      s.pic = K.g(g, {hidden: true});
      K.soil(s.pic, {x: 130, y: SY + 5 * SC, w: 540, layers: [{h: 16 * SC, kind: "sand"}], labels: false});
      s.axis = K.line(s.pic, SX, SY, SX, SY - 64 * SC - 20, {stroke: "rgba(147,161,179,.55)", "stroke-width": 3, "stroke-dasharray": "10 8"});
      s.stick = K.stick(s.pic, {x: SX, y: SY, z: [16, 32, 48, 64].map((z) => z * SC), matW: 64 * SC, matH: 5 * SC, r: 22, slabW: 16 * SC});
      const yR = SY - 64 * SC;
      s.equip = K.g(g, {hidden: true, in: "pop"});
      K.rect(s.equip, SX + 64, yR - 48, 58, 40, {rx: 6, fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.equip, SX + 136, yR - 18, "equipment", {cls: "t-label", color: "var(--ssi)"});
      // RELDISP: X displacement relative to the mat centre (ft), all floors at their maximum together (t = 1.86 s)
      const dft = [0.0155, 0.0349, 0.0559, 0.0766], XR = 200, PXFT = XR * SC;
      s.dft = dft; s.PXFT = PXFT; s.def = 0;
      s.dLab = K.g(g, {hidden: true});
      K.text(s.dLab, SX + dft[3] * PXFT + 16 * SC / 2 + 18, yR + 10, `${(12 * dft[3]).toFixed(2)} in`, {cls: "t-label", color: "var(--ssi)", size: 34});
      K.text(s.dLab, SX, 190, `relative to the mat, t = 1.86 s, displacements × ${XR}`, {cls: "t-small", anchor: "middle"});
      // STRESS: the shear at the base of the stick (6,764 kips, lesson 4)
      s.force = K.g(g, {hidden: true, in: "left"});
      K.arrow(s.force, SX - 220, SY - 50, SX - 26, SY - 50, {color: C.dash, width: 7, head: 24});
      K.text(s.force, SX - 220, SY - 74, "6,800 kips", {cls: "t-label", color: "var(--dash)", size: 34});
      // one plot at a time (right), annotated at the data's own maxima
      const PX = 930, PY = 250, PW = 850, PH = 400;
      s.ta = K.label(s.root, "<span style='font-family:var(--mono)'>MOTION</span> · floor spectra", {x: PX - 30, y: 160, size: 36, color: "var(--ssi)", in: "fade"});
      s.tb = K.label(s.root, "<span style='font-family:var(--mono)'>RELDISP</span> · roof vs mat", {x: PX - 30, y: 160, size: 36, color: "var(--ssi)", in: "fade"});
      s.tc = K.label(s.root, "<span style='font-family:var(--mono)'>STRESS</span> · shear at the base", {x: PX - 30, y: 160, size: 36, color: "var(--dash)", in: "fade"});
      const rs = F("00085TR_X02.RS"), ir = argmax(rs.sa);
      s.pa = K.plot(s.svg, {x: PX, y: PY, w: PW, h: PH, xr: [0.1, 100], xlog: true, yr: [0, 9], xticks: [0.1, 1, 10, 100], yticks: [0, 3, 6, 9],
        xlabel: "equipment frequency (Hz)", ylabel: "roof, 5 % damping (g)", ylabelOffset: 70});
      s.rs = s.pa.line(rs.f, rs.sa, {color: C.ssi, width: 5, draw: true});
      s.rsPk = K.g(s.pa.g, {hidden: true});
      K.circle(s.rsPk, s.pa.X(rs.f[ir]), s.pa.Y(rs.sa[ir]), 11, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.rsPk, s.pa.X(rs.f[ir]) + 26, s.pa.Y(rs.sa[ir]) + 12, `${rs.sa[ir].toFixed(1)} g at ${rs.f[ir].toFixed(1)} Hz`, {cls: "t-label", color: "var(--ssi)", size: 32});
      const thd = F("00085TR_X.THD"), id = argmax(thd.v);
      s.pb = K.plot(s.svg, {x: PX, y: PY, w: PW, h: PH, xr: [0, 24], yr: [-1, 1], xticks: [0, 6, 12, 18, 24], yticks: [-1, -0.5, 0, 0.5, 1],
        xlabel: "time (s)", ylabel: "roof vs mat (in)", ylabelOffset: 70});
      s.thd = s.pb.line(thd.t, thd.v.map((v) => 12 * v), {color: C.ssi, width: 3, draw: true});
      s.thdPk = s.pb.dot(thd.t[id], 12 * thd.v[id], {color: C.ssi, hidden: true, r: 11});
      const fyi = F("BEAMS_002_00001_FYI.THS"), iy = argmax(fyi.v);
      s.pc = K.plot(s.svg, {x: PX, y: PY, w: PW, h: PH, xr: [0, 24], yr: [-7200, 7200], xticks: [0, 6, 12, 18, 24], yticks: [-6000, 0, 6000],
        xlabel: "time (s)", ylabel: "base shear (kips)", ylabelOffset: 112});
      s.fyi = s.pc.line(fyi.t, fyi.v, {color: C.dash, width: 3, draw: true});
      s.fyiPk = s.pc.dot(fyi.t[iy], fyi.v[iy], {color: C.dash, hidden: true, r: 11});
    },
    tick(s) { s.stick.set({bend: s.dft.map((d) => d * s.PXFT * s.def)}); },
    beats: [
      {say: "[m]MOTION gives the floor spectra: how hard equipment on each floor gets shaken.",
        go(k) { k.show([k.s.head, k.s.pic]); },
        m(k) { k.show([k.s.ta, k.s.pa.g]); k.draw(k.s.rs, 1000, {delay: 200}); k.show(k.s.equip, {delay: 400}); }},
      {say: "[r]At the roof, the spectrum peaks at 8 g, [f]right at the new frequency, 3.5 Hz.",
        r(k) { k.show(k.s.rsPk); }, f(k) { k.pulse(k.s.rsPk, {amp: 0.08}); }},
      {say: "[d]RELDISP gives how far each floor moves relative to the mat: [n]0.92 inches at the roof.",
        d(k) {
          const s = k.s;
          k.hide([s.ta, s.pa.g, s.equip]); k.show([s.tb, s.pb.g], {delay: 200}); k.draw(s.thd, 1000, {delay: 300});
          k.tween(s, {def: 1}, 1000, {delay: 300});
        },
        n(k) { k.show([k.s.dLab, k.s.thdPk]); }},
      {say: "[f]STRESS gives the member forces: [n]6,800 kips of shear at the base of the stick.", hold: 600,
        f(k) { const s = k.s; k.hide([s.tb, s.pb.g]); k.show([s.tc, s.pc.g], {delay: 200}); k.draw(s.fyi, 1000, {delay: 300}); },
        n(k) { k.show([k.s.force, k.s.fyiPk]); }},
    ],
  });

  // ---------------------------------------------------------------- 6. why 3.5 Hz: sliding and rocking
  scenes.push({
    id: "springs", title: "Why 3.5 Hz?",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Why 3.5 Hz?", {x: 120, y: 70, size: "h2"});
      s.r = rig(g, {x: 500, y: 520, s: 5, soilW: 700, hsH: 80});
      s.cap = K.text(g, 500, 965, CAP, {cls: "t-small", anchor: "middle", hidden: true});
      s.sch = K.text(g, 500, 928, "spring and dashpot symbols: schematic", {cls: "t-small", anchor: "middle", hidden: true});
      const yR = 520 - 64 * 5;
      s.h20 = K.dim(g, 740, 520, 740, yR, "64 ft", {hidden: true, labelColor: "var(--ink)"});
      // which part of the computed motion is shown
      s.pS = K.pill(s.root, "the mat slides", {x: 1000, y: 180, color: "var(--spring)", style: {fontSize: "32px"}});
      s.pR = K.pill(s.root, "the mat rocks", {x: 1340, y: 180, color: "var(--spring)", style: {fontSize: "32px"}});
      // roof movement per push: sliding against rocking (the run's FOUNSTIF at 3.49 Hz)
      const BX = 1180, BW = 560, U = BW / (FLEX.rock / FLEX.slide);
      s.bars = K.g(g, {hidden: true});
      K.text(s.bars, 1000, 330, "roof movement per push", {cls: "t-label", size: 32, color: "var(--ink)"});
      K.text(s.bars, 1000, 402, "sliding", {cls: "t-label"});
      K.text(s.bars, 1000, 492, "rocking", {cls: "t-label"});
      K.text(s.bars, 1000, 562, "from SASSI's foundation springs at 3.5 Hz", {cls: "t-small", size: 24});
      s.barA = K.rect(g, BX, 372, 0, 40, {rx: 6, fill: C.spring, "fill-opacity": 0.55});
      s.barB = K.rect(g, BX, 462, 0, 40, {rx: 6, fill: C.spring});
      s.barALab = K.text(g, BX + U + 16, 402, "1 ×", {cls: "t-label", hidden: true, color: "var(--ink)"});
      s.barBLab = K.text(g, BX + BW - 14, 494, `${(FLEX.rock / FLEX.slide).toFixed(1)} ×`, {cls: "t-label", anchor: "end", size: 34, color: "#0a111d", hidden: true});
      s.U = U; s.BW = BW; s.bw = {a: 0, b: 0};
      s.res = K.h("div", {x: 1000, y: 650, w: 780, size: 64, in: "up",
        html: "<b style='color:var(--ref)'>5 Hz</b> <span style='color:var(--muted)'>→</span> <b style='color:var(--ssi)'>3.5 Hz</b>"}, s.root);
      s.resL = K.label(s.root, "mostly from rocking", {x: 1002, y: 750, size: 32, color: "var(--muted)", in: "fade"});
      s.w = {slide: 0, rock: 0, bend: 0, ground: 0};
    },
    tick(s, t) {
      pose(s.r, t, s.w);
      s.barA.setAttribute("width", Math.max(0, s.U * s.bw.a));
      s.barB.setAttribute("width", Math.max(0, s.BW * s.bw.b));
    },
    beats: [
      {say: "So why does the building sway at 3.5 Hz, and not 5?",
        go(k) { k.show(k.s.head); k.dim(k.s.r.soil.g, true); }},
      {say: "[s]Remember: the soil acts like a car's suspension, with springs, [d]and shock absorbers.",
        s(k) { k.show([k.s.r.springs, k.s.sch]); }, d(k) { k.show(k.s.r.dashes); }},
      {say: "[w]Springs make the whole system softer, so it sways more slowly.",
        w(k) { k.tween(k.s.w, {slide: 1, rock: 1, bend: 1, ground: 1}, 1000); k.show(k.s.cap); }},
      {say: "[a]On its springs, the mat can do two things: [l]slide sideways, [r]or rock like a seesaw.",
        a(k) { k.tween(k.s.w, {bend: 0}, 800); k.text(k.s.cap, "computed at 3.49 Hz, one part at a time: × 400, slowed down 7 ×"); },
        l(k) { k.show(k.s.pS); k.tween(k.s.w, {slide: 1, rock: 0}, 800); },
        r(k) { k.show(k.s.pR); k.dim(k.s.pS, true); k.tween(k.s.w, {slide: 0, rock: 1}, 800); }},
      {say: "Which matters more? [b]Push on the roof, and compare how far the soil lets it move.",
        b(k) { k.dim(k.s.pS, false); k.show(k.s.bars); k.tween(k.s.bw, {a: 1}, 800); k.show(k.s.barALab, {delay: 400}); }},
      {say: "[x]Rocking lets the roof move 3.5 times more than sliding does.",
        x(k) { k.tween(k.s.bw, {b: 1}, 1000); k.show(k.s.barBLab, {delay: 400}); }},
      {say: "[h]The roof is 64 ft up, so a small tilt at the base becomes a big sway at the top.",
        h(k) { k.show(k.s.h20); k.pulse(k.s.pR, {amp: 0.06}); }},
      {say: "[f]So rocking does most of the softening: that's what takes 5 Hz down to 3.5.", hold: 600,
        f(k) {
          k.tween(k.s.w, {slide: 1, rock: 1, bend: 1}, 1000); k.text(k.s.cap, CAP);
          k.show(k.s.res); k.show(k.s.resL, {delay: 400});
        }},
    ],
  });

  // ---------------------------------------------------------------- 7. why the peak stays high: radiation damping
  scenes.push({
    id: "damping", title: "Why the peak stays high",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Why the peak stays high", {x: 120, y: 70, size: "h2"});
      // the computed motion on the site, waves spreading into the ground (the waves are schematic)
      s.pic = K.g(g, {hidden: true});
      s.r = rig(s.pic, {x: 460, y: 520, s: 5, soilW: 680, hsH: 80});
      s.rip = K.ripples(s.pic, {cx: 460, cy: s.r.ys, r0: 40, r1: 330, n: 4, color: C.dash, width: 5, squash: 0.8});
      s.pic.insertBefore(s.rip.g, s.r.springs);
      K.text(s.pic, 790, s.r.ys + (16 + 39) * 5 + 80 - 27, "waves: schematic", {cls: "t-small", anchor: "end"});   // near the bottom of the soil
      s.cap = K.text(s.pic, 460, 965, CAP, {cls: "t-small", anchor: "middle"});
      s.radLab = K.text(g, 460, 928, "radiation damping", {cls: "t-label", anchor: "middle", color: "var(--dash)", hidden: true});
      s.ra = 0;
      s.w = {slide: 1, rock: 1, bend: 1, ground: 1};
      // FOUNDAMP of the run: damping ratio of sliding and rocking against frequency
      const PX = 980, PY = 190, PW = 800, PH = 440, YT = 25;
      s.p = K.plot(s.svg, {x: PX, y: PY, w: PW, h: PH, xr: [0, 12], yr: [0, YT], xticks: [0, 2, 4, 6, 8, 10, 12], yticks: [0, 5, 10, 15, 20, 25],
        xlabel: "frequency (Hz)", ylabel: "damping ratio (%)", ylabelOffset: 66});
      const dx = imp.dx.map((v) => 100 * v), dy = imp.dyy.map((v) => 100 * v);
      s.lx = s.p.line(imp.f, dx, {color: C.dash, width: 5, draw: true});
      s.ly = s.p.line(imp.f, dy, {color: C.violet, width: 5, draw: true});
      s.leg = s.p.legend([{label: "sliding", color: C.dash}, {label: "rocking", color: C.violet}], {x: PX + 30, y: PY + 36, size: 28, hidden: true});
      // where sliding leaves the frame: its value at 12 Hz, from the data
      const i12 = iF(12), ix = dx.findIndex((v) => v > YT), fx = imp.f[ix - 1] + (imp.f[ix] - imp.f[ix - 1]) * (YT - dx[ix - 1]) / (dx[ix] - dx[ix - 1]);
      s.off = s.p.text(fx, YT, `↑ ${Math.round(dx[i12])} % at 12 Hz`, {cls: "t-small", dx: 12, dy: 26, color: "var(--dash)", hidden: true});
      s.own = s.p.hline(dy[0], {color: C.ref, hidden: true});
      s.ownLab = s.p.text(12, dy[0], "soil's own damping", {cls: "t-label", anchor: "end", dy: 38, color: "var(--ref)", hidden: true});
      const i3 = iF(3.49), i10 = iF(10);
      s.v35 = s.p.vline(imp.f[i3], {color: C.ssi, hidden: true});
      s.v35Lab = s.p.text(imp.f[i3], YT, "3.5 Hz", {cls: "t-label", anchor: "middle", dy: -14, color: "var(--ssi)", hidden: true});
      s.rk = K.g(s.p.g, {hidden: true});
      K.circle(s.rk, s.p.X(imp.f[i3]), s.p.Y(dy[i3]), 12, {fill: "none", stroke: C.violet, "stroke-width": 4});
      K.text(s.rk, s.p.X(imp.f[i3]) + 22, s.p.Y(dy[i3]) + 66, `rocking ${dy[i3].toFixed(1)} %`, {cls: "t-label", color: "var(--violet)"});
      s.r10 = K.g(s.p.g, {hidden: true});
      K.circle(s.r10, s.p.X(imp.f[i10]), s.p.Y(dy[i10]), 12, {fill: "none", stroke: C.violet, "stroke-width": 4});
      K.text(s.r10, s.p.X(imp.f[i10]) - 22, s.p.Y(dy[i10]) - 14, "climbs above 10 Hz", {cls: "t-label", anchor: "end", color: "var(--violet)"});
      s.so = K.card(s.root, {x: PX, y: 770, w: PW, title: "So", size: 32, body: "Little energy escapes: the peak stays high, at a lower frequency."});
    },
    tick(s, t) {
      pose(s.r, t, s.w);
      s.rip.amp = s.ra; s.rip.set(t, TV);
    },
    beats: [
      {say: "Lesson 1 showed the roof peak barely dropped. [q]So why didn't the soil calm it down?",
        go(k) { k.show([k.s.head, k.s.pic]); }},
      {say: "[r]The shock absorbers work by sending waves away into the ground, like ripples from a stone in a pond.",
        r(k) { k.tween(k.s, {ra: 1}, 1000); }},
      {say: "[e]Engineers call that radiation damping. [p]SASSI computes it for each motion, frequency by frequency.",
        e(k) { k.show(k.s.radLab); },
        p(k) { k.show([k.s.p.g, k.s.leg]); k.draw(k.s.lx, 1000); k.draw(k.s.ly, 1000, {delay: 200}); k.show(k.s.off, {delay: 400}); }},
      {say: "[a]At 3.5 Hz, rocking gets about 5 percent, [m]barely more than the soil's own damping.",
        a(k) { k.show([k.s.v35, k.s.v35Lab, k.s.rk]); },
        m(k) { k.show([k.s.own, k.s.ownLab]); }},
      {say: "[t]Why so little? On a layered site, at low frequency, the waves can hardly get away.",
        t(k) { k.tween(k.s, {ra: 0.15}, 1000); k.text(k.s.radLab, "waves can hardly get away"); k.show(k.s.r10, {delay: 400}); }},
      {say: "[k]So the system stays lightly damped: the peak stays high, it just moves to a lower frequency.", hold: 600,
        k(k) { k.show(k.s.so); }},
    ],
  });

  // ---------------------------------------------------------------- 8. recap
  scenes.push({
    id: "recap", title: "Recap",
    build(s) {
      s.head = K.heading(s.root, "Recap", {x: 120, y: 80, size: "h1"});
      s.list = K.bullets(s.root, [
        {t: "SSI model = your ANSYS model", sub: "plus a frequency list and the soil nodes"},
        {t: "The SASSI chain gives design results", sub: "floor spectra, displacements, member forces"},
        {t: "The mat rocks: 5 → 3.5 Hz", sub: "little energy escapes, so the peak stays high"},
      ], {x: 120, y: 230, w: 1000, num: true, size: 44});
      s.q = K.card(s.root, {x: 1180, y: 200, w: 640, kind: "check", title: "Check yourself", size: 34,
        body: "A roof component is tuned to 3.6 Hz. Did SSI raise or lower its demand?"});
      s.a = K.card(s.root, {x: 1180, y: 540, w: 640, title: "Answer", size: 34,
        body: "<b>Raised, about 3 ×:</b> 7.9&nbsp;g with SSI, 2.8&nbsp;g on a fixed base."});
      s.next = K.pill(s.root, "Next · Lesson 5: embedded structures and the flexible-volume method →", {x: 120, y: 880, color: "var(--wave)", style: {fontSize: "30px"}});
    },
    beats: [
      {say: "[a]One: an SSI model is your ANSYS model, plus a frequency list and the nodes that touch the soil.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[b]Two: SASSI's chain turns it into floor spectra, displacements and forces, with the soil included.",
        b(k) { k.show(k.s.list.items[1]); }},
      {say: "[c]Three: the mat rocks on the soil, so the frequency drops; little energy escapes, so the peak stays high.",
        c(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Check yourself: a roof component is tuned to 3.6 Hz. Did SSI raise or lower its demand?", gap: 1500,
        q(k) { k.show(k.s.q); }},
      {say: "[a]It roughly tripled: 7.9 g with SSI, against 2.8 g on a fixed base.",
        a(k) { k.show(k.s.a); }},
      {say: "[n]Next: buildings with basements, and how SASSI handles the soil you dig out.",
        n(k) { k.show(k.s.next); }},
    ],
  });

  SV.video({
    id: "04", n: 4, part: "Fundamentals", lesson: "04-surface-ssi",
    title: "Your first SSI analysis: a stick on a surface mat",
    subtitle: "Build example 1 like an ANSYS model, choose the frequencies, run the chain, read the results, and see why it sways at 3.5 Hz instead of 5.",
    next: {href: "05.html", title: "Embedded structures and the flexible-volume method"},
    pron: [
      [/\bSITE\b/g, "site"], [/\bPOINT\b/g, "point"], [/\bHOUSE\b/g, "house"], [/\bMOTION\b/g, "motion"], [/\bSTRESS\b/g, "stress"],
      [/\bINT\b/g, "I N T"],
    ],
    scenes,
  });
})();
