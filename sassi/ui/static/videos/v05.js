/* Explainer video, lesson 5 (newcomer edition): Embedded structures and the flexible-volume method.
 * The big idea: a basement is buried, so the ground motion changes with depth; SASSI adds the basement and
 * subtracts the soil it replaces (flexible volume, FV); shortcuts with fewer interaction nodes can invent a fake
 * resonance, so they are checked against FV; a stiff basement averages the motion (kinematic interaction).
 * Numbers: sassi/ui/lessons/05_embedded.md; curves: d05.js (the lesson's runs, python -m sassi.ui.video_data);
 * the free field at depth: SV.K.PH.columnWaves (the SHAKE recursion of SITE) for the lesson's site, which
 * reproduces the SITE listing (0.81, 0.32, 0.10 and 0.82 at 16 ft depth).
 * Units: ft, kip, s.  Pictures to scale (the 32 ft x 32 ft x 16 ft box, 3.2 ft levels, 8 ft grid, 16 ft of sand);
 * every moving picture is the computed steady state for 0.04 in of ground motion with one displacement factor per
 * scene, stated on screen. */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["05"], PH = K.PH;
  const TAU = 2 * Math.PI;

  // ---------------------------------------------------------------- data helpers
  /** Complex value of a .TFU file at the computed frequency nearest f. */
  function tfAt(name, f) {
    const o = D[name];
    let i = 0;
    o.f.forEach((v, j) => { if (Math.abs(v - f) < Math.abs(o.f[i] - f)) i = j; });
    return [o.amp[i] * Math.cos(o.ph[i]), o.amp[i] * Math.sin(o.ph[i])];
  }
  const cabs = (z) => Math.hypot(z[0], z[1]);
  const lerp = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t];
  const reph = (z, c, s) => z[0] * c - z[1] * s;                 // Re[z e^{i phi}]
  /** The lesson's site (16 ft sand, 51.2 ft gravel, rock; ft, ft/s, kcf): the free field per unit surface motion
   *  at depth z (ft). */
  const COL = {g: 32.2, layers: [{h: 16, vs: 800, w: 0.120, beta: 0.05}, {h: 51.2, vs: 1300, w: 0.125, beta: 0.04}], hs: {vs: 2600, w: 0.130, beta: 0.02}};
  const freeField = (f) => { const wv = PH.columnWaves(COL, f); return (z) => PH.columnU(wv, COL, z); };
  const TV = 2.4;                                                   // visual period of every moving picture (s)
  const UG = 0.04 / 12;                                             // ground motion of every moving picture: 0.04 in, in ft
  /** Pixels per unit of control motion: 0.04 in of ground motion drawn x times larger, at pxf px per ft. */
  const ampPx = (x, pxf) => x * UG * pxf;
  /** The caption of a moving picture: frequency, 0.04 in of ground motion, the displacement factor, the slow-down. */
  const capOf = (f, x) => `computed at ${Math.round(f)} Hz: ground motion 0.04 in, displacements × ${x}, slowed down ${Math.round(f * TV)} ×`;

  // ---------------------------------------------------------------- drawing helpers
  /** Isometric view of the 32 x 32 x 16 ft excavation mesh (PH.boxNodes, ft): node groups by interaction set.
   *  o: {cx, cy, s (px per ft), ground}. Returns {g, P, on: {base, sides, top, even, odd}, off, all}. */
  function iso(g, o) {
    const a = 0.62, ca = Math.cos(a), sa = Math.sin(a), s = o.s;
    const P = (x, y, z) => [o.cx + s * (x * ca - y * sa), o.cy + s * (-z - 8 + (x * sa + y * ca) * 0.42)];
    const r = {g: K.g(g, {hidden: o.hidden}), P};
    const sq = (z, G) => [P(-G, -G, z), P(G, -G, z), P(G, G, z), P(-G, G, z)];
    if (o.ground !== false) K.path(r.g, K.d(sq(0, o.G || 20.8)) + "Z", {fill: "rgba(200,164,107,.20)", stroke: "rgba(200,164,107,.45)", "stroke-width": 2});
    for (let k = 0; k <= 5; k++) K.path(r.g, K.d(sq(-16 + 3.2 * k, 16)) + "Z", {fill: "none", stroke: k === 5 ? "rgba(223,230,238,.65)" : "rgba(130,160,200,.30)", "stroke-width": k === 5 ? 2.5 : 1.5});
    for (const [x, y] of [[-16, -16], [16, -16], [16, 16], [-16, 16]]) K.line(r.g, ...P(x, y, 0), ...P(x, y, -16), {stroke: "rgba(130,160,200,.45)", "stroke-width": 2});
    const nodes = PH.boxNodes();
    const sorted = nodes.slice().sort((u, v) => (u.x * sa + u.y * ca) - (v.x * sa + v.y * ca) || u.z - v.z);
    const rr = Math.max(4, s * 0.512);
    r.off = K.g(r.g, {});
    sorted.forEach((nd) => K.circle(r.off, ...P(nd.x, nd.y, nd.z), rr * 0.75, {fill: "#0a111d", stroke: C.ref, "stroke-width": 2}));
    r.on = {};
    for (const key of ["base", "sides", "odd", "even", "top"]) r.on[key] = K.g(r.g, {hidden: true});
    sorted.forEach((nd) => {
      const side = nd.i === 0 || nd.j === 0 || nd.i === 4 || nd.j === 4;
      const key = nd.k === 0 ? "base" : side ? "sides" : nd.k === 5 ? "top" : nd.k % 2 === 0 ? "even" : "odd";
      K.circle(r.on[key], ...P(nd.x, nd.y, nd.z), rr, {fill: o.color || C.wave, stroke: "#0a111d", "stroke-width": 1.5});
    });
    r.all = ["base", "sides", "odd", "even", "top"].map((key) => r.on[key]);
    return r;
  }
  /** Section y = 0 of the box (x-z plane) through its 5 x 6 nodes (8 ft apart across, 3.2 ft levels): walls, base
   *  slab, roof slab.  o: {x0 (centre px), y0 (grade px), sx, sz (px per ft), color, dot}. set(disp(i, k) -> [dx, dz] px). */
  function section(g, o) {
    const r = {o, g: K.g(g, {hidden: o.hidden})};
    const col = o.color || C.concrete;
    r.fill = K.path(r.g, "", {fill: "rgba(40,53,70,.95)", stroke: "none"});
    r.base = K.path(r.g, "", {stroke: col, "stroke-width": 11, fill: "none", "stroke-linecap": "round", "stroke-linejoin": "round"});
    r.walls = K.path(r.g, "", {stroke: col, "stroke-width": 7, fill: "none", "stroke-linecap": "round", "stroke-linejoin": "round"});
    r.roof = K.path(r.g, "", {stroke: col, "stroke-width": 7, fill: "none", "stroke-linecap": "round", "stroke-linejoin": "round"});
    r.dots = [];
    const bnd = [];
    for (let i = 0; i <= 4; i++) bnd.push([i, 0], [i, 5]);
    for (let k = 1; k <= 4; k++) bnd.push([0, k], [4, k]);
    bnd.forEach(([i, k]) => r.dots.push({i, k, el: K.circle(r.g, 0, 0, 6, {fill: o.dot || C.ink2, stroke: "#0a111d", "stroke-width": 2})}));
    r.set = (disp) => {
      const pt = (i, k) => { const d = disp ? disp(i, k) : [0, 0]; return [o.x0 + (i - 2) * 8 * o.sx + d[0], o.y0 + (5 - k) * 3.2 * o.sz - d[1]]; };
      const base = [0, 1, 2, 3, 4].map((i) => pt(i, 0)), roof = [0, 1, 2, 3, 4].map((i) => pt(i, 5));
      const left = [0, 1, 2, 3, 4, 5].map((k) => pt(0, k)), right = [0, 1, 2, 3, 4, 5].map((k) => pt(4, k));
      r.base.setAttribute("d", K.d(base));
      r.roof.setAttribute("d", K.d(roof));
      r.walls.setAttribute("d", K.d(left) + K.d(right));
      r.fill.setAttribute("d", K.d(base.concat(right.slice(1), roof.slice().reverse().slice(1), left.slice().reverse().slice(1))) + "Z");
      r.dots.forEach((d) => { const q = pt(d.i, d.k); d.el.setAttribute("cx", q[0]); d.el.setAttribute("cy", q[1]); });
    };
    r.set(null);
    return r;
  }
  /** Motion of the box from the run: X of the base slab (Hb) and roof (Hr) centres (the walls interpolated between
   *  them), Z of the roof edge (Hz, x = +16 ft) for the rocking. */
  const rigid = (Hb, Hr, Hz, A, c, s) => (i, k) => [A * reph(lerp(Hb, Hr, k / 5), c, s), A * reph(Hz, c, s) * (i - 2) / 2];
  /** Soil column with the free field at f: o {x, y (grade), w, pxf (px per ft), depth (ft), px (profile x)};
   *  set(c, s, A): A = px per unit of surface motion. */
  function column(g, o) {
    const r = {o, g: K.g(g, {hidden: o.hidden})};
    r.soil = K.soil(r.g, {x: o.x, y: o.y, w: o.w, labels: o.labels || false, labelSize: 24,
      layers: [{h: 16 * o.pxf, kind: "sand", name: "sand · 16 ft"}, {h: (o.depth - 16) * o.pxf, kind: "gravel", name: "gravel · 51.2 ft"}]});
    r.U = freeField(o.f);
    r.prof = K.g(g, {hidden: true});
    K.line(r.prof, o.px, o.y, o.px, o.y + o.depth * o.pxf, {stroke: "rgba(76,195,255,.35)", "stroke-width": 2, "stroke-dasharray": "8 6"});
    r.curve = K.path(r.prof, "", {stroke: C.wave, "stroke-width": 5, fill: "none", "stroke-linecap": "round"});
    r.zs = K.linspace(0, o.depth, 45);
    r.Uz = r.zs.map((z) => r.U(z));
    r.set = (c, s, A) => {
      r.soil.shear((d) => A * reph(r.U(d / o.pxf), c, s));
      r.curve.setAttribute("d", K.d(r.zs.map((z, j) => [o.px + A * reph(r.Uz[j], c, s), o.y + z * o.pxf]), true));
    };
    r.set(1, 0, 0);
    return r;
  }
  /** Line + dots of a .TFU file (the computed frequencies) on a plot. */
  function tfLine(p, name, o) {
    const t = D[name];
    const line = p.line(t.f, t.amp, Object.assign({width: 4, draw: true}, o));
    const dots = K.g(p.data, {hidden: true});
    t.f.forEach((f, i) => K.circle(dots, p.X(f), p.Y(t.amp[i]), o.r || 6, {fill: o.color, stroke: "#0a111d", "stroke-width": 2}));
    return {line, dots};
  }
  /** A hidden block inside a card body (revealed with k.show). */
  function cardPart(card, html) {
    const d = document.createElement("div");
    d.className = "sv-in sv-fade";
    d.style.marginTop = "8px";
    d.innerHTML = html;
    card.body.appendChild(d);
    return d;
  }
  /** The roof plot shared by the two shortcut scenes: roof X transfer functions, FV and FI-FSIN. */
  function roofPlot(s) {
    s.p = K.plot(s.svg, {x: 1000, y: 190, w: 800, h: 400, xr: [0, 20], yr: [0, 1.2], xticks: [0, 5, 10, 15, 20], yticks: [0, 0.4, 0.8, 1.2],
      xlabel: "frequency (Hz)", ylabel: "roof motion ÷ ground", ylabelOffset: 80});
    s.fv = tfLine(s.p, "ex02/00138TR_X.TFU", {color: C.ref, width: 9, opacity: 0.6, r: 7});
    s.fs = tfLine(s.p, "ex02_fsin/00138TR_X.TFU", {color: C.bad, r: 7});
  }

  // ================================================================== scenes
  const scenes = [];
  const argmax = (a) => a.reduce((b, v, i) => (Math.abs(v) > Math.abs(a[b]) ? i : b), 0);

  // ---------------------------------------------------------------- 0. hook
  scenes.push({
    id: "hook", title: "A buried box",
    build(s) {
      const g = K.g(s.svg);
      const Y0 = 340, PXF = 15.625, F12 = 12.0117;              // 15.625 px per ft: the 32 ft x 16 ft box is 500 x 250 px
      s.col = column(g, {x: 120, y: Y0, w: 760, pxf: PXF, depth: 37, f: F12, px: 790, hidden: true});
      K.text(s.col.g, 144, Y0 + 16 * PXF + 64, "gravel · 51.2 ft", {cls: "t-label", size: 24});
      s.box = K.rect(g, 420 - 16 * PXF, Y0, 32 * PXF, 16 * PXF, {fill: "rgba(10,17,29,.55)", stroke: C.concrete, "stroke-width": 4, "stroke-dasharray": "14 10", hidden: true});
      s.qm = K.text(g, 420, Y0 + 160, "?", {cls: "t-big", anchor: "middle", size: 120, color: "var(--concrete)", hidden: true});
      s.dim5 = K.dim(g, 145, Y0, 145, Y0 + 16 * PXF, "16 ft", {side: "left", hidden: true, labelColor: "var(--ink)"});
      const u5 = Math.hypot(...s.col.U(16)).toFixed(2);                // the free field 16 ft down (0.10)
      s.l1 = K.text(g, 950, Y0 + 16, "1", {cls: "t-label", size: 34, color: "var(--wave)", hidden: true});
      s.l2 = K.text(g, 950, Y0 + 16 * PXF + 12, u5, {cls: "t-label", size: 34, color: "var(--wave)", hidden: true});
      s.slow = K.text(g, 500, 965, "free field: " + capOf(F12, 1000), {cls: "t-small", anchor: "middle", hidden: true});
      s.s1 = K.stat(s.root, {x: 1060, y: 200, w: 340, value: "1", label: "ground surface", color: "var(--wave)", vsize: 120});
      s.s2 = K.stat(s.root, {x: 1430, y: 200, w: 340, value: u5, label: "16 ft down", color: "var(--wave)", vsize: 120});
      s.at = K.label(s.root, "at 12 Hz, on this site", {x: 1060, y: 410, w: 710, size: 30, color: "var(--muted)", in: "fade", style: {textAlign: "center"}});
      s.q = K.heading(s.root, '<span style="color:var(--ssi)">What does a buried box feel?</span>', {x: 1060, y: 610, w: 710, align: "center", size: "h1"});
      s.A = 0; s.AP = ampPx(1000, PXF);
    },
    tick(s, t) { const ph = TAU * t / TV; s.col.set(Math.cos(ph), Math.sin(ph), s.AP * s.A); },
    beats: [
      {say: "Most nuclear buildings sit partly underground, [b]with a basement below grade.",
        go(k) { k.show(k.s.col.g); }, b(k) { k.show([k.s.box, k.s.dim5]); }},
      {say: "[w]And underground, the ground doesn't shake the same at every depth.",
        w(k) { k.show([k.s.col.prof, k.s.slow]); k.tween(k.s, {A: 1}, 1000); }},
      {say: "On this site at 12 Hz, [a]the surface moves 1, [b]but 16 ft down, only 0.10.",
        a(k) { k.show([k.s.s1.el, k.s.l1, k.s.at]); }, b(k) { k.show([k.s.s2.el, k.s.l2]); }},
      {say: "[q]So what motion does a buried basement actually feel?", hold: 600,
        q(k) { k.show([k.s.q, k.s.qm]); }},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 5",
    build(s) {
      s.t = K.titleCard(s.root, {n: 5, part: "Fundamentals", title: "Embedded structures and the flexible-volume method",
        sub: "Dig the hole, subtract the soil, check the shortcuts, and see what a basement feels."});
      s.t.title.style.maxWidth = "1180px";
    },
    beats: [
      {say: "Lesson five. Embedded structures, and the flexible-volume method.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 100}); k.show(t.title, {delay: 200}); k.show(t.sub, {delay: 350}); k.show(t.rule, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. the flexible-volume idea (to scale: 4.375 px per ft)
  scenes.push({
    id: "flexible-volume", title: "Add the box, subtract the soil",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Add the box, subtract the soil", {x: 120, y: 70, size: "h2"});
      const PXF = 4.375, BW = 32 * PXF, BH = 16 * PXF, Y = 300, DZ = 3.2 * PXF;   // the 32 ft x 16 ft box: 140 x 70 px
      const panel = (x, title, color, sub) => {
        const pg = K.g(g, {hidden: true});
        K.rect(pg, x, 180, 340, 340, {rx: 18, fill: "rgba(20,32,51,.6)", stroke: "#2a3b52", "stroke-width": 2});
        K.text(pg, x + 170, 234, title, {cls: "t-label", anchor: "middle", color, size: 32});
        K.text(pg, x + 170, 566, sub, {cls: "t-label", anchor: "middle", color: "var(--ink2)"});
        return pg;
      };
      // 16 ft of sand (the embedment) over the gravel, cut at the panel's frame
      const soilP = (pg, x) => K.soil(pg, {x: x + 25, y: Y, w: 290, layers: [{h: 16 * PXF, kind: "sand"}, {h: 200 - 16 * PXF, kind: "gravel"}], labels: false});
      const boxP = (pg, cx) => {
        K.rect(pg, cx - BW / 2, Y, BW, BH, {fill: "rgba(40,53,70,.95)", stroke: C.concrete, "stroke-width": 5});
        K.line(pg, cx - BW / 2, Y + BH, cx + BW / 2, Y + BH, {stroke: C.concrete, "stroke-width": 9});
      };
      s.p1 = panel(110, "Free field", "var(--wave)", "no hole");
      s.ff = soilP(s.p1, 110);
      s.U = freeField(12.0117);
      s.hole = K.rect(s.p1, 280 - BW / 2, Y, BW, BH, {fill: "none", stroke: C.ink, "stroke-width": 3, "stroke-dasharray": "10 7", hidden: true});
      s.ffCap = K.g(g, {hidden: true});
      K.text(s.ffCap, 280, 604, "free field at 12 Hz, computed", {cls: "t-small", anchor: "middle", size: 22});
      K.text(s.ffCap, 280, 632, "0.04 in × 1000, slowed down 29 ×", {cls: "t-small", anchor: "middle", size: 22});
      s.p2 = panel(560, "Structure", "var(--steel)", "the basement");
      boxP(s.p2, 730);
      // its nodes on the section: 5 across (8 ft), 6 levels (3.2 ft)
      K.nodes(s.p2, [0, 1, 2, 3, 4].flatMap((i) => [[730 - BW / 2 + 0.25 * BW * i, Y], [730 - BW / 2 + 0.25 * BW * i, Y + BH]])
        .concat([1, 2, 3, 4].flatMap((k) => [[730 - BW / 2, Y + DZ * k], [730 + BW / 2, Y + DZ * k]])), {r: 5, color: "var(--ink2)"});
      s.p3 = panel(1010, "Excavated soil", "var(--sand)", "the soil dug out");
      K.rect(s.p3, 1180 - BW / 2, Y, BW, BH, {fill: "rgba(200,164,107,.8)", stroke: C.ink2, "stroke-width": 3, "stroke-dasharray": "10 7"});
      for (let i = 1; i < 4; i++) K.line(s.p3, 1180 - BW / 2 + 0.25 * BW * i, Y, 1180 - BW / 2 + 0.25 * BW * i, Y + BH, {stroke: "rgba(255,240,215,.5)", "stroke-width": 1.5});
      for (let k = 1; k < 5; k++) K.line(s.p3, 1180 - BW / 2, Y + DZ * k, 1180 + BW / 2, Y + DZ * k, {stroke: "rgba(255,240,215,.5)", "stroke-width": 1.5});
      s.p4 = panel(1460, "SSI system", "var(--ssi)", "what SASSI solves");
      soilP(s.p4, 1460); boxP(s.p4, 1630);
      s.ops = ["+", "−", "="].map((t, i) => K.text(g, [535, 985, 1435][i], 376, t, {cls: "t-big", anchor: "middle", size: 70, hidden: true}));
      s.eq = K.eq(s.root, String.raw`\text{SSI system} = \text{free field} + \text{structure} - \htmlClass{t-ex}{\text{excavated soil}}`,
        {x: 160, y: 650, w: 1600, size: 52, align: "center", display: false});
      s.fv = K.pill(s.root, "the flexible-volume method · FV", {x: 700, y: 800, color: "var(--ssi)", style: {fontSize: "36px"}});
      s.A = 0; s.PXF = PXF; s.AP = ampPx(1000, PXF);
    },
    tick(s, t) {
      const ph = TAU * t / TV, c = Math.cos(ph), sn = Math.sin(ph), A = s.AP * s.A;
      s.ff.shear((d) => A * reph(s.U(d / s.PXF), c, sn));
    },
    beats: [
      {say: "[a]It starts from the free field: remember, the ground shaking with no building on it, and no hole.",
        go(k) { k.show(k.s.head); },
        a(k) { k.show([k.s.p1, k.s.ffCap]); k.tween(k.s, {A: 1}, 1000); }},
      {say: "[b]Add the structure: the basement walls and slabs, modelled as you would in ANSYS.",
        b(k) { k.show([k.s.ops[0], k.s.p2], {stagger: 150}); }},
      {say: "[c]But the free field still contains the soil that the basement replaced.",
        c(k) { k.show(k.s.hole); k.pulse(k.s.hole, {amp: 0.08}); }},
      {say: "[c]So SASSI meshes that soil, the excavated soil, [d]only to subtract it.",
        c(k) { k.show([k.s.ops[1], k.s.p3], {stagger: 150}); }, d(k) { k.pulse(k.s.ops[1], {amp: 0.25}); }},
      {say: "[e]Free field, plus structure, minus excavated soil: [f]that's the SSI system.",
        e(k) { k.show(k.s.eq); }, f(k) { k.show([k.s.ops[2], k.s.p4], {stagger: 150}); }},
      {say: "[m]Engineers call it the flexible-volume method: FV for short.", hold: 600,
        m(k) { k.show(k.s.fv); }},
    ],
  });

  // ---------------------------------------------------------------- 3. FV: every excavated node holds on
  scenes.push({
    id: "nodes", title: "Where the soil holds on",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Where the soil holds on", {x: 120, y: 70, size: "h2"});
      s.iso = iso(g, {cx: 540, cy: 580, s: 16.875});
      s.dims = K.text(g, 540, 960, "the excavation: 32 ft × 32 ft, 16 ft deep, to scale", {cls: "t-label", anchor: "middle"});
      s.n = K.stat(s.root, {x: 1080, y: 190, w: 700, value: "0", label: "interaction nodes in FV", color: "var(--wave)", vsize: 130});
      s.exact = K.pill(s.root, "exact: the reference", {x: 1290, y: 420, color: "var(--good)", style: {fontSize: "32px"}});
      s.island = K.label(s.root, "a nuclear island: tens of thousands", {x: 1080, y: 590, w: 700, size: 32, color: "var(--ink)", in: "fade", style: {textAlign: "center"}});
      s.short = K.pill(s.root, "so: shortcuts with fewer nodes", {x: 1225, y: 680, color: "var(--ssi)", style: {fontSize: "32px"}});
    },
    beats: [
      {say: "[i]Where do the basement and the soil hold hands? In FV, at every node of the dug-out volume.",
        go(k) { k.show(k.s.head); },
        i(k) { k.show(k.s.iso.all, {stagger: 150}); k.show(k.s.dims, {delay: 400}); }},
      {say: "[n]That's 150 interaction nodes, for a 32-foot box, 16 ft deep.",
        n(k) { k.show(k.s.n.el); k.count(k.s.n.v, PH.setCount("fv"), {from: 0, ms: 1000}); }},
      {say: "[x]Each one feels the soil, and the free-field motion at its own depth: that's the exact answer.",
        x(k) { k.show(k.s.exact); k.pulse(k.s.iso.all, {amp: 0.03, ms: 900}); }},
      {say: "[r]A real nuclear island can have tens of thousands. [s]So engineers use shortcuts, with fewer nodes.", hold: 600,
        r(k) { k.show(k.s.island); }, s(k) { k.show(k.s.short); }},
    ],
  });

  // ---------------------------------------------------------------- 4. the subtraction method: a fake resonance
  scenes.push({
    id: "shortcut", title: "A shortcut that misleads",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "A shortcut that misleads", {x: 120, y: 70, size: "h2"});
      s.iso = iso(g, {cx: 470, cy: 530, s: 13.75});
      s.isoT = K.text(g, 470, 230, "walls and floor only: FI-FSIN", {cls: "t-label", anchor: "middle", size: 32, color: "var(--bad)", hidden: true});
      s.cnt = K.stat(s.root, {x: 270, y: 800, w: 400, value: String(PH.setCount("fv")), label: "interaction nodes", color: "var(--wave)", vsize: 100});
      roofPlot(s);
      s.leg = s.p.legend([{label: "FV, 150 nodes", color: C.ref}, {label: "shortcut, 105 nodes", color: C.bad}], {x: 1040, y: 430, size: 28, dy: 44, hidden: true});
      // the shortcut's fake peak: its largest excess over FV at the same frequency (the computed points, 11 Hz)
      const fs = D["ex02_fsin/00138TR_X.TFU"], fv = D["ex02/00138TR_X.TFU"], ip = argmax(fs.amp.map((a, i) => a / fv.amp[i]));
      s.pk = K.g(s.p.g, {hidden: true});
      K.circle(s.pk, s.p.X(fs.f[ip]), s.p.Y(fs.amp[ip]), 16, {fill: "none", stroke: C.bad, "stroke-width": 4});
      s.pkL = s.p.text(fs.f[ip], fs.amp[ip], `+${Math.round(100 * (fs.amp[ip] / fv.amp[ip] - 1))} %`, {cls: "t-label", dx: 30, dy: 12, size: 40, color: "var(--bad)", hidden: true});
      const m1 = D["ex02/SHELL_006_00074_MYY.THS"].peak, m2 = D["ex02_fsin/SHELL_006_00074_MYY.THS"].peak;
      s.mom = K.h("div", {x: 1000, y: 690, w: 800, size: 40, in: "up",
        html: `wall moment: <b>${m1.toFixed(2)}</b> <span style='color:var(--muted)'>→</span> <b style='color:var(--bad)'>${m2.toFixed(2)}</b> <span style='font-size:28px;color:var(--muted)'>kip·ft/ft</span> <b style='color:var(--bad)'>+${Math.round(100 * (m2 / m1 - 1))} %</b>`}, s.root);
    },
    beats: [
      {say: "[s]The classic shortcut keeps only the nodes on the walls and the floor: [n]105 instead of 150.",
        go(k) { k.show([k.s.head, k.s.iso.all, k.s.cnt.el]); },
        s(k) { k.hide([k.s.iso.on.top, k.s.iso.on.even, k.s.iso.on.odd]); k.show(k.s.isoT); },
        n(k) { k.count(k.s.cnt.v, PH.setCount("fsin"), {from: PH.setCount("fv"), ms: 1000}); }},
      {say: "[a]At low frequency, the two agree. [k]Then the shortcut shows a peak near 11 Hz.",
        a(k) {
          const s = k.s;
          k.show([s.p.g, s.leg]); k.draw(s.fv.line, 1000); k.draw(s.fs.line, 1000, {delay: 200}); k.show([s.fv.dots, s.fs.dots], {delay: 400});
        },
        k(k) { k.show(k.s.pk); k.pulse(k.s.pk, {amp: 0.2}); }},
      {say: "[v]The roof moves 40 percent too much, [m]and the wall bending moment rises 30 percent.", hold: 600,
        v(k) { k.show(k.s.pkL); }, m(k) { k.show(k.s.mom); }},
    ],
  });

  // ---------------------------------------------------------------- 5. why, the fix, and the rule
  scenes.push({
    id: "check", title: "Check it against FV",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Check it against FV", {x: 120, y: 70, size: "h2"});
      // the two models in section y = 0, to scale (9.375 px per ft), in the free field at 11 Hz
      const PXF = 9.375, XD = 150, F11 = 10.986, Y0 = 330;
      s.U = freeField(F11); s.PXF = PXF; s.AX = ampPx(XD, PXF);
      const band = (x) => K.soil(g, {x: x - 170, y: Y0, w: 340, layers: [{h: 16 * PXF, kind: "sand"}, {h: 9.6 * PXF, kind: "gravel"}], labels: false});
      s.bands = [band(280), band(690)];
      K.text(g, 280, 232, "FV", {cls: "t-label", anchor: "middle", color: "var(--ref)", size: 34});
      K.text(g, 690, 232, "shortcut", {cls: "t-label", anchor: "middle", color: "var(--bad)", size: 34});
      s.fvS = section(g, {x0: 280, y0: Y0, sx: PXF, sz: PXF});
      s.fsS = section(g, {x0: 690, y0: Y0, sx: PXF, sz: PXF, dot: C.bad});
      s.fvL = K.text(g, 280, 612, "follows the ground", {cls: "t-label", anchor: "middle", color: "var(--ref)", hidden: true});
      s.fsL = K.text(g, 690, 612, "roof slab: over 15 ×", {cls: "t-label", anchor: "middle", color: "var(--bad)", hidden: true});
      s.note = K.g(g, {hidden: true});
      K.text(s.note, 485, 652, capOf(F11, XD), {cls: "t-small", anchor: "middle", size: 22});
      K.text(s.note, 485, 682, "soil inside the box not drawn", {cls: "t-small", anchor: "middle", size: 22});
      s.Fv = {b: tfAt("ex02/00013TR_X.TFU", 11), r: tfAt("ex02/00138TR_X.TFU", 11), z: tfAt("ex02/00140TR_Z.TFU", 11)};
      s.Fs = {b: tfAt("ex02_fsin/00013TR_X.TFU", 11), r: tfAt("ex02_fsin/00138TR_X.TFU", 11), z: tfAt("ex02_fsin/00140TR_Z.TFU", 11)};
      s.Fs.u = [s.Fs.z[0] / cabs(s.Fs.z), s.Fs.z[1] / cabs(s.Fs.z)];          // phase of the roof-edge motion (the interior moves with it)
      s.c1 = K.card(s.root, {x: 120, y: 720, w: 780, kind: "warn", title: "A fake vibration", size: 32,
        body: "Roof slab and soil inside, held only at the edges."});
      s.c2 = K.card(s.root, {x: 120, y: 720, w: 780, kind: "check", title: "The fix: FI-EVBN", size: 32,
        body: "9 nodes back on the top face: 114 in all."});
      roofPlot(s);
      s.ev = tfLine(s.p, "ex02_evbn/00138TR_X.TFU", {color: C.good, r: 6});
      s.leg = s.p.legend([{label: "FV, 150", color: C.ref}, {label: "shortcut, 105", color: C.bad}],
        {x: 1040, y: 430, size: 28, dy: 44});
      s.legEv = K.g(s.p.g, {hidden: true});
      K.line(s.legEv, 1040, 518, 1086, 518, {stroke: C.good, "stroke-width": 6, "stroke-linecap": "round"});
      K.text(s.legEv, 1102, 527, `FI-EVBN, ${PH.setCount("evbn")}`, {cls: "t-small", color: "var(--ink2)", size: 28});
      // the largest difference from FV over the computed frequencies (2.4 %), written under the curves
      const ev = D["ex02_evbn/00138TR_X.TFU"], fv = D["ex02/00138TR_X.TFU"];
      const dmax = Math.max(...ev.amp.map((a, i) => Math.abs(a / fv.amp[i] - 1)));
      s.evL = s.p.text(12, 0.2, `within ${(100 * dmax).toFixed(1)} % of FV`, {cls: "t-label", color: "var(--good)", size: 32, hidden: true});
      s.rule = K.card(s.root, {x: 1000, y: 720, w: 800, kind: "check", title: "The rule", size: 32,
        body: "Validate any shortcut against FV first. ASCE 4-16 and SRP 3.7.2 require it."});
      s.A = 0;
    },
    tick(s, t) {
      const ph = TAU * t / TV, c = Math.cos(ph), sn = Math.sin(ph), A = s.AX * s.A;
      s.bands.forEach((b) => b.shear((d) => A * reph(s.U(d / s.PXF), c, sn)));
      s.fvS.set(rigid(s.Fv.b, s.Fv.r, s.Fv.z, A, c, sn));
      const F = s.Fs, ez = reph(F.z, c, sn), uz = reph(F.u, c, sn);
      s.fsS.set((i, k) => {
        const dx = A * reph(lerp(F.b, F.r, k / 5), c, sn);
        let dz = 0;
        if (k === 5) dz = A * [-ez, -15.6 * uz, 0, 15.6 * uz, ez][i];         // roof interior up to 15.6 at 11 Hz (lesson 5, nodes 137, 139)
        else if (i === 0 || i === 4) dz = A * (i === 0 ? -1 : 1) * ez * k / 5;
        return [dx, dz];
      });
    },
    beats: [
      {say: "[a]Why the fake peak? Watch both models shaking at 11 Hz.",
        go(k) {
          const s = k.s;
          k.show([s.head, s.p.g, s.fv.dots, s.fs.dots]); k.draw(s.fv.line, 0); k.draw(s.fs.line, 0);
        },
        a(k) { k.show(k.s.note); k.tween(k.s, {A: 1}, 1000); }},
      {say: "[v]With FV, the box follows the ground smoothly.",
        v(k) { k.show(k.s.fvL); }},
      {say: "[f]With the shortcut, the roof slab see-saws on its own, more than 15 times the ground motion.",
        f(k) { k.show(k.s.fsL); }},
      {say: "[c]Without nodes on the top face, the roof and the soil inside form a vibrating system that doesn't exist.",
        c(k) { k.show(k.s.c1); }},
      {say: "[x]The fix is cheap: put 9 nodes back on the top face. [n]That's FI-EVBN, with 114 nodes.",
        x(k) { k.hide(k.s.c1); k.show(k.s.c2, {delay: 250}); }, n(k) { k.show(k.s.legEv); }},
      {say: "[g]Now it follows FV within 2.4 percent, and the fake peak is gone.",
        g(k) { k.draw(k.s.ev.line, 1000); k.show(k.s.ev.dots, {delay: 400}); k.show(k.s.evL, {delay: 400}); k.dim([k.s.fs.line, k.s.fs.dots], true); }},
      {say: "[r]So the rule is simple, and the codes require it: validate any shortcut against FV first.", hold: 600,
        r(k) { k.show(k.s.rule); }},
    ],
  });

  // ---------------------------------------------------------------- 6. kinematic interaction: the massless box
  scenes.push({
    id: "kinematic", title: "The box averages the ground",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "The box averages the ground", {x: 120, y: 70, size: "h2"});
      const Y0 = 330, PXF = 12.5, F12 = 12.0117;               // 12.5 px per ft: the 32 ft x 16 ft box is 400 x 200 px
      s.col = column(g, {x: 100, y: Y0, w: 700, pxf: PXF, depth: 41.6, f: F12, px: 720});
      s.sec = section(g, {x0: 330, y0: Y0, sx: PXF, sz: PXF, color: C.kin, dot: C.kin, hidden: true});
      s.secL = K.text(g, 330, Y0 - 18, "massless box", {cls: "t-label", anchor: "middle", color: "var(--kin)", hidden: true});
      s.ffL = K.text(g, 720, Y0 - 18, "free field", {cls: "t-label", anchor: "middle", color: "var(--wave)", hidden: true});
      s.l1 = K.text(g, 850, Y0 + 16, "1", {cls: "t-label", size: 34, color: "var(--wave)", hidden: true});
      s.l2 = K.text(g, 850, Y0 + 16 * PXF + 12, Math.hypot(...s.col.U(16)).toFixed(2), {cls: "t-label", size: 34, color: "var(--wave)", hidden: true});
      s.slow = K.text(g, 450, 955, capOf(F12, 1000), {cls: "t-small", anchor: "middle", hidden: true});
      s.Hb = tfAt("ex02_massless/00013TR_X.TFU", 12); s.Hr = tfAt("ex02_massless/00138TR_X.TFU", 12); s.Hz = tfAt("ex02_massless/00140TR_Z.TFU", 12);
      s.kin = K.pill(s.root, "kinematic interaction", {x: 120, y: 172, color: "var(--kin)", style: {fontSize: "32px"}});
      s.p = K.plot(s.svg, {x: 1040, y: 190, w: 760, h: 400, xr: [0, 20], yr: [0, 1.2], xticks: [0, 5, 10, 15, 20], yticks: [0, 0.4, 0.8, 1.2],
        xlabel: "frequency (Hz)", ylabel: "basemat ÷ surface", ylabelOffset: 80});
      s.one = s.p.hline(1, {color: C.wave, dash: "6 6"});
      s.oneL = s.p.text(19.8, 1, "free-field surface", {cls: "t-label", anchor: "end", dy: -14, color: "var(--wave)"});
      s.ml = tfLine(s.p, "ex02_massless/00013TR_X.TFU", {color: C.kin, width: 5, r: 7});
      const ml = D["ex02_massless/00013TR_X.TFU"];
      let i10 = 0; ml.f.forEach((f, j) => { if (Math.abs(f - 10) < Math.abs(ml.f[i10] - 10)) i10 = j; });
      s.mk = K.g(s.p.g, {hidden: true});
      K.circle(s.mk, s.p.X(ml.f[i10]), s.p.Y(ml.amp[i10]), 14, {fill: "none", stroke: C.kin, "stroke-width": 4});
      K.text(s.mk, s.p.X(ml.f[i10]) - 24, s.p.Y(ml.amp[i10]) + 52, `${ml.amp[i10].toFixed(2)} at 10 Hz`, {cls: "t-label", anchor: "end", color: "var(--kin)", size: 34});
      s.fim = K.card(s.root, {x: 1000, y: 760, w: 800, kind: "check", title: "Foundation input motion", size: 32,
        body: "What your building actually receives from the ground."});
      s.A = 0; s.AP = ampPx(1000, PXF);
    },
    tick(s, t) {
      const ph = TAU * t / TV, c = Math.cos(ph), sn = Math.sin(ph), A = s.AP * s.A;
      s.col.set(c, sn, A);
      s.sec.set(rigid(s.Hb, s.Hr, s.Hz, A, c, sn));
    },
    beats: [
      {say: "[m]To see it, take away the box's mass, so only the ground can move it.",
        go(k) { k.show([k.s.head, k.s.col.prof, k.s.ffL]); },
        m(k) { k.show([k.s.sec.g, k.s.secL]); }},
      {say: "[a]At 12 Hz, the surface moves 1 and the base level 0.10. [k]The stiff box follows neither.",
        a(k) { k.tween(k.s, {A: 1}, 1000); k.show([k.s.l1, k.s.l2, k.s.slow]); }, k(k) { k.pulse(k.s.sec.g, {amp: 0.04}); }},
      {say: "It can't follow every wiggle: it averages the motion over its depth, [r]and it rocks.",
        r(k) { k.pulse(k.s.sec.roof, {amp: 0.06}); }},
      {say: "[e]Engineers call this kinematic interaction.",
        e(k) { k.show(k.s.kin); }},
      {say: "[p]Across frequencies, the basemat moves less than the surface: [n]0.39 instead of 1 at 10 Hz.",
        p(k) { k.show(k.s.p.g); k.draw(k.s.ml.line, 1000); k.show(k.s.ml.dots, {delay: 400}); }, n(k) { k.show(k.s.mk); }},
      {say: "[f]That's the foundation input motion: what your building actually receives from the ground.", hold: 600,
        f(k) { k.show(k.s.fim); }},
    ],
  });

  // ---------------------------------------------------------------- 7. what it means for your design
  scenes.push({
    id: "design", title: "What it means for design",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "What it means for your design", {x: 120, y: 70, size: "h2"});
      // left: a clamped box driven by the surface motion; right: the massless box at 10 Hz (computed); 10.625 px per ft
      const PXF = 10.625, Y0 = 290, F10 = 10.01;
      s.PXF = PXF; s.AP = ampPx(1000, PXF); s.U = freeField(F10);
      K.text(g, 510, 236, "fixed base, surface motion", {cls: "t-label", anchor: "middle", color: "var(--ref)", size: 34});
      K.text(g, 1410, 236, "foundation input motion", {cls: "t-label", anchor: "middle", color: "var(--kin)", size: 34});
      s.left = section(g, {x0: 510, y0: Y0, sx: PXF, sz: PXF});
      s.hatch = K.ground(g, 510 - 230, Y0 + 16 * PXF + 6, 460, {color: "var(--ref)"});
      s.band = K.soil(g, {x: 1410 - 280, y: Y0, w: 560, layers: [{h: 16 * PXF, kind: "sand"}, {h: 6.4 * PXF, kind: "gravel"}], labels: false});
      s.right = section(g, {x0: 1410, y0: Y0, sx: PXF, sz: PXF, color: C.kin, dot: C.kin});
      K.text(g, 510, 560, "clamped: moves 1 × the surface motion", {cls: "t-small", anchor: "middle"});
      K.text(g, 1410, 560, `massless SASSI run: ${cabs(tfAt("ex02_massless/00013TR_X.TFU", 10)).toFixed(2)} at the basemat`, {cls: "t-small", anchor: "middle"});
      K.text(g, 960, 600, capOf(F10, 1000), {cls: "t-small", anchor: "middle"});
      s.Hb = tfAt("ex02_massless/00013TR_X.TFU", 10); s.Hr = tfAt("ex02_massless/00138TR_X.TFU", 10); s.Hz = tfAt("ex02_massless/00140TR_Z.TFU", 10);
      s.c1 = K.card(s.root, {x: 120, y: 660, w: 780, kind: "warn", title: "Surface motion as input", size: 32,
        body: "✗ too strong above about 5 Hz"});
      s.c1b = cardPart(s.c1, "✗ misses the rocking");
      s.c2 = K.card(s.root, {x: 1020, y: 660, w: 780, kind: "check", title: "Foundation input motion", size: 32,
        body: "✓ translation and rocking"});
      cardPart(s.c2, "✓ from a massless SASSI run").classList.remove("sv-in", "sv-fade");
      s.A = 0;
    },
    tick(s, t) {
      const ph = TAU * t / TV, c = Math.cos(ph), sn = Math.sin(ph), A = s.AP * s.A;
      s.left.set(() => [A * c, 0]);
      s.hatch.setAttribute("transform", `translate(${A * c} 0)`);
      s.band.shear((d) => A * reph(s.U(d / s.PXF), c, sn));
      s.right.set(rigid(s.Hb, s.Hr, s.Hz, A, c, sn));
    },
    beats: [
      {say: "[a]Drive a fixed-base model with the surface motion, and above about 5 Hz, you overestimate the input.",
        go(k) { k.show(k.s.head); k.tween(k.s, {A: 1}, 1000); },
        a(k) { k.show(k.s.c1); }},
      {say: "[b]You also miss the rocking that the basement picks up from the ground.",
        b(k) { k.show(k.s.c1b); k.pulse(k.s.right.g, {amp: 0.04}); }},
      {say: "[d]The consistent input is the foundation input motion: translation and rocking, from a massless SASSI run.", hold: 600,
        d(k) { k.show(k.s.c2); }},
    ],
  });

  // ---------------------------------------------------------------- 8. recap
  scenes.push({
    id: "recap", title: "Recap",
    build(s) {
      s.head = K.heading(s.root, "Recap", {x: 120, y: 80, size: "h1"});
      s.list = K.bullets(s.root, [
        {t: "Free field + box − excavated soil", sub: "the flexible-volume method, FV"},
        {t: "Check every shortcut against FV", sub: "fewer nodes can invent a fake resonance"},
        {t: "A stiff basement averages the ground", sub: "less motion at high frequency, plus rocking"},
      ], {x: 120, y: 230, w: 1000, num: true, size: 44});
      s.q = K.card(s.root, {x: 1180, y: 200, w: 640, kind: "check", title: "Check yourself", size: 34,
        body: "A shortcut run shows a sharp peak where a stiff, light box has no mode. What do you do?"});
      s.a = K.card(s.root, {x: 1180, y: 580, w: 640, title: "Answer", size: 34,
        body: "<b>Suspect a fake resonance:</b> compare with an FV run before trusting it."});
      s.next = K.pill(s.root, "Next · Lesson 6: seismic input, spectrum-compatible motion and strain-compatible soil →", {x: 120, y: 880, color: "var(--wave)", style: {fontSize: "30px"}});
    },
    beats: [
      {say: "[a]One: SSI equals free field, plus structure, minus the soil you dug out: the flexible-volume method.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[b]Two: shortcuts with fewer nodes can invent a fake resonance, so always check them against FV.",
        b(k) { k.show(k.s.list.items[1]); }},
      {say: "[c]Three: a stiff basement averages the ground motion over its depth, and rocks: kinematic interaction.",
        c(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Check yourself: a shortcut shows a sharp peak where a stiff, light box has no mode. What do you do?", gap: 1500,
        q(k) { k.show(k.s.q); }},
      {say: "[a]Suspect a fake resonance, and compare with an FV run before you trust any number.",
        a(k) { k.show(k.s.a); }},
      {say: "[n]Next: the earthquake itself. A matching ground motion, and soil that softens as it shakes.",
        n(k) { k.show(k.s.next); }},
    ],
  });

  SV.video({
    id: "05", n: 5, part: "Fundamentals", lesson: "05-embedded",
    title: "Embedded structures and the flexible-volume method",
    subtitle: "A basement feels a ground motion that changes with depth: SASSI adds the box and subtracts the soil it replaces, shortcuts are checked against FV, and a stiff box averages the motion.",
    next: {href: "06.html", title: "Seismic input: spectrum-compatible motion and strain-compatible soil"},
    pron: [
      [/\bASCE 4-16\b/g, "A S C E 4 16"], [/\bSRP 3\.7\.2\b/g, "S R P 3 point 7 point 2"],
    ],
    scenes,
  });
})();
