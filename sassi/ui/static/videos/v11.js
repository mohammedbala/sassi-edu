/* Explainer video, lesson 11: An embedded building: the subtraction methods against FV (third edition: pace and accuracy).
 * One idea: a faster shortcut (the subtraction method, FI-FSIN) can invent a resonance of the soil trapped in
 * the basement; FI-EVBN fixes it here; always validate a shortcut against the full FV model.  Closes the course.
 * Numbers: sassi/ui/lessons/11_embedded_building.md; curves: d11.js (the lesson's runs of example 8 with the
 * FV, FI-FSIN and FI-EVBN interaction sets; python -m sassi.ui.video_data).  Units: ft, kip, s.
 * Pictures to scale: the 80 ft x 80 ft building with its 26 ft basement (two levels), walls every 20 ft, the tower over
 * the four central rooms (storeys of 16 ft), true member thicknesses and the ten equipment items where the lesson puts
 * them; the excavated soil as its 9 x 9 x 5 nodes (10 ft x 10 ft x 6.5 ft) with the interaction sets counted from it.
 * Motion: the computed transfer functions (amplitude and phase) of the output nodes at the frequency named, one
 * scale per scene; between the nodes of the excavated soil the lesson's Rayleigh shape (labelled schematic). */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["11"], PH = K.PH;
  const TAU = 2 * Math.PI;

  // ---------------------------------------------------------------- data (the lesson's runs)
  const TF = (dir, node) => D[`${dir}/${String(node).padStart(5, "0")}TR_X.TFU`];
  const RS = (dir, node) => D[`${dir}/${String(node).padStart(5, "0")}TR_X01.RS`];
  const FV = "ex08", SM = "ex08_fsin", MSM = "ex08_evbn";
  const COL = {fv: C.ref, sm: C.bad, msm: C.good};     // FV the reference; FI-FSIN the anomaly; FI-EVBN validated
  // the site of example 8 (lesson 11, step 1; ft, ft/s, kcf): 26 ft at 1,000 ft/s, 42 ft at 1,500, 39 ft at 2,100 ft/s,
  // rock at 5,000 ft/s
  const SITE8 = {g: 32.2, layers: [{h: 26, vs: 1000, w: 0.120, beta: 0.04}, {h: 42, vs: 1500, w: 0.125, beta: 0.03}, {h: 39, vs: 2100, w: 0.130, beta: 0.02}],
    hs: {vs: 5000, w: 0.145, beta: 0.01}};

  // ---------------------------------------------------------------- geometry of example 8 (ft)
  const HB = 40, DE = 26;                                // half the plan size; the embedment depth
  /** Evenly spaced values from a to b (both included), spacing at most about d. */
  const span = (a, b, d) => { const n = Math.max(1, Math.round((b - a) / d)); return Array.from({length: n + 1}, (_, i) => a + (b - a) * i / n); };

  // ---------------------------------------------------------------- complex motion from the computed transfer functions
  const near = (d, f) => d.f.reduce((b, v, i) => (Math.abs(v - f) < Math.abs(d.f[b] - f) ? i : b), 0);
  /** Complex TF value [re, im] of a TFU at the computed frequency closest to f. */
  const cAt = (d, f) => { const i = near(d, f); return [d.amp[i] * Math.cos(d.ph[i]), d.amp[i] * Math.sin(d.ph[i])]; };
  const lerpC = (a, b, t) => [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])];
  /** Output nodes of the building (X): basemat 41, basement slab 675, grade slab 757, floor 446, main roof 567, tower roof 605. */
  const LEVELS = [[-26, 41], [-13, 675], [0, 757], [16, 446], [32, 567], [64, 605]];
  /** The horizontal motion of a run at f: uz(z) of the building (linear between the output levels; the slabs carry
   *  the walls) and soil(x, z) of the excavated soil in section y = 0: the walls and basemat on its boundary, the
   *  computed top-face centre (node 365) and, between them, the Rayleigh shape sin(pi (x + 40)/80) sin(pi (z + 26)/52). */
  function field(dir, f, top) {
    const lv = LEVELS.filter(([z]) => z <= (top === undefined ? 64 : top)).map(([z, n]) => [z, cAt(TF(dir, n), f)]);
    const uz = (z) => {
      if (z <= lv[0][0]) return lv[0][1];
      for (let i = 1; i < lv.length; i++) if (z <= lv[i][0]) return lerpC(lv[i - 1][1], lv[i][1], (z - lv[i - 1][0]) / (lv[i][0] - lv[i - 1][0]));
      return lv[lv.length - 1][1];
    };
    const s0 = cAt(TF(dir, 365), f), g0 = uz(0), b = [s0[0] - g0[0], s0[1] - g0[1]];
    const soil = (x, z) => {
      const w = uz(z), sh = Math.sin(Math.PI * (x + HB) / (2 * HB)) * Math.sin(Math.PI * (z + DE) / (2 * DE));
      return [w[0] + b[0] * sh, w[1] + b[1] * sh];
    };
    return {uz, soil, f: TF(dir, 365).f[near(TF(dir, 365), f)]};
  }
  /** Real part of c e^{i wt}. */
  const re = (c, cw, sw) => c[0] * cw - c[1] * sw;

  // ---------------------------------------------------------------- helpers
  function ring(p, x, y, label, o) {
    o = o || {};
    const g = K.g(p.g, {hidden: true, in: "pop"});
    K.circle(g, p.X(x), p.Y(y), o.r || 22, {fill: "none", stroke: o.color, "stroke-width": 5});
    if (label) K.text(g, p.X(x) + (o.dx === undefined ? 32 : o.dx), p.Y(y) + (o.dy === undefined ? 10 : o.dy), label,
      {cls: "t-label", anchor: o.anchor || "start", color: o.color, size: o.size || 30});
    return g;
  }
  /** Section of the example 8 building (x across, z up, in ft; s px/ft; centre cx, grade gy), to scale: shells at their
   *  mid-surfaces with their true thickness (basemat 6.5 ft, outer basement walls 3.5, outer walls above grade 2.5,
   *  interior and tower walls 2, basement slab 2, grade slab 2.5, floors and roofs 1.5 ft), the hole of the basement,
   *  and the ten equipment items (x, level, weight in kips, slab thickness).
   *  Groups: hole, mat, outer, inner, above, tower, eq.  update(u): u(x, z) -> horizontal offset (px). */
  const EQUIP = [[-30, 0, 330, 2.5], [30, 0, 330, 2.5], [-30, -13, 220, 2], [30, -13, 220, 2], [-10, 16, 175, 1.5], [10, 16, 175, 1.5],
    [-30, 32, 55, 1.5], [30, 32, 55, 1.5], [-10, 48, 130, 1.5], [10, 64, 90, 1.5]];
  function building8(p, o) {
    const s = o.s, cx = o.cx, gy = o.gy;
    const X = (x) => cx + x * s, Z = (z) => gy - z * s;
    const g = K.g(p, {hidden: o.hidden});
    const items = [];
    const grp = {};
    for (const k of ["hole", "mat", "outer", "inner", "above", "tower", "eq"]) grp[k] = K.g(g, {hidden: o.parts !== false && k !== "hole"});
    const hole = K.path(grp.hole, "", {fill: "rgba(8,14,24,.92)", stroke: "none"});
    const brown = "#b07a4f";
    const mem = (gk, v, a, b, c, th, col) => items.push({el: K.path(grp[gk], "", {stroke: col || C.concrete, "stroke-width": th * s, fill: "none", "stroke-linecap": "butt"}), v, a, b, c});
    // v: x, z1, z2 (vertical); h: z, x1, x2 (horizontal)
    mem("mat", false, -26, -41.75, 41.75, 6.5, brown);
    for (const x of [-40, 40]) mem("outer", true, x, -26, 0, 3.5, brown);
    for (const x of [-20, 0, 20]) mem("inner", true, x, -26, 0, 2);
    mem("inner", false, -13, -38.25, 38.25, 2); mem("inner", false, 0, -41.75, 41.75, 2.5);
    for (const x of [-40, 40]) mem("above", true, x, 0, 32, 2.5);
    for (const x of [-20, 0, 20]) mem("above", true, x, 0, 32, 2);
    mem("above", false, 16, -41.25, 41.25, 1.5); mem("above", false, 32, -41.25, 41.25, 1.5);
    for (const x of [-20, 0, 20]) mem("tower", true, x, 32, 64, 2);
    mem("tower", false, 48, -21, 21, 1.5); mem("tower", false, 64, -21, 21, 1.5);
    const eqs = EQUIP.map(([x, z, w, th]) => {
      const a = 0.88 * Math.cbrt(w);                     // side (ft) of a block of the item's weight, for the picture
      return {el: K.rect(grp.eq, 0, 0, a * s, a * s, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 2, rx: 2}), x, z, a, th};
    });
    const obj = {g, grp, X, Z};
    obj.update = function (u) {
      u = u || (() => 0);
      const L = [], R = [];
      for (const z of span(-DE, 0, 3.25)) { L.push([X(-HB) + u(-HB, z), Z(z)]); R.push([X(HB) + u(HB, z), Z(z)]); }
      hole.setAttribute("d", K.d(L) + "L" + R.reverse().map((q) => q.join(",")).join("L") + "Z");
      for (const it of items) {
        const pts = [];
        if (it.v) for (const z of span(it.b, it.c, 3.25)) pts.push([X(it.a) + u(it.a, z), Z(z)]);
        else { const du = u(0, it.a); pts.push([X(it.b) + du, Z(it.a)], [X(it.c) + du, Z(it.a)]); }
        it.el.setAttribute("d", K.d(pts));
      }
      for (const e of eqs) {
        e.el.setAttribute("x", X(e.x) + u(e.x, e.z) - e.a * s / 2);
        e.el.setAttribute("y", Z(e.z) - e.th * s / 2 - e.a * s);
      }
    };
    obj.update();
    return obj;
  }
  /** The excavated soil in section y = 0: x = -40..40 every 10 ft, z = -26..0 every 6.5 ft (9 x 5 nodes).
   *  update(u): u(x, z) -> horizontal offset (px). nodes[]: {x, z, i, k, el}. */
  function excav(p, o) {
    const s = o.s, cx = o.cx, gy = o.gy;
    const X = (x) => cx + x * s, Z = (z) => gy - z * s;
    const g = K.g(p, {hidden: o.hidden});
    const fill = K.path(g, "", {fill: o.fill || "rgba(200,164,107,.34)", stroke: "none"});
    const lines = [];
    for (let i = 0; i <= 8; i++) lines.push({el: K.path(g, "", {stroke: "rgba(232,210,170,.45)", "stroke-width": 2, fill: "none"}), v: true, x: -HB + 10 * i});
    for (let k = 0; k <= 4; k++) lines.push({el: K.path(g, "", {stroke: "rgba(232,210,170,.45)", "stroke-width": 2, fill: "none"}), v: false, z: -DE + 6.5 * k});
    const ng = K.g(g, {});
    const nodes = [];
    for (let k = 0; k <= 4; k++) for (let i = 0; i <= 8; i++) {
      const x = -HB + 10 * i, z = -DE + 6.5 * k;
      nodes.push({x, z, i, k, el: K.circle(ng, X(x), Z(z), o.r || 8, {fill: o.node || "#3a4d66", stroke: "#0a111d", "stroke-width": 2.5})});
    }
    const obj = {g, nodes, X, Z, ng};
    const ZS = span(-DE, 0, 3.25), XS = span(-HB, HB, 5);
    obj.update = function (u) {
      u = u || (() => 0);
      for (const ln of lines) {
        const pts = [];
        if (ln.v) for (const z of ZS) pts.push([X(ln.x) + u(ln.x, z), Z(z)]);
        else for (const x of XS) pts.push([X(x) + u(x, ln.z), Z(ln.z)]);
        ln.el.setAttribute("d", K.d(pts));
      }
      const L = [], R = [];
      for (const z of ZS) { L.push([X(-HB) + u(-HB, z), Z(z)]); R.push([X(HB) + u(HB, z), Z(z)]); }
      const T = XS.map((x) => [X(x) + u(x, 0), Z(0)]);
      fill.setAttribute("d", K.d(L) + "L" + T.map((q) => q.join(",")).join("L") + "L" + R.reverse().map((q) => q.join(",")).join("L") + "Z");
      for (const n of nodes) n.el.setAttribute("cx", X(n.x) + u(n.x, n.z));
    };
    obj.update();
    return obj;
  }
  /** Is a node of the 9 x 9 x 5 excavation (i, j across, k up from the bottom) an interaction node of the set?
   *  (the sets of INTGEN: FV every node; FI-FSIN the lateral faces and the bottom; FI-EVBN also the top face). */
  const BOX = {nx: 9, ny: 9, nz: 5};
  const inSet = (m, i, j, k) => PH.inSet(m, {i, j, k}, BOX);
  /** Colour the nodes of a section (j = 4) by set, staggered by level (instant when seeking). */
  function paint(k, ex, m, col) {
    for (const n of ex.nodes) {
      const on = inSet(m, n.i, 4, n.k);
      const go = () => { n.el.setAttribute("fill", on ? col : "#3a4d66"); n.el.setAttribute("r", on ? 10 : 7); };
      if (k.instant || !on) go(); else k.at(70 * n.k, go);
    }
  }
  /** Hatched support line (a held face), normal (nx, ny) pointing away from the soil. */
  function hatch(pg, x1, y1, x2, y2, nx, ny, col) {
    K.line(pg, x1, y1, x2, y2, {stroke: col, "stroke-width": 5, "stroke-linecap": "round"});
    const L = Math.hypot(x2 - x1, y2 - y1), n = Math.floor(L / 22);
    for (let i = 0; i <= n; i++) {
      const t = i / n, x = x1 + (x2 - x1) * t, y = y1 + (y2 - y1) * t;
      K.line(pg, x, y, x + nx * 16 - ny * 9, y + ny * 16 + nx * 9, {stroke: col, "stroke-width": 2.5});
    }
  }
  /** The basement in section y = 0 for the method scenes: the excavated soil, the basemat and the outer walls.
   *  setField(fn(x, z) -> complex) and move(cw, sw, U) animate it. */
  function basement(p, o) {
    const r = {};
    r.b = building8(p, {s: o.s, cx: o.cx, gy: o.gy, parts: false, hidden: o.hidden});
    r.ex = excav(p, {s: o.s, cx: o.cx, gy: o.gy, r: 7, hidden: o.hidden});         // the soil nodes over the walls they share
    ["hole", "above", "tower", "eq", "inner"].forEach((k) => { r.b.grp[k].style.display = "none"; });
    r.X = (x) => o.cx + x * o.s; r.Z = (z) => o.gy - z * o.s; r.s = o.s;
    r.move = (F, cw, sw, U) => {
      r.ex.update((x, z) => U * re(F.soil(x, z), cw, sw));
      r.b.update((x, z) => U * re(F.uz(z), cw, sw));
    };
    return r;
  }
  function heldFaces(pg, B, col, top) {
    const g = K.g(pg, {hidden: true});
    const off = 2.3 * B.s;                                // just outside the outer walls and the basemat
    hatch(g, B.X(-HB) - off - 6, B.Z(0), B.X(-HB) - off - 6, B.Z(-DE) + off + 14, -1, 0, col);
    hatch(g, B.X(HB) + off + 6, B.Z(0), B.X(HB) + off + 6, B.Z(-DE) + off + 14, 1, 0, col);
    hatch(g, B.X(-HB) - off - 6, B.Z(-DE) + off + 14, B.X(HB) + off + 6, B.Z(-DE) + off + 14, 0, 1, col);
    if (top) hatch(g, B.X(-HB) - off - 6, B.Z(0) - 12, B.X(HB) + off + 6, B.Z(0) - 12, 0, -1, col);
    return g;
  }

  // ================================================================== scenes
  const scenes = [];
  const F65 = field(FV, 6.5);                                  // the building at its SSI peak (FV)
  const FF65 = PH.columnWaves(SITE8, F65.f);                    // the free field at the same frequency
  const ff = (wv, depth) => PH.columnU(wv, SITE8, depth);
  const SLOW65 = Math.round(F65.f * 2.2);                      // the 6.5 Hz motion shown with a period of 2.2 s

  // ---------------------------------------------------------------- 0. hook
  scenes.push({
    id: "hook", title: "Can a shortcut invent a resonance?",
    build(s) {
      const g = K.g(s.svg), S = 6, gy = 550;                   // 6 px/ft
      s.S = S;
      s.soil = K.soil(g, {x: 110, y: gy, w: 880, layers: [{h: 26 * S, kind: "sand"}, {h: 42 * S, kind: "gravel"}], labels: false});
      s.b = building8(g, {s: S, cx: 550, gy, parts: false});
      s.ex = excav(g, {s: S, cx: 550, gy, fill: "rgba(0,0,0,0)", r: 8, hidden: true});
      s.ex.nodes.forEach((n) => n.el.setAttribute("fill", C.ssi));
      s.inner = s.ex.nodes.filter((n) => !inSet("fsin", n.i, 4, n.k)).map((n) => n.el);
      s.note = K.text(g, 110, 990, `computed at ${F65.f.toFixed(1)} Hz (FV) · 1 × ground = 10 px · slowed down ${SLOW65} ×`, {cls: "t-small", color: "var(--muted)"});
      s.big = K.stat(s.root, {x: 1060, y: 190, w: 760, value: "30,000–50,000", label: "points tied to the soil:<br>a nuclear island, full method", color: "var(--ssi)", vsize: 96});
      s.cut = K.pill(s.root, "shortcuts keep only some of them", {x: 1110, y: 470, color: "var(--wave)", style: {fontSize: "32px"}});
      s.q = K.heading(s.root, "Can a shortcut invent a resonance?", {x: 1060, y: 610, w: 760, size: "h1"});
      s.A = 0;
    },
    tick(s, t) {
      const w = TAU * t / 2.2, cw = Math.cos(w), sw = Math.sin(w), U = 10 * s.A;
      s.b.update((x, z) => U * re(F65.uz(z), cw, sw));
      s.soil.shear((d) => U * re(ff(FF65, d / s.S), cw, sw));
    },
    beats: [
      {say: "Real nuclear buildings sit partly buried in the ground, and their SSI models get huge.",
        go(k) { k.show(k.s.b.g); k.tween(k.s, {A: 1}, 1000, {delay: 300}); }},
      {say: "[n]In the full method, a nuclear island can have 30,000 to 50,000 points tied to the soil.",
        n(k) { k.tween(k.s, {A: 0}, 600); k.show(k.s.ex.g); k.hide(k.s.note); k.show(k.s.big.el, {delay: 300}); }},
      {say: "[s]So engineers reach for shortcuts that keep only some of those points.",
        s(k) { k.dim(k.s.inner); k.show(k.s.cut); }},
      {say: "[q]Can a shortcut quietly invent a resonance that isn't there? On this building, yes.",
        q(k) { k.show(k.s.q); }},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 11",
    build(s) {
      s.t = K.titleCard(s.root, {n: 11, part: "Advanced", title: "An embedded building: the subtraction methods against FV",
        sub: "A shortcut against the full model, on a real building. The last lesson of the course."});
    },
    beats: [
      {say: "Lesson eleven, the last of the course: an embedded building, the subtraction methods against FV.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 100}); k.show(t.title, {delay: 200}); k.show(t.sub, {delay: 300}); k.show(t.rule, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. the building
  scenes.push({
    id: "building", title: "Example 8: the building",
    build(s) {
      const g = K.g(s.svg), S = 6, cx = 520, gy = 600;         // 6 px/ft
      s.S = S;
      s.head = K.heading(s.root, "Example 8: a building 26 ft deep", {x: 110, y: 70, size: "h2"});
      s.soil = K.soil(g, {x: 100, y: gy, w: 840, labels: false, layers: [{h: DE * S, kind: "sand"}], hs: {h: 230, kind: "gravel"}});
      s.soilL = K.g(g, {hidden: true});
      K.text(s.soilL, 990, gy + 70, "sand and gravel, 26 ft", {cls: "t-label", size: 28});
      K.text(s.soilL, 990, gy + 102, "Vs = 1,000 ft/s", {cls: "t-small", size: 24});
      K.text(s.soilL, 990, gy + 230, "denser gravels, then rock", {cls: "t-label", size: 28});
      K.text(s.soilL, 990, gy + 262, "Vs = 1,500 → 5,000 ft/s", {cls: "t-small", size: 24});
      s.b = building8(g, {s: S, cx, gy});
      const X = s.b.X, Z = s.b.Z;
      s.d24 = K.dim(g, X(-HB), Z(64) - 34, X(HB), Z(64) - 34, "80 ft", {hidden: true});
      s.d8 = K.dim(g, X(HB) + 46, Z(0), X(HB) + 46, Z(-DE), "26 ft: 2 levels", {hidden: true, labelColor: "var(--ink)"});
      s.d10 = K.dim(g, X(HB) + 46, Z(0), X(HB) + 46, Z(32), "2 × 16 ft", {hidden: true, labelColor: "var(--ink)"});
      s.dT = K.dim(g, X(20) + 40, Z(32), X(20) + 40, Z(64), "2 × 16 ft", {hidden: true, labelColor: "var(--ink)"});
      for (const d of [s.d8, s.d10, s.dT]) { const t = d.querySelector("text"); t.style.paintOrder = "stroke"; t.style.stroke = "#0a111d"; t.style.strokeWidth = "6px"; }
      s.bLab = K.text(g, X(-HB) - 26, Z(-13) + 10, "basement", {cls: "t-label", anchor: "end", color: "var(--ink)", hidden: true});
      s.tow = K.text(g, X(-20) - 22, Z(48) + 10, "tower", {cls: "t-label", anchor: "end", color: "var(--ink)", hidden: true});
      s.eqL = K.g(g, {hidden: true});
      K.text(s.eqL, X(-HB) - 26, Z(16) + 2, "equipment:", {cls: "t-small", anchor: "end", color: "var(--ssi)", size: 24});
      K.text(s.eqL, X(-HB) - 26, Z(16) + 30, "1,780 kips", {cls: "t-small", anchor: "end", color: "var(--ssi)", size: 24});
      // plan: 80 ft x 80 ft, walls every 20 ft (16 rooms), the tower over the four central rooms
      s.plan = K.g(g, {hidden: true});
      const ps = 2.7, px = 1420, py = 520;                     // 2.7 px/ft
      K.rect(s.plan, px, py, 80 * ps, 80 * ps, {fill: "#2b3a4f", stroke: C.concrete, "stroke-width": 5});
      K.rect(s.plan, px + 20 * ps, py + 20 * ps, 40 * ps, 40 * ps, {fill: "rgba(255,183,77,.22)", stroke: "none"});
      for (const t of [20, 40, 60]) {
        K.line(s.plan, px + t * ps, py, px + t * ps, py + 80 * ps, {stroke: C.concrete, "stroke-width": 3});
        K.line(s.plan, px, py + t * ps, px + 80 * ps, py + t * ps, {stroke: C.concrete, "stroke-width": 3});
      }
      K.text(s.plan, px + 40 * ps, py - 18, "plan: 16 rooms", {cls: "t-label", anchor: "middle", size: 26});
      K.text(s.plan, px + 40 * ps, py + 80 * ps + 38, "tower over the 4 central ones", {cls: "t-small", anchor: "middle", color: "var(--ssi)", size: 24});
      s.f65 = K.stat(s.root, {x: 1240, y: 190, w: 560, value: F65.f.toFixed(1), unit: "Hz", label: "where its floors respond most", color: "var(--ssi)", vsize: 110});
      s.note = K.g(g, {hidden: true});
      K.text(s.note, 1520, 432, `computed at ${F65.f.toFixed(1)} Hz (FV) · 1 × ground = 10 px`, {cls: "t-small", anchor: "middle", color: "var(--muted)"});
      K.text(s.note, 1520, 462, `slowed down ${SLOW65} ×`, {cls: "t-small", anchor: "middle", color: "var(--muted)"});
      s.A = 0;
    },
    tick(s, t) {
      const w = TAU * t / 2.2, cw = Math.cos(w), sw = Math.sin(w), U = 10 * s.A;
      s.b.update((x, z) => U * re(F65.uz(z), cw, sw));
      s.soil.shear((d) => U * re(ff(FF65, d / s.S), cw, sw));
    },
    beats: [
      {say: "Meet example 8: a concrete shear-wall building, [p]80 ft square, [m]with a basement 26 ft deep.",
        go(k) { k.show(k.s.head); k.show(k.s.b.grp.above, {delay: 200}); },
        p(k) { k.show(k.s.d24); }, m(k) { k.show([k.s.b.grp.mat, k.s.b.grp.outer, k.s.d8, k.s.bLab]); k.show(k.s.d10, {delay: 300}); }},
      {say: "[i]Walls every 20 ft make sixteen rooms, [t]and a tower rises from the middle.",
        i(k) { k.show([k.s.b.grp.inner, k.s.plan]); }, t(k) { k.show([k.s.b.grp.tower, k.s.tow, k.s.dT]); k.show([k.s.b.grp.eq, k.s.eqL], {delay: 300}); }},
      {say: "[s]It sits in 26 ft of sand and gravel, over denser gravels and rock.",
        s(k) { k.show(k.s.soilL); }},
      {say: "[f]Shaken on this soil, its floors respond most strongly at about 6.5 hertz.",
        f(k) { k.tween(k.s, {A: 1}, 1000); k.show(k.s.f65.el, {delay: 300}); k.show(k.s.note, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 3. the full method and the shortcut
  // the excavated soil in 3D: 9 x 9 x 5 nodes, 10 ft x 10 ft in plan, 6.5 ft vertically (80 ft x 80 ft x 26 ft)
  const P3 = {s: 6.6, cx: 520, cy: 420, az: 30 * Math.PI / 180, el: 28 * Math.PI / 180};     // 6.6 px/ft
  const proj = (x, y, z) => {
    const ca = Math.cos(P3.az), sa = Math.sin(P3.az), x1 = x * ca - y * sa, y1 = x * sa + y * ca;
    return [P3.cx + P3.s * x1, P3.cy - P3.s * (z * Math.cos(P3.el) + y1 * Math.sin(P3.el)), y1];
  };
  const BOXN = [];
  for (let k = 0; k < 5; k++) for (let j = 0; j < 9; j++) for (let i = 0; i < 9; i++) BOXN.push({i, j, k, x: -HB + 10 * i, y: -HB + 10 * j, z: -DE + 6.5 * k});
  const COUNT = (m) => BOXN.filter((n) => inSet(m, n.i, n.j, n.k)).length;
  scenes.push({
    id: "methods", title: "The full method and the shortcut",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "The full method, and a shortcut", {x: 110, y: 70, size: "h2"});
      s.box = K.g(g, {hidden: true});
      const face = (pts, op) => K.path(s.box, K.d(pts.map((q) => proj(...q))) + "Z", {fill: `rgba(200,164,107,${op})`, stroke: "rgba(232,210,170,.55)", "stroke-width": 2});
      const c = (x, y, z) => [x, y, z];
      // back faces and the bottom, then the nodes (far first), then the front faces and the top
      face([c(-HB, HB, -DE), c(HB, HB, -DE), c(HB, HB, 0), c(-HB, HB, 0)], 0.10);
      face([c(HB, -HB, -DE), c(HB, HB, -DE), c(HB, HB, 0), c(HB, -HB, 0)], 0.10);
      face([c(-HB, -HB, -DE), c(HB, -HB, -DE), c(HB, HB, -DE), c(-HB, HB, -DE)], 0.16);
      const order = BOXN.map((n) => Object.assign({p: proj(n.x, n.y, n.z)}, n)).sort((a, b) => b.p[2] - a.p[2] || a.k - b.k);
      s.nodes = order.map((n) => {
        const onFace = n.i === 0 || n.j === 0 || n.i === 8 || n.j === 8 || n.k === 0 || n.k === 4;
        n.r = onFace ? 6 : 4.5;
        n.el = K.circle(s.box, n.p[0], n.p[1], n.r, {fill: "#5a6d86", stroke: "#0a111d", "stroke-width": 1.5});
        return n;
      });
      face([c(-HB, -HB, -DE), c(-HB, HB, -DE), c(-HB, HB, 0), c(-HB, -HB, 0)], 0.06);
      face([c(-HB, -HB, -DE), c(HB, -HB, -DE), c(HB, -HB, 0), c(-HB, -HB, 0)], 0.06);
      face([c(-HB, -HB, 0), c(HB, -HB, 0), c(HB, HB, 0), c(-HB, HB, 0)], 0.05);
      const lab = (q, t, dx, dy, anchor) => { const p = proj(...q); K.text(s.box, p[0] + dx, p[1] + dy, t, {cls: "t-small", anchor: anchor || "middle", size: 24, color: "var(--ink2)"}); };
      lab([0, -HB, -DE], "80 ft", 0, 46);
      lab([-HB, 0, -DE], "80 ft", -30, 40, "end");
      lab([-HB, HB, -13], "26 ft", -18, 8, "end");
      K.text(s.box, P3.cx, 850, "the excavated soil: 9 × 9 × 5 nodes, 10 ft × 10 ft × 6.5 ft", {cls: "t-label", anchor: "middle", size: 28});
      s.minus = K.text(g, P3.cx, 205, "− the block of soil the basement replaces", {cls: "t-label", anchor: "middle", color: "var(--ssi)", size: 32, hidden: true});
      s.n405 = K.stat(s.root, {x: 1080, y: 190, w: 340, value: String(COUNT("fv")), label: "points held:<br>FV, the full method", color: "var(--ssi)", vsize: 110});
      s.n209 = K.stat(s.root, {x: 1460, y: 190, w: 360, value: String(COUNT("fsin")), label: "points held:<br>FI-FSIN, the shortcut", color: "var(--bad)", vsize: 110});
      s.half = K.pill(s.root, "about half", {x: 1555, y: 450, color: "var(--bad)", style: {fontSize: "32px"}});
      s.ref = K.pill(s.root, "FV = the reference", {x: 1120, y: 450, color: "var(--ssi)", style: {fontSize: "32px"}});
      s.faces = K.card(s.root, {x: 1080, y: 600, w: 740, title: "FI-FSIN keeps", size: 32, body: "the four sides and the bottom: 81 + 4 × 32 nodes"});
    },
    beats: [
      {say: "[f]It starts from the ground with no hole, [s]then subtracts the block of soil the basement replaces.",
        go(k) { k.show(k.s.head); }, f(k) { k.show(k.s.box); }, s(k) { k.show(k.s.minus); k.pulse(k.s.minus); }},
      {say: "[n]In the full method, the ground holds that block at every one of its points: 405 here.",
        n(k) {
          for (const n of k.s.nodes) {
            const go = () => { n.el.setAttribute("fill", C.ssi); n.el.setAttribute("r", n.r + 1); };
            if (k.instant) go(); else k.at(110 * n.k, go);
          }
          k.show(k.s.n405.el, {delay: 400});
        }},
      {say: "[v]That's the flexible-volume method, FV: the reference.", v(k) { k.show(k.s.ref); }},
      {say: "[s]The subtraction method, FI-FSIN, holds the block only on its sides and bottom.",
        s(k) {
          for (const n of k.s.nodes) {
            const on = inSet("fsin", n.i, n.j, n.k);
            const go = () => { n.el.setAttribute("fill", on ? C.bad : "#3a4d66"); n.el.setAttribute("r", on ? n.r + 1 : n.r - 1); };
            if (k.instant) go(); else k.at(on ? 110 * n.k : 0, go);
          }
          k.show(k.s.n209.el, {delay: 300}); k.show([k.s.half, k.s.faces], {delay: 400});
        }},
    ],
  });

  // ---------------------------------------------------------------- 4. a resonance that is not there
  scenes.push({
    id: "fsin", title: "A resonance that isn't there",
    build(s) {
      s.head = K.heading(s.root, "Same building, two methods", {x: 110, y: 70, size: "h2"});
      s.p = K.plot(s.svg, {x: 230, y: 190, w: 1040, h: 470, xr: [0, 20], yr: [0, 3], xticks: [0, 4, 8, 10, 12, 16, 20], yticks: [0, 1, 2, 3],
        xlabel: "frequency (Hz)", ylabel: "tower roof amplification", ylabelOffset: 60});
      const dF = TF(FV, 605), dS = TF(SM, 605);
      s.ok = K.rect(s.p.g, s.p.X(0), s.p.box.y, s.p.X(10) - s.p.X(0), s.p.box.h, {fill: "rgba(99,230,164,.08)", hidden: true});
      s.okL = s.p.text(5, 0.3, "agree within 2.2 %", {cls: "t-label", anchor: "middle", color: "var(--good)", hidden: true});
      s.band = K.rect(s.p.g, s.p.X(15.5), s.p.box.y, s.p.X(16) - s.p.X(15.5), s.p.box.h, {fill: "rgba(255,107,107,.22)", hidden: true});
      s.lF = s.p.line(dF.f, dF.amp, {color: COL.fv, width: 12, draw: true, opacity: 0.75});
      s.lS = s.p.line(dS.f, dS.amp, {color: COL.sm, width: 4, draw: true});
      // the shortcut's largest computed value of the 15-17 Hz band (its spurious resonance), against FV there
      const iS = dS.f.reduce((b, f, i) => (f >= 15 && f <= 17 && !(dS.f[b] >= 15 && dS.amp[b] >= dS.amp[i]) ? i : b), 0);
      s.k = ring(s.p, dS.f[iS], dS.amp[iS], `${dS.amp[iS].toFixed(3)} against ${dF.amp[iS].toFixed(3)}`, {color: "var(--bad)", dx: -32, dy: -26, anchor: "end", r: 20});
      s.leg = s.p.legend([{label: "FV, the full method", color: COL.fv}, {label: "FI-FSIN, the shortcut", color: COL.sm}], {x: 860, y: 236, hidden: true, size: 28, dy: 46});
      s.s1 = K.stat(s.root, {x: 1340, y: 200, w: 480, value: F65.f.toFixed(1), unit: "Hz", label: "main peak: both agree", color: "var(--good)", vsize: 100});
      s.s2 = K.stat(s.root, {x: 1340, y: 430, w: 480, value: "15.5–16", unit: "Hz", label: "an extra peak, shortcut only", color: "var(--bad)", vsize: 84});
      s.none = K.card(s.root, {x: 230, y: 820, w: 1590, kind: "warn", title: "No error", size: 32,
        body: "Just a reminder to validate. A check of the main peak alone would pass it."});
    },
    beats: [
      {say: "So we run the building both ways, [r]and compare how much the tower roof amplifies the shaking.",
        go(k) { k.show(k.s.head); k.show(k.s.p.g, {delay: 200}); },
        r(k) { k.draw(k.s.lF, 1000); k.draw(k.s.lS, 1000, {to: 0.5, delay: 300}); k.show(k.s.leg, {delay: 400}); }},
      {say: "[a]Up to 10 hertz, main peak included, the two agree within 2.2 percent.",
        a(k) { k.show([k.s.ok, k.s.okL]); k.show(k.s.s1.el, {delay: 300}); }},
      {say: "[b]Then, between 15.5 and 16 hertz, the shortcut shows a resonance that FV doesn't have.",
        b(k) { k.draw(k.s.lS, 1000, {from: 0.5}); k.show(k.s.band, {delay: 300}); k.show([k.s.k, k.s.s2.el], {delay: 400}); }},
      {say: "[w]No error, just a reminder to validate. A check of the main peak alone would pass it.",
        w(k) { k.show(k.s.none); }},
    ],
  });

  // ---------------------------------------------------------------- 5. why: the soil trapped in the basement
  // the computed motion at 16 Hz, FI-FSIN above and FV below, one scale: 1 x the ground = 4 px
  const F16 = {sm: field(SM, 16, 0), fv: field(FV, 16, 0), msm: field(MSM, 16, 0)};
  const U16 = 4, PER16 = 1.6, SLOW16 = Math.round(F16.sm.f * PER16);
  scenes.push({
    id: "why", title: "Why: the soil trapped in the basement",
    build(s) {
      const g = K.g(s.svg), S = 5.4, cx = 470, gA = 290, gB = 630, xl = cx + HB * S + 44;     // 5.4 px/ft
      s.head = K.heading(s.root, "Why: the soil trapped inside", {x: 110, y: 70, size: "h2"});
      s.A = basement(g, {s: S, cx, gy: gA, hidden: true});
      s.B = basement(g, {s: S, cx, gy: gB, hidden: true});
      s.A.ex.nodes.forEach((n) => { if (inSet("fsin", n.i, 4, n.k)) { n.el.setAttribute("fill", C.bad); n.el.setAttribute("r", 9); } });
      s.B.ex.nodes.forEach((n) => { n.el.setAttribute("fill", C.ssi); n.el.setAttribute("r", 9); });
      s.tA = K.text(g, cx - HB * S, gA - 26, "FI-FSIN, the shortcut", {cls: "t-label", size: 30, color: "var(--bad)", hidden: true});
      s.tB = K.g(g, {hidden: true});
      K.text(s.tB, cx - HB * S, gB - 26, "FV, the full method", {cls: "t-label", size: 30, color: "var(--ink)"});
      K.text(s.tB, xl, gB + 80, "every node held", {cls: "t-label", size: 28, color: "var(--ssi)"});
      s.held = heldFaces(g, s.A, C.ssi, false);
      s.heldL = K.text(g, xl, gA + 130, "held: sides and bottom", {cls: "t-label", size: 28, color: "var(--ssi)", hidden: true});
      s.free = K.text(g, xl, gA + 30, "free: inside and top", {cls: "t-label", size: 28, color: "var(--bad)", hidden: true});
      s.own = K.text(g, cx, gA + DE * S + 78, "its own frequency ≈ 15.6 Hz (eigenvalue analysis)", {cls: "t-small", anchor: "middle", size: 24, color: "var(--ink2)", hidden: true});
      s.note = K.g(g, {hidden: true});
      K.text(s.note, 110, 862, `computed at ${F16.sm.f.toFixed(1)} Hz: basemat, slabs, soil centre · 1 × ground = ${U16} px · slowed down ${SLOW16} ×`,
        {cls: "t-small", size: 22, color: "var(--muted)"});
      K.text(s.note, 110, 890, "the soil between those nodes: schematic shape", {cls: "t-small", size: 22, color: "var(--muted)"});
      s.p = K.plot(s.svg, {x: 1240, y: 200, w: 570, h: 380, xr: [0, 20], yr: [0, 20], xticks: [0, 4, 8, 12, 16, 20], yticks: [0, 5, 10, 15, 20],
        xlabel: "frequency (Hz)", ylabel: "soil inside: amplification", ylabelOffset: 60});
      const sF = TF(FV, 365), sS = TF(SM, 365), i16 = near(sS, 16);
      s.lF = s.p.line(sF.f, sF.amp, {color: COL.fv, width: 6, draw: true});
      s.lS = s.p.line(sS.f, sS.amp, {color: COL.sm, width: 5, draw: true});
      s.k = ring(s.p, sS.f[i16], sS.amp[i16], `${sS.amp[i16].toFixed(1)} ×`, {color: "var(--bad)", dx: -32, dy: 12, anchor: "end", r: 20, size: 38});
      s.kF = ring(s.p, sF.f[i16], sF.amp[i16], "", {color: "var(--ink)", r: 16});
      K.line(s.kF, s.p.X(13.6), s.p.Y(5.4), s.p.X(sF.f[i16]) - 8, s.p.Y(sF.amp[i16]) - 14, {stroke: C.ink, "stroke-width": 2, "stroke-dasharray": "6 5"});
      K.text(s.kF, s.p.X(13.6), s.p.Y(6.2), `FV: ${sF.amp[i16].toFixed(1)} ×`, {cls: "t-label", anchor: "end", size: 30, color: "var(--ink)"});
      s.leg = s.p.legend([{label: "FV", color: COL.fv}, {label: "FI-FSIN", color: COL.sm}], {x: 1280, y: 236, hidden: true, size: 28, dy: 46});
      const b41 = (d) => cAt(TF(d, 41), 16);
      s.drag = K.text(g, cx, gA + DE * S + 110, `basemat ${Math.hypot(...b41(SM)).toFixed(2)} × ground, against ${Math.hypot(...b41(FV)).toFixed(2)} with FV`,
        {cls: "t-small", anchor: "middle", size: 24, color: "var(--ink2)", hidden: true});
      s.name = K.heading(s.root, "A spurious resonance: made by the method", {x: 110, y: 910, w: 1700, size: "h2", color: "var(--bad)"});
      s.M = 0; s.MB = 0;
    },
    tick(s, t) {
      const w = TAU * t / PER16, cw = Math.cos(w), sw = Math.sin(w);
      s.A.move(F16.sm, cw, sw, U16 * s.M);
      s.B.move(F16.fv, cw, sw, U16 * s.MB);
    },
    beats: [
      {say: "Where does it come from? [l]Look at the soil inside the basement.",
        go(k) { k.show(k.s.head); }, l(k) { k.show([k.s.A.ex.g, k.s.A.b.g, k.s.tA]); }},
      {say: "[s]With the shortcut, that soil is held only on its sides and bottom. [t]Its inside and top are free.",
        s(k) { k.show([k.s.held, k.s.heldL]); }, t(k) { k.show(k.s.free); }},
      {say: "[j]Like jelly in a bowl, it can wobble on its own, at its own natural frequency.",
        j(k) { k.show(k.s.own); k.show(k.s.note, {delay: 200}); k.tween(k.s, {M: 1}, 1000); }},
      {say: "[m]Near 16 hertz, it sways about 18 times as much as the ground, [d]and drags the building along.",
        m(k) { k.show(k.s.p.g); k.draw(k.s.lF, 1000); k.draw(k.s.lS, 1000, {delay: 300}); k.show(k.s.leg, {delay: 300}); k.show(k.s.k, {delay: 400}); },
        d(k) { k.show(k.s.drag); k.pulse(k.s.A.b.g, {amp: 0.03}); }},
      {say: "[f]With the full method, the same soil moves only about 1.8 times the ground.",
        f(k) { k.show([k.s.B.ex.g, k.s.B.b.g, k.s.tB]); k.tween(k.s, {MB: 1}, 600); k.show(k.s.kF, {delay: 200}); k.pulse(k.s.lF); }},
      {say: "[n]Engineers call this a spurious resonance: it comes from the method, not from the building.",
        n(k) { k.show(k.s.name); }},
    ],
  });

  // ---------------------------------------------------------------- 6. what it does to the floor spectra
  scenes.push({
    id: "isrs", title: "What it does to your floor spectra",
    build(s) {
      s.head = K.heading(s.root, "Does it matter for design?", {x: 110, y: 70, size: "h2"});
      s.p = K.plot(s.svg, {x: 230, y: 190, w: 1000, h: 470, xr: [0, 25], yr: [0, 1.2], xticks: [0, 5, 10, 15, 20, 25],
        yticks: [0, 0.4, 0.8, 1.2], yfmt: (v) => v.toFixed(1), xlabel: "equipment frequency (Hz)", ylabel: "slab at ground level, 5 % (g)", ylabelOffset: 70});
      const rF = RS(FV, 757), rS = RS(SM, 757);
      // the shortcut's largest excess over FV near the spurious resonance (14.5-17 Hz)
      const i5 = rF.f.reduce((b, f, i) => (f >= 14.5 && f <= 17 && !(rF.f[b] >= 14.5 && rS.sa[b] / rF.sa[b] >= rS.sa[i] / rF.sa[i]) ? i : b), 0);
      const f5 = rF.f[i5], up = 100 * (rS.sa[i5] / rF.sa[i5] - 1);
      s.v155 = s.p.vline(f5, {color: C.bad, hidden: true});
      // the area between the two spectra around the spurious resonance and the dip below it
      const band = (d) => d.f.map((f, i) => [f, d.sa[i]]).filter(([f]) => f >= 13.5 && f <= 18);
      const top = band(rS), bot = band(rF).reverse();
      s.ex = K.path(s.p.g, K.d(s.p.pts(top.map((q) => q[0]), top.map((q) => q[1]))) + "L" +
        s.p.pts(bot.map((q) => q[0]), bot.map((q) => q[1])).map((q) => q.join(",")).join("L") + "Z", {fill: "rgba(255,107,107,.55)", hidden: true});
      s.lF = s.p.line(rF.f, rF.sa, {color: COL.fv, width: 12, draw: true, opacity: 0.75});
      s.lS = s.p.line(rS.f, rS.sa, {color: COL.sm, width: 4, draw: true});
      s.leg = s.p.legend([{label: "FV", color: COL.fv}, {label: "FI-FSIN", color: COL.sm}], {x: 960, y: 236, hidden: true, size: 28, dy: 46});
      s.k = s.p.text(f5, 1.2, `${f5.toFixed(1)} Hz`, {cls: "t-label", anchor: "middle", color: "var(--bad)", dy: -14, hidden: true});
      s.pt = K.g(s.p.g, {hidden: true});
      K.arrow(s.pt, s.p.X(f5) + 70, s.p.Y(rS.sa[i5]) - 70, s.p.X(f5) + 10, s.p.Y(rS.sa[i5]) - 10, {color: C.bad, width: 4, head: 14});
      K.text(s.pt, s.p.X(f5) + 76, s.p.Y(rS.sa[i5]) - 78, `${rS.sa[i5].toFixed(2)} g against ${rF.sa[i5].toFixed(2)} g`, {cls: "t-label", size: 26, color: "var(--bad)"});
      s.s1 = K.stat(s.root, {x: 1320, y: 190, w: 500, value: `+${up.toFixed(0)} %`, label: `slab at ground level, ${f5.toFixed(1)} Hz`, color: "var(--bad)", vsize: 96});
      s.s2 = K.stat(s.root, {x: 1320, y: 400, w: 500, value: "+16 / −10 %", label: "worst of the 11 spectra", color: "var(--bad)", vsize: 80});
      s.s3 = K.stat(s.root, {x: 1320, y: 610, w: 500, value: "≈ 1.5 %", label: "basement wall forces", color: "var(--good)", vsize: 96});
      s.q = K.pill(s.root, "designed for a demand that doesn't exist", {x: 260, y: 760, color: "var(--bad)", style: {fontSize: "30px"}});
      s.tf = K.card(s.root, {x: 230, y: 840, w: 1590, kind: "check", title: "So", size: 32,
        body: "Forces can't validate a shortcut. Compare the transfer functions themselves."});
    },
    beats: [
      {say: "Does it matter for design? [p]Here is the floor spectrum of the slab at ground level.",
        go(k) { k.show(k.s.head); }, p(k) { k.show(k.s.p.g); k.draw(k.s.lF, 1000); k.draw(k.s.lS, 1000, {delay: 300}); k.show(k.s.leg, {delay: 400}); }},
      {say: "[g]Near 15 hertz, the shortcut's spectrum is 12 percent too high. [t]Elsewhere, up to 16 percent, and 10 percent too low just below.",
        g(k) { k.show([k.s.v155, k.s.k, k.s.ex, k.s.pt]); k.show(k.s.s1.el, {delay: 300}); }, t(k) { k.show(k.s.s2.el); }},
      {say: "[q]Equipment tuned near 15.5 hertz would be designed for a demand that doesn't exist.",
        q(k) { k.show(k.s.q); }},
      {say: "[w]Yet the basement wall forces barely change: about 1.5 percent. They come from the 6.5 hertz sway.",
        w(k) { k.show(k.s.s3.el); }},
      {say: "[v]So forces can't validate a shortcut. Compare the transfer functions themselves.",
        v(k) { k.show(k.s.tf); }},
    ],
  });

  // ---------------------------------------------------------------- 7. the fix, and the rule
  scenes.push({
    id: "evbn", title: "The fix, and the rule",
    build(s) {
      const g = K.g(s.svg), S = 6, cx = 500, gy = 290;         // 6 px/ft
      s.head = K.heading(s.root, "The fix: put a lid on it", {x: 110, y: 70, size: "h2"});
      s.B = basement(g, {s: S, cx, gy});
      s.B.ex.nodes.forEach((n) => { if (inSet("fsin", n.i, 4, n.k)) { n.el.setAttribute("fill", C.bad); n.el.setAttribute("r", 9); } });
      s.held = heldFaces(g, s.B, C.ssi, false); s.held.classList.remove("sv-in", "sv-fade");
      s.lid = K.g(g, {hidden: true});
      hatch(s.lid, s.B.X(-HB) - 18, s.B.Z(0) - 14, s.B.X(HB) + 18, s.B.Z(0) - 14, 0, -1, C.good);
      s.lidL = K.text(g, cx, s.B.Z(0) - 52, `top held too: ${COUNT("evbn") - COUNT("fsin")} more points`, {cls: "t-label", anchor: "middle", color: "var(--good)", size: 32, hidden: true});
      s.note = K.text(g, 130, s.B.Z(-DE) + 76, `computed at ${F16.sm.f.toFixed(1)} Hz, the same scale as before: FI-FSIN, then FI-EVBN`,
        {cls: "t-small", size: 22, color: "var(--muted)"});
      s.n258 = K.stat(s.root, {x: 220, y: 560, w: 500, value: String(COUNT("evbn")), label: "points: FI-EVBN, the modified method", color: "var(--good)", vsize: 96});
      s.p = K.plot(s.svg, {x: 1100, y: 190, w: 700, h: 360, xr: [0, 20], yr: [0, 3], xticks: [0, 4, 8, 12, 16, 20], yticks: [0, 1, 2, 3],
        xlabel: "frequency (Hz)", ylabel: "tower roof amplification", ylabelOffset: 58});
      const dF = TF(FV, 605), dS = TF(SM, 605), dE = TF(MSM, 605);
      s.lS = s.p.line(dS.f, dS.amp, {color: COL.sm, width: 3, opacity: 0.35});
      s.lF = s.p.line(dF.f, dF.amp, {color: COL.fv, width: 11, opacity: 0.8});
      s.lE = s.p.line(dE.f, dE.amp, {color: COL.msm, width: 4, draw: true});
      s.leg = s.p.legend([{label: "FV", color: COL.fv}, {label: "FI-EVBN", color: COL.msm}, {label: "FI-FSIN", color: COL.sm}], {x: 1560, y: 226, size: 26, dy: 40});
      s.cut = K.pill(s.root, "its own frequency: 22.7 Hz, above the 20 Hz cut-off", {x: 130, y: 715, color: "var(--good)", style: {fontSize: "28px"}});
      s.s37 = K.stat(s.root, {x: 1200, y: 655, w: 520, value: "5.1 %", label: "worst difference from FV", color: "var(--good)", vsize: 84});
      s.warn = K.card(s.root, {x: 110, y: 810, w: 820, kind: "warn", title: "Not a cure", size: 30, body: "Softer soil or a deeper basement brings the wobble back."});
      s.rule = K.card(s.root, {x: 980, y: 810, w: 840, kind: "check", title: "The rule", size: 30, body: "Validate every shortcut against FV, on the softest soil."});
      s.T = 0;
    },
    tick(s, t) {
      const w = TAU * t / PER16, cw = Math.cos(w), sw = Math.sin(w), m = s.T;
      const F = {uz: (z) => lerpC(F16.sm.uz(z), F16.msm.uz(z), m), soil: (x, z) => lerpC(F16.sm.soil(x, z), F16.msm.soil(x, z), m)};
      s.B.move(F, cw, sw, U16);
    },
    beats: [
      {say: "The fix is surprisingly cheap: [l]hold the top face of the soil block as well.",
        go(k) { k.show(k.s.head); }, l(k) { k.show(k.s.lid); paint(k, k.s.B.ex, "evbn", C.good); k.tween(k.s, {T: 1}, 1000); }},
      {say: "[e]That's 49 more points, 258 in all: the modified subtraction method, FI-EVBN.",
        e(k) { k.show(k.s.lidL); k.show(k.s.n258.el, {delay: 300}); }},
      {say: "[m]With its lid on, the jelly is stiffer: its wobble moves above the 20 hertz cut-off.",
        m(k) { k.show(k.s.cut); }},
      {say: "[r]Now the shortcut follows FV within 5.1 percent, at every frequency and every output.",
        r(k) { k.show(k.s.p.g); k.draw(k.s.lE, 1000, {delay: 300}); k.show(k.s.s37.el, {delay: 400}); }},
      {say: "[w]But it's not a cure: softer soil or a deeper basement brings the wobble back into range.",
        w(k) { k.show(k.s.warn); k.show(k.s.rule, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 8. recap, and the end of the course
  scenes.push({
    id: "recap", title: "Recap, and the end of the course",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Recap", {x: 120, y: 80, size: "h2"});
      s.list = K.bullets(s.root, [
        {t: "A shortcut can invent a resonance", sub: "the soil trapped in the basement resonates"},
        {t: "Floor spectra up to 16 % too high", sub: "with no error, and wall forces unchanged"},
        {t: "Hold the top too, and validate", sub: "against FV, on the softest soil"},
      ], {x: 120, y: 200, w: 980, num: true, size: 40});
      s.list.items.forEach((li) => { li.style.marginBottom = "40px"; });
      s.q = K.card(s.root, {x: 1150, y: 190, w: 680, kind: "check", title: "Check yourself", size: 32,
        body: "The shortcut's wall forces are within about 1&nbsp;% of FV. Can you use it for design?"});
      s.a = K.card(s.root, {x: 1150, y: 470, w: 680, title: "Answer", size: 32,
        body: "<b>No.</b> Near 16&nbsp;Hz its transfer functions are 2.6&nbsp;to&nbsp;5&nbsp;times off."});
      // the course: eleven lessons, lit one after the other
      s.path = K.g(g, {hidden: true});
      const x0 = 200, dx = 152, y = 790;
      K.line(s.path, x0, y, x0 + 10 * dx, y, {stroke: "#2a3b52", "stroke-width": 6, "stroke-linecap": "round"});
      s.lit = K.line(s.path, x0, y, x0, y, {stroke: C.good, "stroke-width": 6, "stroke-linecap": "round"});
      s.dots = [];
      for (let i = 0; i < 11; i++) {
        const c = K.circle(s.path, x0 + i * dx, y, 30, {fill: "#142033", stroke: "#2a3b52", "stroke-width": 4});
        const t = K.text(s.path, x0 + i * dx, y + 10, String(i + 1), {cls: "t-label", anchor: "middle", size: 28, color: "var(--muted)"});
        s.dots.push({c, t});
      }
      s.from = K.text(s.path, x0, y - 54, "rigid ground", {cls: "t-small", anchor: "middle", color: "var(--ink2)", size: 26});
      s.to = K.text(s.path, x0 + 10 * dx, y - 54, "embedded, checked", {cls: "t-small", anchor: "middle", color: "var(--ink2)", size: 26});
      s.motto = K.heading(s.root, "Compute it. Check it. Know its limits.", {x: 120, y: 880, w: 1700, size: "h2", color: "var(--good)"});
      s.L = 0;
    },
    tick(s) {
      const x0 = 200, dx = 152;
      s.lit.setAttribute("x2", x0 + Math.min(10, s.L) * dx);
      s.dots.forEach((d, i) => {
        const on = s.L >= i - 0.01;
        d.c.setAttribute("fill", on ? "rgba(99,230,164,.25)" : "#142033");
        d.c.setAttribute("stroke", on ? "#63e6a4" : "#2a3b52");
        d.t.style.fill = on ? "var(--ink)" : "var(--muted)";
      });
    },
    beats: [
      {say: "[a]A shortcut that holds the basement soil at fewer points can invent a resonance.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[b]Here it pushed floor spectra up to 16 percent too high, with no error.",
        b(k) { k.show(k.s.list.items[1]); }},
      {say: "[c]Hold the top face too, and always validate against FV, on the softest soil.",
        c(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Check yourself. The shortcut's wall forces are within about 1 percent of FV. Can you use it for design?",
        q(k) { k.show(k.s.q); }, gap: 1500},
      {say: "[a]No. Near 16 hertz its transfer functions are 2.6 to 5 times off, and its floor spectra too high.",
        a(k) { k.show(k.s.a); }},
      {say: "[f]Compute it, check it, and know its limits. Thanks for watching, and good luck with your own SSI models.",
        go(k) { k.show(k.s.path); k.tween(k.s, {L: 10}, 1000, {ease: "linear", delay: 200}); }, f(k) { k.show(k.s.motto, {delay: 300}); }},
    ],
  });

  SV.video({
    id: "11", n: 11, part: "Advanced", lesson: "11-embedded-building",
    title: "An embedded building: the subtraction methods against FV",
    subtitle: "A faster shortcut can invent a resonance of the soil trapped in the basement: see it, fix it, and validate against FV.",
    scenes,
  });
})();
