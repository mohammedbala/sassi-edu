/* Explainer video, lesson 9 (third edition: pace and accuracy): Nonlinear soil and cracked concrete by iteration.
 * One big idea: SASSI is a linear method, but strong shaking softens soil and cracks concrete.  The trick is
 * the equivalent-linear iteration (guess the stiffness and damping, run, read the strains, update, repeat
 * until nothing changes), for the soil next to a structure (NLSSIITER) and for cracking walls (Option NON);
 * the result shifts the frequencies and changes the spectra.
 *
 * Accuracy.  Example 6 is drawn to scale in section: the 0.6 m wall (a shell in the plane x = 0, 2 m high,
 * its top at grade) with the 150 t deck lumped on its top, the 4 m x 2 m backfill of 2 x 2 elements of
 * 2 m x 1 m, in 10 m of sand (cut with a fade).  Example 7 is the 12 m box of two 4 m storeys with 0.3 m
 * walls on its 1.5 m mat over 10 m of soil (cut).  The buildings sway with the computed steady state of the
 * lesson (HARMFRAME at the peaks of the roof transfer function): uncracked 8.0 Hz with storey drifts 2.15
 * and 1.86, cracked 7.0 Hz with 2.26 and 1.31 per unit control motion (storeys taken in phase); ground
 * +-1 mm, displacements x 250, slowed down 16 times, relative to the mat.  The soil loop and the soil
 * element are the hyperbolic Masing model of the lesson figure (SV.K.PH: G_max of the 250 m/s sand,
 * reference strain of the library Sand curve); the iteration path is the mean G/G_max of the 8 backfill
 * elements per pass (NLSOIL_CONVERGENCE.TXT) on the library Sand curve, from the start GFAC = 0.4 x the
 * free-field G.  Numbers: sassi/ui/lessons/09_nonlinear.md; curves: d09.js (python -m sassi.ui.video_data). */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["09"], PH = K.PH;
  const TAU = 2 * Math.PI;

  // ---------------------------------------------------------------- helpers
  /** Spectrum {f, sa} -> peak {f, sa}. */
  function rsPeak(r) { let i = 0; r.sa.forEach((v, j) => { if (v > r.sa[i]) i = j; }); return {f: r.f[i], sa: r.sa[i]}; }
  /** Zig-zag crack in local building coordinates (lx px from the centre line, z px above the foundation). */
  function crack(x0, z0, n, dx, dz) {
    const pts = [];
    for (let i = 0; i <= n; i++) pts.push([x0 + i * dx + (i % 2 ? 7 : -7), z0 + i * dz]);
    return pts;
  }
  /** Outline of a building region (lx1..lx2, z1..z2) that follows the building's motion. */
  function quad(b, lx1, lx2, z1, z2) {
    const P = [b.T(lx1, z1), b.T(lx2, z1), b.T(lx2, z2), b.T(lx1, z2)];
    return "M" + P.map((q) => q.join(",")).join("L") + "Z";
  }
  /** A diamond bullet + text line (HTML). */
  const item = (root, x, y, w, color, html, size) => K.h("div", {x, y, w, size: size || 44, in: "left", color: "var(--ink)",
    html: `<span style="display:inline-block;width:22px;height:22px;border-radius:4px;background:${color};margin-right:30px;transform:rotate(45deg)"></span>${html}`}, root);
  /** Linear interpolation (xs ascending). */
  function lerp(xs, ys, x) {
    if (x <= xs[0]) return ys[0];
    for (let i = 1; i < xs.length; i++) if (x <= xs[i]) return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
    return ys[ys.length - 1];
  }

  // ---------------------------------------------------------------- example 7: the computed sway
  const HARM = {el: {f: 8.0, d: [2.15, 1.86]}, cr: {f: 7.0, d: [2.26, 1.31]}};   // lesson 9, part 4
  const UG = 0.001, EXAG = 250, SLOW = 16;
  const HLAB = "example 7, computed: ground ±1 mm, displacements × 250, slowed down 16 ×";
  /** The box of example 7 to scale (S px per metre) on the soil surface ySoil: 12 m, two 4 m storeys, 1.5 m mat. */
  function box(g, x, ySoil, S, o) {
    const b = K.building(g, Object.assign({x, y: ySoil - 1.5 * S, w: 12 * S, storeys: 2, storeyH: 4 * S, windows: false, matH: 1.5 * S}, o));
    return {b, S, x, sh: 4 * S, w: 12 * S, mh: 1.5 * S};
  }
  /** Sway the box: phase ph, crack state cr (0 uncracked ... 1 cracked): storey drifts of the lesson. */
  function sway(bx, ph, cr) {
    const d1 = HARM.el.d[0] + (HARM.cr.d[0] - HARM.el.d[0]) * cr, d2 = HARM.el.d[1] + (HARM.cr.d[1] - HARM.el.d[1]) * cr;
    const k = UG * EXAG * bx.S * Math.cos(ph), sh = bx.sh;
    bx.b.set({drift: (z) => (z <= 0 ? 0 : z <= sh ? k * d1 * z / sh : k * (d1 + d2 * Math.min(1, (z - sh) / sh)))});
    // the mat is 12 m wide like the box (the kit draws an overhang)
    const m = [bx.b.T(-bx.w / 2, 0), bx.b.T(bx.w / 2, 0), bx.b.T(bx.w / 2, -bx.mh), bx.b.T(-bx.w / 2, -bx.mh)];
    bx.b.mat.setAttribute("d", "M" + m.map((q) => q.join(",")).join("L") + "Z");
  }
  /** Advance the sway phase at the frequency of the crack state. */
  const advance = (s, dt) => { if (s.run) s.ph += TAU * (HARM.el.f + (HARM.cr.f - HARM.el.f) * s.Cr) * dt / SLOW; };

  // ---------------------------------------------------------------- example 6 in section (wall, deck, backfill)
  /** The wall of example 6 in section: x = 0 (the wall plane) at sx, ground at sy, m px per metre. */
  function section(g, sx, sy, m, o) {
    const r = {};
    r.soil = K.soil(g, {x: sx - 2.2 * m, y: sy, w: 8 * m, layers: [], hs: {h: o.depth * m, kind: "sand"}, labels: false});
    r.cells = [];
    for (let j = 0; j < 2; j++) for (let i = 0; i < 2; i++) {
      r.cells.push(K.rect(g, sx + i * 2 * m, sy + j * m, 2 * m, m, {fill: "rgba(232,214,170,.55)", stroke: "rgba(255,240,215,.75)", "stroke-width": 2}));
    }
    r.wall = K.rect(g, sx - 0.3 * m, sy, 0.6 * m, 2 * m, {fill: C.concreteFill, stroke: C.concrete, "stroke-width": 3});
    r.deck = K.rect(g, sx - 0.7 * m, sy - 0.45 * m, 1.4 * m, 0.45 * m, {fill: "#3b4a5e", stroke: C.concrete, "stroke-width": 3, rx: 3});
    if (o.labels) {
      r.labs = K.g(g, {hidden: true});
      K.text(r.labs, sx, sy - 0.62 * m, "deck: 150 t", {cls: "t-label", anchor: "middle", size: 28});
      K.text(r.labs, sx + 2 * m, sy + 2.55 * m, "loose backfill", {cls: "t-label", anchor: "middle", color: "var(--ink)", size: 30});
      K.text(r.labs, sx - 0.5 * m, sy + 1.1 * m, "wall", {cls: "t-label", anchor: "end", size: 28});
      K.text(r.labs, sx + 5.7 * m, sy + 3.3 * m, "native sand, 10 m", {cls: "t-label", anchor: "end", size: 28});
      K.dim(r.labs, sx, sy - 0.95 * m, sx + 4 * m, sy - 0.95 * m, "4 m", {labelColor: "var(--ink2)"});
      K.dim(r.labs, sx + 4.25 * m, sy, sx + 4.25 * m, sy + 2 * m, "2 m", {labelColor: "var(--ink2)"});
    }
    return r;
  }

  // ================================================================== scenes
  const scenes = [];

  // ---------------------------------------------------------------- 0. hook
  scenes.push({
    id: "hook", title: "A linear method, a nonlinear earthquake",
    build(s) {
      const g = K.g(s.svg);
      const S = 30, X = 470, YS = 700;
      K.soil(g, {x: 110, y: YS, w: 720, layers: [], hs: {h: 990 - YS, kind: "gravel"}, labels: false});
      s.zone = K.path(g, `M${X - 230},${YS} Q${X},${YS + 250} ${X + 230},${YS} Z`, {fill: "rgba(255,111,165,.24)", stroke: C.dash, "stroke-width": 3, "stroke-dasharray": "10 8", hidden: true});
      s.zoneT = K.text(g, X, YS + 90, "softened soil (illustration)", {cls: "t-label", anchor: "middle", color: "var(--ink)", size: 26, hidden: true});
      s.bx = box(g, X, YS, S);
      s.crk = K.g(g, {hidden: true});
      s.cracks = [crack(-140, 14, 5, 14, 18), crack(-40, 22, 5, 16, 18), crack(60, 10, 5, 12, 19), crack(150, 30, 4, -14, 18)]
        .map((pts) => ({pts, el: K.path(s.crk, "", {stroke: C.bad, "stroke-width": 4, fill: "none", "stroke-linejoin": "round"})}));
      s.lab = K.text(g, 110, 150, HLAB, {cls: "t-small"});
      s.fT = K.text(g, X, YS - 1.5 * S - 8 * S - 40, "sway: 8.0 Hz", {cls: "t-label", anchor: "middle", color: "var(--ink)", size: 30});
      s.head = K.heading(s.root, "SASSI is linear. Earthquakes aren't.", {x: 960, y: 140, w: 860, size: "h2"});
      s.lin = K.pill(s.root, "double the shaking → double the response", {x: 960, y: 290, color: "var(--ref)", size: 32});
      s.i1 = item(s.root, 960, 410, 860, "var(--dash)", "<b>Soil softens</b> in strong shaking");
      s.i2 = item(s.root, 960, 500, 860, "var(--bad)", "<b>Concrete walls crack</b>");
      s.fall = K.pill(s.root, "→ slower sway, shifted floor spectra", {x: 960, y: 610, color: "var(--ssi)", size: 34});
      s.run = 0; s.ph = 0; s.Cr = 0;
    },
    tick(s, t, dt) {
      advance(s, dt);
      sway(s.bx, s.ph, s.Cr);
      for (const c of s.cracks) c.el.setAttribute("d", K.d(c.pts.map(([lx, z]) => s.bx.b.T(lx, z))));
      s.fT.textContent = `sway: ${(HARM.el.f + (HARM.cr.f - HARM.el.f) * s.Cr).toFixed(1)} Hz`;
    },
    beats: [
      {say: "SASSI is a linear method: double the shaking, and every response simply doubles.",
        go(k) { k.s.run = 1; k.show(k.s.head); k.show(k.s.lin, {delay: 400}); }},
      {say: "[s]But real soil softens in strong shaking, [c]and concrete walls crack.",
        s(k) { k.show([k.s.i1, k.s.zone, k.s.zoneT]); },
        c(k) { k.show([k.s.i2, k.s.crk]); k.tween(k.s, {Cr: 1}, 1000); }},
      {say: "[f]Both make the building sway more slowly, and that shifts your floor spectra.", f(k) { k.show(k.s.fall); k.pulse(k.s.fT); }},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 9",
    build(s) {
      s.t = K.titleCard(s.root, {n: 9, part: "Advanced", title: "Nonlinear soil and cracked concrete by iteration",
        sub: "Guess the stiffness, run, check, repeat."});
    },
    beats: [
      {say: "Lesson nine. Nonlinear soil and cracked concrete, by iteration.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 100}); k.show(t.title, {delay: 200}); k.show(t.sub, {delay: 350}); k.show(t.rule, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. equivalent linear
  const rho = 19 / 9.81, G0 = rho * 250 * 250, GR = PH.sandRefStrain();        // the lesson figure's sand (kPa, strain in %)
  const BK = PH.hyperbolic(G0 / 100, GR);
  const ES = 240, EXS = 150;                     // the soil element: 1 m drawn 240 px, strain x 150
  scenes.push({
    id: "eqlin", title: "Soft when shaken hard",
    build(s) {
      const g = K.g(s.svg);
      K.defs();
      s.head = K.heading(s.root, "Soil softens as it shakes harder", {x: 110, y: 70, size: "h2"});
      // the loop
      s.p = K.plot(s.svg, {x: 230, y: 200, w: 640, h: 560, xr: [-0.2, 0.2], yr: [-60, 60], xticks: [-0.2, -0.1, 0, 0.1, 0.2], yticks: [-60, -30, 0, 30, 60],
        xlabel: "strain (%)", ylabel: "stress (kPa)", ylabelOffset: 84});
      s.loopF = K.path(s.p.data, "", {fill: "rgba(255,111,165,.30)", hidden: true});
      s.loop = K.path(s.p.data, "", {stroke: C.ssi, "stroke-width": 6, fill: "none"});
      s.sec = K.path(s.p.data, "", {stroke: C.spring, "stroke-width": 5, "stroke-dasharray": "14 9", fill: "none", hidden: true});
      s.dot = K.circle(s.p.data, 0, 0, 11, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      s.lab1 = K.label(s.root, "<span style='color:var(--spring)'>average slope = stiffness</span>", {x: 250, y: 214, size: 30, in: "fade"});
      s.lab2 = K.label(s.root, "<span style='color:var(--dash)'>loop area = damping</span>", {x: 250, y: 258, size: 30, in: "fade"});
      s.src = K.text(s.p.g, 230, 900, "hyperbolic sand model of the lesson figure (Vs 250 m/s)", {cls: "t-small"});
      // the real object: a 1 m soil element of that sand in cyclic simple shear (the same model)
      const EX0 = 1330, EY0 = 600;
      s.EX0 = EX0; s.EY0 = EY0;
      s.el = K.g(g, {hidden: true});
      K.ground(s.el, EX0 - 60, EY0, ES + 160);
      s.elem = K.path(s.el, "", {fill: "var(--sand)", "fill-opacity": 0.78, stroke: "rgba(255,240,215,.7)", "stroke-width": 3});
      s.elT = K.path(s.el, "", {fill: "url(#sv-pat-sand)"});
      s.tau = K.arrow(s.el, 0, 0, 1, 0, {color: C.ssi, width: 7, head: 20});
      K.text(s.el, EX0 - 60, 290, "stress", {cls: "t-label", color: "var(--ssi)", size: 28});
      K.text(s.el, EX0 + ES / 2, EY0 + 76, "a soil element, shaken in shear", {cls: "t-label", anchor: "middle", size: 30});
      K.text(s.el, EX0 + ES / 2, EY0 + 112, "strain × 150, the same model as the loop", {cls: "t-small", anchor: "middle"});
      // the spring and the damper that replace the loop (schematic), moving with the element's strain
      s.sd = K.g(g, {hidden: true});
      K.line(s.sd, 1150, 330, 1150, 600, {stroke: C.ref, "stroke-width": 5});
      for (let y = 340; y < 600; y += 22) K.line(s.sd, 1150, y, 1130, y + 18, {cls: "hatch"});
      K.line(s.sd, 1150, 600, 1760, 600, {stroke: C.ref, "stroke-width": 4});
      s.blk = K.rect(s.sd, 1500, 360, 200, 220, {rx: 8, fill: "var(--concrete-fill)", stroke: C.concrete, "stroke-width": 3});
      s.spr = K.spring(s.sd, 1150, 420, 1500, 420, {color: C.spring, coils: 7, amp: 20});
      s.dsh = K.dashpot(s.sd, 1500, 520, 1150, 520, {color: C.dash, w: 40, cyl: 110});
      K.text(s.sd, 1320, 382, "stiffness", {cls: "t-label", anchor: "middle", color: "var(--spring)", size: 32});
      K.text(s.sd, 1320, 576, "damping", {cls: "t-label", anchor: "middle", color: "var(--dash)", size: 32});
      K.text(s.sd, 1760, 640, "schematic", {cls: "t-small", anchor: "end"});
      s.eqT = K.pill(s.root, "equivalent-linear", {x: 1300, y: 760, color: "var(--good)", size: 36});
      s.amp = 0.6;
    },
    tick(s, t) {
      const xa = s.amp * GR, lp = PH.masing(BK, xa), p = s.p;
      const pts = [];
      for (let i = 0; i <= 60; i++) { const x = xa - 2 * xa * i / 60; pts.push([p.X(x), p.Y(lp.down(x))]); }
      for (let i = 0; i <= 60; i++) { const x = -xa + 2 * xa * i / 60; pts.push([p.X(x), p.Y(lp.up(x))]); }
      const d = K.d(pts);
      s.loop.setAttribute("d", d); s.loopF.setAttribute("d", d + "Z");
      s.sec.setAttribute("d", `M${p.X(-0.2)},${p.Y(-0.2 * lp.ksec)}L${p.X(0.2)},${p.Y(0.2 * lp.ksec)}`);
      const ph = TAU * t / 2.4, x = xa * Math.cos(ph), F = Math.sin(ph) >= 0 ? lp.down(x) : lp.up(x);
      s.dot.setAttribute("cx", p.X(x)); s.dot.setAttribute("cy", p.Y(F));
      // the element: base fixed, top displaced by the strain; the stress on its top face
      const dx = x / 100 * ES * EXS, x0 = s.EX0, y0 = s.EY0, e = `M${x0},${y0}L${x0 + ES},${y0}L${x0 + ES + dx},${y0 - ES}L${x0 + dx},${y0 - ES}Z`;
      s.elem.setAttribute("d", e); s.elT.setAttribute("d", e);
      const L = F * 2.2, cx = x0 + ES / 2 + dx;
      s.tau.update(cx - L / 2, y0 - ES - 34, cx + L / 2, y0 - ES - 34);
      s.tau.g.style.opacity = String(Math.min(1, Math.abs(L) / 24));
      // the block on its spring and damper, with the same displacement
      s.blk.setAttribute("x", 1500 + dx);
      s.spr.update(1150, 420, 1500 + dx, 420);
      s.dsh.update(1500 + dx, 520, 1150, 520);
    },
    beats: [
      {say: "Soil is like a foam mattress: firm when you press lightly, soft when you press hard.",
        go(k) { k.show(k.s.head); k.show(k.s.el, {delay: 300}); }},
      {say: "[l]Shake it back and forth, and its stress against strain traces a loop, not a straight line.",
        l(k) { k.show(k.s.p.g); }},
      {say: "[b]The harder the shaking, the more the loop leans over: [s]the soil gets softer.",
        b(k) { k.tween(k.s, {amp: 3.1}, 1000); }, s(k) { k.show([k.s.sec, k.s.lab1]); }},
      {say: "[f]And the fatter the loop, the more energy each cycle soaks up: that's damping.",
        f(k) { k.show([k.s.loopF, k.s.lab2]); }},
      {say: "[e]The trick: swap the loop for a spring with that average slope, plus a damper.",
        e(k) { k.hide(k.s.el); k.show(k.s.sd, {delay: 300}); }},
      {say: "[n]Engineers call this equivalent-linear: a linear model, with the right softness for the shaking.",
        n(k) { k.show(k.s.eqT); }},
    ],
  });

  // ---------------------------------------------------------------- 3. the iteration
  // The path of example 6's backfill on the library Sand curve: mean G/G_max of the 8 elements per pass
  // (NLSOIL_CONVERGENCE.TXT), each at the strain where the curve gives it; the start is GFAC = 0.4 x the
  // free-field G of layers L 11 and L 12 (Vs 245 and 230 m/s, 19 kN/m3) over the backfill's G_max (FILE74).
  const F74 = D["ex06/FILE74"].cols;
  const NLC = D["ex06/NLSOIL_CONVERGENCE.TXT"].cols;
  const SANDLG = PH.SAND.g.map(Math.log10);
  const sandG = (g) => lerp(SANDLG, PH.SAND.G, Math.log10(g));
  const sandInv = (r) => Math.pow(10, lerp(PH.SAND.G.slice().reverse(), SANDLG.slice().reverse(), r));
  const R_START = 0.4 * (19 / 9.81) * (245 * 245 + 230 * 230) / 2 / F74[6][0];
  scenes.push({
    id: "iterate", title: "Guess, run, check, repeat",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Guess, run, check, repeat", {x: 110, y: 70, size: "h2"});
      const CX = 640, CY = 590, R = 290;
      s.CX = CX; s.CY = CY; s.R = R;
      // chicken and egg
      s.ce = K.g(g, {hidden: true});
      K.curveArrow(s.ce, CX - 160, 560, CX + 160, 560, {color: C.muted, width: 5, head: 20, bend: -90});
      K.curveArrow(s.ce, CX + 160, 620, CX - 160, 620, {color: C.muted, width: 5, head: 20, bend: -90});
      s.ceP = [K.pill(s.root, "softness", {x: CX - 390, y: 560, color: "var(--spring)", size: 36, in: "fade"}),
        K.pill(s.root, "shaking", {x: CX + 180, y: 560, color: "var(--wave)", size: 36, in: "fade"})];
      // the cycle
      s.cyc = K.g(g, {hidden: true});
      K.circle(s.cyc, CX, CY, R, {fill: "none", stroke: "#2a3b52", "stroke-width": 6});
      [45, 135, 225, 315].forEach((d) => {
        const a = d * Math.PI / 180, px = CX + R * Math.sin(a), py = CY - R * Math.cos(a), tx = Math.cos(a), ty = Math.sin(a);
        K.arrow(s.cyc, px - 16 * tx, py - 16 * ty, px + 16 * tx, py + 16 * ty, {color: C.muted, width: 5, head: 22});
      });
      s.cycDot = K.circle(s.cyc, CX, CY - R, 15, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      const P = [[CX, CY - R], [CX + R, CY], [CX, CY + R], [CX - R, CY]];
      const names = ["1 · guess stiffness and damping", "2 · run the SSI", "3 · read the strains", "4 · update"];
      const cols = ["var(--spring)", "var(--wave)", "var(--ssi)", "var(--dash)"];
      const W = [470, 290, 330, 210];
      s.cycP = names.map((n, i) => K.pill(s.root, n, {x: P[i][0] - W[i] / 2, y: P[i][1] - 28, w: W[i], align: "center", color: cols[i], size: 30}));
      s.until = K.h("div", {x: CX - 130, y: CY - 50, w: 260, align: "center", size: 36, color: "var(--good)", in: "fade",
        html: "<b>until nothing<br>changes</b>"}, s.root);
      // the real iteration: example 6's backfill on its G/Gmax curve
      s.q = K.plot(s.svg, {x: 1200, y: 240, w: 600, h: 420, xr: [0.01, 0.1], xlog: true, yr: [0.4, 1.0], xticks: [0.01, 0.02, 0.05, 0.1], yticks: [0.4, 0.6, 0.8, 1.0],
        xlabel: "shear strain (%)", ylabel: "G / Gmax", ylabelOffset: 76});
      const gs = K.logspace(0.01, 0.1, 40);
      s.q.line(gs, gs.map(sandG), {color: C.spring, width: 5});
      s.q.text(0.035, sandG(0.035), "Sand curve", {cls: "t-label", anchor: "end", color: "var(--spring)", size: 26, dy: 50});
      const rs = NLC[8], gm = rs.map(sandInv), pts = [[0.01, R_START], [gm[0], R_START]];
      rs.forEach((r, i) => { pts.push([gm[i], r]); if (i + 1 < rs.length) pts.push([gm[i + 1], r]); });
      s.path = s.q.line(pts.map((q) => q[0]), pts.map((q) => q[1]), {color: C.wave, width: 5, draw: true});
      s.qL = K.g(s.q.g, {hidden: true});
      K.text(s.qL, s.q.X(0.0105), s.q.Y(R_START) - 14, "guess: run", {cls: "t-label", color: "var(--wave)", size: 26});
      K.text(s.qL, s.q.X(gm[0]) + 10, s.q.Y((R_START + rs[0]) / 2), "update", {cls: "t-label", color: "var(--dash)", size: 26});
      s.conv = K.g(s.q.g, {hidden: true});
      const last = rs.length - 1;
      K.circle(s.conv, s.q.X(gm[last]), s.q.Y(rs[last]), 11, {fill: C.good, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.conv, s.q.X(gm[last]) + 16, s.q.Y(rs[last]) - 22, `converged: ${rs[last].toFixed(2)}`, {cls: "t-label", color: "var(--good)", size: 26});
      s.qT = K.text(s.q.g, 1200, 214, "example 6's backfill, mean of its 8 elements", {cls: "t-small"});
      s.C = 0; s.ang = 0;
    },
    tick(s, t, dt) {
      if (s.C > 0) s.ang = (s.ang + dt * TAU / 4.0) % TAU;
      s.cycDot.setAttribute("cx", s.CX + s.R * Math.sin(s.ang)); s.cycDot.setAttribute("cy", s.CY - s.R * Math.cos(s.ang));
    },
    beats: [
      {say: "But there's a catch: the softness depends on the shaking, and the shaking depends on the softness.",
        go(k) { k.show(k.s.head); k.show([k.s.ce, k.s.ceP], {delay: 300}); }},
      {say: "[d]So SASSI iterates, like setting a shower: try a setting, feel the water, adjust, and repeat.",
        d(k) { k.show(k.s.q.g); k.draw(k.s.path, 1000, {delay: 200}); k.show(k.s.qL, {delay: 300}); }},
      {say: "[a]Guess the stiffness and damping. [b]Run the linear SSI analysis. [c]Read the strains. [u]Update from the soil's curves.",
        go(k) { k.hide([k.s.ce, k.s.ceP]); k.show(k.s.cyc, {delay: 200}); },
        a(k) { k.show(k.s.cycP[0]); }, b(k) { k.show(k.s.cycP[1]); }, c(k) { k.show(k.s.cycP[2]); }, u(k) { k.show(k.s.cycP[3]); }},
      {say: "[r]Then go around again, until the properties stop changing.",
        r(k) { k.tween(k.s, {C: 1}, 300); k.show([k.s.until, k.s.conv]); }},
    ],
  });

  // ---------------------------------------------------------------- 4. the backfill: a nonlinear soil group
  const GRATIO = (e) => F74[3][e - 1] / F74[6][e - 1];                 // converged G / Gmax of backfill element e
  const CELLS = [5, 6, 1, 2];                                          // the section cells (upper row at the wall first)
  const tint = (r, G) => (G > 0 ? `rgba(255,181,71,${(0.12 + G * (0.95 - r) * 1.15).toFixed(3)})` : "rgba(232,214,170,.55)");
  scenes.push({
    id: "backfill", title: "The soil next to the wall",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "The soil next to the wall", {x: 110, y: 70, size: "h2"});
      s.S = section(g, 380, 430, 95, {labels: true, depth: 4.3});
      s.code = K.code(s.root, ["PINGRP,3,1,0.4,1.0,Sand", "NLSSIITER,NLRUN+NLX+NLY+NLZ+NLCOMB,8"], {x: 1010, y: 170, w: 800, size: 28});
      s.cl = [K.label(s.root, "← the backfill: nonlinear sand", {x: 1010, y: 344, size: 28, color: "var(--muted)", in: "fade"}),
        K.label(s.root, "← the loop, up to 8 passes", {x: 1010, y: 388, size: 28, color: "var(--muted)", in: "fade"})];
      const ch = NLC[2].map(Math.abs);
      s.p = K.plot(s.svg, {x: 1150, y: 500, w: 640, h: 250, xr: [-0.6, 6.6], yr: [0, 70], xticks: [0, 1, 2, 3, 4, 5, 6], yticks: [0, 20, 40, 60],
        xlabel: "pass", ylabel: "largest change of G (%)", ylabelOffset: 70});
      s.bars = ch.map((v, i) => K.rect(s.p.data, s.p.X(i) - 30, s.p.Y(v), 60, s.p.Y(0) - s.p.Y(v),
        {fill: i === ch.length - 1 ? C.good : C.ssi, rx: 4, hidden: true, in: "fade"}));
      s.bl = K.g(s.p.g, {hidden: true});
      const last = ch.length - 1;
      K.text(s.bl, s.p.X(0) + 42, s.p.Y(ch[0]) + 26, `${Math.round(ch[0])} %`, {cls: "t-label", color: "var(--ink)", size: 28});
      K.text(s.bl, s.p.X(last), s.p.Y(Math.max(ch[last], 2)) - 16, `${ch[last].toFixed(1)} %`, {cls: "t-label", anchor: "middle", color: "var(--good)", size: 28});
      s.done = s.p.hline(2, {color: C.good, hidden: true});
      s.soft = K.pill(s.root, "a fifth to a quarter as stiff as the native soil", {x: 110, y: 900, color: "var(--ssi)", size: 32});
      s.G = 0;
    },
    tick(s) {
      s.S.cells.forEach((c, i) => c.setAttribute("fill", tint(GRATIO(CELLS[i]), s.G)));
    },
    beats: [
      {say: "SOIL handles the free field. But right next to a building, the soil can strain much more.",
        go(k) { k.show(k.s.head); }},
      {say: "[e]Take example 6: a wall carrying a 150-tonne bridge deck, with loose backfill behind it.",
        e(k) { k.show(k.s.S.labs); k.pulse(k.s.S.cells, {amp: 0.04}); }},
      {say: "[c]In SASSI, one line marks the backfill as nonlinear sand, [n]and NLSSIITER runs the loop for it.",
        c(k) { k.show(k.s.code.el); k.type(k.s.code.lines[0], {cps: 30}); k.show(k.s.cl[0], {delay: 400}); },
        n(k) { k.type(k.s.code.lines[1], {cps: 40}); k.show(k.s.cl[1], {delay: 400}); }},
      {say: "[p]It settles in 6 passes: each change is smaller than the last, down to 1.5 percent.",
        p(k) { k.show(k.s.p.g); k.show(k.s.bars, {stagger: 150, delay: 200}); k.show([k.s.bl, k.s.done], {delay: 400}); k.tween(k.s, {G: 1}, 1000); }},
      {say: "[s]The backfill ends up a fifth to a quarter as stiff as the native soil around it.",
        s(k) { k.show(k.s.soft); }},
    ],
  });

  // ---------------------------------------------------------------- 5. what the wall feels
  scenes.push({
    id: "wall", title: "What the wall feels",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "What the wall feels", {x: 110, y: 70, size: "h2"});
      const SX = 320, SY = 300, M = 70;
      s.S = section(g, SX, SY, M, {depth: 4});
      s.S.cells.forEach((c, i) => c.setAttribute("fill", tint(GRATIO(CELLS[i]), 1)));
      s.n22 = K.circle(g, SX, SY, 13, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3, hidden: true});
      s.n24 = K.circle(g, SX + 4 * M, SY, 13, {fill: C.wave, stroke: "#0a111d", "stroke-width": 3, hidden: true});
      s.st = [
        K.stat(s.root, {x: 110, y: 630, w: 400, value: "0.75", unit: "g", label: "wall top, with the deck", color: "var(--ssi)", vsize: 110}),
        K.stat(s.root, {x: 520, y: 630, w: 400, value: "0.58", unit: "g", label: "behind it: the free field", color: "var(--wave)", vsize: 110}),
      ];
      const ff = D["ex06/RS001_01.RS"], n22 = D["ex06/00022TR_X01.RS"];
      s.q = K.plot(s.svg, {x: 1150, y: 200, w: 660, h: 380, xr: [0.5, 50], xlog: true, yr: [0, 3], xticks: [0.5, 1, 2, 5, 10, 20, 50], yticks: [0, 1, 2, 3],
        xlabel: "frequency (Hz)", ylabel: "spectrum (g), 5 %", ylabelOffset: 62});
      s.lff = s.q.line(ff.f, ff.sa, {color: C.wave, width: 4, draw: true});
      s.l22 = s.q.line(n22.f, n22.sa, {color: C.ssi, width: 6, draw: true});
      s.leg = s.q.legend([{label: "wall top", color: C.ssi}, {label: "free field", color: C.wave}], {x: 1560, y: 236, size: 28, hidden: true});
      s.warn = K.card(s.root, {x: 110, y: 840, w: 860, kind: "warn", title: "Design the deck for the free field?", size: 32,
        body: "You'd miss almost a quarter of its acceleration."});
    },
    beats: [
      {say: "What does this change for the design? [t]The wall top, carrying the deck, now reaches 0.75 g.",
        go(k) { k.show(k.s.head); }, t(k) { k.show([k.s.n22, k.s.st[0].el]); }},
      {say: "[f]Behind the backfill, the ground stays at the free-field 0.58 g.",
        f(k) { k.show([k.s.n24, k.s.st[1].el]); }},
      {say: "[r]And the wall-top spectrum climbs well above the free field's.",
        r(k) { k.show([k.s.q.g, k.s.leg]); k.draw(k.s.lff, 1000); k.draw(k.s.l22, 1000, {delay: 300}); }},
      {say: "[w]Design the deck bearing for the free-field motion, and you'd miss almost a quarter of its acceleration.",
        w(k) { k.show(k.s.warn); }},
    ],
  });

  // ---------------------------------------------------------------- 6. walls that crack: Option NON
  const CRV = D["ex07/Panel0001.crv"].cols;                    // backbone of panel 1: strain, ..., force (kN)
  const EQL = D["ex07/Panel_EQL_Matl_Prop.txt"].cols;          // converged panels: E/E_el (col 5), x_eq (col 10)
  scenes.push({
    id: "walls", title: "Walls that crack: Option NON",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Walls that crack: Option NON", {x: 110, y: 70, size: "h2"});
      s.facts = K.label(s.root, "12 m box · two 4 m storeys · 0.3 m walls · shaken at 0.6 g", {x: 130, y: 172, size: 30, in: "fade"});
      s.lab = K.text(g, 130, 248, HLAB, {cls: "t-small"});
      const S = 40, X = 500, YS = 820;
      K.soil(g, {x: 130, y: YS, w: 740, layers: [], hs: {h: 990 - YS, kind: "gravel"}, labels: false});
      s.bx = box(g, X, YS, S);
      const tw = 0.3 * S;                                          // the end walls, edge-on: 0.3 m
      s.ends = [K.path(g, "", {fill: C.concrete, "fill-opacity": 0.55}), K.path(g, "", {fill: C.concrete, "fill-opacity": 0.55})];
      s.tw = tw;
      s.pan = K.g(g, {hidden: true});
      s.p1 = K.path(s.pan, "", {fill: "rgba(255,181,71,.22)", stroke: C.ssi, "stroke-width": 5});
      s.p2 = K.path(s.pan, "", {fill: "rgba(76,195,255,.14)", stroke: C.wave, "stroke-width": 5});
      s.p1T = K.text(s.pan, 0, 0, "panel · storey 1", {cls: "t-label", anchor: "middle", color: "var(--ssi)", size: 30});
      s.p2T = K.text(s.pan, 0, 0, "panel · storey 2", {cls: "t-label", anchor: "middle", color: "var(--wave)", size: 30});
      s.crk = K.g(g, {hidden: true});
      s.cracks = [crack(-200, 12, 5, 16, 18), crack(-80, 16, 5, 18, 18), crack(40, 10, 5, 14, 19), crack(170, 18, 5, -16, 18)]
        .map((pts) => ({pts, el: K.path(s.crk, "", {stroke: C.bad, "stroke-width": 5, fill: "none", "stroke-linejoin": "round"})}));
      // the backbone of panel 1 (Panel0001.crv: strain, force)
      const xs = CRV[0].map((x) => 100 * x), F = CRV[4].map((v) => v / 1000), kel = CRV[3][0];
      s.p = K.plot(s.svg, {x: 1150, y: 200, w: 660, h: 330, xr: [5e-4, 3], xlog: true, yr: [0, 14], xticks: [0.001, 0.01, 0.1, 1], yticks: [0, 4, 8, 12],
        xfmt: (v) => String(v), xlabel: "wall shear strain (%)", ylabel: "wall shear force (MN)", ylabelOffset: 62});
      const xUn = 100 * 11500 / kel, el = K.logspace(5e-4, xUn, 60);      // the uncracked slope G A_W, up to 11.5 MN
      s.lel = s.p.line(el, el.map((x) => kel * x / 100 / 1000), {color: C.ref, width: 3, draw: true});
      s.lelT = s.p.text(xUn, 11.5, "uncracked", {cls: "t-label", anchor: "end", color: "var(--ref)", size: 28, dx: -14, hidden: true});
      s.lbb = s.p.line(xs, F, {color: C.ssi, width: 7, draw: true});
      let ic = 0;                                                     // cracking: the last point on the elastic slope
      while (ic + 1 < CRV[0].length && CRV[4][ic + 1] / CRV[0][ic + 1] >= 0.99 * kel) ic++;
      s.crP = K.g(s.p.g, {hidden: true});
      K.circle(s.crP, s.p.X(xs[ic]), s.p.Y(F[ic]), 11, {fill: C.bad, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.crP, s.p.X(xs[ic]) - 18, s.p.Y(F[ic]) + 6, "cracks", {cls: "t-label", anchor: "end", color: "var(--bad)", size: 30});
      const xq = 100 * EQL[10][0], Fq = lerp(xs, F, xq);             // converged storey-1 panel: x_eq on the backbone
      s.cvP = K.g(s.p.g, {hidden: true});
      K.circle(s.cvP, s.p.X(xq), s.p.Y(Fq), 11, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.cvP, s.p.X(xq) + 16, s.p.Y(Fq) + 38, "converged, storey 1", {cls: "t-label", color: "var(--ssi)", size: 26});
      s.st = K.stat(s.root, {x: 1050, y: 680, w: 420, value: EQL[5][0].toFixed(2), label: "of the uncracked stiffness:<br>ground-floor walls", color: "var(--ssi)", vsize: 110});
      s.up = K.pill(s.root, "upper walls: uncracked", {x: 1470, y: 734, color: "var(--wave)", size: 30});
      s.run = 0; s.ph = 0; s.Cr = 0;
    },
    tick(s, t, dt) {
      advance(s, dt);
      sway(s.bx, s.ph, s.Cr);
      const b = s.bx.b, hw = s.bx.w / 2, top = 2 * s.bx.sh;
      s.ends[0].setAttribute("d", quad(b, -hw, -hw + s.tw, 0, top));
      s.ends[1].setAttribute("d", quad(b, hw - s.tw, hw, 0, top));
      s.p1.setAttribute("d", quad(b, -hw + 4, hw - 4, 2, s.bx.sh - 2));
      s.p2.setAttribute("d", quad(b, -hw + 4, hw - 4, s.bx.sh + 2, top - 2));
      const a = b.T(0, s.bx.sh - 24), c = b.T(0, top - 24);
      s.p1T.setAttribute("x", a[0]); s.p1T.setAttribute("y", a[1]);
      s.p2T.setAttribute("x", c[0]); s.p2T.setAttribute("y", c[1]);
      for (const k of s.cracks) k.el.setAttribute("d", K.d(k.pts.map(([lx, z]) => b.T(lx, z))));
    },
    beats: [
      {say: "Now the concrete. [b]Example 7: a two-storey shear-wall building, shaken at 0.6 g.",
        go(k) { k.show(k.s.head); k.s.run = 1; }, b(k) { k.show(k.s.facts); }},
      {say: "[p]ACS SASSI's Option NON treats each wall, storey by storey, as a panel.",
        p(k) { k.show(k.s.pan); }},
      {say: "[c]Each panel follows a backbone curve: force against deformation, steep until it cracks, then much flatter.",
        c(k) { k.show(k.s.p.g); k.draw(k.s.lel, 800); k.show(k.s.lelT, {delay: 300}); k.draw(k.s.lbb, 1000, {delay: 300}); k.show(k.s.crP, {delay: 400}); }},
      {say: "[r]After 7 passes, the ground-floor walls keep only 0.35 of their uncracked stiffness. [u]The upper walls never crack.",
        r(k) { k.show([k.s.st.el, k.s.crk, k.s.cvP]); k.tween(k.s, {Cr: 1}, 1000); }, u(k) { k.show(k.s.up); }},
    ],
  });

  // ---------------------------------------------------------------- 7. the spectra move
  scenes.push({
    id: "spectra", title: "The floor spectra move",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "The floor spectra move", {x: 110, y: 70, size: "h2"});
      const e113 = D["ex07/EL_00113TR_X01.RS"], c113 = D["ex07/00113TR_X01.RS"], r0 = rsPeak(e113), r1 = rsPeak(c113);
      s.q = K.plot(s.svg, {x: 240, y: 200, w: 900, h: 440, xr: [1, 50], xlog: true, yr: [0, 8], xticks: [1, 2, 5, 10, 20, 50], yticks: [0, 2, 4, 6, 8],
        xlabel: "frequency (Hz)", ylabel: "roof spectrum (g), 5 %", ylabelOffset: 62});
      s.le = s.q.line(e113.f, e113.sa, {color: C.ref, width: 5, draw: true});
      s.lc = s.q.line(c113.f, c113.sa, {color: C.ssi, width: 7, draw: true});
      s.leg = s.q.legend([{label: "uncracked walls", color: C.ref}, {label: "cracked walls", color: C.ssi}], {x: 870, y: 470, size: 28, hidden: true});
      s.pk0 = K.g(s.q.g, {hidden: true});
      K.circle(s.pk0, s.q.X(r0.f), s.q.Y(r0.sa), 11, {fill: C.ref, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.pk0, s.q.X(r0.f) + 20, s.q.Y(r0.sa) + 6, `${r0.f.toFixed(2)} Hz`, {cls: "t-label", color: "var(--ref)", size: 30});
      s.pk1 = K.g(s.q.g, {hidden: true});
      K.circle(s.pk1, s.q.X(r1.f), s.q.Y(r1.sa), 11, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      K.text(s.pk1, s.q.X(r1.f) - 20, s.q.Y(r1.sa) + 6, `${r1.f.toFixed(2)} Hz`, {cls: "t-label", anchor: "end", color: "var(--ssi)", size: 30});
      const ya = s.q.Y(Math.max(r0.sa, r1.sa)) - 26;
      s.shift = K.arrow(s.q.g, s.q.X(r0.f) - 6, ya, s.q.X(r1.f) + 6, ya, {color: C.ssi, width: 5, head: 18, hidden: true});
      // two buildings: uncracked and cracked, the computed sway at the same scale
      s.bl = K.g(g, {hidden: true});
      const S = 15, YS = 600;
      K.line(s.bl, 1290, YS, 1810, YS, {stroke: C.gravel, "stroke-width": 5});
      s.b1 = box(s.bl, 1420, YS, S, {color: "var(--ref)"});
      s.b2 = box(s.bl, 1680, YS, S);
      s.b2.b.polys[0].style.fill = "rgba(255,181,71,.3)";
      const yt = YS - 1.5 * S - 8 * S - 30;
      K.text(s.bl, 1420, yt, "uncracked", {cls: "t-label", anchor: "middle", color: "var(--ref)", size: 30});
      K.text(s.bl, 1680, yt, "cracked", {cls: "t-label", anchor: "middle", color: "var(--ssi)", size: 30});
      K.text(s.bl, 1420, YS + 50, `faster: ${HARM.el.f.toFixed(1)} Hz`, {cls: "t-label", anchor: "middle", size: 28});
      K.text(s.bl, 1680, YS + 50, `slower: ${HARM.cr.f.toFixed(1)} Hz`, {cls: "t-label", anchor: "middle", size: 28});
      K.text(s.bl, 1810, YS + 96, "computed sway, × 250, slowed down 16 ×", {cls: "t-small", anchor: "end"});
      s.warn = K.card(s.root, {x: 240, y: 790, w: 900, kind: "warn", title: "One linear run", size: 32,
        body: "keeps uncracked walls, though its own strains say they cracked."});
      s.cap = K.card(s.root, {x: 1240, y: 300, w: 570, kind: "warn", title: "One caution: damping", size: 32,
        body: "Cracked walls came out at <b>16.7 %</b>. ASCE 4 accepts <b>7 %</b>: cap it."});
      s.run = 0; s.ph1 = 0; s.ph2 = 0;
    },
    tick(s, t, dt) {
      if (s.run) { s.ph1 += TAU * HARM.el.f * dt / SLOW; s.ph2 += TAU * HARM.cr.f * dt / SLOW; }
      sway(s.b1, s.ph1, 0);
      sway(s.b2, s.ph2, 1);
    },
    beats: [
      {say: "Why should you care? A softer building sways more slowly, [b]so its floor spectra shift.",
        go(k) { k.show(k.s.head); }, b(k) { k.show(k.s.bl); k.s.run = 1; }},
      {say: "[u]Uncracked, the roof spectrum peaks at 7.94 hertz. [c]Cracked, it peaks at 6.76 hertz, and a little lower.",
        u(k) { k.show([k.s.q.g, k.s.leg]); k.draw(k.s.le, 1000); k.show(k.s.pk0, {delay: 400}); },
        c(k) { k.draw(k.s.lc, 1000); k.show([k.s.pk1, k.s.shift], {delay: 400}); }},
      {say: "[l]A single linear run would keep uncracked walls, even though its own strains say they have cracked.",
        l(k) { k.show(k.s.warn); }},
      {say: "[c]One caution: the cracked walls also came out with 16.7 percent damping. ASCE 4 accepts 7, so cap it.",
        c(k) { k.hide(k.s.bl); k.show(k.s.cap, {delay: 300}); }},
    ],
  });

  // ---------------------------------------------------------------- 8. recap
  scenes.push({
    id: "recap", title: "Recap",
    build(s) {
      s.head = K.heading(s.root, "Recap", {x: 110, y: 80, size: "h1"});
      s.list = K.bullets(s.root, [
        {t: "Linear SASSI, nonlinear effects: <b>iterate</b>", sub: "guess, run, check, repeat"},
        {t: "<b>Soil next to a structure</b> may need it", sub: "wall top: 0.75 g, not the free-field 0.58 g"},
        {t: "<b>Option NON</b> cracks the walls", sub: "roof peak: from 7.94 to 6.76 Hz"},
      ], {x: 110, y: 230, w: 980, num: true, size: 42});
      s.q = K.card(s.root, {x: 1170, y: 230, w: 640, kind: "check", title: "Check yourself", size: 34,
        body: "Why does SASSI need several passes, not just one?"});
      s.a = K.card(s.root, {x: 1170, y: 520, w: 640, title: "Answer", size: 34,
        body: "Softness depends on shaking, and shaking on softness: each pass refines the last."});
      s.next = K.pill(s.root, "Next · Lesson 10: good practice, checks and the limits of SSI models →", {x: 110, y: 900, color: "var(--wave)", size: 32});
    },
    beats: [
      {say: "Let's recap. [a]SASSI is linear, so it handles soft soil and cracked walls by iterating: guess, run, check, repeat.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[a]Soil next to a structure may need its own loop: here the wall top reached 0.75 g, not 0.58.",
        a(k) { k.show(k.s.list.items[1]); }},
      {say: "[a]Option NON cracks the walls panel by panel: the roof spectrum peak moved from 7.94 to 6.76 hertz.",
        a(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Check yourself. Why does SASSI need several passes, and not just one?",
        q(k) { k.show(k.s.q); }, gap: 1500},
      {say: "[a]Because the softness depends on the shaking, and the shaking on the softness: each pass refines the last, until they agree.",
        a(k) { k.show(k.s.a); }},
      {say: "[n]Next, in lesson 10: good practice, checks, and the limits of SSI models.", n(k) { k.show(k.s.next); }},
    ],
  });

  SV.video({
    id: "09", n: 9, part: "Advanced", lesson: "09-nonlinear",
    title: "Nonlinear soil and cracked concrete by iteration",
    subtitle: "Equivalent-linear iterations in a linear program: a loose backfill behind a wall (NLSSIITER), and cracking shear walls (Option NON).",
    next: {href: "10.html", title: "Good practice, checks and the limits of SSI models"},
    pron: [[/\bOption NON\b/g, "option non"]],
    scenes,
  });
})();
