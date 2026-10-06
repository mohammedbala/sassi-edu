/* Explainer video, lesson 3 (newcomer edition): Foundation impedance -- soil springs and dashpots.
 * One big idea: the soil under a foundation acts like springs and shock absorbers that change with frequency;
 * SASSI computes them from the interaction nodes instead of taking them from a formula; radiation damping is
 * energy leaving as waves.
 * Numbers: sassi/ui/lessons/03_impedance.md; curves: d03.js (the lesson's runs, python -m sassi.ui.video_data):
 * FOUNSTIF / FOUNDAMP of the 2 m mesh, the Ricker load and the mat's response.
 * Geometry to scale: the 12 m mat of example 3 (7 x 7 interaction nodes at 2 m); the hook's building is the
 * example 1 stick (20 m mat, four 5 m storeys), moving with the one-mode SSI model of lessons 1 and 4
 * (PH.ssiResponse) at its SSI frequency.  Spring and dashpot symbols, waves and the rigid-body motions are
 * labelled schematic / illustration on screen. */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["03"], PH = K.PH;
  const TAU = 2 * Math.PI;
  const ST = D["ex03/FOUNSTIF"].cols, DA = D["ex03/FOUNDAMP"].cols;
  const DISP = D["ex03/00025TR_Z.ACC"], LOAD = D["data/ricker_5hz.th"];
  const FS = ST[0];                                          // the 16 analysis frequencies (Hz)
  const IJ = (r, c) => 1 + 6 * r + c;                        // column of K(r, c) in the FOUN* files (rows X Y Z XX YY ZZ)
  const G = 80000, B = 6;                                    // soil shear modulus (kPa), half-width of the mat (m)
  /** Mode 1 of the example 1 stick: Gamma1 phi1 at the floors (5, 10, 15, 20 m); 1.34 at the roof. */
  const GPHI = [0.2131, 0.5546, 0.9521, 1.3428];

  // ---------------------------------------------------------------- shared pieces
  const re = (c, ph) => c[0] * Math.cos(ph) - c[1] * Math.sin(ph);      // Re(c e^{i ph})
  /** Linear interpolation in a {t, v} history. */
  function interp(h, t) {
    const T = h.t, V = h.v;
    if (t <= T[0]) return V[0];
    if (t >= T[T.length - 1]) return V[V.length - 1];
    let lo = 0, hi = T.length - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (T[m] <= t) lo = m; else hi = m; }
    const a = (t - T[lo]) / ((T[hi] - T[lo]) || 1);
    return V[lo] + a * (V[hi] - V[lo]);
  }
  function halo(e, w) {
    e.style.paintOrder = "stroke"; e.style.stroke = "rgba(10,17,29,.85)"; e.style.strokeWidth = (w || 6) + "px"; e.style.strokeLinejoin = "round";
    return e;
  }
  const note = (p, x, y, str, o) => halo(K.text(p, x, y, str, Object.assign({cls: "t-small", size: 23}, o || {})), 5);
  /** A square mat in plan with n x n nodes over `size` px: {g, pts, dots}. */
  function plan(p, x, y, size, n, o) {
    o = o || {};
    const g = K.g(p, {hidden: o.hidden});
    K.rect(g, x, y, size, size, {fill: "rgba(59,74,94,.75)", stroke: C.concrete, "stroke-width": 3});
    const pts = [];
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) pts.push([x + (size * i) / (n - 1), y + size - (size * j) / (n - 1)]);
    const nd = K.nodes(g, pts, {r: o.r || 8, color: C.wave});
    return {g, pts, dots: nd.dots};
  }
  /** A small fixed wall (vertical hatch) at x from y1 to y2, hatched on the side `dir` (-1 left, +1 right). */
  function wall(p, x, y1, y2, dir) {
    const g = K.g(p, {});
    K.line(g, x, y1, x, y2, {stroke: C.ref, "stroke-width": 5, "stroke-linecap": "round"});
    for (let y = y1 + 6; y < y2; y += 18) K.line(g, x, y, x + 16 * dir, y + 14, {cls: "hatch"});
    return g;
  }
  /** A point of a stick's mat (local lx, ly from the mat centre, ly + down) in its current sway / rock. */
  function matPt(st, x0, y0, matH, lx, ly) {
    const s = st.state, cx = x0 + s.sway, cy = y0 + matH / 2, c = Math.cos(s.rock), sn = Math.sin(s.rock);
    return [cx + lx * c - ly * sn, cy + lx * sn + ly * c];
  }

  // ================================================================== scenes
  const scenes = [];

  // ---------------------------------------------------------------- 0. hook
  // to scale, 14 px per metre: the example 1 building; motion relative to the ground at 3.49 Hz (its SSI
  // frequency), 1 mm of ground motion drawn x 600 = 8.4 px, 5 x slower
  const HP = 14, HA = 8.4, HF = 3.49, MH = 1.5 * HP;
  scenes.push({
    id: "hook", title: "Springs from a formula",
    build(s) {
      const g = K.g(s.svg);
      s.ssi = PH.ssiResponse(PH.ssiModel(1), HF);
      s.head = K.heading(s.root, "Where do soil springs come from?", {x: 120, y: 70});
      const stick = (p, x) => K.stick(p, {x, y: 760 - MH, z: [5 * HP, 10 * HP, 15 * HP, 20 * HP], matW: 20 * HP, matH: MH, r: 15, slabW: 96});
      // left: the ANSYS model on springs and dashpots
      s.L = K.g(g, {hidden: true});
      K.ground(s.L, 250, 840, 400);
      wall(s.L, 250, 728, 790, -1); wall(s.L, 650, 728, 790, 1);
      s.sL = stick(s.L, 450);
      s.vS = [-125, 125].map(() => K.spring(s.L, 0, 0, 0, 0, {coils: 4, amp: 11}));
      s.vD = [-88, 88].map(() => K.dashpot(s.L, 0, 0, 0, 0, {w: 22, cyl: 34}));
      s.hS = K.spring(s.L, 0, 0, 0, 0, {coils: 4, amp: 10});
      s.hD = K.dashpot(s.L, 0, 0, 0, 0, {w: 22, cyl: 34});
      note(s.L, 450, 880, "springs and dashpots, schematic", {anchor: "middle"});
      s.f1 = K.card(s.root, {x: 150, y: 170, w: 600, kind: "ansys", title: "In ANSYS: a textbook formula", body: "<div class='v03-f1'></div>"});
      K.tex(s.f1.querySelector(".v03-f1"), String.raw`K = \frac{8GR}{2-\nu}`, false);
      s.f1.querySelector(".v03-f1").style.fontSize = "44px";
      s.one = K.pill(s.root, "one number, the same at every frequency", {x: 200, y: 910, color: "var(--ref)", size: 30});
      // right: the same building on real soil
      s.R = K.g(g, {hidden: true});
      K.soil(s.R, {x: 1040, y: 760, w: 760, layers: [{h: 100, kind: "sand"}], hs: {h: 130, kind: "sand"}, labels: false});
      s.rip = K.ripples(s.R, {cx: 1420, cy: 762, r0: 30, r1: 330, n: 4, color: C.dash, squash: 0.5, width: 3});
      s.glue = [K.path(s.R, "", {fill: C.sand, "fill-opacity": 0.78}), K.path(s.R, "", {fill: "url(#sv-pat-sand)"})];   // soil stuck to the mat
      s.sR = stick(s.R, 1420);
      K.text(s.R, 1060, 730, "real soil", {cls: "t-label", size: 32, color: "var(--ink)"});
      note(s.R, 1800, 990, "waves leaving the mat, schematic", {anchor: "end"});
      s.q = K.text(g, 1690, 470, "?", {cls: "t-big", size: 170, anchor: "middle", color: "var(--ssi)", hidden: true, in: "pop"});
      s.scale = K.g(g, {hidden: true});
      note(s.scale, 1800, 240, "example 1 building at 3.49 Hz, relative to the ground", {anchor: "end"});
      note(s.scale, 1800, 270, "1 mm of ground motion drawn × 600 · 5 × slower", {anchor: "end"});
      s.A = 0; s.B = 0;
    },
    tick(s, t) {
      const ph = TAU * HF * t / 5, r = s.ssi;
      const pose = (st, A) => st.set({sway: A * re(r.u0, ph), rock: (A / HP) * re(r.th, ph), bend: GPHI.map((gp) => A * gp * re(r.u, ph))});
      pose(s.sL, HA * s.A);
      pose(s.sR, HA * s.B);
      const P = (lx, ly) => matPt(s.sL, 450, 760 - MH, MH, lx, ly);
      [-125, 125].forEach((lx, i) => { const q = P(lx, MH / 2); s.vS[i].update(q[0], q[1], 450 + lx, 840); });
      [-88, 88].forEach((lx, i) => { const q = P(lx, MH / 2); s.vD[i].update(q[0], q[1], 450 + lx, 840); });
      const ql = P(-10 * HP, 0), qr = P(10 * HP, 0);
      s.hS.update(250, ql[1], ql[0], ql[1]);
      s.hD.update(qr[0], qr[1], 650, qr[1]);
      const st = s.sR.state, c = Math.cos(st.rock), sn = Math.sin(st.rock), h = MH / 2;
      const Q = (lx) => [1420 + st.sway + lx * c - h * sn, 760 - h + lx * sn + h * c];
      const a = Q(-10 * HP), b = Q(10 * HP), d = `M${a}L${b}L${b[0]},760L${a[0]},760Z`;
      s.glue.forEach((e) => e.setAttribute("d", d));
      s.rip.amp = s.B;
      s.rip.set(t, 1.6);
    },
    beats: [
      {say: "If you've ever put a building on springs in ANSYS, [f]you took those springs from a textbook formula.",
        go(k) { k.show([k.s.head, k.s.L, k.s.scale]); k.tween(k.s, {A: 1}, 1000); }, f(k) { k.show(k.s.f1); }},
      {say: "[o]One number per direction, the same at every frequency.",
        o(k) { k.show(k.s.one); }},
      {say: "[s]But how does real soil push back on a shaking foundation? [q]And is it really just one number?",
        s(k) { k.show(k.s.R); k.tween(k.s, {B: 1}, 1000); }, q(k) { k.show(k.s.q); }},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 3",
    build(s) {
      s.t = K.titleCard(s.root, {n: 3, part: "Fundamentals", title: "Foundation impedance: soil springs and dashpots",
        sub: "Why the soil acts like a suspension, and how SASSI computes it."});
    },
    beats: [
      {say: "Lesson three: foundation impedance, the soil's springs and dashpots.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 100}); k.show(t.title, {delay: 200}); k.show(t.rule, {delay: 300}); k.show(t.sub, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. springs and dashpots
  // the 12 m mat of example 3 to scale (30 px per metre); the springs and dashpots are symbols
  scenes.push({
    id: "suspension", title: "Springs and dashpots",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Springs and dashpots", {x: 120, y: 70});
      const X0 = 960, Y0 = 470, MT = 18;                       // mat centre x, mat top y, mat thickness (px)
      s.soilG = K.g(g, {hidden: true});
      K.soil(s.soilG, {x: 560, y: Y0 + MT, w: 800, layers: [{h: 230, kind: "sand"}], labels: false}).g.style.opacity = "0.4";
      s.gnd = K.ground(g, 700, 640, 520, {color: C.ref});
      s.mat = K.rect(g, X0 - 180, Y0, 360, MT, {cls: "mat"});
      K.dim(g, X0 - 180, Y0 - 34, X0 + 180, Y0 - 34, "");
      K.text(g, X0, Y0 - 50, "12 m mat", {cls: "t-label", anchor: "middle", size: 30, color: "var(--ink)"});
      const yb = Y0 + MT;
      // one pair, as in a car's suspension
      s.one = K.g(g, {hidden: true});
      s.oneS = K.spring(s.one, X0 - 30, yb, X0 - 30, 640, {coils: 6, amp: 15, hidden: true});
      s.oneD = K.dashpot(s.one, X0 + 30, yb, X0 + 30, 640, {w: 30, cyl: 70, hidden: true});
      s.oneSt = K.text(s.one, X0 - 56, 580, "spring", {cls: "t-label", anchor: "end", size: 32, color: "var(--spring)", hidden: true});
      s.oneDt = K.text(s.one, X0 + 56, 580, "shock absorber", {cls: "t-label", size: 32, color: "var(--dash)", hidden: true});
      // the full set: vertical pairs under the edges (vertical, rocking), a horizontal pair at the mat (sliding)
      s.full = K.g(g, {hidden: true});
      for (const [xs, xd] of [[X0 - 160, X0 - 110], [X0 + 110, X0 + 160]]) {
        K.spring(s.full, xs, yb, xs, 640, {coils: 6, amp: 14});
        K.dashpot(s.full, xd, yb, xd, 640, {w: 26, cyl: 60});
      }
      wall(s.full, 690, Y0 - 30, Y0 + 50, -1); wall(s.full, 1230, Y0 - 30, Y0 + 50, 1);
      K.spring(s.full, 690, Y0 + MT / 2, X0 - 180, Y0 + MT / 2, {coils: 5, amp: 13});
      K.dashpot(s.full, X0 + 180, Y0 + MT / 2, 1230, Y0 + MT / 2, {w: 26, cyl: 30});
      s.labF = K.g(g, {hidden: true});
      halo(K.text(s.labF, X0 - 184, 600, "springs", {cls: "t-label", anchor: "end", size: 30, color: "var(--spring)"}));
      halo(K.text(s.labF, X0 + 184, 600, "dashpots", {cls: "t-label", size: 30, color: "var(--dash)"}));
      s.sch = note(g, X0, 690, "the soil as springs and dashpots: schematic symbols", {anchor: "middle", hidden: true});
      s.term = K.pill(s.root, "together: the foundation's <b>impedance</b>", {x: 1260, y: 170, color: "var(--wave)", size: 32});
      // the three rigid-body motions of the mat, each with its own springs and dashpots (illustration)
      s.mot = ["sliding", "vertical", "rocking"].map((name, i) => {
        const cx = 400 + 560 * i, og = K.g(g, {hidden: true}), mg = K.g(og, {});
        K.ground(og, cx - 120, name === "sliding" ? 878 : 930, 240, {color: C.muted});
        K.rect(mg, cx - 70, 852, 140, 12, {cls: "mat"});
        const parts = [];
        if (name === "sliding") {
          wall(og, cx - 130, 836, 880, -1); wall(og, cx + 130, 836, 880, 1);
          parts.push(K.spring(og, 0, 0, 0, 0, {coils: 3, amp: 9}), K.dashpot(og, 0, 0, 0, 0, {w: 18, cyl: 22}));
          for (const x of [cx - 40, cx + 40]) K.circle(og, x, 870, 6, {fill: "none", stroke: C.muted, "stroke-width": 2.5});
        } else if (name === "vertical") {
          parts.push(K.spring(og, 0, 0, 0, 0, {coils: 4, amp: 9}), K.dashpot(og, 0, 0, 0, 0, {w: 18, cyl: 30}));
        } else {
          parts.push(K.spring(og, 0, 0, 0, 0, {coils: 4, amp: 8}), K.spring(og, 0, 0, 0, 0, {coils: 4, amp: 8}),
            K.dashpot(og, 0, 0, 0, 0, {w: 16, cyl: 28}), K.dashpot(og, 0, 0, 0, 0, {w: 16, cyl: 28}));
        }
        K.text(og, cx, 985, name, {cls: "t-label", anchor: "middle", size: 32, color: "var(--ink)"});
        return {g: og, mg, cx, name, parts};
      });
      s.illu = note(g, 1830, 790, "rigid-body motions: illustration", {anchor: "end", hidden: true});
      s.M = 0;
    },
    tick(s, t) {
      const m = Math.sin(TAU * t / 1.6) * s.M;
      for (const o of s.mot) {
        const cx = o.cx;
        if (o.name === "sliding") {
          const dx = 18 * m;
          o.mg.setAttribute("transform", `translate(${dx},0)`);
          o.parts[0].update(cx - 130, 858, cx - 70 + dx, 858);
          o.parts[1].update(cx + 70 + dx, 858, cx + 130, 858);
        } else if (o.name === "vertical") {
          const dy = 10 * m;
          o.mg.setAttribute("transform", `translate(0,${dy})`);
          o.parts[0].update(cx - 30, 864 + dy, cx - 30, 930);
          o.parts[1].update(cx + 30, 864 + dy, cx + 30, 930);
        } else {
          const a = 0.11 * m, c = Math.cos(a), sn = Math.sin(a);
          o.mg.setAttribute("transform", `rotate(${(a * 180 / Math.PI).toFixed(3)} ${cx} 858)`);
          const P = (lx) => [cx + lx * c - 6 * sn, 858 + lx * sn + 6 * c];
          [[-60, 0], [60, 1]].forEach(([lx, i]) => { const q = P(lx); o.parts[i].update(q[0], q[1], cx + lx, 930); });
          [[-32, 2], [32, 3]].forEach(([lx, i]) => { const q = P(lx); o.parts[i].update(q[0], q[1], cx + lx, 930); });
        }
      }
    },
    beats: [
      {say: "Think of a car's suspension: [sp]a spring, and [d]a shock absorber next to it.",
        go(k) { k.show([k.s.head, k.s.one]); }, sp(k) { k.show([k.s.oneS.g, k.s.oneSt]); }, d(k) { k.show([k.s.oneD.g, k.s.oneDt]); }},
      {say: "[m]The soil under a foundation behaves the same way.",
        m(k) { k.hide(k.s.one); k.show([k.s.soilG, k.s.full, k.s.sch]); }},
      {say: "[n]Engineers call them springs and dashpots. [t]Together, they're the foundation's impedance.",
        n(k) { k.show(k.s.labF); }, t(k) { k.show(k.s.term); }},
      {say: "[x]A mat has a set for each way it can move: [a]sliding, [b]up and down, [c]and rocking.",
        x(k) { k.show(k.s.illu); k.tween(k.s, {M: 1}, 800); }, a(k) { k.show(k.s.mot[0].g); }, b(k) { k.show(k.s.mot[1].g); }, c(k) { k.show(k.s.mot[2].g); }},
    ],
  });

  // ---------------------------------------------------------------- 3. push the mat: radiation damping
  // the 12 m mat to scale (20 px per metre); its computed vertical motion drawn x 10 000, 2.5 x slower
  const RP = 20, RX = 10000;
  scenes.push({
    id: "radiation", title: "Where the energy goes",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Push it: where does the energy go?", {x: 120, y: 70});
      K.soil(g, {x: 120, y: 640, w: 680, layers: [{h: 190, kind: "sand"}], hs: {h: 150, kind: "sand"}, labels: false}).layers[0].edge.style.display = "none";
      // the top 2 m of soil, whose surface follows the mat (welded contact; decay 1/r outside it, as in statics)
      s.top = [K.path(g, "", {fill: C.sand, "fill-opacity": 0.78}), K.path(g, "", {fill: "url(#sv-pat-sand)"})];
      s.soilLab = halo(K.text(g, 140, 900, "uniform soil · Vs 200 m/s · 2 % damping", {cls: "t-label", size: 28, color: "var(--ink)"}));
      // the radiated wave: two fronts leaving the mat after the pulse (schematic, not to speed)
      s.fronts = [0, 1].map(() => K.path(g, "", {stroke: C.dash, "stroke-width": 4, fill: "none"}));
      s.mat = K.rect(g, 460 - 6 * RP, 578, 12 * RP, 22, {cls: "mat"});
      K.text(g, 322, 530, "12 m mat", {cls: "t-label", anchor: "end", size: 28, color: "var(--ink)"});
      s.arr = K.arrow(g, 460, 560, 460, 440, {color: C.ssi, width: 7, head: 22});
      note(g, 800, 960, "mat motion drawn × 10 000 · 2.5 × slower", {anchor: "end"});
      s.wNote = note(g, 800, 990, "waves: schematic, not to speed", {anchor: "end", hidden: true, color: "var(--dash)"});
      // the push and the response (computed)
      const PL = K.plot(s.svg, {x: 1000, y: 160, w: 800, h: 220, xr: [0, 2], yr: [-600, 1200], xticks: [0, 0.5, 1, 1.5, 2], yticks: [0, 1000],
        ylabel: "push (kN)", ylabelOffset: 90});
      const PD = K.plot(s.svg, {x: 1000, y: 500, w: 800, h: 330, xr: [0, 2], yr: [-0.15, 0.32], xticks: [0, 0.5, 1, 1.5, 2], yticks: [-0.1, 0, 0.1, 0.2],
        xlabel: "time (s)", ylabel: "mat moves (mm)", ylabelOffset: 90});
      s.PL = PL; s.PD = PD;
      const lt = LOAD.t.filter((t) => t <= 2), lv = lt.map((t, i) => 1000 * LOAD.v[i]);
      s.lLoad = PL.line(lt, lv, {color: C.ssi, width: 5, draw: true});
      const dt = [], dv = [];
      DISP.t.forEach((t, i) => { if (t <= 2) { dt.push(t); dv.push(1000 * DISP.v[i]); } });
      s.lDisp = PD.line(dt, dv, {color: C.wave, width: 6, draw: true});
      s.stop = K.pill(s.root, "stops almost at once", {x: PD.X(1.0), y: PD.Y(0) - 120, color: "var(--good)", size: 32});
      s.term = K.pill(s.root, "radiation damping", {x: 580, y: 420, color: "var(--dash)", size: 34});
      s.headL = PL.vline(0, {color: C.ink, width: 2, dash: "4 6"});
      s.headD = PD.vline(0, {color: C.ink, width: 2, dash: "4 6"});
      s.ph = 0; s.R = 0; s.Rp = 0;
    },
    tick(s, t, dt) {
      s.ph += (dt || 0) * 0.4 * s.R;
      const td = s.R ? Math.min(2, s.ph % 2.6) : 0;
      const u = 1000 * interp(DISP, td), F = 1000 * interp(LOAD, Math.min(td, 2));
      const up = -u * RP * RX / 1000;                         // mm -> px (x 10 000)
      s.mat.setAttribute("y", 578 + up);
      const pts = [];
      for (let x = 120; x <= 800; x += 10) { const dx = Math.abs(x - 460); pts.push([x, 600 + up * (dx <= 120 ? 1 : 120 / dx)]); }
      const d = "M120,640" + K.d(pts).replace("M", "L") + "L800,640Z";
      s.top.forEach((e) => e.setAttribute("d", d));
      s.arr.update(460, 572 + up, 460, 572 + up - 0.13 * F);
      s.arr.g.style.opacity = String(Math.min(1, Math.abs(F) / 60));
      s.fronts.forEach((f, i) => {
        const age = td - 0.45 - 0.12 * i, r = 30 + 260 * age;
        if (age <= 0 || r > 330) { f.setAttribute("d", ""); return; }
        f.setAttribute("d", `M${460 - r},602A${r},${r * 0.62} 0 0 0 ${460 + r},602`);
        f.style.opacity = String(s.Rp * 0.9 * (1 - r / 330) * (i ? 0.6 : 1));
      });
      [s.headL, s.headD].forEach((h, i) => { const P = i ? s.PD : s.PL; h.setAttribute("x1", P.X(td)); h.setAttribute("x2", P.X(td)); });
    },
    beats: [
      {say: "Take a rigid square mat, 12 metres wide, on uniform soil. [p]Give it one short, sharp push upwards.",
        go(k) { k.show(k.s.head); }, p(k) { k.show(k.s.PL.g); k.draw(k.s.lLoad, 1000); }},
      {say: "[m]The mat jumps up, [s]and then stops almost at once. No long ringing.",
        m(k) { k.show(k.s.PD.g); k.draw(k.s.lDisp, 1000); k.tween(k.s, {R: 1}, 10); }, s(k) { k.show(k.s.stop); }},
      {say: "The soil's own damping is only 2 percent. That can't stop it so fast.",
        go(k) { k.pulse(k.s.soilLab, {amp: 0.08}); k.s.soilLab.style.fill = "var(--ssi)"; }},
      {say: "[r]Instead, the mat sends waves out into the ground, like ripples from a stone in a pond.",
        r(k) { k.s.soilLab.style.fill = "var(--ink)"; k.show(k.s.wNote); k.tween(k.s, {Rp: 1}, 600); }},
      {say: "The waves carry the energy away for good. [n]Engineers call this radiation damping.",
        n(k) { k.show(k.s.term); }},
    ],
  });

  // ---------------------------------------------------------------- 4. how SASSI does it
  // the mat in plan: 7 x 7 interaction nodes at 2 m (400 px = 12 m)
  scenes.push({
    id: "how", title: "How SASSI computes them",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "How SASSI does it", {x: 120, y: 70});
      s.links = K.path(g, "", {stroke: C.spring, "stroke-width": 1.6, fill: "none", hidden: true});
      s.p = plan(g, 160, 230, 400, 7, {hidden: true, r: 9});
      s.pLab = K.g(g, {hidden: true});
      K.text(s.pLab, 360, 735, "12 m mat · 7 × 7 = 49 interaction nodes", {cls: "t-label", anchor: "middle", size: 30, color: "var(--ink)"});
      const [a0, a1] = [s.p.pts[0], s.p.pts[1]];
      K.dim(s.pLab, a0[0], a0[1] + 46, a1[0], a1[1] + 46, "2 m");
      const c = s.p.pts[24];                                  // node 25, the centre
      s.fan = K.path(g, s.p.pts.filter((q, i) => i !== 24).map((q) => `M${c[0]},${c[1]}L${q[0]},${q[1]}`).join(""),
        {stroke: C.ssi, "stroke-width": 2, "stroke-opacity": 0.55, fill: "none", hidden: true});
      let all = "";
      for (let i = 0; i < 49; i++) for (let j = i + 1; j < 49; j++) all += `M${s.p.pts[i][0]},${s.p.pts[i][1]}L${s.p.pts[j][0]},${s.p.pts[j][1]}`;
      s.links.setAttribute("d", all);
      s.links.style.strokeOpacity = "0.22";
      s.p.g.insertBefore(s.links, s.p.g.children[1]);                    // above the mat, below the nodes
      s.rip = K.ripples(g, {cx: c[0], cy: c[1], r0: 14, r1: 290, n: 4, color: C.dash, squash: 1, width: 3, half: false});
      s.ripLab = note(g, 360, 785, "waves from one point load: schematic", {anchor: "middle", hidden: true, color: "var(--dash)"});
      s.push = K.circle(g, c[0], c[1], 15, {fill: C.ssi, stroke: C.ink, "stroke-width": 3, hidden: true, in: "pop"});
      s.list = K.bullets(s.root, [
        {t: "A push at one point moves every point"},
        {t: "Turned around: springs linking every point"},
        {t: "Summed up for the whole mat"},
      ], {x: 760, y: 250, w: 1060, num: true, size: 40});
      s.files = [["FOUNSTIF · springs", "var(--spring)"], ["FOUNDASH · dashpots", "var(--dash)"], ["FOUNDAMP · damping ratios", "var(--good)"]]
        .map(([t, col], i) => K.pill(s.root, t, {x: 812, y: 620 + i * 92, color: col, size: 34}));
      s.Rp = 0;
    },
    tick(s, t) {
      s.rip.amp = s.Rp;
      s.rip.set(t, 1.8);
    },
    beats: [
      {say: "How does SASSI compute all this? [p]It starts with the points where the building and the soil hold hands.",
        go(k) { k.show(k.s.head); }, p(k) { k.show(k.s.p.g); k.show(k.s.p.dots, {stagger: 12}); }},
      {say: "[n]These are the interaction nodes: here, 49 points on a 2 metre grid.",
        n(k) { k.show(k.s.pLab); k.pulse(k.s.p.dots, {amp: 0.4}); }},
      {say: "[k]SASSI works out how a push at any one point moves every other point, through the ground.",
        k(k) { k.show([k.s.push, k.s.fan, k.s.ripLab, k.s.list.items[0]]); k.tween(k.s, {Rp: 1}, 800); }},
      {say: "[x]Turn that around, and you get the soil's springs, linking every point with every other.",
        x(k) { k.hide([k.s.fan, k.s.ripLab]); k.tween(k.s, {Rp: 0}, 600); k.show([k.s.links, k.s.list.items[1]]); }},
      {say: "[g]ANALYS then sums them up for the whole mat, [f]in three files: springs, dashpots, and damping ratios.",
        g(k) { k.show(k.s.list.items[2]); k.dim(k.s.links, true); }, f(k) { k.show(k.s.files, {stagger: 150}); }},
    ],
  });

  // ---------------------------------------------------------------- 5. check against a formula
  scenes.push({
    id: "check", title: "Check it against a formula",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Check it against a formula", {x: 120, y: 70});
      s.side = K.g(g, {hidden: true});
      K.soil(s.side, {x: 140, y: 560, w: 620, layers: [{h: 200, kind: "sand"}], hs: {h: 120, kind: "sand"}, labels: false});
      K.rect(s.side, 330, 538, 240, 22, {cls: "mat"});
      K.dim(s.side, 330, 500, 570, 500, "");
      K.text(s.side, 450, 482, "12 m, rigid", {cls: "t-label", anchor: "middle", size: 30, color: "var(--ink)"});
      halo(K.text(s.side, 160, 840, "uniform soil, Vs 200 m/s", {cls: "t-label", size: 30, color: "var(--ink)"}));
      const P = K.plot(s.svg, {x: 1000, y: 230, w: 780, h: 420, xr: [0, 3], yr: [0, 25], yticks: [0, 10, 20], yfmt: (v) => "+" + v + " %",
        ylabel: "SASSI above the formula", ylabelOffset: 100});
      s.P = P;
      K.text(P.g, 1000, 200, "springs at 0.2 Hz (practically static) · formula: Pais & Kausel", {cls: "t-small", size: 24});
      // FOUNSTIF at 0.195 Hz against Pais & Kausel for a square (nu = 1/3): 5.52 GB, 7.05 GB, 6.00 GB^3
      s.bars = [["sliding", IJ(0, 0), 5.52 * G * B, C.good], ["vertical", IJ(2, 2), 7.05 * G * B, C.good], ["rocking", IJ(4, 4), 6.00 * G * B ** 3, C.ssi]]
        .map(([name, j, pk, col], i) => {
          const v = 100 * (ST[j][0] / pk - 1), bg = K.g(P.g, {hidden: true, in: "fade"}), x = P.X(i + 0.5);
          K.rect(bg, x - 60, P.Y(v), 120, P.Y(0) - P.Y(v), {fill: col, "fill-opacity": 0.85, rx: 4});
          K.text(bg, x, P.Y(v) - 16, "+" + v.toFixed(0) + " %", {cls: "t-label", anchor: "middle", size: 34, color: "var(--ink)"});
          K.text(P.g, x, P.Y(0) + 44, name, {cls: "t-label", anchor: "middle", size: 32, color: "var(--ink)"});
          return bg;
        });
      s.card = K.card(s.root, {x: 1000, y: 790, w: 780, kind: "check", title: "Check the grid like a mesh",
        body: "refine it, and watch the rocking spring"});
    },
    beats: [
      {say: "Can you trust these numbers? [c]Test them where a textbook formula applies: this rigid square mat on uniform soil.",
        go(k) { k.show(k.s.head); }, c(k) { k.show(k.s.side); }},
      {say: "[b]Sliding and vertical land within about 6 percent.",
        b(k) { k.show(k.s.P.g); k.show([k.s.bars[0], k.s.bars[1]], {stagger: 150, delay: 200}); }},
      {say: "[r]Rocking comes out stiffer: partly from the coarse grid of points, partly from the formula itself.",
        r(k) { k.show(k.s.bars[2]); }},
      {say: "[m]So check that grid like any finite-element mesh: refine it, and watch the rocking spring.",
        m(k) { k.show(k.s.card); }},
    ],
  });

  // ---------------------------------------------------------------- 6. springs that change with frequency
  scenes.push({
    id: "frequency", title: "Springs that change with frequency",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Springs that change with frequency", {x: 120, y: 70});
      const iF = FS.map((f, i) => i).filter((i) => FS[i] <= 10.1);
      const fs = iF.map((i) => FS[i]);
      const ratio = (j) => iF.map((i) => ST[j][i] / ST[j][0]);
      const P = K.plot(s.svg, {x: 220, y: 210, w: 780, h: 460, xr: [0, 10], yr: [0.6, 1.1], xticks: [0, 2, 4, 6, 8, 10], yticks: [0.6, 0.8, 1],
        xlabel: "frequency (Hz)", ylabel: "spring ÷ static spring", ylabelOffset: 84});
      s.P = P;
      s.lOne = P.line([0, 10], [1, 1], {color: C.ref, width: 4, draw: true});
      s.oneLab = halo(K.text(P.g, P.X(10) - 10, P.Y(1) - 16, "a formula: one number", {cls: "t-label", anchor: "end", size: 28, color: "var(--ink2)", hidden: true}));
      s.kZ = P.line(fs, ratio(IJ(2, 2)), {color: C.violet, width: 5, draw: true});
      s.kYY = P.line(fs, ratio(IJ(4, 4)), {color: C.ssi, width: 6, draw: true});
      s.leg = P.legend([{label: "vertical", color: C.violet}, {label: "rocking", color: C.ssi}], {x: 260, y: 560, size: 28, hidden: true});
      // the rocking spring at 8 Hz (FOUNSTIF)
      const i8 = FS.findIndex((f) => Math.abs(f - 8) < 0.1), r8 = ST[IJ(4, 4)][i8] / ST[IJ(4, 4)][0];
      s.pk = K.g(P.g, {hidden: true});
      K.circle(s.pk, P.X(FS[i8]), P.Y(r8), 10, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      halo(K.text(s.pk, P.X(FS[i8]) - 24, P.Y(r8) + 50, `−${Math.round(100 * (1 - r8))} % at ${FS[i8].toFixed(0)} Hz`, {cls: "t-big", size: 40, anchor: "end", color: "var(--ssi)"}));
      // a shear wave under the mat: wavelength to scale (20 px per metre, Vs 200 m/s)
      s.W = K.g(g, {hidden: true});
      K.rect(s.W, 1140, 210, 660, 460, {rx: 14, fill: "rgba(20,32,51,.6)", stroke: "#2a3b52", "stroke-width": 2});
      K.text(s.W, 1170, 256, "shear wavelength and the mat, to scale", {cls: "t-label", size: 28});
      K.rect(s.W, 1350, 360, 240, 20, {cls: "mat"});
      K.text(s.W, 1470, 345, "12 m", {cls: "t-label", anchor: "middle", size: 26, color: "var(--ink2)"});
      s.wave = K.path(s.W, "", {stroke: C.wave, "stroke-width": 5, fill: "none"});
      s.lamTxt = K.text(s.W, 1470, 590, "", {cls: "t-label", anchor: "middle", size: 30, color: "var(--wave)"});
      s.fTxt = K.text(s.W, 1470, 634, "", {cls: "t-label", anchor: "middle", size: 32, color: "var(--ink)"});
      s.f = 0.5; s.Wv = 0;
    },
    tick(s, t) {
      // lambda = Vs / f; the crests move at the same speed at every frequency (Vs, slowed down)
      const lam = 20 * 200 / s.f, xc = 160 * t * s.Wv, pts = [];
      for (let x = 1170; x <= 1770; x += 5) pts.push([x, 480 + 40 * Math.sin(TAU * (x - 1470 - xc) / lam)]);
      s.wave.setAttribute("d", K.d(pts));
      s.lamTxt.textContent = "λ = " + (200 / s.f).toFixed(0) + " m";
      s.fTxt.textContent = (s.f < 2 ? "slow shaking" : s.f > 12 ? "fast shaking" : "faster") + " · " + s.f.toFixed(1) + " Hz";
    },
    beats: [
      {say: "Now the big difference from a formula. [c]A car's springs stay the same, however fast you bounce.",
        go(k) { k.show(k.s.head); }, c(k) { k.show(k.s.P.g); k.draw(k.s.lOne, 800); k.show(k.s.oneLab, {delay: 400}); }},
      {say: "[k]The soil's springs don't: the faster the shaking, the softer they get.",
        k(k) { k.draw(k.s.kZ, 1000); k.draw(k.s.kYY, 1000, {delay: 200}); k.show(k.s.leg, {delay: 400}); }},
      {say: "Why? [w]What matters is the size of the mat, compared with the length of the waves.",
        w(k) { k.show(k.s.W); k.tween(k.s, {Wv: 1}, 600); }},
      {say: "[s]Slow shaking makes long waves; [f]fast shaking makes short ones, about the size of the mat.",
        s(k) { k.tween(k.s, {f: 0.5}, 10); }, f(k) { k.tween(k.s, {f: 20}, 2000); }},
      {say: "[r]By 8 hertz, the rocking spring has lost almost 30 percent.",
        r(k) { k.show(k.s.pk); k.pulse(k.s.pk, {amp: 0.08}); }},
    ],
  });

  // ---------------------------------------------------------------- 7. what it means: damping
  scenes.push({
    id: "damping", title: "How much damping the soil adds",
    build(s) {
      s.head = K.heading(s.root, "How much damping does the soil add?", {x: 120, y: 70});
      const iF = FS.map((f, i) => i).filter((i) => FS[i] <= 10.1);
      const fs = iF.map((i) => FS[i]);
      const P = K.plot(s.svg, {x: 230, y: 200, w: 1000, h: 520, xr: [0, 10], yr: [0, 0.8], xticks: [0, 2, 4, 6, 8, 10], yticks: [0, 0.2, 0.4, 0.6, 0.8],
        yfmt: (v) => Math.round(v * 100) + " %", xlabel: "frequency (Hz)", ylabel: "damping ratio", ylabelOffset: 96});
      s.P = P;
      s.dX = P.line(fs, iF.map((i) => DA[IJ(0, 0)][i]), {color: C.wave, width: 6, draw: true});
      s.dYY = P.line(fs, iF.map((i) => DA[IJ(4, 4)][i]), {color: C.ssi, width: 6, draw: true});
      const at = (f) => FS.findIndex((x) => Math.abs(x - f) < 0.1);
      // a computed point (FOUNDAMP) with its value as the label
      const call = (j, i, dx, dy, col, anchor) => {
        const cg = K.g(P.g, {hidden: true});
        const v = DA[IJ(j, j)][i];
        K.circle(cg, P.X(FS[i]), P.Y(v), 10, {fill: col, stroke: "#0a111d", "stroke-width": 3});
        halo(K.text(cg, P.X(FS[i]) + dx, P.Y(v) + dy, Math.round(100 * v) + " %", {cls: "t-big", size: 40, anchor: anchor || "middle", color: col}));
        return cg;
      };
      s.c15 = call(0, at(2), -14, -26, "var(--wave)", "end");
      s.c60 = call(0, at(8), -14, -26, "var(--wave)", "end");
      s.c3 = call(4, at(2), 0, 0, "var(--ssi)");
      { const t3 = s.c3.lastChild, v = DA[IJ(4, 4)][at(2)]; t3.setAttribute("x", P.X(3.3)); t3.setAttribute("y", P.Y(v + 0.09)); t3.setAttribute("text-anchor", "start");
        s.c3.insertBefore(K.line(null, P.X(FS[at(2)]) + 10, P.Y(v) - 8, P.X(3.25), P.Y(v + 0.09) + 4, {stroke: C.ssi, "stroke-width": 2}), t3); }
      s.legX = K.g(s.svg, {hidden: true});
      K.line(s.legX, 1330, 290, 1390, 290, {stroke: C.wave, "stroke-width": 7, "stroke-linecap": "round"});
      K.text(s.legX, 1410, 300, "sliding", {cls: "t-label", size: 34, color: "var(--ink)"});
      s.legY = K.g(s.svg, {hidden: true});
      K.line(s.legY, 1330, 360, 1390, 360, {stroke: C.ssi, "stroke-width": 7, "stroke-linecap": "round"});
      K.text(s.legY, 1410, 370, "rocking", {cls: "t-label", size: 34, color: "var(--ink)"});
      K.text(s.svg, 1330, 440, "radiation + 2 % material damping", {cls: "t-small", size: 24});
      s.typo = K.card(s.root, {x: 1330, y: 480, w: 470, kind: "warn", title: "Not a typo",
        body: "radiation damping can far exceed structural damping"});
      s.card = K.card(s.root, {x: 1330, y: 700, w: 470, title: "For your design",
        body: "the soil's help depends on how your building moves"});
      s.pill = K.pill(s.root, "→ floor spectra, forces, displacements", {x: 230, y: 870, color: "var(--ssi)", size: 30});
    },
    beats: [
      {say: "Now the shock absorbers. [p]How much damping does the soil give a foundation?",
        go(k) { k.show(k.s.head); }, p(k) { k.show(k.s.P.g); }},
      {say: "[x]In sliding, a lot: 15 percent at 2 hertz, [y]and 60 percent at 8.",
        x(k) { k.show(k.s.legX); k.draw(k.s.dX, 1000); k.show(k.s.c15, {delay: 400}); }, y(k) { k.show([k.s.c60, k.s.typo]); }},
      {say: "[r]In rocking, very little at low frequency: 3 percent at 2 hertz, mostly the soil's own damping.",
        r(k) { k.show(k.s.legY); k.draw(k.s.dYY, 1000); k.show(k.s.c3, {delay: 400}); }},
      {say: "[f]And that damping flows straight into your floor spectra, forces and displacements.",
        f(k) { k.show([k.s.card, k.s.pill], {stagger: 150}); }},
    ],
  });

  // ---------------------------------------------------------------- 8. recap
  scenes.push({
    id: "recap", title: "Recap",
    build(s) {
      s.head = K.heading(s.root, "Recap", {x: 120, y: 80, size: "h1"});
      s.list = K.bullets(s.root, [
        {t: "Soil acts like springs and shock absorbers", sub: "the foundation's impedance"},
        {t: "They change with frequency", sub: "SASSI computes them from the interaction nodes"},
        {t: "Radiation damping: energy leaving as waves", sub: "large in sliding, small in rocking"},
      ], {x: 120, y: 230, w: 940, num: true});
      s.q = K.card(s.root, {x: 1120, y: 230, w: 700, kind: "check", title: "Check yourself",
        body: "At 2 Hz, a tall building rocks; a squat one slides. Which gets more soil damping?"});
      s.a = K.card(s.root, {x: 1120, y: 560, w: 700, title: "Answer",
        body: "<b>The squat one</b>: about 15 % in sliding, against 3 % in rocking."});
      s.next = K.pill(s.root, "Next · Lesson 4: your first SSI analysis →", {x: 120, y: 880, color: "var(--wave)", size: 32});
    },
    beats: [
      {say: "Let's recap. [a]The soil under a foundation acts like springs and shock absorbers: its impedance.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[a]They change with frequency, so SASSI computes them at every frequency, from the interaction nodes.",
        a(k) { k.show(k.s.list.items[1]); }},
      {say: "[a]And radiation damping, energy leaving as waves, is large in sliding and small in rocking.",
        a(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Check yourself: at 2 hertz, a tall building rocks and a squat one slides. Which gets more soil damping?",
        q(k) { k.show(k.s.q); }, gap: 1500},
      {say: "[a]The squat one: about 15 percent damping in sliding, against only 3 percent in rocking.",
        a(k) { k.show(k.s.a); }},
      {say: "[n]Next: your first full SSI analysis, a stick on a surface mat. That's lesson four.",
        n(k) { k.show(k.s.next); }},
    ],
  });

  SV.video({
    id: "03", n: 3, part: "Fundamentals", lesson: "03-impedance",
    title: "Foundation impedance: soil springs and dashpots",
    subtitle: "Why the soil under a foundation acts like springs and shock absorbers that change with frequency, how SASSI computes them, and what radiation damping means for your design.",
    next: {href: "04.html", title: "Your first SSI analysis: a stick on a surface mat"},
    pron: [[/\bINT\b/g, "int"]],
    scenes,
  });
})();
