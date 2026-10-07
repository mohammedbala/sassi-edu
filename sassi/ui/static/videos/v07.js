/* Explainer video, lesson 7 (newcomer edition): Three earthquake directions and design ISRS -- three waves rising
 * straight up (SV, SH, P), one ANALYS run for the three directions, an off-centre mass that couples them, SRSS and
 * 100-40-40, then the envelope and the widened peaks of the design floor spectrum.
 * Numbers: sassi/ui/lessons/07_three_components.md; curves: d07.js (the lesson's run, python -m sassi.ui.video_data).
 * Geometry of example 5 to scale (ft): 40 ft x 40 ft basemat 5 ft thick, stick nodes at 13 ft and 26 ft, roof mass
 * at (6.5, 3.25, 26) on a rigid arm; site 16 ft sand, 39 ft gravel, rock.  The model moves with the computed transfer
 * functions of the X input at 5.0 Hz (amplitude and phase); the wave fields are the exact pulse of a uniform
 * half-space of the lesson's sand (PH.halfspacePulse); the widening animates PH.broaden of the SRSS and ends on
 * BROADEN's own result (roof_Y_design.rs). */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["07"], PH = K.PH;
  const TAU = 2 * Math.PI;
  const CX = C.wave, CY = C.violet, CZ = C.spring;           // colour of each direction: X, Y, Z
  const T = (k) => D["ex05/" + k];

  /** The largest y and where it is. */
  const argmax = (xs, ys) => { let i = 0; ys.forEach((y, j) => { if (y > ys[i]) i = j; }); return {x: xs[i], y: ys[i], i}; };
  const pk = (k) => argmax(T(k).f, T(k).sa);
  const g2 = (v) => v.toFixed(2);
  /** A transfer function at the SSI frequency nearest f, as a complex number [re, im] (phase in rad). */
  function tfAt(d, f) {
    let j = 0;
    d.f.forEach((x, i) => { if (Math.abs(x - f) < Math.abs(d.f[j] - f)) j = i; });
    return [d.amp[j] * Math.cos(d.ph[j]), d.amp[j] * Math.sin(d.ph[j])];
  }
  const FREQS = T("X_00030TR_X.TFU").f;                     // the 19 SSI frequencies
  const F5 = FREQS.reduce((a, f) => (Math.abs(f - 5) < Math.abs(a - 5) ? f : a));       // 5.005 Hz
  /** The X input at 5.0 Hz (per unit control motion): roof mass node 30 in X, Y, Z, the stick top's twist
   *  (node 27, rad per ft), the X motion of nodes 13 (basemat centre), 26 and 27.  The arm is rigid, so the
   *  stick top moves as u27 = u30 - theta x r, r = (6.5, 3.25, 0): in Y, u27y = u30y - 6.5 theta_z (lesson: 1.15). */
  const H = {x30: tfAt(T("X_00030TR_X.TFU"), F5), y30: tfAt(T("X_00030TR_Y.TFU"), F5), z30: tfAt(T("X_00030TR_Z.TFU"), F5),
    rz27: tfAt(T("X_00027R_ZZ.TFU"), F5), x13: tfAt(T("00013TR_X.TFU"), F5), x26: tfAt(T("00026TR_X.TFU"), F5), x27: tfAt(T("00027TR_X.TFU"), F5)};
  H.y27 = [H.y30[0] - 6.5 * H.rz27[0], H.y30[1] - 6.5 * H.rz27[1]];
  /** Re(h e^{i phi}). */
  const re = (h, ph) => h[0] * Math.cos(ph) - h[1] * Math.sin(ph);
  const FDRAW = 1;                                           // ft drawn per unit of control motion (as the lesson's animation)
  const SLOW = 10;                                           // 5.0 Hz shown slowed down 10 x

  /** Oblique projection (Y recedes up and to the right at half scale): model (x, y, z) ft -> screen px. */
  const proj = (cx, cy, S) => (x, y, z) => [cx + S * (x + 0.433 * y), cy - S * (z + 0.25 * y)];
  /** Path of a polygon of model points through P. */
  const poly = (P, pts) => "M" + pts.map((q) => P(...q).map((v) => v.toFixed(1)).join(",")).join("L") + "Z";
  /** A box x0..x1, y0..y1, z0..z1 seen in P: its right, front and top faces. */
  function box3(g, P, b, st) {
    const [x0, x1, y0, y1, z0, z1] = b;
    const out = {};
    out.right = K.path(g, poly(P, [[x1, y0, z0], [x1, y1, z0], [x1, y1, z1], [x1, y0, z1]]), Object.assign({}, st, {fill: st.side || st.fill}));
    out.front = K.path(g, poly(P, [[x0, y0, z0], [x1, y0, z0], [x1, y0, z1], [x0, y0, z1]]), st);
    out.top = K.path(g, poly(P, [[x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]]), Object.assign({}, st, {fill: st.topFill || st.fill}));
    return out;
  }
  const box = (root, x, y, w, h, html, color, o) => K.h("div", Object.assign({x, y, w, h, html, in: "pop", align: "center", size: 32,
    style: {border: `3px solid ${color}`, borderRadius: "16px", background: "rgba(15,26,43,.94)", display: "flex", alignItems: "center",
      justifyContent: "center", flexDirection: "column", fontWeight: "650", lineHeight: "1.25"}}, o || {}), root);
  const note = (g, x, y, str, o) => K.text(g, x, y, str, Object.assign({cls: "t-label", size: 24, color: "var(--muted)", hidden: true}, o));
  /** The stick of example 5 at rest in P: basemat (40 ft, 5 ft thick), stick to 26 ft, masses at 13 ft and (6.5, 3.25, 26). */
  function model(g, P, o) {
    K.defs();
    const m = {g: K.g(g, {hidden: o && o.hidden})};
    m.mat = box3(m.g, P, [-20, 20, -20, 20, -2.5, 2.5], {fill: "#2b3a4e", side: "#223044", topFill: "#3b4a5e", stroke: C.concrete, "stroke-width": 2});
    m.stick = K.path(m.g, K.d([P(0, 0, 2.5), P(0, 0, 13), P(0, 0, 26)]), {cls: "struct"});
    m.m1 = K.circle(m.g, ...P(0, 0, 13), o.r1 || 20, {cls: "mass", fill: "url(#sv-grad-mass)"});
    m.arm = K.path(m.g, `M${P(0, 0, 26)}L${P(6.5, 3.25, 26)}`, {stroke: C.ssi, "stroke-width": o.arm || 7, "stroke-linecap": "round"});
    m.top = K.circle(m.g, ...P(0, 0, 26), o.rt || 9, {fill: C.steel});
    m.m2 = K.circle(m.g, ...P(6.5, 3.25, 26), o.r2 || 24, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
    return m;
  }

  const scenes = [];

  // ---------------------------------------------------------------- 0. hook
  scenes.push({
    id: "hook", title: "Three directions at once",
    build(s) {
      K.defs();
      const g = K.g(s.svg);
      // the site of the lesson (16 ft sand, 39 ft gravel, rock) and example 5 on it, one scale: 6.77 px per ft
      const P = proj(580, 360, 22 / 3.25);
      const lay = [[-2.5, -18.5, "var(--sand)", "sv-pat-sand"], [-18.5, -57.5, "var(--gravel)", "sv-pat-gravel"], [-57.5, -70.5, "var(--rock)", "sv-pat-rock"]];
      lay.slice().reverse().forEach(([zt, zb, col, pat], j) => {
        const i = lay.length - 1 - j;                       // rock first, the sand (and the ground surface) last
        const st = {fill: col, "fill-opacity": i === 2 ? 0.5 : 0.78, side: col, stroke: "rgba(255,240,215,.3)", "stroke-width": 2};
        const b = box3(g, P, [-52, 52, -26, 26, zb, zt], st);
        if (i) b.top.remove();
        (i ? [b.right, b.front] : [b.right, b.front, b.top]).forEach((e) => K.path(g, e.getAttribute("d"), {fill: `url(#${pat})`}));
      });
      const fl = P(-52, -26, -2.5), fr = P(52, -26, -2.5), fb = P(52, -26, -57.5);
      s.wf = K.wavefronts(g, {x: fl[0] + 10, w: fr[0] - fl[0] - 20, y0: fb[1] - 6, y1: fl[1] + 14, lambda: 80, color: C.wave, width: 4, hidden: true});
      s.wfn = note(g, fl[0], P(0, -26, -70.5)[1] + 34, "waves rising: schematic · site and building to scale");
      s.mod = model(g, P, {r1: 11, r2: 15, rt: 6, arm: 5});
      const O = [150, 250];
      s.ax = K.arrow(g, O[0], O[1], O[0] + 150, O[1], {color: CX, width: 7, head: 24, hidden: true});
      s.axl = K.text(g, O[0] + 166, O[1] + 14, "X", {cls: "t-big", color: CX, hidden: true});
      s.ay = K.arrow(g, O[0], O[1], O[0] + 0.866 * 120, O[1] - 0.5 * 120, {color: CY, width: 7, head: 24, hidden: true});
      s.ayl = K.text(g, O[0] + 0.866 * 120 + 12, O[1] - 66, "Y", {cls: "t-big", color: CY, hidden: true});
      s.az = K.arrow(g, O[0], O[1], O[0], O[1] - 150, {color: CZ, width: 7, head: 24, hidden: true});
      s.azl = K.text(g, O[0] - 16, O[1] - 120, "Z", {cls: "t-big", anchor: "end", color: CZ, hidden: true});
      // what the equipment engineer wants: one floor spectrum per direction
      s.want = K.heading(s.root, "One floor spectrum per direction?", {x: 1170, y: 150, size: "h3"});
      s.minis = [["X", CX], ["Y", CY], ["Z", CZ]].map(([d, col], i) => {
        const y = 260 + i * 190;
        const grp = K.g(g, {hidden: true});
        K.text(grp, 1180, y + 82, d, {cls: "t-big", color: col});
        K.rect(grp, 1240, y, 440, 130, {rx: 12, fill: "rgba(20,32,51,.7)", stroke: "#2a3b52", "stroke-width": 2});
        K.path(grp, `M1262,${y + 16}L1262,${y + 112}L1660,${y + 112}`, {stroke: "#50657f", "stroke-width": 2, fill: "none"});
        K.text(grp, 1272, y + 36, "SA", {cls: "tick"});
        K.text(grp, 1652, y + 102, "f", {cls: "tick", anchor: "end"});
        K.text(grp, 1460, y + 90, "?", {cls: "t-big", size: 64, anchor: "middle", color: col, weight: 800});
        return grp;
      });
      // the fixed-base recipe
      s.recipe = ["1 · a response per direction", "2 · combine the directions", "3 · envelope the soil cases", "4 · widen the peaks"].map((t, i) =>
        box(s.root, 1170, 230 + i * 130, 640, 96, t, "#33507a", {in: "left"}));
      s.inp = K.pill(s.root, "SSI changes the input", {x: 1170, y: 770, size: 32, color: "var(--ssi)"});
      s.W = 0;
    },
    tick(s, t) {
      if (s.W) s.wf.set(t, 90);
    },
    beats: [
      {say: "[a]An earthquake doesn't shake the ground one way: [b]it pushes sideways in two directions, [c]and up and down.",
        a(k) { k.s.W = 1; k.show([k.s.wf, k.s.wfn]); }, b(k) { k.show([k.s.ax.g, k.s.axl, k.s.ay.g, k.s.ayl], {stagger: 150}); }, c(k) { k.show([k.s.az.g, k.s.azl]); }},
      {say: "[q]Your equipment engineers want one floor spectrum per direction. How do you get there from three shakings?",
        q(k) { k.show(k.s.want); k.show(k.s.minis, {stagger: 150, delay: 200}); }},
      {say: "[g]Good news: it's the recipe you know from fixed-base design. [i]SSI only changes the input.",
        g(k) { k.hide([k.s.want, k.s.minis]); k.show(k.s.recipe, {stagger: 150, delay: 200}); },
        i(k) { k.show(k.s.inp); k.pulse(k.s.wf, {amp: 0.03}); }},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 7",
    build(s) {
      s.t = K.titleCard(s.root, {n: 7, part: "Design applications", title: "Three earthquake directions and design ISRS", sub: "Three waves, one run, and the recipe that turns them into design floor spectra."});
      s.t.title.style.maxWidth = "1180px";
    },
    beats: [
      {say: "Lesson seven: three earthquake directions, and the design floor spectra.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 100}); k.show(t.title, {delay: 200}); k.show(t.sub, {delay: 350}); k.show(t.rule, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. three waves, straight up
  // a uniform half-space of the lesson's sand (Vs 1,000 ft/s, Vp 2,000 ft/s), a 4 Hz pulse, 130 ft shown at 2.92 px/ft;
  // particle motion in model axes (X, Y, Z), drawn in the oblique projection of the whole video
  const WAVES = [
    {name: "X · SV wave", V: 1000, d: [1, 0, 0], color: CX},
    {name: "Y · SH wave", V: 1000, d: [0, 1, 0], color: CY},
    {name: "Z · P wave", V: 2000, d: [0, 0, 1], color: CZ},
  ];
  const UNIT = 56;                                           // px per unit of surface motion (one factor for the scene)
  const scr = (d) => [d[0] + 0.433 * d[1], -(d[2] + 0.25 * d[1])];   // model direction -> screen direction
  scenes.push({
    id: "waves", title: "Three waves, straight up",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Three waves, rising straight up", {x: 120, y: 70, size: "h2"});
      const top = 360, pxm = 380 / 130, Hm = 130, Dm = 40, NX = 6, NZ = 13;     // px per ft, depth and width (ft)
      s.cols = WAVES.map((wv, i) => {
        const x0 = 170 + i * 580, w = 400, P = proj(x0, top, pxm);
        const cg = K.g(g, {hidden: true});
        box3(cg, (x, y, z) => P(x, y, z), [0, w / pxm, 0, Dm, -Hm, 0], {fill: "rgba(212,170,110,.28)", side: "rgba(212,170,110,.18)", topFill: "rgba(212,170,110,.4)", stroke: "rgba(255,240,215,.35)", "stroke-width": 2});
        K.text(cg, x0, top - 96, wv.name, {cls: "t-label", size: 38, color: wv.color, weight: 750});
        const sd = scr(wv.d), L = Math.hypot(sd[0], sd[1]), ax = x0 + w - 40, ay = top - 108;
        K.arrow(cg, ax - (34 * sd[0]) / L, ay - (34 * sd[1]) / L, ax + (34 * sd[0]) / L, ay + (34 * sd[1]) / L, {color: wv.color, width: 5, head: 16, both: true});
        K.arrow(cg, x0 - 30, top + Hm * pxm - 10, x0 - 30, top + Hm * pxm - 110, {color: C.muted, width: 4, head: 14});
        const dots = [];
        for (let kz = 0; kz < NZ; kz++) for (let kx = 0; kx < NX; kx++) {
          const zf = (kz + 0.5) / NZ;
          const x = x0 + (kx + 0.5) * w / NX, y = top + zf * Hm * pxm;
          dots.push({l: K.line(cg, x, y, x, y, {stroke: wv.color, "stroke-width": 3, opacity: 0.55}), c: K.circle(cg, x, y, 7, {fill: wv.color}), x, y, z: zf * Hm});
        }
        // the 40 ft x 40 ft basemat, 5 ft thick, on the surface
        const sl = K.g(g, {hidden: true});
        const c0 = w / pxm / 2;
        box3(sl, P, [c0 - 20, c0 + 20, 0, 40, 0, 5], {fill: "var(--concrete-fill)", side: "#5d6b7c", topFill: "#8796a8", stroke: "var(--concrete)", "stroke-width": 2});
        return {g: cg, dots, wv, slab: sl, sd};
      });
      s.dim = K.dim(g, 110, top, 110, top + Hm * pxm, "130 ft", {side: "left", hidden: true});
      s.note = note(g, 170, top + Hm * pxm + 52, `uniform sand: Vs 1,000 ft/s, Vp 2,000 ft/s · 4 Hz pulse · slowed 20 × · 1 unit of surface motion = ${UNIT} px`);
      s.fast = K.text(g, 1330, top + Hm * pxm + 100, "P: twice as fast here", {cls: "t-label", size: 30, color: "var(--spring)", hidden: true});
      s.unit = K.pill(s.root, "surface motion = 1 unit, in each direction", {x: 170, y: 900, size: 32, color: "var(--ink2)"});
      s.notw = K.pill(s.root, "the whole foundation at once: no twist from the input", {x: 170, y: 900, size: 32, color: "var(--good)"});
      s.W = [0, 0, 0];
    },
    tick(s, t) {
      const fp = 4, T0 = 0.2, tp = ((t * 0.05) % (2 * T0)) - T0;             // physical s; slowed 20 times
      s.cols.forEach((c, i) => {
        const a = s.W[i] * UNIT / 2;                                         // the surface doubles the pulse: peak 2 -> 1 unit
        for (const d of c.dots) {
          const u = a * PH.halfspacePulse(d.z, tp, c.wv.V, fp);
          const px = d.x + u * c.sd[0], py = d.y + u * c.sd[1];      // the dot, and a tail from its rest position
          d.c.setAttribute("cx", px); d.c.setAttribute("cy", py);
          d.l.setAttribute("x2", px); d.l.setAttribute("y2", py);
        }
        const u0 = a * PH.halfspacePulse(0, tp, c.wv.V, fp);
        c.slab.setAttribute("transform", `translate(${(u0 * c.sd[0]).toFixed(2)} ${(u0 * c.sd[1]).toFixed(2)})`);
      });
    },
    beats: [
      {say: "In SASSI, each direction is a wave, rising straight up from the rock to the surface.",
        go(k) { k.show(k.s.head); k.show(k.s.cols.map((c) => c.g), {stagger: 150}); k.show([k.s.note, k.s.dim], {delay: 400}); }},
      {say: "[x]For X, a shear wave that shakes the ground along x. Engineers call it SV.", x(k) { k.tween(k.s.W, {0: 1}, 800); }},
      {say: "[y]For Y, the same kind of wave, shaking along y: SH.", y(k) { k.tween(k.s.W, {1: 1}, 800); }},
      {say: "[z]For up and down, a push-pull wave, like a pulse along a spring: P. It's twice as fast here.",
        z(k) { k.tween(k.s.W, {2: 1}, 800); k.show(k.s.fast, {delay: 400}); }},
      {say: "[u]Each wave is scaled so the ground surface moves by exactly one unit, in its own direction.",
        u(k) { k.show(k.s.unit); }},
      {say: "Each wave reaches the whole foundation at the same moment, [t]so on its own it can't twist a building.",
        go(k) { k.hide(k.s.unit); k.show(k.s.cols.map((c) => c.slab), {stagger: 150}); }, t(k) { k.show(k.s.notw); }},
    ],
  });

  // ---------------------------------------------------------------- 3. one run, three directions
  scenes.push({
    id: "analys", title: "One run, three directions",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "One run, three directions", {x: 120, y: 70, size: "h2"});
      const cols = [CX, CY, CZ], dirs = ["X", "Y", "Z"];
      s.site = K.h("div", {x: 130, y: 172, html: "SITE, three times", size: 32, in: "fade", style: {fontFamily: "var(--mono)", fontWeight: "700", color: "var(--ink2)"}}, s.root);
      s.in = dirs.map((d, i) => box(s.root, 130, 230 + i * 120, 330, 90, `${d} free field`, cols[i]));
      s.mid = box(s.root, 640, 230, 480, 330, "<div style='font:800 46px var(--mono)'>ANALYS</div><div style='font-size:32px;color:var(--ink2);margin-top:14px;font-weight:500'>soil and building,<br>solved once<br>per frequency</div>", "#5f7da6");
      s.out = dirs.map((d, i) => box(s.root, 1310, 230 + i * 120, 400, 90, `response to ${d}`, cols[i]));
      s.aIn = dirs.map((d, i) => K.arrow(g, 470, 275 + i * 120, 628, 275 + i * 120, {color: cols[i], width: 5, head: 18, hidden: true}));
      s.aOut = dirs.map((d, i) => K.arrow(g, 1132, 275 + i * 120, 1298, 275 + i * 120, {color: cols[i], width: 5, head: 18, hidden: true}));
      // the frequency sweep: the 19 SSI frequencies of the lesson
      const FX = (f) => 260 + f * 85, FY = 690;
      s.FX = FX;
      s.fr = K.g(g, {hidden: true});
      K.line(s.fr, FX(0), FY, FX(15), FY, {stroke: "#50657f", "stroke-width": 2});
      [0, 5, 10, 15].forEach((f) => K.text(s.fr, FX(f), FY + 42, f + " Hz", {cls: "tick", anchor: "middle", size: 26}));
      s.fdots = FREQS.map((f) => K.circle(s.fr, FX(f), FY, 9, {fill: "#0a111d", stroke: C.ssi, "stroke-width": 3}));
      K.text(s.fr, FX(0), FY - 34, `${FREQS.length} SSI frequencies, ${FREQS[0].toFixed(1)} to ${Math.round(FREQS[FREQS.length - 1])} Hz`, {cls: "t-label", size: 28});
      s.cur = K.circle(s.fr, FX(FREQS[0]), FY, 14, {fill: C.ssi});
      // the real picture: one matrix per frequency, three right-hand sides
      s.eq = K.eq(s.root, String.raw`\mathbf{C}(f)\;\big[\,\htmlClass{t-x}{U_X}\;\;\htmlClass{t-y}{U_Y}\;\;\htmlClass{t-z}{U_Z}\,\big] = \big[\,\htmlClass{t-x}{Q_X}\;\;\htmlClass{t-y}{Q_Y}\;\;\htmlClass{t-z}{Q_Z}\,\big]`,
        {x: 130, y: 820, size: 46, display: false, box: true});
      s.eq.querySelectorAll(".t-x").forEach((e) => { e.style.color = "var(--wave)"; });
      s.eq.querySelectorAll(".t-y").forEach((e) => { e.style.color = "var(--violet)"; });
      s.eq.querySelectorAll(".t-z").forEach((e) => { e.style.color = "var(--spring)"; });
      s.ansys = K.card(s.root, {x: 1010, y: 810, w: 780, kind: "ansys", title: "In ANSYS terms", size: 32, body: "three load vectors, one factorised matrix"});
      s.F = 0;
    },
    tick(s, t) {
      if (!s.F) return;
      const i = Math.floor(t * 2.2) % FREQS.length;
      if (s.fi !== i) {
        s.fi = i;
        s.cur.setAttribute("cx", s.FX(FREQS[i]));
        s.fdots.forEach((d, j) => d.setAttribute("fill", j <= i ? "rgba(255,181,71,.55)" : "#0a111d"));
      }
      const ph = (t * 2.2) % 1;
      const op = String(0.45 + 0.55 * Math.max(0, Math.sin(Math.PI * ph)));
      [s.aIn, s.aOut].forEach((list) => list.forEach((a) => { a.line.style.opacity = op; a.head.style.opacity = op; }));
    },
    beats: [
      {say: "[a]SITE builds the three free fields: how the ground shakes with no building on it, one per wave.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.site); k.show(k.s.in, {stagger: 150}); }},
      {say: "[b]Then one ANALYS run solves all three directions together.",
        b(k) { k.show(k.s.aIn.map((a) => a.g), {stagger: 100}); k.show(k.s.mid, {delay: 300}); k.show(k.s.fr, {delay: 400}); k.at(400, () => { k.s.F = 1; }); }},
      {say: "[a]In ANSYS terms: three load vectors, one factorised matrix.", a(k) { k.show([k.s.eq, k.s.ansys], {stagger: 150}); }},
      {say: "[r]Out come three sets of transfer functions: how each point answers a unit shake in X, in Y, and in Z.",
        r(k) { k.show(k.s.aOut.map((a) => a.g), {stagger: 100}); k.show(k.s.out, {stagger: 150, delay: 200}); }},
    ],
  });

  // ---------------------------------------------------------------- 4. an off-centre mass couples the directions
  const ISO = proj(560, 800, 16);                            // 16 px per ft
  scenes.push({
    id: "model", title: "An off-centre mass",
    build(s) {
      K.defs();
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "An off-centre roof mass", {x: 120, y: 70, size: "h2"});
      s.mod = model(g, ISO, {hidden: true});
      s.dims = K.g(g, {hidden: true});
      const d0 = ISO(-20, -20, -2.5), d1 = ISO(20, -20, -2.5);
      K.dim(s.dims, d0[0], d0[1] + 44, d1[0], d1[1] + 44, "40 ft");
      const s0 = ISO(-5.2, 0, 0), s1 = ISO(-5.2, 0, 13), s2 = ISO(-5.2, 0, 26);
      K.dim(s.dims, s0[0], s0[1], s1[0], s1[1], "13 ft", {side: "left"});
      K.dim(s.dims, s1[0], s1[1], s2[0], s2[1], "13 ft", {side: "left"});
      s.mnote = note(g, 140, 196, `steady state, X input at ${F5.toFixed(1)} Hz · drawn for ±${FDRAW} ft of ground motion · slowed ${SLOW} ×`);
      const q = ISO(6.5, 3.25, 26);
      s.lab2 = K.g(g, {hidden: true});
      K.text(s.lab2, q[0] + 96, q[1] - 10, "heavy roof mass", {cls: "t-label", size: 32, color: "var(--ssi)", weight: 750});
      K.text(s.lab2, q[0] + 96, q[1] + 30, "6.5 ft and 3.25 ft off-centre", {cls: "t-label", size: 30, color: "var(--ssi)"});
      s.zArr = K.g(g, {hidden: true});
      K.arrow(s.zArr, q[0] + 60, q[1] - 120, q[0] + 60, q[1] - 40, {color: CZ, width: 5, head: 16, both: true});
      K.text(s.zArr, q[0] + 84, q[1] - 74, "up and down", {cls: "t-label", size: 30, color: CZ, weight: 750});
      // plan view: 10.8 px per ft
      s.plan = K.g(g, {hidden: true});
      const PC = [1440, 400], PS = 10.8;
      s.PC = PC; s.PS = PS;
      const PP = (x, y) => [PC[0] + PS * x, PC[1] - PS * y];
      s.PP = PP;
      K.text(s.plan, PC[0], PC[1] - 20 * PS - 30, "seen from above", {cls: "t-label", size: 30, anchor: "middle"});
      s.pMat = K.rect(s.plan, PC[0] - 20 * PS, PC[1] - 20 * PS, 40 * PS, 40 * PS, {fill: "rgba(59,74,94,.6)", stroke: C.concrete, "stroke-width": 2});
      const m0 = PP(6.5, 3.25), a0 = PP(0, 0);
      K.dim(s.plan, a0[0], a0[1] + 70, m0[0], a0[1] + 70, "6.5 ft", {color: C.ink2});
      K.dim(s.plan, m0[0] + 46, a0[1], m0[0] + 46, m0[1], "3.25 ft", {color: C.ink2});
      s.ell = K.path(s.plan, "", {stroke: "rgba(238,243,249,.4)", "stroke-width": 2.5, "stroke-dasharray": "7 7", fill: "none"});
      s.pTop = K.path(s.plan, "", {fill: "rgba(159,176,195,.25)", stroke: C.steel, "stroke-width": 4});
      s.pArm = K.path(s.plan, "", {stroke: C.ssi, "stroke-width": 6});
      s.pMass = K.circle(s.plan, 0, 0, 17, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 3});
      // the ground (control) motion, along X
      s.gy = PC[1] + 20 * PS + 50;
      K.line(s.plan, PC[0] - FDRAW * PS, s.gy, PC[0] + FDRAW * PS, s.gy, {stroke: "rgba(76,195,255,.45)", "stroke-width": 4, "stroke-linecap": "round"});
      s.gnd = K.circle(s.plan, PC[0], s.gy, 9, {fill: CX});
      K.text(s.plan, PC[0] - 20 * PS - 16, s.gy + 10, "ground", {cls: "t-label", size: 28, color: "var(--wave)", anchor: "end"});
      s.yLab = K.text(g, PC[0] - 20 * PS + 18, PC[1] - 20 * PS + 46, "moves in Y too", {cls: "t-label", size: 32, color: CY, weight: 800, hidden: true});
      s.twist = K.text(g, PC[0] - 20 * PS + 18, PC[1] + 20 * PS - 22, "the stick twists", {cls: "t-label", size: 32, color: "var(--ink)", weight: 800, hidden: true});
      s.coup = K.card(s.root, {x: 1100, y: 800, w: 700, kind: "check", title: "Coupling", size: 32, body: "real floors: the mass is rarely centred on the stiffness"});
      s.A = 0;
    },
    tick(s, t) {
      const ph = (TAU * F5 * t) / SLOW, k = FDRAW * s.A;
      // displacements (ft) of the X input at 5.0 Hz
      const x13 = k * re(H.x13, ph), x26 = k * re(H.x26, ph), x27 = k * re(H.x27, ph), y27 = k * re(H.y27, ph);
      const x30 = k * re(H.x30, ph), y30 = k * re(H.y30, ph), z30 = k * re(H.z30, ph), th = k * re(H.rz27, ph);
      // isometric view (the basemat drawn translating with node 13; node 26 in X only, the stick top's
      // vertical motion is axial and negligible)
      const m = s.mod;
      const n0 = ISO(x13, 0, 2.5), n1 = ISO(x26, 0, 13), n2 = ISO(x27, y27, 26), n3 = ISO(6.5 + x30, 3.25 + y30, 26 + z30);
      m.stick.setAttribute("d", K.d([n0, n1, n2], true));
      m.m1.setAttribute("cx", n1[0]); m.m1.setAttribute("cy", n1[1]);
      m.top.setAttribute("cx", n2[0]); m.top.setAttribute("cy", n2[1]);
      m.arm.setAttribute("d", `M${n2}L${n3}`);
      m.m2.setAttribute("cx", n3[0]); m.m2.setAttribute("cy", n3[1]);
      const dx = ISO(x13, 0, 0)[0] - ISO(0, 0, 0)[0];
      [m.mat.right, m.mat.front, m.mat.top].forEach((e) => e.setAttribute("transform", `translate(${dx.toFixed(2)} 0)`));
      // plan view
      const PP = s.PP, PS = s.PS;
      s.pMat.setAttribute("transform", `translate(${(PS * x13).toFixed(2)} 0)`);
      const c = Math.cos(th), sn = Math.sin(th), hw = 2.5;
      const sq = [[-hw, -hw], [hw, -hw], [hw, hw], [-hw, hw]].map(([a, b]) => PP(x27 + a * c - b * sn, y27 + a * sn + b * c));
      s.pTop.setAttribute("d", K.d(sq.concat([sq[0]])) + "Z");
      const tp = PP(x27, y27), mp = PP(6.5 + x30, 3.25 + y30);
      s.pArm.setAttribute("d", `M${tp}L${mp}`);
      s.pMass.setAttribute("cx", mp[0]); s.pMass.setAttribute("cy", mp[1]);
      if (s.ellA !== s.A) {
        s.ellA = s.A;
        const pts = [];
        for (let i = 0; i <= 48; i++) { const a = TAU * i / 48; pts.push(PP(6.5 + k * re(H.x30, a), 3.25 + k * re(H.y30, a))); }
        s.ell.setAttribute("d", K.d(pts) + "Z");
      }
      s.gnd.setAttribute("cx", PP(k * Math.cos(ph), 0)[0]);
    },
    beats: [
      {say: "[a]Now the example: a two-storey stick, on a 40-foot square foundation.",
        go(k) { k.show(k.s.head); }, a(k) { k.show([k.s.mod.g, k.s.dims]); }},
      {say: "[m]But its heavy roof mass sits off-centre: 6.5 feet one way, 3.25 feet the other.",
        m(k) { k.show([k.s.lab2, k.s.plan]); k.pulse(k.s.mod.m2, {amp: 0.3}); }},
      {say: "[p]Shake only in X, and the roof mass moves in X, [y]but also sideways in Y, [r]as the stick twists.",
        p(k) { k.show(k.s.mnote); k.tween(k.s, {A: 1}, 800); }, y(k) { k.show(k.s.yLab); }, r(k) { k.show(k.s.twist); }},
      {say: "[z]It even bobs up and down, because the rocking stick lifts the off-centre mass.",
        z(k) { k.show(k.s.zArr); }},
      {say: "[c]Engineers call this coupling. Real floors do it too: their mass is rarely centred on their stiffness.",
        c(k) { k.show(k.s.coup); }},
    ],
  });

  // ---------------------------------------------------------------- 5. one input, three responses
  const PX = pk("X_00030TR_X01.RS"), PY = pk("X_00030TR_Y01.RS"), PZ = pk("X_00030TR_Z01.RS");
  scenes.push({
    id: "x-input", title: "One input, three responses",
    build(s) {
      s.head = K.heading(s.root, "One input, three responses", {x: 120, y: 70, size: "h2"});
      s.p = K.plot(s.svg, {x: 230, y: 200, w: 980, h: 520, xr: [0.1, 100], xlog: true, yr: [0, 4.5], xticks: [0.1, 1, 10, 100], yticks: [0, 1, 2, 3, 4],
        xlabel: "frequency (Hz)", ylabel: "roof mass: SA (g), 5 %", ylabelOffset: 74});
      s.tag = K.pill(s.root, "the X shaking only", {x: 260, y: 220, size: 30, color: "var(--wave)"});
      const L = (k, col) => s.p.line(T(k).f, T(k).sa, {color: col, width: 6, draw: true});
      s.sx = L("X_00030TR_X01.RS", CX); s.sy = L("X_00030TR_Y01.RS", CY); s.sz = L("X_00030TR_Z01.RS", CZ);
      s.dots = [[PX, CX], [PY, CY], [PZ, CZ]].map(([q, col]) => s.p.dot(q.x, q.y, {color: col, r: 10, hidden: true}));
      s.st = [[PX, "along X", "var(--wave)"], [PY, "along Y", "var(--violet)"], [PZ, "vertical, Z", "var(--spring)"]].map(([q, l, c], i) =>
        K.stat(s.root, {x: 1330, y: 190 + i * 190, w: 470, value: g2(q.y), unit: "g", label: `${l}, at ${q.x.toFixed(1)} Hz`, color: c, vsize: 96}));
      s.warn = K.card(s.root, {x: 230, y: 840, w: 1570, kind: "warn", title: "Why it matters", size: 32,
        body: `Keep only the X response, and you miss ${g2(PY.y)} g in Y and ${g2(PZ.y)} g vertically.`});
    },
    beats: [
      {say: "[s]Here are the floor spectra at the roof mass, for the X shaking alone.",
        go(k) { k.show(k.s.head); }, s(k) { k.show([k.s.p.g, k.s.tag]); }},
      {say: "[a]4.21 g along X, as expected. [b]But also 1.07 g along Y, [c]and 1.45 g vertically.",
        a(k) { k.draw(k.s.sx, 900); k.show([k.s.dots[0], k.s.st[0].el], {delay: 300}); },
        b(k) { k.draw(k.s.sy, 900); k.show([k.s.dots[1], k.s.st[1].el], {delay: 300}); },
        c(k) { k.draw(k.s.sz, 900); k.show([k.s.dots[2], k.s.st[2].el], {delay: 300}); }},
      {say: "[w]Keep only the shaking direction, and your equipment design misses both of those.", w(k) { k.show(k.s.warn); k.pulse([k.s.st[1].el, k.s.st[2].el], {amp: 0.06}); }},
    ],
  });

  // ---------------------------------------------------------------- 6. combine the three directions
  const QY = pk("Y_00030TR_Y01.RS"), QX = pk("X_00030TR_Y01.RS"), QZ = pk("Z_00030TR_Y01.RS");
  const QS = pk("roof_Y_SRSS.rs"), QL = pk("roof_Y_1004040.rs");
  scenes.push({
    id: "combine", title: "Combine the three",
    build(s) {
      s.head = K.heading(s.root, "Combine the three directions", {x: 120, y: 70, size: "h2"});
      s.p = K.plot(s.svg, {x: 230, y: 190, w: 940, h: 540, xr: [0.1, 100], xlog: true, yr: [0, 5], xticks: [0.1, 1, 10, 100], yticks: [0, 1, 2, 3, 4, 5],
        xlabel: "frequency (Hz)", ylabel: "roof mass, Y: SA (g), 5 %", ylabelOffset: 74});
      const L = (k, o) => s.p.line(T(k).f, T(k).sa, Object.assign({width: 5, draw: true}, o));
      s.cy = L("Y_00030TR_Y01.RS", {color: CY});
      s.cx = L("X_00030TR_Y01.RS", {color: CX});
      s.cz = L("Z_00030TR_Y01.RS", {color: CZ});
      s.lin = L("roof_Y_1004040.rs", {color: C.ink, width: 4});
      s.srss = L("roof_Y_SRSS.rs", {color: C.ssi, width: 8});
      s.lg = [[`from the Y shaking: ${g2(QY.y)} g`, CY], [`from X: ${g2(QX.y)} g`, CX], [`from Z: ${g2(QZ.y)} g`, CZ]].map(([l, c], i) => s.p.legend([{label: l, color: c}], {x: 260, y: 230 + i * 42, size: 28, hidden: true}));
      s.lgS = s.p.legend([{label: "SRSS", color: C.ssi}], {x: 260, y: 380, size: 28, hidden: true});
      s.lgL = s.p.legend([{label: "100-40-40", color: C.ink}], {x: 260, y: 422, size: 28, hidden: true});
      s.dS = s.p.dot(QS.x, QS.y, {color: C.ssi, r: 11, hidden: true});
      s.dL = s.p.dot(QL.x, QL.y, {color: C.ink, r: 10, hidden: true});
      s.eq = K.eq(s.root, String.raw`\text{SRSS}=\sqrt{R_X^2+R_Y^2+R_Z^2}`, {x: 1250, y: 200, size: 44});
      s.rule = K.card(s.root, {x: 1250, y: 370, w: 570, title: "100-40-40", size: 32, body: "each one in full, plus 40 % of the other two: keep the worst"});
      s.s1 = K.stat(s.root, {x: 1250, y: 610, w: 270, value: g2(QS.y), unit: "g", label: "SRSS", color: "var(--ssi)", vsize: 76});
      s.s2 = K.stat(s.root, {x: 1550, y: 610, w: 270, value: g2(QL.y), unit: "g", label: `100-40-40 · +${Math.round((QL.y / QS.y - 1) * 100)} %`, color: "var(--ink)", vsize: 76});
    },
    beats: [
      {say: "[a]Now turn it round: the roof's Y spectrum collects something from each of the three inputs.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.p.g); }},
      {say: "[y]4.07 g from the Y shaking, [x]1.07 from X, [z]and 0.34 from the vertical.",
        y(k) { k.draw(k.s.cy, 900); k.show(k.s.lg[0]); }, x(k) { k.draw(k.s.cx, 900); k.show(k.s.lg[1]); }, z(k) { k.draw(k.s.cz, 900); k.show(k.s.lg[2]); }},
      {say: "[e]The SRSS rule: square each one, add them up, take the square root. [s]Here, 4.18 g.",
        e(k) { k.show(k.s.eq); }, s(k) { k.draw(k.s.srss, 1000); k.show(k.s.lgS); k.show([k.s.dS, k.s.s1.el], {delay: 400}); }},
      {say: "[l]The 100-40-40 rule: each one in full, plus 40 percent of the other two, and keep the worst.",
        l(k) { k.show(k.s.rule); }},
      {say: "[m]Here that gives 4.57 g: 9 percent above the SRSS.",
        m(k) { k.draw(k.s.lin, 1000); k.show(k.s.lgL); k.show([k.s.dL, k.s.s2.el], {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 7. widen the peaks: the design ISRS
  const SR = T("roof_Y_SRSS.rs"), DS = T("roof_Y_design.rs");
  const F0 = argmax(SR.f, SR.sa).x;                          // 5.01 Hz
  // the plateau of BROADEN's result: where the design ISRS holds the SRSS peak
  const PLAT = (() => {
    const top = Math.max(...DS.sa), on = DS.f.filter((f, i) => DS.sa[i] >= top - 1e-9);
    return {lo: Math.min(...on), hi: Math.max(...on), top};
  })();
  scenes.push({
    id: "broaden", title: "Widen the peaks",
    build(s) {
      s.head = K.heading(s.root, "The design ISRS: widen the peaks", {x: 120, y: 70, size: "h2"});
      s.p = K.plot(s.svg, {x: 230, y: 190, w: 980, h: 540, xr: [2, 15], xlog: true, yr: [0, 5], xticks: [2, 3, 4, 5, 6, 8, 10, 15], yticks: [0, 1, 2, 3, 4, 5],
        xlabel: "frequency (Hz)", ylabel: "roof mass, Y: SA (g), 5 %", ylabelOffset: 74});
      s.band = K.rect(s.p.data, s.p.X(PLAT.lo), s.p.Y(5), s.p.X(PLAT.hi) - s.p.X(PLAT.lo), s.p.Y(0) - s.p.Y(5), {fill: "rgba(255,181,71,.10)", hidden: true});
      s.sr = s.p.line(SR.f, SR.sa, {color: "rgba(255,181,71,.5)", width: 4, draw: true});
      s.br = K.path(s.p.data, "", {cls: "curve", style: {stroke: C.ssi, strokeWidth: 8}});
      s.br.classList.add("sv-in", "sv-fade");
      s.dsD = K.d(DS.f.map((f, i) => [s.p.X(f), s.p.Y(DS.sa[i])]));
      s.e1 = s.p.vline(F0 * 0.85, {color: C.ink, width: 2, hidden: true});
      s.e2 = s.p.vline(F0 * 1.15, {color: C.ink, width: 2, hidden: true});
      s.lab = K.g(s.p.g, {hidden: true});
      K.text(s.lab, s.p.X(PLAT.lo) - 12, s.p.Y(0.4), `${g2(PLAT.lo)} Hz`, {cls: "t-label", size: 30, anchor: "end", color: "var(--ink)", weight: 750});
      K.text(s.lab, s.p.X(PLAT.hi) + 12, s.p.Y(0.4), `${g2(PLAT.hi)} Hz`, {cls: "t-label", size: 30, color: "var(--ink)", weight: 750});
      s.pk = s.p.text((PLAT.lo + PLAT.hi) / 2, PLAT.top, `${g2(PLAT.top)} g`, {cls: "t-label", size: 34, anchor: "middle", color: "var(--ssi)", weight: 800, dy: -24, hidden: true});
      s.leg = s.p.legend([{label: "SRSS", color: "rgba(255,181,71,.5)"}, {label: "design ISRS, ±15 %", color: C.ssi}], {x: 250, y: 230, size: 28, hidden: true});
      s.steps = K.bullets(s.root, [{t: "Envelope the soil cases", sub: "as in lesson 6"}, {t: "Widen every peak", sub: "uncertain building frequencies"}], {x: 1320, y: 200, w: 500, num: true, size: 36});
      s.code = K.code(s.root, ["BROADEN,9,0,15,4"], {x: 1320, y: 500, w: 500, size: 38});
      s.expl = K.h("div", {x: 1320, y: 660, w: 500, html: "envelope, then ±15 %", size: 32, in: "left", style: {color: "var(--ssi)", fontWeight: "650"}}, s.root);
      s.isrs = K.card(s.root, {x: 230, y: 840, w: 1590, title: "The design ISRS", size: 32, body: "per floor and direction: what your equipment engineers qualify against"});
      s.b = 0;
    },
    tick(s) {
      if (s.bLast === s.b) return;
      s.bLast = s.b;
      if (s.b >= 0.15 - 1e-9) { s.br.setAttribute("d", s.dsD); return; }           // BROADEN's own result
      const B = PH.broaden(SR.f, SR.sa, s.b);
      s.br.setAttribute("d", K.d(B.x.map((f, i) => [s.p.X(f), s.p.Y(B.y[i])])));
    },
    beats: [
      {say: "[a]Last step. With several soil cases, first keep their envelope, as in lesson six.",
        go(k) { k.show([k.s.head, k.s.p.g]); k.draw(k.s.sr, 1000); }, a(k) { k.show(k.s.steps.items[0]); }},
      {say: "[b]Then widen every peak, because your building's own frequencies are uncertain: concrete, mass, modelling.",
        b(k) { k.show(k.s.steps.items[1]); k.pulse(k.s.sr, {amp: 0.02}); }},
      {say: "[c]Your design basis sets how much. Here it's 15 percent either side of every peak.",
        c(k) { k.show([k.s.e1, k.s.e2]); }},
      {say: "[p]So the peak at 5 hertz becomes a plateau, from 4.26 to 5.76 hertz.",
        p(k) { k.show([k.s.br, k.s.leg]); k.tween(k.s, {b: 0.15}, 1500, {done: () => k.show(k.s.lab)}); }, hold: 600},
      {say: "[c]Equipment tuned anywhere in that band is designed for the full peak: 4.18 g.",
        c(k) { k.show([k.s.band, k.s.pk]); }},
      {say: "[k]In SASSI, one BROADEN line does both: the envelope of its sources, then the widening.",
        k(k) { k.show(k.s.code.el); k.type(k.s.code.lines, {cps: 22}); k.show(k.s.expl, {delay: 400}); k.show(k.s.isrs, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 8. recap
  scenes.push({
    id: "recap", title: "Recap",
    build(s) {
      s.head = K.heading(s.root, "Recap", {x: 120, y: 90, size: "h1"});
      s.list = K.bullets(s.root, [
        {t: "Three waves rise straight up", sub: "SV for X, SH for Y, P for Z"},
        {t: "One ANALYS run, three directions", sub: "an off-centre mass mixes them"},
        {t: "Combine, envelope, widen", sub: "SRSS or 100-40-40, then the design ISRS"},
      ], {x: 120, y: 240, w: 1000, num: true, size: 40});
      s.q = K.card(s.root, {x: 1180, y: 240, w: 640, kind: "check", title: "Check yourself", size: 34,
        body: "A symmetric building, shaken by these waves: does it twist?"});
      s.a = K.card(s.root, {x: 1180, y: 540, w: 640, title: "Answer", size: 34,
        body: "No: each wave reaches the whole foundation at once. Twist needs off-centre mass or stiffness."});
      s.next = K.pill(s.root, "Next · Lesson 8: from SSI results to design: forces, displacements, ANSYS →", {x: 120, y: 900, size: 30, color: "var(--wave)"});
    },
    beats: [
      {say: "[a]Let's recap. One: SASSI shakes in three directions, with three waves rising straight up.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[a]Two: one ANALYS run solves all three, and an off-centre mass mixes the directions.", a(k) { k.show(k.s.list.items[1]); }},
      {say: "[a]Three: combine them with SRSS or 100-40-40, envelope the soil cases, then widen the peaks.", a(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Check yourself. A perfectly symmetric building is shaken by these waves. Does it twist?",
        q(k) { k.show(k.s.q); }, gap: 1500},
      {say: "[a]No. Each wave reaches the whole foundation at once, so twisting comes only from off-centre mass or stiffness.",
        a(k) { k.show(k.s.a); }},
      {say: "[n]Next: from SSI results to design. Forces, displacements, and the way back to ANSYS.", n(k) { k.show(k.s.next); }},
    ],
  });

  SV.video({
    id: "07", n: 7, part: "Design applications", lesson: "07-three-components",
    title: "Three earthquake directions and design ISRS",
    subtitle: "Three waves rising straight up, one ANALYS run for all three, an off-centre mass that mixes them, and the recipe for the design floor spectrum.",
    next: {href: "08.html", title: "From SSI results to design: forces, displacements, ANSYS"},
    pron: [[/\b100-40-40\b/g, "100, 40, 40"]],
    scenes,
  });
})();
