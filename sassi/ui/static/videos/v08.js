/* Explainer video, lesson 8 (third edition: pace and accuracy): From SSI results to design -- forces,
 * displacements, ANSYS.
 * One big idea: floor spectra are not the only output.  The design also needs member forces (the soil is
 * inside them; check them by equilibrium), relative displacements (pick the reference for the question),
 * and the SSI motion can be handed to a detailed ANSYS model (Option A, LOADGEN) as long as that model does
 * not change how the foundation moves.
 *
 * Accuracy.  Units: ft, kip, s.  The stick of example 1 is drawn to scale (64 ft x 64 ft mat, 5 ft thick, four
 * 16 ft storeys) on the site of lesson 2 (16 ft sand over 39 ft gravel, cut with a fade).  It moves with the
 * lesson's computed record (d08.js), relative to the free field: the mat centre with 00041TR_X.THD, the floors
 * relative to the mat centre with MAT_000kkTR_X.THD, and the mat tilts with the vertical motion of its edges,
 * which the transfer functions give as ROCK = (H37z - H45z) / (H85x - H41x) = 0.266 times the roof
 * displacement relative to the mat centre, real and constant within 2 % from 0.5 to 6 Hz (0.260 ... 0.265,
 * 0.264 at 3.49 Hz; at t* the D commands of ex01_LGS.inp give 2 x 0.0102 ft against 0.0766 ft).  One factor
 * for every displacement (x 250), the record slowed down 4 times.  Inertia forces: 2,200 kips x the floor's
 * absolute acceleration (000kkTR_X.ACC).  The LOADGEN loads at t* = 1.855 s are the F and D commands of
 * ex01_LGS.inp (the lesson's run).  The "one rule" scene is an illustration computed with the
 * sway-rocking model of lesson 4 (SV.K.PH.ssiModel / ssiResponse).
 * Numbers: sassi/ui/lessons/08_design_outputs.md; curves: d08.js (python -m sassi.ui.video_data). */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["08"], PH = K.PH;
  const TAU = 2 * Math.PI;

  // ---------------------------------------------------------------- data helpers
  /** Largest absolute value of a history {t, v}: {t, v}. */
  function peak(h) {
    let i = 0;
    h.v.forEach((v, j) => { if (Math.abs(v) > Math.abs(h.v[i])) i = j; });
    return {t: h.t[i], v: h.v[i]};
  }
  /** Value of a (thinned, min/max) history at time t, eased between samples (smooth motion). */
  function sample(h, t) {
    const T = h.t, V = h.v;
    if (t <= T[0]) return V[0] * Math.max(0, t) / T[0];
    if (t >= T[T.length - 1]) return V[V.length - 1];
    let a = 0, b = T.length - 1;
    while (b - a > 1) { const m = (a + b) >> 1; if (T[m] <= t) a = m; else b = m; }
    const u = (t - T[a]) / (T[b] - T[a]), e = (1 - Math.cos(Math.PI * u)) / 2;
    return V[a] + (V[b] - V[a]) * e;
  }
  /** The part of a history up to t1, scaled: [ts, vs]. */
  function upto(h, t1, sc) {
    const xs = [], ys = [];
    h.t.forEach((t, i) => { if (t <= t1) { xs.push(t); ys.push(h.v[i] * (sc || 1)); } });
    return [xs, ys];
  }

  // ---------------------------------------------------------------- the computed motion of example 1
  const T_STAR = 1.855;                  // s: LOADGEN's critical time, the largest base shear (lesson)
  const EXAG = 250, SLOW = 4;            // every displacement x 250, the record slowed down 4 times
  const ROCK = 0.266;                    // mat-edge rocking / roof displacement relative to the mat centre (header)
  const W_FLOOR = 2200;                  // kips: floor weight; inertia force = W x a (g)
  const PXK = 0.04;                      // px per kip of the inertia-force arrows
  const T0 = 0.4, T1 = 20;               // the part of the record that plays (s)
  const HM = [82, 83, 84, 85].map((n) => D[`ex01/MAT_000${n}TR_X.THD`]);   // floors relative to the mat centre (ft)
  const H41 = D["ex01/00041TR_X.THD"];                                     // mat centre relative to the free field (ft)
  const HA = [82, 83, 84, 85].map((n) => D[`ex01/000${n}TR_X.ACC`]);      // floor absolute accelerations (g)
  const SHEAR = [6764, 6137, 4847, 2808];           // lesson table: peak storey shear FYI, base to top (kips), at t ≈ 1.85 s
  const FLOORF = SHEAR.map((v, i) => v - (SHEAR[i + 1] || 0));             // floor inertia forces = differences (kips)
  // LOADGEN at t* (ex01_LGS.inp): F at the floors 82-85 (kips, all -X, sum 6,772); D of the interface nodes
  // 37-45 on the x axis (x = -32 ... 32 ft): UX -0.00675 ... -0.00680 ft, UZ below (ft)
  const LG_F = [627, 1295, 2048, 2802];
  const LG_UX = -6.78e-3;
  const LG_UZ = [-10.181, -8.909, -7.327, -4.231, 0, 4.231, 7.327, 8.909, 10.181].map((v) => v * 1e-3);

  /** The state at record time t: u0 mat centre / free field, rel floors / mat centre, th mat tilt (rad,
   *  + = the top moves +X), F inertia forces of the floors (kips, + = +X). */
  function motion(t) {
    const rel = HM.map((h) => sample(h, t));
    return {u0: sample(H41, t), rel, th: ROCK * rel[3] / 64, F: HA.map((h) => -W_FLOOR * sample(h, t))};
  }
  /** Example 1's stick to scale (S px per ft) with its mat top at (x, y); glue: the soil surface under the
   *  mat follows it (a surface mat stays welded to the soil). */
  function exStick(g, x, y, S, glue) {
    const zs = [1, 2, 3, 4].map((i) => 16 * i * S), r = Math.round(Math.max(14, 3.5 * S));
    const gl = glue ? [K.path(g, "", {fill: "var(--sand)", "fill-opacity": 0.78}), K.path(g, "", {fill: "url(#sv-pat-sand)"})] : null;
    const st = K.stick(g, {x, y, z: zs, matW: 64 * S, matH: 5 * S, r, slabW: 16 * S});
    return {st, zs, S, matH: 5 * S, x, y, r, ff: 1, gl};
  }
  /** Pose the stick at the motion m (one factor EXAG for every displacement); ex.ff 1: the mat translation
   *  is shown (relative to the free field), 0: it is removed (relative to the mat centre). */
  function pose(ex, m) {
    const k = ex.S * EXAG, r = m.th * EXAG, c = Math.cos(r), sn = Math.sin(r);
    ex.st.set({sway: m.u0 * k * ex.ff, rock: r, bend: m.rel.map((u, i) => (u * k - (ex.zs[i] + ex.matH / 2) * sn) / c)});
    glue(ex);
    return ex.st.nodes();
  }
  /** The soil surface under the mat follows its underside (where the mat rises; it covers the rest). */
  function glue(ex) {
    if (!ex.gl) return;
    const st = ex.st.state, c = Math.cos(st.rock), sn = Math.sin(st.rock);
    const cx = ex.x + st.sway, cy = ex.y + ex.matH / 2, w = 32 * ex.S, h = ex.matH / 2, ys = ex.y + ex.matH;
    const bl = [cx - w * c - h * sn, Math.min(ys, cy - w * sn + h * c)], br = [cx + w * c - h * sn, Math.min(ys, cy + w * sn + h * c)];
    const d = `M${ex.x - w - 40},${ys}L${bl[0].toFixed(1)},${bl[1].toFixed(1)}L${br[0].toFixed(1)},${br[1].toFixed(1)}L${ex.x + w + 40},${ys}Z`;
    ex.gl.forEach((p) => p.setAttribute("d", d));
  }
  /** The lesson-2 site under a mat: 16 ft sand, 39 ft gravel, rock (S px per ft); cut at yMax with a fade. */
  function site(g, x, y, w, S, yMax) {
    const sand = 16 * S, grav = 39 * S;
    if (y + sand + grav + 40 <= yMax) return K.soil(g, {x, y, w, layers: [{h: sand, kind: "sand"}, {h: grav, kind: "gravel"}], hs: {h: yMax - y - sand - grav, kind: "rock"}, labels: false});
    return K.soil(g, {x, y, w, layers: [{h: sand, kind: "sand"}], hs: {h: yMax - y - sand, kind: "gravel"}, labels: false});
  }
  /** Inertia-force arrows at the floors (F in kips, px per kip). */
  function setForces(arrows, nodes, F, pxK, r) {
    arrows.forEach((a, i) => {
      const q = nodes[i + 1], f = F[i] * pxK, sg = Math.sign(f) || 1, x0 = q[0] + sg * (r + 4);
      a.update(x0, q[1], x0 + f, q[1]);
      a.g.style.opacity = String(Math.min(1, Math.abs(f) / 18));
    });
  }
  /** Advance the record clock of a scene (s.run 1: playing, slowed down SLOW times). */
  function clock(s, dt) {
    if (s.run) { s.tt += dt / SLOW; if (s.tt > T1) s.tt = T0; }
    if (s.clk) s.clk.textContent = !s.run && Math.abs(s.tt - T_STAR) < 1e-6 ? `t = ${T_STAR} s: the largest base shear` : `t = ${s.tt.toFixed(2)} s`;
  }
  const LAB = "example 1, computed: displacements × 250, slowed down 4 ×";

  /** A diamond bullet + text line (HTML). */
  const item = (root, x, y, w, color, html, size) => K.h("div", {x, y, w, size: size || 44, in: "left", color: "var(--ink)",
    html: `<span style="display:inline-block;width:22px;height:22px;border-radius:4px;background:${color};margin-right:30px;transform:rotate(45deg)"></span>${html}`}, root);
  /** A big number with its unit and a label to the right (HTML row). */
  const row = (root, x, y, val, unit, label, color) => K.h("div", {x, y, w: 920, in: "left", style: {display: "flex", alignItems: "center", gap: "34px"},
    html: `<span style="flex:0 0 330px;text-align:right;font-weight:800;font-size:104px;line-height:1;color:${color};font-variant-numeric:tabular-nums">${val}` +
      `<span style="font-size:44px;font-weight:600;color:var(--ink2);margin-left:10px">${unit}</span></span>` +
      `<span style="font-size:34px;line-height:1.25;color:var(--ink2)">${label}</span>`}, root);

  // ================================================================== scenes
  const scenes = [];

  // ---------------------------------------------------------------- 0. hook
  scenes.push({
    id: "hook", title: "Beyond the floor spectra",
    build(s) {
      const g = K.g(s.svg), S = 18 / 3.2;                         // px per ft
      site(g, 110, 600 + 5 * S, 660, S, 990);
      s.ex = exStick(g, 440, 600, S, true);
      s.forces = K.g(g, {hidden: true});
      s.farr = [0, 1, 2, 3].map(() => K.arrow(s.forces, 0, 0, 1, 0, {color: C.ssi, width: 7, head: 20}));
      s.disp = K.arrow(g, 440, 200, 480, 200, {color: C.violet, width: 6, head: 18, both: true, hidden: true});
      s.lab = K.text(g, 110, 150, LAB, {cls: "t-small"});
      s.clk = K.text(g, 110, 182, "", {cls: "t-small"});
      s.kick = K.kicker(s.root, "You have the floor spectra", {x: 960, y: 120});
      s.head = K.heading(s.root, "Can you design the building now?", {x: 960, y: 168, w: 860, size: "h2"});
      s.items = [
        item(s.root, 960, 400, 860, "var(--ssi)", "<b>Forces</b> in walls and slabs"),
        item(s.root, 960, 500, 860, "var(--violet)", "<b>Displacements</b> for gaps and pipes"),
        item(s.root, 960, 600, 860, "var(--steel)", "Your detailed <b>ANSYS</b> model"),
      ];
      s.tt = T0; s.run = 0;
    },
    tick(s, t, dt) {
      clock(s, dt);
      const m = motion(s.tt), n = pose(s.ex, m);
      setForces(s.farr, n, m.F, PXK, s.ex.r);
      const roof = n[4], dx = roof[0] - 440;
      s.disp.update(440, roof[1] - 42, roof[0], roof[1] - 42);
      s.disp.g.style.opacity = String(Math.min(1, Math.abs(dx) / 30));
    },
    beats: [
      {say: "You've run the SSI analysis, and you have your floor spectra. Can you design the building now?",
        go(k) { k.s.run = 1; k.show([k.s.kick, k.s.head], {stagger: 150}); }},
      {say: "Not yet. [f]The walls and slabs need forces, [d]and the gaps and pipes need displacements.",
        f(k) { k.show([k.s.items[0], k.s.forces]); },
        d(k) { k.show([k.s.items[1], k.s.disp]); }},
      {say: "[a]And your detailed design model probably lives in ANSYS, not in SASSI.", a(k) { k.show(k.s.items[2]); }},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 8",
    build(s) {
      s.t = K.titleCard(s.root, {n: 8, part: "Design applications", title: "From SSI results to design",
        sub: "Forces, displacements, and the hand-off to ANSYS."});
    },
    beats: [
      {say: "Lesson eight. From SSI results to design: forces, displacements, and the hand-off to ANSYS.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 100}); k.show(t.title, {delay: 200}); k.show(t.sub, {delay: 350}); k.show(t.rule, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. member forces
  scenes.push({
    id: "forces", title: "Forces, with the soil inside",
    build(s) {
      const g = K.g(s.svg);
      const S = 22 / 3.2, X = 470, Y0 = 720;                       // px per ft
      s.head = K.heading(s.root, "Forces, with the soil inside", {x: 110, y: 70, size: "h2"});
      site(g, 110, Y0 + 5 * S, 720, S, 990);
      s.rip = K.ripples(g, {cx: X, cy: Y0 + 5 * S, r0: 250, r1: 355, n: 3, squash: 0.32, color: C.dash, width: 5, hidden: true});
      s.ripL = K.text(g, 820, 978, "radiation: schematic", {cls: "t-small", anchor: "end", hidden: true});
      s.ex = exStick(g, X, Y0, S, true);
      s.dims = K.g(g, {hidden: true});
      K.dim(s.dims, 745, Y0, 745, Y0 - 64 * S, "4 × 16 ft", {labelColor: "var(--ink2)"});
      K.dim(s.dims, X - 32 * S, Y0 + 5 * S + 64, X + 32 * S, Y0 + 5 * S + 64, "64 ft mat", {labelColor: "var(--ink)"});
      s.forces = K.g(g, {hidden: true});
      s.farr = [0, 1, 2, 3].map(() => K.arrow(s.forces, 0, 0, 1, 0, {color: C.ssi, width: 7, head: 20}));
      s.lab = K.text(g, 110, 182, LAB, {cls: "t-small", hidden: true});
      s.clk = K.text(g, 110, 214, "", {cls: "t-small", hidden: true});
      // at t*: the floor forces (lesson table: differences of the storey shears) and the storey shears
      const nStar = pose(s.ex, motion(T_STAR));
      s.fL = K.g(g, {hidden: true});
      FLOORF.forEach((f, i) => {
        const q = nStar[i + 1];
        K.text(s.fL, q[0] - s.ex.r - 6, q[1] + s.ex.r + 26, `${Math.round(f).toLocaleString("en-US")} kips`, {cls: "t-label", anchor: "end", color: "var(--ssi)", size: 26});
      });
      const X0 = 800, SC = 0.076;                                // px per kip
      s.barT = K.text(g, X0, 250, "peak storey shear, all at t ≈ 1.85 s", {cls: "t-label", color: "var(--ssi)", size: 30, hidden: true});
      s.bars = [];
      const bars = [3, 2, 1, 0].map((e) => ({e, ym: Y0 - (s.ex.zs[e] - 8 * S), len: SHEAR[e] * SC}));
      bars.forEach((b) => s.bars.push(K.path(g, `M${X0},${b.ym}L${X0 + b.len},${b.ym}`, {stroke: C.ssi, "stroke-width": 70, "stroke-opacity": 0.75, fill: "none", draw: true})));
      s.barL = K.g(g, {hidden: true});
      bars.forEach((b) => K.text(s.barL, X0 + b.len - 16, b.ym + 11, `${SHEAR[b.e].toLocaleString("en-US")} kips`, {cls: "t-label", anchor: "end", color: "#1a1206", size: 32, weight: 800}));
      s.card = K.card(s.root, {x: 1360, y: 170, w: 450, title: "STRESS module", size: 30,
        body: "Forces in every member, through the whole earthquake"});
      s.stat = K.stat(s.root, {x: 1360, y: 400, w: 450, value: "0", unit: "kips", label: "base shear =<br>0.77 × weight of the four floors", color: "var(--ssi)", vsize: 120});
      s.inside = K.card(s.root, {x: 1360, y: 670, w: 450, kind: "check", title: "The soil is inside", size: 30,
        body: "Soft support · rocking · energy leaking into the ground"});
      s.tt = T0; s.run = 0; s.R = 0;
      pose(s.ex, motion(T0));
    },
    tick(s, t, dt) {
      clock(s, dt);
      const m = motion(s.tt), n = pose(s.ex, m);
      setForces(s.farr, n, m.F, PXK, s.ex.r);
      s.rip.amp = s.R;
      s.rip.set(t, 1.6);
    },
    beats: [
      {say: "Take the four-storey stick from lesson 4, on its 64-foot mat, shaking on the soil.",
        go(k) { k.s.run = 1; k.show([k.s.head, k.s.lab, k.s.clk]); k.show(k.s.dims, {delay: 300}); }},
      {say: "[a]As it sways, each floor's mass pushes back, like passengers lurching when a bus brakes.",
        a(k) { k.show(k.s.forces); }},
      {say: "[s]Each storey carries the push of every floor above it, so the shear grows toward the base.",
        s(k) {
          k.s.run = 0; k.tween(k.s, {tt: T_STAR}, 600);
          k.hide(k.s.dims); k.show([k.s.barT, k.s.fL]);
          k.s.bars.forEach((b, i) => k.draw(b, 600, {delay: 150 * i}));
          k.show(k.s.barL, {delay: 400});
        }},
      {say: "[b]At the base, the peak shear is 6,764 kips: 0.77 times the weight of the four floors.",
        b(k) { k.show([k.s.card, k.s.stat.el], {stagger: 150}); k.count(k.s.stat.v, 6764, {from: 0, sep: true, ms: 1000}); k.pulse(k.s.bars[3], {amp: 0.05}); }},
      {say: "[i]And the soil is already inside: the soft support, the rocking, and energy leaking into the ground.",
        i(k) { k.show([k.s.inside, k.s.rip, k.s.ripL]); k.tween(k.s, {R: 1}, 800); k.pulse(k.s.ex.st.mat, {amp: 0.04}); }},
    ],
  });

  // ---------------------------------------------------------------- 3. equilibrium check
  scenes.push({
    id: "check", title: "Check it by equilibrium",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Check it like any ANSYS result", {x: 110, y: 70, size: "h2"});
      // free body: the top storey carries the roof (Newton: storey shear = roof mass x roof acceleration)
      const S = 20 / 3.2, x = 330, y0 = 740;                       // px per ft
      s.fb = K.g(g, {hidden: true});
      site(s.fb, 110, y0 + 5 * S, 440, S, 960);
      const ex = exStick(s.fb, x, y0, S), yr = y0 - ex.zs[3];
      s.ring = K.circle(g, x, yr, 44, {fill: "none", stroke: C.wave, "stroke-width": 5, hidden: true});
      s.cut = K.g(g, {hidden: true});
      const yc = y0 - (ex.zs[2] + ex.zs[3]) / 2, L = 130;           // equal arrows: 2,808 kips each
      K.line(s.cut, x - 120, yc, x + 120, yc, {stroke: C.ssi, "stroke-width": 5, "stroke-dasharray": "12 8"});
      K.arrow(s.cut, x + 8, yc + 30, x + 8 + L, yc + 30, {color: C.ssi, width: 7, head: 22});
      K.text(s.cut, x + 20 + L, yc + 40, "storey shear", {cls: "t-label", color: "var(--ssi)", size: 30});
      K.arrow(s.cut, x + 50, yr, x + 50 + L, yr, {color: C.wave, width: 7, head: 22});
      K.text(s.cut, x + 60, yr - 64, "roof mass", {cls: "t-label", color: "var(--wave)", size: 30});
      K.text(s.cut, x + 60, yr - 28, "× acceleration", {cls: "t-label", color: "var(--wave)", size: 30});
      // the two computed histories
      const fyi = D["ex01/BEAMS_002_00004_FYI.THS"], inert = D["ex01/roof_inertia.th"];
      s.p = K.plot(s.svg, {x: 840, y: 250, w: 950, h: 310, xr: [0, 8], yr: [-3400, 3400], xticks: [0, 2, 4, 6, 8], yticks: [-3000, 0, 3000],
        xlabel: "time (s)", ylabel: "top-storey shear (kips)", ylabelOffset: 100});
      const a = upto(fyi, 8), b = upto(inert, 8);
      s.l1 = s.p.line(a[0], a[1], {color: C.ssi, width: 8, draw: true, smooth: true});
      s.l2 = s.p.line(b[0], b[1], {color: C.wave, width: 3, draw: true, smooth: true});
      s.leg = s.p.legend([{label: "storey shear, from STRESS", color: C.ssi}, {label: "roof mass × roof acceleration", color: C.wave}],
        {x: 850, y: 170, size: 28, hidden: true});
      const pk = peak(fyi);
      s.dot = s.p.dot(pk.t, pk.v, {color: C.ssi, r: 11, hidden: true});
      s.st1 = K.stat(s.root, {x: 760, y: 700, w: 520, value: "2,811", unit: "kips", label: "2,200 kips × 1.278 g", color: "var(--wave)", vsize: 96});
      s.st2 = K.stat(s.root, {x: 1290, y: 700, w: 520, value: "2,808", unit: "kips", label: "STRESS", color: "var(--ssi)", vsize: 96});
      s.ok = K.pill(s.root, "✓ only 0.09 % apart", {x: 1120, y: 892, color: "var(--good)", size: 34});
    },
    beats: [
      {say: "Can you trust a number like that? Check it the way you'd check any ANSYS result: by equilibrium.",
        go(k) { k.show(k.s.head); k.show(k.s.fb, {delay: 300}); }},
      {say: "[t]The top storey holds up only the roof, so its shear must equal roof mass times roof acceleration.",
        t(k) { k.show([k.s.ring, k.s.cut]); }},
      {say: "[p]SASSI gives you both histories, and they lie right on top of each other.",
        p(k) { k.show(k.s.p.g); k.draw(k.s.l1, 1000); k.draw(k.s.l2, 1000, {delay: 300}); k.show(k.s.leg, {delay: 300}); }},
      {say: "[n]At the peak, 2,200 kips times 1.278 g gives 2,811 kips. [s]STRESS says 2,808.",
        n(k) { k.show([k.s.dot, k.s.st1.el]); }, s(k) { k.show(k.s.st2.el); k.show(k.s.ok, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 4. relative displacements
  scenes.push({
    id: "displacements", title: "Displacements: compared to what?",
    build(s) {
      const g = K.g(s.svg);
      const S = 20 / 3.2, X = 440, Y0 = 694;                       // px per ft
      s.head = K.heading(s.root, "How far does the roof move?", {x: 110, y: 70, size: "h2"});
      site(g, 110, Y0 + 5 * S, 690, S, 990);
      s.ffl = K.line(g, X, 262, X, Y0 + 5 * S, {stroke: C.wave, "stroke-width": 4, "stroke-dasharray": "12 10"});
      s.ffLab = K.text(g, X - 12, 254, "free field", {cls: "t-label", anchor: "end", color: "var(--wave)", size: 28});
      s.matl = K.line(g, X, 262, X, Y0, {stroke: C.ssi, "stroke-width": 3, "stroke-dasharray": "4 8", hidden: true});
      s.matLab = K.text(g, X + 12, 254, "mat centre", {cls: "t-label", color: "var(--ssi)", size: 28, hidden: true});
      s.ex = exStick(g, X, Y0, S, true);
      s.lab = K.text(g, 110, 176, LAB, {cls: "t-small"});
      s.clk = K.text(g, 110, 208, "", {cls: "t-small"});
      s.drift = K.g(g, {hidden: true});
      s.drA = K.arrow(s.drift, 0, 0, 1, 0, {color: C.ssi, width: 5, head: 15, both: true});
      s.drL = K.line(s.drift, 0, 0, 0, 0, {stroke: C.ssi, "stroke-width": 3, "stroke-dasharray": "6 6"});
      s.drT = K.text(s.drift, 0, 0, "drift", {cls: "t-label", color: "var(--ssi)", size: 30});
      // the answers (lesson numbers: 0.0766, 0.0835 and 0.0074 ft)
      s.rows = [
        row(s.root, 880, 220, "0.92", "in", "roof, relative to <b style='color:var(--ssi)'>its own mat</b>", "var(--ssi)"),
        row(s.root, 880, 400, "1.00", "in", "roof, relative to <b style='color:var(--wave)'>the free field</b>", "var(--wave)"),
        row(s.root, 880, 580, "0.09", "in", "the mat centre itself:<br>it slides, and it rocks", "var(--kin)"),
      ];
      // which reference
      s.c1 = K.card(s.root, {x: 900, y: 170, w: 910, title: "Drift between floors", size: 34, body: "Floor minus floor: <b>either reference</b> works."});
      s.c2 = K.card(s.root, {x: 900, y: 380, w: 910, title: "Gap or pipe to the next building", size: 34, body: "Use the <b>free field</b>."});
      s.pan = K.g(g, {hidden: true});
      K.line(s.pan, 960, 900, 1780, 900, {stroke: C.wave, "stroke-width": 5});
      const b1 = K.building(s.pan, {x: 1170, y: 900, w: 150, storeys: 4, storeyH: 70, windows: false});
      const b2 = K.building(s.pan, {x: 1560, y: 900, w: 130, storeys: 3, storeyH: 70, windows: false});
      const a = b1.T(75, 90), b = b2.T(-65, 90), ga = b1.T(75, 192), gb = b2.T(-65, 192);
      K.path(s.pan, `M${a[0]},${a[1]}C${a[0] + 60},${a[1] + 40} ${b[0] - 60},${b[1] + 40} ${b[0]},${b[1]}`, {stroke: C.violet, "stroke-width": 9, fill: "none", "stroke-linecap": "round"});
      K.arrow(s.pan, ga[0] + 6, ga[1], gb[0] - 6, gb[1], {color: C.ink2, width: 4, head: 14, both: true});
      K.text(s.pan, (ga[0] + gb[0]) / 2, ga[1] - 20, "gap", {cls: "t-label", anchor: "middle", size: 30});
      K.text(s.pan, (a[0] + b[0]) / 2, a[1] + 60, "pipe", {cls: "t-label", anchor: "middle", color: "var(--violet)", size: 30});
      K.text(s.pan, 1780, 940, "free field", {cls: "t-small", anchor: "end", color: "var(--wave)"});
      K.text(s.pan, 960, 940, "schematic", {cls: "t-small"});
      s.tt = T0; s.run = 0;
    },
    tick(s, t, dt) {
      clock(s, dt);
      const m = motion(s.tt), n = pose(s.ex, m), xm = s.ex.x + m.u0 * s.ex.S * EXAG * s.ex.ff;
      s.matl.setAttribute("x1", xm); s.matl.setAttribute("x2", xm);
      // drift between floors 2 and 3
      const p2 = n[2], p3 = n[3], yA = p3[1] - 36;
      s.drL.setAttribute("x1", p2[0]); s.drL.setAttribute("y1", p2[1]); s.drL.setAttribute("x2", p2[0]); s.drL.setAttribute("y2", yA);
      const dx = p3[0] - p2[0];
      if (Math.abs(dx) > 34) s.drA.update(p2[0], yA, p3[0], yA); else s.drA.update(p2[0] - 17, yA, p2[0] + 17, yA);
      s.drA.g.style.opacity = String(Math.min(1, Math.abs(dx) / 34));
      s.drT.setAttribute("x", Math.max(p2[0], p3[0]) + 70); s.drT.setAttribute("y", yA + 10);
    },
    beats: [
      {say: "Next, displacements. How far does the roof move? First ask: compared to what?",
        go(k) { k.show(k.s.head); k.s.run = 1; }},
      {say: "[t]It's like walking on a train: your speed relative to the train, or to the platform?",
        t(k) { k.show([k.s.matl, k.s.matLab]); k.pulse([k.s.ffLab, k.s.matLab]); }},
      {say: "[m]SASSI's RELDISP module answers both. Relative to its own mat, the roof moves 0.92 inches.",
        m(k) { k.show(k.s.rows[0]); k.tween(k.s.ex, {ff: 0}, 600); k.pulse(k.s.matLab); }},
      {say: "[f]Relative to the free field, the ground shaking with no building on it, the roof moves {1.00 inch|a full inch}.",
        f(k) { k.show(k.s.rows[1]); k.tween(k.s.ex, {ff: 1}, 600); k.pulse(k.s.ffLab); }},
      {say: "[r]Why the difference? The mat itself moves on the soil: its centre slides up to 0.09 inches, and it rocks.",
        r(k) { k.show(k.s.rows[2]); k.pulse(k.s.ex.st.mat, {amp: 0.08}); }},
      {say: "So pick the reference for the question. [d]For drift between floors, subtract two floors: either reference works.",
        go(k) { k.hide(k.s.rows); }, d(k) { k.show([k.s.c1, k.s.drift]); }},
      {say: "[g]For the gap to the next building, or a pipe running between buildings, use the free field.",
        g(k) { k.show([k.s.c2, k.s.pan]); }},
    ],
  });

  // ---------------------------------------------------------------- 5. the hand-off (Option A)
  scenes.push({
    id: "handoff", title: "Hand the motion to ANSYS",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "Hand the motion to ANSYS", {x: 110, y: 70, size: "h2"});
      const badge = (p, x, y, n, color) => {
        K.circle(p, x, y, 26, {fill: color});
        K.text(p, x, y + 10, n, {cls: "t-label", anchor: "middle", color: "#04111c", size: 30, weight: 800});
      };
      // step 1: the SSI model (example 1 to scale, computed motion)
      s.L = K.g(g, {hidden: true});
      badge(s.L, 180, 220, "1", C.wave);
      K.text(s.L, 222, 231, "SSI model in SASSI: simple", {cls: "t-label", color: "var(--wave)", size: 32});
      const S1 = 13 / 3.2;                                          // px per ft
      site(s.L, 140, 570 + 5 * S1, 500, S1, 760);
      s.ex = exStick(s.L, 390, 570, S1, true);
      s.lab = K.text(s.L, 140, 272, "computed: × 250, slowed down 4 ×", {cls: "t-small"});
      s.clk = K.text(s.L, 140, 790, "", {cls: "t-small"});
      // step 2: a detailed ANSYS model of the same building: 64 ft x 64 ft walls meshed at 8 ft, slabs every 16 ft
      s.R = K.g(g, {hidden: true});
      badge(s.R, 1220, 220, "2", C.ssi);
      K.text(s.R, 1262, 231, "ANSYS model: detailed", {cls: "t-label", color: "var(--ssi)", size: 32});
      const S2 = 5, MX = 1290, MB = 600, E = 8 * S2;                // px per ft
      K.mesh(s.R, {x: MX, y: MB - 64 * S2, cols: 8, rows: 8, dx: E, dy: E, color: "rgba(185,196,208,.45)"});
      for (let j = 1; j <= 4; j++) K.line(s.R, MX - 8, MB - 16 * j * S2, MX + 64 * S2 + 8, MB - 16 * j * S2, {stroke: C.concrete, "stroke-width": 6, "stroke-linecap": "round"});
      K.rect(s.R, MX, MB, 64 * S2, 5 * S2, {fill: "var(--concrete-fill)", stroke: C.concrete, "stroke-width": 3});
      K.text(s.R, MX + 64 * S2, MB - 64 * S2 - 14, "illustration", {cls: "t-small", anchor: "end"});
      // F: the inertia forces at the floors (LOADGEN, -X); D: the interface displacements (x 600)
      s.inF = K.g(s.R, {hidden: true});
      LG_F.forEach((f, i) => {
        const y = MB - 16 * (i + 1) * S2, x0 = MX + 64 * S2 + 22, len = f * PXK;
        K.arrow(s.inF, x0 + len, y, x0, y, {color: C.ssi, width: 6, head: 18});
      });
      s.ifD = K.g(s.R, {hidden: true});
      const kD = S2 * 600, yb = MB + 5 * S2;
      const dp = LG_UZ.map((uz, i) => [MX + i * E + LG_UX * kD, yb - uz * kD]);
      LG_UZ.forEach((uz, i) => { K.line(s.ifD, MX + i * E, yb, dp[i][0], dp[i][1], {stroke: C.wave, "stroke-width": 2, "stroke-opacity": 0.7}); K.circle(s.ifD, MX + i * E, yb, 5, {fill: C.muted}); });
      K.poly(s.ifD, dp, {stroke: C.wave, "stroke-width": 5, fill: "none"});
      dp.forEach((q) => K.circle(s.ifD, q[0], q[1], 7, {fill: C.wave, stroke: "#0a111d", "stroke-width": 2}));
      // the motion, from step 1 to step 2
      s.arr = K.arrow(g, 690, 420, 1150, 420, {color: C.wave, width: 6, head: 24, dash: "16 12", hidden: true});
      s.motT = K.text(g, 920, 470, "the motion", {cls: "t-label", anchor: "middle", color: "var(--wave)", size: 34, hidden: true});
      s.motT2 = K.text(g, 920, 506, "of the masses and the 81 interface nodes", {cls: "t-small", anchor: "middle", hidden: true});
      // Option A, LOADGEN
      s.card = K.card(s.root, {x: 110, y: 830, w: 1060, kind: "ansys", title: "ACS SASSI: Option A", size: 34,
        body: "The <b>LOADGEN</b> module writes the loads as an ANSYS input file."});
      s.file = K.pill(s.root, "→ ex01_LGS.inp", {x: 1240, y: 900, color: "var(--good)", size: 34});
      s.file.style.fontFamily = "var(--mono)";
      s.loadT = K.label(s.root, `<span style='color:var(--ssi)'>F: mass × acceleration at t = ${T_STAR} s</span>`, {x: 1240, y: 700, size: 28, in: "fade"});
      s.loadT2 = K.label(s.root, "<span style='color:var(--wave)'>D: the 81 interface nodes move (× 600)</span>", {x: 1240, y: 744, size: 28, in: "fade"});
      s.tt = T0; s.run = 0;
    },
    tick(s, t, dt) {
      clock(s, dt);
      pose(s.ex, motion(s.tt));
    },
    beats: [
      {say: "Now the hand-off. Your SSI model is kept simple, while your ANSYS design model has all the detail.",
        go(k) { k.show(k.s.head); k.show(k.s.L, {delay: 200}); k.show(k.s.R, {delay: 400}); }},
      {say: "[a]Step one: SASSI shakes the building on the soil, and works out how it moves.",
        a(k) { k.s.run = 1; k.pulse(k.s.ex.st.g, {amp: 0.03}); }},
      {say: "[b]Step two: it hands that motion, like a baton, to your ANSYS model, which works out the detailed forces.",
        b(k) { k.show(k.s.arr); k.show([k.s.motT, k.s.motT2], {delay: 300}); }},
      {say: "ACS SASSI calls this Option A. [c]Its LOADGEN module writes the loads as an ANSYS input file.",
        c(k) { k.show(k.s.card); k.show(k.s.file, {delay: 400}); }},
      {say: "[l]The loads, frozen at one instant: a push at each floor, its mass times acceleration, [u]plus how the foundation moves.",
        l(k) { k.s.run = 0; k.tween(k.s, {tt: T_STAR}, 600); k.show([k.s.inF, k.s.loadT]); },
        u(k) { k.show([k.s.ifD, k.s.loadT2]); }},
    ],
  });

  // ---------------------------------------------------------------- 6. one instant, and the check
  scenes.push({
    id: "instant", title: "One instant, checked",
    build(s) {
      s.head = K.heading(s.root, "Which instant? Does it add up?", {x: 110, y: 70, size: "h2"});
      const v1 = D["ex01/BEAMS_002_00001_FYI.THS"], c = upto(v1, 12), pk = peak(v1);
      s.p = K.plot(s.svg, {x: 250, y: 220, w: 940, h: 300, xr: [0, 12], yr: [-8000, 8000], xticks: [0, 2, 4, 6, 8, 10, 12], yticks: [-6000, 0, 6000],
        xlabel: "time (s)", ylabel: "base shear (kips)", ylabelOffset: 100});
      s.lv = s.p.line(c[0], c[1], {color: C.ssi, width: 4, draw: true, smooth: true});
      s.glow = s.p.line(c[0], c[1], {color: C.ssi, width: 4, smooth: true, glow: true, hidden: true});
      s.tstar = s.p.vline(pk.t, {color: C.ink, width: 4, hidden: true});
      s.tsl = s.p.text(pk.t, 8000, `largest: ${pk.t} s`, {cls: "t-label", color: "var(--ink)", size: 30, dx: 12, dy: -16, hidden: true});
      s.pLab = s.p.text(12, 8000, "base shear of the stick (STRESS)", {cls: "t-small", anchor: "end", color: "var(--ssi)", size: 26, dy: -16});
      s.code = K.code(s.root, ["LOADGEN,DISPACC,0,0,0,LUMPED,1,FILE8", "LGTIME,V,1", "RUNLOADGEN,STATIC"], {x: 110, y: 700, w: 720, size: 28});
      s.codeL = ["what to write", "which instant", "go"].map((t, i) =>
        K.label(s.root, "← " + t, {x: 850, y: 752 + i * 42 - 2, size: 28, color: "var(--muted)", in: "fade"}));
      s.st1 = K.stat(s.root, {x: 1290, y: 200, w: 520, value: "6,772", unit: "kips", label: "the floor loads, added up", color: "var(--ssi)", vsize: 100});
      s.st2 = K.stat(s.root, {x: 1290, y: 400, w: 520, value: "6,764", unit: "kips", label: "STRESS, same instant", color: "var(--ink)", vsize: 100});
      s.ok = K.pill(s.root, "✓ within 0.12 %", {x: 1440, y: 590, color: "var(--good)", size: 34});
      s.dyn = K.card(s.root, {x: 1110, y: 700, w: 700, title: "The whole earthquake?", size: 32,
        body: "LOADGEN also writes <b>full histories</b> for an ANSYS transient."});
    },
    beats: [
      {say: "Which instant does LOADGEN freeze? [t]By default, the moment of the largest base shear, about 1.85 seconds in.",
        go(k) { k.show(k.s.head); k.show(k.s.p.g); k.draw(k.s.lv, 1000); }, t(k) { k.show([k.s.tstar, k.s.tsl]); }},
      {say: "[c]In SASSI, that's three short lines: what to write, which instant, and go.",
        c(k) { k.show(k.s.code.el); k.type(k.s.code.lines, {cps: 44}); k.show(k.s.codeL, {stagger: 150, delay: 400}); }},
      {say: "[s]The floor loads add up to a base shear of 6,772 kips. [r]STRESS gave 6,764 at that instant: within 0.12 percent.",
        s(k) { k.show(k.s.st1.el); }, r(k) { k.show(k.s.st2.el); k.show(k.s.ok, {delay: 400}); }},
      {say: "[d]Need the whole earthquake, not one instant? LOADGEN can also write full histories, for an ANSYS transient.",
        d(k) { k.show(k.s.dyn); k.dim([k.s.tstar, k.s.tsl], true); k.show(k.s.glow); }},
    ],
  });

  // ---------------------------------------------------------------- 7. the one rule
  // Illustration computed with the sway-rocking model of lesson 4 (PH.ssiModel: mode 1 of the stick, the
  // mat, SASSI's impedance of this mat): steady shaking at the model's SSI peak, free field +-0.08 in,
  // displacements x 250; the ANSYS building made 4 x softer (fixed-base frequency halved).  PU converts ft to
  // the length unit of PH.SSI (its mat width B is the 64 ft mat) for the rocking angle (th per unit u_g).
  const UG = 0.08 / 12, PU = PH.SSI.B / 64;
  const PHI = [0.244, 0.581, 0.956, 1.326];  // floor deformation per unit modal displacement at the modal height h1 (MAT_ at t*, rocking removed)
  scenes.push({
    id: "rule", title: "The one rule",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "The one rule", {x: 110, y: 70, size: "h2"});
      s.rule = K.para(s.root, "Your ANSYS model must <b style='color:var(--ink)'>not change how the foundation moves</b>.", {x: 110, y: 175, w: 1700, size: 46});
      s.dg = K.g(g, {hidden: true});
      const S = 18 / 3.2, X = 960, Y0 = 800;                       // px per ft
      site(s.dg, 660, Y0 + 5 * S, 600, S, 985);
      s.ex = exStick(s.dg, X, Y0, S, true);
      s.ghost = K.path(g, "", {fill: "none", stroke: C.wave, "stroke-width": 4, "stroke-dasharray": "10 7", hidden: true});
      s.ghL = K.text(g, 1270, Y0 + 34, "motion from SASSI", {cls: "t-small", color: "var(--wave)", hidden: true});
      s.lab = K.text(g, X, 392, "same building: same motion", {cls: "t-label", anchor: "middle", color: "var(--good)", size: 32, hidden: true});
      s.lab2 = K.text(g, X, 392, "4 × softer: different motion", {cls: "t-label", anchor: "middle", color: "var(--bad)", size: 32, hidden: true});
      // the model and its SSI peak
      const m1 = PH.ssiModel(1);
      let fpk = 2, best = 0;
      for (let f = 2; f <= 6; f += 0.01) { const a = Math.hypot(...PH.ssiResponse(m1, f).Ht); if (a > best) { best = a; fpk = f; } }
      s.f = fpk; s.r1 = PH.ssiResponse(m1, fpk);
      s.ill = K.text(s.dg, X, 330, `illustration: sway-rocking model of lesson 4, ${fpk.toFixed(1)} Hz, free field ±0.08 in, × 250, slowed down 6 ×`, {cls: "t-small", anchor: "middle"});
      s.c1 = K.card(s.root, {x: 110, y: 420, w: 500, kind: "check", title: "Same building, more detail", size: 32,
        body: "The motion fits: a replay recovers SASSI's result almost exactly."});
      s.c2 = K.card(s.root, {x: 1310, y: 420, w: 500, kind: "warn", title: "Much softer or heavier", size: 32,
        body: "The motion no longer fits: <b>rerun the SSI</b> with the new properties."});
      s.soft = 0; s.run = 0; s.ph = 0;
    },
    tick(s, t, dt) {
      if (s.run) s.ph += TAU * s.f * dt / 6;
      const re = (z) => z[0] * Math.cos(s.ph) - z[1] * Math.sin(s.ph);
      const kap = 1 - 0.75 * s.soft, r = s.soft ? PH.ssiResponse(PH.ssiModel(1, {ffb: PH.SSI.ffb * Math.sqrt(kap)}), s.f) : s.r1;
      const ex = s.ex, kp = ex.S * EXAG * UG;
      ex.st.set({sway: re(r.u0) * kp, rock: re(r.th) * UG * PU * EXAG, bend: PHI.map((p) => re(r.u) * p * kp)});
      glue(ex);
      // the handed-over interface motion (the original building)
      const cx = ex.x + re(s.r1.u0) * kp, cy = ex.y + ex.matH / 2, a = re(s.r1.th) * UG * PU * EXAG, c = Math.cos(a), sn = Math.sin(a);
      const w = 32 * ex.S, h = ex.matH / 2;
      const P = [[-w, h], [w, h], [w, -h], [-w, -h]].map(([lx, ly]) => `${(cx + lx * c - ly * sn).toFixed(1)},${(cy + lx * sn + ly * c).toFixed(1)}`);
      s.ghost.setAttribute("d", `M${P.join("L")}Z`);
    },
    beats: [
      {say: "There's one rule: [r]your ANSYS model must not change how the foundation moves.",
        go(k) { k.show(k.s.head); k.show(k.s.dg); k.s.run = 1; }, r(k) { k.show(k.s.rule); }},
      {say: "[s]Same building, just more detail? Then the motion fits: a replay recovers SASSI's own result almost exactly.",
        s(k) { k.show([k.s.c1, k.s.lab, k.s.ghost, k.s.ghL]); }},
      {say: "[d]But a much softer or heavier ANSYS building would shake the soil differently, so the handed-over motion is wrong.",
        d(k) { k.tween(k.s, {soft: 1}, 1000); k.hide(k.s.lab); k.show(k.s.lab2, {delay: 200}); k.show(k.s.c2, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 8. recap
  scenes.push({
    id: "recap", title: "Recap",
    build(s) {
      s.head = K.heading(s.root, "Recap", {x: 110, y: 80, size: "h1"});
      s.list = K.bullets(s.root, [
        {t: "<b>STRESS</b>: member forces, soil inside", sub: "check them by equilibrium"},
        {t: "<b>Displacements</b> need a reference", sub: "pick the one your question needs"},
        {t: "<b>LOADGEN</b> hands the motion to ANSYS", sub: "while the foundation still moves the same"},
      ], {x: 110, y: 230, w: 980, num: true, size: 42});
      s.q = K.card(s.root, {x: 1170, y: 230, w: 640, kind: "check", title: "Check yourself", size: 34,
        body: "LOADGEN froze the instant of the largest base shear. Is it the worst instant for every member?"});
      s.a = K.card(s.root, {x: 1170, y: 580, w: 640, title: "Answer", size: 34,
        body: "<b>No.</b> Other forces peak at other times: add instants, or use full histories."});
      s.next = K.pill(s.root, "Next · Lesson 9: nonlinear soil and cracked concrete →", {x: 110, y: 900, color: "var(--wave)", size: 32});
    },
    beats: [
      {say: "Let's recap. [a]STRESS gives member forces with the soil already inside; check them by equilibrium.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[a]Displacements need a reference: pick the one your question needs.", a(k) { k.show(k.s.list.items[1]); }},
      {say: "[a]And LOADGEN hands the SSI motion to ANSYS, as long as the foundation still moves the same.",
        a(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Check yourself. LOADGEN froze the instant of the largest base shear. Is that the worst instant for every member?",
        q(k) { k.show(k.s.q); }, gap: 1500},
      {say: "[a]No. A wall stress or a mid-height bending moment can peak at another time. Add instants, or use full histories.",
        a(k) { k.show(k.s.a); }},
      {say: "[n]Next, in lesson 9: what happens when the soil softens and the concrete cracks.", n(k) { k.show(k.s.next); }},
    ],
  });

  SV.video({
    id: "08", n: 8, part: "Design applications", lesson: "08-design-outputs",
    title: "From SSI results to design: forces, displacements, ANSYS",
    subtitle: "Member forces with the soil inside, relative displacements and their reference, and the Option A hand-off to ANSYS with LOADGEN.",
    next: {href: "09.html", title: "Nonlinear soil and cracked concrete by iteration"},
    pron: [[/\bSTRESS\b/g, "stress"]],
    scenes,
  });
})();
