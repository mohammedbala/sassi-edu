/* Explainer video, lesson 10: Good practice, checks and the limits of SSI models (third edition: pace and accuracy).
 * One idea: an SSI run can finish without an error and still be wrong; four quick checks catch it.
 * Numbers: sassi/ui/lessons/10_good_practice.md; curves: d10.js (the lesson's runs: example 1, the coarse
 * frequency set ex01c and the 10 Hz cut-off ex01f; python -m sassi.ui.video_data).
 * Pictures: example 1 to scale (64 ft x 64 ft mat, 9 x 9 nodes at 8 ft, the stick on node 41, four 16 ft storeys,
 * 16 ft sand and 39 ft gravel over rock); the slow-shake motion from the computed transfer functions at 0.098 Hz;
 * every plot annotation from the data it names (CRITFREQ with SV.K.PH.critfreq). */
"use strict";

(function () {
  const K = SV.K, C = K.C, D = SV.DATA["10"], PH = K.PH;
  const TAU = 2 * Math.PI;
  const DF = 1 / (8192 * 0.005);                      // frequency step of example 1 (Hz): f = n df

  // ---------------------------------------------------------------- data (the lesson's runs)
  const R_TFI = D["ex01/00085TR_X.TFI"], R_TFU = D["ex01/00085TR_X.TFU"];                   // roof, 22 SSI frequencies
  const FLOOR_TFU = [41, 82, 83, 84, 85].map((n) => D[`ex01/${String(n).padStart(5, "0")}TR_X.TFU`]);   // mat centre, floors 1-4
  const M_TFI = D["ex01/00041TR_X.TFI"], M_TFU = D["ex01/00041TR_X.TFU"];                   // mat centre, cut-off 20 Hz
  const MF_TFI = D["ex01f/00041TR_X.TFI"];                                                   // mat centre, cut-off 10 Hz
  const CO_TFI = D["ex01c/COARSE_00085TR_X.TFI"], CO_TFU = D["ex01c/COARSE_00085TR_X.TFU"];  // roof, 11 frequencies
  const AD_TFU = D["ex01c/00085TR_X.TFU"];                                                   // 11 + the flagged one
  const RS_R = D["ex01/00085TR_X02.RS"], RS_RF = D["ex01f/00085TR_X02.RS"];                  // roof ISRS, 5 %, cut-off 20 / 10 Hz
  const RS_M = D["ex01/00041TR_X02.RS"], RS_MF = D["ex01f/00041TR_X02.RS"];                  // mat ISRS, 5 %, cut-off 20 / 10 Hz

  // ---------------------------------------------------------------- helpers
  /** Points of a curve with x in [a, b]. */
  function clipXY(xs, ys, a, b) {
    const X = [], Y = [];
    xs.forEach((x, i) => { if (x >= a - 1e-9 && x <= b + 1e-9) { X.push(x); Y.push(ys[i]); } });
    return [X, Y];
  }
  /** Linear interpolation of ys(xs) at x (in log x when log). */
  function lerp(xs, ys, x, log) {
    const t = log ? Math.log10 : (v) => v;
    if (x <= xs[0]) return ys[0];
    for (let i = 1; i < xs.length; i++) {
      if (xs[i] >= x) { const u = (t(x) - t(xs[i - 1])) / (t(xs[i]) - t(xs[i - 1])); return ys[i - 1] + u * (ys[i] - ys[i - 1]); }
    }
    return ys[ys.length - 1];
  }
  const argmax = (a) => a.reduce((b, v, i) => (v > a[b] ? i : b), 0);
  /** Computed SSI frequencies as hollow dots on a plot (hidden, pop in). */
  function dots(p, d, o) {
    o = o || {};
    const g = K.g(p.g, {});
    return d.f.filter((f) => f <= (o.b || 1e9) + 1e-9).map((f, i) =>
      K.circle(g, p.X(f), p.Y(d.amp[i]), o.r || 10, {fill: "#0a111d", stroke: o.color || C.ssi, "stroke-width": 4.5, hidden: true, in: "pop"}));
  }
  /** A ring around a point of a plot, with an optional label beside it (hidden group). */
  function ring(p, x, y, label, o) {
    o = o || {};
    const g = K.g(p.g, {hidden: true, in: "pop"});
    K.circle(g, p.X(x), p.Y(y), o.r || 24, {fill: "none", stroke: o.color, "stroke-width": 5});
    if (label) K.text(g, p.X(x) + (o.dx || 36), p.Y(y) + (o.dy || 10), label, {cls: "t-label", anchor: o.anchor || "start", color: o.color, size: o.size || 30});
    return g;
  }
  /** A program listing (what SASSI-EDU prints): a title and rows of monospace text, each row hidden. */
  function listing(root, o) {
    const el = K.h("div", {x: o.x, y: o.y, w: o.w, in: "fade", cls: "sv-code", style: {fontSize: (o.size || 25) + "px", lineHeight: "1.42"}}, root);
    K.h("div", {cls: "sv-code-h", text: o.title, visible: true}, el).classList.remove("sv-abs");
    const rows = o.rows.map((r) => {
      const d = K.h("div", {html: r, in: "left", style: {whiteSpace: "pre-wrap"}}, el);
      d.classList.remove("sv-abs");
      return d;
    });
    return {el, rows};
  }
  const cmd = (t) => `<span style="color:var(--wave);font-weight:700">${t}</span>`;
  const okT = (t) => `<span style="color:var(--good)">${t}</span>`;
  const warnT = (t) => `<span style="color:var(--ssi)">${t}</span>`;
  /** Readable SVG text over a drawing (dark halo). */
  function halo(e) { e.style.paintOrder = "stroke"; e.style.stroke = "#0a111d"; e.style.strokeWidth = "7px"; e.style.strokeLinejoin = "round"; return e; }

  // ================================================================== scenes
  const scenes = [];
  const MODS = ["SITE", "POINT", "HOUSE", "ANALYS", "MOTION", "STRESS", "RELDISP"];

  // ---------------------------------------------------------------- 0. hook
  scenes.push({
    id: "hook", title: "It ran. Is it right?",
    build(s) {
      const g = K.g(s.svg);
      s.chain = K.chain(s.root, MODS, {x: 130, y: 96, bw: 200, bh: 84, gap: 40, size: 26});
      s.ticks = MODS.map((m, i) => K.h("div", {x: 130 + i * 240, y: 192, w: 200, align: "center", size: 40, color: "var(--good)", html: "✓", in: "pop",
        style: {fontWeight: "800"}}, s.root));
      s.p = K.plot(s.svg, {x: 230, y: 400, w: 660, h: 380, xr: [0.5, 50], xlog: true, yr: [0, 9], xticks: [0.5, 1, 2, 5, 10, 20, 50], yticks: [0, 3, 6, 9],
        xlabel: "equipment frequency (Hz)", ylabel: "roof spectrum, 5 % (g)", ylabelOffset: 70});
      s.rs = s.p.line(RS_R.f, RS_R.sa, {color: C.ssi, width: 6, draw: true});
      s.ok = K.text(g, 680, 520, "looks reasonable…", {cls: "t-label", anchor: "middle", color: "var(--ink2)", hidden: true, size: 32});
      // a plot cursor reading the computed spectrum
      s.cur = K.g(s.p.g, {hidden: true});
      s.curL = K.line(s.cur, 0, s.p.box.y, 0, s.p.box.y + s.p.box.h, {stroke: C.ink2, "stroke-width": 2, "stroke-dasharray": "6 6"});
      s.curD = K.circle(s.cur, 0, 0, 10, {fill: C.ink, stroke: C.ssi, "stroke-width": 4});
      s.curT = K.text(s.cur, 0, s.p.box.y - 14, "", {cls: "t-label", anchor: "middle", size: 26, color: "var(--ink)"});
      s.q = K.heading(s.root, "It ran.<br>Is it right?", {x: 1060, y: 330, size: "h1"});
      s.bad = ["a missing frequency", "a badly chosen cut-off", "a hidden hinge"].map((t, i) =>
        K.pill(s.root, "✗ " + t, {x: 1060, y: 540 + i * 80, color: "var(--bad)", style: {fontSize: "30px"}}));
      s.change = K.heading(s.root, "No error. But your floor spectra change.", {x: 1060, y: 780, w: 760, size: "h3", color: "var(--bad)"});
      s.good = K.pill(s.root, "✓ four quick checks catch them", {x: 1060, y: 910, color: "var(--good)", style: {fontSize: "32px"}});
      s.L = 0;
    },
    tick(s, t) {
      if (!s.L) return;
      const u = 0.5 - 0.5 * Math.cos(TAU * t / 7);
      const f = 0.8 * Math.pow(30 / 0.8, u);
      const sa = lerp(RS_R.f, RS_R.sa, f, true), X = s.p.X(f);
      s.curL.setAttribute("x1", X); s.curL.setAttribute("x2", X);
      s.curD.setAttribute("cx", X); s.curD.setAttribute("cy", s.p.Y(sa));
      s.curT.setAttribute("x", Math.max(s.p.box.x + 110, Math.min(s.p.box.x + s.p.box.w - 110, X)));
      s.curT.textContent = `${f.toFixed(1)} Hz · ${sa.toFixed(2)} g`;
    },
    beats: [
      {say: "Your SSI run just finished. [c]Every module ran to the end, without a single error.",
        go(k) { k.show(k.s.chain.boxes, {stagger: 110}); k.show(k.s.chain.arrows, {stagger: 110, delay: 60}); },
        c(k) {
          const s = k.s;
          s.chain.boxes.forEach((b, i) => k.at(i * 120, () => { b.style.borderColor = "var(--good)"; }));
          k.show(s.ticks, {stagger: 120});
        }},
      {say: "[p]Out comes the roof spectrum: how hard equipment up there gets shaken. It looks reasonable.",
        p(k) { k.show(k.s.p.g); k.draw(k.s.rs, 1000, {delay: 200}); k.show(k.s.ok, {delay: 400}); }},
      {say: "[a]A missing frequency, [b]a badly chosen cut-off, [c]a hidden hinge: none of them stops the run.",
        go(k) { k.show(k.s.q); k.show(k.s.cur); k.tween(k.s, {L: 1}, 10); },
        a(k) { k.show(k.s.bad[0]); }, b(k) { k.show(k.s.bad[1]); }, c(k) { k.show(k.s.bad[2]); }},
      {say: "[x]But every one of them changes your floor spectra.", x(k) { k.show(k.s.change); k.hide(k.s.ok); }},
      {say: "[f]The good news: four quick checks catch them.", f(k) { k.show(k.s.good); k.pulse(k.s.good); }},
    ],
  });

  // ---------------------------------------------------------------- 1. title
  scenes.push({
    id: "title", title: "Lesson 10",
    build(s) {
      s.t = K.titleCard(s.root, {n: 10, part: "Advanced", title: "Good practice, checks and the limits of SSI models",
        sub: "Four quick checks before you trust a floor spectrum."});
    },
    beats: [
      {say: "Lesson ten: good practice, checks, and the limits of SSI models.",
        go(k) { const t = k.s.t; k.show(t.num); k.show(t.part, {delay: 100}); k.show(t.title, {delay: 200}); k.show(t.sub, {delay: 300}); k.show(t.rule, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 2. check 1: before you run
  // example 1 to scale (4.69 px/ft): elevation (mat 64 ft x 5 ft, stick nodes at 16, 32, 48, 64 ft above the mat's
  // mid-plane) and plan (9 x 9 nodes at 8 ft)
  const E1 = {S: 300 / 64, ex: 1135, gy: 680, px: 1625, py: 530};
  scenes.push({
    id: "model-checks", title: "Check 1: before you run",
    build(s) {
      const g = K.g(s.svg), S = E1.S;
      s.head = K.heading(s.root, "Check 1 · Before you run", {x: 110, y: 70, size: "h2"});
      s.list = listing(s.root, {x: 110, y: 190, w: 820, title: "SASSI-EDU · model checks of example 1", rows: [
        cmd("CHECK"),
        "  MODEL: " + okT("0 errors, 0 warnings"),
        "  SITE … RELDISP: " + okT("0 errors, 0 warnings"),
        cmd("FIXEDINT"),
        "  " + okT("no interaction node has a fixed translation"),
        cmd("HINGED"),
        "  " + warnT("1 possible hinge"),
        "  " + warnT("node 41: beam 1 (group 2) meets the coplanar\n  shells out of their plane: drilling (θz)\n  is not transmitted"),
        "  " + okT("→ explained: the rigid spider (group 3)\n    carries θz into the mat"),
      ]});
      // elevation
      const EL = K.g(g, {hidden: true}); s.elev = EL;
      const ex = E1.ex, gy = E1.gy;
      K.text(EL, ex, 318, "elevation", {cls: "t-small", anchor: "middle", size: 26});
      K.soil(EL, {x: ex - 165, y: gy, w: 330, layers: [{h: 40, kind: "sand"}], labels: false});
      s.stick = K.stick(EL, {x: ex, y: gy - 5 * S, z: [16, 32, 48, 64].map((h) => h * S - 2.5 * S), matW: 64 * S, matH: 5 * S, r: 13, slabW: 0});
      K.dim(EL, ex - 32 * S, gy + 72, ex + 32 * S, gy + 72, "64 ft");
      K.dim(EL, ex + 32 * S + 22, gy - 2.5 * S, ex + 32 * S + 22, gy - 2.5 * S - 64 * S, "4 × 16 ft");
      s.elSpider = K.line(g, ex - 8 * S, gy - 2.5 * S, ex + 8 * S, gy - 2.5 * S, {stroke: C.ssi, "stroke-width": 7, "stroke-linecap": "round", hidden: true});
      s.elRing = K.circle(g, ex, gy - 2.5 * S, 20, {fill: "none", stroke: C.ssi, "stroke-width": 4, hidden: true, in: "pop"});
      // plan: the mat edge runs through the edge nodes (nodes at -32 ... 32 ft)
      const P = K.g(g, {hidden: true}); s.plan = P;
      const cx = E1.px, cy = E1.py, dx = 8 * S, x0 = cx - 4 * dx, y0 = cy - 4 * dx;
      K.text(P, cx, 318, "plan", {cls: "t-small", anchor: "middle", size: 26});
      K.rect(P, x0, y0, 8 * dx, 8 * dx, {fill: "#3b4a5e", stroke: C.concrete, "stroke-width": 3});
      for (let i = 0; i <= 8; i++) {
        K.line(P, x0 + i * dx, y0, x0 + i * dx, y0 + 8 * dx, {stroke: "rgba(185,196,208,.35)", "stroke-width": 2});
        K.line(P, x0, y0 + i * dx, x0 + 8 * dx, y0 + i * dx, {stroke: "rgba(185,196,208,.35)", "stroke-width": 2});
      }
      const pts = [];
      for (let j = 0; j <= 8; j++) for (let i = 0; i <= 8; i++) pts.push([x0 + i * dx, y0 + j * dx]);
      s.matNodes = K.nodes(P, pts, {r: 5.5, color: C.concrete});
      K.dim(P, x0, gy + 72, x0 + 8 * dx, gy + 72, "64 ft");
      K.dim(P, x0, y0 - 22, x0 + dx, y0 - 22, "8 ft");
      // the stick base, node 41 (the centre of the 9 x 9 grid)
      K.circle(P, cx, cy, 11, {fill: C.steel, stroke: "#0a111d", "stroke-width": 3});
      s.n41 = K.g(g, {hidden: true, in: "pop"});
      K.circle(s.n41, cx, cy, 20, {fill: "none", stroke: C.ssi, "stroke-width": 4});
      halo(K.text(s.n41, cx + 26, cy - 24, "node 41", {cls: "t-label", size: 26, color: "var(--ssi)"}));
      // the drilling rotation (about the vertical) at node 41
      s.drill = K.g(g, {hidden: true});
      const r = 30, a0 = -150 * Math.PI / 180, a1 = 150 * Math.PI / 180;
      const ex0 = cx + r * Math.cos(a0), ey0 = cy + r * Math.sin(a0), ex1 = cx + r * Math.cos(a1), ey1 = cy + r * Math.sin(a1);
      s.arc = K.path(s.drill, `M${ex0},${ey0}A${r},${r} 0 1 1 ${ex1},${ey1}`, {stroke: C.bad, "stroke-width": 6, fill: "none", "stroke-linecap": "round", draw: true});
      const tx = -Math.sin(a1), ty = Math.cos(a1);
      K.path(s.drill, `M${ex1 + tx * 15},${ey1 + ty * 15}L${ex1 - ty * 12},${ey1 + tx * 12}L${ex1 + ty * 12},${ey1 - tx * 12}Z`, {fill: C.bad});
      halo(K.text(s.drill, cx - 44, cy + 56, "θz", {cls: "t-label", anchor: "end", size: 30, color: "var(--bad)"}));
      // the rigid spider (group 3): node 41 to the 8 surrounding mat nodes (31, 32, 33, 40, 42, 49, 50, 51)
      s.spider = [[3, 3], [4, 3], [5, 3], [3, 4], [5, 4], [3, 5], [4, 5], [5, 5]].map(([i, j]) => {
        const lg = K.g(g, {hidden: true});
        K.line(lg, cx, cy, x0 + i * dx, y0 + j * dx, {stroke: C.ssi, "stroke-width": 7, "stroke-linecap": "round"});
        K.circle(lg, x0 + i * dx, y0 + j * dx, 8, {fill: C.ssi, stroke: "#0a111d", "stroke-width": 2});
        return lg;
      });
      s.spin = K.text(g, 1380, 820, "θz at node 41: the flat shells can't hold it", {cls: "t-label", anchor: "middle", color: "var(--bad)", size: 32, hidden: true});
      s.fine = K.text(g, 1380, 820, "✓ rigid spider: node 41 → 8 mat nodes", {cls: "t-label", anchor: "middle", color: "var(--good)", size: 32, hidden: true});
      s.rule = K.heading(s.root, "Read every warning, then decide. Never just silence it.", {x: 110, y: 890, w: 1700, size: "h3"});
    },
    beats: [
      {say: "[c]The command CHECK reads your input against the manual's list of errors and warnings. [z]Here: zero errors.",
        go(k) { k.show(k.s.head); }, c(k) { k.show(k.s.list.el); k.show(k.s.list.rows[0], {delay: 200}); },
        z(k) { k.show(k.s.list.rows.slice(1, 3), {stagger: 120}); }},
      {say: "[m]Then the model checks hunt for mistakes that run without complaint.",
        m(k) { k.show([k.s.elev, k.s.plan], {stagger: 150}); }},
      {say: "[a]A point where building meets soil, accidentally held fixed. [b]Or a joint that acts like a hinge.",
        a(k) { k.show(k.s.list.rows.slice(3, 5), {stagger: 120}); k.pulse(k.s.matNodes.dots, {amp: 0.6}); },
        b(k) { k.show(k.s.list.rows[5]); }},
      {say: "[h]Our example gets one warning: a possible hinge, where the stick stands on the mat.",
        h(k) { k.show(k.s.list.rows.slice(6, 8), {stagger: 120}); k.show([k.s.n41, k.s.elRing], {delay: 200}); }},
      {say: "[t]The flat mat elements can't hold a twist at that one point, so the stick could spin freely.",
        t(k) { k.show(k.s.drill); k.draw(k.s.arc, 800); k.show(k.s.spin, {delay: 200}); }},
      {say: "[s]But a stiff spider of beams carries the twist into the mat, [o]so here it's fine.",
        s(k) { k.hide([k.s.drill, k.s.spin, k.s.n41]); k.show(k.s.spider, {stagger: 90}); k.show(k.s.elSpider, {delay: 300}); },
        o(k) { k.show(k.s.list.rows[8]); k.show(k.s.fine); }},
      {say: "[r]The rule: read every warning, and decide. Never just silence it.",
        r(k) { k.show(k.s.rule); }},
    ],
  });

  // ---------------------------------------------------------------- 3. check 2: the slow-shake test
  // example 1 to scale (5.31 px/ft): 16 ft sand, 39 ft gravel, rock (cut with a fade); the mat and the 4 x 16 ft stick.
  // Motion: the computed transfer functions (amplitude and phase) at the first SSI frequency, 0.098 Hz.
  const F0 = FLOOR_TFU[0].f[0];
  const H0 = FLOOR_TFU.map((d) => [d.amp[0] * Math.cos(d.ph[0]), d.amp[0] * Math.sin(d.ph[0])]);
  const SLOW = {S: 340 / 64, cx: 450, gy: 590, U: 40, period: 3.4};
  scenes.push({
    id: "low-frequency", title: "Check 2: the slow-shake test",
    build(s) {
      const g = K.g(s.svg), S = SLOW.S, gy = SLOW.gy;
      s.head = K.heading(s.root, "Check 2 · The slow-shake test", {x: 110, y: 70, size: "h2"});
      s.soil = K.soil(g, {x: 130, y: gy, w: 640, labels: false, hidden: true,
        layers: [{h: 16 * S, kind: "sand"}, {h: 39 * S, kind: "gravel"}], hs: {h: 100, kind: "rock"}});
      // layer labels beyond the largest soil displacement (the column moves by up to U)
      [["sand, 16 ft", "Vs = 1,000 ft/s", gy + 8 * S], ["gravel, 39 ft", "Vs = 1,650 ft/s", gy + 35.5 * S], ["rock", "Vs = 3,300 ft/s", gy + 55 * S + 50]].forEach(([n, v, y]) => {
        K.text(s.soil.g, 130 + 640 + SLOW.U + 26, y - 4, n, {cls: "t-label", size: 26});
        K.text(s.soil.g, 130 + 640 + SLOW.U + 26, y + 26, v, {cls: "t-small", size: 22});
      });
      s.stick = K.stick(g, {x: SLOW.cx, y: gy - 5 * S, z: [16, 32, 48, 64].map((h) => h * S - 2.5 * S), matW: 64 * S, matH: 5 * S, r: 16, slabW: 0, hidden: true});
      // the ground motion, to the same scale: +-1 (unit control motion)
      s.gnd = K.g(g, {hidden: true});
      K.arrow(s.gnd, 210 - SLOW.U, gy - 26, 210 + SLOW.U, gy - 26, {color: C.wave, width: 4, head: 14, both: true});
      K.text(s.gnd, 210, gy - 50, "ground ±1", {cls: "t-label", anchor: "middle", color: "var(--wave)", size: 26});
      s.gDot = K.circle(s.gnd, 210, gy - 26, 9, {fill: C.wave, stroke: "#0a111d", "stroke-width": 2});
      s.note = K.text(g, 110, 182, `computed at ${F0.toFixed(3)} Hz · one scale for all · shown ${Math.round(1 / F0 / SLOW.period)} × faster`,
        {cls: "t-small", color: "var(--muted)", hidden: true});
      s.lf = K.text(g, 520, 290, "very slow shaking:", {cls: "t-label", color: "var(--ink)", size: 32, hidden: true});
      s.lf2 = K.text(g, 520, 334, "everything rides together", {cls: "t-label", color: "var(--ink)", size: 32, hidden: true});
      const [xs, ys] = clipXY(R_TFI.f, R_TFI.amp, 0, 2.0);
      s.p = K.plot(s.svg, {x: 1090, y: 190, w: 700, h: 320, xr: [0, 2], yr: [0, 2], xticks: [0, 0.5, 1, 1.5, 2], yticks: [0, 1, 2],
        xlabel: "frequency (Hz)", ylabel: "roof amplification", ylabelOffset: 62});
      s.one = s.p.hline(1, {color: C.good, dash: "12 9", width: 4, hidden: true});
      s.oneL = s.p.text(1.05, 1, "1 = moves like the ground", {cls: "t-label", color: "var(--good)", dy: 44, hidden: true});
      s.line = s.p.line(xs, ys, {color: C.ssi, width: 6, draw: true});
      s.dots = dots(s.p, R_TFU, {b: 2, r: 10});
      s.first = ring(s.p, R_TFU.f[0], R_TFU.amp[0], "", {color: "var(--good)", r: 22});
      s.stat = K.stat(s.root, {x: 1020, y: 600, w: 440, value: "0.103 %", label: "worst miss, at the roof", color: "var(--good)", vsize: 96});
      s.pass = K.pill(s.root, "✓ pass", {x: 1520, y: 640, color: "var(--good)", style: {fontSize: "36px"}});
      // the free-field column at F0 (one-dimensional wave solution of the example 1 site)
      const wv = PH.columnWaves(PH.COLUMN, F0);
      s.colU = (z) => PH.columnU(wv, PH.COLUMN, z);
      s.A = 0;
    },
    tick(s, t) {
      const w = TAU * t / SLOW.period, A = SLOW.U * s.A;
      const u = (H) => A * (H[0] * Math.cos(w) - H[1] * Math.sin(w));          // Re(H e^{iwt}) on screen
      const um = u(H0[0]);
      s.stick.set({sway: um, bend: H0.slice(1).map((H) => u(H) - um)});
      s.soil.shear((d) => u(s.colU(d / SLOW.S)));
      s.gDot.setAttribute("cx", 210 + u([1, 0]));
    },
    beats: [
      {say: "Check two is a thought experiment. [q]What if the ground shakes very, very slowly?",
        go(k) { k.show(k.s.head); },
        q(k) { k.show([k.s.soil.g, k.s.stick.g]); k.show([k.s.gnd, k.s.note], {delay: 300}); k.tween(k.s, {A: 1}, 1000, {delay: 300}); }},
      {say: "[a]Then everything rides along together, like a boat on a long, gentle swell.",
        a(k) { k.show([k.s.lf, k.s.lf2], {stagger: 150}); }},
      {say: "[p]So at the lowest frequency, every point should move exactly like the ground: an amplification of 1.",
        p(k) { k.show(k.s.p.g); k.draw(k.s.line, 1000, {delay: 200}); k.show([k.s.one, k.s.oneL], {delay: 300}); k.show(k.s.dots, {stagger: 100, delay: 400});
          k.show(k.s.first, {delay: 400}); }},
      {say: "[d]SASSI reports the worst miss. Here, at the roof, about 0.1 percent. [g]A pass.",
        d(k) { k.show(k.s.stat.el); }, g(k) { k.show(k.s.pass); }},
    ],
  });

  // ---------------------------------------------------------------- 4. computed points, filled-in curve
  scenes.push({
    id: "interpolation", title: "Computed points, filled-in curve",
    build(s) {
      s.head = K.heading(s.root, "Computed points, filled-in curve", {x: 110, y: 70, size: "h2"});
      s.p = K.plot(s.svg, {x: 240, y: 190, w: 1050, h: 490, xr: [0, 20], yr: [0, 14], xticks: [0, 4, 8, 12, 16, 20], yticks: [0, 4, 8, 12],
        xlabel: "frequency (Hz)", ylabel: "roof amplification", ylabelOffset: 62});
      s.first = s.p.line(R_TFI.f, R_TFI.amp, {color: C.ssi, width: 5, draw: true});
      s.tfi = s.p.line(R_TFI.f, R_TFI.amp, {color: C.ssi, width: 5, draw: true});
      s.dots = dots(s.p, R_TFU, {r: 11});
      s.probe = K.circle(s.p.g, 0, 0, 13, {fill: C.ink, stroke: C.ssi, "stroke-width": 4, hidden: true});
      s.n22 = K.stat(s.root, {x: 1390, y: 200, w: 420, value: "0", label: "computed frequencies", color: "var(--ssi)", vsize: 120});
      const dot = `<span style="display:inline-block;width:22px;height:22px;border:5px solid var(--ssi);border-radius:50%;vertical-align:-2px;margin-right:18px"></span>`;
      const bar = `<span style="display:inline-block;width:44px;height:6px;background:var(--ssi);border-radius:3px;vertical-align:9px;margin-right:18px"></span>`;
      s.leg1 = K.h("div", {x: 1400, y: 470, w: 430, html: `${dot}computed: <b style="color:var(--good)">a result</b>`, size: 34, in: "left"}, s.root);
      s.leg2 = K.h("div", {x: 1400, y: 540, w: 430, html: `${bar}filled in: <b style="color:var(--ssi)">a guess</b>`, size: 34, in: "left"}, s.root);
      s.P = 0;
    },
    tick(s, t) {
      if (!s.P) return;
      const f = 20 * ((t / 9) % 1);
      s.probe.setAttribute("cx", s.p.X(f)); s.probe.setAttribute("cy", s.p.Y(lerp(R_TFI.f, R_TFI.amp, f)));
    },
    beats: [
      {say: "[p]This curve shows how much the roof amplifies the ground shaking, at each frequency.",
        go(k) { k.show(k.s.head); }, p(k) { k.show(k.s.p.g); k.draw(k.s.first, 1000, {delay: 200}); }},
      {say: "Here's the surprise: [d]SASSI solves the full problem at only a few frequencies. [n]Here, 22.",
        d(k) { k.hide(k.s.first); k.show(k.s.dots, {stagger: 70, delay: 300}); },
        n(k) { k.show(k.s.n22.el); k.count(k.s.n22.v, R_TFU.f.length, {from: 0, ms: 1000}); }},
      {say: "[l]Between them, it fills in a smooth curve, like joining a few spot heights into a contour.",
        l(k) { k.draw(k.s.tfi, 1000); k.show(k.s.probe, {delay: 400}); k.tween(k.s, {P: 1}, 10, {delay: 400}); }},
      {say: "[c]Each computed point is a real result. [i]Everything in between is an educated guess.",
        c(k) { k.show(k.s.leg1); k.pulse(k.s.dots, {amp: 0.3}); }, i(k) { k.show(k.s.leg2); }},
    ],
  });

  // ---------------------------------------------------------------- 5. check 3: is every peak backed?
  // CRITFREQ,5,50 on the coarse run, computed here from its TFU and TFI (SV.K.PH.critfreq): the flagged peak,
  // the larger computed neighbour it is compared with, and the difference
  const CF = PH.critfreq(CO_TFU.f, CO_TFU.amp, CO_TFI.f, CO_TFI.amp, 5, 50, DF).filter((r) => r.flag)[0];
  const CF_M = CO_TFU.f.reduce((m, f, i) => (f <= CF.f ? i : m), 0);                         // the bracketing computed points m, m + 1
  const CF_ADD = AD_TFU.f.findIndex((f) => Math.abs(f - CF.f) < DF);                           // the added computed point
  const PK = argmax(R_TFI.amp);                                                                // the SSI peak of the full set (3.47 Hz)
  scenes.push({
    id: "critfreq", title: "Check 3: is every peak backed?",
    build(s) {
      s.head = K.heading(s.root, "Check 3 · Is every peak backed?", {x: 110, y: 70, size: "h2"});
      s.p = K.plot(s.svg, {x: 230, y: 190, w: 960, h: 480, xr: [0, 20], yr: [0, 14], xticks: [0, 4, 8, 12, 16, 20], yticks: [0, 4, 8, 12],
        xlabel: "frequency (Hz)", ylabel: "roof amplification", ylabelOffset: 62});
      const fpk = R_TFI.f[PK];
      s.pk = s.p.vline(fpk, {color: C.muted, hidden: true});
      s.pkL = s.p.text(fpk, 14, `peak: ${fpk.toFixed(2)} Hz`, {cls: "t-label", color: "var(--ink2)", dx: 16, dy: 34, hidden: true});
      s.dots = dots(s.p, CO_TFU, {color: C.violet, r: 11});
      s.n11 = K.pill(s.root, `${CO_TFU.f.length} frequencies`, {x: 820, y: 330, color: "var(--violet)", style: {fontSize: "30px"}});
      const top = Math.max(...CO_TFU.amp);
      s.ceil = s.p.hline(top, {color: C.violet, dash: "10 8", width: 3, hidden: true});
      s.ceilL = s.p.text(19.6, top, `highest computed: ${top.toFixed(2)}`, {cls: "t-label", anchor: "end", color: "var(--violet)", dy: -16, hidden: true});
      s.co = s.p.line(CO_TFI.f, CO_TFI.amp, {color: C.violet, width: 6, draw: true});
      s.top = ring(s.p, CF.f, CF.amp, CF.amp.toFixed(2), {color: "var(--ink)", r: 22, size: 36, dx: -30, dy: 12, anchor: "end"});
      // CRITFREQ: the peak against the larger of its two computed neighbours
      s.cmp = K.g(s.p.g, {hidden: true});
      const fa = CF.f + 1.6, xa = s.p.X(fa);
      K.line(s.cmp, s.p.X(CF.f) + 22, s.p.Y(CF.amp), xa + 12, s.p.Y(CF.amp), {stroke: C.bad, "stroke-width": 3, "stroke-dasharray": "8 7"});
      K.arrow(s.cmp, xa, s.p.Y(CF.ref), xa, s.p.Y(CF.amp), {color: C.bad, width: 5, head: 20, both: true});
      K.text(s.cmp, xa + 22, s.p.Y((CF.amp + CF.ref) / 2) + 12, `${(CF.amp / CF.ref).toFixed(1)} × (+${CF.d.toFixed(0)} %)`, {cls: "t-label", color: "var(--bad)", size: 36});
      s.added = s.p.dot(AD_TFU.f[CF_ADD], AD_TFU.amp[CF_ADD], {color: C.good, r: 14, hidden: true});
      s.code = K.code(s.root, ["CRITFREQ,5,50,00085TR_X,ADDC", "FREQ,1,@ADDC[1]"], {x: 1290, y: 190, w: 540, size: 27});
      s.isrs = K.stat(s.root, {x: 1290, y: 420, w: 540, value: "8.04", unit: "g", label: "roof spectrum peak,<br>resting on that guess", color: "var(--violet)", vsize: 100});
      s.flag = K.pill(s.root, `⚑ flagged: ${CF.amp.toFixed(1)} against ${CF.ref.toFixed(1)}`, {x: 1300, y: 660, color: "var(--bad)", style: {fontSize: "32px"}});
      s.okp = K.pill(s.root, `✓ confirmed: ${AD_TFU.amp[CF_ADD].toFixed(2)} computed`, {x: 1300, y: 660, color: "var(--good)", style: {fontSize: "32px"}});
      s.luck = K.card(s.root, {x: 230, y: 830, w: 1590, title: "Lucky this time", size: 32,
        body: "With several close modes, nothing guarantees the guess. Always back the peak."});
    },
    beats: [
      {say: "Check three. [q]What if the computed points miss a peak entirely?",
        go(k) { k.show(k.s.head); }, q(k) { k.show(k.s.p.g); k.show([k.s.pk, k.s.pkL], {delay: 300}); }},
      {say: "[c]We rerun the example with only 11 frequencies, [n]and none at the 3.47 hertz peak.",
        c(k) { k.show(k.s.dots, {stagger: 100}); k.show(k.s.n11, {delay: 400}); }, n(k) { k.pulse(k.s.pkL); }},
      {say: "[m]The highest computed point reaches only 4.7.",
        m(k) { k.show([k.s.ceil, k.s.ceilL]); }},
      {say: "[i]Yet the filled-in curve still shoots up to 13. [s]And the roof spectrum peak, 8.04 g, rests on it.",
        i(k) { k.hide([k.s.n11, k.s.pkL]); k.draw(k.s.co, 1000); k.show(k.s.top, {delay: 400}); }, s(k) { k.show(k.s.isrs.el); }},
      {say: "[k]That's the job of CRITFREQ. [e]It compares every big peak with its computed neighbours.",
        k(k) { k.show(k.s.code.el); k.type(k.s.code.lines[0], {cps: 30}); },
        e(k) { k.pulse([k.s.dots[CF_M], k.s.dots[CF_M + 1]], {amp: 0.4}); }},
      {say: "[f]This one is about three times higher than anything computed: flagged.",
        f(k) { k.show(k.s.cmp); k.show(k.s.flag, {delay: 300}); }},
      {say: "[a]Add the flagged frequency and run again. [g]Now a computed point says 13: confirmed.",
        a(k) { k.type(k.s.code.lines[1], {cps: 22}); },
        g(k) { k.hide([k.s.cmp, k.s.flag]); k.show(k.s.added); k.show(k.s.okp, {delay: 300}); }},
      {say: "[w]This time the guess was right. With several close modes, nothing guarantees it.",
        w(k) { k.show(k.s.luck); }},
    ],
  });

  // ---------------------------------------------------------------- 6. check 4: the cut-off
  const FN10 = MF_TFI.f[MF_TFI.f.length - 1];                // the last SSI frequency of the 10 Hz model
  // the spectrum frequency nearest 12 Hz (12.02 Hz) and the last one (100 Hz, the ZPA), as the lesson reads them
  const near = (d, f) => d.f.reduce((b, v, i) => (Math.abs(Math.log(v / f)) < Math.abs(Math.log(d.f[b] / f)) ? i : b), 0);
  const F12 = RS_M.f[near(RS_M, 12)];
  const at12 = [RS_M.sa[near(RS_M, 12)], RS_MF.sa[near(RS_MF, 12)]];
  const zpa = [RS_M.sa[near(RS_M, 100)], RS_MF.sa[near(RS_MF, 100)]];
  const fix = (v, d) => (v + 0.5e-6).toFixed(d);             // round half up (the stored values have four digits)
  scenes.push({
    id: "cutoff", title: "Check 4: the cut-off",
    build(s) {
      s.head = K.heading(s.root, "Check 4 · The cut-off", {x: 110, y: 70, size: "h2"});
      // part 1: the computed transfer function of the mat centre, cut-off 20 Hz, then 10 Hz
      const a = K.plot(s.svg, {x: 230, y: 190, w: 1020, h: 470, xr: [0, 25], yr: [0, 2], xticks: [0, 5, 10, 15, 20, 25], yticks: [0, 0.5, 1, 1.5, 2],
        yfmt: (v) => v.toFixed(1), xlabel: "frequency (Hz)", ylabel: "mat centre amplification", ylabelOffset: 70});
      s.a = a;
      const fN = M_TFU.f[M_TFU.f.length - 1], aN = M_TFU.amp[M_TFU.amp.length - 1];
      s.zone = K.rect(a.g, a.X(fN), a.box.y, a.X(25) - a.X(fN), a.box.h, {fill: "rgba(255,107,107,.14)", hidden: true});
      s.zero = a.text((fN + 25) / 2, 1.0, "set to zero", {cls: "t-label", anchor: "middle", color: "var(--bad)", size: 30, hidden: true});
      s.t20 = a.line(M_TFI.f, M_TFI.amp, {color: C.ref, width: 6, draw: true});
      s.tail20 = a.line([fN, fN, 25], [aN, 0, 0], {color: C.ref, width: 6, draw: true});
      s.mdots = dots(a, M_TFU, {color: C.ref, r: 8});
      s.last = ring(a, fN, aN, "last computed", {color: "var(--ink)", r: 18, size: 26, dx: 0, dy: -36, anchor: "middle"});
      s.cut20 = K.g(a.g, {hidden: true});
      K.line(s.cut20, a.X(fN), a.box.y - 6, a.X(fN), a.box.y + a.box.h, {stroke: C.ref, "stroke-width": 4, "stroke-dasharray": "12 8"});
      K.text(s.cut20, a.X(fN), a.box.y - 16, `cut-off ${fN.toFixed(0)} Hz`, {cls: "t-label", anchor: "middle", color: "var(--ink2)", size: 28});
      s.t10 = a.line(MF_TFI.f.concat([FN10, 25]), MF_TFI.amp.concat([0, 0]), {color: C.violet, width: 6, draw: true});
      s.cut10 = K.g(a.g, {hidden: true});
      K.line(s.cut10, a.X(FN10), a.box.y - 6, a.X(FN10), a.box.y + a.box.h, {stroke: C.violet, "stroke-width": 4, "stroke-dasharray": "12 8"});
      K.text(s.cut10, a.X(FN10), a.box.y - 16, `cut-off ${FN10.toFixed(0)} Hz`, {cls: "t-label", anchor: "middle", color: "var(--violet)", size: 28});
      s.lowpass = K.heading(s.root, "a low-pass filter", {x: 270, y: 760, size: "h2", color: "var(--bad)"});
      s.cutDim = K.pill(s.root, "same model: 20 Hz → 10 Hz", {x: 270, y: 870, color: "var(--violet)", style: {fontSize: "32px"}});
      // part 2: the mat spectrum, cut-off 20 Hz against 10 Hz
      s.q = K.plot(s.svg, {x: 230, y: 190, w: 1020, h: 470, xr: [1, 100], xlog: true, yr: [0, 1.6], xticks: [1, 2, 5, 10, 20, 50, 100],
        yticks: [0, 0.4, 0.8, 1.2, 1.6], yfmt: (v) => v.toFixed(1), xlabel: "equipment frequency (Hz)", ylabel: "mat spectrum, 5 % (g)", ylabelOffset: 70});
      s.band = K.rect(s.q.g, s.q.X(10), s.q.box.y, s.q.X(15) - s.q.X(10), s.q.box.h, {fill: "rgba(255,107,107,.18)", hidden: true});
      s.bandL = K.text(s.q.g, (s.q.X(10) + s.q.X(15)) / 2, s.q.box.y + 40, "10–15 Hz", {cls: "t-label", anchor: "middle", color: "var(--bad)", hidden: true});
      s.m20 = s.q.line(RS_M.f, RS_M.sa, {color: C.ref, width: 5, draw: true});
      s.m10 = s.q.line(RS_MF.f, RS_MF.sa, {color: C.violet, width: 6, draw: true});
      s.qleg = s.q.legend([{label: "cut-off 20 Hz", color: C.ref}, {label: "cut-off 10 Hz", color: C.violet}], {x: 930, y: 240, hidden: true, size: 28, dy: 46});
      s.mk = K.g(s.q.g, {hidden: true});
      K.arrow(s.mk, s.q.X(F12), s.q.Y(at12[0]) + 10, s.q.X(F12), s.q.Y(at12[1]) - 8, {color: C.bad, width: 5, head: 18});
      K.circle(s.mk, s.q.X(F12), s.q.Y(at12[0]), 10, {fill: C.ref, stroke: "#0a111d", "stroke-width": 3});
      K.circle(s.mk, s.q.X(F12), s.q.Y(at12[1]), 10, {fill: C.violet, stroke: "#0a111d", "stroke-width": 3});
      const pk20 = Math.max(...RS_R.sa), pk10 = Math.max(...RS_RF.sa);
      s.s1 = K.stat(s.root, {x: 1340, y: 190, w: 480, value: pk20.toFixed(2) === pk10.toFixed(2) ? pk20.toFixed(2) : `${pk20.toFixed(2)} · ${pk10.toFixed(2)}`,
        unit: "g", label: "roof peak, both runs", color: "var(--ssi)", vsize: 90});
      s.s2 = K.stat(s.root, {x: 1340, y: 400, w: 480, value: `${fix(zpa[0], 3)} · ${fix(zpa[1], 3)}`, unit: "g", label: "mat ZPA, 20 · 10 Hz", color: "var(--ink)", vsize: 62});
      s.s3 = K.stat(s.root, {x: 1340, y: 560, w: 480, value: `−${Math.round(100 * (1 - at12[1] / at12[0]))} %`,
        label: `at ${F12.toFixed(0)} Hz: ${fix(at12[0], 2)} → ${fix(at12[1], 2)} g`, color: "var(--bad)", vsize: 110});
      s.warn = K.card(s.root, {x: 230, y: 820, w: 960, kind: "warn", title: "Equipment at 10–15 Hz", size: 32, body: "designed for too little, and nothing warns you"});
      s.rule = K.card(s.root, {x: 1240, y: 820, w: 580, kind: "check", title: "The rule", size: 32, body: "cut-off above your band, with a margin"});
    },
    beats: [
      {say: "Check four. [q]What happens above the last frequency SASSI computes?",
        go(k) { k.show(k.s.head); },
        q(k) { k.show(k.s.a.g); k.draw(k.s.t20, 1000, {delay: 200}); k.show(k.s.mdots, {stagger: 40, delay: 200}); k.show(k.s.last, {delay: 400}); }},
      {say: "[z]Nothing. Above that cut-off, every response is set to zero, like a radio with the treble switched off.",
        z(k) { k.hide(k.s.last); k.show([k.s.cut20, k.s.zone]); k.draw(k.s.tail20, 600); k.show(k.s.zero, {delay: 300}); k.show(k.s.lowpass, {delay: 400}); }},
      {say: "[c]Let's stop the same model at 10 hertz instead of 20.",
        c(k) { k.show([k.s.cutDim, k.s.cut10]); k.draw(k.s.t10, 1000, {delay: 200}); }, hold: 600},
      {say: "[a]The big roof peak doesn't move: 8.03 g in both runs. [z]The ZPA barely changes either.",
        a(k) {
          const s = k.s;
          k.hide([s.a.g, s.lowpass, s.cutDim]);
          k.show(s.q.g, {delay: 200}); k.draw(s.m20, 1000, {delay: 200}); k.draw(s.m10, 1000, {delay: 400}); k.show(s.qleg, {delay: 400});
          k.show(s.s1.el, {delay: 300});
        },
        z(k) { k.show(k.s.s2.el); }},
      {say: "[r]But look at the mat, between 10 and 15 hertz.",
        r(k) { k.show([k.s.band, k.s.bandL]); }},
      {say: "[d]At 12 hertz, the spectrum drops from 0.82 g to 0.48: [p]41 percent less.",
        d(k) { k.show(k.s.mk); }, p(k) { k.show(k.s.s3.el); }},
      {say: "[w]Equipment there would be designed for too little, and neither the peak nor the ZPA warns you.",
        w(k) { k.show(k.s.warn); k.show(k.s.rule, {delay: 400}); }},
    ],
  });

  // ---------------------------------------------------------------- 7. the limits of the tool
  scenes.push({
    id: "limits", title: "What SASSI-EDU is, and is not",
    build(s) {
      const g = K.g(s.svg);
      s.head = K.heading(s.root, "What SASSI-EDU is, and is not", {x: 110, y: 70, size: "h2"});
      s.grid = K.g(g, {hidden: true});
      s.sq = [];
      for (let i = 0; i < 74; i++) {
        const c = i % 10, r = Math.floor(i / 10);
        s.sq.push(K.rect(s.grid, 150 + c * 50, 210 + r * 50, 40, 40, {rx: 6, fill: "#1a2a40", stroke: "#2a3b52", "stroke-width": 2}));
      }
      s.n74 = K.stat(s.root, {x: 150, y: 640, w: 490, value: "0", label: "verification problems,<br>all passing", color: "var(--good)", vsize: 120});
      s.cards = [
        K.card(s.root, {x: 780, y: 260, w: 1040, kind: "warn", title: "A teaching tool", size: 32, body: "Not qualified for licensing work: that needs qualified software."}),
        K.card(s.root, {x: 780, y: 520, w: 1040, title: "Reconstructed details", size: 32, body: "Where the manual is silent, results can differ from ACS SASSI."}),
      ];
      s.F = 0;
    },
    tick(s) {
      s.sq.forEach((e, i) => {
        const on = i < s.F;
        e.setAttribute("fill", on ? "rgba(99,230,164,.85)" : "#1a2a40");
        e.setAttribute("stroke", on ? "#63e6a4" : "#2a3b52");
      });
    },
    beats: [
      {say: "[v]SASSI-EDU is tested against 74 verification problems, from textbook formulas to published benchmarks. [p]All of them pass.",
        go(k) { k.show(k.s.head); k.show(k.s.grid, {delay: 200}); },
        v(k) { k.show(k.s.n74.el); k.count(k.s.n74.v, 74, {from: 0, ms: 1000, ease: "linear"}); k.tween(k.s, {F: 74}, 1000, {ease: "linear"}); },
        p(k) { k.pulse(k.s.n74.el); }},
      {say: "[q]But it's a teaching tool. It is not qualified for licensing work: that needs qualified software.",
        q(k) { k.show(k.s.cards[0]); }},
      {say: "[r]And where the manual is silent, it reconstructs, so details can differ from ACS SASSI.",
        r(k) { k.show(k.s.cards[1]); }},
    ],
  });

  // ---------------------------------------------------------------- 8. recap
  scenes.push({
    id: "recap", title: "Recap",
    build(s) {
      s.head = K.heading(s.root, "Recap", {x: 120, y: 80, size: "h2"});
      s.list = K.bullets(s.root, [
        {t: "Check the model before you run", sub: "and the slow shake: amplification 1"},
        {t: "Back every peak with a computed point", sub: "CRITFREQ flags the guesses"},
        {t: "Cut-off above your equipment band", sub: "with a margin"},
      ], {x: 120, y: 210, w: 960, num: true, size: 44});
      s.list.items.forEach((li) => { li.style.marginBottom = "52px"; });
      s.q = K.card(s.root, {x: 1150, y: 190, w: 680, kind: "check", title: "Check yourself", size: 32,
        body: "Equipment at 10–30 Hz. The analysis stopped at 15 Hz. Can you use the spectrum?"});
      s.a = K.card(s.root, {x: 1150, y: 530, w: 680, title: "Answer", size: 32,
        body: "Only well below 15 Hz. Extend the frequencies, or prove the band is unaffected."});
      s.next = K.pill(s.root, "Next · Lesson 11: an embedded building, the subtraction methods against FV →", {x: 120, y: 900, color: "var(--wave)"});
    },
    beats: [
      {say: "[a]Check the model before you run, and check that everything rides with the ground at the lowest frequency.",
        go(k) { k.show(k.s.head); }, a(k) { k.show(k.s.list.items[0]); }},
      {say: "[b]Back every peak that matters with a computed point. CRITFREQ finds the guesses.",
        b(k) { k.show(k.s.list.items[1]); }},
      {say: "[c]And put the cut-off above your equipment's frequencies, with a margin.",
        c(k) { k.show(k.s.list.items[2]); }},
      {say: "[q]Check yourself. Your equipment sits at 10 to 30 hertz; the analysis stopped at 15. Can you use the spectrum?",
        q(k) { k.show(k.s.q); }, gap: 1500},
      {say: "[a]Only well below 15 hertz. For the rest, extend the frequencies, or prove the band isn't affected.",
        a(k) { k.show(k.s.a); }},
      {say: "[n]Next, the final lesson: a real embedded building, and a shortcut that can invent a resonance.",
        n(k) { k.show(k.s.next); }},
    ],
  });

  SV.video({
    id: "10", n: 10, part: "Advanced", lesson: "10-good-practice",
    title: "Good practice, checks and the limits of SSI models",
    subtitle: "An SSI run can finish without an error and still be wrong: four quick checks, and the honest limits of the tool.",
    next: {href: "11.html", title: "An embedded building: the subtraction methods against FV"},
    pron: [[/\bCHECK\b/g, "check"], [/\bFIXEDINT\b/g, "fixed-int"], [/\bHINGED\b/g, "hinged"]],
    scenes,
  });
})();
