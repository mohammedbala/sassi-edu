/* SASSI-EDU explainer videos -- the motion-graphics kit (SV.K).
 *
 * Builders for the stage (1920 x 1080 logical pixels).  HTML builders take the scene root (s.root) and
 * place an absolutely positioned element; SVG builders take an SVG parent (s.svg or a group).  Every
 * builder returns its element(s) so that beats can reveal, move or draw them with the k helpers of
 * player.js.  HTML pieces start hidden (class sv-in) unless {visible: true}; set the entrance with
 * {in: "up" | "left" | "right" | "down" | "fade" | "pop"}.  SVG pieces are visible unless created with
 * {hidden: true} (or class "sv-in sv-fade"); strokes with {draw: true} start undrawn (class sv-draw).
 *
 * Colours: use the semantic CSS variables of player.css through K.C (K.C.wave, K.C.ssi ...).
 * Authoring guide: docs/internal/explainer_videos.md. */
"use strict";

(function () {
  const SV = window.SV;
  const NS = SV.NS;
  const K = (SV.K = {});
  const TAU = 2 * Math.PI;

  /** Semantic colours (the CSS variables of player.css, usable as SVG fill / stroke values). */
  K.C = {
    ink: "var(--ink)", ink2: "var(--ink2)", muted: "var(--muted)", faint: "var(--faint)", line: "var(--line)",
    wave: "var(--wave)", ssi: "var(--ssi)", ref: "var(--ref)", spring: "var(--spring)", dash: "var(--dash)", kin: "var(--kin)",
    good: "var(--good)", bad: "var(--bad)", violet: "var(--violet)",
    sand: "var(--sand)", gravel: "var(--gravel)", rock: "var(--rock)", deep: "var(--deep)",
    concrete: "var(--concrete)", concreteFill: "var(--concrete-fill)", steel: "var(--steel)",
  };
  /** The physics of the lesson figures (sassi/ui/static/figures.js), when the page loads it. */
  K.PH = (typeof SASSI !== "undefined" && SASSI.Figures && SASSI.Figures.physics) || null;

  // ================================================================== element builders
  function setAttrs(e, attrs) {
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v === undefined || v === null || v === false) continue;
      if (k === "text") e.textContent = v;
      else if (k === "html") e.innerHTML = v;
      else if (k === "cls") e.setAttribute("class", ((e.getAttribute("class") || "") + " " + v).trim());
      else if (k === "style" && typeof v === "object") Object.assign(e.style, v);
      else if (k === "hidden" || k === "draw" || k === "in" || k === "visible") continue;
      else e.setAttribute(k, v);
    }
    return e;
  }
  function entrance(e, o) {
    o = o || {};
    if (o.visible) return e;
    e.classList.add("sv-in");
    const how = o.in || "up";
    if (how !== "up") e.classList.add("sv-" + how);
    return e;
  }
  function svgFlags(e, o) {
    if (!o) return e;
    if (o.hidden) e.classList.add("sv-in", o.in && o.in !== "up" ? "sv-" + o.in : "sv-fade");
    if (o.draw) e.classList.add("sv-draw");
    return e;
  }
  /** SVG element. */
  K.s = function (tag, attrs, parent) {
    const e = document.createElementNS(NS, tag);
    setAttrs(e, attrs);
    svgFlags(e, attrs);
    if (parent) parent.appendChild(e);
    return e;
  };
  /** HTML element at stage coordinates: o.x, o.y, o.w, o.h (px) and the entrance o.in. */
  K.h = function (tag, o, parent) {
    o = o || {};
    const e = document.createElement(tag);
    setAttrs(e, {cls: o.cls, text: o.text, html: o.html, style: o.style});
    e.classList.add("sv-abs");
    if (o.x !== undefined) e.style.left = o.x + "px";
    if (o.y !== undefined) e.style.top = o.y + "px";
    if (o.w !== undefined) e.style.width = o.w + "px";
    if (o.h !== undefined) e.style.height = o.h + "px";
    if (o.align) e.style.textAlign = o.align;
    if (o.color) e.style.color = o.color;
    if (o.size) e.style.fontSize = o.size + "px";
    entrance(e, o);
    if (parent) parent.appendChild(e);
    return e;
  };
  K.g = (parent, attrs) => K.s("g", attrs || {}, parent);
  K.line = (p, x1, y1, x2, y2, a) => K.s("line", Object.assign({x1, y1, x2, y2}, a), p);
  K.rect = (p, x, y, w, h, a) => K.s("rect", Object.assign({x, y, width: w, height: h}, a), p);
  K.circle = (p, cx, cy, r, a) => K.s("circle", Object.assign({cx, cy, r}, a), p);
  K.path = (p, d, a) => K.s("path", Object.assign({d}, a), p);
  /** Polyline path through [[x, y], ...] (smooth: Catmull-Rom). */
  K.d = function (pts, smooth) {
    if (!pts.length) return "";
    if (!smooth || pts.length < 3) return "M" + pts.map((q) => q[0].toFixed(1) + "," + q[1].toFixed(1)).join("L");
    let d = `M${pts[0][0].toFixed(1)},${pts[0][1].toFixed(1)}`;
    for (let i = 0; i < pts.length - 1; i++) {
      const p0 = pts[i - 1] || pts[i], p1 = pts[i], p2 = pts[i + 1], p3 = pts[i + 2] || p2;
      const c1 = [p1[0] + (p2[0] - p0[0]) / 6, p1[1] + (p2[1] - p0[1]) / 6];
      const c2 = [p2[0] - (p3[0] - p1[0]) / 6, p2[1] - (p3[1] - p1[1]) / 6];
      d += `C${c1[0].toFixed(1)},${c1[1].toFixed(1)} ${c2[0].toFixed(1)},${c2[1].toFixed(1)} ${p2[0].toFixed(1)},${p2[1].toFixed(1)}`;
    }
    return d;
  };
  K.poly = (p, pts, a) => K.path(p, K.d(pts, a && a.smooth), a);
  /** SVG text; o.anchor start | middle | end, o.cls (t-label t-small t-big t-mono), o.size, o.color, o.weight. */
  K.text = function (p, x, y, str, o) {
    o = o || {};
    const e = K.s("text", {x, y, cls: o.cls || "t-label", "text-anchor": o.anchor || "start", hidden: o.hidden, in: o.in,
      "dominant-baseline": o.baseline || "auto"}, p);
    if (o.size) e.style.fontSize = o.size + "px";
    if (o.color) e.style.fill = o.color;
    if (o.weight) e.style.fontWeight = o.weight;
    if (o.mono) e.style.fontFamily = "var(--mono)";
    const lines = String(str).split("\n");
    if (lines.length === 1) e.textContent = str;
    else lines.forEach((ln, i) => K.s("tspan", {x, dy: i === 0 ? 0 : "1.25em", text: ln}, e));
    return e;
  };
  /** Arrow from (x1, y1) to (x2, y2): {g, line, head, update(x1, y1, x2, y2)}. */
  K.arrow = function (p, x1, y1, x2, y2, o) {
    o = o || {};
    const color = o.color || K.C.ink2, w = o.width || 4, hs = o.head || 16;
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const line = K.line(g, x1, y1, x2, y2, {stroke: color, "stroke-width": w, "stroke-linecap": "round", draw: o.draw, "stroke-dasharray": o.dash});
    const head = K.path(g, "", {fill: color});
    const a = {g, line, head};
    a.update = function (ax, ay, bx, by) {
      const L = Math.hypot(bx - ax, by - ay) || 1, ux = (bx - ax) / L, uy = (by - ay) / L;
      const ex = bx - ux * hs * 0.9, ey = by - uy * hs * 0.9;
      line.setAttribute("x1", ax); line.setAttribute("y1", ay); line.setAttribute("x2", ex); line.setAttribute("y2", ey);
      head.setAttribute("d", `M${bx},${by}L${bx - ux * hs - uy * hs * 0.55},${by - uy * hs + ux * hs * 0.55}L${bx - ux * hs + uy * hs * 0.55},${by - uy * hs - ux * hs * 0.55}Z`);
      if (o.both) {
        if (!a.head2) a.head2 = K.path(g, "", {fill: color});
        a.head2.setAttribute("d", `M${ax},${ay}L${ax + ux * hs - uy * hs * 0.55},${ay + uy * hs + ux * hs * 0.55}L${ax + ux * hs + uy * hs * 0.55},${ay + uy * hs - ux * hs * 0.55}Z`);
        line.setAttribute("x1", ax + ux * hs * 0.9); line.setAttribute("y1", ay + uy * hs * 0.9);
      }
    };
    a.update(x1, y1, x2, y2);
    return a;
  };
  /** Curved arrow (quadratic, bend = sideways offset of the control point in px). */
  K.curveArrow = function (p, x1, y1, x2, y2, o) {
    o = o || {};
    const color = o.color || K.C.ink2, w = o.width || 4, hs = o.head || 16, bend = o.bend === undefined ? 60 : o.bend;
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const mx = (x1 + x2) / 2, my = (y1 + y2) / 2, L = Math.hypot(x2 - x1, y2 - y1) || 1;
    const cx = mx - ((y2 - y1) / L) * bend, cy = my + ((x2 - x1) / L) * bend;
    const tx = x2 - cx, ty = y2 - cy, tl = Math.hypot(tx, ty) || 1, ux = tx / tl, uy = ty / tl;
    const ex = x2 - ux * hs * 0.9, ey = y2 - uy * hs * 0.9;
    const line = K.path(g, `M${x1},${y1}Q${cx},${cy} ${ex},${ey}`, {stroke: color, "stroke-width": w, fill: "none", "stroke-linecap": "round", draw: o.draw, "stroke-dasharray": o.dash});
    const head = K.path(g, `M${x2},${y2}L${x2 - ux * hs - uy * hs * 0.55},${y2 - uy * hs + ux * hs * 0.55}L${x2 - ux * hs + uy * hs * 0.55},${y2 - uy * hs - ux * hs * 0.55}Z`, {fill: color});
    return {g, line, head};
  };
  /** Engineering dimension line with end ticks and a label. */
  K.dim = function (p, x1, y1, x2, y2, label, o) {
    o = o || {};
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const c = o.color || K.C.muted;
    K.line(g, x1, y1, x2, y2, {stroke: c, "stroke-width": 2});
    const L = Math.hypot(x2 - x1, y2 - y1) || 1, nx = -(y2 - y1) / L * 10, ny = (x2 - x1) / L * 10;
    K.line(g, x1 - nx, y1 - ny, x1 + nx, y1 + ny, {stroke: c, "stroke-width": 2});
    K.line(g, x2 - nx, y2 - ny, x2 + nx, y2 + ny, {stroke: c, "stroke-width": 2});
    if (label) {
      const vertical = Math.abs(x2 - x1) < Math.abs(y2 - y1);
      const tx = (x1 + x2) / 2 + (vertical ? (o.side === "left" ? -14 : 14) : 0), ty = (y1 + y2) / 2 + (vertical ? 8 : -14);
      K.text(g, tx, ty, label, {cls: "t-small", anchor: vertical ? (o.side === "left" ? "end" : "start") : "middle", color: o.labelColor});
    }
    return g;
  };
  /** Fixed support: a line with hatching underneath (the fixed base of a structural model). */
  K.ground = function (p, x, y, w, o) {
    o = o || {};
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    K.line(g, x, y, x + w, y, {stroke: o.color || K.C.ref, "stroke-width": 5, "stroke-linecap": "round"});
    for (let xi = x + 6; xi < x + w; xi += 22) K.line(g, xi, y + 4, xi - 18, y + 24, {cls: "hatch"});
    return g;
  };

  // ================================================================== shared defs (patterns, gradients)
  K.defs = function () {
    if (K._defs) return K._defs;
    const svg = document.createElementNS(NS, "svg");
    svg.setAttribute("width", "0"); svg.setAttribute("height", "0");
    svg.style.position = "absolute";
    document.body.appendChild(svg);
    const defs = K.s("defs", {}, svg);
    const pat = (id, w, hgt, draw) => { const p = K.s("pattern", {id, width: w, height: hgt, patternUnits: "userSpaceOnUse"}, defs); draw(p); };
    pat("sv-pat-sand", 26, 22, (p) => {
      for (const [x, y] of [[4, 5], [16, 3], [10, 13], [22, 15], [3, 19]]) K.circle(p, x, y, 1.6, {fill: "rgba(255,240,210,.35)"});
    });
    pat("sv-pat-gravel", 44, 36, (p) => {
      for (const [x, y, r] of [[9, 9, 5], [30, 7, 4], [21, 24, 6], [39, 28, 3.5], [6, 30, 3]]) K.circle(p, x, y, r, {fill: "none", stroke: "rgba(255,235,200,.32)", "stroke-width": 1.6});
    });
    pat("sv-pat-rock", 30, 30, (p) => {
      K.path(p, "M0,30L30,0M-8,8L8,-8M22,38L38,22", {stroke: "rgba(255,255,255,.14)", "stroke-width": 2});
    });
    pat("sv-pat-deep", 30, 30, (p) => {
      K.path(p, "M0,30L30,0M-8,8L8,-8M22,38L38,22M0,0L30,30", {stroke: "rgba(255,255,255,.08)", "stroke-width": 2});
    });
    const lg = K.s("linearGradient", {id: "sv-grad-mass", x1: 0, y1: 0, x2: 0, y2: 1}, defs);
    K.s("stop", {offset: "0", "stop-color": "#f4f7fb"}, lg); K.s("stop", {offset: "1", "stop-color": "#9fb0c3"}, lg);
    const lg2 = K.s("linearGradient", {id: "sv-grad-fade", x1: 0, y1: 0, x2: 0, y2: 1}, defs);
    K.s("stop", {offset: "0", "stop-color": "#000", "stop-opacity": 0}, lg2); K.s("stop", {offset: "1", "stop-color": "#000", "stop-opacity": 0.55}, lg2);
    K._defs = defs;
    return defs;
  };

  // ================================================================== HTML components
  /** Lesson title card: {n, part, title, sub} -> {el, num, part, title, sub, rule} (reveal them in turn). */
  K.titleCard = function (root, o) {
    const el = K.h("div", {cls: "sv-title-card", visible: true}, root);
    const r = {el};
    r.num = K.h("div", {cls: "sv-num", text: String(o.n).padStart(2, "0"), in: "fade"}, el);
    r.num.classList.remove("sv-abs");
    r.part = K.h("div", {cls: "sv-part", text: o.part || "", in: "left"}, el);
    r.title = K.h("div", {cls: "sv-t", text: o.title}, el);
    r.sub = K.h("div", {cls: "sv-st", text: o.sub || ""}, el);
    r.rule = K.h("div", {cls: "sv-rule", in: "left"}, el);
    for (const k of ["part", "title", "sub", "rule"]) r[k].classList.remove("sv-abs");
    r.num.style.position = "absolute";
    return r;
  };
  K.kicker = (root, text, o) => K.h("div", Object.assign({cls: "sv-kicker", text}, o), root);
  /** Heading: o.size "h1" | "h2" | "h3" (default h2). */
  K.heading = (root, text, o) => K.h("div", Object.assign({cls: "sv-" + ((o && o.size) || "h2"), html: text}, o), root);
  K.para = (root, html, o) => K.h("div", Object.assign({cls: "sv-body", html}, o), root);
  K.label = (root, html, o) => K.h("div", Object.assign({cls: "sv-label", html}, o), root);
  /** Bullet list: items are strings (HTML) or {t, sub}; returns {el, items}. Reveal items one by one. */
  K.bullets = function (root, items, o) {
    o = o || {};
    const el = K.h("ul", {cls: "sv-bullets" + (o.num ? " num" : ""), x: o.x, y: o.y, w: o.w, visible: true}, root);
    const lis = items.map((it) => {
      const li = document.createElement("li");
      if (typeof it === "string") li.innerHTML = it;
      else { li.innerHTML = it.t; if (it.sub) { const s = document.createElement("span"); s.className = "sv-sub"; s.innerHTML = it.sub; li.appendChild(s); } }
      if (o.size) li.style.fontSize = o.size + "px";
      entrance(li, {in: o.in || "left", visible: o.visible});
      el.appendChild(li);
      return li;
    });
    if (o.color) lis.forEach((li) => li.style.setProperty("--b", o.color));
    if (o.dot) lis.forEach((li) => li.style.setProperty("--dot", o.dot));
    return {el, items: lis};
  };
  /** Card: {x, y, w, h, title, body (HTML), kind: "ansys" | "check" | "warn"} -> element. */
  K.card = function (root, o) {
    const el = K.h("div", Object.assign({}, o, {cls: "sv-card" + (o.kind ? " " + o.kind : "")}), root);
    if (o.title) K.h("div", {cls: "sv-card-t", text: o.title, visible: true}, el).classList.remove("sv-abs");
    const b = K.h("div", {cls: "sv-card-b", html: o.body || "", visible: true}, el);
    b.classList.remove("sv-abs");
    if (o.size) b.style.fontSize = o.size + "px";
    el.body = b;
    return el;
  };
  /** Big number: {x, y, w, value, unit, label, color} -> {el, v} (k.count(v, ...)). */
  K.stat = function (root, o) {
    const el = K.h("div", Object.assign({}, o, {cls: "sv-stat", in: o.in || "pop"}), root);
    const row = document.createElement("div");
    const v = document.createElement("span"); v.className = "sv-stat-v"; v.textContent = o.value === undefined ? "0" : o.value;
    if (o.color) v.style.color = o.color;
    if (o.vsize) v.style.fontSize = o.vsize + "px";
    row.appendChild(v);
    if (o.unit) { const u = document.createElement("span"); u.className = "sv-stat-u"; u.textContent = o.unit; row.appendChild(u); }
    el.appendChild(row);
    if (o.label) { const l = document.createElement("div"); l.className = "sv-stat-l"; l.innerHTML = o.label; el.appendChild(l); }
    return {el, v};
  };
  /** SASSI command block: lines (comment lines start with *); {x, y, w, title} -> {el, lines}.
   *  Lines start empty and visible only when typed (k.type(lines)) or shown (k.show). */
  K.code = function (root, lines, o) {
    o = o || {};
    const el = K.h("div", Object.assign({}, o, {cls: "sv-code"}), root);
    if (o.title !== false) { const t = document.createElement("div"); t.className = "sv-code-h"; t.textContent = o.title || "Command Entry"; el.appendChild(t); }
    if (o.size) el.style.fontSize = o.size + "px";
    const lns = lines.map((ln) => {
      const d = document.createElement("div");
      d.className = "ln" + (/^\s*\*/.test(ln) ? " c" : "");
      d.setAttribute("data-text", ln);
      d.setAttribute("data-cmd", "1");
      if (o.typed === false || o.static) {
        const i = ln.indexOf(",");
        if (/^\s*\*/.test(ln)) d.textContent = ln;
        else { const sp = document.createElement("span"); sp.className = "cmd"; sp.textContent = i < 0 ? ln : ln.slice(0, i); d.appendChild(sp); if (i >= 0) d.appendChild(document.createTextNode(ln.slice(i))); }
      }
      el.appendChild(d);
      return d;
    });
    return {el, lines: lns};
  };
  /** Equation (KaTeX; \htmlClass{name}{...} marks terms to highlight with k.cls). {x, y, size, box, color}. */
  K.eq = function (root, tex, o) {
    o = o || {};
    const el = K.h("div", Object.assign({}, o, {cls: "sv-eq" + (o.box ? " sv-eq-box" : "")}), root);
    if (o.size) el.style.fontSize = o.size + "px";
    K.tex(el, tex, o.display !== false);
    return el;
  };
  K.tex = function (el, tex, display) {
    if (typeof katex !== "undefined") {
      try { katex.render(tex, el, {displayMode: !!display, throwOnError: false, trust: true, strict: false}); return el; } catch (e) { /* fall through */ }
    }
    el.textContent = tex;
    return el;
  };
  /** Table: {x, y, head: [...], rows: [[...]], colW: [...], size} -> {el, rows: [tr]} (rows hidden). */
  K.table = function (root, o) {
    const el = K.h("table", Object.assign({}, o, {cls: "sv-table", visible: true}), root);
    if (o.size) el.style.fontSize = o.size + "px";
    if (o.head) {
      const tr = document.createElement("tr");
      o.head.forEach((c, i) => { const th = document.createElement("th"); th.innerHTML = c; if (o.colW && o.colW[i]) th.style.width = o.colW[i] + "px"; tr.appendChild(th); });
      el.appendChild(tr);
    }
    const rows = o.rows.map((r) => {
      const tr = document.createElement("tr");
      r.forEach((c) => { const td = document.createElement("td"); td.innerHTML = c; tr.appendChild(td); });
      if (!o.visible) entrance(tr, {in: "fade"});
      el.appendChild(tr);
      return tr;
    });
    return {el, rows};
  };
  /** Pill label (HTML) at x, y: {color}. */
  K.pill = (root, text, o) => K.h("div", Object.assign({cls: "sv-pill", html: text, in: (o && o.in) || "pop"}, o), root);
  /** The SASSI module chain as boxes with arrows: names [...] -> {el, boxes, arrows}. */
  K.chain = function (root, names, o) {
    o = o || {};
    const x = o.x || 120, y = o.y || 400, bw = o.bw || 200, bh = o.bh || 96, gap = o.gap || 56;
    const el = K.h("div", {x, y, w: names.length * (bw + gap), h: bh + 80, visible: true}, root);
    const boxes = [], arrows = [];
    names.forEach((n, i) => {
      const b = K.h("div", {x: i * (bw + gap), y: 0, w: bw, h: bh, in: o.in || "up",
        style: {border: "2px solid #33507a", borderRadius: "14px", background: "rgba(20,32,51,.9)", display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "center",
          fontFamily: "var(--mono)", fontWeight: "700", fontSize: (o.size || 30) + "px", color: "var(--ink)"}}, el);
      b.innerHTML = typeof n === "string" ? n : `${n.t}${n.sub ? `<span style="font: 500 19px var(--ui); color: var(--muted); margin-top: 4px">${n.sub}</span>` : ""}`;
      boxes.push(b);
      if (i < names.length - 1) {
        const a = K.h("div", {x: i * (bw + gap) + bw + 6, y: bh / 2 - 18, w: gap - 12, h: 36, in: "fade", html: `<svg viewBox="0 0 44 36" width="${gap - 12}" height="36"><path d="M2 18H36M26 8l12 10-12 10" stroke="#5f7da6" stroke-width="4" fill="none" stroke-linecap="round" stroke-linejoin="round"/></svg>`}, el);
        arrows.push(a);
      }
    });
    return {el, boxes, arrows};
  };

  // ================================================================== soil
  const LAYER_STYLE = {sand: ["var(--sand)", "sv-pat-sand"], gravel: ["var(--gravel)", "sv-pat-gravel"], rock: ["var(--rock)", "sv-pat-rock"], deep: ["var(--deep)", "sv-pat-deep"]};
  /** Layered soil block: {x, y (ground surface), w, layers: [{h (px), name, vs, kind: sand|gravel|rock|deep, color}],
   *  hs: {h (px), name, vs, kind} (half-space, fades out), labels: "right" | "left" | "inside" | false}
   *  -> {g, layers: [{g, poly, label}], ys, bottom, shear(u)}: shear(u) offsets the column by u(depth px) (px). */
  K.soil = function (p, o) {
    K.defs();
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const x = o.x, y = o.y, w = o.w;
    const list = o.layers.slice();
    if (o.hs) list.push(Object.assign({halfspace: true}, o.hs, {h: o.hs.h || 160}));
    const res = {g, layers: [], ys: [y], x, w};
    let top = y;
    const NSEG = 12;
    list.forEach((L, i) => {
      const kind = L.kind || ["sand", "gravel", "rock", "deep"][Math.min(i, 3)];
      const [col, pat] = LAYER_STYLE[kind] || LAYER_STYLE.sand;
      const lg = K.g(g, {});
      const poly = K.path(lg, "", {fill: L.color || col, "fill-opacity": L.halfspace ? 0.55 : 0.78});
      const tex = K.path(lg, "", {fill: `url(#${pat})`});
      const edge = K.path(lg, "", {stroke: "rgba(255,240,215,.35)", "stroke-width": 2, fill: "none"});
      const item = {g: lg, poly, tex, edge, top, h: L.h, L};
      if (L.halfspace) {
        const fade = K.rect(lg, x - 2, top + L.h * 0.35, w + 4, L.h * 0.65 + 2, {fill: "url(#sv-grad-fade)"});
        item.fade = fade;
      }
      if (o.labels !== false && (L.name || L.vs)) {
        const lx = o.labels === "left" ? x - 24 : o.labels === "inside" ? x + 24 : x + w + 24;
        const anchor = o.labels === "left" ? "end" : "start";
        const ly = top + Math.min(L.h, 120) / 2 + (L.halfspace ? 4 : 0);
        const lab = K.g(lg, {});
        K.text(lab, lx, ly - (L.vs ? 4 : -8), L.name || "", {cls: "t-label", anchor, size: o.labelSize || 26});
        if (L.vs) K.text(lab, lx, ly + 28, typeof L.vs === "number" ? `Vs = ${L.vs} m/s` : L.vs, {cls: "t-small", anchor, size: (o.labelSize || 26) - 4});
        item.label = lab;
      }
      res.layers.push(item);
      top += L.h;
      res.ys.push(top);
    });
    res.bottom = top;
    res.shear = function (u) {
      u = u || (() => 0);
      for (const it of res.layers) {
        const left = [], right = [];
        for (let j = 0; j <= NSEG; j++) {
          const yy = it.top + (it.h * j) / NSEG, du = u(yy - y);
          left.push([x + du, yy]); right.push([x + w + du, yy]);
        }
        const d = "M" + left.map((q) => q.join(",")).join("L") + "L" + right.reverse().map((q) => q.join(",")).join("L") + "Z";
        it.poly.setAttribute("d", d);
        it.tex.setAttribute("d", d);
        const du0 = u(it.top - y);
        it.edge.setAttribute("d", `M${x + du0},${it.top}L${x + w + du0},${it.top}`);
      }
    };
    res.shear();
    return res;
  };

  // ================================================================== structures
  /** Lumped-mass stick on a mat: {x (centre), y (top of the mat), z: [floor heights px above the mat],
   *  matW, matH, r (mass radius), slabW (floor bar width, 0 = none), color}
   *  -> {g, set({sway, rock, bend: [px per floor] | fn(z)}), nodes(): [[x, y] base..roof]}.
   *  sway: mat translation (px); rock: rotation (rad, + = clockwise on screen) about the mat centre. */
  K.stick = function (p, o) {
    K.defs();
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const x0 = o.x, y0 = o.y, zs = o.z, matW = o.matW || 260, matH = o.matH || 26, r = o.r || 20, slabW = o.slabW === undefined ? 120 : o.slabW;
    const mat = K.path(g, "", {cls: "mat"});
    const slabs = zs.map(() => (slabW ? K.path(g, "", {stroke: "rgba(223,230,238,.55)", "stroke-width": 6, "stroke-linecap": "round", fill: "none"}) : null));
    const stick = K.path(g, "", {cls: "struct", style: o.color ? {stroke: o.color} : null});
    const masses = zs.map(() => K.circle(g, 0, 0, r, {cls: "mass", fill: "url(#sv-grad-mass)"}));
    const obj = {g, mat, stick, masses, slabs, state: {sway: 0, rock: 0, bend: null}};
    const pt = (z, ub) => {
      // rotate about the mat centre (x0 + sway, y0 + matH / 2)
      const st = obj.state, cx = x0 + st.sway, cy = y0 + matH / 2;
      const lx = ub, ly = -z - matH / 2;
      const c = Math.cos(st.rock), s = Math.sin(st.rock);
      return [cx + lx * c - ly * s, cy + lx * s + ly * c];
    };
    obj.nodes = function () {
      const st = obj.state;
      const bend = (i) => (typeof st.bend === "function" ? st.bend(zs[i]) : st.bend ? st.bend[i] || 0 : 0);
      return [pt(0, 0)].concat(zs.map((z, i) => pt(z, bend(i))));
    };
    obj.set = function (st) {
      Object.assign(obj.state, st || {});
      const n = obj.nodes();
      stick.setAttribute("d", K.d(n, true));
      n.slice(1).forEach((q, i) => {
        masses[i].setAttribute("cx", q[0]); masses[i].setAttribute("cy", q[1]);
        if (slabs[i]) {
          const c = Math.cos(obj.state.rock), s = Math.sin(obj.state.rock);
          slabs[i].setAttribute("d", `M${q[0] - slabW / 2 * c},${q[1] - slabW / 2 * s}L${q[0] + slabW / 2 * c},${q[1] + slabW / 2 * s}`);
        }
      });
      const corners = [[-matW / 2, matH / 2], [matW / 2, matH / 2], [matW / 2, -matH / 2], [-matW / 2, -matH / 2]].map(([lx, ly]) => {
        const cx = x0 + obj.state.sway, cy = y0 + matH / 2, c = Math.cos(obj.state.rock), s = Math.sin(obj.state.rock);
        return [cx + lx * c - ly * s, cy + lx * s + ly * c];
      });
      mat.setAttribute("d", K.d(corners.concat([corners[0]])) + "Z");
      return obj;
    };
    obj.set({});
    return obj;
  };

  /** Shear-wall building (storeys of boxes, optional basement below the surface):
   *  {x (centre), y (ground surface), w, storeys, storeyH, basement: levels below grade, color}
   *  -> {g, set({sway, rock, drift: fn(z px above the foundation) -> px}), parts}. */
  K.building = function (p, o) {
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const x0 = o.x, w = o.w || 300, sh = o.storeyH || 90, ns = o.storeys || 3, nb = o.basement || 0;
    const base = o.y + nb * sh;                       // underside of the foundation (y px)
    const bh = o.matH || 22;
    const polys = [];
    const fill = o.fill || "var(--concrete-fill)", stroke = o.color || "var(--concrete)";
    const levels = nb + ns;
    for (let i = 0; i < levels; i++) {
      const pg = K.path(g, "", {fill: i < nb ? "rgba(40,53,70,.95)" : fill, stroke, "stroke-width": 3, "stroke-linejoin": "round"});
      polys.push(pg);
    }
    const slabs = [];
    for (let i = 0; i <= levels; i++) slabs.push(K.path(g, "", {stroke, "stroke-width": i === 0 ? 0 : 7, fill: "none", "stroke-linecap": "round"}));
    const mat = K.path(g, "", {cls: "mat"});
    const windows = [];
    if (o.windows !== false) for (let i = nb; i < levels; i++) for (let j = 0; j < 3; j++) windows.push({i, j, el: K.path(g, "", {fill: "rgba(76,195,255,.18)", stroke: "rgba(76,195,255,.35)", "stroke-width": 2})});
    const obj = {g, polys, slabs, mat, state: {sway: 0, rock: 0, drift: null}, base};
    const T = (lx, z) => {
      // lx: offset from the centre line, z: height above the underside of the foundation (px)
      const st = obj.state, cx = x0 + st.sway, cy = base;
      const u = st.drift ? st.drift(z) : 0;
      const X = lx + u, Y = -z;
      const c = Math.cos(st.rock), s = Math.sin(st.rock);
      return [cx + X * c - Y * s, cy + X * s + Y * c];
    };
    obj.T = T;
    obj.set = function (st) {
      Object.assign(obj.state, st || {});
      const seg = (z1, z2) => {
        const L = [], R = [];
        for (let k = 0; k <= 4; k++) { const z = z1 + ((z2 - z1) * k) / 4; L.push(T(-w / 2, z)); R.push(T(w / 2, z)); }
        return "M" + L.map((q) => q.join(",")).join("L") + "L" + R.reverse().map((q) => q.join(",")).join("L") + "Z";
      };
      for (let i = 0; i < levels; i++) polys[i].setAttribute("d", seg(i * sh, (i + 1) * sh));
      for (let i = 0; i <= levels; i++) { const a = T(-w / 2 - 6, i * sh), b = T(w / 2 + 6, i * sh); slabs[i].setAttribute("d", `M${a}L${b}`); }
      const m = [T(-w / 2 - 20, 0), T(w / 2 + 20, 0), T(w / 2 + 20, -bh), T(-w / 2 - 20, -bh)];
      mat.setAttribute("d", "M" + m.map((q) => q.join(",")).join("L") + "Z");
      for (const wd of windows) {
        const z1 = wd.i * sh + sh * 0.3, z2 = wd.i * sh + sh * 0.72, x1 = -w / 2 + (w * (wd.j + 0.22)) / 3, x2 = x1 + (w / 3) * 0.56;
        const q = [T(x1, z1), T(x2, z1), T(x2, z2), T(x1, z2)];
        wd.el.setAttribute("d", "M" + q.map((v) => v.join(",")).join("L") + "Z");
      }
      return obj;
    };
    obj.set({});
    return obj;
  };

  // ================================================================== springs, dashpots
  /** Zig-zag spring path from (x1, y1) to (x2, y2). */
  K.springD = function (x1, y1, x2, y2, o) {
    o = o || {};
    const n = o.coils || 6, amp = o.amp || 14;
    const L = Math.hypot(x2 - x1, y2 - y1) || 1, ux = (x2 - x1) / L, uy = (y2 - y1) / L, nx = -uy, ny = ux;
    const lead = Math.min(o.lead || 14, L * 0.2), body = L - 2 * lead;
    const pts = [[x1, y1], [x1 + ux * lead, y1 + uy * lead]];
    for (let i = 0; i < 2 * n; i++) {
      const t = lead + (body * (i + 0.5)) / (2 * n), s = i % 2 ? -amp : amp;
      pts.push([x1 + ux * t + nx * s, y1 + uy * t + ny * s]);
    }
    pts.push([x2 - ux * lead, y2 - uy * lead], [x2, y2]);
    return K.d(pts);
  };
  K.spring = function (p, x1, y1, x2, y2, o) {
    o = o || {};
    const el = K.path(p, K.springD(x1, y1, x2, y2, o), {cls: "spring", hidden: o.hidden, in: o.in, style: o.color ? {stroke: o.color} : null});
    return {el, g: el, update: (a, b, c, d) => el.setAttribute("d", K.springD(a, b, c, d, o))};
  };
  /** Dashpot from (x1, y1) (rod end) to (x2, y2) (cylinder end). */
  K.dashpot = function (p, x1, y1, x2, y2, o) {
    o = o || {};
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const st = o.color ? {stroke: o.color} : null;
    const rod = K.path(g, "", {cls: "dashpot", style: st});
    const cyl = K.path(g, "", {cls: "dashpot", style: st});
    const pist = K.path(g, "", {cls: "dashpot", style: st});
    const wd = o.w || 30;
    const update = (ax, ay, bx, by) => {
      const L = Math.hypot(bx - ax, by - ay) || 1, ux = (bx - ax) / L, uy = (by - ay) / L, nx = -uy, ny = ux;
      const cl = Math.min(L * 0.55, o.cyl || 70);                   // cylinder length from b
      const c0 = [bx - ux * cl, by - uy * cl];
      const hw = wd / 2;
      cyl.setAttribute("d", `M${c0[0] + nx * hw},${c0[1] + ny * hw}L${bx - ux * 6 + nx * hw},${by - uy * 6 + ny * hw}L${bx - ux * 6 - nx * hw},${by - uy * 6 - ny * hw}L${c0[0] - nx * hw},${c0[1] - ny * hw}M${bx - ux * 6},${by - uy * 6}L${bx},${by}`);
      const pp = [c0[0] + ux * cl * 0.45, c0[1] + uy * cl * 0.45];
      pist.setAttribute("d", `M${pp[0] + nx * (hw - 7)},${pp[1] + ny * (hw - 7)}L${pp[0] - nx * (hw - 7)},${pp[1] - ny * (hw - 7)}`);
      rod.setAttribute("d", `M${ax},${ay}L${pp[0]},${pp[1]}`);
    };
    update(x1, y1, x2, y2);
    return {g, update};
  };

  // ================================================================== waves
  /** Vertically propagating wave fronts between y0 (bottom) and y1 (top): {x, w, y0, y1, lambda (px),
   *  color, n} -> {g, set(t, speed px/s)}. */
  K.wavefronts = function (p, o) {
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const lam = o.lambda || 120, n = Math.ceil((o.y0 - o.y1) / lam) + 1;
    const lines = [];
    for (let i = 0; i < n; i++) lines.push(K.path(g, "", {stroke: o.color || K.C.wave, "stroke-width": o.width || 4, fill: "none", "stroke-linecap": "round"}));
    const obj = {g, set(t, speed) {
      const v = speed === undefined ? 120 : speed;
      lines.forEach((ln, i) => {
        const yy = o.y0 - (((i * lam + t * v) % (n * lam)));
        if (yy < o.y1 || yy > o.y0) { ln.setAttribute("d", ""); return; }
        const f = (o.y0 - yy) / (o.y0 - o.y1);
        const pts = [];
        for (let k = 0; k <= 30; k++) { const xx = o.x + (o.w * k) / 30; pts.push([xx, yy + Math.sin(k / 30 * TAU * 1.5 + t * 2) * (o.ripple || 0)]); }
        ln.setAttribute("d", K.d(pts));
        ln.style.opacity = String(Math.min(1, 4 * f, 4 * (1 - f)) * (o.opacity || 0.8));
      });
    }};
    obj.set(0);
    return obj;
  };
  /** Radiating ripples (radiation damping) from (cx, cy): {r0, r1, n, color, half: true (below cy only)}
   *  -> {g, set(t, period s)}. */
  K.ripples = function (p, o) {
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const n = o.n || 4, r0 = o.r0 || 40, r1 = o.r1 || 400;
    const arcs = [];
    for (let i = 0; i < n; i++) arcs.push(K.path(g, "", {stroke: o.color || K.C.dash, "stroke-width": o.width || 4, fill: "none"}));
    const obj = {g, amp: 1, set(t, period) {
      const T = period || 2.2;
      arcs.forEach((a, i) => {
        const ph = ((t / T + i / n) % 1);
        const r = r0 + (r1 - r0) * ph;
        const ry = r * (o.squash || 0.8);
        a.setAttribute("d", o.half === false ? `M${o.cx - r},${o.cy}A${r},${ry} 0 1 0 ${o.cx + r},${o.cy}A${r},${ry} 0 1 0 ${o.cx - r},${o.cy}` : `M${o.cx - r},${o.cy}A${r},${ry} 0 0 0 ${o.cx + r},${o.cy}`);
        a.style.opacity = String((1 - ph) * 0.9 * obj.amp);
      });
    }};
    obj.set(0);
    return obj;
  };
  /** A sine / signal trace in a box: {x, y, w, h, fn(t in 0..1) -> -1..1, color, n} -> path. */
  K.signal = function (p, o) {
    const n = o.n || 300, pts = [];
    for (let i = 0; i <= n; i++) pts.push([o.x + (o.w * i) / n, o.y + o.h / 2 - (o.h / 2) * o.fn(i / n)]);
    return K.path(p, K.d(pts), {stroke: o.color || K.C.wave, "stroke-width": o.width || 4, fill: "none", "stroke-linejoin": "round", draw: o.draw, hidden: o.hidden, cls: o.cls});
  };

  // ================================================================== plots
  /** Axes: {x, y, w, h (px box), xr: [a, b], yr: [a, b], xlog, ylog, xticks, yticks, xfmt, yfmt, xlabel, ylabel,
   *  grid (default true), hidden (default true: reveal with k.show(plot.g))}
   *  -> {g, X(v), Y(v), line(xs, ys, o), area(xs, ys, o), dot(x, y, o), vline(x, o), hline(y, o), text(x, y, s, o), legend(items, o)}. */
  K.plot = function (p, o) {
    const g = K.g(p, {hidden: o.hidden !== false, in: o.in || "fade"});
    const {x, y, w, h} = o;
    const lx = (v) => Math.log10(v);
    const X = (v) => x + (w * ((o.xlog ? lx(v) : v) - (o.xlog ? lx(o.xr[0]) : o.xr[0]))) / ((o.xlog ? lx(o.xr[1]) : o.xr[1]) - (o.xlog ? lx(o.xr[0]) : o.xr[0]));
    const Y = (v) => y + h - (h * ((o.ylog ? lx(v) : v) - (o.ylog ? lx(o.yr[0]) : o.yr[0]))) / ((o.ylog ? lx(o.yr[1]) : o.yr[1]) - (o.ylog ? lx(o.yr[0]) : o.yr[0]));
    const fmtv = (v) => (Math.abs(v) >= 1000 ? v.toLocaleString("en-US") : String(+v.toPrecision(4)));
    const axes = K.g(g, {});
    const grid = K.g(axes, {});
    if (o.xticks) for (const t of o.xticks) {
      if (o.grid !== false) K.line(grid, X(t), y, X(t), y + h, {cls: "grid"});
      K.text(axes, X(t), y + h + 34, o.xfmt ? o.xfmt(t) : fmtv(t), {cls: "tick", anchor: "middle"});
    }
    if (o.yticks) for (const t of o.yticks) {
      if (o.grid !== false) K.line(grid, x, Y(t), x + w, Y(t), {cls: "grid"});
      K.text(axes, x - 14, Y(t) + 8, o.yfmt ? o.yfmt(t) : fmtv(t), {cls: "tick", anchor: "end"});
    }
    K.path(axes, `M${x},${y}L${x},${y + h}L${x + w},${y + h}`, {cls: "ax"});
    if (o.xlabel) K.text(axes, x + w / 2, y + h + 76, o.xlabel, {cls: "axlab", anchor: "middle"});
    if (o.ylabel) {
      const t = K.text(axes, x - (o.ylabelOffset || 86), y + h / 2, o.ylabel, {cls: "axlab", anchor: "middle"});
      t.setAttribute("transform", `rotate(-90 ${x - (o.ylabelOffset || 86)} ${y + h / 2})`);
    }
    const clipId = "sv-clip-" + Math.random().toString(36).slice(2, 9);
    const cp = K.s("clipPath", {id: clipId}, g);
    K.rect(cp, x - 2, y - 30, w + 4, h + 32);
    const data = K.g(g, {"clip-path": o.clip === false ? null : `url(#${clipId})`});
    const P = {g, axes, data, X, Y, box: {x, y, w, h}};
    P.pts = (xs, ys) => xs.map((v, i) => [X(v), Y(Math.max(o.ylog ? 1e-9 : -Infinity, ys[i]))]);
    /** Curve: o.color, o.width, o.dash, o.draw (start undrawn), o.smooth, o.glow, o.hidden. */
    P.line = (xs, ys, lo) => {
      lo = lo || {};
      const e = K.path(data, K.d(P.pts(xs, ys), lo.smooth), {cls: "curve" + (lo.glow ? " glow" : ""), draw: lo.draw, hidden: lo.hidden, "stroke-dasharray": lo.dash,
        style: {stroke: lo.color || K.C.wave, strokeWidth: lo.width || 5, color: lo.color || K.C.wave, opacity: lo.opacity}});
      return e;
    };
    P.area = (xs, ys, lo) => {
      lo = lo || {};
      const pts = P.pts(xs, ys);
      const base = Y(lo.base !== undefined ? lo.base : (o.ylog ? o.yr[0] : Math.max(o.yr[0], 0)));
      const d = K.d(pts) + `L${pts[pts.length - 1][0]},${base}L${pts[0][0]},${base}Z`;
      return K.path(data, d, {fill: lo.color || K.C.wave, "fill-opacity": lo.opacity || 0.15, hidden: lo.hidden, in: lo.in});
    };
    P.dot = (vx, vy, lo) => { lo = lo || {}; return K.circle(data, X(vx), Y(vy), lo.r || 9, {fill: lo.color || K.C.ink, hidden: lo.hidden, in: lo.in || "pop", stroke: lo.stroke || "#0a111d", "stroke-width": 3}); };
    P.vline = (vx, lo) => { lo = lo || {}; return K.line(g, X(vx), y, X(vx), y + h, {stroke: lo.color || K.C.muted, "stroke-width": lo.width || 3, "stroke-dasharray": lo.dash || "10 8", hidden: lo.hidden, draw: lo.draw}); };
    P.hline = (vy, lo) => { lo = lo || {}; return K.line(g, x, Y(vy), x + w, Y(vy), {stroke: lo.color || K.C.muted, "stroke-width": lo.width || 3, "stroke-dasharray": lo.dash || "10 8", hidden: lo.hidden, draw: lo.draw}); };
    P.text = (vx, vy, s, lo) => K.text(g, X(vx) + ((lo && lo.dx) || 0), Y(vy) + ((lo && lo.dy) || 0), s, lo);
    /** Legend: items [{label, color, dash}] at o.x, o.y (px) -> group. */
    P.legend = (items, lo) => {
      lo = lo || {};
      const lg = K.g(g, {hidden: lo.hidden, in: lo.in});
      items.forEach((it, i) => {
        const yy = (lo.y !== undefined ? lo.y : y + 20) + i * (lo.dy || 42), xx = lo.x !== undefined ? lo.x : x + w - 330;
        K.line(lg, xx, yy, xx + 46, yy, {stroke: it.color, "stroke-width": 6, "stroke-dasharray": it.dash, "stroke-linecap": "round"});
        K.text(lg, xx + 62, yy + 9, it.label, {cls: "t-small", color: lo.textColor || "var(--ink2)", size: lo.size || 24});
      });
      return lg;
    };
    return P;
  };
  /** n log-spaced or linear values. */
  K.logspace = (a, b, n) => Array.from({length: n}, (_, i) => a * Math.pow(b / a, i / (n - 1)));
  K.linspace = (a, b, n) => Array.from({length: n}, (_, i) => a + ((b - a) * i) / (n - 1));
  /** |H(f)| of a damped single mode (hysteretic damping, participation P): a schematic transfer function. */
  K.modeTF = (f, f0, beta, P) => {
    const r = f / f0, re = 1 - r * r, im = 2 * beta;
    const base = Math.hypot(1, 2 * beta) / Math.hypot(re, im);
    return 1 + (P === undefined ? 1 : P) * (base - 1);
  };

  // ================================================================== meshes and nodes
  /** Dots at points [[x, y]] -> {g, dots}. */
  K.nodes = function (p, pts, o) {
    o = o || {};
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const dots = pts.map((q) => K.circle(g, q[0], q[1], o.r || 7, {fill: o.color || K.C.wave, stroke: o.stroke || "#0a111d", "stroke-width": o.sw || 2.5}));
    return {g, dots};
  };
  /** Rectangular mesh {x, y, cols, rows, dx, dy, color, nodes: true} -> {g, lines, pts}. */
  K.mesh = function (p, o) {
    const g = K.g(p, {hidden: o.hidden, in: o.in});
    const c = o.color || "rgba(185,196,208,.55)";
    const lines = K.g(g, {});
    for (let i = 0; i <= o.cols; i++) K.line(lines, o.x + i * o.dx, o.y, o.x + i * o.dx, o.y + o.rows * o.dy, {stroke: c, "stroke-width": o.width || 2});
    for (let j = 0; j <= o.rows; j++) K.line(lines, o.x, o.y + j * o.dy, o.x + o.cols * o.dx, o.y + j * o.dy, {stroke: c, "stroke-width": o.width || 2});
    const pts = [];
    for (let j = 0; j <= o.rows; j++) for (let i = 0; i <= o.cols; i++) pts.push([o.x + i * o.dx, o.y + j * o.dy, i, j]);
    return {g, lines, pts};
  };
})();
