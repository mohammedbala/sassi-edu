/* SASSI-EDU GUI -- plot tabs drawn with Plotly from the plot data of sassi.plotting.state
 * (GET /api/plot/<id>), the toolbars (requirements 5.9) and the Results browser.
 *
 * Families: 2d (SPECPLOT, THPLOT), layer (LAYERPLOT), soilprop (SOILPROPPLOT), 3d (MODELPLOT,
 * NODEPLOT, CUTPLOT) and anim (BUBBLEPLOT, VECTORPLOT, CONTOURPLOT, DEFORMPLOT).  Plot state lives
 * on the server; every toolbar button submits its command text (L17) and the 'update' event of the
 * plot state triggers the redraw. */
"use strict";

(function (S) {
  const el = S.el;
  const P = {};
  S.P = P;
  // no cloud buttons (plotly.js 3 "Share chart..." uploads to Plotly Cloud): the GUI is local and offline
  const PLOT_CONFIG = {displaylogo: false, responsive: true,
    modeBarButtonsToRemove: ["sendDataToCloud", "sendChartToCloud", "lasso2d", "select2d"], toImageButtonOptions: {format: "png"}};
  const THREE_D = ["MODELPLOT", "NODEPLOT", "CUTPLOT", "BUBBLEPLOT", "VECTORPLOT", "CONTOURPLOT", "DEFORMPLOT"];
  const ANIMS = ["BUBBLEPLOT", "VECTORPLOT", "CONTOURPLOT", "DEFORMPLOT"];
  const ELEMENT_PLOTS = ["MODELPLOT", "CUTPLOT", "CONTOURPLOT", "DEFORMPLOT"];
  /* setting command -> plot kinds it acts on (sassi.plotting.state.CAPABILITY) */
  const CAP = {
    ELECOLOR: ["MODELPLOT"], SHRINK: ["MODELPLOT"], WIREFRAME: ["MODELPLOT"],
    ELENUM: ELEMENT_PLOTS, GROUPNUM: ELEMENT_PLOTS, NODENUM: THREE_D,
    SHOWDOF: ["MODELPLOT", "NODEPLOT"], SHOWMASS: ["MODELPLOT", "NODEPLOT"],
    CNGVIEW: THREE_D, RSTVIEW: THREE_D, CNGCENTER: THREE_D, RSTCENTER: THREE_D, DEBUG: THREE_D, PAUSE: ANIMS,
  };
  P.CAP = CAP;

  // ------------------------------------------------------------------ helpers
  function plotRec(pid) {
    const ps = (S.state && S.state.plots && S.state.plots.plots) || [];
    return ps.find((p) => p.id === pid);
  }
  P.activePlot = () => S.state && S.state.plots && S.state.plots.active ? plotRec(S.state.plots.active) : null;
  P.activeKind = () => { const p = P.activePlot(); return p ? p.kind : null; };
  P.activeIs3D = () => THREE_D.includes(P.activeKind());
  function upsertPlotRec(rec) {
    if (!S.state.plots) S.state.plots = {plots: [], active: null};
    const ps = S.state.plots.plots;
    const i = ps.findIndex((p) => p.id === rec.id);
    if (i >= 0) ps[i] = rec; else ps.push(rec);
  }
  function jet(t) {
    const r = (x) => Math.max(0, Math.min(1, 1.5 - Math.abs(4 * t - x)));
    return `rgb(${Math.round(255 * r(3))},${Math.round(255 * r(2))},${Math.round(255 * r(1))})`;
  }
  const JET = Array.from({length: 17}, (_, k) => [k / 16, jet(k / 16)]);
  function rgbStr(c) { return `rgb(${Math.round(255 * c[0])},${Math.round(255 * c[1])},${Math.round(255 * c[2])})`; }
  const log10 = (v) => Math.log(v) / Math.LN10;

  // ------------------------------------------------------------------ plot tabs
  P.tabs = {};       // pid -> tab
  function tabId(pid) { return "plot" + pid; }
  function makePlotTab(rec) {
    const pid = rec.id;
    if (P.tabs[pid]) return P.tabs[pid];
    const host = el("div", {class: "plot-host"});
    const barInfo = el("span", {class: "grow"});
    const bar = el("div", {class: "pane-bar"}, barInfo,
      el("button", {class: "btn small", text: "Window Settings", onclick: () => S.D.windowSettings()}),
      el("button", {class: "btn small", text: "Export Image", onclick: () => S.D.exportImage()}));
    const body = el("div", {class: "pane-body"}, host);
    const pane = el("div", {class: "pane"}, bar, body);
    const t = S.addTab({id: tabId(pid), title: rec.caption + (rec.title ? ` - ${rec.title}` : ""), kind: "plot", pane, plotId: pid});
    t.host = host;
    t.bar = bar;
    t.barInfo = barInfo;
    t.rec = rec;
    t.dirty = true;
    t.onActivate = () => { if (t.dirty) scheduleRender(t); else resize(t); };
    t.onClose = () => { stopAnim(t); if (t.observer) t.observer.disconnect(); try { Plotly.purge(t.plotDiv || host); } catch (e) { /* ignore */ } };
    // the plot follows its tab's size (window resize, View > Toolbars, the status bar wrapping ...)
    if (window.ResizeObserver) {
      let pending = false;
      t.observer = new ResizeObserver(() => {
        if (pending || S.activeTab !== t) return;
        pending = true;
        requestAnimationFrame(() => { pending = false; resize(t); });
      });
      t.observer.observe(host);
    }
    P.tabs[pid] = t;
    return t;
  }
  function resize(t) {
    const div = t.plotDiv;
    if (!div || !div.data) return;
    try {
      Plotly.Plots.resize(div);
      // 3D: the orthographic fit depends on the width / height of the plot (sceneGeometry) -- unless the
      // learner has rotated or zoomed the view with the mouse: then the view is kept as it is
      const d = t.data;
      if (d && (d.family === "3d" || d.family === "anim") && d.camera && !t.userView) {
        t.ownRelayout = true;
        const lay = sceneLayout(d, null, d.family === "anim", div);
        t.fitAspectX = lay.scene.aspectratio.x;
        Promise.resolve(Plotly.relayout(div, {scene: lay.scene}))
          .finally(() => { t.ownRelayout = false; });
      }
    } catch (e) { /* hidden tab */ }
  }
  window.addEventListener("resize", () => { const t = S.activeTab; if (t && t.plotId) resize(t); });
  function scheduleRender(t) {
    t.dirty = true;
    if (S.activeTab !== t) return;            // rendered when the tab is brought to the front
    clearTimeout(t._timer);
    t._timer = setTimeout(() => render(t), 60);
  }
  async function render(t) {
    t.dirty = false;
    let d;
    try {
      d = await S.get(`/api/plot/${t.plotId}`);
    } catch (e) {
      t.host.innerHTML = "";
      t.host.appendChild(el("div", {class: "empty", text: `plot ${t.plotId}: ${e.message}`}));
      return;
    }
    t.data = d;
    S.setTabTitle(t.id, d.caption + (d.title ? ` - ${d.title}` : ""));
    t.barInfo.textContent = `plot ${d.id}  ${d.kind}` + (d.model !== null && d.model !== undefined ? `  model ${d.model}` : "");
    try {
      if (d.family === "2d") render2D(t, d);
      else if (d.family === "layer") renderLayer(t, d);
      else if (d.family === "soilprop") renderSoilProp(t, d);
      else if (d.family === "3d") render3D(t, d);
      else if (d.family === "anim") renderAnim(t, d);
    } catch (e) {
      console.error(e);
      t.host.innerHTML = "";
      t.host.appendChild(el("div", {class: "empty", text: `plot ${t.plotId} could not be drawn: ${e.message}`}));
    }
    S.updateToolbars();
  }
  P.render = render;
  function plotDiv(t, cls) {
    if (!t.plotDiv || !t.host.contains(t.plotDiv)) {
      try { Plotly.purge(t.plotDiv); } catch (e) { /* none */ }
      t.host.innerHTML = "";
      t.plotDiv = el("div", {class: cls || "plot-host"});
      t.host.appendChild(t.plotDiv);
    }
    return t.plotDiv;
  }

  // ------------------------------------------------------------------ 2D line plots
  function axisRange(lo, hi, log) {
    if (lo === null || hi === null || lo === undefined || hi === undefined) return undefined;
    if (log) { if (!(lo > 0 && hi > 0)) return undefined; return [log10(lo), log10(hi)]; }
    if (lo === hi) return [lo - 1, hi + 1];
    return [lo, hi];
  }
  function render2D(t, d) {
    const s = d.settings;
    const traces = d.lines.map((L) => ({
      type: "scattergl", mode: L.markers ? "lines+markers" : "lines", name: L.name || `line ${L.number}`,
      x: L.x, y: L.y, line: {color: L.color, width: 1.6, dash: L.dash === "dashdot" ? "dashdot" : L.dash},
      marker: {size: 4, color: L.color}, hovertemplate: `${S.esc(L.name || "line " + L.number)}<br>x=%{x}<br>y=%{y}<extra>line ${L.number}</extra>`,
    }));
    const ext = d.extent || [];
    const layout = {
      title: {text: S.esc(d.title || ""), font: {size: 14}},
      margin: {l: 70, r: 20, t: d.title ? 40 : 18, b: 50},
      showlegend: true, legend: {x: 1, xanchor: "right", y: 1, bgcolor: "rgba(255,255,255,.85)", bordercolor: "#999", borderwidth: 1},
      xaxis: {title: {text: S.esc(s.xtitle || (d.kind === "SPECPLOT" ? "Frequency (Hz)" : "Time (s)"))}, type: s.log_x ? "log" : "linear",
        showgrid: s.major_x, gridcolor: "#bbb", gridwidth: 1, minor: {showgrid: s.minor_x, gridcolor: "#e4e4e4"},
        range: axisRange(ext[0], ext[1], s.log_x), autorange: axisRange(ext[0], ext[1], s.log_x) ? false : true,
        showline: true, mirror: true, linecolor: "#000", zeroline: false, exponentformat: "e"},
      yaxis: {title: {text: S.esc(s.ytitle || "")}, type: s.log_y ? "log" : "linear",
        showgrid: s.major_y, gridcolor: "#bbb", minor: {showgrid: s.minor_y, gridcolor: "#e4e4e4"},
        range: axisRange(ext[2], ext[3], s.log_y), autorange: axisRange(ext[2], ext[3], s.log_y) ? false : true,
        showline: true, mirror: true, linecolor: "#000", zeroline: true, zerolinecolor: "#999", exponentformat: "e"},
      plot_bgcolor: "#fff", paper_bgcolor: "#fff", uirevision: JSON.stringify(s),
    };
    if (!traces.length) layout.annotations = [{text: "no defined line (lines deleted?)", showarrow: false}];
    Plotly.react(plotDiv(t), traces, layout, PLOT_CONFIG);
  }

  // ------------------------------------------------------------------ soil layer plot
  const GREYS = ["#e8e8e8", "#c8c8c8", "#a8a8a8", "#888888", "#686868", "#d8d8d8", "#b8b8b8", "#989898"];
  function renderLayer(t, d) {
    const tb = d.table;
    const show = d.show || {};
    try { Plotly.purge(t.plotDiv); } catch (e) { /* none */ }
    t.host.innerHTML = "";
    const left = el("div", {class: "left"});
    const right = el("div", {class: "right"});
    t.host.appendChild(el("div", {class: "plot-split"}, left, right));
    t.plotDiv = left;
    const rows = tb.layers;
    const total = rows.reduce((a, r) => a + r.thick, 0) || 1;
    const hsH = tb.halfspace ? 0.15 * total : 0;
    const top0 = rows.length ? rows[0].top : 0;
    const shapes = [], ann = [];
    rows.forEach((r, k) => {
      const y0 = r.top - top0, y1 = y0 + r.thick;
      shapes.push({type: "rect", x0: 0, x1: 1, y0, y1, fillcolor: GREYS[k % GREYS.length], line: {color: "#444", width: 1}});
      ann.push({x: -0.04, y: (y0 + y1) / 2, xanchor: "right", text: String(r.index), showarrow: false, font: {size: 12}});
      ann.push({x: 0.5, y: (y0 + y1) / 2, text: `L${r.layer}`, showarrow: false, font: {size: 10, color: "#222"}});
    });
    if (tb.halfspace) {
      const y0 = total, y1 = total + hsH;
      shapes.push({type: "rect", x0: 0, x1: 1, y0, y1, fillcolor: "#5c5c5c", line: {color: "#444", width: 1}});
      ann.push({x: 0.5, y: (y0 + y1) / 2, text: "Halfspace", showarrow: false, font: {size: 11, color: "#fff"}});
    }
    const layout = {margin: {l: 40, r: 10, t: 24, b: 10}, shapes, annotations: ann, showlegend: false,
      xaxis: {visible: false, range: [-0.25, 1.05]}, yaxis: {autorange: "reversed", title: {text: "depth"}, range: [total + hsH, 0], zeroline: false},
      title: {text: "Soil Layers", font: {size: 13}}, plot_bgcolor: "#fff"};
    Plotly.react(left, [], layout, PLOT_CONFIG);
    const cols = [["thick", "Thickness"], ["weight", "Unit Weight"], ["vp", "P-Wave Velocity"], ["vs", "S-Wave Velocity"],
      ["pdamp", "P-Wave Damping Ratio"], ["sdamp", "S-Wave Damping Ratio"]].filter(([k]) => show[k] !== false && (k !== "thick" || show.thick));
    const tbl = el("table", {class: "grid"});
    tbl.appendChild(el("thead", {}, el("tr", {}, el("th", {text: "Layer"}), el("th", {text: "L"}), ...cols.map(([, h]) => el("th", {text: h})))));
    const body = el("tbody");
    for (const r of rows) body.appendChild(el("tr", {}, el("td", {text: r.index}), el("td", {text: r.layer}), ...cols.map(([k]) => el("td", {text: S.fmt(r[k])}))));
    if (tb.halfspace) body.appendChild(el("tr", {}, el("td", {class: "l", text: "Halfspace"}), el("td", {text: tb.halfspace.layer}),
      ...cols.map(([k]) => el("td", {text: k === "thick" ? "" : S.fmt(tb.halfspace[k])}))));
    tbl.appendChild(body);
    right.appendChild(el("div", {class: "opt-note", text: `Source: ${tb.source}${tb.source === "L" ? " (no TOPL list: every L layer except the half-space)" : ""}; ` +
      `layers ${tb.start}..${tb.end} of ${tb.n_layers}` + (tb.omitted_above ? `; ${tb.omitted_above} layer(s) above not shown` : "") +
      (tb.omitted_below ? `; ${tb.omitted_below} layer(s) below not shown` : "") + (tb.missing && tb.missing.length ? `; TOPL layers without L: ${tb.missing.join(", ")}` : "")}));
    right.appendChild(tbl);
  }

  // ------------------------------------------------------------------ soil property plot
  function renderSoilProp(t, d) {
    const cv = d.curves, s = d.settings, show = (d.params && d.params.show) || {};
    const traces = [];
    if (show.damping !== false) traces.push({type: "scatter", mode: "lines+markers", name: "Damping Ratio/Shear Strain", x: cv.d_strain, y: cv.d,
      line: {color: "#e00000"}, marker: {size: 6}, yaxis: "y"});
    if (show.modulus !== false) traces.push({type: "scatter", mode: "lines+markers", name: "Shear Modulus/Shear Strain", x: cv.g_strain, y: cv.g,
      line: {color: "#00a000"}, marker: {size: 6}, yaxis: "y2"});
    const layout = {title: {text: S.esc(d.title || cv.name || ""), font: {size: 14}}, margin: {l: 70, r: 70, t: 40, b: 50},
      xaxis: {title: {text: S.esc(s.xtitle || "Shear Strain %")}, type: s.log_x ? "log" : "linear", showline: true, mirror: true, linecolor: "#000",
        showgrid: s.major_x, minor: {showgrid: s.minor_x}, exponentformat: "power"},
      yaxis: {title: {text: S.esc(s.ytitle || "Damping Ratio")}, color: "#c00000", type: s.log_y ? "log" : "linear", showline: true, linecolor: "#000",
        showgrid: s.major_y, minor: {showgrid: s.minor_y}},
      yaxis2: {title: {text: S.esc(s.ytitle2 || "Shear Modulus")}, overlaying: "y", side: "right", range: [0, 1.05], color: "#008000", showgrid: false,
        tickmode: "linear", tick0: 0, dtick: 0.2},
      legend: {x: 0.02, y: 0.02, bgcolor: "rgba(255,255,255,.85)", bordercolor: "#999", borderwidth: 1}, plot_bgcolor: "#fff",
      uirevision: JSON.stringify(s)};
    Plotly.react(plotDiv(t), traces, layout, PLOT_CONFIG);
  }

  // ------------------------------------------------------------------ 3D scenes
  /** Scene geometry of a 3D plot (computed by sceneLayout before cameraOf).
   *
   *  Plotly maps a data point p to normalised scene coordinates n = (p - mid) / range * aspectratio
   *  (the plot box spans +-aspectratio/2 about the box centre), and an orthographic view shows +-1
   *  normalised unit vertically and +-W/H horizontally whatever the eye distance (measured on the
   *  bundled plotly.js; its orthographic zoom works through the aspect ratio).  With explicit axis
   *  ranges and aspectratio = range * lam, a model length L is lam * L normalised units, so
   *  lam = min(1, W/H) / half shows exactly the screen window of the server's camera (GET /api/plot
   *  camera.half: fit of the rotated model about the rotation centre, divided by the zoom). */
  function sceneGeometry(d, fixed, div) {
    const b = (d.scene && d.scene.bbox) || [0, 1, 0, 1, 0, 1];
    const span = Math.max(b[1] - b[0], b[3] - b[2], b[5] - b[4], 1e-9);
    const m = (fixed ? 0.25 : 0.02) * span;          // animations: room for the deformed shape
    const r = [[b[0] - m, b[1] + m], [b[2] - m, b[3] + m], [b[4] - m, b[5] + m]];
    const W = (div && div.clientWidth) || 800, H = (div && div.clientHeight) || 600;
    const ex = d.camera.extent || [-1, 1, -1, 1];
    const half = d.camera.half || Math.max((ex[1] - ex[0]) / 2, (ex[3] - ex[2]) / 2, 1e-9);
    let lam = Math.min(1, W / H) / Math.max(half, 1e-12);
    // Fit the whole axis box (with the tick labels) rather than the model alone: the half extents of
    // the box corners about the rotation centre along the screen axes, 10 % margin, times the zoom.
    const B = d.camera.basis, c = d.camera.center;
    if (B && c) {
      let hx = 1e-12, hy = 1e-12;
      for (const x of r[0]) for (const y of r[1]) for (const z of r[2]) {
        const q = [x - c[0], y - c[1], z - c[2]];
        hx = Math.max(hx, Math.abs(q[0] * B[0][0] + q[1] * B[0][1] + q[2] * B[0][2]));
        hy = Math.max(hy, Math.abs(q[0] * B[1][0] + q[1] * B[1][1] + q[2] * B[1][2]));
      }
      const zoom = Number(d.view && d.view.zoom) > 0 ? Number(d.view.zoom) : 1;
      lam = 0.9 * zoom * Math.min((W / H) / hx, 1 / hy);
    }
    return {ranges: r, mid: r.map(([lo, hi]) => 0.5 * (lo + hi)), lam,
      aspect: {x: (r[0][1] - r[0][0]) * lam, y: (r[1][1] - r[1][0]) * lam, z: (r[2][1] - r[2][0]) * lam}};
  }
  /** Mouse rotation that keeps the model upright: when the view has no roll (its screen "up" is the
   *  projection of global Z) the camera's up vector is global Z and the scene rotates as a turntable
   *  about it; a rolled view (CNGVIEW about the view axis) or a view along Z keeps its own up vector and
   *  free (orbit) rotation. */
  function viewUp(B) {
    const dir = B[2], up = B[1];
    const zp = [-dir[2] * dir[0], -dir[2] * dir[1], 1 - dir[2] * dir[2]];
    const n = Math.hypot(zp[0], zp[1], zp[2]);
    if (n > 0.15 && (zp[0] * up[0] + zp[1] * up[1] + zp[2] * up[2]) / n > 0.999) {
      return {up: {x: 0, y: 0, z: 1}, dragmode: "turntable"};
    }
    return {up: {x: up[0], y: up[1], z: up[2]}, dragmode: "orbit"};
  }
  /** Plotly camera of a 3D plot: rotation (basis rows right / up / toward the viewer), the point
   *  looked at = rotation centre (CNGCENTER, else the box centre: camera.center) shifted by the pan
   *  (CNGVIEW <px>, <py> along the screen axes), and the zoom through sceneGeometry (D-UI-06). */
  function cameraOf(d) {
    const g = d._geo, B = d.camera.basis, c = d.camera.center, v = d.view;
    const px = Number(v.px) || 0, py = Number(v.py) || 0;
    const look = [0, 1, 2].map((i) => c[i] + px * B[0][i] + py * B[1][i]);
    const n = look.map((x, i) => (x - g.mid[i]) * g.lam);
    // the eye only gives the direction in an orthographic view; keep it outside the plot box
    const k = 2 + Math.hypot(g.aspect.x, g.aspect.y, g.aspect.z) + Math.hypot(n[0], n[1], n[2]);
    return {eye: {x: n[0] + B[2][0] * k, y: n[1] + B[2][1] * k, z: n[2] + B[2][2] * k}, up: viewUp(B).up,
      center: {x: n[0], y: n[1], z: n[2]}, projection: {type: "orthographic"}};
  }
  /** The commanded view of a 3D plot: the view values, the view counter (CNGVIEW / RSTVIEW / ... ) and the
   *  extents of what is drawn (a hide / show / display-volume change re-fits).  A new key replaces a view
   *  rotated or zoomed with the mouse. */
  function viewKey(d) {
    const v = (d && d.view) || {}, b = (d && d.scene && d.scene.bbox) || [];
    return JSON.stringify([v.rx, v.ry, v.rz, v.zoom, v.px, v.py, v.center, v.view_rev, b]);
  }
  /** sceneLayout keeping the learner's mouse view across re-renders of the same commanded view: a toggle
   *  such as NODENUM or the next animation frame re-renders the plot, and the rotation (the live camera)
   *  and the orthographic zoom (the aspect ratio) would otherwise go back to the commanded view.  Plotly's
   *  uirevision does not do it reliably (it compares camera values that differ in the last digits). */
  function layoutKeepingZoom(t, d, extra, fixed, div) {
    const lay = sceneLayout(d, extra, fixed, div);
    const a = lay.scene.aspectratio, fl = div && div._fullLayout && div._fullLayout.scene;
    const v = d.view || {};
    const key = viewKey(d);
    if (t.userView && t.viewKey === key && fl && fl.aspectratio && t.fitAspectX) {
      const f = fl.aspectratio.x / t.fitAspectX;
      if (isFinite(f) && f > 0) lay.scene.aspectratio = {x: a.x * f, y: a.y * f, z: a.z * f};
      try { if (fl._scene && fl._scene.getCamera) lay.scene.camera = fl._scene.getCamera(); } catch (e) { /* keep the commanded camera */ }
    }
    t.fitAspectX = a.x;
    return lay;
  }
  function sceneLayout(d, extra, fixed, div) {
    const v = d.view;
    d._geo = sceneGeometry(d, fixed, div);
    const g = d._geo;
    const key = viewKey(d);
    const ax = (t, r) => ({title: {text: t}, showbackground: false, gridcolor: "#e6e6e6", zerolinecolor: "#ccc", showspikes: false,
      range: r, autorange: false});
    const scene = {aspectmode: "manual", aspectratio: g.aspect, camera: cameraOf(d), xaxis: ax("X", g.ranges[0]), yaxis: ax("Y", g.ranges[1]),
      zaxis: ax("Z", g.ranges[2]), uirevision: key, dragmode: viewUp(d.camera.basis).dragmode};
    return Object.assign({
      margin: {l: 0, r: 0, t: d.title ? 30 : 6, b: 0}, title: {text: S.esc(d.title || ""), font: {size: 14}},
      showlegend: false, paper_bgcolor: "#fff", uirevision: key, scene,
    }, extra || {});
  }
  /** Face polygons (node index lists) of the scene; shrink toward the element centroid (5.7). */
  function faceList(sc, opts) {
    const out = [];
    const F = sc.faces, fe = sc.face_elem, fb = sc.face_boundary;
    for (let f = 0; f < F.length; f++) {
      if (!opts.all && !fb[f]) continue;
      if (opts.only && !opts.only(fe[f])) continue;
      const idx = F[f].filter((i) => i >= 0);
      out.push({f, e: fe[f], idx});
    }
    return out;
  }
  function polyCoords(sc, face, shrink, xyz) {
    const P0 = xyz || sc.xyz;
    const pts = face.idx.map((i) => P0[i]);
    if (!shrink) return pts;
    const c = sc.elem_centroid[face.e];
    return pts.map((p) => [c[0] + (1 - shrink) * (p[0] - c[0]), c[1] + (1 - shrink) * (p[1] - c[1]), c[2] + (1 - shrink) * (p[2] - c[2])]);
  }
  function meshTrace(sc, faces, opts) {
    const x = [], y = [], z = [], I = [], J = [], K = [], fc = [], tri2e = [];
    for (const face of faces) {
      const pts = polyCoords(sc, face, opts.shrink, opts.xyz);
      const b = x.length;
      for (const p of pts) { x.push(p[0]); y.push(p[1]); z.push(p[2]); }
      for (let j = 1; j < pts.length - 1; j++) {
        I.push(b); J.push(b + j); K.push(b + j + 1);
        fc.push(opts.color ? opts.color(face) : sc.elem_color[face.e]);
        tri2e.push(face.e);
      }
    }
    const tr = {type: "mesh3d", x, y, z, i: I, j: J, k: K, flatshading: true, hoverinfo: "skip", showscale: false,
      lighting: {ambient: 0.75, diffuse: 0.45, specular: 0.05, roughness: 0.9, fresnel: 0.05}, opacity: opts.opacity || 1, meta: {role: "mesh"}};
    if (opts.intensity) {
      tr.intensity = opts.intensity(x.length, faces, sc);
      tr.colorscale = JET; tr.cmin = 0; tr.cmax = 1; tr.intensitymode = "vertex";
      tr.showscale = true;
      tr.colorbar = opts.colorbar;
    } else tr.facecolor = fc;
    tr._tri2e = tri2e;
    return tr;
  }
  function outlineTrace(sc, faces, opts) {
    const x = [], y = [], z = [];
    for (const face of faces) {
      const pts = polyCoords(sc, face, opts.shrink, opts.xyz);
      for (const p of pts.concat([pts[0]])) { x.push(p[0]); y.push(p[1]); z.push(p[2]); }
      x.push(null); y.push(null); z.push(null);
    }
    return {type: "scatter3d", mode: "lines", x, y, z, line: {color: opts.color || "#202020", width: opts.width || 1.5}, hoverinfo: "skip"};
  }
  function lineElementTraces(sc, xyz) {
    const P0 = xyz || sc.xyz;
    const byColor = {};
    sc.edges.forEach((e, k) => {
      const c = sc.elem_color[sc.edge_elem[k]] || "#008000";
      const b = byColor[c] || (byColor[c] = {x: [], y: [], z: [], text: []});
      for (const i of e) { b.x.push(P0[i][0]); b.y.push(P0[i][1]); b.z.push(P0[i][2]); }
      b.x.push(null); b.y.push(null); b.z.push(null);
    });
    const tr = Object.entries(byColor).map(([c, b]) => ({type: "scatter3d", mode: "lines", x: b.x, y: b.y, z: b.z, line: {color: c, width: 5}, hoverinfo: "skip"}));
    if (sc.points.length) {
      tr.push({type: "scatter3d", mode: "markers", x: sc.points.map((i) => P0[i][0]), y: sc.points.map((i) => P0[i][1]), z: sc.points.map((i) => P0[i][2]),
        marker: {size: 5, symbol: "diamond", color: sc.point_elem.map((e) => sc.elem_color[e])}, hoverinfo: "skip"});
    }
    return tr;
  }
  function markerTrace(sc, idx, opts, xyz) {
    const P0 = xyz || sc.xyz;
    return {type: "scatter3d", mode: "markers", x: idx.map((i) => P0[i][0]), y: idx.map((i) => P0[i][1]), z: idx.map((i) => P0[i][2]),
      marker: Object.assign({size: 4}, opts.marker), name: opts.name || "",
      text: idx.map((i) => `Node ${sc.node_id[i]}`), hovertemplate: `%{text}${opts.hover || ""}<br>(%{x:.4g}, %{y:.4g}, %{z:.4g})<extra></extra>`,
      customdata: idx.map((i) => sc.node_id[i]), meta: {role: "nodes"}};
  }
  function labelTraces(d, sc, xyz) {
    const v = d.view, P0 = xyz || sc.xyz, out = [];
    const vis = sc.node_visible;
    if (v.node_labels) {
      const idx = sc.node_id.map((_, i) => i).filter((i) => vis[i]);
      out.push({type: "scatter3d", mode: "text", x: idx.map((i) => P0[i][0]), y: idx.map((i) => P0[i][1]), z: idx.map((i) => P0[i][2]),
        text: idx.map((i) => String(sc.node_id[i])), textfont: {size: 10, color: "#003c9c"}, textposition: "top right", hoverinfo: "skip"});
    }
    if (v.elem_labels || v.group_labels) {
      const C = sc.elem_centroid;
      out.push({type: "scatter3d", mode: "text", x: C.map((c) => c[0]), y: C.map((c) => c[1]), z: C.map((c) => c[2]),
        text: (v.elem_labels ? sc.elem_id : sc.elem_group).map(String), textfont: {size: 10, color: v.elem_labels ? "#8a3b00" : "#a00060"}, hoverinfo: "skip"});
    }
    return out;
  }
  function nodeMarkers(d, sc, nodePlot) {
    const out = [];
    const vis = sc.node_visible;
    const all = sc.node_id.map((_, i) => i).filter((i) => vis[i]);
    if (nodePlot) out.push(markerTrace(sc, all, {marker: {size: 3, color: "#000"}}));
    else out.push(markerTrace(sc, all, {marker: {size: 2, color: "rgba(0,0,0,0.25)"}}));
    if (sc.interaction.length) out.push(markerTrace(sc, sc.interaction, {marker: {size: nodePlot ? 4 : 3.5, color: "#e00000"}, hover: " (interaction)"}));
    if (sc.fixed.length) out.push(markerTrace(sc, sc.fixed, {marker: {size: 7, color: "rgba(0,0,0,0)", line: {color: "#00b000", width: 2}, symbol: "square-open"}, hover: " (fixed DOF)"}));
    if (sc.mass.length) {
      const dirs = [["#e00000", 0], ["#00b000", 1], ["#0000e0", 2]];
      for (const [c, j] of dirs) {
        const idx = sc.mass.filter((_, k) => sc.mass_dir[k][j]);
        if (idx.length) out.push(markerTrace(sc, idx, {marker: {size: 6, color: c, symbol: "diamond"}, hover: ` (mass ${"XYZ"[j]})`}));
      }
      const rot = sc.mass.filter((_, k) => sc.mass_rot[k] && !sc.mass_dir[k].some(Boolean));
      if (rot.length) out.push(markerTrace(sc, rot, {marker: {size: 6, color: "#9000c0", symbol: "diamond-open"}, hover: " (rotational mass)"}));
    }
    if (sc.selected.length) out.push(markerTrace(sc, sc.selected, {marker: {size: 8, color: "rgba(0,0,0,0)", line: {color: "#0060ff", width: 2}, symbol: "square-open"}, hover: " (selected)"}));
    return out;
  }
  function debugBox(t, d, extra) {
    t.host.querySelectorAll(".debug-box").forEach((b) => b.remove());
    if (!d.view.debug) return;
    const v = d.view;
    const txt = `view rX ${v.rx} rY ${v.ry} rZ ${v.rz}  pan ${v.px} ${v.py}  zoom ${v.zoom}\ncentre ${(d.camera.center || []).map((c) => S.fmt(c, 4)).join(", ")}` + (extra ? "\n" + extra : "");
    t.host.appendChild(el("div", {class: "debug-box", text: txt}));
  }
  /** A rotation / zoom with the mouse marks the view as the learner's: resize keeps it (see resize).  A
   *  new view from the server (CNGVIEW, RSTVIEW, CNGCENTER: another view key) clears the mark. */
  function trackView(t, div) {
    const v = (t.data && t.data.view) || {};
    const key = viewKey(t.data || {});
    if (t.viewKey !== key) { t.viewKey = key; t.userView = false; }
    if (!div || div._sassiView) return;
    div._sassiView = true;
    div.on("plotly_relayout", (ev) => {
      if (!t.ownRelayout && ev && Object.keys(ev).some((k) => k.startsWith("scene.camera") || k.startsWith("scene.aspectratio"))) t.userView = true;
    });
  }
  function attachClick(t) {
    const div = t.plotDiv;
    trackView(t, div);
    if (div._sassiClick) return;
    div._sassiClick = true;
    div.on("plotly_click", (ev) => {
      const pt = ev && ev.points && ev.points[0];
      if (!pt || !t.data) return;
      const tr = pt.data;
      if (tr.meta && tr.meta.role === "nodes" && pt.customdata !== undefined) {
        S.local("INFO", `Node ${pt.customdata}: (${S.fmt(pt.x)}, ${S.fmt(pt.y)}, ${S.fmt(pt.z)})`);
      } else if (tr.meta && tr.meta.role === "mesh" && tr._tri2e && pt.pointNumber !== undefined) {
        const sc = t.data.scene, e = tr._tri2e[pt.pointNumber];
        if (e !== undefined) S.local("INFO", `Element ${sc.elem_id[e]}, Group ${sc.elem_group[e]} (material ${sc.elem_mat[e]}, property ${sc.elem_prop[e]})`);
      }
    });
  }
  function render3D(t, d) {
    const sc = d.scene, v = d.view;
    const traces = [];
    if (!sc.node_id.length) {
      const div0 = plotDiv(t);
      Plotly.react(div0, [], sceneLayout(d, {annotations: [{text: "nothing to draw (all elements hidden?)", showarrow: false}]}, false, div0), PLOT_CONFIG);
      return;
    }
    if (d.kind === "MODELPLOT") {
      const shrink = v.shrink ? (d.shader.shrink || 0.06) : 0;
      const faces = faceList(sc, {all: !!shrink});
      if (!v.wireframe) traces.push(meshTrace(sc, faces, {shrink}));
      traces.push(outlineTrace(sc, v.wireframe ? faceList(sc, {all: true}) : faces, {shrink, color: v.wireframe ? "#303030" : "#202020", width: v.wireframe ? 1.2 : 1.5}));
      traces.push(...lineElementTraces(sc));
      traces.push(...nodeMarkers(d, sc, false));
    } else if (d.kind === "NODEPLOT") {
      traces.push(outlineTrace(sc, faceList(sc, {}), {color: "#c8c8c8", width: 1}));
      traces.push(...lineElementTraces(sc).map((tr) => Object.assign(tr, {line: {color: "#c8c8c8", width: 2}})));
      traces.push(...nodeMarkers(d, sc, true));
    } else if (d.kind === "CUTPLOT") {
      traces.push(outlineTrace(sc, faceList(sc, {}), {color: "#404040", width: 1}));
      const inCut = faceList(sc, {all: true, only: (e) => sc.elem_in_cut[e]});
      if (inCut.length) traces.push(meshTrace(sc, inCut, {color: () => "#e00000"}));
      traces.push(...lineElementTraces(sc));
    }
    traces.push(...labelTraces(d, sc));
    const div = plotDiv(t);
    Plotly.react(div, traces, layoutKeepingZoom(t, d, null, false, div), PLOT_CONFIG);
    attachClick(t);
    debugBox(t, d);
  }

  // ------------------------------------------------------------------ animations (5.8, basic)
  function animTraces(d, fr) {
    const sc = d.scene, v = d.view, traces = [];
    const kind = d.kind;
    if (kind === "DEFORMPLOT") {
      const faces = faceList(sc, {});
      if (v.show_undeformed) {        // the undeformed shape in grey: face outlines and line elements (sticks)
        traces.push(outlineTrace(sc, faces, {color: "#c8c8c8", width: 1}));
        traces.push(...lineElementTraces(sc).filter((tr) => tr.mode === "lines").map((tr) => Object.assign(tr, {line: {color: "#c8c8c8", width: 3}})));
      }
      traces.push(meshTrace(sc, faces, {xyz: fr.xyz, color: () => "#c00000", opacity: 0.85}));
      traces.push(outlineTrace(sc, faces, {xyz: fr.xyz, color: "#600000", width: 1}));
      traces.push(...lineElementTraces(sc, fr.xyz));
    } else if (kind === "CONTOURPLOT") {
      const faces = faceList(sc, {});
      const anim = d.animation;
      const lab = fr.colorbar || [];
      traces.push(meshTrace(sc, faces, {intensity: (n, fcs) => {
        const out = [];
        for (const face of fcs) for (const i of face.idx) out.push(fr.t[i] === null ? 0 : fr.t[i]);
        return out;
      }, colorbar: {title: {text: `col ${anim.col}`}, tickvals: [1, 0.75, 0.5, 0.25, 0], ticktext: lab, len: 0.8}}));
      traces.push(outlineTrace(sc, faces, {color: "#202020", width: 1}));
    } else {
      traces.push(outlineTrace(sc, faceList(sc, {}), {color: "#c8c8c8", width: 1}));
      traces.push(...lineElementTraces(sc).map((tr) => Object.assign(tr, {line: {color: "#c8c8c8", width: 2}})));
      const found = sc.node_id.map((_, i) => i).filter((i) => fr.found[i] && sc.node_visible[i]);
      if (kind === "BUBBLEPLOT") {
        traces.push({type: "scatter3d", mode: "markers", x: found.map((i) => sc.xyz[i][0]), y: found.map((i) => sc.xyz[i][1]), z: found.map((i) => sc.xyz[i][2]),
          marker: {size: found.map((i) => fr.size[i]), color: found.map((i) => rgbStr(fr.rgb[i])), line: {width: 0}},
          text: found.map((i) => `Node ${sc.node_id[i]}: ${S.fmt(fr.value[i])}`), hovertemplate: "%{text}<extra></extra>"});
        traces.push({type: "scatter3d", mode: "markers", x: [null], y: [null], z: [null], marker: {colorscale: JET, cmin: 0, cmax: 1, color: [0], showscale: true,
          colorbar: {tickvals: [1, 0.75, 0.5, 0.25, 0], ticktext: fr.colorbar || [], len: 0.8}}, hoverinfo: "skip"});
      } else if (kind === "VECTORPLOT") {
        const cols = ["#e00000", "#00b000", "#0000e0"];
        const dirs = v.direction === "ALL" ? [0, 1, 2] : [{X: 0, Y: 1, Z: 2}[v.direction] || 0];
        for (const dd of dirs) {
          const x = [], y = [], z = [];
          for (const i of found) {
            const p = sc.xyz[i], w = fr.vectors[i][dd];
            x.push(p[0], p[0] + w[0], null); y.push(p[1], p[1] + w[1], null); z.push(p[2], p[2] + w[2], null);
          }
          traces.push({type: "scatter3d", mode: "lines", x, y, z, line: {color: cols[dd], width: 4}, hoverinfo: "skip"});
        }
      }
    }
    traces.push(...labelTraces(d, sc, kind === "DEFORMPLOT" ? fr.xyz : null));
    return traces;
  }
  function stopAnim(t) { if (t.anim && t.anim.timer) { clearTimeout(t.anim.timer); t.anim.timer = null; } }
  function renderAnim(t, d) {
    stopAnim(t);
    const a = d.animation;
    const frames = a.frames || [];
    t.anim = t.anim || {cache: {}};
    t.anim.cache = {};
    t.anim.cache[d.frame.frame] = d.frame;
    t.anim.k = Math.max(0, frames.indexOf(a.current));
    // animation controls in the tab bar
    t.bar.querySelectorAll(".anim-bar").forEach((b) => b.remove());
    const label = el("span", {text: ""});
    const slider = el("input", {type: "range", min: 0, max: Math.max(frames.length - 1, 0), value: t.anim.k, "aria-label": "animation frame"});
    // Pause keeps the frame on the screen: PAUSE + WINDOWSETTINGS,FRAME,<k> (the replayed session and
    // CAPTUREPLOT show the same frame, L17); Start is PAUSE alone
    const togglePause = () => {
      if (!d.view.paused && frames.length) S.command(["PAUSE", `WINDOWSETTINGS,FRAME,${frames[t.anim.k]}`]);
      else S.command("PAUSE");
    };
    t.anim.togglePause = togglePause;
    const play = el("button", {class: "btn small", text: d.view.paused ? "Start" : "Pause", title: "PAUSE (toggle)", onclick: togglePause});
    const step = (s) => {
      if (!d.view.paused) return;
      const k = (t.anim.k + s + frames.length) % frames.length;
      S.command(`WINDOWSETTINGS,FRAME,${frames[k]}`);
    };
    const bar = el("span", {class: "anim-bar"}, play,
      el("button", {class: "btn small", text: "−", title: "previous frame (paused)", onclick: () => step(-1), disabled: !d.view.paused}),
      el("button", {class: "btn small", text: "+", title: "next frame (paused)", onclick: () => step(1), disabled: !d.view.paused}),
      slider, label);
    t.bar.insertBefore(bar, t.bar.children[1]);
    const show = async (k) => {
      const fk = frames[k];
      let fr = t.anim.cache[fk];
      if (!fr) {
        try { fr = await S.get(`/api/plot/${t.plotId}?frame=${fk}`); } catch (e) { S.status(e.message); return; }
        t.anim.cache[fk] = fr;
        const keys = Object.keys(t.anim.cache);
        if (keys.length > 400) delete t.anim.cache[keys[0]];
      }
      t.anim.k = k;
      slider.value = k;
      label.textContent = `frame ${fk} (${k + 1}/${frames.length}) ${fr.label || ""}`;
      const div = plotDiv(t);
      Plotly.react(div, animTraces(d, fr), layoutKeepingZoom(t, d, null, true, div), PLOT_CONFIG);
      trackView(t, div);
      debugBox(t, d, `frame ${fk}  ${fr.label || ""}  ${d.kind}  ${a.buffer_dir}`);
    };
    slider.addEventListener("input", () => show(Number(slider.value)));
    // a paused animation: the chosen frame is the plot's frame (WINDOWSETTINGS,FRAME)
    slider.addEventListener("change", () => { if (d.view.paused && frames.length) S.command(`WINDOWSETTINGS,FRAME,${frames[Number(slider.value)]}`); });
    const pause = Math.max(Number(d.params.frame_pause || 33), 60);
    const loop = async () => {
      if (!t.anim || d.view.paused || S.activeTab !== t || !frames.length) { t.anim.timer = null; return; }
      await show((t.anim.k + 1) % frames.length);
      t.anim.timer = setTimeout(loop, pause);
    };
    show(t.anim.k).then(() => { if (!d.view.paused) t.anim.timer = setTimeout(loop, pause); });
    const prevActivate = t.onActivate;
    t.onActivate = () => { prevActivate(); if (!t.anim.timer && !d.view.paused && t.data === d) t.anim.timer = setTimeout(loop, pause); };
  }

  // ------------------------------------------------------------------ plot events (sassi.plotting contract)
  P.onEvent = function (ev) {
    const kind = ev.event, pid = ev.plot_id, data = ev.data || {};
    if (kind === "open") {
      upsertPlotRec(data.plot);
      S.state.plots.active = pid;
      const t = makePlotTab(data.plot);
      if (!S.keepFocus) S.selectTab(t.id);
      scheduleRender(t);
    } else if (kind === "update") {
      upsertPlotRec(data.plot);
      const t = P.tabs[pid] || makePlotTab(data.plot);
      scheduleRender(t);
    } else if (kind === "close") {
      S.state.plots.plots = S.state.plots.plots.filter((p) => p.id !== pid);
      S.state.plots.active = data.active;
      const t = P.tabs[pid];
      delete P.tabs[pid];
      if (t) S.removeTab(t.id);
      if (data.active && P.tabs[data.active] && !S.keepFocus) S.selectTab(P.tabs[data.active].id);
    } else if (kind === "activate") {
      upsertPlotRec(data.plot);
      S.state.plots.active = pid;
      const t = P.tabs[pid] || makePlotTab(data.plot);
      if (S.activeTab !== t && !S.keepFocus) S.selectTab(t.id);
    } else if (kind === "lines") {
      if (S.D && S.D.onLinesChanged) S.D.onLinesChanged(data.numbers || []);
    } else if (kind === "capture") {
      S.status(`plot ${pid} saved to ${data.path}`);
    } else if (kind === "progress") {
      if (data.fraction >= 1) S.progress(null); else S.progress(data.fraction, `${data.command || ""} ${data.item || ""}`);
    } else if (kind === "dialog") {
      if (S.D && S.D.openNamedDialog) S.D.openNamedDialog(data.dialog, data.context || {});
    }
    S.updateToolbars();
  };
  P.restorePlots = function (ps) {
    if (!ps || !ps.plots) return;
    for (const rec of ps.plots) makePlotTab(rec);
    if (ps.active && P.tabs[ps.active]) S.selectTab(P.tabs[ps.active].id);
  };
  /** The model changed (a command ran): redraw the visible 3D plot of the changed models. */
  P.onModelChanged = function () {
    for (const t of Object.values(P.tabs)) {
      if (["3d", "anim"].includes(t.rec.family) || THREE_D.includes(t.rec.kind) || t.rec.kind === "LAYERPLOT" || t.rec.kind === "SOILPROPPLOT") {
        if (S.activeTab === t) scheduleRender(t); else t.dirty = true;
      }
    }
    P.refreshResults();
  };
  /** Keyboard controls of the 3D plots (spec 06 section 1.3): rotate 5 degrees per press. */
  P.onKey = function (ev) {
    const t = S.activeTab;
    if (!t || !t.plotId || !t.data || !THREE_D.includes(t.data.kind)) return;
    const v = t.data.view;
    const rot = {Insert: [5, 0, 0], Delete: [-5, 0, 0], Home: [0, 5, 0], End: [0, -5, 0], PageUp: [0, 0, 5], PageDown: [0, 0, -5]}[ev.key];
    if (rot) {
      ev.preventDefault();
      S.command(`CNGVIEW,${S.fmt(v.rx + rot[0])},${S.fmt(v.ry + rot[1])},${S.fmt(v.rz + rot[2])}`);
    } else if (ev.key === "Pause" && ANIMS.includes(t.data.kind)) {
      ev.preventDefault();
      P.togglePause();
    }
  };
  /** Toolbar Pause/Start, the Pause key and the tab's Pause button: see renderAnim. */
  P.togglePause = function () {
    const t = S.activeTab;
    if (t && t.anim && t.anim.togglePause) t.anim.togglePause();
    else S.command("PAUSE");
  };

  // ------------------------------------------------------------------ Results browser (extension)
  P.openResults = function () {
    let t = S.tab("results");
    if (t) { S.selectTab("results"); P.refreshResults(); return; }
    const list = el("div", {class: "pane-body"});
    const dirLabel = el("span", {class: "grow"});
    const filter = el("select", {}, ...[["all", "all files"], ["result", "results (TF, RS, histories)"], ["spectrum", "spectra / TF"], ["history", "time histories"], ["deck", "input decks"], ["text", "text / listings"]]
      .map(([v, l]) => el("option", {value: v, text: l})));
    const pane = el("div", {class: "pane"}, el("div", {class: "pane-bar"}, dirLabel, filter,
      el("button", {class: "btn small", text: "Refresh", onclick: () => P.refreshResults()})), list);
    t = S.addTab({id: "results", title: "Results", kind: "results", pane});
    t.list = list; t.dirLabel = dirLabel; t.filter = filter; t.dir = null;
    filter.addEventListener("change", () => P.refreshResults());
    S.selectTab("results");
    P.refreshResults();
  };
  P.refreshResults = async function (dir) {
    const t = S.tab("results");
    if (!t) return;
    if (dir !== undefined) t.dir = dir;
    let d;
    try { d = await S.get("/api/files" + (t.dir ? `?dir=${encodeURIComponent(t.dir)}` : "")); } catch (e) {
      t.list.innerHTML = ""; t.list.appendChild(el("div", {class: "empty", text: e.message})); return;
    }
    t.dir = d.dir;
    t.dirLabel.textContent = d.dir;
    const f = t.filter.value;
    const keep = (x) => f === "all" || (f === "result" ? ["spectrum", "history"].includes(x.kind) : x.kind === f);
    const tbl = el("table", {class: "grid results-list", style: {width: "100%"}});
    tbl.appendChild(el("thead", {}, el("tr", {}, el("th", {class: "l", text: "Name"}), el("th", {class: "l", text: "Type"}), el("th", {text: "Size"}),
      el("th", {text: "Modified"}), el("th", {text: ""}))));
    const body = el("tbody");
    if (d.parent) body.appendChild(el("tr", {ondblclick: () => P.refreshResults(d.parent)}, el("td", {class: "l name", text: ".."}), el("td"), el("td"), el("td"), el("td")));
    for (const sd of d.dirs) body.appendChild(el("tr", {ondblclick: () => P.refreshResults(sd.path)}, el("td", {class: "l name", text: sd.name + "/"}),
      el("td", {class: "l", text: "folder"}), el("td"), el("td"), el("td")));
    for (const x of d.files.filter(keep)) {
      const actions = el("td", {class: "l"});
      if (x.plot) actions.appendChild(el("button", {class: "btn small", text: "Plot", onclick: () => S.D.quickPlot(x.path, x.plot)}));
      if (["text", "deck", "spectrum", "history", "other"].includes(x.kind) || x.type === "PRE" || x.type === "ERR" || x.type === "OUT")
        actions.appendChild(el("button", {class: "btn small", text: "Open", onclick: () => S.openEditor(x.path)}));
      body.appendChild(el("tr", {ondblclick: () => x.plot ? S.D.quickPlot(x.path, x.plot) : (x.kind !== "binary" && x.kind !== "image" ? S.openEditor(x.path) : null)},
        el("td", {class: "l name", text: x.name}), el("td", {class: "l"}, el("span", {class: "chip " + x.kind, text: x.type || x.kind}), " ", x.label),
        el("td", {text: x.size.toLocaleString()}), el("td", {text: x.mtime}), actions));
    }
    tbl.appendChild(body);
    t.list.innerHTML = "";
    t.list.appendChild(tbl);
    if (!d.files.length && !d.dirs.length) t.list.appendChild(el("div", {class: "empty", text: "no files"}));
  };

  // ------------------------------------------------------------------ toolbars (5.9, spec 06 section 11)
  const ICON = {
    new: '<path d="M4 1.5h5.5l3 3V14.5H4z" fill="#fff" stroke="#555"/><path d="M9.5 1.5v3h3" fill="none" stroke="#555"/>',
    open: '<path d="M1.5 4.5h4l1.5 1.5h7v8h-12.5z" fill="#f2c94c" stroke="#8a6d1e"/><path d="M8 13V8M5.8 10.2 8 8l2.2 2.2" stroke="#1f5fbf" stroke-width="1.6" fill="none"/>',
    save: '<path d="M1.5 4.5h4l1.5 1.5h7v8h-12.5z" fill="#f2c94c" stroke="#8a6d1e"/><path d="M8 7.5v5M5.8 10.3 8 12.5l2.2-2.2" stroke="#1f5fbf" stroke-width="1.6" fill="none"/>',
    outpre: '<path d="M1.5 2.5h9l-3.5 4.5v4l-2 1.5V7z" fill="#9db9e6" stroke="#2c4f86"/><path d="M9 8h5.5v6.5H9z" fill="#fff" stroke="#555"/><path d="M10.3 10h3M10.3 12h3" stroke="#555"/>',
    outansys: '<path d="M1.5 2.5h9l-3.5 4.5v4l-2 1.5V7z" fill="#9db9e6" stroke="#2c4f86"/><text x="9" y="15" font-size="7.5" font-weight="700" fill="#b03000">A</text>',
    H: '<rect x="1.5" y="1.5" width="13" height="13" rx="2" fill="#fff" stroke="#555"/><text x="4.2" y="12.2" font-size="10" font-weight="700" fill="#2c4f86">H</text>',
    A: '<rect x="1.5" y="1.5" width="13" height="13" rx="2" fill="#fff" stroke="#555"/><text x="4.2" y="12.2" font-size="10" font-weight="700" fill="#b03000">A</text>',
    Sx: '<rect x="1.5" y="1.5" width="13" height="13" rx="2" fill="#fff" stroke="#555"/><text x="4.5" y="12.2" font-size="10" font-weight="700" fill="#555">S</text>',
    camera: '<rect x="1.5" y="4.5" width="13" height="9" rx="1.5" fill="#666" stroke="#333"/><circle cx="8" cy="9" r="2.7" fill="#cfe0f7" stroke="#fff"/><path d="M5 4.5l1-2h4l1 2" fill="#666" stroke="#333"/>',
    element: '<path d="M2 5l6-3 6 3v7l-6 3-6-3z" fill="#f6d55c" stroke="#7c6514"/><path d="M2 5l6 3 6-3M8 8v7" fill="none" stroke="#7c6514"/>',
    node: '<g fill="#333"><circle cx="3.5" cy="3.5" r="1.4"/><circle cx="8" cy="3.5" r="1.4"/><circle cx="12.5" cy="3.5" r="1.4"/><circle cx="3.5" cy="8" r="1.4"/><circle cx="8" cy="8" r="1.4" fill="#e00000"/><circle cx="12.5" cy="8" r="1.4"/><circle cx="3.5" cy="12.5" r="1.4"/><circle cx="8" cy="12.5" r="1.4"/><circle cx="12.5" cy="12.5" r="1.4"/></g>',
    cut: '<path d="M2 5l6-3 6 3v7l-6 3-6-3z" fill="#fff" stroke="#555"/><path d="M2 9.5l12-3" stroke="#e07000" stroke-width="2"/>',
    spectrum: '<path d="M1.5 14.5V1.5M1.5 14.5h13" stroke="#333" fill="none"/><path d="M2.5 13 5 9l2 1.5 2.5-6 2 4 2.5-2" fill="none" stroke="#e00000" stroke-width="1.6"/>',
    th: '<circle cx="8" cy="8" r="6.3" fill="#fff" stroke="#c62828" stroke-width="1.6"/><path d="M8 4v4.3l2.8 1.7" stroke="#c62828" stroke-width="1.5" fill="none"/>',
    layers: '<path d="M1.5 3h13v2.5h-13z" fill="#e0e0e0" stroke="#555"/><path d="M1.5 5.5h13v3h-13z" fill="#bdbdbd" stroke="#555"/><path d="M1.5 8.5h13v2.5h-13z" fill="#9e9e9e" stroke="#555"/><path d="M1.5 11h13v3h-13z" fill="#6d6d6d" stroke="#555"/>',
    soilprop: '<path d="M1.5 14.5V1.5M1.5 14.5h13" stroke="#333" fill="none"/><path d="M2.5 3c4 0 6 1 8 4s2.5 6 3.5 7" stroke="#00a000" fill="none" stroke-width="1.5"/><path d="M2.5 13c4 0 6-1 8-4s2.5-5 3.5-6" stroke="#e00000" fill="none" stroke-width="1.5"/>',
    frames: '<circle cx="8" cy="8" r="6.5" fill="#555"/><g fill="#ddd"><circle cx="8" cy="4.3" r="1.5"/><circle cx="8" cy="11.7" r="1.5"/><circle cx="4.3" cy="8" r="1.5"/><circle cx="11.7" cy="8" r="1.5"/></g>',
    bubble: '<circle cx="4.5" cy="10.5" r="3" fill="#e53935"/><circle cx="10.5" cy="10.5" r="2.3" fill="#43a047"/><circle cx="8" cy="5" r="2.8" fill="#1e88e5"/>',
    vector: '<circle cx="4" cy="12" r="2" fill="#222"/><path d="M4.8 11.2 13 3M13 3H8.5M13 3v4.5" stroke="#222" stroke-width="1.6" fill="none"/>',
    contour: '<circle cx="8" cy="8" r="6.5" fill="#1e5fd8"/><circle cx="8" cy="8" r="4.5" fill="#f6d55c"/><circle cx="8" cy="8" r="2.3" fill="#e53935"/>',
    deformed: '<path d="M4 15V8q0-5 7-6" stroke="#555" stroke-width="2.4" fill="none"/><circle cx="11.5" cy="2.5" r="1.7" fill="#e53935"/><path d="M4 15V2" stroke="#bbb" stroke-dasharray="2 1.5"/>',
    options: '<circle cx="8" cy="8" r="3" fill="#f6d55c" stroke="#555"/><path d="M8 1.5v2.5M8 12v2.5M1.5 8H4M12 8h2.5M3.4 3.4l1.8 1.8M10.8 10.8l1.8 1.8M3.4 12.6l1.8-1.8M10.8 5.2l1.8-1.8" stroke="#555" stroke-width="1.6"/>',
    afwrite: '<circle cx="8" cy="8" r="3" fill="#cfd8e3" stroke="#555"/><path d="M8 1.5v2.5M8 12v2.5M1.5 8H4M12 8h2.5M3.4 3.4l1.8 1.8M10.8 10.8l1.8 1.8M3.4 12.6l1.8-1.8M10.8 5.2l1.8-1.8" stroke="#555" stroke-width="1.6"/><path d="M9 11h6v4H9z" fill="#fff" stroke="#1f5fbf"/>',
    check: '<path d="M8 1.5l7 12.5H1z" fill="#ffd54f" stroke="#8a6d1e"/><path d="M8 6v4" stroke="#333" stroke-width="1.8"/><circle cx="8" cy="12" r="1" fill="#333"/>',
    results: '<path d="M2.5 1.5h11v13h-11z" fill="#fff" stroke="#555"/><path d="M4.5 4.5h7M4.5 7h7M4.5 9.5h7M4.5 12h4" stroke="#1f5fbf"/>',
    view: '<circle cx="6.5" cy="6.5" r="4.3" fill="#e3edf9" stroke="#333" stroke-width="1.4"/><path d="M9.7 9.7l4.5 4.5" stroke="#333" stroke-width="2"/>',
    resetview: '<circle cx="6.5" cy="6.5" r="4.3" fill="#e3edf9" stroke="#333" stroke-width="1.4"/><path d="M9.7 9.7l4.5 4.5" stroke="#333" stroke-width="2"/><path d="M4.5 4.5l4 4M8.5 4.5l-4 4" stroke="#1f5fbf" stroke-width="1.5"/>',
    center: '<circle cx="8" cy="8" r="5.2" fill="none" stroke="#333" stroke-width="1.3"/><circle cx="8" cy="8" r="2" fill="#e00000"/><path d="M8 0.5v4M8 11.5v4M0.5 8h4M11.5 8h4" stroke="#333"/>',
    resetcenter: '<circle cx="8" cy="8" r="5.2" fill="none" stroke="#333" stroke-width="1.3"/><path d="M8 0.5v4M8 11.5v4M0.5 8h4M11.5 8h4" stroke="#333"/><path d="M5.8 5.8l4.4 4.4M10.2 5.8l-4.4 4.4" stroke="#1f5fbf" stroke-width="1.5"/>',
    wire: '<path d="M2.5 5.5h7v7h-7z M6.5 2.5h7v7h-7z" fill="none" stroke="#333"/><path d="M2.5 5.5l4-3M9.5 5.5l4-3M9.5 12.5l4-3M2.5 12.5l4-3" stroke="#333"/>',
    shrink: '<rect x="1.5" y="1.5" width="5.5" height="5.5" fill="#e57373"/><rect x="9" y="1.5" width="5.5" height="5.5" fill="#64b5f6"/><rect x="1.5" y="9" width="5.5" height="5.5" fill="#64b5f6"/><rect x="9" y="9" width="5.5" height="5.5" fill="#e57373"/>',
    G: '<g fill-opacity=".8"><circle cx="3" cy="13" r="2" fill="#e53935"/><circle cx="8" cy="13.5" r="2" fill="#43a047"/><circle cx="13" cy="13" r="2" fill="#1e88e5"/></g><text x="4" y="10" font-size="10" font-weight="700" fill="#222">G</text>',
    M: '<g fill-opacity=".8"><circle cx="3" cy="13" r="2" fill="#e53935"/><circle cx="8" cy="13.5" r="2" fill="#43a047"/><circle cx="13" cy="13" r="2" fill="#1e88e5"/></g><text x="3.3" y="10" font-size="10" font-weight="700" fill="#222">M</text>',
    Pp: '<g fill-opacity=".8"><circle cx="3" cy="13" r="2" fill="#e53935"/><circle cx="8" cy="13.5" r="2" fill="#43a047"/><circle cx="13" cy="13" r="2" fill="#1e88e5"/></g><text x="4.3" y="10" font-size="10" font-weight="700" fill="#222">P</text>',
    nodenum: '<circle cx="4" cy="12" r="2" fill="#333"/><text x="6" y="8.5" font-size="8" font-weight="700" fill="#1f5fbf">N1</text>',
    elemnum: '<rect x="2" y="6" width="7" height="8" fill="#9db9e6" stroke="#2c4f86"/><text x="8.5" y="8" font-size="8" font-weight="700" fill="#e07000">E</text>',
    groupnum: '<rect x="2" y="6" width="7" height="8" fill="#9db9e6" stroke="#2c4f86"/><text x="8.5" y="8" font-size="8" font-weight="700" fill="#c62828">G</text>',
    showdof: '<path d="M8 3.5 3.5 10h9z" fill="#c8e6c9" stroke="#2e7d32"/><path d="M2.5 12.5h11" stroke="#2e7d32" stroke-width="1.5"/><circle cx="8" cy="3.5" r="1.6" fill="#333"/>',
    showmass: '<circle cx="6" cy="10" r="2" fill="#333"/><path d="M11 2l1 2.3 2.5.2-1.9 1.6.6 2.4-2.2-1.3-2.2 1.3.6-2.4L7.5 4.5 10 4.3z" fill="#9000c0"/>',
    pause: '<path d="M2 2.5v11l6-5.5z" fill="#2e7d32"/><path d="M10 3h1.8v10H10zM13 3h1.8v10H13z" fill="#c62828"/>',
    learn: '<path d="M8 2.2 0.8 5.6 8 9l7.2-3.4z" fill="#1f5fbf" stroke="#123c7a" stroke-width=".7"/><path d="M3.6 7.3v3.3c0 1.3 2 2.4 4.4 2.4s4.4-1.1 4.4-2.4V7.3L8 9.4z" fill="#8fb3ea" stroke="#123c7a" stroke-width=".7"/><path d="M14.4 6v4.6" stroke="#123c7a" stroke-width="1.1"/><circle cx="14.4" cy="11.3" r="1" fill="#123c7a"/>',
    debug: '<ellipse cx="8" cy="9" rx="3.5" ry="4.5" fill="#8d5524"/><path d="M8 4.5v9M3 6l2 1.3M13 6l-2 1.3M2.5 10h2M11.5 10h2M3 13.5l2-1.3M13 13.5l-2-1.3" stroke="#3e2410" stroke-width="1.2"/><circle cx="8" cy="4" r="1.8" fill="#3e2410"/>',
  };
  function icon(name) { return `<svg viewBox="0 0 16 16" xmlns="http://www.w3.org/2000/svg" aria-hidden="true">${ICON[name] || ""}</svg>`; }
  const D = () => S.D;
  const MAIN_TB = [
    ["new", "Create a New Model (ACTM)", () => S.command(`ACTM,${S.lowestUnusedModel()}`)],
    ["open", "Open the Model Database", () => D().loadModel()],
    ["save", "Save a Model (SAVE)", () => S.command("SAVE")],
    ["outpre", "Output model to .pre (WRITE)", () => D().output()],
    ["outansys", "Output model to ANSYS (ANSYS)", () => D().exportAnsys()],
    ["H", "use .hou converter (CONVERT,SSI)", () => D().converter("SSI")],
    ["A", "use ANSYS converter (CONVERT,ANSYS)", () => D().converter("ANSYS")],
    ["Sx", "use GT STRUDL converter (not included)", null],
    "|",
    ["camera", "capture a plot image (CAPTUREPLOT)", () => D().exportImage(), () => !!(S.state && S.state.plots.active)],
    "|",
    ["element", "Element Plot (MODELPLOT)", () => S.command("MODELPLOT")],
    ["node", "Node plot (NODEPLOT)", () => S.command("NODEPLOT")],
    ["cut", "Cut Plot (CUTPLOT)", () => D().cutPlot()],
    ["spectrum", "Spectrum Plot (READSPEC + SPECPLOT)", () => D().lineSelection("SPECPLOT")],
    ["th", "Time History Plot (READTH + THPLOT)", () => D().lineSelection("THPLOT")],
    ["layers", "Soil Layer Plot (LAYERPLOT)", () => S.command("LAYERPLOT")],
    ["soilprop", "Soil Properties Plot (SOILPROPPLOT)", () => D().soilProperty()],
    "|",
    ["frames", "Process Animation Frames (PROCFRAME)", () => D().procFrame()],
    ["bubble", "Bubble Plot (BUBBLEPLOT)", () => D().loadFrameData("BUBBLEPLOT")],
    ["vector", "Vector Plot (VECTORPLOT)", () => D().loadFrameData("VECTORPLOT")],
    ["contour", "Contour Plot (CONTOURPLOT)", () => D().loadFrameData("CONTOURPLOT")],
    ["deformed", "Deformed Shape Plot (DEFORMPLOT)", () => D().loadFrameData("DEFORMPLOT")],
    "|",
    ["options", "Change Analysis Options", () => D().optionsDialog("ANALYSIS")],
    ["afwrite", "Afwrite current model (AFWRITE)", () => S.command("AFWRITE")],
    ["check", "Check Errors window (CHECK)", () => D().checkErrors()],
    ["results", "Results browser", () => S.openResults()],
    "|",
    ["learn", "Learn: guided course, examples and command explainer", () => SASSI.Learn.openStart()],
  ];
  const PLOT_TB = [
    ["view", "Change View of Current Plot (CNGVIEW)", () => D().changeView(), "CNGVIEW"],
    ["resetview", "Reset current plot view to default (RSTVIEW)", () => S.command("RSTVIEW"), "RSTVIEW"],
    ["center", "Change center of current plot (CNGCENTER)", () => D().changeCenter(), "CNGCENTER"],
    ["resetcenter", "Reset Center of current plot to default (RSTCENTER)", () => S.command("RSTCENTER"), "RSTCENTER"],
    "|",
    ["wire", "Wireframe view (WIREFRAME)", () => S.command("WIREFRAME"), "WIREFRAME", (v) => v.wireframe],
    ["shrink", "Shrink faces of elements (SHRINK)", () => S.command("SHRINK"), "SHRINK", (v) => v.shrink],
    "|",
    ["G", "Colors of elements determined by group (ELECOLOR,1)", () => S.command("ELECOLOR,1"), "ELECOLOR", (v) => v.color_by === 1],
    ["M", "Colors of elements determined by material (ELECOLOR,2)", () => S.command("ELECOLOR,2"), "ELECOLOR", (v) => v.color_by === 2],
    ["Pp", "Colors of elements determined by property (ELECOLOR,3)", () => S.command("ELECOLOR,3"), "ELECOLOR", (v) => v.color_by === 3],
    "|",
    ["nodenum", "Show node labels (NODENUM)", () => S.command("NODENUM"), "NODENUM", (v) => v.node_labels],
    ["elemnum", "Show element labels (ELENUM)", () => S.command("ELENUM"), "ELENUM", (v) => v.elem_labels],
    ["groupnum", "Show group labels (GROUPNUM)", () => S.command("GROUPNUM"), "GROUPNUM", (v) => v.group_labels],
    ["showdof", "show fixed degrees of freedom (SHOWDOF)", () => D().showDof(), "SHOWDOF", (v) => v.show_dof && v.show_dof.length > 0],
    ["showmass", "show lumped masses (SHOWMASS)", () => S.command("SHOWMASS"), "SHOWMASS", (v) => v.show_mass],
    "|",
    ["pause", "Pause Start Animation (PAUSE)", () => P.togglePause(), "PAUSE", (v) => !v.paused],
    ["debug", "Show Debug info (DEBUG)", () => S.command("DEBUG"), "DEBUG", (v) => v.debug],
  ];
  const Toolbars = {buttons: []};
  SASSI.Toolbars = Toolbars;
  Toolbars.build = function () {
    const mk = (bar, defs, plot) => {
      for (const d of defs) {
        if (d === "|") { bar.appendChild(el("span", {class: "tb-sep"})); continue; }
        const b = el("button", {class: "tb", title: d[1], "aria-label": d[1], html: icon(d[0])});
        if (!d[2]) b.disabled = true;
        else b.addEventListener("click", () => d[2]());
        Toolbars.buttons.push({b, def: d, plot});
        bar.appendChild(b);
      }
    };
    const main = S.$("#toolbar-main"), pl = S.$("#toolbar-plot");
    main.innerHTML = ""; pl.innerHTML = "";
    mk(main, MAIN_TB, false);
    pl.appendChild(el("span", {class: "tb-label", text: "3D Plot"}));
    mk(pl, PLOT_TB, true);
  };
  Toolbars.update = function () {
    const rec = P.activePlot();
    for (const {b, def, plot} of Toolbars.buttons) {
      if (!def[2]) continue;
      if (plot) {
        const cmd = def[3];
        const ok = !!rec && (CAP[cmd] || THREE_D).includes(rec.kind);
        b.disabled = !ok;
        b.classList.toggle("on", !!(ok && def[4] && rec.view && def[4](rec.view)));
      } else if (def[3]) {
        b.disabled = !def[3]();
      }
    }
  };
})(SASSI);
