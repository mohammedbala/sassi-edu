/* SASSI-EDU GUI -- concept figures of the guided course (the ```figure blocks of the lessons).
 *
 * A lesson writes  ```figure  /  <name> [key=value ...]  /  caption (Markdown)  /  ```  (format:
 * docs/internal/lesson_format.md).  sassi/ui/learn.py renders the block as a placeholder
 * <div class="lesson-figure" data-figure-block="k"> and the caption as HTML; learn.js calls S.Figures.fill,
 * which mounts the registered figure there.  Every figure is drawn on a canvas in vanilla JavaScript and
 * computes what it shows from the formulas named in its formula strip (the physics functions below are
 * pure and are also run by tests/unit/test_lesson_figures.py in Node against the Python modules):
 *
 *   fixed-base-vs-ssi    one-mode stick on a rigid base and on sway-rocking springs and dashpots (closed-form
 *                        3-DOF harmonic solution)
 *   kinematic-inertial   SSI = kinematic interaction (massless) + inertial interaction (schematic)
 *   substructuring       free field + structure - excavated soil = SSI system (manual Eq. 2.1, build-up)
 *   soil-column          vertically propagating SH wave in a layered column on a half-space (exact SHAKE
 *                        recursion with the complex modulus of Theory Section 2)
 *   impedance-ellipse    K(w) = k + i w c under harmonic displacement: force-displacement ellipse
 *   tf-to-isrs           control motion -> H(f) -> floor motion -> oscillators (Nigam-Jennings) -> ISRS
 *   interaction-sets     FV, FI-FSIN, FI-EVBN and FFV on the excavation mesh of lesson 5
 *   wave-types           vertically propagating SV, SH and P waves (X, Y, Z input)
 *   isrs-broadening      envelope of soil cases and +-b peak broadening (BROADEN window maximum)
 *   hysteresis           backbone, Masing loop, secant stiffness, xi = E_D / (4 pi E_S)
 *   freq-interpolation   computed vs interpolated transfer function (Tajirian windows, option 1), CRITFREQ
 *
 * Animated figures have Play / Pause; they stop while off screen (IntersectionObserver), start paused
 * under prefers-reduced-motion, and stop for good when their lesson page is replaced.  Colours come from
 * the CSS variables of styles.css (.lfig).  No innerHTML: captions arrive as nodes from learn.js. */
"use strict";

(function (S) {
  const FIG = {};            // registry: name -> figure definition
  const PH = {};             // physics (pure functions, no DOM)
  const TAU = 2 * Math.PI;

  // ================================================================== complex arithmetic ([re, im])
  const C = {
    add: (a, b) => [a[0] + b[0], a[1] + b[1]],
    sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
    mul: (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]],
    div: (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; },
    scale: (a, s) => [a[0] * s, a[1] * s],
    abs: (a) => Math.hypot(a[0], a[1]),
    arg: (a) => Math.atan2(a[1], a[0]),
    exp: (a) => { const e = Math.exp(a[0]); return [e * Math.cos(a[1]), e * Math.sin(a[1])]; },
    conj: (a) => [a[0], -a[1]],
    sqrt: (a) => {
      const r = Math.hypot(a[0], a[1]);
      const re = Math.sqrt(Math.max(0, (r + a[0]) / 2));
      let im = Math.sqrt(Math.max(0, (r - a[0]) / 2));
      if (a[1] < 0) im = -im;
      return [re, im];
    },
  };
  PH.C = C;
  /** SASSI complex-modulus factor c(beta) = 1 - 2 beta^2 + 2 i beta sqrt(1 - beta^2) (Theory Eq. 2.1). */
  PH.cfac = (beta) => [1 - 2 * beta * beta, 2 * beta * Math.sqrt(1 - beta * beta)];

  /** Solve the complex linear system A x = b (Gaussian elimination, partial pivoting).  Returns
   *  {x, pivot} with pivot = smallest |pivot| / largest |entry| (0 when singular). */
  PH.csolve = function (A0, b0) {
    const n = b0.length;
    const A = A0.map((r) => r.map((v) => [v[0], v[1]]));
    const b = b0.map((v) => [v[0], v[1]]);
    let big = 0;
    for (const r of A) for (const v of r) big = Math.max(big, C.abs(v));
    let minp = Infinity;
    for (let c = 0; c < n; c++) {
      let p = c, pv = C.abs(A[c][c]);
      for (let r = c + 1; r < n; r++) { const v = C.abs(A[r][c]); if (v > pv) { pv = v; p = r; } }
      minp = Math.min(minp, pv);
      if (pv === 0) return {x: null, pivot: 0};
      if (p !== c) { [A[p], A[c]] = [A[c], A[p]]; [b[p], b[c]] = [b[c], b[p]]; }
      for (let r = c + 1; r < n; r++) {
        const f = C.div(A[r][c], A[c][c]);
        if (f[0] === 0 && f[1] === 0) continue;
        for (let k = c; k < n; k++) A[r][k] = C.sub(A[r][k], C.mul(f, A[c][k]));
        b[r] = C.sub(b[r], C.mul(f, b[c]));
      }
    }
    const x = new Array(n);
    for (let r = n - 1; r >= 0; r--) {
      let s = b[r];
      for (let k = r + 1; k < n; k++) s = C.sub(s, C.mul(A[r][k], x[k]));
      x[r] = C.div(s, A[r][r]);
    }
    return {x, pivot: big > 0 ? minp / big : 0};
  };

  // ================================================================== FFT (radix 2, in place)
  /** In-place complex FFT of (re, im), length a power of 2; sign -1 forward (numpy convention), +1 inverse
   *  (without the 1/N). */
  PH.fft = function (re, im, sign) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = sign * TAU / len, wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < len / 2; k++) {
          const a = i + k, b = a + len / 2;
          const tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
          re[b] = re[a] - tr; im[b] = im[a] - ti;
          re[a] += tr; im[a] += ti;
          const ncr = cr * wr - ci * wi;
          ci = cr * wi + ci * wr; cr = ncr;
        }
      }
    }
  };
  /** rfft of a real signal zero-padded to nfft: {re, im} of the bins 0..nfft/2. */
  PH.rfft = function (a, nfft) {
    const re = new Float64Array(nfft), im = new Float64Array(nfft);
    for (let i = 0; i < Math.min(a.length, nfft); i++) re[i] = a[i];
    PH.fft(re, im, -1);
    return {re: re.slice(0, nfft / 2 + 1), im: im.slice(0, nfft / 2 + 1)};
  };
  /** Inverse of rfft (with 1/N): the real signal of length nfft. */
  PH.irfft = function (R, I, nfft) {
    const re = new Float64Array(nfft), im = new Float64Array(nfft);
    for (let k = 0; k <= nfft / 2; k++) { re[k] = R[k]; im[k] = I[k]; }
    for (let k = 1; k < nfft / 2; k++) { re[nfft - k] = R[k]; im[nfft - k] = -I[k]; }
    im[0] = 0; im[nfft / 2] = 0;
    PH.fft(re, im, +1);
    for (let i = 0; i < nfft; i++) re[i] /= nfft;
    return re;
  };

  // ================================================================== 1. SDOF on sway-rocking springs
  /** One-mode idealisation of the stick of example 1 (lessons 1 and 4): mode 1 of the fixed-base stick
   *  (M1 = 3063 t, h1 = 15.6 m, 5.05 Hz, beta 5 %; Gamma1 phi1 = 1.34 at the roof), the 20 m x 20 m mat
   *  (1468 t at grade), and SASSI's impedance of that mat at 3.49 Hz (lesson 4, FOUNSTIF / FOUNDAMP):
   *  Kx = 1.48e7 kN/m, Ktheta = 1.70e9 kN m/rad, damping ratios 7.7 % (sliding) and 4.9 % (rocking), of
   *  which 4.8 % is the soil's material damping (hysteretic) and the rest radiation (a dashpot fitted at
   *  3.49 Hz).  Units: kN, m, t, s. */
  PH.SSI = {M: 3062.6, h: 15.59, ffb: 5.054, beta: 0.05, m0: 1468, B: 20, Kx: 1.48e7, Kt: 1.70e9,
    xim: 0.048, xix: 0.077, xit: 0.049, fref: 3.49};
  /** The model for a velocity factor s of the site (springs x s^2, dashpots x s). */
  PH.ssiModel = function (s, p) {
    p = Object.assign({}, PH.SSI, p || {});
    s = s || 1;
    const wfb = TAU * p.ffb, wref = TAU * p.fref;
    return {p, s, k: wfb * wfb * p.M, I0: p.m0 * p.B * p.B / 12,
      Kx: p.Kx * s * s, Kt: p.Kt * s * s,
      cx: 2 * (p.xix - p.xim) * p.Kx / wref * s, ct: 2 * (p.xit - p.xim) * p.Kt / wref * s};
  };
  /** Steady state at f (Hz) per unit free-field displacement u_g: u (stick deformation), u0 (mat sway
   *  relative to the free field), th (rocking, rad/m of u_g), Ht = total motion of the mass, Hfb = total
   *  motion of the mass on a fixed base, prad = radiated power per unit |u_g|^2. */
  PH.ssiResponse = function (m, f) {
    const p = m.p, w = TAU * f, w2 = w * w;
    const ks = C.scale(PH.cfac(p.beta), m.k);
    const Sx = [m.Kx, 2 * p.xim * m.Kx + w * m.cx];
    const St = [m.Kt, 2 * p.xim * m.Kt + w * m.ct];
    const M = p.M, h = p.h;
    const A = [
      [C.sub(ks, [w2 * M, 0]), [-w2 * M, 0], [-w2 * M * h, 0]],
      [[-w2 * M, 0], C.sub(Sx, [w2 * (M + p.m0), 0]), [-w2 * M * h, 0]],
      [[-w2 * M * h, 0], [-w2 * M * h, 0], C.sub(St, [w2 * (M * h * h + m.I0), 0])],
    ];
    const b = [[w2 * M, 0], [w2 * (M + p.m0), 0], [w2 * M * h, 0]];
    const [u, u0, th] = PH.csolve(A, b).x;
    const Ht = C.add([1, 0], C.add(u0, C.add(C.scale(th, h), u)));
    const Hfb = C.div(ks, C.sub(ks, [w2 * M, 0]));
    const prad = 0.5 * w2 * (m.cx * (u0[0] ** 2 + u0[1] ** 2) + m.ct * (th[0] ** 2 + th[1] ** 2));
    return {u, u0, th, Ht, Hfb, prad};
  };
  /** Peak and half-power damping of |H(f)| sampled on fs. */
  PH.peakInfo = function (fs, amp) {
    let i = 0;
    for (let k = 1; k < amp.length; k++) if (amp[k] > amp[i]) i = k;
    const half = amp[i] / Math.SQRT2;
    const cross = (k0, dk) => {
      for (let k = k0; k + dk >= 0 && k + dk < amp.length; k += dk) {
        if (amp[k + dk] < half) return fs[k] + (fs[k + dk] - fs[k]) * (amp[k] - half) / (amp[k] - amp[k + dk]);
      }
      return NaN;
    };
    const f1 = cross(i, -1), f2 = cross(i, 1);
    return {f: fs[i], amp: amp[i], zeta: (f2 - f1) / (2 * fs[i]), i};
  };

  // ================================================================== 2. layered soil column (SHAKE)
  /** Lesson 2 site (example 1): 5 m sand over 12 m gravel over weathered rock; w = unit weight kN/m3. */
  PH.COLUMN = {g: 9.81, layers: [{h: 5, vs: 300, w: 19.0, beta: 0.05, name: "sand"}, {h: 12, vs: 500, w: 20.0, beta: 0.04, name: "gravel"}],
    hs: {vs: 1000, w: 21.0, beta: 0.02, name: "weathered rock"}};
  /** Complex velocity V* = V (sqrt(1 - beta^2) + i beta) (Theory Eq. 2.2). */
  PH.vstar = (v, beta) => [v * Math.sqrt(1 - beta * beta), v * beta];
  /** Up- and down-going amplitudes E_m, F_m of every layer (the half-space last) at f, with E_1 = F_1 = 1
   *  (Theory Eq. 12.2).  Returns {E, F, vs (complex V*), rho, h}. */
  PH.columnWaves = function (col, f) {
    const L = col.layers.concat([col.hs]);
    const n = L.length, w = TAU * f;
    const rho = L.map((l) => l.w / col.g), vs = L.map((l) => PH.vstar(l.vs, l.beta));
    const E = [[1, 0]], F = [[1, 0]];
    for (let m = 0; m < n - 1; m++) {
      const alpha = C.div(C.scale(vs[m], rho[m]), C.scale(vs[m + 1], rho[m + 1]));
      const ikh = C.mul([0, w * L[m].h], C.div([1, 0], vs[m]));
      const ep = C.exp(ikh), em = C.exp(C.scale(ikh, -1));
      const onep = C.add([1, 0], alpha), onem = C.sub([1, 0], alpha);
      E.push(C.scale(C.add(C.mul(C.mul(E[m], onep), ep), C.mul(C.mul(F[m], onem), em)), 0.5));
      F.push(C.scale(C.add(C.mul(C.mul(E[m], onem), ep), C.mul(C.mul(F[m], onep), em)), 0.5));
    }
    return {E, F, vs, rho, h: L.map((l) => l.h), w};
  };
  /** Displacement at depth z (m, from the surface) per unit surface motion. */
  PH.columnU = function (wv, col, z) {
    let top = 0, m = 0;
    while (m < col.layers.length && z > top + col.layers[m].h) { top += col.layers[m].h; m++; }
    const k = C.div([wv.w, 0], wv.vs[m]);
    const zm = z - top;
    const u = C.add(C.mul(wv.E[m], C.exp([-k[1] * zm, k[0] * zm])),
      C.mul(wv.F[m], C.exp([k[1] * zm, -k[0] * zm])));
    return C.scale(u, 0.5);              // the surface motion E_1 + F_1 = 2
  };
  /** Surface / outcrop (2E) and surface / within (E + F) at the top of the half-space. */
  PH.columnRatios = function (col, f) {
    const wv = PH.columnWaves(col, f);
    const n = wv.E.length - 1;
    return {outcrop: C.div([2, 0], C.scale(wv.E[n], 2)), within: C.div([2, 0], C.add(wv.E[n], wv.F[n])), wv};
  };

  // ================================================================== 3. synthetic motion, oscillators
  /** Deterministic synthetic acceleration (g): random phases (seeded), Kanai-Tajimi / Clough-Penzien
   *  amplitude spectrum (fg = 3.5 Hz, zeta_g = 0.6, high-pass 0.5 Hz), envelope rise 1.5 s, strong 4 s, decay;
   *  scaled to the peak `pga`.  Illustrative only. */
  PH.synthMotion = function (o) {
    o = Object.assign({dt: 0.01, dur: 10, nfft: 2048, pga: 0.3, seed: 11975, fg: 3.5, zg: 0.6, ff: 0.5, fmax: 25}, o || {});
    let st = o.seed >>> 0;
    const rnd = () => { st = (Math.imul(1664525, st) + 1013904223) >>> 0; return st / 4294967296; };
    const n = o.nfft, df = 1 / (n * o.dt);
    const R = new Float64Array(n / 2 + 1), I = new Float64Array(n / 2 + 1);
    for (let k = 1; k < n / 2; k++) {
      const f = k * df;
      const phase = TAU * rnd();
      if (f > o.fmax) continue;
      const r = f / o.fg, q = f / o.ff;
      const kt = Math.sqrt((1 + 4 * o.zg * o.zg * r * r) / ((1 - r * r) ** 2 + 4 * o.zg * o.zg * r * r));
      const hp = (q * q) / Math.sqrt((1 - q * q) ** 2 + 4 * o.zg * o.zg * q * q);
      R[k] = kt * hp * Math.cos(phase); I[k] = kt * hp * Math.sin(phase);
    }
    const a = PH.irfft(R, I, n);
    const nd = Math.round(o.dur / o.dt);
    let peak = 0;
    for (let i = 0; i < n; i++) {
      const t = i * o.dt;
      const env = i >= nd ? 0 : t < 1.5 ? (t / 1.5) ** 2 : t < 5.5 ? 1 : Math.exp(-0.6 * (t - 5.5));
      a[i] *= env;
      peak = Math.max(peak, Math.abs(a[i]));
    }
    for (let i = 0; i < n; i++) a[i] *= o.pga / peak;
    return {a, dt: o.dt, nfft: n, n: nd};
  };
  /** Nigam-Jennings coefficients (Theory Section 11.3) for x'' + 2 z w x' + w^2 x = -a(t), a linear in a step. */
  PH.njCoefficients = function (w, z, dt) {
    const sq = Math.sqrt(1 - z * z), wd = w * sq;
    const E = Math.exp(-z * w * dt), Sn = Math.sin(wd * dt), Cs = Math.cos(wd * dt);
    const a11 = E * (z / sq * Sn + Cs), a12 = E * Sn / wd, a21 = -w / sq * E * Sn, a22 = E * (Cs - z / sq * Sn);
    const w2 = w * w, w3 = w2 * w, t1 = (2 * z * z - 1) / (w2 * dt), t2 = 2 * z / (w3 * dt);
    const b11 = E * ((t1 + z / w) * Sn / wd + (t2 + 1 / w2) * Cs) - t2;
    const b12 = -E * (t1 * Sn / wd + t2 * Cs) - 1 / w2 + t2;
    const b21 = E * ((t1 + z / w) * (Cs - z / sq * Sn) - (t2 + 1 / w2) * (wd * Sn + z * w * Cs)) + 1 / (w2 * dt);
    const b22 = -E * (t1 * (Cs - z / sq * Sn) - t2 * (wd * Sn + z * w * Cs)) - 1 / (w2 * dt);
    return [a11, a12, a21, a22, b11, b12, b21, b22];
  };
  /** Oscillator (f Hz, damping z) under base acceleration a (sampled at dt), the input linear between the
   *  samples and sub-divided so that a period has at least 64 steps (as MOTION does).  Returns the relative
   *  displacement x at the samples, the running maximum of the absolute acceleration |2 z w x' + w^2 x| up to
   *  each sample (between-sample peaks included) and the final maximum (the spectral acceleration). */
  PH.oscillator = function (a, dt, f, z, nmax) {
    const n = Math.min(nmax || a.length, a.length), w = TAU * f;
    const m = Math.max(1, Math.min(64, Math.ceil(64 * f * dt)));
    const [a11, a12, a21, a22, b11, b12, b21, b22] = PH.njCoefficients(w, z, dt / m);
    const xs = new Float64Array(n), run = new Float64Array(n);
    let x = 0, v = 0, mx = 0;
    for (let i = 0; i < n - 1; i++) {
      const a0 = a[i], da = (a[i + 1] - a0) / m;
      for (let j = 0; j < m; j++) {
        const p = a0 + j * da, q = p + da;
        const xn = a11 * x + a12 * v + b11 * p + b12 * q;
        v = a21 * x + a22 * v + b21 * p + b22 * q;
        x = xn;
        const aa = Math.abs(2 * z * w * v + w * w * x);
        if (aa > mx) mx = aa;
      }
      xs[i + 1] = x; run[i + 1] = mx;
    }
    return {x: xs, run, sa: mx};
  };
  /** Response spectrum (absolute acceleration) at the frequencies fs. */
  PH.spectrum = function (a, dt, fs, z, nmax) { return fs.map((f) => PH.oscillator(a, dt, f, z, nmax).sa); };
  /** Log-spaced frequencies. */
  PH.logspace = (f1, f2, n) => Array.from({length: n}, (_, i) => f1 * Math.pow(f2 / f1, i / (n - 1)));
  /** Total-motion transfer function of one hysteretic mode (f0, beta) with participation P at the output:
   *  H = 1 + P r^2 / (c(beta) - r^2), r = f/f0 (P = 1: Theory Eq. 2.4). */
  PH.modalTF = function (f, modes) {
    let H = [1, 0];
    for (const md of modes) {
      const r2 = (f / md.f) ** 2;
      H = C.add(H, C.div([md.P * r2, 0], C.sub(PH.cfac(md.beta), [r2, 0])));
    }
    return H;
  };
  /** Floor acceleration a(t) = IFFT[H(f) A(f)] (Theory Eq. 11.1) for a TF function of f. */
  PH.convolve = function (mo, tf) {
    const n = mo.nfft, df = 1 / (n * mo.dt);
    const A = PH.rfft(mo.a, n);
    const R = new Float64Array(n / 2 + 1), I = new Float64Array(n / 2 + 1);
    for (let k = 0; k <= n / 2; k++) {
      const H = tf(k * df);
      R[k] = A.re[k] * H[0] - A.im[k] * H[1];
      I[k] = A.re[k] * H[1] + A.im[k] * H[0];
    }
    return PH.irfft(R, I, n);
  };

  // ================================================================== 4. peak broadening (BROADEN)
  function interpLin(X, Y, x) {
    if (x <= X[0]) return Y[0];
    if (x >= X[X.length - 1]) return Y[Y.length - 1];
    let lo = 0, hi = X.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (X[mid] <= x) lo = mid; else hi = mid; }
    return Y[lo] + (Y[hi] - Y[lo]) * (x - X[lo]) / (X[hi] - X[lo]);
  }
  PH.interpLin = interpLin;
  /** Envelope of several spectra on the union of their frequencies. */
  PH.envelope = function (lines) {
    const X = Array.from(new Set(lines.flatMap((l) => l.x))).sort((a, b) => a - b);
    return {x: X, y: X.map((x) => Math.max(...lines.map((l) => interpLin(l.x, l.y, x))))};
  };
  /** Peak broadening by +-b: B(f) = max{E(f') : f/(1+b) <= f' <= f/(1-b)} (BROADEN step 2, the window
   *  maximum of sassi/plotting/lines.py), evaluated on X, X(1-b), X(1+b). */
  PH.broaden = function (X, E, b) {
    if (!(b > 0)) return {x: X.slice(), y: E.slice()};
    const G = Array.from(new Set(X.concat(X.map((x) => x * (1 - b)), X.map((x) => x * (1 + b)))))
      .filter((x) => x >= X[0] && x <= X[X.length - 1]).sort((p, q) => p - q);
    const y = G.map((f) => {
      const lo = f / (1 + b), hi = f / (1 - b);
      let m = Math.max(interpLin(X, E, lo), interpLin(X, E, hi));
      for (let i = 0; i < X.length; i++) if (X[i] > lo && X[i] < hi && E[i] > m) m = E[i];
      return m;
    });
    return {x: G, y};
  };

  // ================================================================== 5. backbones and Masing loops
  /** Hyperbolic backbone tau = G0 g / (1 + |g|/gr) (Theory Eq. 18.1, beta = s = 1), its integral. */
  PH.hyperbolic = (G0, gr) => ({
    kel: G0, xcr: 0,
    F: (x) => G0 * x / (1 + Math.abs(x) / gr),
    I: (x) => { const z = Math.abs(x) / gr; return G0 * gr * gr * (z - Math.log1p(z)); },
  });
  /** Piecewise-linear odd backbone through the origin and (xs, ys), constant beyond the last point. */
  PH.polyline = function (xs, ys) {
    const X = [0].concat(xs), Y = [0].concat(ys);
    const area = [0];
    for (let i = 1; i < X.length; i++) area.push(area[i - 1] + (X[i] - X[i - 1]) * (Y[i] + Y[i - 1]) / 2);
    const Fp = (x) => (x >= X[X.length - 1] ? Y[Y.length - 1] : interpLin(X, Y, x));
    return {
      kel: ys[0] / xs[0], xcr: xs[0],
      F: (x) => Math.sign(x) * Fp(Math.abs(x)),
      I: (x) => {
        x = Math.abs(x);
        if (x >= X[X.length - 1]) return area[area.length - 1] + (x - X[X.length - 1]) * Y[Y.length - 1];
        let i = 1;
        while (X[i] < x) i++;
        const y = Fp(x);
        return area[i - 1] + (x - X[i - 1]) * (Y[i - 1] + y) / 2;
      },
    };
  };
  /** Stabilised symmetric Masing loop at amplitude a (Theory Section 17, GMR): E_D = 8 I(a) - 4 a F(a),
   *  E_S = a F(a) / 2, xi_h = E_D / (4 pi E_S), K_sec = F(a)/a; the branches of the loop. */
  PH.masing = function (bb, a) {
    const Fa = bb.F(a);
    const ED = Math.max(0, 8 * bb.I(a) - 4 * a * Fa), ES = 0.5 * a * Fa;
    return {Fa, ED, ES, xi: ES > 0 ? ED / (4 * Math.PI * ES) : 0, ksec: Fa / a, ratio: Fa / a / bb.kel,
      down: (x) => Fa + 2 * bb.F((x - a) / 2),           // unloading from +a
      up: (x) => -Fa + 2 * bb.F((x + a) / 2)};            // reloading from -a
  };
  /** SHAKE91 library 'Sand' (sassi/data/dynp_library.pre: Seed & Idriss 1970 upper range, Idriss 1990):
   *  strain %, G/Gmax, damping %. */
  PH.SAND = {g: [1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1, 3, 10],
    G: [1, 1, 0.99, 0.96, 0.85, 0.64, 0.37, 0.18, 0.08, 0.05, 0.035],
    D: [0.24, 0.42, 0.8, 1.4, 2.8, 5.1, 9.8, 15.5, 21, 25, 28]};
  /** Reference strain (%) where the library curve has G/Gmax = 0.5 (linear in log10 strain, SHAKE91 rule). */
  PH.sandRefStrain = function () {
    const s = PH.SAND;
    for (let i = 1; i < s.g.length; i++) {
      if (s.G[i] <= 0.5) {
        const t = (s.G[i - 1] - 0.5) / (s.G[i - 1] - s.G[i]);
        return Math.pow(10, Math.log10(s.g[i - 1]) + t * (Math.log10(s.g[i]) - Math.log10(s.g[i - 1])));
      }
    }
    return NaN;
  };
  /** BBCGEN backbone (Theory Section 17): cracking (Vcr/GA, Vcr), 20 points up to the yield (gy, Vu) on
   *  V = Vcr + (Vu - Vcr)[1 - (1 - xi)^2], failure (0.02, 1.02 Vu).  Lesson 9 wall: Vu = 12 472 kN,
   *  Vcr = 0.3 Vu, G A_W = 4.5e7 kN (Vcr / gamma_cr with gamma_cr = 8.31e-5). */
  PH.bbcgen = function (vu, vcr, ga, gy, gf) {
    gy = gy || 0.004; gf = gf || 0.02;
    const gcr = vcr / ga, xs = [gcr], ys = [vcr];
    for (let k = 1; k <= 20; k++) { const xi = k / 20; xs.push(gcr + xi * (gy - gcr)); ys.push(vcr + (vu - vcr) * (1 - (1 - xi) ** 2)); }
    xs.push(gf); ys.push(1.02 * vu);
    return {xs, ys};
  };
  PH.PANEL = {vu: 12472, vcr: 3742, gcr: 8.31e-5};

  // ================================================================== 6. TF interpolation (Tajirian)
  function quadRoots(a, b, c) {
    // roots of a x^2 + b x + c (complex), cancellation-free; missing roots -> null
    if (C.abs(a) === 0) { if (C.abs(b) === 0) return [null, null]; return [C.div(C.scale(c, -1), b), null]; }
    let disc = C.sqrt(C.sub(C.mul(b, b), C.scale(C.mul(a, c), 4)));
    if ((b[0] * disc[0] + b[1] * disc[1]) < 0) disc = C.scale(disc, -1);
    const q = C.scale(C.add(b, disc), -0.5);
    return [C.div(q, a), C.abs(q) === 0 ? [0, 0] : C.div(c, q)];
  }
  const evalRat = (Cs, x) => C.div(C.add(C.scale(C.add(C.scale(Cs[0], x), Cs[1]), x), Cs[2]),
    C.add(C.scale(C.add([x, 0], Cs[3]), x), Cs[4]));
  /** Fit the 2-DOF hysteretic form H = (C1 x^2 + C2 x + C3)/(x^2 + C4 x + C5), x = (w/w_max)^2, to five
   *  points (Theory Eq. 10.1-10.2) and decide the spurious-pole guard (Theory Section 10.1). */
  PH.fitWindow = function (f5, H5) {
    const scale = TAU * f5[4];
    const x = f5.map((f) => (TAU * f / scale) ** 2);
    const A = x.map((xp, p) => [[xp * xp, 0], [xp, 0], [1, 0], C.scale(H5[p], -xp), C.scale(H5[p], -1)]);
    const b = x.map((xp, p) => C.scale(H5[p], xp * xp));
    let sol = PH.csolve(A, b);
    let Cs = sol.x;
    const hmax = Math.max(...H5.map(C.abs)) || 1;
    const nodeErr = (Cv) => Math.max(...x.map((xp, p) => {
      const num = C.add(C.scale(C.add(C.scale(Cv[0], xp), Cv[1]), xp), Cv[2]);
      const den = C.add(C.scale(C.add([xp, 0], Cv[3]), xp), Cv[4]);
      return C.abs(C.sub(num, C.mul(H5[p], den))) / C.abs(den);
    })) / hmax;
    let err = Cs ? nodeErr(Cs) : Infinity;
    if (!Cs || !isFinite(err) || sol.pivot < 1e-13) {
      // (near) rank deficient: regularised normal equations, close to the minimum-norm solution
      const AH = (i, j) => x.reduce((s, _, p) => C.add(s, C.mul(C.conj(A[p][i]), A[p][j])), [0, 0]);
      const N = [0, 1, 2, 3, 4].map((i) => [0, 1, 2, 3, 4].map((j) => AH(i, j)));
      const tr = N.reduce((s, r, i) => s + r[i][0], 0);
      for (let i = 0; i < 5; i++) N[i][i] = C.add(N[i][i], [1e-14 * tr, 0]);
      const rhs = [0, 1, 2, 3, 4].map((i) => x.reduce((s, _, p) => C.add(s, C.mul(C.conj(A[p][i]), b[p])), [0, 0]));
      const s2 = PH.csolve(N, rhs);
      if (s2.x) { const e2 = nodeErr(s2.x); if (!(e2 >= err)) { Cs = s2.x; err = e2; } }
    }
    let guard = !(err <= 1e-6);
    if (!guard) {
      // trigger: min|D| < 1e-3 max|D| on [x0, 1]; extrema at the ends and the real roots of d|D|^2/dx
      const p4 = Cs[3], p5 = Cs[4], x0 = x[0];
      const D = (t) => C.abs(C.add(C.scale(C.add([t, 0], p4), t), p5));
      const dd = (t) => 2 * t ** 3 + 3 * p4[0] * t * t + (p4[0] ** 2 + 2 * p5[0] + p4[1] ** 2) * t + (p4[0] * p5[0] + p4[1] * p5[1]);
      const cand = [x0, 1];
      const NS = 400;
      let prev = dd(x0);
      for (let i = 1; i <= NS; i++) {
        const t = x0 + (1 - x0) * i / NS, cur = dd(t);
        if (prev === 0 || prev * cur < 0) {
          let lo = t - (1 - x0) / NS, hi = t;
          for (let it = 0; it < 60; it++) { const mid = (lo + hi) / 2; if (dd(lo) * dd(mid) <= 0) hi = mid; else lo = mid; }
          cand.push((lo + hi) / 2);
        }
        prev = cur;
      }
      const Ds = cand.map(D);
      const trigger = Math.min(...Ds) < 1e-3 * Math.max(...Ds);
      if (trigger) {
        const P = quadRoots([1, 0], p4, p5), Z = quadRoots(Cs[0], Cs[1], Cs[2]);
        const dist = (p) => C.abs(C.sub(p, [Math.min(1, Math.max(x0, p[0])), 0]));
        const reach = (p) => Math.max(C.abs(C.sub(p, [x0, 0])), C.abs(C.sub(p, [1, 0])));
        const close = (p, z) => p && z && C.abs(C.sub(p, z)) <= 1e-3 * dist(p);
        const pa = [close(P[0], Z[0]), close(P[1], Z[1])], pb = [close(P[0], Z[1]), close(P[1], Z[0])];
        const canc = (pb[0] + pb[1] > pa[0] + pa[1]) ? pb : pa;
        const allZero = Cs.slice(0, 3).every((c) => c[0] === 0 && c[1] === 0);
        guard = P.some((p, i) => {
          if (!p) return false;
          const physical = Math.sin(C.arg(p) / 2) >= 1e-4;
          const remote = dist(p) >= Math.sqrt(1e-3) * reach(p);
          return !(canc[i] || allZero || physical || remote);
        });
      }
    }
    return {C: Cs, scale, guard, err};
  };
  /** Interpolate computed TFs (fs strictly increasing, H complex) to fout with the window scheme of
   *  MOTION option 1 (SASSI 1982: windows 1, 5, 9, ..., the last five points for the tail; Theory Section 10.2).
   *  Below f_1 linear from H(0) = h0 to H_1, above f_N zero, exact at the computed frequencies. */
  PH.interpTF = function (fs, H, fout, h0) {
    const N = fs.length;
    h0 = h0 || [1, 0];
    const fits = {};
    const fitAt = (s) => fits[s] || (fits[s] = PH.fitWindow(fs.slice(s, s + 5), H.slice(s, s + 5)));
    return fout.map((f) => {
      if (f > fs[N - 1] * (1 + 1e-9)) return [0, 0];
      if (f < fs[0]) { const t = f / fs[0]; return C.add(C.scale(h0, 1 - t), C.scale(H[0], t)); }
      for (let j = 0; j < N; j++) if (Math.abs(f - fs[j]) <= 1e-9 * fs[N - 1]) return H[j];
      if (N < 5) {      // complex linear (fewer than 5 points: the lesson never gets here)
        let j = 0; while (fs[j + 1] < f) j++;
        const t = (f - fs[j]) / (fs[j + 1] - fs[j]);
        return C.add(C.scale(H[j], 1 - t), C.scale(H[j + 1], t));
      }
      let j = 0; while (j < N - 2 && fs[j + 1] <= f) j++;          // interval [f_j, f_j+1], 0-based
      const j1 = j + 1;
      const m1 = Math.min(j1 < 1 ? 1 : 1 + 4 * Math.floor((j1 - 1) / 4), N - 4);
      const s = m1 - 1, fit = fitAt(s);
      if (!fit.guard) return evalRat(fit.C, (TAU * f / fit.scale) ** 2);
      // guard: the complex cubic through the 4 window points nearest to f
      const st = Math.abs(f - fs[s]) <= Math.abs(f - fs[s + 4]) ? s : s + 1;
      let out = [0, 0];
      for (let a = 0; a < 4; a++) {
        let Lw = 1;
        for (let b = 0; b < 4; b++) if (a !== b) Lw *= (f - fs[st + b]) / (fs[st + a] - fs[st + b]);
        out = C.add(out, C.scale(H[st + a], Lw));
      }
      return out;
    });
  };
  /** CRITFREQ (spec 09 section 2.4): peaks of |TFI| above (1 - minfilter/100) max not supported by the larger
   *  computed amplitude of the two bracketing SSI frequencies within tol %. */
  PH.critfreq = function (fu, au, fi, ai, tol, minfilter, df) {
    const amax = Math.max(...ai), thr = (1 - minfilter / 100) * amax, out = [];
    for (let j = 0; j < ai.length; j++) {
      const l = j > 0 ? ai[j - 1] : -Infinity, r = j < ai.length - 1 ? ai[j + 1] : -Infinity;
      if (!(ai[j] > l && ai[j] >= r) || ai[j] < thr) continue;
      let m = -1;
      for (let q = 0; q < fu.length; q++) if (fu[q] <= fi[j]) m = q;
      if (m < 0 || m >= fu.length - 1) continue;
      const ref = Math.max(au[m], au[m + 1]);
      if (!(ref > 0)) continue;
      const d = 100 * Math.abs(ai[j] - ref) / ref;
      out.push({f: fi[j], n: Math.round(fi[j] / df), amp: ai[j], ref, d, flag: d > tol});
    }
    return out;
  };

  // ================================================================== 7. the excavation mesh of lesson 5
  /** Nodes of the 10 m x 10 m x 5 m excavation (5 x 5 per level at 2.5 m, 6 levels at 1 m, z = -5 ... 0) and
   *  the interaction sets of Theory Section 3.3 (INTGEN 1, 3, 2, 5 with skip 2). */
  PH.boxNodes = function (o) {
    o = Object.assign({nx: 5, ny: 5, nz: 6, dx: 2.5, dz: 1}, o || {});
    const out = [];
    for (let k = 0; k < o.nz; k++) for (let j = 0; j < o.ny; j++) for (let i = 0; i < o.nx; i++) {
      out.push({i, j, k, x: (i - (o.nx - 1) / 2) * o.dx, y: (j - (o.ny - 1) / 2) * o.dx, z: -(o.nz - 1) * o.dz + k * o.dz});
    }
    out.o = o;
    return out;
  };
  PH.inSet = function (method, nd, o, skip) {
    const side = nd.i === 0 || nd.j === 0 || nd.i === o.nx - 1 || nd.j === o.ny - 1;
    const fsin = nd.k === 0 || side;
    if (method === "fv") return true;
    if (method === "fsin") return fsin;
    if (method === "evbn") return fsin || nd.k === o.nz - 1;
    if (method === "ffv") return fsin || nd.k === o.nz - 1 || (nd.k > 0 && nd.k < o.nz - 1 && nd.k % (skip || 2) === 0);
    return false;
  };
  PH.setCount = function (method) {
    const nodes = PH.boxNodes();
    return nodes.filter((nd) => PH.inSet(method, nd, nodes.o)).length;
  };

  // ================================================================== 8. waves in a uniform half-space
  /** Ricker wavelet of central frequency fp, peak 1 at t = 0. */
  PH.ricker = (t, fp) => { const a = (Math.PI * fp * t) ** 2; return (1 - 2 * a) * Math.exp(-a); };
  /** Uniform half-space with a free surface: incident up-going pulse plus its reflection,
   *  u(z, t) = p(t + z/V) + p(t - z/V) (z depth, the incident pulse reaches the surface at t = 0). */
  PH.halfspacePulse = (z, t, V, fp) => PH.ricker(t + z / V, fp) + PH.ricker(t - z / V, fp);

  // ================================================================== UI infrastructure (DOM only here)
  const STATE = new Map();          // figure state kept across re-renders of the lesson page
  const reducedMotion = () => !!(typeof matchMedia === "function" && matchMedia("(prefers-reduced-motion: reduce)").matches);
  function el(tag, attrs, ...kids) {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v === null || v === undefined || v === false) continue;
      if (k === "text") e.textContent = v;
      else if (k === "class") e.className = v;
      else if (k.startsWith("on") && typeof v === "function") e.addEventListener(k.slice(2), v);
      else e.setAttribute(k, v === true ? "" : v);
    }
    for (const c of kids.flat()) if (c !== null && c !== undefined && c !== false) e.append(c);
    return e;
  }
  const fmt = (v, d) => (Number.isFinite(v) ? v.toFixed(d === undefined ? 2 : d) : "–");
  const fmtSig = (v, n) => {
    if (!Number.isFinite(v)) return "–";
    if (v === 0) return "0";
    const e = Math.floor(Math.log10(Math.abs(v)));
    if (e >= 5 || e <= -4) return v.toExponential(Math.max(0, (n || 3) - 1)).replace("e+", "e");
    return v.toFixed(Math.max(0, (n || 3) - 1 - e));
  };
  const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
  /** A <span> of text with _{subscript} / ^{superscript} parts as <sub> / <sup>. */
  function richSpan(s) {
    const out = el("span");
    const re = /([_^])\{([^}]*)\}/g;
    let last = 0, m;
    while ((m = re.exec(s))) {
      if (m.index > last) out.append(s.slice(last, m.index));
      out.append(el(m[1] === "_" ? "sub" : "sup", {text: m[2]}));
      last = re.lastIndex;
    }
    out.append(s.slice(last));
    return out;
  }
  const lerp = (a, b, t) => a + (b - a) * t;

  /** Colours of the figure from the CSS variables of .lfig (styles.css). */
  function palette(node) {
    const cs = getComputedStyle(node);
    const v = (n, d) => (cs.getPropertyValue(n) || "").trim() || d;
    return {ink: v("--fig-ink", "#1c222a"), muted: v("--fig-muted", "#5d6774"), grid: v("--fig-grid", "#e3e6ea"),
      axis: v("--fig-axis", "#9aa3ae"), bg: v("--fig-bg", "#ffffff"), s1: v("--fig-s1", "#1f5fbf"), s2: v("--fig-s2", "#c25a00"),
      s3: v("--fig-s3", "#7a4fb5"), s4: v("--fig-s4", "#00897b"), ref: v("--fig-ref", "#8a94a0"), hi: v("--fig-hi", "#c62828"),
      soil1: v("--fig-soil1", "#efe3c6"), soil2: v("--fig-soil2", "#e2cf9f"), rock: v("--fig-rock", "#cfc8bb"),
      soilLine: v("--fig-soil-line", "#b59d6c"), struct: v("--fig-struct", "#6b7684"), structFill: v("--fig-struct-fill", "#d9dee5"),
      ok: v("--fig-ok", "#2e7d32"), font: v("--fig-font", "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif")};
  }
  function withAlpha(hex, a) {
    const m = /^#([0-9a-f]{6})$/i.exec(hex || "");
    if (!m) return hex;
    const n = parseInt(m[1], 16);
    return `rgba(${n >> 16},${(n >> 8) & 255},${n & 255},${a})`;
  }

  /** Drawing helper over a 2D context (CSS pixels). */
  class Pen {
    constructor(g, col) { this.g = g; this.c = col; }
    font(size, weight) { return `${weight || 400} ${size || 11}px ${this.c.font}`; }
    line(pts, color, w, dash) {
      const g = this.g;
      if (pts.length < 2) return;
      g.save(); g.beginPath(); g.strokeStyle = color; g.lineWidth = w || 1.5; g.lineJoin = "round"; g.lineCap = "round";
      if (dash) g.setLineDash(dash);
      g.moveTo(pts[0][0], pts[0][1]);
      for (let i = 1; i < pts.length; i++) g.lineTo(pts[i][0], pts[i][1]);
      g.stroke(); g.restore();
    }
    poly(pts, fill, stroke, w) {
      const g = this.g;
      g.save(); g.beginPath(); g.moveTo(pts[0][0], pts[0][1]);
      for (let i = 1; i < pts.length; i++) g.lineTo(pts[i][0], pts[i][1]);
      g.closePath();
      if (fill) { g.fillStyle = fill; g.fill(); }
      if (stroke) { g.strokeStyle = stroke; g.lineWidth = w || 1; g.lineJoin = "round"; g.stroke(); }
      g.restore();
    }
    rect(x, y, w, h, fill, stroke, lw) { this.poly([[x, y], [x + w, y], [x + w, y + h], [x, y + h]], fill, stroke, lw); }
    circle(x, y, r, fill, stroke, w) {
      const g = this.g;
      g.save(); g.beginPath(); g.arc(x, y, r, 0, TAU);
      if (fill) { g.fillStyle = fill; g.fill(); }
      if (stroke) { g.strokeStyle = stroke; g.lineWidth = w || 1; g.stroke(); }
      g.restore();
    }
    text(s, x, y, o) {
      o = o || {};
      if (/[_^]\{/.test(s)) { this.rich(s, x, y, o); return; }
      const g = this.g;
      g.save();
      g.font = this.font(o.size || 11, o.weight);
      g.textAlign = o.align || "left"; g.textBaseline = o.base || "middle";
      if (o.rot) { g.translate(x, y); g.rotate(o.rot); x = 0; y = 0; }
      if (o.halo !== false) { g.lineWidth = 3; g.strokeStyle = o.haloColor || this.c.bg; g.lineJoin = "round"; g.strokeText(s, x, y); }
      g.fillStyle = o.color || this.c.ink;
      g.fillText(s, x, y);
      g.restore();
    }
    /** Text with _{subscript} and ^{superscript} parts. */
    rich(s, x, y, o) {
      const size = o.size || 11, sub = Math.round(size * 0.78), segs = [];
      const re = /([_^])\{([^}]*)\}/g;
      let last = 0, m;
      while ((m = re.exec(s))) {
        if (m.index > last) segs.push([s.slice(last, m.index), 0]);
        segs.push([m[2], m[1] === "_" ? 1 : -1]);
        last = re.lastIndex;
      }
      if (last < s.length) segs.push([s.slice(last), 0]);
      const widths = segs.map(([t, k]) => this.measure(t, k ? sub : size, o.weight));
      const total = widths.reduce((a, b) => a + b, 0);
      let x0 = o.align === "center" ? x - total / 2 : o.align === "right" ? x - total : x;
      segs.forEach(([t, k], i) => {
        this.text(t, x0, y + k * size * 0.3, Object.assign({}, o, {size: k ? sub : size, align: "left"}));
        x0 += widths[i];
      });
    }
    measure(s, size, weight) { this.g.save(); this.g.font = this.font(size, weight); const w = this.g.measureText(s).width; this.g.restore(); return w; }
    arrow(x1, y1, x2, y2, color, w, head) {
      const L = Math.hypot(x2 - x1, y2 - y1);
      if (L < 0.5) return;
      head = Math.min(head || 7, L * 0.6);
      const ux = (x2 - x1) / L, uy = (y2 - y1) / L;
      this.line([[x1, y1], [x2 - ux * head * 0.6, y2 - uy * head * 0.6]], color, w || 1.5);
      this.poly([[x2, y2], [x2 - ux * head - uy * head * 0.45, y2 - uy * head + ux * head * 0.45],
        [x2 - ux * head + uy * head * 0.45, y2 - uy * head - ux * head * 0.45]], color);
    }
    /** Zig-zag spring from (x1, y1) to (x2, y2). */
    spring(x1, y1, x2, y2, color, coils, width) {
      const L = Math.hypot(x2 - x1, y2 - y1), ux = (x2 - x1) / L, uy = (y2 - y1) / L, nx = -uy, ny = ux;
      const n = (coils || 4) * 2, lead = Math.min(6, L * 0.15), wz = width || 5, pts = [[x1, y1], [x1 + ux * lead, y1 + uy * lead]];
      for (let i = 1; i < n; i++) {
        const s = lead + (L - 2 * lead) * i / n, sgn = i % 2 ? 1 : -1;
        pts.push([x1 + ux * s + nx * wz * sgn, y1 + uy * s + ny * wz * sgn]);
      }
      pts.push([x2 - ux * lead, y2 - uy * lead], [x2, y2]);
      this.line(pts, color, 1.3);
    }
    /** Dashpot (piston in a cylinder) from (x1, y1) to (x2, y2). */
    dashpot(x1, y1, x2, y2, color, width) {
      const L = Math.hypot(x2 - x1, y2 - y1), ux = (x2 - x1) / L, uy = (y2 - y1) / L, nx = -uy, ny = ux, w = width || 6;
      const a = L * 0.32, p = L * 0.55, b = L * 0.75;
      const P = (s, o) => [x1 + ux * s + nx * o, y1 + uy * s + ny * o];
      this.line([P(0, 0), P(p, 0)], color, 1.3);                          // piston rod
      this.line([P(p, -w * 0.7), P(p, w * 0.7)], color, 1.8);             // piston
      this.line([P(a, -w), P(b, -w), P(b, w), P(a, w)], color, 1.3);      // cylinder, open towards the rod
      this.line([P(b, 0), [x2, y2]], color, 1.3);
    }
    /** Hatched ground strip (x, y) .. (x + w), hatches below y. */
    ground(x, y, w, color) {
      this.line([[x, y], [x + w, y]], color, 1.5);
      for (let s = x + 4; s < x + w; s += 7) this.line([[s, y], [s - 5, y + 6]], color, 1);
    }
  }

  // ---------------------------------------------------------------- charts
  function niceStep(span, n) {
    const raw = span / Math.max(1, n), p = Math.pow(10, Math.floor(Math.log10(raw))), r = raw / p;
    return (r < 1.5 ? 1 : r < 3.5 ? 2 : r < 7.5 ? 5 : 10) * p;
  }
  function linTicks(a, b, n) {
    const st = niceStep(b - a, n), out = [];
    for (let v = Math.ceil(a / st - 1e-9) * st; v <= b + 1e-9 * st; v += st) out.push(Math.abs(v) < 1e-12 * st ? 0 : v);
    return out;
  }
  function logTicks(a, b) {
    const out = [];
    for (let e = Math.floor(Math.log10(a)); e <= Math.ceil(Math.log10(b)); e++) {
      for (const m of [1, 2, 5]) { const v = m * Math.pow(10, e); if (v >= a * 0.999 && v <= b * 1.001) out.push(v); }
    }
    return out;
  }
  const tickLabel = (v) => (Math.abs(v) >= 1e5 || (Math.abs(v) < 0.01 && v !== 0) ? v.toExponential(0).replace("e+", "e") : String(+v.toPrecision(4)));
  /** Axes in box {x, y, w, h}: grid, ticks, labels; returns the scales {X, Y, ix, iy, box}. */
  function axes(p, box, o) {
    const xl = o.xlog, yl = o.ylog;
    const tx = (v) => (xl ? Math.log10(v) : v), ty = (v) => (yl ? Math.log10(v) : v);
    const X = (v) => box.x + (tx(v) - tx(o.x0)) / (tx(o.x1) - tx(o.x0)) * box.w;
    const Y = (v) => box.y + box.h - (ty(v) - ty(o.y0)) / (ty(o.y1) - ty(o.y0)) * box.h;
    const ix = (px) => { const t = tx(o.x0) + (px - box.x) / box.w * (tx(o.x1) - tx(o.x0)); return xl ? Math.pow(10, t) : t; };
    const iy = (py) => { const t = ty(o.y0) + (box.y + box.h - py) / box.h * (ty(o.y1) - ty(o.y0)); return yl ? Math.pow(10, t) : t; };
    const c = p.c;
    const xt = o.xt || (xl ? logTicks(o.x0, o.x1) : linTicks(o.x0, o.x1, Math.max(2, Math.floor(box.w / 60))));
    const yt = o.yt || (yl ? logTicks(o.y0, o.y1) : linTicks(o.y0, o.y1, Math.max(2, Math.floor(box.h / 34))));
    p.rect(box.x, box.y, box.w, box.h, c.bg);
    for (const v of xt) p.line([[X(v), box.y], [X(v), box.y + box.h]], c.grid, 1);
    for (const v of yt) p.line([[box.x, Y(v)], [box.x + box.w, Y(v)]], c.grid, 1);
    p.line([[box.x, box.y + box.h], [box.x + box.w, box.y + box.h]], c.axis, 1);
    p.line([[box.x, box.y], [box.x, box.y + box.h]], c.axis, 1);
    let lastx = -1e9;
    const decadesOnly = xl && !o.xt && box.w / Math.max(1, Math.log10(o.x1 / o.x0)) < 150;
    for (const v of xt) {
      if (decadesOnly && Math.abs(Math.log10(v) - Math.round(Math.log10(v))) > 1e-9) continue;
      const s = tickLabel(v), x = X(v), w = p.measure(s, 10);
      if (x - w / 2 < lastx + 4) continue;
      p.text(s, x, box.y + box.h + 9, {size: 10, color: c.muted, align: "center", halo: false});
      lastx = x + w / 2;
    }
    for (const v of yt) p.text(tickLabel(v), box.x - 4, Y(v), {size: 10, color: c.muted, align: "right", halo: false});
    if (o.xlabel) p.text(o.xlabel, box.x + box.w, box.y + box.h + 22, {size: 10.5, color: c.muted, align: "right", halo: false});
    if (o.ylabel) p.text(o.ylabel, box.x, box.y - 8, {size: 10.5, color: c.muted, align: "left", halo: false});
    if (o.title) p.text(o.title, box.x + box.w, box.y - 8, {size: 11, weight: 600, color: c.ink, align: "right", halo: false});
    return {X, Y, ix, iy, box, o};
  }
  /** Polyline of (xs, ys) through the scales, clipped to the box. */
  function plot(p, ax, xs, ys, color, w, dash) {
    const g = p.g, b = ax.box;
    g.save(); g.beginPath(); g.rect(b.x - 1, b.y - 1, b.w + 2, b.h + 2); g.clip();
    const pts = [];
    for (let i = 0; i < xs.length; i++) if (Number.isFinite(ys[i])) pts.push([ax.X(xs[i]), ax.Y(ys[i])]);
    p.line(pts, color, w || 2, dash);
    g.restore();
  }
  const inBox = (b, x, y) => x >= b.x && x <= b.x + b.w && y >= b.y && y <= b.y + b.h;

  // ---------------------------------------------------------------- the figure frame
  /** Mount figure `name` into host.  opts: {params, caption (node), key}. */
  function mount(host, name, opts) {
    opts = opts || {};
    const def = FIG[name];
    if (!def) { host.append(el("p", {class: "lfig-missing", text: `(unknown figure "${name}")`})); return null; }
    const key = opts.key || name;
    const saved = STATE.get(key);
    const st = {};
    for (const [k, d] of Object.entries(def.params || {})) {
      let v = d.def;
      const given = opts.params && opts.params[k];
      if (given !== undefined) v = d.values ? String(given) : clamp(Number(given), d.min, d.max);
      if (d.values && !d.values.includes(v)) v = d.def;
      if (!d.values && !Number.isFinite(v)) v = d.def;
      st[k] = v;
    }
    if (saved) Object.assign(st, saved.st);
    const canvas = el("canvas", {class: "lfig-canvas", role: "img", "aria-label": def.alt || def.title});
    const stage = el("div", {class: "lfig-stage"}, canvas);
    const playBtn = def.animated ? el("button", {class: "btn small lfig-play", type: "button"}) : null;
    const head = el("div", {class: "lfig-head"}, el("span", {class: "lfig-title", text: def.title}), el("span", {class: "lfig-grow"}), playBtn);
    const ctrls = el("div", {class: "lfig-ctrls"});
    const read = el("div", {class: "lfig-read", "aria-live": "polite"});
    const tex = (def.tex || []).length ? el("div", {class: "lfig-tex"}, ...def.tex.map((t) => el("span", {class: "lfig-f", "data-tex": t}))) : null;
    const note = def.note ? el("div", {class: "lfig-note"}, richSpan(def.note)) : null;
    const root = el("figure", {class: "lfig", "data-figure": name}, head, stage, ctrls, read, tex, note, opts.caption || null);
    host.replaceChildren(root);
    const fig = {def, name, root, stage, canvas, st, key, t: saved ? saved.t : 0, W: 0, H: 0, visible: true, dirty: true,
      playing: def.animated ? (saved ? saved.playing : !reducedMotion()) : false, controls: {}, regions: [], data: {}};
    fig.col = palette(root);
    fig.save = () => STATE.set(key, {st: Object.assign({}, fig.st), t: fig.t, playing: fig.playing});
    fig.read = (parts) => { read.replaceChildren(...parts.map((s) => (typeof s === "string" ? richSpan(s) : s))); };
    fig.redraw = () => { fig.dirty = true; if (!raf && fig.inited) draw(); };
    fig.set = (k, v, from) => {
      fig.st[k] = v;
      const c = fig.controls[k];
      if (c && from !== c) c.sync();
      if (def.change) def.change(fig, k);
      fig.save();
      fig.redraw();
    };
    // ---- controls
    fig.slider = function (k, label, o) {
      o = o || {};
      const d = def.params[k];
      const log = !!o.log, N = 1000;
      const toPos = (v) => (log ? Math.log(v / d.min) / Math.log(d.max / d.min) : (v - d.min) / (d.max - d.min)) * N;
      const toVal = (pos) => { const t = pos / N; let v = log ? d.min * Math.pow(d.max / d.min, t) : d.min + t * (d.max - d.min); if (d.step) v = Math.round(v / d.step) * d.step; return clamp(v, d.min, d.max); };
      const inp = el("input", {type: "range", min: 0, max: N, step: 1, "aria-label": label.replace(/[_^]\{([^}]*)\}/g, "$1")});
      const out = el("output", {class: "lfig-val"});
      const c = {el: inp, sync() { inp.value = String(Math.round(toPos(fig.st[k]))); out.textContent = (o.fmt || ((v) => fmt(v)))(fig.st[k]); }};
      inp.addEventListener("input", () => { fig.set(k, toVal(Number(inp.value)), c); out.textContent = (o.fmt || ((v) => fmt(v)))(fig.st[k]); });
      fig.controls[k] = c;
      c.sync();
      ctrls.append(el("label", {class: "lfig-ctl"}, el("span", {class: "lfig-lab"}, richSpan(label)), inp, out));
      return c;
    };
    fig.choice = function (k, label, options) {
      const btns = options.map(([v, text, tip]) => el("button", {class: "lfig-seg", type: "button", title: tip || null, text, onclick: () => fig.set(k, v)}));
      const c = {sync() { btns.forEach((b, i) => { const on = options[i][0] === fig.st[k]; b.classList.toggle("on", on); b.setAttribute("aria-pressed", on ? "true" : "false"); }); }};
      fig.controls[k] = c;
      c.sync();
      ctrls.append(el("div", {class: "lfig-ctl lfig-choice", role: "group", "aria-label": label}, label ? el("span", {class: "lfig-lab", text: label}) : null, ...btns));
      return c;
    };
    fig.button = function (text, fn, tip) { const b = el("button", {class: "btn small", type: "button", title: tip || null, text, onclick: fn}); ctrls.append(b); return b; };
    /** Pointer on a chart box sets the value of parameter k (drag or click). */
    fig.drag = (boxFn, fn) => fig.regions.push({boxFn, fn});
    // ---- play / pause
    const paintPlay = () => {
      if (!playBtn) return;
      playBtn.textContent = fig.playing ? "❚❚ Pause" : "▶ Play";
      playBtn.setAttribute("aria-pressed", fig.playing ? "true" : "false");
      playBtn.title = fig.playing ? "pause the animation" : "play the animation";
    };
    if (playBtn) playBtn.addEventListener("click", () => { fig.playing = !fig.playing; paintPlay(); fig.save(); kick(); });
    paintPlay();
    // ---- layout and drawing
    let raf = 0, last = 0, observer = null, resizer = null;
    function layout() {
      const W = Math.max(260, Math.floor(stage.clientWidth || root.clientWidth || 600));
      if (W === fig.W && fig.H) return false;
      fig.W = W;
      fig.H = def.layout(fig, W);
      const dpr = Math.min(2.5, (typeof devicePixelRatio === "number" && devicePixelRatio) || 1);
      canvas.width = Math.round(W * dpr); canvas.height = Math.round(fig.H * dpr);
      canvas.style.width = W + "px"; canvas.style.height = fig.H + "px";
      fig.g = canvas.getContext("2d");
      fig.g.setTransform(dpr, 0, 0, dpr, 0, 0);
      fig.pen = new Pen(fig.g, fig.col);
      return true;
    }
    function draw() {
      if (!fig.g) layout();
      if (!fig.g) return;
      fig.g.clearRect(0, 0, fig.W, fig.H);
      fig.pen.rect(0, 0, fig.W, fig.H, fig.col.bg);
      def.draw(fig, fig.pen);
      fig.dirty = false;
    }
    function destroy() {
      if (raf) cancelAnimationFrame(raf);
      raf = 0;
      if (observer) observer.disconnect();
      if (resizer) resizer.disconnect();
      fig.dead = true;
    }
    function frame(now) {
      raf = 0;
      if (!root.isConnected) { destroy(); return; }
      const dt = Math.min(0.1, Math.max(0, (now - last) / 1000));
      last = now;
      if (fig.playing && fig.visible) {
        fig.t += dt;
        if (def.tick) def.tick(fig, dt);
        draw();
        raf = requestAnimationFrame(frame);
      } else if (fig.dirty) draw();
    }
    function kick() {
      if (fig.dead || !fig.inited) return;
      if (fig.playing && fig.visible && !raf) { last = performance.now(); raf = requestAnimationFrame(frame); }
      else if (!fig.playing) { fig.save(); draw(); }
    }
    fig.kick = kick;
    // pointer regions (charts that set a parameter)
    const pointer = (ev) => {
      const r = canvas.getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top;
      for (const rg of fig.regions) {
        const b = rg.boxFn();
        if (b && inBox(b, x, y)) { rg.fn(x, y, ev); ev.preventDefault(); return true; }
      }
      return false;
    };
    let dragging = false;
    canvas.addEventListener("pointerdown", (ev) => { if (pointer(ev)) { dragging = true; canvas.setPointerCapture(ev.pointerId); } });
    canvas.addEventListener("pointermove", (ev) => { if (dragging) pointer(ev); });
    canvas.addEventListener("pointerup", () => { dragging = false; });
    canvas.addEventListener("pointercancel", () => { dragging = false; });
    // the definition builds its controls and state
    if (def.init) def.init(fig);
    if (typeof IntersectionObserver === "function") {
      observer = new IntersectionObserver((entries) => {
        for (const e of entries) fig.visible = e.isIntersecting;
        if (!root.isConnected) { destroy(); return; }
        kick();
      }, {threshold: 0.05});
      observer.observe(root);
    }
    if (typeof ResizeObserver === "function") {
      resizer = new ResizeObserver(() => {
        if (!root.isConnected) { destroy(); return; }
        if (!fig.inited) { if (stage.clientWidth > 0) first(); return; }
        if (layout()) { if (def.change) def.change(fig, "_layout"); draw(); }
      });
      resizer.observe(stage);
    }
    // first layout as soon as the stage has a width: the lesson page is attached right after fill(), so
    // try again after the current task (timers also run in a hidden tab, animation frames do not)
    const first = () => { if (fig.inited || fig.dead) return; fig.inited = true; layout(); if (def.change) def.change(fig, "_init"); draw(); kick(); };
    const tryFirst = (n) => {
      if (fig.inited || fig.dead) return;
      if (stage.clientWidth > 0) first();
      else if (!root.isConnected && n > 40) destroy();
      else if (n < 60) setTimeout(() => tryFirst(n + 1), n < 5 ? 0 : 100);   // later: the ResizeObserver
    };
    Promise.resolve().then(() => tryFirst(0));
    root._fig = fig;
    return fig;
  }

  /** Mount the figures of `blocks` (server-rendered lesson blocks of kind "figure") into the placeholders
   *  <div class="lesson-figure" data-figure-block="k"> under root.  caption(b) returns the caption node. */
  function fill(root, blocks, caption, keyPrefix) {
    const byK = {};
    for (const b of blocks || []) if (b.kind === "figure") byK[b.k] = b;
    root.querySelectorAll(".lesson-figure[data-figure-block]").forEach((ph) => {
      const b = byK[Number(ph.dataset.figureBlock)];
      if (!b) { ph.remove(); return; }
      const cap = b.caption_html && caption ? caption(b) : null;
      const wrap = cap ? el("figcaption", {class: "lfig-cap"}, cap) : null;
      mount(ph, b.figure, {params: b.params || {}, caption: wrap, key: `${keyPrefix || ""}/${b.k}/${b.figure}`});
      ph.classList.add("mounted");
    });
  }

  // ================================================================== shared drawing bits
  const VIS = 0.6;                         // animation cycles per second (motion is shown slowed down)
  /** Wrap text into lines no wider than maxw. */
  function wrap(p, s, maxw, size, weight) {
    const words = String(s).split(/\s+/), out = [];
    let cur = "";
    for (const w of words) {
      const t = cur ? cur + " " + w : w;
      if (cur && p.measure(t, size, weight) > maxw) { out.push(cur); cur = w; } else cur = t;
    }
    if (cur) out.push(cur);
    return out;
  }
  function title(p, s, x, y, maxw, o) {
    o = o || {};
    const lines = wrap(p, s, maxw, o.size || 11.5, o.weight || 600);
    lines.forEach((ln, i) => p.text(ln, x, y + i * ((o.size || 11.5) + 3), {size: o.size || 11.5, weight: o.weight || 600, color: o.color || p.c.ink, align: o.align || "left", base: "top", halo: false}));
    return lines.length * ((o.size || 11.5) + 3);
  }
  /** Panel frame (light border) with optional highlight. */
  function frame(p, b, on) {
    p.rect(b.x + 0.5, b.y + 0.5, b.w - 1, b.h - 1, null, on ? p.c.s1 : p.c.grid, on ? 2 : 1);
  }
  /** Soil block with layer bands; free-field lines displaced by disp(zFrac) (px). */
  function soilBlock(p, x, y, w, h, bands, disp, lines) {
    let yy = y;
    for (const bd of bands) { p.rect(x, yy, w, h * bd.f, bd.color); yy += h * bd.f; }
    const n = lines || Math.max(4, Math.round(w / 16));
    const g = p.g;
    g.save(); g.beginPath(); g.rect(x, y, w, h); g.clip();
    for (let i = 0; i <= n; i++) {
      const x0 = x + (i + 0.5) * w / (n + 1), pts = [];
      for (let k = 0; k <= 16; k++) { const zf = k / 16; pts.push([x0 + (disp ? disp(zf) : 0), y + zf * h]); }
      p.line(pts, withAlpha(p.c.soilLine, 0.55), 1);
    }
    g.restore();
    p.line([[x, y], [x + w, y]], p.c.soilLine, 1.5);
  }
  /** A stick with a lumped mass: base (bx, by), height hpx, lateral displacement d(zFrac) in px. */
  function stick(p, bx, by, hpx, d, o) {
    o = o || {};
    const pts = [];
    for (let k = 0; k <= 12; k++) { const zf = k / 12; pts.push([bx + d(zf), by - zf * hpx]); }
    p.line(pts, o.color || p.c.struct, o.w || 3, o.dash);
    const top = pts[pts.length - 1];
    p.circle(top[0], top[1], o.r || 8, o.fill === undefined ? p.c.struct : o.fill, o.ring || null, 1.5);
    return top;
  }
  const shapeCant = (zf) => zf * zf * (3 - zf) / 2;          // a cantilever-like deflected shape

  // ================================================================== figure: fixed base vs SSI
  FIG["fixed-base-vs-ssi"] = {
    title: "Fixed base or on the soil: the same stick, two systems",
    alt: "Animated comparison of a stick with a lumped mass on a rigid base and on sway and rocking springs and dashpots, with the amplitude of the mass motion against frequency for both.",
    animated: true,
    params: {f: {def: 3.5, min: 0.2, max: 10, step: 0.01}, s: {def: 1, min: 0.5, max: 4, step: 0.01}},
    tex: ["k^* = \\omega_1^2 M_1\\,c(\\beta),\\qquad S_x = K_x\\,(1 + 2i\\beta_g) + i\\omega c_x,\\qquad S_\\theta = K_\\theta\\,(1 + 2i\\beta_g) + i\\omega c_\\theta",
      "\\left(\\tilde f / f_1\\right)^2 = 1\\big/\\left(1 + k/K_x + k\\,h^2/K_\\theta\\right)"],
    note: "Computed in closed form: the steady state of the three-degree-of-freedom system (stick deformation, mat sway, rocking) per unit free-field motion. Mode 1 of the example 1 stick (M₁ = 3063 t at h₁ = 15.6 m, 5.05 Hz, β = 5 %), the 1468 t mat, and SASSI's impedance of this mat at 3.49 Hz (lesson 4: K_{x} = 1.48e7 kN/m, K_{θ} = 1.70e9 kN m/rad, 4.8 % material damping β_{g}, radiation dashpots for 7.7 % sliding and 4.9 % rocking) held constant with frequency. The second formula is the Veletsos-Meek estimate without the mat mass. Motion slowed down; wave fronts scaled by the power the dashpots radiate.",
    init(fig) {
      fig.slider("f", "Frequency f", {fmt: (v) => fmt(v, 2) + " Hz"});
      fig.slider("s", "Soil velocities × s", {log: true, fmt: (v) => `× ${fmt(v, 2)} (sand Vs ${Math.round(300 * v)} m/s)`});
      fig.drag(() => fig.data.ax && fig.data.ax.box, (x) => fig.set("f", clamp(+fig.data.ax.ix(x).toFixed(2), 0.2, 10)));
      fig.data.pref = null;
    },
    change(fig, k) {
      const d = fig.data;
      if (!d.curve || k === "s" || k === "_init") {
        const m = PH.ssiModel(fig.st.s);
        const fs = [], ht = [], hb = [];
        for (let f = 0.1; f <= 10.0001; f += 0.005) { const r = PH.ssiResponse(m, f); fs.push(f); ht.push(C.abs(r.Ht)); hb.push(C.abs(r.Hfb)); }
        d.model = m; d.curve = {fs, ht, hb}; d.pk = PH.peakInfo(fs, ht); d.pkfb = PH.peakInfo(fs, hb);
        if (!d.pref) { const m1 = PH.ssiModel(1), pk = PH.peakInfo(fs, ht); d.pref = PH.ssiResponse(m1, pk.f).prad; }
      }
      const m = d.model, f = fig.st.f, r = PH.ssiResponse(m, f), w = TAU * f;
      d.r = r;
      fig.read([`f = ${fmt(f, 2)} Hz`, `|H| fixed base ${fmt(C.abs(r.Hfb), 2)}, on soil ${fmt(C.abs(r.Ht), 2)}`,
        `peak on soil ${fmt(d.pk.f, 2)} Hz (fixed base ${fmt(d.pkfb.f, 2)} Hz), system damping ≈ ${fmt(100 * d.pk.zeta, 1)} % (half-power)`,
        `foundation radiation damping at f: sliding ${fmt(100 * w * m.cx / (2 * m.Kx), 1)} %, rocking ${fmt(100 * w * m.ct / (2 * m.Kt), 1)} %`]);
    },
    layout(fig, W) {
      const pad = 8, sh = W >= 600 ? 250 : 230, ch = W >= 600 ? 180 : 165;
      const sw = (W - 3 * pad) / 2;
      fig.L = {a: {x: pad, y: pad, w: sw, h: sh}, b: {x: 2 * pad + sw, y: pad, w: sw, h: sh},
        chart: {x: 44, y: sh + pad + 26, w: W - 44 - 14, h: ch - 26 - 30}};
      return sh + ch + pad;
    },
    draw(fig, p) {
      const L = fig.L, d = fig.data, c = p.c, r = d.r;
      if (!r) return;
      const phi = TAU * VIS * fig.t, e = [Math.cos(phi), Math.sin(phi)];
      const re = (z) => z[0] * e[0] - z[1] * e[1];                   // Re[z e^{i phi}]
      const ug = e[0];
      const maxH = Math.max(...d.curve.ht, ...d.curve.hb);
      const kap = Math.min(3.2, 0.3 * L.a.w / maxH);
      // ---- fixed base
      const A = L.a, gy = A.y + Math.round(0.6 * A.h), hp = gy - A.y - 58, cxa = A.x + A.w / 2;
      frame(p, A, false);
      title(p, A.w > 300 ? "Fixed base: the base moves with the ground" : "Fixed base", A.x + 8, A.y + 7, A.w - 16, {size: 11});
      const db = kap * ug;
      p.rect(A.x + 10 + db, gy, A.w - 20, 16, c.rock);
      p.ground(A.x + 10 + db, gy, A.w - 20, c.struct);
      const dm = kap * re(r.Hfb);
      const topA = stick(p, cxa + db, gy, hp, (zf) => (dm - db) * shapeCant(zf));
      p.text("M₁", topA[0] + 12, topA[1] - 2, {size: 10.5, color: c.muted});
      p.text(`|H| = ${fmt(C.abs(r.Hfb), 1)}`, A.x + 10, gy + 28, {size: 10.5, color: c.s2, weight: 600});
      // ---- on soil
      const B = L.b, gyb = B.y + Math.round(0.6 * B.h), cxb = B.x + B.w / 2;
      frame(p, B, false);
      const thB = title(p, B.w > 300 ? "On the soil: sway and rocking springs and dashpots" : "On soil springs and dashpots", B.x + 8, B.y + 7, B.w - 16, {size: 11});
      const soilTop = gyb, soilH = B.y + B.h - 1 - soilTop;
      soilBlock(p, B.x + 1, soilTop, B.w - 2, soilH, [{f: 1, color: c.soil1}], () => kap * ug, Math.round(B.w / 18));
      // radiated wave fronts (lower half circles), strength ~ sqrt(radiated power)
      const strength = clamp(Math.sqrt(r.prad / d.pref), 0, 3);
      const g = p.g;
      g.save(); g.beginPath(); g.rect(B.x + 1, soilTop, B.w - 2, soilH); g.clip();
      const per = 1 / VIS, speed = 34, rmax = Math.hypot(B.w / 2, soilH);
      for (let k = 0; k < 8; k++) {
        const age = ((fig.t % per) + k * per);
        const rad = 18 + speed * age;
        if (rad > rmax) break;
        const a = clamp(0.55 * strength * (1 - rad / rmax), 0, 0.9);
        if (a < 0.02) continue;
        g.beginPath(); g.arc(cxb, soilTop, rad, 0.05, Math.PI - 0.05);
        g.strokeStyle = withAlpha(c.s1, a); g.lineWidth = 1.6 + strength; g.stroke();
      }
      g.restore();
      // mat: sway u0 and rocking theta (rigid rotation drawn consistently with the stick)
      const hpb = hp - 7, mw = Math.min(0.4 * B.w, 110), mt = 8;
      const d0 = kap * (ug + re(r.u0));
      const rot = kap * re(r.th) * PH.SSI.h / hpb;                   // visual rotation (rad)
      const mx = cxb + d0;
      p.poly([[mx - mw / 2, gyb - mt - rot * mw / 2], [mx + mw / 2, gyb - mt + rot * mw / 2], [mx + mw / 2, gyb + rot * mw / 2], [mx - mw / 2, gyb - rot * mw / 2]], c.structFill, c.struct, 1.2);
      // rocking springs and dashpots under the mat edges, sway pair at the left
      const anchor = kap * ug;
      for (const sx of [-0.36, 0.36]) {
        const ex = mx + sx * mw, ey = gyb + rot * sx * mw;
        p.spring(ex - 7, ey + 1, cxb + sx * mw - 7 + anchor, gyb + 40, c.struct, 4, 4);
        p.dashpot(ex + 7, ey + 1, cxb + sx * mw + 7 + anchor, gyb + 40, c.struct, 4.5);
        p.ground(cxb + sx * mw - 15 + anchor, gyb + 40, 30, c.struct);
      }
      const ax0 = cxb - mw / 2 - Math.min(52, (B.w - mw) / 2 - 12) + anchor, mlx = mx - mw / 2;
      p.line([[mlx, gyb - mt - rot * mw / 2], [mlx, gyb + 11]], c.struct, 2);          // connector on the mat side
      p.spring(ax0, gyb - 4, mlx, gyb - 4, c.struct, 4, 4);
      p.dashpot(ax0, gyb + 8, mlx, gyb + 8, c.struct, 4.5);
      p.line([[ax0, gyb - 12], [ax0, gyb + 15]], c.struct, 2);
      p.text("K_{x}, c_{x}", ax0 - 2, gyb - 20, {size: 10.5, color: c.muted});
      p.text("K_{θ}, c_{θ}", cxb + anchor, gyb + 54, {size: 10.5, color: c.muted, align: "center"});
      // stick: rigid part (sway + rocking) plus its own deformation
      const du = kap * re(r.u);
      const topB = stick(p, mx, gyb - mt, hpb, (zf) => rot * zf * hpb + du * shapeCant(zf));
      p.text("M₁", topB[0] + 12, topB[1] - 2, {size: 10.5, color: c.muted});
      p.text(`|H| = ${fmt(C.abs(r.Ht), 1)}`, B.x + 10, B.y + thB + 20, {size: 10.5, color: c.s1, weight: 600});
      // ---- chart
      const cv = d.curve, ymax = Math.max(12, Math.ceil(maxH * 1.12));
      const ax = axes(p, L.chart, {x0: 0, x1: 10, y0: 0, y1: ymax, xlabel: "f (Hz)", ylabel: "|H| = mass motion / free-field motion"});
      d.ax = ax;
      plot(p, ax, cv.fs, cv.hb, c.s2, 2);
      plot(p, ax, cv.fs, cv.ht, c.s1, 2);
      p.text("fixed base", ax.X(d.pkfb.f) + 6, ax.Y(d.pkfb.amp) - 2, {size: 10.5, color: c.ink});
      p.text("on soil", ax.X(d.pk.f) - 6, ax.Y(d.pk.amp) - 2, {size: 10.5, color: c.ink, align: "right"});
      const f = fig.st.f, X = ax.X(f);
      p.line([[X, ax.box.y], [X, ax.box.y + ax.box.h]], c.ink, 1, [3, 3]);
      p.circle(X, ax.Y(C.abs(r.Hfb)), 4.5, c.s2, c.bg, 2);
      p.circle(X, ax.Y(C.abs(r.Ht)), 4.5, c.s1, c.bg, 2);
      if (ax.box.w > 420) p.text("drag on the plot to change f", ax.box.x + ax.box.w, ax.box.y + 8, {size: 9.5, color: c.muted, align: "right"});
    },
  };

  // ================================================================== figure: kinematic + inertial
  /** Free field of a uniform soil column per unit surface motion: cos(2 pi z / lambda) (undamped standing
   *  wave of a vertically incident wave in a half-space with a free surface). */
  const ffUniform = (zOverLambda) => Math.cos(TAU * zOverLambda);
  /** Rigid embedded box following the free field over its depth (schematic: least-squares translation at
   *  mid-depth and rotation); D/lambda = r.  Returns {u (at mid-depth), b (per unit depth fraction)}. */
  function boxAverage(r) {
    let a = 0, b = 0;
    const n = 64;
    for (let i = 0; i < n; i++) { const zf = (i + 0.5) / n, u = ffUniform(r * zf); a += u / n; b += 12 * (zf - 0.5) * u / n; }
    return {u: a, b};
  }
  FIG["kinematic-inertial"] = {
    title: "SSI = kinematic interaction + inertial interaction",
    alt: "Schematic animation: the SSI response equals the kinematic response of a massless foundation and structure plus the inertial response of the masses on the soil springs and dashpots driven by the foundation input motion.",
    animated: true,
    params: {case: {def: "embedded", values: ["surface", "embedded"]}, r: {def: 0.15, min: 0, max: 0.3, step: 0.005}},
    tex: ["\\text{1: } \\left[K^*_s - C^e + X_{ff}\\right] U_k = X_{ff}\\,U'_f",
      "\\text{2: } \\left[K^*_s - \\omega^2 M_s - C^e + X_{ff}\\right] U_i = \\omega^2 M_s\\,U_k",
      "U = U_k + U_i"],
    note: "Schematic, amplitudes illustrative. C^{e} = K*_{e} − ω²M_{e} is the excavated soil (zero for a surface mat). Adding equations 1 and 2 gives the full SSI equation, so the split is exact for a linear system. The embedded box is drawn following an average of the free field cos(2πz/λ) over its depth; lesson 5 computes the real foundation input motion with a massless basement.",
    init(fig) {
      fig.choice("case", "Foundation", [["surface", "surface mat"], ["embedded", "embedded box"]]);
      fig.slider("r", "Embedment / wavelength D/λ", {fmt: (v) => fmt(v, 3)});
    },
    change(fig) {
      const emb = fig.st.case === "embedded";
      const k = emb ? boxAverage(fig.st.r) : {u: 1, b: 0};
      fig.data.k = k;
      // inertial part: the stick of figure fixed-base-vs-ssi at 0.85 of its SSI frequency (illustrative)
      if (!fig.data.inr) { const m = PH.ssiModel(1); fig.data.inr = PH.ssiResponse(m, 0.85 * 3.905); }
      fig.read(emb ? [`D/λ = ${fmt(fig.st.r, 3)}`, `foundation input motion: translation ${fmt(k.u, 2)} × free-field surface motion (mid-depth), plus rocking`]
        : ["surface mat, vertically incident waves: the foundation input motion is the free-field motion (lesson 1, massless run), so all of the SSI effect is inertial"]);
    },
    layout(fig, W) {
      const pad = 8;
      if (W >= 640) {
        const ow = 22, pw = (W - 2 * pad - 2 * ow) / 3, h = 240;
        fig.L = {wide: true, p: [0, 1, 2].map((i) => ({x: pad + i * (pw + ow), y: pad, w: pw, h})), ops: [{x: pad + pw + ow / 2, y: pad + h / 2, s: "="}, {x: pad + 2 * pw + 1.5 * ow, y: pad + h / 2, s: "+"}]};
        return h + 2 * pad;
      }
      const h = 200, oh = 20;
      fig.L = {wide: false, p: [0, 1, 2].map((i) => ({x: pad, y: pad + i * (h + oh), w: W - 2 * pad, h})), ops: [{x: W / 2, y: pad + h + oh / 2, s: "="}, {x: W / 2, y: pad + 2 * h + 1.5 * oh, s: "+"}]};
      return 3 * h + 2 * oh + 2 * pad;
    },
    draw(fig, p) {
      const L = fig.L, c = p.c, emb = fig.st.case === "embedded", r = emb ? fig.st.r : 0, k = fig.data.k, ri = fig.data.inr;
      if (!k || !ri) return;
      const phi = TAU * VIS * fig.t, cs = Math.cos(phi), e = [cs, Math.sin(phi)];
      const re = (z) => z[0] * e[0] - z[1] * e[1];
      const kap = 9;
      const titles = ["SSI analysis: one run, total motion", "1  Kinematic: massless foundation and structure", "2  Inertial: the masses on the soil springs and dashpots, shaken by the foundation input motion"];
      for (const o of L.ops) p.text(o.s, o.x, o.y, {size: 22, weight: 600, color: c.muted, align: "center", halo: false});
      L.p.forEach((b, i) => {
        frame(p, b, false);
        const th = title(p, titles[i], b.x + 8, b.y + 7, b.w - 16, {size: 11});
        const gy = b.y + Math.max(th + 18, b.h * 0.5), soilH = b.y + b.h - 1 - gy, cx = b.x + b.w / 2;
        const D = emb ? Math.min(0.45 * soilH, 46) : 0, bw = Math.min(0.42 * b.w, 92);
        // free field per unit surface motion at depth zpx (px below grade): cos(2 pi z / lambda), z/lambda = r zpx/D
        const ff = (zpx) => (emb ? ffUniform(r * zpx / D) : 1);
        if (i < 2) soilBlock(p, b.x + 1, gy, b.w - 2, soilH, [{f: 1, color: c.soil1}], (zf) => kap * ff(zf * soilH) * cs, Math.round(b.w / 18));
        else p.rect(b.x + 1, gy, b.w - 2, soilH, withAlpha(c.soil1, 0.45));
        if (i < 2) {           // the vertically incident wave
          const wx = b.x + b.w - 22;
          p.arrow(wx, b.y + b.h - 22, wx, b.y + b.h - 22 - Math.min(40, soilH - 30), c.s1, 1.6, 6);
          p.text("wave", wx - 6, b.y + b.h - 12, {size: 9.5, color: c.muted, align: "right"});
        }
        // kinematic (foundation input motion): lateral motion at grade and rotation (px per px of height)
        const xgK = emb ? kap * cs * (k.u - 0.5 * k.b) : kap * cs;
        const rotK = emb ? -0.4 * kap * cs * k.b / D : 0;                     // drawn at 0.4 scale (schematic)
        // inertial (panels 0 and 2): mat sway, rocking and stick deformation driven by the input motion
        const hs = Math.min(gy - b.y - th - 26, 92);
        const xgI = i === 1 ? 0 : kap * k.u * re(ri.u0);
        // the rocking is drawn at 0.4 of its size; the rest goes into the stick so that its top moves right
        const rotFull = i === 1 ? 0 : kap * k.u * re(ri.th) * PH.SSI.h / Math.max(hs, 1);
        const rotI = 0.4 * rotFull;
        const defI = i === 1 ? 0 : kap * k.u * re(ri.u) + 0.6 * rotFull * hs;
        const xg = cx + xgK + xgI, rot = rotK + rotI;
        const fill = i === 1 ? withAlpha(c.structFill, 0.6) : c.structFill;
        // rigid (small) rotation: lateral motion xg - rot * depth, vertical motion rot * horizontal offset
        if (emb) p.poly([[xg - bw / 2, gy - rot * bw / 2], [xg + bw / 2, gy + rot * bw / 2], [xg - rot * D + bw / 2, gy + D + rot * bw / 2], [xg - rot * D - bw / 2, gy + D - rot * bw / 2]], fill, c.struct, 1.3);
        else p.poly([[xg - bw / 2, gy - 7 + rot * (-bw / 2)], [xg + bw / 2, gy - 7 + rot * (bw / 2)], [xg + bw / 2, gy + rot * (bw / 2)], [xg - bw / 2, gy - rot * (bw / 2)]], fill, c.struct, 1.2);
        if (i === 2) {           // springs and dashpots whose far ends move with the foundation input motion
          const by = gy + D, ay = by + 36, anc = cx + xgK;
          for (const sx of [-0.32, 0.32]) {
            p.spring(xg - rot * D + sx * bw - 5, by + rot * sx * bw, anc + sx * bw - 5, ay, c.struct, 3, 3);
            p.dashpot(xg - rot * D + sx * bw + 5, by + rot * sx * bw, anc + sx * bw + 5, ay, c.struct, 3.5);
            p.ground(anc + sx * bw - 12, ay, 24, c.struct);
          }
        }
        const foot = ["free field: " + (emb ? "varies with depth" : "uniform over the mat"), "massless: no inertia forces", "base of the springs moves with U_{k}"][i];
        p.text(foot, b.x + 8, b.y + b.h - 10, {size: 10, color: c.muted});
        // superstructure
        stick(p, xg, gy - (emb ? 0 : 7), hs, (zf) => rot * zf * hs + defI * shapeCant(zf),
          {fill: i === 1 ? c.bg : c.struct, ring: i === 1 ? c.struct : null, dash: i === 1 ? [4, 3] : null, w: i === 1 ? 2 : 3, r: 7});
      });
    },
  };

  // ================================================================== figure: substructuring (Eq. 2.1)
  const SUB_TERMS = ["\\bigl[(K^*_s-\\omega^2 M_s)", "-\\,(K^*_e-\\omega^2 M_e)", "+\\,X_{ff}\\bigr]\\,U", "=\\,X_{ff}\\,U'_f"];
  const SUB_ON = [[0, 1, 2, 3], [2, 3], [0], [1], [0, 1, 2, 3]];
  const SUB_TEXT = {
    embedded: [
      "The SSI system: a building with a basement in a layered site. SASSI splits it into three parts and adds them up.",
      "Free field: the layered site without any excavation (SITE, POINT). At the interaction nodes it is represented by its dynamic stiffness X_ff (springs and dashpots) and its motion U'_f.",
      "+ Structure: the finite-element model of the building, basement included (HOUSE: K*_s, M_s), with no supports.",
      "− Excavated soil: an FE model of the soil the basement replaces, with the free-field properties. It is subtracted because the free field already contains it (HOUSE: K*_e, M_e).",
      "= SSI system: at every frequency ω, ANALYS solves the equation below for the total motion U of every node. With all excavated nodes as interaction nodes (FV) it is exact.",
    ],
    surface: [
      "The SSI system: a building on a surface mat. SASSI splits it into the free field and the structure.",
      "Free field: the layered site (SITE, POINT). At the interaction nodes under the mat it is represented by its dynamic stiffness X_ff and its motion U'_f.",
      "+ Structure: the finite-element model of the building and its mat (HOUSE: K*_s, M_s), with no supports.",
      "− Excavated soil: none for a surface mat (K*_e = M_e = 0).",
      "= SSI system: at every frequency ω, ANALYS solves the equation below for the total motion U of every node.",
    ],
  };
  FIG.substructuring = {
    title: "Substructuring: free field + structure − excavated soil = SSI system",
    alt: "Step-by-step build-up of the SASSI substructuring equation from the free field, the structure and the excavated soil.",
    animated: true,
    params: {case: {def: "embedded", values: ["surface", "embedded"]}, stage: {def: 0, min: 0, max: 4, step: 1}},
    texParts: SUB_TERMS,
    note: "Manual Eq. 2.1 (Theory Section 3). The highlighted terms belong to the highlighted part. Play steps through the parts; ◀ ▶ step by hand.",
    init(fig) {
      fig.choice("case", "Foundation", [["surface", "surface mat"], ["embedded", "embedded"]]);
      fig.button("◀", () => { fig.playing = false; fig.kick(); fig.set("stage", (fig.st.stage + 4) % 5); }, "previous part");
      fig.button("▶", () => { fig.playing = false; fig.kick(); fig.set("stage", (fig.st.stage + 1) % 5); }, "next part");
      const terms = SUB_TERMS.map((t) => el("span", {class: "lfig-eqt", "data-tex": t}));
      fig.data.terms = terms;
      fig.root.insertBefore(el("div", {class: "lfig-eq", "aria-label": "substructuring equation"}, ...terms), fig.root.querySelector(".lfig-read"));
      fig.data.acc = 0;
    },
    tick(fig, dt) {
      fig.data.acc += dt;
      if (fig.data.acc > 3.2) { fig.data.acc = 0; fig.st.stage = (fig.st.stage + 1) % 5; fig.save(); paintSub(fig); }
    },
    change(fig) { paintSub(fig); },
    layout(fig, W) {
      const pad = 8;
      if (W >= 700) {
        const ow = 20, pw = (W - 2 * pad - 3 * ow) / 4, h = 200;
        fig.L = {p: [0, 1, 2, 3].map((i) => ({x: pad + i * (pw + ow), y: pad, w: pw, h})),
          ops: ["+", "−", "="].map((s, i) => ({s, x: pad + (i + 1) * pw + i * ow + ow / 2, y: pad + h / 2}))};
        return h + 2 * pad;
      }
      const ow = 20, lm = 16, pw = (W - 2 * pad - ow - lm) / 2, h = 170, oh = 16, x0 = pad + lm;
      fig.L = {p: [0, 1, 2, 3].map((i) => ({x: x0 + (i % 2) * (pw + ow), y: pad + Math.floor(i / 2) * (h + oh), w: pw, h})),
        ops: [{s: "+", x: x0 + pw + ow / 2, y: pad + h / 2}, {s: "−", x: pad + lm / 2, y: pad + h + oh + h / 2}, {s: "=", x: x0 + pw + ow / 2, y: pad + h + oh + h / 2}]};
      return 2 * h + oh + 2 * pad;
    },
    draw(fig, p) {
      const c = p.c, emb = fig.st.case === "embedded", stg = fig.st.stage;
      const lit = [[1, 1, 1, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, 1, 1, 1]][stg];
      const phi = TAU * VIS * fig.t, cs = Math.cos(phi);
      for (const o of fig.L.ops) p.text(o.s, o.x, o.y, {size: 20, weight: 600, color: c.muted, align: "center", halo: false});
      const names = ["free field", "structure", "excavated soil", "SSI system"];
      fig.L.p.forEach((b, i) => {
        const g = p.g;
        g.save();
        g.globalAlpha = lit[i] ? 1 : 0.32;
        frame(p, b, lit[i] && stg > 0 && stg < 4);
        p.text(names[i], b.x + 8, b.y + 12, {size: 11.5, weight: 600, color: c.ink});
        const gy = b.y + b.h * 0.42, cx = b.x + b.w / 2, soilH = b.y + b.h - 1 - gy;
        const bw = Math.min(0.5 * b.w, 86), D = emb ? Math.min(0.42 * soilH, 44) : 0;
        const bands = [{f: 0.3, color: c.soil1}, {f: 0.45, color: c.soil2}, {f: 0.25, color: c.rock}];
        const wave = (zf) => 3.5 * Math.cos(TAU * 0.45 * zf) * cs;
        const nodes = [];
        if (emb) { for (let a = 0; a <= 4; a++) for (let q = 0; q <= 3; q++) nodes.push([cx - bw / 2 + a * bw / 4, gy + q * D / 3]); }
        else for (let a = 0; a <= 4; a++) nodes.push([cx - bw / 2 + a * bw / 4, gy]);
        const building = (dx) => {
          if (emb) p.rect(cx - bw / 2 + dx, gy, bw, D, c.structFill, c.struct, 1.4);
          else p.rect(cx - bw / 2 + dx, gy - 6, bw, 6, c.structFill, c.struct, 1.2);
          const top = gy - (emb ? 0 : 6), hh = Math.min(gy - b.y - 30, 70);
          p.rect(cx - bw * 0.3 + dx, top - hh, bw * 0.6, hh, null, c.struct, 1.6);
          for (let fl = 1; fl < 3; fl++) p.line([[cx - bw * 0.3 + dx, top - hh * fl / 3], [cx + bw * 0.3 + dx, top - hh * fl / 3]], c.struct, 1.2);
        };
        if (i === 0) {
          soilBlock(p, b.x + 1, gy, b.w - 2, soilH, bands, wave, Math.round(b.w / 16));
          if (emb) p.rect(cx - bw / 2, gy, bw, D, null, c.ink, 1);
          for (const [x, y] of nodes) {
            const u = 3.5 * Math.cos(TAU * 0.45 * (y - gy) / soilH) * cs;
            p.circle(x + u, y, 2.8, c.s1);
            if (stg === 1) p.arrow(x + u, y, x + u + 9 * Math.sign(cs || 1), y, c.s1, 1.2, 4);
          }
          p.text("X_{ff}, U′_{f} at the interaction nodes", b.x + 8, b.y + b.h - 10, {size: 9.5, color: c.muted});
        } else if (i === 1) {
          p.line([[b.x + 6, gy], [b.x + b.w - 6, gy]], withAlpha(c.soilLine, 0.6), 1, [4, 3]);
          building(0);
          for (const [x, y] of nodes) p.circle(x, y, 2.4, c.struct);
          p.text("K*_{s}, M_{s} (HOUSE)", b.x + 8, b.y + b.h - 10, {size: 9.5, color: c.muted});
        } else if (i === 2) {
          p.line([[b.x + 6, gy], [b.x + b.w - 6, gy]], withAlpha(c.soilLine, 0.6), 1, [4, 3]);
          if (emb) {
            p.rect(cx - bw / 2, gy, bw, D * 0.45, c.soil1); p.rect(cx - bw / 2, gy + D * 0.45, bw, D * 0.55, c.soil1);
            for (let a = 1; a < 4; a++) p.line([[cx - bw / 2 + a * bw / 4, gy], [cx - bw / 2 + a * bw / 4, gy + D]], c.soilLine, 1);
            for (let q = 1; q < 3; q++) p.line([[cx - bw / 2, gy + q * D / 3], [cx + bw / 2, gy + q * D / 3]], c.soilLine, 1);
            p.rect(cx - bw / 2, gy, bw, D, null, c.soilLine, 1.4);
            for (const [x, y] of nodes) p.circle(x, y, 2.4, c.soilLine);
            p.text("K*_{e}, M_{e} (HOUSE)", b.x + 8, b.y + b.h - 10, {size: 9.5, color: c.muted});
          } else {
            p.text("none for a surface mat", cx, gy + 22, {size: 10.5, color: c.muted, align: "center"});
          }
        } else {
          soilBlock(p, b.x + 1, gy, b.w - 2, soilH, bands, wave, Math.round(b.w / 16));
          const dx = 2.5 * cs;
          if (emb) p.rect(cx - bw / 2 + dx, gy, bw, D, c.bg);
          building(dx);
          for (const [x, y] of nodes) p.circle(x + dx, y, 2.6, c.s1);
          p.text("solved for U at each ω (ANALYS)", b.x + 8, b.y + b.h - 10, {size: 9.5, color: c.muted});
        }
        g.restore();
      });
    },
  };
  function paintSub(fig) {
    const on = SUB_ON[fig.st.stage] || [];
    const emb = fig.st.case === "embedded";
    (fig.data.terms || []).forEach((t, i) => { t.classList.toggle("on", on.includes(i)); t.classList.toggle("zero", !emb && i === 1); });
    fig.read([`${fig.st.stage === 0 ? "" : `Step ${fig.st.stage} of 4. `}${SUB_TEXT[fig.st.case][fig.st.stage]}`]);
  }

  // ================================================================== figure: layered soil column
  FIG["soil-column"] = {
    title: "A vertically propagating shear wave in the layered site",
    alt: "Animated displacement profile of the lesson 2 soil column over a half-space at the chosen frequency, with the amplification of the surface over the outcrop and the within motion at the top of the rock.",
    animated: true,
    params: {f: {def: 7.3, min: 0.2, max: 25, step: 0.01}, vhs: {def: 1000, min: 400, max: 100000, step: 10}, damp: {def: 1, min: 0.25, max: 3, step: 0.01}},
    tex: ["G^* = G\\,c(\\beta),\\quad V^* = V\\left(\\sqrt{1-\\beta^2} + i\\beta\\right),\\quad k^* = \\omega/V^*",
      "u_m(z) = E_m\\,e^{ik^*_m z} + F_m\\,e^{-ik^*_m z}",
      "\\text{within} = E + F,\\qquad \\text{outcrop} = 2E"],
    note: "Exact one-dimensional solution in each layer with the transfer of the up- and down-going amplitudes E, F across the interfaces (the SHAKE recursion of Theory Section 12.1), normalised to unit motion at the surface (the control point of lesson 2). Site of lesson 2: 5 m sand (Vs 300 m/s, 5 %), 12 m gravel (500 m/s, 4 %), rock half-space (1000 m/s, 2 %). Profile scaled to its largest amplitude; motion slowed down.",
    init(fig) {
      fig.slider("f", "Frequency f", {fmt: (v) => fmt(v, 2) + " Hz"});
      fig.slider("vhs", "Rock Vs", {log: true, fmt: (v) => `${Math.round(v)} m/s` + (v >= 50000 ? " (≈ rigid)" : "")});
      fig.slider("damp", "Soil damping ×", {fmt: (v) => `× ${fmt(v, 2)} (${fmt(5 * v, 1)} %, ${fmt(4 * v, 1)} %)`});
      fig.drag(() => fig.data.ax && fig.data.ax.box, (x) => fig.set("f", clamp(+fig.data.ax.ix(x).toFixed(2), 0.2, 25)));
    },
    change(fig, k) {
      const d = fig.data, st = fig.st;
      const col = {g: PH.COLUMN.g, layers: PH.COLUMN.layers.map((l) => Object.assign({}, l, {beta: Math.min(0.45, l.beta * st.damp)})),
        hs: Object.assign({}, PH.COLUMN.hs, {vs: st.vhs})};
      if (!d.curve || k === "vhs" || k === "damp" || k === "_init") {
        const fs = [], oc = [], wi = [];
        for (let f = 0.05; f <= 25.0001; f += 0.025) { const rr = PH.columnRatios(col, f); fs.push(f); oc.push(C.abs(rr.outcrop)); wi.push(C.abs(rr.within)); }
        const first = (a) => { let i = 0; for (let q = 1; q < fs.length && fs[q] < 11; q++) if (a[q] > a[i]) i = q; return i; };
        d.curve = {fs, oc, wi, io: first(oc), iw: first(wi)};
        d.col = col;
      }
      const rr = PH.columnRatios(col, st.f);
      d.wv = rr.wv; d.rr = rr;
      const zs = [];
      for (let z = 0; z <= 22.0001; z += 0.25) zs.push(z);
      d.prof = zs.map((z) => [z, PH.columnU(rr.wv, col, z)]);
      d.pmax = Math.max(...d.prof.map((q) => C.abs(q[1])));
      const cv = d.curve, n = rr.wv.E.length - 1;
      const vbar = (5 * 300 + 12 * 500) / 17;
      fig.read([`f = ${fmt(st.f, 2)} Hz`, `surface / outcrop ${fmt(C.abs(rr.outcrop), 2)}`, `surface / within (17 m) ${fmt(C.abs(rr.within), 2)}`,
        `first peaks: outcrop ${fmt(cv.oc[cv.io], 2)} at ${fmt(cv.fs[cv.io], 1)} Hz, within ${fmt(cv.wi[cv.iw], 1)} at ${fmt(cv.fs[cv.iw], 1)} Hz`,
        `V̄s/4H = ${fmt(vbar, 0)}/(4 × 17) = ${fmt(vbar / 68, 1)} Hz`, `at 17 m per unit surface motion: within |E+F|/2 = ${fmt(C.abs(C.add(rr.wv.E[n], rr.wv.F[n])) / 2, 2)}, outcrop |E| = ${fmt(C.abs(rr.wv.E[n]), 2)}`]);
    },
    layout(fig, W) {
      const pad = 8;
      if (W >= 620) {
        const cw = Math.min(300, Math.round(0.4 * W)), h = 300;
        fig.L = {col: {x: pad, y: pad, w: cw, h: h - pad}, chart: {x: cw + 2 * pad + 40, y: pad + 22, w: W - cw - 3 * pad - 40 - 8, h: h - 2 * pad - 22 - 28}};
        return h;
      }
      const h1 = 260, h2 = 200;
      fig.L = {col: {x: pad, y: pad, w: W - 2 * pad, h: h1 - pad}, chart: {x: 44, y: h1 + 22, w: W - 44 - 14, h: h2 - 22 - 30}};
      return h1 + h2;
    },
    draw(fig, p) {
      const d = fig.data, c = p.c, b = fig.L.col;
      if (!d.prof) return;
      const phi = TAU * VIS * fig.t, e = [Math.cos(phi), Math.sin(phi)];
      const re = (z) => z[0] * e[0] - z[1] * e[1];
      frame(p, b, false);
      const top = b.y + 26, H = b.h - 34, zmax = 22, Z = (z) => top + z / zmax * H;
      const cw = Math.min(110, 0.36 * b.w), cx = b.x + 12 + cw / 2 + 18;
      p.rect(cx - cw / 2, Z(0), cw, Z(5) - Z(0), c.soil1);
      p.rect(cx - cw / 2, Z(5), cw, Z(17) - Z(5), c.soil2);
      p.rect(cx - cw / 2, Z(17), cw, Z(zmax) - Z(17), c.rock);
      for (const z of [0, 5, 10, 17, 22]) p.text(`${z} m`, cx - cw / 2 - 4, Z(z), {size: 9.5, color: c.muted, align: "right", halo: false});
      const lx = cx + cw / 2 + 10, lw = b.x + b.w - lx - 4;
      const st = fig.st;
      const lab = (y, a, bb) => { p.text(a, lx, y - 6, {size: 10.5, color: c.ink, halo: false}); p.text(bb, lx, y + 7, {size: 9.5, color: c.muted, halo: false}); };
      if (lw > 70) {
        lab(Z(2.5), "sand", `Vs 300 m/s, β ${fmt(5 * st.damp, 1)} %`);
        lab(Z(11), "gravel", `Vs 500 m/s, β ${fmt(4 * st.damp, 1)} %`);
        lab(Z(19.5), "rock half-space", `Vs ${Math.round(st.vhs)} m/s, 2 %`);
      }
      // displacement profile and its envelope
      const A = 0.42 * cw / d.pmax;
      const env = d.prof.map(([z, u]) => [cx + A * C.abs(u), Z(z)]), env2 = d.prof.map(([z, u]) => [cx - A * C.abs(u), Z(z)]);
      p.line(env, withAlpha(c.s1, 0.35), 1, [3, 3]); p.line(env2, withAlpha(c.s1, 0.35), 1, [3, 3]);
      p.line([[cx, Z(0)], [cx, Z(zmax)]], withAlpha(c.ink, 0.35), 1);
      const prof = d.prof.map(([z, u]) => [cx + A * re(u), Z(z)]);
      p.line(prof, c.s1, 2.2);
      for (const [z, u] of d.prof.filter((q, i) => i % 8 === 0)) p.circle(cx + A * re(u), Z(z), 2.4, c.s1);
      p.line([[cx - cw / 2, Z(17)], [cx + cw / 2, Z(17)]], c.ink, 1.2);
      p.text("control point (surface)", cx, Z(0) - 9, {size: 9.5, color: c.muted, align: "center"});
      // chart
      const cv = d.curve;
      const ax = axes(p, fig.L.chart, {x0: 0, x1: 25, y0: 0.2, y1: 100, ylog: true, yt: [0.2, 0.5, 1, 2, 5, 10, 20, 50, 100], xlabel: "f (Hz)", ylabel: "amplification of the surface motion"});
      d.ax = ax;
      const vbar = (5 * 300 + 12 * 500) / 17, fq = vbar / 68;
      p.line([[ax.X(fq), ax.box.y], [ax.X(fq), ax.box.y + ax.box.h]], c.ref, 1, [2, 3]);
      p.text("V̄s/4H", ax.X(fq) - 3, ax.box.y + ax.box.h - 8, {size: 9.5, color: c.muted, align: "right"});
      plot(p, ax, cv.fs, cv.wi, c.s2, 2);
      plot(p, ax, cv.fs, cv.oc, c.s1, 2);
      p.text("surface / within at 17 m (E+F)", ax.X(cv.fs[cv.iw]) + 6, ax.Y(Math.min(90, cv.wi[cv.iw])) + 2, {size: 10.5, color: c.ink});
      p.text("surface / outcrop (2E)", ax.X(cv.fs[cv.io]) + 6, ax.Y(cv.oc[cv.io]) - 9, {size: 10.5, color: c.ink});
      const X = ax.X(st.f);
      p.line([[X, ax.box.y], [X, ax.box.y + ax.box.h]], c.ink, 1, [3, 3]);
      p.circle(X, ax.Y(clamp(C.abs(d.rr.within), 0.2, 100)), 4.5, c.s2, c.bg, 2);
      p.circle(X, ax.Y(clamp(C.abs(d.rr.outcrop), 0.2, 100)), 4.5, c.s1, c.bg, 2);
    },
  };

  // ================================================================== figure: impedance ellipse
  const IMP = {k: 2.82e6, U: 1e-3};       // lesson 3: Kx of the 12 m mat (2 m mesh, 0.2 Hz) in kN/m; U = 1 mm
  FIG["impedance-ellipse"] = {
    title: "A spring and a dashpot under harmonic motion: K(ω) = k + iωc",
    alt: "A footing on a spring and a dashpot moved harmonically, the force and displacement histories with their phase lag, and the force-displacement ellipse whose area is the energy dissipated per cycle.",
    animated: true,
    params: {f: {def: 5, min: 0.2, max: 12, step: 0.01}, c: {def: 5.76e4, min: 0, max: 2e5, step: 100}},
    tex: ["K(\\omega) = k + i\\omega c,\\qquad F = k\\,u + c\\,\\dot u,\\quad u = U\\sin\\omega t",
      "\\varphi = \\arctan\\frac{\\omega c}{k},\\qquad E_D = \\oint F\\,du = \\pi c\\,\\omega\\,U^2,\\qquad E_S = \\tfrac12 k U^2",
      "\\xi = \\frac{E_D}{4\\pi E_S} = \\frac{\\omega c}{2k}"],
    note: "Defaults from lesson 3: k = 2.82e6 kN/m (Kx of the 12 m × 12 m mat at 0.2 Hz) and c = ρVsA = 5.76e4 kN s/m (the plane-wave dashpot the FOUNDASH values approach); U = 1 mm (illustrative). Motion slowed down.",
    init(fig) {
      fig.slider("f", "Frequency f", {fmt: (v) => fmt(v, 2) + " Hz"});
      fig.slider("c", "Dashpot c", {fmt: (v) => fmtSig(v, 3) + " kN s/m"});
    },
    change(fig) {
      const w = TAU * fig.st.f, c = fig.st.c, k = IMP.k, U = IMP.U;
      const d = fig.data;
      d.w = w; d.phi = Math.atan2(w * c, k); d.ED = Math.PI * c * w * U * U; d.ES = 0.5 * k * U * U; d.xi = w * c / (2 * k);
      d.Fmax = U * Math.hypot(k, w * c);
      fig.read([`ω = ${fmt(w, 1)} rad/s`, `φ = atan(ωc/k) = ${fmt(d.phi * 180 / Math.PI, 1)}°`, `E_{D} = πcωU² = ${fmt(d.ED, 2)} kJ`,
        `E_{S} = ½kU² = ${fmt(d.ES, 2)} kJ`, `ξ = E_{D}/(4πE_{S}) = ωc/(2k) = ${fmt(d.xi, 3)}`, `|K| = ${fmtSig(Math.hypot(k, w * c), 3)} kN/m`]);
    },
    layout(fig, W) {
      const pad = 8;
      if (W >= 640) {
        const w1 = Math.round(0.24 * W), w2 = Math.round(0.36 * W), h = 280;
        fig.L = {s: {x: pad, y: pad, w: w1 - pad, h: h - 2 * pad}, t: {x: w1 + 40, y: pad + 22, w: w2 - 52, h: h - 2 * pad - 22 - 28},
          e: {x: w1 + w2 + 44, y: pad + 22, w: W - w1 - w2 - 44 - 12, h: h - 2 * pad - 22 - 28}};
        return h;
      }
      const h1 = 150, h2 = 160, h3 = Math.min(280, W * 0.75);
      fig.L = {s: {x: pad, y: pad, w: W - 2 * pad, h: h1 - pad}, t: {x: 40, y: h1 + 22, w: W - 40 - 12, h: h2 - 22 - 28},
        e: {x: 52, y: h1 + h2 + 22, w: W - 52 - 12, h: h3 - 22 - 30}};
      return h1 + h2 + h3;
    },
    draw(fig, p) {
      const d = fig.data, c = p.c, L = fig.L;
      if (!d.w) return;
      const psi = TAU * VIS * fig.t;
      const u = Math.sin(psi), F = Math.sin(psi + d.phi);           // u/U and F/Fmax
      // ---- schematic
      const s = L.s;
      frame(p, s, false);
      const wx = s.x + 18, cy = s.y + s.h / 2 - 2, bx = s.x + s.w * 0.62 + 16 * u;
      p.line([[wx, cy - 34], [wx, cy + 34]], c.struct, 2);
      for (let y = cy - 30; y < cy + 34; y += 7) p.line([[wx, y], [wx - 6, y + 5]], c.struct, 1);
      p.spring(wx, cy - 14, bx - 20, cy - 14, c.struct, 4, 5);
      p.dashpot(wx, cy + 14, bx - 20, cy + 14, c.struct, 6);
      p.rect(bx - 20, cy - 28, 40, 56, c.structFill, c.struct, 1.3);
      p.text("k", (wx + bx - 20) / 2, cy - 27, {size: 11, color: c.muted, align: "center"});
      p.text("c", (wx + bx - 20) / 2, cy + 30, {size: 11, color: c.muted, align: "center"});
      p.arrow(bx + 20, cy, bx + 20 + 34 * F, cy, c.s2, 2, 7);
      p.text("F(t)", bx + 24, cy - 12, {size: 10.5, color: c.ink});
      p.line([[s.x + s.w * 0.62, cy + 40], [s.x + s.w * 0.62, cy + 46]], c.muted, 1);
      p.arrow(s.x + s.w * 0.62, cy + 43, bx, cy + 43, c.s1, 1.6, 6);
      p.text("u(t)", s.x + s.w * 0.62, cy + 54, {size: 10.5, color: c.ink, align: "center"});
      p.text("footing on the soil", s.x + 8, s.y + 12, {size: 11, weight: 600, color: c.ink});
      // ---- histories over two periods
      const at = axes(p, L.t, {x0: 0, x1: 2, y0: -1.4, y1: 1.4, yt: [-1, 0, 1], xt: [0, 0.5, 1, 1.5, 2], xlabel: "t / T", ylabel: "u / U,  F / F_{max}"});
      const xs = [], yu = [], yf = [];
      for (let i = 0; i <= 200; i++) { const x = 2 * i / 200, ps = TAU * x; xs.push(x); yu.push(Math.sin(ps)); yf.push(Math.sin(ps + d.phi)); }
      plot(p, at, xs, yu, c.s1, 2);
      plot(p, at, xs, yf, c.s2, 2);
      const tn = (psi / TAU) % 2;
      p.line([[at.X(tn), at.box.y], [at.X(tn), at.box.y + at.box.h]], c.ink, 1, [3, 3]);
      p.circle(at.X(tn), at.Y(Math.sin(TAU * tn)), 4, c.s1, c.bg, 1.5);
      p.circle(at.X(tn), at.Y(Math.sin(TAU * tn + d.phi)), 4, c.s2, c.bg, 1.5);
      const xu = 0.25, xf = 0.25 - d.phi / TAU;
      p.line([[at.X(xf), at.Y(1.16)], [at.X(xu), at.Y(1.16)]], c.ink, 1);
      p.line([[at.X(xf), at.Y(1.1)], [at.X(xf), at.Y(1.22)]], c.ink, 1); p.line([[at.X(xu), at.Y(1.1)], [at.X(xu), at.Y(1.22)]], c.ink, 1);
      p.text("φ", at.X(xu) + 8, at.Y(1.2), {size: 11, color: c.ink});
      p.text("u", at.X(1.25) + 4, at.Y(1) + 8, {size: 10.5, color: c.ink});
      p.text("F", at.X(1.25 - d.phi / TAU) - 4, at.Y(1) - 8, {size: 10.5, color: c.ink, align: "right"});
      // ---- force-displacement ellipse
      const Fm = d.Fmax, U = IMP.U, kN = (v) => v;              // forces in kN, displacement in mm
      const ae = axes(p, L.e, {x0: -1.25, x1: 1.25, y0: -1.25 * Fm, y1: 1.25 * Fm, xlabel: "u (mm)", ylabel: "F (kN)"});
      const ell = [];
      for (let i = 0; i <= 160; i++) { const ps = TAU * i / 160; ell.push([ae.X(Math.sin(ps)), ae.Y(kN(Fm * Math.sin(ps + d.phi)))]); }
      p.poly(ell, withAlpha(c.s2, 0.14), null);
      p.poly([[ae.X(0), ae.Y(0)], [ae.X(1), ae.Y(0)], [ae.X(1), ae.Y(IMP.k * U)]], withAlpha(c.s3, 0.18), null);
      p.line([[ae.X(-1.2), ae.Y(-1.2 * IMP.k * U)], [ae.X(1.2), ae.Y(1.2 * IMP.k * U)]], c.s3, 1.5, [5, 3]);
      p.line(ell, c.s2, 2);
      p.circle(ae.X(u), ae.Y(Fm * F), 5, c.s2, c.bg, 2);
      p.text("spring k·u", ae.X(-0.62) + 4, ae.Y(-0.62 * IMP.k * U) + 12, {size: 10, color: c.ink});
      p.text("E_{S}", ae.X(0.72), ae.Y(0.3 * IMP.k * U), {size: 10, color: c.ink, align: "center"});
      p.text("area = E_{D}", ae.X(-0.05), ae.Y(-0.08 * Fm), {size: 10.5, color: c.ink, align: "right"});
    },
  };

  // ================================================================== shared: synthetic motion (cached)
  let MOTION = null;
  const motion = () => MOTION || (MOTION = PH.synthMotion());

  // ================================================================== figure: from TF to ISRS
  const BANK = [1, 1.6, 2.5, 4, 6.3, 10, 16, 25];
  FIG["tf-to-isrs"] = {
    title: "From transfer function to in-structure response spectrum",
    alt: "Animated pipeline: a control motion, the floor transfer function, the floor motion, a bank of oscillators responding to it and the response spectrum building up point by point.",
    animated: true,
    params: {f0: {def: 3.5, min: 1, max: 12, step: 0.05}, beta: {def: 0.05, min: 0.01, max: 0.2, step: 0.005}},
    tex: ["a(t) = \\operatorname{IFFT}\\left[H(f)\\,A(f)\\right],\\qquad H(f) = 1 + \\Gamma_1\\phi_1\\,\\frac{r^2}{c(\\beta) - r^2},\\quad r = f/f_0",
      "\\ddot x + 2\\zeta\\omega\\dot x + \\omega^2 x = -a(t),\\qquad S_a(f) = \\max_t\\,\\lvert 2\\zeta\\omega\\dot x + \\omega^2 x\\rvert"],
    note: "Synthetic control motion (random phases, Kanai-Tajimi spectrum, 10 s, PGA 0.3 g): illustrative. The floor transfer function is one hysteretic mode at f₀ with damping β and Γ₁φ₁ = 1.34 (the roof of lesson 1). Oscillators of 5 % damping integrated exactly between samples (Nigam-Jennings, at least 64 steps per period), as MOTION does. Played at 0.75 × speed; oscillator deflections scaled each to its own maximum.",
    init(fig) {
      fig.slider("f0", "Floor mode f₀", {fmt: (v) => fmt(v, 2) + " Hz"});
      fig.slider("beta", "Mode damping β", {fmt: (v) => fmt(100 * v, 1) + " %"});
      fig.data.speed = 0.75;
    },
    change(fig, k) {
      const d = fig.data;
      if (k === "_layout") return;
      clearTimeout(d.timer);
      const run = () => {
        const mo = motion(), nT = 1200, dt = mo.dt;
        const modes = [{f: fig.st.f0, beta: fig.st.beta, P: 1.34}];
        const fl = PH.convolve(mo, (f) => PH.modalTF(f, modes));
        const FS = PH.logspace(0.5, 50, 100);
        d.nT = nT; d.dt = dt; d.ag = mo.a; d.af = fl; d.FS = FS;
        d.runs = FS.map((f) => PH.oscillator(fl, dt, f, 0.05, nT).run);
        d.bank = BANK.map((f) => { const o = PH.oscillator(fl, dt, f, 0.05, nT); let mx = 0; for (const v of o.x) mx = Math.max(mx, Math.abs(v)); return {f, x: o.x, run: o.run, mx}; });
        if (!d.ground) d.ground = FS.map((f) => PH.oscillator(mo.a, dt, f, 0.05, nT).sa);
        d.final = d.runs.map((r) => r[nT - 1]);
        d.Hf = []; for (let f = 0; f <= 20.0001; f += 0.05) d.Hf.push([f, C.abs(PH.modalTF(f, modes))]);
        d.agmax = Math.max(...Array.from(mo.a.slice(0, nT), Math.abs)); d.afmax = Math.max(...Array.from(fl.slice(0, nT), Math.abs));
        d.ready = true;
        fig.redraw();
      };
      if (k === "_init" || !d.ready) run(); else d.timer = setTimeout(run, 40);
    },
    layout(fig, W) {
      const pad = 8;
      if (W >= 700) {
        const r1 = 128, r2 = 220, gap = 26, w3 = (W - 2 * pad - 2 * gap) / 3;
        fig.L = {wide: true, ag: {x: pad + 34, y: pad + 20, w: w3 - 40, h: r1 - 20 - 26},
          H: {x: pad + w3 + gap + 34, y: pad + 20, w: w3 - 40, h: r1 - 20 - 26},
          af: {x: pad + 2 * (w3 + gap) + 34, y: pad + 20, w: w3 - 40, h: r1 - 20 - 26},
          bank: {x: pad, y: r1 + 2 * pad, w: Math.round(0.38 * W), h: r2 - 2 * pad},
          rs: {x: Math.round(0.38 * W) + 2 * pad + 40, y: r1 + 2 * pad + 20, w: W - Math.round(0.38 * W) - 3 * pad - 40 - 6, h: r2 - 2 * pad - 20 - 28},
          arrows: [[pad + w3 + 4, pad + r1 / 2 - 8], [pad + 2 * w3 + gap + 4, pad + r1 / 2 - 8]]};
        return r1 + r2 + pad;
      }
      const h = [98, 110, 98, 150, 210];
      let y = pad;
      const box = (hh, lx) => { const b = {x: 40, y: y + 20, w: W - 40 - 12, h: hh - 20 - 26}; y += hh; return b; };
      const ag = box(h[0]), H = box(h[1]), af = box(h[2]);
      const bank = {x: pad, y: y, w: W - 2 * pad, h: h[3] - 8}; y += h[3];
      const rs = box(h[4]);
      fig.L = {wide: false, ag, H, af, bank, rs, arrows: []};
      return y + pad;
    },
    draw(fig, p) {
      const d = fig.data, c = p.c, L = fig.L;
      if (!d.ready) { p.text("computing …", 12, 20, {color: c.muted}); return; }
      const total = d.nT * d.dt, hold = 3, cyc = total / d.speed + hold;
      const tt = Math.min(total, ((fig.t % cyc)) * d.speed);
      const idx = !fig.playing && fig.t === 0 ? d.nT - 1 : Math.max(1, Math.min(d.nT - 1, Math.floor(tt / d.dt)));
      const tnow = idx * d.dt;
      const ts = [], A = [], Fl = [];
      for (let i = 0; i < d.nT; i += 2) { ts.push(i * d.dt); A.push(d.ag[i]); Fl.push(d.af[i]); }
      const ym = Math.max(d.agmax, d.afmax) * 1.1;
      const a1 = axes(p, L.ag, {x0: 0, x1: total, y0: -ym, y1: ym, yt: linTicks(-ym, ym, 3), title: "control motion a_{g}(t) (g)", xlabel: "t (s)"});
      plot(p, a1, ts, A, c.s4, 1.2);
      p.line([[a1.X(tnow), a1.box.y], [a1.X(tnow), a1.box.y + a1.box.h]], c.ink, 1);
      const Hmax = Math.max(...d.Hf.map((q) => q[1]));
      const a2 = axes(p, L.H, {x0: 0, x1: 20, y0: 0, y1: Math.ceil(Hmax * 1.1), title: "floor TF |H(f)|", xlabel: "f (Hz)"});
      plot(p, a2, d.Hf.map((q) => q[0]), d.Hf.map((q) => q[1]), c.s3, 2);
      const a3 = axes(p, L.af, {x0: 0, x1: total, y0: -ym, y1: ym, yt: linTicks(-ym, ym, 3), title: "floor motion a(t) (g)", xlabel: "t (s)"});
      plot(p, a3, ts, Fl, c.s1, 1.2);
      p.line([[a3.X(tnow), a3.box.y], [a3.X(tnow), a3.box.y + a3.box.h]], c.ink, 1);
      for (const [x, y] of L.arrows) p.text("→", x + 8, y + 10, {size: 18, color: c.muted, align: "center", halo: false});
      // oscillator bank
      const b = L.bank;
      frame(p, b, false);
      p.text("oscillators on the floor (5 %)", b.x + 8, b.y + 12, {size: 11, weight: 600, color: c.ink});
      const fy = b.y + b.h - 26, n = BANK.length, sp = (b.w - 20) / n, hh = Math.min(b.h - 64, 90);
      p.rect(b.x + 8, fy, b.w - 16, 6, c.structFill, c.struct, 1);
      d.bank.forEach((o, i) => {
        const x0 = b.x + 10 + sp * (i + 0.5), hs = hh * (0.55 + 0.45 * (1 - i / (n - 1)));
        const defl = 0.32 * sp * (o.x[idx] / (o.mx || 1));
        stick(p, x0, fy, hs, (zf) => defl * shapeCant(zf), {color: c.struct, w: 2, r: 5, fill: c.s1});
        p.text(o.f >= 10 ? String(o.f) : String(o.f), x0, fy + 14, {size: 9.5, color: c.muted, align: "center", halo: false});
      });
      p.text("Hz", b.x + b.w - 10, fy + 14, {size: 9.5, color: c.muted, align: "right", halo: false});
      // spectrum building up
      const cur = d.runs.map((r) => r[idx]);
      const smax = Math.max(...d.final, ...d.ground) * 1.12;
      const a4 = axes(p, L.rs, {x0: 0.5, x1: 50, xlog: true, y0: 0, y1: smax, xlabel: "f (Hz)", title: "ISRS, 5 % (g)"});
      plot(p, a4, d.FS, d.ground, c.ref, 1.5, [4, 3]);
      plot(p, a4, d.FS, cur, c.s1, 2);
      d.bank.forEach((o) => { p.circle(a4.X(o.f), a4.Y(o.run[idx]), 4, c.s1, c.bg, 1.5); });
      p.text("control motion", a4.X(30), a4.Y(d.ground[d.FS.findIndex((f) => f >= 30)]) - 9, {size: 10, color: c.muted, align: "center"});
      let ip = 0; cur.forEach((v, i) => { if (v > cur[ip]) ip = i; });
      if (d.lastRead !== undefined && Math.abs(idx - d.lastRead) < 40 && idx !== d.nT - 1) return;
      d.lastRead = idx;
      fig.read([`t = ${fmt(tnow, 1)} s`, `floor peak acceleration so far ${fmt(Math.max(...Array.from(d.af.slice(0, idx + 1), Math.abs)), 2)} g (control motion ${fmt(Math.max(...Array.from(d.ag.slice(0, idx + 1), Math.abs)), 2)} g)`,
        `ISRS peak so far ${fmt(cur[ip], 2)} g at ${fmt(d.FS[ip], 2)} Hz`]);
    },
  };

  // ================================================================== figure: interaction-node sets
  const SET_INFO = {
    fv: ["FV", "every node of the excavated soil is an interaction node: the exact substructuring equation, the reference."],
    fsin: ["FI-FSIN", "only the lateral and bottom faces (subtraction method). The interior and the roof-slab nodes hang on the excavated soil alone: spurious peaks at the frequencies of that sub-system (near 6 Hz in lesson 5)."],
    evbn: ["FI-EVBN", "FSIN plus the ground-surface face (modified subtraction). The surface nodes raise the frequencies of the unphysical sub-system."],
    ffv: ["FFV", "EVBN plus the interior of every second internal level (z = −3 m and −1 m): small, stiff sub-volumes, close to FV."],
  };
  FIG["interaction-sets"] = {
    title: "Interaction-node sets on the excavation mesh of lesson 5",
    alt: "Rotating isometric view of the 5 by 5 by 6 node excavation mesh with the interaction nodes of the selected method highlighted and the node counts.",
    animated: true,
    params: {method: {def: "fv", values: ["fv", "fsin", "evbn", "ffv"]}},
    tex: ["\\text{memory of } X_{ff} = (3N_{\\text{int}})^2 \\times 16\\ \\text{bytes},\\qquad \\text{time} \\propto N_{\\text{int}}^3"],
    note: "10 m × 10 m × 5 m excavation, nodes every 2.5 m horizontally and 1 m vertically (150 nodes, numbered bottom-up). The sets are those of INTGEN 1, 3, 2 and 5 (skip 2); the counts are what INTCOUNT reports in the lesson. Drag the view to turn it.",
    init(fig) {
      fig.choice("method", "Interaction set", ["fv", "fsin", "evbn", "ffv"].map((m) => [m, SET_INFO[m][0]]));
      fig.data.nodes = PH.boxNodes();
      fig.data.counts = {};
      for (const m of ["fv", "fsin", "evbn", "ffv"]) fig.data.counts[m] = PH.setCount(m);
      fig.data.rot = 0;
      fig.drag(() => ({x: 0, y: 0, w: fig.W, h: fig.H}), (x, y, ev) => {
        if (ev.type === "pointerdown") { fig.data.dx = x; return; }
        fig.data.rot += (x - fig.data.dx) * 0.01; fig.data.dx = x; fig.redraw();
      });
    },
    change(fig) {
      const m = fig.st.method, n = fig.data.counts[m];
      fig.read([["fv", "fsin", "evbn", "ffv"].map((q) => `${SET_INFO[q][0]} ${fig.data.counts[q]}`).join(" · "),
        `${SET_INFO[m][0]}: ${n} interaction nodes, ${3 * n} DOFs, X_{ff} ${fmt((3 * n) ** 2 * 16 / 1e6, 2)} MB`, `${SET_INFO[m][0]}: ${SET_INFO[m][1]}`]);
    },
    layout(fig, W) { return W >= 600 ? 340 : 280; },
    draw(fig, p) {
      const c = p.c, W = fig.W, H = fig.H, nodes = fig.data.nodes, o = nodes.o, m = fig.st.method;
      const a = 0.62 + 0.22 * Math.sin(0.35 * fig.t) + fig.data.rot, ca = Math.cos(a), sa = Math.sin(a);
      // projection fitted to the canvas (true vertical scale), above the legend
      const G = 7, raw = (x, y, z) => [x * ca - y * sa, -z + (x * sa + y * ca) * 0.42];
      const x1 = G * Math.SQRT2, x0 = -x1, y0 = -0.42 * G * Math.SQRT2, y1 = 5 + 0.42 * 5 * Math.SQRT2;   // any rotation
      const avW = W - 40, avH = H - 56, s = Math.min(avW / (x1 - x0), avH / (y1 - y0));
      const ox = (W - s * (x1 + x0)) / 2, oy = 10 + (avH - s * (y1 - y0)) / 2 - s * y0;
      const P = (x, y, z) => { const r = raw(x, y, z); return [ox + s * r[0], oy + s * r[1]]; };
      const depth = (x, y) => x * sa + y * ca;
      // ground surface around the excavation
      p.poly([P(-G, -G, 0), P(G, -G, 0), P(G, G, 0), P(-G, G, 0)], withAlpha(c.soil1, 0.8), c.soilLine, 1);
      const hx = 5;
      p.poly([P(-hx, -hx, 0), P(hx, -hx, 0), P(hx, hx, 0), P(-hx, hx, 0)], c.bg, c.ink, 1);
      // grid lines of the mesh (faint)
      for (let k = 0; k < o.nz; k++) {
        const z = -(o.nz - 1) * o.dz + k * o.dz;
        p.poly([P(-hx, -hx, z), P(hx, -hx, z), P(hx, hx, z), P(-hx, hx, z)], null, withAlpha(c.axis, 0.7), 0.8);
      }
      for (const [x, y] of [[-hx, -hx], [hx, -hx], [hx, hx], [-hx, hx]]) p.line([P(x, y, 0), P(x, y, -5)], c.axis, 1);
      const sorted = nodes.map((nd) => ({nd, d: depth(nd.x, nd.y)})).sort((u, v) => u.d - v.d || u.nd.z - v.nd.z);
      for (const {nd} of sorted) {
        const [x, y] = P(nd.x, nd.y, nd.z), on = PH.inSet(m, nd, o);
        if (on) p.circle(x, y, 4.3, c.s1, c.bg, 1.4); else p.circle(x, y, 3, c.bg, c.ref, 1.3);
      }
      const lab = (t, q, dy) => p.text(t, Math.min(q[0] + 8, W - 4 - p.measure(t, 10)), q[1] + dy, {size: 10, color: c.muted});
      lab("z = 0 (grade)", P(hx, -hx, 0), -6);
      lab("z = −5 m", P(hx, -hx, -5), 0);
      p.circle(14, H - 30, 4.3, c.s1, c.bg, 1.4); p.text(`interaction node (${fig.data.counts[m]})`, 24, H - 30, {size: 10.5, color: c.ink});
      p.circle(14, H - 13, 3, c.bg, c.ref, 1.3); p.text(`excavated-soil node only (${150 - fig.data.counts[m]})`, 24, H - 13, {size: 10.5, color: c.ink});
    },
  };

  // ================================================================== figure: SV, SH and P waves
  const WAVES = [
    {name: "SV wave · X input", sub: ["WAVE,2", "V_{s} = 300 m/s"], V: 300, dir: [1, 0], axis: "x"},
    {name: "SH wave · Y input", sub: ["WAVE,4", "V_{s} = 300 m/s"], V: 300, dir: [0.62, -0.36], axis: "y"},
    {name: "P wave · Z input", sub: ["WAVE,3", "V_{p} = 600 m/s"], V: 600, dir: [0, -1], axis: "z"},
  ];
  FIG["wave-types"] = {
    title: "Vertically propagating SV, SH and P waves: the three input directions",
    alt: "Three soil columns with particles: a pulse travels up, doubles at the free surface and reflects; particles move in x for SV, in y for SH and vertically for P.",
    animated: true,
    params: {fp: {def: 4, min: 2, max: 10, step: 0.1}},
    tex: ["u(z,t) = p\\!\\left(t + \\tfrac{z}{V}\\right) + p\\!\\left(t - \\tfrac{z}{V}\\right),\\qquad V = V_s\\ (\\text{SV, SH}),\\ V_p\\ (\\text{P})"],
    note: "A Ricker pulse p(t) in a uniform half-space with a free surface (the sand of lesson 7: Vs = 300 m/s, Vp = 600 m/s), z the depth; each incident pulse reaches the surface at the same instant, where incident and reflected waves add (twice the incident motion). Particle motion along x (SV), y (SH) or z (P), the wave travels vertically. 40 m of soil shown; time slowed down about 12 times.",
    init(fig) { fig.slider("fp", "Pulse frequency", {fmt: (v) => fmt(v, 1) + " Hz"}); },
    change(fig) {
      const f = fig.st.fp;
      fig.read([`wavelength at ${fmt(f, 1)} Hz: S ${fmt(300 / f, 0)} m, P ${fmt(600 / f, 0)} m`, "surface motion = incident + reflected = 2 × incident (free surface)"]);
    },
    layout(fig, W) { return W >= 600 ? 360 : 340; },
    draw(fig, p) {
      const c = p.c, W = fig.W, H = fig.H, fp = fig.st.fp, Zm = 40;
      const T0 = Zm / 300 + 1.4 / fp, tphys = ((fig.t * 0.08) % (2 * T0)) - T0;
      const colw = (W - 28) / 3, cw = Math.min(150, colw - 40), dep = Math.min(42, cw * 0.38);
      const th = Math.max(...WAVES.map((wv) => wrap(p, wv.name, colw - 8, 11, 600).length)) * 14;
      const top = 6 + th + 30 + 50, bot = H - 22, hh = bot - top;
      WAVES.forEach((wv, i) => {
        const cx0 = 14 + i * colw, x0 = cx0 + (colw - cw - dep * 0.8) / 2;
        const ox = dep * 0.8, oy = -dep * 0.45;
        title(p, wv.name, cx0 + 4, 6, colw - 8, {size: 11});
        wv.sub.forEach((t, q) => p.text(t, cx0 + 4, 6 + th + 8 + 13 * q, {size: 10, color: c.muted}));
        p.poly([[x0, top], [x0 + cw, top], [x0 + cw + ox, top + oy], [x0 + ox, top + oy]], withAlpha(c.soil1, 0.9), c.soilLine, 1);
        p.poly([[x0 + cw, top], [x0 + cw + ox, top + oy], [x0 + cw + ox, bot + oy], [x0 + cw, bot]], withAlpha(c.soil2, 0.7), c.soilLine, 1);
        p.rect(x0, top, cw, hh, withAlpha(c.soil1, 0.6), c.soilLine, 1);
        const nx = cw > 120 ? 5 : 4, nz = 15, amp = Math.min(7, cw / 14);
        for (let kz = 0; kz < nz; kz++) {
          const zf = (kz + 0.5) / nz, z = zf * Zm, u = PH.halfspacePulse(z, tphys, wv.V, fp);
          for (let kx = 0; kx < nx; kx++) {
            const x = x0 + (kx + 0.5) * cw / nx, y = top + zf * hh;
            p.circle(x + amp * u * wv.dir[0] * 1.6, y + amp * u * wv.dir[1] * 1.6, 2.6, c.s1);
          }
        }
        const ax = x0 + cw / 2 + ox / 2, ay = top + oy - 14;
        p.arrow(ax, ay, ax + 16 * wv.dir[0], ay + 16 * wv.dir[1], c.ink, 1.3, 5);
        p.arrow(ax, ay, ax - 16 * wv.dir[0], ay - 16 * wv.dir[1], c.ink, 1.3, 5);
        p.text(wv.axis, ax + 20 * wv.dir[0] + (wv.axis === "z" ? 6 : 4), ay + 20 * wv.dir[1] - 2, {size: 10.5, color: c.ink});
        p.arrow(x0 - 8, bot - 4, x0 - 8, bot - 34, c.muted, 1.2, 5);
      });
      p.text("wave travels up ↑", 10, H - 8, {size: 10, color: c.muted});
    },
  };

  // ================================================================== figure: envelope and peak broadening
  const CASES = [{f: 4.4, name: "lower-bound soil", col: "s3"}, {f: 5.01, name: "best estimate", col: "s4"}, {f: 5.6, name: "upper-bound soil", col: "s2"}];
  FIG["isrs-broadening"] = {
    title: "Design ISRS: envelope of the soil cases and ±b peak broadening",
    alt: "Animated construction of a design floor spectrum: the spectra of the soil cases, their envelope, and the envelope broadened by plus or minus b around its peaks.",
    animated: true,
    params: {b: {def: 0.15, min: 0, max: 0.25, step: 0.005}, cases: {def: "1", values: ["1", "3"]}},
    tex: ["B(f) = \\max\\left\\{E(f') : \\frac{f}{1+b} \\le f' \\le \\frac{f}{1-b}\\right\\}",
      "\\text{peak at } f_p \\;\\rightarrow\\; \\text{plateau } \\bigl[f_p(1-b),\\ f_p(1+b)\\bigr]"],
    note: "Illustrative spectra (5 %): the synthetic motion of the transfer-function figure through one floor mode (β = 5 %, Γ₁φ₁ = 1.34) at 5.01 Hz, the frequency of lesson 7's SRSS peak (the spectrum of this motion peaks a little above it), and with three cases also at 4.4 and 5.6 Hz for lower- and upper-bound soil (illustrative shifts). BROADEN takes the envelope of its sources first, then broadens it by ±Smooth2 % (the window maximum above).",
    init(fig) {
      fig.slider("b", "Broadening ±b", {fmt: (v) => `±${fmt(100 * v, 1)} %`});
      fig.choice("cases", "Soil cases", [["1", "one"], ["3", "three"]]);
    },
    change(fig, k) {
      const d = fig.data;
      if (!d.spec) {
        const mo = motion();
        d.FS = PH.logspace(0.5, 50, 200);
        d.spec = CASES.map((cs) => { const fl = PH.convolve(mo, (f) => PH.modalTF(f, [{f: cs.f, beta: 0.05, P: 1.34}])); return PH.spectrum(fl, mo.dt, d.FS, 0.05, 1200); });
      }
      const use = fig.st.cases === "3" ? [0, 1, 2] : [1];
      d.use = use;
      d.env = use.length > 1 ? PH.envelope(use.map((i) => ({x: d.FS, y: d.spec[i]}))) : {x: d.FS.slice(), y: d.spec[1].slice()};
      d.br = PH.broaden(d.env.x, d.env.y, fig.st.b);
      let ip = 0; d.env.y.forEach((v, i) => { if (v > d.env.y[ip]) ip = i; });
      d.fp = d.env.x[ip]; d.pk = d.env.y[ip];
      const b = fig.st.b;
      fig.read([`envelope peak ${fmt(d.pk, 2)} g at ${fmt(d.fp, 2)} Hz`, `±${fmt(100 * b, 1)} %: plateau ${fmt(d.fp * (1 - b), 2)} to ${fmt(d.fp * (1 + b), 2)} Hz`,
        use.length > 1 ? "three cases: envelope first, then broadened" : "one case: broadened directly"]);
      void k;
    },
    layout(fig, W) {
      const h = W >= 600 ? 270 : 240;
      fig.L = {chart: {x: 44, y: 26, w: W - 44 - 14, h: h - 26 - 32}};
      return h;
    },
    draw(fig, p) {
      const d = fig.data, c = p.c;
      if (!d.br) return;
      const cyc = 10, ph = fig.playing || fig.t > 0 ? fig.t % cyc : cyc;
      const showEnv = ph > 2 || !fig.playing, bnow = !fig.playing ? fig.st.b : clamp((ph - 3.5) / 3.5, 0, 1) * fig.st.b;
      const ymax = Math.max(...[].concat(...d.use.map((i) => d.spec[i]))) * 1.15;
      const ax = axes(p, fig.L.chart, {x0: 0.5, x1: 50, xlog: true, y0: 0, y1: ymax, xlabel: "f (Hz)", ylabel: "spectral acceleration (g), 5 %"});
      const br = bnow === fig.st.b ? d.br : PH.broaden(d.env.x, d.env.y, bnow);
      if (bnow > 0) {
        p.rect(ax.X(d.fp * (1 - bnow)), ax.box.y, ax.X(d.fp * (1 + bnow)) - ax.X(d.fp * (1 - bnow)), ax.box.h, withAlpha(c.s1, 0.08));
      }
      for (const i of d.use) plot(p, ax, d.FS, d.spec[i], c[CASES[i].col], d.use.length > 1 ? 1.5 : 2);
      if (d.use.length > 1) {
        d.use.forEach((i, q) => {
          let ip = 0; d.spec[i].forEach((v, j) => { if (v > d.spec[i][ip]) ip = j; });
          p.text(CASES[i].name, ax.X(d.FS[ip]) + (q === 0 ? -6 : 6), ax.Y(d.spec[i][ip]) + (q === 1 ? -10 : 12), {size: 10, color: c.ink, align: q === 0 ? "right" : "left"});
        });
      }
      if (d.use.length === 1) {
        const i0 = d.FS.findIndex((f) => f >= d.fp * 0.8);
        p.text("computed", ax.X(d.FS[i0]) - 6, ax.Y(d.spec[1][i0]), {size: 10, color: c.ink, align: "right"});
      }
      if (showEnv && d.use.length > 1) plot(p, ax, d.env.x, d.env.y, c.ink, 1.2, [4, 3]);
      if (bnow > 0 || !fig.playing) plot(p, ax, br.x, br.y, c.s1, 2.4);
      if (bnow > 0) {
        p.line([[ax.X(d.fp * (1 - bnow)), ax.Y(d.pk) - 10], [ax.X(d.fp * (1 + bnow)), ax.Y(d.pk) - 10]], c.ink, 1);
        p.text(`±${fmt(100 * bnow, 0)} %`, ax.X(d.fp), ax.Y(d.pk) - 18, {size: 10.5, color: c.ink, align: "center"});
        p.text("broadened", ax.X(d.fp * (1 + bnow)) + 8, ax.Y(d.pk * 0.82), {size: 10.5, color: c.ink});
      }
      const stage = !fig.playing ? "" : ph < 2 ? "spectra of the cases" : ph < 3.5 ? (d.use.length > 1 ? "their envelope" : "the spectrum") : ph < 7 ? "peak broadening" : "design spectrum";
      if (stage) p.text(stage, ax.box.x + 8, ax.box.y + 10, {size: 10.5, weight: 600, color: c.muted});
    },
  };

  // ================================================================== figure: hysteresis and equivalent damping
  function hystModel(mode) {
    if (mode === "panel") {
      const P = PH.PANEL, bb = PH.bbcgen(P.vu, P.vcr, P.vcr / P.gcr);
      const bk = PH.polyline(bb.xs, bb.ys);
      return {bk, xref: P.gcr, xlo: P.gcr / 10, xhi: 0.02, xunit: (x) => fmt(100 * x, 4) + " %", xs: 100, xlabel: "shear strain γ (%)", ylabel: "panel shear V (kN)", ref: "γ_{cr}", F: "V"};
    }
    const rho = 19 / 9.81, G0 = rho * 250 * 250, gr = PH.sandRefStrain();       // kPa, strain in %
    const bk = PH.hyperbolic(G0 / 100, gr);                                     // tau (kPa) of strain in %
    return {bk, xref: gr, xlo: 1e-4, xhi: 1, xunit: (x) => fmtSig(x, 3) + " %", xs: 1, xlabel: "shear strain γ (%)", ylabel: "shear stress τ (kPa)", ref: "γ_{r}", F: "τ", G0, gr};
  }
  FIG.hysteresis = {
    title: "Backbone, hysteresis loop, secant stiffness and equivalent damping",
    alt: "A backbone curve with the stabilised Masing loop at the chosen amplitude, its secant line and areas, and the secant-stiffness ratio and equivalent damping against the amplitude.",
    animated: true,
    params: {mode: {def: "soil", values: ["soil", "panel"]}, a: {def: 2, min: 0.05, max: 200}},
    tex: ["F = F_r + 2F_{bb}\\!\\left(\\frac{x - x_r}{2}\\right),\\qquad K_\\text{sec} = \\frac{F(x_a)}{x_a}",
      "\\xi_h = \\frac{E_D}{4\\pi E_S},\\qquad E_D = 8\\int_0^{x_a} F_{bb}\\,dx - 4x_a F(x_a),\\qquad E_S = \\tfrac12 x_a F(x_a)",
      "\\tau = \\frac{G_{\\max}\\,\\gamma}{1 + \\gamma/\\gamma_r}"],
    init(fig) {
      fig.choice("mode", "Material", [["soil", "soil (hyperbolic)"], ["panel", "shear-wall panel (BBCGEN)"]]);
      fig.slider("a", "Amplitude x_{a} / reference", {log: true, fmt: (v) => "× " + fmtSig(v, 3)});
      fig.data.noteEl = el("div", {class: "lfig-note"});
      fig.root.insertBefore(fig.data.noteEl, fig.root.querySelector(".lfig-cap"));
    },
    change(fig) {
      const d = fig.data, M = hystModel(fig.st.mode);
      d.M = M;
      const xa = clamp(fig.st.a * M.xref, M.xlo, M.xhi);
      d.xa = xa; d.loop = PH.masing(M.bk, xa);
      d.curve = PH.logspace(M.xlo, M.xhi, 90).map((x) => { const l = PH.masing(M.bk, x); return [x, l.ratio, 100 * l.xi]; });
      const parts = [`x_{a} = ${M.xunit(xa)} (${fmtSig(xa / M.xref, 3)} ${M.ref})`, `K_{sec}/K_{el} = ${fmt(d.loop.ratio, 3)}`,
        fig.st.mode === "soil" ? `E_{D} = ${fmtSig(0.01 * d.loop.ED, 3)} kJ/m³, E_{S} = ${fmtSig(0.01 * d.loop.ES, 3)} kJ/m³`
          : `E_{D} = ${fmtSig(100 * d.loop.ED, 3)} kN·%, E_{S} = ${fmtSig(100 * d.loop.ES, 3)} kN·%`,
        `ξ_{h} = E_{D}/(4πE_{S}) = ${fmt(100 * d.loop.xi, 1)} %`];
      if (fig.st.mode === "soil") {
        const s = PH.SAND, lg = Math.log10(xa);
        const at = (arr) => (xa <= s.g[0] ? arr[0] : xa >= s.g[s.g.length - 1] ? arr[arr.length - 1] : PH.interpLin(s.g.map(Math.log10), arr, lg));
        parts.push(`library Sand at x_{a}: G/G_{max} ${fmt(at(s.G), 2)}, D ${fmt(at(s.D), 1)} %`);
      }
      fig.read(parts);
      d.noteEl.replaceChildren(richSpan(fig.st.mode === "soil"
        ? `Hyperbolic backbone with G_{max} of the lesson's sand (Vs = 250 m/s, 19 kN/m³: ${fmt(M.G0 / 1000, 0)} MPa) and γ_{r} = ${fmt(M.gr, 3)} %, where the SHAKE91 library Sand curve (points) has G/G_{max} = 0.5; Masing loops. SOIL's equivalent-linear iterations read G/G_{max} and damping from the library curves instead; Masing damping overshoots the measured curve at large strains.`
        : "The BBCGEN backbone of the lesson's wall panels (ACI 318-08: V_{u} = 12 472 kN, cracking at 0.3 V_{u} = 3742 kN, γ_{cr} = 8.31e-5, yield at 0.4 %), with loops by the Masing rule (GMR, used by Option NON for springs). The panels of lesson 9 use the Cheng-Mertz shear rules instead: their pinched loops dissipate much less (ξ_{h} about 10-12 % between γ = 2e-4 and 1e-2, against 28-50 % for these Masing loops; sassi/core/hysteresis.py)."));
    },
    layout(fig, W) {
      const pad = 8;
      if (W >= 640) {
        const w1 = Math.round(0.55 * W), h = 300;
        fig.L = {loop: {x: 60, y: pad + 22, w: w1 - 60 - 10, h: h - 2 * pad - 22 - 30},
          k: {x: w1 + 46, y: pad + 22, w: W - w1 - 46 - 12, h: (h - 2 * pad) / 2 - 22 - 26},
          xi: {x: w1 + 46, y: pad + (h - 2 * pad) / 2 + 22, w: W - w1 - 46 - 12, h: (h - 2 * pad) / 2 - 22 - 30}};
        return h;
      }
      const h1 = Math.min(300, W * 0.85), h2 = 130;
      fig.L = {loop: {x: 60, y: pad + 22, w: W - 60 - 12, h: h1 - pad - 22 - 30},
        k: {x: 46, y: h1 + 22, w: W - 46 - 12, h: h2 - 22 - 28}, xi: {x: 46, y: h1 + h2 + 22, w: W - 46 - 12, h: h2 - 22 - 30}};
      return h1 + 2 * h2;
    },
    draw(fig, p) {
      const d = fig.data, c = p.c, M = d.M, L = fig.L;
      if (!d.loop) return;
      const xa = d.xa, lp = d.loop, bk = M.bk, xs = M.xs;
      const xm = 1.2 * xa, fm = 1.25 * Math.max(Math.abs(bk.F(xm)), lp.Fa);
      const ax = axes(p, L.loop, {x0: -xm * xs, x1: xm * xs, y0: -fm, y1: fm, xlabel: M.xlabel, ylabel: M.ylabel});
      const X = (x) => ax.X(x * xs);
      // backbone
      const bb = [];
      for (let i = 0; i <= 120; i++) { const x = -xm + 2 * xm * i / 120; bb.push([X(x), ax.Y(bk.F(x))]); }
      p.line(bb, c.ref, 1.5);
      // loop: unloading from +xa to -xa, reloading back
      const loop = [];
      for (let i = 0; i <= 100; i++) { const x = xa - 2 * xa * i / 100; loop.push([X(x), ax.Y(lp.down(x))]); }
      for (let i = 0; i <= 100; i++) { const x = -xa + 2 * xa * i / 100; loop.push([X(x), ax.Y(lp.up(x))]); }
      p.poly(loop, withAlpha(c.s1, 0.16), null);
      p.poly([[X(0), ax.Y(0)], [X(xa), ax.Y(0)], [X(xa), ax.Y(lp.Fa)]], withAlpha(c.s2, 0.2), null);
      p.line(loop, c.s1, 2);
      p.line([[X(-xm), ax.Y(-lp.ksec * xm)], [X(xm), ax.Y(lp.ksec * xm)]], c.s2, 1.6, [5, 3]);
      const psi = TAU * 0.35 * fig.t, x = xa * Math.cos(psi), F = Math.sin(psi) >= 0 ? lp.down(x) : lp.up(x);
      p.circle(X(x), ax.Y(F), 5, c.s1, c.bg, 2);
      p.text("backbone", X(-0.98 * xm), ax.Y(bk.F(-0.98 * xm)) + 12, {size: 10, color: c.muted});
      p.text("secant K_{sec}", X(0.55 * xm), ax.Y(lp.ksec * 0.55 * xm) - 12, {size: 10, color: c.ink, align: "right"});
      p.text("E_{D} = loop area", X(-0.05 * xa), ax.Y(0.45 * lp.Fa), {size: 10, color: c.ink, align: "right"});
      p.text("E_{S}", X(0.78 * xa), ax.Y(0.22 * lp.Fa), {size: 10, color: c.ink, align: "center"});
      // amplitude curves (small multiples, one quantity each)
      const cx0 = d.curve[0][0], cx1 = d.curve[d.curve.length - 1][0];
      const ak = axes(p, L.k, {x0: cx0 * xs, x1: cx1 * xs, xlog: true, y0: 0, y1: 1.05, yt: [0, 0.5, 1], title: "K_{sec} / K_{el}", xlabel: ""});
      plot(p, ak, d.curve.map((q) => q[0] * xs), d.curve.map((q) => q[1]), c.s2, 2);
      const xmax = Math.max(30, ...d.curve.map((q) => q[2])) * 1.05;
      const axi = axes(p, L.xi, {x0: cx0 * xs, x1: cx1 * xs, xlog: true, y0: 0, y1: xmax, title: "ξ_{h} (%)", xlabel: M.xlabel});
      plot(p, axi, d.curve.map((q) => q[0] * xs), d.curve.map((q) => q[2]), c.s1, 2);
      if (fig.st.mode === "soil") {
        const s = PH.SAND;
        s.g.forEach((g, i) => {
          if (g < cx0 || g > cx1) return;
          p.circle(ak.X(g), ak.Y(s.G[i]), 3.2, c.bg, c.s3, 1.5);
          p.circle(axi.X(g), axi.Y(s.D[i]), 3.2, c.bg, c.s3, 1.5);
        });
        p.text("○ library Sand", ak.box.x + 6, ak.box.y + ak.box.h - 10, {size: 9.5, color: c.s3});
      }
      for (const [a2, v] of [[ak, lp.ratio], [axi, 100 * lp.xi]]) {
        p.line([[a2.X(xa * xs), a2.box.y], [a2.X(xa * xs), a2.box.y + a2.box.h]], c.ink, 1, [3, 3]);
        p.circle(a2.X(xa * xs), a2.Y(v), 4, a2 === ak ? c.s2 : c.s1, c.bg, 1.5);
      }
    },
  };

  // ================================================================== figure: frequency interpolation and CRITFREQ
  const DF = 1 / (8192 * 0.005);
  const FREQ_SETS = {
    lesson: [4, 20, 41, 61, 82, 102, 123, 143, 164, 184, 205, 225, 246, 287, 328, 369, 410, 492, 573, 655, 737, 819],
    coarse: [4, 41, 82, 123, 164, 205, 287, 369, 492, 655, 819],
  };
  const TRUE_MODES = {
    "2": [{f: 3.49, beta: 0.052, P: 1.34}, {f: 13.0, beta: 0.10, P: -0.34}],
    "3": [{f: 3.49, beta: 0.052, P: 1.0}, {f: 4.6, beta: 0.04, P: 0.34}, {f: 13.0, beta: 0.10, P: -0.34}],
  };
  FIG["freq-interpolation"] = {
    title: "Computed and interpolated transfer functions, and CRITFREQ",
    alt: "A true transfer function, the values computed at the SSI frequencies, the transfer function MOTION interpolates between them, and the peaks CRITFREQ flags.",
    animated: false,
    params: {modes: {def: "2", values: ["2", "3"]}, set: {def: "coarse", values: ["lesson", "coarse"]}},
    tex: ["H(\\omega) = \\frac{C_1\\,\\omega^4 + C_2\\,\\omega^2 + C_3}{\\omega^4 + C_4\\,\\omega^2 + C_5}",
      "100\\,\\frac{\\lvert A_I(f_p) - A_\\text{ref}\\rvert}{A_\\text{ref}} > \\text{tol}"],
    note: "True transfer function (grey): 1 + Σ P r²/(c(β) − r²) with the SSI mode of lesson 1 (3.49 Hz, 5.2 %, P = 1.34) and an illustrative mode at 13 Hz (P = −0.34); 'three modes' adds an illustrative close mode at 4.6 Hz. The dots are its values at the SSI frequencies f = nΔf (Δf = 0.0244 Hz); the blue curve is MOTION's interpolation, option 1 (five-point windows of the 2-DOF form, with the spurious-pole guard of Theory Section 10.1). CRITFREQ,5,50: peaks above 50 % of the maximum are compared with the larger of the two bracketing computed values. Click the plot to add a frequency there.",
    init(fig) {
      fig.choice("set", "Frequencies", [["lesson", "lesson 1 set (22)"], ["coarse", "coarse set (11, lesson 10)"]]);
      fig.choice("modes", "True TF", [["2", "two modes"], ["3", "three modes"]]);
      fig.button("Add the flagged frequencies (CRITFREQ)", () => {
        const add = (fig.data.crit || []).filter((q) => q.flag).map((q) => q.n);
        if (add.length) fig.set("added", Array.from(new Set((fig.st.added || []).concat(add))));
      }, "append the frequency numbers CRITFREQ flags, as lesson 10 does with FREQ,1,@ADDC[1]");
      fig.button("Reset", () => fig.set("added", []), "back to the chosen frequency set");
      if (!Array.isArray(fig.st.added)) fig.st.added = [];
      fig.drag(() => fig.data.ax && fig.data.ax.box, (x, y, ev) => {
        if (ev.type !== "pointerdown") return;
        const n = Math.round(fig.data.ax.ix(x) / DF);
        if (n >= 1 && n <= 819 && !fig.data.ns.includes(n)) fig.set("added", (fig.st.added || []).concat([n]));
      });
    },
    change(fig, k) {
      const d = fig.data;
      if (k === "set" || k === "modes") fig.st.added = [];
      const modes = TRUE_MODES[fig.st.modes];
      const ns = Array.from(new Set(FREQ_SETS[fig.st.set].concat(fig.st.added || []))).filter((n) => n >= 1 && n <= 819).sort((a, b) => a - b);
      d.ns = ns;
      const fu = ns.map((n) => n * DF), Hu = fu.map((f) => PH.modalTF(f, modes));
      const fi = [];
      for (let q = 0; q <= 819; q++) fi.push(q * DF);
      const Hi = PH.interpTF(fu, Hu, fi, [1, 0]);
      d.fu = fu; d.au = Hu.map(C.abs); d.fi = fi; d.ai = Hi.map(C.abs);
      d.ft = []; d.at = [];
      for (let f = 0; f <= 20.0001; f += 0.01) { d.ft.push(f); d.at.push(C.abs(PH.modalTF(f, modes))); }
      d.crit = PH.critfreq(d.fu, d.au, d.fi, d.ai, 5, 50, DF);
      const lines = [`${ns.length} computed frequencies` + ((fig.st.added || []).length ? ` (added: ${fig.st.added.join(", ")})` : "")];
      for (const q of d.crit) lines.push(`peak ${fmt(q.f, 2)} Hz (n = ${q.n}): interpolated ${fmt(q.amp, 2)}, computed neighbours ≤ ${fmt(q.ref, 2)}, ${fmt(q.d, 1)} % → ${q.flag ? "flagged" : "supported"}`);
      if (!d.crit.length) lines.push("no interpolated peak above 50 % of the maximum");
      fig.read(lines);
    },
    layout(fig, W) {
      const h = W >= 600 ? 300 : 270;
      fig.L = {chart: {x: 40, y: 26, w: W - 40 - 14, h: h - 26 - 40}};
      return h;
    },
    draw(fig, p) {
      const d = fig.data, c = p.c;
      if (!d.fi) return;
      const ymax = Math.ceil(Math.max(...d.at, ...d.ai) * 1.12);
      const ax = axes(p, fig.L.chart, {x0: 0, x1: 20, y0: 0, y1: ymax, xlabel: "f (Hz)", ylabel: "|H| (roof)"});
      d.ax = ax;
      plot(p, ax, d.ft, d.at, c.ref, 1.6);
      plot(p, ax, d.fi, d.ai, c.s1, 2);
      for (let i = 0; i < d.fu.length; i++) {
        p.circle(ax.X(d.fu[i]), ax.Y(d.au[i]), 4, c.s1, c.bg, 1.8);
        p.line([[ax.X(d.fu[i]), ax.box.y + ax.box.h], [ax.X(d.fu[i]), ax.box.y + ax.box.h - 6]], c.s1, 1.5);
      }
      for (const q of d.crit) {
        p.circle(ax.X(q.f), ax.Y(q.amp), 7, null, q.flag ? c.hi : c.ok, 2);
        p.text(q.flag ? `flagged: +${fmt(q.d, 0)} %` : "supported", ax.X(q.f) + 10, ax.Y(q.amp) - 4, {size: 10.5, color: q.flag ? c.hi : c.ok, weight: 600});
      }
      p.text("true", ax.X(15.5), ax.Y(PH.interpLin(d.ft, d.at, 15.5)) - 10, {size: 10, color: c.muted});
      p.text("● computed   — interpolated (TFI)   — true (grey)", ax.box.x + ax.box.w, ax.box.y + 8, {size: 9.5, color: c.muted, align: "right"});
    },
  };

  /** Public interface: learn.js calls fill(); mount() draws one figure into a host element; the registry and
   *  the physics functions are read by tests/unit/test_lesson_figures.py (Node). */
  S.Figures = Object.assign(S.Figures || {}, {registry: FIG, physics: PH, fill, mount, names: () => Object.keys(FIG),
    reducedMotion});
})(typeof SASSI !== "undefined" ? SASSI : (globalThis.SASSI = globalThis.SASSI || {}));
