/* SASSI-EDU explainer videos -- the player.
 *
 * A video is a list of scenes; a scene is a list of beats; a beat is one narrated sentence plus the
 * animation cues that go with it (authoring contract: docs/internal/explainer_videos.md).
 *
 *   SV.video({id, n, part, title, subtitle, next, pron, scenes: [
 *     {id, title, build(s), tick(s, t, dt), beats: [
 *       {say: "The springs [sp]lower the frequency.", go(k) {...}, sp(k) {...}, hold: ms, gap: ms}, ...]}]})
 *
 * The narration is spoken by the browser (Web Speech API): no audio files.  A beat starts its cues
 * (go), speaks its sentence and fires the [marker] cues when the voice reaches them (word-boundary
 * events, or the calibrated speaking rate when a voice sends none).  The next beat starts when the
 * sentence has been spoken and the beat's hold time has passed.  With the voice off the same timeline
 * runs on estimated durations (captions only).
 *
 * Seeking rebuilds the scene and replays the cues of the earlier beats instantly (k.instant = true), so
 * every cue must be able to jump to its end state: use the k helpers, which do.  The clock that drives
 * tweens, k.at timers and scene ticks stops while paused.  Vanilla JavaScript, no modules (works from
 * file:// and from the GUI server).  Debug: ?debug shows scene.beat; SV.player is the running player. */
"use strict";

(function () {
  const SV = (window.SV = window.SV || {});
  const W = 1920, H = 1080;
  const NS = "http://www.w3.org/2000/svg";
  SV.W = W; SV.H = H; SV.NS = NS;

  // ================================================================== small helpers
  const store = {
    get(key, dflt) {
      try { const v = window.localStorage.getItem("sassi-edu.video." + key); return v === null ? dflt : JSON.parse(v); } catch (e) { return dflt; }
    },
    set(key, value) {
      try { window.localStorage.setItem("sassi-edu.video." + key, JSON.stringify(value)); } catch (e) { /* storage unavailable */ }
    },
  };
  function h(tag, attrs, parent) {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v === undefined || v === null || v === false) continue;
      if (k === "text") e.textContent = v;
      else if (k === "html") e.innerHTML = v;
      else if (k === "style" && typeof v === "object") Object.assign(e.style, v);
      else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
      else e.setAttribute(k, v === true ? "" : v);
    }
    if (parent) parent.appendChild(e);
    return e;
  }
  const ICON = {
    play: '<svg viewBox="0 0 24 24"><path d="M7 4.5v15l13-7.5z"/></svg>',
    pause: '<svg viewBox="0 0 24 24"><path d="M6 4h4.5v16H6zM13.5 4H18v16h-4.5z"/></svg>',
    prev: '<svg viewBox="0 0 24 24"><path d="M6 5h2.5v14H6zM20 5v14L9.5 12z"/></svg>',
    next: '<svg viewBox="0 0 24 24"><path d="M15.5 5H18v14h-2.5zM4 5v14l10.5-7z"/></svg>',
    cc: '<svg viewBox="0 0 24 24"><path d="M3 5h18v14H3zm2 2v10h14V7zm2.5 2.5h4v1.6h-2.4v1.8h2.4v1.6h-4zm5.5 0h4v1.6h-2.4v1.8h2.4v1.6h-4z"/></svg>',
    list: '<svg viewBox="0 0 24 24"><path d="M4 6h16v2H4zm0 5h16v2H4zm0 5h10v2H4z"/></svg>',
    gear: '<svg viewBox="0 0 24 24"><path d="M12 8.5a3.5 3.5 0 1 0 0 7 3.5 3.5 0 0 0 0-7zm8.4 4.6.1-1.1-.1-1.1 2.1-1.6-2-3.5-2.5 1a7.6 7.6 0 0 0-1.9-1.1L15.7 3h-4l-.4 2.7a7.6 7.6 0 0 0-1.9 1.1l-2.5-1-2 3.5L7 10.9l-.1 1.1.1 1.1-2.1 1.6 2 3.5 2.5-1c.6.5 1.2.8 1.9 1.1l.4 2.7h4l.4-2.7c.7-.3 1.3-.6 1.9-1.1l2.5 1 2-3.5z"/></svg>',
    full: '<svg viewBox="0 0 24 24"><path d="M4 4h6v2H6v4H4zm10 0h6v6h-2V6h-4zM4 14h2v4h4v2H4zm14 0h2v6h-6v-2h4z"/></svg>',
    voice: '<svg viewBox="0 0 24 24"><path d="M4 9h4l5-4v14l-5-4H4zm12.5-1.5a6 6 0 0 1 0 9l-1.2-1.3a4.2 4.2 0 0 0 0-6.4z"/></svg>',
    mute: '<svg viewBox="0 0 24 24"><path d="M4 9h4l5-4v14l-5-4H4zm11.6.6 1.3-1.3 2 2 2-2 1.3 1.3-2 2 2 2-1.3 1.3-2-2-2 2-1.3-1.3 2-2z"/></svg>',
  };

  // ================================================================== easing
  const EASE = {
    linear: (t) => t,
    in: (t) => t * t * t,
    out: (t) => 1 - Math.pow(1 - t, 3),
    inOut: (t) => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2),
    back: (t) => { const c1 = 1.70158, c3 = c1 + 1; return 1 + c3 * Math.pow(t - 1, 3) + c1 * Math.pow(t - 1, 2); },
    sine: (t) => -(Math.cos(Math.PI * t) - 1) / 2,
  };
  SV.EASE = EASE;

  // ================================================================== narration text
  /** Spoken forms of the abbreviations and symbols of the course (captions keep the written form).
   *  A video adds its own with `pron`; a one-off is written {shown|spoken} in the text. */
  const PRON = [
    [/\bSASSI-EDU\b/g, "sassy E D U"], [/\bACS SASSI\b/g, "A C S sassy"], [/\bSASSI\b/g, "sassy"],
    [/\bSSI\b/g, "S S I"], [/\bISRS\b/g, "I S R S"], [/\bZPA\b/g, "Z P A"], [/\bSRSS\b/g, "S R S S"], [/\bSRP\b/g, "S R P"],
    [/\bFI-FSIN\b/g, "F I, F S I N"], [/\bFI-EVBN\b/g, "F I, E V B N"], [/\bFFV\b/g, "F F V"], [/\bFV\b/g, "F V"],
    [/\bMSM\b/g, "M S M"], [/\bSM\b/g, "S M"], [/\bSV\b/g, "S V"], [/\bSH\b/g, "S H"],
    [/\bANALYS\b/g, "analys"], [/\bRELDISP\b/g, "rel-disp"], [/\bEQUAKE\b/g, "E-quake"], [/\bLOADGEN\b/g, "load-gen"],
    [/\bINTGEN\b/g, "int-gen"], [/\bCRITFREQ\b/g, "crit-freq"], [/\bHARMFRAME\b/g, "harm-frame"], [/\bAFWRITE\b/g, "A F write"],
    [/\bNLSSIITER\b/g, "N L S S I iter"], [/\bSITEX\b/g, "site X"], [/\bFOUNSTIF\b/g, "found-stiff"], [/\bFOUNDASH\b/g, "found-dash"],
    [/\bFOUNDAMP\b/g, "found-damp"], [/\bFILE(\d+)\b/g, "file $1"], [/\bTFU\b/g, "T F U"], [/\bTFI\b/g, "T F I"], [/\bINP\b/g, "I N P"],
    [/\bAPDL\b/g, "A P D L"], [/\bCDB\b/g, "C D B"], [/\bNFFT\b/g, "N F F T"], [/\bDOFs\b/g, "degrees of freedom"], [/\bDOF\b/g, "degree of freedom"],
    [/\bVs\b/g, "V s"], [/\bVp\b/g, "V p"], [/\b1D\b/g, "one-D"], [/\b2D\b/g, "two-D"], [/\b3D\b/g, "three-D"],
    [/\bHz\b/g, "hertz"], [/\bft\/s\b/g, "feet per second"], [/\bkip·ft\b/g, "kip-feet"], [/\bksf\b/g, "K S F"], [/\bksi\b/g, "K S I"],
    [/\bkcf\b/g, "K C F"], [/\b1 ft\b/g, "1 foot"], [/(\d) ?ft\b/g, "$1 feet"], [/\bm\/s\b/g, "metres per second"], [/\bkN\b/g, "kilonewtons"], [/\bMN\b/g, "meganewtons"],
    [/\bMB\b/g, "megabytes"], [/\bGB\b/g, "gigabytes"], [/(\d) ?g\b/g, "$1 gee"], [/\b1 m\b/g, "1 metre"], [/(\d) ?m\b/g, "$1 metres"],
    [/(\d) ?mm\b/g, "$1 millimetres"], [/(\d) ?s\b/g, "$1 seconds"], [/\be\.g\./g, "for example"], [/\bi\.e\./g, "that is"], [/\bvs\./g, "versus"],
    [/Δf/g, "delta f"], [/Δ/g, "delta "], [/ω/g, "omega"], [/β/g, "beta"], [/ξ/g, "xi"], [/γ/g, "gamma"], [/λ/g, "lambda"], [/θ/g, "theta"],
    [/ = /g, " equals "], [/≈/g, " about "], [/→/g, " to "], [/±/g, " plus or minus "], [/×/g, " times "], [/−/g, "minus "], [/(\d)–(\d)/g, "$1 to $2"], [/[–—]/g, ", "],
  ];
  SV.PRON = PRON;
  /** FNV-1a hash (8 hex digits) of a spoken text: the key of its recorded clip in SV.AUDIO[id].beats. */
  SV.hash = function (text) {
    let x = 0x811c9dc5;
    for (let i = 0; i < text.length; i++) { x ^= text.charCodeAt(i); x = Math.imul(x, 0x01000193) >>> 0; }
    return ("0000000" + x.toString(16)).slice(-8);
  };
  const MARK = /\[([A-Za-z_]\w*)\]/g;
  const ALT = /\{([^{}|]*)\|([^{}]*)\}/g;
  /** Parse the say text of a beat: caption text, spoken text, markers (name + spoken char offset). */
  function parseSay(text, pron) {
    text = text || "";
    const parts = [];
    let last = 0, m;
    MARK.lastIndex = 0;
    while ((m = MARK.exec(text))) { parts.push({t: text.slice(last, m.index)}, {mark: m[1]}); last = m.index + m[0].length; }
    parts.push({t: text.slice(last)});
    let caption = "", spoken = "";
    const marks = [];
    for (const p of parts) {
      if (p.mark) { marks.push({name: p.mark, at: spoken.length, cap: caption.length}); continue; }
      caption += p.t.replace(ALT, "$1");
      let s = p.t.replace(ALT, "$2");
      for (const [re, rep] of pron) s = s.replace(re, rep);
      spoken += s;
    }
    caption = caption.replace(/\s+/g, " ").trim();
    return {caption, spoken: spoken.replace(/\s+/g, " "), marks};
  }
  SV.parseSay = parseSay;

  // ================================================================== voices
  const NOVELTY = /Albert|Bad News|Bahh|Bells|Boing|Bubbles|Cellos|Good News|Jester|Organ|Superstar|Trinoids|Whisper|Wobble|Zarvox|Fred|Junior|Ralph|Kathy|Grandma|Grandpa|Rocko|Shelley|Flo\b|Eddy|Reed|Sandy/i;
  function voiceScore(v) {
    if (!/^en([-_]|$)/i.test(v.lang)) return -1000;
    const n = v.name;
    let s = 0;
    if (NOVELTY.test(n)) s -= 500;
    if (/Natural|Neural|Premium|Enhanced/i.test(n)) s += 60;
    if (/Online/i.test(n)) s += 30;
    if (/Google US English|Google UK English/i.test(n)) s += 35;
    if (/Ava|Samantha|Allison|Susan|Zoe|Evan|Nathan|Tom|Aria|Jenny|Guy|Libby|Sonia|Ryan|Daniel|Serena|Karen|Moira|Alex|Matilda|Lee/i.test(n)) s += 20;
    if (/en[-_]US/i.test(v.lang)) s += 8; else if (/en[-_]GB/i.test(v.lang)) s += 6; else if (/en[-_](AU|CA|IE)/i.test(v.lang)) s += 3;
    if (v.localService) s += 4;
    return s;
  }
  const synth = window.speechSynthesis || null;
  function listVoices() {
    if (!synth) return [];
    return synth.getVoices().filter((v) => voiceScore(v) > -100).sort((a, b) => voiceScore(b) - voiceScore(a));
  }

  // ================================================================== the player
  const DEFAULT_CPS = 16.5;           // spoken characters per second at rate 1 before calibration
  const GAP = 120;                    // ms between beats
  const SCENE_GAP = 250;              // ms after the last beat of a scene
  const DEFAULT_RATE = 1.15;          // playback speed of the voice (pitch preserved)
  const TWEEN = 0.8;                  // every animation runs at this fraction of its written duration
  const DELAY = 0.85;                 // k.at delays, staggers and holds (they follow the faster voice)

  class Player {
    constructor(def) {
      this.def = def;
      this.pron = (def.pron || []).concat(PRON);
      this.scenes = def.scenes.map((sc, si) => ({
        def: sc, si,
        beats: sc.beats.map((b, bi) => Object.assign({si, bi, def: b}, parseSay(b.say, this.pron))),
      }));
      this.flat = [];
      this.scenes.forEach((sc) => sc.beats.forEach((b) => this.flat.push(b)));
      this.flat.forEach((b, i) => { b.i = i; });
      // recorded narration (aNN.js, web/voice_videos.mjs): one clip per beat, found by the hash of its text
      this.audio = (SV.AUDIO && SV.AUDIO[def.id]) || null;
      this.flat.forEach((b) => { const t = b.spoken.trim(); b.h = SV.hash(t); b.clip = (this.audio && t && this.audio.beats[b.h]) || null; });
      this.hasNatural = this.flat.some((b) => b.clip);
      this.natural = this.hasNatural && store.get("natural", true) !== false;
      this.clock = 0;                 // seconds, advances only while playing
      this.playing = false;
      this.cur = null;                // the running scene instance
      this.beat = null;               // the running beat state
      this.rate = store.get("rate2", DEFAULT_RATE);
      this.cps = store.get("cps", DEFAULT_CPS);
      this.captions = store.get("captions", true);
      this.voiceOn = store.get("voiceOn", true) && (!!synth || this.hasNatural);
      this.voiceName = store.get("voice", "");
      this.voice = null;
      this.ended = false;
      this.debug = /[?&]debug\b/.test(location.search);
      this.shot = /[?&]shot\b/.test(location.search);       // screenshots: the stage fills the window
    }

    // ---------------------------------------------------------------- durations
    estimate(b) {
      const clip = this.natural && b.clip;
      const speak = !b.spoken.trim() ? 0 : clip ? (clip.d / this.rate) * 1000 + 60 : (b.spoken.length / (this.cps * this.rate)) * 1000 + 250;
      return Math.max(speak, (b.def.hold || 0) * DELAY) + (b.def.gap !== undefined ? b.def.gap : GAP)
        + (b.bi === this.scenes[b.si].beats.length - 1 ? SCENE_GAP : 0);
    }
    timeline() {
      let t = 0;
      this.starts = this.flat.map((b) => { const s = t; t += this.estimate(b); return s; });
      this.total = t;
      return t;
    }

    // ---------------------------------------------------------------- DOM
    mount(host) {
      document.body.classList.add("sv-page");
      document.title = `${String(this.def.n).padStart(2, "0")} · ${this.def.title} — SASSI-EDU explainer`;
      const root = h("div", {class: "sv-player"}, host);
      const main = h("div", {class: "sv-main"}, root);
      this.viewport = h("div", {class: "sv-viewport"}, main);
      this.frame = h("div", {class: "sv-frame"}, this.viewport);
      this.stage = h("div", {class: "sv-stage", "aria-hidden": "true"}, this.frame);
      h("div", {class: "sv-brand", text: "SASSI-EDU · Lesson " + this.def.n}, this.stage);
      this.chapEl = h("div", {class: "sv-chap"}, this.stage);
      this.transcript = h("aside", {class: "sv-transcript", "aria-label": "Transcript"}, main);
      this.caption = h("div", {class: "sv-caption" + (this.captions ? "" : " off"), "aria-live": "off"}, root);
      this.capSpan = h("span", {}, this.caption);
      this.buildControls(root);
      this.buildOverlays();
      this.buildTranscript();
      if (this.debug) this.dbg = h("div", {class: "sv-debug"}, this.frame);
      this.root = root;
      new ResizeObserver(() => this.fit()).observe(this.viewport);
      this.fit();
      this.timeline();
      this.initVoices();
      this.updateTime();
      document.addEventListener("keydown", (ev) => this.key(ev));
      // frames stall in a hidden tab: a timer keeps the clock (and the narration's beats) going
      const at = this.parseHash();
      this.seek(at.si, at.bi, false);
      if (this.debug || this.shot) this.overlayStart.hidden = true;
      if (this.shot) root.classList.add("sv-shot");
      requestAnimationFrame((t) => this.frameLoop(t));
      setInterval(() => { const now = performance.now(); if (now - (this.lastNow || 0) > 220) this.advanceClock(now); }, 250);
      return this;
    }
    parseHash() {
      const m = /^#(\d+)(?:\.(\d+))?$/.exec(location.hash || "");
      if (m) {
        const si = Math.min(+m[1], this.scenes.length - 1);
        return {si, bi: Math.min(+(m[2] || 0), this.scenes[si].beats.length - 1)};
      }
      const id = (location.hash || "").slice(1);
      const si = this.scenes.findIndex((s) => s.def.id === id);
      return {si: Math.max(0, si), bi: 0};
    }
    fit() {
      const r = this.viewport.getBoundingClientRect();
      const pad = r.width > 700 && !this.shot ? 16 : 0;
      const s = Math.max(0.05, Math.min((r.width - pad * 2) / W, (r.height - pad * 2) / H));
      this.scale = s;
      this.frame.style.width = W * s + "px";
      this.frame.style.height = H * s + "px";
      this.stage.style.transform = `scale(${s})`;
    }
    buildControls(root) {
      const c = h("div", {class: "sv-controls"}, root);
      const B = (icon, title, fn) => { const b = h("button", {class: "sv-btn", title, "aria-label": title, html: ICON[icon], onclick: fn}, c); return b; };
      this.bPlay = B("play", "Play (space)", () => this.toggle());
      B("prev", "Previous sentence (←); previous chapter (shift ←)", () => this.step(-1));
      B("next", "Next sentence (→); next chapter (shift →)", () => this.step(1));
      this.timeEl = h("span", {class: "sv-time"}, c);
      const bar = h("div", {class: "sv-bar", title: ""}, c);
      h("div", {class: "sv-bar-track"}, bar);
      this.fill = h("div", {class: "sv-bar-fill"}, bar);
      this.ticks = h("div", {}, bar);
      this.knob = h("div", {class: "sv-bar-knob"}, bar);
      this.tip = h("div", {class: "sv-bar-tip"}, bar);
      this.bar = bar;
      bar.addEventListener("click", (ev) => {
        const r = bar.getBoundingClientRect();
        const t = ((ev.clientX - r.left) / r.width) * this.total;
        let i = 0;
        while (i + 1 < this.flat.length && this.starts[i + 1] <= t) i++;
        const b = this.flat[i];
        this.seek(b.si, b.bi, this.playing);
      });
      bar.addEventListener("mousemove", (ev) => {
        const r = bar.getBoundingClientRect();
        const x = ev.clientX - r.left, t = (x / r.width) * this.total;
        let i = 0;
        while (i + 1 < this.flat.length && this.starts[i + 1] <= t) i++;
        this.tip.style.left = x + "px";
        this.tip.textContent = `${fmt(t)} · ${this.scenes[this.flat[i].si].def.title || ""}`;
      });
      this.bVoice = B(this.voiceOn ? "voice" : "mute", "Narration voice on/off (m)", () => this.setVoiceOn(!this.voiceOn));
      this.bCC = B("cc", "Captions (c)", () => this.setCaptions(!this.captions));
      this.bCC.classList.toggle("on", this.captions);
      this.bTr = B("list", "Transcript and chapters (t)", () => this.toggleTranscript());
      this.bSet = B("gear", "Voice and speed", () => this.menu.classList.toggle("open"));
      B("full", "Full screen (f)", () => this.fullscreen());
      // settings menu
      const m = h("div", {class: "sv-menu"}, root);
      this.menu = m;
      h("label", {text: "Narration voice"}, m);
      this.voiceSel = h("select", {onchange: () => {
        const v = this.voiceSel.value;
        this.pickVoice(v);
        if (v !== "__natural__") store.set("voice", v);
        this.timeline(); this.drawTicks(); this.restartBeatSpeech();
      }}, m);
      h("label", {text: "Speed"}, m);
      this.rateSel = h("select", {onchange: () => { this.rate = +this.rateSel.value; store.set("rate2", this.rate); this.timeline(); this.drawTicks(); this.restartBeatSpeech(); }}, m);
      for (const r of [1, 1.15, 1.25, 1.4, 1.6]) h("option", {value: r, text: r + "×", selected: r === this.rate}, this.rateSel);
      this.voiceNote = h("div", {class: "sv-note"}, m);
      document.addEventListener("click", (ev) => { if (!m.contains(ev.target) && !this.bSet.contains(ev.target)) m.classList.remove("open"); });
    }
    drawTicks() {
      this.ticks.textContent = "";
      let i = 0;
      for (const sc of this.scenes) {
        if (i > 0) h("div", {class: "sv-bar-tick", style: {left: (100 * this.starts[i]) / this.total + "%"}}, this.ticks);
        i += sc.beats.length;
      }
    }
    buildOverlays() {
      const d = this.def;
      const o = h("div", {class: "sv-overlay"}, this.frame);
      h("div", {class: "sv-o-kicker", text: `Lesson ${d.n} · ${d.part || ""} · explainer`}, o);
      h("div", {class: "sv-o-title", text: d.title}, o);
      if (d.subtitle) h("div", {class: "sv-o-sub", text: d.subtitle}, o);
      h("button", {class: "sv-o-play", "aria-label": "Play", html: ICON.play, onclick: () => this.play()}, o);
      this.oLen = h("div", {class: "sv-o-voice"}, o);
      this.oVoice = h("div", {class: "sv-o-voice"}, o);
      this.overlayStart = o;
      const e = h("div", {class: "sv-overlay", hidden: true}, this.frame);
      h("div", {class: "sv-o-kicker", text: `End of lesson ${d.n} explainer`}, e);
      h("div", {class: "sv-o-title", text: d.next ? `Up next: ${d.next.title}` : "That's the course."}, e);
      const row = h("div", {class: "sv-o-row"}, e);
      if (d.next) h("a", {class: "sv-o-btn primary", href: d.next.href, text: `Lesson ${d.n + 1} explainer →`}, row);
      h("button", {class: "sv-o-btn", text: "↺ Replay", onclick: () => { e.hidden = true; this.seek(0, 0, true); }}, row);
      h("a", {class: "sv-o-btn", href: "index.html", text: "All explainers"}, row);
      this.overlayEnd = e;
    }
    buildTranscript() {
      const t = this.transcript;
      this.trLines = [];
      h("div", {class: "sv-kicker", style: {fontSize: "12px", letterSpacing: ".12em"}, text: "Transcript"}, t);
      for (const sc of this.scenes) {
        h("h3", {text: sc.def.title || sc.def.id}, t);
        for (const b of sc.beats) {
          if (!b.caption) continue;
          const p = h("p", {text: b.caption, onclick: () => this.seek(b.si, b.bi, true)}, t);
          this.trLines[b.i] = p;
        }
      }
    }
    toggleTranscript() {
      const on = !this.transcript.classList.contains("open");
      this.transcript.classList.toggle("open", on);
      this.bTr.classList.toggle("on", on);
      setTimeout(() => this.fit(), 0);
    }
    setCaptions(on) {
      this.captions = on;
      store.set("captions", on);
      this.caption.classList.toggle("off", !on);
      this.bCC.classList.toggle("on", on);
      setTimeout(() => this.fit(), 0);
    }
    setVoiceOn(on) {
      if (on && !synth && !this.hasNatural) return;
      this.voiceOn = on;
      store.set("voiceOn", on);
      this.bVoice.innerHTML = on ? ICON.voice : ICON.mute;
      this.restartBeatSpeech();
      this.describeVoice();
    }
    fullscreen() {
      const r = this.root;
      if (document.fullscreenElement) document.exitFullscreen();
      else if (r.requestFullscreen) r.requestFullscreen().catch(() => {});
      else if (r.webkitRequestFullscreen) r.webkitRequestFullscreen();
    }
    key(ev) {
      if (ev.target && /INPUT|SELECT|TEXTAREA/.test(ev.target.tagName)) return;
      const k = ev.key;
      if (k === " " || k === "k") { ev.preventDefault(); this.toggle(); }
      else if (k === "ArrowRight" || k === "l") { ev.preventDefault(); ev.shiftKey ? this.chapter(1) : this.step(1); }
      else if (k === "ArrowLeft" || k === "j") { ev.preventDefault(); ev.shiftKey ? this.chapter(-1) : this.step(-1); }
      else if (k === "c") this.setCaptions(!this.captions);
      else if (k === "m") this.setVoiceOn(!this.voiceOn);
      else if (k === "t") this.toggleTranscript();
      else if (k === "f") this.fullscreen();
    }

    // ---------------------------------------------------------------- voices
    initVoices() {
      const load = () => {
        this.voices = listVoices();
        this.voiceSel.textContent = "";
        if (this.hasNatural) h("option", {value: "__natural__", text: `${this.audio.voice || "Recorded"} · natural voice (ElevenLabs)`}, this.voiceSel);
        for (const v of this.voices) h("option", {value: v.name, text: `${v.name} (${v.lang})${v.localService ? "" : " · online"} · browser voice`}, this.voiceSel);
        this.pickVoice(this.natural ? "__natural__" : this.voiceName);
        this.describeVoice();
      };
      load();
      if (synth && "onvoiceschanged" in synth) synth.addEventListener("voiceschanged", load);
    }
    pickVoice(name) {
      const vs = this.voices || [];
      // the recorded voice, or a browser voice (also the fallback for a beat without a recording)
      this.natural = name === "__natural__" && this.hasNatural;
      if (this.hasNatural) store.set("natural", this.natural);
      this.voice = vs.find((v) => v.name === (this.natural ? this.voiceName : name)) || vs[0] || null;
      this.voiceSel.value = this.natural ? "__natural__" : (this.voice ? this.voice.name : "");
      if (!this.natural) this.voiceName = this.voice ? this.voice.name : "";
      this.describeVoice();
    }
    describeVoice() {
      const len = `About ${Math.round(this.total / 60000)} minutes · ${this.scenes.length} chapters`;
      if (this.oLen) this.oLen.textContent = len;
      let msg;
      if (this.natural && this.voiceOn) msg = `Narrated by ${this.audio.voice || "a recorded voice"} (ElevenLabs). Voice and speed under ⚙.`;
      else if (!synth && !this.hasNatural) msg = "This browser has no text-to-speech: the video plays with captions.";
      else if (!this.voiceOn) msg = "Narration off: captions only (press m to turn the voice on).";
      else if (!this.voice) msg = "Looking for an English voice…";
      else msg = `Narrated by your browser's voice “${this.voice.name}”. Change it under ⚙.`;
      if (this.oVoice) this.oVoice.textContent = msg;
      if (this.voiceNote) {
        this.voiceNote.innerHTML = (this.hasNatural ? "The narration was recorded with ElevenLabs text-to-speech; you can also choose one of your browser's voices. " : "The narration is read by your browser's text-to-speech. ")
          + "Among browser voices, the ones marked Natural, Neural, Premium or Enhanced sound best (Edge and Chrome offer some; "
          + "on a Mac, System Settings › Accessibility › Spoken Content › System Voice › Manage Voices adds Premium voices).";
      }
    }

    // ---------------------------------------------------------------- transport
    toggle() { this.playing ? this.pause() : this.play(); }
    play() {
      this.overlayStart.hidden = true;
      this.menu.classList.remove("open");
      if (this.ended) { this.overlayEnd.hidden = true; this.seek(0, 0, true); return; }
      if (this.playing) return;
      this.playing = true;
      this.bPlay.innerHTML = ICON.pause;
      this.bPlay.title = "Pause (space)";
      if (this.beat && !this.beat.started) this.startBeat();
      else if (this.beat) this.resumeBeat();
    }
    pause() {
      if (!this.playing) return;
      this.playing = false;
      this.bPlay.innerHTML = ICON.play;
      this.bPlay.title = "Play (space)";
      if (this.beat) this.suspendSpeech();
    }
    step(d) {
      const i = this.beat ? this.beat.b.i : 0;
      const j = Math.max(0, Math.min(this.flat.length - 1, i + d));
      const b = this.flat[j];
      this.seek(b.si, b.bi, this.playing);
    }
    chapter(d) {
      const si = this.beat ? this.beat.b.si : 0;
      const sj = Math.max(0, Math.min(this.scenes.length - 1, si + d));
      this.seek(d < 0 && this.beat && this.beat.b.bi > 0 ? si : sj, 0, this.playing);
    }

    /** Show scene si at the start of beat bi (earlier beats applied instantly); play it if `play`. */
    seek(si, bi, play) {
      si = Math.max(0, Math.min(this.scenes.length - 1, si | 0));
      bi = Math.max(0, Math.min(this.scenes[si].beats.length - 1, bi | 0));
      this.ended = false;
      this.overlayEnd.hidden = true;
      this.stopSpeech();
      this.enterScene(si, false);
      const sc = this.scenes[si];
      this.stage.classList.add("sv-instant");
      for (let j = 0; j < bi; j++) this.runBeatInstant(sc.beats[j]);
      this.cur.tickNow();
      void this.stage.offsetWidth;
      this.stage.classList.remove("sv-instant");
      this.prepareBeat(sc.beats[bi]);
      if (play) { if (!this.playing) this.play(); else this.startBeat(); }
      else this.renderBeatInfo();
    }
    enterScene(si, fade) {
      const old = this.cur;
      if (old) {
        old.dispose();
        if (fade) { old.root.classList.add("sv-leaving"); setTimeout(() => old.root.remove(), 500); }
        else old.root.remove();
      }
      this.cur = new SceneRun(this, this.scenes[si]);
      this.chapEl.textContent = this.scenes[si].def.title || "";
      if (fade) {
        this.cur.root.style.opacity = "0";
        void this.cur.root.offsetWidth;
        this.cur.root.style.opacity = "";
      }
    }
    /** Debug / screenshots: scene si as it is at the end of beat bi (every cue applied, no animation). */
    snap(si, bi) {
      this.seek(si, Math.min(bi, this.scenes[si].beats.length - 1), false);
      this.stage.classList.add("sv-instant");
      this.runBeatInstant(this.beat.b);
      this.cur.tickNow();
      void this.stage.offsetWidth;
      this.stage.classList.remove("sv-instant");
      return {scenes: this.scenes.length, beats: this.scenes[si].beats.length};
    }
    runBeatInstant(b) {
      const k = this.cur.ctx(true, b);
      if (b.def.go) b.def.go.call(b.def, k);
      for (const m of b.marks) if (typeof b.def[m.name] === "function") b.def[m.name].call(b.def, k);
      this.cur.flushTimers();
    }
    prepareBeat(b) {
      this.beat = {b, started: false, fired: new Set(), spokeTo: 0, speechDone: false, t0: 0, token: 0, utter: null};
      if (this.natural && b.clip && !(this.preloaded && this.preloaded.h === b.h)) {
        const a = new Audio(b.clip.f);                 // ready before Play (or the next beat) asks for it
        a.preload = "auto";
        this.preloaded = {h: b.h, a};
      }
      this.showCaption(b);
      this.highlightTranscript(b);
    }
    renderBeatInfo() {
      this.updateTime();
      if (this.dbg && this.beat) this.dbg.textContent = `${this.beat.b.si}.${this.beat.b.bi}  ${this.scenes[this.beat.b.si].def.id}`;
    }
    showCaption(b) {
      this.capSpan.textContent = b.caption || "";
    }
    highlightTranscript(b) {
      if (this.trNow) this.trNow.classList.remove("now");
      this.trNow = this.trLines[b.i] || null;
      if (this.trNow) {
        this.trNow.classList.add("now");
        if (this.transcript.classList.contains("open")) this.trNow.scrollIntoView({block: "nearest"});
      }
    }
    startBeat() {
      const st = this.beat, b = st.b;
      st.started = true;
      st.t0 = this.clock;
      st.k = this.cur.ctx(false, b);
      this.renderBeatInfo();
      if (b.def.go) b.def.go.call(b.def, st.k);
      this.speak(0);
    }
    fire(name) {
      const st = this.beat;
      if (!st || st.fired.has(name)) return;
      st.fired.add(name);
      const fn = st.b.def[name];
      if (typeof fn === "function") fn.call(st.b.def, st.k);
    }
    /** Fire the markers up to spoken offset `upTo` (inclusive), in order. */
    fireUpTo(upTo) {
      for (const m of this.beat.b.marks) if (m.at <= upTo) this.fire(m.name);
    }

    // ---------------------------------------------------------------- speech
    speak(from) {
      const st = this.beat, b = st.b;
      const text = b.spoken.slice(from);
      st.token++;
      const token = st.token;
      st.from = from;
      st.speakStart = this.clock;
      st.speechDone = false;
      // marker fallback timers on the clock (a voice without boundary events, or captions only)
      st.markTimes = b.marks.filter((m) => m.at >= from && !st.fired.has(m.name))
        .map((m) => ({name: m.name, t: this.clock + Math.max(0, m.at - from) / (this.cps * this.rate)}));
      st.boundaries = false;
      if (!text.trim()) { st.speechDone = true; st.estEnd = this.clock; this.fireUpTo(Infinity); return; }
      if (this.voiceOn && this.natural && b.clip && !b.clipFailed) { this.playClip(st); return; }
      st.estEnd = this.clock + text.length / (this.cps * this.rate) + 0.25;
      if (!this.voiceOn || !synth || !this.voice) return;   // timed by the clock
      const u = new SpeechSynthesisUtterance(text);
      u.voice = this.voice;
      u.lang = this.voice.lang;
      u.rate = this.rate;
      u.onstart = () => { if (token === st.token && this.beat === st) { st.spoke = true; st.realStart = performance.now(); } };
      u.onboundary = (ev) => {
        if (token !== st.token || this.beat !== st) return;
        st.boundaries = true;
        st.spokeTo = from + ev.charIndex;
        this.fireUpTo(st.spokeTo + Math.max(1, ev.charLength || 1));
        // markers just ahead of the word being spoken
        this.fireUpTo(st.spokeTo);
      };
      const done = () => {
        if (token !== st.token || this.beat !== st) return;
        if (st.realStart && text.length > 40) {
          const cps = text.length / ((performance.now() - st.realStart) / 1000) / this.rate;
          if (cps > 6 && cps < 30) { this.cps = 0.75 * this.cps + 0.25 * cps; store.set("cps", this.cps); }
        }
        st.speechDone = true;
        this.fireUpTo(Infinity);
      };
      u.onend = done;
      u.onerror = (ev) => { if (ev && ev.error !== "interrupted" && ev.error !== "canceled") done(); };
      st.utter = u;            // keep a reference (Chrome drops the events of collected utterances)
      st.spoke = false;
      try {
        if (synth.speaking || synth.pending) synth.cancel();
        synth.speak(u);
        if (synth.paused) synth.resume();
      } catch (e) { st.utter = null; }
    }
    /** Play the recorded clip of the beat (from its start); its markers fire at their recorded times. */
    playClip(st) {
      const b = st.b, clip = b.clip;
      const token = st.token;
      st.markTimes = null;
      st.clipMarks = b.marks.filter((m) => !st.fired.has(m.name)).map((m) => ({name: m.name,
        t: clip.m && clip.m[m.name] !== undefined ? clip.m[m.name] : (m.at / Math.max(1, b.spoken.length)) * clip.d}));
      st.estEnd = this.clock + clip.d / this.rate + 0.3;
      const a = this.takeAudio(b);
      st.audio = a;
      a.playbackRate = this.rate;
      a.preservesPitch = true;
      const fail = () => {
        if (token !== st.token || this.beat !== st) return;
        st.audio = null;
        b.clipFailed = true;                         // this beat falls back to the browser voice
        this.speak(0);
      };
      a.onended = () => { if (token === st.token && this.beat === st && st.audio === a) { st.speechDone = true; this.fireUpTo(Infinity); } };
      a.onerror = fail;
      try { if (a.currentTime) a.currentTime = 0; } catch (e) { /* not loaded yet */ }
      const p = a.play();
      if (p && p.catch) p.catch((e) => {
        if (token !== st.token || this.beat !== st || st.audio !== a) return;
        if (e && e.name === "NotAllowedError") { st.audio = null; this.pause(); this.overlayStart.hidden = false; }   // autoplay blocked
        else if (e && e.name !== "AbortError") fail();
      });
      this.preloadNext(b);
    }
    takeAudio(b) {
      if (this.preloaded && this.preloaded.h === b.h) { const a = this.preloaded.a; this.preloaded = null; return a; }
      const a = new Audio(b.clip.f);
      a.preload = "auto";
      return a;
    }
    preloadNext(b) {
      const n = this.flat[b.i + 1];
      if (!n || !n.clip || n.clipFailed) return;
      const a = new Audio(n.clip.f);
      a.preload = "auto";
      this.preloaded = {h: n.h, a};
    }
    suspendSpeech() {
      const st = this.beat;
      if (!st || !st.started) return;
      if (st.audio) { st.audio.pause(); st.resumeFrom = st.speechDone ? null : "audio"; return; }
      st.token++;
      if (synth) synth.cancel();
      st.resumeFrom = st.speechDone ? null : this.resumePoint(st);
    }
    resumePoint(st) {
      const b = st.b;
      let at = st.from || 0;
      if (st.boundaries) at = st.spokeTo;
      else if (st.utter && st.realStart) at = Math.min(b.spoken.length, (st.from || 0) + Math.floor(((performance.now() - st.realStart) / 1000) * this.cps * this.rate));
      else at = Math.min(b.spoken.length, (st.from || 0) + Math.floor((this.clock - st.speakStart) * this.cps * this.rate));
      // back up to the start of the sentence part (after the last comma or full stop), or of the word
      const s = b.spoken;
      let k = at;
      while (k > (st.from || 0) && !/[,.;:]/.test(s[k - 1])) k--;
      if (at - k > 60) { k = at; while (k > 0 && s[k - 1] !== " ") k--; }
      return k;
    }
    resumeBeat() {
      const st = this.beat;
      if (st.resumeFrom === null || st.resumeFrom === undefined) return;   // the speech had finished
      if (st.resumeFrom === "audio" && st.audio) {
        st.audio.playbackRate = this.rate;
        st.estEnd = this.clock + Math.max(0, (st.audio.duration || st.b.clip.d) - st.audio.currentTime) / this.rate + 0.3;
        const p = st.audio.play();
        if (p && p.catch) p.catch(() => {});
        return;
      }
      this.speak(st.resumeFrom);
    }
    /** After a change of voice, speed or voice on/off in the middle of a beat. */
    restartBeatSpeech() {
      const st = this.beat;
      if (!st || !st.started || st.speechDone) return;
      const wantClip = this.voiceOn && this.natural && st.b.clip && !st.b.clipFailed;
      if (st.audio && wantClip) { st.audio.playbackRate = this.rate; return; }          // a new speed only
      const at = st.audio ? 0 : this.resumePoint(st);
      if (st.audio) { st.audio.pause(); st.audio = null; }
      st.token++;
      if (synth) synth.cancel();
      if (!this.playing) { st.resumeFrom = at; return; }
      setTimeout(() => { if (this.beat === st && this.playing) this.speak(at); }, 60);
    }
    stopSpeech() {
      if (this.beat) {
        this.beat.token++;
        if (this.beat.audio) { this.beat.audio.pause(); this.beat.audio = null; }
      }
      if (synth && (synth.speaking || synth.pending)) synth.cancel();
    }

    // ---------------------------------------------------------------- frame loop
    frameLoop(now) {
      this.advanceClock(now);
      requestAnimationFrame((t) => this.frameLoop(t));
    }
    advanceClock(now) {
      const dt = this.lastNow ? Math.min(0.3, Math.max(0, now - this.lastNow) / 1000) : 0;
      this.lastNow = now;
      if (this.playing) {
        this.clock += dt;
        this.cur.update(dt);
        this.checkBeat();
        this.updateTime();
      }
    }
    checkBeat() {
      const st = this.beat;
      if (!st || !st.started) return;
      if (st.audio && !st.speechDone) {
        // recorded clip: markers at their recorded times; a stalled clip is given up after a while
        const t = st.audio.currentTime;
        for (const m of st.clipMarks) if (!st.fired.has(m.name) && t >= m.t - 0.04) this.fire(m.name);
        if (this.clock > st.estEnd + 5) { st.speechDone = true; this.fireUpTo(Infinity); }
        return;
      }
      // marker timers (always armed: a voice may skip boundary events for some words)
      if (st.markTimes) for (const m of st.markTimes) {
        if (this.clock >= m.t && !st.fired.has(m.name)) {
          // with live boundary events, wait for the voice unless it is far behind the estimate
          if (!st.boundaries || this.clock >= m.t + 1.5) this.fire(m.name);
        }
      }
      const usingVoice = this.voiceOn && synth && this.voice && st.utter;
      if (!st.speechDone) {
        if (!usingVoice && this.clock >= st.estEnd) { st.speechDone = true; this.fireUpTo(Infinity); }
        // watchdog: a voice that never starts or never ends
        else if (usingVoice && !st.spoke && this.clock - st.speakStart > 3) { st.utter = null; }
        else if (usingVoice && this.clock > st.estEnd + Math.max(6, (st.estEnd - st.speakStart) * 1.5)) { st.speechDone = true; this.fireUpTo(Infinity); }
      }
      if (!st.speechDone) return;
      if (!st.doneAt) st.doneAt = this.clock;
      const d = st.b.def;
      const hold = ((d.hold || 0) * DELAY) / 1000;
      const last = st.b.bi === this.scenes[st.b.si].beats.length - 1;
      const gap = (d.gap !== undefined ? d.gap : GAP) / 1000 + (last ? SCENE_GAP / 1000 : 0);
      if (this.clock >= Math.max(st.t0 + hold, st.doneAt + gap)) this.advance();
    }
    advance() {
      const b = this.beat.b;
      const sc = this.scenes[b.si];
      if (b.bi + 1 < sc.beats.length) {
        this.prepareBeat(sc.beats[b.bi + 1]);
        this.startBeat();
        return;
      }
      if (b.si + 1 < this.scenes.length) {
        this.enterScene(b.si + 1, true);
        this.prepareBeat(this.scenes[b.si + 1].beats[0]);
        this.startBeat();
        return;
      }
      this.ended = true;
      this.playing = false;
      this.bPlay.innerHTML = ICON.play;
      this.overlayEnd.hidden = false;
      this.capSpan.textContent = "";
    }
    updateTime() {
      if (!this.total) return;
      let t = 0;
      const st = this.beat;
      if (st) {
        t = this.starts[st.b.i];
        if (st.started) t += Math.min(this.estimate(st.b), (this.clock - st.t0) * 1000);
      }
      const f = Math.min(1, t / this.total);
      this.fill.style.width = 100 * f + "%";
      this.knob.style.left = 100 * f + "%";
      this.timeEl.textContent = `${fmt(t)} / ${fmt(this.total)}`;
      if (!this.ticksDrawn) { this.drawTicks(); this.ticksDrawn = true; }
    }
  }
  function fmt(ms) {
    const s = Math.max(0, Math.round(ms / 1000));
    return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, "0")}`;
  }

  // ================================================================== a running scene
  class SceneRun {
    constructor(player, sc) {
      this.player = player;
      this.sc = sc;
      this.t = 0;                         // scene time (s)
      this.tweens = [];
      this.timers = [];
      this.loops = [];
      this.root = h("div", {class: "sv-scene"});
      player.stage.insertBefore(this.root, player.chapEl);
      this.svg = document.createElementNS(NS, "svg");
      this.svg.setAttribute("class", "sv-svg");
      this.svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
      this.root.appendChild(this.svg);
      this.s = {root: this.root, svg: this.svg, scene: this, $: (sel) => this.root.querySelector(sel), $$: (sel) => Array.from(this.root.querySelectorAll(sel))};
      if (sc.def.build) sc.def.build.call(sc.def, this.s);
      this.prepDraw(this.root);
    }
    /** Elements with class sv-draw start undrawn (stroke dash = length). */
    prepDraw(root) {
      for (const e of root.querySelectorAll(".sv-draw")) {
        if (e.__svLen) continue;
        let len = 0;
        try { len = e.getTotalLength(); } catch (err) { len = 0; }
        len = Math.ceil(len) + 2;
        e.__svLen = len;
        e.style.strokeDasharray = `${len} ${len}`;
        e.style.strokeDashoffset = String(len);
      }
    }
    dispose() { this.tweens = []; this.timers = []; this.loops = []; this.disposed = true; }
    update(dt) {
      this.t += dt;
      const clock = this.player.clock;
      // timers (k.at)
      if (this.timers.length) {
        const due = this.timers.filter((tm) => clock >= tm.t);
        this.timers = this.timers.filter((tm) => clock < tm.t);
        for (const tm of due) tm.fn();
      }
      this.stepTweens(clock);
      for (const fn of this.loops) fn(this.t, dt);
      if (this.sc.def.tick) this.sc.def.tick.call(this.sc.def, this.s, this.t, dt);
    }
    tickNow() {
      this.stepTweens(this.player.clock);
      for (const fn of this.loops) fn(this.t, 0);
      if (this.sc.def.tick) this.sc.def.tick.call(this.sc.def, this.s, this.t, 0);
    }
    stepTweens(clock) {
      if (!this.tweens.length) return;
      for (const tw of this.tweens.slice()) {
        if (tw.killed || tw.finished) continue;
        const u = tw.dur <= 0 ? 1 : Math.min(1, Math.max(0, (clock - tw.t0) / tw.dur));
        if (u <= 0 && clock < tw.t0) continue;          // delayed, not started
        tw.apply(tw.ease(u), u);
        if (u >= 1) { tw.finished = true; if (tw.done) tw.done(); }
      }
      this.tweens = this.tweens.filter((tw) => !tw.killed && !tw.finished);
    }
    flushTimers() {
      // instant mode: run every pending timer now, in time order (they may add more)
      let guard = 0;
      while (this.timers.length && guard++ < 1000) {
        this.timers.sort((a, b) => a.t - b.t);
        const tm = this.timers.shift();
        tm.fn();
      }
      let guard2 = 0;
      while (this.tweens.length && guard2++ < 50) {
        const list = this.tweens;
        this.tweens = [];
        for (const tw of list) if (!tw.killed) { tw.apply(1, 1); if (tw.done) tw.done(); }
      }
    }
    /** The cue context of a beat. */
    ctx(instant, b) {
      const run = this, player = this.player;
      const els = (target) => {
        if (!target) return [];
        if (typeof target === "string") return Array.from(run.root.querySelectorAll(target));
        if (target instanceof Element) return [target];
        if (Array.isArray(target) || target instanceof NodeList) return Array.from(target).flatMap((t) => els(t));
        if (target.g instanceof Element) return [target.g];       // kit objects with a group
        if (target.el instanceof Element) return [target.el];
        return [];
      };
      const k = {
        instant, s: run.s, beat: b, scene: run,
        $: run.s.$, $$: run.s.$$,
        /** Schedule fn after ms (scene clock); instantly in instant mode. */
        at(ms, fn) {
          if (k.instant) { fn(); return; }
          run.timers.push({t: player.clock + (ms * DELAY) / 1000, fn});
        },
        /** Reveal elements (classes sv-in + sv-left / sv-right / sv-down / sv-fade / sv-pop set at build). */
        show(target, o) {
          o = o || {};
          const list = els(target);
          const stagger = o.stagger || 0, delay = o.delay || 0;
          list.forEach((e, i) => {
            const go = () => {
              e.classList.remove("sv-gone", "sv-dim");
              if (!e.classList.contains("sv-in")) return;
              e.classList.add("sv-on");
            };
            const d = delay + i * stagger;
            if (d > 0 && !k.instant) k.at(d, go); else go();
          });
        },
        hide(target, o) {
          o = o || {};
          els(target).forEach((e, i) => {
            const go = () => e.classList.add("sv-gone");
            const d = (o.delay || 0) + i * (o.stagger || 0);
            if (d > 0 && !k.instant) k.at(d, go); else go();
          });
        },
        /** Dim (true) or undim (false) elements. */
        dim(target, on) { els(target).forEach((e) => e.classList.toggle("sv-dim", on !== false)); },
        cls(target, c, on) { els(target).forEach((e) => e.classList.toggle(c, on !== false)); },
        attr(target, attrs) { els(target).forEach((e) => { for (const [a, v] of Object.entries(attrs)) e.setAttribute(a, v); }); },
        text(target, t) { els(target).forEach((e) => { e.textContent = t; }); },
        /** Tween numeric properties of a plain object (scene state); `update` is called each step. */
        tween(obj, props, ms, o) {
          o = o || {};
          const ease = typeof o.ease === "function" ? o.ease : EASE[o.ease || "inOut"];
          const from = {};
          for (const p of Object.keys(props)) from[p] = obj[p] === undefined ? 0 : obj[p];
          const apply = (e) => { for (const p of Object.keys(props)) obj[p] = from[p] + (props[p] - from[p]) * e; if (o.update) o.update(obj); };
          const delay = ((o.delay || 0) * DELAY) / 1000;
          if (k.instant || !ms) { apply(1); if (o.done) o.done(); return; }
          ms *= TWEEN;
          for (const tw of run.tweens) if (tw.obj === obj && Object.keys(props).some((p) => tw.props.includes(p))) tw.killed = true;
          run.tweens.push({obj, props: Object.keys(props), t0: player.clock + delay, dur: ms / 1000, ease, apply, done: o.done});
        },
        /** Move / fade an element: {x, y, scale, rotate, opacity} (CSS individual transforms: an SVG
         *  transform attribute is kept). */
        move(target, props, ms, o) {
          for (const e of els(target)) {
            const st = e.__svT || (e.__svT = {x: 0, y: 0, scale: 1, rotate: 0, opacity: e.style.opacity === "" ? 1 : +e.style.opacity});
            const upd = () => {
              e.style.translate = `${st.x}px ${st.y}px`;
              e.style.scale = String(st.scale);
              e.style.rotate = `${st.rotate}deg`;
              e.style.opacity = String(st.opacity);
            };
            if (props.scale !== undefined || props.rotate !== undefined) { e.style.transformBox = "fill-box"; e.style.transformOrigin = (o && o.origin) || "center"; }
            k.tween(st, props, ms, Object.assign({}, o, {update: upd}));
          }
        },
        /** Draw strokes (elements of class sv-draw are prepared undrawn): from/to are fractions. */
        draw(target, ms, o) {
          o = o || {};
          for (const e of els(target)) {
            run.prepDraw(e.parentNode || run.root);
            if (!e.__svLen) { let len = 0; try { len = e.getTotalLength(); } catch (err) { len = 0; } e.__svLen = Math.ceil(len) + 2; e.style.strokeDasharray = `${e.__svLen} ${e.__svLen}`; }
            const len = e.__svLen;
            const st = {f: o.from !== undefined ? o.from : (1 - (parseFloat(e.style.strokeDashoffset) || 0) / len)};
            e.style.opacity = "";
            e.classList.remove("sv-gone");
            k.tween(st, {f: o.to !== undefined ? o.to : 1}, ms === undefined ? 1200 : ms,
              {ease: o.ease || "inOut", delay: o.delay, update: () => { e.style.strokeDashoffset = String(len * (1 - st.f)); }, done: o.done});
          }
        },
        /** Count a number up in an element's text. */
        count(target, to, o) {
          o = o || {};
          const dec = o.dec || 0;
          const fmtN = (v) => (o.pre || "") + (o.sep ? Number(v.toFixed(dec)).toLocaleString("en-US", {minimumFractionDigits: dec, maximumFractionDigits: dec}) : v.toFixed(dec)) + (o.post || "");
          for (const e of els(target)) {
            const st = {v: o.from !== undefined ? o.from : (parseFloat(e.textContent.replace(/[^\d.\-]/g, "")) || 0)};
            k.tween(st, {v: to}, o.ms || 1200, {ease: o.ease || "out", update: () => { e.textContent = fmtN(st.v); }});
          }
        },
        /** Type the data-text of elements (code lines) one after the other. */
        type(target, o) {
          o = o || {};
          const cps = o.cps || 38;
          let t = o.delay || 0;
          for (const e of els(target)) {
            const full = e.getAttribute("data-text") || "";
            const isCmd = e.getAttribute("data-cmd");
            const set = (n) => {
              const s = full.slice(0, n);
              if (isCmd && !e.classList.contains("c")) {
                const i = s.indexOf(",");
                e.innerHTML = "";
                const head = i < 0 ? s : s.slice(0, i);
                const sp = document.createElement("span"); sp.className = "cmd"; sp.textContent = head; e.appendChild(sp);
                if (i >= 0) e.appendChild(document.createTextNode(s.slice(i)));
              } else e.textContent = s;
            };
            if (k.instant) { set(full.length); e.classList.remove("typing"); e.classList.add("sv-on"); continue; }
            const ms = (full.length / cps) * 1000;
            k.at(t, () => {
              e.classList.add("typing", "sv-on");
              const st = {n: 0};
              k.tween(st, {n: full.length}, ms, {ease: "linear", update: () => set(Math.round(st.n)), done: () => e.classList.remove("typing")});
            });
            t += ms + (o.pause || 180);
          }
        },
        /** A continuous animation for the rest of the scene: fn(t, dt) every frame. */
        loop(fn) { run.loops.push(fn); fn(run.t, 0); },
        /** Attention pulse (scale) on an element. */
        pulse(target, o) {
          if (k.instant) return;
          for (const e of els(target)) {
            e.style.transformBox = "fill-box"; e.style.transformOrigin = "center";
            const st = {p: 0};
            k.tween(st, {p: 1}, (o && o.ms) || 700, {ease: "linear", update: () => { e.style.scale = String(1 + ((o && o.amp) || 0.12) * Math.sin(Math.PI * st.p)); }});
          }
        },
      };
      return k;
    }
  }

  // ================================================================== entry point
  /** Define and mount a video (or register it when SV.catalog is set: the index page). */
  SV.video = function (def) {
    if (SV.catalog) { SV.catalog.push(def); return; }
    const start = () => { SV.player = new Player(def).mount(document.getElementById("sv-root") || document.body); };
    if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
    else start();
  };
  SV.Player = Player;
  SV.estimateMinutes = function (def) {
    const p = Object.create(Player.prototype);
    p.def = def; p.pron = (def.pron || []).concat(PRON); p.rate = 1; p.cps = DEFAULT_CPS;
    p.scenes = def.scenes.map((sc, si) => ({def: sc, si, beats: sc.beats.map((b, bi) => Object.assign({si, bi, def: b}, parseSay(b.say, p.pron)))}));
    p.flat = [];
    p.scenes.forEach((sc) => sc.beats.forEach((b) => p.flat.push(b)));
    p.audio = (SV.AUDIO && SV.AUDIO[def.id]) || null;           // recorded clips: their real durations
    p.flat.forEach((b) => { const t = b.spoken.trim(); b.clip = (p.audio && t && p.audio.beats[SV.hash(t)]) || null; });
    p.natural = p.flat.some((b) => b.clip);
    return p.timeline() / 60000;
  };
  SV.h = h;
  SV.store = store;
})();
