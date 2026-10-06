// Node test of the explainer videos (sassi/ui/static/videos): every beat of every video in headless Chrome.
//
//   node web/test_videos.mjs                         # all videos: every beat applied, errors reported
//   node web/test_videos.mjs 01 03                   # some videos
//   node web/test_videos.mjs 01 --shots DIR          # also a JPEG of the end of every scene (DIR/01_s03.jpg)
//   node web/test_videos.mjs 01 --shots DIR --all    # ... of the end of every beat (DIR/01_s03_b02.jpg)
//   node web/test_videos.mjs 01 --play 20            # also play 20 s from the start (voice off) and check the clock
//   node web/test_videos.mjs 01 --audio 15           # also play 15 s with the recorded narration (muted output)
//
// Uses the Chrome DevTools protocol over Node's built-in WebSocket (Node 22+), no packages.  The pages are
// opened from file:// in ?shot mode (the 1920 x 1080 stage fills the window).  Exit code 1 on any error.
import {spawn} from "node:child_process";
import {mkdtempSync, rmSync, writeFileSync, mkdirSync, readdirSync, existsSync} from "node:fs";
import {tmpdir} from "node:os";
import {join, resolve, dirname} from "node:path";
import {fileURLToPath, pathToFileURL} from "node:url";

const HERE = dirname(fileURLToPath(import.meta.url));
const VIDEOS = resolve(HERE, "../sassi/ui/static/videos");
const CHROME = process.env.CHROME || "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const args = process.argv.slice(2);
const opt = (name) => { const i = args.indexOf(name); return i >= 0 ? args[i + 1] : null; };
const shots = opt("--shots");
const all = args.includes("--all");
const play = opt("--play");
const audio = opt("--audio");
const scale = +(opt("--scale") || 0.5);
const ids = args.filter((a, i) => /^\d\d$/.test(a) && !["--shots", "--play", "--scale", "--audio"].includes(args[i - 1]));
const list = ids.length ? ids : readdirSync(VIDEOS).filter((f) => /^\d\d\.html$/.test(f)).map((f) => f.slice(0, 2)).sort();

const prof = mkdtempSync(join(tmpdir(), "sv-chrome-"));
const port = 9300 + Math.floor(Math.random() * 500);
const chrome = spawn(CHROME, ["--headless=new", `--remote-debugging-port=${port}`, `--user-data-dir=${prof}`, "--hide-scrollbars",
  "--mute-audio", "--autoplay-policy=no-user-gesture-required", "--no-first-run", "--disable-extensions", "--allow-file-access-from-files", "--window-size=1920,1080", "about:blank"], {stdio: "ignore"});
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
let failed = 0;

async function target() {
  for (let i = 0; i < 80; i++) {
    try {
      const r = await fetch(`http://127.0.0.1:${port}/json/list`);
      const t = (await r.json()).find((x) => x.type === "page");
      if (t) return t.webSocketDebuggerUrl;
    } catch (e) { /* not up yet */ }
    await sleep(150);
  }
  throw new Error("Chrome did not start (set CHROME=/path/to/chrome)");
}

function cdp(url) {
  const ws = new WebSocket(url);
  let id = 0;
  const pending = new Map(), handlers = [];
  ws.onmessage = (ev) => {
    const m = JSON.parse(ev.data);
    if (m.id && pending.has(m.id)) { const {res, rej} = pending.get(m.id); pending.delete(m.id); m.error ? rej(new Error(m.error.message)) : res(m.result); }
    else if (m.method) handlers.forEach((h) => h(m));
  };
  const ready = new Promise((r) => { ws.onopen = r; });
  return {
    ready, on: (h) => handlers.push(h), close: () => ws.close(),
    send: (method, params) => new Promise((res, rej) => { const i = ++id; pending.set(i, {res, rej}); ws.send(JSON.stringify({id: i, method, params: params || {}})); }),
  };
}

try {
  const c = cdp(await target());
  await c.ready;
  const errors = [];
  c.on((m) => {
    if (m.method === "Runtime.exceptionThrown") errors.push(m.params.exceptionDetails.exception?.description || m.params.exceptionDetails.text);
    if (m.method === "Runtime.consoleAPICalled" && (m.params.type === "error" || m.params.type === "warning"))
      errors.push("console." + m.params.type + ": " + m.params.args.map((a) => a.value ?? a.description).join(" "));
    if (m.method === "Log.entryAdded" && m.params.entry.level === "error") errors.push("log: " + m.params.entry.text + " " + (m.params.entry.url || ""));
  });
  await c.send("Runtime.enable");
  await c.send("Log.enable");
  await c.send("Page.enable");
  await c.send("Emulation.setDeviceMetricsOverride", {width: 1920, height: 1080, deviceScaleFactor: 1, mobile: false});
  const ev = async (expr) => {
    const r = await c.send("Runtime.evaluate", {expression: expr, awaitPromise: true, returnByValue: true});
    if (r.exceptionDetails) throw new Error(r.exceptionDetails.exception?.description || r.exceptionDetails.text);
    return r.result.value;
  };
  if (shots) mkdirSync(shots, {recursive: true});
  for (const v of list) {
    const file = join(VIDEOS, `${v}.html`);
    if (!existsSync(file)) { console.log(`${v}: no ${v}.html`); failed++; continue; }
    errors.length = 0;
    await c.send("Page.navigate", {url: pathToFileURL(file).href + "?shot"});
    let ok = false;
    for (let i = 0; i < 60 && !ok; i++) { await sleep(100); try { ok = await ev("!!(window.SV && SV.player && SV.player.scenes)"); } catch (e) { /* loading */ } }
    if (!ok) { console.log(`${v}: the player did not start ${errors.join(" | ")}`); failed++; continue; }
    await ev("document.fonts.ready.then(() => true)");
    const info = await ev("({scenes: SV.player.scenes.map((s) => ({id: s.def.id, n: s.beats.length})), total: SV.player.total})");
    const beatErrors = [];
    for (let si = 0; si < info.scenes.length; si++) {
      const n = info.scenes[si].n;
      for (let bi = 0; bi < n; bi++) {
        try { await ev(`SV.player.snap(${si}, ${bi}), true`); }
        catch (e) { beatErrors.push(`scene ${si} (${info.scenes[si].id}) beat ${bi}: ${e.message.split("\n")[0]}`); continue; }
        if (shots && (all || bi === n - 1)) {
          await sleep(60);
          const r = await c.send("Page.captureScreenshot", {format: "jpeg", quality: 72, clip: {x: 0, y: 0, width: 1920, height: 1080, scale}});
          writeFileSync(join(shots, all ? `${v}_s${String(si).padStart(2, "0")}_b${String(bi).padStart(2, "0")}.jpg` : `${v}_s${String(si).padStart(2, "0")}.jpg`), Buffer.from(r.data, "base64"));
        }
      }
    }
    // a short real-time play (voice off: the clock runs on estimated durations)
    if (play) {
      await ev("SV.player.setVoiceOn(false), SV.player.seek(0, 0, true), true");
      await sleep(+play * 1000);
      const st = await ev("({i: SV.player.beat.b.i, clock: SV.player.clock})");
      if (st.i < 1) beatErrors.push(`play: still on the first beat after ${play} s`);
      console.log(`${v}: played ${play} s -> beat ${st.i}, clock ${st.clock.toFixed(1)} s`);
      await ev("SV.player.pause(), true");
    }
    // a real-time play with the recorded narration: the clips load, play and drive the beats
    if (audio) {
      const has = await ev("SV.player.hasNatural");
      if (!has) beatErrors.push("audio: no recorded narration (aNN.js)");
      else {
        await ev("SV.player.setVoiceOn(true), SV.player.pickVoice('__natural__'), SV.player.timeline(), SV.player.seek(0, 0, true), true");
        const t0 = Date.now();
        let maxT = 0, beats = new Set();
        while (Date.now() - t0 < +audio * 1000) {
          await sleep(500);
          const st = await ev("({i: SV.player.beat.b.i, t: SV.player.beat.audio ? SV.player.beat.audio.currentTime : -1, failed: SV.player.flat.filter((b) => b.clipFailed).length})");
          beats.add(st.i); maxT = Math.max(maxT, st.t);
          if (st.failed) { beatErrors.push(`audio: ${st.failed} clip(s) failed to play`); break; }
        }
        if (beats.size < 2 || maxT <= 0) beatErrors.push(`audio: the clips did not drive the beats (beats ${[...beats]}, max t ${maxT})`);
        console.log(`${v}: ${audio} s of recorded narration -> beats ${[...beats].join(",")}`);
        await ev("SV.player.pause(), true");
      }
    }
    const errs = beatErrors.concat(errors.splice(0));
    const beats = info.scenes.reduce((a, s) => a + s.n, 0);
    console.log(`${v}: ${info.scenes.length} scenes, ${beats} beats, ~${(info.total / 60000).toFixed(1)} min${errs.length ? "  FAILED" : "  ok"}`);
    for (const e of errs) console.log("   " + e);
    if (errs.length) failed++;
  }
  c.close();
} catch (e) {
  console.error(e.message);
  failed++;
} finally {
  chrome.kill();
  await sleep(300);
  rmSync(prof, {recursive: true, force: true});
}
process.exit(failed ? 1 : 0);
