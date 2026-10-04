/* SASSI-EDU in the browser -- Node test of the built bundle in Pyodide (no browser needed).
 *
 *   cd web && npm install            # once: the pinned pyodide package (314.0.7, as the site loads)
 *   python web/build.py              # the site, with web/site/web/sassi-edu.zip
 *   node web/test_pyodide.mjs        # [--site DIR] [--quick] [--lessons 01-why-ssi,02-free-field]
 *
 * It loads web/site/web/sassi-edu.zip into Pyodide exactly as web/worker.js does and drives
 * sassi.web.bridge.Bridge the way the page does:
 *
 *   1. start: GET /api/state, /api/events (must not block), /api/lessons, Help with a figure, the explainer;
 *   2. lesson 01 as a learner runs it: open (fresh workspace), every step's commands one per request, the
 *      RUN<MODULE> jobs by run_pending (their listing pushed while they run, then GET /api/jobs/<id>), the
 *      step actions (plots: SPECPLOT / MODELPLOT data through /api/plot/<id>; listings through /api/file);
 *   3. File > Upload to Workspace / Download (POST and GET /api/file);
 *   4. unless --quick: every lesson headless (sassi.ui.lessons.run_lesson_headless), with its time.
 *
 * Prints one "RESULT {...}" line (JSON) at the end; exit status 1 when a check failed. */
import { loadPyodide } from "pyodide";
import { readFileSync, existsSync } from "node:fs";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const HERE = dirname(fileURLToPath(import.meta.url));
const args = process.argv.slice(2);
const opt = (name, dflt) => { const i = args.indexOf(name); return i >= 0 ? args[i + 1] : dflt; };
const SITE = resolve(opt("--site", join(HERE, "site")));
const QUICK = args.includes("--quick");
const ONLY = (opt("--lessons", "") || "").split(",").filter(Boolean);
const ROOT = "/home/pyodide/sassi-edu";
const WORK = "/home/pyodide/work";
const LESSON = "01-why-ssi";

const failures = [];
const timings = {};
function check(cond, what) { if (!cond) { failures.push(what); console.log(`FAIL ${what}`); } }
const now = () => performance.now() / 1000;
const zip = join(SITE, "web", "sassi-edu.zip");
if (!existsSync(zip)) { console.error(`${zip} not found: run python web/build.py first`); process.exit(2); }

// ------------------------------------------------------------------ Pyodide as web/worker.js sets it up
let t = now();
const py = await loadPyodide();
await py.loadPackage(["numpy", "scipy"], {messageCallback: () => {}});
py.unpackArchive(new Uint8Array(readFileSync(zip)), "zip", {extractDir: ROOT});   // as fetch().arrayBuffer() gives it
py.runPython(`
import os, sys
sys.path.insert(0, ${JSON.stringify(ROOT)})
os.makedirs(${JSON.stringify(WORK)}, exist_ok=True)
os.chdir(${JSON.stringify(WORK)})
`);
const pushes = [];
const bridge = py.pyimport("sassi.web.bridge").Bridge((text) => pushes.push(JSON.parse(text)));
timings.load = +(now() - t).toFixed(1);
console.log(`Pyodide ${py.version} + numpy + scipy + bundle: ${timings.load} s`);

/** One request as window.SASSI_TRANSPORT sends it; returns {status, payload}. */
function api(method, path, body) {
  const i = path.indexOf("?");
  const text = bridge.handle(method, i < 0 ? path : path.slice(0, i), i < 0 ? "" : path.slice(i + 1),
                             body === undefined ? "" : JSON.stringify(body));
  return JSON.parse(text);
}
const get = (p) => api("GET", p);
const post = (p, b) => api("POST", p, b || {});

/** The page's job view: listing lines from pushes and answers, de-duplicated by index (static/app.js applyJob). */
const shown = {};
function applyJob(d) {
  const j = shown[d.id] || (shown[d.id] = {next: 0, lines: []});
  const first = d.next - d.messages.length;
  for (const m of d.messages.slice(Math.max(0, j.next - first))) j.lines.push(m.text);
  j.next = Math.max(j.next, d.next);
}
/** What web/worker.js does after posting an answer: run the queued jobs. */
function runJobs() {
  let n = 0;
  while (bridge.pending) { bridge.run_pending(); n++; }
  return n;
}

// ------------------------------------------------------------------ 1. start-up requests
t = now();
let r = get("/api/state");
check(r.status === 200 && r.payload.web === true && r.payload.cwd === WORK, "GET /api/state (web mode, workspace)");
let t1 = now();
r = get("/api/events?since=0&timeout=20");
check(r.status === 200 && now() - t1 < 1.0, "GET /api/events does not wait (timeout ignored)");
r = get("/api/lessons");
const lessonIds = r.payload.parts.flatMap((p) => p.lessons.map((l) => l.id));
check(lessonIds.length === 10 && lessonIds[0] === LESSON, `10 lessons (got ${lessonIds.length})`);
r = get("/api/help");
check(r.status === 200 && r.payload.docs.some((d) => d.id === "docs/user/GUI.md"), "Help lists the documents");
check(!r.payload.docs.some((d) => d.id.startsWith("docs/spec/") || d.id.startsWith("docs/internal/")), "no internal documents");
r = get("/api/help/doc?name=" + encodeURIComponent("docs/verification/VERIFICATION_MANUAL.md"));
const img = (r.payload.html || "").match(/<img src="(docs\/[^"]+)"/);
check(img && existsSync(join(SITE, decodeURIComponent(img[1]))), "Help figure with a relative docs/ URL that the site has");
r = get(`/api/lessons/${LESSON}`);
check(JSON.stringify(r.payload).includes("data-tex"), "lesson text with formulas for KaTeX");
r = get("/api/explain?line=" + encodeURIComponent("RUNSITE"));
check(r.status === 200 && r.payload.name === "RUNSITE", "the command explainer");
r = get("/api/about");
check(r.payload.web === true && /Pyodide/.test(r.payload.build), "Help > About in web mode");
timings.startup = +(now() - t).toFixed(1);

// ------------------------------------------------------------------ 2. lesson 01 through the bridge
t = now();
r = post(`/api/lessons/${LESSON}/open`, {});
check(r.status === 200 && r.payload.commands && r.payload.workspace.startsWith(WORK), "open lesson 01");
for (const ln of r.payload.commands) check(post("/api/command", {lines: [ln]}).payload.ok, `setup ${ln}`);
const lesson = get(`/api/lessons/${LESSON}`).payload;
const jobs = [];
let plotsDrawn = 0;
let animations = 0;
for (const step of lesson.steps) {
  const ts = now();
  for (const ln of step.commands) {
    if (!ln.trim()) continue;
    const res = post("/api/command", {lines: [ln]}).payload;
    if (res.job) {
      check(res.job.state === "starting", `${ln}: the job is queued (starting) when the request answers`);
      const before = pushes.length;
      runJobs();
      for (const p of pushes.slice(before)) if (p.job) applyJob(p.job);
      const d = get(`/api/jobs/${res.job.id}?since=0`).payload;
      applyJob(d);
      check(d.state === "done" && d.ok, `${ln}: job ${d.state}`);
      check(shown[d.id].lines.length === d.next, `${ln}: listing shown once (${shown[d.id].lines.length} of ${d.next})`);
      jobs.push({line: ln, messages: d.next, pushes: pushes.length - before, s: +d.elapsed.toFixed(2)});
    } else {
      check(res.ok, `step ${step.index}: ${ln}`);
    }
  }
  post(`/api/lessons/${LESSON}/progress`, {step: step.index});
  for (const a of step.actions) {
    const ar = post(`/api/lessons/${LESSON}/action`, {verb: a.verb, args: a.args});
    if (ar.status !== 200) { check(false, `action ${a.verb} ${a.args}: ${ar.payload.error}`); continue; }
    const act = ar.payload;
    if (act.kind === "commands") {
      const res = post("/api/command", {lines: act.lines}).payload;
      check(res.ok, `action ${a.verb} ${a.args}`);
      const st = get("/api/plots").payload;
      const pd = get(`/api/plot/${st.active}`);
      check(pd.status === 200 && pd.payload && Object.keys(pd.payload).length > 0, `plot data of ${a.verb}`);
      if (pd.status === 200 && pd.payload.family === "anim") {      // an animation: the player fetches every frame
        const fr = pd.payload.animation.frames;
        const f2 = get(`/api/plot/${st.active}?frame=${fr[Math.min(1, fr.length - 1)]}`);
        check(f2.status === 200 && f2.payload.label && (f2.payload.xyz || f2.payload.value || f2.payload.vectors),
              `animation frame of ${a.verb} ${a.args}`);
        animations++;
      }
      plotsDrawn++;
    } else if (act.kind === "file") {
      const f = get("/api/file?name=" + encodeURIComponent(act.path));
      check(f.status === 200 && f.payload.text.length > 0, `action ${a.verb} ${a.args}: file`);
    }
  }
  console.log(`  step ${step.index} ${step.title}: ${(now() - ts).toFixed(1)} s`);
}
const hist = get("/api/state").payload.history;
check(["RUNSITE", "RUNPOINT", "RUNHOUSE", "RUNANALYS", "RUNMOTION"].every((x) => hist.includes(x)), "replay history holds the module runs (L17)");
check(jobs.some((j) => j.pushes > 0), "job output pushed while a job ran");
const evPushes = pushes.filter((p) => p.events && p.events.length).length;
check(evPushes > 0, "events pushed while Python was busy");
timings.lesson01 = +(now() - t).toFixed(1);

// ------------------------------------------------------------------ 3. upload / download
const text = "* uploaded by the test\nN,1,0,0,0\n";
r = post("/api/file", {name: `${WORK}/uploaded.pre`, text});
check(r.status === 200 && r.payload.ok, "File > Upload to Workspace (POST /api/file)");
r = get("/api/file?name=" + encodeURIComponent(`${WORK}/uploaded.pre`));
check(r.status === 200 && r.payload.text === text, "File > Download (GET /api/file)");

// ------------------------------------------------------------------ 4. every lesson headless
const headless = [];
if (!QUICK) {
  const ids = ONLY.length ? ONLY : lessonIds;
  for (const id of ids) {
    const out = JSON.parse(py.runPython(`
import json, tempfile, time
from sassi.ui import lessons as L
les = L.lesson_by_id(${JSON.stringify(id)})
t = time.time()
rep = L.run_lesson_headless(les, tempfile.mkdtemp(dir=${JSON.stringify(WORK)}))
json.dumps({"id": les.id, "ok": rep.ok, "s": round(time.time() - t, 1),
            "errors": [str(e)[:200] for e in rep.errors[:3]], "missing": [str(m) for m in rep.missing[:3]]})
`));
    headless.push(out);
    check(out.ok, `lesson ${id} headless: ${JSON.stringify(out.errors.concat(out.missing))}`);
    console.log(`  headless ${out.id}: ${out.ok ? "ok" : "FAILED"} ${out.s} s`);
  }
}

check(animations > 0, "lesson 01 animations (animate actions) drawn");
const result = {ok: failures.length === 0, failures, timings, jobs, plots: plotsDrawn, animations, pushes: pushes.length,
                event_pushes: evPushes, headless};
console.log("RESULT " + JSON.stringify(result));
process.exit(failures.length ? 1 : 0);
