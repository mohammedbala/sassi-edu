/* SASSI-EDU in the browser (GitHub Pages) -- the Web Worker that runs Python.
 *
 * Loads Pyodide 314.0.7 (CPython 3.14 compiled to WebAssembly) from the jsDelivr CDN with NumPy and SciPy,
 * unpacks web/sassi-edu.zip (the sassi package, the examples and the public documentation, written by
 * web/build.py) into the in-memory file system, and answers the requests of the page with
 * sassi.web.bridge.Bridge -- the same API (sassi/ui/api.py) the local server sassi-gui serves.
 *
 * One message at a time: a request is answered at once (the answer is posted), then the module runs that
 * the request queued (RUN<MODULE>, VERIFY) are run by Bridge.run_pending.  While Python is busy the worker
 * cannot answer, so the bridge pushes the new events and the job output to the page ({push}).
 *
 * Python kept between visits: the files Pyodide fetches (the runtime's WebAssembly, the standard library, NumPy,
 * SciPy and their libraries) are read from and stored in the Cache Storage of the site ("sassi-edu-pyodide-<version>",
 * the cache web/sw.js fills too) by keptFetch, so a later visit downloads none of them -- also in a browser where
 * the service worker does not run.  info.cache reports what came from this computer and what was downloaded.
 *
 * Files kept between visits: the workspace and the settings folder are IndexedDB file systems (Emscripten
 * IDBFS): restored before the session starts, saved at most every SAVE_DELAY ms after a request that may
 * write (not GET) or a job, and when the page is hidden or closed ({type: "flush"}).  File > Clear Saved Files
 * (web/boot.js) deletes the two databases.
 *
 * Protocol: see web/boot.js.  A module worker ({type: "module"}); the build fills in the build id. */

const PYODIDE_VERSION = "314.0.7";
const PYODIDE_URL = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;
const BUILD = "__SASSI_BUILD__";
const ZIP_V = "__ZIP_V__";                     // web/build.py: the version (content hash) of web/sassi-edu.zip
const PY_CACHE = `sassi-edu-pyodide-${PYODIDE_VERSION}`;
const ROOT = "/home/pyodide/sassi-edu";        // the unpacked bundle (sassi/, docs/, examples/, README.md)
const WORK = "/home/pyodide/work";             // the workspace (sassi.web.bridge.WORK_DIR)
const SETTINGS = "/home/pyodide/.sassi-edu";   // SASSIini.xml, SASSIdb.xml, SASSIani.xml (bridge.SETTINGS_DIR)
const KEPT = [WORK, SETTINGS];                 // IndexedDB file systems: kept between visits
const SAVE_DELAY = 1000;

let bridge = null;
let pyFS = null;                               // Pyodide's file system once the kept folders are mounted
const early = [];                              // requests that arrived before Python was ready
let timer = null;

function status(text, fraction) { self.postMessage({type: "status", text, fraction}); }

// ------------------------------------------------------------------ Python kept between visits (Cache Storage)
const fetched = {kept: 0, keptBytes: 0, downloaded: 0, downloadedBytes: 0, storable: typeof caches !== "undefined"};
const networkFetch = self.fetch.bind(self);
const sizeOf = (res) => Number(res.headers.get("content-length")) || 0;
/** fetch for the Pyodide files: the stored copy when there is one, else the network (and the copy is stored).
 *  Other requests, and a browser without Cache Storage (a private window may refuse it), use the network. */
async function keptFetch(input, init) {
  const url = typeof input === "string" ? input : (input && input.url) || String(input);
  const method = (init && init.method) || (input && input.method) || "GET";
  if (!fetched.storable || method !== "GET" || !url.startsWith(PYODIDE_URL)) return networkFetch(input, init);
  let cache = null;
  try {
    cache = await caches.open(PY_CACHE);
    const hit = await cache.match(url);
    if (hit) { fetched.kept += 1; fetched.keptBytes += sizeOf(hit); return hit; }
  } catch (err) {
    fetched.storable = false;
  }
  const res = await networkFetch(input, init);
  fetched.downloaded += 1;
  fetched.downloadedBytes += sizeOf(res);
  if (cache && res.ok && res.status === 200) cache.put(url, res.clone()).catch(() => {});
  return res;
}
self.fetch = keptFetch;
/** Number of Pyodide files an earlier visit stored (0 without Cache Storage). */
async function storedCount() {
  try { return fetched.storable && (await caches.has(PY_CACHE)) ? (await (await caches.open(PY_CACHE)).keys()).length : 0; }
  catch (err) { return 0; }
}

self.onmessage = (ev) => {
  const m = ev.data || {};
  if (m.type === "flush") { save(); return; }
  if (m.id === undefined) return;
  if (!bridge) { early.push(m); return; }
  answer(m);
};

/** Answer one request, then let the jobs it queued run (after any request already waiting). */
function answer(m) {
  let text;
  try {
    text = bridge.handle(m.method, m.path, m.query || "", m.body || "");
  } catch (err) {
    if (err && err.name !== "PythonError") { fatal(err); return; }
    text = JSON.stringify({status: 500, payload: {error: `Python error: ${err.message || err}`}});
  }
  self.postMessage({id: m.id, answer: text});
  if (m.method !== "GET") saveSoon();          // a GET reads (the event poll runs every 250 ms)
  schedule();
}

function schedule() {
  if (timer === null && bridge && bridge.pending) timer = setTimeout(runJobs, 0);
}

/** One queued job, to its end (Bridge.run_pending pushes its output on the way); then the next one. */
function runJobs() {
  timer = null;
  let more = false;
  try {
    more = bridge.run_pending();
  } catch (err) {
    if (err && err.name !== "PythonError") { fatal(err); return; }
    console.error(err);
  }
  if (more) schedule();
  saveSoon();
}

// ------------------------------------------------------------------ files kept between visits (IndexedDB)
let saveTimer = null, saving = false, saveAgain = false;
function saveSoon() {
  if (pyFS && saveTimer === null) saveTimer = setTimeout(save, SAVE_DELAY);
}
function save() {
  if (saveTimer !== null) { clearTimeout(saveTimer); saveTimer = null; }
  if (!pyFS) return;
  if (saving) { saveAgain = true; return; }
  saving = true;
  pyFS.syncfs(false, (err) => {
    saving = false;
    if (err) console.warn("SASSI-EDU: the workspace could not be saved in this browser:", err);
    if (saveAgain) { saveAgain = false; saveSoon(); }
  });
}
/** Mount the kept folders on IndexedDB and restore what an earlier visit saved; false without IndexedDB
 *  (a private window may refuse it): the files then live in this tab only. */
async function restoreKept(py) {
  const FS = py.FS;
  if (typeof indexedDB === "undefined" || !FS || !FS.filesystems || !FS.filesystems.IDBFS) return false;
  try {
    for (const d of KEPT) {
      FS.mkdirTree(d);
      FS.mount(FS.filesystems.IDBFS, {}, d);
    }
    await new Promise((res, rej) => FS.syncfs(true, (err) => (err ? rej(err) : res())));
    pyFS = FS;
    return true;
  } catch (err) {
    console.warn("SASSI-EDU: files are not kept between visits:", err);
    return false;
  }
}

function fatal(err, what) {
  self.postMessage({type: "fatal", error: `${what || "The Python engine stopped"}: ${(err && err.message) || err}`});
  bridge = null;
}

async function start() {
  const t0 = performance.now();
  const stored = await storedCount();
  const where = stored >= 4 ? " from this computer" : "";
  status(stored >= 4 ? `Starting the Python runtime (Pyodide ${PYODIDE_VERSION}) from this computer ...`
    : `Downloading the Python runtime (Pyodide ${PYODIDE_VERSION}; first visit) ...`, 0.05);
  let loadPyodide;
  try {
    ({loadPyodide} = await import(`${PYODIDE_URL}pyodide.mjs`));
  } catch (err) {
    throw new Error(`Pyodide could not be loaded from ${PYODIDE_URL} (${err.message || err})`);
  }
  const py = await loadPyodide({indexURL: PYODIDE_URL});
  status(stored >= 4 ? `Loading NumPy and SciPy${where} ...` : "Downloading NumPy and SciPy (first visit) ...", 0.3);
  await py.loadPackage(["numpy", "scipy"], {messageCallback: () => {}});
  status("Loading SASSI-EDU ...", 0.85);
  const url = new URL(`sassi-edu.zip?v=${ZIP_V}`, self.location.href);
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url.pathname}: HTTP ${r.status}`);
  py.unpackArchive(await r.arrayBuffer(), "zip", {extractDir: ROOT});
  status("Restoring your files ...", 0.9);
  const kept = await restoreKept(py);
  status("Starting the SASSI-EDU session ...", 0.93);
  py.runPython(`
import os, sys
sys.path.insert(0, ${JSON.stringify(ROOT)})
os.makedirs(${JSON.stringify(WORK)}, exist_ok=True)
os.chdir(${JSON.stringify(WORK)})
`);
  const mod = py.pyimport("sassi.web.bridge");
  bridge = mod.Bridge((text) => self.postMessage({push: text}));
  const info = JSON.parse(py.runPython(`
import json, os, sys, numpy, scipy, sassi
n = size = 0
for base, _dirs, names in os.walk(${JSON.stringify(WORK)}):
    for name in names:
        try:
            size += os.path.getsize(os.path.join(base, name))
            n += 1
        except OSError:
            pass
json.dumps({"python": sys.version.split()[0], "numpy": numpy.__version__, "scipy": scipy.__version__,
            "sassi": sassi.__version__, "kept": {"files": n, "bytes": size}})
`));
  info.kept.on = kept;
  info.cache = Object.assign({storedBefore: stored}, fetched);
  info.pyodide = PYODIDE_VERSION;
  info.build = BUILD;
  info.seconds = (performance.now() - t0) / 1000;
  // what the worker downloaded (bytes over the network: 0 when the browser cache had the file)
  info.downloads = performance.getEntriesByType("resource").map((e) => ({
    name: e.name, transfer: e.transferSize, encoded: e.encodedBodySize, decoded: e.decodedBodySize,
    ms: Math.round(e.duration)}));
  self.postMessage({type: "ready", info});
  while (early.length && bridge) answer(early.shift());
}

start().catch((err) => fatal(err, "Loading failed"));
