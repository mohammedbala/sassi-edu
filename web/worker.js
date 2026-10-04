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
 * Protocol: see web/boot.js.  A module worker ({type: "module"}); the build fills in the build id. */

const PYODIDE_VERSION = "314.0.7";
const PYODIDE_URL = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;
const BUILD = "__SASSI_BUILD__";
const ROOT = "/home/pyodide/sassi-edu";        // the unpacked bundle (sassi/, docs/, examples/, README.md)
const WORK = "/home/pyodide/work";             // the workspace (sassi.web.bridge.WORK_DIR)

let bridge = null;
const early = [];                              // requests that arrived before Python was ready
let timer = null;

function status(text, fraction) { self.postMessage({type: "status", text, fraction}); }

self.onmessage = (ev) => {
  const m = ev.data || {};
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
}

function fatal(err, what) {
  self.postMessage({type: "fatal", error: `${what || "The Python engine stopped"}: ${(err && err.message) || err}`});
  bridge = null;
}

async function start() {
  const t0 = performance.now();
  status(`Downloading the Python runtime (Pyodide ${PYODIDE_VERSION}) ...`, 0.05);
  let loadPyodide;
  try {
    ({loadPyodide} = await import(`${PYODIDE_URL}pyodide.mjs`));
  } catch (err) {
    throw new Error(`Pyodide could not be loaded from ${PYODIDE_URL} (${err.message || err})`);
  }
  const py = await loadPyodide({indexURL: PYODIDE_URL});
  status("Downloading NumPy and SciPy ...", 0.3);
  await py.loadPackage(["numpy", "scipy"], {messageCallback: () => {}});
  status("Downloading SASSI-EDU ...", 0.85);
  const url = new URL(`sassi-edu.zip?v=${BUILD}`, self.location.href);
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url.pathname}: HTTP ${r.status}`);
  py.unpackArchive(await r.arrayBuffer(), "zip", {extractDir: ROOT});
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
import json, sys, numpy, scipy, sassi
json.dumps({"python": sys.version.split()[0], "numpy": numpy.__version__, "scipy": scipy.__version__,
            "sassi": sassi.__version__})
`));
  info.pyodide = PYODIDE_VERSION;
  info.seconds = (performance.now() - t0) / 1000;
  // what the worker downloaded (bytes over the network: 0 when the browser cache had the file)
  info.downloads = performance.getEntriesByType("resource").map((e) => ({
    name: e.name, transfer: e.transferSize, encoded: e.encodedBodySize, decoded: e.decodedBodySize,
    ms: Math.round(e.duration)}));
  self.postMessage({type: "ready", info});
  while (early.length && bridge) answer(early.shift());
}

start().catch((err) => fatal(err, "Loading failed"));
