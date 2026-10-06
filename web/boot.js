/* SASSI-EDU in the browser (GitHub Pages) -- main-thread boot, loaded before the GUI scripts.
 *
 *  * starts the Web Worker (web/worker.js) that runs the sassi package with Pyodide (Python compiled to
 *    WebAssembly) -- there is no server: everything runs on the visitor's computer;
 *  * shows the loading overlay (what is being downloaded, elapsed time, the disclaimer) until Python is
 *    ready, and a clear error when WebAssembly, module workers or the Pyodide CDN are not available;
 *  * sets window.SASSI_TRANSPORT(method, path, body) -> Promise<{status, payload}>, which static/app.js uses
 *    instead of fetch(): the request goes to the worker by postMessage and the answer comes back the same way;
 *  * hands the pushes of the worker (events and job output sent while Python is busy) to SASSI.push;
 *  * asks the worker to save the workspace in this browser when the page is hidden or closed, and offers
 *    SASSI_WEB.clearSaved() (File > Clear Saved Files: delete the saved workspace and start afresh).
 *
 * Protocol (web/worker.js): page -> worker {id, method, path, query, body (JSON text)}; worker -> page
 * {id, answer (JSON text {status, payload})}, {push (JSON text)}, {type: "status" | "ready" | "fatal", ...}.
 * Vanilla JavaScript, no build step (web/build.py only fills in the build id). */
"use strict";

(function () {
  const PYODIDE_VERSION = "314.0.7";
  const BUILD = "__SASSI_BUILD__";            // web/build.py: a hash of the bundle (cache busting)
  const DISCLAIMER = "SASSI-EDU is an educational re-implementation of the ACS SASSI soil-structure-interaction " +
    "methodology. It is not affiliated with or endorsed by ACS SASSI or its vendor, and it is not qualified for " +
    "design or licensing work.";
  const t0 = Date.now();

  // ------------------------------------------------------------------ loading overlay
  const css = `
#web-overlay { position: fixed; inset: 0; z-index: 1000; display: flex; align-items: center; justify-content: center;
  background: rgba(38, 49, 63, .55); font: 13px/1.45 var(--ui, sans-serif); color: var(--text, #1c222a); }
#web-overlay.gone { opacity: 0; pointer-events: none; transition: opacity .35s; }
#web-overlay .card { width: min(560px, calc(100vw - 32px)); background: var(--panel, #fff); border: 1px solid var(--border, #c5cbd3);
  border-radius: 6px; box-shadow: 0 12px 36px rgba(0,0,0,.25); padding: 20px 24px; }
#web-overlay h1 { font-size: 18px; margin: 0 0 4px; }
#web-overlay .sub { color: var(--muted, #5d6774); margin: 0 0 14px; }
#web-overlay .stage { font-weight: 600; margin: 0 0 6px; }
#web-overlay .bar { height: 8px; background: var(--border-2, #dde1e6); border-radius: 4px; overflow: hidden; }
#web-overlay .bar > div { height: 100%; width: 0; background: var(--accent, #1f5fbf); transition: width .4s; }
#web-overlay .time { color: var(--muted, #5d6774); font-size: 12px; margin: 4px 0 12px; }
#web-overlay .note { color: var(--muted, #5d6774); font-size: 12px; margin: 8px 0 0; }
#web-overlay .disclaimer { font-size: 12px; margin: 12px 0 0; padding: 8px 10px; background: var(--panel-2, #f6f7f9);
  border-left: 3px solid var(--warn, #a86b00); }
#web-overlay .error { display: none; margin: 10px 0 0; padding: 8px 10px; color: var(--danger, #c62828);
  background: #fdecea; border: 1px solid #f3b8b3; border-radius: 4px; white-space: pre-wrap; overflow-wrap: anywhere; }
#web-overlay.failed .error { display: block; }
#web-overlay.failed .bar > div { background: var(--danger, #c62828); }
`;
  const style = document.createElement("style");
  style.textContent = css;
  document.head.appendChild(style);

  function node(tag, cls, text) {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text) e.textContent = text;
    return e;
  }
  const overlay = node("div");
  overlay.id = "web-overlay";
  overlay.setAttribute("role", "status");
  overlay.setAttribute("aria-live", "polite");
  const card = node("div", "card");
  const stage = node("p", "stage", "Starting ...");
  const bar = node("div", "bar");
  const fill = node("div");
  bar.appendChild(fill);
  const time = node("p", "time", "");
  const error = node("div", "error");
  const note = node("p", "note", "");
  const FIRST_VISIT = "This first visit downloads about 27 MB (the Python runtime, NumPy, SciPy, Plotly and SASSI-EDU) " +
    "and keeps it on this computer: the next visits start without downloading anything.";
  const FROM_CACHE = "Loading from this computer: an earlier visit stored the Python runtime, NumPy, SciPy and " +
    "SASSI-EDU, so nothing is downloaded. Starting Python still takes a few seconds.";
  const NO_SW = "The first visit downloads about 27 MB (the Python runtime, NumPy, SciPy, Plotly and SASSI-EDU); the " +
    "browser caches it for the next visits.";
  card.append(node("h1", "", "SASSI-EDU — soil-structure interaction in your browser"),
    node("p", "sub", "The guided course, the examples, the plots and Help of the SASSI-EDU GUI, with Python running in this tab."),
    stage, bar, time, note,
    node("p", "note", "Everything runs on your computer: nothing is sent to a server. Your files are kept in this " +
      "browser between visits (File > Clear Saved Files removes them); File > Download saves one on your computer."),
    node("p", "disclaimer", "Disclaimer. " + DISCLAIMER), error);
  overlay.appendChild(card);
  document.body.appendChild(overlay);

  let fraction = 0;
  const tick = setInterval(() => { time.textContent = `${((Date.now() - t0) / 1000).toFixed(0)} s`; }, 500);
  function setStage(text, frac) {
    stage.textContent = text;
    if (typeof frac === "number") fraction = Math.max(fraction, Math.min(1, frac));
    fill.style.width = `${(100 * fraction).toFixed(0)}%`;
  }
  function fail(message) {
    clearInterval(tick);
    overlay.classList.add("failed");
    overlay.classList.remove("gone");
    stage.textContent = "SASSI-EDU could not start";
    error.textContent = message + "\n\nSASSI-EDU in the browser needs a current desktop browser (Chrome, Edge, Firefox or " +
      "Safari) with WebAssembly and module Web Workers, and access to the Pyodide files on cdn.jsdelivr.net. " +
      "Reload the page to try again.";
  }
  function done() {
    clearInterval(tick);
    setStage(`Ready in ${((Date.now() - t0) / 1000).toFixed(1)} s`, 1);
    window.SASSI_WEB.readyTime = (Date.now() - t0) / 1000;
    setTimeout(() => { overlay.classList.add("gone"); setTimeout(() => overlay.remove(), 400); }, 300);
  }

  // ------------------------------------------------------------------ the worker and the transport
  let readyResolve, readyReject;
  const ready = new Promise((res, rej) => { readyResolve = res; readyReject = rej; });
  ready.catch(() => {});                       // a failed start is shown by the overlay
  const pending = new Map();
  let seq = 0;
  let worker = null;
  window.SASSI_WEB = {version: PYODIDE_VERSION, build: BUILD, ready, readyTime: null, info: null};
  /** The IndexedDB databases of the kept folders (web/worker.js KEPT: Emscripten names them by mount point). */
  const KEPT_DBS = ["/home/pyodide/work", "/home/pyodide/.sassi-edu"];
  /** How much this site stores in the browser (the app files and the saved workspace), when the browser says. */
  window.SASSI_WEB.storageUse = async function () {
    try { return navigator.storage && navigator.storage.estimate ? await navigator.storage.estimate() : null; }
    catch (e) { return null; }
  };
  /** File > Clear Saved Files: stop Python, delete the saved workspace and settings, start afresh (a reload).
   *  The course progress (localStorage) and the stored app files (sw.js) are kept. */
  window.SASSI_WEB.clearSaved = async function () {
    if (worker) worker.terminate();
    dead = "the saved files are being cleared";
    await Promise.all(KEPT_DBS.map((name) => new Promise((res) => {
      try {
        const r = indexedDB.deleteDatabase(name);
        r.onsuccess = r.onerror = r.onblocked = () => res();
      } catch (e) { res(); }
    })));
    location.reload();
  };
  // save the workspace when the visitor leaves or hides the tab (the worker also saves after every change)
  const flush = () => { if (worker && !dead) worker.postMessage({type: "flush"}); };
  window.addEventListener("pagehide", flush);
  document.addEventListener("visibilitychange", () => { if (document.visibilityState === "hidden") flush(); });

  let dead = null;                             // the message once the worker has stopped
  function stop(message) {
    dead = message;
    fail(message);
    readyReject(new Error(message));
    for (const p of pending.values()) p.reject(new Error(message));
    pending.clear();
  }

  /** The service worker (sw.js, from web/sw.js) keeps the files of the app on this computer: a visit after the
   *  first one downloads nothing.  On the first visit Python starts once the service worker controls the page,
   *  so that the one download is stored (at most a few seconds' wait).  Resolves true when it controls the page. */
  async function serviceWorker() {
    if (!("serviceWorker" in navigator) || !window.isSecureContext) return false;
    try {
      await navigator.serviceWorker.register("sw.js", {scope: "./"});
    } catch (e) {
      console.warn("SASSI-EDU: no service worker (files are cached by the browser only):", e);
      return false;
    }
    if (navigator.serviceWorker.controller) return true;
    await new Promise((res) => {
      const t = setTimeout(res, 4000);
      navigator.serviceWorker.addEventListener("controllerchange", () => { clearTimeout(t); res(); }, {once: true});
    });
    // ask the browser to keep the stored files under storage pressure (Firefox would ask the visitor: not there)
    if (navigator.storage && navigator.storage.persist && !/Firefox\//.test(navigator.userAgent)) {
      navigator.storage.persist().catch(() => {});
    }
    return !!navigator.serviceWorker.controller;
  }
  /** Number of Pyodide files (runtime, NumPy, SciPy ...) an earlier visit stored (sw.js: its Pyodide cache). */
  async function storedFiles() {
    try {
      const name = `sassi-edu-pyodide-${PYODIDE_VERSION}`;
      if (!window.caches || !(await caches.has(name))) return 0;
      return (await (await caches.open(name)).keys()).length;
    } catch (e) {
      return 0;
    }
  }

  if (typeof WebAssembly !== "object" || typeof Worker !== "function") {
    stop("This browser has no WebAssembly or no Web Workers.");
  } else {
    setStage("Preparing ...", 0.01);
    (async () => {
      const [sw, stored] = await Promise.all([serviceWorker(), storedFiles()]);
      note.textContent = !sw ? NO_SW : stored >= 4 ? FROM_CACHE : FIRST_VISIT;
      window.SASSI_WEB.storage = {serviceWorker: sw, storedFiles: stored};
      setStage("Starting the Python engine (Web Worker) ...", 0.02);
      try {
        worker = new Worker(`web/worker.js?v=${BUILD}`, {type: "module"});
      } catch (e) {
        stop(`The Web Worker could not be created: ${e.message || e}`);
        return;
      }
      connect(worker);
    })();
  }
  function connect(worker) {
    worker.onerror = (ev) => {
      ev.preventDefault();
      stop(`The Python engine stopped (${ev.message || "the worker script failed to load"}).`);
    };
    worker.onmessage = (ev) => {
      const m = ev.data || {};
      if (m.id !== undefined) {                       // the answer of a request
        const p = pending.get(m.id);
        if (!p) return;
        pending.delete(m.id);
        try { p.resolve(JSON.parse(m.answer)); } catch (e) { p.reject(e); }
      } else if (m.push !== undefined) {              // events / job output sent while Python is busy
        let data;
        try { data = JSON.parse(m.push); } catch (e) { return; }
        // static/app.js declares SASSI as a global const (a binding of the global scope, not a window property)
        const gui = typeof SASSI !== "undefined" ? SASSI : null;
        if (gui && gui.push) {
          try { gui.push(data); } catch (e) { console.error(e); }
        }
      } else if (m.type === "status") {
        setStage(m.text, m.fraction);
      } else if (m.type === "ready") {
        window.SASSI_WEB.info = m.info || null;
        readyResolve();
        done();
      } else if (m.type === "fatal") {
        stop(m.error || "The Python engine stopped.");
      }
    };
  }

  /** The GUI's transport (static/app.js S.api): one request to the Python session of the worker. */
  window.SASSI_TRANSPORT = function (method, path, body) {
    return ready.then(() => new Promise((resolve, reject) => {
      if (dead) { reject(new Error(dead)); return; }
      const id = ++seq;
      pending.set(id, {resolve, reject});
      const i = path.indexOf("?");
      worker.postMessage({id, method, path: i < 0 ? path : path.slice(0, i), query: i < 0 ? "" : path.slice(i + 1),
        body: body === undefined ? "" : JSON.stringify(body)});
    }));
  };
})();
