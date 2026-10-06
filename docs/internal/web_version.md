# SASSI-EDU in the browser (GitHub Pages) -- design contract (internal)

Goal (the user's request): "do the github pages in browser version" -- a public, free, static web site
on GitHub Pages where anyone can open SASSI-EDU (the GUI with the guided course, the examples, the
plots, Help) with **no server**: the Python package runs in the visitor's browser with
[Pyodide](https://pyodide.org) (CPython compiled to WebAssembly) in a **Web Worker**.

## Feasibility (measured)

* Pyodide **314.0.7** (CPython 3.14.2, numpy 2.4.6, scipy 1.18.0) runs the package unchanged: lesson
  01 runs headlessly in Node (`npm install pyodide`) in 16.9 s (12 s native). Only numpy and scipy are
  needed (matplotlib is imported lazily, plotly only server-side for `plotly.min.js`).
* What cannot run in a browser: worker **processes** (`subprocess`, `sassi/ui/jobs.py` +
  `sassi/ui/worker.py`), **threads** (`threading.Thread.start` fails in Pyodide), blocking long-polls
  (`/api/events?timeout=20` would freeze the worker), external module executables, the HTTP server.

## Architecture

```
index.html (static)  --  static/app.js dialogs.js plots.js learn.js (unchanged GUI)  --  S.api(method, path, body)
        |                                        |  web mode: window.SASSI_TRANSPORT(method, path, body) -> Promise<json>
web/boot.js  (main thread: starts the worker, loading overlay, transport, push handling)
        |  postMessage {id, method, path, query, body}  /  {id, status, payload}  /  {push: ...}
web/worker.js (module worker: loadPyodide from the jsDelivr CDN, numpy+scipy, unpack web/sassi-edu.zip,
        |      import sassi.web.bridge, answer requests, then run queued jobs)
sassi/web/bridge.py  -- Bridge(push): GuiSession(web mode) + route() + run_pending()
```

* **Same code, two transports.** The GUI's JS and the API (`sassi/ui/api.py`, `GuiSession.route`) are
  shared by the local server and the web version; web mode is selected by `window.SASSI_TRANSPORT`
  (set by `web/boot.js`) on the JS side and by `GuiSession(..., web=True)` (or an injected job manager) on
  the Python side. The local `sassi-gui` must keep working exactly as now (all existing tests pass).
* **Relative URLs.** A GitHub Pages project site lives under `https://<user>.github.io/<repo>/`, so the
  page must not use root-absolute URLs: `index.html` and the JS use `static/...` (not `/static/...`); help
  figures and document links use `docs/...` (helpdocs resolver; `markdown.SAFE_HREF` accepts it). The
  local server keeps serving these (the page is at `/`, so relative = absolute there).
* **Jobs inline.** `sassi/ui/jobs.py` gets an `InlineJobManager` (same interface as `JobManager`:
  `jobs`, `running()`, `get()`, `start()`, `cancel()`, `shutdown()`), used in web mode. `start()` only
  queues the job (state `starting`) and returns, so the request that started it answers at once (the GUI
  opens the job's output tab); `run_pending()` (called by the worker right after it has posted the
  response) runs it **synchronously in-process** with the same semantics as `sassi/ui/worker.py`: a fresh
  interpreter holding a copy of the session's models (`to_json`/`from_json`), the AFWRITE state, the
  progress wrapper of `run_module` -- refactor `worker.main` into a reusable `run_spec(spec, send)` used
  by both. Messages, progress and the finish go through the same callbacks as now (`_message`,
  `on_progress`, `_finish`), so the Command History, the replay history and the commands deferred after a
  `RUN<MODULE>` behave as in the local GUI. Cancel of a running inline job is not possible (documented);
  a cancel before it starts works.
* **Events without blocking.** In web mode `/api/events` is answered with `timeout=0`; the JS poll loop in
  web mode polls without a long-poll (e.g. every 250 ms when idle, at once after a push). While a job runs
  the worker is busy; **pushes** (`postMessage` from Python via a `push` callback given to the Bridge) carry
  the job's streamed messages and progress so the job's output tab and the status bar update live
  (de-duplicate against the later `/api/jobs/<id>?since=` answer by message index).
* **Files.** The workspace lives in the worker's in-memory file system (`/home/pyodide/work`, start
  directory; course workspaces under it). Add to the GUI (web mode only or both): **File > Upload to
  Workspace...** (text files from the visitor's computer -> `POST /api/file`) and **File > Download...**
  (a text file of the workspace via `GET /api/file` -> Blob download), and a short note that files live in
  this browser tab only (a reload starts fresh; course progress in localStorage persists).
* **Not in the web version** (hidden or explained): Model > Exit, Modules > Location (external
  executables), Cancel of a running module, anything needing a process. VERIFY runs inline (slow, allowed).
* **Bundle.** `web/build.py` writes the static site to `web/site/` (git-ignored):
  `index.html` (from `sassi/ui/static/index.html`: token removed, `web/boot.js` added), `static/*` (incl.
  `katex/`), `static/plotly.min.js` (from the installed plotly package), `web/boot.js`, `web/worker.js`,
  `web/sassi-edu.zip` (the `sassi` package without `__pycache__`, `README.md`, `examples/` without
  `examples/sassi-course/`, the public docs: `docs/index.md`, `docs/user`, `docs/theory`,
  `docs/verification`, `docs/reference`), the help figures under `docs/...` (PNG/JPEG/GIF), `.nojekyll`.
  **Never** bundle `reference/` (the copyrighted ACS SASSI manual text), `docs/spec/`, `docs/internal/`,
  `tests/`, `.venv`. The build prints the site size; keep the zip small (a few MB).
* **Pyodide version** pinned (314.0.7, CDN `https://cdn.jsdelivr.net/pyodide/v314.0.7/full/`), the same
  version the Node test uses.
* **Loading screen.** First visit downloads ~30 MB (Pyodide, numpy, scipy; kept by the service worker sw.js,
  afterwards): show a progress overlay with what is loading, and a clear error if WebAssembly/workers are
  unavailable. Show the disclaimer (educational, not affiliated with ACS SASSI, not for design) on the
  overlay and in Help > About.
* **Deployment.** `.github/workflows/pages.yml`: on push to `main`, set up Python, `pip install plotly`
  (for `plotly.min.js`), run `python web/build.py`, upload `web/site` with `actions/upload-pages-artifact`
  and deploy with `actions/deploy-pages`. `.gitignore` for the repository: `.venv/`, `__pycache__/`,
  `*.egg-info/`, `web/site/`, `web/node_modules/`, `reference/acs-sassi.txt`, `examples/sassi-course/`,
  `sassi-course/`. `web/README.md`: how to build, test locally (`python -m http.server -d web/site`), and
  publish (create the repository, push, Settings > Pages > Source: GitHub Actions).

## Tests

* Python (always): `InlineJobManager` and the bridge with a fake push -- a `RUN<MODULE>` started by
  `POST /api/command` returns a `starting` job, `run_pending()` runs it, the messages/progress were pushed,
  `/api/jobs/<id>` reports `done`, deferred commands ran, the replay history is right; `/api/events` never
  blocks in web mode; relative `docs/` links; the build script (bundle contents: no `reference/`,
  `docs/spec`, `docs/internal`, `tests`, `examples/sassi-course`; every file the page references exists).
* Pyodide (when `web/node_modules/pyodide` exists, `npm install` in `web/`): `web/test_pyodide.mjs` loads
  the built zip into Pyodide in Node and drives the Bridge like the browser: open lesson 01, run its steps
  (commands one per request, jobs via `run_pending`), action plots via the API; and runs every lesson
  headlessly (`run_lesson_headless`). Report the times.
* Browser: serve `web/site` locally and use it in a real browser (Pyodide from the CDN): start page, a
  lesson step with module runs (live output), plots, Help with KaTeX, the explainer, upload/download.

## As built (2026-10-04)

* Python: `sassi/web/bridge.py` (`Bridge(push)`: `handle(method, path, query, body) -> JSON text
  {status, payload}`, `pending`, `run_pending()`, `flush()`); `GuiSession(web=True)`;
  `jobs.InlineJobManager` (shares `_Jobs` with `JobManager`); `worker.run_spec(spec, send)` with
  `progress_reporting()`, which restores `run_module` / the P2 runners after the job;
  `EventLog.subscribe()` (the bridge pushes at most every 0.1 s while a request or a job runs; the first
  push of a request not before 0.1 s).
* Push format `{"events", "last", "reset", "job"}`; the page skips events by `seq` and job messages by
  index (`app.js handleEvents / applyJob`); the bridge tracks the newest event it has sent (answers of
  `/api/events` and pushes), so pushes have no gap.
* JS: `S.transport` / `S.web`, `webApi`, `S.push`, `webWait` (250 ms poll, woken by a push); menu items
  `web: false` (Model > Exit, Modules > Location) and `webOnly` (File > Upload to Workspace...,
  Download...); job Cancel hidden; File > Export Image = `Plotly.downloadImage` in web mode.
  `boot.js` reaches the GUI through the global `SASSI` binding (a `const`, not a `window` property).
* Bundle: `sassi/ui/static/` is not in the zip (the site has it); the Help figures are both in the zip
  (`image_file()` checks them) and in the site. `index.html` gets a CSP meta tag and `?v=<build id>` on
  every local URL (hash of the content; also in `web/worker.js?v=` and `sassi-edu.zip?v=`).
* Measured: bundle zip 1.93 MB, site 7.7 MB (3.9 MB gzip), CDN 23.1 MB (brotli); Node: Pyodide +
  NumPy + SciPy + bundle 3.5 s, lesson 01 through the bridge 15.9 s, the ten lessons headless 128 s
  (native 1.2-1.6 x faster); browser: ready in 3-8 s with a warm cache. Lesson 11 (2026-10-05): 236 s
  headless (93 s native: 2.5 x; three ANALYS runs of 41 frequencies with 627 to 1215 interaction DOFs);
  the eleven lessons together take about 6 min.

## Files kept between visits, files of the visitor's computer (2026-10-06)

* `web/worker.js` mounts `/home/pyodide/work` and `/home/pyodide/.sassi-edu` on Emscripten IDBFS (one
  IndexedDB database per mount point) before the session starts and restores them (`syncfs(true)`); it
  saves (`syncfs(false)`) at most 1 s after a request that is not a GET (the 250 ms event poll does not
  save) and after every job, and at once on `{type: "flush"}`, which `web/boot.js` sends on `pagehide` and
  when the page is hidden. Without IndexedDB (a private window may refuse it) the files stay in memory.
  `info.kept` reports the files restored. File > Clear Saved Files: `SASSI_WEB.clearSaved()` terminates the
  worker, deletes the two databases and reloads; the course progress (localStorage) and the stored app files
  (`sw.js`) stay.
* Every file dialog (`D.pickFile`, web version, not folder pickers or Download) has **From your
  computer…**: text files (<= 8 MB, no NUL in the first 4 kB) copied into the folder shown (`POST
  /api/file`), the first selected. Model > Output and Export Table download the file they wrote
  (`GET /api/file`), as Export to ANSYS already did.
