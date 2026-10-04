# SASSI-EDU in the browser (GitHub Pages)

This folder builds a **static web site** of the SASSI-EDU GUI: the guided course, the examples, the
plots, Help and the Command Entry, with **no server**. Python runs in the visitor's browser with
[Pyodide](https://pyodide.org) (CPython 3.14 compiled to WebAssembly, with NumPy and SciPy), inside a
Web Worker. The same front end (`sassi/ui/static`) and the same API (`sassi/ui/api.py`) serve the local
`sassi-gui` and the web version.

| File | Role |
|---|---|
| `build.py` | writes the site to `web/site/` (git-ignored) |
| `boot.js` | page side: starts the worker, loading overlay with the disclaimer, `window.SASSI_TRANSPORT`, pushes |
| `worker.js` | the Web Worker: Pyodide 314.0.7 from the jsDelivr CDN, NumPy, SciPy, the bundle, `sassi.web.bridge` |
| `test_pyodide.mjs` | Node test: the built bundle in Pyodide, driven like the page (lesson 01, all lessons headless) |
| `package.json` | the pinned `pyodide` npm package for that test |
| `../sassi/web/bridge.py` | the Python side: requests, inline module runs, pushes |
| `../.github/workflows/pages.yml` | builds and deploys the site on every push to `main` |

The design is described in `docs/internal/web_version.md` (not published).

## Build

From the project root (needs Python 3.9+ and the `plotly` package, which provides Plotly.js; the project
virtual environment has both):

```bash
.venv/bin/python web/build.py
```

It prints the size of the site (about 7.7 MB: Plotly.js 4.8 MB, the Python bundle `web/sassi-edu.zip`
1.9 MB, KaTeX 0.6 MB; about 3.9 MB over the network, as GitHub Pages compresses text files). Visitors
also download Pyodide, NumPy and SciPy from cdn.jsdelivr.net at their first visit (23 MB compressed,
32 MB unpacked), which their browser then caches: about 27 MB in all.

The build never bundles `reference/` (the ACS SASSI manual text), `docs/spec/`, `docs/internal/`,
`tests/`, `examples/sassi-course/` (your local course runs) or virtual environments; it stops with an
error if one of them would be written. `tests/unit/test_web.py` checks the bundle as well.

## Test locally

Any static file server will do. A Web Worker cannot start from a `file://` page, so do not open
`index.html` directly:

```bash
.venv/bin/python -m http.server 8811 -d web/site
```

Then open <http://127.0.0.1:8811/>. Python is ready after about 5 to 15 s at the first start (downloads,
then Python starts; less once the browser has cached the files); the overlay shows what is loading. The
browser needs internet access for the Pyodide files.

The Node test runs the built bundle in Pyodide without a browser:

```bash
cd web && npm install && cd ..          # once: the pinned pyodide package (about 30 MB, git-ignored)
node web/test_pyodide.mjs               # lesson 01 through the bridge, then every lesson headless (~2.5 min)
node web/test_pyodide.mjs --quick       # lesson 01 only (~20 s); pytest runs this one
```

`npm install` gets the Pyodide runtime; NumPy and SciPy are downloaded from the CDN at the first run and
kept in `web/node_modules/pyodide/`.

## Publish on GitHub Pages

The site is free on GitHub Pages for a **public** repository (a private repository needs a paid plan).

1. **Check what the repository will contain.** `.gitignore` keeps out `.venv/`, caches, `web/site/`,
   `web/node_modules/`, `reference/acs-sassi.txt` and the course workspaces (`examples/sassi-course/`).
   It does *not* keep out `docs/spec/` and `docs/internal/`: the specifications digest the ACS SASSI
   manual. They are not part of the web site, but a public repository shows them. If you do not want
   them public, add `docs/spec/` and `docs/internal/` to `.gitignore` before the first commit (the
   build and the site do not need them).
2. **Create the repository** on GitHub (for example `sassi-edu`), without a README.
3. **Push** the project from its root:

   ```bash
   git init -b main
   git add .
   git status                      # check: no reference/, .venv/, web/site/, examples/sassi-course/
   git commit -m "SASSI-EDU"
   git remote add origin https://github.com/<user>/sassi-edu.git
   git push -u origin main
   ```

4. **Turn on Pages:** on GitHub, *Settings > Pages > Build and deployment > Source: GitHub Actions*.
5. The workflow *Deploy the web version to GitHub Pages* (Actions tab) builds and deploys the site;
   re-run it once if the first push happened before step 4. The site is at
   `https://<user>.github.io/sassi-edu/`.

Every later push to `main` publishes a new version. Each build has its own id in the URLs of the page,
so a browser never mixes files of two builds.

## What the web version cannot do

* **Files live in the browser tab.** The workspace (`/home/pyodide/work`) is in memory: a reload or
  closing the tab starts afresh. Your course progress is kept (browser storage). Use **File > Upload to
  Workspace...** and **File > Download...** to copy text files in and out; binary files (FILE1 ... FILE8,
  `.sdb`) cannot be transferred.
* **No Cancel of a running module.** Module runs execute in the Python engine of the page, one at a time;
  their listing and progress are shown live, and other requests wait until the run ends.
* **Not available:** Model > Exit, Modules > Location (external module executables), anything that needs a
  process. File > Export Image saves the active plot as PNG with Plotly instead of `CAPTUREPLOT`.
* **Speed.** The lessons take 1.2 to 1.6 times as long as with native Python (lesson 1: 14.6 s against
  11.4 s); Help > Verification runs, but slowly.
