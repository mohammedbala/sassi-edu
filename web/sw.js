/* SASSI-EDU in the browser (GitHub Pages) -- the service worker: keeps the files of the app on this computer, so
 * that a visit after the first one downloads nothing (and the site also works offline).
 *
 * Two caches (Cache Storage of the site's origin):
 *  - "sassi-edu-pyodide-<version>": Pyodide (the Python runtime), NumPy, SciPy and their libraries from the
 *    jsDelivr CDN.  The URLs carry the pinned Pyodide version and never change: cache first, and the cache is kept
 *    across deployments of the site (a new SASSI-EDU build does not download Python again);
 *  - "sassi-edu-site": the files of this site, kept across deployments.  Files whose URL carries a version
 *    (?v=<hash of the file>: the scripts, the style sheets, the SASSI-EDU bundle web/sassi-edu.zip) are served
 *    from the cache first -- a new deployment changes the version of the files that changed only, so only those
 *    are downloaded again; the MANIFEST of this build lists the versioned URLs it uses, and the others (files of
 *    earlier builds) are deleted when this service worker takes over.  The page itself goes to the network
 *    first (a new deployment shows at once; the cache answers offline); other files (fonts, pictures, help
 *    figures) come from the cache at once and are refreshed in the background.
 *
 * The Python engine (web/worker.js) also reads and fills the Pyodide cache itself, so the packages are kept in a
 * browser where this service worker cannot run or does not see the worker's requests.
 *
 * Registered by web/boot.js (scope: the site folder).  web/build.py fills in the build and the Pyodide version. */

const BUILD = "__SASSI_BUILD__";
const PYODIDE_VERSION = "__PYODIDE_VERSION__";
const SITE_CACHE = "sassi-edu-site";
const MANIFEST = new Set(__SASSI_MANIFEST__);      // the versioned files of this build: "<path>?v=<hash>"
const PY_CACHE = `sassi-edu-pyodide-${PYODIDE_VERSION}`;
const CDN = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/`;
const SCOPE = new URL("./", self.location.href).pathname;

// install: store the page and every versioned file of this build that the cache does not have yet (after a new
// deployment: only the files that changed), so that the next visit needs no network at all -- on the first
// visit the page loaded its scripts before this service worker ran.  A file that cannot be fetched now is
// stored at its first use instead.
self.addEventListener("install", (ev) => {
  self.skipWaiting();
  ev.waitUntil((async () => {
    try {
      const cache = await caches.open(SITE_CACHE);
      const have = new Set((await cache.keys()).map((r) => r.url));
      const page = new URL("./", self.location.href).href;
      const want = [page, ...[...MANIFEST].map((p) => new URL(p, self.location.href).href)];
      await Promise.all(want.filter((u) => u === page || !have.has(u)).map(async (u) => {
        try { keep(cache, u, await fetch(u)); } catch (err) { /* stored at its first use */ }
      }));
    } catch (err) { /* no Cache Storage: the service worker only passes the requests on */ }
  })());
});

self.addEventListener("activate", (ev) => {
  ev.waitUntil((async () => {
    // other Pyodide versions and the per-build caches of earlier service workers
    for (const k of await caches.keys()) {
      if (k.startsWith("sassi-edu-") && k !== SITE_CACHE && k !== PY_CACHE) await caches.delete(k);
    }
    // versioned files this build no longer uses
    const site = await caches.open(SITE_CACHE);
    for (const req of await site.keys()) {
      const u = new URL(req.url);
      if (u.searchParams.has("v") && !MANIFEST.has(u.pathname.slice(SCOPE.length) + u.search)) await site.delete(req);
    }
    await self.clients.claim();               // control the open page at once (web/boot.js waits for it)
  })());
});

function keep(cache, req, res) {
  // only complete, successful answers; a full disk (quota) must not break the page
  if (res && res.ok && res.status === 200) cache.put(req, res.clone()).catch(() => {});
}

async function cacheFirst(req, name) {
  const cache = await caches.open(name);
  const hit = await cache.match(req);
  if (hit) return hit;
  const res = await fetch(req);
  keep(cache, req, res);
  return res;
}

async function networkFirst(req, name) {
  const cache = await caches.open(name);
  try {
    const res = await fetch(req);
    keep(cache, req, res);
    return res;
  } catch (err) {
    const hit = await cache.match(req);
    if (hit) return hit;
    throw err;
  }
}

async function staleWhileRevalidate(ev, req, name) {
  const cache = await caches.open(name);
  const hit = await cache.match(req);
  const update = fetch(req).then((res) => { keep(cache, req, res); return res; });
  if (hit) {
    ev.waitUntil(update.catch(() => {}));
    return hit;
  }
  return update;
}

self.addEventListener("fetch", (ev) => {
  const req = ev.request;
  if (req.method !== "GET") return;
  const url = new URL(req.url);
  if (url.href.startsWith(CDN)) {
    ev.respondWith(cacheFirst(req, PY_CACHE));
    return;
  }
  if (url.origin !== self.location.origin || !url.pathname.startsWith(SCOPE)) return;
  if (req.mode === "navigate" || url.pathname.endsWith("/") || url.pathname.endsWith(".html")) {
    ev.respondWith(networkFirst(req, SITE_CACHE));
  } else if (url.searchParams.has("v")) {        // a versioned file never changes under its URL
    ev.respondWith(cacheFirst(req, SITE_CACHE));
  } else {
    ev.respondWith(staleWhileRevalidate(ev, req, SITE_CACHE));
  }
});
