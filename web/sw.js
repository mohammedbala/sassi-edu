/* SASSI-EDU in the browser (GitHub Pages) -- the service worker: keeps the files of the app on this computer, so
 * that a visit after the first one downloads nothing (and the site also works offline).
 *
 * Two caches (Cache Storage of the site's origin):
 *  - "sassi-edu-pyodide-<version>": Pyodide (the Python runtime), NumPy, SciPy and their libraries from the
 *    jsDelivr CDN.  The URLs carry the pinned Pyodide version and never change: cache first, and the cache is kept
 *    across deployments of the site (a new SASSI-EDU build does not download Python again);
 *  - "sassi-edu-site-<build>": the files of this site.  Files whose URL carries this build (?v=<build>: the
 *    scripts, the style sheets, the SASSI-EDU bundle web/sassi-edu.zip) are served from the cache first; the page
 *    itself goes to the network first (a new deployment shows at once; the cache answers offline); other files
 *    (fonts, pictures, help figures) come from the cache at once and are refreshed in the background.
 *    The caches of other builds are deleted when a new service worker takes over.
 *
 * Registered by web/boot.js (scope: the site folder).  web/build.py fills in the build and the Pyodide version. */

const BUILD = "__SASSI_BUILD__";
const PYODIDE_VERSION = "__PYODIDE_VERSION__";
const SITE_CACHE = `sassi-edu-site-${BUILD}`;
const PY_CACHE = `sassi-edu-pyodide-${PYODIDE_VERSION}`;
const CDN = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/`;
const SCOPE = new URL("./", self.location.href).pathname;

self.addEventListener("install", () => { self.skipWaiting(); });

self.addEventListener("activate", (ev) => {
  ev.waitUntil((async () => {
    for (const k of await caches.keys()) {
      if (k.startsWith("sassi-edu-") && k !== SITE_CACHE && k !== PY_CACHE) await caches.delete(k);
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
  } else if (url.searchParams.get("v") === BUILD) {
    ev.respondWith(cacheFirst(req, SITE_CACHE));
  } else {
    ev.respondWith(staleWhileRevalidate(ev, req, SITE_CACHE));
  }
});
