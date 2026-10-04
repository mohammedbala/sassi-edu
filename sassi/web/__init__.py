"""SASSI-EDU in the browser: the static web version of the GUI (GitHub Pages).

The GUI of :mod:`sassi.ui` normally talks to the local server ``sassi-gui``.  The web version has no
server: the page (``web/boot.js``) starts a Web Worker (``web/worker.js``) that runs this package with
Pyodide (CPython compiled to WebAssembly, with NumPy and SciPy) and passes it the requests of the page.
:mod:`sassi.web.bridge` answers them with one :class:`sassi.ui.api.GuiSession` in web mode (module runs
inline, no long poll) and pushes the session's events to the page while Python is busy.

``web/build.py`` writes the static site (``web/site/``); ``web/README.md`` tells how to test and publish it.
"""
