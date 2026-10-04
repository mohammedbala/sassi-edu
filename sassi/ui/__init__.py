"""SASSI-EDU graphical user interface: a local web application (ARCHITECTURE section 9).

The system Tk on macOS is 8.5.9 and deprecated, so the GUI of requirements section 5 is a small
web application served by the Python standard library and shown in the user's browser
(decision D-GEN-05 revised).  Everything runs on the user's machine: the server binds to
``127.0.0.1`` only and Plotly.js is served from the installed ``plotly`` package (works offline).

Start it with ``sassi-gui`` (or ``python -m sassi.ui.server``)::

    sassi-gui [--port 8765] [--no-browser] [--model-dir DIR]

The same front end and API also run with no server at all, in the visitor's browser (Python compiled
to WebAssembly by Pyodide, in a Web Worker): the static web version for GitHub Pages, built by
``web/build.py`` -- see :mod:`sassi.web.bridge` and ``web/README.md``.

Modules of this package
-----------------------
``server.py``   ``main()``, the HTTP request handler, static files, security checks.
``api.py``      :class:`~sassi.ui.api.GuiSession` -- one interpreter per server session and the JSON
                API of ARCHITECTURE section 9 (commands, state, model, options, files, plots ...).
``dialogs.py``  the Options dialogs (Model, Write, Check, Analysis with its 12 tabs) and the Modules-menu
                dialogs ANSYS Eq. Static Load / ANSYS Dynamic Load: form layout built from
                :mod:`sassi.prep.options`, current values, and the command text a dialog commit
                submits (rule L17) after validation with the CHECK rules (UI-06).
``cmdrecords.py`` the dialog bindings of the command records outside the option records (Option NON
                EQL / P / S / BBC, SOIL-NON NLSOIL / NLSLAYER, Option A LOADGEN ...), validated by a
                dry run of their commands on a copy of the model.
``jobs.py``     background module runs (``RUN<MODULE>``) and ``VERIFY`` in a worker process with the
                listing streamed and a Cancel that terminates the process (UI-04); inline jobs of
                the browser version.
``worker.py``   the worker process (``python -m sassi.ui.worker``) and ``run_spec``, the job itself.
``events.py``   the event log the browser polls (messages, plot events, progress, jobs).
``files.py``    model-directory listing with result types and path safety.
``modeldata.py`` model JSON for 3D plotting (nodes, elements by group, flags, masses, fixities).
``settings.py`` ``SASSIini.xml`` (module locations, extensions, Command Display) and ``SASSIdb.xml``
                (Load Model tree) persistence (requirements 5.10, UI-07, D-UI-07).
``helpdocs.py`` Help > Help: the documentation set (``docs/``, ``examples/README.md``), links between
                documents, search, the figures served at ``/docs/<path>``.
``markdown.py`` the safe Markdown -> HTML renderer of the Help pages (no raw HTML passes through).
``lessons.py``  Learn: the lesson files of the guided course (``lessons/*.md``; parser, workspace, headless
                runner, validator; format ``docs/internal/lesson_format.md``).
``learn.py``    Learn: course outline, lesson rendering, lesson / example workspaces, step actions and
                the examples gallery (``/api/lessons``, ``/api/examples``).
``explain.py``  Learn: the command explainer (``/api/explain``) with the curated argument table
                ``command_args.json``.
``static/``     the front end: ``index.html``, ``styles.css``, ``app.js`` (main window, menus,
                command entry, tabs), ``dialogs.js`` (dialogs), ``plots.js`` (Plotly drawing),
                ``learn.js`` (Learn tab, lesson panel, command explainer, examples).

Rule L17 / UI-01: every GUI action that changes the model or the plots submits command text to
the one :class:`sassi.prep.Interpreter` of the session, so a GUI session can be replayed as a
``.pre`` file (the Command History is that replay).  Only the session settings that have no
command in ACS SASSI (Options > Write, Options > Check, Command Display filters, module locations
and extensions) are set directly.
"""
