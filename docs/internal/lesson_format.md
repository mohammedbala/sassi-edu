# SASSI-EDU guided course: lesson file format (contract)

The guided course ("Learn" in the GUI) is a set of Markdown lesson files in `sassi/ui/lessons/`,
parsed by `sassi/ui/lessons.py` (`parse_lesson`, `load_lessons`) and validated by actually running
them (`run_lesson_headless`, used by `tests/unit/test_lessons.py`). The GUI renders the lessons
with its Markdown renderer (`sassi/ui/markdown.py`) and executes their commands through the session
interpreter (rule L17), so everything a lesson does is visible in the Command History and replayable.

## Audience and tone

Structural engineers who design nuclear (and other critical) facilities: they know seismic design
(ASCE 4 / ASCE 43 / ACI 349 / SRP 3.7 practice), response spectra and in-structure response spectra
(ISRS), fixed-base finite-element models (mostly ANSYS), modal analysis, damping, member design. They do
**not** know SSI theory or SASSI. Write for them: start from what they know (fixed-base model, base
input, springs and dashpots, harmonic analysis in ANSYS), explain what SSI changes and why it matters for
their design deliverables (ISRS, member forces, relative displacements, soil pressures), and give the
technical basis precisely (equations, assumptions, limits). Never pad; never oversell. Every number you
state must come from running the lesson (or from the cited source).

## File layout

`sassi/ui/lessons/NN_slug.md` — NN orders the files. UTF-8 Markdown with a front-matter header:

```
---
id: 04-surface-ssi
title: Your first SSI analysis: a stick on a surface mat
part: Fundamentals
order: 4
minutes: 25
example: ex01_surface_stick
summary: One sentence shown in the course outline.
objectives: [Build a stick model, Run the SASSI chain, Read ISRS and compare with fixed base]
prerequisites: [02-free-field, 03-impedance]
---
```

* `id` (required, unique, `[0-9a-z-]+`), `title` (required), `part` (`Fundamentals`, `Design
  applications`, `Advanced`), `order` (int), `minutes` (estimated time), `summary`, `objectives` and
  `prerequisites` (lists in `[a, b]` form, items without commas), `example` (optional: name of an
  `examples/*.pre` file without extension; its `.pre` and the `examples/data/` folder are copied into the
  lesson workspace).
* Values are plain text after `key:`; lists use `[ ... ]` with comma-separated items.

## Body

* Text before the first `##` heading is the lesson **introduction** (shown on the lesson's first page).
* Every `## ` heading starts a **step**. The step's text up to its first `###` heading is the step
  **narrative**.
* Inside a step, `### ` headings start **sections**. Recognised titles (case-insensitive; the GUI renders
  each as a labelled panel):
  - `What this does` — the commands of the step explained argument by argument where useful.
  - `Why it matters` — the design/engineering significance (ISRS, member forces, licensing practice).
  - `Technical basis` — equations, assumptions, method details; cite the theory manual section
    (`docs/theory/THEORY_MANUAL.md#...`) and the verification problem that verifies it (VP-xx).
  - `In ANSYS terms` — the analogy for an ANSYS user (optional).
  - `Try this` — a variation the learner can make and re-run (optional).
  - `Check yourself` — a question; the line starting with `Answer:` and everything after it is the
    answer (the GUI hides it behind "Show answer").
  Other `###` titles are allowed and rendered as plain sections.

## Fenced blocks

| Block | Meaning |
|---|---|
| ```` ```sassi ```` | Commands the step runs when the learner presses **Run step**, in order (one command per line; `*` comment lines are echoed as comments). A step may have several; they run in document order. |
| ```` ```sassi-show ```` | Commands or syntax shown but not run (e.g. `SITE,<opmode>,<mode1>,...` templates or a variation). |
| ```` ```sassi-setup ```` | Only in the introduction: commands run automatically when the lesson is opened, after the workspace is prepared (e.g. to restore the state a lesson builds on). Keep them fast. |
| ```` ```action ```` | GUI actions offered as buttons after the step has run, one per line, `verb: arguments` (below). |
| any other fence | Displayed as code (equations in plain text, APDL, file excerpts). |

Action verbs (paths are relative to the lesson workspace):

| Verb | Arguments | Effect |
|---|---|---|
| `plot-model` | — | Plot ▸ Model ▸ Elements of the active model |
| `plot-nodes` | — | Plot ▸ Model ▸ Nodes |
| `plot-layers` | — | Plot ▸ Soil Layers (LAYERPLOT) |
| `plot-soilprops` | `<DYNP label>` | Plot ▸ Soil Properties |
| `plot-spectrum` | `file1, file2, ... [| log]` | READSPEC each file and SPECPLOT them (log frequency axis with `| log`) |
| `plot-history` | `file1, file2, ...` | READTH each file and THPLOT them |
| `open-file` | `path` | open a text file (listing, deck, result) in an editor tab |
| `open-listing` | `MODULE` | open `<model dir>/<model>_<MODULE>.out` |
| `open-dialog` | `ANALYSIS[/TAB]` or `MODEL`/`WRITE`/`CHECK` | open Options ▸ Analysis at a tab, etc. |
| `open-doc` | `docs/...md[#anchor]` | open a document in Help |
| `explain` | `<command line>` | open the command explainer for that line |

## Equations (LaTeX)

Every equation and every formula in the prose is written in LaTeX and typeset by KaTeX in the GUI
(`sassi/ui/static/katex`; the renderer is `sassi/ui/markdown.py`):

* inline: `$\omega^2 M_s$`, `$a_0 = \omega B / V_s$` (no space after the opening or before the closing
  `$`; a literal dollar is `\$`);
* display: a fenced block with the language `math` (preferred), or `$$ ... $$` on lines of their own.
  Several lines go in one `aligned` (or `gathered`) environment:

  ````
  ```math
  \begin{aligned}
  \left[(K^*_s-\omega^2 M_s)-(K^*_e-\omega^2 M_e)+X_{ff}\right]U &= X_{ff}\,U'_f\\
  X_{ff} &= F_{ff}^{-1}
  \end{aligned}
  ```
  ````
* notation (as the Theory Manual): $K^*_s, M_s$ structure, $K^*_e, M_e$ excavated soil, $F_{ff}$
  flexibility and $X_{ff}$ impedance at the interaction DOFs, $U'_f$ free-field motion, $\omega$, $f$,
  $\beta$ damping ratio, $V_s, V_p$, $\rho$, $G$, $\nu$, $a_0 = \omega B/V_s$, $H(f)$ transfer
  function; units in `\text{...}` with a thin space (`300\,\text{m/s}`);
* KaTeX supports most of LaTeX math (no `\usepackage`, no `\newcommand` across formulas); inside a
  table cell write `\lvert x \rvert` instead of `|x|` (a bare `|` ends the cell);
* command names, file names and SASSI command lines stay in code spans / command blocks, never in math;
* the definitions of the symbols of a display equation go in a sentence or a short list after it
  (not inside the equation as text);
* check: `.venv/bin/python -m sassi.ui.mathcheck sassi/ui/lessons/*.md` (runs the GUI's KaTeX in
  Node.js; `tests/unit/test_lessons.py` runs the same check).

## Workspace and state

When a lesson is opened the GUI creates `<workspace root>/<lesson id>/` (default workspace root:
`<GUI start directory>/sassi-course/`), copies the example files when `example` is set, changes the
interpreter's working directory there (`CD`), starts a fresh model (the learner is asked before an
unsaved model is replaced) and runs the `sassi-setup` block. Lessons are self-contained: a lesson that
needs the results of another runs what it needs in `sassi-setup` (e.g. `INP` of the example) — keep
the total runtime of any lesson below about two minutes on a laptop.

## Validation (`tests/unit/test_lessons.py`)

* the front matter is complete and ids are unique; prerequisites exist;
* every section title is known or deliberate, every action verb is known, every runnable command's
  name exists in the interpreter's command catalogue;
* running each lesson headlessly (setup, then every step's `sassi` blocks in order, in a fresh
  temporary workspace) produces **no error messages**, and every file named by an action exists after
  its step;
* every link to `docs/` resolves.
