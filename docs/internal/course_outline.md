# SASSI-EDU guided course: outline and style guide (internal)

Lesson format: `docs/internal/lesson_format.md`. Files: `sassi/ui/lessons/NN_slug.md`.

## Learner

A structural engineer who designs nuclear (or other safety-related) facilities: comfortable with
seismic design practice (design response spectra, ISRS, ASCE 4 / ASCE 43, ACI 349, NRC SRP 3.7),
fixed-base FE models (typically ANSYS), modal and harmonic analysis, damping, member design. New to
SSI theory and to SASSI. They want to know: what SASSI does at each step, why, what the numbers mean,
what to check, and how it changes their design deliverables compared with a fixed-base model.

## Style

* Start each lesson from what the learner already knows (fixed-base model, base input, springs).
* One idea per step; 3–8 commands per step; explain arguments that matter (not every zero).
* Sections per step: narrative, `What this does`, `Why it matters`, `Technical basis` (equations as plain
  text code blocks, the assumption behind them, the theory-manual section and the VP that verifies it),
  optionally `In ANSYS terms`, `Try this`, `Check yourself` (question + `Answer:`).
* Quote numbers only from what the lesson actually produces when run (read the listings/results), or
  from the cited source. Never invent code clauses or results; when citing the ACS SASSI manual,
  ASCE 4, ASCE 43 or SRP 3.7.x, cite only what appears in `reference/acs-sassi.txt` or `docs/spec/`.
* Link documents with repository-root paths, e.g. `[Theory §3](docs/theory/THEORY_MANUAL.md#3-flexible-volume-substructuring)`.
* Keep every lesson's total run time under about two minutes on a laptop.
* Use actions (`plot-model`, `plot-spectrum: ex01/00085TR_X01.RS, ex01/00041TR_X01.RS | log`,
  `open-listing: ANALYS`, `open-dialog: ANALYSIS/SITE`, ...) so the learner sees results right away.

## Outline (ids, titles, examples are fixed; content authors fill them)

| # | id | Part | Title | Example | Content |
|---|---|---|---|---|---|
| 1 | `01-why-ssi` | Fundamentals | From fixed base to SSI: what changes and why | ex01_surface_stick | fixed-base vs SSI; kinematic vs inertial interaction; frequency domain; substructuring; the SASSI module chain (manual Fig. 1.1) mapped to the engineer's workflow; a first look at a model and its results |
| 2 | `02-free-field` | Fundamentals | The site: layers, half-space and free-field motion | ex01_surface_stick | L/TOPL layering, sublayer thickness rule (Vs/(5h)), half-space by the variable-depth method, control point and within vs outcrop motion, SITE Mode 1/2, free-field transfer functions |
| 3 | `03-impedance` | Fundamentals | Foundation impedance: soil springs and dashpots | ex03_forced_vibration | interaction nodes, POINT and the central-zone radius, rigid mat, FORCE/ANALYS vibration analysis and global impedance (FOUNSTIF/FOUNDASH), static stiffness vs closed forms, frequency dependence and radiation damping |
| 4 | `04-surface-ssi` | Fundamentals | Your first SSI analysis: a stick on a surface mat | ex01_surface_stick | structure (BEAMS, masses, rigid links), frequency set (df, NFFT, cut-off), HOUSE, ANALYS (Eq. 2.1), MOTION (interpolation, ISRS), RELDISP, STRESS; frequency shift and radiation damping vs fixed base |
| 5 | `05-embedded` | Fundamentals | Embedded structures and the flexible-volume method | ex02_embedded_box | excavated soil (ETYPE), FV vs FI-FSIN vs FI-EVBN vs FFV interaction sets (INTGEN), kinematic interaction, the subtraction-method anomaly, validation of reduced methods against FV |
| 6 | `06-seismic-input` | Design applications | Seismic input: spectrum-compatible motion and strain-compatible soil | ex04_site_response | EQUAKE (target RS, SRP 3.7.1 checks, PSD), SOIL (SHAKE equivalent-linear iterations, curves, strain-compatible properties, FILE88), SITEX, best-estimate/lower/upper-bound soil cases |
| 7 | `07-three-components` | Design applications | Three earthquake directions and design ISRS | ex05_xyz_simultaneous | SV/SH/P vertically propagating waves for X/Y/Z, simultaneous cases, directional combination of ISRS (SRSS / 100-40-40 as the manual or ASCE 4 state), envelope and peak broadening (BROADEN), the design ISRS deliverable |
| 8 | `08-design-outputs` | Design applications | From SSI results to design: forces, displacements, ANSYS | ex01_surface_stick | STRESS element forces (EOUT), RELDISP relative displacements, section cuts, Option A LOADGEN to ANSYS (equivalent static and dynamic loads) |
| 9 | `09-nonlinear` | Advanced | Nonlinear soil and cracked concrete by iteration | ex06_nonlinear_soil | primary vs secondary soil nonlinearity, near-field soil iterations (.pin, FILE74/78), Option NON concrete cracking with panels and backbone curves (ex07, described and/or run within the time budget), when these are needed |
| 10 | `10-good-practice` | Advanced | Good practice, checks and the limits of SSI models | ex01_surface_stick | frequency selection and interpolation checks (TFU vs TFI, CRITFREQ), cut-off and mesh rules, interaction-node count and method choice, CHECK/FIXEDINT/EXCSTRCHK/HINGED, how SASSI-EDU is verified (Verification Manual), incoherency, 2D, SYMM, what SASSI-EDU does not do |
| 11 | `11-embedded-building` | Advanced | An embedded building: the subtraction methods against FV | ex08_embedded_building | a shear-wall building with a two-level basement, a grid of interior walls and a tower (walls, slabs, roofs, equipment masses, concrete properties), the excavated soil and the interior structure on separate nodes (manual rule 11, EXCSTRCHK, RMVUNUSED/NCOM), the element plot with the soil hidden, shown again with the interaction nodes, and a cutaway (WINDOWSETTINGS), the backfill question, INTGEN sets and their counts, memory and run time; FV, FI-FSIN (SM) and FI-EVBN (MSM) runs; transfer functions, ISRS and basement-wall forces against FV; the spurious resonance of SM animated (HARMFRAME) and explained by the enclosed-soil frequency; CRITFREQ and what a validation must document; FFV and when each method is acceptable (stiff and soft sites) |
