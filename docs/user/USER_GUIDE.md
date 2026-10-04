# SASSI-EDU User Guide

SASSI-EDU is an educational Python re-implementation of the **SASSI flexible-volume method** of
soil-structure interaction (SSI) analysis (Lysmer et al., UC Berkeley, 1981), with the modules,
commands, option codes, file names and messages of the **ACS SASSI Version 3** user manual. This
guide is written for a structural engineer who builds models in ANSYS and is learning SSI. It
explains how to use the program, why each step exists and what can go wrong.

> **Disclaimer.** SASSI-EDU is a learning tool. It is not affiliated with, endorsed by, or a
> substitute for ACS SASSI or its vendor, and it is not qualified (10 CFR 50 Appendix B, ASME
> NQA-1) for licensing or design work. Its results are verified against published solutions
> (see the [Verification Manual](../verification/VERIFICATION_MANUAL.md)), but you are responsible
> for any use you make of them.

**Companion documents.** The [Theory Manual](../theory/THEORY_MANUAL.md) derives the equations;
the [Verification Manual](../verification/VERIFICATION_MANUAL.md) shows computed versus reference
values for every verification problem; the [Command Reference](../reference/COMMAND_REFERENCE.md)
lists every command. Five guides cover parts of the program in depth: [GUI.md](GUI.md) the browser
interface, [ANSYS.md](ANSYS.md) the ANSYS model converters, [OPTION_A.md](OPTION_A.md) SSI loads for a
second-step ANSYS analysis (module LOADGEN), [OPTION_NON.md](OPTION_NON.md) nonlinear structures
(module NONLINEAR) and [WATER.md](WATER.md) water in pools and tanks. The
[documentation index](../index.md) links everything.

## Contents

1. [Getting started](#1-getting-started)
2. [How an SSI analysis works](#2-how-an-ssi-analysis-works)
3. [Units, axes and sign conventions](#3-units-axes-and-sign-conventions)
4. [The command language](#4-the-command-language)
5. [Building the structural model](#5-building-the-structural-model)
6. [ANSYS models and the two-step approach](#6-ansys-models-and-the-two-step-approach)
7. [Site and soil](#7-site-and-soil)
8. [Structure, excavation and interaction nodes](#8-structure-excavation-and-interaction-nodes)
9. [Frequencies and the FFT](#9-frequencies-and-the-fft)
10. [Running ANALYS](#10-running-analys)
11. [Post-processing](#11-post-processing)
12. [Advanced options](#12-advanced-options): [incoherency](#121-incoherency-wave-passage-and-multiple-excitation), [2D](#122-two-dimensional-plane-strain-ssi), [symmetry](#123-symmetry-planes-half-and-quarter-models), [nonlinear soil and SOIL-NON](#124-nonlinear-soil-soil-near-field-iterations-and-soil-non), [Option NON](#125-nonlinear-structures-option-non), [water](#126-water-in-pools-and-tanks), [Option A](#127-option-a-and-option-aa), [TSHELL and the node optimizer](#128-thick-shells-tshell-and-the-node-optimizer), [AFWRBAT](#129-frequency-split-runs-afwrbat), [not available](#1210-what-is-not-available)
13. [Common errors and how to fix them](#13-common-errors-and-how-to-fix-them)
14. [Engineering guidance](#14-engineering-guidance)
15. [The tutorial examples](#15-the-tutorial-examples)
16. [Verification](#16-verification)
17. [ACS SASSI terms and their SASSI-EDU equivalents](#17-acs-sassi-terms-and-their-sassi-edu-equivalents)
18. [Limitations](#18-limitations)

---

## 1. Getting started

### 1.1 Installation

SASSI-EDU is pure Python (no compiled extensions). It needs Python 3.9 or newer with NumPy, SciPy,
matplotlib and plotly. From the project root:

```bash
python3 -m venv .venv
.venv/bin/pip install -e .[test]        # installs the 'sassi' and 'sassi-gui' commands
.venv/bin/python -m pytest -q           # optional: unit tests and verification problems
```

On Windows replace `.venv/bin/` by `.venv\Scripts\`. Everything runs locally; the GUI serves its
page on `127.0.0.1` only and needs no internet connection.

### 1.2 Three ways to drive the program

All three drive **one command interpreter**, so a model gives the same results whichever you use
(requirement UI-01, rule L17).

| Way | Command | Use it for |
|---|---|---|
| Console | `.venv/bin/sassi` | typing commands one by one, as in the ACS SASSI Command Entry line |
| Batch | `.venv/bin/sassi run model.pre` (`--quiet`, `--cwd DIR`) or `sassi -c "LINE" -c "LINE"` | running a complete `.pre` command file; exit status 0 = no error |
| Browser GUI | `.venv/bin/sassi-gui` (or `python -m sassi.ui.server`) | menus, Options dialogs, 3D model plots, result plots |

In the console, `exit` or Ctrl-D quits; the Up arrow recalls earlier lines. In the GUI, every menu
item, dialog OK and toolbar button is turned into command text and shown in the Command History,
so a GUI session can be saved and replayed as a `.pre` file. See [GUI.md](GUI.md) for the menus
(Model, File, Plot, Modules, Options, View, Help) and the Options > Analysis dialog, whose tabs
(EQUAKE, SOIL, SITE, POINT, HOUSE, FORCE, ANALYS, MOTION, STRESS, RELDISP, NONLINEAR, AFWRITE)
correspond one to one to the option commands described in this guide.

### 1.3 A first analysis in five minutes

Run the first tutorial example from the `examples` directory:

```bash
cd examples
../.venv/bin/sassi run ex01_surface_stick.pre
```

(or, from the project root, `.venv/bin/sassi --cwd examples run ex01_surface_stick.pre`). The file
builds a four-storey stick on a stiff surface mat, checks it, writes the module input decks and runs
SITE, POINT, HOUSE, ANALYS, MOTION, STRESS and RELDISP in about four seconds. The results are in
`examples/ex01/`:

* `00085TR_X.TFU` / `.TFI`: the roof transfer function at the computed and at the interpolated
  frequencies (it is 1.0 at low frequency and peaks near 3.5 Hz with SSI);
* `00085TR_X01.RS`: the roof in-structure response spectrum (ISRS) for the first damping value;
* `ex01_MOTION.out`, `ex01_STRESS.out`: listings with the maxima.

Open the file `ex01_surface_stick.pre` in a text editor while you read section 15.1: every line is
commented.

### 1.4 Files of a model

A model has a **name** and a **directory** (`MDL,<name>,<path>`). Every module reads and writes its
files in that directory:

| Kind | Examples | Written by |
|---|---|---|
| command file | `ex01.pre` | you, or `WRITE` |
| module input decks | `ex01.sit`, `.poi`, `.hou`, `.anl`, `.mot`, `.str`, `.rdi`, `.soi`, `.equ`, `.frc`, `.eql` (NONLINEAR) | `AFWRITE` |
| side files of the decks | `.pin` (nonlinear soil groups), `.nls` (SOIL-NON layers), `HOUSE.opt` (EDUOPT switches of HOUSE) | `AFWRITE` |
| LOADGEN deck | `ex01.lgn` | `RUNLOADGEN` |
| listings | `ex01_SITE.out`, `ex01_ANALYS.out`, ... | each module |
| binary inter-module files | `FILE1`, `FILE2`, `FILE3`, `ex01.N4` (FILE4), `COOSK`, `COOSM`, `FILE77` (incoherency factors), `FILE8`, `FILE9` | SITE, POINT, HOUSE, ANALYS, FORCE |
| iteration files (text) | `FILE78`, `FILE74` (nonlinear soil); `PANEL.NON`, `Panel_EQL_Matl_Prop.txt`, `ex07_new.hou` (Option NON) | HOUSE, STRESS; NONLINEAR |
| results | `nnnnnTR_X.TFU/.TFI/.ACC`, `nnnnnTR_X01.RS`, `.THD`, `BEAMS_002_00001_MZI.THS` | MOTION, RELDISP, STRESS |
| ANSYS loads (APDL) | `ex01_LGS.inp` (static), `ex01_LGD.inp` (dynamic) | LOADGEN |
| check report | `ex01.err` | `CHECK` |
| saved database | `ex01.sdb` | `SAVE` |

The binary files keep the legacy names of ACS SASSI (FILE1, FILE8 ...) but are NumPy `.npz`
containers without an extension: they can be opened in Python with
`sassi.io.container.read_container(path)`.

---

## 2. How an SSI analysis works

### 2.1 What problem SASSI solves

In ANSYS you would analyse soil-structure interaction with a large box of soil elements around the
structure, absorbing boundaries on its sides and the earthquake applied at its base, in the time
domain. SASSI does something different, and much cheaper:

* the analysis is in the **frequency domain**: the structure and soil are linear (equivalent-linear),
  damping is hysteretic (complex moduli), and the system is solved at a few tens to a few hundred
  frequencies; time histories are obtained afterwards by interpolation and inverse FFT;
* the **soil is not meshed** outside the structure. The horizontally layered site (the *free field*)
  is treated semi-analytically by the thin-layer method: its dynamic stiffness at the nodes where
  the structure touches the soil (*interaction nodes*) is computed exactly for a layered medium
  that extends to infinity, including radiation damping;
* the soil that the basement replaces is handled by **substructuring**: the finite-element model
  contains the structure *plus* the excavated soil, the excavated soil is *subtracted*, and the
  free-field impedance is *added* at the interaction nodes. This is the **flexible-volume method**:

```
   SSI system  =  free field (layered soil, no excavation)
               +  structure (FE model, basement included)
               -  excavated soil (FE model of the removed soil, free-field properties)
```

The equation solved at each frequency ω is (manual Eq. 2.1; [Theory Manual §3](../theory/THEORY_MANUAL.md#3-flexible-volume-substructuring))

```
[ (K*s - ω² Ms) - (K*e - ω² Me) + Xff ] U = Xff U'f
```

where `K*s, Ms` are the structure matrices (HOUSE), `K*e, Me` the excavated-soil matrices (HOUSE),
`Xff` the free-field impedance at the interaction degrees of freedom (POINT + ANALYS) and `U'f` the
free-field motion at the interaction nodes (SITE). `U` are total displacements per unit control
motion: the **transfer functions** stored in FILE8.

### 2.2 The module chain (manual Fig. 1.1)

```
 INPUT                    ANALYSIS                                RESULTS
 FORCE --(FILE9)------------\
 EQUAKE                      \
   | control motion           \
 SOIL --(FILE73)----------------\-------------------------------\
   | FILE88 (strain-compatible) \                                \
 SITE --(FILE1)---------------> ANALYS --(FILE8)--> COMBIN --+--> MOTION --> .TFU .TFI .ACC .RS
   | FILE2                      /                            |      | complex .TFI
 POINT --(FILE3)-------------> /                             |    RELDISP --> .TFD .THD
 HOUSE --(FILE4, COOSK, COOSM,/                              +--> STRESS --> .TFU .TFI .THS
   ^      FILE77 incoherency)
   ^-------------------------(FILE74, near-field soil iterations)--- STRESS
   ^-------------------------(<model>_new.hou, Option NON)---------- NONLINEAR <-- .THD (X+Y+Z)
                       MOTION .ACC + RELDISP .THD (or FILE8) --> LOADGEN --> APDL for ANSYS (Option A)
```

| Module | What it computes | Reads | Writes |
|---|---|---|---|
| EQUAKE | spectrum-compatible acceleration histories, SRP 3.7.1 checks | `.equ`, target spectrum `.rsi` | `.acc .vel .dis .rso .psd .fft` |
| SOIL | 1D site response: equivalent linear (SHAKE) or nonlinear in the time domain (SOIL-NON) | `.soi`, [`.nls`], control motion | `ACCxxx.TH`, spectra, FILE73, FILE88 |
| SITE | Mode 1: Rayleigh/Love modes of the layered site; Mode 2: free-field motion for unit control motion | `.sit` | FILE2, FILE1 |
| POINT | point-load solutions of the layered site (POINT3 axisymmetric central zone; POINT2 plane strain) | `.poi`, FILE2 | FILE3 |
| HOUSE | stiffness and mass of the structure and of the excavated soil; incoherency factors; near-field soil properties | `.hou`, `.sit`, [`.pin`, FILE74] | `<model>.N4` (FILE4), COOSK, COOSM, DOFSMAP, FILE90, FILE91, [FILE77, FILE78] |
| FORCE | load vectors of external harmonic loads | `.frc` | FILE9 |
| ANALYS | impedance `Xff = Fff⁻¹`, flexible-volume equation, transfer functions | `.anl`, FILE1 or FILE9, FILE3, FILE4, [FILE77] | FILE8, restart files, FOUN* impedance tables |
| COMBIN | merge two FILE8 computed for different frequency sets | FILE81, FILE82 | FILE8 |
| MOTION | interpolation, convolution with the control motion, ISRS | `.mot`, FILE8, control motion | `.TFU .TFI .ACC .RS` |
| RELDISP | relative displacements from complex `.TFI` | `.rdi`, `.TFI` | `.TFD .THD` |
| STRESS | element stress, strain and force transfer functions and histories | `.str`, FILE4, FILE8 | `.TFU .TFI .THS`, maxima, [FILE74] |
| NONLINEAR | Option NON: hysteresis of wall panels and springs, equivalent-linear properties | `.eql`, `.hou`, RELDISP `.THD` | `<model>_new.hou`, `Panel*` / `SPRING*` files, convergence history |
| LOADGEN | Option A: SSI motions as ANSYS loads (APDL) | `.lgn`, FILE4, COOSM, MOTION `.ACC`, RELDISP `.THD` (or FILE8) | `<model>_LGS.inp`, `<model>_LGD.inp` |

Brackets mark the files of the optional analyses (incoherency, nonlinear soil).

Each module is an independent program. The interpreter writes its input deck with `AFWRITE`, and
`RUN<MODULE>` (for example `RUNSITE`) runs it in the model directory. A module can also be run in
batch exactly as the original executables: `python -m sassi.modules.site < SITE.inp`, where
`SITE.inp` holds three lines (model name, deck file, listing file).

### 2.3 The standard seismic workflow

Paraphrasing the 13 steps of manual §4.1.1:

1. **Control motion.** Generate spectrum-compatible acceleration histories (EQUAKE) or use given ones.
2. **Know your structure.** Run a fixed-base modal analysis (in ANSYS or with HOUSE) to find the
   important frequencies.
3. **Cut-off frequency.** Choose the highest SSI frequency from the input content and the structure
   (typically 30-40 Hz for soil sites, 60-70 Hz for rock sites); the soil layers and the excavation
   mesh must pass it (section 7.2).
4. **Free-field soil properties.** Run SOIL to obtain strain-compatible shear-wave velocities and
   damping for the input motion.
5. **Site model.** Use the SOIL properties in SITE (`SITEX,1`) and for the excavated soil layers.
6. **Structure model.** Build the FE model (it may include near-field soil or backfill).
7. **Excavation model.** Build the excavated soil for an embedded structure.
8. **Method.** Choose the interaction-node set: FV, FI or FFV (section 8.3).
9. **SSI frequencies.** Choose 40-80 frequencies for stick models, 100-300 for large FE models.
10. **Initiation run.** SITE (Mode 1 and 2), POINT, HOUSE, ANALYS.
11. **Post-processing.** MOTION, RELDISP, STRESS; compare computed (`.TFU`) and interpolated (`.TFI`)
    transfer functions.
12. **Refine.** Add frequencies where the interpolation is not supported by computed points
    (CRITFREQ), run ANALYS for them and merge with COMBIN.
13. **Restarts.** Re-analyse a changed structure or a new seismic environment with the ANALYS
    restart modes (section 10.3).

For **forced vibration** (machine foundations, impedance studies) the control motion is replaced by
a reference load history, SOIL is not used, SITE runs Mode 1 only, and FORCE produces FILE9, which
replaces FILE1 in ANALYS (`ANALYS <type>` = 1).

### 2.4 What runs in what order

| Situation | Modules |
|---|---|
| Seismic initiation | [EQUAKE] → [SOIL] → SITE → POINT → HOUSE → ANALYS → [COMBIN] → MOTION → RELDISP; STRESS |
| X, Y and Z input in one run | SITE ×3 (copy FILE1 to FILE1X/Y/Z with FCOPY) → POINT → HOUSE → ANALYS (`<simul>` = 1) → MOTION per direction |
| Forced vibration | SITE (Mode 1) → POINT → FORCE (copy FILE9 to FILE9001 ...) → HOUSE → ANALYS (`<type>` = 1) → MOTION / STRESS |
| New time history only | MOTION (STRESS, RELDISP) |
| New structure, same soil and interaction nodes | HOUSE → ANALYS restart `<mode>` = 1 → MOTION / STRESS |
| New seismic environment | SITE → ANALYS restart `<mode>` = 2 |
| More frequencies | ANALYS on the new frequencies → FMOVE old/new FILE8 to FILE81/FILE82 → COMBIN → MOTION |
| Frequency-split (parallel) run | `AFWRBAT,<n>` → `run_k` in each of the n folders (SITE POINT HOUSE ANALYS) → `combine` (COMBIN, then MOTION / STRESS / RELDISP) |
| Incoherent stochastic simulation | SITE ×3 → POINT → HOUSE (FILE77001 ... FILE77Ns) → ANALYS (`<simul>` = Ns) → MOTION per sample → average of the spectra |
| Near-field soil nonlinearity | SOIL → SITE ×3 → POINT → HOUSE → ANALYS (`<save>` = 1) → STRESS ×3 → COMBXYZSTRAIN, then repeat HOUSE → ANALYS `<mode>` = 1 → STRESS ×3 → COMBXYZSTRAIN (`NLSSIITER`) |
| Nonlinear structure (Option NON) | elastic SSI with restart files → MOTION + RELDISP ×3 → COMBXYZTHD → NONLINEAR, then repeat HOUSE (`<model>_new.hou`) → ANALYS `<mode>` = 1 → MOTION + RELDISP ×3 → COMBXYZTHD → NONLINEAR (`NONLINITER`) |
| ANSYS second step (Option A) | a finished SSI run with MOTION and RELDISP → LOADGEN (`RUNLOADGEN,STATIC` or `DYNAMIC`) → ANSYS |

MOTION and STRESS are independent of each other; RELDISP needs MOTION's complex `.TFI` files.

---

## 3. Units, axes and sign conventions

* **Units.** Any consistent set. Weights are specific weights (force/length³) and masses are entered
  as weights or masses (`MUNITS`); the program divides by the **gravity** you give
  (`GRAVITY,9.81` for m, kN, t, s; `GRAVITY,32.2` for ft, kip, s). New models start with the
  manual's default 32.2: set `GRAVITY` (or HOUSE argument 1) first when you work in SI. SOIL has its
  own gravity (SOIL argument 2). Control motions are always in **g**.
* **Axes.** Global right-handed Cartesian, **Z up**. Two-dimensional models lie in the X-Z plane.
  Depths are measured down from the ground elevation `GROUNDELEV` (HOUSE `<gelev>`).
* **Damping.** Hysteretic (frequency independent), applied as a complex modulus
  `G* = G(1 - 2β² + 2iβ√(1-β²))` to every material, layer and spring (`CMODFORM,1` selects
  `G(1 + 2iβ)` for benchmarking only). An ANSYS model with structural damping `g` reproduces it
  with `E_ANSYS = E(1-2β²)` and `g = 2β√(1-β²)/(1-2β²)` ([ANSYS.md](ANSYS.md) section 4).
* **Time dependence.** `exp(+iωt)`; transfer functions follow the FFT convention of `numpy.fft`.
* **Transfer functions.** A seismic FILE8 holds `H = U/U_cp`, the ratio of the total motion of a
  node to the control motion. It is dimensionless and the same for displacement and acceleration.
* **Element results.** BEAMS end forces are the forces exerted *on* the element at I and J in local
  axes; SHELL membrane outputs are stresses and bending outputs moments per unit length in the local
  shell axes; SOLID stresses are centroid stresses in global axes.

---

## 4. The command language

### 4.1 Syntax

A command line is

```
Keyword, p1, p2, ..., pn
```

| Rule | Example |
|---|---|
| Fields are separated by commas; the first comma may be a blank | `N,1,0,0,0` or `N 1,0,0,0` |
| Names are case-insensitive; file names keep their case | `ngen,3,5,1,5,1,0,0,10` |
| A blank or omitted field takes the documented default | `MOTION,0,0,0,20` |
| A line starting with `*` is a comment; there are no inline comments | `* soil profile` |
| Numbers are free format, Fortran exponents allowed | `M,1,3.0E7,0.2,24.0,0.05,0.05,1`, `1.0D9` |
| A text argument in the last position takes the rest of the line, commas included | `TIT,Stick model, 4 storeys` |
| `"..."` embeds commas in any field (SASSI-EDU extension) | `VAR,CMD,"EDUOPT,TFFILE,FILE8X"` |
| Node and element lists accept ranges | `NOUT,1,1,1,0,0,1,1,41,82-85` |
| Only the full name or the documented 4-character abbreviation is accepted | `AFWR` = `AFWRITE`; `AFW` is `Command not found` |

An unknown command prints `<X> Command not found` and processing continues; a failing command prints
its error and the next line still runs. The [Command Reference](../reference/COMMAND_REFERENCE.md)
lists every command with its syntax, abbreviation, tier and implementation status. The few manual
commands this build does not implement (the binary result databases, GETENV/SETENV, FREAD/MREAD and a
handful of others, section 18) are still accepted, so `.pre` files of the original program load
without loss: option commands are stored and written back by `WRITE`, and every such command prints
`<CMD> is not available in this build (tier Pn)` or the manual's own message.

### 4.2 Models in memory

The interpreter can hold several numbered models (manual §5.2). The active model receives the
commands.

```
* name and directory of the active model (also the working directory)
MDL,ex02,ex02
* copy the active model to model 2
CPMODEL,2
* make model 2 active
ACTM,2
* give the copy its own name and directory
MDL,ex02fsin,../ex02_fsin
MODELLIST                     * list the models in memory
```

### 4.3 Variables, loops and macros

Variables hold lists; `FOREACH` repeats a command for each item (`#` is the 1-based loop index):

```
VAR,FLOOR,82,83,84,85
* 9810 kN = 1000 t on each floor
FOREACH,FLOOR,MT,@FLOOR[#],9810,9810,9810
VAR,LFIRST,1,26,51,76,101,126
VAR,LLAST,25,50,75,100,125,150
* interaction nodes level by level
FOREACH,LFIRST,INT,@LFIRST[#],@LLAST[#],1,1
```

`@X` is the counter of variable X, `@X++` increments it, `@X[i]` is item i. `LOADMACRO,<name>,<file>`
loads a macro file whose `$1$`, `$2$` ... placeholders `MACRO,<name>,<a1>,<a2>,...` replaces.
`RND`, `ADDRND`, `RNDSEED` and `REDUCESET` build random or sorted lists (useful for parametric and
sensitivity studies).

### 4.4 Saving and reading models

| Command | Effect |
|---|---|
| `WRITE,[file]` | writes the complete model, including every analysis option, as commands (`<model>.pre`); `INP` of that file rebuilds an identical model |
| `INP,<file>` | executes a `.pre` file; nested INP is allowed |
| `SAVE` / `RESUME` | binary snapshot `<model>.sdb` of the active model (needs `MDL`) |
| `STATUS`, `NLIST`, `ELIST`, `MLIST`, `LLIST`, `LFREQ` ... | list model data |

`WRITE` is the best way to keep a model: the `.pre` file is readable, diffable and can be edited.

### 4.5 CHECK, AFWRITE and RUN

```
* modules to check and write: SITE POINT HOUSE ANALYS MOTION STRESS RELDISP
AOPT,0,0,0,1,1,1,0,0,1,0,1,1,1,0
CHECK                              * errors and warnings per module -> <model>.err
AFWRITE                            * CHECK, then the deck of every enabled module without errors
RUNSITE
RUNPOINT
...
```

`AOPT` has one flag per module in the order EQUAKE, SOIL, (DEP), SITE, POINT, HOUSE, (DEP), FORCE,
ANALYS, COMBIN, MOTION, STRESS, RELDISP, NONLINEAR; the two DEP slots must be 0. LOADGEN (Option A)
has no AOPT flag: `RUNLOADGEN` writes its deck `<model>.lgn` and runs it. A module with any
CHECK **error** gets no deck and cannot be run; **warnings** never block. Run `AFWRITE` again after
every change of the model or of the options: the decks are snapshots. `RUN<MODULE>` checks that the
module's input files exist (for example ANALYS needs FILE1, FILE3 and FILE4) and names the module
that produces a missing file.

`FCOPY,<src>,<dst>` and `FMOVE,<src>,<dst>` copy or rename files in the model directory (FILE1 to
FILE1X, FILE8 to FILE81 ...), so that multi-run workflows can be scripted completely.

---

## 5. Building the structural model

### 5.1 Nodes and coordinate systems

```
* node 1
N,1,-10,-10,0
* node 9
N,9,10,-10,0
* nodes 2..8 equally spaced between 1 and 9
FILL,1,9
* 8 copies of nodes 1..9, numbers +9, shifted 2.5 in y each time
NGEN,8,9,1,9,1,0,2.5,0
```

`NGEN,<itim>,<step>,<n1>,<n2>,<inc>,<dx>,<dy>,<dz>` is the workhorse for regular meshes. `NMED`,
`LMOVE`, `NMOVE` and `NSCALE` create or move nodes; `NDEL` deletes them. Local Cartesian systems are
defined by origin and angles (`LOC`) or by three nodes (`LOCAL`) and activated with `CSYS`; nodes keep
the system they were defined in, and `GLOBAL` converts them.

### 5.2 Groups and element types

Elements live in **groups**; a group has one element type. `GROUP,<ng>,<type>` creates or activates
a group, then `E,<ne>,<n1>,...` defines its elements and `EGEN` copies element patterns.

| Type | Name | Nodes | DOF/node | Formulation | Mass |
|---|---|---|---|---|---|
| 1 | SOLID | 8 (prism/pyramid by repeated nodes) | 3 | trilinear hexahedron; 9 incompatible modes for structural solids (`MOPT`) | ½ lumped + ½ consistent |
| 2 | BEAMS | I, J, K (K orients the section, no DOF) | 6 | 3D Timoshenko frame, shear areas, end releases `KI`/`KJ` | consistent |
| 3 | SHELL | I, J, K, L (L omitted or repeated = triangle) | 6 | flat thin shell: plane-stress membrane + DKQ/DKT Kirchhoff plate; **no drilling stiffness** | lumped |
| 4 | PLANE | I, J, K, L in the X-Z plane | 2 | plane strain, unit thickness | ½ + ½ |
| 5 | TSHELL | I, J, K, L | 6 | thick (Mindlin-Reissner) shell, MITC4 transverse shear (no shear locking), small automatic drilling stiffness (section 12.8) | lumped, with rotary inertia |
| 7 | SPRING | I, J | 6 | six uncoupled global springs `SC` with one damping ratio | none |
| 9 | GENERAL | I, J (global) or I, J, K (local axes) | 6 | user 12×12 complex stiffness `MXR + i MXI` and mass `MXM` | user |

Per-element attributes are set by element ranges: `ETYPE` (structure or excavated soil, section 8.1),
`EINT` (SOLID integration order 0/1/2; TSHELL 0 reduced / 1 selective), `THICK` (shell thickness),
`MSET` (material, or soil layer for excavated solids), `RSET` (beam section, spring or matrix
property) and `KI`/`KJ` (beam end releases P1 P2 P3 M1 M2 M3). `MACT` and `RACT` set the material
and property index of the elements created next.

**Beam local axes** (like ANSYS with a K node): axis 1 runs from I to J, axis 2 points toward K in
the plane I-J-K, axis 3 = 1 × 2. Bending in the 1-2 plane uses `I3` and `As2`; bending in the 1-3
plane uses `I2` and `As3`. **Shell local axes**: x' joins the mid-points of sides L-I and J-K, z'
is the normal (I → J → K → L counter-clockwise seen from +z'), y' = z' × x'.

### 5.3 Materials, sections, springs and matrices

```
* type 1 (E, nu), 2 (M, G), 3 (Vp, Vs)
M,<nm>,<val1>,<val2>,<weight>,<pdamp>,<sdamp>,<type>
* concrete, 5 % damping, 24 kN/m3
M,1,3.0E7,0.2,24.0,0.05,0.05,1
* beam section (As = 0: no shear deformation)
R,<nm>,<A>,<As2>,<As3>,<J>,<I2>,<I3>
* spring constants in global axes
SC,<nm>,<kx>,<ky>,<kz>,<kxx>,<kyy>,<kzz>,<damp>
* GENERAL matrices, upper-triangle rows
MXR,<p>,<row>,<t1>,...   MXI,...   MXM,...
```

The specific weight gives the mass density `ρ = weight / gravity`. A material with weight 0 is
massless (useful for sticks whose masses are lumped at the floors). Damping `pdamp` applies to the
constrained (P-wave) modulus and `sdamp` to the shear modulus; for beams and shells use equal values.
`MLIST`, `RLIST`, `SCLIST` and `MXLIST` print the tables with derived quantities.

### 5.4 Masses and loads

```
* nodes 82..85: masses given as weights (default 1 = weight)
MUNITS,82,85,1,1
* translational mass (weight units -> divided by gravity)
MT,82,9810,9810,9810
* rotational inertia about X, Y, Z
MR,30,1.2E6,1.2E6,2.4E6
* force factor 1 in X at node 25 (FORCE, vibration analysis)
F,25,1,0,0
* moment factor about Z
MM,25,0,0,1
```

`F` and `MM` take factors (and arrival times in arguments 5-7) that multiply one reference load
history; a moving load is a set of loads with increasing arrival times. `MOPT,<incomp>,<matrix>,<mass>,<force>`
controls the incompatible modes, the GENERAL mass units and whether repeated mass/force definitions
add or overwrite.

### 5.5 Boundary conditions and unstiffened rotations

`D,<n1>,<n2>,<inc>,<val>,<labels>` fixes (1) or frees (0) degrees of freedom (labels UX UY UZ ROTX
ROTY ROTZ, DISP, ROT, ALL). The active DOFs of a node are those of the elements attached to it:
the rotations of a node connected only to solids are removed automatically. A DOF that an element
*defines* but does not *stiffen* must be restrained by you, because the system would be singular:

* **shell drilling rotations** (rotation about the shell normal): `FIXSHLROT,[k]` adds soft springs,
  `FIXROT` fixes axis-parallel shells with `D` and puts springs on oblique ones;
* rotations of nodes shared by solids and springs: `FIXSPRROT`;
* rotations of solid-only nodes: `FIXSLDROT`.

CHECK warns about unrestrained drilling rotations (EDU-06). Do not fix a rotation that a beam or the
soil should resist: in example 1 the stick base keeps its torsion free so that the torsional soil
impedance acts on the mat.

### 5.6 Model checks

| Command | Reports |
|---|---|
| `EXCSTRCHK` | interior excavation nodes shared with structural elements (a serious modelling error) |
| `FIXEDINT` | interaction nodes with a fixed translation |
| `HINGED` | 6-DOF elements meeting solids at a single node (unintended hinges), drilling joints |
| `KINT` | beam orientation nodes that are interaction nodes |
| `FREESPRING` | unrestrained nodes connected only to springs |
| `INTCOUNT` | number of interaction nodes and the ANALYS memory estimate |
| `USED` | fixes all DOFs of unused nodes (modifies the model) |

Use the GUI (Plot > Model, Plot > Nodes) or `MODELPLOT` to look at the mesh, the interaction nodes
and the fixities before you run anything.

---

## 6. ANSYS models and the two-step approach

There are two ways to combine SASSI-EDU with ANSYS: move the *model* between the programs (section
6.1), and use the SSI *results* as loads of a refined ANSYS model, the two-step approach of ACS SASSI
Option A (section 6.2).

### 6.1 Moving models between ANSYS and SASSI-EDU

Most engineers build the structure in ANSYS. SASSI-EDU reads an ANSYS `CDWRITE` file and writes an
APDL file back:

```
* ANSYS -> SASSI-EDU (gravity in model units is required)
CONVERT,ANSYS,,mybuilding.cdb,9.81
ANSYS                                 * SASSI-EDU -> <model>.inp (APDL), read with /INPUT in ANSYS
```

The converter maps the common ANSYS elements (SOLID45/65/185/186/187/95, SHELL63/181/281,
BEAM4/44/188/189, PIPE288, LINK8/180, COMBIN14, MASS21, MATRIX27, PLANE42/182), the materials, real
constants and sections; it writes an element map `<model>_cdb.map` and warns about everything it
cannot map (for example, SHELL181 becomes a thin Kirchhoff SHELL and mid-side nodes are dropped). After a conversion you still have to: set `FIXROT`, the ground elevation,
`ETYPE`/`ETYPEGEN` and the interaction nodes, and add the excavated soil. Compare the fixed-base
frequencies of the two programs before the SSI run (verification problem VP-A1 does this for a
BEAM188 cantilever). The complete description, including the damping mapping and the beam-axis
rules, is in **[ANSYS.md](ANSYS.md)**.

### 6.2 The two-step approach: SSI motions as ANSYS loads (Option A, LOADGEN)

Step 1 is the SSI analysis in SASSI-EDU; step 2 is a refined ANSYS model of the structure (finer mesh,
other element types, nonlinear materials, contact, a refined soil model for wall pressures) loaded by
the motion that step 1 computed. Module **LOADGEN** writes the step-2 loads as APDL. The decomposition
behind it is `u = u_r + ι u_g`: a rigid-body translation of the ground motion strains nothing, so the
structure with its interface displacements prescribed obeys
`M_ss ü_r + C u̇_r + K_ss u_r = -M_ss ι a_g(t) - K_sb u_r,b(t)` ([Theory Manual §20](../theory/THEORY_MANUAL.md#20-option-a-ssi-loads-for-a-second-step-ansys-analysis)).

* **Equivalent static loads** (`LOADGEN`, "ANSYS Eq. Static Load"): at critical times t* (by default the
  largest base shear) the nodal inertia forces `F = -m a(t*)` of the absolute accelerations (MOTION
  `.ACC`) and, with "Disp. and Accel.", the interface displacements relative to the free field
  (RELDISP `.THD`); with "Acceleration" the interface is fixed (the classical fixed-base equivalent
  static analysis). The masses are generated from the HOUSE mass matrix or read from your own file.
* **Dynamic loads** (`LOADGENDYN`, "ANSYS Dynamic Load"): TABLE arrays of the ground acceleration
  (`ACEL`, the acceleration of the reference frame) and of the interface displacements relative to it,
  Rayleigh damping and the transient solution commands.

A finished SSI run with MOTION (complex TFs saved) and RELDISP (free-field reference) on the mass and
interface nodes is all LOADGEN needs (or FILE8 alone with `<source>` = FILE8):

```
* Disp. and Accel.; lumped masses generated from COOSM
LOADGEN,DISPACC,0,0,0,LUMPED,1
* one critical time: the largest base shear
LGTIME,V,1
* <model>.lgn -> <model>_LGS.inp
RUNLOADGEN,STATIC
* Rayleigh damping 5 % at 1 and 20 Hz
LOADGENDYN,,,REL,0,RESULTS,0,0,0.05,1,20
* -> <model>_LGD.inp
RUNLOADGEN,DYNAMIC
ANSYS                                      * the structural model as APDL (<model>.inp)
```

In ANSYS: `/INPUT,<model>,inp`, then `/INPUT,<model>_LGS,inp` (or `_LGD`); after the static step the
total base reaction equals the SSI base shear (VP-LA1 checks the inertia forces against the STRESS base
shear and moment to 0.04 %). Option A is exact for a linear structure with the same mass, stiffness and
damping as the SSI model (VP-LA2 replays the exported second step and recovers the SSI accelerations to
5e-9) and approximate otherwise: ANSYS Rayleigh damping is exact at two frequencies only, the static step
is one instant without damping forces, and anything you add in ANSYS must not change the motion of the
foundation-soil interface. The commands, the files, the node numbering (`LGMAP`) and a step-by-step
example are in **[OPTION_A.md](OPTION_A.md)**. Option AA (SSI with ANSYS matrices) is not available.

For a **nonlinear** structure, the converged cracked moduli and spring constants of an Option NON
analysis (section 12.5) are the secant properties to use in an ANSYS model of the same building.

---

## 7. Site and soil

### 7.1 Soil layers and the layer list

```
* one soil layer type
L,<nm>,<thick>,<weight>,<Vp>,<Vs>,<pdamp>,<sdamp>
L,1,0.5,19.0,600,300,0.05,0.05
L,2,1.0,20.0,1000,500,0.04,0.04
* half-space properties (thickness not used)
L,3,1.0,21.0,2000,1000,0.02,0.02
* the site, top layer first (append)
TOPL,1,1,1,1,1,1,1,1,1,1
TOPL,2,2,2,2,2,2,2,2,2,2
TOPL,2,2
```

`L` defines layer *types*; `TOPL` lists the layers of the site from the ground surface down (it
appends; `TOPL,0` clears). The **interfaces** are the tops of the TOPL layers: interface 1 is the
ground surface, interface `i` the top of TOPL layer `i`, and interface `nTOPL+1` the top of the
half-space. `LAYERPLOT` draws the column.

### 7.2 Layering rules

* **1/5-wavelength rule.** Displacements vary linearly inside a thin layer, so each layer must be
  thin compared with the shortest wavelength: `h ≤ Vs/(5 f_cut)`, i.e. the layer's *passing
  frequency* `f_pass = Vs/(5h)` must exceed the cut-off frequency. Use the strain-compatible Vs
  (it is lower). CHECK warns (EDU-10) when a layer does not pass the cut-off.
* **Enough layers.** The manual recommends more than 20 top layers; too few layers degrade the
  Rayleigh and Love modes on which POINT is built.
* **Interaction nodes on interfaces.** Every interaction node must lie on a layer interface at or
  below the ground surface (EDU-01). For an embedded structure the excavation mesh levels must
  therefore coincide with TOPL interfaces, and each embedment layer should have its own L number
  (EDU-28), so that excavated elements can refer to the layer they replace.
* **Poisson's ratio.** Avoid ν > 0.47 (saturated soils with high Vp): the thin-layer solution can
  become unstable at isolated frequencies (EDU-11). Look at the transfer functions, and remove such
  frequencies (`REMOVEFREQ`) if needed.

### 7.3 Half-space or rigid base

SITE `<nl>` sets how the bottom of the profile is modelled:

* `<nl>` = 0: **rigid base** at the bottom of the last TOPL layer (bedrock much stiffer than the soil);
* `<nl>` = 4 ... 20 (10-20 recommended, 20 for best accuracy): a **visco-elastic half-space** with
  the properties of layer `<hs>`. SITE adds `<nl>` sublayers of half-space material with a total
  depth of 1.5 shear wavelengths, `1.5 Vs_hs/f` (the depth changes with frequency), and puts viscous
  dashpots `ρVs` and `ρVp` at their base ("variable-depth method + viscous boundary").

The sublayer thickness law is `EDUOPT,HSLAW,UNIFORM` (default, lead decision D-W1-01, verified by
VP-02b), `GEOMETRIC` or `LINEAR`.

### 7.4 Wave field and control point (SITE options)

```
SITE,<opmode>,<mode1>,<fstep>,<nl>,<hs>,<mode2>,<wopt>,<freq1>,<freq2>,<cl>,<cm>,<delt>,<nft>,<freq>
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.005,8192,1
WAVE,<type>,<opt>,<ratio1>,<ratio2>,<angle>
* vertically incident SV wave
WAVE,2,1,1,1,0
```

| Argument | Meaning |
|---|---|
| `<mode1>`, `<mode2>` | run Mode 1 (layer eigen-solutions, FILE2) and Mode 2 (free-field motion, FILE1) |
| `<nl>`, `<hs>` | generated half-space sublayers; L number of the half-space |
| `<wopt>` | 0 in-plane waves (P, SV, Rayleigh), 1 anti-plane waves (SH, Love) |
| `<freq1>`, `<freq2>` | frequency numbers at which the WAVE ratios are given (linear in between) |
| `<cl>`, `<cm>` | control point: top of TOPL layer `<cl>`, direction 0 x', 1 y', 2 z' |
| `<delt>`, `<nft>`, `<fstep>`, `<freq>` | time step, NFFT, frequency step (0 = 1/(Δt·NFFT)), frequency set (section 9) |

**Wave types** (`WAVE <type>`): 1 Rayleigh, 2 SV, 3 P, 4 SH, 5 Love. Standard practice is vertically
propagating SV for the X input, SH for Y and P for Z (incidence angle 0). Inclined body waves
(`<angle>` > 0) and surface waves are available too. With several WAVE records the free field is
their mixture with the given ratios (they must sum to 1).

The free field is **normalised so that the control point moves with unit amplitude** in the control
direction: FILE1 and FILE8 are transfer functions per unit control motion. The control motion is
the **within** (in-column) motion at the top of layer `<cl>`. If your design motion is a rock
*outcrop* motion, run SOIL first (it deconvolves the outcrop motion through the column) and use the
SOIL surface motion as control motion at layer 1, as examples 4 and 6 do.

### 7.5 Strain-compatible properties: SOIL and SITEX

SOIL performs the SHAKE equivalent-linear analysis of the free-field column:

```
* SHAKE91 modulus and damping curves (labels Clay, Sand, Rock)
INP,../sassi/data/dynp_library.pre
* SOIL sublayer -> L properties + curve label (last SPRO = half-space)
SPRO,<layer>,<L>,<label>
* values, gravity, header lines, outcrop, save FILE88, 8 iterations, ratio 0.65
SOIL,4000,9.81,1,1,1,8,0.65,1,0
* horizontal input, factor 1, control at sublayer 23
SOILX,0,1,0,23
* save the surface acceleration ACC001.TH
SACC,1,2,0
RUNSOIL
* SITE uses the FILE88 properties ("Non-Linear Soil")
SITEX,1
```

The effective strain ratio is typically 0.5-0.7 (SHAKE91 suggests (M-1)/10 for magnitude M). Eight
iterations are usually enough; the listing prints the change of G and damping per iteration. For
an embedded structure, copy the strain-compatible properties into the L layers used by the
excavated soil as well (example 6 shows how): the excavated soil must have the free-field properties.
For strong shaking SOIL can also integrate the column in the time domain with a hyperbolic soil
model (SOIL-NON, `NLSOIL,1`, section 12.4); its FILE88 then holds the equivalent properties of the
nonlinear run.

---

## 8. Structure, excavation and interaction nodes

### 8.1 Structure or excavated soil: ETYPE

SOLID and PLANE elements are either structure or excavated soil:

| `ETYPE` | Meaning | Properties |
|---|---|---|
| 0 (default) | implicit: excavated if all nodes are at or below the ground elevation and the centroid below it, structure otherwise | as resolved |
| 1 | structure (including near-field soil and backfill) | M table (`MSET` = material number) |
| 2 | excavated soil (SOLID/PLANE) or buried shell (SHELL/TSHELL) | L table (`MSET` = soil layer number) |

`ETYPEGEN,<type>` sets ETYPE for all elements at once (0 resolves the implicit rule into explicit
1/2). HOUSE lists, for every excavated element, the layer it uses and warns (EDU-08) when that is not
the TOPL layer at the element's depth.

### 8.2 Interaction nodes

Interaction nodes are where the free-field impedance acts and the free-field motion is applied:

```
* set = 1 sets, 0 resets; code 0 = interaction node
INT,<n1>,<n2>,<inc>,<set>,[<code>]
* nodes 1..81 are interaction nodes
INT,1,81,1,1
INTLIST                               * list them
```

Rules (manual §4.1.2 item 13; checked by CHECK, HOUSE and the model-check commands):

* interaction nodes belong to the **excavated soil**, not to the structure; nodes shared by the
  structure and the outer surface of the excavation are interaction nodes; no other structural node
  may be one;
* no **interior** excavation node may be connected to a structural element (`EXCSTRCHK`);
* every interaction node lies **on a soil-layer interface at or below the ground surface** (EDU-01);
* interaction nodes of embedded models are numbered **bottom-up** in ascending order (EDU-21);
* an interaction node must not be fixed (`FIXEDINT`, Error 124).

For a **surface foundation** the interaction nodes are the foundation nodes at grade; there is no
excavated soil and FV, FI and FFV are the same.

### 8.3 Choosing the interaction-node set: FV, FI and FFV

The flexible-volume method is exact when *every* node of the excavated soil is an interaction node.
Fewer interaction nodes are cheaper, because the impedance matrix is dense
(`(3 N_int)² × 16` bytes):

| Method | Interaction nodes | `INTGEN` type | Remarks |
|---|---|---|---|
| **FV** (flexible volume, "direct") | all excavated-soil nodes | 1 | reference method |
| **FFV** (fast flexible volume) | FI-EVBN + internal horizontal levels every `skip` levels | 5 | close to FV at a fraction of the cost |
| **FI-EVBN** (modified subtraction, MSM) | lateral, bottom **and top** faces of the excavation | 2 | much better than FI-FSIN |
| **FI-FSIN** (subtraction, SM) | lateral and bottom faces (the soil-foundation interface) | 3 | spurious resonances at and above the excavated-volume frequencies |
| surface | foundation nodes at grade | 4 | surface foundations |

`INTGEN,<type>,[skip]` generates the set on an excavation mesh (types are additive; 0 clears). HOUSE
`<imp>` records the method (0 FV, 1 FFV, 2 FI); the computation is the same, only the set differs.
ASCE 4-16 and SRP 3.7.2 require that a non-FV method be **validated against FV** (CHECK warning
EDU-12): compare the transfer functions at the common structure/excavation nodes. Example 2 shows the
spurious peak of FI-FSIN at 6 Hz (roof amplitude 1.34 against 0.85 with FV, +58 %) and its removal by
FI-EVBN. Example 8 does the same comparison on a full shear-wall building with a two-level basement:
FI-FSIN resonates at 16 Hz (2 to 3 times the FV transfer functions), FI-EVBN stays within 2.4 % of FV.

### 8.4 The central-zone radius

POINT spreads each unit point load over a small axisymmetric **central zone** of radius R0, because a
true point load has an infinite displacement under it. R0 calibrates the diagonal terms of the
flexibility matrix to the mesh:

| Interaction-node mesh | Rule |
|---|---|
| 3D square mesh of size h | R0 = 0.90 h |
| 3D triangular (circular foundation) mesh of size h | R0 = 0.85 h |
| 2D mesh of size h | R0 = h |

```
POINT,<opmode>,<layer>,<rad>
* surface foundation (no embedded layers), R0 = 0.9 x 2.5 m
POINT,0,0,2.25
* embedment over 5 TOPL layers: point loads at interfaces 1..6
POINT,0,5,2.25
```

For non-uniform meshes `RADIUS,<scale>,<file>` lists the equivalent radius of every excavated
element and their average; run sensitivity studies with the minimum, average and maximum. `<layer>`
must cover the deepest interaction node.

### 8.5 Mesh-size rules for the excavation

* vertical element size ≤ Vs/(5 f_cut), as for the soil layers;
* horizontal size: equal to the vertical size unless a sensitivity study justifies more (1.2-1.5,
  sometimes 2 times the vertical size is often acceptable for vertically propagating waves);
* the excavated element sizes are set by the interaction-node spacing; keep basement nodes that
  belong to internal structures separate from the excavation nodes.

### 8.6 Generation tools

| Command | Use |
|---|---|
| `EXCAV,<model>,[delta]` | generate the excavation volume from the lowest basement grid |
| `SOILMESH,...` | near-field soil mesh around a basement |
| `MERGE`, `MERGESOIL` | join a structure model and an excavation model (coincident nodes merged, or stiff springs) |
| `WELD`, `NCOM`, `GCOM`, `RMVUNUSED` | merge coincident nodes, compress numbers, remove unused nodes |
| `ROTATE`, `TRANSLATE` | move the whole model |
| `ETYPEGEN`, `INTGEN` | ETYPE and interaction nodes for the whole model |

Verification problems VP-50 and VP-51 check these tools on small models.

---

## 9. Frequencies and the FFT

### 9.1 Time step, NFFT and the frequency step

The analysis frequencies are integer **frequency numbers** `n` times the frequency step

```
Δf = 1 / (Δt · NFFT)          (or SITE <fstep> when it is > 0)
f  = n · Δf
```

with Δt the time step of the control motion and NFFT the number of FFT points (a power of 2, at most
32,768 for the input modules and 65,536 for MOTION/RELDISP/STRESS). The Fourier period `NFFT·Δt`
must be longer than the record plus a **quiet zone** of zeros in which the response decays; otherwise
the free vibration wraps around to the start of the record (EDU-19 warns). Example: Δt = 0.005 s,
NFFT = 8192 gives Δf = 0.0244 Hz and a 41 s period for a 20 s record.

### 9.2 Choosing the SSI frequencies

```
* append frequency numbers to set <set>; FREQ,<set>,0 deletes it
FREQ,<set>,<n1>,...,<n10>
LFREQ                          * list the sets with their values in Hz
```

* 40-80 frequencies for stick models, 100-300 for large FE models, at least 200 from the start for
  incoherent analyses (manual §4.1.2 item 10); at most 500 per ANALYS run.
* Put more frequencies where the response changes fast: around the SSI frequencies of the
  structure (often 0.5-0.7 times the fixed-base frequencies) and around the soil column resonances.
* Start low (about 0.1-0.5 Hz): below the first SSI frequency a seismic transfer function is
  interpolated linearly to its rigid-body value (1 in the input direction) at f = 0.
* The last frequency is the **cut-off**: transfer functions are zero above it, so the cut-off acts
  as a low-pass filter on every result. Typical values are 30-40 Hz for soil sites and 60-70 Hz for
  rock sites; the layers and the excavation mesh must pass it.

### 9.3 Interpolation between SSI frequencies

MOTION and STRESS interpolate the complex transfer functions from the SSI frequencies to every
Fourier frequency. Options 0-5 fit the transfer function of a **two-degree-of-freedom system with
hysteretic damping** through windows of five consecutive SSI frequencies (Tajirian's method); they
reproduce the computed values exactly and resolve resonance peaks between computed points:

| Option | Windows | Combination |
|---|---|---|
| 0 | SASSI2000: all overlapping windows | weighted average |
| 1 | original SASSI (1982): non-overlapping windows | single window (default) |
| 2 | all overlapping windows | average |
| 3 | the three closest overlapping windows | average |
| 4 / 5 | non-overlapping windows shifted by one / two points | single window |
| 6 | complex not-a-knot cubic splines of Re H and Im H | needs denser frequencies |

`MOTION <smo>` (smoothing parameter S) damps interpolated peaks that leave the band of the two
neighbouring computed values; keep S = 0 for coherent analyses (EDU-14 warns), because genuine
resonance peaks are clipped too. `MOTION <pzadj>` = 1 reduces the phase toward zero ("phase
adjustment", an upper-bound approach of incoherent analyses).

**Always compare** the computed `.TFU` and the interpolated `.TFI` transfer functions (GUI: Plot >
Transfer Functions). An interpolated peak in the middle of an interval, much higher than the
computed neighbours, means that a frequency is missing.

### 9.4 Adding frequencies

```
* frequencies where TFI peaks exceed TFU -> variable
CRITFREQ,<tol>,<minfilter>,<TF>,<Var>
```

Run ANALYS for the added frequencies only (they must exist in FILE1/FILE9 and FILE3, otherwise run
SITE and POINT for them too), then merge:

```
* the first run
FMOVE,FILE8,FILE81
...                           * ANALYS on the new frequency set writes FILE8
FMOVE,FILE8,FILE82
RUNCOMBIN                     * FILE81 + FILE82 -> FILE8 (duplicate frequencies are an error)
```

`REMOVEFREQ,<infile>,<outfile>,<n1>,...` removes spurious frequencies from a FILE8.

---

## 10. Running ANALYS

```
ANALYS,<opmode>,<type>,<mode>,<save>,<prnt>,<fopt>,<ang>,<xc>,<yc>,<zc>,<impe>,[simul]
* seismic initiation run, control point at the origin
ANALYS,0,0,0,0,1,0,0,0,0,0,0
```

| Argument | Values |
|---|---|
| `<type>` | 0 seismic (FILE1), 1 foundation vibration (FILE9) |
| `<mode>` | 0 initiation; restarts: 1 New Structure, 2 New Seismic Environment, 3 New Dynamic Loading |
| `<save>` | 1 saves the restart files COOXqqq (impedance) and COOTKqqq (factorised system) |
| `<prnt>` | 1 amplitudes only in the listing, 0 real and imaginary parts |
| `<fopt>` | 1 solves every frequency of FILE1/FILE9 instead of the frequency set |
| `<ang>` | angle from the SITE x' axis to the global x axis (degrees, counter-clockwise about z) |
| `<xc>`, `<yc>`, `<zc>` | control point: phase origin of inclined/surface waves and reference of the global impedance |
| `<impe>` | global (unconstrained) foundation impedance: 0 none, 1 diagonal, 2 full 6×6 (3D models without SYMM planes) |
| `<simul>` | simultaneous cases: X, Y and Z input, load cases, or incoherent stochastic samples (section 10.1) |

`ANALYSX,<ffm>,<delrst>` adds two dialog fields: how incoherency factors are applied (0 free-field
load, the default, or 1 free-field motion; section 12.1) and the deletion of the restart files after a
successful restart run.

At every SSI frequency ANALYS computes the flexibility `Fff` of the free field at the interaction
DOFs from FILE3, inverts it to the impedance `Xff`, assembles the flexible-volume equation with the
HOUSE matrices, builds the load (free field of FILE1 or forces of FILE9) and solves it by a Schur
complement on the interaction DOFs (sparse LU of the rest, dense LU of the condensed system). FILE8
then holds the transfer function of every degree of freedom. The listing prints a per-frequency
summary and a **low-frequency check**: at the first frequency every node should move with the
control motion (|ATF| close to 1 in the input direction, G-19). A large deviation points to a
modelling error (missing interaction nodes, a structure not connected to the soil, an
ill-conditioned model).

### 10.1 Simultaneous cases

* **X, Y and Z input in one run** (`<simul>` = 1): run SITE three times (SV with `<cm>` = 0, SH with
  `<cm>` = 1, P with `<cm>` = 2, angle 0) and copy each FILE1 to FILE1X, FILE1Y, FILE1Z with `FCOPY`.
  ANALYS factorises the system once per frequency and writes FILE8X, FILE8Y and FILE8Z. Run MOTION
  once per direction with `EDUOPT,TFFILE,FILE8X` (Y, Z). See example 5.
* **Load cases** (vibration, `<simul>` = Nl ≥ 2): run FORCE once per load case and copy FILE9 to
  FILE9001 ... FILE9Nl; ANALYS writes FILE8001 ... FILE8Nl. See example 3.
* **Incoherent stochastic samples** (`<simul>` = Ns ≤ 50, with HOUSE samples FILE77001 ... FILE77Ns):
  each sample is applied to FILE1X, FILE1Y and FILE1Z; ANALYS writes `FILE8{3(s-1)+d}`, i.e. FILE8001,
  FILE8002, FILE8003 for sample 1 in X, Y, Z, then FILE8004 ... (section 12.1).

### 10.2 Restarts

| Change | Restart | Re-run | Reused |
|---|---|---|---|
| only the control-motion history | "New Time History" | MOTION / STRESS / RELDISP | FILE8 |
| structure or near-field soil; same interaction nodes and soil | `<mode>` = 1 New Structure | HOUSE → ANALYS | impedance (COOX) |
| wave field or control point | `<mode>` = 2 New Seismic Environment | SITE → ANALYS | impedance and factorised system (COOX, COOTK) |
| external load pattern | `<mode>` = 3 New Dynamic Loading | FORCE → ANALYS | impedance and factorised system |

The initiation run must save the restart files (`<save>` = 1). A restart checks that the
interaction nodes, the soil layering and (for Mode 3) the structure are unchanged (hashes in FILE90)
and refuses mismatching files. A restart gives the same FILE8 as a fresh run (VP-22). New Structure
restarts are the basis of the rigid-basement study recommended by the manual (multiply the basement
modulus by 10⁴ and 10⁵), of the near-field soil iterations (section 12.4) and of the Option NON
iterations (section 12.5).

### 10.3 Foundation impedance

With `<impe>` = 1 or 2 ANALYS condenses the impedance on the interaction nodes to a 6×6 matrix about
the control point, `KG = Tᵀ Xff T` (T = rigid-body transformation), and writes

| File | Content |
|---|---|
| `FOUNSTIF` | Re KG (stiffness) |
| `FOUNDASH` | Im KG / ω (dashpot coefficients) |
| `FOUNDAMP` | Im KG / (2 \|Re KG\|) (equivalent damping ratio) |
| `FOUNIMPD` | \|KG\| (modulus) |

This is the *unconstrained* impedance of the soil at the interaction nodes; it equals the impedance
of a rigid foundation only for a surface foundation whose nodes all move as a rigid body. Example 3
compares it with the impedance obtained from a rigid mat loaded by unit forces and moments.

### 10.4 Memory and run time

The dense impedance matrix needs `(3 N_int)² × 16` bytes: 576 MB for 2,000 interaction nodes, 14.4
GB for 10,000. `INTCOUNT` prints the estimate; ANALYS checks it against the available memory before
it starts (EDU-20). Run time grows with the cube of the number of interaction nodes, so FFV or
FI-EVBN instead of FV, a coarser excavation mesh (after a sensitivity study) or symmetry planes
(`SYMM`: a half model keeps a little more than half of the interaction nodes, a quarter model a little
more than a quarter; section 12.3) are the main levers. For large models split the frequency set into several ANALYS runs and merge
the FILE8 with COMBIN; `AFWRBAT` prepares such a split run in separate folders that can run in
parallel or on other computers (section 12.9).

---

## 11. Post-processing

### 11.1 MOTION: transfer functions, accelerations and ISRS

```
MOTION,<opmode>,<out>,<step>,<dur>,<res>,<freq1>,<freq2>,<fstep>,<mult>,<max>,<rec1>,<rec2>,
       <fopt>,<bl>,<smo>,<cplx>,<cnvrt>,<pzadj>,<interp>
MOTION,0,0,0,20,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1
* ISRS damping ratios (MOTION uses up to 5)
DAMP,0.02,0.05
* control motion in g (first line dt, then values)
THFILE,../data/rg160h_030g.acc
* X: TF, history, ISRS, maxima of nodes 41, 82..85
NOUT,1,1,1,0,0,1,1,41,82-85
* Z: TF and maxima of nodes 37 and 45
NOUT,3,1,0,0,0,0,1,37,45
```

| Argument | Meaning |
|---|---|
| `<dur>` | duration of interest; the `.ACC` histories cover 1.2 × `<dur>` (0 = the whole Fourier period) |
| `<freq1>`, `<freq2>`, `<fstep>` | ISRS frequency range and number of log-spaced points (≥ 301 recommended, EDU-15) |
| `<mult>`, `<max>` | scale the control motion by a factor, or to a peak value (exactly one non-zero) |
| `<rec1>`, `<rec2>`, `<fopt>` | records used; file format (0: dt then one value per line, 1: time-value pairs) |
| `<bl>` | 1 = Hudson-Housner baseline correction of the output histories |
| `<smo>`, `<pzadj>`, `<interp>` | smoothing, phase adjustment, interpolation option (section 9.3) |
| `<cplx>` | 1 = save complex `.TFI` (needed by RELDISP) |

`NOUT,<dir>,<six flags>,<nodes>` requests output for direction 1 X, 2 Y, 3 Z, 4 XX, 5 YY, 6 ZZ;
the flags are (TF print, save history, plot history, plot RS, save RS, print maximum). Output files:

| File | Content |
|---|---|
| `00041TR_X.TFU` | computed TF at the SSI frequencies (frequency, amplitude, phase) |
| `00041TR_X.TFI` | interpolated TF on the Fourier grid up to the cut-off |
| `00041TR_X.ACC` | acceleration history in g (vibration analysis: displacement, velocity or acceleration per `MOTIONX`) |
| `00041TR_X01.RS` | response spectrum (absolute acceleration, g) for damping number 01 |
| `<model>_MOTION.out` | maxima (ZPA), input summary |

The seismic response is `a(t) = IFFT[H(f) A(f)]` with A the FFT of the control acceleration; for a
vibration analysis the reference load history (THFILE) replaces the control motion.

### 11.2 RELDISP: relative displacements

Relative displacements are computed analytically from the complex `.TFI` files, which avoids the
drift of double integration:

```
* save the complex .TFD; no all-node output
RELD,1,0,0
* reference node and direction (empty: the free field)
RELFILE,00041TR_X.TFI
* node 85, X component
RDND,85,1,0,0,0,0,0
```

`d(t) = IFFT[(H_node - H_ref) U_g]` with `U_g = -g A/ω²` the control displacement spectrum. Results:
`00085TR_X.THD` (history) and `.TFD` (transfer function, length per g). One DOF per run (the DOF of
the reference file).

### 11.3 STRESS: element forces and stresses

```
STRESS,<opmode>,<iter>,<save>,<itran>,<interopt>
STRESS,0,0,1,1,1
* per component: 0 none, 1 maximum, 2 maximum + history
EOUT,<12 codes>,<group>,<elements>
* BEAMS group 2: histories of FYI and MZI
EOUT,1,2,1,1,1,2,1,1,1,1,1,1,2,1-4
```

| Element | Components (EOUT order) |
|---|---|
| SOLID | SXX SYY SZZ SXY SXZ SYZ SOCT (octahedral shear stress) |
| BEAMS | FXI FYI FZI MXI MYI MZI FXJ FYJ FZJ MXJ MYJ MZJ (local axes, forces on the element) |
| SHELL | FXX FYY FXY (membrane stresses) MXX MYY MXY (moments per unit length), local axes |
| TSHELL | NXX NYY NXY QXZ QYZ MXX MYY MXY (forces and moments per unit length; `THSHLSTR,1` adds the face stresses N/t ± 6M/t² and strains) |
| PLANE | SXX SZZ TXZ |
| SPRING | FX FY FZ MXX MYY MZZ |

STRESS interpolates the *stress* transfer functions (not the nodal ones) and convolves them with
the control displacement spectrum. Files: `BEAMS_002_00001_MZI.THS` (history), `.TFU`/`.TFI` with
`<itran>` = 1, the maxima in `<model>_STRESS.out`, `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`, and with
`SECDATAOPT,1` the frames `ESTRESS_n.ess` used by the section cuts.

### 11.4 Plots

In the GUI, the Plot menu draws the model (3D), transfer functions (TFU against TFI), histories,
spectra, soil layers and soil-property curves; plots can be exported as PNG. The same plots are
available as commands: `MODELPLOT`, `NODEPLOT`, `LAYERPLOT`, `SOILPROPPLOT,<label>`, `SPECPLOT,<lines>`,
`THPLOT,<lines>`, with `PLOTTITLE`, `XTITLE`, `AXES`, `PLOTRANGE`, `LINENAME`, `MARKERS` and
`CAPTUREPLOT,<file.png>`. Line objects are created by `READSPEC` (spectrum files) and `READTH`
(histories).

Animations play frame files on the model: the restart frames of MOTION, RELDISP and STRESS, or the
steady-state motion at one SSI frequency that `HARMFRAME` (SASSI-EDU extension) writes from FILE8,
u(t) = Re(H e^{iωt}) per unit control motion over one period:

```
* the SSI mode of example 1: every node at 3.49 Hz, 24 frames, drawn 0.2 m per unit of control motion
HARMFRAME,FILE8,3.49,HARM_SSI
PROCFRAME,HARM_SSI,HARM_SSI_ani,SSI system at 3.49 Hz,3
DEFORMPLOT,HARM_SSI_ani,1,24,1,0.2
WINDOWSETTINGS,UNDEFORMED,1
```

`HARMFRAME,<Src>,<Freq>,<OutDir>,[NFrames],[Ref]` uses the computed frequency closest to `<Freq>` and
prints it; `<Ref>` 0 gives the motion relative to the free field, a node number the motion relative to
that node ([GUI §6.3](GUI.md#63-animations-requirements-58)).

### 11.5 Line mathematics and spectrum broadening

Line operations work on the union of the abscissas of their sources:

```
* line 1 = roof ISRS from the X run
READSPEC,00085TR_X01.RS,1,1
* line 2 = from the Y run
READSPEC,00085TR_Y01.RS,1,2
* line 3 = SRSS of lines 1 and 2
SRSS,3,1,2
* line 4 = envelope, +-15 % peak broadening, peak bridging
BROADEN,4,15,15,3
WRITESPEC,roof_design.rs,4
```

`ADDITION`, `SUBTRACTION`, `LINECOMBIN` (coefficients), `AVERAGE` (e.g. the mean of several
simulations) and `SRSS` combine lines; `BROADEN,<dest>,<Smooth1>,<Smooth2>,<sources>` forms the
envelope, broadens the peaks by ±`Smooth2` % (Regulatory Guide 1.122 practice) and bridges valleys
between peaks within `Smooth1` %.

### 11.6 Section cuts

To obtain the resultant forces on a cut through a solid or shell structure (for example the base
shear and overturning moment of a shear wall), let STRESS write the element-centre stresses of
selected time steps (`SECDATAOPT,1` and a `Frames.txt` list) and then, as in manual §5.8.2:

```
* element stresses of one time step onto the model
READSTR,NSTRESS/ESTRESS_00101.ess
* cut 3 = the elements of one wall (box selection)
CUTVOL,3,52.5,52.8,-320,320,2.53,45.22
* cross-section model 1 through the plane x = 52.65
CSECT,1,3,52.65,0,4,1,0,0
ACTM,1
* area, centroid, inertias and the six resultants
CALCPAR,1,0,0,0,1,0,10
ACTM,0
* resultant histories
CALCSECTHIST,NSTRESS/ESTRESS.lst,3,52.65,0,4,1,0,0,0,1,0,10,0.01,wall.csv
```

`CUTADD`/`CUTRMV` (element lists) and `SLICE` (elements crossing a plane) build cuts too;
`CALCMOI`, `CALCC` and `CALCM` give section inertias, centroids and masses. `CALCSECTHIST` writes a
7-column CSV (time, three forces, three moments) with a final row of signed absolute maxima. VP-49
checks the algebra on two small models.

---

## 12. Advanced options

Each subsection says what the option does, gives a short command sequence, points to the theory and
names the verification problems. Section 2.4 lists the module sequences of these analyses.

### 12.1 Incoherency, wave passage and multiple excitation

Two points of a foundation a distance D apart do not receive exactly the same motion: the waves arrive
at different times (*wave passage*) and are partly uncorrelated (*incoherence*, measured by the
coherency γ(f, D), which decays with frequency and distance). For large stiff foundations, typically
nuclear islands on rock, incoherency reduces the high-frequency translational response and adds
rotations. HOUSE computes complex **incoherency factors** for every SSI frequency, direction and
interaction node (FILE77): it builds the coherency matrix of the interaction nodes, decomposes it into
*incoherent spatial modes* and combines them; ANALYS multiplies the free-field load by the factors.
The equations are in the [Theory Manual §14](../theory/THEORY_MANUAL.md#14-incoherency-wave-passage-and-multiple-excitation).

| Command | Meaning |
|---|---|
| `HOUSE,<gravity>,<gelev>,<opmode>,<dim>,<imp>,<coh>,<wpass>,<me>,<cmplxspec>` | arguments 6-9: 1 = incoherent motion, wave passage, multiple excitation, complex amplification ratios |
| `WPASS,<appv>,<ang>,<cohf>` | apparent velocity and angle of Line D (the direction of wave passage); coherency model 1-7 |
| `INCOH,<gammax>,<gammay>,<gammaz>,<alpha>,<ngp>,<ipr>,<nmodes>,<met>,<HSeed>,<VSeed>,<RandPhz>` | Luco-Wong parameters per component (model 1, ≥ 0.1); directionality factor α (models 2-7, 0.5 = isotropic) or mean Vs (model 1); embedded levels; I N C O listing; modes (0 all, k the first k, -k mode k only); units (0 British, 1 SI); seeds and random phase angle (degrees) |
| `HOUSEX,<optimize>,<supmode>,<nsim>` | superposition 0 Linear / 1 Quadratic (SRSS TF); number of stochastic samples (≤ 50) |
| `ME,<no>,<nfirst>,<nlast>,...`, `AMP,<no>,<a1>,...` | multiple-excitation zone (interaction nodes nfirst ... nlast) and its spectral amplification ratios, one per SSI frequency |
| `ANALYSX,<ffm>` | 0 free-field load (FFL, default), 1 free-field motion (FFM, surface foundations) |
| `BUILDFILE77,<out>,<in1>,...` | combine the FILE77 of per-level HOUSE runs (deeply embedded foundations) |
| `EDUOPT,INCOHSIGN,ADJUST\|RAW`, `EDUOPT,INCOHMERGE,1` | sign convention of the modes (ADJUST: f → 0 gives coherent motion); merge coincident plan positions |

**Coherency models** (`WPASS <cohf>`). The coefficients of the Abrahamson models were transcribed from
the EPRI reports into `sassi/data/coherency/*.json`, with the report, equation, table and page of every
number; a model whose source could not be retrieved refuses to run with *coefficients not available*.

| `<cohf>` | Model | In SASSI-EDU |
|---|---|---|
| 1 | Luco and Wong (1986), `exp[-(γ ω D/Vs)²]` | available |
| 2 | Abrahamson (1993), all soil types | refused (coefficients not available; use 3 or 7) |
| 3 | Abrahamson (2005), all sites, surface foundations | available (EPRI 1012968) |
| 4 | Abrahamson (2006), all sites, embedded foundations | refused (coefficients not available; EPRI recommends 5 for embedded foundations) |
| 5 | Abrahamson (2007), hard rock | available (EPRI 1015110) |
| 6 | Abrahamson (2007), soil sites, surface foundations | available (EPRI 1015110) |
| 7 | user tables `COHXUSER`, `COHYUSER`, `COHZUSER` on `FREQCOH` × `DISTCOH` | available |

Models 2-7 are applied with the wave-passage option on (`<wpass>` = 1); give `WPASS,1e9,...` when you
want no delay.

**Approaches.** The manual's reference is the stochastic simulation; the deterministic approaches are
approximate and validated for stick models on rigid mats only (EDU-13).

| Approach | Input | Runs and results |
|---|---|---|
| stochastic simulation ("Simulation Mean") | non-zero seeds and RandPhz (180 = random phases on the full circle), HOUSEX `<nsim>` = Ns, Linear | one HOUSE run writes FILE77001 ... FILE77Ns; ANALYS `<simul>` = Ns writes FILE8001 ... FILE8(3Ns); MOTION per sample; the mean of the spectra (`AVERAGE`) |
| algebraic sum (AS) | zero seeds or RandPhz = 0, Linear | one FILE77, one FILE8 per direction |
| SRSS TF (Quadratic) | `<nmodes>` = -k, `<supmode>` = 1 | one HOUSE + ANALYS run per mode k; MOTIONX `<srss>` with `SRSSTF.txt` combines the modal transfer functions |
| SRSS FRS | `<nmodes>` = -k, `<supmode>` = 0 | one run per mode; SRSS of the modal spectra (`SRSS` line operation) |

A stochastic simulation of a surface foundation with the 2007 hard-rock model (X, Y and Z input as in
example 5, SI units):

```
* 3D, incoherent motion with wave passage
HOUSE,9.81,0,0,2,0,1,1,0,0
* no delay (V_app = 1e9), Line D along X, model 5
WPASS,1e9,0,5
* alpha 0.5, I N C O listing, all modes, SI, seeds, +-180 deg
INCOH,0.1,0.1,0.2,0.5,1,1,0,1,1975,2026,180
* Linear superposition, 20 samples
HOUSEX,0,0,20
* 20 samples x FILE1X/Y/Z -> FILE8001 ... FILE8060
ANALYS,0,0,0,0,1,0,0,0,0,0,0,20
* MOTION of sample 1, X input (then FILE8004, FILE8007 ...)
EDUOPT,TFFILE,FILE8001
```

For each sample copy the spectra under their own names (`FCOPY`), read them with `READSPEC` and form the
mean with `AVERAGE`. Rules of the manual: 3D models only and no SYMM planes (CHECK: EDU-26),
interaction nodes numbered bottom-up (EDU-21 is an error for incoherent analyses), at least 200 SSI
frequencies from the start, interpolation option 6, and stochastic simulation without phase adjustment
(`MOTION <pzadj>` = 0) as the reference approach.

**Wave passage** alone (`<coh>` = 0, `<wpass>` = 1) delays the coherent motion along Line D by
`τ = d/V_app`, measured from the ANALYS control point. **Multiple excitation** (`<me>` = 1, which needs
`<wpass>` = 1) scales the motion of each foundation zone by its spectral amplification ratios.

*Verification:* VP-27 (decomposition, f → 0 limit, distances, wave-passage phase), VP-I1 (rigid mat
under Luco-Wong coherency against an independent rigid-foundation average, FFL and FFM, 1e-8), VP-I2
(stochastic samples: the mean of s sᴴ converges to the coherency matrix as 1/√Ns), VP-26 (FFL versus
FFM: identical where the theory makes them identical; lead decision D-W3-05 explains why the manual's
"identical results" holds only for the generalised force in the input direction and for decoupled
foundation DOFs).

### 12.2 Two-dimensional (plane-strain) SSI

`HOUSE <dim>` = 1 makes the model a plane-strain slice of unit thickness in the X-Z plane: PLANE
elements for the soil and plane structures, POINT2 line-load solutions with transmitting boundaries,
interaction DOFs UX and UZ, and ANALYS with the in-plane impedance. Masses, stiffnesses and loads are
per unit length.

```
* <dim> 1: 2D plane strain (POINT2)
HOUSE,9.81,0,0,1,0,0,0,0,0
* excavated soil (ETYPE 2) and plane structures, nodes in the X-Z plane
GROUP,1,PLANE
* 2D rule: central-zone radius R0 = h, the interaction-node spacing
POINT,0,2,1.0
* <simul> 1: the in-plane cases X and Z (FILE1X, FILE1Z -> FILE8X, FILE8Z)
ANALYS,0,0,0,0,1,0,0,0,0,0,0,1
```

Restrictions: in-plane input only (SV, P or Rayleigh waves; the anti-plane SH/Love case is not
implemented), the coordinate transformation angle 0 or 180°, no global impedance and no incoherency.
BEAMS, SPRING and GENERAL nodes keep six DOFs: fix the out-of-plane ones (UY, ROTX, ROTZ). The manual
does not recommend 2D SSI for design (it overestimates radiation damping, which is unconservative); use
it for sensitivity studies, for example of non-horizontal layering. Theory:
[Theory Manual §15.1](../theory/THEORY_MANUAL.md#151-plane-strain-ssi-house-dim--1). *Verification:*
VP-T1 (POINT2 far field against Kausel's line-load solution, 0.33 %), VP-T3 (2D zero-SSI identity), VP-42
(rigid strip on a layer against Jakub and Roesset, 5 %).

### 12.3 Symmetry planes: half and quarter models

A structure symmetric about one or two vertical planes, under a loading that is symmetric or
antisymmetric about them, can be analysed as a half or a quarter model. `SYMM,<no>,<type>,<n1>,<n2>,<n3>`
defines plane `<no>` by three nodes of the model that lie in it (two in 2D) with `<type>` 0 symmetric or
1 antisymmetric loading; the planes must be parallel to the XZ or YZ planes. HOUSE fixes the DOFs the
symmetry makes zero on the plane nodes, and ANALYS forms the soil impedance of the reduced set of
interaction nodes from their mirror images
([Theory Manual §15.2](../theory/THEORY_MANUAL.md#152-symmetry-and-antisymmetry-planes-symm)).

```
* quarter model x >= 0, y >= 0 under the X input (vertical SV, x' = X)
* plane x = 0: antisymmetric (the X input reverses in the mirror image)
SYMM,1,1,101,105,141
* plane y = 0: symmetric
SYMM,2,0,101,121,141
```

(nodes 101 and 141 lie on the vertical line x = y = 0, node 105 in the plane x = 0 and node 121 in the
plane y = 0). For the Y input swap the types; for the Z (P-wave) input both planes are symmetric. The
plane types belong to the HOUSE model, so each input direction is a separate run. Elements, masses and
loads lying *in* a plane get 1/2 of their full-model values (1/4 on the intersection of two planes).
ANALYS checks that the free field has the declared symmetry and stops otherwise. Not allowed with
SYMM: incoherency, wave passage, multiple excitation and the global impedance. *Verification:* VP-T2
(half and quarter models of a surface mat and of an embedded basement reproduce the full model to 4e-11).

### 12.4 Nonlinear soil: SOIL, near-field iterations and SOIL-NON

| Kind | In SASSI-EDU |
|---|---|
| Primary (free-field) soil nonlinearity, equivalent linear | **SOIL** (SHAKE, section 7.5) |
| Primary soil nonlinearity, nonlinear in the time domain | **SOIL-NON**: `NLSOIL`, `NLSLAYER`, `DELNLS` (below) |
| Secondary (near-field) soil nonlinearity: equivalent-linear SSI iterations on soil elements of the structure model | `PIN`, `PINGRP`, `PINMAT` (the `.pin` input), `HOUSEX,0,0,1,1`, `NLSSIRESET`, `NLSSIITER`, `COMBXYZSTRAIN`; example 6 |

**Near-field soil iterations** (example 6). The soil next to the structure (a loose backfill, the soil
under a heavy neighbouring building) is modelled with SOLID (PLANE in 2D) elements of the *structure*
model (ETYPE 1) at the nodes of the excavated soil, whose shear modulus and damping are iterated element
by element:

1. Give the near-field elements a material with the **low-strain** soil properties (G_max): the curve
   softens it. `PIN,<esf>` sets the effective-strain factor; `PINGRP,<igrp>,<istr>,<gfac>,<dfac>,<curve>`
   declares a nonlinear group (strain measure 0 = largest shear-strain component, 1 = octahedral; GFAC
   and DFAC scale the free-field G and damping at the element for the first run, 1 = the free field; the
   soil curve, e.g. `Sand`); `HOUSEX,0,0,1,1` switches on non-linear SSI in HOUSE.
2. The initiation run (HOUSE, ANALYS with `<save>` = 1, three directions) is followed by STRESS with
   `<iter>` = 1 for X, Y and Z (effective strains → FILE74X/Y/Z) and `COMBXYZSTRAIN` (SRSS → FILE74).
3. `NLSSIITER,<variables>,<maxit>` repeats HOUSE (new G and β from the strains and the curves) → ANALYS
   New Structure restart → STRESS ×3 → COMBXYZSTRAIN until max |ΔG/G| < 2 % and max |Δβ| < 0.5 %.

```
* effective strain = 0.65 x maximum strain
PIN,0.65
* group 3: octahedral strain, start at 0.4 G_free-field, curve Sand
PINGRP,3,1,0.4,1.0,Sand
* HOUSE non-linear SSI on (AFWRITE writes <model>.pin)
HOUSEX,0,0,1,1
NLSSIRESET                               * a new analysis: the next HOUSE run takes the .pin properties
...                                      * initiation run, STRESS x 3, COMBXYZSTRAIN (example 6)
* <mode> 1 = New Structure restart for the iterations
ANALYS,0,0,1,0,1,0,0,0,0,0,0,1
* HOUSE, ANALYS, STRESS x 3, COMBXYZSTRAIN until converged
NLSSIITER,NLRUN+NLX+NLY+NLZ+NLCOMB,8
```

Each iteration is a New Structure restart, 2-4 times faster than an initiation run. The history is in
`NLSOIL_CONVERGENCE.TXT` and the HOUSE listing. Theory:
[Theory Manual §16](../theory/THEORY_MANUAL.md#16-nonlinear-soil-ssi-iterations). *Verification:* VP-N1
(a column of near-field soil identical to the free field reproduces SOIL within 1.3 %), VP-N2 (restart
consistency).

**SOIL-NON** replaces the equivalent-linear iterations of SOIL by a nonlinear time-domain analysis of the
column (DEEPSOIL-type): a modified hyperbolic backbone `τ = G0 γ / (1 + β(γ/γ_r)^s)` per sublayer with
Masing unloading and reloading, small-strain viscous damping, a rigid or an elastic (Joyner-Chen) base,
implicit Newmark integration with Newton iterations
([Theory Manual §18](../theory/THEORY_MANUAL.md#18-soil-non-nonlinear-time-domain-site-response)). The
input must be the motion at bedrock (control layer = the last SPRO sublayer), as an outcrop motion with
the elastic base or a within motion with the rigid base.

```
* SOIL-NON on; flexible sub-steps; default tolerances; elastic base; damping type 1
NLSOIL,1,0,0,0,0,1,1,0,0
* sublayer 1: fit the hyperbolic model to its DYNP G/Gmax curve
NLSLAYER,1,1
* sublayer 2: beta 1, s 0.9, reference strain 0.03 %, no viscosity
NLSLAYER,2,0,1.0,0.9,0.03,0
...                              * one set per soil sublayer (Error 125), none for the half-space
AFWRITE                          * writes <model>.soi and the side file <model>.nls
RUNSOIL
```

The output files are those of SOIL-EQL (`ACCxxx.TH`, `SNxxx.TH`, `SSxxx.TH`, spectra, FILE88 with the
equivalent properties of the nonlinear run). *Verification:* VP-SN1 (linear limit = SOIL within 0.4 %),
VP-SN2 (Masing loop damping), VP-SN3 (SOIL-NON against SOIL-EQL on the SHAKE91 sample, informative).

### 12.5 Nonlinear structures: Option NON

Option NON (module NONLINEAR) includes **cracking reinforced-concrete wall panels** and **nonlinear
springs** (isolators, pile-soil interfaces) by equivalent-linear SSI iterations: the deformation history
of each element from RELDISP (X + Y + Z, combined by COMB_XYZ_THD) is run through a hysteresis model
(Cheng-Mertz shear for panels, General Masing Rule for springs; Takeda experimental), the panel modulus or
spring constant is set to the secant value at 0.8 × the peak deformation and the damping to the loop-area
damping, and the SSI analysis is repeated (New Structure restarts) until E changes by less than 2 % and
the damping by less than 0.5 %.

```
PNLGEN                                      * one panel per vertical shell group (each with its own material)
* shear capacities: f'c 30 MPa, f_y 420 MPa, 0.5 % web steel (kN/m2)
SHEAR,0,30000,420000,0.005,0
* 22-point backbones from the ACI 318-08 capacity, cracking at 0.3 V_u
BBCGEN,0,1,30000,420000,0.005,0,0,0,0.3
* EDF 0.8, panels, no damping cut-off, elastic damping added
EQL,0.8,1,0,0,1
NONLINMOTDISP                               * panel corners into the MOTION and RELDISP requests
...                                         * elastic SSI with restart files; MOTION + RELDISP per direction
* writes the generic script <model>_NONLINBAT.pre
NONLINBAT,1
* iterate until |dE/E| < 2 % and |d xi| < 0.5 %
NONLINITER,NLHOUSE+NLDX+NLDY+NLDZ+NLNON,10
```

Use smooth backbone curves and enough SSI frequencies around the structural frequencies, which fall as
the walls crack; a change that grows from one iteration to the next means divergence. The complete guide,
with the hysteresis rules and their sources, the capacity equations and a convergence discussion, is
**[OPTION_NON.md](OPTION_NON.md)**; the theory is in
[Theory Manual §17](../theory/THEORY_MANUAL.md#17-option-non-nonlinear-structures-by-equivalent-linearisation).
*Verification:* VP-45 (Masing loops, state machine), VP-46 (SHEAR/BBCGEN), VP-NON1 (SDOF through the
whole chain against a closed-form fixed-point iteration); example 7.

### 12.6 Water in pools and tanks

Water adds hydrodynamic mass to pool and tank walls. `FILLPOOL` fills a pool sub-model with water SOLIDs
(bulk modulus 2.2 GPa, shear modulus 10⁻⁸ K) attached to the walls by springs along the wall normals;
`REFINEMODEL` refines the wall mesh first, `LISTPOOLINTER` checks the interface nodes and `MERGEPOOL`
imports the water into the building model:

```
* the pool walls and floor (group 7) ...
CUTADD,1,7,RANGE,1,240
* ... copied to sub-model 5
CUT2SUB,1,5
ACTM,5
REFINEMODEL                      * optional: a finer wall mesh gives a finer water mesh
* springs 1e9 along the normals, 1 empty level, water nodes from 12001
FILLPOOL,1E9,0,1,-1,12000,0
ACTM,0
* every interface node found in the building
LISTPOOLINTER,5
* water and springs (and area shells) into the building
MERGEPOOL,5
FIXROT
```

The water SOLIDs reproduce the **impulsive** mass of Westergaard and Housner (VP-W1: within 1.5-1.9 % of
the exact potential flow); they do not reproduce sloshing (the convective mass), their response below
about 0.5 Hz has no meaning, and they need incompatible modes (FILLPOOL sets `MOPT,0`). Why the shear
modulus is 10⁻⁸ K and not ν = 0.49 (lead decision D-W3-07): with ν = 0.49 the water moves rigidly with the
tank below its shear modes and gives the total instead of the impulsive mass. The guide is
**[WATER.md](WATER.md)**; theory in
[Theory Manual §19](../theory/THEORY_MANUAL.md#19-water-hydrodynamic-mass-by-soft-solids).

### 12.7 Option A and Option AA

Option A, the two-step approach (SSI motions as loads of a refined ANSYS model, module LOADGEN), is
described in section 6.2 and in [OPTION_A.md](OPTION_A.md). Option AA (SSI with ANSYS dynamic matrices,
modules HOUSEFSA/ANALYSFSA and SSI2ANSYS) is not available.

### 12.8 Thick shells (TSHELL) and the node optimizer

* **TSHELL** (GROUP type 5) is a Mindlin-Reissner shell for walls and slabs whose thickness is not small
  compared with the span or the wavelength: transverse shear deformation is included and the MITC4
  assumed-strain field prevents shear locking, so the same element also works for thin plates. `EINT` 0
  (default) uses reduced bending integration with hourglass control, 1 selective (2 × 2) integration. A
  small drilling stiffness is added automatically: TSHELL nodes need no FIXROT. `THSHLSTR,1` adds face
  stresses to the STRESS output. Theory:
  [Theory Manual §13.2](../theory/THEORY_MANUAL.md#132-the-tshell-thick-shell-mindlin-reissner);
  *verification:* VP-38T (NAFEMS Test 21 thick plate within 2 % on 16 × 16 and 32 × 32 meshes), VP-TS1
  (patch tests, no locking), VP-54T (face stresses).
* **HOUSE node optimizer** (`HOUSEX,1`): reverse Cuthill-McKee renumbering started from the interaction
  nodes, which keep their bottom-up order; writes `<model>.hounew` and the `old new` pairs `<model>.map`.
  FILE4 and FILE8 then use the new numbers: as the manual warns, give the MOTION and RELDISP output nodes
  with the new numbers (look them up in the `.map`); ANALYS translates the FORCE loads and LOADGEN
  translates back to the model numbers. The results do not change, because SASSI-EDU's sparse solver
  orders the equations itself; the optimizer exists for fidelity with the manual and to renumber badly
  numbered models (VP-O1).

### 12.9 Frequency-split runs (AFWRBAT)

`AFWRBAT,<n>` splits the SSI frequency set into n contiguous blocks and writes, for block k, a folder
`<model>_<k>/` with the decks of SITE, POINT, HOUSE, FORCE and ANALYS for that block, a command file
`run_<k>.pre` and launchers `run_<k>.sh` / `run_<k>.bat`. The folders are self-contained and can run in
parallel or on other computers. In the model folder it writes the MOTION, STRESS and RELDISP decks and
`combine.pre` (`.sh`, `.bat`), which collects the FILE8 of every folder, merges them with COMBIN (also
FILE8X/Y/Z, load cases and incoherent samples) and runs the post-processors. Run EQUAKE and SOIL before
AFWRBAT; a FILE88 is copied into every folder.

```
AOPT,0,0,0,1,1,1,0,0,1,0,1,1,1,0
* four folders ex01_1 ... ex01_4 and combine.pre
AFWRBAT,4
```

Then run `run_1.sh` in folder `ex01_1` ... `run_4.sh` in `ex01_4` (in any order, at the same time if you
like; `.bat` on Windows) and finally `combine.sh` in the model folder.

### 12.10 What is not available

Binary result databases (`BINOUT` products and the `*DB` commands), Option AA, the vendor's New Load
Vector restart (ANALYS Mode 6), nonlinear beams in Option NON (`<NonLinOpts>` 4) and the Cheng-Mertz
bending model, the Abrahamson 1993 and 2006 coherency models (coefficients not available), anti-plane SH
analysis in 2D, sloshing of water, lower-bound incoherent spectra (`LBINCORS`: the manual gives no
algorithm), `SOILREDEF`, `GETENV`/`SETENV`, `FREAD`/`MREAD`, `CALCSECTHISTDB`, `BEAMPILE` and
`DCOUPLEBEAM`. The commands are accepted (stored for `WRITE`) and print a message; the
[Command Reference](../reference/COMMAND_REFERENCE.md) marks them. See also section 18.

---

## 13. Common errors and how to fix them

CHECK lists **errors** (the module's deck is not written) and **warnings** per module, with the
numbers of manual Chapter 10. SASSI-EDU adds checks numbered EDU-nn so the manual's numbering stays
faithful. The full catalogue is in `docs/spec/11_commands_water_nonlin_binary_tshell_errors_refs.md`
section 5. The ones you will meet most often:

| Message | Cause | Fix |
|---|---|---|
| Error 42 No Nodes Defined / 43 No Groups Defined | empty model (CHECK stops) | define nodes and groups |
| Error 13 Material `<m>` Is not Defined / 19 Soil Layer Is not Defined | an element refers to a missing M or L entry; for excavated solids MSET is an **L** number | `MSET`, `MACT`, `ETYPE` |
| Error 41 Node Is Not Defined | element node missing | `N` |
| Error 9 Nodes ... Are Collinear | collinear beam K node, degenerate shell | check the nodes |
| Error 10 Improper Release Code | the same component released at both beam ends | `KI`/`KJ` |
| Error 19 (SITE) half-space layer `<hs>` is not defined | SITE `<hs>` not an L number while `<nl>` > 0 | set `<hs>` |
| Error 44 Frequency Set Is Not Defined | SITE `<freq>` names a missing set | `FREQ` |
| Error 47 Illegal Number of Layers for Halfspace Simulation | `<nl>` not 0 or 4-20 | SITE `<nl>` |
| Error 55 Illegal Sum of Wave Ratios | WAVE ratios do not sum to 1 | `WAVE` |
| Error 57 Illegal Radius of Central Zone | POINT `<rad>` ≤ 0 | POINT, section 8.4 |
| Error 64 No Nodal Output Request | MOTION/RELDISP enabled without NOUT/RDND | `NOUT`, `RDND` |
| Error 73 Acceleration Time History File Does Not Exist | THFILE path wrong (a relative path is looked for in the model directory, then the working directory) | `THFILE`, `MDL`, run with `--cwd` |
| Errors 77 / 78 Multiplication Factor and Maximum Value ... | both zero or both non-zero | MOTION `<mult>`/`<max>` |
| Error 124 Node Is a Fixed Interaction Node | interaction node with fixed translations | `D` or `INT` |
| EDU-01 Interaction Node Is Not on a Soil Layer Interface | node elevation does not match a TOPL interface | adjust layers or mesh levels |
| EDU-21 Interaction Nodes Not Numbered Bottom-Up (warning; error for incoherent analyses) | embedded model numbered top-down | renumber the nodes, or let the HOUSE optimizer do it (`HOUSEX,1`) |
| EDU-23 Duplicate Frequency Numbers | the same number twice in FREQ | `FREQ,<set>,0` and redefine |
| EDU-03 Rounded Number of Fourier Components Is Shorter Than the Records | NFFT not a power of 2 and rounded below the record length | larger NFFT |
| Warning 9 Number of Values for Fourier Transform Is Not Power of 2 | NFFT rounded to the nearest power of 2 | set a power of 2 |
| EDU-06 Unrestrained Shell Drilling Rotation | shell nodes without drilling stiffness | `FIXROT` / `FIXSHLROT` |
| EDU-10 Passing Frequency Below the Cut-Off Frequency | layer or excavation element too thick | split the layers |
| EDU-11 Poisson Ratio > 0.47 | near-incompressible soil layer | check Vp; watch for unstable frequencies |
| EDU-12 Non-FV Method Selected | FI or FFV in use | validate against FV |
| EDU-19 Quiet Zone Too Short | record fills the Fourier period | larger NFFT |
| EDU-26 Option Combination or Value Not Allowed | for example incoherency in a 2D model or with SYMM planes, coherency model 2-7 or multiple excitation without wave passage, a FILE3 of the wrong dimension (POINT2/POINT3), TSHELL EINT not 0/1 | full 3D model; HOUSE `<wpass>` = 1 (`WPASS,1e9,...` = no delay); re-run POINT after changing HOUSE `<dim>` |
| coherency model 2 / 4: coefficients not available | the Abrahamson 1993 and 2006 coefficients could not be transcribed from a primary source | coherency model 3, 5, 6 or 7 (section 12.1) |
| Errors 58-60, 113-119 (incoherency, wave passage, multiple excitation) | Luco-Wong parameter < 0.1, mean Vs ≤ 0, embedded levels, apparent velocity ≤ 0, model number, ME zones whose first/last node is not an interaction node, amplification ratios out of [0, 10] or not one per SSI frequency | `INCOH`, `WPASS`, `ME`, `AMP` |
| Errors 121-128, EDU-44 (Option NON) | no panels or springs, a missing material or group, a force or displacement option this version does not support, further `.eql` findings | `P`, `S`, `BBC`, `GROUPMAT`; `EDUOPT,NONEXT,1` unlocks Takeda and bending panels (experimental) |
| Error 125 Not enough nonlinear soil properties | SOIL-NON (`NLSOIL,1`) without an NLSLAYER set for every soil sublayer | `NLSLAYER` for each SPRO sublayer above the half-space |
| ANALYS: the free field does not have the declared symmetry | SYMM types that do not match the input direction | X input: antisymmetric (1) about x = const, symmetric (0) about y = const; Z input: symmetric about both |
| W1 Gap Found at Node / W4 Unused Node | missing node numbers, unused nodes | harmless (fixed in the analysis files); clean with `RMVUNUSED`/`NCOM` |

**Runtime messages.** A module that stops prints the reason in its listing and on the screen, for
example *FILE1 missing (run SITE Mode 2)*, *Frequency n not in FILE1/FILE3*, *time history longer
than Fourier period*, *DF of the control motion differs from FILE8* (NFFT or Δt changed after
ANALYS). After any change, run `AFWRITE` again before `RUN<MODULE>`.

**Results that look wrong.**

* *Transfer functions far from 1 at the first frequency*: interaction nodes missing or not connected
  to the structure, wrong `<cm>`/`<ang>`, a mechanism (unrestrained DOFs).
* *Sharp spikes in the interpolated TF*: missing SSI frequencies (section 9.4) or an isolated
  unstable frequency (ν > 0.47); compare `.TFU` and `.TFI`.
* *Response that does not decay at the end of the history*: quiet zone too short (increase NFFT).
* *Spurious peaks with FI-FSIN*: the subtraction-method anomaly; use FI-EVBN, FFV or FV.

---

## 14. Engineering guidance

A paraphrase of the "Engineering Considerations" of manual §4.1.2, with the numbers you need.

1. **Basement.** Model the basement with its real (flexible) properties: a rigid basement saves no
   computation. To study the rigid limit, run New Structure restarts with the basement modulus × 10⁴
   and × 10⁵.
2. **Embedment.** Treating an embedded structure as a surface structure saves a lot of run time but
   is not conservative at every frequency: forces are usually conservative, ISRS can show new peaks.
3. **2D or 3D.** Prefer 3D; 2D SSI overestimates radiation damping. 2D is for sensitivity studies
   (section 12.2).
4. **Incoherency.** Only for 3D models without symmetry planes; stochastic simulation is the
   reference approach; deterministic approaches are approximate and validated for stick models on
   rigid mats only (section 12.1).
5. **Several structures (SSSI).** Put them, with their excavations, in one model; SSSI affects
   mainly ISRS and wall pressures.
6. **Non-uniform input.** Foundation zones with their own amplification (multiple excitation).
7. **Symmetry.** Half and quarter models are valid for symmetric structures under symmetric or
   antisymmetric loading, but not with incoherency or the impedance options (section 12.3).
8. **Half-space.** 10-20 generated layers (20 recommended); 0 means a rigid base.
9. **Cut-off frequency.** Limited by the input content, the system frequencies and the time step;
   30-40 Hz for soil sites, 60-70 Hz for rock sites; it limits the element sizes.
10. **SSI frequencies.** 40-80 for sticks, 100-200 (up to 300) for complex FE models, at least 200
    for incoherent analyses; refine where TFI peaks are not supported by TFU points.
11. **Frequency numbers.** `f = n Δf`, `Δf = 1/(Δt NFFT)`; at most 500 per run.
12. **Soil layers.** h ≤ Vs/(5 f_cut); more than 20 layers; beware ν > 0.47 and very soft deep
    deposits (isolated unstable frequencies: remove them). Model backfill as near-field soil, or use
    springs between duplicate perimeter nodes (INTGEN generates them) to estimate wall pressures.
13. **Structure.** Interaction nodes on the excavated soil only; below grade on layer interfaces;
    ascending (bottom-up) numbering; run EXCSTRCHK, FIXEDINT and HINGED.
14. **Excavated soil.** Solids (3D) or planes (2D) sized by the interaction-node spacing; keep
    internal basement structures separate from the excavation nodes; results are sensitive to the
    excavation mesh.
15. **Method.** FV is the reference; validate SM (FI-FSIN), MSM (FI-EVBN) and FFV against FV before
    production use (ASCE 4-16, SRP 3.7.2). Quarter models can hide instabilities of SM/MSM.
16. **Adding frequencies.** Possible as long as they exist in FILE1 (or FILE9) and FILE3; merge with
    COMBIN.
17. **Impedance.** `K + iωD = (f + ig)⁻¹`, a dense matrix: memory `(3 N_int)² × 16` bytes.
18. **Soil nonlinearity.** Primary (free field, SOIL) and secondary (near field, SSI iterations);
    the near field needs explicit soil elements and New Structure restarts (section 12.4). Nonlinear
    walls and springs of the structure are iterated the same way by Option NON (section 12.5).
19. **Mesh size.** Vertical ≤ Vs/(5 f_cut) (valid with the ½ lumped + ½ consistent mass);
    horizontal ≈ vertical unless a sensitivity study justifies 1.2-2 times more.
20. **Half-space simulation.** Variable depth (1.5 λ) plus viscous boundary, both in SITE.

---

## 15. The tutorial examples

The examples in `examples/` are ordinary `.pre` files, commented line by line. Run them from the
`examples` directory (`../.venv/bin/sassi run exNN_....pre`) or with `--cwd examples` from the
project root; the results go to `examples/exNN/`. [examples/README.md](../../examples/README.md)
lists the expected results, which `tests/integration/test_examples.py` (examples 1-5 and 8),
`tests/unit/test_nlsoil_example.py` (example 6) and `tests/unit/test_nonlinear_example.py` (example 7)
check.

| Example | Topic | Modules | Run time |
|---|---|---|---|
| ex01 | stick on a rigid surface mat: ISRS, relative displacement, beam forces | SITE POINT HOUSE ANALYS MOTION STRESS RELDISP | ~4 s |
| ex02 | embedded basement: FV, FI-FSIN and FI-EVBN | SITE POINT HOUSE ANALYS MOTION STRESS (3 models) | ~5 s |
| ex03 | forced vibration of a rigid surface foundation: compliance and impedance | SITE POINT HOUSE FORCE ANALYS MOTION | ~1 s |
| ex04 | free-field chain EQUAKE → SOIL → SITE | EQUAKE SOIL SITE | ~15 s |
| ex05 | X, Y and Z input in one ANALYS run | SITE ×3 POINT HOUSE ANALYS MOTION ×3 RELDISP | ~4 s |
| ex06 | loose backfill behind a wall: near-field soil iterations | SOIL SITE ×3 POINT HOUSE ANALYS STRESS ×3, then 6 iterations HOUSE ANALYS STRESS ×3; MOTION | ~5 s |
| ex07 | Option NON: shear-wall building with cracking wall panels | SITE ×3 POINT HOUSE ANALYS MOTION ×3 RELDISP ×3 NONLINEAR, then 7 iterations; MOTION | 30-120 s |
| ex08 | embedded shear-wall building: FV against FI-FSIN (SM) and FI-EVBN (MSM) | SITE POINT HOUSE ANALYS MOTION STRESS (3 models) | ~60 s |

### 15.1 Example 1: stick on a surface mat (the basic workflow)

*Model.* A 20 m × 20 m SHELL mat (9 × 9 nodes, 1.5 m thick, ten times stiffer than concrete) carries
a four-storey BEAMS stick (5 m storeys, 1000 t per floor given as weights with `MT` in a `FOREACH`
loop). A rigid "spider" of beams connects the stick base to the surrounding mat nodes. The site is
5 m of sand (Vs 300 m/s) and 12 m of gravel (Vs 500 m/s) on rock (Vs 1000 m/s), discretised in 22
TOPL layers of 0.5-1 m. Input: vertically incident SV wave, control motion in X at the free surface;
22 SSI frequencies from 0.1 Hz to the 20 Hz cut-off.

*Walk-through.* `L`/`TOPL` build the site; `FREQ` gives frequency numbers (Δf = 1/(8192 × 0.005) =
0.0244 Hz); `SITE`/`WAVE` set the wave field; `N`/`FILL`/`NGEN` and `E`/`EGEN` build the mesh;
`FIXROT` restrains the shell drilling rotations of shell-only nodes; `INT,1,81,1,1` makes the 81 mat
nodes interaction nodes; `POINT,0,0,2.25` sets R0 = 0.9 × 2.5 m; `HOUSE`, `ANALYS`, `MOTION`, `NOUT`,
`STRESS`, `EOUT`, `RELD`/`RELFILE`/`RDND` set the analyses; `AOPT`, `CHECK`, `AFWRITE` and the `RUN`
commands run them; `WRITE` saves the complete model.

*Results to look at.* The fixed-base frequency is 5.05 Hz; with SSI the roof transfer function peaks
at 3.49 Hz (|ATF| = 13.1): the soil springs lower the frequency and radiation damping limits the
peak. At 0.1 Hz every transfer function is 1.000. The ZPA grows from 0.36 g at the mat to 1.29 g at
the roof (input PGA 0.324 g). The roof drift relative to the mat is 23.2 mm. The top-storey shear
equals the roof mass times the roof acceleration within 2 %. The rocking motions of the two mat
edges are equal and opposite. VP-E1 runs this example and checks that AFWRITE writes exactly the
intended decks.

### 15.2 Example 2: embedded basement, FV versus FI-FSIN and FI-EVBN

*Model.* A 10 m × 10 m × 5 m concrete box (SHELL walls and slabs) embedded in 5 m of sand. The
excavated soil is 80 SOLID elements with `ETYPE` 2 in five groups, one per 1 m embedment layer; each
embedment layer has its own L number and each excavated group uses it (`MACT,k`). Interaction nodes
are set level by level with `INT` in `FOREACH` loops. `CPMODEL`, `ACTM` and `MDL` create two more
models whose interaction sets are FI-FSIN (lateral and bottom faces) and FI-EVBN (plus the top
face); `FCOPY` reuses FILE1 and FILE3 of the first run.

*Results to look at.* With FV the roof X transfer function stays at or below 1 (kinematic
interaction filters the motion). FI-FSIN shows a spurious resonance at 6.0 Hz (roof amplitude +58 % above FV): the
roof slab rests on top-face nodes that are not interaction nodes, and the slab and the subtracted
soil form an unphysical oscillator. FI-EVBN follows FV within 2.5 % in amplitude at a fraction of the
interaction nodes. This is the subtraction-method problem reported by DOE in 2011, in miniature.

### 15.3 Example 3: forced vibration, compliance and impedance

*Model.* A rigid massless 12 m × 12 m mat (7 × 7 interaction nodes) on a uniform half-space, tied to
its centre node by stiff beams. Six unit load cases (`F`/`MM` at the centre) are run through FORCE
and copied to FILE9001 ... FILE9006; one ANALYS run (`<type>` = 1, `<simul>` = 6, `<impe>` = 2) solves
all six and writes the global impedance (FOUNSTIF, FOUNDASH, FOUNDAMP, FOUNIMPD). MOTION
(`EDUOPT,TFFILE,FILE8003`, displacement response) gives the vertical displacement under a 1000 kN
Ricker pulse.

*Results to look at.* The inverse of the 6×6 compliance from FILE8001-8006 equals the global
impedance within 0.2 %. Compared with Pais and Kausel (1988) the 2 m mesh is 6 % too stiff in
translation and about 20 % in rocking and torsion: the impedance converges with mesh refinement
(VP-14 shows the convergence). The Ricker pulse response peaks at 0.205 mm, below the static value
0.279 mm, and radiation damping stops the motion almost at once.

### 15.4 Example 4: the free-field chain EQUAKE → SOIL → SITE

*Model.* 10 m of sand on 12 m of clay on rock, in 1 m sublayers, no structure. EQUAKE matches a 20 s
record to the RG 1.60 spectrum (`RSIN`/`RSOUT`/`ACCOUT`); SOIL runs the equivalent-linear analysis
with the record as a rock-outcrop motion, using the SHAKE91 curves of `sassi/data/dynp_library.pre`
assigned with `SPRO` in `FOREACH` loops; `SITEX,1` makes SITE read FILE88.

*Results to look at.* EQUAKE: every SRP 3.7.1 acceptance check passes (listing `ex04_EQUAKE.out`).
SOIL: the sand softens from Vs 200 m/s to 74 m/s at 10 m depth, damping 2 % to 18 %; the
surface/base amplification peaks at 4.65 at 1.83 Hz (68 at 3.10 Hz with the low-strain properties).
SITE with FILE88 reproduces the SOIL amplification within 0.5 % up to 7 Hz; the difference grows as
(kh)², the thin-layer discretisation error, which tells you how thin the layers must be.

### 15.5 Example 5: X, Y and Z input with simultaneous cases

*Model.* A two-storey stick on a 12 m × 12 m mat with an eccentric roof mass. SITE runs three times
(SV along x', SH along y', P along z'; FILE1 copied to FILE1X/Y/Z); ANALYS with `<simul>` = 1 solves the
three cases with one factorisation per frequency (FILE8X/Y/Z); MOTION runs once per direction
(`EDUOPT,TFFILE,FILE8X` ...; vertical motion scaled to 2/3), RELDISP for Y with the free-field
reference.

*Results to look at.* Each case has |ATF| = 1 in its own direction at low frequency. The eccentric
mass couples the directions: the X input produces a Y response of the roof mass of up to 1.85 and a
torsion of the stick. The directional results are combined afterwards (SRSS or 100-40-40), outside
the SASSI modules.

### 15.6 Example 6: near-field soil nonlinearity

*Model.* A 0.6 m reinforced-concrete wall (SHELL, 4 m long, 2 m high) carrying a 150 t bridge deck,
with a 4 m × 4 m × 2 m block of loose backfill (Vs 160 m/s at low strain) behind it, in 10 m of
medium-dense sand on rock. SOIL first computes the strain-compatible free field for a 0.30 g rock
outcrop motion (0.58 g at the surface). The backfill is SOLID elements of the *structure* (ETYPE 1)
at the nodes of the excavated soil (FV); `PIN`/`PINGRP` declare it a nonlinear soil group and
`NLSSIITER` runs the iterations HOUSE → ANALYS (New Structure restart) → STRESS (X, Y, Z) →
COMBXYZSTRAIN.

*Walk-through.* Part 1 runs SOIL alone (`SPRO` with the SHAKE91 Sand curve per sublayer, rock outcrop
scaled to 0.30 g with `SOILX`, `SACC` saves the surface motion `ACC001.TH`). Part 2 copies the
strain-compatible layers of FILE88 into L 11-16 (each embedment layer with its own L number), builds
the excavated soil (two groups, ETYPE 2, `MACT` = the layer), the backfill (group 3, ETYPE 1, a
low-strain material `M,2` with Vs 160 m/s), the wall (SHELL) and the deck masses, and makes every
excavated node an interaction node (FV). Part 3 declares the nonlinear group (`PIN,0.65`,
`PINGRP,3,1,0.4,1.0,Sand`, `HOUSEX,0,0,1,1`). Part 4 is the initiation: three SITE runs (FILE1X/Y/Z),
POINT, HOUSE, ANALYS with restart files, then STRESS (`<iter>` = 1) per direction and
`COMBXYZSTRAIN`, written as command lists in variables (`VAR`) and run with `FOREACH`. Part 5 switches
ANALYS to New Structure (`<mode>` 1) and runs `NLSSIITER` with the same command lists. Part 6 computes
the spectra of the wall top with the converged backfill.

*Results to look at.* The iterations converge in six passes (max |ΔG/G| = 66, 22, 14, 9.4, 4.0, 2.9,
1.5 %); the converged backfill has G/Gmax = 0.52-0.57 and 6.4-7.2 % damping, and the elements next
to the wall strain more than those at the back. The wall top reaches 0.75 g against 0.58 g in the
free field. See `ex06_HOUSE.out` (convergence tables), `NLSOIL_CONVERGENCE.TXT`, FILE74 and FILE78.
The expected values are tested by `tests/unit/test_nlsoil_example.py`.

### 15.7 Example 7: nonlinear shear walls (Option NON)

*Model.* A 12 m × 12 m two-storey reinforced-concrete box (storeys of 4 m, 0.3 m walls meshed 3 m × 2
m, 0.4 m slabs with 400 t of equipment each) on a 1.5 m stiff surface mat, on 10 m of stiff soil (Vs
400 m/s, 25 sublayers) over rock. Each wall of each storey is a SHELL group with its own material: 8
panels (h/l = 1/3, low-rise, shear governed). The RG 1.60 record of the examples is scaled to 0.6 g in X
and Y and 0.4 g in Z (Demo 9 of the manual is an RC shear-wall building at 0.60 g).

*Walk-through.* `PNLGEN` makes one panel per vertical shell group; `SHEAR` lists the ACI 318-08, Wood,
Barda and Gulec-Whittaker capacities; `BBCGEN,0,1,...,0.3` builds the 22-point backbone of every panel
from the ACI capacity with cracking at 0.3 V_u; `EQL,0.8,1,0,0,1` sets EDF 0.8 and adds the elastic 4 %
damping; `NONLINMOTDISP` adds the panel corners to the MOTION and RELDISP requests. The elastic SSI
analysis (three SITE runs, POINT, HOUSE, ANALYS with restart files), MOTION + RELDISP per direction,
`COMBXYZTHD` and the elastic NONLINEAR run follow; `NONLINITER` then repeats HOUSE (with
`ex07_new.hou`) → ANALYS New Structure → MOTION + RELDISP ×3 → COMBXYZTHD → NONLINEAR.

*Results to look at.* Convergence after 7 iterations (max |ΔE/E| = 28.0, 22.9, 17.1, 11.5, 6.6, 4.2,
2.4, 1.3 %). The four storey-1 walls crack: E/E_el = 0.35 and 16.7 % damping (4 % elastic + 12.7 %
hysteretic), peak shear strain 3.7e-4, i.e. 4.4 times the cracking strain (F_μ = 3.3); the storey-2
walls stay elastic (μ = 0.84). With the converged model the X acceleration is 0.64 g at the mat centre,
0.85 g at the floor and 1.12 g at the roof. See `ex07_NONLINEAR.out`, `NONLINEAR_CONVERGENCE.TXT`,
`Panel_EQL_Matl_Prop.txt`, `Panel.fmu` and `Panel000k.thd` / `.ths` (plot one against the other to see
the hysteresis loops). [OPTION_NON.md §11](OPTION_NON.md#11-tutorial-example-7) discusses the variants
(a 7 % damping cut-off, the Gulec-Whittaker backbone). Tested by `tests/unit/test_nonlinear_example.py`.

### 15.8 Example 8: an embedded shear-wall building, FV against the subtraction methods

*Model.* A 24 m × 24 m reinforced-concrete shear-wall building embedded 8 m: a 2 m basemat, 1.0 m outer
basement walls, 0.6 m interior walls on the axes x = 0 and y = 0, a 0.6 m basement slab and a 0.8 m grade
slab, three storeys of 5 m above grade (0.8 m outer walls, 0.6 m floors, 0.5 m roof), 830 t of equipment
(`MT`); 15,063 t in all. Site: 8 m of sand and gravel (Vs 300 m/s) over dense gravels and rock. The
excavated soil (256 SOLID elements, 3 m × 3 m × 2 m: 20 Hz by the λ/5 rule) shares nodes with the
basemat and the outer walls only; the interior walls and slabs have nodes of their own (manual rule 11,
`EXCSTRCHK`), generated as "excavation node + 1000" and renumbered by `RMVUNUSED` and `NCOM`. `INTGEN`
builds the FV (405 nodes), FI-FSIN (209) and FI-EVBN (258) sets on copies of the model (`CPMODEL`,
`FCOPY` of FILE1 and FILE3).

*Results to look at.* FV: roof transfer function 2.44 at 6.5 Hz, 5 % ISRS ZPA 0.277 g at the basemat,
0.322 g at grade and 0.420 g at the roof. FI-FSIN follows FV within 2 % below 10 Hz but resonates at
16 Hz (roof 0.88 against 0.32, the soil enclosed by the basement 12.4 against 1.8) with a dip at 15 Hz:
the excavated soil inside the interface, which carries no impedance, has its own natural frequency
there. Its ISRS differ by up to +9.5 % / −8 % at 14-16 Hz and the wall forces by less than 2 %, because
the resonance is narrow and far above the 6.5 Hz SSI mode. FI-EVBN follows FV within 2.4 % (ISRS 0.7 %).
On this small model the reduced sets save little time (ANALYS 18 s FV, 13 s FI-FSIN, 15 s FI-EVBN): the
structure's equations cost as much as the impedance. Lesson 11 of the guided course builds the example
step by step, animates the spurious resonance, runs FFV and discusses the validation against FV.

---

## 16. Verification

Every physics module is verified by verification problems (VPs) that compare computed values with
exact, published or independently derived references and state the tolerance. The
[Verification Manual](../verification/VERIFICATION_MANUAL.md) shows, for each VP, the model, the
computed and reference values, the errors, the tolerances and the result.

```bash
.venv/bin/python -m pytest -q tests/verification            # all VPs as pytest cases
.venv/bin/python -m sassi.verify.report                     # regenerate the Verification Manual
.venv/bin/python -m sassi.verify.report --only VP-01,VP-30   # partial report: ./VERIFICATION_REPORT_partial.md
```

In the console or the GUI command line:

```
* list the problems
VERIFY,LIST
* run one problem and print its table
VERIFY,VP-15
* all P0 problems (default)
VERIFY,P0
* run the problems that are not slow and write a Markdown report
VERIFYREPORT,FAST,my_report.md
```

Some published references proved to be approximate (for example the Veletsos-Verbic rocking
damping); the lead decisions of requirements sections 7.16 to 7.18 replaced them by rigorous
references and keep the published values as *informative* comparisons. The Verification Manual
explains each case. The problems of the later features are listed in requirements §6.3: incoherency
VP-26, VP-27, VP-I1, VP-I2; 2D and symmetry VP-42, VP-T1, VP-T2, VP-T3; nonlinear soil VP-N1, VP-N2;
TSHELL and the optimizer VP-38T, VP-54T, VP-TS1, VP-O1; Option NON VP-45, VP-46, VP-NON1; Option A
VP-LA1, VP-LA2; SOIL-NON VP-SN1, VP-SN2, VP-SN3; water VP-54W, VP-W1.

---

## 17. ACS SASSI terms and their SASSI-EDU equivalents

| ACS SASSI term | SASSI-EDU command / file | Notes |
|---|---|---|
| UI (User Interface), PREP | `sassi` console, `sassi run`, `sassi-gui`; `sassi.prep.Interpreter` | one interpreter for all |
| Command Entry, Command History | console prompt; GUI Command Entry and Command History tab | GUI actions emit commands |
| Options > Analysis dialog tabs | option commands `EQUAKE`, `SOIL`, `SITE`, `POINT`, `HOUSE`, `FORCE`, `ANALYS`, `MOTION`, `STRESS`, `RELD` and the X-commands `SITEX`, `SOILX`, `HOUSEX`, `ANALYSX`, `MOTIONX`, `STRESSX`, `RELDX` | X-commands carry dialog fields without a manual argument |
| AFWRITE (write analysis files) | `AFWRITE` → `<model>.sit .poi .hou .anl .mot .str .rdi .soi .equ .frc .eql` (and `.pin`, `.nls`) | decks are text keyword files |
| CHECK, Check Errors window | `CHECK` → `<model>.err` | Errors 1-128, Warnings 1-11, EDU-nn |
| Run Module (Modules menu) | `RUNSITE`, `RUNPOINT`, `RUNHOUSE`, `RUNANALYS` ... | or `python -m sassi.modules.<name> < NAME.inp` |
| HOUSEFS / ANALYSFS (fast solvers) | HOUSE / ANALYS | one implementation (sparse + Schur complement) |
| Model database, numbered models | `ACTM`, `CPMODEL`, `DMODEL`, `MODELLIST`, `SAVE`/`RESUME` (`.sdb`) | |
| `.pre` file | `WRITE` / `INP` | exact round trip |
| Flexible Volume (FV), Flexible Interface (FI-FSIN, FI-EVBN), Fast Flexible Volume (FFV) | interaction-node set by `INT` / `INTGEN`; HOUSE `<imp>` | section 8.3 |
| Interaction node | `INT,...,1,0`; `INTLIST`, `INTCOUNT` | |
| Excavated soil | SOLID/PLANE with `ETYPE` 2 and `MSET` = L layer | |
| Radius of Central Zone | POINT `<rad>`; `RADIUS` | R0 = 0.9 h (square), 0.85 h (triangular), h (2D) |
| Number of generated layers / halfspace layer | SITE `<nl>`, `<hs>`; `EDUOPT,HSLAW` | |
| Control point layer / direction | SITE `<cl>`, `<cm>` | within motion |
| Frequency set, frequency numbers | `FREQ`, `LFREQ`; SITE `<freq>` | `f = n Δf` |
| Simultaneous Cases | ANALYS `<simul>`; FILE1X/Y/Z, FILE9001 ... | `FCOPY` to create them |
| Save Restart Files; New Structure / New Seismic Environment / New Dynamic Loading | ANALYS `<save>`, `<mode>` = 1/2/3; COOXqqq, COOTKqqq, COOXI, COOTKI, DOFSMAP, FILE90, FILE91 | |
| Global (unconstrained) impedance | ANALYS `<impe>`; FOUNSTIF, FOUNDASH, FOUNDAMP, FOUNIMPD, FILE11 | |
| Interpolation Option 0-6, Smoothing, Phase Adjustment | MOTION `<interp>`, `<smo>`, `<pzadj>`; STRESS `<interopt>`, STRESSX | |
| Save Complex TF (.TFI) | MOTION `<cplx>` = 1 | needed by RELDISP |
| `Remove_Frequencies_from_FILE8.exe` | `REMOVEFREQ` | |
| `COMBIN_XYZ_STRAIN.exe` | `COMBXYZSTRAIN` | SRSS of the X, Y, Z effective strains (FILE74) |
| `Build_FILE77.exe` | `BUILDFILE77` | per-level incoherency factors combined into one FILE77 |
| `COMB_XYZ_THD` | `COMBXYZTHD`; `python -m sassi.modules.nonlinear COMB_XYZ_THD COMB_XYZ_THD.inp` | X + Y + Z relative-displacement histories for Option NON |
| Unlagged coherency model, Line D, apparent velocity | `WPASS,<appv>,<ang>,<cohf>`; `INCOH` | section 12.1; coefficients in `sassi/data/coherency/` |
| Incoherent modes, Simulation Mean, Superposition Mode (Linear / Quadratic), Number of Simulations | `INCOH <nmodes>`, `<HSeed>`, `<VSeed>`, `<RandPhz>`; `HOUSEX <supmode>`, `<nsim>`; FILE77, FILE77001 ... | stochastic simulation, AS, SRSS TF / SRSS FRS |
| Free-Field Load / Free-Field Motion | `ANALYSX,<ffm>` | FFL default |
| Multiple excitation, spectral amplification ratios | `ME`, `AMP`; HOUSE `<me>`, `<cmplxspec>` | |
| Symmetry of the System | `SYMM,<no>,<type>,<n1>,<n2>,<n3>` | half and quarter models, section 12.3 |
| 2D analysis (POINT2) | HOUSE `<dim>` = 1 | plane strain in the X-Z plane, section 12.2 |
| Optimize Model | `HOUSEX,1`; `<model>.hounew`, `<model>.map` | section 12.8 |
| Non-Linear SSI (HOUSE), `.pin`, `.liq` | `PIN`, `PINGRP`, `PINMAT`, `HOUSEX,0,0,1,1`, `NLSSIRESET`, `NLSSIITER` | near-field soil iterations, example 6 |
| SOIL-NON | `NLSOIL`, `NLSLAYER`, `DELNLS`; side file `<model>.nls` | section 12.4 |
| Option NON, NONLINEAR ("PANEL") module | `EQL`, `P`, `S`, `BBC`, `BBCX`/`BBCY`/`BBCP`, `SHEAR`, `BBCGEN`, `PNLGEN`, `NONLINMOTDISP`, `RUNNONLINEAR`, `NONLINITER`, `NONLINBAT`; `.eql` | [OPTION_NON.md](OPTION_NON.md), example 7 |
| Option A, LOADGEN, ANSYS Eq. Static Load / ANSYS Dynamic Load | `LOADGEN`, `LOADGENDYN`, `LGFILE`, `LGNODE`, `LGTIME`, `LGMAP`, `LGOPT`, `RUNLOADGEN`; `.lgn` | [OPTION_A.md](OPTION_A.md) |
| Water Modeling commands | `FILLPOOL`, `REFINEMODEL`, `LISTPOOLINTER` (and `MERGEPOOL`, `POOLDATA` of SASSI-EDU) | [WATER.md](WATER.md) |
| TSHELL, thick-shell stresses | GROUP type 5, `EINT`, `THSHLSTR` | section 12.8 |
| AFWRBAT (frequency-split batch run) | `AFWRBAT,<n>`; `run_<k>.pre`, `combine.pre` | section 12.9 |
| file copy/rename in the model directory | `FCOPY`, `FMOVE` | SASSI-EDU extension |
| complex modulus form | `CMODFORM` | 0 SASSI form (default), 1 = 1 + 2iβ |
| algorithm switches | `EDUOPT,<key>,<value>` | requirements §7 decisions |
| V&V Manual problems | `VERIFY`, `VERIFYREPORT`, `python -m sassi.verify.report` | Verification Manual |
| Demo 1-12 (vendor examples) | `examples/ex01` ... `ex08` | similar workflows (ex07 as Demo 9; ex02 and ex08 as the embedded Demo 5) |
| Binary files FILE1 ... FILE91 | same names, NumPy `.npz` containers without extension | readable with `sassi.io.container` |

---

## 18. Limitations

SASSI-EDU reproduces the methodology and the nomenclature of ACS SASSI, not its code. Main
differences and missing features in this build:

* **Not available:** binary result databases (BINOUT products and the `*DB` commands), Option AA
  (SSI with ANSYS matrices), the vendor's New Load Vector restart (ANALYS Mode 6), the coherency models
  Abrahamson 1993 and 2006 (no verifiable coefficients), anti-plane SH analysis in 2D, nonlinear beams
  and the Cheng-Mertz bending model in Option NON, sloshing of water, lower-bound incoherent spectra
  (LBINCORS) and a few utilities (section 12.10). The commands are accepted and stored or print
  *not available in this build*.
* **Restrictions of the available options.** Incoherency, wave passage and multiple excitation need a
  full 3D model (no SYMM, no 2D); 2D analysis is plane strain in the X-Z plane with in-plane input; the
  global impedance is 3D only and not available with SYMM; FFM is meant for surface foundations; Takeda
  and bending panels of Option NON are experimental (`EDUOPT,NONEXT,1`); the water SOLIDs give the
  impulsive mass only.
* **Reconstructions.** Where the manual is silent the implementation follows published theory and
  recorded decisions (requirements §7): the half-space sublayer law, the central-zone discretisation
  (one radial element), the weighting of interpolation option 0, the window shifts of options 4/5, the
  inclined-wave half-space treatment, the incompatible-mode formulations, the sign and basis of the
  incoherent modes, the CMS and Takeda hysteresis rules (from Cheng and Mertz's INRESB-3D-SUP listing),
  the SOIL-NON damping and time-stepping details (DEEPSOIL literature), the water shear modulus, and
  the LOADGEN deck and APDL layout (the vendor's ANSYS integration manual was not available; the APDL
  has been checked by a parser, not run in ANSYS). Results can differ from ACS SASSI in the details.
* **Files.** Binary files are NumPy containers, not the vendor's formats; the `.pre` language is that
  of the manual, so `.pre` files can be exchanged, but legacy fixed-format decks cannot be read.
* **Performance.** Pure Python and NumPy/SciPy; the target is about 2,000 interaction nodes and
  20,000 DOFs per frequency on a laptop.
* **Quality assurance.** No NQA-1 program; use the verification problems as a learning tool and to
  check changes, not as a qualification.

For a command-by-command list of what is implemented, see the
[Command Reference](../reference/COMMAND_REFERENCE.md).
