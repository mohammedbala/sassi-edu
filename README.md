# SASSI-EDU

An educational, open re-implementation of the **SASSI flexible-volume method** of soil-structure
interaction (SSI) analysis in Python. It follows the principles, module structure, command language,
option codes, file names and messages of the **ACS SASSI Version 3** user manual, so that what you
learn here carries over to the commercial code. It is written for structural engineers who know
finite elements (ANSYS) and want to learn frequency-domain SSI.

> **Disclaimer.** SASSI-EDU is a learning tool written from the published methodology (Lysmer et al.
> 1981; SASSI2000; Kausel's thin-layer method; SHAKE). It is **not affiliated with, endorsed by, or a
> substitute for ACS SASSI** or its vendor, and it is not qualified (10 CFR 50 Appendix B,
> ASME NQA-1) for licensing or design work.

## What it does

The same chain of modules as ACS SASSI (manual Fig. 1.1), each an independent program that reads an
input deck and the binary files of the modules before it:

| Module | Purpose | Reads | Writes |
|---|---|---|---|
| EQUAKE | spectrum-compatible acceleration histories with the SRP 3.7.1 checks | `.equ` | `.acc .vel .dis .rso .psd .fft` |
| SOIL | 1D free-field site response: equivalent linear (SHAKE) or nonlinear in the time domain (SOIL-NON: hyperbolic soil, Masing rules) | `.soi`, `.nls` | `ACCxxx.TH`, spectra, FILE73, FILE88 |
| SITE | layered free field by the thin-layer method: Rayleigh/Love modes, free-field motion (vertical and inclined body waves, surface waves) | `.sit` | FILE2, FILE1 |
| POINT | point-load flexibility of the layered site (POINT3 axisymmetric central zone; POINT2 plane strain) | `.poi`, FILE2 | FILE3 |
| HOUSE | finite-element matrices of the structure and the excavated soil (SOLID, BEAMS, SHELL, TSHELL, PLANE, SPRING, GENERAL), node optimizer, incoherency factors, near-field nonlinear soil | `.hou`, `.sit`, `.pin` | FILE4 (`<model>.N4`), COOSK, COOSM, FILE77, FILE78 |
| FORCE | external load vectors (vibration analysis) | `.frc` | FILE9 |
| ANALYS | impedance `X = F⁻¹`, flexible-volume equation solved by a Schur complement, restarts, simultaneous cases, global foundation impedance, incoherent input (FFL/FFM), 2D plane strain, symmetry planes | `.anl`, FILE1/FILE9, FILE3, FILE4, FILE77 | FILE8 |
| COMBIN | merge two FILE8 computed for different frequency sets | FILE81, FILE82 | FILE8 |
| MOTION | transfer-function interpolation (options 0-6), convolution, in-structure response spectra | `.mot`, FILE8 | `.TFU .TFI .ACC .RS` |
| RELDISP | relative displacements | `.rdi`, `.TFI` | `.TFD .THD` |
| STRESS | element stresses, forces and strains; near-field soil strains | `.str`, FILE4, FILE8 | `.TFU .TFI .THS`, FILE74 |
| NONLINEAR | Option NON: hysteresis of wall panels and springs, equivalent-linear properties for the next iteration | `.eql`, `.hou`, `.THD` | `<model>_new.hou`, `Panel*` files |
| LOADGEN | Option A: SSI motions as equivalent static or dynamic ANSYS loads (APDL) | `.lgn`, FILE4, `.ACC`, `.THD` or FILE8 | `<model>_LGS.inp`, `<model>_LGD.inp` |

Around the modules:

* the **ACS SASSI command language** (`N`, `NGEN`, `E`, `EGEN`, `GROUP`, `M`, `L`, `R`, `SC`, `SITE`,
  `HOUSE`, `ANALYS`, `MOTION`, `CHECK`, `AFWRITE`, `RUNANALYS`, variables, `FOREACH`, macros ...),
  from `.pre` files, the console or the browser GUI, with `WRITE`/`INP` round trips;
* **CHECK** with the manual's error and warning catalogue plus SASSI-EDU checks;
* model **generation and checking** tools (INTGEN, EXCAV, MERGESOIL, WELD, EXCSTRCHK, FIXEDINT,
  HINGED, FIXROT ...);
* **post-processing**: plots, line mathematics, spectrum broadening, section cuts;
* **ANSYS**: `.cdb` import and APDL export of models, and the **Option A** two-step approach (SSI
  motions as ANSYS loads, module LOADGEN);
* **incoherent motion, wave passage and multiple excitation**: Luco-Wong, Abrahamson 2005 and 2007
  and user coherency models, stochastic simulation, algebraic-sum and SRSS approaches, free-field load
  or motion;
* **2D plane-strain SSI** and **symmetry planes** (half and quarter models);
* **nonlinearity**: near-field soil SSI iterations (equivalent linear, New Structure restarts),
  **Option NON** for cracking shear walls and nonlinear springs, **SOIL-NON** time-domain site response;
* **water** in pools and tanks (FILLPOOL), the **TSHELL** thick shell, the HOUSE **node optimizer** and
  **AFWRBAT** frequency-split runs;
* **verification problems** for every module, run by pytest, by the `VERIFY` command and by the
  report generator that writes the Verification Manual.

## Quick start

```bash
python3 -m venv .venv && .venv/bin/pip install -e .[test]
.venv/bin/sassi --cwd examples run ex01_surface_stick.pre    # batch run of the first tutorial model
.venv/bin/sassi                                              # interactive command console
.venv/bin/sassi-gui                                          # browser GUI; opens on Learn, the guided course (10 lessons)
.venv/bin/python -m pytest -q                                # unit tests and verification problems
.venv/bin/python -m sassi.verify.report                      # regenerate the Verification Manual and the Command Reference
.venv/bin/python -m sassi.verify.report --commands-only      # regenerate the Command Reference only
```

**No installation:** the same GUI, course and examples also run entirely in a web browser, with Python
in the page (Pyodide); `python web/build.py` builds that static site for GitHub Pages
([web/README.md](web/README.md), [GUI.md §16](docs/user/GUI.md#16-web-version-github-pages)).

**New to SSI or SASSI?** Start the GUI: it opens on the **Learn** tab, a guided course of ten lessons
for structural engineers who design with fixed-base models (fundamentals, design applications,
advanced). Each step runs real commands next to their plots and listings and explains what they do,
why they matter for ISRS, member forces and displacements, their technical basis and the ANSYS
analogy; any command line can be explained argument by argument ([GUI.md §13](docs/user/GUI.md)).

The shape of a model file (`.pre`, or the same lines typed in the console):

```
* model name and directory
MDL,demo,demo
* SI units: m, kN, t, s
GRAVITY,9.81
* soil layer: thickness, weight, Vp, Vs, damping
L,1,1.0,19,500,250,0.05,0.05
* half-space
L,2,1.0,21,2000,1000,0.02,0.02
* ten 1 m layers on the half-space
TOPL,1,1,1,1,1,1,1,1,1,1
* SSI frequency numbers: f = n df, df = 1/(0.005 x 8192)
FREQ,1,4,20,41,82,123,164,205,246,328,410
* ... up to the 20 Hz cut-off
FREQ,1,492,573,655,737,819
SITE,0,1,0,20,2,1,0,1,4096,1,0,0.005,8192,1
* vertically incident SV wave
WAVE,2,1,1,1,0
...                                   * structure, interaction nodes, POINT, HOUSE, ANALYS, MOTION
AOPT,0,0,0,1,1,1,0,0,1,0,1,0,0,0
CHECK
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
RUNMOTION
```

The seven tutorial examples in [`examples/`](examples/README.md) are complete, commented models: a stick
on a surface mat, an embedded box (FV, FI-FSIN, FI-EVBN), forced vibration and foundation impedance,
the EQUAKE-SOIL-SITE chain, X+Y+Z input in one run, near-field soil nonlinearity (a loose backfill behind
a wall) and Option NON (a shear-wall building with cracking walls).

## Documentation

Start at [`docs/index.md`](docs/index.md).

| Document | Content |
|---|---|
| [User Guide](docs/user/USER_GUIDE.md) | how to build and run models, the advanced options, engineering guidance, the examples, ACS SASSI term mapping, limitations |
| [Theory Manual](docs/theory/THEORY_MANUAL.md) | the equations as implemented, with the code and the verification problems of each part |
| [Verification Manual](docs/verification/VERIFICATION_MANUAL.md) | computed versus reference values of every verification problem (generated) |
| [Command Reference](docs/reference/COMMAND_REFERENCE.md) | every command: syntax, abbreviation, tier, implemented or not (generated) |
| [GUI](docs/user/GUI.md), [ANSYS](docs/user/ANSYS.md) | the browser interface; the ANSYS converters |
| [Option A](docs/user/OPTION_A.md), [Option NON](docs/user/OPTION_NON.md), [Water](docs/user/WATER.md) | SSI loads for ANSYS (LOADGEN); nonlinear structures (NONLINEAR); water in pools and tanks (FILLPOOL) |
| [Architecture](docs/ARCHITECTURE.md), [requirements](docs/spec/00_requirements.md), [theory note R1](docs/spec/R1_sassi_theory.md), [benchmarks R2](docs/spec/R2_benchmarks.md) | specifications for developers |

## Limitations

* **Not available in this build:** binary result databases (BINOUT products and the `*DB` commands),
  Option AA (SSI with ANSYS matrices), the vendor's New Load Vector restart (ANALYS Mode 6), the
  Abrahamson 1993 and 2006 coherency models (their coefficients could not be verified from a primary
  source), anti-plane SH analysis in 2D, nonlinear beams in Option NON, water sloshing (the convective
  mass) and a few utilities. The commands are accepted (option commands are stored and written back)
  and print *not available in this build*. See [User Guide §18](docs/user/USER_GUIDE.md#18-limitations)
  and the [Command Reference](docs/reference/COMMAND_REFERENCE.md).
* Where the manual does not document an algorithm, SASSI-EDU follows the published theory and records
  its choice as a decision (requirements §7); results can differ from ACS SASSI in details. The LOADGEN
  input and APDL layouts are SASSI-EDU's own (the vendor's ANSYS integration manual was not available),
  and the APDL has been checked by a parser, not run in ANSYS.
* Binary files keep the legacy names (FILE1 ... FILE8) but are NumPy containers; legacy fixed-format
  SASSI2000/PREP decks cannot be read.
* Pure Python with NumPy/SciPy: the design target is about 2,000 interaction nodes and 20,000 DOFs per
  frequency on a laptop.

## Sources and data

The SHAKE91 sample data and soil curves in `sassi/data/` are public domain (see the README files
there). The coherency-model coefficients in `sassi/data/coherency/` were transcribed from the EPRI and
NRC reports cited in each file (report, equation, table and page). The method follows the open
literature listed in the [Theory Manual](docs/theory/THEORY_MANUAL.md#23-references).
