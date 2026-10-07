# SASSI-EDU tutorial examples

Nine small models that run the complete ACS SASSI workflow the way a user runs it: a `.pre` file
of ACS SASSI V3 commands builds the model, `CHECK` lists the errors and warnings of the model
and of the enabled modules, `AFWRITE` writes one input deck per module, and `RUN<MODULE>` runs the
modules one after the other in the model directory. Each file is commented line by line; read it
next to the results.

| file | topic | modules run | run time |
|---|---|---|---|
| `ex01_surface_stick.pre` | stick on a rigid surface mat, vertical SV (X): ISRS, relative displacement, beam forces | SITE POINT HOUSE ANALYS MOTION STRESS RELDISP | ~4 s |
| `ex02_embedded_box.pre` | embedded basement: FV, FI-FSIN and FI-EVBN interaction sets | SITE POINT HOUSE ANALYS MOTION STRESS (3 models) | ~5 s |
| `ex03_forced_vibration.pre` | forced vibration of a rigid surface foundation: compliance, impedance, global impedance | SITE POINT HOUSE FORCE (x6) ANALYS MOTION | ~1 s |
| `ex04_site_response.pre` | free-field chain: EQUAKE, SOIL (equivalent linear), SITEX 1 + SITE | EQUAKE SOIL SITE | ~15 s |
| `ex05_xyz_simultaneous.pre` | X, Y and Z input in one ANALYS run (simultaneous cases) | SITE (x3) POINT HOUSE ANALYS MOTION (x3) RELDISP | ~4 s |
| `ex06_nonlinear_soil.pre` | loose backfill behind an embedded wall: near-field soil nonlinearity by equivalent-linear SSI iterations | SOIL SITE (x3) POINT HOUSE ANALYS STRESS (x3), then 4 iterations of HOUSE ANALYS (New Structure) STRESS (x3); MOTION | ~5 s |
| `ex07_option_non.pre` | Option NON: two-storey shear-wall building with cracking wall panels | SITE (x3) POINT HOUSE ANALYS MOTION (x3) RELDISP (x3) NONLINEAR, then 7 iterations; MOTION | 30-120 s |
| `ex08_embedded_building.pre` | embedded shear-wall building with a two-level basement and a tower: FV against the subtraction (FI-FSIN) and modified subtraction (FI-EVBN) methods | SITE POINT HOUSE ANALYS MOTION STRESS (3 models) | 60-90 s |
| `ex09_steel_frame.pre` | three-storey braced steel frame (AISC W and HSS shapes, moment frames, composite slabs, equipment) on a 4 ft RC mat: SSI against a fixed base, ISRS, brace and column forces, storey drifts | SITE POINT HOUSE ANALYS MOTION STRESS RELDISP (x3), 2 models | ~30 s |

`data/` holds the input data:

* `rg160h_030g.rsi`: target spectrum for EQUAKE, the RG 1.60 horizontal design spectrum anchored
  to 0.30 g, 5 % damping, 27 frequency / spectral acceleration (g) pairs;
* `rg160h_030g.acc`: an acceleration history matched to it (EQUAKE, seed 11975, dt 0.005 s,
  20 s, PGA 0.324 g; first line dt, then one value in g per line: MOTION `<fopt>` = 0). It is the
  control motion of examples 1, 2, 5, 8 and 9, the rock-outcrop motion of example 6 (scaled to 0.30 g) and the
  input of example 7 (scaled to 0.6 g horizontally, 0.4 g vertically). Example 4 generates the same
  kind of record;
* `ricker_5hz.th`: a Ricker wavelet (5 Hz, peak 1 at t = 0.5 s), the load history of example 3.

## How to run

The model directory (`MDL,ex01,ex01`) and the data paths (`THFILE,../data/...`) are relative to
the **working directory**, which must be the `examples` directory. `MDL` also makes the model
directory the working directory, so `../data` is `examples/data`.

```bash
cd examples
../.venv/bin/sassi run ex01_surface_stick.pre
```

or, from the repository root, give the working directory with `--cwd`:

```bash
.venv/bin/sassi --cwd examples run ex01_surface_stick.pre
```

Do not run `sassi run examples/ex01_surface_stick.pre` from the repository root without
`--cwd examples`: the model directory would be created in the repository root (`./ex01`), the
control-motion file would not be found (CHECK Error 73), and MOTION would be blocked. The results of
example `exNN` are written to `examples/exNN/` (examples 2, 8 and 9 also write to `examples/ex02_fsin/`,
`examples/ex02_evbn/`, `examples/ex08_fsin/`, `examples/ex08_evbn/` and `examples/ex09_fixed/`). Example 6 also reads the soil curves `../sassi/data/dynp_library.pre`, so keep
the repository layout. Example 7 leaves about 105 MB in `examples/ex07/` (72 MB of binary inter-module
files such as the ANALYS restart files `COOTKnnn`, which can be deleted after the run); example 8 about 60 MB
in its three folders, example 9 about 58 MB in its two. Each example ends with `WRITE`, which saves the complete model as
commands (`examples/exNN/<model>.pre`). `INP` of that file rebuilds the same model.

The same commands can be typed in the console (`sassi`) or the browser GUI (`sassi-gui`); both
drive the same interpreter. In the GUI the module options are the tabs of Options > Analysis.

## How the examples map to the ACS SASSI workflow

ACS SASSI separates the user interface, which holds the model and writes the module input files,
from the analysis modules, which are batch programs that exchange binary files. SASSI-EDU keeps
this split:

```
.pre commands ──> model database ──> CHECK ──> AFWRITE ──> <model>.sit .poi .hou .anl .mot .str .rdi ...
                                                              │
   RUNSITE ──> FILE1 (free field), FILE2 (layer eigen-solutions)
   RUNPOINT ──> FILE3 (point-load solutions)                    ──> RUNANALYS ──> FILE8 (transfer functions)
   RUNHOUSE ──> FILE4 = <model>.N4, COOSK, COOSM (structure + excavated soil)        │
   RUNFORCE ──> FILE9 (loads, vibration analysis)                                   │
                                                                  RUNMOTION / RUNSTRESS / RUNRELDISP
```

`AOPT` selects the modules whose decks AFWRITE writes. A module with a CHECK error gets no deck
and cannot be run. Re-run AFWRITE after every change of the model or the options. `FCOPY` copies
files between runs, as the manual's file-management commands do (FILE1 to FILE1X, FILE9 to
FILE9001, ...).

The vendor's Demos 1-12 are not available to this project (requirements §1, traceability only);
the examples cover the same workflows:

* ex01: the basic surface-model workflow (Demo 1: `Demo1_X` input, control motion in X);
* ex02: embedded models (Demo 5 is the manual's reference embedded model), with the
  interaction-node sets of the manual (FV, FI-FSIN = subtraction method, FI-EVBN = modified
  subtraction method);
* ex03: foundation vibration (`ANALYS <type>` = 1, FORCE) and the global impedance option
  (`ANALYS <impe>`);
* ex04: EQUAKE with a target spectrum file `.rsi` (as in Demo 1) and the SOIL -> SITE chain
  (Non-Linear Soil option `SITEX,1`, FILE88);
* ex05: X, Y and Z inputs (manual §1.5.4, Demos 9 and 10) with simultaneous cases
  (`ANALYS <simul>` = 1);
* ex06: the secondary (near-field) soil nonlinearity of manual §1.5.4 and §4.1.2 item 18: the
  non-linear soil group of HOUSE (`.pin` input), New Structure restarts, STRESS effective strains and
  COMB_XYZ_STRAIN (`COMBXYZSTRAIN`), iterated by `NLSSIITER`;
* ex07: the Option NON workflow of manual Fig. 1.2 (Demo 9 is a reinforced-concrete shear-wall
  building at 0.60 g): wall panels, SHEAR/BBCGEN backbones, the NONLINEAR module and COMB_XYZ_THD
  (`COMBXYZTHD`), iterated by `NONLINITER`;
* ex08: a realistic embedded building (Demo 5 is the manual's reference embedded model): the separate
  meshes of the basement interior and the excavated soil (manual 1.5.1 rule 11, `EXCSTRCHK`), `INTGEN`
  for the FV, FI-FSIN and FI-EVBN sets on one model, and the validation of the reduced sets against FV
  that the manual's ASCE 4-16 warning asks for;
* ex09: a steel-framed structure (the BEAMS element with real sections, K nodes and end releases) on a
  surface mat, with the fixed-base reference run of lesson 1 and storey drifts by RELDISP.

## Example 1: lumped-mass stick on a rigid surface mat

**Model.** Units ft, kip, s (g = 32.2 ft/s²). A 64 ft x 64 ft SHELL mat of 8 ft elements, 5 ft thick
and ten times stiffer than concrete (concrete: E = 576,000 ksf = 4,000 ksi, 0.150 kcf; 9 x 9
nodes, `N`/`FILL`/`NGEN`, `E`/`EGEN`). On it stands a four-storey BEAMS stick (16 ft storeys,
2,200 kips per floor, given as weights with `MT` in a `FOREACH` loop). A rigid "spider" of BEAMS
elements connects the stick base to the eight surrounding mat nodes. The site is 16 ft of sand
(Vs 1,000 ft/s) over 39 ft of gravel (Vs 1,650 ft/s) over rock (Vs 3,300 ft/s). The input is a
vertically incident SV wave with the control motion in X at the free surface. `FIXROT` restrains
the shell drilling rotation only at the nodes connected to shells alone. The stick base and the
spider ends keep it, so the soil and not a fixity resists the torsion of the foundation. HOUSE
weighs 11,872 kips in all: 8,800 kips of floor masses and the 3,072 kip mat.

**Look at.**
* `ex01/FILE8` or `00085TR_X.TFU/.TFI`: the transfer function (ATF) of the roof. It is 1 at
  low frequency and has one SSI peak.
* `000nnTR_X01.RS` and `...02.RS`: the ISRS at the mat centre (41) and the floors (82-85), for
  2 % and 5 % damping.
* `000nnTR_X.THD`: the relative displacements of the floors with respect to the mat centre.
  The reference is `RELFILE,00041TR_X.TFI`.
* `BEAMS_002_0000k_FYI.THS` and `..._MZI.THS`: the storey shear and the bending moment at the
  bottom of each storey. `ex01_STRESS.out` lists their maxima.
* `00037TR_Z.TFU` and `00045TR_Z.TFU`: the vertical motion of the two mat edges on the x axis
  (rocking). The two are equal and opposite.

**Expected results** (integration tests in parentheses):
* The fixed-base frequency is 4.97 Hz (second mode 18.5 Hz; from the HOUSE matrices with the mat
  clamped, between 4.8 and 5.2 Hz). With SSI the roof ATF peaks at **3.47 Hz** with |ATF| = 13.1
  (between 3.2 Hz and 0.8 times the fixed-base frequency). At 0.1 Hz every ATF is 1.000 (within 1 %).
  The mat is stiff but not rigid: with the mat 1000 times as stiff as concrete the peak moves to
  3.83 Hz (|ATF| = 12.7), the value of a rigid-foundation estimate (lesson 4, Try this).
* The ISRS zero-period acceleration of the mat centre is **0.36 g** for an input PGA of 0.324 g
  (ratio 1.10; between 0.8 and 1.3). The ZPA grows up the stick: 0.47, 0.65, 0.93 and
  **1.28 g** at the roof. The roof ISRS peaks at 3.55 Hz with 15.4 g (2 %) and 8.0 g (5 %).
* Relative displacements of the floors: 0.0155, 0.0349, 0.0559 and **0.0766 ft** at the roof
  (0.19, 0.42, 0.67 and 0.92 in; between 0.015 and 0.3 ft). The mat centre relative to itself is
  exactly 0.
* Beam forces: the base shear is 6,760 kips and the base moment 3.29E5 kip ft. The storey shear
  decreases upwards (6,760, 6,140, 4,850 and 2,810 kips). The top-storey shear equals the roof
  weight times the roof acceleration (2,200 kips x 1.28 g), within 2 %, because the stick is
  massless.
* The rocking response of nodes 37 and 45 is antisymmetric to 1E-8.

**VP-E1** (`sassi/verify/problems/vp_examples.py`) runs this example through the interpreter up
to `RUNANALYS`. It then writes the same model by hand as module decks, with the primitives of
`sassi.verify.builders`, and runs the modules directly. The decks are identical, and the FILE8
of the two runs is identical (relative difference 0, tolerance 1E-12): AFWRITE writes exactly
the intended decks.

## Example 2: embedded box, FV versus FI-FSIN and FI-EVBN

**Model.** A 32 ft x 32 ft x 16 ft reinforced-concrete basement: SHELL walls (1.6 ft), a base slab
(3.2 ft) and a roof slab (1.6 ft). It is embedded in 16 ft of sand (Vs 800 ft/s) over gravel
(1,300 ft/s) over rock (2,600 ft/s). The excavated soil is 80 SOLID elements with `ETYPE` 2, in
five groups, one per 3.2 ft embedment layer.
The manual's rules for embedded models are followed:
* each embedment layer has its own L number in TOPL (`TOPL,1,2,3,4,5,...`). L numbers may be
  repeated for identical layers, but not for embedment layers;
* each excavated group uses the L number of its layer (`MACT,k` before group k).

CHECK warns about violations: EDU-28 for an L number repeated among the embedment layers, EDU-08
for excavated soil whose L differs from the TOPL layer at its depth. Interaction nodes are set
with `INT` in `FOREACH` loops over the first and last node of each level. Drilling rotations are
fixed face by face with `D` in loops. `FIXROT` cannot be used here, because the wall nodes are
also nodes of excavated solids.

Three analyses share the site and the embedment. `CPMODEL` copies the model, and `ACTM` and
`MDL` give the copy its own name and directory. `FCOPY` reuses FILE1 and FILE3 of the first run.
* model 0 `ex02`: FV, where every excavated node is an interaction node;
* model 2 `ex02fsin`: FI-FSIN, the lateral and bottom faces only;
* model 3 `ex02evbn`: FI-EVBN, the FI-FSIN set plus the top face.

**Look at.** The FILE8 of each directory, or `00138TR_X.TFU` (roof centre, X), `00013TR_X.TFU`
(base-slab centre, X) and `00136TR_Z.TFU` / `00140TR_Z.TFU` (roof edges, vertical = rocking). The
ISRS `000nnTR_X01.RS` (5 %), and the wall bending moments `SHELL_006_000nn_M*.THS` of the
STRESS run.

**Expected results.**
* All three sets give |ATF| = 1.000 at 0.1 Hz (within 1 %). The embedded box rocks under the
  vertically incident wave, and the rocking is antisymmetric (1E-8).
* FV: the roof X transfer function stays at or below 1 over the whole range. The stiff box
  follows the ground and kinematic interaction filters the motion. The roof-edge vertical TF
  grows to 0.24 at 20 Hz.
* **FI-FSIN** has a spurious resonance near **11 Hz** (the 22 frequencies include 11.0 Hz to
  resolve it): at 11.0 Hz roof |ATF_x| = 0.89 against 0.63 with FV (+40 % in amplitude; the complex
  difference |H - H_FV| is 64 % of the FV peak), base slab 0.55 against 0.34 and roof-edge vertical
  0.47 against 0.11. The roof slab rests on top-face nodes that are not
  interaction nodes. The slab and the subtracted soil under it then form an unphysical
  oscillator. This is the known weakness of the subtraction method.
* **FI-EVBN** follows FV within **1.5 %** of the FV peak at every frequency (complex difference; the
  largest amplitude difference is 2.4 %, at 20 Hz; tested: < 3 %). Adding the top face
  removes the spurious mode at a fraction of the FV interaction-node count.

## Example 3: forced vibration, compliance and impedance

**Model.** A rigid, massless 39 ft x 39 ft mat (7 x 7 interaction nodes at 6.5 ft) on a uniform
half-space (Vs 650 ft/s, nu 1/3, unit weight 0.125 kcf, G 1640 ksf, 2 % damping). 48 stiff
massless BEAMS ("rigid links", generated with `FOREACH` over a node list) tie the mat to its
centre node 25. For each of six unit load cases (`F` / `MM` at node 25), the example runs
AFWRITE and RUNFORCE, then `FCOPY,FILE9,FILE900k`. One ANALYS run (`<type>` = 1,
`<simul>` = 6, `<impe>` = 2) then solves all six cases and writes FILE8001 ... FILE8006 and the
global impedance K_G = T' X_ff T (FILE11, FOUNSTIF, FOUNDASH, FOUNDAMP, FOUNIMPD). MOTION
(`EDUOPT,TFFILE,FILE8003`, `MOTIONX` response = displacement) gives the vertical displacement
under a 225 kip Ricker pulse.

**Look at.** `FOUNSTIF` (Re K, 6 x 6 per frequency), `FOUNDASH` and `FOUNDAMP`. The 6 x 6
compliance of node 25 from FILE8001-8006, whose inverse is the impedance. The displacement
history `00025TR_Z.ACC` (the extension is `.ACC` whatever the response type).

**Expected results.**
* The two routes agree: the inverse of the 6 x 6 compliance equals K_G within **0.14-0.21 %**
  at every frequency. The links are stiff but not infinitely rigid. K is symmetric, with
  K_x = K_y and K_xx = K_yy (1E-8).
* Static stiffnesses at 0.2 Hz compared with Pais & Kausel (1988):

  | | computed | Pais-Kausel | difference |
  |---|---|---|---|
  | K_x = K_y (kip/ft) | 1.88E5 | 1.77E5 | +6.3 % |
  | K_z (kip/ft) | 2.39E5 | 2.25E5 | +5.8 % |
  | K_xx = K_yy (kip-ft/rad) | 8.75E7 | 7.30E7 | +19.8 % |
  | K_zz (kip-ft/rad) | 1.19E8 | 1.01E8 | +18.1 % |

  The 6.5 ft mesh (6 elements across the mat) over-predicts the stiffness, the moments most.
  VP-14 shows how the error decreases with mesh refinement.
* Ricker pulse, 225 kips: the peak vertical displacement is **0.00069 ft (0.0083 in)** at 0.525 s.
  The static value is 225/K_z = 0.00094 ft (0.0113 in). Because |K| grows with frequency, the
  dynamic peak is smaller. Radiation damping stops the motion almost at once: from t = 1 s (0.5 s
  after the pulse peak) it stays below 0.03 % of the peak.

## Example 4: the free-field chain EQUAKE -> SOIL -> SITE

**Model.** 30 ft of sand (Vs 650 ft/s) on 42 ft of clay (Vs 1,000 ft/s) on rock (Vs 4,000 ft/s),
in 3 ft sublayers in the sand and 3.5 ft sublayers in the clay. There is no structure. Gravity
32.2 ft/s² (units ft, kip, s; EQUAKE writes velocity in in/s and displacement in in).
1. EQUAKE matches a 20 s record to the RG 1.60 spectrum of `data/rg160h_030g.rsi`
   (`RSIN`/`RSOUT`/`ACCOUT`, seed 11975).
2. SOIL runs the equivalent-linear (SHAKE) analysis with the record as a rock-outcrop motion. It
   reads the strain-dependent curves of SHAKE91 with `INP` of `sassi/data/dynp_library.pre`
   (DYNP labels Sand and Clay) and assigns them per sublayer with `SPRO` in `FOREACH` loops. It
   writes the strain-compatible properties to FILE88.
3. `SITEX,1` (Non-Linear Soil) makes SITE read FILE88 instead of the low-strain L properties.

**Look at.** `ex04_EQUAKE.out` (spectral match and the SRP 3.7.1 acceptance table),
`rg160h_eq.rso`, `ex04_SOIL.out` (iterations, strain-compatible properties),
`FILE88` (text), `ACC001.TH` (surface motion), `RS001_01.RS` / `RS023_01.RS` (spectra at the
surface and of the rock outcrop), `SN005.TH` / `SN015.TH` (strains), and `SAF001W_023W.TFU`
(surface / base amplification).

**Expected results.**
* EQUAKE: PGA 0.324 g, PGV 22.5 in/s, PGD 16.0 in. The ratio to the target is 0.967 to 1.166,
  and every acceptance check passes: no point more than 10 % below or 30 % above, at most 9
  adjacent points below, 100 points per decade, 20 s. The strong-motion duration is 9.0 s.
* SOIL: 8 iterations. In the sand, the strain-compatible Vs falls from 650 ft/s to 620 ft/s at
  the surface and to **242 ft/s** at the bottom of the sand (sublayer 10, 27-30 ft), with damping
  2 % to 18 %. In the clay, Vs is 843-903 ft/s and damping 6-8 %. The surface / base (within)
  amplification peaks at **4.96 at 1.97 Hz**. With the low-strain properties (1 % damping) it
  would be 69 at 3.22 Hz. The surface PGA is 0.627 g for the 0.324 g rock-outcrop input.
* SITE with FILE88 reproduces the SOIL amplification: the difference is below **1 %** up to
  7.8 Hz and grows to 3.4 % at 13.7 Hz (3.6 % at 12.7 Hz). This is the thin-layer discretisation
  error, which grows as (k h)^2 (tested: < 1 % below 7.5 Hz, < 6 % up to the 13.7 Hz cut-off).
  The frequency set stops at 13.7 Hz because the softened sand (Vs 240-260 ft/s) passes
  Vs/(5h) = 16-17 Hz with 3 ft sublayers.
* For an SSI analysis, the control motion would now be the SOIL surface motion `ACC001.TH` at
  the top of layer 1 (within motion), consistent with FILE1.

## Example 5: X, Y and Z input with simultaneous cases

**Model.** A 40 ft x 40 ft basemat (5 ft thick) on the site of example 1. A two-storey concrete
stick (13 ft storeys) stands on it, and an eccentric 660 kip roof mass (node 30, 7.27 ft off the
axis) is tied to the stick top by a rigid arm. Masses are given as weights in kips with `MT`, and
the slab rotary inertias with `MR` (weight units: I x g, kip ft2); HOUSE divides by g = 32.2 ft/s2.
The three-component workflow of the manual (§1.5.4):
1. SITE three times: SV with control motion along x' (`FCOPY,FILE1,FILE1X`), SH along y'
   (`FILE1Y`) and P along z' (`FILE1Z`). The second and third runs reuse FILE2 (SITE Mode 1 off).
2. POINT and HOUSE once. ANALYS with `<simul>` = 1 solves the three cases with one impedance
   and one factorisation per frequency, and writes FILE8X, FILE8Y and FILE8Z.
3. MOTION once per direction (`EDUOPT,TFFILE,FILE8X` ...). The vertical motion is 2/3 of the
   horizontal one (MOTION `<mult>` = 0.6667). RELDISP runs for Y with the free-field reference
   (no RELFILE).

AFWRITE writes the control direction of the case (x', y', z') into the MOTION, STRESS and RELDISP
decks for FILE8X/Y/Z, whatever the last SITE run was. MOTION, STRESS and RELDISP also read it
back from the FILE8 itself.

MOTION names its files by node and DOF only. To keep the coupled response of the X run, the
example copies `00030TR_Y.*` and `00027R_ZZ.TFU` to `X2Y_*` names before the Y run writes its
own `00030TR_Y.*`.

**Look at.** `X2Y_00030TR_Y.TFU` and `X2Y_00030TR_Y01.RS` (Y response of the roof mass to the X
input), `X2Y_00027R_ZZ.TFU` (torsion of the stick top under X input), `00030TR_Y.TFU` (direct Y
response), `00027TR_Y.THD` (roof drift relative to the free field), and the `00030TR_*01.RS`
spectra.

**Expected results.**
* Each case moves with its own control motion at low frequency: |ATF| = 1.000 in X for FILE8X,
  in Y for FILE8Y and in Z for FILE8Z.
* The roof-mass ATF peaks near 5 Hz: 4.7 in X and 6.2 in Y. The mass is 6.5 ft off the stick
  axis in x and 3.25 ft in y, so the Y input twists the stick twice as much as the X input
  (stick-top torsion 0.225 against 0.112 rad per ft of control motion), and the twist adds to
  the Y motion of the mass. The vertical ATF peaks at 1.7 at 12 Hz (1.6 at 10 Hz).
* Coupling through the eccentric mass: the X input gives a Y response of the roof mass of up to
  **1.78** (0 at low frequency) and a torsion of the stick top. The direct Y response (6.2) is
  3.5 times larger. The 5 % ISRS of the roof mass peaks at 4.21 g (X), 4.07 g (Y) and 1.07 g
  for the Y response to the X input.
* RELDISP: the peak roof drift relative to the ground in Y is **0.0160 ft (0.19 in)**.

The three directional responses are combined afterwards (SRSS or 100-40-40). This
post-processing is outside the SASSI modules.

## Example 6: near-field soil nonlinearity (loose backfill behind a wall)

**Model.** A 2 ft reinforced-concrete wall (SHELL, 13 ft long, 6.5 ft high, its top at the ground
surface) carries 330 kips of a bridge deck at its top edge, like an abutment. Behind it (x = 0 to 13 ft)
the soil is a 13 ft x 13 ft x 6.5 ft block of loose backfill (Vs 500 ft/s at low strain) in a 32.5 ft
deposit of medium-dense sand (Vs 800 ft/s) on rock (Vs 2,600 ft/s). SSI has two kinds of soil nonlinearity (manual §1.5.4):
* **primary** (free field): SOIL iterates the layer properties of the site. Part 1 runs SOIL alone
  with the SHAKE91 Sand curves (`INP,../sassi/data/dynp_library.pre`, `SPRO` per sublayer) and the
  rock-outcrop motion scaled to 0.30 g (`SOILX`); `SACC` saves the surface motion `ACC001.TH`, the SSI
  control motion;
* **secondary** (near field): the backfill strains more than the free field. It is modelled with
  SOLID elements of the *structure* (group 3, ETYPE 1) at the nodes of the excavated soil, and their
  shear modulus and damping are iterated element by element.

The strain-compatible layers of FILE88 are copied into L 11-16 (each embedment layer has its own L
number, and the excavated groups use it with `MACT`). Every excavated node is an interaction node (FV).
The backfill material `M,2,1000,500,0.120,0.02,0.02,3` holds the **low-strain** properties (G_max = 932 ksf):
the curve softens it, so it must not be given strain-compatible values. `PIN,0.65` sets the effective
strain factor, `PINGRP,3,1,0.4,1.0,Sand` declares group 3 nonlinear (octahedral strain, curve Sand,
start at 0.4 x the free-field G: the native soil has G = 2,270 ksf, so the start is close to the backfill's
own G_max) and `HOUSEX,0,0,1,1` switches on non-linear SSI in HOUSE. The initiation run (three SITE runs
for FILE1X/Y/Z, POINT, HOUSE, ANALYS with restart files) is followed by STRESS with `<iter>` = 1 for each
direction (FILE74X/Y/Z) and `COMBXYZSTRAIN` (SRSS into FILE74). The command lists are kept in variables
(`VAR,NLX,...`): `FOREACH` runs them once, then `NLSSIITER,NLRUN+NLX+NLY+NLZ+NLCOMB,8` runs HOUSE
(new properties from FILE74) → ANALYS New Structure restart (`<mode>` 1) → STRESS x 3 → COMBXYZSTRAIN
until max |dG/G| < 2 % and max |d beta| < 0.5 %. MOTION finally gives the spectra of the wall top with
the converged backfill.

**Look at.** `ex06_HOUSE.out` (tables "Non-linear soil SSI", "Non-linear soil elements" and
"Non-linear soil convergence"), `NLSOIL_CONVERGENCE.TXT`, FILE78 (the properties used), FILE74X/Y/Z
and FILE74 (effective strains), `ex06_STRESS.out`, and the 5 % spectra `00022TR_X01.RS` (wall top) and
`00024TR_X01.RS` (back of the backfill).

**Expected results** (`tests/unit/test_nlsoil_example.py`):
* SOIL: rock outcrop 0.30 g → 0.58 g at the surface; the native Vs falls from 781 to 467 ft/s and the
  damping rises from 1.5 % to 10.7 % from the top to the bottom of the sand.
* The iterations converge in **4 passes**: max |dG/G| = 74, 26, 27, 9.8 and 0.7 % (iterations 0 to 4;
  iteration 2 overshoots and iteration 3 goes back). The converged backfill has G/G_max = 0.49-0.53 and
  7.0-7.7 % damping (effective octahedral strain 0.049-0.059 %); the elements next to the wall
  (x = 0 to 6.5 ft) strain and soften more than those at the back, clearly in the lower level.
* X input: the wall top (node 22, carrying the deck) reaches **0.77 g**, the back of the backfill
  (node 24) 0.58 g, the free-field value.

## Example 7: Option NON, a shear-wall building with cracking walls

**Model.** A 40 ft x 40 ft two-storey box (storeys of 13 ft) with 1 ft reinforced-concrete walls (meshed
10 ft x 6.5 ft) on a 5 ft stiff surface mat; 1.25 ft floor and roof slabs carry 875 kips of equipment each
(`MT`, `MTGEN`). Site: 32.5 ft of stiff soil (Vs 1,300 ft/s) in 25 sublayers on rock. Input: the RG 1.60 record
scaled to 0.6 g in X and Y and 0.4 g in Z. Each wall of each storey is a SHELL group with a material of
its own (`M,11` ... `M,18`): NONLINEAR changes E panel by panel (D-NON-12). The eight panels are low-rise
(h/l = 13/40) and shear governed.

The Option NON data: `PNLGEN` makes one panel per vertical shell group (Disp. Opt 1 = shear strain from
the four corner nodes, Force Opt 1 = Cheng-Mertz shear hysteresis); `SHEAR,0,4,60,0.005,0` lists
the ACI 318-08, Wood, Barda and Gulec-Whittaker capacities (f'c = 4 ksi, f_y = 60 ksi, 0.5 % web
steel; with the gravity 32.2 SHEAR and BBCGEN take ksi and kips); `BBCGEN,0,1,...,0.3` builds the 22-point backbone of every panel from the ACI capacity V_u
(cracking at 0.3 V_u, yield at 0.4 % strain, failure at 2 %); `EQL,0.8,1,0,0,1` sets the equivalent-linear
displacement factor 0.8 and adds the elastic 4 % damping to the hysteretic damping; `NONLINMOTDISP` adds
the panel corners to the MOTION and RELDISP requests. The workflow of manual Fig. 1.2:
1. elastic SSI: three SITE runs (FILE1X/Y/Z), POINT, HOUSE, ANALYS with restart files and X/Y/Z cases;
2. MOTION + RELDISP for each direction (`NONLINTHD,X` keeps the corner THD files as `X_*`), COMB_XYZ_THD
   (`COMBXYZTHD,COMB_XYZ_THD.inp`, the X + Y + Z sum per time step) and the elastic NONLINEAR run;
   `NONLINSAVE,elastic` keeps its results;
3. iterations: `NONLINITER,NLHOUSE+NLDX+NLDY+NLDZ+NLNON,10` copies `ex07_new.hou` to `ex07.hou` and runs
   HOUSE → ANALYS New Structure → MOTION + RELDISP x 3 → COMBXYZTHD → NONLINEAR until max |dE/E| < 2 % and
   max |d xi| < 0.5 % (D-NON-06). `NONLINBAT,1` writes the same loop as a generic script
   (`ex07_NONLINBAT.pre`).

**Look at.** `ex07_NONLINEAR.out` (panels, backbone curves, equivalent-linear properties, ductility,
convergence), `NONLINEAR_CONVERGENCE.TXT`, `Panel_EQL_Matl_Prop.txt` (the converged E and damping of
every panel), `Panel.fmu` (ductility mu and force-reduction factor F_mu), `Panel000k.thd` / `.ths`
(strain and force histories: plot one against the other for the hysteresis loops), `Panel000k.crv`
(E/E_el and damping versus strain) and the floor spectra `00063TR_X01.RS` of the converged model.

**Expected results** (`tests/unit/test_nonlinear_example.py`, marked `slow`):
* Convergence after **7 iterations**: max |dE/E| = 25.5, 20.9, 16.1, 11.1, 7.4, 4.0, 2.6 and 1.7 %
  (analyses 0 to 7).
* The four storey-1 walls crack: E/E_el = **0.37**, damping 16.5 % (4 % elastic + 12.5 % hysteretic),
  peak shear strain 3.6e-4, i.e. mu = 4.1 times the cracking strain, F_mu = 3.1. The storey-2 walls stay
  elastic (mu = 0.82, E/E_el = 1, 4 %).
* Converged model, X input: peak acceleration 0.65 g at the mat centre (node 13), 0.86 g at the floor
  centre (node 63) and 1.14 g at the roof centre (node 113).

Try `EQL,0.8,1,7,0,1` (damping cut-off 7 %, the ASCE 4 level for cracked concrete) or BBCGEN with
`<ShearModel>` 4 (Gulec-Whittaker) and see how the convergence changes; `docs/user/OPTION_NON.md`
sections 8 and 11 discuss both.

## Example 8: an embedded shear-wall building, FV against the subtraction methods

**Model.** A reinforced-concrete shear-wall building, 80 ft x 80 ft in plan, embedded 26 ft. Shear walls
on x, y = -40, -20, 0, 20, 40 ft make sixteen rooms of 20 ft x 20 ft below and above grade; the main block
rises two storeys of 16 ft and its four central rooms continue as a tower two storeys higher (roof at
64 ft). Units: ft, kip, s.

| Member (group, colour in the 3D view) | Thickness | Elements |
|---|---|---|
| basemat at z = -26 ft (5, purple) | 6.5 ft | 64 SHELL (10 ft x 10 ft) on the bottom face of the excavation |
| outer basement walls (12, brown) | 3.5 ft | 128 SHELL (10 ft x 6.5 ft) on its lateral faces |
| interior walls on x, y = -20, 0, 20 ft, basement, main block and tower (10, teal) | 2 ft | 304 SHELL |
| floor slabs: basement z = -13 ft, grade, z = 16 ft, tower z = 32 and 48 ft (9, pink) | 2, 2.5, 1.5 ft | 221 SHELL |
| roofs: main roof z = 32 ft around the tower, tower roof z = 64 ft (11, lavender) | 1.5 ft | 64 SHELL |
| outer walls above grade: main block (2.5 ft), tower (2 ft) (13, cream) | | 96 SHELL (10 ft x 16 ft) |

The group numbers are chosen for the colours of the 3D view (group g is drawn in colour g of the element
palette). Concrete E = 576,000 ksf (4,000 ksi), nu = 0.2, 0.150 kcf, 4 % damping (ASCE 4 response level 1
as the manual's Damping Cutoff note quotes it: uncracked stiffness, 4 %); slabs and roofs carry 0.025 kcf
more for distributed equipment, piping and live load. A 10 ft x 10 ft stair opening passes the basement
slab, the grade slab and the first floor (`EDEL`, then `ECOMPR`, because CHECK treats a missing element
number as an error). Ten equipment items of 55 to 330 kips (1,780 kips) sit at room centres (`MT`, as
weights). CALCM: 34,239 kips of elements, 36,019 kips with the equipment (1,118.6 kip-s2/ft); the
excavated soil weighs 19,968 kips.

The site: 26 ft of sand and gravel (Vs 1,000 ft/s, four 6.5 ft embedment layers with their own L
numbers), 42 ft of dense gravel (Vs 1,500 ft/s, 3.5 ft sublayers), 39 ft of very dense gravel
(Vs 2,100 ft/s, 6.5 ft sublayers) on rock (Vs 5,000 ft/s): 22 TOPL layers. The excavated soil is 256 SOLID
elements (10 ft x 10 ft x 6.5 ft) in four groups, one per embedment layer. The mesh passes
Vs/(5h) = 1000/(5 x 10) = 20 Hz in plan and 31 Hz vertically, the cut-off of the 41 SSI frequencies
(0.1 Hz, then every 0.5 Hz from 0.5 to 20 Hz). The input is the RG 1.60 record as a vertically incident
SV wave (control motion in X at the surface).

The building touches the excavated soil only on the foundation-soil interface (basemat and outer
walls). The interior walls and the two slabs inside the basement have nodes of their own at the
coordinates of the excavated-soil nodes (manual 1.5.1 rule 11; `EXCSTRCHK` finds no shared interior
node): they are generated as "excavation node + 1000", then `RMVUNUSED` and `NCOM` number the 781 nodes
without gaps (1-405 excavation, 406-567 main block, 568-617 tower, 618-781 basement interior). `FIXROT`
restrains the shell drilling rotations of the shell-only nodes, `D` in `FOREACH` loops those of the
face-interior nodes shared with the soil. `WINDOWSETTINGS,HIDEGROUP` hides the four soil groups in the 3D
views (display only), so the element plot and the gallery picture show the building closed, with the
interaction nodes on its basement walls; `WINDOWSETTINGS,SHOWGROUP,1` ... shows the soil again and
`WINDOWSETTINGS,VOLUME,,,3` (y >= 3 ft) cuts the building open north of y = 0 (lesson 11). `INTGEN`
builds the interaction sets on the same model and `CPMODEL` / `FCOPY` reuse the free field and the
point-load solutions, as in example 2:
* model 0 `ex08`: FV, 405 interaction nodes (1215 DOFs);
* model 2 `ex08fsin`: FI-FSIN (subtraction method, SM), 209 nodes;
* model 3 `ex08evbn`: FI-EVBN (modified subtraction method, MSM), 258 nodes.

FFV (`INTGEN,5`, 307 nodes) is not run here: on this model it costs as much as FV (lesson 11 has it
as a Try this).

**Look at.** The FILE8 of each folder, or the transfer functions `00605TR_X.TFU` (tower roof centre),
`00567TR_X.TFU` (main roof corner), `00757TR_X.TFU` (grade slab centre), `00041TR_X.TFU` (basemat centre),
`00045TR_Z.TFU` (basemat edge, rocking), `00739TR_Z.TFU` (grade slab under a 330-kip item) and
`00365TR_X.TFU` (the excavated soil at the centre of the excavation top face, inside the basement); the
5 % ISRS `000nnTR_X01.RS` and `000nnTR_Z01.RS`; the wall histories `SHELL_012_00004_FXY.THS` (south wall,
in-plane shear) and `SHELL_012_00069_MYY.THS` (west wall, bending); `ex08fsin.err` (EDU-12); the ANALYS
listings (interaction DOFs, memory estimate, time per frequency).

**Model checks.**
* CHECK: no error and no warning for FV; FI-FSIN and FI-EVBN only EDU-12 (Non-FV Method Selected:
  Validate Against FV). `EXCSTRCHK`, `FIXEDINT` and `HINGED` find nothing.
* At 0.1 Hz every transfer function is 1.000 (ANALYS low-frequency check: largest deviation 0.017 %).
* Fixed-base frequencies (an eigenvalue analysis of the HOUSE structural matrices `COOSK`/`COOSM`, done
  outside the program, the massless rotations condensed): with the basemat clamped the first lateral
  modes are at 12.6 Hz (X and Y), 64 % of the mass each, the next lateral ones at 22.5 Hz, the first
  strong vertical mode at 30.4 Hz; with the whole basement clamped, about 17 Hz (three close modes at
  16.8-17.4 Hz). These are in the range of squat shear-wall buildings.

**Expected results** (`tests/integration/test_examples.py`, marked `slow`):
* Run time about 60-90 s on a laptop: ANALYS 26 s (FV), 21 s (FI-FSIN) and 25 s (FI-EVBN) for 41
  frequencies, measured with other applications loading the machine. The dense impedance is 23.6 MB per matrix for FV and 6.3 MB for FI-FSIN, but on a
  model this small the 3930 equations of the structure and the excavation cost as much as the
  impedance: FV takes only 1.25 times as long as FI-FSIN. The manual's tens of times apply to excavations with thousands of nodes.
  About 60 MB in the three folders.
* FV: the tower roof X transfer function peaks at **2.78 at 6.5 Hz** (the SSI frequency, about half the
  fixed-base one), the main roof corner at 1.73; the basemat never moves more than the free surface.
  5 % ISRS zero-period accelerations: 0.283 g at the basemat (kinematic interaction: below the 0.324 g
  control motion), 0.321 g at grade, 0.361 g at the main roof and 0.451 g at the tower roof; the tower
  roof ISRS peaks at 2.23 g at 6.8 Hz. South basement wall FXY 7.20 ksf (50 psi), west wall MYY
  10.4 kip-ft/ft (bottom panels).
* **FI-FSIN** follows FV within 2.2 % (complex difference, relative to the FV peak) below 10 Hz and
  within 5.6 % from 10 to 14 Hz, then has a **spurious resonance near 16 Hz**: at 16 Hz basemat 0.923
  against 0.288, grade slab 0.982 against 0.332, tower roof 1.809 against 0.565, and the soil inside the
  basement (node 365) 18.4 against 1.8; a dip precedes it (the tower roof 0.118 against 0.603 at
  15.5 Hz, the main roof corner 0.039 against 0.116 at 14.5 Hz). `CRITFREQ` flags the interpolated SM
  peaks near 15.85 Hz as not supported by computed values. With 0.25 Hz steps (a scratch run) the computed peak is at
  15.75 Hz (basemat 1.11, tower roof 1.35, enclosed soil 21.6, against 0.30, 0.58 and 1.5 with FV); at
  15.84 Hz, the frequency CRITFREQ proposes, the basemat moves 2.10, the tower roof 3.28 and the
  enclosed soil 41.7 (FV: 0.29, 0.58 and 1.6).
* The 5 % ISRS of FI-FSIN differ from FV by at most 2.3 % below 10 Hz, are **up to 16 % too high at
  15.1-15.9 Hz** (basemat X +11 %, grade slab X +12 %, main roof X +13 %, basemat edge Z +16 %) and up to
  10 % too low at 13.8-14.8 Hz, and differ by up to 16 % above 17 Hz; the peak wall forces change by
  about 1.5 % (FXY -1.2 %, MYY +1.5 %): they come from the 6.5 Hz SSI response.
* **FI-EVBN** follows FV within **5.1 %** at every computed frequency and output (the largest difference
  is the grade slab Z under the 330-kip item at 20 Hz, where the FI-EVBN set approaches its own
  enclosed-soil frequency, 22.7 Hz; X within 1.1 %), the ISRS within 1.2 %, the wall forces within
  1.3 %. FFV (lesson 11, Try this) gives 4.3 % at the same place, X within 1.4 %, ISRS within 1.1 %.

**Why the subtraction method resonates.** In FI-FSIN the 196 excavated nodes inside the basement and
on the excavation top face carry neither impedance nor free-field load: their rows keep only
-(K_e - w2 M_e). No structure rests on them (the interior structure has its own nodes), so the anomaly
comes from the soil alone: the enclosed soil, held at the walls and the basemat and free at the top,
has natural frequencies of its own, and near the first one the subtracted soil resonates through the
interface. An eigenvalue analysis of the HOUSE matrices `K_e`, `M_e` with the interaction DOFs fixed
(done outside the program) gives 15.6 Hz for the FI-FSIN set, 22.7 Hz for FI-EVBN (above the 20 Hz
cut-off: no anomaly) and 38.5 Hz for FFV; they depend on the excavation only, not on the building. How
much of the anomaly reaches the ISRS depends on the building and on the input, not only on the soil:
here the floor spectra take 9-16 % of error at 15.1-15.9 Hz, while the transfer functions are off by a
factor of about 3. The frequency scales with the soil velocity. The same model with other embedment
soils (L 1-4 changed, everything else equal; scratch runs):
* Vs 2,000 ft/s (stiff): FI-FSIN follows FV within 1.5 % up to 20 Hz, its ISRS within 0.9 %; FI-EVBN
  within 1.1 %. This is the manual's statement that FI-FSIN agrees with FV on stiff soil and rock sites.
* Vs 650 ft/s (soft): FI-FSIN resonates at 10.5 Hz (basemat 0.93 against 0.59; enclosed-soil frequency
  10.1 Hz), ISRS +39 % / -20 %; FI-EVBN departs from FV by up to 46 % near 19.5 Hz (the grade slab Z
  under the 330-kip item; X up to 19 %), its own enclosed-soil frequency having fallen to 14.8 Hz, ISRS
  +1.4 % / -5.1 %. (The 10 ft mesh passes only 13 Hz for this soil: the comparison of the methods
  holds, the absolute results above 13 Hz do not.)

Lesson 11 of the guided course builds this example step by step, shows the soil and the interaction
nodes and a cutaway of the building, animates the FV and FI-FSIN motion at 16 Hz (`HARMFRAME`) and
discusses the validation a design engineer must document.

## Example 9: a braced steel frame on a reinforced-concrete mat

**Model.** A three-storey steel building of the kind that houses auxiliary systems: 60 ft x 40 ft in
plan (column lines 1-4 at x = -30, -10, 10, 30 ft and A-C at y = -20, 0, 20 ft, 20 ft bays), storeys of
16, 15 and 15 ft (roof at 46 ft), on a 4 ft mat of 70 ft x 50 ft (5 ft beyond the outer column lines)
on the ground surface. The lateral system:

* X: concentric X-braces in the end bays of lines A and C, a welded moment frame on line B;
* Y: X-braces in both bays of lines 1 and 4, welded moment frames on lines 2 and 3;
* the middle bays are free of braces; every other girder and the floor beams have shear-tab connections.

| Member (group, colour in the 3D view) | Shape (AISC Shapes Database v15.0) | A | Ix | Iy | J | Elements |
|---|---|---|---|---|---|---|
| columns, storeys 1-2 (3, blue) | W14x132 | 38.8 in2 | 1530 in4 | 548 in4 | 12.3 in4 | 24 BEAMS |
| columns, storey 3 (3) | W14x90 | 26.5 | 999 | 362 | 4.06 | 12 |
| moment-frame girders, lines B, 2, 3 (1, red) | W24x76 | 22.4 | 2100 | 82.5 | 2.68 | 42 |
| girders of lines A and C (4, orange) | W24x76 | | | | | 36 |
| edge girders of lines 1 and 4, floor beams on x = -20, 0, 20 (4) | W21x44 | 13.0 | 843 | 20.7 | 0.770 | 60 |
| braces, storeys 1-2 (2, green) | HSS8x8x1/2 | 13.5 | 125 (0 in the model) | 125 (0) | 204 | 32 |
| braces, storey 3 (2) | HSS6x6x1/2 | 9.74 | 48.3 (0) | 48.3 (0) | 81.1 | 16 |
| composite slabs, floors 1-2 (9) and roof (11) | 6 in (0.5 ft) SHELL | | | | | 72 |
| mat (13, cream) | 4 ft SHELL, 5 ft mesh | | | | | 140 |

The R table gets the AISC values in ft units (in2 / 144 = ft2, in4 / 20,736 = ft4) with the shear
areas As2 = d tw (web) and As3 = 5/6 x 2 bf tf (flanges). **Orientation.** Local axis 2 of a BEAMS
element points to its K node and I3 goes with bending in the 1-2 plane (spec 08 4.5), so the K node is
put on the web side of every W shape and I3 = Ix: the columns of line B have their webs along X (strong
axis in the X moment frame), the others along Y (strong axis in the moment frames of lines 2 and 3),
and the girders and beams have K 15 ft above (web vertical; the roof beams use a level of K-only nodes,
which are not drawn and which AFWRITE fixes). Check: the element matrices of a 16 ft W14x132
cantilever with K in +X give a tip deflection of 1.100 in in X and 2.983 in in Y under 20 kips, equal
to P L3/(3 E Ix) + P L/(G d tw) and P L3/(3 E Iy) + P L/(G As3); a W24x76 girder along X with K above
is stiff vertically (0.210 in on a 10 ft cantilever) and flexible sideways (4.84 in). **Connections.**
Shear tabs release M3 at the supported end of each girder or beam (two 10 ft elements per member, `KI`
on the odd and `KJ` on the even elements). The four corner columns are pinned at the base (`KI`
releases M2 and M3); the other eight are moment-frame columns fixed into the mat (the mat's rotational
stiffness at a node, 1.3E7 kip ft/rad with its translations held, is about 170 times 4EI/L of a
W14x132, 7.7E4 kip ft/rad, so no rigid spider is needed). The braces are pin-ended axial members:
`KI`/`KJ` cannot release the same moment at both ends of an element (CHECK Error 10), so their R rows
have I2 = I3 = 0 (a truss member, as LINK180 in ANSYS); the two diagonals of an X-brace cross without a
node. **Slabs.** The 6 in SHELL slabs (the average concrete thickness of a 3 in deck with 4.5 in of
concrete above the ribs) share the beam nodes, so they are the diaphragms but not composite with the
beams: the floor vertical frequencies are lower bounds. `FIXROT` fixes the drilling rotation of the
shell-only mat nodes and `D` that of the twelve column bases (without it, local torsion modes of these
nodes appear at 6.8-13.6 Hz); the floor nodes get their stiffness about Z from the beams (only M3 is
released). `WINDOWSETTINGS,HIDEGROUP` hides the slabs in the 3D views (display only), so the picture
shows the steel frame on the mat.

Materials: mat concrete E = 4,000 ksi = 576,000 ksf (f'c 5 ksi), slabs 3,600 ksi = 518,400 ksf (f'c
4 ksi), steel E = 29,000 ksi = 4,176,000 ksf, nu = 0.3, 0.490 kcf. Damping is hysteretic and per
element, from ASCE 4-16 Table 3-2 at response level 1 (the level of the concrete values the manual
quotes in its Damping Cutoff note, as in ex08): concrete 4 %, welded steel 2 % (columns and
moment-frame girders), bearing-bolted steel 4 % (braces, shear-tab girders and beams). The slab weights
carry the superimposed loads: floors 135 psf (slab and deck 80, services 30, 25 % of a 100 psf live
load), roof 115 psf (slab and deck 80, roofing and services 20, 75 % of a 20 psf snow load). Equipment
(`MT`, weights in kips, placed symmetrically): two 65 kip heat exchangers on the line B girders of
floor 1, four 22 kip cabinets on the girders of lines 2 and 3 of floor 2, two 25 kip air-handling
units on the roof. CALCM and HOUSE: **3,494 kips** (mass 108.5 kip s2/ft; mat 2,100 kips, slabs 924,
steel 202, equipment 268); 1,373 kips act at the three floor levels (527, 477 and 369 kips).

The site: 19.5 ft of medium-dense sand (Vs 800 ft/s), 39 ft of dense sand and gravel (Vs 1,300 ft/s),
39 ft of very dense gravel (Vs 2,000 ft/s) on rock (Vs 5,000 ft/s), 24 TOPL layers. The 165 mat nodes
are the interaction nodes (495 DOFs); the 5 ft mesh passes Vs/(5h) = 800/(5 x 5) = 32 Hz, above the
30 Hz cut-off. 97 SSI frequencies: 0.1 Hz, every 0.5 Hz to 2.5 Hz, every 0.098 Hz from 3.0 to 5.0 Hz
(the sharp peaks of a frame with 2-4 % damping), every 0.25 Hz to 14 Hz and every 0.5 Hz to 30 Hz.
Input: the RG 1.60 record as a vertically incident SV wave, control motion in X at the free surface.
Model 2 (`ex09fb`, results in `ex09_fixed/`) is the same frame on a practically rigid site (all layers
Vs = 33,000 ft/s), the fixed-base reference of lesson 1.

**Look at.**
* `ex09/00253TR_X.TFI` against `ex09_fixed/00253TR_X.TFI`: the roof transfer function with and without
  SSI; `00183TR_X` and `00218TR_X` the floors, `00083TR_X` the mat centre.
* `00nnnTR_X01.RS` (2 %) and `...02.RS` (5 %): the ISRS at the mat centre (83), the floor centres (183,
  218), the roof (253), the heat exchanger (185) and the air-handling unit (255).
* `00076TR_Z` and `00090TR_Z` (mat edges, rocking), `00017TR_Z` and `00021TR_Z` (the bases of columns A1
  and A2, a braced bay), `00185TR_Z` and `00255TR_Z` (vertical response at the equipment).
* `BEAMS_002_0000k_FXI.THS`: the axial force of the storey-1 braces of line A (1-4) and of one brace in
  storeys 2 (17) and 3 (33); `BEAMS_003_00001_FXI.THS` / `00002` (axial force of the braced-frame
  columns A1 and A2), `BEAMS_003_00005_FYI` / `_MZI.THS` (shear and strong-axis moment at the base of
  the moment-frame column B1); `ex09_STRESS.out` lists the maxima.
* `00183TR_X.THD`, `00218TR_X.THD`, `00253TR_X.THD`: the inter-storey drifts (each floor relative to the
  floor below: three RELDISP runs with `RELFILE` = the `.TFI` of the floor below).

**Model checks.**
* CHECK: no error and no warning for either model. `INTCOUNT`: 165 interaction nodes. `FIXEDINT` and
  `KINT` find nothing. `HINGED` lists four nodes where a single column meets coplanar shells (the bases
  and roof tops of B2 and B3): the column's torsion cannot pass into a shell. The bases have ROTZ fixed,
  at the roof the moment-frame girders carry it; nothing is hinged. CPMODEL warns that the copy has
  the name of the original until `MDL` renames it (as in lesson 1).
* At 0.1 Hz every transfer function is 1.000 (largest deviation 0.08 %), the vertical ones 0.
* Fixed-base frequencies (an eigenvalue analysis of the HOUSE matrices `COOSK`/`COOSM` with the mat
  nodes clamped, done outside the program): X 4.38 Hz (84 % of the mass above the mat), Y 4.84 Hz (88 %),
  torsion 8.32 Hz, floor vertical modes 7.1-9.6 Hz, second X mode 11.7 Hz (12 %); 99.9 % of the X mass
  is in modes below 20 Hz. The rigid-site run peaks at 4.37 Hz, the computed frequency next to 4.38 Hz.
* The model is doubly symmetric: braces 1 and 4 (and 2 and 3) have the same forces, the four corner
  columns the same axial force, and the mat edges move vertically in antiphase (sum 3E-14 of either).
* **Brace forces against a hand estimate.** The floor-centre accelerations times the level masses give
  a base shear of 1,344 kips with SSI (0.98 times the weight above the mat) and 1,276 kips on the fixed
  base. If the 8 storey-1 X-braces carried all of it, each would carry V/(8 cos theta) with cos theta =
  20/25.6 = 0.781: 215 and 204 kips. The computed maxima are 208 and 198 kips (97 %): at the instant of
  the peak shear the braces carry 95 % of it, the moment frames the rest. 208 kips is about 0.4 of the
  design compression strength of an HSS8x8x1/2 with a 12.8 ft buckling length (AISC 360, Fy 50 ksi:
  504 kips), consistent with response level 1.

**Expected results** (`tests/integration/test_examples.py`, marked `slow`):
* Run time about 30 s for the two models (ANALYS 6 s and MOTION 7 s each); about 58 MB in `ex09/` and
  `ex09_fixed/` (FILE2 alone is 16 MB per model).
* **Frequency shift.** Roof X transfer function: **20.2 at 4.37 Hz** on the fixed base, **18.3 at
  3.98 Hz** with SSI (-9 % in frequency and in amplitude). The mat moves 1.48 times the free field at
  3.88 Hz, and it rocks: its edges move 0.68 vertically per unit control motion at 3.98 Hz. The
  steel's 2-4 % damping gains little radiation damping, because rocking radiates little at 4 Hz.
* **ISRS** (X; SSI against fixed base). Zero-period accelerations: mat 0.363 g against 0.325 g (the
  free field, PGA 0.324 g), floor 1 0.654 / 0.645 g, floor 2 1.008 / 1.002 g, roof **1.623 / 1.577 g**.
  5 % peaks: floor 1 3.77 / 3.63 g, floor 2 7.10 / 6.87 g, roof **10.27 g at 4.07 Hz / 10.01 g at
  4.27 Hz** (+3 to +4 %). 2 % peaks: 7.25 / 5.94 g, 13.96 / 10.88 g, roof **20.27 / 15.52 g** (+22 to
  +31 %): the 2 % spectrum of this record has a peak at 4.07 Hz (1.54 g against 1.35 g at 4.38 Hz)
  onto which the SSI frequency moves. The equipment nodes follow their floors (heat exchanger 3.79 g,
  air-handling unit 10.29 g at 5 %).
* **Vertical response under X input.** Fixed base: 0.29 at the heat exchanger (12.3 Hz) and 0.62 at
  the air-handling unit (8.7 Hz), from the overturning of the frame. With SSI the mat rocking excites
  the floor-1 vertical mode: 0.34 at 7.25 Hz at the heat exchanger (5 % ISRS 0.18 g against 0.12 g),
  while the roof unit drops to 0.45 (0.23 g against 0.27 g).
* **Member forces** (maxima, SSI / fixed base): storey-1 braces **208 / 198 kips** (+5 %), storey 2
  160 / 156 kips, storey 3 88.5 / 86.5 kips; braced-frame columns A2 404 / 400 kips and A1 328 /
  327 kips axial; moment-frame column B1 151 / 130 kips axial, 10.9 / 9.6 kips shear and 109 /
  98 kip·ft at the base; B2 132 / 111 kip·ft; the line B girder at B1 9.7 / 8.1 kips and 104 /
  93 kip·ft. The braces carry axial force only; the pinned bases carry no moment.
* **Drifts** (RELDISP, each floor relative to the floor below): **0.304, 0.309 and 0.297 in** with SSI
  (drift ratios 0.16-0.17 %) against 0.250, 0.266 and 0.255 in on the fixed base (+22, +16, +17 %),
  while the braces carry only 2-5 % more. The difference is foundation rocking: the bases of the braced
  bay A1-A2 move 0.036 and 0.017 in vertically, more than the rigid rocking of the mat would give there
  (0.030 and 0.010 in: the 4 ft mat bends under the braced-frame columns); the rotation of the bay
  times 16 ft is 0.042 in, and the storey-1 drift without it is 0.262 in, 5 % above the fixed base like
  the brace forces. Pipes and ducts between floors see the total drift, the braces only the racking.

**What it shows.** For a light steel frame on a heavy mat, SSI lowers the frequency by about 10 % and
hardly reduces the peaks, because the structure has little damping and rocking radiates little energy.
Whether the ISRS and member forces go up or down then depends on the input spectrum at the two
frequencies, which is why design ISRS are enveloped over soil cases and broadened. The drifts include
rocking of the foundation that the fixed-base model does not have.

## Tests

* `tests/integration/test_examples.py` runs examples 1-5, 8 and 9 in a temporary copy, through the
  interpreter as `sassi --cwd examples run` does. It asserts that every module finishes with
  status OK, checks the key results above (with physically motivated bounds), checks that
  `WRITE` -> `INP` gives back the same model for every model of every example, and checks the
  documented `--cwd` command line. The example 4, 8 and 9 tests are marked `slow` (EQUAKE; three ANALYS
  runs; two models of 97 frequencies).
* `tests/unit/test_nlsoil_example.py` runs example 6: a clean run, the convergence history, the
  wall-top spectrum and the `WRITE` -> `INP` round trip.
* `tests/unit/test_nonlinear_example.py` runs example 7 (marked `slow`): a clean run, the panels and
  backbone curves, the convergence, the cracked first storey, the next HOUSE deck and the round trip.
* `tests/verification/test_vp_examples.py` runs **VP-E1**.

```bash
.venv/bin/python -m pytest -q tests/integration tests/verification/test_vp_examples.py
.venv/bin/python -m pytest -q tests/unit/test_nlsoil_example.py tests/unit/test_nonlinear_example.py
```
