# SASSI-EDU tutorial examples

Eight small models that run the complete ACS SASSI workflow the way a user runs it: a `.pre` file
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
| `ex06_nonlinear_soil.pre` | loose backfill behind an embedded wall: near-field soil nonlinearity by equivalent-linear SSI iterations | SOIL SITE (x3) POINT HOUSE ANALYS STRESS (x3), then 6 iterations of HOUSE ANALYS (New Structure) STRESS (x3); MOTION | ~5 s |
| `ex07_option_non.pre` | Option NON: two-storey shear-wall building with cracking wall panels | SITE (x3) POINT HOUSE ANALYS MOTION (x3) RELDISP (x3) NONLINEAR, then 7 iterations; MOTION | 30-120 s |
| `ex08_embedded_building.pre` | embedded shear-wall building with a two-level basement and a tower: FV against the subtraction (FI-FSIN) and modified subtraction (FI-EVBN) methods | SITE POINT HOUSE ANALYS MOTION STRESS (3 models) | 60-90 s |

`data/` holds the input data:

* `rg160h_030g.rsi`: target spectrum for EQUAKE, the RG 1.60 horizontal design spectrum anchored
  to 0.30 g, 5 % damping, 27 frequency / spectral acceleration (g) pairs;
* `rg160h_030g.acc`: an acceleration history matched to it (EQUAKE, seed 11975, dt 0.005 s,
  20 s, PGA 0.324 g; first line dt, then one value in g per line: MOTION `<fopt>` = 0). It is the
  control motion of examples 1, 2, 5 and 8, the rock-outcrop motion of example 6 (scaled to 0.30 g) and the
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
example `exNN` are written to `examples/exNN/` (examples 2 and 8 also write to `examples/ex02_fsin/`,
`examples/ex02_evbn/`, `examples/ex08_fsin/` and `examples/ex08_evbn/`). Example 6 also reads the soil curves `../sassi/data/dynp_library.pre`, so keep
the repository layout. Example 7 leaves about 95 MB in `examples/ex07/` (60 MB of binary inter-module
files such as the ANALYS restart files `COOTKnnn`, which can be deleted after the run); example 8 about 57 MB
in its three folders. Each example ends with `WRITE`, which saves the complete model as
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
  that the manual's ASCE 4-16 warning asks for.

## Example 1: lumped-mass stick on a rigid surface mat

**Model.** A 20 m x 20 m SHELL mat, 1.5 m thick and ten times stiffer than concrete (9 x 9
nodes, `N`/`FILL`/`NGEN`, `E`/`EGEN`). On it stands a four-storey BEAMS stick (5 m storeys,
1000 t per floor, given as weights with `MT` in a `FOREACH` loop). A rigid "spider" of BEAMS
elements connects the stick base to the eight surrounding mat nodes. The site is 5 m of sand
(Vs 300 m/s) over 12 m of gravel (Vs 500 m/s) over rock (Vs 1000 m/s). The input is a vertically
incident SV wave with the control motion in X at the free surface. `FIXROT` restrains the shell
drilling rotation only at the nodes connected to shells alone. The stick base and the spider ends
keep it, so the soil and not a fixity resists the torsion of the foundation.

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
* The fixed-base frequency is 5.05 Hz (second mode 19 Hz). With SSI the roof ATF peaks at
  **3.49 Hz** with |ATF| = 13.1 (between 3 and 4 Hz). At 0.1 Hz every ATF is 1.000 (within 1 %).
  The mat is stiff but not rigid: with the mat 1000 times as stiff as concrete the peak moves to
  3.86 Hz (|ATF| = 12.6), the value of a rigid-foundation estimate (lesson 4, Try this).
* The ISRS zero-period acceleration of the mat centre is **0.36 g** for an input PGA of 0.324 g
  (ratio 1.11; between 0.8 and 1.3). The ZPA grows up the stick: 0.48, 0.66, 0.94 and
  **1.29 g** at the roof. The roof ISRS peaks at 3.5 Hz with 15.7 g (2 %) and 8.2 g (5 %).
* Relative displacements of the floors: 4.7, 10.6, 16.9 and **23.2 mm** at the roof. The mat
  centre relative to itself is exactly 0.
* Beam forces: the base shear is 30.3 MN and the base moment 4.6E5 kN m. The storey shear
  decreases upwards. The top-storey shear (12.6 MN) equals the roof mass times the roof
  acceleration (1000 t x 1.29 g), within 2 %, because the stick is massless.
* The rocking response of nodes 37 and 45 is antisymmetric to 1E-8.

**VP-E1** (`sassi/verify/problems/vp_examples.py`) runs this example through the interpreter up
to `RUNANALYS`. It then writes the same model by hand as module decks, with the primitives of
`sassi.verify.builders`, and runs the modules directly. The decks are identical, and the FILE8
of the two runs is identical (relative difference 0, tolerance 1E-12): AFWRITE writes exactly
the intended decks.

## Example 2: embedded box, FV versus FI-FSIN and FI-EVBN

**Model.** A 10 m x 10 m x 5 m reinforced-concrete basement: SHELL walls (0.5 m), a base slab
(1.0 m) and a roof slab (0.5 m). It is embedded in 5 m of sand over gravel over rock. The
excavated soil is 80 SOLID elements with `ETYPE` 2, in five groups, one per 1 m embedment layer.
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
* **FI-FSIN** has a spurious resonance at **6.0 Hz**: roof |ATF_x| = 1.34 against 0.85 with FV (+58 % in
  amplitude; the complex difference |H - H_FV| is 63 % of the FV peak) and roof-edge vertical 0.58
  against 0.03. The roof slab rests on top-face nodes that are not
  interaction nodes. The slab and the subtracted soil under it then form an unphysical
  oscillator. This is the known weakness of the subtraction method.
* **FI-EVBN** follows FV within **1.5 %** of the FV peak at every frequency (complex difference; the
  largest amplitude difference is 2.5 %, at 20 Hz; tested: < 3 %). Adding the top face
  removes the spurious mode at a fraction of the FV interaction-node count.

## Example 3: forced vibration, compliance and impedance

**Model.** A rigid, massless 12 m x 12 m mat (7 x 7 interaction nodes at 2 m) on a uniform
half-space (Vs 200 m/s, nu 1/3, rho 2 t/m3, 2 % damping). 48 stiff massless BEAMS ("rigid
links", generated with `FOREACH` over a node list) tie the mat to its centre node 25. For each of
six unit load cases (`F` / `MM` at node 25), the example runs AFWRITE and RUNFORCE, then
`FCOPY,FILE9,FILE900k`. One ANALYS run (`<type>` = 1, `<simul>` = 6, `<impe>` = 2) then solves
all six cases and writes FILE8001 ... FILE8006 and the global impedance K_G = T' X_ff T (FILE11,
FOUNSTIF, FOUNDASH, FOUNDAMP, FOUNIMPD). MOTION (`EDUOPT,TFFILE,FILE8003`, `MOTIONX` response =
displacement) gives the vertical displacement under a 1000 kN Ricker pulse.

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
  | K_x = K_y (kN/m) | 2.82E6 | 2.65E6 | +6.3 % |
  | K_z (kN/m) | 3.58E6 | 3.38E6 | +5.8 % |
  | K_xx = K_yy (kN m/rad) | 1.24E8 | 1.04E8 | +19.6 % |
  | K_zz (kN m/rad) | 1.69E8 | 1.44E8 | +17.9 % |

  The 2 m mesh (6 elements across the mat) over-predicts the stiffness, the moments most.
  VP-14 shows how the error decreases with mesh refinement.
* Ricker pulse, 1000 kN: the peak vertical displacement is **0.205 mm** at 0.525 s. The static
  value is 1000/K_z = 0.279 mm. Because |K| grows with frequency, the dynamic peak is
  smaller. Radiation damping stops the motion almost at once: from t = 1 s (0.5 s after the
  pulse peak) it stays below 0.1 % of the peak.

## Example 4: the free-field chain EQUAKE -> SOIL -> SITE

**Model.** 10 m of sand (Vs 200 m/s) on 12 m of clay (Vs 300 m/s) on rock (Vs 1200 m/s), in
1 m sublayers. There is no structure.
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
* EQUAKE: PGA 0.324 g. The ratio to the target is 0.967 to 1.166, and every acceptance check
  passes: no point more than 10 % below or 30 % above, at most 9 adjacent points below, 100
  points per decade, 20 s. The strong-motion duration is 9.0 s.
* SOIL: 8 iterations. In the sand, the strain-compatible Vs falls from 200 m/s to 191 m/s at the
  surface and to **74 m/s** at 10 m depth, with damping 2 % to 18 %. In the clay, Vs is 255-268
  m/s and damping 6-8 %. The surface / base (within) amplification peaks at **4.65 at 1.83 Hz**.
  With the low-strain properties (1 % damping) it would be 68 at 3.10 Hz.
* SITE with FILE88 reproduces the SOIL amplification: the difference is below **0.5 %** up to
  7 Hz and grows to 4.7 % at 13.7 Hz. This is the thin-layer discretisation error, which grows
  as (k h)^2 (tested: < 1 % below 7.5 Hz, < 6 % up to the 13.7 Hz cut-off). The frequency set
  stops at 13.7 Hz because the softened sand (Vs ~ 75 m/s) passes Vs/(5h) = 15 Hz with 1 m
  sublayers.
* For an SSI analysis, the control motion would now be the SOIL surface motion `ACC001.TH` at
  the top of layer 1 (within motion), consistent with FILE1.

## Example 5: X, Y and Z input with simultaneous cases

**Model.** A 12 m x 12 m basemat on the site of example 1. A two-storey stick stands on it, and
an eccentric 300 t roof mass (node 30, 2.24 m off the axis) is tied to the stick top by a rigid
arm. Masses are given with `MT`, and the slab rotary inertias with `MR` (weight units: I x g).
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
* The roof-mass ATF peaks near 5 Hz: 5.0 in X and 6.5 in Y. The mass is 2 m off the stick axis
  in x and 1 m in y, so the Y input twists the stick twice as much as the X input (stick-top
  torsion 0.74 against 0.37 rad/m), and the twist adds to the Y motion of the mass. The vertical
  ATF peaks at 1.7 near 10 Hz.
* Coupling through the eccentric mass: the X input gives a Y response of the roof mass of up to
  **1.85** (0 at low frequency) and a torsion of the stick top. The direct Y response (6.5) is
  3.5 times larger. The 5 % ISRS of the roof mass peaks at 4.27 g (X), 4.25 g (Y) and 1.07 g
  for the Y response to the X input.
* RELDISP: the peak roof drift relative to the ground in Y is **5.2 mm**.

The three directional responses are combined afterwards (SRSS or 100-40-40). This
post-processing is outside the SASSI modules.

## Example 6: near-field soil nonlinearity (loose backfill behind a wall)

**Model.** A 0.6 m reinforced-concrete wall (SHELL, 4 m long, 2 m high, its top at the ground surface)
carries 150 t of a bridge deck at its top edge, like an abutment. Behind it (x = 0 to 4 m) the soil is a
4 m x 4 m x 2 m block of loose backfill (Vs 160 m/s at low strain) in a 10 m deposit of medium-dense
sand (Vs 250 m/s) on rock (Vs 800 m/s). SSI has two kinds of soil nonlinearity (manual §1.5.4):
* **primary** (free field): SOIL iterates the layer properties of the site. Part 1 runs SOIL alone
  with the SHAKE91 Sand curves (`INP,../sassi/data/dynp_library.pre`, `SPRO` per sublayer) and the
  rock-outcrop motion scaled to 0.30 g (`SOILX`); `SACC` saves the surface motion `ACC001.TH`, the SSI
  control motion;
* **secondary** (near field): the backfill strains more than the free field. It is modelled with
  SOLID elements of the *structure* (group 3, ETYPE 1) at the nodes of the excavated soil, and their
  shear modulus and damping are iterated element by element.

The strain-compatible layers of FILE88 are copied into L 11-16 (each embedment layer has its own L
number, and the excavated groups use it with `MACT`). Every excavated node is an interaction node (FV).
The backfill material `M,2,320,160,18.5,0.02,0.02,3` holds the **low-strain** properties (G_max = 48 MPa):
the curve softens it, so it must not be given strain-compatible values. `PIN,0.65` sets the effective
strain factor, `PINGRP,3,1,0.4,1.0,Sand` declares group 3 nonlinear (octahedral strain, curve Sand,
start at 0.4 x the free-field G: the native soil has G = 116 MPa, so the start is close to the backfill's
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
* SOIL: rock outcrop 0.30 g → 0.58 g at the surface; the native Vs falls from 245 to 148 m/s and the
  damping rises from 1.4 % to 10.3 % from the top to the bottom of the sand.
* The iterations converge in **6 passes**: max |dG/G| = 66, 22, 14, 9.4, 4.0, 2.9 and 1.5 % (iterations
  0 to 6). The converged backfill has G/G_max = 0.52-0.57 and 6.4-7.2 % damping (effective octahedral
  strain 0.041-0.053 %); the elements next to the wall (x = 0 to 2 m) strain and soften more than those
  at the back.
* X input: the wall top (node 22, carrying the deck) reaches **0.75 g**, the back of the backfill
  (node 24) 0.58 g, the free-field value.

## Example 7: Option NON, a shear-wall building with cracking walls

**Model.** A 12 m x 12 m two-storey box (storeys of 4 m) with 0.3 m reinforced-concrete walls (meshed
3 m x 2 m) on a 1.5 m stiff surface mat; 0.4 m floor and roof slabs carry 400 t of equipment each (`MT`,
`MTGEN`). Site: 10 m of stiff soil (Vs 400 m/s) in 25 sublayers on rock. Input: the RG 1.60 record
scaled to 0.6 g in X and Y and 0.4 g in Z. Each wall of each storey is a SHELL group with a material of
its own (`M,11` ... `M,18`): NONLINEAR changes E panel by panel (D-NON-12). The eight panels are low-rise
(h/l = 1/3) and shear governed.

The Option NON data: `PNLGEN` makes one panel per vertical shell group (Disp. Opt 1 = shear strain from
the four corner nodes, Force Opt 1 = Cheng-Mertz shear hysteresis); `SHEAR,0,30000,420000,0.005,0` lists
the ACI 318-08, Wood, Barda and Gulec-Whittaker capacities (f'c = 30 MPa, f_y = 420 MPa, 0.5 % web
steel); `BBCGEN,0,1,...,0.3` builds the 22-point backbone of every panel from the ACI capacity V_u
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
* Convergence after **7 iterations**: max |dE/E| = 28.0, 22.9, 17.1, 11.5, 6.6, 4.2, 2.4 and 1.3 %
  (analyses 0 to 7).
* The four storey-1 walls crack: E/E_el = **0.35**, damping 16.7 % (4 % elastic + 12.7 % hysteretic),
  peak shear strain 3.7e-4, i.e. mu = 4.4 times the cracking strain, F_mu = 3.3. The storey-2 walls stay
  elastic (mu = 0.84, E/E_el = 1, 4 %).
* Converged model, X input: peak acceleration 0.64 g at the mat centre (node 13), 0.85 g at the floor
  centre (node 63) and 1.12 g at the roof centre (node 113).

Try `EQL,0.8,1,7,0,1` (damping cut-off 7 %, the ASCE 4 level for cracked concrete) or BBCGEN with
`<ShearModel>` 4 (Gulec-Whittaker) and see how the convergence changes; `docs/user/OPTION_NON.md`
sections 8 and 11 discuss both.

## Example 8: an embedded shear-wall building, FV against the subtraction methods

**Model.** A reinforced-concrete shear-wall building, 24 m x 24 m in plan, embedded 8 m. Shear walls on
x, y = -12, -6, 0, 6, 12 m make sixteen rooms of 6 m x 6 m below and above grade; the main block rises
two storeys of 5 m and its four central rooms continue as a tower two storeys higher (roof at 20 m).

| Member (group, colour in the 3D view) | Thickness | Elements |
|---|---|---|
| basemat at z = -8 m (5, purple) | 2.0 m | 64 SHELL (3 m x 3 m) on the bottom face of the excavation |
| outer basement walls (12, brown) | 1.0 m | 128 SHELL (3 m x 2 m) on its lateral faces |
| interior walls on x, y = -6, 0, 6 m, basement, main block and tower (10, teal) | 0.6 m | 304 SHELL |
| floor slabs: basement z = -4 m, grade, z = 5 m, tower z = 10 and 15 m (9, pink) | 0.6, 0.8, 0.5 m | 221 SHELL |
| roofs: main roof z = 10 m around the tower, tower roof z = 20 m (11, lavender) | 0.5 m | 64 SHELL |
| outer walls above grade: main block (0.8 m), tower (0.6 m) (13, cream) | | 96 SHELL (3 m x 5 m) |

The group numbers are chosen for the colours of the 3D view (group g is drawn in colour g of the element
palette). Concrete E = 30 GPa, nu = 0.2, 24 kN/m3, 4 % damping (ASCE 4 response level 1 as the manual's
Damping Cutoff note quotes it: uncracked stiffness, 4 %); slabs and roofs carry 4 kN/m3 more for
distributed equipment, piping and live load. A 3 m x 3 m stair opening passes the basement slab, the grade
slab and the first floor (`EDEL`, then `ECOMPR`, because CHECK treats a missing element number as an error).
Ten equipment items of 25 to 150 t (810 t) sit at room centres (`MT`). CALCM: 1.545E5 kN of elements,
16,556 t with the equipment; the excavated soil weighs 8,925 t.

The site: 8 m of sand and gravel (Vs 300 m/s, four 2 m embedment layers with their own L numbers), 12 m
of dense gravel (Vs 450 m/s, 1 m sublayers), 12 m of very dense gravel (Vs 650 m/s, 2 m sublayers) on
rock (Vs 1500 m/s): 22 TOPL layers. The excavated soil is 256 SOLID elements (3 m x 3 m x 2 m) in four
groups, one per embedment layer. The mesh passes Vs/(5h) = 300/(5 x 3) = 20 Hz in plan and 30 Hz
vertically, the cut-off of the 41 SSI frequencies (0.1 Hz, then every 0.5 Hz from 0.5 to 20 Hz). The
input is the RG 1.60 record as a vertically incident SV wave (control motion in X at the surface).

The building touches the excavated soil only on the foundation-soil interface (basemat and outer
walls). The interior walls and the two slabs inside the basement have nodes of their own at the
coordinates of the excavated-soil nodes (manual 1.5.1 rule 11; `EXCSTRCHK` finds no shared interior
node): they are generated as "excavation node + 1000", then `RMVUNUSED` and `NCOM` number the 781 nodes
without gaps (1-405 excavation, 406-567 main block, 568-617 tower, 618-781 basement interior). `FIXROT`
restrains the shell drilling rotations of the shell-only nodes, `D` in `FOREACH` loops those of the
face-interior nodes shared with the soil. `WINDOWSETTINGS,HIDEGROUP` hides the four soil groups in the 3D
views (display only), so the element plot and the gallery picture show the building closed, with the
interaction nodes on its basement walls; `WINDOWSETTINGS,SHOWGROUP,1` ... shows the soil again and
`WINDOWSETTINGS,VOLUME,,,1` cuts the building open north of y = 0 (lesson 11). `INTGEN` builds the
interaction sets on the same model and `CPMODEL` / `FCOPY` reuse the free field and the point-load
solutions, as in example 2:
* model 0 `ex08`: FV, 405 interaction nodes (1215 DOFs);
* model 2 `ex08fsin`: FI-FSIN (subtraction method, SM), 209 nodes;
* model 3 `ex08evbn`: FI-EVBN (modified subtraction method, MSM), 258 nodes.

FFV (`INTGEN,5`, 307 nodes) is not run here: on this model it costs as much as FV (lesson 11 has it as
a Try this).

**Look at.** The FILE8 of each folder, or the transfer functions `00605TR_X.TFU` (tower roof centre),
`00567TR_X.TFU` (main roof corner), `00757TR_X.TFU` (grade slab centre), `00041TR_X.TFU` (basemat centre),
`00045TR_Z.TFU` (basemat edge, rocking), `00739TR_Z.TFU` (grade slab under a 150 t item) and
`00365TR_X.TFU` (the excavated soil at the centre of the excavation top face, inside the basement); the
5 % ISRS `000nnTR_X01.RS` and `000nnTR_Z01.RS`; the wall histories `SHELL_012_00004_FXY.THS` (south wall,
in-plane shear) and `SHELL_012_00069_MYY.THS` (west wall, bending); `ex08fsin.err` (EDU-12); the ANALYS
listings (interaction DOFs, memory estimate, time per frequency).

**Model checks.**
* CHECK: no error and no warning for FV; FI-FSIN and FI-EVBN only EDU-12 (Non-FV Method Selected:
  Validate Against FV). `EXCSTRCHK`, `FIXEDINT` and `HINGED` find nothing.
* At 0.1 Hz every transfer function is 1.000 (ANALYS low-frequency check: largest deviation 0.033 %).
* Fixed-base frequencies (an eigenvalue analysis of the HOUSE structural matrices `COOSK`/`COOSM`, done
  outside the program): with the basemat clamped the first lateral modes are at 12.4 Hz (Y) and 12.5 Hz
  (X), 65 % of the mass each, the next lateral ones at 21-23 Hz, the first strong vertical mode at
  31 Hz; with the whole basement clamped, 16.4 Hz. These are in the range of squat shear-wall buildings.

**Expected results** (`tests/integration/test_examples.py`, marked `slow`):
* Run time about 60-90 s on a laptop: ANALYS 26 s (FV), 21 s (FI-FSIN) and 25 s (FI-EVBN) for 41
  frequencies, measured with other applications loading the machine. The dense impedance is 23.6 MB per
  matrix for FV and 6.3 MB for FI-FSIN, but on a model this small the 3930 equations of the structure
  and the excavation cost as much as the impedance: FV takes only 1.25 times as long as FI-FSIN. The
  manual's tens of times apply to excavations with thousands of nodes. About 57 MB in the three folders.
* FV: the tower roof X transfer function peaks at **2.82 at 6.5 Hz** (the SSI frequency, about half the
  fixed-base one), the main roof corner at 1.71; the basemat never moves more than the free surface.
  5 % ISRS zero-period accelerations: 0.280 g at the basemat (kinematic interaction: below the 0.324 g
  control motion), 0.320 g at grade, 0.364 g at the main roof and 0.464 g at the tower roof; the tower
  roof ISRS peaks at 2.27 g at 6.8 Hz. South basement wall FXY 374 kPa, west wall MYY 41.6 kN m/m
  (bottom panels).
* **FI-FSIN** follows FV within 2.4 % (complex difference, relative to the FV peak) below 10 Hz and
  within 6.7 % from 10 to 14 Hz, then has a **spurious resonance at 15.5-16 Hz**: basemat 0.675 against
  0.325 and grade slab 0.545 against 0.338 at 15.5 Hz, tower roof 1.223 against 0.598 at 16 Hz, and the
  soil inside the basement (node 365) 12.6 against 1.6; a dip precedes it (the main roof corner 0.011
  against 0.095 at 14.5 Hz). `CRITFREQ` flags the interpolated SM peaks near 15.7 Hz as not supported by
  computed values. With 0.25 Hz steps (a scratch run) the computed peak is at 15.75 Hz: basemat 2.10,
  tower roof 3.33, enclosed soil 41.7, against 0.31, 0.61 and 1.8 with FV.
* The 5 % ISRS of FI-FSIN differ from FV by at most 4.9 % below 10 Hz, are **up to 31 % too high at
  15.5 Hz** (basemat X +25 %, grade slab X +23 %, grade slab Z under the 150 t item +31 %; the tower roof
  only +3 % / -9 %) and differ by up to 15 % above 17 Hz; the peak wall forces change by about 1 %
  (FXY -1.2 %, MYY +0.9 %): they come from the 6.5 Hz SSI response.
* **FI-EVBN** follows FV within **3.7 %** at every computed frequency and output (the largest difference
  is the grade slab Z under the 150 t item at 20 Hz; X within 1.3 %), the ISRS within 1.0 %, the wall
  forces within 1.2 %. FFV (lesson 11, Try this) gives the same 3.7 % at the same place, X within 1.5 %.

**Why the subtraction method resonates.** In FI-FSIN the 196 excavated nodes inside the basement and
on the excavation top face carry neither impedance nor free-field load: their rows keep only
-(K_e - w2 M_e). No structure rests on them (the interior structure has its own nodes), so the anomaly
comes from the soil alone: the enclosed soil, held at the walls and the basemat and free at the top,
has natural frequencies of its own, and near the first one the subtracted soil resonates through the
interface. An eigenvalue analysis of the HOUSE matrices `K_e`, `M_e` with the interaction DOFs fixed
(done outside the program) gives 15.45 Hz for the FI-FSIN set, 22.35 Hz for FI-EVBN (above the 20 Hz
cut-off: no anomaly) and 37.7 Hz for FFV; they depend on the excavation only, not on the building. How
much of the anomaly reaches the ISRS depends on the building and on the input, not only on the soil:
here the floor spectra take 20-31 % of error at 15.5 Hz. The frequency scales with the soil velocity. The same model with other
embedment soils (L 1-4 changed, everything else equal; scratch runs):
* Vs 600 m/s (stiff): FI-FSIN follows FV within 1.3 % up to 20 Hz, its ISRS within 1.3 %; FI-EVBN
  within 0.9 %. This is the manual's statement that FI-FSIN agrees with FV on stiff soil and rock sites.
* Vs 200 m/s (soft): FI-FSIN resonates at 10.5 Hz (basemat 1.31 against 0.57), ISRS +16 % / -16 %;
  FI-EVBN departs from FV by up to 24 % near 19 Hz (its own enclosed-soil frequency falls to about
  15 Hz), ISRS +3.1 % / -9.3 %. (The 3 m mesh passes only 13.3 Hz for this soil: the comparison of the
  methods holds, the absolute results above 13 Hz do not.)

Lesson 11 of the guided course builds this example step by step, shows the soil and the interaction
nodes and a cutaway of the building, animates the FV and FI-FSIN motion at 16 Hz (`HARMFRAME`) and
discusses the validation a design engineer must document.

## Tests

* `tests/integration/test_examples.py` runs examples 1-5 and 8 in a temporary copy, through the
  interpreter as `sassi --cwd examples run` does. It asserts that every module finishes with
  status OK, checks the key results above (with physically motivated bounds), checks that
  `WRITE` -> `INP` gives back the same model for every model of every example, and checks the
  documented `--cwd` command line. The example 4 and 8 tests are marked `slow` (EQUAKE; three ANALYS
  runs).
* `tests/unit/test_nlsoil_example.py` runs example 6: a clean run, the convergence history, the
  wall-top spectrum and the `WRITE` -> `INP` round trip.
* `tests/unit/test_nonlinear_example.py` runs example 7 (marked `slow`): a clean run, the panels and
  backbone curves, the convergence, the cracked first storey, the next HOUSE deck and the round trip.
* `tests/verification/test_vp_examples.py` runs **VP-E1**.

```bash
.venv/bin/python -m pytest -q tests/integration tests/verification/test_vp_examples.py
.venv/bin/python -m pytest -q tests/unit/test_nlsoil_example.py tests/unit/test_nonlinear_example.py
```
