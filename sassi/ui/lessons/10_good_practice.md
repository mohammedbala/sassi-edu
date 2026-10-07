---
id: 10-good-practice
title: Good practice, checks and the limits of SSI models
part: Advanced
order: 10
minutes: 40
example: ex01_surface_stick
summary: The checks that make an SSI result trustworthy: model checks before the run, the low-frequency check, computed against interpolated transfer functions and CRITFREQ, the cut-off frequency, model size and method choice, how SASSI-EDU is verified, and what it does not do.
objectives: [Run CHECK and the model checks and read their findings, Use the low-frequency check and TFU against TFI to judge a frequency set, Detect a missing SSI frequency with CRITFREQ and fix it, See the cut-off frequency act as a low-pass filter on ISRS, Choose the interaction-node set and estimate the cost of a model, Know how SASSI-EDU is verified and where its limits are]
prerequisites: [04-surface-ssi, 05-embedded, 07-three-components]
---
An SSI analysis can run to the end and still be wrong: a missing frequency, a cut-off that filters out
the band of your equipment, a fixed interaction node, a hinge in the mesh. None of these stops a
module; all of them change the ISRS. A reviewer of a nuclear SSI calculation looks for exactly these
checks, and the ACS SASSI manual's "Engineering Considerations" (paraphrased in
[User Guide §14](docs/user/USER_GUIDE.md#14-engineering-guidance)) lists them.

This lesson goes through them on example 1 (the stick on a surface mat). Opening the lesson runs the
complete example (about 4 s); two experiments then run modified copies of it, one with a frequency set
that misses the SSI peak and one with a cut-off frequency of 10 Hz instead of 20 Hz. The steps run
in about 10 s.

```sassi-setup
INP,ex01_surface_stick.pre
```

## Before you run: CHECK and the model checks

`CHECK` checks the model and every module enabled by `AOPT` against the error and warning list of the
manual (Chapter 10) plus SASSI-EDU's EDU-nn checks. The model-check commands look for the modelling
errors that CHECK cannot see from the input alone.

```sassi
CHECK
* interior excavation nodes shared with the structure
EXCSTRCHK
* interaction nodes with a fixed translation
FIXEDINT
* beams or shells meeting solids or shells at one node (unintended hinges)
HINGED
* beam orientation (K) nodes that are interaction nodes
KINT
* nodes held only by springs
FREESPRING
```

### What this does
* `CHECK` writes `ex01/ex01.err` and reports per module: here 0 errors and 0 warnings for the model,
  SITE, POINT, HOUSE, ANALYS, MOTION, STRESS and RELDISP. A module with an error gets no input deck
  from AFWRITE and cannot run.
* `EXCSTRCHK`, `FIXEDINT`, `KINT` and `FREESPRING` find nothing in this model.
* `HINGED` reports one *possible* hinge: node 41, where the stick (beam 1 of group 2) meets the mat
  shells out of their plane, so the shells do not transmit the rotation about their normal
  (drilling).

### Why it matters
HINGED is right that a flat shell has no drilling stiffness, but it looks at the stick and the shells
only: the rigid "spider" of group 3 connects node 41 to the eight surrounding mat nodes, and its beams
carry the stick's base torsion and moments into the mat. That is exactly the remedy the message suggests
("penetrate ... with massless beams or shells"). Read every warning and decide; do not silence it.
The other checks catch errors that change results without stopping a run:

| Check | Typical cause | Consequence |
|---|---|---|
| EXCSTRCHK | an internal basement wall or slab sharing interior excavation nodes | the excavated soil is coupled to the structure where it should not be: a serious modelling error |
| FIXEDINT | a support condition left on an interaction node | the soil impedance at that node is short-circuited |
| HINGED | a beam or shell meeting solids at one node | a mechanism or an artificial pin at the joint |
| KINT | a beam orientation (K) node placed on an interaction node | incorrect results when the node only orients beams (manual) |

### Technical basis
CHECK applies the manual's Errors 1-128 and Warnings 1-11 and the EDU-nn checks
([User Guide §13](docs/user/USER_GUIDE.md#13-common-errors-and-how-to-fix-them) lists the frequent
ones with their fixes); the model checks are described in
[User Guide §5.6](docs/user/USER_GUIDE.md#56-model-checks) and verified with planted faults by
[VP-52](docs/verification/VERIFICATION_MANUAL.md#vp-52). Shell drilling rotations need `FIXROT`
(EDU-06), which example 1 applies only at the nodes connected to shells alone.

### In ANSYS terms
The equivalent of the element shape checks, `CHECK` of the solution and a look at the reaction
forces for unintended supports, before you trust a model.

```action
open-file: ex01/ex01.err
open-dialog: CHECK
plot-nodes
explain: HINGED
```

## The low-frequency check and computed against interpolated transfer functions

ANALYS solves the SSI equation only at the 22 SSI frequencies of set 1; MOTION and STRESS interpolate
the transfer functions to the 820 Fourier frequencies up to the cut-off. Two checks belong to every
run: the low-frequency check of the ANALYS listing, and the comparison of the computed (`.TFU`) and
interpolated (`.TFI`) transfer functions. `CRITFREQ` automates the second.

```sassi
LFREQ
* CRITFREQ,<tol %>,<minfilter %>,<TF file without extension>,<variable>
CRITFREQ,5,50,00085TR_X,ADDF
* the same check on the mat centre
CRITFREQ,5,50,00041TR_X,ADDM
```

### What this does
* `LFREQ` lists the frequency numbers of set 1 with their values in Hz: 0.098 Hz, then every 0.5 Hz
  from 0.49 to 6.0 Hz, then 7 to 10 Hz in 1 Hz and 12 to 20 Hz in 2 Hz steps.
* `CRITFREQ,5,50,00085TR_X,ADDF` looks at the peaks of the interpolated roof transfer function
  `00085TR_X.TFI` that are higher than 50 % of its maximum, compares each with the larger of the two
  computed values that bracket it in `00085TR_X.TFU`, and puts the frequency numbers of the peaks that
  differ by more than 5 % into the variable ADDF. Here: one peak, 3.47 Hz (number 142),
  $\lvert\mathrm{TFI}\rvert = 13.11$, supported by the computed value 13.02 at 3.49 Hz (number 143,
  difference 0.7 %); ADDF stays empty.
* On the mat centre (`00041TR_X`) CRITFREQ finds two peaks. The one at 8.89 Hz is supported (0.0 %);
  the one at 3.37 Hz (number 138) is not: $\lvert\mathrm{TFI}\rvert = 1.70$ against computed
  neighbours of at most 1.42, 20.2 %, so ADDM = 138. This is the unsupported mat peak that lesson 04
  found by eye.

### Why it matters
* **Low-frequency check.** As $f \to 0$ the whole system moves with the free field, so every transfer
  function in the input direction tends to 1. The ANALYS listing reports the largest deviation at the
  first frequency (0.098 Hz): 0.103 % at node 85. A large deviation points to a modelling error:
  interaction nodes missing or not connected to the structure, a mechanism, a wrong control direction.
* **TFU against TFI.** An interpolated peak much higher than its computed neighbours means the
  frequency set does not resolve that peak; a spike between two computed points can be an artefact.
  Either way the ISRS at that frequency rests on interpolation, not on a solution. The check is per
  node and direction: the same frequency set resolves the roof but not the mat centre. Lesson 04's
  Try this added five frequencies around 3.4 Hz, among them number 139 (3.39 Hz) next to the flagged
  138, and computed 1.696 there, the value the interpolation predicted: the interpolation was right
  this time, but only the added computed points show it. The next two steps make this check routine.

### Technical basis
MOTION and STRESS fit, in windows of five consecutive SSI frequencies, the transfer function of a
two-degree-of-freedom system with hysteretic damping (Tajirian 1981),

```math
H(\omega) = \frac{C_1\,\omega^4 + C_2\,\omega^2 + C_3}{\omega^4 + C_4\,\omega^2 + C_5}
```

with five complex constants $C_1, \dots, C_5$ per window, and reproduce the computed values exactly
([Theory §10](docs/theory/THEORY_MANUAL.md#10-transfer-function-interpolation); verified by
[VP-28](docs/verification/VERIFICATION_MANUAL.md#vp-28) and, for CRITFREQ, VP-53). The low-frequency
check is G-19 ([Theory §9.4](docs/theory/THEORY_MANUAL.md#94-simultaneous-cases-and-the-low-frequency-check),
VP-41). Recommended numbers of SSI frequencies: 40-80 for stick models, 100-300 for large FE models,
at least 200 for incoherent analyses, at most 500 per ANALYS run
([User Guide §9.2](docs/user/USER_GUIDE.md#92-choosing-the-ssi-frequencies)).

```action
open-listing: ANALYS
plot-spectrum: ex01/00085TR_X.TFU, ex01/00085TR_X.TFI
plot-spectrum: ex01/00041TR_X.TFU, ex01/00041TR_X.TFI
explain: CRITFREQ,5,50,00085TR_X,ADDF
```

## Experiment 1: a frequency set that misses the SSI peak

A copy of the model (model 1, directory `ex01c`) with 11 frequencies spaced 1 to 4 Hz apart and
none at the SSI frequency of the stick (3.47 Hz). SITE, POINT, HOUSE, ANALYS and MOTION run again.

```sassi
CPMODEL,1
ACTM,1
MDL,ex01c,../ex01c
* 0.1, 1, 2, 3, 4, 5, 7, 9, 12, 16 and 20 Hz
FREQ,1,0
FREQ,1,4,41,82,123,164,205,287,369,492,655
FREQ,1,819
AOPT,0,0,0,1,1,1,0,0,1,0,1,0,0,0
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
RUNMOTION
CRITFREQ,5,50,00085TR_X,ADDC
* keep this run's roof transfer functions
FCOPY,00085TR_X.TFU,COARSE_00085TR_X.TFU
FCOPY,00085TR_X.TFI,COARSE_00085TR_X.TFI
```

### What this does
`CPMODEL,1` copies example 1 to model 1 and warns that the copy still has the original's name and
directory; `ACTM,1` and `MDL,ex01c,../ex01c` give it its own (lesson 06). `FREQ,1,0` deletes set 1;
the two `FREQ` lines define the coarse set. `AOPT` enables SITE to MOTION
(STRESS and RELDISP are not needed here). `CRITFREQ` checks the new roof transfer function and stores
the flagged frequency numbers in ADDC; `FCOPY` keeps the result before the next step overwrites it.

### Why it matters
With 11 frequencies the computed roof transfer function never sees the resonance: its largest computed
value is 4.69 at 3.0 Hz. The interpolated one still peaks at 13.12 at 3.47 Hz (13.11 with the full set
of example 1), and the 5 % roof ISRS hardly changes (8.04 g against 8.03 g, at 3.55 Hz): the two-degree-of-
freedom interpolant reconstructed the resonance from points on its flanks. CRITFREQ flags it anyway,
interpolated 13.12 against computed neighbours of at most 4.69 (180 %), and stores 142 in ADDC. It is
right to: a peak almost three times higher than anything computed is a prediction, not a result, until a
computed point confirms it.

### Technical basis
The interpolant is exact for a two-degree-of-freedom hysteretic system; a stick on a mat is close to
one around its SSI mode, which is why the peak was recovered. Nothing guarantees it for a building
with several close modes, soil-column resonances and coupled directions. The CRITFREQ criterion:
flag a TFI peak at $f_p$ when

```math
100\,\frac{\lvert A_I(f_p) - A_\text{ref}\rvert}{A_\text{ref}} > \text{tol}
```

with $A_I(f_p)$ the interpolated amplitude at the peak, $A_\text{ref}$ the larger computed amplitude
of the two bracketing SSI frequencies and tol the tolerance in percent (`<tol %>`).

```figure
freq-interpolation set=coarse modes=2
This experiment with a stand-in transfer function. On the coarse set no computed point is near
the SSI peak (about 3.5 Hz), yet the two-degree-of-freedom interpolant recovers the peak, and
CRITFREQ flags it because no computed value supports it; add the flagged frequency and the flag
clears. With three modes the interpolant misses the close mode, and it takes a few rounds.
```

```action
plot-spectrum: ex01c/COARSE_00085TR_X.TFU, ex01c/COARSE_00085TR_X.TFI
explain: FREQ,1,4,41,82,123,164,205,287,369,492,655
```

## Add the flagged frequency and confirm

CRITFREQ left the frequency number to add in the variable ADDC. Add it to the set, run again and
check again.

```sassi
* append the flagged frequency number (142 = 3.47 Hz) to set 1
FREQ,1,@ADDC[1]
LFREQ
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
RUNMOTION
CRITFREQ,5,50,00085TR_X,ADDC
```

### What this does
`@ADDC[1]` is the first item of the variable written by CRITFREQ (frequency number 142). `LFREQ`
shows the 12 frequencies in ascending order. After the new run CRITFREQ finds the peak supported by a
computed value: 13.11 at 3.47 Hz, difference 0.0 %, nothing flagged.

### Why it matters
This is the loop the ACS SASSI workflow prescribes (initiation run, review of the computed against the
interpolated transfer functions, added frequencies, merge): you add frequencies until every peak that
matters for the ISRS is supported by a computed point. For a large model you do not rerun everything:
run SITE and POINT for the new frequencies (they must exist in FILE1 and FILE3), ANALYS for the new
frequencies only, and merge the two FILE8 with COMBIN:

```sassi-show
FMOVE,FILE8,FILE81
* ... frequency set with the added numbers only, then AFWRITE and RUNANALYS
FMOVE,FILE8,FILE82
* FILE81 + FILE82 -> FILE8 (duplicate frequencies are an error)
RUNCOMBIN
```

### Technical basis
[User Guide §9.4](docs/user/USER_GUIDE.md#94-adding-frequencies); COMBIN is verified by
[VP-25](docs/verification/VERIFICATION_MANUAL.md#vp-25) (merging odd and even frequency sets gives
the complete FILE8 exactly).

### Check yourself
Why does the roof transfer function need many frequencies around 3.5 Hz but few between 10 and 20 Hz?

Answer: the response changes fast only near resonances: around the SSI frequency of the stick
(3.47 Hz) the computed amplitude rises from 4.7 at 3.0 Hz to 13 and falls back to 4.3 at 4.0 Hz, so
the interpolation needs computed points close to the peak. Between 10 and 20 Hz the computed roof
transfer function varies smoothly between 0.76 and 0.88 (`ex01/00085TR_X.TFU`), and a few points
describe it. Put frequencies where the structure's and the soil column's resonances are,
and verify with CRITFREQ.

```action
plot-spectrum: ex01c/COARSE_00085TR_X.TFI, ex01c/00085TR_X.TFU, ex01c/00085TR_X.TFI
```

## Experiment 2: the cut-off frequency is a low-pass filter

The last SSI frequency is the **cut-off**: above it every transfer function is zero. Model 2 repeats
example 1 with the frequency set stopped at 10 Hz instead of 20 Hz.

```sassi
ACTM,0
CPMODEL,2
ACTM,2
MDL,ex01f,../ex01f
* the frequencies of example 1 up to 10 Hz (numbers 4 ... 410)
FREQ,1,0
FREQ,1,4,20,41,61,82,102,123,143,164,184
FREQ,1,205,225,246,287,328,369,410
AOPT,0,0,0,1,1,1,0,0,1,0,1,0,0,0
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
RUNMOTION
ACTM,0
```

### What this does
`ACTM,0` returns to example 1 and `CPMODEL,2` copies it; the copy keeps the 17 frequencies of
example 1 up to number 410 (10.0 Hz) and runs in its own directory `ex01f`. `ACTM,0` at the end makes
example 1 the active model again.

### Why it matters
The SSI peak is untouched (roof ISRS 8.03 g at 3.55 Hz in both runs) and the zero-period
accelerations barely move (mat centre 0.356 and 0.364 g, roof 1.281 and 1.282 g at 100 Hz). Between the
cut-off and the rigid range, however, the spectra drop: the 5 % ISRS of the mat centre falls at
10 Hz from 0.96 g to 0.79 g (-18 %), at 12 Hz, the largest drop, from 0.82 g to 0.48 g (-41 %), and
at 15 Hz from 0.51 g to 0.41 g (-20 %); the roof at 10 Hz from 1.77 g to 1.63 g. Below about 8.7 Hz
the mat-centre and roof spectra of the two runs agree within 5 %. Equipment on the mat with frequencies of 10-15 Hz would be qualified for too little, and
neither the peak nor the ZPA would warn you. The transfer-function plot shows why: the TFI of the
10 Hz model stops at 10 Hz.

### Technical basis
Above the last SSI frequency $f_N$ MOTION sets every transfer function to zero
([Theory §10.2](docs/theory/THEORY_MANUAL.md#102-window-schemes-options-0-6)), so $f_N$ is the corner of
an ideal low-pass filter on every result. Its choice (User Guide §14 item 9 and spec
[03 §2.1](docs/spec/03_guidelines.md)): the frequency content of the control motion, the structural
modes (the 90 % cumulative modal mass criterion of ASCE 4 Section 3), typically 30-40 Hz for soil and
60-70 Hz for rock sites, and never above the frequency the mesh passes: soil sublayers and excavation
elements no thicker than $V_s/(5 f_\text{cut})$ (example 1: 1.6 ft sand sublayers pass 125 Hz, 3.25 ft
gravel sublayers 102 Hz; EDU-10 warns otherwise).

### Check yourself
Your equipment is qualified in the 10-30 Hz band and the SSI analysis stopped at 15 Hz to save run
time. Which part of the ISRS can you use?

Answer: only the part well below 15 Hz, after a check. Between the cut-off and the rigid range the
ISRS is built from the motion below the cut-off: it is not conservative in general (here a 10 Hz
cut-off lowered the mat ISRS by 18 % at 10 Hz and 41 % at 12 Hz while the ZPA hardly changed, and
only below about 8.7 Hz, 0.87 times the cut-off, did the mat spectra of the two runs agree within 5 %). Extend the frequency set
(and check that the soil layers and the mesh pass the new cut-off) or demonstrate by a sensitivity run
that the band is not affected.

```action
plot-spectrum: ex01/00041TR_X02.RS, ex01f/00041TR_X02.RS | log
plot-spectrum: ex01/00085TR_X02.RS, ex01f/00085TR_X02.RS | log
plot-spectrum: ex01/00085TR_X.TFI, ex01f/00085TR_X.TFI
```

## Model size, method choice and the advanced options

The cost of an SSI run is set by the number of interaction nodes: the soil impedance is a dense matrix
on their degrees of freedom, inverted at every frequency.

```sassi
INTCOUNT
```

### What this does
`INTCOUNT` prints the number of interaction nodes and the memory estimate of ANALYS: here
$N_\text{int} = 81$ nodes (243 DOFs), a dense complex matrix of
$(3N_\text{int})^2 \times 16\,\text{B} = 0.00094\,\text{GB}$, about 0.0028 GB for ANALYS.

### Why it matters
For 2,000 interaction nodes the matrix alone needs 576 MB, for 10,000 it needs 14.4 GB, and the run
time grows with the cube of the node count. The levers:

* **Interaction-node set** (lesson 05). **FV** (every excavated node) is the reference. FI-FSIN
  (subtraction method) and FI-EVBN (modified subtraction) use the faces of the excavation, FFV a
  reduced volume; ASCE 4-16 and SRP 3.7.2 require a validation against FV before such a method is used
  in production (the ACS SASSI manual warns that quarter models can hide the instabilities of the
  subtraction method). Example 2 (lesson 05) shows FI-FSIN with a spurious 11 Hz resonance and
  FI-EVBN within 2.4 % of FV.
* **Symmetry planes** (`SYMM`): a half or quarter model of a symmetric structure under symmetric or
  antisymmetric loading, one run per input direction; not with incoherency, wave passage or the global
  impedance.
* **Frequency-split runs** (`AFWRBAT`) and **restarts** (New Structure, New Seismic Environment)
  reuse or distribute the expensive part.
* **2D plane-strain models** are for sensitivity studies: 2D SSI overestimates radiation damping.

Advanced input options change the answer: **incoherency** reduces the high-frequency
translational response of large stiff foundations (typically on rock) and adds rotations; the
stochastic simulation is the reference approach, and the deterministic AS and SRSS approaches were
validated only for stick models on rigid mats (EDU-13 warns).

```sassi-show
* half model about the plane x = 0 for the X input (antisymmetric), nodes 101, 105, 141 in the plane
SYMM,1,1,101,105,141
* incoherent motion with wave passage: 2007 hard-rock coherency, no delay, 20 stochastic samples
HOUSE,32.2,0,0,2,0,1,1,0,0
WPASS,1e9,0,5
INCOH,0.1,0.1,0.2,0.5,1,1,0,1,1975,2026,180
HOUSEX,0,0,20
```

### Technical basis
$(3N_\text{int})^2 \times 16$ bytes for the dense impedance
([User Guide §10.4](docs/user/USER_GUIDE.md#104-memory-and-run-time)); the method variants
([Theory §3.3](docs/theory/THEORY_MANUAL.md#33-method-variants-the-choice-of-the-interaction-set);
VP-16 checks the exactness of FV with a zero-SSI identity, VP-50 the interaction sets INTGEN
generates), symmetry ([Theory §15.2](docs/theory/THEORY_MANUAL.md#152-symmetry-and-antisymmetry-planes-symm),
VP-T2: half and quarter models reproduce the full model to 4e-11), incoherency
([Theory §14](docs/theory/THEORY_MANUAL.md#14-incoherency-wave-passage-and-multiple-excitation),
VP-26, VP-27, VP-I1, VP-I2), 2D ([Theory §15.1](docs/theory/THEORY_MANUAL.md#151-plane-strain-ssi-house-dim--1)).

```action
open-doc: docs/user/USER_GUIDE.md#14-engineering-guidance
explain: INTCOUNT
```

## How SASSI-EDU is verified, and what it does not do

Every module is verified by verification problems (VPs) that compare computed values with exact,
published or independently derived references, with stated tolerances. `VERIFY` runs them from the
command line.

```sassi
* low-frequency transfer functions = 1 at all nodes (G-19)
VERIFY,VP-41
* transfer-function interpolation exactness and ISRS from interpolated transfer functions
VERIFY,VP-28
* a command of the manual that this build does not implement
SOILREDEF
```

### What this does
`VERIFY,<VP id>` builds the problem's model, runs the modules and prints one row per criterion
(computed, reference, error, tolerance, pass/fail). VP-41 passes: over 102 nodes of its model the
largest deviation of $\lvert\mathrm{ATF}\rvert$ from 1 at the first frequency is 0.0014 (tolerance
0.05). VP-28 passes: the interpolation reproduces an exact two-degree-of-freedom hysteretic transfer
function to 2e-13 (options 0-5), the interpolated values equal the computed ones at the SSI
frequencies exactly, and the ISRS from 58 SSI frequencies differs from that of a dense solution by at
most 0.07 % (options 0-5; 1.2 % for the spline option 6 with 152 frequencies), against a 2 %
tolerance. `SOILREDEF` is accepted and answers "not available in this build" with a warning, like
every command of the manual that SASSI-EDU does not implement.

### Why it matters
The [Verification Manual](docs/verification/VERIFICATION_MANUAL.md) reports 74 problems, all passing,
from closed-form checks (a hysteretic SDOF, layer transfer functions, response-spectrum closed forms)
through published benchmarks (SHAKE91 sample problem, Wong's surface Green functions, Tajirian and
Tabatabaie's disk on a layer, NAFEMS plates) to independent reconstructions. Some published references
proved approximate; the lead decisions replaced them by rigorous references and keep the published
values as informative comparisons. Use the VPs to understand the methods and to check changes.

**What SASSI-EDU does not do** ([User Guide §18](docs/user/USER_GUIDE.md#18-limitations)):

* it is an educational reimplementation, **not an NQA-1 qualified program**: results for licensing
  work come from qualified software;
* not available: binary result databases, Option AA (SSI with ANSYS matrices), the New Load Vector
  restart, the Abrahamson 1993 and 2006 coherency models, anti-plane SH analysis in 2D, nonlinear beams
  and the Cheng-Mertz bending model in Option NON, water sloshing, lower-bound incoherent spectra,
  Option PRO (probabilistic site response and SSI);
* reconstructions where the manual is silent (the half-space sublayer law, interpolation option 0
  weights, hysteresis rules, the LOADGEN file layout, ...): results can differ from ACS SASSI in the
  details, and the APDL of LOADGEN has been checked by a parser, not run in ANSYS;
* performance: pure Python and NumPy/SciPy, about 2,000 interaction nodes and 20,000 DOFs per
  frequency on a laptop.

### Technical basis
[Verification Manual §1](docs/verification/VERIFICATION_MANUAL.md#1-how-the-verification-works)
explains the criteria (relative, absolute, logical, informative); the full set runs with
`.venv/bin/python -m pytest -q tests/verification` or `VERIFY,ALL`
([User Guide §16](docs/user/USER_GUIDE.md#16-verification)).

### Check yourself
A colleague shows you ISRS from an SSI run and the input file. Name five things you check before
looking at the spectra.

Answer: (1) CHECK and the model checks are clean, warnings explained; (2) the low-frequency check of
ANALYS is close to 0 % in the input direction; (3) TFU against TFI (CRITFREQ) at the floors of interest
shows no unsupported peaks; (4) the cut-off frequency covers the input content and the equipment band,
and the soil sublayers and the mesh pass it ($f_\text{pass}$ in the SITE listing, EDU-10); (5) the
interaction-node set is FV, or the reduced method was validated against FV; and the inputs: a
spectrum-compatible record that passes SRP 3.7.1, strain-compatible soil with the right control
motion (within or outcrop), the soil cases, three directions, and the combination and broadening
rules used.

```action
open-doc: docs/verification/VERIFICATION_MANUAL.md#2-summary
open-doc: docs/user/USER_GUIDE.md#18-limitations
```
