---
id: 01-why-ssi
title: From fixed base to SSI: what changes and why
part: Fundamentals
order: 1
minutes: 25
example: ex01_surface_stick
summary: Run a complete SSI analysis of a stick on a surface mat, then compare it with the same stick on a fixed base and on a massless foundation.
objectives: [Run a complete SASSI analysis and find its results, Map the SASSI modules to your seismic design workflow, Read a transfer function and an in-structure response spectrum, Measure how SSI shifts the frequency and changes the ISRS, Separate kinematic from inertial interaction]
prerequisites: []
---
You already know how to analyse a building for an earthquake on a **fixed base**: an ANSYS model
clamped at the foundation level, the design motion applied as a base acceleration, a modal or
time-history analysis, and then in-structure response spectra (ISRS), member forces and relative
displacements. That model assumes that the ground under the foundation is rigid and that the
foundation moves exactly with the free-field ground motion.

On a soil site neither assumption holds. Soil-structure interaction (SSI) has two parts:

* **Kinematic interaction**: a stiff foundation cannot follow the free-field motion point by point.
  It averages the motion over its area and depth, so even a *massless* foundation moves differently
  from the free field (the "foundation input motion").
* **Inertial interaction**: the inertia forces of the building load the soil, which deforms. The
  soil acts as frequency-dependent springs and dashpots under the foundation. The springs lower the
  frequencies of the system; the dashpots represent energy radiated into the ground (radiation
  damping) on top of the soil's material damping.

SASSI computes both in one analysis, in the **frequency domain**, by **substructuring**: the
horizontally layered site is solved semi-analytically, the building is a finite-element model, and
the two are joined at the *interaction nodes*. In this lesson you run a complete SASSI analysis of
a four-storey stick on a 64 ft × 64 ft surface mat, look at its results, and then run the same stick
on a practically rigid site (the fixed-base reference) and on a massless foundation. The numbers
you see in this lesson come from these runs.

## Run a complete SSI analysis

The example file `ex01_surface_stick.pre` builds the model, checks it, writes the module input
files and runs every module. One command runs it all.

```sassi
* read and execute the commented example file (about 4 s)
INP,ex01_surface_stick.pre
```

```action
plot-model
open-file: ex01_surface_stick.pre
open-listing: MOTION
```

### What this does
`INP,<file>` executes the commands of a `.pre` file, exactly as if you typed them; the Command
History shows each one. This file holds 80 commands:

* `MDL,ex01,ex01` names the model `ex01` and makes `ex01/` its directory and the working
  directory. Every module writes its files there.
* `L` and `TOPL` define the site: 16 ft of dense sand ($V_s = 1000\,\text{ft/s}$) over 39 ft of
  gravel ($V_s = 1650\,\text{ft/s}$) over weathered rock ($V_s = 3300\,\text{ft/s}$). The model is in
  ft, kip and s: unit weights in kcf, moduli in ksf, the floor masses given as weights in kips.
* `FREQ` and `SITE` define the 22 analysis frequencies between 0.1 and 20 Hz.
* `N`, `E` and their generators build the structure: a 5 ft SHELL mat (81 nodes, ten times
  stiffer than concrete: stiff, but not rigid, as lesson 4 shows), a four-storey BEAMS stick with
  2,200 kips per floor (weights; HOUSE divides them by $g = 32.2\,\text{ft/s}^2$),
  and a rigid "spider" of beams that clamps the stick to the mat.
* `INT,1,81,1,1` makes the 81 mat nodes **interaction nodes**: the nodes where the soil acts.
* the module options (`POINT`, `HOUSE`, `ANALYS`, `MOTION`, `STRESS`, `RELD` ...), then `CHECK`,
  `AFWRITE` (one input deck per module) and the runs `RUNSITE` ... `RUNRELDISP`.

Lessons 2 to 4 go through these commands one by one. Here you only look at what comes out.

### Why it matters
A complete SSI analysis of this model takes a few seconds. The workflow, the files and the checks
are the same for a full nuclear island model; only the run time grows. The design deliverables you
know (ISRS at the floors, member forces, relative displacements) all come out of this one chain.

### Technical basis
At every analysis frequency $\omega$ SASSI solves one complex linear system, manual Eq. 2.1
([Theory §3](docs/theory/THEORY_MANUAL.md#3-flexible-volume-substructuring)):

```math
\left[(K^*_s - \omega^2 M_s) - (K^*_e - \omega^2 M_e) + X_{ff}\right] U = X_{ff}\,U'_f
```

where

* $K^*_s$, $M_s$ are the complex stiffness (damping) and the mass of the structure (HOUSE);
* $K^*_e$, $M_e$ the same for the excavated soil, none for a surface mat (HOUSE);
* $X_{ff}$ the dynamic stiffness of the layered site at the interaction nodes: springs (real part)
  and dashpots (imaginary part) (POINT + ANALYS);
* $U'_f$ the free-field motion at the interaction nodes (SITE);
* $U$ the total motion of every node per unit control motion (FILE8).

```figure
substructuring case=surface
The three parts of Eq. 2.1 for this model. The free field is the layered site of SITE and POINT,
which the structure sees only through $X_{ff}$ and $U'_f$ at the 81 interaction nodes; a surface mat
replaces no soil, so $K^*_e = M_e = 0$.
```

The verification problem VP-E1 runs this very example and checks that the decks written by the
interpreter are exactly the intended ones (FILE8 identical to 1e-12)
([Verification Manual, VP-E1](docs/verification/VERIFICATION_MANUAL.md#vp-e1)).

### In ANSYS terms
`INP` is `/INPUT`; the `.pre` file is your APDL input file. The difference is that SASSI's
analysis runs as separate modules that exchange files, not inside one solver.

## The module chain is your seismic workflow

Each `RUN<module>` command ran one program. Together they do what you do by hand around a
fixed-base analysis: site response, soil springs, the structural model, the dynamic analysis and
the post-processing.

```sassi
* global model information and the size of the soil impedance matrix
STATUS
INTCOUNT
```

```action
open-listing: SITE
open-listing: POINT
open-listing: HOUSE
open-listing: ANALYS
open-doc: docs/user/USER_GUIDE.md#22-the-module-chain-manual-fig-11
```

### What this does
`STATUS` lists the model: 90 nodes of which 81 are interaction nodes, 3 groups with 76 elements, 4
lumped masses, one frequency set of 22 frequencies with $\Delta f = 0.0244\,\text{Hz}$. `INTCOUNT`
counts the interaction nodes (81, i.e. 243 interaction degrees of freedom) and estimates the memory
of the dense impedance matrix. The listings (buttons) are the printed output of the modules:

| Module | Computes | Writes | In your fixed-base workflow |
|---|---|---|---|
| SITE | free-field motion of the layered site for a unit control motion, at every layer interface and frequency | FILE1, FILE2 | site response analysis (SHAKE) |
| POINT | displacements of the layered site under unit point loads | FILE3 | the soil-spring calculation |
| HOUSE | stiffness and mass matrices of the structure (and of the excavated soil) | `ex01.N4`, COOSK, COOSM | building the ANSYS model |
| ANALYS | soil impedance $X_{ff}$ at the interaction nodes, the SSI equation, transfer functions of every DOF | FILE8 | harmonic analysis on springs and dashpots |
| MOTION | accelerations and ISRS from the transfer functions and the control motion | `.TFU .TFI .ACC .RS` | time-history analysis and ISRS generation |
| STRESS | element force (stress) histories | `.THS`, maxima | member forces |
| RELDISP | relative displacements between nodes | `.THD` | seismic gaps, piping and cable displacements |

### Why it matters
Because the modules exchange files, you rerun only what changed. A new earthquake record needs
only MOTION, STRESS and RELDISP; a changed structure on the same soil, with the same interaction
nodes, needs HOUSE and ANALYS (a "New Structure" restart that reuses the soil impedance) and then
the post-processing. On large models this is what makes parametric SSI studies (soil cases,
cracked and uncracked concrete) affordable.

### Technical basis
At every frequency ANALYS forms the dynamic stiffness of the structure, inverts the soil flexibility
at the interaction nodes to get the impedance $X_{ff} = F_{ff}^{-1}$ (a full complex matrix), adds
it to the structure, forms the seismic load $X_{ff}\,U'_f$ and solves
([Theory §9](docs/theory/THEORY_MANUAL.md#9-solution-of-the-ssi-equation-analys)). The impedance
matrix is dense: $(3N_{\text{int}})^2 \times 16$ bytes, about 1 MB here, 14.4 GB for 10,000
interaction nodes ([User Guide §10.4](docs/user/USER_GUIDE.md#10-running-analys)). This is why the
number of interaction nodes, not the number of structural nodes, governs the cost of an SSI model.

### Check yourself
You receive a new set of spectrum-compatible records for the same site and soil properties. Which
modules do you rerun?

Answer: only MOTION (and STRESS and RELDISP if you need forces and displacements). The transfer
functions in FILE8 depend on the structure, the soil and the wave field, not on the record
([User Guide §2.4](docs/user/USER_GUIDE.md#24-what-runs-in-what-order)).

## Transfer functions and in-structure response spectra

ANALYS does not produce accelerations. It produces **transfer functions**: the complex response of
every degree of freedom to a harmonic control motion of unit amplitude, at each analysis
frequency. MOTION turns them into histories and spectra.

```sassi
* the analysis frequencies: frequency numbers n and f = n x df
LFREQ
* the steady-state motion of every node at the SSI frequency: 24 frames over one period
HARMFRAME,FILE8,3.49,HARM_SSI
```

```action
plot-spectrum: ex01/00085TR_X.TFU, ex01/00085TR_X.TFI, ex01/00041TR_X.TFI
plot-spectrum: ex01/00085TR_X02.RS, ex01/00084TR_X02.RS, ex01/00041TR_X02.RS | log
animate: ex01/HARM_SSI | deformed 0.6 front | SSI system at 3.49 Hz
open-listing: MOTION
explain: NOUT,1,1,1,0,0,1,1,41,82-85
explain: HARMFRAME,FILE8,3.49,HARM_SSI
```

### What this does
`LFREQ` lists the 22 analysis frequencies (0.098 ... 19.995 Hz). The first plot shows the X
transfer function of the roof (node 85): `.TFU` at the 22 computed frequencies and `.TFI`
interpolated by MOTION at all 820 Fourier frequencies from 0 to 20 Hz, with the mat centre
(node 41). The second plot shows the 5 % ISRS of the roof, the floor below it (node 84, 48 ft) and
the mat centre. From the results:

* at 0.098 Hz every transfer function is 1.000: at low frequency the whole system moves with the
  ground (ANALYS checks this; the largest deviation is 0.10 %, at the roof);
* the roof transfer function peaks at **3.47 Hz** with $\lvert H\rvert = \mathbf{13.1}$ (the
  interpolated `.TFI`; the computed frequency next to it, 3.49 Hz, has 13.0): the first mode of the
  soil-structure system;
* the mat centre is not 1: $\lvert H\rvert = 1.42$ at 3.49 Hz and 0.43 at 20 Hz. The foundation does
  not move with the free field;
* above 20 Hz (the last analysis frequency) every transfer function is zero;
* the roof ISRS (5 %) peaks at 8.0 g at 3.55 Hz; the peak accelerations (the zero-period
  accelerations of the ISRS; MOTION listing, "Maximum requested response") grow from 0.36 g at the
  mat to 1.28 g at the roof, for a control motion with a peak of 0.324 g.

`HARMFRAME,FILE8,3.49,HARM_SSI` turns the transfer functions of every node at 3.49 Hz into the
steady-state motion under a harmonic control motion of unit amplitude, $u(t) = \operatorname{Re}\left[H\,e^{i\omega t}\right]$,
sampled in 24 frames over one period (folder `ex01/HARM_SSI`; it prints the frequency it used, the
computed one closest to the request). The **Animate** button stores the frames (`PROCFRAME`) and plays
them on the model (`DEFORMPLOT`, displacements drawn 0.6 ft per unit of control motion, the undeformed
model in grey, seen from the front: `CNGVIEW`). You see the soil-structure mode: the roof swings 13.0 times the ground motion,
94° (about a quarter period) behind it; the mat slides 1.42 times the ground motion, 54° behind it, and
rocks, its edges 32 ft from the centre moving up and down by 1.58. The ground itself, which moves with
amplitude 1 in phase with frame 1, is not drawn.

### Why it matters
The transfer function is where SSI shows itself most clearly: its peak frequency is the SSI
frequency that your ISRS peaks will sit on, and its peak height reflects the system damping. A
design engineer reviewing an SSI analysis looks at transfer functions first, before any spectrum.

### Technical basis
A seismic transfer function is $H(f) = U(f)/U_{\text{cp}}(f)$: the total motion of a DOF per unit
motion of the control point. It is dimensionless and the same for accelerations and displacements.
MOTION interpolates $H$ to every Fourier frequency, multiplies by the Fourier transform of the
control acceleration and transforms back:

```math
a(t) = \operatorname{IFFT}\left[H(f)\,A(f)\right]
```

with $A(f)$ the FFT of the control acceleration (in g), then computes the response spectrum of
$a(t)$ (Nigam-Jennings integration)
([Theory §11](docs/theory/THEORY_MANUAL.md#11-convolution-response-spectra-relative-displacements-and-stresses)).

```figure
tf-to-isrs f0=3.5 beta=0.05
What MOTION does for a floor: the control motion is multiplied by $H(f)$ frequency by frequency,
the floor motion is transformed back, and every oscillator of the spectrum keeps its largest
absolute acceleration. The default mode (3.5 Hz, 5 %) stands for the roof mode of this run; move
$f_0$ and the ISRS peak follows it.
```

The low-frequency check is VP-41, the interpolation VP-28, the convolution VP-29 and the response
spectra VP-30 ([Verification Manual](docs/verification/VERIFICATION_MANUAL.md#vp-28)).

### In ANSYS terms
A harmonic analysis (`ANTYPE,HARMIC`) of a fixed-base model driven by a unit base acceleration
gives the same kind of function: the response of every node per unit input, frequency by
frequency. SASSI keeps it for every DOF in FILE8 and then does the convolution with the actual
record itself; in ANSYS that would be a separate transient analysis.

### Check yourself
Why is every transfer function equal to 1 at 0.1 Hz, whatever the building?

Answer: at very low frequency the wavelengths are miles long (10,000 ft in the sand at 0.1 Hz),
inertia forces vanish and the
soil-structure system moves as a rigid body with the ground. A value far from 1 at the first
frequency points to a modelling error (missing interaction nodes, a disconnected structure, a
mechanism).

## The reference: the same stick on a fixed base

To see what SSI changes, run the same model on a practically rigid site. SASSI drives the
structure through the interaction nodes; with a very stiff soil the mat follows the free field
exactly, which is a fixed base with base excitation.

```sassi
* copy the model, give the copy its own name and directory
CPMODEL,2
ACTM,2
MDL,ex01fb,../ex01_fixed
TIT,Ex01 - the same stick on a practically rigid site (fixed-base reference)
* the same three layers with Vs = 33,000 ft/s (Vp = 2 Vs as before); damping unchanged
L,1,1.6,0.120,66000,33000,0.05,0.05
L,2,3.25,0.125,66000,33000,0.04,0.04
L,3,3.25,0.130,66000,33000,0.02,0.02
* SITE POINT HOUSE ANALYS MOTION only
AOPT,0,0,0,1,1,1,0,0,1,0,1,0,0,0
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
RUNMOTION
* the fixed-base stick at its first mode
HARMFRAME,FILE8,5.0,HARM_FB
```

```action
plot-spectrum: ex01/00085TR_X.TFI, ex01_fixed/00085TR_X.TFI
plot-spectrum: ex01/00085TR_X02.RS, ex01_fixed/00085TR_X02.RS | log
plot-spectrum: ex01/00041TR_X02.RS, ex01_fixed/00041TR_X02.RS | log
animate: ex01_fixed/HARM_FB | deformed 0.6 front | Fixed base at 5.0 Hz
animate: ex01/HARM_SSI | deformed 0.6 front | SSI system at 3.49 Hz
```

### What this does
`CPMODEL,2` copies the active model (everything, options included) to model 2 and `ACTM,2` makes
it active. The copy still has the name and directory of the original (the warning says so), so
`MDL,ex01fb,../ex01_fixed` gives it its own; the path is relative to the current working
directory `ex01/`. `L` redefines the three soil layer types in place. `AOPT` selects the modules
that `AFWRITE` writes decks for. The run takes about 3 s.

The comparison (computed values, 5 % ISRS):

| | Fixed base (rigid site) | SSI |
|---|---|---|
| first-mode peak of the roof transfer function | 13.4 at 4.93 Hz | 13.1 at 3.47 Hz |
| roof ISRS peak, 5 % | 8.03 g at 4.9 Hz | 8.03 g at 3.55 Hz |
| roof ISRS peak, 2 % | 13.4 g at 4.8 Hz | 15.4 g at 3.55 Hz |
| roof peak acceleration (ZPA) | 1.18 g | 1.28 g |
| mat ISRS peak, 5 % | 1.01 g at 2.6 Hz (the free field) | 1.39 g at 3.2 Hz |
| mat peak acceleration (ZPA) | 0.324 g | 0.355 g |

`HARMFRAME,FILE8,5.0,HARM_FB` writes the steady-state motion of the fixed-base model at 5.005 Hz,
the computed frequency closest to its first-mode peak. Play the two animations one after the other
(same scale, one tab each): on the rigid site the mat moves exactly with the ground (1.00, in phase,
no rocking: its edges move vertically by 0.002) and the stick bends in its first mode, the roof at
12.7, 107° behind the ground (5.005 Hz is just above the 4.97 Hz mode, so a little more than a quarter
period). On the soil the same stick sits on a mat that slides (1.42) and rocks (±1.58 at the edges),
and the whole system vibrates at 3.49 Hz instead of 5.0 Hz.

### Why it matters
SSI lowered the frequency of the first mode by 30 %, from 5.0 Hz to 3.5 Hz. Every ISRS peak
moved with it: equipment qualified to a fixed-base ISRS peak near 5 Hz would see its demand peak
somewhere else. The peaks did **not** get lower here (the 5 % roof peak is the same, 8.03 g in both
runs; the 2 % peak and the peak accelerations are higher): the system damping stayed at about 5 %
(Technical basis below; lesson 4 explains why), the input is as strong at 3.5 Hz as at 5 Hz (the
fixed-base mat ISRS, which is the free-field spectrum, has 0.92 g at 3.5 Hz and 0.91 g at 5 Hz), and
the foundation itself now shakes at 3.2-3.5 Hz. SSI is neither automatically beneficial nor
automatically conservative to ignore; it has to be computed, and the uncertainty of the soil
properties (lower, best and upper estimates) moves the SSI frequency further. This is why design
ISRS are enveloped over soil cases and broadened.

### Technical basis
The fixed-base peak at 4.93 Hz (interpolated) is the stick's first mode (4.97 Hz,
[examples/README.md](examples/README.md#example-1-lumped-mass-stick-on-a-rigid-surface-mat)). For a
single mode with hysteretic damping $\beta$ the transfer-function peak is
$1/\bigl(2\beta\sqrt{1-\beta^2}\bigr) = 10.0$ for $\beta = 5\,\%$
([Theory §2](docs/theory/THEORY_MANUAL.md#2-complex-modulus-damping)); for the roof of a
multi-storey stick it is multiplied by the participation of mode 1 at the roof (participation factor
times mode-shape value, $\Gamma_1\phi_1$), here $13.4/10.0 = 1.34$. VP-01 checks the formula through
the whole SASSI chain on a rigid site ([VP-01](docs/verification/VERIFICATION_MANUAL.md#vp-01)). The
half-power width of the SSI peak (in the interpolated `.TFI`) corresponds to about 5 % damping as
well (5.2 %, against 5.0 % on the rigid site): the soil adds damping, but the structure's own 5 %
counts less once the period lengthens, because a smaller share of the system's deformation energy is
in the structure. Lesson 4 takes this apart with the foundation springs and dashpots.

```figure
fixed-base-vs-ssi f=3.5
The stick of this example reduced to its first mode, on a fixed base and on the springs and
dashpots of its mat. The springs move the peak from 4.98 Hz to about 3.9 Hz (SASSI's full model
gives 3.47 Hz because the real mat also bends; lesson 4) while the system damping stays near 5 %: on this layered site little energy
radiates at that frequency, so the wave fronts stay faint. Raise the soil velocities and the peak
moves back towards the fixed-base one.
```

### In ANSYS terms
This run is your usual fixed-base model: base nodes clamped, acceleration applied at the base,
5 % material damping (`MP,DMPR`). Mode 1 at 4.97 Hz is what a modal analysis of the stick gives.

### Try this
Make another copy of model 0 (`ACTM,0`, `CPMODEL,4`, `ACTM,4`, `MDL,ex01stiff,../ex01_stiff`) with
twice the shear-wave velocity in the sand and the gravel (four times the shear modulus), run it
like the fixed-base model (`AOPT`, `AFWRITE`, `RUNSITE` ... `RUNMOTION`) and plot the three roof
transfer functions together. The SSI frequency moves towards the fixed-base value as the soil gets
stiffer (the interpolated roof peak moves from 3.5 Hz to about 4.1 Hz).

```sassi-show
L,1,1.6,0.120,4000,2000,0.05,0.05
L,2,3.25,0.125,6600,3300,0.04,0.04
```

## Kinematic or inertial interaction? A massless foundation

Now remove the masses. With no inertia, any difference between the mat motion and the free-field
motion is kinematic interaction.

```sassi
* copy the original SSI model (model 0) to model 3
ACTM,0
CPMODEL,3
ACTM,3
MDL,ex01ml,../ex01_massless
TIT,Ex01 - massless stick and mat (kinematic interaction only)
* delete the floor masses and give the mat material zero weight
MTDEL,82,85
M,1,5.76E6,0.2,0.0,0.05,0.05,1
* same site and frequencies: reuse the free field and the point-load solutions of the first run
FCOPY,../ex01/FILE1,FILE1
FCOPY,../ex01/FILE3,FILE3
* HOUSE ANALYS MOTION (SITE stays enabled because HOUSE reads the layer table of its deck)
AOPT,0,0,0,1,0,1,0,0,1,0,1,0,0,0
AFWRITE
RUNHOUSE
RUNANALYS
RUNMOTION
```

```action
plot-spectrum: ex01_massless/00041TR_X.TFI, ex01/00041TR_X.TFI, ex01_massless/00085TR_X.TFI
open-listing: MOTION
```

### What this does
`MTDEL` deletes the four floor masses; `M,1,...` redefines the mat material with zero unit weight.
`FCOPY` copies FILE1 (free field) and FILE3 (point-load solutions) from the first run: the site and
the frequencies are unchanged, so SITE and POINT need not run again. HOUSE builds the new mass
matrix and ANALYS solves.

The result: the transfer functions of the mat centre and of the roof are 1.000 at every
frequency, and the vertical motion of the mat edges (the rocking) is zero to round-off. Every node
moves with the free field. The peak acceleration is 0.3235 g at all nodes; the control motion has
0.3239 g, the difference being the content above the 20 Hz cut-off.

### Why it matters
For a surface foundation under vertically propagating waves, kinematic interaction is zero: the free
field moves every point of the ground surface in phase, and a rigid massless mat follows it exactly.
So everything you saw in the SSI run (mat $\lvert H\rvert = 1.42$ at 3.49 Hz, 0.43 at 20 Hz, mat ZPA
0.36 g instead of 0.324 g) is **inertial** interaction: the building and the mat push on the soil
springs. For an **embedded** structure (lesson 5) the massless foundation no longer follows the free
field (the free-field motion varies with depth), and the foundation input motion is reduced at
higher frequencies. Inclined or incoherent waves also create kinematic effects for surface mats;
those are advanced topics.

```figure
kinematic-inertial case=surface
A linear SSI response splits exactly into a kinematic problem (no masses: the foundation input
motion $U_k$) and an inertial one (the masses on the soil springs and dashpots, shaken by $U_k$).
For this surface mat under vertically incident waves $U_k$ is the free-field motion, as the
massless run shows, so all of the change you saw is inertial; switch to the embedded box for the
case of lesson 5.
```

### Technical basis
With zero mass the SSI equation becomes $\left[K^*_s + X_{ff}\right] U = X_{ff}\,U'_f$. With
vertically incident waves $U'_f$ is the same at every node of the ground surface, a rigid-body
translation, and $K^*_s$ times a rigid-body motion is zero, so $U = U'_f$ exactly. VP-15 checks this
case (rigid massless surface foundation under vertical SV and SH: translation 1, rotation 0)
([VP-15](docs/verification/VERIFICATION_MANUAL.md#vp-15)). The "massless foundation" analysis is
also how kinematic interaction is isolated in practice, and the manual recommends it for checking
reduced interaction-node sets against the flexible-volume method (lesson 5).

### Where to go next
Lesson 2 builds the site and the free field (SITE), lesson 3 the soil springs and dashpots (POINT,
ANALYS), lesson 4 steps through this model command by command, and lesson 5 treats embedded
structures, where kinematic interaction matters.

### Check yourself
In the SSI run the mat-centre ISRS has a peak of 1.39 g at 3.2 Hz where the free-field spectrum has
about 1.0 g. Is this kinematic or inertial interaction, and what does it mean for equipment
anchored to the basemat?

Answer: inertial interaction. The massless run shows that the mat would otherwise follow the free
field. The building's response feeds back into the foundation at the SSI frequency, so basemat
equipment sees an amplified spectrum near 3-3.5 Hz; using the free-field (design) spectrum at the
basemat would be unconservative there.
