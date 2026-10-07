---
id: 04-surface-ssi
title: Your first SSI analysis: a stick on a surface mat
part: Fundamentals
order: 4
minutes: 35
example: ex01_surface_stick
summary: Build example 1 command by command, run the SASSI chain, and read the transfer functions, ISRS, relative displacements and member forces, then explain the SSI frequency with the foundation springs and dashpots.
objectives: [Build a stick model on a mat with BEAMS SHELL masses and rigid links, Choose the frequency set (df and NFFT and cut-off), Run SITE POINT HOUSE and ANALYS and check the solution, Produce ISRS with MOTION and check the interpolation, Get relative displacements with RELDISP and member forces with STRESS, Explain the frequency shift and the system damping with the foundation impedance]
prerequisites: [02-free-field, 03-impedance]
---
This lesson builds example 1 command by command: the four-storey stick on a 64 ft × 64 ft mat that
lesson 1 ran in one go. The structure is modelled the way you would model it in ANSYS (beams,
shells, lumped masses, a rigid region); what is new is everything around it: the frequency set,
the soil impedance on the interaction nodes, and the way results are produced from transfer
functions.

The site is the one built in lesson 2 (16 ft of sand, $V_s = 1000\,\text{ft/s}$, on 39 ft of gravel,
$V_s = 1650\,\text{ft/s}$, on weathered rock, $V_s = 3300\,\text{ft/s}$); it is restored automatically
when the lesson opens:

```sassi-setup
MDL,ex01,ex01
TIT,Ex01 - 4-mass stick on a rigid surface mat, vertical SV (X) input
L,1,1.6,0.120,2000,1000,0.05,0.05
L,2,3.25,0.125,3300,1650,0.04,0.04
L,3,3.25,0.130,6600,3300,0.02,0.02
TOPL,1,1,1,1,1,1,1,1,1,1
TOPL,2,2,2,2,2,2,2,2,2,2
TOPL,2,2
WAVE,2,1,1,1,0
```

Units are ft, kip and s (US customary: masses in kip·s²/ft, entered as weights in kips, and
$g = 32.2\,\text{ft/s}^2$); the vertical axis is Z, up; the input is a vertically propagating SV wave
with the control motion in X at the ground surface.

## The basemat

```sassi
GRAVITY,32.2
* mat nodes 1..81: a 9 x 9 grid at grade, 8 ft spacing, numbered row by row (x fastest)
N,1,-32,-32,0
N,9,32,-32,0
FILL,1,9
NGEN,8,9,1,9,1,0,8,0
* M,<nm>,<E>,<nu>,<weight>,<pdamp>,<sdamp>,<type 1>: 10 x concrete stiffness, 0.150 kcf
M,1,5.76E6,0.2,0.150,0.05,0.05,1
GROUP,1,SHELL
GTIT,1,basemat
MACT,1
E,1,1,2,11,10
EGEN,7,1,1
EGEN,7,9,1,8
THICK,1,64,1,5.0
```

```action
plot-model
explain: NGEN,8,9,1,9,1,0,8,0
explain: EGEN,7,9,1,8
```

### What this does
`N` defines nodes 1 and 9 at the two ends of the first row and `FILL` puts nodes 2 to 8 between
them; `NGEN,8,9,1,9,1,0,8,0` copies that row 8 times, adding 9 to the node numbers and 8 ft to Y
each time: 81 nodes. `M,1` is the mat material ($E = 5.76 \times 10^6\,\text{ksf}$, i.e. 40,000 ksi,
ten times concrete, $\nu = 0.2$, unit weight 0.150 kcf, 5 % damping). `GROUP,1,SHELL` opens a shell
group; `E,1,1,2,11,10` is the first 8 ft × 8 ft element, `EGEN` copies it 7 times along X (node
increment 1) and that row 7 times along Y (increment 9): 64 elements. `THICK` gives them 5 ft.

### Why it matters
The mat is deliberately almost rigid, so that the results can be compared with a rigid-foundation
idealisation. In a real model the basemat has its real stiffness: the manual recommends a flexible
basement and, if the rigid assumption must be studied, a New Structure restart with the modulus
multiplied by $10^4$ and $10^5$ ([User Guide §14](docs/user/USER_GUIDE.md#14-engineering-guidance)).

### Technical basis
SHELL is a flat facet element: a plane-stress membrane with incompatible modes plus a
discrete-Kirchhoff plate (DKQ), lumped translational mass, and **no stiffness for the rotation
about its normal** (the drilling rotation, restrained in step 3). The complex modulus
$E^* = E\left(1 - 2\beta^2 + 2i\beta\sqrt{1-\beta^2}\right)$ gives the 5 % hysteretic damping
([Theory §13.1](docs/theory/THEORY_MANUAL.md#13-element-formulations-house); VP-38 checks the
NAFEMS free-vibration benchmarks).

### In ANSYS terms
`N`, `FILL`, `NGEN` and `EGEN` have the same argument order as in APDL, with one trap: SASSI's
`<itim>` counts the copies, APDL's `ITIME` counts the original as well, so SASSI `NGEN,8,9,...` is
APDL `NGEN,9,9,...`. SASSI's `E` gives the element number first (`E,<ne>,<n1>,...`), APDL's `E`
only the nodes. A SHELL group is a set of shell elements with a section thickness, like
`SHELL181` but without its drilling stiffness (step 3).

## The stick and the floor masses

```sassi
* stick nodes 82..85 above the mat centre (node 41), 16 ft storeys
N,82,0,0,16
NGEN,3,1,82,82,1,0,0,16
* beam orientation (K) nodes 86..89 beside each storey in +X
N,86,32,0,8
NGEN,3,1,86,86,1,0,0,16
* stick: concrete, massless (the masses are lumped at the floors)
M,2,5.76E5,0.2,0.0,0.05,0.05,1
* R,<nm>,<A>,<As2>,<As3>,<J>,<I2>,<I3> (ft2, ft4): the shear-wall core
R,1,215.0,108.0,108.0,46000.0,23000.0,23000.0
GROUP,2,BEAMS
GTIT,2,stick
MACT,2
RACT,1
E,1,41,82,86
E,2,82,83,87
EGEN,2,1,2
* 2,200 kips per floor, entered as weights
VAR,FLOOR,82,83,84,85
FOREACH,FLOOR,MT,@FLOOR[#],2200,2200,2200
```

```action
plot-model
explain: R,1,215.0,108.0,108.0,46000.0,23000.0,23000.0
explain: FOREACH,FLOOR,MT,@FLOOR[#],2200,2200,2200
```

### What this does
Four BEAMS elements from the mat centre (node 41) to the roof (node 85), each with a K node that
fixes the orientation of the local axes (local axis 2 points to +X, so bending in X is about local
axis 3). `R` gives the section: $A = 215\,\text{ft}^2$, shear areas 108 ft²,
$J = 46{,}000\,\text{ft}^4$, $I = 23{,}000\,\text{ft}^4$. The material
is massless (`<weight>` = 0) because the floor masses are lumped: `VAR` lists the floor nodes and
`FOREACH` runs `MT,<node>,2200,2200,2200` for each, i.e. 2,200 kips of weight (a mass of
68.3 kip·s²/ft) in X, Y and Z.

### Why it matters
This is the classical lumped-mass stick of nuclear practice. Its fixed-base first mode is at
4.97 Hz and the second at 18.5 Hz ([examples/README.md](examples/README.md#example-1-lumped-mass-stick-on-a-rigid-surface-mat));
lesson 1 reproduced the first on a rigid site. Masses are entered as weights by default
(`MUNITS`); forgetting this, or the gravity, is the classic unit error.

### Technical basis
BEAMS is a 3D Timoshenko beam (shear areas $A_{s2}$, $A_{s3}$) with a consistent mass; the stick
here is massless and the masses are nodal
([Theory §13](docs/theory/THEORY_MANUAL.md#13-element-formulations-house); element closed forms are
checked by VP-37).

## Rigid links, drilling rotations and interaction nodes

```sassi
* K node of the horizontal rigid links (any point off their axes)
N,90,0,0,8
M,3,5.76E8,0.2,0.0,0.0,0.0,1
R,2,64.0,53.0,53.0,576.0,341.0,341.0
GROUP,3,BEAMS
GTIT,3,rigid links
MACT,3
RACT,2
* a rigid spider from the stick base (41) to the 8 surrounding mat nodes
VAR,SPIDER,31,32,33,40,42,49,50,51
FOREACH,SPIDER,E,#,41,@SPIDER[#],90
* restrain the shell drilling rotation where only shells meet
FIXROT
* the 81 mat nodes are interaction nodes
INT,1,81,1,1
INTCOUNT
CALCM
```

```action
plot-model
plot-nodes
explain: FIXROT
```

### What this does
Eight stiff, massless beams ($E = 5.76 \times 10^8\,\text{ksf}$, 1000 times concrete) connect the
stick base to the eight mat nodes around it: a single shell node cannot take the concentrated base
moment of the stick without a large, mesh-dependent local rotation. `FIXROT` fixes the rotation
about the shell normal (ROTZ) at the 72 mat nodes that are connected to shells only; the stick base
and the spider ends keep it, because the beams give it stiffness and because fixing it there would
clamp the torsion of the whole mat to the ground. `INT,1,81,1,1` makes the 81 mat nodes interaction
nodes. `CALCM` checks the mass: 95.4 kip·s²/ft of mat (element weight 3,072 kips) plus
273.3 kip·s²/ft of floors (4 × 2,200 kips), 368.7 kip·s²/ft in total, a weight of 11,872 kips.

### Why it matters
For a surface foundation the interaction nodes are simply the foundation nodes at grade; every
one of them gets the soil impedance. A fixed interaction node would cut the soil out at that node
(`FIXEDINT` finds them), and an unrestrained drilling rotation makes the stiffness matrix singular
(CHECK warning EDU-06). Checking the total mass is the first thing to do with any dynamic model,
SSI or not.

### Technical basis
A flat shell has no stiffness for the rotation about its normal; FIXROT applies `D` to that DOF at
shell-only nodes ([User Guide §5.5](docs/user/USER_GUIDE.md#55-boundary-conditions-and-unstiffened-rotations)).
Interaction nodes carry translations only; the rotations of the mat are transmitted to the soil
by the translations of its nodes
([User Guide §8.2](docs/user/USER_GUIDE.md#82-interaction-nodes)).

### In ANSYS terms
The spider is `CERIG` (or `MPC184` rigid beams); `FIXROT` is the `D,...,ROTZ` you would add for a
shell element type without drilling stiffness.

## The frequency set: df, NFFT and the cut-off

```sassi
* frequency numbers n, f = n x df: 22 frequencies from 0.1 to 20 Hz
FREQ,1,4,20,41,61,82,102,123,143,164,184
FREQ,1,205,225,246,287,328,369,410,492,573,655
FREQ,1,737,819
* SITE: half-space (20 sublayers of layer 3), surface control point, dt 0.005 s, NFFT 8192
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.005,8192,1
LFREQ
```

```action
explain: FREQ,1,4,20,41,61,82,102,123,143,164,184
open-dialog: ANALYSIS/SITE
open-doc: docs/user/USER_GUIDE.md#92-choosing-the-ssi-frequencies
```

### What this does
The control motion has a time step $\Delta t = 0.005\,\text{s}$; with $\text{NFFT} = 8192$ points
the Fourier period is $8192 \times 0.005 = 40.96\,\text{s}$ (the 20 s record plus a 21 s quiet zone
of zeros) and the frequency step is $\Delta f = 1/40.96 = 0.0244\,\text{Hz}$. The analysis
frequencies are integer multiples of $\Delta f$: `FREQ` lists the frequency numbers, from $n = 4$
(0.098 Hz) to $n = 819$ (19.995 Hz): about every 0.5 Hz up to 6 Hz, every 1 Hz up to 10 Hz and every
2 Hz above. `LFREQ` prints them in Hz.

### Why it matters
The frequency set is the main accuracy and cost decision of an SSI analysis:

* the **last frequency is the cut-off**: transfer functions are zero above it, so it acts as a
  low-pass filter on every result. It must cover the input's frequency content and the important
  modes; the manual quotes 30-40 Hz for soil sites and 60-70 Hz for rock sites. The 20 Hz of this
  tutorial keeps it fast; it removes almost nothing from the peak of this record (0.3239 g becomes
  0.3235 g, lesson 1), but it sits just above the stick's second fixed-base mode (18.5 Hz), so a
  production analysis would go higher;
* points must be **dense where the response changes fast**: around the SSI frequencies, which are
  often well below the fixed-base ones (here 3.5 Hz against 4.97 Hz). The manual recommends 40-80
  frequencies for stick models and 100-300 for large FE models; this tutorial uses 22 to stay fast;
* the **quiet zone** must be long enough for the free vibration to decay, or it wraps around to the
  start of the record (EDU-19).

### Technical basis
$\Delta f = 1/(\Delta t \cdot \text{NFFT})$, $f = n\,\Delta f$; MOTION interpolates the transfer
functions from the analysis frequencies to every Fourier frequency $k\,\Delta f$ up to the cut-off
([User Guide §9](docs/user/USER_GUIDE.md#9-frequencies-and-the-fft),
[Theory §10](docs/theory/THEORY_MANUAL.md#10-transfer-function-interpolation)). The sublayers must
pass the cut-off (lesson 2: 102 and 125 Hz here).

### Check yourself
Your record has 4000 points at $\Delta t = 0.01\,\text{s}$ and you choose $\text{NFFT} = 4096$. What
is wrong?

Answer: the Fourier period is 40.96 s and the record lasts 40 s, leaving a quiet zone of under
1 s; the free vibration of the system would wrap around to the start of the record. Use
$\text{NFFT} = 8192$ ($\Delta f = 0.0122\,\text{Hz}$).

## Run the SSI solution: SITE, POINT, HOUSE, ANALYS

```sassi
* POINT: surface foundation, R0 = 0.9 x 8 ft
POINT,0,0,7.2
* HOUSE,<gravity>,<gelev>,<opmode>,<dim>,<imp>,...: g = 32.2 ft/s2, ground at z = 0, 3D, FV
HOUSE,32.2,0,0,2,0,0,0,0,0
* ANALYS: seismic initiation run, control point at the origin; <impe> = 2 also writes the
* 6 x 6 foundation impedance of the mat (FOUNSTIF ...), used in the last step
ANALYS,0,0,0,0,1,0,0,0,0,0,2
AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0
CHECK
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
```

```action
open-listing: HOUSE
open-listing: ANALYS
explain: ANALYS,0,0,0,0,1,0,0,0,0,0,2
open-dialog: ANALYSIS/ANALYS
```

### What this does
`POINT,0,0,7.2`: point loads at the surface only, central zone
$R_0 = 0.9 \times 8\,\text{ft} = 7.2\,\text{ft}$. `HOUSE` sets gravity 32.2 ft/s² (the listing
reports British units), ground elevation 0, a 3D model and the FV method (`<imp>` = 0; for a
surface foundation all methods coincide). `ANALYS,0,0,0,...`: solution, seismic (`<type>` 0),
initiation run (`<mode>` 0), no restart files, amplitudes in the listing; `<impe>` = 2 adds the
global impedance of lesson 3. `CHECK` and `AFWRITE` check the model and write the four decks;
AFWRITE fixes the 5 K-only nodes (86-90) that carry no element (Warning 4). The four modules take
about 0.7 s:

* HOUSE: 438 equations, 72 DOFs fixed by FIXROT, 81 interaction nodes on interface 1;
* ANALYS: at each frequency the structure is condensed onto the 243 interaction DOFs (a Schur
  complement, the same static condensation as an ANSYS superelement), the dense soil impedance is
  added there, and the condensed system is solved; the **low-frequency check** finds every
  transfer function within about 0.1 % of 1 at 0.098 Hz (the largest deviation, 0.103 %, is at the
  roof).

### Why it matters
This is the SSI analysis proper; MOTION, STRESS and RELDISP afterwards only post-process FILE8. The
low-frequency check is your first acceptance test of any SSI run: a model that does not move rigidly
with the ground at 0.1 Hz has a modelling error.

### Technical basis
ANALYS solves, at each frequency, manual Eq. 2.1 (with no excavated soil here):

```math
\left[K^*_s - \omega^2 M_s + X_{ff}\right] U = X_{ff}\,U'_f, \qquad X_{ff} = F_{ff}^{-1}
```

with $X_{ff}$ and $F_{ff}$ on the 243 interaction DOFs.

The structure DOFs are condensed onto the interaction DOFs (sparse LU), the dense condensed system
is solved (dense LU), and the transfer function of every DOF is written to FILE8
([Theory §9](docs/theory/THEORY_MANUAL.md#9-solution-of-the-ssi-equation-analys); VP-41 low-frequency
check, VP-E1 this example, VP-17 inertial SSI of an SDOF against the exact 3-DOF solution).

### In ANSYS terms
A harmonic analysis of the structure with a frequency-dependent, fully coupled spring-dashpot matrix
at the 81 base nodes, driven by forces $X_{ff}\,U'_f$ (the springs pulled by the moving ground),
repeated for 22 frequencies.

## MOTION: transfer functions, accelerations and ISRS

```sassi
* MOTION: 20 s of output, ISRS 0.1-100 Hz at 301 points, factor 1 (motion in g),
* complex .TFI saved for RELDISP, interpolation option 1
MOTION,0,0,0,20,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1
DAMP,0.02,0.05
THFILE,../data/rg160h_030g.acc
THTIT,RG 1.60 horizontal spectrum-compatible motion (EQUAKE, seed 11975)
* NOUT,<dir>,<TF print>,<save TH>,<plot TH>,<plot RS>,<save RS>,<print max>,<nodes>
NOUT,1,1,1,0,0,1,1,41,82-85
NOUT,3,1,0,0,0,0,1,37,45
AOPT,0,0,0,1,1,1,0,0,1,0,1,0,0,0
AFWRITE
RUNMOTION
```

```action
plot-spectrum: ex01/00085TR_X.TFU, ex01/00085TR_X.TFI
plot-spectrum: ex01/00041TR_X.TFU, ex01/00041TR_X.TFI
plot-spectrum: ex01/00041TR_X02.RS, ex01/00082TR_X02.RS, ex01/00083TR_X02.RS, ex01/00084TR_X02.RS, ex01/00085TR_X02.RS | log
plot-history: ex01/00085TR_X.ACC
open-listing: MOTION
explain: MOTION,0,0,0,20,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1
```

### What this does
`MOTION` (19 arguments; the dialog tab shows them): `<dur>` 20 s of output, ISRS from
`<freq1>` 0.1 to `<freq2>` 100 Hz at `<fstep>` 301 log-spaced frequencies, `<mult>` 1 (the record
is already in g), `<cplx>` 1 saves the complex interpolated transfer functions (RELDISP needs them),
`<interp>` 1 the interpolation scheme. `DAMP` gives the ISRS damping ratios (file suffix 01 = 2 %,
02 = 5 %). `THFILE` is the control motion: an RG 1.60 spectrum-compatible record with a peak of
0.324 g. `NOUT,1,...` requests X output (transfer function, history, spectrum, maximum) at the mat
centre and the four floors, `NOUT,3,...` the vertical transfer function of the two mat edges on the
X axis (their motion is the rocking of the mat).

Results (MOTION listing and files):

* roof transfer function: 1.00 at 0.1 Hz, **13.0 at 3.49 Hz** (the computed frequency closest to
  the peak; the interpolated peak is 13.1 at 3.47 Hz); mat centre 1.42 at 3.49 Hz; the two mat
  edges move vertically by 1.58 at 3.49 Hz, in opposite directions: the mat rocks;
* peak accelerations: 0.36 g (mat), 0.47, 0.65, 0.93 and **1.28 g** (roof);
* 5 % ISRS: roof 8.03 g at 3.55 Hz, mat 1.39 g at 3.2 Hz.

### Why it matters
These are your design ISRS (for one direction, one soil case, before broadening). Before using them,
compare the computed (`.TFU`) and interpolated (`.TFI`) transfer functions: the interpolation is
trusted only where computed points support it. The mat-centre plot shows an interpolated peak of
1.70 at 3.37 Hz between computed values of 1.35 (3.00 Hz) and 1.42 (3.49 Hz): exactly the situation
in which the manual recommends adding an analysis frequency (Try this).

### Technical basis
Between analysis frequencies MOTION fits, window by window, the transfer function of a
two-degree-of- freedom system with hysteretic damping, which reproduces the computed values exactly
and resolves resonance peaks between them
([Theory §10](docs/theory/THEORY_MANUAL.md#10-transfer-function-interpolation)); then
$a(t) = \operatorname{IFFT}\left[H(f)\,A(f)\right]$ and the response spectrum of $a(t)$
([Theory §11](docs/theory/THEORY_MANUAL.md#11-convolution-response-spectra-relative-displacements-and-stresses);
VP-28 interpolation, VP-30 spectra).

### Try this
Add five frequencies around the mat peak, rerun the chain and look at the mat-centre TFU again. The
computed value at 3.394 Hz (number 139) is 1.696, the value the interpolation predicted there, and
the ISRS do not change (8.03 g at the roof): here the interpolation was right. Lesson 10 shows how
`CRITFREQ` finds such frequencies automatically.

```sassi-show
FREQ,1,131,135,139,147,151
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
RUNMOTION
```

## RELDISP: relative displacements

```sassi
* RELD,<RelDisOutput>,<RelDispSAll>,<RelDispNumFiles>: save the complex .TFD
RELD,1,0,0
* reference: the mat centre (its interpolated complex transfer function)
RELFILE,00041TR_X.TFI
* RDND,<node>,<X>,<Y>,<Z>,<XX>,<YY>,<ZZ>: X displacement of the mat centre and the floors
RDND,41,1,0,0,0,0,0
RDND,82,1,0,0,0,0,0
RDND,83,1,0,0,0,0,0
RDND,84,1,0,0,0,0,0
RDND,85,1,0,0,0,0,0
AOPT,0,0,0,1,1,1,0,0,1,0,1,0,1,0
AFWRITE
RUNRELDISP
```

```action
plot-history: ex01/00085TR_X.THD, ex01/00083TR_X.THD
open-listing: RELDISP
```

### What this does
RELDISP computes the X displacement of each requested node relative to the reference given by
`RELFILE`, here the mat centre. The maxima (RELDISP listing, in ft): 0.0155, 0.0349, 0.0559 and
**0.0766 ft** (0.19, 0.42, 0.67 and 0.92 in) at the four floors, all at $t = 1.86\,\text{s}$, and
exactly 0 for the mat centre itself.

### Why it matters
Relative displacements are design quantities of their own: seismic gaps between buildings, the
displacement demand on piping, cable trays and ducts crossing between floors or structures,
storey drifts. Referenced to the mat centre they contain the deformation of the stick **and** the
rigid-body rocking of the foundation (the mat rotation times the height). Referenced to the free
field instead (no `RELFILE`), they also include the sliding of the mat relative to the ground,
which matters for components that run from the building into the ground (buried pipes, duct
banks).

### Technical basis
Displacements are not integrated from accelerations (which drifts); they are obtained from the
complex transfer functions and the control displacement spectrum:

```math
D(f) = \left(H_{\text{node}}(f) - H_{\text{ref}}(f)\right) U_g(f), \qquad
U_g = -\frac{g\,A(f)}{\omega^2}, \qquad
d(t) = \operatorname{IFFT}\left[D(f)\right]
```

with $H_{\text{node}}$ and $H_{\text{ref}}$ the complex transfer functions of the node and of the
reference, $A(f)$ the FFT of the control acceleration (in g), $U_g$ the control displacement
spectrum and $d(t)$ the relative displacement history.

([Theory §11.4](docs/theory/THEORY_MANUAL.md#114-relative-displacements-reldisp); VP-34.)

### Try this
Watch the earthquake response in time. `RELDX,0,1` (*Restart For Frame Generation*) makes RELDISP
write one frame per time step in `ex01/THD/`: the X displacements of the five requested nodes relative
to the mat centre, 4800 frames for the 24 s of output, about 40 MB on disk once `PROCFRAME` has stored
them as well (a second or two each to write and to store). `DEFORMPLOT` then plays the first 6 s,
every fourth frame (0.02 s apart), with the displacements drawn 100 times larger: the stick sways back
and forth, its roof reaching 0.077 ft (0.92 in) at 1.86 s. Only the five frame nodes move; the mat
has no frame data and stays at rest.
Delete the two folders afterwards (Learn > Free Disk Space removes the whole workspace).

```sassi-show
RELDX,0,1
AFWRITE
RUNRELDISP
PROCFRAME,THD,THD_ani,RELDISP X relative to the mat centre,3
DEFORMPLOT,THD_ani,1,1200,4,100
WINDOWSETTINGS,UNDEFORMED,1
RELDX,0,0
```

## STRESS: member forces

```sassi
* STRESS,<opmode>,<iter>,<save>,<itran>,<interopt>: save histories and the force transfer functions
STRESS,0,0,1,1,1
* EOUT,<12 codes FXI..MZJ>,<group>,<elements>: 2 = maximum and history (FYI and MZI), 1 = maximum
EOUT,1,2,1,1,1,2,1,1,1,1,1,1,2,1-4
AOPT,0,0,0,1,1,1,0,0,1,0,1,1,1,0
AFWRITE
RUNSTRESS
* save the complete model as commands (ex01/ex01.pre)
WRITE
```

```action
plot-history: ex01/BEAMS_002_00001_FYI.THS, ex01/BEAMS_002_00004_FYI.THS
open-listing: STRESS
open-file: ex01/ex01.pre
```

### What this does
`EOUT` asks for the end forces of the four stick elements (group 2): the maximum of all twelve
components, plus the history of the shear FYI and the moment MZI at the bottom of each storey. The
STRESS listing gives the storey shears 6,764, 6,137, 4,847 and 2,808 kips (bottom to top) and the
base moment $3.29 \times 10^5\,\text{kip\,ft}$. `WRITE` saves the whole model as a `.pre` command
file that `INP` can read back.

### Why it matters
These are the seismic member forces for one direction and one soil case, the input to the
design of the walls and the basemat. A quick equilibrium check: the top-storey shear (2,808 kips)
equals the roof weight times the roof acceleration in g,
$2{,}200\,\text{kips} \times 1.278 = 2{,}811\,\text{kips}$ (both peak at $t = 1.85\,\text{s}$), because
the stick itself is massless.

### Technical basis
STRESS forms the **force** transfer functions element by element (the recovery operator times the
element displacements, minus the rigid-body part), interpolates those (not the nodal transfer
functions) and convolves them with the control displacement spectrum
([Theory §11.5](docs/theory/THEORY_MANUAL.md#115-stresses-and-forces-stress); VP-35 recovery, VP-S1
the base moment of a cantilever stick).

### In ANSYS terms
The BEAMS end forces correspond to the member end forces of a `BEAM188` element in its element
coordinate system (`SMISC` results). SASSI reports them as forces exerted on the element at I and J,
and sign conventions differ between programs, so compare magnitudes. Here they are histories over
the whole record, so the maximum of each component is a true time-history maximum (not a modal
combination).

## Why 3.5 Hz: the foundation springs and dashpots

The stick has its first fixed-base mode at 4.97 Hz; on the soil its first mode is at 3.47 Hz (the
interpolated roof peak; 3.49 Hz is the computed frequency next to it), and its peak is not damped
much more than on a fixed base (lesson 1). The impedance written by ANALYS
explains both.

```sassi
* the soil moduli, for a hand estimate of the springs
LLIST,1,3
* the motion at 3.49 Hz over one period: total, and relative to the mat centre (node 41)
HARMFRAME,FILE8,3.49,HARM_349
HARMFRAME,FILE8,3.49,HARM_349R,24,41
```

```action
open-file: ex01/FOUNSTIF
open-file: ex01/FOUNDAMP
open-listing: MOTION
animate: ex01/HARM_349 | deformed 0.6 front | Total motion at 3.49 Hz
animate: ex01/HARM_349R | deformed 0.6 front | Relative to the mat centre at 3.49 Hz
```

### What this does
`LLIST` gives $G = 3{,}726.7\,\text{ksf}$ for the sand and 10,568.7 ksf for the gravel. FOUNSTIF and
FOUNDAMP hold the 6 × 6 impedance of the mat about its centre at each frequency. At 0.098 Hz:
$K_x = 1.17 \times 10^6\,\text{kip/ft}$ and $K_{yy} = 1.23 \times 10^9\,\text{kip\,ft/rad}$ (rocking
about Y); at 3.49 Hz $K_x = 1.01 \times 10^6\,\text{kip/ft}$,
$K_{yy} = 1.18 \times 10^9\,\text{kip\,ft/rad}$, with damping ratios 7.6 % in sliding and 4.9 % in
rocking (FOUNDAMP).

The two `HARMFRAME` commands write the steady-state motion at 3.49 Hz in 24 frames over one period,
per unit harmonic control motion: the total motion, and the motion relative to the mat centre
(`<Ref>` = 41 subtracts the X, Y and Z motion of node 41 from every node). The animations show the
three parts of conclusion 1 below. In the total motion the mat slides and rocks under the swinging
stick. Relative to the mat centre the sliding is gone: the mat turns about its centre (edges ±1.58)
and the stick tilts with it while it bends; the roof moves 12.0 relative to the mat centre instead
of 13.0 in total.

### Why it matters
Three conclusions an engineer can take from these numbers:

1. **Rocking governs the soil flexibility.** For a force at the roof ($h = 64\,\text{ft}$) the soil
   adds a sliding flexibility of $1/K_x = 0.0119\,\text{in}$ per 1,000 kips and a rocking flexibility
   of $h^2/K_{yy} = 0.0416\,\text{in}$ per 1,000 kips, 3.5 times more. The transfer functions at
   3.49 Hz (MOTION listing) agree: the mat edges, 32 ft either side of the centre, move vertically
   by ±1.58, so the rocking alone moves the roof (64 ft up) by $2 \times 1.58 = 3.2$ times the
   control motion, against a mat translation of 1.42. The three parts are not in phase (compare the
   phases in the listing), so they add as complex numbers: subtracting the mat translation and the
   rocking from the roof transfer function (13.0) leaves about 8.8 for the bending and shear of the
   stick.
2. **Rocking gets almost no radiation damping here.** At 3.49 Hz the rocking damping ratio (4.9 %)
   is hardly more than the material damping of the soil (4.8 % at 0.1 Hz, where nothing radiates);
   sliding has about 3 % of radiation damping on top of it (7.6 %). On a layered site a foundation
   radiates little energy below the first resonance of the soil column, and that threshold depends
   on the motion. Horizontal motion radiates mainly shear waves: the threshold is the shear
   resonance, about 7.3 Hz here ($V_s/(4H)$, lesson 2). Vertical and rocking motion radiate mainly
   compression waves: the threshold is higher, towards the compression resonance $V_p/(4H)$, which
   is twice the shear one here because $V_p = 2V_s$. FOUNDAMP shows it: the sliding damping climbs
   from about 5 Hz (12 % at 5 Hz, 27 % at 7 Hz), the rocking damping only from about 10 Hz (9.8 %
   at 10 Hz, 20 % at 12 Hz). On a uniform half-space (lesson 3) the same $a_0$ would already give
   some radiation damping in rocking. A rocking-dominated SSI mode on a layered site therefore stays
   lightly damped, which is why the peak barely dropped.
3. **Hand springs need the layering.** Pais and Kausel with the sand modulus (half-width
   $B = 32\,\text{ft}$) give $K_x = 5.52\,GB = 6.6 \times 10^5\,\text{kip/ft}$ and
   $K_{yy} = 6.0\,GB^3 = 7.3 \times 10^8\,\text{kip\,ft/rad}$, 44 % and 40 % below the SASSI static
   values: the gravel at 16 ft depth, only $B/2$ below the mat, stiffens the foundation (the coarse
   interaction-node mesh adds a few percent, lesson 3). Half-space formulas on a layered site are a
   first estimate only.

### Technical basis
For an SDOF of stiffness $k$ at height $h$ on a massless rigid foundation with springs $K_x$ and
$K_\theta$, the period lengthens as (Veletsos and Meek, [R2 F.1](docs/spec/R2_benchmarks.md)):

```math
\left(\frac{\tilde{T}}{T}\right)^2 = 1 + \frac{k}{K_x} + \frac{k\,h^2}{K_\theta}
```

with $T$ the fixed-base period and $\tilde{T}$ the period on the soil springs. The system damping
combines the foundation damping with the structural damping reduced by $(\tilde{T}/T)^n$ ($n = 2$
for hysteretic damping): the longer the period, the smaller the share of the structure's own
damping. The $h^2$ term is why rocking dominates for tall, stiff structures. The formula ignores the
mass of the mat (95.4 kip·s²/ft here, a weight of 3,072 kips) and the sliding-rocking coupling of
a welded mat ($K_{x,yy} = -7.1 \times 10^5\,\text{kip}$ at 3.49 Hz), which SASSI includes; VP-17
checks SASSI against the exact 3-DOF solution of an SDOF on a rigid disk and against Veletsos-Meek
([VP-17](docs/verification/VERIFICATION_MANUAL.md#vp-17)). The layer resonances that set the onset
of radiation damping are checked by VP-13 on a disk on a layer over a rigid base, where below them
no energy radiates at all: the horizontal ones fall at $(2n-1)\,V_s/(4H)$ within 1.3 %; for vertical
motion radiation sets in slightly below $V_p/(4H)$ ([R2 C.7](docs/spec/R2_benchmarks.md),
[VP-13](docs/verification/VERIFICATION_MANUAL.md#vp-13)).

```figure
fixed-base-vs-ssi f=3.49
The estimate computed exactly for one mode: $M_1 = 209.8\,\text{kip\,s}^2/\text{ft}$ (a weight of
6,754 kips) at $h_1 = 49.8\,\text{ft}$ and the mat on the springs and dashpots that FOUNSTIF and
FOUNDAMP give at 3.49 Hz. The peak lands near 3.9 Hz, the value of Try this, against SASSI's
3.47 Hz (Try this explains the difference); the rocking dashpot adds almost nothing to the material
damping, so the system damping stays near 5 %.
```

### Try this
Use the formula with the effective mass and height of mode 1 from a modal analysis of the stick in
ANSYS ($k = \omega_1^2 M_1$) and $K_x$, $K_{yy}$ at 3.5 Hz from FOUNSTIF. The estimate lands above
the SASSI value (about 3.9 Hz against 3.47 Hz). The mat mass and the sway-rocking coupling are not
the reason: the whole stick and mat on the same 6 x 6 springs and dashpots, frequency-dependent and
complex, peaks at 3.83 Hz. The reason is the mat itself. $K_G$ assumes that the 81 interaction nodes
move as one rigid body, but the spider clamps the stick to the nine central mat nodes only, and the
5 ft mat, ten times as stiff as concrete, still bends under the rocking moment, which softens the
rocking spring. Make the mat 1000 times as stiff as concrete (`M,1,5.76E8,...`) and SASSI's roof
peak moves to 3.83 Hz (|H| = 12.7). A foundation that looks rigid in a static check can be flexible
for SSI; this is why SASSI keeps the real mat and all its interaction nodes rather than a rigid-body
spring set. Then repeat lesson 1's Try this ($V_s$ of the sand and the gravel doubled): the frequency
moves back towards 5 Hz (about 4.1 Hz).

### Check yourself
Your fixed-base ISRS of this building peaks at 5 Hz. After the SSI analysis the peak is at 3.5 Hz
with about the same height. A component on the roof has its own frequency at 3.6 Hz. What happened
to its demand, and what else would you check before accepting the SSI result?

Answer: its demand went up about threefold: at 3.6 Hz the 5 % roof ISRS is 7.9 g with SSI against
2.8 g on the low-frequency flank of the fixed-base peak (lesson 1's two runs). Before accepting,
check the soil cases (lower and upper bound shift the peak), the interpolation around the peak
(TFU against TFI), the cut-off and the quiet zone, and envelope and broaden the spectra over the
cases.
