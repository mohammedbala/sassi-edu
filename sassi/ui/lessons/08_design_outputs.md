---
id: 08-design-outputs
title: From SSI results to design: forces, displacements, ANSYS
part: Design applications
order: 8
minutes: 35
example: ex01_surface_stick
summary: Extract the quantities a structural design needs from an SSI run (member forces, basemat section resultants, relative displacements) and hand the SSI motion to a refined ANSYS model with the Option A two-step approach (LOADGEN).
objectives: [Read member-force histories and check them by equilibrium, Choose the right reference for relative displacements, Get section resultants of a shell structure with a section cut, Write equivalent static and dynamic ANSYS loads with LOADGEN and know their assumptions]
prerequisites: [04-surface-ssi]
---
ISRS are only one of the deliverables of a seismic analysis. The structural design itself needs
**member forces** (storey shears and moments, wall and slab forces), **section resultants** of walls
and mats, and **relative displacements** (drifts, gaps, seismic anchor motions of piping that spans
between floors or buildings). In a fixed-base ANSYS design you get these from a response-spectrum or a
transient analysis of one model. With SSI there are two routes:

1. directly from the SSI model, with **STRESS** (element forces) and **RELDISP** (relative
   displacements), both computed from the transfer functions and the record; or
2. with the **two-step approach** (ACS SASSI Option A): the SSI model gives the motion, and a refined
   ANSYS model of the structure, loaded by that motion, gives the forces. Module **LOADGEN** writes the
   ANSYS loads.

This lesson does both on example 1, the four-storey stick on a 64 ft × 64 ft surface mat of lesson 04.
Opening the lesson runs the complete example (about 4 s): the setup below reads
`ex01_surface_stick.pre`, which builds the model and runs SITE, POINT, HOUSE, ANALYS, MOTION, STRESS
and RELDISP. The steps then run in about 3 s.

```sassi-setup
INP,ex01_surface_stick.pre
```

## Member forces of the stick: storey shears and moments

The setup's STRESS run computed the end forces of the four stick elements. Here you load two of the
histories it saved and check them against the inertia of the roof, the kind of equilibrium check
you would do on any ANSYS result.

```sassi
* the shear at the bottom of the top storey (element 4, FYI) and the roof acceleration (g)
READTH,BEAMS_002_00004_FYI.THS,0,1
READTH,00085TR_X.ACC,0,2
* roof mass x absolute roof acceleration = (2200 kips / g) x a(g)
LINECOMBIN,3,2,2200
LINENAME,3,roof mass x roof acceleration (kips)
WRITETH,roof_inertia.th,3
```

### What this does
* The example requested the forces with `EOUT,1,2,1,1,1,2,1,1,1,1,1,1,2,1-4`: twelve codes for the
  end forces FXI FYI FZI MXI MYI MZI FXJ ... MZJ (0 none, 1 maximum, 2 maximum and history), then the
  group (2, the stick) and the elements (1-4). FYI (shear along local axis 2, which points to the
  K node, i.e. +X) and MZI (moment about local axis 3 at the bottom end I) are saved as histories
  `BEAMS_002_0000k_FYI.THS` and `..._MZI.THS`.
* `READTH,<file>,0,<line>` loads a history (format 0: time step first, then values),
  `LINECOMBIN,3,2,2200` scales the roof acceleration (in g) by the roof weight (2,200 kips), and
  `WRITETH` writes the result as a history.

### Why it matters
The STRESS listing gives the maxima of the storey shear FYI and moment MZI (bottom of each storey,
kips and kip·ft):

| Storey | Max shear FYI (kips) | Max moment MZI (kip·ft) |
|---|---|---|
| 1 (base, element 1) | 6,764 | 328,742 |
| 2 | 6,137 | 220,554 |
| 3 | 4,847 | 122,489 |
| 4 (top, element 4) | 2,808 | 44,935 |

All maxima occur at $t \approx 1.85\,\text{s}$. The two histories of the first plot lie on top of each
other: the stick is massless, so the top storey carries exactly the roof mass times its absolute
acceleration (peak $2{,}200\,\text{kips} \times 1.278\,\text{g} = 2{,}811\,\text{kips}$ against
$\mathrm{FYI} = 2{,}808\,\text{kips}$, 0.09 % apart: MOTION interpolates the nodal transfer functions,
STRESS the force transfer functions). FYI is the force exerted on the element at its lower end, so it
has the sign of $m\,a$. The shear grows down the stick to 6,764 kips at the base, 0.77 times the weight
of the four floors (8,800 kips). These forces contain the SSI effects directly: the foundation
flexibility, the rocking and the radiation damping of the soil are in them.

### Technical basis
HOUSE stores per element a complex recovery operator $S$ (for a beam the local end forces $k_L T u$
exerted on the element); STRESS forms the **stress transfer functions** $S\,(U_e - r_e)$ at the SSI
frequencies ($r_e$ is the rigid-body part), interpolates them, and convolves them with the control
displacement spectrum. Interpolating the force transfer functions, not the nodal ones, keeps forces
accurate between SSI frequencies.
[Theory §11.5](docs/theory/THEORY_MANUAL.md#115-stresses-and-forces-stress); verified by VP-35 (rigid
body gives zero force, beam equilibrium) and [VP-S1](docs/verification/VERIFICATION_MANUAL.md#vp-s1)
(base moment of a cantilever stick against an independent calculation).

### In ANSYS terms
The member forces of a transient analysis (ETABLE / element results in POST26). The difference: they
are computed per frequency and convolved, so a maximum is a true time-domain maximum, not a modal
combination.

```action
open-listing: STRESS
plot-history: ex01/BEAMS_002_00004_FYI.THS, ex01/roof_inertia.th
plot-history: ex01/BEAMS_002_00001_FYI.THS, ex01/BEAMS_002_00004_FYI.THS
plot-history: ex01/BEAMS_002_00001_MZI.THS
explain: EOUT,1,2,1,1,1,2,1,1,1,1,1,1,2,1-4
```

## Relative displacements: to the basemat or to the free field?

The setup's RELDISP run gave the displacements of the floors relative to the mat centre
(`RELFILE,00041TR_X.TFI`). This step keeps them under new names and runs RELDISP again with the
free field as reference.

```sassi
* keep the displacements relative to the mat centre as MAT_*
VAR,THD,00041TR_X.THD,00082TR_X.THD,00083TR_X.THD,00084TR_X.THD,00085TR_X.THD
FOREACH,THD,FCOPY,@THD[#],MAT_@THD[#]
* no reference file: displacements relative to the free-field ground motion
RELFILE
AOPT,0,0,0,0,0,0,0,0,0,0,0,0,1,0
AFWRITE
RUNRELDISP
```

### What this does
* `RELFILE` without a file name clears the reference: RELDISP then subtracts the free-field motion
  (unit amplitude, zero phase, in the input direction). The node requests of the example
  (`RDND,41,1,...` and `RDND,82` ... `85`) stay; `AOPT` enables RELDISP only.
* The `.THD` files of nodes 41 (mat centre) and 82-85 (floors) are overwritten with the new
  reference; the `MAT_` copies keep the old one.

### Why it matters
The roof moves 0.0766 ft (0.92 in) relative to the mat centre (floors 0.0155, 0.0349, 0.0559 and
0.0766 ft, all at $t \approx 1.86\,\text{s}$) but 0.0835 ft (1.00 in) relative to the free field; the
mat centre itself moves up to 0.0074 ft (0.09 in) relative to the free field (at 6.2 s): the foundation
translates and rocks on the soil. A
fixed-base model would report only the deformation of the stick.

Which one to use:

* **interstorey drift** for the structure: the difference of two floors (either reference gives the
  same difference);
* **gap to an adjacent building** or a **seismic anchor motion** of piping between buildings: the
  displacement relative to the free field (or, better, the difference of the two buildings' motions
  in one SSI model);
* displacements relative to the basemat include the **rigid-body rocking** of the foundation times the
  height, which moves equipment but does not strain the structure.

### Technical basis
```math
D(f) = \left[H_\text{node}(f) - H_\text{ref}(f)\right] U_g(f), \qquad U_g(f) = -\frac{g\,A(f)}{\omega^2}, \qquad d(t) = \mathrm{IFFT}[D]
```
with $H_\text{node}$ and $H_\text{ref}$ the transfer functions of the node and of the reference, $A(f)$
the FFT of the control acceleration (in g) and $U_g(f)$ the control displacement; computed from the
complex `.TFI` transfer functions (MOTION `<cplx>` = 1), not by double integration
of accelerations, which drifts. The FFT solution is periodic: its zero-frequency term is undetermined,
so the histories are relative displacements with zero mean over the Fourier period.
[Theory §11.4](docs/theory/THEORY_MANUAL.md#114-relative-displacements-reldisp); verified by
[VP-34](docs/verification/VERIFICATION_MANUAL.md#vp-34).

```action
plot-history: ex01/MAT_00085TR_X.THD, ex01/00085TR_X.THD, ex01/00041TR_X.THD
open-listing: RELDISP
explain: RDND,85,1,0,0,0,0,0
```

## Basemat forces and a section cut

Walls and mats are shells or solids: the design needs **resultants** over a section (shear, normal
force, moment), not element stresses. STRESS writes the element-centre stresses of every time step,
and a section cut integrates them over a plane. Here: the basemat across its full width at
$x = 12\,\text{ft}$, between the stick's base spider and the mat edge.

```sassi
* SHELL forces of the 64 basemat elements (FXX FYY FXY MXX MYY MXY), maxima only
EOUT,0
EOUT,1,1,1,1,1,1,0,0,0,0,0,0,1,1-64
* no transfer-function files for these 384 components (STRESS <itran> = 0)
STRESS,0,0,1,0,1
* element-centre stress frames (NSTRESS/ESTRESS_nnnnn.ess), every 5th time step
SECDATAOPT,1
STRESSX,0,0,5
AOPT,0,0,0,0,0,0,0,0,0,0,0,1,0,0
AFWRITE
RUNSTRESS
* cut 1 = the 8 elements between x = 8 and 16 ft; resultants on the plane x = 12 ft
CUTVOL,1,7.9,16.1,-32.1,32.1,-0.1,0.1
CALCSECTHIST,NSTRESS/ESTRESS.lst,1,12,0,0,1,0,0,0,1,0,0,0.025,mat_cut_x12.csv
* back to the example's STRESS settings for later runs
STRESS,0,0,1,1,1
SECDATAOPT,0
STRESSX,0,0,1
```

### What this does
* `EOUT,0` clears the beam request; the new request asks for the six shell components (membrane
  forces per unit length FXX FYY FXY, moments per unit length MXX MYY MXY, local axes) of group 1.
* `STRESS,0,0,1,0,1` switches off `<itran>` (Output Transfer Function): with it on, STRESS would also
  write a `.TFU` and a `.TFI` file for each of the $64 \times 6 = 384$ requested components (about 20
  MB) that this step does not use. The last block restores the example's `STRESS,0,0,1,1,1`.
* `SECDATAOPT,1` makes STRESS write a frame of element-centre stresses per output step;
  `STRESSX,0,0,5` (`<skip>` = 5) keeps every fifth step ($\Delta t = 0.025\,\text{s}$, 960 frames, 7.5
  MB) to keep the files small.
* `CUTVOL,1,...` puts into cut 1 the elements whose nodes all lie in the box ($x$ from 8 to 16 ft, the
  full width in $y$). `CALCSECTHIST,<list>,<cut>,<px,py,pz>,<nx,ny,nz>,<rx,ry,rz>,<sysno>,<ts>,<csv>`
  integrates them on the plane through (12, 0, 0) with normal X, local $x$ axis along global Y, for
  every frame of the list, and writes a CSV of time, $F_x$, $F_y$, $F_z$ (normal force), $M_x$, $M_y$,
  $M_z$ (torsion)
  with a final row of signed absolute maxima. `<ts>` = 0.025 s is the time between frames.

### Why it matters
The CSV ends with the maxima: a bending moment $M_x =$ **82,740 kip·ft** across the 64 ft width at
$t = 1.85\,\text{s}$ (the time of the peak base moment), about 1,290 kip·ft per foot on average, and a
normal (membrane) force $F_z = 2{,}592\,\text{kips}$ at the same time, with an almost equal peak
(2,591 kips) at $t = 6.2\,\text{s}$, the time when the mat slides most relative to the free field. The
section moment is the sum of the shell moments MXX of the eight cut elements times their 8 ft width
(480 to 2,230 kip·ft/ft in the STRESS listing, largest in the middle strips under the stick). The
transverse shear $F_y$ is zero: the Kirchhoff SHELL element reports no transverse shear forces; model
the mat with TSHELL (Mindlin, QXZ and QYZ output) when you need them.

The same commands give the base shear and overturning moment of a shear wall, the forces across a
construction joint or the moment in a slab strip: the resultants you size reinforcement for.

### Technical basis
Each cut element contributes the integral of its stresses over its intersection with the plane;
resultants are taken about the area centroid of the section in the local axes ($x$ along the "right"
vector projected on the plane, $z$ along the normal). The frame step `<skip>` samples the history, so
the maxima of a skipped history can miss the true peak by a little: use `<skip>` 1 for final design
values. [User Guide §11.6](docs/user/USER_GUIDE.md#116-section-cuts); verified by
[VP-49](docs/verification/VERIFICATION_MANUAL.md#vp-49) (a solid and a shell wall).

### In ANSYS terms
`FSUM` over the elements on one side of the cut and the nodes on the cut (`ESEL`/`NSEL` first),
repeated for every result set in a `*DO` loop over `SET`; in Workbench, a construction surface with
force and moment reaction probes.

```action
open-file: ex01/mat_cut_x12.csv
open-listing: STRESS
explain: CALCSECTHIST,NSTRESS/ESTRESS.lst,1,12,0,0,1,0,0,0,1,0,0,0.025,mat_cut_x12.csv
```

## Option A: equivalent static loads for ANSYS (LOADGEN)

The two-step approach. Split the motion of the structure into a rigid-body motion and a relative one,
$u = u_r + \iota\,u_g$: a rigid translation strains nothing, so the structure with its foundation
interface prescribed obeys

```math
M_{ss}\,\ddot{u}_{r,s} + C\,\dot{u}_{r,s} + K_{ss}\,u_{r,s} = -M_{ss}\,\iota_s\,a_g(t) - K_{sb}\,u_{r,b}(t)
```

with $s$ the structural and $b$ the prescribed interface DOFs, $\iota\,u_g$ the rigid-body translation
with the ground motion $u_g$, $u_r$ the motion relative to it and $a_g$ the ground acceleration. The
**equivalent static** version freezes this at a critical time $t^*$ and drops damping:
$K_{ss}\,u_s = -M_{ss}\,a_s(t^*) - K_{sb}\,u_b(t^*)$, i.e. the inertia forces of the absolute
accelerations plus the interface displacements, applied statically to an ANSYS model of the structure.

```sassi
* LOADGEN,<data>,<multi>,<rotdisp>,<rotacc>,<masstype>,<genmass>,<source>
LOADGEN,DISPACC,0,0,0,LUMPED,1,FILE8
* critical time: the largest base shear along the input direction
LGTIME,V,1
RUNLOADGEN,STATIC
```

### What this does
* `DISPACC` ("Disp. and Accel."): nodal inertia forces $F = -m\,a(t^*)$ at every mass node, plus the
  displacements of the interface nodes (by default the 81 interaction nodes) relative to the free
  field at $t^*$. `ACC` would fix the interface instead (the classical fixed-base equivalent static
  analysis).
* `LUMPED,1`: generate the lumped masses (`ex01.masl`) from the HOUSE mass matrix (row sums per
  direction). `FILE8`: compute the accelerations and displacements from FILE8 and the record with the
  MOTION and RELDISP algorithms, so no `.ACC` file is needed for the 81 mat nodes.
* `LGTIME,V,1`: one critical time, the largest base shear (sum of the inertia forces) along the input
  direction. `RUNLOADGEN,STATIC` writes the deck `ex01.lgn`, runs LOADGEN and writes `ex01_LGS.inp`.

### Why it matters
The LOADGEN listing reports: masses of 368.7 kip·s²/ft per direction, i.e. 11,872 kips of weight
($4 \times 2{,}200\,\text{kips}$ of floors plus 3,072 kips of mat; the mat nodes are interface nodes,
so their inertia goes straight to the reactions), the critical time **$t = 1.855\,\text{s}$**, and a
base shear of **6,772 kips** (sum of the 12 inertia forces of the floors). The STRESS base shear of the
stick at that time is 6,764 kips: the two routes agree to 0.12 %. The moment of the inertia forces
about the base is 329,113 kip·ft (STRESS: 328,742 kip·ft, 0.11 % apart). `ex01_LGS.inp` holds the 12
`F` commands and 243 `D` commands (81 interface nodes × 3 translations).

The listing also warns that the interface rotations are not prescribed (`<rotdisp>` 0): with 81
interface nodes the translations restrain the ANSYS model; set `<rotdisp>` 1 for a single interface
node (a stick on one support).

In ANSYS: `/INPUT,ex01,inp` (the structure, written in the next step), then `/INPUT,ex01_LGS,inp`
(loads, then SOLVE); in POST1 `PRRSOL,F` must show a total X reaction equal to the sum of the applied
inertia forces with the opposite sign (6,772 kips): by equilibrium the reactions balance the applied
forces whatever the ANSYS mesh. Dead load is not included: combine with a gravity load case as your
design code requires.

### Technical basis
Masses: $m_{i,d} = \sum_j M_{(i,d),(j,d)}$ (row sums of the HOUSE mass matrix per direction, exact for
lumped masses). Critical times: the largest local maxima of $\lvert V(t)\rvert$. The relative
displacements are started at rest (the zero-mean offset of the periodic FFT solution is removed).
[Theory §20](docs/theory/THEORY_MANUAL.md#20-option-a-ssi-loads-for-a-second-step-ansys-analysis),
[OPTION_A.md §5](docs/user/OPTION_A.md#5-equivalent-static-loads-ansys-eq-static-load); verified by
[VP-LA1](docs/verification/VERIFICATION_MANUAL.md#vp-la1): the APDL inertia forces of a stick sum to
the STRESS base shear and moment within 0.04 %.

### Check yourself
LOADGEN wrote the loads at one instant. Is the ANSYS stress state at that instant the maximum of every
member force?

Answer: No. It is the state at the time of the largest base shear; the moment at mid-height, a slab
shear or a wall stress can peak at other times. Ask for more critical times (`LGTIME,V,3`, or
`LGTIME,MX|MY,...` for overturning, `LGTIME,ACC,...` for a node) with `<multi>` = 1 to get one load
file per time, and envelope the ANSYS results; or use the dynamic route of the next step.

```action
open-file: ex01/ex01_LGS.inp
open-listing: LOADGEN
open-file: ex01/ex01.masl
explain: LOADGEN,DISPACC,0,0,0,LUMPED,1,FILE8
```

## Option A: dynamic loads and the structural model in APDL

The **dynamic** second step keeps the whole equation: an ANSYS transient analysis with `ACEL` = the
ground acceleration and the interface displacements relative to the ground prescribed as tables.
The `ANSYS` command writes the structural model itself.

```sassi
* LOADGENDYN,<alpha>,<beta>,<method>,<refnode>,<source>,<rotdisp>,<rotacc>,<zeta>,<f1>,<f2>
LOADGENDYN,,,REL,0,FILE8,0,0,0.05,3.5,20
* LGOPT,<digits>: 6 significant digits in the APDL tables instead of 12 (a smaller file)
LGOPT,6
RUNLOADGEN,DYNAMIC
* the structure as APDL with BEAM188 / SHELL181 (ex01.inp)
EDUOPT,ANSYSMODERN,1
ANSYS
```

### What this does
* Method `REL` (the manual's): `ACEL` = the control motion, `D` = TABLE arrays of the interface
  displacements relative to it. `<refnode>` 0 = the free field. Rayleigh damping from $\zeta = 5\,\%$ at
  `<f1>` = 3.5 Hz (the SSI frequency of the stick, lesson 04) and `<f2>` = 20 Hz.
* `RUNLOADGEN,DYNAMIC` writes `ex01_LGD.inp`: 244 TABLE arrays of 4800 values (the ACEL history and
  243 interface displacement histories), the transient solution commands and `SOLVE`. The file
  grows with the interface DOFs times the samples; `LGOPT,6` writes 6 significant digits instead of
  the default 12, which is ample for loads (a relative precision of 1e-6) and keeps the file at about
  17 MB instead of 24 MB. Keep 12 digits when you want to replay the SSI solution exactly (VP-LA2).
  Read the LOADGEN listing rather than the file.
* `EDUOPT,ANSYSMODERN,1` + `ANSYS` write `ex01.inp`: nodes, BEAM188 and SHELL181 elements, materials
  with `DENS` = $\text{weight}/g$, MASS21 masses, and the interaction nodes as the component `SSI_INT`.

### Why it matters
The listing prints the Rayleigh damping that ANSYS will apply: 15.1 % at 1 Hz, 7.9 % at 2 Hz,
4.0 % at 5 Hz, 3.6 % at 10 Hz and 5.0 % at 20 Hz, against the 5 % of the SSI model at every
frequency. For this linear stick that damping model is what separates the ANSYS transient from the
SSI solution; to a lesser degree so do the Newmark time integration, the element formulations
(SHELL181 is a Mindlin shell, SASSI's SHELL a Kirchhoff one) and the mass matrix (SASSI lumps the
shell mass; ANSYS uses consistent mass unless you set `LUMPM,ON`).

**What Option A assumes** (and when to redo the SSI): the ANSYS model must not change the motion of the
foundation-soil interface. It is exact for a linear structure with the same mass, stiffness and
damping as the SSI model (VP-LA2 replays the exported second step and recovers the SSI accelerations
to 5e-9). It becomes approximate when ANSYS uses Rayleigh damping instead of SASSI's
frequency-independent damping (exact at two frequencies only, "a significant limitation" in the
manual's words), when the static step represents one instant, and when the refined model is much
softer, heavier or more nonlinear than the SSI model: then the interface motion it would produce is
different, and the SSI analysis must be repeated with the refined properties.

### Technical basis
Rayleigh damping:

```math
\begin{aligned}
C &= \alpha M + \beta K, \qquad \alpha = \frac{2\zeta\,\omega_1\omega_2}{\omega_1 + \omega_2}, \qquad \beta = \frac{2\zeta}{\omega_1 + \omega_2}\\
\zeta(\omega) &= \frac{\alpha}{2\omega} + \frac{\beta\,\omega}{2}
\end{aligned}
```

with $\omega_1 = 2\pi f_1$ and $\omega_2 = 2\pi f_2$: $\zeta(\omega)$ is equal to $\zeta$ at $f_1$ and
$f_2$, lower in between, higher outside. Radiation damping is not in $C$: it is in the prescribed
interface motion, which is why the D tables matter.
[OPTION_A.md §6](docs/user/OPTION_A.md#6-dynamic-loads-ansys-dynamic-load),
[ANSYS.md §3](docs/user/ANSYS.md#3-what-the-ansys-command-writes); verified by
[VP-LA2](docs/verification/VERIFICATION_MANUAL.md#vp-la2). The APDL has been checked by a parser in
SASSI-EDU, not run in ANSYS.

### In ANSYS terms
`ANTYPE,TRANS` with `TRNOPT,FULL`, a constant `DELTIM` equal to the SSI time step, `ALPHAD`/`BETAD`,
`ACEL,%LG_ACX%,0,0` and `D,node,UX,%LGD_node_UX%`. ANSYS displacements are then relative to the
moving ground; stresses and internal forces are the actual ones; absolute accelerations = ANSYS
acceleration + ACEL.

### Try this
Change `<f1>` to 1 Hz (`LOADGENDYN,,,REL,0,FILE8,0,0,0.05,1,20`) and rerun: the listing's damping
at 5 Hz drops. Which choice of $f_1$ and $f_2$ brings the Rayleigh damping closest to 5 % over the band
that carries the response of this stick?

```action
open-listing: LOADGEN
open-file: ex01/ex01.inp
open-doc: docs/user/OPTION_A.md#2-why-it-works-the-theory-in-one-page
```
