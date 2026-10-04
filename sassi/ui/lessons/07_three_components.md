---
id: 07-three-components
title: Three earthquake directions and design ISRS
part: Design applications
order: 7
minutes: 35
example: ex05_xyz_simultaneous
summary: Apply X, Y and Z input as SV, SH and P waves in one ANALYS run, see how an eccentric mass couples the directions, combine the co-directional responses (SRSS, 100-40-40) and turn them into a broadened design ISRS.
objectives: [Set up the three vertically propagating wave fields of a 3D seismic input, Solve the three directions in one ANALYS run (simultaneous cases), Read direct and cross-coupled responses of an eccentric structure, Combine co-directional ISRS by SRSS and 100-40-40, Envelope and broaden an ISRS into the design deliverable]
prerequisites: [04-surface-ssi, 06-seismic-input]
---
A fixed-base design applies three earthquake components, two horizontal and one vertical, and for
every response quantity (an ISRS at a floor in one direction, a member force) combines the
contributions of the three components: SRSS or the 100-40-40 rule. Then the ISRS of the soil
cases are enveloped and their peaks broadened. The design ISRS you hand to the equipment qualification
engineer is that broadened envelope.

SSI changes how the input is defined, not this logic. In SASSI each input direction is a **wave
field**: a vertically propagating SV wave for X, an SH wave for Y and a P wave for Z, each normalised
to a unit motion at the control point. SITE builds the three free fields, ANALYS solves the three
cases with one factorisation per frequency, and MOTION convolves each case with its own record. The
directional combination is post-processing, done with the line operations of the Plot menu.

Example 5 is built to show coupling: a two-storey stick on a 12 m × 12 m basemat with its 300 t roof
mass 2.24 m off the stick axis, so an X input also moves the roof in Y and twists the stick. The lesson
runs in about 6 s.

## The site and the frequencies

The site of lesson 04 (sand over gravel over rock) and a set of 19 SSI frequencies from 0.1 to 15 Hz.

```sassi
MDL,ex05,ex05
TIT,Ex05 - eccentric two-storey stick, X + Y + Z input with ANALYS simul = 1
L,1,0.5,19.0,600,300,0.05,0.05
L,2,1.0,20.0,1000,500,0.04,0.04
L,3,1.0,21.0,2000,1000,0.02,0.02
TOPL,1,1,1,1,1,1,1,1,1,1
TOPL,2,2,2,2,2,2,2,2,2,2
TOPL,2,2
* df = 1/(8192 x 0.005 s) = 0.0244 Hz; numbers 4 ... 614 = 0.1 ... 15 Hz
FREQ,1,4,20,41,61,82,102,123,143,164,184
FREQ,1,205,246,287,328,369,410,492,573,614
```

### What this does
Ten sublayers of 0.5 m of sand ($V_s = 300\,\text{m/s}$), twelve of 1 m of gravel
($V_s = 500\,\text{m/s}$), rock half-space ($V_s = 1000\,\text{m/s}$); `FREQ` appends 19 frequency
numbers to set 1. The damping of the layers (5 %, 4 %, 2 %) stands for strain-compatible values: in a
project these layers come from SOIL (lesson 06).

### Why it matters
The frequency set must resolve the peaks of every case: the X and Y responses of this structure peak
near 5 Hz, the vertical one near 10 Hz, so the set extends to 15 Hz. One set serves the three
directions, because ANALYS solves them together.

### Technical basis
$f = n\,\Delta f$, $\Delta f = 1/(\Delta t\cdot\mathrm{NFFT})$
([User Guide §9](docs/user/USER_GUIDE.md#9-frequencies-and-the-fft)); the interpolation between them
is lesson 10's subject.

```action
plot-layers
```

## The structure: a stick with an eccentric roof mass

A 12 m × 12 m basemat of 4 × 4 shells (3 m mesh, 1.5 m thick, ten times stiffer than concrete)
carries a two-storey stick (4 m storeys). The roof mass, node 30 at (2, 1, 8), is tied to the stick
top (node 27) by a rigid arm.

```sassi
* basemat nodes 1..25 (3 m grid), centre node 13; stick nodes 26 (z = 4) and 27 (z = 8)
N,1,-6,-6,0
N,5,6,-6,0
FILL,1,5
NGEN,4,5,1,5,1,0,3,0
N,26,0,0,4
N,27,0,0,8
N,28,6,0,2
N,29,6,0,6
* the eccentric roof mass and the K node of the rigid links
N,30,2,1,8
N,31,0,0,2
M,1,3.0E8,0.2,24.0,0.05,0.05,1
M,2,3.0E7,0.2,0.0,0.05,0.05,1
M,3,3.0E10,0.2,0.0,0.0,0.0,1
R,1,10.0,5.0,5.0,10.0,5.0,5.0
R,2,9.0,7.5,7.5,11.4,6.75,6.75
GROUP,1,SHELL
MACT,1
E,1,1,2,7,6
EGEN,3,1,1
EGEN,3,5,1,4
THICK,1,16,1,1.5
GROUP,2,BEAMS
MACT,2
RACT,1
E,1,13,26,28
E,2,26,27,29
* rigid spider at the stick base and the rigid arm to the roof mass
GROUP,3,BEAMS
MACT,3
RACT,2
VAR,SPIDER,7,8,9,12,14,17,18,19
FOREACH,SPIDER,E,#,13,@SPIDER[#],31
E,9,27,30,31
FIXROT
* 300 t per storey (weights 2943 kN) and the rotary inertia of each 12 m x 12 m slab (I x g)
MT,26,2943,2943,2943
MT,30,2943,2943,2943
MR,26,35316,35316,70632
MR,30,35316,35316,70632
INT,1,25,1,1
POINT,0,0,2.7
HOUSE,9.81,0,0,2,0,0,0,0,0
```

### What this does
* `MT,30,2943,2943,2943` puts 300 t (as a weight, kN) at node 30 in X, Y and Z; `MR,30,...` adds the
  rotary inertia of the slab about X, Y and Z ($I \times g$: 3600 t m² about X and Y, 7200 t m² about
  Z).
* `E,9,27,30,31` is the rigid arm (material 3, 1000 × concrete stiffness, massless) from the stick
  top to the mass: the mass sits 2 m off the axis in x and 1 m in y.
* `INT,1,25,1,1` makes the 25 basemat nodes interaction nodes; `POINT,0,0,2.7` is a surface foundation
  with central-zone radius $0.9 \times 3\,\text{m}$.

### Why it matters
Real buildings have eccentric masses: the centre of mass of a floor rarely coincides with the centre
of rigidity. An eccentricity turns a horizontal input into torsion and a response in the other
horizontal direction. A design that combines directions only "in the input direction" misses this;
the combination rules exist precisely to collect these cross terms.

### Technical basis
The slab inertias follow from a uniform 12 m × 12 m slab of 300 t:

```math
\begin{aligned}
I_Z &= \frac{m\,(a^2 + b^2)}{12} = \frac{300 \times 288}{12} = 7200\,\text{t\,m}^2\\
I_X = I_Y &= \frac{m\,a^2}{12} = 3600\,\text{t\,m}^2
\end{aligned}
```

with $m = 300\,\text{t}$ the mass and $a = b = 12\,\text{m}$ the sides of the slab. `MT`/`MR` take
weights ($I \times g$ with the default `MUNITS` 1, divided by $g$ in HOUSE). The stick and the arm are
2-node 3D Timoshenko frames (shear areas of `R`), the basemat flat shells
([Theory §13.1](docs/theory/THEORY_MANUAL.md#131-the-element-library)).

### In ANSYS terms
`MT` and `MR` are a MASS21 element with all six entries; the rigid arm plays the part of a rigid
link (CERIG or MPC184). The stick is BEAM188 with an arbitrary section (`SECTYPE,...,BEAM,ASEC`,
`SECDATA` for the area and inertias of `R,1`, `SECCONTROL` for its shear areas), which is what the
`ANSYS` export command writes (lesson 08).

```action
plot-model
explain: MR,30,35316,35316,70632
```

## Three wave fields: SV for X, SH for Y, P for Z

SITE runs three times. Each run computes the free field of one vertically propagating wave, normalised
to a unit motion in one direction at the control point (the free surface), and `FCOPY` keeps it under
its own name.

```sassi
* SITE,<opmode>,<mode1>,<fstep>,<nl>,<hs>,<mode2>,<wopt>,<freq1>,<freq2>,<cl>,<cm>,<delt>,<nft>,<freq>
AOPT,0,0,0,1,0,0,0,0,0,0,0,0,0,0
* run 1: SV wave, control motion along x' (cm 0)
SITE,0,1,0,20,3,1,0,1,4096,1,0,0.005,8192,1
WAVE,2,1,1,1,0
CHECK
AFWRITE
RUNSITE
FCOPY,FILE1,FILE1X
* run 2: SH wave, control motion along y' (cm 1); Mode 1 off, FILE2 of run 1 is reused
SITE,0,0,0,20,3,1,1,1,4096,1,1,0.005,8192,1
WAVE,2,0
WAVE,4,1,1,1,0
AFWRITE
RUNSITE
FCOPY,FILE1,FILE1Y
* run 3: P wave, control motion along z' (cm 2)
SITE,0,0,0,20,3,1,0,1,4096,1,2,0.005,8192,1
WAVE,4,0
WAVE,3,1,1,1,0
AFWRITE
RUNSITE
FCOPY,FILE1,FILE1Z
```

### What this does
* Run 1: `<wopt>` 0 (in-plane wave family), `<cm>` 0 ($x'$), `WAVE,2,1,1,1,0` = SV only, vertical
  incidence (angle 0). Mode 1 computes the layer eigen-solutions (FILE2) once.
* Run 2: `<mode1>` 0 reuses FILE2 (it holds both the Rayleigh and the Love modes), `<wopt>` 1
  (anti-plane family), `<cm>` 1 ($y'$); `WAVE,2,0` removes the SV wave and `WAVE,4,...` sets SH.
* Run 3: `<cm>` 2 ($z'$), `WAVE,3,...` = P wave.

### Why it matters
This is the standard three-component input of ASCE and USNRC practice as the ACS SASSI manual
describes it: vertically propagating SV, SH and P waves, which reproduce a uniform free-field motion in
each direction over the plan of the foundation. A vertically propagating wave arrives at every point of
the basemat at the same time, so it excites translation and (through the structure's response) rocking
and torsion, but no torsion from the input itself; inclined and surface waves, incoherency and wave
passage are separate options (lesson 10).

### Technical basis
For vertical incidence the SV and SH free fields are horizontal shear waves in $x'$ and $y'$, the P
field a compression wave in $z'$; each is scaled so that the motion at the control point is 1 in its
direction. [Theory §6.1](docs/theory/THEORY_MANUAL.md#61-vertically-propagating-body-waves) and
[§6.2](docs/theory/THEORY_MANUAL.md#62-normalisation-to-the-control-point);
[User Guide §10.1](docs/user/USER_GUIDE.md#101-simultaneous-cases).

```figure
wave-types fp=4
The three wave fields of this step in a uniform column of the lesson's sand. Each pulse travels
vertically and doubles at the free surface; only the particle motion differs: $x$ for SV, $y$ for
SH, $z$ for P, which travels twice as fast here ($V_p = 2V_s$). Every point of a horizontal plane
moves in phase, which is why these waves give no torsional input (Check yourself).
```

### Check yourself
With these vertically propagating waves, does a doubly symmetric building (no eccentricity) twist?

Answer: No. A vertically propagating wave moves every point of the basemat in phase, so it applies
no rotational input about the vertical axis; torsion comes only from eccentricities of mass or
stiffness. In this model the X and Y inputs twist the stick only because of the eccentric roof mass,
and the Z input produces no torsion at all (`Z_00027R_ZZ.TFU` is zero to round-off after step 6).
So this input does not cover the accidental torsion that design practice adds for the uncertain
location of mass and stiffness; incoherency or wave passage, which do produce rotational input, are
separate options (lesson 10).

```action
open-listing: SITE
explain: SITE,0,0,0,20,3,1,1,1,4096,1,1,0.005,8192,1
```

## One ANALYS for three directions (simultaneous cases)

POINT and HOUSE run once. `ANALYS ... <simul>` = 1 reads FILE1X, FILE1Y and FILE1Z and solves the
three cases as three right-hand sides of one factorisation per frequency.

```sassi
* ANALYS,<opmode>,<type>,<mode>,<save>,<prnt>,<fopt>,<ang>,<xc>,<yc>,<zc>,<impe>,<simul>
ANALYS,0,0,0,0,1,0,0,0,0,0,0,1
AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0
AFWRITE
RUNPOINT
RUNHOUSE
RUNANALYS
```

### What this does
`<simul>` = 1 is the "Simultaneous Cases" option for seismic input: the result is FILE8X, FILE8Y and
FILE8Z, the transfer functions of every degree of freedom for a unit control motion in X, Y and Z.
Everything else is the seismic initiation run of lesson 04.

### Why it matters
The expensive part of an SSI run is the impedance of the soil at the interaction nodes and the
factorisation of the system, frequency by frequency. They do not depend on the input direction, so
three directions cost little more than one. The ANALYS listing's **low-frequency check** should show
$\lvert\mathrm{ATF}\rvert$ close to 1 in the input direction of each case: at $f \to 0$ the whole
system moves with the free field. Here the largest deviation at 0.098 Hz is 0.028 % (FILE8X), 0.030 %
(FILE8Y) and 0.006 % (FILE8Z), all at node 30.

### Technical basis
Per frequency ANALYS solves $[C]\,U = Q$ with the flexible-volume matrix $C$ (structure plus soil
impedance) and one load vector $Q$ per case; VP-23 checks that each simultaneous case equals a single
run to 1e-12. [Theory §9.4](docs/theory/THEORY_MANUAL.md#94-simultaneous-cases-and-the-low-frequency-check).

### In ANSYS terms
Like solving several load cases with one factorised stiffness matrix (multiple load vectors in one
solve), here in a harmonic analysis at each frequency.

```action
open-listing: ANALYS
explain: ANALYS,0,0,0,0,1,0,0,0,0,0,0,1
```

## The X input: direct and coupled response

MOTION convolves the transfer functions of FILE8X with the record. Besides the X response of the
basemat, floor, stick top and roof mass, the requests ask for the **Y and Z response of the roof mass**
and the torsion of the stick top under the X input.

```sassi
DAMP,0.05
THFILE,../data/rg160h_030g.acc
THTIT,RG 1.60 horizontal spectrum-compatible motion (EQUAKE, seed 11975)
AOPT,0,0,0,0,0,0,0,0,0,0,1,0,0,0
MOTION,0,0,0,20,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1
EDUOPT,TFFILE,FILE8X
* NOUT,<dir>,<TF print>,<save TH>,<plot TH>,<plot RS>,<save RS>,<print max>,<nodes>
NOUT,1,1,1,0,0,1,1,13,26,27,30
NOUT,2,1,0,0,0,1,1,30
NOUT,3,1,0,0,0,1,1,30
NOUT,6,1,0,0,0,0,1,27
AFWRITE
RUNMOTION
* MOTION names files by node and direction only: keep this run's files under X_ names
VAR,R30,00030TR_X01.RS,00030TR_Y01.RS,00030TR_Z01.RS,00030TR_X.TFU,00030TR_Y.TFU,00030TR_Z.TFU,00027R_ZZ.TFU
FOREACH,R30,FCOPY,@R30[#],X_@R30[#]
* the steady-state motion under the X input at 5.0 Hz, 24 frames over one period
HARMFRAME,FILE8X,5.0,HARM_X5
```

### What this does
* `EDUOPT,TFFILE,FILE8X` makes MOTION read the X case; AFWRITE writes the control direction x' into
  the MOTION deck from the file name.
* `NOUT,2,...,30` and `NOUT,3,...,30` request the Y and Z response of node 30 under this X input;
  `NOUT,6,...,27` the rotation about Z (torsion) of the stick top.
* `00030TR_Y01.RS` means "node 30, translation Y, damping 1" whatever the input direction, so the
  next run would overwrite it: `FOREACH ... FCOPY` keeps copies with the prefix `X_`.
* `HARMFRAME,FILE8X,5.0,HARM_X5` writes the steady-state motion of every node under a harmonic X
  input of unit amplitude at 5.005 Hz (24 frames over one period). In the animation (0.3 m per unit of
  control motion) the roof mass does not move along X only: X 5.00 (39° behind the control motion),
  Y 1.85 (72° ahead of it) and Z 1.25. Out of phase, X and Y make it run round an ellipse in plan, and
  the rigid arm swings about the stick top (Y 1.22): the stick twists. Turn the view with the mouse
  to look at it from above.

### Why it matters
The X input moves the roof mass mostly in X (transfer function peak 5.00 at 5.0 Hz, 1.000 at
0.1 Hz), but also in Y (1.85 at 5.0 Hz, zero at low frequency) and vertically (1.49 at 6.0 Hz: the
mass is off the axis, so the rocking of the stick moves it up and down), and it twists the stick
top (0.37 rad per metre of control motion at 5.0 Hz). The 5 % spectra of the roof mass under this
one input peak at **4.27 g in X, 1.07 g in Y and 1.47 g in Z**. A design that kept only the response
in the input direction would miss a horizontal spectrum of 1.07 g and a vertical one of 1.47 g at the
roof mass.

### Technical basis
$a(t) = \mathrm{IFFT}[H(f)\,A(f)]$ with $H$ the interpolated transfer function and $A$ the FFT of the
record
([Theory §11.1](docs/theory/THEORY_MANUAL.md#111-convolution)); spectra by the Nigam-Jennings
recurrence at 301 frequencies from 0.1 to 100 Hz, the range and number the ACS SASSI manual recommends
for nuclear ISRS per SRP 3.7.1
([Theory §11.3](docs/theory/THEORY_MANUAL.md#113-response-spectra-nigam-jennings), VP-30).

```action
plot-spectrum: ex05/X_00030TR_X01.RS, ex05/X_00030TR_Y01.RS, ex05/X_00030TR_Z01.RS | log
plot-spectrum: ex05/X_00030TR_Y.TFU, ex05/X_00027R_ZZ.TFU
animate: ex05/HARM_X5 | deformed 0.3 | X input at 5.0 Hz
open-listing: MOTION
```

## The Y and Z inputs, and a relative displacement

The same requests for the Y case (FILE8Y), with RELDISP for the drift of the stick top relative to the
free field, and for the Z case (FILE8Z) with the vertical motion scaled to 2/3 of the horizontal one.

```sassi
* ---- Y input
EDUOPT,TFFILE,FILE8Y
NOUT,0
NOUT,1,1,0,0,0,1,1,30
NOUT,2,1,1,0,0,1,1,13,26,27,30
NOUT,3,1,0,0,0,1,1,30
NOUT,6,1,0,0,0,0,1,27
* RELDISP: Y displacement of node 27 relative to the free field (no RELFILE)
RELD,1,0,0
RDND,27,0,1,0,0,0,0
AOPT,0,0,0,0,0,0,0,0,0,0,1,0,1,0
AFWRITE
RUNMOTION
RUNRELDISP
FOREACH,R30,FCOPY,@R30[#],Y_@R30[#]
* ---- Z input: MOTION <mult> = 0.6667
EDUOPT,TFFILE,FILE8Z
NOUT,0
NOUT,1,1,0,0,0,1,1,30
NOUT,2,1,0,0,0,1,1,30
NOUT,3,1,1,0,0,1,1,13,26,27,30
NOUT,6,1,0,0,0,0,1,27
MOTION,0,0,0,20,0,0.1,100,301,0.6667,0,1,0,0,0,0,1,0,0,1
AOPT,0,0,0,0,0,0,0,0,0,0,1,0,0,0
AFWRITE
RUNMOTION
FOREACH,R30,FCOPY,@R30[#],Z_@R30[#]
```

### What this does
* `NOUT,0` clears the previous requests. Each case asks for the X, Y and Z response of node 30, so
  that the three co-directional contributions of every response direction exist after the three runs.
* `RELD,1,0,0` with no `RELFILE` uses the free-field reference; `RDND,27,0,1,0,0,0,0` requests the
  Y component of node 27. Result: `00027TR_Y.THD`.
* `MOTION ... <mult>` = 0.6667 scales the record for the vertical case.

### Why it matters
* **Y input.** The Y transfer function of the roof mass peaks at 6.53 at 5.0 Hz, higher than the X one
  (5.00): the mass is 2 m off the axis in x and only 1 m in y, so the Y input twists the stick twice
  as much (0.74 against 0.37 rad/m) and the twist adds to the Y motion of the mass. The cross terms
  are equal: X response to Y input = Y response to X input = 1.85. Its 5 % spectrum peaks at
  4.25 g at 5.0 Hz.
* **Z input.** The vertical transfer function of the roof mass peaks at 1.67 at 10.0 Hz; through the
  eccentric mass the vertical input also produces X (1.03) and Y (0.51) motion at 6.0 Hz, and no
  torsion at all (the stick-top rotation about Z is zero to round-off).
* **RELDISP.** The stick top drifts up to **5.2 mm** in Y relative to the free field
  ($t = 3.6\,\text{s}$):
  the interstorey and building-to-building displacements for gaps and for distribution systems that
  span between structures come from RELDISP, not from double integration of the accelerations.

The example uses **one record** for all three directions (and 2/3 of it vertically) to keep the lesson
short. A design analysis uses three statistically independent records (EQUAKE generates them and
reports their cross-correlation, required $\lvert\rho\rvert \le 0.16$), and the vertical record
matches the vertical design spectrum. The SRSS combination of the next step assumes that independence.

### Technical basis
$d(t) = \mathrm{IFFT}[(H_\text{node} - H_\text{ref})\,U_g]$ with $U_g = -g\,A/\omega^2$ and the free
field as reference (unit amplitude in the input direction): the drift is computed from the transfer
functions, without double integration
([Theory §11.4](docs/theory/THEORY_MANUAL.md#114-relative-displacements-reldisp), VP-34).

```action
plot-spectrum: ex05/X_00030TR_Y01.RS, ex05/Y_00030TR_Y01.RS, ex05/Z_00030TR_Y01.RS | log
plot-history: ex05/00027TR_Y.THD
open-listing: RELDISP
```

## Combine the three directions: SRSS and 100-40-40

For the Y response of the roof mass there are now three spectra: due to the X, the Y and the Z input.
The ACS SASSI manual offers two combinations of co-directional ISRS: the SRSS and a weighted linear
combination. The 100-40-40 rule is the weighted linear combination that the LINECOMBIN specification
gives as its example ($1.0\,X + 0.4\,Y + 0.4\,Z$); which rule applies is set by your design basis.

```sassi
READSPEC,X_00030TR_Y01.RS,1,1
READSPEC,Y_00030TR_Y01.RS,1,2
READSPEC,Z_00030TR_Y01.RS,1,3
LINENAME,1,Y from X input
LINENAME,2,Y from Y input
LINENAME,3,Y from Z input
* SRSS of the three contributions
SRSS,4,1,2,3
LINENAME,4,SRSS
* 100-40-40: each contribution in turn at 100 %, the others at 40 %, then the envelope
LINECOMBIN,5,2,1.0,1,0.4,3,0.4
LINECOMBIN,6,1,1.0,2,0.4,3,0.4
LINECOMBIN,7,3,1.0,1,0.4,2,0.4
BROADEN,8,0,0,5,6,7
LINENAME,8,100-40-40
WRITESPEC,roof_Y_SRSS.rs,4
WRITESPEC,roof_Y_1004040.rs,8
```

### What this does
* `SRSS,4,1,2,3`: line 4, $y_4 = \sqrt{y_1^2 + y_2^2 + y_3^2}$, at every frequency (on the union of the
  frequencies), with $y_k$ the ordinate of line $k$.
* `LINECOMBIN,5,2,1.0,1,0.4,3,0.4`: line 5, $y_5 = 1.0\,y_2 + 0.4\,y_1 + 0.4\,y_3$; lines 6 and 7 put
  the other two contributions at 100 %; `BROADEN,8,0,0,5,6,7` (no broadening, no bridging) is their
  envelope.
* `WRITESPEC` writes the two combined spectra to `ex05/`.

### Why it matters
The Y response of the roof mass collects 4.25 g from the Y input, 1.07 g from the X input and 0.34 g
from the Z input (spectral peaks). The SRSS peaks at **4.36 g** at 5.0 Hz, 2.6 % above the Y
contribution alone; the 100-40-40 envelope peaks at **4.76 g** at 5.1 Hz, 9 % above the SRSS. At
100 Hz (the ZPA) they give 0.736 g and 0.798 g. The cross-direction terms are not negligible for an
eccentric structure, and the choice of rule is a design-basis decision: state it with the ISRS.

### Technical basis
Both rules combine peak responses that do not occur at the same time; they are justified when the
three input components are statistically independent. SRSS estimates the peak of the sum of three
uncorrelated responses; 100-40-40 is a linear alternative that usually gives somewhat larger values
(here 9 % at the peak). For time histories the manual offers the algebraic sum of the
co-directional histories, time step by time step: the response to the three components acting
together, which needs no combination rule (Option NON does this in lesson 09). Line operations:
[User Guide §11.5](docs/user/USER_GUIDE.md#115-line-mathematics-and-spectrum-broadening) and
[spec 10 §4.1](docs/spec/10_commands_module_calc_plot_prog.md) (SRSS, LINECOMBIN), verified by
[VP-47](docs/verification/VERIFICATION_MANUAL.md#vp-47).

### Check yourself
When does 100-40-40 exceed SRSS the most, and when do the two rules agree?

Answer: with the largest contribution $R_1$ and the two others fractions $a$ and $b$ of it,
$R_\text{SRSS} = R_1\sqrt{1 + a^2 + b^2}$ and $R_\text{100-40-40} = R_1\,(1 + 0.4a + 0.4b)$. They
agree when one direction acts alone ($a = b = 0$), differ by 4 % when the three are equal (1.80
against 1.73), and by up to 15 % when the two others are each about 40 % of the dominant one (1.32
against 1.15). For the Y response of the roof mass the X input gives about 25 % and the Z input 8 % of
the Y input's peak, hence the 9 %. Try the X response yourself with `X_00030TR_X01.RS`,
`Y_00030TR_X01.RS` and `Z_00030TR_X01.RS`.

```action
plot-spectrum: ex05/X_00030TR_Y01.RS, ex05/Y_00030TR_Y01.RS, ex05/Z_00030TR_Y01.RS, ex05/roof_Y_SRSS.rs, ex05/roof_Y_1004040.rs | log
explain: LINECOMBIN,5,2,1.0,1,0.4,3,0.4
```

## The design ISRS: envelope and peak broadening

The design spectrum of a floor is the combined spectrum, enveloped over the soil cases and with its
peaks broadened to cover the uncertainty of the structure's frequencies. `BROADEN` does all of it:
envelope of its sources, $\pm\text{Smooth2}\,\%$ peak broadening, valley bridging within Smooth1 %.

```sassi
* BROADEN,<dest>,<Smooth1 = bridging %>,<Smooth2 = broadening %>,<sources>
BROADEN,9,0,15,4
LINENAME,9,roof mass Y design ISRS 5 % (+-15 %)
WRITESPEC,roof_Y_design.rs,9
```

### What this does
`BROADEN,9,0,15,4` broadens the SRSS spectrum (line 4) by ±15 % without bridging and writes it to
`ex05/roof_Y_design.rs`. With several soil cases you would list the SRSS line of every case as
sources (`BROADEN,9,0,15,4,14,24`): the envelope is taken first, then broadened.

### Why it matters
The broadened spectrum keeps the SRSS peak, 4.36 g, over a plateau from 4.26 to 5.76 Hz
($5.01\,\text{Hz} \pm 15\,\%$), and leaves the ZPA practically unchanged (0.737 g against 0.736 g). A
component or a piping mode anywhere in that band is qualified for the peak. Broadening covers the
uncertainty of the structure's frequencies (material properties, modelling, mass), which the soil
cases alone do not.

This file, for every floor, direction and damping value, is what the equipment and piping engineers
qualify against. Typical deliverables use several damping values (`DAMP,0.02,0.03,0.05,...`, up to
ten): run MOTION once with the full list.

### Technical basis
A peak at $f_0$ becomes a plateau over $[f_0(1 - b),\ f_0(1 + b)]$, $b = \text{Smooth2}/100$; the
flanks are shifted by $(1 - b)$ and $(1 + b)$, evaluated on an augmented grid so that the plateau
edges are exact. The ACS SASSI manual names the BROADEN parameters but gives neither a formula nor a
broadening value; the BROADEN specification
([spec 10 §4.1.3](docs/spec/10_commands_module_calc_plot_prog.md)) lists, as typical values from
outside the manual, ±15 % for ASCE 4 ISRS practice and a frequency-dependent broadening of at least
±10 % for RG 1.122. The ±15 % used here follows that note; your design basis sets the value. The
manual warns that fewer than 301 spectrum frequencies reduce the accuracy of the broadening (EDU-15
warns). [User Guide §11.5](docs/user/USER_GUIDE.md#115-line-mathematics-and-spectrum-broadening);
verified by [VP-48](docs/verification/VERIFICATION_MANUAL.md#vp-48) (T-B1 to T-B3).

```figure
isrs-broadening b=0.15 cases=1
`BROADEN` as a picture: the envelope of the sources is taken first, then every peak becomes a
plateau from $f_p(1-b)$ to $f_p(1+b)$ and the flanks move outwards. In this step the SRSS peak at
5.01 Hz becomes a plateau from 4.26 to 5.76 Hz ($b = 15\,\%$); switch the figure to three soil
cases to see the envelope come first.
```

### In ANSYS terms
The same as the ISRS post-processing you do after an ANSYS transient run (or a spectrum-to-spectrum
method): only here the floor motion already contains the soil's frequency shift and radiation
damping.

### Try this
Make the design spectrum of the roof mass in X: `READSPEC` the three `?_00030TR_X01.RS` files into
lines 11-13, `SRSS,14,11,12,13`, `BROADEN,15,0,15,14`, `WRITESPEC,roof_X_design.rs,15`. Then try a
bridging of 10 % (`BROADEN,16,10,15,14`) and see which valleys it fills.

```action
plot-spectrum: ex05/roof_Y_SRSS.rs, ex05/roof_Y_design.rs | log
open-doc: docs/user/USER_GUIDE.md#115-line-mathematics-and-spectrum-broadening
```
