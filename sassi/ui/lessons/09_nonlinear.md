---
id: 09-nonlinear
title: Nonlinear soil and cracked concrete by iteration
part: Advanced
order: 9
minutes: 50
example: ex06_nonlinear_soil
summary: The frequency-domain SSI solution is linear; soil and concrete nonlinearity enter through equivalent-linear iterations. Run the free-field (primary) and near-field (secondary) soil iterations of example 6, then the cracking shear walls of example 7 with Option NON.
objectives: [Tell primary from secondary soil nonlinearity and know when each matters, Model a near-field soil zone as a nonlinear soil group and iterate it with New Structure restarts, Define wall panels with capacities and backbone curves for Option NON, Run the Option NON iterations and read the cracked state and its convergence, Judge convergence and the limits of equivalent linearisation]
prerequisites: [05-embedded, 06-seismic-input]
---
SASSI solves the SSI equation frequency by frequency with complex stiffnesses: it is a linear method.
You know the consequence from ANSYS: a harmonic analysis cannot crack concrete or yield soil. Yet both
happen in a design-basis earthquake, and both change the result:

* **soil** softens and dissipates energy as it strains. In the **free field** (primary nonlinearity)
  SOIL takes care of it (lesson 06). Next to the structure (**secondary**, near-field nonlinearity)
  the soil can strain much more than the free field: loose backfill against a wall, soil under a heavy
  neighbouring building;
* **reinforced-concrete walls** crack: the effective stiffness of a shear wall can drop well below its
  uncracked value and its damping rises; the building's frequencies fall and the ISRS shift.

SASSI handles both the way SHAKE handles a soil column: replace each nonlinear element by an
**equivalent-linear** one whose stiffness and damping depend on how much it deforms, run the SSI
analysis, update the properties from the computed deformations, and repeat until they stop changing.
Each repetition is a cheap **New Structure restart**: the soil impedance of the first run is reused.

The lesson runs example 6 (a loose backfill block behind an abutment wall, about 5 s) and then builds
and runs example 7 (a two-storey shear-wall building with eight cracking wall panels at 0.6 g, about
15-60 s depending on the computer). The iterations need the ANALYS restart files, about 55 MB for
example 7; the last step deletes them once the iterations are over, and the workspace keeps about
50 MB.

## Primary nonlinearity: the free-field column (SOIL)

Example 6 sits in 10 m of medium-dense sand ($V_s = 250\,\text{m/s}$ at low strain) on rock. SOIL runs
the equivalent-linear site response with the RG 1.60 record scaled to a 0.30 g rock outcrop.

```sassi
MDL,ex06,ex06
TIT,Ex06 - loose backfill behind an embedded wall (equivalent-linear SSI iterations)
* low-strain sand (sublayers of 1, 1, 2, 2, 2, 2 m) on rock
L,1,1.0,19.0,500,250,0.02,0.02
L,2,1.0,19.0,500,250,0.02,0.02
L,3,2.0,19.0,500,250,0.02,0.02
L,7,1.0,21.0,1600,800,0.01,0.01
SITE,0,1,0,20,7,1,0,1,4096,1,0,0.005,8192,1
WAVE,2,1,1,1,0
HOUSE,9.81,0,0,2,0,0,0,0,0
INP,../sassi/data/dynp_library.pre
SPRO,1,1,Sand
SPRO,2,2,Sand
VAR,SANDL,3,4,5,6
FOREACH,SANDL,SPRO,@SANDL[#],3,Sand
SPRO,7,7
SOIL,4000,9.81,1,1,1,8,0.65,1,0
* SOILX,<indir>,<mult>,<max>,<cl>,<file>: rock outcrop scaled to 0.30 g at the top of sublayer 7
SOILX,0,0,0.3,7,../data/rg160h_030g.acc
EDUOPT,SOILCUTOFF,15
* the surface motion ACC001.TH (the SSI control motion) and its 5 % spectrum
SACC,1,2,0
SRS,1,1,0
DAMP,0.05
AOPT,0,1,0,0,0,0,0,0,0,0,0,0,0,0
CHECK
AFWRITE
RUNSOIL
```

### What this does
The chain of lesson 06 in short: `SPRO` assigns the Sand curve of the SHAKE91 library to the six
sublayers, `SOILX,0,0,0.3,7,...` scales the record to a 0.30 g peak (`<max>`) and applies it as an
outcrop motion at the top of the rock, `EDUOPT,SOILCUTOFF,15` removes Fourier components above
15 Hz (the cut-off of the SSI run that follows), `SACC,1,2,0` saves the surface motion `ACC001.TH`.
`SRS,1,1,0` (with `DAMP,0.05`) adds its 5 % spectrum `RS001_01.RS` for comparison in step 6.

### Why it matters
The 0.30 g rock outcrop becomes **0.58 g** at the surface. The strain-compatible sand softens with
depth, from $V_s = 245\,\text{m/s}$ ($G/G_{\max} = 0.96$, 1.4 % damping) at 0.5 m to 148 m/s
($G/G_{\max} = 0.35$, 10.3 %) at 9 m; the 5 % surface spectrum peaks at 2.20 g at 3.98 Hz. These
properties are the free field of the SSI model of the next step, and `ACC001.TH` is its control
motion. This is the primary nonlinearity: it is in every SSI analysis whose soil properties come from
SOIL.

### Technical basis
Equivalent-linear SHAKE iterations, $G = G_{\max}\,(G/G_{\max})(\gamma_\text{eff})$,
$\beta = D(\gamma_\text{eff})$, $\gamma_\text{eff} = 0.65\,\gamma_{\max}$
([Theory §12](docs/theory/THEORY_MANUAL.md#12-the-shake-equivalent-linear-method-soil), VP-04).

```action
open-listing: SOIL
plot-history: ex06/ACC001.TH
```

## The SSI model: strain-compatible layers and the excavated soil

The SSI site takes the strain-compatible properties of SOIL (FILE88, rounded) as L layers 11-16.
The wall's excavation is a 4 m × 4 m × 2 m block behind the wall, meshed with 2 × 2 × 2 SOLID
elements; every node of it is an interaction node (flexible volume).

```sassi
* strain-compatible site from FILE88 (each embedment layer its own L number)
L,11,1.0,19.0,489,245,0.014,0.014
L,12,1.0,19.0,459,230,0.029,0.029
L,13,2.0,19.0,410,205,0.048,0.048
L,14,2.0,19.0,361,181,0.072,0.072
L,15,2.0,19.0,324,162,0.089,0.089
L,16,2.0,19.0,297,148,0.103,0.103
TOPL,11,12,13,14,15,16
FREQ,1,4,20,41,61,82,102,123,143,164,184
FREQ,1,205,225,246,266,287,307,328,369,410,451
FREQ,1,492,533,573,614
* nodes at x = 0, 2, 4; y = -2, 0, 2; z = -2, -1, 0, numbered bottom-up
N,1,0,-2,-2
N,3,4,-2,-2
FILL,1,3
NGEN,2,3,1,3,1,0,2,0
NGEN,2,9,1,9,1,0,0,1
* excavated soil, one group per embedment layer (MACT = its L number)
GROUP,1,SOLID
GTIT,1,excavated soil 0 to -1 m
MACT,11
E,1,10,11,14,13,19,20,23,22
EGEN,1,1,1
EGEN,1,3,1,2
ETYPE,1,4,1,2
GROUP,2,SOLID
GTIT,2,excavated soil -1 to -2 m
MACT,12
E,1,1,2,5,4,10,11,14,13
EGEN,1,1,1
EGEN,1,3,1,2
ETYPE,1,4,1,2
INT,1,27,1,1
```

### What this does
* L 11-16 are the SOIL results rounded (for example L 11: $V_s = 245\,\text{m/s}$, 1.4 % damping;
  L 16: 148 m/s, 10.3 %). Each embedment layer has its own L number, and each excavated group uses the L
  number of its layer (`MACT,11`, `MACT,12`).
* `ETYPE,1,4,1,2` marks elements 1-4 of the active group as **excavated soil** (type 2): they are the
  soil that the flexible-volume method removes; their properties must equal the free field.
* `INT,1,27,1,1`: all 27 nodes of the block are interaction nodes (FV).

### Why it matters
The excavated soil and the free field must be the same material, at the same strain level: if the
excavation kept the low-strain properties while the free field softened, the "removed" soil would not
be the soil that SITE and POINT computed, and the substructuring would be inconsistent. The example
copies the FILE88 values into the L table by hand (SITEX,1 would update the free field but not the L
table the excavated elements use).

### Technical basis
Flexible volume: the dynamic stiffness of the structure minus the excavated soil,
$(K^*_s - \omega^2 M_s) - (K^*_e - \omega^2 M_e)$, plus the free-field impedance $X_{ff}$ at the
interaction nodes ([Theory §3](docs/theory/THEORY_MANUAL.md#3-flexible-volume-substructuring),
lesson 05). Excavated SOLIDs take the L-table layer of their `MSET` index
([Theory §13.1](docs/theory/THEORY_MANUAL.md#131-the-element-library)).

```action
open-file: ex06/FILE88
plot-layers
```

## The backfill as a nonlinear soil group, and the wall

The loose backfill ($V_s = 160\,\text{m/s}$ at low strain) occupies the same volume as the excavated
native soil. It is modelled as SOLID elements of the **structure** (type 1) at the same nodes: group
3, whose shear modulus and damping HOUSE will iterate. The wall is a 0.6 m SHELL in the plane $x = 0$
carrying 150 t of a bridge deck at its top.

```sassi
* the backfill: LOW-STRAIN properties (G_max of the curve): Vp 320, Vs 160 m/s, 2 % damping
M,2,320,160,18.5,0.02,0.02,3
GROUP,3,SOLID
GTIT,3,backfill (non-linear soil)
MACT,2
E,1,1,2,5,4,10,11,14,13
EGEN,1,1,1
EGEN,1,3,1,2
EGEN,1,9,1,4
ETYPE,1,8,1,1
* the wall (concrete, 0.6 m) in the plane x = 0
M,1,3.0E7,0.2,24.0,0.05,0.05,1
GROUP,4,SHELL
GTIT,4,wall
MACT,1
E,1,1,4,13,10
EGEN,1,3,1
EGEN,1,9,1,2
THICK,1,4,1,0.6
D,1,25,3,1,ROTX
* 150 t of the deck on the three top nodes of the wall
MT,19,490,490,490
MT,22,490,490,490
MT,25,490,490,490
* the non-linear soil input (.pin): ESF 0.65; group 3, octahedral strain, start at 0.4 G_ff, curve Sand
PIN,0.65
PINGRP,3,1,0.4,1.0,Sand
HOUSEX,0,0,1,1
PINLIST
```

### What this does
* `M,2,320,160,18.5,0.02,0.02,3` (type 3: $V_p$, $V_s$) holds the **low-strain** backfill:
  $G_{\max} = \rho V_s^2 = 48\,\text{MPa}$ is the reference of the $G/G_{\max}$ curve. Never enter
  strain-compatible values here: the curve would soften them a second time.
* `PIN,0.65` sets the effective-strain factor; `PINGRP,3,1,0.4,1.0,Sand` declares group 3 nonlinear
  with the octahedral shear strain (ISTR 1), starting at $\mathrm{GFAC} = 0.4$ times the free-field
  $G$ at the element and $\mathrm{DFAC} = 1$ times its damping, iterated on the Sand curve.
  `HOUSEX,0,0,1,1` switches on non-linear SSI in HOUSE; `PINLIST` lists the input.
* `D,1,25,3,1,ROTX` fixes the drilling rotation (about the wall normal X) at the wall nodes (every
  third node from 1 to 25). Nothing else restrains it: the flat shell has no drilling stiffness and
  the solids that share these nodes have no rotations. `FIXROT` does this only at nodes connected to
  shells alone, and these nodes also belong to the excavated and backfill solids, hence the explicit
  `D`.

### Why it matters
Near-field nonlinearity is the soil effect a free-field site response cannot see. A loose backfill
behind a wall softens and damps more than the native soil and changes what the wall and its deck
feel: wall pressures, the deck acceleration, the relative displacement of wall and backfill. Without
a model of it you would analyse the wall in the free-field soil.

### Technical basis
HOUSE sets the properties of the nonlinear elements as

```math
\begin{aligned}
&\text{iteration 0:} && G = \mathrm{GFAC}\cdot G_\text{layer}, && \beta = \mathrm{DFAC}\cdot\beta_\text{layer}\\
&\text{later iterations:} && G = G_{\max}\,(G/G_{\max})(\gamma_\text{eff}), && \beta = D(\gamma_\text{eff})
\end{aligned}
```

with $G_\text{layer}$, $\beta_\text{layer}$ the free field at the element's mid-depth and $G_{\max}$
from the M table. The strain of ISTR 1 is the octahedral shear strain

```math
\gamma_\text{oct} = \frac{2}{3}\sqrt{(\varepsilon_x-\varepsilon_y)^2 + (\varepsilon_y-\varepsilon_z)^2 + (\varepsilon_z-\varepsilon_x)^2 + 1.5\,(\gamma_{xy}^2 + \gamma_{yz}^2 + \gamma_{xz}^2)}
```

which gives $\sqrt{2/3}\,\gamma$ in simple shear $\gamma$. GFAC only chooses the start: the converged
state is set by $G_{\max}$, the curve and the motion.
[Theory §16](docs/theory/THEORY_MANUAL.md#16-nonlinear-soil-ssi-iterations).

### In ANSYS terms
A group of SOLID185 elements whose material is updated between linear runs, like a secant-stiffness
iteration that you would script around a harmonic analysis.

```action
plot-model
explain: PINGRP,3,1,0.4,1.0,Sand
```

## Iteration 0: the initiation run and the effective strains

The initiation run is an ordinary three-direction SSI analysis (lesson 07) that saves its restart
files. STRESS then computes the effective strain of every backfill element for each input direction,
and `COMBXYZSTRAIN` combines them.

```sassi
POINT,0,2,1.8
THFILE,ACC001.TH
THTIT,SOIL surface motion (rock outcrop 0.30 g)
MOTION,0,0,0,20,0,0.1,100,301,1,0,1,5000,0,0,0,1,0,0,1
* STRESS,<opmode>,<iter>,...: <iter> 1 = strains of the nonlinear soil elements
STRESS,0,1,0,0,1
AOPT,0,0,0,1,0,0,0,0,0,0,0,0,0,0
CHECK
AFWRITE
RUNSITE
FCOPY,FILE1,FILE1X
SITE,0,0,0,20,7,1,1,1,4096,1,1,0.005,8192,1
WAVE,2,0
WAVE,4,1,1,1,0
AFWRITE
RUNSITE
FCOPY,FILE1,FILE1Y
SITE,0,0,0,20,7,1,0,1,4096,1,2,0.005,8192,1
WAVE,4,0
WAVE,3,1,1,1,0
AFWRITE
RUNSITE
FCOPY,FILE1,FILE1Z
* initiation with Save Restart Files (<save> 1) and X/Y/Z cases (<simul> 1)
ANALYS,0,0,0,1,1,0,0,0,0,0,0,1
NLSSIRESET
AOPT,0,0,0,1,1,1,0,0,1,0,0,1,0,0
CHECK
AFWRITE
RUNPOINT
RUNHOUSE
RUNANALYS
* the STRESS commands of each direction, kept in variables for the iterations
VAR,NLX,"EDUOPT,TFFILE,FILE8X","MOTION,0,0,0,20,0,0.1,100,301,1,0,1,5000,0,0,0,1,0,0,1",AFWRITE,RUNSTRESS,"FCOPY,FILE74,FILE74X"
VAR,NLY,"EDUOPT,TFFILE,FILE8Y",AFWRITE,RUNSTRESS,"FCOPY,FILE74,FILE74Y"
VAR,NLZ,"EDUOPT,TFFILE,FILE8Z","MOTION,0,0,0,20,0,0.1,100,301,0.6667,0,1,5000,0,0,0,1,0,0,1",AFWRITE,RUNSTRESS,"FCOPY,FILE74,FILE74Z"
VAR,NLCOMB,"COMBXYZSTRAIN,FILE74X,FILE74Y,FILE74Z,FILE74"
AOPT,0,0,0,0,0,0,0,0,0,0,0,1,0,0
FOREACH,NLX,@NLX[#]
FOREACH,NLY,@NLY[#]
FOREACH,NLZ,@NLZ[#]
FOREACH,NLCOMB,@NLCOMB[#]
```

### What this does
* The control motion is the SOIL surface motion `ACC001.TH` (8192 values; records 1-5000, the 20 s
  motion plus 5 s, leave a quiet zone). `POINT,0,2,1.8`: two embedded layers, central-zone radius
  $0.9 \times 2\,\text{m}$.
* Three SITE runs (SV, SH, P), `ANALYS ... <save>` 1 (restart files `COOXqqq` with the soil
  impedance) and `<simul>` 1; `NLSSIRESET` starts a new nonlinear analysis (`ex06.liq` = 0), so this
  HOUSE run uses the `.pin` starting properties.
* `VAR,NLX,"...",...` stores command lists (quoted items contain commas); `FOREACH,NLX,@NLX[#]` runs
  them: STRESS for the X case writes the effective strains `FILE74`, copied to `FILE74X`; the same for
  Y and Z (vertical motion 2/3). `COMBXYZSTRAIN` combines them into `FILE74`.

### Why it matters
The nonlinear behaviour must be driven by the simultaneous three-component input, not by each direction
alone: an element that strains in X and Y at the same time is softer than either direction suggests.
`COMBXYZSTRAIN` takes the SRSS of the directional effective strains, the ACS SASSI companion program
COMB_XYZ_STRAIN.

Read the warnings of this step (lesson 10), three of them, and decide:

* SITE, G-05: the deepest sand sublayer (2 m, $V_s = 148\,\text{m/s}$) passes 14.8 Hz, just below the
  15.0 Hz cut-off; SITE also notes that the manual recommends more than 20 soil layers for accurate
  Rayleigh and Love modes (here 6). Both are shortcuts that keep the example small; a design model
  would split the sublayers.
* HOUSE, EXCSTRCHK: interior excavation node 14 is shared with structural SOLIDs. That is the backfill,
  built on the excavation mesh on purpose, and SASSI-EDU accepts it with a warning: for near-field
  soil on the excavation mesh the FV result is close (exact if the elements had the properties of the
  excavated soil), but separate structure and excavation nodes are the correct model (manual §1.5.1
  rule 11). A wall or slab on interior excavation nodes would be a serious modelling error.

### Technical basis
$\gamma_\text{eff} = \mathrm{ESF}\cdot\max_t\lvert\gamma(t)\rvert$ per element and direction (from the
element strain transfer functions, interpolated and convolved like stresses), then
$\gamma_\text{eff} = \sqrt{\gamma_X^2 + \gamma_Y^2 + \gamma_Z^2}$
([Theory §16](docs/theory/THEORY_MANUAL.md#16-nonlinear-soil-ssi-iterations)).

```action
open-listing: STRESS
open-file: ex06/FILE74
```

## Iterate until the backfill is strain-compatible (NLSSIITER)

Each iteration: HOUSE reads FILE74 and the curve, writes the new properties (FILE78) into the structure
matrices; ANALYS restarts with **New Structure** (the soil impedance of the initiation run is reused);
STRESS three times; COMBXYZSTRAIN.

```sassi
* ANALYS <mode> 1 = New Structure restart (X_ff read from COOXqqq)
ANALYS,0,0,1,0,1,0,0,0,0,0,0,1
AOPT,0,0,0,0,0,1,0,0,1,0,0,1,0,0
AFWRITE
VAR,NLRUN,RUNHOUSE,RUNANALYS
* NLSSIITER,<command variables>,<max iterations>
NLSSIITER,NLRUN+NLX+NLY+NLZ+NLCOMB,8
```

### What this does
`NLSSIITER,NLRUN+NLX+NLY+NLZ+NLCOMB,8` runs the commands of the five variables in order, at most
8 times, and stops when the largest change of $G$ is below 2 % and of the damping below 0.5 %
(absolute). Its history is in `NLSOIL_CONVERGENCE.TXT` and in the HOUSE listing.

### Why it matters
The iterations converge in **6 passes** (`NLSOIL_CONVERGENCE.TXT`): $\max\lvert\Delta G/G\rvert$ = 66,
22, 14, 9.4, 4.0, 2.9 and 1.5 % (iterations 0-6), and the damping changes by 4.7, 1.8, 1.3, 0.8, 0.4,
0.3 and 0.1 points. The first change is the largest because $\mathrm{GFAC} = 0.4$ started the backfill
near its own low-strain modulus. The converged backfill (FILE74) has effective octahedral strains of
0.041-0.053 %, $G/G_{\max} = 0.51\text{–}0.57$ ($G = 25\text{–}27\,\text{MPa}$, a fifth to a quarter
of the 102-116 MPa of the native soil at the same depth) and 6.4-7.3 % damping. The four elements next
to the wall ($x = 0$ to 2 m, odd numbers) strain and soften more than the four at the back.

### Technical basis
Row $k$ of `NLSOIL_CONVERGENCE.TXT` compares the properties used in iteration $k$ (FILE78) with the
strain-compatible properties of its response (FILE74):
$\Delta G/G = (G_\text{new} - G_\text{used})/G_\text{new}$,
$\Delta\beta = \beta_\text{new} - \beta_\text{used}$. A New Structure restart reuses the soil
impedance and is, in the manual's guidance, 2-4 times faster than an initiation run
([Theory §16](docs/theory/THEORY_MANUAL.md#16-nonlinear-soil-ssi-iterations)); with unchanged
properties it reproduces the initiation FILE8 exactly
([VP-N2](docs/verification/VERIFICATION_MANUAL.md#vp-n2)).
[VP-N1](docs/verification/VERIFICATION_MANUAL.md#vp-n1) checks the method: a near-field column
identical to the free field, iterated from the low-strain properties, converges to the SHAKE profile
(within 1.3 % on $G/G_{\max}$).

### Check yourself
The convergence test looks only at the change made by one iteration. Why can a slowly converging
sequence pass it before it reaches the strain-compatible state, and how would you check?

Answer: when each iteration changes the properties only a little but always in the same direction, the
change per iteration can fall below 2 % while the sum of the remaining changes is larger. Check the
trend of the history (monotone and shrinking fast, as here, or slow), and compare the last properties
with the curve at the last effective strain (FILE74 gives $\gamma_\text{eff}$, $G$ and $\beta$ per
element).

```action
open-file: ex06/NLSOIL_CONVERGENCE.TXT
open-listing: HOUSE
open-file: ex06/FILE78
```

## The wall top with the converged backfill

MOTION on the X case of the converged model: the deck node at the wall top (22) and the back of the
backfill (24).

```sassi
EDUOPT,TFFILE,FILE8X
MOTION,0,0,0,20,0,0.1,100,301,1,0,1,5000,0,0,0,1,0,0,1
NOUT,1,1,1,0,0,1,1,22,24
AOPT,0,0,0,0,0,0,0,0,0,0,1,0,0,0
AFWRITE
RUNMOTION
```

### What this does
The X-direction transfer functions of the last iteration (FILE8X) convolved with `ACC001.TH`;
5 % spectra (`DAMP,0.05` of step 1) of nodes 22 and 24.

### Why it matters
With the converged backfill the wall top (node 22, carrying the deck) reaches **0.75 g**, while the
back of the backfill (node 24) stays at the free-field 0.58 g. The 5 % spectrum of the wall top peaks
at 2.67 g at 4.27 Hz, against 2.29 g at 3.98 Hz at node 24 and 2.20 g at 3.98 Hz in the free field.
Designing the deck bearing for the free-field motion would underestimate its peak acceleration by
23 %.

**When do you need near-field iterations?** When the soil next to the structure is clearly different
from the free field (backfill, an excavation refilled with engineered fill, soft pockets) or strains
much more than the free field (heavy structures on soft soil, adjacent buildings), and the quantities
that depend on it (wall pressures, the response of components attached to the wall, building
interaction) matter for the design. For a stiff surface mat on competent soil the free-field
strain-compatible properties are usually enough.

### Technical basis
The near-field elements are part of the structure model, so the result depends on their mesh like any
FE model (lesson 10); the iterations use the same curve interpolation as SOIL (linear in $\log\gamma$).
[User Guide §12.4](docs/user/USER_GUIDE.md#124-nonlinear-soil-soil-near-field-iterations-and-soil-non).

```action
plot-spectrum: ex06/RS001_01.RS, ex06/00022TR_X01.RS, ex06/00024TR_X01.RS | log
open-listing: MOTION
```

## Option NON, part 1: a two-storey shear-wall building

Example 7: a 12 m × 12 m box, two storeys of 4 m, 0.3 m reinforced-concrete walls meshed 3 m × 2 m, on
a 1.5 m stiff surface mat, on 10 m of stiff soil ($V_s = 400\,\text{m/s}$) over rock. A second model
in memory (model 1) keeps it apart from example 6.

```sassi
ACTM,1
MDL,ex07,../ex07
TIT,Ex07 - Option NON: two-storey RC shear-wall building, nonlinear wall panels
L,1,0.4,19.6,800,400,0.05,0.05
L,2,1.0,22.0,2400,1200,0.02,0.02
TOPL,1,1,1,1,1,1,1,1,1,1
TOPL,1,1,1,1,1,1,1,1,1,1
TOPL,1,1,1,1,1
* 28 SSI frequencies, 0.1-20 Hz, dense where the frequencies fall as the walls crack
FREQ,1,4,20,41,61,82,102,123,143,164,184
FREQ,1,205,225,246,266,287,307,328,348,369,389
FREQ,1,410,451,492,532,573,655,737,819
* node(i, j, k) = 25 k + 5 j + i + 1 at x = 3i - 6, y = 3j - 6, z = 2k
N,1,-6,-6,0
N,5,6,-6,0
FILL,1,5
NGEN,4,5,1,5,1,0,3,0
NGEN,4,25,1,25,1,0,0,2
* interior nodes of the mid-storey levels are not used
NDEL,32,34
NDEL,37,39
NDEL,42,44
NDEL,82,84
NDEL,87,89
NDEL,92,94
* mat (10 x concrete), slab concrete, and one material per wall panel (M,11 ... M,18): 4 % damping
M,1,3.0E8,0.2,24.0,0.04,0.04,1
M,2,3.0E7,0.2,24.0,0.04,0.04,1
VAR,PM,11,12,13,14,15,16,17,18
FOREACH,PM,M,@PM[#],3.0E7,0.2,24.0,0.04,0.04,1
GROUP,1,SHELL
GTIT,1,foundation mat
MACT,1
E,1,1,2,7,6
EGEN,3,1,1
EGEN,3,5,1,4
THICK,1,16,1,1.5
* storey-1 walls: groups 2-5 (South, East, North, West), 4 x 2 shells each
GROUP,2,SHELL
GTIT,2,panel 1: storey 1 south wall
MACT,11
E,1,1,2,27,26
EGEN,3,1,1
EGEN,1,25,1,4
THICK,1,8,1,0.3
GROUP,3,SHELL
GTIT,3,panel 2: storey 1 east wall
MACT,12
E,1,5,10,35,30
EGEN,3,5,1
EGEN,1,25,1,4
THICK,1,8,1,0.3
GROUP,4,SHELL
GTIT,4,panel 3: storey 1 north wall
MACT,13
E,1,21,22,47,46
EGEN,3,1,1
EGEN,1,25,1,4
THICK,1,8,1,0.3
GROUP,5,SHELL
GTIT,5,panel 4: storey 1 west wall
MACT,14
E,1,1,6,31,26
EGEN,3,5,1
EGEN,1,25,1,4
THICK,1,8,1,0.3
```

### What this does
The commands come in four groups.

* **Site and frequencies.** `ACTM,1` and `MDL,ex07,../ex07` start a second model in its own
  directory. The soil is 25 sublayers of 0.4 m ($V_s = 400\,\text{m/s}$, 5 % damping) on rock
  ($V_s = 1200\,\text{m/s}$); the frequency set has 28 frequencies up to 20 Hz.
* **Nodes.** `N,1` and `N,5` with `FILL` make the five nodes of the line $y = -6$ at $z = 0$ (3 m
  apart); `NGEN,4,5,1,5,1,0,3,0` copies them four times 3 m apart in $y$ (the 25 nodes of the mat,
  numbered with step 5); `NGEN,4,25,1,25,1,0,0,2` copies that level four times 2 m apart in $z$
  (levels $z$ = 2, 4, 6 and 8 m, nodes 26-125, numbered with step 25). The six `NDEL` remove the nine
  interior nodes of each mid-storey level, $z = 2\,\text{m}$ (32-34, 37-39, 42-44) and
  $z = 6\,\text{m}$ (82-84, 87-89, 92-94): only the walls use those levels (CHECK later lists the 18
  missing numbers as gaps, Warning 1, and AFWRITE writes them with all DOFs fixed: harmless). The
  floor ($z = 4\,\text{m}$) and roof ($z = 8\,\text{m}$) levels keep all 25 nodes for the slabs.
* **Materials.** `FOREACH,PM,M,@PM[#],...` defines materials 11 to 18 in one line: identical concrete
  ($E = 30\,\text{GPa}$, $\nu = 0.2$, 4 % damping), but **one material per wall panel**, because
  Option NON changes $E$ panel by panel. Material 1 (the mat) is ten times stiffer, material 2 is the
  slab concrete.
* **Elements.** Each wall is one SHELL group with its own material, built the same way: `E` defines
  the first 3 m × 2 m shell by its four corner nodes, `EGEN,3,1,1` (or `EGEN,3,5,1` for the walls
  along y) copies it three times along the wall, and `EGEN,1,25,1,4` copies those four elements
  one level (node increment 25) up. Storey 1: group 2 (south, $y = -6$), 3 (east), 4 (north),
  5 (west); `THICK,1,8,1,0.3` gives the 8 elements of the group a thickness of 0.3 m. The mat
  (group 1) is meshed like the mat of lesson 07.

### Why it matters
The panel is the unit of nonlinearity: one wall (or wall segment) of one storey, deforming mostly in
shear, with one backbone curve and one equivalent-linear modulus. The mesh and the groups must be
prepared for it: each panel a group of coplanar shells in a vertical plane, with a material of its
own. Low-rise walls ($\text{height}/\text{length} = 4/12$ here) are shear governed, the case the
Cheng-Mertz shear model is made for.

### Technical basis
A panel uses only its four corner nodes: the shear strain is

```math
\gamma = \frac{1}{2}\,\frac{(u_{TL} - u_{BL}) + (u_{TR} - u_{BR})}{H} + \frac{1}{2}\,\frac{(w_{BR} - w_{BL}) + (w_{TR} - w_{TL})}{L}
```

with $u$ and $w$ the in-plane horizontal and the vertical displacement of the corners BL, BR, TR and
TL (bottom/top, left/right) and $L$, $H$ the width and the height of the panel; it is zero for rigid
translations and in-plane rotations
([Theory §17](docs/theory/THEORY_MANUAL.md#17-option-non-nonlinear-structures-by-equivalent-linearisation),
[OPTION_NON.md §3.3](docs/user/OPTION_NON.md#33-wall-panels-p-pnlgen-plist)).

```action
plot-model
```

## Option NON, part 2: storey 2, slabs, masses and the panel data

The storey-2 walls (groups 6-9), the floor and roof slabs (group 10) with 400 t of equipment each,
then the Option NON data: panels, capacities, backbone curves and the equivalent-linear options.

```sassi
GROUP,6,SHELL
GTIT,6,panel 5: storey 2 south wall
MACT,15
E,1,51,52,77,76
EGEN,3,1,1
EGEN,1,25,1,4
THICK,1,8,1,0.3
GROUP,7,SHELL
GTIT,7,panel 6: storey 2 east wall
MACT,16
E,1,55,60,85,80
EGEN,3,5,1
EGEN,1,25,1,4
THICK,1,8,1,0.3
GROUP,8,SHELL
GTIT,8,panel 7: storey 2 north wall
MACT,17
E,1,71,72,97,96
EGEN,3,1,1
EGEN,1,25,1,4
THICK,1,8,1,0.3
GROUP,9,SHELL
GTIT,9,panel 8: storey 2 west wall
MACT,18
E,1,51,56,81,76
EGEN,3,5,1
EGEN,1,25,1,4
THICK,1,8,1,0.3
GROUP,10,SHELL
GTIT,10,floor and roof slabs
MACT,2
E,1,51,52,57,56
EGEN,3,1,1
EGEN,3,5,1,4
EGEN,1,50,1,16
THICK,1,32,1,0.4
FIXROT
* 400 t on each slab: 16 t (157 kN) at each of its 25 nodes
MT,51,157,157,157
MTGEN,24,1,51
MTGEN,1,50,51,75
INT,1,25,1,1
POINT,0,0,2.7
HOUSE,9.81,0,0,2,0,0,0,0,0
* ---- Option NON data
PNLGEN
PLIST
* SHEAR,<panel>,<fc>,<fy>,<rho>,<Nu>: f'c 30 MPa, fy 420 MPa, 0.5 % web steel, no axial force
SHEAR,0,30000,420000,0.005,0
* BBCGEN,<Panel>,<ShearModel>,<fc>,<fy>,<Pn>,<Nu>,<bre>,<bys>,<CrackingForceLevel>
BBCGEN,0,1,30000,420000,0.005,0,0,0,0.3
* EQL,<disp>,<NonLinOpts>,<dampCutoff %>,<dampScale>,<ElasicD>
EQL,0.8,1,0,0,1
NONLINMOTDISP
```

### What this does
* **The rest of the structure.** Groups 6-9 are the storey-2 walls, built like those of storey 1 one
  level higher. Group 10 holds the floor ($z = 4\,\text{m}$) and roof ($z = 8\,\text{m}$) slabs: the 3
  × 3 m floor mesh is made like the mat, and `EGEN,1,50,1,16` copies its 16 elements 50 node numbers
  up to the roof (32 elements, 0.4 m thick). `FIXROT` fixes the drilling rotation where only coplanar
  shells meet (51 nodes inside the mat, the slabs and the walls), since a flat shell has no drilling
  stiffness. `MT,51,...` with the two `MTGEN` puts 157 kN (16 t) on each of the 25 floor nodes and
  copies them to the 25 roof nodes; the concrete's own weight comes from the materials. `INT`, `POINT`
  and `HOUSE` are those of a surface mat (lesson 07).
* `PNLGEN` makes one panel per vertical shell group (panels 1-8 = groups 2-9; displacement option 1
  = shear strain from the corners, force option 1 = Cheng-Mertz shear hysteresis); `PLIST` lists
  them with their corners, $L$, $H$, $t$ and material.
* `SHEAR,0,30000,420000,0.005,0` prints the shear capacity of every panel by ACI 318-08, Wood 1990,
  Barda 1977 and Gulec-Whittaker 2009 (kN/m², kN).
* `BBCGEN,0,1,...,0.3` builds a 22-point backbone curve for every panel from the ACI 318-08 capacity
  $V_u$: cracking at $0.3\,V_u$ (point 1, on the elastic slope $G A_W$), a smooth rise to $V_u$ at
  0.4 % shear strain, failure at 2 %.
* `EQL,0.8,1,0,0,1`: equivalent-linear displacement factor 0.8, panels only, no damping cut-off,
  damping scale 1, elastic damping (4 %) added to the hysteretic damping. `NONLINMOTDISP` adds the
  panel corners to the MOTION and RELDISP requests.

### Why it matters
`SHEAR` prints, for each of the identical walls ($h_W = 4\,\text{m}$, $l_W = 12\,\text{m}$,
$t_W = 0.3\,\text{m}$): ACI 318-08 **12 472 kN**; Wood 1990 1890 kN from its formula, raised to its
lower bound $6\sqrt{f'_c}\,A_W =$ **9824 kN**; Barda 1977 **11 576 kN**; Gulec-Whittaker 2009
**7527 kN**. The capacities differ by a factor of 1.7: the choice of the shear-strength model weighs
as much as the analysis. `BBCGEN` uses ACI: $V_u = 12\,472\,\text{kN}$, cracking at
$V_\text{cr} = 3742\,\text{kN}$ and $\gamma_\text{cr} = 8.31 \times 10^{-5}$ (on the elastic slope
$G A_W = 4.5 \times 10^{7}\,\text{kN}$), yield at $\gamma = 0.4\,\%$, failure at
$(0.02,\ 12\,721\,\text{kN})$. `NONLINMOTDISP` adds the 12 corner nodes (36 MOTION
requests, 12 RELDISP nodes).

### Technical basis
```math
\begin{aligned}
x_\text{eq} &= \mathrm{EDF}\cdot\max_t\lvert\gamma(t)\rvert\\
E_\text{new} &= E_\text{el}\,\frac{K_\text{sec}(x_\text{eq})}{K_\text{el}}, \qquad K_\text{sec} = \frac{F_\text{bb}(x_\text{eq})}{x_\text{eq}}\\
\xi &= \min\left(\text{cutoff},\ \text{scale}\cdot\xi_h + \xi_\text{el}\right), \qquad \xi_h = \frac{E_D}{4\pi E_S}
\end{aligned}
```
with $\mathrm{EDF} = 0.8$ here, $F_\text{bb}$ the backbone curve (BBC), $K_\text{el}$ its first slope,
$E_\text{el}$ and $\xi_\text{el}$ the elastic modulus and damping of the panel, and $\xi_h$ the
hysteretic damping of the stabilised Cheng-Mertz loop at $x_\text{eq}$ ($E_D$ the energy dissipated
per cycle, $E_S$ the strain energy); cutoff and scale are `<dampCutoff %>` and `<dampScale>` of `EQL`.
The cracking strength $0.3\,V_u$ is a choice of this example (`CrackingForceLevel` 0.10-0.50; 0 uses
the ASCE 4-17 cracking stress $3\sqrt{f'_c}$). On the choice of the capacity model the ACS SASSI
manual says that Wood and Gulec-Whittaker come close to the median strengths of squat-wall tests
(Gulec-Whittaker being sensitive to the panel aspect ratio), and that ACI 318-08 and, for typical
nuclear shear walls, Barda can overestimate the strength. This example uses ACI; a design would
justify its choice
([OPTION_NON.md §5](docs/user/OPTION_NON.md#5-shear-capacities-and-backbone-generation-shear-bbcgen);
SHEAR and BBCGEN verified by [VP-46](docs/verification/VERIFICATION_MANUAL.md#vp-46)).

```action
plot-model
explain: BBCGEN,0,1,30000,420000,0.005,0,0,0,0.3
open-doc: docs/user/OPTION_NON.md#5-shear-capacities-and-backbone-generation-shear-bbcgen
```

## Option NON, part 3: the elastic SSI analysis and the elastic NONLINEAR run

Analysis 0 is the uncracked building: three SITE runs, POINT, HOUSE and ANALYS with restart files;
then MOTION and RELDISP for each direction (the record scaled to 0.6 g horizontally and 0.4 g
vertically), the sum of the three directions and the first NONLINEAR run.

```sassi
AOPT,0,0,0,1,0,0,0,0,0,0,0,0,0,0
SITE,0,1,0,20,2,1,0,1,4096,1,0,0.005,8192,1
WAVE,2,1,1,1,0
CHECK
AFWRITE
RUNSITE
FCOPY,FILE1,FILE1X
SITE,0,0,0,20,2,1,1,1,4096,1,1,0.005,8192,1
WAVE,2,0
WAVE,4,1,1,1,0
AFWRITE
RUNSITE
FCOPY,FILE1,FILE1Y
SITE,0,0,0,20,2,1,0,1,4096,1,2,0.005,8192,1
WAVE,4,0
WAVE,3,1,1,1,0
AFWRITE
RUNSITE
FCOPY,FILE1,FILE1Z
ANALYS,0,0,0,1,1,0,0,0,0,0,0,1
AOPT,0,0,0,0,1,1,0,0,1,0,0,0,0,0
AFWRITE
RUNPOINT
RUNHOUSE
RUNANALYS
THFILE,../data/rg160h_030g.acc
THTIT,RG 1.60 spectrum-compatible motion
DAMP,0.05
RELD,1,0,0
* writes ex07_NONLINBAT.pre and COMB_XYZ_THD.inp; NONLINRESET starts a new nonlinear analysis
NONLINBAT,1
NONLINRESET
AOPT,0,0,0,0,0,0,0,0,0,0,1,0,1,0
VAR,NLDX,"AOPT,0,0,0,0,0,0,0,0,0,0,1,0,1,0","EDUOPT,TFFILE,FILE8X","MOTION,0,0,0,20,0,0.1,100,301,0,0.6,1,0,0,0,0,1,0,0,1",AFWRITE,RUNMOTION,RUNRELDISP,"NONLINTHD,X"
VAR,NLDY,"EDUOPT,TFFILE,FILE8Y",AFWRITE,RUNMOTION,RUNRELDISP,"NONLINTHD,Y"
VAR,NLDZ,"EDUOPT,TFFILE,FILE8Z","MOTION,0,0,0,20,0,0.1,100,301,0,0.4,1,0,0,0,0,1,0,0,1",AFWRITE,RUNMOTION,RUNRELDISP,"NONLINTHD,Z"
VAR,NLNON,"COMBXYZTHD,COMB_XYZ_THD.inp",RUNNONLINEAR
FOREACH,NLDX,@NLDX[#]
FOREACH,NLDY,@NLDY[#]
FOREACH,NLDZ,@NLDZ[#]
FOREACH,NLNON,@NLNON[#]
NONLINSAVE,elastic
* for comparison: floor spectra of the uncracked building (X input), kept as EL_*
NOUT,0
NOUT,1,1,0,0,0,1,1,13,63,113
AOPT,0,0,0,0,0,0,0,0,0,0,1,0,0,0
EDUOPT,TFFILE,FILE8X
MOTION,0,0,0,20,0,0.1,100,301,0,0.6,1,0,0,0,0,1,0,0,1
AFWRITE
RUNMOTION
VAR,FL,00013TR_X01.RS,00063TR_X01.RS,00113TR_X01.RS
FOREACH,FL,FCOPY,@FL[#],EL_@FL[#]
* restore the panel-corner requests for the iterations
NOUT,0
NONLINMOTDISP
```

### What this does
* The SSI part is that of lesson 07 with `ANALYS ... <save>` 1 (restart files).
* `NONLINBAT,1` writes `ex07_NONLINBAT.pre`, a generic script of the whole X/Y/Z Option NON loop
  (this lesson runs the same commands itself), and `COMB_XYZ_THD.inp`, the list of X/Y/Z history
  triplets that COMBXYZTHD sums. `NONLINRESET` deletes the state files of any earlier Option NON
  analysis, so the next NONLINEAR run is the elastic one.
* The four variables hold the per-direction command lists, as in example 6. `NLDX`: MOTION and
  RELDISP on FILE8X, with `MOTION ... <max>` = 0.6 scaling the record to a 0.6 g peak (0.4 g for Z
  in `NLDZ`). RELDISP (free-field reference) gives the corner displacement histories, and
  `NONLINTHD,X` keeps them as `X_*.THD`.
* `COMBXYZTHD` adds the X, Y and Z histories time step by time step (ACS SASSI's COMB_XYZ_THD);
  `RUNNONLINEAR` without a `PANEL.NON` file is the **elastic run**: it computes the panel strains of
  the uncracked building, the equivalent-linear properties, and writes the HOUSE deck of the next
  analysis (`ex07_new.hou`). `NONLINSAVE,elastic` keeps its results under `_elastic` names.
* The last block runs MOTION once more on the elastic X case for the mat, floor and roof centres
  (nodes 13, 63, 113), keeps the spectra as `EL_*`, and restores the corner requests with `NOUT,0`
  and `NONLINMOTDISP`.

### Why it matters
In the uncracked building the storey-1 walls reach a peak shear strain of 1.51e-4, 1.8 times their
cracking strain ($\mu = 1.81$ in `Panel_elastic.fmu`); the storey-2 walls stay below it
($\mu = 0.94$). The elastic run therefore proposes, for the storey-1 walls, $E/E_\text{el} = 0.72$ and
11.4 % damping (4 % elastic + 7.4 % hysteretic) for the next analysis, while the storey-2 walls keep
their elastic properties. A linear SSI analysis would stop here, with walls that its own results show
to be cracked.

### Technical basis
Panel deformations come from RELDISP's relative displacements (free-field reference, which cancels in
every difference the panel strain forms), summed over X, Y and Z per time step because the walls see
the three components at once. The Cheng-Mertz shear model (transcribed from INRESB-3D-SUP) gives the
loop at $x_\text{eq}$; the state files `PANEL.NON` and `Panel_EQL_Matl_Prop.txt` drive the
elastic/iteration state machine of the manual
([OPTION_NON.md §2](docs/user/OPTION_NON.md#2-the-workflow-manual-fig-12)).

```action
open-listing: NONLINEAR
open-file: ex07/Panel_elastic.fmu
open-file: ex07/Panel_EQL_Matl_Prop_elastic.txt
```

## Option NON, part 4: iterate to the cracked state and read the result

The iterations are New Structure restarts with the HOUSE deck NONLINEAR wrote. `NONLINITER` runs
them until the panel moduli change by less than 2 % and the damping by less than 0.5 %.

```sassi
ANALYS,0,0,1,0,1,0,0,0,0,0,0,1
AOPT,0,0,0,0,0,0,0,0,1,0,0,0,0,0
AFWRITE
VAR,NLHOUSE,"FCOPY,ex07_new.hou,ex07.hou",RUNHOUSE,RUNANALYS
NONLINITER,NLHOUSE+NLDX+NLDY+NLDZ+NLNON,10
* spectra of the converged building: mat, floor and roof centres (nodes 13, 63, 113), X input
EDUOPT,TFFILE,FILE8X
MOTION,0,0,0,20,0,0.1,100,301,0,0.6,1,0,0,0,0,1,0,0,1
NOUT,0
NOUT,1,1,1,0,0,1,1,13,63,113
AOPT,0,0,0,0,0,0,0,0,0,0,1,0,0,0
AFWRITE
RUNMOTION
* done: one more New Structure run (the same FILE8) that deletes the restart files (ANALYSX <delrst>)
ANALYSX,0,1
AOPT,0,0,0,0,0,0,0,0,1,0,0,0,0,0
AFWRITE
RUNANALYS
ANALYSX,0,0
```

### What this does
* `ANALYS,0,0,1,...`: New Structure (`<mode>` 1); the ANALYS deck is written once, here, and AFWRITE
  inside the loop only writes MOTION and RELDISP decks, so it never overwrites the iterated HOUSE
  deck.
* `NONLINITER,NLHOUSE+NLDX+NLDY+NLDZ+NLNON,10`: per pass, copy `ex07_new.hou` to `ex07.hou`, HOUSE,
  ANALYS, MOTION + RELDISP × 3, COMBXYZTHD, NONLINEAR; at most 10 passes. (`FCOPY` warns in every pass
  that it replaces `ex07.hou`: that is the intended step.)
* The next block computes the 5 % spectra of the converged building for the X input.
* The last block is housekeeping. The restart files hold the soil impedance (`COOXqqq`) and the
  factorised system (`COOTKqqq`) of every SSI frequency: 54 MB here, and, for a model with 2,000
  interaction nodes, more than 500 MB per frequency for the impedance alone (lesson 10). Once the
  iterations are over they are not needed. `ANALYSX,0,1` sets *Delete Restart Files* (`<delrst>`), a
  New Structure run without `<save>` reproduces the converged FILE8 and then deletes them, and
  `ANALYSX,0,0` restores the default so that a later restart (the Try this below) keeps its files.

### Why it matters
The iterations converge after **7 passes** (`NONLINEAR_CONVERGENCE.TXT`): $\max\lvert\Delta E/E\rvert$
= 28.0, 22.9, 17.1, 11.5, 6.6, 4.2, 2.4 and 1.3 % (analyses 0-7), while the damping changes fall
faster (7.4 to 0.01 points). The converged state (`Panel_EQL_Matl_Prop.txt`, `Panel.fmu`):

* storey-1 walls: **$E/E_\text{el} = 0.35$**, damping **16.7 %** (4 % + 12.7 % hysteretic), peak
  shear strain $3.69 \times 10^{-4} = 4.4\,\gamma_\text{cr}$ ($\mu = 4.43$), far below the 0.4 % yield
  strain; $F_\mu = 3.34$ (the initial stiffness at that strain would give 16 584 kN, the backbone gives
  4962 kN, 40 % of $V_u$);
* storey-2 walls: $\mu = 0.84$, elastic ($E/E_\text{el} = 1$, 4 %).

Cracking softens the building: the X spectrum of the floor centre peaks at 5.17 g at 6.76 Hz instead
of 5.39 g at 7.94 Hz uncracked, the roof centre at 6.68 g at 6.76 Hz instead of 7.39 g at 7.94 Hz.
Peak accelerations of the converged model: 0.64 g (mat centre), 0.85 g (floor) and 1.12 g (roof). The
16.7 % damping is above the 7 % that ASCE 4 accepts for cracked concrete (OPTION_NON.md §11): a design
analysis would cap it with `<dampCutoff>`.

**When is Option NON warranted?** When the walls are expected to crack under the design motion (high
shear stress relative to the cracking strength), and when the frequency shift of the cracked building
matters for the ISRS or the demand: low-rise shear-wall buildings near the peak of the input spectrum,
beyond-design-basis and margin assessments, base isolators or nonlinear springs (sliding,
pile-soil interfaces). For design-basis analyses the usual alternative is a linear analysis with the
stiffness and damping of the expected response level (the manual quotes the ASCE 4-17 levels:
uncracked with 4 %, or 0.5 × the uncracked shear and bending stiffness with 7 %); Option NON computes
the cracked stiffness and damping panel by panel instead of assuming them.

### Technical basis
Convergence test (D-NON-06):

```math
\max\frac{\lvert E_\text{new} - E_\text{old}\rvert}{E_\text{old}} < 2\,\%, \qquad \max\,\lvert\xi_\text{new} - \xi_\text{old}\rvert < 0.5\,\%
```

In force control (a stiff wall whose shear is set by the inertia), each iteration reduces the distance
to the fixed point only by the factor $1 - x\,F_\text{bb}'/F_\text{bb}$, close to 1 just after
cracking where the backbone is flat: hence slow convergence and the manual's advice to use smooth
backbones and enough SSI frequencies around the falling structural frequencies. The method is verified
on an SDOF with a closed-form fixed point by
[VP-NON1](docs/verification/VERIFICATION_MANUAL.md#vp-non1);
[Theory §17](docs/theory/THEORY_MANUAL.md#17-option-non-nonlinear-structures-by-equivalent-linearisation),
[OPTION_NON.md §8](docs/user/OPTION_NON.md#8-convergence-what-to-expect).

### In ANSYS terms
The converged E of each panel is the cracked (secant) modulus to use in an ANSYS model of the same
building, for example a modal or harmonic analysis with cracked walls; the damping ratios map to
material damping as for the elastic model.

### Try this
Add the damping cut-off of 7 % that ASCE 4 accepts for cracked concrete (`EQL,0.8,1,7,0,1`), then
`NONLINRESET` and rerun parts 3 and 4: OPTION_NON.md §11 reports that the walls then deform more
($E/E_\text{el} = 0.18$, $\mu = 12.5$) and the iterations do not converge in 10 passes.

```action
open-file: ex07/NONLINEAR_CONVERGENCE.TXT
open-file: ex07/Panel_EQL_Matl_Prop.txt
open-file: ex07/Panel.fmu
plot-spectrum: ex07/EL_00063TR_X01.RS, ex07/00063TR_X01.RS, ex07/EL_00113TR_X01.RS, ex07/00113TR_X01.RS | log
plot-history: ex07/Panel0001.thd, ex07/Panel0001_elastic.thd
open-doc: docs/user/OPTION_NON.md#11-tutorial-example-7
```
