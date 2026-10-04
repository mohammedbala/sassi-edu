---
id: 03-impedance
title: Foundation impedance: soil springs and dashpots
part: Fundamentals
order: 3
minutes: 30
example: ex03_forced_vibration
summary: Compute the springs and dashpots of a rigid square mat on a half-space, compare them with the closed-form solutions you know, and see radiation damping at work.
objectives: [Explain interaction nodes and the central-zone radius, Run a foundation-vibration analysis with FORCE and ANALYS, Read the global impedance files FOUNSTIF FOUNDASH and FOUNDAMP, Compare static stiffness with Pais-Kausel and judge the mesh, Separate radiation damping from material damping]
prerequisites: [02-free-field]
---
If you have ever put a building on springs in ANSYS, you computed the springs from closed-form
formulas: $K = 8GR/(2-\nu)$ for sliding of a disk, the Pais and Kausel or Gazetas formulas for a
rectangle, perhaps with dashpots $\rho V_s A$ and $\rho V_p A$. SASSI computes the same quantity,
the **impedance** of the soil at the foundation, numerically for any layering and any foundation
shape, as a full frequency-dependent complex matrix: the real part is the spring, the imaginary part
the dashpot.

In this lesson you compute it for a case where the closed forms apply: a rigid, massless 12 m × 12 m
mat (half-width $B = 6\,\text{m}$) on a uniform half-space with $V_s = 200\,\text{m/s}$,
$\nu = 1/3$, $\rho = 2.0\,\text{t/m}^3$ ($G = 80{,}000\,\text{kPa}$) and 2 % material damping. You
load the mat with unit forces and moments (a foundation-vibration analysis, as for a machine
foundation), let ANALYS condense the soil impedance to the six rigid-body degrees of freedom,
compare the static values with Pais and Kausel, look at how springs and dashpots change with
frequency, and refine the mesh of interaction nodes.

## A uniform half-space as thin layers

```sassi
MDL,ex03,ex03
TIT,Ex03 - rigid square surface foundation, forced vibration and global impedance
* 24 sublayers of 1 m (4 B deep) of the same soil as the half-space (L 2)
L,1,1.0,19.62,400,200,0.02,0.02
L,2,1.0,19.62,400,200,0.02,0.02
TOPL,1,1,1,1,1,1,1,1,1,1
TOPL,1,1,1,1,1,1,1,1,1,1
TOPL,1,1,1,1
* dt 0.005 s, NFFT 4096: df = 0.0488 Hz; 16 frequencies 0.2 .. 20 Hz
FREQ,1,4,20,41,61,82,102,123,143,164,184
FREQ,1,205,246,287,328,369,410
* Mode 1 only (eigen-solutions for POINT): no free field, no WAVE
SITE,0,1,0,20,2,0,0,1,2048,1,0,0.005,4096,1
```

```action
plot-layers
explain: SITE,0,1,0,20,2,0,0,1,2048,1,0,0.005,4096,1
```

### What this does
The uniform half-space is modelled as 24 m of 1 m sublayers (four half-widths deep) of the same soil
as the half-space (layer type 2), with 20 generated half-space sublayers below. The sublayers pass
$V_s/(5h) = 40\,\text{Hz}$. The 16 frequencies run from 0.195 to 20.0 Hz. `SITE` runs **Mode 1
only** (`<mode2>` = 0): a vibration analysis needs the eigen-solutions of the layered soil for the
point-load problem, but no seismic free field.

### Why it matters
Impedance functions are usually plotted against the dimensionless frequency $a_0 = \omega B/V_s$.
Here $a_0 = 2\pi f \times 6/200 = 0.188\,f$, so the frequency set covers $a_0 = 0.04$ to 3.8
($a_0 = 1$ at 5.3 Hz). The same mat on a stiffer soil would reach the same $a_0$ only at a higher
frequency: what
matters is the ratio of the foundation size to the wavelength.

### Technical basis
Thin-layer discretisation and half-space simulation as in lesson 2
([Theory §4](docs/theory/THEORY_MANUAL.md#4-the-thin-layer-method-site-mode-1),
[§5](docs/theory/THEORY_MANUAL.md#5-half-space-variable-depth-and-viscous-boundary)). For an
impedance study the half-space simulation is essential: the dashpots of the foundation represent
energy carried away by waves into an unbounded medium, which a rigid base would reflect back.

## The foundation: interaction nodes and a rigid mat

```sassi
* 7 x 7 interaction nodes at 2 m, numbered row by row; node 25 is the centre
N,1,-6,-6,0
N,7,6,-6,0
FILL,1,7
NGEN,6,7,1,7,1,0,2,0
INT,1,49,1,1
* K node of the links, above the centre
N,50,0,0,5
* rigid links: E = 1E5 x G of the soil, massless, 2 m x 2 m section
M,1,8.0E9,0.25,0.0,0.0,0.0,1
R,1,4.0,3.333,3.333,2.25,1.333,1.333
GROUP,1,BEAMS
GTIT,1,rigid links
MACT,1
RACT,1
* one link from the centre (25) to each of the other 48 mat nodes
VAR,RIM,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49
FOREACH,RIM,E,#,25,@RIM[#],50
```

```action
plot-nodes
plot-model
explain: INT,1,49,1,1
```

### What this does
`N`, `FILL` and `NGEN` make a 7 × 7 grid of nodes at 2 m spacing on the ground surface ($z = 0$).
`INT,<n1>,<n2>,<inc>,<set>` sets the **interaction** flag on nodes 1 to 49. `FOREACH` creates one
BEAMS element from the centre node 25 to each of the other 48 nodes (the variable `RIM` lists them,
`#` is the element number); with $E = 8 \times 10^9\,\text{kPa}$, $10^5$ times the soil modulus, and
zero weight, the mat is rigid and massless, and a load at node 25 moves the whole mat.

### Why it matters
**Interaction nodes are where the soil acts on the structure.** Only their translations couple to
the soil, and the soil impedance is a full matrix on them, so every interaction node interacts with
every other one through the ground. The mat itself carries no soil springs: its stiffness comes
entirely from the soil at these 49 points. The spacing of the interaction nodes is a mesh: too
coarse a grid makes the soil too stiff (last step of this lesson).

### Technical basis
The impedance matrix on the interaction nodes is $X_{ff} = F_{ff}^{-1}$, where $F_{ff}$ is the
flexibility of the layered soil: the displacement at node $i$ due to a unit force at node $j$, for
every pair ([Theory §8](docs/theory/THEORY_MANUAL.md#8-flexibility-and-impedance-assembly)). It has
$3 \times 49 = 147$ rows and columns. The rules for interaction nodes (on layer interfaces, at or
below the surface, not fixed, belonging to the soil side) are in
[User Guide §8.2](docs/user/USER_GUIDE.md#82-interaction-nodes).

### In ANSYS terms
The links are a rigid region (`CERIG`) or a spider of very stiff beams. Where you would attach a
`COMBIN14` spring and dashpot at each base node, SASSI attaches the full matrix $X_{ff}$, coupling
all nodes with each other.

## Point-load solutions: POINT and the central zone

```sassi
* POINT,<opmode>,<layer>,<rad>: surface foundation, central-zone radius R0 = 0.9 x 2 m
POINT,0,0,1.8
HOUSE,9.81,0,0,2,0,0,0,0,0
FORCE,0
* SITE, POINT, HOUSE
AOPT,0,0,0,1,1,1,0,0,0,0,0,0,0,0
CHECK
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
```

```action
open-listing: POINT
explain: POINT,0,0,1.8
open-doc: docs/user/USER_GUIDE.md#84-the-central-zone-radius
```

### What this does
`POINT,<opmode>,<layer>,<rad>`: `<layer>` = 0 embedded layers (loads only at the ground surface),
`<rad>` = 1.8 m central-zone radius. `HOUSE,9.81,...` sets gravity (SI units), the ground elevation
0 and a 3D analysis. `FORCE,0` only initialises the FORCE options. The three modules run in a
fraction of a second. The POINT listing prints the self-flexibility on the load axis: at 0.195 Hz
a unit horizontal load moves its own point $4.57 \times 10^{-6}\,\text{m/kN}$, a vertical load
$3.76 \times 10^{-6}\,\text{m/kN}$; the
imaginary parts grow with frequency (energy carried away by waves).

### Why it matters
POINT is where the soil springs are born. It solves the layered soil once per load depth for a unit
point load, and ANALYS later evaluates that solution at the distance between any two interaction
nodes. A true point load would have an infinite displacement under it, so the load is spread over a
small central zone of radius $R_0$. **$R_0$ calibrates the diagonal terms of the flexibility to the
interaction-node mesh**: the manual's rule is $R_0 = 0.9\,h$ for a square mesh of size $h$
($0.85\,h$ for a triangular mesh, $h$ in 2D). For a non-uniform mesh, `RADIUS` gives the average
equivalent radius, and the manual asks for sensitivity runs with the minimum, average and maximum
radius.

### Technical basis
The near field ($r < R_0$) is an axisymmetric finite element core with the thin-layer discretisation
in depth; outside it the displacements are an exact superposition of outgoing Rayleigh and Love
modes (Hankel functions), joined by a transmitting boundary
([Theory §7](docs/theory/THEORY_MANUAL.md#7-point-load-solutions-central-zone-and-transmitting-boundary-point),
[§7.5](docs/theory/THEORY_MANUAL.md#75-the-central-zone-radius)). VP-07 checks the static limit
against Boussinesq and Cerruti, VP-08 the dynamic surface Green functions of a half-space (Wong),
VP-09 the far field against the exact thin-layer series
([VP-07](docs/verification/VERIFICATION_MANUAL.md#vp-07),
[VP-09](docs/verification/VERIFICATION_MANUAL.md#vp-09)).

## Six unit load cases (FORCE)

A foundation-vibration analysis replaces the earthquake by external loads. Here each load case is a
unit force or moment at the centre node; the result is the dynamic **compliance** of the mat.

```sassi
* FORCE only
AOPT,0,0,0,0,0,0,0,1,0,0,0,0,0,0
* case 1: unit force in X at node 25 -> FILE9001
F,25,1,0,0
AFWRITE
RUNFORCE
FCOPY,FILE9,FILE9001
FDEL,25
* case 2: Y
F,25,0,1,0
AFWRITE
RUNFORCE
FCOPY,FILE9,FILE9002
FDEL,25
* case 3: Z
F,25,0,0,1
AFWRITE
RUNFORCE
FCOPY,FILE9,FILE9003
FDEL,25
* case 4: unit moment about X
MM,25,1,0,0
AFWRITE
RUNFORCE
FCOPY,FILE9,FILE9004
MMDEL,25
* case 5: moment about Y
MM,25,0,1,0
AFWRITE
RUNFORCE
FCOPY,FILE9,FILE9005
MMDEL,25
* case 6: moment about Z
MM,25,0,0,1
AFWRITE
RUNFORCE
FCOPY,FILE9,FILE9006
```

```action
open-listing: FORCE
explain: F,25,1,0,0
```

### What this does
The 30 commands are one five-command group repeated for six load cases:

1. `F,25,1,0,0` (or `MM,25,...`) defines the load. `F,<node>,<fx>,<fy>,<fz>` gives force
   **factors** (and optionally arrival times) at a node; `MM` does the same for moments about X, Y
   and Z.
2. `AFWRITE` writes the FORCE deck (`AOPT` enables FORCE only).
3. `RUNFORCE` turns the factors into a load vector at every analysis frequency (FILE9).
4. `FCOPY,FILE9,FILE900k` keeps that load vector as load case $k$.
5. `FDEL,25` or `MMDEL,25` deletes the load, so that the next case starts clean (the last case
   keeps its moment; the last step of this lesson deletes it).

The FORCE listing of the last case shows the load: a unit moment about Z at node 25, the same at
every frequency.

### Why it matters
This is how a machine foundation, a vibration test or any harmonic load is analysed in SASSI. The
factor multiplies a reference load history given later in MOTION; an arrival time delays the load,
which is how moving loads are modelled. Several load cases are solved in one ANALYS run, with one
factorisation per frequency.

### Technical basis
$P_k(\omega) = a_k \exp(-i\omega t_k)$ for a factor $a_k$ and an arrival time $t_k$
([Theory §9.1](docs/theory/THEORY_MANUAL.md#91-per-frequency)). With `<simul>` = 6, ANALYS reads
FILE9001 ... FILE9006 and writes FILE8001 ... FILE8006 (VP-23 checks that simultaneous cases equal
single runs to 1e-12, [VP-23](docs/verification/VERIFICATION_MANUAL.md#vp-23)).

## ANALYS: compliance and the global impedance

```sassi
* ANALYS,<opmode>,<type>,<mode>,<save>,<prnt>,<fopt>,<ang>,<xc>,<yc>,<zc>,<impe>,<simul>
ANALYS,0,1,0,0,0,0,0,0,0,0,2,6
* SITE (layer table for HOUSE), HOUSE, ANALYS
AOPT,0,0,0,1,0,1,0,0,1,0,0,0,0,0
AFWRITE
RUNANALYS
```

```action
open-listing: ANALYS
open-file: ex03/FOUNSTIF
explain: ANALYS,0,1,0,0,0,0,0,0,0,0,2,6
open-dialog: ANALYSIS/ANALYS
```

### What this does
`<type>` = 1 is a foundation-vibration analysis (loads from FILE9 instead of the free field);
`<prnt>` = 0 prints real and imaginary parts; `<xc>,<yc>,<zc>` = (0, 0, 0) is the reference point,
here the mat centre; `<impe>` = 2 requests the **full 6 × 6 global impedance**; `<simul>` = 6 solves
the six load cases. ANALYS writes:

* FILE8001 ... FILE8006: the displacements of every node per unit load, the compliance;
* `FOUNSTIF` ($\operatorname{Re} K_G$, the springs), `FOUNDASH` ($\operatorname{Im} K_G/\omega$, the
  dashpot coefficients), `FOUNDAMP` ($\operatorname{Im} K_G/(2\lvert\operatorname{Re} K_G\rvert)$,
  the equivalent damping ratios) and `FOUNIMPD` ($\lvert K_G\rvert$), one 6 × 6 matrix per frequency
  in rows X Y Z XX YY ZZ. The ANALYS listing prints the same matrices.

### Why it matters
These files are the soil springs and dashpots of the foundation, ready to compare with hand
calculations, to use in a simplified model, or to check an SSI model. Two independent routes give
the same impedance: inverting the 6 × 6 compliance of node 25 from FILE8001-8006 reproduces $K_G$ to
about 0.2 % at every frequency (largest difference 0.14-0.21 % of the largest term), the difference
being the finite stiffness of the "rigid" links
([examples/README.md](examples/README.md#example-3-forced-vibration-compliance-and-impedance)).

### Technical basis
The global impedance condenses the soil impedance on the interaction nodes to the six rigid-body
motions of the foundation about the reference point
([Theory §8.3](docs/theory/THEORY_MANUAL.md#83-global-unconstrained-foundation-impedance)):

```math
\begin{aligned}
K_G(\omega) &= T^{\mathsf{T}} X_{ff}\, T, \qquad X_{ff} = F_{ff}^{-1}\\[4pt]
T_j &= \begin{bmatrix}
1 & 0 & 0 & 0 & d_z & -d_y\\
0 & 1 & 0 & -d_z & 0 & d_x\\
0 & 0 & 1 & d_y & -d_x & 0
\end{bmatrix}
\end{aligned}
```

where $T_j$ is the block of $T$ for interaction node $j$ and $d = (d_x, d_y, d_z)$ the position of
node $j$ relative to the reference point. $X_{ff} = F_{ff}^{-1}$ is manual Eq. 4.3,
$K + i\omega D = (f + ig)^{-1}$.

For a surface foundation whose interaction nodes all move as a rigid body this is the impedance of
the rigid foundation; for an embedded foundation it is an "unconstrained" impedance of the soil only.
VP-10, VP-11, VP-13 and VP-14 compare it with the classical solutions
([VP-14](docs/verification/VERIFICATION_MANUAL.md#vp-14)).

### In ANSYS terms
For one frequency, FOUNSTIF and FOUNDASH are what you would enter in two `MATRIX27` elements
between the mat reference node and a fixed ground node, one with `KEYOPT(3) = 4` (stiffness) and
one with `KEYOPT(3) = 5` (damping), or in six `COMBIN14` springs and dashpots if you drop the
coupling terms.

## Static stiffness: compare with the closed forms

```sassi
* the soil modulus for the hand calculation
LLIST,1,2
```

```action
open-file: ex03/FOUNSTIF
open-listing: ANALYS
```

### What this does
`LLIST` confirms $G = 80{,}000\,\text{kPa}$ ($\rho = 2.0\,\text{t/m}^3$). The first row of FOUNSTIF
(0.195 Hz, $a_0 = 0.04$) is practically the static stiffness. Compared with Pais and Kausel (1988)
for a square, $\nu = 1/3$ ($GB = 4.8 \times 10^5\,\text{kN/m}$,
$GB^3 = 1.728 \times 10^7\,\text{kN\,m}$):

| | SASSI, 2 m mesh | Pais & Kausel | difference | welded BEM (VP-14) | difference |
|---|---|---|---|---|---|
| $K_x = K_y$ (kN/m) | $2.82 \times 10^6$ | $5.52\,GB = 2.65 \times 10^6$ | +6.3 % | $5.586\,GB$ | +5.0 % |
| $K_z$ (kN/m) | $3.58 \times 10^6$ | $7.05\,GB = 3.38 \times 10^6$ | +5.8 % | $7.057\,GB$ | +5.7 % |
| $K_{xx} = K_{yy}$ (kN m/rad) | $1.24 \times 10^8$ | $6.00\,GB^3 = 1.04 \times 10^8$ | +19.8 % | $6.457\,GB^3$ | +11.4 % |
| $K_{zz}$ (kN m/rad) | $1.70 \times 10^8$ | $8.31\,GB^3 = 1.44 \times 10^8$ | +18.1 % | $8.588\,GB^3$ | +14.2 % |

The matrix also has an off-diagonal term $K_{x,yy} = -1.44 \times 10^6\,\text{kN}$: a horizontal
force at the mat centre also rocks the mat.

### Why it matters
This is the check you can always make: the low-frequency impedance of a surface mat must be close
to the closed-form static stiffness. Two lessons from the numbers:

* the coarse mesh (6 interaction-node spacings across the mat) is **too stiff**, a few percent in
  translation, 11-14 % in rocking and torsion. Like a coarse FE mesh, the discretised soil is
  stiffer than the continuum and softens as the mesh is refined (here and in VP-14); the rotations,
  which depend on the concentration of the contact stresses towards the edges, suffer most. The
  last step refines the mesh;
* closed forms are fits with their own errors: the welded-contact boundary-element solution of
  VP-14 shows that the Pais-Kausel rocking coefficient is about 7 % low for welded contact. A
  difference from a formula is not automatically an error of the model.

### Technical basis
The static stiffnesses and their accuracy are collected in
[docs/spec/R2_benchmarks.md](docs/spec/R2_benchmarks.md) section B.2 (Pais and Kausel, Gazetas).
VP-14 runs 16, 24 and 32 cells per side and extrapolates to within 1.1 % ($K_x$), 0.04 % ($K_z$) of
Pais and Kausel and within 0.7 % of the welded BEM in rocking
([VP-14](docs/verification/VERIFICATION_MANUAL.md#vp-14)). The sliding-rocking coupling term of a
welded mat (about $-0.085\,K_x B$ here) is absent from the closed-form formulas, which give only the
six diagonal terms; VP-17 shows its effect on an SSI response
([VP-17](docs/verification/VERIFICATION_MANUAL.md#vp-17)).

### Check yourself
Your ANSYS model of the same mat uses Pais-Kausel springs. Which degree of freedom would you expect
to differ most from a SASSI analysis with a fine interaction-node mesh, and why?

Answer: rocking (and torsion). The Pais-Kausel rocking fit is about 7 % below the rigorous welded
solution, and the sliding-rocking coupling of a welded mat is not in the formulas. Translation
agrees within a few percent.

## Frequency dependence and radiation damping

A spring and a dashpot that change with frequency are hard to picture. Load the mat with a short
vertical pulse and watch it in the time domain, then read the springs and dashpots.

```sassi
* MOTION on the vertical load case (FILE8003), response = displacement
EDUOPT,TFFILE,FILE8003
MOTIONX,0,0
MOTION,0,0,0,5,0,0.1,100,301,1000,0,1,0,0,0,0,1,0,0,1
* the load history: a 5 Hz Ricker wavelet, peak 1 at t = 0.5 s, times the factor 1000 (kN)
THFILE,../data/ricker_5hz.th
THTIT,Ricker wavelet 5 Hz, peak 1 at t = 0.5 s
NOUT,3,1,1,0,0,0,1,25
AOPT,0,0,0,1,0,0,0,0,0,0,1,0,0,0
AFWRITE
RUNMOTION
```

```action
plot-history: ex03/00025TR_Z.ACC
open-file: ex03/FOUNDASH
open-file: ex03/FOUNDAMP
open-listing: MOTION
```

### What this does
`EDUOPT,TFFILE,FILE8003` makes MOTION read the vertical load case; `MOTIONX,0,0` asks for
displacements; `MOTION` takes the load history from `THFILE` with the factor `<mult>` = 1000, so the
mat receives a 1000 kN vertical Ricker pulse; `NOUT` requests the Z history of node 25 (the file is
called `.ACC` whatever the response type). The load acts upwards (+Z) and so does the mat: it
moves **0.205 mm** at 0.525 s, just after the load peak at 0.5 s, and stops almost at once: from
$t = 1\,\text{s}$ on the motion stays below 0.03 % of the peak. The static displacement under
1000 kN would be $1000/K_z = 0.279\,\text{mm}$.

The springs and dashpots across the frequency range (from FOUNSTIF, FOUNDAMP):

| $f$ (Hz) | $a_0$ | $K_z/K_{z,\text{static}}$ | $K_{yy}/K_{yy,\text{static}}$ | damping ratio X | Z | rocking YY |
|---|---|---|---|---|---|---|
| 2.0 | 0.38 | 0.95 | 0.95 | 0.15 | 0.20 | 0.03 |
| 5.0 | 0.94 | 0.88 | 0.82 | 0.37 | 0.52 | 0.12 |
| 8.0 | 1.51 | 0.79 | 0.71 | 0.60 | 0.99 | 0.30 |

### Why it matters
* **Radiation damping** is large for translation (here 15-60 % equivalent damping in sliding and
  20-99 % vertically between 2 and 8 Hz) and small for rocking at low frequency (3 % at 2 Hz, of
  which 2 % is the soil's material damping). This is why the pulse dies out at once, and why a
  rocking-dominated building gets much less damping from the soil than a sliding one.
* The dashpot coefficients (FOUNDASH) become nearly constant above about 4 Hz:
  $5.95\text{--}6.32 \times 10^4\,\text{kN\,s/m}$ in X against the plane-wave value
  $\rho V_s A = 5.76 \times 10^4$, $1.05\text{--}1.24 \times 10^5$ in Z against
  $\rho V_p A = 1.15 \times 10^5$. These are the dashpots of the classical lumped models.
* The springs soften with frequency (rocking by 29 % at 8 Hz). A single constant spring is a
  compromise; SASSI uses the correct value at each frequency.
* The dynamic peak (0.205 mm) is below the static value (0.279 mm): the pulse has its energy around
  5 Hz, where $\lvert K\rvert$ is larger than the static stiffness because of the damping term.

On a layered site (lesson 4) the picture changes: a soft layer over much stiffer soil traps waves.
Below the first resonance of the layer, $V_s/(4H)$ for horizontal motion and about $V_p/(4H)$ for
vertical and rocking motion, little energy radiates, and over a rigid base none at all (VP-13).

### Technical basis
At 8 Hz ($a_0 = 1.51$) the Pais-Kausel dynamic modifiers and radiation damping ratios for a square
are $\alpha_z = 0.89$, $\alpha_{yy} = 0.71$, $\beta_x = 0.55$, $\beta_z = 0.96$, $\beta_{yy} = 0.27$
(formulas of [R2 C.3](docs/spec/R2_benchmarks.md) evaluated at $a_0 = 1.51$), to which the 2 %
material damping is added; SASSI gives 0.79, 0.71, 0.58, 0.97 and 0.28 (FOUNDAMP minus the 2 %). The
radiation damping agrees within about 5 %; the vertical spring softens about 10 % more in SASSI than
in the fit. The high-frequency dashpots $\rho V_s A$ and $\rho V_p A$ follow from plane-wave
radiation (R2 C.1). VP-11 checks the disk against Veletsos-Verbic, VP-13 a disk on a layer over a
rigid base ([VP-11](docs/verification/VERIFICATION_MANUAL.md#vp-11),
[VP-13](docs/verification/VERIFICATION_MANUAL.md#vp-13)). The response is
$u(t) = \operatorname{IFFT}\left[H(f)\,F(f)\right]$ with $H$ the compliance and $F$ the load
spectrum ([Theory §11.1](docs/theory/THEORY_MANUAL.md#111-convolution)).

### In ANSYS terms
FOUNDAMP is the damping ratio $c\,\omega/(2k)$ of a `COMBIN14` spring-dashpot pair at that
frequency. A 60 % "damping ratio" is not a typo: the radiation damping of a foundation in
translation is often far above any structural damping.

```figure
impedance-ellipse f=5
One `COMBIN14` pair at one frequency: under $u = U\sin\omega t$ the force leads the displacement by
$\varphi$, and the loop it traces encloses the energy the dashpot dissipates per cycle. With the
sliding values of this lesson ($c = \rho V_s A$) the damping ratio at 5 Hz is 0.32; FOUNDAMP gives
0.37 there, with the 2 % material damping and the softened spring.
```

### Check yourself
A tall reactor building on this site responds mainly in rocking at about 2 Hz; a squat, heavy
building responds mainly in sliding at the same frequency. Which one gets more help from the soil
damping?

Answer: the squat building. At 2 Hz the sliding damping ratio is 15 %, the rocking one 3 %, of
which 2 % is the material damping of the soil. Radiation damping in rocking only becomes large at
higher $a_0$ (30 % at 8 Hz here).

## The mesh of interaction nodes

The 2 m mesh was too stiff. Halve the spacing and compute the impedance again; only the soil
impedance is needed, so a seismic run with the global impedance option is enough.

```sassi
* a copy of the model, cleared of the links, nodes and loads
CPMODEL,2
ACTM,2
MDL,ex03f,../ex03_fine
TIT,Ex03 - 12 m x 12 m mat, 1 m mesh of interaction nodes
GDEL,1
NDEL,1,50
FDEL,25
MMDEL,25
* 13 x 13 interaction nodes at 1 m, a stiff massless SHELL mat (the structure only closes the model)
N,1,-6,-6,0
N,13,6,-6,0
FILL,1,13
NGEN,12,13,1,13,1,0,1,0
INT,1,169,1,1
GROUP,1,SHELL
GTIT,1,rigid mat
MACT,1
E,1,1,2,15,14
EGEN,11,1,1
EGEN,11,13,1,12
THICK,1,144,1,2.0
FIXROT
* R0 = 0.9 x 1 m; a free field (Mode 2, vertical SV) so that ANALYS can run a seismic case
POINT,0,0,0.9
SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,4096,1
WAVE,2,1,1,1,0
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
plot-nodes
open-file: ex03_fine/FOUNSTIF
open-listing: ANALYS
```

### What this does
The 32 commands form five groups:

1. **Copy and clear** (`CPMODEL` ... `MMDEL`): copy the model to a new directory; `GDEL` deletes
   the links (group and elements), `NDEL` the 50 nodes, `FDEL`/`MMDEL` the last load.
2. **New interaction nodes** (`N`, `FILL`, `NGEN`, `INT`): a 13 × 13 grid at 1 m replaces the
   7 × 7 grid.
3. **A mat to close the model** (`GROUP` ... `FIXROT`): 144 stiff, massless SHELL elements 2 m
   thick (material 1, `E`/`EGEN` as in lesson 4), with `FIXROT` for the drilling rotations.
4. **Module options**: $R_0$ becomes 0.9 m (`POINT`). Because the global impedance depends only on
   the soil at the interaction nodes, a seismic run (SITE Mode 2 plus a vertical SV wave) with
   `<impe>` = 2 gives it without load cases.
5. **Run** SITE, POINT, HOUSE and ANALYS: 169 interaction nodes (507 interaction DOFs) take
   about 1 s.

The static stiffnesses (first row of `ex03_fine/FOUNSTIF`), compared with the welded BEM of VP-14:

| | 2 m mesh | 1 m mesh | welded BEM |
|---|---|---|---|
| $K_x$ | +5.0 % | +4.6 % | $2.68 \times 10^6\,\text{kN/m}$ |
| $K_z$ | +5.7 % | +5.1 % | $3.39 \times 10^6\,\text{kN/m}$ |
| $K_{xx}$ | +11.4 % | +7.1 % | $1.12 \times 10^8\,\text{kN\,m/rad}$ |
| $K_{zz}$ | +14.2 % | +10.1 % | $1.48 \times 10^8\,\text{kN\,m/rad}$ |

### Why it matters
Halving the interaction-node spacing reduces the rocking and torsion errors by a third; the
translations, already close, change less. What remains comes from the vertical discretisation
(1 m soil sublayers) as well: see Try this. In a real model the interaction-node spacing is set by
the basemat or excavation mesh, so check it the way you check any FE mesh: by refinement, on the
quantities that matter (here the rocking stiffness, which governs the SSI frequency of a tall
building).

### Technical basis
The central-zone radius is tied to the mesh ($R_0 = 0.9\,h$), so refining the mesh also shrinks the
zone over which each point load is spread. The convergence is roughly first order in $h$; VP-14
uses 16 to 32 cells per side and extrapolates
([VP-14](docs/verification/VERIFICATION_MANUAL.md#vp-14)).

### Try this
Also halve the soil sublayers (48 sublayers of 0.5 m to the same 24 m depth) and rerun SITE, POINT,
HOUSE and ANALYS in this model: the errors against the welded BEM drop to +2.6 % ($K_x$), +3.8 %
($K_z$), +4.1 % ($K_{xx}$) and +5.8 % ($K_{zz}$). Both the horizontal (interaction nodes) and the
vertical (sublayers) discretisation stiffen the soil.

```sassi-show
L,1,0.5,19.62,400,200,0.02,0.02
TOPL,0
VAR,HALF,1,2,3,4,5,6,7,8
FOREACH,HALF,TOPL,1,1,1,1,1,1
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
```
