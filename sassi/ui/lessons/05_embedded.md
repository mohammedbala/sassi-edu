---
id: 05-embedded
title: Embedded structures and the flexible-volume method
part: Fundamentals
order: 5
minutes: 40
example: ex02_embedded_box
summary: Model an embedded basement with its excavated soil, compare the FV, FI-FSIN, FI-EVBN and FFV interaction sets, and isolate kinematic interaction with a massless basement.
objectives: [Model the excavated soil and its embedment layers, Apply the interaction-node rules and the model checks, Explain the flexible-volume equation for an embedded structure, Compare FV with the reduced interaction sets and recognise the subtraction-method anomaly, Validate a reduced method against FV, Measure kinematic interaction with a massless basement]
prerequisites: [04-surface-ssi]
---
Most nuclear buildings are embedded: the basement walls and the basemat are below grade. In a
fixed-base model you would clamp the basemat and perhaps add soil springs on the walls. In SASSI the
embedment is modelled explicitly: the finite-element model contains the basement **and** the soil it
displaces (the *excavated soil*). The excavated soil is subtracted, because the free field already
contains it, and the soil impedance and the free-field motion act at the interaction nodes: every
node of the excavated volume in the reference method, or a subset of them in the reduced methods.
This is the **flexible-volume** substructuring method:

```math
\text{SSI system} = \text{free field} + \text{structure} - \text{excavated soil}
```

with

* the free field: the layered soil, no excavation;
* the structure: the FE model, basement included;
* the excavated soil: an FE model with the free-field properties.

```figure
substructuring case=embedded
The three parts for this lesson's box, one at a time. The excavated soil is meshed only to be
subtracted; with every one of its nodes an interaction node (FV) the equation is exact, and the
reduced sets of steps 7 to 9 drop some of these nodes.
```

Embedment brings in two effects the surface mat of lessons 1 and 4 did not have: **kinematic
interaction** (the basement cannot follow a free-field motion that varies with depth) and the choice
of **interaction nodes**, which governs both the cost and the accuracy. In this lesson you build
example 2, a 10 m × 10 m × 5 m reinforced-concrete box in 5 m of sand ($V_s = 250\,\text{m/s}$) over
gravel ($V_s = 400\,\text{m/s}$) over rock ($V_s = 800\,\text{m/s}$), run it with the reference FV
interaction set and with three reduced sets, and run a massless basement.

## The site: one layer number per embedment layer

```sassi
MDL,ex02,ex02
TIT,Ex02 - embedded basement 10 x 10 x 5 m, FV interaction set
GRAVITY,9.81
* L 1-5: the five 1 m embedment layers (sand, Vs 250 m/s), one L number each
L,1,1.0,19.0,500,250,0.05,0.05
L,2,1.0,19.0,500,250,0.05,0.05
L,3,1.0,19.0,500,250,0.05,0.05
L,4,1.0,19.0,500,250,0.05,0.05
L,5,1.0,19.0,500,250,0.05,0.05
* L 6: 16 m of dense gravel (Vs 400 m/s) in 1 m sublayers; L 7: the half-space (Vs 800 m/s)
L,6,1.0,20.0,800,400,0.04,0.04
L,7,1.0,21.0,1600,800,0.02,0.02
TOPL,1,2,3,4,5
TOPL,6,6,6,6,6,6,6,6,6,6
TOPL,6,6,6,6,6,6
* 21 frequencies from 0.1 to 20 Hz (df = 0.0244 Hz)
FREQ,1,4,20,41,61,82,102,123,143,164,184
FREQ,1,205,246,287,328,369,410,492,573,655,737
FREQ,1,819
SITE,0,1,0,20,7,1,0,1,4096,1,0,0.005,8192,1
WAVE,2,1,1,1,0
```

```action
plot-layers
```

### What this does
Seven layer types: five identical sand layers (L 1-5), the gravel (L 6) and the rock half-space (L
7). `TOPL,1,2,3,4,5` makes the first five 1 m sublayers the five embedment layers, each with its own
L number; 16 gravel sublayers follow: 21 sublayers on the half-space. The sublayers pass
$V_s/(5h) = 50\,\text{Hz}$ (sand) and 80 Hz (gravel). The frequencies, SITE and WAVE are as in
lessons 2 and 4 (vertically incident SV, control motion in X at the ground surface), except that 5.5
Hz is left out (21 frequencies) and the half-space is layer type 7 (`<hs>` = 7).

### Why it matters
The embedment depth (5 m) **coincides with a layer interface**, and every basement level that will
carry interaction nodes lies on an interface: interaction nodes must lie on layer interfaces. Each
embedment layer has **its own L number** because the excavated soil of that layer must have the
free-field properties of that layer, and the manual's rule is that L numbers may be repeated in
`TOPL` for identical layers but not for embedment layers. After a SOIL analysis (lesson 6) the five
sand layers would no longer be identical, and each excavated group follows its own layer.

### Technical basis
The excavated soil is subtracted from the free field, so it must have exactly the free-field
properties; otherwise the subtraction leaves a spurious "difference" soil. CHECK warns when an
L number is repeated among embedment layers (EDU-28) and HOUSE checks each excavated element
against the layer at its depth (EDU-08)
([User Guide §7.2](docs/user/USER_GUIDE.md#72-layering-rules)).

## Nodes on the layer interfaces

```sassi
* a 5 x 5 grid (2.5 m) on each of the 6 interface levels z = -5, -4, ..., 0, numbered bottom-up:
*   node(i, j, k) = 25 k + 5 j + i + 1   (i, j = 0..4 along x, y; k = 0..5, z = -5 + k)
N,1,-5,-5,-5
N,5,5,-5,-5
FILL,1,5
NGEN,4,5,1,5,1,0,2.5,0
NGEN,5,25,1,25,1,0,0,1
```

```action
plot-nodes
```

### What this does
Nodes 1-5 along x at the bottom level ($z = -5\,\text{m}$), `NGEN` copies the row four times in y
(2.5 m) and the whole 25-node level five times upwards (1 m): 150 nodes on the six interfaces
$z = -5, \ldots, 0$. The numbering runs from the bottom up.

### Why it matters
These nodes are shared by the basement (shells) and the excavated soil (solids). Interaction nodes
of embedded models must be numbered in ascending order from the bottom up (a warning, EDU-21, and an
error for incoherent analyses). The mesh spacing (2.5 m horizontally, 1 m vertically) is the
excavation mesh: its vertical size must pass the cut-off (1 m passes 50 Hz in the sand); the
horizontal size may be somewhat larger than the vertical one, 1.2 to 1.5 and sometimes 2 times
according to the manual, after a sensitivity study. This tutorial mesh is coarser than that (2.5 m
against 1 m) to stay small; 2.5 m still passes $250/(5 \times 2.5) = 20\,\text{Hz}$, the cut-off,
but a production model would justify such a ratio with a sensitivity study.

### Technical basis
The manual's mesh rules for the excavated soil
([User Guide §8.5](docs/user/USER_GUIDE.md#85-mesh-size-rules-for-the-excavation)): vertical element
size $\le V_s/(5\,f_{\text{cut}})$ (the same one-fifth-wavelength rule as the sublayers), horizontal
size close to the vertical one unless a sensitivity study justifies more, and a mesh as uniform as
possible (the central-zone radius $R_0 = 0.9 \times 2.5\,\text{m}$ assumes a uniform 2.5 m mesh).

## The excavated soil

```sassi
* groups 1-5: 16 SOLID elements each, one group per embedment layer (top down);
* ETYPE 2 = excavated soil; MACT = the L number of the layer the group replaces
GROUP,1,SOLID
GTIT,1,excavated soil 0 to -1 m
MACT,1
E,1,101,102,107,106,126,127,132,131
EGEN,3,1,1
EGEN,3,5,1,4
ETYPE,1,16,1,2
GROUP,2,SOLID
GTIT,2,excavated soil -1 to -2 m
MACT,2
E,1,76,77,82,81,101,102,107,106
EGEN,3,1,1
EGEN,3,5,1,4
ETYPE,1,16,1,2
GROUP,3,SOLID
GTIT,3,excavated soil -2 to -3 m
MACT,3
E,1,51,52,57,56,76,77,82,81
EGEN,3,1,1
EGEN,3,5,1,4
ETYPE,1,16,1,2
GROUP,4,SOLID
GTIT,4,excavated soil -3 to -4 m
MACT,4
E,1,26,27,32,31,51,52,57,56
EGEN,3,1,1
EGEN,3,5,1,4
ETYPE,1,16,1,2
GROUP,5,SOLID
GTIT,5,excavated soil -4 to -5 m
MACT,5
E,1,1,2,7,6,26,27,32,31
EGEN,3,1,1
EGEN,3,5,1,4
ETYPE,1,16,1,2
```

```action
plot-model
explain: ETYPE,1,16,1,2
```

### What this does
The 35 commands are one seven-command group repeated for the five embedment layers; only the group
number, the title, the L number and the node numbers of the first element change:

1. `GROUP,k,SOLID` and `GTIT` open group $k$ and name it.
2. `MACT,k` makes L number $k$ the active material. Because the elements will be excavated soil,
   their material index is an **L number**: group $k$ takes the properties of embedment layer $k$
   ($V_s = 250\,\text{m/s}$, 5 %).
3. `E,1,...` defines the first eight-node SOLID element: the bottom face counter-clockwise seen
   from above, then the top face.
4. `EGEN,3,1,1` copies it three times along x (node increment 1), `EGEN,3,5,1,4` copies that row
   three times along y (increment 5): $4 \times 4 = 16$ elements.
5. `ETYPE,1,16,1,2` declares elements 1 to 16 of the group **excavated soil** (type 2).

In all, 80 elements fill the 10 m × 10 m × 5 m excavation.

### Why it matters
The excavated soil is not a structure: it is the soil that the free field already contains and that
the basement replaces. It must coincide with the excavation volume, use the free-field properties
layer by layer, and touch the structure only on the excavation boundary (`EXCSTRCHK`, next step).
One group per embedment layer also lets STRESS recover soil pressures layer by layer.

### Technical basis
HOUSE builds two sets of matrices: $K^*_s$, $M_s$ for the structure and $K^*_e$, $M_e$ for the
excavated soil, the latter with the complex moduli of the L layer ($G = \rho V_s^2$,
$M = \rho V_p^2$) and never with
incompatible modes ([Theory §13.1](docs/theory/THEORY_MANUAL.md#13-element-formulations-house)).
The HOUSE listing checks every excavated element against the layer at its mid-depth and prints its
passing frequency (50 Hz here); VP-H2 checks the excavated mass $\sum \rho V$ and stiffness
([VP-H2](docs/verification/VERIFICATION_MANUAL.md#vp-h2)).

### In ANSYS terms
There is no ANSYS counterpart: in a direct (box-of-soil) ANSYS model the excavation is simply
absent. In SASSI the soil outside the excavation is never meshed, so the excavated volume is meshed
in order to be taken away.

## The basement and its drilling rotations

```sassi
* reinforced concrete: E = 3E7 kPa, nu 0.2, 24 kN/m3, 5 % damping
M,1,3.0E7,0.2,24.0,0.05,0.05,1
GROUP,6,SHELL
GTIT,6,basement
MACT,1
* base slab (z = -5, elements 1-16) and roof slab (z = 0, 17-32)
E,1,1,2,7,6
EGEN,3,1,1
EGEN,3,5,1,4
EGEN,1,125,1,16
* walls y = -5 (33-52), y = +5 (53-72), x = -5 (73-92), x = +5 (93-112)
E,33,1,2,27,26
EGEN,3,1,33
EGEN,4,25,33,36
E,53,21,22,47,46
EGEN,3,1,53
EGEN,4,25,53,56
E,73,1,6,31,26
EGEN,3,5,73
EGEN,4,25,73,76
E,93,5,10,35,30
EGEN,3,5,93
EGEN,4,25,93,96
THICK,1,16,1,1.0
THICK,17,112,1,0.5
* fix the rotation about the face normal at the face-interior nodes (rows of 3 nodes)
VAR,ZS,7,12,17,132,137,142
VAR,ZE,9,14,19,134,139,144
FOREACH,ZS,D,@ZS[#],@ZE[#],1,1,ROTZ
VAR,YS,27,52,77,102,47,72,97,122
VAR,YE,29,54,79,104,49,74,99,124
FOREACH,YS,D,@YS[#],@YE[#],1,1,ROTY
VAR,XS,31,56,81,106,35,60,85,110
VAR,XE,41,66,91,116,45,70,95,120
FOREACH,XS,D,@XS[#],@XE[#],5,1,ROTX
```

```action
plot-model
explain: FOREACH,ZS,D,@ZS[#],@ZE[#],1,1,ROTZ
```

### What this does
112 SHELL elements on the faces of the excavation: a 1.0 m base slab, a 0.5 m roof slab and four
0.5 m walls, sharing the nodes of the excavated soil. The 31 commands form four groups:

1. **Material and group** (`M` ... `MACT`): concrete with 24 kN/m³ and 5 % damping, group 6.
2. **Slabs** (`E,1` ... `EGEN,1,125,1,16`): the 16 base-slab elements at $z = -5$, copied once with
   a node increment of 125 (five levels of 25 nodes) to make the roof slab at $z = 0$.
3. **Walls** (`E,33` ... `EGEN,4,25,93,96`): for each wall one element, `EGEN` along the wall (node
   increment 1 or 5) and then four times upwards (increment 25): $4 \times 5 = 20$ elements per
   wall. `THICK` gives the slabs and walls their thickness.
4. **Drilling rotations** (`VAR`, `FOREACH`, `D`): each `FOREACH` loops over pairs of variables
   (first and last node of a row) and runs `D,<first>,<last>,<inc>,1,<label>`, which fixes the
   rotation about the face normal at the face-interior nodes (66 DOFs in all, HOUSE listing).

`FIXROT` cannot be used here: it acts on nodes connected to shells only, and these nodes also
belong to excavated solids.

### Why it matters
The basement is a real, flexible structure; the manual recommends modelling it with its actual
properties. The drilling rotations need care: a flat shell has no stiffness about its normal, and
the solids that share the nodes have no rotations at all, so an unrestrained drilling DOF would be
singular. Fixing it inside a face is harmless; on a wall-slab edge the other face's bending
restrains it.

### Technical basis
The structure's nodes on the excavation boundary are shared with the excavated soil, which is the
condition of the flexible-volume method without near-field soil: structure and excavation are
connected only at the foundation-soil interface ([User Guide §8.2](docs/user/USER_GUIDE.md#82-interaction-nodes)).

## Interaction nodes (FV) and the model checks

```sassi
* FV: every node of the excavated volume, level by level (first and last node of each level)
VAR,LFIRST,1,26,51,76,101,126
VAR,LLAST,25,50,75,100,125,150
FOREACH,LFIRST,INT,@LFIRST[#],@LLAST[#],1,1
INTCOUNT
* the checks the manual strongly recommends before any production run
EXCSTRCHK
FIXEDINT
HINGED
CALCM
* POINT: 5 embedded layers (point loads at interfaces 1..6), R0 = 0.9 x 2.5 m
POINT,0,5,2.25
```

```action
plot-nodes
plot-soil
open-doc: docs/user/USER_GUIDE.md#83-choosing-the-interaction-node-set-fv-fi-and-ffv
explain: POINT,0,5,2.25
```

### What this does
`INT` flags all 150 nodes of the excavated volume as interaction nodes: the **flexible-volume (FV)**
set, 450 interaction DOFs (`INTCOUNT`: 3.2 MB per dense matrix). The checks:

* `EXCSTRCHK`: no interior node of the excavation is connected to the structure;
* `FIXEDINT`: no interaction node has a fixed translation;
* `HINGED`: no unintended hinge between 6-DOF and 3-DOF elements;
* `CALCM`: the basement weighs 6000 kN (612 t); the excavated soil it replaces 9500 kN (968 t). The
  basement is lighter than the soil it displaces.

`POINT,0,5,2.25` now has `<layer>` = 5: point loads at the six interfaces 1 to 6 (0 to 5 m), because
interaction nodes exist at all of them.

### Why it matters
In FV every excavated node feels the soil impedance and the free-field motion, which is the exact
form of the substructuring equation and the reference against which every reduced method is judged.
Its cost is the dense impedance on all excavated nodes: $(3N_{\text{int}})^2$ complex numbers, and a
solution time that grows with $N_{\text{int}}^3$. For a nuclear island the FV set can reach tens of
thousands of nodes, which is why the reduced sets exist. The four checks catch the modelling errors
that produce wrong results without any error message.

### Technical basis
Manual Eq. 2.1 with the DOF sets $s$ (structure only), $i$ (shared), $w$ (excavation only) and $f$
(interaction):

```math
\begin{bmatrix}
C^s_{ss} & C^s_{si} & 0\\
C^s_{is} & C^s_{ii} - C^e_{ii} + X_{ii} & -C^e_{iw} + X_{iw}\\
0 & -C^e_{wi} + X_{wi} & -C^e_{ww} + X_{ww}
\end{bmatrix}
\begin{bmatrix} U_s\\ U_i\\ U_w \end{bmatrix}
=
\begin{bmatrix} 0\\ X_{ii}\,U'_i + X_{iw}\,U'_w\\ X_{wi}\,U'_i + X_{ww}\,U'_w \end{bmatrix}
```

where $C = K^* - \omega^2 M$ is the dynamic stiffness of the structure ($C^s$) or of the excavated
soil ($C^e$), and $X$ the soil impedance on the interaction DOFs (zero where a DOF is not one).

It is exact when every excavated node is an interaction node (FV). A check of the whole chain: if
the structure were identical to the excavated soil, $U = U'$ everywhere, the free field; VP-16
verifies this zero-SSI identity to 2e-15
([Theory §3](docs/theory/THEORY_MANUAL.md#3-flexible-volume-substructuring),
[VP-16](docs/verification/VERIFICATION_MANUAL.md#vp-16)).

### In ANSYS terms
A direct ANSYS model would mesh a soil island around the basement, several basement widths across
and deep, with absorbing (viscous) boundaries on its sides and base. SASSI meshes no soil outside the
excavation: the layered site is horizontally infinite, and its whole effect enters through the
impedance $X$ on these 150 interaction nodes. **Plot the model in its soil** (`SHOWSOIL`, the soil
button of the 3D Plot toolbar) draws the free-field layers around the basement as such an island, with
the quarter facing you cut away, for orientation only: hover over a layer for its depth and
properties.

## Run the FV analysis

```sassi
HOUSE,9.81,0,0,2,0,0,0,0,0
ANALYS,0,0,0,0,1,0,0,0,0,0,0
MOTION,0,0,0,20,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1
DAMP,0.05
THFILE,../data/rg160h_030g.acc
THTIT,RG 1.60 horizontal spectrum-compatible motion (EQUAKE, seed 11975)
* X at the base-slab centre (13) and the roof centre (138); Z at the roof edges x = -5 / +5
NOUT,1,1,1,0,0,1,1,13,138
NOUT,3,1,0,0,0,0,1,136,140
* STRESS: maxima of the four bottom panels of wall x = -5 (group 6, elements 73-76), histories of MXX and MYY
STRESS,0,0,1,0,1
EOUT,1,1,1,2,2,1,0,0,0,0,0,0,6,73-76
AOPT,0,0,0,1,1,1,0,0,1,0,1,1,0,0
CHECK
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
RUNMOTION
RUNSTRESS
WRITE
* the steady-state motion at 6.0 Hz, 24 frames over one period (compared with FI-FSIN in the next step)
HARMFRAME,FILE8,6.0,HARM6
```

```action
plot-spectrum: ex02/00138TR_X.TFU, ex02/00013TR_X.TFU
plot-spectrum: ex02/00140TR_Z.TFU, ex02/00136TR_Z.TFU
animate: ex02/HARM6 | deformed 0.3 | FV at 6.0 Hz
open-listing: SITE
open-listing: ANALYS
open-listing: STRESS
```

### What this does
The usual chain; the run takes a few seconds. HOUSE assembles 726 equations (excavated soil mass
968 t); ANALYS solves with 150 interaction nodes on six interfaces; the low-frequency check finds all
transfer functions within 0.008 % of 1 at 0.1 Hz. `EOUT` asks for the maxima of the six shell
components (membrane stresses FXX FYY FXY in kPa, moments MXX MYY MXY per unit length) and the
histories of the two bending moments. Results:

* the roof centre (node 138) transfer function in X stays **at or below 1**: 0.85 at 6 Hz, 0.68 at
  10 Hz, 0.49 at 20 Hz; the base-slab centre (node 13) 0.81 at 5 Hz, 0.29 at 12 Hz, 0.45 at 20 Hz;
* the roof edges move vertically in opposite directions, up to 0.24 at 20 Hz: the box **rocks**
  under a horizontal, vertically incident wave;
* the bending moment MYY of the wall panel 74 reaches 20.9 kN m/m (STRESS listing).

`HARMFRAME,FILE8,6.0,HARM6` writes the steady-state motion of every node at 6.006 Hz (the computed
frequency closest to 6 Hz), per unit harmonic control motion, in 24 frames over one period. The
animation draws it 0.3 m per unit of control motion: the box translates almost as a rigid body, the
roof centre at 0.85 and the base-slab centre at 0.74, nearly in phase with the control motion (3°
and 5° ahead of it), and hardly rocks (roof edges 0.03). The drawing includes the excavated soil
inside the box (groups 1-5), whose nodes it shares.

### Why it matters
Unlike the surface stick, this stiff, light box amplifies nothing: it follows the ground at low
frequency and moves **less** than the free-field surface at high frequency. Compare with the free
field: the SITE listing shows the free-field motion at 5 m depth (interface 6) falling to 0.81 at
5 Hz, 0.32 at 10 Hz and 0.10 at 12 Hz, then rising to 0.82 at 20 Hz. The box cannot follow a motion
that changes so much over its height: it averages it, and it rocks. This is **kinematic
interaction** of an embedded foundation; the last step separates it from the inertial effects.

### Technical basis
ANALYS evaluates the free field $U'$ at every interaction node from FILE1 at the node's depth
([Theory §6.4](docs/theory/THEORY_MANUAL.md#64-evaluation-at-the-interaction-nodes)). The
free-field notch at 12 Hz at 5 m depth is the quarter-wavelength resonance of the 5 m of sand,
$V_s/(4H) = 250/20 = 12.5\,\text{Hz}$, one of the screening frequencies quoted for embedded models
([Theory §3.3](docs/theory/THEORY_MANUAL.md#33-method-variants-the-choice-of-the-interaction-set)).

## FI-FSIN: the subtraction method

Now reduce the interaction set to the foundation-soil interface: the lateral faces and the bottom
face of the excavation. This is FI-FSIN (flexible interface, foundation-soil interface nodes), the
"subtraction method" (SM) of ASCE 4.

```sassi
CPMODEL,2
ACTM,2
MDL,ex02fsin,../ex02_fsin
TIT,Ex02 - embedded basement 10 x 10 x 5 m, FI-FSIN interaction set
* reset the interior nodes (i, j = 1..3) of the levels above the bottom: 15 rows of 3 nodes
VAR,RFIRST,32,37,42,57,62,67,82,87,92,107,112,117,132,137,142
VAR,RLAST,34,39,44,59,64,69,84,89,94,109,114,119,134,139,144
FOREACH,RFIRST,INT,@RFIRST[#],@RLAST[#],1,0
* HOUSE <imp> = 2 (FI): informative, the interaction set defines the method
HOUSE,9.81,0,0,2,2,0,0,0,0
* same site, frequencies and embedment: reuse the free field and the point-load solutions
FCOPY,../ex02/FILE1,FILE1
FCOPY,../ex02/FILE3,FILE3
AOPT,0,0,0,1,0,1,0,0,1,0,1,1,0,0
CHECK
AFWRITE
RUNHOUSE
RUNANALYS
RUNMOTION
RUNSTRESS
WRITE
* the steady-state motion at 6.0 Hz, as for FV
HARMFRAME,FILE8,6.0,HARM6
```

```action
plot-nodes
plot-spectrum: ex02/00138TR_X.TFU, ex02_fsin/00138TR_X.TFU
plot-spectrum: ex02/00140TR_Z.TFU, ex02_fsin/00140TR_Z.TFU
plot-history: ex02/SHELL_006_00074_MYY.THS, ex02_fsin/SHELL_006_00074_MYY.THS
animate: ex02/HARM6 | deformed 0.3 | FV at 6.0 Hz
animate: ex02_fsin/HARM6 | deformed 0.3 | FI-FSIN at 6.0 Hz
open-file: ex02_fsin/ex02fsin.err
```

### What this does
`INT,<n1>,<n2>,<inc>,0` (`<set>` = 0) resets the interaction flag of the 45 interior nodes of the
five upper levels: 105 interaction nodes remain. Because HOUSE `<imp>` = 2 declares a non-FV
method, CHECK now warns **EDU-12: Non-FV Method Selected: Validate Against FV** (the `.err` file).
CHECK reads `<imp>` only; it does not inspect the interaction set, so keep `<imp>` consistent with
the set you build. FILE1 and FILE3 are copied because the site, the frequencies and the embedment
depth have not changed. Compared with FV (computed points):

| | FV | FI-FSIN |
|---|---|---|
| roof X transfer function at 6.0 Hz | 0.85 | **1.34** (+58 %) |
| roof-edge vertical transfer function at 6.0 Hz | 0.029 | **0.58** |
| wall panel 74, moment MYY | 20.9 kN m/m | **43.0 kN m/m** |

`HARMFRAME` writes the FI-FSIN motion at 6.006 Hz as it did for FV. Play the two animations (same
scale, 0.3 m per unit of control motion). In FV the box translates as a nearly rigid body. In FI-FSIN
the walls and the base slab, held by the interaction nodes, move at most 1.38 times the control
motion, but the non-interaction nodes vibrate in a mode of their own: the roof slab see-saws about its
centre line x = 0, its interior nodes moving vertically by up to 15.6 (nodes 137 and 139, in opposite
directions), and the soil inside the box sways by up to 4.1 (node 113). This is the spurious
sub-system of the Technical basis below, made visible.

### Why it matters
FI-FSIN creates a resonance near 6 Hz that does not exist: the roof motion is overestimated by
58 %, the rocking twentyfold, and the wall moment is doubled. Away from 6 Hz the method is close
to FV at low frequency and deviates at high frequency. Nothing in the run flags it except the
comparison with FV: no error, no warning beyond the EDU-12 reminder. This is, in miniature, the
subtraction-method problem identified by the US Department of Energy in 2011 (operating-experience
report OE-3 2011-02, Theory §3.3); it is why, as the ACS SASSI manual states (§4.1.2 item 15),
ASCE 4-16 and SRP 3.7.2 require reduced methods to be **validated against FV** before production
runs. Lesson 11 repeats this comparison on a full shear-wall building with a two-level basement.

### Technical basis
In FI-FSIN the rows of the non-interaction excavated DOFs (the $w$ set) keep only $-C^e_{ww}$ (plus
the structure where it shares those nodes): no impedance and no free-field load act there. Here the
non-interaction nodes are the interior of the excavation and the nine roof-slab nodes on its top
face. Held only by the interface nodes, they form a sub-system "roof slab minus excavated soil"
with no counterpart in reality. At its natural frequencies its dynamic stiffness $C^s - C^e$ is
singular, and the whole system is near-singular there: spurious peaks at and above the first such
frequency, here near 6 Hz
([Theory §3.3](docs/theory/THEORY_MANUAL.md#33-method-variants-the-choice-of-the-interaction-set)).
The manual adds that FI-FSIN is expected to agree with FV on stiff soil and rock sites and can
become numerically unstable at higher frequencies, depending on the surrounding soil stiffness and
the excavation.

### Check yourself
Without the FV run, what in the FI-FSIN results alone could have warned you?

Answer: a sharp peak near 6 Hz in transfer functions that should be smooth (a stiff, light box
has no mode of its own near 6 Hz on this site; its FV roof transfer function decreases steadily
from 1 to 0.49), a large vertical motion of the roof edges at that frequency, and wall membrane
stresses far above what the box's light loading suggests (883 kPa for FYY of panel 74 in the
STRESS listing; FV gives 30 kPa). Such symptoms are hints, not proof; only the comparison with FV
is a validation.

## FI-EVBN: the modified subtraction method

Add the top face of the excavation (the ground-surface nodes) back to the interaction set. This is
FI-EVBN (flexible interface, excavated-volume boundary nodes), the "modified subtraction method"
(MSM).

```sassi
CPMODEL,3
ACTM,3
MDL,ex02evbn,../ex02_evbn
TIT,Ex02 - embedded basement 10 x 10 x 5 m, FI-EVBN interaction set
* set the interaction flag again on the interior nodes of the top level (z = 0)
VAR,TFIRST,132,137,142
VAR,TLAST,134,139,144
FOREACH,TFIRST,INT,@TFIRST[#],@TLAST[#],1,1
FCOPY,../ex02/FILE1,FILE1
FCOPY,../ex02/FILE3,FILE3
CHECK
AFWRITE
RUNHOUSE
RUNANALYS
RUNMOTION
RUNSTRESS
WRITE
```

```action
plot-spectrum: ex02/00138TR_X.TFU, ex02_fsin/00138TR_X.TFU, ex02_evbn/00138TR_X.TFU
plot-spectrum: ex02/00013TR_X.TFU, ex02_evbn/00013TR_X.TFU
open-listing: STRESS
```

### What this does
Model 3 is a copy of the FI-FSIN model with the 9 interior nodes of the top face flagged again: 114
interaction nodes, the **excavated-volume boundary nodes**. The spurious peak is gone: the roof X
transfer function follows FV within 2.5 % at every frequency (the largest difference is at 20 Hz,
0.479 against 0.491), the base slab within 1.3 %, and the wall moment is 20.8 kN m/m (FV 20.9).

### Why it matters
Nine more interaction nodes (114 instead of 105, against 150 for FV) remove the anomaly. FI-EVBN is
what the manual describes as typically accurate and reasonably fast for shallowly embedded nuclear
islands, but not guaranteed for deeply embedded ones. The validation you just did, transfer
functions at the nodes shared by structure and excavation, is the comparison the manual describes
for the ASCE 4-16 validation; add the response quantities you design for (here the wall moment).
The manual recommends doing it on the full model, or on a massless-foundation model, rather than
on a quarter model, which can hide the instability.

### Technical basis
The ground-surface interaction nodes carry the free-field surface motion and the impedance at the
top of the excavation; they raise the frequencies of the unphysical sub-system above the range of
interest ([Theory §3.3](docs/theory/THEORY_MANUAL.md#33-method-variants-the-choice-of-the-interaction-set)).
The equations and the code are the same for every method: only the interaction set differs.

## FFV and the interaction sets generated by INTGEN

Writing interaction sets by hand is error-prone. `INTGEN` generates them from the excavation mesh.
Start from a copy of the FV model and generate each set in turn, ending with FFV (fast flexible
volume).

```sassi
ACTM,0
CPMODEL,4
ACTM,4
MDL,ex02ffv,../ex02_ffv
TIT,Ex02 - embedded basement 10 x 10 x 5 m, FFV interaction set
* INTGEN,<type>: 0 clear, 1 FV, 2 FI-EVBN, 3 FI-FSIN, 4 surface, 5 FFV (types add up)
INTGEN,0
INTGEN,1
INTCOUNT
INTGEN,0
INTGEN,3
INTCOUNT
INTGEN,0
INTGEN,2
INTCOUNT
* FFV: FI-EVBN plus every second internal level of the excavation
INTGEN,0
INTGEN,5
INTCOUNT
* HOUSE <imp> = 1 (FFV)
HOUSE,9.81,0,0,2,1,0,0,0,0
FCOPY,../ex02/FILE1,FILE1
FCOPY,../ex02/FILE3,FILE3
AOPT,0,0,0,1,0,1,0,0,1,0,1,1,0,0
CHECK
AFWRITE
RUNHOUSE
RUNANALYS
RUNMOTION
RUNSTRESS
```

```action
plot-nodes
plot-spectrum: ex02/00138TR_X.TFU, ex02_evbn/00138TR_X.TFU, ex02_ffv/00138TR_X.TFU
plot-spectrum: ex02/00140TR_Z.TFU, ex02_evbn/00140TR_Z.TFU, ex02_ffv/00140TR_Z.TFU
open-listing: STRESS
```

### What this does
Each `INTGEN,0` clears the interaction flags, so that the next type starts from an empty set (the
types add up otherwise), and `INTCOUNT` reports the result: **FV 150, FI-FSIN 105, FI-EVBN 114, FFV
132** interaction nodes, the same sets you built by hand. The FFV set is FI-EVBN plus the 9 interior
nodes of every second internal level ($z = -3\,\text{m}$ and $z = -1\,\text{m}$). The last set
generated, FFV, is the one that runs (`HOUSE` `<imp>` = 1 records the method). The FFV run follows
FV within 0.44 % for the roof X transfer function and 1.5 % for the roof-edge vertical one; the wall
moment is 21.0 kN m/m.

```figure
interaction-sets method=ffv
The four sets of this step on the excavation mesh, with the counts INTCOUNT reports. The dense
impedance grows with the square of the number of interaction nodes and the solution time with its
cube.
```

### Why it matters
On this small box ANALYS solves every set in about half a second, so the cost argument is invisible;
the memory and time scale as $N_{\text{int}}^2$ and $N_{\text{int}}^3$, so on a basement with tens
of thousands of excavated nodes the reduced sets decide whether the analysis is feasible. FFV is the
manual's choice for deeply embedded structures such as small modular reactors, where FI-EVBN may
fail. With HOUSE `<imp>` = 1 (FFV) or 2 (FI), CHECK reminds you (EDU-12) that the method must be
validated against FV.

### Technical basis
The internal levels cut the excavated volume into thin, stiff sub-volumes whose frequencies lie far
above the frequencies of interest (Ghiocel 2013;
[Theory §3.3](docs/theory/THEORY_MANUAL.md#33-method-variants-the-choice-of-the-interaction-set)).
VP-50 checks the INTGEN counts on a box mesh
([VP-50](docs/verification/VERIFICATION_MANUAL.md#vp-50)); `INTGEN,5,<skip>` sets the spacing of the
internal levels.

### Check yourself
A colleague proposes FI-FSIN for a 20 m deep reactor building on a soft-soil site to save run time,
and shows that the transfer functions of a quarter model agree with FV. Is that enough?

Answer: no. FI-FSIN is the method most prone to spurious resonances on soft sites, the manual warns
that quarter models can look more stable than the full model, and for deep embedment it recommends
FFV or FV. Validate on the full (or massless-foundation) model, at the nodes shared by structure and
excavation and on the design quantities.

## Kinematic interaction: a massless basement

Finally, isolate kinematic interaction: give the basement zero mass (the excavated soil keeps its
mass, it is still subtracted) and run FV again.

```sassi
ACTM,0
CPMODEL,5
ACTM,5
MDL,ex02ml,../ex02_massless
TIT,Ex02 - massless basement, FV (kinematic interaction)
M,1,3.0E7,0.2,0.0,0.05,0.05,1
FCOPY,../ex02/FILE1,FILE1
FCOPY,../ex02/FILE3,FILE3
AOPT,0,0,0,1,0,1,0,0,1,0,1,0,0,0
AFWRITE
RUNHOUSE
RUNANALYS
RUNMOTION
* the foundation input motion at 12 Hz, 24 frames over one period
HARMFRAME,FILE8,12,HARM12
```

```action
plot-spectrum: ex02_massless/00013TR_X.TFU, ex02/00013TR_X.TFU, ex02_massless/00138TR_X.TFU
plot-spectrum: ex02_massless/00140TR_Z.TFU, ex02/00140TR_Z.TFU
plot-spectrum: ex02_massless/00013TR_X01.RS, ex02/00013TR_X01.RS | log
animate: ex02_massless/HARM12 | deformed 3 | Massless box at 12 Hz
```

### What this does
The basement material gets zero unit weight; everything else is the FV model. The base-slab centre
of the massless box moves 0.78 at 5 Hz, 0.39 at 10 Hz, 0.29 at 12 Hz and 0.46 at 20 Hz; the roof
0.41 at 20 Hz; the roof edges move vertically up to 0.16 (rocking). With its real mass the box
gives 0.81, 0.41, 0.29 and 0.45 at the base slab: here inertia changes little, because the basement
is a light (612 t, less than the 968 t of soil it replaces), very stiff box with no superstructure.
A real building on this basement would add its inertia, and the difference would be larger.

`HARMFRAME,FILE8,12,HARM12` writes the motion of the massless box at 12.01 Hz, where the free field
at the base level (5 m deep) moves only 0.10 while the ground surface moves 1. The animation (3 m per
unit of control motion) shows the foundation input motion: the box neither follows the surface nor
the base level. The base slab moves 0.29 and the roof 0.52, 63° and 34° ahead of the control motion
(top and bottom do not move together), and the roof edges move up and down by ±0.10: translation and
rocking from the free field alone, without any inertia of the box.

### Why it matters
The motion of the massless foundation is the **foundation input motion**: the motion your
structure actually receives from the ground, before its own inertia acts. For this box it is well
below the free-field surface motion above 5 Hz and contains a rocking component that a fixed-base
analysis driven by the surface motion would never see. This is one reason why embedded structures on
soil sites have lower high-frequency ISRS than a fixed-base analysis driven by the surface motion
predicts. Recall lesson 1: for a surface mat under vertically incident waves the same analysis gave
exactly the free-field motion. A massless-foundation model is also what the manual suggests for
validating reduced interaction sets against FV.

```figure
kinematic-inertial case=embedded r=0.2
The split for the embedded box: the massless run of this step is problem 1, and its foundation
input motion (translation and rocking) drives problem 2. $D/\lambda$ is the embedment over the
shear wavelength, $f/(50\,\text{Hz})$ for 5 m of sand at 250 m/s: the larger it is, the more the box
averages a free field that varies over its depth.
```

### Technical basis
With zero structural mass the equation keeps the excavated-soil terms,
$\left[K^*_s - C^e + X\right] U = X\,U'$,
so the result depends on the geometry and stiffness of the basement and on the variation of the free
field with depth. Under vertically incident waves this is an **embedment** effect: the walls and the
base slab average a free field that varies over the 5 m depth, and the unequal motion along the
walls rocks the box. (Averaging over the plan of the base slab, "base-slab averaging", needs
inclined or incoherent waves; with vertical coherent waves every point of a horizontal plane moves
in phase.) ([Theory §3](docs/theory/THEORY_MANUAL.md#3-flexible-volume-substructuring).) The full
SSI response then adds inertial interaction on top of it.

### Check yourself
You have a fixed-base ANSYS model of this basement and the design response spectrum at the ground
surface. What does the massless run tell you about the base input you are using?

Answer: the free-field surface motion overestimates the translational input to the basemat above
about 5 Hz (0.39 instead of 1 at 10 Hz) and ignores the rocking the box picks up from the
depth-varying free field. The foundation input motion (translation and rocking at the basemat from
the massless run) is the consistent input; it accounts for kinematic interaction only, so the
inertial part still needs the soil springs and dashpots (lesson 3) or the full SSI run.

### In ANSYS terms
If you must use a simplified ANSYS model of an embedded building, the input should be the
foundation input motion (translation **and** rocking) at the basemat, not the free-field surface
motion, applied through the foundation springs and dashpots of the embedded foundation; a
massless-foundation SASSI run is how that motion is obtained. A clamped base driven by the
foundation input motion still ignores inertial interaction.
