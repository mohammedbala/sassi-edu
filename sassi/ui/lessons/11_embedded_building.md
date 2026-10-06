---
id: 11-embedded-building
title: An embedded building: the subtraction methods against FV
part: Advanced
order: 11
minutes: 50
example: ex08_embedded_building
summary: Build a reinforced-concrete shear-wall building with a two-level basement and a tower in its excavated soil, run it with the FV, FI-FSIN (subtraction) and FI-EVBN (modified subtraction) interaction sets, find where the subtraction method goes wrong and why, and validate the reduced sets against FV.
objectives: [Model a multi-storey embedded building with its excavated soil and an interior structure on nodes of its own, Generate and count the interaction-node sets with INTGEN and weigh their cost, Compare FV with the subtraction and modified subtraction methods on transfer functions ISRS and wall forces, Explain the spurious resonance of the subtraction method from the excavated-soil equations, Document a validation against FV and judge when SM MSM or FFV is acceptable]
prerequisites: [05-embedded]
---
Lesson 5 showed the subtraction-method anomaly on a bare 10 m basement box whose roof slab rested on
excavation nodes. A real nuclear building is different: thick outer walls on the soil, a grid of
interior shear walls and floor slabs inside the basement, heavy equipment, storeys above grade. Does
the subtraction method (SM, FI-FSIN in the ACS SASSI manual) still go wrong on such a building, where,
by how much, and does it matter for the ISRS and the member forces you design with? This lesson
answers with example 8:

* a reinforced-concrete shear-wall building, 24 m × 24 m in plan, with a 2 m basemat 8 m below grade
  and two basement levels; shear walls every 6 m in both directions make sixteen rooms; the main
  block rises two storeys above grade and its four central rooms continue as a tower two storeys
  higher; 810 t of equipment on its floors;
* a site of 8 m of sand and gravel ($V_s = 300\,\text{m/s}$) over dense gravels and rock;
* the same model run with three interaction-node sets: FV (the reference), FI-FSIN (SM) and FI-EVBN
  (the modified subtraction method, MSM), and in the last step FFV.

The analyses run in about one and a half minutes on a laptop, about four in the browser version. The
workspace keeps about 60 MB of results.

```figure
substructuring case=embedded
The three parts of the SSI system. For this building the "structure" is the whole shear-wall
building, basement included, and the "excavated soil" is the 24 m × 24 m × 8 m of soil its basement
replaces; the interaction set decides on which of its nodes the free field acts.
```

## The site, the frequencies and the mesh rule

```sassi
MDL,ex08,ex08
TIT,Ex08 - embedded shear-wall building 24 x 24 m, FV interaction set
GRAVITY,9.81
* L 1-4: the four 2 m embedment layers (sand and gravel, Vs 300 m/s), one L number each
L,1,2.0,19.0,600,300,0.04,0.04
L,2,2.0,19.0,600,300,0.04,0.04
L,3,2.0,19.0,600,300,0.04,0.04
L,4,2.0,19.0,600,300,0.04,0.04
* L 5: dense gravel (1 m sublayers), L 6: very dense gravel (2 m sublayers), L 7: rock half-space
L,5,1.0,20.0,900,450,0.03,0.03
L,6,2.0,21.0,1300,650,0.02,0.02
L,7,1.0,23.0,3000,1500,0.01,0.01
TOPL,1,2,3,4,5,5,5,5,5,5
TOPL,5,5,5,5,5,5,6,6,6,6
TOPL,6,6
* 41 frequencies: 0.1 Hz, then every 0.5 Hz from 0.5 to 20 Hz (df = 0.0244 Hz)
FREQ,1,4,20,41,61,82,102,123,143,164,184
FREQ,1,205,225,246,266,287,307,328,348,369,389
FREQ,1,410,430,451,471,492,512,532,553,573,594
FREQ,1,614,635,655,676,696,717,737,758,778,799
FREQ,1,819
SITE,0,1,0,20,7,1,0,1,4096,1,0,0.005,8192,1
WAVE,2,1,1,1,0
```

```action
plot-layers
```

### What this does
Seven layer types and 22 sublayers on a rock half-space: the four 2 m embedment layers (0 to 8 m,
$V_s = 300\,\text{m/s}$, 4 % damping), 12 m of dense gravel in 1 m sublayers
($V_s = 450\,\text{m/s}$), 12 m of very dense gravel in 2 m sublayers ($V_s = 650\,\text{m/s}$),
then rock ($V_s = 1500\,\text{m/s}$). The 41 SSI frequencies are 0.1 Hz and every 0.5 Hz from
0.5 Hz to the 20 Hz cut-off. SITE and WAVE are those of lessons 4 and 5: a vertically incident SV
wave, control motion in X at the free surface.

### Why it matters
The embedment depth (8 m) falls on a layer interface, and so does every level of the excavation mesh
that will carry interaction nodes. Each embedment layer has its own L number, so that each layer of
excavated soil can follow its own free-field layer (lesson 5); 22 sublayers meet the manual's advice
of more than 20 layers for the surface waves. The stiffness of the embedment soil will also decide,
later in the lesson, at which frequency the subtraction method goes wrong.

### Technical basis
The one-fifth-wavelength rule sets both the sublayers and the excavation mesh: an element of size
$h$ passes waves up to

```math
f_{\max} = \frac{V_s}{5\,h}
```

The excavation will be meshed 3 m in plan and 2 m vertically: $300/(5 \times 3) = 20\,\text{Hz}$
horizontally, $300/(5 \times 2) = 30\,\text{Hz}$ vertically, a plan size 1.5 times the vertical one,
within the ratio the manual accepts after a sensitivity study
([User Guide §8.5](docs/user/USER_GUIDE.md#85-mesh-size-rules-for-the-excavation),
[Theory §4](docs/theory/THEORY_MANUAL.md#4-the-thin-layer-method-site-mode-1)).

### Check yourself
The geotechnical report gives a lower-bound $V_s$ of 220 m/s for the embedment soil. Is the 3 m
mesh still good to 20 Hz for that soil case?

Answer: no: $220/(5 \times 3) = 14.7\,\text{Hz}$. Either refine the plan mesh to 2.2 m (and the
interaction-node count grows by about $(3/2.2)^2 \approx 1.9$) or justify a lower cut-off for that
case. The soil cases of a design (lower, best, upper estimate) can need different meshes.

## The excavated soil and the basement on its boundary

```sassi
* 9 x 9 grid (3 m) on the levels z = -8, -6, -4, -2, 0 (excavation) and z = 5, 10 (main block)
*   node(i, j, k) = 81 k + 9 j + i + 1,  x = -12 + 3 i,  y = -12 + 3 j
N,1,-12,-12,-8
N,9,12,-12,-8
FILL,1,9
NGEN,8,9,1,9,1,0,3,0
NGEN,4,81,1,81,1,0,0,2
NGEN,2,81,325,405,1,0,0,5
* concrete: M 1 walls and basemat; M 2 slabs and roofs (+4 kN/m3 for equipment, piping, live load)
M,1,3.0E7,0.2,24.0,0.04,0.04,1
M,2,3.0E7,0.2,28.0,0.04,0.04,1
* groups 1-4: the excavated soil, one group per embedment layer from the surface down (ETYPE 2)
GROUP,1,SOLID
GTIT,1,excavated soil 0 to -2 m
MACT,1
E,1,244,245,254,253,325,326,335,334
EGEN,7,1,1
EGEN,7,9,1,8
ETYPE,1,64,1,2
GROUP,2,SOLID
GTIT,2,excavated soil -2 to -4 m
MACT,2
E,1,163,164,173,172,244,245,254,253
EGEN,7,1,1
EGEN,7,9,1,8
ETYPE,1,64,1,2
GROUP,3,SOLID
GTIT,3,excavated soil -4 to -6 m
MACT,3
E,1,82,83,92,91,163,164,173,172
EGEN,7,1,1
EGEN,7,9,1,8
ETYPE,1,64,1,2
GROUP,4,SOLID
GTIT,4,excavated soil -6 to -8 m
MACT,4
E,1,1,2,11,10,82,83,92,91
EGEN,7,1,1
EGEN,7,9,1,8
ETYPE,1,64,1,2
* group 5: the basemat (2.0 m) on the bottom face
GROUP,5,SHELL
GTIT,5,basemat
MACT,1
E,1,1,2,11,10
EGEN,7,1,1
EGEN,7,9,1,8
THICK,1,64,1,2.0
* group 12: the outer basement walls (1.0 m) on the lateral faces: south, north, west, east
GROUP,12,SHELL
GTIT,12,outer basement walls
MACT,1
E,1,1,2,83,82
EGEN,7,1,1
EGEN,3,81,1,8
E,33,73,74,155,154
EGEN,7,1,33
EGEN,3,81,33,40
E,65,1,10,91,82
EGEN,7,9,65
EGEN,3,81,65,72
E,97,9,18,99,90
EGEN,7,9,97
EGEN,3,81,97,104
THICK,1,128,1,1.0
```

```action
plot-model
explain: NGEN,4,81,1,81,1,0,0,2
```

### What this does
* **Nodes.** One row of 9 nodes along x, copied 8 times in y and 4 times upwards (2 m): the 405
  nodes of the excavation, numbered bottom-up. The top level (nodes 325-405, $z = 0$) is copied twice
  5 m upwards for the two storeys of the main block (nodes 406-567).
* **Materials.** Concrete with $E = 30\,\text{GPa}$, $\nu = 0.2$, 24 kN/m³ and 4 % damping. The
  slabs and roofs get 4 kN/m³ more: the distributed equipment, piping and live load (2 kPa on a
  0.5 m slab).
* **Excavated soil.** Four groups of 8 × 8 SOLID elements (256 in all), group $k$ in embedment layer
  $k$ and using its L number, declared excavated soil by `ETYPE`.
* **Basemat and outer basement walls.** 64 shells of 2.0 m on the bottom face and 4 × 32 shells of
  1.0 m on the lateral faces, on the same nodes as the excavated soil.

The structural groups are numbered for the colours of the 3D view (colour $g$ of the element
palette for group $g$): purple for the basemat, brown for the walls in the ground; floors, interior
walls, roofs and the walls above grade follow as groups 9, 10, 11 and 13.

### Why it matters
The basemat and the outer walls are where the building touches the soil: the foundation-soil
interface. They share their nodes with the excavated soil, and only there are the two models
connected. A 2 m basemat and 1 m walls are flexible members with their own stiffness and mass, not a
rigid box; SASSI analyses them as they are, and STRESS will give their forces.

### Technical basis
In the notation of the flexible-volume equation the shared nodes are the $i$ set (structure and
excavation), the excavated nodes not shared with the structure the $w$ set
([Theory §3.2](docs/theory/THEORY_MANUAL.md#32-derivation)). HOUSE assembles the excavated soil with
the complex moduli of its L layer, and its mass is $\sum \rho V$
([VP-H2](docs/verification/VERIFICATION_MANUAL.md#vp-h2)).

### In ANSYS terms
The walls and slabs are shells at their mid-surfaces, as SHELL181 elements would be; the basemat
mid-plane sits at the bottom of the excavation. The excavated soil has no counterpart in an ANSYS
model: it is meshed only to be subtracted.

## The interior structure on nodes of its own

Shear walls every 6 m divide the basement into sixteen rooms. These walls and the two slabs inside
the basement stand where the excavated soil is. They must not be connected to it.

```sassi
* separate nodes for the structure inside the excavation: node(i, j, k) + 1000 at the same
* coordinates, interior positions i, j = 1..7 of the levels k = 1..4 (z = -6, -4, -2, 0)
N,1092,-9,-9,-6
N,1098,9,-9,-6
FILL,1092,1098
NGEN,6,9,1092,1098,1,0,3,0
NGEN,3,81,1092,1152,1,0,0,2
* group 10: interior shear walls (0.6 m) on x = -6, 0, 6 and y = -6, 0, 6; in the basement they stand
* on the basemat and end on the outer walls: 1-96 walls x = -6, 0, 6; 97-192 walls y = -6, 0, 6
GROUP,10,SHELL
GTIT,10,interior walls
MACT,1
E,1,3,12,1093,84
E,2,12,21,1102,1093
EGEN,5,9,2
E,8,66,75,156,1147
E,9,84,1093,1174,165
E,10,1093,1102,1183,1174
EGEN,5,9,10
E,16,1147,156,237,1228
EGEN,2,81,9,16
EGEN,2,2,1,32
E,97,19,20,1101,100
E,98,20,21,1102,1101
EGEN,5,1,98
E,104,26,27,108,1107
E,105,100,1101,1182,181
E,106,1101,1102,1183,1182
EGEN,5,1,106
E,112,1107,108,189,1188
EGEN,2,81,105,112
EGEN,2,18,97,128
THICK,1,192,1,0.6
* group 9: floor slabs; here the basement slab z = -4 (0.6 m) and the grade slab z = 0 (0.8 m)
GROUP,9,SHELL
GTIT,9,floor slabs
MACT,2
E,1,163,164,1173,172
E,2,164,165,1174,1173
EGEN,5,1,2
E,8,170,171,180,1179
E,9,172,1173,1182,181
EGEN,5,9,9
E,15,1173,1174,1183,1182
EGEN,5,1,15
EGEN,5,9,15,20
E,51,1179,180,189,1188
EGEN,5,9,51
E,57,226,1227,236,235
E,58,1227,1228,237,236
EGEN,5,1,58
E,64,1233,234,243,242
EGEN,1,162,1,64
THICK,1,64,1,0.6
THICK,65,128,1,0.8
```

```action
plot-model
explain: EGEN,2,18,97,128
```

### What this does
* The separate nodes are the excavation nodes plus 1000, on the 7 × 7 interior positions of the four
  upper levels. The edge nodes of these members are the shared nodes of the outer walls and the
  basemat.
* **Group 10**, the interior basement walls (0.6 m): the wall $x = -6\,\text{m}$ row by row from the
  basemat (the first and last element of a row end on an outer wall, the six between are generated
  with `EGEN`, `EGEN,2,81,...` copies a row twice upwards), then `EGEN,2,2,1,32` copies the whole wall
  to $x = 0$ and $x = 6\,\text{m}$ (node increment 2), and the same for the walls along x
  (`EGEN,2,18,97,128`); 192 shells.
* **Group 9**, the basement slab: corners, edge strips and the 6 × 6 interior block (0.6 m), then
  `EGEN,1,162,1,64` copies it two levels up (node increment 162 for both kinds of node) as the grade
  slab (0.8 m); 128 shells.

### Why it matters
The ACS SASSI manual makes this a basic modelling rule: the excavated soil and the basement share
nodes only on the foundation-soil interface; every other excavated node stays independent of the
structure, even at the same coordinates, so that the two models vibrate independently except where
they really touch. `EXCSTRCHK` checks it, and SASSI-EDU's CHECK treats a structural element on an
interior excavation node as an error. The manual warns that ignoring the rule gave very poor results
in past practice, especially with the reduced methods, and that FV hides it only for stiff walls and
floors. The grade slab follows the rule too: its interior lies on the top face of the excavation, but
that face touches no soil in the real building (it is the open top of the hole), so its interior
nodes are separate as well.

### Technical basis
Rules 4 and 11 of manual §1.5.1 and §9.8.1 (`EXCSTRCHK`); SASSI-EDU's audit measured what sharing
does: a flexible floor slab on interior FV nodes lost its resonance entirely, and an FI-EVBN model
changed by up to 478 % ([Audit report, F-01](docs/verification/AUDIT_REPORT.md)). The checks of
[VP-52](docs/verification/VERIFICATION_MANUAL.md#vp-52) plant such faults and find them.

### The backfill question
In the field the excavation is larger than the building and is backfilled after the walls are cast.
If the backfill is as stiff as the native soil, the model above is the usual idealisation. If it is
softer (or strains more), model it as near-field soil: SOLID elements of the *structure* (`ETYPE` 1)
with the backfill properties, inside a larger excavated volume whose boundary is again the
interaction surface (example 6 and lesson 9 do this for a backfill block). The interaction set is
then built on the larger excavation, and its count grows.

## The superstructure, the tower, the masses and the model checks

```sassi
* the tower nodes (x, y = -6 .. 6) at z = 15 and 20, same numbering formula
N,588,-6,-6,15
N,592,6,-6,15
FILL,588,592
NGEN,4,9,588,592,1,0,3,0
NGEN,1,81,588,628,1,0,0,5
* group 10 continued (MACT is global: set it again): interior walls of the main block (z = 0 to 10)
* and inside the tower (to z = 20)
GROUP,10
MACT,1
E,193,327,1336,417,408
E,194,1336,1345,426,417
EGEN,5,9,194
E,200,1390,399,480,471
E,201,408,417,498,489
EGEN,7,9,201
EGEN,2,2,193,208
E,241,343,1344,425,424
E,242,1344,1345,426,425
EGEN,5,1,242
E,248,1350,351,432,431
E,249,424,425,506,505
EGEN,7,1,249
EGEN,2,18,241,256
E,289,509,518,599,590
EGEN,3,9,289
EGEN,1,81,289,292
E,297,525,526,607,606
EGEN,3,1,297
EGEN,1,81,297,300
THICK,193,304,1,0.6
* group 9 continued: floor z = 5, tower floors z = 10 and 15 (0.5 m), and a stair opening
GROUP,9
MACT,2
E,129,406,407,416,415
EGEN,7,1,129
EGEN,7,9,129,136
E,193,507,508,517,516
EGEN,3,1,193
EGEN,3,9,193,196
EGEN,1,81,193,208
THICK,129,224,1,0.5
EDEL,60
EDEL,124
EDEL,188
ECOMPR
* group 11: the roofs (0.5 m): the main roof at z = 10 around the tower, the tower roof at z = 20
GROUP,11,SHELL
GTIT,11,roofs
MACT,2
E,1,487,488,497,496
EGEN,7,1,1
EGEN,1,9,1,8
E,17,505,506,515,514
EGEN,1,1,17
EGEN,3,9,17,18
E,25,511,512,521,520
EGEN,1,1,25
EGEN,3,9,25,26
E,33,541,542,551,550
EGEN,7,1,33
EGEN,1,9,33,40
E,49,669,670,679,678
EGEN,3,1,49
EGEN,3,9,49,52
THICK,1,64,1,0.5
* group 13: the outer walls above grade: main block (0.8 m) and tower (0.6 m), one shell per storey
GROUP,13,SHELL
GTIT,13,outer walls above grade
MACT,1
E,1,325,326,407,406
EGEN,7,1,1
EGEN,1,81,1,8
E,17,397,398,479,478
EGEN,7,1,17
EGEN,1,81,17,24
E,33,325,334,415,406
EGEN,7,9,33
EGEN,1,81,33,40
E,49,333,342,423,414
EGEN,7,9,49
EGEN,1,81,49,56
E,65,507,508,589,588
EGEN,3,1,65
EGEN,1,81,65,68
E,73,543,544,625,624
EGEN,3,1,73
EGEN,1,81,73,76
E,81,507,516,597,588
EGEN,3,9,81
EGEN,1,81,81,84
E,89,511,520,601,592
EGEN,3,9,89
EGEN,1,81,89,92
THICK,1,64,1,0.8
THICK,65,96,1,0.6
* remove the separate nodes no element uses and close the numbering gaps
RMVUNUSED
NCOM
* drilling rotations: FIXROT for shell-only nodes, D on the face-interior nodes shared with the soil
FIXROT
VAR,BS,11,29,47,65
VAR,BE,17,35,53,71
FOREACH,BS,D,@BS[#],@BE[#],2,1,ROTZ
VAR,YS,83,245,155,317
VAR,YE,89,251,161,323
FOREACH,YS,D,@YS[#],@YE[#],2,1,ROTY
VAR,XS,91,253,99,261
VAR,XE,145,307,153,315
FOREACH,XS,D,@XS[#],@XE[#],18,1,ROTX
* equipment as weights (kN): grade 2 x 150 t, basement slab 2 x 100 t, z = 5 2 x 80 t,
* main roof 2 x 25 t, tower floor z = 15 60 t, tower roof 40 t
MT,733,1471.5,1471.5,1471.5
MT,739,1471.5,1471.5,1471.5
MT,693,981,981,981
MT,699,981,981,981
MT,418,784.8,784.8,784.8
MT,420,784.8,784.8,784.8
MT,551,245.25,245.25,245.25
MT,557,245.25,245.25,245.25
MT,584,588.6,588.6,588.6
MT,611,392.4,392.4,392.4
* the 3D view: hide the excavated soil, so that the building shows whole (display only)
WINDOWSETTINGS,HIDEGROUP,1
WINDOWSETTINGS,HIDEGROUP,2
WINDOWSETTINGS,HIDEGROUP,3
WINDOWSETTINGS,HIDEGROUP,4
CALCM
```

```action
plot-model
plot-soil
explain: NCOM
```

### What this does
* **Superstructure.** The main block has two storeys of 5 m: outer walls of 0.8 m, the interior walls
  of the basement continued (0.6 m), a 0.5 m floor at $z = 5\,\text{m}$ and a 0.5 m roof at
  $z = 10\,\text{m}$. The four central rooms continue as a tower to $z = 20\,\text{m}$: outer walls of
  0.6 m on the lines $x, y = \pm 6\,\text{m}$, one interior wall each way, floors at 10 and 15 m and a
  roof. One shell per storey (3 m × 5 m). The groups: 304 interior walls (group 10), 221 floor slabs
  (group 9), 64 roof shells (group 11), 96 outer walls above grade (group 13).
* **Stair opening.** `EDEL` removes one 3 m × 3 m slab element against the north wall
  ($x = -3 \ldots 0$, $y = 9 \ldots 12\,\text{m}$) in the basement slab, the grade slab and the floor
  at 5 m, and `ECOMPR` renumbers the slab elements without gaps (CHECK treats a missing element number
  as an error).
* **Numbering.** `RMVUNUSED` deletes the 32 separate nodes that no element uses (the room centres of
  $z = -6$ and $-2\,\text{m}$). `NCOM` then numbers the nodes 1 to 781 without gaps: nodes 1-567 keep
  their numbers, the tower nodes become 568-617 (at $z = 20\,\text{m}$ node $593 + 5(j-2) + (i-2)$)
  and the separate basement nodes 618-781 (at $z = -4\,\text{m}$ node $643 + 7j + i$, at grade
  $725 + 7j + i$). The masses below use these numbers.
* **Drilling rotations.** `FIXROT` fixes the rotation about the normal of the nodes connected to
  coplanar shells only, and the rotations of the solid-only nodes. Where shells share nodes with the
  excavated soil (basemat, outer basement walls), `D` in three `FOREACH` loops fixes it at the
  face-interior nodes that no wall or slab restrains (every second node of 4 rows on the basemat and
  of 8 rows on the walls).
* **Masses.** `MT` gives each equipment item as a weight (the default `MUNITS` 1) at the centre of a
  room: 810 t in all.
* **Display.** `WINDOWSETTINGS,HIDEGROUP` hides the four groups of excavated soil in the 3D views:
  the element plot now shows the building, closed, with the brown walls in the ground (the example's
  picture also shows the interaction nodes on them, set in the next step).
* **CALCM**: the elements weigh $1.545 \times 10^5\,\text{kN}$ (15,746 t), the equipment 810 t:
  16,556 t in all. The excavated soil weighs $8.76 \times 10^4\,\text{kN}$ (8,925 t).

### Why it matters
The building is almost twice as heavy as the soil it replaces, and much of that mass is above grade:
inertial interaction will matter, unlike for lesson 5's box, which was lighter than its soil. The
concrete properties are the uncracked stiffness with 4 % damping, the combination the manual's Damping
Cutoff note quotes for ASCE 4 response level 1; a building that cracks under the design earthquake
would use the reduced stiffness and the higher damping of the next level (lesson 9). A fixed-base
eigenvalue analysis of these structural matrices with the basemat clamped puts the first lateral modes
at 12.4 and 12.5 Hz ([examples/README.md](examples/README.md)), in the range of squat nuclear
shear-wall buildings.

### Technical basis
A flat shell element has no stiffness about its normal, and the excavated SOLID elements that share
its nodes carry no rotations: the rotation is singular unless something restrains it
([User Guide §5.5](docs/user/USER_GUIDE.md#55-boundary-conditions-and-unstiffened-rotations)); CHECK
reports any such node (EDU-06). `RMVUNUSED` and `NCOM` remap every node reference (elements, masses,
fixities, output requests) ([VP-51](docs/verification/VERIFICATION_MANUAL.md#vp-51)). Masses given as
weights are divided by the gravity of `HOUSE`.

### In ANSYS terms
The equipment items are MASS21 elements; the extra slab density is a non-structural mass smeared over
the slab. `NCOM` is what `NUMCMP,NODE` does; the hidden groups are what `ESEL,U` before `EPLOT` would
give.

### Try this
Look inside: `WINDOWSETTINGS,VOLUME,,,1` keeps only the elements with a node at $y \ge 1\,\text{m}$
(blank fields keep the model's extent); open a new element plot (`MODELPLOT`) and it shows the
building cut open just north of the wall $y = 0$: the basemat, the two basement levels, the floors and
the interior walls of the rooms, the tower. `WINDOWSETTINGS,VOLUME` without values
restores the whole building.

See the building in its site: **Plot the model in its soil** (`MODELPLOT`, then `SHOWSOIL,1`) draws
the 22 free-field layers and the half-space around the basement, the quarter facing you cut away so
that the basement walls show; the softer layers are lighter. `SHOWSOIL,1,0` closes the cut and
`SHOWSOIL,0` takes the soil away. The picture changes nothing in the analysis: SASSI's layers are
horizontally infinite.

### Check yourself
Why are the basement nodes of the interior structure numbered from 1092 upwards before `NCOM`, and
not given their final numbers directly?

Answer: with "excavation node + 1000" every element of the interior structure can be generated with
the same `EGEN` node increments as the excavation grid (1 along x, 9 along y, 81 or 162 upwards, 2 or
18 from one wall to the next), and the number tells where the node is. `NCOM` afterwards removes the
gaps that CHECK would otherwise list one by one (Warning 1, Gap Found at Node).

## Interaction-node sets with INTGEN: counts, memory and time

```sassi
* show the excavated soil again (the interaction nodes are its nodes), and cut the south and east
* outer walls away so that it shows inside the basement (display only; the copies made later keep it)
WINDOWSETTINGS,SHOWGROUP,1
WINDOWSETTINGS,SHOWGROUP,2
WINDOWSETTINGS,SHOWGROUP,3
WINDOWSETTINGS,SHOWGROUP,4
WINDOWSETTINGS,HIDEELEM,13,1-16
WINDOWSETTINGS,HIDEELEM,13,49-64
WINDOWSETTINGS,HIDEELEM,12,1-32
WINDOWSETTINGS,HIDEELEM,12,97-128
* INTGEN,<type>: 0 clear, 1 FV, 2 FI-EVBN (MSM), 3 FI-FSIN (SM), 5 FFV; the types add up
INTGEN,1
INTCOUNT
INTGEN,0
INTGEN,3
INTCOUNT
INTGEN,0
INTGEN,2
INTCOUNT
INTGEN,0
INTGEN,5
INTCOUNT
* the reference set for the first analysis: FV
INTGEN,0
INTGEN,1
EXCSTRCHK
FIXEDINT
HINGED
```

```action
plot-nodes
plot-model
open-doc: docs/user/USER_GUIDE.md#83-choosing-the-interaction-node-set-fv-fi-and-ffv
```

```figure
interaction-sets method=fsin
The four sets on the small excavation of lesson 5, where they are easy to see. On this building's
9 × 9 × 5 node excavation INTGEN counts FV 405, FI-FSIN 209, FI-EVBN 258 and FFV 307 nodes.
```

### What this does
`WINDOWSETTINGS,SHOWGROUP` brings the excavated soil back into the 3D views and
`WINDOWSETTINGS,HIDEELEM` cuts the south and east outer walls (above and below grade) away: the element
plot shows the four layers of excavated soil filling the basement, the node plot its 405 nodes with the
interaction nodes in red. `INTGEN` builds each set from the excavated-soil elements;
`INTGEN,0` clears the flags so that each type starts from an empty set. `INTCOUNT` reports:

| Set | Interaction nodes | DOFs | One dense matrix $(3N)^2 \times 16\,\text{B}$ |
|---|---|---|---|
| FV (all excavated nodes) | 405 | 1215 | 23.6 MB |
| FI-FSIN, SM (lateral and bottom faces) | 209 | 627 | 6.3 MB |
| FI-EVBN, MSM (FI-FSIN + top face) | 258 | 774 | 9.6 MB |
| FFV (FI-EVBN + the level $z = -4\,\text{m}$) | 307 | 921 | 13.6 MB |

The model ends with the FV set. `EXCSTRCHK` finds no interior excavation node shared with the
structure, `FIXEDINT` no fixed interaction translation, `HINGED` no hinge.

### Why it matters
ANALYS inverts the free-field flexibility on the interaction DOFs and solves with it at every
frequency: memory grows with $N^2$ and time with up to $N^3$ in the interaction count $N$. For a
nuclear island with 30,000 to 50,000 FV interaction nodes, which the manual quotes, this decides
whether the analysis can be run at all. FI-FSIN keeps about half the nodes here, and a smaller
fraction of a larger excavation, whose faces grow more slowly than its volume. On this small model
the dense part is not yet dominant: the 3,930 equations of the structure and the excavation
cost as much, and you will see FV take only about 1.25 times as long as SM.

### Technical basis
The sets and their consequences are the table of
[Theory §3.3](docs/theory/THEORY_MANUAL.md#33-method-variants-the-choice-of-the-interaction-set):
the equations are the same for every set, only $f$ changes. For an $n_x \times n_y \times n_z$ box of
elements FV has $(n_x+1)(n_y+1)(n_z+1)$ nodes, FI-EVBN removes the $(n_x-1)(n_y-1)(n_z-1)$ interior
ones and FI-FSIN also the $(n_x-1)(n_y-1)$ interior nodes of the top face: here 405, 405 − 147 = 258
and 258 − 49 = 209 ([VP-50](docs/verification/VERIFICATION_MANUAL.md#vp-50)). The memory estimate
is that of [User Guide §10.4](docs/user/USER_GUIDE.md#104-memory-and-run-time).

## The reference: FV

```sassi
* POINT: 4 embedded layers (interfaces 1..5), R0 = 0.9 x 3 m
POINT,0,4,2.7
HOUSE,9.81,0,0,2,0,0,0,0,0
ANALYS,0,0,0,0,1,0,0,0,0,0,0
MOTION,0,0,0,20,0,0.1,100,301,1,0,1,0,0,0,0,1,0,0,1
DAMP,0.05
THFILE,../data/rg160h_030g.acc
THTIT,RG 1.60 horizontal spectrum-compatible motion (EQUAKE, seed 11975)
* X: basemat 41, basement slab 675, grade slab 757, floor 446, main roof corner 567, tower roof 605, soil 365
NOUT,1,1,1,0,0,1,1,41,675,757,446,567,605,365
* Z: basemat edge 45, grade edge 369, grade slab under a 150 t item 739, main roof corner 567, tower roof corner 617
NOUT,3,1,1,0,0,1,1,45,369,739,567,617
* wall forces: south wall FXY (elements 4-5), west wall MXX and MYY (68-69), bottom panels (group 12)
STRESS,0,0,1,0,1
EOUT,1,1,2,1,1,1,0,0,0,0,0,0,12,4-5
EOUT,1,1,1,2,2,1,0,0,0,0,0,0,12,68-69
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
* the steady-state motion at 16 Hz, 24 frames over one period (compared with SM in the next step)
HARMFRAME,FILE8,16,HARM16
```

```action
plot-spectrum: ex08/00605TR_X.TFU, ex08/00757TR_X.TFU, ex08/00041TR_X.TFU
plot-spectrum: ex08/00605TR_X01.RS, ex08/00757TR_X01.RS, ex08/00041TR_X01.RS | log
open-listing: ANALYS
open-listing: HOUSE
```

### What this does
The chain of lesson 5 with the outputs a design needs: transfer functions, histories and 5 % ISRS at
the basemat, the basement slab, the grade slab, the floor, the main roof and the tower roof (X), the
vertical response of the basemat and grade edges, of the grade slab under a 150 t item and of the two
roof corners (Z), and the bottom panels of two outer basement walls in STRESS. Node 365 is not
structure: it is the excavated soil at the centre of the top face, inside the basement. CHECK reports
no error and no warning. HOUSE assembles 3,930 equations; ANALYS takes about 26 s for the
41 frequencies (the listing gives the time of each), and finds every transfer function within
0.033 % of 1 at 0.1 Hz.

Results:

* the tower roof transfer function in X peaks at **2.82 at 6.5 Hz**, the corner of the main roof at
  1.71 at 6.5 Hz, the grade slab at 1.15 at 4.5 Hz; the basemat never exceeds 1;
* the 5 % ISRS zero-period accelerations are 0.280 g at the basemat, 0.320 g at grade, 0.364 g at
  the main roof and 0.464 g at the tower roof, for a control motion of 0.324 g; the tower roof ISRS
  peaks at 2.27 g at 6.8 Hz;
* the south basement wall carries an in-plane shear stress FXY of 374 kPa in its bottom panels, the
  west wall a bending moment MYY of 41.6 kN m/m (STRESS listing).

### Why it matters
The building responds at 6.5 Hz on this site, against 12.4-12.5 Hz on a clamped basemat: SSI moves
the floor-spectrum peaks down by about half in frequency. The basemat moves less than the free-field
surface (ZPA 0.280 g against 0.324 g): kinematic interaction of the 8 m embedment (lesson 5). This FV
run is the reference every reduced set must be validated against.

### Technical basis
FV solves the flexible-volume equation exactly in the sense of
[Theory §3](docs/theory/THEORY_MANUAL.md#3-flexible-volume-substructuring): every excavated node
carries the free-field impedance and load; the zero-SSI identity of
[VP-16](docs/verification/VERIFICATION_MANUAL.md#vp-16) verifies the formulation. The low-frequency
check is that of [Theory §9.4](docs/theory/THEORY_MANUAL.md#94-simultaneous-cases-and-the-low-frequency-check).

## The subtraction method: FI-FSIN

```sassi
CPMODEL,2
ACTM,2
MDL,ex08fsin,../ex08_fsin
TIT,Ex08 - embedded shear-wall building 24 x 24 m, FI-FSIN interaction set (SM)
* lateral and bottom faces of the excavation: 209 nodes
INTGEN,0
INTGEN,3
INTCOUNT
* HOUSE <imp> = 2 (FI) records the method
HOUSE,9.81,0,0,2,2,0,0,0,0
* same site, frequencies and embedment: reuse the free field and the point-load solutions
FCOPY,../ex08/FILE1,FILE1
FCOPY,../ex08/FILE3,FILE3
AOPT,0,0,0,1,0,1,0,0,1,0,1,1,0,0
CHECK
AFWRITE
RUNHOUSE
RUNANALYS
RUNMOTION
RUNSTRESS
WRITE
HARMFRAME,FILE8,16,HARM16
```

```action
plot-spectrum: ex08/00605TR_X.TFU, ex08_fsin/00605TR_X.TFU
plot-spectrum: ex08/00041TR_X.TFU, ex08_fsin/00041TR_X.TFU
plot-spectrum: ex08/00365TR_X.TFU, ex08_fsin/00365TR_X.TFU
animate: ex08/HARM16 | deformed 0.3 front | FV at 16 Hz
animate: ex08_fsin/HARM16 | deformed 0.3 front | FI-FSIN at 16 Hz
open-file: ex08_fsin/ex08fsin.err
```

### What this does
The FV model is copied (with the display settings of step 5: soil shown, south and east walls cut
away), its interaction set replaced by the FI-FSIN set (209
nodes), and HOUSE, ANALYS, MOTION and STRESS run again on the free field of the FV run. CHECK warns
EDU-12: a non-FV method must be validated against FV. ANALYS takes about 21 s. Computed transfer
functions, FV against FI-FSIN:

| | 15.0 Hz | 15.5 Hz | 16.0 Hz |
|---|---|---|---|
| tower roof centre X (605) | 0.661 / 0.374 | 0.630 / 0.465 | 0.598 / **1.223** |
| main roof corner X (567) | 0.093 / 0.088 | 0.106 / **0.452** | 0.122 / **0.350** |
| floor $z = 5\,\text{m}$ X (446) | 0.213 / 0.145 | 0.230 / **0.539** | 0.244 / **0.512** |
| grade slab centre X (757) | 0.332 / 0.204 | 0.338 / **0.545** | 0.339 / **0.641** |
| basemat centre X (41) | 0.344 / 0.254 | 0.325 / **0.675** | 0.304 / **0.581** |
| grade slab under 150 t, Z (739) | 0.353 / 0.231 | 0.353 / **0.498** | 0.353 / **0.615** |
| excavated soil, top-face centre X (365) | 0.456 / 3.691 | 1.613 / **12.61** | 1.958 / **11.09** |

Below 10 Hz FI-FSIN follows FV within 2.4 % of the FV peak at every output (complex difference), from
10 to 14 Hz within 6.7 %; then comes a resonance between 15.5 and 16 Hz that FV does not have,
preceded by a dip (the main roof corner moves 0.011 at 14.5 Hz, against 0.095).

`HARMFRAME` writes the steady-state motion at 15.99 Hz of both models. Play the two animations (same
scale, 0.3 m per unit of control motion, seen from the south with the south walls cut away). In FV
the building moves at most 0.62 times the control motion horizontally (the tower roof) and 0.46
vertically (the basement slab under the 100 t items), and the excavated soil inside the basement at
most 2.0. In FI-FSIN the soil inside the basement sways back and forth by up to 11.1 times the control
motion and heaves by up to 2.3, and it drags the building along: the tower roof moves 1.26 and the
basement slab 0.80, twice their FV motion.

### Why it matters
The subtraction method creates a resonance of this building near 15.5 Hz that does not exist, and a
dip just below it, in the range where equipment and distribution systems on the floors of a stiff
building often have their own frequencies. Nothing in the run flags it: no error, only the EDU-12
reminder. Note where it is *not*: below 10 Hz, around the 6.5 Hz SSI peak that governs the building's
overall response, SM is as good as FV. A validation limited to the SSI peak would have passed it.

### Technical basis
In FI-FSIN the 196 excavated nodes inside the basement and on its top face carry neither impedance nor
load: their rows of Eq. (3.2) keep only $-C^e_{ww}$. That block is singular at the natural frequencies
of the enclosed soil with its lateral and bottom faces held, and near the first one the subtracted soil
resonates through the interface ([Theory §3.3](docs/theory/THEORY_MANUAL.md#33-method-variants-the-choice-of-the-interaction-set)).
Unlike lesson 5's box, no structure rests on these nodes here: the anomaly comes from the soil alone.
For a block of plan size $B$ and depth $D$, held on its sides and bottom and free on top, a Rayleigh
estimate of that frequency with the shape $\sin(\pi x/B)\sin(\pi y/B)\sin(\pi z/2D)$ is

```math
f_{EV} \approx \frac{1}{2\pi}\sqrt{\left(V_p^2+V_s^2\right)\left(\frac{\pi}{B}\right)^2+V_s^2\left(\frac{\pi}{2D}\right)^2}
```

an upper bound: $B = 24\,\text{m}$, $D = 8\,\text{m}$, $V_s = 300\,\text{m/s}$,
$V_p = 600\,\text{m/s}$ give 16.8 Hz; the run puts the resonance between 15.5 and 16 Hz. The screening
frequency $V_s/(4D) = 9.4\,\text{Hz}$ of the soil column alone is a lower bound. The anomaly moves with
$V_s$: the stiffer the embedment soil, the higher the frequency, which is why the manual expects
FI-FSIN to agree with FV on stiff soil and rock sites.

### Check yourself
You only have the FI-FSIN results. What in them could make you suspicious?

Answer: a sharp resonance of the basemat and of every floor at the same high frequency, far from the
building's own modes (the basemat moves 0.675 at 15.5 Hz, against 0.254 at 15 Hz and 0.284 at
17 Hz); a dip of the floor motion to almost nothing just below it (0.011 at the main roof corner at
14.5 Hz); and, if you output it, the excavated soil inside the basement moving 12 times the ground.
These are hints. Only the comparison with FV is a validation.

## The modified subtraction method: FI-EVBN

```sassi
CPMODEL,3
ACTM,3
MDL,ex08evbn,../ex08_evbn
TIT,Ex08 - embedded shear-wall building 24 x 24 m, FI-EVBN interaction set (MSM)
* lateral, bottom and top faces of the excavation: 258 nodes
INTGEN,0
INTGEN,2
INTCOUNT
FCOPY,../ex08/FILE1,FILE1
FCOPY,../ex08/FILE3,FILE3
CHECK
AFWRITE
RUNHOUSE
RUNANALYS
RUNMOTION
RUNSTRESS
WRITE
```

```action
plot-spectrum: ex08/00605TR_X.TFU, ex08_fsin/00605TR_X.TFU, ex08_evbn/00605TR_X.TFU
plot-spectrum: ex08/00739TR_Z.TFU, ex08_fsin/00739TR_Z.TFU, ex08_evbn/00739TR_Z.TFU
open-listing: ANALYS
```

### What this does
Model 3 is a copy of the FI-FSIN model with the FI-EVBN set: the 49 interior nodes of the top face
become interaction nodes again, 258 in all. ANALYS takes about 25 s. FI-EVBN follows FV within
3.7 % of the FV peak at every computed frequency and output (the largest difference is the vertical
motion of the grade slab under the 150 t item at 20 Hz); the horizontal transfer functions within
1.3 %.

### Why it matters
Forty-nine more interaction nodes remove the anomaly up to the 20 Hz cut-off, at 64 % of the FV
interaction count. FI-EVBN is the method the manual describes as typically accurate and reasonably
fast, and appropriate for shallowly embedded nuclear islands. It is not a cure: it moves the problem
up in frequency.

### Technical basis
Holding the top face as well raises the natural frequencies of the enclosed soil: in the Rayleigh
estimate $\pi/2D$ becomes $\pi/D$, which gives 23.4 Hz for this block, above the cut-off. The remaining
non-interaction block ($-C^e_{ww}$ of the 147 interior nodes) still has its own frequencies; a softer
site or a deeper excavation brings them back into the range of interest (last step).

## ISRS and basement-wall forces

```sassi
* the 5 % ISRS of the grade slab centre (X) of the three runs, and their differences
READSPEC,../ex08/00757TR_X01.RS,1,1
READSPEC,../ex08_fsin/00757TR_X01.RS,1,2
READSPEC,../ex08_evbn/00757TR_X01.RS,1,3
LINENAME,1,FV
LINENAME,2,FI-FSIN (SM)
LINENAME,3,FI-EVBN (MSM)
SUBTRACTION,4,2,1
SUBTRACTION,5,3,1
LINENAME,4,SM - FV
LINENAME,5,MSM - FV
WRITESPEC,grade_isrs_methods.rs,1,2,3,4,5
```

```action
plot-spectrum: ex08/00757TR_X01.RS, ex08_fsin/00757TR_X01.RS, ex08_evbn/00757TR_X01.RS | log
plot-spectrum: ex08/00739TR_Z01.RS, ex08_fsin/00739TR_Z01.RS, ex08_evbn/00739TR_Z01.RS | log
plot-spectrum: ex08/00605TR_X01.RS, ex08_fsin/00605TR_X01.RS, ex08_evbn/00605TR_X01.RS | log
plot-history: ex08/SHELL_012_00069_MYY.THS, ex08_fsin/SHELL_012_00069_MYY.THS
open-listing: STRESS
```

### What this does
`READSPEC` loads the three grade-slab spectra as lines 1-3, `SUBTRACTION` forms SM − FV and MSM − FV,
and `WRITESPEC` writes the five lines to `ex08_evbn/grade_isrs_methods.rs` for the calculation
package. Over all eleven 5 % ISRS of the outputs:

| | FI-FSIN (SM) against FV | FI-EVBN (MSM) against FV |
|---|---|---|
| the 11 ISRS below 10 Hz | within 4.9 % | within 1.0 % at every frequency |
| the 11 ISRS from 13 to 17 Hz | highs of +3.0 % to +31.2 % (15.5-16.6 Hz), lows down to −8.9 % (13.2-15.1 Hz) | |
| the 11 ISRS above 17 Hz | up to 14.8 % | |
| basemat centre X | +24.9 % at 15.5 Hz | within 0.5 % |
| grade slab centre X | +23.0 % at 15.5 Hz | within 0.2 % |
| grade slab under the 150 t item, Z | +31.2 % at 15.5 Hz | within 1.0 % |
| tower roof centre X | +3.0 % at 16.6 Hz, −8.9 % at 15.1 Hz | within 0.5 % |
| south wall FXY, bottom panels (FV 374 kPa) | 370 kPa (−1.2 %) | 372 kPa (−0.4 %) |
| west wall MYY, bottom panels (FV 41.6 kN m/m) | 42.0 kN m/m (+0.9 %) | 42.1 kN m/m (+1.2 %) |

### Why it matters
Here the spurious resonance reaches the floor spectra: +20 % to +31 % at 15.5 Hz at the basemat, the
grade slab, the floor and the main roof, where an item with its frequency near 15.5 Hz would be
qualified against a demand that does not exist; the tower roof, whose own response at 15.5 Hz is
larger, changes by less than 9 %. The peak wall forces hardly move (about 1 %): they come from the
6.5 Hz SSI mode, far below the anomaly. The size of the ISRS error depends on how much the building
itself moves at the anomaly frequency, so it cannot be judged from another building: the transfer
functions are what the validation compares.

### Technical basis
The ISRS is the response spectrum of the floor history obtained from the interpolated transfer
function and the control motion
([Theory §11](docs/theory/THEORY_MANUAL.md#11-convolution-response-spectra-relative-displacements-and-stresses));
a 5 % oscillator at 15.5 Hz integrates the floor motion over a band of about
$2\beta f = 1.6\,\text{Hz}$, which the spurious resonance fills. STRESS recovers the wall forces from
the stress transfer functions
([Theory §11.5](docs/theory/THEORY_MANUAL.md#115-stresses-and-forces-stress)). Line operations:
[User Guide §11.5](docs/user/USER_GUIDE.md#115-line-mathematics-and-spectrum-broadening).

### Check yourself
The wall forces of the subtraction method are within about 1 % of FV. Can you use SM for the
production runs of this building?

Answer: no. The validation the manual and ASCE 4-16 ask for compares transfer functions, and those
differ by a factor of 2 to 4 near 15.5 Hz; the floor spectra are 20 % to 31 % too high there. Member
forces governed by the low-frequency response do not reveal the anomaly. FI-EVBN passes the same
comparison within 3.7 % and costs little more.

## Validation against FV: what to check and what to document

```sassi
* is every peak of the interpolated transfer functions supported by computed values?
CRITFREQ,5,50,../ex08/00605TR_X,CFV
CRITFREQ,5,50,../ex08_fsin/00605TR_X,CFS
CRITFREQ,5,50,../ex08_fsin/00041TR_X,CFB
CRITFREQ,5,50,00605TR_X,CFE
* the interaction-node count and memory estimate of the active model (FI-EVBN)
INTCOUNT
```

```action
open-file: ex08_fsin/ex08fsin_ANALYS.out
open-file: ex08_evbn/ex08evbn.err
```

### What this does
`CRITFREQ` (lesson 10) compares each peak of an interpolated transfer function (`.TFI`) with the
computed values around it (`.TFU`). The FV tower-roof peak at 6.67 Hz is supported (1.6 %), and so is
the FI-EVBN one (1.5 %). In FI-FSIN it flags a peak at 15.70 Hz on the tower roof (|TFI| 3.80
against computed neighbours of at most 1.22, 211 %) and at 15.72 Hz on the basemat (3.04 against
0.675, 351 %): the 0.5 Hz frequency step does not resolve the spurious resonance, and the
interpolation reconstructs it between 15.5 and 16 Hz.

### Why it matters
A validation is a calculation like any other and must be documented so that a reviewer can repeat
it. For a reduced interaction set the file should contain:

* the method (SM, MSM or FFV), the interaction-node counts of the reduced set and of FV, and why the
  reduced set is needed (memory, run time);
* the validation model: the full SSI model, or a massless-foundation model for the kinematic part as
  the manual suggests; not a quarter model, which the manual reports can hide the instability of an
  SM or MSM model;
* the comparison of the transfer functions over the whole frequency range at the nodes shared by the
  structure and the excavation (here the basemat 41 and 45, the grade edge 369) and at the critical
  locations inside the building (floors, equipment supports), with the largest difference and its
  frequency;
* the soil cases: the softest case governs (the anomaly moves down with $V_s$);
* the frequency set and evidence that it resolves every peak (CRITFREQ, added frequencies), so that
  a spurious peak is not missed between two computed frequencies;
* the acceptance criterion, stated before the comparison, and the design quantities compared (ISRS,
  forces, displacements).

### Technical basis
The validation requirement is stated in the manual's warnings (§1.5.1 item 6 and §4.1.2 item 15):
ASCE 4-16 and SRP 3.7.2 require a validation study against FV before SM, MSM or FFV is used for
production runs; the comparison is of acceleration transfer functions at the common nodes of the
structure and the excavated soil, and the most complete study uses the full SSI model with ATF at
critical locations in the building. The manual also warns that numerical instabilities can occur at
isolated frequencies with any method and recommends inspecting the ATF at several nodes and adding
frequencies near suspect ones. CRITFREQ:
[VP-53](docs/verification/VERIFICATION_MANUAL.md#vp-53); symmetry models:
[VP-T2](docs/verification/VERIFICATION_MANUAL.md#vp-t2).

### Try this
Resolve the anomaly: add frequency number 645 (15.75 Hz) to set 1 of the FI-FSIN model and rerun
SITE, POINT, HOUSE, ANALYS and MOTION there (`ACTM,2`, `FREQ,1,645`, `AOPT,0,0,0,1,1,1,0,0,1,0,1,0,0,0`,
`AFWRITE` and the `RUN` commands; the copied FILE1 and FILE3 do not have the new frequency). The computed
values at 15.75 Hz are 2.10 at the basemat and 3.33 at the tower roof, against 0.314 and 0.614 with
FV: the resonance is real for this model, and the interpolated peaks above (3.04 and 3.80) even
overstate it.

## FFV, and when each method is acceptable

```sassi
ACTM,0
CPMODEL,4
ACTM,4
MDL,ex08ffv,../ex08_ffv
TIT,Ex08 - embedded shear-wall building 24 x 24 m, FFV interaction set
* FI-EVBN plus every second internal level of the excavation: z = -4 m, 307 nodes
INTGEN,0
INTGEN,5
INTCOUNT
* HOUSE <imp> = 1 (FFV); the run itself is the Try this below
HOUSE,9.81,0,0,2,1,0,0,0,0
FCOPY,../ex08/FILE1,FILE1
FCOPY,../ex08/FILE3,FILE3
AOPT,0,0,0,1,0,1,0,0,1,0,1,0,0,0
```

```sassi-show
* the FFV run (about as long as the FV run of this lesson)
CHECK
AFWRITE
RUNHOUSE
RUNANALYS
RUNMOTION
```

```action
plot-nodes
```

### What this does
FFV adds the 49 interior nodes of the level $z = -4\,\text{m}$ to FI-EVBN: 307 interaction nodes
(`INTGEN,5` reports the five node levels and the internal level it keeps). The step prepares the FFV
model; its run is the Try this below. ANALYS takes about as long as for FV (about 30 s against 26 s),
and FFV follows FV within 3.7 % of the FV peak (the vertical motion of the grade slab under the 150 t
item at 20 Hz), the horizontal transfer functions within 1.5 % and the ISRS within 0.8 %.

### Why it matters
On this building FFV buys nothing over FI-EVBN: it is not closer to FV (3.7 % for both, at the
same node and frequency; 1.5 % against 1.3 % horizontally), and it costs as much as FV because the
model is small and its cost is in the structure. Its place is the deep embedment of small modular
reactors, where the manual reports that FI-EVBN and FI-FSIN may fail and where the FV count is
prohibitive. What decides between the methods is the frequency of the enclosed soil against the
cut-off, for the softest soil case. The same building on other sites (rerun the lesson with other
`L,1` to `L,4` lines, Try this below):

| Embedment soil | FI-FSIN against FV | FI-EVBN against FV |
|---|---|---|
| $V_s = 600\,\text{m/s}$ (stiff) | within 1.3 % up to 20 Hz; ISRS within 1.3 % | within 0.9 %; ISRS within 1.1 % |
| $V_s = 300\,\text{m/s}$ (this lesson) | resonance at 15.5-16 Hz, floors 1.9 to 4.3 times FV; ISRS up to +31 % | within 3.7 %; ISRS within 1.0 % |
| $V_s = 200\,\text{m/s}$ (soft) | resonance at 10.5 Hz, basemat 1.31 against 0.57; ISRS +16 % / −16 % | up to 24 % near 19 Hz; ISRS +3.1 % / −9.3 % |

(For $V_s = 200\,\text{m/s}$ the 3 m mesh passes only 13.3 Hz: the comparison of the methods on the
same mesh still holds, the absolute results above 13 Hz do not.)

### Technical basis
The manual: FI-FSIN can become numerically unstable in the higher frequency range depending on the
surrounding soil stiffness and the excavation, and is expected to coincide with FV and FI-EVBN on
stiff soil and rock sites; FI-EVBN is typically accurate and reasonably fast, appropriate for
shallowly embedded nuclear islands, and may fail for deeply embedded models such as SMRs, for which
FFV or FV are appropriate. The resonance frequency scales with $V_s$ (with 0.25 Hz steps the
resonance of this lesson's site is at 15.75 Hz; 15.75 Hz × 200/300 = 10.5 Hz, run: 10.5 Hz), and the
Rayleigh estimate shows how the depth $D$ and
the plan size $B$ enter
([Theory §3.3](docs/theory/THEORY_MANUAL.md#33-method-variants-the-choice-of-the-interaction-set)).

### Try this
Make the embedment soil soft: in the first step use `L,1,2.0,19.0,400,200,0.04,0.04` to
`L,4,2.0,19.0,400,200,0.04,0.04` and rerun the lesson. FI-FSIN now resonates at 10.5 Hz, and FI-EVBN
departs from FV by up to 24 % near 19 Hz: its own enclosed-soil frequency, 22.4 Hz on this lesson's
site (an eigenvalue analysis in examples/README.md), falls to about 15 Hz. With
$V_s = 600\,\text{m/s}$ (`L,k,2.0,19.0,1200,600,0.04,0.04`) FI-FSIN follows FV within 1.3 %.

Run FFV as well: the commands of the second box of this step, on the model the step prepared
(about 35 s).

### Check yourself
A reactor building on a soft-soil site is embedded 20 m. The project proposes FI-EVBN, validated on
the best-estimate soil case. What do you ask for?

Answer: the validation on the softest (lower-bound) soil case, because the enclosed-soil frequencies
fall with $V_s$ and grow lower with depth; transfer functions over the full range at the common
nodes and at the equipment locations, on the full model (not a quarter model); evidence that the
frequency set resolves every peak; and, for deep embedment, a comparison with FFV, which the manual
recommends there.
