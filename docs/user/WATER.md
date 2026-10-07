# Water modelling in SASSI-EDU (FILLPOOL, REFINEMODEL, LISTPOOLINTER)

Pools and tanks such as spent-fuel pools, cooling-water basins and elevated tanks add **hydrodynamic
mass** to a structure: when the walls accelerate, part of the water is pushed along with them and
pushes back on the walls with a dynamic pressure. ACS SASSI models this water with ordinary SOLID
elements that have the properties of water. Three commands (manual section 9.16) support this:

| Command | What it does |
|---|---|
| `FILLPOOL,<Stiff>,<Sensitivity>,<EmptyLevels>,<ShellArea>,<offset>,<stiff2>` | Fills a pool sub-model with water SOLIDs, interface springs and, optionally, interface-area shells |
| `REFINEMODEL` | Splits every quadrilateral shell into 4 and every hexahedral solid into 8 |
| `LISTPOOLINTER,<Pool>` | Lists the pool interface nodes that the original model has with the same number and position |

SASSI-EDU adds two commands that are not in the manual:

| Command | What it does |
|---|---|
| `MERGEPOOL,<Pool>` | Imports the water, the springs and the area shells of a filled pool into the active (original) model, keeping the node numbers |
| `POOLDATA,...` | Storage record. WRITE writes it so that LISTPOOLINTER and MERGEPOOL still work after WRITE → INP. You do not type it |

Where the code lives:

* algorithms: `sassi/prep/water_lib.py`;
* commands: `sassi/prep/commands/water.py`;
* verification: `sassi/verify/problems/vp_water.py` (VP-54W, VP-W1);
* tests: `tests/unit/test_water*.py`, `tests/verification/test_vp_water.py`.

The specification is in requirements §3.4.N, spec 11 §1 and decision D-WAT-01. D-WAT-01 is amended
here: see §5.

> **In one sentence:** the water SOLIDs reproduce the *impulsive* hydrodynamic mass of Westergaard and
> Housner. They do **not** reproduce sloshing (the *convective* mass), they need incompatible-mode
> SOLIDs (`MOPT,0`, which FILLPOOL sets for you), and their response below about 0.5 Hz has no meaning.

---

## 1. Quick start

### 1.1 A pool inside a building (the manual's workflow)

```
* model 0 is the building; group 7 holds the pool walls and floor (SHELLs)
ACTM,0
* the pool's elements, copied to sub-model 5 (same node, group and element numbers)
CUTADD,1,7,RANGE,1,240
CUT2SUB,1,5
ACTM,5
* optional: a finer wall mesh gives a finer water mesh
REFINEMODEL
* offset 12000 >= the building's last node number
FILLPOOL,1E8,0,1,-1,12000,0
ACTM,0
* every interface node must be found (same number and position)
LISTPOOLINTER,5
* water, springs and area shells copied into the building
MERGEPOOL,5
* as usual for shell models
FIXROT
```

(Comments in a `.pre` file are whole lines that start with `*`; there are no inline comments.)

* `FILLPOOL` fills the pool. Here the first Z-level below the wall tops stays empty (`EmptyLevels` = 1,
  freeboard). The interface springs have a stiffness of 10⁸ along the wall normals (model units: kip/ft
  here, §4.4) and zero along the walls. No interface-area shells are made (`ShellArea` = −1). The water
  nodes are numbered from 12001.
* `LISTPOOLINTER` checks that the wall nodes the springs are attached to still exist in the building,
  with the same number and position.
* `MERGEPOOL` appends the FILLPOOL groups after the building's last group. It renumbers their materials
  and spring properties after the building's, and keeps all node numbers.

### 1.2 A tank on its own (as in VP-W1)

```
* gravity first: it sets the water density (kip, ft, s)
HOUSE,32.2,0,0,2,0,0,0,0,0
* ... N / E / THICK commands of the tank walls and floor (SHELL), in ft ...
* fill up to 2 levels below the wall tops; interface springs 1E8 kip/ft
FILLPOOL,1E8,0,2
FIXROT
* interaction nodes: the tank floor only (see 4.8)
INT,1,63,1,1
* ... SITE / POINT / ANALYS / FREQ / AOPT, CHECK, AFWRITE, RUNSITE ... RUNANALYS
```

VP-W1 builds the same kind of model in SI units (`HOUSE,9.81`, springs 10⁹ kN/m): see `tank_lines()` in
`sassi/verify/problems/vp_water.py`.

---

## 2. The commands

### 2.1 FILLPOOL

`FILLPOOL,<Stiff>,<Sensitivity>,<EmptyLevels>,<ShellArea>,<offset>,<stiff2>`. The defaults are 10⁶, 0,
0, −1, −1 and 0.

The active model must be a **pool sub-model**: it holds only the walls and floor of one pool, made of
SHELL/TSHELL or SOLID elements. Make it with CUT2SUB (or CPMODEL + GDEL).

| Argument | Meaning |
|---|---|
| `Stiff` > 0 | Stiffness (force/length) of each interface spring along the wall normal |
| `Sensitivity` ≥ 0 | Allowed variation of Z among the nodes of one Z-level (as EXCAV `delta`). Also the coincidence distance between a water node and a wall node |
| `EmptyLevels` ≥ 0 | Number of Z-levels, counted down from the highest, left without water (freeboard) |
| `ShellArea` | `1` creates the interface-area shells. Any other value creates none |
| `offset` | The water nodes are numbered from `offset + 1`. A value ≤ 0 uses the largest node number of the pool model. A positive value below it is an error. To import the water into the original model, use a value ≥ the original model's last node |
| `stiff2` ≥ 0 | Interface stiffness *along* the wall. Use 0 for water (frictionless). Use > 0 for a pool material that carries shear |

**How the pool is filled.** FILLPOOL uses the EXCAV algorithm:

1. The elevations of all wall and floor nodes are grouped into **Z-levels**. A new level starts when the
   gap to the previous one exceeds `Sensitivity`. Each level's elevation is the mean of its nodes.
2. The **floor template** is the plan mesh of the water. It is the lowest band of horizontal SHELLs and
   of free (unshared), upward-facing SOLID faces. For a SOLID pool, that is the top of the floor slab
   inside the walls.
3. Each template cell is extruded between consecutive levels, from the floor up to the level
   `EmptyLevels` below the top. A quadrilateral cell gives an 8-node SOLID; a triangular cell gives a
   prism with repeated nodes. The water nodes are new nodes, numbered bottom-up, level by level. They
   are separate from the wall nodes even where the two coincide.

**What is added.** The new groups go after the model's last group:

| Group | Contents |
|---|---|
| `FILLPOOL water` | SOLID, ETYPE 1 (always structure, even below grade), the new water material |
| `FILLPOOL water-wall interface springs` | One zero-length SPRING per water node of the wetted boundary that coincides with a wall/floor node: node I = wall, node J = water |
| `FILLPOOL interface areas` (ShellArea = 1) | Dummy SHELLs on the wall-side nodes of every wetted water face |

The water nodes at the floor level and on the side boundary are attached. The free surface is not.

**Spring constants.** SASSI springs (SC) are uncoupled constants in global axes.

* At each attached node the springs act along every wetted-face normal. A floor node gets Z, a node on
  a wall parallel to the YZ plane gets X, a vertical corner gets X and Y, a wall foot gets X and Z.
* Along a normal the constant is `Stiff`; in the other directions it is `stiff2`. The rotational
  constants are 0 and the damping is 0.
* Normals of neighbouring boundary edges that differ by less than 30° are averaged (a polygonal
  "circular" wall counts as smooth); larger angles are corners.
* For walls that are **not** parallel to a global plane, see §4.7.

**Water material.** A new type-3 M entry (Vp, Vs) holds the D-WAT-01 water in the units implied by the
HOUSE gravity: kip-ft when g > 20 (g = 32.2), kN-m otherwise.

* K = 45,950 ksf (2.2 GPa; 2.2 × 10⁶ kPa in kN-m);
* unit weight 0.0624 kcf (62.4 pcf; 9.81 kN/m³ in kN-m);
* G = 10⁻⁸ K, so Vs = 10⁻⁴ Vp;
* damping 0.5 %;
* Vp ≈ 4,870 ft/s (1,483 m/s).

FILLPOOL prints the material number. For other units, redefine it with
`M,<n>,<Vp>,<Vs>,<weight>,<Dp>,<Ds>,3`.

**Other changes FILLPOOL makes:**

* **Rotations.** SPRING elements carry six DOFs per node, and nothing restrains the rotations they
  add at the water nodes, so FILLPOOL fixes them with D. It does the same for the interface nodes of
  SOLID walls. At flat SHELL-wall interface nodes it fixes the drilling rotation, as FIXROT does: once
  the springs exist, FIXROT no longer treats those nodes as shell-only and would skip them. Run FIXROT
  for the rest of the model as usual.
* **Incompatible modes.** MOPT `<incomp>` is set to 0 (with a warning; see §4.3).
* **Pool data.** The `POOLDATA` record is stored (§2.4).

**What FILLPOOL prints:**

* the Z-levels and the water depth;
* the water material;
* the volume, mass and wetted area;
* the spring properties;
* the "frequency of the water mass on the interface springs" (§4.4);
* an estimate of the **sloshing** that is *not* modelled (§4.1);
* warnings for any water boundary nodes that have no wall node (non-conforming wall mesh), more than
  one wall node (unwelded junction) or a wall inside the water, and for a roof that touches the water
  surface (the surface stays free: use `EmptyLevels` ≥ 1).

### 2.2 REFINEMODEL

`REFINEMODEL` takes no arguments.

* Every quadrilateral SHELL, TSHELL and PLANE element (4 distinct nodes) is split into 4.
* Every hexahedral SOLID (8 distinct nodes) is split into 8.
* Triangles, prisms and pyramids (repeated nodes), BEAMS, SPRING and GENERAL elements are not changed.

The splits use the edge midpoints, plus the face centres and the body centre that a split needs. All
are shared between neighbours, so a conforming mesh stays conforming. Each centre is the average of
its corners, which is the bilinear or trilinear centre. Every child is therefore an exact piece of its
parent: the area and volume are preserved, and so is the node order (orientation, positive Jacobian).

Children inherit the material, property, ETYPE, EINT, THICK and the beam releases.

New nodes:

* are numbered after the last node, in order of creation;
* get an INT code, or a D fixity, only when **every** parent corner node has it (an edge's two ends, a
  face's four corners, a hexahedron's eight);
* that become interaction nodes may lie between soil-layer interfaces, and REFINEMODEL then warns
  (EDU-01: refine the TOPL layers too).

The elements of a refined group are renumbered 1…n, with the children of each parent consecutive. EOUT
requests and the session cuts are updated to list the children. Nodal loads and masses are **not**
redistributed (a warning). Where an unrefined triangle or prism borders a refined element there are
hanging nodes (a warning).

Refine a pool's walls **before** FILLPOOL; the interface springs are not regenerated.

### 2.3 LISTPOOLINTER

`LISTPOOLINTER,<Pool>`. Run it with the **original** model active (ACTM) and the filled pool sub-model
`<Pool>` in memory. The pool's interface nodes are the wall nodes of its FILLPOOL springs.

* The command lists those that the active model has with the **same number and the same position**
  (within the geometric tolerance).
* Interface nodes that are missing from the active model, or that lie somewhere else in it, are
  reported as warnings. These are the nodes to fix before importing the water.
* It works only for a pool filled with FILLPOOL.

### 2.4 MERGEPOOL and POOLDATA (SASSI-EDU)

The manual says to "import the water and spring group back into the original model". The general MERGE
command cannot do that, because it always renumbers the second model's nodes (D-MDL-12): the water
would be joined to copies of the walls.

`MERGEPOOL,<Pool>` imports the water without renumbering nodes. It first checks:

* the LISTPOOLINTER condition (every interface node found, same number and position);
* that no water node number is already in use. Fill with `offset` ≥ the original model's last node,
  or check that the numbers are free.

It then:

* appends the FILLPOOL groups with their element numbers;
* gives the water and shell materials and the spring properties new numbers after the active model's;
* copies the water nodes with their fixed rotations;
* fixes the rotations of interface wall nodes that nothing restrains in the active model (the same
  rules as FILLPOOL);
* sets MOPT `<incomp>` = 0.

A second MERGEPOOL of the same pool is refused.

`POOLDATA` is the record FILLPOOL stores. It holds:

* the water, spring and shell group numbers;
* the FILLPOOL arguments;
* the water material;
* the floor and surface elevations.

WRITE writes it after the groups (with a comment line), so INP restores it. `POOLDATA` with no
arguments deletes it.

---

## 3. Why it works: water as a "soft solid"

### 3.1 From a solid to an inviscid fluid

An isotropic elastic solid with bulk modulus K and shear modulus G obeys

&nbsp;&nbsp;&nbsp;&nbsp;ρ ü = (K + 4G/3) ∇(∇·u) − G ∇×(∇×u)

If G → 0, the shear term disappears and the only stress is the pressure p = −K ∇·u:

&nbsp;&nbsp;&nbsp;&nbsp;ρ ü = −∇p,  so  p̈ = c² ∇²p with c = √(K/ρ) ≈ 4,870 ft/s.

This is the **linear acoustic equation** of an inviscid compressible fluid, written in displacements. A
SOLID with the bulk modulus of water and almost no shear modulus therefore behaves like water, as long
as nothing excites its (spurious) shear modes; see §4.2.

### 3.2 The incompressible limit: Westergaard and Housner

Under seismic excitation the walls move at frequencies far below the compression frequencies of the
water. The first is f_a = c/(4H), about 76 Hz for H = 16 ft of water. In that range the water is
practically incompressible, ∇²p = 0, with these boundary conditions:

* on a wall or the floor, the water follows the wall's normal acceleration: ∂p/∂n = −ρ aₙ (the
  interface springs impose this);
* the water is free to slide along the walls (no shear: `stiff2` = 0);
* at the free surface p = 0. Without gravity, the surface can heave freely.

This is the potential-flow problem of Westergaard (1933) and Housner (1963). The water pushed by the
walls moves with them and adds the **impulsive mass** m_i. The rest of the water lags behind, which
the free surface allows.

For a rigid rectangular tank (length 2L in the direction of the motion, water depth H), the series
solution is

&nbsp;&nbsp;&nbsp;&nbsp;m_i/m = 2/(L H²) Σₙ tanh(λₙ L)/λₙ³,  λₙ = (2n+1)π/(2H)

(`water_lib.impulsive_ratio_exact`). Its properties:

* it is exactly ½ for L = H;
* it tends to 1 for narrow tanks;
* it tends to 0.5428 H/L for long ones (Westergaard's dam);
* it satisfies r(L/H) + r(H/L) = 1.

Two classical closed forms approximate it:

* **Housner** (1963): tanh(√3 L/H)/(√3 L/H). It is 6–11 % conservative.
* **Westergaard** (1933): (7/12) H/L, for long reservoirs only.

The resultant of the wall pressures acts at 0.40 H above the floor (exact series; Housner: 3/8 H,
excluding the base pressure).

### 3.3 What VP-W1 shows

VP-W1 is a rigid tank 10 m × 2 m filled to 5 m (L/H = 1), on a practically rigid site. It runs the full
chain SITE → POINT → HOUSE → ANALYS from commands. The table gives the base shear, which is the sum of
the X forces of the interface springs, divided by the water mass times the input acceleration:

| Frequency | F/(m a), FE (0.5 m mesh) | exact m_i/m | Housner | Westergaard |
|---|---|---|---|---|
| 2 Hz | 0.5073 | 0.5000 | 0.5423 | 0.5833 |
| 4 Hz | 0.5077 | 0.5000 | 0.5423 | 0.5833 |
| 8 Hz | 0.5093 | 0.5000 | 0.5423 | 0.5833 |

A second tank with SOLID walls 0.5 m thick on a 0.5 m SOLID slab, with the interface-area shells, gives
the same values to 10⁻⁵ (0.50729, 0.50771, 0.50929). It has the same water mesh, and this checks the SOLID-wall
path of FILLPOOL: the floor template from the slab top and the rotation fixes.

The FE value is 1.5–1.9 % above the exact potential flow. That is the mesh discretisation error (§4.5)
plus a small compressibility effect that grows with frequency (§4.6). Housner's design formula is
6.4 % higher. The resultant of the wall pressures acts at 0.4096 H (exact 0.4047 H).

Exact properties are checked to round-off:

* the water mass in the HOUSE mass matrix is ρV;
* the spring forces balance the water's inertia;
* the tank follows the control motion.

---

## 4. Limits and modelling advice

### 4.1 Sloshing (the convective mass) is not modelled

Surface waves need the restoring force ρg at the free surface, which an elastic solid does not have.
The **convective** part of the water and its sloshing frequencies (a fraction of a hertz for pools) are
therefore missing.

For seismic analysis of the structure this is usually acceptable. The convective mass responds at
very low frequency, where spectral accelerations are small, and it is only loosely coupled to the
structure. It matters for freeboard (wave height) and for the pool walls at low frequencies.

FILLPOOL prints an estimate from the plan extents, using Housner's rectangular-tank model
(`water_lib.housner_convective`, with k = √(5/2)):

&nbsp;&nbsp;&nbsp;&nbsp;m_c/m = (k/3)(L/H) tanh(k H/L),  ω_c² = (k g/L) tanh(k H/L)

The exact first sloshing frequency is ω² = (πg/2L) tanh(πH/2L) (`water_lib.sloshing_frequency_exact`).

If the convective response matters, add it yourself as Housner's spring-mass: a mass m_c attached to
the walls at the height h_c by a spring of stiffness m_c ω_c². Reduce the water accordingly, or accept
a small double counting.

### 4.2 Spurious low-frequency modes

G is small but not zero, so the water has many "shear" modes at very low frequency. With G = 10⁻⁸ K
(Vs ≈ 0.5 ft/s) they lie below about 0.05 Hz for pools more than a few feet deep. Far above them the
water behaves as an inviscid fluid (§5).

* Do not interpret the water response below about 0.5 Hz.
* Do not use f = 0: the static problem is nearly singular.

### 4.3 Volumetric locking: the water needs incompatible modes (MOPT,0)

A trilinear SOLID integrated with 2×2×2 Gauss points enforces the volume constraint of a nearly
incompressible material at too many points. It *locks*: it can hardly deform at constant volume.
Locked water moves with the walls or resonates at wrong frequencies (§5, last row).

The Wilson–Taylor incompatible modes (D-ELM-02) remove the locking on regular meshes. FILLPOOL and
MERGEPOOL therefore set MOPT `<incomp>` = 0. This also switches the incompatible modes on for every
other structural SOLID of the model; they are then more accurate in bending.

Hexahedral water elements behave best. Triangles in the floor template give prisms, which are less
accurate.

### 4.4 Interface spring stiffness

The springs impose the normal coupling. They must be stiff compared with the water itself:

* FILLPOOL prints, per direction, the frequency at which the whole water mass would vibrate on the
  springs alone. Keep it far above the frequencies of interest.
* The default `Stiff` = 10⁶ is in model units. In kip-ft it is about 10 times the water's own bulk
  stiffness K × (element size) for 2 ft elements (K h = 45,950 ksf × 2 ft ≈ 9 × 10⁴ kip/ft), which is
  stiff enough (table below); in kN-m it is only comparable to K h ≈ 10⁶ kN/m (0.5 m elements).
* The examples of §1 use 10⁸ kip/ft (about 1,000 K h for 2 ft elements); VP-W1 uses 10⁹ kN/m.

At low frequency the result is not very sensitive. In the 2D slice of §5 with a 2 ft mesh
(L = H = 16 ft) and an interface spring at every wetted node, F/(m a) changes as follows:

| Stiff / (K h) | 1 Hz | 2 Hz | 5 Hz | 10 Hz |
|---|---|---|---|---|
| 0.1 | 0.5106 | 0.5114 | 0.5164 | 0.5352 |
| 1 | 0.5104 | 0.5107 | 0.5118 | 0.5157 |
| 10 | 0.5104 | 0.5106 | 0.5113 | 0.5138 |
| ≥ 100 | 0.5104 | 0.5106 | 0.5113 | 0.5136 |

### 4.5 Mesh

The FE solution converges to the exact potential flow from above, at about O(h²). At 2 Hz, with
L = H = 16 ft (the 2D slice of §5, square elements):

| element size | F/(m a) | error |
|---|---|---|
| 4 ft (4 over the depth) | 0.5333 | +6.7 % |
| 2 ft (8 over the depth) | 0.5106 | +2.1 % |
| 1 ft (16 over the depth) | 0.5033 | +0.65 % |

Use about 8–10 water elements over the depth for 2 %. Refine the pool walls with REFINEMODEL before
FILLPOOL if needed.

### 4.6 Compressibility

The incompressible references hold only well below the first compression frequency c/(4H) of the
water column. As that frequency is approached, the effective mass grows (dynamic amplification):

* for H = 16 ft (c/(4H) = 76 Hz), F/(m a) is 0.6 % higher at 10 Hz than at 2 Hz, 2.4 % higher at
  20 Hz and 5.8 % higher at 30 Hz (the §5 slice, G = 10⁻⁸ K row);
* deep pools have lower compression frequencies (40 ft of water: 30 Hz).

### 4.7 Oblique walls (D-WAT-01 / OQ-2)

For a wall that is not parallel to a global plane, the exact interface is the coupled 3×3 matrix
`Stiff P_N + stiff2 (I − P_N)`, where P_N projects on the wall normals.

D-WAT-01 proposed zero-length 2-node GENERAL elements for these walls. CHECK rejects them with Error 8
("Has 0 Length", D-CHK-07), so they cannot reach AFWRITE. FILLPOOL therefore uses the **diagonal** of
that matrix in a SPRING (OQ-2 (b)) and warns.

The diagonal keeps the normal restraint but adds a tangential one (for a 45° wall, `Stiff`/2 in both
directions). The water partly **sticks** to oblique walls. Because G ≈ 0, this affects only the first
row of water elements, and the hydrodynamic mass is slightly overestimated.

Corner nodes, where the normals span the whole horizontal plane, and floor nodes stay exact.

### 4.8 Interaction nodes

`INTGEN,4` (surface foundation) takes every SOLID/SHELL node at the ground elevation. That includes
the bottom water nodes of a tank standing at grade. Define the interaction nodes on the tank floor with
INT (as in VP-W1), or reset the INT flag of the water nodes.

### 4.9 Pressures, hydrostatics and units

* The analysis is dynamic. The hydrostatic pressure is not part of it.
* The hydrodynamic pressure in a water SOLID is p = −(SXX + SYY + SZZ)/3 (STRESS). With the
  interface-area shells, the spring forces divided by the shell areas give the wall pressures (for
  example, for an ANSYS model).
* The water material is defined for kip-ft or kN-m. With other units, redefine it with `M`. FILLPOOL
  warns when the gravity is neither ~32.2 nor ~9.81.

---

## 5. Choice of the shear modulus (D-WAT-01 amended)

D-WAT-01 proposed ν = 0.49, which is G = 0.0201 K. The table shows why FILLPOOL uses G = 10⁻⁸ K. It
gives F/(m a) for the 2D slice of tests/unit/test_water_physics.py, here in kip-ft with L = H = 16 ft
(16 × 16 elements, exact value 0.500):

| water model | 0.2 Hz | 0.5 Hz | 1 Hz | 2 Hz | 5 Hz | 10 Hz | 20 Hz | 30 Hz |
|---|---|---|---|---|---|---|---|---|
| ν = 0.49 (G = 0.0201 K) | 1.000 | 1.000 | 1.002 | 1.007 | 1.045 | 1.272 | 0.417 | 0.046 |
| G = 10⁻⁵ K | 1.190 | 0.704 | 1.139 | 0.498 | 0.504 | 0.509 | 0.520 | 0.537 |
| G = 10⁻⁶ K | 0.686 | 0.475 | 0.490 | 0.503 | 0.508 | 0.510 | 0.520 | 0.537 |
| G = 10⁻⁷ K | 0.498 | 0.503 | 0.506 | 0.507 | 0.508 | 0.511 | 0.520 | 0.537 |
| **G = 10⁻⁸ K (FILLPOOL)** | **0.503** | **0.507** | **0.507** | **0.508** | **0.508** | **0.511** | **0.520** | **0.537** |
| G = 10⁻⁸ K, MOPT,1 (locked) | 1.002 | 1.011 | 1.048 | 1.286 | 0.501 | 0.457 | 0.527 | 0.542 |

What the rows show:

* With ν = 0.49 the "water" is an elastic solid. Its shear modes lie in the seismic band (the first
  resonance of the slice is near 15 Hz for 16 ft of water), and below them it moves rigidly with the
  walls: the whole mass instead of half of it.
* The smaller G is, the lower the spurious modes, and the wider the band in which the water behaves
  as an inviscid fluid.
* G = 10⁻⁸ K gives that behaviour from 0.2 Hz up.

A very small G is harmless numerically. At every non-zero frequency the inertia term controls the
near-zero-stiffness shear deformations, and K/G = 10⁸ is far from the double-precision limit.

---

## 6. Verification

| ID | What is verified | Result |
|---|---|---|
| VP-54W | REFINEMODEL: 1 quad → 4 (5 new nodes), 2×2 → 16 quads / 25 nodes, hexahedron → 8 (19 new nodes), area and volume preserved (1e-12), node order kept, triangles and prisms untouched, hanging nodes reported. FILLPOOL: volume a·b·(H − m·h) (1e-12), one spring per coincident pair, `Stiff` along normals and `stiff2` along walls at walls, floor, corners, wall foot and free surface, offset error with the model unchanged, water K = 2.2 GPa and G/K = 10⁻⁸, mass ρV, area shells = wetted area, MOPT,0, WRITE → INP round trip; SOLID pool. CUT2SUB → FILLPOOL → LISTPOOLINTER (all found) → MERGEPOOL (node numbers kept, second import refused) | 54 checks pass |
| VP-W1 | Rigid tank (SHELL walls; and SOLID walls on a slab with area shells), full chain: water mass (error 3 × 10⁻¹⁶), equilibrium (≤ 10⁻⁸), rigid tank (max \|H − 1\| = 1.0 × 10⁻⁴, criterion 10⁻³), F/(m a) = exact m_i/m within 5 % (observed +1.5 % … +1.9 %); Housner, Westergaard and h_i reported | passes |

The unit tests (`tests/unit/test_water.py`, `test_water_refine.py`, `test_water_physics.py`) cover:

* every argument and error of the commands;
* oblique and triangular pools, jittered levels with `Sensitivity`, British units;
* the rotation fixes, POOLDATA, LISTPOOLINTER mismatches and MERGEPOOL;
* the analytical references;
* the physics of §4–§5.

---

## 7. Decisions and differences from the manual

* **D-WAT-01 (amended):** G = 10⁻⁸ K instead of ν = 0.49 (§5). The other values are unchanged: K, unit
  weight, 0.5 % damping, spring damping 0, Z springs at floor nodes, dummy area shells.
  * The area shells use E = 10⁻⁶ K_water, ν = 0.3, zero weight and zero damping, and a thickness of
    1/100 of the mean wetted-face edge.
  * Delete the area-shell group (GDEL) if you do not need it.
* **Oblique walls:** a diagonal SPRING instead of a GENERAL element (§4.7; CHECK Error 8).
* **Rotations and MOPT:** FILLPOOL and MERGEPOOL fix the rotations the springs add and set MOPT,0. The
  manual is silent on both; without them the system is singular or the water locks.
* **offset:** applies to node numbers. The elements of the new groups are numbered from 1, as in every
  SASSI group.
* **MERGEPOOL, POOLDATA:** SASSI-EDU additions (§2.4).
* **REFINEMODEL (OQ-6):** also refines PLANE quadrilaterals. It inherits INT codes and fixities by the
  all-corners rule, renumbers refined groups and remaps EOUT and cuts.
* **LISTPOOLINTER:** lists *nodes*, as the manual's detailed description says (its first sentence says
  "elements"), and also reports the mismatches.

---

## 8. References

* ACS SASSI Version 3 User Manual, section 9.16 (Water Modeling Commands).
* Westergaard, H. M. (1933). Water pressures on dams during earthquakes. *Transactions of the ASCE*,
  98, 418–433.
* Housner, G. W. (1963). The dynamic behavior of water tanks. *Bulletin of the Seismological Society of
  America*, 53(2), 381–387. Also: U.S. AEC, *Nuclear Reactors and Earthquakes*, TID-7024 (1963),
  chapter 6 and appendix F.
* Wilson, E. L., Taylor, R. L., Doherty, W. P. and Ghaboussi, J. (1973). Incompatible displacement
  models. In *Numerical and Computer Methods in Structural Mechanics*, Academic Press, 43–57.
* Taylor, R. L., Beresford, P. J. and Wilson, E. L. (1976). A non-conforming element for stress
  analysis. *International Journal for Numerical Methods in Engineering*, 10, 1211–1219.
