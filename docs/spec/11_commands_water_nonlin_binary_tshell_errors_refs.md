# Spec 11: Water Modeling, Nonlinear & Panel, Binary Database and Thick Shell Commands; CHECK Errors & Warnings; References

**Source:** ACS SASSI Version 3 User Manual, text lines 13261–15558 (`reference/acs-sassi.txt`).
This covers printed pages 312–357 (PDF pages 314–359). For these pages, PDF page = printed page + 2.

**Sections covered:**
* 9.16 Water Modeling Commands
* 9.17 Nonlinear & Panel Commands
* 9.18 Binary Commands
* 9.19 Thick Shell Commands
* Chapter 10 Check Errors and Warnings
* Chapter 11 References

The equation and the dialog-like listings were checked against the PDF page images (pp. 314–336 and 337–353). The
NLSLAYER equation is the only equation in Ch. 9.16–9.19, and it was transcribed from the image. The SHEAR
equations, which BBCGEN depends on, were transcribed from PDF pp. 291–292 for cross-reference.

**Tags:**
* **[CORE]** needed for computational correctness
* **[IO]** file or input format
* **[UI]** user interface, plotting or convenience (model-generation utilities are included here)
* **[ADV]** advanced options: Option NON (nonlinear), incoherency, Options A/AA/PRO

**Related specs:**
* `05d_options_stress_nonlinear_settings.md`: NONLINEAR tab, Option NON algorithm, `.eql`, hysteretic models,
  BBC conventions, RELDISP/STRESS binary options.
* `05a_options_equake_soil_site_point.md`: SOIL tab including the Nonlinear Soil group, SITE and POINT checks.
* `05b_options_house_force.md`: HOUSE incoherency and multiple-excitation checks.
* `05c_options_analys_motion.md`: ANALYS and MOTION checks.

This spec uses the same codes and terms as those files.

---

## 0. Conventions common to all commands in this range

### 0.1 Command syntax [IO]

* Commands are comma-separated: `NAME,arg1,arg2,…`.
* `<arg>` is required and `[arg]` is optional.
* A blank argument (for example `BINOUT,,1`) means "use the default". For `BINOUT` it means "leave the current
  value unchanged".
* Commands work on models held in **ACS SASSI UI memory**:
  * several models can be loaded at once, each identified by a **model number**;
  * the **active model** is selected with `ACTM,<model>`.
  * Commands that take a model number (`LISTPOOLINTER`, `MERGEPANEL`, `SOILREDEF`) need that model to be loaded.
* Names are case-insensitive, as in the rest of the UI. The manual does not state this, but the examples assume it.

### 0.2 Range/stride deletion and listing convention [UI]

`DELBM`, `DELNLS`, `DELSPR`, `PDEL`, `PLIST` and `DELBBC` all take `start,[end],[stride]`:

* Items handled: `start, start+stride, start+2·stride, …` while the item is ≤ `end`.
* `end` defaults to `start`. For `PLIST`, `end` defaults to the largest panel number in memory and `start`
  defaults to 1.
* `stride` defaults to 1. If `stride` is **blank or < 1**, the command uses 1 and prints a **warning** that the
  default stride was used. This is stated for DELBM, DELNLS/DELSPR and PDEL.
* Numbers in the range that do not exist are skipped. The manual does not say this; it is the implementer's
  choice (see OQ-11).

### 0.3 Status flags printed in the manual (must be reproduced in the UI help/status) [UI]

| Command | Manual status |
|---|---|
| B | "Not usable in this version" |
| BEAMPILE | "Not applicable in this version" |
| DCOUPLEBEAM | "Not usable in this version" (body text: "Not applicable") |
| SOILREDEF | "NOT VALIDATED in this version" |
| THSHLSMH | "NON VALIDATED in this version" |
| BINOUT `reldisp=1` (TFD binary) | "Not used in this version" |
| NLSOIL `BedInt=1` (viscoelastic) | "Currently Disabled" |

Recommendation for the re-implementation:
* Commands marked not usable or not applicable should parse their arguments and store them, then print a "not
  available in this version" message.
* Commands marked not validated should run but print a "not validated" warning.

---

## 1. Water Modeling Commands (Manual 9.16) [UI]

| Command | Action | Section |
|---|---|---|
| FILLPOOL | Fill a shell or solid pool with water solids | 9.16.1 |
| REFINEMODEL | Subdivide all shell and solid elements in a model | 9.16.2 |
| LISTPOOLINTER | List the interface nodes of a pool in the original model | 9.16.3 |

### 1.1 FILLPOOL

**Syntax:** `FILLPOOL,<Stiff>,<Sensitivity>,<EmptyLevels>,<ShellArea>,<offset>,<stiff2>`

**Purpose.** FILLPOOL works on the active model, which must be a pool sub-model. It adds three things:
1. a **group of SOLID "water" elements** that fill the pool volume;
2. a **spring interface** connecting the water to the pool walls and floor;
3. optionally, a **SHELL group** that approximates the interface areas. These shells are useful for later
   processing of the pool in ANSYS, for example to convert spring forces to pressures.

**Fill algorithm.** The fill uses the same algorithm as `EXCAV`. The pool floor gives the plan template, and
the **Z-levels of the wall nodes** give the layers of solids.

**Preconditions [UI]:**
* Create a **sub-model containing only the walls and floor of a single pool**, then activate it.
* Walls and floor must be made of **SHELL elements or SOLID elements only**.

**Arguments:**

| Arg | Type | Default | Units | Valid | Meaning | Tag |
|---|---|---|---|---|---|---|
| `Stiff` | real | **1.0e6** | force/length (model units) | > 0 | Stiffness of the water–wall interface springs **along the wall normal** | [CORE] |
| `Sensitivity` | real | 0 | length | ≥ 0 | Allowed variation in Z for nodes on the **same Z-level**. Same role as EXCAV `delta`. Use > 0 when the floor or level Z values vary slightly. | [UI] |
| `EmptyLevels` | int | 0 | – | 0 ≤ n < number of levels | Number of Z-levels, counted down from the **highest** level, that are **not** filled with water. Models freeboard, i.e. a water surface below the top of the wall. | [UI] |
| `ShellArea` | int | −1 | – | any | `1` creates the interface-area shell group. Any other value (≠ 1) does not. | [UI] |
| `offset` | int | −1 | – | ≤ 0, or ≥ the max node number of the pool walls | Number from which the new pool nodes and elements are numbered. **≤ 0** uses the maximum node number of the pool-wall model. A **positive value smaller than that maximum is an ERROR**. To import the water and springs back into the original model, use a value ≥ the original model's last node number. | [IO] |
| `stiff2` | real | 0 | force/length | ≥ 0 | Interface stiffness **perpendicular to the normal** (tangential). Used when the pool holds a non-fluid material that can carry shear. 0 means a frictionless fluid. | [CORE] |

**Algorithm.** The manual describes only the outcome. The steps below are the proposed implementation.

1. **Validate.** The active model has only SHELL and/or SOLID elements. Determine `Nmax_pool` = max node number.
   Resolve `offset`:
   * if `offset ≤ 0`: `offset = Nmax_pool`;
   * else if `offset < Nmax_pool`: stop with an error.
2. **Z-levels.** Collect the distinct Z coordinates of all wall nodes. Sort them ascending and merge values that
   lie within `Sensitivity` into one level. Use a representative value per cluster, for example the mean (OQ-1).
3. **Floor template (plan mesh).**
   * Shell floor: the horizontal floor shells, i.e. normal ∥ Z, at the lowest level.
   * Solid floor: the **free (unshared) upward-facing faces** at the floor-top level. These are the faces that
     will touch the water.
   * The template is the 2D quad or triangle mesh made of these faces.
4. **Fill levels.** Let the levels be L0 (floor or floor-top) < L1 < … < Ln (top of the wall). Fill each
   interval [Lk, Lk+1] for k = 0 … n − 1 − `EmptyLevels`. Each interval is extruded from the template:
   * quad → 8-node SOLID;
   * triangle → prism, i.e. a SOLID with repeated nodes.
5. **Water nodes.**
   * Create new nodes at every template (x, y) for every filled level, numbered from `offset + 1`.
   * Water nodes are **separate** from the wall nodes, even where they coincide in space.
6. **Interface springs.** Each water node that coincides with a wall or floor node (within tolerance) gets a
   2-node SPRING between the wall node and the water node:
   * normal stiffness = `Stiff`;
   * tangential stiffness = `stiff2`;
   * rotational stiffness = 0.
   * SASSI spring constants (`SC`) are in **global** xyz and uncoupled. So for walls parallel to the global
     planes, the normal is the global X, Y or Z axis. Oblique walls need an approximation (OQ-2).
   * Damping of the interface springs is not specified (OQ-3).
7. **Water material.** Not specified by the manual (OQ-4). See the proposal in OQ-4.
8. **Shell-area group** (if `ShellArea = 1`). Add a SHELL group on the wetted wall and floor faces, built from
   the wall-side nodes, so each wetted face has an area. Thickness, material and mass are not specified (OQ-5).
   Proposed: massless, negligible stiffness.
9. **Append.** Add the new groups to the end of the pool model's group list:
   1. water SOLIDs;
   2. interface SPRING(s);
   3. optional SHELL areas.

**Workflow [UI]:**
1. Make a sub-model of the pool and fill it with `FILLPOOL`.
2. `ACTM` back to the original model and check node correspondence with `LISTPOOLINTER,<pool>`.
3. Merge the water and spring groups back into the original model. The general `MERGE` command can do this with
   zero offsets.

**Warning to reproduce [UI]:** an error is reported for a positive `offset` below the pool-wall maximum node number.

### 1.2 REFINEMODEL

**Syntax:** `REFINEMODEL` (no arguments)

**Purpose.** REFINEMODEL subdivides every **quad SHELL** and every **non-prism SOLID** (hexahedron) in the active
model:
* each quad shell → **4 shells**, split at the **edge midpoints**;
* each hexahedral solid → **8 solids**, split at the **edge midpoints**.

**Triangular shells and prism solids are unchanged.**

**Algorithm (proposed) [UI]:**
1. Edge midpoints are shared between neighbouring elements. Key each edge by its sorted node pair so that
   adjacent elements reuse the same midpoint node (keeps the mesh conforming).
2. Splitting a quad into four quads also needs a **face-centre node**. Splitting a hexahedron into eight needs
   **6 face-centre nodes and 1 body-centre node**. The manual mentions only edge midpoints, but these extra nodes
   are geometrically required. Place them at the average of the corner nodes (bilinear/trilinear centre).
3. Quad split, corners n1…n4 with edge midpoints m12, m23, m34, m41 and centre c:
   `(n1,m12,c,m41) (m12,n2,m23,c) (c,m23,n3,m34) (m41,c,m34,n4)`. This keeps the node order and normal
   direction.
4. Hexahedron split: the standard 2×2×2 subdivision, keeping the node order.
5. Children inherit group, material, thickness, ETYPE and element options.
6. Renumber elements in each group contiguously, so there are no gaps (Error 6).
7. Append new nodes after the current max node number.

**Element classification:**
* A **triangle** is a 4-node shell with 3 distinct nodes.
* A **prism** is a SOLID with 6 distinct nodes.
* PLANE, BEAMS, SPRING and GENERAL elements are not mentioned and are left unchanged (OQ-6).

**Warnings to add (implementer):**
* Hanging nodes appear where a refined element borders an unrefined triangle or prism.
* New nodes on interaction surfaces may need interaction status. Advise re-running `INTGEN`/`INT`.
* Output requests, loads and masses are not redistributed.

### 1.3 LISTPOOLINTER

**Syntax:** `LISTPOOLINTER,<Pool>`

* `<Pool>`: model number of the pool sub-model.

**Behaviour [UI]:**
* The original model and the pool sub-model must both be loaded.
* The **original model must be active** (`ACTM`).
* The command compares the pool's **interface nodes** (the wall/floor nodes connected to water by FILLPOOL
  springs) with the original model's nodes. It lists the nodes that have **both the same node number and the
  same position**.
* It works **only for a pool filled with FILLPOOL**, which needs the stored interface metadata.

**Output:** a list of nodes printed to the command window. Proposal: also report interface nodes that do not
match, since those are the actual problem cases.

Note: the manual's first sentence says it lists "elements". The detailed description says it lists nodes.
Implement nodes.

---

## 2. Nonlinear & Panel Commands (Manual 9.17): Option NON model preparation [ADV]

### 2.0 Command index

| Command | Action | Section | Status |
|---|---|---|---|
| B | Nonlinear beam definition | 9.17.1 | not usable |
| BBC | Add a backbone curve (BBC) from a file | 9.17.2 | |
| BBCGEN | Generate BBC curves for defined panels | 9.17.3 | |
| BBCI | Set BBC information | 9.17.4 | |
| BBCP | Define a single point on a BBC | 9.17.5 | |
| BBCX | BBC definition, X values | 9.17.6 | |
| BBCY | BBC definition, Y values | 9.17.7 | |
| BEAMPILE | Separate beam piles and reconnect them with springs at the pile interface | 9.17.8 | not applicable |
| DCOUPLEBEAM | Add springs at beam intersections | 9.17.9 | not usable |
| DELBBC | Delete backbone curves | 9.17.10 | |
| DELBM | Delete nonlinear beam | 9.17.11 | |
| DELNLS | Delete nonlinear soil layers | 9.17.12 | text mis-printed (see §2.11) |
| DELSPR | Delete nonlinear spring | 9.17.13 | text mis-printed (see §2.11) |
| DGRDFLR | Modify Young's modulus of floor panels | 9.17.14 | |
| EDGE | Split a panel based on its edges | 9.17.15 | |
| EDGEMODEL | Apply EDGE to all wall groups in the model | 9.17.16 | |
| EQL | Nonlinear analysis options (`.eql` header) | 9.17.17 | |
| MERGEPANEL | Merge a panel model with the solids and beams of the original model | 9.17.18 | |
| NLSLAYER | Add a nonlinear soil layer definition | 9.17.19 | |
| NLSOIL | Set parameters for the nonlinear soil option | 9.17.20 | |
| NONLINBAT | Create a generic nonlinear batch run | 9.17.21 | |
| NONLINMOTDISP | Add panel corner nodes to the output request lists | 9.17.22 | |
| P | Add a panel to the active model | 9.17.23 | |
| PANELIZE | Separate shells in the model into panels | 9.17.24 | |
| PDEL | Delete panel(s) from the active model | 9.17.25 | |
| PLIST | List or check panel(s) in the active model | 9.17.26 | |
| PNLGEN | Create panel definitions for shell groups (called `PANELGEN` in Ch. 6) | 9.17.27 | |
| S | Nonlinear spring definition | 9.17.28 | |
| SOILREDEF | Redefine excavation-volume layers from a 2D soil model | 9.17.29 | not validated |
| SOLIDPILE | Separate solid piles and reconnect them with springs at the pile interface | 9.17.30 | |
| UNIPNL | Create a unique group for each element of a group | 9.17.31 | |
| WALLFLR | Separate shell walls and floors into separate groups | 9.17.32 | |

### 2.1 Nonlinear data model (UI memory, written to `.eql` by AFWRITE) [IO]

The `.eql` record layout is not documented (see 05d OQ-N1). The data model below is required.

```
BBC        { num:int>=1, type:int, npoints:int, yield:int (1-based index of yield point),
             X[npoints]:real, Y[npoints]:real }        # origin (0,0) is NOT stored; point 1 = cracking
Panel      { num:int, group:int, bbc:int, disp:int, force:int }
NLSpring   { num:int, group:int, elem:int, bbc:int, disp:int, force:int }
NLBeam     { num:int, group:int, spgroup:int, bbc:int, force:int, end1:int, end2:int }   # not usable
EQLHeader  { disp:real (EDF), NonLinOpts:int, dampCutoff:real(%), dampScale:real, ElasicD:0|1 }
NLSoilLayer{ num:int, curvefit:0|1, B:real, S:real, refStrain:real, Vis:real }        # SOIL module
NLSoilGlobal{ Opt:0|1, NSTimeSunInc:int, DispConv:real, ForceConv:real, EqualIt:int,
              BedInt:0|1, NLDampType:int, MMmult:real, SMmult:real }                   # SOIL module
```

**Code tables** (from 05d and the CHECK rules):

| Field | Codes | Supported in this version (CHECK) |
|---|---|---|
| BBC `type` | 1 CMS (Cheng–Mertz Shear), 2 CMB (Cheng–Mertz Bending, not included), 3 TAK (Takeda), 4 GMR (General Masing Rule, springs) | used for titles in the printout |
| Panel `disp` | 1 shear strain (average rotation of the two vertical edges); 2 bending rotation/curvature | **only 1** (Error 128) |
| Panel `force` | 1 CMS, 2 CMB, 3 TAK | **only 1** (Error 126) |
| Spring `disp` | 1 X, 2 Y, 3 Z translation | 1–3 |
| Spring `force` | 4 GMR | **only 4** (Error 127) |

### 2.2 Backbone curve definition commands [ADV][IO]

**BBC: load from a file**
* Syntax: `BBC,<num>,<type>,<points>,<yield>,<file>`
* The file is a list of **X Y pairs**: whitespace separated, one pair per line. X is the strain or displacement
  and Y is the force or moment.
* The command tries to read exactly `<points>` pairs. The user must know the number of points beforehand.
* Arguments: `num` curve number; `type` curve type (codes above); `points` number of points; `yield` index of
  the yield point; `file` file name.
* Implementation:
  * fewer pairs than `points` → error, curve not created;
  * extra lines → ignored with a warning.

**BBCI: set curve information**
* Syntax: `BBCI,<num>,<yield>,<type>`
* Sets the yield-point index and the type of curve `num`. Note the argument order differs from BBC.

**BBCP: set a single point**
* Syntax: `BBCP,<num>,<point>,<X>,<Y>`
* Sets point number `point` (1-based order along the curve) to (X, Y).
* The dialog grid cannot add or remove points (05d), so require 1 ≤ point ≤ npoints. Proposal: allow
  `point = npoints + 1` to append (OQ-7).

**BBCX / BBCY: set the X or Y vector**
* Syntax: `BBCX,<num>,<points>,<yield>,<X1>,…,<Xn>` and `BBCY,<num>,<points>,<yield>,<Y1>,…,<Yn>`
* Each sets the full X (or Y) vector, together with `points` and `yield`.
* Creates the curve if it does not exist.
* If a later BBCX or BBCY gives a different `points` or `yield`, the last command wins, with a warning.
* Before AFWRITE, validate that both vectors exist and have length `points`.

**DELBBC: delete curves**
* Syntax: `DELBBC,<start>,<end>,<stride>`
* Deletes curves from UI memory using the §0.2 convention.

**BBC validity checks** (implementer; the manual has no CHECK error for these):
* npoints ≥ 2.
* 1 ≤ yield ≤ npoints.
* X strictly increasing and > 0.
* Y > 0 (monotonic up to the yield point is expected).
* Point 1 is the cracking point and the origin is implicit.
* Warn when the slope Y1/X1 differs from the elastic stiffness of the referencing panel or spring by more than
  a tolerance. See 05d §3.7: the BBC must be consistent with G·A (shear) or E·I (bending).

### 2.3 BBCGEN: automatic panel BBC generation [ADV][CORE for Option NON]

**Syntax:** `BBCGEN,<Panel>,<ShearModel>,[fc],[fy],[Pn],[Nu],[bre],[bys],[CrackingForceLevel]`

**Behaviour.** BBCGEN builds a BBC for one panel, or for all panels when `Panel = 0` (useful with sub-models to
treat subsets of panels).
* The ultimate shear strength is computed **exactly as the SHEAR command does**, using the selected model.
* Every optional argument left blank is set to **0**.
* The command assumes the panel's group is non-empty and has a uniform material and thickness. Checking this is
  the user's responsibility. Proposal: the implementation should check it and warn.

| Arg | Meaning | Units |
|---|---|---|
| `Panel` | panel number. 0 = all panels with the same parameters. | – |
| `ShearModel` | 1 ACI 318-08, 2 Wood 1990, 3 Barda 1977, 4 Gulec–Whittaker 2009 | – |
| `fc` | f'c, concrete compressive strength | ksi or kN/m² |
| `fy` | reinforcement yield strength | ksi or kN/m² |
| `Pn` | web reinforcement ratio: ρH for ACI, ρV for Wood, Barda and G–W | – |
| `Nu` | axial force | kips or kN |
| `bre`, `bys` | Gulec–Whittaker only. Syntax names `[bre],[bys]`; description names "Fvw", "Fbe". See OQ-8: interpret as boundary-element reinforcement area A_BE and yield stress f_y,BE, so F_BE = A_BE·f_y,BE, as in SHEAR `<Abe>`, `<Fybe>`. | ft² or m²; ksi or kN/m² |
| `CrackingForceLevel` | 0 (default): cracking stress = **3√f'c** (ASCE 4-17 §C.3.3.2). A non-zero value is the cracking/ultimate shear ratio, **allowed only in [0.10, 0.50]**; other non-zero values are rejected. | – |

**BBC construction [CORE].** The BBC always has **22 points**:
1. **Point 1 = cracking point** (γ_cr, V_cr):
   * `CrackingForceLevel = 0`: V_cr = 3√f'c · A_w, in psi units (§2.3.1);
   * otherwise V_cr = CFL · V_u.
   * γ_cr = V_cr / (G·A_w) with G = E/(2(1+ν)) from the panel material. The shear-area choice A_w versus
     (5/6)A_w is not stated (OQ-9).
2. **Points 2–21:** 20 points **equally spaced in strain** from the cracking point to the **yield point**.
   * Point 21 is the yield point, so **yield index = 21**.
   * Yield force = V_u, the computed ultimate shear.
   * The yield strain γ_y and the curve shape between cracking and yield are **not given** (OQ-10).
   * Proposed smooth monotonic shape:
     `V(γ) = V_cr + (V_u − V_cr)·[1 − (1 − ξ)²]` with `ξ = (γ − γ_cr)/(γ_y − γ_cr)`.
     This has zero slope at yield. Check that its initial slope does not exceed G·A_w.
   * Proposed configurable default: γ_y = 0.004.
3. **Point 22 = failure point** at **shear strain = 0.02 (2 %)** and **force = 1.02 × V_u**.

X is shear strain (decimal) and Y is shear force (kips or kN). The BBC type should be 1 (CMS). The destination
curve number is not stated (OQ-10). Proposed: the panel's `bbc` field; if that is 0, use the panel number and
update the panel. This matches PNLGEN, where bbc = panel number.

**Bounds [CORE]:**
* If the ACI 318-08 or Wood 1990 result falls outside its equation's bounds, the bound is used.
* Wood: clamp to [6√f'c·A_w, 10√f'c·A_w].
* ACI: clamp to ≤ 10√f'c·A_w.

#### 2.3.1 SHEAR equations (cross-reference, Manual 9.13.9; transcribed from PDF pp. 291–292) [CORE]

Panel geometry is taken automatically from the panel's shells:
* h_w = panel height;
* l_w = horizontal length;
* t_w = thickness;
* A_w = t_w·l_w (web area).

The unit system is detected from the **gravity value**: ft for British (g ≈ 32.2), m for SI (g ≈ 9.81).
Forces are kips or kN; stresses are ksi or kN/m².

```
Barda et al. 1977 (axial included):
  V = ( 8·sqrt(f'c) − 2.5·sqrt(f'c)·(h_w/l_w) + N_u/(4·l_w·t_w) + ρ_V·f_y ) · t_w · (0.6·l_w)

Wood 1990 (no axial):
  6·sqrt(f'c)·A_w  ≤  V = ρ_V·A_w·f_y / 4  ≤  10·sqrt(f'c)·A_w

ACI 318-08 (no axial):
  V = ( α_c·sqrt(f'c) + ρ_H·f_y ) · A_w  ≤  10·sqrt(f'c)·A_w

Gulec–Whittaker 2009, eq. 6-9 (axial included):
  V = ( 1.5·sqrt(f'c)·A_w + 0.25·F_VW + 0.20·F_BE + 0.40·N_U ) / sqrt(h_w/l_w)
  with F_VW = ρ_V·A_w·f_y ,  F_BE = A_BE·f_y,BE
```

**SHEAR output.** For each panel, six columns:
1. panel number;
2. ACI 318-08 upper bound;
3. Wood 1990 upper bound;
4. Wood 1990 lower bound;
5. Barda 1977;
6. Gulec–Whittaker 2009.

A single-panel call also prints a title line with the equation names.

**Manual guidance:** Barda applies only to barbell walls with heavy flanges and may overestimate the capacity of
nuclear shear walls, so avoid it. Gulec–Whittaker is sensitive to h_w/l_w and grows large for long panels.

**Implementation note (unit trap, OQ-9).** The empirical coefficients (8, 2.5, 6, 10, 1.5, 3, α_c) assume
√f'c in **psi** with the result in psi. With f'c in ksi or kN/m², convert internally:
* US: k·√f'c[psi] psi = k·0.031623·√f'c[ksi] ksi. Areas in ft² must be converted to in² (×144) when the
  stress is in ksi.
* SI: k·√f'c[psi] psi ≈ k·0.08304·√f'c[MPa] MPa = k·2.626·√f'c[kPa] kPa.

The manual does not define α_c. Standard ACI 318-08 §21.9.4.1: α_c = 3.0 for h_w/l_w ≤ 1.5, 2.0 for
h_w/l_w ≥ 2.0, and linear in between.

### 2.4 Panel definition and output support [ADV]

**P: define a wall panel**
* Syntax: `P,<num>,<group>,<bbc>,<disp>,<force>`
* Defines wall panel `num` as SHELL group `group` with backbone curve `bbc`, displacement type `disp` and force
  option `force`. Codes are in §2.1.
* It creates no groups or elements and changes no linear model data.
* **All shells in the group should be coplanar** (vertical plane). WALLFLR creates coplanar groups.

**PNLGEN: generate panels**
* Syntax: `PNLGEN` (no arguments)
* For each shell group whose plane is **vertical**, it creates a panel record. Orientation is taken from the
  **first element of the group**: normal · Z ≈ 0.
* Panel numbers are assigned in ascending order of group number: the first vertical group gets panel 1, and so on.
* Each panel gets `bbc = panel number`, `disp = 1` and `force = 1`.
* The BBC data must then be provided with BBCI + BBCP, BBC (file), BBCX/BBCY or BBCGEN.
* Intended for a panel model, i.e. shells regrouped at least with WALLFLR.

**PDEL: delete panels**
* Syntax: `PDEL,<start>,[end],[stride]`
* Deletes nonlinear panel records only (§0.2). Elements, linear data and other nonlinear data are kept.
* The manual's argument text says "backbone curve" by copy error; it means panels.

**PLIST: list panels**
* Syntax: `PLIST,[start],[end],[stride]`
* Lists panel records. Defaults: 1, largest panel number, 1.
* The index table also describes it as "check". Proposal: show for each panel the group, element count,
  coplanarity residual, BBC presence, and whether disp and force are supported.

**NONLINMOTDISP: panel corner nodes to output requests [CORE for Option NON]**
* Syntax: `NONLINMOTDISP` (no arguments)
* Finds the **4 corner nodes** of every panel defined with P and **adds them to the MOTION and RELDISP nodal
  output request lists**. Ch. 6 says to run it **before AFWRITE**.
* **Algorithm (node-connection counting):** within the panel group, count how many elements use each node.
  * In a structured quad mesh, corners are used by exactly **1** element, edge nodes by 2 and interior nodes
    by 4.
  * Select the nodes with count 1. If the result is not exactly 4, report a failure for that panel.
  * The count fails for triangular meshes and concave or non-rectangular panels. The user must check the panels
    visually.
  * Guidance: keep panels rectangular, or close to rectangular, with low aspect ratio for shear-dominated panels.
* Implementation:
  * do not add a node that is already in a list (that would trigger Error 66);
  * order the corners as BL, BR, TR, TL in panel-local axes and store that order for the NONLINEAR module.

### 2.5 Panel-model construction commands [ADV][UI]

**WALLFLR: separate walls and floors**
* Syntax: `WALLFLR`
* Works on the **active model**:
  1. **Deletes all non-shell elements.** So run it on a copy or sub-model of the original.
  2. Groups shells by a **coplanarity test**: shells with the same plane (normal and offset within tolerance;
     OQ-12). Any coplanar set with **≥ 5 elements** becomes a new group.
  3. Shells that cannot form a wall or floor of at least 5 coplanar elements go to **group 1**.
* Each new group gets a title through `GTIT`, based on orientation: parallel to a global plane (wall along X,
  wall along Y, floor) or **oblique**. The exact title strings are not given.
* After WALLFLR the model is a **panel model**. It can be refined further with PANELIZE.

**PANELIZE: split at intersections**
* Syntax: `PANELIZE`
* Use only on a panel model made by WALLFLR.
* Finds where the wall and floor groups **intersect** and subdivides each group along those intersection lines.
* Proposed algorithm:
  1. For each group, the cut lines are its intersections with other groups' planes, restricted to shared nodes
     and edges.
  2. Remove element adjacency across edges lying on a cut line.
  3. Each resulting connected component becomes a group. The first keeps the original number; the others are
     appended.
* It may produce too many groups. Combine them with `MERGEGROUP` (same element type; then `GCOM` to compress
  group numbers). Use `SPLITGROUP` to split along a wall that does not intersect the group.

**EDGE: split a panel at its edges**
* Syntax: `EDGE,<panel>,[X],[Y],[Z]`
* `panel` is a **group number**, despite the name.
* Splits a panel group using edge detection: the external panel edges plus the edges of every opening.
* Use it after WALLFLR. New groups are appended to the end of the model.
* Flags X, Y, Z: `0` = use edges parallel to that global axis (default); `1` = ignore them.
* Proposed algorithm, which reproduces the manual's Fig. 1.5 example:
  1. Boundary edges = element edges used by exactly one element of the group. These include the outer perimeter
     and the opening perimeters.
  2. Classify each boundary edge as parallel to X, Y or Z (angular tolerance). Drop the classes whose flag is 1.
     Oblique edges are always kept (OQ-13).
  3. Each kept edge defines an infinite **cut line** in the panel plane. Merge collinear lines within a
     tolerance.
  4. Assign each element to a cell by the side of each cut line its centroid lies on. Non-empty cells become
     groups: the first keeps `panel`, the others are new groups.
* Fig. 1.5 example (Ch. 1):
  * `EDGE,1,0,0,1` ignores vertical (Z) edges, so only the horizontal top and bottom edges of the openings cut
    the wall. The result is 3 horizontal strips: group 1 plus new groups 2 and 3.
  * `EDGE,2` on the strip with the openings then cuts along the vertical opening edges. The result is the piers:
    groups 2, 4 and 5. The opening cells are empty and are dropped.
  * A BBC must then be defined per sub-panel.

**EDGEMODEL: EDGE on every wall group**
* Syntax: `EDGEMODEL,[x],[y],[z]`
* Applies `EDGE,<g>,x,y,z` to every **wall** (vertical) group. Flags as for EDGE; default 0 = use edges.

**UNIPNL: one group per element**
* Syntax: `UNIPNL,<group>`
* Gives each element of `group` its own group. Only group membership and element numbers change.
* Use it for panel models, for example curved walls such as a containment, where each shell becomes its own
  panel. See 05d: adjacent panels should differ in angle by less than 15°.
* Proposal: element 1 stays in `group`, and the other elements go to new appended groups, each numbered element 1.

**GROUPMAT (cross-reference, 9.7.18)**
* Creates **one material per group**, copied from the material of the first element in the group.
* The old material list is deleted.
* Required before DGRDFLR, and in general so that NONLINEAR can update each panel's E independently.

**DGRDFLR: scale floor-panel stiffness**
* Syntax: `DGRDFLR,<scale>`
* Builds the list of materials used by **Z-parallel floor panels**, i.e. horizontal shell groups with normal
  ∥ Z, in a panel model.
* Multiplies their **Young's modulus E** by `scale`.
* Assumes WALLFLR and GROUPMAT have been run. Implementation: warn if a floor material is also used by a
  non-floor group.

**MERGEPANEL: put the panel model back into the original model**
* Syntax: `MERGEPANEL,<Panel>`
* Steps for the user:
  1. `ACTM` the **original** model.
  2. **Delete its shell groups**, because the panel model replaces them.
  3. Run `MERGEPANEL,<panel model number>`.
* The panel model's **groups and materials are appended to the end** of the original's group and material
  lists.
* Proposed node handling:
  * panel-model nodes keep the original numbering, because they came from the same model;
  * verify that coordinates match and reuse those nodes;
  * add any node that is missing, with a warning (OQ-14).
* Element material references are remapped to the new material numbers.

**Recommended panel workflow [UI]** (assembled from Ch. 1, Ch. 6 and §9.17; the order is an implementer
proposal, OQ-15):
1. Make a copy or sub-model of the original model and activate it.
2. `WALLFLR`.
3. `PANELIZE`, then `MERGEGROUP`/`SPLITGROUP` as needed.
4. `EDGEMODEL`/`EDGE` for openings; `UNIPNL` for curved walls.
5. `GROUPMAT`.
6. Optionally `DGRDFLR`.
7. `ACTM` the original model, delete its shell groups, and run `MERGEPANEL`.
8. `PNLGEN` (or `P`).
9. `SHEAR` to review capacities, then `BBCGEN` (or BBC/BBCX/BBCY/BBCP).
10. `EQL`.
11. `S` for nonlinear springs.
12. `NONLINMOTDISP`.
13. `AFWRITE` (CHECK errors 121–128).
14. `NONLINBAT,1`.

### 2.6 Nonlinear springs and pile interface [ADV]

**S: nonlinear spring**
* Syntax: `S,<num>,<group>,<elem>,<bbc>,<disp>,<force>`
* Adds nonlinear properties to an **existing SPRING element** (defined with `E`) in group `group`, element
  `elem`. No elements are added.
* `disp`: 1/2/3 = X/Y/Z translation.
* `force`: must be 4 (GMR) (Error 127).
* Springs with several translational stiffnesses should be **split into groups of 1D springs**, so each DOF has
  its own stiffness and damping.

**SOLIDPILE: separate solid piles**
* Syntax: `SOLIDPILE,<group>,[stiff],[soft],[stiff2]`
* `group` is the pile group and must be **SOLID** elements.
* Steps:
  1. Find the **interface nodes**: nodes of the pile group that are also used by elements outside it.
  2. Duplicate each interface node and reconnect the **pile elements** to the duplicate.
  3. Add springs between each original node and its duplicate.
* Each call adds **4 spring groups**:

  | Spring group | Location | Direction (1D) | Stiffness | Default |
  |---|---|---|---|---|
  | 1 | pile side-wall interface nodes | X | `stiff` | **1.0e7** |
  | 2 | side-wall | Y | `stiff` | 1.0e7 |
  | 3 | side-wall | Z | `soft` | **10** |
  | 4 | pile bottom (tip) interface nodes | Z | `stiff2` | **1.0e7** |

* **All springs have 4 % damping** (the SC `damp` field = 0.04).
* The soft vertical side-wall springs allow shaft slip. These springs are the intended targets for `S` with
  GMR to model a nonlinear pile–soil interface.
* Classification of nodes on the bottom perimeter (side-wall vs bottom) is not specified (OQ-16).

**BEAMPILE (not applicable in this version)**
* Syntax: `BEAMPILE`
* Same idea for **beam** piles, i.e. BEAMS with **ETYPE = 2**.
* Duplicates the beam nodes and reconnects the beams to the duplicates.
* Adds **one** spring group per call with translational stiffness **1.0e7** and 4 % damping.

**DCOUPLEBEAM (not usable in this version)**
* Syntax: `DCOUPLEBEAM`
* At beam ends and beam intersections perpendicular to the Z axis, found from connectivity, it adds a new node
  and a spring of stiffness **1.0e7 in all directions**.
* The springs go in a new group, **each with its own spring constant**, so each can be modified independently.

### 2.7 Nonlinear beams (not usable in this version) [ADV]

**B: nonlinear beam**
* Syntax: `B,<num>,<group>,<spgroup>,<bbc>,<force>,<end1>,<end2>`
* Links existing beam element(s) to a BBC. Arguments:
  * `num`: beam number;
  * `group`: beam group;
  * `spgroup`: spring group at the beam end;
  * `bbc`: backbone curve;
  * `force`: force option (dialog: 4 = General Masing Rule);
  * `end1`, `end2`: elements at the two beam ends.
* No elements, materials or properties change.

**DELBM: delete nonlinear beams**
* Syntax: `DELBM,<start>,[end],[stride]`
* Removes the nonlinear-beam flags only, not elements (§0.2).

### 2.8 Global nonlinear options and batch [ADV]

**EQL: nonlinear options**
* Syntax: `EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>`
* Sets the header of the NONLINEAR module input file (`*.eql`), in the same way that HOUSE, ANALYS and STRESS set
  their module headers. Same fields as the Options→Analysis **NONLINEAR** tab (05d §3.8).

| Arg | Meaning | Notes |
|---|---|---|
| `disp` | Equivalent-linear displacement factor (EDF) | typical 0.7–0.9, best ≈ 0.80 |
| `NonLinOpts` | Nonlinear options: use panels, springs, beams | encoding not documented (05d OQ-N10) |
| `dampCutoff` | damping cut-off value (%) | ASCE 4 / NRC limits |
| `dampScale` | damping scale factor | |
| `ElasicD` [sic] | include elastic damping | 0 = don't include, 1 = include |

**NONLINBAT: batch file**
* Syntax: `NONLINBAT,<Sel>`
* Writes a **generic batch run file** for the Option NON analysis, using the active model name and the module
  executable locations configured in the UI.
* `Sel`: `0` = single direction; `1` = three directions, which includes the **COMB_XYZ_THD** step. The user
  must provide `COMB_XYZ_THD.inp`.
* The user can run the file as is or edit it.
* File name and content are not specified (OQ-17). Proposed content (see 05d §3.2):
  * an initial elastic pass: SITE, POINT, HOUSE, ANALYS, MOTION, RELDISP [×3 directions], COMB_XYZ_THD,
    NONLINEAR;
  * then an iteration loop: HOUSE, ANALYS (restart), MOTION, RELDISP, COMB_XYZ_THD, NONLINEAR, with per-iteration
    file renaming.

### 2.9 Nonlinear soil (SOIL module, time-domain hyperbolic model) [ADV]

**NLSLAYER: nonlinear soil layer**
* Syntax: `NLSLAYER,<Num>,[curvefit],[B],[S],[refStrain],[Vis]`
* Defines a **nonlinear soil layer data set** for the SOIL module. **Every argument defaults to 0.**
* `curvefit = 1`: the parameters are fitted automatically and the others are ignored.
* `curvefit = 0`: the user must give all curve parameters.

The model is a modified hyperbolic stress–strain law with viscosity. Transcribed from PDF p. 323 (printed 321):

```
            G0 · γ
  τ  =  ───────────────────  +  η · γ̇
         1 + β · (γ / γr)^s
```

Symbols and manual typical values:
* τ: shear stress.
* γ: shear strain; γ̇: shear strain rate.
* G0: elastic (small-strain) shear modulus. Proposed: ρ·Vs² of the layer.
* γr: reference shear strain, typically about **0.03 %**.
* β: constant, typically about **1 to 1.4**.
* s: exponent, typically about **0.7 to 0.9**.
* η: soil viscosity.

The implied backbone is G/G0 = 1/(1 + β(γ/γr)^s). This is the MKZ form used by DEEPSOIL; the manual cites
DEEPSOIL for SOIL-NON.

| Arg | Meaning |
|---|---|
| `Num` | soil layer number |
| `curvefit` | 0 = no curve fit; 1 = fit β, s, γr automatically (to the layer's G/Gmax curve; proposal: least squares on log-strain points) |
| `B` | β |
| `S` | s exponent |
| `refStrain` | γr. Whether it is entered in % or decimal is not stated (OQ-18). |
| `Vis` | η, viscosity (used with damping type 2) |

Unloading and reloading rules (Masing or non-Masing) are not stated (OQ-18).

**NLSOIL: global nonlinear soil options**
* Syntax: `NLSOIL,<Opt>,<NSTimeSunInc>,<DispConv>,<ForceConv>,<EqualIt>,<BedInt>,<NLDampType>,<MMmult>,<SMmult>`
* Global options for the SOIL module. It does **not** change the linear (equivalent-linear) soil model.

| Arg | Meaning | Values |
|---|---|---|
| `Opt` | use the nonlinear soil time-domain analysis | 0 no, 1 yes |
| `NSTimeSunInc` [sic] | sub-increments per time step | int ≥ 1 (implementer) |
| `DispConv` | displacement convergence error (tolerance) | real > 0 |
| `ForceConv` | force convergence error (tolerance) | real > 0 |
| `EqualIt` | maximum equilibrium iterations | int ≥ 1 |
| `BedInt` | bedrock interface | 0 Rigid; 1 Viscoelastic (**currently disabled**) |
| `NLDampType` | small-strain damping type | 1 frequency-independent, 2 visco-elastic (uses Vis), 3 Rayleigh (uses MMmult, SMmult) (codes from 05a §6.6) |
| `MMmult` | mass matrix multiplier α in C = αM + βK | real |
| `SMmult` | stiffness matrix multiplier β | real |

CHECK **Error 125**: there are fewer NLSLAYER definitions than soil layers.

The time-integration scheme (Newmark, Newton–Raphson on the 1D shear column) is not described here (OQ-19).
Validation is V&V Problem 49, against DEEPSOIL.

**SOILREDEF: layers from a 2D soil model (not validated)**
* Syntax: `SOILREDEF,<soil>,<dir>`
* Redefines the **excavation-volume layering** of the active 3D model using a 2D soil model.
* Both models must be in UI memory, with the 3D model active. **Only the layer list and the excavation-volume
  layer assignments** of the 3D model change.
* The 2D model must lie in the **global X–Z plane**.
  * `dir = 0`: the 2D x coordinate maps to 3D **X**.
  * `dir = 1`: the 2D model is moved from the X–Z plane to the **Y–Z** plane, so the 2D x maps to 3D Y.
* Every excavation-volume element must be **explicitly** flagged with `ETYPE`/`ETYPEGEN` (type 2). Only
  global-axis alignment is supported; use `ROTATE` first if needed.
* Proposed algorithm: for each excavation SOLID, project its centroid to (x or y, z), find the 2D soil element
  that contains it, and assign that element's soil layer and properties. Build the new layer list from the
  distinct assignments.

### 2.10 Manual text errors in §9.17 to correct in the re-implementation

| Item | Manual | Correct behaviour (inferred) |
|---|---|---|
| DELNLS (9.17.12) | Syntax line shows `DELBBC,…`; text "Delete backbone curves" | `DELNLS,<start>,[end],[stride]` deletes NLSLAYER data sets. It does **not** touch L, SOIL or any other linear soil data. |
| DELSPR (9.17.13) | Syntax line shows `DELNLS,…`; its text (pp. 319) describes NLSLAYER deletion | `DELSPR,<start>,[end],[stride]` deletes nonlinear spring records made by `S`. The springs stay as linear elements. |
| PDEL args | "backbone curve" | panels |
| EDGE / EDGEMODEL flag lists | second option printed without "1" | `1` = ignore edges |
| BBCGEN args | syntax `[bre],[bys]` vs text "Fvw", "Fbe"; `CrackingForceLevel` vs `CrackForceLevel` | see OQ-8 |
| Ch. 6 names | `PANELGEN`, `NONLINMODISP` | aliases of `PNLGEN`, `NONLINMOTDISP` (accept both) |

---

## 3. Binary Database Commands (Manual 9.18) [IO][UI]

In the manual's index table, every row points to "Section 9.17.29". This is a copy error; the real sections are
9.18.1–9.18.15.

### 3.0 Binary databases: producers, names and memory model [IO]

| DB code | Content | Produced by | Enabled by | Default file name(s) (Ch. 3 / Ch. 6) |
|---|---|---|---|---|
| **ACC** | nodal acceleration histories | MOTION | `BINOUT,<mot>=1` (MOTION tab "Save Binary Database") | `Modelname_ACC.bin` (one per input direction run) |
| **THS** | element stress, force and moment histories | STRESS | `BINOUT,,<str>=1` | `Modelname_STRESS.bin` |
| **DISP** | nodal relative-displacement histories | RELDISP | `BINOUT,,,<reldisp>=2` (THD binary) | per response direction: `Modelname_TR_X_THD.bin`, `Modelname_TR_Y_THD.bin`, `Modelname_TR_Z_THD.bin`. Combined: `Modelname_THD.bin` |

* **Only one database of each type** (ACC, DISP, THS) can be loaded in UI memory at a time. Loading another of
  the same type **replaces** it.
* The binary layout is **not documented** (OQ-20). A re-implementation must define its own container while
  keeping the file names. A proposed schema is in §3.6.
* Databases are "compressed" according to Ch. 3. Stress databases may also be in **ANSYS format**
  (`LOADTHSDB` sel = 1), for Option A/AA.

### 3.1 BINOUT: binary output flags [IO]

**Syntax:** `BINOUT,[mot],[str],[reldisp]`

Sets the binary-output flags written into the SSI module inputs. **A blank argument leaves its flag unchanged.**

| Arg | 0 | 1 | 2 |
|---|---|---|---|
| `mot` | MOTION writes no ACC database | write the ACC database | – |
| `str` | STRESS writes no THS database | write the THS database | – |
| `reldisp` | RELDISP writes no database | TFD binary (**not used in this version**) | write the **THD** (time-history) database |

### 3.2 Load and delete [UI]

| Command | Syntax | Behaviour |
|---|---|---|
| LOADACCDB (heading misprinted as "LOADACCDBANI") | `LOADACCDB,<file>` | Load a MOTION ACC database (full path). Replaces the loaded ACC DB. |
| LOADDISPDB | `LOADDISPDB,<file>` | Load a RELDISP DISP database. Replaces the loaded DISP DB. |
| LOADTHSDB | `LOADTHSDB,<file>,[sel]` | Load a STRESS THS database. `sel`: **0 = ACS SASSI format (default)**, 1 = ANSYS® format. |
| DELDB | `DELDB,[sel]` | Remove database(s) from UI memory. `sel`: **ALL (default)**, ACC, THS, DISP. Files on disk are not deleted. |

### 3.3 Combine directional databases [IO][CORE for post-processing]

All four commands have the form `<CMD>,<Xfile>,<Yfile>,<Zfile>,<Comb>`, with full paths.

| Command | Inputs | Output |
|---|---|---|
| COMBACCDB | ACC databases from the X-, Y- and Z-input runs | one combined ACC DB |
| COMBDISPDB | DISP databases for the X, Y and Z **input** directions (each already made by COMBDISPDIR) | one combined DISP DB |
| COMBDISPDIR | **RELDISP only.** For one input direction, the three **response-direction** databases: the principal component (for example X–X) and the coupling responses (for example Y–X, Z–X) | one "single input direction" DISP DB. **Use it before COMBDISPDB**, once per input direction. |
| COMBTHSDB | THS databases for the X, Y and Z inputs | one combined THS DB |

**RELDISP 3-direction workflow (Ch. 3):**
1. RELDISP is run **three times per input direction**, giving `…TR_X_THD.bin`, `…TR_Y_THD.bin` and
   `…TR_Z_THD.bin`.
2. `COMBDISPDIR` merges them into `Modelname_THD.bin` for that input direction.
3. `COMBDISPDB` combines the three input directions.

**Combination rule (not stated; OQ-21).** For COMBACCDB, COMBDISPDB and COMBTHSDB, the proposal is
**algebraic summation, time step by time step**, of co-directional responses. This gives the response to the
three simultaneous input components and is consistent with ASCE 4 time-history combination.
* Inputs must share the node/element sets, number of steps and time step. Otherwise it is an error.
* COMBDISPDIR does **not** sum. It assembles response components from separate files into one record per node.

### 3.4 Frame and table output [IO][UI]

**BINFRAMEOUT: ASCII frames from a database**
* Syntax: `BINFRAMEOUT,<db>,<frame>,<TS>,[Split],<dir>`
* Writes ASCII frame file(s) in the **UI frame format** for one database at one time step. These are the same
  frames MOTION, RELDISP and STRESS write (Table 3.2), so FRAMECOMBIN and MODFRAMES can use them.

| Arg | Values |
|---|---|
| `db` | `THS` stress, `DISP` displacement, `ACC` acceleration |
| `frame` | frame number (integer order) **or** frame time. `frame = −1` with `TS ≤ 0` writes the **maximum-component frame**. |
| `TS` | `TS ≤ 0`: `frame` is an integer frame number. `TS > 0`: `frame` is a time and TS is the database time step. Proposed: index = round(time/TS), plus 1 if frames are 1-based (OQ-22). |
| `Split` | 0 = one file; 1 = translations and rotations in **separate files** (ACC and DISP only; ignored for THS) |
| `dir` | output directory |

Frame names follow Table 3.2, for example `ACC_<time>_<fnum>`, `THD_<time>_<fnum>` and
`stress_<time>_<fnum>_<comp>`. Maximum frames follow the `*_ABS_MAX_*` pattern.

**MAXDBFRAME: maximum frame**
* Syntax: `MAXDBFRAME,<Type>,[dir]`
* Finds the **maximum value** of each component in the loaded database (`THS` | `DISP` | `ACC`) and stores it as
  a **single-frame animation**: binary file(s) in `dir` (default: current working directory), registered in the
  animation database.
* Proposed: per DOF or component, the max |x(t)| with its sign kept, as in BINSTRTBL.

**BINSTRTBL: stress table to CSV**
* Syntax: `BINSTRTBL,<group>,<EVar>,[step],<file>`
* Writes a **CSV** table from the loaded THS database.
  * Line 1: column labels.
  * Each following line: `Group, Element, comp1, comp2, …` for one element.
* `EVar` is the name of a UI **variable** (`VAR,<Name>,<X1>…` or `LOADVAR`, §9.15) that must already hold the
  element numbers, all positive integers.
* `step`:
  * blank or **−1 (default)**: for each component, the value with the largest absolute value, written **with its
    sign**;
  * otherwise an integer time-step index. If it is greater than the last step, an error is posted in the command
    window.
* Elements not in the database are written with all components **0**.
* `file`: full path.

### 3.5 Animation creation [UI]

| Command | Syntax | Source DB | Result |
|---|---|---|---|
| ACCDBANI (heading misprinted "ACCANIDB") | `ACCDBANI,<dir>,[label]` | loaded ACC | SSI model animation (deformed/bubble type) |
| DISPDBANI | `DISPDBANI,<dir>,[label]` | loaded **DISP**. The manual text says "acceleration" by copy error. | deformed-shape animation |
| THSDBANI | `THSDBANI,<dir>,[label]` (text calls it `desc`) | loaded THS | element **contour-plot** animation |

* `dir`: buffer directory for the animation frame files. They are in the UI binary animation format and should
  not be edited.
* `label`: description stored in the animation database (`SASSIani.xml`, Ch. 7). It is used to pick the
  animation when reloading.
* These appear in the animation list as the "(binary)" items.

### 3.6 Proposed database schema (implementation decision, not from the manual) [IO]

Use one self-describing container per file. HDF5 (`h5py`) or NumPy `.npz` works well. Keep the manual's file
names. Contents:

```
/meta      : db_type ("ACC"|"DISP"|"THS"), model_name, input_dir ("X"|"Y"|"Z"|"XYZ"),
             response_dirs, dt, nsteps, units, format_version, source_module
/time      : float64[nsteps]
ACC/DISP   : /nodes int32[nn]; /data float32[nn, ndof(6: TX,TY,TZ,RX,RY,RZ), nsteps]
THS        : /groups int32[]; per group: /etype str, /elements int32[ne], /components str[nc],
             /data float32[ne, nc, nsteps]
```

Add an exporter for the "ANSYS format" stress database that LOADTHSDB sel=1 reads. Its layout is also unknown
(OQ-20).

---

## 4. Thick Shell Commands (Manual 9.19) [IO]

The index table again points to "Section 9.17.29" by copy error.

### 4.1 THSHLSTR: TSHELL face stress/strain output flag [IO][CORE for STRESS output]

**Syntax:** `THSHLSTR,<flag>`

Affects **only the values computed by the STRESS module** for TSHELL (Mindlin–Reissner) elements.

**flag = 0 (default):** write only the **8 basic element components**, in the local element system:

| # | Component | Symbol | Units |
|---|---|---|---|
| 1 | Force XX | Nxx | F/L |
| 2 | Force YY | Nyy | F/L |
| 3 | Force XY | Nxy | F/L |
| 4 | Transverse shear XZ | Qxz | F/L |
| 5 | Transverse shear YZ | Qyz | F/L |
| 6 | Moment XX | Mxx | F·L/L |
| 7 | Moment YY | Myy | F·L/L |
| 8 | Moment XY (printed "ZZ-Direction (Mxy)") | Mxy | F·L/L |

**flag = 1:** also write **face (top/bottom) stresses and strains**, computed from the **maximum values of the
basic components** (not from time-history combination).

The face formulas are not printed. Proposed, using thickness t and the local system:

```
σxx(±) = Nxx/t ± 6·Mxx/t² ;  σyy(±) = Nyy/t ± 6·Myy/t² ;  τxy(±) = Nxy/t ± 6·Mxy/t²
τxz,max = 1.5·Qxz/t , τyz,max = 1.5·Qyz/t  (mid-surface, parabolic; zero at faces)
εxx = (σxx − ν·σyy)/E ; εyy = (σyy − ν·σxx)/E ; γxy = τxy/G ,  G = E/(2(1+ν))
```

How the maximum values are combined into face values (sign permutations, as STRESS does for SHELL; see 05d
§2) is OQ-23.

Note: Ch. 3 warns that **no text frame files are generated for TSHELL** elements.

### 4.2 THSHLSMH: smoothing of TSHELL transverse shear (not validated) [UI]

**Syntax:** `THSHLSMH,<passes>,[type],[workdir]`

* A "smear" function applied to TSHELL shear output components, mainly the local **transverse shear forces**
  (Qxz, Qyz). Element-by-element values are smoothed using the adjacent elements.
* Arguments:
  * `passes`: number of adjacent element layers (rings of neighbours) included.
  * `type`: **0 = weighted average with a Parzen blending window (default)**; 1 = simple average (constant
    window).
  * `workdir`: directory for the output frame (default: working directory).
* Proposed algorithm:
  1. Build element adjacency through shared nodes or edges and compute the ring distance k (0 = self) up to
     L = passes.
  2. Smoothed value = Σ w(k)·Q_e / Σ w(k).
  3. For type 1, w = 1.
  4. For type 0, use the Parzen weights with M = L + 1:
     * w(k) = 1 − 6(k/M)² + 6(k/M)³ for k ≤ M/2;
     * w(k) = 2(1 − k/M)³ for M/2 < k ≤ M.
* Restrict smoothing to elements of the same group or plane, so that values are not averaged across wall/slab
  junctions (OQ-24).

---

## 5. CHECK Errors and Warnings (Manual Chapter 10): the model checker [CORE]

### 5.1 Behaviour [CORE][UI]

* `CHECK` checks every enabled module's data and lists errors and warnings in the **Check Errors** window.
  * It **simulates writing the analysis files**, so analysis options must be set first.
  * It can run on a partial model.
* The output is like legacy PREP, with **module headers** that group the messages by module. Proposed line format
  (01_capabilities §): `Errors and Warnings for <MODULE>`, then lines `Error <n> : <text>` and
  `Warning <n> : <text>`.
* `AFWRITE` runs CHECK first. **A module with any error does not have its input file written.** Warnings never
  block.
* Model-level errors (1–43, 83, 124) are attributed to HOUSE. HOUSE then fails, and the modules that depend on
  it are also not runnable. The exact attribution is inferred; see the module column below.
* **Fatal errors:** Errors **42** (no nodes) and **43** (no groups) stop CHECK immediately.
* Placeholders: `<i>` index; `<e>` element; `<g>` group; `<m>` material; `<l>` soil layer; `<p>` property;
  `<n>` node; `<s>` frequency set; `<w>` wave; `<k>` panel or spring.
* Remedy text: most option errors tell the user to "set/correct … by selecting the Options / Analysis command".
  Model errors name the command to use (ECOMPR, KI/KJ, SC, MTDEL/MRDEL, N, FDEL/MMDEL, MXR/MXI/MXM, FREQ, P, M,
  GROUP, INT, D).

**Tolerance policy (implementer, OQ-25).** The manual gives no numeric tolerances. Use relative tolerances based
on the model size L_ref (bounding-box diagonal):
* coincident nodes: |Δx| < 1e-9·L_ref;
* collinear: |a×b| < 1e-9·|a||b|;
* zero area: A < 1e-12·L_ref²;
* warping: the measure in Errors 12 and W3 below.

### 5.2 Error catalogue

The message titles are exact, because they are the nomenclature. Conditions are given as implementable
predicates. The Module column is the header group (inferred from the related option tabs in specs 05a–05d).

| # | Message | Check condition (error when…) | Module | Remedy |
|---|---|---|---|---|
| 1 | Illegal Acceleration of Gravity | gravity g ≤ 0 | SITE/HOUSE (analysis options) | Options/Analysis |
| 2 | Illegal Node for Symmetry Plane/Line `<i>` | a SYMM node of plane/line i is an illegal number (≤ 0) or undefined | HOUSE | redefine via SYMM |
| 3 | Illegal Number of Nodes for Symmetry Plane/Line `<i>` | plane/line i has only one node | HOUSE | SYMM |
| 4 | Nodes from Symmetry Line `<i>` Are Equal | line nodes are the same id or have identical coordinates | HOUSE | SYMM |
| 5 | Nodes from Symmetry Plane `<i>` Are Collinear | the 3 plane nodes are collinear | HOUSE | SYMM |
| 6 | Undefined Element `<e>`, Group `<g>` | gap in element numbering: e < max element number of g and e not defined | HOUSE | define e, or ECOMPR |
| 7 | Element `<e>` from Group `<g>` Has Too Few Nodes | fewer nodes defined than the group type requires | HOUSE | define nodes / change type (GROUP) |
| 8 | Element `<e>`, Group `<g>` Has 0 Length | two-node element with the same node twice or coincident nodes | HOUSE | check nodes |
| 9 | Nodes from Element `<e>`, Group `<g>` Are Collinear | 3- or 4-node element whose nodes are the same, coincident or collinear | HOUSE | check nodes |
| 10 | Element `<e>`, Group `<g>` Has Improper Release Code | BEAMS: for some component c ∈ {P1,P2,P3,M1,M2,M3}, KI[c] = 1 **and** KJ[c] = 1 | HOUSE | KI / KJ |
| 11 | Zero Area for Element `<e>`, Group `<g>` | 4-node element whose area ≈ 0 (some nodes collinear or the same) | HOUSE | check nodes |
| 12 | Warped Element `<e>`, Group `<g>` | 4-node element whose nodes are not in one plane (beyond the error tolerance) | HOUSE | check nodes |
| 13 | Material `<m>` Is not Defined | an element references undefined material m | HOUSE | check indices / ETYPE / define M |
| 14 | Elasticity Modulus from Material `<m>` Is Illegal | E ≤ 0 | HOUSE | M |
| 15 | Poisson Coefficient from Material `<m>` Is Illegal | ν ≤ 0 | HOUSE | M |
| 16 | Specific Weight from Material `<m>` Is Illegal | γ < 0 | HOUSE | M |
| 17 | P-Wave Damping Ratio from Material `<m>` Is Illegal | Dp < 0 | HOUSE | M |
| 18 | S-Wave Damping Ratio from Material `<m>` Is Illegal | Ds < 0 | HOUSE | M |
| 19 | Soil Layer `<l>` Is not Defined | an (excavated-soil) element references undefined layer l | HOUSE | check indices / ETYPE / define L |
| 20 | Thickness from Soil Layer `<l>` Is Illegal | thickness ≤ 0 | HOUSE/SITE | L |
| 21 | Specific Weight from Soil Layer `<l>` Is Illegal | γ < 0 | HOUSE/SITE | L |
| 22 | P-Wave Velocity from Soil Layer `<l>` Is Illegal | Vp < 0 | HOUSE/SITE | L |
| 23 | S-Wave Velocity from Soil Layer `<l>` Is Illegal | Vs < 0 | HOUSE/SITE | L |
| 24 | P-Wave Damping Ratio from Soil Layer `<l>` Is Illegal | Dp < 0 | HOUSE/SITE | L |
| 25 | S-Wave Damping Ratio from Soil Layer `<l>` Is Illegal | Ds < 0 | HOUSE/SITE | L |
| 26 | Property `<p>` Is not Defined | a BEAMS element references undefined property p | HOUSE | define BEAM property |
| 27 | Axial Area from Property `<p>` Is Illegal | A ≤ 0 | HOUSE | property |
| 28 | Shear Area 2 from Property `<p>` Is Illegal | As2 < 0 | HOUSE | property |
| 29 | Shear Area 3 from Property `<p>` Is Illegal | As3 < 0 | HOUSE | property |
| 30 | Torsion Inertia Moment from Property `<p>` Is Illegal | J ≤ 0 | HOUSE | property |
| 31 | Flexural Inertia Moment 2 from Property `<p>` Is Illegal | I2 < 0 | HOUSE | property |
| 32 | Flexural Inertia Moment 3 from Property `<p>` Is Illegal | I3 < 0 | HOUSE | property |
| 33 | Spring Property `<p>` Is not Defined | a SPRING element references undefined spring property p | HOUSE | SC |
| 34 | Spring Constant X from Spring Property `<p>` Is Illegal | kx < 0 | HOUSE | SC |
| 35 | Spring Constant Y from Spring Property `<p>` Is Illegal | ky < 0 | HOUSE | SC |
| 36 | Spring Constant Z from Spring Property `<p>` Is Illegal | kz < 0 | HOUSE | SC |
| 37 | Spring Constant XX from Spring Property `<p>` Is Illegal | kxx < 0 | HOUSE | SC |
| 38 | Spring Constant YY from Spring Property `<p>` Is Illegal | kyy < 0 | HOUSE | SC |
| 39 | Spring Constant ZZ from Spring Property `<p>` Is Illegal | kzz < 0 | HOUSE | SC |
| 40 | Mass out of Defined Nodes Range | a translational or rotational mass is on a node number outside [1, max node] | HOUSE | define nodes or MTDEL / MRDEL |
| 41 | Node `<n>` Is Not Defined - Group `<g>` Element `<e>` | element e of group g uses undefined node n | HOUSE | define N |
| 42 | No Nodes Defined | node count = 0. **FATAL: CHECK stops.** | HOUSE | – |
| 43 | No Groups Defined | group count = 0. **FATAL: CHECK stops.** | HOUSE | – |
| 44 | Frequency Set `<s>` Is Not Defined | the frequency set selected in the analysis options does not exist | SITE (also ANALYS/RELDISP) | FREQ / options |
| 45 | Mode 1 and Mode 2 Are Both Deselected | SITE Mode 1 and Mode 2 both off | SITE | options |
| 46 | No Top Layers | SITE top-layer list empty | SITE | TOPL / options |
| 47 | Illegal Number of Layers for Halfspace Simulation | NL ≠ 0 and NL ∉ [4, 20] | SITE | options |
| 48 | Illegal Frequency Step | Δf < 0 | SITE | options |
| 49 | Illegal Time Step of Control Motion | Δt < 0 | SITE (shared Δt) | options |
| 50 | Illegal Number of Values for Fourier Transform | NFFT < 0 | SITE (shared NFFT) | options |
| 51 | All Wave Fields Are Deselected | no wave type (P, SV, SH, R, L) selected | SITE | options |
| 52 | Illegal Incident Angle of Wave `<w>` | angle ∉ (0, 360) degrees. Implement [0, 360), because 0 = vertical incidence is the normal default (OQ-26). | SITE | options |
| 53 | Illegal Value for Frequency `<i>` | ratio-curve frequency i ≤ 0 | SITE | options |
| 54 | Illegal Value for Wave `<w>` Ratio at Frequency `<i>` | ratio ≤ 0 or > 1 | SITE | options |
| 55 | Illegal Sum of Wave Ratios at Frequency `<i>` | Σ ratios of the selected waves ≠ 1 (tolerance, e.g. 1e-6) | SITE | options |
| 56 | Illegal Last Layer Number in Near Field Zone | last near-field (embedment) layer < 0 | POINT | options |
| 57 | Illegal Radius of Central Zone | radius ≤ 0 | POINT | options |
| 58 | Illegal Coherence Parameter | coherence parameter < 0.1 (Luco–Wong γ) | HOUSE (incoherency) | options |
| 59 | Illegal Mean Soil Shear Wave Velocity | mean Vs ≤ 0 | HOUSE (incoherency) | options |
| 60 | Illegal Number of Mesh Points / Embedment Level | value ≤ 0 (when incoherency is active; OQ-27) | HOUSE (incoherency) | options |
| 61 | No Forces Defined | FORCE module active and no F / MM defined | FORCE | F, MM |
| 62 | Node for Force / Moment `<i>` is not defined | force or moment i applied to an undefined node | FORCE | N, or FDEL / MMDEL |
| 63 | Illegal Coordinate Transformation Angle | angle ∉ (0, 360). Implement [0, 360) (OQ-26). | ANALYS | options |
| 64 | No Nodal Output Request | no node lists defined for output | MOTION (and RELDISP) | options |
| 65 | Illegal Nodal Output Request: `<n>` | a node list contains an illegal or undefined node n, explicitly or within a range | MOTION / RELDISP | options |
| 66 | Nodal Output Request Defined More Than Once: `<n>` | node n appears more than once over all lists, ranges included | MOTION / RELDISP | options |
| 67 | Illegal Output Time History Step | output TH step < 0 | MOTION / STRESS | options |
| 68 | Illegal Total Duration To Be Plotted | duration < 0 | MOTION | options |
| 69 | Illegal First Frequency for RS Analysis | f_first < 0 | MOTION | options |
| 70 | Illegal Last Frequency for RS Analysis | f_last < 0 | MOTION | options |
| 71 | Illegal Number of Frequency Steps For RS Analysis | n_steps < 0 | MOTION | options |
| 72 | Illegal Damping Ratio For RS Analysis | a damping ratio in the list ∉ (0, 1) | MOTION | options |
| 73 | Acceleration Time History File Does Not Exist | THFILE path missing or misspelled | MOTION / STRESS / SOIL | options / THFILE |
| 74 | Illegal First Record Number | first record < 0 or > number of records in the file | MOTION / STRESS | options |
| 75 | Illegal Last Record Number | last record < 0 | MOTION / STRESS | options |
| 76 | First Record Number Larger Than Last Record Number | first > last | MOTION / STRESS | options |
| 77 | Multiplication Factor and Maximum Value Of Time History Are Both Zero | MF = 0 and Max = 0 | MOTION / STRESS | options |
| 78 | Multiplication Factor and Maximum Value Of Time History Are Both Non-Zero | MF ≠ 0 and Max ≠ 0 (exactly one must be non-zero) | MOTION / STRESS | options |
| 79 | No Element Output Request | no element lists defined | STRESS | options |
| 80 | Illegal Group For Output Request: `<g>` | an element list names an illegal or nonexistent group g | STRESS | options |
| 81 | Illegal Element Output Request: `<e>`, Group `<g>` | an element list for group g contains an illegal element e (≤ 0 or > count), explicitly or within a range | STRESS | options |
| 82 | Element Output Request Defined More Than Once: `<e>`, Group `<g>` | (g, e) appears more than once | STRESS | options |
| 83 | Matrix Property `<p>` Is not Defined | a GENERAL element references undefined matrix property p | HOUSE | MXR, MXI, MXM |
| 84 | No RS Input Files Specified | all three spectrum input files blank | EQUAKE | options |
| 85 | RS Input File `<i>` Does Not Exist | spectrum input file i (non-blank) not found | EQUAKE | options |
| 86 | Invalid RS Output File `<i>` | spectrum output file i blank (for an active component) | EQUAKE | options |
| 87 | Invalid Acceleration Output File `<i>` | acceleration output file i blank (for an active component) | EQUAKE | options |
| 88 | Invalid Acceleration Input File `<i>` | acceleration-input option selected and file i blank | EQUAKE | set it, or deselect the option |
| 89 | Number of Frequencies Does Not Match RS Input File `<i>` | number of records in spectrum file i ≠ number of frequencies | EQUAKE | options |
| 90 | Illegal Initial Random Number | seed ≤ 0 | EQUAKE | options |
| 91 | Illegal Number of Frequencies | n_freq ≤ 0 | EQUAKE | options |
| 92 | Illegal Duration | total duration ≤ 0 | EQUAKE | options |
| 93 | No Correlation Factors Defined | correlation option on and no factors | EQUAKE | define or deselect |
| 94 | Illegal Correlation Factor | a factor > 1 (proposal: \|ρ\| > 1) | EQUAKE | options |
| 95 | No Dynamic Soil Properties Assigned | no dynamic property assigned to any SOIL sublayer | SOIL | options |
| 96 | Too Many Dynamic Soil Properties | more than **15** distinct dynamic properties used by the sublayers | SOIL | options |
| 97 | Dynamic property `<p>` has no shear modulus curve | the G/Gmax–γ curve of p is empty | SOIL | options |
| 98 | Dynamic property `<p>` has illegal shear modulus values | some G/Gmax < 0 or > 1 | SOIL | options |
| 99 | Dynamic property `<p>` has no damping curve | the D–γ curve of p is empty | SOIL | options |
| 100 | Number of Acceleration Values Is Illegal | n_values ≤ 0 | SOIL | options |
| 101 | Cut-Off Frequency Is Illegal | f_cut < 0 | SOIL | options |
| 102 | Illegal Reading Format | the reading-format specification cannot be parsed or used | SOIL | options |
| 103 | Illegal Number of Header Lines | header lines < 0 | SOIL | options |
| 104 | Illegal Control Layer Number | control layer ∉ [1, number of layers] (inferred) | SOIL | options |
| 105 | Illegal Number of Iterations | iterations < 0 | SOIL | options |
| 106 | Illegal Strain Ratio | equivalent/max strain ratio ∉ (0, 1) | SOIL | options |
| 107 | No Damping Ratios Defined | the RS damping list is empty | SOIL | DAMP / options |
| 108 | Illegal Multiplier for Acceleration of Gravity | multiplier ≤ 0 | SOIL | options |
| 109 | Illegal Second Layer Number for Layer `<i>` | spectral-amplification second layer of layer i ∉ [1, number of layers] (inferred) | SOIL | options |
| 110 | Illegal Frequency Step for Layer `<i>` | frequency step ≤ 0 (for an active SSAF request) | SOIL | options |
| 111 | Illegal Number of Smoothings for Layer `<i>` | smoothings < 0 | SOIL | options |
| 112 | Illegal Number of Values to Be Saved for Layer `<i>` | values to save < 0 | SOIL | options |
| 113 | Illegal Apparent Velocity for Line D | apparent velocity ≤ 0 (wave passage on) | HOUSE (incoherency / wave passage) | options |
| 114 | Illegal Directional Coherence Factor | factor < 0 | HOUSE | options |
| 115 | No Multiple Excitation Data Defined | multiple-excitation option on and no data | HOUSE | define or deselect |
| 116 | Illegal First Node Number for Motion `<i>` | first node of motion i ≤ 0, undefined, or not an interaction node (inferred) | HOUSE | options |
| 117 | Illegal Last Node Number for Motion `<i>` | last node illegal (as in 116) or < first node (inferred) | HOUSE | options |
| 118 | Illegal Spectral Amplification Ratio for Motion `<i>` | an amplification factor ∉ [0, 10] | HOUSE | options / AMP |
| 119 | Spectral Amplification Ratios for Motion `<i>` Do Not Match Frequencies | count of ratios ≠ number of frequencies in the selected set | HOUSE | options / AMP / FREQ |
| 120 | Frequency set `<i>` is Empty. Change the active frequency set in the SITE tab of the Analysis Options | the active frequency set (used for relative displacement) has no frequencies | RELDISP (set chosen in SITE) | SITE tab, FREQ |
| 121 | No Panels Specified | panel (NONLINEAR) module enabled and no P records | NONLINEAR | P |
| 122 | Material `<i>` referenced in panel `<k>` does not exist | the material of panel k's group elements is undefined | NONLINEAR | P, or define M |
| 123 | Group `<i>` referenced in panel `<k>` does not exist | panel k's group i does not exist | NONLINEAR | P, or GROUP |
| 124 | Node `<i>` is a fixed interaction node | node i is an interaction node **and** fully fixed (proposal: all translational DOFs fixed; OQ-28) | HOUSE/ANALYS | INT or D |
| 125 | Not enough nonlinear soil properties have been defined for the Nonlinear soil model. | NLSOIL Opt = 1 and the number of NLSLAYER sets < the number of soil layers | SOIL | NLSLAYER |
| 126 | Force Option `<i>` not supported for panel `<k>` in this version | panel force ≠ 1 | NONLINEAR | set 1 |
| 127 | Force Option `<i>` not supported for spring `<k>` | spring force ≠ 4 | NONLINEAR | set 4 |
| 128 | Displacement Option `<i>` not supported for panel `<k>` | panel disp ≠ 1 | NONLINEAR | set 1 |

**Element-specific implementation notes:**
* **Error 7 node counts:** BEAMS need I, J (+K orientation node); SHELL/TSHELL 4 (a triangle is a repeated or
  degenerate 4th node); PLANE 4; SOLID 8 (prisms by repeated nodes); SPRING 2; GENERAL 2 or 3.
* **Error 8 must not be raised for SPRING elements.** Zero-length springs between coincident nodes are
  intentional; FILLPOOL and SOLIDPILE create them. Apply Error 8 to BEAMS (I–J) and to 2-node GENERAL elements
  only (OQ-29).
* **Error 9 for BEAMS:** a K node collinear with I–J leaves the local axes undefined. Treat it as a "three-node"
  collinearity error.
* **Error 11 and triangles:** remove consecutive duplicate nodes before computing the area, so that legitimate
  triangles entered as 4-node shells with a repeated node do not trigger Error 11.
* **Errors 12 / Warning 3 warping measure** (proposal): w = |n̂·(x4 − x1)| / sqrt(A), with n̂ from the
  diagonals' cross product.
  * Warning 3 when w > 1e-3.
  * Error 12 when w > 5e-2.
  * Make both thresholds configurable.

### 5.3 Warning catalogue (warnings never block AFWRITE)

| # | Message | Condition | AFWRITE / analysis consequence |
|---|---|---|---|
| W1 | Gap Found at Node `<n>` | node number n is missing (n < max node, undefined) | AFWRITE generates node n with **fixed DOFs**, in the analysis files only, not in the model |
| W2 | Distorted Element `<e>`, Group `<g>`, Face `<f>` | SHELL/PLANE (3 or 4 nodes): smallest interior angle too small compared with the largest. SOLIDs: every face is checked. Threshold not given (OQ-30); proposal: warn if min angle < 15° or max angle > 165°. | none; results may be poor |
| W3 | Warped Element `<e>`, Group `<g>` | 4-node SHELL/PLANE nodes "not quite" coplanar (between the W3 and Error 12 thresholds) | none |
| W4 | Unused Node `<n>` | node not used by any element. The manual excepts K-nodes of BEAMS; interpret as "K-nodes are counted as used" (OQ-31). | AFWRITE fixes its DOFs in the analysis files |
| W5 | Translational Mass in Node `<n>` Is on Fixed DOF | translational mass on a fixed DOF | mass ignored by the analysis |
| W6 | Rotational Mass in Node `<n>` Is on Fixed DOF | rotational mass on a fixed DOF | mass ignored |
| W7 | Group `<g>` Has no Elements | empty group | AFWRITE skips the group |
| W8 | Too Many Top Layers | SITE top-layer list > 100 | only the first 100 are written (conflicts with the 200 limit stated elsewhere; 05a) |
| W9 | Number of Values for Fourier Transform Is Not Power of 2 | NFFT not 2^m | the nearest power of 2 is written (05a: 3000 → 2048) |
| W10 | Force in Node `<n>` Is on Fixed DOF | force on a fixed DOF | ignored |
| W11 | Moment in Node `<n>` Is on Fixed DOF | moment on a fixed DOF | ignored |

### 5.4 Suggested additional checks (not in the manual; label as "ACS-edu" checks so numbering stays faithful)

* Poisson ratio ν ≥ 0.5 (an incompressibility singularity for SOLIDs).
* Mesh passing frequency: element size h > Vs/(5·f_max).
* Interaction-node bottom-up numbering rules (05b).
* Option NON panels with ≠ 4 detected corners.
* BBC monotonicity and consistency (§2.2).
* Zero-length SPRING with all constants 0.

---

## 6. References (Manual Chapter 11)

Role flags:
* **KEY-SSI**: core theory for SASSI's flexible-volume/substructuring, site and impedance methodology.
* **KEY-SOIL**: free-field and equivalent-linear soil theory.
* **KEY-INCOH**: motion incoherency (Option [ADV]).
* **KEY-NON**: nonlinear structure and wall capacity (Option NON, SHEAR/BBCGEN).
* **STD**: codes and standards.
* **BKG**: background (Ghiocel and others: applications, probabilistic/Option PRO, validation studies).

| # | Citation | Topic | Role |
|---|---|---|---|
| 1 | Abrahamson, N. (1993). Spatial variation of multiple support inputs. 1st US Seminar on Seismic Evaluation & Retrofit of Steel Bridges, UC Berkeley | coherency (1993 model, Type 2) | KEY-INCOH |
| 2 | Abrahamson, N. (2005). Spatial Coherency for Soil-Structure Interaction. EPRI Report 1012968 | 2005 coherency model (Type 3) | KEY-INCOH |
| 3 | Abrahamson, N. (2006). Program on Technology Innovation: Spatial Coherency for SSI. EPRI/DOE Report 1014101 | 2006 model, embedded foundations (Type 4) | KEY-INCOH |
| 4 | Abrahamson, N. (2007). Hard Rock Coherency Functions Based on the Pinyon Flat Data. EPRI/DOE Report 1015110 | 2007 hard-rock model (Type 5) | KEY-INCOH |
| 5 | ACI 349-06 (2006). Code Requirements for Nuclear Safety-Related Concrete Structures | concrete design | STD |
| 6 | ASCE 43-05 (2005). Seismic Design Criteria for SSCs in Nuclear Facilities | inelastic absorption factors, limit states | STD / KEY-NON |
| 7 | ASCE 4-17 (2017). Seismic Analysis of Safety-related Nuclear Structures and Commentary | cracking stress 3√f'c (§C.3.3.2), damping/stiffness response levels | STD / KEY-NON |
| 8 | ATC PEER/ATC 72-1 (2010). Modeling and Acceptance Criteria for Seismic Design and Analysis of Tall Buildings | cyclic degradation of shear capacities | KEY-NON |
| 9 | Barda, F., Hanson, J.M., Corley, G.W. (1977). Shear Strength of Low-Rise Walls with Boundary Elements. ACI SP 53-8 | Barda equation (ShearModel 3) | KEY-NON |
| 10 | Gergely, P. (1984). Seismic Fragility of RC Structures and Components for Nuclear Facilities. NUREG/CR-4123 | shear-dominated low-rise walls | KEY-NON |
| 11 | Ghiocel, D.M. (2015a). SASSI Flexible Volume Substructuring Methods for Deeply Embedded Structures; Selection of Excavated Soil Interaction Nodes and Element Meshing. SMiRT22 [sic], Manchester | FV, FI, MSM interaction-node selection (INTGEN) | KEY-SSI |
| 12 | Ghiocel, D.M. (2015b). Seismic Motion Incoherency Effects of SSI and SSSI of Nuclear Structures for Different Soil Site Conditions. SMiRT22 [sic], Manchester | incoherent SSI/SSSI | KEY-INCOH |
| 13 | Ghiocel, D.M. (2015c). Fast Nonlinear Seismic SSI Analysis Using a Hybrid Time-Complex Frequency Approach for Low-Rise Nuclear Concrete Shearwall Buildings. SMiRT22 [sic], Manchester | **Option NON methodology** | KEY-NON |
| 14 | Ghiocel, D.M. (2014a). Effects of Seismic Motion Incoherency on SSI and SSSI Responses… DOE NPH Meeting, Germantown | incoherency | KEY-INCOH |
| 15 | Ghiocel, D.M. (2014b). SASSI Methodology-Based Sensitivity Studies for Deeply Embedded Structures such as SMRs. DOE NPH Meeting | deeply embedded SASSI | KEY-SSI |
| 16 | Ghiocel, D.M., Yue, D., Fuyama, H., Kitani, T., McKenna, M. (2013a). Validation of Modified Subtraction Method for Seismic SSI Analysis of Large-Size Embedded Nuclear Islands. SMiRT22, San Francisco | **MSM (FI-EVBN) validation** | KEY-SSI |
| 17 | Ghiocel, D.M. (2013b). Comparative Studies on Seismic Incoherent SSI Analysis Methodologies. SMiRT22, San Francisco | incoherency methods comparison | KEY-INCOH |
| 18 | Ghiocel, D.M., Todorovski, L., Fuyama, H., Mitsuzawa, D. (2011a). Seismic Incoherent SSI Analysis of a Reactor Building Complex on a Rock Site. SMiRT21, New Delhi, Paper 825 | incoherency application | BKG |
| 19 | Ghiocel, D.M., Lee, I. (2011b). Seismic SSI Effects for Large-Size Surface and Embedded Nuclear Facility Structures. ASEM11, Seoul | SSI application | BKG |
| 20 | Ghiocel, D.M. (2010a). Some Insights and Brief Guidance for Application of Subtraction/Flexible Interface Method to Seismic SSI Analysis of Embedded Nuclear Facilities. GPT-TIR-01-0930-2010 | **subtraction / flexible-interface guidance** | KEY-SSI |
| 21 | Ghiocel, D.M., Todorovski, L., Fuyama, H. (2010b). Seismic SSI Response of Reactor Building Structures. OECD NEA/IAEA SSI Workshop, Ottawa | application | BKG |
| 22 | Ghiocel, D.M., Short, S., Hardy, G. (2010c). Seismic Motion Incoherency Effects for Nuclear Complex Structures on Different Soil Site Conditions. OECD NEA/IAEA Workshop | incoherency | KEY-INCOH |
| 23 | Ghiocel, D.M., Li, D., Brown, N., Zhang, J.J. (2010d). EPRI AP1000 NI Model Studies on Seismic SSSI Effects. OECD NEA Workshop | SSSI | BKG |
| 24 | Ghiocel, D.M., Stoyanov, G., Adhikari, S., Aziz, T. (2010e). Seismic Motion Incoherency Effects for CANDU Reactor Building Structure. OECD NEA/IAEA Workshop | incoherency application | BKG |
| 25 | Ghiocel, D.M., Li, D., Tunon-Sanjur, L. (2009a). Seismic Motion Incoherency Effects for AP1000 Nuclear Island Complex. SMiRT20, Helsinki, Paper 1852 | incoherency application | BKG |
| 26 | Ghiocel, D.M., Short, S., Hardy, G. (2009b). Seismic Motion Incoherency Effects on SSI Response of Nuclear Islands with Significant Mass Eccentricities and Different Embedment Levels. SMiRT20, Paper 1853 | incoherency with eccentricity | KEY-INCOH |
| 27 | Ghiocel, D.M. (2007a). Stochastic and Deterministic Approaches for Incoherent Seismic SSI Analysis as Implemented in ACS SASSI. Appendix C, EPRI/DOE TR-1015110 | **ACS SASSI incoherency algorithms** (stochastic simulation, SRSS, AS, …) | KEY-INCOH |
| 28 | Ghiocel, D.M., Ostadan, F. (2007b). Seismic Ground Motion Incoherency Effects on SSI Response. SMiRT19, Toronto, Paper K05/4 | incoherency | KEY-INCOH |
| 29 | Ghiocel, D.M. (2004). Stochastic Simulation in Engineering Predictions. Ch. 20, CRC Engineering Design Reliability Handbook | stochastic simulation (Option PRO) | BKG |
| 30 | Ghiocel, D.M., Wang, L. (2004). Seismic Motion Incoherency Effects on Structures. 3rd UJNR SSI Workshop, Menlo Park | incoherency | BKG |
| 31 | Ghiocel, D.M., Ghanem, R. (2002). Stochastic Finite-Element Analysis of Seismic SSI. ASCE J. Eng. Mech. 128(1) | stochastic FE / probabilistic SSI (Option PRO) | BKG |
| 32 | Ghiocel, D.M. (1998). Uncertainties of Seismic SSI Analysis: Significance, Modeling and Examples. US-Japan SSI Workshop, USGS Menlo Park | uncertainty | BKG |
| 33 | Ghiocel, D.M. et al. (1996a). Seismic Motion Incoherency Effects on Dynamic Response. 7th ASCE Specialty Conf. Probabilistic Mechanics & Structural Reliability, Worcester | incoherency (looks like a duplicate of #35) | BKG |
| 34 | Ghiocel, D.M. et al. (1996b). On SSI Issues for Deep Foundation Structures. 6th Symp. NPP Structures, Equipment & Piping, Raleigh | deep foundations | BKG |
| 35 | Ghiocel, D.M. et al. (1996c). Seismic Motion Incoherency Effects on Dynamic Response. 7th ASCE EMD/STD Joint Specialty Conf., Worcester | incoherency | BKG |
| 36 | Ghiocel, D.M. et al. (1996d). Probabilistic Seismic Analysis Including SSI. same conference | probabilistic SSI | BKG |
| 37 | Ghiocel, D.M. et al. (1996e). Effects of Random Field Modeling of Seismic Ground Motion on Structural Dynamic Response. 37th AIAA/ASME/ASCE SDM Conf., Salt Lake City | random fields | BKG |
| 38 | Ghiocel, D.M. et al. (1995). Seismic SSI Effects on Probabilistic Floor Response Spectra. ASME PVP, Honolulu | probabilistic FRS | BKG |
| 39 | Ghiocel, D.M. et al. (1991). Evaluation of Seismic SSI by Different Approaches. 2nd Int. Conf. Recent Advances in Geotech. EQ Eng., St. Louis | SSI methods comparison | BKG |
| 40 | Ghiocel, D.M. et al. (1990a). Seismic SSI Effects for Shearwall Buildings on Soft Clays. 10th ECEE, Moscow | application | BKG |
| 41 | Ghiocel, D.M. et al. (1990b). Evaluation of SSI and SSSI Effects on Seismic Response of Nuclear Heavy Buildings by Different Approaches. 10th ECEE | SSSI | BKG |
| 42 | Ghiocel, D.M. (1986a). Probabilistic Seismic SSI Analysis. 8th ECEE, Lisbon | probabilistic SSI | BKG |
| 43 | Ghiocel, D.M. et al. (1986b). Effects of Spatial Character of Seismic Random Excitations on Structures. 8th ECEE | spatial variation | BKG |
| 44 | Ghiocel, D.M. (1986c). PRELAMOS: Computer Program for Stochastic Parameter Estimation of Accelerograms. 5th Natl. Symp. Informatics in Civil Eng., Sibiu | accelerogram stochastic parameters | BKG |
| 45 | Ghiocel, D.M. (1985a). Probabilistic Seismic Analysis of the Containment Structure of a Nuclear Reactor Building. J. Civil Eng. 11, Bucharest | probabilistic | BKG |
| 46 | Ghiocel, D.M. et al. (1985b). Seismic Risk Evaluation for Buildings Including SSI. Sci. Bull. Civil Eng., ICB Bucharest | risk | BKG |
| 47 | Ghiocel, D.M. et al. (1983a). Actual Tendencies in Seismic Analysis of NPP Structures. Sci. Bull. Civil Eng. vol. 2, ICB | review | BKG |
| 48 | Ghiocel, D.M. et al. (1983b). Advanced Computational Methods in Seismic Analysis of NPPs. Sci. Bull. Civil Eng. vol. 1, ICB | methods | BKG |
| 49 | Ghiocel, D.M. et al. (1983c). Structural Reliability Analysis of NPP Subjected to Seismic Load Using Advanced Numerical Methods. Natl. Symp. Cybernetic Applications in Industry, Bucharest | reliability | BKG |
| 50 | Gulec, C.K., Whittaker, A.S. (2009). Performance-Based Assessment and Design of Squat Reinforced Concrete Shear Walls. MCEER-09-0010 | **G–W equation 6-9; source compilation for Barda/Wood/ACI** (ShearModels 1–4) | KEY-NON |
| 51 | Hardin, B.O., Drnevich, V.P. (1972). Shear Modulus and Damping in Soils: Design Equations and Curves. JSMFD ASCE 98(SM7), 667–692 | G/Gmax and damping curves | KEY-SOIL |
| 52 | Johnson, J.J., Short, S.A., Hardy, G.S. (2007). Modeling Seismic Incoherence Effects on NPP Structures: Unifying CLASSI and SASSI Approaches. SMiRT19, K05/5 | incoherency (CLASSI vs SASSI) | KEY-INCOH |
| 53 | Reed, J.W., Kennedy, R.P. (1994). Methodology for Developing Seismic Fragilities. EPRI TR-103959 | fragility; shear-dominated walls | KEY-NON / BKG |
| 54 | Luco, J.E., Mita, A. (1987). Response of Circular Foundation to Spatially Random Ground Motion. J. Eng. Mech. ASCE 113(1), 1–15 | coherency, rigid foundations | KEY-INCOH |
| 55 | Luco, J.E., Wong, H.L. (1986). Response of a Rigid Foundation Subjected to a Spatially Random Ground Motion. EESD 14, 891–908 | **Luco–Wong coherency model** (Type 1; Error 58) | KEY-INCOH |
| 56 | Lysmer, J., Kuhlemeyer, R.L. (1969). Finite Dynamic Model for Infinite Media. J. Eng. Mech. Div. ASCE 95(EM4), 859–877 | **viscous (transmitting) boundaries**, used for the half-space base dashpots | KEY-SSI |
| 57 | Lysmer, J., Udaka, T., Tsai, C.-F., Seed, H.B. (1975). FLUSH: A Computer Program for Approximate 3-D Analysis of SSI Problems. EERC 75-30 | equivalent-linear FE SSI, transmitting boundaries | KEY-SSI |
| 58 | Lysmer, J. (1978). Analytical Procedures in Soil Dynamics. EERC 78/29 | **thin-layer / layered-media eigenproblem; procedures behind SITE/POINT** | KEY-SSI |
| 59 | Lysmer, J., Tabatabaie-Raissi, M., Tajirian, F., Vahdani, S., Ostadan, F. (1981). SASSI: A System for Analysis of Soil-Structure Interaction. Report UCB/GT 81-02 | **THE primary SASSI theory and user reference** (modules SITE, POINT, HOUSE, ANALYS, MOTION, STRESS) | KEY-SSI |
| 60 | Nie, J., Braverman, J., Costantino, M. (2013). Seismic SSI Analyses of a Deeply Embedded Model Reactor: SASSI Analyses. BNL 102434-2013 | FV vs subtraction verification for deep embedment | KEY-SSI / BKG |
| 61 | Ostadan, F., Tseng, W.S., Lilhanand, K. (1987). Application of Flexible Volume Method to SSI Analysis of Flexible and Embedded Foundation. 9th SMiRT, Lausanne | flexible volume method application | KEY-SSI |
| 62 | Seed, H.B., Idriss, I.M. (1969). The Influence of Soil Conditions on Ground Motion during Earthquakes. JSMFD ASCE (SM1), 99–137 | **site response / equivalent-linear basis** | KEY-SOIL |
| 63 | Seed, H.B., Idriss, I.M. (1970). Soil Moduli and Damping Factors for Dynamic Response Analysis. EERC 70-10 | **standard G/Gmax and damping curves** (Seed–Idriss iterative equivalent-linear) | KEY-SOIL |
| 64 | Short, S.A., Hardy, G.S., Merz, K.L., Johnson, J.J. (2006). Effect of Seismic Wave Incoherence on Foundation and Building Response. EPRI TR-1013504 | **EPRI incoherence study** | KEY-INCOH |
| 65 | Short, S.A., Hardy, G.S., Merz, K.L., Johnson, J.J. (2007). Validation of CLASSI and SASSI to Treat Seismic Wave Incoherence in SSI Analysis of NPP Structures. EPRI TR-1015110 | **EPRI validation of SASSI incoherency** | KEY-INCOH |
| 66 | Tabatabaie-Raissi, M. (1982). Flexible Volume Method for Dynamic SSI Analysis. Ph.D. Diss., UC Berkeley | **flexible volume method (FV) theory** | KEY-SSI |
| 67 | Tajirian, F. (1981). Impedance Matrices and Interpolation Techniques for 3-D Interaction Analysis by the Flexible Volume Method. Ph.D. Diss., UC Berkeley | **impedance matrices (POINT/ANALYS) and transfer-function interpolation (MOTION)** | KEY-SSI |
| 68 | Tseng, W.S., Lilhanand, K. (1997). SSI Analysis Incorporating Spatial Incoherence of Ground Motions. EPRI TR-102631 2225 | incoherent SSI | KEY-INCOH |
| 69 | US NRC ISG-01 (2008). Seismic Issues Associated with High Frequency Ground Motion in DC and COL Applications | regulatory basis for incoherency and high-frequency motion | STD |
| 70 | Vahdani, S. (1984). Impedance Matrices for SSI Analysis by the Flexible Volume Method. Ph.D. Diss., UC Berkeley | **impedance matrix computation** | KEY-SSI |
| 71 | Wass [sic: Waas], G. (1972). Earth Vibration Effects and Abatement for Military Facilities: Analysis Method for Footing Vibrations through Layered Media. Tech. Rep. S-71-14, USAE WES, Vicksburg | **origin of the thin-layer / transmitting-boundary formulation for layered media (SITE)** | KEY-SSI |
| 72 | Wood, S. (1990). Shear Strength of Low-Rise Reinforced Concrete Walls. ACI Struct. J. 87, 99–107 | Wood equation (ShearModel 2) | KEY-NON |

### 6.1 Key theory references by module (for the re-implementation)

| Module / feature | Primary references in the list |
|---|---|
| SITE (layered free field, thin-layer eigenproblem, half-space simulation by sublayers + viscous base) | 58 Lysmer 1978; 71 Waas 1972; 56 Lysmer–Kuhlemeyer 1969; 59 SASSI 1981 |
| POINT (point-load Green's functions, central zone) | 59; 58; 67 Tajirian 1981; 70 Vahdani 1984 |
| HOUSE/ANALYS (flexible volume and subtraction substructuring, impedance) | 59; 66 Tabatabaie-Raissi 1982; 67; 70; 61 Ostadan et al. 1987; 11, 15, 16, 20 Ghiocel (FV/FI/MSM); 60 Nie et al. 2013 |
| MOTION (TF interpolation, time histories) | 67 Tajirian 1981; 59 |
| SOIL (equivalent-linear SHAKE methodology) | 62, 63 Seed–Idriss; 51 Hardin–Drnevich; 57 FLUSH |
| Incoherency (HOUSE/ANALYS [ADV]) | 55 Luco–Wong; 54 Luco–Mita; 1–4 Abrahamson; 64, 65 Short et al. EPRI; 52 Johnson et al.; 68 Tseng–Lilhanand; 27 Ghiocel 2007a (ACS SASSI algorithms); 17, 28; 69 NRC ISG-01 |
| Option NON / SHEAR / BBCGEN | 13 Ghiocel 2015c; 50 Gulec–Whittaker; 72 Wood; 9 Barda; 6 ASCE 43-05; 7 ASCE 4-17; 8 ATC 72-1; 10 Gergely; 53 Reed–Kennedy |
| Option PRO (probabilistic) | 29, 31, 32, 36, 38, 42 Ghiocel |

### 6.2 Important references *not* in the manual's list (the implementer will need them; external, flagged)

The manual relies on the following methods but does not cite them. The citations below are standard literature
and do not come from the manual.

* **SHAKE methodology** (SOIL module): Schnabel, P.B., Lysmer, J., Seed, H.B. (1972). *SHAKE: A computer
  program for earthquake response analysis of horizontally layered sites*. EERC 72-12, UC Berkeley. Also
  SHAKE91 (Idriss & Sun, 1992), which the manual names as a verification code. **Schnabel is not cited in the
  manual.**
* **Kausel** layered-media and thin-layer Green's functions, relevant to SITE/POINT: e.g. Kausel & Roesset
  (1981), *Stiffness matrices for layered soils*, BSSA 71(6); Kausel & Peek (1982), *Dynamic loads in the
  interior of a layered stratum: an explicit solution*, BSSA 72(5). **Kausel is not cited in the manual.**
* **DEEPSOIL / MKZ hyperbolic model** (NLSLAYER/NLSOIL): Hashash & Park (2001), Eng. Geology 62;
  Matasovic & Vucetic (1993), J. Geotech. Eng. 119(11).
* **Hysteretic models** (Option NON): Takeda, Sozen & Nielsen (1970), J. Struct. Div. ASCE 96(ST12), for TAK;
  Cheng & Mertz (1989) (University of Missouri-Rolla) for CMS/CMB; Masing (1926) for GMR.
* **ACI 318-08** (SHEAR/BBCGEN ShearModel 1, α_c), plus the related cracking-stress basis.
* For FILLPOOL validation: Westergaard (1933) hydrodynamic pressure and Housner (1963) impulsive/convective
  pool models.

### 6.3 Citation inconsistencies spotted (keep the manual's text in UI help; correct in the docs)

* Refs 11–13 (2015, Manchester) are labelled "SMiRT22". SMiRT22 was San Francisco 2013 (refs 16–17); Manchester
  2015 was SMiRT23.
* Refs 33 and 35 have the same title and conference, so they are probably duplicates.
* Report number 1015110 is given for three different documents: refs 4, 27 and 65. Spec 05b cites
  TR-1015111 for the 2007 validation, so verify which number is correct.
* Ref 62 gives "Vol. 94 … December" for a 1969 SM1 paper. Verify the volume and month.
* The text cites "ACI 349-08", "ASCE 43-17" and "ACI 318-08", but the list has ACI 349-06 and ASCE 43-05, and
  ACI 318-08 is missing.
* Ref 71 author is spelled "Wass"; the correct spelling is Waas.

---

## 7. Verification tests (for a verifiable re-implementation)

1. **REFINEMODEL**
   * A single quad becomes 4 quads with 5 new nodes (4 midpoints + 1 centre). A 2×2 quad mesh becomes 16 quads
     with **shared** midpoints, giving 25 nodes in total.
   * A single hexahedron becomes 8 with 19 new nodes.
   * The total area or volume is unchanged to round-off.
   * Triangles and prisms are untouched.
2. **FILLPOOL**
   * Rectangular shell pool, plan a×b, wall height H, levels every h, `EmptyLevels` = m: the water solids'
     volume equals a·b·(H − m·h).
   * The number of springs equals the number of coincident wall-water node pairs.
   * Spring constants are `Stiff` along the normal and `stiff2` tangentially.
   * `offset` 0 < offset < Nmax raises an error.
3. **BBCGEN**
   * The BBC has 22 points, yield index 21, point 22 = (0.02, 1.02·V_u), and points 2–21 are equally spaced in
     strain.
   * CFL = 0.3 gives Y1 = 0.3·V_u.
   * CFL = 0.6 is rejected.
   * The Wood clamp triggers for low or high ρ_V.
   * Hand calculations of the four SHEAR equations agree within 1e-6 relative.
4. **NONLINMOTDISP**
   * A 4×3 quad panel gives exactly its 4 corner nodes.
   * The node lists show no duplicates afterwards (Error 66 not triggered).
5. **EDGE**
   * Reproduce Fig. 1.5: `EDGE,1,0,0,1` gives 3 groups; then `EDGE,2` gives groups 2, 4 and 5.
   * The element count is conserved.
6. **Binary DB**
   * Write DBs for X, Y and Z, then COMB*: each value equals the algebraic sum of the three inputs.
   * BINFRAMEOUT at index k equals the slice at k.
   * MAXDBFRAME equals the signed abs-max.
   * BINSTRTBL with step −1 equals the signed abs-max per component, and gives 0-rows for missing elements.
7. **THSHLSTR**
   * Pure membrane Nxx with t = 1: σ_top = σ_bot = Nxx.
   * Pure bending Mxx: σ = ±6M/t².
8. **CHECK**: one unit test per error and warning number, building a minimal model or option set that violates
   exactly that rule. Also test:
   * fatal-stop behaviour for 42 and 43;
   * module gating in AFWRITE (a module with an error is not written; warnings do not block);
   * no Error 8 for zero-length springs created by SOLIDPILE.
9. **NLSLAYER** curve fit: fit to synthetic G/Gmax points generated from known (β, s, γr) and recover them
   within tolerance.

---

## 8. Open questions / ambiguities (implementer decisions)

* **OQ-1** FILLPOOL and EXCAV: the rule for merging nearby Z values (mean, min or first) and the tolerance
  semantics of `Sensitivity`. Proposed: single-linkage clustering with gap ≤ Sensitivity, using the mean value.
* **OQ-2** FILLPOOL springs on **oblique** walls. SC constants are global and uncoupled, so a true normal spring
  cannot be represented. Options:
  * (a) restrict to axis-aligned walls and error out otherwise;
  * (b) project: kx = Stiff·nx² + stiff2·(1 − nx²), etc. This is an approximation;
  * (c) use GENERAL (matrix) elements with a full 3×3 coupling.

  Corner nodes shared by two walls, or by a wall and the floor, need both normals. Proposed: sum the
  normal-direction stiffnesses. Also undecided: whether floor nodes get springs (proposed: yes, with Z normal).
* **OQ-3** FILLPOOL spring damping value. Not given; SOLIDPILE uses 4 %. Proposed: 0, or user-set later with SC.
* **OQ-4** Water SOLID material: not specified. Proposal:
  * unit weight 62.4 pcf (9.81 kN/m³);
  * ν = 0.49;
  * E chosen so that the bulk modulus K = E/(3(1 − 2ν)) ≈ 2.2 GPa (≈ 4.6e4 ksf), giving E ≈ 1.3e5 kPa at
    ν = 0.49;
  * small damping, about 0.5 %.

  The resulting small but non-zero G (shear) is an approximation of a fluid. This should be documented as an
  educational simplification. Sloshing (convective mode) is not captured by these solids.
* **OQ-5** Shell-area group: element thickness, material, mass, whether it is built on wall-side or water-side
  nodes, and whether it should be excluded from the dynamic analysis. Proposed: dummy shells (ETYPE = 1, very
  low E, zero mass), recorded for ANSYS export only.
* **OQ-6** REFINEMODEL scope:
  * PLANE quads, TSHELL quads (proposed: refine TSHELL like SHELL; refine PLANE in 2D models);
  * handling of interaction-node status of new nodes;
  * loads, masses, output requests and node restraints D on new edge nodes (proposed: inherit restraint if both
    end nodes share it).
* **OQ-7** BBCP outside the existing point range: grow the curve or error out?
* **OQ-8** BBCGEN arguments 7–8 (`bre`, `bys`) versus the description's `Fvw`, `Fbe`. SHEAR also prints
  `[Fvw],[Fbe]` in its syntax but describes `<Abe>`, `<Fybe>`. Proposed: A_BE and f_y,BE, giving
  F_BE = A_BE·f_y,BE, with F_VW computed internally as ρ_V·A_w·f_y.
* **OQ-9** Unit consistency of the √f'c equations (psi-based coefficients vs ksi or kN/m² input). Also α_c
  (ACI) is undefined, and the shear area used for cracking strain (A_w vs 5/6·A_w) is not stated.
* **OQ-10** BBCGEN:
  * yield strain γ_y;
  * curve shape between cracking and yield;
  * whether the yield force = V_u;
  * which BBC number is written;
  * the BBC type assigned (CMS = 1 proposed);
  * the order of operations for `Panel = 0` (every panel needs uniform properties).
* **OQ-11** Range deletions with non-existent IDs: skip silently or warn? Renumbering after deletion: proposed
  none, IDs are user keys.
* **OQ-12** WALLFLR:
  * coplanarity tolerances (angle and offset);
  * whether coplanar but **disconnected** shell sets (two separate walls in the same plane) go into one group;
  * the exact GTIT title strings;
  * how "floor" (normal ∥ Z) and "oblique" are classified.
* **OQ-13** EDGE:
  * angle tolerance for "parallel to an axis";
  * treatment of oblique edges and of walls oblique to the global axes. Proposed: classify by in-plane local
    horizontal/vertical, and map the X/Y flags to the "horizontal" direction for oblique walls.
* **OQ-14** MERGEPANEL node policy: reuse nodes by number, or by coordinates; behaviour for mismatches.
* **OQ-15** The order of PNLGEN relative to MERGEPANEL. Group numbers change on merge, so panels must reference
  the final group numbers. Proposed: run PNLGEN after MERGEPANEL.
* **OQ-16** SOLIDPILE:
  * how bottom-perimeter nodes are classified;
  * per-pile versus whole-group "bottom" detection;
  * spring-property (SC) creation and numbering;
  * whether duplicate nodes copy restraints;
  * how the four new groups are named.
* **OQ-17** NONLINBAT: batch file name, shell (DOS `.bat` vs portable script), iteration count and convergence
  test (see 05d OQ-N5).
* **OQ-18** NLSLAYER:
  * units of `refStrain` (% or decimal);
  * the meaning of `Num` (free-field SOIL sublayer vs `L` layer number);
  * the curve-fit objective and point weighting;
  * hysteresis rule (Masing / non-Masing / DEEPSOIL MRDF);
  * how `Vis` interacts with NLDampType.
* **OQ-19** NLSOIL time-domain solver details: Newmark parameters, sub-incrementation, convergence norms (DispConv
  and ForceConv in absolute or relative units), behaviour when EqualIt is exceeded, and the input/output file
  layout.
* **OQ-20** Binary database layouts (ACC, DISP, THS, and the ANSYS-format stress DB) and the "compression"
  scheme are undocumented. A schema must be defined (§3.6).
* **OQ-21** COMB* combination rule. Proposed: algebraic time-wise sum. Also to decide: phasing, time-step
  alignment, and handling of databases with different node or element sets.
* **OQ-22** BINFRAMEOUT time-to-index mapping (round vs floor, 0- vs 1-based frames) and the exact ASCII frame
  header (rows/cols header per FRAMECOMBIN/MODFRAMES).
* **OQ-23** THSHLSTR face formulas and how the "maximum basic components" combine (sign permutations vs
  simultaneous values). Also the face-strain definitions.
* **OQ-24** THSHLSMH:
  * the definition of "adjacent element layers" (shared node vs shared edge);
  * the Parzen window width;
  * whether to smooth across groups or planes;
  * the output frame format and name.
* **OQ-25** Numeric tolerances for all geometric CHECK tests (coincident, collinear, zero area, warping,
  distortion).
* **OQ-26** Errors 52 and 63 state the open interval (0, 360). Literal reading makes 0° illegal, which conflicts
  with vertical incidence and the zero-angle default. Proposed: [0, 360).
* **OQ-27** Error 60 (mesh points/embedment level ≤ 0) conflicts with a screenshot value of 0 (05b). Proposed:
  apply only when incoherency is active with an embedded-foundation coherency model.
* **OQ-28** Error 124: does "fully fixed" mean all 6 DOFs, or all 3 translations? The interaction DOFs are
  translations. Proposed: all three translational DOFs fixed.
* **OQ-29** Error 8 scope: it must exclude zero-length SPRING (and possibly GENERAL) elements. The manual says
  "two-node elements" without exception.
* **OQ-30** Warning 2 distortion criterion: "smallest angle too small compared with the greatest angle". Proposed
  thresholds: min angle < 15°, max angle > 165°, or ratio min/max < 0.2. Make them configurable.
* **OQ-31** Warning 4 "except K-nodes from BEAMS": either K-nodes are excluded from the warning, or nodes used
  only as K-nodes are still reported. The bubble-plot note in Ch. 7 counts K-nodes among the unused nodes.
  AFWRITE fixing a K-node's DOFs is harmless because K-nodes only define orientation. Proposed: do not warn for
  K-nodes, but still fix them in the analysis files.
* **OQ-32** Errors 116/117/104/109/102: the manual says only "illegal". The exact predicates (range limits,
  interaction-node membership, format grammar) are inferred.
* **OQ-33** Module attribution of Errors 1, 20–25, 44, 49, 50, 64–67 and 73–78. The variables are shared across
  tabs (SITE/SOIL/MOTION/STRESS). Proposed: report the error under **every** module that uses the variable, so
  that each such module is gated, because AFWRITE skips modules with errors.
* **OQ-34** Warning 8 (100 top layers) conflicts with the documented 200-layer limit (05a).
* **OQ-35** DELNLS/DELSPR syntax and text are swapped or mis-printed in the manual (§2.10). The intended commands
  are inferred.
