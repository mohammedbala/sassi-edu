# 08 — Command Reference: Node, Element and Load Commands (Ch. 9.3 – 9.5)

**Source:** ACS SASSI V3 User Manual, Chapter 9, §9.3 *Node Commands*, §9.4 *Element Commands*,
§9.5 *Load Commands* (extracted-text lines 9717–10687; printed pages 237–259; PDF pages 239–261).
I viewed all of PDF pages 239–261 as images. Figures 9.1–9.4, the GENERAL-element DOF/matrix
figure, the beam section-property figure and all tables are transcribed below. The page images
were also used to check which part of each command name the manual underlines (that is the
accepted abbreviation, see §1.3).

Context taken from other parts of the manual is marked *(xref §x / line n)*. Those facts are
specified in full in the other spec files: 01_capabilities, 02_theory_modules, 05a (MOPT dialog),
05b (HOUSE/FORCE options) and the Chapter-10 CHECK messages.

**Tags:** **[CORE]** needed for computational correctness · **[IO]** file/input format or
command syntax · **[UI]** interface, listing or convenience · **[ADV]** advanced or legacy
(incoherency, Option A/AA/NON, deprecated fields).

**Provenance markers**
- *(manual)*: stated in the manual.
- *(figure)*: read from a figure or table image.
- *(inferred)*: my reconstruction to fill a gap, or a design decision for this implementation.
  Every inferred item is repeated in §12 *Open questions*.

---

## 0. Quick index

| § | Command group | Commands |
|---|---|---|
| 3 | Node (9.3) | CSYS, D, FILL, GLOBAL, INT, INTLIST, LMOVE, LOC, LOCAL, N, NDEL, NGEN, NLIST, NMED, NMOVE, NSCALE, SDEL, SLIST |
| 4 | Element library | SOLID, BEAMS, SHELL, TSHELL, PLANE, SPRING, GENERAL |
| 5 | Element (9.4) | DELL, DELM, DELR, DELSC, E, ECOMPR, EDEL, EGEN, EINT, ELIST, ETYPE, GDEL, GLIST, GROUP, GTIT, KI, KJ, L, LLIST, M, MACT, MLIST, MSET, MTYPE, MXDEL, MXI, MXLIST, MXM, MXR, R, RACT, RLIST, RSET, SC, SCLIST, THICK |
| 6 | Property tables | material conversion, soil-layer use, beam sections, springs, matrix properties |
| 7 | Load (9.5) | F, FDEL, FLIST, FSCALE, MM, MMDEL, MMLIST, MR, MRGEN, MRDEL, MRSCALE, MSCALE, MT, MTDEL, MTGEN, MTLIST, MTSCALE, MUNITS |
| 8 | CHECK rules linked to these commands | Errors 6–43, 61, 62, 83, 124; Warnings 1–7, 10, 11 |
| 9 | ANSYS cross-walk | for users coming from ANSYS |
| 10 | Worked input example | |
| 11 | Verification tests | |
| 12 | Open questions | |

---

## 1. Conventions shared by all commands

### 1.1 Syntax notation [IO] *(manual, §9 intro, line 8693)*
- The manual writes commands in Backus-Naur style.
  - `<value>` is replaced by a number or a string.
  - Items in `[ ]` are optional.
  - Every other character (command name, commas) is typed literally.
- Arguments are **comma-separated** and **positional**.
- **Command names are not case-sensitive.**
- *(inferred)* Parsing rules, to be shared with the general parser spec:
  - An empty field (`,,`) or an omitted trailing field takes the documented default. A field with no stated default takes 0 (numeric) or "" (string).
  - Numbers are read in free format. Accept Fortran-style exponents: `1.0E5`, `1.0e+5`, `1.0D5`.
  - Integer arguments given as reals (e.g. `5.0`) are rounded to the nearest integer. Warn if the fractional part is non-zero.
  - Surrounding blanks are trimmed.
  - String arguments (GROUP type, D labels, GTIT title) are case-insensitive keywords. GTIT titles are free text, case preserved.
  - Extra trailing arguments beyond those documented are ignored with a warning, except where §3.10 (N) says they are legacy fields.

### 1.2 Command execution model [IO][UI]
- Commands act on the **active model**.
- They can be typed at the UI command line or read from a `.pre` batch file (INP, or Model > Input).
- *(xref line 4130)* An invalid command name prints `<X> Command not found`. Processing of the `.pre` file then **continues**.
- *(inferred)* Every command in this file has the same effect whether typed or read from `.pre`.
- *(inferred)* WRITE (xref §9.2.47) must be able to regenerate an equivalent command stream for all data defined here.

### 1.3 Accepted abbreviations [IO] *(manual line 8697 + figure underlines)*
- Old PREP looked at only the **first 4 letters** of a command name. The new UI reads the full name.
- In the manual, a command longer than 4 characters has its **accepted abbreviation underlined**. A user may type either the full name or exactly the underlined part. Any other variant is rejected.
- Underlines seen on the page images for this section:

| Full name | Underlined (accepted short form) | Full name | Underlined |
|---|---|---|---|
| GLOBAL | GLOB | MXDEL | MXDE |
| LOCAL | LOCA | MXLIST | MXLI |
| NLIST | NLIS | RLIST | RLIS |
| NMOVE | NMOV | SCLIST | SCLI |
| NSCALE | NSCA | THICK | THIC |
| ECOMPR | ECOM | FLIST | FLIS |
| ELIST | ELIS | FSCALE | FSCA |
| GLIST | GLIS | MMDEL | MMDE |
| GROUP | GROU | MMLIST | MMLI |
| LLIST | LLIS | MRGEN | MRGE |
| MLIST | MLIS | MRDEL | MRDE |
| MTYPE | MTYP | MRSCALE | MRSC |
| MSCALE | MSCA | MTDEL | MTDE |
| MTGEN | MTGE | MTLIST | MTLI |
| MTSCALE | MTSC | MUNITS | MUNI |

- Commands of 4 characters or fewer are typed in full: CSYS, D, FILL, INT, LOC, N, NDEL, NGEN, NMED, SDEL, DELL, DELM, DELR, E, EDEL, EGEN, EINT, GDEL, GTIT, KI, KJ, L, M, MACT, MSET, MXI, MXM, MXR, R, RACT, RSET, SC, F, FDEL, MM, MR, MT.
- **No underline is printed** for INTLIST, LMOVE, SLIST, DELSC and ETYPE, although all five are longer than 4 characters and are listed as legacy PREP commands. *(inferred)* Also accept INTL, LMOV, SLIS, DELS and ETYP (§12, Q1).
- *(inferred)* Lookup algorithm:
  1. Upper-case the token.
  2. If it equals a full command name (legacy or UI-only), use that command. This rule handles GROUPMAT vs GROUP, ETYPEGEN vs ETYPE, and INTGEN/INTCOUNT vs INT/INTLIST.
  3. Otherwise, if the token is exactly 4 characters and equals the abbreviation of a legacy command, use that command.
  4. Otherwise, reject the token.
  - There are no 4-letter clashes inside this command set. LOC, MXI and INT are exact 3-letter names and do not collide with LOCA, MXLI or INTL.

### 1.4 Ranges `<a>,<b>,<inc>` [IO] *(inferred where not stated)*
- Many commands act on a numeric range: start, end, step. It is used for nodes, elements, groups, materials, layers, properties and coordinate systems.
- Meaning: the inclusive set {a, a+inc, a+2·inc, …} up to and including b.
- Defaults are given per command. The most common pattern is: end = start, step = 1.
- *(inferred)*
  - If `inc ≤ 0`, use 1.
  - If `b < a`, swap them.
  - IDs in the range that are **not defined are skipped silently** by modify/delete/list commands.
  - Report how many items were affected.

### 1.5 "Last defined" bookkeeping [IO] *(inferred, needed by FILL/NGEN/FSCALE/...)*
Several defaults refer to the most recently defined entity:
- FILL: "latest two consecutively defined nodes".
- NGEN: "second to last node" and "last node to be defined".
- FSCALE, MSCALE, MTSCALE, MRSCALE: "last 2 defined" forces, moments or masses.
- LLIST, MLIST, RLIST, SCLIST: `<m2>` defaults to the "last defined".

Keep an ordered **definition history** per entity kind:
- nodes;
- forces, moments, translational masses and rotational masses (each with its own list);
- materials, soil layers, real properties, spring properties, matrix properties.

A definition (or redefinition) moves the ID to the end of its list. Deleting an entity removes it from the list. For the table lists, "last defined" means the **highest index defined**. See §12, Q2 for the alternative reading.

### 1.6 Units [CORE] *(xref 05b §3.4)*
- The program does **not** convert units. Use one consistent unit system for the whole model.
- The acceleration of gravity `g` (HOUSE `<gravity>` or GRAVITY command) turns weights into masses:
  - mass density ρ = γ / g, where γ is the specific weight in force per volume;
  - nodal or matrix mass = weight / g.
- Typical g: 32.2 ft/s² (British) or 9.81 m/s² (SI).
- *(inferred)* Store raw user values. Do all weight→mass and modulus conversions when the analysis files are written (AFWRITE/CHECK), because g may be set after the materials.

### 1.7 In-memory data model [CORE][IO] *(inferred structure; field meanings are from the manual)*

```
Model
  csys_active: int = 0                       # CSYS
  coord_systems: {ns:int -> CoordSys}        # LOC, LOCAL; ns >= 1 (0 is global, implicit)
  nodes: {nd:int -> Node}
  groups: {ng:int -> Group}; group_active: int|None
  materials:   {nm -> Material}              # M
  soil_layers: {nm -> SoilLayer}             # L
  real_props:  {nm -> RealProp}              # R      (BEAMS)
  spring_props:{nm -> SpringProp}            # SC     (SPRING)
  matrix_props:{p  -> MatrixProp}            # MXR/MXI/MXM (GENERAL)
  mact: int = 1  (inferred default)          # MACT
  ract: int = 1  (inferred default)          # RACT
  forces:  {n -> NodalLoad}                  # F
  moments: {n -> NodalLoad}                  # MM
  tmass:   {n -> (mx,my,mz)}                 # MT
  rmass:   {n -> (mxx,myy,mzz)}              # MR
  mass_units: {n -> 0|1}                     # MUNITS (0=mass, 1=weight)
  options (from MOPT, xref 05a): incomp, matrix_units, overwrite_mass, overwrite_force

CoordSys:   ns, type(0=Cartesian), origin[3] (global), R[3x3] (columns = local x,y,z unit vectors in global)
Node:       nd, xyz[3] (in csys), csys:int, bc[6] (0 free / 1 fixed; order UX,UY,UZ,ROTX,ROTY,ROTZ),
            flags: {interaction, intermediate, interface, internal} (booleans),
            legacy: pile, control (ignored)
Group:      ng, type in {1 SOLID,2 BEAMS,3 SHELL,4 PLANE,5 TSHELL,7 SPRING,9 GENERAL}, title:str,
            elements: {ne -> Element}  (ne = 1..N, gap-free when written)
Element:    ne, nodes[<=8] (0 = not given), mat:int, prop:int, etype:int(0/1/2),
            eint:int (SOLID: 0/1/2 ; TSHELL: 0/1), thick:float, ki[6], kj[6]
Material:   nm, val1, val2, weight, pdamp, sdamp, mtype(1/2/3)
SoilLayer:  nm, thick, weight, vp, vs, pdamp, sdamp
RealProp:   nm, axial, shear2, shear3, tors, flex2, flex3
SpringProp: nm, scx, scy, scz, scxx, scyy, sczz, damp
MatrixProp: p, KR[12x12], KI[12x12], M[12x12]    (symmetric, stored as upper triangle)
NodalLoad:  n, factor[3], arrival[3]
```

---

## 2. Coordinate systems and transformations [CORE]

### 2.1 Conventions *(manual, §9.3.10)*
- The global system is Cartesian and right-handed, with **Z vertical and positive upward**.
- X and Y follow the right-hand rule.
- For **1D and 2D analyses**, the program **ignores the Y coordinate** of the nodes. A 2D model lies in the **X–Z plane**.
- **Local coordinate systems** are numbered `ns ≥ 1`. `ns = 0` is the global system.
- Only **Cartesian** local systems (`type = 0`) are implemented. Cylindrical (`type = 1`) is "not included in this version".

### 2.2 Node storage and transformation *(manual §9.3.4 + inferred)*
- A node is stored with its coordinates **in the system that was active when it was defined**, plus that system's number.
- GLOBAL turns stored coordinates into global ones and sets the node's csys to 0.
- The same transformation is applied **automatically when the analysis files are written** (AFWRITE).
- GLOBAL **does not deactivate** the active local system. A later N still uses the active CSYS.
- Local to global, for a node with local coordinates `x_l` in system s:
  `X_g = O_s + R_s · x_l`
- Global to local: `x_l = R_sᵀ · (X_g − O_s)`
- *(inferred)* Generation commands (NGEN, FILL, LMOVE, NMOVE, NSCALE, NMED) work in the **active** coordinate system, as ANSYS does:
  - convert each source node to the active system;
  - apply the operation;
  - store the result in the active system.
  - Because only Cartesian systems exist, FILL and NMED give the same geometric result in any system. NGEN, LMOVE, NMOVE and NSCALE do not, because the increments or factors are applied along the active axes. See §12, Q3.

---

## 3. Node commands (§9.3)

### 3.1 CSYS — activate a coordinate system [CORE][IO]
`CSYS,<ns>`

| # | Arg | Req | Default | Meaning |
|---|---|---|---|---|
| 1 | ns | yes | — | Local system to activate. **0 = global.** |

- Sets `csys_active = ns`. All later N commands interpret coordinates in this system.
- *(inferred)* It is an error if `ns ≠ 0` and system ns is not defined.
- Example: `CSYS,1` (activate local system 1), then `CSYS,0` (back to global).

### 3.2 D — boundary conditions (fixed/free DOFs) [CORE][IO]
`D,<n1>,<n2>,[<inc>],[<val>],<label1>,[<label2>,…<label6>]`

| # | Arg | Req | Default | Meaning |
|---|---|---|---|---|
| 1 | n1 | yes | — | first node |
| 2 | n2 | yes (syntax) | *(inferred lenient: n1)* | last node |
| 3 | inc | no | **1** | node step |
| 4 | val | no | **0** | **0 = free DOF, 1 = fixed DOF** |
| 5–10 | label1…label6 | ≥1 | — | DOF labels: `UX`, `UY`, `UZ`, `DISP`, `ROTX`, `ROTY`, `ROTZ`, `ROT`, `ALL` |

- Label expansion *(inferred meaning of the group labels)*:
  - `DISP` = UX, UY, UZ
  - `ROT` = ROTX, ROTY, ROTZ
  - `ALL` = all six DOFs
- Action: for each existing node in the range and each DOF in the expanded labels, set `bc[dof] = val`.
- D is meant to be applied **after** the nodes exist. Node generation can therefore ignore boundary conditions. *(inferred)* NGEN, FILL, LMOVE and NMOVE do **not** copy the bc codes of their source nodes.
- **Only 0/1 codes are supported.** There are no prescribed non-zero displacements; support motions come from the seismic input.
- *(inferred)* DOF labels always refer to **global** axes. There are no nodal rotated coordinate systems in SASSI.
- Effect on loads and masses *(manual)*:
  - Concentrated masses (MT, MR), forces (F) and moments (MM) must sit on **free** DOFs.
  - Values on a fixed DOF are **ignored** (CHECK Warnings 5, 6, 10, 11).
  - A fixed DOF has no translation or rotation.
- Fixing unwanted DOFs makes the equation system smaller.
- **Element DOF table** *(figure, §9.3.2)*. Bullets show the DOFs each element type defines at its nodes:

| Element type | X | Y | Z | XX | YY | ZZ |
|---|---|---|---|---|---|---|
| SOLID | • | • | • | | | |
| BEAM(S) | • | • | • | • | • | • |
| SHELL / TSHELL | • | • | • | • | • | • |
| PLANE | • | | • | | | |
| SPRING | • | • | • | • | • | • |
| GENERAL | • | • | • | • | • | • |

- A node shared by several element types gets the **union** of their DOFs. Example: a node shared by a BEAM and a SOLID has all 6 DOFs ("beam governs").
- *(inferred, recommended)* When assembling, compute each node's active DOFs as: (union of DOFs of the attached elements) minus (DOFs fixed by D). Then:
  - **Auto-fix**, with an informational message, any DOF that no attached element defines (for example the rotations of SOLID-only nodes).
  - For DOFs that are defined but carry no stiffness, such as the SHELL drilling rotation, keep the manual's workflow: the user runs FIXROT, FIXSHLROT, FIXSLDROT or FIXSPRROT (UI commands, xref §9.7.7–9.7.10). Warn if this is not done.
  - See §12, Q4.
- Symmetric structures under symmetric loading can be cut to a half or quarter model by fixing the DOFs on the symmetry planes (see SYMM, xref §9.2.39).
- **Interaction nodes:** the soil impedance acts only on their **translational** DOFs. Structural rotations are passed to the soil through the nodal translations (xref line 805; 01_capabilities §5.2).
  - *(inferred)* A shell or beam attached to an interaction node keeps its rotational DOFs there as ordinary structural DOFs. They are simply not coupled to the soil.
  - A node that is **fully fixed and flagged as interaction is CHECK Error 124**.
- Examples:
  - `D,1,121,1,1,ALL` fixes all 6 DOFs of nodes 1..121.
  - `D,200,300,10,1,ROT` fixes the rotations of nodes 200, 210, …, 300.
  - `D,5,5,,1,UX,UZ` fixes UX and UZ of node 5.
  - `D,5,5,,0,UX` frees UX of node 5 again.
  - Note that `D,5,5,,,ALL` **frees** all DOFs, because val defaults to 0.

### 3.3 FILL — generate a line of nodes by interpolation [CORE][IO]
`FILL,[<n1>],[<n2>],[<nr>]`

| # | Arg | Default | Meaning |
|---|---|---|---|
| 1 | n1 | second-to-last node of the definition history | start node (must exist) |
| 2 | n2 | last node of the definition history | end node (must exist) |
| 3 | nr | `n2 − n1 − 1` | number of nodes to create between n1 and n2 |

- Positions are equally spaced on the straight line: `X_k = X_1 + k/(nr+1)·(X_2 − X_1)`, k = 1..nr.
- Numbering *(inferred, ANSYS-like)*:
  - node increment `d = (n2 − n1)/(nr + 1)`;
  - new nodes are `n1 + k·d`, k = 1..nr;
  - if d is not an integer, use `n1 + k` and warn;
  - it is an error if any new ID ≥ n2 or collides with a node that exists outside this FILL.
  - With the default nr, d = 1, which gives nodes n1+1 … n2−1.
- Example: `N,1,0,0,0` / `N,5,40,0,0` / `FILL` creates nodes 2, 3, 4 at x = 10, 20, 30.

### 3.4 GLOBAL — transform nodes to global coordinates [CORE][IO]
`GLOBAL,<n1>,<n2>,<inc>` (abbr. `GLOB`)

| # | Arg | Default (inferred lenient) | Meaning |
|---|---|---|---|
| 1 | n1 | first node | start |
| 2 | n2 | last node | end |
| 3 | inc | 1 | step |

- For each node in the range with `csys ≠ 0`: `xyz ← O + R·xyz`, then `csys ← 0`.
- AFWRITE does the same for **all** nodes.
- GLOBAL does **not** change `csys_active`.
- Example: `GLOBAL,1,500,1`.

### 3.5 INT — set or reset node classification flags [CORE][IO]
`INT,<n1>,<n2>,[<inc>],<set>,[<code>]`

| # | Arg | Req | Default | Meaning |
|---|---|---|---|---|
| 1 | n1 | yes | — | first node |
| 2 | n2 | yes | — | last node |
| 3 | inc | no | 1 *(inferred)* | step |
| 4 | set | yes | — | **1 = set** the flag, **0 = reset** (clear) it |
| 5 | code | no | **0** *(inferred)* | **0 = interaction**, **1 = intermediate**, **2 = interface**, **3 = internal** *(inferred: the manual text and the PDF page are both cut off after "internal (<")* |

- The flags are independent booleans, so a node can carry several. INTLIST filters on them separately.
- **Interaction nodes [CORE]** are the excavated-soil nodes where the soil impedance is applied (xref 05b §1.1, line 5565 and lines 3008–3024). Rules:
  - They belong to the **excavated soil** model, not to the structure. The only exception is nodes shared by the structure basement and the excavation boundary.
  - Every interaction node below the ground surface must lie **on a soil-layer interface**.
  - For embedded models they must be numbered in **ascending order from the deepest excavation level up to the ground surface**.
  - **FV:** all excavated-volume nodes are interaction nodes. **FI and FFV:** only a subset. **Surface foundations:** the ground-surface nodes only.
  - The UI command **INTGEN** (xref §9.7.19) generates them automatically. INT is the manual way.
  - All nodes on the excavation / far-field interface must be interaction nodes (xref line 884).
- **Intermediate / interface / internal [ADV, legacy]:** the manual gives no computational meaning for these in V3 (§12, Q5). *(inferred)* Store them, list them and write them out. They have no effect on the solution in this implementation.
- Examples:
  - `INT,1,121,1,1,0` (or `INT,1,121,,1`) flags nodes 1..121 as interaction nodes.
  - `INT,50,50,1,0,0` clears the interaction flag of node 50.

### 3.6 INTLIST — list classified nodes [UI]
`INTLIST,[<n1>],[<n2>],[<step>],[<c1>],[<c2>],[<c3>],[<c4>]` (no underline printed; see §1.3)

| # | Arg | Default | Meaning |
|---|---|---|---|
| 1 | n1 | first node | start |
| 2 | n2 | last node | end |
| 3 | step | 1 | step |
| 4 | c1 | 0 | 1 = list interaction nodes |
| 5 | c2 | 0 | 1 = list intermediate nodes |
| 6 | c3 | 0 | 1 = list interface nodes |
| 7 | c4 | 0 | 1 = list internal nodes |

- If all of c1..c4 are 0, every node that carries any of the four flags is listed.
- Example: `INTLIST,,,,1` lists all interaction nodes.

### 3.7 LMOVE — generate a node list by translation [CORE][IO]
`LMOVE,[<dx>],[<dy>],[<dz>],<nd>,<l1>,[<l2>,…,<l15>]` (no underline printed)

| # | Arg | Default | Meaning |
|---|---|---|---|
| 1–3 | dx, dy, dz | 0 | translation, in the active system *(inferred)* |
| 4 | nd | — | first new node number |
| 5… | l1 … l15 | ≥ 1 required | source nodes (an explicit list, not a range) |

- New node `nd + (i−1)` = node `l_i` + (dx, dy, dz), for i = 1..k.
- The syntax line allows **15** list entries, but the text says `<l9>`. *(inferred)* Accept up to 15, or more (§12, Q6).
- The list must contain at least one node.
- *(inferred)* Error if a source node does not exist. Overwrite an existing target node, with a warning.
- Example: `LMOVE,0,0,10,101,1,2,3,4` creates nodes 101..104 = nodes 1..4 shifted +10 in Z.

### 3.8 LOC — local coordinate system from Euler angles [CORE][IO]
`LOC,<ns>,<type>,<x0>,<y0>,<z0>,<txy>,<tyz>,<txz>`

| # | Arg | Meaning |
|---|---|---|
| 1 | ns | system number (≥ 1) |
| 2 | type | 0 = Cartesian; 1 = cylindrical, z is the rotation axis (**not included in this version**, so reject it) |
| 3–5 | x0, y0, z0 | origin in **global** coordinates |
| 6–8 | txy, tyz, txz | rotation angles in **degrees** |

- *(inferred, ANSYS LOCAL convention; the manual does not give the order or signs)* Apply three successive rotations about the moving axes:
  1. `txy` about Z, turning X toward Y.
  2. `tyz` about the new X, turning Y toward Z.
  3. `txz` about the new Y, turning Z toward X.
- The rotation matrix is `R = Rz(txy) · Rx(tyz) · Ry(txz)`. Its columns are the local axes expressed in global coordinates. With `c = cos` and `s = sin`:
  - `Rz(a) = [[c,−s,0],[s,c,0],[0,0,1]]`
  - `Rx(b) = [[1,0,0],[0,c,−s],[0,s,c]]`
  - `Ry(g) = [[c,0,s],[0,1,0],[−s,0,c]]`
- Defining an existing ns again replaces that system. *(inferred)* Nodes already stored in system ns keep their local numbers, so their global position moves. Warn about this.
- See LOCAL for the node-based definition.
- Example: `LOC,1,0,10,0,0,45,0,0`, then `CSYS,1`, `N,100,5,0,0` puts node 100 at global (13.536, 3.536, 0).

### 3.9 LOCAL — local coordinate system from three nodes [CORE][IO]
`LOCAL,<ns>,<type>,<n1>,<n2>,<n3>` (abbr. `LOCA`)

| # | Arg | Meaning |
|---|---|---|
| 1 | ns | system number |
| 2 | type | 0 = Cartesian (1 = cylindrical is not implemented) |
| 3 | n1 | origin node |
| 4 | n2 | node on the +x axis (x = direction n1→n2) |
| 5 | n3 | node in the positive x–y quadrant (local y > 0) |

Algorithm, using global coordinates of n1, n2, n3:
```
a = X(n2) - X(n1);  b = X(n3) - X(n1)
ex = a/|a|
ez = (a x b)/|a x b|          # error if n1,n2,n3 collinear or coincident
ey = ez x ex                  # (manual text says "X x Z"; that gives -Y. ez x ex is the
                              #  right-handed choice that puts n3 at y>0; see §12 Q7)
O = X(n1);  R = [ex ey ez]    # columns
```
- Example: n1 = (0,0,0), n2 = (0,1,0), n3 = (−1,0,0) gives ex = +Y, ey = −X, ez = +Z. This is the same system as `LOC,ns,0,0,0,0,90,0,0`.

### 3.10 N — define a node [CORE][IO]
`N,<nd>,[<x>],[<y>],[<z>]`

| # | Arg | Default | Meaning |
|---|---|---|---|
| 1 | nd | — | node number (positive integer) |
| 2–4 | x, y, z | 0 | coordinates in the **active** system |
| 5–6 | pile, control node | 0 | **deprecated** legacy fields [ADV]. Accept them and ignore them. |

- An existing nd is redefined: its coordinates and csys are replaced *(inferred)*. Its bc codes and flags are kept *(inferred)*.
- **Node numbering rules [CORE]** *(manual)*:
  - Node numbers may be arbitrary.
  - However, to minimize storage and block operations, to allow restart with a **new superstructure**, and for the **incoherence** option, the nodes **at or below ground surface must be numbered first, layer by layer, starting from the bottom**.
  - Z must be vertical upward.
  - For 1D/2D analysis Y is ignored.
  - *(xref line 5608)* Recommended practice is bottom-up numbering from the foundation level, then layer by layer to the surface.
- *(xref Warning 1)* Gaps in node numbering are allowed in the model. AFWRITE writes the missing nodes to the analysis files with all DOFs fixed.
- Example: `N,1,0,0,-20` / `N,2,10,0,-20`.

### 3.11 NDEL — delete nodes [IO]
`NDEL,<n1>,[<n2>],[<inc>]`. Defaults: n2 = n1, inc = 1.
- Removes the nodes in the range.
- *(inferred)*
  - Elements that reference a deleted node are kept. CHECK then reports Error 41.
  - Forces, moments and masses on deleted nodes are deleted together with the node.
  - Interaction flags go with the node.
  - Warn with a count of the dangling references.
- Example: `NDEL,10,20,2`.

### 3.12 NGEN — generate nodes by copying a pattern [CORE][IO]
`NGEN,[itim],[step],[n1],[n2],[inc],[dx],[dy],[dz]`

| # | Arg | Default | Meaning |
|---|---|---|---|
| 1 | itim | 1 | number of new sets, not counting the original pattern |
| 2 | step | `n2 − n1 + 1` | node-number increment between sets |
| 3 | n1 | second-to-last node defined | first pattern node |
| 4 | n2 | last node defined | last pattern node |
| 5 | inc | 1 | step inside the pattern |
| 6–8 | dx, dy, dz | 0 | coordinate increment per set (active system, *inferred*) |

Algorithm:
```
for k in 1..itim:
  for n in range(n1, n2, inc) if n exists:
     new = n + k*step
     xyz(new) = xyz_active(n) + k*(dx,dy,dz);  csys(new) = csys_active
```
- *(inferred)*
  - bc codes and flags are not copied.
  - An existing target node is overwritten, with a warning.
  - "Last node defined" means the most recent entry in the definition history (§1.5).
- Example: after nodes 1..5 at z = 0, `NGEN,3,5,1,5,1,0,0,10` creates 6–10 at z = 10, 11–15 at z = 20 and 16–20 at z = 30.

### 3.13 NLIST — list nodes [UI]
`NLIST,[<n1>],[<n2>],[<inc>]` (abbr. `NLIS`). Defaults: first node, last node, 1.
- *(inferred)* Columns: node number, x, y, z (as stored), csys, the six bc codes, and the flags (I/M/F/N).

### 3.14 NMED — node at the average of other nodes [CORE][IO]
`NMED,<nd>,<n1>,[<n2>,…,<n8>]`
- Defines node nd at the arithmetic mean of the coordinates of 1 to 8 listed nodes. At least one node is required.
- *(inferred)* Average the **global** coordinates. Store the result in the active system.
- Example: `NMED,500,1,2,7,6` puts node 500 at the centroid of nodes 1, 2, 7 and 6.

### 3.15 NMOVE — generate a node list by scaling [CORE][IO]
`NMOVE,[<dx>],[<dy>],[<dz>],<nd>,<l1>,[<l2>,…,<l15>]` (abbr. `NMOV`)
- New node `nd + (i−1)` has coordinates `(x·dx, y·dy, z·dz)` of node `l_i`. **The defaults for dx, dy and dz are 1.**
- Scaling is about the origin of the active system *(inferred)*.
- The list length issue is the same as for LMOVE (15 in the syntax, 9 in the text).
- *(inferred)* A factor of exactly 0 is used as given; the manual does not apply the "0 → 1" rule here. Warn, because 0 collapses the coordinate. See §12, Q8.
- Example: `NMOVE,2,1,1,201,1,2,3` creates nodes 201..203 with doubled x.

### 3.16 NSCALE — scale node coordinates in place [CORE][IO]
`NSCALE,<n1>,<n2>,[<inc>],[<sfx>],[<sfy>],[<sfz>]` (abbr. `NSCA`)
- For each node in the range: `x ← x·sfx`, `y ← y·sfy`, `z ← z·sfz`. inc defaults to 1.
- **A factor of 0.0 (including a blank one) is replaced by 1.0.**
- *(inferred)* Coordinates are scaled in the active system.
- Example: `NSCALE,1,500,1,0.3048,0.3048,0.3048` converts ft to m.

### 3.17 SDEL — delete coordinate systems [IO]
`SDEL,<s1>,[<s2>],[<inc>]`. Defaults: s2 = s1, inc = 1.
- *(inferred)*
  - Before a system is deleted, every node stored in it is converted to global coordinates.
  - If the active system is deleted, `csys_active` becomes 0.
  - System 0 cannot be deleted.

### 3.18 SLIST — list coordinate systems [UI]
`SLIST,[<s1>],[<s2>],[<inc>]` (no underline printed). Defaults: first system, last system, 1.
- *(inferred)* Columns: ns, type, origin, the 3 Euler angles (for LOC) or the defining nodes (for LOCAL), and the axis unit vectors.

---

## 4. Finite-element library [CORE]

### 4.1 Group types *(figure, GROUP table, §9.4.14)*

| Type no. | Type string | Nodes given in `E` | Description |
|---|---|---|---|
| 1 | `SOLID` | 8 (7, 6 or 5 distinct nodes, obtained by repeating node numbers) | 3D solid |
| 2 | `BEAMS` | 3 (I, J, K) | 3D beam |
| 3 | `SHELL` | 4 or 3 | 3D plate/shell, thin, **Kirchhoff** theory |
| 5 | `TSHELL` | 4 or 3 | 3D plate/shell, thick, **Mindlin–Reissner** theory |
| 4 | `PLANE` | 4 or 3 | 2D **plane-strain** solid |
| 7 | `SPRING` | 2 | 3D spring (translations and/or rotations) |
| 9 | `GENERAL` | 3 (I, J, K: local axes) or 2 (I, J: global axes) | 3D general stiffness/mass ("GM") element |

- Type numbers 6 and 8 are not used.
- Only SOLID (3D) and PLANE (2D) can model the **excavated soil**.
- *(xref line 5570; §1.5.3)* Further facts about the library:
  - SOLID: 8 nodes; optional **9 incompatible displacement modes** for structural SOLIDs (global MOPT `<incomp>`: 0 include, 1 suppress). Never used for excavated-soil SOLIDs.
  - TSHELL: HOUSE adds the drilling stiffness automatically.
  - SHELL: has no drilling stiffness. Use FIXROT or FIXSHLROT.
  - GENERAL: 6 DOFs per node. The super-element matrices are complex.
  - Option AA converts ANSYS MATRIX50 elements to GENERAL elements.
- *(xref line 965)* Mass matrices:
  - **50% lumped + 50% consistent** for all elements,
  - **except beams (consistent)** and **plates/shells (lumped)**.
- *(xref line 969)* Damping uses **complex moduli**, which gives frequency-independent damping that can differ from element to element. The complex factor convention is in §6.6.

### 4.2 Element data assignable per group *(figure, §9.4.14 second table)*

| Group type | Element data | Command |
|---|---|---|
| SOLID | material / soil-layer index | MSET (or MACT) |
| | element type | ETYPE |
| | *(also)* integration order | EINT |
| BEAMS | material index | MSET |
| | real-property index | RSET (or RACT) |
| | I-node release code | KI |
| | J-node release code | KJ |
| SHELL / TSHELL | material index | MSET |
| | thickness | THICK |
| | element type | ETYPE |
| | *(TSHELL)* integration type | EINT |
| PLANE | material / soil-layer index | MSET |
| | element type | ETYPE |
| SPRING | spring-property index | RSET |
| GENERAL | matrix-property index | RSET |

### 4.3 Grouping rules [CORE][IO] *(manual WARNING, §9.4.14)*
- All elements must be **grouped by type**: one group holds one element type.
- **No gaps in element numbering** inside a group are allowed in the analysis model. Use ECOMPR (CHECK Error 6).
- Element numbers **start at 1 in each group**.
- One element type may use several groups. Example: structural bricks in one group and excavated-soil bricks in another.
- *(xref line 5625; recommended practice, a strict rule if stress contour plots are wanted)*
  - Build the excavation volume from **one group per embedment layer**: a horizontal layer of elements tied to one far-field soil layer.
  - Number these groups **from the surface down to the foundation level**.
  - Number the excavation elements top-down.
  - Each excavation group uses the **L layer number equal to its far-field layer**.
- Empty groups are skipped by AFWRITE (Warning 7).
- *(inferred)* Gaps in group numbers are allowed. The GCOM UI command removes them.

### 4.4 SOLID (type 1) [CORE]
- **Nodes:** 8 node numbers are always given.
  - Nodes 1-2-3-4 are one face, numbered in sequence around it. Nodes 5-6-7-8 are the opposite face, with node 4+i opposite node i.
  - **Degenerate shapes (Figure 9.1)** come from repeating node numbers:
    - 7 distinct nodes: repeat 7 and 8 (`n7 = n8`).
    - 6 distinct nodes (wedge with the top edge collapsed): `n5 = n6` and `n7 = n8`.
    - 5 distinct nodes (pyramid): `n5 = n6 = n7 = n8`.
- **Natural axes (Figure 9.1, inferred from the arrows):**
  - ξ points toward face 2-3-7-6.
  - ζ points toward face 5-6-7-8.
  - η is drawn pointing out of face 1-2-6-5.
  - *(inferred)* Use the standard isoparametric mapping with nodes at (ξ,η,ζ) = (−1,−1,−1), (1,−1,−1), (1,1,−1), (−1,1,−1), (−1,−1,1), (1,−1,1), (1,1,1), (−1,1,1). The internal orientation does not affect the results; only the face/opposite-face topology does.
  - Degenerate elements work through the repeated nodes in the trilinear mapping. The Jacobian must stay non-singular at the Gauss points, which holds for the collapsed patterns above.
- **DOFs:** UX, UY, UZ at each node.
- **Material:**
  - Structural (ETYPE 1, or ETYPE 0 above ground): **M** table.
  - Excavated soil (ETYPE 2, or ETYPE 0 below ground): **L** table (§6.2).
- **Integration:** set by EINT. 0 = rectangular, 1 = skewed, 2 = extremely distorted. The default is "2×2" (2×2×2 Gauss).
  - *(inferred)* 0 → 2×2×2, 1 → 3×3×3, 2 → 4×4×4 points (SAP IV heritage) (§12, Q9).
- **Incompatible modes:** 9 Wilson modes when MOPT `<incomp>` = 0 (include) and the element is structural. The 05a screenshot default is Suppress.
- **Mass:** 50% lumped + 50% consistent, with ρ = γ/g.
- **Stiffness:** linear isotropic elasticity with complex Lamé constants from M\* and G\* (§6.6).

### 4.5 BEAMS (type 2) [CORE]
- **Nodes:** I, J (end nodes) and **K**, a geometric reference point.
  - K may be any other node, but it **must not lie on local axis 1** (CHECK Error 9 if I, J, K are collinear).
  - K is a reference only. It contributes no DOFs. If K is used by nothing else it is reported as an unused node (Warning 4) and fixed by AFWRITE.
  - **KINT** (xref §9.8.6) warns when K-nodes are also interaction nodes.
- **Local axes (Figure 9.4):**
  ```
  e1 = (X_J - X_I)/|X_J - X_I|                    # axis 1 along I->J
  v  = X_K - X_I
  e2 = (v - (v·e1) e1)/|v - (v·e1) e1|             # axis 2 in plane I-J-K, toward K
  e3 = e1 x e2                                      # axis 3 completes the right-handed triad
  ```
- **DOFs:** 6 per node (UX, UY, UZ, ROTX, ROTY, ROTZ).
- **Section (R table):** axial area A, shear areas As2 and As3, torsion constant J, inertias I2 and I3 (§6.3).
  - **I2 is about local axis 2. I3 is about local axis 3.**
  - *(inferred, standard)*
    - Bending in the 1–2 plane (deflection along axis 2) uses E·I3 and shear area As2.
    - Bending in the 1–3 plane uses E·I2 and As3.
    - Shear-deformation parameters: `φ2 = 12·E·I3/(G·As2·L²)` and `φ3 = 12·E·I2/(G·As3·L²)`.
    - If a shear area is 0, shear deformation is ignored in that plane (φ = 0). The manual recommends this when shear deformation is not wanted.
- **End releases (KI, KJ):** six 0/1 digits for the end forces **P1, P2, P3, M1, M2, M3** (local axes) at node I or J.
  - **1 = released** (the force is zero: hinge or roller).
  - CHECK **Error 10** if the same component is released at both I and J.
  - *(inferred)* Remove released components by static condensation of the 12×12 local stiffness and mass before transforming to global.
- **Material:** M table, using E, G (or ν), γ and damping (§6.1, §6.6).
- **Mass:** consistent. *(inferred)*
  - Standard 12×12 consistent mass from ρ·A·L, with polar inertia ρ·J for torsion (ρ·(I2+I3) if J is not a polar value; §12, Q10).
  - Rotary inertia is neglected.
- **Transformation:** `K_g = Tᵀ K_l T`, with `T = blockdiag(Λ,Λ,Λ,Λ)`, where the rows of Λ are e1, e2 and e3.

### 4.6 SHELL (type 3) and TSHELL (type 5) [CORE]
- **Nodes:** I, J, K, L, in sequence and **counter-clockwise** around the element. Omit L for a triangle.
- **Local axes (Figure 9.3):**
  ```
  m_LI = (X_L+X_I)/2; m_JK = (X_J+X_K)/2; m_IJ = (X_I+X_J)/2; m_KL = (X_K+X_L)/2
  x' = (m_JK - m_LI)/|...|
  z' = (x' x (m_KL - m_IJ))/|...|        # normal to x' and to the line IJ-KL
  y' = z' x x'                            # right-handed; used for the resultant forces
  ```
  - *(inferred, triangles)* Use L ≡ K in these formulas, so the K-L midpoint is node K.
  - With I→J→K→L counter-clockwise as seen from +z', z' points toward the viewer.
- **DOFs:** 6 per node.
- **Thickness:** set per element with THICK. *(inferred)* Error at CHECK if the thickness is ≤ 0 or was never set.
- **Material:** M table. **`<pdamp>` must equal `<sdamp>`** for materials used by SHELL elements *(manual, §9.4.20)*. *(inferred)* The same applies to TSHELL; CHECK should warn otherwise.
- **SHELL formulation:** thin plate, Kirchhoff theory. It has **no drilling stiffness**, so FIXROT or FIXSHLROT must be used (default soft spring stiffness 10, xref §9.7.8, 9.7.10).
  - *(inferred)* Flat facet = plane-stress membrane + Kirchhoff (DKQ/DKT) bending (§12, Q11).
- **TSHELL formulation:** thick plate, Mindlin–Reissner theory, with drilling stiffness added automatically by HOUSE.
  - **EINT** for TSHELL: 0 = **Reduced** (1×1 Gauss for both plate bending and transverse shear; the default), 1 = **Selective** (2×2 for bending, 1×1 for transverse shear).
- **Mass:** lumped, `ρ·t·A/n_nodes` on each translational DOF. *(inferred)* Rotational lumped inertia is ρ·t³/12·A/n_nodes, or zero (§12, Q10).
- **ETYPE:**
  - 0 or 1: structural shell.
  - 2: **embedded/buried shell**, a shell inside or on the soil whose nodes become interaction nodes through INTGEN. These shells are **not** part of the excavation volume (xref §9.7.6, 9.7.19).
- *(xref §9.19)* THSHLSTR and THSHLSMH control TSHELL stress output (STRESS module).

### 4.7 PLANE (type 4) [CORE]
- **Nodes:** I, J, K, L counter-clockwise (Figure 9.2). Omit L for a triangle.
  - The element lies in the **X–Z plane**, because Y is ignored in 2D.
  - *(inferred)* The manual does not say which side of the X–Z plane "counter-clockwise" is viewed from. Accept either orientation. If the signed area (the Jacobian) is negative, reorder the nodes internally and warn (§12, Q12).
- **DOFs:** UX, UZ.
- **Formulation:** 2D **plane strain**. *(inferred)* Use **unit thickness**, i.e. per unit length out of plane, consistent with 2D SSI (POINT2).
- **Material:**
  - Structural: M table.
  - Excavated soil (ETYPE 2, or ETYPE 0 below ground): L table.
- **Integration:** *(inferred)* 2×2 Gauss. EINT applies to SOLID and TSHELL only.
- **Mass:** 50% lumped + 50% consistent.

### 4.8 SPRING (type 7) [CORE]
- **Nodes:** I, J.
- **Property:** SC table: **six uncoupled spring constants in the global directions** (kx, ky, kz, kxx, kyy, kzz) and one damping ratio.
- Element stiffness. For each global DOF d in {UX..ROTZ}:
  ```
  K_e[d_I,d_I] += k_d* ; K_e[d_J,d_J] += k_d* ; K_e[d_I,d_J] -= k_d* ; K_e[d_J,d_I] -= k_d*
  k_d* = k_d · c(β)      (complex factor, §6.6)
  ```
- No mass.
- The constants are added **directly to the global stiffness matrix**, so they must be given in **global XYZ**. The local system and node orientation are irrelevant.
- *(inferred)* A grounded (one-node) spring is made by connecting J to a node with all DOFs fixed. Whether `J = 0` is accepted is not stated (§12, Q13).
- *(xref §1.5.3)* For Option NON nonlinear springs, split 3D springs into 1D spring groups [ADV].

### 4.9 GENERAL (type 9) [CORE][ADV]
- **Nodes:**
  - **I and J only:** the matrix terms are in **global** coordinates.
  - **I, J and K:** the terms are in the **local** system defined by I, J, K (figure).
  - K is a reference point only.
- **DOF numbering (figure):**
  - 1, 2, 3 = translations of I along axes 1/2/3 (or X/Y/Z). 4, 5, 6 = rotations of I about those axes.
  - 7, 8, 9 = translations of J. 10, 11, 12 = rotations of J.
  - Matrix property: three symmetric **12×12** matrices, entered as the **upper triangle, row by row** (MXR, MXI, MXM):
    - K_R = real stiffness;
    - K_I = imaginary stiffness;
    - M = mass, or weight when MOPT `<matrix>` = 1.
- Element complex stiffness: `K* = K_R + i·K_I`. *(inferred hint)* Hysteretic damping β is entered as `K_I = 2β·K_R` (or with the SASSI factor of §6.6).
- Mass: `M_e = M/g` if MOPT `<matrix>` = 1 (weight units, applies to all GENERAL elements). Otherwise `M_e = M`. The 05a screenshot default is mass units.
- *(inferred)* Local axes for 3 nodes: the same construction as BEAMS (e1 along I→J, e2 in plane IJK toward K, e3 = e1 × e2). Transform with `T = blockdiag(Λ,Λ,Λ,Λ)`. The figure shows parallel triads at I and J, consistent with this (§12, Q14).
- Each GENERAL element couples only **two** nodes. A larger super-element is assembled from several GENERAL elements, one per node pair. *(inferred)* The user or the converter must split the diagonal blocks so they are not counted twice.

---

## 5. Element commands (§9.4)

### 5.1 DELL — delete soil layers [IO]
`DELL,<m1>,[<m2>],[<inc>]`. Defaults: m2 = m1, inc = 1.
- Elements that still reference a deleted layer give CHECK Error 19.

### 5.2 DELM — delete materials [IO]
`DELM,<m1>,[<m2>],[<inc>]`. Defaults: m2 = m1, inc = 1. Dangling references give Error 13.

### 5.3 DELR — delete real properties [IO]
`DELR,<r1>,[<r2>],[<inc>]`. Defaults: r2 = r1, inc = 1. Dangling references give Error 26.

### 5.4 DELSC — delete spring properties [IO]
`DELSC,<r1>,[<r2>],[<inc>]` (no underline printed). Defaults: r2 = r1, inc = 1. Dangling references give Error 33.

### 5.5 E — define an element [CORE][IO]
`E,<ne>,<n1>,<n2>,…,<n8>`

| # | Arg | Meaning |
|---|---|---|
| 1 | ne | element number **within the active group** (starts at 1 in each group) |
| 2… | n1…n8 | node numbers. The count depends on the group type (§4.1). |

- Node order per type:

| Group | Node order |
|---|---|
| SOLID | n1..n8 per Figure 9.1, repeating numbers for prism or pyramid shapes |
| BEAMS | I, J, K |
| SHELL / TSHELL / PLANE | I, J, K, L (counter-clockwise). Omit L for a triangle. |
| SPRING | I, J |
| GENERAL | I, J (global input) or I, J, K (local input) |

- **A group must be active** (GROUP).
- The new element gets `mat = mact` and `prop = ract` (MACT, RACT), `etype = 0`, `eint = 0`, `thick = 0` and releases all 0 *(inferred defaults)*.
- An existing ne is overwritten *(inferred)*.
- CHECK:
  - Error 7: too few nodes for the group type.
  - Error 41: an undefined node.
  - Error 8: a 2-node element of zero length.
  - Error 9: collinear or coincident 3- or 4-node elements.
  - Error 11: zero area.
  - Error 12: warped quadrilateral. Warning 3: slightly warped.
  - Warning 2: distorted angles.
- Examples:
  - `GROUP,1,SOLID` / `E,1,1,2,7,6,26,27,32,31`
  - `GROUP,2,SHELL` / `E,1,1,2,7,6` / `E,2,2,3,8` (triangle)
  - `GROUP,3,BEAMS` / `E,1,21,22,99` (K = 99)

### 5.6 ECOMPR — compress element numbering [IO]
`ECOMPR` (abbr. `ECOM`)
- Renumbers the elements of the **active group** as 1..N, keeping their relative order. This removes numbering gaps.
- *(inferred)* Element output requests (EOUT) that refer to elements of this group are remapped, or flagged.

### 5.7 EDEL — delete elements [IO]
`EDEL,<e1>,[<e2>],[<inc>]`. Defaults: e2 = e1, inc = 1. Acts on the active group.
- The gaps this leaves must be compressed before AFWRITE (Error 6).

### 5.8 EGEN — generate elements by copying a pattern [CORE][IO]
`EGEN,[<itim>],<ninc1>,<e1>,[<e2>],[<inc>],[<ee>]`

| # | Arg | Default | Meaning |
|---|---|---|---|
| 1 | itim | 1 | number of new sets, not counting the pattern |
| 2 | ninc1 | — | node-number increment per set |
| 3 | e1 | — | first pattern element |
| 4 | e2 | e1 | last pattern element |
| 5 | inc | 1 | pattern step |
| 6 | ee | max element number + 1 | number of the first generated element |

Algorithm *(inferred numbering and attribute copying)*:
```
pattern = [e for e in range(e1,e2,inc) if e exists in active group]
next = ee
for k in 1..itim:
  for p in pattern:
     new.nodes = [n + k*ninc1 if n != 0 else 0 for n in p.nodes]   # includes the BEAMS/GENERAL K node (see §12 Q15)
     copy p.mat, p.prop, p.etype, p.eint, p.thick, p.ki, p.kj
     element[next] = new; next += 1
```
- Example (a 5 × 4 node grid in X–Z with nodes 1..20, 5 per row):
  - `E,1,1,2,7,6` then `EGEN,3,1,1` creates elements 2..4.
  - `EGEN,2,5,1,4` then creates elements 5..12, two more rows with nodes +5 and +10.

### 5.9 EINT — integration order [CORE][IO]
`EINT,<e1>,<e2>,[<inc>],<order>`. inc defaults to 1. Acts on the active group.

| Group | order | Meaning |
|---|---|---|
| SOLID | 0 | rectangular elements. Default integration is "2×2". |
| | 1 | skewed elements |
| | 2 | extremely distorted elements. Avoid distorted elements for stress analysis wherever possible. |
| TSHELL | 0 | **Reduced** (default): 1×1 for plate bending and transverse shear |
| | 1 | **Selective**: 2×2 for plate bending, 1×1 for transverse shear |

- *(inferred)* The Gauss rule per SOLID order is in §4.4 (§12, Q9). EINT on other group types is ignored with a warning.
- Example: `EINT,1,240,1,1`.

### 5.10 ELIST — list elements [UI]
`ELIST,[<e1>],[<e2>],[<inc>]` (abbr. `ELIS`). With no arguments, all elements of the active group are listed.
- *(inferred)* Columns: ne, nodes, mat, prop, etype, eint, thick, KI/KJ.

### 5.11 ETYPE — structural / excavated / embedded flag [CORE][IO]
`ETYPE,<e1>,<e2>,[<inc>],<type>` (no underline printed). inc defaults to 1. Acts on the active group.

| type | SOLID / PLANE | SHELL / TSHELL |
|---|---|---|
| 0 (default) | excavated soil if **below the ground surface**, otherwise structural | structural |
| 1 | structural (including near-field backfill modelled as structure) | structural |
| 2 | **excavated soil** (MSET index → L table) | **embedded / buried shell** (nodes are interaction nodes, not part of the excavation volume) |

- ETYPE has no effect on BEAMS (xref §9.7.19).
- For ETYPE 0, the classification uses the **ground elevation** (HOUSE `<gelev>` / GROUNDELEV) and is resolved **when AFWRITE writes the `.hou` file**. WRITE keeps the implicit 0 in the `.pre` file (xref §9.7.6).
- The test for "below ground surface" is not specified (xref 05b §3.5). The recommendation there: every node is at or below g_elev + tol, and the centroid is strictly below.
- INTGEN options 1–3 and 5 need ETYPE set **explicitly**; ETYPEGEN does this for all elements.
- Option A (ANSYS interface) ignores elements that are not explicitly structural.
- Example: `GROUP,10,SOLID` / `ETYPE,1,400,1,2`.

### 5.12 GDEL — delete groups [IO]
`GDEL,<g1>,[<g2>],[<inc>]`. Defaults: g2 = g1, inc = 1.
- The groups and all their elements are removed. *(inferred)* If the active group is deleted, no group is active.

### 5.13 GLIST — list groups [UI]
`GLIST,[<g1>],[<g2>],[<inc>]` (abbr. `GLIS`). With no arguments, all groups are listed.
- *(inferred)* Columns: ng, type (number and string), title, element count, and whether the group is active.

### 5.14 GROUP — create or activate a group [CORE][IO]
`GROUP,<ng>,<type>` (abbr. `GROU`)
- If group ng does not exist, it is created with the given type and becomes active.
- If it exists, it becomes active.
- `<type>` is a number (1, 2, 3, 4, 5, 7, 9) or a string (SOLID, BEAMS, SHELL, PLANE, TSHELL, SPRING, GENERAL), case-insensitive.
- *(inferred)*
  - `<type>` may be omitted when activating an existing group.
  - If a different type is given for an existing group, raise an error and point to MTYPE.
  - Also accept `BEAM` as an alias of `BEAMS`.
- Example: `GROUP,1,SOLID`, `GROUP,2,5` (TSHELL), `GROUP,1` (re-activate).

### 5.15 GTIT — group title [UI][IO]
`GTIT,[gr],<title>`. gr defaults to the active group.
- *(inferred)* The title is the rest of the line after the second comma, with commas kept.
- Example: `GTIT,3,Reactor building walls`.

### 5.16 KI — beam end-release code at node I [CORE][IO]
`KI,<e1>,[<e2>],[<inc>],<k1>,<k2>,…,<k6>`. Defaults: e2 = e1, inc = 1.
- **The active group must be of type BEAMS** (error otherwise).
- k1..k6 ∈ {0, 1} correspond to **P1, P2, P3, M1, M2, M3** at node I, in local axes 1/2/3.
- **1 = the end force is zero (released)**: a hinge or roller.
- *(inferred)* Values other than 0 or 1 are an error.
- Example: `KI,1,10,1,0,0,0,0,1,1` releases the two bending moments at I (a pinned I-end).

### 5.17 KJ — beam end-release code at node J [CORE][IO]
`KJ,<e1>,[<e2>],[<inc>],<k1>,…,<k6>`: the same as KI, for node J.
- CHECK **Error 10** if any component is released at both I and J.

### 5.18 L — define a soil layer [CORE][IO]
`L,<nm>,<thick>,<weight>,<pveloc>,<sveloc>,<pdamp>,<sdamp>`

| # | Arg | Meaning | CHECK |
|---|---|---|---|
| 1 | nm | layer number | |
| 2 | thick | layer thickness | Error 20 if ≤ 0 |
| 3 | weight | specific weight γ | Error 21 if < 0 |
| 4 | pveloc | P-wave velocity Vp | Error 22 if < 0 |
| 5 | sveloc | S-wave velocity Vs | Error 23 if < 0 |
| 6 | pdamp | P-wave damping ratio βp (fraction, e.g. 0.02) | Error 24 if < 0 |
| 7 | sdamp | S-wave damping ratio βs | Error 25 if < 0 |

- Uses:
  1. Properties of **excavated-soil** SOLID or PLANE elements, through MSET (§6.2).
  2. Soil layering and halfspace properties for the **SITE** module.
- For embedded models, the layers must be defined **in the same order as in the TOPL command**, and embedment layers must **not be repeated** in TOPL (xref §9.2.44).
- *(inferred)* Damping ratios are fractions, not percent. Warn if a value is > 0.5.
- Example: `L,1,10.0,0.120,2000.0,1000.0,0.02,0.02`.

### 5.19 LLIST — list soil layers [UI]
`LLIST,<m1>,[<m2>],[<step>]` (abbr. `LLIS`). Defaults: m2 = last defined, step = 1. *(inferred)* No arguments lists all layers.

### 5.20 M — define a material [CORE][IO]
`M,<nm>,<val1>,<val2>,<weight>,<pdamp>,<sdamp>,<type>`

| # | Arg | type = 1 | type = 2 | type = 3 |
|---|---|---|---|---|
| 1 | nm | material number | | |
| 2 | val1 | Young's modulus E | **constrained** modulus M = λ + 2G | P-wave velocity Vp |
| 3 | val2 | Poisson's ratio ν | shear modulus G | S-wave velocity Vs |
| 4 | weight | specific weight γ (all types) | | |
| 5 | pdamp | P-wave damping ratio βp | | |
| 6 | sdamp | S-wave damping ratio βs | | |
| 7 | type | 1, 2 or 3 | | |

- *(inferred)* If `<type>` is blank, use 1.
- **For SHELL elements, pdamp must equal sdamp.**
- CHECK:
  - Error 14: E ≤ 0.
  - Error 15: ν ≤ 0. Note that a **zero or negative ν is rejected**; derived values are checked too.
  - Error 16: γ < 0.
  - Errors 17 and 18: negative damping.
  - *(inferred)* Also flag ν ≥ 0.5, or Vp ≤ √2·Vs, since these give ν ≤ 0.
- The conversions are in §6.1.
- Example (concrete, kip–ft): `M,1,519120,0.17,0.150,0.04,0.04,1`.

### 5.21 MACT — active material / soil-layer index [IO]
`MACT,<index>`
- Every element defined after this command, until the next MACT, gets mat = index.
- *(inferred)* The initial value is 1. The index is applied to all groups.
- Example: `MACT,3`.

### 5.22 MLIST — list materials [UI]
`MLIST,<m1>,[<m2>],[<step>]` (abbr. `MLIS`). Defaults: m2 = last defined, step = 1.
- *(inferred)* List the raw values and the derived E, ν, G, M, Vp, Vs, ρ once g is known.

### 5.23 MSET — assign a material / soil-layer index [CORE][IO]
`MSET,<e1>,[<e2>],[<inc>],<index>`. Defaults: e2 = e1, inc = 1. Acts on the active group.
- How the index is read:
  - In a SOLID or PLANE group, if the element is excavated soil (by default or explicitly), the index is a **soil-layer (L)** number.
  - Otherwise it is a **material (M)** number.
- *(inferred)*
  - Store the raw index and interpret it at CHECK/AFWRITE time, so that a later ETYPE change reinterprets it correctly.
  - MSET on SPRING or GENERAL groups is ignored with a warning.
- Example: `GROUP,11,SOLID` / `MSET,1,120,1,4` (layer 4 for excavation group 11).

### 5.24 MTYPE — change a group's type [IO]
`MTYPE,[<gr>],<type>` (abbr. `MTYP`). gr defaults to the active group. type is a number or a string (as in GROUP).
- *(inferred)* Warn if existing elements do not have the node count the new type needs (CHECK Error 7). Data that the new type does not use is kept but ignored.

### 5.25 MXDEL — delete matrix properties [IO]
`MXDEL,<p1>,[<p2>],[<step>]` (abbr. `MXDE`). Defaults: p2 = p1, step = 1.
- The manual's description uses `<m1>`, `<m2>`, `<inc>` for these same arguments.
- Dangling references give Error 83.

### 5.26 MXI — imaginary stiffness terms of a matrix property [CORE][IO]
`MXI,<p>,<row>,<t1>,<t2>,…<t12>`
- Sets row `<row>` (1..12) of the **imaginary part of the stiffness matrix** of matrix property p.
- The property is **created** if it does not exist.
- Only the **upper triangle** is entered:
  - row r takes `13 − r` terms, from the diagonal to the right;
  - `t1 = A[r,r]`, `t2 = A[r,r+1]`, …, `t_{13−r} = A[r,12]`;
  - row 1 uses t1..t12, row 2 uses t1..t11, and row 12 uses only t1.
- The lower triangle is filled by symmetry.
- *(inferred)* Missing trailing terms are set to 0. Terms beyond `13 − r` are an error.
- Input coordinates: **global** if the element has only nodes I and J; **local** (I, J, K) if K is given.
- Example (two-node axial dashpot-like term on X, β = 0.02, k = 1000): `MXI,1,1,40,0,0,0,0,0,-40` and `MXI,1,7,40`.

### 5.27 MXLIST — list a matrix property [UI]
`MXLIST,<p>` (abbr. `MXLI`). Lists the real stiffness, imaginary stiffness and mass/weight matrices of property p.

### 5.28 MXM — mass/weight terms of a matrix property [CORE][IO]
`MXM,<p>,<row>,<t1>,…<t12>`
- The same row and upper-triangle rules as MXI, for the **mass / weight matrix**.
- Enter the terms in **mass or weight units, as chosen by MOPT `<matrix>`** (Options / Model). The choice applies to all GENERAL elements.
- The global/local rule is the same as for MXI.

### 5.29 MXR — real stiffness terms of a matrix property [CORE][IO]
`MXR,<p>,<row>,<t1>,…<t12>`
- The same rules, for the **real part of the stiffness matrix**.
- Example (a k = 1000 axial X-spring between I and J):
  - `MXR,1,1,1000,0,0,0,0,0,-1000` (row 1: K11 = 1000, K17 = −1000)
  - `MXR,1,7,1000` (row 7: K77)

### 5.30 R — define a real (beam section) property [CORE][IO]
`R,<nm>,<axial>,<shear2>,<shear3>,<tors>,<flex2>,<flex3>`

| # | Arg | Symbol | Meaning | CHECK |
|---|---|---|---|---|
| 1 | nm | | property number | |
| 2 | axial | A | axial area | Error 27 if ≤ 0 |
| 3 | shear2 | As2 | shear area for local axis 2 | Error 28 if < 0 |
| 4 | shear3 | As3 | shear area for local axis 3 | Error 29 if < 0 |
| 5 | tors | J | torsional inertia moment | Error 30 if ≤ 0 |
| 6 | flex2 | I2 | flexural inertia about local axis 2 | Error 31 if < 0 |
| 7 | flex3 | I3 | flexural inertia about local axis 3 | Error 32 if < 0 |

- Set shear2 = shear3 = 0 to exclude shear deformation.
- Section formulas from the manual's figure are in §6.3.
- Example (rectangle b = 0.5 along axis 3, h = 1.0 along axis 2): `R,1,0.5,0.41667,0.41667,0.02861,0.010417,0.041667`.

### 5.31 RACT — active real / spring / matrix property index [IO]
`RACT,<index>`
- Every element defined after this command, until the next RACT, gets prop = index.
- *(inferred)* The initial value is 1.

### 5.32 RLIST — list real properties [UI]
`RLIST,<r1>,[<r2>],[<step>]` (abbr. `RLIS`). Defaults: r2 = last defined, step = 1.

### 5.33 RSET — assign a real / spring / matrix property index [CORE][IO]
`RSET,<e1>,[<e2>],[<inc>],<index>`. Defaults: e2 = e1, inc = 1. Acts on the active group.
- How the index is read:
  - BEAMS → real-property table (R).
  - SPRING → spring-property table (SC).
  - GENERAL → matrix-property table (MXR/MXI/MXM).
- *(inferred)* Ignored, with a warning, for other group types.

### 5.34 SC — define a spring property [CORE][IO]
`SC,<nm>,<scx>,<scy>,<scz>,<scxx>,<scyy>,<sczz>,<damp>`

| # | Arg | Meaning | CHECK |
|---|---|---|---|
| 1 | nm | spring property number | |
| 2–4 | scx, scy, scz | translational constants (force/length) in global X, Y, Z | Errors 34–36 if < 0 |
| 5–7 | scxx, scyy, sczz | rotational constants (moment/radian) about global X, Y, Z | Errors 37–39 if < 0 |
| 8 | damp | damping ratio of the spring (fraction) | *(inferred)* error if < 0 |

- The constants are added directly to the global stiffness. They are **global and uncoupled**.
- Example: `SC,1,1.0E5,1.0E5,2.0E5,0,0,0,0.05`.

### 5.35 SCLIST — list spring properties [UI]
`SCLIST,<r1>,[<r2>],[<step>]` (abbr. `SCLI`). Defaults: r2 = last defined, step = 1.

### 5.36 THICK — shell thickness [CORE][IO]
`THICK,<e1>,<e2>,[<inc>],<thick>` (abbr. `THIC`). inc defaults to 1.
- Sets the thickness of elements e1..e2 in the active **SHELL** group. *(inferred)* TSHELL groups too.
- *(inferred)* Error if the active group is not SHELL/TSHELL, or if thick ≤ 0.
- Example: `THICK,1,12,1,1.5`.

---

## 6. Property tables and conversions [CORE]

### 6.1 Material (M) conversions *(inferred standard elasticity; definitions from the manual)*
Let ρ = γ/g.

| type | Given | Derived |
|---|---|---|
| 1 | E, ν | G = E/(2(1+ν)); M = E(1−ν)/((1+ν)(1−2ν)); Vs = √(G/ρ); Vp = √(M/ρ) |
| 2 | M, G | ν = (M−2G)/(2(M−G)); E = 2G(1+ν) = G(3M−4G)/(M−G) |
| 3 | Vp, Vs | G = ρ·Vs²; M = ρ·Vp²; ν = (Vp²−2Vs²)/(2(Vp²−Vs²)); E = 2G(1+ν) |

- Lamé constant: λ = M − 2G.
- Type 3 needs g, so convert lazily (§1.6).
- βp goes with the constrained/P-wave (dilatational) modulus M. βs goes with the shear modulus G.

### 6.2 Soil layer (L) used by excavated-soil elements
- ρ = γ/g, G = ρ·Vs², M = ρ·Vp², with damping βp and βs.
- *(xref 02_theory)* The excavated-soil element properties should match the free-field layer at the same depth. That is why MSET points to the same L table the SITE module uses.
  - Use the strain-compatible values from SOIL if the profile was iterated (the user's responsibility).
  - Thickness is not used by element formulation; it is used by SITE.
- *(xref line 5644)* WARNING: always check, in the HOUSE output file, that excavation SOLIDs got the right soil-layer numbers.

### 6.3 Beam section formulas *(figure, §9.4.30)*
The shear correction (form) factor f gives the shear area `As = A/f`.

| Section | I2 | I3 | f2 = f3 | A | J |
|---|---|---|---|---|---|
| Solid circle, radius r | πr⁴/4 | πr⁴/4 | 10/9 | πr² | πr⁴/2 |
| Solid rectangle: b along axis 3, h along axis 2 | h·b³/12 | b·h³/12 | 6/5 | b·h | [1/3 − 0.21·(b/h)·(1 − b⁴/(12h⁴))]·h·b³ |

- The rectangle torsion formula is printed as `0.21(bh)`. That is a typo for `0.21(b/h)` (Roark), valid for h ≥ b. If b > h, swap b and h in the J formula only (§12, Q16).

### 6.4 Spring (SC)
- Six uncoupled global constants and one damping ratio. See §4.8.

### 6.5 Matrix property (MXR / MXI / MXM)
- Three 12×12 symmetric matrices. See §4.9 and §5.26–5.29.
- Unset terms are 0.

### 6.6 Complex-modulus damping factor *(inferred; xref 01_capabilities §5, Q1)*
- Use the same factor for all materials, layers, springs and structural elements:
  - default (SASSI/SHAKE form): `c(β) = 1 − 2β² + 2iβ√(1−β²)`;
  - selectable alternative: `c(β) = 1 + 2iβ`.
- Applying the factor:
  - Solid and plane elements: `M* = M·c(βp)`, `G* = G·c(βs)`, `λ* = M* − 2G*`.
  - Beams: complex E\* and G\* derived from M\* and G\* using the type-2 formulas, so that E\* = E·c(β) when βp = βs.
  - Shells: βp = βs, so E\* = E·c(β) and ν stays real.
  - Springs: `k* = k·c(damp)`.
  - GENERAL: `K_R + i·K_I`, as entered.
- Mass is always real.

---

## 7. Load commands (§9.5)

### 7.1 Common semantics [CORE]
- **Forces and moments** (F, MM) are used **only by the FORCE module**, for external-force or vibration analysis (xref 05b §FORCE, line 6315).
  - All loads share **one reference time history** f_ref(t), which has the maximum reference amplitude and **zero arrival time**.
  - Each DOF load has a **factor** a (its maximum amplitude divided by the reference maximum) and an **arrival time** t0 (a time lag, which allows moving loads).
  - *(xref 05c)* With the e^{+iωt} convention: `P_dof(ω) = a·F_ref(ω)·exp(−iω·t0)`.
  - A factor of 0 means no load on that DOF.
  - Seismic and external-force analyses cannot be combined in one run.
- **Masses** (MT, MR) are added to the structure's mass matrix on the node DOFs.
  - Translational masses act on UX/UY/UZ. Rotational masses act on ROTX/ROTY/ROTZ.
  - Masses and loads on **fixed** DOFs are ignored (Warnings 5, 6, 10, 11).
- **Units** of nodal masses are set **per node** by MUNITS: 0 = mass, 1 = weight. Weight is converted to mass by dividing by g.
  - Rotational "weight" units are weight × length².
  - Masses on GENERAL elements (MXM) use MOPT `<matrix>` instead.
  - *(xref 5.1, line 3407)* Recommended workflow: define nodal masses with MT/MR, then set their units with MUNITS.
- **Overwrite flags** (MOPT, xref 05a):
  - `<mass>`: 0 = add to an existing nodal mass, 1 = set (replace).
  - `<force>`: the same for forces and moments.
  - The 05a screenshot defaults are set (1) for both.
  - *(inferred)* The flags act when MT/MR/F/MM run. In add mode for F/MM, factors are added and the arrival time of the new command replaces the old one, with a warning if they differ (§12, Q17).
- CHECK:
  - Error 40: a mass on a node beyond the defined node range.
  - Error 61: FORCE requested but no F or MM defined.
  - Error 62: a force or moment on an undefined node.

### 7.2 F — define a nodal force [CORE][IO]
`F,<n>,<fx>,<fy>,<fz>,<tx>,<ty>,<tz>`

| # | Arg | Meaning |
|---|---|---|
| 1 | n | node |
| 2–4 | fx, fy, fz | force factors in global X, Y, Z |
| 5–7 | tx, ty, tz | arrival times (s) for the X, Y, Z components |

- Defaults *(inferred)*: 0.
- Example: `F,50,1.0,0,0` (the reference load in X at node 50), then `F,51,0.5,0,0,0.1` (half the amplitude, arriving 0.1 s later).

### 7.3 FDEL — delete forces [IO]
`FDEL,<n1>,[<n2>],[<inc>]`. Defaults: n2 = n1, inc = 1.

### 7.4 FLIST — list forces [UI]
`FLIST,[<n1>],[<n2>],[<inc>]` (abbr. `FLIS`). inc defaults to 1. With no arguments, all forces are listed.

### 7.5 FSCALE — scale force factors [CORE][IO]
`FSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` (abbr. `FSCA`)
- Multiplies the factors fx, fy, fz of the forces on nodes n1..n2 (step inc) by sx, sy, sz.
- **A factor of 0.0 is replaced by 1.0.**
- If n1 and n2 are omitted, they are set from the **last two defined nodal forces**. *(inferred)* The range runs between the node numbers of the two most recent F entries, low to high.
- inc defaults to 1.
- Arrival times are unchanged *(inferred)*.
- It works like NSCALE.

### 7.6 MM — define a nodal moment [CORE][IO]
`MM,<n>,<fxx>,<fyy>,<fzz>,<txx>,<tyy>,<tzz>`
- Moment factors about global X, Y, Z, and their arrival times. The semantics are the same as F.

### 7.7 MMDEL — delete moments [IO]
`MMDEL,<n1>,[<n2>],[<inc>]` (abbr. `MMDE`). Defaults: n2 = n1, inc = 1.

### 7.8 MMLIST — list moments [UI]
`MMLIST,[<n1>],[<n2>],[<inc>]` (abbr. `MMLI`). With no arguments, all moments are listed.

### 7.9 MR — define rotational masses [CORE][IO]
`MR,<n>,<mxx>,<myy>,<mzz>`
- Rotational masses (mass moments of inertia) about global X, Y, Z at node n.
- The manual says these are "in weight units". **Conflict:** MUNITS also sets the units of rotational masses per node (§12, Q18). *(inferred)* Use the node's MUNITS flag.

### 7.10 MRGEN — generate rotational masses [CORE][IO]
`MRGEN,[<itim>],[<ninc>],<n1>,<n2>,[<inc>],[<mxx>],[<myy>],[<mzz>]` (abbr. `MRGE`)

| # | Arg | Default | Meaning |
|---|---|---|---|
| 1 | itim | 1 | number of new sets |
| 2 | ninc | n2 − n1 + 1 | node increment per set |
| 3–4 | n1, n2 | — | pattern node range |
| 5 | inc | 1 | pattern step |
| 6–8 | mxx, myy, mzz | 0 | **increment** added to the pattern masses for each set |

```
for k in 1..itim:
  for n in range(n1,n2,inc) if rmass[n] exists:
     rmass[n + k*ninc] = rmass[n] + k*(mxx,myy,mzz)    # subject to MOPT overwrite flag
     mass_units[n + k*ninc] = mass_units[n]              # (inferred)
```
- It works like NGEN. *(inferred)* Pattern nodes without a rotational mass are skipped.

### 7.11 MRDEL — delete rotational masses [IO]
`MRDEL,<n1>,[<n2>],[<inc>]` (abbr. `MRDE`). Defaults: n2 = n1, inc = 1.

### 7.12 MRSCALE — scale rotational masses [IO]
`MRSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` (abbr. `MRSC`)
- Multiplies mxx, myy, mzz by sx, sy, sz.
- A factor of 0.0 becomes 1.0.
- If n1 and n2 are omitted, they are set from the last two defined rotational masses. inc defaults to 1.

### 7.13 MSCALE — scale moment factors [IO]
`MSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` (abbr. `MSCA`)
- The same as FSCALE, for moments. If n1 and n2 are omitted, they are set from the last two defined moments.

### 7.14 MT — define translational masses [CORE][IO]
`MT,<n>,<mx>,<my>,<mz>`
- Translational masses in global X, Y, Z at node n. Units follow the node's MUNITS flag.
- Different values per direction are allowed.
- Example: `MT,100,2.5,2.5,2.5`.

### 7.15 MTDEL — delete translational masses [IO]
`MTDEL,<n1>,[<n2>],[<inc>]` (abbr. `MTDE`). Defaults: n2 = n1, inc = 1.

### 7.16 MTGEN — generate translational masses [CORE][IO]
`MTGEN,[<itim>],[<ninc>],<n1>,<n2>,[<inc>],[<mx>],[<my>],[<mz>]` (abbr. `MTGE`)
- The same algorithm as MRGEN, for translational masses.
- Example: `MTGEN,3,10,100,105` copies the masses on nodes 100..105 to 110..115, 120..125 and 130..135.

### 7.17 MTLIST — list masses [UI]
`MTLIST,[<n1>],[<n2>],[<inc>]` (abbr. `MTLI`)
- Lists **both translational and rotational** masses. With no arguments, all are listed.
- *(inferred)* Also show the units flag (mass or weight) and the converted mass.

### 7.18 MTSCALE — scale translational masses [IO]
`MTSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` (abbr. `MTSC`)
- The same as MRSCALE, for translational masses. If n1 and n2 are omitted, they are set from the last two defined translational masses.

### 7.19 MUNITS — units of the nodal masses [CORE][IO]
`MUNITS,<n1>,[<n2>],[<step>],<units>` (abbr. `MUNI`). Defaults: n2 = n1, step = 1.
- units: **0 = mass units**, **1 = weight units**. The flag applies to both the translational and rotational masses of each node in the range.
- *(inferred)* The default for nodes never given MUNITS is unknown (§12, Q18). Recommended: 1 (weight units), which matches the MR wording. Always show the flag in MTLIST and in the HOUSE echo.
- Example: `MUNITS,100,135,1,0` (mass units for nodes 100..135).

---

## 8. CHECK rules tied to this section [CORE][IO] *(xref Chapter 10)*

| Code | Condition | Source command |
|---|---|---|
| Err 6 | element numbering gap in a group | E, EDEL (fix with ECOMPR) |
| Err 7 | element has fewer nodes than its group type needs | E, MTYPE |
| Err 8 | 2-node element with coincident nodes | E (BEAMS I = J, SPRING, GENERAL) |
| Err 9 | 3- or 4-node element with coincident or collinear nodes | E (includes beam I-J-K) |
| Err 10 | the same release component is 1 at both I and J | KI, KJ |
| Err 11 | zero area (4-node) | E |
| Err 12 | warped 4-node element (error level) | E |
| Err 13 | material not defined | MSET / MACT, M |
| Err 14–18 | E ≤ 0; ν ≤ 0; γ < 0; βp < 0; βs < 0 | M |
| Err 19 | soil layer not defined | MSET on excavated elements, L |
| Err 20–25 | thick ≤ 0; γ < 0; Vp < 0; Vs < 0; βp < 0; βs < 0 | L |
| Err 26 | real property not defined (BEAMS) | RSET, R |
| Err 27–32 | A ≤ 0; As2 < 0; As3 < 0; J ≤ 0; I2 < 0; I3 < 0 | R |
| Err 33 | spring property not defined | RSET, SC |
| Err 34–39 | a spring constant < 0 | SC |
| Err 40 | mass on a node beyond the defined node range | MT, MR, MTGEN, MRGEN |
| Err 41 | element uses an undefined node | E, NDEL |
| Err 42 / 43 | no nodes / no groups (fatal) | N / GROUP |
| Err 61 / 62 | no F/MM for FORCE; F/MM on an undefined node | F, MM |
| Err 83 | matrix property not defined | RSET, MXR/MXI/MXM |
| Err 124 | node fully fixed **and** interaction | D, INT |
| Warn 1 | gap in node numbering (AFWRITE writes the missing nodes, fixed) | N |
| Warn 2 | distorted element (small angle relative to the largest angle); all faces checked for SOLID | E |
| Warn 3 | slightly warped 4-node SHELL/PLANE | E |
| Warn 4 | unused node, including nodes used only as BEAMS K-nodes (AFWRITE fixes it) | N, E |
| Warn 5 / 6 | translational / rotational mass on a fixed DOF (ignored) | MT/MR + D |
| Warn 7 | group with no elements (skipped) | GROUP |
| Warn 10 / 11 | force / moment on a fixed DOF (ignored) | F/MM + D |

- *(inferred)* The numeric thresholds for Warnings 2/3 and Error 12 are not given. Proposals:
  - Warning 2: min angle / max angle < 0.2 (or an interior angle < 15° or > 165°).
  - Warning 3: warp = out-of-plane distance / mean diagonal > 1e−3.
  - Error 12: warp > 0.05.
  - See §12, Q19.
- AFWRITE runs CHECK first. A module that fails CHECK does not get its input file written.

---

## 9. ANSYS cross-walk [UI] (informational, for users coming from ANSYS)

| ACS SASSI | Nearest ANSYS | Key differences |
|---|---|---|
| `N,nd,x,y,z` | `N,NODE,X,Y,Z` | Same. Z must be vertical up. |
| `NGEN,itim,step,n1,n2,inc,dx,dy,dz` | `NGEN,ITIME,INC,NODE1,NODE2,NINC,DX,DY,DZ` | Same argument order. The defaults use the last two defined nodes. |
| `FILL,n1,n2,nr` | `FILL,NODE1,NODE2,NFILL` | |
| `CSYS,ns` | `CSYS,KCN` | Only Cartesian local systems. |
| `LOC,ns,type,x0,y0,z0,txy,tyz,txz` | `LOCAL,KCN,KCS,XC,YC,ZC,THXY,THYZ,THZX` | **The names are swapped:** SASSI `LOCAL` is node-based, like ANSYS `CS`. |
| `D,n1,n2,inc,val,labels` | `D,NODE,Lab,VALUE,,NEND,NINC` | val is a **0/1 fixity code**, not a displacement value. The range comes first. |
| `GROUP,ng,SOLID` + `E` | `ET` + `TYPE` + `E` | Element numbers restart at 1 in each group. Elements are grouped by type. |
| SOLID / BEAMS / SHELL / TSHELL / PLANE / SPRING / GENERAL | SOLID45/185; BEAM4/188; SHELL63; SHELL181; PLANE42 (plane strain); COMBIN14; MATRIX27/50 | Damping comes from complex moduli per material, not from Rayleigh α/β. |
| `M,nm,E,ν,γ,βp,βs,1` | `MP,EX/PRXY/DENS` + `MP,DMPR` | Specific **weight**, not density. Two damping ratios (P and S). |
| `R,nm,A,As2,As3,J,I2,I3` | `SECTYPE/SECDATA` or BEAM4 real constants | Shear **areas**, not shear coefficients. |
| `MT`, `MR` + `MUNITS` | MASS21 | Mass or weight units are set per node. |
| `F,n,fx,fy,fz,tx,ty,tz` | `F` + `TABLE` | The values are factors on a single reference time history, with arrival times. |

- *(xref line 4200)* The UI `.cdb` converter maps BEAM4/44 (with a K node, RBLOCK), COMBIN14, MASS21, SOLID45/185, SHELL63/181 and BEAM188/PIPE288 onto these types.

---

## 10. Worked input example (illustrative, not from the manual) [IO]
A surface basemat on soil: 4 × 4 SHELL elements in the X–Y plane at z = 0, and a 5 × 5 ft concrete
column from the mat centre (node 13) to a lumped top mass at node 100. Units are kip–ft–s.
```
GRAVITY,32.2
N,1,0,0,0
N,5,40,0,0
FILL                           (nodes 2..4)
NGEN,4,5,1,5,1,0,10,0          (nodes 6..25 at y = 10..40)
INT,1,25,1,1,0                 (all mat nodes are interaction nodes: surface foundation)
D,1,12,1,1,ROTZ                (drilling DOF of the flat mat; node 13 keeps ROTZ for column torsion)
D,14,25,1,1,ROTZ
M,1,519120,0.17,0.150,0.04,0.04,1
GROUP,1,SHELL
GTIT,1,Basemat
MACT,1
E,1,1,2,7,6
EGEN,3,1,1                     (elements 2..4)
EGEN,3,5,1,4                   (elements 5..16)
THICK,1,16,1,2.0
N,100,20,20,30                 (column top)
NMED,101,1,5,21,25             (helper node at the plan centroid, (20,20,0); only to show NMED)
NDEL,101
N,102,-10,20,0                 (beam K reference node: off the column axis, not an interaction node)
D,102,102,1,1,ALL
GROUP,2,BEAMS
R,1,25,20.83,20.83,88.02,52.08,52.08
MACT,1
RACT,1
E,1,13,100,102                 (I = 13, J = 100, K = 102: axis 1 vertical, axis 2 toward -X)
MUNITS,100,100,1,1
MT,100,500,500,500             (500 kip weight at the top)
L,1,20,0.120,2000,1000,0.02,0.02
```
- Text in parentheses is explanation only. It is not part of the syntax; the comment syntax belongs to the general parser spec.
- `D,…,ROTZ` assumes the flat mat lies in the X–Y plane. A sloped shell would need FIXSHLROT instead. Node 13 keeps ROTZ because the column's torsion acts on it.
- The R values are those of a 5 × 5 ft square: A = 25, As = A/1.2, I = 5⁴/12, and J from §6.3.

---

## 11. Verification tests (to implement) [CORE]
1. **LOC/CSYS:** `LOC,1,0,10,0,0,45,0,0`, `CSYS,1`, `N,100,5,0,0`, `GLOBAL,100,100,1` gives (13.5355, 3.5355, 0). `NLIST` shows csys 0 after GLOBAL, and `CSYS` is still 1.
2. **LOCAL vs LOC:** n1 = (0,0,0), n2 = (0,1,0), n3 = (−1,0,0). `LOCAL,2,0,1,2,3` gives the same R as `LOC,3,0,0,0,0,90,0,0`, to 1e−12.
3. **Euler order:** `LOC,4,0,0,0,0,0,90,0` maps local y to global +Z.
4. **FILL / NGEN defaults:** the §3.3 and §3.12 examples give the expected node IDs and coordinates. A FILL with no arguments right after NGEN uses the last two defined nodes.
5. **NSCALE zero rule:** `NSCALE,1,1,,0,2,0` changes only y (×2).
6. **D expansion:** `D,1,1,,1,DISP` fixes UX/UY/UZ only. `D,1,1,,,ALL` frees all six.
7. **EGEN:** the §5.8 grid example gives 12 elements with the expected connectivity. ECOMPR after `EDEL,5` gives elements 1..11 with the order kept.
8. **Material conversions:** E = 30000, ν = 0.25, γ = 0.15, g = 32.2 in types 1, 2 and 3 give the same (G, M, Vs, Vp) to 1e−10.
9. **Beam element:**
   - Cantilever, L = 10, tip load P along axis 2, As2 = 0: tip deflection = PL³/(3EI3).
   - With As2: tip deflection = PL³/(3EI3) + PL/(G·As2).
   - Rotating the beam with different K nodes leaves the response along e2 unchanged.
10. **Release:** a beam along X with I fully fixed, and J free only in UY (all other J DOFs fixed with D), has local axis 2 = +Y, chosen through K.
    - Without releases the lateral stiffness is 12·E·I3/L³.
    - With `KJ,1,1,1,0,0,0,0,0,1` (M3 released at J) it is 3·E·I3/L³.
    - Releasing M3 at both ends triggers CHECK Error 10.
11. **Section table:** r = 1 gives I2 = 0.785398, J = 1.570796. b = 0.5, h = 1 gives J = 0.028610.
12. **SPRING:** `SC,1,1000,0,0,0,0,0,0.05` between a free node (MT = 1, mass units) and a fixed node.
    - Undamped frequency √1000/(2π) = 5.033 Hz.
    - Resonant amplification of the complex SDOF ≈ 1/(2β) = 10, for c = 1 + 2iβ.
13. **GENERAL equivalence:** a GM element with `MXR` rows equal to the spring's 12×12 matrix (global, 2 nodes) gives the same response as test 12, and `MXI = 2β·MXR` gives the same damped response. Using 3 nodes with a rotated I-J line, with the local matrix entered along axis 1, reproduces an axial spring along I→J.
14. **SOLID patch test:** a single distorted 8-node brick under a uniform-strain field gives an exact constant stress, with and without incompatible modes. With incompatible modes, a 1-element-thick cantilever approaches the beam theory value.
15. **Shell/TSHELL:** square plate, simply supported, uniform load. The central deflection converges to the Navier solution. Reduced and Selective EINT settings both converge.
16. **Mass units:**
    - MT = 32.2 with MUNITS = 1 equals MT = 1 with MUNITS = 0.
    - MXM in weight units with MOPT `<matrix>` = 1 gives the same result as mass units.
17. **Load phase:** FORCE with F factor 0.5 and t0 = 0.1 s gives an output time history equal to 0.5·f_ref(t − 0.1) (circular shift).
18. **CHECK:** build one model with each fault in §8 and assert the error/warning codes.

---

## 12. Open questions / ambiguities

1. **Abbreviations:** INTLIST, LMOVE, SLIST, DELSC and ETYPE have no underlined short form in the manual. Should INTL, LMOV, SLIS, DELS and ETYP be accepted? *Decision:* yes, as legacy 4-character names.
2. **"Last defined"** in FILL/NGEN/xSCALE/xLIST defaults: definition order or highest number? *Decision:* definition order for nodes and loads; highest index for the property tables.
3. **Coordinate system for generation commands:** are NGEN dx/dy/dz, LMOVE translations, NMOVE/NSCALE factors and FILL interpolation applied in the active system, in the source node's system, or globally? *Decision:* the active system (ANSYS behaviour).
4. **Automatic DOF elimination:** does HOUSE/AFWRITE remove DOFs that no element defines (for example rotations at SOLID-only nodes), or must the user run FIXROT/FIXSLDROT? *Decision:* auto-eliminate the undefined DOFs. Require the FIX* commands for defined-but-unstiffened DOFs (shell drilling), with a warning.
5. **INT codes:** the manual text is cut off after "internal (<". Assume code 3 = internal and default code 0. The computational meaning of intermediate (1), interface (2) and internal (3) nodes in V3 is not given: perhaps legacy PINT/pile or subtraction-method features. *Decision:* store them, list them and pass them through; they have no effect on the solution.
6. **LMOVE/NMOVE list length:** the syntax says 15, the text says 9. *Decision:* accept up to 15.
7. **LOCAL Y axis:** the text says "Y local from X x Z", which is left-handed. *Decision:* Y = Z × X, so that n3 is in the positive x–y quadrant.
8. **NMOVE zero factors:** NSCALE turns 0 into 1, but NMOVE's rule is not stated (its default is 1). *Decision:* a blank field gives 1. An explicit 0 is used as given, with a warning. An alternative is to apply 0 → 1 as in NSCALE.
9. **SOLID EINT mapping:** what Gauss rule do order 0, 1 and 2 use? *Decision:* 2, 3 and 4 points per direction (SAP IV heritage). The manual says only "default 2×2".
10. **Mass details:** beam torsional inertia (ρJ or ρ(I2+I3)) and rotary inertia; whether shells get lumped rotational inertia; whether excavated-soil elements use the same 50/50 mass scheme. *Decision:* torsional ρ·(I2+I3) (polar inertia), no rotary inertia, no shell rotational lumped mass, and the 50/50 scheme for both structure and excavated soil. Note that the excavated soil must use the **same** formulation as any consistent free-field treatment.
11. **SHELL membrane/bending formulation:** the manual only names Kirchhoff (SHELL) and Mindlin–Reissner (TSHELL). Is there an in-plane drilling formulation, and are the plate and membrane coupled for warped elements? *Decision:* flat facet = plane-stress membrane (Q4 with incompatible modes) + DKQ/DKT for SHELL; MITC4 or the manual's reduced/selective Mindlin Q4 for TSHELL. The automatic TSHELL drilling stiffness is small, e.g. 1e−4 of the minimum membrane diagonal.
12. **PLANE orientation:** "counter-clockwise" is ambiguous in the X–Z plane. *Decision:* accept either orientation and fix the node order internally if the Jacobian is negative.
13. **Grounded springs:** can a SPRING or GENERAL element have J = 0 (to ground)? *Decision:* not allowed. The user connects to a fully fixed node.
14. **GENERAL local axes:** how I, J, K define the local 1-2-3 axes is not written out. *Decision:* the same as BEAMS.
15. **EGEN details:** whether the BEAMS/GENERAL K-node is incremented; whether the material, property, thickness, release and ETYPE attributes are copied from the pattern or taken from the active MACT/RACT. *Decision:* increment K, and copy the attributes from the pattern.
16. **Rectangular J formula:** the figure prints `0.21(bh)`. *Decision:* use b/h with b ≤ h (Roark).
17. **Add-mode overwrite for F/MM:** with MOPT `<force>` = 0, how do two arrival times combine? *Decision:* add the factors and use the latest arrival time, with a warning. A physically exact superposition of two different arrival times on one DOF cannot be stored as one (factor, time) pair.
18. **Nodal mass units default:** MR says "in weight units", but MUNITS sets the units per node. *Decision:* default per-node flag = 1 (weight). Confirm against a reference `.hou` file or the HOUSE manual if one becomes available.
19. **CHECK thresholds** for distortion/warping (Warning 2/3, Error 12) are not given. Proposals are in §8.
20. **Redefinition semantics:** redefining N, E, M, L, R or SC on an existing number overwrites it. Should LOC/LOCAL redefinition move the nodes stored in that system? *Decision:* overwrite, and warn; nodes keep their local coordinates.
21. **Deletion cascades:** NDEL (loads and masses on the node), GDEL (output requests), DELM/DELL (dangling indices). *Decision:* delete the node-attached loads and masses; leave dangling indices for CHECK to report.
22. **"Below ground surface" test** for ETYPE 0 (centroid vs all nodes) and the tolerance. *Decision:* follow 05b §3.5.
23. **MACT/RACT initial values** and their scope (global vs per group). *Decision:* both start at 1 and are global.
24. **Thickness for PLANE:** unit thickness (plane strain per unit length) is assumed. There is no command to set it.
25. **TSHELL with THICK:** the manual says THICK sets thickness "from the SHELL active group". *Decision:* it also applies to TSHELL groups.
