# 09 — Command Reference: UI General, Model Checking, Model Generation and Combination, Cut and Submodeling, File Conversion

**Source:** ACS SASSI V3 User Manual, §9.6–§9.11 (extracted-text lines 10688–11783; printed pages
259–282; PDF pages 261–284). I viewed every PDF page in the range as an image. Those pages are all
text with simple tables: there are no dialog screenshots, figures or equations in §9.6–§9.11. The
image check corrected one extraction defect: the MERGESOIL `<Stiff>` default is **10^7**, not "107".

Where the commands depend on facts stated elsewhere in the manual, I looked at those passages and
cite them as *(cross-ref §x / line n)*:
- §1.5.1 (lines 723–963): modelling rules, FIXROT, the FV/FI/FFV methods, EXCSTRCHK
- §2.1 (lines 1477–1590): FV substructuring
- §3 and §4 (lines 2040–2290, 2960–3080): TFU/TFI files, frame files, interaction-node rules
- §5.7–§5.9 (lines 3657–3890): MERGE and MERGESOIL examples, section-cut example, variables
- §6.1.5–§6.1.9 (lines 4136–4345): converter dialogs, which I viewed as images (pp. 91–92)
- §6.5.3 (Check options), the HOUSE options (radius figure, p. 128), §9.2–§9.4 (D, INT, ETYPE, HOUSE, SC, MX*)

**Tags:** **[CORE]** needed for computational correctness · **[IO]** file or input format ·
**[UI]** interface, plotting or convenience · **[ADV]** advanced (Option A/AA, NON, incoherency,
fast solver).

**Provenance markers:**
- *(manual)*: stated explicitly in the manual.
- *(inferred)*: my reconstruction where the manual gives only intent. The implementer must confirm it
  or pick an alternative. Every inferred algorithm is listed again in §11 *Open questions*.

---

## 0. Scope: commands covered

| § | Group | Commands |
|---|---|---|
| 9.7 | UI general | ACTM, AFWRBAT, CPMODEL, CRITFREQ, DMODEL, ETYPEGEN, FIXSLDROT, FIXSHLROT, FIXSPRROT, FIXROT, FRAMECOMBIN, FRAMESEL, GCOM, GETENV, GLB2LOC, GRAVITY, GROUNDELEV, GROUPMAT, INTGEN, MDL, MDLNAME, NCOM, MODELLIST, MODFRAMES, RADIUS, RMVUNUSED, SETENV, TPSD |
| 9.8 | Model checking | EXCSTRCHK, FIXEDINT, FREESPRING, HINGED, INTCOUNT, KINT, USED |
| 9.9 | Generation and combination | EXCAV, MERGE, MERGEGROUP, MERGEPANEL, MERGESOIL, ROTATE, SOILMESH, TRANSLATE, WELD |
| 9.10 | Cut and submodeling | CSECT, CUT2SUB, CUTADD, CUTCLR, CUTRMV, CUTVOL, EXTRACTEXCAV, SLICE, SPLITGROUP, TRANELEM, TRANVOL |
| 9.11 | File conversion | ANSYS, ANSYSMODELTYPE, ANSYSREFORMAT, CONVERT, GENMATRIXDAMP |

TRANELEM and TRANVOL (§9.10.10–9.10.11) fall inside the assigned line range, so they are covered here
as well.

---

## 1. Global conventions these commands depend on

### 1.1 Command syntax [IO][UI]
- A command is one line of **comma-separated fields**: `NAME,arg1,arg2,...`. *(manual)*
- In the manual, `<x>` marks a required argument and `[x]` an optional one.
- An **empty field** (nothing, or only white space, between commas) means "use the default".
  Example: `cutvol,3,,,,,2.53` (the manual writes this as `<blank>`). *(manual, line 3790–3800)*
- Command names are **case-insensitive**. The manual's examples use `Cut2Sub`, `CutAdd`,
  `MergeSoil`, `etypegen` and `Actm`. *(manual)*
- **Full names required:** every command in §9.6 onward is valid **only in the ACS SASSI UI** (not
  the legacy PREP) and must be typed with its **full name**. *(manual, §9.6)*
  - Legacy PREP commands (§9.2–§9.5) accept abbreviations (the underlined prefix, e.g. `MTSC` for
    `MTSCALE`). The UI commands do not.
  - Implementation: register the §9.6+ commands with `allow_abbrev=False`.
- Variables (§5.9, cross-ref): `VAR,name,v1,v2,...` defines a **list of strings plus an integer
  counter**.
  - Reference syntax: `@name`, `@name[i]` (1-based), `@name++`, `@name+5`, `@name=-15`.
  - Variable names are **not case-sensitive**.
  - CRITFREQ and FRAMESEL **write their results into such a variable**. Writing redefines the
    variable, so its counter resets to 0, as `VAR` does.

### 1.2 Model registry and the active model [CORE][UI]
- The UI holds **several models in memory**, each keyed by an **integer reference number**. *(manual)*
- Exactly one model is the **active model**. Commands operate on it unless they take explicit
  model-number arguments.
- A model holds:
  - nodes: number → (x, y, z), coordinate system, BC flags for 6 DOF, interaction code
  - local coordinate systems
  - groups: number → type, title, and element table (element number → node list, MSET, RSET, ETYPE, THICK, KI/KJ)
  - materials, soil layers, real properties (R), spring properties (SC), matrix properties (MXR, MXI, MXM)
  - masses, forces
  - analysis settings, including gravity and ground elevation from `HOUSE`
  - model name, path and title (TIT)
- Element numbers start at 1 in each group. Gaps are not allowed when files are written (ECOMPR
  compresses them). *(cross-ref §9.4.14)*
- **Cuts** (§9.10) are numbered element selections: sets of (group number, element number).
  - The Cut Plot takes a cut number **and** a model number, and is "the only plot that can plot a
    model that is not the current active model" (cross-ref line 8217).
  - So cuts appear to be stored **globally in UI memory**, not inside a model, and are interpreted
    against whichever model uses them. *(inferred; see Q-C1)*

### 1.3 Element types, DOF per node, ETYPE [CORE]
DOF active per element type (cross-ref §9.3.2, D command):

| Group type (number/string) | Nodes | UX | UY | UZ | ROTX | ROTY | ROTZ |
|---|---|---|---|---|---|---|---|
| 1 SOLID | 8 (7, 6, 5 by repeating nodes) | • | • | • | | | |
| 2 BEAMS | 3 (I, J, K) | • | • | • | • | • | • |
| 3 SHELL | 4 or 3 | • | • | • | • | • | • |
| 5 TSHELL | 4 or 3 | • | • | • | • | • | • |
| 4 PLANE (2D) | 4 or 3 | • | | • | | | |
| 7 SPRING | 2 | • | • | • | • | • | • |
| 9 GENERAL | 2 (global axes) or 3 (local axes) | • | • | • | • | • | • |

- A node shared by several element types gets the **union** of their DOF ("beam governs").
- SOLID node order (Figure 9.1, viewed):
  - nodes 1-2-3-4 form the bottom face, counter-clockwise seen from above;
  - nodes 5-6-7-8 sit directly above 1-2-3-4;
  - prisms, pyramids and tetrahedra are made by repeating nodes (e.g. 7=8; 5=6, 7=8; 5=6=7=8).
- BEAMS node K is a geometric reference point only. *(cross-ref §9.4.5)*

**ETYPE** (cross-ref §9.4.11) is a per-element attribute for SOLID, PLANE, SHELL and TSHELL:

| ETYPE | Meaning |
|---|---|
| 0 (default, *implicit*) | SOLID/PLANE: excavation if below the ground surface, otherwise structure. SHELL: structure. |
| 1 | structure (SOLID, PLANE, SHELL/TSHELL). May include near-field or backfill soil modelled as structure. |
| 2 | SOLID/PLANE: excavated soil. SHELL/TSHELL: **buried, or sitting on the soil**. The nodes of an ETYPE=2 shell are **interaction nodes** and are not part of the excavation volume. ETYPE=2 has no effect on beams for INTGEN. |

- Excavated-soil elements' `MSET` index points to the **soil-layer table** (L command), not to the
  material table. *(cross-ref line 5706, §9.4.18)*
- Ground elevation is set by `HOUSE,<gravity>,<gelev>,...` or by GROUNDELEV, and decides the
  implicit ETYPE=0 classification.
- `AFWRITE` resolves implicit ETYPE=0 when it writes `.hou`.
- `WRITE` keeps implicit ETYPE=0 as-is in `.pre`.
- **INTGEN (options 1–3, 5), EXTRACTEXCAV, MERGESOIL and the soil-layer redefinition commands need
  ETYPE defined explicitly** (1 or 2, not 0).

### 1.4 Interaction-node rules the generators must respect [CORE]
These come from cross-ref §1.5.1 items 5–11 and §4 item 13:
- Interaction nodes belong to the **excavated-soil model**. Structural nodes are interaction nodes
  only where the basement and the excavation **share** nodes on the excavation's lateral and bottom
  surfaces.
- Every interaction node lies **at or below the ground surface** and **on a soil-layer interface**,
  and has translational DOF only.
- Interaction nodes should be **numbered in ascending order**. For incoherent analysis, numbering
  should run from the bottom layer up to the ground surface.
- Interaction nodes **always** include all nodes on the excavation / far-field interface. Skipping
  any of them can badly degrade accuracy.
- Excavation **interior** nodes must not be shared with structure basement nodes. EXCSTRCHK checks
  this.
- Nodes of buried BEAM/SHELL elements outside the excavation volume should be interaction nodes. For
  shells this is automated by ETYPE=2 together with INTGEN.
- `INT,<n1>,<n2>,[<inc>],<set>,[<code>]` sets (`set`=1) or resets (`set`=0) a node set as:
  - interaction (`code`=0)
  - intermediate (`code`=1)
  - interface (`code`=2)
  - internal (code text truncated in the manual; presumably 3)
  
  INTGEN works on the *interaction* code (0). *(cross-ref §9.3.5)*

### 1.5 Geometric tolerance [CORE] *(inferred)*
- The manual never gives a coincidence tolerance (WELD, MERGESOIL, EXCAV, INTGEN, SLICE).
- Use one model-level tolerance `tol = 1e-6 × (model bounding-box diagonal)`, with a minimum of 1e-9,
  and allow a per-command override.
- EXCAV's `delta` argument is a separate **z-level clustering** tolerance.

### 1.6 Size limits [CORE]
- The maximum node number is **99,999** for the IKTR9 build, and 650,000 nodes / 2,500,000 DOF for
  IKTR9_650K (cross-ref §1.2).
- Commands that renumber or offset (MERGE, MERGESOIL, SOILMESH, CONVERT from ANSYS) can exceed this.
  NCOM is the cure.
- An implementation should warn when the max node number is > 99,999 while the target is the
  IKTR9-compatible writer.

---

## 2. §9.7 UI General commands

### 2.1 ACTM — change the active model [CORE][UI]
`ACTM,<Model>`
- `<Model>`: integer model reference number.
- Makes `<Model>` the active model. All later commands act on it.
- If the model exists, its data is **unchanged**.
- If it does not exist, ACTM **creates an empty model**: no nodes, no elements, and **default
  simulation settings**. *(manual)*
- Example (§5.7.1): `Actm,1` / `Inp,Model1.pre` / `Actm,2` / `Inp,Model2.pre` / `Actm,3` / `Merge,1,2,0,0,0`.

### 2.2 AFWRBAT — split an SSI run by frequency subsets [IO][ADV]
`AFWRBAT,<splits>`
- `<splits>`: number of subsets into which the **SSI frequency set** is split.

What it does *(manual)*:
- Creates `<splits>` copies of the model, **each with a separate frequency subset, each in its own
  folder**. The folders can be moved to different computers.
- Writes **batch files** that run the necessary modules in each folder.
- Writes **another batch file** that combines the per-folder results afterwards.
- Folder names and locations come from the model name and path set by **MDL**. MDL (or the
  database) is required.
- Module executable paths come from the UI **module-location settings**.
- The ACS SASSI installation path must be the same on every machine, or the user edits the batch
  files.

Recommended implementation *(inferred)*:
1. Read the active model's SSI frequency list (the frequency set used by ANALYS with
   `<fopt>`=0). Split it into `<splits>` **contiguous** blocks whose sizes differ by at most 1.
2. For k = 1..splits, create the folder `<Path>/<Model>_<k>/`. Run AFWRITE into it with only block k
   in the frequency set. File names stay `<Model>.<ext>`.
3. Write `run_<k>.bat` (and a `.sh` twin for portability) that runs, in order, the modules enabled in
   AOPT: SITE → POINT → HOUSE → ANALYS. Each folder must regenerate FILE1/FILE3 for its own
   frequencies, or copy them. HOUSE output is frequency-independent and may be copied.
4. Write `combine.bat` in `<Path>`. It gathers each `FILE8` and calls COMBIN repeatedly: COMBIN takes
   exactly two inputs, renamed `FILE81` and `FILE82`, and writes `FILE8` (cross-ref §3, COMBIN). Fold
   pairwise until one `FILE8` holds every frequency, then optionally run MOTION, STRESS and RELDISP.
- See Q-G1 for the open points.

### 2.3 CPMODEL — copy the active model [UI]
`CPMODEL,<Mdl>`
- `Mdl`: destination model number.
- Copies the **active model** into model `<Mdl>`. Nothing is said about overwriting; overwrite with
  a full deep copy *(inferred)*.
- The active model does not change.
- Copy the name and path too. A later AFWRITE from the copy would then overwrite the original's
  files, so warn about this (Q-G2).

### 2.4 CRITFREQ — find frequencies where the interpolated ATF is wrong [CORE][IO]
`CRITFREQ,<tol>,<minfilter>,<TF>,<Var>`

| Arg | Meaning *(manual)* |
|---|---|
| `<tol>` | Percentage difference between TFU and TFI that causes a frequency to be added to the result. |
| `<minfilter>` | Percentage below the global maximum within which TFU–TFI differences are ignored. Peaks too far below the maximum are not considered. |
| `<TF>` | Full path of the transfer-function file **without** the `.TFU`/`.TFI` extension. Both `<TF>.TFU` (computed ATF at the SSI frequencies) and `<TF>.TFI` (interpolated ATF at the Fourier frequencies) are read. |
| `<Var>` | Name of the UI variable that receives the list of identified frequencies. |

Purpose:
- Automatically identifies frequencies where an **interpolated** ATF peak differs significantly from
  the **computed** ATF values near that peak's frequency.
- The identified frequencies should be **added to the SSI analysis frequency set**. This is the
  "automatic selection of additional SSI frequencies" capability (cross-ref line 430, item 10 at
  line 2900).

Recommended algorithm *(inferred; Q-G3)*:
1. Read the TFI amplitudes A_I(f_j) at the Fourier frequencies and the TFU amplitudes A_U(F_m) at the
   SSI frequencies. If the files are complex (amplitude and phase columns, MOTION `<cplx>`=1), use
   the amplitude column.
2. Let A_max = max_j A_I(f_j).
3. Find the local maxima of A_I: A_I(f_j) ≥ A_I(f_{j−1}) and A_I(f_j) > A_I(f_{j+1}).
4. Discard any peak with A_I(f_p) < (1 − minfilter/100)·A_max.
5. For each remaining peak f_p:
   - find the bracketing SSI frequencies F_m ≤ f_p ≤ F_{m+1};
   - take the reference A_ref = max(A_U(F_m), A_U(F_{m+1}));
   - if 100·|A_I(f_p) − A_ref| / A_ref > tol, add f_p to the result.
   
   Peaks that coincide with a computed frequency (|f_p − F_m| < ½ Δf_Fourier) differ by nothing and
   are naturally skipped.
6. Store the frequencies, ascending and formatted as strings, in `<Var>`. The user can then pass
   `@Var[i]` to frequency-set commands.

### 2.5 DMODEL — delete a model from memory [UI]
`DMODEL,<Mdl>`
- Removes model `<Mdl>` from UI memory **only**. Saved files and databases are not deleted.
  *(manual)*
- If `<Mdl>` is the active model, make the active slot an empty model with the same number, as ACTM
  would *(inferred; Q-G2)*.

### 2.6 ETYPEGEN — set element types explicitly for the whole model [CORE]
`ETYPEGEN,<type>`

| `<type>` | Effect *(manual)* |
|---|---|
| 0 | Element type defined by location relative to the ground-surface level. |
| 1 | Change **all** elements to Structure (ETYPE=1). |
| 2 | Change **all** elements to Excavated Soil (ETYPE=2). |

- Applies to SOLID and PLANE elements (structure or excavation), and to SHELL/TSHELL (1 = structure,
  2 = buried or sitting on soil, so its nodes become interaction nodes).
- Leave beams, springs and GENERAL elements untouched. ETYPE means nothing for them in INTGEN
  *(inferred)*.
- **Option A:** every element **not defined as structural is ignored** by the ANSYS interface. Run
  ETYPEGEN,1 on structure-only models before Option A. *(manual)*
- Interpretation of `<type>`=0 *(inferred; Q-G4)*: resolve the implicit classification into explicit
  values. Every SOLID/PLANE element lying entirely at or below `gelev` (all node z ≤ gelev + tol)
  gets ETYPE=2. Every other SOLID/PLANE gets ETYPE=1, and every SHELL/TSHELL gets ETYPE=1. This
  matches the AFWRITE rule. Provide `ETYPEGEN,0,reset` to restore the implicit 0 if wanted.
- **Caution:** `ETYPEGEN,2` on a model that contains shells turns them into buried shells. INTGEN
  then makes all their nodes interaction nodes.
- Typical use (§5.7.2): `etypegen,1` on the structure model and `etypegen,2` on the excavation model,
  then MERGESOIL.

### 2.7 FIXSLDROT — fix rotations at solid-only nodes [CORE]
`FIXSLDROT` (no arguments)
- Every node connected **only to SOLID elements** gets ROTX, ROTY and ROTZ fixed, i.e. the equivalent
  of `D,n,n,1,1,ROT`. *(manual)*
- Together with FIXSHLROT and FIXSPRROT it is an alternative to FIXROT. It better approximates the
  stiffness matrix of models generated in ANSYS.
- It also saves disk space and run time (cross-ref line 742).
- Extension *(inferred, optional)*: for 2D models, nodes connected only to PLANE elements get UY,
  ROTX, ROTY and ROTZ fixed.

### 2.8 FIXSHLROT — soft drilling springs at coplanar-shell nodes [CORE]
`FIXSHLROT,[stiff]`
- `stiff`: stiffness of the rotational springs added to shells to remove singularities. Default
  **10**.
- Adds a **soft rotational spring about the shell normal** (the in-plane or "drilling" rotation) at
  every node connected to **coplanar shells**. *(manual)*
- Why: the thin-plate SHELL (Kirchhoff) has no in-plane rotational stiffness, so the drilling DOF at
  a node whose shells all share one normal is singular. Oblique shells cannot be handled with D, so
  springs are used (cross-ref line 733–740).
- Recommended practice: the spring stiffness should be **no more than 10% of the shell element
  bending stiffness**. *(manual, line 739)*
- TSHELL already gets a small rotational stiffness inside HOUSE, so it needs no fix (cross-ref
  line 759).

Algorithm *(inferred where marked)*:
1. For each node n whose connected elements are all SHELL, compute the unit normals of those shells:
   n_e = (x_IJ... cross product of the two diagonals, normalised).
2. Call the node **coplanar** if every |n_e · n_ref| ≥ 1 − ε (use ε = 1e-4). The orientation sign
   is ignored.
3. Add a rotational spring K_rot = stiff · n nᵀ on the rotation DOF, connected to a **new coincident
   node with every DOF fixed**. The WELD note in the manual confirms these springs use coincident
   nodes (§4.9 below).
   - With global-axis uncoupled SPRING elements (`SC` constants are global and uncoupled; cross-ref
     §9.4.34), use scxx = stiff·n_x², scyy = stiff·n_y², sczz = stiff·n_z² and zero translational
     constants *(inferred; Q-G5)*. This equals k·nnᵀ exactly for axis-aligned shells.
   - A Python engine that supports coupled nodal stiffness should instead add the exact
     stiff·nnᵀ.
4. Put the springs in a new SPRING group with a new SC property. Use damping 0 unless configured.

### 2.9 FIXSPRROT — fix rotations of translational-spring nodes [CORE]
`FIXSPRROT` (no arguments). The manual also spells it FIXSPROT (line 742); accept it as an alias.
- Fixes the rotational DOF of translational spring nodes that are connected **only to springs**, or
  **only to springs and solids**. *(manual)*
- Node connected only to springs: the fixed DOF follow **the spring stiffness**. Fix every DOF whose
  summed spring constant over the attached springs is zero *(inferred: any of the 6, but at least
  the rotations)*.
- Node connected to springs and solids: the rotational DOF follow **the springs' rotational
  stiffness**. Fix ROTX, ROTY or ROTZ wherever the summed scxx, scyy or sczz is zero.

### 2.10 FIXROT — automatic rotation fixing [CORE]
`FIXROT,[Stiff]`
- `Stiff`: stiffness of the springs added to shells to remove singularities. Default **10**.
- Fixes rotations at nodes of **solids connected only to solids**, and the **in-plane rotations of
  shells connected only to shells**, by applying the D command automatically. *(manual)*
- **Shell faces not parallel to a global coordinate plane cannot be fixed with D.** For those, soft
  rotational springs are added in the shell plane with stiffness `Stiff`. *(manual)*
- Equivalent in effect to FIXSLDROT + FIXSHLROT + FIXSPRROT (cross-ref line 746).
- **Strongly recommended for thin-shell SSI models, especially those with oblique shells.**

Algorithm:
- Solid-only nodes: D on ROT, as in FIXSLDROT.
- Coplanar shell-only nodes with common normal n:
  - if |n · e_k| ≥ 1 − ε for some global axis k, fix the rotation about axis k (D on ROTk);
  - otherwise add the soft spring as in FIXSHLROT.
- Spring-only nodes: as in FIXSPRROT.
- Nodes where non-coplanar shells meet (wall–floor corners) are left alone: each shell's drilling DOF
  is restrained by the other shell's bending.

### 2.11 FRAMECOMBIN — combine animation frame files [UI][IO]
`FRAMECOMBIN,<op>,<num>,<InFile1>,…,<InFileX>,<Outfile>`

| Arg | Meaning |
|---|---|
| `op` | **0** = SRSS, **1** = sum, **2** = average |
| `num` | Number of input frames to combine. X = num. |
| `InFile1…InFileX` | Full paths of the input ASCII frame files |
| `Outfile` | Full path of the output frame file |

- Combines ASCII frame files from MOTION, RELDISP or STRESS (Table 3.2: `\ACC\ACC_time_fnum`,
  `\THD\...`, `\NSTRESS\stress_time_fnum_comp`, `\SOILPRES\press_...`, `\TFU\TFU_freq_fnum`,
  `\RS\RS##_freq_fnum`) into a new frame for animations.
- Each frame header must state the **number of rows and columns**. Current ACS SASSI writes this by
  default. Legacy frames must first be converted with **MODFRAMES**.

Algorithm:
- Check that all `num` frames have identical rows × cols. Copy the first column(s) that identify a
  node or element from the first file *(inferred; see Q-G6 on which columns are IDs)*.
- Combine each data cell (i, j) element-wise:
  - SRSS: sqrt(Σ_k v_k²)
  - sum: Σ_k v_k
  - average: (1/num)·Σ_k v_k
- Write the header and data in the same layout.

### 2.12 FRAMESEL — pick critical animation frames [UI]
`FRAMESEL,<tol>,<Acc>,<Var>`
- `tol`: percentage of the global maximum **below which local maxima are ignored** as critical frames.
- `Acc`: acceleration time-history file to process (the MOTION `.ACC` format).
- `Var`: name of the variable that stores the critical-frame list.

Algorithm *(inferred details)*:
- Take a = acceleration column and A_max = max|a|.
- Collect the indices j of local extrema of a: maxima where a_j ≥ a_{j±1}, and minima where
  a_j ≤ a_{j±1}. Keep those with |a_j| ≥ (tol/100)·A_max.
- Store their **1-based frame numbers**, matching the `fnum` of frame files (`ACC_00.000_00001`), in
  `<Var>`.
- Q-G7: frame number versus time value.

### 2.13 GCOM — compress group numbers [UI]
`GCOM` (no arguments)
- Renumbers groups 1..G with no gaps. **Relative order is kept.**
- Element numbers inside groups are **not** compressed (ECOMPR does that, per group).
- Update every reference to a group number: cuts, nonlinear panel data, and so on. Group titles
  move with their groups *(inferred)*.

### 2.14 GETENV — show fast-solver environment variables [UI][ADV]
`GETENV` (no arguments)
- Prints the fast-solver environment variables and their values. See SETENV.

### 2.15 GLB2LOC — transform nodes to a local Cartesian system [CORE]
`GLB2LOC,<Start>,<End>,<Stride>,<Sysno>`
- `Start`, `End`, `Stride`: node range. `Sysno`: local coordinate-system number (LOC or LOCAL).
- Only **Cartesian** local systems are supported. *(manual)*
- Any requested node already defined in another local system is **first converted to global**, then
  to system `Sysno`. *(manual)*

Formula *(standard)*:
- Let system s have origin O_s and rotation R_s, whose columns are the local unit axes expressed in
  global coordinates (the LOC Euler angles or LOCAL three-node definition, cross-ref §9.3.8/9.3.9).
- First x_g = O_old + R_old·x_old (skip this if the node is global).
- Then x_loc = R_sᵀ (x_g − O_s).
- Store x_loc as the node coordinates and tag the node with system s.
- Node numbers and connectivity do not change.

### 2.16 GRAVITY — set the gravity constant [CORE]
`GRAVITY,<grav>`
- Sets only the model's gravity constant: the first `HOUSE` argument and the HOUSE dialog's
  "Acceleration of Gravity" (the same value is shown in SITE).
- Other HOUSE options are untouched. *(manual)*
- Units must be consistent: 32.2 ft/s² for British units, 9.81 m/s² for SI.

### 2.17 GROUNDELEV — set the ground elevation [CORE]
`GROUNDELEV,<elev>`
- Sets only the ground-surface elevation: the second `HOUSE` argument (`gelev`). *(manual)*
- This value drives the implicit ETYPE, AFWRITE classification, EXCAV, INTGEN option 4 and
  MERGESOIL welding.
- Set it **before** MERGESOIL (§5.7.2).

### 2.18 GROUPMAT — one material per group [ADV]
`GROUPMAT` (no arguments)
- Creates a **new material for every group** and assigns it to all elements of that group. *(manual)*
- Each new material copies the original material of the **first element of the group**.
- The **original material list is deleted**. The new list holds only the group materials.
- Intended for the **Option NON (nonlinear) modelling** process.
- Number new material g' consecutively in group order, starting at 1 *(inferred)*.
- Groups with no material (SPRING, GENERAL) are skipped *(inferred)*.
- Excavation groups use soil layers, not materials, so skip them too *(inferred; Q-G8)*.

### 2.19 INTGEN — automatic interaction-node generation [CORE]
`INTGEN,<type>,[level skip]`

| `<type>` | Method *(manual)* |
|---|---|
| 0 | Set **all** nodes to non-interaction nodes (clears every interaction flag). |
| 1 | Embedded foundation, **Flexible Volume (FV)** (the "Direct" method) |
| 2 | Embedded foundation, Flexible Interface with Excavation Volume Boundary Nodes: **FI-EVBN**, also called the **Modified Subtraction Method (MSM)** (sometimes ESM) |
| 3 | Embedded foundation, Flexible Interface with Foundation–Soil Interface Nodes: **FI-FSIN**, also called the **Subtraction Method (SM)** |
| 4 | **Surface foundation**: interaction nodes only at the ground-surface level |
| 5 | Embedded foundation, **Fast FV (FFV)** with multiple layers of internal interaction nodes |

- `[level skip]`: **option 5 only**. Number of levels skipped between interaction-node levels.
  Default **1**.

Rules *(manual)*:
- Options 1–5 **never remove** interaction nodes defined earlier. The result is the **union** with
  the existing set. Use option 0 first for a clean regeneration.
- Options 1, 2, 3 and 5 need the excavation volume **explicitly** defined by ETYPE or ETYPEGEN. With
  default ETYPE=0 the command **does not work**: raise an error.
- **Buried shells** (SHELL/TSHELL with ETYPE=2): their nodes are interaction nodes.
- ETYPE=2 has **no effect for beam elements**.

Definitions used by the algorithm *(inferred from §1.5.1 and §2.1 definitions)*:
- **Excavation elements X**: SOLID elements (3D) or PLANE elements (2D) with ETYPE=2.
- **Excavation nodes N_X**: every node of an element in X.
- **Boundary faces B**:
  - 3D: build all faces of every element in X (a hex has 6 quad faces). Drop repeated nodes, so
    degenerate faces of prisms and pyramids become triangles or vanish if they have fewer than 3
    unique nodes. Key each face by its sorted unique node tuple. A face used by **exactly one**
    element of X is a boundary face.
  - 2D: use element edges the same way.
- **Boundary nodes N_B**: nodes of faces in B.
- **Top (ground-surface) faces B_top**: boundary faces whose nodes all satisfy |z − gelev| ≤ tol.
- **FSIN nodes N_FSIN**: nodes of the faces in B \ B_top, i.e. the lateral and bottom surface of the
  excavation, which is the foundation–soil or near-field–far-field interface.
- **z-levels**: the sorted distinct z values of N_X, clustered with tol.

Generation per option:

| Option | Interaction set added |
|---|---|
| 1 FV | N_X (all excavation nodes) ∪ buried-shell nodes |
| 2 FI-EVBN / MSM | N_B (lateral + bottom + ground-surface top face) ∪ buried-shell nodes |
| 3 FI-FSIN / SM | N_FSIN (lateral + bottom only; top-face-only nodes excluded) ∪ buried-shell nodes |
| 4 Surface | Nodes with \|z − gelev\| ≤ tol that belong to a structural element. No ETYPE needed. For surface foundations, FV = FI = FFV. |
| 5 FFV | N_B ∪ {nodes of N_X on the selected internal z-levels} ∪ buried-shell nodes |

Option 5 level selection *(inferred; Q-G9)*:
- Number the internal levels from the bottom of the excavation, j = 1..L−2 between bottom (j=0) and
  top (j=L−1).
- Select levels with j mod (skip+1) == 0. So skip=1 gives every other internal level and skip=0
  gives every level (identical to FV).

After generation:
- Set the INT code to interaction (`INT,n,n,1,1,0`).
- Warn about any interaction node that is not on a soil-layer interface, or that lies above
  `gelev` + tol.
- Warn if interaction-node numbering is not ascending bottom-up (this matters for the incoherency
  option).
- **Do not renumber** nodes.

Closed-form check for a regular box excavation of n_x × n_y × n_z hexes (useful as a verification
test):

| Option | Count |
|---|---|
| FV | (n_x+1)(n_y+1)(n_z+1) |
| EVBN | (n_x+1)(n_y+1)(n_z+1) − (n_x−1)(n_y−1)(n_z−1) |
| FSIN | EVBN − (n_x−1)(n_y−1) |
| FFV | EVBN + (n_x−1)(n_y−1) × (number of selected internal levels) |
| Surface | (n_x+1)(n_y+1) for a surface mat of the same plan grid |

Guidance (cross-ref §1.5.1):
- Runtime grows roughly with the 2nd–3rd power of the number of interaction nodes. Keep it below
  about 20,000.
- **ASCE 4-16 / SRP 3.7.2:** a validation against FV is required before SM, MSM or FFV is used in
  production.
- SM can become unstable at high frequency.
- MSM suits shallowly embedded nuclear islands. FFV or FV suits deeply embedded structures such as
  SMRs.
- Quarter models are **not** recommended for qualifying SM/MSM.

### 2.20 MDL — set model name and path [IO]
`MDL,<Model>,<Path>`
- `<Model>`: model **name** (string). `<Path>`: model directory.
- Replaces the old PREP database mechanism. The database stays optional.
- Commands that need a name and path (AFWRITE, AFWRBAT, RUN*, ANSYS default) fail without it.
- Every file written by WRITE or AFWRITE is `<Path>/<Model>.<ext>`.

### 2.21 MDLNAME — change the model name only [IO]
`MDLNAME,<name>`
- Changes the model name. The path and the title (TIT command) are unchanged.
- WRITE and AFWRITE then use `<name>.<ext>`.

### 2.22 NCOM — compress node numbers [CORE]
`NCOM` (no arguments)
- Renumbers nodes 1..N with **no gaps**, keeping their **relative order**.
- Updates element connectivity, beam K nodes included. Also update every other node reference
  *(inferred, required for correctness)*: D fixities, interaction codes, masses (MT/MR), forces
  (F/MM), LOCAL coordinate-system definition nodes, multiple-support foundation node ranges, output
  node lists, nonlinear data.
- Needed when the 99,999 node limit is exceeded, typically after an ANSYS conversion or
  MERGE / MERGESOIL / SOILMESH.
- Recommended practice: on a finished model run **RMVUNUSED first, then NCOM**.
- Optionally write a mapping file (old → new), as `.map` *(inferred convenience)*.

### 2.23 MODELLIST — list models in memory [UI]
`MODELLIST`: lists each model's reference number and name (if defined).

### 2.24 MODFRAMES — fix legacy frame headers [IO][UI]
`MODFRAMES,<cols>,<framelist>`
- `<cols>`: number of columns in the frame files.
- `<framelist>`: full path of the animation frame-list file. This is the same list used to load
  animation frames: extensions `zpani`, `thiani`, `contani`, `thani`; elsewhere also `.dispani`,
  `.tfiani`, `.impani`.

Behaviour:
- Reads every frame file named in `<framelist>` and replaces the **second number of each frame
  header** with `<cols>`.
- The rest of each frame file, and the frame-list file itself, are unchanged.
- The frame-list file's first line is a header and is skipped when listing frames (cross-ref §7.3.1).

### 2.25 RADIUS — excavation element radii for POINT [CORE][IO]
`RADIUS,<Scale>,<FileName>`
- `<Scale>`: scale factor for the radius file. `<FileName>`: full path of the output file.
- Writes the **radius of every excavation element**, plus the **average radius**, to `<FileName>`.
- Purpose: choose the POINT module's **"Radius of Central Zone"** for irregular excavation meshes.

Cross-ref POINT guidance (figure, printed p.128, viewed):
- uniform rectangular mesh of spacing h: r = **0.90 h**
- uniform triangular mesh: r = **0.85 h**
- 2D: r = **h**
- non-uniform meshes: use the **average** radius, and run sensitivity studies with the min, average
  and max values, enveloping responses if they differ.

Algorithm *(inferred; Q-G10)*:
- For each excavation SOLID, take the characteristic plan size h_e = sqrt(A_plan,e). A_plan,e is the
  element footprint area projected on the XY plane: the area of the convex hull of its nodes' (x, y).
- For PLANE elements, h_e = the element's horizontal width.
- Radius r_e = Scale · h_e.
- Write one line per element (`group element r_e`), then a final line with the average
  r̄ = (1/N)Σ r_e. Also print min and max (convenience).
- The recommended Scale is 0.9 for 3D quad meshes and 1.0 for 2D.

### 2.26 RMVUNUSED — delete unused nodes [CORE]
`RMVUNUSED` (no arguments)
- Finds the nodes used by elements (all element node lists, including beam K nodes *(inferred)*) and
  the interaction nodes.
- Removes every node that is **unused and non-interaction**.
- Does **not** compress numbers and does **not** change connectivity. Pair it with NCOM.
- Also drop masses, forces and fixities attached to removed nodes *(inferred)*.

### 2.27 SETENV — configure the fast solver [ADV]
`SETENV,<mem>`
- `<mem>`: memory limit for the fast solver, in **MB**.
- Recommended: **90–95% of physical RAM**. **More than 100% has been shown to give incorrect
  results.** *(manual)*
- Sets **3 environment variables** in the user's registry. They are per-user and persistent.
- Run once per user account after installation, and again only if physical RAM changes.
- The manual does not name the 3 variables. A Python re-implementation should keep a user-level
  config entry (e.g. `fast_solver.mem_mb`) and expose it to the solver (Q-G11).
- Validation: refuse values > 100% of detected RAM, and warn above 95%.

### 2.28 TPSD — target PSD for EQUAKE [IO][ADV]
`TPSD,<num>,<file>`
- `<num>`: input-file number. It matches the numbering of the other EQUAKE commands ACCIN, ACCOUT,
  RSIN and RSOUT.
- `<file>`: full path of the target PSD file.
- Used when EQUAKE's `<tpsd>` option is 1 (cross-ref §9.2, `EQUAKE,…,[tpsd]`, and §9.2.42, which is
  the same command).

---

## 3. §9.8 Model checking commands

None of these commands changes the model, **except USED**. Each prints to the **command history**.
The number of entries listed is capped by the Check Options **break message number** (§6.5.3:
messages after the break are not shown, but totals are still counted; the option is not saved
between sessions).

### 3.1 EXCSTRCHK — excavation interior nodes shared with the structure [CORE]
`EXCSTRCHK`
- Checks whether excavation **interior** nodes are also structure **basement** nodes. That is wrong
  for SASSI: structure and excavation must vibrate independently except at their common interface
  nodes.
- Prints the list of shared nodes. Any reported node **may cause incorrect SSI results**.
- Called "a very important SSI model check" in §1.5.1.

Algorithm *(inferred)*:
- interior = N_X \ N_B (same definitions as INTGEN).
- Report every interior node that also belongs to an element with ETYPE ≠ 2, or to any BEAMS,
  SPRING or GENERAL element.
- Needs explicit ETYPE. If it is implicit, resolve it as in ETYPEGEN,0 for the check only.

### 3.2 FIXEDINT — fixed interaction nodes [CORE]
`FIXEDINT`
- Lists every interaction node with **any fixed translational DOF** (UX, UY or UZ), or reports that
  there are none.
- The user must fix the model.
- **Strongly recommended** before production runs.

### 3.3 FREESPRING — free spring nodes [CORE]
`FREESPRING`
- Warns about every **unconstrained** node (no fixed DOF) that is connected **only to a spring**.
- Such a node has a singular stiffness unless it carries mass.

### 3.4 HINGED — possible unintended hinges [CORE]
`HINGED`
- Warns about places where an element type with **more active DOF** (BEAMS, SHELL, TSHELL: 6) meets
  one with fewer (SOLID: 3, or PLANE).
- Rotations are not transmitted at such a joint. The 6-DOF elements should **penetrate the solids
  along an edge or face** (extra massless beams or shells), so moments pass as force couples (cross-ref
  §1.5.1 item 3).
- The user decides whether to act.

Algorithm *(inferred)*:
- For each 6-DOF element e, count its nodes that are also SOLID or PLANE nodes.
- If **exactly one** of e's nodes touches solids (a single-point connection), flag that node.
- Also flag beam–shell joints where a beam meets a shell node and the beam's rotation about the shell
  normal (drilling) is unresisted: the beam axis is not in the shell plane and no other 6-DOF element
  at the node is out of the plane. These need a "tripod" of beams, per §1.5.1.
- **Strongly recommended** before production runs.

### 3.5 INTCOUNT — count interaction nodes [UI]
`INTCOUNT`
- Prints the number of interaction nodes. That number strongly drives run time.
- Interaction nodes are defined or deleted with INT and INTGEN.

### 3.6 KINT — beam K-nodes that are interaction nodes [CORE]
`KINT`
- Lists every beam K-node (orientation node) that is also an interaction node.
- These can give **incorrect results**, especially when the node is used only as a K-node.
- Report only. The model is unchanged.

### 3.7 USED — fix unused nodes [CORE]
`USED`
- Finds the nodes not used by any element and **fixes them**, the same as `D,n,n,1,1,ALL` (all 6
  DOF fixed).
- Interaction nodes not attached to any element should be reported rather than fixed: fixing them
  would trip FIXEDINT *(inferred)*.

---

## 4. §9.9 Model generation and combination commands

### 4.1 EXCAV — generate an excavation volume from a structure model [CORE]
`EXCAV,<model>,[delta]`
- `<model>`: model number in which the generated excavation model is stored. The source is the
  **active** model.
- `[delta]`: allowed z variation within a single level. Default **0**. It must be a **positive
  floating-point number**, otherwise the default is used.

What it does *(manual)*:
- Uses the **lowest z-level grid** of the active model as a template.
- From it, builds a **homogeneous mesh up to the ground surface**.
- Parts of the basement that "outcrop" below grade but do not reach the **bottom** z-level (for
  example a shallower wing) get **no** excavation volume.
- The ground elevation must be set correctly in the source model.
- Use `delta > 0` for floors whose nodes have slightly different z. Otherwise each small variation
  creates a separate level.
- Used together with **MERGESOIL** to join the excavation and structure models.
- Precondition (§5.7.2): the basemat geometry has an **identical node-layer mesh at every underground
  elevation**.

Algorithm *(inferred; Q-M1)*:
1. Gather the source nodes with z ≤ gelev + tol.
2. Cluster their z values into levels: sort the z values and start a new level whenever
   z − z_level_start > delta. Each level's elevation is the mean, or the minimum, of its cluster.
   Add gelev as the top level if it is missing.
3. **Template:** the faces at the lowest level z₀.
   - The 2D cells (quads or triangles) are taken from structural element faces lying at z₀: SHELL
     elements at z₀, or bottom faces of SOLIDs at z₀.
   - The template node set is those faces' nodes, using their (x, y).
4. For each pair of consecutive levels (z_k, z_{k+1}), for k = 0..L−2, and for each template cell,
   create one SOLID: an 8-node hex for a quad, a 6-node prism for a triangle using repeated nodes per
   Figure 9.1.
   - Bottom face nodes 1–4: (x, y, z_k).
   - Top face nodes 5–8: (x, y, z_{k+1}).
   - Keep the counter-clockwise order seen from above.
5. Create **new** nodes in the destination model for every (template node, level). Number them
   bottom-up, level by level, as the N command recommends for embedded models.
   - Coordinates are snapped to the level elevation, so delta clustering merges near-equal z.
   - MERGESOIL later welds coincident nodes to the structure.
6. Put all elements into one SOLID group with ETYPE=2.
7. Set MSET of each element to the soil-layer index whose depth range contains the element mid-height,
   if soil layers are defined. Otherwise leave MSET=1 and warn. Copy the source model's soil layers
   and ground elevation to the destination.
8. Report the number of levels and elements. Warn if a level is not a soil-layer interface.

### 4.2 MERGE — combine two models [CORE]
`MERGE,<Mdl1>,<Mdl2>,<X>,<Y>,<Z>`
- The result is stored in the **active** model. In the example, `Actm,3` comes first.
- Mdl1 data is kept as-is: same nodes, elements and so on.
- Mdl2 data is **appended with renumbering offsets**: node numbers, group numbers, material numbers,
  "etc." (§5.7.1).
- Mdl2 nodes are **translated** by (X, Y, Z). Rotate first with ROTATE if needed.
- Typical use: building **SSSI** (structure-soil-structure) models.

Algorithm (offsets *inferred* as "max ID in Mdl1"):
- Node offset = max node number of Mdl1. Apply the same rule to:
  - groups (group offset), with element numbers kept within each group;
  - materials, soil layers, real properties, spring properties, matrix properties;
  - local coordinate systems.
- Remap every reference in Mdl2's copy accordingly: connectivity, MSET/RSET, fixities, interaction
  codes, masses, forces.
- Coincident nodes are **not** merged. Use WELD or MERGESOIL for that.
- Analysis settings (gravity, gelev, frequencies, ...) come from Mdl1 *(inferred; Q-M2)*.

### 4.3 MERGEGROUP — merge groups [UI][ADV]
`MERGEGROUP,<dest>,[G1],…,[G10]`
- Merges up to **11 groups** (dest plus G1–G10) into `<dest>`.
- Dest elements keep their numbers.
- Each added group's elements are numbered after the destination's current maximum element number,
  as they are appended in argument order.
- G1–G10 are **removed** afterwards, which leaves gaps in group numbering. Recommend GCOM afterwards.
- **All groups must have the same element type.** Error otherwise.
- Used in panel modelling for Option NON (WALLFLR → PANELIZE → MERGEGROUP).

### 4.4 MERGEPANEL — merge a panel model back into the original [ADV]
`MERGEPANEL,<Panel>`
- `<Panel>`: model number of the **panel model**, created with WALLFLR/PANELIZE from a copy of the
  original.

Procedure:
1. `ACTM` to the original model.
2. Delete its shell group(s) (GDEL), since the panel shells replace them.
3. `MERGEPANEL,<panel>`.

Result:
- The original model's non-shell elements, with the panel model's **groups and materials appended
  to the end** of the original group and material lists.
- Nodes are shared, because the panel model came from the same node set. Copy the panel model's
  groups referencing the existing node numbers, and do not offset nodes *(inferred)*.

### 4.5 MERGESOIL — merge the structure and excavation models [CORE]
`MERGESOIL,<Struct>,<Soil>,[Mode],[Stiff],[Stiff2],[SepLevel],[Mapping]`

| Arg | Meaning | Default |
|---|---|---|
| `<Struct>` | model number of the structure (elements ETYPE=1) | required |
| `<Soil>` | model number of the excavation volume (elements ETYPE=2) | required |
| `[Mode]` | how interface nodes are joined: **0** = fully unbonded interface on all sides; **1** = fully bonded on all sides by **merging nodes**; **2** = **stiff spring** connection on all sides; **3** = stiff springs **below** `SepLevel` and **soft** springs above it | **1** |
| `[Stiff]` | stiffness of the stiff springs (Modes 2, 3) | **10^7** |
| `[Stiff2]` | stiffness of the soft springs above the separation level (Mode 3) | **10** |
| `[SepLevel]` | global Z of the soil separation level (Mode 3) | none |
| `[Mapping]` | mapping file of excavation-model nodes (old node, new node) | none |

The manual's Mode 0 text says "unbounded"; this means *unbonded*.

Behaviour *(manual)*:
- The result goes into the **active model**.
- **All soil-model elements become soil (excavation) elements** (ETYPE=2), and their **materials
  become soil layers**.
- When nodes are merged, **the higher-numbered node is no longer used**. Remove it with RMVUNUSED.
- Set the **ground elevation before** merging, so interface nodes weld properly.
- **Option AA:** MERGESOIL is **required** to join ANSYS-built structure and excavation models
  (converted `.cdb` → `.pre`). The mapping file must be named **`modelname_Excv.map`** (see Demo 7).

Algorithm *(inferred where marked)*:
1. Copy Struct into the active model unchanged.
2. Append Soil with node, group and other offsets as in MERGE, without translation. Convert Soil
   materials into soil-layer entries, the material → L layer mapping preserving order, and set the
   MSET of its elements to the new layer indices.
3. Interface pairs: for each soil node s (renumbered), find a structure node t at the same location
   (‖x_s − x_t‖ ≤ tol) with z ≤ gelev + tol. Use a spatial hash.
4. By Mode:
   - **0**: no connection. The duplicate nodes stay separate.
   - **1**: replace s by t in every soil-element connectivity, since t has the lower number. The
     interaction flag of s moves to t *(inferred)*. s becomes unused.
   - **2**: add a SPRING element (s, t) per pair with translational constants scx = scy = scz =
     Stiff, rotational constants 0, in a new SPRING group and a new SC property. Damping is
     unspecified; use 0 or configure it (Q-M3).
   - **3**: as Mode 2, but pairs with z_s > SepLevel get Stiff2, and pairs with z_s ≤ SepLevel get
     Stiff. Mode 3 needs SepLevel; error if it is missing.
5. If `[Mapping]` is given, write a text file with one `old_soil_node  new_node` line per excavation
   node. new_node is t for merged nodes in Mode 1, otherwise the offset number.

Manual example inconsistency:
- The §5.7.2 example is `MergeSoil,1,2,1,modelname_Excv.map`, which puts the map file in argument 4
  (`[Stiff]`).
- Accept a non-numeric string in any trailing optional position as `[Mapping]` (lenient parsing).
  Document the canonical form `MergeSoil,1,2,1,,,,modelname_Excv.map`.

### 4.6 ROTATE — rotate the model about a point [CORE]
`ROTATE,<x>,<y>,<z>,<rxy>,<ryz>,<rzx>`
- `<x>,<y>,<z>`: centre of rotation.
- `<rxy>`: rotation about **Z** in degrees. `<ryz>`: about **X**. `<rzx>`: about **Y**.
- Rotates node coordinates. Used with MERGE to orient models for SSSI.

Convention *(inferred; Q-M4)*:
- Right-hand rule, positive counter-clockwise seen from the positive axis.
- Apply in the order Rz(rxy), then Rx(ryz), then Ry(rzx): x' = c + Ry·Rx·Rz·(x − c).
- Beam K nodes rotate with all other nodes, so beam orientation is preserved.
- **Global-axis spring constants (SC) and global-axis GENERAL matrices are NOT rotated.** Warn when
  the model has anisotropic springs or 2-node GENERAL elements.

### 4.7 SOILMESH — generate a near-field soil mesh for soil-pressure models [CORE][ADV]
`SOILMESH,<dest>,<sX>,<sY>,<hori>,<vert>,<xAdj>,<yAdj>,<Zdepth>,<contact>,<rNum>`

| Arg | Meaning *(manual)* |
|---|---|
| `dest` | destination model number |
| `sX` | float: percentage growth in X of each level |
| `sY` | float: percentage growth in Y of each level |
| `hori` | number of new horizontal layers |
| `vert` | number of vertical layers |
| `xAdj` | float: centroid correction in X |
| `yAdj` | float: centroid correction in Y |
| `Zdepth` | float: thickness of each new level in Z |
| `contact` | 0 = no contact surfaces. Otherwise contact surfaces are included. |
| `rNum` | real (spring) property number for the contact surface |

- Builds a **soil-field mesh** around the active model's basement, used to compute **pressure on the
  excavation volume**. This is the near-field soil SOLID layer that STRESS needs for soil-pressure
  frames.
- The basement must be a **continuous convex hull with a uniform external wall at all layers**.

Note *(manual)*:
- SOILMESH predates MERGE. Recommended practice: **WRITE the SOILMESH model, then load it (INP) over
  the original model**.
- MERGE also works, but both MERGE and SOILMESH apply numbering offsets (groups, elements, nodes,
  ...). Using both **doubles the offset** and leaves gaps.
- So the destination model contains only the new soil-mesh entities, already numbered after the
  original model's maxima.

Algorithm *(inferred; Q-M5)*:
1. Take the basement's outer wall perimeter at each underground level: the convex hull of the
   active model's nodes at that level. Take the centroid c = mean(perimeter) + (xAdj, yAdj).
2. Ring i = 1..hori around the perimeter: each perimeter node p maps to
   - x = c_x + (p_x − c_x)(1 + i·sX/100)
   - y = c_y + (p_y − c_y)(1 + i·sY/100)
   
   Growth may compound instead, (1+sX/100)^i (Q-M5). This gives concentric rings joined into hex
   elements between rings i−1 and i, at the existing basement level elevations.
3. vert layers of thickness Zdepth: extend the mesh below the basemat by `vert` levels at
   z_bottom − j·Zdepth, j = 1..vert. Cover the basemat footprint as well, so that soil exists under
   the mat.
4. New nodes, groups and elements are numbered from (original max + 1). The soil elements are
   SOLIDs with ETYPE=1, since near-field soil is modelled as part of the structure, with materials
   from the soil layers at the same elevation.
5. If `contact` ≠ 0: do not share nodes between the soil and the basement walls/mat. Create duplicate
   soil-side nodes and connect each pair with SPRING elements using property `rNum`, the contact
   surface. Spring axial forces divided by tributary area estimate pressures (cross-ref line 2969).
   If `contact` = 0, share the basement nodes.

### 4.8 TRANSLATE — translate every node [CORE]
`TRANSLATE,<x>,<y>,<z>`
- Moves **all** nodes by (x, y, z). Each defaults to **0**.
- No other model data changes. The ground elevation is not shifted, so warn if the model is moved
  vertically.

### 4.9 WELD — merge coincident nodes [CORE]
`WELD`
- Searches the node list for nodes at the **same location** and builds a **substitution table**
  (keep the lowest node number *(inferred)*).
- Rewrites element connectivity with that table.
- Designed for **excavation-only or structure-only models** before MERGESOIL.

Notes *(manual)*:
- Substituted nodes stay in the node list but are unattached. Clean up with RMVUNUSED.
- Masses, boundary conditions and interaction status are **not** changed. The user may need to fix
  them.
- WELD changes connectivity at **all** coincident nodes, even intended ones:
  - separate structure and excavation meshes that share coordinates;
  - the **zero-length springs created by FIXROT**, which become degenerate.
  
  Implementation: warn if any element would end up with repeated nodes when that is not a legal
  degenerate SOLID/SHELL form, and refuse to collapse SPRING elements unless `force` is given.

---

## 5. §9.10 Cut and submodeling commands

### 5.1 Cut data structure [CORE][UI]
- A cut is a numbered set of (group, element) pairs.
- CUTADD and CUTVOL **create** the cut if it does not exist.
- CUTCLR deletes it. After that, any plot or submodel use of the cut is an **error**.
- SLICE now **adds** to an existing cut. In the older SUBMODELER module it replaced the cut.
- Workflow (§5.8):
  1. Select elements with CUTADD, CUTVOL or SLICE.
  2. Review with Plot > Cuts / CUTPLOT.
  3. Build a submodel (CUT2SUB) or a cross-section (CSECT).
  4. Run calculations (CALCPAR, CALCM, CALCSECTHIST in §9.13) after loading element stresses (READSTR /
     `.ess` frames).
- Transfer commands (TRANELEM, TRANVOL) bypass the review step. The manual recommends cuts instead.

### 5.2 CSECT — cross-section model [CORE]
`CSECT,<dest>,<cutnum>,<pointx>,<pointy>,<pointz>,<normalx>,<normaly>,<normalz>`
- The infinite plane passes through point P = (pointx, pointy, pointz) with normal n = (normalx,
  normaly, normalz). n need not be unit length: normalise it.
- Takes the cross-section of **every element in the cut** with that plane, then **expands it to a
  unit thickness** to form a 3D model that can be plotted. This model is stored in `<dest>`.

Details (§5.8.2):
- The cross-section model keeps the **same group and element numbers** as the original. Node numbers
  and locations are new.
- Each element's extent is trimmed to its intersection with the plane.
- Shell thickness is **adjusted when the shell is oblique** to the cut plane.
- Load the `.ess` element stresses **before** creating the cross-section.
- Example: `readstr,stressfile.txt` / `cutvol,3,52.5,52.8,-320,320,2.53,45.22` /
  `csect,1,3,55.6,0,4,1,0,0` / `actm,1` / `calcpar,1,0,0,0,1,0` / `calcm`.

Algorithm *(inferred; Q-C2)*:
- For each element compute the signed distances d_i = n̂·(x_i − P) of its nodes. Keep the element
  only if min d < 0 < max d, or if |d| ≤ tol at some node.
- **Shells:** the shell mid-surface meets the plane in a segment AB, found by interpolating along
  edges whose ends have d values of opposite sign.
  - Build a 4-node SHELL with nodes A − ½n̂, B − ½n̂, B + ½n̂, A + ½n̂ (unit thickness along n̂).
  - Its thickness is t_eff = t / sinθ, where θ is the angle between the shell plane and the cut plane
    (sinθ = ‖n̂ × n_shell‖). This is the true width of the cut through the shell.
- **Solids:** intersect with the plane to get a polygon of 3–6 points. Extrude it ±½ along n̂ and
  split it into hex or prism SOLIDs. Keep the original element number for the first sub-element
  (Q-C2 covers multiple sub-elements).
- **Beams:** the intersection point Q gives a BEAMS element from Q − ½n̂ to Q + ½n̂, keeping the
  K-node orientation.
- Copy materials, properties and the element stress data with the elements.

### 5.3 CUT2SUB — submodel from a cut [CORE]
`CUT2SUB,<cutnum>,<dest>,[solid]`
- Builds a submodel in `<dest>` from the active model's elements in the cut: those elements with the
  same group and element numbers, every node they reference, and the materials and properties they
  use.
- `[solid]`: if ≥ 1, **convert all shells into solid elements** ("thick shells"), by extruding ±t/2
  along the shell normal. Default **−1** (no conversion).
- The thick-shell solids are **for plotting only**: the generator cannot check whether the new solids
  intersect.
- Example: `cutvol,3,,,,,2.53` / `cut2sub,3,1` / `actm,1` / `write,model above ground.pre,<dir>`
  writes the elements above a ground surface at 2.53.

### 5.4 CUTADD — add elements to a cut [UI]
Two forms:
- `CUTADD,<cutnum>,<group num>,<elem 1>,…,<elem N>`: an explicit, non-sequential list.
- `CUTADD,<cutnum>,<group num>,RANGE,<elem start>,[elem end],[stride]`: a range.
  - `elem end` defaults to `elem start`. `stride` defaults to 1.
  - The literal keyword `RANGE` (case-insensitive) in field 3 selects this form.
- All elements come from one group. Creates the cut if needed. Ignore duplicates; warn on elements
  that do not exist.

### 5.5 CUTCLR — clear cuts [UI]
`CUTCLR,<first cut>,[last cut],[step]`
- `last cut` defaults to `first cut`. `step` defaults to 1.
- Deletes the cuts entirely. They can be rebuilt by adding elements again.

### 5.6 CUTRMV — remove elements from a cut [UI]
Same two forms as CUTADD:
- `CUTRMV,<cutnum>,<group num>,<e1>,…,<eN>`
- `CUTRMV,<cutnum>,<group num>,RANGE,<start>,[end],[stride]`

Defaults: end = start, stride = 1.

### 5.7 CUTVOL — add the elements in a box to a cut [UI]
`CUTVOL,<cutnum>,[Xmin],[Xmax],[Ymin],[Ymax],[Zmin],[Zmax]`
- Any blank bound is replaced by the active model's corresponding minimum or maximum.
- Selection criterion *(inferred; Q-C3)*: an element is selected if **all its nodes** lie inside the
  closed box, with tol.
  - Check: the manual's `cutvol,3,,,,,2.53` selects the elements above ground at 2.53.
  - Check: `cutvol,3,52.5,52.8,-320,320,2.53,45.22` isolates one wall lying in a thin X slab.

### 5.8 EXTRACTEXCAV — submodel of the excavation volume [CORE]
`EXTRACTEXCAV,<Model>`
- Copies every **excavation** element (ETYPE=2 SOLID/PLANE, which must be **explicit**; use
  ETYPEGEN), with its nodes, soil layers and interaction flags, into model `<Model>`.
- The active model is unchanged.

### 5.9 SLICE — add the elements crossing a plane to a cut [UI]
`SLICE,<cutnum>,<pointx>,<pointy>,<pointz>,<normalx>,<normaly>,<normalz>`
- Adds to the cut every element that **crosses** the infinite plane through the point with that
  normal.
- Test: the signed node distances have both signs, or one node lies within tol of the plane.

### 5.10 SPLITGROUP — split a group by a plane [ADV]
`SPLITGROUP,<group>,<split>,[dir]`
- `group`: the group to split. It can be any group of any element type.
- Without `dir`:
  - `split` is the number of a **SHELL group** whose **first shell** defines the splitting plane.
  - The second group must be shells.
  - Designed to work with **WALLFLR** groups.
  - Uses the same algorithm PANELIZE uses to split groups.
- With `dir` = `X`, `Y` or `Z` (added in UI 1.0.2):
  - `split` is a float giving the plane's location in global coordinates;
  - the plane is perpendicular to the chosen axis.
- Only `<group>` is split.

Algorithm *(inferred)*:
- Classify each element by the sign of its **centroid's** signed distance to the plane.
- Elements on the negative side stay in `group`. Elements on the positive side move to a new group,
  number = max group + 1, of the same type, with the title suffixed `_B`. Renumber them 1..k in their
  original order.
- Elements whose centroid lies on the plane (|d| ≤ tol) stay.
- Warn if either side is empty.

### 5.11 TRANELEM — transfer a list of elements to another model [UI]
`TRANELEM,<dest>,<group>,<begin>,<end>,<stride>`
- Copies elements begin..end (step `stride`) of `group` in the active model **directly** into model
  `<dest>`.
- Destination elements with the same numbers are **overwritten**.
- The destination group's **type is changed** to the source type.
- **Every node referenced is overwritten** in the destination as well.
- The active model is unchanged.

### 5.12 TRANVOL — transfer the elements in a box to another model [UI]
`TRANVOL,<dest>,[Xmin],[Xmax],[Ymin],[Ymax],[Zmin],[Zmax]`
- Transfers every element inside the box (global coordinates) to `<dest>`. Blank bounds mean the
  model extents.
- Same overwrite semantics as TRANELEM: destination group types, elements and nodes change.
- Use the same inside test as CUTVOL.

---

## 6. §9.11 File conversion commands

### 6.1 ANSYS — export to ANSYS APDL [IO]
`ANSYS,[FileName],[Dir]`
- Writes the active model as an ANSYS **APDL** input file. The manual writes "ADPL".
- Default name and location: `<modelname>.inp` in the model directory (cross-ref §6.1.9). With no
  name or path defined, the file is called `.inp` in a platform-dependent working directory.
- Not every model translates directly. **Use ANSYSREFORMAT first** to avoid most problems.
- The export uses elements compatible with **ANSYS V11–15**, which newer ANSYS treats as legacy, so
  output may fail in the latest versions.
- Element mapping *(inferred, the inverse of the importer list in §6.1.6; Q-F1)*:

| ACS SASSI | ANSYS | Notes |
|---|---|---|
| SOLID | SOLID45 | |
| SHELL / TSHELL | SHELL63 | real constant = thickness |
| BEAMS | BEAM4 | releases through KEYOPT 7/8, hence ANSYSREFORMAT |
| SPRING | COMBIN14 | one element per non-zero direction; KEYOPT(2)=1..6 |
| node masses | MASS21 | |
| GENERAL | MATRIX27 | |

  Write MPDATA for materials, R/RBLOCK for real constants, D for constraints.

### 6.2 ANSYSMODELTYPE — Option AA model type [ADV]
`ANSYSMODELTYPE,<type>`
- `1` = embedded operation mode. `2` = surface operation mode.
- Same as the HOUSE dialog's "ANSYS Model Input" check box plus "Ansys Model Type" (Embedded or
  Surface). The manual recommends the dialog to avoid input mistakes.
- Option AA works only with the fast solver (HOUSEFSA / ANALYSFA).
- In Option AA, the `.hou` written by AFWRITE is **fake**: it holds topology only and cannot be run
  without Option AA.

### 6.3 ANSYSREFORMAT — regroup beams for ANSYS export [IO]
`ANSYSREFORMAT,<Org>,<Map>`
- `<Org>`: model number to reformat. `<Map>`: file name for the mapping file between old and new
  beam groups.
- Run it with an **empty active model**; the reformatted copy goes there.
- Reason: ACS SASSI sets beam end releases **per element** (KI/KJ). ANSYS sets them **per element
  type** (KEYOPTs), so beam groups have to be regrouped.

Algorithm *(inferred)*:
- Copy Org into the active model.
- For each BEAMS group, partition its elements by the key (KI code, KJ code), and create one new
  BEAMS group per distinct key. A group with a single key keeps its number. Elements are renumbered
  1..k within each new group.
- Write `<Map>` with lines `old_group old_elem new_group new_elem`.

### 6.4 CONVERT — run a file converter without the GUI [IO]
`CONVERT,<ConSel>,<model>,<filename>,…`

| `<ConSel>` | Converter |
|---|---|
| `SSI` | SASSI `*.hou` → model. Also reads `*.sit` and `*.poi` with the same base name in the same folder. |
| `ANSYS` | ANSYS `*.cdb` → model. Extra argument: `CONVERT,ANSYS,<model>,<filename>,<gravity>`. |
| `STRUDL` | GT STRUDL → model. **Not applicable in this version.** Extra arguments: `CONVERT,STRUDL,<model>,<jointfile>,<group>,<member>,<fematrib>,<prop>,<const>,<input>,<tie>` |

- `<model>`: destination model number. `<filename>`: full path of the file to convert.
- The STRUDL extra files are the group database, member database, FEM attribute file, properties,
  constraints, mass input and spring properties.
  - The joint file goes in `filename` with its full path.
  - The others are names only, in the **same directory**.
  - Very limited testing. **Not for NQA use.**

`CONVERT` with no arguments opens the converter dialog (fields read from the screenshots on printed
pp. 91–92):

| Dialog | Field | Notes |
|---|---|---|
| SASSI .hou to .pre Converter | Input File Name | required (browse "<<") |
| | Output .pre File Name | optional; if given, a `.pre` is written |
| | Save Converted Data to Model Number | non-negative integer; blank means the active model |
| | Buttons | Convert, Cancel |
| ANSYS .cdb to .pre Converter | Input File Name | required |
| | Output .pre File Name | optional |
| | Save Converted Data to Model Number | blank means the active model |
| | Enter Value for Gravity | **required**; the example shows 32.2 |

Both dialogs carry the disclaimer that the converter has had limited testing, so check the model.

- Without an output name, the model lives in memory only. Use WRITE later.
- `.cdb` conversion limits (cross-ref §6.1.6), summarised:
  - materials only via MPDATA;
  - BEAM4/44 need a K node and RBLOCK with 6, 8, 10, 12 or 19–24 fields; releases via KEYOPT 7/8;
  - COMBIN14 needs KEYOPT 2 (1–6) or KEYOPT 3 (1–2), with KEYOPT 2 taking precedence, and stiffness in
    RBLOCK;
  - MASS21 needs RBLOCK mass with KEYOPT 1 and 2 = 0;
  - SOLID45 and SOLID185 (KEYOPT4 = 0);
  - SHELL63 (thickness in RBLOCK) and SHELL181 (thickness via section commands);
  - BEAM188/PIPE288: K node, ASEC or RECT sections, no end releases, no pipe offsets.
- In Option AA, display-only element types are also accepted: TRUSS180 and PIPE16/18 as beams,
  MPC184 as a spring, FLUID80 as a solid. The converted model is flagged as not runnable by plain
  HOUSE.

### 6.5 GENMATRIXDAMP — complex stiffness for converted MATRIX27 [IO] (NOT USABLE IN THIS VERSION)
`GENMATRIXDAMP,<begin>,[end],[stride],<damp>`
- `<begin>`: first matrix property to convert.
- `[end]`: last matrix property. Blank or default (−1) means end = begin.
- `[stride]`: default 1.
- `<damp>`: user damping value.

Behaviour:
- ACS SASSI reads MATRIX27 connectivity and mass from `.cdb`. ANSYS has no complex stiffness, so this
  command builds one from the real stiffness and `damp`.
- Formula *(inferred; Q-F2)*: Im K = 2·damp·Re K, i.e. K* = K(1 + 2iβ), with damp a ratio (not %).
  Write it into MXI for each property in the range. Keep the project-wide hysteretic-damping
  convention used for element materials.
- In this release the command should print "not available in this version" and do nothing, matching
  the manual.

---

## 7. Typical workflows (from §5 cross-refs) [UI]

**Build an SSSI model:**
```
ACTM,1
INP,Model1.pre
ACTM,2
INP,Model2.pre
ROTATE,...        (if needed, on model 2 while active)
ACTM,3
MERGE,1,2,dx,dy,dz
```

**Structure plus excavation (including Option AA):**
```
ACTM,1
INP,struct.pre
ETYPEGEN,1
GROUNDELEV,g
ACTM,2
INP,soil.pre      (or EXCAV,2 from model 1)
ETYPEGEN,2
GROUNDELEV,g
ACTM,3
MERGESOIL,1,2,1,,,,modelname_Excv.map
RMVUNUSED
NCOM
INTGEN,0
INTGEN,<method>
FIXROT
```
Then run the model checks: EXCSTRCHK, FIXEDINT, HINGED, KINT, FREESPRING, INTCOUNT.

**Section cut:**
```
INP,model.pre
READSTR,file
CUTVOL,...        (or CUTADD / SLICE)
CSECT,dest,cut,P,n
ACTM,dest
CALCPAR,...
CALCM
```

**Frequency refinement:**
1. Run MOTION.
2. `CRITFREQ,tol,minf,<path>/xxxxxTR_X,FR`
3. Add `@FR[i]` to the frequency set.
4. Run ANALYS for the new frequencies.
5. COMBIN.

---

## 8. Verification tests to implement (for "verifiable")

1. **INTGEN counts** on a generated box excavation of n_x × n_y × n_z hexes: compare against the
   closed forms in §2.19 for options 1, 2, 3, 5 (skip = 0, 1, 2) and 4. Check that options 1–5 are
   additive and that option 0 clears.
2. **EXCAV:** a basement with a 4×4 basemat grid and 3 embedment levels with z jitter of ±0.001.
   - delta = 0 gives extra levels.
   - delta = 0.01 gives 48 hexes and 75 nodes.
   - Every node lies in the template (x, y) set.
3. **MERGESOIL Mode 1:** coincident interface nodes collapse to the lower number. The higher numbers
   become unused. The map file lists every excavation node.
4. **MERGESOIL Modes 2 and 3:** the spring count equals the interface-pair count, and stiffness
   splits correctly at SepLevel.
5. **NCOM after RMVUNUSED:** node numbers become 1..N, order is kept, and connectivity, fixities,
   masses and interaction flags are remapped. A nodal mass at the old node 57 must end up at the
   right new node.
6. **ROTATE/TRANSLATE:** (1, 0, 0) rotated about the origin by rxy=90 gives (0, 1, 0). A
   round-trip rotate of +θ then −θ restores the coordinates to 1e-12.
7. **GLB2LOC:** a node in system 1 (origin (1, 2, 3), rotated 30° about Z) transformed to system 2
   matches the analytic coordinates.
8. **FIXROT:** in a model of one flat shell plate in the XY plane, every shell-only node gets ROTZ
   fixed and nothing else. The same plate rotated 30° about X gets springs, not D.
9. **WELD:** two coincident 8-node cubes become a 12-node connected mesh, and 4 nodes become unused.
10. **CRITFREQ:** synthetic TFU with a sharp peak between two SSI frequencies plus a TFI peak
    overshoot of 40%. tol = 20 flags the peak; tol = 50 does not.
11. **CUTVOL/SLICE/CSECT:** on a 1×1×1 cube of 2×2×2 solids, `SLICE` through x = 0.5 selects all 8;
    CSECT gives a 2×2 patch of unit-thickness solids with the original element numbers.

---

## 9. Manual errata and inconsistencies noted

| Location | Issue |
|---|---|
| §9.7 table | GCOM is listed as "Section 9.7.11"; it is actually 9.7.13 (9.7.11 is FRAMECOMBIN). |
| §9.8 table | INTCOUNT is listed as 9.8.4 (actually 9.8.5); KINT as 9.8.5 (actually 9.8.6). |
| §9.9 table | MERGEGROUP is described as "Merge 2 groups", but it merges up to 11. |
| §9.10 table | TRANELEM is listed as "Section 0", a broken cross-reference; it is actually 9.10.10. |
| §9.9.5 | MERGESOIL Mode 0 says "unbounded"; this means *unbonded*. The Stiff default prints as "107" in the extracted text; the PDF image shows **10^7**. |
| §5.7.2 | The example `MergeSoil,1,2,1,modelname_Excv.map` puts the mapping file in the 4th slot, while the syntax puts it in the 7th. |
| §5.8.2 | The CSECT example plane x = 55.6 does not fall inside the CUTVOL box X ∈ [52.5, 52.8]. Probably a typo; do not use it as a test. |
| §1.5.1 | "FIXSPROT" is used for FIXSPRROT. |
| §9.7.19 | "Type of iteration node generation" means *interaction*. |
| §9.11.1 | "ADPL" means APDL. |

---

## 10. Implementation-critical summary

- **ETYPE must be explicit** (via ETYPEGEN or ETYPE) for INTGEN 1/2/3/5, EXTRACTEXCAV and MERGESOIL.
  Raise a clear error when it is implicit.
- **INTGEN is additive.** Use option 0 to clear.
  - FV = all excavation nodes.
  - EVBN = all excavation boundary nodes, top face included.
  - FSIN = lateral + bottom boundary nodes.
  - FFV = EVBN + selected internal levels (skip default 1).
  - Surface = ground-level nodes.
  - Nodes of buried shells (ETYPE=2) are always added.
- **MERGESOIL:**
  - default Mode 1: merge nodes, keeping the lower number;
  - Stiff default 1e7, Stiff2 default 10, SepLevel in global Z;
  - soil materials → soil layers, soil elements → ETYPE=2;
  - the map file is `modelname_Excv.map` for Option AA;
  - set gelev first.
- **MERGE:** Mdl1 unchanged; Mdl2 offset by Mdl1 maxima and translated; result in the active model;
  no welding.
- **RMVUNUSED then NCOM** to respect the 99,999-node limit. NCOM must remap **every** node
  reference.
- **FIXROT family:**
  - D on rotations at solid-only nodes;
  - drilling DOF fixed via D when the shell normal is parallel to a global axis;
  - otherwise a soft spring (default 10, ≤ 10% of the shell bending stiffness) to a fixed coincident
    node;
  - TSHELL needs no fix.
- **Model checks** (EXCSTRCHK, FIXEDINT, HINGED, KINT, FREESPRING) are read-only and capped by the
  Check break number. USED fixes unused nodes with D ALL.
- Commands in §9.6 onward need **full names**. Blank fields mean defaults.

---

## 11. Open questions / ambiguities

**General (§9.7)**
- **Q-G1 (AFWRBAT):** The manual leaves these unspecified, so the implementer must choose:
  - folder naming;
  - which modules each batch runs (SITE/POINT/HOUSE/ANALYS or ANALYS only);
  - whether the split is contiguous or interleaved;
  - how the combine batch chains COMBIN (pairwise FILE81/FILE82).
  
  Recommended: contiguous blocks, `<name>_<k>` folders, full SITE→ANALYS per folder, pairwise COMBIN
  fold.
- **Q-G2 (CPMODEL/DMODEL):** Unspecified:
  - whether CPMODEL overwrites an existing destination;
  - whether the model name and path are copied (this risks file collisions);
  - what happens when DMODEL deletes the active model.
- **Q-G3 (CRITFREQ):** Unspecified:
  - the exact peak detection;
  - the meaning of `minfilter` (threshold (1 − m/100)·A_max, or m/100·A_max);
  - the "vicinity" comparison (bracketing SSI frequencies, max or interpolated);
  - TFU/TFI column layout (amplitude only, or amplitude + phase in radians when `<cplx>`=1);
  - the output format (frequencies or frequency indices).
- **Q-G4 (ETYPEGEN,0):** Does it reset elements to implicit 0, or explicitly assign 1/2 by location?
  What is the "below ground" test: centroid, or all nodes? Does ETYPEGEN touch beams?
- **Q-G5 (FIXSHLROT/FIXROT springs):** SC springs are global and uncoupled, so how is a spring
  "about the normal" represented for oblique shells?
  - Options: diagonal k·n_i², k on all three rotations, or a coupled matrix.
  - Is a duplicate fixed node created, or does the spring connect to ground some other way?
  - What damping do the springs carry?
- **Q-G6 (FRAMECOMBIN/MODFRAMES):** The exact frame header (the "number of rows and columns"; MODFRAMES
  changes the *second* header number) and which columns are IDs versus data are not given.
  - Do SRSS and average apply to the ID columns?
  - The frame text formats (Table 3.2) give file naming only, not content.
- **Q-G7 (FRAMESEL):** Is the output frame numbers or times? Is it 1-based? Do minima count by
  absolute value? Which column of the `.ACC` file holds acceleration?
- **Q-G8 (GROUPMAT):** How are groups without materials (SPRING, GENERAL, excavation soil layers)
  handled, and how are new materials numbered?
- **Q-G9 (INTGEN):**
  - Option 5: the exact meaning of `level skip` (default 1 = every other level?) and the counting
    origin (bottom or top).
  - Option 3: is FSIN the geometric lateral + bottom boundary, or strictly the nodes shared with
    structural elements? They are the same when no backfill is modelled.
  - Option 4: which nodes qualify (all nodes at gelev, or only structural foundation nodes)?
  - Does INTGEN set code 0 (interaction) only?
- **Q-G10 (RADIUS):** The element-radius formula and the role of `Scale` are not given. The file
  format is not given either.
- **Q-G11 (SETENV/GETENV):** The names and meanings of the 3 fast-solver environment variables are
  not given. A Python port should define its own configuration keys.

**Model generation (§9.9)**
- **Q-M1 (EXCAV):**
  - how the template grid is extracted (from shells at the lowest level, or from solid bottom faces);
  - how level elevations are chosen inside a delta cluster;
  - node numbering;
  - soil-layer (MSET) assignment;
  - whether gelev must coincide with a structure node level.
- **Q-M2 (MERGE):**
  - offset rule (max ID, or max ID + 1);
  - which entity tables are offset;
  - which model supplies the analysis settings, soil layers and frequency set;
  - whether duplicate soil-layer tables are merged.
- **Q-M3 (MERGESOIL):**
  - spring damping;
  - which node keeps the interaction flag in Mode 1/2/3;
  - whether welding is limited to z ≤ gelev;
  - the boundary case z = SepLevel;
  - the mapping-file format;
  - argument-position leniency, given the §5.7.2 example.
- **Q-M4 (ROTATE):** The rotation order and sign convention are not stated. Should local coordinate
  systems also rotate?
- **Q-M5 (SOILMESH):** The algorithm is described only by its parameters:
  - linear versus compounded growth;
  - meaning of `vert` versus `hori`;
  - whether the mesh extends under the basemat;
  - ETYPE and materials of the generated soil;
  - contact-surface element type (springs with RSET `rNum`) and node duplication;
  - how the "convex hull with uniform external wall" precondition is checked.
- **WELD tolerance**, and which node survives, are unspecified.

**Cut and submodel (§9.10)**
- **Q-C1:** Are cuts global (UI-level) or per model? The Cut Plot accepting any model number suggests
  global.
- **Q-C2 (CSECT):**
  - how solids are cut and extruded, and how multiple sub-elements keep the "same element numbers";
  - the beam cross-section representation;
  - whether the unit thickness is centred on the plane or one-sided;
  - the shell oblique-thickness formula.
- **Q-C3 (CUTVOL/TRANVOL):** Is "inside the box" all nodes, the centroid, or any node? Are bounds
  inclusive?
- **SPLITGROUP:** the new group number for the split-off part, the side classification (centroid or
  all nodes), and the treatment of elements straddling the plane are unspecified.
- **TRANELEM/TRANVOL:** Do materials, properties and soil layers referenced by the elements also
  transfer? Are destination elements beyond the transferred range kept?

**File conversion (§9.11)**
- **Q-F1 (ANSYS export):** The element type mapping, KEYOPT settings, unit handling and file layout
  of the `.inp` are not specified. Only "V11–15 compatible legacy elements" is stated.
- **Q-F2 (GENMATRIXDAMP):** The complex-stiffness formula, (1 + 2iβ) versus SASSI's
  (1 − 2β² + 2iβ√(1 − β²)), is not stated. Is `damp` a ratio or a percentage? The command is marked
  not usable in this version.
- **ANSYSREFORMAT:** the regrouping key and the map-file format are unspecified.
- **CONVERT,SSI:** the `.hou` fixed-format layout is defined in the HOUSE input specification, not
  here. The SASSI2000 legacy preparation rules (add zeros on the 3rd option line, remove `$`
  comments, insert a new first line with 0 in column 5 of `.sit`) must be honoured or automated.
