# Option A: SSI results as ANSYS loads (module LOADGEN)

Option A of ACS SASSI is the **two-step approach**: step 1 is the soil-structure interaction (SSI)
analysis; step 2 is a refined ANSYS model of the structure, loaded by the motion that step 1 computed
at the foundation-soil interface. Step 2 can use what SASSI cannot: finer meshes, other element types,
nonlinear materials, contact, cracking, a refined model of the soil for the pressures on basement walls.
The module that writes the ANSYS loads is **LOADGEN** (manual section 6.4.15). This guide shows how to
run it, what it writes, and how to use the result in ANSYS Mechanical APDL.

Related documents: [ANSYS.md](ANSYS.md) (moving models between ANSYS and SASSI-EDU), the
[User Guide](USER_GUIDE.md) (the SSI chain), the code in `sassi/modules/loadgen.py` and
`sassi/prep/commands/loadgen_cmds.py`, and the verification problems VP-LA1 and VP-LA2
(`sassi/verify/problems/vp_loadgen.py`).

> **Assumption of Option A (manual, Introduction).** Whatever you add in the ANSYS model must not change
> the motion of the foundation-soil interface. Option A is exact for a linear structure with the same
> mass, stiffness and damping as the SSI model (VP-LA2 shows this), and an approximation otherwise. If
> the second-step model is much softer or heavier than the SSI model, redo the SSI analysis with it.

> **Source of the details.** ACS SASSI documents LOADGEN fully in a separate manual ("ACS SASSI-ANSYS
> Integration Capability", Options A and AA), which is not available. The dialog fields below are those of
> the V3 manual; the input deck, the file layouts and the APDL are SASSI-EDU's own design. Every
> interpretation is listed in section 10 (decisions D-LGN-01 ... D-LGN-08).

---

## 1. Quick start

You need a finished SSI run: SITE, POINT, HOUSE, ANALYS, then MOTION and RELDISP with the outputs that
LOADGEN reads (section 3). Then, in the same session (or after `INP` of the model):

```
* ---- equivalent static loads at the time of the largest base shear
* Disp. and Accel.; lumped masses generated from the HOUSE matrices
LOADGEN,DISPACC,0,0,0,LUMPED,1
* one critical time: largest base shear along the input direction
LGTIME,V,1
* writes <model>.lgn, runs LOADGEN -> <model>_LGS.inp
RUNLOADGEN,STATIC

* ---- dynamic loads for a transient analysis
* Rayleigh damping 5 % at 1 Hz and 20 Hz
LOADGENDYN,,,REL,0,RESULTS,0,0,0.05,1,20
* -> <model>_LGD.inp
RUNLOADGEN,DYNAMIC
```

In ANSYS (the structural model first, for example the file the SASSI-EDU `ANSYS` command writes):

```
/INPUT,mymodel,inp        ! the structure (nodes, elements, materials, masses)
/INPUT,mymodel_LGS,inp    ! the equivalent static loads; solves one load step
/POST1
SET,LAST
PRRSOL,F                  ! reaction forces; TOTAL VALUES = - sum of the inertia forces (the SSI base shear)
```

or `/INPUT,mymodel_LGD,inp` for the transient analysis. Every LOADGEN file starts with a header that
says what it contains, in which units and axes, and how to run it.

LOADGEN is fast (well under a second for models with a few thousand result histories). Its listing is
`<model>_LOADGEN.out` in the model directory. In the GUI, use **Modules ▸ ANSYS Eq. Static Load** or
**Modules ▸ ANSYS Dynamic Load**: the dialogs hold the fields of the commands below, **OK** emits the
changed LOADGEN/LGFILE/LGNODE/LGTIME/LGMAP/LGOPT commands and **Run** also runs `RUNLOADGEN` as a job with
its listing in an output tab (see `docs/user/GUI.md`). The commands can also be typed in the Command Entry.

---

## 2. Why it works (the theory in one page)

Split the total displacement of the structure into the motion of a reference frame and the motion
relative to it:

```
u(t) = u_r(t) + iota u_g(t)
```

with `iota u_g` a rigid-body translation (the free-field ground motion, or the motion of a reference
node). A rigid-body translation strains nothing (`K iota = 0`), so the equations of the structural
degrees of freedom `s`, with the interface degrees of freedom `b` prescribed, are

```
M_ss u_r,s'' + C u_r,s' + K_ss u_r,s = -M_ss iota_s a_g(t) - K_sb u_r,b(t)
```

* **Dynamic second step.** This is an ANSYS transient analysis with `ACEL = a_g(t)` (ANSYS applies
  `-M ACEL` to every mass: ACEL is the acceleration of the reference frame) and the interface
  displacements relative to the ground motion prescribed with `D` (RELDISP's `.THD` with the free-field
  reference). This is how the manual describes the ANSYS Dynamic Load: "the ANSYS dynamic load generator
  uses all the relative displacement data files with respect to the free-field motion" plus a ground
  acceleration file. ANSYS then computes `u_r`; stresses are those of the total motion.
* **Equivalent static second step.** Freeze the equation at a critical time `t*` and drop the damping and
  the velocity terms: `K_ss u_s = -M_ss a_s(t*) - K_sb u_b(t*)`. The loads are the nodal **inertia forces
  of the absolute accelerations**, `F = -m a(t*)` (MOTION's `.ACC`), and the interface displacements at
  `t*` (RELDISP's `.THD`). Relative and total interface displacements differ by a rigid translation,
  which produces no stress. The dynamic equilibrium of the free body above the interface says that the
  sum of the inertia forces is minus the SSI base shear: VP-LA1 checks it against STRESS.

What is approximate:

* the damping of the dynamic step: SASSI uses frequency-independent hysteretic damping (the complex
  modulus), ANSYS direct integration uses Rayleigh damping `C = alpha M + beta K`, whose damping ratio
  `zeta(w) = alpha/(2 w) + beta w/2` is exact at two frequencies only. The manual warns that this is a
  significant limitation of the dynamic second step;
* the equivalent static step leaves out the damping forces and represents one instant of the response;
* anything in the ANSYS model that changes the interface motion (Option A assumption).

---

## 3. What the SSI run must provide

LOADGEN reads the SSI results directory (the model directory by default):

| needed for | file | how to get it |
|---|---|---|
| everything | `<model>.N4` (FILE4), `<model>.hou` | HOUSE, AFWRITE |
| generated masses | COOSM | HOUSE |
| inertia forces (static), check tables, method ACC | MOTION `.ACC` of the mass / master / check / reference nodes, every direction that carries mass | MOTION: NOUT flag 2 (save history), or MOTIONX `<saveacc>` (Save ACC in All Points); rotations with MOTIONX `<saverot>` ("Save Rotation for Ansys") or NOUT on the rotational DOFs |
| interface displacements | RELDISP `.THD` of the interface nodes, free-field reference (RELFILE blank), or MOTION's complex `.TFI` | RELDISP with RDND on the interface nodes; MOTION must save their complex TFs (NOUT flag 1 or MOTIONX `<savetf>`, `<cplx>` = 1) |
| ACEL of the REL method | the control motion | the THFILE of the RELDISP (or MOTION) deck, or LGFILE,GROUND |

If a per-node file is missing, LOADGEN looks for the MOTION / RELDISP **frames** (`ACC/`, `ACCR/`,
`THD/`, `THDR/`), which is what the ACS SASSI dialog calls "frame data". A relative displacement without
`.THD` or frame is computed from MOTION's **complex `.TFI`** of the node with RELDISP's free-field formula
(the same numbers as RELDISP's `.THD`), so MOTION with complex TFs saved is enough. If nothing is found
the listing names the missing node/DOFs.

**Source FILE8** (`LOADGEN,...,FILE8` or `LOADGENDYN,...,FILE8`) skips all of this: LOADGEN computes the
same histories itself from FILE8 and the control motion, with the MOTION deck's settings (interpolation,
smoothing, phase adjustment, Hudson-Housner baseline) and RELDISP's formula `u = irfft((H - H_ref) U_g)`.
It only needs FILE8 and the MOTION deck (`<model>.mot`, written by AFWRITE). The two sources agree to
1e-9 (VP-LA1, VP-LA2). Use FILE8 for large models: no thousands of history files.

**Units.** SASSI-EDU uses the model's consistent units throughout. MOTION's accelerations are in g;
LOADGEN multiplies them by the HOUSE gravity, so the ANSYS file holds length/s², forces in the model's
force unit (mass = weight / g) and displacements in the model's length unit.

---

## 4. The commands (the two dialogs)

The manual has no LOADGEN command ("Interactive only"): two Modules-menu dialogs write the LOADGEN input
and start it. SASSI-EDU keeps the dialog fields in option records so that sessions replay from `.pre`
files and WRITE/INP keep them (they are written under "Other stored options"). `RUNLOADGEN` is the
dialog's Ok.

### 4.1 ANSYS Static Load Converter: `LOADGEN`

```
LOADGEN,<data>,<multi>,<rotdisp>,<rotacc>,<masstype>,<genmass>,<source>
```

| field | dialog control | values (default first) |
|---|---|---|
| `<data>` | *Data to Add From ACS SASSI* | `2`/`ACC` Acceleration, `1`/`DISP` Displacement, `3`/`DISPACC` Disp. and Accel., `4`/`SOILDISP` Disp. for Soil Module |
| `<multi>` | Use Multiple File List Inputs | `0` one APDL file (the first critical time), `1` one file per critical time |
| `<rotdisp>` | Rotational Disp. | `0`, `1` = prescribe the interface rotations too |
| `<rotacc>` | Rotational Accel. | `0`, `1` = rotary inertia moments / master rotations |
| `<masstype>` | Mass Type | `1`/`LUMPED` Lumped Mass, `2`/`MASTER` Master Node Mass |
| `<genmass>` | Generate Mass Data | `0` read the mass file, `1` write (overwrite) it from the HOUSE mass matrix |
| `<source>` | (SASSI-EDU) | `RESULTS` (MOTION/RELDISP files or frames), `FILE8` |

### 4.2 ANSYS Dynamic Load Converter: `LOADGENDYN`

```
LOADGENDYN,<alpha>,<beta>,<method>,<refnode>,<source>,<rotdisp>,<rotacc>,<zeta>,<f1>,<f2>,<gfopt>,<gmult>
```

| field | meaning |
|---|---|
| `<alpha>`, `<beta>` | Rayleigh damping coefficients (dialog "Raleigh Damping Coeff."), `ALPHAD` and `BETAD` |
| `<method>` | `REL` (default, the manual's): ACEL = ground motion, D = interface displacements relative to it. `ACC`: interface fixed, ACEL = absolute acceleration of `<refnode>` (a fixed-base model driven by the SSI foundation motion; valid for a rigid foundation, the listing reports how rigid it is) |
| `<refnode>` | node whose motion is the reference frame; `0` = the control motion (free field). A node is useful for an embedded foundation (the manual's "kinematic SSI acceleration"): the D tables are then relative to that node (source FILE8, or RELDISP run with that node's `.TFI` as RELFILE) |
| `<source>`, `<rotdisp>`, `<rotacc>` | as for LOADGEN; `<rotacc>` adds rotational check tables |
| `<zeta>`, `<f1>`, `<f2>` | Rayleigh damping from a damping ratio `zeta` at two frequencies (Hz): `alpha = 2 zeta w1 w2/(w1 + w2)`, `beta = 2 zeta/(w1 + w2)`; used when `zeta > 0` |
| `<gfopt>`, `<gmult>` | format (`0` dt then values, `1` (t, a) pairs) and factor of the Ground Acceleration File |

### 4.3 Files and paths: `LGFILE`

`LGFILE,<key>,<name>` — the name is the rest of the line (spaces allowed); a blank name restores the
default.

| key | dialog box | default |
|---|---|---|
| `SSIPATH` | SASSI Model and Results Input: Path | model directory |
| `HOUSE` | HOUSE Module Input | `<model>.hou` |
| `DISP`, `DISPROT` | Displacement Results (frames, translations / rotations) | `THD`, `THDR` |
| `ACC`, `ACCROT` | Acceleration Results (frames) | `ACC`, `ACCR` |
| `ANSYSPATH` | ANSYS Model and Data Input: Path (where the outputs go) | model directory |
| `LUMPED`, `MASTER` | Lumped node / Master Node Mass file | `<model>.masl`, `<model>.masm` |
| `APDL`, `APDLDYN` | ANSYS Output File (static / dynamic) | `<model>_LGS.inp`, `<model>_LGD.inp` |
| `GROUND` | Ground Acceleration File (in g) | the control motion of the RELDISP / MOTION deck |

### 4.4 Node lists, critical times, numbering, output options

* `LGNODE,<kind>,<nodes>` adds nodes (model numbers, ranges `a-b`); `LGNODE,<kind>,0` clears the list.
  Kinds: `D` the nodes that receive `D` (default: the interaction nodes; for "Disp. for Soil Module" also
  the nodes of the excavated soil elements), `M` the master nodes for generating a Master Node Mass file,
  `A` the nodes whose SSI absolute accelerations are exported as check tables (dynamic).
* `LGTIME,<crit>,...` chooses the critical times of the static loads (section 5.3).
* `LGMAP,<mode>,<tol>,<file>` the ANSYS node numbering (section 7).
* `LGOPT,<digits>,<opmode>,<rest>`: significant digits of the APDL numbers (default 12), the data-check
  mode (`1`: check everything, write no file) and `<rest>` (`1`, default: every relative-displacement
  history starts at rest, section 6.3; `0`: as computed by MOTION/RELDISP).
* `LGLIST,[STATIC|DYNAMIC]` lists the settings and what RUNLOADGEN will do.
* `RUNLOADGEN,[STATIC|DYNAMIC],[model]` writes `<model>.lgn` (the LOADGEN deck) and runs LOADGEN in the
  model directory (the RUN<MODULE> order `RUNLOADGEN,<model>,[STATIC|DYNAMIC]` works too). Batch:
  `python -m sassi.modules.loadgen < LOADGEN.inp` with the three lines model, deck, listing (manual
  section 3.2).

---

## 5. Equivalent static loads ("ANSYS Eq. Static Load")

### 5.1 What each "Data to add" writes

| data | F (inertia forces) | D | use it for |
|---|---|---|---|
| 2 Acceleration | `F = -m a(t*)` at every mass node except the interface nodes | `D = 0` at the interface nodes (every DOF they have): a **fixed base** | the classical equivalent static analysis of a structure on a stiff (rigid) foundation |
| 3 Disp. and Accel. | the same | the interface displacements relative to the free field at `t*` | the complete quasi-static decomposition: foundation flexibility, rocking, differential motion of the interface |
| 1 Displacement | none | the relative displacements at the `D` nodes | driving a model (or a sub-model) by the SSI deformation only |
| 4 Disp. for Soil Module | none | the relative displacements at the interaction nodes and the nodes of the excavated soil | the boundary displacements of a refined ANSYS soil model (soil pressures) |

The inertia of the interface nodes goes straight to the reactions when their displacements are
prescribed; the listing shows that mass separately. Supports: the interface `D` must restrain every
rigid-body motion of the ANSYS model — with several interface nodes the translations are enough; a single
interface node (a stick) needs its rotations (`<rotdisp>` 1) or a fixity in your model.

### 5.2 Masses

**Lumped Mass** (`.masl`, one line per node: `node mx my mz mxx myy mzz`, model numbers, mass units). With
Generate Mass Data LOADGEN builds it from the HOUSE structure mass matrix `Ms` (COOSM):

```
m_(i,d) = sum_j M[(i,d), (j,d)]       (translations: row sums over the same direction)
I_(i,r) = M[(i,r), (i,r)]             (rotations: diagonal)
```

so lumped MT/MR masses are reproduced exactly and the total mass per direction of a consistent (beam,
solid, shell) mass matrix is preserved. Without COOSM the MT/MR masses of the HOUSE deck are used
(warning). You can also give your own file — for example the nodal masses of your ANSYS model, if its
nodes coincide with the SSI nodes. A line `node m` means `mx = my = mz = m`.

**Master Node Mass** (`.masm`, lines `node master x y z mx my mz [mxx myy mzz]`): each **load node** — an
ANSYS node, numbered as in ANSYS, at `x y z` (SASSI-EDU global axes) with its own masses — moves with the
rigid-body motion of its **master node** of the SSI model:

```
a_k = a_M + alpha_M x (x_k - x_M)          F_k = -m_k a_k,   M_k = -I_k alpha_M
```

This is the load transfer from a coarse SSI model (one master node per floor of a stick) to a refined
ANSYS model (every node of the floor). The rotational accelerations `alpha_M` need `<rotacc>` = 1 (and
the rotational `.ACC` files or source FILE8); without them `alpha_M = 0` (warning). Generate Mass Data
writes a template from the SSI model: every node with mass becomes a load node, its master is the master
node (`LGNODE,M,...`) closest in elevation, then in plan. Replace the load nodes and masses by those of
your ANSYS model.

### 5.3 Critical times

`LGTIME` (default `LGTIME,V,1`):

| form | critical times |
|---|---|
| `LGTIME,V,<n>,<tsep>` | the n largest peaks of the base shear along the input (control) direction |
| `LGTIME,VX\|VY\|VZ,<n>,<tsep>` | of a base-shear component |
| `LGTIME,MX\|MY\|MZ,<n>,<tsep>,<x0>,<y0>,<z0>` | of the overturning / torsional moment about (x0, y0, z0); blank = centroid of the interface nodes |
| `LGTIME,ACC\|DISP,<n>,<tsep>,<node>,<dof>` | of a node's absolute acceleration / relative displacement |
| `LGTIME,TIME,<t1>,<t2>,...` | given times (s; the nearest sample is used, with a warning) |
| `LGTIME,STEP,<k1>,<k2>,...` | given time steps (1-based; the frame numbers of ACC/THD frames) |

V and M are the resultants of the inertia forces `F = -m a` of the loaded nodes, in SASSI-EDU global
axes. "Peaks" are local maxima of the absolute value, at least `tsep` seconds apart. With `<multi>` = 0
only the first (largest) critical time is written — the manual's single seismic load file; with
`<multi>` = 1 every critical time gets its own file `<apdl>_01.inp`, `<apdl>_02.inp`, ...

The resultant histories are also written to `<apdl>_res.txt` (columns `t VX VY VZ MX MY MZ`), handy to
plot the base shear and to see what the critical times are.

### 5.4 An equivalent static file

```
! SASSI-EDU 0.1.0  LOADGEN (Option A): ANSYS equivalent static seismic load, file 1 of 3
! Load step  : SSI time t = 6.8 s (sample 341 of 1024, dt 0.02 s); critical time 1 by base shear
!              along the control direction (sum of the inertia forces)
! Data       : Disp. and Accel. ...
! Units      : the SSI model's consistent units, SI-like (m, s): g = 9.81 length/s^2 ...
LG_SOLVE = 1                ! 1 = solve at the end of this file; 0 = apply the loads only
LG_TSSI = 6.8               ! SSI time of this load step (s)
/SOLU
ANTYPE,STATIC
TIME,1                      ! load step 1 <- SSI time 6.8 s
! ---- inertia forces F = -m a(t) (9 components)
F,10,FX,-9.40897475678E+00
F,11,FX,-1.81924664364E+01
F,12,FX,-2.16485514706E+01
...
! ---- SSI displacements relative to the free field (27 components)
D,1,UX,-1.66704335837E-06
D,1,UZ,-4.07154261628E-06
...
*IF,LG_SOLVE,EQ,1,THEN
SOLVE
*ENDIF
FINISH
```

Components that are numerically zero (for example `FY` of an X-input) are written too: they are part of
the SSI solution and do no harm. Several files of one run can be read one after the other in one ANSYS
session: each re-applies all its `F` and `D` (same nodes, new values) and solves a new load step
(`TIME,1`, `TIME,2`, ...). Dead load is not included: combine with a gravity load case (`ACEL,0,0,g` for
Z up) as your design code requires.

**Check in ANSYS:** `PRRSOL,F` in POST1 after the solution: the TOTAL VALUES of the reactions are minus
the sum of the `F` commands (the header's "Resultant"), which is the SSI base shear of the structure above
the interface (the STRESS base element force) at that time.

---

## 6. Dynamic loads ("ANSYS Dynamic Load")

### 6.1 The file

```
LG_SOLVE = 1                ! 1 = SOLVE at the end of this file, 0 = define the loads only
LG_DT = 0.02                ! time step of the SSI histories (s)
LG_NT = 1024                ! samples per table
LG_TEND = 20.46             ! end time (s) = (LG_NT - 1) LG_DT
! ---- TABLE arrays: column 0 = time (s), column 1 = value
*DIM,LG_ACX,TABLE,1024,1,1,TIME   ! reference-frame acceleration X (length/s^2)
*VFILL,LG_ACX(1,0),RAMP,0,LG_DT
LG_ACX(1,1)=0,-1.18632218893E-04,-7.85195784617E-05, ... (10 values per line)
...
*DIM,LGD_1_UX,TABLE,1024,1,1,TIME   ! node 1 UX relative to the reference motion
*VFILL,LGD_1_UX(1,0),RAMP,0,LG_DT
LGD_1_UX(1,1)=...
...
/SOLU
ANTYPE,TRANS
TRNOPT,FULL
TIMINT,ON
AUTOTS,OFF
KBC,0
DELTIM,LG_DT                ! constant time step = SSI time step
TIME,LG_TEND
ALPHAD,0.5711986642890534
BETAD,0.0014468631190172306
OUTRES,ALL,ALL              ! every substep (reduce for large models, e.g. OUTRES,ALL,10)
ACEL,%LG_ACX%,0,0           ! acceleration of the reference frame = ground motion
D,1,UX,%LGD_1_UX%
D,1,UY,%LGD_1_UY%
...
*IF,LG_SOLVE,EQ,1,THEN
SOLVE
*ENDIF
FINISH
```

* One 1-column TABLE per history, primary variable TIME; the time column is filled by `*VFILL ... RAMP`,
  the values by assignments of at most 10 values. ANSYS interpolates the tables at every substep.
  Names: `LG_ACX/Y/Z` (ACEL), `LGD_<node>_<label>` (D), `LGA_<node>_<label>` (check tables).
* One load step, constant time step = the SSI time step, full transient (Newmark). The structure starts
  at rest (zero initial conditions).
* Set `LG_SOLVE = 0` to add your own options (NLGEOM, contact, OUTRES, TINTP ...) before `SOLVE`.
* ANSYS counts the parameters: with several thousand interface DOFs the number of TABLE parameters may
  exceed what your release allows (LOADGEN warns above 4000).

### 6.2 Reading the results

With method REL the ANSYS displacements, velocities and accelerations are **relative to the reference
frame**; stresses and internal forces are the actual ones. The absolute acceleration of a node is the
ANSYS acceleration plus the ACEL table. The check tables (`LGNODE,A,...`) hold the SSI absolute
accelerations of chosen nodes, so that you can compare: with the same damping they would coincide
(VP-LA2 replays the second step and reproduces the SSI accelerations to 5e-9); the differences you see in
ANSYS come from Rayleigh damping, the integration scheme and your model changes.

### 6.3 Rayleigh damping and other practical points

* Choose `zeta` (`LOADGENDYN,...,<zeta>,<f1>,<f2>`) equal to the structural damping of the SSI model, `f1`
  near the fundamental SSI frequency and `f2` near the highest frequency that matters (20-33 Hz for
  in-structure response). Between f1 and f2 the damping is lower than `zeta`, outside higher; the
  listing prints the damping ratio at 1, 2, 5, 10 and 20 Hz. Rayleigh damping cannot reproduce the
  radiation damping of the soil: that is in the prescribed interface motion, which is why the D tables
  matter.
* **Start at rest.** MOTION and RELDISP work with the FFT, whose zero-frequency displacement is
  undetermined and set to zero: their relative displacements are the physical ones minus their mean over
  the Fourier period, so they do not start at zero. LOADGEN subtracts the initial value of every relative
  displacement history (the listing gives the size of that offset); this is exact when the quiet zone
  (NFFT dt minus the record length) lets the response decay. The static D values get the same
  correction. Accelerations are not changed.
* ANSYS's default Newmark parameters add a little numerical damping (TINTP); keep the time step equal to
  the SSI time step (the histories have no information between samples).
* **Gravity** is not in the file. ANSYS keeps one ACEL per load step, so with `LG_SOLVE = 0` put it in the
  Z field of the ACEL line (Z up; for example a TABLE holding the constant g, defined like the others) and
  solve a static gravity step first (`TIMINT,OFF`, a short `TIME`, `SOLVE`, then `TIMINT,ON` and the
  transient) so that the dead load is not applied suddenly.
* The mass of the ANSYS model is what ACEL loads: the `ANSYS` command exports MASS21 elements and DENS
  from the SSI model; SASSI uses ½ lumped + ½ consistent solid mass and lumped shell mass, ANSYS consistent
  mass unless you set `LUMPM,ON` (see ANSYS.md §3).

---

## 7. Node numbering

Three numberings meet in Option A:

1. **FILE4 / FILE8 numbering** — the numbers in the MOTION and RELDISP file names. They are the model
   numbers unless the HOUSE node optimizer renumbered the model (HOUSEX `<optimize>` = 1). LOADGEN
   translates them back automatically (FILE4 `x_node_old_id`), so you never see them.
2. **Model numbering** — the numbers of your `.pre` / `.hou`. All LOADGEN inputs (LGNODE, LGTIME nodes,
   `<refnode>`, lumped-mass files, masters) use them.
3. **ANSYS numbering** (`LGMAP`):
   * `LGMAP,IDENTITY` (default): ANSYS node = model node. True when the ANSYS model was written by the
     SASSI-EDU `ANSYS` command or the SSI model came from `CONVERT,ANSYS` (both keep the node numbers).
   * `LGMAP,PAIRS,,<file>`: lines `sassi_node ansys_node`.
   * `LGMAP,COORD,<tol>,<file.cdb or APDL>`: each SSI node takes the ANSYS node at the same position
     (within `tol`, default 1e-6 of the model size) of the ANSYS file (`CDWRITE` output or `N` commands);
     the map found is written to `<model>_LG.map`. Nodes of a refined ANSYS mesh that lie between SSI
     nodes get no load: use a Master Node Mass file for the inertia forces of a refined mesh.

A load or interface node without an ANSYS node stops the run with the list of those nodes. LOADGEN
writes global components; nodes with rotated nodal coordinate systems in ANSYS (`NROTAT`) would take them
in their nodal axes (COORD warns when the `.cdb` has such nodes).

**2-D models** (PLANE groups) are rotated exactly as the `ANSYS` command does: the SASSI X-Z plane becomes
the ANSYS X-Y plane (`UZ -> UY`, `FZ -> FY`, `ROTY -> -ROTZ`).

---

## 8. Step-by-step example (a stick on a rigid mat)

The model of VP-LA1: a three-mass stick (2, 2 and 1.5 t at 4, 8 and 12 m) on a 4 m × 4 m rigid mat of 9
interaction nodes, uniform site Vs = 1000 m/s, earthquake-like record (PGA 0.3 g).

```
* ... model, SITE, POINT, HOUSE, ANALYS as usual, then
MOTION,0,0,0,0,0,0.1,100,3,1,0,1,0,0,0,0,1,0,0,1
THFILE,eq.acc
* TF + history of the stick nodes, X ...
NOUT,1,1,1,0,0,0,1,10-12
* ... Y ...
NOUT,2,1,1,0,0,0,1,10-12
* ... and Z (the lumped masses act in all directions)
NOUT,3,1,1,0,0,0,1,10-12
* Save TF in All Points (RELDISP of the mat nodes needs them)
MOTIONX,0,2,0,1
* RELFILE blank = free-field reference
RELD,1,0,0
* relative displacements X, Y, Z of mat node 1 ...
RDND,1,1,1,1,0,0,0
...                                 ! ... one RDND line per mat node, up to RDND,9,1,1,1,0,0,0
AFWRITE
RUNSITE
RUNPOINT
RUNHOUSE
RUNANALYS
RUNMOTION
RUNRELDISP
* Disp. and Accel., generate the lumped masses
LOADGEN,DISPACC,0,0,0,LUMPED,1
LGTIME,V,1
RUNLOADGEN,STATIC
* the structural model in ANSYS format (ANSYS.md)
ANSYS
```

(VP-LA1 builds this model with `sassi.verify.builders`; the lines above are the equivalent commands.)
LOADGEN finds the largest base shear at t = 6.8 s; the inertia forces sum to −49.11 kN, the STRESS base
shear of the stick is 49.13 kN (VP-LA1: 0.04 % apart). In ANSYS, after `/INPUT,model,inp` and
`/INPUT,model_LGS,inp`, `PRRSOL,F` in POST1 must show a total base reaction of +49.11 kN in X: by
equilibrium the reactions balance the applied forces exactly, whatever the ANSYS mesh.

---

## 9. Files

| file | written when | content |
|---|---|---|
| `<model>.lgn` | RUNLOADGEN | the LOADGEN deck (SASSI-EDU deck format, `SASSI-EDU LOADGEN DECK v1`) |
| `<model>_LOADGEN.out` | every run | listing: options, inputs, masses, critical times, loads, files, warnings |
| `<model>_LGS.inp` (or `_01`, `_02`...) | static | equivalent static load step(s) |
| `<model>_LGS_res.txt` | static with masses | resultant histories `t VX VY VZ MX MY MZ` of the inertia forces |
| `<model>_LGD.inp` | dynamic | tables, ACEL, D and the transient solution |
| `<model>.masl` / `.masm` | Generate Mass Data | lumped / master-node mass file |
| `<model>_LG.map` | LGMAP,COORD | the node map found (`sassi ansys` pairs) |

---

## 10. SASSI-EDU interpretations (decisions D-LGN-01 ... D-LGN-08)

| ID | subject | decision |
|---|---|---|
| D-LGN-01 | input | the dialogs are option records (LOADGEN, LOADGENDYN, LGFILE, LGNODE, LGTIME, LGMAP, LGOPT); RUNLOADGEN writes the LOADGEN deck `<model>.lgn` and runs the module (the dialog's Ok). The manual's "no command" rule is kept in spirit: the records only describe the dialog |
| D-LGN-02 | static data types | Acceleration = inertia forces + fixed interface; Disp. and Accel. = inertia forces + relative interface displacements; Displacement and Disp. for Soil Module = displacements only (default nodes: interaction nodes; soil module: + excavated-soil nodes) |
| D-LGN-03 | masses | Lumped = row sums per direction of the HOUSE `Ms` (diagonal rotary); Master Node Mass = load nodes moving rigidly with a master node; generated master = closest in elevation, then in plan |
| D-LGN-04 | dynamic | REL = ACEL of the reference motion + D of the interface displacements relative to it (manual); ACC = fixed interface + ACEL of a node's absolute acceleration; Rayleigh damping; ACEL is the reference-frame acceleration (`-M ACEL` on the masses) |
| D-LGN-05 | numbering | FILE4 (optimizer) -> model -> ANSYS by identity, pairs file or coordinates |
| D-LGN-06 | histories | RESULTS = MOTION `.ACC` / RELDISP `.THD`, else frames, else (displacements) MOTION's complex `.TFI` with RELDISP's formula; FILE8 = MOTION's and RELDISP's algorithms; relative displacements started at rest |
| D-LGN-07 | critical times | largest local maxima of the absolute criterion, at least `tsep` apart; or given times / steps; one file unless Use Multiple File List Inputs |
| D-LGN-08 | APDL | model units, global axes (2-D rotation of the ANSYS export), 12 significant digits, ≤ 10 values per assignment, names ≤ 32 characters, header with units/axes/instructions, `LG_SOLVE` switch, nothing of the user model redefined except the interface supports of the fixed-base variants |

---

## 11. Verification

* **VP-LA1** (`VERIFY,VP-LA1`; `pytest tests/verification/test_vp_loadgen.py`): the stick of section 8
  through the whole chain. Generated masses = MT masses (1e-12); critical time = time of the largest STRESS
  base shear (±1 sample); sum of the APDL `FX` = − STRESS base shear `FYI` and the moment about the base =
  − STRESS base moment `MZI` at that time, both within 1 % (observed 0.037 % and 0.032 %: MOTION
  interpolates nodal TFs, STRESS the element force TF); interface `D` = RELDISP `.THD` started at rest
  (2e-12); source FILE8 = source RESULTS (exact).
* **VP-LA2** (`VERIFY,VP-LA2`): the same stick on a mat on soil springs, FILE8 exact at every Fourier
  frequency. The APDL tables reproduce MOTION `.ACC` (ACC method's ACEL and the check tables), the control
  motion (REL method's ACEL) and RELDISP `.THD` to 1e-10 (observed ≤ 5e-12); FILE8 and RESULTS sources
  agree to 1e-9; the APDL parses (`*DIM`, assignments, `*VFILL`, `*IF`, D, F, ACEL, solution commands).
  Finally the second step is **replayed**: the structure with the exported D tables and ACEL, solved with
  the HOUSE matrices at every frequency, reproduces the SSI absolute accelerations of the stick to 5e-9
  (criterion 1e-6) — the exactness of the decomposition of section 2.
* Unit tests: `tests/unit/test_loadgen.py` (deck, axes, numbering, masses, peaks, APDL, every data type,
  mass type, source, method and numbering) and `tests/unit/test_loadgen_cmds.py` (commands, WRITE/INP,
  RUNLOADGEN).

---

## 12. Limitations

* The exact ACS SASSI LOADGEN formats are unknown: files from the commercial program are not read and the
  APDL differs from the commercial one (same purpose, documented here).
* No interpolation of the SSI motion to ANSYS nodes between SSI nodes: interface nodes must coincide (use
  a refined SSI interface, or the Master Node Mass transfer for inertia).
* Dynamic: one control direction per run (as the SSI run); the large-mass method for absolute
  accelerations at many supports is not generated (the D tables of relative displacements are exact and
  preferable).
* One APDL TABLE per interface DOF: very large interfaces may exceed the parameter limit of your ANSYS
  release.
