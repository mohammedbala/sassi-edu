# 10 - Command Reference: Module, Calculation, Plot and Programming Commands (Manual 9.12-9.15)

Source: ACS SASSI V3 User Manual, Chapter 9, Sections 9.12 (Module Commands), 9.13 (Calculation
Commands), 9.14 (Plot Commands) and 9.15 (Programming Commands). Text lines 11784-13260 of
`reference/acs-sassi.txt`; printed pages 282-312 (PDF pages 284-314).

PDF pages viewed as images and transcribed: 285-305 (equations of SHEAR, page 291-292; SHOWDOF
"Boundary Conditions" window, page 304). For dialogs these commands open, the Chapter 7 screenshots
were also viewed: "Line Selection" and "Graph Plot Options" (PDF 200, 202), "Parse Frame Data" and
"Load Frame Data" (PDF 209-210), "Window Options" and "Shader Options" (PDF 194-195), "Select Cut
to Display" (PDF 198), "Soil Layer Windows Setting", "Select Dynamic Soil Property", "Soil
Properties Window Settings" (PDF 206-208), Tables 3.1/3.2 (PDF 49-50).

Cross-referenced manual text outside this range (only to make these commands implementable):
5.6 Macros, 5.8 Section cuts, 5.9 Variables and loops, 3.x STRESS/MOTION outputs, 6.4 Modules menu,
6.5.6 Colors, 6.5.8 Shader options, 7.1-7.3 Plotting, 9.3 CSYS/LOC/LOCAL, 9.7.16 GRAVITY, 9.7.20 MDL,
9.10 CSECT/cuts, 9.17.3 BBCGEN, 9.17.23 P.

Companion spec: `04_model_ui_modules.md` (workflow layer for macros, variables, cuts, modules).
This file is the **authoritative per-command spec** for 9.12-9.15. Where the two differ, see
Open Questions Q-30.

## Tag legend

- **[CORE]** needed for computational correctness
- **[IO]** file or input format
- **[UI]** user interface, plotting, convenience
- **[ADV]** advanced option (incoherency, nonlinear Option NON, Option A/AA/PRO)
- **(derived)** the rule is inferred by this spec. The manual does not state it. Treat it as an
  implementer decision and verify it.
- **(ext. ref.)** the rule comes from an external standard, not from the manual. Verify it.

## Notation

- Backus-Naur style as in the manual: `<x>` required, `[x]` optional. Arguments are
  **comma-separated**. A blank argument (nothing, or only whitespace, between two commas) means
  "use the default".
- Command names are **not case sensitive**. None of the commands in 9.12-9.15 shows an underlined
  abbreviation in the PDF, so they must be typed **in full** (the PREP 4-letter rule does not
  apply to these UI-only commands).
- "Active model" means the model selected with `ACTM`. A `[model]` argument of `-1` means the
  active model.
- Lines starting with `*` in a `.pre` file are comments (manual examples).

---

## 0. Command index for this file

| Group | Commands | Manual section |
|---|---|---|
| Module run | RUNANALYS, RUNCOMBIN, RUNEQUAKE, RUNFORCE, RUNHOUSE, RUNMOTION, RUNPOINT, RUNRELDISP, RUNSOIL, RUNSTRESS | 9.12.1-9.12.10 |
| Calculation | CALCC, CALCSECTHIST, CALCSECTHISTDB, CALCM, CALCMOI, CALCPAR, READSTR, SECDATAOPT, SHEAR | 9.13.1-9.13.9 |
| Plot (52) | ADDITION, AVERAGE, AXES, BROADEN, BUBBLEPLOT, CAPTUREPLOT, CLOSEPLOT, CNGCENTER, CNGVIEW, COLOR, CONTOURPLOT, CUTPLOT, DEBUG, DEFORMPLOT, ELECOLOR, ELENUM, GROUPNUM, LBINCORS, LAYERPLOT, LINECOMBIN, LINENAME, MARKERS, MODELPLOT, NODENUM, NODEPLOT, NODESEL, PAUSE, PLOTRANGE, PLOTTITLE, PROCFRAME, READSPEC, READTH, RSTCENTER, RSTVIEW, SHADEROPTIONS, SHOWDOF, SHOWMASS, SHRINK, SOILPROPPLOT, SPECPLOT, SRSS, STIPPLE, SUBTRACTION, THPLOT, VECTORPLOT, WINDOWSETTINGS, WIREFRAME, WRITESPEC, WRITETH, XTITLE, YTITLE, YTITLE2 | 9.14.1-9.14.52 |
| Programming | ADDRND, CD, FOREACH, LOADMACRO, LOADVAR, MACRO, MACROLIST, MKDIR, REDUCESET, RND, RNDSEED, SETVAR, SHOWVAR, VAR, VARLIST | 9.15.1-9.15.15 |

Manual numbering quirks (reproduce the commands, not the typos): the 9.13 summary table lists
CALCM as 9.13.5, CALCMOI 9.13.6, CALCPAR 9.13.7, READSTR 9.13.8, and both SECDATAOPT and SHEAR as
9.13.9. The actual headings are CALCM 9.13.4, CALCMOI 9.13.5, CALCPAR 9.13.6, READSTR 9.13.7,
SECDATAOPT 9.13.8, SHEAR 9.13.9. The table calls CALCC "center of area". The body text says
"volume centroid".

### 0.1 Quick map for ANSYS users [UI]

The closest ANSYS APDL equivalents, given only as orientation. Behaviour is **not** identical.

| ACS SASSI UI | Closest ANSYS idea |
|---|---|
| RUNxxx | running a solver step / batch job |
| CALCPAR, CALCSECTHIST(DB) | section forces: `SPOINT` + `FSUM` on a selected cut, or SURF surface operations (`SUCR`/`SUMAP`/`SUEVAL`); Mechanical "Surface" construction geometry + force/moment reaction probes |
| CALCM, CALCC, CALCMOI | mass and centroid from `*GET`/`PRECISION` mass summary; section properties |
| READSPEC/READTH + ADDITION/LINECOMBIN/SRSS | POST26 variables and `ADD`/`PROD`/`SQRT` operations; `*VREAD`/`*VWRITE` |
| BROADEN | no direct APDL command (user spectrum smoothing) |
| AXES, PLOTRANGE, XTITLE/YTITLE | `/GRID`, `/XRANGE`/`/YRANGE`, `/AXLAB` |
| CAPTUREPLOT | `/SHOW,PNG` |
| VAR / FOREACH / MACRO / LOADMACRO | `*DIM`/`*SET`, `*DO...*ENDDO`, `*USE`, `*CREATE` |
| CD, MKDIR | `/CWD`, `/MKDIR` |

---

## 1. Shared definitions used by many commands

### 1.1 Model reference numbers and the active model [CORE][UI]

- Models in memory are addressed by integer reference numbers (see ACTM, MODELLIST in Section 9.7).
- A command argument `[model]` or `<model>` with value `-1` (the default) means the active model.
- A model needs a **name and path** (MDL command, or the model database) before any module can be
  run on it, because module files are named `<modelname>.<ext>` in `<path>`.

### 1.2 Line objects [CORE][UI]

The manual calls plotted numeric data **line objects**:

- They are loaded by commands such as READSPEC and READTH, or created by line-math commands
  (ADDITION, AVERAGE, BROADEN, LINECOMBIN, SRSS, SUBTRACTION, LBINCORS).
- They are stored **for the session only** and cleared when the UI closes.
- They are addressed by an integer **reference number** chosen by the user. Writing to a number that
  already holds a line **overwrites** it.

Data model (derived):

```
LineObject:
    number : int              # reference number (key in the session line table)
    name   : str              # legend text; default set by the creating command
    x      : float[n]         # abscissa: frequency (Hz) for spectra/TF, time (s) for histories
    y      : float[n]         # ordinate
    markers: bool = False     # global property (MARKERS)
    color  : palette entry    # SpecLines / THLines palette by line number (COLOR)
```

- `name` and `markers` are **global line properties**. Changing them affects every graph that
  shows the line, not only the active graph.
- Default names given by creating commands (exact strings, including the period):

| Creating command | Default name |
|---|---|
| ADDITION, SUBTRACTION, LINECOMBIN | `Linear Combin.` |
| AVERAGE | `Average Line` |
| BROADEN | `Envelope` |
| SRSS | `SRSS Line` |
| READSPEC, READTH | file name (derived from the Line Selection screenshot, which lists `1: Test2.rs`) |
| LBINCORS | not stated (Open Questions) |

### 1.3 Union-grid resampling (common to all multi-line operations) [CORE]

ADDITION, AVERAGE, BROADEN, LINECOMBIN, SRSS, SUBTRACTION and WRITESPEC all use the same rules,
which the manual states identically for each command:

1. The result has **the X values of all the source line objects** (union of abscissas).
2. **Linear interpolation** is used to get a line's value between its own data points.
3. **Extrapolation** outside a line's range uses that line's **first or last defined y-value**
   (constant extrapolation, not zero, not linear).

Implementation (derived where noted):

```
def union_grid(lines):
    xs = sorted(set(x for L in lines for x in L.x))
    return merge_close(xs, rtol=1e-9)           # (derived) merge values equal within tolerance

def value_on(L, xq):                             # piecewise linear, constant extrapolation
    if xq <= L.x[0]:  return L.y[0]
    if xq >= L.x[-1]: return L.y[-1]
    i = index with L.x[i] <= xq < L.x[i+1]
    return L.y[i] + (L.y[i+1]-L.y[i]) * (xq-L.x[i]) / (L.x[i+1]-L.x[i])
```

- Interpolation is linear in **linear x and linear y**, even when the plot axes are logarithmic
  (derived; the manual says only "linear interpolation").
- Each source line must have strictly increasing x. If a file gives duplicate x values, keep the
  first occurrence and warn (derived).
- WARNING (derived): constant extrapolation of a time history past its last point repeats the last
  sample value. For acceleration histories of different durations, that is not the same as
  zero-padding. Tell the user this in the command history when the source ranges differ.

### 1.4 Cut-plane local coordinate system (CALCSECTHIST, CALCSECTHISTDB, CALCMOI, CALCPAR) [CORE]

Inputs: plane normal `n = (normalx, normaly, normalz)` in global coordinates. "Right" vector
`r = (rightx, righty, rightz)` in global coordinates, which the manual defines as the direction of
the **local system's positive X**. System number `sysno`.

Construction (derived, consistent with "r defines local +X"):

```
ez = n / |n|                                   # error if |n| = 0
r' = r - (r . ez) ez                           # project r onto the cut plane (Gram-Schmidt)
if |r'| < 1e-8 * |r|: error "right vector is parallel to the plane normal"
ex = r' / |r'|                                 # local X (in the plane)
ey = ez x ex                                   # local Y (in the plane), right-handed: ex x ey = ez
origin = area centroid C of the section        # manual: centroid is the origin of the local system
```

- Store the system as Cartesian local coordinate system number `sysno` in the same table used by
  `LOC`/`LOCAL`/`CSYS` (Section 9.3). Do **not** activate it (derived: storing must not change the
  active coordinate system for later `N` commands).
- `sysno = 0` is the global system and must not be overwritten. If `sysno` is 0 or blank, compute
  everything but do not store the system, and print an informational message (derived; the manual
  example `calcpar,1,0,0,0,1,0` leaves sysno out).
- Local components reported (naming and order are derived; see 3.6.3):
  - in-plane shear forces `Fx` (along ex) and `Fy` (along ey),
  - normal (axial) force `Fz` (along ez, positive = tension),
  - bending moments `Mx` (about ex) and `My` (about ey),
  - torsion `Mz` (about ez).

### 1.5 Unit system detection from gravity [CORE]

The SHEAR command (and by extension CALCM) uses the model's **gravity constant** (`GRAVITY`
command, or the value set by HOUSE options) to decide the unit system of the FE model:

| Unit system | Geometry | Force | Stress | Typical g |
|---|---|---|---|---|
| British | ft | kips | ksi (SHEAR inputs); ksf elsewhere in the model | 32.2 ft/s^2 |
| International | m | kN | kN/m^2 | 9.81 m/s^2 |

Rule (derived): British if `g > 20`, otherwise International. Warn if `g` is not within 5% of
32.174 or 9.80665.

### 1.6 Command history messages [UI]

Every command writes to the Command History, using the message types of manual 6.6.3: echo,
acceptance, information, warning, error, comment. Errors stop only the failing command. A `.pre`
or macro keeps executing with the next line (derived; consistent with batch use).

---

## 2. Module commands (9.12) [CORE][IO][UI]

### 2.1 Common syntax and behaviour

| Command | Syntax | Runs module |
|---|---|---|
| RUNANALYS | `RUNANALYS,[model]` | ANALYS |
| RUNCOMBIN | `RUNCOMBIN,[model]` | COMBIN |
| RUNEQUAKE | `RUNEQUAKE,[model]` | EQUAKE |
| RUNFORCE | `RUNFORCE,[model]` | FORCE |
| RUNHOUSE | `RUNHOUSE,[model]` | HOUSE |
| RUNMOTION | `RUNMOTION,[model]` | MOTION |
| RUNPOINT | `RUNPOINT,[model]` | POINT (POINT2 / POINT3) |
| RUNRELDISP | `RUNRELDISP,[model]` | RELDISP |
| RUNSOIL | `RUNSOIL,[model]` | SOIL |
| RUNSTRESS | `RUNSTRESS,[model]` | STRESS |

Argument: `<model>` is the user selection of the model the module uses. Default `-1` = active model.

Manual rules (all ten commands are identical apart from the module name):

1. The model must be in UI memory.
2. The user must first run **AFWRITE** to write the module input file. AFWRITE runs CHECK first and
   skips any analysis file that has errors.
3. The model name and path must be defined by **MDL** or by the model database.
4. With the default argument, the module runs on the active model.

Execution algorithm (derived from 6.4 and Chapter 3; see also spec 04, section 15.3):

```
RUN<MOD>(model=-1):
    m = active_model if model == -1 else models[model]       # error if not in memory
    require m.name and m.path                                  # else error "use MDL or database"
    inp = m.path / (m.name + ext_in[MOD])                      # ext table: Modules > Extension
    require exists(inp)                                        # else error "run AFWRITE first"
    require predecessor files (table 2.2)                      # else error naming the module to run
    out = m.path / (m.name + ext_out[MOD])                     # e.g. modelname_house.out
    exe = module_location[MOD]                                 # Modules > Location (SASSIini.xml)
    cwd = m.path                                               # all generated files go to the model dir
    open a new UI tab streaming the module console output
    run synchronously; feed stdin lines: m.name, inp.name, out.name   # batch protocol of the modules
    report exit status; on failure show the tail of the output file
```

- Run synchronously while a `.pre` file or macro is executing, so that sequences such as
  `RUNHOUSE` then `RUNANALYS` work in batch (derived).
- The input-file extensions are user-configurable (Modules > Extension, stored in
  `SASSIini.xml`). Defaults: `.equ .soi .sit .poi .hou .frc .anl .mot .str .rdi` (spec 04,
  section 15.2).
- Python implementation (decision): `module_location[MOD]` may be an external executable (for
  cross-checking against the commercial code) or the built-in Python module entry point. Default to
  the built-in one.

### 2.2 Module prerequisites (needed to give useful errors) [CORE]

From manual Chapter 3 and 6.4 (details in specs 02 and 04):

| Module | Input file (AFWRITE) | Also requires | Produces (main) |
|---|---|---|---|
| EQUAKE | `.equ` | - | accelerograms compatible with target spectra |
| SOIL | `.soi` | - | strain-compatible soil properties |
| (SITE, menu only) | `.sit` | - | FILE1, FILE2 |
| POINT | `.poi` | FILE2 (SITE) | FILE3 |
| HOUSE | `.hou` | (`.sit` data for embedded models) | FILE4 (`modelname.N4`), FILE77 (incoherent) |
| FORCE | `.frc` | - | FILE9 |
| ANALYS | `.anl` | FILE1, FILE3, FILE4; + FILE9 (external loads); + FILE77 (incoherent) | FILE8 |
| COMBIN | (not stated) | FILE81, FILE82 (two renamed FILE8 files) | FILE8 |
| MOTION | `.mot` | FILE8 | TFU/TFI/ACC/RS text files, FILE13, frames |
| STRESS | `.str` | FILE4, FILE8 | FILE14, FILE15, `.THS/.TFU/.TFI`, `ESTRESS_<frame>.ess` (if SECDATAOPT,1) |
| RELDISP | `.rdi` | `.TFI` from MOTION | `.TFD`, `.THD`, frames |

Typical linear SSI order: (EQUAKE) -> (SOIL) -> SITE -> POINT -> HOUSE -> (FORCE) -> ANALYS ->
(COMBIN) -> MOTION -> RELDISP / STRESS.

NOTE: 9.12 has **no RUNSITE** command. Section 6.4.5 says the SITE menu item "performs the same
functionality as the RUNSITE command". Implement `RUNSITE,[model]` with the same semantics so that
a complete analysis can be scripted (Open Questions Q-1). There is also no RUN command for the
NONLINEAR (PANEL) module; it runs in batch mode (spec 05d).

---

## 3. Calculation commands (9.13)

### 3.1 Data needed by the section-cut commands [CORE][IO]

Workflow (manual 5.8.2):

1. `READSTR` loads element stresses onto the **original** model. It must come **before** CSECT.
2. Build a cut: `CUTADD`, `CUTVOL`, `SLICE` (Section 9.10).
3. `CSECT,<dest>,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>` builds the cross-section model in
   `<dest>` (Section 9.10.1).
4. `ACTM,<dest>`, then `CALCPAR` (one time) / `CALCMOI` / `CALCC` / `CALCM`.
5. For time histories: `CALCSECTHIST` (ASCII `.ess` frames) or `CALCSECTHISTDB` (binary stress
   database). These work from the original model, the cut and the plane.

Element stress components available for section cuts (from manual 6.5 STRESS options and Chapter
3):

| Element type | 6 components, in order | Reference axes | Units |
|---|---|---|---|
| SOLID | Sxx, Syy, Szz, Sxy, Sxz, Syz | **global** | F/L^2 |
| SHELL | Sx'x', Sy'y', Sx'y' (membrane), Mx'x', My'y', Mx'y' (plate bending) | **local element** x', y' | membrane F/L^2; moments F*L/L |
| TSHELL | not available | - | the manual says no text ESTRESS frames are generated for TSHELL |
| BEAM, SPRING | not in `.ess` frames (nodal forces) | - | excluded from section integration (derived) |

`.ess` file format (derived): the same layout as `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT` (manual
Chapter 3), with signed values at one time step:

```
<number of groups>
<etype> <group#> <ordered group#> <#elements>       (one line per group)
...
<etype> <group#> <ordered group#>                   (block header, per group)
<element#> <c1> <c2> <c3> <c4> <c5> <c6>            (one line per element)
...
```

Frames are named `ESTRESS_<framenumber>.ess` and written to the `.\NSTRESS\` subfolder when
SECDATAOPT = 1 (STRESS "Section Cut Options - Save Time History").

### 3.2 Section-cut integration algorithm (used by CALCPAR, CALCSECTHIST, CALCSECTHISTDB) [CORE] (derived)

The manual gives no formulas. This algorithm follows the manual's description of CSECT: pieces
clipped to the plane, shell thickness adjusted for obliquity, stresses taken at element centres.

**Section pieces** (one per element of the cut that the plane actually crosses):

- **SOLID**: compute signed node distances `d_k = (X_k - P) . ez`, with tolerance
  `tol = 1e-6 * model size`.
  - The element is crossed if `min d < -tol` and `max d > +tol`.
  - Intersection points: element nodes with `|d| <= tol`, plus edge points where `d` changes sign.
  - Order the points by angle about their mean, in the (ex, ey) plane, to get a convex polygon.
  - Piece data: polygon area `A`, centroid `c`, second moments (exact polygon formulas below).
- **SHELL** (mid-surface quad or triangle):
  - The intersection of the shell polygon with the plane is a segment of length `L`, unit
    direction `d` and midpoint `c`.
  - Effective section width `w = t / |m x ez|`, where `m` is the shell unit normal and `t` the
    shell thickness. This is the manual's "thickness adjusted when the shell face is oblique".
  - Piece = rectangle `L x w` with area `A = L w`.
  - If the shell lies in the plane (`|m x ez| < 1e-6`): skip it and warn.
- **Coincident faces/edges**: the plane may lie on a face shared by two solids, or an edge shared
  by two shells (for example a cut exactly at a floor level). Count each face/edge **once**: key it
  by its node set, and keep the element on the `-ez` side if it is in the cut, otherwise the `+ez`
  side element. Print an information message advising the user to move the plane to mid-height of
  an element row.

**Area properties** in local plane coordinates `(x', y')` relative to the section centroid `C`:

```
A = sum A_i ;   C = sum A_i c_i / A                       (c_i are global points on the plane)
Ixx = integral y'^2 dA ;  Iyy = integral x'^2 dA ;  Ixy = integral x' y' dA ;  Izz = Ixx + Iyy (polar)
polygon (vertices (x_k, y_k), k=0..n-1, closed, CCW):
    a_k = x_k y_{k+1} - x_{k+1} y_k
    A   = 1/2 sum a_k
    Cx  = 1/(6A) sum (x_k + x_{k+1}) a_k ;  Cy = 1/(6A) sum (y_k + y_{k+1}) a_k
    Iyy0= 1/12 sum (x_k^2 + x_k x_{k+1} + x_{k+1}^2) a_k     (about the local origin)
    Ixx0= 1/12 sum (y_k^2 + y_k y_{k+1} + y_{k+1}^2) a_k
    Ixy0= 1/24 sum (x_k y_{k+1} + 2 x_k y_k + 2 x_{k+1} y_{k+1} + x_{k+1} y_k) a_k
    then shift to C with the parallel-axis theorem
shell strip (2-D unit vectors in the plane: d along the segment, nu perpendicular to d):
    J_c = (A/12) (L^2 d d^T + w^2 nu nu^T)                    (second-moment tensor about the strip centroid)
    add A (c_i - C)(c_i - C)^T  (parallel axis), then read Ixx, Iyy, Ixy
```

**Resultants** (global first, then rotated to local):

```
SOLID piece:  S = global 3x3 stress tensor from (Sxx,Syy,Szz,Sxy,Sxz,Syz)
              f_i = A_i (S . ez)                          # traction times area
SHELL piece:  e1, e2 = shell local axes x', y' (same convention as the STRESS module output)
              S = Sx'x' e1e1 + Sy'y' e2e2 + Sx'y' (e1e2 + e2e1)   # plane-stress tensor, global
              f_i = A_i (S . ez)                          # = L t (S . nu_s): exact for oblique cuts
              plate bending (option, default ON; Q-12):
                  nu_s = in-shell unit normal to the cut line, oriented with nu_s . ez > 0
                  Mb   = Mx'x' e1e1 + My'y' e2e2 + Mx'y' (e1e2 + e2e1)
                  mb_i = L * ( m x (Mb . nu_s) )          # couple vector; convention M = integral z*sigma dz
F = sum f_i
M = sum (c_i - C) x f_i  +  sum mb_i                      # about the section centroid
Fx = F.ex ; Fy = F.ey ; Fz = F.ez ; Mx = M.ex ; My = M.ey ; Mz = M.ez
```

Sign convention: `F` and `M` are the actions of the **+ez side** of the structure on the **-ez
side** (traction with outward normal `+ez` of the -ez part). Positive `Fz` = tension.

Elements without an `.ess` record (beams, springs, missing data) are skipped with a warning listing
group/element numbers. If no piece remains, the command errors out.

### 3.3 CALCC (9.13.1) [CORE]

- **Syntax:** `CALCC` (no arguments).
- **Action:** computes the **volume centroid** of the active model in global coordinates. The
  same calculation is part of CALCPAR's output.
- **Algorithm** (derived):
  - `Xc = sum(V_e c_e) / sum(V_e)`.
  - SOLID volume: exact, by splitting into tetrahedra (5 or 6 for a hexahedron, 3 for a wedge, 1
    for a tetrahedron). `c_e` is the volume centroid.
  - SHELL: `V_e = area x thickness`, centroid of the mid-surface polygon.
  - BEAM: `V_e = length x section area` (from the R property), if the property is defined.
  - Springs, lumped masses and general elements are excluded.
  - On a CSECT cross-section model, which is extruded to unit thickness, the result must equal the
    area centroid. So extrude the plotting model symmetrically, **+/-0.5 along n** (derived).
- **Output:** `Xc = ..., Yc = ..., Zc = ...` and total volume.

### 3.4 CALCSECTHIST (9.13.2) [CORE][IO]

- **Syntax:**
  `CALCSECTHIST,<infile>,<cutnum>,<pointx>,<pointy>,<pointz>,<normalx>,<normaly>,<normalz>,<rightx>,<righty>,<rightz>,<sysno>,[ts],<outfile>`
- **Purpose:** batch the section-force calculation of one cut over every time step. Input is a
  history of element-centre stresses stored as one `.ess` frame per step. The result goes to
  `<outfile>`.

| Arg | Meaning |
|---|---|
| infile | Full path of the stress-history **list file**. "Similar format to the animation list file for PREP": one header line (ignored), then one `.ess` frame path per line (derived) |
| cutnum | cut to use (must exist; see Section 9.10) |
| pointx, pointy, pointz | a point on the section plane (global) |
| normalx, normaly, normalz | normal of the original cut plane, global |
| rightx, righty, rightz | local "right" (+X) direction, global |
| sysno | number under which the local system is stored (section 1.4) |
| ts | time step of each calculation. If `ts <= 0` (or blank), column 1 is the **integer position** in the infile |
| outfile | full path of the output force-history file |

- **Preconditions:** the original (full) model is active and the cut is defined. The stress frames
  belong to that model.
- **Algorithm:**
  ```
  frames = read_list(infile)                 # skip header line and blank lines; relative paths
                                             # resolve against infile's folder, then the CWD (derived)
  pieces = build_section(active_model, cut, P, n)   # geometry is the same for all steps
  for k, f in enumerate(frames, start=1):
      stress = read_ess(f)                   # 3.1 layout
      Fx..Mz = integrate(pieces, stress, ex, ey, ez, C)   # 3.2
      t = (k-1)*ts if ts > 0 else k          # frame 1 = time 0.000 (frame naming stress_00.000_00001)
      write row
  write MAX row; store local system sysno
  ```
- **Output format:** same as CALCSECTHISTDB (3.5). The manual does not say this explicitly
  (Q-8).

### 3.5 CALCSECTHISTDB (9.13.3) [CORE][IO]

- **Syntax:**
  `CALCSECTHISTDB,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sys>,<ts>,[start],[end],[Stride],<outfile>`
- **Purpose:** same as CALCSECTHIST, but the stresses come from the **binary stress database**
  currently loaded in UI memory (`Modelname_STRESS.bin`, written by STRESS when the "Save Binary
  Database" option / BINOUT str flag is set, and loaded with the Section 9.18 database commands).

| Arg | Meaning | Default |
|---|---|---|
| cutnum | cut number to use for the section history | required |
| px, py, pz | point on the plane of the cross section | required |
| nx, ny, nz | normal of the cross-section plane | required |
| rx, ry, rz | vector defining the local system's positive X | required |
| sys | system number for storing the local cut system | required |
| ts | time step of the stress database. If `0`, column 1 = step number | required |
| start | first step included in the table | 1 |
| end | last step included in the table | last time step in the database |
| Stride | number of time steps advanced between reported calculations | 1 |
| outfile | full path of the results file | required |

- **Preconditions (manual):**
  - The stress database is loaded into the UI.
  - The associated model is in memory **and is the active model**.
  - `cutnum` is defined.
- **Steps:** `k = start, start+Stride, ..., <= end`. Time `t = (k-1)*ts` (derived; same rule as
  3.4), or `k` when `ts = 0`. Also treat negative `ts` as 0 (derived).
- **Output format (manual):** a **7-column CSV ASCII table**.
  - Column 1 = simulation time, or step number if `ts = 0`.
  - Columns 2-7 = the six section resultants. Order (derived): `Fx, Fy, Fz, Mx, My, Mz` in the
    local cut system.
  - **Last line:** the word `MAX` in column 1, then for each component the value with the
    **largest absolute value, sign kept**. Ties go to the first occurrence (derived).
  - Header row (derived, configurable): `Time,Fx,Fy,Fz,Mx,My,Mz` (or `Step,...`). Number format
    `%.6E` (derived).

  ```
  Time,Fx,Fy,Fz,Mx,My,Mz
  0.000000E+00,1.2E+01,...
  ...
  MAX,-3.4E+02,5.1E+01,...
  ```

### 3.6 CALCM (9.13.4) [CORE]

- **Syntax:** `CALCM` (no arguments).
- **Action:**
  - Computes the total mass of the active model's elements from material properties and element
    volumes.
  - Computes the **total lumped masses** of the active model.
- **Algorithm** (derived):
  - Element mass `m_e = (w / g) V_e`, where `w` is the material **specific weight** (`M` command,
    `<weight>` field) and `g` the model gravity. Volumes are computed as in CALCC.
  - Excavated-soil elements (`ETYPE` = excavation) are reported separately and are not added to
    the structural mass (derived).
  - Lumped masses: sum of the nodal masses defined by `MT` (translational, per X/Y/Z) and `MR`
    (rotational), converted with the `MUNITS` option. If MUNITS says the values are weights,
    divide by `g`.
- **Output:**
  - element mass and weight,
  - lumped mass per direction (Mx, My, Mz) and rotational inertia sums,
  - total mass per translational direction.

### 3.7 CALCMOI (9.13.5) [CORE]

- **Syntax:** `CALCMOI,<normalx>,<normaly>,<normalz>,<rightx>,<righty>,<rightz>,<sysno>`
- **Action:** moments of inertia of the active model. Meant for **cut-plane (CSECT)
  models**.

| Arg | Meaning |
|---|---|
| normalx, normaly, normalz | normal of the original cut plane, global coordinates |
| rightx, righty, rightz | local right (+X), global coordinates |
| sysno | system number where the local system is stored |

- **Algorithm:**
  - Build the local system (1.4) with its origin at the area centroid.
  - Compute `A, Ixx, Iyy, Ixy, Izz` from the section pieces (3.2).
  - The pieces are recovered from the cross-section model: polygons and strips stored by CSECT.
    Alternatively, intersect the extruded unit-thickness elements with their mid-plane (derived).

### 3.8 CALCPAR (9.13.6) [CORE]

- **Syntax:** `CALCPAR,<normalx>,<normaly>,<normalz>,<rightx>,<righty>,<rightz>,<sysno>,[verbose]`
  (written `CalcPar` in the manual).
- **Action:** prints the **area centroid**, **moments of inertia** and the **current stress
  calculation** (six section resultants) of the active model. Meant for cross-sections built by
  CSECT.
  - The centroid is the origin of the local system used for the MOI and the stress resultants.
  - Stress data must be loaded onto the **original** model with READSTR **before** the
    cross-section is built. CSECT copies the element stresses into the cross-section model.
- **Args:** as CALCMOI, plus `verbose`, the output format flag. Default `-1`.
  - `verbose = -1` (default): a small table, one `<parname> = <value>` line per parameter.
  - Any other value: all parameters on **one line**, in the same order, without names (good for
    collecting many sections into one file).
- **Parameter list and order** (derived; the manual gives no names):

| # | Name | Meaning |
|---|---|---|
| 1 | `Area` | section area |
| 2-4 | `Xc`, `Yc`, `Zc` | area centroid, global |
| 5 | `Ixx` | integral y'^2 dA about local X through C |
| 6 | `Iyy` | integral x'^2 dA about local Y through C |
| 7 | `Ixy` | integral x'y' dA |
| 8 | `Izz` | polar moment Ixx + Iyy |
| 9 | `Fx` | shear force along local X |
| 10 | `Fy` | shear force along local Y |
| 11 | `Fz` | normal force (+ tension) |
| 12 | `Mx` | moment about local X |
| 13 | `My` | moment about local Y |
| 14 | `Mz` | torsion about local Z (normal) |

- The manual also says the calculation commands can be used on a model that is not a cross-section,
  for example the global base forces and moments of a whole building. Implement this through a
  cut that holds the whole structure and a CSECT plane at the base (derived; Q-11).
- **Manual example (5.8.2)**, a workflow pattern only (Q-10: it looks geometrically
  inconsistent):
  ```
  inp,model.pre
  readstr,stressfile.txt
  cutvol,3,52.5,52.8,-320,320,2.53,45.22
  csect,1,3,55.6,0,4,1,0,0
  actm,1
  calcpar,1,0,0,0,1,0
  calcm
  ```

### 3.9 READSTR (9.13.7) [IO]

- **Syntax:** `READSTR,<Filename>,[Dir]`
- **Action:** reads the stress results of the **current (active) model** from an `.ess` file and
  attaches them to the elements. Afterwards cross-sections can be built and section stresses
  computed (CALCPAR on the CSECT model; CALCSECTHIST).
- **Args:**
  - `Filename`: name of the `.ess` file to load. May be a full path.
  - `Dir`: directory that holds the file. If given, the path is `Dir/Filename`, otherwise
    `Filename` relative to the CWD (derived).
- **Rules (derived):**
  - Match records by **group number + element number**.
  - Warn for elements without data. Error if the group or element type in the file does not match
    the model.
  - The manual example reads `stressfile.txt`, so any extension is accepted.
- WARNING (derived): `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT` has the same layout, but its values are
  per-component **absolute maxima** that do not occur at the same time. Section resultants computed
  from it are not physically consistent. Warn when the file name matches.

### 3.10 SECDATAOPT (9.13.8) [IO]

- **Syntax:** `SECDATAOPT,<flag>`
- **Action:** sets the `*.ess` output request in the **STRESS module input file** (`.str`). This
  is the same as the STRESS options "Section Cut Options - Save Time History" check box.

| flag | Meaning |
|---|---|
| 0 | no `.ess` files written |
| 1 | save the `.ess` files for the **entire time history**: `ESTRESS_<framenumber>.ess` in `.\NSTRESS\` |

- Takes effect at the next AFWRITE / RUNSTRESS.

### 3.11 SHEAR (9.13.9) [CORE][ADV]

- **Syntax (as printed):** `SHEAR,<panel>,[fc],[fy],[P],[Nu],[Fvw],[Fbe]`
- **Purpose:** peak (ultimate) shear strength of one wall panel, or of all wall panels (Option NON
  panels defined with `P`/`PNLGEN`). Four empirical models are used (Gulec and Whittaker, 2009,
  give the details): ACI 318-08, Wood (1990), Barda et al. (1977) and Gulec-Whittaker (2009).
  The lower bound of Wood and the upper bound shared by Wood and ACI are reported too.
- **Arguments.** The argument list in the manual text differs from the syntax line (Q-14). This
  spec follows the descriptive list:

| # | Syntax name | Description name | Meaning | Units (British / SI) |
|---|---|---|---|---|
| 1 | panel | panel | panel number. `0` = all panels, one row per panel in a space-separated table | - |
| 2 | fc | fc | f'c, concrete compressive strength | ksi / kN/m^2 |
| 3 | fy | fy | reinforcement yield strength | ksi / kN/m^2 |
| 4 | P | P | rho_V or rho_H, web reinforcement ratio (one value is used for both) | - |
| 5 | Nu | Nu | axial force (compression positive, derived) | kips / kN |
| 6 | Fvw | **Abe** | A_BE, boundary-element vertical reinforcement area (Gulec-Whittaker only) | ft^2 / m^2 |
| 7 | Fbe | **Fybe** | f_y,BE, boundary-element vertical reinforcement yield stress (Gulec-Whittaker only) | ksi / kN/m^2 |

- Any optional argument left blank is **0**.
- Unit system detection follows 1.5. Results are reported in **kips** (British) or **kN** (SI).
- The panel geometry is taken from the model automatically for each panel: height `h_W`, length
  `l_W`, thickness `t_W` and web area `A_W = l_W t_W` (derived; Q-15).
  - `h_W` = vertical (global Z) extent of the panel group nodes.
  - `l_W` = horizontal extent in the panel plane.
  - `t_W` = shell thickness. BBCGEN requires uniform thickness and material in a panel.
- Using `panel = 0` is valid only if the concrete and reinforcement-ratio inputs are the same for
  all panels. This is not true for Gulec-Whittaker, which has more inputs.

**Equations.** Transcribed from the images on PDF pages 291-292. They are in the US-customary
form of Gulec and Whittaker (2009): `sqrt(f'c)` is in **psi**, areas in in^2, forces in lb.

```
Barda et al. 1977 (G&W eq. 2-7 / 4-7) - axial force included:
    V_Barda = ( 8 sqrt(f'c) - 2.5 sqrt(f'c) (h_W/l_W) + N_U/(4 l_W t_W) + rho_V f_y ) * t_W * (0.6 l_W)

Wood 1990 (G&W eq. 2-8) - axial force not included:
    6 sqrt(f'c) A_W  <=  V_Wood = rho_V A_W f_y / 4  <=  10 sqrt(f'c) A_W

ACI 318-08 (G&W eq. 4-1, "Chapter 11") - axial force not included:
    V_ACI = ( alpha_c sqrt(f'c) + rho_H f_y ) A_W  <=  10 sqrt(f'c) A_W

Gulec-Whittaker 2009 (G&W eq. 6-9) - axial force included:
    V_GW = ( 1.5 sqrt(f'c) A_W + 0.25 F_VW + 0.20 F_BE + 0.40 N_U ) / sqrt(h_W / l_W)
    F_VW = rho_V A_W f_y           (vertical web reinforcement strength)
    F_BE = A_BE f_y,BE             (boundary element vertical reinforcement strength)
```

- `alpha_c` is not defined in the manual. Use ACI 318-08 21.9.4.1 (ext. ref.): `alpha_c = 3.0` for
  `h_W/l_W <= 1.5`, `2.0` for `h_W/l_W >= 2.0`, linear in between (Q-16).

**Unit handling** [CORE] (derived; required because the coefficients only work with psi):

```
British (geometry ft, inputs ksi/kips): L_in = 12 L_ft ; A_in2 = 144 A_ft2 ; f_psi = 1000 f_ksi ; N_lb = 1000 N_kip
SI (geometry m, inputs kN/m^2, kN):     L_in = 39.3700787 L_m ; A_in2 = 1550.0031 A_m2
                                        f_psi = 0.1450377 f_kPa ; N_lb = 224.808943 N_kN
compute every V in lb, then convert:    kips = lb/1000 ;  kN = lb * 0.00444822162
```

- Equivalent SI-MPa check: `k sqrt(f'c[psi])` psi = `0.08304 k sqrt(f'c[MPa])` MPa.

**Output columns.** The manual text is ambiguous (Q-17). This spec uses 6 columns per panel:

| Col | Content (default interpretation) |
|---|---|
| 1 | panel number |
| 2 | ACI 318-08, capped at the upper bound: `min(V_ACI, 10 sqrt(f'c) A_W)` |
| 3 | Wood 1990, capped at the upper bound: `min(V_Wood, 10 sqrt(f'c) A_W)` |
| 4 | Wood 1990 lower bound `6 sqrt(f'c) A_W` |
| 5 | Barda 1977 |
| 6 | Gulec-Whittaker 2009 |

- A single-panel command prints a title row with the equation names first.
- `panel = 0` prints one space-separated row per panel.
- Also print, in an info line, the uncapped values and the unit system used, to support
  verification (derived).
- Recommended practice (manual):
  - Barda can **overestimate** the strength of typical nuclear shear walls. It applies to squat
    walls with heavily reinforced flanges (barbells).
  - ACI 318-08 may overestimate.
  - Wood and Gulec-Whittaker are close to the median of squat-wall tests.
  - Gulec-Whittaker is sensitive to the aspect ratio: for long panels it rises to Barda/ACI values
    or higher.
- Related: `BBCGEN` (9.17.3) uses the same strengths. Its ShearModel numbers are 1 ACI 318-08,
  2 Wood 1990, 3 Barda 1977, 4 Gulec-Whittaker 2009. The upper/lower bound applies when ACI or Wood
  falls outside its bounds.

**Numeric verification example** (computed for this spec; British):
- Inputs: `h_W` = 20 ft, `l_W` = 30 ft, `t_W` = 2 ft, f'c = 5 ksi, f_y = 60 ksi, rho = 0.005,
  N_U = 1000 kips, A_BE = 0.
- Intermediate values: `A_W` = 8640 in^2, sqrt(f'c) = 70.7107 psi, alpha_c = 3.0.

| Quantity | kips |
|---|---|
| V_ACI (raw = capped) | 4424.82 |
| upper bound 10 sqrt(f'c) A_W | 6109.40 |
| V_Wood raw | 648.00 |
| Wood lower bound 6 sqrt(f'c) A_W | 3665.64 |
| V_Barda | 4026.77 |
| V_GW (F_VW = 2592 kips, F_BE = 0) | 2405.90 |
| V_GW with A_BE = 0.1 ft^2, f_y,BE = 60 ksi (F_BE = 864 kips) | 2617.54 |

- The same panel in SI (6.096 m x 9.144 m x 0.6096 m, f'c = 34473.8 kN/m^2, f_y = 413685 kN/m^2,
  N_U = 4448.22 kN) must give V_ACI = 19682.6 kN, V_Barda = 17912.0 kN, V_GW = 10702.0 kN.

---

## 4. Plot commands (9.14)

### 4.0 Plot-state model [UI] (derived)

- The UI has a set of **plot tabs**; one of them is the **active plot**.
- Plot kinds:
  - 2D line plots: Spectrum Plot, Time History Plot, Soil Layer Plot, Soil Property Plot.
  - 3D model plots: Model (element) Plot, Node Plot, Cut Plot.
  - 3D animations: Bubble, Vector, Contour, Deformed Shape.
- A command meant for another plot kind is **ignored with a warning**, for example AXES on a 3D
  plot or SHRINK on a node plot. The manual says commands act only "if the plot is capable".
- Persistent, model-level display state survives across plots: NODESEL selections, hide requests.
- Toggle-argument conventions differ by command. Reproduce them exactly:

| Command | toggle value | on | off | default |
|---|---|---|---|---|
| ELENUM, GROUPNUM, NODENUM, SHOWMASS | -1 | 1 | 0 | -1 (toggle) |
| SHRINK, STIPPLE, WIREFRAME | -1 | 1 | 0 | -1 (toggle) |
| DEBUG | **2** | 1 | 0 | **2** (toggle) |
| PAUSE | -1 | **0 = start** | **1 = stop** | -1 (toggle) |

- Python implementation suggestion: matplotlib for 2D; VTK/pyvista (or OpenGL) for 3D.
  `CAPTUREPLOT` must render off-screen in batch mode.

### 4.1 Line-math commands (spectrum / time-history post-processing) [CORE]

All of them obey 1.3 (union grid, linear interpolation, constant extrapolation). Source and
destination numbers are line reference numbers. The destination may be one of the sources: compute
first, then overwrite (derived). Undefined source numbers are an **error** (derived; the plotting
commands ignore them, but silently dropping a term from a combination is unsafe). The same
operations are available as buttons in the 2D "Graph Plot Options" window (4.4.4).

#### 4.1.1 ADDITION (9.14.1)

- **Syntax:** `ADDITION,<dest>,<source1>,...,<source100>`
- `dest` = integer destination line. `source1..100` = integer input lines (up to 100).
- Result: `y = sum y_i`. This is LINECOMBIN with all coefficients 1. Default name
  `Linear Combin.`
- Use: algebraic summation of co-directional time histories (manual 1.x item xiii).

#### 4.1.2 AVERAGE (9.14.2)

- **Syntax:** `AVERAGE,<dest>,<source1>,...,<source100>`
- Result: arithmetic mean `y = (1/N) sum y_i` on the union grid. Default name `Average Line`.
- Use: averaging ISRS from several input-motion sets or analyses.

#### 4.1.3 BROADEN (9.14.4) [CORE]

- **Syntax:** `BROADEN,<Dest>,<Smooth1>,<Smooth2>,<source1>,...,<source100>`

| Arg | Meaning | Dialog label (Graph Plot Options, "Broadening and Enveloping Spectra") | Dialog default |
|---|---|---|---|
| Dest | destination line | (Post Processing Results - Line number) | - |
| Smooth1 | **peak bridging percentage** | `Peak Difference(%)` | 0 |
| Smooth2 | **peak broadening percentage** | `Broaden(%)` | 0 |
| source1..100 | input lines (up to 100) | checked lines | - |

- Result name `Envelope`. The command **envelopes** the sources, then broadens. It uses the union
  grid and constant extrapolation (1.3).
- Percentages are entered as percent: `15` means 15%.

**Algorithm** (derived; the manual names the parameters but gives no formula; ext. ref. ASCE 4
ISRS smoothing/peak-broadening practice):

```
BROADEN(dest, p_bridge%, p_broad%, sources):
    b = p_broad/100 ; p = p_bridge/100              # require 0 <= b < 1 ; 0 <= p <= 1
    G = union_grid(sources)
    E(f) = max_i value_on(source_i, f)               # 1) envelope (piecewise linear on G)
    # 2) peak broadening: every point of E is spread over [f(1-b), f(1+b)]
    G2 = G  U  {g(1-b), g(1+b) for g in G}, clipped to [min G, max G]     # extra breakpoints
    B(f) = max{ E(f') : f/(1+b) <= f' <= f/(1-b) }  for f in G2
           # exact for piecewise-linear E: max over E at both window ends and at all vertices of
           # G strictly inside the window
    # 3) peak bridging (only if p > 0): fill shallow valleys between adjacent peaks
    repeat until no change:
        peaks = local maxima of B (a flat plateau counts as one peak)
        for each adjacent pair (P_i, P_j), from low to high frequency:
            A_low = min(A_i, A_j) ; V = min of B between them
            if (A_low - V) / A_low <= p:                    # valley-depth criterion (default, Q-19)
                B(f) = max(B(f), A_low) for f in (f_i, f_j)  # insert the crossing points of
                                                            # B = A_low so the fill is exact
    store dest = B on G2 (or on G, if the "ACS-compatible grid" option is set; see below)
```

- **Step 2 is the standard plus/minus b% peak broadening.** A sharp peak at `f0` becomes a flat
  plateau over `[f0(1-b), f0(1+b)]`. The left flank shifts by `(1-b)` and the right flank by
  `(1+b)`.
- Grid choice:
  - The default output grid is the augmented `G2`, which gives exact plateau edges.
  - Option "ACS-compatible grid" (decision): evaluate only on `G`, using the discrete form
    `B(f_k) = max{E(f_j): f_j in [f_k/(1+b), f_k/(1-b)]}`.
- The manual warns that **fewer frequency steps reduce the accuracy of the UI spectrum broadening
  algorithm**. For nuclear safety-related ISRS, use 0.1-100 Hz with **at least 301** frequency
  steps (US NRC SRP 3.7.1). This warning suggests the original works on the discrete points. Print
  it when a source has fewer than 301 points.
- Typical values (ext. ref.):
  - ASCE 4 ISRS practice: broaden peaks by +/-15%.
  - US NRC RG 1.122: frequency-dependent broadening, at least +/-10%.
  - `Smooth1 = 0` means no bridging (the dialog default).
- Algorithm-order decision: envelope, then broaden, then bridge.

**Verification tests for BROADEN** (computed with the algorithm above):

- **T-B1 (broadening).**
  - Source: `x = [1,5,9,10,11,15,20]`, `y = [1,1,1,5,1,1,1]`.
  - Command: `BROADEN,2,0,15,1`.
  - Expected `B = 5` on `[8.5, 11.5]`.
  - Linear ramp from (7.65, 1) to (8.5, 5), and linear ramp from (11.5, 5) to (12.65, 1).
  - `B = 1` elsewhere.
  - Grid points include 7.65, 8.5, 9.35, 10.35, 11.5, 12.65.
- **T-B2 (bridging).**
  - Source: `x = [1,4,5,6,6.5,7,8,9,12,20]`, `y = [0.5,2,4,3,3.5,3.0,3.8,2,1,0.5]`.
  - Command: `BROADEN,2,15,0,1`.
  - Expected: values unchanged up to x = 5 (4.0), then 3.8 at x = 5.2 (crossing), **flat 3.8 from
    5.2 to 8**, then unchanged.
  - The same result is obtained with the alternative "peak amplitude difference" criterion (Q-19).
  - With `Smooth1 = 10` the two criteria differ: the valley-depth criterion gives no change, while
    the amplitude criterion bridges. Use this case to pin down the original behaviour.
- **T-B3 (envelope).** Two lines with peaks at different frequencies, `b = 0`, `p = 0`: the result
  must be the pointwise maximum on the union grid.

#### 4.1.4 LBINCORS (9.14.18) [ADV]

- **Syntax:** `LBINCORS,<out>,<in>`
- `out` = reference number where the output line is stored. `in` = reference number of the input
  line.
- Summary-table action: "Calculate the Lower Bound Incoherent Response Spectra".
- The manual gives **no algorithm**, default name or parameters. Implement the parser and the
  line-store plumbing. Leave the transformation as a pluggable rule that raises "not specified"
  until the rule is confirmed (Q-20). Do not invent a reduction rule silently.

#### 4.1.5 LINECOMBIN (9.14.20)

- **Syntax:** `LINECOMBIN,<Dest>,<Line1>,<Coeff1>,...,[Line100],[Coeff100]`
- Result: `y = sum c_i y_i`. Each line's y values are scaled by its coefficient and the scaled
  lines are added.
- A coefficient **must** be given for every line. An odd number of arguments after Dest is an
  error.
- Up to 100 pairs. Default name `Linear Combin.`
- The window button handles only up to 3 lines ("Spectra 1/2/3" coefficients, default 1). The
  command has no such limit.
- Use: weighted linear combination of co-directional ISRS, e.g. `1.0 X + 0.4 Y + 0.4 Z` (100-40-40
  rule).

#### 4.1.6 SRSS (9.14.41)

- **Syntax:** `SRSS,<dest>,<source1>,...,[source100]`
- Result: `y = sqrt(sum y_i^2)` on the union grid. Default name `SRSS Line`.
- Use: SRSS combination of co-directional ISRS (X, Y, Z input directions). The result is
  non-negative, so it is meant for spectra and amplitudes, not signed histories.

#### 4.1.7 SUBTRACTION (9.14.43)

- **Syntax:** `SUBTRACTION,<dest>,<source1>,...,<source100>`
- Result: `y = y_1 - sum_{i>=2} y_i`. Up to 99 lines are subtracted from the first. This is
  LINECOMBIN with coefficients (1, -1, -1, ...). Default name `Linear Combin.`

**Line-math verification test T-L1:**
- Lines: `A: x=[0,1,2], y=[0,1,2]` and `B: x=[0.5,1.5], y=[10,20]`.
- Union grid: `[0,0.5,1,1.5,2]`. On that grid B is `[10,10,15,20,20]`.

| Command | Expected y |
|---|---|
| ADDITION | [10, 10.5, 16, 21.5, 22] |
| AVERAGE | [5, 5.25, 8, 10.75, 11] |
| SUBTRACTION A-B | [-10, -9.5, -14, -18.5, -18] |
| LINECOMBIN 2A+0.5B | [5, 6, 9.5, 13, 14] |
| SRSS | [10, 10.01249, 15.03330, 20.05617, 20.09975] |

### 4.2 Line file I/O [IO]

#### 4.2.1 READSPEC (9.14.31)

- **Syntax:** `READSPEC,<SpecFile>,<numLines>,<Line1>,...,<LineN>`
- `SpecFile` = full path of the spectrum file. `numLines` = number of data columns to read.
  `Line_i` = reference number for each line read.
- Format: column 1 is the **frequency** (x). Data column `i` (i = 1..numLines, counted after the
  frequency column) goes to `Line_i`.
- **Do not count the frequency column** in numLines. For TFU, TFI, TFD and RS files the frequency
  column is excluded. Example: an RS file has frequency + acceleration, so `numLines = 1`.
- Columns are read **from the first data column**; you cannot skip columns. To get column 3 you
  must read columns 1-3 (manual 7.2.1).
- Parsing (derived): whitespace- or comma-separated numbers. Lines that do not parse fully as
  numbers (headers) are skipped. Fortran `D` exponents are accepted. Error if a row has fewer than
  `1 + numLines` numbers.
- Line name = file name (derived). Existing lines are overwritten.
- Manual example: `READSPEC,$1$,1,1` inside a macro.

#### 4.2.2 READTH (9.14.32)

- **Syntax:** `READTH,<THFile>,<Pair>,<Num>`

| Pair | Format |
|---|---|
| 0 | one-column time history (**ACS SASSI output format**) |
| 1 | two-column time history, time/acceleration pairs on each line |

- `Num` = reference number for the line.
- Pair 0 parsing (derived from the MOTION "File Contains Pairs" option text): the first value is
  the **time step** dt, followed by one acceleration value per record. Accept several values per
  line. `t_k = (k-1) dt`, starting at t = 0 (Q-22).
- Pair 1: each line is `t a`.

#### 4.2.3 WRITESPEC (9.14.48)

- **Syntax:** `WRITESPEC,<SpecFile>,<Num1>,...,[Num50]`
- Writes up to **50** lines into one spectrum-format file. Every line is resampled onto the union
  of all their x-values (1.3).
- Output (derived): one row per grid point: `x y_1 y_2 ... y_n`, whitespace-separated, `%14.6E`,
  no header. This lets READSPEC read it back exactly.

#### 4.2.4 WRITETH (9.14.49)

- **Syntax:** `WRITETH,<THFile>,<Num>`
- Writes one line in the **single-column** time-history format. The line must have a **constant
  time interval**, which is computed from the first two points: `dt = x2 - x1`.
- Output (derived): dt on the first line, then one y value per line (mirror of READTH Pair 0).
- Check (derived): warn if any interval differs from dt by more than 1e-6 dt.

### 4.3 2D plot creation and formatting [UI]

| Command | Syntax | Behaviour |
|---|---|---|
| SPECPLOT (9.14.40) | `SPECPLOT,<Line1>,...,<Line50>` | New spectrum plot of up to 50 lines. The default extent is the min/max of the plotted lines. Numbers with no line are **ignored**. A value of **-1 ends the list**: anything after it is ignored. Same as Plot > Spectrum TFU-TFI |
| THPLOT (9.14.44) | `THPLOT,<Line1>,...,<Line50>` | Same rules for a time-history plot. Same as Plot > Time History |
| LAYERPLOT (9.14.19) | `LAYERPLOT` | Soil layer plot of the active model: layer column and table (4.4.8) |
| SOILPROPPLOT (9.14.39) | `SOILPROPPLOT,<PropName>` | Soil property plot (G/Gmax and damping vs shear strain %) for the named dynamic soil property. `PropName` is **case-sensitive** and must match a property in memory, or no plot appears. Without it, the "Select Dynamic Soil Property" window opens (4.4.9) |
| AXES (9.14.3) | `AXES,<MaxTickX>,<MaxTickY>,<MinTickX>,<MinTickY>,<LogX>,<LogY>` | Axis control on the active 2D plot. Each argument is 0/1: MaxTick = thick (major) grid lines on X/Y; MinTick = thin (minor) grid lines on X/Y; LogX/LogY: 0 standard, 1 logarithmic |
| PLOTRANGE (9.14.28) | `PLOTRANGE,<Xmin>,<Xmax>,<Ymin>,<Ymax>` | Set the axis ranges of the active 2D plot. Error if min >= max; error on a log axis if min <= 0 (derived) |
| PLOTTITLE (9.14.29) | `PLOTTITLE,<Title>` | Title of the active plot. Works for 2D line plots and **all 3D plots**. **Ignored by the soil layer plot** |
| XTITLE (9.14.50) | `XTITLE,<label>` | X-axis title of the current graph |
| YTITLE (9.14.51) | `YTITLE,<label>` | **Left** Y-axis title |
| YTITLE2 (9.14.52) | `YTITLE2,<label>` | **Right** Y-axis title (e.g. "Shear Modulus" on the soil property plot) |
| LINENAME (9.14.21) | `LINENAME,<Num>,<Name>` | Rename a line object. Global, so all graphs show the new name |
| MARKERS (9.14.22) | `MARKERS,<Mark>,<Ln1>,...,<Ln50>` | Data-point markers 0 off / 1 on for up to 50 lines. Global line property |
| STIPPLE (9.14.42) | `STIPPLE,<switch>` | Dash pattern on the lines of the active plot, for black-and-white prints. 0 off, 1 on, -1 (default) toggle. Each line gets a distinct dash pattern (derived) |

- Text arguments (titles, names) take the rest of the argument verbatim, trimmed. A title cannot
  contain a comma unless quoted (derived; Q-24).

### 4.4 Dialogs opened by plot commands (transcribed) [UI]

#### 4.4.1 "Line Selection" (Plot > Spectrum / Time History, before plotting; PDF 200)

| Group | Field | Type | Notes |
|---|---|---|---|
| (top) | line list | check-list | entries `n: name`, e.g. `1: Test2.rs` |
| | Ok / Cancel | buttons | |
| Title_Axis Labels | Title, X-Label, Y-Label | text | |
| X Axis Options | Min, Max | text | |
| | Logarithmic, Show Ticks | check | Show Ticks = minor ticks |
| Y Axis Options | Min, Max, Logarithmic, Show Ticks | as X | |
| Input Line File | File Name | text + `<<` browse | |
| | Starting Number | text | first line reference number |
| | Lines in file | text | number of data columns to read (numLines) |
| | Add Line(s) | button | reads the file into consecutive numbers starting at Starting Number |

#### 4.4.2 "Graph Plot Options" (WINDOWSETTINGS on a 2D plot; PDF 202)

| Group | Field | Default (screenshot) | Command equivalent |
|---|---|---|---|
| line list | check-list (`1: Test.rs`) | - | lines shown |
| Title_Axis Labels | Title, X-Label, Y-Label | blank | PLOTTITLE, XTITLE, YTITLE |
| X Axis Options | Min, Max, Logarithmic, Show Ticks | 0.10, 100.00, off, off | PLOTRANGE, AXES |
| Y Axis Options | Min, Max, Logarithmic, Show Ticks | 0.02, 4.12, off, off | PLOTRANGE, AXES |
| Line Options | Data Points (highlighted line) | off (grey) | MARKERS |
| | Line Stippling (all lines) | off | STIPPLE |
| Spectra Analysis Postprocessing | Spectral Average: `Average` button | - | AVERAGE |
| | SRSS Combination: `SRSS` button | - | SRSS |
| Broadening and Enveloping Spectra | Peak Difference(%) | 0 | BROADEN Smooth1 |
| | Broaden(%) | 0 | BROADEN Smooth2 |
| | `Broaden` button | - | BROADEN |
| Linear Combination and Coefficients | Spectra 1, Spectra 2, Spectra 3 | 1, 1, 1 | LINECOMBIN (up to 3 lines) |
| | `Linear Combin` button | - | LINECOMBIN |
| Temporal Analysis Post Processing | `Addition`, `Subtraction` buttons | - | ADDITION, SUBTRACTION |
| Post Processing Results | Line number | blank | destination line (derived: blank = next free number) |
| | Ok, Cancel | | |

- The buttons act on the **checked** lines. The new line is appended to the list.
- The buttons do **not** write files. Use WRITESPEC / WRITETH for that.

#### 4.4.3 "Window Options" (WINDOWSETTINGS on a 3D plot; PDF 194)

| Group | Fields | Default / notes |
|---|---|---|
| Model Display Volume | X, Y, Z: Min, Max | model bounding box. Geometry outside is hidden |
| Colormap Value Range (Bubble and Contour only) | Min, Max | Values below Min show the min colour (dark blue); above Max, the max colour (dark red) |
| Output Direction (Vector only) | X, Y, Z, All (radio) | X |
| Show Element Group | check-list `Group <n> <TYPE> - <#elements>` | all checked. Space bar toggles the highlighted rows |
| Hide/Show Elements (Nodes on node plots) | Hide/Show radio; Group; Elem. Numbers; `Hide/Show Elem.` button | The list is space-separated, or one range `a-b`. One group at a time. **Applied immediately** when the button is pressed, even without OK. Hide requests are stored **in the model** |
| Animation Options | Scale Factor (1.00), Frame Pause (ms) (33) | Scale applies to vector and deformed plots. Frame pause applies to all graphs |
| | Show Undeformed Shape | Deformed Shape plot only |
| | Title | plot title |
| | OK, Cancel | |

- The manual notes that a future version should take arguments for WINDOWSETTINGS. Optional
  extension: `WINDOWSETTINGS,<field>,<value>` (derived, not required).

#### 4.4.4 "Shader Options" (SHADEROPTIONS with no arguments; Options > Shader Options; PDF 195)

| Field | Default | SHADEROPTIONS argument |
|---|---|---|
| Node/Bubble Node Size | 10.00 | `points` (max point size, node/bubble plots; also vector plots per 7.1.1) |
| Vector/Displacement Scale Factor | 1.00 | `scale` (animation scale factor) |
| Element Outline Thickness (% of element) | 0.02 | `linew` (outline width of solid elements: the fraction of the element edge drawn in the outline colour, for Element, Cut, Contour and Deformed plots) |
| Element Shrink (% of element) | 0.06 | `shrink` (face shrink in the Element plot) |
| Cancel / OK | | OK updates the active plot |

#### 4.4.5 "Load Frame Data" (animation plots with no arguments; PDF 210)

| Group | Field | Notes |
|---|---|---|
| Select From Database | list with columns Description, Animation Directory, Type, Frames | entries come from `SASSIani.xml` |
| | Remove Animation | Deletes the entry **and** the processed frame files |
| Animation Control - Frame Selection | Start (1), End (number of frames), Stride (1) | defaults come from SASSIani.xml |
| Data Scale | Scale Factor (1.000000) | deformed/vector plots. Bubble/contour plots show a colormap Min/Max range instead |
| | Ok, Cancel | |

#### 4.4.6 "Parse Frame Data" (PROCFRAME with blank AniFile or BufferDIR; PDF 209)

| Field | Notes |
|---|---|
| List File Name + `<<` | animation frame list file |
| Frame Storage Dir + `<<` | processed (binary) frame directory |
| Data Description | free text stored in the database |
| Plot Type (radio) | Bubble, Vector, Contour, Time History, Stress DB (Binary), ACC DB (Binary), RelDisp DB (Binary) |
| Ok, Cancel | |

- The last three types build animations from the SSI response **binary databases**.

#### 4.4.7 "Select Cut to Display" (CUTPLOT with no arguments; PDF 198)

- Fields: Cut Number (default 1), Model Number (default 1), Ok, Cancel.

#### 4.4.8 "Soil Layer Windows Setting" (WINDOWSETTINGS on the soil layer plot; PDF 206)

| Field | Default |
|---|---|
| Start Layer | 1 (surface level) |
| EndLayer | -1 (deepest layer in the model) |
| Show Thickness | off |
| Show Specific Weight | on |
| Show P-Wave Velocity | on |
| Show S-Wave Velocity | on |
| Show P-Wave Damping Ratio | on |
| Show S-Wave Damping Ratio | on |

- The plot shows the layers as a column with relative thickness, coloured from the SoilLayer
  palette.
- Table columns: Layer, Unit Weight, P-Wave Velocity, S-Wave Velocity, P-Wave Damping Ratio,
  S-Wave Damping Ratio. The last row is `Halfspace`.
- The data comes from the `L` and `TOPL` commands.

#### 4.4.9 "Select Dynamic Soil Property" (SOILPROPPLOT with no name; PDF 208)

- A property name list (e.g. Clay, Rock, Sand) with buttons New, Edit, Delete, Ok, Cancel.
- A table with columns `Strain | Mod. Red. | Strain | Damp`, and a Title field.
- Plot: damping ratio on the left Y axis, shear modulus (reduction) on the right Y axis, shear
  strain % on X.

#### 4.4.10 "Soil Properties Window Settings" (WINDOWSETTINGS on the soil property plot)

- G,D Axis group: Logarithmic, Show Ticks.
- Shear Strain Axis group: Logarithmic, Show Ticks.
- Show Shear Modulus Line (on), Show Damping Line (on), Ok, Cancel.

### 4.5 3D model plots, animations and frame processing [UI][IO]

| Command | Syntax | Behaviour |
|---|---|---|
| MODELPLOT (9.14.23) | `MODELPLOT` | Element plot of the active model, coloured by group/material/property (ELECOLOR) |
| NODEPLOT (9.14.25) | `NODEPLOT` | Node plot of the active model. Only nodes connected to elements are shown. Non-interaction nodes are black, interaction nodes red |
| CUTPLOT (9.14.12) | `CUTPLOT,[Cut],[Model]` | Wireframe of the whole model with the cut's elements filled. With no arguments, the "Select Cut to Display" window opens. This is the **only plot that can show a non-active model** |
| BUBBLEPLOT (9.14.5) | `BUBBLEPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>` | Bubble (ZPA / nodal value) plot. Bubble size and colour follow the nodal data and the colour bar |
| CONTOURPLOT (9.14.11) | `CONTOURPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>` | Contour plot on element faces. Colours are interpolated linearly across faces from nodal data |
| DEFORMPLOT (9.14.14) | `DEFORMPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>` | Animated deformed shape from nodal displacement frames. The undeformed wireframe is optional (Window Options) |
| VECTORPLOT (9.14.45) | `VECTORPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>` | Three vectors per used node: red = x, green = y, blue = z. A complex component tilts the vector out of its axis. Manual example: x-data 5+0.3i is drawn as x = 5, y = z = 0.3 |
| PROCFRAME (9.14.30) | `PROCFRAME,[AniFile],[BufferDIR],[Data],[anitype]` | Converts an animation frame list into the binary frame buffer (see below) |

Animation argument meanings (BUBBLEPLOT, CONTOURPLOT, DEFORMPLOT, VECTORPLOT):

| Arg | Meaning |
|---|---|
| BufferDir | processed frame storage directory (output of PROCFRAME) |
| MnF | starting frame number |
| MxF | end frame number |
| ST | stride |
| MnR, MxR | minimum / maximum data value of the colormap range (bubble, contour) |
| Col | data column selection in the frame files (bubble, contour) |
| Scale | scalar that magnifies the deformation in the animation (deformed, vector) |

- With **no arguments**, these commands open the "Load Frame Data" window.
- Partial argument entry may give **inconsistent results**. Recommended (derived): require all
  arguments, or none. If some are missing, fill them from the SASSIani.xml defaults and warn.
- "Used nodes" means nodes connected to elements. K-nodes and orphan nodes are not drawn.

PROCFRAME details:

| Arg | Meaning |
|---|---|
| AniFile | full path of the animation frame list file |
| BufferDIR | processed frame storage directory |
| Data | label for the animation, stored in the database |
| anitype | integer tag for sorting the database only. All frame data is processed the same way: 0 Bubble (default), 1 Vector, 2 Contour, 3 Time History |

- If AniFile or BufferDIR is blank, the "Parse Frame Data" window opens.
- Processing can take long. A progress bar appears at the bottom of the screen, and **user input
  is blocked** until processing finishes.
- Frame list file: the same format as the PREP frame listing files (`*.dispani`, `*.tfiani`,
  `*.impani`; also `zpani`, `thiani`, `contani`, `thani`). The **first header line is ignored**.
  Every other line names one frame file.
- Frame files: produced by MOTION/RELDISP/STRESS restart options in `\TFU`, `\RS`, `\ACC`,
  `\ACCR`, `\THD`, `\NSTRESS`, `\SOILPRES`. Naming examples from Table 3.2: `RS01_000.10_00001`,
  `TFU_000.02_00001`, `ACC_00.000_00001`, `THD_00.000_00001`, `stress_00.000_00001_sig`,
  `pres_00.000_00001_nod`. Maximum-value frames: `stress_ABS_MAX_<comp>`, `pres_ABS_MAX_<type>`.
- The header's second number is the number of columns (see MODFRAMES, 9.7.24).
- Output: binary frame files in BufferDIR (users must not edit them, only delete them), plus an
  entry in `SASSIani.xml`.
  - SASSIani.xml lives in a system-dependent default location.
  - Entry fields: description, directory, type, frame count, default start/end/stride, scale,
    colormap range (derived).
- Python implementation (decision): one `.npy` (or one `.npz`/HDF5 file) per frame. Keep the
  `SASSIani.xml` name for nomenclature, holding the same fields.

HARMFRAME (SASSI-EDU extension, requirements §3.4.R) — steady-state harmonic frames at one frequency:

`HARMFRAME,<Src>,<Freq>,<OutDir>,[NFrames],[Ref]`

| Arg | Meaning |
|---|---|
| Src | transfer-function source: a FILE8-type file (FILE8, FILE8X, FILE81 …) or the `TFU` frame folder / list file of MOTION *Restart for TF* (Re/Im of X, Y, Z) |
| Freq | frequency in Hz; the computed (SSI) frequency closest to it is used and printed, with a warning when it lies outside the computed range or more than 5 % away |
| OutDir | frame folder (created; earlier `HARM_*` frames in it are deleted; other files are reported because PROCFRAME of the folder would mix them) |
| NFrames | frames per period, 4 … 360 (default 24) |
| Ref | blank: total motion per unit control motion; 0: relative to the free field (unit control motion, zero phase, in the control direction `cm`/`ang` of the FILE8, as RELDISP without RELFILE; FILE8 sources only; a vibration FILE8 has no free field); n: relative to node n |

- Frame k (1-based) holds, for every node with a translational DOF in the source, the displacement
  u_k = Re(H e^{iφ_k}) = Re(H) cos φ_k − Im(H) sin φ_k at φ_k = 2π(k−1)/N, with H the complex transfer
  function of X, Y and Z at the chosen frequency (fixed DOFs 0). The time factor e^{iωt} is that of the
  convolution a(t) = IFFT[H A]: a DOF with H = |H| e^{iθ} moves as |H| cos(ωt + θ) while the control
  motion moves as cos ωt.
- Frame names `HARM_<φ_k in degrees, 000.0>_<k>` (Table 3.2 pattern), header `nrows 4`, rows
  `node ux uy uz` (D-FIL-03): PROCFRAME stores them with layout `xyz`, DEFORMPLOT draws x + Scale u, and
  the frame label is `ωt = <φ>°`.
- HOUSE optimizer: when `<model>.map` lies beside the source, the FILE8 (new) node numbers are written
  in the model numbering, which DEFORMPLOT draws.
- Output: the frequency used (and its frequency number or TFU frame), the period, the reference, the
  number of frames and nodes, and the largest amplitude with its node, direction and phase.

### 4.6 3D view and display controls [UI]

| Command | Syntax | Behaviour |
|---|---|---|
| CNGCENTER (9.14.8) | `CNGCENTER,<X>,<Y>,<Z>` | Set the centre of rotation of the 3D plot (global coordinates) |
| RSTCENTER (9.14.33) | `RSTCENTER` | Reset the centre to the default: the centre of the bounding box of all nodes **used by elements** (global) |
| CNGVIEW (9.14.9) | `CNGVIEW,<rX>,<rY>,<rZ>,<px>,<py>,<zoom>` | Set the view: rotations about X/Y/Z (degrees, derived), horizontal/vertical screen pan, zoom constant. The current values are shown when DEBUG mode is on. Used to view two plots from exactly the same viewpoint. Same as the toolbar "Change View" button |
| RSTVIEW (9.14.34) | `RSTVIEW` | Reset to the default view: **top-down** view at first opening. Zoom and location fit the bounding box of element-connected nodes |
| DEBUG (9.14.13) | `DEBUG,[switch]` | Debug overlay with view angles and animation info. 0 off, 1 on, **2 (default) toggle** |
| ELECOLOR (9.14.15) | `ELECOLOR,<val>` | Element colouring: **1 = Group, 2 = Material, 3 = Property** (beam/spring property). Elements without a material/property use a default colour. Same as toolbar "Group Colors / Material Colors / Property Colors". Colours come from the ElemPalette (up to 128 entries; index = number mod 128) |
| ELENUM (9.14.16) | `ELENUM,[opt]` | Element number labels: -1 toggle (default), 0 off, 1 on. **Mutually exclusive with GROUPNUM**: turning on ELENUM hides group numbers |
| GROUPNUM (9.14.17) | `GROUPNUM,[opt]` | Group number labels: -1 toggle (default), 0 off, 1 on. Mutually exclusive with ELENUM |
| NODENUM (9.14.24) | `NODENUM,[opt]` | Node number labels: -1 toggle (default), 0 off, 1 on (the manual text wrongly says "element numbers") |
| NODESEL (9.14.26) | `NODESEL,<N1>,...,<N20>` | Highlight up to 20 nodes on **all current and future** 3D plots (shown as a blue square border on node plots). Nodes are not deselected between commands. Calling NODESEL with an already-selected node **deselects** it (toggle per node). Repeat the command for more than 20 nodes |
| SHOWDOF (9.14.36) | `SHOWDOF,[label1],...,[label6]` | Mark nodes with a **fixed** DOF in the requested directions (green border on node plots). Labels: `X, Y, Z, XX, YY, ZZ, DISP` (X,Y,Z), `ROT` (XX,YY,ZZ), `ALL`. With no labels, the "Boundary Conditions" window opens: check boxes UX, UY, UZ, ROTX, ROTY, ROTZ and Ok / Cancel (Figure 9.7) |
| SHOWMASS (9.14.37) | `SHOWMASS,[opt]` | Lumped-mass markers: -1 toggle (default), 0 off, 1 on. On the element plot, marker sections show the direction of the mass: **red = x, green = y, blue = z**. Other plots only show that a mass exists (purple border on node plots) |
| SHRINK (9.14.38) | `SHRINK,[switch]` | Shrink element faces (Element plot only). 1 on, 0 off, -1 (default) toggle. Amount set by SHADEROPTIONS shrink |
| WIREFRAME (9.14.47) | `WIREFRAME,<Switch>` | Wireframe mode of the element plot. 0 off, 1 on, -1 (default) toggle |
| SHADEROPTIONS (9.14.35) | `SHADEROPTIONS,[points],[linew],[shrink],[scale]` | Change only the arguments that are filled in. With **no arguments**, the Shader Options window opens (4.4.4). Note that the argument order differs from the window's field order |
| PAUSE (9.14.27) | `PAUSE,[pz]` | Animation: -1 toggle start/stop (default), **0 start, 1 stop**. Keyboard: Pause key; `+`/`-` step frames while paused |
| COLOR (9.14.10) | `COLOR,<Palette>,<Num>,<R>,<G>,<B>` | Set colour `Num` of palette `Palette` to RGB, each 0-255. Same as Options > Colors |
| CAPTUREPLOT (9.14.6) | `CAPTUREPLOT,<FileName>` | Save the active plot image. **PNG** if `.png` appears **anywhere** in FileName (case-insensitive, derived), otherwise **Bitmap (BMP)**. Error if no plot is active. Same as File > Export Image |
| CLOSEPLOT (9.14.7) | `CLOSEPLOT` | Close the active plot window. Warn if there is none |
| WINDOWSETTINGS (9.14.46) | `WINDOWSETTINGS` | Open the settings window of the active plot: 2D "Graph Plot Options", 3D "Window Options", or the soil layer / soil property settings |

- SHOWDOF with labels sets the displayed set to exactly those directions (derived).
  `SHOWDOF,NONE` turns it off (derived extension; Q-25).
- COLOR palette names (Options > Colors tabs, manual 6.5.6):

| Palette | Controls |
|---|---|
| UI | command history text colours by message type |
| Element | element plot colours other than element fills |
| ElemPalette | element fills by group/material/property (max 128 colours; index = num mod 128) |
| Node | node plot |
| Cut | cut plot |
| Bubble | bubble plot |
| Vector | vector plot |
| Contour | contour plot |
| Deformed | deformed shape plot |
| Spec | spectrum plot axes |
| SpecLines | spectrum line colours, by line number |
| TimeHist | time-history plot axes |
| THLines | time-history line colours, by line number |
| SoilLayer | soil layers plot |
| SoilProp | soil properties plot |

- Palette names are case-insensitive. `Num` is 1-based within the palette (derived; Q-26).
- Mouse controls on 3D plots (manual 7.1.1):
  - right-drag: rotate;
  - left-drag vertically: zoom;
  - middle-drag: pan;
  - Shift + middle-drag: zoom box;
  - left double-click: identify a node in the command history;
  - right double-click: identify an element.
- Keyboard controls: Insert/Delete rotate about X, Home/End about Y, PageUp/PageDown about Z.

---

## 5. Programming commands (9.15) [CORE for the interpreter][UI]

These commands rest on the variable and macro machinery of manual 5.6 and 5.9 (summarised in spec
04, sections 5-6). The rules are repeated here where they are needed per command.

### 5.0 Variable model and substitution (recap, exact)

- A variable = **name** (case-insensitive), **list of strings** (1-based), **integer counter**
  (initially 0).
- Reference syntax in any command argument **except FOREACH's own `<Var>`**, which has no `@`:

| Expression | Counter effect | Substituted text |
|---|---|---|
| `@X` | none | counter |
| `@X+k` / `@X-k` | += k / -= k | new counter |
| `@X++` / `@X--` | +1 / -1 | new counter (pre-increment) |
| `@X=k` | = k | k |
| `@X[i]` | none | i-th element (1-based) |
| `@X[#]` | none | element at the current index of the FOREACH loop **over X** |
| `#` (with postfix operators) | - | loop counter (1-based) |

- Evaluate left to right, once per execution of a command.
- Macro `$k$` placeholders are substituted **before** `@` processing (spec 04, 5.3).

### 5.1 ADDRND (9.15.1)

- **Syntax:** `ADDRND,<var>,<numsamples>,<dist>,<param1>,...,<paramN>`
- Appends random numbers to the end of `var`, using the same algorithms as RND. It does **not**
  delete the current list.
- `numsamples` = "number of values the variable will have **after** this command". So it appends
  `max(0, numsamples - len(var))` values. Warn if nothing is appended (literal reading; Q-27).
- If `var` does not exist, create it (derived).
- `dist` and parameters: as RND (5.10). The counter is unchanged (derived).

### 5.2 CD (9.15.2)

- **Syntax:** `CD,<dir>`
- Changes the UI working directory. `dir` is absolute or relative to the current working
  directory.
- The directory **must exist**, otherwise error.
- All later relative paths resolve against the new CWD (INP, READSPEC, LOADMACRO, MKDIR, ...).
- The MDL command also changes the working directory (manual 9.15.8).

### 5.3 FOREACH (9.15.3)

- **Syntax:** `FOREACH,<Var>,<Command>`
- Runs **one** command once per element of `Var`. `Var` is the name without `@`.
- `<Command>` = everything after the second comma, verbatim. It may itself be FOREACH (nesting)
  or MACRO (for multi-command bodies).
- Substitution of `@`/`#` in the body is deferred to each iteration.
- Nested loops **must not reuse the same variable**: this must be an **error**, and the loop does
  not run.
- Regression test (manual 5.9): `Var,NNUM` / `Var,X,1,2,3,4,5` / `Var,Y,...` / `Var,Z,...` /
  `ForEach,Z,ForEach,Y,ForEach,X,N,@NNUM++,@X[#],@Y[#],@Z[#]` creates nodes 1-125 on a 5x5x5 grid.
  Node `25(iz-1)+5(iy-1)+ix` is at `(ix,iy,iz)`.

### 5.4 LOADMACRO (9.15.4)

- **Syntax:** `LOADMACRO,<Name>,<MACROFILE>`
- Binds a macro name to a file (full or relative path).
- The name is converted to **upper case**, maximum **50 characters**.
- The file is read and cached at load time. The manual advises against calling LOADMACRO inside
  macros, for performance.
- Reloading a name replaces its definition. A missing file is an error.

### 5.5 LOADVAR (9.15.5)

- **Syntax:** `LOADVAR,<filename>`
- Loads a variable from an ASCII file, **one item per line**. Blank lines are skipped and items
  trimmed (derived).
- The manual does not say how the variable is **named** (Q-28). Decision: the variable name is the
  file stem (`nodes.txt` -> `NODES`). Accept an optional extension `LOADVAR,<filename>,[name]`.
- The counter is reset to 0.

### 5.6 MACRO (9.15.6)

- **Syntax:** `MACRO,<Name>,<1>,<2>,...,<N>`
- Runs a loaded macro, replacing every `$k$` (substring-level, e.g. `Node$1$x.rs`) with argument k.
- Limits:
  - 3000 characters per substituted value and per resulting line;
  - no limit on the number of variables;
  - macros may call macros (add a recursion guard, e.g. depth 64).
- A macro may share a command's name. It does not override the command, because macros are reached
  only through MACRO.
- A missing argument for `$k$` is undefined in the manual. Decision: substitute an empty string
  and warn.

### 5.7 MACROLIST (9.15.7)

- **Syntax:** `MACROLIST`
- Lists every macro name in UI memory with its associated file name.

### 5.8 MKDIR (9.15.8)

- **Syntax:** `MKDIR,<dir>`
- Creates a directory, absolute or relative to the CWD.
  - In the original, the CWD defaults to the UI install directory. Python decision: the launch
    directory.
  - CD and MDL change the CWD.
- Reports an **error** in the command history if the directory was not created, and a **warning**
  if it already exists.
- Creates intermediate directories (derived).

### 5.9 REDUCESET (9.15.9)

- **Syntax:** `REDUCESET,<var>,[sorttype]`
- The manual heading is truncated: "(Not usable in this" (version?). See Q-29.
- Sorts the list ascending and **removes duplicates**. Values are stored as strings.

| sorttype | Ordering |
|---|---|
| `STRING` (default) | lexicographic |
| `INT` | integer order: elements converted before reduction, so "1" and "01" are duplicates |
| `FLOAT` | double-precision order: "1" and "1.0" are duplicates |

- After a numeric reduction the values are stored back as strings: `str(int)`, or the
  shortest round-trip `repr` for floats (derived).
- INT or FLOAT on non-numeric items may cause errors or inconsistent sorting. Decision: error, and
  leave the variable unchanged.

### 5.10 RND (9.15.10)

- **Syntax:** `RND,<var>,<numsamples>,<dist>,<param1>,...,<paramX>`
- Fills `var` with `numsamples` random values, **overwriting** its current contents. The counter
  is reset to 0 (derived).

| dist | Distribution | param1 | param2 |
|---|---|---|---|
| `UNI` | uniform, double | minimum | maximum |
| `UNIINT` (manual also writes `UNINT`) | uniform, integer, inclusive range (derived) | minimum | maximum |
| `NORM` | normal | mean | standard deviation |
| `LOGNORM` | lognormal | mean of the distribution | standard deviation of the distribution |
| `POISSON` (manual also writes `POSSION`) | Poisson, integer | mean | - |

- Accept the misspellings `UNINT` and `POSSION` as aliases.
- LOGNORM parameters are the mean and standard deviation of the lognormal variable itself
  (derived; Q-27): `s^2 = ln(1 + (sd/mean)^2)`, `mu = ln(mean) - s^2/2`.
- Store values as strings: floats via `repr`, integers via `str`.
- Generator (decision): numpy `PCG64` (or Python `random.Random`), seeded by RNDSEED, otherwise
  from OS entropy at start-up. The streams will not match the commercial code. For verifiable runs,
  always call RNDSEED first.

### 5.11 RNDSEED (9.15.11)

- **Syntax:** `RNDSEED,<seed>`
- Seeds or reseeds the generator used by RND and ADDRND. `seed` must be a **positive integer**,
  otherwise error.

### 5.12 SETVAR (9.15.12)

- **Syntax:** `SETVAR,...`
- A placeholder with **no effect** on the model or data. Its arguments are still substituted, so
  `@X=5`, `@X++` etc. update counters ("set or check variable counters").
- It echoes the substituted line (derived), so the user can check counter values.

### 5.13 SHOWVAR (9.15.13)

- **Syntax:** `SHOWVAR,<varname>`
- Shows the contents of the variable: all elements with their 1-based index, and the counter
  (derived format).
- An unknown variable gives a **"variable not defined"** error.

### 5.14 VAR (9.15.14)

- **Syntax:** `VAR,<Name>,<X1>,...,<Xn>`
- Defines or overwrites a variable. The list may be empty (`Var,NNUM`).
- The counter is set or reset to 0 on every VAR.
- Names are case-insensitive: redefining with different capitalisation overwrites.
- According to the manual, the only way to change a list is to redefine it with VAR. ADDRND, RND,
  REDUCESET and LOADVAR were added later and also modify lists.

### 5.15 VARLIST (9.15.15)

- **Syntax:** `VARLIST`
- Lists every variable name in memory. For a single-element variable, prints that value. Otherwise
  prints the **number of elements**. Use SHOWVAR to see multi-value contents.
- Format (derived): `NAME = value` or `NAME : n elements`, plus the counter.

---

## 6. Verification checklist for this section

| ID | Test | Expected |
|---|---|---|
| T-L1 | line math on lines A and B (4.1) | table in 4.1.7 |
| T-B1..B3 | BROADEN (4.1.3) | as listed |
| T-IO1 | WRITESPEC lines 1-3, then READSPEC,file,3,11,12,13 | lines 11-13 equal lines 1-3 resampled on the union grid |
| T-IO2 | WRITETH, then READTH Pair 0 | identical x and y (dt from x2 - x1) |
| T-S1 | CALCPAR on a two-hex block section (below) | values below |
| T-S2 | shell wall section (below) | values below |
| T-SH1 | SHEAR numeric example (3.11) | table in 3.11 |
| T-P1 | FOREACH 5x5x5 node grid | 125 nodes, numbering rule |
| T-P2 | `VAR,X,1.23,2.83,3`; `SETVAR,@X+5` -> 5; `SETVAR,@X++` -> 6; `SETVAR,@X=-15` -> -15; `VAR,X,9` -> 0 | counters |
| T-P3 | `VAR,S,10,9,1,9`; `REDUCESET,S` -> `1,10,9`; `REDUCESET,S,INT` -> `1,9,10` | |
| T-P4 | `RNDSEED,7`; `RND,R,1000,UNI,2,4` | all values in [2,4]; mean ~3; the same seed reproduces the stream |
| T-M1 | manual macro examples 5.6.1-5.6.3 | `N,1,13.52,15,100.25` etc. |

**T-S1 (solid section).**
- Model: two 1x1x1 hexahedra, element 1 at x in [0,1] and element 2 at x in [1,2], with
  y in [0,1] and z in [0,1].
- Plane: through (0,0,0.5) with n = (0,0,1); r = (1,0,0). So ex = X, ey = Y, ez = Z.
- Expected section properties: A = 2, C = (1, 0.5, 0.5), Ixx = 2*1^3/12 = 0.166667,
  Iyy = 1*2^3/12 = 0.666667, Ixy = 0.

| Stress state (both elements unless noted) | Expected |
|---|---|
| Szz = +10 (el. 1), -10 (el. 2) | Fz = 0, **My = +10**, Mx = Mz = 0 |
| plus Sxz = 3 | Fx = 6 |
| Syz = +1 (el. 1), -1 (el. 2) | Mz = -1 |

**T-S2 (shell wall).**
- Model: two vertical shells in the plane y = 0, each 2 wide (x in [0,2] and [2,4]) and 2 high
  (z in [0,2]), thickness 0.5, local x' = X and y' = Z.
- Plane: z = 1, n = (0,0,1), r = (1,0,0).
- Expected section properties: A = 2, Iyy = 0.5*4^3/12 = 2.6667.

| Stress state | Expected |
|---|---|
| Sy'y' = +100 (el. 1), -100 (el. 2) | Fz = 0, My = +200 |
| Sy'y' = +100 in both | Fz = 200 |

---

## 7. Open questions / ambiguities

| ID | Topic | Issue | Decision taken in this spec |
|---|---|---|---|
| Q-1 | RUNSITE | 6.4.5 refers to a RUNSITE command that 9.12 does not list | Implement `RUNSITE,[model]` with the RUN semantics |
| Q-2 | RUN execution | Synchronous or asynchronous is not stated. The UI only "opens a tab with the module output" | Synchronous in `.pre`/macro execution, with streamed output |
| Q-3 | COMBIN input | Extension and AFWRITE input of COMBIN not stated (it needs FILE81/FILE82) | See spec 04 / 05 |
| Q-4 | `.ess` layout | Single-frame `.ess` format not given | Assume the ELEMENT_CENTER_ABS_MAX_STRESSES.TXT layout, signed |
| Q-5 | Section integration | No formulas in the manual (piece geometry, sign convention, components) | Section 3.2 algorithm. Verify against Demo 8 / Verification Problem 47 |
| Q-6 | Local axes | "right" = local +X. The normal's role (local Z?) is not stated | ez = n, ex = projected r, ey = ez x ex |
| Q-7 | Output component order and names | CALCPAR parameter names/order and CSV column order not given | Table 3.8; CSV `Fx,Fy,Fz,Mx,My,Mz` |
| Q-8 | CALCSECTHIST output | Format not described (only CALCSECTHISTDB's is) | Same 7-column CSV with MAX row |
| Q-9 | Time column | Whether step k maps to (k-1)ts or k*ts | (k-1)ts (frame 1 is time 0.000 in the frame naming) |
| Q-10 | Manual example 5.8.2 | CUTVOL x-range [52.5, 52.8] with the CSECT plane at x = 55.6 looks inconsistent | Use it as a workflow pattern only |
| Q-11 | Whole-model base forces | The manual says the calculation commands work on non-cross-section models, but the force integration needs a plane | Use a cut holding the whole structure plus CSECT at the base |
| Q-12 | Shell plate bending in section cuts | Not stated whether shell bending moments (Mx'x', My'y', Mx'y') enter the resultants, nor their sign convention | Include them by default (option), with M = integral z*sigma dz about the shell normal; verify |
| Q-13 | Coincident faces | Not addressed | Count once; prefer the element on the -n side |
| Q-14 | SHEAR arguments 6-7 | Syntax line `[Fvw],[Fbe]` vs description `<Abe>` (area) and `<Fybe>` (yield stress). BBCGEN describes Fvw/Fbe as forces, while the SHEAR text computes F_VW from rho, A_W, f_y | Args 6, 7 = A_BE, f_y,BE; F_VW computed. Offer a flag to read 6, 7 as forces F_VW, F_BE directly |
| Q-15 | SHEAR geometry | How h_W, l_W, t_W, A_W are derived from a panel group is not stated | Bounding extents in the panel plane, with Z vertical; A_W = l_W t_W |
| Q-16 | alpha_c | Not defined (and "Chapter 11" is cited while the form matches ACI 318-08 eq. 21-7) | ACI 318-08 21.9.4.1 interpolation |
| Q-17 | SHEAR columns | "six columns ... panel number, upper bound of ACI and Wood, lower bound of Wood, Barda, G-W" can be read several ways | Panel, ACI(capped), Wood(capped), Wood LB, Barda, G-W. Also print raw values |
| Q-18 | SHEAR sign of N_U | Not stated | Compression positive (raises strength, as in Barda/G-W) |
| Q-19 | BROADEN Smooth1 ("peak bridging", dialog "Peak Difference(%)") | Criterion not defined: valley depth relative to the lower peak, or amplitude difference between adjacent peaks | Valley-depth criterion by default. Alternative selectable. Test T-B2 with Smooth1 = 10 tells them apart |
| Q-20 | LBINCORS | No algorithm, no parameters | Parser only; transformation pluggable; raise "not specified" |
| Q-21 | BROADEN grid | The original likely works only on the discrete points (hence the 301-point warning) | Augmented grid by default; "ACS-compatible" discrete option |
| Q-22 | READTH Pair 0 / WRITETH formats | "ACS SASSI output format" single column: header lines? dt on the first line? values per line? | First numeric = dt, then values; t starts at 0 |
| Q-23 | READSPEC headers | Header-line handling of RS/TFU/TFI files not stated | Skip non-numeric lines |
| Q-24 | Text arguments with commas | Titles/labels containing commas | Allow double-quoted arguments (extension) |
| Q-25 | SHOWDOF off | No documented way to turn the markers off except the window | `SHOWDOF,NONE` extension; unchecking all boxes in the window |
| Q-26 | COLOR numbering | Index base of `<Num>` and the exact palette spelling | 1-based; case-insensitive names from 6.5.6 |
| Q-27 | ADDRND numsamples / LOGNORM parameters | "values the variable will have after" (total) vs "number appended"; lognormal mean/sd of X or of ln X | Total count; mean/sd of X |
| Q-28 | LOADVAR variable name | The syntax has no name argument | Name = file stem, optional extension argument |
| Q-29 | REDUCESET availability | Heading says "(Not usable in this ..." (truncated) | Implement it anyway; mark experimental |
| Q-30 | Spec 04 vs this file | Spec 04 section 8.3 sketched local axes with e1 = n. This file uses ex = r, ez = n, following the manual's "r = local +X" | This file is authoritative for 9.13. Update spec 04 |
| Q-31 | CNGVIEW units | Rotation units and the meaning of the zoom constant | Degrees; zoom = camera scale factor (1 = default fit) |
| Q-32 | Partial animation arguments | "Partial argument entry may cause inconsistent results" | Fill from SASSIani.xml defaults and warn |
| Q-33 | CALCM mass of excavated soil, MUNITS | Whether excavation elements count, and the weight/mass units of MT/MR | Report excavation separately; convert by MUNITS |
| Q-34 | Error vs ignore for undefined lines in line math | Not stated (only the plotting commands say "ignored") | Error for line math; ignore for SPECPLOT/THPLOT |
