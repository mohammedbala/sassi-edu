# Spec 01 — Program Capabilities, Architecture, Element Library and Modeling Rules

**Manual source:** ACS SASSI Version 3 (IKTR9) User Manual, Chapter 1 "Introduction" (Sections 1, 1.1 to 1.7).
Extracted text lines 89 to 1472 of `reference/acs-sassi.txt`. Printed pages 1 to 32 (PDF pages about 3 to 36).
Figures 1.1 to 1.5 and the CHECK window screenshot were read from the PDF images and are transcribed below.

**Purpose of this document.** Chapter 1 sets out the program's scope. It covers the module architecture, the run sequences, the file flow between modules, size limits, the finite element (FE) library, the modeling rules and checks, loading, the nonlinear iteration schemes, the interpolation options, units, the command input syntax and data checking. Later chapters give the per-module input details. This spec gives a developer the frame and the rules every other spec has to fit into.

**Tags.**
- **[CORE]** Needed for computational correctness.
- **[IO]** File or input format, file naming, or data exchange between modules.
- **[UI]** User interface, plotting or convenience.
- **[ADV]** Advanced options: incoherency, nonlinear (Option NON), Option A/AA, Option PRO, FFV and similar.
- **[INFO]** Background only. Nothing to implement.
- **(Impl. note)** Guidance this spec adds where the manual is silent. Each one is listed again under Open Questions.

---

## 0. Executive summary for implementers

1. ACS SASSI is a frequency-domain substructuring program for linear soil-structure interaction (SSI) analysis, following the SASSI Flexible Volume family of methods. Nonlinear behavior is handled by iterative equivalent linearization: each iteration is a linear complex-frequency solution, followed by a time-domain post-processing step that updates the properties. **[CORE]**
2. The program is a chain of separate **modules**: EQUAKE, SOIL, SITE, POINT2/POINT3, FORCE, HOUSE, ANALYS, COMBIN, MOTION, RELDISP and STRESS. Option modules add NONLINEAR, COMB_XYZ_THD, COMB_XYZ_STRAIN, LOADGEN and SSI2ANSYS. Modules exchange data through numbered files (FILE1, FILE2, FILE3, FILE4, FILE8, FILE9, FILE73, FILE74, FILE77, FILE78, FILE88 and others) and through text result files (`.TFU`, `.TFI`, `.ACC`, `.RS`, `.THD`, `.TFD`, `.THS`). A **UI** layer holds the model database, writes each module's input file (AFWRITE), runs CHECK, runs the modules and post-processes the results. **[CORE][IO]**
3. Linear run order: **SITE → POINT → HOUSE → ANALYS → MOTION → RELDISP → STRESS**. EQUAKE (spectrum-compatible motion generation) and SOIL (equivalent-linear or nonlinear site response) are optional and run before SITE. FORCE replaces the seismic load for external-force problems. COMBIN merges FILE8 results from two ANALYS runs. **[CORE]**
4. The element library has SOLID (8 nodes), BEAM (3 nodes), SHELL (4 nodes, Kirchhoff), TSHELL (4 nodes, Mindlin), PLANE (4 nodes, plane strain), SPRING (2 nodes) and General Matrix (GM, 2 or 3 nodes). The excavated soil is modeled with SOLID (3D) or PLANE (2D) elements. **[CORE]**
5. Interaction nodes have translational DOFs only and must lie on soil layer interfaces. The substructuring method (FV, FI-FSIN/SM, FI-EVBN/MSM, FFV) is chosen by which nodes are made interaction nodes. **[CORE]**
6. Time-domain inputs go through the FFT. The number of points must be a power of 2, the time step must be uniform, and the record needs a trailing "quiet zone" of zeros. **[CORE]**
7. Damping is hysteretic, applied as a complex modulus (frequency independent). The mass matrix is 50% lumped + 50% consistent, except BEAM (consistent) and plates/shells (lumped). **[CORE]**
8. Transfer functions are computed at a sparse set of SSI frequencies and interpolated in the complex domain onto the FFT grid. There are 7 interpolation options, 0 to 6. Option 6 is a complex bicubic spline. **[CORE]**
9. Input is free-format commands, `Keyword, p1, p2, ..., pn`, separated by commas and case-insensitive. They can be typed, generated from menus, or read from `.pre` files with `INP,<fname>`. CHECK validates the model, and AFWRITE runs CHECK before writing any module input file. **[IO][UI]**

---

## 1. Product overview and option set (Manual §1, printed pp. 1–5)

### 1.1 Scope statements **[INFO]/[CORE]**

| Item | Specification |
|---|---|
| Analysis class | 3D (and 2D) linear and nonlinear SSI for shallow, embedded, deeply embedded and buried structures. Seismic input may be coherent or incoherent. |
| Methodology | The linearized SASSI complex-frequency substructuring methodology. It is extended to nonlinear soil (equivalent-linear iterations), nonlinear structures (Option NON) and probabilistic analysis (Option PRO, separate manual and out of scope here). |
| Run modes | **Interactive** for a single SSI model. **Batch** for one or many models. Implement both: a GUI/CLI session, and a batch runner that executes a module sequence from a script. **[UI]** |
| Translators | (a) ANSYS → ACS SASSI: the manual calls the ANSYS file a "CBD file". This is presumably ANSYS `.cdb`, see Open Questions. (b) University SASSI2000 fixed-format input files → ACS SASSI. (c) ACS SASSI → ANSYS as APDL input. **[IO]** Very relevant to this user, who works in ANSYS. |
| Resource management | The UI manages all files, directories and links between modules automatically. **[UI]** |
| Demos | 12 interactive demos (Demo 1 to Demo 12), meant to be worked in ascending order. They cover Options A, AA and NON, post-processing macros, frequency selection and animation frame selection. Demo 3 covers automatic frequency addition, Demo 9 the nonlinear RC shearwall building at 0.60g, and Demo 10 the base-isolated building. **[INFO]** These suggest tutorial and regression cases for our implementation. |

### 1.2 Multiple load cases in a single run **[CORE]**

- **Coherent seismic analysis.** The X, Y and Z input directions are solved together in memory as three right-hand sides of the same factorized system at each frequency. No restart files are needed to get the three directions.
- **External force / impedance problems.** Up to **500 separate load cases** in one run, limited by RAM. Typical use is the many unit harmonic loads needed to build the impedance of a flexible foundation.
- **Incoherent seismic analysis.** Up to **50 stochastic simulations** in one run, each with X, Y and Z components.
- **(Impl. note)** Build the ANALYS solver around a per-frequency factorization of the complex dynamic stiffness, reused for many right-hand sides. The limits 500 and 50 are the manual's practical RAM limits. In Python they can be configurable soft limits.

### 1.3 Restart analysis types **[CORE]**

Restart files are optional for coherent 3-direction seismic runs. They remain useful for two kinds of rerun:

- **"New Structure"** restart. The soil (site and impedance) is unchanged and the structure or near-field properties change. This is used for nonlinear soil iterations and Option NON iterations. Reusing the soil impedance matrix already computed in the initiation run cuts the time per iteration by a factor of 2 to 5, depending on the embedment size.
- **"New Seismic Environment"** restart. The structure is unchanged and the seismic load vector changes, as in incoherent simulations.

Restart file names (from Figure 1.1, with details in the Chapter 3 spec) are `COOXxxx` and `COOTKxxx`, where `xxx` is the 3-digit frequency order number, for example `COOX001`. The other restart files are the index files `COOXI` and `COOTKI`, plus `DOFSMAP`, `FILE90` and `FILE91`. **[IO]**

### 1.4 Option modules **[ADV]**

| Option | What it does | Modules | Notes |
|---|---|---|---|
| **Option A** ("ANSYS") | Two-step approach. Step 1 is the global SSI or SSSI analysis in ACS SASSI. Step 2 is a refined ANSYS stress analysis, quasi-static or dynamic, with boundary conditions taken from the SSI responses: motions at all time steps or at selected critical steps. Step 2 can also compute soil pressure on basement walls and slabs, including soil plasticity, separation and sliding, using refined ANSYS soil models. | **LOADGEN** (exports loads and BCs to the ANSYS model) | Assumes that nonlinearity in the ANSYS model does **not** change the SSI soil motion at the foundation-soil interface. Follows ASCE 4-16 §4.6. Has a separate manual. |
| **Option AA** ("Advanced ANSYS") | Runs the SSI analysis directly with the ANSYS FE model's dynamic matrices, with no conversion to ACS SASSI elements. Only the model **topology** is transferred, not materials or real constants. Supports advanced ANSYS element types, pipes, shells with shear flexibility, coupled nodes, constraint equations, MPC elements, FLUID80 and MATRIX50 super-elements. | **SSI2ANSYS** (imports the ANSYS FE model). HOUSE and ANALYS are replaced by **HOUSEFSA** and **ANALYSFSA**. The manual also writes "ANALYSFA", treated here as a typo. | **WARNING:** for Option AA runs, HOUSEFS and ANALYSFS must be swapped for HOUSEFSA and ANALYSFSA. In the UI this is done with **Modules/Location**, which changes the executable paths. MATRIX50 super-elements can either stay inside the Option AA ANSYS model or be converted automatically to ACS SASSI **GM** elements. Option AA always uses the HOUSE node-renumbering optimizer. To get stresses, use Option A to transfer the motions to ANSYS. |
| **Option NON** | Nonlinear structure SSI by piecewise (iterative) equivalent linearization: a hybrid time / complex-frequency approach. Covers RC shearwall cracking and post-cracking (in-plane shear or in-plane bending), base-isolator hysteresis (nonlinear springs), pile-soil interface slip and potential foundation sliding. | **NONLINEAR**, **COMB_XYZ_THD** | Detailed in §8 below. Validated for low-rise RC shearwalls dominated by in-plane shear. **Out-of-plane** concrete nonlinearity is **not** included. Nonlinear springs are translational only (X, Y, Z). **No nonlinear rotational springs.** |
| **Option PRO** | Probabilistic site response (PSRA) and probabilistic SSI (PSSIA) by Latin Hypercube Sampling (LHS). Randomizes the response-spectrum shape, the Vs and D profiles, the G/Gmax and D curves, and the effective structural stiffness and damping per element group. | separate | **Not described in this manual. Out of scope** for this re-implementation unless added later. |

Fast solver: the standard SSI modules HOUSE and ANALYS exist as fast-solver executables **HOUSEFS** and **ANALYSFS**. All other modules are the same. **[INFO]** In Python there is a single implementation of each, and these names are kept only as aliases in the module-location table.

---

## 2. Basic capabilities, Manual §1.1 items (i) to (xiii) → feature requirements

Each row is a feature our program must provide or explicitly defer. The "Owner" column names the module or spec where the details live.

| # | Capability | Key requirements and nomenclature | Tag | Owner |
|---|---|---|---|---|
| i | **Spectrum-compatible input generation** | Generate three-component acceleration histories that match a target design ground response spectrum, **with or without time-varying correlation** between components. Alternatively generate them by complex **Fourier phasing of selected "seed records"** (ASCE 4-16 term). Provide **baseline correction**. Compute the **PSD** and the peak ground acceleration, velocity and displacement, so the user can check the USNRC SRP 3.7.1 requirements. | [CORE] | EQUAKE (Ch. 3/6 spec) |
| ii | **Motion incoherency and wave passage** | Both isotropic (radial) and anisotropic (directional) coherency models. **Six incoherent SSI approaches:** two simplified deterministic ones, **AS** and **SRSS** (EPRI-benchmarked), three alternate deterministic ones, and the rigorous stochastic **"Simulation Mean"** approach. **Seven plane-wave coherency models:** Luco-Wong 1986 (theoretical, not validated), Abrahamson 1993, 2005 (all sites, surface foundations), 2006 (all sites, embedded foundations), 2007a (rock sites, all foundations), 2007b (soil sites, surface foundations), plus **user-defined** models. Directional (anisotropic) Abrahamson variants are included. User-defined models may differ between the two horizontal principal directions, for example from 2D site response along the maximum and minimum layer-slope directions. | [ADV] | HOUSE (random field decomposition, FILE77), ANALYS |
| ii-W | **WARNING (SRSS):** SRSS is implemented **for benchmarking only**. It was validated by EPRI for stick models with rigid basemats. It has no clear convergence criterion for the number of incoherent spatial modes, which for flexible foundations can run to tens or hundreds. It can be over-conservative in the mid-frequency range (sometimes above the coherent result) and unconservative at high frequencies. It is not recommended for design-basis FE models with elastic foundations. AS and SRSS should be used only for rigid-basemat stick models. For FE models with flexible foundations the stochastic simulation approach is recommended, because it does not intrude on the SSI system dynamics. AS and SRSS on flexible foundations can be badly biased, especially in the vertical direction (as also noted in ASCE 4-16). | [ADV][UI] | Show this as a UI warning whenever SRSS is selected |
| iii | **Nonlinear hysteretic soil** | Seed-Idriss iterative equivalent-linear procedure for **global** soil nonlinearity (free-field wave propagation, done in SOIL with the SHAKE methodology) and **local** soil nonlinearity (SSI effects, using near-field PLANE or SOLID elements in the HOUSE model). Uses fast SSI restart (reusing the impedance from the initiation run, 2 to 5 times faster per iteration). In batch mode the simultaneous X, Y and Z effect at each iteration comes from the **COMB_XYZ_STRAIN** auxiliary program. | [CORE][ADV] | SOIL, HOUSE, STRESS (§8.1) |
| iv | **Nonlinear hysteretic structure (Option NON)** | Iterative equivalent linearization for RC wall cracking and post-cracking (shell elements), rubber-bearing base isolators (nonlinear springs), and local pile-soil interface slip in the vertical direction (nonlinear springs). | [ADV] | NONLINEAR (§8.2) |
| v | **Nonuniform soil motion / multiple seismic input** | Variable-amplitude seismic input. For continuous foundations the free-field complex amplitude varies in the horizontal plane following frequency-dependent spectral patterns, given by **complex amplification factors (relative transfer functions)** relative to a reference motion. Can be combined with incoherency and wave passage. Multiple support excitation of isolated foundations is also supported, with up to 5,000 foundation zones. See the AMP command (`AMP,<no>,<a1>...<a100>`) in the Chapter 9 spec. | [ADV] | HOUSE / ANALYS |
| vi | **Seven complex interpolation schemes** | Options 0 to 6, see §9. Option 6 (bicubic complex spline) is recommended for incoherent analysis with a dense frequency set (more than 200 frequencies). Comparative plots of computed TF (`.TFU`) against interpolated TF (`.TFI`) must be available. | [CORE][UI] | MOTION, STRESS |
| vii-a | **Fast Flexible Volume (FFV)** | For deeply embedded structures (DES, such as SMRs). The interaction nodes are the outer excavation surface nodes **plus internal horizontal node layers** inside the excavation. Up to tens of times faster than FV. Generated with **INTGEN**. | [ADV][CORE] | UI INTGEN, ANALYS |
| vii-b | **Automatic selection of additional SSI frequencies** | Add frequencies where the interpolated ATF or STF shows peaks not captured by the computed frequencies. Done with UI commands and macros (Demo 3). | [UI] | UI macros, COMBIN |
| viii | **Complex TF vector visualization** | At a selected SSI frequency, animate colored vectors at the nodes: **red = X, green = Y, blue = Z**. Vector length is the TF amplitude and vector orientation is the TF phase. | [UI] | Plot spec |
| ix | **Amplitude TF / spectral acceleration field plots** | At a selected SSI frequency or damping value, show amplitude TF or spectral acceleration over the whole model as an animated deformed shape (with controllable frame speed, so it can also be shown static) or as bubble plots. ZPA or spectral amplitude at resonant frequencies can be shown as a deformed shape. | [UI] | Plot spec |
| x | **Time-history deformed shapes** | Acceleration and relative-displacement histories shown as static deformed shapes (at a selected time or at the maximum) or as animations. | [UI] | Plot spec |
| xi | **Nodal stress contours** | Average nodal stresses (all **6 components in global coordinates**) from **element-center** stresses of SHELL and SOLID elements. Node stress is the average of the center stresses of the adjacent elements, with **no shape-function extrapolation**. Both maximum and time-varying values. Static (max) or animated contours with automatic frame selection. Maximum element-center stresses are also written to a text file. | [CORE][UI] | STRESS + UI |
| xii | **Seismic soil pressure on walls** | Uses near-field SOLID elements. Nodal pressure is the average of the adjacent element-center pressures. Gives maximum and time-varying values. Can be **combined with the static soil bearing pressure** and plotted. Static or animated contours with automatic frame selection. | [CORE][UI] | STRESS + UI |
| xiii | **Post-SSI combination** | Time histories (acceleration, displacement, stress) are combined by **algebraic summation** of co-directional effects. In-structure response spectra (ISRS) are combined by (i) **weighted linear combination** or (ii) **SRSS** of the X, Y and Z co-directional effects. Also provides the **envelope**, **broadening** and **average** of multiple spectral curves, and maximum stresses, forces, moments and soil pressures, with or without the static bearing component. | [UI][CORE] | UI calculation commands (Ch. 9.13) |

**QA / V&V context [INFO].** The NQA version follows 10 CFR50 App. B, 10 CFR21 and ASME NQA-1 2008/2009a. IKTR9 ships **56 V&V problems**, benchmarked against analytical solutions, SHAKE91, SASSI2000 and ANSYS. Problem 49 checks SOIL-NON against DEEPSOIL (Ch. 3). Problem 51 checks Option NON for a base-isolated building against ANSYS V15. The RC nonlinear building is checked against PERFORM3D. **(Impl. note)** Our verification suite should mirror these categories: shallow, embedded and buried foundations; rigid and flexible foundations; piles; surface and body waves; incoherency; wave passage; multiple support; linear and nonlinear.

---

## 3. Model size restrictions (Manual §1.2) **[CORE][IO]**

### 3.1 Global limits

| Build | Max nodes | Max equations / DOF |
|---|---|---|
| IKTR9 | 100,000 nodes. HOUSE lists 99,999. | n/a |
| IKTR9_650K (also written IKTR8_650K) | 650,000 nodes | 2,500,000 equations |

Tested with up to 625,000 nodes and 35,000 interaction nodes on 128 to 512 GB RAM workstations.

Hardware guidance [INFO]:

| Machine | Supported problem |
|---|---|
| 16 GB RAM | Up to 100,000 nodes with up to 8,000 interaction nodes runs efficiently. |
| 32 to 192 GB | More than 10,000 and up to 20,000 interaction nodes. |
| 512 GB | More than 15,000 to 20,000 interaction nodes. |

Practical recommendation: keep the number of interaction nodes **below 20,000**. If needed, use FI (MSM/ESM) or FFV instead of FV.

**Model-size effect on file formats [IO]:**
- For models with **more than 99,999 nodes**, the integer fields in the **fixed-format** module input files written by **AFWRITE** are widened. This applies to the FORCE, HOUSE, MOTION, RELDISP and STRESS inputs. The UI `.pre` file is **not** affected by model size.
- Result file names: node numbers use **5 digits** for models below 100,000 nodes and **6 digits** above (Ch. 3 cross-reference, for example `xxxxxTR_y.TFU`).

### 3.2 Per-module limits (enforced as validation checks in a "compatibility" profile)

| Module | Quantity | Limit |
|---|---|---|
| EQUAKE | Time steps of simulated acceleration histories | 32,768 |
| SOIL | Time steps of simulated acceleration histories | 32,768 |
| SOIL | Soil material curves | 100 |
| SOIL | Data points per soil curve | 11 |
| SOIL | Soil layers | 200 |
| SITE | Soil layers | 200 |
| SITE | Half-space layers | 20 |
| SITE | Analysis frequencies | 500 |
| POINT | Soil layers | 200 |
| POINT | Half-space layers | 20 |
| POINT | Analysis frequencies | 500 |
| POINT | Embedment layers | 50 |
| FORCE | Analysis frequencies | 500 |
| HOUSE | Nodes | 99,999 (IKTR9) / 650,000 nodes or 2,500,000 DOF (650K build) |
| HOUSE | Interaction nodes | 100,000 (impractical: needs about 4 to 5 TB RAM) |
| HOUSE | Materials or cross-section geometries | 9,999 |
| HOUSE | Analysis frequencies | 500 |
| HOUSE | Embedment layers | 50 |
| HOUSE | Multiple supports (foundation zones) | 5,000 |
| HOUSE | Interaction nodes per embedment layer for incoherent analysis | "20,0000" as printed, read as 20,000. See Open Questions. |
| ANALYS | Analysis frequencies | 500 |
| MOTION | Analysis frequencies | 1,500 |
| MOTION | Time steps / Fourier frequencies | 65,536 |
| MOTION | Damping values for response spectra | 5 |
| RELDISP | Analysis frequencies | 1,500 |
| RELDISP | Time steps / Fourier frequencies | 65,536 |
| STRESS | Analysis frequencies | 1,500 |
| STRESS | Time steps / Fourier frequencies | 65,536 |
| STRESS | Elements per group | 100,000 |
| STRESS | Element groups per model | 10,000 |
| NONLINEAR | Nodes | 300,000 |
| NONLINEAR | Data points per backbone curve | 100 |

**(Impl. note)** A Python implementation does not need fixed array sizes. Provide a `limits` profile with `"acs_sassi_v3_iktr9"` as the default, which **warns** (does not abort) when these limits are exceeded, and an `"unlimited"` profile. Two limits drive physics and numerics and should stay enforced as **errors**, because they keep results comparable with the reference program:
- the power-of-2 FFT length;
- consistency of the analysis frequency count across SITE, POINT, HOUSE and ANALYS.

Note the inconsistency in the maximum FFT points: §1.5.2 says 32,768, while the MOTION, RELDISP and STRESS limits say 65,536. See Open Questions.

---

## 4. Modular structure configuration (Manual §1.3, Figure 1.1)

### 4.1 Module list **[CORE]**

The main SSI modules ("eleven" in the manual) are **EQUAKE, SOIL, SITE, POINT (POINT2 / POINT3), FORCE, HOUSE, ANALYS, MOTION, STRESS, RELDISP, COMBIN**. POINT2 is for 2D SSI and POINT3 for 3D SSI.

| Module | Role (one line) | Input file ext. (written by AFWRITE) | Key outputs |
|---|---|---|---|
| EQUAKE | Simulate spectrum-compatible 3-component acceleration histories (optional) | `.equ` | Acc files (text; `.acc/.vel/.dis`) |
| SOIL | 1D site response: SHAKE-type equivalent linear (SOIL-EQL), or time-domain hyperbolic nonlinear (SOIL-NON, DEEPSOIL-like, no water table) | `.soi` | FILE88 (iterated soil properties → SITE), FILE73 (soil curves → STRESS), `.TH` files |
| SITE | Free-field site response for layered soil. Computes the free-field motion information (FILE1) and transmitting-boundary data (FILE2). Adds half-space layers and the viscous boundary. Holds the control point and wave composition. | `.sit` | FILE1, FILE2 |
| POINT2/POINT3 | Information to form the frequency-dependent soil flexibility matrix at the interaction nodes, using the axisymmetric point-load model with a user "radius" | `.poi` | FILE3 (needs FILE2) |
| HOUSE | Element mass and stiffness matrices of the structure and near-field soil, plus the excavated soil. Also does incoherency random-field decomposition, multiple-support data and nonlinear-soil property updates. | `.hou` | FILE4 (`.n4`), COOSK/COOSM (fast solver), FILE77, FILE78 |
| FORCE | Load vector for external force cases (not used for seismic, except foundation impedance) | `.frc` | FILE9 (per case: FILE9001 … FILE9500) |
| ANALYS | At each SSI frequency, forms the soil flexibility, inverts it to the impedance at the interaction nodes, forms the seismic or force load vector, and solves for the complex transfer functions | `.anl` | FILE8; restart files; FILE11 (global impedance option) |
| COMBIN | Merges two FILE8 solutions computed at different frequency sets (inputs renamed FILE81 and FILE82) | (none listed) | new FILE8 |
| MOTION | Interpolates the TFs in the complex frequency domain, then computes the response acceleration, velocity and displacement histories and ISRS at the requested nodes | `.mot` | `.TFU`, `.TFI`, `.ACC`, `.RS`, FILE12, FILE13, complex TFI (for RELDISP) |
| RELDISP | Relative displacements from the interpolated complex ATF | `.rdi` | `.TFD`, `.THD`, max relative displacements |
| STRESS | Element stress, strain and force histories and peaks. Nonlinear-soil strain evaluation. | `.str` | `.TFU`, `.TFI`, `.THS`, FILE14, FILE15, FILE74 |

Option modules [ADV]: **LOADGEN**, **SSI2ANSYS** (Options A/AA), **NONLINEAR**, **COMB_XYZ_THD** (Option NON). Auxiliary programs: **COMB_XYZ_STRAIN** (nonlinear soil), **Remove_Frequencies_from_FILE8.exe** (§6, item 11).

Name spellings in the manual vary: `COMBIN_XYZ_THD`, `COMBINE_XYZ_THD`, `COMB_XYZ_THD`, `COMB_XYZ`; also `COMB_XYZ_STRAIN` and `COMBIN_XYZ_STRAIN.exe`. **(Impl. note)** Use `COMB_XYZ_THD` and `COMB_XYZ_STRAIN` as canonical and accept the others as aliases.

UI module [UI]: the "ACS SASSI UI". It builds and checks models, runs analyses, and reviews and post-processes results, including animations. It has a parametric macro language. It contains all the old MAIN, PREP and SUBMODELER commands.

### 4.2 Run sequences **[CORE]**

```
LINEAR SSI (seismic):
  [EQUAKE] -> [SOIL] -> SITE -> POINT(2|3) -> HOUSE -> ANALYS -> [COMBIN] -> MOTION -> RELDISP -> STRESS
LINEAR SSI (external forces / impedance):
  SITE -> POINT -> HOUSE -> FORCE -> ANALYS -> MOTION/RELDISP/STRESS
NONLINEAR SOIL (equivalent-linear SSI iterations, local soil nonlinearity):
  initiation: SITE -> POINT -> HOUSE -> ANALYS(save restart) -> STRESS
  loop k=2..: HOUSE -> ANALYS(restart "New Structure") -> STRESS [-> COMB_XYZ_STRAIN]   until converged
NONLINEAR STRUCTURE (Option NON):
  iteration 1: SITE -> POINT -> HOUSE -> ANALYS -> MOTION -> RELDISP -> COMB_XYZ_THD -> NONLINEAR
  loop k=2..: HOUSE -> ANALYS(restart) -> MOTION -> RELDISP -> COMB_XYZ_THD -> NONLINEAR  until converged
```

Batch convention (Ch. 3 cross-reference) **[IO]**:
- Each module runs as `MODULE.exe < MODULE.inp`.
- The `MODULE.inp` file has exactly three lines: `modelname`, `modelname.<ext_input>`, `modelname_<MODULE>.out`.
- **(Impl. note)** Support the same three-line driver file in the Python CLI, for example `python -m acssassi.site < SITE.inp`, so the original batch scripts remain readable.

Iteration flag files (Ch. 3 cross-reference) **[IO][ADV]**:
- HOUSE writes a `.liq` file. If it holds `1`, the next HOUSE run reads FILE74 for nonlinear-soil materials.
- `PANEL.NON` or `SPRING.NON` plays the same role for NONLINEAR. If the file is missing, NONLINEAR treats the run as the linear 1st iteration.
- The `.pin` file holds the nonlinear soil inputs for HOUSE.

Constraints:
- The seismic option and the external force option **cannot be combined in one run**. Run them separately and superpose linearly (§7).
- HOUSE can run without SITE and POINT, **except for embedded models**. HOUSE still needs the `.sit` file in the working directory (Ch. 3 cross-reference).
- SITE must run before POINT, because POINT reads FILE2.

### 4.3 Figure 1.1 data-flow (transcribed) **[IO]**

Figure legend: red outline = **text** file, green outline = **binary** file. Blue module names belong to the options.

| From | File / data | To | Figure colour |
|---|---|---|---|
| EQUAKE | "Acc file" (acceleration history) | SOIL | red (text) |
| SOIL | **FILE88** (equivalent soil properties) | SITE | red (text) |
| SOIL | **FILE73** (soil material curves) | STRESS | green in figure. §1.5.4 calls it a text file, see Open Questions. |
| SITE | **FILE1** (free-field motion data) | ANALYS | green |
| SITE | **FILE2** (transmitting boundary data) | POINT | green |
| POINT | **FILE3** (flexibility data) | ANALYS | green |
| HOUSE | **FILE4**, **FILE77** (structure matrices; incoherency / multiple support) | ANALYS. FILE4 also goes to STRESS (Ch. 3). | green |
| FORCE | **FILE9** (load vectors) | ANALYS | green |
| ANALYS ↔ ANALYS | Restart files **COOXxxx, COOTKxxx, DOFSMAP, FILE90, FILE91** (written on initiation, read on restart) | ANALYS | (green text in figure) |
| ANALYS | **FILE8** (complex TF solution) | COMBIN → MOTION, STRESS (and RELDISP) | green |
| MOTION | `*.tfu, *.tfi, *.acc, *.rs` | Results / LOADGEN | red |
| MOTION | **FILE12, FILE13** | Results | red |
| MOTION | "Complex TFI" (interpolated complex ATF) | RELDISP | red |
| RELDISP | `*.thd, *.tfd` | Results, COMB_XYZ → NONLINEAR, LOADGEN | red |
| STRESS | `*.tfu, *.tfi, *.ths` | Results / LOADGEN | red |
| STRESS | **FILE14, FILE15** | Results | red |
| HOUSE | **FILE78** | STRESS (nonlinear soil) | red |
| STRESS | **FILE74** (effective strains) | HOUSE (next nonlinear-soil iteration) | red |
| COMB_XYZ + NONLINEAR (Option NON) | updated HOUSE model (`Modelname_new.hou`) | HOUSE | — |
| ANSYS model ("Import FE Model") | → SSI2ANSYS (Option AA) | HOUSE | — |
| Results (all) | → LOADGEN (Option A) | "Export Load" → ANSYS model | — |

ASCII rendering of Figure 1.1:

```
 INPUT                    ANALYSIS                                  RESULTS
 FORCE --(9)-----------\
 EQUAKE                 \        Restart: COOXxxx, COOTKxxx,
   |(Acc file)           \       DOFSMAP, FILE90, FILE91 (in/out)
 SOIL --(73)--------------\---------------------------\
   |(88)                   \                           \
 SITE --(1)-------------> ANALYS --(8)--> COMBIN --+--> MOTION --> *.tfu *.tfi *.acc *.rs ; (12,13)
   |(2)                   /                         |      | (Complex TFI)
 POINT --(3)-----------> /                          |    RELDISP --> *.thd *.tfd --> COMB_XYZ -> NONLINEAR --> HOUSE
 HOUSE --(4,77)-------> /                           +--> STRESS --> *.tfu *.tfi *.ths ; (14,15)
   ^  \---(78)----------------------------------------> STRESS
   |<-----------------------------------------(74)----- STRESS
   |<-- SSI2ANSYS <-- ANSYS Model (Import FE Model)          [Option AA]
 Results --> LOADGEN --> ANSYS Model (Export Load)            [Option A]
```

**(Impl. note)** Keep the numbered file names (FILE1, FILE2, … FILE91) as the on-disk names of the inter-module data, for nomenclature fidelity. Their internal format can be our own (HDF5 or NPZ), documented in the IO specs. Text result files must follow the manual's naming and column formats, which are defined in the Chapter 3 spec (Tables 3.1 and 3.2).

---

## 5. Finite element library (Manual §1.4 and §1.5.3) **[CORE]**

### 5.1 Element types

| Type keyword (manual) | Dim. | Nodes | DOF per node | Formulation | Allowed use | Nonlinear capability | Notes |
|---|---|---|---|---|---|---|---|
| **SOLID** | 3D | 8 (degenerate pyramids allowed) | 3 translations (ux, uy, uz) | Isoparametric hexahedron. Optionally **9 incompatible displacement modes** when the element models the **structure**. | Structure, near-field soil, **excavated soil** | Equivalent-linear soil nonlinearity via HOUSE + STRESS iterations | No rotational DOF, so rotations are not transmitted from BEAM or SHELL elements |
| **BEAM** (listed as "BEAMS" in §1.4 and in result file names) | 3D | 3 | 6 (3 translations + 3 rotations) | 3D frame element | Structure, buried elements | Linear only | Consistent mass. **(Impl. note)** Two end nodes I, J and a third orientation node K defining the local axes. Confirm in the element-command spec. |
| **SHELL** | 3D | 4 (triangles allowed) | 6 | Thin plate/shell, **Kirchhoff** theory | Structure, buried (ETYPE = 2) | **Option NON** wall panels (in-plane shear or in-plane bending) | **Zero in-plane rotational (drilling) stiffness.** Use FIXROT / FIXSHLROT. Lumped mass. |
| **TSHELL** | 3D | 4 (triangles allowed) | 6 | Thick plate/shell, **Mindlin-Reissner** theory | Structure, buried | Option NON panels (as SHELL) | HOUSE adds a small drilling stiffness automatically, so FIXROT is not needed. Lumped mass. |
| **PLANE** | 2D | 4 (triangles allowed) | 2 translations | Plane strain | 2D structure, near-field soil, **excavated soil** in 2D | Equivalent-linear soil nonlinearity via HOUSE + STRESS | |
| **SPRING** | 3D | 2 | 6 (3 translational + 3 rotational stiffnesses) | Discrete spring | Structure, isolators, pile-soil interface | **Option NON** translational only, as **1D spring groups** | For nonlinear use, a 3D spring with up to 3 translational stiffnesses must be **split into groups of 1D nonlinear springs**, so that each DOF gets its own stiffness and damping. |
| **General Matrix (GM)** | 3D | 2 (I, J: input in **global** coordinates) or 3 (I, J, K: input in **local** coordinates) | 6 | User complex stiffness and mass matrices (super-element) | Super-elements, equipment, fluid effects, external subsystems | — | Input by the **MXI** command, `MXI,<p>,<row>,<t1>,<t2>,…<t12>`, one 12-term row per call, which implies a 12×12 matrix for two 6-DOF nodes (Ch. 9.4.26). Option AA converts ANSYS **MATRIX50** to GM automatically. |

**External loads** that can be applied: **nodal forces, nodal moments, nodal translational masses, nodal rotational masses.**

**Excavated soil volume** element types: **SOLID** (3D) and **PLANE** (2D) only.

### 5.2 DOF management **[CORE]**

- Each structural node can have up to 6 DOFs. The user can **delete (fix) DOFs** to reduce problem size. This is the D-command family (Ch. 9).
- **Interaction nodes carry only translational DOFs**, 3 in 3D or 2 in 2D. Structural rotations are passed to the soil through the nodal translations (§6, rule 5).
- **(Impl. note)** Build the active-DOF map from the union of DOFs the attached elements activate, minus DOFs fixed by the user. Warn about DOFs that are active but have zero stiffness: free rotations at SOLID-only nodes, and drilling at SHELL-only nodes. That is the purpose of FIXROT, FIXSLDROT, FIXSHLROT and FIXSPRROT.

### 5.3 Mass matrix convention (§1.5.1 item 12) **[CORE]**

- Default element mass matrix: **M = 0.5·M_lumped + 0.5·M_consistent**, for SOLID and PLANE (structure and soil).
- **BEAM**: consistent mass only.
- **SHELL / TSHELL (plate elements)**: lumped mass only.
- SPRING: massless. GM: user mass matrix. Nodal masses: added to the diagonal (translational and rotational).
- Ch. 4 cross-reference: this 50/50 mix is the reason the **1/5-wavelength** mesh-size rule is enough.

### 5.4 Damping convention (§1.5.1 item 13) **[CORE]**

- Material damping uses **complex moduli**. This gives effective damping ratios that are **independent of frequency** and can vary from element to element.
- **(Impl. note)** Implement the SASSI/SHAKE91 convention as the default: `E* = E·(1 − 2β² + 2iβ·sqrt(1 − β²))`. This is equivalent to the complex wave velocity `v* = v·(sqrt(1−β²) + iβ)`. Apply the same factor to G and to spring stiffnesses. Also offer the simpler `E·(1 + 2iβ)` as a selectable alternative for benchmarking. Confirm against the Chapter 2/3/6 specs. The complex element stiffness is `K* = Σ_e K_e(E*_e)`, assembled with each element's own β_e.

---

## 6. Modeling capabilities and limitations (Manual §1.5.1, items 1–13) **[CORE]**

Each item is phrased as a rule, followed by what the program must do.

**Rule 1. Soil layering.** The site is semi-infinite, horizontally layered viscoelastic soil, resting on either (a) a rigid base or (b) a semi-infinite elastic or viscoelastic half-space. SITE represents the half-space by **adding soil layers automatically** with a **viscous boundary** at the bottom. Soil nonlinearity is represented by the Seed-Idriss iterative equivalent-linear method.
- **WARNING (layer count).** Use **more than 20 soil layers**. This matters especially for uniform deposits, because too few layers spoil the Rayleigh and Love modes computed by SITE. User-defined **half-space layers: 10 to 20** computational layers.
  - Program action: emit a CHECK **warning** if the layer count is 20 or fewer, or if the half-space layer count is outside 10 to 20.
- **WARNING (equivalent properties).** The layer Vs and damping for seismic SSI should be the strain-compatible (equivalent-linear) values from SOIL. Depending on the SITE input selection, the program can pass the SOIL results (FILE88) to SITE automatically, and to HOUSE for the excavated soil layers.
  - Program action: offer a "use SOIL iterated properties" switch in SITE and HOUSE. **[IO]**
- SOIL also has a time-domain nonlinear option for soft soils: a hyperbolic hysteretic model like DEEPSOIL, without water-table effects. [ADV]

**Rule 2. Structure discretization.** Standard 2D and 3D FE connected at nodes, up to 6 DOFs per node, with user DOF deletion. GM elements carry super-element matrices. Option AA converts MATRIX50 to GM.

- **FIXROT / FIXSHLROT** (thin SHELL models). These improve numerical conditioning by removing singularities caused by the **zero drilling stiffness** of the Kirchhoff element.
  - Fixing a rotation by a D-constraint only works for shell faces parallel to a global plane. For **oblique** shells the commands add a **small torsional (drilling) spring** about the shell normal.
  - **Default spring stiffness = 10.** The user may change it, but it should be **no more than 10% of the shell bending stiffness**.
  - Syntax: `FIXROT,[Stiff]` and `FIXSHLROT,[stiff]`.
- **FIXSLDROT** fixes the unused free rotational DOFs at nodes connected only to SOLID elements.
- **FIXSPRROT** fixes the unused rotational DOFs of springs. Chapter 1 spells it "FIXSPROT"; the command reference uses **FIXSPRROT**. Accept both.
  - These save storage and run time.
- **FIXROT ≈ FIXSHLROT + FIXSLDROT + FIXSPRROT.** The combination of the three individual commands is said to approximate ANSYS-generated stiffness more closely.
- **WARNING.** For thin-shell models, especially oblique shells, use of FIXROT or FIXSHLROT is **highly recommended**. TSHELL gets an automatic small drilling stiffness in HOUSE, so FIXROT is not needed for TSHELL.
- **(Impl. note)** Define the drilling spring in the shell's local frame: `k_θn·(θ_n)²/2`, added at each node about the element normal and transformed to global axes. Its units are moment per radian in the model's units. "10" is an absolute value, so it is unit-dependent; see Open Questions.

**Rule 3. Mixed-DOF connections (hinges).**
- Where BEAM or SHELL/TSHELL nodes are shared with SOLID nodes, the **rotations are not transmitted** to the SOLIDs, which have 3 DOFs. To pass moments, add **massless** BEAM or SHELL/TSHELL elements along the SOLID edges or faces. These carry the moment into the solid as force couples at the SOLID nodes.
- BEAM-to-SHELL connections: the shell has no drilling stiffness, so the beam's in-plane rotation is not passed to the shell. Build a local **"tripod" of BEAMs** to spread the rotation to neighbouring nodes.
- Model checks:
  - **FIXEDINT** finds interaction nodes that are mistakenly fixed (any fixed translational DOF).
  - **HINGED** finds potential unintended hinges between SOLID, SHELL/TSHELL and BEAM elements, that is, a 6-DOF element meeting a 3-DOF element without penetrating along an edge or face.
- **WARNING.** Running FIXEDINT and HINGED before production runs is **strongly recommended**.
- **(Impl. note)** FIXEDINT: for each interaction node, report it if any of ux, uy, uz is constrained. HINGED: for each node shared by a 6-DOF element (BEAM, SHELL, TSHELL) and a SOLID, report it if the 6-DOF element connects to the solid mesh through only that single node (beam) or only one edge-collinear set of nodes, so the rotation is free. The exact algorithm is in the Chapter 9.8 spec.

**Rule 4. Excavated soil modeling.**
- The excavation is modeled with PLANE elements (2D) or SOLID elements (3D).
- **Without near-field soil:** excavated soil elements connect to the structure **only at the foundation-soil interface nodes**.
- **With near-field soil** (backfill, irregular layering): excavated soil elements connect **only to near-field soil elements**, never to structural elements. The connection is at the nodes on the boundary between near-field soil and far-field layering.
- **WARNING.** **EXCSTRCHK** is the key model check. It finds excavation-volume **internal** nodes that are wrongly shared with structure basement nodes.
  - Syntax: `EXCSTRCHK`, no arguments. It does not modify the model.
  - Output goes to the command history, and the list length follows the Check Options break number.
- **(Impl. note)** "Internal" excavation nodes are nodes of excavation elements (ETYPE = 2) that are **not** on the outer boundary surface of the excavation volume. A boundary face is an excavation element face not shared with another excavation element. Report any internal node that is also referenced by a structural element (ETYPE = 1) or a near-field soil element.

**Rule 5. Interaction nodes.** All interaction nodes must lie **on soil layer interfaces** and have **only translational DOFs**.
- Program action: CHECK error if an interaction node's elevation does not match a layer interface depth within tolerance. This means the excavation mesh's horizontal node levels must coincide with the SITE layer boundaries, so layering and mesh must be built together. Ignore rotations at interaction nodes in the soil coupling. **[CORE]**

**Rule 6. Substructuring methods.** These are defined by the choice of interaction nodes. **[CORE][ADV]**

| Method | Also called | Interaction nodes | Typical use / accuracy |
|---|---|---|---|
| **FV**, Flexible Volume | "Direct" method; the reference method | **All** nodes of the excavated soil volume | Reference accuracy. Costliest. |
| **FI-FSIN**, Flexible Interface, Foundation-Soil-Interface Nodes | **Subtraction Method (SM)** | Only the foundation-soil interface nodes | Fastest. Can become **numerically unstable at higher frequencies**, depending on soil stiffness and excavation geometry. Accurate for stiff soil and rock sites. |
| **FI-EVBN**, Flexible Interface, Excavated-Volume-Boundary-Nodes | **Modified Subtraction Method (MSM)**, sometimes **Extended Subtraction Method (ESM)** | All nodes on the excavation volume boundary | Several times faster than FV and only a few times slower than FI-FSIN. Good for **shallowly embedded** nuclear islands. **May fail for deeply embedded** models such as SMRs. |
| **FFV**, Fast Flexible Volume | — | The FI-EVBN nodes **plus a reduced number of internal horizontal node layers** in the excavation | For deeply embedded structures and SMRs. Highly accurate and up to tens of times faster than FV. Useful when FV would need 30,000 to 50,000 interaction nodes. |
| Surface (INTGEN type 4) | — | Nodes at ground-surface level only | Surface foundations |

- **INTGEN** generates interaction nodes automatically: `INTGEN,<type>,[level skip]`. From Ch. 9.7.19:
  - type 0 = set all nodes to non-interaction;
  - type 1 = FV;
  - type 2 = FI-EVBN (MSM);
  - type 3 = FI-FSIN (SM);
  - type 4 = Surface foundation;
  - type 5 = FFV with multiple internal levels. `[level skip]` (default 1) is the number of node levels skipped between internal interaction levels.
  - Types 1 to 5 do not remove existing interaction nodes.
  - The excavation must be explicitly flagged with ETYPE or ETYPEGEN for types 1 to 3 and 5.
- Theory of FV, SM and MSM is in Manual §2.1 (Chapter 2 spec).
- **WARNING (ASCE 4-16).** Using any method other than FV (SM, MSM or FFV) requires a **preliminary validation study against FV**.
  - ASCE 4-16 allows an "excavated-soil-only" ("swimming pool") model, compared by ATF at the nodes shared between structure and excavated soil.
  - The manual instead suggests a **simplified massless foundation** model for kinematic SSI, or best, the **full SSI model** compared by ATF at critical in-structure locations.
  - For deep soft deposits and DES, the excavated-soil-only model can be poorly conditioned and should be avoided.
- **WARNING (quarter models).** ASCE 4-16 suggests quarter models for the SM/MSM vs FV comparison. The manual **does not recommend** quarter models for qualifying SM/MSM, because the symmetric and antisymmetric boundary conditions can hide instabilities that the full model shows.

**Rule 7. Method trade-off.** Soil impedance cost grows with the **power 2 to 3 of the number of interaction nodes**. FV can take tens of times longer than FI. FI-EVBN is the usual accuracy/speed compromise. Sensitivity studies comparing FI with FV are always recommended for embedded structures.
- **(Impl. note)** Have the UI print a cost estimate, about N_int³ per frequency, before the run.

**Rule 8. Excavation boundary completeness.** Interaction nodes must **always include every node on the excavation / far-field soil interface**. Skipping any of them can badly hurt accuracy, so do sensitivity studies if any are skipped.
- Program action: CHECK warning if a node on the excavation boundary is not an interaction node, for methods FV, MSM and FFV.

**Rule 9. Mesh uniformity and POINT radius.**
- Keep the excavation mesh **as uniform as possible**. Non-uniform meshes give non-uniform local impedances. Use transition meshes between the excavation mesh and the basement (BNL / Nie et al. 2013). An automatic transition-mesh feature is planned but not available.
- **Axisymmetric-model radius** (used by POINT to compute the soil layer flexibility):
  - For non-uniform meshes there is no exact rule. Practical rule: use the **average radius over all excavation-volume solid elements**.
  - Run sensitivity cases with the **maximum, average and minimum** radius. If the results differ significantly, **envelope** them.
  - Non-uniform meshes can amplify torsion and rocking, especially under incoherent input.
- **(Impl. note)** Define the element "radius" as the equivalent radius of the element's plan (horizontal) area, `r_e = sqrt(A_plan,e / π)`. Provide a UI helper that reports the min, average and max r_e over the excavation solids. See Open Questions.

**Rule 10. Interaction nodes outside the excavation.** Buried BEAM or SHELL/TSHELL elements (which have no volume) may reach outside the basement envelope. Their nodes must also be interaction nodes.
- For shells this is automatic when **ETYPE = 2** is combined with **INTGEN**. ETYPE = 2 has no effect on beams.
- ETYPE values (Ch. 9.4.11):
  - 0 = default: SOLIDs below grade are excavation, otherwise structural; SHELLs are structural;
  - 1 = structural;
  - 2 = excavated soil (SOLID/PLANE) or embedded shells (SHELL/TSHELL).

**Rule 11. Separate meshes for structure and excavation.**
- Excavated soil nodes must be **different from** structural basement nodes, **except** at the foundation-soil interface. There they share nodes, or are tied with rigid or very stiff springs.
- With backfill, connection happens only on the outer surface of the excavated soil (the near-field / far-field interface).
- All other internal excavated-soil nodes must be **independent** of the basement nodes, **even when they have the same coordinates**.
- Breaking this rule can give very poor results, especially with FI methods. FV can tolerate it for stiff basements, but **flexible basement subsystems** such as piping and equipment are affected significantly.
- The EXCSTRCHK warning applies again here.
- **WARNING (isolated-frequency instability).** Embedded models can show numerical instability at isolated frequencies with any method, FV included.
  - Inspect the ATFs at several nodes.
  - Exclude spurious frequencies from the TF interpolation in MOTION and STRESS.
  - Run sensitivity checks at adjacent frequencies.
  - The tool **Remove_Frequencies_from_FILE8.exe** (installed zipped under `C:\ACSV300\`) deletes chosen frequencies from FILE8.
  - **(Impl. note)** Implement this as a command, for example `FILE8 remove-frequencies <list>`, that writes a new FILE8 without the chosen frequency records and keeps a backup. **[IO][UI]**

**Rule 12. Mass matrix** — see §5.3.

**Rule 13. Damping** — see §5.4.

---

## 7. Dynamic loading for seismic and vibration analysis (Manual §1.5.2) **[CORE]**

1. **Seismic environment.** Any 3D superposition of **inclined body waves and surface waves** (SITE wave-composition input).
2. **Control motion.**
   - An acceleration time history, the **control motion**, assigned to one of the **three global directions**.
   - It is applied at the **control point**, which **lies on a soil layer interface** (SITE input).
   - Program action: CHECK that the control-point depth is a layer interface.
3. **Coherent and incoherent input**, with or without directional wave passage.
   - Use different incoherency models for rock and soil sites.
   - Up to 50 stochastic simulations, each with X, Y and Z components, can run in one SSI run.
4. **External forces and moments** (impact, wave forces, rotating machinery) are applied at nodes.
   - All loads share **one time-history shape** f(t). Each nodal load can have its own **maximum amplitude a_j** and **arrival time t_j**, which lets the user model **moving loads**.
   - **(Impl. note)** In the frequency domain, `F_j(ω) = a_j · F(ω) · exp(−iωt_j)`, where F(ω) is the FFT of the common shape.
   - Up to 500 load cases per run.
5. **FFT-based transient analysis.**
   - Inputs must be sampled at a **uniform Δt**.
   - The **Fourier period** `T = N·Δt` must be **much longer than the excitation**.
   - The points after the excitation must be **zero**. This is the "quiet zone" (also called "trailing zero part"), during which the SSI system vibrates freely. Its length depends on the **lowest damping** in the system.
   - **N must be a power of 2.** Maximum **32,768** points here (see the §3.2 note on 65,536).
   - Frequency step: `Δf = 1/(N·Δt)` (Ch. 4 cross-reference).
   - **(Impl. note)** CHECK warning if the quiet zone is shorter than `T_q ≥ ln(1/ε)/(2π·f_min·β_min)`, with ε = 0.01 by default, where f_min is the lowest significant system frequency. Pad with zeros to the next power of 2 if the user asks; Chapter 10 Warning 9 says PREP writes the nearest power of 2.
6. **Seismic and external-force analyses cannot be in the same run.** Run them separately and combine the final results by linear superposition.

---

## 8. Nonlinear soil and structure hysteretic behaviour (Manual §1.5.4)

### 8.1 Equivalent-linear iteration framework **[CORE][ADV]**

- The substructuring solution in complex frequency is **linear**. Nonlinear soil and structure behaviour is handled by **iterative equivalent linearization**. Each iteration is a fast **restart** of ANALYS using **"New Structure"**.
- **Local soil nonlinearity** is modeled with 3D SOLID elements (3D models, typically next to the foundation) or 2D PLANE elements (2D models).
  - Only HOUSE and STRESS run differently. ANALYS runs as a "New Structure" restart. The loop is **HOUSE → ANALYS(restart) → STRESS**.
- **Effective strain measure** (user choice):
  - (i) the effective soil **shear strain component**;
  - (ii) the **maximum shear strain** in 2D soil elements;
  - (iii) the **octahedral shear strain** in 3D soil elements.
  - **(Impl. note)** Standard definitions:
    - `γ_max,2D = sqrt((ε_x − ε_y)² + γ_xy²)`
    - `γ_oct = (2/3)·sqrt((ε_x−ε_y)² + (ε_y−ε_z)² + (ε_z−ε_x)² + 1.5(γ_xy² + γ_yz² + γ_zx²))`
    - Effective strain = `R_γ · max_t |γ(t)|`, with the strain ratio R_γ set by the STRESS / SOIL options (SHAKE commonly uses 0.65). Confirm in the Chapter 6 STRESS spec.
- **Directional combination.**
  - STRESS writes the effective shear strain for each input direction (X, Y, Z) to **FILE74**.
  - **COMB_XYZ_STRAIN** combines them at each iteration, by **SRSS** of the directional effective strains (Ch. 4 cross-reference).
  - HOUSE reads the combined FILE74 to update G and D for the next iteration, and writes **FILE78**, which STRESS uses.
- **Soil curves.** The user enters G/Gmax(γ) and D(γ) curves in SOIL. SOIL writes them to **FILE73**, which STRESS reads during the SSI iterations.
- **1D (horizontal layer) soil.** SOIL alone handles the hysteresis, either with the Seed-Idriss equivalent-linear method (SHAKE) or with the DEEPSOIL-like true nonlinear model.

**(Impl. note) Iteration algorithm, nonlinear soil.**
```
props_0 = low-strain (or SOIL-iterated) G, D for each nonlinear soil element group
for k = 1..k_max:
    HOUSE(props_{k-1}) -> FILE4/FILE78     # only the structure + near-field matrices change
    ANALYS restart "New Structure"          # reuse soil impedance per frequency (COOX…)
    STRESS -> γ_eff per element per direction -> FILE74(X), FILE74(Y), FILE74(Z)
    COMB_XYZ_STRAIN: γ_eff = sqrt(γx² + γy² + γz²)
    props_k = curves(γ_eff)  (G = Gmax·G/Gmax(γ_eff), D = D(γ_eff))
    stop if max_e |G_k − G_{k−1}|/G_k < tol_G and |D_k − D_{k−1}| < tol_D
```
The convergence tolerances are not given in Chapter 1; see Open Questions.

### 8.2 Option NON: nonlinear structure **[ADV]**

**Element choices.**
- (i) Nonlinear **RC wall panels** ("macroshell" elements), which are groups of SHELL elements.
- (ii) **Nonlinear springs** with translational DOFs.

**Applications.**
- RC cracking and post-cracking in **low-rise shearwalls**.
- Rubber base-isolators.
- Pile-soil interface nonlinearity.
- Limited foundation sliding.
- Better inter-building displacement predictions for gap sizing. Incoherency can significantly reduce the predicted gaps.

**Wall panel definition.**
- A panel is a subdivision of a wall assumed to deform in **uniform shear or uniform bending**.
- It is a set of shell elements in a **vertical plane**.
- Its shape should be **rectangular**, though triangular shells are allowed inside.
- It can have **any orientation in plan**, not necessarily parallel to the global axes.
- Defined with the `P,<num>,<group>,<bbc>,<disp>,<force>` command (Ch. 9.17.23). All shells in the panel group must be coplanar.
- `NONLINMOTDISP` finds the panel corner nodes and adds them to the MOTION and RELDISP output requests.

**Computational steps (NONLINEAR module).**
1. Initial **linear elastic (uncracked)** SSI analysis with the initial elastic properties of the nonlinear elements.
2. Compute the local nonlinear-element behaviour **in the time domain** from the local relative displacements. Use it to calibrate a linearized hysteretic model per element, in complex frequency.
3. Run a new SSI iteration with a fast **restart** in complex frequency, using the linearized properties from step 2.
4. Check convergence and either stop or go back to step 2.

**WARNING (3-directional combination).** After each iteration, combine the co-directional effects so that nonlinear behaviour is driven by the **three-component** input, not by each direction separately.
- **COMB_XYZ_STRAIN** does this for nonlinear SOLID/PLANE soil, from the STRESS strains.
- **COMB_XYZ_THD** does this for Option NON shell walls and springs, from the RELDISP relative displacements.
- COMB_XYZ_THD is put into the batch file automatically by **`NONLINBAT,1`** (`NONLINBAT,<Sel>`: 0 = single direction, 1 = three directions).
- Its input file **`COMB_XYZ_THD.inp`** must be written by the user. An example is on the installation DVD; the format is not given in Chapter 1, see Open Questions.

**Recommended practice (axial).**
- There is **no hysteretic model for axial deformation** in panels.
- ASCE 4-17 reduces only shear and bending stiffness for cracking and keeps axial stiffness unchanged.
- So treat the structure as **nonlinear under the two horizontal components and linear under the vertical component**. The panel corner displacements used each iteration must still include the combined effect of all three components.

**NONLINEAR outputs.**
- Nonlinear structural displacements, accelerations and element forces.
- **Ductility ratios** relative to the initial stiffness (that is, relative to cracking for concrete).
- **Inelastic absorption reduction factors**, defined as initial elastic force ÷ effective nonlinear force. These depend on the strain level and correspond to the ASCE 43-05 inelastic energy absorption factors **F_μ** for the different limit states. The likely file is `*.FMU`.

**Figure 1.2 (transcribed): nonlinear SSI flow and files.** Use the string "PANEL" for panels and "SPRING" for springs.

```
Initial Elastic SSI Analysis                         Iterative SSI Analyses (loop)
Step 1: SITE, POINT, HOUSE, ANALYS, MOTION, RELDISP   -> HOUSE, ANALYS (Restart), MOTION, RELDISP
Step 2: COMBIN_XYZ_THD, NONLINEAR                      -> COMBIN_XYZ_THD, NONLINEAR -> CONVERGENCE?
Step 3: copy & save output files (elastic run):          Result files for iteration #:
   Panelxxxx.THD          -> Panelxxxx_elastic.THD          Panelxxxx_it#.THD and Panelxxxx_it#.THS
   Panelxxxx.THS          -> Panelxxxx_elastic.THD (sic;    Panelxxxx_it#.crv
                              presumably _elastic.THS)      Panel_it#.FMU
   Panel.FMU              -> Panelxxxx_elastic.FMU          Panel_EQL_Matl_Prop_It#.txt
   Modelname.hou          -> Modelname_elastic.hou          Modelname_it#.hou
   Modelname_new.hou      -> Modelname.hou   (updated model becomes the next HOUSE input)
   FILE8                  -> FILE8_It1
   PANEL_EQL_MATL-PROP.TXT-> PANEL_EQL_MATL_PROP_IT1.TXT
```

(Impl. note) Probable file meanings:

| File | Probable content |
|---|---|
| `.THD` | Panel deformation (strain or curvature) histories |
| `.THS` | Panel force histories |
| `.crv` | Hysteresis loop (force–deformation) data |
| `.FMU` | F_μ inelastic absorption factors |
| `*_EQL_MATL_PROP*` | Equivalent-linear E and damping for each panel |
| `Modelname_new.hou` | HOUSE input regenerated by NONLINEAR with updated panel and spring properties |

**Hysteretic models (low-rise RC shearwalls).**

| No. | Model | Status | Use |
|---|---|---|---|
| 1 | **Cheng-Mertz Shear (CMS)** | implemented | Panels dominated by in-plane **shear**. Typical low-rise walls, height < width. Most appropriate for nuclear buildings. |
| 2 | **Cheng-Mertz Bending (CMB)** | **not included in this version** | Bending |
| 3 | **Takeda (TAK)** | implemented | Panels dominated by in-plane **bending** (manual Ch. 6 says both bending and shear) |

**Figure 1.3** compares CMS, CMB and TAK force–strain loops for the same strain history:
- The BBC is a smooth curve rising to a plateau, peaking at about 5000 force units near 0.0002 to 0.0004 strain, then flat out to about 1.2e-3.
- Panel titles give two cases. **Dc = 0.000062, Vc = 2470** ("higher cracking point", no pinching in CMS). **Dc = 0.00002, Vc = 795** ("low cracking point", which gives **pinched** CMS loops).
- **Dc / Vc are the cracking deformation and cracking force** that define the first BBC segment.
- The detailed CMS and TAK loop rules are not in Chapter 1; see Open Questions and the Chapter 6/9.17 specs.

**Panel kinematics (Figure 1.4).** **[ADV][CORE]**
- Shear strain and bending curvature are computed from the displacements of the **four corner nodes**.
- They are expressed in the **local wall-plane axes**, after **removing the panel rigid-body motion**.
- **Shear strain** is the **average rotation of the two vertical edges**: γ = (γ₁ + γ₂)/2, where γᵢ is the angle between the deformed and undeformed vertical edge i.
- **Bending rotational strain (curvature)** comes from the **relative rotation of the horizontal edges**, θ₁ and θ₂ in the figure, measured between the deformed and undeformed top edge.
- The **uniform vertical axial strain** is also computed, so the user can get the seismic axial forces. Axial effects are usually small for low-rise walls.

(Impl. note) Recommended formulas:
- Local axes: `e_h` is the unit horizontal vector along the panel bottom edge, `e_v` is global Z, `e_n = e_h × e_v`.
- Corners: B1 bottom-left, B2 bottom-right, T2 top-right, T1 top-left. L is the width and H the height.
- In-plane displacements: `u = d·e_h` and `w = d·e_v`.
- Formulas:
```
γ   = ½[(u_T1−u_B1)+(u_T2−u_B2)]/H + ½[(w_B2−w_B1)+(w_T2−w_T1)]/L     (rigid-rotation invariant)
θ_t = (w_T2 − w_T1)/L ;  θ_b = (w_B2 − w_B1)/L ;  κ = (θ_t − θ_b)/H       (bending curvature)
ε_v = ½[(w_T1−w_B1)+(w_T2−w_B2)]/H                                       (axial strain)
```
These are evaluated at each time step, from the 3-direction combined corner displacement histories.

**Property update.**
- After each iteration, the equivalent-linear properties for each panel are applied as a **reduction of the elastic modulus E** of the panel's shells.
- **Poisson's ratio is kept constant**, so shear, axial and bending stiffness degrade **equally** (shear and bending are fully coupled).
- This is acceptable for low-rise, shear-dominated walls. Gergely (NUREG/CR-4123) and Reed & Kennedy (EPRI 1994) note that flexural effects play a negligible role there.

**Axial and gravity effects.**
- If axial force matters, the panel shear and bending capacities, and therefore the BBC, must include the combined seismic and gravity vertical force.
- Get the gravity forces by converting to ANSYS and running a static analysis there, or by a static analysis inside ACS SASSI.
- **Static gravity trick in ACS SASSI** (frequency-domain code). Apply a **vertical acceleration of amplitude 1g** as a **one-period harmonic** whose period is **exactly the Fourier period** T = N·Δt. Example: N = 32,768 and Δt = 0.005 s give T = 163.84 s.
- For combined analysis, superpose the vertical seismic motion on that slow harmonic, ideally around the middle of the **first half-period**, with the half-period defined with a **negative sign** so it acts in the gravity direction.
- **(Impl. note)** Provide a helper that generates this record: `a(t) = −g·sin(2πt/T)`, with the seismic record inserted around t ≈ T/4.

**Backbone curves (BBC).**
- The user defines the material backbone curves and the hysteretic model type. **Defining the BBC is the analyst's responsibility.**
- BBC quantities:

| Element | BBC relationship | Initial slope |
|---|---|---|
| Spring | Force vs relative displacement | Spring stiffness k |
| Shear panel | In-plane shear force vs shear strain | G·A_shear |
| Bending panel | In-plane bending moment vs curvature | E·I |

- BBC input:
  - from an external data file, with the **BBC** command;
  - directly as point vectors, with **`BBCX,<num>,<points>,<yield>,<X1>…<Xn>`** for strain and **`BBCY,<num>,<points>,<yield>,<Y1>…<Yn>`** for force or moment. Up to 100 points (NONLINEAR limit).
- Ways to build the BBC:
  1. From published capacity equations for the peak or ultimate in-plane capacity, together with the elastic stiffness.
  2. From a separate static nonlinear pushover analysis, done outside ACS SASSI.
- The **SHEAR** command (`SHEAR,<panel>,[fc],[fy],[P],[Nu],[Fvw],[Fbe]`) computes and compares shear capacities from different equations: Gulec & Whittaker 2009, Wood 1990, ACI 349-08, Barda et al. 1977.
- **BBCGEN** (`BBCGEN,<Panel>,<ShearModel>,[fc],[fy],[Pn],[Nu],[bre],[bys],[CrackingForceLevel]`) generates smooth BBCs for many panels at once.
- **Avoid the Barda equation** for typical nuclear shearwalls. It applies only to barbell walls with heavy flanges, can overestimate capacity, and has been removed from the ASCE 43-17 draft. ATC 72-1 Option 3 (2010) gives guidance for reducing the peak capacity for cyclic degradation.
- If it is unclear whether a panel is governed by shear or bending, run preliminary comparison analyses with each assumption.

**Walls with large openings (Figure 1.5).** Use **EDGE** to split panels into subpanels around openings: `EDGE,<panel>,[X],[Y],[Z]`.

Worked example:
- Start with a single panel in group 1.
- `EDGE 1,0,0,1` splits it horizontally into 3 subpanels. Group 1 is the top band, group 2 the band containing the openings, and group 3 the bottom band.
- `EDGE, 2` then splits group 2 into groups 2, 4 and 5, which are the piers between and beside the openings.
- The result is 5 subpanels: 1 (top), 2 (left pier), 4 (middle pier), 5 (right pier), 3 (bottom).
- Each subpanel needs its own BBC, from its geometry and reinforcement.

Figure 1.5 (lower right) shows a shear stress/force vs shear strain BBC:
- solid wall: higher peak, "global wall failure";
- wall with openings: lower peak with faster softening, "local wall failure" in the piers next to the openings.

**Damping options for nonlinear elements.**
- Total damping is either **hysteretic only**, or **viscous-elastic + hysteretic**. The UI option is "Include Elastic Damping" (Ch. 6). If hysteretic only, it should not be smaller than the initial viscous damping.
- A **damping cut-off value** keeps the effective damping at or below the code maximum. For example, a **7% cut-off** for cracked-concrete design-basis analysis, per ASCE 4/43 and USNRC.
- A **damping scale (reduction) factor** can calibrate the damping to test data.

**WARNINGS.**
- **BBC shape.** BBCs must be smooth. **Bilinear or trilinear BBCs with sharp corners can badly hurt convergence**, giving oscillation or divergence.
- **BBC initial slope.** The first-segment slope of each BBC **must equal the elastic stiffness in the initial HOUSE input (`.hou`)**:
  - for shells, the elastic modulus E;
  - for springs, the spring stiffness.
  - NONLINEAR **normalizes each BBC by its first-segment slope** and applies the normalized curve to the initial E (or k) from the `.hou` file.
  - For concrete panels the cracking-point slope must equal **G·A_shear** for shear BBCs and **E·I** for bending BBCs.
  - Include vertical axial force effects in the BBC if they are significant.
  - Program action: CHECK warning if `|BBC initial slope − (G·A_shear or E·I or k from .hou)|/ref > tolerance`, for example 1%.
- **Nonlinear springs.** Split them into 1D spring groups, one translational DOF per group.
- **Frequencies.** Nonlinear SSI needs **more SSI frequencies** than linear SSI. Make the frequency grid denser where spectral peaks are expected to move down as the structure softens.

**(Impl. note) Equivalent-linearization of a nonlinear element**, not specified in Chapter 1:
```
input: combined deformation history δ(t) (γ, κ or Δu), BBC F_bb(δ), hysteretic model rule set (CMS|TAK|spring)
1. run hysteretic model over δ(t) → F(t); record loops
2. δ_eff = max|δ| (or R·max|δ|; option), F_eff from loop at δ_eff
3. secant stiffness ratio  r_k = (F_eff/δ_eff) / K0,  K0 = BBC first-segment slope
4. hysteretic damping  ξ_h = E_D / (4π·E_S),  E_D = loop area, E_S = ½·F_eff·δ_eff (averaged over cycles)
5. ξ = (ξ_h [+ ξ_0 if Include Elastic Damping]) × scale_factor;  ξ = min(ξ, cutoff)
6. E_new = r_k · E_0  (shell panels) ; k_new = r_k · k_0 (springs);  write Modelname_new.hou
7. converge when max |E_k − E_{k−1}|/E_k and |ξ_k − ξ_{k−1}| are below the Displacement/Force convergence errors
```
Ductility = δ_max/δ_cr. F_μ = F_elastic(δ_max)/F_nonlinear(δ_max) = (K0·δ_max)/F(δ_max).

---

## 9. SSI solution interpolation in frequency (Manual §1.1(vi), §1.5.4 item 7) **[CORE]**

The complex TFs are computed at NF selected SSI frequencies, integer multiples of Δf. They are interpolated onto every Fourier frequency up to the cut-off, for both nodal acceleration TFs (**ATF**) and stress/force TFs (**STF**). Seven options exist. The names come from the Ch. 6 MOTION/STRESS dialogs and are used identically in both modules:

| Option | Name / description | Windowing | Averaging |
|---|---|---|---|
| **0** | SASSI2000 scheme: dense overlapping windows | moving window, overlapping | **weighted** averaging |
| **1** | Original SASSI 1982 scheme | non-overlapping windows | none |
| **2** | Dense overlapping windows | overlapping | averaging |
| **3** | Only three overlapping windows | overlapping (3) | averaging |
| **4** | Non-overlapping windows with **one** position shift | non-overlapping | none |
| **5** | Non-overlapping windows with **two** position shift | non-overlapping | none |
| **6** | **Complex bicubic spline** interpolation | none; needs a denser frequency grid | n/a |

Notes and warnings:
- Options 0 to 5 can produce **spurious narrow-band peaks and valleys** when there are too few frequencies. This is common for complex nuclear models with several structures on a common mat and dense mid/high-frequency modes. The fix is to add frequencies.
- Option 6 is the **best choice for incoherent SSI**, with **at least about 200 frequencies**, because it does not overshoot. It should be used **only when the frequency set is dense enough that spectral peaks are not clipped or smoothed**.
- Recommended number of SSI frequencies:
  - **40 to 80** for simple stick models;
  - **100 to 200** for complex flexible FE models with many local modes (coherent);
  - **200 to 300** for incoherent analysis;
  - more for nonlinear SSI.
  - Details are in Manual §4.1.2 (Chapter 4 spec).
- Accuracy check: compare the computed TF (`.TFU`) with the interpolated TF (`.TFI`) in "Spectrum TFU-TFI" plots. **[UI]**
- The interaction between the smoothing parameter and option 6 (smoothing must be 0 with the spline) is covered in the Chapter 6 MOTION spec.

**(Impl. note) Recommended algorithms** (the manual does not give formulas):
- **Options 0 to 5.** Classic SASSI rational, SDOF-like interpolant through **3 consecutive computed frequencies** (a "window"):
  - Form: `H(ω) = (A + B·ω²)/(1 + C·ω²)`, with complex A, B, C solved from the 3 points as a 3×3 complex linear system. This form reproduces exactly the TF of an SDOF oscillator with hysteretic damping.
  - **Non-overlapping** (option 1): windows {f₁,f₂,f₃}, {f₃,f₄,f₅}, … Each window evaluates the Fourier frequencies inside its span.
  - **Options 4 and 5**: the same, with the window start shifted by 1 or 2 computed points. The leading intervals use the first available window.
  - **Overlapping** (options 0, 2, 3): every window {f_k, f_{k+1}, f_{k+2}} that covers the interval contributes, and the results are averaged.
    - Option 2: simple average over all covering windows.
    - Option 3: at most 3 windows (the centered one and its two neighbours).
    - Option 0: weighted average, with the weight decreasing with the distance from the target frequency to the window center.
  - Below f₁: use the window that contains f = 0. For seismic ATF the value at 0 is 1, the rigid-body limit; for force TFs it is the static value. Above f_NF: TF = 0 (cut-off).
- **Option 6.** Cubic spline fitted **separately to Re(H) and Im(H)** against frequency (not-a-knot end conditions), evaluated at the FFT frequencies up to the cut-off.
- Validate with tests that (a) reproduce an analytic SDOF TF with 3 to 5 points per half-power bandwidth, and (b) compare against SASSI2000 benchmark outputs when available.

---

## 10. System of units (Manual §1.5.5) **[CORE][IO]**

- Any system of units may be used, **as long as it is consistent across all modules**.
- Gravity acceleration: **g = 32.2 ft/s²** for British Units (**BU**) or **g = 9.81 m/s²** for International Units (**IU**).
- **(Impl. note)** Store a model-level unit system flag, `BU` or `IU` (optionally `USER` with a user-set g). Use g wherever accelerations are given or reported in g units (EQUAKE, MOTION, ISRS, the gravity trick). Never convert other quantities internally. Echo the units in every module output header.

---

## 11. User data input (Manual §1.6) **[IO][UI]**

Input can be entered in three ways, mixed freely and in any order:
1. The **instruction (command) line**.
2. **Menu commands** and dialog boxes. Each dialog must produce the equivalent command, so the session can be replayed.
3. **Import of input text files** with **`INP,<fname>`** (Manual §9.2.18).
   - `INP` switches input to the file and returns to the keyboard at EOF.
   - The default path is the active model's directory.
   - The file is a `*.pre` file, written by the user or by the `WRITE` command.

**Command grammar** (PREP instructions, free format):
```
command   := keyword { "," param }
keyword   := alphanumeric token, case-insensitive (e.g. "afwrite" == "AFWRITE")
param     := number | string | empty            (empty → default value of that argument)
```
- The general form is `Keyword, p1, p2, ..., pn`, with **commas** as delimiters.
- Instruction names are **not case-sensitive**.
- Data may be entered **in any order**, alternating the command line, dialogs, menus and data files. **(Impl. note)** The model is therefore a declarative database, and module input files are generated only at AFWRITE time.
- **(Impl. note)** Points the manual does not define here:
  - Trim whitespace around tokens.
  - Allow a space after the keyword in place of the first comma. Figure 1.5's example `EDGE 1,0,0,1` uses a space.
  - Treat missing trailing optional parameters (shown as `[x]` in the Command Reference) as defaults.
  - Comment syntax and variable/macro syntax are defined in Ch. 5.6/5.9/9.15 (other specs).
- The `.pre` format does not depend on model size (§3.1).

---

## 12. Data checking (Manual §1.7) **[UI][IO][CORE]**

Data checking methods:
- **Data lists** generated while entering data (list commands such as GLIST, Ch. 9).
- **Interactive plots** of the whole model or selected parts, soil layers, time-history files, response spectra, and G/Gmax-strain and damping-strain curves.
- **Colour coding** of elements by group, property, material and so on.
- **Plot Info window**: **right-click** on an element to show all its properties.
- The **CHECK** command (Manual §9.2.7) verifies the data and shows errors and warnings in the **Check window**.
  - CHECK **simulates writing the analysis files**, so the analysis parameters must be set first.
  - It can run on a partial model.
  - Results are written to the model's **`.err`** file.
- **AFWRITE** (Manual §9.2.3) writes the requested module input files: `model.equ`, `.soi`, `.sit`, `.poi`, `.hou`, `.frc`, `.anl`, `.mot`, `.rdi`, `.str`, placed in the model directory.
  - It **always runs CHECK first**.
  - **If CHECK finds errors (warnings excluded), the affected analysis file(s) are not written.** Other modules' files are still written.
  - Model name and path must be defined, by the database or by `MDL,<Model>,<Path>`.
- **Message limits.** The **Options/Check** menu sets the maximum number of messages. Chapter 1 cites §6.6.1, but the dialog is actually §6.5.3; §6.6.1 is View → "Check Errors", which displays the `.err` file.
  - Options: show or hide errors and warnings; a per-type, per-module "break" number.
  - CHECK still counts **all** messages, but displays only up to the break number.
  - These options are global, **not saved**, and reset to defaults each time the UI starts.

**Check window format** (transcribed from the screenshot; window title "CHECK: Errors and Warning for - <modelname>"):
```
Errors and Warnings for EQUAKE
Error 85 : RS Input File 1 Does Not Exist
Error 85 : RS Input File 2 Does Not Exist
Error 85 : RS Input File 3 Does Not Exist
Errors and Warnings for SOIL
Error 73 : Acceleration Time History File Does Not Exist
Errors and Warnings for SITE
Errors and Warnings for POINT
Errors and Warnings for HOUSE
Errors and Warnings for ANALYS
Errors and Warnings for MOTION
Error 73 : Acceleration Time History File Does Not Exist
```
**(Impl. note)** Format each message line as `Error <n> : <text>` or `Warning <n> : <text>`, grouped under `Errors and Warnings for <MODULE>` headers in run order. Use the numbered error and warning catalogue in Manual Chapter 10 (spec 10) so numbers match the original, for example Error 73, Error 85 and Warning 9 "Number of Values for Fourier Transform Is Not Power of 2".

**Model-check commands named in Chapter 1** (full specs in Ch. 9.8): **EXCSTRCHK**, **FIXEDINT**, **HINGED**. Related commands in the same chapter: FREESPRING, INTCOUNT, USED.

**(Impl. note) CHECK rules derived from Chapter 1.** Add these to the CHECK catalogue as warnings unless noted:

| # | Rule | Severity |
|---|---|---|
| C1 | Soil layers ≤ 20 | warning |
| C2 | Half-space layers outside 10–20 | warning |
| C3 | Interaction node not on a layer interface | error |
| C4 | Interaction node with a fixed translation (FIXEDINT) | error |
| C5 | Internal excavation node shared with structure (EXCSTRCHK) | error |
| C6 | Excavation boundary node not an interaction node (for FV/MSM/FFV) | warning |
| C7 | SHELL model without FIXROT/FIXSHLROT and with oblique shells | warning |
| C8 | FFT length not a power of 2 | Warning 9, round to the nearest power of 2 |
| C9 | Quiet zone too short | warning |
| C10 | Seismic and force loads both requested in one run | error |
| C11 | Number of SSI frequencies below the recommended count for the selected interpolation option, incoherency or nonlinear mode (option 6 with fewer than 200) | warning |
| C12 | BBC first slope ≠ HOUSE elastic stiffness | warning |
| C13 | Sharp-cornered (bilinear/trilinear) BBC | warning |
| C14 | 3D nonlinear spring not split into 1D groups | error |
| C15 | SRSS incoherency selected for a flexible-foundation FE model | warning |
| C16 | Non-FV method selected (ASCE 4-16 validation reminder) | info |
| C17 | Size limits of §3.2 | warning (in the compatibility profile) |

---

## 13. Verification hooks derived from Chapter 1 **[CORE]**

These are suggested unit and regression tests that make the implementation verifiable against the statements in this chapter:
1. **Mass matrix.** For SOLID, PLANE, BEAM and SHELL, the total translational mass equals ρ·V (or the beam/shell equivalent). Check the 50/50 mix: the diagonal of a SOLID's M must equal 0.5·lumped + 0.5·consistent.
2. **Complex damping.** An SDOF with hysteretic β: peak |H| ≈ 1/(2β) for small β. With E·(1−2β²+2iβ√(1−β²)), the free-field layer amplification must match SHAKE91 for a uniform layer on rigid rock: `|H| = 1/|cos(ωH/v*)|`.
3. **Interpolation.** Options 0 to 6 recover an analytic SDOF TF with enough points. Deliberately sparse grids reproduce the "spurious peak" behaviour that the TFU/TFI comparison is meant to catch.
4. **FFT rules.** Power-of-2 enforcement. Δf = 1/(NΔt). Zero quiet-zone padding. Moving loads by phase shift exp(−iωt_j) match a time-shifted input.
5. **Gravity trick.** The one-period 1g harmonic (N = 32,768, Δt = 0.005 s, T = 163.84 s) on a fixed-base model returns quasi-static displacements equal to a static 1g analysis, within about 1% at T/4.
6. **Model checks.** Small models with planted errors (a shared interior excavation node, a fixed interaction node, a beam meeting a solid at one node, an oblique Kirchhoff shell without FIXROT) must trigger EXCSTRCHK, FIXEDINT, HINGED and the FIXROT warning respectively.
7. **Substructuring consistency.** For a surface rigid massless foundation on a layered site, FV, MSM, SM and FFV must give ATFs consistent with each other, and with published impedance solutions for the V&V categories listed in §2.
8. **Option NON.** A spring-only base-isolation model compared against an independent nonlinear time-history analysis (mirroring V&V Problem 51, which used ANSYS V15). Since the user uses ANSYS, an APDL export of the same model is the natural benchmark.

---

## 14. Cross-reference index (names introduced in Chapter 1)

| Name | Kind | Where specified |
|---|---|---|
| EQUAKE, SOIL, SITE, POINT2/POINT3, FORCE, HOUSE(FS/FSA), ANALYS(FS/FSA), COMBIN, MOTION, RELDISP, STRESS | modules | Ch. 3, Ch. 6.4–6.5 specs |
| NONLINEAR, COMB_XYZ_THD, COMB_XYZ_STRAIN, LOADGEN, SSI2ANSYS | option / aux modules | Ch. 3, 6, 9.17; Option A-AA manual (separate) |
| FILE1, FILE2, FILE3, FILE4 (`.n4`), FILE8, FILE81/FILE82, FILE9 (FILE9001…), FILE11, FILE12, FILE13, FILE14, FILE15, FILE73, FILE74, FILE77, FILE78, FILE88, FILE90, FILE91, COOXxxx, COOTKxxx, COOXI, COOTKI, DOFSMAP, COOSK, COOSM | inter-module files | Ch. 3 spec |
| `.pre`, `.err`, `.equ`, `.soi`, `.sit`, `.poi`, `.hou`, `.hounew`, `.map`, `.frc`, `.anl`, `.mot`, `.rdi`, `.str`, `.TFU`, `.TFI`, `.ACC`, `.RS`, `.TFD`, `.THD`, `.THS`, `.TH`, `.crv`, `.FMU` | file extensions | Ch. 3, 6, 9 specs |
| `Remove_Frequencies_from_FILE8.exe`, `COMB_XYZ_THD.inp` | aux tool / input | not specified in manual (see OQ) |
| INP, CHECK, AFWRITE, MDL, WRITE | general commands | Ch. 9.2 / 9.7 |
| FIXROT, FIXSHLROT, FIXSLDROT, FIXSPRROT (alias FIXSPROT), INTGEN, ETYPE, ETYPEGEN, MXI, D | model commands | Ch. 9.4 / 9.7 |
| EXCSTRCHK, FIXEDINT, HINGED | model checking | Ch. 9.8 |
| P, NONLINMOTDISP, NONLINBAT, BBC, BBCX, BBCY, BBCGEN, SHEAR, EDGE | Option NON commands | Ch. 9.13 / 9.17 |
| AMP | multiple-support amplification | Ch. 9.2.4 |
| Options/Check, View/Check Errors, Modules/Location, Options/Analysis | UI menus | Ch. 6.5.3, 6.6.1, 6.4.1, 6.5.4 |
| Interpolation options 0–6, smoothing parameter, phase adjustment | MOTION/STRESS inputs | Ch. 6 (MOTION, STRESS) |

---

## 15. Open questions / ambiguities

1. **Complex modulus form.** Chapter 1 says only "complex moduli". Is it `(1 − 2β² + 2iβ√(1−β²))` (SASSI/SHAKE91) or `(1 + 2iβ)`? Is it the same for soil layers (SITE), soil elements and structural elements (HOUSE)? Recommendation: the former, made selectable.
2. **Interpolation algorithms 0 to 5.** The manual gives names only: no window size, interpolant form, weights, or overlap/shift rules. The rational SDOF form `(A+Bω²)/(1+Cω²)` over 3-point windows is a reconstruction of classic SASSI and must be confirmed or calibrated against SASSI2000 output. Also undefined: what "dense" means in options 0 and 2, the weighting law in option 0, and how option 3 picks its "three windows".
3. **Option 6 "bicubic" spline.** Is it a cubic spline on Re/Im, on amplitude/phase, or truly bicubic (2D)? Which end conditions? How is it extrapolated below the first frequency and above the last? Recommendation: separate natural or not-a-knot cubic splines on Re and Im.
4. **Maximum FFT points.** 32,768 (§1.5.2, EQUAKE, SOIL) or 65,536 (MOTION, RELDISP, STRESS)? Recommendation: allow 65,536 in the response modules and 32,768 for input generation.
5. **Analysis frequency limit.** It is 500 for SITE, POINT, HOUSE, ANALYS and FORCE but 1,500 for MOTION, RELDISP and STRESS. This presumably allows several FILE8 runs to be combined (COMBIN). Confirm whether COMBIN output can exceed 500.
6. **FILE73 type.** It is drawn as a binary file (green) in Figure 1.1, but §1.5.4 calls it a "text file". The implementer must pick one and document it. Recommendation: text, per the prose.
7. **"20,0000" interaction nodes per embedment layer for incoherent analysis.** This is a typo, read here as 20,000.
8. **Build names.** IKTR8_650K vs IKTR9_650K. This is informational only.
9. **Translator input "CBD file".** Assumed to mean the ANSYS `.cdb` (CDWRITE) archive. The supported ANSYS element and keyword subset is not specified in Chapter 1.
10. **Option AA matrix import format.** The ANSYS matrix files used by HOUSEFSA (for example `.full`/HBMAT/MMF) are not specified. Option AA has its own separate manual, which is not available.
11. **FIXROT spring stiffness default "10".** This is an absolute number with model-unit dependence (force·length/rad). Should it be scaled to the units, or to the shell bending stiffness D = Et³/12(1−ν²)? The manual gives only the guidance "≤ 10% of the shell bending stiffness".
12. **Incompatible modes in SOLID.** What switches the 9 incompatible modes on: automatic for structural SOLIDs (ETYPE = 1), or a HOUSE option? Are they also used in PLANE?
13. **BEAM third node.** Presumably an orientation node (K). Confirm in the element-command spec. Also: is BEAM shear deformation included (Timoshenko) or not?
14. **POINT "radius".** The rule "average radius over all excavation solid elements" needs a definition of an element's radius. Equivalent plan-area radius `sqrt(A/π)` is assumed. Should it use only bottom-layer elements or all of them?
15. **Interaction-node-on-layer-interface tolerance.** Not given. Suggestion: 1e-4 × the minimum layer thickness.
16. **Nonlinear soil strain choices.** For "effective shear strain component": which component (γ_xz, γ_yz, or the maximum horizontal-vertical)? What is the effective strain ratio (0.65 or a magnitude-based value)? Is the SRSS combination of X, Y and Z effective strains really SRSS (per Ch. 4) for all three options? What are the convergence tolerances and maximum number of iterations for nonlinear soil SSI?
17. **Option NON hysteretic rules.** The CMS and TAK loop rules (unloading and reloading stiffness, pinching parameters), the definition of "effective" strain for linearization, how equivalent damping is evaluated from the loops, and the convergence norm ("Displacement/Force Convergence Error", Ch. 6) are not in Chapter 1. They must come from the Chapter 6/9.17 specs or the cited Ghiocel et al. 2017 paper. Without them, Option NON cannot be reproduced to match the original.
18. **Panel bending measure.** Figure 1.4 draws θ₁ and θ₂ on the top edge only. It is unclear whether curvature is (θ_top − θ_bottom)/H, (θ₁ + θ₂)/2 per unit height, or a rotation rather than a curvature. The units of the bending BBC (moment vs curvature, slope E·I) suggest curvature, as recommended in §8.2.
19. **Panel shear definition.** "Average vertical edge rotation after rigid-body extraction". How is the rigid rotation removed: from the bottom edge, or from a least-squares fit of all 4 corners? The recommended invariant formula in §8.2 adds the horizontal-edge rotation term. The two agree only if the rigid rotation is defined consistently.
20. **Figure 1.2 file names.** "Panelxxxx.THS → Panelxxxx_elastic.THD" looks like a typo for `_elastic.THS`. "Panel.FMU → Panelxxxx_elastic.FMU" mixes the per-panel and global naming. The exact contents of `.THD`, `.THS`, `.crv`, `.FMU` and `*_EQL_MATL_PROP*` for Option NON are not specified.
21. **COMB_XYZ_THD.inp format.** The user must write this file, but Chapter 1 gives no format; it points only to an example on the installation DVD. It is not the generic three-line batch `<module>.inp` file (model name / input file / output file) described in Chapter 3. Its content must be defined from the Chapter 3/4/6/9.17 specs or designed by us.
22. **Remove_Frequencies_from_FILE8 behaviour.** Interactive DOS tool with an unspecified interface. Our version needs a defined CLI.
23. **Quiet zone length rule.** The manual gives a qualitative rule only (it depends on the lowest damping). The decay criterion in §7 is our own.
24. **Gravity trick sign and placement.** "Best over the middle of the first half-period variation defined with a negative sign" is ambiguous about whether the first half-period is negative, that is a(t) = −g·sin(2πt/T). Recommendation: as in §8.2, and verify with test 5 of §13.
25. **Options/Check reference.** Chapter 1 cites §6.6.1, but the settings dialog is §6.5.3 and §6.6.1 is the error-viewer window. The default break numbers are not given.
26. **Command syntax details.** Chapter 1 does not say whether a space may replace the first comma (`EDGE 1,0,0,1`), whether quoted strings and comments are allowed, or how blank parameters are handled. See the Chapter 5/9 specs.
27. **AFWRITE partial writing.** "Affected file(s) are not written" implies a per-module error check. Must an existing stale file of the same name be deleted, to avoid running modules on outdated input? Recommendation: delete or rename it and warn.
28. **Seven vs six interpolation options.** The text says the "first six options (0–5)" cover SASSI 1982, SASSI2000 and four new schemes. The Ch. 6 dialog maps 0 = SASSI2000 and 1 = SASSI 1982. Use the Ch. 6 mapping as authoritative.
29. **Nonlinear spring "limited foundation sliding"** and **pile slipping in the vertical direction.** No friction or gap model is described. Is the BBC alone (with a plateau at the friction capacity) the intended model?
30. **Hardware, parallel and GPU statements, and RAM-driven case limits (500 load cases, 50 simulations).** These are not functional requirements for an educational Python version. Treat them as soft limits.
