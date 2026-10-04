# 05b — Analysis Options: HOUSE and FORCE modules

**Source:** ACS SASSI V3 User Manual, §6.5.4 *Analysis Options* — "HOUSE Module Options" and
"FORCE Module Options" (extracted-text lines 5547–6332; printed pages 129–146; PDF pages 131–148).
I viewed the two dialog screenshots (HOUSE: printed p.129, FORCE: printed p.145), the node-numbering
figure (p.131), equations 6.1–6.4 (p.136) and the `.pin` example window (p.143) as images and
transcribed them below.

Where this section depends on facts stated elsewhere in the manual, I cite them as
"(cross-ref §x / line n)". Those are context only; they are specified in full in other spec files.

**Tags:** **[CORE]** needed for computational correctness · **[IO]** file/input format ·
**[UI]** interface/plotting/convenience · **[ADV]** advanced (incoherency, nonlinear, multiple
excitation, Option AA).

**Provenance markers** used in this file:
- *(manual)*: stated explicitly in the manual.
- *(screenshot)*: read from a dialog image.
- *(inferred)*: my reconstruction to fill a gap. The implementer must confirm it or pick an
  alternative. Every inferred item is repeated in §13 *Open questions*.

---

## 1. HOUSE module: purpose, inputs, outputs

### 1.1 Purpose [CORE]
- HOUSE builds the **static stiffness and mass matrices** for two FE models that are combined in the SSI model:
  1. the **structure** FE model, which can include near-field or backfill soil modelled as structure, and
  2. the **excavated soil volume** FE model, i.e. the soil removed from the free field, as in the flexible-volume substructuring idea.
- Steps: read the nodes, nodal types (interaction or not), soil layer properties and element data
  for the structural and excavated-soil elements; form the element K and M; assemble them; store
  them in **compact format** ready for the SSI solution (ANALYS).
- **Interaction nodes belong to the excavated soil volume.**
  - **FV** (Flexible Volume): every node of the excavated volume is an interaction node.
  - **FI** (Flexible Interface, either FI-EVBN or FI-FSIN) and **FFV** (Fast Flexible Volume): only part of the excavation nodes are interaction nodes.
  - **Surface foundation:** no excavation volume, interaction nodes only at ground level, so FV = FI = FFV.
- Structural nodes are **not** interaction nodes, with one exception: nodes shared by the
  structure basement and the excavated-soil model on the lateral and bottom excavation surfaces.
- **Incoherent analysis:** HOUSE also builds the free-field **coherency matrix**, computes its
  **eigen-solution** (spectral factorization) and from it computes the free-field incoherent motion
  variations for **each frequency**. These are saved to **FILE77**.
- **Nonuniform / multiple-support excitation:** the data also goes to **FILE77** (cross-ref line 1919).

### 1.2 Module variants [CORE][ADV]
| Name | Use |
|---|---|
| `HOUSE` | baseline (standard solver) |
| `HOUSEFS` | fast-solver version; the node-numbering optimizer is **optional** |
| `HOUSEFSA` | Option AA (Advanced ANSYS) version; the optimizer is used **automatically** |

### 1.3 Files [IO]
| File | Dir | Content |
|---|---|---|
| `modelname.hou` | in | HOUSE input text file written by `AFWRITE`. If line 1, column 1 holds `1`, the node optimizer is enabled. |
| `modelname.sit` | in | SITE input; must be in the working directory (needed for embedded models: layer data) (cross-ref line 1929) |
| `modelname.pin` | in | Nonlinear-soil initial-property file (§6). Edited through the "Input Data" button. |
| `modelname.liq` | in/out | Written by HOUSE on the initial run. If it is non-empty and contains `1`, HOUSE takes model materials from FILE74 on the next iteration (cross-ref line 2350). |
| `FILE74` | in | Effective strains and updated properties from STRESS (nonlinear iterations) |
| `COHXUSER`, `COHYUSER`, `COHZUSER` | in | User-defined coherency tables for X, Y, Z (no extension) (§3.16) |
| `FREQCOH`, `DISTCOH` | in | Frequency vector and distance vector for the user tables (no extension) |
| `COOSK`, `COOSM` | out | Binary assembled stiffness and mass matrices (compact storage) |
| `FILE4` = `modelname.N4` | out | Topology: node coordinates and element connectivity (plus K/M for the baseline solver, cross-ref line 1914) |
| `FILE77` | out | Incoherency eigen-solution / incoherent motion variations, or multiple-excitation data. One file per simulation: `FILE77001` … `FILE77050`. Each file holds X, Y and Z. |
| `FILE78` | out | Non-empty only for nonlinear SSI; read by STRESS |
| `modelname.hounew` | out | Renumbered (optimized) model input. The manual also misspells it ".hownew". |
| `modelname.map` | out | Pairs of old and new node numbers |
| HOUSE output listing | out | Echo of the model, soil-layer assignment of excavation elements, and the incoherent mode contribution table (search string `I N C O`) |

### 1.4 FE library handled by HOUSE [CORE] (see GROUP command)
| Type | Description (manual) | DOF/node (cross-ref §1.5.3) |
|---|---|---|
| `SOLID` | 3D solid, 8 nodes (degenerate pyramids allowed), with or without incompatible modes (9 incompatible modes, structure only) | 3 translations |
| `BEAMS` | 3D beam, 3 nodes (3rd = orientation node) | 3 translations + 3 rotations |
| `SHELL` | 3D thin plate, Kirchhoff theory, 4 nodes (triangles allowed) | 6 |
| `TSHELL` | 3D thick plate, Mindlin–Reissner theory, 4 nodes | 6 |
| `PLANE` | 2D plane-strain solid, 4 nodes (triangles allowed) | 2 translations |
| `SPRING` | 3D spring, 2 nodes | 6 |
| `GENERAL` | 3D general stiffness/mass matrix element. 2 nodes = global coordinates; 3 nodes = local coordinates (see MXI) | 6 |

The table marks SHELL with `*` and TSHELL with `**`, but the page has no footnote text (see §13).

Cross-referenced conventions HOUSE must honour [CORE]:
- Mass matrix: 50 % lumped + 50 % consistent, except beams (consistent) and plates (lumped)
  (line 965).
- Material damping: complex moduli, which gives frequency-independent damping that can differ from
  element to element (line 969). Inputs are the P-wave and S-wave damping ratios of `M` and `L`.
- For `SHELL` materials, `pdamp` must equal `sdamp` (M command note).

### 1.5 HOUSE modelling remarks, recast as rules for the implementation
| # | Rule | Tag | Implementation action |
|---|---|---|---|
| R1 | Thin `SHELL`: the drilling rotation (about the plate normal) is singular. The user should apply `FIXROT` or `FIXSHLROT` (small rotational stiffness). `TSHELL` handles this automatically inside HOUSE. | CORE | TSHELL element: add a small drilling stiffness internally. SHELL: warn if there is no FIXROT/FIXSHLROT. |
| R2 | SOLID incompatible modes can be switched on from the UI; they never apply to excavated-soil SOLIDs. TSHELL integration (Reduced or Selective) is chosen per element with `EINT` (0 = Reduced 1×1 bending+shear [default]; 1 = Selective 2×2 bending / 1×1 shear). | CORE | Force incompatible modes OFF for excavated-soil SOLIDs. |
| R3 | The excavated soil is modelled only with 3D `SOLID` or 2D `PLANE`. **Every interaction node below ground must lie on a soil-layer interface.** | CORE | Check: each sub-surface interaction node elevation equals a layer interface elevation (within tolerance). Error if not. |
| R4 | Suggested **node numbering is bottom-up**: start at foundation level and go layer by layer to the surface (figure p.131: two cubes with numbering 1…90 and 1…88 ascending from bottom). **Interaction nodes must always be in ascending order from the deepest excavation level to the ground surface** (WARNING). | CORE | Validate interaction-node order against elevation. Error/warning if a lower level has higher numbers than an upper level. The optimizer (§3.7) may renumber. |
| R5 | Element numbering: there is no restriction without embedment, but it must be continuous and bandwidth-friendly. **With embedment, excavation elements are numbered top-down (surface → base)**, so their associated far-field layer numbers increase. | CORE | Check that layer index is monotone non-decreasing with element number in excavation groups. |
| R6 | Strongly recommended, and a **strict rule if stress contours are needed**: one excavation **group per embedment layer**, each a horizontal slab matching one far-field layer. Number of excavation groups = number of embedment layers. Groups numbered surface → foundation. Each group's soil layer is set with `L` using the same number as the far-field layer. | CORE/UI | Validation warning; strict error when STRESS contour output is requested. |
| R7 | Always review the HOUSE output. In particular, check that the excavation SOLIDs got the correct soil-layer numbers before running the SSI analysis. | UI | HOUSE output must list, per excavation element (or group), the assigned layer number and its properties. |
| R8 | Demo 5 (`Demo5.pre`) is the reference example for an embedded model. | — | Use as a regression case if it becomes available. |

---

## 2. HOUSE dialog: transcription [UI]

Window "Analysis Options". Tabs: `EQUAKE | SOIL | SITE | POINT | HOUSE | FORCE | ANALYS | MOTION | STRESS | RELDISP | NONLINEAR | AFWRITE`. Buttons: **Ok**, **Cancel**.
Values below are those **shown in the screenshot**. They are an example state, **not necessarily defaults** (see §13).

### 2.1 Left column
| Group box | Control (exact label) | Widget | Screenshot value | Command mapping |
|---|---|---|---|---|
| Operation Mode | `Solution` / `Data Check` | radio | Solution | `HOUSE <opmode>`: 0 = complete solution, 1 = data check only |
| Dimension of Analysis | `1D` / `2D` / `3D` | radio | 3D | `HOUSE <dim>`: 0 = 1D (not available), 1 = 2D, 2 = 3D |
| Flexible Volume Method | `Flexible Volume(FV)` / `Fast Flexible Volume(FFV)` / `Flexible Interface(FI)` | radio | Flexible Interface(FI) | `HOUSE <imp>`: 0 = FV, 1 = FFV, 2 = FI |
| — | `Acceleration of Gravity` | text | 32.2 | `HOUSE <gravity>` (also `GRAVITY` command) |
| — | `Ground Elevation` | text | -10 | `HOUSE <gelev>` (also `GROUNDELEV`) |
| — | `Non-Linear SSI` + button `Input Data` | button | — | opens `modelname.pin` in a text editor (§6) |
| — | `Optimize Model` | checkbox | unchecked | writes `1` at `.hou` line 1, column 1 |

### 2.2 Middle column
| Group box | Control (exact label) | Widget | Screenshot value | Command mapping |
|---|---|---|---|---|
| Soil Motion | `Coherent` / `Incoherent` | radio | Incoherent | `HOUSE <coh>`: 0 = coherent, 1 = incoherent |
| — | `Coherence Parameter X Dir` | text | 0.1 | `INCOH <gammax>` |
| — | `Coherence Parameter Y Dir` | text | 0.1 | `INCOH <gammay>` |
| — | `Coherence Parameter Z Dir` | text | 0.2 | `INCOH <gammaz>` |
| — | `Alpha Directionality Factor` | text | 1000 | `INCOH <alpha>` (and possibly `WPASS <cohf>`, see §13) |
| — | `Number of Embedded Layers` | text | 0 | `INCOH <ngp>` ("number of mesh points per each embedment level", see §13) |
| — | `Time Step of Sesmic Motion` (sic) | text | 0.005 | mirrors SITE `<delt>` |
| — | `Nr. of Fourier Components` | text | 4096 | mirrors SITE `<nft>` |
| — | `Frequency Set Number` | text | 1 | mirrors SITE `<freq>` |
| — | `Number of Incoh. Modes` | text | 0 | `INCOH <nmodes>` |
| — | `Print Coherency Matrix` | checkbox | unchecked | `INCOH <ipr>` (the text calls it "Incoherent Mode Contributions") |

### 2.3 Right column
| Group box | Control (exact label) | Widget | Screenshot value | Command mapping |
|---|---|---|---|---|
| Multiple Excitation | `Use Multiple Excitation` | checkbox | unchecked | `HOUSE <me>`: 0/1 |
| Multiple Excitation | `Input Motion Number` | spin box (▲▼) | 1 | `ME <no>` (1…10) / `AMP <no>` |
| Multiple Excitation | `First Foundation Node` | text | 1 | `ME <nfirst>` |
| Multiple Excitation | `Last Foundation Node` | text | 69 | `ME <nlast>` |
| Multiple Excitation | `X Coord. of Control Point` | text | 0 | `ME <xc>` (not used) |
| Multiple Excitation | `Y Coord. of Control Point` | text | 0 | `ME <yc>` (not used) |
| Multiple Excitation | `Z Coord. of Control Point` | text | 0 | `ME <zc>` (not used) |
| Spectral Amplifcation (sic) | multi-line text box with scroll bar | text | `1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1` | `AMP <no>,<a1>,…` |
| (below the box) | `Use Complex Spectral Amp.` | checkbox | unchecked | `HOUSE <cmplxspec>`: 0/1 |
| (below the box) | `Non-Uniform Motion` | checkbox | **disabled (greyed)** | none documented |
| (below the box) | `Non-Uniform Soil` | checkbox | **disabled (greyed)** | none documented |

### 2.4 Lower area
| Group box | Control (exact label) | Widget | Screenshot value | Command mapping |
|---|---|---|---|---|
| Wave Passage | `Use Wave Passage` | checkbox | **checked** | `HOUSE <wpass>`: 0/1 |
| Wave Passage | `Apparent Velocity for Line D` | text | 1e+008 | `WPASS <appv>` |
| Wave Passage | `Angle Line D with X Axis` | text | 0 | `WPASS <ang>` (degrees) |
| Wave Passage | `Unlagged Coherency Model` | text (integer) | 3 | probably `WPASS <cohf>` *(inferred)* |
| Motion Incoherency Simulation | `Deterministic (Median) Incoherency Input` / `Stochastically Simulated Incoherency Inupe` (sic) | radio | Deterministic | implied by `INCOH <HSeed>,<VSeed>,<RandPhz>` (all zero = deterministic) |
| Superposition Mode | `Linear` / `Quadratic` | radio | Linear | no documented argument (see §13) |
| — | `Ansys Model Input` | checkbox | unchecked | enables Option AA (`HOUSEFSA`) |
| Ansys Model Type | `Embedded` / `Surface` | radio, disabled unless Ansys Model Input is checked | Embedded | `ANSYSMODELTYPE <type>`: 1 = embedded, 2 = surface |

The text also mentions stochastic-simulation inputs that are **not visible** in the screenshot:
a Horizontal SEED, a Vertical SEED, a Random Phase angle and the Number of Simulations.
An implementation should show them in the "Motion Incoherency Simulation" group, enabled only when
"Stochastically Simulated" is selected [UI] *(inferred layout)*.

### 2.5 UI enable/disable logic (recommended; partly inferred) [UI]
- Coherent selected: disable every incoherency field (coherence parameters, alpha, embedded layers,
  incoherent modes, print, simulation, superposition). Wave Passage and Multiple Excitation can still be used *(inferred)*.
- Dimension = 2D: force Coherent and disable Incoherent (incoherency needs a full 3D model, §3.8).
- Unlagged Coherency Model ∈ {2…7}: requires `Use Wave Passage` (§3.15). Either auto-check it or
  block OK with a message.
- Use Multiple Excitation requires `Use Wave Passage` (§3.19).
- Ansys Model Type is enabled only when `Ansys Model Input` is checked.
- Time step, NFFT and frequency set mirror the SITE tab. Keep a single source of truth: edit in one
  place and reflect the value in both.

---

## 3. HOUSE fields: semantics, ranges, rules

### 3.1 Operation Mode [CORE]
- `Solution`: full run.
- `Data Check`: read and check the input only.
- `HOUSE <opmode>`: 0 = complete solution, 1 = data check only.

### 3.2 Dimension of Analysis [CORE]
- Allowed values: 2D or 3D. 1D is shown in the dialog but is not available.
- 2D: the POINT module must be **POINT2**, not POINT3. Keep the HOUSE and POINT choices consistent and check them together.
- `<dim>`: 0 = 1D, 1 = 2D, 2 = 3D.

### 3.3 Flexible Volume Method [CORE]
- Choices: **FV**, **FFV** or **FI**. FI covers both FI-EVBN and FI-FSIN; which one depends on how the interaction nodes were generated.
- `<imp>`: 0 = FV, 1 = FFV, 2 = FI.
- Surface foundation: all three methods are identical.
- Interaction nodes are generated automatically with `INTGEN,<type>,[level skip]` (cross-ref §9.7.19):
  - 0 = clear
  - 1 = FV
  - 2 = FI-EVBN (MSM)
  - 3 = FI-FSIN (SM)
  - 4 = surface
  - 5 = FFV with multiple internal interaction-node layers; `level skip` default 1

  The older `INT` command still works.
- **WARNING (regulatory):** whenever FI or FFV is used, sensitivity studies against FV are required
  (ASCE 4-16, USNRC SRP 3.7.2).
  [UI] Show a non-blocking warning when imp ≠ 0 and the model is embedded.

### 3.4 Acceleration of Gravity [CORE]
- Same value as in the SITE tab.
- Typical values: 32.2 ft/s² (British units) or 9.81 m/s² (SI) (cross-ref §1.5.5).
- **Error 1** if ≤ 0.
- HOUSE uses it to turn specific weights into mass densities: ρ = γ_w / g. Masses can be given in weight units (MR, MUNITS).
- Unit system: the whole model must use consistent units. HOUSE does not convert anything except
  in the coherency models (see `met`, §3.9).

### 3.5 Ground Elevation and element classification [CORE]
- Elevation of the ground surface, in model length units.
- HOUSE also uses it to classify SOLID and PLANE elements as **structure** or **excavated soil**, together with `ETYPE` (cross-ref §9.4.11).

Classification algorithm:
```
for each element e in SOLID/PLANE groups:
    if ETYPE(e) == 1: structure (material index -> material table M)
    elif ETYPE(e) == 2: excavated soil (MSET index -> far-field SOIL LAYER table L)
    else (ETYPE 0, default):
        if e lies below ground elevation: excavated soil (index -> L table)
        else: structure (index -> M table)
SHELL/TSHELL: default structure; ETYPE 2 = embedded shell (its nodes become interaction nodes via INTGEN)
BEAMS: ETYPE 2 has no effect
```
- `ETYPE,1` is only for element groups that belong to the structure FE model. That model can contain a near-field backfill-soil part.
- What "below ground surface" means is unspecified: centroid below, all nodes ≤ g_elev, or something else (§13). Recommended: every node at or below g_elev + tol, with the centroid strictly below.

### 3.6 Non-Linear SSI "Input Data" button [UI][ADV]
- Opens the `.pin` file for editing (§6).
- The manual text also describes an "Input Data" button for the **pile interface**, with the note "PINT module is not available". The dialog shows the button next to the "Non-Linear SSI" label, so treat it as the `.pin` editor (§13).

### 3.7 Optimize Model (node-numbering optimizer) [CORE][IO]
- Checkbox. It writes `1` in line 1, column 1 of the `.hou` file.
- Optional in HOUSEFS. Always on in HOUSEFSA (Option AA).
- Purpose: faster and less RAM for embedded SSI.
- The manual "strongly recommends" it for any larger embedded model.
- Outputs:
  - `modelname.hounew`: the renumbered model. FILE4, COOSK and COOSM are built from it.
  - `modelname.map`: the old-to-new node number pairs.
- **Post-processing consequence:**
  - MOTION and RELDISP node requests must use the **new** numbers.
  - STRESS reports use the new numbering.
  - FILE8 is based on the renumbered model.
- Validation case: V&V Problem 31.
- The renumbering algorithm is not specified (§13). It must keep the interaction-node bottom-up ascending order (R4).

### 3.8 Soil Motion: Coherent / Incoherent [CORE][ADV]
Incoherence applies **only** to 3D models with **no symmetry** (full models):
- not 2D,
- not half or quarter models (`SYMM` must be absent).

The interaction-node numbering must start at the bottom layer (baserock side) and go up to the ground surface.

The six incoherent SSI approaches (five deterministic, one stochastic), and how each is selected:

| # | Approach | HOUSE selection | Downstream | EPRI 2007 status |
|---|---|---|---|---|
| 1 | **Stochastic Simulation (SS)**, "Simulation Mean" | Stochastically Simulated; non-zero HSeed and VSeed; RandPhz = 180; Nsim ≤ 50; all modes (`nmodes=0`) | ANALYS runs all simulations (FILE8xxx); MOTION averages | reference approach (validated) |
| 2 | **AS / Linear with ATF phase adjustment** | Deterministic; Linear; `nmodes` = 0 or n | single ANALYS run | validated |
| 3 | AS / Linear without phase adjustment | Deterministic; Linear | single run | not used by EPRI |
| 4 | **SRSS TF, zero phase** (Quadratic) | Deterministic; Quadratic; one HOUSE+ANALYS run per mode with `nmodes = −k` | MOTION with `SRSSTF.txt` phase flag 0 | validated |
| 5 | SRSS TF with coherent phase | as #4 | `SRSSTF.txt` flag 1; first FILE8 listed = coherent run | not validated |
| 6 | **SRSS FRS** | Deterministic; **Linear**; one run per mode (`nmodes = −k`) | MOTION/STRESS per mode; the user SRSS-combines the end results (ATF, ISRS, ZPA) | not validated, not endorsed by ISG-01; benchmarks only |

How the AS phase adjustment is switched on or off is not documented (§13).

Guidance to put in the UI help and as warnings:
- Deterministic approaches should be limited to **rigid** foundations (stick models with rigid basemats).
- **Flexible** foundations should use **stochastic simulation**.
- SRSS approaches need tens to hundreds of modes for flexible mats. That makes them impractical; they exist for benchmarking.
- Before using SRSS, establish the number of modes with the 90 % cumulative-variance criterion (§3.14).
- The EPRI zero-phase approaches may be over-conservative in the mid-frequency range and unconservative at high frequency when only a few modes are used. In the vertical direction, 40–50 modes may still be too few.

### 3.9 Coherence Parameters X / Y / Z [CORE][ADV]
- Used by the **Luco–Wong** model (Type 1). There is one parameter per motion component: X and Y horizontal, Z vertical.
- Typical range: **0.10–0.30**. Higher values give upper bounds on incoherence effects.
- The horizontal field can be isotropic (γx = γy) or anisotropic.
- **Error 58:** the coherence parameter must be **≥ 0.1**.
- `INCOH` also has `<met>` (0 = British, ft; 1 = SI, m). It matters because the Abrahamson models are empirical in metres and Hz. The dialog shows no unit control (§13).

Luco–Wong model *(standard published form, not printed in this section)*:
`γ_LW(ω, D) = exp[ −( γ_c · ω · D / Vs )² ]`, with ω = 2πf (rad/s), D = separation distance,
Vs = mean soil shear-wave velocity (the "Alpha" field when Type 1 is used, see 3.10),
and γ_c = γx, γy or γz depending on the motion component.

### 3.10 Alpha Directionality Factor [CORE][ADV]
- **Abrahamson and user-defined models (Types 2–7):** α ∈ [0, 1].
  - α = 0.50: isotropic (radial) coherency, circular in plan.
  - α ≠ 0.50: anisotropic (directional), elliptical.
- Directional distance between interaction nodes i and j (manual equation, "homotopic relationship"):

  `D_ij = sqrt( 2 · ( α · DX_ij² + (1 − α) · DY_ij² ) )`

  - DX and DY are the directional (horizontal) separations, measured along the principal axes of the coherency ellipse.
  - The ellipse can be rotated so its axes line up with Line D (the angle in §3.15) *(the manual says "could be rotated")*.
  - α = 0.5 gives `D = sqrt(DX² + DY²)`, the Euclidean distance.
  - α = 0.1 gives `D = sqrt(0.2 DX² + 1.8 DY²)`, so Y distances are weighted 3× more than X distances (√(1.8/0.2) = 3).
  - α = 0.9 is the reverse.
- **Luco–Wong (Type 1):** the same field holds the **mean soil shear-wave velocity Vs**, in model velocity units.
  - **Error 59** if ≤ 0.
  - The screenshot value 1000 fits a Vs value.
- **Error 114** ("Illegal Directional Coherence Factor") if < 0.
- The manual gives two contradictory readings of α = 0.1:
  - the dialog description says the ellipse "decays three times faster in X";
  - the paragraph with the formula says "Y-distances weighted 3× more", which means faster decay in Y.

  **The formula is authoritative** (§13).
- EPRI TR-1015111 (2007) validation used the **2005 Abrahamson isotropic** model. The 2007 models came out after those studies.

### 3.11 Number of Embedded Layers [CORE][ADV]
- The number of embedment layers defined by the set of interaction nodes.
- FI-EVBN or FI-FSIN: internal nodes must **not** be interaction nodes. FV: all excavation-volume nodes are interaction nodes.
- The batch argument is `INCOH <ngp>`, described as "number of mesh points per each embedment level".
- **Error 60**: "Illegal Number of Mesh Points / Embedment Level" if ≤ 0.
- The screenshot shows 0, which conflicts with Error 60 if the check were unconditional. Its exact role in the algorithm is unspecified (§13).

### 3.12 Time Step of Seismic Motion / Nr. of Fourier Components / Frequency Set Number [CORE]
- All three repeat the SITE tab values. A single value should be shared.
- **Δt** in seconds. **Error 49** if < 0.
- **NFFT** must be a power of 2; otherwise the UI rounds to the nearest power of 2 and warns (SITE rule). Maximum 32 768. **Error 50** if < 0.
- **Frequency Set Number** is a FREQ set index.
  - Frequencies are `f_j = n_j · Δf` with `Δf = 1/(Δt·NFFT)`.
  - The frequency numbers n_j are positive integers, sorted ascending; duplicates stop the run (cross-ref SITE, line 5351–5364).
  - **Error 44** if the set is undefined. **Error 120** if the set is empty.
- HOUSE needs the frequency list because the coherency matrix and its eigen-solution depend on frequency. They are computed at **each SSI frequency**.

### 3.13 Number of Incoherency Modes (`nmodes`) [CORE][ADV]
| Input | Meaning |
|---|---|
| `0` | **all** modes (N = number of interaction nodes). This is the default and is recommended for SS and AS. |
| `n > 0` | use the first n modes only (n < N) |
| `−n < 0` | use **only mode n**. This is used with the Quadratic/SRSS approaches, which need one SSI run per mode. |

- Using all modes gives the best accuracy at negligible extra run time.
- Quadratic (SRSS TF) also needs `SRSSTF.txt` for MOTION.
- Stochastic simulation always uses all modes (one per interaction node per direction).
- Mode ordering for "first n": by **decreasing eigenvalue** *(inferred, see §13)*.

### 3.14 Print Coherency Matrix / Incoherent Mode Contributions (`ipr`) [CORE][UI]
- When checked, HOUSE prints, for every frequency, the percent contribution of each incoherent mode to rebuilding the free-field coherency matrix.
- Purpose: check the accuracy of the random-field decomposition.

Equations (transcribed from p.136). N = number of interaction nodes, λ_j = eigenvalues of the coherency matrix Σ (Σ = Φ Λ Φᵀ). Because the diagonal of Σ is 1 (zero distance means coherency 1), `trace(Λ) = trace(Σ) = N`.

```
(6.1)   Σ_{j=1..N} λ_j = N                         (all modes: total variance recovered)
(6.2)   Σ_{j=1..m} λ_j < N        for m < N        (truncated)
(6.3)   υ_j = (λ_j / N) · 100                      (percent contribution of mode j)
(6.4)   Σ_{j=1..m} υ_j = Σ_{j=1..m} (λ_j / N) · 100   (cumulative percent, m modes)
```
- Each λ_j is the variance carried by mode shape j. N is the total variance of the incoherent amplitude-variation field at that frequency.
- **Convergence criterion** (recommended): pick m so the cumulative contribution is ≥ **90 %** at each frequency. This is the analogue of the 90 % cumulative modal-mass rule. It may be conservative for extremely rigid foundations.
- **Accuracy check WARNING:** verify that (6.1) holds at every frequency. If it does not, the coherency matrix is numerically ill-conditioned, most likely because of the spatial layout of the interaction nodes, e.g. near-coincident horizontal projections.
- The manual text calls these equations "5.1" and "5.4", which are typos for 6.1 and 6.4.
- Output format [UI][IO]: the table header must contain the string **`I N C O`**, the documented search key. For each frequency, list j, λ_j, υ_j and the cumulative Συ.

### 3.15 Wave Passage group [CORE][ADV]
- **Use Wave Passage** (`HOUSE <wpass>` 0/1). **Required** for:
  - the Abrahamson models (Types 2–6),
  - the user-defined model (Type 7),
  - Multiple Excitation.

  The Luco–Wong model works with or without it.
  WARNING: the Abrahamson and user-defined models are applied **only** when wave passage is checked.
- **Apparent Velocity for Line D** (`WPASS <appv>`): apparent horizontal propagation speed along Line D.
  - **Error 113** if ≤ 0.
  - To make wave passage negligible, use a very large value, e.g. **1.0e9** (EPRI 2007 rock-site practice). The screenshot shows 1e+008.
- **Angle Line D with X Axis** (`WPASS <ang>`): angle in degrees of the horizontal Line D from global X.
  - Line D sets the wave-passage direction.
  - It also sets the principal axes of the directional coherency models (Types 2–7).
- **Unlagged Coherency Model** (integer 1–7):

| Type | Model |
|---|---|
| 1 | 1986 Luco–Wong (theoretical, not validated against records) |
| 2 | 1993 Abrahamson, all soil types (similar to 2005) |
| 3 | 2005 Abrahamson, all sites, surface foundations |
| 4 | 2006 Abrahamson, all sites, embedded foundations |
| 5 | 2007 Abrahamson, hard-rock sites, surface and embedded foundations |
| 6 | 2007 Abrahamson, soil sites, surface foundations (also shallowly embedded) |
| 7 | User-defined coherency functions for X, Y, Z (files, §3.16) |

  - Types 2–7 can be rotated so their principal axes follow Line D.
  - Directional variants use α ≠ 0.5 (e.g. 0.1 or 0.9). Radial variants use α = 0.5.
  - The Abrahamson models are empirical: no parameters beyond f, D and α.
  - **Type 6 WARNING:** valid for surface or shallowly embedded foundations only. Deeply embedded structures in soil (e.g. SMRs) need a **multilevel** incoherent analysis (§7).

Wave-passage phase *(inferred, standard; the manual gives no formula)*:
- With the e^{+iωt} convention, a node at plan position (x_i, y_i) gets the arrival delay `τ_i = (x_i cosθ + y_i sinθ)/V_app` and the phase factor `exp(−iωτ_i)`.
- Lagged coherency: `γ_ij,lagged = γ_ij,unlagged · exp(−iω(τ_i − τ_j))`.
- The delay reference point (τ = 0) is not specified (§13).

### 3.16 User-defined coherency files (Type 7) [IO][ADV]
| File | Content |
|---|---|
| `FREQCOH` | frequency points vector (Hz), length NF |
| `DISTCOH` | relative distance points vector, length ND (model length units, see `met`) |
| `COHXUSER` | coherency matrix for X motion, NF × ND |
| `COHYUSER` | coherency matrix for Y motion, NF × ND |
| `COHZUSER` | coherency matrix for Z motion, NF × ND |

- No file extensions.
- Default and UI-imposed size: **100 × 100** (100 frequencies × 100 distances).
- Smaller sizes are allowed if the integer sizes are written into the `.hou` file. V&V Problem 39 uses 61 frequencies × 81 distances. The UI cannot change the size.
- The frequency and distance ranges must cover the analysis frequency range and the foundation size.
- Text layout, the orientation of rows (frequencies) versus columns (distances), and interpolation are unspecified (§13). Recommended: bilinear interpolation in (f, D), clamped at the ends, with the coherency forced to 1 at D = 0.

### 3.17 Motion Incoherency Simulation [CORE][ADV]
- **Deterministic (Median/Mean) Incoherency Input**: zero phase angles between the spatial wavelength components (modes). Selected when **HSeed = VSeed = 0 and RandPhz = 0**.
- **Stochastically Simulated Incoherency Input**: random phase angles uniform in **[−RandPhz, +RandPhz]**, with RandPhz = **180°**. Requires:
  - arbitrary non-zero **SEED** integers for Horizontal (`HSeed`, used for X and Y) and Vertical (`VSeed`, used for Z),
  - a **Number of Simulations**: 1 to **≤ 50**, all run in a **single HOUSE run**.
- Outputs: `FILE77001 … FILE77050`. ANALYS then produces up to 150 FILE8s (3 directions × 50). See the ANALYS spec, "Simultaneous Cases"; the number of simulations must match between HOUSE and ANALYS.
- Final SS response: statistical mean of the per-sample SSI responses (done downstream).

### 3.18 Superposition Mode [CORE][ADV]
- **Linear**: algebraic sum (AS) of the scaled incoherent spatial modes. Also used for SRSS FRS one mode at a time.
- **Quadratic**: SRSS of the modal complex ATF amplitudes (SRSS TF). One run per mode is needed (`nmodes = −k`).
- The SRSS TF phase option (0 = zero phase, 1 = coherent phase) is set in MOTION's `SRSSTF.txt`, not in HOUSE.
  - The manual names this file "STRSSTF.txt" in one place. The correct name is `SRSSTF.txt` (cross-ref MOTION, line 6696).
  - Format: line 1 `[# of modes] [phase option]`, then one FILE8 name per line. If phase option = 1, the first name is the coherent FILE8.

### 3.19 Multiple Excitation (nonuniform input) [CORE][ADV]
- **Use Multiple Excitation** (`HOUSE <me>`).
  - **Requires Wave Passage** to be on.
  - Applies to a single continuous foundation divided into zones, or to several separate foundations.
  - Nonuniform amplitude is introduced as a variable free-field Fourier amplitude under the foundation.
  - Can be combined with incoherence and wave passage.
  - A single foundation can be partitioned into up to **5000 zones**. `ME <no>` allows **1–10**, which is a conflict (§13).
- **Input Motion Number** (spin box): which zone or motion the fields below refer to.
- **First / Last Foundation Node**: the range of **interaction nodes** belonging to that foundation or zone.
  - Each range must be **contiguous with unit increment**. Do not skip numbers between foundations.
  - The bottom-up numbering rule still applies.
  - **Error 116 / 117**: illegal first or last node for motion i.
  - **Error 115**: option checked but no data defined.
- **X/Y/Z Coord. of Control Point**: reserved, **not used**.
- **Spectral Amplification Ratios (SAR)** (text box; `AMP,<no>,<a1>,…<a100>`; `a1 = 0` deletes the list):
  - Definition, at each SSI frequency: `SAR(ω_j) = ATF_local(ω_j) / ATF_ref(ω_j)`. ATF_local is the complex ATF of the "local" motion for that foundation or zone; ATF_ref is the "reference" motion from the single-control-motion SSI analysis.
  - **Count = number of SSI frequencies** in the selected set (**Error 119** otherwise).
  - Separators: blank or comma.
  - **Error 118** if any value is outside **[0, 10]**.
  - Real values by default. With **Use Complex Spectral Amp.** (`<cmplxspec>` = 1), complex SAR carry the phase effects of differential motion. The complex input syntax is unspecified (§13).
  - SAR usually differ for X, Y and Z inputs because of the layering under the foundation. In that case, write **separate `.hou` files per direction** and run HOUSE + ANALYS (restart) per direction.
  - V&V Problem 35 is the example.
- Application *(inferred)*: for interaction nodes in zone k,
  `U_ff,i(ω_j) ← SAR_k(ω_j) · U_ff,i(ω_j)`, applied after any incoherent synthesis and wave-passage phase. Nodes in no zone keep SAR = 1.
- `Non-Uniform Motion` / `Non-Uniform Soil` checkboxes: disabled in this version and undocumented. Implement as disabled placeholders [UI].

### 3.20 ANSYS Model Input / ANSYS Model Type (Option AA) [ADV][IO]
- Checkbox: the ANSYS structural model is used directly for SSI (Option AA, **fast-solver version only**, module `HOUSEFSA`).
- Then choose **Embedded** or **Surface**, either here or with `ANSYSMODELTYPE,<type>` (1 = embedded, 2 = surface). The UI selection is recommended.
- All ANSYS model files must be copied into the SSI working directory. See the separate "ACS SASSI-ANSYS Integration Capability" manual.
- The node optimizer is applied automatically, and `modelname.map` is written.
- **WARNING:** with Option AA, the `.hou` from `AFWRITE` is a **placeholder**: topology and connectivity only, with fake materials and constants. It is **not runnable** by regular HOUSE. The implementation must mark such files (e.g. a header flag) and refuse to run them in non-AA mode.
- Related (cross-ref §1.5.3): Option AA can convert ANSYS MATRIX50 super-elements to GENERAL (GM) elements.

---

## 4. Incoherent free-field computation in HOUSE: algorithm [CORE][ADV]

The manual describes this at concept level only. The steps below are a coherent reconstruction.
Each step is labelled manual-stated (M) or inferred (I).

1. (M) Gather the N interaction nodes (excavation-volume nodes per FV/FI/FFV, or the surface nodes).
   Use their **horizontal projections** (x_i, y_i). The coherency matrix is built on projections
   (cross-ref line 2750).
2. (M/I) For each SSI frequency f_j (from the frequency set) and each motion component c ∈ {X, Y, Z}:
   1. Distances: rotate (x, y) into the Line-D frame (θ = angle of Line D). Then
      `D_ik = sqrt(2(α ΔX'² + (1−α) ΔY'²))` for Types 2–7. For Type 1 (Luco–Wong), use the plain
      horizontal distance and Vs = Alpha field *(I)*.
   2. Unlagged coherency matrix: `Σ_ik = γ_model,c(f_j, D_ik)`. Σ is real, symmetric and positive
      semi-definite, with `Σ_ii = 1`. Abrahamson models have separate horizontal and vertical forms:
      the horizontal form is used for X and Y, the vertical form for Z *(I)*. Type 7 uses
      COHXUSER, COHYUSER or COHZUSER.
   3. Eigen-solution `Σ = Φ Λ Φᵀ`, with λ sorted descending *(I)* and negative round-off clipped to 0 *(I)*.
   4. Check `Σλ = N` (6.1). If `ipr`, print υ_j and the cumulative (6.3, 6.4).
   5. Keep the modes given by `nmodes` (0 = all; n = first n; −n = only mode n).
3. (M/I) Synthesize the incoherent amplitude-variation vector `s_c(f_j)` (length N) and write it to FILE77:
   - **Stochastic** (sample r = 1…Nsim): `s = Φ_m Λ_m^{1/2} η`, with `η_k = exp(iθ_k)` and θ_k ~ U[−RandPhz, +RandPhz] (degrees → radians).
     - HSeed seeds X and Y; VSeed seeds Z *(I: whether X and Y share a stream)*.
     - Independent θ per frequency? The manual is silent (§13).
     - With uniform phases on the full circle, `E[s sᴴ] = Φ_m Λ_m Φ_mᵀ`, so the coherency is recovered in the ensemble mean.
   - **Deterministic Linear (AS):** `s = Σ_k sqrt(λ_k) φ_k` (θ = 0).
     - Mode sign ambiguity: eigenvectors are defined only up to ±1, and the manual warns that results depend on this sign.
     - Recommended convention *(I)*: choose the sign of φ_k so that `Σ_i φ_ik ≥ 0`. Then the f → 0 case (Σ = 11ᵀ, λ₁ = N, φ₁ = 1/√N) gives `s = 1`, which reproduces coherent motion. This is consistent with the manual's check that the zero-frequency ATF in the input direction should be ≈ 1.00 (cross-ref line 2761).
     - The "ATF phase adjustment" variant is unspecified (§13).
   - **Quadratic / single mode** (`nmodes = −k`): `s = sqrt(λ_k) φ_k`.
4. (I) Wave passage: `s_i ← s_i · exp(−iω τ_i)`, with τ_i from §3.15.
5. (I) Multiple excitation: `s_i ← s_i · SAR_zone(i)(ω)`.
6. (M) Downstream use: ANALYS turns the FILE77 variations into the seismic load vector. Two options:
   - **FFL**: randomize the free-field *load* at the interaction nodes. EPRI-validated; slightly conservative.
   - **FFM**: randomize the free-field *motion*, then multiply by the coherent soil impedance. Not recommended for embedded structures.

   In both cases the coherent free-field motion of each interaction node comes from SITE (FILE1), so the incoherent free-field motion is `u_i,inc(ω) = s_i(ω) · u_i,coh(ω)` *(I for FFM; see ANALYS spec)*.

Ill-conditioning guard (recommended, inferred):
- Nodes at different depths but the same plan position give identical rows in Σ.
- Option: build Σ on **unique plan positions** (merge within a tolerance), decompose that matrix, then broadcast the modes back to all nodes. N in (6.1) is then the number of unique positions.
- Otherwise, report the condition and the trace check. The manual's alternative is the "per level" approach (§7).

---

## 5. HOUSE batch commands (input-file equivalents of the dialog) [IO]
| Command | Arguments (exact order) | Dialog fields |
|---|---|---|
| `HOUSE` | `<gravity>,<gelev>,<opmode>,<dim>,<imp>,<coh>,<wpass>,<me>,<cmplxspec>` | gravity, ground elev., op. mode, dimension (0/1/2), method (0 FV/1 FFV/2 FI), coherent(0)/incoherent(1), wave passage 0/1, multiple excitation 0/1, complex SAR 0/1 |
| `INCOH` | `<gammax>,<gammay>,<gammaz>,<alpha>,<ngp>,<ipr>,<nmodes>,<met>,<HSeed>,<VSeed>,<RandPhz>` | coherence X/Y/Z, alpha (or Vs), embedded layers / mesh points per level, print flag 0/1, modes (0 / k / −k), units 0 = ft / 1 = m, horizontal seed, vertical seed, random phase angle in degrees |
| `WPASS` | `<appv>,<ang>,<cohf>` | apparent velocity, angle of Line D, `<cohf>` "directional coherence factor" (§13) |
| `ME` | `<no>,<nfirst>,<nlast>,<xc>,<yc>,<zc>` | motion number (1–10), first/last foundation node, control point (not used) |
| `AMP` | `<no>,<a1>,…,<a100>` | SAR list for motion `<no>`; `<a1>=0` deletes it |
| `ANSYSMODELTYPE` | `<type>` | 1 = embedded, 2 = surface |
| `GRAVITY`, `GROUNDELEV` | value | gravity, ground elevation (UI commands, §9.7.16/17) |
| `INTGEN` | `<type>,[level skip]` | interaction-node generation (§3.3) |
| `ETYPE` | `<e1>,<e2>,[<inc>],<type>` | 0 default / 1 structure / 2 excavated or embedded shell |
| `EINT` | `<e1>,<e2>,[<inc>],<order>` | SOLID: 0 rectangular, 1 skewed, 2 extremely distorted (default 2×2). TSHELL: 0 Reduced (default), 1 Selective |

- PREP syntax is free format: `Keyword, p1, p2, …`, comma-delimited, keywords not case-sensitive (cross-ref §1.6).
- `WRITE` saves all dialog state into `modelname.pre`.
- `AFWRITE` writes `.hou` (and `.frc` for FORCE). `AFWRITE` runs `CHECK` first and does not write the affected file if there are errors.

---

## 6. Nonlinear SSI input file `.pin` [IO][ADV]

- Opened through the HOUSE dialog's "Input Data" button.
- Free format; the example uses commas.
- Defines the **initial** properties of the near-field nonlinear soil element groups.

Layout:
```
Line 1  (3 items): NGRP, ESF, NCURV
          NGRP  = number of nonlinear soil element groups
          ESF   = effective strain factor (e.g. 0.60)
          NCURV = number of soil material curves defined in SOIL (G-γ / D-γ pairs in FILE73)
Repeat for each nonlinear group g = 1..NGRP:
  Line 2 (4 items): IGRP, NMAT, NELEM, ISTR
          IGRP  = group number of the nonlinear soil group
          NMAT  = number of materials in the group (need not equal number of layers)
          NELEM = number of solid elements in the group
          ISTR  = effective shear strain flag:
                  0 = maximum component shear strain among X, Y, Z
                  1 = octahedral shear strain (3D)  /  maximum shear strain (2D)
  Line 3 (3 items) x NMAT lines: GFAC, DFAC, ICURVE
          GFAC   = initial shear-modulus reduction factor (1.00 = same G as free field)
          DFAC   = initial damping-ratio factor (1.00 = same damping as free field)
          ICURVE = soil material curve order number (pair of G and D curves in FILE73)
```
The manual's wording ("block after 1st line for all groups") means lines 2 and 3 repeat for each group.

Example (screenshot `C:\SSI\Demo1\XDIR\RBX.pin`): one group, ESF 0.60, 2 curves; group 2 has 5
materials and 180 elements, with octahedral strain (flag 1):
```
1, 0.60, 2
2, 5, 180, 1
1.0, 1.0, 1
1.0, 1.0, 1
1.0, 1.0, 1
1.0, 1.0, 1
1.0, 1.0, 1
```

Iteration workflow (context, cross-ref §1.5.4 and §4.2.4):
1. Run HOUSE. On the first run it reads `.pin` and writes `.liq`. On later runs, if `.liq` contains 1, it reads FILE74.
2. Run ANALYS as restart **Mode 2 "New Structure"**.
3. Run STRESS: it computes effective strains, takes the updated G and D from the FILE73 curves, and writes FILE74. FILE78 from HOUSE is used here.
4. For 3-directional input, run each of X, Y and Z, then **`COMBIN_XYZ_STRAIN`**, which SRSS-combines the three directional FILE74 strains into one FILE74.
5. Repeat until converged.

Notes:
- The same mechanism works for any nonlinear material modelled with solids (massive concrete, rubber, …). Add a new G-γ / D-γ pair to FILE73 (or define it in SOIL) and reference it with ICURVE.
- Equivalent-linear convention (standard SHAKE practice, inferred): `γ_eff = ESF · γ_max`, then interpolate G/G_max and D from the curve at γ_eff.

---

## 7. Multilevel ("per level") incoherent analysis and BuildFILE77 [ADV]
- Needed for deeply embedded foundations in soil, where coherency at foundation level differs
  from that at the surface (Model 6 is a surface model). Also needed when the coherency matrix of
  an embedded mesh is ill-conditioned.
- Procedure:
  1. Run HOUSE **incoherent** separately for each interaction-node level, defining only that level's interaction nodes. Each run gives one FILE77.
  2. Run a final HOUSE **coherent** run with all interaction nodes to build the correct FILE4 for ANALYS.
  3. Combine the per-level FILE77s with the auxiliary `BuildFILE77.exe` (default folder `C:\ACSV300\Build_FILE77`) into one FILE77 covering all interaction nodes.
- Accurate **only** if the eigenmodes keep the same node-numbering pattern at every level, so mode signs do not flip between levels. This needs mode-by-mode checking.
- The implementation needs a deterministic eigenvector sign convention (§4, step 3) and a merge tool `build_file77`.

---

## 8. HOUSE checks to implement (CHECK / data-check mode) [CORE][UI]
| Error | Condition |
|---|---|
| 1 | gravity ≤ 0 |
| 44 / 120 | frequency set undefined / empty |
| 49 / 50 | Δt < 0 / NFFT < 0 (NFFT also: power of 2, ≤ 32 768) |
| 58 | coherence parameter < 0.1 |
| 59 | Luco–Wong mean Vs ≤ 0 |
| 60 | mesh points per embedment level ≤ 0 |
| 113 | apparent velocity ≤ 0 |
| 114 | directional coherence factor < 0 |
| 115 | multiple excitation on, no data |
| 116 / 117 | illegal first / last node for motion i |
| 118 | SAR outside [0, 10] |
| 119 | SAR count ≠ number of frequencies |

Additional rule checks (warnings unless noted):
- R3 interface rule: **error**.
- R4 bottom-up interaction ordering: **error when incoherent**.
- R5, R6 layering rules.
- Incoherent analysis with 2D or SYMM: **error**.
- Types 2–7 or multiple excitation without wave passage: **error**.
- FI/FFV without an FV sensitivity study: warning.
- SHELL without FIXROT/FIXSHLROT: warning.
- `nmodes` < 0 with the Stochastic option: **error** *(inferred)*.
- Nsim > 50: **error**.

---

## 9. FORCE module [CORE][IO]

### 9.1 Purpose
- For **each specified frequency**, FORCE builds the load-vector entries for **external forces**
  acting directly on the structure: impact, rotating machinery, wave forces, moving loads, or unit
  forces used to compute flexible-foundation impedance.
- Results go to **FILE9**.
- Input file: `modelname.frc`, written by `AFWRITE`.
- FORCE is not used in seismic problems, except for computing foundation impedances.
- **Seismic and external-force analyses cannot be combined in one run.** Run them separately and superpose linearly (cross-ref §1.5.2).

### 9.2 FORCE dialog (screenshot, printed p.145)
Same tab strip as HOUSE; FORCE tab active. Buttons Ok/Cancel.
| Group / Control (exact label) | Widget | Screenshot value | Meaning / mapping |
|---|---|---|---|
| Operation Mode: `Solution` / `Data Check` | radio | Solution | `FORCE,<opmode>`: 0 = complete solution, 1 = data check only |
| `Acceleration of Gravity` | text | 32.2 | same as SITE |
| `Frequency Step` | text | 0 | Hz; same as SITE. 0 or blank means `Δf = 1/(Δt·NFFT)` is used |
| `Time Step of Motion Control` (text calls it "Time Step of Control Motion") | text | 0.005 | sec; same as SITE |
| `Nr. of Fourier Components` | text | 4096 | power of 2; same as SITE |
| `Frequency Set number` | text | 1 | FREQ set; same as SITE |

Frequency rules are the same as SITE:
- Time-history analysis: give Δt and NFFT; Δf may be left blank.
- Single harmonic analysis: Δf must be given; Δt and NFFT are unused.
- `f_j = n_j · Δf`.
- Errors 44, 48 (Δf < 0), 49, 50, 120 apply.

### 9.3 Load definition [CORE][IO]
All forces and moments defined in the UI are written to the FORCE input file:

| Command | Syntax | Meaning |
|---|---|---|
| `F` | `F,<n>,<fx>,<fy>,<fz>,<tx>,<ty>,<tz>` | force **factors** and **arrival times** for X, Y, Z at node n |
| `MM` | `MM,<n>,<fxx>,<fyy>,<fzz>,<txx>,<tyy>,<tzz>` | moment factors and arrival times about X, Y, Z at node n |
| `FDEL` / `MMDEL` | `<n1>,[<n2>],[<inc>]` | delete forces / moments |
| `FLIST` / `MMLIST` | `[<n1>],[<n2>],[<inc>]` | list forces / moments |
| `FSCALE` / `MSCALE` | `[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` | scale factors. A scale value of 0.0 is treated as 1.0. If n1/n2 are omitted, the last 2 defined entries are used. |
| `FREAD` / `MREAD` | (named in the manual; syntax not documented in §9.5) | read forces / moments from file (§13) |

Load model (manual):
- All nodal loads share **one reference time history** with the maximum reference amplitude and **zero arrival time t₀**: the reference load starts at t = 0.
- Each nodal load component is defined relative to the reference by:
  - a **load factor** = (maximum amplitude of that load) / (maximum reference amplitude),
  - an **arrival time** (time lag), which becomes a complex phase in frequency.
- This allows **moving loads**.
- Up to **500 load cases** per ANALYS batch run, RAM permitting.

Frequency-domain load vector *(inferred, standard, consistent with the manual's "time lag that corresponds to complex force phasing")*:
```
P_k(ω_j) = a_k · exp(−i ω_j t_k)        k = loaded DOF (Fx,Fy,Fz,Mx,My,Mz at node n)
```
- e^{+iωt} convention, ω_j = 2π f_j.
- FILE9 stores `P(ω_j)` for every SSI frequency and is used by ANALYS (in place of FILE1) to get transfer functions per unit reference load.
- The response time history is then `IFFT[ H(ω) · F_ref(ω) ]`, where F_ref is the FFT of the reference load history. Where F_ref is supplied (MOTION input) is outside this section.

### 9.4 Load-case files [IO]
- Multiple load cases:
  1. Run FORCE once per case.
  2. Copy each resulting `FILE9` to `FILE9001`, `FILE9002`, … `FILE9500`.
  3. ANALYS then runs them as "Simultaneous Cases" (restart Mode 3 "New Dynamic Loading" is also available).
- The manual also says "Only one digit load case number can be appended to FILE9". This contradicts the 3-digit names (§13).

### 9.5 FORCE checks
- **Error 61**: no forces or moments defined.
- **Error 62**: force or moment on an undefined node (define it with `N`, or delete with FDEL/MMDEL).
- Frequency errors as in 9.2.

---

## 10. Recommended practice values (summary) [UI help text]
| Item | Value |
|---|---|
| Luco–Wong coherence parameter | 0.10–0.30 (≥ 0.1 enforced) |
| Alpha, isotropic | 0.50 |
| Alpha, directional | 0.1 or 0.9 (3:1 distance weighting) |
| Apparent velocity to suppress wave passage | 1.0e9 |
| Number of incoherent modes, SS and AS | 0 (all) |
| Mode-truncation criterion (SRSS) | ≥ 90 % cumulative variance per frequency |
| Rigid-mat stick models, SRSS | ~10 modes typical (cross-ref line 2711) |
| Flexible mats | tens to hundreds of modes, so use SS instead |
| Stochastic simulations | random seeds, phase 180°, Nsim ≤ 50 (e.g. 20) |
| Node optimizer | ON for large embedded models |
| Effective strain factor (example) | 0.60 |
| FI/FFV use | always with an FV sensitivity study |
| Multiple-excitation zones | ≤ 5000; contiguous node ranges |

---

## 11. Verification tests suggested for this section (to make the implementation verifiable)
1. **Trace test:** for every frequency and direction, `|Σλ − N|/N < 1e-8` (eq. 6.1); Σ is symmetric with unit diagonal.
2. **Coherent limit:**
   - As f → 0 (or with V_app = 1e9 and coherency ≈ 1), Σ → 11ᵀ, so λ₁ = N and the others ≈ 0.
   - The AS synthesis gives s = 1, so the zero-frequency ATF in the input direction is ≈ 1.00 (manual check).
3. **Ensemble test (SS):** the mean of `s sᴴ` over Nsim samples converges to Σ (error ∝ 1/√Nsim). With seed 0 and phase 0, the result is the deterministic field.
4. **Mode count:** `nmodes = n` gives the cumulative printout equal to Σ_{j≤n} υ_j. `nmodes = −k` uses only mode k.
5. **Directional distance:** two nodes 1 apart along X' and Y': D = sqrt(2α) and sqrt(2(1−α)). For α = 0.1 the ratio is 1/3.
6. **Wave passage:** the phase difference between two nodes a distance L apart along Line D is ωL/V_app.
7. **SAR identity:** all SAR = 1 reproduces the uniform-input results. A complex SAR e^{iφ} on a zone rotates that zone's free-field phase by φ.
8. **Optimizer invariance:** results mapped back through `.map` match the un-optimized run to round-off. Bandwidth/profile does not increase.
9. **FORCE:** for a single load with factor a and arrival t₀, `IFFT(P·F_ref)` equals a·f_ref(t − t₀) (circular shift within the Fourier period).
10. Manual V&V cross-references: Problem 31 (optimizer), Problem 35 (multiple excitation), Problem 39 (user coherency 61 × 81).

---

## 12. Output naming summary (HOUSE and FORCE) [IO]
- `FILE4` / `modelname.N4`, `COOSK`, `COOSM`, `FILE77` (or `FILE77001…FILE77050`), `FILE78`, `.liq`, `.hounew`, `.map`
- `FILE9` (or `FILE9001…FILE9500` after renaming)
- Downstream: `FILE8` / `FILE8xxx`, `SRSSTF.txt` (MOTION), `FILE74` (STRESS), `FILE73` (SOIL)

---

## 13. Open questions / ambiguities (implementer must decide)
1. **Dialog defaults.** The screenshots show one example state (FI selected, Incoherent, Alpha = 1000, Model 3, wave passage on, Last Foundation Node 69). No defaults are documented. Proposed defaults: Solution, 3D, FV, Coherent, γx = γy = 0.1, γz = 0.2, α = 0.5, nmodes 0, Linear, Deterministic, WP off, V_app = 1e9, θ = 0, Model 1, ME off.
2. **Alpha = 1000 with Model 3** in the screenshot. This is meaningful only as a Luco–Wong Vs. Treat it as a leftover value. Validate α ∈ [0, 1] for Types 2–7 and α > 0 for Type 1.
3. **α = 0.1 direction contradiction** (§3.10). The formula says Y distances are weighted 3× (faster decay in Y); the dialog description says faster decay in X. Use the formula.
4. **Rotation of the distance frame.** Are DX, DY measured in the Line-D frame always, only for Types 2–7, or only when wave passage is on? Proposed: Line-D frame for Types 2–7; global horizontal distance for Type 1.
5. **Which field carries the model type and which carries α.** `WPASS <cohf>` is described as a "directional coherence factor" (Error 114: < 0), while `INCOH <alpha>` is the directionality parameter. The dialog places "Unlagged Coherency Model" in the Wave Passage group, which suggests `<cohf>` = model type. Proposed: `WPASS <cohf>` = Unlagged Coherency Model (1–7), `INCOH <alpha>` = α/Vs. This needs a decision.
6. **Superposition Mode and the stochastic Number of Simulations have no documented command arguments.**
   - Proposed: extend `INCOH` with trailing optional arguments `<supmode>` (0 = Linear, 1 = Quadratic) and `<nsim>`.
   - Alternatively, infer stochastic from non-zero seeds and Quadratic from negative nmodes.
7. **AS "with/without ATF phase adjustment".** The switch and the algorithm are unspecified. Proposed: phase adjustment = deterministic eigenvector sign normalisation (Σ_i φ_ik ≥ 0) per mode and frequency. The "without" variant leaves the raw eigensolver signs. Confirm against EPRI TR-1015111.
8. **Random phase details.** Are phases independent per frequency, per direction and per sample? Do X and Y share HSeed? Proposed: independent per frequency, mode, direction and sample; HSeed seeds separate streams for X and Y; VSeed for Z. Use a documented PRNG (e.g. PCG64) for reproducibility.
9. **"Number of Embedded Layers" vs INCOH `<ngp>` ("number of mesh points per each embedment level")**, Error 60 (≤ 0 illegal) vs screenshot value 0. The algorithmic role is unknown. It may group interaction nodes per level for a per-level coherency build. Proposed: 0 = surface/auto; the computation uses unique horizontal projections (§4 guard); store the value for fidelity.
10. **Mode ordering.** Eq. 6.4 says modes are "ordered with increasing contributions". Proposed: descending eigenvalue, so the cumulative sum increases fastest; "first n modes" are the n largest.
11. **Units flag `met`** (ft vs m) is not in the dialog. Proposed: derive it from gravity (≈ 32.2 → ft; ≈ 9.81 → m), with an explicit override. Convert distances to metres before evaluating Abrahamson models.
12. **Abrahamson model equations (Types 2–6) are not given in the manual.**
    - The implementer must source them: Abrahamson 1993 (US Seminar on Seismic Evaluation…), EPRI 2005 report (all sites, surface), EPRI 2006 report (embedded), 2007 hard-rock (Pinyon Flat, EPRI TR-1015110), 2007 soil sites; see manual references 1–4.
    - The commonly published 2005/2007 plane-wave form is
      `γ_pw(f,ξ) = [1 + (f·tanh(a3 ξ)/(a1·fc(ξ)))^{n1}]^{-1/2} · [1 + (f·tanh(a3 ξ)/a2)^{n2}]^{-1/2}`
      with separate horizontal and vertical coefficient sets. Coefficients must be transcribed from the original reports, or from USNRC DC/COL-ISG-01 where reproduced. **Do not use from memory.**
    - Recommended architecture: a pluggable coherency-function registry with coefficient tables in data files and unit tests against published plots. Type 7 lets any tabulated function be checked independently.
13. **Luco–Wong form.** Assumed `exp[−(γ ω D/Vs)²]` with ω in rad/s (per Luco & Wong 1986). Confirm whether ACS SASSI uses ω or f.
14. **Wave-passage sign and reference point.** Is τ = 0 at the global origin, the first interaction node, or the foundation centroid? A constant shift does not change relative motion, but it does change absolute phases in FILE8. Proposed: global origin, with the e^{+iωt} convention.
15. **User coherency file format.** Row/column orientation (frequency rows × distance columns assumed), delimiter, header lines, how the non-default sizes (e.g. 61 × 81) are encoded in `.hou`, and interpolation and extrapolation rules.
16. **Complex SAR input syntax** (`Use Complex Spectral Amp.`) is undocumented. Proposed: pairs `re im` per frequency, or `(re,im)` tokens. Also: does the [0, 10] check apply to the modulus? Proposed: yes.
17. **Multiple excitation limits.** The text allows up to 5000 zones, while `ME <no>` is limited to 1–10. Also undefined: SAR for interaction nodes outside every zone (proposed 1.0); how ME combines with incoherence (proposed multiplicative, after synthesis); what the X/Y/Z control points are for (unused).
18. **`Non-Uniform Motion` / `Non-Uniform Soil` checkboxes** are disabled and undocumented. Implement as placeholders.
19. **"Input Data" button.** Pile interface (PINT, not available) or `.pin` nonlinear input? Proposed: `.pin` editor.
20. **"Below ground surface" test for ETYPE 0** (centroid vs all nodes) and the tolerance.
21. **Node-renumbering algorithm** is unspecified. Proposed: reverse Cuthill–McKee or Sloan on the structural plus excavation graph, constrained so interaction nodes stay bottom-up ascending. Write `.hounew` and `.map` (old new pairs, one per line).
22. **`.pin` material ordering.** How do the NMAT lines map to the group's materials: ascending material index, order of first appearance, or layer order? Proposed: ascending material/layer index used in the group.
23. **`.liq` semantics** beyond "contains 1 means use FILE74"; the format of FILE74 and FILE78 (STRESS spec).
24. **FILE77, FILE4, COOSK, COOSM, FILE9 binary layouts** are proprietary and not documented. Define our own versioned formats (e.g. NumPy `.npz` or HDF5 with the legacy names as filenames), keeping the legacy names.
25. **FORCE: where is the reference load time history supplied** (MOTION input?), and what exactly does FILE9 contain (normalized vectors per frequency assumed)? `FREAD`/`MREAD` syntax is not documented. Load-case naming: "only one digit" vs `FILE9001…FILE9500`; proposed 3-digit suffix.
26. **FE table footnotes** `*` (SHELL) and `**` (TSHELL) have no footnote text.
27. **Typos to normalise:** ".hownew" → `.hounew`; "STRSSTF.txt" → `SRSSTF.txt`; "equation 5.1/5.4" → 6.1/6.4; dialog labels "Sesmic", "Inupe", "Amplifcation" should be corrected in our UI. Keep the exact field names otherwise.
28. **Deterministic label:** the dialog says "(Median)", the text says "(Mean)". Use "Deterministic (Median) Incoherency Input" to match the dialog.
