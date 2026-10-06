# 00 — Consolidated Requirements: Educational Re-implementation of ACS SASSI (SASSI-EDU)

**Status:** Normative master requirements document. It drives the implementation.
**Date:** 2026-10-01.
**Scope:** an educational Python re-implementation of the ACS SASSI Version 3 (IKTR9) methodology. It keeps
the program's modules, file names, command language, option codes and UI nomenclature, and it must be
technically accurate and verifiable against published solutions.

**How this document relates to the other specs.** Specs 01–11, R1 and R2 under `docs/spec/` record what the
manual says (01–11), the underlying SASSI theory (R1) and the reference solutions (R2). This document
consolidates them, sets priority tiers and resolves every conflict and open question that affects
implementation (§7). **When this document and a numbered spec disagree, this document wins.** Where this
document says "see spec NN §x", that section is normative for the detail it covers.

| Spec | Content | Role here |
|---|---|---|
| `01_capabilities.md` | Ch. 1: capabilities, limits, module chain, element library, modelling rules, nonlinear frameworks, interpolation names, units, input, CHECK | architecture frame |
| `02_theory_modules.md` | Ch. 2 theory (Eq. 2.1–2.3), Ch. 3 modules, Tables 3.1/3.2, restarts | file flow, naming |
| `03_guidelines.md` | Ch. 4: 13-step procedure, 20 engineering considerations, run types, validators G-01…G-29 | orchestration, validators |
| `04_model_ui_modules.md` | Ch. 5, 6.1–6.4: workflow, models in memory, macros, variables, MERGE/MERGESOIL, cuts, menus, converters, module runs | interpreter, UI, IO |
| `05a`–`05d` | Ch. 6.5–6.7: every Options dialog field and its command mapping; module algorithms | module options |
| `06_plotting_toolbar.md` | Ch. 7–8: plots, line maths, animations, toolbars | plotting |
| `07`–`11` | Ch. 9–11: every command, CHECK catalogue (Errors 1–128, Warnings 1–11), references | command language |
| `R1_sassi_theory.md` | SASSI/SHAKE formulations (TLM, POINT core, Eq. 2.1, Tajirian interpolation, damping) with numerical checks in `R1_checks/` | **normative theory** |
| `R2_benchmarks.md` | Reference solutions with tolerances; data in `R2_benchmark_data/` | verification |

**Contents.** §0 Conventions (tiers, provenance, global conventions, environment) · §1 Glossary · §2 Module
data flow · §3 Command language (rules, abbreviations, inventory by category with tiers, extension commands) ·
§4 Computational methodology per module (with the conflict register §4.16) · §5 User-interface requirements ·
§6 Verification plan · §7 Decisions · §8 Implementation order.

---

## 0. Conventions used in this document

### 0.1 Priority tiers

| Tier | Meaning | Content (summary) |
|---|---|---|
| **P0** | Core linear SSI chain. Must work end-to-end and be verified before anything else ships. | Model definition (nodes, elements, materials, layers, properties, masses, BCs, interaction nodes by INT); EQUAKE (frequency-domain matching, SRP checks); SOIL (SHAKE equivalent-linear); SITE; POINT3; HOUSE (3D, coherent, FV/FI/FFV through the interaction set); FORCE; ANALYS (incl. restarts, simultaneous cases, global impedance); COMBIN; MOTION; RELDISP; STRESS; AFWRITE / CHECK / RUNxxx; `.pre` I/O (INP, WRITE, SAVE/RESUME); model checks (EXCSTRCHK, FIXEDINT, HINGED…); FIXROT family |
| **P1** | Important engineering features | ANSYS `.cdb` import and APDL export; interaction-node and excavation generation (INTGEN, ETYPEGEN, EXCAV, SOILMESH, MERGE, MERGESOIL, WELD, NCOM…); macros and variables; cuts, submodels and section-cut calculations; line maths (envelope, broadening, SRSS, combinations); plotting commands and animations; incoherency, wave passage, multiple excitation; 2D SSI (POINT2, PLANE); symmetry (SYMM); nonlinear-soil SSI iterations (`.pin`, FILE74/78, COMB_XYZ_STRAIN); EQUAKE time-domain (AB) refinement, X–Y correlation and multiple random-seed trials; frequency-refinement tools (CRITFREQ, FILE8 frequency removal); AFWRBAT |
| **P2** | Advanced or cosmetic | Option NON (NONLINEAR, panels, BBC, springs, SHEAR/BBCGEN, COMB_XYZ_THD); SOIL-NON (NLSOIL/NLSLAYER); water modelling (FILLPOOL…); binary databases (BINOUT products, COMB*DB, *DBANI…); TSHELL element and thick-shell commands; Option A (LOADGEN), Option AA (SSI2ANSYS, HOUSEFSA/ANALYSFSA, ANSYSMODELTYPE); cosmetic UI (COLOR, SHADEROPTIONS, STIPPLE, DEBUG, fonts) |
| **OOS** | Out of scope | Option PRO (separate manual); GT-STRUDL converter (not included in V3); LIQUEF, PINT (not in V3); New Load Vector (ANALYS Mode 6, vendor-only); legacy byte-exact PREP/SASSI2000 deck formats |

A command or feature in a lower tier must still be **parsed and stored** in P0 (so `.pre` files written by the
original UI load without loss) and must print `<CMD> is not available in this build (tier Pn)` when executed.

### 0.2 Provenance markers

| Marker | Meaning |
|---|---|
| **[M]** | Stated in the ACS SASSI V3 manual (via specs 01–11) |
| **[S]** | Confirmed by an external source cited in R1/R2 |
| **[R]** | Reconstruction or engineering decision (recorded in §7 with an ID) |
| **[V]** | Checked numerically (R1_checks or R2) |

### 0.3 Global technical conventions (normative; details in §4.0 and §7)

| Item | Convention | Decision |
|---|---|---|
| Time factor | `e^{+iωt}`; outgoing waves `e^{-ikx}`, `H^{(2)}` Hankel functions | D-CNV-01 |
| FFT | forward `A_k = Σ a_n e^{-2πikn/N}`, inverse with `1/N` (`numpy.fft.rfft/irfft`) | D-CNV-02 |
| Complex modulus | `G* = G(1 − 2β² + 2iβ√(1−β²))` (SASSI/SHAKE91), applied to G with βs and to M = λ+2G with βp; `(1+2iβ)` selectable for benchmarks only | D-CNV-03 |
| Mass matrices | ½ lumped + ½ consistent (SOLID, PLANE, SITE/POINT layers, excavated soil); BEAMS consistent; SHELL/TSHELL lumped | D-CNV-05 |
| Transfer functions | FILE8 holds `H = U/U_cp`, dimensionless (equal for displacement and acceleration) for seismic; DTF per unit load factor for vibration | D-CNV-06 |
| Units | any consistent set; gravity g sets weight→mass; British (g≈32.2 ft/s²) or SI (g≈9.81 m/s²) detected from g | D-CNV-08 |
| Geometry | global right-handed Cartesian, **Z up**; 2D models in the X–Z plane (Y ignored) | [M] |
| Frequencies | integer frequency numbers × Δf; `Δf = 1/(Δt·NFFT)` unless `fstep > 0` | [M] |
| Numbering | node ids user-defined; element numbers start at 1 per group, gap-free when written | [M] |

### 0.4 Implementation environment (as found in the project venv)

Python 3.9.6, numpy 2.0.2, scipy 1.13.1, matplotlib 3.9.4, plotly 7.1, pytest 8.4, tkinter (Tk 8.5).
h5py, VTK/pyvista and Qt are **not** installed. Requirements below therefore use NumPy `.npz` containers,
tkinter + matplotlib for the GUI, matplotlib/plotly for plots, and pytest for verification (D-GEN-02…05).

---

## 1. Glossary of ACS SASSI nomenclature

One-line meanings. "Tier" gives the tier of the feature the term belongs to. Section references point to
the normative spec.

### 1.1 Products, builds, options and auxiliary programs

| Term | Meaning | Tier |
|---|---|---|
| ACS SASSI | Commercial SSI program (Ghiocel Predictive Technologies, Version 3, build IKTR9) whose methodology and nomenclature are reproduced | — |
| SASSI | *System for Analysis of Soil-Structure Interaction* (Lysmer et al. 1981, UC Berkeley); SASSI2000 is the 1999 Berkeley update | — |
| SASSI-EDU | Name of this educational re-implementation (Python package `sassi`); not affiliated with or endorsed by the vendor | — |
| PREP | Legacy pre-processor (≤ IKTR4). Its `.pre` command language is a subset of the UI language | P0 |
| SUBMODELER, MAIN | Legacy submodeling and main programs, merged into the UI | — |
| UI | "ACS SASSI User Interface": command interpreter, model database, CHECK/AFWRITE, module launcher, plots | P0 |
| IKTR9 / IKTR9_650K (IKTR8_650K) / IKTR8_300K | Builds limited to 99,999 nodes / 650,000 nodes and 2.5 M equations / 300,000 nodes (extended integer fields) | info |
| Fast solver | HOUSEFS / ANALYSFS executables; one implementation here, names kept as aliases | P0 |
| HOUSEFSA, ANALYSFSA (also "ANALYSFA") | Option AA versions of HOUSE / ANALYS that use ANSYS matrices | P2 |
| NQA | Nuclear-QA version (10 CFR 50 App. B, ASME NQA-1); 56 V&V problems (not available) | info |
| Option A | Two-step SSI: SASSI global SSI, then ANSYS stress analysis with SSI motions as BCs (module LOADGEN); implemented (D-W3-11, `docs/user/OPTION_A.md`) | P2 |
| Option AA | "Advanced ANSYS": SSI with ANSYS dynamic matrices directly; topology-only `.hou` (module SSI2ANSYS); not implemented (ANSYSMODELTYPE is stored only) | P2 |
| Option NON | Nonlinear-structure SSI by iterative equivalent linearisation (NONLINEAR, COMB_XYZ_THD); implemented (D-W3-09, `docs/user/OPTION_NON.md`) | P2 |
| Option PRO | Probabilistic site response / SSI by Latin Hypercube Sampling | OOS |
| LOADGEN | Option A module writing ANSYS load / BC (APDL) files. The manual has two dialogs and no command; SASSI-EDU keeps the dialog fields as option records (LOADGEN, LOADGENDYN, LGFILE, LGNODE, LGTIME, LGMAP, LGOPT; LGLIST lists them) and RUNLOADGEN,[STATIC\|DYNAMIC] writes the deck `<model>.lgn` and runs the module (D-W3-11, D-LGN-01) | P2 |
| SSI2ANSYS ("SASSIANSYS Module") | Option AA super-element utility: MATRIX50 → GENERAL elements, or add SE matrices to ANSYS matrices; not implemented | P2 |
| NONLINEAR ("PANEL module") | Option NON module: hysteretic models, equivalent-linear properties, new `.hou`; implemented as a module (deck `.eql` written by AFWRITE, run by RUNNONLINEAR; iterations by NONLINITER; D-W3-09) | P2 |
| COMB_XYZ_THD (aliases COMBIN_XYZ_THD, COMBINE_XYZ_THD, COMB_XYZ) | Combines co-directional X/Y/Z relative-displacement histories for Option NON (algebraic sum per time step, D-NON-08); implemented as the command COMBXYZTHD and as `python -m sassi.modules.nonlinear COMB_XYZ_THD <inpfile>` | P2 |
| COMB_XYZ_STRAIN (alias COMBIN_XYZ_STRAIN) | SRSS of the X, Y, Z directional FILE74 effective strains; command COMBXYZSTRAIN (D-NLS-04) | P1 |
| Build_FILE77 (BuildFILE77) | Assembles per-level FILE77 incoherency files into one; command BUILDFILE77 | P1 |
| Remove_Frequencies_from_FILE8 | Deletes chosen (spurious) frequencies from FILE8; command REMOVEFREQ | P1 |
| LIQUEF, PINT | Modules not in V3 (greyed; the `<DEP>` slots of AOPT) | OOS |
| SOIL-EQL / SOIL-NON | SOIL sub-options: SHAKE equivalent-linear / DEEPSOIL-like time-domain hyperbolic (SOIL-NON implemented: `NLSOIL,1` with NLSLAYER sets in `<model>.nls`, D-W3-10) | P0 / P2 |
| Demo 1–12, V&V Problem n | Vendor examples cited by the manual (not available; used only for traceability) | info |

### 1.2 Modules

| Module | One-line meaning | Input deck | Tier |
|---|---|---|---|
| EQUAKE | Spectrum-compatible acceleration histories (LW frequency-domain + AB time-domain matching, baseline correction, SRP 3.7.1 checks, PSD/FFT/Arias) | `.equ` | P0 (AB, correlation, multiple seed trials P1) |
| SOIL | 1D free-field nonlinear site response: SHAKE equivalent-linear (SOIL-EQL) or hyperbolic time domain (SOIL-NON) | `.soi` (SOIL-NON: + `.nls`) | P0 / P2 |
| SITE | Linearised layered free field: Mode 1 Rayleigh/Love eigenproblems (FILE2); Mode 2 free-field motion for unit control motion (FILE1) | `.sit` | P0 |
| POINT (POINT2 2D, POINT3 3D) | Point-load solutions of the layered soil column → frequency-dependent flexibility data (FILE3) | `.poi` | POINT3 P0, POINT2 P1 |
| HOUSE (HOUSEFS, HOUSEFSA) | Element K and M of structure and excavated soil (FILE4/COOSK/COOSM), incoherency decomposition (FILE77), nonlinear-soil updates | `.hou` | P0 (incoherency P1) |
| FORCE | External-load vectors per frequency (FILE9) | `.frc` | P0 |
| ANALYS (ANALYSFS, ANALYSFSA) | Per frequency: impedance X_ff = F_ff⁻¹, assemble Eq. 2.1, load vector, solve → FILE8; restarts; global impedance | `.anl` | P0 |
| COMBIN | Merges two FILE8 solutions (FILE81 + FILE82 → FILE8) | (none) | P0 |
| MOTION | Interpolates nodal TFs (options 0–6), convolves with the control motion, outputs histories, ISRS, TFU/TFI, frames | `.mot` | P0 |
| RELDISP | Relative-displacement TFs and histories from complex `.TFI` files | `.rdi` | P0 |
| STRESS | Element stress/strain/force TFs, interpolation, convolution, maxima, histories, frames, FILE74 | `.str` | P0 |
| NONLINEAR | Option NON equivalent-linearisation module (implemented, D-W3-09) | `.eql` | P2 |
| LOADGEN | Option A: ANSYS loads and boundary conditions (APDL) from the SSI results (implemented, D-W3-11) | `.lgn` (written by RUNLOADGEN) | P2 |

### 1.3 SSI concepts, methods and analysis terms

| Term | Meaning |
|---|---|
| Flexible volume substructuring | SASSI method: SSI system = free field + structure − excavated soil, coupled through the free-field impedance at interaction nodes (Eq. 2.1) |
| Free field / free-field system | Horizontally layered soil on a half-space or rigid base, without structure or excavation |
| Structural system | FE model of structure + basement + optional near-field soil |
| Excavated soil (excavation volume) | FE model of the soil removed by the embedded structure, with free-field layer properties (SOLID 3D / PLANE 2D, ETYPE 2) |
| Near-field soil / backfill | Soil included in the structural FE model (ETYPE 1), e.g. backfill or irregular zones |
| Far-field soil | The layered free field outside the FE model, represented by the impedance |
| Basement / superstructure | Structure part in contact with excavated soil / part above or not connected to it |
| DOF sets s, i, w, f | superstructure; shared structure/excavation (basement interface); excavation-only; interaction DOFs |
| Interaction node | Excavated-soil node where X_ff acts and U′_f is applied; translational DOFs only; must lie on a layer interface at or below grade |
| Intermediate / interface / internal node | Legacy INT codes 1/2/3; stored and listed only, no computational effect |
| FV (Flexible Volume, "Direct") | Interaction set = all excavated-soil nodes; reference method |
| FI (Flexible Interface) | Interaction set = a subset on the excavation surface (FI-FSIN or FI-EVBN) |
| FI-FSIN / SM (Subtraction Method) | Interaction set = foundation–soil interface nodes (lateral + bottom); fastest; spurious peaks near excavated-volume frequencies |
| FI-EVBN / MSM (Modified Subtraction; ESM in some texts) | FSIN + ground-surface (top) face nodes of the excavation |
| FFV (Fast Flexible Volume) | EVBN + selected internal horizontal node levels ("level skip") |
| Surface foundation | No excavation; interaction nodes at ground level only; FV = FI = FFV |
| DES / SMR | Deeply embedded structure / small modular reactor (motivation for FFV) |
| Impedance X_ff (K + iωD) | Dynamic stiffness of the free field at the interaction DOFs; `X_ff = F_ff⁻¹` (Eq. 4.3) |
| Flexibility / compliance F_ff (f + ig) | Free-field displacement at interaction DOFs due to unit harmonic point loads |
| Dynamic stiffness C(ω) | `C = K* − ω²M` with complex K* (Eq. 2.2) |
| ATF / DTF / STF | Acceleration (total, per unit control motion) / displacement (per unit load) / stress transfer function |
| TFU / TFI | Computed (uninterpolated) TF at SSI frequencies / interpolated TF on the Fourier grid |
| U′_f | Free-field motion amplitudes at the interaction nodes for unit control motion (from FILE1) |
| Control motion / control point / control point layer / control direction | The input acceleration history; specified at the top of the control layer (an interface) in direction X′, Y′ or Z′ |
| x′y′z′ | SITE wave axes: z′ up, x′ in the vertical propagation plane, y′ = z′ × x′ |
| Coordinate transformation angle | ANALYS angle (deg) from x′ to global x, counter-clockwise about z |
| Transmitting boundary | Consistent (Waas–Lysmer–Kausel) boundary of a layered semi-infinite region, from SITE Mode-1 eigen-solutions |
| Thin-layer method (TLM) | Layered-soil discretisation with displacements linear through each sublayer (Lysmer, Waas, Kausel) |
| Half-space layer / generated layers | L layer giving the half-space properties / `n` buffer sublayers (0 = rigid base, 4–20, recommended 10–20) of total depth 1.5 λ_s(f) |
| Variable-depth method / viscous boundary | Frequency-dependent buffer depth / Lysmer–Kuhlemeyer dashpots (ρVs, ρVp) at the base |
| Central zone / R0 | Axisymmetric (plane-strain in 2D) FE core around the point load; "Radius of Central Zone" (0.90h square mesh, 0.85h triangular, h in 2D) |
| Frequency step Δf / frequency number / frequency set | `Δf = 1/(Δt·NFFT)`; SSI frequencies are integers n with `f = nΔf`; a numbered list (FREQ) of frequency numbers |
| NFFT / Fourier period / quiet zone | FFT length (power of 2) / `T = NFFT·Δt` / trailing zeros allowing free decay |
| Cut-off frequency f_cut | Highest SSI frequency; TF is zero above it |
| Passing frequency / 1/5-wavelength rule | `f_pass = Vs/(5h)`; layer and vertical element sizes `h ≤ Vs/(5 f_cut)` |
| Interpolation options 0–6 | Complex TF interpolation schemes (§1.6, §4.8) |
| Smoothing parameter S | Damps spurious interpolated peaks/valleys (options 0–5); 0 for coherent runs |
| Phase adjustment | Reduces the differential Fourier phase toward zero ("minimum-delay" upper bound) for incoherent runs |
| Initiation / restart | ANALYS full run (Mode 1) / reuse of saved impedance or factorisation |
| New Structure (Mode 2) | Restart with new FILE4, same interaction nodes and soil; reuses X_ff (COOX) |
| New Seismic Environment (Mode 3) | Restart with new FILE1 and/or FILE77; reuses X_ff and the factorised system (COOX + COOTK) |
| New Dynamic Loading (Mode 3) | Same as above with new FILE9 |
| New Time History | Only the control-motion history changes: re-run MOTION/STRESS/RELDISP only |
| New Load Vector (Mode 6) | Batch-only vendor restart; out of scope (D-ANL-09) |
| Simultaneous cases | Several right-hand sides per factorisation: X/Y/Z coherent, ≤ 50 incoherent simulations, ≤ 500 load cases |
| Global (unconstrained) impedance | `K_G = Tᵀ X_ff T` (6×6) about the control point; FOUNSTIF/FOUNDASH/FOUNDAMP/FOUNIMPD |
| FFL / FFM | Incoherency applied to the free-field load `X_ff U′_f` (default) / to the free-field motion U′_f (surface only) |
| Incoherency / coherency | Spatial variation of free-field motion / its correlation function γ(f, D) between points |
| Unlagged / lagged coherency | Without / with wave-passage phase `exp(−iω(τ_i−τ_j))` |
| Wave passage / Line D / apparent velocity | Horizontal propagation delay along a horizontal Line D at velocity V_app |
| Directionality factor α / coherence parameter γ | Anisotropy weight in `D = √(2(αΔX² + (1−α)ΔY²))` / Luco–Wong parameter (≥ 0.1) |
| Incoherent spatial mode | Eigenvector of the coherency matrix Σ = ΦΛΦᵀ |
| SS ("simulation mean") | Stochastic simulation with random modal phases; mean of ≤ 50 samples |
| AS | Algebraic sum of scaled incoherent modes (deterministic, Linear superposition) |
| SRSS TF / SRSS FRS | SRSS of modal ATF amplitudes (Quadratic) / of modal end results (ISRS, ZPA…) |
| Deterministic (median) input | Zero modal phases (HSeed = VSeed = RandPhz = 0) |
| Per-level approach | Coherency decomposition per embedment level, merged by Build_FILE77 |
| Multiple excitation / foundation zone / SAR | Non-uniform input: interaction-node ranges driven by complex spectral amplification ratios |
| SSSI | Structure–soil–structure interaction (several structures in one model) |
| Symmetry / antisymmetry plane | SYMM half/quarter models (not allowed with incoherency) |
| ISRS / ARS / FRS / RS | In-structure / acceleration / floor response spectrum |
| ZPA | Zero-period acceleration (peak absolute acceleration) |
| PSD | One-sided power spectral density of the strong-motion part, ±20 % band averaged |
| Arias intensity / strong-motion duration | `I_A(t) ∝ ∫a²dt` / time between 5 % and 75 % of I_A |
| Seed record | Recorded history whose Fourier phasing is kept by EQUAKE |
| LW / AB | Levy–Wilkinson frequency-domain / Abrahamson (RspMatch-type) time-domain spectral matching |
| Baseline correction | Removal of polynomial drift (EQUAKE: FLUSH complex-frequency + time-domain polynomials; MOTION: Hudson–Housner) |
| Outcrop / within motion | `2A_m` (twice the up-going wave) / `A_m + B_m` (in-column) at a layer top |
| Equivalent-linear method | Seed–Idriss/SHAKE iteration of strain-compatible G and β |
| Effective strain / strain ratio R_γ / ESF | `γ_eff = R_γ·γ_max`; SOIL `<ratio>` (0.60–0.70); `.pin` effective-strain factor |
| Dynamic soil property | Named pair of curves G/Gmax(γ) and D(γ) (DYNP), strain and damping in percent |
| Strain-compatible properties | Iterated Vs, β from SOIL (FILE88) |
| Primary / secondary soil nonlinearity | Free-field (SOIL) / local SSI-induced (near-field iterations) |
| Relative displacement / free-field reference | Displacement relative to a reference node DOF / synthetic unit-amplitude zero-phase reference TFI |
| Panel (macroshell) | Option NON wall panel: coplanar vertical SHELL group deforming in shear or bending |
| BBC / cracking point / yield num | Backbone curve (origin omitted; point 1 = cracking) / 1-based yield-point index |
| CMS / CMB / TAK / GMR | Cheng–Mertz Shear / Cheng–Mertz Bending (not in V3) / Takeda / General Masing Rule hysteresis (BBC types 1/2/3/4) |
| EDF ("Disp. Factor") | Equivalent-linear displacement factor, `x_eq = EDF·max\|x\|` (≈ 0.8) |
| Ductility / Fμ | `max\|x\|/x_cracking` / initial elastic force ÷ effective nonlinear force (inelastic absorption factor) |
| Gravity trick | One-period 1 g harmonic over the Fourier period to obtain static gravity effects in a frequency-domain code |

### 1.4 Wave types (SITE)

| Code (WAVE `<type>`) | Wave | Family (`SITE <wopt>`) | Motion | `<opt>` values |
|---|---|---|---|---|
| 1 | R (Rayleigh) | 0 (R/SV/P, in-plane x′z′, 2 DOF/node) | elliptical, x′z′ | 0 none, 1 shortest wavelength, 2 least decay |
| 2 | SV | 0 | x′z′, perpendicular to propagation | 0 none, 1 field |
| 3 | P | 0 | along propagation | 0 none, 1 field |
| 4 | SH | 1 (SH/L, anti-plane y′, 1 DOF/node) | y′ | 0 none, 1 field |
| 5 | L (Love) | 1 | y′ | 0 none, 1 field |

Standard nuclear practice [M]: vertically propagating SV for X, SH for Y, P for Z (incident angle 0).
Wave ratios are given at Frequency 1 and Frequency 2, interpolated linearly, each in (0, 1], summing to 1.

### 1.5 Elements and element data

| GROUP type | Name | Nodes in `E` | DOF/node | Formulation | Mass | Tier |
|---|---|---|---|---|---|---|
| 1 | SOLID | 8 (prism/pyramid by repeated nodes 7=8; 5=6,7=8; 5=6=7=8) | 3 | trilinear hexahedron; optional 9 incompatible modes for structural solids | ½ lumped + ½ consistent | P0 |
| 2 | BEAMS | I, J, K (K = orientation node, no DOFs) | 6 | 3D Timoshenko frame (shear areas; 0 = no shear deformation); end releases | consistent | P0 |
| 3 | SHELL | I, J, K, L counter-clockwise (L omitted = triangle) | 6 | flat thin shell, Kirchhoff plate + membrane, **no drilling stiffness** | lumped | P0 |
| 4 | PLANE | I, J, K, L (X–Z plane) | 2 (UX, UZ) | plane strain, unit thickness | ½ + ½ | P1 |
| 5 | TSHELL | I, J, K, L | 6 | thick shell, Mindlin–Reissner; automatic small drilling stiffness; EINT 0 reduced / 1 selective | lumped | P2 |
| 7 | SPRING | I, J | 6 | six uncoupled global spring constants + one damping ratio | none | P0 |
| 9 | GENERAL (GM) | I, J (global input) or I, J, K (local input) | 6 | user 12×12 complex stiffness (MXR + i·MXI) and mass (MXM) | user | P0 |

| Term | Meaning |
|---|---|
| ETYPE | Per-element flag: 0 implicit (SOLID/PLANE excavated if below grade, else structure), 1 structure, 2 excavated soil (SOLID/PLANE) or buried shell (SHELL/TSHELL) |
| ETYPEGEN | Sets ETYPE for all elements (0 by location, 1 all structure, 2 all excavation) |
| EINT | Integration order: SOLID 0 rectangular / 1 skewed / 2 extremely distorted; TSHELL 0 reduced / 1 selective |
| MSET / MACT | Assign / set active material index (soil-layer L index for excavated SOLID/PLANE) |
| RSET / RACT | Assign / set active property index (R for BEAMS, SC for SPRING, MX for GENERAL) |
| KI / KJ | Beam end-release digits for P1, P2, P3, M1, M2, M3 at node I / J (1 = released) |
| THICK | Shell thickness per element |
| M material | type 1 (E, ν), 2 (constrained M, G), 3 (Vp, Vs); plus specific weight, βp, βs |
| L soil layer | thickness, specific weight, Vp, Vs, βp, βs; used by SITE/POINT/SOIL and excavated elements |
| R real property | A, As2, As3, J, I2, I3 (beam local axes) |
| SC spring property | kx, ky, kz, kxx, kyy, kzz (global) + damping ratio |
| MXR / MXI / MXM | Upper-triangle rows of the GENERAL element real stiffness / imaginary stiffness / mass (or weight) |
| MT / MR / MUNITS | Nodal translational / rotational masses; per-node units flag 0 mass, 1 weight |
| F / MM | Nodal force / moment factors with arrival times on one reference load history |
| D | Fixity codes 0 free / 1 fixed for UX…ROTZ (labels DISP, ROT, ALL) |
| Local axes | BEAMS: e1 = I→J, e2 toward K in plane IJK, e3 = e1×e2. SHELL: x′ through edge midpoints LI→JK, z′ normal, y′ = z′×x′ |

### 1.6 Option-code quick reference

| Item | Codes |
|---|---|
| `HOUSE <dim>` | 0 1D (n/a), 1 2D (POINT2), 2 3D (POINT3) |
| `HOUSE <imp>` | 0 FV, 1 FFV, 2 FI |
| `INTGEN <type>` | 0 clear, 1 FV, 2 FI-EVBN (MSM), 3 FI-FSIN (SM), 4 surface, 5 FFV (+ `[level skip]`, default 1) |
| `ANALYS <mode>` | 0 Initiation (Mode 1), 1 New Structure (Mode 2), 2 New Seismic Environment (Mode 3), 3 New Dynamic Loading (Mode 3) |
| `ANALYS <type>` | 0 seismic (FILE1), 1 foundation vibration (FILE9) |
| `ANALYS <impe>` | 0 none, 1 six diagonal terms, 2 full 6×6 |
| Interpolation option (MOTION `<interp>`, STRESS `<interopt>`) | 0 SASSI2000 dense overlapping windows, weighted averaging; 1 original SASSI 1982 non-overlapping windows; 2 dense overlapping, averaging; 3 three overlapping windows, averaging; 4 non-overlapping, one-position shift; 5 non-overlapping, two-position shift; 6 complex bicubic spline |
| Unlagged coherency model (`WPASS <cohf>`, D-INC-01) | 1 Luco–Wong 1986; 2 Abrahamson 1993; 3 Abrahamson 2005 (surface); 4 Abrahamson 2006 (embedded); 5 Abrahamson 2007 hard rock; 6 Abrahamson 2007 soil surface; 7 user-defined (COHXUSER…) |
| `INCOH <nmodes>` | 0 all modes; n > 0 first n (largest λ); −n only mode n |
| Superposition mode | Linear (AS) / Quadratic (SRSS TF) |
| `MERGESOIL [Mode]` | 0 unbonded, 1 merge coincident nodes (default), 2 stiff springs (1e7), 3 stiff below SepLevel / soft (10) above |
| `MOPT` | `<incomp>` 0 include / 1 suppress; `<matrix>` 0 mass / 1 weight; `<mass>`,`<force>` 0 add / 1 overwrite |
| EOUT code digit | 0 no request, 1 print maximum, 2 print maximum + save history |
| SACC `<opt>` | 0 none, 1 maximum, 2 maximum + history |
| BINOUT `reldisp` | 0 none, 1 TFD (unused), 2 THD |
| BBC type / panel force / spring force | 1 CMS, 2 CMB (n/a), 3 TAK, 4 GMR; CHECK requires panel force 1, disp 1, spring force 4 |
| NLSOIL `<NLDampType>` | 1 frequency independent, 2 visco-elastic, 3 Rayleigh |
| FRAMECOMBIN `<op>` | 0 SRSS, 1 sum, 2 average |
| PROCFRAME `<anitype>` | 0 bubble, 1 vector, 2 contour, 3 time history (sorting tag only) |
| ELECOLOR | 1 group, 2 material, 3 property |
| Toggles | ELENUM/GROUPNUM/NODENUM/SHOWMASS/SHRINK/STIPPLE/WIREFRAME: −1 toggle, 0 off, 1 on; DEBUG 2 toggle; PAUSE −1 toggle, 0 start, 1 stop |

### 1.7 Files

"Text"/"binary" follows Figure 1.1. Binary formats are proprietary in the original; here every binary file is a
versioned `.npz`-format container written **under the legacy file name without an added extension**
(D-FIL-01).

**1.7.1 Module input decks (written by AFWRITE as `<model>.<ext>` in the model directory)**

| File | Module | Notes |
|---|---|---|
| `.equ` | EQUAKE | output listing `<model>_equake.out` |
| `.soi` | SOIL | |
| `.sit` | SITE | also required by HOUSE in the working directory |
| `.poi` | POINT2/3 | |
| `.hou` | HOUSE | line 1 col 1 = 1 enables the node optimizer; Option AA `.hou` is a topology-only placeholder |
| `.frc` | FORCE | |
| `.anl` | ANALYS | |
| `.mot` | MOTION | |
| `.str` | STRESS | |
| `.rdi` | RELDISP | |
| `.eql` | NONLINEAR | keyword deck (D-NON-01) written by AFWRITE (D-W3-09); P2 |
| `.lgn` | LOADGEN | written by RUNLOADGEN from the LOADGEN records (D-W3-11); P2 |
| `.nls` | SOIL (SOIL-NON) | NLSOIL / NLSLAYER side file next to `.soi`, written by AFWRITE (D-W3-10); P2 |
| `.pin` | HOUSE (nonlinear soil) | user-edited initial nonlinear-soil properties (P1) |
| `.liq` | HOUSE | iteration flag; "1" → read FILE74 next run (P1) |
| `<MODULE>.inp` | any module (batch) | exactly three lines: `modelname`, `modelname.<ext>`, `modelname_<MODULE>.out` |

**1.7.2 Inter-module numbered files**

| File | Kind | Producer → consumer(s) | Content |
|---|---|---|---|
| FILE1 (FILE1X, FILE1Y, FILE1Z) | bin | SITE Mode 2 → ANALYS | Free-field motion at all interfaces per frequency for unit control motion; wave data |
| FILE2 | bin | SITE Mode 1 → SITE Mode 2, POINT | Generated sublayers, Rayleigh/Love eigenvalues and mode shapes per frequency |
| FILE3 | bin | POINT → ANALYS | Point-load solutions (mode amplitudes, axis values, R0) per frequency, load interface 1…L+1, harmonic μ = 0, 1 |
| FILE4 = `<model>.N4` | bin | HOUSE → ANALYS, STRESS | Topology, DOF map, element data (recovery operators), excavation flags |
| COOSK, COOSM | bin | HOUSE → ANALYS | Assembled sparse stiffness (complex) and mass, split into structure and excavated parts |
| DOFSMAP, FILE90, FILE91 | bin | HOUSE → ANALYS (restart) | DOF map; model and interaction-set hashes; restart metadata |
| FILE77 (FILE77001…FILE77050) | bin | HOUSE → ANALYS | Incoherency / multiple-excitation factors per frequency and direction (X, Y, Z in each file); one per simulation |
| FILE78 | text | HOUSE → STRESS | Nonlinear-soil iteration data (properties used) |
| FILE74 | text | STRESS (COMB_XYZ_STRAIN) → HOUSE | Effective strains and strain-compatible G, β per nonlinear soil element |
| FILE73 | text | SOIL → STRESS | Soil material curves (G/Gmax–γ, D–γ) |
| FILE88 | text | SOIL → SITE (and HOUSE excavated layers via UI) | Strain-compatible layer properties |
| FILE9 (FILE9001…FILE9500) | bin | FORCE → ANALYS | Load vectors per SSI frequency; one per load case |
| FILE8 (FILE8X/Y/Z; FILE8001…FILE8150 incoherent; FILE8001…FILE8500 load cases) | bin | ANALYS, COMBIN → MOTION, STRESS, COMBIN | Complex TFs at all DOFs and SSI frequencies; DOF map; frequency numbers |
| FILE81, FILE82 | bin | user rename of FILE8 → COMBIN | The two solutions merged |
| COOXxxx, COOTKxxx, COOXI, COOTKI | bin | ANALYS → ANALYS restart | xxx = 3-digit frequency order number; COOX = X_ff; COOTK = factorised condensed system; index files |
| FILE11 | bin | ANALYS | Global unconstrained impedance data (large) |
| FOUNSTIF, FOUNDASH, FOUNDAMP, FOUNIMPD | text | ANALYS | Re K; Im K/ω; Im K/(2\|Re K\|); \|K\| per frequency |
| FILE12 | text | MOTION | Acceleration histories and RS (optional) |
| FILE13 | text | MOTION | Baseline-corrected `time acc vel disp` per node/DOF (optional) |
| FILE14, FILE15 | text | STRESS | Beam force/moment TFs; stress time histories (listings) |

**1.7.3 Result and auxiliary text files**

| File / pattern | Producer | Content |
|---|---|---|
| `.rsi` (any name) | user | EQUAKE target RS: 2 columns (Hz, amplitude), one damping |
| `.acc`, `.vel`, `.dis` | EQUAKE | generated histories (accel in g) |
| `.rso` | EQUAKE | RS of generated motion; also `name.RSO` from MOTION external-file conversion |
| `.psd` | EQUAKE | 2 columns: Hz, PSD (cm²/s³ SI, in²/s³ British) |
| `.fft` | EQUAKE | 3 columns: Hz, Re, Im (positive frequencies, strong-motion part) |
| `ACCxxx.TH`, `SNxxx.TH`, `SSxxx.TH` | SOIL | layer acceleration (top), strain (mid), stress (mid); xxx = 3-digit layer number from the surface |
| `nnnnnTR_d.TFU` / `.TFI` / `.ACC` | MOTION | node (5 digits, 6 if > 99,999 nodes), DOF tag TR_X/TR_Y/TR_Z or R_XX/R_YY/R_ZZ |
| `nnnnnTR_dzz.RS` | MOTION | RS for damping order number zz (01, 02, …) |
| `nnnnnTR_d.TFD` / `.THD` | RELDISP | relative-displacement TF / history |
| `etype_ggg_eeeee_comp.TFU` / `.TFI` / `.THS` | STRESS | element STF and stress history, e.g. `BEAMS_003_00045_MXJ.THS` |
| `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT` | STRESS | maximum element-centre components (groups header then element rows) |
| `STATIC_SOIL_PRESSURES.TXT`, `pres_max_ele`, `pres_max_nod` | STRESS | static pressures to add; element and averaged nodal maxima |
| `.sig`, `.tau`, `.bdsig`, `.bdtau` | STRESS | all-element maxima / histories by component class |
| `ESTRESS_<frame>.ess` | STRESS (`SECDATAOPT,1`) | element-centre stresses of one time step (section cuts) |
| `Frames.txt` | user / STRESS dialog | `[# frames]`, frame numbers, `[# soil-pressure groups]`, groups on one line |
| `SRSSTF.txt` | user / MOTION dialog | `[# modes] [phase 0/1]`, then FILE8 names (coherent first if phase = 1) |
| `CONTTRS.txt` | user / MOTION dialog | list of external `.ACC` files for RS conversion |
| `ACC_max.txt` | MOTION | ZPA frame |
| `COHXUSER`, `COHYUSER`, `COHZUSER`, `FREQCOH`, `DISTCOH` | user | user coherency tables (default 100 × 100) |
| `<model>.hounew`, `<model>.map` | HOUSE optimizer | renumbered deck; `old new` node pairs |
| `<model>_Excv.map` | MERGESOIL | excavation node map (Option AA name) |
| `PANEL.NON`, `SPRING.NON` | NONLINEAR | state flag "1" |
| `Panelxxxx.crv/.thd/_AXIAL.thd/.ths`, `Panel_EQL_Matl_Prop.txt`, `Panel.fmu` (SPRING… analogues), `Modelname_new.hou` | NONLINEAR | P2 outputs |
| `COMB_XYZ_THD.inp` | user | P2 combination control (format defined in D-NON-08) |

**1.7.4 Frame directories (post-processing "Restart" options)**

| Directory / pattern | Producer | Content |
|---|---|---|
| `\TFU\TFU_<fff.ff>_<nnnnn>` | MOTION Restart for TF | complex TF of all active nodal DOFs at one frequency |
| `\ACC\ACC_<tt.ttt>_<nnnnn>`, `\ACCR\…` | MOTION Restart for ACC | accelerations (rotations) at one time step |
| `\RS\RS<dd>_<fff.ff>_<nnnnn>` | MOTION Restart for RS | RS values at one frequency, damping dd |
| `\THD\THD_<tt.ttt>_<nnnnn>`, `\THDR\…` | RELDISP Restart | relative displacements (rotations) at one step |
| `\NSTRESS\stress_<tt.ttt>_<nnnnn>_<sig\|tau\|bdsig\|bdtau>`, `stress_ABS_MAX_<comp>` | STRESS | averaged nodal stresses (plotting only) |
| `\SOILPRES\pres_<tt.ttt>_<nnnnn>_<ele\|nod>`, `pres_ABS_MAX_<type>` | STRESS | soil pressures |
| `*.dispani`, `*.tfiani`, `*.impani` (legacy `zpani`, `thiani`, `contani`, `thani`) | user/UI | animation frame list (first line ignored) |

**1.7.5 UI, persistence and exchange files**

| File | Meaning |
|---|---|
| `<model>.pre` | Command (PREP/UI) text file; written by WRITE, read by INP |
| `<model>-Sim.pre` | Simulation commands written when Extended Write "Simulation Commands" is set |
| `<model>.err` | CHECK output (`Errors and Warnings for <MODULE>`, `Error n : text`) |
| `<model>.sdb` | SAVE/RESUME binary model database (D-MDL-02) |
| `SASSIini.xml` | module locations, extensions, display settings |
| `SASSIdb.xml` | group → model tree (name and location only) |
| `SASSIani.xml` | processed-animation registry |
| `<model>.inp` (APDL) | ANSYS export (`ANSYS` command) — distinct from the 3-line batch `.inp` |
| `.cdb` | ANSYS CDWRITE file (the manual's "CBD file") read by CONVERT,ANSYS |
| `Modelname_ACC.bin`, `Modelname_STRESS.bin`, `Modelname_TR_[X\|Y\|Z]_THD.bin`, `Modelname_THD.bin` | binary history databases (P2) |

### 1.8 Naming formulas (normative)

```
nodal result     f"{node:0{w}d}{tag}.{ext}"   w = 5 (6 if max node > 99,999); tag ∈ TR_X TR_Y TR_Z R_XX R_YY R_ZZ
nodal RS         f"{node:0{w}d}{tag}{idamp:02d}.RS"
element result   f"{etype}_{group:03d}_{elem:05d}_{comp}.{ext}"      etype ∈ SOLID BEAMS SHELL TSHELL PLANE SPRING
soil layer       f"ACC{layer:03d}.TH"  f"SN{layer:03d}.TH"  f"SS{layer:03d}.TH"
incoherent FILE8 f"FILE8{3*(s-1)+d:03d}"      s = 1..50, d = 1 X, 2 Y, 3 Z
load-case FILE8  f"FILE8{l:03d}"  ; FILE9 f"FILE9{l:03d}" ; FILE77 f"FILE77{s:03d}"
restart          f"COOX{q:03d}"  f"COOTK{q:03d}"   q = 1-based order in the sorted frequency list
frames           freq f"{f:06.2f}"  time f"{t:06.3f}"  frame f"{n:05d}"
```
Element component codes (D-FIL-04): SOLID `SXX SYY SZZ SXY SXZ SYZ SOCT`; BEAMS `FXI FYI FZI MXI MYI MZI FXJ FYJ
FZJ MXJ MYJ MZJ` (local 1/2/3 written as X/Y/Z); SHELL `FXX FYY FXY MXX MYY MXY`; TSHELL `NXX NYY NXY QXZ QYZ MXX
MYY MXY`; PLANE `SXX SZZ TXZ`; SPRING `FX FY FZ MXX MYY MZZ`.

### 1.9 CHECK and messaging terms

| Term | Meaning |
|---|---|
| Error n / Warning n | Numbered messages of Manual Ch. 10 (Errors 1–128, Warnings 1–11), titles exact (spec 11 §5) |
| Fatal error | Errors 42 (no nodes) and 43 (no groups) stop CHECK |
| Break Check at N | Maximum printed messages per type per module (default 100); totals still counted |
| AFWRITE gating | A module with any CHECK error does not get its deck written; warnings never block |
| EDU-nn | Additional checks of this implementation, numbered separately so manual numbering stays faithful |
| Message classes | command echo, confirmation, comment, information (cannot be hidden), warning, error |

---

## 2. Module data flow

### 2.1 Module input/output table (normative)

All files live in the model directory (the module's working directory). "Req." files must exist or the
module stops with an error naming the producer module (D-RUN-03).

| Module | Deck | Reads (required) | Reads (conditional) | Writes | Consumed by |
|---|---|---|---|---|---|
| EQUAKE | `.equ` | RSIN target spectra (≤ 3) | ACCIN seed/external histories; TPSD target PSD | ACCOUT `.acc/.vel/.dis`; RSOUT `.rso`; `.psd`; `.fft`; listing | SOIL, MOTION, STRESS, RELDISP (as THFILE) |
| SOIL | `.soi` | THFILE control motion; DYNP curves (in deck) | — | listing; `ACCxxx/SNxxx/SSxxx.TH`; FILE73; FILE88 (if `<save>`=1) | SITE (FILE88), STRESS (FILE73), user |
| SITE | `.sit` | — | FILE2 (Mode 2 only, if Mode 1 off); FILE88 (Non-Linear Soil flag) | FILE2 (Mode 1); FILE1 (Mode 2); listing | POINT, ANALYS (user copies FILE1 → FILE1X/Y/Z) |
| POINT2/3 | `.poi` | FILE2 | — | FILE3; listing | ANALYS |
| HOUSE | `.hou` | `.sit` (layer table; always) | FILE74 + `.liq`=1 (nonlinear soil); `.pin`; COHXUSER…/FREQCOH/DISTCOH (model 7) | FILE4 (`<model>.N4`), COOSK, COOSM, DOFSMAP, FILE90, FILE91; FILE77/FILE77sss (incoherent or ME); FILE78 (nonlinear soil); `.liq`; `.hounew`, `.map` (optimizer); listing (incl. `I N C O` table) | ANALYS, STRESS |
| FORCE | `.frc` | — | — | FILE9; listing | ANALYS (user copies to FILE9lll) |
| ANALYS | `.anl` | FILE3, FILE4, COOSK, COOSM, DOFSMAP; seismic: FILE1 (or FILE1X/Y/Z); vibration: FILE9 (or FILE9001…) | FILE77/FILE77sss; restart: COOXqqq, COOTKqqq, COOXI, COOTKI, FILE90, FILE91 | FILE8 (or FILE8X/Y/Z, FILE8nnn); COOX/COOTK (if `<save>`=1); FOUN* + FILE11 (if `<impe>`>0); listing | COMBIN, MOTION, STRESS |
| COMBIN | (none; D-CMB-02) | FILE81, FILE82 | — | FILE8; listing | MOTION, STRESS |
| MOTION | `.mot` | FILE8 (or the file named in the deck); THFILE | SRSSTF.txt + modal FILE8s; CONTTRS.txt + `.ACC` files | `.TFU/.TFI/.ACC/.RS`; FILE12/FILE13 (optional); `.RSO`; frames `\TFU \ACC \ACCR \RS`, `ACC_max.txt`; `Modelname_ACC.bin` (P2); listing | RELDISP (`.TFI`), UI plots, LOADGEN |
| RELDISP | `.rdi` | complex `.TFI` of reference and output nodes (MOTION `<cplx>`=1); THFILE | — | `.TFD`, `.THD`; frames `\THD \THDR`; `Modelname_TR_d_THD.bin` (P2); listing | UI, COMB_XYZ_THD |
| STRESS | `.str` | FILE4, FILE8, THFILE | FILE78 + FILE73 (nonlinear soil); `Frames.txt`; `STATIC_SOIL_PRESSURES.TXT` | `.TFU/.TFI/.THS`; FILE14, FILE15; FILE74; `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`; `.sig/.tau/.bdsig/.bdtau`; `\NSTRESS`, `\SOILPRES`, `pres_max_*`; `ESTRESS_n.ess`; `Modelname_STRESS.bin` (P2); listing | HOUSE (FILE74), UI |
| COMB_XYZ_STRAIN | (command, D-NLS-04) | three directional FILE74s | — | combined FILE74 | HOUSE |
| NONLINEAR | `.eql` | `.sit`, `.hou`, combined `.THD` histories | `PANEL.NON`/`SPRING.NON`, `*_EQL_MATL_PROP.txt` | `.crv/.thd/.ths/.fmu`, `*_EQL_Matl_Prop.txt`, `Modelname_new.hou`, `.NON` | HOUSE |
| COMB_XYZ_THD | `COMB_XYZ_THD.inp` | X/Y/Z-input `.THD` sets | — | combined `.THD` | NONLINEAR |

### 2.2 Figure 1.1 (transcribed)

```
 INPUT                    ANALYSIS                                  RESULTS
 FORCE --(9)-----------\
 EQUAKE                 \        Restart: COOXqqq, COOTKqqq, COOXI, COOTKI,
   |(acc file)           \       DOFSMAP, FILE90, FILE91 (in/out)
 SOIL --(73)--------------\---------------------------\
   |(88)                   \                           \
 SITE --(1)-------------> ANALYS --(8)--> COMBIN --+--> MOTION --> .TFU .TFI .ACC .RS ; FILE12,13
   |(2)                   /                         |      | (complex .TFI)
 POINT --(3)-----------> /                          |    RELDISP --> .THD .TFD --> COMB_XYZ_THD -> NONLINEAR --> HOUSE
 HOUSE --(4,77)-------> /                           +--> STRESS --> .TFU .TFI .THS ; FILE14,15
   ^  \---(78)----------------------------------------> STRESS
   |<-----------------------------------------(74)----- STRESS
   |<-- SSI2ANSYS <-- ANSYS model                                    [Option AA]
 Results --> LOADGEN --> ANSYS model                                 [Option A]
```

### 2.3 File producer/consumer matrix (W = writes, R = reads, r = reads conditionally)

| File | EQUAKE | SOIL | SITE | POINT | HOUSE | FORCE | ANALYS | COMBIN | MOTION | RELDISP | STRESS |
|---|---|---|---|---|---|---|---|---|---|---|---|
| THFILE history | W | R | | | | | | | R | R | R |
| FILE88 | | W | r | | (UI) | | | | | | |
| FILE73 | | W | | | | | | | | | r |
| FILE2 | | | W/r | R | | | | | | | |
| FILE1, FILE1X/Y/Z | | | W | | | | R | | | | |
| FILE3 | | | | W | | | R | | | | |
| `.sit` | | | R | | R | | | | | | |
| FILE4, COOSK, COOSM, DOFSMAP, FILE90/91 | | | | | W | | R | | | | R (FILE4) |
| FILE77 / FILE77sss | | | | | W | | r | | | | |
| FILE78 | | | | | W | | | | | | r |
| FILE74 | | | | | r | | | | | | W |
| FILE9 / FILE9lll | | | | | | W | r | | | | |
| COOXqqq, COOTKqqq, COOXI, COOTKI | | | | | | | W/r | | | | |
| FILE8 family | | | | | | | W | W | R | | R |
| FILE81, FILE82 | | | | | | | | R | | | |
| `.TFI` (complex) | | | | | | | | | W | R | |

### 2.4 Run sequences

```
LINEAR SEISMIC (initiation):
  [EQUAKE] → [SOIL] → SITE(Mode 1, Mode 2) → POINT → HOUSE → ANALYS(Mode 1) → [COMBIN] → MOTION → RELDISP ; STRESS
  (MOTION and STRESS are independent; RELDISP requires MOTION .TFI with <cplx>=1)
LINEAR SEISMIC, X+Y+Z coherent in one ANALYS run:
  SITE ×3 (SV x′ / SH y′ / P z′, angle 0) → copy FILE1 → FILE1X, FILE1Y, FILE1Z ; POINT ; HOUSE ;
  ANALYS(simul=1, ang=0) → FILE8X/Y/Z ; MOTION/STRESS per direction (FILE8 = FILE8d)
FORCED VIBRATION:
  SITE(Mode 1 only) → POINT → FORCE (×Nl, copy FILE9 → FILE9lll) → HOUSE → ANALYS(type 1) → MOTION/STRESS
  (SOIL not used; low-strain properties; FILE9 replaces FILE1)
FREQUENCY REFINEMENT:
  MOTION → CRITFREQ (TFU vs TFI) → add frequencies (must exist in FILE1/FILE9 and FILE3, else re-run SITE/POINT/FORCE)
  → ANALYS (new subset) → rename old/new FILE8 → FILE81/FILE82 → COMBIN → MOTION/STRESS   (repeat)
RESTARTS: see 2.5
INCOHERENT, STOCHASTIC (P1):
  SITE ×3 → FILE1X/Y/Z ; POINT ; HOUSE(coh=1, Ns ≤ 50 → FILE77001…) ; ANALYS(simul=Ns) → FILE8001…FILE8(3Ns) ;
  MOTION per FILE8 → average responses (AVERAGE)
INCOHERENT SRSS TF (P1, benchmarking only):
  for k = 1..m: HOUSE(nmodes = −k) ; ANALYS Mode 3 → rename FILE8_k ; then MOTION with SRSSTF.txt
NONLINEAR SOIL SSI ITERATIONS (P1):
  initiation SITE → POINT → HOUSE(.pin) → ANALYS(save) → STRESS(X,Y,Z) → COMB_XYZ_STRAIN
  loop: HOUSE(.liq=1 reads FILE74) → ANALYS(Mode 2, X,Y,Z) → STRESS(X,Y,Z) → COMB_XYZ_STRAIN → converged?
OPTION NON (P2):
  SITE → POINT → HOUSE → ANALYS(save) → MOTION → RELDISP (×dir) → COMB_XYZ_THD → NONLINEAR (writes PANEL.NON)
  loop: HOUSE(Modelname_new.hou) → ANALYS(Mode 2) → MOTION → RELDISP → COMB_XYZ_THD → NONLINEAR → converged?
AFWRBAT SPLIT (P1):
  k folders <name>_k with contiguous frequency blocks: SITE→POINT→HOUSE→ANALYS each; combine.bat folds COMBIN pairwise
```

### 2.5 Restart and re-analysis types

| Change | Name | Modules to re-run | ANALYS `<mode>` (manual Mode) | Reused | Preconditions (CHECK at run time) |
|---|---|---|---|---|---|
| Control-motion history or spectrum only, same wave field | New Time History | MOTION (STRESS, RELDISP) | — | FILE8 | FILE8 present; Δf of new NFFT·Δt equals FILE8 Δf |
| Structure / near-field properties or geometry; interaction nodes and soil unchanged | New Structure | HOUSE → ANALYS → MOTION/STRESS | 1 (Mode 2) | X_ff (COOXqqq) | COOX present for every requested frequency; FILE90 hash of interaction-node coordinates/numbering and layer table unchanged |
| Wave type / incidence / control location; new incoherent sample | New Seismic Environment | SITE (or HOUSE if incoherent) → ANALYS | 2 (Mode 3) | X_ff + factorised system (COOX + COOTK) | as above, plus FILE4 hash unchanged |
| New external-load pattern | New Dynamic Loading | FORCE → ANALYS → MOTION | 3 (Mode 3) | X_ff + factorisation | as above |
| Only the load time history changes | (New Time History) | MOTION | — | FILE8 | — |
| Special seismic load vector | New Load Vector | — | Mode 6 (batch only) | — | OOS |

Every restart result must equal a fresh initiation run to round-off (VP-22).

### 2.6 Pipeline-runner dependency rules (normative)

| Before running | Require | Message if missing |
|---|---|---|
| any RUNxxx | model name and path (MDL or database); deck `<model>.<ext>` written by AFWRITE; no CHECK error for that module | "Model name/path not defined — use MDL" / "Run AFWRITE first" |
| POINT | FILE2 | "FILE2 missing — run SITE Mode 1" |
| SITE Mode 2 alone | FILE2 | as above |
| HOUSE | `.sit` in the model directory | "HOUSE needs <model>.sit (run AFWRITE with SITE enabled)" |
| ANALYS seismic | FILE1 (or FILE1X/Y/Z when simul = 1), FILE3, FILE4 family | name the producing module |
| ANALYS vibration | FILE9 (or FILE9001…FILE9Nl), FILE3, FILE4 family | idem |
| ANALYS incoherent / ME | FILE77 (or FILE77001…FILE77Ns) | idem |
| ANALYS frequency survey | every requested frequency number present in FILE1/FILE9 **and** FILE3 | "Frequency n not in FILE1/FILE3 — re-run SITE/POINT" (stop) |
| COMBIN | FILE81, FILE82 with identical DOF maps and Δf | error |
| MOTION / STRESS | FILE8 (or named FILE8 file); analysis type matches FILE8 | error |
| RELDISP | complex `.TFI` for the reference and every requested node/DOF | "Run MOTION with Save Complex TF for nodes …" |
| STRESS | FILE4 matching FILE8 DOF map | error |
| HOUSE optimizer used | MOTION/RELDISP node requests are in the new numbering | warning listing nodes absent from `.map` targets |

---

## 3. The `.pre` command language

### 3.1 Lexical and syntax rules (normative; consolidates specs 01 §11, 04 §4–6, 07 §1, 08 §1, 09 §1, 10 §1)

| # | Rule | Source / decision |
|---|---|---|
| L1 | A command line is `Keyword, p1, p2, …, pn`, with commas as delimiters. A single run of blanks or tabs may replace the **first** comma only (`EDGE 1,0,0,1`) | [M] / D-PAR-02 |
| L2 | Whitespace around every token is trimmed | [M] |
| L3 | Command names are case-insensitive. String arguments keep their case. DYNP labels and SOILPROPPLOT names are case-sensitive. Macro names are upper-cased. Variable names are case-insensitive | [M] |
| L4 | A line whose first non-blank character is `*` is a comment. Blank lines are ignored. There are no inline comments (a `*` inside a line is data) | [M] / D-PAR-03 |
| L5 | An empty field (nothing or only blanks between commas) and an omitted trailing field mean "use the documented default". An omitted **required** numeric field is read as 0 with a warning (PREP blank-is-zero convention) | [M] / D-PAR-05 |
| L6 | Numbers are free format; Fortran exponents `E`/`D` accepted (`1.0D9`, `1e+008`). Integer fields given as reals are rounded (warning if the fractional part ≠ 0) | D-PAR-07 |
| L7 | Text arguments in the **last** position (TIT, THTIT, EQTIT, GTIT, SSAF title, PLOTTITLE, XTITLE, YTITLE, YTITLE2, LINENAME name, DYNP label, PROCFRAME data when last) take the rest of the line, commas included. A double-quoted token `"…"` may be used anywhere to embed commas (extension) | D-PAR-06 |
| L8 | Node and element lists (NOUT, EOUT, CUTADD list form, hide lists) are trailing tokens of integers or inclusive ranges `a-b` (`a ≤ b`), separated by commas, blanks, tabs or `;` | D-PAR-08 |
| L9 | Legacy SOIL decks: a final parenthesised token `(…)` is one token (inner commas do not split). The command is flagged legacy, its format string kept, and a warning (Error 102 context) is issued; arguments are not auto-mapped | [M] / D-SOL-10 |
| L10 | Command lookup: exact match of the full name, else exact match of the documented 4-character abbreviation (§3.2). No prefix matching. Full names take priority (GROUPMAT vs GROUP, ETYPEGEN vs ETYPE) | [M] / D-PAR-04 |
| L11 | UI commands (manual §9.6 onward) and extension commands (§3.4.R) must be typed in full | [M] |
| L12 | An unknown name prints `<X> Command not found` and processing continues. A failing command prints its error and processing continues; `.pre` execution never aborts except on an unreadable file. A summary (commands, warnings, errors) is printed at EOF | [M] / D-PAR-09 |
| L13 | Per line, substitution order is: (1) macro `$k$` placeholders as plain text; (2) variable expressions `@NAME…` and `#`, left to right, each evaluated once; (3) tokenisation; (4) dispatch. A FOREACH body is substituted afresh at each iteration | [M] / D-PAR-10 |
| L14 | Maximum line length after substitution 3000 characters; macro names ≤ 50 characters; macro nesting allowed (depth guard 64); INP nesting allowed (depth guard 32) | [M] / D-PAR-11 |
| L15 | Relative paths resolve in this order: absolute; active model path (if MDL set); current working directory (CD/MDL); directory of the calling file | [M] / D-PAR-12 |
| L16 | While INP/macro input runs, command echo and confirmations are suppressed; comments, warnings, errors and information are shown. At EOF print `INPUT FILE REACHED EOF, INPUT SWITCHED TO KEYBOARD` | [M] |
| L17 | Every GUI action emits the equivalent command text through the same interpreter (logged in the Command History), so a session can be replayed exactly | [M] |
| L18 | Commands that open a dialog when called without arguments (CONVERT, CUTPLOT, SOILPROPPLOT, PROCFRAME, animation plots, SHADEROPTIONS, SHOWDOF, WINDOWSETTINGS) raise an error in batch/headless mode | D-UI-08 |

**Variables and loops (P1)** [M]: `VAR,<Name>,<X1>…` defines a list of strings and resets the counter to 0.
`@X` gives the counter; `@X+k`, `@X-k`, `@X++`, `@X--`, `@X=k` change the counter first and substitute the
new value (pre-increment); `@X[i]` is the 1-based i-th item; inside FOREACH, `#` is the loop index (1-based)
and `@V[#]` uses the index of the loop over V. `FOREACH,<Var>,<command…>` runs one command per item;
nesting a FOREACH over the same variable is an error and the loop does not run. Regression test: Loop.pre
creates 125 nodes, node `25(iz−1)+5(iy−1)+ix` at `(ix,iy,iz)` (UT-10).

**Macros (P1)** [M]: `LOADMACRO,<Name>,<file>` caches the file; `MACRO,<Name>,<a1>,…` replaces every `$k$`
(also inside words, e.g. `Node$1$x.rs`) before variable processing. A missing argument substitutes an empty
string with a warning (D-PAR-13). Macros live in their own namespace (a macro named like a command does not
override it).

### 3.2 Accepted abbreviations

| Full | Abbr | Full | Abbr | Full | Abbr | Full | Abbr |
|---|---|---|---|---|---|---|---|
| ACCIN | ACCI | HOUSE | HOUS | THTIT | THTI | MLIST | MLIS |
| ACCOUT | ACCO | INCOH | INCO | WPASS | WPAS | MTYPE | MTYP |
| AFWRITE | AFWR | LFREQ | LFRE | WRITE | WRIT | MXDEL | MXDE |
| ANALYS | ANAL | MOTION | MOTI | GLOBAL | GLOB | MXLIST | MXLI |
| CHECK | CHEC | POINT | POIN | LOCAL | LOCA | RLIST | RLIS |
| EQTIT | EQTI | RESUME | RESU | NLIST | NLIS | SCLIST | SCLI |
| EQUAKE | EQUA | RSOUT | RSOU | NMOVE | NMOV | THICK | THIC |
| FORCE | FORC | STATUS | STAT | NSCALE | NSCA | FLIST | FLIS |
| ECOMPR | ECOM | ELIST | ELIS | GLIST | GLIS | FSCALE | FSCA |
| GROUP | GROU | LLIST | LLIS | MMDEL | MMDE | MMLIST | MMLI |
| MRGEN | MRGE | MRDEL | MRDE | MRSCALE | MRSC | MSCALE | MSCA |
| MTDEL | MTDE | MTGEN | MTGE | MTLIST | MTLI | MTSCALE | MTSC |
| MUNITS | MUNI | THFILE | THFI | | | | |

- No abbreviation: RELFILE (whole name underlined) and every command of 4 characters or fewer.
- Legacy aliases accepted by decision D-PAR-04 (no underline printed in the manual): INTLIST → INTL,
  LMOVE → LMOV, SLIST → SLIS, DELSC → DELS, ETYPE → ETYP, STRESS → STRE.
- Name aliases: FIXSPROT → FIXSPRROT; PANELGEN → PNLGEN; NONLINMODISP → NONLINMOTDISP; ACCANIDB → ACCDBANI;
  LOADACCDBANI → LOADACCDB; COMDISPDB → COMBDISPDB.

### 3.3 Command classes (state behaviour)

| Class | Commands | Behaviour |
|---|---|---|
| Record setter | ANALYS, AOPT, EQUAKE, FORCE, HOUSE, INCOH, MOPT, MOTION, POINT, SITE, SOIL, STRESS, WPASS, RELD, EQL, NLSOIL, BINOUT (blank = unchanged), extension X-commands | replaces the module's option record (every positional argument stored) |
| Indexed setter | ACCIN, ACCOUT, RSIN, RSOUT, TPSD, CORR, DYNP, ME, SYMM, WAVE, SPRO, SACC, SRS, SSTR, SSAF, SFOU, L, M, R, SC, MX*, N, E, LOC, LOCAL, BBC*, P, S, NLSLAYER | creates or overwrites the entry at its key |
| List append / delete | FREQ (per set), DAMP, TOPL, AMP (per motion) | non-zero first value appends; first value 0 deletes the list |
| Request append | NOUT, EOUT, RDND | appends; `NOUT,0` / `EOUT,0` / `RDND,0` clears the list (D-PAR-14) |
| String setter | TIT, THTIT, EQTIT, THFILE, RELFILE, GTIT | sets a string |
| Action | INP, WRITE, SAVE, RESUME, CHECK, AFWRITE, RUNxxx, model checks, generators, plots, calculations | performs an operation |

### 3.4 Command inventory by category

Columns: **Tier** as defined in §0.1. "Spec" is the normative detail. Syntax uses `<required>` and
`[optional]`.

#### 3.4.A Session, files, models and global options

| Command | Abbr | Syntax | Tier | Meaning | Spec |
|---|---|---|---|---|---|
| INP | — | `INP,<filename>` | P0 | Execute commands from a `.pre` file | 07 §9.2.18 |
| WRITE | WRIT | `WRITE,[<file>],[<path>]` | P0 | Write the model (incl. analysis options) as commands; default `<model>.pre` in the model path | 07 §9.2.47 |
| SAVE | — | `SAVE` | P0 | Save the active model (binary, `<model>.sdb`) | 07 §9.2.29 |
| RESUME | RESU | `RESUME` | P0 | Reload the last SAVE | 07 §9.2.25 |
| MDL | — | `MDL,<Model>,<Path>` | P0 | Set model name and directory (also sets the working directory) | 09 §2.20 |
| MDLNAME | — | `MDLNAME,<name>` | P0 | Change the model name only | 09 §2.21 |
| TIT | — | `TIT,<title>` | P0 | Model title | 07 §9.2.43 |
| STATUS | STAT | `STATUS` | P0 | Global model information | 07 §9.2.37 |
| ACTM | — | `ACTM,<Model>` | P0 | Activate model number (create empty if absent) | 09 §2.1 |
| CPMODEL | — | `CPMODEL,<Mdl>` | P0 | Copy active model to `<Mdl>` (name/path copied, warning) | 09 §2.3 |
| DMODEL | — | `DMODEL,<Mdl>` | P0 | Remove a model from memory only | 09 §2.5 |
| MODELLIST | — | `MODELLIST` | P0 | List models in memory | 09 §2.23 |
| CD | — | `CD,<dir>` | P0 | Change working directory (must exist) | 10 §5.2 |
| MKDIR | — | `MKDIR,<dir>` | P0 | Create directory (warning if it exists) | 10 §5.8 |
| MOPT | — | `MOPT,<incomp>,<matrix>,<mass>,<force>` | P0 | Incompatible modes 0 include/1 suppress; GM mass 0 mass/1 weight; mass and force 0 add/1 overwrite | 07 §9.2.21 |
| GRAVITY | — | `GRAVITY,<grav>` | P0 | Set HOUSE `<gravity>` only | 09 §2.16 |
| GROUNDELEV | — | `GROUNDELEV,<elev>` | P0 | Set HOUSE `<gelev>` only | 09 §2.17 |
| AOPT | — | `AOPT,<EQUAKE>,<SOIL>,<DEP>,<SITE>,<POINT>,<HOUSE>,<DEP>,<FORCE>,<ANALYS>,<COMBIN>,<MOTION>,<STRESS>,<RELDISP>,<PANEL>` | P0 | Modules processed by AFWRITE and CHECK (DEP slots must be 0) | 07 §9.2.6 |
| CHECK | CHEC | `CHECK` | P0 | Validate model and enabled modules; write `<model>.err` | 07 §9.2.7, 11 §5 |
| AFWRITE | AFWR | `AFWRITE` | P0 | CHECK, then write the deck of each enabled module without errors | 07 §9.2.3 |
| AFWRBAT | — | `AFWRBAT,<splits>` | P1 | Split the frequency set into folders with batch files and a COMBIN batch | 09 §2.2 |
| SETENV | — | `SETENV,<mem>` | P2 | Solver memory limit (MB) | 09 §2.27 |
| GETENV | — | `GETENV` | P2 | Show solver settings | 09 §2.14 |

#### 3.4.B Frequency sets

| Command | Abbr | Syntax | Tier | Meaning | Spec |
|---|---|---|---|---|---|
| FREQ | — | `FREQ,<ndx>,<f1>,…,<f10>` | P0 | Append frequency numbers to set `ndx`; `f1 = 0` deletes the set | 07 §9.2.15 |
| LFREQ | LFRE | `LFREQ,[start],[end],[step]` | P0 | List sets with numbers and Hz values | 07 §9.2.19 |

#### 3.4.C Nodes and coordinate systems (Manual §9.3)

| Command | Abbr | Syntax | Tier | Meaning |
|---|---|---|---|---|
| N | — | `N,<nd>,[<x>],[<y>],[<z>]` | P0 | Define node in the active system (extra legacy fields ignored) |
| NDEL | — | `NDEL,<n1>,[<n2>],[<inc>]` | P0 | Delete nodes (attached loads/masses removed) |
| NGEN | — | `NGEN,[itim],[step],[n1],[n2],[inc],[dx],[dy],[dz]` | P0 | Copy node pattern itim times |
| FILL | — | `FILL,[<n1>],[<n2>],[<nr>]` | P0 | Interpolate nodes between two nodes |
| NMED | — | `NMED,<nd>,<n1>,[<n2>,…,<n8>]` | P0 | Node at the mean of 1–8 nodes |
| LMOVE | LMOV* | `LMOVE,[<dx>],[<dy>],[<dz>],<nd>,<l1>,[…,<l15>]` | P0 | Copy a node list with translation |
| NMOVE | NMOV | `NMOVE,[<dx>],[<dy>],[<dz>],<nd>,<l1>,[…,<l15>]` | P0 | Copy a node list with coordinate scaling (defaults 1) |
| NSCALE | NSCA | `NSCALE,<n1>,<n2>,[<inc>],[<sfx>],[<sfy>],[<sfz>]` | P0 | Scale coordinates in place (0 → 1) |
| NLIST | NLIS | `NLIST,[<n1>],[<n2>],[<inc>]` | P0 | List nodes |
| D | — | `D,<n1>,<n2>,[<inc>],[<val>],<label1>,[…,<label6>]` | P0 | Fixity 0 free (default) / 1 fixed; UX UY UZ ROTX ROTY ROTZ DISP ROT ALL |
| INT | — | `INT,<n1>,<n2>,[<inc>],<set>,[<code>]` | P0 | Set (1) / reset (0) code 0 interaction (default), 1 intermediate, 2 interface, 3 internal |
| INTLIST | INTL* | `INTLIST,[<n1>],[<n2>],[<step>],[<c1>],[<c2>],[<c3>],[<c4>]` | P0 | List classified nodes |
| CSYS | — | `CSYS,<ns>` | P0 | Activate coordinate system (0 = global) |
| LOC | — | `LOC,<ns>,<type>,<x0>,<y0>,<z0>,<txy>,<tyz>,<txz>` | P0 | Cartesian local system by origin and angles (type 1 not available) |
| LOCAL | LOCA | `LOCAL,<ns>,<type>,<n1>,<n2>,<n3>` | P0 | Local system from three nodes |
| GLOBAL | GLOB | `GLOBAL,<n1>,<n2>,<inc>` | P0 | Convert stored coordinates to global |
| SDEL | — | `SDEL,<s1>,[<s2>],[<inc>]` | P0 | Delete systems (nodes converted to global first) |
| SLIST | SLIS* | `SLIST,[<s1>],[<s2>],[<inc>]` | P0 | List systems |

(* = legacy alias by D-PAR-04.)

#### 3.4.D Groups, elements and properties (Manual §9.4)

| Command | Abbr | Syntax | Tier | Meaning |
|---|---|---|---|---|
| GROUP | GROU | `GROUP,<ng>,<type>` | P0 | Create/activate group; type 1 SOLID, 2 BEAMS, 3 SHELL, 4 PLANE, 5 TSHELL, 7 SPRING, 9 GENERAL (number or name) |
| MTYPE | MTYP | `MTYPE,[<gr>],<type>` | P0 | Change group type |
| GTIT | — | `GTIT,[gr],<title>` | P0 | Group title |
| GDEL | — | `GDEL,<g1>,[<g2>],[<inc>]` | P0 | Delete groups |
| GLIST | GLIS | `GLIST,[<g1>],[<g2>],[<inc>]` | P0 | List groups |
| E | — | `E,<ne>,<n1>,…,<n8>` | P0 | Element in the active group (node count by type) |
| EGEN | — | `EGEN,[<itim>],<ninc1>,<e1>,[<e2>],[<inc>],[<ee>]` | P0 | Copy element pattern (K node incremented, attributes copied) |
| EDEL | — | `EDEL,<e1>,[<e2>],[<inc>]` | P0 | Delete elements |
| ECOMPR | ECOM | `ECOMPR` | P0 | Renumber active-group elements 1..N (EOUT remapped) |
| ELIST | ELIS | `ELIST,[<e1>],[<e2>],[<inc>]` | P0 | List elements |
| ETYPE | ETYP* | `ETYPE,<e1>,<e2>,[<inc>],<type>` | P0 | 0 implicit / 1 structure / 2 excavated soil or buried shell |
| EINT | — | `EINT,<e1>,<e2>,[<inc>],<order>` | P0 (TSHELL part P2) | SOLID 0/1/2 integration; TSHELL 0 reduced / 1 selective |
| THICK | THIC | `THICK,<e1>,<e2>,[<inc>],<thick>` | P0 | Shell thickness |
| KI / KJ | — | `KI,<e1>,[<e2>],[<inc>],<k1>,…,<k6>` | P0 | Beam end releases at I / J (1 = released) |
| MSET | — | `MSET,<e1>,[<e2>],[<inc>],<index>` | P0 | Material (M) or, for excavated SOLID/PLANE, soil-layer (L) index |
| MACT | — | `MACT,<index>` | P0 | Active material index for new elements |
| RSET | — | `RSET,<e1>,[<e2>],[<inc>],<index>` | P0 | Property index (R beams, SC springs, MX general) |
| RACT | — | `RACT,<index>` | P0 | Active property index |
| M | — | `M,<nm>,<val1>,<val2>,<weight>,<pdamp>,<sdamp>,<type>` | P0 | Material: type 1 (E, ν), 2 (M, G), 3 (Vp, Vs) |
| L | — | `L,<nm>,<thick>,<weight>,<pveloc>,<sveloc>,<pdamp>,<sdamp>` | P0 | Soil layer |
| R | — | `R,<nm>,<axial>,<shear2>,<shear3>,<tors>,<flex2>,<flex3>` | P0 | Beam section |
| SC | — | `SC,<nm>,<scx>,<scy>,<scz>,<scxx>,<scyy>,<sczz>,<damp>` | P0 | Spring constants (global) and damping ratio |
| MXR / MXI / MXM | — | `MXR,<p>,<row>,<t1>,…,<t12>` | P0 | GENERAL real stiffness / imaginary stiffness / mass rows (upper triangle, row r has 13−r terms) |
| DELM / DELL / DELR / DELSC | DELS* (DELSC) | `DELM,<m1>,[<m2>],[<inc>]` | P0 | Delete materials / layers / beam properties / spring properties |
| MXDEL | MXDE | `MXDEL,<p1>,[<p2>],[<step>]` | P0 | Delete matrix properties |
| MLIST / LLIST / RLIST / SCLIST | MLIS / LLIS / RLIS / SCLI | `MLIST,<m1>,[<m2>],[<step>]` | P0 | List tables |
| MXLIST | MXLI | `MXLIST,<p>` | P0 | List a matrix property |

#### 3.4.E Loads and masses (Manual §9.5)

| Command | Abbr | Syntax | Tier | Meaning |
|---|---|---|---|---|
| F | — | `F,<n>,<fx>,<fy>,<fz>,<tx>,<ty>,<tz>` | P0 | Force factors and arrival times (global X, Y, Z) |
| MM | — | `MM,<n>,<fxx>,<fyy>,<fzz>,<txx>,<tyy>,<tzz>` | P0 | Moment factors and arrival times |
| FDEL / MMDEL | — / MMDE | `FDEL,<n1>,[<n2>],[<inc>]` | P0 | Delete |
| FLIST / MMLIST | FLIS / MMLI | `FLIST,[<n1>],[<n2>],[<inc>]` | P0 | List |
| FSCALE / MSCALE | FSCA / MSCA | `FSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` | P0 | Scale factors (0 → 1; default range = last two defined) |
| MT / MR | — | `MT,<n>,<mx>,<my>,<mz>` / `MR,<n>,<mxx>,<myy>,<mzz>` | P0 | Nodal translational / rotational masses |
| MTGEN / MRGEN | MTGE / MRGE | `MTGEN,[<itim>],[<ninc>],<n1>,<n2>,[<inc>],[<mx>],[<my>],[<mz>]` | P0 | Copy masses with increments |
| MTDEL / MRDEL | MTDE / MRDE | `MTDEL,<n1>,[<n2>],[<inc>]` | P0 | Delete masses |
| MTSCALE / MRSCALE | MTSC / MRSC | `MTSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` | P0 | Scale masses (0 → 1) |
| MTLIST | MTLI | `MTLIST,[<n1>],[<n2>],[<inc>]` | P0 | List masses (both kinds, with units flag) |
| MUNITS | MUNI | `MUNITS,<n1>,[<n2>],[<step>],<units>` | P0 | Per-node mass units 0 mass / 1 weight (÷ g) |

#### 3.4.F Module analysis options

| Command | Abbr | Syntax | Tier | Meaning | Spec |
|---|---|---|---|---|---|
| EQUAKE | EQUA | `EQUAKE,<accopt>,<nrfreq>,<rand>,<damp>,<dur>,<corr>,<seeds>,[tpsd]` | P0 | accopt 0 none / 1 seed record / 2 external (D-EQK-01); nrfreq = records per RSIN file | 07, 05a §5 |
| RSIN / RSOUT | — / RSOU | `RSIN,<no>,<file>` | P0 | Target / output spectrum file for spectrum 1–3 | 07 |
| ACCIN / ACCOUT | ACCI / ACCO | `ACCIN,<no>,<file>` | P0 | Seed-or-external input / generated output acceleration file | 07 |
| TPSD | — | `TPSD,<num>,<file>` | P0 | Target PSD file | 07, 09 §2.28 |
| CORR | — | `CORR,<no>,<time>,<val>` | P1 | Time-varying X–Y correlation pair | 07 |
| EQTIT | EQTI | `EQTIT,<title>` | P0 | Spectra title | 07 |
| SOIL | — | `SOIL,<nrval>,<grav>,<header>,<outcrop>,<save>,<iter>,<ratio>,<gravmult>,<cof>` | P0 | SHAKE options (cof forced to Nyquist) | 07, 05a §6 |
| SPRO | — | `SPRO,<layer>,<prop>,<dynprop>` | P0 | SOIL sublayer → L property + dynamic property label (last SPRO = half-space, D-SOL-02) | 07 |
| DYNP | — | `DYNP,<no>,<sg>,<g>,<sd>,<d>,<label>` | P0 | Curve point: strain % → G/Gmax; strain % → damping % | 07 |
| SACC / SRS / SSTR | — | `SACC,<layer>,<opt>,<outcrop>` / `SRS,<layer>,<save>,<outcrop>` / `SSTR,<layer>,<opt1>,<opt2>,<opt3>,<opt4>` | P0 | SOIL output requests | 07 |
| SSAF | — | `SSAF,<layer>,<save>,<outcrop1>,<outcrop2>,<layer2>,<freqstep>,<title>` | P0 | Spectral amplification factor RS(layer)/RS(layer2) | 07 |
| SFOU | — | `SFOU,<layer>,<out>,<save>,<outcrop>,<smooth>,<nrval>` | P0 (parse only) | "Not usable in this version" | 07 |
| DAMP | — | `DAMP,<d1>,…,<d10>` | P0 | RS damping list (fractions; 0 clears; MOTION uses ≤ 5) | 07 |
| THFILE / THTIT | THFI / THTI | `THFILE,<file>` / `THTIT,<title>` | P0 | Control-motion history file / title (shared) | 07 |
| SITE | — | `SITE,<opmode>,<mode1>,<fstep>,<nl>,<hs>,<mode2>,<wopt>,<freq1>,<freq2>,<cl>,<cm>,<delt>,<nft>,<freq>` | P0 | SITE options; owns Δt, NFFT, Δf, frequency set, control layer | 07, 05a §7 |
| TOPL | — | `TOPL,<l1>,[l2],…,[l200]` | P0 | Free-field layer list (L numbers, top first); 0 clears | 07 |
| WAVE | — | `WAVE,<type>,<opt>,<ratio1>,<ratio2>,<angle>` | P0 | Wave field: type 1 R … 5 L | 07 |
| POINT | POIN | `POINT,<opmode>,<layer>,<rad>` | P0 | Embedment layers (0 surface), central-zone radius | 07, 05a §8 |
| HOUSE | HOUS | `HOUSE,<gravity>,<gelev>,<opmode>,<dim>,<imp>,<coh>,<wpass>,<me>,<cmplxspec>` | P0 (coh/wpass/me P1) | HOUSE options; owns SSI gravity and ground elevation | 07, 05b |
| INCOH | INCO | `INCOH,<gammax>,<gammay>,<gammaz>,<alpha>,<ngp>,<ipr>,<nmodes>,<met>,<HSeed>,<VSeed>,<RandPhz>` | P1 | Incoherency parameters | 07, 05b §3 |
| WPASS | WPAS | `WPASS,<appv>,<ang>,<cohf>` | P1 | Wave passage; `<cohf>` = unlagged coherency model 1–7 (D-INC-01) | 07 |
| ME | — | `ME,<no>,<nfirst>,<nlast>,<xc>,<yc>,<zc>` | P1 | Multiple-excitation zone node range (control point unused) | 07 |
| AMP | — | `AMP,<no>,<a1>,…,<a100>` | P1 | Spectral amplification ratios (append; 0 clears) | 07 |
| SYMM | — | `SYMM,<no>,[<type>],[<node1>],[<node2>],[<node3>]` | P1 (parse P0) | Symmetry (0) / antisymmetry (1) plane; ≤ 2 planes | 07, D-ANL-12 |
| FORCE | FORC | `FORCE,<opmode>` | P0 | FORCE options (rest shared) | 07 |
| ANALYS | ANAL | `ANALYS,<opmode>,<type>,<mode>,<save>,<prnt>,<fopt>,<ang>,<xc>,<yc>,<zc>,<impe>,[simul]` | P0 | ANALYS options | 07, 05c A |
| MOTION | MOTI | `MOTION,<opmode>,<out>,<step>,<dur>,<res>,<freq1>,<freq2>,<fstep>,<mult>,<max>,<rec1>,<rec2>,<fopt>,<bl>,<smo>,<cplx>,<cnvrt>,<pzadj>,<interp>` | P0 | MOTION options | 07, 05c B |
| NOUT | — | `NOUT,<dir>,<code1>,…,<code6>,<node list>` | P0 | Nodal output request (dir 1 x … 6 zz; six 0/1 flags in dialog order) | 07 |
| STRESS | STRE* | `STRESS,<opmode>,<iter>,<save>,<itran>,<interopt>` | P0 | STRESS options | 07, 05d §1 |
| EOUT | — | `EOUT,<code1>,…,<code12>,<group>,<element list>` | P0 | Element output request (0/1/2 per component) | 07 |
| RELD | — | `RELD,<RelDisOutput>,<RelDispSAll>,<RelDispNumFiles>` | P0 | RELDISP options (D-RDP-03) | 07, 05d §2 |
| RELFILE | — | `RELFILE,<FileName>` | P0 | Reference-node complex `.TFI` | 07 |
| RDND | — | `RDND,<NodeNum>,<X>,<Y>,<Z>,<XX>,<YY>,<ZZ>` | P0 | RELDISP node/DOF request (≥ 1 = on) | 07 |

#### 3.4.G Module run commands (Manual §9.12)

| Command | Syntax | Tier | Module |
|---|---|---|---|
| RUNEQUAKE, RUNSOIL, RUNSITE†, RUNPOINT, RUNHOUSE, RUNFORCE, RUNANALYS, RUNCOMBIN, RUNMOTION, RUNRELDISP, RUNSTRESS | `RUN<MOD>,[model]` (default −1 = active) | P0 | the named module; synchronous in batch; output streamed to a UI tab; stdin protocol of three lines |
| RUNNONLINEAR† | `RUNNONLINEAR,[model]` | P2 | NONLINEAR |

† RUNSITE is referenced by Manual §6.4.5 but not listed in §9.12 (D-RUN-01). RUNNONLINEAR is an extension.

#### 3.4.H Model checking (Manual §9.8; read-only except USED)

| Command | Syntax | Tier | Reports |
|---|---|---|---|
| EXCSTRCHK | `EXCSTRCHK` | P0 | Excavation interior nodes shared with structure/beam/spring/GM elements |
| FIXEDINT | `FIXEDINT` | P0 | Interaction nodes with any fixed translation |
| FREESPRING | `FREESPRING` | P0 | Unconstrained nodes connected only to springs |
| HINGED | `HINGED` | P0 | Single-point 6-DOF/3-DOF connections (unintended hinges); beam–shell drilling joints |
| INTCOUNT | `INTCOUNT` | P0 | Number of interaction nodes (and memory/cost estimate, D-ANL-10) |
| KINT | `KINT` | P0 | Beam K nodes that are interaction nodes |
| USED | `USED` | P0 | Fixes all DOFs of unused nodes (modifies model) |

#### 3.4.I Model conditioning and generation (Manual §9.7, §9.9)

| Command | Syntax | Tier | Meaning |
|---|---|---|---|
| FIXSLDROT | `FIXSLDROT` | P0 | D on rotations of solid-only nodes |
| FIXSHLROT | `FIXSHLROT,[stiff]` | P0 | Drilling springs (default 10) at coplanar-shell nodes |
| FIXSPRROT | `FIXSPRROT` | P0 | Fix unstiffened DOFs of spring-only / spring+solid nodes |
| FIXROT | `FIXROT,[Stiff]` | P0 | Combined: D for axis-parallel shells, springs for oblique shells, solids, springs |
| ETYPEGEN | `ETYPEGEN,<type>` | P1 | Set ETYPE for all elements (0 by location → explicit 1/2) |
| INTGEN | `INTGEN,<type>,[level skip]` | P1 | Generate interaction nodes (0 clear, 1 FV, 2 EVBN, 3 FSIN, 4 surface, 5 FFV); additive |
| RADIUS | `RADIUS,<Scale>,<FileName>` | P1 | Excavation element radii r_e = Scale·√A_plan and average (D-PNT-05) |
| GLB2LOC | `GLB2LOC,<Start>,<End>,<Stride>,<Sysno>` | P1 | Convert nodes to a local system |
| GROUPMAT | `GROUPMAT` | P2 | One material per group (Option NON preparation) |
| RMVUNUSED | `RMVUNUSED` | P1 | Remove unused non-interaction nodes |
| NCOM | `NCOM` | P1 | Compress node numbers 1..N, remapping every reference |
| GCOM | `GCOM` | P1 | Compress group numbers |
| WELD | `WELD` | P1 | Merge coincident nodes (keep lowest number; springs protected) |
| MERGE | `MERGE,<Mdl1>,<Mdl2>,<X>,<Y>,<Z>` | P1 | Combine two models into the active model (offset + translate Mdl2) |
| MERGESOIL | `MERGESOIL,<Struct>,<Soil>,[Mode],[Stiff],[Stiff2],[SepLevel],[Mapping]` | P1 | Join structure and excavation models (Mode 0–3) |
| MERGEGROUP | `MERGEGROUP,<dest>,[G1],…,[G10]` | P1 | Merge same-type groups |
| ROTATE | `ROTATE,<x>,<y>,<z>,<rxy>,<ryz>,<rzx>` | P1 | Rotate model about a point (deg) |
| TRANSLATE | `TRANSLATE,<x>,<y>,<z>` | P1 | Translate all nodes |
| EXCAV | `EXCAV,<model>,[delta]` | P1 | Generate excavation volume from the lowest basement grid |
| SOILMESH | `SOILMESH,<dest>,<sX>,<sY>,<hori>,<vert>,<xAdj>,<yAdj>,<Zdepth>,<contact>,<rNum>` | P1 | Near-field soil mesh around the basement |
| CRITFREQ | `CRITFREQ,<tol>,<minfilter>,<TF>,<Var>` | P1 | Frequencies where TFI peaks deviate from TFU → variable |
| FRAMECOMBIN | `FRAMECOMBIN,<op>,<num>,<InFile1>,…,<InFileX>,<Outfile>` | P1 | Combine frame files (0 SRSS, 1 sum, 2 average) |
| FRAMESEL | `FRAMESEL,<tol>,<Acc>,<Var>` | P1 | Critical frame numbers → variable |
| MODFRAMES | `MODFRAMES,<cols>,<framelist>` | P1 | Fix legacy frame headers |

#### 3.4.J Cuts, submodels and section calculations (Manual §9.10, §9.13)

| Command | Syntax | Tier | Meaning |
|---|---|---|---|
| CUTADD / CUTRMV | `CUTADD,<cutnum>,<group>,<e1>,…,<eN>` or `CUTADD,<cutnum>,<group>,RANGE,<start>,[end],[stride]` | P1 | Add / remove elements |
| CUTVOL | `CUTVOL,<cutnum>,[Xmin],[Xmax],[Ymin],[Ymax],[Zmin],[Zmax]` | P1 | Add elements fully inside a box |
| SLICE | `SLICE,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>` | P1 | Add elements crossing a plane |
| CUTCLR | `CUTCLR,<first>,[last],[step]` | P1 | Delete cuts |
| CUT2SUB | `CUT2SUB,<cutnum>,<dest>,[solid]` | P1 | Submodel from a cut (solid ≥ 1: shells → plotting solids) |
| CSECT | `CSECT,<dest>,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>` | P1 | Cross-section model (unit thickness), same group/element numbers |
| EXTRACTEXCAV | `EXTRACTEXCAV,<Model>` | P1 | Submodel of explicit excavation elements |
| SPLITGROUP | `SPLITGROUP,<group>,<split>,[dir]` | P1 | Split a group by a shell plane or axis plane |
| TRANELEM | `TRANELEM,<dest>,<group>,<begin>,<end>,<stride>` | P1 | Transfer elements (overwrite) |
| TRANVOL | `TRANVOL,<dest>,[Xmin],…,[Zmax]` | P1 | Transfer elements in a box |
| READSTR | `READSTR,<Filename>,[Dir]` | P1 | Load `.ess` element stresses onto the active model |
| SECDATAOPT | `SECDATAOPT,<flag>` | P1 | STRESS writes `ESTRESS_n.ess` (1) or not (0) |
| CALCPAR | `CALCPAR,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sysno>,[verbose]` | P1 | Area, centroid, inertias, six resultants |
| CALCMOI | `CALCMOI,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sysno>` | P1 | Section moments of inertia |
| CALCC | `CALCC` | P1 | Volume centroid |
| CALCM | `CALCM` | P1 | Element mass and lumped masses |
| CALCSECTHIST | `CALCSECTHIST,<infile>,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sysno>,[ts],<outfile>` | P1 | Section-force history from `.ess` frames (7-column CSV + MAX row) |
| CALCSECTHISTDB | `CALCSECTHISTDB,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sys>,<ts>,[start],[end],[Stride],<outfile>` | P2 | Same from the binary stress DB |
| SHEAR | `SHEAR,<panel>,[fc],[fy],[P],[Nu],[Fvw],[Fbe]` (args 6–7 read as A_BE, f_y,BE — D-NON-10) | P2 | Wall shear capacities (ACI, Wood, Barda, Gulec–Whittaker) |

#### 3.4.K File conversion (Manual §9.11)

| Command | Syntax | Tier | Meaning |
|---|---|---|---|
| CONVERT | `CONVERT,SSI,<model>,<filename>` / `CONVERT,ANSYS,<model>,<filename>,<gravity>` / `CONVERT,STRUDL,…` | P1 (STRUDL OOS) | `.hou` (+ `.sit`, `.poi`) or `.cdb` → model |
| ANSYS | `ANSYS,[FileName],[Dir]` | P1 | Export APDL (`<model>.inp`) |
| ANSYSREFORMAT | `ANSYSREFORMAT,<Org>,<Map>` | P1 | Regroup beams by release pattern before export |
| ANSYSMODELTYPE | `ANSYSMODELTYPE,<type>` | P2 | Option AA: 1 embedded, 2 surface |
| GENMATRIXDAMP | `GENMATRIXDAMP,<begin>,[end],[stride],<damp>` | P2 | "Not usable in this version" (message only) |

#### 3.4.L Plotting and line mathematics (Manual §9.14)

| Command | Syntax | Tier | Meaning |
|---|---|---|---|
| READSPEC | `READSPEC,<SpecFile>,<numLines>,<Line1>,…,<LineN>` | P1 | Load spectrum-format columns (frequency column not counted) |
| READTH | `READTH,<THFile>,<Pair>,<Num>` | P1 | Load history (0 one-column dt-header, 1 pairs) |
| WRITESPEC | `WRITESPEC,<SpecFile>,<Num1>,…,[Num50]` | P1 | Write ≤ 50 lines on the union grid |
| WRITETH | `WRITETH,<THFile>,<Num>` | P1 | Write a one-column history (constant dt) |
| ADDITION | `ADDITION,<dest>,<source1>,…,<source100>` | P1 | Sum ("Linear Combin.") |
| SUBTRACTION | `SUBTRACTION,<dest>,<source1>,…,<source100>` | P1 | First minus the rest |
| LINECOMBIN | `LINECOMBIN,<Dest>,<Line1>,<Coeff1>,…,[Line100],[Coeff100]` | P1 | Σ c_i y_i |
| AVERAGE | `AVERAGE,<dest>,<source1>,…,<source100>` | P1 | Mean ("Average Line") |
| SRSS | `SRSS,<dest>,<source1>,…,[source100]` | P1 | √Σy² ("SRSS Line") |
| BROADEN | `BROADEN,<Dest>,<Smooth1>,<Smooth2>,<source1>,…,<source100>` | P1 | Envelope, peak broadening (Smooth2 %), peak bridging (Smooth1 %) ("Envelope") |
| LBINCORS | `LBINCORS,<out>,<in>` | P2 | Lower-bound incoherent RS — no algorithm in the manual; "not specified" error |
| SPECPLOT / THPLOT | `SPECPLOT,<Line1>,…,<Line50>` | P1 | 2D spectrum / history plot (−1 ends list; unknown lines ignored) |
| LAYERPLOT | `LAYERPLOT` | P1 | Soil layer column and table |
| SOILPROPPLOT | `SOILPROPPLOT,<PropName>` | P1 | G/Gmax and damping vs strain % |
| MODELPLOT / NODEPLOT | `MODELPLOT` / `NODEPLOT` | P1 | Element / node plot of the active model |
| CUTPLOT | `CUTPLOT,[Cut],[Model]` | P1 | Wireframe + filled cut elements (any model) |
| PROCFRAME | `PROCFRAME,[AniFile],[BufferDIR],[Data],[anitype]` | P1 | Process frame list → frame store + SASSIani.xml |
| BUBBLEPLOT / CONTOURPLOT | `BUBBLEPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>` | P1 | Animated bubble / contour |
| VECTORPLOT / DEFORMPLOT | `VECTORPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>` | P1 | Animated complex-TF vectors / deformed shape |
| AXES | `AXES,<MaxTickX>,<MaxTickY>,<MinTickX>,<MinTickY>,<LogX>,<LogY>` | P1 | Grid and log axes |
| PLOTRANGE | `PLOTRANGE,<Xmin>,<Xmax>,<Ymin>,<Ymax>` | P1 | 2D extent |
| PLOTTITLE / XTITLE / YTITLE / YTITLE2 | `PLOTTITLE,<Title>` | P1 | Titles (soil layer plot ignores PLOTTITLE) |
| LINENAME / MARKERS | `LINENAME,<Num>,<Name>` / `MARKERS,<Mark>,<Ln1>,…,<Ln50>` | P1 | Global line name / markers |
| CAPTUREPLOT / CLOSEPLOT / WINDOWSETTINGS | `CAPTUREPLOT,<FileName>` | P1 | PNG if `.png` in name else BMP / close / settings dialog |
| CNGVIEW / RSTVIEW / CNGCENTER / RSTCENTER | `CNGVIEW,<rX>,<rY>,<rZ>,<px>,<py>,<zoom>` | P1 | 3D view (degrees) and rotation centre |
| ELECOLOR / ELENUM / GROUPNUM / NODENUM | `ELECOLOR,<val>` / `ELENUM,[opt]` | P1 | Colouring (1 group, 2 material, 3 property); labels (ELENUM ⟂ GROUPNUM) |
| NODESEL | `NODESEL,<N1>,…,<N20>` | P1 | Toggle node selection (persistent) |
| SHOWDOF / SHOWMASS | `SHOWDOF,[label1],…,[label6]` / `SHOWMASS,[opt]` | P1 | Fixity / mass markers |
| SHRINK / WIREFRAME / PAUSE | `SHRINK,[switch]` / `WIREFRAME,<Switch>` / `PAUSE,[pz]` | P1 | Element plot modes; animation start/stop |
| COLOR / SHADEROPTIONS / STIPPLE / DEBUG | `COLOR,<Palette>,<Num>,<R>,<G>,<B>` / `SHADEROPTIONS,[points],[linew],[shrink],[scale]` / `STIPPLE,<switch>` / `DEBUG,[switch]` | P2 | Cosmetic |

#### 3.4.M Programming (Manual §9.15)

| Command | Syntax | Tier | Meaning |
|---|---|---|---|
| VAR | `VAR,<Name>,<X1>,…,<Xn>` | P1 | Define variable (counter reset) |
| SETVAR | `SETVAR,…` | P1 | No-op that evaluates counter expressions |
| SHOWVAR / VARLIST | `SHOWVAR,<varname>` / `VARLIST` | P1 | Inspect variables |
| FOREACH | `FOREACH,<Var>,<Command>` | P1 | Loop |
| LOADMACRO / MACRO / MACROLIST | `LOADMACRO,<Name>,<MACROFILE>` / `MACRO,<Name>,<1>,…,<N>` / `MACROLIST` | P1 | Macros |
| LOADVAR | `LOADVAR,<filename>,[name]` | P1 | Variable from file (name = file stem unless given; extension arg) |
| RND / ADDRND | `RND,<var>,<numsamples>,<dist>,<p1>,…` / `ADDRND,…` | P1 | Random lists: UNI, UNIINT (UNINT), NORM, LOGNORM, POISSON (POSSION) |
| RNDSEED | `RNDSEED,<seed>` | P1 | Seed (positive integer) |
| REDUCESET | `REDUCESET,<var>,[sorttype]` | P1 | Sort + de-duplicate (STRING, INT, FLOAT) |

#### 3.4.N Water modelling (Manual §9.16)

| Command | Syntax | Tier | Meaning |
|---|---|---|---|
| FILLPOOL | `FILLPOOL,<Stiff>,<Sensitivity>,<EmptyLevels>,<ShellArea>,<offset>,<stiff2>` | P2 | Water solids + interface springs (defaults 1e6, 0, 0, −1, −1, 0) |
| REFINEMODEL | `REFINEMODEL` | P2 | Split quads 1→4, hexes 1→8 |
| LISTPOOLINTER | `LISTPOOLINTER,<Pool>` | P2 | Pool interface nodes matching the original model |

#### 3.4.O Option NON and nonlinear soil (Manual §9.17) — all P2

| Command | Syntax | Meaning |
|---|---|---|
| EQL | `EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>` | NONLINEAR header (EDF, options bit mask D-NON-02, cut-off %, scale, include elastic damping) |
| P | `P,<num>,<group>,<bbc>,<disp>,<force>` | Wall panel |
| S | `S,<num>,<group>,<elem>,<bbc>,<disp>,<force>` | Nonlinear spring (force 4 GMR) |
| B | `B,<num>,<group>,<spgroup>,<bbc>,<force>,<end1>,<end2>` | Nonlinear beam (not usable; message) |
| BBC | `BBC,<num>,<type>,<points>,<yield>,<file>` | BBC from X–Y file |
| BBCI / BBCP | `BBCI,<num>,<yield>,<type>` / `BBCP,<num>,<point>,<X>,<Y>` | BBC info / single point |
| BBCX / BBCY | `BBCX,<num>,<points>,<yield>,<X1>,…,<Xn>` | BBC vectors |
| BBCGEN | `BBCGEN,<Panel>,<ShearModel>,[fc],[fy],[Pn],[Nu],[bre],[bys],[CrackingForceLevel]` | 22-point BBC from shear capacity |
| DELBBC / DELBM / DELNLS / DELSPR / PDEL / PLIST | `PDEL,<start>,[end],[stride]` | Delete / list nonlinear records (DELNLS deletes NLSLAYER sets; DELSPR deletes S records) |
| PNLGEN (PANELGEN) / PANELIZE / WALLFLR / EDGE / EDGEMODEL / UNIPNL / DGRDFLR / MERGEPANEL | `EDGE,<panel>,[X],[Y],[Z]`; `EDGEMODEL,[x],[y],[z]`; `UNIPNL,<group>`; `DGRDFLR,<scale>`; `MERGEPANEL,<Panel>` | Panel-model construction |
| NONLINMOTDISP | `NONLINMOTDISP` | Add panel corner nodes to MOTION/RELDISP requests |
| NONLINBAT | `NONLINBAT,<Sel>` | Batch script (0 one direction, 1 three with COMB_XYZ_THD) |
| SOLIDPILE | `SOLIDPILE,<group>,[stiff],[soft],[stiff2]` | Pile interface springs (1e7, 10, 1e7; 4 % damping) |
| BEAMPILE / DCOUPLEBEAM | `BEAMPILE` / `DCOUPLEBEAM` | Not applicable / not usable (message) |
| SOILREDEF | `SOILREDEF,<soil>,<dir>` | Excavation layers from a 2D soil model (not validated: warning) |
| NLSOIL | `NLSOIL,<Opt>,<NSTimeSunInc>,<DispConv>,<ForceConv>,<EqualIt>,<BedInt>,<NLDampType>,<MMmult>,<SMmult>` | SOIL-NON global options |
| NLSLAYER | `NLSLAYER,<Num>,[curvefit],[B],[S],[refStrain],[Vis]` | SOIL-NON layer: τ = G0γ/(1+β(γ/γr)^s) + ηγ̇ |

#### 3.4.P Binary databases (Manual §9.18) — P2 except BINOUT parsing (P0)

| Command | Syntax | Meaning |
|---|---|---|
| BINOUT | `BINOUT,[mot],[str],[reldisp]` | Flags (blank = unchanged); reldisp 2 = THD |
| LOADACCDB / LOADDISPDB / LOADTHSDB | `LOADTHSDB,<file>,[sel]` | Load one DB of each type (sel 0 SASSI, 1 ANSYS) |
| DELDB | `DELDB,[sel]` | ALL, ACC, THS, DISP |
| COMBACCDB / COMBDISPDB / COMBDISPDIR / COMBTHSDB | `<CMD>,<Xfile>,<Yfile>,<Zfile>,<Comb>` | Combine directional DBs (algebraic sum; COMBDISPDIR assembles components) |
| BINFRAMEOUT | `BINFRAMEOUT,<db>,<frame>,<TS>,[Split],<dir>` | ASCII frames from a DB (−1 = max frame) |
| BINSTRTBL | `BINSTRTBL,<group>,<EVar>,[step],<file>` | CSV stress table (−1 = signed abs-max) |
| MAXDBFRAME | `MAXDBFRAME,<Type>,[dir]` | Single-frame maximum animation |
| ACCDBANI / DISPDBANI / THSDBANI | `ACCDBANI,<dir>,[label]` | Animations from DBs |

#### 3.4.Q Thick shell (Manual §9.19) — P2

| Command | Syntax | Meaning |
|---|---|---|
| THSHLSTR | `THSHLSTR,<flag>` | 0 eight basic components; 1 adds face stresses/strains |
| THSHLSMH | `THSHLSMH,<passes>,[type],[workdir]` | Transverse-shear smoothing (not validated: warning) |

#### 3.4.R Extension commands of SASSI-EDU (dialect; full names only)

These carry dialog fields that have no argument in the manual's commands, and the auxiliary programs that
the manual runs as separate executables. WRITE emits an X-command only when a value differs from its
default, directly after the base command, so default models produce manual-identical `.pre` files.

| Command | Syntax (defaults) | Tier | Carries |
|---|---|---|---|
| SITEX | `SITEX,<soilmode 0>` | P0 | 0 Linear Soil (L table) / 1 Non-Linear Soil (FILE88 strain-compatible properties) |
| SOILX | `SOILX,<indir 0>,<mult 0>,<max 0.1>,[<cl 0>],[<file>]` | P0 | SOIL Input Direction (0 horizontal, 1 vertical); SOIL's own scaling (exactly one non-zero); optional SOIL-only control layer (0 = SITE `<cl>`) and input history (blank = THFILE) (D-SOL-13) |
| HOUSEX | `HOUSEX,<optimize 0>,<supmode 0>,<nsim 1>,<nlssi 0>,<ansys 0>` | P0 parse (optimize P1, supmode/nsim P1, nlssi P1, ansys P2) | Optimize Model; Superposition Linear/Quadratic; stochastic simulations; Non-Linear SSI (.pin); ANSYS Model Input |
| ANALYSX | `ANALYSX,<ffm 0>,<delrst 0>` | P0 (ffm P1) | FFL (0) / FFM (1); delete restart files after success |
| MOTIONX | `MOTIONX,<f1213 0>,<resp 2>,<srss 0>,<savetf 0>,<saveacc 0>,<savers 0>,<saverot 0>,<rsttf 0>,<rstacc 0>,<rstrs 0>` | P0 (frames P1) | Save FILE12/13 (0/1 FILE13/2 FILE12); vibration response 0 disp / 1 vel / 2 acc; use SRSSTF.txt; post-processing Save/Restart flags |
| STRESSX | `STRESSX,<pzadj 0>,<smo 0>,<skip 0>,<savemax 0>,<saveth 0>,<rstns 0>,<rstsp 0>` | P0 (frames P1) | Phase adjustment, smoothing, output step, post-processing and restart flags |
| RELDX | `RELDX,<saverot 0>,<rstframes 0>` | P1 | Save Rotations for ANSYS; Restart for Frame Generation |
| CMODFORM | `CMODFORM,<form 0>` | P0 | Complex-modulus form 0 SASSI `1−2β²+2iβ√(1−β²)` / 1 `1+2iβ` (benchmarks) |
| EDUOPT | `EDUOPT,<key>,<value>` | P0 | Algorithm switches and tolerances listed in §7 (e.g. `HSLAW`, `NFFTROUND`, `LIMITS`, `BROADENGRID`, `GEOMTOL`) |
| FCOPY / FMOVE | `FCOPY,<src>,<dst>` | P0 | Copy / rename files in the model directory (FILE1 → FILE1X, FILE8 → FILE81 …) so workflows are scriptable |
| REMOVEFREQ | `REMOVEFREQ,<infile>,<outfile>,<n1>,…` | P1 | Remove_Frequencies_from_FILE8 (backup kept) |
| HARMFRAME | `HARMFRAME,<Src>,<Freq>,<OutDir>,[NFrames 24],[Ref]` | P1 | Steady-state harmonic frames `HARM_<ωt°>_<k>` (node, X, Y, Z) of every node at the computed frequency closest to `<Freq>`: u = Re(H e^{iωt}) over one period, from a FILE8-type file or the TFU restart frames; `<Ref>` blank total motion, 0 relative to the free field, n relative to node n (spec 10 §4.5) |
| COMBXYZSTRAIN | `COMBXYZSTRAIN,<FILE74x>,<FILE74y>,<FILE74z>,<out>` | P1 | COMB_XYZ_STRAIN (SRSS of effective strains) |
| BUILDFILE77 | `BUILDFILE77,<out>,<in1>,…,<inN>` | P1 | Build_FILE77 for per-level incoherency |
| COMBXYZTHD | `COMBXYZTHD,<inpfile>` | P2 | COMB_XYZ_THD |
| VERIFY | `VERIFY,[id\|P0\|ALL]` | P0 | Run verification problems of §6 and print pass/fail |
| SHOWSOIL | `SHOWSOIL,[opt -1],[cut],[margin],[depth]` | P1 | The free-field soil drawn around the foundation in the element and node plots (display only; §7.20) |

---

## 4. Computational methodology per module

This section is the implementation contract for the numerics. Equations marked [S]/[V] are taken from R1/R2
and are normative. Where the manual is silent, the [R] choices are fixed by §7 decisions (IDs given). Full
matrix listings that are only referenced here are normative in R1 (section numbers given).

### 4.0 Cross-module conventions

**4.0.1 Time and frequency.**
- Harmonic factor `e^{+iωt}`, ω = 2πf. Horizontal propagation `e^{−ikx}`, outgoing waves have Im k < 0
  (or Re k > 0 when k is real). Cylindrical outgoing waves use `H^{(2)}_μ` (D-CNV-01).
- FFT: `A_k = Σ_{n=0}^{N−1} a_n e^{−2πikn/N}`, `a_n = (1/N) Σ_k A_k e^{+2πikn/N}`; implemented with
  `numpy.fft.rfft/irfft`, k = 0…N/2, `f_k = kΔf` (D-CNV-02).
- `Δf = fstep` if `fstep > 0`, else `Δf = 1/(Δt·NFFT)`. SSI frequencies `f_i = n_i Δf`, n_i positive
  integers, sorted ascending; duplicates are an error; `n_i ≤ NFFT/2` (warning above) [M].
- NFFT must be a power of 2. Otherwise AFWRITE writes the nearest power of 2 with Warning 9, and adds EDU-03
  if that value is shorter than the record (D-CNV-10). Maximum NFFT: 32,768 for EQUAKE/SOIL/SITE/HOUSE/FORCE,
  65,536 for MOTION/RELDISP/STRESS (D-CNV-11).
- Maximum SSI frequencies: 500 per ANALYS run; 1,500 in COMBIN output and MOTION/RELDISP/STRESS input [M].

**4.0.2 Material damping (complex moduli)** [S][V] (D-CNV-03, D-CNV-04).
```
c(β) = 1 − 2β² + 2iβ√(1−β²)            (|c| = 1; V* = V(√(1−β²) + iβ))      default
c(β) = 1 + 2iβ                          selectable by CMODFORM,1 (benchmark only)
G* = G·c(βs)      M* = (λ+2G)·c(βp)      λ* = M* − 2G*                  solids, planes, layers
E*, ν*  from (M*, G*):  ν* = (M*−2G*)/(2(M*−G*)),  E* = 2G*(1+ν*)      beams (E* = E·c(β) when βp = βs)
E* = E·c(β), ν real                     shells (βp must equal βs; CHECK warning otherwise)
k* = k·c(damp)                          springs (SC)
K* = MXR + i·MXI                        GENERAL (as entered)
```
Mass is always real. β ≥ 0.5 is rejected (EDU-04). For ANSYS cross-checks, the equivalent ANSYS structural
damping is `E_ANSYS = E(1−2β²)`, `g = 2β√(1−β²)/(1−2β²)` (R1 §1.2).

**4.0.3 Units and gravity.** ρ = γ/g with g = HOUSE `<gravity>` (SSI modules) or SOIL `<grav>` (SOIL). Masses
with MUNITS = 1 and GM masses with MOPT `<matrix>` = 1 are divided by g. Accelerations of control motions are
in g and converted with the module's g. British/SI is detected from g (g > 20 → British) for unit-dependent
outputs (EQUAKE velocity/displacement/PSD units, SHEAR) (D-CNV-08).

**4.0.4 Geometry tolerance.** `tol = max(1e-6·L_ref, 1e-9)`, with L_ref the bounding-box diagonal of the
active model. Interaction-node elevation matches an interface if `|z − z_iface| ≤ max(tol, 1e-4·h_min)`
(D-GEN-06).

**4.0.5 Transfer-function normalisation.** SITE normalises the free field to unit control motion, so FILE8 holds
`H_k = U_k/U_cp` (dimensionless, the same for total acceleration and total displacement). A response quantity
linear in displacement (stress, strain, relative displacement) is convolved with the **control displacement
spectrum** `U_g(ω) = −g·A(ω)/ω²` (A = FFT of the control acceleration in g; the f = 0 term is set to 0).
Accelerations are convolved with `g·A(ω)` (or with A when output in g) (D-CNV-06).

### 4.1 Element formulations (HOUSE) [M]+[R]

| Element | Stiffness | Mass | Notes / decisions |
|---|---|---|---|
| SOLID | 8-node isoparametric hexahedron, isotropic complex D*; Gauss rule by EINT 0/1/2 → 2×2×2 / 3×3×3 / 4×4×4 (D-ELM-01). Incompatible modes (MOPT `<incomp>` = 0, structural elements only): 9 Wilson–Taylor modes (1−ξ², 1−η², 1−ζ² per direction) with the Taylor centroid-Jacobian correction, statically condensed from K* (D-ELM-02) | M = ½ M_lumped (row sum) + ½ M_consistent | Degenerate prism/pyramid by repeated nodes; positive Jacobian required at all Gauss points (EDU-05) |
| PLANE | 4-node plane strain, unit thickness, X–Z plane, 2×2 Gauss; triangles by repeated node; same incompatible-mode option for structural elements (D-ELM-03) | ½ lumped + ½ consistent | node order normalised to positive Jacobian with a warning |
| BEAMS | 2-node 3D Timoshenko (12×12), axis 1 = I→J, axis 2 toward K, axis 3 = 1×2; bending in 1–2 plane uses E·I3 and As2, in 1–3 plane E·I2 and As3; `φ = 12EI/(G·As·L²)`, As = 0 → φ = 0; releases by static condensation of the local matrices (D-ELM-04) | consistent (translation ρAL; torsion ρ(I2+I3)L/3 coupling as standard; no rotary inertia) (D-ELM-05) | K node carries no DOF; K collinear with I–J → Error 9 |
| SHELL | Flat facet: membrane = plane-stress Q4 with incompatible modes (CST for triangles); bending = DKQ (quad) / DKT (triangle) Kirchhoff plate; **zero drilling stiffness** (D-ELM-06) | lumped: ρtA/n per translational DOF, zero rotational inertia (D-ELM-07) | warped quads projected on the mean plane (W3/E12 thresholds D-CHK-05) |
| TSHELL (P2) | MITC4 Mindlin–Reissner; EINT 0 reduced (1-point bending and shear, hourglass-stabilised), 1 selective (2×2 bending, 1-point shear); drilling stiffness `1e-4 × min membrane diagonal × element area` (units of a rotational stiffness; D-ELM-08) | lumped | THSHLSTR face stresses |
| SPRING | six uncoupled global constants k_d·c(damp) between I and J | none | zero length allowed (Error 8 excluded) |
| GENERAL | 12×12 `K_R + iK_I`, mass M (÷ g if weight units); 3 nodes → local axes as BEAMS, transformed `T = blockdiag(Λ,Λ,Λ,Λ)` | user | each element couples two nodes only |
| Nodal masses | MT/MR on diagonal (÷ g if MUNITS = 1); ignored on fixed DOFs (W5/W6) | | default MUNITS = 1 (D-MDL-08) |

**DOF management.** Active DOFs of a node = union of the DOFs of attached element types, minus D-fixed DOFs.
DOFs no element defines (e.g. rotations of SOLID-only nodes) are eliminated automatically with an
information message; DOFs defined but unstiffened (SHELL drilling) are kept and need FIXROT/FIXSHLROT
(warning EDU-06 otherwise) (D-ELM-09). Gap and unused nodes are fixed in the analysis files (W1/W4).
Interaction DOFs are the translations of interaction nodes only [M].

### 4.2 SITE — layered free field (thin-layer method) [S: R1 §2]

**Inputs.** TOPL list of L layers (top first), half-space layer `<hs>`, `<nl>` generated sublayers, frequency
set, Δf, control layer `<cl>` and direction `<cm>`, wave family `<wopt>`, WAVE records, Frequency 1/2
(frequency numbers, D-SIT-05), SITEX soil mode.

**Layer matrices.** For each sublayer, the Kausel-form TLM matrices A_m, B_m, G_m and the mixed mass
`M_m = ½M_cons + ½M_lump` (R1 §2.2, normative; DOF order `{u_x, u_y, ũ_z = i·u_z}` per interface, z up).
Assembly by overlapping interface blocks. Interface relation `P̄ = (A k² + B k + G − ω²M) Ū`.

**Half-space (variable-depth + viscous boundary).** For `nl > 0` add `nl` sublayers with half-space properties,
total thickness `H_hs(f) = 1.5·Vs_hs/f`, thickness non-decreasing with depth, every sublayer `≤ λ_s/8`
(law D-SIT-02), and add `+iω·ρ_hs·V*` dashpots at the bottom interface (`c_s` on x, y; `c_p` on z, using
complex V* of the half-space). `nl = 0` → rigid base (bottom DOFs removed) [M]+[R].

**Mode 1 (FILE2).** Per frequency:
- Rayleigh: linearised generalised eigenproblem in k² of size 2N_f:
  `(k² [A_x 0; B_xzᵀ A_z] + [C_x B_xz; 0 C_z]) {φ_x; kφ_z} = 0`, `C = G − ω²M` (+ dashpots).
  Root `k_j = √(k_j²)` with Im k_j < 0 (Re k_j > 0 if real). Normalise
  `φ_xᵀA_xφ_x + φ_zᵀA_zφ_z + k_j⁻¹φ_zᵀB_xzᵀφ_x = 1`.
- Love: `(A_y k² + C_y)φ_y = 0`, N_f modes, `φ_yᵀA_yφ_y = 1`.
- FILE2 stores: generated sublayer thicknesses/properties, {k_j^R, φ_x, φ_z} and {k_l^L, φ_y} per frequency.
  One FILE2 serves 2D and 3D.

**Mode 2 (FILE1).** Unit-normalised free field at every interface for each frequency:
- Vertically incident SV/SH (shear column, `G_s` with G*) and P (`G_p` with M*):
  `(G − ω²M + iωc e_N e_Nᵀ) u = 2iωc E_b e_N` (dashpot base) or `u_N` prescribed (rigid base);
  normalise `Û(z) = u(z)/u_d(z_cp)` at the **top of the control layer**, control motion taken as the
  **within** (in-column) motion (D-SIT-07).
- Inclined body waves (P1): `k = ω sinθ/V_hs` (V = Vs for SV/SH, Vp for P), base load from the incident plane
  wave, exact Kausel–Roesset half-space stiffness for θ > 0 when `nl > 0` (D-SIT-08); field
  `Û(z) e^{−ik(x′−x′_cp)}`.
- Surface waves (P1): mode j selected by `<opt>` (1 shortest wavelength = largest Re k among propagating modes;
  2 least decay = smallest |Im k|; Love uses the shortest-wavelength rule) (D-SIT-09); field
  `φ_j(z) e^{−ik_j(x′−x′_cp)}/φ_{j,d}(z_cp)`, physical `u_z = −iφ_z`.
- Wave mix: `U = Σ_w r_w(n) Û_w`, ratios interpolated linearly between Frequency 1 and 2 (frequency
  numbers), held constant outside, each in (0,1], sum = 1 (Errors 53–55) [M].
- SITEX soil mode 1: SOIL FILE88 strain-compatible Vs, βs (and Vp, βp per D-SOL-06) replace the TOPL layer
  properties position by position; the layer counts must match (EDU-07).
- FILE1 stores per frequency: Û at every interface (3 complex components), wave numbers, wave types, ratios,
  control point, `<cm>`.

**Rules.** Interaction nodes must lie on interfaces (Error EDU-01 if not, in HOUSE/ANALYS); control point at an
interface [M]; 1/5-wavelength rule `h ≤ Vs/(5 f_cut)` warning G-05; ν > 0.47 warning G-09; ≤ 20 layers
warning; FILE1X/Y/Z must be produced with angle 0 (SV x′, SH y′, P z′) for simultaneous cases [M].

### 4.3 POINT — point-load solutions (POINT3 P0, POINT2 P1) [S: R1 §3]

- **Central zone (3D):** one radial axisymmetric element of radius R0 (`<rad>`), vertical discretisation =
  the SITE sublayers (including generated half-space sublayers and base dashpots), same mixed mass. Fourier
  harmonics μ = 0 (vertical load) and μ = 1 (horizontal load). Axis constraints: μ = 0 → ũ_ρ = ũ_θ = 0;
  μ = 1 → ũ_ρ = ũ_θ tied, ũ_z = 0.
- **Exterior:** `ũ(ρ) = Ψ_μ(ρ) α` with outgoing-mode matrix Ψ_μ built from FILE2 modes and `H^{(2)}_μ`;
  boundary stiffness `R_μ = −𝒯_μ Ψ_μ(R0)⁻¹` (R1 §3.3 (ii)–(iii), normative).
- **Solve** for unit loads at interfaces n = 1 … L+1 (L = `<layer>`, 0 = surface only):
  `[K_c − ω²M_c + R_μ] Ũ = F̃` (F̃ = P/π for μ = 1, P/(2π) for μ = 0); then `α = Ψ_μ(R0)⁻¹ Ũ_boundary`.
- **FILE3** stores per frequency, load interface and μ: α (3N_f complex), axis displacements at interfaces
  1…L+1, and R0. ANALYS evaluates `Ψ_μ(ρ)α` exactly at any distance (D-PNT-02).
- **POINT2:** strip |x| ≤ R0 (2 elements) with Waas–Lysmer consistent boundaries at ±R0 (R1 §2.6, §3.5); far
  field `Σ α_j φ_j e^{−ik_j(|x|−R0)}`, vertical response to horizontal load odd in x.
- POINT2 vs POINT3 is selected from HOUSE `<dim>` (1 → POINT2, 2 → POINT3) (D-PNT-01).

### 4.4 HOUSE — model matrices, classification, incoherency, nonlinear soil

1. **Classification** of SOLID/PLANE: ETYPE 1 structure (M table), 2 excavated soil (L table via MSET), 0 →
   excavated if every node `z ≤ gelev + tol` and the centroid is strictly below gelev, else structure
   (D-HOU-01). SHELL/TSHELL ETYPE 2 = buried shell (structure, nodes are interaction nodes). Incompatible
   modes never on excavated elements.
2. **Assembly:** frequency-independent complex sparse matrices `K*_s, M_s` (structure, incl. near field) and
   `K*_e, M_e` (excavated soil), on one global DOF map. Excavated elements use the L layer of their MSET
   index; EDU-08 warns if that layer is not the TOPL layer at the element mid-depth.
3. **Interaction set:** INT/INTGEN flags; translational DOFs only. Rules enforced: interface elevations
   (EDU-01, error), at or below grade (error), bottom-up ascending numbering (error when incoherent, warning
   otherwise), no structural node is an interaction node unless shared with the excavation boundary (G-14),
   EXCSTRCHK condition (error).
4. **Outputs:** FILE4 (`<model>.N4`: topology, DOF map, element recovery data), COOSK/COOSM, DOFSMAP,
   FILE90 (hashes of interaction-node coordinates/numbering, layer table, FILE4), FILE91 (metadata), listing
   with the excavation-element layer assignment table [M].
5. **Optimize Model (P1):** reverse Cuthill–McKee renumbering with interaction nodes kept in their original
   relative (bottom-up) order; writes `.hounew` and `.map`; FILE4/FILE8 in the new numbering (D-HOU-03).
   The internal sparse solver uses its own fill-reducing ordering regardless.
6. **Incoherency (P1)** per SSI frequency and direction c ∈ {X, Y, Z} (D-INC-*):
   - horizontal projections of the N interaction nodes; distances in the Line-D frame (angle `WPASS <ang>`);
     `D_ij = √(2(αΔX′² + (1−α)ΔY′²))` for models 2–7; plain horizontal distance for model 1.
   - Σ_ij = γ_model,c(f, D_ij) (unit diagonal, symmetric); Luco–Wong `γ = exp[−(γ_c ω D / Vs)²]` with Vs = INCOH
     `<alpha>`; Abrahamson models from coefficient data files (D-INC-04); user tables bilinear in (f, D),
     clamped, γ(D = 0) = 1.
   - `Σ = ΦΛΦᵀ`, λ descending, negative round-off clipped to 0; check `|Σλ − N|/N < 1e-8` (Eq. 6.1); print the
     `I N C O` table (υ_j = 100λ_j/N, cumulative) if `<ipr>` = 1.
   - Synthesis of factors s (length N): stochastic `s = Φ_m Λ_m^{½} e^{iθ}`, θ ~ U[−RandPhz, +RandPhz],
     independent per frequency, mode, direction and sample (PCG64 streams seeded from HSeed for X and Y, VSeed
     for Z); deterministic AS `s = Σ_k √λ_k φ_k` with sign convention `Σ_i φ_ik ≥ 0`; single mode
     (`nmodes = −k`) `s = √λ_k φ_k`.
   - Wave passage `s_i ← s_i·exp(−iωτ_i)`, `τ_i = ((x_i−x_c)cosθ_D + (y_i−y_c)sinθ_D)/V_app`, reference =
     ANALYS control point (D-INC-06). Multiple excitation `s_i ← s_i·SAR_zone(i)(ω)` (SAR = 1 outside zones).
   - FILE77 (or FILE77sss) stores s per frequency, direction and interaction node.
   - Allowed only for 3D, no SYMM, and models 2–7 / ME require wave passage (errors).
7. **Nonlinear soil (P1):** `.pin` defines nonlinear groups, ESF, ISTR, and per-material GFAC, DFAC, ICURVE.
   First run: G = GFAC·G_layer, β = DFAC·β_layer; writes `.liq` = 1. Later runs (`.liq` = 1): read FILE74
   (combined effective strains), set `G = G_max·(G/G_max)(γ_eff)`, `β = D(γ_eff)` from FILE73 curve ICURVE,
   one internal material per element; write FILE78 (properties used) (D-NLS-*).

### 4.5 FORCE — load vectors [M]+[R]

For each SSI frequency and loaded DOF k (F: global X/Y/Z forces, MM: moments) with factor a_k and arrival time
t_k: `P_k(ω_j) = a_k·exp(−iω_j t_k)`. Repeated definitions follow MOPT `<force>` (overwrite or add; in add mode
factors add and the latest arrival time is kept with a warning, D-FRC-02). Loads on fixed DOFs are ignored
(W10/W11). FILE9 stores P per frequency on the FILE4 DOF map. The response is
`IFFT[H(ω)·F_ref(ω)]` with F_ref the FFT of the reference load history (MOTION THFILE) [R].

### 4.6 ANALYS — impedance and SSI solution [M]+[S: R1 §4]

Per SSI frequency (order number q, ω = 2πf_q):

1. **Frequency survey:** match integer frequency numbers in the requested set (or all of FILE1/FILE9 when
   `<fopt>` = 1) against FILE1/FILE9 and FILE3; stop on any miss [M].
2. **Flexibility F_ff** (3n_f × 3n_f) from FILE3 [S]: for interaction node i on interface m and j on
   interface n, r = horizontal distance, θ = atan2(y_i−y_j, x_i−x_j), with POINT values for load at n observed
   at m (μ = 1: ũ, ṽ, w̃; μ = 0: p̃, q̃):
   ```
   F_ij = [ ũc²+ṽs²     (ũ−ṽ)sc    p̃c ]
          [ (ũ−ṽ)sc     ũs²+ṽc²    p̃s ]      c = cosθ, s = sinθ
          [ w̃c          w̃s         q̃  ]
   ```
   r = 0: `diag(ũ_ax, ũ_ax, q̃_ax)`; 0 < r < R0: linear interpolation in r between axis and R0 values;
   symmetrise `F ← ½(F + Fᵀ)` (D-ANL-02). 2D: the 2×2 (P-SV) or 1×1 (SH) block with odd terms × sgn(x_i−x_j).
3. **Impedance** `X_ff = F_ff⁻¹` (dense complex LU; complex-symmetric, never Cholesky/Hermitian).
4. **System** (general assembly, which reduces to Eq. 2.1):
   `𝐂(ω) = 𝒜_sᵀ(K*_s − ω²M_s)𝒜_s − 𝒜_eᵀ(K*_e − ω²M_e)𝒜_e + 𝒜_fᵀ X_ff 𝒜_f`.
5. **Load:** seismic `𝐛 = 𝒜_fᵀ X_ff U′_f` with U′_f from FILE1 at each interaction node's interface, rotated
   `[u_x,u_y] = R(a)[u_x′,u_y′]`, oblique phase `e^{−ik x′_j}` with `x′_j = (x_j−x_c)cos a + (y_j−y_c)sin a`,
   incoherency (FFL: `P ← s ⊙ P`; FFM: `U′ ← s ⊙ U′` before X_ff; the factor s_i of the control-motion
   direction multiplies all three components at node i) and SAR; vibration `𝐛` = FILE9 vector.
6. **Solve** by Schur complement (D-ANL-03): order DOFs as non-interaction n and interaction f; sparse LU of
   `𝐂_nn`; `S = 𝐂_ff − 𝐂_fn𝐂_nn⁻¹𝐂_nf` (where 𝐂_ff excludes X_ff); dense LU of `(S + X_ff)`;
   `U_f = (S+X_ff)⁻¹(𝐛_f − 𝐂_fn𝐂_nn⁻¹𝐛_n)`, `U_n = 𝐂_nn⁻¹(𝐛_n − 𝐂_nf U_f)`. All simultaneous right-hand sides
   share the factorisations. Fallback: if 𝐂_nn is singular or ill-conditioned at a frequency (e.g. undamped
   structure with the interaction DOFs clamped; estimated rcond < 1e-12), solve the full system by sparse LU
   with the dense X_ff block inserted, and warn; the restart record then stores that factorisation.
7. **Output:** H for every DOF (fixed DOFs = 0) to FILE8 (or FILE8X/Y/Z, FILE8nnn); listing amplitudes or
   Re/Im per `<prnt>`.
8. **Restart files** (`<save>` = 1): COOXqqq = X_ff; COOTKqqq = sparse LU factors of 𝐂_nn (L, U, row/column
   permutations) + dense LU of (S+X_ff); every record stores f and the FILE90 hashes; COOXI/COOTKI index the
   records by frequency number (lookup by value, so subsets work) (D-ANL-05). Mode 2 reuses COOX; Mode 3
   reuses COOX and COOTK.
9. **Global impedance** (`<impe>` 1/2, 3D only): `K_G(ω) = Tᵀ X_ff T` over translational interaction DOFs,
   `T_j = [[1,0,0,0,dz,−dy],[0,1,0,−dz,0,dx],[0,0,1,dy,−dx,0]]`, d = r_j − (x_c, y_c, z_c);
   FOUNSTIF = Re K, FOUNDASH = Im K/ω, FOUNDAMP = Im K/(2|Re K|), FOUNIMPD = |K| (D-ANL-07). Warning:
   unconstrained impedance; equals the rigid-foundation impedance only for surface foundations.
10. **Resource check:** memory ≈ (3N_int)²·16 B per dense matrix (× ≈ 3 working copies) compared with available
    RAM before the run (G-22, D-ANL-10).

Physical checks built in: C^s = C^e on the excavation (FV) gives U = U′ (VP-16); SM is near-singular at the
excavated-volume frequencies (expected; VP-19).

### 4.7 COMBIN [M]+[R]

FILE81 ∪ FILE82 by frequency number, sorted. DOF maps and Δf must be identical (error). A frequency present in
both is an error unless `EDUOPT,COMBINDUP,PREFER82` (D-CMB-01). Output may exceed 500 frequencies (≤ 1,500).

### 4.8 MOTION — interpolation, convolution, spectra [M]+[S: R1 §5]+[R]

1. **Control motion:** read THFILE (`<fopt>` 0: first value dt then one value per record; 1: (t, a) pairs);
   records rec1…rec2 (1-based samples after the header); scale by `<mult>` or to `<max>` (exactly one non-zero,
   Errors 77/78); error if record length > NFFT; zero-pad; `A = rfft(a)`; DF must equal FILE8 Δf (1e-6 rel.).
2. **Interior interpolation** of H from the SSI frequencies f_1 < … < f_N to f_k ≤ f_N:
   - **Options 0–5:** local 2-DOF hysteretic rational form in ω² (Tajirian) on 5-point windows
     `H(ω) = (C1ω⁴ + C2ω² + C3)/(ω⁴ + C4ω² + C5)` with the five complex constants from
     `[ω_p⁴, ω_p², 1, −ω_p²H_p, −H_p]·C = ω_p⁴H_p`, p = 1…5, ω scaled by the window maximum, solved by
     SVD/least squares (rcond 1e-10) to handle rank-deficient (SDOF-like) windows. Guard: if
     `min|ω⁴+C4ω²+C5| < 1e-3·max|·|` on the window span, use the complex cubic through the 4 nearest points
     (D-MOT-01; supersedes spec 05c B.6.4).
   - **Window schemes** (start index m; window = points m…m+4; candidate starts for interval j:
     `max(1, j−3) ≤ m ≤ min(j, N−4)`):

     | Option | Windows | Combination |
     |---|---|---|
     | 1 (SASSI 1982) | starts 1, 5, 9, … sharing end points; last window = last 5 points | single |
     | 4 | starts 2, 6, 10, …; leading interval uses m = 1 | single |
     | 5 | starts 3, 7, 11, …; leading intervals use m = 1 | single |
     | 2 | all candidate windows | arithmetic mean |
     | 0 (SASSI2000) | all candidate windows | weighted mean, `w_m = max(1e-3, 1 − \|f − c_m\|/h_m)` |
     | 3 | 3 candidates with the smallest `\|f − c_m\|` | arithmetic mean |

     Every scheme returns H_j exactly at computed frequencies. N < 5 → option 6; N < 3 → complex linear.
   - **Option 6:** independent not-a-knot cubic splines of Re H and Im H through the computed points plus the
     f = 0 anchor (D-MOT-02).
3. **Outside the computed range:** `f > f_N` → H = 0. Seismic `0 ≤ f < f_1` → complex linear between the
   anchor H(0) and H_1, H(0) = projection of the control direction on the output DOF (1 parallel translation,
   cos a / sin a for rotated horizontals, 0 otherwise and for rotations). Vibration → H = H_1. Nyquist bin real
   (D-MOT-03).
4. **Smoothing** (options 0–5, S > 0): for f in [f_j, f_{j+1}], linear reference L(f),
   `r = max(0, |H|/A_max − 1) + max(0, 1 − |H|/max(A_min, 1e-12A_max))`, `H_s = L + (H − L)/(1 + S·r)`
   (D-MOT-04). S = 0 → identity. Ignored with a warning for option 6.
5. **Phase adjustment** (`<pzadj>` = 1): `φ = unwrap(arg H)` from f = 0; `H′ = |H|·e^{iρφ}`, ρ = 1/(1+S) for
   options 0–5, ρ = 0 for option 6 (D-MOT-05).
6. **SRSS TF** (MOTIONX srss = 1 with SRSSTF.txt): interpolate each modal FILE8 TF, `|H| = √Σ|H_k|²`, phase 0
   (p = 0) or coherent phase (p = 1) (D-MOT-06).
7. **Convolution:** `R_k = H_k·A_k`; `r = irfft(R, NFFT)`; output length `min(NFFT, round(1.2·dur/Δt))`
   (dur = 0 → NFFT). Vibration: displacement `H·F`, velocity `iωHF`, acceleration `−ω²HF` (MOTIONX resp).
8. **Response spectra:** exact piecewise-linear Nigam–Jennings recursion (R2 G.1 coefficients, normative) over
   the full NFFT record; frequencies log-spaced from f1 to f2 inclusive with `<fstep>` points (D-MOT-07); SA =
   max|absolute acceleration| (`−(2ζωẋ + ω²x)`) written to `.RS`; SV = max|ẋ| (relative) in the listing; for
   f > 0.1/Δt the record is FFT-upsampled ×4 before integration. Damping list from DAMP (≤ 5 for MOTION).
9. **Baseline correction** (`<bl>` = 1, D-MOT-08): trapezoidal integration; least-squares fit
   `v ≈ c1t + c2t²/2`; subtract `c1 + c2t` from the acceleration; re-integrate; report final displacement;
   FILE13 (`t acc vel disp`) only when baseline correction is on and MOTIONX f1213 = 1.
10. **Outputs:** `.TFU` (computed f, |H| [, phase rad]), `.TFI` (Fourier grid 0…f_N), `.ACC`, `.RS`, FILE12/13,
    frames, `ACC_max.txt`, CONTTRS `.RSO` (formats D-FIL-02). Duplicate NOUT nodes are merged (flags OR'ed) at
    AFWRITE with a warning (D-MOT-10).

### 4.9 RELDISP [M]+[R]

`H_rel(f) = H_node,d(f) − H_ref,d(f)` on the Fourier grid, from complex `.TFI` files (RELFILE = reference;
free-field reference = synthetic TFI with amplitude 1, phase 0 in the input direction). Seismic:
`D(f) = H_rel(f)·(−g·A(f)/ω²)` (f = 0 → 0), `d(t) = irfft(D)`. Vibration: `D = H_rel·F_ref`. `.TFD` stores the
complex relative-displacement TF per unit control acceleration (length per g) (D-RDP-02); `.THD` the history;
the listing the maxima. One DOF per run (the reference TFI's DOF tag). Frames and all-node output per RELD and
RELDX flags.

### 4.10 STRESS [M]+[R]

1. **STF at SSI frequencies:** element nodal TFs from FILE8, `STF_c(f_j) = S_c·U_e(f_j)`:
   SOLID/PLANE `σ = D*B u` and `ε = B u` at the centroid in global axes (incompatible modes recovered from the
   condensation); SHELL membrane stress `D*_m B_m u` and moments per unit length `D*_b κ` in local x′y′;
   TSHELL N, Q, M (P2); BEAMS local end forces `f_L = k_L T u` (forces exerted on the element at I and J, local
   axes 1/2/3, D-STR-05); SPRING `f_d = k*_d (u_J,d − u_I,d)` (global).
2. **Interpolate** each STF with the same option, smoothing and phase-adjustment logic as MOTION (STRESS
   `<interopt>`, STRESSX).
3. **Convolve** with the control displacement spectrum `U_g = −g·A/ω²` (seismic) or F_ref (vibration); inverse
   FFT. Derived quantities (octahedral shear stress `τ_oct = ⅓√[(σxx−σyy)² + (σyy−σzz)² + (σzz−σxx)² +
   6(τxy²+τyz²+τxz²)]`, shear strains) are computed in the time domain; the six components are computed
   whenever octahedral output is requested.
4. **Outputs:** maxima (listing), `.THS`, `.TFU/.TFI` (`<itran>` = 1, all element types, D-STR-06), FILE14/15,
   `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT` (all-element option), `.sig/.tau/.bdsig/.bdtau`, frames (`Frames.txt`
   frame numbers = 1-based time-step indices, D-STR-08), nodal averages = plain mean of adjacent element
   centres (plotting only), soil pressure = −σ_n on the solid face shared with the structure (compression
   positive) (D-STR-09); `STATIC_SOIL_PRESSURES.TXT` added algebraically; `ESTRESS_n.ess` (SECDATAOPT = 1).
5. **Nonlinear soil (P1):** for nonlinear groups, γ(t) per ISTR (0 max of γxy, γxz, γyz; 1 octahedral
   `(2/3)√[(εx−εy)² + (εy−εz)² + (εz−εx)² + 1.5(γxy²+γyz²+γxz²)]` in 3D, `√[(εx−εz)² + γxz²]` in 2D);
   `γ_eff = ESF·max|γ(t)|`; write FILE74 (element, γ_eff, G, β from FILE73).

### 4.11 SOIL — SHAKE equivalent-linear site response [S: R1 §6, SHAKE91]

- Column = SPRO sublayers 1…N; the **last SPRO entry is the half-space** (D-SOL-02); ρ = γ/g_SOIL,
  G_max = ρVs².
- Recursion (e^{iωt}, z_m down from the layer top): `u_m = E_m e^{ik*z} + F_m e^{−ik*z}`,
  ```
  E_{m+1} = ½E_m(1+α_m)e^{ik_m h_m} + ½F_m(1−α_m)e^{−ik_m h_m}
  F_{m+1} = ½E_m(1−α_m)e^{ik_m h_m} + ½F_m(1+α_m)e^{−ik_m h_m}
  α_m = ρ_m V*_m /(ρ_{m+1} V*_{m+1}),   E_1 = F_1
  ```
  within motion at a layer top `E + F`, outcrop `2E`. Input at the top of the control layer `<cl>` as outcrop
  (`<outcrop>` = 1) or within; TFs scaled so the control motion equals the input.
- Strain at mid-sublayer `γ(ω) = (E_a e^{ik*h/2} − F_a e^{−ik*h/2})/(iωV*)` with acceleration amplitudes;
  `γ_max = max_t|γ(t)|`; `γ_eff = R_γ·γ_max` (`<ratio>`).
- Update `G = G_max·(G/G_max)(γ_eff)`, `β = D(γ_eff)` by **linear interpolation in log10(γ)** with constant
  extrapolation; DYNP strain and damping in percent (D-SOL-03). Run exactly `<iter>` iterations and report
  ΔG %, Δβ % per iteration (D-SOL-04).
- Input direction 1 (SOILX): Vp, βp, constrained modulus, no iterations [M].
- Frequencies above the cut-off are removed: `<cof>` is fixed to Nyquist [M]; `EDUOPT,SOILCUTOFF,<Hz>`
  reproduces SHAKE91 runs (25 Hz in the sample problem) (D-SOL-08).
- Outputs: `ACCxxx.TH` (top of layer, within or outcrop per SACC), `SNxxx.TH`/`SSxxx.TH` (mid-layer strain %,
  stress), RS per DAMP at layer tops (SA in g × `<gravmult>`, D-SOL-07), SSAF ratios, FILE73 curves, FILE88
  (per sublayer: thickness, γ_eff, G, Vs, βs, Vp, βp) with Vp/βp per D-SOL-06.

### 4.12 EQUAKE — spectrum-compatible motions [M]+[R]

1. **Target:** RSIN two-column spectrum (Hz, SA in g) at damping `<damp>`, log-log interpolated; check record
   count = `<nrfreq>` (Error 89).
2. **Initial motion** (D-EQK-02): random phases (PCG64 seeded by `<rand>`) with a stationary-part envelope
   (Saragoni–Hart, strong-motion part ≥ 6 s within `<dur>`), or seed-record Fourier phases (`<accopt>` = 1);
   Fourier amplitudes from the target by the SIMQKE (Gasparini–Vanmarcke) PSD estimate.
3. **LW frequency-domain matching (P0):** iterate `|A(f)| ← |A(f)|·SA_target(f)/SA_comp(f)` (5 % RS by
   Nigam–Jennings at 100 points/decade), phases kept, envelope re-applied, up to 20 iterations or until the
   SRP band is met.
4. **AB time-domain matching (P1):** RspMatch-type tapered-cosine wavelet adjustments (Al Atik & Abrahamson
   2010 form) to remove residual misfit without drift.
5. **Baseline correction:** complex-frequency correction (zero mean velocity: remove f = 0 content) followed by
   a time-domain polynomial (order ≤ 3) fit to the displacement, subtracting its second derivative, so final
   velocity and displacement ≈ 0 (D-EQK-03).
6. **Checks (warnings in output and screen):** duration ≥ 20 s; Δt ≤ 0.005 s (manual) and Nyquist ≥ 50 Hz (SRP
   Rev 4); RS at ≥ 100 log points per decade over the target range: never > 10 % below, never > 30 % above, ≤ 9
   adjacent points below; strong-motion (Arias 5–75 %) ≥ 6 s; |ρ| ≤ 0.16 between components; PSD ≥ 80 % of
   TPSD target over 0.3–24 Hz when given (D-EQK-04).
7. **Outputs:** `.acc` (g), `.vel`, `.dis` (in/s, in for British; cm/s, cm for SI), `.rso`, `.psd`
   (one-sided `S0 = 2|F|²/(2πT_D)` of the strong-motion window, ±20 % band average; in²/s³ or cm²/s³),
   `.fft` (Hz, Re, Im; strong-motion window, f ≥ 0), listing with PGA/PGV/PGD, V/A, AD/V², strong-motion
   duration, stationary and 2-s moving-window correlations.
8. **Seeds:** `<seeds>` trials (0 → 1) with seeds rand, rand+1, …; keep the trial with the smallest maximum
   deviation from target (D-EQK-05).
9. **Correlation (P1):** target ρ(t) piecewise linear from CORR pairs; in each 2-s window mix normalised
   components `y′ = ρx + √(1−ρ²)y`, then rerun LW/AB with the result as seed (manual warning) (D-EQK-06).

### 4.13 Line mathematics and spectrum broadening (P1) [M]+[R]

- All multi-line operations work on the **union of abscissas**, with linear interpolation inside each line and
  constant extrapolation outside (`numpy.interp(X, xs, ys, left=ys[0], right=ys[-1])`) [M].
- ADDITION Σy; SUBTRACTION y1 − Σy_{i≥2}; LINECOMBIN Σc_i y_i; AVERAGE mean; SRSS √Σy² [M]. Undefined source
  lines are errors (plots ignore them) (D-LIN-02).
- BROADEN,Dest,Smooth1,Smooth2 (D-LIN-01): (1) envelope E = max_i y_i; (2) peak broadening with b = Smooth2/100:
  `B(f) = max{E(f′) : f/(1+b) ≤ f′ ≤ f/(1−b)}` on the grid augmented with f(1±b) points; (3) peak bridging
  with p = Smooth1/100: raise valleys between adjacent peaks to min(P_k, P_{k+1}) when
  `(min(P_k,P_{k+1}) − valley)/min(P_k,P_{k+1}) ≤ p`, repeat until stable; name "Envelope". Warning if a source
  has < 301 points.

### 4.14 Section-cut integration (P1) [R: spec 10 §3.2, normative]

Local axes `ez = n/|n|`, `ex = (r − (r·ez)ez)/|·|`, `ey = ez × ex`, origin = area centroid, stored as
coordinate system `sysno` (not activated; 0/blank = not stored). SOLID pieces: plane–element polygon,
`f_i = A_i(S·ez)`; SHELL pieces: segment of length L, width `w = t/|m × ez|`, membrane traction plus plate-bending
couples (default on, D-SEC-03); `F = Σf_i`, `M = Σ(c_i − C) × f_i + Σmb_i`; reported Fx, Fy, Fz (+ tension), Mx,
My, Mz. Faces lying on the plane are counted once. CALCSECTHIST(DB) write a 7-column CSV with a final signed
abs-max `MAX` row; time `(k−1)·ts` or step k when ts ≤ 0.

### 4.15 Option NON summary (P2) [M]+[R]

Per nonlinear element and iteration: x(t) from combined 3-direction THD (panel shear strain from the four corner
nodes, rigid-body invariant:
`γ = ½[(u_TL−u_BL)+(u_TR−u_BR)]/H + ½[(w_BR−w_BL)+(w_TR−w_TL)]/L`; bending curvature `κ = (θ_t − θ_b)/H`;
axial `ε_v`); `x_eq = EDF·max|x|`; hysteresis model (CMS, TAK from Cheng–Mertz 1989 / Takeda 1970; GMR = Masing
with factor 2) gives F(t); `E_new = E_el·K_sec/K_el` with K_el = BBC first slope Y1/X1 (must equal G·A_shear,
E·I or k, warning EDU-09); damping `ξ = min(cutoff, scale·ξ_h + [ξ_el])` (scale 0 → 1, cutoff 0 → none);
ductility `max|x|/x_cr`; `Fμ = K_el·max|x| / F(max|x|)`; convergence when max relative change of E < 2 % and of
ξ < 0.5 % (absolute), at most 10 iterations (D-NON-*). BBCGEN builds 22 points (cracking, 20 equal strain steps
to yield at index 21 with V = V_u, failure at (0.02, 1.02V_u)); SHEAR equations per spec 10 §3.11 in psi units
internally.

### 4.16 Conflict register (where sources disagree; resolution is binding)

| ID | Conflict | Sources | Resolution |
|---|---|---|---|
| C-01 | TF interpolant: rational in shifted linear x vs 2-DOF form in ω² | 05c B.6.4 vs R1 §5.1 (MHI transcription, Tajirian) | R1 ω² form (D-MOT-01) |
| C-02 | Options 4/5 averaging with the unshifted solution vs single window | R1 §5.2 vs 01 §9 / 05c ("no averaging") | single window, no averaging (manual table) |
| C-03 | Damping form `1+2iβ` (SHAKE 1972, some spec background) vs Udaka–Lysmer form | 05a §6.7 vs R1 §1, R2 §0 | SASSI form default; `1+2iβ` via CMODFORM (D-CNV-03) |
| C-04 | Half-space sublayer law: linear `H(2j−1)/n²`, geometric from thin layer, uniform | 03 item 20, R1 §2.4 (uniform most accurate in V4), 05a OQ12 | geometric, non-decreasing, each ≤ λ/8, uniform fallback; gate by VP-02b (D-SIT-02) |
| C-05 | NFFT rounding "nearest" vs "round up" | Warning 9 text, 05a (3000 → 2048) vs 07 OQ26 | nearest (manual) + EDU-03 if the record is truncated (D-CNV-10) |
| C-06 | Max FFT points 32,768 vs 65,536 | §1.5.2 vs §1.2 limits | 32,768 input generation, 65,536 response modules (D-CNV-11) |
| C-07 | Generated layers "10–20" vs Error 47 "0 or 4–20" | 05a OQ1 | accept 0 or 4–20; EDU warning below 10 (D-SIT-03) |
| C-08 | TOPL max 200 vs Warning 8 (> 100 truncated) | 05a OQ2, 11 OQ34 | 200 allowed; W8 only in legacy-limits profile (D-SIT-04) |
| C-09 | SITE `<opmode>` (solution/check) vs dialog Linear/Non-Linear Soil | 07 OQ19 | keep `<opmode>`; SITEX carries soil mode (D-SIT-06) |
| C-10 | MOTION `<bl>` "0 time / 1 frequency domain" vs dialog No/With correction | 07 OQ17 | 0 none, 1 Hudson–Housner time domain (D-MOT-08) |
| C-11 | ANALYS `[simul]` default 0 vs dialog 1 | 05c OQ11 | 0 single case; 1 = X/Y/Z coherent; Ns, Nl as counts (D-ANL-06) |
| C-12 | FILE9 "one digit" vs FILE9001…FILE9500 | 02 OQ12 | 3 digits |
| C-13 | Restart suffix "xxxxx" vs 3 digits | 03 OQ7 | 3-digit order number |
| C-14 | Smoothing range "1–1000" vs "0 for coherent" | 05c OQ2 | 0 = off; 1–1000 valid |
| C-15 | Errors 52/63 interval (0, 360) vs default 0° | 11 OQ26 | [0, 360) (D-CHK-03) |
| C-16 | Error 60 (≤ 0 illegal) vs dialog value 0 | 05b OQ9 | only when incoherent and embedded (D-CHK-04) |
| C-17 | Dynamic properties: 15 used (E96) vs 100 curves / 11 points | 05a OQ16 | E96 above 15 in legacy profile; 100 curves; > 11 points warning (D-SOL-05) |
| C-18 | Section-cut axes e1 = n vs ez = n, ex = r | 04 §8.3 vs 10 Q-30 | spec 10 |
| C-19 | RADIUS: `√(A/π)` vs `Scale·√A_plan` | 01 §6 rule 9 vs 09 Q-G10 | `Scale·√A_plan` (matches 0.9h) (D-PNT-05) |
| C-20 | Quarter models for SM/MSM validation: required (ASCE/SRP) vs not recommended (ACS) | 03 OQ16 | show both statements; enforce neither |
| C-21 | SRP criteria: manual (Δt ≤ 0.005 s) vs SRP Rev 4 (Δt ≤ 0.010 s, strong motion ≥ 6 s, \|ρ\| ≤ 0.16, PSD) | 02/05a vs R2 G.2 | check both, labelled separately (D-EQK-04) |
| C-22 | α = 0.1 decays faster in X (dialog text) vs in Y (formula) | 05b OQ3 | formula |
| C-23 | MERGESOIL mapping file in argument 4 (example) vs 7 (syntax) | 04, 09 | lenient: trailing non-numeric token = Mapping (D-MDL-12) |
| C-24 | Panel force 3/disp 2 documented vs Errors 126/128 | 05d OQ-N12 | keep Errors 126/128 (fidelity); EDUOPT NONEXT unlocks (D-NON-05) |
| C-25 | SOIL unit weight lb/ft³ vs kcf | 07 OQ31 | one model unit system; kcf/kN/m³ (D-SOL-09) |
| C-26 | PSD units cm/s³ vs cm²/s³ | 05a OQ9 | cm²/s³ / in²/s³ |
| C-27 | Interpolation options 0/1 naming (Ch. 1 vs Ch. 6) | 01 OQ28 | Ch. 6 mapping (0 SASSI2000, 1 SASSI 1982) |
| C-28 | SDOF test expectation 1/(2β) | 08 §11 test 12 | SASSI form gives `1/(2β√(1−β²))` = 10.0125 at β = 0.05 |
| C-29 | 05c test C8 viscous SDOF exact under options 0–5 | 05c C8 vs R1 §5.1 | use hysteretic 2-DOF/SDOF |
| C-30 | Table 3.1 examples SN001 = layer 2, SS001 = layer 3 | 02 OQ13 | typo: xxx = layer number |
| C-31 | RELDISP DB "Modulename_STRESS.bin" | 02 OQ16 | `Modelname_TR_d_THD.bin`, `Modelname_THD.bin` |
| C-32 | SHEAR/BBCGEN args 6–7 forces (syntax) vs area and yield stress (text) | 10 Q-14, 11 OQ8 | A_BE and f_y,BE (D-NON-10) |
| C-33 | INT code 3 and default code (text cut off) | 08 Q5 | 3 internal, default 0 |
| C-34 | LOCAL Y = X × Z (left-handed) | 08 Q7 | Y = Z × X |

---

## 5. User-interface requirements

### 5.1 Architecture requirements

| ID | Requirement | Tier |
|---|---|---|
| UI-01 | One command interpreter serves everything: Command Entry, INP files, macros, FOREACH, menus, dialogs and toolbar buttons. GUI actions only build command text and submit it (L17). GUI and scripted runs therefore give identical results | P0 |
| UI-02 | Three front ends share the same core library (`sassi`): (a) interactive GUI `sassi-gui` (tkinter + embedded matplotlib), (b) text console `sassi` (prompt with history; also `sassi run file.pre` for batch), (c) per-module batch entry points `python -m sassi.modules.<module> < <MODULE>.inp` that read the three-line protocol | P0 |
| UI-03 | All plotting works headless (matplotlib Agg); CAPTUREPLOT works in batch | P1 |
| UI-04 | Long operations (module runs, PROCFRAME) show progress in the status bar and stream output to a tab. The GUI stays responsive: modules run in a worker process; a Cancel button terminates them. During `.pre` execution RUNxxx is synchronous | P0 |
| UI-05 | The UI holds several numbered models (model 0 at start). Plot tabs and dialogs act on the active model, except the Cut plot | P0 |
| UI-06 | Every dialog commit is validated with the same rules as CHECK for the fields it owns; invalid values are refused with the Chapter 10 message text | P0 |
| UI-07 | Settings persistence: SASSIini.xml (module locations, extensions, Command Display filters, colours, fonts, shader options); SASSIdb.xml (group/model tree); SASSIani.xml (animations). Check Options and toolbar visibility are **not** persisted (manual fidelity) | P0 (ani P1) |
| UI-08 | Educational aids (extensions, clearly labelled): INTCOUNT cost/memory estimate before ANALYS; passing-frequency report; "Verification" menu running VERIFY; tooltips quoting the manual WARNING texts listed in specs 01–05d | P1 |

### 5.2 Main window [M]

Title `SASSI-EDU User Interface (ACS SASSI V3 methodology)`. Menu bar `Model | File | Plot | Modules | Options |
View | Help`. Two toolbars (Main, 3D Plot; §5.9). Tabbed document area: `Command History` tab first, then plot
tabs (`Model n - Model Plot`, …), editor tabs (`File Editor - <path>`) and module-output tabs. `Command Entry`
line at the bottom (Up/Down recall). Status bar with the progress bar at the lower right (current line / total
lines during INP).

### 5.3 Menu tree (each item calls the command shown)

| Menu | Item | Command / action | Tier |
|---|---|---|---|
| Model | New | `ACTM,<lowest unused number>` | P0 |
| | Open (Ctrl+O) | Load Model dialog (SASSIdb.xml tree; Open, Add/Remove Group, Add/Remove Model with two-step delete prompts) → `MDL` + `RESUME` | P0 |
| | Save | `SAVE` | P0 |
| | Input | file dialog → `INP,<file>` | P0 |
| | Converters ▸ SASSI .hou / ANSYS .cdb / GT-STRUDL | `CONVERT,SSI…` (P1) / `CONVERT,ANSYS…` (P1) / greyed (OOS) | P1 |
| | Output | `WRITE,<file>,<path>` | P0 |
| | Export to ANSYS / Export to STRUDL | `ANSYS,…` (P1) / greyed | P1 |
| | Exit | save SASSIini.xml; prompt for unsaved models (extension) | P0 |
| File | Open | text editor (create if missing; File ▸ Save; Input ▸ Connect to Command Entry — one window at a time appends every accepted command) | P0 |
| | Export Image | `CAPTUREPLOT,<file>` (.bmp default; .png by extension) | P1 |
| | Export Table | CSV of the active 2D plot (x column, one column per line, header of names) | P1 |
| Plot | Model ▸ Elements / Nodes | `MODELPLOT` / `NODEPLOT` | P1 |
| | Cuts | `CUTPLOT` (dialog: Cut Number, Model Number) | P1 |
| | Spectrum TFU-TFI / Time History | Line Selection dialog → `READSPEC`/`READTH` + `SPECPLOT`/`THPLOT` | P1 |
| | Soil Layers / Soil Properties | `LAYERPLOT` / `SOILPROPPLOT` | P1 |
| | Non Uniform Soil Field | greyed | OOS |
| | Process Animation Frame List | `PROCFRAME` (Parse Frame Data dialog) | P1 |
| | Bubble / Vector / Contour / Deformed Shape | Load Frame Data dialog → `BUBBLEPLOT`/`VECTORPLOT`/`CONTOURPLOT`/`DEFORMPLOT` | P1 |
| Modules | Location | Module Directories dialog: per module "built-in" (default) or an external executable path | P0 |
| | Extension | File Extension Options (module, input ext, output ext; OK saves all edits) | P0 |
| | EQUAKE … RELDISP | `RUN<MODULE>` on the active model (incl. SITE → `RUNSITE`) | P0 |
| | LIQUEF, PINT | greyed | OOS |
| | NONLINEAR | `RUNNONLINEAR` | P2 |
| | ANSYS Eq. Static Load / ANSYS Dynamic Load / ANSYS Super Element Utilities | LOADGEN / SSI2ANSYS dialogs (spec 04 §15.6–15.7) | P2 |
| Options | Model | Model Options (MOPT) | P0 |
| | Write | Extended Write Options (WRITE/AFWRITE) | P0 |
| | Check | Check Options (session only) | P0 |
| | Analysis | Analysis Options (tabs, §5.4) | P0 |
| | Windows Settings / Colors / Font / Shader Options / Reset Plot | `WINDOWSETTINGS` (P1) / Select Colors (P2) / Font Selection (P2) / `SHADEROPTIONS` (P2) / `RSTVIEW` (P1) | P1–P2 |
| View | Check Errors | window `CHECK: Errors and Warning for - <model>` showing `<model>.err` | P0 |
| | Command Window | show/reopen Command History (reopening clears its output) | P0 |
| | Command Display ▸ Command Echo / Output Confirmation / Comments / Warnings & Errors | toggles (Info always shown); persisted | P0 |
| | Toolbars ▸ Main / Plot | show/hide (not persisted) | P1 |
| Help | Help (F1) | local HTML help (manual-derived topic pages + spec cross-links) | P1 |
| | About | version, build, environment (Python, numpy, scipy versions) | P0 |
| | Verification (extension) | `VERIFY` with result table | P0 |

### 5.4 Options dialogs

**Options ▸ Model** (MOPT): Incompatible Modes (Include / **Suppress**), General Matrix (**Mass** / Weight), Overwrite
Mass (☑), Overwrite Force (☑).

**Options ▸ Write** (Extended Write Options): MDL command in *.Pre (☐); Extend integer fields (☐, accepted and
recorded, no effect on SASSI-EDU decks, D-AFW-03); AFWR Command in *.Pre (☐); Simulation Commands (☐) with
location `Write to *-Sim.Pre` (default) / `Write to *.Pre` — writes RUN commands for AOPT-enabled modules in
run order (D-AFW-04).

**Options ▸ Check**: Show Warnings ☑, Show Errors ☑, Suppress Error Window ☐, Break Check at [100] messages (per
type, per module). Session only.

**Options ▸ Analysis** — one window, tabs in this order:
`EQUAKE | SOIL | SITE | POINT | HOUSE | FORCE | ANALYS | MOTION | STRESS | RELDISP | NONLINEAR | AFWRITE`, buttons
Ok/Cancel. Shared fields are one stored value shown on several tabs (Δt, NFFT, frequency set and Δf owned by
SITE; control point layer by SITE; SSI gravity and ground elevation by HOUSE; analysis type by ANALYS;
history-file data by MOTION; coherency/wave passage/ME by HOUSE). Defaults below are the **new-model defaults**
(D-UI-03); screenshot values in the manual are examples only.

| Tab | Field (dialog label) | Type | Default | Stored in |
|---|---|---|---|---|
| EQUAKE | Spectrum Number | 1–3 | 1 | index of RSIN/RSOUT/ACCIN/ACCOUT/TPSD |
| | Spectrum Input File / Spectrum Output File / Acceleration Output File (+ Edit, + Library) | path | — (blank: built-in default, §7.19) | RSIN, RSOUT, ACCOUT |
| | Accel. Record / External Accel | check (exclusive) | off | EQUAKE `<accopt>` 1 / 2 |
| | Acceleration Input File | path | — | ACCIN |
| | Use Target PSD / PSD File | check / path | off | EQUAKE `[tpsd]`, TPSD |
| | Number of Frequencies | int > 0 | (records of RSIN 1) | `<nrfreq>` |
| | Initial Random SEED | int > 0 | 11975 | `<rand>` |
| | Damping Value | real | 0.05 | `<damp>` |
| | Time Step | real | 0.005 | SITE `<delt>` |
| | Total Duration | real > 0 | 20 | `<dur>` |
| | Number of SEEDs | int | 1 | `<seeds>` |
| | Correlated + Correlation grid (Time, Corr.) | check / table | off | `<corr>`, CORR (P1) |
| | Spectra Title | text | "" | EQTIT |
| SOIL | Number of Fourier Components / Time Step of Input Motion | shared | 4096 / 0.005 | SITE `<nft>`, `<delt>` |
| | Number of Values | int > 0 | — | `<nrval>` |
| | Multiplication Factor / Max Value for Time History (g) | real (exactly one ≠ 0) | 0 / 0.1 | SOILX `<mult>`, `<max>` |
| | Gravity Accel. (free field) | real > 0 | 32.2 | SOIL `<grav>` |
| | Number of Header Lines | int ≥ 0 | 0 | `<header>` |
| | Input Direction | 0/1 | 0 | SOILX `<indir>` |
| | Control Point Layer | int | 1 | SITE `<cl>` (SOIL-only override SOILX `<cl>`) |
| | File (+ Library) | path | — (blank: built-in record, §7.19) | THFILE (SOIL-only override SOILX `<file>`) |
| | Assign as Outcrop Motion | check | on | `<outcrop>` |
| | Save Strain-Compatible Soil Properties | check | on | `<save>` |
| | Number of Iterations | int ≥ 0 | 8 | `<iter>` |
| | Equiv. Uniform / Max Strain | (0,1) | 0.65 | `<ratio>` |
| | Layer Number / Property Number / Dynamic Soil Property (…) | spin / int / label + dialog | 1 / — / — | SPRO, DYNP |
| | Accelerations (No / Maximum / Maximum + TH; Outcropping) | radio / check | Maximum | SACC |
| | Response Spectrum (Save, Outcropping, Multiplier for g, Damping Ratios) | | on, off, 1, "0.02,0.05" | SRS, SOIL `<gravmult>`, DAMP |
| | Stresses/Strains (4 checks) | | compute both, save none | SSTR |
| | Spectral Amplification Factor (save, outcrops, second layer, step, title) | | off | SSAF |
| | Fourier Spectrum | disabled | — | SFOU |
| | Nonlinear Soil group (P2) | see spec 05a §6.6 | off | NLSOIL, NLSLAYER |
| SITE | Linear Soil / Non-Linear Soil | radio | Linear | SITEX `<soilmode>` |
| | Mode 1 / Mode 2 | checks | on / on | `<mode1>`, `<mode2>` |
| | Gravity Accel. | real | 32.2 | HOUSE `<gravity>` |
| | Frequency Step | real ≥ 0 | 0 | `<fstep>` |
| | Time Step Control Motion / Nr. of Fourier Component | real / power of 2 | 0.005 / 4096 | `<delt>`, `<nft>` |
| | Frequency Set Number | int | 1 | `<freq>` |
| | Number of Generated Layers | 0 or 4–20 | 20 | `<nl>` |
| | Halfspace Layer | L number | — | `<hs>` |
| | Top Layers | list (≤ 200; blank, tab, `,`, `;`, Enter) | — | TOPL |
| | R-, SV- and P-Waves / SH- and L-Waves | radio | R/SV/P | `<wopt>` |
| | Wave pages (field option, Wave Ratio 1, Wave Ratio 2, Incident Angle on body-wave pages) | | SV field, 1, 1, 0 | WAVE |
| | Frequency 1 / Frequency 2 (frequency numbers; Hz shown alongside) | int | 1 / NFFT/2 | `<freq1>`, `<freq2>` |
| | Control Point Layer / Direction (X, Y, Z) | int / radio | 1 / X | `<cl>`, `<cm>` |
| POINT | Operation Mode | radio | Solution | `<opmode>` |
| | Number of Embedment Soil Layers | int ≥ 0 | 0 | `<layer>` |
| | Point Load Central Zone Radius (+ "From mesh" helper → RADIUS) | real > 0 | — | `<rad>` |
| HOUSE | Operation Mode / Dimension (2D, 3D; 1D greyed) / Flexible Volume Method (FV, FFV, FI) | radio | Solution / 3D / FV | `<opmode>`, `<dim>`, `<imp>` |
| | Acceleration of Gravity / Ground Elevation | real | 32.2 / 0 | `<gravity>`, `<gelev>` |
| | Non-Linear SSI + Input Data (`.pin` editor) / Optimize Model | check / button / check | off | HOUSEX `<nlssi>`, `<optimize>` |
| | Soil Motion Coherent / Incoherent | radio | Coherent | `<coh>` |
| | Coherence Parameter X / Y / Z Dir | real ≥ 0.1 | 0.1 / 0.1 / 0.2 | INCOH 1–3 |
| | Alpha Directionality Factor (Vs for model 1) | real | 0.5 | INCOH `<alpha>` |
| | Number of Embedded Layers | int | 1 | INCOH `<ngp>` |
| | Number of Incoh. Modes / Print Coherency Matrix | int / check | 0 / off | INCOH `<nmodes>`, `<ipr>` |
| | Use Multiple Excitation; Input Motion Number; First/Last Foundation Node; X/Y/Z Coord. (unused) | | off | HOUSE `<me>`, ME |
| | Spectral Amplification (text) / Use Complex Spectral Amp. | text / check | — / off | AMP, `<cmplxspec>` |
| | Non-Uniform Motion / Non-Uniform Soil | disabled | — | — |
| | Use Wave Passage; Apparent Velocity for Line D; Angle Line D with X Axis; Unlagged Coherency Model | | off; 1e9; 0; 1 | `<wpass>`, WPASS |
| | Deterministic (Median) / Stochastically Simulated; Horizontal SEED; Vertical SEED; Random Phase Angle; Number of Simulations | | Deterministic; 0; 0; 0 (180 when stochastic); 1 | INCOH 9–11, HOUSEX `<nsim>` |
| | Superposition Mode Linear / Quadratic | radio | Linear | HOUSEX `<supmode>` |
| | Ansys Model Input / Ansys Model Type (P2) | check / radio | off / Embedded | HOUSEX `<ansys>`, ANSYSMODELTYPE |
| FORCE | Operation Mode; Gravity; Frequency Step; Time Step; NFFT; Frequency Set | | Solution; shared | FORCE `<opmode>`; shared |
| ANALYS | Operation Mode / Type of Analysis / Mode of Analysis | radio | Solution / Seismic / Initiation | `<opmode>`, `<type>`, `<mode>` |
| | Simultaneous Cases | int | 0 | `[simul]` |
| | Save Restart Files / Delete Restart Files | check | off / off | `<save>`, ANALYSX `<delrst>` |
| | Take Frequency Numbers from File1/File9 / Frequency Set Number | check / int | off / 1 | `<fopt>`, SITE `<freq>` |
| | X-, Y-, Z-Coordinate of Control Point / Coordinate Transformation Angle | real | 0 | `<xc>`, `<yc>`, `<zc>`, `<ang>` |
| | Coherent / Incoherent; Wave Passage Effects Included (mirrors of HOUSE) | | — | HOUSE |
| | Free-Field Load / Free-Field Motion (enabled when incoherent) | radio | FFL | ANALYSX `<ffm>` |
| | Print Amplitude Only | check | on | `<prnt>` |
| | Multiple Excitation block (mirror of HOUSE) | | | HOUSE, ME |
| | Global Impedance: None / Diagonal / Full 6×6 | radio | None | `<impe>` |
| MOTION | Operation Mode / Type of Analysis | radio | Solution / (ANALYS) | `<opmode>`, ANALYS `<type>` |
| | Baseline Correction No / With | radio | No | `<bl>` |
| | Response Spectrum: First / Last Frequency / Total Number of Freq. Steps / Damping Ratios | real / int / list | 0.1 / 100 / 301 / "0.02 0.05" | `<freq1>`, `<freq2>`, `<fstep>`, DAMP |
| | Output Only Transfer Functions / Save Complex Transfer Functions | check | off / on | `<out>`, `<cplx>` |
| | Save FILE 12 or FILE 13 / Total Duration to be Plotted | int / real | 0 / 0 | MOTIONX `<f1213>`, `<dur>` |
| | Response type for Foundation Vibration (Displacement / Velocity / Acceleration) — shown only for vibration | radio | Acceleration | MOTIONX `<resp>` |
| | Incoherent SSI: Input (SRSSTF.txt editor) / use SRSS | button / check | off | MOTIONX `<srss>` |
| | Interpolation Option / Phase Adjustment / Smoothing Parameter | int 0–6 / 0–1 / real ≥ 0 | 1 / 0 / 0 | `<interp>`, `<pzadj>`, `<smo>` |
| | Node List (Add/Edit/Delete; "1, 3-6 10") × Direction (X…ZZ) × six flags | list | — | NOUT |
| | Nr. of Fourier Components / Time Step (shared) / Multiplication Factor / Max Value / First Record / Last Record / Title / File / File Contains Pairs | | shared / 1 / 0 / 1 / last / — / — (blank: built-in record or load pulse, §7.19) / off | MOTION args, THFILE, THTIT |
| | Convert TH to RS: Select External Files / Input Time History Files | check / button | off | `<cnvrt>`, CONTTRS.txt |
| | Post Processing: Save TF / ACC / RS in All Points; Save Rotation for Ansys; Restart for TF / ACC / RS | checks | off | MOTIONX flags (P1) |
| | Save Binary Database | check | off | BINOUT `[mot]` (P2) |
| STRESS | Operation Mode / Type of Analysis | radio | Solution / (ANALYS) | `<opmode>` |
| | Auto Computation of Strains in Soil El. / Save Stress Time Histories / Output Transfer Function | checks | off / on / off | `<iter>`, `<save>`, `<itran>` |
| | Phase Adjustment / Interpolation Option / Smoothing Option / Skip Time History Steps | | 0 / 1 / 0 / 0 | STRESSX, `<interopt>` |
| | Acceleration Time History Data (shared with MOTION) | | | MOTION, THFILE |
| | Element Output Data list (Group, Element List, Output Code) + Components radio + Component Request radio | list | — | EOUT |
| | Post Processing: Save Max Value; Save Time History; Restart for Nodal Stress Contours; Restart for Soil Pressure Contours; Frame Selection (Frames.txt editor) | | off | STRESSX (P1) |
| | Section Cut Options: Save Time History | check | off | SECDATAOPT (P1) |
| | Save Binary Database | check | off | BINOUT `[str]` (P2) |
| RELDISP | Complex TF File Name | path | — | RELFILE |
| | Save Rel Disp Complex TF | check | on | RELD `<RelDisOutput>` |
| | Acceleration Time History Data (shared) | | | MOTION |
| | Nodal Output table (Node, X…ZZ) + Add/Edit/Delete | | — | RDND |
| | Save Relative Displacement in All Nodes / Save Rotations for ANSYS V11.0 / Restart For Frame Generation | checks | off | RELD `<RelDispSAll>`, RELDX |
| | Binary Disp. Option: No / TFD (greyed) / THD | radio | No | BINOUT `[reldisp]` (P2) |
| NONLINEAR (P2) | Disp. Factor; Damping Cutoff %; Damping Scale Factor; Material Parameter (greyed); Use Non-linear Panels / Springs / Beams (greyed); Include Elastic Damping; BBC grid (view/edit only); Panel, Spring, Beam (greyed) data | | 0.8; 0; 1; —; off; off | EQL, BBC*, P, S, B |
| AFWRITE | EQUAKE, SOIL, LIQUEF (greyed), SITE, POINT, HOUSE, PINT (greyed), FORCE, ANALYS, COMBIN, MOTION, STRESS, RELDISP, NONLINEAR | checks | SITE, POINT, HOUSE, ANALYS, MOTION on; others off | AOPT |

**Enable/disable logic** (D-UI-04): Coherent → incoherency fields disabled; 2D → Incoherent disabled; Unlagged
model 2–7 or Multiple Excitation → Wave Passage forced on (with message); FFL/FFM enabled only when
incoherent; smoothing disabled with option 6; vibration response radio only for Foundation Vibration; ANSYS
type only with ANSYS Model Input.

### 5.5 Check Errors window and `.err` format [M]

```
CHECK: Errors and Warning for - <model>
Errors and Warnings for EQUAKE
Error 85 : RS Input File 1 Does Not Exist
Errors and Warnings for SITE
Warning 9 : Number of Values for Fourier Transform Is Not Power of 2
...
<MODULE>: <n> errors, <m> warnings   (totals line, extension)
```
Module headers appear in run order for every AOPT-enabled module (and a `MODEL` header for model-level checks,
D-CHK-01). EDU checks print as `Warning EDU-nn : <text>`.

### 5.6 Command Entry, Command History, File Editor [M]

- Command Entry: single line, Up/Down history, Enter submits.
- Command History: read-only, selectable/copyable; six message classes with user colours (echo black,
  confirmation blue, error red, comment green, info purple, warning amber); View ▸ Command Display filters
  (info cannot be hidden); during INP echo and confirmation are suppressed.
- File Editor: open/create any text file; File ▸ Save; Input ▸ Connect to Command Entry (one window at a time
  appends each accepted command); a "Run" action executes the buffer through INP (extension).

### 5.7 Plots (P1 unless noted; normative detail in spec 06)

| Plot | Content | Interaction / options |
|---|---|---|
| Element (MODELPLOT) | shaded elements coloured by group / material / property (ElemPalette, 128 colours, index `((n−1) mod 128)+1`); beams/springs as lines; mass markers coloured by direction | wireframe, shrink, labels (node / element ⟂ group), SHOWDOF, SHOWMASS, hide/show, display volume |
| Node (NODEPLOT) | element-connected nodes only: black ordinary, red interaction, green border fixed DOF (SHOWDOF), purple border mass, blue square selected | NODESEL, labels |
| Cut (CUTPLOT) | wireframe of a model with the cut's elements filled; any model number | dialog Cut Number / Model Number (default active model) |
| Spectrum (SPECPLOT) / Time History (THPLOT) | up to 50 line objects; legend; generic XY axes | AXES, PLOTRANGE, titles, MARKERS, STIPPLE; Graph Plot Options dialog with Average, SRSS, Broaden, Linear Combin (3 lines), Addition, Subtraction (no file writing) |
| Soil Layer (LAYERPLOT) | layer column proportional to thickness + table (Unit Weight, Vp, Vs, Dp, Ds; Thickness optional), Halfspace row | Start Layer 1, EndLayer −1 |
| Soil Properties (SOILPROPPLOT) | damping (left axis, red) and G/Gmax (right axis, green) vs Shear Strain % (log axis default) | Select Dynamic Soil Property dialog (New/Edit/Delete; 4-column table) |
| Bubble / Contour | jet colour map clamped (dark blue below MnR, dark red above MxR), 5 colour-bar labels; contour interpolates nodal scalars across faces | Load Frame Data (Start, End, Stride, range), PAUSE, ± step |
| Vector | three vectors per used node (X red, Y green, Z blue); complex `a+ib` drawn as (a, b, b) for X (and analogues) × Scale | Output Direction X/Y/Z/All |
| Deformed | `x′ = x + Scale·u`; optional undeformed wireframe | Scale, Frame Pause 33 ms |

3D controls: right-drag rotate, left-drag zoom, middle-drag pan, Shift+middle zoom box, left double-click
node id, right double-click element id; Insert/Delete, Home/End, PageUp/PageDown rotate 5° per press;
orthographic projection; axis triad. Implementation: matplotlib `mplot3d` with a vectorised poly collection
for models up to ~50k faces; a plotly HTML export for large models and animations (D-UI-06).

### 5.8 Animation pipeline (P1)

Frame files from module Restart options → PROCFRAME (list file, first line ignored) → frame store (one `.npy`
per frame plus `index.json`) + SASSIani.xml entry (description, directory, type, frames, defaults) → animation
commands. Frame-file layout (D-FIL-03): header `nrows ncols`, then rows `node v1 … v(ncols−1)` (node id counted
as column 1); complex frames store Re/Im pairs.

### 5.9 Toolbars [M] (spec 06 §11)

Main: New, Open DB, Save, Output .pre, Output ANSYS, .hou converter, ANSYS converter, STRUDL (greyed), Capture,
Element plot, Node plot, Cut plot, Spectrum, Time History, Soil Layer, Soil Properties, Process Frames, Bubble,
Vector, Contour, Deformed, Analysis Options (opens §5.4), AFWRITE.
3D Plot: Change View, Reset View, Change Centre, Reset Centre, Wireframe, Shrink, Colour by Group/Material/
Property, Node/Element/Group labels, Show DOF, Show Mass, Pause/Start, Debug. Buttons not applicable to the
active plot are disabled; each button submits its command text.

### 5.10 Persistence and locations

| File | Location (D-UI-07) | Written |
|---|---|---|
| SASSIini.xml | `~/Library/Application Support/SASSI-EDU/` (macOS), `%APPDATA%\SASSI-EDU\`, `~/.config/sassi-edu/` | OK of Location/Extension/Colors/Font/Shader dialogs; Exit |
| SASSIdb.xml | same directory | Load Model dialog changes |
| SASSIani.xml | same directory | PROCFRAME, *DBANI, Remove Animation |
| `<model>.sdb`, `.pre`, decks, results | model directory (MDL path) | SAVE, WRITE, AFWRITE, modules |

---

## 6. Verification plan

### 6.1 Principles

1. Every verification problem (VP) is an automated pytest case under `tests/verification/vp_NN_*.py`. It is
   also runnable from the UI with `VERIFY,VP-NN`. Where modules are involved, the VP is built as an
   ordinary `.pre` model plus a comparer, so a student can open, run and inspect it in the GUI.
2. Reference values come from R2 (published P, exact E, approximate A, derived D) and R1 checks [V]. Data files
   are in `docs/spec/R2_benchmark_data/`. A VP fails if any quantity is outside its tolerance. Tolerances
   marked *provisional* must be calibrated once (record the observed error) and then frozen.
3. A tier is complete only when every VP of that tier passes. P0 release gate: all P0 VPs and unit tests pass.
4. Each VP report states the inputs, the computed and reference values, the error and the tolerance. This
   is the "verifiable" requirement of the user request.

### 6.2 Interpreter, data-model and I/O tests (UT)

| ID | Test | Expected | Tier |
|---|---|---|---|
| UT-01 | Every full name and abbreviation of §3.2 resolves; `ACC`, `ANALY`, `RELF` → `<X> Command not found` and later lines still run; names case-insensitive, paths case-kept | as stated | P0 |
| UT-02 | Blank fields take defaults (`cutvol,3,,,,,2.53`); titles with commas; legacy SOIL `(5X,F10.4)` stays one token and raises the legacy warning | as stated | P0 |
| UT-03 | WRITE → INP into an empty model gives deep-equal state (nodes, elements, tables, options, FREQ, DAMP, TOPL, AMP, NOUT, EOUT, RDND, DYNP labels, SYMM, X-commands) | identical | P0 |
| UT-04 | List semantics: `DAMP,0.02` then `DAMP,0.05` → [0.02, 0.05]; `DAMP,0` clears; same for FREQ (per set), TOPL, AMP | as stated | P0 |
| UT-05 | AFWRITE gating: invalid SITE `<nl>` = 3 blocks `.sit` only; other decks written; existing stale `.sit` renamed `.sit.bak` with a warning | as stated | P0 |
| UT-06 | Shared variables: changing SITE `<delt>` changes Δt seen by SOIL, MOTION, STRESS, RELDISP, EQUAKE decks; `GRAVITY,9.81` changes only HOUSE argument 1 | as stated | P0 |
| UT-07 | MOTION scaling: mult/max exclusivity (Errors 77/78) and `a·max/max\|a\|` on a synthetic history | exact | P0 |
| UT-08 | CHECK catalogue: one minimal model per Error 1–128 and Warning 1–11 that violates exactly that rule; fatal stop for 42/43; no Error 8 for zero-length springs; break count honoured | each message number appears | P0 |
| UT-09 | Macros (manual 5.6.1–5.6.3): `MACRO,Node,13.52,15,100.25` → `N,1,13.52,15,100.25`; graphing and nested drivers produce identical file names | as stated | P1 |
| UT-10 | Loop.pre: `ForEach,Z,ForEach,Y,ForEach,X,N,@NNUM++,@X[#],@Y[#],@Z[#]` → 125 nodes, node `25(iz−1)+5(iy−1)+ix` at (ix,iy,iz); counter test `VAR,X,1.23,2.83,3`; `SETVAR,@X+5` → 5; `@X++` → 6; `@X=-15` → −15; `VAR,X,9` → 0; nested same-variable FOREACH → error, not run | as stated | P1 |
| UT-11 | REDUCESET: `VAR,S,10,9,1,9`; `REDUCESET,S` → `1,10,9`; `REDUCESET,S,INT` → `1,9,10`; `RNDSEED,7` + `RND,R,1000,UNI,2,4` reproducible, all in [2,4], mean ≈ 3 | as stated | P1 |
| UT-12 | Coordinate systems: `LOC,1,0,10,0,0,45,0,0`, `CSYS,1`, `N,100,5,0,0`, `GLOBAL,100,100,1` → (13.5355, 3.5355, 0), CSYS still 1; `LOCAL,2,0,1,2,3` with n1 = (0,0,0), n2 = (0,1,0), n3 = (−1,0,0) equals `LOC,3,0,0,0,0,90,0,0` to 1e-12; `LOC,4,0,0,0,0,0,90,0` maps local y to global +Z | as stated | P0 |
| UT-13 | Generation: `N,1,0,0,0` `N,5,40,0,0` `FILL` → nodes 2–4 at x = 10, 20, 30; `NGEN,3,5,1,5,1,0,0,10` → nodes 6–20 at z = 10, 20, 30; `EGEN` grid example → 12 elements; `NSCALE,1,1,,0,2,0` scales y only; `D,1,1,,,ALL` frees all DOFs | as stated | P0 |
| UT-14 | Materials: E = 30000, ν = 0.25, γ = 0.15, g = 32.2 entered as types 1, 2, 3 give the same (G, M, Vs, Vp) to 1e-10; circle r = 1: I2 = 0.785398, J = 1.570796; rectangle b = 0.5, h = 1: J = 0.028610 | as stated | P0 |
| UT-15 | Mass units: MT = 32.2 with MUNITS = 1 equals MT = 1 with MUNITS = 0; MXM weight units with MOPT `<matrix>` = 1 equals mass units | exact | P0 |
| UT-16 | Frequency grid: Δt = 0.005, NFFT = 4096 → Δf = 0.048828125 Hz; f_cut = 25 Hz → NFREQ 512; Δt = 0.01, NFFT = 2048 → 0.048828125 Hz; Δf = 0.9765625 Hz with numbers {1,3,5,7,9,11,13,15,16,18} → {0.977, 2.930, 4.883, 6.836, 8.789, 10.742, 12.695, 14.648, 15.625, 17.578} Hz; NFFT = 3000 → 2048 with Warning 9 (+EDU-03 if records > 2048); duplicates → error | exact | P0 |
| UT-17 | Naming: incoherent sim 17 dir Y → `FILE8050`, sim 50 dir Z → `FILE8150`; node 12 X RS damping 1 → `00012TR_X01.RS`; 6-digit nodes above 99,999; `BEAMS_003_00045_MXJ.THS` | exact | P0 |
| UT-18 | Memory estimate: N_int = 10,000 → 14.4 GB per dense matrix; 20,000 → 57.6 GB | exact | P0 |
| UT-19 | Half-space buffer: Vs_hs = 1500 m/s → H = 450 m at 5 Hz, 75 m at 30 Hz; generated thicknesses sum to H, are non-decreasing, each ≤ λ/8 | exact | P0 |
| UT-20 | Strain SRSS: γ = 0.10, 0.08, 0.02 % → 0.1296 % | exact | P1 |
| UT-21 | Section-cut piece geometry and resultant algebra (VP-49 inputs) | see VP-49 | P1 |

### 6.3 Physics verification problems

| ID | Problem | Modules verified | Tier | Type | Setup | Expected reference | Tolerance | Source |
|---|---|---|---|---|---|---|---|---|
| VP-01 | Fixed-base hysteretic SDOF transfer-function peak | HOUSE (SPRING/mass), ANALYS (no soil: fixed base), MOTION TFU, damping convention | P0 | E | k, m, β = 0.02…0.20, fine frequency set | peak `1/(2β√(1−β²))` = 25.0050, 10.0125, 5.0252, 3.3715, 2.5516 at `ω0√(1−2β²)` (ω/ω0 = 0.99960, 0.99750, 0.98995, 0.97724, 0.95917); CMODFORM,1 gives 25.0200, 10.0499, 5.0990, 3.4801, 2.6926; recovered damping 1/(2\|H\|max) = 4.99/9.95/14.83/19.60 % (SASSI2000 published 5.0/9.9/14.7/19.5 %) | 1e-6 rel. (frequency step permitting) | R2 §0 |
| VP-02a | Uniform layer TFs, closed form | SOIL (0 iterations) | P0 | E | H = 30 m, Vs = 200 m/s, ρ = 2.0, ξ = 5 %; rock Vr = 1000 m/s, ρ = 2.2, ξ = 1 % | \|surf/base-within\| (SASSI form) at 0.5, 1.0, 1.5, 1.6667, 2.0, 3.0, 5.0, 8.3333 Hz: 1.1216, 1.6931, 5.7710, 12.7152, 3.1156, 1.0411, 4.2038, 2.4815; \|surf/outcrop\|: 1.1150, 1.6207, 3.4116, 3.8327, 2.4115, 1.0087, 2.3554, 1.6702; CMODFORM,1: 1.1209, 1.6878, 5.6735, 12.7633, 3.1590, 1.0436, 4.2202, 2.4918 and 1.1143, 1.6161, 3.3916, 3.8344, 2.4313, 1.0110, 2.3610, 1.6756; any sub-division of the layer identical | 1e-8 rel. | R2 A.2 |
| VP-02b | Same site in SITE (TLM + variable-depth half-space) | SITE Mode 1/2, sublayer law | P0 | D/E | layer split to h ≤ λ/20; nl = 20 and 10 | within TF converges to VP-02a; outcrop (half-space top) TF within 1 % (nl = 20) and 2 % (nl = 10) of the exact `1/[cos k*H + iα* sin k*H]`; R1 V4 achieved ≤ 0.7 % / ≤ 1.5 % with uniform sublayers | 1 % / 2 % (provisional); failure switches D-SIT-02 default to uniform | R1 §2.4, V4 |
| VP-03 | Discrete soil column dispersion | HOUSE SOLID/PLANE mass, SITE layer mass | P0 | D | fixed–free column of Ne equal elements; plane wave kh | ω_h/ω for consistent / lumped / 50-50 at N = 5 per λ: 1.06632 / 0.93549 / 0.99451; column Ne = 2, 4, 8 (50-50): ω1 ratio 0.99919, 0.99995, 1.00000; ω2 ratio (Ne = 2) 0.92712 | 1e-6 | R2 A.3 |
| VP-04 | SHAKE91 sample problem | SOIL (complete) | P0 | P | `R2_benchmark_data/shake91_example/` (INP.DAT, DIAM.ACC): 16 sublayers + rock, g = 32.2, NFFT 4096, Δt 0.02, 1900 values scaled to 0.10 g, cut-off 25 Hz (`EDUOPT,SOILCUTOFF,25`), outcrop input at sublayer 17, 8 iterations, ratio 0.50, log-strain interpolation | surface PGA 0.19037 g; base within 0.07617 g; sublayer 1: γ_eff 0.00077 %, G 3851.5 ksf, G/Gmax 0.992, D 0.007; sublayer 9: 0.01356 %, 5402.8, 0.792, 0.034; sublayer 16: 0.00865 %, 11292.4, 0.863, 0.026 (full Table B-2 in R2 H.1); \|surf/outcrop\| 1.3270 (1 Hz), 3.3181 (2 Hz), 3.4369 (2.125 Hz), 1.7386 (3 Hz), 2.0450 (5 Hz), 1.6690 (10 Hz); max outcrop amplification 3.44 at 2.12 Hz; final surface/within-base max 20.47 at 2.11 Hz; 5 % SA (g) at T = 0.10, 0.40, 1.0 s: 0.3977, 0.7832, 0.1351 | PGA 0.5 %; strain, G, D 1 %; SA 2 % | R2 H.1 |
| VP-05 | Rayleigh phase velocity | SITE Mode 1 | P0 | E/D | deep uniform stratum, h = λ/20 | c_R/Vs = 0.919402 (ν = 0.25), 0.932526 (1/3), 0.948960 (0.45), 0.954074 (0.49); R1 V2: 0.9202 mixed mass at ν = 0.25 | < 0.5 % at ≥ 10 sublayers per λ; O(h²) convergence | R2 D.5, R1 V2 |
| VP-06 | Love-wave dispersion | SITE Mode 1 (SH) | P0 | E/D | H = 10 m, β1 = 200, β2 = 400 m/s, equal ρ | cut-offs 0, 11.547, 23.094 Hz; fundamental c = 387.536 (2 Hz), 300.026 (5 Hz), 224.175 (10 Hz), 205.949 (20 Hz) m/s; mode 1 c = 279.222 m/s at 20 Hz | 0.5 %, convergence from above | R2 D.5 |
| VP-07 | Static point-load Green functions and reciprocity | POINT3, ANALYS F_ff | P0 | E | quasi-static (f → 0.01 Hz), deep graded stratum | Boussinesq `u_z = P(1−ν)/(2πμr)`, `u_r = −P(1−2ν)/(4πμr)`; Cerruti `u_x = Q/(2πμr)[(1−ν)+νx²/r²]`, `u_z = +Q(1−2ν)x/(4πμr²)`; F_ff symmetric | 2 % at r = 1–5 m (R1 V3: 0.1–1.7 %); `‖F−Fᵀ‖/‖F‖ < 1e-3` before symmetrisation | R2 D.1, R1 V3/V10 |
| VP-08 | Dynamic surface Green functions (Wong 1975) | POINT3 | P0 | P | near-elastic half-space (β ≤ 1e-4, many thin sublayers + 20 generated layers), ν = 0.33, r0 = ωr/Vs = 0…5.5 | μru/P (Re, Im), e.g. r0 = 0: U_r0 −0.027, U_z0 0.106, U_r1 0.159, U_θ1 −0.106; r0 = 1.0: (−.033,.025), (.037,−.102), (.112,−.105), (−.045,.099); r0 = 3.0: (.074,.035), (−.114,.053), (−.064,−.090), (.132,−.020); r0 = 5.0: (.005,−.120), (.127,.064), (−.054,.004), (−.088,−.129); full table R2 D.3 | ±0.005 absolute | R2 D.3 |
| VP-09 | POINT far field vs exact Kausel point load | POINT3 | P0 | D | 28 layers, f = 5 Hz, R0 = 0.9 m | agreement with R1 §3.2 series | ≤ 2.3 % at 3.3R0; ≤ 0.8 % at 6.7R0 (dominant components); ≤ 0.25 % beyond 13R0 | R1 V8 |
| VP-10 | Rigid disk static stiffness on a half-space | POINT, ANALYS global impedance (`<impe>` = 2), surface foundation | P0 | E | rigid massless disk (rigid shell mat), fine interaction mesh, low frequency, ν = 1/3 | K_v = 4GR/(1−ν) = 6.0000 GR; K_h = 8GR/(2−ν) = 4.8000 GR; K_r = 8GR³/(3(1−ν)) = 4.0000 GR³; K_t = 16GR³/3 = 5.3333 GR³ | 5 % after mesh extrapolation (welded vs relaxed contact) | R2 B.1 |
| VP-11 | Rigid disk dynamic impedance | ANALYS global impedance | P0 | A/E | as VP-10, a0 = ωR/Vs ≤ 1.5 and a0 ≈ 6–8 | Veletsos–Verbic ν = 1/3: horizontal k = 1, c = 0.65; at a0 = 0.5 / 1.0 / 1.5: k_r 0.9529 / 0.8400 / 0.7120, c_r 0.0235 / 0.0800 / 0.1440, k_v 0.9517 / 0.8634 / 0.7934, c_v 0.7886 / 0.8593 / 0.9152; high-frequency c(∞): horizontal π(2−ν)/8 = 0.654, vertical π(1−ν)Vp/(4Vs) = 1.047, rocking 0.393, torsion 0.295 | ±10 % (a0 ≤ 1.5); ±5 % asymptotes | R2 C.1, C.2 |
| VP-12 | Embedded rigid cylinder | ANALYS (FV embedded), global impedance | P1 | P/E | cylinder radius a, depth h·a, ν = 1/4 | high-frequency dashpots (ρVs a² or ρVs a⁴): h = 0: C_HH 3.142, C_VV 5.441, C_TT 1.571, C_MM 1.360; h = 0.5: 7.433, 8.583, 4.712, 3.289; h = 1: 11.725, 11.725, 7.854, 7.363; h = 2: 20.308, 18.008, 14.137, 30.532; static torsion K_t/(Ga³) = 5.33, 13.23, 19.89, 32.75 (h = 0, 0.5, 1, 2) | ±5 % | R2 B.3, C.1 |
| VP-13 | Tajirian & Tabatabaie (1985) disk on a layer over rigid base (first end-to-end benchmark) | SITE, POINT, ANALYS | P0 | P/E | r = 0.5, H = 3, Vs = 1, Vp = 2, ρ = 1, β = 5 % and 15 %, 24 sublayers of 0.125 | resonances: vertical A0 = 0.524, 1.571, 2.618; horizontal A0 = 0.262, 0.785, 1.309; static `K_h ≈ 8GR/(2−ν)(1+R/2H)`, `K_r ≈ 8GR³/(3(1−ν))(1+R/6H)`; compliance curves vs digitised Figs. 4–11 (action item) | resonance A0 ±2 %; statics ±5 % | R2 C.7, B.4 |
| VP-14 | Square rigid footing statics | ANALYS global impedance | P0 | A | L = B, ν = 1/3 | Pais–Kausel (GB, GB³): K_z 7.0500, K_x = K_y 5.5200, K_zz 8.3100, K_xx = K_yy 6.0000 (exactly symmetric); equal-area circle K_z 6.7703 | ±5 % | R2 B.2 |
| VP-15 | Surface rigid massless foundation under vertical coherent SH/SV | SITE, POINT, HOUSE, ANALYS, MOTION | P0 | E | rigid massless mat on layered soil, vertical SV (X) and SH (Y) | foundation ATF H_u = 1 and rocking H_θ = 0 at every frequency | \|H_u − 1\| < 1 %, \|H_θ\|·L < 1 % | R2 E6 |
| VP-16 | Zero-SSI identity (FV) | HOUSE, ANALYS | P0 | E | embedded model whose "structure" equals the excavated soil (same mesh and layer properties), FV, vertical SV | U = U′_f at every node; ATF equals SITE free-field TF; control point ATF = 1 | 1e-8 rel. | 02 Part 6, R1 §4.4 |
| VP-17 | Inertial SSI 3-DOF | HOUSE, ANALYS | P0 | E/D | SDOF on rigid massless disk; F.2 equations fed with the code's own impedances | identity with F.2 solution; worked example (R = 10 m, Vs = 100 m/s, ρ = 2 t/m³, ν = 1/3, h = 10 m, m = 942.48 t, T = 0.5 s): Veletsos–Meek T̃/T = 1.15805 with static springs (K_x = 9.60e5 kN/m, K_θ = 8.00e7 kN·m); with Veletsos–Verbic impedances resonance 1.717 Hz, peak \|ü/ü_g\| = 6.71 | 1e-8 (F.2 identity); ±5 % vs VV-based values | R2 F.1, F.2 |
| VP-18 | Embedment kinematic notches (NUREG/CR-6896 soils) | SITE, ANALYS (FV), MOTION | P1 | E/D | cylinder DOB 11.5–46 m in columns A, B, C, bedrock at 80 m | notches near V/(4·DOB): A 5.4, 2.7, 1.8, 1.4 Hz; B 21.7, 10.9, 7.3, 5.4 Hz; C 13.5, 7.5, 5.3, 4.2 Hz (DOB 11.5, 23, 34.5, 46 m) | ±5 % | R2 E2 |
| VP-19 | FV vs FFV vs FI-EVBN vs FI-FSIN on an embedded box | INTGEN, HOUSE, ANALYS | P1 | D | Mertz-type 120×120×30 ft excavation, uniform 6 ft bricks | SM deviates at/above f_EV (excavated volume fixed at interaction nodes; ≈ 10 Hz in Mertz 2011); MSM closer; FFV skip 1 within 5 % of FV below cut-off | FFV ≤ 5 %; anomaly located within ±10 % of computed f_EV | R2 E3, R1 §4.5 |
| VP-20 | FV vs ANSYS full-harmonic soil box (external) | ANALYS, ANSYS export | P1 | D | Anderson & Ostadan (2014) coarse model: Vs = 2500 fps, Vp = 4677 fps, 5 %, foundation −140 ft, bedrock −240 ft, 600×600 ft box | ATFs agree | 5 % | R2 E4 |
| VP-21 | Kinematic interaction approximations | MOTION | P1 | A | embedded rigid foundation | `H_u = cos(Dω/Vs)` (Dω/Vs < 1.1), `H_θL/u_g = 0.26[1 − cos(Dω/Vs)]` | qualitative ±20 % | R2 E6 |
| VP-22 | Restart equivalence | ANALYS modes 1/2/3, COOX/COOTK | P0 | E | Mode 1 with save; Mode 2 with identical FILE4; Mode 3 with identical FILE1; Mode 3 on a frequency subset | FILE8 identical to the initiation run | 1e-10 rel. | 05c C3 |
| VP-23 | Simultaneous cases | ANALYS | P0 | E | simul = 1 with FILE1X/Y/Z; vibration with 3 load cases | FILE8X equals a single run with FILE1 = FILE1X (Y, Z alike); FILE8001…3 equal single runs | 1e-12 | 05c C4 |
| VP-24 | Coordinate transformation angle | ANALYS | P0 | E | axisymmetric model, x′ input at angle a | `ATF_x = cos a·ATF(a=0)`, `ATF_y = sin a·ATF(a=0)` | 1e-8 | 05c C5 |
| VP-25 | COMBIN merge | COMBIN | P0 | E | ANALYS on odd and even frequency numbers separately | merged FILE8 equals the full run; duplicate frequency → error | exact | 02 Part 6 |
| VP-26 | FFL vs FFM for a rigid massless surface foundation on rock | SITE, POINT, HOUSE, ANALYS (incoherent input) | P1 | E | 20 m × 20 m rigid massless mat (stiff massless links; rigid limit by Richardson extrapolation in the link stiffness) on uniform rock (Vs = 2000 m/s); X/Y/Z input; Abrahamson 2007 hard-rock coherency (isotropic, V_app = 1e9), deterministic AS factors | FFL (`s ⊙ X U′`) and FFM (`X (s ⊙ U′)`) differ by the commutator [S, X]; `Tᵀ[S, X]T` is antisymmetric (D-W3-05), so they are identical for: the generalised force (K_G q)_d in the input direction; the vertical response of the doubly symmetric mat (decoupled by K_G); coherent factors s = 1 (both equal the coherent analysis). The exact relation `K_G(q_FFL − q_FFM) = Tᵀ[S, X]T t_d` is evaluated independently from FILE11 and FILE77. The horizontal ATF differences (≈ 1e-3 X, 3e-2 Y) and the rotations are informative | 1e-8 (identities, rigid limit); 1e-6 (commutator relation); 1e-12 (s = 1) | 05c C6; requirements 4.6 item 5; D-W3-05 |
| VP-27 | Incoherency decomposition | HOUSE incoherency, ANALYS | P1 | E/D | surface grids of interaction nodes; Luco-Wong (model 1) and user tables (model 7); 20 stochastic samples; Line D at 30°, α = 0.1; wave passage V_app = 400 m/s | `\|Σλ − N\|/N < 1e-8` (Eq. 6.1; FILE77 and an independent eigvalsh); f → 0 (0.0122 Hz): λ1/N = 100 % (0.001 % abs.), AS factors \|s − 1\| ≤ 0.01 and ∝ f (ratio 2 at 2f, 2 %), zero-frequency ATF 1 within 5 %; ensemble mean of s sᴴ → Σ (mean squared error / exact expectation within 30 %; error ratio 2 between Ns = 5 and 20, 25 %); directional distances √(2α), √(2(1−α)), ratio 1/3 at α = 0.1 (1e-9); wave-passage phase −ωL/V_app (1e-12) with \|s\| = 1; model 7 bilinear tables vs the tabulated function (λ1/N within 0.05 % abs.) | as stated | 05b §11; requirements 4.4 item 6; D-INC-03, D-W3-06 |
| VP-I1 | Rigid surface foundation under incoherent vertical SH (Luco-Wong, AS) | SITE, POINT, HOUSE, ANALYS | P1 | E/D | 10 m × 10 m rigid massless mat (5 × 5 interaction nodes, stiff links, Richardson extrapolation) on a uniform half-space (Vs = 250 m/s); vertical SH; Luco-Wong γ = 0.5, Vs = 250 m/s; deterministic AS factors (D-INC-03); FFL and FFM | FILE77 = independent coherency matrix and eigen-decomposition (1e-10 abs.); centre ATF_y (rigid limit) = the rigid-foundation average `K_G⁻¹ Tᵀ b` evaluated from FILE3, FILE1 and FILE77, for FFL and FFM; \|ATF_y\| → 1 as f → 0 (1e-4); \|ATF_y\| decreasing with frequency; raw results converge as 1/k (ratio 2, 5 %) | 1e-8 (ATF) | Luco & Wong (1986); requirements 4.4 item 6, 4.6 item 5 |
| VP-I2 | Stochastic incoherency simulation: mean of s sᴴ → coherency matrix | HOUSE | P1 | E | 4 × 4 surface grid (4 m), Luco-Wong γ = 0.3, Vs = 200 m/s, 10 frequencies; HSeed 11, VSeed 17, RandPhz 180°, 50 samples (FILE77001 … FILE77050) | the mean squared error of the Ns-sample mean of s sᴴ equals the exact estimator variance `Σ_ij Var_ij / Ns`, `Var_ij = 1 − Σ_k λ_k² φ_ik² φ_jk²` (Ns = 2, 5, 10, 20, 50); log-log slope of the RMS error −0.5; samples distinct; a re-run with the same seeds identical | 25 % (error ratio); ±0.1 (slope); 0 (re-run) | requirements 4.4 item 6; D-INC-07; 05b §11 test 3 |
| VP-28 | TF interpolation exactness | MOTION, STRESS interpolation | P0 | E | exact 2-DOF hysteretic TF sampled at SSI frequencies | options 0–5 reproduce it to 1e-10 (R1 V5: 1e-14); hysteretic SDOF windows rank-deficient but exact; option 6 exact on cubic polynomial data and O(Δf⁴) otherwise; TFI = TFU at computed frequencies for every option and S; ISRS from interpolated vs dense TF within 2 % | as stated | R1 V5, 05c C8/C9/C14 |
| VP-29 | Identity convolution | MOTION | P0 | E | H ≡ 1 up to Nyquist | output equals input | 1e-12 | 05c C7 |
| VP-30 | Response spectra (Nigam–Jennings) | MOTION RS, EQUAKE RS | P0 | E | step, impulse, harmonic inputs | step: PSA/a0 = 1 + e^{−ζπ/√(1−ζ²)} = 1.854468 (ζ = 5 %); impulse PSV/ΔV = 0.92669; harmonic at resonance PSA → A/(2ζ) = 10.000 (5 %), 25.000 (2 %); max absolute acceleration A√(1+4ζ²)/(2ζ) = 10.0499; SA(f → ∞) → PGA | 1e-5 rel. (Δt ≤ T/100); harmonic 1e-3 | R2 G.1 |
| VP-31 | Baseline correction | MOTION | P0 | E | input with a constant acceleration offset | corrected velocity trend ≈ 0; final displacement ≈ 0 | 1e-6 rel. | 05c C11 |
| VP-32 | SRSS TF | MOTION | P1 | E | one modal FILE8 equal to the coherent FILE8 | phase option 1 → coherent TF; option 0 → \|H_coh\| | exact | 05c C12 |
| VP-33 | Phase adjustment | MOTION | P1 | E | option 6 + pzadj 1; option 2, S = 1000 | option 6 gives \|H\| (zero phase) exactly; option 2 phases scaled by 1/(1+S) | exact / 1e-12 | 05c C13 |
| VP-34 | Relative displacement | RELDISP | P0 | E | reference = node itself; free-field unit reference | zero; equals irfft((H−1)·U_g) and double integration of the acceleration difference | 1e-10; 1e-3 vs integration | 05d §8 |
| VP-35 | Stress recovery | STRESS | P0 | E | rigid-body nodal TF; pure shear and uniaxial states; cantilever beam at low frequency | rigid body → zero STF; τ_oct = (√6/3)τ = 0.8165τ; uniaxial (√2/3)σ; beam end forces in equilibrium (F_I + F_J = 0 unloaded) and M3 = k_L u | 1e-10 | 05d §8 |
| VP-36 | Spectrum-compatible motion (RG 1.60) | EQUAKE | P0 | P | RG 1.60 horizontal 5 % at 1 g; Δt 0.005 s, 20 s | target points reproduced (A 1.0 at 33 Hz, B 2.61 at 9 Hz, C 3.13 at 2.5 Hz, D 2.05 (displ.) at 0.25 Hz → SA 0.4716 g; vertical C 2.98 at 3.5 Hz, D 1.37 → 0.3152 g); generated motion passes: ≥ 100 pts/decade, none > 10 % below, ≤ 9 adjacent below, none > 30 % above, duration ≥ 20 s, strong motion ≥ 6 s; three components \|ρ\| ≤ 0.16; PSD ≥ 80 % of the RG 1.60 App. A target over 0.3–24 Hz | criteria pass | R2 G.2, G.3 |
| VP-37 | Element closed forms | HOUSE elements | P0 | E | cantilever, plates, bars | Euler–Bernoulli β_nL = 1.875104, 4.694091, 7.854757, 10.995541, 14.137168; tip `PL³/(3EI)` (+ `PL/(G·As)` with shear area); released end: lateral 12EI/L³ → 3EI/L³; KJ/KI same component → Error 10; SS Kirchhoff plate `f_mn = (π/2)[(m/a)²+(n/b)²]√(D/ρt)`; bar `f_n = (2n−1)√(E/ρ)/(4L)`; patch test exact constant stress | 1 % (modes, converged mesh); 1e-10 (patch) | R2 I.1, 08 §11 |
| VP-38 | NAFEMS free vibration | HOUSE SHELL, SOLID (TSHELL Test 21: VP-38T) | P0 | P | E = 200 GPa, ν = 0.3, ρ = 8000 kg/m³ | FV12 free plate 10×10×0.05 m: 1.622, 2.360, 2.922, 4.233, 4.233, 7.416 Hz; FV16 cantilever plate: 0.421, 1.029, 2.582, 3.306, 3.753, 6.555 Hz; FV32 membrane (SHELL in-plane) [VERIFY geometry]: 44.623, 130.03, 162.70, 246.05, 379.90, 391.44 Hz; Test 21 thick plate (TSHELL, P2): 45.897, 109.44 (×2), 167.89, 204.51 (×2), 256.50 (×2) Hz | ≤ 2 % first 4 modes (NAFEMS mesh) | R2 I.2 |
| VP-38T | NAFEMS Test 21 thick plate and thin-limit FV12/FV16 (TSHELL) | HOUSE TSHELL | P2 | P | simply supported plate 10 × 10 × 1 m (E = 200 GPa, ν = 0.3, ρ = 8000 kg/m³; hard supports), EINT 0 and 1; quadrilaterals 8 × 8, 16 × 16, 32 × 32 and two triangle patterns on 32 × 32; FV12/FV16 (t = 0.05 m) on 32 × 32 | 45.897, 109.44 (×2), 167.89, 204.51 (×2), 256.50 (×2) Hz (the Mindlin closed form with κ = π²/12 and rotary inertia reproduces them to 2e-4); observed convergence order 2; FV12/FV16 NAFEMS values and within 0.5 % of the Kirchhoff SHELL (DKQ) on the same mesh; the NAFEMS 8 × 8 rows are informative (D-W3-01: −5.4 % on mode 4 is the O(h²) error of the lumped mass) | 2 % (first four modes on 16 × 16, all eight on 32 × 32); order 2 ± 15 %; 0.5 % vs DKQ | R2 I.2; D-ELM-08, D-W1-05, D-W3-01 |
| VP-39 | Spring–mass SDOF and GENERAL equivalence | HOUSE SPRING, GENERAL | P0 | E | `SC,1,1000,0,0,0,0,0,0.05`, MT = 1 (mass units); GM with `MXR = k(1−2β²)`, `MXI = 2kβ√(1−β²)` in the spring's 12×12 pattern | f = 5.033 Hz; ATF peak `1/(2β√(1−β²))` = 10.0125; GM identical (2-node global and 3-node rotated local) | 1e-8 | 08 §11 (corrected C-28) |
| VP-40 | Moving load phase | FORCE, ANALYS, MOTION | P0 | E | factor 0.5, arrival 0.1 s | response = 0.5·f_ref(t − 0.1) response (circular shift) | 1e-10 | 05b §11 |
| VP-41 | Low-frequency ATF ≈ 1 regression | ANALYS, MOTION | P0 | E | any coherent model, X input | \|ATF_x\| at the first SSI frequency ≈ 1 at all nodes | 5 % warning threshold (G-19) | 03 §7 |
| VP-42 | 2D rigid strip on a layer over a rigid base: static K_x and K_φ | SITE, POINT (POINT2), HOUSE (PLANE/BEAMS), FORCE, ANALYS (2D) | P1 | A | rigid massless strip of half-width B = 1 m (stiff links) on a uniform layer (ν = 0.30) of thickness H on a rigid base, B/H = 1/2, 1/4, 1/8 (and 1/16, 1/32), a0 = 0.01; 8/16/32 interaction elements across the strip; three FORCE load cases; R0 = h | Jakub & Roesset (1977): `K_x/G = 1.175(1 + 2.15B/H)`, `K_φ/(GB²) = 2.394(1 + 0.17B/H)` (1/8 ≤ B/H ≤ 1/2) for the Richardson value 2K(32) − K(16) (O(h) convergence: mesh-difference ratio 2 ± 15 %); FORCE + ANALYS = the direct rigid-strip impedance `Tᵀ X_ff T` (1e-4); K_x decreases as H grows, `d(G/K_x)/d ln H` → (1 − ν)/π (5 %); finest-mesh values and half-plane punch solutions informative | ±5 % | R2 B.4; D-PNT-01, D-W3-13 |
| VP-T1 | POINT2 far field vs Kausel's exact line-load Green functions | SITE, POINT (POINT2) | P1 | E/D | R1 V7 column, f = 4 Hz, R0 = 1 m, load at interface 4, observation at the surface | the exact discrete-layer line-load Green functions (Kausel 1981, Eq. 47-48; independent modal sum `greens_tlm.line_load`); parity (u_z \| P_x and u_x \| P_z odd in x, the others even; 1e-12); the core with the prototype's consistent mass = the printed R1 V7 values (1e-9) | 0.33 % at \|x\| = 2R0 (D-W3-12), 0.3 % beyond | R1 2.6, 3.5, V7; Kausel (1981); D-W3-12 |
| VP-T2 | SYMM half and quarter models give the transfer functions of the full model | SITE, POINT, HOUSE, ANALYS | P1 | E | doubly symmetric flexible surface mat with a stick, and embedded FV basement (SOLID excavation, SHELL walls and slabs) with a stick; vertical SV (X) and P (Z); full, half (x ≥ 0) and quarter (x, y ≥ 0) models; X input antisymmetric about x = 0 (type 1) and symmetric about y = 0 (type 0), Z symmetric about both | the reduced models reproduce the full model at every DOF; the DOFs the symmetry constrains vanish in the full model; FILE4 carries the planes and ANALYS lists the image terms `F_red(i,j) = Σ_S s_S F(i, j_S) P_S` (D-W3-13) | 1e-8 | D-ANL-12, D-W3-13; spec 07 9.2.39 |
| VP-T3 | 2D zero-SSI identity: embedded PLANE excavation with structure = soil | SITE, POINT, HOUSE, ANALYS (2D) | P1 | E | plane-strain excavation 6 m wide through the two upper layers (PLANE, 4 columns, FV) whose structure is PLANE elements with the layer properties; vertical SV (X case), P (Z case) and inclined SV (30°) | U = U′_f at every node (UX, UZ) and frequency; control-point ATF = 1; 2D FILE4 `int_dofs` = [1, 3]; `<simul>` = 1 writes FILE8X and FILE8Z only | 1e-8 | R1 4.4 (built-in check); 02 Part 6 |
| VP-43 | Thin-layer monotone convergence and Lamb pulse | SITE, POINT | P1 | E/D | stratum with N = 4, 6, 12 sublayers; inverse FFT of the surface vertical Green function | compliance increases toward the exact value as N doubles; Lamb (ν = 1/4) check values w·πμr/P = −0.018957 (τ = 0.6), −0.010234 (0.8), −0.028501 (0.9), +0.375 (τ ≥ 1.0877) [VERIFY sign] | monotone; 2 % | R2 D.4, D.5 |
| VP-44 | Gravity trick | ANALYS, MOTION/STRESS | P2 | E | fixed-base model, 1 g one-period harmonic, N = 32,768, Δt = 0.005 s (T = 163.84 s) | quasi-static displacement equals a static 1 g solution at T/4 | 1 % | 01 §13 |
| VP-45 | Option NON element models | NONLINEAR | P2 | E | elastic-perfectly-plastic BBC cycled at x = 1.25, 1.5, 2, 5, 10, 50 x_y (GMR/Masing); the module on a relative-displacement history; the state machine | closed loop and peak forces ±F_y (1e-12); `ξ_h = 2(x − x_y)/(πx)`; `K_sec = F_y/x`; `F_μ = K_el max\|x\|/F` = 3; state machine (.NON absent → elastic run + file; = 1 → iteration; the properties of `*_EQL_Matl_Prop.txt` used; deleting both restarts) | 1e-6 | 05d §8 items 7-8; D-NON-03, D-NON-04 |
| VP-46 | SHEAR / BBCGEN | SHEAR, BBCGEN | P2 | D | h_W = 20 ft, l_W = 30 ft, t_W = 2 ft, f′c = 5 ksi, f_y = 60 ksi, ρ = 0.005, N_U = 1000 kips, A_BE = 0; the SI version 6.096 m × 9.144 m × 0.6096 m | V_ACI 4424.82; 10√f′c·A_W 6109.40; V_Wood raw 648.00; Wood lower bound 3665.64; V_Barda 4026.77; V_GW 2405.90 (2617.54 with A_BE = 0.1 ft², f_y,BE = 60 ksi) kips; SI version 19682.6, 17912.0, 10702.0 kN (printed to 0.1 kN); BBCGEN 22 points, yield index 21, point 22 = (0.02, 1.02V_u), first slope G·A_W | 1e-6 rel. (British; SI path vs British × 4.44822162 kN/kip); 0.05 kN abs. vs the printed SI values (D-W3-08); BBCGEN 1e-9 | 10 §3.11, 11 §7; D-NON-09, D-W3-08 |
| VP-NON1 | Option NON SDOF: nonlinear spring iterated through the whole chain | HOUSE, ANALYS, MOTION, RELDISP, NONLINEAR | P2 | E/D | 1 t mass on a GMR spring (k_el = 1000 kN/m; bilinear BBC: yield 10 kN at 0.01 m, 10 % hardening; 2 % elastic damping) on one interaction node of a practically rigid site; harmonic control acceleration 2.5 g at 8.008 Hz with an integer number of cycles in the NFFT window; EDF = 1; New Structure restarts | every SSI response = the closed-form steady state of the equivalent-linear SDOF; amplitude and k of every iteration = an independent closed-form fixed-point iteration (same number of iterations, 5); converged k and ξ = backbone secant and Masing loop damping at the converged amplitude; convergence per D-NON-06 | 1e-6 | requirements 4.15; D-NON-04, D-NON-06 |
| VP-N1 | Near-field soil iterations reproduce SOIL (SHAKE) for a uniform column | SOIL, SITE, POINT, HOUSE, ANALYS, STRESS | P1 | D | laterally uniform column of nonlinear SOLIDs (ETYPE 1, low-strain material = G_max) in an FV excavation of a uniform sand deposit on rock (10 layers of 1 m), vertical SV; free field and excavated soil at the SOIL strain-compatible properties of the same column, motion and curves; (a) start from the low-strain properties, (b) start from the free field (GFAC = DFAC = 1) | (a) monotone convergence within the D-NLS-02 limit to the SHAKE G/G_max and damping profile; (b) converged at iteration 0 with the same profile | 5 % (G/G_max), 10 % (damping) per layer | requirements 4.4 item 7, 4.10 item 5; R1 §6; D-NLS-02, D-NLS-03 |
| VP-N2 | Restart consistency of the near-field soil iterations | HOUSE, ANALYS, STRESS | P1 | E | (a) an iteration whose properties do not change (flat curves); (b) HOUSE repeated with the same FILE74 strains | (a) the New Structure restart (`<mode>` 1, X_ff from COOXqqq) reproduces the initiation FILE8; (b) the same FILE78 properties, FILE8 and FILE74 | 1e-10 (FILE8, FILE74); 1e-12 (properties) | requirements 2.5 (VP-22), 4.4 item 7; D-NLS-01 |
| VP-TS1 | TSHELL patch tests, shear locking, rigid-body modes | HOUSE TSHELL | P2 | E | distorted patches (quadrilaterals and triangles, EINT 0 and 1, oblique planes); cantilever plates t/L = 1/1000 and 1/10; simply supported plates t/a = 1/1000 and 1/10 under uniform load with two triangulations; single elements | membrane, bending and transverse-shear patch tests exact; tip deflection = Kirchhoff / Timoshenko; centre deflection = Navier at O(h²); exactly six zero-energy modes and no near-zero triangle mode; `K* = c(β) K₀`; the plain-MITC3 and locking-element rows informative | 1e-10 (patches); 1 % (8 × 2 cantilever), 0.1 % (32 × 4); 2 % (16 × 16) and 0.5 % (32 × 32) plates; order 2 ± 15 % | requirements 1.5, 4.1 TSHELL; D-ELM-08; MacNeal & Harder (1985) |
| VP-54T | THSHLSTR TSHELL face stresses | HOUSE, STRESS | P2 | E | TSHELL plate under a pure membrane and a pure bending state (synthetic FILE8), EINT 0 and 1 | STF NXX … MXY = D ε, D κ; maxima = closed form × max\|u_g\|; face stresses N/t on both faces and ±6M/t² for the four sign permutations, TXY (++), principal stresses and plane-stress strains (D-TSH-01) | 1e-9 | spec 11 §4.1, §7 item 7; spec 05d test 6; D-TSH-01, D-W3-02 |
| VP-O1 | HOUSE node-numbering optimizer | HOUSE, ANALYS | P1 | E | embedded TSHELL basement on a layered site, deliberately badly numbered, with and without HOUSEX,1 | FILE8 mapped back through `<model>.map` identical to the reference; same number of equations; the map a bijection; bandwidth and profile reduced (profile 0.59 × the original, informative); interaction nodes ascending in their bottom-up order; `.hounew` reproduces the optimised FILE4 (1e-14) | 1e-10 | requirements 4.4 item 5; D-HOU-03, D-FIL-09; spec 05b test 8 |
| VP-LA1 | LOADGEN equivalent static loads | SITE, POINT, HOUSE, ANALYS, MOTION, RELDISP, STRESS, LOADGEN | P2 | E/D | three-mass BEAMS stick (2, 2, 1.5 t) on a rigid 4 m × 4 m surface mat (9 interaction nodes), Vs = 1000 m/s, earthquake-like record (PGA 0.3 g); LOADGEN Disp. and Accel., lumped masses generated from COOSM, critical time = largest base shear | sum of the APDL inertia forces = − STRESS base shear FYI, their moment about the base = − STRESS base moment MZI at the critical time; critical time = time of the largest \|FYI\| (±1 sample); generated masses = MT masses; interface D = RELDISP `.THD` started at rest; source FILE8 = source RESULTS | 1 % (resultants; observed 0.04 %); 1e-12 (masses); 1e-10 (D); 1e-9 (sources) | manual 6.4.15; spec 04 §15.6; D-W3-11 |
| VP-LA2 | LOADGEN dynamic loads and replay of the second step | HOUSE, MOTION, RELDISP, LOADGEN | P2 | E | the same stick on a 3 × 3 mat on frequency-independent soil springs (FILE8 exact at every Fourier frequency); LOADGENDYN with the REL and ACC methods and Rayleigh damping | TABLE time columns and ALPHAD/BETAD exact; ACEL table = control motion (REL) or the mat-centre `.ACC` (ACC); check tables = MOTION `.ACC`; D tables = RELDISP `.THD` − `.THD`(0); source FILE8 = RESULTS; the APDL parses; the exported second step, solved with the HOUSE matrices at every Fourier frequency, reproduces the SSI absolute accelerations | 1e-10 (tables); 1e-9 (sources); 1e-6 (replay; observed 5e-9) | manual 6.4.15; spec 04 §15.6; D-W3-11 |
| VP-SN1 | SOIL-NON linear limit equals the linear SOIL (SHAKE) solution | SOIL | P2 | E/D | γ_r = 1e6 % (linear backbone), same profile and motion as a SOIL run with 0 iterations; elastic base (Joyner-Chen) with outcrop input and rigid base with within input; damping type 1 (frequency independent); Rayleigh and visco-elastic damping informative | surface PGA and 5 % SA at T = 0.05 … 2 s of SOIL-NON = SOIL-EQL; small-strain amplification peak and its frequency of the discrete column = the continuum (1 %) | 2 % (observed ≤ 0.4 %) | R1 §6; Joyner & Chen (1975); Phillips & Hashash (2009); D-W3-10 |
| VP-SN2 | SOIL-NON single-element Masing loop | SOIL | P2 | E | one MKZ element cycled at γ_a = 0.1, 1, 10 γ_r for (β, s) = (1, 1), (1.3, 0.85), (1, 0.7) | loop secant modulus = backbone secant; loop-area damping = the analytical Masing damping (closed form for s = 1, R2 H.2; hypergeometric integral otherwise); the FILE88 damping function = analytical | 1e-3 (loops); 1e-6 (FILE88 damping) | Masing (1926); Kramer (1996) 6.4.3; R2 H.2; D-W3-10 |
| VP-SN3 | SOIL-NON vs SOIL-EQL at moderate shaking (informative) | SOIL | P2 | D | SHAKE91 sample profile and motion (0.1 g rock outcrop); MKZ fitted to the same curves; SOIL-EQL with 8 iterations | both runs complete; surface and depth PGA, maximum strains and spectra compared, with the shear-strain index I_γ = PGV/V_S30 of Kim et al. (2016) (observed 0.04 %; surface PGA 12 % lower with SOIL-NON) | informative (no criterion) | R2 H.1; Kim et al. (2016); D-W3-10 |
| VP-54W | Water modelling and refinement (REFINEMODEL, FILLPOOL, LISTPOOLINTER, MERGEPOOL) | UI (model generation), HOUSE material | P2 | E | REFINEMODEL of one quadrilateral, a 2 × 2 mesh and a hexahedron (+ prism); rectangular SHELL pool 4 × 3 m, 3 m walls, 0.5 m levels, EmptyLevels 2; a SOLID pool; a CUT2SUB pool of a building | 1 quad → 4 (5 new nodes), 2 × 2 → 16 quads / 25 nodes, hexahedron → 8 (19 new nodes), area and volume preserved, attributes kept, triangles and prisms unchanged, hanging nodes reported; water volume a·b·(H − m·h); one spring per coincident wall/water pair, Stiff along the normals and stiff2 along the walls; offset error with the model unchanged; water K = 2.2 GPa and G/K = 1e-8 (D-W3-07), mass ρV; area shells = wetted area; MOPT,0; WRITE → INP; LISTPOOLINTER finds every interface node; MERGEPOOL keeps the node numbers and refuses a second import | exact counts; 1e-12 (areas, volumes, masses); 1e-8 (K) | requirements 6.4 VP-54; spec 11 §1, §7 items 1-2; D-WAT-01, D-W3-07 |
| VP-W1 | Hydrodynamic (impulsive) mass of a rigid rectangular tank filled by FILLPOOL | UI (FILLPOOL), SITE, POINT, HOUSE, ANALYS | P2 | E/A | rigid tank 10 m (direction of motion) × 2 m, walls 6 m, water depth H = 5 m (0.5 m levels), on a practically rigid site (Vs = 10 km/s); SHELL walls, and SOLID walls on a SOLID slab with interface-area shells; harmonic SV (X) input at 2, 4 and 8 Hz | water mass in the HOUSE mass matrix = ρ·2L·2B·H (1e-10); base shear = inertia of the water (1e-6); the tank follows the control motion, \|H − 1\| ≤ 1e-3; F/(m a) = the exact potential-flow impulsive mass ratio m_i/m = 0.5 for L = H (observed 0.507-0.509); Housner (0.542), Westergaard (0.583) and the resultant height (0.4047 H) informative | 5 % (impulsive mass) | Westergaard (1933); Housner (1963); exact potential-flow series; D-WAT-01, D-W3-07 |

### 6.4 Utilities, post-processing and generation tests

| ID | Test | Expected | Tier |
|---|---|---|---|
| VP-47 | Line maths T-L1: A x=[0,1,2], y=[0,1,2]; B x=[0.5,1.5], y=[10,20] | grid [0,0.5,1,1.5,2]; ADDITION [10,10.5,16,21.5,22]; AVERAGE [5,5.25,8,10.75,11]; SUBTRACTION A−B [−10,−9.5,−14,−18.5,−18]; LINECOMBIN 2A+0.5B [5,6,9.5,13,14]; SRSS [10,10.01249,15.03330,20.05617,20.09975]; spec 06 case L1 = {(1,1),(3,3)}, L2 = {(2,10),(4,20)} → ADDITION {11,12,18,23} | P1 |
| VP-48 | BROADEN T-B1: x=[1,5,9,10,11,15,20], y=[1,1,1,5,1,1,1], `BROADEN,2,0,15,1` | B = 5 on [8.5, 11.5]; ramps (7.65,1)→(8.5,5) and (11.5,5)→(12.65,1); B = 1 elsewhere. T-B2: x=[1,4,5,6,6.5,7,8,9,12,20], y=[0.5,2,4,3,3.5,3.0,3.8,2,1,0.5], `BROADEN,2,15,0,1` → flat 3.8 from 5.2 to 8. T-B3: two-line envelope = pointwise max | P1 |
| VP-49 | Section cut T-S1: two unit hexes x∈[0,1],[1,2], plane z = 0.5, n = (0,0,1), r = (1,0,0) | A = 2, C = (1,0.5,0.5), Ixx = 0.166667, Iyy = 0.666667, Ixy = 0; Szz = +10/−10 → My = +10, Fz = 0; + Sxz = 3 → Fx = 6; Syz = +1/−1 → Mz = −1. T-S2: two shells y = 0, 2×2 each, t = 0.5, plane z = 1: A = 2, Iyy = 2.6667; Sy′y′ = ±100 → My = +200; Sy′y′ = 100 both → Fz = 200 | P1 |
| VP-50 | INTGEN counts on an n_x × n_y × n_z box of hexes | FV (n_x+1)(n_y+1)(n_z+1); EVBN FV − (n_x−1)(n_y−1)(n_z−1); FSIN EVBN − (n_x−1)(n_y−1); FFV EVBN + (n_x−1)(n_y−1)·(selected internal levels); surface (n_x+1)(n_y+1); options 1–5 additive; 0 clears | P1 |
| VP-51 | Generation tools | MERGE offsets/translation; MERGESOIL Mode 1 (n coincident nodes collapse to the lower number, map has every excavation node), Mode 2 (n springs, 1e7), Mode 3 split at SepLevel (1e7 / 10); EXCAV 4×4 basemat, 3 levels, jitter ±0.001: delta 0 → extra levels, delta 0.01 → 48 hexes and 75 nodes; WELD two coincident cubes → 12 connected nodes, 4 unused; NCOM after RMVUNUSED remaps masses, BCs, interaction flags; ROTATE (1,0,0) by rxy = 90 → (0,1,0); FIXROT flat XY plate → ROTZ fixed only; rotated 30° about X → springs | P1 |
| VP-52 | Model checks | planted faults: shared interior excavation node → EXCSTRCHK; fixed interaction node → FIXEDINT and Error 124 (all translations fixed); beam meeting solids at one node → HINGED; K node interaction → KINT; free spring node → FREESPRING; oblique Kirchhoff shell without FIXROT → EDU-06 | P0 |
| VP-53 | Frame and file formats | WRITESPEC/READSPEC and WRITETH/READTH round trips exact; READSPEC column counting (1 + 3 columns, numLines = 3); FRAMECOMBIN SRSS/sum/average; CRITFREQ flags a 40 % TFI overshoot at tol 20, not at 50 | P1 |
| VP-54 | Water, refinement, databases, thick shell | REFINEMODEL: 1 quad → 4 (5 new nodes), 2×2 mesh → 16 quads / 25 nodes, hex → 8 (19 new nodes), area/volume preserved; FILLPOOL water volume a·b·(H − m·h); binary COMB*DB = algebraic sum; BINSTRTBL step −1 = signed abs-max; THSHLSTR pure bending ±6M/t²; implemented as VP-54T (thick shell) and VP-54W (water, refinement) in §6.3, the binary-database part is not implemented | P2 |

### 6.5 Module coverage matrix

| Module / area | VPs and UTs |
|---|---|
| Interpreter, model data, CHECK, AFWRITE | UT-01…UT-19, VP-52 |
| EQUAKE | VP-30, VP-36, UT-16 |
| SOIL | VP-02a, VP-04 (SOIL-NON P2: VP-SN1, VP-SN2, VP-SN3) |
| SITE | VP-02b, VP-05, VP-06, VP-15, UT-19 |
| POINT | VP-07, VP-08, VP-09, VP-13, VP-42, VP-43, VP-T1 |
| HOUSE elements | VP-01, VP-03, VP-37, VP-38, VP-39, UT-14, UT-15 (P1: VP-O1; P2 TSHELL: VP-38T, VP-TS1, VP-54T) |
| HOUSE / ANALYS SSI | VP-10, VP-11, VP-13, VP-14, VP-15, VP-16, VP-17, VP-22, VP-23, VP-24, VP-41 (P1: VP-12, VP-18, VP-19, VP-20, VP-26, VP-27, VP-I1, VP-I2, VP-42, VP-T2, VP-T3, VP-N1, VP-N2) |
| FORCE | VP-40, VP-23 |
| COMBIN | VP-25 |
| MOTION | VP-01, VP-15, VP-28, VP-29, VP-30, VP-31 (P1: VP-32, VP-33, VP-21) |
| RELDISP | VP-34 |
| STRESS | VP-35, VP-28 (nonlinear-soil strains: VP-N1, VP-N2; TSHELL face stresses: VP-54T) |
| Line maths, plots, cuts, generation (P1) | VP-47, VP-48, VP-49, VP-50, VP-51, VP-53 |
| NONLINEAR (Option NON, P2) | VP-45, VP-46, VP-NON1 |
| LOADGEN (Option A, P2) | VP-LA1, VP-LA2 |
| Water modelling (P2) | VP-54W, VP-W1 |
| P2 | VP-38T, VP-44, VP-45, VP-46, VP-54T, VP-54W, VP-NON1, VP-LA1, VP-LA2, VP-SN1, VP-SN2, VP-SN3, VP-TS1, VP-W1 |

Planned but not implemented as verification problems in this build: VP-12, VP-18 to VP-21 and VP-44 (the registered problems are listed by `python -m sassi.verify.report --list`).

### 6.6 Cross-checks against ANSYS (for the user's own verification; optional)

1. Export each element benchmark (VP-37/38) with `ANSYS` and compare ANSYS modal frequencies with SASSI-EDU
   fixed-base eigenvalues to 0.1 % (isolates converter and element differences from mesh effects).
2. Convert ANSYS `.cdb` models (BEAM188 B ≠ H rectangle, SHELL181, SOLID185, COMBIN14, MASS21) and compare
   HOUSE stiffness against ANSYS static solutions (beam-axis mapping test, D-ANS-03).
3. VP-20 ANSYS full-harmonic soil box with damping mapped by `E_ANSYS = E(1−2β²)`, `g = 2β√(1−β²)/(1−2β²)`.

### 6.7 Reference data still to be acquired (action items)

| Item | Needed for | Action |
|---|---|---|
| Tajirian & Tabatabaie (1985) Figs. 4–11 | VP-13 curves | digitise the open Bechtel/MTR scan |
| Luco–Westmann 1971, Veletsos–Wei 1971, Wong–Luco 1978/1985 tables | VP-11 as published data | obtain papers |
| Abrahamson coherency coefficients (1993, 2005, 2006, 2007 rock/soil) | incoherency models 2–6 | **partly done (wave 3):** models 3 (2005), 5 (2007 hard rock) and 6 (2007 soil) transcribed with citations into `sassi/data/coherency/*.json` and plot-tested; models 2 (1993) and 4 (2006 embedded) still need a confirmed source (they are refused with "coefficients not available") |
| NAFEMS TNSB figures (FV2, FV4, FV32 geometry; FV52 conflicting targets) | VP-38 | obtain TNSB Rev. 3 |
| Pekeris 1955 Fig. 1 sign | VP-43 | confirm |
| Kausel 1974 layer vertical/torsion formulas | VP-13 statics | confirm |

---

## 7. Decisions (resolution of open questions)

Each decision is binding for the implementation. "Source" cites the open question(s) it resolves (spec file
and OQ number). Decisions marked **(gated)** must be confirmed by the named verification problem; if it fails,
the stated fallback applies. All switches named `EDUOPT,<KEY>` are runtime options with the default shown.

### 7.1 General and platform (GEN)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-GEN-01 | Product identity | "SASSI-EDU", Python package `sassi`. Reproduce module names, file names, commands, option codes, dialog labels and CHECK messages. Byte compatibility with vendor binaries or decks is a non-goal | user request (same features and nomenclature); vendor formats are proprietary |
| D-GEN-02 | Language and dependencies | Python ≥ 3.9 (venv 3.9.6); numpy, scipy for numerics; matplotlib (plots), plotly (optional HTML export), tkinter (GUI), pytest (tests); no compiled extensions | installed environment |
| D-GEN-03 | Binary container (01 OQ, 02 OQ7, 05b OQ24) | NumPy `.npz` written to the legacy name without extension; entry `meta.json` holds format name, version, model hash, units, creation module | no h5py; self-describing; inspectable by students |
| D-GEN-04 | Module deck formats (05a OQ18, 05d OQ-S12, OQ-N1) | Keyword-structured text decks with header `SASSI-EDU <MODULE> DECK v1`, one record per line, documented in the IO spec; written only by AFWRITE; read by the module and by CONVERT,SSI | vendor fixed formats unavailable; readable decks help learning |
| D-GEN-05 | GUI toolkit | tkinter with embedded matplotlib (TkAgg); Agg for headless | only toolkit available |
| D-GEN-06 | Geometric tolerance (09 §1.5, 11 OQ25) | `tol = max(1e-6·L_ref, 1e-9)`; layer-interface match `max(tol, 1e-4·h_min)`; `EDUOPT,GEOMTOL,<value>` overrides | scale-independent |
| D-GEN-07 | Size limits (01 §3.2) | `EDUOPT,LIMITS,ACS` (default) warns (EDU-02) when manual V3 limits are exceeded; `PREP` adds the legacy limits (100 top layers, 15 dynamic properties); `UNLIMITED` silences. Hard errors only where physics needs them (power-of-2 NFFT, frequency consistency) | limits are RAM artefacts of the original |
| D-GEN-08 | Random numbers (05b OQ8, 10 Q-27) | numpy PCG64 everywhere (EQUAKE, incoherency, RND); documented seeding | reproducible; streams differ from the vendor's by necessity |
| D-GEN-09 | Determinism | Same input → bit-identical output on one platform (fixed ordering, no nondeterministic parallel reductions) | needed for restart-equivalence tests |
| D-GEN-10 | Performance targets | P0 must handle ≥ 2,000 interaction nodes and ≥ 20,000 DOF per frequency on a laptop; dense X_ff in core; no out-of-core solver | educational scale |

### 7.2 Conventions (CNV)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-CNV-01 | Time factor (02 OQ, 03 OQ18) | `e^{+iωt}`, outgoing `e^{−ikx}`, `H^{(2)}` | SHAKE, Kausel, Waas, SASSI damping sign (R1 §0) |
| D-CNV-02 | FFT convention | numpy rfft/irfft as in §4.0.1 | TFs under e^{+iωt} apply directly |
| D-CNV-03 | Complex-modulus form (01 OQ1, 02 OQ1, 05a OQ11, 08 Q, 09 Q-F2) | `G(1−2β²+2iβ√(1−β²))` default; `CMODFORM,1` selects `G(1+2iβ)` | SHAKE91 code, TLUSH, SASSI2000 SDOF evidence [V] (R1 §1, R2 §0) |
| D-CNV-04 | Application of damping (R1 §1, 08 §6.6) | βs on G, βp on M = λ+2G; beams from (M*, G*); shells require βp = βs; springs k·c(damp); GM as entered; identical factor in SOIL, SITE, POINT, HOUSE | consistency between free field and excavated soil (VP-16) |
| D-CNV-05 | Mass matrices (01 §5.3, R1 §2.1) | ½ lumped + ½ consistent for SOLID, PLANE, SITE/POINT layers and excavated soil; BEAMS consistent; SHELL/TSHELL lumped | manual; 4th-order dispersion accuracy (R2 A.3: 0.55 % at λ/5) |
| D-CNV-06 | TF normalisation and convolution (02 OQ5, 05d OQ-S1) | FILE8 H = U/U_cp dimensionless; stresses, strains and relative displacements use `U_g = −g·A/ω²` (f = 0 term 0); accelerations use `g·A` (output in g) | ratio is the same for displacement and acceleration; units explicit |
| D-CNV-07 | ATF below f_1 and above the cut-off (05c OQ5) | see D-MOT-03 | rigid-body limit; cut-off acts as low-pass |
| D-CNV-08 | Units (01 §10, 07 OQ31) | one consistent model unit system; only weight→mass conversions (÷ g); British if g > 20 (warn if g differs > 5 % from 32.174 or 9.80665) | manual; detection needed for g-unit outputs |
| D-CNV-09 | Rotations at interaction nodes (08 §3.2) | interaction coupling uses translations only; rotational DOFs of shells/beams at interaction nodes stay ordinary structural DOFs | manual rule 5 |
| D-CNV-10 | NFFT not a power of 2 (05a, 07 OQ26) | write the **nearest** power of 2 (Warning 9; tie → up); EDU-03 is raised as an **error** for the modules that read the history (SOIL, MOTION, STRESS, RELDISP) when the rounded value is shorter than the records used, recommending the next power of 2 | manual Warning 9 text, with a safety net |
| D-CNV-11 | Maximum FFT points (01 OQ4) | 32,768 for EQUAKE, SOIL, SITE, HOUSE, FORCE; 65,536 for MOTION, RELDISP, STRESS | manual §1.2 per-module limits |
| D-CNV-12 | Frequency numbers above NFFT/2 | warning (G-02); MOTION/STRESS interpolate only up to `min(f_N, f_Nyq)` | harmonic-only runs may legitimately exceed |

### 7.3 Command language and interpreter (PAR)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-PAR-01 | One code path for GUI and scripts (04 §9) | single interpreter; dialogs emit commands | verifiability |
| D-PAR-02 | Space before first argument (01 OQ26) | allowed for the first separator only (`EDGE 1,0,0,1`) | manual example |
| D-PAR-03 | Comment syntax (07 §1.2) | `*` as first non-blank character; no inline comments | manual examples; `*` may occur in data |
| D-PAR-04 | Abbreviations (07 OQ1, 08 Q1) | exact full name or documented abbreviation; legacy aliases INTL, LMOV, SLIS, DELS, ETYP, STRE; spelling aliases of §3.2; UI commands full name only; full-name match first | old PREP decks load; no ambiguity |
| D-PAR-05 | Blank / missing arguments (07 §1.4) | blank or omitted → documented default; missing required numeric → 0 with a warning | PREP convention |
| D-PAR-06 | Text with commas (07 OQ3, 10 Q-24) | last-position text takes the rest of the line; `"…"` quoting anywhere | titles in practice contain commas |
| D-PAR-07 | Numeric parsing | free format, `E`/`D` exponents, integer fields rounded with warning | Fortran-era decks |
| D-PAR-08 | Lists in NOUT/EOUT (07 OQ2) | integers and `a-b` ranges separated by `,` blank tab `;`; descending range → error | dialog list syntax |
| D-PAR-09 | Error behaviour (04 §4.1) | `<X> Command not found`; failing commands report and processing continues; EOF summary | manual |
| D-PAR-10 | Substitution order (04 §5.3) | `$k$` → `@`/`#` → tokenise → dispatch; FOREACH body re-substituted per iteration | manual examples |
| D-PAR-11 | Recursion and length guards | line ≤ 3000 chars; macro depth 64; INP depth 32 | runaway protection |
| D-PAR-12 | Relative paths (04 OQ5, 07 OQ7) | absolute → model path → CWD → calling file's folder | manual says model path; others are fallbacks |
| D-PAR-13 | Missing macro argument (04 OQ5, 10 §5.6) | empty string and warning | non-fatal like the rest of the language |
| D-PAR-14 | Deleting output requests (07 OQ9) | `NOUT,0`, `EOUT,0`, `RDND,0` clear their lists; WRITE emits the clear first | consistent with FREQ/DAMP |
| D-PAR-15 | `@` inside longer tokens and `#` grammar (04 OQ6) | `@` + longest defined variable name (+ optional `[i]` or operator) is substituted anywhere; unknown name left literally with a warning; `#` substituted only inside FOREACH bodies (innermost loop index; inside `@V[#]` the loop over V); `#+k` evaluates without mutation | supports `Node@X[#].rs` naming |
| D-PAR-16 | `@X[i]` out of range | error; the command is skipped | silent wrong data is worse |
| D-PAR-17 | WRITE order (07 §9.2.47) | canonical order of spec 07 §9.2.47; X-commands only when non-default; byte-stable output | reviewable `.pre`; round trip |

### 7.4 Model data and generation (MDL)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-MDL-01 | Model numbering (04 OQ1–2, 06 OQ20) | 0-based; Model ▸ New and Model ▸ Open use the lowest unused number and activate it; Cut dialog default = active model | screenshots show "Model 0" |
| D-MDL-02 | SAVE/RESUME format (04 OQ3) | `<model>.sdb` (npz container of the full state incl. hide sets, versioned); RESUME refuses unknown major versions | round trip; forward safety |
| D-MDL-03 | Commands without MDL (04 OQ4) | SAVE, RESUME, AFWRITE, AFWRBAT, RUNxxx refuse with an error; ANSYS export writes `unnamed.inp` in the CWD with a warning | the manual's "may run anywhere" is unsafe |
| D-MDL-04 | Coordinate systems (08 Q3, Q7, Q20) | nodes stored in their defining system; converted at AFWRITE/GLOBAL; LOC = intrinsic Rz(txy)·Rx(tyz)·Ry(txz) (ANSYS); LOCAL Y = Z × X; generation commands work in the active system; redefining a system does not move stored local coordinates (warning) | ANSYS-familiar; right-handed |
| D-MDL-05 | "Last defined" (08 Q2) | definition order for nodes and loads; highest index for property tables | matches manual examples |
| D-MDL-06 | D defaults | n2 defaults to n1; val defaults to 0 (free) | manual + lenient |
| D-MDL-07 | INT codes (08 Q5) | 0 interaction (default), 1 intermediate, 2 interface, 3 internal; only code 0 affects the solution | truncated manual text |
| D-MDL-08 | Nodal mass units and MOPT timing (05a OQ20, 08 Q17–18) | MUNITS default 1 (weight); MOPT flags act when MT/MR/F/MM execute; add mode for F/MM adds factors, keeps the latest arrival time, warns | MR text says weight units |
| D-MDL-09 | EGEN, NMOVE, LMOVE (08 Q6, Q8, Q15) | EGEN increments K and copies attributes; NMOVE blank = 1, explicit 0 kept with warning; lists up to 15 nodes | manual syntax line |
| D-MDL-10 | ETYPE 0 and MSET binding (08 Q22) | raw MSET index stored; classification resolved at CHECK/AFWRITE (D-HOU-01) | later ETYPE changes reinterpret correctly |
| D-MDL-11 | Deletion cascades (08 Q21) | NDEL removes loads, masses, BCs and flags of deleted nodes; dangling element references remain for CHECK (Error 41); GDEL drops that group's EOUT requests (warning); ECOMPR remaps EOUT | predictable |
| D-MDL-12 | MERGE / MERGESOIL details (04 OQ8–9, 09 Q-M2, Q-M3) | MERGE: Mdl2 node, group, material, R, SC, MX and coordinate-system ids offset by the Mdl1 maxima (new ids start at max+1); **L table not offset** (far-field layers are shared; mismatching same-number layers warn); analysis options and frequency sets from Mdl1. MERGESOIL: soil-model material m → L layer m (existing L kept if different, warning); coincidence = geometric tol, only z ≤ gelev + tol; Mode 1 keeps the lower number and moves interaction flags to it; springs translational only, damping 0; z = SepLevel counts as below; a trailing non-numeric optional token is the Mapping file | excavation layers must equal far-field layer numbers; lenient example parsing |
| D-MDL-13 | Cuts (04 OQ11–12, 09 Q-C1, Q-C3) | cuts are session-global sets of (group, element); CUTVOL selects elements with all nodes inside the closed box (tol); CUT2SUB copies elements, their nodes, referenced materials/properties/layers, BCs, masses and interaction flags | Cut plot accepts any model |
| D-MDL-14 | ROTATE (04 OQ10, 09 Q-M4) | Rz(rxy), then Rx(ryz), then Ry(rzx), right-hand positive, about (x,y,z); local systems rotate too; warn when global SC or 2-node GM properties exist | argument order |
| D-MDL-15 | WELD (09 §4.9) | keep the lowest number; never collapse SPRING elements unless `WELD,FORCE` | protects FIXROT/MERGESOIL springs |
| D-MDL-16 | GROUP type mismatch | error, use MTYPE; `BEAM` accepted for BEAMS | explicit |
| D-MDL-17 | MACT/RACT initial values (08 Q23) | 1, global | manual silence |
| D-MDL-18 | Option AA flag (04 OQ16) | model attribute `ansys_aa`; AFWRITE writes a topology-only `.hou` marked `OPTION AA PLACEHOLDER`; RUNHOUSE refuses it outside Option AA | manual warning |
| D-MDL-19 | EXCAV and SOILMESH algorithms (09 Q-M1, Q-M5) | as spec 09 §4.1 and §4.7 (template from lowest-level faces; level clustering by delta, level value = cluster mean; nodes numbered bottom-up; MSET from the layer containing the element mid-height); SOILMESH linear growth `1 + i·s/100`, soil ETYPE 1, contact springs with property rNum when contact ≠ 0 | reconstruction documented and testable (VP-51) |

### 7.5 Elements (ELM)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-ELM-01 | SOLID EINT rules (08 Q9) | 0 → 2×2×2, 1 → 3×3×3, 2 → 4×4×4 Gauss | manual default "2×2"; SAP heritage |
| D-ELM-02 | Incompatible modes (01 OQ12) | Wilson–Taylor 9 modes with centroid-Jacobian correction; only for structural SOLIDs and only if MOPT `<incomp>` = 0; condensed from K* | manual |
| D-ELM-03 | PLANE incompatible modes | same MOPT switch for structural PLANE elements; never for excavated soil | 2D/3D consistency |
| D-ELM-04 | Beam theory (01 OQ13, 08 §4.5) | Timoshenko with shear areas (0 → Euler–Bernoulli); releases by static condensation of the local stiffness; the condensed mass uses the same transformation | manual R table |
| D-ELM-05 | Beam mass (08 Q10) | consistent translational mass ρAL; torsional ρ(I2+I3)L/3 (polar from I2+I3); no rotary inertia | standard |
| D-ELM-06 | SHELL formulation (08 Q11) | flat facet: Q4 plane-stress membrane with incompatible modes (CST for triangles) + DKQ/DKT bending; no drilling | Kirchhoff per manual |
| D-ELM-07 | Shell lumped mass | ρtA/n on translations; zero rotational inertia | manual "plates lumped" |
| D-ELM-08 | TSHELL (P2) | MITC4; EINT 0 reduced (hourglass-stabilised), 1 selective; drilling `1e-4 × min membrane diagonal × element area` | manual + stability |
| D-ELM-09 | DOF elimination (08 Q4) | DOFs no attached element defines are removed automatically (info); defined but unstiffened DOFs (SHELL drilling) need FIXROT/FIXSHLROT (EDU-06 warning) | avoids singular systems without hiding modelling issues |
| D-ELM-10 | FIXROT/FIXSHLROT spring representation (01 OQ11, 09 Q-G5) | drilling spring = 2-node GENERAL element between the node and a new fully fixed coincident node, `K = stiff·n nᵀ` on the rotational block (exact about the normal); FIXROT uses D instead when the normal is parallel to a global axis; default stiff 10 (model units); EDU warning if stiff > 0.1 × the smallest bending-rotation diagonal at that node | exact for oblique shells; "≤ 10 % of bending stiffness" made testable |
| D-ELM-11 | Gap and unused nodes | fixed in the decks only (W1/W4); K-node-only nodes are not warned but fixed (D-CHK-08) | manual |
| D-ELM-12 | Excavated-soil element properties | same complex moduli, densities and mixed mass as the SITE layer they represent | required for the FV identity (VP-16) |

### 7.6 SITE and POINT (SIT, PNT)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-SIT-01 | Layer discretisation | Kausel TLM matrices of R1 §2.2 with mixed mass; user layers used as entered (no automatic sub-division; G-05 warns) | interaction nodes must coincide with user interfaces |
| D-SIT-02 **(gated by VP-02b)** | Half-space sublayer law (03 OQ11, 05a OQ12, R1 §10.1) | geometric series, non-decreasing, `h1 = min(h_last·Vs_hs/Vs_last, 1.5λ/n)`, ratio q ≥ 1 by bisection so Σh = 1.5λ, each h ≤ λ/8; uniform fallback when infeasible. `EDUOPT,HSLAW,GEOMETRIC\|UNIFORM\|LINEAR`. If VP-02b fails, default becomes UNIFORM (R1 V4: ≤ 0.7 % at n = 20). No low-frequency depth cap (optional `EDUOPT,HSMAXDEPTH`) | manual "increasing with depth" with an accuracy safeguard |
| D-SIT-03 | Generated layers range (05a OQ1) | 0 or 4–20 (Error 47); EDU warning below 10 | CHECK text vs recommendation |
| D-SIT-04 | Top-layer count (05a OQ2) | ≤ 200; Warning 8 only with `LIMITS,PREP` | V3 limit |
| D-SIT-05 | Frequency 1/2 units (05a OQ3, 07 OQ20) | frequency numbers (Hz shown in the UI); ratios constant outside the range | defaults 1 / 4000 and "≥ NFFT/2" |
| D-SIT-06 | Linear/Non-Linear Soil storage (05a OQ4, 07 OQ19) | `SITEX,<soilmode>`; SITE reads FILE88 at run time when 1 | `<opmode>` keeps its documented meaning |
| D-SIT-07 | Control motion definition (R1 §2.7a) | within (in-column) motion at the top of the control layer; outcrop control motions are converted with SOIL or the control point is put at the surface | SASSI convention; avoids the SHAKE "2E" ambiguity |
| D-SIT-08 | Inclined body waves and surface waves (05a OQ22, R1 §10.3) | vertical body waves P0; inclined body waves and R/L wave fields P1. For θ > 0 the half-space is represented by the exact Kausel–Roesset half-space stiffness at the apparent wavenumber (generated layers ignored, info) | dashpot base is poor for oblique incidence |
| D-SIT-09 | R/L mode selection (05a OQ13) | shortest wavelength = largest Re k among propagating modes (Im k ≤ 0); least decay = smallest \|Im k\|; Love = shortest wavelength | option names |
| D-SIT-10 | 1/5-wavelength rule velocity (03 OQ14) | strain-compatible Vs (FILE88 when soil mode 1); Vp checked only for P-wave input | governing wavelength |
| D-SIT-11 | FILE1X/Y/Z production | SITE writes FILE1; users copy with FCOPY; ANALYS warns if FILE1X is not SV/x′, FILE1Y not SH/y′, FILE1Z not P/z′, or angle ≠ 0 | manual workflow, checked |
| D-PNT-01 | POINT2 vs POINT3 (05a OQ21) | from HOUSE `<dim>` (1 → POINT2, 2 → POINT3); CHECK error if inconsistent | manual |
| D-PNT-02 | FILE3 content (R1 §10.5) | mode amplitudes α and axis values; exact far-field evaluation in ANALYS | no tabulation error |
| D-PNT-03 | Central-zone discretisation (R1 §10.2) | one radial axisymmetric element, vertical mesh = SITE sublayers, mixed mass | only discretisation R0 fully defines; [V] R1 V8 |
| D-PNT-04 | Load interfaces (05a OQ23) | interfaces 1 … L+1 included (L = 0 → interface 1 only) | base of embedment needed |
| D-PNT-05 | RADIUS formula (01 OQ14, 09 Q-G10) | `r_e = Scale·√A_plan` (A_plan = area of the convex hull of the element's plan projection; PLANE: horizontal width); output `group element r_e`, then min, average, max; Scale default 0.9 | reproduces the manual's 0.90h rule |

### 7.7 HOUSE, incoherency and multiple excitation (HOU, INC)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-HOU-01 | "Below ground" test for ETYPE 0 (05b OQ20) | all nodes `z ≤ gelev + tol` and centroid strictly below | robust for elements touching grade |
| D-HOU-02 | `.sit` requirement (02 OQ20) | always required (layer table, units) | manual |
| D-HOU-03 | Node optimizer (05b OQ21) | P1; RCM keeping interaction nodes in original relative order; `.hounew` + `.map` (`old new` per line); FILE4/FILE8 in new numbering; solver ordering independent | nomenclature kept; sparse solver orders internally |
| D-HOU-04 | Excavation numbering/grouping rules R5/R6 | warnings; error only if STRESS soil-pressure or nodal-contour output is requested | manual "strict if contours needed" |
| D-HOU-05 | Interaction-node order | bottom-up ascending: error when incoherent, warning otherwise | manual |
| D-INC-01 | Model-type argument (05b OQ5, 07 OQ23) | `WPASS <cohf>` = unlagged coherency model 1–7; `INCOH <alpha>` = α (models 2–7) or mean Vs (model 1) | dialog grouping |
| D-INC-02 | Superposition mode, simulations (05b OQ6) | HOUSEX `<supmode>`, `<nsim>`; stochastic iff (HSeed ≠ 0 or VSeed ≠ 0) and RandPhz > 0 | no manual argument |
| D-INC-03 | Mode order and sign (05b OQ7, OQ10) | λ descending; sign normalised so `Σ_i φ_ik ≥ 0` ("with ATF phase adjustment"); `EDUOPT,INCOHSIGN,RAW` = "without" | f → 0 reproduces coherent motion |
| D-INC-04 | Abrahamson coefficients (05b OQ12) | from versioned data files `data/coherency/*.json` transcribed from EPRI/NRC sources with plot regression tests; models 2–6 refuse to run while the file is marked unverified | do not use from memory |
| D-INC-05 | Luco–Wong form (05b OQ13) | `exp[−(γ_c ω D/Vs)²]`, ω in rad/s | Luco & Wong 1986 |
| D-INC-06 | Wave-passage reference (05b OQ14) | τ = 0 at the ANALYS control point (x_c, y_c) | same origin as the oblique-wave phase |
| D-INC-07 | Random phases (05b OQ8) | independent per frequency, mode, direction and sample; X and Y from two spawned HSeed streams, Z from VSeed | documented reproducibility |
| D-INC-08 | `<ngp>` / Error 60 (05b OQ9, 07 OQ14) | stored; must be ≥ 1 only when incoherent with embedded interaction nodes; used to validate the number of distinct interaction-node levels; 0 = auto | conflict with dialog value 0 |
| D-INC-09 | Ill-conditioned projections (03 §3 item 4) | default: no merging; EDU-17 warning when projections from different levels are closer than 0.1 × mesh size; `EDUOPT,INCOHMERGE,1` merges unique plan positions | fidelity first; remedy available |
| D-INC-10 | Multiple excitation (05b OQ16–17, 07 OQ8, OQ15) | ME numbers ≥ 1 (warning above 10, ≤ 5000); AMP appends; complex SAR as consecutive (Re, Im) pairs; [0, 10] checked on the modulus; SAR = 1 outside zones; applied after incoherency and wave passage | manual limits reconciled |
| D-INC-11 | User coherency files (05b OQ15) | whitespace text: FREQCOH (NF values), DISTCOH (ND values), COH?USER (NF rows × ND columns); sizes from the files; bilinear in (f, D), clamped, γ(D=0) = 1 | simplest consistent layout |
| D-INC-12 | FFL/FFM argument (05c OQ9) | `ANALYSX,<ffm>`; FFL default; FFM with embedded interaction nodes → EDU warning | manual recommendation |

### 7.8 FORCE, ANALYS, COMBIN (FRC, ANL, CMB)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-FRC-01 | Load phasing (05b OQ25) | `P = a·e^{−iωt0}` per DOF; response `IFFT[H·F_ref]`, F_ref = FFT of the MOTION/STRESS THFILE (force history) | manual "time lag = complex phasing" |
| D-FRC-02 | Add mode (08 Q17) | factors add, latest arrival time kept, warning | one (a, t0) pair per DOF |
| D-FRC-03 | Load-case files (02 OQ12) | `FILE9001…FILE9500` (3 digits) | examples |
| D-FRC-04 | FREAD / MREAD | not implemented (undocumented syntax): "not available" message | manual names only |
| D-ANL-01 | Eq. 2.1 generalisation (02 OQ2) | general assembly (§4.6 step 4) | handles structure on w nodes |
| D-ANL-02 | F_ff details (R1 §10.5) | 3×3 rotation blocks; r = 0 axis values; 0 < r < R0 linear in r; symmetrise; 2D sign rule | R1 §4.1 |
| D-ANL-03 | Solver | Schur complement: sparse LU (scipy `splu`) of C_nn, dense LU (`lu_factor`) of S + X_ff; never Cholesky | complex symmetric, non-Hermitian |
| D-ANL-04 | Dense impedance | `X_ff` formed explicitly by LU solve of F_ff with identity (needed for COOX and global impedance) | manual stores X_ff |
| D-ANL-05 | Restart files (02 OQ6, 03 OQ6, 05c §A.2) | COOXqqq = X_ff; COOTKqqq = (L, U, perm_r, perm_c of C_nn) + (lu, piv of S+X_ff); records keyed by frequency value with FILE90 hashes; COOXI/COOTKI index; subsets allowed | true reuse; verifiable equivalence |
| D-ANL-06 | Simultaneous cases (05c OQ11, 07 OQ10) | `[simul]` 0 = single FILE1/FILE9 → FILE8; seismic coherent 1 = FILE1X/Y/Z → FILE8X/Y/Z; seismic incoherent Ns = FILE77001… (+FILE1X/Y/Z) → FILE8001…; vibration Nl = FILE9001… → FILE8001…; simul = 1 with angle ≠ 0 → error | manual semantics |
| D-ANL-07 | Global impedance outputs (05c OQ14) | translational interaction DOFs only; files per frequency row: `f` + 6 diagonal (option 1) or 36 row-major values (option 2); FOUNDAMP = Im/(2\|Re\|); FILE11 holds X_ff and T per frequency | definitions explicit |
| D-ANL-08 | Frequency survey | match integer frequency numbers (\|f − nΔf\| < 1e-6Δf) | exact |
| D-ANL-09 | Mode 6 (05c OQ22) | out of scope; error message | vendor-only |
| D-ANL-10 | Memory guard (02 §3.7) | estimate `3·(3N_int)²·16 B` + sparse factor estimate; refuse above 80 % of physical RAM (`EDUOPT,MEMLIMIT`) | clear failure instead of "access violation" |
| D-ANL-11 | Angle range (05c OQ15) | [0, 360); negative normalised with warning | 0 is the default |
| D-ANL-12 | SYMM (07 OQ22) | P1: structure BCs on symmetry planes (symmetric: normal translation and in-plane rotations fixed; antisymmetric: in-plane translations and normal rotation fixed); soil flexibility by image superposition (formula corrected by D-W3-13); global impedance and incoherency forbidden with SYMM; VP compares half and full models | manual rules; standard image method |
| D-ANL-13 | Output of fixed DOFs | H = 0; rotations included for all 6-DOF nodes | MOTION rotational output |
| D-CMB-01 | Duplicate frequencies (03 OQ4) | error; `EDUOPT,COMBINDUP,PREFER82` uses FILE82 | explicit |
| D-CMB-02 | COMBIN deck (04 OQ18) | no deck; RUNCOMBIN reads FILE81/FILE82; the AOPT COMBIN flag only enables CHECK of their presence (warnings) | manual lists none |

### 7.9 MOTION, RELDISP, STRESS (MOT, RDP, STR)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-MOT-01 | Interpolation algorithms 0–5 (01 OQ2, 05c OQ1, R1 §10.4) | Tajirian 2-DOF ω² rational form on 5-point windows (§4.8), lstsq for rank deficiency, Lagrange guard; windows per the §4.8 table; options 4/5 single window (no averaging) | sourced formula (MHI RAI, TLUSH); manual "no averaging" |
| D-MOT-02 | Option 6 (01 OQ3, 05c OQ4) | separate not-a-knot cubic splines on Re and Im, including the f = 0 anchor | "bicubic" = cubic on both parts |
| D-MOT-03 | Range ends (05c OQ5) | above f_N: 0; seismic below f_1: complex linear from H(0) (rigid-body projection) to H_1; vibration: H_1; Nyquist bin real | rigid-body limit (manual low-frequency ATF ≈ 1 check) |
| D-MOT-04 | Smoothing algorithm (05c OQ2) | band filter of 05c B.6.7; S = 0 off; ignored for option 6 with warning | manual purpose; no-op for coherent |
| D-MOT-05 | Phase adjustment (05c OQ3) | unwrapped phase scaled by ρ = 1/(1+S) (options 0–5), ρ = 0 (option 6); reference f = 0 | only reading consistent with manual guidance |
| D-MOT-06 | SRSS order (05c OQ17) | interpolate then combine (default); `EDUOPT,SRSSORDER,COMBINEFIRST` alternative | SRSS of amplitudes is not interpolable linearly |
| D-MOT-07 | RS definitions (05c OQ7) | log spacing, endpoints included, `<fstep>` = number of points; SA = absolute acceleration; SV = relative velocity; Nigam–Jennings exact; full record | 301 points over 0.1–100 Hz = 100/decade |
| D-MOT-08 | Baseline codes (05c OQ6, 07 OQ17) | `<bl>` 0 none, 1 Hudson–Housner (time domain); FILE13 only when on | dialog wins |
| D-MOT-09 | Duration (05c OQ18) | dur = 0 → full Fourier period; output 1.2·dur clipped to NFFT·Δt; RS uses the full record | manual 20 % rule |
| D-MOT-10 | Duplicate output nodes (07 OQ18) | AFWRITE merges per (node, direction) with OR of flags, warning EDU-16 (Error 66 reported by CHECK as warning in this case) | manual says both error and auto-delete |
| D-MOT-11 | Which FILE8 is read (05c B.2) | default `FILE8`; `EDUOPT,TFFILE,<name>` sets the TF file for MOTION and STRESS (e.g. FILE8X) | avoids copying in batch while keeping the default |
| D-MOT-12 | Records (05c OQ19, 07 OQ28) | a record = one acceleration sample after the header/time-step line; pair format: one record per line, Δt checked against the time column | explicit |
| D-MOT-13 | Vibration output control (05c OQ8) | MOTIONX `<resp>` 0 displacement, 1 velocity, 2 acceleration (default); RS from acceleration | footnotes reference it |
| D-MOT-14 | MOTION NFFT/Δt vs FILE8 (05c OQ21) | must reproduce FILE8 Δf (1e-6) else error | interpolation grid defined by Δf |
| D-MOT-15 | RS damping count | > 5 values: EDU warning, all computed | manual MOTION limit 5 |
| D-RDP-01 | RELDISP algorithm (05d §2.1) | `H_rel = H_node − H_ref` on complex TFI; seismic `D = H_rel·(−gA/ω²)`; vibration `D = H_rel·F_ref` | standard form |
| D-RDP-02 | `.TFD` content (05d OQ-R2) | complex relative-displacement TF per unit control acceleration (length per g), amplitude + phase (rad) | directly usable |
| D-RDP-03 | RELD codes (05d OQ-R3, 07 OQ24) | `<RelDisOutput>` bit 0 = complex TFD (default 1), bit 1 = use TFU instead of TFI (discouraged, warning); `<RelDispSAll>` 1 = all nodes + frames; `<RelDispNumFiles>` written = RDND count, ignored on input | round-trip storage of raw values |
| D-RDP-04 | Free-field reference (05d §2.1) | RELFILE keyword `FREEFIELD` uses a synthetic unit-amplitude, zero-phase reference in the input direction | manual workflow made scriptable |
| D-STR-01 | Where interpolation applies (05d §1.2) | interpolate the STF, not nodal TFs | manual |
| D-STR-02 | Convolution (05d OQ-S1, OQ-S2) | seismic `U_g = −gA/ω²`; vibration F_ref | D-CNV-06 |
| D-STR-03 | Output code digits (05d OQ-S4) | digit k = component slot k of spec 05d §1.6; 0/1/2 | screenshot and EOUT |
| D-STR-04 | SOLID stress vs strain (05d OQ-S5) | components 1–6 output stresses; strains additionally with `EDUOPT,STRAINOUT,1` (files suffixed `E`, e.g. `SOLID_001_00001_EXX`) | stresses are the design quantity |
| D-STR-05 | Sign conventions (05d OQ-S8) | beam end forces = forces exerted on the element at I and J in local axes; spring `F = k*(u_J − u_I)` per global component (manual warns sign ≠ tension) | standard stiffness convention |
| D-STR-06 | `<itran>` scope (05d OQ-S13) | all element types | uniform behaviour |
| D-STR-07 | Skip Time History Steps (05d OQ-S3) | STRESSX `<skip>` = output every k-th sample (0/1 = all), same meaning as MOTION `<step>` print step | Error 67 shared |
| D-STR-08 | Frames.txt numbers (05d OQ-S11) | 1-based time-step indices | frame naming `_00001` = t 0 |
| D-STR-09 | Soil pressure (05d OQ-S10) | element pressure = −n·σ·n on the SOLID face shared with the structure (compression positive); nodal = plain average | manual "normal stress on the face" |
| D-STR-10 | Units (05d OQ-S7) | SHELL membrane outputs are stresses (F/L²), moments F·L/L; TSHELL N, Q per unit length | manual |
| D-STR-11 | Mixed-group beams (05d §1.5 warning) | identical results whether alone or mixed (regression test) | correctness |
| D-STR-12 | `.ess` layout (10 Q-4) | ELEMENT_CENTER layout with signed values at one step | consistent readers |

### 7.10 SOIL and EQUAKE (SOL, EQK)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-SOL-01 | SHAKE conventions | SHAKE91 recursion, outcrop 2E, mid-layer strain, SASSI complex modulus | VP-04 |
| D-SOL-02 | Half-space in SOIL (05a OQ14) | last SPRO entry; CHECK warning if its L differs from SITE `<hs>` | SHAKE convention |
| D-SOL-03 | Curve units and interpolation (05a OQ15, 07 OQ5) | DYNP strain and damping in percent; linear in log10 γ; constant extrapolation | SHAKE91 code |
| D-SOL-04 | Convergence (05a §6.7) | exactly `<iter>` iterations; report ΔG and Δβ; optional `EDUOPT,SOILTOL,<pct>` early stop | SHAKE91 behaviour |
| D-SOL-05 | Curve limits (05a OQ16) | 100 curves, > 11 points warned; Error 96 (> 15 used) only with LIMITS,PREP | V3 limits |
| D-SOL-06 | Vp/βp after iteration (R1 §6.3) | default constant Poisson's ratio (Vp scales with Vs) and βp = βs; `EDUOPT,VPPOLICY,VP` keeps Vp (saturated soils) | sources silent; both options needed |
| D-SOL-07 | `<gravmult>` (05a OQ17) | multiplies RS ordinates (SA in g × gravmult) | unit scaling |
| D-SOL-08 | Cut-off | Nyquist (manual); `EDUOPT,SOILCUTOFF,<Hz>` for SHAKE91 comparisons | VP-04 needs 25 Hz |
| D-SOL-09 | SOIL units (07 OQ31) | British kip–ft–ksf–kcf, SI kN–m–kN/m²–kN/m³; lb/ft³ input is not converted (user responsibility) | one model unit system |
| D-SOL-10 | Legacy SOIL argument order (07 OQ21) | detected by the parenthesised token; warning; not auto-mapped | order unknown |
| D-SOL-11 | SOIL scaling fields (05a OQ6) | SOIL has its own mult/max (SOILX), separate from MOTION's | screenshots show different values |
| D-SOL-12 | FILE88 → SITE matching | by position: SOIL sublayer k ↔ TOPL entry k (excluding the half-space); count mismatch → EDU-07 error | unambiguous |
| D-SOL-13 | Shared control layer and history file (05a §6.2, 07 §3) | shared with SITE/MOTION as in the manual; SOILX optional `<cl>` and `<file>` override them for SOIL only | allows rock-outcrop SOIL input with a surface SSI control point in one model |
| D-EQK-01 | `<accopt>` encoding (05a OQ7, 07 OQ4) | 0 none, 1 seed record, 2 external (RS/PSD/FFT only) | two dialog options |
| D-EQK-02 | Initial motion | SIMQKE-type PSD-based amplitudes, random or seed phases, Saragoni–Hart envelope | standard practice |
| D-EQK-03 | Baseline correction | zero-mean velocity in frequency + order-≤ 3 polynomial displacement fit | manual "complex-frequency + polynomial" |
| D-EQK-04 | Acceptance checks (C-21) | manual criteria and SRP 3.7.1 Rev 4 criteria, reported separately | both are cited |
| D-EQK-05 | Seeds (05a OQ8) | 0 treated as 1; best trial by smallest maximum deviation | dialog default 0 |
| D-EQK-06 | Correlation (05a OQ10) | P1; per-2-s-window mixing then re-matching as seed | manual warning |
| D-EQK-07 | Output units | .acc in g; .vel/.dis in in/s, in (British) or cm/s, cm (SI); PSD in²/s³ or cm²/s³ | manual PSD units |
| D-EQK-08 | History file format | one-column with Δt on the first line (= THFILE `<fopt>` 0 and `.ACC`) | interchangeable files |

### 7.11 CHECK, AFWRITE, module runs and files (CHK, AFW, RUN, FIL)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-CHK-01 | Module attribution (11 OQ33) | model-level errors under a `MODEL` header and gate HOUSE and every module needing FILE4; errors on shared variables are reported under every module using them | AFWRITE gating correctness |
| D-CHK-02 | Break count (05a OQ19) | per message type per module | manual "per module" |
| D-CHK-03 | Errors 52/63 (11 OQ26) | valid [0, 360) | 0 is the default |
| D-CHK-04 | Error 60 (11 OQ27) | only when incoherent with embedded interaction nodes | dialog value 0 |
| D-CHK-05 | Geometry thresholds (08 Q19, 11 OQ30) | W2 if an interior angle < 15° or > 165°; warp w = \|n̂·(x4−x1)\|/√A: W3 > 1e-3, E12 > 5e-2; configurable | manual gives none |
| D-CHK-06 | Error 124 (11 OQ28) | interaction node with all three translations fixed | interaction DOFs are translations |
| D-CHK-07 | Error 8 scope (11 OQ29) | BEAMS (I = J) and 2-node GENERAL only; never SPRING | FILLPOOL/SOLIDPILE zero-length springs |
| D-CHK-08 | Warning 4 and K nodes (11 OQ31) | no warning for K-only nodes; still fixed in decks | harmless |
| D-CHK-09 | Additional checks | EDU-01 interaction node not on an interface (error); EDU-02 size limit; EDU-03 rounded NFFT shorter than the records used (error for SOIL/MOTION/STRESS/RELDISP); EDU-04 damping ≥ 0.5 (error); EDU-05 non-positive Jacobian (error); EDU-06 unrestrained shell drilling; EDU-07 FILE88 layer mismatch (error); EDU-08 excavated element layer ≠ TOPL layer at its depth; EDU-09 BBC slope ≠ elastic stiffness; EDU-10 passing frequency < f_cut (G-05/G-06/G-08); EDU-11 ν > 0.47 (G-09); EDU-12 non-FV method without FV validation (info, G-21); EDU-13 AS/SRSS on a flexible foundation (warning); EDU-14 smoothing ≠ 0 with coherent input or option 6; EDU-15 RS points < 301; EDU-16 duplicate output nodes merged; EDU-17 near-coincident coherency projections (G-18); EDU-18 post-run \|ATF\| at f_1 deviates > 5 % from 1 (G-19); EDU-19 quiet zone shorter than `ln(100)/(2π f_min β_min)`; EDU-20 memory estimate (G-22); EDU-21 interaction nodes not bottom-up ascending; EDU-22 excavation-boundary node not an interaction node (FV/MSM/FFV) | specs 01 C1–C17, 03 G-01…G-29 consolidated |
| D-AFW-01 | AFWRITE with errors (01 OQ27) | the module's deck is not written; an existing stale deck is renamed `<model>.<ext>.bak` with a warning | prevents running outdated input |
| D-AFW-02 | Deck formats | D-GEN-04 | — |
| D-AFW-03 | "Extend integer fields" | accepted and stored; no effect (free-format decks) | fidelity of the option set |
| D-AFW-04 | "Simulation Commands" (05a OQ18) | RUN commands for AOPT-enabled modules in run order (and FCOPY lines for simultaneous-case files when configured), to `<model>-Sim.pre` or appended to `.pre` | most likely meaning; scriptable |
| D-AFW-05 | Default AOPT | SITE, POINT, HOUSE, ANALYS, MOTION enabled; others off | avoids spurious errors in new models |
| D-AFW-06 | Write-time fix-ups | ETYPE 0 resolved; coordinates to global; gap/unused nodes fixed; empty groups skipped; NFFT rounded; MOTION duplicates merged | manual §9.2.3 |
| D-RUN-01 | RUNSITE (10 Q-1) | implemented with the RUN semantics | manual §6.4.5 references it |
| D-RUN-02 | Module locations | built-in Python modules by default; an external executable is allowed only if it accepts SASSI-EDU decks | vendor executables cannot read these decks |
| D-RUN-03 | Prerequisite checks | §2.6 rules before launching | clear errors |
| D-RUN-04 | Batch protocol | three stdin lines (or `--model/--input/--output` arguments); exit code 0 success, 1 input error, 2 numerical failure | legacy batch files stay readable |
| D-RUN-05 | Synchronous runs (10 Q-2) | synchronous in `.pre`/macros; GUI uses a worker process with streaming and Cancel | sequencing |
| D-FIL-01 | Binary file names | legacy names, npz content, no extension | nomenclature |
| D-FIL-02 | Text result formats (05c OQ16, 06 OQ1, 10 Q-22) | `.TFU/.TFI/.TFD`: optional `#` header line, then `f amp [phase_rad]`; `.ACC/.THD/.THS`: first line Δt, then one value per line; `.RS/.RSO`: `f SA`; non-numeric lines are skipped by readers | READSPEC column counting works; ACC reusable as THFILE |
| D-FIL-03 | Frame files (06 OQ16, 09 Q-G6) | header `nrows ncols`; rows `node v1 … v(ncols−1)` (node id is column 1); complex values as Re/Im column pairs; list files: first line ignored, then one path per line (relative to the list file) | MODFRAMES edits the second number |
| D-FIL-04 | Element component codes (05d OQ-S6) | as §1.8 | complete naming |
| D-FIL-05 | Node digits | 5, or 6 when the maximum node number > 99,999 | manual |
| D-FIL-06 | ELEMENT_CENTER file (02 OQ17) | six columns: SOLID Sxx Syy Szz Sxy Sxz Syz; SHELL Sx′x′ Sy′y′ Sx′y′ Mx′x′ My′y′ Mx′y′; SPRING Fx Fy Fz Mxx Myy Mzz; PLANE Sxx Szz Txz 0 0 0; "ordered group #" = 1-based order of the group among groups of the same type | manual example layout |
| D-FIL-07 | SRSSTF.txt, CONTTRS.txt, Frames.txt | manual layouts; case-insensitive lookup; CONTTRS one path per line | manual |
| D-FIL-08 | FILE74 / FILE78 (05b OQ23) | text: header (iteration, ESF), rows `group element γ_eff G β curve` | readable convergence check |
| D-FIL-09 | `.map` files (05b OQ21) | `old new` per line | manual |
| D-FIL-10 | FILE73 / FILE88 type (01 OQ6) | text files (curves; strain-compatible layer table) | §1.5.4 prose; inspectable |

### 7.12 UI, plotting, line maths and section cuts (UI, LIN, SEC)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-UI-01 | Interpreter sharing | UI-01 | verifiability |
| D-UI-02 | Toolkit | D-GEN-05 | environment |
| D-UI-03 | Dialog defaults (05b OQ1, 05c OQ20) | §5.4 "Default" column; manual screenshot values are examples only. Exceptions to note: Save Complex TF on (RELDISP needs it), Print Amplitude Only on, Equiv. strain ratio 0.65, EQUAKE duration 20 s | recommended practice |
| D-UI-04 | Enable/disable logic (05b §2.5) | §5.4 rules | prevents invalid combinations |
| D-UI-05 | Load Model "Add Model" fields (04 OQ2) | group, model name, model directory (created if missing), title | minimal useful set |
| D-UI-06 | 3D rendering | matplotlib mplot3d poly collections (≤ ~50k faces), plotly HTML export for large models and animations; orthographic projection | no VTK available |
| D-UI-07 | Persistence (05d OQ-U2) | colours, fonts, shader options and Command Display persisted in SASSIini.xml; Check Options and toolbars not persisted | manual + convenience |
| D-UI-08 | Dialog commands in batch | error ("requires arguments in batch mode") | headless safety |
| D-UI-09 | ElemPalette index (05d OQ-U1) | colour = `((n−1) mod 128) + 1` | 1-based palette |
| D-UI-10 | Image format (04 OQ22) | CAPTUREPLOT: PNG if `.png` in the name, else BMP; File ▸ Export Image defaults to `.bmp`, PNG selectable | manual |
| D-UI-11 | Hide/show and clipping (06 OQ10–11) | hide sets stored in the model (saved by SAVE, not WRITE); element clipped only when all its nodes are outside the display volume; one list or one range per press | manual statements |
| D-UI-12 | View conventions (06 OQ9, 10 Q-31) | CNGVIEW angles in degrees, zoom = scale factor (1 = fit); keys rotate 5°; default view isometric (rX −60°, rZ −45°) | screenshots are oblique |
| D-UI-13 | Vector rule for Y and Z (06 OQ19) | Y: (b, a, b), Z: (b, b, a) × Scale | generalisation of the manual's X example |
| D-UI-14 | Bubble size (06 OQ15) | size = S_max·max(0.2, t), t = colour-map position | visible minimum |
| D-UI-15 | SASSIani.xml (06 OQ17) | `<Animation description directory type frames start end stride scale cmin cmax listfile/>` | complete for reload |
| D-UI-16 | Line names (06 OQ3) | READSPEC column 1 → file name; column k > 1 → `<file>[k]` | distinguishable legends |
| D-UI-17 | Export Table (06 OQ27) | spectrum, time history, soil properties; CSV with header (x name, line names) on the union grid | spreadsheet-ready |
| D-LIN-01 | BROADEN (06 OQ5, 10 Q-19, Q-21) | envelope → broaden ±b on the augmented grid → bridge valleys (valley-depth criterion); `EDUOPT,BROADENGRID,ACS` evaluates on source points only; `EDUOPT,BRIDGE,AMPLITUDE` alternative criterion | ASCE 4 practice; T-B tests pin behaviour |
| D-LIN-02 | Undefined lines (10 Q-34) | error in line maths; ignored in SPECPLOT/THPLOT | silent omission unsafe |
| D-LIN-03 | LBINCORS (10 Q-20) | parsed; "algorithm not specified in the manual" error | no invented rule |
| D-LIN-04 | READTH Pair 0 (10 Q-22) | first number = Δt, then values (several per line allowed); t_k = (k−1)Δt | matches D-FIL-02 |
| D-SEC-01 | Local axes (10 Q-6, Q-30) | ez = n, ex = projected r, ey = ez × ex; origin at area centroid; sysno 0/blank not stored | manual "r = local +X" |
| D-SEC-02 | Output order (10 Q-7, Q-8) | CALCPAR: Area, Xc, Yc, Zc, Ixx, Iyy, Ixy, Izz, Fx, Fy, Fz, Mx, My, Mz; CSV `Time,Fx,Fy,Fz,Mx,My,Mz` + `MAX` row for both CALCSECTHIST and CALCSECTHISTDB | explicit |
| D-SEC-03 | Shell bending in resultants (10 Q-12) | included by default (`EDUOPT,SECTBEND,0` excludes) | physically complete |
| D-SEC-04 | Coincident faces (10 Q-13) | counted once (element on the −n side preferred) | no double counting |
| D-SEC-05 | Time column (10 Q-9) | (k−1)·ts; step k when ts ≤ 0 | frame 1 = time 0 |
| D-SEC-06 | READSTR on maxima file (10 §3.9) | warning that maxima are not simultaneous | physical consistency |

### 7.13 ANSYS interfaces (ANS)

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-ANS-01 | `.cdb` records (04 OQ14) | NBLOCK, EBLOCK, ET, KEYOPT, MPDATA (EX; PRXY or NUXY; DENS → γ = DENS·g; DMPR → β = DMPR/2 for both βp and βs; others ignored with warning), RLBLOCK/RBLOCK, SECTYPE/SECBLOCK/SECDATA, D with value 0 (converted to fixities), nodal rotation angles / CE / CP / ENDRELEASE (warning); the gravity argument also sets the model GRAVITY | supported subset per manual; damping mapping explicit |
| D-ANS-02 | Grouping (04 §11.6) | one group per ANSYS element type number; elements renumbered 1…n; node numbers kept; map file `<model>_cdb.map` (`ansys_elem group elem`) | manual group rules |
| D-ANS-03 | Beam axes and releases (04 OQ15) | computed per element by matching ANSYS and SASSI local frames from I, J, K (no hard-coded table); releases mapped through the same frame; unit test with B ≠ H section against ANSYS | robust to convention differences |
| D-ANS-04 | Oblique COMBIN14 (04 §11.3) | axial spring along a non-axis I–J → 2-node GENERAL element `k·e eᵀ` | SC is global and uncoupled |
| D-ANS-05 | Masses | MASS21 → MT/MR with MUNITS 0 (mass) | ANSYS masses are masses |
| D-ANS-06 | Export (04 OQ17, 09 Q-F1) | legacy V11–15 elements by default (SOLID45, SHELL63, BEAM44 releases by KEYOPT, COMBIN14 per non-zero component, MASS21, MATRIX27, PLANE42 plane strain); `EDUOPT,ANSYSMODERN,1` writes SOLID185/SHELL181/BEAM188; ETYPE 0 above grade counts as structural; excavated soil and layers not exported; interaction nodes as component `SSI_INT`; damping as `MP,DMPR,2β` | manual + modern option |
| D-ANS-07 | SOLID185 condition (04 OQ14) | accept KEYOPT(3) = 0 (non-layered); warn otherwise | manual's "nonuniform materials" refers to the layered option |
| D-ANS-08 | CONVERT,SSI scope | P1 reads SASSI-EDU `.hou`/`.sit`/`.poi` (round trip); legacy SASSI2000 fixed format OOS until its record layout is available | format unknown |

### 7.14 Nonlinear soil SSI (NLS, P1) and Option NON (NON, P2) and other P2 items

| ID | Question (source) | Decision | Rationale |
|---|---|---|---|
| D-NLS-01 | Iteration driver | no new command; shipped macro `nlsoil_iter.pre` loops RUNHOUSE → RUNANALYS (Mode 2, ×3 via TFFILE) → RUNSTRESS ×3 → COMBXYZSTRAIN | manual batch practice |
| D-NLS-02 | Convergence (01 OQ16, 03 OQ19) | max \|ΔG/G\| < 2 % and max \|Δβ\| < 0.5 % absolute; default limit 8 iterations; report per iteration | SHAKE-style practice |
| D-NLS-03 | Strain measures | per ISTR as §4.10 item 5; ESF from `.pin` | manual |
| D-NLS-04 | COMB_XYZ_STRAIN | `COMBXYZSTRAIN` SRSS of γ_eff per element; output FILE74 | manual formula |
| D-NON-01 | `.eql` format (05d OQ-N1) | keyword text deck (D-GEN-04) carrying EQL, BBC, P, S records | none documented |
| D-NON-02 | `<NonLinOpts>` (05d OQ-N10) | bit mask 1 panels, 2 springs, 4 beams | three check boxes |
| D-NON-03 | Hysteresis rules (05d OQ-N3) | CMS per Cheng & Mertz (1989), TAK per Takeda et al. (1970) with unloading exponent 0.4, GMR = Masing factor 2 on the BBC | published models |
| D-NON-04 | Damping order (05d OQ-N6) | `ξ = min(cutoff, scale·ξ_h + [ξ_el])`; scale 0 → 1; cutoff 0 → none | screenshot zeros |
| D-NON-05 | Unsupported options (05d OQ-N12) | keep Errors 126–128; `EDUOPT,NONEXT,1` unlocks TAK and bending panels as experimental | fidelity first |
| D-NON-06 | Convergence (05d OQ-N5) | max relative change of E < 2 %, of ξ < 0.5 % absolute, ≤ 10 iterations | explicit |
| D-NON-07 | File renaming (05d OQ-N14) | performed by the NONLINBAT script, with the corrected names (`_elastic.THS`) | Fig. 1.2 typo |
| D-NON-08 | `COMB_XYZ_THD.inp` (01 OQ21) | line 1 `n`; then n lines `outfile fileX fileY fileZ`; algebraic sum per time step | simplest complete form |
| D-NON-09 | BBCGEN shape (11 OQ10) | V(γ) = V_cr + (V_u − V_cr)[1 − (1 − ξ)²], γ_y default 0.004, γ_cr = V_cr/(G·A_w), type 1 (CMS), written to the panel's BBC number (panel number if 0) | smooth, zero slope at yield |
| D-NON-10 | SHEAR/BBCGEN args 6–7 (10 Q-14, 11 OQ8) | A_BE and f_y,BE (F_BE = A_BE·f_y,BE; F_VW = ρ_V A_W f_y); `EDUOPT,SHEARFORCEARGS,1` reads them as forces | descriptive text |
| D-NON-11 | K_el (05d OQ-N7) | BBC first slope Y1/X1; EDU-09 warns if it differs > 1 % from G·A_shear, E·I or k | manual normalisation rule |
| D-NON-12 | Panel materials (05d OQ-N13) | each panel must have its own material (GROUPMAT/UNIPNL); shared → error | per-panel E updates |
| D-BDB-01 | Binary DB layout and combination (11 OQ20–22) | npz schema of spec 11 §3.6; COMB*DB algebraic time-wise sum; one DB per type in memory | spec 11 proposal |
| D-WAT-01 | FILLPOOL (11 OQ1–5) | water ν = 0.49, E from K = 2.2 GPa, 62.4 pcf / 9.81 kN/m³, 0.5 % damping; oblique-wall springs as coupled GENERAL elements; spring damping 0; floor nodes get Z springs; shell-area group dummy (very low E, zero mass) | documented educational approximation |
| D-TSH-01 | THSHLSTR faces (11 OQ23) | σ = N/t ± 6M/t² with the four sign permutations of (N, M) maxima; τ_xy with `++` only; principal values from those | manual description |

### 7.15 Items that remain open (tracked, not blocking P0)

1. Abrahamson coherency coefficients: models 3, 5 and 6 are implemented with cited coefficients (wave 3);
   models 2 (Abrahamson 1993) and 4 (2006 embedded) remain unavailable until their coefficients are
   confirmed from a primary source.
2. Tajirian & Tabatabaie (1985) compliance curves need digitising for a quantitative VP-13.
3. ~~The half-space sublayer law default is confirmed only after VP-02b (D-SIT-02).~~ Resolved: UNIFORM
   (D-W1-01).
4. Kausel layer vertical/torsion formulas and the Pekeris sign (VP-13, VP-43) need primary-source
   confirmation.
5. ~~CMS/TAK hysteresis details for Option NON must be taken from the original papers (P2).~~ Done (wave 3):
   CMS follows HYST04 of INRESB-3D-SUP (Cheng & Mertz 1989), TAK the Takeda/Otani rules (HYST06); see
   `docs/user/OPTION_NON.md`. ACS SASSI's own implementation is unpublished, so loop details may differ.
6. LBINCORS and Mode 6 remain unimplemented until the algorithms are documented.
7. Legacy SASSI2000/PREP fixed-format deck import requires the original record layouts.

### 7.16 Lead decisions after implementation wave 1 (binding; supersede the cited rows)

| ID | Subject | Decision | Evidence |
|---|---|---|---|
| D-W1-01 | D-SIT-02 half-space law | VP-02b failed with GEOMETRIC (1.10 % at nl = 20) and passed with UNIFORM (0.20 % / 0.95 % at nl = 20 / 10); per the D-SIT-02 fallback the default is **UNIFORM** (`EDUOPT,HSLAW` and the SITE deck default) | VP-02b |
| D-W1-02 | VP-05 criterion | Linear thin layers lock volumetrically for ν ≥ 0.45 (1.24 % at ν = 0.45, 2.95 % at ν = 0.49 with 10 sublayers per λ). Criterion: 0.5 % at 10 and 20 sublayers per λ with order ≈ 2 for ν ≤ 1/3; 0.5 % at 40 sublayers per λ with monotone convergence for ν ≥ 0.45 (users are warned by G-09/EDU-11 for ν > 0.47) | VP-05 |
| D-W1-03 | VP-06 criteria | (a) The variable-depth half-space cannot represent the grazing field near a Love cut-off; discrete cut-offs are 0.85 / 12.40 / 23.92 Hz vs 0 / 11.547 / 23.094 Hz. Cut-offs are checked with a **provisional 10 %** tolerance (absolute 10 % of 11.547 Hz for the zero cut-off). (b) "Convergence from above" holds for consistent mass only; with the mixed mass of D-CNV-05 the requirement is monotone convergence of the 20 Hz phase velocities. Phase velocities keep 0.5 % | VP-06 |
| D-W1-04 | VP-09 tolerance at 3.3 R0 | 2.35 % (R1 V8's prototype value 2.310 % had been rounded to 2.3 %; SASSI-EDU gives 2.308 %) | VP-09 |
| D-W1-05 | VP-38 mesh | The lumped translational SHELL mass (D-ELM-07) converges from below at O(h²) and is 2.4–6.7 % low on the NAFEMS 8×8 mesh. FV12/FV16 are checked on a 32×32 mesh (2 %) plus an observed order ≈ 2 (8/16/32 meshes); FV32 on the NAFEMS 16×8 mesh. The FV16 fundamental converges to λ₁ = 3.4710, the accurate superposition and Legendre–Ritz value; NAFEMS 0.421 Hz is +0.73 % high. FV12: NAFEMS mode 4 (4.233 Hz) is +1.02 % above the exact (Leissa) value and mode 5 +0.82 %, so the 32×32 SHELL mode 4 is −0.37 % from exact but −1.38 % from NAFEMS | VP-38 |
| D-W1-06 | VP-04 spectral accelerations | The printed SHAKE91 spectra are spectra of a record corrupted by SHAKE91 subroutine DRCTSP (samples that set a new running maximum are not converted from g). Emulating that defect on SOIL's surface motion reproduces all 151 printed values within 0.6 %; this is the VP-04 spectrum criterion. The correct 5 % spectrum of the motion is reported in the VP notes (≈ 0.249 / 0.833 / 0.190 g at T = 0.1 / 0.4 / 1.0 s) | VP-04 |
| D-W1-07 | Response-spectrum sampling (D-MOT-07 amended) | In addition to the 4× FFT up-sampling above 0.1/Δt, the piecewise-linear input is sub-divided so every oscillator period has ≥ 64 integration steps (`sassi.core.spectra.STEPS_PER_PERIOD`), resolving response peaks between samples (previously up to −1.2 % at 20 samples per period). `substep=False` reproduces the plain sampling | spectra unit tests, EQUAKE bank tests |
| D-W1-08 | D-SIT-09 mode selection | "Propagating" = \|arg k\| ≤ atan(0.5) + atan(loss) with loss = tan(δ/2) of the column (damping rotates all wavenumbers); least-decay ties within 1e-8·max\|k\| are broken by the largest Re k. FILE3 carries `x_loss` | SITE/POINT review |
| D-W1-09 | D-EQK-03 baseline | Degree-5 displacement polynomial (two free coefficients fitted by constrained least squares; cubic acceleration correction) after the zero-mean-velocity correction | EQUAKE review |
| D-W1-10 | D-MOT-01 guard | The literal guard (min\|D\| < 1e-3·max\|D\|) destroys exact fits on lightly damped poles; the guard triggers only on spurious poles (no cancelling zero, damping < 1e-4 or remote), see `sassi.core.interp` | VP-28 |
| D-W1-11 | D-MOT-05 phase adjustment | `H' = \|H\|·exp(i(φ0 + ρ(φ − φ0)))` with φ0 the f = 0 (or first non-zero) phase, so negative rigid-body anchors are preserved; ρ = 0 when options 0–5 fall back to spline/linear | VP-33 |
| D-W1-12 | COMBIN | Also rejects different `case`, `type`, `cm` and `ang` meta values | VP-25 |
| D-W1-13 | DOFSMAP / FILE90 / FILE91 / COOX | Defined in `sassi.io.files.SCHEMAS`; HOUSE writes DOFSMAP, FILE90 (hashes) and FILE91 | contract |
| D-W1-14 | SHELL recovery units | Membrane outputs FXX FYY FXY are stresses (D-STR-10); 'force' per unit length is available as an element option | elements review |
| D-W1-15 | SOIL deck | `sstr` columns are opt1 compute stress, opt2 save SSxxx.TH, opt3 compute strain, opt4 save SNxxx.TH; `ssaf.freqstep` is a real (Hz); `vppolicy` carries `EDUOPT,VPPOLICY` | SOIL review |

### 7.17 Lead decisions after implementation wave 2 (binding)

Approximate published references that were superseded by rigorous ones remain in the VPs as
**informative** comparisons (`VPResult.inform`, kind `info`): they are printed with their
deviation but are not pass/fail criteria. The superseding references are independent
implementations validated against exact solutions (see `docs/verification/impedance_study.md`):
a static BEM with Boussinesq/Cerruti kernels (reproduces Mossakovskii's welded vertical stiffness
6.1299 GR, Reissner–Sagoci torsion and the relaxed 4GR/(1−ν), 8GR³/(3(1−ν)) to 0.1 %) and a dynamic
BEM with exact Lamb kernels (reproduces the Wong/Apsel table of R2 D.3 to 0.001).

| ID | Subject | Decision | Evidence |
|---|---|---|---|
| D-W2-01 | VP-11 rocking damping | The Veletsos–Verbic c_r (R2 C.2) is 28–53 % below the rigorous relaxed solution and the exact low-frequency limit c_r = 0.240 a0² (ν = 1/3); criterion: relaxed c_r from F_ff vs the exact-kernel relaxed BEM within 10 % (observed ≤ 1.2 %). VV c_r is informative; all other VV coefficients and the high-frequency dashpots stay criteria | VP-11 |
| D-W2-02 | VP-13 vertical resonances | The vertical compliance peaks of a disk on a layer are not at (2n−1)Vp/(4H): the first sits below the zero-group-velocity onset of the second P-SV mode (A0 = 0.509, continuum), the second is masked by 5 % damping. Criteria: SITE ZGV onset vs continuum dispersion (1 %), β = 5 % peak below and within 6 % of the onset, β = 1 % second peak vs 3Vp/(4H) (2 %). The 1-D rule is informative; horizontal resonances and statics remain criteria | VP-13 |
| D-W2-03 | VP-14 square rocking | The Pais–Kausel rocking fit (6.0 GB³) is 7 % below the welded BEM (6.457) and 4 % below the relaxed BEM (6.238); published fits disagree by 10 %. Criterion: K_xx = K_yy vs the welded BEM within 5 % (observed 0.7 %); K_z, K_x, K_y, K_zz keep Pais–Kausel | VP-14 |
| D-W2-04 | VP-17 peak | R2 F.1's 6.71 uses Veletsos–Verbic impedances (low rocking damping) and omits the sliding–rocking coupling of a welded mat. Criterion: peak and resonance vs F.2 fed with the rigorous welded impedance (5 %; observed 4e-5 and 0.2 %); the 6.71 value is informative. The F.2 identity (1e-8), statics, T̃/T and the 1.717 Hz resonance remain criteria | VP-17 |
| D-W2-05 | VP-51 EXCAV | "48 hexes and 75 nodes" is internally inconsistent (a 4×4 template on L node levels gives 25L nodes and 16(L−1) hexes); the 3-level basement gives 75 nodes and **32** hexes | VP-51 |
| D-W2-06 | D-FRC-02 wording | In add mode the factors are summed and the arrival time is that of the most recent definition with a non-zero factor (a zero-factor row never changes it), as in the interpreter (spec 08 §9.5) | FORCE review |
| D-W2-07 | FILE4 meta | `int_dofs` ([1, 3] in 2D, [1, 2, 3] in 3D; `int_eq` = −1 marks non-interaction DOFs) and `excav_layering` (R5/R6 findings); STRESS stops with an error when soil-pressure or nodal-contour output is requested and R5/R6 findings exist (D-HOU-04) | HOUSE review |
| D-W2-08 | Section cuts | D-SEC-07: CALCPAR/CALCMOI use stored CSECT pieces, else a recognised unit-thickness extrusion, else the mid-plane; D-SEC-08: BEAMS are drawn by CSECT but are never section pieces (information/warning lists) | cuts review |
| D-W2-09 | ACTIVATEPLOT | New extension command `ACTIVATEPLOT,<id>` (P1, full name) selects the active plot (GUI tabs emit it, L17) | post-processing |
| D-W2-10 | Memory estimate (D-ANL-10) | inputs + 3 dense copies (3.25 with `<save>` = 1, 2 in Mode 3) + bounded work + 30·nnz·16 B sparse fill + results/FILE11 buffers | ANALYS review |
| D-W2-11 | Temporary files | `run_problem` removes its temporary directory unless `keep=True`; pytest keeps only failed `tmp_path` directories | disk use |

### 7.18 Lead decisions during implementation wave 3 (binding)

| ID | Subject | Decision | Evidence |
|---|---|---|---|
| D-W3-01 | VP-38 TSHELL (NAFEMS Test 21) | As D-W1-05: the lumped TSHELL mass converges from below at O(h²) (mode 4: −5.4 % on 8×8, −1.4 % on 16×16, −0.35 % on 32×32 vs the element's own κ = 5/6 solution). The NAFEMS 8×8 rows are informative; the criteria are 2 % on 16×16 and 32×32 and an observed order ≈ 2 | VP-38T |
| D-W3-02 | THSHLSTR transport | The STRESS deck carries `thshlstr` (AFWRITE writes 0/1 from the model record, authoritative); −1 = not set (decks built without AFWRITE fall back to `THSHLSTR.opt`) | TSHELL review |
| D-W3-03 | HOUSEX messages | The node optimizer and non-linear soil SSI are implemented; HOUSEX no longer reports them as unavailable | nonlinear-soil / TSHELL reviews |
| D-W3-05 | VP-26 FFL vs FFM | FFL (`s ⊙ X U′`) and FFM (`X (s ⊙ U′)`) differ by the commutator [S, X], whose rigid-body projection Tᵀ[S, X]T is antisymmetric: they are identical only for the generalised force in the input direction, K_G-decoupled DOFs and coherent factors. Criteria: those identities plus the exact relation K_G(q_FFL − q_FFM) = Tᵀ[S, X]T t_d computed independently from FILE11 and FILE77; the horizontal ATF differences (≈ 1e-3 to 3e-2) are informative | VP-26 |
| D-W3-06 | D-INC-03 extension | With repeated eigenvalues of the coherency matrix (e.g. symmetric meshes) the eigenvectors are put in a canonical basis before the sign convention, so AS factors are reproducible and independent of solver details | incoherency |
| D-W3-07 | D-WAT-01 water model (amended) | FILLPOOL water is SOLID with bulk modulus K and G = 1e-8 K (not ν = 0.49, which makes the water move rigidly with the tank below ~5 Hz and gives the total instead of the impulsive mass); incompatible modes are required (FILLPOOL/MERGEPOOL set MOPT,0). Oblique wall interfaces use the diagonal projection of the coupled spring matrix (zero-length GENERAL elements are rejected by Error 8, D-CHK-07). New extension commands MERGEPOOL (P2 action) and POOLDATA (P2 record). Verified by VP-54W and VP-W1 (impulsive base shear +1.5 % to +1.9 % from exact potential flow at 2–8 Hz on the 0.5 m mesh, criterion 3 %; mesh convergence in `tests/unit/test_water_physics.py`) | water package |
| D-W3-08 | VP-46 SI references | The SI capacities of spec 10 §3.11 are printed to 0.1 kN; they are compared within half their last digit (0.05 kN). The SI computation path is checked at 1e-6 against the British results × 4.44822162 kN/kip | Option NON |
| D-W3-09 | Option NON integration | NONLINEAR is a module (`sassi.modules.base.MODULE_NAMES`, deck `.eql` written by AFWRITE from `build_eql`); CHECK applies the numbered Errors 121–128 (121 when neither panels nor springs exist; 126/128 lifted by `EDUOPT,NONEXT,1` for force 3/disp 2) and reports the further `build_eql` findings as EDU-44 (error) / EDU-45 (warning) | Option NON |
| D-W3-10 | SOIL-NON | Implemented per the decisions SN-1 … SN-12 recorded in the `sassi/core/soilnon.py` docstring (backbone fit with β = 1 and s ∈ [0.2, 1]; extended Masing rules; damping types 1–3; Newmark average acceleration with Newton iterations and step halving; rigid base or Joyner–Chen elastic base). NLSOIL/NLSLAYER travel in `<model>.nls` written by AFWRITE. Error 125 requires one NLSLAYER set per soil sublayer, not for the half-space. Verified by VP-SN1 (linear limit vs SHAKE ≤ 0.4 %), VP-SN2 (Masing loop damping 2e-7), VP-SN3 (informative) | SOIL-NON |
| D-W3-11 | Option A (LOADGEN) | Implemented as module LOADGEN (deck `<model>.lgn`, commands LOADGEN, LOADGENDYN, LGFILE, LGNODE, LGTIME, LGMAP, LGOPT, LGLIST, RUNLOADGEN,[STATIC\|DYNAMIC]) per the interpretations D-LGN-01 … D-LGN-08 of `docs/user/OPTION_A.md` §10 (the vendor's Integration manual is unavailable). Static: inertia forces −M a(t*) at critical times plus interface displacements; dynamic: APDL TABLE arrays, ACEL of the reference frame and relative interface displacements (initial value removed). Verified by VP-LA1 (base shear and moment vs STRESS ≤ 4e-4) and VP-LA2 (tables vs MOTION/RELDISP to 1e-11; replay of the exported second step with the HOUSE matrices reproduces the SSI accelerations to 5e-9). The APDL is checked by a parser; it has not been run in ANSYS | Option A |
| D-W3-12 | VP-T1 tolerance | 0.33 % at \|x\| = 2 R0 (the R1 V7 prototype prints 0.32 % for u_z due to P_z there; R1 rounded it to 0.3 %); 0.3 % at larger distances | 2D package |
| D-W3-13 | D-ANL-12 (corrected) and 2D scope | The image-superposition formula is `F_red(i,j) = Σ_S s_S F(i, j_S) P_S` over the 2^k plane combinations (P_S reflection, s_S = ±1 per symmetric/antisymmetric plane; equivalently `s P F(i′, j)`); the earlier text `F(i,j) ± P F(i,j′) P` was wrong. Symmetry-constrained translations are removed before inversion and ANALYS checks that the free field has the declared symmetry. Wave passage, multiple excitation, incoherency and global impedance are refused with SYMM and in 2D (manual §2.7); 2D is plane strain in X–Z (UX, UZ); anti-plane SH analysis is not implemented. Verified by VP-T2 (half/quarter = full to 1e-11), VP-T3 (2D zero-SSI identity 1e-15), VP-42 (strip statics within 5 % after extrapolation) | 2D package |
| D-W3-14 | EXCSTRCHK severity (refines §4.4 item 3 and G-15) | One rule shared by CHECK and HOUSE (`sassi.core.house_lib.excstrchk_kind`): an excavation-interior node shared with the structure is an **error**, unless the node is an interaction node and every structural element at it is a SOLID/PLANE on the excavation mesh (near-field soil, solid basement, the VP-16 same-mesh identity), which is a **warning**. Evidence (final audit F-01): shared vs separate nodes changes the response by 0.5–2.5 % with FV on solid blocks, 5 % to 478 % with FI-EVBN, and removes the resonance of a flexible floor slab on interior FV nodes (peak 1.0 instead of 9.7); manual §9.8.1, §1.5.1 rule 11 and guideline 13a | final audit |
| D-W3-04 | VERIFY work files | VERIFY removes its temporary work directory after a clean run and keeps it (and reports its path) only when a problem fails; informative comparisons are printed as `info` | disk use |


### 7.19 Lead decisions: built-in inputs and defaults of blank inputs (wave 5, binding)

ACS SASSI requires the input files of an analysis: a blank THFILE is Error 73, no RSIN file Error 84, blank
RSOUT / ACCOUT Errors 86 / 87, a SOIL profile without curves Errors 95 / 97 / 99. The user asked that the
standard inputs have defaults instead of having to be loaded. SASSI-EDU therefore ships a built-in input
library and fills these blanks from it, **reports every default it uses** so that a design engineer always sees
what input was assumed, and keeps the ACS behaviour one switch away (`EDUOPT,DEFAULTS,OFF`). The policy is one
module, `sassi/prep/defaults.py`, used by CHECK, AFWRITE, the RUN commands, `LIBRARY,DEFAULTS` and the GUI.

| ID | Subject | Decision | Evidence |
|---|---|---|---|
| D-W5-01 | Built-in input library | `sassi/data/library/` (package data, also in the browser bundle) holds the RG 1.60 0.30 g record `rg160h_030g.acc` and spectrum `rg160h_030g.rsi`, the RG 1.60 H and V 1 g spectra, the SRP 3.7.1 App. A target PSDs (1 g, and 0.30 g = 0.09 × 1 g; cm²/s³ and in²/s³), the 5 Hz Ricker pulse `ricker_5hz.th` and the SHAKE91 curves `dynp_library.pre` (README.txt: source and units of each). A file is named `@<file>` (case-insensitive), resolved from the package's own location in a source tree, a pip install and Pyodide (`/home/pyodide/sassi-edu/sassi`). The library is read-only: an `@` output name is refused (command error; EQUAKE Errors 86 / 87). New command `LIBRARY,[kind\|@name\|DEFAULTS]` (P0) lists it | user request; `tests/unit/test_input_defaults.py` (copies identical to their originals) |
| D-W5-02 | `@` resolution | Every reader of a file name resolves `@name`: the interpreter (`resolve_path`: INP, READSPEC, READTH ...), CHECK (`resolve_file`), AFWRITE (`deck_file_name` keeps the `@` name), the modules (`sassi.io.library.module_path`: EQUAKE, SOIL, MOTION, STRESS, RELDISP, LOADGEN) and the GUI file API (read-only). Variable substitution (L13) leaves `@<library file>` as it is: library names contain a '.', which variable names cannot, and the library name wins over a variable of the same stem | interpreter tests |
| D-W5-03 | THFILE blank, seismic | MOTION, STRESS and RELDISP of a seismic analysis (ANALYS `<type>` ≠ 1) read `@rg160h_030g.acc` (RG 1.60 H record matched by EQUAKE, PGA 0.324 g, 20 s, Δt 0.005 s) | ACS: Error 73 |
| D-W5-04 | THFILE blank, forced vibration | With ANALYS `<type>` = 1 (the `type` of the MOTION / STRESS / RELDISP decks) the load history is `@ricker_5hz.th` (5 Hz Ricker wavelet, peak 1, 2 s) | ACS: Error 73 |
| D-W5-05 | When the THFILE default applies | Only with MOTION `<fopt>` = 0 (the files hold Δt on line 1) and SITE `<delt>` = 0.005 s (relative 1e-6); otherwise there is no default and Error 73 names the reason. MOTION with `<out>` = 1 and no TH-to-RS conversion reads no history: no default is written or reported | a record at another Δt would be silently stretched |
| D-W5-06 | SOIL input motion | SOILX `<file>`, else THFILE, else `@rg160h_030g.acc` (SOIL is a seismic site response whatever ANALYS `<type>`), with the Δt condition of D-W5-05. A built-in record holds Δt on line 1: SOIL with `<header>` = 0 skips that line (read as 1) and checks that it equals `<delt>`; the rule is reported like a default | SOIL module tests |
| D-W5-07 | EQUAKE target spectrum | With no RSIN file at all and `<accopt>` ≠ 2, RSIN 1 = `@rg160h_030g.rsi` (RG 1.60 H, 0.30 g, 5 %, 27 rows); a blank Number of Frequencies takes its 27 rows (the existing "records of RSIN 1" rule); a given `<nrfreq>` ≠ 27 stays Error 89 with a hint. RSIN 2 / 3 given with RSIN 1 blank: no default (the user chose the components) | ACS: Error 84 |
| D-W5-08 | EQUAKE output names | For every spectrum that runs, a blank RSOUT is `<model>_eq<i>.rso` and (`<accopt>` ≠ 2) a blank ACCOUT `<model>_eq<i>.acc`, in the model folder. ACCIN of `<accopt>` 1 / 2 (Error 88) and TPSD have no default | ACS: Errors 86 / 87 |
| D-W5-09 | DYNP labels from the library | A SPRO label that no model DYNP defines resolves to the library curve of that label (Clay, Sand, Rock; case-sensitive) without INP; a model DYNP of the same label always wins (an incomplete one keeps Errors 97–99). AFWRITE writes the curve to the SOIL deck in FILE73 order; PIN ICURVE labels, SOILPROPPLOT and the GUI follow | ACS: Errors 97 / 99 |
| D-W5-10 | Default SOIL profile | With **no** SPRO entry at all: sublayer k = TOPL layer k, labelled Rock when its Vs ≥ 760 m/s (2493.4 ft/s; unit system of SOIL `<grav>`, or the HOUSE gravity while SOIL is not given) and Sand otherwise, then the SITE half-space `<hs>` without a label (linear) as the last sublayer. It needs TOPL layers and a defined `<hs>`; otherwise Error 95 names the reason. SPRO entries given without labels stay Error 95. 760 m/s is the ASCE 7 site class B/C boundary: a reported starting point that a site-specific profile replaces | ACS: Error 95 |
| D-W5-11 | Reporting (never silent) | Every default used is reported: CHECK **Warning EDU-29** "Built-In Default Input Used" under each module that uses it, naming the blank input and what replaces it (a warning: it never blocks AFWRITE); AFWRITE notes; RUN`<MODULE>` writes the defaults after the listing header (and so on the screen); the modules note every built-in input they read (`built-in input @x: ...`); `LIBRARY,DEFAULTS` lists them per AOPT module | `tests/unit/test_input_defaults.py` |
| D-W5-12 | ACS behaviour | `EDUOPT,DEFAULTS,ON\|OFF` (default ON). OFF: no defaults, blank inputs are Errors 73, 84, 86, 87, 95, 97, 99 as in ACS SASSI; `@` names still work. The CHECK catalogue cases of these errors run with OFF | `tests/unit/test_check_catalogue.py` |
| D-W5-13 | Decks | AFWRITE writes the default into the deck (Shown = used, like `<nrfreq>` and SITE `<freq2>`); a built-in file is written as its portable `@` name, which every module resolves. A module run from a hand-written deck with a blank input keeps the ACS module errors: the defaults are a preprocessor policy | module tests unchanged |
| D-W5-14 | GUI | A blank file field shows its default as placeholder ("built-in: RG 1.60, 0.30 g (@rg160h_030g.acc)"), computed with the rules above on the values being edited (type, `<fopt>`, Δt, `<accopt>`, spectrum number); file fields have a Library button listing the built-in files that fit (records for THFILE / ACCIN / SOILX, loads for THFILE, spectra for RSIN, PSDs for TPSD); the soil-curve field and Select Dynamic Soil Property list the library curves; a model without SPRO shows the default profile on the layer pages, stored with the first page changed (as the implicit wave field) | `tests/unit/test_input_defaults.py` |


### 7.20 Lead decisions: display aids, the soil island and deformed shapes from the results (wave 6, binding)

The user expected to see a soil island. SASSI has none: the layered site (TOPL / L layers on the SITE
half-space) is horizontally infinite and enters the analysis through the impedance at the interaction nodes;
the only soil elements of a model are the excavated soil. An engineer used to a direct (box-of-soil) model
still wants to see the model in its site, so SASSI-EDU draws the free-field profile around the foundation as a
picture, and says on the plot that it is one. The user also found that deformed shapes did not load: an
animation plays frames, the modules write frames only on a restart option and no example asks for them, so
after an analysis Load Frame Data was empty; a deformed shape needed `HARMFRAME` and `PROCFRAME` typed by
hand (only the lessons' Animate buttons did it).

| ID | Subject | Decision | Evidence |
|---|---|---|---|
| D-W6-01 | Command | New extension command `SHOWSOIL,[opt],[cut],[margin],[depth]` (P1, full name, class ui): `opt` -1 toggle (default), 0 off, 1 on; `cut` -1 automatic (cut when the foundation is embedded), 0 none, 1 cut; `margin`, `depth` lengths, 0 = automatic. A blank field keeps its value. A 3D view setting like `SHOWMASS` (`View3D.show_soil`, `soil_cut`, `soil_margin`, `soil_depth`): it acts on the active MODELPLOT / NODEPLOT, else becomes the default of new 3D plots; the 3D Plot toolbar has a button. With a plot it prints what is drawn and that it is a display aid | user request; `tests/unit/test_soil_island.py` |
| D-W6-02 | Geometry | `sassi.plotting.state.soil_island`: the foundation is the plan box of the interaction nodes (else of the element nodes at or below grade, else of all element nodes) and its lowest point; the top of the soil is the ground elevation (HOUSE `<gelev>`), each layer of `layer_table` at its depth, the half-space below. Automatic extent: margin = max(B/2, 0.3 × depth drawn, embedment); depth = the whole profile plus a half-space band of a quarter of the profile depth (a profile deeper than max(2.5 B, D + 1.5 B) is cut there, the plot says how many layers are not drawn); B = plan width, D = embedment. The plan box down to the lowest point is left open for the model (excavated soil, basement); the soil top under a surface foundation is not drawn (the model's face is there). The cut removes the quarter facing the viewer (sign of the view direction in X and Y; +X, -Y in the default view) through the centre of the foundation, down to the bottom. Only the boundary faces of the soil cells are drawn, with their outward normals; outline and crease lines, and the layer interfaces as thin lines | `tests/unit/test_soil_island.py` (face count and closure, opening, cut, depths, extents) |
| D-W6-03 | Colours | One colour per layer from light sand (smallest Vs of the profile) to dark brown (largest), neighbouring layers alternating 7 % in shade so that equal sublayers stay distinct; the half-space grey. Hover (GUI) gives the layer, its depths, thickness, Vs, Vp, unit weight and damping; a click prints them | `sassi/ui/static/plots.js` `soilTraces` |
| D-W6-04 | Fit and rendering | The plot box, the rotation centre and the camera fit include the soil (toggling re-fits). GUI: one plotly mesh per layer; the node plot draws the soil see-through. Headless renderer (CAPTUREPLOT, batch): the soil faces that face the viewer join the model faces and the beams (as screen-plane ribbons) in one depth-sorted collection; interaction markers inside the soil are not drawn | `tests/unit/test_soil_island.py` (PNG rendering) |
| D-W6-05 | Deformed shape from the results | Load Frame Data of Plot > Deformed Shape also offers the steady-state motion at one computed frequency of the active model's FILE8-type files (`FILE8[A-Z0-9_]*` in the model folder that read as FILE8 containers), total or relative to the free field (seismic only). **Animate** submits command text (L17): `HARMFRAME,<file>,<f>,HARM_<f>[,,0]` (folder `HARM_<f with p>[R]`), then the lines of a lesson's `animate` action for that folder (`PROCFRAME` into `<folder>_ani`, `DEFORMPLOT` with the automatic scale 0.15 × model size / largest displacement, the undeformed shape, the title). API: `GET /api/harmonic`, `POST /api/harmonic/plan`, `POST /api/harmonic/show` (`sassi/ui/harmonic.py`) | user report; `tests/unit/test_harmonic_dialog.py` |
| D-W6-06 | Preselected frequency | The computed frequency of the largest node-to-node spread of the motion, max over nodes of the absolute value of H_n - mean_n H, in the control direction (seismic) or per translation (vibration). The largest absolute H would pick the rigid translation at the lowest frequency for an embedded box that never exceeds the free field; the spread picks 3.49 Hz for example 1 (the lessons' SSI frequency), 3.98 / 4.37 Hz for example 9 SSI / fixed base (its first frequencies) and 12 Hz for example 2 (FV) | `sassi/ui/harmonic.py` `_summary` |
---

## 8. Implementation order and acceptance gates (informative)

| Milestone | Content | Gate |
|---|---|---|
| M0 Core | interpreter (§3.1–3.3), model data (§3.4.A–E), FREQ, CHECK catalogue, AFWRITE decks, WRITE/INP/SAVE/RESUME, RUN framework and batch protocol, file containers | UT-01…UT-19, VP-52 (model checks) |
| M1 Free field | SOIL (SHAKE), SITE Mode 1/2 (vertical waves), EQUAKE (LW + checks), NJ spectra | VP-02a, VP-02b, VP-04, VP-05, VP-06, VP-30, VP-36 |
| M2 Elements | HOUSE element library (SOLID, BEAMS, SHELL, SPRING, GENERAL, masses), DOF management, FIXROT family | VP-01, VP-03, VP-37, VP-38 (P0 part), VP-39 |
| M3 Impedance and SSI | POINT3, ANALYS (F_ff, X_ff, assembly, Schur solve, restarts, simultaneous cases, global impedance), FORCE, COMBIN | VP-07…VP-11, VP-13…VP-17, VP-22…VP-25, VP-40, VP-41 |
| M4 Post-processing | MOTION (options 0–6, smoothing, phase adjustment, RS, baseline), RELDISP, STRESS | VP-28…VP-31, VP-34, VP-35 |
| M5 GUI shell | main window, menus, Options dialogs, Check window, command history, editor | UI-01…UI-07 behaviours, round trip through dialogs |
| P1 | INTGEN/ETYPEGEN/EXCAV/MERGESOIL…, macros and variables, cuts and section calculations, plots and animations, line maths, ANSYS import/export, incoherency and ME, 2D, SYMM, nonlinear-soil iterations, CRITFREQ, AFWRBAT, EQUAKE AB/correlation | VP-12, VP-18…VP-21, VP-26, VP-27, VP-32, VP-33, VP-42, VP-43, VP-47…VP-51, VP-53, UT-09…UT-11, UT-20, UT-21 |
| P2 | Option NON, SOIL-NON, water, binary databases, TSHELL, Option A/AA, cosmetic UI | VP-44…VP-46, VP-54 |

Each milestone ships with its VP reports (`VERIFY,P0` output) archived in `docs/verification/`.
