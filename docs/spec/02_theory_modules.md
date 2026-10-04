# Spec 02 — Theoretical Basis (Ch. 2) and SSI Module Descriptions (Ch. 3)

**Source:** ACS SASSI Version 3 User Manual, extracted text lines 1473–2432
(`reference/acs-sassi.txt`), printed pages 33–54 (PDF pages ≈ 35–56 in the
Read-tool numbering; printed page = PDF page − 2).
Images transcribed by viewing the PDF: Figure 2.1, Figure 2.2, Eq. (2.1), Table 3.1,
Table 3.2. Figure 1.1 (module chart, printed p. 13, outside this range) was also viewed
because Chapter 3 depends on it; its file-flow arrows are summarised in §3.0.

**Tag legend**

| Tag | Meaning |
|---|---|
| [CORE] | Needed for computational correctness |
| [IO] | File / input-output format or naming |
| [UI] | User interface, plotting, convenience |
| [ADV] | Advanced: incoherency, nonlinear (Option NON), Option A/AA, Option PRO |
| *(impl. note)* | Not stated in this section of the manual. It is the implementer's interpretation or standard SASSI theory (Lysmer et al. 1981; Tabatabaie‑Raissi 1982; Tajirian 1981) and has to be checked against later sections. |

---

## Part 1 — Chapter 2: Theoretical Basis for Linearized SSI Analysis

### 1.1 Flexible Volume (FV) substructuring: concept [CORE]

- The SSI solution is computed for a **coupled structure + excavated-soil system**.
  Its dynamic stiffness is the *difference* between the complex dynamic stiffness of
  the structure (including basement) and that of the excavated soil, which is
  *subtracted*.
- Physical interpretation, which also gives a verification test (see §5): if the
  structure's complex dynamic stiffness equals the excavated soil's at every
  excavated-soil node, the difference is zero and there is **no SSI effect**. The
  larger the difference, the larger the SSI effect on an embedded structure.
- The SSI system is split into **three coupled substructures**:
  1. the **free-field system** (horizontally layered soil on a half-space, with no
     structure and no excavation);
  2. the **structural system** (FE model: structure + basement + optional near-field soil);
  3. the **excavated soil system** (FE model of the soil volume that the embedded
     structure displaces, with the free-field soil properties).
- The "FE model part" is made of two FE models: (i) the structural system and
  (ii) the excavated soil system.

**Figure 2.1 "Flexible Volume Method Concept"** (transcribed from the image). Title
in the figure: *Flexible Volume SSI Substructuring*. Text bullets in the figure:
- "No wave scattering analysis."
- "Impedance problem is trivial; reduced to a free-field problem."
- "Structural dynamic problem slightly more complex since includes a coupled excavated soil."

Other labels: "No Mass Structure" (a grid of interaction points in the layered soil)
with the note "Each load case solved fast using axisymmetric soil model". The right
side shows a "Structure Minus Excavated Soil" FE model attached by springs and dashpots
at the interaction nodes, the "Excavated Soil Model", and the seismic load term
written as **X_f ü_f** (impedance times free-field motion).
→ Each interaction-node point-load case is solved with an axisymmetric (3D) or
plane (2D) soil column model. That is the POINT module (§3.4).

### 1.2 DOF sets and Figure 2.2 [CORE]

Subscripts used in Eq. (2.1):

| Subscript | Node set | Meaning |
|---|---|---|
| `s` | superstructure | Structure nodes that are not shared with the excavated soil (above grade, or internal structure nodes not on excavated-soil nodes) |
| `i` | basement | Nodes shared by the structure (basement) and the excavated soil |
| `w` | excavated soil | Excavated-soil nodes that are *not* connected to the structure (internal nodes of the excavated volume) |
| `f` | interaction nodes | The nodes where the impedance X_ff is applied and the free-field motion U′_f is evaluated. Which nodes are in `f` depends on the method (§1.3). |

**Figure 2.2 "The Flexible Volume SSI Substructuring Approaches"** (transcribed):
- **(a) The FV Method (Direct):** *Free Field Problem* (layered soil; **all**
  excavated-volume nodes, both `w` interior and `i` boundary, are drawn as
  interaction nodes) **+** *Structure* (mesh with `s` nodes; `i` nodes on the
  embedded basement walls and base) **−** *Excavated Soil* (`w` interior nodes and
  `i` boundary nodes, all interaction nodes) **=** *SSI Problem (FV)*. The seismic
  input enters through the free-field problem.
- **(b) The FI-EVBN Method (Modified Subtraction):** the same decomposition, but
  in the free-field problem only the **boundary nodes** of the excavated volume are
  interaction nodes. This boundary includes the top row at the ground surface; the
  interior is drawn as a dotted, non-interaction grid. The excavated soil shows
  interior `w` nodes as plain mesh intersections, which are not interaction nodes,
  and a purple boundary ring of `i` nodes.

### 1.3 Method variants [CORE]

| Method name (ACS SASSI) | ASCE 4-17 name | Interaction node set `f` | Accuracy / cost |
|---|---|---|---|
| **FV** (Flexible Volume, reference, "Direct") | Direct | **All** excavated-soil nodes (`i` ∪ `w`) | Reference solution; most expensive (largest flexibility and impedance matrices) |
| **FI-FSIN** (Flexible Interface, Foundation-Soil Interface Nodes) | Subtraction | Only nodes on the **foundation–soil interface** (embedded walls and base) | Fastest; least accurate. *(impl. note: the literature reports spurious frequency spikes)* |
| **FI-EVBN** (Flexible Interface, Excavation Volume Boundary Nodes) | Modified Subtraction | Foundation–soil interface nodes **plus the nodes on the ground-surface (top) face of the excavated soil** | Much more accurate than FSIN. The ground-surface nodes are "of key importance" for capturing scattering of surface waves. |
| **FFV** (Fast Flexible Volume, Ghiocel 2013a) | — | EVBN set **plus one or more internal layers** of excavated-soil nodes | More accurate than the FI methods for deeply embedded structures (e.g. SMRs). The internal layers have to be chosen case by case, using expert judgment and/or sensitivity runs with alternative selections. |

Rules:
- The FI methods approximate FV so that the soil impedance calculation costs less for
  embedded structures. Their interaction nodes lie only on the excavation lateral
  surface (FSIN or EVBN), so they are much faster than FV.
- The **SSI matrix formulation is identical** for FV, FI and FFV. Only the excavated-soil
  equations differ, depending on which nodes are interaction nodes.
- For FI methods, the equations of motion of the **internal** (non-interaction)
  excavated-soil nodes are **approximate**: they omit the soil impedance terms and the
  free-field seismic load terms (X = 0 and no U′ for those nodes).
- The size of the flexibility and impedance matrices (§1.6) depends on the method:
  n_f × n_f, where n_f = number of interaction DOFs.

### 1.4 The SSI matrix equation — Eq. (2.1) (exact transcription) [CORE]

$$
\begin{bmatrix}
\mathbf{C}^{s}_{ii}-\mathbf{C}^{e}_{ii}+\mathbf{X}_{ii} & -\mathbf{C}^{e}_{iw}+\mathbf{X}_{iw} & \mathbf{C}^{s}_{is}\\[2pt]
-\mathbf{C}^{e}_{wi}+\mathbf{X}_{wi} & -\mathbf{C}^{e}_{ww}+\mathbf{X}_{ww} & \mathbf{0}\\[2pt]
\mathbf{C}^{s}_{si} & \mathbf{0} & \mathbf{C}^{s}_{ss}
\end{bmatrix}
\begin{Bmatrix}\mathbf{U}_i\\ \mathbf{U}_w\\ \mathbf{U}_s\end{Bmatrix}
=
\begin{Bmatrix}
\mathbf{X}_{ii}\mathbf{U}'_i+\mathbf{X}_{iw}\mathbf{U}'_w\\
\mathbf{X}_{wi}\mathbf{U}'_i+\mathbf{X}_{ww}\mathbf{U}'_w\\
\mathbf{0}
\end{Bmatrix}
\qquad (2.1)
$$

- Superscript `s`: structure. Superscript `e`: excavated soil.
- **U** = complex nodal displacement amplitudes (total motion). **U′** = free-field
  displacement amplitudes at the same nodes, from the SITE module.
- **X** (also written **X_ff**) = the frequency-dependent **impedance matrix**. It is
  the dynamic boundary of the foundation at the interaction nodes, i.e. the dynamic
  stiffness of the free-field layered soil at those nodes.
- Eq. (2.1) gives the final **total** motions of the structure.

**Eq. (2.2), dynamic stiffness:**

$$\mathbf{C}(\omega)=\mathbf{K}-\omega^{2}\mathbf{M}\qquad(2.2)$$

**M** = total mass matrix. **K** = total **complex** stiffness matrix, where material
damping enters as complex moduli (§1.5.1 item 13 of the manual: frequency-independent
damping that can vary from element to element). ω = circular frequency (rad/s).

#### 1.4.1 Specialisation by method *(impl. note, derived from §1.3)*

- **FV:** `f = i ∪ w`. X is a full n_f × n_f matrix over all `i` and `w` DOFs. U′ is
  needed at every `i` and `w` node.
- **FI-FSIN / FI-EVBN / FFV:** interaction DOFs `f` ⊂ `i ∪ w`. For every DOF not in
  `f`, the X rows and columns are zero and its U′ term is dropped. Eq. (2.1) then
  becomes the excavated-soil subtraction with impedance only at `f`. If some `w`
  nodes are not interaction nodes, their rows read
  `−C^e_wi U_i − C^e_ww U_w = 0` (+ X terms only for those in `f`).
- Eq. (2.1) as printed has no structural stiffness on `w` nodes (zero `s`–`w`
  coupling). Real models may have internal slabs or walls passing through
  excavated-soil nodes. The **general assembly** below covers both cases, and the
  implementation should use it.

#### 1.4.2 General assembly form (recommended implementation) *(impl. note)*

For each frequency ω, over the global DOF vector of the HOUSE model:

```
A(ω) = Σ_struct-elems (K_e* − ω² M_e)          # structure + basement + near-field soil
     − Σ_excavated-soil-elems (K_e* − ω² M_e)   # excavated soil, free-field properties
A[f,f] += X_ff(ω)                                # impedance on interaction DOFs only
P_seismic[f] = X_ff(ω) · U'_f(ω)                 # all other entries 0
P_force      = nodal external forces (FILE9)      # forced vibration
A · U = P                                         # complex symmetric system
```

Permuting the DOFs into (i, w, s) order and assuming no structure on `w` gives exactly
Eq. (2.1). The matrices are **complex symmetric**, not Hermitian. Use complex
symmetric (LDLᵀ) or general LU solvers, never Cholesky or Hermitian solvers.

### 1.5 Three main computational steps, per frequency [CORE]

1. **Site response problem:** compute the free-field motion **U′_f** in the embedded
   part of the structure (SITE).
2. **Impedance problem:** compute the free-field impedance matrix **X_ff**
   (POINT + ANALYS).
3. **SSI problem:** form the full SSI complex stiffness and load vectors and solve
   Eq. (2.1) for the final SSI displacements (HOUSE + ANALYS).

### 1.6 Site response analysis (Section 2.2) [CORE]

- Site model: **horizontal soil layers over a uniform half-space**. All materials are
  **viscoelastic**.
- Layer shear stiffness and hysteretic damping are normally the **iterated
  equivalent-linear values** from a nonlinear free-field analysis.
  - The *linearized* site response for SSI is done in **SITE**.
  - The *nonlinear* (equivalent-linear iteration) site response that produces those
    properties is done in **SOIL**, using the SHAKE equivalent-linearization method.
- Only the free-field displacements at the **layer interfaces where the structure
  is connected** are needed. *(impl. note: every interaction node must lie at a soil
  sublayer interface elevation. The excavated-soil mesh has to be compatible with the
  SITE layering. Validate this and report an error otherwise.)*
- **Eq. (2.3):**

$$u'_f(x)=U'_f\,\exp\!\left[i(\omega t-kx)\right]\qquad(2.3)$$

  - **U′_f** = "mode shape" vector of the interface amplitudes at and below the
    **control point** (x = 0).
  - **k** = complex horizontal wave number. It sets how fast the wave propagates and
    decays in the horizontal x direction. *(impl. note: for vertically propagating
    body waves k = 0, so the motion is uniform in x. For inclined body waves
    k = ω sin θ / V, which is constant through the layers by Snell's law. For surface
    waves, k comes from the layered-medium eigenproblem.)*
- Wave types, with mode shapes and wave numbers found by "effective discrete and
  iterative eigen solution methods" (*impl. note: e.g. the thin-layer method of
  Lysmer/Waas/Kausel*), for control motion at **any soil layer interface**:
  - inclined **P**, **SV**, **SH** body waves;
  - **Rayleigh** and **Love** surface waves.
- Two wave combinations exist:
  1. **SV + P + Rayleigh:** in-plane, **2 DOF per node**.
  2. **SH + Love:** out-of-plane, **1 DOF per node**.
- A complete recorded wave field contains all five types. For typical seismic SSI
  (ASCE, USNRC), use **only vertically propagating SV, SH and P**.
- **Recommended practice:** **SV waves for X**, **SH waves for Y**, **P waves for Z**
  (to rebuild a 3D field from vertically propagating waves).
- [ADV] Incoherent motion option (seismic random field) and **wave passage** effects
  can be included. References: EPRI 2006–2007 reports (Short et al. 2006, 2007;
  Ghiocel 2007a) and Ghiocel 1998, 2007b, 2009a/b, 2013b, 2014a, 2015b.

### 1.7 Impedance analysis (Section 2.3) [CORE]

- The soil impedance matrix is the dynamic stiffness of the free-field layering at the
  interaction nodes.
- **X_ff = F_ff⁻¹**: the inverse of the free-field dynamic soil **flexibility
  (compliance)** matrix at the interaction nodes.
  *(Cross-ref: manual Eq. 4.3, printed p. ~70: K + iD = (f + ig)⁻¹.)*
- The matrix size depends on the method chosen (FV, FI or FFV).
- Section 2.5 step 3: the flexibility matrix is **inverted in place** with a routine
  for **symmetric** matrices, and the impedance is stored in the same (symmetric)
  form.

### 1.8 Structural analysis (Section 2.4) [CORE]

- The FE part can include the **structure**, optionally some surrounding **backfill
  soil**, and the **excavated soil**.
- **Near-field soil** = backfill soil included in the FE part.
- **Far-field soil** = soil layering outside the FE part. It connects to the
  excavated-soil FE part at the outer-surface interaction nodes.
- The backfill zone can be modelled in **2D with PLANE** elements or in **3D with
  SOLID** elements.

### 1.9 Per-frequency solution algorithm (Section 2.5) [CORE]

For each frequency ω, form the impedance matrix, the load vector, and the dynamic
stiffness C = K − ω²M for (a) the structure, including basement and any irregular
soil zone, and (b) the excavated soil. Large 3D embedded problems must be handled with
the available RAM and disk in mind (out-of-core where needed).

Operations per SSI frequency step:

1. **Form complex dynamic stiffness of structure** from its total stiffness and
   mass matrices.
2. **Form complex dynamic stiffness of excavated soil** from its total stiffness and
   mass matrices.
3. **Form impedance matrix.** For the flexible volume or interface methods, build the
   full flexibility matrix for the interaction nodes, then invert it in place
   (symmetric routine) to get the impedance matrix, stored in the same form.
4. **Form total stiffness of the soil–structure system.** Add the impedance matrix
   and **subtract** the excavated-soil dynamic stiffness.
5. **Form load vector.**
   - *Seismic:* impedance matrix × free-field motion vector (P_f = X_ff U′_f).
   - *Forced vibration:* taken directly from the given nodal external forces.
6. **Solve the linear complex system.** The result is the **transfer function
   matrix**:
   - **ATF**, absolute acceleration transfer functions, for seismic analysis;
   - **DTF**, displacement transfer functions, for forced vibration.

*(impl. note — ATF normalisation.)* SITE computes U′_f for a **unit control-motion
amplitude** in the selected direction, so the solved U is the transfer function from
control motion to total motion. Because ü = −ω²U for both the response and the control
point, the displacement ratio equals the **absolute acceleration** ratio, which is the
ATF. For forced vibration, U per unit force is the DTF. Several right-hand sides can
be solved at once: three directions, up to 50 incoherent simulations, or up to 500
load cases (§3.7).

---

## Part 2 — Chapter 3: Description of SSI Modules

Options A and AA are documented in the separate "ACS SASSI-ANSYS Integration Capability"
manual. Option PRO (probabilistic site response and SSI) has its own manual. Option
NON modules are included in this chapter.

### 3.0 Module inventory, file flow and run order

**Main SSI modules:** EQUAKE, SOIL, SITE, POINT (POINT2/POINT3), HOUSE, FORCE, ANALYS,
MOTION, STRESS, RELDISP, COMBIN.
**Option NON:** NONLINEAR, COMB_XYZ_THD (the manual also spells it COMBIN_XYZ_THD and
COMBINE_XYZ_THD). Auxiliary: COMB_XYZ_STRAIN.
**Executable variants:** HOUSEFS / ANALYSFS (fast solver); HOUSEFSA / ANALYSFA
(Option AA).

**Figure 1.1 file flow** (cross-ref, printed p. 13). Red ovals are text files and
green ovals are binary files in the original figure.

| From | File | To |
|---|---|---|
| EQUAKE | acceleration file (.acc) | SOIL (and MOTION, RELDISP, STRESS as control motion) |
| SOIL | FILE88 | SITE |
| SOIL | FILE73 | STRESS (nonlinear soil curves) |
| SITE | FILE2 | POINT |
| SITE | FILE1 | ANALYS |
| POINT | FILE3 | ANALYS |
| HOUSE | FILE4 (modelname.n4), FILE77 | ANALYS (FILE4 also STRESS) |
| HOUSE | FILE78 | STRESS |
| STRESS | FILE74 | HOUSE |
| FORCE | FILE9 | ANALYS |
| ANALYS ↔ disk | restart files COOXxxx, COOTKxxx, DOFSMAP, FILE90, FILE91 | ANALYS (restart) |
| ANALYS | FILE8 | COMBIN, MOTION, RELDISP, STRESS |
| MOTION | FILE12, FILE13; *.tfu, *.tfi, *.acc, *.rs | user |
| MOTION | complex TFI (.TFI) | RELDISP |
| RELDISP | *.thd, *.tfd | user, NONLINEAR, COMB_XYZ |
| STRESS | FILE14, FILE15; *.tfu, *.tfi, *.ths | user |
| NONLINEAR/COMB_XYZ | → | HOUSE (next iteration) |
| LOADGEN [ADV, Option A] | export load to ANSYS model | — |
| SSI2ANSYS [ADV, Option AA] | import ANSYS FE model | HOUSE |

**Module ↔ file I/O matrix** (Ch. 3; ✱ = cross-ref from other chapters):

| Module | Input deck (by AFWRITE) | Reads | Writes |
|---|---|---|---|
| EQUAKE | `.equ` | target RS `.rsi` (2 cols); optional seed records; optional target PSD; or external accelerograms | `.acc`, `.vel`, `.dis`, `.rso`, `.psd`, `.fft`, output file |
| SOIL | `.soi` | input accelerogram | output; `ACCxxx.TH`, `SNxxx.TH`, `SSxxx.TH`; `FILE73`; `FILE88` |
| SITE | `.sit` | (FILE88 if strain-compatible/nonlinear option) | `FILE1`, `FILE2`, output |
| POINT2/POINT3 | `.poi` | `FILE2` | `FILE3` |
| HOUSE (HOUSEFS/HOUSEFSA) | `.hou` | `.sit` (embedded models); `FILE74` (nonlinear iter.); `.liq`, `.pin` (nonlinear soil); ANSYS matrix files (Option AA) | `FILE4` (=`modelname.n4`), `COOSK`, `COOSM`, `FILE77`/`FILE77xxx`, `FILE78`, `.hounew`, `.map`, `.liq`; ✱`DOFSMAP`, `FILE90`, `FILE91` |
| FORCE | `.frc` | — | `FILE9` (user copies to `FILE9001`…`FILE9500`) |
| ANALYS | `.anl` | `FILE1` (or `FILE1X/Y/Z`), `FILE3`, `FILE4`(+`COOSK`,`COOSM`); `FILE9`/`FILE9xxx` (forces); `FILE77`/`FILE77xxx` (incoherent); restart: `COOXxxx`, `COOTKxxx`, `COOXI`, `COOTKI`, `DOFSMAP`, `FILE90`, `FILE91` | `FILE8` (or `FILE8X/Y/Z`, `FILE8xxx`); restart files (if "Saving Restart Files"); `FILE11` (global unconstrained impedance) |
| MOTION | `.mot` | `FILE8`; control motion; `CONTTRS` list (external RS mode) | output; `FILE12` (opt.), `FILE13`; `*.TFU`, `*.TFI`, `*.ACC`, `*.RS`; frames `\TFU`, `\RS`, `\ACC`, `\ACCR`; `ACC_max.txt`; `Modelname_ACC.bin`; `*.RSO`; ✱`SRSSTF.txt` |
| STRESS | `.str` | `FILE4`, `FILE8`; `FILE78` (nonlinear); ✱`FILE73`; `STATIC_SOIL_PRESSURES.TXT` | output; `*.TFU`, `*.TFI`, `*.THS`; `Modelname_STRESS.bin`; `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`; frames `\NSTRESS`, `\SOILPRES`; `STATIC_SOIL_PRESSURES.TXT`; `pres_max_ele`, `pres_max_nod`; `ESTRESS_framenumber.ess`; `FILE74`; ✱`FILE14`, `FILE15` |
| RELDISP | `.rdi` | `.TFI` files from MOTION | output (max relative displacements); `*.TFD`, `*.THD`; frames `\THD`, `\THDR`; `Modelname_TR_y_THD.bin` |
| COMBIN | — | `FILE81`, `FILE82` | `FILE8` |
| NONLINEAR [ADV] | — | `PANEL.NON` / `SPRING.NON`; previous-iteration results | (see Option NON spec) |
| COMB_XYZ_STRAIN [ADV] | — | directional `FILE74` | combined `FILE74` (SRSS) |
| COMB_XYZ_THD [ADV] | `COMB_XYZ_THD.inp` ✱ | X/Y/Z displacement histories | combined co-directional histories |

---

### 3.1 Module EQUAKE — spectrum-compatible accelerogram generation

**Purpose [CORE]:** generate acceleration time histories compatible with design ground
response spectra (RS), and characterise both simulated and user-supplied histories.

**Algorithm [CORE]:**
1. **Frequency-domain matching**, Levy–Wilkinson (**LW**) algorithm, run first to get
   an approximate RS match.
2. **Time-domain matching**, Abrahamson (**AB**) algorithm (as in RspMatch), run next
   to bring the computed RS closer to the target.
3. **Phasing:** either (a) **random phases**, or (b) **"seed records"** (ASCE 4-16
   approach), where the simulated motions **keep the Fourier phasing of the seed
   record** components for X, Y and Z.
4. **Baseline correction:** the FLUSH-code complex-frequency algorithm, followed by
   extra **polynomial corrections in the time domain**.

**SRP 3.7.1 compliance checks** (single time history, Option 1, Approach 2) [CORE]:

| Criterion | Rule | On failure |
|---|---|---|
| Duration | Total duration ≥ **20 s** | Warning on screen and in output file |
| Time step / Nyquist | Time step **Δt ≤ 0.005 s**, i.e. Nyquist ≥ 100 Hz (the manual's wording is garbled; see open questions) | Warning on screen and in output file if Δt > 0.005 s |
| RS check frequencies | At least **100 points per frequency decade**, between the lowest and highest frequency in the target RS input file | — |
| Matching band | The computed **5 %-damping** RS must nowhere be **more than 10 % below** the target and nowhere **more than 30 % above** it | — |
| Low points | **No more than 9 adjacent frequency points** below the target | — |

**Feature parameters computed** [CORE]:
- Strong-motion duration = time between **5 % and 75 % Arias intensity**.
  *(impl. note: I_A(t) ∝ ∫₀ᵗ a² dτ.)*
- **V/A** and **AD/V²**, where A = PGA, V = PGV, D = PGD.
- With more than one component: **stationary** cross-correlation (whole duration) and
  **nonstationary** cross-correlation (**2-second moving window**) for each pair.
- If a target **PSD** is given, check the generated history against the PSD
  requirement.

**Outputs and formats** [IO]:
- Input deck: **`.equ`**, written by AFWRITE.
- Target RS input: a text file (extension **`.rsi`** in Demo 1) with **2 columns:
  frequency, amplitude**, for one damping ratio.
- Generated histories: text files **`.acc`**, **`.vel`**, **`.dis`**. The user names
  the acceleration file, and the same base name is used for the velocity and
  displacement files.
- **`.rso`**: computed response spectra.
- **`.psd`**: 2 columns, frequency and PSD amplitude. PSD uses **±20 % frequency
  averaging intervals** (ASCE 4, USNRC). Units depend on the chosen gravity units:
  ft/s² → **in²/s³**; m/s² → **cm²/s³**.
- **`.fft`**: 3 columns, frequency, Re(FT), Im(FT). **Positive frequencies only.**
- FFT and PSD are computed from the **strong-motion part only** (5–75 % Arias window).
- Output file contents: input data echo, parameters of the generated histories,
  stationary pair correlations (full duration), and nonstationary correlations
  (2-s window).
- EQUAKE also computes RS, PSD and FFT for **external** acceleration histories
  supplied by the user.

**Uses** [CORE]: the generated accelerograms feed SOIL (site response), MOTION and
STRESS (SSI responses).

**Nonstationary-correlation option** [ADV]: the nonstationary correlations of
recorded components (NS, EW, V) can be used to simulate histories with the same
correlation patterns (an alternative to using recorded phasing). They also give
insight into incoming wave patterns and principal axes.
> **WARNING:** This option is not in the ASCE standards or USNRC guides. If it is
> used, restore RS compatibility of the correlated components by feeding the
> generated correlated histories back as **"seed" records** and running EQUAKE again.

*(impl. note: SRP 3.7.1 also limits the correlation coefficient between components,
typically |ρ| ≤ 0.16. This section gives no threshold.)*

### 3.2 Module SOIL — free-field nonlinear site response

- [CORE] **SOIL-EQL:** nonlinear site response for **vertically propagating S waves**,
  using the **Seed–Idriss iterative equivalent-linear** model (SHAKE methodology plus
  later improvements). The equivalent soil properties it computes can then be used in
  the SSI analysis.
- [IO] Input deck **`.soi`** (AFWRITE).
- [IO] Plotting time histories, text files with extension **`.TH`**:
  - `ACCxxx.TH`: acceleration time history of soil layer xxx;
  - `SNxxx.TH`: strain time history of layer xxx;
  - `SSxxx.TH`: stress time history of layer xxx;
  - `xxx` = 3-digit free-field layer number, counted **from the ground surface
    downward** (e.g. `ACC001.TH` = layer 1).
- [IO] **`FILE73`** (text): material soil curves (G/Gmax and damping versus strain)
  used by **STRESS** for nonlinear SSI.
- [IO] **`FILE88`**: iterated equivalent-linear (effective) soil properties used by
  **SITE** when the user selects the nonlinear (strain-compatible) SSI option.
- [ADV] **SOIL-NON:** nonlinear **time-domain** site response with a **hyperbolic
  hysteretic** soil model, based on DEEPSOIL theory (University of Illinois at
  Urbana-Champaign). It needs minimal input (lower part of the SOIL input window).
  Output file names are the same as for the equivalent-linear option. Validated
  against DEEPSOIL in **V&V Problem 49**.
  > **WARNING:** Not for direct use on nuclear licensing projects unless the analyst
  > fully justifies it. It is not part of standard US practice for SSI inputs. It is a
  > benchmark tool for soft soils with strong nonlinearity, to compare against the
  > standard equivalent-linear procedure.
  >
  > **WARNING:** Valid only for **nonlinear convolution** analysis with the input
  > motion defined as the compatible motion **at bedrock**.

### 3.3 Module SITE — linearized site response

- [IO] Input deck **`.sit`** (AFWRITE). It must define the **control point** (layer
  interface) and the **wave composition** of the control motion.
- [CORE] Solves the site response problem and writes to **`FILE1`** the data needed
  to compute the **free-field displacement vector** (U′ at the interaction nodes, for
  unit control motion).
- [CORE] Writes **`FILE2`**: data for the **transmitting boundary** calculations,
  used by POINT.
- The control-motion **time history is not needed** in SITE. MOTION uses it later.
- ✱ Cross-ref (Ch. 4/6): SITE has **Mode 1**, which forms and solves the
  transmitting-boundary eigenproblem (Rayleigh/Love) and saves it to FILE2, and
  **Mode 2**, which reads FILE2 + `.sit` and solves the linearized site response,
  saving it to FILE1. Both are normally run.
- ✱ For simultaneous three-direction coherent ANALYS runs: run SITE three times and
  copy FILE1 to **FILE1X** (SV, x′ direction, angle 0), **FILE1Y** (SH, y′
  direction, angle 0) and **FILE1Z** (P, z direction, angle 0). The coordinate
  transformation angle in `.anl` must be 0.
- The manual's sentence "In addition to the output and binary files FILE1 and FILE2."
  is truncated. SITE outputs are taken to be: the output listing + FILE1 + FILE2.

### 3.4 Module POINT (POINT2 / POINT3) — soil flexibility data

- Two programs: **POINT2** (2D SSI) and **POINT3** (3D SSI).
- [IO] Input deck **`.poi`** (AFWRITE).
- [CORE] Computes the data needed to form the **frequency-dependent flexibility
  matrix**, i.e. point-load solutions of the layered soil column at the interaction
  nodes. 2D uses plane-strain elements and 3D uses axisymmetric elements
  (✱ Ch. 4 item 17). Results go to **`FILE3`**.
- [CORE] Requires **`FILE2`** from SITE, so **SITE must run before POINT**.
- ✱ Cross-ref: `.poi` also specifies the maximum embedment and the radius of the
  central zone (point-load radius).

### 3.5 Module HOUSE — structural element matrices

- [CORE] Forms the **element mass and stiffness matrices** of all elements (structure,
  near-field soil, excavated soil). The FE model can be the structure alone, or
  include near-field soil, especially where **irregular soil zones** exist.
- Names: **HOUSEFS** (fast-solver version), **HOUSEFSA** (Option AA, which reads the
  **ANSYS model matrix files** directly).
- [IO] Input deck **`.hou`** (AFWRITE).
- [IO] Outputs:
  - Baseline code with the standard solver: stiffness and mass go to **`FILE4`**
    (= **`modelname.n4`**).
  - Structure and basement stiffness and mass matrices are also stored in
    **`COOSK`** and **`COOSM`** (fast solver).
  - ✱ `DOFSMAP`, `FILE90`, `FILE91`, which ANALYS restart needs (Ch. 6, "Save
    Restart Files").
- [ADV] **Incoherency:** the coherence kernel random-field decomposition is done
  here. Results for incoherent SSI or **nonuniform/multiple support excitation** go to
  **`FILE77`**, used by ANALYS. HOUSE can print the incoherent mode contributions to
  check the decomposition accuracy. Several stochastic simulations produce several
  FILE77s (✱ `FILE77001`…`FILE77050`, up to 50). **Each FILE77 holds X, Y and Z.**
- [ADV] **Nonlinear SSI:** reads **`FILE74`** (from STRESS). Writes **`FILE78`**,
  which is non-empty **only** when the nonlinear SSI option is used and is read by
  STRESS during nonlinear iterations.
- [CORE] HOUSE can run independently of SITE and POINT **except for embedded
  models**. In all cases it requires the **`.sit`** file in the working directory.
- [CORE/UI] **Node-numbering optimizer** (bandwidth/profile reduction):
  - Strongly recommended for large embedded models. HOUSEFSA (Option AA) always
    uses it.
  - Selected in **UI Options/Analysis/HOUSE** input window.
  - After a run, HOUSE writes **`modelname.hounew`** (new HOUSE input with optimised
    numbering) and **`modelname.map`** (original node → new node pairs).
  > **WARNING:** With the optimizer, use the **new** node numbers to select MOTION and
  > RELDISP output nodes. The `.hounew` model is what generates FILE4 (`.n4`), COOSK
  > and COOSM. **FILE8** (used by MOTION, RELDISP, STRESS) is in the **renumbered**
  > numbering.
  >
  > **WARNING:** Option AA always uses the optimizer and writes `modelname.map` with
  > new/old node pairs.

### 3.6 Module FORCE — external load vectors

- [IO] Input deck **`.frc`** (AFWRITE).
- [CORE] Forms the load force vector for external load cases: impact, rotating
  machinery, or **unit forces** used to get the impedance of a flexible foundation.
  Not used for seismic problems except for computing foundation impedances.
- [IO] Output **`FILE9`**.
- [IO] **Multiple load cases** in one ANALYS run, **up to 500** if RAM allows (this
  saves a lot of run time): after each FORCE run, copy FILE9 to `FILE9` + load case
  number, i.e. **`FILE9001`, `FILE9002`, …, `FILE9500`** (3-digit zero-padded
  suffix; see open question on the "one digit" wording).

### 3.7 Module ANALYS — SSI solution per frequency

- [IO] Input deck **`.anl`** (AFWRITE).
- [IO] Always required: **`FILE1`**, **`FILE3`**, **`FILE4`**. Also **`FILE9`** for
  external load cases and **`FILE77`** for incoherent analysis.
- [CORE] Computational steps:
  1. Form the soil **flexibility** matrix (from FILE3 + FILE2 data).
  2. Compute the soil **impedance** matrix at the interaction nodes (inversion).
  3. Form the **seismic load vector**, including incoherency effects when they are
     considered.
  4. **Solve** the system at each frequency step and get the complex response
     **transfer functions** for every DOF.
- [CORE/IO] Output **`FILE8`**: complex TFs, either control motion → final (total)
  motions (seismic) or external loads → total displacements (forced). Used by
  **MOTION, RELDISP, STRESS**. Interpolation of the TFs and other output requirements
  are handled downstream (MOTION/STRESS).
- [IO] **Restart:**
  - Run the initiation with the ANALYS option **"Saving Restart Files"** checked.
  - Restart database files: **`COOXxxx`** and **`COOTKxxx`**, where `xxx` = 3-digit
    **frequency order number**. Example: `COOX001`, `COOTK001` are needed to restart
    the 1st frequency.
  - Index files **`COOXI`**, **`COOTKI`** must also be in the restart working
    directory, along with **`DOFSMAP`**, **`FILE90`**, **`FILE91`**.
  - Because of the naming, restart works only with the **full sequential set** of SSI
    frequencies. To restart on a subset, rename the COOX/COOTK files to match.
  - ✱ Ch. 6: "New Structure" restart needs only the COOXxxx files.
  - *(impl. note, inferred: COOX = impedance per frequency; COOTK = factored total
    system stiffness per frequency.)*
- [ADV/IO] **Global "unconstrained" soil impedance option** writes **`FILE11`**,
  which is very large.
  > **WARNING:** Use global soil impedances only for **surface stick models with rigid
  > basemats**. Accuracy is crude for elastic foundations and interpretation is hard.
  > For embedded models the unconstrained global impedances are correct **only with
  > FI-FSIN**.
- [CORE] **Simultaneous cases** in one run, with no restart needed:
  - coherent seismic: all **three X, Y, Z** inputs at once
    (✱ writes `FILE8X`, `FILE8Y`, `FILE8Z`);
  - incoherent seismic: **up to 50 simulations** per run
    (✱ writes `FILE8xxx`, up to 150 = 50 × 3 directions; FILE8001/002/003 = sim 1
    X/Y/Z);
  - vibration/external forces: **up to 500 load cases** (✱ writes `FILE8xxx`,
    xxx ≤ 500).
  - The limit depends on model size and RAM. If RAM runs out, ANALYS fails with
    **"access violation"**. *(impl. note: the Python version should estimate memory
    first and raise a clear error.)*

### 3.8 Module MOTION — interpolation and nodal responses

- [IO] Input deck **`.mot`** (AFWRITE). Requires **only `FILE8`** as SSI input.
- [CORE] Reads the TFs from FILE8, **interpolates them in frequency with a complex
  domain scheme**, and uses the interpolated TFs to compute SSI response motions at
  the user-selected nodes. *(impl. note: the details of the interpolation and the
  FFT convolution with the control motion are in the MOTION input chapter.)*
- [CORE] Accelerations, velocities, displacements and response spectra can be
  requested at chosen nodes and DOFs (with the baseline-correction option for
  integrated quantities). Relative displacements obtained this way are **much more
  approximate** than those from **RELDISP**.
- [IO] **`FILE13`** (text): nodal motions, **4 columns** per time step: accumulated
  time, absolute acceleration, absolute velocity, absolute displacement.
- [IO] **`FILE12`** (optional): computed acceleration histories and response spectra.
- [IO] The output listing can be very large if time histories are printed.
- [IO] Post-processing text files (three translational DOFs):
  - **`.TFU`**: computed (uninterpolated) TF;
  - **`.TFI`**: interpolated TF;
  - **`.ACC`**: acceleration time history;
  - **`.RS`**: in-structure response spectra for the selected damping values.
  - Naming: **`xxxxxTR_y.ext`**, where `xxxxx` = node number (5 digits for models
    with < 100,000 nodes, 6 digits otherwise), `y` ∈ {X, Y, Z}, `ext` ∈ {TFU, TFI,
    ACC}.
  - RS naming: **`xxxxxTR_yzz.RS`**, where `zz` = 2-digit damping-ratio order number
    (e.g. `01` → 0.02 and `02` → 0.05 if those two ratios were chosen).
- [UI/IO] **MOTION post-processing restart option** writes frame files to
  subdirectories **`\TFU`**, **`\RS`**, **`\ACC`** (plus **`\ACCR`** for rotational
  accelerations). Each frame holds the values for **all active nodal DOFs** at one
  frequency step or time step. The UI/PREP uses them for bubble plots, TF vector
  plots, contour plots and deformed-shape animation (Table 3.2). MOTION also writes
  **`ACC_max.txt`**, the maximum-acceleration (ZPA) frame.
- [IO] Compressed binary database of acceleration histories:
  **`Modelname_ACC.bin`** (see **BINOUT** command). **COMBACCDB** combines the X, Y,
  Z input-direction databases.
- [IO] **External RS mode:** MOTION can compute response spectra for external
  acceleration histories listed in the **`CONTTRS`** file. These files must use the
  `.ACC` format and the **`.ACC` extension**. Results go to files with the same name
  and the **`.RSO`** extension.

### 3.9 Module STRESS — element stresses, strains, forces

- [IO] Input deck **`.str`** (AFWRITE). Requires **`FILE4`** and **`FILE8`**.
- [CORE] Computes stress, strain and force **time histories and peak values** in the
  structural elements.
- [ADV] Nonlinear SSI: writes **`FILE74`** after each SSI iteration (effective shear
  modulus and damping of near-field soil elements, read by HOUSE next iteration) and
  reads **`FILE78`** from HOUSE. **COMB_XYZ_STRAIN** combines the X, Y, Z
  directional shear strains by **SRSS** into a new FILE74.
- [IO] Text files for post-processing:
  - **`.TFU`**: element stress TF, computed;
  - **`.TFI`**: element stress TF, interpolated;
  - **`.THS`**: stress time histories.
  - Naming **`etype_gnum_enum_comp`** + extension. Example
    `BEAMS_003_00045_MXJ` = moment MX at end J of BEAM element 45 in group 3.
    `gnum` is 3 digits and `enum` is 5 digits.
- [IO] Binary database of element output histories **`Modelname_STRESS.bin`**
  (BINOUT). **COMBTHSDB** combines the X, Y, Z input-direction databases
  (Demo 9).
- [IO] **`ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`**: maximum absolute stress components
  at element **centres**, for all elements.
  > **WARNING:** It is written only if **all element stresses in all elements** are
  > selected with the post-processing check box.
  >
  > **WARNING:** Not usable for **beam** elements, which need end forces and moments
  > at the nodes. It is meaningful for **shell, solid and spring** elements.

  Format for a model with n groups:
  ```
  [# of groups]
  [group-1 element type] [group #] [ordered group #] [# of elements in group]
  ... (one line per group, n lines)
  [group-1 element type] [group #] [ordered group #]
  [element #] [comp1] [comp2] [comp3] [comp4] [comp5] [comp6]   (one line per element)
  [group-2 element type] [group #] [ordered group #]
  ...
  ```
  Example (from the manual):
  ```
  3
  SOLID           1    1       3
  SHELL           2    1       1
  SOLID           3    1       5
  SOLID           1    1
                1 0.69504 0.61290 0.93326 0.21454 1.36011 0.45008
                2 0.82394 0.70086 0.68225 0.20217 0.65360 0.34301
                3 1.63296 1.09535 1.41437 0.49395 1.40079 0.39915
                SHELL    2    1
                1 6.98477 12.51727 9.28106 0.87311 0.58065 0.44423
                SOLID    3    1
                1 1.07909 1.12969 1.91468 0.10359 1.03872 0.25992
                ... (elements 2–5)
  ```
  Notes: "ordered group #" is 1 for every group in the example, so its meaning is
  unclear (see open questions). Element numbers restart at 1 within each group. Six
  components are always listed.
- [UI/IO] **STRESS post-processing restart option** writes frame files to
  **`\NSTRESS`** for node-stress contour plots (static at a chosen time, maximum, or
  animated).
  - Only **SOLID and SHELL** elements of 3D models are handled.
  - With both SOLID and SHELL present, frames hold only **average nodal membrane
    stresses**.
  - SHELL-only models also get separate **average nodal bending-stress** frames
    (`bd` in the name).
  - Naming is in Table 3.2.
- [UI/IO] **Soil-pressure frames** apply when near-field soil elements are adjacent
  to foundation walls. They are saved in **`\SOILPRES`**: seismic pressure frames
  at each time step plus **one frame of maximum pressures**.
  - **`STATIC_SOIL_PRESSURES.TXT`** is generated when soil-pressure frames are
    requested. It has the same format as `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT` but
    **one value column** (element soil pressure).
  - On the first run it holds **zeros**. If the user enters non-zero static (bearing,
    geological) pressures and runs the STRESS soil-pressure restart again, total
    pressure = static + seismic (**algebraic sum**) is written to `\SOILPRES`.
  - Two more files: **`pres_max_ele`** (maximum element soil pressures) and
    **`pres_max_nod`** (average nodal soil pressures, for plotting only), for the
    SOLID elements modelling adjacent near-field soil.
  > **WARNING:** No text frame files for **TSHELL** elements.
  >
  > **WARNING:** Frame nodal stresses and pressures are **averages** of element-centre
  > values, with no shape-function extrapolation. Nodal stress is taken as equal to
  > the element-centre stress, and soil pressure is the normal stress on the SOLID
  > face. Use them for **plotting only**. Design values are the **element-centre**
  > values (STRESS output, `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT`, `pres_max_ele`).
- [IO] **Element output frames** at all time steps, **`ESTRESS_framenumber.ess`**,
  enabled in **UI Options/Analysis/STRESS**. Section-cut calculations need them
  (Demo 8, Verification Problem 47).
  > **WARNING:** No ESTRESS frames for **TSHELL** elements.

### 3.10 Module RELDISP — relative displacements

- [IO] Input deck **`.rdi`** (AFWRITE).
- [CORE] Uses the **interpolated complex ATF from MOTION** (`.TFI` files) to compute
  **relative displacements** between selected nodes. **MOTION must run first** for
  every node of interest.
- [IO] Output listing of **maximum nodal relative displacements**; **`.TFD`** (complex
  relative-displacement TF) and **`.THD`** (relative-displacement time histories),
  named as in Table 3.1.
- [UI/IO] **RELDISP post-processing restart option** writes frames to **`\THD`** (and
  **`\THDR`** for rotational displacements), with all active nodal DOFs per time step,
  for deformed-shape animation.
- [IO] Binary databases (BINOUT). The manual calls it "Modulename_STRESS.bin" here,
  which is apparently a copy-paste error (see open questions). **COMBTHDDB** combines
  X/Y/Z input-direction component histories.
  - For relative-displacement X/Y/Z combination, run RELDISP **three times per input
    direction** (principal plus coupled responses), giving for the X input:
    **`Modelname_TR_X_THD.bin`**, **`TR_Y_THD.bin`**, **`TR_Z_THD.bin`**.
  - **COMBDISPDIR** merges these three into **`Modelname_THD.bin`**. Run it once for
    each input direction.

### 3.11 Module COMBIN — merge frequency sets

- [CORE/IO] Combines results for **different frequency sets from two ANALYS runs**,
  e.g. when extra frequencies turn out to be needed.
- Input: two FILE8-type solutions renamed **`FILE81`** and **`FILE82`**. Output: a new
  **`FILE8`**. *(impl. note: merge by frequency, sort ascending, reject or resolve
  duplicates, and require identical DOF maps.)*

### 3.12 Module NONLINEAR [ADV, Option NON]

Details are in manual §1.5.4 and Demos 9–10. Steps:
1. Run an initial **linear SSI** analysis with the initial elastic properties of the
   nonlinear elements.
2. Compute the **local time-domain behaviour** of each nonlinear element from the local
   relative displacements, and calibrate an equivalent **linearized hysteretic model**
   (complex frequency) for each element.
3. Run a new SSI iteration with a **fast SSI restart** in the complex frequency domain,
   using the step-2 linearized models.
4. Check convergence, then stop or continue.

**COMB_XYZ_THD** (companion) combines **co-directional** displacement histories for the
X, Y and Z inputs (§1.5.4, Demos 9 and 10).

---

## Part 3 — Tables 3.1 and 3.2 (full transcription) [IO]

### Table 3.1 — Useful text files for result verification and post-processing

Caption note: node numbers need **6 digits** for models with **more than 99,999
nodes**. All character positions below assume 5-digit node numbers.

| Extension / file | Description |
|---|---|
| **RS** | Response spectra data files generated by the MOTION module |
| | *Naming scheme for TFU, TFI, TFD, ACC files* (as printed under RS; these positions apply to `.RS` names): Characters **1–5** = node number; Characters **6–9** = translation (`TR`) or rotational (`R`) DOF (e.g. `TR_X`); Characters **10–11** = damping ratio number |
| **TFU** | Uninterpolated acceleration transfer functions written by MOTION, and stress transfer functions (by STRESS) |
| **TFI** | Interpolated acceleration transfer functions written by MOTION, and stress transfer functions written by STRESS |
| **TFD** | Displacement transfer functions generated by RELDISP |
| **THD** | Displacement time history written by RELDISP |
| **ACC** | Acceleration time history written by MOTION |
| | *Naming scheme for acceleration TFU, acceleration TFI, TFD, THD and ACC files*: Characters **1–5** = node number; Characters **6–9** = translation (`TR`) or rotational (`R`) DOF |
| **TH** | Soil time history for layers (SOIL) |
| | *Naming scheme*: `ACC***` = acceleration time history for soil layer *** (e.g. `ACC001.TH` → layer 1); `SN***` = strain time history for soil layer *** (manual's example: "`SN001.TH` is the strain time history for soil layer 2", which looks like a typo); `SS***` = stress time history for soil layer *** (manual's example: "`SS001.TH` … soil layer 3", same typo) |
| **THS** | Stress time history written by STRESS |
| | *Naming scheme for THS, stress TFU and stress TFI*: `etype_gnum_enum_comp`, e.g. `BEAMS_012_00001_FXI.THS`. `etype` = element type; `gnum` = group number; `enum` = element number; `comp` = stress component |
| **Frames.txt** | Post-processing frames for stress and motion |
| **ELEMENT_CENTER_ABS_MAX_STRESSES.TXT** | List of maximum stresses for each element |
| **STATIC_SOIL_PRESSURES.TXT** | Defines extra soil pressure (geological pressure) to include in soil-pressure frames |
| **SRSSTF.txt** | SRSS option in MOTION |

Concrete name patterns *(impl. note, from the examples)*:
- nodal TF and history: `f"{node:05d}TR_{dof}.{ext}"`, with dof ∈ {X, Y, Z} and ext ∈
  {TFU, TFI, ACC, TFD, THD}. Use `:06d` when the model has more than 99,999 nodes.
- nodal RS: `f"{node:05d}TR_{dof}{idamp:02d}.RS"`.
- element: `f"{etype}_{group:03d}_{elem:05d}_{comp}.{ext}"`, with ext ∈ {THS, TFU,
  TFI}; e.g. `BEAMS_012_00001_FXI.THS`.
- soil layers: `f"ACC{layer:03d}.TH"`, `f"SN{layer:03d}.TH"`, `f"SS{layer:03d}.TH"`.

### Table 3.2 — Post-processing frame files produced by MOTION, RELDISP and STRESS

| Frame type | Pattern | Example (as printed) | Fields |
|---|---|---|---|
| RS frames | `RS##_freq_filenum` | `\RS\RS01_000.10_00001` | `##` = damping number; `freq` = frequency; `fnum` = frame number |
| TFU frames | `TFU_freq_filenum` | `\TFU\TFU_000.02_00001` | `freq` = frequency; `fnum` = frame number |
| ACC frames | `ACC_time_filenum` | `\ACC\ACC_00.000_00001` | `time` = time; `fnum` = frame number |
| THD frames | `THD_time_filenum` | `\THD\THD_00.000_00001` | `time` = time; `fnum` = frame number |
| Stress frames | `stress_time_fnum_comp` | `\NTRESS\stress_00.000_00001_sig` (typo for `\NSTRESS`) | `time`; `fnum`; `comp` = stress component (below) |
| Soil pressure frames | `press_time_fnum_type` | `\SOILPRES\pres_00.000_00001_nod` (prefix printed as `pres`) | `time`; `fnum`; `type` = `ele` (element values) or `nod` (nodal values) |

Stress component codes (`comp`):

| Code | Solids | Shells |
|---|---|---|
| `sig` | Normal stress | Membrane stress |
| `tau` | Shear stress | Membrane shear |
| `bdsig` | — | Bending stress (shell elements only) |
| `bdtau` | — | Bending shear (shell elements only) |

**Maximum-value frames:**

| Frame | Pattern | Example | Fields |
|---|---|---|---|
| Stress | `stress_ABS_MAX_comp` | `\NSTRESS\stress_ABS_MAX_sig` | `comp` ∈ {sig, tau, bdsig, bdtau}, as above |
| Soil pressure | `press_ABS_MAX_type` | `\SOILPRES\pres_ABS_MAX_nod` | `type` ∈ {ele, nod} |
| Acceleration (text in §3.8) | `ACC_max.txt` | — | ZPA frame |

Numeric field formats in the examples *(impl. note)*:
- frequency `000.10` → `f"{freq:06.2f}"`;
- time `00.000` → `f"{t:06.3f}"`;
- frame number `00001` → `f"{n:05d}"`.

No file extension is shown in the frame names.

---

## Part 4 — Section 3.2: Initiation and restart analyses

### 4.1 Working directory and run sequences [CORE/IO]

- Keep **all text input files written by AFWRITE in the same working directory**.
- **Linear seismic SSI:** EQUAKE (if needed) → SOIL (if needed) → **SITE → POINT →
  HOUSE → ANALYS** → MOTION (if needed) → STRESS (if needed). ✱ §1.3 adds RELDISP
  after MOTION.
- **Linear external-force (vibration) SSI:** SITE → POINT → **FORCE** → HOUSE →
  ANALYS → MOTION (if needed) → STRESS (if needed).
- ✱ Nonlinear soil (from §1.3): iterate HOUSE → ANALYS (restart) → STRESS
  (+ COMB_XYZ_STRAIN).
- ✱ Nonlinear structure: SITE-POINT-HOUSE-ANALYS-MOTION-RELDISP-COMB_XYZ_THD-NONLINEAR
  first, then HOUSE-ANALYS(restart)-MOTION-RELDISP-COMB_XYZ_THD-NONLINEAR.

**Dependencies to enforce** (for a pipeline runner):

| Module | Needs before it runs |
|---|---|
| POINT | SITE (FILE2) |
| HOUSE | `.sit` present (always). Embedded models: SITE/POINT. Nonlinear iteration: FILE74 |
| ANALYS | FILE1, FILE3, FILE4 (+FILE9 / FILE77) |
| MOTION | FILE8 |
| RELDISP | MOTION `.TFI` for its nodes |
| STRESS | FILE4, FILE8 |
| COMBIN | FILE81, FILE82 |

### 4.2 Nonlinear-analysis control files [ADV/IO]

- Nonlinear soil SSI needs two extra input text files: **`.liq`** and **`.pin`** (see
  HOUSE input for the nonlinear option).
- **`.liq`** is written by HOUSE in the initial run. If it is **non-empty and contains
  `1`**, HOUSE takes the model materials from **FILE74** in the next iteration.
- NONLINEAR uses the analogous files **`PANEL.NON`** (nonlinear panels) or
  **`SPRING.NON`** (nonlinear springs). If the `.NON` file **exists and contains
  `1`**, NONLINEAR uses the previous iteration's SSI results and runs a nonlinear
  analysis for the new iteration. If it **does not exist**, NONLINEAR runs a
  **linear** analysis (1st iteration).

### 4.3 Interactive and batch execution [UI/IO]

- **Interactive:** run one module at a time from the UI **"Modules"** menu, or run a set
  of modules automatically with **Options\Write "Running SSI Modules"**.
- **Batch (DOS window):** one command per module:
  ```
  SSI_module_name.exe < SSI_module_name.inp
  ```
  e.g. `SITE.exe < SITE.inp`. Executables are installed by default in
  **`C:\ACSV300\EXEB`** (also on the installation DVD).
- Each **`.inp`** file has exactly **three lines**:
  ```
  modelname
  modelname.ext_input              (ext_input = extension written by AFWRITE, e.g. .sit)
  modelname_SSI_module_name.out    (output listing)
  ```
- NQA version: the V&V runs are set up for batch mode and serve as examples for
  building batch files.
- *(impl. note: the Python CLI should accept the same 3-line `.inp`, e.g.
  `python -m sassi.site < SITE.inp`, so that batch files can be reproduced.)*

### 4.4 Restart and reanalysis types [CORE]

| Change | Name in manual | Modules to re-run | ✱ ANALYS mode (Ch. 4/6) |
|---|---|---|---|
| a. New control motion time history or RS, **same wave field type** | "New Time History" (MOTION restart) | **MOTION only** | — (no ANALYS) |
| b. New seismic environment (e.g. SV/SH → vertical P or Rayleigh, same control-point motion) | "New Seismic Environment" (ANALYS restart) | **SITE + ANALYS**. With incoherency (it changes the seismic loads): **HOUSE + ANALYS** | Mode 3 (needs COOXxxx + COOTKxxx) |
| c. New dynamic loads applied on the structure | "New Dynamic Loading" (ANALYS restart) | **FORCE + ANALYS + MOTION**. If only the load **time history** changes and the spatial pattern does not: **MOTION only** | Mode 3 |
| d. Structure or near-field soil changed, **interaction nodes unchanged** | "New Structure" (ANALYS restart) | **HOUSE + ANALYS + MOTION + STRESS** | Mode 2 (needs COOXxxx) |

For nonlinear-soil SSI iterations (type d), STRESS computes the effective near-field
soil shear modulus and damping and passes them to HOUSE through **FILE74**.
✱ Initiation = ANALYS Mode 1. Mode 6 ("New Load Vector") is batch-only and rarely used.

---

## Part 5 — Master file inventory (this section) [IO]

| File | Kind | Producer | Consumer(s) | Content / notes |
|---|---|---|---|---|
| `.equ` `.soi` `.sit` `.poi` `.hou` `.frc` `.anl` `.mot` `.str` `.rdi` | text input decks | UI AFWRITE | EQUAKE, SOIL, SITE, POINT, HOUSE, FORCE, ANALYS, MOTION, STRESS, RELDISP | Module inputs (formats in later specs) |
| `.rsi` | text | user | EQUAKE | Target RS: freq, amplitude (one damping) |
| `.acc` `.vel` `.dis` | text | EQUAKE | SOIL, MOTION, STRESS, RELDISP | Generated histories |
| `.rso` `.psd` `.fft` | text | EQUAKE (`.RSO` also MOTION) | user | RS; PSD (2 col); FFT (3 col, f ≥ 0) |
| `ACCxxx.TH` `SNxxx.TH` `SSxxx.TH` | text | SOIL | user | Layer acc / strain / stress histories |
| FILE73 | text | SOIL | STRESS | Soil material curves for nonlinear SSI |
| FILE88 | text | SOIL | SITE | Iterated equivalent-linear soil properties |
| FILE1 (FILE1X/Y/Z) | binary | SITE | ANALYS | Free-field motion data (U′) |
| FILE2 | binary | SITE | POINT (and SITE Mode 2) | Transmitting-boundary eigen data |
| FILE3 | binary | POINT2/3 | ANALYS | Point-load / flexibility data |
| FILE4 = `modelname.n4` | binary | HOUSE | ANALYS, STRESS | Structure + excavated-soil K, M |
| COOSK, COOSM | binary | HOUSE | ANALYS | Stiffness and mass matrices (fast solver) |
| `.hounew`, `.map` | text | HOUSE (optimizer) | user, HOUSE | Renumbered model; old ↔ new node map |
| FILE77 / FILE77xxx | binary | HOUSE | ANALYS | Incoherency / multi-support decomposition, X+Y+Z |
| FILE78 | — | HOUSE | STRESS | Non-empty only for nonlinear SSI |
| FILE74 | text ✱ | STRESS (COMB_XYZ_STRAIN) | HOUSE | Iterated near-field soil G, D |
| `.liq`, `.pin` | text | HOUSE / user | HOUSE | Nonlinear soil control |
| PANEL.NON, SPRING.NON | text | — | NONLINEAR | Iteration flag ("1") |
| FILE9 / FILE9001…FILE9500 | binary | FORCE (user renames) | ANALYS | Load vectors |
| FILE8 / FILE8X/Y/Z / FILE8xxx | binary | ANALYS, COMBIN | MOTION, RELDISP, STRESS, COMBIN | Complex TF solution, all DOFs, all frequencies |
| FILE81, FILE82 | binary | user rename of FILE8 | COMBIN | Two solutions to merge |
| COOXxxx, COOTKxxx, COOXI, COOTKI | binary | ANALYS | ANALYS restart | Per-frequency restart database + index |
| DOFSMAP, FILE90, FILE91 | binary | HOUSE ✱ | ANALYS restart | DOF map and auxiliary data |
| FILE11 | binary | ANALYS | user | Global unconstrained impedances (large) |
| FILE12 | — | MOTION | user | Optional acc histories + RS |
| FILE13 | text | MOTION | user | time, abs acc, abs vel, abs disp |
| FILE14, FILE15 ✱ | text | STRESS | user | Beam force TFs; stress histories |
| `*.TFU` `*.TFI` `*.ACC` `*.RS` | text | MOTION | user, RELDISP (TFI) | Nodal results |
| `*.TFD` `*.THD` | text | RELDISP | user | Relative-displacement TF and history |
| `*.THS` (+ stress `.TFU`/`.TFI`) | text | STRESS | user | Element results |
| `\TFU \RS \ACC \ACCR \THD \THDR \NSTRESS \SOILPRES` | text dirs | MOTION, RELDISP, STRESS restarts | UI | Frames (Table 3.2) |
| `ACC_max.txt` | text | MOTION | UI | ZPA frame |
| `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT` | text | STRESS | user | Max element-centre stresses |
| `STATIC_SOIL_PRESSURES.TXT` | text | STRESS / user edits | STRESS | Static pressures to add |
| `pres_max_ele`, `pres_max_nod` | text | STRESS | user, UI | Max soil pressures |
| `ESTRESS_framenumber.ess` | text | STRESS | section-cut tools | Element output frames |
| `Modelname_ACC.bin`, `Modelname_STRESS.bin`, `Modelname_TR_y_THD.bin`, `Modelname_THD.bin` | binary | MOTION / STRESS / RELDISP (BINOUT) | COMBACCDB, COMBTHSDB, COMBTHDDB, COMBDISPDIR | Compressed history databases |
| `CONTTRS` | text | user | MOTION | List of external `.ACC` files for RS |
| `Frames.txt`, `SRSSTF.txt` | text | MOTION/STRESS | — | Frame list; SRSS option |
| `SSI_module_name.inp` | text | user | batch | 3-line batch input |

---

## Part 6 — Verification hooks derived from this section *(impl. note)*

1. **Zero-SSI identity (FV):** fill the excavation with "structure" solids that have
   exactly the free-field layer properties, so C^s = C^e on all excavated nodes, with
   no superstructure. Using the FV method, the ATF at every node must equal the
   free-field TF from SITE (U = U′), to solver precision. This is the manual's §2.1
   statement and checks the assembly signs, the impedance inversion and the
   load-vector formation.
2. **Surface rigid massless mat under vertical SV:** with no excavation (Eq. 2.1 reduces
   to the s/i blocks), the mat ATF must be about 1 at all frequencies. Kinematic
   interaction vanishes for a surface foundation under vertical waves.
3. **FV vs FFV vs FI-EVBN vs FI-FSIN** on one embedded model: FV is the reference.
   Expected accuracy order is FFV > EVBN > FSIN, and FSIN may show spurious peaks.
4. **Impedance symmetry:** X_ff must be complex symmetric (X = Xᵀ). Check
   ‖X − Xᵀ‖/‖X‖ < 1e-10.
5. **EQUAKE SRP checks:** unit-test each criterion (−10 % / +30 % band, ≤ 9 adjacent
   low points, ≥ 100 pts/decade, duration ≥ 20 s, Δt ≤ 0.005 s warning).
6. **COMBIN:** merging FILE81 (odd frequencies) and FILE82 (even frequencies) must
   reproduce the FILE8 of a single run with all frequencies.
7. **Restart equivalence:** an ANALYS Mode 3 run with a changed FILE1 must equal a
   fresh initiation run (Mode 1) with that FILE1. A Mode 2 run with a changed FILE4
   (same interaction nodes) must equal a fresh run.

---

## Open questions / ambiguities

1. **Complex-modulus damping form.** The manual says only that damping uses complex
   moduli. Options are K(1 + 2iβ) or the SASSI/SHAKE form
   K(1 − 2β² + 2iβ√(1 − β²)). Check the HOUSE/SITE chapters, or choose one and
   document it (recommend the SASSI form as default, with a switch).
2. **`w` nodes with structural stiffness.** Eq. (2.1) has no C^s on `w` and no `s`–`w`
   coupling. Real models can have internal slabs or walls through excavated-soil
   nodes. Recommendation: use general global assembly (§1.4.2), which reduces to
   Eq. (2.1).
3. **Definition of FSIN vs EVBN node sets** for automatic selection: are the
   ground-surface nodes of the excavation top face (EVBN) all nodes at grade inside the
   footprint? Are the FFV internal layers chosen by the user only (by elevation or
   node list)? The UI selection mechanism is not in this section.
4. **Interaction nodes and layer interfaces.** The text implies that free-field
   motion is known only at layer interfaces. The tolerance for node elevations versus
   layer interfaces, and the error behaviour, are not specified.
5. **ATF normalisation.** Is the FILE1 free-field vector normalised to unit control
   *displacement* or unit control *acceleration*? Both give the same ratio, but the
   storage convention in FILE8 matters for STRESS and RELDISP (relative
   displacements need TF/(−ω²) scaling).
6. **COOX / COOTK contents.** The manual gives only the names. Inferred:
   COOX = impedance and COOTK = factored total stiffness. The index files COOXI/COOTKI
   and DOFSMAP/FILE90/FILE91 have unspecified content. Python can use its own
   internal format but should keep the names.
7. **Binary file formats** (FILE1, FILE2, FILE3, FILE4, FILE8, FILE9, FILE77, FILE11,
   COOSK/COOSM) are proprietary Fortran unformatted files. Choose a documented
   container (e.g. NPZ/HDF5) under the same names. Byte-compatibility is not
   possible from this manual.
8. **EQUAKE Nyquist wording:** "The Nyquist frequency is not higher than 100 Hz" alongside
   "warning if Δt > 0.005 s" is contradictory. Δt = 0.005 s gives Nyquist = 100 Hz,
   so the intended rule is presumably Nyquist ≥ 100 Hz (Δt ≤ 0.005 s). Implement the
   Δt check as stated.
9. **EQUAKE PSD normalisation** (one-sided vs two-sided, per Hz vs per rad/s, the
   2π T_D factor) and the exact meaning of "±20 % averaging". The units stated are
   in²/s³ (for ft/s² gravity units, so the length unit is **inches**, not ft) and
   cm²/s³ (for m/s²). Check the EQUAKE input chapter. SRP 3.7.1 Appendix A is the
   likely reference.
10. **EQUAKE cross-correlation threshold** (SRP: |ρ| ≤ 0.16) is not stated here.
11. **LW and AB algorithm parameters** (iterations, tolerances, wavelet type for AB) are
    not given here. They are presumably in the EQUAKE input chapter.
12. **FILE9 multi-case suffix:** the text says "Only one digit load case number can be
    appended" but the examples use 3 digits (`FILE9001`…`FILE9500`). Use 3-digit
    zero-padded suffixes.
13. **Table 3.1 TH examples** say `SN001.TH` = layer 2 and `SS001.TH` = layer 3,
    against the stated `***` = layer number rule. Treat as a typo: `SNxxx` and `SSxxx`
    = layer xxx.
14. **Table 3.1 character positions 6–9** ("TR or R") versus the name `xxxxxTR_y`.
    For rotations, is it `xxxxxR_y` (3 chars) or `xxxxxRR_y`? Positions shift by
    one with 6-digit nodes. The table lists "Characters 10–11 = damping ratio
    number" under the TFU/TFI/TFD/ACC heading, but it only applies to `.RS`.
15. **Table 3.2 inconsistencies:** pattern `press_…` versus example `pres_…`;
    `\NTRESS` versus `\NSTRESS`; the term "filenum" versus "fnum". Frame file
    extensions are not shown (and `ACC_max.txt` has `.txt`). Field widths are inferred
    from the examples only.
16. **RELDISP binary database name** is printed as "Modulename_STRESS.bin", likely a
    copy-paste error. Use `Modelname_THD.bin` (as COMBDISPDIR output) and
    `Modelname_TR_{X|Y|Z}_THD.bin`. Also "Modulename" versus "Modelname".
17. **"ordered group #"** in ELEMENT_CENTER_ABS_MAX_STRESSES.TXT is 1 for every group in
    the example. Its meaning (order within element type? sequence?) is unclear. The
    meaning and order of the six stress components per element type (SOLID σx, σy, σz,
    τxy, τyz, τzx? SHELL?) are not given here.
18. **MOTION sentence on baseline correction / FILE13** is garbled. It is unclear
    whether FILE13 is written only with baseline correction on, and which DOFs or
    nodes it covers.
19. **SITE outputs sentence is truncated** ("In addition to the output and binary files
    FILE1 and FILE2."). There may be more outputs.
20. **HOUSE `.sit` requirement:** the text says HOUSE runs independently of SITE/POINT
    except for embedded models, and then that it "requires" `.sit` in the working
    directory. Does a surface model also need `.sit`? Recommend: always require it
    (it supplies the soil layering for excavated-soil properties and units).
21. **COMB_XYZ_THD naming variants** (COMBIN_XYZ_THD, COMBINE_XYZ_THD). Choose
    `COMB_XYZ_THD` and accept aliases.
22. **Restart with frequency subsets:** the manual says the COOX/COOTK files would need
    to be renamed. A Python version could store the frequency value with each
    record, which removes this limitation while keeping the names.
23. **MOTION "complex domain" interpolation scheme** is not defined here. The SASSI
    classic approach (Tajirian 1981) fits a rational, 2-DOF-like transfer function
    between solved frequencies. Confirm in the MOTION chapter.
24. **Mass matrix assumptions** (✱ §1.5.1: 50 % lumped + 50 % consistent, except beams
    consistent and plates lumped) affect C = K − ω²M. They belong to the HOUSE spec but
    must match the theory here.
