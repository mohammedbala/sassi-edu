# 03 — Application Guidelines (ACS SASSI Manual, Chapter 4)

**Source range:** `reference/acs-sassi.txt` lines 2433–3384, which is manual Chapter 4 "Application Guidelines". This is printed pages 55–73, or PDF pages 57–75. All of those PDF pages were also viewed as images.
**Content type:** the chapter is about procedure and methodology. It has no dialog boxes or tables, and only three numbered equations (4.1–4.3). Those equations, plus the inline formulas for the half-space depth and the rigid-basement factor, were checked against the PDF images and are transcribed below.
**Cross-references:** where Chapter 4 points to other chapters (1.5, 2.1, 3, 6, 9), a short `[xref]` note is added. The note gives enough context for the reader, but the full spec for that material belongs to the spec file for that chapter.

Tag legend:

| Tag | Meaning |
|---|---|
| **[CORE]** | Needed for computational correctness (algorithms, formulas, numerical rules) |
| **[IO]** | File / input-format / file-naming behaviour |
| **[UI]** | User interface, plotting, convenience, checks presented to the user |
| **[ADV]** | Advanced options (incoherency, nonlinear, Option A / AA / NON / PRO, multiple excitation) |
| **[DECISION]** | Not specified by the manual; a proposed implementer choice (see §9 Open questions) |

---

## 1. Chapter scope

Chapter 4 has two parts:

* **§4.1 SSI Analysis Procedure**
  * §4.1.1 is a 13-step recipe for a *linearized* (equivalent-linear) seismic SSI analysis, with changes for forced-vibration problems.
  * §4.1.2 gives 20 "Engineering Considerations": modelling rules, mesh and frequency rules, damping/half-space treatment, incoherency guidance, and so on. These add to the capabilities and limitations in Manual §1.5.
* **§4.2 ACS SASSI Runs** says which modules to run, and in what order, for each of the three *run types*:
  1. Initiation SSI solution runs
  2. Post-processing runs
  3. Restart (re-analysis) solution runs. These are the ANALYS restart modes "New Structure", "New Seismic Environment" and "New Dynamic Loading".

For our Python re-implementation, this chapter defines four things:
1. **The run orchestration / pipeline state machine** (§2, §4 below).
2. **The model-checking and parameter-validation rules** that the UI should enforce or warn about (§3, collected in §6).
3. **The core numerical rules** (§5): the frequency grid (Eq. 4.1/4.2), the 1/5-wavelength mesh rule, the half-space buffer depth, and impedance as the inverse of compliance (Eq. 4.3).
4. **The file naming conventions** for multi-direction, multi-simulation and multi-load-case runs (§8).

---

## 2. §4.1.1 Computational steps for a linearized SSI analysis

### 2.1 Step table (seismic analysis)

| Step | Action | Module(s) | Inputs | Outputs | Tag |
|---|---|---|---|---|---|
| 1 | Generate a spectrum-compatible acceleration time-history set for the **control motion**. It must meet ASCE 4 and USNRC SRP 3.7.1. | **EQUAKE** | Target response spectra | Acceleration histories, used by SOIL, MOTION, RELDISP, STRESS | [CORE] |
| 2 | Find the important frequency ranges from the dominant frequencies of the structure. Use a **fixed-base modal analysis** in a standard FE code (for example ANSYS). The UI has automatic **ANSYS ⇄ ACS SASSI** model converters. After the first SSI run, add frequencies as needed to capture all the structural modes. Solve the added frequencies separately, then merge the result databases with the **COMBIN** module. (The manual calls it "COMBINED" in this step. That is a typo for COMBIN.) | external FE / UI converters, COMBIN | FE model | Frequency ranges of interest | [CORE]/[UI] |
| 3 | Set the **SSI cut-off frequency**. Base it on (a) the frequency content of the control motion and (b) the **90 % cumulative modal mass** criterion (ASCE 4, Section 3). Constraint: **f_cut ≤ maximum wave passing frequency** of both the soil-layer mesh and the FE model mesh (see §3 items 9, 12, 19). | — | Steps 1–2 | f_cut | [CORE] |
| 4 | Assume **vertically propagating waves**. Specify the **location of the control motion** and compute the free-field iterated **strain-compatible soil properties**. | **SOIL** (SHAKE-type equivalent linear) | Control motion, low-strain profile, degradation curves | Iterated Vs, damping per layer | [CORE] |
| 5 | Use the linearized (iterated) soil profile from SOIL for the site-response analysis in **SITE**. This can be automatic: choose the **"Nonlinear"** operation-mode option in the SITE input dialog. For embedded models, the UI also copies the iterated soil properties to the **HOUSE** model input for the **excavated soil layers**. | SITE, UI | SOIL results | SITE layer data, HOUSE excavated-soil materials | [CORE]/[UI] |
| 6 | Build the **structure FE model**, or several models for multiple structures. It may include **backfill soil** as a **near-field soil zone**. | UI / HOUSE input | — | Structure model | [CORE] |
| 7 | Build the **excavated soil FE model** for each embedded structure. For adjacent structures, the excavated soil volumes may be merged into one larger excavation volume. | UI / HOUSE input | — | Excavation model | [CORE] |
| 8 | Choose the substructuring method used for the impedance matrix and the SSI solution: **FV, FI or FFV**. | HOUSE option `<imp>` [xref §9.2.16: 0 = FV, 1 = FFV, 2 = FI] | — | Method flag | [CORE] |
| 9 | Choose the **SSI frequencies** at which the site-response and point-load problems are solved. Recommendation: **at least 40–80** for simple stick models and **100–300** for complex FE models. Base the choice on Steps 1–2, and add more later if needed to make the interpolated transfer functions more accurate. | FREQ set [xref §9.2.15] | Steps 1–2 | Frequency set (integer frequency numbers, Eq. 4.1) | [CORE] |
| 10 | **Initiation SSI analysis**: compute the complex acceleration transfer functions (ATF) at all nodes. Sub-steps a–d are in §2.2. | SITE, POINT, HOUSE, ANALYS | — | FILE8 | [CORE] |
| 11 | Compute the time-domain response from the transfer functions: (a) **COMBIN** if needed, to merge FILE8s from different frequency sub-ranges; (b) **MOTION, RELDISP, STRESS**. Check the computed vs interpolated **ATFs** (from MOTION) and **STFs** (from STRESS) to see whether more frequencies are needed. | COMBIN, MOTION, RELDISP, STRESS | FILE8 | Response histories, ARS, stresses | [CORE] |
| 12 | Add new SSI frequencies if the ATF/STF review in Step 11 shows they are needed. Loop back to ANALYS for the new frequencies, then COMBIN. | ANALYS, COMBIN | — | Augmented FILE8 | [CORE] |
| 13 | **Restart SSI analysis** when the superstructure changes ("**New Structure**" restart) or the seismic environment changes ("**New Seismic Environment**" restart). Forced-vibration substitutions are in §2.4. | HOUSE/SITE/FORCE + ANALYS | — | New FILE8 | [CORE] |

### 2.2 Step 10 detail: initiation-run module/file flow [CORE][IO]

**a. SITE: Mode 1, then Mode 2.** Uses the data from Steps 5 and 9 and the specified control-motion location.
  * *Mode 1* forms and solves the soil-layer eigenproblems that **POINT** needs to build the **transmitting boundary**. It writes **FILE2**.
  * *Mode 2* solves the site-response problem. It reads **FILE2** and the **`.sit`** input file, which describes the seismic wave field (wave types, control point). It writes the **free-field motions** to **FILE1**.

**b. POINT** (POINT2 for 2D / POINT3 for 3D) reads **FILE2** and the **`.poi`** input file. The `.poi` file gives the **maximum embedment of the structure** and the **radius of the point load** (the "Radius of Central Zone" [xref §6 POINT dialog]). POINT writes the **point-load solution** to **FILE3**.

**c. HOUSE(FS)**, the "fast-solver" variant of HOUSE, uses the data from Steps 6, 7 and 8. It writes:
  * **FILE4**: the complex stiffness and mass matrices of the structure and the excavated soil. FILE4 is **renamed `modelname.N4`**.
  * **COOSK** and **COOSM**, the FE model matrix files (stiffness and mass).
  * **Option AA [ADV]:** **HOUSE(FSA)** also reads **ANSYS FE model matrix files** (see the separate "ACS SASSI-ANSYS Integration Capability" manual). In that case the ACS SASSI model only supplies the **topology** (nodes and element connectivity). The topology is needed for **equation (DOF) mapping** and for UI post-processing.

**d. ANALYS(FS)**, the fast-solver ANALYS, reads:
  * **FILE1** (seismic) **or FILE9** (external-force option),
  * **FILE3**,
  * **FILE4** (`modelname.n4`).

  It produces the **soil impedance matrices**, the **triangularized** (factorized) **total-system stiffness**, and the **complex ATFs** for every requested SSI frequency and every nodal point. These are stored in **FILE8**.
  * Recommended practice for large embedded models: split the frequency range into **sub-ranges** and run ANALYS separately for each one. This cuts run time and scratch-file size. Merge the results afterwards with COMBIN (Step 11a).

Context [xref Ch. 2, Eq. 2.1–2.2]. ANALYS solves, for each frequency ω, the partitioned complex system built from the dynamic stiffness `C(ω) = K − ω² M` (K is complex, so it carries the damping), plus the free-field impedance `X_ff` added at the interaction DOFs. The load is `X_ff · u'_f`, where `u'_f` is the free-field motion at the interaction nodes taken from FILE1. Per Ch. 3, ANALYS does four things for each frequency: (1) forms the soil flexibility matrix, (2) computes the impedance at the interaction nodes, (3) forms the seismic load vector, including incoherency when it is active, and (4) solves for the complex response transfer functions at every DOF.

### 2.3 Steps 11–12: the frequency refinement loop [CORE][UI]

```
repeat
    ANALYS (new frequencies only; FILE1/FILE9, FILE3 must contain them; FILE4 is frequency-independent)
    COMBIN (FILE81 = old FILE8, FILE82 = new FILE8  ->  FILE8)
    MOTION / STRESS -> compare computed (.TFU) vs interpolated (.TFI) transfer functions
    identify interpolated peaks not supported by computed points (CRITFREQ command [xref §9.7.4])
until no significant unsupported interpolated peaks
```

### 2.4 Step 13: substitutions for **forced foundation vibration** problems [CORE]

When the excitation is an external dynamic force instead of an earthquake:

| Item | Change |
|---|---|
| a | The control-motion time history of Step 1 is replaced by the **reference time history of the external dynamic forces**. |
| b | **Step 4 (SOIL) is not used.** |
| c | The iterated soil properties of Step 5 are replaced by the **initial (low-strain) soil properties**. |
| d | The **site-response problem is dropped** from Step 9 and from the second part of Step 10a. That is, SITE runs **Mode 1 only**, and no FILE1 is made. |
| e | A **FORCE** module run produces the **load vector in FILE9**, and **FILE9 replaces FILE1** as ANALYS input. (The manual's sub-item letters here, "item d" and "item e of Step 10", do not match the a–d list above. See Open Questions.) |
| f | In Step 13, the "dynamic environment" is replaced by the **external dynamic forces**. The relevant restart is therefore **"New Dynamic Loading"**, not "New Seismic Environment". |

[xref Ch. 3.2] The module order for linear vibration analysis is **SITE → POINT → FORCE → HOUSE → ANALYS → MOTION (if needed) → STRESS (if needed)**.

### 2.5 Module / file dependency graph [IO][CORE]

```
EQUAKE ──(acc. histories)──► SOIL ──(iterated soil props)──► SITE ( & HOUSE excavated layers, via UI)
                                         │
                       SITE Mode 1 ──► FILE2 ──► SITE Mode 2 (+ .sit) ──► FILE1  (seismic only)
                                         └──────► POINT2/POINT3 (+ .poi) ──► FILE3
FORCE (+ .frc) ──► FILE9                              (vibration only)
HOUSE (+ .hou [, FILE74] [, ANSYS matrices: Opt AA]) ──► FILE4 = modelname.N4, COOSK, COOSM,
                                                          FILE77 / FILE77sss (incoherent),
                                                          FILE78 (nonlinear), DOFSMAP, FILE90, FILE91
ANALYS (+ .anl): FILE1|FILE1X/Y/Z|FILE9|FILE9lll + FILE3 + FILE4 [+ FILE77] [+ restart files]
          ──► FILE8 | FILE8X/Y/Z | FILE8xxx ; restart files COOXxxx, COOTKxxx, COOXI, COOTKI (if "Save Restart Files")
COMBIN: FILE81 + FILE82 ──► FILE8
MOTION: FILE8 (+ time history) ──► .TFU (computed ATF), .TFI (interpolated ATF), accelerations, ARS
RELDISP: .TFI (from MOTION) ──► relative displacements
STRESS: FILE8 [+ FILE78] ──► stresses/forces, STF ; FILE74 (nonlinear soil)
COMB_XYZ_STRAIN: FILE74 (X), FILE74 (Y), FILE74 (Z) ──► combined FILE74
Build_FILE77 (aux): per-level FILE77s ──► one FILE77
```

---

## 3. §4.1.2 Engineering considerations (items 1–20)

These items add to the modelling capabilities and limitations in Manual §1.5. Each item below has a paraphrase of the requirement, the numbers it uses, and what the implementation has to do.

### Item 1 — Rigid vs. flexible basement [CORE][UI]
* The basement is **always modelled as flexible** in the FE model. Assuming a rigid basement **saves no computation**, so model the basement with its actual visco-elastic properties.
* To study a rigid basement, run a **restart** (New Structure, Mode 2) with the basement **elastic modulus multiplied by 10⁴ and by 10⁵**, i.e. `E_rigid = 10⁴·E` and `E_rigid = 10⁵·E`. Running both values shows whether the result has converged to the rigid limit.
* Implementation: provide a UI helper that scales E for the selected groups by a factor. Run it as a New Structure restart.

### Item 2 — Surface vs. embedded structure [CORE]
* An embedded structure **can be treated as a surface structure** when embedment effects are small for the responses of interest. This gives a large run-time saving because there are far fewer interaction nodes.
* Caveat: ignoring embedment is **not conservative at every frequency**. It is broadly conservative for **structural forces**, but **ISRS** can show amplification and new peaks in narrow frequency bands.
* Implementation: no algorithm. The UI and documentation should state the caveat.

### Item 3 — 2D vs. 3D SSI [CORE]
* **2D SSI is not recommended**, especially at soil sites. It cannot represent the complex dynamic foundation stiffness correctly, and it usually has **too much radiation damping**, which makes it **unconservative**.
* 2D is acceptable for **sensitivity studies** of 2D layering effects: non-uniform or sloped layers, inclined baserock, topography.
* ACS SASSI can run **nonlinear (equivalent-linear) site response on 2D soil profiles**. This is also useful for building site-specific **directional incoherency models** [ADV].
* Implementation: 2D analysis uses POINT2 and plane-strain PLANE elements. Incoherency is **blocked** in 2D (see Item 4).

### Item 4 — Coherent vs. incoherent motion [ADV]
**Validity limits [CORE]**
* Incoherent motion is allowed **only for 3D models**. It is **not permitted** for quarter models, half models or 2D models. The validator must reject `coh = 1` when SYMM planes are active or when dim ≠ 3D.
* Wave-passage (directional) effects can be combined with incoherency (HOUSE `<wpass>` flag [xref §9.2.16]).

**Approaches.** There are seven incoherent SSI approaches. The breakdown below combines this chapter with xref §6 HOUSE/MOTION dialogs.

| # | Approach | Type | Notes |
|---|---|---|---|
| 1 | **Stochastic simulation, no phase adjustment** ("Simulation Mean", the reference approach in the 2007 EPRI study) | Stochastic | **Strongly recommended**, especially for large elastic foundations. It is accurate and efficient, and it gives response statistics. It is the "theoretically exact" Monte-Carlo approach recommended for SSSI. |
| 2 | Stochastic simulation **with phase adjustment** | Stochastic | Approximate upper bound. Phase adjustment is a MOTION option [xref]. |
| 3 | **AS** (algebraic sum of scaled incoherent modes) **with** ATF phase adjustment | Deterministic | The version used in EPRI 2007. |
| 4 | AS **without** phase adjustment | Deterministic | |
| 5 | **SRSS TF**, zero phase | Deterministic | Validated in EPRI 2007. Combines the **ATF amplitudes** of the incoherent modes by SRSS, assuming **zero phase** (as for a static response). ISRS are then computed from the combined ATF. |
| 6 | SRSS TF with phase = coherent-response phase | Deterministic | Not used by EPRI. |
| 7 | **SRSS FRS** | Deterministic | From the 1997 EPRI report (Tseng & Lilhanand). **Not validated** in 2007. Applies SRSS to the **end results** of each incoherent mode (ATF amplitudes, ISRS, maximum displacements, ZPA). |

Formulas (k = incoherent spatial mode index, N_m = number of modes used):
* SRSS TF: `|ATF_SRSS(f)| = sqrt( Σ_k |ATF_k(f)|² )`, with phase = 0. That combined complex ATF is then used to compute the ISRS.
* SRSS FRS: `R_SRSS = sqrt( Σ_k R_k² )`, where R is the final quantity (an ISRS ordinate, a ZPA, a peak displacement, or an ATF amplitude).

Guidance and warnings to implement as UI messages:
* AS and SRSS are **approximate**. They are validated **only for stick models on rigid mats**. Using them for flexible-foundation FE models needs a preliminary sensitivity study showing that ISRS at **floor corner locations** are accurate or conservative compared with stochastic simulation.
* SRSS needs a **separate SSI analysis for each incoherent spatial mode**. That makes it harder to use and slower than stochastic simulation. It exists mainly for benchmarking.
* Number of modes: **about 10 modes** is usually enough for **stick models on rigid mats** (EPRI 2007). **Large elastic mats can need tens to hundreds**. The vertical direction is the most critical because the mat is more flexible there.
* SRSS ignores the **coupling between closely spaced incoherent modal responses**. Deterministic SRSS also cannot choose the **sign (±) of each eigen-mode** uniquely, because the sign depends on the interaction-node numbering and changes the modal complex ATF.
* WARNING: for **structural forces**, SRSS TF would need STF modal responses. That is **not validated and not implemented** in ACS SASSI.
* **AS and SRSS are not recommended for licensing projects.**

Coherency models (plane-wave). These are listed here for nomenclature only; their equations belong to the HOUSE spec:

| Model # | Name | Applicability |
|---|---|---|
| 1 | Luco–Wong (1986) | Theoretical, not validated in practice |
| 2 | Abrahamson 1993 | All sites, surface foundations |
| 3 | Abrahamson 2005 | All sites, surface foundations |
| 4 | Abrahamson 2006 | All sites, embedded foundations |
| 5 | Abrahamson 2007 | Hard-rock sites, all foundations |
| 6 | Abrahamson 2007 | Soil sites, surface foundations |
| 7 | User-defined coherency model | — |

**Coherency matrix conditioning [CORE][ADV]**
* ACS SASSI builds the free-field **coherency matrix from the projections of all interaction nodes onto the horizontal plane**. If interaction nodes at different embedment levels project to nearly the same (x, y), the coherency matrix becomes **ill-conditioned**. A typical cause is a backfill mesh whose perimeter node planes lean slightly away from vertical.
* In such cases use the **"per level" approach**: compute the incoherent modes separately for each embedment level of interaction nodes. Then assemble the coherency eigenvectors for all interaction nodes with the auxiliary program **`Build_FILE77.exe`**. In the original program it is installed in `C:\ACSV300\Build_FILE77`; the manual also spells it "BuildFILE77".
  * Limitation: this is accurate **only if the incoherent modes at every level keep the same node-numbering pattern**, so that mode signs do not flip between levels. Checking this takes effort for every mode at every level.
  * [xref §6 HOUSE] After the per-level incoherent HOUSE runs, a final **coherent** HOUSE run with all interaction nodes is needed to build FILE4 correctly.
* "Per level" is also **required for deeply embedded foundations at soil sites**, where the coherence functions differ between the surface and the foundation level.
* WARNING / verification check: for multilevel embedment,
  1. request the extra HOUSE output on **incoherent mode contributions** (INCOH `<ipr>` = 1 [xref §9.2.17]), and
  2. check the **ANALYS ATF at the zero frequency (frequency number 1)**. In the input direction its amplitude must be **close to 1.00**. A large deviation means the coherency matrix is ill-conditioned.
* Implementation:
  * an automatic **"low-frequency ATF ≈ 1" check** on FILE8 [UI][CORE];
  * a **projected-node proximity check**: flag interaction-node pairs with different z and horizontal distance below a tolerance [DECISION: tolerance];
  * reporting of the eigenvalue spectrum or condition number of the coherency matrix.

### Item 5 — Isolated vs. multiple structures (SSSI) [CORE][ADV]
* **Structure-soil-structure interaction (SSSI)** can be analyzed in 3D by putting several structures (and their excavations) in one model (see Step 7).
* SSSI mainly affects **local responses**: ISRS, and seismic soil pressures on foundation walls and mats. **Different foundation levels** of neighbouring buildings usually amplify SSSI and raise wall pressures.
* Incoherency can increase SSSI effects at soft sites and can shrink inter-building gaps, possibly to impact. At soil sites the **minimum gap computed with incoherent motion is about 2–3 times the coherent value**.
* Incoherent SSSI should use **stochastic simulation with no phase adjustment**.

### Item 6 — Non-uniform seismic motion input [ADV]
* For sites that are strongly non-uniform horizontally (for example, boreholes with different soil columns), the input can vary across the foundation.
* The foundation, or each of several foundations, is divided into **zones**. Each zone gets its own excitation via **frequency-dependent amplitude amplification factors**. Incoherency and wave passage can be added on top.
* Limit: **up to 5,000 zones** under a foundation. This matches the HOUSE limit "Number of multiple support (foundation zones) = 5,000" [xref §1 limits]. The HOUSE `<me>` flag and the ME command [xref §9.2.20] control it. Note that the ME command text says 1–10 input motions (see Open Questions).

### Item 7 — Symmetry of the system [CORE]
* Symmetric structures under **symmetric or antisymmetric loading** can be analyzed as half or quarter models, which greatly cuts run time.
* The **SYMM** command's parameter gives the **number of symmetry planes**. [xref §9.2.39] The syntax is `SYMM,<no>,[<type>],[<node1>],[<node2>],[<node3>]`, with at most 2 planes. `type` 0 = symmetry, 1 = antisymmetry, both relative to the loading. In 3D the planes must be parallel to xz or yz; in 2D the line must be parallel to z.
* WARNING: **incoherent motion and impedance evaluation options do not work with half or quarter models.** The validator must block them.

### Item 8 — Rigid base rock vs. half-space [CORE]
* A visco-elastic **half-space** below the user-defined layers is simulated by having SITE add extra **computational half-space sublayers** automatically. Their total depth is

  `H_hs(f) = 1.5 · Vs_hs / f`   (= 1.5 λ_s, where λ_s is the half-space shear wavelength at frequency f)

  Vs_hs is the half-space shear-wave velocity and f is the analysis frequency.
* **Viscous dashpots** are also placed at the bottom of the extended deposit.
* The user input is the **number of generated half-space layers**:
  * recommended **10–20**;
  * the limit is **≤ 20** (SITE/POINT limit);
  * **0 or blank means no half-space simulation**, i.e. a **rigid base**.

  (See Item 20 for the method details.)

### Item 9 — SSI cut-off frequency [CORE]
* The cut-off frequency sets the **upper limit of the analysed frequencies**. It also controls **model size**, because the maximum allowable element sizes follow from it, and those set the size of the K and M matrices.
* It is governed by:
  (a) the frequency content of the input motion,
  (b) the dominant frequencies of the whole dynamic system,
  (c) the **time increment** of the input time history. Implied: the Nyquist frequency is `f_Nyq = 1/(2·DT)`.
* **Typical values: 30–40 Hz for soil sites; 60–70 Hz for rock sites.**
* Implementation, as a validator/advisor:
  * require `f_cut ≤ f_pass,min`, the minimum passing frequency over all soil layers and all excavation/near-field soil elements (Items 12, 19);
  * require `f_cut ≤ f_Nyq`;
  * advise `f_cut ≥ f_90%`, the frequency at which the cumulative fixed-base modal mass reaches 90 % in each direction (Step 3), when modal data have been imported.

### Item 10 — Selection of SSI frequencies [CORE][UI]
* What sets the choice: the frequency content of the input and of the system response, how narrow and how large the spectral peaks are, and how close together they are.
  * A first estimate comes from **fixed-base modal analysis**: the frequencies, the mode participation factors, and the shapes of the fixed-base transfer functions.
  * **Warning:** SSI can change the effective stiffness and damping a lot and shift the frequencies. Be careful when extrapolating fixed-base results.
* SSI tends to **flatten sharp peaks**, and frequency interpolation can **miss** structural peaks.
* Recommended number of SSI frequencies:

| Model | Coherent | Incoherent |
|---|---|---|
| Simple structures / stick models | **40–80** (no more is usually needed) | — |
| Complex FE models | **100–200** (Step 9 says 100–300) | **200–300** (**at least 200 from the start**) |
| Strategy | Can start with fewer (e.g. **80**) and refine | Start with ≥ 200 |

* The complex transfer functions (nodal acceleration ATF, and element stress/force STF) are **interpolated in the complex frequency domain** to **every Fourier frequency up to the cut-off**. ACS SASSI has **seven interpolation schemes**. [xref §6 MOTION "Interpolation Option"] They are:
  * 0 = SASSI2000 dense overlapping windows (weighted averaging),
  * 1 = original SASSI 1982 non-overlapping windows,
  * 2 = dense overlapping windows (averaging),
  * 3 = three overlapping windows,
  * 4 = non-overlapping windows with one position shift,
  * 5 = non-overlapping windows with two position shifts,
  * 6 = complex bicubic spline (needs a denser grid; recommended for incoherent analysis).
* **Frequency-addition rule [CORE]:** if an **interpolated peak lies near the midpoint** between two consecutive computed frequencies, and its amplitude is **much higher than the neighbouring computed amplitudes**, add a frequency there.
  * UI commands and macros can automate this. See the **CRITFREQ** command [xref §9.7.4: `CRITFREQ,<tol>,<minfilter>,<TF>,<Var>`, which compares `.TFI` peaks with `.TFU` values] and Demo 3.
* **TFU–TFI plot** [UI]: overlay the **computed ATF (`.TFU` files)** and the **interpolated ATF (`.TFI` files)**.

### Item 11 — Computation of SSI frequency points (NFREQ) and steps [CORE][IO]
* Transfer functions are computed only at frequencies that are **integer multiples of the frequency step DF**.
* **Maximum total number of SSI analysis frequencies = 500.**
* General (time-history) analysis:

  `DF = 1 / (DT · NFFT)`

  where DT is the time step of the input history and NFFT is the number of points in the Fourier transform. [xref §6 SITE: NFFT should be a power of 2; the UI rounds it to the nearest power of 2 and warns.]
* **Single-harmonic forced-vibration analysis:** no time history is needed, so the user gives **DF directly**.
* Integer frequency numbers:

  `NFREQ_i = f_i / DF`,  i = 1, 2, …, NF     (Eq. 4.1)

  NF is the total number of selected frequency points. The text says "according to item 6"; that reference is stale and Item 10 is meant.
* The largest frequency number comes from the cut-off:

  `NFREQ_NF = f_NF / DF`     (Eq. 4.2)

  where f_NF is the cut-off frequency.
* Implementation:
  * store the frequency set as **positive integers** (NFREQ) and compute `f_i = NFREQ_i · DF`;
  * [xref §6 SITE] sort the set ascending and **stop with an error on duplicates**;
  * enforce `NF ≤ 500`, and `NFREQ_i ≤ NFFT/2` (Nyquist);
  * when the user enters frequencies in Hz, convert with `round(f/DF)` and report the snapped value [DECISION].

### Item 12 — Soil deposit modelling [CORE]
* The soil is a stack of **horizontal, laterally semi-infinite elastic/visco-elastic layers** on **rigid base rock** or a **visco-elastic half-space**.
* **1/5-wavelength layer rule:** each layer's thickness must be **at most one fifth of the wavelength** at the highest analysis frequency:

  `h_j ≤ λ_min / 5 = Vs_j / (5 · f_cut)`,  equivalently  `f_pass,j = Vs_j / (5 · h_j) ≥ f_cut`.

  The reason is the **thin-layer method** assumption that displacements **vary linearly within each layer**. The geotechnical profile is therefore split into computational **sublayers** that meet the rule.
  * [xref §1.5.1] Use **more than 20 soil layers**; too few layers degrade the Rayleigh/Love modes from SITE.
* WARNING — numerical instability:
  * Deep soft deposits that are non-uniform with depth may need many layers to pass high frequencies for vertically propagating waves.
  * If properties are **non-uniform with depth**, or **Poisson's ratio > 0.47**, the free-field solution can become unstable at isolated frequencies.
  * Find those frequencies by inspecting the ATFs at several nodes. **Remove them from the SSI analysis**, and exclude them from the transfer-function interpolation in MOTION and STRESS. [xref §1.5.1: the aux tool `Remove_Frequencies_from_FILE8.exe`.]
  * Sensitivity runs at adjacent frequencies are highly recommended.
  * Implementation: warn when ν > 0.47 in any layer; provide a frequency-removal utility for FILE8.
* **Backfill / near-field soil:**
  * Irregular backfill next to the structure can go into the FE model.
  * To get **seismic soil pressures** on walls and mat, and to include local nonlinear hysteretic behaviour, add **one or a few rows of adjacent soil elements** around the walls and under the mat: **solids in 3D, plane elements in 2D**. These can **slow the run by up to about 2×**.
  * **Alternative:** put **springs** between the excavation-volume perimeter interaction nodes and the basement perimeter nodes. These are **duplicate nodes at the same coordinates**, generated automatically by the **INTGEN** command. The spring **axial force divided by the tributary (afferent) area** estimates the wall soil pressure: `p ≈ N_spring / A_trib`.
  * **Option A [ADV]** handles seismic pressure including **nonlinear foundation–soil separation and sliding**.
  * When **secondary (local) soil nonlinearity** is expected to matter (often the case for deeply embedded foundations), adjacent **nonlinear near-field soil elements** are **required**, together with **equivalent-linear SSI iterations** (Item 18).
* **Piles [ADV]:**
  * **Option NON** can model pile–soil interface slip with **nonlinear springs**.
  * The **SOLIDPILE** command generates the springs between the piles and the backfill automatically [xref §9.17.30: `SOLIDPILE,<group>,[stiff],[soft],[stiff2]`].

### Item 13 — Structural modelling [CORE]
* Model the structure with 2D or 3D finite elements that follow its geometry. Element library: 3D **BEAM, SHELL, SOLID, SPRING**, 2D **PLANE**, and 3D **General Matrix (GM)** elements. GM elements carry substructure or super-element matrices, or fluid inertial effects.
* **Option AA [ADV]** converts **ANSYS MATRIX50** super-elements (`.sub` files) to ACS SASSI GM elements (`.pre` files).
* **Interaction node rules [CORE]**. The validator must enforce all of these:
  * a. Interaction nodes are defined **on the excavated soil model, not on the structure**.
    * If the outer surface of the excavated volume shares nodes with the structure (as it should when the FE model has no surrounding backfill), those **common nodes are interaction nodes**.
    * **No other structural node may be an interaction node**, including basement nodes.
    * **No internal excavation-volume node may be connected to a structural node.**
  * b. Every interaction node must be **below the ground surface** and **lie on a soil-layer interface**. (For surface foundations "below" has to mean "at or below"; see Open Questions.)
  * c. For embedded models, **interaction nodes must be numbered in ascending order**.
* WARNING: run **EXCSTRCHK** [xref §9.8.1]. It lists excavation **interior** nodes that are wrongly shared with structure basement nodes. It is a very important model check. [xref §1.5.1] The **FIXEDINT** check (interaction nodes with fixed translations) and the **HINGED** check are also recommended.

### Item 14 — Excavated soil volume modelling [CORE]
* For embedded structures, the excavated soil volume is modelled with **solid elements (3D)** or **plane elements (2D)**. If there is no backfill zone, it connects to the structure at the **perimeter interaction nodes**.
* The excavated element sizes are set by the **spacing of the interaction nodes**.
* **Internal basement nodes** that belong to internal structures vibrating independently of the soil must be **kept separate** from the SSI interaction nodes. In other words, use separate meshes for the basement and the excavation [xref §1.5.1 item 11].
* SSI results for embedded models are **sensitive to the excavation mesh size**:
  * **vertical** size: the 1/5-wavelength rule;
  * **horizontal** size: the 1/5 rule can be too strict (see Item 19).

### Item 15 — SSI substructuring methods [CORE]
* There are three main flexible-volume substructuring methods: **FV** (flexible volume / direct), **FI** (flexible interface: **SM** = FI-FSIN subtraction, **MSM** = FI-EVBN modified subtraction), and **FFV** (fast flexible volume). [xref §1.5.1, §2.1]
* WARNING: **ASCE 4-2016 and USNRC SRP 3.7.2** require a **preliminary validation against FV** before **SM, MSM or FFV** are used in production.
  * The validation compares **ATFs at the common nodes** between the structure and the excavated soil.
  * The simplest form is an "excavated soil model" with no structure (the "swimming pool model"). ACS SASSI suggests a **simplified massless foundation model** instead, to capture kinematic SSI.
  * For **deeply embedded SMR-type models in deep soft soil**, the excavation-only model can be **poorly conditioned** and should be avoided.
  * The most complete validation is the **full SSI model**, checking ATFs at critical locations.
* WARNING: ASCE 4-16 and SRP 3.7.2 call for this validation with **quarter models of the excavation**. However, **quarter models are much more numerically stable** than full models for unstable SM/MSM embedded models, because of the kinematic symmetry/antisymmetry conditions. A quarter model can therefore **hide** an instability. [xref §1.5.1: ACS SASSI does *not* recommend quarter models for qualifying SM/MSM.]
* Implementation: provide a "method validation" workflow. Run the same model with FV and with the chosen method, then plot and compare the ATFs at the selected (common) nodes [UI].

### Item 16 — Addition of SSI frequencies [CORE][IO]
* Frequencies can be added or combined with **COMBIN**, **as long as the frequencies are in the computed FILE8 data files**. FILE4 (`modelname.N4`) is **frequency-independent**, so it never needs regenerating for added frequencies. The frequencies to be added must already be in the frequency-dependent inputs, i.e. **FILE1 (or FILE9) and FILE3**, so SITE/POINT (or FORCE) must have been run for them.
* Manual example (paraphrased):
  * SITE and POINT were run for 10 frequencies: 0.98, 2.93, 4.88, 6.84, 8.79, 10.74, 12.70, 14.65, 15.62, 17.58 Hz.
  * ANALYS was run for only 5 of them.
  * Two more frequencies (stated as "2.93 and 15.66 Hz") are then solved with ANALYS, and COMBIN merges them into a new FILE8.
  * Check: these 10 values are exactly `k · 0.9765625 Hz` for k = 1, 3, 5, 7, 9, 11, 13, 15, 16, 18, i.e. DF = 1/1.024 Hz. "15.66" is evidently 15.62(5).
* COMBIN input/output convention [xref §4.2.2]: the two old FILE8s must be **renamed FILE81 and FILE82**. The output is a new **FILE8**.
* COMBIN algorithm [DECISION, consistent with the manual]: take the union of the frequency sets, sorted ascending. The manual does not say how a frequency present in both files is resolved (see Open Questions).

### Item 17 — Flexibility and impedance matrices of the foundation [CORE]
* The **free-field flexibility (compliance) matrix** comes from the **point-load solution of a soil column**. The column is a FE model of **plane-strain elements (2D)** or **axisymmetric elements (3D)**, with transmitting boundaries from the SITE Mode 1 eigen-solution (POINT2/POINT3).
* The **impedance matrix** C, with real part K and imaginary part ωD, is the **inverse of the compliance matrix**, frequency by frequency:

  `K + i ω D = (f + i g)⁻¹`     (Eq. 4.3)

  Here ω is the circular frequency of analysis (rad/s, ω = 2πf_Hz), i = √−1, and **f** and **g** are the real and imaginary parts of the free-field compliance matrix.
  * Naming clash: this `f` is the compliance matrix, not the frequency. In code, use names like `F_re` and `F_im`.
  * In Ch. 2 notation, this C is `X_ff`.
* **Columns of the compliance matrix** come from applying a **unit harmonic force (or moment), one at a time,** at **each translational DOF of the interaction nodes**. Interaction nodes have translational DOFs only [xref §1.5.1 item 5]. The compliance is then **inverted at every frequency**.
* The compliance and impedance matrices are **fully populated (dense)**, so the inversion has to be efficient and parallel.
* **Memory:** about **14.4 GB for 10,000 interaction nodes**, growing with the **square** of the node count. This is consistent with a dense complex double-precision matrix of size 3N × 3N:

  `Mem ≈ (3 · N_int)² × 16 bytes`. For N_int = 10,000 that is 9·10⁸ × 16 B = 1.44·10¹⁰ B = 14.4 GB.

  The implementation should show this estimate before running.
* Implementation notes:
  * use complex128;
  * solve or factorize (LU, or complex-symmetric LDLᵀ) instead of forming the explicit inverse where possible. Mathematically the result must equal Eq. 4.3.
  * The matrix should be symmetric by reciprocity. Check this, and symmetrize if needed [DECISION].

### Item 18 — Nonlinear hysteretic soil behaviour [CORE][ADV]
* The complex frequency-domain solution only works for **linear (linearized) systems**. Nonlinear soil is approximated with an **equivalent-linear iterative procedure (SHAKE methodology)** for both free-field and SSI analyses.
* There are two kinds of soil nonlinearity:
  * **Primary (global):** free-field wave propagation. This is handled by SOIL → SITE, and is usually the only one considered.
  * **Secondary (local):** caused by SSI. It matters when the model has a **soft backfill zone** or **several heavy neighbouring buildings (SSSI)**, and it affects how soil pressure is distributed on embedded walls.
* Secondary nonlinearity needs an **extended near-field soil zone** in the FE model. The effective shear strains in that zone are found iteratively, following the equivalent-linear procedure.
* **Each SSI iteration is a "New Structure" restart.** This is about **2–4 times faster** than an initiation run.
* After the X, Y and Z directional runs of an iteration, the **three-directional effective strain** is computed with the auxiliary **`COMBIN_XYZ_STRAIN.exe`**, also called **COMB_XYZ_STRAIN** elsewhere. It applies the **SRSS rule** to the effective shear strains that STRESS wrote to **FILE74** for each direction:

  `γ_eff,e = sqrt( γ_eff,e,X² + γ_eff,e,Y² + γ_eff,e,Z² )`  for each nonlinear soil element e.

  (The iteration loop is in §4.5.)

### Item 19 — Excavated soil and structure FE model discretization [CORE]
* **Vertical mesh size:** the largest dimension of each element must be **≤ 1/5 of the shortest wavelength**. That wavelength corresponds to the **cut-off frequency**:

  `Δz_max ≤ Vs / (5 · f_cut)`.

  The rule is acceptable because the mass matrix is **50 % lumped + 50 % consistent** [xref §1.5.1 item 12: except beams, which use consistent mass, and plates, which use lumped mass].
* The 1/5 rule governs the sizing of:
  * the **far-field soil layers**,
  * the **near-field soil**,
  * the **excavated soil elements in the vertical direction**.

  It is the right rule for vertical sizes when the input is mainly **vertically propagating waves**.
* **Horizontal mesh size of the excavation.** Inside the excavation, scattering produces a mix of incident and scattered body waves plus **surface waves**, which travel horizontally at close to Vs. In principle the horizontal size should therefore be **close to or equal to the vertical size**. In practice:
  * for realistic layered deposits, **about 1.2–1.5 times, sometimes up to 2 times, the vertical size** is often acceptable, case by case;
  * if the wave field were **purely 1D vertical propagation**, the horizontal size would be **unrestricted**, and could be several times the vertical size;
  * the **larger the mismatch** between the complex dynamic stiffness of the basement and that of the excavated soil, the **more scattering**, and the **finer** the horizontal mesh must be.
* WARNING: always run **mesh sensitivity studies** to justify the horizontal size. Larger horizontal sizes cut run time a lot. Such studies can use **quarter FE models of the excavated volume with no structure**. **Without a sensitivity study, take the horizontal size ≈ the vertical size.**
* Implementation, as validator rules:
  * vertical element size ≤ Vs/(5 f_cut), using the **strain-compatible Vs** of the layer that holds the element [DECISION];
  * warn if horizontal size > 1.0 × vertical (default ratio limit 1.0; user-adjustable up to 2.0) [DECISION];
  * report the element **passing frequency** `f_pass = Vs/(5·h_max)`.

### Item 20 — Half-space simulation [CORE]
Two techniques are used together, both in **SITE**:
1. **Variable-depth method.**
   * Up to **20** extra half-space computational layers are added under the profile. They have the **half-space properties** and a **total thickness of 1.5 λ**, where `λ = Vs_hs / f` is the half-space shear wavelength. This depth therefore **changes with frequency**.
   * Why 1.5 λ: the fundamental Rayleigh mode in a half-space decays with depth and is essentially gone at 1.5 λ.
   * The 1.5 λ depth is split into **n layers whose thickness increases with depth**. **n = 20** is suggested for best accuracy.
   * As a result, the sublayer thickness **increases with depth and decreases with frequency**.
2. **Viscous boundary.** The rigid base of the extended layer system is replaced by a **viscous boundary**, with **dashpots in the horizontal and vertical directions**.

* Implementation [DECISION, not specified in the manual]:
  * Sublayer thicknesses, with H = 1.5 Vs_hs / f and n sublayers: `h_j = H · (2j − 1) / n²`, j = 1..n. These increase linearly with depth and sum exactly to H.
    * For n = 20 the deepest sublayer is 0.146 λ, which meets the 1/5 rule.
    * The deepest sublayer is `1.5(2n−1)/n² · λ`. It meets the 1/5 rule only for n ≥ 15 (n = 15 gives 0.193 λ; n = 14 gives 0.207 λ), so for n = 10–14 it breaks the rule. A geometric progression with a capped last layer is an alternative.
  * Base dashpots per unit area: Lysmer–Kuhlemeyer type, `c_h = ρ·Vs_hs` and `c_v = ρ·Vp_hs`.
  * Both choices must be documented and verified against analytical half-space solutions (§7).

---

## 4. §4.2 ACS SASSI runs (run orchestration)

### 4.1 Run types [CORE][UI]
First decide which run type applies:
* **a. Initiation** SSI solution runs
* **b. Post-processing** runs (after FILE8 exists)
* **c. Restart / re-analysis** SSI solution runs (ANALYS restart)

[xref Ch. 3.2] There are also changes that need **no ANALYS** run:
* "**New Time History**": only the control-motion time history or response spectrum changes, and the wave field stays the same. Re-run **MOTION** only.
* A new force-history shape with the same load pattern: MOTION only.

### 4.2 §4.2.1 Initiation run for a seismic linearized SSI solution [CORE]
Modules run in sequence:

| # | Module run | Produces | Consumers |
|---|---|---|---|
| 1 | **EQUAKE** | Spectrum-compatible seismic input accelerations | SOIL, MOTION, RELDISP, STRESS |
| 2 | **SOIL** | Iterated effective soil properties from the nonlinear free-field (SHAKE) analysis | SITE |
| 3 | **SITE** (Mode 1 and 2) | Free-field soil wave modes and soil-layering information (FILE2, FILE1) | POINT, HOUSE, ANALYS |
| 4 | **POINT** (POINT2 or POINT3) | Soil-layer flexibility/compliance results (FILE3) | ANALYS |
| 5 | **HOUSE** | Structural complex K and M matrices (FILE4/`.N4`, COOSK, COOSM), and the **eigen-solution for incoherent motion** (FILE77) | ANALYS, STRESS |
| 6 | **ANALYS** (Mode 1) | **Seismic:** complex **ATF** (acceleration transfer function) solution. **Vibration / external force:** **DTF** (displacement transfer function) solution. Stored in FILE8. | MOTION, RELDISP, STRESS |

* Seismic ANALYS runs can cover **one direction (X)** or **three directions (X, Y, Z)** at once.
* Vibration ANALYS runs can cover **1 load case or up to 500 load cases**, depending on model size and RAM.
* "ANALYS initiation and restart are done in the same way" when it comes to multiple directions and cases.
* [xref §6 ANALYS dialog] When the RAM is insufficient, the original program fails with an "access violation". The re-implementation should estimate memory beforehand and refuse the run with a clear message [UI].

### 4.3 §4.2.2 Post-processing runs [CORE][IO]
Once FILE8 exists:

| # | Run | Purpose |
|---|---|---|
| 1 | **COMBIN** | Merge FILE8s with different SSI frequency sets into a new FILE8. Needed **only** when new frequencies are added. The inputs **must be renamed FILE81 and FILE82**. |
| 2 | **MOTION** | Accelerations, ATF (interpolated `.TFI`, computed `.TFU`), ARS (acceleration response spectra) |
| 3 | **RELDISP** | Relative displacements. **Needs MOTION run first**, because it reads the complex ATF `.TFI` files for the nodes of interest. |
| 4 | **STRESS** | Element stresses and forces, STF |

MOTION and STRESS are **independent** of each other.
Dependency rules for the orchestrator:
* COMBIN (optional) runs before MOTION and STRESS;
* MOTION runs before RELDISP;
* STRESS can run after FILE8 alone.

### 4.4 §4.2.3 ANALYS restart modes [CORE]

| ANALYS "Mode of Analysis" | Name | When to use | What changes | What is re-used |
|---|---|---|---|---|
| **Mode 1** | Initiation | New problem | — | — |
| **Mode 2** | **New Structure** | Structure or near-field soil properties/geometry change; **interaction (embedment) nodes and soil layering unchanged** | New FILE4 (HOUSE re-run) | Soil impedances (restart files) and the old FILE1/FILE9 |
| **Mode 3** | **New Seismic Environment** (seismic) | Type or location of the seismic input changes, or a new incoherent simulation | New FILE1 (SITE re-run) and/or new FILE77 (HOUSE re-run) | Factorized system + impedances (restart files) |
| **Mode 3** | **New Dynamic Loading** (external loads) | New external force set | New FILE9 (FORCE re-run) | Factorized system + impedances |
| *Mode 6* | *New Load Vector* [xref §6 ANALYS; batch only, not in GUI] | Special changes to the seismic load vector (e.g. per-level incoherent embedded analyses) | — | — |

The manual says "each of the above modes involves only two computer runs", but the lists that follow have three runs (see Open Questions).

Theory of re-use [CORE, inferred from Ch. 2 + Ch. 6 ANALYS dialog]:
* The impedance `X_ff(ω)` depends only on the **soil layering and the interaction-node geometry and numbering**.
  * **Mode 2** can reuse X_ff and re-assemble and refactor only with the new structure matrices `C(ω) = K − ω²M`.
  * The ANALYS dialog says the "Save Restart Files" option keeps **only COOXxxx** for New Structure. That suggests COOXxxx stores the impedance-related matrix at frequency xxx.
* For Mode 3, the total system matrix is unchanged and **already factorized ("triangularized")**. Only the **right-hand side** changes:
  * seismic: `X_ff · u'_f` with new free-field motions, or a new incoherent field;
  * external: the FILE9 load vector.

  That suggests COOTKxxx stores the triangularized total stiffness at frequency xxx.

### 4.5 §4.2.4 New Structure (or near-field strain-dependent soil properties) restart [CORE]
Used for **seismic and vibration** problems when the FE model changes **without changing the embedment (interaction) nodes**.

Runs:
1. **HOUSE** (new FE model) → new **FILE4** (`modelname.n4`)
2. **ANALYS** restart, **"New Structure", Mode 2**. Inputs: the new FILE4, the **old FILE1** (seismic) or **old FILE9** (vibration), and the restart files. Output: new **FILE8**.
3. **STRESS** (stresses in the new FE model)

**Nonlinear soil SSI (equivalent-linear SSI iterations) [ADV]** use this same restart. Each iteration runs **HOUSE, ANALYS and STRESS**:

```
iteration k = 0: initiation run (linear, initial near-field soil properties)
loop k = 1, 2, ...
    for dir in (X, Y, Z):
        STRESS(dir)            -> FILE74(dir)   (effective shear strain / strain-compatible G, D per nonlinear soil element)
    COMB_XYZ_STRAIN            -> FILE74 (SRSS-combined over X, Y, Z; Item 18 formula)
    HOUSE (reads FILE74 for nonlinear soil groups; writes FILE78) -> new FILE4
    ANALYS Mode 2 (X, Y, Z)    -> new FILE8X/Y/Z
    check convergence using FILE74 and FILE78 text files
until converged
```

* **FILE74** is a text file written by STRESS and read by HOUSE. **FILE78** is a text file written by HOUSE and read by STRESS. Both are **written at every iteration** and are used to **check convergence**.
* [xref Ch. 3.2] HOUSE creates a **`.liq`** file on the initial run. When that file contains **1**, HOUSE takes its material input for the next iteration from **FILE74**. A **`.pin`** file is also part of the nonlinear input.
* The manual does not give a convergence criterion or a strain ratio (the effective/maximum strain factor) in this chapter (see Open Questions).

### 4.6 §4.2.5 New Seismic Environment, or incoherent simulation [CORE][ADV][IO]
Seismic problems only. Used when the **type of seismic input** (wave field composition) or **its location** changes, or for incoherent simulations.

Runs:
1. **SITE**, for the new seismic environment, or once each for the X, Y and Z directions
2. **HOUSE**, for a new incoherent stochastic simulation run for X, Y, Z. This is needed only when incoherency is considered, because incoherency changes the seismic loads [xref Ch. 3.2 b].
3. **ANALYS** restart, **"New Seismic Environment", Mode 3**

The FILE1 from SITE is the ANALYS input. ANALYS writes a new FILE8, or a set of FILE8s.

**Coherent, three directions at once:**
* Before running ANALYS, create **FILE1X, FILE1Y, FILE1Z** by running SITE once per direction and copying or renaming FILE1.
* ANALYS then produces **FILE8X, FILE8Y, FILE8Z** in one run.
* [xref §6 ANALYS]
  * Set "Simultaneous Cases" = 1 for coherent analysis.
  * Wave types per direction in `.sit`: X uses SV waves (x' direction, angle 0), Y uses SH waves (y' direction, angle 0), Z uses P waves (z direction, angle 0).
  * The coordinate transformation angle in `.anl` must be 0.

**Incoherent simulations:**
* **HOUSE can run up to 50 incoherent motion simulations** in one run. It creates up to **50 FILE77s**, named **FILE77001 … FILE77050**.
* **ANALYS can run up to 50 simulations × 3 directions in one run**, producing **up to 150 FILE8s** depending on model size and RAM. They are named **FILE8001, FILE8002, FILE8003** for simulation 1 (X, Y, Z), FILE8004–FILE8006 for simulation 2, and so on:

  `FILE8{3(s−1)+d:03d}`, where s = 1..50 is the simulation and d = 1 (X), 2 (Y), 3 (Z).

* [xref §6 ANALYS] The number of simulations given to ANALYS must equal the number used in HOUSE.

**Restart files required by ANALYS** (written during the initiation run when **"Save Restart Files"** is checked):

| File | Meaning |
|---|---|
| **COOXxxx** | Per-frequency restart database. xxx = **3-digit frequency order number** (COOX001 is the 1st frequency). |
| **COOTKxxx** | Per-frequency restart database ("TK" = triangularized stiffness, inferred) |
| **COOXI**, **COOTKI** | Index files for the above |
| **DOFSMAP** | DOF map (produced by HOUSE) |
| **FILE90**, **FILE91** | Additional restart files (produced by HOUSE) |

* The text in this section says "xxxxx" (five characters), but Ch. 3 and Ch. 6 say three digits. **Use 3 digits**, which is enough for 500 frequencies.
* Because the files are numbered by order, a restart has to use the **full sequential frequency set**. Using a subset means renaming COOXxxx/COOTKxxx accordingly [xref Ch. 3].

### 4.7 §4.2.6 New Dynamic Loading [CORE][IO]
Vibration (external force) problems only.

Runs:
1. **FORCE** → new **FILE9**
2. **ANALYS** restart, **"New Dynamic Load", Mode 3** → new **FILE8**

**Multiple load cases:**
* ANALYS can run **up to 500 load cases** in one restart.
* All the FILE9s must exist before ANALYS runs. They are named **FILE9001, FILE9002, …, FILE9500**: `FILE9{l:03d}`, l = 1..500.
* [xref §6 ANALYS] ANALYS writes one **FILE8xxx** per load case.

---

## 5. Formula quick reference [CORE]

| # | Formula | Variables / units |
|---|---|---|
| (4.1) | `NFREQ_i = f_i / DF`, i = 1..NF | f_i in Hz; NFREQ_i a positive integer |
| (4.2) | `NFREQ_NF = f_NF / DF` | f_NF = cut-off frequency (Hz) |
| — | `DF = 1 / (DT · NFFT)` | DT time step (s); NFFT a power of 2 |
| — | `NF ≤ 500` | Maximum number of SSI frequencies (SITE/POINT/HOUSE/ANALYS) |
| (4.3) | `K + iωD = (f + i g)⁻¹` | f, g = real and imaginary parts of the free-field compliance at the interaction-node translational DOFs; ω = 2πf_Hz |
| — | `h ≤ Vs/(5·f_cut)`; `f_pass = Vs/(5·h)` | Layer thickness, or vertical element size |
| — | `h_H ≈ h_V` (default); `≤ 1.2–1.5 (max ~2) × h_V` once sensitivity studies justify it | Horizontal excavation element size |
| — | `H_hs(f) = 1.5·Vs_hs/f = 1.5 λ` | Total half-space buffer depth; n = 10–20 sublayers (≤ 20); n = 0 means rigid base |
| — | `E_rigid = 10⁴·E` and `10⁵·E` | Rigid-basement sensitivity (New Structure restart) |
| — | `Mem ≈ (3 N_int)² · 16 B` (= 14.4 GB at N_int = 10,000) | Dense complex impedance or compliance matrix |
| — | `γ_eff = sqrt(γ_X² + γ_Y² + γ_Z²)` | COMB_XYZ_STRAIN, per nonlinear soil element |
| — | `|ATF_SRSS| = sqrt(Σ_k |ATF_k|²)`, phase 0 | SRSS TF incoherent approach |
| — | `R = sqrt(Σ_k R_k²)` | SRSS FRS incoherent approach (end results) |
| — | FILE8 index = `3(s−1)+d` | Incoherent multi-simulation naming |

Recommended practice values:

| Parameter | Recommended value |
|---|---|
| Cut-off frequency | 30–40 Hz (soil), 60–70 Hz (rock) |
| SSI frequencies | 40–80 (sticks); 100–200 (complex FE, coherent; Step 9 says 100–300); 200–300 (incoherent, ≥ 200 from the start); coherent may start at about 80 |
| Half-space sublayers | 10–20; 20 suggested for best accuracy |
| Soil layers | > 20 [xref §1.5.1] |
| SRSS incoherent modes | About 10 for stick models on rigid mats; tens to hundreds for flexible mats |
| Poisson's ratio warning threshold | > 0.47 |
| Backfill element run-time penalty | Up to about 2× |
| Nonlinear SSI iteration run time | Each New Structure restart is about 2–4× faster than the initiation run |
| Incoherent SSSI gap | About 2–3× the coherent minimum gap (soil sites) |
| Non-uniform input zones | Up to 5,000 |
| Incoherent simulations per run | ≤ 50 (HOUSE and ANALYS) |
| Load cases per ANALYS run | ≤ 500 |

---

## 6. Validator and model-check catalogue derived from this chapter [UI][CORE]

The program should implement these checks. "E" means error (block the run) and "W" means warning.

| ID | Check | Level | Source item |
|---|---|---|---|
| G-01 | Number of SSI frequencies NF ≤ 500; frequency numbers are positive integers, unique, and sorted ascending | E | 11, [xref SITE] |
| G-02 | `max(NFREQ)·DF ≈ f_cut`; `NFREQ ≤ NFFT/2` | W/E | 11 |
| G-03 | NFFT is a power of 2 (round and warn) | W | [xref SITE] |
| G-04 | Recommended frequency count for the model type: < 40 (stick), < 100 (FE coherent), < 200 (incoherent) | W | 9, 10 |
| G-05 | Each soil layer: `Vs/(5h) ≥ f_cut` (layer passing frequency) | W (offer auto-sublayering) | 12, 19 |
| G-06 | Each excavation / near-field soil element: vertical size ≤ Vs/(5 f_cut) | W | 14, 19 |
| G-07 | Horizontal-to-vertical excavation size ratio > user limit (default 1.0, max advised 2.0) | W | 19 |
| G-08 | f_cut ≤ min passing frequency of the mesh and of the layering, and ≤ 1/(2DT) | W | 3, 9 |
| G-09 | Poisson's ratio > 0.47 in any layer | W | 12 |
| G-10 | Number of half-space layers is 0 (rigid base) or 10–20; > 20 is an error | W/E | 8, 20 |
| G-11 | 20 or fewer soil layers (the manual recommends more than 20) | W | [xref §1.5.1] |
| G-12 | Interaction nodes: at or below the ground surface, at a soil-layer interface elevation (within tolerance) | E | 13b |
| G-13 | Interaction nodes in ascending order (embedded models) | E | 13c |
| G-14 | No structural node is an interaction node unless it is a common node on the outer surface of the excavation | E | 13a |
| G-15 | **EXCSTRCHK**: no excavation interior node is shared with a structure element | E | 13 |
| G-16 | FIXEDINT: no interaction node has fixed translations | W | [xref §1.5.1] |
| G-17 | Incoherency (`coh = 1`) only with 3D and no SYMM planes; impedance-evaluation options not allowed with SYMM | E | 4, 7 |
| G-18 | Projected horizontal near-coincidence of interaction nodes at different elevations (coherency ill-conditioning); suggest "per level" | W | 4 |
| G-19 | Post-run: \|ATF\| at the first (lowest) frequency in the input direction ≈ 1.00 at all nodes | W | 4 |
| G-20 | Post-run: interpolated (TFI) peak near a mid-interval and much higher than the neighbouring computed (TFU) points; suggest added frequencies (CRITFREQ) | W | 10 |
| G-21 | Non-FV method (SM/MSM/FFV/FI) selected without a recorded FV validation | W | 15 |
| G-22 | Memory estimate for the dense impedance matrix vs available RAM | W/E | 17 |
| G-23 | Restart prerequisites: Mode 2 needs the restart files with the interaction nodes and soil unchanged (compare a hash of the interaction-node coordinates/numbering and the layering); Mode 3 needs COOXxxx/COOTKxxx, COOXI, COOTKI, DOFSMAP, FILE90, FILE91, and a new FILE1 or FILE9 | E | 4.2.3–4.2.6 |
| G-24 | COMBIN: FILE81 and FILE82 exist and come from the same model and DF | E | 16 |
| G-25 | RELDISP requested without MOTION `.TFI` files | E | 4.2.2 |
| G-26 | Simultaneous X/Y/Z coherent run: FILE1X/FILE1Y/FILE1Z exist | E | 4.2.5 |
| G-27 | Number of incoherent simulations ≤ 50, and the same in HOUSE and ANALYS | E | 4.2.5 |
| G-28 | FILE9001…FILE9nnn exist for the requested number of load cases (≤ 500) | E | 4.2.6 |
| G-29 | Forced-vibration analysis: SITE Mode 2 skipped; SOIL not used; initial soil properties used | W | 13 |

---

## 7. Verification examples (worked numbers for unit tests) [CORE]

1. **Frequency step:** DT = 0.005 s, NFFT = 4096 → DF = 1/(0.005·4096) = 0.048828125 Hz. With f_cut = 25 Hz, NFREQ_NF = 25/0.048828125 = **512**.
2. **Item 16 frequency grid:** DF = 0.9765625 Hz (e.g. DT = 0.001 s, NFFT = 1024). Frequency numbers {1, 3, 5, 7, 9, 11, 13, 15, 16, 18} give {0.977, 2.930, 4.883, 6.836, 8.789, 10.742, 12.695, 14.648, 15.625, 17.578} Hz, which are the manual's values after rounding.
3. **1/5 rule:** Vs = 300 m/s, f_cut = 33.3 Hz → λ = 9.0 m → h_max = 1.8 m. Conversely, a 2 m layer with Vs = 400 m/s has f_pass = 400/(5·2) = 40 Hz.
4. **Half-space buffer:** Vs_hs = 1500 m/s.
   * At f = 5 Hz: H = 1.5·1500/5 = 450 m.
   * At f = 30 Hz: H = 75 m.
   * With n = 20 and the [DECISION] linear-increase rule, the sublayers at 30 Hz run from 75/400 = 0.1875 m (top) to 75·39/400 = 7.31 m (bottom). The bottom value is ≤ λ/5 = 10 m.
5. **Impedance memory:** N_int = 10,000 → (30,000)²·16 B = 14.4 GB. N_int = 20,000 → 57.6 GB (scales with the square of N).
6. **Strain SRSS:** γ_X = 0.10 %, γ_Y = 0.08 %, γ_Z = 0.02 % → γ_eff = sqrt(0.01 + 0.0064 + 0.0004) % = 0.1296 %.
7. **Incoherent FILE8 naming:** simulation 17, direction Y → index 3·16 + 2 = 50 → `FILE8050`. Simulation 50, direction Z → `FILE8150`.
8. **Low-frequency ATF check:** for any coherent, correctly built model, the X-input ATF in X at all nodes at the first frequency is 1.00 ± a small tolerance. Use this as a regression test.
9. **Recommended additional benchmarks (not from the manual)** for checking Eq. 4.3 and Item 20:
   * static-limit impedance of a rigid circular surface footing on a uniform half-space: `K_v = 4GR/(1−ν)`, `K_h = 8GR/(2−ν)`;
   * frequency-dependent impedance against published Luco/Veletsos values;
   * free-field 1D amplification of a uniform layer on rigid base, with peaks at `f_n = (2n−1)Vs/(4H)`.

---

## 8. File naming conventions from this chapter [IO]

| Name | Producer → consumer | Description |
|---|---|---|
| `.sit` | user/UI → SITE | Seismic wave-field description, control motion location |
| `.poi` | user/UI → POINT | Maximum embedment, point-load (central zone) radius |
| FILE2 | SITE Mode 1 → SITE Mode 2, POINT | Soil-layer eigen-solutions (transmitting boundary) |
| FILE1 (FILE1X/Y/Z) | SITE Mode 2 → ANALYS | Free-field motions (per direction) |
| FILE3 | POINT → ANALYS | Point-load (compliance) solution |
| FILE9 (FILE9001–FILE9500) | FORCE → ANALYS | External load vectors (per load case) |
| FILE4 = `modelname.N4` | HOUSE → ANALYS | Complex K and M of structure and excavated soil (frequency-independent) |
| COOSK, COOSM | HOUSE → ANALYS | FE model stiffness and mass matrix files |
| FILE77 (FILE77001–FILE77050) | HOUSE → ANALYS (STRESS) | Incoherent-motion eigen-solution / simulations |
| FILE74 | STRESS → HOUSE, COMB_XYZ_STRAIN | Nonlinear soil effective strains / properties (text) |
| FILE78 | HOUSE → STRESS | Nonlinear iteration data (text) |
| DOFSMAP, FILE90, FILE91 | HOUSE → ANALYS (restart) | Restart support files |
| COOXxxx, COOTKxxx, COOXI, COOTKI | ANALYS → ANALYS (restart) | Per-frequency restart databases plus their index files |
| FILE8 (FILE8X/Y/Z; FILE8001–FILE8150; FILE8xxx per load case) | ANALYS → COMBIN, MOTION, RELDISP, STRESS | SSI transfer-function database |
| FILE81, FILE82 | user rename → COMBIN | The two FILE8s to merge |
| `.TFU` | MOTION | Computed (uninterpolated) ATF |
| `.TFI` | MOTION → RELDISP | Interpolated ATF at the Fourier frequencies |
| `Build_FILE77.exe` | aux | Assembles per-level FILE77s |
| `COMBIN_XYZ_STRAIN.exe` / COMB_XYZ_STRAIN | aux | SRSS combination of directional FILE74s |
| `Remove_Frequencies_from_FILE8.exe` | aux [xref §1.5.1] | Removes spurious frequencies from FILE8 |

Implementation: keep these exact names (case-insensitive on input; the manual uses both `.N4` and `.n4`), for nomenclature fidelity. Wrap them in a Python "run directory" abstraction that resolves names such as `file8(sim=s, dir='Y')`.

---

## 9. Open questions / ambiguities

1. **Step 13e lettering.** The text says "Item d of Step 10 is replaced by … FORCE … FILE9, which replaces FILE1 in item e of Step 10", but Step 10 has only items a–d. *Interpretation:* FORCE replaces the site-response part (SITE Mode 2), and FILE9 replaces FILE1 as an ANALYS input (Step 10d). Ch. 3.2 gives the order SITE → POINT → FORCE → HOUSE → ANALYS.
2. **Item 11 references "item 6 of this section"** for choosing NF. This is stale and Item 10 is meant.
3. **Item 16 example** lists "2.93 and 15.66 Hz" as *new* frequencies, but 2.93 Hz is already one of the 10 SITE/POINT frequencies, and 15.66 Hz is not on the DF grid (15.625 Hz is). It also says "FILE9 and FILE3" for what looks like a seismic case; FILE1 would be expected. *Decision:* treat the example as illustrative. The requirement is that added frequencies exist in the frequency-dependent inputs (FILE1 or FILE9, and FILE3).
4. **COMBIN duplicates:** the manual does not say what happens when the same frequency is in both FILE81 and FILE82. *Proposal:* error by default, with an option to prefer FILE82 (the newer file).
5. **"Each of the above modes involves only two computer runs"** (§4.2.3) conflicts with the three-run lists in §4.2.4 (HOUSE, ANALYS, STRESS) and §4.2.5 (SITE, HOUSE, ANALYS). *Interpretation:* the core is "regenerate the changed input + ANALYS restart". HOUSE in §4.2.5 is needed only for incoherent analysis, and STRESS in §4.2.4 is post-processing.
6. **What COOXxxx and COOTKxxx contain** is not documented. The inference (COOX = impedance-related per-frequency matrix; COOTK = triangularized total stiffness) comes from "Save Restart Files" keeping only COOXxxx for New Structure. The re-implementation can use its own binary format (for example HDF5/NPZ) under these names.
7. **Restart file numbering:** "xxxxx" (§4.2.5) vs three digits (Ch. 3/6). *Decision:* 3-digit zero-padded frequency **order** number (1-based position in the sorted frequency list, not NFREQ).
8. **FILE9 numbering:** Ch. 3 says "Only one digit load case number can be appended to FILE9", but the names given are FILE9001…FILE9500 (3 digits). *Decision:* 3-digit zero-padded.
9. **Non-uniform input zones:** up to 5,000 zones (this chapter and the HOUSE limit), but the ME command allows `<no>` between 1 and 10 input motions. Clarify how zones relate to motions. It is likely many zones but at most 10 distinct motions, or the ME text is outdated.
10. **Interaction nodes "below the ground surface"** (Item 13b): surface foundations have interaction nodes *at* the ground surface (z = ground elevation, which is the top layer interface). *Decision:* interpret as z ≤ ground elevation, within a tolerance.
11. **Half-space sublayer distribution:** the manual says only "increasing thickness with depth" for the n layers over 1.5λ. The exact SITE rule, and the dashpot coefficients, are not given. *Proposal:* see Item 20 [DECISION]. Verify against analytical half-space solutions.
12. **Half-space depth at very low frequency:** H = 1.5 Vs/f grows without limit as f → 0. At the first frequency (e.g. DF = 0.049 Hz with Vs = 1500 m/s, H ≈ 46 km) this may be numerically awkward. Should there be a cap?
13. **"Zero-frequency ATF"** (Item 4 warning): frequency numbers must be positive integers, so the "zero frequency" is the first (lowest) frequency. The pass/fail tolerance for "close to 1.00" is not given. *Proposal:* warn if the deviation is > 5 %.
14. **Which wave velocity sets the 1/5 rule:** the manual says "wavelength" without naming S or P. *Decision:* use the **strain-compatible Vs** (the governing, shortest wavelength), and optionally check Vp for P-wave input. It is also unclear whether low-strain or iterated Vs should be used; iterated Vs is lower and therefore conservative.
15. **Horizontal mesh ratio:** given only as 1.2–1.5 (sometimes 2) "case by case". The default warning threshold is an implementer choice.
16. **Quarter models for method validation:** this chapter (ASCE/SRP requirement) and Ch. 1.5.1 (ACS SASSI does not recommend quarter models for qualifying SM/MSM) say different things. The UI should show both statements and not enforce either.
17. **Seven incoherent approaches:** this chapter names five deterministic approaches plus "stochastic simulation". The count of seven assumes stochastic with and without phase adjustment, taken from Ch. 6. Confirm against the HOUSE/MOTION spec.
18. **Complex-modulus damping convention** (how `K + iωD` relates to the hysteretic damping ratio β, e.g. `G* = G(1 − 2β² + 2iβ√(1−β²))` in classic SASSI vs `G(1 + 2iβ)`), and the **time-harmonic sign convention** (`e^{+iωt}` vs `e^{−iωt}`). Neither is stated in this chapter; take them from the Ch. 2 / SITE / HOUSE specs. Ch. 1.5.1 says only that damping is applied through complex moduli and is frequency-independent.
19. **Nonlinear SSI iteration convergence criterion** (strain or modulus tolerance, maximum iterations) and the effective-strain ratio used by STRESS are not given here. Take them from the STRESS/HOUSE nonlinear spec.
20. **90 % cumulative modal mass criterion** (Step 3) is cited only through ASCE 4 Section 3. How it combines with the input motion's frequency content to set f_cut is left to the engineer. The tool should only advise.
21. **Option AA wording** in Step 10c is garbled. *Interpretation:* with Option AA, K and M come from ANSYS matrix files, and the ACS SASSI model supplies only the topology for DOF mapping and post-processing.
22. **The "seven smart interpolation schemes"** are defined only by name in Ch. 6. Their window, averaging and smoothing algorithms belong to the MOTION/STRESS spec and need separate clarification.
