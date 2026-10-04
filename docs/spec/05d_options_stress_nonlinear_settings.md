# 05d — Analysis Options: STRESS, RELDISP, NONLINEAR, AFWRITE/CHECK tabs; Options-menu settings (Windows Settings, Colors, Font, Shader, Reset Plot); View and Help submenus

Source: ACS SASSI V3 User Manual, Section 6.5.4 (STRESS, RELDISP, NONLINEAR and AFWRITE tabs of the
*Analysis Options* dialog), 6.5.5–6.5.9, 6.6, 6.7.
Text lines 6986–8043 of `reference/acs-sassi.txt` (printed pages 162–190; Read-tool PDF pages 164–192).
I also checked these related parts of the manual so the spec is complete. They are cited inline and are not
re-specified here:

| Topic | Where in the manual | Used for |
|---|---|---|
| Section 1.5.4 nonlinear soil and structure theory, Figs 1.2/1.3/1.4 | lines 1062–1268, printed pp. 22–28 | NONLINEAR workflow, panel kinematics |
| Section 3 module descriptions (STRESS #9, RELDISP #10, NONLINEAR #12) and Tables 3.1/3.2 | lines 2090–2330, printed pp. 46–52 | file names, frame naming |
| MOTION options (Phase Adjustment, Interpolation, Smoothing, time-history data) | lines 6621–6985 | shared fields |
| HOUSE "Nonlinear SSI Analysis Input" (.pin file) | lines 6191–6245 | nonlinear soil via STRESS |
| Command reference: STRESS, EOUT, RELD, RELFILE, RDND, AOPT, AFWRITE, CHECK, EQL, BBC/BBCI/BBCP/BBCX/BBCY, P, S, B, NONLINBAT, NONLINMOTDISP, SECDATAOPT, READSTR, BINOUT, COMBTHSDB, COMBDISPDIR, COMBDISPDB, THSHLSTR | Section 9 | command equivalents |
| Chapter 10 CHECK errors 49, 50, 67, 73–82, 121–128, Warning 9 | lines 14377+ | validation rules |

Tags: **[CORE]** needed for computational correctness · **[IO]** file/input format · **[UI]** user interface,
plotting or convenience · **[ADV]** advanced (incoherency, nonlinear/Option NON, etc.).
"(inferred)" means the manual does not state the item outright. I derived it from screenshots or from
consistency with other sections, or I propose it as standard mechanics. Everything marked "(inferred)" also
appears in the Open Questions list.

---

## 0. Architecture notes for the implementer

* **[UI]** Each module's options live on one tab of a single modal **Analysis Options** dialog
  (Options → Analysis). The tab order shown in screenshots is:
  `EQUAKE | SOIL | SITE | POINT | HOUSE | FORCE | ANALYS | MOTION | STRESS | RELDISP | NONLINEAR | AFWRITE`.
  Buttons: **Ok**, **Cancel**.
* **[CORE]** Several fields are *shared model variables*. The manual repeatedly says "the value is the same as
  set in ..." another tab. Implement them as one source of truth that is displayed on several tabs:
  * *Type of Analysis* (Seismic / Foundation Vibration) belongs to ANALYS.
  * *Nr. of Fourier Components*, *Time Step of Control Motion* and *Frequency Set Number* belong to SITE.
  * *Multiplication Factor*, *Max. Value for Time History*, *First Record*, *Last Record*, *Title*,
    *Acceleration History File* and *File Contains Pairs* belong to MOTION. STRESS and RELDISP show them too.
* **[IO]** `AFWRITE` writes one input file per enabled module into the model directory, named
  `<modelname><postfix>.<ext>`. This section uses: STRESS → `.str`, RELDISP → `.rdi`, NONLINEAR → `.eql`.
  (The other modules: `.equ`, `.soi`, `.sit`, `.poi`, `.hou`, `.frc`, `.anl`, `.mot`.) Extensions are
  user-configurable (Options → Extension window) and are stored in `SASSIini.xml`.
* **[CORE]** `AFWRITE` always runs `CHECK` first. A module whose CHECK has errors (warnings are fine) does not
  get its input file written.

---

## 1. STRESS Module Options (STRESS tab)

### 1.1 Purpose and data flow [CORE]

STRESS computes **maximum** (absolute-peak) stresses, strains, forces or moments in requested elements. It can
also print and save their **time histories**.

* **Inputs:**
  * `FILE8` holds the complex nodal transfer functions (TFs) at the solved frequencies, from ANALYS. The manual
    calls them "acceleration (displacement) transfer functions".
  * `FILE4` holds element information (connectivity, properties, local axes), from HOUSE.
  * The `.str` input file comes from AFWRITE.
  * The control-motion time-history file.
  * For nonlinear soil it also reads `FILE78` (from HOUSE) and `FILE73` (soil G/γ and D/γ curves, from SOIL).
* **Per requested element, the processing is:**
  1. compute the stress/force/moment component TFs at every solved frequency. These are the **stress transfer
     functions (STF)**;
  2. **interpolate** the STF in complex frequency to the full FFT frequency grid, using the MOTION
     interpolation schemes;
  3. **convolve** with the control motion, i.e. multiply by its Fourier transform;
  4. return to the time domain by **inverse FFT**.
* **Outputs:**
  * the module output listing (maxima);
  * `.TFU` / `.TFI` / `.THS` text files per element component;
  * optional post-processing files (`.sig`, `.tau`, `.bdsig`, `.bdtau`), frames in `\NSTRESS` and `\SOILPRES`,
    and section-cut frames `ESTRESS_<frame>.ess`;
  * the binary database `Modelname_STRESS.bin`;
  * `FILE74` for nonlinear soil iterations.
* **[CORE] WARNING (directional inputs):** STRESS handles **one input direction per run**. To combine
  co-directional effects for nonlinear soil, use the auxiliary DOS batch program **COMB_XYZ_STRAIN**. The
  manual also spells it COMBIN_XYZ_STRAIN. It runs after each SSI iteration, after the X, Y and Z STRESS runs,
  and SRSS-combines the directional effective shear strains from the directional `FILE74` files into a new
  `FILE74`.

### 1.2 Computational algorithm (implementation recipe) [CORE]

Notation: N = *Nr. of Fourier Components* (power of 2), Δt = *Time Step of Control Motion*,
f_m = m/(N·Δt) for m = 0…N/2, f_k = solved frequencies stored in FILE8.

1. **Read the control motion** a_g(t). The file format is either:
   * first line = time step, then one acceleration value per line; or
   * if *File Contains Pairs Time Step–Accel.* is checked, "time, acceleration" pairs on each line.

   Keep records *First Record* … *Last Record*.

   Scale: if *Multiplication Factor* ≠ 0, multiply by it. If *Max. Value for Time History* ≠ 0, scale so that
   max|a_g| equals that value (in g). Exactly one of the two must be non-zero (CHECK errors 77/78).

   Zero-pad to N samples.
2. **FFT:** A_g(f_m) = FFT{a_g}.
3. **STF at the solved frequencies:** for element e, gather the complex nodal DOF transfer functions
   U_e(f_k) from FILE8. Then STF_c(f_k) = S_c · U_e(f_k), where S is the element stress-recovery operator.

   | Element type | Operator S |
   |---|---|
   | SOLID, PLANE, SHELL membrane/bending, TSHELL | D·B evaluated at the element centroid |
   | BEAMS | local stiffness k_L · transformation T, giving end forces |
   | SPRING | spring constants × relative DOF displacements |

   Strains use B alone. **Rigid-body motion produces no strain**, so total TFs can be used directly.
4. **Interpolate** STF_c from f_k to every f_m with the selected *Interpolation Option* (0–6). Apply the
   *Smoothing* parameter and *Phase Adjustment* exactly as in MOTION (see the MOTION spec). Interpolation is
   applied to the **stress TF**, not to the nodal TFs: the manual says it "computes the stress ... components
   at each frequency, performs interpolation".
5. **Convolution:** R_c(f_m) = STF_c(f_m) · G(f_m). G is the Fourier transform of the control motion
   expressed in the kinematic quantity that matches the TF normalisation.

   For seismic analysis the FILE8 TF ratio is the same for total acceleration and total displacement
   (Ü/Ü_g = U/U_g). So stresses need the control **displacement** spectrum:
   G = −g·A_g(f)/(2πf)², with g = acceleration of gravity from the model and the f = 0 term set to 0.
   (inferred, see OQ-S1.)

   For Foundation Vibration, G is the Fourier transform of the load history (inferred, OQ-S2).
6. **Inverse FFT** with Hermitian symmetry gives r_c(t), t = n·Δt.
7. **Derived components** (octahedral shear stress or strain, maximum shear strain) are computed **in the time
   domain** from the component histories, then their peak is taken. STRESS therefore also computes and prints
   the maxima of any components it needs internally (for example all 6 stress components for octahedral
   stress). It saves their histories only if they were explicitly requested.

   The manual says "MOTION" and "FILE12" in this passage. It is the STRESS run, and "not saved" means no
   time-history file is written for the auxiliary components.
8. **Output:**
   * print max|r_c| for requests 1 and 2;
   * save r_c(t) for request 2 (`.THS`), honouring *Skip Time History Steps*;
   * if *Output Transfer Function* is set, save the TF at the solved frequencies (`.TFU`) and the interpolated
     TF (`.TFI`).

### 1.3 STRESS tab — dialog fields

Layout as in the screenshot (left column top to bottom, right column top to bottom). The defaults/example
column shows the values in the manual's screenshot.

| Group box | Field (label in dialog) | Text name in manual | Type | Values / meaning | Example in screenshot | Command arg | Tag |
|---|---|---|---|---|---|---|---|
| Operation Mod(e) | Solution / Data Check | Operation Mode | radio | Solution = 0 (complete solution), Data Check = 1 | Solution | `STRESS,<opmode>` | [CORE] |
| Type of Analysis | Seismic / Foundation Vibration | Type of Analysis | radio | shared with the ANALYS tab (0 = seismic, 1 = foundation vibration) | Seismic | (ANALYS `<type>`) | [CORE] |
| Output Control | Auto Computation of Strains in Soil El. | Auto Computation of Strains in Soil Elements | check | 1 = STRESS automatically computes strains in **all soil elements** for nonlinear soil iterations. Use it together with the TF save options. | off | `STRESS,,<iter>` (1/0) | [ADV] |
| Output Control | Save Stress Time Histories | Save Stress Time Histories | check | 1 = save stress time histories in `.ths` files. Use together with per-element component requests. | on | `STRESS,,,<save>` | [IO] |
| Output Control | Output Transfer Function | Save Transfer Functions | check | 1 = write the TFs of element components. The text says "for beam element nodal forces and moments". Writes `.TFU` (computed) and `.TFI` (interpolated). | off | `STRESS,,,,<itran>` | [IO] |
| Output Control | Phase Adjustment | Phase Adjustment | int | 0 = none, 1 = phase adjustment (defined as in MOTION) | 0 | (shared semantics with MOTION) | [ADV] |
| Output Control | Interpolation Option | Interpolation Option | int 0–6 | see the table in 1.4 | 1 | `STRESS,,,,,<interopt>` | [CORE] |
| Output Control | Smoothing Option | Smoothing Parameter | real ≥ 0 | 0 for coherent analysis. Typically 10–1000 for incoherent analysis. Must be 0 with option 6. | 0 | (inferred to be in `.str`) | [ADV] |
| (text only) | Skip Time History Steps | Skip Time History Steps | int ≥ 0 | number of steps skipped in the output time history (CHECK error 67 if negative). **Not visible in the 2017 screenshot.** | — | — | [IO] |
| Acceleration Time History Data | Nr. of Fourier Components | same | int, power of 2 | shared with SITE. Warning 9: if not a power of 2, the UI writes the nearest power of 2. | 8192 | — | [CORE] |
| " | Time Step of Control Motion | same | real > 0 (s) | shared with SITE (error 49) | 0.005 | — | [CORE] |
| (text only) | Frequency Set Number | same | int | shared with SITE (error 44 if undefined). **Not in the 2017 screenshot.** | — | — | [CORE] |
| " | Multiplication Factor | same | real | shared with MOTION. Use only if Max Value = 0. | 1 | — | [CORE] |
| " | Max Value for Time History | Max. Value for Time History | real (g) | shared with MOTION. Use only if Mult. Factor = 0. | 0 | — | [CORE] |
| " | First Record | same | int ≥ 1 | shared with MOTION (error 74) | 1 | — | [IO] |
| " | Last Record | same | int ≥ First | shared with MOTION (errors 75, 76) | 5000 | — | [IO] |
| " | Title | same | text | shared with MOTION | acc_X_8192 | `THTIT` | [UI] |
| " | File | Acceleration History File | path | shared with MOTION (error 73 if missing) | C:/test/tshell/acc_X_8192.acc | `THFILE` | [IO] |
| " | File Contains Pairs Time Step – Accel. | same | check | shared with MOTION | off | — | [IO] |
| Binary Processing Option | Save Binary Database | Save Binary Database | check | writes `Modelname_STRESS.bin` (see `BINOUT`) | on | `BINOUT,,<str>` | [IO] |
| Element Output Data | list (Group / Element List / Output Code) + **Add**, **Edit**, **Delete** | Element Output Data | list | see 1.5 | `10 | 1-28 | 000000000000` | `EOUT` | [IO] |
| Components | radio set (depends on the element type of the selected list) | Components | radio | see 1.6 | Force NXX (TSHELL) | — | [UI] |
| Component Request | No Request / Print Only Maximum / Print Maximum and Save Time History Response | — | radio | 0/1/2 for the selected component (see 1.7) | No Request | `EOUT <codeN>` | [IO] |
| Post Processing Options | Save Max Value; Save Time History; Restart for Nodal Stress Contours; Restart for Soil Pressure Contours; **Frame Selection** button | — | checks + button | see 1.9 | all off | (inferred to be in `.str`) | [IO]/[UI] |
| Section Cut Options | Save Time History | — | check | writes `ESTRESS_<frame>.ess` files (see 1.10) | off | `SECDATAOPT,<flag>` | [IO] |

**[UI] Label mismatches to keep in mind.** In each pair the dialog label comes first, then the manual text:
"Operation Mod" (truncated) / "Operation Mode"; "Output Transfer Function" / "Save Transfer Functions";
"Smoothing Option" / "Smoothing Parameter". Use the dialog labels on screen and the text names in
documentation and tooltips.

### 1.4 Interpolation Option values (same as MOTION) [CORE]

| Value | Scheme |
|---|---|
| 0 | SASSI2000: dense overlapping windows (weighted averaging) |
| 1 | Original SASSI 1982: non-overlapping windows (no averaging) |
| 2 | Dense overlapping windows (averaging) |
| 3 | Only three overlapping windows (averaging) |
| 4 | Non-overlapping windows with reduced shift (no averaging). MOTION says "one position shift". |
| 5 | Non-overlapping windows with large shift (no averaging). MOTION says "two position shift". |
| 6 | Complex bicubic spline interpolation (no windowing; needs a denser frequency grid) |

* **WARNING [ADV]:** spline (6) is recommended for incoherent SSI because it does not overshoot. First confirm,
  on the coherent analysis, that the frequency set is dense enough that the spline does not smooth or clip
  significant STF peaks. MOTION guidance: more than about 200 frequencies is adequate for incoherent analysis.
* **WARNING [ADV]:** the smoothing parameter must be **0 for coherent analysis** and **non-zero for incoherent
  analysis**, except with Option 6, where it must be 0 (smoothing has no effect with the spline). Run
  sensitivity studies before production runs. The 2007 EPRI study found essentially identical ISRS for
  smoothing values of 10–500.
* **[ADV] Phase Adjustment** (MOTION definition):
  * 0 = keep the physical complex phasing.
  * 1 = reduce the differential phase between Fourier components to nearly zero, depending on the smoothing
    parameter. This gives a conservative "minimum-delay" motion and is used with the stochastic or AS
    incoherent approaches.

  The manual notes there is no literature supporting or rejecting phase adjustment for stresses. With
  stochastic simulation and no phase adjustment, the result is exactly Monte Carlo simulation, whose mean is
  the statistical mean of the incoherent response.

### 1.5 Element Output Data list and the Add/Edit Element List dialog [IO]/[UI]

* **List columns:** `Group | Element List | Output Code`. Example row: `10 | 1-28 | 000000000000`.
* **Buttons:** **Add** adds a list, **Edit** edits the elements of the selected list, **Delete** deletes it.
* **Add/Edit Element List** dialog (title "Add/Edit Element List"):
  * **Element List** text box, for example `1-20`. Separators are blank, tab, `,` or `;`. A range is written
    `a-b`. **Element numbers in a group list must be in ascending order.**
  * **Groups** list box with columns `No | TYPE | Elem. | Title`. Example rows: `1 SOLID 88`, `2 BEAMS 18`,
    `3 BEAMS 32`. "Elem." is the number of elements in the group. The user selects the group here.
  * **Ok** / **Cancel**.
* Each list belongs to exactly **one group**, so it has one element type. Each list stores one request code
  per component of that type.
* **Output Code** (inferred from the screenshot and `EOUT,<code1>,…<code12>,<group>,<element list>`) is a fixed
  string of **12 digits**, one per component slot in the order listed in 1.6. Each digit is 0/1/2 (see 1.7).
  Slots beyond the element type's component count are 0. TSHELL uses 8 slots and BEAMS uses all 12.
* **Command equivalent:** `EOUT,<code1>,…,<code12>,<group>,<element list>` adds an element output request for
  STRESS.
* **[CORE] Nonlinear soil:** if the nonlinear soil SSI option is used, the user must include the
  **nonlinear soil element group** in the output data: all SOLID elements in 3D models, or PLANE elements in
  2D. The nonlinear SSI results are saved at each iteration in **FILE74**.
* **WARNING [CORE] (mixed models):** when the output requests mix BEAMS groups with SOLID and/or SHELL groups,
  check the BEAM results against a separate STRESS run with BEAMS-only requests. There is an *unconfirmed*
  report of different printed maximum axial forces for the first elements of BEAMS groups in mixed runs. A
  re-implementation should produce identical results either way. Use this as a regression test.

### 1.6 Element output components per element type [CORE]

Slot order = Output Code digit order (inferred). Location, axes and units are as stated by the manual.

**SOLID** (3D solid) — 7 components. Computed at the **element centroid**, in **global axes** (X, Y, Z).

| Slot | Dialog label |
|---|---|
| 1 | Stress / Strain XX Direction |
| 2 | Stress / Strain YY Direction |
| 3 | Stress / Strain ZZ Direction |
| 4 | Stress / Strain XY Direction |
| 5 | Stress / Strain XZ Direction |
| 6 | Stress / Strain YZ Direction |
| 7 | Octahedral Shear Stress |

* The figure shows an 8-node hexahedron with nodes 1–8 and centroid 0, and a stress cube with σxx, σyy, σzz,
  τxy, τxz, τyz.
* Octahedral shear stress (standard definition, inferred):
  `τ_oct = (1/3)·sqrt[(σxx−σyy)² + (σyy−σzz)² + (σzz−σxx)² + 6(τxy² + τyz² + τxz²)]`.
  Compute it at each time step from the six component histories. **Requesting it forces all 6 components to
  be computed.**
* For nonlinear soil, the strain measures (HOUSE `.pin` flag) are:
  * flag 0: max component shear strain among γxy, γxz, γyz;
  * flag 1: octahedral shear strain (3D), standard engineering-strain form
    `γ_oct = (2/3)·sqrt[(εxx−εyy)² + (εyy−εzz)² + (εzz−εxx)² + 1.5(γxy² + γyz² + γxz²)]`;
  * maximum shear strain (2D): `γ_max = sqrt[(εxx−εzz)² + γxz²]` (inferred standard forms).

**BEAMS** — 12 components. End forces and moments at nodes I and J, in **local beam axes 1, 2, 3**.

| Slots | Labels |
|---|---|
| 1–6 | Force 1-Direction – Node I, Force 2-Direction – Node I, Force 3-Direction – Node I, Moment 1-Direction – Node I, Moment 2-Direction – Node I, Moment 3-Direction – Node I |
| 7–12 | the same six labels for Node J |

* Figure: axis **1** runs along the member from I to J. Node **K** (orientation node) defines the 1–2 plane,
  with axis **2** pointing toward K's side. Axis **3** completes a right-handed system.
* Recovery: f_L = k_L · T · u_e in local axes (stiffness-based end forces, no member loads).
* File-name component codes from Table 3.1 examples: `FXI`, `MXJ`. So slots map to `FXI FYI FZI MXI MYI MZI
  FXJ FYJ FZJ MXJ MYJ MZJ`, with local 1/2/3 written as X/Y/Z (inferred mapping).

**SHELL** (thin shell / plate) — 6 components. Membrane stresses and plate bending moments at the element
centroid ("infinitesimal element"), in the **local element system x′, y′, z′**.

| Slot | Label | Quantity |
|---|---|---|
| 1 | Force XX-Direction (S x′x′) | membrane normal stress |
| 2 | Force YY-Direction (S y′y′) | membrane normal stress |
| 3 | Force XY-Direction (S x′y′) | membrane shear, S x′y′ = S y′x′ |
| 4 | Moment XX-Direction (M x′x′) | bending moment per unit length |
| 5 | Moment YY-Direction (M y′y′) | bending moment per unit length |
| 6 | Moment ZZ-Direction (M x′y′) | twisting moment, M x′y′ = M y′x′. The label says "ZZ" but it is the twisting moment. |

* Units: the "forces" are **force/length/length (F/L²), i.e. stress**. The moments are **moment/length
  (F·L/L)**.
* The figure shows a quad I-J-K-L and a triangle I-J-K(L), with local z′ normal to the element. The local
  x′/y′ definition belongs to the HOUSE shell element spec.

**TSHELL** (thick plate/shell, the newer element) — 8 components. Membrane forces, transverse shears and
bending moments at the centroid, in **local element axes**.

| Slot | Label |
|---|---|
| 1 | Force XX-Direction (Nxx) |
| 2 | Force YY-Direction (Nyy) |
| 3 | Force XY-Direction (Nxy) |
| 4 | Force XZ-Direction (Qxz) |
| 5 | Force YZ-Direction (Qyz) |
| 6 | Moment XX-Direction (Mxx) |
| 7 | Moment YY-Direction (Myy) |
| 8 | Moment ZZ-Direction (Mxy) |

* Dialog radio labels: `Force NXX, Force NYY, Force NXY, Force QXZ, Force QYZ, Moment MXX, Moment MYY,
  Moment MXY`.
* Units: forces **F/L**. The manual adds "or force over area", which is ambiguous (OQ-S7). Moments **F·L/L**.
* The figure shows a plate with z, w up, θx about y… and resultants Qxz, Qyz, Mxx, Myy, Mxy, Myx on its edges.
* **[ADV] Face stresses (command `THSHLSTR,<flag>`):**
  * flag 0 (default) writes only the 8 basic components.
  * flag 1 also writes top/bottom **face stresses and strains**, built from the maximum component values. For
    the in-plane normal stresses, take the four sign permutations of the maximum force and maximum moment
    contributions: `++`, `−−`, `+−`, `−+` for (Nxx, Mxx) and for (Nyy, Myy). For in-plane shear, use only
    `++` for (Nxy, Mxy), assuming the in-plane shear is the largest. Then compute **principal face stresses
    and strains** from these maximum component stresses.
  * Standard plate formulas, inferred: σ_face = N/t ± 6M/t², τ_face = Nxy/t + 6Mxy/t²,
    σ1,2 = (σx+σy)/2 ± sqrt(((σx−σy)/2)² + τ²). Face strains follow from plane-stress Hooke's law with E, ν.
* `THSHLSMH,<passes>,[type],[workdir]` (not validated in this version) smooths transverse shear forces of
  adjacent elements. type 0 = weighted average with a Parzen window (default), 1 = plain average.
* **No NSTRESS frame and no ESTRESS `.ess` frame generation for TSHELL.**

**PLANE** (2D plane element) — 3 components. At the **element center**, **global axes** (X–Z plane):

| Slot | Label |
|---|---|
| 1 | Sigma XX |
| 2 | Sigma ZZ |
| 3 | Tau XZ |

The figure shows a quad I-J-K-L with center 0, in the global x (horizontal) and z (vertical) axes.

**SPRING** — 6 components, in the global X/Y/Z directions:

| Slot | Label |
|---|---|
| 1 | Force X-Direction |
| 2 | Force Y-Direction |
| 3 | Force Z-Direction |
| 4 | Moment XX-Direction |
| 5 | Moment YY-Direction |
| 6 | Moment ZZ-Direction |

* **WARNING [CORE]:** the sign of a STRESS spring force does **not** indicate tension or compression. Get the
  sign from the relative displacement between the spring end nodes, computed with RELDISP.

### 1.7 Component Request codes [IO]

| Code | Label in text | Label in dialog |
|---|---|---|
| 0 | No Request | No Request |
| 1 | Print Only Maximum Response | Print Only Maximum |
| 2 | Print Maximum and Save Time History | Print Maximum and Save Time History Response |

The UI interaction: select a list, select a component radio button, then select the request radio button. The
request is stored in that component's digit of the list's Output Code.

### 1.8 Nonlinear soil via STRESS (equivalent-linear soil iterations) [ADV]/[CORE]

* For nonlinear soil SSI, STRESS computes the **shear strain time history in the time domain** at each time
  step for every soil element in the nonlinear group.
* Effective strain = effective strain factor × max|γ(t)|. The factor is typically 0.6–0.7 and is set in the
  HOUSE `.pin` file. The effective strain gives new strain-compatible properties through the FILE73 G/γ and
  D/γ curves. Results are written to **FILE74**, which HOUSE reads in the next iteration.
* Loop: `HOUSE → ANALYS (restart, New Structure) → STRESS (X, Y, Z) → COMB_XYZ_STRAIN (SRSS of the X/Y/Z
  FILE74)`, repeated until convergence.
* The strain measure (component / octahedral / max-2D) is chosen in the HOUSE "Nonlinear SSI" input (`.pin`).
* The **Auto Computation of Strains in Soil El.** checkbox (`<iter>` = 1) enables this. The manual says to use
  it together with the save-TF options.

### 1.9 Post Processing Options (STRESS) [IO]/[UI]

These save stresses in **all elements** for contour plotting, and provide **restart** runs that generate frame
files after the response has already been computed.

| Check box | Effect | Output |
|---|---|---|
| Save Max Value(s) | save only the maximum values in all elements | `.sig`, `.tau`, `.bdsig`, `.bdtau` files |
| Save Time History (in All Elements) | save only the time histories in all elements | `.sig`, `.tau`, `.bdsig`, `.bdtau` files |
| Restart for Nodal Stress Contours | compute and save frames for static or animated contour plots | `\NSTRESS\` subdirectory |
| Restart for Soil Pressure Contours | compute and save nodal soil-pressure frames for static or animated contour plots | `\SOILPRES\` subdirectory |
| **Frame Selection** (button) | opens/edits the text file **`Frames.txt`**, which the restart options require | `Frames.txt` |

* "Save" options save responses in all DOFs and can be combined with the Element Output Data requests.
  "Restart" options generate frame files for animation.
* The component suffixes, from Table 3.2:
  * `sig` = normal stress (solids) or membrane normal stress (shells);
  * `tau` = shear stress (solids) or membrane shear (shells);
  * `bdsig` = bending stress (shells only);
  * `bdtau` = bending shear (shells only).
* **Frame naming (Table 3.2) [IO]:**
  * Stress frames: `\NSTRESS\stress_<time>_<fnum>_<comp>`, for example `\NSTRESS\stress_00.000_00001_sig`.
    Maximum frame: `\NSTRESS\stress_ABS_MAX_<comp>`.
  * Soil pressure frames: `\SOILPRES\pres_<time>_<fnum>_<type>` with type `ele` (element values) or `nod`
    (nodal values), for example `\SOILPRES\pres_00.000_00001_nod`. Maximum frame:
    `\SOILPRES\pres_ABS_MAX_<type>`.
  * `<time>` is formatted like `00.000` and `<fnum>` is a 5-digit frame number.
* **Element-type limits [CORE]:**
  * Nodal stress frames handle only **SOLID and SHELL** elements in 3D models.
  * If the model has both SOLID and SHELL elements, the frames contain only average nodal **membrane**
    stresses.
  * If the model has SHELL elements only, separate bending frames (the `bd` files) are also generated.
  * No TSHELL frames.
* **Nodal values [CORE]:** nodal stress = average of the adjacent element-center stresses. No shape-function
  extrapolation is used, and the values are for plotting only. Design values are the element-center values:
  the STRESS output, `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT` and `pres_max_ele`.
* **Soil pressures [CORE]:**
  * They are computed for foundation walls and mat from the **near-field soil SOLID groups adjacent** to the
    structure. These groups must be defined and listed in `Frames.txt`.
  * Element soil pressure = the element-center normal stress on the solid face adjacent to the wall
    (σ_n = nᵀσn; the face/normal rule is inferred, OQ-S10).
  * In addition to the per-time-step frames, one **maximum** pressure frame is written.
  * Also written: `pres_max_ele` (element maxima) and `pres_max_nod` (averaged nodal values, for plotting
    only).
* **Static + seismic soil pressure workflow [IO]:**
  1. The first soil-pressure restart writes **`STATIC_SOIL_PRESSURES.TXT`** filled with zeros, plus
     seismic-only frames. The file has the same layout as `ELEMENT_CENTER_ABS_MAX_STRESSES.TXT` but only one
     value column.
  2. The user combines the three directional seismic pressure frames, in time or by SRSS.
  3. The user enters the non-zero static (geological) pressures.
  4. Rerunning the restart adds static + seismic algebraically, giving total-pressure frames.
* **`ELEMENT_CENTER_ABS_MAX_STRESSES.TXT` [IO]:**
  * It is written only when **all element stresses in all elements** are selected through the
    post-processing check box.
  * It is not usable for BEAMS, whose nodal end forces are needed; it is useful for SHELL, SOLID and SPRING.
  * Format:
    ```
    [# of groups]
    [elem type] [group #] [ordered group #] [# elements in group]      (one line per group)
    ...
    [elem type] [group #] [ordered group #]                            (block header per group)
    [element #] [comp1] [comp2] [comp3] [comp4] [comp5] [comp6]        (one line per element)
    ```
  * Manual example: `3`, then `SOLID 1 1 3 / SHELL 2 1 1 / SOLID 3 1 5`, then blocks.
* **`Frames.txt` format [IO]:**
  ```
  [# of Frames]
  List of Frame Numbers            (one per line)
  [# of Soil Pressure Groups]
  List of Soil Pressure Groups     (all on the same line)
  ```
  The manual's example is `10`, then `1`…`10` on separate lines, then `3`, then `8 9 10`. It saves stress
  frames 1–10 and computes soil-pressure frames for near-field SOLID groups 8, 9 and 10.

### 1.10 Section Cut Options [IO]

* **Save Time History**: STRESS saves frames of **all element-center stresses at all time steps** as
  `estress_<framenumber>.ess` in `.\NSTRESS\`. The manual also writes `ESTRESS_framenumber.ess`.
* The UI then computes section-cut force and moment time histories. Relevant commands:
  * `READSTR,<Filename>,[Dir]` loads an `.ess` file;
  * `CALCSECTHIST` and `CALCSECTHISTDB` compute the histories;
  * `CSECT`, `CUTADD` and related commands define cuts.

  See Demo 8 and Verification Manual Problem 47.
* Command equivalent: `SECDATAOPT,<flag>` (0 = no `.ess` files; 1 = save `.ess` files for the entire time
  history).
* Not available for TSHELL.

### 1.11 Binary Processing Option [IO]

* **Save Binary Database** writes element stress/force component time histories to
  **`Modelname_STRESS.bin`** for fast post-processing. Command: `BINOUT,[mot],[str],[reldisp]` with
  `str` = 1 (write) or 0 (don't). A blank argument keeps the current flag.
* Combine the X, Y, Z directional databases with **`COMBTHSDB,<Xfile>,<Yfile>,<Zfile>,<Comb>`** (full path
  names).
* Related commands:
  * `LOADTHSDB,<file>,[sel]` (sel 0 = ACS SASSI format, the default; 1 = ANSYS format);
  * `THSDBANI` (contour animation);
  * `BINSTRTBL` (CSV table; step −1 = signed absolute maximum);
  * `MAXDBFRAME,THS`;
  * `BINFRAMEOUT,THS,…`.

  These are specified in the binary-database spec.

### 1.12 STRESS output text files [IO]

* Per element component: `etype_gnum_enum_comp.<ext>`. `etype` ∈ {`SOLID`, `BEAMS`, `SHELL`, `TSHELL`,
  `PLANE`, `SPRING`}. `gnum` is 3-digit and `enum` is 5-digit, both zero-padded. `comp` is the component code.
  * Examples: `BEAMS_003_00045_MXJ` (MX moment at node J, beam 45, group 3) and `BEAMS_012_00001_FXI.THS`.
  * Extensions: `.TFU` (uninterpolated, i.e. computed, stress TF), `.TFI` (interpolated stress TF), `.THS`
    (stress time history).
* `FILE74`: nonlinear soil effective strains and properties per iteration.

### 1.13 STRESS command equivalents [IO]

```
STRESS,<opmode>,<iter>,<save>,<itran>,<interopt>
   opmode   0 complete solution | 1 data check only
   iter     1 automatic computation of strains in all soil elements | 0 otherwise
   save     1 save stress time histories in .ths files | 0 do not save
   itran    1 output transfer functions | 0 no TF output
   interopt interpolation option for stress (0..6)
EOUT,<code1>,…,<code12>,<group>,<element list>
SECDATAOPT,<flag>        (0 | 1)
BINOUT,[mot],[str],[reldisp]
THSHLSTR,<flag>          (0 | 1)
THFILE,<file> ; THTIT,<title>     (shared time-history file/title)
RUNSTRESS,[model]        (runs the module)
```

### 1.14 CHECK rules for STRESS (Chapter 10) [CORE]

| Code | Condition |
|---|---|
| Error 79 | No Element Output Request: there are no element lists |
| Error 80 | Illegal Group For Output Request: `<g>` |
| Error 81 | Illegal Element Output Request: `<e>`, Group `<g>` (element not in the group, explicitly or within a range) |
| Error 82 | Element Output Request Defined More Than Once: `<e>`, Group `<g>` |
| Error 49 | Illegal Time Step of Control Motion (negative) |
| Error 50 | Illegal Number of Values for Fourier Transform (negative) |
| Warning 9 | Number of Values for Fourier Transform Is Not Power of 2 (the nearest power of 2 is written) |
| Error 44 / 120 | Frequency Set not defined / empty |
| Error 67 | Illegal Output Time History Step (negative), i.e. Skip Time History Steps |
| Error 73 | Acceleration Time History File Does Not Exist |
| Errors 74 / 75 / 76 | Illegal First Record; Illegal Last Record; First > Last |
| Errors 77 / 78 | Multiplication Factor and Max Value are both zero / both non-zero |

---

## 2. RELDISP Options (RELDISP tab)

### 2.1 Purpose and methodology [CORE]

* RELDISP computes **relative displacement** time histories, and their maxima, at selected nodes and DOFs
  with respect to a **reference location and direction** (reference node/DOF).
* It uses the **interpolated complex acceleration TFs (`.TFI`)** that MOTION produced for the output nodes.
  MOTION must have been run first with complex TFI output at those nodes.
* The options are written to `<model>.rdi`.
* Outputs:
  * the maximum relative displacements;
  * `.TFD` (relative-displacement complex TF) and `.THD` (relative-displacement time history) per node/DOF;
  * optional frames in `\THD` (and `\THDR` for rotations);
  * optional binary databases.
* **Algorithm (inferred standard form; the manual calls it a "refined, analytical approach" based on complex
  ATFs):**
  1. H_rel(f) = H_node,d(f) − H_ref,d(f), using the complex interpolated ATFs on the full FFT grid for the
     same DOF d.
  2. Relative displacement spectrum: D(f) = H_rel(f) · [−g·A_g(f)/(2πf)²], with the f = 0 term set to 0.
  3. d(t) = IFFT{D}.
  4. Report max|d(t)|. Save the `.TFD` (H_rel/(−ω²) or H_rel; OQ-R2) and `.THD` files.
* **One DOF per run:** the reference TFI defines the DOF, so **RELDISP must run six times** to get all six DOF
  relative displacements. Reference TFI file names should carry the DOF tag: `TR_X`, `TR_Y`, `TR_Z`, `R_XX`,
  `R_YY`, `R_ZZ`, as in the naming `nnnnnTR_X.TFI` (characters 1–5 = node number; characters 6–9 =
  translation `TR` or rotation `R` plus the direction).
* **Relative to the free-field input:** the user builds a synthetic reference TFI with amplitude 1.0 in the
  input direction and 0 in the other directions, at every frequency.
  * Coherent input: phase 0 in all directions at all frequencies.
  * Incoherent stochastic simulation: the reference TFI must include the zero-frequency phase.
* Recommended over MOTION's baseline-correction option for relative displacements.

### 2.2 RELDISP tab — dialog fields

| Group box | Field (dialog label) | Type | Meaning | Example in screenshot | Command | Tag |
|---|---|---|---|---|---|---|
| Reference Location and Direction | Complex TF File Name | path | the reference node/DOF `.TFI` file | `C:/test/tshell/00415X.TFI` | `RELFILE,<FileName>` | [IO] |
| Output Control | Save Rel Disp Complex TF | check | write the **complex** (rather than amplitude) relative-displacement TF. The text also says this group selects whether `.TFI` or `.TFU` files are used for the output. | disabled in the screenshot | `RELD,<RelDisOutput>,…` | [IO] |
| Acceleration Time History Data | Nr. of Fourier Components, Time Step of Control Motion, Multiplication Factor, Max Value for Time History, First Record, Last Record, Title, File, File Contains Pairs Time Step – Accel. | — | same as MOTION/STRESS (shared) | 8192, 0.005, 1, 0, 1, 5000, acc_X_8192, C:/test/tshell/acc_X_8192.acc, off | — | [CORE] |
| Nodal Output Data | table `Node Num… | X | Y | Z | XX | YY | ZZ` (✓ marks) + **Add**, **Edit**, **Delete** | list | node/DOF output requests | `1 ✓ ✓ ✓` | `RDND,<NodeNum>,<X>,<Y>,<Z>,<XX>,<YY>,<ZZ>` | [IO] |
| Post Processing Options | Save Relative Displacement in All Nodes | check | compute and save relative displacements, with respect to the reference node/DOF, in **all nodes**, in THD files | off | `RELD,,<RelDispSAll>` | [IO] |
| " | Save Rotations for ANSYS V11.0 | check | same, for the **rotational** DOFs in all nodes (THD files) | off | — | [IO] |
| " | Restart For Frame Generation | check | compute and save relative-displacement frames in `\THD` (and `\THDR` for rotations) for deformed-shape animation | off | — | [UI] |
| Binary Disp. Option | No Binary / TFD Binary / THD Binary | radio | THD Binary writes the relative-displacement history database. TFD Binary is not used in this version (disabled). | THD Binary | `BINOUT,,,<reldisp>` (0 / 1 TFD (unused) / 2 THD) | [IO] |

* **Add/Edit Node** dialog: **Node Number** text box (example `15`), check boxes `X Y Z XX YY ZZ` (example
  X, Z and YY checked), **Ok** / **Cancel**.
* `RDND`: a DOF argument < 1 means the DOF is ignored; ≥ 1 means it is output.
* `RELD,<RelDisOutput>,<RelDispSAll>,<RelDispNumFiles>`:
  * `RelDisOutput` = requested output type;
  * `RelDispSAll` = flag that overrides the node list with all node components and makes RELDISP generate
    displacement frames;
  * `RelDispNumFiles` = number of file names in the output node list.
* THD frame naming: `\THD\THD_<time>_<fnum>`, for example `\THD\THD_00.000_00001`.
* **Binary databases:** the default names for the three response directions are
  `Modelname_TR_X_THD.bin`, `Modelname_TR_Y_THD.bin` and `Modelname_TR_Z_THD.bin`.
  1. For each input direction, combine the principal and coupling response databases into one database
     (`Modelname_THD.bin`) with `COMBDISPDIR,<Xfile>,<Yfile>,<Zfile>,<Comb>`.
  2. Then combine the X-, Y- and Z-input databases with `COMBDISPDB,<Xfile>,<Yfile>,<Zfile>,<Comb>`. The
     manual's text misspells this as "COMDISPDB".
* Supporting commands: `LOADDISPDB`, `DISPDBANI`, `MAXDBFRAME,DISP`, `BINFRAMEOUT,DISP,…`.

---

## 3. NONLINEAR Module Options (NONLINEAR tab) — Option NON [ADV]

### 3.1 Scope

* Option NON consists of the NONLINEAR module plus the COMB_XYZ_THD auxiliary program. It covers nonlinear
  **structure** SSI by iterative **equivalent linearisation** in the complex frequency domain.
* Applications: RC low-rise shear-wall cracking and post-cracking, rubber base isolators (LRB), pile–soil
  interface, and limited foundation sliding.
* Two element choices:
  1. nonlinear **RC wall panels** (macro-shell: groups of coplanar SHELL elements in a vertical plane);
  2. nonlinear **springs** with **translational** DOFs.
* Nonlinear **beams** are shown in the dialog but are **not available** in this version.
* Demos: Demo 9 (RC shear-wall building, nonlinear panels, 0.60 g), Demo 10 (base-isolated building, LRB as
  nonlinear shear springs), Demo 12 (containment pushover calibrated against SNL test NUREG/CR-6783).
* The Verification Manual validates against PERFORM3D.

### 3.2 Iterative workflow (Fig. 1.2) [CORE]

**Initial elastic SSI analysis**

* Step 1: run `SITE, POINT, HOUSE, ANALYS, MOTION, RELDISP`. MOTION/RELDISP output requests must include all
  panel corner nodes and spring end nodes.
* Step 2: run `COMB_XYZ_THD` (to combine X/Y/Z), then `NONLINEAR`. PANEL.NON is absent, so this is an elastic
  run that creates PANEL.NON = 1.
* Step 3: copy and save all output files:

  | From | To |
  |---|---|
  | `Panelxxxx.THD` | `Panelxxxx_elastic.THD` |
  | `Panelxxxx.THS` | `Panelxxxx_elastic.THD` (as printed; probably `_elastic.THS`) |
  | `Panel.FMU` | `Panelxxxx_elastic.FMU` |
  | `Modelname.hou` | `Modelname_elastic.hou` |
  | `Modelname_new.hou` | `Modelname.hou` |
  | `FILE8` | `FILE8_It1` |
  | `PANEL_EQL_MATL_PROP.TXT` | `PANEL_EQL_MATL_PROP_IT1.TXT` |

**Iterative SSI analyses**

* Run `HOUSE, ANALYS (Restart, "New Structure"), MOTION, RELDISP`, then `COMB_XYZ_THD`, then `NONLINEAR`.
* Then copy and save the per-iteration files: `Panelxxxx_It#.THD`, `Panelxxxx_It#.THS`, `Panelxxxx_It#.crv`,
  `Panel_It#.FMU`, `Panel_EQL_Matl_Prop_It#.txt`, `Modelname_It#.hou`.
* Repeat until **CONVERGENCE**. The convergence criterion is unspecified (OQ-N5).

**Supporting tools**

* `NONLINBAT,<Sel>` generates a generic batch file for the whole analysis, using the active model and the
  module locations: Sel 0 = single direction, 1 = three directions, which adds COMB_XYZ_THD. The user writes
  `COMB_XYZ_THD.inp` (an example ships on the install DVD).
* The panel strings in file names are replaced by **SPRING** for nonlinear springs (`SPRINGxxxx.*`,
  `SPRING.NON`, …).

The module's computational steps:

1. Initial linear elastic (uncracked) SSI analysis with the elastic properties.
2. Time-domain local behaviour of each nonlinear element, from its local relative displacements. This
   calibrates the linearised hysteretic model (secant stiffness plus damping) in complex frequency.
3. A new SSI iteration by fast restart (ANALYS "New Structure") with the updated properties.
4. Check convergence, then stop or continue.

### 3.3 Files [IO]

* **Inputs:** `modelname.sit`, `modelname.hou`, `modelname.eql` (written by the UI) and the nodal
  displacement histories `xxxxxTR_[X,Y,Z].THD`, `xxxxxR_[XX,YY,ZZ].THD` for the panel and spring nodes, from
  RELDISP and combined by COMB_XYZ_THD.
* **Outputs (panels):**
  * `Panelxxxx.crv`: equivalent-linear **elastic modulus–damping curves vs. panel strain amplitude**, after
    rigid-body removal. Computed during the initial run.
  * `Panelxxxx.thd`: panel equivalent-linear strain history.
  * `Panelxxxx_AXIAL.thd`: uniform vertical axial strain history.
  * `Panelxxxx.ths`: nonlinear force history. Together with `.thd` it allows plotting the hysteresis loops.
  * `Panel_EQL_Matl_Prop.txt`: new material properties of the walls.
  * `Panel.fmu`: ductility μ and force-reduction factor Fμ per wall.
  * A general output text file named by the user, containing nearly all input and output.
  * The updated HOUSE input (`Modelname_new.hou`) for the next iteration.
* **Outputs (springs):** `SPRINGxxxx.crv`, `SPRINGxxxx.thd`, `SPRINGxxxx.ths`, `SPRING_EQL_Matl_Prop.txt`,
  `Spring.fmu`.
* **State files [CORE]:**
  * `PANEL.NON` / `SPRING.NON`:
    * absent → NONLINEAR runs a **linear elastic** analysis and creates the `.NON` file containing `1`;
    * present and containing `1` → NONLINEAR uses the previous iteration's SSI results and runs a nonlinear
      iteration.
  * `PANEL_EQL_MATL_PROP.txt` / `SPRING_EQL_MATL_PROP.txt`: if present, the next iteration uses the
    equivalent-linear stiffness and damping in this file.
  * **To restart from the elastic analysis, delete both the `.NON` file and the `*_EQL_MATL_PROP.txt`
    file.**

### 3.4 Kinematics of nonlinear elements [CORE]

**Wall panels** (Figs 1.4 and 1.3 references):

* A panel is a set of **coplanar shell elements in a vertical plane**.
  * It should be rectangular, though it may contain triangular shells.
  * Its plane may have any orientation in the horizontal plane.
  * It is assumed to deform uniformly in shear or in bending.
* Only the **four corner nodes** are used. `NONLINMOTDISP` finds them with a node-connection-counting
  algorithm and adds them to the MOTION and RELDISP output requests. **Run it before AFWRITE.** Warped,
  non-4-cornered, triangular-meshed or concave panels can make it, or NONLINEAR, fail.
* Local panel axes (inferred):
  * h = unit horizontal vector along the bottom edge;
  * v = global Z;
  * corner displacements are projected onto h and v.
* Remove the panel's rigid-body motion first. Then:
  * **Shear strain** (Disp. Opt 1) = the **average rotation of the two vertical (lateral) edges**, i.e. the
    average relative story drift of the two lateral edges / height: γ = (γ1 + γ2)/2 (Fig. 1.4 left).
  * **Bending rotational strain / curvature** (Disp. Opt 2) = the **relative rotation of the top and bottom
    (horizontal) edges**, θ1 − θ2 (Fig. 1.4 right). It is divided by height when expressed as curvature
    (units "radian/length").
  * **Uniform vertical axial strain** is also computed (`Panelxxxx_AXIAL.thd`), but no axial hysteretic model
    exists.
* Proposed rigid-body-invariant formulas (inferred, OQ-N2), with corners BL, BR, TL, TR, width L and height H:
  ```
  γ   = ½[(u_h,TL − u_h,BL) + (u_h,TR − u_h,BR)]/H  +  ½[(u_v,BR − u_v,BL) + (u_v,TR − u_v,TL)]/L
  θ_b = (u_v,BR − u_v,BL)/L ;  θ_t = (u_v,TR − u_v,TL)/L ;  κ = (θ_t − θ_b)/H
  ε_v = ½[(u_v,TL − u_v,BL) + (u_v,TR − u_v,BR)]/H
  ```
* **Curved walls** (for example a containment cylinder): refine the panel mesh so that the angle between
  adjacent panels is **< 15°**. One shell per panel is allowed; `UNIPNL` creates a unique group per element.

**Nonlinear springs:**

* Use the relative displacement between the spring end nodes along the active DOF. **Disp. Type / Dof.:**
  1 = X, 2 = Y, 3 = Z translation.
* Nonlinear springs with up to three translational stiffnesses should be split into groups of **1D**
  nonlinear springs, so that each DOF gets its own stiffness and damping.
* The **GMR** hysteretic model applies under shear or axial deformation.

### 3.5 Equivalent-linearisation algorithm [CORE]

For each nonlinear element (panel or spring), at each iteration:

1. Get x(t): the panel strain (γ or κ) or the spring relative displacement, from the combined 3-direction
   THD histories.
2. **Equivalent-linear amplitude:** x_eq = EDF · max_t|x(t)|. **EDF = Disp. Factor**, typically 0.7–0.9 with a
   best fit near **0.80**. EDF = 1.0 makes elements too soft; 0.6 makes them too stiff.
3. Run the selected **hysteretic model** (Force Opt) on x(t) in the time domain to get the nonlinear force
   history F(t). Write `.thd` (x) and `.ths` (F).
4. **Secant stiffness:** K_sec(x_eq) from the BBC (and/or the hysteresis loop). Normalise it by the elastic
   stiffness K_el and scale the elastic modulus: **E_new = E_elastic · K_sec/K_el**.
   * Shear BBC: K_el is in G·A_shear units.
   * Bending BBC: K_el is in E·I units.
   * Spring: K_el is the elastic spring constant.

   Under the isotropy assumption, one modulus reduction degrades the **shear, axial and bending** stiffness
   equally, with **Poisson's ratio constant**. This is reasonable only for shear-governed low-rise walls.
   ASCE 4-17 reduces only shear and bending, which is why the vertical direction is suggested to stay
   elastic.
5. **Equivalent damping:** hysteretic damping ξ_h(x_eq), interpolated from the `.crv` E–D curves.
   * Then apply *Damping Scale Factor*.
   * If **Include Elastic Damping** is set, ξ = ξ_h + ξ_elastic (the initial viscous damping). Otherwise
     ξ = ξ_h, which should not be smaller than the initial elastic damping.
   * Cap ξ at *Damping Cutoff %* when that is non-zero.
   * The order of these operations is unspecified (OQ-N6).
6. Write the new E and damping into `Modelname_new.hou` and `*_EQL_Matl_Prop.txt`.
7. **Ductility** μ = max|x_nonlinear| / x_cracking, where x_cracking is the displacement where the elastic
   BBC region ends, i.e. BBC point 1 (not the yield point). **Force-reduction (inelastic absorption) factor**
   Fμ = initial linear-elastic force / final nonlinear force. Write both to `Panel.fmu`.

   These factors correspond to the ASCE 43-05 inelastic absorption factors for shear-wall limit states.

### 3.6 Hysteretic models [ADV]

| Code (BBC Type / Force Opt) | Model | Used for | Status in this version |
|---|---|---|---|
| 1 | Cheng-Mertz Shear (CMS) | panel shear deformation. Best for low-rise walls (height < width). Pinched loops when the cracking point is low. | available for panels |
| 2 | Cheng-Mertz Bending (CMB) | panel bending deformation | **not included** |
| 3 | Takeda (TAK) | bending and shear. Best for bending-governed panels. | listed. **WARNING:** only usable with equivalent damping = the constant elastic damping. The yield point force/moment is the panel's **peak capacity** in Takeda, so a lower yield point gives a lower capacity. |
| 4 | General Massing Rule (GMR) | nonlinear springs (shear or axial) | springs only |

* The manual gives no constitutive equations for these models. Fig. 1.3 compares the CMS, CMB and TAK loops
  for the same strain history. Implementers must use the published formulations: Cheng & Mertz (1989), Takeda
  et al. (1970), and the Masing rule (unloading/reloading branch F − F_r = 2·f((x − x_r)/2)). See OQ-N3.
* **CHECK in this UI version:**
  * Error 126: panel Force Option must be **1**;
  * Error 127: spring Force Option must be **4**;
  * Error 128: panel Displacement Option must be **1**.

### 3.7 Backbone curve (BBC) conventions [CORE]

* The BBC is given by the user as force–displacement, shear force–shear strain, or moment–rotation/curvature
  points.
* **The origin (0,0) is omitted.** **Point 1 = cracking point** (cracking displacement and force).
* **Yield Num.** is the 1-based index of the yield point within the BBC points, and it must be defined.
* Units are user-consistent and never converted.
  * Disp. Opt 1 (shear strain) needs a BBC in *shear strain (decimal) – force*.
  * Disp. Opt 2 needs *section rotation/curvature (rad/length) – moment*.
  * The BBC must be consistent with the elastic stiffness (G·A_shear for shear, E·I for bending).
* Generation commands:
  * `BBCGEN,<Panel>,<ShearModel>,[fc],[fy],[Pn],[Nu],[bre],[bys],[CrackingForceLevel]` always makes a
    **22-point** BBC: point 1 = cracking, then 20 points equally spaced in strain up to yield, then failure
    at (strain 2 %, 1.02 × ultimate shear).
    * ShearModel: 1 ACI 318-08, 2 Wood 1990, 3 Barda 1977, 4 Gulec-Whittaker 2009.
    * CrackingForceLevel: 0 = the ASCE 4-17 C.3.3.2 cracking stress 3√f′c; otherwise 0.10–0.50 × ultimate.
  * `SHEAR,<panel>,[fc],[fy],[P],[Nu],[Fvw],[Fbe]` estimates ultimate shear.
* Pushover or test data, if available, can be used to adjust the elementary panel BBCs (Demo 12).

### 3.8 NONLINEAR tab — dialog fields

**Global Modeling Options** — command `EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>`. This
sets the `.eql` header.

| Field (dialog) | Text name | Type | Meaning / recommended | Example in screenshot | EQL arg | Tag |
|---|---|---|---|---|---|---|
| Disp. Factor | Equivalent-Linear Displacement Factor (EDF) | real | 0.7–0.9, best ≈ 0.80. Gives x_eq = EDF·max. | 0.8 | `disp` | [CORE] |
| Damping Cutoff % | Damping Cutoff | real (%) | user cut-off on the damping value, used to match ASCE 4 / USNRC limits (see the note below) | 0 | `dampCutoff` | [CORE] |
| Damping Scale Factor | Damping Scale Factor | real | scales the damping to calibrate against test data | 0 | `dampScale` | [CORE] |
| Material Parameter | Material Parameter | — | **disabled in this version** | (greyed) | — | [UI] |
| Use Non-linear Panels | Use Nonlinear (Element) Options | check | include the nonlinear wall panels | off | part of `NonLinOpts` | [CORE] |
| Use Non-linear Springs | " | check | include the nonlinear springs | on | part of `NonLinOpts` | [CORE] |
| Use Non-linear Beams | " | check | **disabled** (beams not available) | greyed | part of `NonLinOpts` | [UI] |
| Include Elastic Damping | Include Elastic Damping | check | ξ_total = ξ_hysteretic + ξ_elastic. If unchecked, ξ = ξ_hysteretic only (which should be ≥ ξ_elastic). | off | `ElasicD` (0 = don't include, 1 = include) | [CORE] |

Damping Cutoff note. The manual states the ASCE 4-17 Section C.3.3.2 response levels as:

| Level | Stiffness | Damping |
|---|---|---|
| 1 | uncracked | 4 % |
| 2 | 0.50 × uncracked shear/bending | 7 % |
| 3 | 0.50 × uncracked shear/bending | 7 % |

The manual gives Levels 2 and 3 identically; this is probably a typo (OQ-N9).

**Backbone Curve Data:**

| Field | Meaning | Example | Command | Tag |
|---|---|---|---|---|
| Backbone Curve (spin box) | number of the displayed BBC | 1 | — | [UI] |
| Type | BBC hysteretic type, used **only for the model titles in the printout**: 1 CMS, 2 CMB, 3 TAK, 4 GMR (seen in the screenshot) | 4 | `BBCI,<num>,<yield>,<type>` / `BBC,<num>,<type>,…` | [IO] |
| Yield Num. | index of the yield point among the BBC points | 11 | `BBCX/BBCY … <yield>` | [CORE] |
| X / Y grid (rows 1…n) | BBC points. **View/edit only**: the grid cannot add or remove points. Curve size and points must first be created with the BBC commands. | rows: (0.01, 100), (0.0223, 220), (0.0232, 226), (0.0244, 232), (0.0265, 238), (0.0302, 244), (0.0374, 251), … | `BBCX,<num>,<points>,<yield>,<X1>…<Xn>`; `BBCY,<num>,<points>,<yield>,<Y1>…<Yn>`; `BBCP,<num>,<point>,<X>,<Y>`; `BBC,<num>,<type>,<points>,<yield>,<file>` (file of X Y pairs) | [IO] |

**Panel Data** (wall panels only). Command `P,<num>,<group>,<bbc>,<disp>,<force>`. See also `PANELGEN` /
`PNLGEN`.

| Field (dialog) | Text name | Meaning | Example |
|---|---|---|---|
| Panel (spin) | Panel | wall panel number in the nonlinear FE model | 1 |
| Group Num. | Group Num. | the SHELL group that forms the panel (all shells coplanar) | 0 |
| BBC Num. | BBC Num. | BBC assigned to the panel | 0 |
| Disp Type | Disp. Opt. | 1 = shear strain (average relative story drift of the two lateral edges after rigid-body removal); 2 = bending section rotation (relative rotation of the top and bottom edges) | 0 |
| Force Opt | Force Opt. | hysteretic model, same codes as the BBC type: 1 CMS, 2 CMB, 3 TAK | 0 |

**Spring Data** (nonlinear springs). Command `S,<num>,<group>,<elem>,<bbc>,<disp>,<force>`. The springs must
already exist (E command).

| Field (dialog) | Text name | Meaning | Example |
|---|---|---|---|
| Spring (spin) | Spring | nonlinear spring number | 1 |
| Group Num. | Group Num. | spring group number | 6 |
| Elem Num. | Elem. Num. | element number in the group | 1 |
| BBC Num. | BBC Num. | BBC used | 1 |
| Dof. | Disp. Type | 1 X, 2 Y, 3 Z translation | 1 |
| Force Opt | Force Opt. | 1/2/3 not applicable (CMS/CMB/TAK); **4 = GMR** | 4 |

**Beam Data** (**NOT AVAILABLE IN THIS VERSION**; all fields greyed). Command
`B,<num>,<group>,<spgroup>,<bbc>,<force>,<end1>,<end2>`.

| Field | Meaning | Example |
|---|---|---|
| Beam (spin) | beam number | 1 |
| Group Num. | nonlinear beam group | 0 |
| Spring Gr. | spring group at the beam ends | 0 |
| BBC Num | BBC | 0 |
| Force Opt | 1–3 not applicable; 4 = General Massing Rule | 0 |
| Beam End 1 | element defining the node at one end | 0 |
| Beam End 2 | element defining the node at the other end | 0 |

**Related commands [UI]:**

* Panelisation: `MERGEPANEL`, `WALLFLR`, `PANELIZE`, `MERGEGROUP`, `SPLITGROUP`, `EDGEMODEL`, `EDGE`, `UNIPNL`,
  `PANELGEN`/`PNLGEN`.
* Data: `EQL`, `P`, `S`, `BBCX`, `BBCY`, `SHEAR`, `BBCGEN`.
* Lists and deletion: `PLIST`, `PDEL`, `DELBBC`, `DELSPR`, `DELBM`.
* For models with many panels, prefer commands over the dialog (Demo 9).

### 3.9 CHECK rules (NONLINEAR) [CORE]

| Code | Condition |
|---|---|
| Error 121 | No Panels Specified |
| Error 122 | Material `<i>` referenced in panel `<k>` does not exist |
| Error 123 | Group `<i>` referenced in panel `<k>` does not exist |
| Error 126 | Force Option `<i>` not supported for panel `<k>` (set 1) |
| Error 127 | Force Option `<i>` not supported for spring `<k>` (set 4) |
| Error 128 | Displacement Option `<i>` not supported for panel `<k>` (set 1) |

---

## 4. AFWRITE and CHECK Options (AFWRITE tab) [IO]/[UI]

* One checkbox per module enables or disables it for `AFWRITE` and `CHECK`. Screenshot states, in this order:

  | Module | State |
  |---|---|
  | EQUAKE | ✓ |
  | SOIL | ✓ |
  | LIQUEF | disabled |
  | SITE | ✓ |
  | POINT | ✓ |
  | HOUSE | ✓ |
  | PINT | disabled |
  | FORCE | ☐ |
  | ANALYS | ✓ |
  | COMBIN | ☐ |
  | MOTION | ✓ |
  | STRESS | ✓ |
  | RELDISP | ✓ |
  | NONLINEAR | ☐ |

* Generated files go in the active model's directory, named `<modelname><postfix>.<ext>` with the postfix and
  extension set in the UI.
* Command: `AOPT,<EQUAKE>,<SOIL>,<DEP>,<SITE>,<POINT>,<HOUSE>,<DEP>,<FORCE>,<ANALYS>,<COMBIN>,<MOTION>,<STRESS>,<RELDISP>,<PANEL>`.
  * Each argument is 1 = include, 0 = skip.
  * The two `<DEP>` slots are modules not in this version and must be 0. They correspond to the greyed
    LIQUEF and PINT boxes.
  * `<PANEL>` is the NONLINEAR flag; the manual's text says "RELDISP", a typo.
* `AFWRITE` runs `CHECK` first. A module with errors (warnings excluded) is not written. Modules unchecked here
  are ignored by both.

---

## 5. Options menu: 6.5.5 – 6.5.9 [UI]

### 5.1 Windows Settings (6.5.5)

Opens the window-settings dialog for the **active plot**. The dialogs are specified in the Plotting chapter
(Section 7, "Windows Settings" / "Windows Option"). For 3D plots, hide requests for groups, elements and nodes
are stored in the model data.

### 5.2 Colors (6.5.6)

* Dialog **"Select Colors"**. It is a tabbed palette for all plot and text colors in one place; older versions
  scattered them across windows and allowed only the active window.
* Click the colour swatch button next to a name to open the standard colour picker. On selection the picker
  closes and the swatch updates. **Ok** / **Cancel**.

| Tab | Controls |
|---|---|
| **UI** | Command-history text colours by message type. Rows and defaults as in the screenshot: Command Echo Color (black), Conformation Color (blue), Error Color (red), Comment Color (green), Info Color (purple), Warning Color (yellow/amber). |
| **Element** | element-plot colours other than the element fill colours |
| **ElemPalette** | element fill colours by group, material or beam/spring property. Filled as elements are plotted. **Up to 128 colours**; with more than 128 items, colour index = item number mod 128 (`Elemnum ≡ Colnum (mod 128)`). One shared palette serves groups, materials and properties. |
| **Node** | node plot |
| **Cut** | cut plot |
| **Bubble** | bubble plot |
| **Vector** | vector plot |
| **Contour** | contour plot |
| **Deformed** | deformed-shape plot |
| **Spec** | spectrum plot axis colours (line colours are in SpecLines) |
| **SpecLines** | spectrum line colours, filled as line numbers are plotted |
| **TimeHist** | time-history plot axis colours (line colours are in THLines) |
| **THLines** | time-history line colours, filled as line numbers are plotted |
| **SoilLayer** | Soil Layers plot |
| **SoilProp** | Soil Properties plot |

### 5.3 Font (6.5.7)

* Dialog **"Font Selection"**. Each font-description button opens the platform's font dialog. **OK** applies
  the fonts to the UI, including plots already open; **Cancel** discards.

| Group | Button | Default shown |
|---|---|---|
| Command History | Command History Font | Courier New, 9 |
| 3D Plots | Title Font | Courier New, 9 |
| 3D Plots | Notation Font | Courier New, 9 |
| 2D Plots | Title Font | Courier New, 9 |
| 2D Plots | Notation Font | Courier New, 9 |
| 2D Plots | Legend Font | Courier New, 9 |

* The manual also says fonts apply to the command entry window (OQ-U2).

### 5.4 Shader Options (6.5.8)

Opens the Shader Options window. There is no screenshot, and defaults and ranges are unspecified. The window
sets:

* the maximum point size in node and bubble plots;
* the outline width of solids (the text is truncated to "soli");
* the amount of shrink in element plots;
* the scale factor for animations.

### 5.5 Reset Plot (6.5.9)

Resets the 3D plot view to its default view and zoom.

---

## 6. View submenu (6.6) [UI]

The View menu items are: **Check Errors**, **Command Window**, **Command Display ▸**, **Toolbars ▸**.

### 6.1 Check Errors (6.6.1)

* Opens the error listing window, titled `CHECK: Errors and Warning for - <model>`, for the active model. It
  shows the contents of the `.err` file produced by CHECK.
* If the model has no name or path, or has not been checked, the window is blank or a "file not found"
  message box appears.
* `.err` layout (from the screenshot):
  * a header line per module: `Errors and Warnings for <MODULE>`;
  * then lines `Error <n> : <text>` / `Warning <n> : <text>`.
  * Example: `Errors and Warnings for EQUAKE` / `Error 85 : RS Input File 1 Does Not Exist` / … /
    `Error 73 : Acceleration Time History File Does Not Exist`.
* The Options → Check settings (6.5.3) control the `.err` content: show or hide errors and warnings, and a
  per-module limit (break number) on how many messages of each type are printed. Totals are still counted.
  These settings are not saved and reset to defaults when the UI opens.

### 6.2 Command Window (6.6.2)

* Opens the command window, or makes its tab active (tab "Command History" above a "Command Entry" box).
* Reopening the command history **clears its output**. This does not affect models in memory.

### 6.3 Command Display (6.6.3)

* A submenu of check items that toggle echo of message types to the command history: **Command Echo**,
  **Output Conformation** [sic], **Comments**, **Warnings & Errors**.
* There are 6 message types: echo, confirmation, comment, error, warning, info. Five can be turned off;
  **Info** (default colour purple) cannot.
* All items are checked by default for a new user account. Clicking toggles. Selections persist in
  **`SASSIini.xml`**.
* While reading a `.pre` file, command echo and command acceptance (confirmation) are **always suppressed**,
  regardless of these settings. This greatly speeds parsing. Error and comment output during `.pre` input
  still follows the settings.
* Purpose: copy and paste clean command text from the history into a custom `.pre`.

### 6.4 Toolbars (6.6.4)

* A submenu with **Main Toolbar** and **Plot Toolbar** check items. Checked = visible; both are shown by
  default.
* **Not saved** to `SASSIini.xml`, so they return to the defaults on each start.

---

## 7. Help submenu (6.7) [UI]

* **Help** (shortcut **F1**) opens the online help Overview page. Help moved to a wiki format to avoid
  proprietary help formats.
* **About** shows the UI version number, release date and other information.
* The main window title in the screenshots is "SSI Submodeler". Menus: `Model, File, Plot, Modules, Options,
  View, Help`.

---

## 8. Verification hooks (suggested tests) [CORE]

1. **Octahedral stress:** pure shear τxy = τ gives τ_oct = (√6/3)τ ≈ 0.8165τ. Uniaxial σ gives
   τ_oct = (√2/3)σ.
2. **Beam end forces:** a cantilever beam with a quasi-static, very low-frequency input. The STF of
   M3 at node I should equal k_L·u. Check equilibrium: F_I + F_J = 0 for unloaded members.
3. **Rigid-body check:** a uniform rigid translation TF in all element nodes must give zero STF for solids,
   shells and planes.
4. **Interpolation path:** with option 6 and smoothing 0 on a dense grid, the `.TFI` must reproduce the `.TFU`
   at the solved frequencies.
5. **RELDISP:** reference TFI = the node's own TFI gives zero. With the free-field unit-amplitude reference,
   the result must equal the IFFT of (H − 1)·U_g. Compare with double integration of the acceleration
   difference.
6. **TSHELL face stress:** with N = 0 and M ≠ 0, face stress = ±6M/t².
7. **GMR / Masing with an elastic–perfectly-plastic BBC:** cycling at amplitude x > x_y should give a closed
   loop with ξ_h = 2(x − x_y)/(πx) and K_sec = F_y/x.
8. **NONLINEAR state machine:**
   * no `.NON` file → elastic run and `.NON = 1` is created;
   * `.NON = 1` → iteration;
   * `*_EQL_MATL_PROP.txt` present → its values are used.
9. **`Frames.txt` parser:** the manual example gives frames [1..10] and soil groups [8, 9, 10].
10. **ElemPalette:** item 129 maps to colour 1, or to colour 0 with 0-based indexing (OQ-U1).
11. **CHECK:** duplicate elements across lists raise Error 82; MF = 0 together with Max = 0 raises Error 77.

---

## 9. Open questions / ambiguities

**STRESS**

* **OQ-S1** Whether control-motion scaling to displacement (−g·A_g/ω²) is applied to the STF, or whether FILE8
  holds displacement TFs per unit control acceleration. The manual says only "acceleration (displacement)
  transfer functions". The units of g also need deciding. Recommend: the FILE8 TF is dimensionless and STRESS
  multiplies by the control displacement spectrum.
* **OQ-S2** For Foundation Vibration: which time history is convolved (force history), and the TF units.
* **OQ-S3** Exact semantics of *Skip Time History Steps*: output every (k+1)-th step, or skip the first k
  steps? It is not in the 2017 dialog. Error 67 only requires ≥ 0.
* **OQ-S4** The 12-digit Output Code mapping (digit order = component order, values 0/1/2) is inferred from
  the screenshot and EOUT. The manual does not state it.
* **OQ-S5** SOLID "Stress / Strain" components: does one request output both stress and strain, or is there a
  separate switch? Is octahedral *strain* also output?
* **OQ-S6** Component codes in file names (`etype_gnum_enum_comp`) are given only for BEAMS (FXI, MXJ). The
  codes for SOLID, SHELL, TSHELL, PLANE and SPRING are unspecified. A 6-digit node/element padding may be
  needed above 99,999.
* **OQ-S7** TSHELL force units: the manual says "F/L or force over area". SHELL forces are stresses (F/L²).
  This must be defined per element.
* **OQ-S8** Sign conventions for beam end forces (forces on the element versus on the node) and for springs.
  The manual warns that spring signs do not indicate tension or compression.
* **OQ-S9** Content and format of `.sig`, `.tau`, `.bdsig`, `.bdtau`, the NSTRESS/SOILPRES frame files, `.ess`,
  and `Modelname_STRESS.bin` are not given in this section.
* **OQ-S10** Soil pressure definition: which solid face is used (the one shared with the structure), the sign
  (compression positive?), and how pres_max_nod averages.
* **OQ-S11** Whether a frame number in `Frames.txt` is a time-step index (1-based) or a frame index after
  decimation.
* **OQ-S12** The `.str` file layout: the STRESS command gives only 5 fields. Where phase adjustment,
  smoothing, post-processing flags, the frame option and binary flags are stored is unspecified.
* **OQ-S13** The text says *Save Transfer Functions* applies to "beam element nodal forces and moments". Does
  it also apply to other element types?
* **OQ-S14** How *Phase Adjustment* is implemented for the STF. It is defined only in MOTION, and the manual
  notes there is no literature support for its use on stresses.

**RELDISP**

* **OQ-R1** The `.TFI` and `.TFD` file formats (frequency and complex value columns) are specified in the
  MOTION spec, not here.
* **OQ-R2** Whether `.TFD` stores the relative *displacement* TF (H_rel/(−ω²)·g) or the relative ATF.
* **OQ-R3** The "Output Control" semantics: the dialog shows only *Save Rel Disp Complex TF*, but the text says
  it also chooses TFI or TFU input. The meaning of `RELD,<RelDisOutput>` values is unknown.
* **OQ-R4** The incoherent "zero-frequency phase" in the reference TFI is not defined.

**NONLINEAR**

* **OQ-N1** The `.eql` file format is not documented. It must be designed from EQL, BBC*, P, S and B data.
* **OQ-N2** The exact rigid-body removal and the panel shear, bending and axial strain formulas. Fig. 1.4 is
  only qualitative; the rigid-body-invariant formulas above are proposed. Treatment of the rotational DOFs in
  the THD files is unknown.
* **OQ-N3** The constitutive rules of the CMS, CMB, TAK and "General" Massing models are not given. The
  implementer must choose published versions: Cheng & Mertz 1989, Takeda 1970 (unloading-stiffness exponent),
  and the Masing factor 2 or a generalised factor.
* **OQ-N4** How `.crv` curves are built (cyclic runs at which amplitudes, or from the actual history), and
  whether K_sec comes from the BBC or from the loop peak.
* **OQ-N5** No convergence criterion or tolerance is given (on E, on x_eq, and the maximum number of
  iterations).
* **OQ-N6** The order of damping scale, elastic-damping addition and cutoff. The meaning of value 0 for
  *Damping Scale Factor*: 0 is shown in the screenshot and probably means "not applied", i.e. 1.0. The
  meaning of *Damping Cutoff %* = 0 is probably "no cutoff".
* **OQ-N7** K_el used for normalisation: the slope to BBC point 1 (Y1/X1), or G·A_shear / E·I computed from
  model properties?
* **OQ-N8** The definitions of the "initial linear elastic force" and "final nonlinear force" in Fμ (peak
  values from which iterations?).
* **OQ-N9** The ASCE 4-17 damping levels as quoted (Level 2 and Level 3 are identical in the manual).
* **OQ-N10** The encoding of `EQL <NonLinOpts>` (bit flags for panels, springs and beams?).
* **OQ-N11** The BBC Type field accepts 4 (GMR) in the screenshot, but the text lists only 1–3.
* **OQ-N12** CHECK errors 126 and 128 force panel Force Opt = 1 and Disp Opt = 1, which contradicts the
  documented TAK (3) and bending (2) options. Decide whether to allow them as [ADV] in the re-implementation.
* **OQ-N13** How updated properties map onto materials or spring properties in `.hou`. Each panel or spring
  presumably needs its own material or property (UNIPNL helps).
* **OQ-N14** The file rename table in Fig. 1.2 has an apparent typo (`Panelxxxx.THS → Panelxxxx_elastic.THD`).
  Also, is the renaming done by NONLINEAR or by the batch file?
* **OQ-N15** The `COMB_XYZ_THD.inp` format is not described (it is on the install DVD).

**UI**

* **OQ-U1** The base of ElemPalette indexing (0 or 1) for the mod-128 rule.
* **OQ-U2** Whether colours, fonts and shader options persist in `SASSIini.xml`. Only Command Display is
  stated to persist, and toolbars are explicitly not saved. The Font dialog shows no entry for the Command
  Entry window although the text mentions it.
* **OQ-U3** Shader option defaults and ranges.
* **OQ-U4** The "Windows Settings" dialog contents are deferred to the plotting chapter (7.1.1).
