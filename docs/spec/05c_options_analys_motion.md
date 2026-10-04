# Spec 05c — Analysis Options: ANALYS and MOTION modules (§6.5.4)

**Source:** ACS SASSI Version 3 User Manual, extracted text lines 6333–6985
(`reference/acs-sassi.txt`). This is printed pages 147–161 (PDF pages 149–163; printed page = PDF page − 2).
Both dialog screenshots were viewed in the PDF and transcribed: the ANALYS tab (printed p. 147), the
coordinate-system figure (p. 151), the MOTION tab (p. 154), the "Output Motion Node List"
sub-dialog (p. 159) and the SRSSTF example (p. 156).

**Cross-references used** (outside the assigned range, consulted for consistency only):
- ANALYS / MOTION / NOUT / DAMP / THFILE / THTIT / ME / HOUSE / INCOH / WPASS batch commands
  (§9.2.5, 9.2.9, 9.2.16, 9.2.17, 9.2.20, 9.2.22, 9.2.23, 9.2.40, 9.2.41, 9.2.46).
- BINOUT / COMBACCDB (§9.18.3, 9.18.5) and CRITFREQ (§9.7.4).
- Module descriptions in Ch. 3 (lines 1975–2090), restart runs in §3.2 / §4.2 and
  Eq. (2.1) (printed p. 35).
- Check errors 44, 63–78 and 113–120 (Ch. 10).
- The STRESS-tab text that repeats the interpolation options (lines 7041–7080).

**Tag legend**

| Tag | Meaning |
|---|---|
| [CORE] | Needed for computational correctness |
| [IO] | File / input-output format or naming |
| [UI] | User interface, plotting, convenience |
| [ADV] | Advanced: incoherency, wave passage, multiple excitation, Option A/AA/PRO, batch-only restarts |
| *(impl. decision)* | The manual does **not** specify this. The rule given is the implementer's design choice and is listed again in §E (Open questions). |

The tabs of the **Analysis Options** dialog appear in this order (both screenshots):
`EQUAKE | SOIL | SITE | POINT | HOUSE | FORCE | ANALYS | MOTION | STRESS | RELDISP | NONLINEAR | AFWRITE`.
Each tab has **Ok / Cancel** buttons. The AFWRITE command writes the tab contents to the module input deck
`modelname.anl` (ANALYS) or `modelname.mot` (MOTION). It runs CHECK first, and if CHECK finds errors the
affected deck is not written.

---

## Part A — ANALYS module options

### A.1 Role of ANALYS [CORE]

- ANALYS is "the heart" of the code. It computes the soil-layer flexibility and impedance matrices,
  assembles every SSI system matrix, and solves the SSI problem for the complex response
  amplitudes (transfer functions) at each selected frequency.
- More than 90 % of total run time is spent in ANALYS. Total run time can be estimated beforehand as
  (time for one frequency) × (number of frequencies).
- ANALYS also controls the restart modes. There is **one initiation mode and three restart modes** in the
  GUI, plus one batch-only restart (Mode 6).
- Per-frequency computational steps (from Ch. 2/3; restated here because the options act on them):
  1. Form the soil flexibility (compliance) matrix at the interaction nodes (from POINT results, FILE3).
  2. Invert it to get the free-field impedance matrix: `X_ff(ω) = K + iD = (f + i g)^(-1)` (Eq. 4.3).
  3. Form the load vector. **Seismic:** `P = X_ff · U'_f` (free-field motion at the interaction nodes,
     including incoherency / wave passage / multiple-excitation modifications). **Vibration:**
     nodal forces read from FILE9.
  4. Assemble and solve the SSI system (Eq. 2.1, transcribed from printed p. 35):

     ```
     [ C^s_ii − C^e_ii + X_ii    −C^e_iw + X_iw    C^s_is ] {U_i}   { X_ii U'_i + X_iw U'_w }
     [ −C^e_wi + X_wi           −C^e_ww + X_ww     0      ] {U_w} = { X_wi U'_i + X_ww U'_w }
     [ C^s_si                    0                 C^s_ss ] {U_s}   { 0                     }
     ```
     Here `C(ω) = K − ω² M` (Eq. 2.2), with `K` complex (material damping enters through complex moduli),
     `s` = superstructure DOFs, `i` = basement (structure–soil interface) DOFs,
     `w` = excavated-soil DOFs (interaction nodes inside the excavation), `e` = excavated soil and
     `X` = impedance at the interaction DOFs. The system is solved by **LU decomposition plus
     back-substitution** (§6.4.9 text).
  5. Store the complex transfer functions for **all SSI DOFs** in **FILE8**. For seismic analysis these are
     the **ATF** (total acceleration TF relative to the control motion). For external-force analysis they
     are the **DTF** (total displacement per unit reference load).
- FILE8 is read by **MOTION, STRESS and RELDISP**. It is also merged by **COMBIN**, which reads two FILE8s
  renamed `FILE81` and `FILE82`.

### A.2 ANALYS file inventory [IO]

| File | Direction | Produced by | Content / notes |
|---|---|---|---|
| `modelname.anl` | in | AFWRITE | ANALYS input deck (this tab) |
| `FILE1` | in (seismic) | SITE (Mode 2) | Free-field wave modes / motions at the selected frequencies |
| `FILE1X`, `FILE1Y`, `FILE1Z` | in (seismic, simultaneous cases) | SITE run 3×, renamed by the user | One FILE1 per input direction |
| `FILE9` | in (vibration) | FORCE | Load vector |
| `FILE9001` … `FILE9500` | in (vibration, multiple load cases) | FORCE run once per case, renamed by the user | One per load case |
| `FILE3` | in | POINT | Point-load (compliance) solution for the soil layering |
| `FILE4` (= `modelname.N4`) | in | HOUSE | Structure + excavated soil FE data |
| `COOSK`, `COOSM` | in | HOUSE | Structure/basement stiffness and mass matrices (fast solver) |
| `FILE77` / `FILE77001` … `FILE77050` | in (incoherent, multiple excitation) | HOUSE | Incoherent motion field / multiple-excitation info for X, Y and Z, one file per stochastic simulation |
| `DOFSMAP`, `FILE90`, `FILE91` | in (restart) | HOUSE | Additional files that restart runs need |
| `FILE8` | out | ANALYS | SSI solution database: complex TF at all frequencies for all DOFs |
| `FILE8X`, `FILE8Y`, `FILE8Z` | out | ANALYS (coherent seismic, simultaneous) | One per input direction |
| `FILE8001` … `FILE8150` | out | ANALYS (incoherent, ≤ 50 sims × 3 directions) | See A.5.4 for the numbering |
| `FILE8001` … `FILE8500` | out | ANALYS (vibration, ≤ 500 load cases) | One per load case |
| `COOXxxx`, `COOTKxxx` | out/in (restart) | ANALYS | `xxx` = 3-digit **frequency order number** (e.g. `COOX001`). Inferred content: COOX = impedance matrix, COOTK = factorized total stiffness. |
| `COOXI`, `COOTKI` | out/in (restart) | ANALYS | Index files for the COOX/COOTK sets |
| `FOUNSTIF`, `FOUNDASH`, `FOUNDAMP`, `FOUNIMPD` | out (global impedance option) | ANALYS | Text files: dynamic stiffness, viscous damping coefficient, effective damping ratio, absolute value of impedance |
| `FILE11` | out (global impedance option) | ANALYS | Large file written when the global "unconstrained" impedance option is on (Ch. 3) |
| `modelname_ANALYS.out` (listing) | out | ANALYS | Printed ATF (amplitude only, or Re/Im; see A.5.6) |

Restart files carry the frequency order number in their names. Because of this, only the **full sequential set**
of SSI frequencies can be reused in a restart unless the user renames the files. *(impl. decision: store
the frequency value inside each restart record and look records up by value, but keep the
`COOXnnn`/`COOTKnnn` names.)*

### A.3 ANALYS dialog transcription [UI]

Values shown are those in the manual screenshot (printed p. 147). Use them as the UI defaults unless noted.

| Group box | Control | Type | Screenshot value | Notes |
|---|---|---|---|---|
| Operation Mode | Solution / Data Check | radio | **Solution** | |
| Type of Analysis | Seismic / Foundation Vibration | radio | **Seismic** | Shared with MOTION/STRESS "Type of Analysis" |
| Mode Of Analysis | Initiation / New Structure / New Seismic Environment / New Dynamic Loading | radio | **Initiation** | |
| (bottom left) | Simultaneous Cases | integer edit | `1` | Value shown in light grey in the screenshot. Command default is 0 (see A.5.4). |
| (bottom left) | Save Restart Files | checkbox | ☐ | |
| Frequency Numbers | Take Frequency Numbers from File1 / File9 | checkbox | ☐ | |
| Frequency Numbers | Frequency Set Number | integer edit | `1` | Same value as the SITE tab |
| Control Motion Foundation Reference Point | X-Coordinate of Control Point | real | `0` | |
| 〃 | Y-Coordinate of Control Point | real | `0` | |
| 〃 | Z-Coordinate of Control Point | real | `0` | |
| 〃 | Coordinate Transformation Angle | real (deg) | `0` | |
| (incoherency box) | Coherent / Incoherent | radio | **Incoherent** (in the screenshot) | Mirrors HOUSE `<coh>`. *(impl. decision: default Coherent for a new model.)* |
| 〃 | Wave Passage Effects Included | checkbox | ☐ | Mirrors HOUSE `<wpass>` |
| 〃 | Free-Field Load / Free-Field Motion | radio | **Free-Field Load** | Active only when Incoherent |
| (bottom middle) | Delete Restart Files | checkbox | ☐ (disabled) | Not described in the text (see §E) |
| (bottom middle) | Print Amplitude Only | checkbox | ☐ | |
| Multiple Excitation | Use Multiple Excitation | checkbox | ☐ | Mirrors HOUSE `<me>` |
| 〃 | Input Motion Number | integer spinner | `1` | 1–10 (ME command) |
| 〃 | First Foundation Node | integer | `0` | Per active motion |
| 〃 | Last Foundation Node | integer | `0` | Per active motion |
| 〃 | X Coord. of Control Point | real | `0` | Per active motion |
| 〃 | Y Coord. of Control Point | real | `0` | Per active motion |
| 〃 | Z Coord. of Control Point | real | `0` | Per active motion |
| Global Impedance Calculations | No Impedance Calculations / Only Decoupled (Diagonal) Impedances / Full Rigid Body Impedance Matrix 6X6 | radio | **No Impedance Calculations** | |

UI binding rule: the Coherent/Incoherent, Wave Passage and Multiple Excitation controls (including the
motion table) show **the same model data** as the HOUSE tab. They must be one source of truth, not
copies.

### A.4 Batch command `ANALYS` (§9.2.5) ↔ dialog mapping [IO]

`ANALYS,<opmode>,<type>,<mode>,<save>,<prnt>,<fopt>,<ang>,<xc>,<yc>,<zc>,<impe>,[simul]`

| Arg | Values | Dialog field |
|---|---|---|
| `<opmode>` | 0 = complete solution, 1 = data check only | Operation Mode |
| `<type>` | 0 = seismic, 1 = foundation vibration | Type of Analysis |
| `<mode>` | 0 = initiation, 1 = new structure, 2 = new seismic environment, 3 = new dynamic loading | Mode of Analysis. Internal solver mode numbers are 1, 2, 3, 3 (see A.5.3). |
| `<save>` | 0 = do not save restart files, 1 = save | Save Restart Files |
| `<prnt>` | 0 = print complex TF, 1 = print amplitude only | Print Amplitude Only |
| `<fopt>` | 0 = frequencies from the frequency set, 1 = from File1 (or File9) | Take Frequency Numbers from File1/File9 |
| `<ang>` | degrees | Coordinate Transformation Angle |
| `<xc>,<yc>,<zc>` | length | X/Y/Z-Coordinate of Control Point (foundation reference / control point) |
| `<impe>` | 0 = none, 1 = six diagonal terms, 2 = full 6×6 rigid-body matrix | Global Impedance Calculations |
| `[simul]` | integer, default **0** | Simultaneous Cases |

Not in the ANALYS command:
- **Frequency Set Number** comes from the SITE command `<freq>`.
- **Coherent/Incoherent, Wave Passage, Multiple Excitation** come from `HOUSE,…,<coh>,<wpass>,<me>,…`
  and `ME,<no>,<nfirst>,<nlast>,<xc>,<yc>,<zc>`.
- **FFL/FFM** has **no command argument anywhere**. *(impl. decision: add an extra optional trailing
  argument `[ffm]` to our ANALYS command: 0 = FFL (default), 1 = FFM.)*

### A.5 Field-by-field specification

#### A.5.1 Operation Mode [CORE]
- **Solution**: full run.
- **Data Check**: read and verify the inputs (file presence, frequency survey, DOF maps, option
  consistency) and write the listing, but do not factorize or solve.

#### A.5.2 Type of Analysis [CORE]
- **Seismic**: the load vector comes from FILE1 (free-field) and the impedance. FILE8 holds the **ATF**
  (total acceleration per unit control-point acceleration).
- **Foundation Vibration**: the load vector comes from FILE9 (FORCE module). FILE8 holds the **DTF**
  (total displacement per unit reference force amplitude). FORCE builds the load phasing from the load
  factor and arrival time t0, i.e. complex phasing `factor · e^{−iω t0}` relative to the reference history.
- Seismic and vibration cannot be mixed in one run. Combine them afterwards by superposition (Ch. 1.5.2).

#### A.5.3 Mode of Analysis (initiation and restart) [CORE]

| GUI choice | Cmd `<mode>` | Solver mode | Typical trigger | Inputs needed | What is recomputed | What is reused |
|---|---|---|---|---|---|---|
| **Initiation** | 0 | Mode 1 | New problem | FILE1 (or FILE9), FILE3, FILE4/`modelname.N4`, COOSK, COOSM (+FILE77 if incoherent / ME) | Everything: impedance, assembly, factorization, load, solve | — |
| **New Structure** | 1 | Mode 2 | Structural properties/geometry changed, **same interaction nodes and same soil layering**. Also used for every nonlinear (equivalent-linear) SSI iteration. | New FILE4 (re-run HOUSE), old FILE1/FILE9, `COOXxxx`, COOXI, DOFSMAP, FILE90, FILE91 | Structure matrices, assembly, factorization, load, solve | Impedance matrices `X_ff` (COOX) |
| **New Seismic Environment** | 2 | Mode 3 | Seismic environment changed (wave type, incidence, control location; or a new incoherent simulation from HOUSE). FILE1 changes, restart files unchanged. | New FILE1 (and/or FILE77), `COOXxxx`, `COOTKxxx`, COOXI, COOTKI, DOFSMAP, FILE90, FILE91 | Load vector and back-substitution only | `X_ff` (to form `X_ff·U'_f`) and the factorized system (COOTK) |
| **New Dynamic Loading** | 3 | Mode 3 (same as above) | New set of external forces | New FILE9 (or FILE9xxx), restart files | Load vector and back-substitution | Factorized system |
| *New Load Vector* | — (batch only) | **Mode 6** | Special change of the seismic load vector, e.g. per-level incoherent analysis of embedded models | Requires vendor tech support | — | — [ADV] |

Notes:
- Mode 2 needs only the `COOXxxx` files. Modes 3/3 need both `COOXxxx` and `COOTKxxx` (see A.5.5).
- A time-history change alone (same wave field) does **not** need ANALYS. Re-run only MOTION (and STRESS)
  (Ch. 3.2a). A time-history change of external loads with an unchanged load pattern also needs only MOTION.
- Mode 6 is not exposed in the GUI. *(impl. decision: implement it as "user-supplied load vector per
  frequency" for advanced use. The manual does not specify it.)*

#### A.5.4 Simultaneous Cases [CORE]/[IO]
Number of dynamic load cases solved in one ANALYS run, reusing one factorization per frequency for
many right-hand sides.

| Situation | Value to enter | Inputs read | FILE8 outputs |
|---|---|---|---|
| Classic single case | 0 (command default) | `FILE1` or `FILE9` | `FILE8` |
| Coherent seismic, X+Y+Z in one run | **1** | `FILE1X`, `FILE1Y`, `FILE1Z` | `FILE8X`, `FILE8Y`, `FILE8Z`. The listing contains the ATF for all three directions. |
| Incoherent seismic (stochastic) | `Ns` = number of simulations, **≤ 50**, **equal to the HOUSE value** | `FILE77001…FILE77Ns` (+ FILE1X/Y/Z) | `FILE8001 … FILE8(3·Ns)`, up to `FILE8150` |
| External force / vibration | `Nl` = number of load cases, **≤ 500** | `FILE9001 … FILE9Nl` | `FILE8001 … FILE8Nl` |

- Incoherent numbering (Ch. 4.2.5): simulation 1 → `FILE8001` (X), `FILE8002` (Y), `FILE8003` (Z),
  simulation 2 → `FILE8004…FILE8006`, and so on. In general, simulation `s` (1-based), direction `d`
  (1 = X, 2 = Y, 3 = Z) → `FILE8{3(s−1)+d:03d}`.
- **Requirements for seismic simultaneous cases** (WARNING): run SITE three times beforehand and copy each
  FILE1 to `FILE1X`, `FILE1Y`, `FILE1Z`.
  - X: **SV** waves, control direction **x'**, incidence angle 0.
  - Y: **SH** waves, control direction **y'**, angle 0.
  - Z: **P** waves, control direction **z**, angle 0.
  - The ANALYS **Coordinate Transformation Angle must be 0**.
- **Vibration**: run FORCE once per case and rename each output `FILE9001…FILE9500`.
- Capacity limit WARNING: the feasible number of cases depends on model size and RAM. The original code
  fails with "access violation". *(impl. decision: estimate memory per right-hand side up front and raise
  a clear error before solving.)*

#### A.5.5 Save Restart Files / Delete Restart Files [IO]
- **Save Restart Files** (checked, `<save>=1`): keep `COOXxxx` and `COOTKxxx` (+ `COOXI`, `COOTKI`)
  for later "New Seismic Environment" / "New Dynamic Loading" restarts, or only `COOXxxx` when
  only "New Structure" restarts are needed. Restart runs also need `DOFSMAP`, `FILE90`, `FILE91` from
  HOUSE.
- Unchecked: the COOX/COOTK files are not kept.
- **Delete Restart Files** appears (disabled) in the dialog but is not documented. *(impl. decision:
  enabled only for restart modes; deletes COOX*/COOTK* after a successful run.)*

#### A.5.6 Print Amplitude Only [IO]
- Checked (`<prnt>=1`): the listing prints the ATF **amplitude** for all nodes at every solved frequency.
- Unchecked (`<prnt>=0`): it prints the **real and imaginary parts** separately for **all 6 DOFs per node**.
- Complex amplitude/phase output for chosen nodes is better obtained from MOTION (Save Complex TF).

#### A.5.7 Frequency Numbers [CORE]
- SSI frequencies are integer multiples of the frequency step: `f_i = NFREQ_i · DF`, with
  `DF = 1/(DT · NFFT)` (Eq. 4.1). For harmonic-only vibration analysis DF can be given directly. **At most
  500 SSI frequencies.**
- **Take Frequency Numbers from File1/File9** checked (`<fopt>=1`): solve every frequency stored in
  FILE1 (seismic) or FILE9 (vibration).
- Unchecked (`<fopt>=0`): solve the frequencies of the **Frequency Set Number** (same set as the SITE tab;
  sets are defined with `FREQ,<ndx>,…` and listed with `LFREQ`).
- **Frequency survey** (required behaviour): before solving, check that every requested frequency exists in
  the input files (FILE1/FILE9 and FILE3). If any is missing, **ANALYS stops** with an error.
  *(impl. decision: match on the integer frequency number `NFREQ = round(f/DF)`, i.e. |f − NFREQ·DF| <
  1e-6·DF.)*
- **Subsets**: the full set can be split into subsets run separately (same or different computers) and
  merged afterwards with COMBIN (FILE81 + FILE82 → FILE8). A new frequency set can be solved later and
  merged with old results. All subsets must use frequencies that exist in the same FILE1/FILE3 (or FILE9).
- Recommended numbers (cross-ref §1.5.4 / §4.1.2): 40–80 for simple stick models; 100–200 for complex FE
  models (coherent); **200–300 for incoherent**. Use a denser grid where peaks are expected in nonlinear
  analysis.

#### A.5.8 Control Motion Foundation Reference Point and Coordinate Transformation Angle [CORE]
- **X-, Y-Coordinate of Control Point**: horizontal location of the control point for **oblique
  (inclined) input motion**, i.e. the origin of the horizontal phase term. It is also the foundation
  reference point for global rigid-body impedance (A.5.10).
- **Z-Coordinate of Control Point**: used **only** as the z of the foundation reference point for the
  global impedance calculation.
- **Coordinate Transformation Angle** `a` (degrees): angle between the SITE local axis x' and the global
  (HOUSE) axis x. Figure on printed p. 151: x' is rotated by `a` **counter-clockwise** from x about the
  vertical axis, y' is perpendicular to x', and z' coincides with z. SITE defines the control motion in
  x'y'z'. ANALYS transforms it to xyz:

  ```
  [u_x]   [cos a  −sin a  0] [u_x']
  [u_y] = [sin a   cos a  0] [u_y']
  [u_z]   [  0       0    1] [u_z']
  ```
  The horizontal propagation coordinate of node j for inclined waves is
  `x'_j = (x_j − x_c) cos a + (y_j − y_c) sin a`. The free-field motion has the form
  `u'_f(x') = U'_f exp[i(ωt − k x')]` (Eq. 2.3), so the node factor is `exp(−i k x'_j)`, with k the complex
  horizontal wave number from SITE. For vertically propagating waves k = 0 and the control point x/y do not
  matter.
- Valid range: Check Error 63 requires the angle to be in (0, 360) degrees. *(impl. decision: accept
  [0, 360).)*
- Must be **0** when Simultaneous Cases is used for seismic analysis (A.5.4).

#### A.5.9 Coherent/Incoherent, Wave Passage, Free-Field Load (FFL) vs Free-Field Motion (FFM) [ADV]
- **Coherent/Incoherent** and **Wave Passage Effects Included**: the same settings as on the HOUSE tab.
  HOUSE produces the incoherent-field information (FILE77). Wave passage uses the WPASS data:
  `<appv>` = apparent velocity of Line D, `<ang>` = angle of Line D with x, `<cohf>` = directional
  coherence factor. *(impl. note from the HOUSE tab: delay factor `exp(−iω s_j / V_app)`, with
  `s_j = (x_j − x_c) cos α_D + (y_j − y_c) sin α_D`. V_app = 1e9 effectively disables it.)*
- **Free-Field Load (FFL)** — the **recommended default**, the "official" approach accepted by EPRI and
  USNRC and validated in the 2007 EPRI studies (Short et al., 2007). The incoherency randomization is
  applied to the **seismic load vector** at the interaction nodes:
  `P_coh = X_ff · U'_f` (coherent), then `P_incoh,j = S_j(ω) · P_coh,j` per interaction DOF, where
  `S_j` is the incoherency factor/realization for node j from FILE77 *(impl. decision on the exact form of
  S_j; it belongs to the HOUSE spec)*.
  - Slightly conservative for elastic foundations. It reduces the coupling between forces at different
    plan locations, which makes the incoherency effect look slightly smaller.
  - The conservatism also covers uncertainty in using the isotropic Abrahamson coherency model, which is
    the same in all directions.
  - Valid from surface to deeply embedded structures.
- **Free-Field Motion (FFM)**: randomize the **free-field motion** at the interaction nodes, then form
  `P = X_ff · (S ⊙ U'_f)`, i.e. multiply the coherent impedance by the incoherent motion vector.
  - Gives results identical to FFL for **surface rigid foundations on rock**, and slightly lower for
    elastic foundations.
  - May be slightly more accurate than FFL at low frequency for surface structures. In some cases FFL can
    bias a few low-frequency incoherent ATF amplitudes slightly above the coherent ones.
  - WARNING: **not recommended for embedded structures**. Mixing incoherent motion with the coherent soil
    impedance creates artificial wave scattering in the excavation cavity. FFM is accurate only for surface
    structures.
- WARNING (phase adjustment, also relevant to MOTION B.5.6): adjusting the complex response phase is a
  departure from physics-based SSI. In the EPRI studies, stochastic simulation **with** phase adjustment
  matched the deterministic SRSS (zero-phase) results and was conservative for rigid basemats. It is
  **not guaranteed** to be conservative everywhere for elastic foundations (Ghiocel 2016a).
- *(impl. decision: FFL/FFM is greyed out unless Incoherent is selected.)*

#### A.5.10 Global Impedance Calculations [CORE] (diagnostic)
Options (`<impe>`):

| Value | Dialog label | Output |
|---|---|---|
| 0 | No Impedance Calculations | none |
| 1 | Only Decoupled (Diagonal) Impedances | Six diagonal terms per frequency → `FOUNSTIF`, `FOUNDASH`, `FOUNDAMP`, `FOUNIMPD` (+ `FILE11`) |
| 2 | Full Rigid Body Impedance Matrix 6X6 | Full 6×6 complex matrix per frequency (+ the same four files, *impl. decision*) |

Algorithm (rigid-body integration of nodal impedances; reference point = (X, Y, Z-Coordinate of Control
Point)):
1. For each interaction node j with offset `d_j = r_j − r_0 = (dx, dy, dz)`, form the rigid-body
   kinematic matrix (3×6, for translations u0 and small rotations θ about r_0) from
   `u_j = u_0 + θ × d_j`:
   ```
   T_j = [ 1 0 0   0    dz  −dy ]
         [ 0 1 0  −dz   0    dx ]
         [ 0 0 1   dy  −dx   0  ]
   ```
2. Stack the matrices: `T = [T_1; T_2; …; T_n]` (3n × 6), ordered like the X_ff interaction DOFs.
3. `K_G(ω) = Tᵀ · X_ff(ω) · T` (6×6 complex). DOF order is (X, Y, Z, XX, YY, ZZ).
4. Report per frequency `f` (ω = 2πf) for each term `K_lm`:
   - dynamic stiffness `k = Re K` → **FOUNSTIF**
   - viscous damping coefficient `c = Im K / ω` → **FOUNDASH**
   - effective damping ratio `ξ = Im K / (2 Re K)` → **FOUNDAMP** *(impl. decision: this definition; take
     the absolute value in the denominator if Re K < 0)*
   - absolute value `|K|` → **FOUNIMPD**

Limitations and warnings:
- Applies only to **3D** SSI models, with or without embedment.
- These are **"unconstrained"** soil impedances: only the distributed soil stiffness (X_ff) is included, as
  if the foundation were infinitely flexible. **For surface foundations they equal the rigid-foundation
  impedances.** For embedded models they differ from the rigid-foundation impedances and are computed
  correctly only with the **FI-FSIN** method.
- Ch. 3 WARNING: apply global impedances only to **surface stick models with rigid basemats**. Accuracy is
  crude for elastic foundations. For a two-step analysis use **Option A** instead [ADV, separate manual].

#### A.5.11 Multiple Excitation [ADV]
- **Use Multiple Excitation**: same as on the HOUSE tab. Activates non-uniform (variable-amplitude) input
  at different parts of the foundation. HOUSE requires Wave Passage to be selected as well.
- **Input Motion Number**: active motion, 1–10 (ME command). All fields below apply to that motion.
- **First Foundation Node / Last Foundation Node**: node range of the interaction nodes driven by this
  motion.
- **X/Y/Z Coord. of Control Point**: location where this input motion applies. The ME command says these are
  "not used in this version" (see §E).
- Spectral amplification ratios (`AMP,<no>,<a1>…<a100>`) scale the reference motion for motion `<no>`
  (complex if HOUSE `<cmplxspec>=1`). There is one ratio per frequency of the selected set (Error 119), each
  in [0, 10] (Error 118). *(impl. note: the free-field motion at the interaction nodes in
  [first, last] of motion n is multiplied by a_n(f).)*

### A.6 ANALYS reference algorithm (pseudo-code) [CORE]

```
read .anl; resolve frequency list F (set or all-in-file); survey F against FILE1/FILE9 and FILE3 -> stop if missing
determine load cases L (A.5.4)
for each f in F (order number q = 1..NF):
    ω = 2π f
    if mode == 1 (initiation) or not exists COOX[q]:
        F_soil = compliance at interaction DOFs from FILE3        # POINT results
        Xff = inv(F_soil)                                          # Eq. 4.3
        if save_restart: write COOX[q]
    else: Xff = read COOX[q]
    if mode in (1, 2):
        Cs = Ks − ω² Ms  (structure, from FILE4/COOSK/COOSM; complex Ks carries damping)
        Ce = Ke − ω² Me  (excavated soil)
        A  = assemble(Cs) − assemble_iw(Ce) + assemble_iw(Xff)    # Eq. 2.1
        LU = factorize(A); if save_restart: write COOTK[q]
    else:  LU = read COOTK[q]                                      # mode 3
    for each load case ℓ in L:
        if seismic: U'f = free-field motion at interaction DOFs from FILE1[ℓ]
                      -> rotate by angle a, oblique phase exp(-i k x'), wave passage, multiple excitation
                      -> P = Xff·U'f (FFL: then apply incoherency factors to P; FFM: apply to U'f first)
        else:       P = FILE9[ℓ] load vector at f
        U = LU.solve(P)                     # TF for all DOFs (ATF seismic / DTF vibration)
        append (f, U) to FILE8[ℓ]
    if impe > 0: K_G = Tᵀ Xff T -> FOUN* files
    print listing (amplitude only or Re/Im for 6 DOF/node)
```

### A.7 ANALYS validation checks (UI CHECK / module) [IO]
- Error 44: Frequency Set `<s>` Is Not Defined. Error 120: Frequency set `<i>` is empty.
- Error 63: Illegal Coordinate Transformation Angle (outside (0, 360) degrees).
- Errors 115–119 (multiple excitation): no ME data defined; illegal first/last node for motion `<i>`;
  amplification ratio outside 0–10; number of ratios ≠ number of frequencies in the set.
- Errors 113/114 (wave passage): apparent velocity ≤ 0; directional coherence factor < 0.
- Runtime: requested frequency missing from FILE1/FILE9/FILE3 → stop. Simultaneous cases above the limits
  (1 coherent; ≤ 50 incoherent; ≤ 500 vibration) → error. Restart mode without the required restart files
  → error naming the missing file.

---

## Part B — MOTION module options

### B.1 Role and pipeline of MOTION [CORE]
MOTION computes SSI **acceleration time histories, transfer functions and in-structure response spectra
(ISRS)** at the selected nodes. Seismic pipeline:
1. Read the control-motion acceleration history and transform it to the frequency domain (FFT).
2. Read the computed complex TFs from **FILE8** for the selected output nodes/DOFs.
3. **Interpolate** them to every Fourier frequency, using the selected complex-frequency scheme.
4. Multiply ("convolve") the TF with the control-motion Fourier transform.
5. Inverse FFT back to the time domain.
6. Output the time histories directly and/or convert them to response spectra.

MOTION needs only FILE8 as SSI input, plus the time-history file. Its input deck is `modelname.mot`
(AFWRITE). Ch. 6.4.11 says the original interpolation uses "a two SDOF transfer function model with five
parameters". This is the basis of B.6.4.

### B.2 MOTION file inventory [IO]

| File | Dir. | Content |
|---|---|---|
| `modelname.mot` | in | MOTION input deck |
| `FILE8` (or `FILE8X/Y/Z`, `FILE8nnn`) | in | SSI TF database. *(impl. decision: let the .mot name the FILE8 to read; the original expects `FILE8`, so the user renames.)* |
| Control-motion file (e.g. `acc_X_8192.acc`) | in | Acceleration history (B.5.10) |
| `SRSSTF.txt` | in (optional) | SRSS incoherent combination control (B.5.4) |
| `CONTTRS.txt` | in (optional) | List of external `.ACC` files for RS-only conversion (B.5.11) |
| `xxxxxTR_y.TFU` | out | Computed (uninterpolated) TF at the SSI frequencies; y ∈ {X, Y, Z} |
| `xxxxxTR_y.TFI` | out | Interpolated TF at the Fourier frequencies (RELDISP input, which needs complex values) |
| `xxxxxTR_y.ACC` | out | Acceleration time history |
| `xxxxxTR_yzz.RS` | out | Response spectrum for the damping value with order number zz (01, 02, …) |
| `xxxxxR_yy.ACC` etc. | out | Rotational DOFs: `R_XX`, `R_YY`, `R_ZZ` (naming confirmed by the RELDISP tab) |
| `name.RSO` | out | RS of an external `.ACC` file listed in CONTTRS |
| `FILE13` | out (optional) | Baseline-corrected absolute acc/vel/disp, 4 columns: time, acc, vel, disp |
| `FILE12` | out (optional) | Acceleration histories and response spectra |
| `\TFU\TFU_<freq>_<fnum>` | out (restart frames) | TF frame at one frequency, e.g. `\TFU\TFU_000.02_00001` |
| `\ACC\ACC_<time>_<fnum>`, `\ACCR\…` | out (restart frames) | Acceleration frame at one time step, e.g. `\ACC\ACC_00.000_00001` (ACCR = rotations) |
| `\RS\RS##_<freq>_<fnum>` | out (restart frames) | RS frame for damping ## at one frequency, e.g. `\RS\RS01_000.10_00001` |
| `ACC_max.txt` | out | Maximum-acceleration (ZPA) frame |
| `Frames.txt` | out | Post-processing frame list (Table 3.1) |
| `Modelname_ACC.bin` | out (optional) | Binary acceleration database (BINOUT `[mot]=1`). Combine X/Y/Z with `COMBACCDB,<Xfile>,<Yfile>,<Zfile>,<Comb>`. |
| `modelname_MOTION.out` | out | Listing (printed TF plots, maxima, optional tables) |

Node number field: 5 digits zero-padded for models with fewer than 100,000 nodes, 6 digits otherwise.
Name layout (Table 3.1):
- characters 1–5: node number;
- characters 6–9: DOF tag `TR_X`, `TR_Y`, `TR_Z`, `R_XX`, `R_YY`, `R_ZZ`;
- characters 10–11 (RS only): damping order number.

Example: `00012TR_X01.RS`.

### B.3 MOTION dialog transcription [UI]

Values are those in the manual screenshot (printed p. 154). They are example values. Proposed defaults
are in the last column.

| Group box | Control | Type | Screenshot | Proposed default / notes |
|---|---|---|---|---|
| Operation Mode | Solution / Data Check | radio | Solution | Solution |
| Type of Analysis | Seismic / Foundation Vibration | radio | Seismic | Mirrors ANALYS |
| Baseline Correction | No Correction / With Correction | radio | No Correction | No Correction |
| Response Spectrum Data | First Frequency | real (Hz) | `1` | **0.1** (SRP 3.7.1 practice) |
| 〃 | Last Frequency | real (Hz) | `8192` | **100** |
| 〃 | Total Number of Freq. Steps | integer | `32` | **301** |
| 〃 | Damping Ratios | text list | (empty) | e.g. `0.02 0.05`. Fractions; separators blank / tab / `,` / `;` |
| Output Control | Output Only Transfer Functions | checkbox | ☐ | ☐ |
| 〃 | Save Complex Transfer Functions | checkbox | ☐ | ☐ (☑ if RELDISP will be run) |
| 〃 | Save FILE 12 or FILE 13 | integer | `0` | 0 (none), 1 (FILE13), 2 (FILE12) |
| 〃 | Total Duration to be Plotted | real (s) | `0` | 0 = full record *(impl. decision)* |
| 〃 | Incoherent SSI — **Input** button | button | — | Opens the SRSSTF.txt editor (B.5.4) |
| 〃 | Interpolation Option | integer 0–6 | `1` | **1** coherent; **6** incoherent |
| 〃 | Phase Adjustment | integer 0–1 | `1` | **0** coherent |
| 〃 | Smoothing Parameter | real | `1` | **0** coherent |
| Nodal Output | Node List | list box (entries such as `1-30`) + **Add / Edit / Delete** | `1-30` | Each entry is a node list with its own flags |
| 〃 | Direction | radio X / Y / Z / XX / YY / ZZ | X | Flags below are stored per (list, direction) |
| 〃 | Printed Plot of Transfer Function | checkbox | ☑ | |
| 〃 | Save Time History of Requested Response | checkbox | ☐ | |
| 〃 | Plot Time History of Requested Response | checkbox | ☑ | |
| 〃 | Plot Acceleration and Velocity R.S. | checkbox | ☑ | Deprecated ("does not work in newer versions; will be taken out") |
| 〃 | Save Acceleration and Velocity R.S. | checkbox | ☑ | |
| 〃 | Print Maximum Requested Response | checkbox | ☑ | |
| Acceleration Time History Data | Nr. of Fourier Components | integer (power of 2) | `8192` | = SITE value |
| 〃 | Time Step of Control Motion | real (s) | `0.005` | = SITE value |
| 〃 | Multiplication Factor | real | `1` | 1 |
| 〃 | Max Value for Time History | real | `0` | 0 (exactly one of Factor/Max may be non-zero) |
| 〃 | First Record | integer | `1` | 1 |
| 〃 | Last Record | integer | `5000` | last record in the file |
| 〃 | Title | text | `acc_X_8192` | THTIT |
| 〃 | File | path | `C:/test/tshell/acc_X_8192.acc` | THFILE |
| 〃 | File Contains Pairs Time Step - Accel. | checkbox | ☐ | ☐ |
| Convert Time History to Response Spectrum | Select External Files | checkbox | ☐ | |
| 〃 | Input Time History Files | button | — | Browse; writes CONTTRS.txt |
| Post Processing Options | Save TF in All Points | checkbox | ☐ | |
| 〃 | Save ACC in All Points | checkbox | ☐ | |
| 〃 | Save RS in All Points | checkbox | ☐ | |
| 〃 | Save Rotation for Ansys 11.0 | checkbox | ☐ | The text calls it "Save Rotations for ANSYS" |
| 〃 | Restart for TF | checkbox | ☐ | |
| 〃 | Restart for ACC | checkbox | ☐ | |
| 〃 | Restart for RS | checkbox | ☐ | |
| Binary Output Option | Save Binary Database | checkbox | ☐ | = BINOUT `[mot]` |

**Output Motion Node List** sub-dialog (opened by Add/Edit): a single **List** text field (example
`1, 3-6 10`) with **Ok / Cancel**.

### B.4 Batch commands ↔ dialog mapping [IO]

`MOTION,<opmode>,<out>,<step>,<dur>,<res>,<freq1>,<freq2>,<fstep>,<mult>,<max>,<rec1>,<rec2>,<fopt>,<bl>,<smo>,<cplx>,<cnvrt>,<pzadj>,<interp>`

| Arg | Values / meaning | Dialog field |
|---|---|---|
| `<opmode>` | 0 = complete solution, 1 = data check | Operation Mode |
| `<out>` | 0 = full output, 1 = only transfer functions | Output Only Transfer Functions |
| `<step>` | output time-history print step: 0 = print only table; >1 = print every `<step>`-th point | *(not in dialog; default 0; Error 67 if < 0)* |
| `<dur>` | total duration of time histories to be plotted (s) | Total Duration to be Plotted |
| `<res>` | not used (Data Check) | — |
| `<freq1>`, `<freq2>` | first/last RS frequency (Hz) | First / Last Frequency |
| `<fstep>` | total number of RS frequency steps | Total Number of Freq. Steps |
| `<mult>` | multiplication factor | Multiplication Factor |
| `<max>` | target peak value | Max Value for Time History |
| `<rec1>`, `<rec2>` | first/last record (default last = last in file) | First / Last Record |
| `<fopt>` | 0 = time step on the first line, then one acceleration per line; 1 = pairs per line | File Contains Pairs… |
| `<bl>` | **"0 – time domain, 1 – frequency domain"** per §9.2.22 | Baseline Correction (dialog: No/With). **Conflict, see §E.** |
| `<smo>` | smoothing parameter ("between 1 and 1000") | Smoothing Parameter |
| `<cplx>` | 0 = TFU/TFI amplitude only; 1 = amplitude **and phase (radians)** | Save Complex Transfer Functions |
| `<cnvrt>` | 0 = no RS for external files; 1 = compute RS for the CONTTRS files | Select External Files |
| `<pzadj>` | phase adjustment 0/1 | Phase Adjustment |
| `<interp>` | interpolation option 0–6 | Interpolation Option |

Related commands:
- `NOUT,<dir>,<code1>,…,<code6>,<node list>`: adds a nodal output request. `<dir>`: 1 = x, 2 = y, 3 = z,
  4 = xx, 5 = yy, 6 = zz. *(impl. decision: `<code1>…<code6>` = the six checkboxes in dialog order —
  Printed Plot of TF, Save TH, Plot TH, Plot RS, Save RS, Print Max — each 0/1.)*
- `DAMP,<d1>,…,<d10>`: adds non-zero damping ratios to the RS list. `<d1> = 0` deletes the list. Up to 10
  per command.
- `THFILE,<file>`: time-history file. `THTIT,<title>`: time-history title.
- `BINOUT,[mot],[str],[reldisp]`: `mot` 0/1 writes `Modelname_ACC.bin`. A blank argument leaves the flag
  unchanged.
- The Save FILE12/13 flag, Baseline (with/without), Phase adjustment, the post-processing flags and the
  SRSS input have no separate documented command beyond the above. *(impl. decision: add `MOTPP,…`
  for the post-processing flags in our command set. Name it so it cannot be confused with original
  commands.)*

### B.5 Field-by-field specification

#### B.5.1 Operation Mode / Type of Analysis [CORE]
- Data Check: validate the inputs (B.8) and echo them; no FFT and no output files.
- **Type of Analysis** must match the FILE8 being read (seismic ATF vs vibration DTF). The value is the same
  as on the ANALYS tab.
  - **Seismic**: the requested response is **acceleration**, and TFs are total-acceleration TFs.
  - **Foundation Vibration**: TFs are **total-displacement** TFs, and the requested response type is set by
    the "Output Control" selection (displacement / velocity / acceleration; see §E — that control is not
    visible in the seismic-mode screenshot).

#### B.5.2 Baseline Correction [CORE]
- **No Correction / With Correction**. "With Correction" applies the **classical Hudson–Housner
  time-domain** baseline correction to the computed response acceleration histories before integrating to
  velocity and displacement (algorithm in B.6.13).
- It gives only **approximate absolute displacements**. The user should check that the ground displacement
  is about zero at the end of the motion. Relative displacements between locations are better obtained as
  differences of absolute displacements, and much better from RELDISP.
- WARNING: use the **RELDISP** module, not baseline correction, for relative displacements. RELDISP works
  analytically in complex frequency from the ATFs.
- FILE13 is produced only when baseline correction is on (B.5.3).

#### B.5.3 Output Control fields [IO]/[CORE]
- **Output Only Transfer Functions** (`<out>=1`): extract only TFs at the selected nodes, i.e. write
  TFU/TFI. Skip convolution, time histories and RS.
- **Save Complex Transfer Function** (`<cplx>=1`): TFU and TFI files store **amplitude and phase
  (radians)**. Otherwise they store amplitude only. Not recommended for routine use because the Fourier
  phase is hard to interpret. **However, RELDISP requires complex TFI files**, so turn it on when RELDISP
  follows.
- **Save FILE13 or FILE12** (integer):
  - `0` = save neither.
  - `1` = save the **baseline-corrected absolute acceleration, velocity and displacement** histories in
    **FILE13**. No FILE13 is written if baseline correction is off. FILE13 can become extremely large for
    full post-processing of big models, so it is not recommended routinely.
  - `2` = save the acceleration histories and response spectra in **FILE12**.
- **Total Duration to be Plotted / for Extracted Time Histories** (`<dur>`, s): duration of the output
  histories. The extracted length is **20 % longer** than the given duration, to include part of the free
  vibration: `T_out = 1.2·dur`, clipped to `NFFT·DT`. *(impl. decision: dur = 0 → T_out = NFFT·DT.)*
  Error 68 if negative.

#### B.5.4 Incoherent SRSS (SRSSTF.txt) [ADV]
- The **Input** button creates/edits **`SRSSTF.txt`**. MOTION reads it to apply the **SRSS TF**
  deterministic incoherent approach. The modal FILE8s come from separate SSI analyses, one per incoherency
  mode. These are usually New Seismic Environment restarts after HOUSE is re-run for each mode (HOUSE
  `<nmodes> = −k` → only mode k).
- File format (free-format text):
  ```
  [# of modes n] [phase option p]
  [FILE8 for coherent analysis]     <- present ONLY if p = 1
  [FILE8 for mode 1]
  [FILE8 for mode 2]
  ...
  [FILE8 for mode n]
  ```
  When p = 1 the first file line is the coherent FILE8. For p = 0 that line is **omitted** and the list
  starts with the mode-1 FILE8.
- Phase options:
  - `0` = SRSS TF with **zero phase** (validated in the 2007 EPRI studies);
  - `1` = SRSS TF with the **coherent phase** (not validated by EPRI).
- The user must rename FILE8 after each ANALYS run. Manual example: 10 modes, phase option 1:
  ```
  10 1
  FILE8_coh
  FILE8_01
  ...
  FILE8_10
  ```
- Combination (per node/DOF and frequency): `|H_SRSS(f)| = sqrt( Σ_{k=1..n} |H_k(f)|² )`.
  - p = 0: `H(f) = |H_SRSS(f)|` (real, non-negative).
  - p = 1: `H(f) = |H_SRSS(f)| · exp(i·arg H_coh(f))`.
  - Algorithm detail in B.6.9.
- WARNING: the ACS SASSI SRSS follows the basic SRSS theory but may not reproduce all the "artifacts" of
  SASSI2010 / EPRI INCOH, so results may differ from those codes.
- The manual spells this file `SRSSTF.txt`, `SRSSTF.TXT` and (typo) `STRSSTF.txt`. *(impl. decision:
  case-insensitive lookup of `SRSSTF.txt`.)*

#### B.5.5 Interpolation Option [CORE]
Seven complex-frequency interpolation schemes for the ATF (the same codes are used for the STF in STRESS):

| Code | MOTION name | STRESS-tab name | Averaging |
|---|---|---|---|
| 0 | SASSI2000, dense overlapping windows | same | **weighted** averaging |
| 1 | Original SASSI 1982, non-overlapping windows | same | none |
| 2 | Dense overlapping windows | same | averaging |
| 3 | Only three overlapping windows | same | averaging |
| 4 | Non-overlapping windows with **one** position shift | "reduced shift" | none |
| 5 | Non-overlapping windows with **two** position shift | "large shift" | none |
| 6 | Complex bicubic spline interpolation | same | no windowing; needs a **denser** frequency grid |

- The algorithms are in B.6.4–B.6.6.
- WARNING: **option 6 is recommended for incoherent SSI** because it avoids interpolation overshoot. First
  check that the frequency grid is dense enough (usually > 200 frequencies) so the spline does not clip or
  smooth genuine incoherent ATF peaks.
- Check interpolation quality with the UI **Spectrum TFU-TFI** plot (computed vs interpolated). The
  **CRITFREQ** command (`CRITFREQ,<tol>,<minfilter>,<TF>,<Var>`) automatically finds frequencies where an
  interpolated peak departs from the nearby computed values by more than `tol` %, ignoring levels more
  than `minfilter` % below the global maximum. Add those frequencies to the SSI run.

#### B.5.6 Phase Adjustment [ADV]
- Intended for **incoherent** SSI with the **stochastic simulation** approach ("Simulation Mean" in EPRI
  2007) or the **deterministic linear AS** approach. It gives an approximate **upper-bound** solution by
  zeroing the ATF phasing of the incoherent responses.
- `0` — no adjustment. Keeps the physics-based complex ATF phase. This is the "theoretically exact" Monte
  Carlo approach. It was not used in EPRI 2007 because it gave lower ISRS than the other industry SRSS
  approaches.
- `1` — adjustment. It **strongly reduces, to nearly zero, the differential phase between neighbouring
  Fourier components**, by an amount that **depends on the smoothing parameter**. This avoids cancellation
  between neighbouring frequencies and creates a critical "minimum delay" motion. It is generally
  conservative, close to SRSS TF for stick models with rigid basemats, and was the EPRI-2007 "consensus"
  approach. Algorithm in B.6.8.
- WARNING: phase adjustment (like SRSS) distorts the **cross-correlation** between responses at different
  locations. Do **not** use the resulting motions for **multiple time history analysis** of secondary
  systems.
- STRESS uses the same parameter definition.

#### B.5.7 Smoothing Parameter [CORE]
- Filters spurious sharp spectral peaks and valleys that the classical SASSI interpolation schemes
  (**options 0–5 only**) can introduce into the interpolated complex ATF. It has **no effect for option 6**.
- Typical range **10–1000**. The batch command text says 1–1000. **0 = no smoothing** (B.6.7).
- WARNING: use **0 for coherent analysis** and a **non-zero value for incoherent analysis** (except with
  option 6, where it must be 0). Run sensitivity studies before fixing a production value. EPRI 2007 found
  ISRS practically unchanged for values from 10 to 500.
- WARNING: do not smooth when the spline (option 6) is used. Set it to zero.
- The parameter also controls the strength of the phase adjustment (B.5.6, B.6.8).

#### B.5.8 Nodal Output Data [UI]/[IO]
- The **Node List** holds lists of nodes that share the same output request. Use **Add** (new list),
  **Edit** (selected list) or **Delete** (selected list).
- List syntax: separators blank, tab, `,` or `;`. `a-b` adds the inclusive range a…b.
  Example `1, 3-6 10` → {1, 3, 4, 5, 6, 10}. *(impl. decision: allow spaces around `-`; reject
  descending ranges with an error.)*
- For each list **and each direction** (X, Y, Z, XX, YY, ZZ), six flags are stored:

  | # | Flag | Meaning |
  |---|---|---|
  | 1 | Printed Plot of Transfer Functions (*) | TF output (listing printer-plot; TFU/TFI files) |
  | 2 | Save Time History of Requested Response (**) | Write the history (`.ACC` for seismic) |
  | 3 | Plot Time History of Requested Response (**) | UI plot request |
  | 4 | Plot Acceleration and Velocity R.S. (***) | Deprecated; accept and ignore with a notice |
  | 5 | Save Acceleration and Velocity R.S. (***) | Write `.RS` files per damping |
  | 6 | Print Maximum Requested Response (**) | Peak value in the listing |

  (*) TF = total acceleration TF for seismic, total displacement TF for vibration.
  (**) Requested response = acceleration for seismic. For vibration it is the type selected in Output
  Control.
  (***) RS are computed independently of the Output Control selection, so **displacement response
  spectra cannot be requested**.
- Output requests on **constrained (fixed) DOFs are ignored**. *(impl. decision: emit a warning.)*
- WARNING: defining an output node **twice** gives incorrect results. The UI **removes duplicate nodes**
  when AFWRITE writes the MOTION deck. CHECK Error 66 reports duplicates.
  *(impl. decision: de-duplicate per (node, direction), merging flags by logical OR, and warn.)*
- If HOUSE node renumbering (optimizer) was used, the output nodes must use the **new** numbers
  (`.hounew` / `.map`), because FILE8 is in the renumbered system.

#### B.5.9 Response Spectrum Data [CORE]
- **First Frequency** and **Last Frequency** (Hz) of the RS.
- **Total Number of Freq. Steps**: number of RS frequency points. For nuclear safety work (US NRC SRP 3.7.1)
  use **0.1–100 Hz with at least 301 steps**.
  - WARNING: fewer steps degrade the UI **spectrum broadening** algorithm (BROADEN command). Use 301.
  - *(impl. decision: logarithmic spacing `f_i = f1·(f2/f1)^((i−1)/(n−1))`, i = 1…n. With 301 points over
    0.1–100 Hz this gives exactly 100 points per decade.)*
- **Damping Ratios**: list of critical-damping fractions, e.g. `0.02, 0.05`. Separators blank / tab / `,`
  / `;`. Each must be in (0, 1) (Error 72). The order in the list sets the `zz` file suffix
  (01, 02, …). DAMP allows up to 10 per command. *(impl. decision: maximum 10 damping values in total.)*

#### B.5.10 Acceleration Time History Data [CORE]/[IO]
- **Nr. of Fourier Components** (`NFFT`): FFT length. Must be a **power of 2**, ≤ **32,768**, equal to the
  SITE value. The UI rounds to the nearest power of 2 with Warning 9.
- **Time Step of Control Motion** (`DT`, s): equal to the SITE value. `DF = 1/(NFFT·DT)` must equal the DF of
  FILE8. *(impl. decision: error if the relative mismatch exceeds 1e-6.)*
- **Multiplication Factor** and **Max. Value for Time History**: scale the motion. Use **exactly one**.
  - Factor ≠ 0 and Max = 0: `a ← factor · a`.
  - Max ≠ 0 and Factor = 0: `a ← a · (Max / max|a|)`.
  - Both zero → Error 77. Both non-zero → Error 78. The manual's wording "use only if the other is blank"
    means blank = 0.
- **First Record / Last Record**: 1-based record range used from the file. Defaults are the first and last
  records. Errors 74–76.
- **Title**: label of the history (THTIT).
- **Acceleration History File**: full path (THFILE). It can be viewed with Plot / Time History. Error 73 if
  missing.
- **File Contains Pairs Time Step–Acceleration** (`<fopt>`):
  - checked: each line holds a (time, acceleration) pair *(impl. decision: DT is checked against the time
    column; non-uniform steps are an error)*;
  - unchecked: **line 1 = time step**, followed by one acceleration value per line.
- A "record" is one acceleration value (one line after the header) *(impl. decision)*.
- Requirement (Ch. 1.5.2): the Fourier period `NFFT·DT` must be much longer than the excitation. The
  trailing part ("quiet zone") is zero-padded and covers the free vibration. Its required length grows as
  the system's lowest damping decreases.

#### B.5.11 Convert Time History to Response Spectrum [UI]
- Convenience feature, independent of the SSI solution. It computes RS for **external** acceleration
  histories, e.g. histories obtained by algebraic summation outside the code.
- Check **Select External Files** (`<cnvrt>=1`) and browse with **Input Time History Files**. The chosen file
  names go into **`CONTTRS.txt`** (user-editable, one path per line *(impl. decision)*).
- The external files must have the **same one-column format, time step and duration** as the
  MOTION-computed accelerations (no shorter, no longer), and **the `.ACC` extension**. Output goes to the
  same base name with extension **`.RSO`**.
- The RS settings (frequencies, dampings, number of time steps) are the nodal-output ones (B.5.9).

#### B.5.12 Post Processing Options [UI]/[IO]
These generate data for **bubble, vector-TF and deformed-shape** plots at all nodes. They can be combined
with the Node Output requests.

| Option | Effect |
|---|---|
| Save TF in All Points | Save the computed **and** interpolated acceleration TF for all **translational** DOFs (`.TFU`, plus `.TFI` *impl. decision*) |
| Save ACC in All Points | Save acceleration histories for all translational DOFs (`.ACC`) |
| Save RS in All Points | Save acceleration RS for the selected dampings for all translational DOFs (`.RS`) |
| Save Rotation(s) for ANSYS (11.0) | Save acceleration histories for all **rotational** DOFs (`.ACC`, `R_XX/R_YY/R_ZZ`) |
| Restart for TF | Compute and save frames for animated complex-TF **vector** plots → `\TFU` |
| Restart for ACC | Frames for animated acceleration deformed shape or static bubble plots → `\ACC` (and `\ACCR` for rotations). Also `ACC_max.txt`. |
| Restart for RS | Frames for RS deformed-shape (animated) or bubble (static) plots → `\RS` |

- "Save" options write responses for all DOFs. "Restart" options write one **frame file per frequency step
  or time step** containing the values for **all active nodal DOFs**.
- Frame names (Table 3.2):
  - `RS##_freq_fnum`: `##` = damping number, `freq` formatted `000.10`, `fnum` 5-digit frame number;
  - `TFU_freq_fnum`, e.g. `TFU_000.02_00001`;
  - `ACC_time_fnum`, e.g. `ACC_00.000_00001`.
- The frame header gives the number of rows/columns (FRAMECOMBIN note). *(impl. decision: header line
  `nrows ncols`, then rows `node  vx vy vz`.)*

#### B.5.13 Save Binary Database [IO]
Stores all acceleration histories in one binary database **`Modelname_ACC.bin`** for fast
post-processing (BINOUT). Combine the X, Y and Z direction databases with **COMBACCDB**.
*(impl. decision: HDF5/NPZ container under the same name.)*

### B.6 Algorithms in implementable detail [CORE]

#### B.6.1 Conventions
- Time dependence is **e^{+iωt}** (Eq. 2.3, `exp[i(ωt − kx)]`). Forward transform
  `A_k = Σ_{n=0}^{NFFT−1} a_n e^{−2πi kn/NFFT}`; inverse `a_n = (1/NFFT) Σ_k A_k e^{+2πi kn/NFFT}`. This is
  the numpy `rfft`/`irfft` convention, and a TF computed by ANALYS under e^{+iωt} applies directly.
- Fourier frequencies: `f_k = k·DF`, k = 0…NFFT/2. Nyquist `f_Ny = 1/(2DT)`.
- Computed SSI frequencies from FILE8: `f_1 < … < f_N` (N ≤ 500), with complex values `H_j`. The cut-off
  frequency is `f_c = f_N`.
- Phase is reported in **radians** (cplx = 1). Damping is a fraction of critical. Frequency is in Hz,
  time in s, angles in degrees. Accelerations keep the input-file units; the ATF is dimensionless.

#### B.6.2 Control-motion preprocessing
1. Read the file (B.5.10), select records rec1…rec2 → `n` samples.
2. Scale (factor or max).
3. If `n > NFFT` → error ("time history longer than Fourier period"). Otherwise zero-pad to NFFT.
4. `A = rfft(a)` (length NFFT/2 + 1).

#### B.6.3 TF definition outside the computed range *(impl. decision)*
- `f_k > f_N`: `H(f_k) = 0`. Response content above the SSI cut-off is discarded, so the cut-off acts as a
  low-pass filter.
- `0 ≤ f_k < f_1`, **seismic**: linear complex interpolation between an anchor `H(0)` and `H_1`. `H(0)` is
  the rigid-body value: the projection of the control-motion direction onto the output DOF, i.e. 1 for the
  translational DOF parallel to the input after the coordinate transformation, `cos a`/`sin a` for oblique
  horizontal components, and 0 for the others and for rotations.
- `0 ≤ f_k < f_1`, **vibration**: `H(f_k) = H_1` (constant extrapolation; static compliance approximation).
- At the Nyquist bin, set `H` to its real part so `irfft` is consistent (usually `f_N < f_Ny`, so it is 0).
- The interior interpolation (B.6.4–B.6.6) uses only the computed points, plus the f = 0 anchor if
  `include_zero_anchor = true` (default true for option 6, false for the window schemes).

#### B.6.4 Local "two-SDOF, five-parameter" rational model (basis of options 0–5) *(impl. decision on exact form)*
For a window W of **p = 5 consecutive computed points** `{f_m, …, f_{m+4}}`, use local normalized frequency
`x = (f − c_W)/h_W`, with `c_W = (f_m + f_{m+4})/2` and `h_W = (f_{m+4} − f_m)/2`. Fit:

```
R_W(x) = (β1 + β2 x + β3 x²) / (1 + β4 x + β5 x²)        (5 complex parameters)
```
The fit is exact at the 5 window points. It comes from the linear system (p = 1…5):
```
β1 + β2 x_p + β3 x_p² − H_p β4 x_p − H_p β5 x_p² = H_p
```
Solve the 5×5 complex system (LU with partial pivoting). This quadratic-over-quadratic form can represent the
total-acceleration TF of a damped SDOF exactly (and approximately the near-resonance behaviour of two modes).
Safeguards:
- If the system is singular (cond > 1e12), or `min |1 + β4 x + β5 x²| < 1e-3` on `x ∈ [−1, 1]` (a near-real
  pole inside the window), replace R_W in that window by the complex cubic **Lagrange** polynomial through
  the 4 points nearest the target, and log a warning.
- If `N < 5`, all options 0–5 fall back to option 6, or to complex linear interpolation if N < 3.

#### B.6.5 Window schemes for options 0–5 *(impl. decision; derived from the option names)*
Let j be the interval index with `f_j ≤ f < f_{j+1}`. The candidate windows containing interval j have start
indices `M_j = {m : max(1, j−3) ≤ m ≤ min(j, N−4)}`.

| Opt | Window selection for interval j | Combination |
|---|---|---|
| **1** (SASSI 1982) | Fixed partition with stride 4: starts m ∈ {1, 5, 9, …}. Neighbouring windows share only end points. Interval j uses `m = 1 + 4⌊(j−1)/4⌋`, clipped to ≤ N−4 (the tail window is `N−4…N`). | single window |
| **4** (one-position shift) | Same as opt 1 but with starts {2, 6, 10, …}. Interval 1 uses m = 1. | single window |
| **5** (two-position shift) | Starts {3, 7, 11, …}. Intervals 1–2 use m = 1. | single window |
| **2** (dense, averaging) | All m ∈ M_j (up to 4 windows) | arithmetic mean of the complex values R_m(f) |
| **0** (SASSI2000, weighted) | All m ∈ M_j | weighted mean `Σ w_m R_m / Σ w_m`, `w_m = max(ε, 1 − |f − c_m|/h_m)`, ε = 1e-3 (a window counts most when f is near its centre) |
| **3** (three windows) | The 3 windows in M_j with smallest `|f − c_m|` (ties → smaller m). All of M_j if fewer than 3. | arithmetic mean |

At computed frequencies every scheme returns `H_j` exactly, because each R_m interpolates its nodes.

#### B.6.6 Option 6 — complex bicubic spline
- Fit independent **cubic splines** to `Re H` and `Im H` against f over the computed points (plus the f = 0
  anchor, B.6.3). *(impl. decision: not-a-knot end conditions, e.g. `scipy.interpolate.CubicSpline`
  default.)* Evaluate at every `f_k ≤ f_N`.
- No windowing and no smoothing. Ignore the smoothing parameter and warn if it is non-zero.
- Requires a dense grid (more than 200 frequencies for incoherent analysis).

#### B.6.7 Smoothing filter (options 0–5, S > 0) *(impl. decision)*
Goal: suppress interpolated excursions that are not supported by the computed neighbours, while leaving
values at the computed points unchanged. For f in interval j:
```
L(f)   = H_j + (H_{j+1} − H_j)(f − f_j)/(f_{j+1} − f_j)      # complex linear reference
A_max  = max(|H_j|, |H_{j+1}|),  A_min = min(|H_j|, |H_{j+1}|)
r(f)   = max(0, |H(f)|/A_max − 1) + max(0, 1 − |H(f)|/max(A_min, 1e-12·A_max))
g(f)   = 1 / (1 + S · r(f))
H_s(f) = L(f) + g(f) · (H(f) − L(f))
```
- S = 0 → identity (no smoothing).
- A large S pulls overshoots (peaks) and undershoots (valleys) toward the linear reference.
- Inside the band [A_min, A_max] nothing changes.
- This is why S must be 0 for coherent runs: genuine resonant peaks between computed points would be
  clipped.

#### B.6.8 Phase adjustment (option 1) *(impl. decision)*
Apply after interpolation (and after smoothing) on the Fourier grid k = 0…K, where `f_K ≤ f_N`:
1. `φ_k = unwrap(arg H_k)`, starting from `φ_0 = arg H(0)` (= 0 for a positive rigid-body anchor).
2. Residual phase factor ρ:
   - **Options 0–5:** `ρ = 1/(1 + S)`. S = 10 gives ρ ≈ 0.09; S = 1000 gives ρ ≈ 0.001 ("almost zero").
     S = 0 gives ρ = 1, i.e. no reduction. Warn that phase adjustment needs S > 0.
   - **Option 6:** `ρ = 0`, i.e. **complete zeroing**, regardless of S. The manual says S has no effect on
     the spline scheme and must be 0 there, yet it recommends both phase adjustment and the spline for
     incoherent analysis. The only consistent reading is that the spline path zeroes the phase completely.
3. `H'_k = |H_k| · exp(i · ρ · φ_k)`. Every differential phase `Δφ_k = φ_k − φ_{k−1}` is scaled by ρ, i.e.
   "reduced to almost zero depending on the smoothing parameter".
- A zero-phase TF gives an acausal (time-symmetric) response kernel. With the periodic FFT, part of the
  response wraps to the end of the Fourier period. This is expected. Keep the quiet zone long, and extract
  `T_out ≤ NFFT·DT`.
- If phase adjustment = 1 is used with options 0–5 and 0 < S < 10 (weak adjustment), warn.

#### B.6.9 SRSS TF combination (when SRSSTF.txt is active)
For each output node/DOF:
1. Read `H_k(f_j)` from each modal FILE8 (k = 1…n) and, if p = 1, `H_coh(f_j)`.
2. *(impl. decision, default "interpolate-then-combine")*: interpolate each modal complex TF to the
   Fourier grid with the selected option. Then `A_k = sqrt(Σ_k |H_k(f_k)|²)`. Then:
   - p = 0: `H = A`;
   - p = 1: `H = A·exp(i·arg H_coh,interp)`.
   - Alternative "combine-then-interpolate" (combine at the computed f_j, then interpolate the real
     amplitude) is available as a flag for benchmarking.
3. Convolve as in B.6.10. Phase adjustment should be 0 with SRSS (warn otherwise).

All modal FILE8s must share the same frequency list and DOF map. Otherwise → error.

#### B.6.10 Convolution and inverse transform
```
R_k = H_k · A_k   (k = 0…NFFT/2)
r   = irfft(R, NFFT)                 # response history, length NFFT, step DT
output r[0 : n_out], n_out = min(NFFT, round(1.2·dur/DT)) (or NFFT if dur = 0)
peak = max |r[0:n_out]|              # "Print Maximum Requested Response"
```
- **Identity check:** with H ≡ 1 for all k, r = a exactly. For the free-field DOF at the control point,
  r should reproduce the input up to f_c.

#### B.6.11 Foundation-vibration specifics
- `A` = FFT of the **reference load time history** (the history file is the force history), and `H` = DTF
  (displacement per unit load).
- Displacement `R = H·A`; velocity `R = iω_k·H·A`; acceleration `R = −ω_k²·H·A`, as selected in Output
  Control. RS (if requested) are computed from the acceleration response.

#### B.6.12 Response spectra (ISRS)
For each damping ζ in the list and each RS frequency f_i (ω = 2πf_i), integrate the base-excited SDOF
`ẍ + 2ζωẋ + ω²x = −r(t)` over the response history r (taken as piecewise-linear between samples).
*(impl. decision: use exact first-order-hold state transition, which is equivalent to Nigam–Jennings.)*
- `z = [x, ẋ]ᵀ`, `F = [[0, 1], [−ω², −2ζω]]`, `G = [0, −1]ᵀ`.
- Build `M = [[F, G, 0], [0, 0, 1], [0, 0, 0]]` (4×4) and `E = expm(M·Δt) = [[Φ, Γa, Γb], …]`.
- Step: `z_{n+1} = Φ z_n + Γa r_n + Γb (r_{n+1} − r_n)/Δt`, starting at rest.

Outputs:
- **SA** = max_n |ω² x_n + 2ζω ẋ_n| (absolute acceleration). This is written to `.RS` as two columns,
  frequency and SA.
- **SV** = max |ẋ_n| (relative velocity) for the "velocity R.S." listing *(impl. decision; pseudo-velocity
  ω·SD as an optional variant)*.
- At high frequency, SA → ZPA = max|r|.
- *(impl. decision: if f_i > 0.1/Δt, also evaluate on an FFT-upsampled (×4) history to reduce sampling
  error in the peak.)*

Run the integration over the full NFFT record (including the quiet zone), not only T_out, so the free
vibration is captured.

#### B.6.13 Baseline correction (Hudson–Housner, time domain) *(impl. decision on the exact variant)*
For each output acceleration history `r_n` (n = 0…M−1, step Δt):
1. Integrate (trapezoidal rule) to get `v` (v_0 = 0), then `d` (d_0 = 0).
2. Least-squares fit `v(t) ≈ c1·t + c2·t²/2`, i.e. the velocity effect of an acceleration baseline
   `c1 + c2·t`. This minimizes `Σ (v_n − c1 t_n − c2 t_n²/2)²`, the classical minimum mean-square-velocity
   criterion.
3. Corrected values: `r_c = r − (c1 + c2 t)`, `v_c = v − (c1 t + c2 t²/2)`, `d_c = ∫ v_c` (re-integrate).
4. Write to FILE13 (if Save = 1) per node/DOF: columns `time, acc, vel, disp`.
5. Report `d_c(T_end)` in the listing so the user can check that it is ≈ 0.

Units: if accelerations are in g, multiply by the gravity constant before integrating. *(impl. decision:
add a "units of acceleration" option, g or length/s², default g with the SITE/HOUSE gravity.)*

### B.7 Output file content *(impl. decision; the manual gives names and column semantics only)*
- `.TFU`: header line `# node dof analysis=seismic|vibration complex=0|1`, then rows
  `f  |H|` or `f  |H|  phase_rad` at the **computed** frequencies.
- `.TFI`: the same at **Fourier** frequencies `0 … f_N`.
- `.ACC`: one-column format, the same as the input non-pair format: line 1 = DT, then one value per line, so
  `.ACC` files can be read back as control motions and in CONTTRS.
- `.RS` / `.RSO`: two columns `frequency  SA` (the UI line-object reader counts one data column per RS file).
- `FILE13`: per node/DOF, a header line, then `t  acc  vel  disp`.

### B.8 MOTION validation checks [IO]
Each check's error number:
- **64**: no nodal output request. *(impl. decision: waived if any "…in All Points" option or
  `<cnvrt>=1` is set.)*
- **65**: illegal node in a list.
- **66**: node defined more than once.
- **67**: negative output time-history step.
- **68**: negative total duration.
- **69 / 70 / 71**: negative first RS frequency / last RS frequency / number of RS steps.
  *(impl. decision: also require 0 < f1 < f2 and n ≥ 2.)*
- **72**: damping not in (0, 1).
- **73**: time-history file missing.
- **74 / 75 / 76**: bad first record / bad last record / first > last.
- **77 / 78**: multiplication factor and max value both zero / both non-zero.
- **Warning 9**: NFFT not a power of 2.

Additional runtime checks (implementation):
- NFFT ≤ 32,768;
- DF consistency with FILE8;
- Type of Analysis matches FILE8;
- interpolation option in 0…6;
- phase adjustment in {0, 1};
- S ≥ 0;
- warn if S > 0 with option 6, or S > 0 in a coherent run;
- SRSS file list consistent.

---

## Part C — Verification tests (recommended for the Python implementation)

| # | Module | Test | Expected |
|---|---|---|---|
| C1 | ANALYS | **Free-field reproduction:** embedded model whose "structure" has exactly the excavated-soil properties (C^s = C^e) | U = U'_f at all interaction DOFs, so the ATF equals the SITE free-field TF and is 1 at the control point (no SSI effect) |
| C2 | ANALYS | **Surface rigid foundation global impedance** (option 2) on a uniform half-space, circular footing of radius a: low-frequency limits | Re K_vert → 4Ga/(1−ν), K_horiz → 8Ga/(2−ν), K_rock → 8Ga³/(3(1−ν)), K_tors → 16Ga³/3, within mesh-discretization tolerance. Frequency trends follow Veletsos/Luco–Westmann. |
| C3 | ANALYS | **Restart equivalence:** Mode 1 with save restart, then Mode 3 with the identical FILE1, and Mode 2 with the identical FILE4 | FILE8 identical to round-off |
| C4 | ANALYS | **Simultaneous cases:** coherent run with simul = 1 | FILE8X equals a single-case run with FILE1 = FILE1X (same for Y and Z) |
| C5 | ANALYS | **Coordinate angle:** x' input at angle a | ATF_x = cos a · ATF(x-input, a = 0), ATF_y = sin a · …, for a symmetric model |
| C6 | ANALYS | **FFL vs FFM** for a surface rigid foundation on rock | Identical results (as stated in the manual) |
| C7 | MOTION | **Identity convolution:** H ≡ 1 up to Nyquist | Output equals input to round-off |
| C8 | MOTION | **Interpolation exactness:** sample an analytic SDOF ATF `(ω0² + 2iζω0ω)/(ω0² − ω² + 2iζω0ω)` at the SSI frequencies | Options 0–5 reproduce it exactly (the rational model contains this form). Option 6 error decreases as O(Δf⁴). |
| C9 | MOTION | **Spline on cubic data:** H(f) = complex cubic polynomial | Option 6 is exact |
| C10 | MOTION | **RS harmonic check:** sinusoidal base acceleration at f0, long duration, ζ = 0.05 | SA at f0 ≈ steady-state amplitude ratio `sqrt(1+(2ζ)²)/(2ζ)` × input amplitude; SA(100 Hz) ≈ ZPA |
| C11 | MOTION | **Baseline:** input with a constant acceleration offset added | Corrected velocity has ≈ zero trend; d_c(T_end) ≈ 0 |
| C12 | MOTION | **SRSS:** one mode identical to the coherent FILE8 | p = 1 gives the coherent TF; p = 0 gives |H_coh| |
| C13 | MOTION | **Phase adjustment:** option 2 with S = 1000 vs option 6 (any S) on the same TF | Option 6 gives exactly `|H|` (zero phase). Option 2 / S = 1000 gives phases scaled by ≈ 1e-3, so its RS should be close to the zero-phase RS (set the tolerance in the test plan). Options 0–5 with S = 0 leave H unchanged. |
| C14 | MOTION | **Smoothing invariance:** S = 0 | Bitwise-identical TFI to the no-smoothing path; at computed f, TFI = TFU for every S |

---

## Part D — Recommended-practice summary (from the manual)
- FFL is the default incoherent randomization (EPRI / USNRC). Use FFM only for surface structures.
- Coherent seismic simultaneous cases = 1 (X, Y, Z in one run). Incoherent = number of HOUSE simulations
  (≤ 50). Vibration ≤ 500 load cases.
- Global impedances: use only for surface stick models with rigid mats. For embedded models only with
  FI-FSIN. Use Option A for two-step analysis.
- Interpolation:
  - coherent: options 0–5 with **S = 0**, plus TFU-TFI checks and CRITFREQ-driven frequency additions;
  - incoherent: **option 6** (≥ 200 frequencies) with **S = 0**, or options 0–5 with **S = 10–1000**
    (10–500 gives practically equal ISRS) after sensitivity studies.
- Phase adjustment only for incoherent (stochastic or AS). Never use phase-adjusted or SRSS motions for
  multi-support time-history analysis of secondary systems.
- RS: 0.1–100 Hz, 301 frequencies. Dampings as needed (e.g. 0.02, 0.05).
- Do not save FILE13 for large models. Use RELDISP for relative displacements.
- Turn on Save Complex TF when RELDISP will be run (it needs complex TFI).

---

## Part E — Open questions / ambiguities

1. **Exact interpolation algorithms for options 0–5.** The manual gives only names, plus the 6.4.11 hint
   "two SDOF transfer function model with five parameters". B.6.4–B.6.5 (5-point quadratic/quadratic
   rational fit; stride-4 partitions; averaging rules) is a reconstruction. The window size, the weights
   for option 0, the selection rule for option 3, and how "one/two position shift" is applied (start offset
   vs per-interval shift) are all assumptions. Compare against Tajirian (1981) and the SASSI2000 MOTION
   source/theory if available.
2. **Smoothing-parameter algorithm** (B.6.7) is not documented. The manual gives only the purpose, the
   typical range 10–1000 and the EPRI insensitivity for 10–500. The command reference says 1–1000, but the
   GUI text says 0 for coherent runs. 0 is treated here as "off".
3. **Phase-adjustment algorithm** (B.6.8): "reduces to almost zero the differential phase… depending on the
   smoothing parameter" is implemented as unwrapped-phase scaling by ρ = 1/(1+S) for options 0–5, and as
   full zeroing (ρ = 0) for option 6. The manual says S has no effect for the spline, recommends S = 0 with
   it, and still recommends phase adjustment for incoherent analysis. Two things remain unconfirmed: the
   true dependence on S, and whether the phase is measured from the f = 0 anchor or relative to the
   coherent/free-field phase.
4. **Option 6 spline details:** end conditions, whether an f = 0 anchor is used, whether Re/Im or
   amplitude/phase are splined. The meaning of "bicubic" is unclear (taken here as "cubic on Re and on
   Im").
5. **TF below the first and above the last computed frequency** (B.6.3): H(0) anchor value, and
   zero-filling above the cut-off, are not stated in this section.
6. **Baseline-correction variant:** the dialog says No/With Correction, but MOTION `<bl>` documents
   "0 = time domain, 1 = frequency domain". Also: which polynomial order is used in the Hudson–Housner
   procedure, whether it is applied to every output history or only to FILE13 output, and whether
   velocity/displacement RS can then be requested. Ch. 3 says baseline correction allows "acceleration,
   velocity, or displacement response spectra", but the footnote says displacement RS cannot be requested.
7. **"Acceleration and Velocity R.S."**: is the velocity spectrum relative velocity or pseudo-velocity?
   Is SA absolute acceleration or pseudo-acceleration? The RS frequency spacing (log vs linear), and whether
   "Total Number of Freq. Steps" counts points or intervals, are also not stated (301 points over 0.1–100 Hz
   strongly suggests log spacing, 100 points per decade).
8. **Foundation-vibration Output Control**: the footnotes refer to a displacement/velocity/acceleration
   selection "from the Output Control group box", but the seismic screenshot does not show such a control.
   Its UI form and default are unknown.
9. **FFL/FFM** have no batch-command argument, and the exact form of the incoherency factors applied to the
   load (FFL) or to the motion (FFM) depends on the FILE77 content (HOUSE spec).
10. **Delete Restart Files** checkbox (disabled in the screenshot) is undocumented.
11. **Simultaneous Cases semantics:** command default 0 vs GUI value 1. Whether incoherent simultaneous
    runs also need FILE1X/Y/Z. Mapping of the coherent "1" to three directions vs one. The file numbering for
    incoherent cases is given only by example (FILE8001–003 for simulation 1).
12. **FILE9xxx naming:** Ch. 3 says "only one digit load case number can be appended to FILE9", which
    contradicts `FILE9001…FILE9500` elsewhere. The 3-digit form is used here.
13. **Multiple-excitation control-point coordinates:** the ANALYS dialog says they define the motion's
    application location, but the ME command says "not used in this version".
14. **Global impedance outputs for the 6×6 option:** file layout of FOUNSTIF/FOUNDASH/FOUNDAMP/FOUNIMPD and
    FILE11 contents, the effective damping-ratio definition (Im/(2Re) assumed), and whether rotational
    interaction DOFs (shell/beam interaction nodes) are included in T.
15. **Coordinate Transformation Angle range:** Error 63 says (0, 360). Whether 0 is legal is unclear (the
    screenshot default is 0, and 0 is required for simultaneous cases). [0, 360) is accepted here.
16. **File formats** of TFU/TFI/ACC/RS/RSO/FILE12/FILE13/frames/`Modelname_ACC.bin` are not specified beyond
    names and column semantics. Python formats are defined in B.7. Whether `.ACC` includes a DT header line
    is unknown.
17. **SRSS combination order** (interpolate-then-combine vs combine-then-interpolate) and the handling of
    the coherent-phase FILE8 when frequency lists differ.
18. **Total Duration = 0** meaning (assumed: whole Fourier period). Whether the extracted 20 % extension is
    also applied to RS integration (assumed: RS uses the full record).
19. **Record definition** for First/Last Record in pair-format files (assumed: one record = one data line).
20. **Screenshot example values** (RS 1–8192 Hz, 32 steps, smoothing 1, phase adjustment 1) contradict the
    recommended practice (0.1–100 Hz, 301 steps, S = 0 and phase adjustment 0 for coherent). They are
    treated as examples, not defaults.
21. **Relationship of the MOTION NFFT/DT to SITE**: the manual says they are "the same as SITE". Whether
    MOTION may use a different NFFT (with interpolation onto a different DF grid) is not stated. A strict
    match is required here.
22. **Mode 6 (New Load Vector)** is batch-only and needs vendor information, so its behaviour is undefined.
