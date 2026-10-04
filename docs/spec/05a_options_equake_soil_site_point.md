# 05a: Options Submenu (Model, Write, Check) and Analysis Options for EQUAKE, SOIL, SITE, POINT

**Source:** ACS SASSI V3 User Manual, Section 6.5 "Options Submenu", 6.5.1 to 6.5.4 (EQUAKE, SOIL, SITE and POINT tabs only).
Extracted text lines 4967 to 5546. Printed pages 113 to 129 (PDF pages 115 to 131).
All dialog screenshots on these pages were viewed and transcribed. Some details come from other parts of the manual. Each one is marked "xref" with its section:

- Ch.1.5 limits (text ~L536-560)
- Ch.3.1 module descriptions (~L1727-1895)
- Ch.4 guidelines (~L2329-2540, L2829-2848, L3104-3222)
- Ch.9 command reference: EQUAKE, SOIL, SITE, POINT, WAVE, TOPL, SPRO, SACC, SRS, SSTR, SSAF, SFOU, DAMP, DYNP, CORR, RSIN, RSOUT, ACCIN, ACCOUT, TPSD, THFILE, MOPT, AOPT, AFWRITE, WRITE, NLSOIL, NLSLAYER, RADIUS, GRAVITY
- Ch.10 Check errors and warnings

Tags: **[CORE]** needed for computational correctness · **[IO]** file or input format · **[UI]** user interface or convenience · **[ADV]** advanced option (nonlinear, incoherency, Option A/AA/PRO). Text marked **BACKGROUND (not in manual)** is standard SASSI or SHAKE theory. It is included so the module can be implemented, and it must be verified by the implementer. It is *not* a statement of what the manual says.

---

## 0. Overview of the Options submenu  [UI]

Main window menu bar: `Model  File  Plot  Modules  Options  View  Help`. The **Options** dropdown (screenshot, p.113) contains, in order:

| Menu item | Purpose (xref Ch.6 menu table, L3989-4000) | In this spec |
|---|---|---|
| Model | Sets the model options (dialog "Model Options") | yes, §1 |
| Write | Modifies file-writing options for the AFWRITE/WRITE commands ("Extended Write Options") | yes, §2 |
| Check | Sets the check options ("Check Options") | yes, §3 |
| Analysis | Sets SSI analysis options for the ACS SASSI modules ("Analysis Options", tabbed) | yes, §4 to §8 |
| *(separator)* Windows Settings | Options for the active window | no (other spec) |
| Colors, Font | Plot appearance | no |
| *(separator)* Shader Options, Reset Plot | 3D plotting | no |

Each dialog has an `OK` (also spelled `Ok`) and a `Cancel` button. **OK** commits the values to the active model's in-memory data, exactly as if the matching commands (MOPT, EQUAKE, SOIL, SITE, ...) had been typed. **Cancel** discards the edits. An implementation should map every dialog commit to the equivalent command, so that a `.pre` written by WRITE reproduces the state (see §2).

---

## 1. Options > Model ("Model Options" dialog)  (6.5.1)

The dialog changes the model options of the **active model**. The command equivalent is `MOPT,<incomp>,<matrix>,<mass>,<force>` (xref 9.2.21). Note that the manual's syntax line spells the first argument `<incomp>` while its description uses `<icomp>`.

| Group / control | Type | Default (screenshot) | MOPT arg | Values | Meaning | Tag |
|---|---|---|---|---|---|---|
| Incompatible Modes: "Include Incompatible Modes" / "Suppress Incompatible Modes" | radio | **Suppress** | `<incomp>` | 0 = include, 1 = suppress | Whether 8-node SOLID elements use incompatible (bubble) modes. Never applied to excavated-soil SOLIDs (xref HOUSE remarks, p.130). | [CORE] |
| General Matrix: "Mass Matrix" / "Weight Matrix" | radio | **Mass Matrix** | `<matrix>` | 0 = mass units, 1 = weight units | Units of the mass matrix supplied for GENERAL (general stiffness/mass matrix) elements. With weight units, the mass is entered as weight and must be divided by g. | [CORE] |
| Overwrite Mass | checkbox | **checked** | `<mass>` | 0 = add, 1 = set | A nodal mass defined again on a node that already has one either replaces it (checked, 1) or is added to it (0). | [CORE] |
| Overwrite Force | checkbox | **checked** | `<force>` | 0 = add, 1 = set | Same rule for nodal forces. | [CORE] |

**Implementation notes**
- The mass and force flags act when the mass/force definition commands run. Changing them later does not alter masses or forces that are already stored. This is an interpretation; see Open questions.
- For "Weight Matrix", use the model gravity (HOUSE/GRAVITY) when converting to mass.

---

## 2. Options > Write ("Extended Write Options" dialog)  (6.5.2)  [IO]

**Purpose.** This dialog selects added behaviours for the `WRITE` command (writes the model as a `.pre` command file) and the `AFWRITE` command (writes the module input files). The UI was first designed as a drop-in replacement for the older PREP, with WRITE and AFWRITE unchanged. These options were added later.

**Baseline rule.** If no option in this window is selected, AFWRITE and WRITE must write *exactly* the same files as the previous PREP `.pre` workflow (backward compatibility).

| Control (exact label) | Type | Default | Effect | Tag |
|---|---|---|---|---|
| `MDL command in *.Pre` | checkbox | unchecked | WRITE puts an `MDL,<Model>,<Path>` command (model name and path) into the `.pre`. A `.pre` read later with `INP` then knows its name and path, so SAVE, RESUME and AFWRITE work (xref 5.4, 9.7.20). | [IO] |
| `Extend integer fields (for Models with more than 100000 nodes)` | checkbox | unchecked | Widens the fixed-format integer fields in the module input files written by AFWRITE. This is needed when the SSI model has more than 99,999 nodes. Per xref Ch.1 (L143-146) it affects the **FORCE, HOUSE, MOTION, RELDISP and STRESS** inputs. The `.pre` file is never affected. It applies **only** to the enhanced-solver build that handles up to 300,000 nodes (`IKTR8_300K`; the intro also mentions `IKTR8_650K`). | [IO] |
| `AFWR Command in *.Pre` | checkbox | unchecked | WRITE adds an `AFWRITE` command to the `.pre`, so that `INP` of the file also regenerates the analysis input files (useful for batch work). | [IO] |
| `Simulation Commands` | checkbox | unchecked | WRITE also writes the "simulation" commands that support batch analysis. The manual never defines these; see Open questions. The most likely meaning is the analysis-option and module-run commands (RUNEQUAKE, RUNSOIL, RUNSITE, RUNPOINT, RUNHOUSE, RUNANALYS, ...) for the modules enabled by AOPT. | [IO] |
| Group `Simulation Command Location`: `Write to *-Sim.Pre` / `Write to *.Pre` | radio (greyed out unless Simulation Commands is checked) | **Write to \*-Sim.Pre** | Where the simulation commands go: a separate file `<modelname>-Sim.pre`, or appended to `<modelname>.pre`. | [IO] |

**Related facts** (xref 9.2.3, 9.2.47, 3.2, 10.1):
- `AFWRITE` (no arguments) writes the module input files that the user requests, into the model directory, as `<modelname>.<ext>` with one extension per module. In this range: `.equ` (EQUAKE), `.soi` (SOIL), `.sit` (SITE), `.poi` (POINT). Others: `.hou` (HOUSE), `.anl` (ANALYS), etc. AFWRITE **runs CHECK first**. A module with any **error** (warnings do not count) does not get its input file written. Which modules are processed is set by `AOPT` (the AFWRITE tab, other spec).
- `WRITE,[<file>],[<path>]` writes all model data as instruction lines. The default name is `<modelname>.pre` in the model path. Reloading uses `INP,<file>`.
- The manual (3.2) also mentions running a set of SSI modules automatically through the Options\Write "Running SSI Modules" menu. Batch mode runs each module as `SSI_module_name.exe < SSI_module_name.inp`. The `.inp` file has three lines: `modelname`, `modelname.ext_input`, `modelname_SSI_module_name.out`.

---

## 3. Options > Check ("Check Options" dialog)  (6.5.3)  [UI]

These settings control what the CHECK function writes to the model's **`.err` file** and shows in the Check Errors window.

| Control (exact label) | Type | Default | Meaning |
|---|---|---|---|
| `Show Warnings` | checkbox | checked | Include warning messages in the output |
| `Show Errors` | checkbox | checked | Include error messages in the output |
| `Suppress Error Window` | checkbox | unchecked | Do not pop up the Check Errors window (results still go to `.err`) |
| `Break Check at [N] Messages` | integer | **100** | Most messages of each type (errors and warnings) printed **per module** in `.err` |

Rules:
1. CHECK still checks everything, and reports the **total** number of errors and warnings for every module selected in the command window. Only the *printed* messages are cut off. Messages after the break number are not reported.
2. The options are **global**: they apply to the next check of *any* model.
3. They are **not saved**. They reset to the defaults every time the UI starts. The reimplementation may keep this behaviour for fidelity, or persist the settings as a UI convenience.
4. Since the UI version, the `.err` file has headers that show which module each message belongs to (xref 10.1).

---

## 4. Options > Analysis ("Analysis Options" window): general  (6.5.4)  [UI]

- This sets the analysis options of the **active model**, for the SSI modules that will be used.
- It is a single window with one **tab per module**, in this order: `EQUAKE | SOIL | SITE | POINT | HOUSE | FORCE | ANALYS | MOTION | STRESS | RELDISP | NONLINEAR | AFWRITE`. This spec covers EQUAKE, SOIL, SITE and POINT. The window has `Ok` / `Cancel` buttons.
- **Shared variables:** some fields appear on several tabs and are one stored variable. Editing one updates all of them. From the screenshots and text:

| Shared quantity | EQUAKE tab | SOIL tab | SITE tab | HOUSE tab (other spec) | Command field(s) |
|---|---|---|---|---|---|
| Time step of motion Δt (s) | `Time Step` | `Time Step of Input Motion` (text: "Time Step of Control Motion") | `Time Step Control Motion` | `Time Step of Sesmic Motion` [sic] | `SITE <delt>` (EQUAKE and SOIL commands have no Δt argument) |
| Number of FFT components NFFT | — | `Number of Fourier Components` | `Nr. of Fourier Component` | `Nr. of Fourier Components` | `SITE <nft>` |
| Control point layer | — | `Control Point Layer` | `Control Point Layer` | — | `SITE <cl>` |
| Frequency set number | — | — | `Frequency Set Number` | `Frequency Set Number` | `SITE <freq>` |
| Gravity (free-field, SOIL) | — | `Gravity Accel. (ft/s^2 or m/s^s) (used for free-fixed analysis)` [sic] | — | — | `SOIL <grav>` |
| Gravity (SSI) | — | — | `Gravity Accel. ...` | `Acceleration of Gravity` | `HOUSE <gravity>` / `GRAVITY,<grav>` (inferred; SITE command has no gravity argument) |

  **WARNING:** the SOIL gravity and the SITE (SSI) gravity are **independent**. Changing one does not change the other (stated explicitly in both field descriptions).
- In a command (`.pre`) workflow, every field below maps to a command argument (Ch.9). The note on each command says it is "provided for input files"; interactive users use this window.

---

## 5. EQUAKE tab: "EQUAKE Module Options"  (pp.115-117)

### 5.1 Purpose and methodology (xref 3.1, 6.4.3)  [CORE]
EQUAKE generates acceleration time histories that match given design ground response spectra (RS).

- **Algorithm:** a frequency-domain match by **Levy-Wilkinson (LW)** gives a first approximation. A time-domain match by the **Abrahamson (AB)** algorithm (as in RspMatch) then improves the fit to the target RS.
- **Phasing:** either random phases, or the Fourier phasing of a "seed record" (ASCE 4-16). With seed records, each simulated component keeps the phasing of the seed's X, Y and Z components.
- **Baseline correction:** complex-frequency method (as in FLUSH) plus polynomial corrections in the time domain.
- **Outputs** (text files): `.acc`, `.vel`, `.dis` histories. The user gives the acceleration filename and the velocity/displacement files reuse its base name. Also `.rso` (computed RS), `.psd` (2 columns: freq, PSD), and `.fft` (3 columns: freq, Re, Im, positive frequencies only).
- **PSD:** averaged over ±20 % frequency windows (ASCE 4 / USNRC). PSD and FFT are computed over the **strong-motion duration** only, defined as the time between 5 % and 75 % Arias intensity.
- **Feature parameters reported:** strong-motion duration, V/A, AD/V² (A = PGA, V = PGV, D = PGD). For more than one component, stationary cross-correlation (whole record) and nonstationary cross-correlation (2-s moving window).
- **SRP 3.7.1 Option 1 Approach 2 checks:**
  - Duration ≥ 20 s; otherwise a warning on screen and in the output.
  - Nyquist frequency ≤ 100 Hz; warn if Δt > 0.005 s.
  - At least 100 frequency points per decade.
  - The 5 %-damped RS is never more than 10 % below the target and never more than 30 % above it.
  - No more than 9 adjacent frequency points below the target.
- **Input file:** `.equ`, written by AFWRITE.
- **Spectrum input file format** (e.g. `.rsi` in Demo 1): two columns, frequency and amplitude, for the given damping ratio.
- **Up to three spectra** (three translational directions), indexed by *Spectrum Number* 1..3.

### 5.2 Dialog fields  (screenshot p.115; descriptions p.116-117)

Group **Spectrum Files** (one set per Spectrum Number):

| Label (dialog) | Text name | Type | Default (screenshot) | Command | Meaning / rules | Tag |
|---|---|---|---|---|---|---|
| `Spectrum Number` | Spectrum Number | spin 1..3 | 1 | `<no>` index of RSIN/RSOUT/ACCIN/ACCOUT/TPSD | Selects which spectrum (direction) the file fields below edit. At most 3. | [UI] |
| `Edit` | Edit | button | (disabled when no file) | — | Opens the selected spectrum input file in an editor | [UI] |
| `Spectrum Input File` + `<<` browse | Spectrum Input File | path | e.g. `...\Demo1s\R...` | `RSIN,<no>,<file>` | Target design RS file (2 cols: freq, amplitude) | [IO] |
| `Spectrum Output File` + `<<` | Spectrum Output File | path | `...\Demo1s\XI...` | `RSOUT,<no>,<file>` | Output RS of the simulated motion (`.rso`) | [IO] |
| `Acceleration Output File` + `<<` | Acceleration Output File | path | `...\Demo1s\R...` | `ACCOUT,<no>,<file>` | Simulated acceleration file (`.acc`; `.vel`, `.dis` share its name) | [IO] |

Group **Optional Spectrum Files** (the text calls it "Optional Acceleration Files"):

| Label | Type | Default | Command | Meaning | Tag |
|---|---|---|---|---|---|
| `Accel. Record` | checkbox | unchecked | `EQUAKE <accopt>` | Use "seed record" acceleration input files. The simulated motion keeps the seed's Fourier phasing. | [CORE] |
| `External Accel` | checkbox | unchecked | `EQUAKE <accopt>` (see Open questions) | No simulation. Compute RS (and PSD/FFT) of an **external** acceleration history, which must have the same time step and number of data points. | [CORE] |
| `Acceleration Input File` + `<<` | path | `...\Demo1s` | `ACCIN,<no>,<file>` | Acceleration input file used by whichever option above is selected | [IO] |

Group **Target PSD**:

| Label | Type | Default | Command | Meaning | Tag |
|---|---|---|---|---|---|
| `Use Target PSD` | checkbox | (checked in example) | `EQUAKE ... ,[tpsd]` (0 disabled, 1 enabled; default 0) | Check the generated motion against a target PSD | [CORE] |
| `PSD File` + `<<` | path | blank | `TPSD,<num>,<file>` (`<num>` = spectrum number) | Target PSD text file: 2 columns, frequency and PSD amplitude. Units are **cm²/s³** (IS) or **in²/s³** (BS), depending on the gravity unit system. The text prints "cm/sec^3"; the PSD output description uses cm²/s³. | [IO] |

Ungrouped fields:

| Label (dialog) | Text name | Type | Default | Command arg | Meaning / validation | Tag |
|---|---|---|---|---|---|---|
| `Number of Frequencies` | Number of Frequencies | int | 24 | `<nrfreq>` | Number of frequency points that define the design spectrum. **Must equal the number of records in the spectrum input file** (Error 89). Must be > 0 (Error 91). | [CORE] |
| `Inital Random SEED` [sic] | Initial Random Number | int, 5 digits | 11975 | `<rand>` | Initial seed of the random-phase generator. Must be > 0 (Error 90). | [CORE] |
| `Damping Value` | Damping Value | real | 0.05 | `<damp>` | Damping **ratio** (decimal) of the input design spectrum | [CORE] |
| `Time Step` | Time Step | real (s) | 0.005 | shared Δt (`SITE <delt>`) | Time step of the generated history. Should equal the Δt used by the other SSI modules (shared). | [CORE] |
| `Total Duration` | Total Duration | real (s) | 15 | `<dur>` | Total motion duration. Must be > 0 (Error 92). Warn if < 20 s (SRP). | [CORE] |
| `Number of SEEDs` | Number of Random SEEDs | int | 0 | `<seeds>` ("number of random seeds per Acceleration file") | Number of random-seed trials. EQUAKE keeps the history with the best RS fit. The PSD criterion is **not** checked for these trials; only the `.psd` file is produced. The text says it "should be a non-zero integer" but the default shown is 0 (see Open questions). | [CORE] |
| `Correlated` | Correlated | checkbox | unchecked | `<corr>` 0/1 | Generate correlated X and Y components. When checked, enter time/correlation pairs in the grid. | [CORE] |
| Grid `Correlation` (cols `Time`, `Corr.`; rows 1..n) | — | table | rows: (0,0),(1,0),(3,0),(5,0),(6,0),(8,0),(9,0),... | `CORR,<no>,<time>,<val>` per row | Time-varying target correlation coefficient between X and Y: pairs (time s, coefficient). At least one pair is needed when Correlated (Error 93). A coefficient > 1 is illegal (Error 94; implement as \|ρ\| ≤ 1). | [CORE] |
| `Spectra Title` | Spectra Title | text | "Title for spectra" | `EQTIT,<title>` | Title | [IO] |

Command form: `EQUAKE,<accopt>,<nrfreq>,<rand>,<damp>,<dur>,<corr>,<seeds>,[tpsd]`.

### 5.3 Units WARNING (verbatim intent)  [CORE]
- EQUAKE takes the **gravity acceleration units** to set the units of velocity, displacement and PSD.
- The gravity value is defined in the **SOIL, SITE or HOUSE** window. It must be in m/s² (IS) or ft/s² (BS).
- The simulated **acceleration** history is in **g**. The acceleration PSD amplitude is in **cm²/s³** (IS) or **in²/s³** (BS) (xref L1800-1802). The text here prints "cm/sec^2 / ft/sec^2", which is inconsistent.
- PGA/PGV/PGD-based quantities in the output use the units printed in **NUREG/CR-6728 Tables 3-5 and 3-6**.

### 5.4 EQUAKE check rules (Ch.10)  [CORE]
| Error | Condition |
|---|---|
| 84 | All three spectrum input files are blank |
| 85 | Spectrum input file *i* does not exist |
| 86 | Spectrum output file *i* is blank (for a defined spectrum) |
| 87 | Acceleration output file *i* is blank |
| 88 | Acceleration input file *i* is blank while the acceleration-input option is on |
| 89 | Records in spectrum input file *i* ≠ Number of Frequencies |
| 90 | Initial random number ≤ 0 |
| 91 | Number of frequencies ≤ 0 |
| 92 | Total duration ≤ 0 |
| 93 | Correlated selected but no correlation factors |
| 94 | A correlation factor > 1 |
| 1 | Gravity ≤ 0 (the global gravity check also governs EQUAKE units) |

Limit (xref 1.5): at most **32,768 time steps** per simulated history.

---

## 6. SOIL tab: "SOIL Module Analysis Options"  (pp.118-121)

### 6.1 Purpose and methodology (xref 3.1, 6.4.4, 1.5)  [CORE]
- SOIL runs a **1D free-field site response** for **vertically propagating S-waves** (or P-waves if Input Direction = 1).
- Default model: the **Seed-Idriss iterative equivalent-linear** model, implemented as SHAKE (submodule "SOIL-EQL").
- Optional model: a **nonlinear time-domain** hyperbolic hysteretic model, DEEPSOIL-like, without water-table effects (submodule "SOIL-NON") [ADV].
- **Input file:** `.soi`, written by AFWRITE.
- **Outputs:**
  - Output listing (with an input echo).
  - Text `.TH` time histories: `ACCxxx` (accelerations), `SNxxx` (strains), `SSxxx` (stresses). `xxx` is the free-field layer number, counted from the ground surface downward.
  - **FILE73**: soil material curves, used by STRESS for nonlinear SSI.
  - **FILE88**: iterated strain-compatible (effective) soil properties. SITE uses it when the Non-Linear Soil option is selected. The UI also passes them to the HOUSE excavated-soil layers (xref Step 5, L2473-2484).
- **Units:**
  - BS: ft (thickness), ksf (shear modulus), kcf (unit weight; the gravity text also says lb/ft³), ft/s (velocities).
  - IS: m, kN/m², kN/m³, m/s.
  - **WARNING:** tonnes with metres is not permitted for IS.
- **Limits:**
  - ≤ 32,768 time steps.
  - ≤ 100 soil material curves.
  - ≤ 11 data points per soil curve.
  - ≤ 200 soil layers.
  - ≤ 15 *distinct* dynamic soil properties assigned to sublayers (Error 96).
- **WARNING (nonlinear time domain):** this option is not for direct use in nuclear licensing site response unless fully justified. It is a benchmarking tool for soft soils. It is valid **only** for nonlinear convolution with the input motion defined as the compatible motion **at bedrock**. It was validated against DEEPSOIL in V&V Problem 49.

### 6.2 Dialog fields: group "Input Motion"

| Label (dialog) | Type | Default (screenshot) | Command arg | Meaning / rules | Tag |
|---|---|---|---|---|---|
| `Number of Fourier Components` | int | 4096 | shared `SITE <nft>` | Number of values in the FFT. Same variable as on the SITE tab. Should be a power of 2 (Warning 9: rounded to the nearest power of 2). Must be ≥ Number of Values (implementer: pad with zeros). | [CORE] |
| `Time Step of Input Motion` | real (s) | 0.005 | shared `SITE <delt>` | Time step of the control motion. Same variable as SITE. | [CORE] |
| `Number of Values` | int | 3000 | `SOIL <nrval>` | Number of acceleration values read from the time-history file. Must be > 0 (Error 100). | [IO] |
| `Multiplication Factor` | real | 0 | (see Open questions; MOTION has `<mult>`) | Scale factor for the history. Use **only** if Max Value is blank or 0. | [CORE] |
| `Max Value for Time History` (text: "Max. Value for Time History (in g's)") | real (g) | 0.1 | (see Open questions; MOTION has `<max>`) | Scale the history so its peak \|a\| equals this value (g). Use **only** if Multiplication Factor is blank or 0. Exactly one of the two must be non-zero (Errors 77, 78). | [CORE] |
| `Gravity Accel. (ft/s^2 or m/s^s) (used for free-fixed analysis)` | real | 32.2 | `SOIL <grav>` | Gravity for the **free-field (SOIL)** analysis. ft/s² for BS (unit weight then in lb/ft³ per text, or kcf per units note) or m/s² for IS (kN/m³). Independent of the SITE (SSI) gravity. Must be > 0 (Error 1). | [CORE] |
| `Number of Header Lines` | int | 0 | `SOIL <header>` | Header lines skipped at the start of the history file. Must be ≥ 0 (Error 103). | [IO] |
| `Input Direction` | int flag | 0 | (no documented command argument; see Open questions) | **0 = horizontal** (wave velocity Vs, S-wave damping). **1 = vertical** (velocity **Vp**, P-wave damping). For vertical, **no iterations** should be run (set Number of Iterations = 0). Check the SOIL echo. | [CORE] |
| `Control Point Layer` | int | 1 | shared `SITE <cl>` | Layer number of the control point. The motion is applied at the **top** of this layer (1 = ground surface). Same variable as SITE. Must be a valid layer (Error 104). | [CORE] |
| `File` + `<<` | path | `C:\Users\Owner\Desktop\...` | `THFILE,<file>` (title via `THTIT,<title>`) | Full path of the acceleration time-history input file (in g). View with Plot/Time History. Must exist (Error 73). | [IO] |
| `Assign as Outcrop Motion` | checkbox | **checked** | `SOIL <outcrop>` 0/1 | Checked: the input is an **outcrop** motion. Unchecked: an **in-column (within)** motion at the control point. | [CORE] |

### 6.3 Group "Iteration Parameters" (text: "Equivalent-Linear Soil Behavior")

| Label | Type | Default | Command arg | Meaning | Tag |
|---|---|---|---|---|---|
| `Save Strain-Compatible Soil Properties` | checkbox | checked | `SOIL <save>` (0 skip, 1 save) | Save the iterated properties in **FILE88** for the SSI analysis | [IO] |
| `Number of Iterations` | int | 8 | `SOIL <iter>` | Number of equivalent-linear iterations. **Recommended 8.** Must be ≥ 0 (Error 105). Use 0 for vertical input. | [CORE] |
| `Equiv. Uniform / Max Strain` | real | 0.6 | `SOIL <ratio>` | Ratio R_γ = effective (equivalent uniform) strain / maximum strain. The same value for all layers. **Recommended 0.60 to 0.70.** Must lie in (0,1) (Error 106: "not between 0 and 1"). | [CORE] |

SOIL command form: `SOIL,<nrval>,<grav>,<header>,<outcrop>,<save>,<iter>,<ratio>,<gravmult>,<cof>`.
- `<cof>` is the cut-off frequency. It is **disabled** and is always set to the Nyquist frequency 1/(2Δt) for accuracy. Error 101 if < 0.
- Older PREP versions of the SOIL command end with a format flag in parentheses (e.g. `(…)`). The reader must accept and ignore it, and warn (see Error 102 "Illegal Reading Format").

### 6.4 Group "Soil Profile" (per sublayer; the following groups also apply to the active layer)

| Label | Type | Default | Command | Meaning | Tag |
|---|---|---|---|---|---|
| `Layer Number` | spin int | 1 | `SPRO,<layer>,...` and the first arg of SACC/SRS/SSTR/SSAF/SFOU | Selects the **active** SOIL sublayer (1 = top). All controls below edit this layer. | [UI] |
| `Property Number` | int | 2 | `SPRO <prop>` | Soil layer property number (defined by the `L` command: thickness, unit weight, Vp, Vs, P- and S-damping) for this sublayer | [CORE] |
| `Dynamic Soil Property` + `...` | text + button | "Sand" | `SPRO <dynprop>` (label) | Label of the strain-dependent curves (normalized shear modulus G/Gmax and damping vs effective shear strain) for this layer. `...` opens the "Select Dynamic Soil Property" dialog (xref 7.2.6) to choose, create (New), rename (Edit) or delete a property and edit its table. | [CORE] |

**Dynamic soil property data** (xref DYNP 9.2.10, Fig. p.206)  [CORE][IO]
- `DYNP,<no>,<sg>,<g>,<sd>,<d>,<label>` sets the `<no>`-th point of two curves for property `<label>`: shear strain `<sg>` vs G/Gmax `<g>`, and shear strain `<sd>` vs damping `<d>`.
- The editor table has columns `Strain | Mod. Red. | Strain | Damp`, plus a `Title` field.
- Units, from the screenshot (Rock example: strains 0.0001...; damping 0.4, 0.8, 1.5, 3, 4.6, 6; plot axis "Shear Strain %"): **strain in percent**, **damping in percent**, G/Gmax dimensionless in [0,1].
- Checks:
  - Error 97: no modulus curve.
  - Error 98: a G/Gmax value < 0 or > 1.
  - Error 99: no damping curve.
  - Error 95: no dynamic properties assigned to sublayers.

### 6.5 Output-request groups (per active layer)

| Group | Controls (exact labels) | Defaults | Command | Semantics | Tag |
|---|---|---|---|---|---|
| `Accelerations` | radio `No Computation` / `Compute  Maximum` / `Compute Maximum _Time History`; checkbox `Outcropping` | Compute Maximum; Outcropping off | `SACC,<layer>,<opt>,<outcrop>` with opt 0 = no computation, 1 = compute maximum, 2 = compute maximum and save time history | Acceleration at the **top** of the layer. Outcropping = report the outcrop motion (2 × upgoing wave) instead of the within motion. The saved history is `ACCxxx.TH`. | [IO] |
| `Response Spectrum` | `Save Response Spectrum`; `Outcropping`; `Multiplier for Acceleration of Gravity` [1]; `Damping Ratios` [0.02,0.05] | Save on; Outcrop off; 1; "0.02,0.05" | `SRS,<layer>,<save>,<outcrop>`; `SOIL <gravmult>`; `DAMP,<d1>,…,<d10>` | ARS at the top of the layer, for each damping ratio in the list. Damping ratios are **decimal** fractions (0.02 = 2 %), at most 10 (DAMP; `<d1>`=0 deletes the list). Error 107: no damping ratios defined. `gravmult` must be > 0 (Error 108). It scales the RS output units (see Open questions). | [IO] |
| `Stresses _Strains` | `Compute Stresses`; `Save Stress Time History`; `Compute Strains`; `Save Strain Time History` | on; off; on; off | `SSTR,<layer>,<opt1>,<opt2>,<opt3>,<opt4>` | Stresses and strains at the **centre (mid-depth)** of the layer: compute, and optionally save the histories (`SSxxx.TH`, `SNxxx.TH`) | [IO] |
| `Spectral Amplificarion Factor` [sic] | `Save Spectral Amplification Factor`; `Outcropping of Second Layer`; `Outcropping of First Layer`; `Second Layer Number` [0]; `Frequency Step` [0]; `Title` | all off/0/blank | `SSAF,<layer>,<save>,<outcrop1>,<outcrop2>,<layer2>,<freqstep>,<title>` | Spectral amplification factor: ratio of the response spectrum at this layer to that at `<layer2>`, each optionally as outcrop, sampled at `<freqstep>`. Error 109: illegal second layer number. Error 110: frequency step ≤ 0 (when enabled). | [IO] |
| `Fourier Spectrum` | `Compute Fourier Spectrum`; `Save to File`; `Outcropping`; `Nr. of Smoothings` [0]; `Nr. of Values to be Saved` [0] | off/0 | `SFOU,<layer>,<out>,<save>,<outcrop>,<smooth>,<nrval>` | **Do not use** ("Not usable in this version"). Compute Fourier spectra with EQUAKE instead. The UI may grey it out. Errors 111/112 if negative. | [UI] |

### 6.6 Group "Nonlinear Soil" (time-domain hyperbolic model)  [ADV]
Command forms (xref 9.17.19 to 9.17.20):
- Global options: `NLSOIL,<Opt>,<NSTimeSunInc>,<DispConv>,<ForceConv>,<EqualIt>,<BedInt>,<NLDampType>,<MMmult>,<SMmult>`
- Per layer: `NLSLAYER,<Num>,[curvefit],[B],[S],[refStrain],[Vis]` (all arguments default 0). Error 125: fewer NL soil layers defined than soil layers.

| Label (dialog) | Text name | Default | Command arg | Meaning |
|---|---|---|---|---|
| `Nonlinear Time Domain` | Nonlinear Time Domain | off | `NLSOIL <Opt>` 0/1 | Use the hyperbolic model in the time domain (SOIL-NON) instead of equivalent-linear |
| `Subincrements per Timestep` | Subincrements Per Timestep | 0 | `<NSTimeSunInc>` | Number of sub-increments within each input time step |
| `Dispacment Convergence Error` [sic] | Displacement Convergence Error | 0 | `<DispConv>` | Allowed error norm of displacements for convergence |
| `Force Convergence Error` | Force Convergence Error | 0 | `<ForceConv>` | Allowed error norm of forces (residuals) |
| `Equillibrium Iterations` [sic] | Equilibrium Iterations | 0 | `<EqualIt>` | Most Newton-Raphson equilibrium iterations per step before the program stops |
| `Bedrock Interface`: `Rigid` / `Viscoelastic` | Bedrock interface | Rigid | `<BedInt>` 0 = Rigid, 1 = Viscoelastic (**disabled**, "not applicable to this version") | 0: the bottom layer is rigidly tied to the half-space, so both move with the input acceleration. 1: a viscous dashpot with half-space properties (not available). |
| `Damping Type (1,2,3)` | Damping Type | 0 | `<NLDampType>` | Small-strain damping: **1** frequency independent, **2** visco-elastic (needs Viscosity), **3** Rayleigh (needs the two multipliers) |
| `Mass Matrix Mult.` | Mass Matrix Mult. | 0 | `<MMmult>` | Rayleigh mass multiplier α (C = αM + βK) |
| `Stiff Matirix Mult.` [sic] | Stiff Matrix Mult. | 0 | `<SMmult>` | Rayleigh stiffness multiplier β |
| `Curve fit Hyprbolic Parameters` [sic] | Curve Fit Hyperbolic Parameters | off | `NLSLAYER [curvefit]` | **0** = parameters entered by the user; **1** = parameters fitted to the layer's G/Gmax curve by error minimization. With curve fit, Beta, S exponent, Reference Strain and Viscosity are **ignored**. |
| `Beta` | Beta | 0 | `[B]` | β of the hyperbolic model for the layer (typically about 1 to 1.4) |
| `S exponent` | S exponent | 0 | `[S]` | s exponent (typically about 0.7 to 0.9) |
| `Reference Strain` | Reference Strain ("yield shear strain") | 0 | `[refStrain]` | γ_r (typically about 0.03 %, so in percent) |
| `Viscosity` | Viscosity | 0 | `[Vis]` | Soil viscosity η |

**Hyperbolic constitutive law** (xref 9.17.19, equation image verified, PDF p.323):

$$\tau = \frac{G_0\,\gamma}{1+\beta\left(\gamma/\gamma_r\right)^{s}} + \eta\,\dot\gamma$$

G₀ = elastic (small-strain) shear modulus, γ = shear strain, γ_r = reference strain, γ̇ = strain rate.

BACKGROUND (not in manual, DEEPSOIL practice):
- Unloading and reloading follow the extended Masing rules, giving G/Gmax(γ) = 1/(1+β(γ/γ_r)^s).
- The curve-fit option minimizes the error between this expression and the user G/Gmax curve over (β, s, γ_r).
- The integration uses Newmark-β (average acceleration) with Newton-Raphson equilibrium iterations, and the time step is split into sub-increments.

### 6.7 Equivalent-linear algorithm (implementation outline)  [CORE]
Manual-defined pieces: vertically propagating waves, a layered column above a half-space (or rigid base), the Seed-Idriss/SHAKE iteration, the strain ratio R_γ, the number of iterations, outcrop or within input at the top of the control layer, and the per-layer output locations (acceleration at the top of the layer, stress and strain at the centre, RS at the top).

BACKGROUND (SHAKE methodology, implementer to verify):
1. Build the column from the SPRO list: sublayer *m* → L property → (h_m, γ_unit,m, Vs_m, Vp_m, ξ_s,m, ξ_p,m). ρ_m = γ_unit,m / g_SOIL. G_max,m = ρ_m Vs_m². The last layer is the half-space (see Open questions).
2. Read the motion: skip `<header>` lines and read `<nrval>` values (g). Scale: if mult ≠ 0, a = mult·a_file; otherwise a = a_file·(max / max\|a_file\|). Convert to length units (× g_SOIL) if needed, zero-pad to NFFT, then FFT with Δt.
3. For each iteration: complex modulus G\*_m = G_m(1 + 2iξ_m). SHAKE91 instead uses G_m(1 − 2ξ² + 2iξ√(1−ξ²)); the choice is listed in Open questions. Use Vs\*_m = √(G\*_m/ρ_m) and k\*_m = ω/Vs\*_m. With α\*_m = ρ_m Vs\*_m / (ρ_{m+1} Vs\*_{m+1}):
   - A_{m+1} = ½A_m(1+α\*_m)e^{ik\*_m h_m} + ½B_m(1−α\*_m)e^{−ik\*_m h_m}
   - B_{m+1} = ½A_m(1−α\*_m)e^{ik\*_m h_m} + ½B_m(1+α\*_m)e^{−ik\*_m h_m}
   - Start from A_1 = B_1 (free surface).
4. Control motion at the top of control layer *c*: outcrop motion = 2A_c; within motion = A_c + B_c. Scale all A, B by (input spectrum / control-point transfer value).
5. Strain at the mid-depth of layer *m*: γ(z) = ∂u/∂z = i k\*_m (A_m e^{ik\*z} − B_m e^{−ik\*z}) in displacement form (u = ü/(−ω²)). Take the inverse FFT, γ_max = max\|γ(t)\|, and γ_eff = R_γ·γ_max.
6. Update G/G_max and ξ by interpolating the curves at γ_eff (strain in %; semi-log interpolation in strain is customary). Repeat for `<iter>` iterations. A convergence report (change in G and ξ per iteration) is recommended.
7. For Input Direction = 1 (vertical), use Vp, constrained modulus and P-wave damping, and do not iterate.
8. Outputs per §6.5. FILE88 holds the final strain-compatible Vs (or G) and damping per layer. FILE73 holds the curves.

### 6.8 SOIL check rules (Ch.10)
| Error | Condition |
|---|---|
| 1 | Gravity ≤ 0 |
| 73 | Time-history file missing |
| 77 | Multiplication factor and max value both zero |
| 78 | Both non-zero |
| 95 | No dynamic soil properties assigned |
| 96 | More than 15 dynamic properties used |
| 97 | Property has no modulus curve |
| 98 | Illegal modulus values (not in [0,1]) |
| 99 | Property has no damping curve |
| 100 | Number of values ≤ 0 |
| 101 | Cut-off frequency < 0 |
| 102 | Illegal reading format |
| 103 | Header lines < 0 |
| 104 | Illegal control layer |
| 105 | Iterations < 0 |
| 106 | Strain ratio not in (0,1) |
| 107 | No damping ratios |
| 108 | Gravity multiplier ≤ 0 |
| 109 | Illegal second layer (SAF) |
| 110 | SAF frequency step ≤ 0 |
| 111 | Smoothings < 0 |
| 112 | Values to save < 0 |
| 125 | Not enough NL soil layers |
| Warning 9 | NFFT not a power of 2 (the nearest power of 2 is written) |

---

## 7. SITE tab: "SITE Module Options"  (pp.122-127)

### 7.1 Methodology  [CORE]
SITE has two basic operation modes, enabled by the **Mode 1** and **Mode 2** check boxes. At least one must be on (Error 45).

**Mode 1: form and solve the transmitting-boundary eigenvalue problem (for POINT).**
- Read the soil layer properties. For **each analysis frequency**, form the transmitting-boundary submatrices for the **Rayleigh** (in-plane, 2 DOF per interface) and **Love** (out-of-plane, 1 DOF per interface) cases.
- Solve the two eigenvalue problems to get the eigenvalues (wavenumbers) and eigenvectors (mode shapes) of the layered site.
- **Half-space simulation:** SITE generates a set of sublayers below the user layers, with thickness that **varies with frequency**, and attaches **viscous dashpots** at their base.
- Results go to **FILE2**. FILE2 is input to Mode 2 and to POINT. The eigenproblem of a horizontally layered 3D site is the same as that of a plane-strain model, so **FILE2 serves both 2D and 3D**.

**Half-space details** (xref 4.1 items 8 and 20, 1.5):
- *Variable depth method:* up to **20** extra layers with half-space properties are added, with total thickness **1.5 λ**, where λ = Vs_hs / f (the half-space shear wavelength at frequency f). The fundamental Rayleigh mode essentially vanishes at 1.5λ depth.
- The 1.5λ is split into *n* sublayers whose thickness **increases with depth** (and decreases with frequency).
- *Viscous boundary method:* horizontal and vertical dashpots replace the rigid base under the extended system.
- *n* = 0 means no half-space simulation: **rigid base**.
- Recommended *n* = 10 to 20. **20** gives the best accuracy.
- The total number of user soil layers should be **> 20** for accurate Rayleigh and Love modes, especially for uniform deposits.

**Mode 2: solve the linearized site response problem.**
- Recover the layer properties and the Rayleigh/Love eigen-solutions from FILE2.
- For each wave type present, compute mode shapes and wave numbers in the defined (x′y′z′) coordinate system.
- Scale and superimpose all wave types according to the wave composition and the nature of the control motion.
- Store the result in **FILE1** (free-field motion at the layer interfaces per frequency, for unit control motion), used later by ANALYS for seismic analysis.
- **FILE1 is not generated for foundation-vibration analysis.** If the seismic environment is the same, FILE1 serves both 2D and 3D.
- The actual control-motion time history is *not* needed by SITE. It is used later in MOTION.

Displacement form (xref 2.2, eq. 2.3): u′_f(x) = U′_f exp[i(ωt − kx)].
- U′_f is the vector of interface amplitudes at and below the control point (x = 0).
- k is the complex horizontal wavenumber, which sets propagation (Re) and decay (Im) in x.

**BACKGROUND (thin-layer method of Lysmer & Waas 1972 / Kausel 1981, standard SASSI; not in manual):**
- Each layer is discretized with displacements varying linearly through its thickness. For a layer with complex modulus G\* = G(1+2iξ) (see Open questions on the damping form), complex λ\*, density ρ and thickness h:
  - **Love (SH) case**, DOF u_y at the top and bottom of the layer:
    - A_L = (G\*h/6)[[2,1],[1,2]]
    - G_L = (G\*/h)[[1,−1],[−1,1]]
    - M_L = (ρh/6)[[2,1],[1,2]] (consistent mass)
    - Assemble and solve (k²A + G − ω²M)φ = 0 for N eigenpairs. Keep the root with Re(k) > 0 and Im(k) ≤ 0 for e^{i(ωt−kx)} outgoing, decaying waves.
  - **Rayleigh (P-SV) case**, DOF (u_x, u_z) per interface: (k²A + ikB + G − ω²M)φ = 0. This is quadratic in k and has 2N roots. Kausel's change of variables (scale the vertical components by i) reduces it to a linear eigenproblem in k². The implementer must take the exact A, B, G matrices and sign conventions from Kausel (1981) / Waas (1972).
- Base dashpots (viscous boundary, Lysmer-Kuhlemeyer) per unit area: c_x = c_y = ρ_hs Vs_hs and c_z = ρ_hs Vp_hs. They are added as iω·c to the bottom-interface diagonal of G.
- Body waves (SV, P, SH) at incident angle α: apparent horizontal wavenumber k = ω sin α / V_hs (V = Vs for SV/SH, Vp for P). Vertical incidence (α = 0) gives k = 0, i.e. a 1D column response.
- Each wave field is normalized so that its motion at the control point, in the control direction, equals 1. The total field is then U = Σ_w r_w(f) · U_w / U_w(control, dir). That is why the ratios at each frequency must sum to 1.

### 7.2 Dialog fields  (screenshot p.122)

Group **Operation Mode**:

| Label | Type | Default | Command arg | Meaning | Tag |
|---|---|---|---|---|---|
| `Linear Soil` / `Non-Linear Soil` | radio | **Non-Linear Soil** (example) | not in the SITE command (see Open questions; the command's `<opmode>` is 0 complete solution / 1 data check) | Linear: SITE uses the L-command soil layer properties as entered. Non-Linear: the primary soil nonlinearity is included through the SOIL equivalent-linear results. SITE uses the strain-compatible properties from **FILE88**, and the UI also passes them to the HOUSE excavated-soil layers. With Non-Linear, the SITE and SOIL data must use **consistent units**: BS (kip, ft, ksf, kcf) or IS (kN, m, kN/m², kN/m³). | [CORE] |
| `Mode 1` | checkbox | checked | `SITE <mode1>` 0 skip / 1 write | Enable Mode 1 (FILE2) | [CORE] |
| `Mode 2` | checkbox | checked | `SITE <mode2>` 0 skip / 1 write | Enable Mode 2 (FILE1) | [CORE] |

Group **Mode 1**:

| Label | Text name | Type | Default | Command arg | Meaning / rules | Tag |
|---|---|---|---|---|---|---|
| `Gravity Accel. (ft/s^2 or m/s^s) (used for free-fixed analysis)` (label copied from SOIL) | Acceleration of Gravity | real | 32.2 | model gravity (HOUSE `<gravity>` / `GRAVITY`; inferred) | g for the **SSI** analysis: ft/s² (BS) or m/s² (IS). Converts unit weight to mass density, ρ = γ/g. Independent of the SOIL gravity. | [CORE] |
| `Frequency Step` | Frequency Step | real (Hz) | 0 | `SITE <fstep>` | Δf. Leave blank for a time-history analysis. Required for a single-harmonic analysis. Must be ≥ 0 (Error 48). | [CORE] |
| `Time Step Control Motion` | Time Step of Control Motion | real (s) | 0.005 | `SITE <delt>` (shared) | Δt of the control motion. Must be ≥ 0 (Error 49). | [CORE] |
| `Nr. of Fourier Component` | Nr. of Fourier Components | int | 4096 | `SITE <nft>` (shared) | NFFT. Must be a power of 2; otherwise the UI uses the **closest** power of 2 and warns (Warning 9). Must be ≥ 0 (Error 50). | [CORE] |
| `Frequency Set Number` | Frequency Set Number | int | 1 | `SITE <freq>` | Index of a frequency set built with `FREQ,<ndx>,<f1>..<f10>`. It holds **positive integer frequency numbers**. SITE sorts them ascending and **stops if two are equal**. Analysis frequency f_i = n_i·Δf. The set must exist (Error 44) and not be empty (Error 120). At most 500 analysis frequencies. | [CORE] |
| `Number of Generated Layers` | Number of Generated Layers (used for half-space simulation) | int | 20 | `SITE <nl>` | Number of half-space sublayers. Text says 10 to 20, Check (Error 47) allows 0 or 4 to 20, maximum 20. Blank or 0 = no half-space simulation (rigid base). Otherwise the frequency-dependent sublayers and the base viscous boundary for body-wave radiation are added. | [CORE] |
| `Halfspace Layer` | Half-Space Layer | int | 2 | `SITE <hs>` | **L-command soil layer property number** used for the half-space or bedrock | [CORE] |
| `Top Layers` (multi-line text) | Top Layers | int list | `2,2,2,2,2,2,1,2,2,2` | `TOPL,<l1>,[l2],…,[l200]` (`<l1>`=0 deletes the list) | The free-field layering above the half-space or bedrock, as a list of **L property numbers**, **first = topmost**. Separators: blank, tab, `,`, `;`, Enter. At most **200** entries. Layers with identical properties may reuse the same L number. Must not be empty (Error 46). View with Plot/Layers. SOIL and SITE properties must be consistent for nonlinear SSI. | [CORE] |

**WARNING (Top Layers, embedded models):** the L layer numbers must be the embedment soil layer numbers, and they must equal the far-field layer numbers. Per TOPL/L: the L layers are defined in the same order as in TOPL, and embedment layers are **not repeated** in TOPL. Otherwise the UI may generate wrong L materials for the excavated soil layers in the HOUSE input (`.hou`).

**Frequency step logic**  [CORE]
- *Time-history analysis:* Δt and NFFT must both be defined. Δf may be blank. SITE computes **Δf = 1 / (Δt · NFFT)** and uses it with the frequency numbers.
  - Example: Δt = 0.005 s, NFFT = 4096 gives Δf = 0.048828125 Hz.
  - Nyquist frequency = 100 Hz, at frequency number NFFT/2 = 2048.
- *Single-harmonic analysis:* Δf must be given and is used directly. Δt and NFFT are not used and may be blank.
- Recommended practice (xref Step 9): 40 to 80 frequencies for simple stick models, 100 to 300 for complex FE models.

Group **Mode 2**:

| Label | Type | Default | Command | Meaning | Tag |
|---|---|---|---|---|---|
| `R-,SV-,and P-Waves` / `Sh- amd L-Waves` [sic] | radio | R-,SV-,P- | `SITE <wopt>` 0 = R-, SV- and P-waves; 1 = SH- and L-waves | Wave family. In-plane (x′z′), **2 DOF per node**; or out-of-plane (y′), **1 DOF per node**. Picks which wave tab pages appear. | [CORE] |
| Tab pages `R-Wave`, `SV-Wave`, `P-Wave` (family 0) or `SH-Wave`, `L-Wave` (family 1), scrolled with ◂ ▸ | tabs | R-Wave shown | `WAVE,<type>,<opt>,<ratio1>,<ratio2>,<angle>`, one per wave type (type 1 R, 2 SV, 3 P, 4 SH, 5 L) | One page per wave type | [UI] |
| R-Wave page radio: `No Wave Field` / `Shortest Wavelength` / `Least Decay` | radio | No Wave Field | `WAVE,1,<opt>`: 0 none, 1 shortest wavelength, 2 least decay (R-waves only) | How the Rayleigh mode is chosen from the Mode-1 eigen-solutions. BACKGROUND interpretation: "Shortest Wavelength" = the mode with the largest Re(k) (slowest phase velocity, usually the fundamental mode); "Least Decay" = the mode with the smallest attenuation Im(k) (or \|Im k\|/Re k). | [CORE] |
| SV-Wave page: `No SV-Wave Field` / `SV-Wave Field` | radio | — | `WAVE,2,<opt 0/1>` | Include SV | [CORE] |
| P-Wave page: `No P-Wave Field` / `P-Wave Field` | radio | — | `WAVE,3,<opt 0/1>` | Include P | [CORE] |
| SH-Wave page: `No SH-Wave Field` / `SH-Wave Field` | radio | — | `WAVE,4,<opt 0/1>` | Include SH | [CORE] |
| L-Wave page: `No L-Wave Field` / `L-Wave Field` | radio | — | `WAVE,5,<opt 0/1>` | Include Love waves (BACKGROUND: shortest-wavelength mode) | [CORE] |
| `Wave Ratio 1` (per page) | real | 0 | `WAVE <ratio1>` | Participation ratio of this wave type at Frequency 1 | [CORE] |
| `Wave Ratio 2` (per page) | real | 0 | `WAVE <ratio2>` | Participation ratio at Frequency 2 | [CORE] |
| `Incident Angle` (per body-wave page; not visible on the R-Wave page shown) | real (deg) | 0 | `WAVE <angle>` | Angle between the **propagation direction and the z′ axis**. 0 = vertical propagation. Check range (0,360) (Error 52); implement as [0, 360). | [CORE] |
| `Frequency 1` | real | 1 | `SITE <freq1>` ("lowest for defining wave contribution (= 1)") | First frequency of the ratio curve, **common to all wave types** | [CORE] |
| `Frequency 2` | real | 4000 | `SITE <freq2>` ("highest ... (= > NFFT/2)") | Second frequency of the ratio curve. The command note suggests these are **frequency numbers** (index × Δf), not Hz (see Open questions). | [CORE] |
| `Control Point Layer` | int | 1 | `SITE <cl>` (shared with SOIL) | Layer whose **top** is the control point, where the control motion is specified (1 = surface) | [CORE] |
| `Direction`: `X` / `Y` / `Z` | radio | X | `SITE <cm>` 0 = X, 1 = Y, 2 = Z | Direction of the control motion in **x′y′z′**. ANALYS transforms x′y′z′ to the structure's global xyz using the coordinate transformation angle `ANALYS <ang>`. | [CORE] |

**Wave-ratio rules**  [CORE]
- Several wave types: give each one's participation ratio at **two frequencies** (Frequency 1, Frequency 2) that **cover the whole analysis range**.
- Ratios at intermediate frequencies come from **simple (linear) interpolation**. They need not be given at the exact solution frequencies.
- A single wave type needs only two frequencies (one at the start and one at the end of the range) with ratio = **1**.
- All ratios are positive decimals ≤ 1. The ratios of all participating waves must **sum to 1 at each frequency**.
- Checks:
  - Error 51: all wave fields deselected.
  - Error 53: a frequency ≤ 0.
  - Error 54: a selected wave ratio ≤ 0 or > 1.
  - Error 55: sum ≠ 1 (use a tolerance such as 1e-6).
- Linear interpolation of ratios that sum to 1 at both ends keeps the sum equal to 1 in between.

**Wave types and coordinate system**  [CORE]
- Body waves: P and S. Surface waves are generated when body waves reach the free surface or layer interfaces: R (Rayleigh) and L (Love).
- P moves along the propagation direction. S moves perpendicular to it: SV in the vertical plane, SH horizontal.
- R gives horizontally propagating **elliptical** motion in the vertical plane. L gives horizontal motion perpendicular to the horizontal propagation direction.
- SITE frame: z′ is **vertical up**. x′ lies in the vertical plane of propagation. y′ = z′ × x′ (right-hand rule).
- P, SV and R move in **x′z′**. SH and L move along **y′**.
- The figure on p.125 shows the incident wave from below at angle α to z′, in the vertical x′z′ propagation plane. The arrows mark P, R in plane, SV, R in plane, and SH, L along y′.

**WARNING: "Simultaneous Cases"** (ANALYS option)  [CORE][IO]
If ANALYS "Simultaneous Cases" is selected, run SITE **three times** before ANALYS and copy FILE1 after each run:
- **FILE1X**: SV waves, direction **x′**, angle 0.
- **FILE1Y**: SH waves, direction **y′**, angle 0.
- **FILE1Z**: P waves, direction **z′**, angle 0.

The `.anl` coordinate transformation angle must be **0**. The manual (2.2) also recommends SV for X, SH for Y and P for Z, all vertically propagating, as the standard ASCE/USNRC 3D input.

### 7.3 SITE check rules (Ch.10)
| Error / Warning | Condition |
|---|---|
| Error 1 | Gravity ≤ 0 |
| Error 44 | Frequency set not defined |
| Error 45 | Mode 1 and Mode 2 both off |
| Error 46 | No top layers |
| Error 47 | Half-space layers not 0 and not in 4..20 |
| Error 48 | Δf < 0 |
| Error 49 | Δt < 0 |
| Error 50 | NFFT < 0 |
| Error 51 | All wave fields off |
| Error 52 | Incident angle outside (0,360) |
| Error 53 | Ratio-curve frequency ≤ 0 |
| Error 54 | Ratio ≤ 0 or > 1 |
| Error 55 | Ratio sum ≠ 1 |
| Error 120 | Frequency set empty |
| Errors 19 to 25 | L-layer data: undefined layer, thickness ≤ 0, negative unit weight, negative Vp, negative Vs, negative damping |
| Warning 8 | More than 100 top layers: only the first 100 are written (conflicts with the stated maximum of 200) |
| Warning 9 | NFFT not a power of 2 |

Limits: 200 soil layers, 20 half-space layers, 500 analysis frequencies.

---

## 8. POINT tab: "POINT Module Options"  (pp.127-129)

### 8.1 Methodology  [CORE]
- POINT has two programs: **POINT2** (2D, plane strain) and **POINT3** (3D, axisymmetric).
- It computes the data needed to form the **frequency-dependent free-field flexibility (compliance) matrix** at the interaction nodes, and saves it in **FILE3**.
- It needs **FILE2** from SITE, so SITE must run first.
- **Input file:** `.poi`, written by AFWRITE.
- Per xref 4.1 item 17: the free-field flexibility follows the **point-load solution of a soil column**. The column is modelled with plane-strain elements (2D) or axisymmetric elements (3D). ANALYS later inverts the compliance, K + iD = (f + ig)⁻¹, to get the impedance.

BACKGROUND (SASSI POINT, not in manual):
- A unit harmonic load in each translational direction is spread over a small **central zone** (disk of radius r in 3D, strip in 2D) at each layer interface from the surface down to the last embedment interface.
- The central zone is discretized (axisymmetric or plane-strain elements). Outside it, the far field is represented by the consistent transmitting boundary built from the FILE2 eigen-solutions.
- The resulting interface displacements as functions of horizontal distance are stored in FILE3. ANALYS evaluates them at the actual distances between interaction nodes.
- Interaction nodes lie only on layer interfaces (xref 1.5.1 item 5). With L embedment layers, loads and responses are needed at interfaces 1..L+1.

### 8.2 Dialog fields  (screenshot p.127)

| Label (dialog) | Text name | Type | Default (screenshot) | Command arg | Meaning / rules | Tag |
|---|---|---|---|---|---|---|
| `Operation Mode`: `Solution` / `Data Check` | Operation Mode | radio | **Solution** | `POINT <opmode>` 0 = complete solution, 1 = data check only | Run fully, or check the input only | [CORE] |
| `Number of Embedment Soil Layers` | Last Layer Number in Near Field Zone | int | 0 | `POINT <layer>` ("last layer number in near field zone, equal to the number of the embedment layers") | The deepest layer (counted from the surface in the SITE top-layer list) that holds the structure and any irregular near-field soil zone. Smaller means less data in FILE3 and less ANALYS work. It must be large enough that the **excavated soil never extends below this layer**. **0 or blank** for surface structures with no irregular soil zone. Must be ≥ 0 (Error 56). At most 50 embedment layers. | [CORE] |
| `Point Load Central Zone Radius` | Radius of Central Zone | real (length) | 13.8 | `POINT <rad>` | Radius of the central zone of the 2D or axisymmetric point-load solution. Must be **> 0** (Error 57). It depends on the FE mesh of the foundation or excavation. | [CORE] |

### 8.3 Central-zone radius rules (figures p.128)  [CORE]
For a **uniform** foundation mesh with element size *h*:

| Mesh | Figure geometry | Rule |
|---|---|---|
| 3D rectangular (quadrilateral) mesh | Plan L = 6h, W = 4h, square elements of side h. The circle of radius r is centred on a node and covers the 2×2 element block around it; the hatched corners lie outside the circle. | **r = 0.90 h** |
| 3D triangular mesh (circular foundation, equilateral triangles of side h) | Circular plan of diameter 6h = L. The circle covers the hexagon of 6 triangles around a node. | **r = 0.85 h** |
| 2D (plane strain) | Strip of width L = 6h with square elements. r is shown as one element width. | **r = h** |

- **Non-uniform excavation meshes:** use an average value. The `RADIUS,<Scale>,<FileName>` command (xref 9.7.25) writes the per-element radii of the excavation and their average.
- Uniform meshes are recommended wherever possible, with transition meshes between the structure and the excavated soil so that the excavation volume keeps a regular mesh (Nie et al., 2013).
- Sensitivity studies are always recommended. Run with the **minimum, average and maximum** radius from the excavated mesh sizes, and envelope the responses if they differ significantly (xref 1.5.1 item 9).
- **WARNING:** SSI results should be fairly insensitive to the radius at low and mid frequencies, as long as r is about the size of the interaction-volume elements. At **higher frequencies** the approximation can deviate more from the exact solution.

### 8.4 POINT check rules
| Error | Condition |
|---|---|
| 56 | Last layer number < 0 |
| 57 | Radius ≤ 0 |

---

## 9. Data-flow summary for these modules  [IO]

| Module | Input written by AFWRITE | Reads | Writes |
|---|---|---|---|
| EQUAKE | `<model>.equ` | target RS files (RSIN), optional seed or external acceleration (ACCIN), optional target PSD (TPSD) | `.acc`, `.vel`, `.dis` (ACCOUT name), `.rso` (RSOUT), `.psd`, `.fft`, output listing |
| SOIL | `<model>.soi` | acceleration time history (THFILE), soil curves (DYNP) | output listing, `ACCxxx`/`SNxxx`/`SSxxx` `.TH`, **FILE73** (curves), **FILE88** (strain-compatible properties) |
| SITE | `<model>.sit` | L layers, TOPL, FREQ set, **FILE88** (Non-Linear Soil) | **FILE2** (Mode 1 eigen-solutions), **FILE1** (Mode 2 free-field; copied to FILE1X/Y/Z for simultaneous cases) |
| POINT | `<model>.poi` | **FILE2** | **FILE3** (point-load or flexibility data) |

Run order for linear seismic SSI (xref 3.2): EQUAKE (if needed), SOIL (if needed), SITE, POINT, HOUSE, ANALYS, MOTION, STRESS. For vibration analysis: SITE, POINT, FORCE, HOUSE, ANALYS, ...

Note: HOUSE needs the `.sit` file to be present (xref L1930).

---

## 10. Suggested verification tests (not from the manual; for "verifiable" goal)
1. **Δf arithmetic.** Δt = 0.005, NFFT = 4096 should give Δf = 0.048828125 Hz. The frequency set {1, 2, …} maps to f = nΔf. Also test duplicate detection and the sort.
2. **NFFT rounding.** 3000 should round to 2048 (nearest power of 2, |3000−2048| = 952 < 1096) with Warning 9. Use the same rounding rule in the UI and the writer.
3. **SOIL, uniform damped layer on rigid base, linear (0 iterations), within input at base.** The surface/base amplification should be \|H(f)\| = 1/\|cos(ωH/Vs\*)\|. The first peak is at f₁ = Vs/(4H), with amplitude ≈ 2/(πξ) for small ξ (with G\* = G(1+2iξ)).
4. **SOIL outcrop vs within.** For a uniform elastic half-space, the outcrop control motion is 2A and the surface motion equals the outcrop motion.
5. **SOIL equivalent-linear iteration.** Compare against SHAKE91 or DEEPSOIL (EQL mode) for a published example. The manual's V&V suite uses SHAKE91; V&V Problem 49 covers SOIL-NON against DEEPSOIL.
6. **SITE Rayleigh mode, homogeneous half-space** with 20 generated half-space layers and fine top layers. The phase velocity should approach the Rayleigh velocity c_R ≈ 0.9194 Vs (ν = 0.25), and be independent of frequency.
7. **SITE vertical SV in a uniform column.** The free-field FILE1 amplitudes should match the 1D SH/SV transfer function from test 3.
8. **Wave ratios.** With R at 0.3 → 0.5 and SV at 0.7 → 0.5 between Frequency 1 and 2, the interpolated sum is 1 everywhere.
9. **POINT static limit.** For a uniform vertical pressure q over a disk of radius a on a homogeneous elastic half-space, the centre displacement is w₀ = (1−ν)qa/G. The low-frequency POINT3 surface vertical compliance should approach it.

---

## 11. Open questions / ambiguities (implementer must decide)
1. **Number of Generated Layers range.** The text says 10 to 20 (maximum 20); Check Error 47 accepts 0 or 4 to 20. Proposal: accept 0 or 4 to 20, and warn below 10.
2. **Top Layers maximum.** The text and limits say 200; Warning 8 says more than 100 are truncated to the first 100. Proposal: allow 200 (Version 3 limit) and drop Warning 8, or keep it configurable.
3. **Frequency 1 / Frequency 2 units.** The dialog defaults (1, 4000) and the SITE command note ("= 1", "≥ NFFT/2") suggest **frequency numbers** (multiples of Δf). The field text says "frequency". Proposal: treat them as frequency numbers and convert by × Δf, but show the Hz value in the UI.
4. **SITE Linear/Non-Linear Soil flag.** It is not an argument of `SITE,...` (whose `<opmode>` is complete/data-check). The SITE dialog has no Data Check option. It is unclear where the flag is stored, and whether FILE88 is read by SITE at run time or merged into `.sit` by AFWRITE. Proposal: add a model-level flag (e.g. an extra trailing SITE argument with default 0) and have SITE read FILE88 when it is set.
5. **SITE gravity.** The SITE command has no gravity argument. The SITE gravity is presumably the model or HOUSE gravity (`GRAVITY`, `HOUSE <gravity>`). The SITE dialog label says "used for free-fixed analysis" (a copy of the SOIL label) while the text says it is for SSI. Proposal: SITE gravity = model (HOUSE) gravity, separate from SOIL gravity.
6. **SOIL fields without a documented command argument:** Input Direction, Multiplication Factor, Max Value, NFFT, Δt, Control Point Layer. Proposal: NFFT, Δt and control layer come from SITE (shared). Mult and Max are stored per model (possibly shared with MOTION `<mult>`, `<max>`; decide whether SOIL and MOTION share them). Input Direction needs a new stored field (e.g. an extra SOIL argument).
7. **EQUAKE `<accopt>`** is documented as 0/1, but the dialog has two independent checkboxes (Accel. Record and External Accel). Proposal: 0 = none, 1 = seed record, 2 = external acceleration (RS only); verify against `.equ` writer needs.
8. **EQUAKE Number of SEEDs.** The text requires a non-zero integer but the screenshot default is 0. Proposal: treat 0 as 1 (a single trial).
9. **EQUAKE PSD units.** The target PSD text says "cm/sec^3" and "inch^2/sec^3". The WARNING says acceleration PSD is in "cm/sec^2 or ft/sec^2". Ch.3 says cm²/s³ or in²/s³. Proposal: cm²/s³ (IS) and in²/s³ (BS).
10. **EQUAKE correlation grid.** Semantics assumed: piecewise-linear target correlation ρ(t) between X and Y; time in s. Valid range assumed \|ρ\| ≤ 1 (the manual only flags > 1). How the correlation is imposed (e.g. Cholesky mixing of the phase or white-noise processes per window) is not specified.
11. **Damping definition** for complex moduli in SOIL and SITE (G(1+2iξ) vs SHAKE91 G(1−2ξ²+2iξ√(1−ξ²))) is not given in this section. It must be consistent across SOIL, SITE, POINT and HOUSE.
12. **Sublayer thickness distribution** of the generated half-space layers. Only "total 1.5λ, increasing with depth, up to 20 layers" is given. Proposal: a geometric progression, with the first sublayer thickness tied to the bottom user-layer thickness. Document the choice and test it against the homogeneous half-space Rayleigh velocity.
13. **R-wave mode selection** ("Shortest Wavelength" vs "Least Decay"). The exact criterion is not defined. Proposal: max Re(k) vs min \|Im(k)\| (or min \|Im k\|/Re k) among modes with Re(k) > 0. The same question applies to which Love mode is used.
14. **SOIL half-space.** It is not stated whether the last SPRO sublayer is the half-space (SHAKE convention), or whether SITE's `Halfspace Layer` is reused. Proposal: last SPRO entry = half-space (elastic, damped), unless the Nonlinear rigid bedrock option applies.
15. **SOIL dynamic-curve units.** The screenshot shows strain and damping in **percent** in the property table, while the RS Damping Ratios are **decimal**. Curve interpolation (linear vs semi-log) is not specified. Proposal: semi-log interpolation in strain, constant extrapolation outside the table.
16. **Too many dynamic properties.** The limit is 15 used (Error 96) vs 100 curves (limits list). There are 11 points per curve (limits list) vs a typical 20 in SHAKE. Proposal: enforce 100 defined, 15 used, and allow more than 11 points with a warning.
17. **`Multiplier for Acceleration of Gravity`** (`SOIL <gravmult>`, RS group) has an unstated role. Proposal: RS output = computed spectral acceleration in g × gravmult (unit scaling of the RS output).
18. **Extended Write "Simulation Commands"** is undefined: analysis-option commands, RUNxxx module commands, or both. The `-Sim.pre` file content and order are also unspecified. The extended integer field widths are not given. Original fixed formats are typically I5 for node numbers; the extension could be I6 to I7.
19. **Check Options "Break Check at".** The limit is described as "per module"; whether it applies to each message type separately or to the combined count is unclear. Proposal: each type separately, per module.
20. **Model Options mass/force flags.** It is unclear whether Overwrite Mass / Overwrite Force act only on later definitions. Proposal: they affect only subsequent mass/force commands.
21. **POINT2 vs POINT3 selection.** It presumably follows HOUSE "Dimension of Analysis" (2D/3D). This is not stated in this section.
22. **Incident Angle field location.** It is described in the text but not visible on the R-Wave page screenshot. Proposal: show it on the SV, P and SH pages only. R and L are surface waves and the angle does not apply, though the WAVE command still carries it (write 0).
23. **Last Layer Number semantics.** The dialog label says "Number of Embedment Soil Layers" and the text says "Last Layer Number in Near Field Zone". These are equal when counting from the surface. Confirm that interface L+1 (the base of the embedment) is included in FILE3.
