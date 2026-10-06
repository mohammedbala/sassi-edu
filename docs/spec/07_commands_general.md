# 07 — Command Reference, Part 1: §9.1–9.2 "General Commands Compatible with Previous PREP Versions"

**Source:** ACS SASSI V3 User Manual, Chapter 9 intro, §9.1 and §9.2.1–9.2.50.
Text lines 8691–9716 of `reference/acs-sassi.txt`, printed pages 216–237. With the Read tool the PDF page is printed + 2, so PDF pp. 218–239.

**Cross-checked against** (for argument meaning, dialog fields, limits and check rules):
- §6.5.1–6.5.4, the Options ▸ Model / Write / Check / Analysis dialogs, printed pp. 113–184. I viewed the screenshots.
- §1.2 (size limits), §1.6 (input syntax), §1.7 (data checking), §5.4–5.6 and §5.9 (`.pre` handling, MDL, macros, comments), §6.1.4 (Input).
- §7.2.6, the Select Dynamic Soil Property table (units for DYNP).
- §9.7.16/17/20/28 (GRAVITY, GROUNDELEV, MDL, TPSD) and Chapter 10 (Check errors and warnings).

Companion specs: `05a…05d` (dialog-level specs of the Options ▸ Analysis tabs), `02` (module file flow) and `03` (guidelines). **This file is the authoritative definition of the `.pre` command syntax for these 50 commands.** Where 05a–05d propose an interpretation, this file uses the same one.

**Tags:**
- **[CORE]** needed for computational correctness.
- **[IO]** file or input format.
- **[UI]** user interface or convenience.
- **[ADV]** advanced options (incoherency, nonlinear, Option A/AA/PRO).

Statements marked *(inferred)* are not written in the manual. They are reconstructed from dialog screenshots, Check-error texts or engineering logic, and they need a decision (see the Open questions).

---

## 0. Where these commands sit in the program

- The ACS SASSI UI (successor of the old **PREP** pre-processor) keeps a model state: nodes, elements, materials, layers, and so on. It also keeps an **"Analysis Options"** state for every SSI module (EQUAKE, SOIL, SITE, POINT, HOUSE, FORCE, ANALYS, MOTION, STRESS, RELDISP, NONLINEAR, AFWRITE).
- Almost every §9.2 command ends with the boilerplate note *"This instruction is provided for input files…"*. **The command is the text serialization of an Options ▸ Analysis (or Options ▸ Model) dialog setting.** The dialog and the command must edit the *same* state variables. [CORE][IO]
- **WRITE** serializes the whole model, analysis options included, into `modelname.pre`. **INP** replays a `.pre`. **AFWRITE** turns the state into the module input decks (`.equ .soi .sit .poi .hou .frc .anl .mot .str .rdi .eql`). **CHECK** validates the state.
- **WARNING (manual, p. 218):** the dialogs are the recommended way to set these options. After that, WRITE dumps all inputs to `.pre` in a logical order, for review. [UI]

---

## 1. The `.pre` command language: lexical and syntax rules  [IO][CORE]

### 1.1 Notation used in the manual (Backus-Naur style)
| Notation | Meaning |
|---|---|
| `<value>` | Replace with an actual number or string. |
| `[item]` or `[<item>]` | Optional. It may be omitted, usually trailing. |
| anything else (keywords, commas) | Typed literally. |
| `…` | Repetition up to the stated maximum (e.g. `<a1>…<a100>`). |

### 1.2 Instruction form
- General form: `Keyword, p1, p2, …, pn`. Free format, with **comma** as the parameter delimiter. Spaces after commas occur in the manual (e.g. `EQUAKE, <accopt>,…`, `RSOUT, <no>,<file>`), so the parser **must trim whitespace** around tokens.
- **Command names are not case-sensitive.** File names, titles and dynamic-property labels keep their case. `SOILPROPPLOT` states that a soil-property name is case-sensitive.
- Data may come in any order. Instruction line, dialogs, menus and data files can be mixed (§1.6).
- **Comment lines start with `*`** (all `.pre` examples in §5.6 and §5.9). *(inferred: a `*` in column 1, after optional leading blanks, comments out the whole line.)* Blank lines are ignored.
- **Unknown command:** the UI prints `"<X> Command not found"` and **keeps processing** the file. Only later commands that depend on the missing one are affected (§6.1.4). [UI]
- Text substitution comes before the line is parsed as a command:
  - Macro arguments `$n$` (§5.6).
  - Variables `@NAME`, `@NAME[i]`, `@NAME++`, `@NAME+k`, `@NAME=k`, and `#` in FOREACH (§5.9).
  - Limit: 3000 characters per macro command line.

### 1.3 Command-name abbreviation rule  [IO]
- Old PREP read **only the first 4 letters** of an instruction name. The new UI reads the whole name.
- For names longer than 4 characters, the manual **underlines** an accepted abbreviation. The **full name** and the **exact underlined part** are both valid. **No other variant is valid**: no other prefix length, and nothing longer than the name.
- Names of 4 characters or fewer must be typed in full.
- Implementation: use an **exact-match dictionary** that maps both full name and abbreviation to the handler. Do **not** use prefix matching.

I measured the underlines from the PDF vector graphics (PyMuPDF line extraction):

| Full name | Accepted abbreviation | Full name | Accepted abbreviation |
|---|---|---|---|
| ACCIN | **ACCI** | MOTION | **MOTI** |
| ACCOUT | **ACCO** | POINT | **POIN** |
| AFWRITE | **AFWR** | RESUME | **RESU** |
| ANALYS | **ANAL** | RSOUT | **RSOU** |
| CHECK | **CHEC** | STATUS | **STAT** |
| EQTIT | **EQTI** | THFILE | **THFI** |
| EQUAKE | **EQUA** | THTIT | **THTI** |
| FORCE | **FORC** | WPASS | **WPAS** |
| HOUSE | **HOUS** | WRITE | **WRIT** |
| INCOH | **INCO** | RELFILE | *whole name underlined → no abbreviation* |
| LFREQ | **LFRE** | STRESS | *not underlined in the manual (see Open Q 1)* |

Names of 4 characters or fewer (full name only): `AMP AOPT CORR DAMP DYNP EOUT FREQ INP ME MOPT NOUT RSIN SACC SAVE SFOU SITE SOIL SPRO SRS SSAF SSTR SYMM TPSD TIT TOPL WAVE RELD RDND`.

### 1.4 Token types and parsing conventions *(implementation decisions where the manual is silent)*
| Token kind | Rule |
|---|---|
| integer flag/option | Parse as int. Accept `1.0` → 1 (Fortran-era decks). |
| real | Accept Fortran exponents: `1e+008`, `1.0E9`, `1.0D9`. |
| file path | Rest-of-token. Windows paths with `\` and drive letters. If relative, resolve against the model directory (as INP does). |
| title / label (TIT, THTIT, EQTIT, SSAF `<title>`, DYNP `<label>`) | Take the **remainder of the line** as the title, commas included, for the last-position title args (TIT, THTIT, EQTIT, SSAF). DYNP `<label>` is last too. *(See Open Q 3.)* |
| node/element list (NOUT, EOUT) | Last-position list. Separators: blank, tab, `,`, `;`. `a-b` is an inclusive range, as in the dialog list boxes (pp. 159, 165). Range expansion must reject descending ranges. |
| omitted trailing args | Optional args (`[simul]`, `[tpsd]`, LFREQ, SYMM, WRITE, TOPL tail) take their documented defaults. Omitted **required** args: set to 0 and warn. The PREP free-format convention reads a blank field as 0, and the dialogs say "leave blank (or type 0)". |
| parenthesised token | Only for the **legacy SOIL** command: a trailing `( … )` Fortran-format flag. Commas inside the parentheses must **not** split tokens (e.g. `(5X,F10.4)`). See SOIL. |

### 1.5 Semantic classes of commands (state behaviour)  [CORE]
| Class | Commands | Behaviour |
|---|---|---|
| **Record setter** | ANALYS, AOPT, EQUAKE, FORCE, HOUSE, INCOH, MOPT, MOTION, POINT, SITE, SOIL, STRESS, WPASS, RELD | Replaces the whole option record of that module. Every positional arg is stored. |
| **Indexed setter** | ACCIN, ACCOUT, RSIN, RSOUT, TPSD (index = spectrum no.); CORR (pair no.); DYNP (pair no. within label); ME (motion no.); SYMM (plane no.); WAVE (wave type); SPRO, SACC, SRS, SSTR, SSAF, SFOU (sublayer no.) | Creates or overwrites the entry at that index or key. |
| **List append / delete** | FREQ (per set), DAMP, TOPL, AMP (per motion) | A non-zero first value **appends** non-zero values. A first value of **0 deletes** the list (or set). |
| **Request append** | NOUT, EOUT, RDND | Appends an output request. *(Deletion syntax not documented, Open Q 9.)* |
| **String setter** | TIT, THTIT, EQTIT, THFILE, RELFILE | Sets one string. |
| **Action** | AFWRITE, CHECK, INP, LFREQ, RESUME, SAVE, STATUS, WRITE | Performs an operation. Does not change analysis options. |

---

## 2. Master table of §9.2 commands

| § | Command | Abbr. | Module / scope | Purpose | Class | Tag |
|---|---|---|---|---|---|---|
| 9.2.1 | ACCIN | ACCI | EQUAKE | acceleration input file (seed record / external) for spectrum *no* | indexed | IO |
| 9.2.2 | ACCOUT | ACCO | EQUAKE | simulated acceleration output file for spectrum *no* | indexed | IO |
| 9.2.3 | AFWRITE | AFWR | all | CHECK, then write module input decks | action | IO |
| 9.2.4 | AMP | — | HOUSE | spectral amplification ratios (SAR) for input motion *no* | list | ADV |
| 9.2.5 | ANALYS | ANAL | ANALYS | ANALYS options | record | CORE |
| 9.2.6 | AOPT | — | AFWRITE/CHECK | which modules are written/checked | record | IO |
| 9.2.7 | CHECK | CHEC | all | validate model and options | action | CORE/UI |
| 9.2.8 | CORR | — | EQUAKE | time-varying correlation pair *no* | indexed | ADV |
| 9.2.9 | DAMP | — | SOIL, MOTION (RS) | RS damping-ratio list | list | CORE |
| 9.2.10 | DYNP | — | SOIL | G/Gmax–γ and D–γ curve points of dynamic property *label* | indexed | CORE |
| 9.2.11 | EOUT | — | STRESS | element output request | request | IO |
| 9.2.12 | EQTIT | EQTI | EQUAKE | spectra title | string | UI |
| 9.2.13 | EQUAKE | EQUA | EQUAKE | EQUAKE options | record | CORE |
| 9.2.14 | FORCE | FORC | FORCE | FORCE options | record | CORE |
| 9.2.15 | FREQ | — | SITE/ANALYS… | add/delete frequency numbers of a frequency set | list | CORE |
| 9.2.16 | HOUSE | HOUS | HOUSE | HOUSE options (incl. global gravity, ground elevation) | record | CORE |
| 9.2.17 | INCOH | INCO | HOUSE | incoherency options | record | ADV |
| 9.2.18 | INP | — | UI | read commands from file | action | IO |
| 9.2.19 | LFREQ | LFRE | UI | list frequency sets | action | UI |
| 9.2.20 | ME | — | HOUSE | multiple-excitation motion *no* (node range) | indexed | ADV |
| 9.2.21 | MOPT | — | model | model options | record | CORE |
| 9.2.22 | MOTION | MOTI | MOTION | MOTION options | record | CORE |
| 9.2.23 | NOUT | — | MOTION | nodal output request | request | IO |
| 9.2.24 | POINT | POIN | POINT | POINT options | record | CORE |
| 9.2.25 | RESUME | RESU | UI | reload last SAVE | action | UI |
| 9.2.26 | RSIN | — | EQUAKE | target RS input file for spectrum *no* | indexed | IO |
| 9.2.27 | RSOUT | RSOU | EQUAKE | RS output file for spectrum *no* | indexed | IO |
| 9.2.28 | SACC | — | SOIL | acceleration output for sublayer | indexed | IO |
| 9.2.29 | SAVE | — | UI | save model (binary) | action | UI |
| 9.2.30 | SFOU | — | SOIL | Fourier spectrum output (**not usable in this version**) | indexed | UI |
| 9.2.31 | SITE | — | SITE | SITE options (incl. Δt, NFFT, Δf, frequency set) | record | CORE |
| 9.2.32 | SOIL | — | SOIL | SOIL options | record | CORE |
| 9.2.33 | SPRO | — | SOIL | soil profile: sublayer → layer property + dynamic property | indexed | CORE |
| 9.2.34 | SRS | — | SOIL | response-spectrum output for sublayer | indexed | IO |
| 9.2.35 | SSAF | — | SOIL | spectral amplification factor output | indexed | IO |
| 9.2.36 | SSTR | — | SOIL | stress/strain output for sublayer | indexed | IO |
| 9.2.37 | STATUS | STAT | UI | list global info | action | UI |
| 9.2.38 | STRESS | (none shown) | STRESS | STRESS options | record | CORE |
| 9.2.39 | SYMM | — | model/HOUSE | symmetry / anti-symmetry plane or line | indexed | CORE |
| 9.2.40 | THFILE | THFI | SOIL/MOTION/STRESS/RELDISP | control-motion acceleration file | string | IO |
| 9.2.41 | THTIT | THTI | same | its title | string | UI |
| 9.2.42 | TPSD | — | EQUAKE | target PSD file for spectrum *num* | indexed | ADV |
| 9.2.43 | TIT | — | model | model title | string | UI |
| 9.2.44 | TOPL | — | SITE | free-field top-layer list (L-layer numbers) | list | CORE |
| 9.2.45 | WAVE | — | SITE | wave-type field and participation | indexed | CORE |
| 9.2.46 | WPASS | WPAS | HOUSE | wave-passage data | record | ADV |
| 9.2.47 | WRITE | WRIT | UI | write model to `.pre` | action | IO |
| 9.2.48 | RELD | — | RELDISP | RELDISP options | record | IO |
| 9.2.49 | RELFILE | — | RELDISP | reference-node `.TFI` file | string | IO |
| 9.2.50 | RDND | — | RELDISP | add node + DOF flags to output list | request | IO |

Commands documented elsewhere but tied to this group: `MDL,<Model>,<Path>` (§9.7.20, needed by SAVE/RESUME/AFWRITE when the model came from INP), `MDLNAME`, `GRAVITY,<grav>` and `GROUNDELEV,<elev>` (§9.7.16–17, partial setters of HOUSE values), `L` (§9.4.18, soil layers referenced by TOPL/SITE/SPRO), `INTGEN`/`INT` (interaction nodes), `BINOUT` (§9.18.3), `ANSYSMODELTYPE` (§9.11.2), `AFWRBAT` (§9.7.2), `RUNxxx` (§9.12).

---

## 3. Shared analysis variables  [CORE]

The Analysis dialog tabs share variables: changing one in any tab changes it everywhere (p. 115). Several dialog fields therefore have **no argument in their own module's command**. They are owned by another command. Implement **one** storage location per variable.

| Shared variable | Owning command (arg) | Also shown / used in | Notes |
|---|---|---|---|
| Acceleration of gravity for **SSI** (`g`) | `HOUSE <gravity>` (also `GRAVITY,<grav>`) | SITE, FORCE tabs; EQUAKE units | ft/s² (BS) or m/s² (IS). Error 1 if ≤ 0. Separate from the SOIL gravity. |
| Acceleration of gravity for **SOIL** free-field | `SOIL <grav>` | SOIL tab | Independent of the SSI value (p. 119). |
| Ground elevation | `HOUSE <gelev>` (also `GROUNDELEV`) | — | Decides structure vs excavated soil for SOLID/PLANE with ETYPE 0/2 below grade. |
| Time step of control motion `Δt` (s) | `SITE <delt>` | SOIL, HOUSE (incoh.), FORCE, MOTION, STRESS, RELDISP, EQUAKE ("Time Step") | — |
| Number of Fourier components `NFFT` | `SITE <nft>` | SOIL, HOUSE, FORCE, MOTION, STRESS, RELDISP | Power of 2. Otherwise the nearest power of 2 is written with Warning 9. |
| Frequency step `Δf` (Hz) | `SITE <fstep>` | FORCE | 0/blank → Δf = 1/(Δt·NFFT). |
| Frequency set number | `SITE <freq>` | HOUSE, FORCE, ANALYS (when `<fopt>`=0), STRESS, RELDISP | Set defined by FREQ. |
| Control point layer | `SITE <cl>` | SOIL "Control Point Layer" | Top of that layer (1 = surface). |
| Analysis type (seismic / vibration) | `ANALYS <type>` | MOTION, STRESS "Type of Analysis" | — |
| Coherent/incoherent, wave passage, multiple excitation | `HOUSE <coh>,<wpass>,<me>` | ANALYS tab | — |
| Multiple-excitation motions | `ME`, `AMP` | ANALYS tab | — |
| Time-history scaling: multiplication factor, max value, first/last record, pairs-format flag | `MOTION <mult>,<max>,<rec1>,<rec2>,<fopt>` | STRESS, RELDISP (SOIL shows Mult./Max fields too, Open Q 6) | — |
| Acceleration history file and title | `THFILE`, `THTIT` | SOIL, MOTION, STRESS, RELDISP | — |
| RS damping list | `DAMP` | SOIL (RS group), MOTION (RS group) | — |

---

## 4. Command-by-command specification

### 9.2.1 ACCIN — acceleration input file (EQUAKE)  [IO]
**Syntax:** `ACCIN,<no>,<file>` (also `ACCI`).

| # | Arg | Meaning | Values |
|---|---|---|---|
| 1 | `<no>` | file/spectrum number | 1–3. The EQUAKE dialog Spectrum Number spinner allows at most three spectra, one per translational direction. |
| 2 | `<file>` | full path of the acceleration input file | text |

- **Semantics:** sets the acceleration input history used by EQUAKE for spectrum `<no>`. Dialog "Optional Spectrum Files ▸ Acceleration Input File". It is used when **Accel. Record** ("seed record": the simulated motion keeps the seed's Fourier phasing) or **External Accel.** (RS of an external history with the same Δt and number of points) is selected. Enabled by `EQUAKE <accopt>` = 1.
- **Checks:** Error 88 if the acceleration-input option is on and the file name is blank.
- **Example:** `ACCIN,1,C:\SSI\Demo1\seedNS.acc`

### 9.2.2 ACCOUT — acceleration output file (EQUAKE)  [IO]
**Syntax:** `ACCOUT,<no>,<file>` (also `ACCO`).

| # | Arg | Meaning |
|---|---|---|
| 1 | `<no>` | spectrum number (1–3) |
| 2 | `<file>` | path of the simulated acceleration history (`.acc`). `.vel` and `.dis` take the same base name. |

- **Checks:** Error 87 if blank for a defined spectrum.
- **SASSI-EDU:** a blank ACCOUT of a spectrum that runs is `<model>_eq<no>.acc` in the model folder, reported as
  Warning EDU-29 (requirements §7.19, D-W5-08; `EDUOPT,DEFAULTS,OFF` restores Error 87).
- **Example:** `ACCOUT,1,C:\SSI\Demo1\H1.acc`

### 9.2.3 AFWRITE — write analysis (module input) files  [IO][CORE]
**Syntax:** `AFWRITE` (also `AFWR`). No arguments.

- Writes the input decks for every module enabled by **AOPT** (Options ▸ Analysis ▸ AFWRITE tab) into the **model directory**. The file name is the **model name + the module's postfix/extension** (Modules ▸ Extension window).
- Default extensions (spec 02 / 05d): EQUAKE `.equ` (output `_equake.out`, from the Extension-window screenshot), SOIL `.soi`, SITE `.sit`, POINT `.poi`, HOUSE `.hou`, FORCE `.frc`, ANALYS `.anl`, MOTION `.mot`, STRESS `.str`, RELDISP `.rdi`, NONLINEAR `.eql`.
- **Algorithm:**
  1. Run **CHECK** (§9.2.7) for the enabled modules.
  2. For each enabled module, write its deck **only if CHECK found no *errors* for it**. Warnings do not block. Message: "Any module that does not pass CHECK will not have the associated module input written" (§10.1).
  3. Apply the write-time fix-ups that CHECK warnings describe:
     - **W1/W4:** a gap node or unused node is written with all DOF fixed, in the analysis file only. The model is not changed.
     - **W7:** an empty group is skipped.
     - **W8:** only the first 100 top layers are written (see Open Q 12).
     - **W9:** the nearest power of 2 is written for NFFT.
     - MOTION: duplicate output nodes are deleted (p. 159).
- **Prerequisite:** the model name and path must be known, from the database or `MDL,<Model>,<Path>`. After a bare INP the commands "may run but there is no way to tell where the models will be written" (§5.4). The implementation should raise an error if no model path is set.
- The manual notes that setting up an analysis is complex. The UI dialogs are recommended.
- **Option AA warning [ADV]:** with ANSYS Model Input, AFWRITE writes a `.hou` that holds only topology, with fake materials. It is **not** runnable by standard HOUSE.
- `Options ▸ Write` adds extended behaviour (see §5.3 below). With nothing checked there, AFWRITE writes exactly what old PREP wrote.

### 9.2.4 AMP — spectral amplification ratios for multiple excitation (HOUSE)  [ADV]
**Syntax:** `AMP,<no>,<a1>,<a2>,…,<a100>`

| # | Arg | Meaning | Values |
|---|---|---|---|
| 1 | `<no>` | reference input motion number (same index as ME) | 1–10 (ME range) |
| 2…101 | `<a1>…<a100>` | SAR values, one per SSI frequency of the active frequency set | real. 0 < a ≤ 10 (Error 118 if not in [0,10]). `<a1>` = 0 → delete the list. |

- **Meaning [CORE/ADV]:** at each SSI frequency, the SAR is the ratio of the complex ATF of the "local" motion of an isolated foundation (or zone) to the "reference" motion computed for the single control motion. It introduces non-uniform amplitude input. A complex SAR also carries phasing (`HOUSE <cmplxspec>`=1).
- **Count rule:** the number of ratios must equal the number of frequencies in the selected set (Error 119).
- Up to 100 values per line. Sets can have up to 500 frequencies, so **repeated AMP lines for the same `<no>` append** *(inferred from FREQ/DAMP/TOPL semantics, Open Q 8)*.
- Dialog: HOUSE tab ▸ Spectral Amplification text box ("1,1,1,…"). Separators: blank or comma.
- SAR usually differ per input direction (X, Y, Z), so separate `.hou` files and HOUSE+ANALYS restart runs per direction are needed (p. 142).
- **Example:** `AMP,2,1.00,1.02,1.05,1.10,1.18,1.25,1.20,1.12,1.06,1.03`

### 9.2.5 ANALYS — ANALYS module options  [CORE]
**Syntax:** `ANALYS,<opmode>,<type>,<mode>,<save>,<prnt>,<fopt>,<ang>,<xc>,<yc>,<zc>,<impe>,[simul]` (also `ANAL`).

| # | Arg | Meaning | Values |
|---|---|---|---|
| 1 | `<opmode>` | operation mode | 0 = complete solution, 1 = data check only |
| 2 | `<type>` | analysis type | 0 = seismic (load from FILE1), 1 = foundation vibration / external forces (load from FILE9) |
| 3 | `<mode>` | analysis (restart) mode | 0 = initiation, 1 = new structure, 2 = new seismic environment, 3 = new dynamic loading |
| 4 | `<save>` | restart file save | 0 = do not save, 1 = save restart files |
| 5 | `<prnt>` | transfer-function print option | 0 = print complex TF (Re and Im separately, 6 DOF/node), 1 = amplitude only |
| 6 | `<fopt>` | frequency source | 0 = frequency set (SITE `<freq>`), 1 = all frequencies in FILE1 (seismic) or FILE9 (vibration) |
| 7 | `<ang>` | coordinate transformation angle (deg) | angle from the SITE local x′ axis to the global x axis. Range [0,360) (Error 63). |
| 8–10 | `<xc>,<yc>,<zc>` | foundation reference / control point | coordinates (model length units) |
| 11 | `<impe>` | global foundation impedance calculation | 0 = none, 1 = six diagonal terms, 2 = full 6×6 rigid-body matrix |
| 12 | `[simul]` | number of simultaneous cases | integer, **default 0** |

- **Mode codes vs manual "Mode n" names.** The command uses 0-based codes. The dialog text says "Initiation (Mode 1), New Structure (Mode 2), New Seismic Environment (Mode 3), New Dynamic Loading (Mode 3)". Mapping: `<mode>` 0 ↔ Mode 1, 1 ↔ Mode 2, 2 ↔ Mode 3, 3 ↔ Mode 3 (identical to 2 but with a new FILE9). "New Load Vector (Mode 6)" exists only in batch, is undocumented and **has no command code**. [CORE]
- **Restart files:**
  - `<save>`=1 keeps `COOXxxx` and `COOTKxxx` (xxx = 3-digit frequency order). New Structure needs only the COOX files.
  - Restart also needs `DOFSMAP`, `FILE90` and `FILE91` from HOUSE (plus the index files `COOXI`, `COOTKI`, see spec 02).
- **Frequency survey:** ANALYS checks that every requested frequency is present in the input files and **stops** if one is missing. A set can be split into subsets, run separately and merged with COMBIN.
- **Control point** (`xc,yc,zc`): reference for oblique input motion and/or for the rigid-body global impedance integration. `zc` is used only for the impedance reference.
- **Global impedance (`<impe>`)** [CORE/ADV]:
  - Rigid-body integration of the nodal impedances. **3D only**.
  - For embedded models it is valid only with **FI-FSIN**. These are "unconstrained" impedances: an infinitely flexible foundation, equal to the rigid-foundation values only for surface models.
  - Option 1 writes the text files `FOUNSTIF` (dynamic stiffness), `FOUNDASH` (viscous dashpot coefficients), `FOUNDAMP` (effective damping ratios) and `FOUNIMPD` (|impedance|).
  - For a true two-step analysis the manual recommends Option A instead.
- **Simultaneous cases `[simul]`** [CORE]:
  - Coherent seismic: 1. Solves X, Y and Z together and writes `FILE8X/Y/Z`. Needs `FILE1X/Y/Z` from three SITE runs (SV in x′ at 0°, SH in y′ at 0°, P in z at 0°) and **`<ang>` = 0**.
  - Incoherent: number of stochastic simulations, ≤ 50, equal to the HOUSE value. Writes `FILE8xxx`, up to 150.
  - Vibration: up to 500 load cases (`FILE9xxx` → `FILE8xxx`).
  - The default 0 means a classic single-case run *(inferred, Open Q 10)*.
  - Too many cases for the RAM gives an "access violation".
- **Dialog-only items with no argument here:** Coherent/Incoherent, "Wave Passage Effects Included", Multiple Excitation (all from HOUSE/ME). **Free-Field Load (FFL) vs Free-Field Motion (FFM)** incoherency randomization has no documented argument (Open Q 11). The manual recommends FFL, and FFM only for surface structures. "Delete Restart Files" is greyed out.
- **Checks:** Error 44 (frequency set not defined), Error 63 (angle).
- **Example:** initiation, seismic, save restart, amplitude print, set from SITE, no rotation, no impedance:
  `ANALYS,0,0,0,1,1,0,0.0,0.0,0.0,0.0,0,0`

### 9.2.6 AOPT — AFWRITE/CHECK module selection  [IO]
**Syntax:** `AOPT,<EQUAKE>,<SOIL>,<DEP>,<SITE>,<POINT>,<HOUSE>,<DEP>,<FORCE>,<ANALYS>,<COMBIN>,<MOTION>,<STRESS>,<RELDISP>,<PANEL>` — **14 flags**.

| # | Flag | Module | Dialog check box (AFWRITE tab, p. 183) |
|---|---|---|---|
| 1 | `<EQUAKE>` | EQUAKE | EQUAKE |
| 2 | `<SOIL>` | SOIL | SOIL |
| 3 | `<DEP>` | module not in this version, **always 0** | LIQUEF (greyed) *(position match)* |
| 4 | `<SITE>` | SITE | SITE |
| 5 | `<POINT>` | POINT | POINT |
| 6 | `<HOUSE>` | HOUSE | HOUSE |
| 7 | `<DEP>` | not in this version, **always 0** | PINT (greyed) *(position match)* |
| 8 | `<FORCE>` | FORCE | FORCE |
| 9 | `<ANALYS>` | ANALYS | ANALYS |
| 10 | `<COMBIN>` | COMBIN | COMBIN |
| 11 | `<MOTION>` | MOTION | MOTION |
| 12 | `<STRESS>` | STRESS | STRESS |
| 13 | `<RELDISP>` | RELDISP | RELDISP |
| 14 | `<PANEL>` | manual says "flag for the RELDISP module" (apparent copy error) | NONLINEAR *(position match: panel = nonlinear wall-panel module)* |

- A flag of 0 excludes the module from both AFWRITE and CHECK. Non-zero (1) includes it.
- **Example**, matching the screenshot (FORCE, COMBIN, NONLINEAR off): `AOPT,1,1,0,1,1,1,0,0,1,0,1,1,1,0`

### 9.2.7 CHECK — data check  [CORE][UI]
**Syntax:** `CHECK` (also `CHEC`).

- Checks all data (model and analysis options of the AOPT-enabled modules) for **errors and warnings**. It works on partial models too.
- It **simulates AFWRITE**, so analysis options must be set first.
- Results go to the **Check Errors window** and to the model's **`.err` file**, with headers per module.
- **Options ▸ Check dialog** (p. 114): `Show Warnings` ☑, `Show Errors` ☑, `Suppress Error Window` ☐, `Break Check at [100] Messages` (per message type, per module).
  - Totals are always counted, but messages after the break number are not printed.
  - These are global session settings, not saved, and reset on each UI start. [UI]
- Fatal errors stop the check: Error 42 (no nodes), Error 43 (no groups).
- §6 below lists the check rules that apply to the commands in this file.

### 9.2.8 CORR — spectra correlation pairs (EQUAKE)  [ADV]
**Syntax:** `CORR,<no>,<time>,<val>`

| # | Arg | Meaning |
|---|---|---|
| 1 | `<no>` | pair index (row of the dialog Correlation grid) |
| 2 | `<time>` | time (s) |
| 3 | `<val>` | correlation coefficient at that time, ≤ 1 (Error 94) |

- Defines the time-varying correlation between the horizontal X and Y components. It is active only when `EQUAKE <corr>`=1 (Error 93 if the option is on and no pairs are defined). Indexed set: it overwrites row `<no>`.
- *(inferred)* Between rows the correlation is interpolated linearly in time. The manual does not say.
- **Example:** `CORR,1,0.0,0.0` / `CORR,2,5.0,0.1` / `CORR,3,15.0,0.1`

### 9.2.9 DAMP — damping ratios for response spectra  [CORE]
**Syntax:** `DAMP,<d1>,…,<d10>`

- Adds the **non-zero** values to the RS damping-ratio list. `<d1>` = 0 deletes the whole list. At most 10 values per line.
- Values are **decimal fractions** (0.05 = 5 %), each in (0,1) (Error 72).
- The list is used by **SOIL** response spectra (Error 107 if empty when SOIL RS is requested) and by **MOTION** ISRS.
- MOTION supports at most **5** damping values (§1.2). CHECK should warn above 5 *(inferred)*.
- RS file naming uses the damping **order number** (`xxxxxTR_yzz.RS`, zz = 01, 02…), so the list order matters. [IO]
- **Example:** `DAMP,0` then `DAMP,0.02,0.05`

### 9.2.10 DYNP — dynamic soil property curves (SOIL)  [CORE]
**Syntax:** `DYNP,<no>,<sg>,<g>,<sd>,<d>,<label>`

| # | Arg | Meaning | Units / range |
|---|---|---|---|
| 1 | `<no>` | point (pair) number within the curve set | 1… (max 11 data per curve, §1.2) |
| 2 | `<sg>` | shear strain of the G-curve point | **percent** (dialog axis "Shear Strain %") |
| 3 | `<g>` | normalized shear modulus G/Gmax ("Mod. Red.") | 0–1 (Error 98) |
| 4 | `<sd>` | shear strain of the damping-curve point | percent |
| 5 | `<d>` | damping ratio | **percent** (table values 0.4, 0.8, 1.5, 3, 4.6…) *(inferred from p. 206; see Open Q 5)* |
| 6 | `<label>` | name of the dynamic soil property (e.g. `Sand`, `Clay`, `Rock`) | text, case-sensitive |

- The two curves keep separate strain abscissae (`sg` vs `sd`) but share the index `<no>`.
- The "Select Dynamic Soil Property" editor has columns `Strain | Mod. Red. | Strain | Damp` plus a `Title`, and New/Edit/Delete buttons.
- Labels are referenced by `SPRO <dynprop>`.
- **Checks:**
  - Error 97: label has no G curve.
  - Error 99: label has no D curve.
  - Error 98: G outside [0,1].
  - Error 96: more than **15** dynamic properties used by the sublayers. §1.2 says 100 curves, so there is a conflict (Open Q 13).
- **Algorithm hint [CORE]:** SOIL interpolates G/Gmax(γ) and D(γ) at the effective strain. The usual practice is **linear in log10(γ)** *(inferred, SHAKE practice)*.
- **Example** (values from the manual screenshot, property "Rock"):
  ```
  DYNP,1,0.0001,1.0,0.0001,0.4,Rock
  DYNP,2,0.001,1.0,0.0003,0.8,Rock
  DYNP,3,0.01,0.9875,0.001,1.5,Rock
  DYNP,4,0.1,0.9525,0.003,3.0,Rock
  DYNP,5,1.0,0.90,0.01,4.6,Rock
  ```

### 9.2.11 EOUT — element output request (STRESS)  [IO]
**Syntax:** `EOUT,<code1>,…,<code12>,<group>,<element list>`

| # | Arg | Meaning |
|---|---|---|
| 1–12 | `<code1>…<code12>` | output request per component (component order per element type below). **0 = No Request, 1 = Print Only Maximum, 2 = Print Maximum and Save Time History** *(from the dialog radio group and the 12-character "Output Code" column `000000000000`, p. 162)* |
| 13 | `<group>` | element group number |
| 14… | `<element list>` | element numbers in that group. Ranges `a-b` allowed (Open Q 2). Must be ascending. |

- Component order (code index → component), from §6.5.4 STRESS. Unused trailing codes must be 0.

| Type | Components (code 1 → n) |
|---|---|
| SOLID (7) | σxx, σyy, σzz, τxy, τxz, τyz (stress/strain, global axes, at centroid), octahedral shear stress |
| BEAMS (12) | F1, F2, F3, M1, M2, M3 at node I; F1, F2, F3, M1, M2, M3 at node J (local beam axes) |
| SHELL (6) | Sx′x′, Sy′y′, Sx′y′ (F/L/L); Mx′x′, My′y′, Mx′y′ (FL/L), local axes |
| TSHELL (8) | Nxx, Nyy, Nxy, Qxz, Qyz (F/L); Mxx, Myy, Mxy (FL/L) |
| PLANE (3) | σxx, σzz, τxz (global axes, centre) |
| SPRING (6) | Fx, Fy, Fz, Mxx, Myy, Mzz |

- STRESS may also print components needed internally (e.g. all six stresses for octahedral). Their time histories are saved only if requested.
- The nonlinear-soil SSI option requires that the nonlinear soil group (all its SOLID/PLANE elements) be requested.
- **Checks:**
  - Error 79: no element request.
  - Error 80: illegal group.
  - Error 81: illegal element.
  - Error 82: element requested more than once.
- **Example:** SOLID group 3, elements 1–20, max of the six stresses plus a saved octahedral history:
  `EOUT,1,1,1,1,1,1,2,0,0,0,0,0,3,1-20`

### 9.2.12 EQTIT — EQUAKE spectra title  [UI]
**Syntax:** `EQTIT,<title>` (also `EQTI`). Sets the title written to the `.equ` deck (dialog "Spectra Title").
**Example:** `EQTIT,Design spectra 0.3g, 5% damping` (the rest of the line is the title).

### 9.2.13 EQUAKE — EQUAKE options  [CORE]
**Syntax:** `EQUAKE,<accopt>,<nrfreq>,<rand>,<damp>,<dur>,<corr>,<seeds>,[tpsd]` (also `EQUA`).

| # | Arg | Meaning | Values / notes | Dialog field (screenshot value) |
|---|---|---|---|---|
| 1 | `<accopt>` | acceleration input files option | 0 disabled, 1 enabled (seed records / ACCIN files; see Open Q 4 on External Accel.) | Accel. Record (☐) |
| 2 | `<nrfreq>` | number of frequencies in the spectrum files | > 0 (Error 91). **Must equal the number of records of each RSIN file** (Error 89). | Number of Frequencies (24) |
| 3 | `<rand>` | initial random number (seed) | > 0 (Error 90). "5 digits" recommended. | Initial Random SEED (11975) |
| 4 | `<damp>` | damping of the target design spectrum | ratio, e.g. 0.05 | Damping Value (0.05) |
| 5 | `<dur>` | total duration (s) | > 0 (Error 92). SRP 3.7.1 wants ≥ 20 s. | Total Duration (15) |
| 6 | `<corr>` | correlated X–Y spectra option | 0/1. Needs CORR pairs. | Correlated (☐) |
| 7 | `<seeds>` | number of random seeds (trials) per acceleration file | the best-fitting history is kept. Should be a non-zero integer (UI text). The screenshot shows 0. | Number of SEEDs (0) |
| 8 | `[tpsd]` | use target PSD files | 0 (default) / 1. Needs TPSD. | Use Target PSD (☑) |

- **Time step:** not in this command. It is the shared Δt (`SITE <delt>`). The text says it "should be the same" as for the other modules.
- **Units warning:** EQUAKE derives velocity, displacement and PSD units from the gravity value set in SOIL/SITE/HOUSE. Generated acceleration is in g.
- **Checks:** Error 84 (no RS input file), Error 85 (RS input file missing), Error 86 (RS output blank), Error 87 (acceleration output blank), Errors 88–94.
- **Example:** `EQUAKE,0,24,11975,0.05,20.0,0,5,0`

### 9.2.14 FORCE — FORCE options  [CORE]
**Syntax:** `FORCE,<opmode>` (also `FORC`). `<opmode>`: 0 = complete solution, 1 = data check only.

- All other FORCE dialog fields are **shared**: gravity (HOUSE), Δf, Δt, NFFT, frequency set (SITE).
- The loads come from the F/MM/FREAD/MREAD commands (§9.5). Each load is a reference time history scaled by a load factor and delayed by an arrival time t0 (phase lag).
- **Checks:** Error 61 (no forces), Error 62 (force on an undefined node). Warnings 10/11 (force/moment on a fixed DOF is ignored).
- **Example:** `FORCE,0`

### 9.2.15 FREQ — frequency set editing  [CORE]
**Syntax:** `FREQ,<ndx>,<f1>,…,<f10>`

| # | Arg | Meaning |
|---|---|---|
| 1 | `<ndx>` | frequency set number |
| 2–11 | `<f1>…<f10>` | **frequency numbers** (positive integers) added to the set. `<f1>` = 0 deletes set `<ndx>`. At most 10 per line; repeat the line to add more. |

- **Frequency computation [CORE]:** SITE sorts the set ascending and **stops on duplicates**. Then f_k = n_k·Δf with:
  ```
  Δf = SITE<fstep>                if fstep > 0   (harmonic analysis: fstep required)
  Δf = 1 / (SITE<delt> · SITE<nft>)  otherwise    (time-history analysis: delt & nft required)
  ```
  Example: Δt = 0.005 s, NFFT = 4096 → Δf = 0.048828 Hz. n = 512 → 25.0 Hz.
- Limits: at most **500** SSI analysis frequencies (SITE/POINT/HOUSE/FORCE/ANALYS, §1.2). *(inferred)* Each frequency number should satisfy n ≤ NFFT/2 (Nyquist).
- **Checks:** Error 44 (set used but not defined), Error 120 (set empty, RELDISP).
- **Example:**
  ```
  FREQ,1,0
  FREQ,1,2,4,6,8,10,12,14,16,18,20
  FREQ,1,24,28,32,36,40,48,56,64,80,96
  FREQ,1,112,128,160,192,224,256,320,384,448,512
  ```

### 9.2.16 HOUSE — HOUSE options  [CORE]
**Syntax:** `HOUSE,<gravity>,<gelev>,<opmode>,<dim>,<imp>,<coh>,<wpass>,<me>,<cmplxspec>` (also `HOUS`).

| # | Arg | Meaning | Values | Dialog (screenshot) |
|---|---|---|---|---|
| 1 | `<gravity>` | acceleration of gravity, the **model/SSI value** | > 0 (Error 1). 32.2 ft/s² or 9.81 m/s². | Acceleration of Gravity (32.2) |
| 2 | `<gelev>` | ground surface elevation | length | Ground Elevation (−10) |
| 3 | `<opmode>` | operation mode | 0 complete, 1 data check | Solution |
| 4 | `<dim>` | analysis dimension | 0 = 1D (*not available*), 1 = 2D (use POINT2), 2 = 3D (use POINT3) | 3D |
| 5 | `<imp>` | impedance method | 0 = Flexible Volume (FV), 1 = Fast Flexible Volume (FFV), 2 = Flexible Interface (FI) | FI |
| 6 | `<coh>` | soil motion | 0 coherent, 1 incoherent | Incoherent |
| 7 | `<wpass>` | wave passage | 0/1 | ☑ |
| 8 | `<me>` | multiple excitation | 0/1 | ☐ |
| 9 | `<cmplxspec>` | complex spectral amplification ratios | 0/1 | ☐ |

- **Ground elevation [CORE]:** SOLID/PLANE elements below grade are treated as excavated soil when they have ETYPE 0 (default) or ETYPE 2. With ETYPE 2 the MSET material is read as the free-field layer. ETYPE 1 marks structure, possibly including near-field backfill.
- **Method:** for surface foundations FV, FFV and FI are identical. With FFV or FI, sensitivity studies against FV are required (ASCE 4-16, SRP 3.7.2). The interaction nodes come from INTGEN/INT.
- **Incoherency** (`<coh>`=1) [ADV]:
  - 3D full models only (no symmetry). Interaction nodes are numbered bottom-up.
  - Abrahamson and user-defined coherency models need `<wpass>`=1.
  - `<me>`=1 also needs `<wpass>`=1.
  - Parameters come from INCOH, WPASS, ME and AMP.
- **Dialog-only items:** Optimize Model (node renumbering → `.hounew`, `.map`), Non-Linear SSI Input Data (`.pin`), Superposition Mode (Linear/Quadratic), ANSYS Model Input/Type. **None has an argument here** (Open Q 11).
- Partial setters: `GRAVITY,<grav>` and `GROUNDELEV,<elev>` change args 1 and 2 only (§9.7.16–17).
- **Example:** BS units, ground at z = 0, 3D FV, coherent: `HOUSE,32.2,0.0,0,2,0,0,0,0,0`

### 9.2.17 INCOH — incoherency options (HOUSE)  [ADV]
**Syntax:** `INCOH,<gammax>,<gammay>,<gammaz>,<alpha>,<ngp>,<ipr>,<nmodes>,<met>,<HSeed>,<VSeed>,<RandPhz>` (also `INCO`).

The syntax line spells the 3rd argument `<gammmaZ>`; the argument list spells it `<gammaz>`.

| # | Arg | Meaning | Values / guidance | Dialog (screenshot) |
|---|---|---|---|---|
| 1 | `<gammax>` | coherence parameter, X | Luco-Wong. ≥ 0.1 (Error 58). Typical 0.10–0.30. | Coherence Parameter X Dir (0.1) |
| 2 | `<gammay>` | coherence parameter, Y | as above | (0.1) |
| 3 | `<gammaz>` | coherence parameter, Z (vertical) | as above | (0.2) |
| 4 | `<alpha>` | directionality factor α (Abrahamson / user models), **or** mean Vs (Luco-Wong) | α ∈ [0,1]: 0.5 isotropic, 0.1/0.9 strongly directional. For Luco-Wong, Vs > 0 (Error 59). | Alpha Directionality Factor (1000) |
| 5 | `<ngp>` | "number of mesh points per each embedment level" | > 0 (Error 60). The dialog label is "Number of Embedded Layers" (screenshot 0, Open Q 14). | Number of Embedded Layers |
| 6 | `<ipr>` | print incoherent mode contributions per frequency | 0/1 | Print Coherency Matrix |
| 7 | `<nmodes>` | number of incoherent modes | 0 = all N; k > 0 = modes 1…k (k < N); −k = only mode k | Number of Incoh. Modes (0) |
| 8 | `<met>` | unit flag | 0 = British (ft), 1 = SI (m) | (no field visible) |
| 9 | `<HSeed>` | horizontal seed | integer. 0 → deterministic | Motion Incoherency Simulation |
| 10 | `<VSeed>` | vertical seed | integer | 〃 |
| 11 | `<RandPhz>` | random phase angle (deg) | 180 for stochastic, 0 for deterministic | 〃 |

- **Distance metric [CORE/ADV]:** for directional models the effective horizontal distance between interaction nodes i, j is

  `D = [ 2 ( α·DX² + (1−α)·DY² ) ]^(1/2)`

  where DX and DY are the separations along the axes of Line D (rotated by WPASS `<ang>`). α = 0.5 gives D = √(DX² + DY²).
- **Spectral factorization check:** for each frequency the coherency matrix Σ (N×N, unit diagonal) is decomposed into eigenpairs λ_j.
  - With all modes, Σ_{j=1}^{N} λ_j = N (6.1). With m < N modes the sum is < N (6.2).
  - Percent contribution υ_j = 100·λ_j/N (6.3). Cumulative contribution Σ_{j=1}^{m} υ_j (6.4).
  - Recommended criterion: 90 % cumulative. Check (6.1) to detect ill-conditioning. The HOUSE output contains the "I N C O" string.
- **Deterministic vs stochastic:**
  - `HSeed` = `VSeed` = `RandPhz` = 0 → deterministic (median, zero phase between spatial modes).
  - Non-zero seeds and RandPhz = 180 → stochastic, with random phases in [−180°, 180°].
  - Up to **50** simultaneous simulations in one HOUSE run (FILE77001…). The count argument is not in this command (Open Q 11).
- **Recommendations:**
  - All modes (`nmodes` = 0) for stochastic simulation and for linear AS.
  - `nmodes` = −k (one SSI run per mode) for SRSS quadratic.
  - Deterministic approaches are suitable for rigid foundations only. Use stochastic simulation for flexible foundations.
- **Example:** Abrahamson-type isotropic, deterministic, all modes, metres:
  `INCOH,0.1,0.1,0.2,0.5,1,1,0,1,0,0,0`

### 9.2.18 INP — switch input to a file  [IO]
**Syntax:** `INP,<filename>`

- Reads and executes commands from `<filename>`, a `.pre` written by the user or by WRITE. **At EOF, input returns to the keyboard** (command line).
- Relative paths are taken from the **active model's path**. The menu equivalent is Model ▸ Input, which shows a progress bar (current line / total lines).
- Only comments, warnings, errors and model information are echoed, which keeps loading fast.
- *(inferred)* A nested INP inside a `.pre` pushes a new input source and returns to the caller at EOF. Guard against recursion depth.
- After INP alone the model has no name or path; use MDL before SAVE, RESUME or AFWRITE.
- **Example:** `INP,Demo5.pre`

### 9.2.19 LFREQ — list frequency sets  [UI]
**Syntax:** `LFREQ,[start],[end],[step]` (also `LFRE`). Lists the frequency sets from `start` (default 1) to `end` (default the last set) in steps of `step` (default 1).
- *(suggested output)* Set number, count, sorted frequency numbers and their Hz values using the current Δf.

### 9.2.20 ME — multiple excitation motion definition (HOUSE)  [ADV]
**Syntax:** `ME,<no>,<nfirst>,<nlast>,<xc>,<yc>,<zc>`

| # | Arg | Meaning |
|---|---|---|
| 1 | `<no>` | input motion number, 1–10 |
| 2 | `<nfirst>` | first interaction (foundation) node of this motion |
| 3 | `<nlast>` | last foundation node of this motion |
| 4–6 | `<xc>,<yc>,<zc>` | control point of this motion. **Not used in this version.** |

- The nodes of each foundation or zone must form one **continuous, unit-increment node range**. Bottom-up interaction-node numbering still applies. Up to 5000 zones (§1.2, p. 141). The command range is 1–10 (Open Q 15).
- Active only with `HOUSE <me>`=1, which needs `<wpass>`=1.
- **Checks:**
  - Error 115: ME on but no data.
  - Error 116/117: illegal first/last node.
  - Errors 118/119: AMP values.
- **Example:** `ME,1,1,69,0,0,0` / `ME,2,70,138,0,0,0`

### 9.2.21 MOPT — model options  [CORE]
**Syntax:** `MOPT,<incomp>,<matrix>,<mass>,<force>`. The syntax spells the 1st argument `<incomp>`; the list spells it `<icomp>`. Dialog: **Options ▸ Model**.

| # | Arg | Meaning | 0 | 1 | Dialog (screenshot) |
|---|---|---|---|---|---|
| 1 | `<incomp>` | incompatible modes for SOLID elements | include | suppress | Suppress Incompatible Modes ◉ |
| 2 | `<matrix>` | units of GENERAL element matrices (MXM) | mass units | weight units | Mass Matrix ◉ |
| 3 | `<mass>` | repeated nodal mass definitions | add | set (overwrite) | Overwrite Mass ☑ |
| 4 | `<force>` | repeated nodal force definitions | add | set (overwrite) | Overwrite Force ☑ |

- Incompatible modes never apply to excavated-soil SOLIDs.
- With `<matrix>`=1, divide the MXM matrix entries by g (the HOUSE gravity) to get mass *(inferred)*.
- `<mass>`/`<force>` define what happens when MT/MR/F/MM hit a node that already has a value.
- **Example:** `MOPT,1,0,1,1`

### 9.2.22 MOTION — MOTION options  [CORE]
**Syntax (19 args):**
`MOTION,<opmode>,<out>,<step>,<dur>,<res>,<freq1>,<freq2>,<fstep>,<mult>,<max>,<rec1>,<rec2>,<fopt>,<bl>,<smo>,<cplx>,<cnvrt>,<pzadj>,<interp>` (also `MOTI`).

| # | Arg | Meaning | Values / rules | Dialog field |
|---|---|---|---|---|
| 1 | `<opmode>` | operation mode | 0 complete, 1 data check | Operation Mode |
| 2 | `<out>` | output option | 0 full output, 1 transfer functions only | Output Only Transfer Functions |
| 3 | `<step>` | output time-history print step | 0 = print only table; > 1 = print every step-th point. ≥ 0 (Error 67). | (no visible field; Open Q 16) |
| 4 | `<dur>` | total duration of histories to extract/plot (s) | ≥ 0 (Error 68). The output is 20 % longer, to include free vibration. | Total Duration to be Plotted |
| 5 | `<res>` | "not used for Data Check in this version" | write 0 | — |
| 6 | `<freq1>` | first RS frequency (Hz) | ≥ 0 (Error 69) | First Frequency |
| 7 | `<freq2>` | last RS frequency (Hz) | ≥ 0 (Error 70) | Last Frequency |
| 8 | `<fstep>` | total number of RS frequency steps | ≥ 0 (Error 71). **Recommended 0.1–100 Hz with ≥ 301 steps** (SRP 3.7.1). Fewer steps hurt spectrum broadening. | Total Number of Freq. Steps |
| 9 | `<mult>` | time-history multiplication factor | exactly one of mult/max must be non-zero (Errors 77, 78) | Multiplication Factor |
| 10 | `<max>` | target peak of the time history (history rescaled to this max) | see above | Max Value for Time History |
| 11 | `<rec1>` | first record written to the analysis file | ≥ 0, ≤ number of records (Error 74) | First Record |
| 12 | `<rec2>` | last record (**default: last record in file**) | ≥ 0 (Error 75). rec1 ≤ rec2 (Error 76). | Last Record |
| 13 | `<fopt>` | history file format | 0 = 1st line holds the time step, then one acceleration per line; 1 = pairs (time step, acceleration) per line | File Contains Pairs Time Step – Accel. |
| 14 | `<bl>` | baseline correction | 0 = time domain, 1 = frequency domain (per §9.2.22; the dialog shows "No Correction / With Correction", Open Q 17) | Baseline Correction |
| 15 | `<smo>` | smoothing parameter | 1–1000 per §9.2.22. Typical 10–1000. **0 for coherent analysis; must be 0 with `interp` = 6.** Only affects interpolation options 0–5. | Smoothing Parameter |
| 16 | `<cplx>` | TFU/TFI content | 0 amplitude only, 1 complex (amplitude + phase in **radians**). **Needs 1 if RELDISP will be run.** | Save Complex Transfer Functions |
| 17 | `<cnvrt>` | RS of external acceleration files | 0 no, 1 compute RS for the files listed in `CONTTRS.txt` (`.ACC` format, same Δt/length) → `.RSO` | Convert TH to RS ▸ Select External Files |
| 18 | `<pzadj>` | phase adjustment (incoherent) | 0 none; 1 phase adjustment ("minimum delay", near-zero differential phase, EPRI 2007 practice) | Phase Adjustment |
| 19 | `<interp>` | complex TF interpolation scheme | 0 SASSI2000 dense overlapping windows (weighted averaging); 1 original SASSI 1982 non-overlapping windows; 2 dense overlapping windows (averaging); 3 only three overlapping windows (averaging); 4 non-overlapping, one-position shift; 5 non-overlapping, two-position shift; 6 complex bicubic spline (no windowing; needs a dense frequency grid, > 200 frequencies for incoherent) | Interpolation Option |

- **Shared inputs:**
  - NFFT and Δt come from SITE.
  - Damping list: DAMP. Title: THTIT. File: THFILE.
  - Type of analysis: `ANALYS <type>`. For vibration problems the TF is displacement and the requested response depends on Output Control.
- **Scaling algorithm [CORE]:** read records rec1…rec2 (after the format header), then:
  ```
  if mult ≠ 0:  a(t) ← mult · a(t)
  else:         a(t) ← a(t) · max / max|a(t)|
  ```
  Acceleration is in g.
- **Dialog-only items not carried by MOTION:**
  - Save FILE12 or FILE13 (0/1/2).
  - Incoherent SSI ▸ Input (writes `SRSSTF.txt`: `[#modes] [phase option 0|1]`, then one FILE8 name per line, with the coherent FILE8 first when the phase option is 1).
  - Post-processing Save/Restart options, Save Binary Database (BINOUT).
  - These are carried by other commands or files (see 05c).
- **Warnings (manual):**
  - Spline (6) is recommended for incoherent analysis.
  - With phase adjustment or SRSS, do not use the multi-location motions for multiple-time-history analysis of secondary systems.
  - Use RELDISP rather than baseline correction for relative displacements.
- **Example:** full output, 25 s, RS 0.1–100 Hz / 301 steps, scale ×1, records 1–4096, single-column file, no baseline, coherent, complex TF for RELDISP, interpolation 0:
  `MOTION,0,0,0,25.0,0,0.1,100.0,301,1.0,0.0,1,4096,0,0,0,1,0,0,0`

### 9.2.23 NOUT — nodal output request (MOTION)  [IO]
**Syntax:** `NOUT,<dir>,<code1>,…,<code6>,<node list>`

| # | Arg | Meaning |
|---|---|---|
| 1 | `<dir>` | DOF direction: 1 = x, 2 = y, 3 = z, 4 = xx, 5 = yy, 6 = zz |
| 2–7 | `<code1>…<code6>` | 0/1 per output option, in dialog order *(inferred)*: (1) Printed Plot of Transfer Functions, (2) Save Time History of Requested Response, (3) Plot Time History of Requested Response, (4) Plot Acceleration and Velocity R.S. (deprecated, "does not work in newer versions"), (5) Save Acceleration and Velocity R.S., (6) Print Maximum Requested Response |
| 8… | `<node list>` | node numbers, ranges allowed (`1, 3-6 10` in the dialog) |

- Seismic: TF = total acceleration TF, and the requested response is acceleration. Vibration: TF = total displacement TF. RS do not depend on Output Control, so displacement RS cannot be requested.
- Requests on constrained DOFs are ignored.
- **With the HOUSE node optimizer, use the new node numbers** (from `.hounew` / `.map`).
- **Checks:**
  - Error 64: no request.
  - Error 65: illegal node.
  - Error 66: node requested more than once.
  - Duplicates are deleted at AFWRITE (p. 159). The manual thus says both "error" and "auto-fix" (Open Q 18).
- One NOUT per direction per list. WRITE emits up to six NOUT lines per dialog list.
- **Example:** `NOUT,1,1,1,0,0,1,1,1-30`

### 9.2.24 POINT — POINT options  [CORE]
**Syntax:** `POINT,<opmode>,<layer>,<rad>` (also `POIN`).

| # | Arg | Meaning | Values / guidance |
|---|---|---|---|
| 1 | `<opmode>` | operation mode | 0 complete, 1 data check |
| 2 | `<layer>` | last layer number of the near-field zone = number of embedment layers | ≥ 0 (Error 56). **0 for surface structures.** It must reach the bottom of the excavated region. A smaller value means a smaller FILE3 and less work in ANALYS. Max 50 embedment layers. |
| 3 | `<rad>` | radius of the central zone of the point-load (2D or axisymmetric) solution | > 0 (Error 57). It should be of the order of the excavation element size h. |

- **Recommended central-zone radius [CORE]**, from the figures on p. 128:
  - 3D uniform rectangular mesh: r = 0.90 h.
  - 3D circular/triangular mesh: r = 0.85 h.
  - 2D: r = h.
  - Non-uniform meshes: use the average and run sensitivity studies with min/avg/max values.
  - The `RADIUS,<Scale>,<FileName>` command reports element radii.
- Dialog labels: "Number of Embedment Soil Layers" (0) and "Point Load Central Zone Radius" (13.8).
- HOUSE `<dim>` selects POINT2 (2D) or POINT3 (3D).
- **Example:** 3 embedment layers, 3 ft square elements: `POINT,0,3,2.7`

### 9.2.25 RESUME — reload the saved model  [UI]
**Syntax:** `RESUME` (also `RESU`).

- Reloads the binary model data written by the last **SAVE**. This undoes the commands run since that save.
- Needs the model name/path (database or MDL).
- On model open the data loads automatically.
- PREP binary data is **not** compatible with the UI. Move models between them as `.pre` files (WRITE).

### 9.2.26 RSIN — target response-spectrum input file (EQUAKE)  [IO]
**Syntax:** `RSIN,<no>,<file>`. `<no>` = spectrum number (1–3). `<file>` = target RS file: 2 columns (frequency, spectral acceleration) at the damping given by `EQUAKE <damp>`; record count = `<nrfreq>`.
- **Checks:** Error 84 (all three blank), Error 85 (file missing), Error 89 (record count ≠ nrfreq).
- **SASSI-EDU:** with all three blank, RSIN 1 is the built-in RG 1.60 H spectrum `@rg160h_030g.rsi` (0.30 g, 5 %,
  27 rows), reported as Warning EDU-29 (requirements §7.19, D-W5-07); `@<file>` names a built-in file (LIBRARY).
- **Example:** `RSIN,1,C:\SSI\Demo1\H1.rsi`

### 9.2.27 RSOUT — response-spectrum output file (EQUAKE)  [IO]
**Syntax:** `RSOUT,<no>,<file>` (also `RSOU`). Output RS of the generated motion (`.rso`). Error 86 if blank
(SASSI-EDU: `<model>_eq<no>.rso`, Warning EDU-29, requirements §7.19, D-W5-08).
**Example:** `RSOUT,1,C:\SSI\Demo1\H1.rso`

### 9.2.28 SACC — acceleration output for a SOIL sublayer  [IO]
**Syntax:** `SACC,<layer>,<opt>,<outcrop>`

| # | Arg | Meaning | Values |
|---|---|---|---|
| 1 | `<layer>` | sublayer number. Acceleration is at the **top** of the layer. | 1 = surface |
| 2 | `<opt>` | output | 0 none, 1 compute maximum, 2 maximum + save time history (`ACCxxx.TH`) |
| 3 | `<outcrop>` | outcrop motion | 0 within-column, 1 outcrop |

- **Outcrop (impl.) [CORE]:** an outcropping motion is twice the upgoing wave amplitude, 2A, at that interface. A within-column motion is A + B.
- **Example:** `SACC,1,2,0`

### 9.2.29 SAVE — save the model  [UI]
**Syntax:** `SAVE`. Saves the active model in a binary format for RESUME or a later session. Needs the model name/path. PREP and UI binaries are incompatible; exchange models as `.pre`.
- *(impl.)* The format is free to choose (e.g. versioned JSON/HDF5 under the model directory). It must round-trip exactly with RESUME.

### 9.2.30 SFOU — Fourier spectrum output (SOIL) — **not usable in this version**  [UI]
**Syntax:** `SFOU,<layer>,<out>,<save>,<outcrop>,<smooth>,<nrval>`

| # | Arg | Meaning |
|---|---|---|
| 1 | `<layer>` | sublayer number |
| 2 | `<out>` | 0 no computation, 1 computation |
| 3 | `<save>` | 0/1 save |
| 4 | `<outcrop>` | 0/1 outcrop |
| 5 | `<smooth>` | number of smoothing passes, ≥ 0 (Error 111) |
| 6 | `<nrval>` | number of values saved, ≥ 0 (Error 112) |

- Parse and store it (for round-trip WRITE), but **do not act on it**. Grey out the UI. The manual says to compute Fourier spectra with EQUAKE instead.

### 9.2.31 SITE — SITE options  [CORE]
**Syntax (14 args):** `SITE,<opmode>,<mode1>,<fstep>,<nl>,<hs>,<mode2>,<wopt>,<freq1>,<freq2>,<cl>,<cm>,<delt>,<nft>,<freq>`

| # | Arg | Meaning | Values / rules | Dialog (screenshot) |
|---|---|---|---|---|
| 1 | `<opmode>` | operation mode | 0 complete, 1 data check (the dialog shows "Linear Soil / Non-Linear Soil", Open Q 19) | Operation Mode (Non-Linear Soil) |
| 2 | `<mode1>` | Mode 1: transmitting-boundary eigenproblem → FILE2 | 0 skip, 1 write | Mode 1 ☑ |
| 3 | `<fstep>` | frequency step Δf (Hz) | ≥ 0 (Error 48). 0 → computed from Δt and NFFT. | Frequency Step (0) |
| 4 | `<nl>` | number of generated sublayers for half-space simulation | 0 = rigid base, else **4–20** (Error 47). The text recommends 10–20, max 20. | Number of Generated Layers (20) |
| 5 | `<hs>` | half-space layer: L-command soil layer number with the base-rock properties | defined L layer | Halfspace Layer (2) |
| 6 | `<mode2>` | Mode 2: site response → FILE1 | 0 skip, 1 write | Mode 2 ☑ |
| 7 | `<wopt>` | wave combination | 0 = R-, SV-, P-waves (in-plane, 2 DOF/node); 1 = SH- and L-waves (anti-plane, 1 DOF/node) | R-,SV-,and P-Waves |
| 8 | `<freq1>` | "frequency 1", lowest for the wave-ratio curve | recommended 1. > 0 (Error 53). | Frequency 1 (1) |
| 9 | `<freq2>` | "frequency 2", highest for the wave-ratio curve | recommended ≥ NFFT/2 | Frequency 2 (4000) |
| 10 | `<cl>` | layer number of the control point (top of layer) | 1 = ground surface | Control Point Layer (1) |
| 11 | `<cm>` | control motion direction in x′y′z′ | 0 = X(x′), 1 = Y(y′), 2 = Z(z′) | Direction X |
| 12 | `<delt>` | time step of the seismic motion (s) | ≥ 0 (Error 49) | Time Step Control Motion (0.005) |
| 13 | `<nft>` | number of Fourier components | power of 2 (Warning 9) | Nr. of Fourier Component (4096) |
| 14 | `<freq>` | frequency set number (FREQ) | defined (Error 44) | Frequency Set Number (1) |

- **Rules:**
  - Error 45 if `<mode1>` and `<mode2>` are both 0.
  - Error 46 if there are no top layers (TOPL).
  - Error 51 if all wave fields are off.
- **Time history vs harmonic:** time-history analysis needs Δt and NFFT (Δf may be blank). Single-harmonic analysis needs Δf (Δt and NFFT may be blank).
- **Coordinates:** SITE works in x′y′z′, with z′ vertical up and x′ in the vertical plane of propagation.
  - P, SV and R motions lie in x′z′. SH and L motions lie along y′.
  - ANALYS `<ang>` rotates x′ to the global x.
- **Half-space simulation [CORE]:** SITE generates `nl` sublayers whose thickness varies with frequency, with viscous dashpots at the base (body-wave radiation), below the fixed top layers.
- **Wave-ratio curve:** each wave type's participation ratio is given at frequency 1 and 2 (WAVE `<ratio1>`, `<ratio2>`) and interpolated linearly between them. The two frequencies must span the analysis range.
- `<freq1>`/`<freq2>` are probably **frequency numbers** (multiples of Δf), not Hz, given the recommended values "1" and "≥ NFFT/2" and the screenshot 1/4000 (Open Q 20).
- **Simultaneous-case FILE1X/Y/Z:** run SITE three times with (`cm`=0, SV, 0°), (`cm`=1, SH with `wopt`=1, 0°) and (`cm`=2, P, 0°). Copy FILE1 after each run.
- **Example:** vertically incident SV, X control at the surface, 20 half-space sublayers:
  ```
  SITE,0,1,0,20,5,1,0,1,2048,1,0,0.005,4096,1
  WAVE,2,1,1.0,1.0,0.0
  ```

### 9.2.32 SOIL — SOIL options  [CORE]
**Syntax (9 args):** `SOIL,<nrval>,<grav>,<header>,<outcrop>,<save>,<iter>,<ratio>,<gravmult>,<cof>`

| # | Arg | Meaning | Values / guidance | Dialog (screenshot) |
|---|---|---|---|---|
| 1 | `<nrval>` | number of acceleration values read from the history file | > 0 (Error 100). ≤ 32 768 (§1.2). | Number of Values (3000) |
| 2 | `<grav>` | gravity for SOIL free-field analysis | ft/s² (BS, weights in lb/ft³ → see units note) or m/s² (IS, kN/m³) | Gravity Accel. (32.2) |
| 3 | `<header>` | number of header lines in the history file | ≥ 0 (Error 103) | Number of Header Lines (0) |
| 4 | `<outcrop>` | input motion is outcrop (1) or within-column (0) | 0/1 | Assign as Outcrop Motion ☑ |
| 5 | `<save>` | save strain-compatible properties to **FILE88** (read by SITE in nonlinear mode) | 0 skip, 1 save | Save Strain-Compatible Soil Properties ☑ |
| 6 | `<iter>` | number of equivalent-linear iterations | ≥ 0 (Error 105). **Recommended 8.** 0 for vertical input. | Number of Iterations (8) |
| 7 | `<ratio>` | effective (equivalent uniform) / maximum strain | in (0,1) (Error 106). **Recommended 0.60–0.70.** | Equiv. Uniform / Max Strain (0.6) |
| 8 | `<gravmult>` | multiplier for the acceleration of gravity (RS group) | > 0 (Error 108) | Multiplier for Acceleration of Gravity (1) |
| 9 | `<cof>` | cut-off frequency | **disabled**: always set to the Nyquist frequency 1/(2Δt). ≥ 0 (Error 101). | — |

- **Legacy form [IO]:**
  - The argument list changed from old PREP. The old SOIL command has a **format flag as its last argument, in parentheses** (a Fortran read format such as `(8F10.6)`).
  - Parser: if the last token is parenthesised, mark the command as legacy, keep the format string and warn the user to review the arguments.
  - Error 102 (Illegal Reading Format) applies to it.
  - The old argument order is not given (Open Q 21).
- **Units note (p. 119):** BS = ft, ksf/kcf (the dialog text also says lb/ft³ for weight density). IS = m, kN/m², kN/m³. "Tons with meters" is not allowed.
- **Dialog fields with no SOIL argument:**
  - Input Direction (0 horizontal using Vs; 1 vertical using Vp, no iterations).
  - Control Point Layer (= SITE `<cl>`), NFFT and Δt (= SITE).
  - File (= THFILE), Multiplication Factor / Max Value (Open Q 6).
  - Nonlinear Time Domain block [ADV].
- **Equivalent-linear iteration [CORE]** (SHAKE methodology, for implementers):
  1. For each sublayer j, compute γ_eff,j = `ratio`·γ_max,j.
  2. Update G_j = Gmax,j·(G/Gmax)(γ_eff,j) and D_j = D(γ_eff,j) from the DYNP curves of the SPRO label.
  3. Repeat `iter` times (or until converged).
- **Example:** `SOIL,3000,32.2,0,1,1,8,0.65,1.0,0`

### 9.2.33 SPRO — SOIL profile assignment  [CORE]
**Syntax:** `SPRO,<layer>,<prop>,<dynprop>`

| # | Arg | Meaning |
|---|---|---|
| 1 | `<layer>` | SOIL sublayer number (1 = top) |
| 2 | `<prop>` | soil layer property number (an `L` layer: thickness, unit weight, Vp, Vs, damping) |
| 3 | `<dynprop>` | label of the dynamic soil property (DYNP) |

- Error 95 if no dynamic properties are assigned. Error 19 if the L layer is undefined.
- For nonlinear SSI, the SOIL and SITE soil properties must be consistent.
- **SASSI-EDU** (requirements §7.19, D-W5-09/10, Warning EDU-29): the labels Clay, Sand and Rock resolve to built-in
  SHAKE91 curves when no DYNP of the model defines them; with no SPRO entry at all SOIL uses a default profile
  (the TOPL layers with Sand / Rock by Vs, the SITE half-space last).
- **Example:** `SPRO,1,2,Sand` / `SPRO,2,2,Sand` / `SPRO,3,3,Clay`

### 9.2.34 SRS — SOIL response-spectrum output  [IO]
**Syntax:** `SRS,<layer>,<save>,<outcrop>`. `<layer>` = layer (ARS at the **top** of the layer). `<save>` 0/1. `<outcrop>` 0/1.
- RS are computed for each DAMP ratio.
- **Example:** `SRS,1,1,0`

### 9.2.35 SSAF — spectral amplification factor output (SOIL)  [IO]
**Syntax:** `SSAF,<layer>,<save>,<outcrop1>,<outcrop2>,<layer2>,<freqstep>,<title>`

| # | Arg | Meaning |
|---|---|---|
| 1 | `<layer>` | first sublayer (numerator) |
| 2 | `<save>` | 0/1 save |
| 3 | `<outcrop1>` | outcrop option for the first sublayer |
| 4 | `<outcrop2>` | outcrop option for the second sublayer |
| 5 | `<layer2>` | second sublayer number (Error 109 if illegal) |
| 6 | `<freqstep>` | frequency step, > 0 (Error 110) |
| 7 | `<title>` | title of the SAF output (rest of line) |

- *(inferred)* SAF(f) = RS_layer(f) / RS_layer2(f) at each DAMP ratio, evaluated on a grid with spacing `freqstep`.
- **Example:** `SSAF,1,1,0,1,12,0.1,Surface to bedrock outcrop SAF`

### 9.2.36 SSTR — SOIL stress/strain output  [IO]
**Syntax:** `SSTR,<layer>,<opt1>,<opt2>,<opt3>,<opt4>`. Stresses and strains are computed at the **centre** of the layer.
- `<opt1>`: compute stress (0/1). `<opt2>`: save stress history `SSxxx.TH` (0/1). `<opt3>`: compute strain (0/1). `<opt4>`: save strain history `SNxxx.TH` (0/1).
- **Example:** `SSTR,5,1,0,1,1`

### 9.2.37 STATUS — general information  [UI]
**Syntax:** `STATUS` (also `STAT`). Shows the global variables and general model information, including the active symmetry planes/lines (see SYMM).
- *(suggested content)* Model name, path and title. Counts of nodes, elements, groups, materials and layers. Units/gravity. Frequency sets. AOPT flags. Analysis-option summary. SYMM entries.

### 9.2.38 STRESS — STRESS options  [CORE]
**Syntax:** `STRESS,<opmode>,<iter>,<save>,<itran>,<interopt>`. No abbreviation is underlined (Open Q 1).

| # | Arg | Meaning | Values | Dialog |
|---|---|---|---|---|
| 1 | `<opmode>` | operation mode | 0 complete, 1 data check | Operation Mode |
| 2 | `<iter>` | automatic strains in all soil elements (iterative nonlinear soil SSI) | 1 yes, 0 otherwise | Auto Computation of Strains in Soil El. |
| 3 | `<save>` | save stress time histories in **`.ths`** files | 1 save, 0 no | Save Stress Time Histories |
| 4 | `<itran>` | output transfer functions (beam forces/moments → `.TFU`/`.TFI`) | 1 yes, 0 no | Output Transfer Function |
| 5 | `<interopt>` | interpolation option for stress TF | 0–6, same schemes as MOTION `<interp>` (§9.2.22) | Interpolation Option |

- **Not carried by this command:**
  - The STRESS tab's own Phase Adjustment and Smoothing values (separate from MOTION's in the screenshots) and Skip Time History Steps.
  - Post-processing/frame options, Section Cut "Save Time History", Save Binary Database.
  - See Open Q 16. The acceleration-history data is shared with MOTION, and the analysis type with ANALYS.
- STRESS works for one input direction at a time. COMB_XYZ_STRAIN combines X/Y/Z strains (SRSS) between nonlinear iterations [ADV].
- **Example:** `STRESS,0,0,1,0,0`

### 9.2.39 SYMM — symmetry / anti-symmetry plane or line  [CORE]
**Syntax:** `SYMM,<no>,[<type>],[<node1>],[<node2>],[<node3>]`

| # | Arg | Meaning |
|---|---|---|
| 1 | `<no>` | plane/line number, **1 or 2** (max 2) |
| 2 | `<type>` | 0 = symmetry, 1 = anti-symmetry, **with respect to the loading** |
| 3–5 | `<node1>…<node3>` | nodes defining the plane (3 non-collinear nodes, 3D) or line (2 nodes, 2D). `<node1>` = 0 resets plane `<no>`. |

- **Geometry rules:**
  - 3D: the planes must be parallel to the **xz or yz** plane.
  - 2D: the line must be parallel to the **z axis**.
  - Any combination of at most two planes/lines is allowed.
- **Checks:**
  - Error 2: node illegal or undefined.
  - Error 3: only one node.
  - Error 4: line nodes equal or coincident.
  - Error 5: plane nodes collinear.
- **Incoherency is not allowed** on symmetric (half or quarter) models (§6.5.4).
- **Implementation [CORE]**, standard symmetry boundary conditions *(inferred)* for a plane with normal n = x (a yz-parallel plane):
  - symmetric loading: u_x = 0, θ_y = θ_z = 0;
  - anti-symmetric loading: u_y = u_z = 0, θ_x = 0.
  - For a plane with normal n = y, use the analogous conditions.
  - The soil flexibility of a half-space cut by a symmetry plane needs image-source treatment (Open Q 22).
- **Example:** `SYMM,1,0,1,2,5` (plane through nodes 1, 2, 5) / `SYMM,1,0,0` (reset)

### 9.2.40 THFILE — acceleration time-history file  [IO]
**Syntax:** `THFILE,<file>` (also `THFI`). Sets the control-motion acceleration file used by SOIL, MOTION, STRESS and RELDISP ("File" in those tabs).
- Error 73 if it does not exist.
- Format per MOTION `<fopt>` and SOIL `<header>`/`<nrval>`. Values in g.
- **SASSI-EDU** (requirements §7.19, D-W5-03 to D-W5-06): `@<file>` names a built-in file; a blank THFILE is
  `@rg160h_030g.acc` (seismic, SOIL) or `@ricker_5hz.th` (ANALYS `<type>` 1) when `<fopt>` = 0 and Δt = 0.005 s,
  reported as Warning EDU-29; `EDUOPT,DEFAULTS,OFF` restores Error 73.
- **Example:** `THFILE,C:\SSI\Demo1\H1.acc`

### 9.2.41 THTIT — time-history title  [UI]
**Syntax:** `THTIT,<title>` (also `THTI`). Title of the acceleration history (shown as "Title" in MOTION/STRESS/RELDISP).
**Example:** `THTIT,acc_X_8192`

### 9.2.42 TPSD — target PSD file (EQUAKE)  [ADV]
**Syntax:** `TPSD,<num>,<file>` (duplicated as §9.7.28).
- `<num>` = spectrum number, the same index as ACCIN/ACCOUT/RSIN/RSOUT.
- `<file>` = full path of the target PSD.
- **File format:** two columns, frequency and target PSD amplitude. Units follow the gravity units: cm²/s³ (IS) or in²/s³ (BS). The p. 117 text prints "cm/sec^3, inch^2/sec^3"; see 05a.
- Active when `EQUAKE [tpsd]`=1. The PSD criterion is **not** used to accept seeds; only the `.psd` file is generated.
- **Example:** `TPSD,1,C:\SSI\Demo1\H1_target.psd`

### 9.2.43 TIT — model title  [UI]
**Syntax:** `TIT,<title>`. Sets the active model title. MDLNAME changes the file name, not the title.
**Example:** `TIT,Reactor building SSI model - BE soil`

### 9.2.44 TOPL — SITE top-layer list  [CORE]
**Syntax:** `TOPL,<l1>,[l2],…,[l200]`

- Appends the non-zero L-layer numbers to the free-field layer list above the half-space. **The first entry is the topmost layer.** `<l1>` = 0 deletes the list.
- The same L number may be repeated for identical layers, **but not for embedment layers.**
- **Embedded models:** the L numbers of the embedment layers must be the far-field layer numbers, in the same order as the TOPL list. The excavated-soil L material written to `.hou` depends on this.
- Dialog separators: blank, tab, `,`, `;`, Enter.
- **Limits:**
  - Max 200 entries (UI and §1.2).
  - Warning 8 says only the first **100** are written (Open Q 12).
  - Error 46 if the list is empty. Error 19 if a layer is undefined.
- **Example:** `TOPL,0` then `TOPL,1,2,3,4,5,6,6,6,7,7`

### 9.2.45 WAVE — wave field data (SITE)  [CORE]
**Syntax:** `WAVE,<type>,<opt>,<ratio1>,<ratio2>,<angle>`

| # | Arg | Meaning | Values |
|---|---|---|---|
| 1 | `<type>` | wave type | 1 R, 2 SV, 3 P, 4 SH, 5 L |
| 2 | `<opt>` | wave field option | 0 none; 1 wave field (**for R: shortest wavelength method**); 2 least-decay method (**R only**) |
| 3 | `<ratio1>` | participation fraction at frequency 1 | 0 < r ≤ 1 (Error 54) |
| 4 | `<ratio2>` | participation fraction at frequency 2 | 0 < r ≤ 1 |
| 5 | `<angle>` | incident angle (deg), from the propagation direction to the z′ axis (0 = vertical) | Error 52 if outside "(0,360)". 0 must be allowed. |

- **Rules:**
  - Types 1–3 apply with `SITE <wopt>`=0. Types 4–5 apply with `<wopt>`=1.
  - At each of the two frequencies, the ratios of all participating waves must **sum to 1** (Error 55).
  - For a single wave type, use ratios 1 and 1.
  - Intermediate frequencies use linear interpolation *(clamped outside the range, inferred)*.
- **Example:** `WAVE,2,1,1.0,1.0,0.0` (vertical SV only), or `WAVE,4,1,1.0,1.0,0.0` with `wopt`=1 (vertical SH)

### 9.2.46 WPASS — wave passage data (HOUSE)  [ADV]
**Syntax:** `WPASS,<appv>,<ang>,<cohf>` (also `WPAS`).

| # | Arg | Meaning | Values |
|---|---|---|---|
| 1 | `<appv>` | apparent (horizontal) velocity along Line D | > 0 (Error 113). Use ~1.0e9 to switch off wave passage effects (EPRI rock-site practice). |
| 2 | `<ang>` | angle of Line D with the x axis (deg) | also rotates the axes of directional coherency models 2–7 |
| 3 | `<cohf>` | "directional coherence factor" | ≥ 0 (Error 114). **Proposed:** the dialog "Unlagged Coherency Model" type 1–7 (Open Q 23) |

- Unlagged Coherency Model types (dialog):
  1. 1986 Luco-Wong
  2. 1993 Abrahamson, all soils (similar to 2005)
  3. 2005 Abrahamson, all sites, surface foundations
  4. 2006 Abrahamson, all sites, embedded foundations
  5. 2007 Abrahamson, hard rock, surface/embedded
  6. 2007 Abrahamson, soil sites, surface foundations
  7. user-defined, from the files `COHXUSER`, `COHYUSER`, `COHZUSER`, `FREQCOH`, `DISTCOH` (default 100×100)
- Lagged coherency: γ_lagged(ω, i, j) = γ_unlagged · exp(−iω(τ_i − τ_j)), with τ = (distance along Line D)/appv *(standard, see 05b)*.
- **Example:** `WPASS,1.0E9,0.0,3`

### 9.2.47 WRITE — write the model to a `.pre` input file  [IO]
**Syntax:** `WRITE,[<file>],[<path>]` (also `WRIT`).

- Writes every existing datum as instruction lines, so that INP rebuilds the model.
- Default file name: **model name + `.pre`**. Default path: the model path.
- This is an ASCII file of PREP instructions. It is the archive or exchange format between PREP, the UI and versions. Menu equivalent: Model ▸ Output.
- **Options ▸ Write ("Extended Write Options", p. 114):**
  - `MDL command in *.Pre`.
  - `Extend integer fields (for Models with more than 100000 nodes)`, IKTR8_300K solver.
  - `AFWR Command in *.Pre`.
  - `Simulation Commands`, with `Simulation Command Location`: `Write to *-Sim.Pre` (default) or `Write to *.Pre`.
  - With nothing checked, the output is identical to old PREP.
- **Round-trip requirement [CORE]:** WRITE followed by INP into an empty model must reproduce an identical state. This is a key test.
- **Suggested emission order** *(inferred, "logical sequence")*:
  1. `TIT`, optional `MDL`, then `MOPT`.
  2. Coordinate systems, nodes, BCs, INT; groups, elements; M, L, R, SC, MX*; MSET/RSET.
  3. Masses, forces; SYMM.
  4. FREQ sets.
  5. DYNP, SPRO, SACC, SRS, SSTR, SSAF, SFOU, DAMP.
  6. EQTIT, EQUAKE, RSIN, RSOUT, ACCIN, ACCOUT, TPSD, CORR; THFILE, THTIT.
  7. SOIL; SITE, TOPL, WAVE; POINT.
  8. HOUSE, INCOH, WPASS, ME, AMP; FORCE; ANALYS.
  9. MOTION, NOUT; STRESS, EOUT; RELD, RELFILE, RDND.
  10. AOPT, and an optional trailing `AFWR`.
- **Example:** `WRITE` / `WRITE,Archive_rev2.pre,C:\SSI\Archive`

### 9.2.48 RELD — RELDISP options  [IO]
**Syntax:** `RELD,<RelDisOutput>,<RelDispSAll>,<RelDispNumFiles>`. The syntax spells the 1st argument `RelDisOutput` (sic).

| # | Arg | Meaning |
|---|---|---|
| 1 | `<RelDisOutput>` | user-requested output type. Dialog "Output Control": complex vs amplitude relative-displacement TF, and whether `.TFI` or `.TFU` is used (codes undocumented, Open Q 24) |
| 2 | `<RelDispSAll>` | flag. When the `.rdi` is written, it overrides the output node list with **all node components** and tells RELDISP to **make displacement frames**. Dialog: Save Relative Displacement in All Nodes / Restart for Frame Generation. |
| 3 | `<RelDispNumFiles>` | number of entries ("file names") in the output node list, i.e. the count of RDND entries |

- The reference TF comes from RELFILE. NFFT/Δt/scaling/file are shared with MOTION. RELDISP needs **complex** `.TFI` from MOTION (`MOTION <cplx>`=1).
- For a free-field reference, the user builds a reference TFI with amplitude 1 in the input direction, 0 otherwise, and zero phase (coherent).
- One RELDISP run gives one DOF (TR_X, TR_Y, TR_Z, R_XX, R_YY or R_ZZ), so six runs give all DOFs.
- **Example:** `RELD,0,0,2`

### 9.2.49 RELFILE — reference-node TFI file (RELDISP)  [IO]
**Syntax:** `RELFILE,<FileName>`. **Only the full name is accepted**: the whole name is underlined.
- Sets the reference node's complex `.TFI` (dialog "Reference Location and Direction ▸ Complex TF File Name").
- Naming follows MOTION: `xxxxxTR_X.TFI` etc.
- **Example:** `RELFILE,C:\SSI\M1\00415TR_X.TFI`

### 9.2.50 RDND — RELDISP output node  [IO]
**Syntax:** `RDND,<NodeNum>,<X>,<Y>,<Z>,<XX>,<YY>,<ZZ>`
- Adds node `<NodeNum>` to the RELDISP output list.
- Each DOF flag: **< 1 → ignored; ≥ 1 → considered.** Dialog "Add/Edit Node" check boxes: X, Y, Z, XX, YY, ZZ.
- **Example:** `RDND,15,1,0,1,0,1,0` (X, Z, YY)

---

## 5. Related dialogs that define command defaults and behaviour

### 5.1 Options ▸ Model (MOPT)
Groups: `Incompatible Modes` (◯ Include / ◉ Suppress), `General Matrix` (◉ Mass Matrix / ◯ Weight Matrix), `Overwrite Mass` ☑, `Overwrite Force` ☑. Buttons OK, Cancel. ◉ = shown selected in the screenshot.

### 5.2 Options ▸ Check (CHECK)
`Show Warnings` ☑, `Show Errors` ☑, `Suppress Error Window` ☐, `Break Check at [100] Messages`. Session-only.

### 5.3 Options ▸ Write (WRITE/AFWRITE)
`MDL command in *.Pre` ☐, `Extend integer fields (for Models with more than 100000 nodes)` ☐, `AFWR Command in *.Pre` ☐, `Simulation Commands` ☐; `Simulation Command Location` ◉ `Write to *-Sim.Pre` / ◯ `Write to *.Pre` (greyed until Simulation Commands is checked).

### 5.4 Options ▸ Analysis ▸ AFWRITE tab (AOPT)
Check boxes: EQUAKE, SOIL, LIQUEF (greyed), SITE, POINT, HOUSE, PINT (greyed), FORCE, ANALYS, COMBIN, MOTION, STRESS, RELDISP, NONLINEAR.

### 5.5 Field → command map, condensed (full tab specs in 05a–05d)
| Tab | Field | Command.arg |
|---|---|---|
| EQUAKE | Spectrum Number / Input / Output / Acceleration Output File | `<no>` of RSIN / RSOUT / ACCOUT |
| EQUAKE | Accel. Record; Acceleration Input File | EQUAKE `<accopt>`; ACCIN |
| EQUAKE | Use Target PSD; PSD File | EQUAKE `[tpsd]`; TPSD |
| EQUAKE | Number of Frequencies, Initial Random SEED, Damping Value, Total Duration, Number of SEEDs, Correlated | EQUAKE args 2–7 |
| EQUAKE | Time Step | SITE `<delt>` |
| EQUAKE | Correlation grid; Spectra Title | CORR; EQTIT |
| SOIL | Number of Values, Gravity, Header Lines, Assign as Outcrop, Save Strain-Compatible, Iterations, Equiv. Uniform/Max Strain, Multiplier for g | SOIL args 1–8 |
| SOIL | Layer/Property/Dynamic Property | SPRO |
| SOIL | Accelerations / Response Spectrum / Stresses_Strains / SAF / Fourier groups | SACC / SRS (+DAMP) / SSTR / SSAF / SFOU |
| SITE | Mode 1/2, Frequency Step, Δt, NFFT, Frequency Set, Generated Layers, Halfspace Layer, wave option, Frequency 1/2, Control Point Layer, Direction | SITE args |
| SITE | Top Layers; wave tabs | TOPL; WAVE |
| SITE | Gravity | HOUSE `<gravity>` |
| POINT | Operation Mode, Number of Embedment Soil Layers, Point Load Central Zone Radius | POINT |
| HOUSE | Operation Mode, Dimension, Method, Gravity, Ground Elevation, Soil Motion, Wave Passage, Multiple Excitation, Complex Spectral Amp. | HOUSE |
| HOUSE | Coherence X/Y/Z, Alpha, Embedded Layers, Incoh. Modes, Print, simulation radio | INCOH |
| HOUSE | Apparent Velocity, Angle Line D, Unlagged Coherency Model | WPASS |
| HOUSE | Input Motion Number, First/Last Node, X/Y/Z Coord.; Spectral Amplification | ME; AMP |
| FORCE | Operation Mode | FORCE (the rest is shared) |
| ANALYS | all except FFL/FFM and coherence/WP/ME mirrors | ANALYS |
| MOTION | see 9.2.22 table; node lists | MOTION; NOUT |
| STRESS | Operation Mode, Auto Strains, Save Stress TH, Output TF, Interpolation; element lists | STRESS; EOUT |
| RELDISP | Complex TF File Name; Output Control; Save in All Nodes/Frames; node list | RELFILE; RELD; RDND |
| AFWRITE | module check boxes | AOPT |

---

## 6. Validation catalogue for this command set (CHECK)  [CORE][UI]

These are the Chapter 10 messages that apply to the §9.2 commands. Every module check runs only for modules enabled by AOPT.

| No. | Condition | Source command(s) |
|---|---|---|
| E1 | gravity ≤ 0 | HOUSE `<gravity>` (and SOIL `<grav>`) |
| E2–E5 | SYMM node illegal/undefined; one node only; line nodes equal; plane nodes collinear | SYMM |
| E44 | frequency set not defined | SITE `<freq>`, FREQ |
| E45 | SITE Mode 1 and Mode 2 both off | SITE |
| E46 | no top layers | TOPL |
| E47 | half-space layers not 0 and not in 4…20 | SITE `<nl>` |
| E48/E49/E50 | Δf / Δt / NFFT negative | SITE |
| E51 | all wave fields off | WAVE |
| E52 | incident angle outside (0,360) | WAVE `<angle>` |
| E53 | ratio-curve frequency ≤ 0 | SITE `<freq1>,<freq2>` |
| E54 | wave ratio ≤ 0 or > 1 | WAVE |
| E55 | sum of ratios ≠ 1 at a frequency | WAVE |
| E56 / E57 | POINT layer < 0 / radius ≤ 0 | POINT |
| E58 | coherence parameter < 0.1 | INCOH |
| E59 | Luco-Wong mean Vs ≤ 0 | INCOH `<alpha>` |
| E60 | mesh points per embedment level ≤ 0 | INCOH `<ngp>` |
| E61 / E62 | no forces / force on undefined node | FORCE (F, MM) |
| E63 | coordinate transformation angle outside (0,360) | ANALYS `<ang>` |
| E64–E66 | no nodal request / illegal node / duplicate node | NOUT |
| E67 / E68 | output step < 0 / duration < 0 | MOTION |
| E69–E71 | RS first/last frequency or number of steps < 0 | MOTION |
| E72 | RS damping not in (0,1) | DAMP |
| E73 | acceleration file missing | THFILE |
| E74–E76 | first record illegal / last record < 0 / first > last | MOTION `<rec1>,<rec2>` |
| E77 / E78 | mult and max both zero / both non-zero | MOTION |
| E79–E82 | element request missing / illegal group / illegal element / duplicate | EOUT |
| E84–E89 | EQUAKE RS files (none, missing, output blank, acceleration output blank, acceleration input blank, record count ≠ nrfreq) | RSIN, RSOUT, ACCOUT, ACCIN, EQUAKE |
| E90–E92 | seed ≤ 0 / nrfreq ≤ 0 / duration ≤ 0 | EQUAKE |
| E93 / E94 | correlation on but no pairs / factor > 1 | CORR |
| E95–E99 | no dynamic properties assigned / > 15 used / label lacks G curve / G outside [0,1] / no D curve | SPRO, DYNP |
| E100–E108 | SOIL: nrval ≤ 0, cof < 0, illegal format, header < 0, illegal control layer, iterations < 0, ratio outside (0,1), no damping ratios, gravmult ≤ 0 | SOIL, DAMP, SITE `<cl>` |
| E109–E112 | SSAF layer2 illegal / freqstep ≤ 0; SFOU smoothings < 0 / values < 0 | SSAF, SFOU |
| E113 / E114 | apparent velocity ≤ 0 / directional coherence factor < 0 | WPASS |
| E115–E119 | ME on with no data / first or last node illegal / SAR outside [0,10] / SAR count ≠ number of frequencies | ME, AMP |
| E120 | frequency set empty (RELDISP) | FREQ |
| W1, W4 | gap or unused node → fixed at AFWRITE | (model) |
| W7 | empty group skipped | (model) |
| W8 | > 100 top layers → first 100 written | TOPL |
| W9 | NFFT not a power of 2 → nearest written | SITE `<nft>` |
| W10/W11 | force/moment on a fixed DOF ignored | FORCE |

**Recommended extra checks** *(inferred, flag as warnings)*:
- More than 5 DAMP values (MOTION limit).
- EQUAKE `<dur>` < 20 s or Δt > 0.005 s (SRP 3.7.1).
- MOTION RS steps < 301.
- `smo` ≠ 0 with `interp` = 6.
- `smo` ≠ 0 with `<coh>` = 0.
- Simultaneous seismic cases with `<ang>` ≠ 0.
- `<coh>`=1 with SYMM defined or `<dim>`≠2.
- Abrahamson/user models or `<me>`=1 with `<wpass>`=0.
- RELDISP enabled with MOTION `<cplx>`=0.
- Frequency numbers > NFFT/2.
- More than 500 SSI frequencies.

---

## 7. Worked example: analysis-option block of a `.pre` file (illustrative)

```
* ---- model identification
TIT,Demo surface mat - coherent, BS units
MDL,DemoMat,C:\SSI\DemoMat
MOPT,1,0,1,1
* ---- frequency set 1 (df = 1/(0.005*4096) = 0.048828 Hz)
FREQ,1,0
FREQ,1,2,4,6,8,10,12,14,16,20,24
FREQ,1,28,32,40,48,56,64,80,96,128,160
FREQ,1,192,256,320,384,448,512
* ---- free-field (SOIL) equivalent-linear site response
DYNP,1,0.0001,1.0,0.0001,0.5,Sand
DYNP,2,0.001,0.98,0.001,0.8,Sand
DYNP,3,0.01,0.85,0.01,2.5,Sand
DYNP,4,0.1,0.45,0.1,9.0,Sand
DYNP,5,1.0,0.10,1.0,20.0,Sand
SPRO,1,1,Sand
SPRO,2,2,Sand
SACC,1,2,0
SRS,1,1,0
SSTR,2,1,0,1,0
DAMP,0
DAMP,0.02,0.05
THFILE,C:\SSI\DemoMat\H1.acc
THTIT,H1 0.3g
SOIL,4000,32.2,0,1,1,8,0.65,1.0,0
* ---- linearized site response (SITE) and point loads (POINT)
SITE,0,1,0,20,3,1,0,1,2048,1,0,0.005,4096,1
TOPL,0
TOPL,1,2
WAVE,2,1,1.0,1.0,0.0
POINT,0,0,2.7
* ---- structure (HOUSE) and SSI solution (ANALYS)
HOUSE,32.2,0.0,0,2,0,0,0,0,0
ANALYS,0,0,0,1,1,0,0.0,0.0,0.0,0.0,0,0
* ---- post-processing
MOTION,0,0,0,25.0,0,0.1,100.0,301,1.0,0.0,1,4096,0,0,0,1,0,0,0
NOUT,1,1,1,0,0,1,1,101-105
NOUT,3,1,1,0,0,1,1,101-105
STRESS,0,0,1,0,0
EOUT,1,1,1,1,1,1,1,0,0,0,0,0,2,1-40
RELFILE,C:\SSI\DemoMat\00101TR_X.TFI
RDND,105,1,0,0,0,0,0
RELD,0,0,1
* ---- write and check decks
AOPT,0,1,0,1,1,1,0,0,1,0,1,1,1,0
AFWRITE
```
The DYNP values above are illustrative placeholders, not a published curve.

---

## 8. Implementation notes and verification hooks

1. **Parser unit tests [IO]:**
   - Every full name and abbreviation in §1.3 resolves. Every other prefix (e.g. `ACC`, `ANALY`, `RELF`) is rejected with "Command not found" and processing continues.
   - Names are case-insensitive; file paths keep their case.
   - Parentheses grouping works for legacy SOIL.
   - Titles may contain commas.
2. **Round trip [IO]:** WRITE then INP gives an identical state (deep equality of the analysis-option records and the lists FREQ, DAMP, TOPL, AMP, NOUT, EOUT and RDND).
3. **Frequency derivation [CORE]:** Δf = 1/(Δt·NFFT); FREQ numbers are sorted; duplicates are rejected by SITE; f = n·Δf. Test: Δt = 0.01, NFFT = 2048 → Δf = 0.048828125 Hz exactly.
4. **List semantics:** `DAMP,0.02` then `DAMP,0.05` gives [0.02, 0.05]. `DAMP,0` clears. Same for FREQ (per set), TOPL and AMP (per motion).
5. **AFWRITE gating:** a deliberately invalid SITE (e.g. `<nl>`=3) blocks `.sit` only. The other enabled decks are still written.
6. **Scaling:** MOTION mult/max exclusivity, and the max-scaling formula, are unit tested on a synthetic history.
7. **Shared variables:** changing `SITE <delt>` changes the Δt seen by SOIL, MOTION, STRESS, RELDISP and EQUAKE writers. `GRAVITY,9.81` changes only HOUSE arg 1.

---

## Open questions / ambiguities

1. **STRESS abbreviation.** No underline is printed for STRESS, although the name has more than 4 letters. Old PREP decks may contain `STRE`. Proposal: accept `STRESS` and `STRE`, and log a compatibility note for `STRE`.
2. **List syntax inside commands.** Ranges (`1-20`) and blank/`;` separators are documented only for dialog list boxes. Proposal: in NOUT/EOUT the tokens after the fixed args may be single numbers or `a-b` ranges, separated by commas or blanks.
3. **Titles with commas** (TIT, THTIT, EQTIT, SSAF `<title>`). The manual does not say. Proposal: the rest of the line is the title. For SSAF, `<title>` is everything after the 6th comma.
4. **EQUAKE `<accopt>`** has values 0/1, but the dialog has two independent options: Accel. Record (seed) and External Accel. (RS only). Proposal (same as 05a): 0 none, 1 seed record, 2 external. Otherwise add an argument.
5. **DYNP damping units.** The table on p. 206 suggests damping in percent and strain in percent (SHAKE convention). DAMP uses fractions. Confirm and document the units in the UI.
6. **SOIL fields without arguments.** Input Direction (0/1), Multiplication Factor, Max Value (in g), File, Control Point Layer, NFFT, Δt. Proposal: NFFT, Δt and control layer come from SITE; File from THFILE. Mult/Max: decide whether SOIL shares MOTION's or keeps its own pair (the screenshots show different values: 0/0.1 in SOIL, 1/0 in MOTION). Input Direction needs a new stored field. If so, add optional trailing SOIL arguments in our dialect and keep legacy compatibility.
7. **INP nesting and path resolution.** The manual does not say whether an INP inside a `.pre` returns to the calling file or to the keyboard, or which directory relative paths inside a nested file use. Proposal: use a stack of input sources, return to the caller at EOF, and resolve relative paths against the active model path (as documented). Add a recursion-depth guard.
8. **AMP append vs replace**, and the **complex SAR syntax** when `<cmplxspec>`=1 (real/imag pairs? amplitude/phase?). This is undocumented. Proposal: append semantics; with cmplxspec=1 read values as consecutive (Re, Im) pairs, so 2·Nfreq numbers.
9. **Deleting output requests** (NOUT, EOUT, RDND) by command is not documented. Proposal: `NOUT,0` / `EOUT,0` / `RDND,0` clears that request list (consistent with FREQ/DAMP), and WRITE emits the clear first.
10. **ANALYS `[simul]` default 0 vs dialog 1.** Proposal: treat 0 and 1 the same for vibration and incoherent runs. For coherent seismic, interpret ≥ 1 as the three-direction FILE1X/Y/Z run only when those files exist? This needs a decision. The manual says coherent simultaneous should be 1. Safer: add an explicit flag in our UI.
11. **Dialog options with no command argument.**
    - HOUSE: Optimize Model, Superposition Mode (Linear/Quadratic), number of stochastic simulations, Non-Linear SSI `.pin`, ANSYS Model Input.
    - ANALYS: FFL/FFM.
    - MOTION: Save FILE12/13, post-processing and binary options.
    - STRESS: Phase Adjustment, Smoothing, Skip TH Steps.
    - Some are probably set by other §9 commands (BINOUT, ANSYSMODELTYPE…). The rest need our own extension commands. **Decision needed on dialect extensions** (e.g. optional trailing args).
12. **TOPL limit:** 200 (UI, §1.2) vs Warning 8 "only the first 100 are written". Proposal: allow 200 and warn above 100 only when writing in "legacy" format.
13. **Dynamic property count:** Error 96 says > 15 used is illegal; §1.2 says 100 soil material curves, and 11 data points per curve. Proposal: enforce 100 curves and 11 points, and warn above 15 (legacy).
14. **INCOH `<ngp>`:** "mesh points per embedment level" vs dialog "Number of Embedded Layers" (screenshot 0, but Error 60 forbids ≤ 0). Its algorithmic role is unclear. Proposal: store the value. Use it as the number of interaction-node levels for per-level coherency; 0 → auto (surface), and relax E60.
15. **ME motion count:** the command says 1–10; the limits say up to 5000 foundation zones. Proposal: allow ≥ 1 with no hard maximum of 10, and warn above 10.
16. **MOTION `<step>` and `<res>`** have no visible dialog field. STRESS "Skip Time History Steps" may share `<step>`. Proposal: `<res>` is reserved (write 0). `<step>` is the print step shared with STRESS.
17. **MOTION `<bl>`:** the command says 0 = time domain, 1 = frequency domain; the dialog says No Correction / With Correction (Hudson-Housner time domain), and FILE13 exists only with correction. Proposal: 0 = no correction, 1 = correction (time domain, Hudson-Housner), so the dialog wins. Flag it in the docs.
18. **Duplicate NOUT nodes:** Error 66 vs automatic deletion at AFWRITE. Proposal: CHECK reports a warning and AFWRITE de-duplicates (merging code flags by OR).
19. **SITE `<opmode>`:** the command says complete/data check; the dialog says Linear Soil / Non-Linear Soil (Non-Linear reads FILE88). Proposal: keep `<opmode>` = complete/check, and add our own flag for Linear/Non-Linear soil. Alternatively, find in the `.sit` format how FILE88 use is selected.
20. **SITE `<freq1>,<freq2>`:** Hz or frequency numbers? The recommended "1" and "≥ NFFT/2" and the screenshot 1/4000 imply frequency numbers (index × Δf). Proposal: interpret them as frequency numbers.
21. **Old (legacy) SOIL argument order.** It is not given; the manual says only that the arguments changed and that the old form ends with a parenthesised format flag. Legacy decks can be detected but not mapped reliably. Proposal: warn and refuse to auto-map, or allow a user-supplied mapping.
22. **Symmetry handling in the soil (POINT/ANALYS):** the manual gives only the geometric rules. It does not say how the layered-soil flexibility of a half or quarter model is formed (image superposition). This must be designed and verified against a full model.
23. **WPASS `<cohf>`:** "directional coherence factor" (E114: ≥ 0) vs the dialog field "Unlagged Coherency Model" (1–7) in the same group. Proposal (same as 05b): `<cohf>` = model type. α stays in INCOH `<alpha>`.
24. **RELD `<RelDisOutput>` codes** (complex vs amplitude, TFI vs TFU) and the exact effect of `<RelDispSAll>` (all nodes **and** frames) are not enumerated. Proposal: bit 0 = complex output, bit 1 = use TFU instead of TFI. Store the raw integers for round trip.
25. **Angle range checks (E52, E63)** say "(0,360)", but 0° is the normal vertical / no-rotation value. Proposal: valid range is [0,360). Negative angles are normalized with a warning.
26. **NFFT "closest/nearest power of 2"** rounding rule (up, down or nearest) is unspecified. Proposal: round **up** to the next power of 2, which keeps all records, and warn.
27. **RS frequency spacing** for MOTION `<freq1>…<freq2>` with `<fstep>` steps (log or linear) is unspecified. Proposal: logarithmic spacing (standard ISRS practice), with the endpoints included.
28. **Record definition** for `<rec1>/<rec2>` (one value vs one line) when `<fopt>`=0 or 1. Proposal: record = one acceleration sample after the header/time-step line.
29. **EQUAKE `<seeds>`**: the UI text requires a non-zero integer, but the screenshot shows 0. Proposal: 0 is treated as 1.
30. **Interpolation between CORR pairs**, and between and beyond the WAVE ratio frequencies, is not specified. Proposal: linear interpolation, held constant outside the range.
31. **SOIL British units for unit weight:** p. 119 says lb/ft³, while §6.4.4 and SITE say kip/ft/ksf/kcf. The SITE text also requires SOIL and SITE units to be consistent for nonlinear SSI. Proposal: one model-wide unit system (BS = kip, ft, ksf, kcf; IS = kN, m, kN/m², kN/m³). Treat lb/ft³ input as an error unless it is explicitly converted.
