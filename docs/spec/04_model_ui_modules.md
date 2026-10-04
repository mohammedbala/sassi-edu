# 04 - Model Construction & Manipulation; UI Menus (Model, File, Plot, Modules)

Source: ACS SASSI V3 User Manual, Chapter 5 (all) and Chapter 6 intro + Sections 6.1-6.4
(printed pages 74-112, PDF pages 76-114; text lines 3385-4966 of `reference/acs-sassi.txt`).
Dialog screenshots on PDF pages 76, 88-91, 93-94, 99-103, 110-113 were viewed and transcribed.
Command signatures quoted here for context were checked against Chapter 9 (Sections 9.2, 9.4,
9.7, 9.9-9.15). The authoritative spec for each command is the Chapter 9 spec file; this file
fixes the *workflow / semantics / UI* layer that ties them together.

Tag legend:

- **[CORE]** needed for computational correctness
- **[IO]** file or input format
- **[UI]** user interface, plotting, convenience
- **[ADV]** advanced option (Option A / AA / NON / PRO, incoherency, nonlinear)
- **(derived)** means the implementer rule is inferred from the manual or from ANSYS element
  documentation, not stated word-for-word in the ACS SASSI manual. Treat it as a
  recommendation and verify it.

Command notation follows Chapter 9 (Backus-Naur style): `<x>` = required, `[x]` = optional,
arguments are comma separated, command names are case-insensitive.

---

## 0. Big picture: how this chapter fits the program

```
            +-------------------- ACS SASSI UI (the "pre/post processor") --------------------+
            |  command interpreter  <- Command Entry line, .pre files (INP), macros, loops       |
            |  N models in memory (model numbers 0,1,2,...; one "active" model, ACTM)            |
            |  per model: name+path (MDL), nodes, groups/elements, materials, layers, props,     |
            |             masses, BCs, interaction nodes, model options, analysis options, freqs |
            |  session-wide: variables (VAR), macros (LOADMACRO), cuts (CUTADD/CUTVOL...)        |
            |  persistence: SAVE/RESUME (binary), WRITE (.pre text), SASSIdb.xml, SASSIini.xml   |
            |  converters: .hou -> model, ANSYS .cdb -> model ; exporter: model -> ANSYS .inp    |
            |  AFWRITE -> module input files (.equ .soi .sit .poi .hou .frc .anl .mot .str .rdi)  |
            +-------------------------------------|--------------------------------------------+
                                                  v  RUNxxx / Modules menu (exe path from SASSIini.xml)
   EQUAKE  SOIL  SITE -> FILE1,FILE2   POINT(FILE2) -> FILE3   HOUSE -> FILE4 (.N4) [+FILE77]
   FORCE -> FILE9   ANALYS(FILE1,FILE3,FILE4[,FILE9][,FILE77]) -> FILE8   COMBIN(FILE81,FILE82)->FILE8
   MOTION(FILE8) -> TFU/TFI/ACC/RS [FILE13]   RELDISP(.TFI) -> TFD/THD frames   STRESS(FILE4,FILE8) -> FILE14, FILE15
   Option NON: NONLINEAR   Option A: LOADGEN   Option AA: SSI2ANSYS, HOUSEFSA, ANALYSFA
```

---

## 1. Section 5.1: Recommended workflow for building a new SSI model [CORE][UI]

The manual gives this ordered checklist. An implementation should support it end-to-end and may
offer it as a guided "New model wizard" or checklist panel (UI choice).

| Step | Action | Commands / menus |
|---|---|---|
| 1 | Define structure nodes, **starting from the bottom** (recommended practice) | `N` |
| 2 | Set interaction nodes (automatic or manual); add fixed boundary conditions if needed | `INTGEN` or `INT`; `D` |
| 3 | Define groups and elements | `GROUP`, `E` |
| 4 | Check elements visually | Plot > Model > Elements (`MODELPLOT`) |
| 5 | Define materials, soil layers, beam real properties | `M`, `L`, `R` |
| 6 | Assign materials/properties to elements | `MSET`, `RSET` |
| 7 | Check assignments: color elements by material/property and compare with the color legend | Options > Colors or `COLOR`; plot toolbar color-by-group/material/property (`ELECOLOR`) |
| 8 | Define nodal masses and mass units | `MT`, `MR`, `MUNITS`. If General (GM) elements carry the mass via `MXM`, the units are set by `MOPT` (Options > Model "General Matrix: Mass Matrix / Weight Matrix") |
| 9 | Check masses | Plot > Model > Nodes (node plot shows lumped masses, `SHOWMASS`) |
| 10 | Set model options | Options > Model (`MOPT`) |
| 11 | Define analysis frequency set | `FREQ` |
| 12 | Set analysis options for the modules | Options > Analysis |
| 13 | Write module input files | `AFWRITE` (needs model name/path, see Section 3) |
| 14 | Run the SSI modules | Modules menu or `RUNxxx` commands |

The order of steps 5-6 vs 3 is not strict. For example `MACT`/`RACT` set "active" indices so that
elements defined later inherit them.

---

## 2. Models in memory, active model, Model > New [CORE][UI]

- The UI holds **several models in memory at once**, each identified by a non-negative integer
  **model number**. At startup the active model is **model 0**. The first plot tab in the
  screenshots is titled `Model 0 - Model Plot`.
- `ACTM,<Model>` makes `<Model>` active. If it does not exist, an empty model is created with
  no nodes or elements and default simulation settings. If it exists, its data is unchanged.
  All later commands act on the active model.
- **Model > New** [UI]: creates an empty model with the **lowest model number not currently in
  use** and makes it active. The manual only says "a number not currently in use"; "lowest" is
  a (derived) choice. It never deletes or changes other models in memory.
- Related commands (Chapter 9): `CPMODEL,<Mdl>` copies the active model to `<Mdl>`;
  `DMODEL,<Mdl>` removes a model from memory only (files are untouched); `MODELLIST` lists the
  models.
- Plot tabs are per model ("Model N - ..."). Commands that build a model into a destination
  (MERGE, MERGESOIL, CSECT, CUT2SUB, EXCAV, converters) take model numbers as arguments.

---

## 3. Sections 5.2-5.4 and 6.1.2-6.1.3: Model identity, saving/loading, database [IO][UI]

### 3.1 Model identity: name + path (`MDL`) [IO]

- A model has a **model name** and a **model directory path**. These determine where SAVE,
  RESUME, AFWRITE, WRITE (default name), RUNxxx and Export to ANSYS read and write files.
- A model that came only from a `.pre` file (INP / Model > Input) or a converter has **no
  name/path**. In that state SAVE, RESUME and AFWRITE "may run" but the output location is
  undefined. The implementation should **refuse with a clear error** (recommended) or use a
  documented fallback directory.
- `MDL,<Model>,<Path>` sets the name and path. It may be issued before or after the model is
  input. After MDL, SAVE/RESUME/AFWRITE work normally. Files written by WRITE/AFWRITE are named
  `<Model>.<ext>` in `<Path>`.
- `MDLNAME,<name>` changes only the name (path and title `TIT` unchanged).
- Options > Write offers "MDL command in *.Pre". When checked, WRITE emits an `MDL` line so the
  `.pre` file is self-locating (see the Options spec, manual 6.5.2).

### 3.2 SAVE / RESUME (binary) [IO]

- `SAVE` writes the active model's full data in a **binary database format** to the model
  directory, for later `RESUME`. It requires the name/path.
- `RESUME` reloads the binary data from the last SAVE. It discards changes made after the last
  save and requires the name/path.
- These do **not** update the database tree (SASSIdb.xml).
- The PREP (pre-IKTR4) binary format is **not** compatible with the UI binary format. Models move
  between the two programs only through `.pre` files (WRITE then INP).
- Implementation decision (open): the binary container. Recommendation: one versioned file such
  as `<name>.sdb` (zip of JSON + NumPy arrays, or HDF5) with a format-version header, so that
  RESUME can reject or migrate old versions.

### 3.3 The two-level database (Model > Open) [UI][IO]

- One database per user account, cached in a **system-determined location**. The user never has
  to browse for it. It is a convenience and is **never required**.
- It is stored in **`SASSIdb.xml`** at a default, OS-dependent location. It holds **only** the
  tree: group names and, for each model, its **name and location**. No model data is stored.
- Structure: **Group** (level 1, like a project folder) -> **Model** (level 2). Every model must
  be in a group. A model added with no group selected is **not** stored.
- **Load Model** dialog (Model > Open, shortcut **Ctrl+O**), transcribed from the screenshot:

| Element | Behavior |
|---|---|
| Tree view (left) | Groups as expandable nodes, with models as children. Example: group `Trial` holding `Test`, `teststruct`. It is empty the first time it is used. |
| `Open` | Load the highlighted model and close the dialog. Model files are loaded as if by RESUME, and the name and path are set. |
| `Cancel` | Close without action. |
| `Add Group` | Popup asking for a group name. The group is appended at the **bottom** of the tree. |
| `Remove Group` | Needs a highlighted group. (1) Confirm removal; if yes, the group leaves the tree. (2) A second popup asks whether to **delete all data of that group's models from disk**: yes removes the models' files and folders, no leaves them on disk. |
| `Add Model` | Needs a highlighted group, or a highlighted model, in which case the new model goes into that model's group. Opens a "model information" popup (fields not shown in the manual; see Open Questions). The model may be new (not yet saved), may already exist in another group, or may have been made elsewhere or on another computer. |
| `Remove Model` | Needs a highlighted model. Uses the same two-step confirm / delete-files popups as Remove Group. |
| Info pane (right, below buttons) | Text about the selection, e.g. `Group Name: Trial`. |

- Opening a **newly defined** model (never saved) produces "file not found" errors. This is
  **expected and harmless**: the database does not create model files, and the files appear at
  the first SAVE. The implementation may suppress these errors for new entries (UI choice) but
  must still allow opening.
- To register a model that was saved without the database (via MDL + SAVE), use Add Model with
  its location.

### 3.4 Working without the database [UI]

The manual's recommended path for a one-off analysis: `INP,<file.pre>` (or Model > Input, or a
converter), then `MDL,<name>,<path>`, then `AFWRITE`, then `RUNxxx`. To continue an earlier
saved model without the database: `MDL,<name>,<path>` then `RESUME`.

### 3.5 Model > Save (6.1.3) [UI]

The same as `SAVE`. Name and path must be defined first; otherwise show an error that suggests
MDL or the database.

### 3.6 Echo policy when reading `.pre` files [UI]

While a `.pre` file is executing, **do not echo every command**. Echo only **comments**
(lines starting with `*`), **warnings**, **errors**, and **model information text**. The manual
reports this cut load time from about 20 min to under 1 min. When the file ends, print
`INPUT FILE REACHED EOF, INPUT SWITCHED TO KEYBOARD` (seen in the screenshot).

---

## 4. Sections 5.5, 6.1.4, 6.2.1: PREP compatibility, .pre input, command entry, text editor

### 4.1 `.pre` compatibility [IO]

- Old PREP `.pre` files (up to IKTR4) must load in the new UI **unchanged**. The UI language is a
  superset of PREP.
- Any `.pre` written by WRITE is compatible with earlier UI modules, but a hand-written `.pre`
  may use new commands that PREP / SUBMODELER do not know.
- Unknown command: print `<X> Command not found` (where `<X>` is the name), **continue** with the
  next line, and do not abort. Later commands fail only if they depend on the missing one.
- PREP matched only the first 4 letters of a command name. The UI accepts either the full name
  or its documented abbreviation (the underlined part in Chapter 9) and nothing else. (Chapter 9
  rule, repeated here because INP must honor it.)
- Comment lines start with `*`. An empty argument (nothing, or whitespace only, between commas)
  means "use the default" (the `<blank>` convention in Section 5.8).

### 4.2 Model > Input (6.1.4) [UI][IO]

- A file dialog picks a `.pre` file, which is executed **on the current (active) model** as a
  batch of commands. Any valid UI command may appear, and a `.pre` can change any part of the
  model.
- Progress bar in the lower-right status bar: **current line / total lines** (not a time
  estimate). At the end, show a "loading completed" message.
- Same as `INP,<filename>`. The default path for `<filename>` is the active model's path. At EOF,
  input switches back to the keyboard.

### 4.3 Structure of a WRITE-generated `.pre` (from screenshots) [IO]

The header (File Editor screenshot, PDF p.99) is a banner of asterisks, then two comment lines
saying the file was written by the ACS SASSI SUBMODELER and can be reloaded with
`INP,<this file>`, then another asterisk banner. Section comment lines follow in this order (from
the command-history screenshot, PDF p.76):

```
* Nodes                 (N,...)
* Boundary Conditions   (D,...)
* Interaction Nodes     (INT,...)
* Material Table        (M,...)
* Soil Layer Table      (L,...)
* Real Property Table   (R,... SC,... MX*,...)
* Groups and Elements   (GROUP,..., E,..., MSET/RSET/KI/KJ/THICK/ETYPE,...)
* Masses                (MT,..., MR,..., MUNITS,...)
* Model Options         (MOPT,...)
* Analysis Options      (HOUSE, SITE, ANALYS, MOTION, ... option commands)
* Frequencies           (FREQ,...)
```

Example node lines: `N,1,27.65,-3,-16`. The section order shown is authoritative. The commands
listed in parentheses under each heading are (derived).

### 4.4 Command Entry and Command History [UI]

- **Command Entry**: a one-line input at the bottom of the main window. **Up/Down arrows** recall
  earlier commands while the caret is in the box.
- **Command History**: a tab showing output. It is copyable (select + copy). Message types each
  have a **user-selectable color**: command echo, output confirmation, comments, information,
  warnings and errors. View > Command Display filters them with toggles: Command Echo, Output
  Confirmation, Comments, Warning & Errors.

### 4.5 Text editor window (File > Open, 6.2.1; replaces PREP's read-only viewer) [UI]

- File > Open shows a file picker. The user picks a file **or types a new name, which creates the
  file**. A simple editor opens with title `File Editor - <full path>` and menus `File`
  (incl. **Save**) and `Input`.
- **Input > Connect to Command Entry**: while connected and open, every command typed in Command
  Entry that **does not produce an entry error** is appended to the editor text. **Only one window
  can be connected at a time.** This is how a user records a session into a reusable `.pre` or
  macro.

---

## 5. Section 5.6: Macros [CORE for the interpreter][UI]

### 5.1 Definition

- A **macro** is a specially formatted `.pre` file whose text holds **macro variables**: a
  positive integer between dollar signs, `$1$`, `$2$`, ..., `$10293$`.
- `LOADMACRO,<Name>,<MACROFILE>` binds a macro name to a file (full or relative path).
- `MACRO,<Name>,<1>,<2>,...,<N>` runs it, replacing `$k$` with the k-th argument after the name
  (1-based).
- `MACROLIST` lists the loaded macro names and their files.

### 5.2 Rules and limits (exact)

| Rule | Value |
|---|---|
| Name case | Names are converted to **UPPER CASE** (`Node` and `NODE` are the same macro) |
| Max macro name length | **50 characters** |
| Max length of the text substituted for one variable | **3000 characters** |
| Max length of one command line in a macro | **3000 characters** |
| Number of variables per macro | unlimited |
| Content | anything, **including calls to other macros** (nesting) |
| Name collision with a command | allowed. It does **not** override the command: macros live in their own namespace, reached only through `MACRO,<Name>` |
| Advice | do not call LOADMACRO inside a macro. Reloading the same file repeatedly hurts performance |

### 5.3 Execution algorithm (derived, consistent with all examples)

```
LOADMACRO(name, file):
    key = upper(name); check len(key) <= 50
    text = read_file(resolve(file))            # load once and cache (performance advice)
    macros[key] = (file, list_of_lines(text))  # reloading replaces the definition

MACRO(name, args[1..N]):
    m = macros[upper(name)] or error "macro not defined"
    for each line L in m.lines:
        L2 = replace every token  \$([1-9][0-9]*)\$  in L with args[k]   # plain text substitution,
                                                                       # also inside words
        check len(L2) <= 3000
        execute_command_line(L2)    # normal interpreter: comments, @variables, FOREACH, nested MACRO
```

- Substitution is **textual and substring-level**. The nested example uses `Node$1$x.rs`, which
  becomes `Node17x.rs`.
- Apply `$k$` substitution **before** `@variable` processing of the resulting line (derived).
- A placeholder with no matching argument: undefined in the manual. Recommendation: error
  "macro argument k missing" and skip the line, or substitute an empty string with a warning.
  This is an implementer choice; see Open Questions.
- Add a recursion depth guard (e.g. 64) to catch a macro that calls itself.

### 5.4 Manual examples (use as regression tests)

**Basic** (5.6.1): the file `Node-Macro.pre` holds `N,1,$1$,$2$,$3$`. Then
```
LoadMacro,Node,.\Node-Macro.pre
Macro,Node,13.52,15,100.25          -> executes N,1,13.52,15,100.25
```

**Graphing** (5.6.2): the file `SRSS-macro.pre`:
```
READSPEC,$1$,1,1
READSPEC,$2$,1,2
READSPEC,$3$,1,3
SRSS,4,1,2,3
WRITESPEC,$4$,4
SPECPLOT,1,2,3,4
CAPTUREPLOT,$4$.png
CLOSEPLOT
```
Driver: `LOADMACRO,SRSS,SRSS-macro.pre` then one call per node, e.g.
`MACRO,SRSS,Node1x.rs,Node1y.rs,Node1z.rs,Node1srss.rs`, and so on through node 5214. Each call
reads three spectra into plot lines 1-3, SRSS-combines them into line 4, writes line 4 to the file
named by arg 4, plots, saves a PNG named `<arg4>.png`, and closes the plot.

**Nested** (5.6.3): the file `Nested-macro.pre` holds
`MACRO,SRSS,Node$1$x.rs,Node$1$y.rs,Node$1$z.rs,Node$1$srss.rs`. The driver loads both macros and
then calls `MACRO,NESTED,1` ... `MACRO,NESTED,5214`. Result is identical to the graphing example.

Demo 3 of the manual has more post-processing macros.

---

## 6. Section 5.9: Variables and loops [CORE for the interpreter]

### 6.1 Variable model

- A variable = **(name, list of strings, integer counter)**.
- Created or overwritten with `VAR,<Name>,<X1>,...,<Xn>`. The list may be empty (`Var,NNUM`).
- Names are **case-insensitive**. Redefining with a different case overwrites the first.
- The list is **immutable** after declaration and changes only when VAR redefines the variable.
  Every VAR resets the **counter to 0**. Other list producers in Chapter 9: `LOADVAR`, `RND`,
  `ADDRND`, `REDUCESET`, `CRITFREQ`.
- The counter can change only when it is referenced inside a valid command. `SETVAR,...` is a
  no-op placeholder command used purely to evaluate counter expressions (set or check counters)
  without touching the model.
- `SHOWVAR,<name>` and `VARLIST` inspect variables.

### 6.2 Reference syntax inside any command argument (except FOREACH's own variable argument)

Using `VAR,X,1.23,2.83,3` as the running example:

| Expression | Effect on counter | Replaced by |
|---|---|---|
| `@X` | none | current counter (default 0) |
| `@X+k` | counter += k | the **new** counter value |
| `@X-k` | counter -= k | the new counter value |
| `@X++` | counter += 1 (same as `@X+1`) | the new value |
| `@X--` | counter -= 1 | the new value |
| `@X=k` | counter = k (k may be negative, e.g. `@X=-15`) | k |
| `@X[i]` | none | i-th list item, **1-based** (`@X[1]` -> `1.23`, `@X[2]` -> `2.83`) |
| `#` (inside a FOREACH body) | - | the loop counter (1-based iteration index). The manual says `#` may also take the postfix operators |
| `@X[#]` | none | list item at the current loop index (see 6.3 for which loop) |

Semantics are **pre-increment**: the counter changes first, then the new value is substituted.
This follows from "@X+5 will add 5 ... and that value will replace the entire statement" and from
the loop example, where node numbers start at 1 when NNUM starts at 0.

Evaluate the arguments **left to right**. Each `@` expression is evaluated exactly once per
command execution.

### 6.3 FOREACH

`FOREACH,<Var>,<Command...>`. `<Var>` has no `@`. Everything after it is **one command line**
(which may itself be `FOREACH,...` or `MACRO,...`).

```
FOREACH(var, body_text):
    v = vars[lower(var)] or error
    if var is already the loop variable of an enclosing FOREACH: error "nested loops cannot reuse
       the same variable" and do not execute       # manual: would loop forever, must error
    push loop frame (var, index=0)
    for i in 1..len(v.list):
        frame.index = i
        execute_command_line(body_text)   # substitution is deferred: @/# are resolved now,
                                          # at each iteration, not when FOREACH is parsed
    pop frame
```

- Iteration count = number of items in the variable's list (`FOREACH,X,...` runs 3 times).
- Only one command per loop. For more, loop over a `MACRO` call.
- `#` resolution (derived from the nested example): inside `@V[#]`, `#` is the index of the
  FOREACH loop **over V**. A bare `#` elsewhere means the innermost loop index.

**Manual example (Loop.pre), a regression test:**
```
Var,NNUM
Var,X,1,2,3,4,5
Var,Y,1,2,3,4,5
Var,Z,1,2,3,4,5
ForEach,Z,ForEach,Y,ForEach,X,N,@NNUM++,@X[#],@Y[#],@Z[#]
```
Expected: 125 nodes. Z is outermost and X innermost, so node
`n = 25*(iz-1) + 5*(iy-1) + ix` sits at `(ix, iy, iz)`. Node 1 = (1,1,1), node 2 = (2,1,1),
node 6 = (1,2,1), node 125 = (5,5,5).

**Counter unit test** (derived): `VAR,X,1.23,2.83,3`. Then `SETVAR,@X+5` -> counter 5;
`SETVAR,@X++` -> 6; `SETVAR,@X=-15` -> -15; `VAR,X,9` -> counter 0.

---

## 7. Section 5.7: Merging models (SSSI and excavation) [CORE]

### 7.1 MERGE, for building SSSI (multi-structure) models

`MERGE,<Mdl1>,<Mdl2>,<X>,<Y>,<Z>`. The result goes into the **active model** (in the manual
example, model 3, set by `Actm,3` first).

- Model 1's data is copied **unchanged**.
- Model 2's data is appended with **numbering offsets** so nothing collides: node numbers, group
  numbers, material numbers, "etc." (derived: also real/spring/matrix property indices). The
  offset rule is (derived): `offset = max id used in Mdl1` for each table, and every reference
  inside Mdl2 (element node lists, MSET/RSET indices, mass/BC/interaction node ids) is remapped.
  Element numbers are local to each group and stay as they are because groups are new.
- Model 2's nodes are **translated** by the constants: `x += X`, `y += Y`, `z += Z`.
- **No coincident-node merging** happens in MERGE (use WELD or MERGESOIL for that).
- Rotation: use `ROTATE,<x>,<y>,<z>,<rxy>,<ryz>,<rzx>` on one model **before** MERGE.
  `(x,y,z)` is the rotation center. The angles are in degrees: `rxy` about Z, `ryz` about X,
  `rzx` about Y. The order in which multiple rotations apply is not stated (Open Questions).
  `TRANSLATE,<x>,<y>,<z>` (defaults 0) moves all nodes and nothing else.
- Undefined in the manual: whether model options, analysis options, frequencies and the soil
  layer table come from Mdl1 (recommended) and whether Mdl2's soil layers are appended.

Manual example:
```
Actm,1
Inp,Model1.pre
Actm,2
Inp,Model2.pre
Actm,3
Merge,1,2,0,0,0          (models already positioned; zero offset)
```

Warning (from the SOILMESH reference): SOILMESH also offsets ids. Do not chain SOILMESH and MERGE,
or the ids will have doubled gaps.

### 7.2 MERGESOIL, to combine a structure FE model with its excavation-volume model

`MERGESOIL,<Struct>,<Soil>,[Mode],[Stiff],[Stiff2],[SepLevel],[Mapping]`. The result goes into the
**active model**.

Preconditions:

- Both models use the **same coordinate system**.
- The structure model is all structural elements (`ETYPEGEN,1`).
- The excavation model is all excavation elements (`ETYPEGEN,2`).
- **Ground elevation must be set** (`GROUNDELEV,<elev>`, or through the HOUSE options) before the
  merge, so the interface can be welded correctly.
- The excavation model may be generated with `EXCAV,<model>,[delta]` (Section 7.3) if the basemat
  has an identical node-layer mesh at every embedded elevation.

Behavior:

- Every element of `<Soil>` becomes a **soil (excavation) element**, and its materials are
  **converted into soil layers** (L table).
- Interface treatment by `[Mode]`:

| Mode | Meaning |
|---|---|
| 0 | Fully **unbonded** interface on all sides: no node merge, no springs |
| **1** (default) | Fully **bonded** on all sides by **merging coincident nodes** |
| 2 | **Stiff springs** link every coincident structure/excavation node pair, all sides |
| 3 | Stiff springs **below** `SepLevel` (global Z) and **soft springs above** it (soil separation near the surface) |

| Arg | Default | Meaning |
|---|---|---|
| `Stiff` | **1.0E7** (the text extraction shows "107"; the PDF shows 10^7) | spring constant for stiff springs (Mode 2 and 3) |
| `Stiff2` | **10** | spring constant for soft springs above SepLevel (Mode 3) |
| `SepLevel` | - | global Z of the separation level |
| `Mapping` | - | output text file listing, for each excavation-model node, its **old number and new number** |

- When nodes merge, the **higher-numbered node** of a coincident pair stops being used. It stays
  in the node list and can be removed with `RMVUNUSED`.
- Option AA: MERGESOIL is **required** to join the structure and excavation models that were
  each converted from ANSYS `.cdb`. In this case the mapping file **must be named
  `modelname_Excv.map`** (Demo 7).

Derived algorithm:

```
S = copy(Struct); E = copy(Soil) with node/group/material ids offset past S
mark all E elements ETYPE=2; for each E material -> new L entry
    (Vs = sqrt(G/rho), Vp = sqrt(M/rho), rho = weight/g, G = E/(2(1+nu)),
     M = E(1-nu)/((1+nu)(1-2nu)); damping = material damping; layer thickness: Open Question)
pairs = {(s,e) : |x_s - x_e| <= tol, s in S nodes, e in E nodes}     # tol: Open Question
Mode 1: substitute max(s,e) -> min(s,e) in all connectivity; write Mapping (old,new)
Mode 2/3: create SPRING group; for each pair add SPRING element (s,e) with property
          SC(k,k,k, kr,kr,kr, damp), k = Stiff (Mode 2, or Mode 3 with z < SepLevel) else Stiff2.
          Rotational components and damping: Open Question (recommend translational only, damp 0)
Mode 0: no connection
```

Manual example:
```
Actm,1
Inp,struct.pre
etypegen,1
* Adjust Ground Elev. before merge
Actm,2
Inp,Soil.pre
* Adjust Ground Elev. before merge
etypegen,2
Actm,3
MergeSoil,1,2,1,modelname_Excv.map
```

The example puts the mapping file as the **4th** argument, which by the formal syntax is the
`[Stiff]` slot (see Open Questions). Recommended parser rule: if a numeric slot from position 4 on
holds a token that is not a number and it is the **last** argument, treat it as `[Mapping]`.

### 7.3 Supporting commands (signatures only; full spec in the Chapter 9 file)

| Command | Purpose |
|---|---|
| `ETYPEGEN,<type>` | Sets all elements of the active model: 0 = implicit (by location relative to ground level), 1 = all structural, 2 = all excavated soil. Buried shells with ETYPE 2 have interaction nodes that are not part of the excavation volume. Non-structural elements are ignored by the ANSYS interface (Option A). |
| `EXCAV,<model>,[delta]` | Builds an excavation-volume model in `<model>`. It uses the lowest z-level node grid as the template for a homogeneous mesh up to the ground surface. Areas that do not reach the bottom level are not generated. `delta` (default 0, positive float) groups slightly different z values into one level. |
| `WELD` | Merges coincident nodes by building a substitution table and rewriting element connectivity. Substituted nodes stay in the list unconnected. Masses, BCs and interaction flags are **not** changed. It affects **all** coincident nodes, even intended ones (e.g. FIXROT springs, structure/excavation pairs). Meant for structure-only or excavation-only models before MERGESOIL. |
| `RMVUNUSED` | Removes nodes not referenced by any element. |
| `GROUNDELEV,<elev>` | Sets the ground elevation without resetting other HOUSE options. |
| `GRAVITY,<grav>` | Sets the gravity constant. |

---

## 8. Section 5.8: Submodels, section cuts, section-cut calculations [CORE][UI]

### 8.1 Cuts

- A **cut** is a numbered, user-defined **set of elements** (each a group number + element
  number) used for submodeling and for engineering calculations. It is viewed with Plot > Cuts
  (`CUTPLOT`).
- Builders (all **add** to an existing cut, which is created if missing):
  - `CUTADD,<cut>,<group>,<e1>,...,<eN>` adds a list.
  - `CUTADD,<cut>,<group>,RANGE,<start>,[end],[stride]` adds a range (end defaults to start,
    stride to 1).
  - `CUTVOL,<cut>,[Xmin],[Xmax],[Ymin],[Ymax],[Zmin],[Zmax]` adds the elements inside a box.
    Every blank bound defaults to the active model's min/max extent.
  - `SLICE,<cut>,<px>,<py>,<pz>,<nx>,<ny>,<nz>` adds elements that cross an infinite plane. In the
    old SUBMODELER it replaced the cut; now it adds.
- Editors: `CUTRMV` (same two forms as CUTADD) removes elements. `CUTCLR,<first>,[last],[step]`
  deletes cuts; using a cleared cut is an error until it is rebuilt.
- Box membership test for CUTVOL is not specified (Open Questions). Recommended: all element
  nodes within the closed box, with a small tolerance. Example 1 below ("elements above ground,
  Zmin = 2.53") is consistent with this.

### 8.2 Submodels from cuts and transfer commands

- `CUT2SUB,<cutnum>,<dest>,[solid]` copies the cut's elements from the **active** model into model
  `<dest>`, with the nodes they reference. If `solid >= 1` (default -1), shells are turned into
  solid "thick shell" elements **for plotting only**: there is no self-intersection check.
- Transfer commands write straight into another model and **overwrite** existing elements with the
  same numbers. Group types and referenced nodes in the destination are overwritten too; the
  active model is unchanged. The manual recommends **cuts plus CUT2SUB instead**, because the
  selection can be reviewed first.
  - `TRANELEM,<dest>,<group>,<begin>,<end>,<stride>`
  - `TRANVOL,<dest>,[Xmin],[Xmax],[Ymin],[Ymax],[Zmin],[Zmax]`
- Other related commands: `EXTRACTEXCAV,<Model>` (submodel of the explicitly typed excavation
  elements) and `SPLITGROUP,<group>,<split>,[dir]`.

**Manual example 1** (submodel of everything above ground at z = 2.53):
```
inp,model.pre
cutvol,3,,,,,2.53            (blank = default extent; <blank> may be empty or whitespace)
cut2sub,3,1
actm,1
write,model above ground.pre,C:/users/user/.../
```
Note that file names may contain spaces, because arguments are comma-delimited.

### 8.3 Cross-section model and section-cut calculations

Workflow:

1. Load element stresses onto the original model: `READSTR,<file.ess>,[Dir]`. The `.ess` files
   come from STRESS when `SECDATAOPT,1` is set (frames named `ESTRESS_<frame>.ess`). This must
   happen **before** CSECT.
2. Build a cut (CUTADD/CUTVOL/SLICE).
3. `CSECT,<dest>,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>` intersects every element in the cut with
   the infinite plane through `p` with normal `n`, and stores a **cross-section model** in `<dest>`.
4. `ACTM,<dest>`, then run the calculations:
   - `CALCPAR,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sysno>,[verbose]` at one time (the loaded stress
     state).
   - `CALCSECTHIST` (from a list of stress files) or `CALCSECTHISTDB` (from the binary stress
     database) for **time histories** of section forces.
   - Also available: `CALCM` (total mass), `CALCC` (centroid), `CALCMOI` (moments of inertia).

Cross-section model rules (manual):

- It keeps the **same group and element numbers** as the original. Node numbers and locations are
  **new**.
- Each element's extent is **clipped to its intersection** with the plane.
- **Shells**: the shell/plane intersection is found and the in-section **thickness is adjusted
  when the shell is oblique** to the plane. Derived formula, with shell mid-surface unit normal
  `m`, section normal `n`, and shell thickness `t`:
  `w = t / |m x n|` = `t / sin(alpha)`, where alpha is the angle between `m` and `n`.
  A wall perpendicular to the section plane gives `w = t`.
- The 2D intersection footprint is **extruded to unit thickness along n** to make a 3D model,
  used for plotting instead of a flat 2D one.
- The calculation commands also work on a whole model that is not a cross-section (e.g. global
  base forces and moments of a building).

CALCPAR outputs: area centroid (the origin of the local system), moments of inertia, and the
**six section resultants** (axial, two shears, torsion, two bending moments) in a local system
defined by the section normal `n` and the "right" vector `r`, stored as system `sysno`. Default
`verbose = -1` prints a labeled table (`<parname> = value`). Any other value prints all
parameters on one line, unlabeled, in the same order.

Derived resultant algorithm (verify against the manual's Demo 8 / Verification Problem 47):

```
for each section piece i: area A_i, centroid c_i, global stress tensor S_i (from .ess, element center)
C   = sum(A_i c_i) / sum(A_i)
t_i = S_i . n                      # traction on the section plane
F   = sum(t_i A_i);  M = sum((c_i - C) x t_i A_i)
local axes: e1 = n (normal), e2 = r projected onto the plane and normalized, e3 = e1 x e2
report P = F.e1, V2 = F.e2, V3 = F.e3, T = M.e1, M2 = M.e2, M3 = M.e3 ; I22, I33 about C
```

For SHELL elements, `.ess` may hold stress resultants rather than stresses. The implementer must
map them (Open Questions).

**Manual example 2** (wall section):
```
inp,model.pre
readstr,stressfile.txt
cutvol,3,52.5,52.8,-320,320,2.53,45.22
csect,1,3,55.6,0,4,1,0,0
actm,1
calcpar,1,0,0,0,1,0
calcm
```
Caution: this example looks internally inconsistent (Open Questions), so do not use it as a
numeric regression test.

---

## 9. Chapter 6 overview: complete menu tree [UI]

Main window (screenshots, PDF p.89):

- Title `ACS-SASSI User Interface` (older builds: `SSI Submodeler`).
- Menu bar `Model | File | Plot | Modules | Options | View | Help`.
- Two toolbars: **Main** and **3D Plot**; their icon maps are in the Chapter 8 spec.
- Tabbed document area: a `Command History` tab plus plot and editor tabs; module runs open their
  own output tabs.
- `Command Entry` single-line box at the bottom; a status bar with the progress bar at the lower
  right.

| Menu | Item | Sub-item | Function | Command equivalent |
|---|---|---|---|---|
| Model | New | | empty model, unused number, made active | (`ACTM` on an unused number) |
| | Open (Ctrl+O) | | open from database | (database + RESUME) |
| | Save | | save model (binary database format) | `SAVE` |
| | Input | | load a `.pre` file | `INP` |
| | Converters | SASSI .hou | model from HOUSE `.hou` (+ `.sit`, `.poi`) | `CONVERT,SSI,...` |
| | | ANSYS .cdb | model from ANSYS `.cdb` | `CONVERT,ANSYS,...` |
| | | GT-STRUDL Database | (not included) | `CONVERT,STRUDL,...` (N/A) |
| | Output | | write model to `.pre` | `WRITE` |
| | Export to ANSYS | | write ANSYS APDL input | `ANSYS` |
| | Export to STRUDL | | (not included, greyed) | - |
| | Exit | | close UI, save settings to SASSIini.xml | - |
| File | Open | | open/create ASCII file in editor | - |
| | Export Image | | save plot image (.bmp) | ~`CAPTUREPLOT` |
| | Export Table | | export 2D-plot data to CSV | - |
| Plot | Model | Elements / Nodes | element plot / node plot | `MODELPLOT` / `NODEPLOT` |
| | Cuts | | plot a section cut of the active model | `CUTPLOT` |
| | Spectrum TFU-TFI | | load spectrum-format files, 2D line plot | `SPECPLOT` (+`READSPEC`) |
| | Time History | | load time-history files, 2D line plot | `THPLOT` (+`READTH`) |
| | Soil Layers | | soil layer plot | `LAYERPLOT` |
| | Soil Properties | | soil property plot | `SOILPROPPLOT` |
| | Non Uniform Soil Field | | (not active, greyed) | - |
| | Process Animation Frame List | | process frames for animation | `PROCFRAME` |
| | Bubble / Vector / Contour / Deformed Shape | | animations | `BUBBLEPLOT` / `VECTORPLOT` / `CONTOURPLOT` / `DEFORMPLOT` |
| Modules | Location | | module executable locations | - |
| | Extension | | module input/output extensions | - |
| | EQUAKE, SOIL, SITE, POINT, HOUSE, FORCE, ANALYS, COMBIN, MOTION, STRESS, RELDISP | | run module on active model | `RUNEQUAKE`, `RUNSOIL`, `RUNSITE`, `RUNPOINT`, `RUNHOUSE`, `RUNFORCE`, `RUNANALYS`, `RUNCOMBIN`, `RUNMOTION`, `RUNSTRESS`, `RUNRELDISP` |
| | LIQUEF, PINT | | not included (greyed) | - |
| | NONLINEAR | | run NONLINEAR (the table says "Run the PANEL module") | (none; see `NONLINBAT`) |
| | ANSYS Eq. Static Load | | LOADGEN, static | none (dialog only) |
| | ANSYS Dynamic Load | | LOADGEN, dynamic | none (dialog only) |
| | ANSYS Super Element Utilities | | SSI2ANSYS (MATRIX50 SE) | none (dialog only) |
| Options | Model / Write / Check / Analysis | | model options, extended write options, check options, analysis options | `MOPT`, ... (see Options spec, manual 6.5) |
| | Windows Settings / Colors / Font / Shader Options / Reset Plot | | display settings | `WINDOWSETTINGS`, `COLOR`, -, `SHADEROPTIONS`, (`RSTVIEW`/`RSTCENTER`, approximate) |
| View | Check Errors | | error-checking window | - |
| | Command Window | | show/open command history | - |
| | Command Display | Command Echo / Output Confirmation / Comments / Warning & Errors | filters | - |
| | Toolbars | Main Toolbar / Plot Toolbar | show/hide | - |
| Help | Help | | online help in browser | - |
| | About | | about box | - |

Menu separators (screenshots):

- Model: New, Open, Save | Input, Converters | Output, Export to Ansys, Export to Strudl | Exit.
- Plot: Model ... Non Uniform Soil Field | Process Animation Frame List | Bubble ... Deformed
  Shape.
- Modules: Location, Extension | EQUAKE ... NONLINEAR | three ANSYS items.
- Options: Model, Write, Check, Analysis | Windows Settings, Colors, Font | Shader Options, Reset
  Plot.

For a Python reimplementation, every menu action should call the **same command handler** as its
command equivalent, so that GUI and script results match exactly. This is a verifiability
requirement.

---

## 10. Section 6.1.5: SASSI fixed-format (.hou) converter [IO]

Dialog **"SASSI .hou to .pre Converter"** (PDF p.91):

| Field / button | Type | Default / example | Required | Meaning |
|---|---|---|---|---|
| Input File Name | text + `<<` browse | e.g. `C:/test.hou` | **yes** | HOUSE input file to convert |
| Output .pre File Name | text + `<<` browse | empty | no | if given, the converted model is also written as `.pre` |
| Save Converted Data to Model Number | integer | e.g. `1` | no | non-negative integer means store into that model; blank means the **active** model |
| Convert | button | | | run conversion |
| Cancel | button | | | close |
| (static text) | | | | disclaimer: limited testing, possibly inaccurate data, check every model before simulating |

Behavior:

- Parses the HOUSE fixed-format input. This is the inverse of the AFWRITE HOUSE writer; the
  record layout is in the HOUSE input spec.
- If `<base>.sit` (SITE) and `<base>.poi` (POINT) exist **in the same folder with the same base
  name**, they are read too and their data (soil layers, site and point options) is added to the
  model. The manual claims "no limitation" for this converter.
- Legacy SASSI2000 inputs need manual preparation:
  1. rename the HOUSE file to end in `.hou`;
  2. append a few zeros to the **3rd general-information (option) line** of the HOUSE file;
  3. delete any `$` comment markers;
  4. in the SITE file, **insert a new first line holding a zero in column 5**.

  An implementation may automate these steps in a "legacy SASSI2000" mode (UI choice).
- Command form: `CONVERT,SSI,<model>,<filename>`.

---

## 11. Section 6.1.6: ANSYS .cdb converter (ANSYS to ACS SASSI) [IO][CORE]

### 11.1 Dialog and command

Dialog **"ANSYS .cdb to .pre Converter"** (PDF p.92):

| Field / button | Type | Default / example | Required | Meaning |
|---|---|---|---|---|
| Input File Name | text + `<<` | e.g. `C:/test.cdb` | **yes** | ANSYS CDWRITE file |
| Output .pre File Name | text + `<<` | empty | no | if blank, the model lives only in memory; use WRITE later |
| Save Converted Data to Model Number | integer | blank | no | blank means the **active model** |
| Enter Value for Gravity | float | **32.2** | **yes** | gravity in the model's length/time units |
| Convert / Cancel | buttons | | | |
| (static text) | | | | same disclaimer as the .hou converter |

- Opened from the menu or by `CONVERT` with no arguments. Scripted form:
  `CONVERT,ANSYS,<model>,<filename>,<gravity>`.
- Gravity is needed to turn ANSYS mass-density data into ACS SASSI **specific weight** (the M and L
  commands take weight per unit volume) and to set the model gravity constant.
  - Derived formula: `weight = DENS * g`.
  - Units warning: 32.2 is ft/s^2. For lbf-in-s models use 386.4; for SI (N, m, kg) use 9.81.
- **WARNING (manual)**: only a subset of ANSYS elements is converted, and not every KEYOPT
  variant. **Any option not listed here is unsupported.** The conversion may still "succeed" but
  produce an **incorrect** ACS SASSI model. The implementation must emit a warning for every
  unsupported element type, KEYOPT, real-constant layout or command it skips.

### 11.2 Global data rules

| Item | Rule |
|---|---|
| Materials | read **only from `MPDATA`** records. Any other way of defining or redefining material data (e.g. TB tables, later MP edits not written as MPDATA) is ignored |
| Real constants | read from the real-constant block. The manual writes "RBLOCK"; CDWRITE files write `RLBLOCK`. Accept both |
| Sections | read only from combinations of `SECTYPE`, `SECBLOCK`, `SECDATA`. Other section commands are not transferred |
| Unsupported elements | ignored, i.e. not in the model, with a warning (except the Option AA display list, 11.5) |

### 11.3 Supported ANSYS element types (structural conversion)

| ANSYS element | SASSI target (GROUP type) | Conditions (manual) | Implementer mapping (derived; verify against ANSYS element docs) |
|---|---|---|---|
| **BEAM4**, **BEAM44** (legacy from ANSYS V15) | BEAMS (2) | Must use the **K node** (3-node definition). Properties from the real-constant block, which must hold **6, 8, 10, 12 or 19-24 fields**. End releases come from **KEYOPT(7)** (node I) and **KEYOPT(8)** (node J). Section commands (SECTYPE/SECBLOCK/SECDATA) are the only other recognized property source. | `E,ne,I,J,K`. BEAM4 RLBLOCK order: AREA, IZZ, IYY, TKZ, TKY, THETA, ISTRN, IXX, SHEARZ, SHEARY, SPIN, ADDMAS (6/8/10/12 fields). BEAM44 adds node-J properties (fields 7-12), offsets DX1..DZ2 (13-18), SHEARZ, SHEARY (19-20), TKZT1..TKYT2 (21-24). Axis mapping: see 11.4. |
| **COMBIN14** | SPRING (7) | **KEYOPT(2) or KEYOPT(3)** must be defined for the type. **KEYOPT(2) wins** if both are set. Spring constant from the real-constant block. Supported: **KEYOPT(2) = 1..6**; **KEYOPT(3) = 1, 2**. | KEYOPT(2)=1..6 means a 1-D spring on UX, UY, UZ, ROTX, ROTY, ROTZ, so `SC,n,k,0,0,0,0,0,damp` with k (= R1) placed in component 1..6. KEYOPT(3)=1 (torsional) and =2 (2-D longitudinal) are axial along I-J. This is exact only if I-J is aligned with a global axis. Otherwise use a 2-node GENERAL element with K = k*(e e^T) in global axes (Open Questions). COMBIN14 viscous damping CV1/CV2 has no ACS SASSI counterpart (SC damp is a ratio), so warn. |
| **MASS21** | nodal masses (`MT`, `MR`) | Masses must come from the real-constant block. **KEYOPT(1) and KEYOPT(2) = 0 or undefined**. Mass units then follow `MUNITS`. | KEYOPT(1)=0 means the constants are masses/inertias; KEYOPT(2)=0 means global axes. Fields: MASSX, MASSY, MASSZ, IXX, IYY, IZZ (KEYOPT(3)=0); a single MASS (KEYOPT(3)=2,3). Emit `MT,n,mx,my,mz`, `MR,n,ixx,iyy,izz`, and `MUNITS,n,,,0` (mass units). Sum several MASS21 at one node. |
| **SOLID45** (legacy from V15) | SOLID (1) | "Fully convertible"; **some highly distorted pyramids** may not convert. | ANSYS I,J,K,L,M,N,O,P map one-to-one to SASSI nodes 1-8 (SASSI Fig. 9.1: 1-4 bottom counter-clockwise, 5-8 above them). Degenerate forms keep their repeated nodes. A pyramid (M=N=O=P) matches the SASSI form 5=6=7=8. A tetrahedron (K=L, M=N=O=P) is allowed. Check positive Jacobian and reorder if needed. |
| **SHELL63** (legacy from V15) | SHELL (3) (same FE formulation, Kirchhoff) | Thickness **must** come from the real-constant block. | Thickness = TK(I) (R1). If TK(J..L) are non-zero and differ, use the average with a warning. Triangle (K=L): write 3 nodes and omit L, as SASSI requires for triangles. Nodes must be counter-clockwise; ANSYS ordering is compatible. |
| **SHELL181** | SHELL (3) (thin, similar to SHELL63) | Thickness **must** come from section commands. | Thickness = sum of layer thicknesses of the SECTYPE,SHELL section (SECBLOCK/SECDATA). Material from the element MAT (or the single layer MAT). Multi-layer composites are unsupported, so warn. Note that ANSYS's thick-shell (Mindlin) behavior is **not** carried over; mapping to TSHELL is an optional implementer extension. |
| **SOLID185** | SOLID (1) (similar to SOLID45) | KEYOPT(4) ("nonuniform materials" option, as the manual calls it) undefined or 0. | Same node mapping as SOLID45. |
| **BEAM188**, **PIPE288** | BEAMS (2) | Pipes become **equivalent straight beams**. **K node required**. **End releases are not possible**: ANSYS ENDRELEASE adds nodes and couplings that ACS SASSI cannot represent. Properties **must** come from section commands, of type **ASEC or RECT**. **No section offset** for pipes. | ASEC SECDATA: A, Iyy, Iyz, Izz, Iw, J, CGy, CGz, SHy, SHz, TKz, TKy, giving axial=A, tors=J, I's per 11.4; SHy/SHz are shear-center offsets, **not** shear areas. RECT SECDATA B, H: A = BH, Iyy = BH^3/12, Izz = HB^3/12, J by the rectangle torsion formula, shear areas 5A/6 (optional). PIPE section (Do, tw): A = pi(Do^2-Di^2)/4, I = pi(Do^4-Di^4)/64, J = 2I. |

### 11.4 Beam local-axis and release mapping (derived; verify)

- SASSI BEAMS (manual Fig. 9.4): axis 1 runs from I to J. **Axis 2 points toward the K node**
  (K lies in the 1-2 plane). Axis 3 = 1 x 2.
- ANSYS BEAM4/44/188 (ANSYS documentation): element x runs from I to J. The **K node lies in the
  x-z plane**.
- With the same K node: SASSI axis 2 = ANSYS z and SASSI axis 3 = ANSYS -y. Hence:

| SASSI `R` field | ANSYS source |
|---|---|
| `<axial>` | AREA / A |
| `<shear2>` (shear along axis 2 = ANSYS z) | AREA/SHEARZ (BEAM4/44, if SHEARZ > 0), else 0 |
| `<shear3>` | AREA/SHEARY (if > 0), else 0 |
| `<tors>` | IXX (BEAM4/44; if 0, ANSYS uses IYY+IZZ) / J (ASEC) |
| `<flex2>` (inertia about axis 2) | IZZ / Izz |
| `<flex3>` (inertia about axis 3) | IYY / Iyy |

- Robust method (recommended): do not hard-code the table. For each element, compute both local
  frames from (I, J, K) and map every inertia and shear term by matching axes (with sign). Run a
  unit test with an asymmetric rectangle (B is not H) and compare the HOUSE stiffness against
  ANSYS.
- Release codes:
  - BEAM44 KEYOPT(7)/(8) are 6-digit patterns, read as digits for UX UY UZ ROTX ROTY ROTZ
    (1 = released).
  - SASSI `KI`/`KJ` take k1..k6 for P1, P2, P3, M1, M2, M3 in the local axes.
  - Mapping: k1 = UX, k2 = UZ, k3 = UY, k4 = ROTX, k5 = ROTZ, k6 = ROTY.
- Tapered BEAM44 (node-J properties differ): SASSI beams are prismatic. Use the node-I values (or
  the average) and warn. Non-zero offsets DX..DZ cannot be represented, so warn.
- ADDMAS (BEAM4 added mass per length) can be lumped as `ADDMAS*L/2` at each end node via MT
  (optional extension, warn).

### 11.5 Option AA (topology-only) conversion [ADV]

- In Option AA the structural matrices come from ANSYS (COOSK*/COOSM* files). HOUSE uses the
  `.hou` file **only for node coordinates and element connectivity**. The converted model is
  therefore only a **topology/display model**.
- The model is **flagged** as possibly containing features incompatible with ACS SASSI FE
  modeling. The `.hou` that AFWRITE writes for it is **not runnable by HOUSE(FS)**, and a
  **warning** appears at conversion time. Only HOUSEFSA/ANALYSFA (the Option AA modules) use it.
- Extra element types converted **for display only** (no limits; only connectivity is needed):

| ANSYS | Displayed as SASSI |
|---|---|
| TRUSS180 | Beam |
| MPC184 | Spring |
| PIPE16 | Beam |
| PIPE18 | Beam (curved elbow shown as a straight beam (derived); K is ANSYS's curvature node) |
| FLUID80 | Solid |

- All other unlisted ANSYS elements are **ignored**.
- `ANSYSMODELTYPE,<type>` sets the Advanced ANSYS operation mode: 1 = embedded, 2 = surface.
- The structure and excavation `.cdb` models are converted **separately** and then joined with
  MERGESOIL. The mapping file must be named `modelname_Excv.map` (Section 7.2).

### 11.6 `.cdb` parsing notes for the implementer (derived from the ANSYS CDWRITE format; verify on real files)

```
ET,<itype>,<ename#>                 e.g. ET,1,185   (4=BEAM4 44=BEAM44 14=COMBIN14 21=MASS21 45=SOLID45
                                                      63=SHELL63 181 185 188 288 180 184 16 18 80 50=MATRIX50)
KEYOPT,<itype>,<knum>,<value>
NBLOCK,<nfields>,SOLID,<max>,<count>
(3i9,6e21.13e3)                     <- Fortran format line: honor the field widths given
 node  solid_ent  line_loc  X  Y  Z  [THXY THYZ THZX]    (nodal rotation angles: ignore + warn if nonzero)
N,R5.3,LOC,-1,                      <- end of NBLOCK
EBLOCK,19,SOLID,<max>,<count>
(19i9)                              <- (19i8)/(19i10) also occur
 mat type real secnum esys birth solidref shape nnodes excl  elem#  n1..n8   (continuation line if nnodes > 8)
-1                                  <- end of EBLOCK
MPTEMP,R5.0,...                     <- temperature table (use the first temperature only)
MPDATA,R5.0,<len>,<Lab>,<mat>,<stloc>,<v1>,...    Lab: EX, EY, EZ, PRXY/NUXY, GXY, DENS, DMPR, DAMP ...
RLBLOCK,<nsets>,<maxset>,<maxrl>,<nperline>
(2i8,6g16.9)                        <- first line of each set: set#, nvalues, up to 6 values
(7g16.9)                            <- continuation lines, 7 values each
SECTYPE,<id>,<BEAM|SHELL|PIPE>,<subtype RECT|ASEC|...>,<name>,<refine>
SECOFFSET,...                       <- for PIPE288 must be absent or centroid (no offset)
SECDATA,<v1>,<v2>,...               <- beam section data
SECBLOCK,<nlayers>                  <- shell layer lines follow: thickness, mat, theta, nip
D,<node>,<UX|UY|UZ|ROTX|ROTY|ROTZ|ALL>,<value>   <- optional: map value 0 to SASSI D (not in manual)
CE / CP / ENDRELEASE                <- unsupported: warn
```

- SASSI requires one group per element type and **element numbers starting at 1 with no gaps**.
  Recommended grouping: one SASSI group per ANSYS element-type number (`itype`), elements
  renumbered 1..n in ANSYS order, with an (ANSYS elem# -> group, local elem#) map file. Node
  numbers are kept.
- Converted elements get ETYPE 0 (implicit). The user then runs `ETYPEGEN,1` or `ETYPEGEN,2`
  (Option AA / MERGESOIL workflow).
- Material: `M,<mat>,EX,PRXY,DENS*g,<pdamp>,<sdamp>,1`. Damping is not addressed in the manual
  (Open Questions). Recommendation: if DMPR exists, the ratio is DMPR/2 (structural damping
  coefficient = 2*zeta); else use a user default with a warning.

---

## 12. Sections 6.1.7-6.1.10: GT-STRUDL, Output, Export to ANSYS, Exit

### 12.1 GT-STRUDL Database Converter [IO] (not included)

Converts one specific GT-STRUDL database format, and breaks on any format change. Command form
(N/A in this version): `CONVERT,STRUDL,<model>,<jointfile>,<group>,<member>,<fematrib>,<prop>,<const>,<input>,<tie>`.
Only the joint file takes a full path; the others are bare names in the same directory. Do not
implement it (keep the menu item greyed).

### 12.2 Model > Output [UI][IO]

The same as `WRITE,[<file>],[<path>]`: a save dialog asks for name and path and writes the
current model as `.pre`. Default name: `<modelname>.pre` in the model path.

### 12.3 Model > Export to ANSYS [IO]

- Writes an **ANSYS APDL** (the manual spells it "ADPL") input file named **`<modelname>.inp`**
  in the **model directory**.
- With no name/path defined, the file is still written, but with an **empty base name (`.inp`)**
  in a platform-dependent working directory. The implementation should warn.
- Command: `ANSYS,[FileName],[Dir]` (defaults from the model info).
- **WARNING (manual)**: run **`ANSYSREFORMAT`** first. Group and element structures differ: ACS
  SASSI sets beam end releases **per element**, while ANSYS sets them **per element type via
  KEYOPT**.
  - `ANSYSREFORMAT,<Org>,<Map>` reads model `<Org>`, regroups the beams so that each new group has
    one release pattern, and writes `<Map>`, a mapping of old to new beam groups.
  - Run it with an **empty active model**; the reformatted model goes into the active model.
- The exporter uses elements compatible with **ANSYS V11-15** (legacy elements). Newer ANSYS
  versions may reject them. Recommended (derived) target mapping, with an optional "modern
  elements" switch:

| SASSI | Legacy (default) | Modern option | Notes |
|---|---|---|---|
| SOLID | SOLID45 (KEYOPT(1)=0 includes extra shapes, =1 suppresses; set from MOPT incompatible-mode option) | SOLID185 | node order 1-8 = I..P |
| SHELL | SHELL63, R1 = THICK | SHELL181 + SECTYPE,SHELL | |
| TSHELL | SHELL181 | SHELL181 | |
| BEAMS | BEAM4 (no releases) or BEAM44 with KEYOPT(7),(8) from KI/KJ | BEAM188 (+ ENDRELEASE warning) | invert the 11.4 axis mapping; real constants or section from R |
| SPRING | COMBIN14 with KEYOPT(2)=1..6, one element per non-zero SC component | same | SC damp ratio has no direct counterpart |
| GENERAL | MATRIX27 (stiffness KEYOPT(3)=4, mass KEYOPT(3)=2) | same | 3-node (local axes) GM must be rotated to global first |
| PLANE | PLANE42 with KEYOPT(3)=2 (plane strain) | PLANE182 | 2D models only |
| MT/MR | MASS21 KEYOPT(3)=0 | same | divide by g if MUNITS = weight units |
| M | MP,EX / PRXY / DENS = weight/g (/ DMPR = 2*sdamp) | same | |
| D | D commands | same | |
| interaction nodes | node component (e.g. `CM,SSI_INT,NODE`) | same | convenience |
| excavation elements (ETYPE 2), soil layers | **not exported** | | the Option A ANSYS interface ignores non-structural elements |

### 12.4 Model > Exit [UI]

- Closes the UI and saves the environment settings to **`SASSIini.xml`**. This file holds the
  module locations, extensions and display options; toolbar visibility is *not* saved.
- **Models and plots are NOT saved**. Prompt the user if there are unsaved changes (recommended
  improvement).

---

## 13. Section 6.2: File submenu [UI]

| Item | Behavior |
|---|---|
| Open | Section 4.5 (text editor, create-if-missing, File > Save, Input > Connect to Command Entry) |
| Export Image | Asks for a file name and saves the **active plot** as **.bmp**. Behaves like `CAPTUREPLOT` (which picks the format from the extension, e.g. `.png` in the macro example) |
| Export Table | Exports the data of certain **2D plots** (spectrum, time history) as **CSV**. Recommended layout: first column the abscissa (frequency/period/time), one column per plotted line, header row with line names |

---

## 14. Section 6.3: Plot submenu [UI]

Every item points to Chapter 7 (Plotting) and its spec file:

| Item | Manual section |
|---|---|
| Element Plot | 7.1.2 |
| Node Plot | 7.1.3 |
| Cut Plot | 7.1.4 |
| Spectrum TFU-TFI | 7.2.3 |
| Time History Plot | 7.2.4 |
| Soil Layer Plot | 7.2.5 |
| Soil Properties Plot | 7.2.6 |
| Processing Animation | 7.3.1 |
| Bubble | 7.3.4 |
| Vector | 7.3.5 |
| Contour | 7.3.6 |
| Deformed Shape | 7.3.7 |

Non Uniform Soil Field is not active.

---

## 15. Section 6.4: Modules submenu

### 15.1 Location window: "Module Directories" (6.4.1) [UI][IO]

One row per module: a label, a full-path text box and a `<<` browse button. `Ok` and `Cancel`
sit at the right.

| Row label | State | Example default (screenshot) |
|---|---|---|
| EQUAKE Module | enabled | `...\SSI Installer V3\FastSolver\EXEB\Equakeb.exe` |
| SOIL Module | enabled | `C:\SSI Files\SSI Installer V3\FastSolve...` |
| LIQUEF Module | **disabled** | - |
| SITE Module | enabled | same folder |
| POINT Module | enabled | same folder |
| HOUSE Module | enabled | same folder |
| PINT Module | **disabled** | - |
| FORCE, ANALYS, COMBIN, MOTION, STRESS, RELDISP Module | enabled | same folder |
| LOADGEN Module | enabled | empty |
| NONLINEAR Module | enabled | empty |
| SASSIANSYS Module (the SSI2ANSYS executable) | enabled | empty |

- Stored in **`SASSIini.xml`** (OS-dependent default location) and reloaded at startup. The paths
  are also used by `AFWRITE`-related batch generators such as `AFWRBAT`.
- Python implementation (decision): rows may point either to external executables (for
  cross-checking against the commercial modules) or to the built-in Python module entry points.
  Default to the built-in ones.

### 15.2 Extension window: "File Extension Options" (6.4.2) [UI][IO]

- Fields: `Module` (drop-down list of modules), `Input File Extension` (text), `Output File
  Extension` (text), `Ok`, `Cancel`. Screenshot example: EQUAKE has input `.equ` and output
  `_equake.out`.
- Picking another module shows its two extensions. Edits for several modules are kept until `Ok`
  saves **all** of them to SASSIini.xml. `Cancel` drops **all** edits made in this dialog.
- Default table: EQUAKE is shown in the screenshot. The other input extensions come from manual
  Chapter 3. Output extensions follow the batch rule `modelname_<SSI_module_name>.out` (derived).

| Module | Input ext | Output ext |
|---|---|---|
| EQUAKE | `.equ` | `_equake.out` |
| SOIL | `.soi` | `_soil.out` |
| SITE | `.sit` | `_site.out` |
| POINT | `.poi` | `_point.out` |
| HOUSE | `.hou` | `_house.out` |
| FORCE | `.frc` | `_force.out` |
| ANALYS | `.anl` | `_analys.out` |
| COMBIN | (not stated; Open Questions) | `_combin.out` |
| MOTION | `.mot` | `_motion.out` |
| STRESS | `.str` | `_stress.out` |
| RELDISP | `.rdi` | `_reldisp.out` |
| NONLINEAR | `.eql` (written by the UI; see the Option NON spec) | `_nonlinear.out` |

### 15.3 Common run semantics for every module menu item / RUNxxx command [CORE][IO][UI]

- `RUN<MODULE>,[model]`: default `-1` means the **active model**. Choosing the menu item does the
  same thing for the active model.
- Preconditions:
  - Model name and path are defined (MDL or database).
  - **AFWRITE has already written** the module input file `<name><inext>` in the model directory.
    AFWRITE first runs CHECK and skips writing analysis files that have errors.
  - Input files from earlier modules exist (dependency table below).
- **Working directory = model directory.** All files the module generates go there.
- The UI opens a **new tab** that streams the module's console output.
- Batch protocol of the original executables (manual 3.x): each module reads three lines from
  standard input:
  ```
  modelname
  modelname.<input ext>
  modelname_<SSI_module_name>.out
  ```
  run as `<MODULE>.exe < <MODULE>.inp`. The Python modules should accept the same three values
  (stdin or CLI arguments) so batch files stay interchangeable.

### 15.4 Per-module summary (6.4.3-6.4.13) [CORE]

| Module | What it does (paraphrase) | Requires | Produces |
|---|---|---|---|
| **EQUAKE** | Generates earthquake accelerograms compatible with given ground (design) spectra. A time-varying correlation between the horizontal components can be specified. Results feed SOIL, MOTION, RELDISP and STRESS. | `.equ` | accelerograms and spectra files (see EQUAKE spec) |
| **SOIL** | Nonlinear site response with an **equivalent-linear** model of soil hysteresis. The iterated (effective) properties can then be used in the SSI analysis. Units are **British** (ft thickness, **ksf** shear modulus, **kcf** unit weight, ft/s velocity) or **SI** (m, **kN/m^2**, **kN/m^3**, m/s). | `.soi` | iterated soil properties, FILE73, FILE88 (Chapter 3) |
| **SITE** | Site-response part of the problem. The **control point and wave composition** of the control motion must be defined. Computes and stores what is needed for the free-field displacement vector. The control-motion time history is **not** needed here (it is used in MOTION). Incoherency is added later, in HOUSE. | `.sit` | **FILE1** (free-field data); **FILE2** (transmitting-boundary data) |
| **POINT** (POINT2 for 2-D, POINT3 for 3-D) | Computes what is needed to form the **frequency-dependent flexibility matrix**. SITE must run first. | `.poi`, **FILE2** | **FILE3** |
| **HOUSE** | Forms the element **mass and stiffness matrices** of the discretized model (structure only, or structure plus irregular soil zone). Can run without SITE/POINT for **surface** models. Embedded models need the `.sit` (Chapter 3). Does the **random-field decomposition** for incoherent motion. | `.hou` (and `.sit` when embedded) | **FILE4 = `modelname.N4`**; **FILE77** (incoherent) |
| **FORCE** | Builds the load vector for **external load cases** (impact, rotating machinery, unit forces for flexible-foundation impedance). Of no use for seismic SSI. When foundation flexibility does not matter, ANALYS's direct impedance option is more efficient. | `.frc` | **FILE9** |
| **ANALYS** | Solves the problem at every frequency step. Steps: (1) soil flexibility matrix at the interaction nodes; (2) soil **impedance** matrix at the interaction nodes; (3) external or seismic load vectors, including incoherency; (4) solve the SSI equations at each frequency by **LU decomposition + back-substitution**, giving transfer functions for every DOF: **ATF** (seismic, from the control motion to the final motion) or **DTF** (external force to total displacement). | `.anl`; seismic: **FILE1, FILE3, FILE4**; external: + **FILE9**; incoherent: + **FILE77** | **FILE8** (binary, complex TFs) |
| **COMBIN** | Merges the frequency sets of two solutions, e.g. when more frequencies are added after a run. Inputs are two FILE8 files renamed **FILE81** and **FILE82**. | FILE81, FILE82 | new **FILE8** |
| **MOTION** | Reads the TFs from FILE8 and interpolates efficiently in the frequency domain with a **two-SDOF transfer-function model with five parameters** (other interpolation schemes are in the MOTION spec). Computes the final response at user-selected nodes. Acceleration, velocity or displacement **response spectra** can be requested per location and DOF. | `.mot`, **FILE8** only | TFU/TFI/ACC/RS text files; optional **FILE13** (formatted) with baseline-corrected nodal motions |
| **RELDISP** | **Relative displacements**. Writes the displacements of all points in **frame format** for animation. Needs MOTION's `.TFI` (Chapter 3/4). | `.rdi`, `.TFI` | TFD/THD files, frame files |
| **STRESS** | Stress, strain and force **time histories and peaks** in structural elements. | `.str`, **FILE4 (`modelname.n4`)**, **FILE8** | **FILE15** (stress time histories, formatted); **FILE14** (beam section force/moment TFs, formatted); `.ess` frames if `SECDATAOPT,1` |

- Linear SSI run order (manual Chapter 1): `SITE -> POINT -> HOUSE -> ANALYS -> MOTION -> RELDISP
  -> STRESS`. EQUAKE and SOIL are optional pre-steps. COMBIN is optional after ANALYS.
- The implementation should check FILEn presence before launching and report the missing
  predecessor module (recommended).

### 15.5 NONLINEAR module, Option NON (6.4.14) [ADV]

- Part of the Option NON capability. It is meant for an **iterative equivalent-linearization**
  procedure and should run in **batch mode**. Each iteration runs
  `HOUSE -> ANALYS (restart) -> MOTION -> RELDISP -> COMB_XYZ_THD -> NONLINEAR` until convergence.
  The first iteration starts with SITE -> POINT (Chapter 1).
- `NONLINBAT` writes a generic batch file for this loop (spec in the Nonlinear/Panel chapter).
- The menu item runs NONLINEAR once for the active model, using the same location/extension
  rules.
- Demo 9 (concrete structure) and Demo 10 (base-isolated structure) show its use.

### 15.6 LOADGEN module, Option A (6.4.15) [ADV]

- **Interactive only. There is no command-line command**, because of the number of interdependent
  options.
- Two Modules-menu items open dialogs. On `Ok`, the UI **writes the LOADGEN input file**, then
  **launches LOADGEN**. Outputs go to the model directory.
- Full details are in the separate "ACS SASSI-ANSYS Integration Capability (Options A and AA)"
  manual, Revision 4. Demos 5 and 6 show it.

**(a) "ANSYS Eq. Static Load"** (2nd step = equivalent-static stress analysis in ANSYS). Its data
comes from AFWRITE, an SSI run, and an ANSYS model. Dialog **"ANSYS Static Load Converter"**
(PDF p.110):

| Group / field | Control | Example | Meaning |
|---|---|---|---|
| *Data to Add From ACS SASSI to the ANALYS model* (manual text: "Mass Data to add from the ACS SASSI to the ANALYS module") | radio: `Displacement`, `Acceleration` (selected in screenshot), `Disp. and Accel.`, `Disp. for Soil Module` | Acceleration | sets the LOADGEN analysis-type flag and **disables** the input boxes it does not need |
| `Use Multiple File List Inputs` | checkbox | off | make several seismic load files for several chosen critical time steps in one run |
| *SSI Model and Results Input*: `Path` | text | `C:\data\model1` | base folder; the other inputs in this group are resolved against it when they are not full paths |
| `HOUSE Module Input` | text + `<<` | `C:\data\model1\Demo1_X.hou` | the `.hou` written by AFWRITE |
| `Displacement Results` | 2 texts + `<<` each; checkbox `Rotational Disp.` | (disabled for Acceleration) | 1st box: translational displacement frame file. 2nd box (enabled by the checkbox): rotational |
| `Trans. Acceleration Results` | 2 texts + `<<` each; checkbox `Rotational Accel.` | `...model1\Demo1_X.pre` | 1st box: translational acceleration frame data. 2nd box (enabled by the checkbox): rotational |
| *Ansys Model and Data Input*: `Path` | text | `C:\data\model1\` | ANSYS input data folder |
| *Mass Data for Internal Load (Ignore for Displacement)*: `Mass Type` | radio `Lumped Mass` (default) / `Master Node Mass` | Lumped | which mass representation |
| `Generate Mass Data` | checkbox | off | if checked, the selected mass-file box names an **output** file that is created or **overwritten** |
| *For Lumped Mass*: `Lumped node` | text + `<<` | `C:\data\model1\demo1_x.masl` | lumped-mass file used or created |
| *For Master Mass*: `Master Node Mass` | text + `<<` | (disabled for Lumped) | master-mass file used or created |
| *ANSYS Output File*: `ADPL File` | text + `<<` | `C:\data\model1\demo1_x.inp` | full path of the APDL file LOADGEN writes |
| `Ok` / `Cancel` | buttons | | |

**(b) "ANSYS Dynamic Load"** (2nd step = dynamic stress analysis in ANSYS by time-domain direct
integration). Dialog **"ANSYS Dynamic Load Converter"** (PDF p.111):

| Group / field | Control | Example | Meaning |
|---|---|---|---|
| *SASSI Model and Results Input*: `Path` | text | `C:/data/model1/` | full path of the folder holding the **displacement frames** |
| `HOUSE Module Input` | text + `<<` | `C:/data/model1/demo1_x.hou` | `.hou` file |
| `Ground Acceleration File` | text + `<<` | `C:/data/model1/equake.acc` | **ground** acceleration (surface foundation) or **kinematic SSI** acceleration (embedded foundation), prepared by the user |
| *ANSYS Model and Data Input*: `Path` | text | `C:/data/model1/` | ANSYS files folder |
| *Raleigh Damping Coeff.*: `Alpha`, `Beta` | 2 floats | `.5`, `.03` | Rayleigh damping `C = alpha*M + beta*K`, so `zeta(omega) = alpha/(2 omega) + beta*omega/2`. The screenshot values are **placeholders, not recommendations** |
| *ANSYS Output File*: `ADPL File` | text + `<<` | `C:/data/model1/demo1_x.inp` | APDL file to load into ANSYS |
| `Ok` / `Cancel` | buttons | | |

Notes for the dynamic load generator:

- It uses **all relative-displacement data files (relative to the free field)** from the SSI
  results. Copy them into the folder given in the "Path" box of "SASSI Model and Results Input".
- **WARNING (manual)**: Rayleigh damping is a significant limitation of the ANSYS 2nd-step dynamic
  stress analysis.

### 15.7 SSI2ANSYS module: ANSYS MATRIX50 super-elements, Option AA (6.4.16) [ADV]

- Interactive only (no command). Background:
  - An ANSYS **MATRIX50** super-element (SE) condenses many finite elements into one element
    defined by its **K, M, C** matrices (ANSYS substructuring).
  - SSI2ANSYS turns an ANSYS model that contains MATRIX50 into either a complete ANSYS model or an
    ACS SASSI-runnable model.
- Two approaches:
  - **Option 1 "SE-to-GE"**: the MATRIX50 matrices are converted into ACS SASSI **General Matrix
    (GM) elements** (stiffness/mass), giving an ACS SASSI-runnable model.
  - **Option 2 "Add-SE"**: the MATRIX50 K, M, C (extracted into the SE `.sub` files by what the manual calls the
    "GENERAL PASS" option, i.e. the ANSYS substructure generation pass) are **added to the main ANSYS structure matrices outside ANSYS**, for Option AA.
- **WARNING (manual)**: assembling MATRIX50 SEs directly inside ANSYS does **not** give accurate
  results for every problem.
- Dialog **"Super Element Utility"** (PDF p.113):

| Field | Control | Used by | Meaning |
|---|---|---|---|
| *ANSYS MATRIX50 Super Element Operation* | radio: `Convert ANSYS SE Matrices to SASSI General Elements` (SE-to-GE, default) / `Assemble SE Matrices into ANSYS Main Structure Matrices (Option AA)` (Add-SE) | both | picks the approach. Fields not needed by the chosen one are **disabled** |
| `SE Matrix Folder` | text | both | folder with the MATRIX50 `.sub` files. The SE utility input file is also written here |
| `Main Structure Matrix Folder` | text | Add-SE only (disabled for SE-to-GE) | folder with the main-structure matrices |
| `Number of Super Elements` | integer | both | number of `.sub` files to process |
| `General Matrix ID Start` | integer | SE-to-GE | first matrix property id to use in the output |
| `Element Group ID Start` | integer | SE-to-GE | first group number to use in the output |
| `Input SE Files Names (.sub) One by One:` | text + `Add`; list + `Remove` | both | `Add` appends the text box entry to the list. `Remove` deletes the highlighted entry |
| `General Element Output Folder` | text | SE-to-GE | output folder |
| `General Element Output File (.pre)` | text | SE-to-GE | output `.pre` name **without** the `.pre` extension |
| `Ok` / `Cancel` | buttons | | `Ok` writes the SE-utility input file, then launches SSI2ANSYS |

- Outputs:
  - SE-to-GE: a `.pre` file of GENERAL elements (GROUP type 9, `MXR`/`MXI`/`MXM` terms).
  - Add-SE: modified **`COOSKI_r`, `COOSK_r`, `COOSMI_r`, `COOSM_r`, `COOCI_r`, `COOSC_r`** and
    **`Node2Equ_Stru.map`**, ready for HOUSE(FSA) / ANALYS(FSA).
- The manual text names an input called "SE Utility File Name", but the dialog has no such field
  (Open Questions).
- Demo 11 is the worked example.
- SE-to-GE algorithm (derived):
  - SASSI GM elements have only 2 nodes (12x12, global axes) or 3 nodes (local axes).
  - An SE with master nodes {1..N} becomes N(N-1)/2 two-node GM elements, one per node pair
    (i, j).
  - Each element gets the coupling blocks K_ij and K_ji, plus a **share** of the diagonal blocks.
    For example, put K_ii in full on the first element that holds node i and zero elsewhere, so
    that assembly reproduces K exactly. Do the same for M.
  - SASSI GM stiffness is complex (MXR real, MXI imaginary). With hysteretic damping ratio zeta,
    `MXI = 2*zeta*MXR`. A viscous C cannot be represented frequency-independently, so warn.
  - Masses follow the `MOPT` mass/weight option.
  - Reading the binary ANSYS `.sub` file needs the ANSYS binary format. The Python version may
    instead accept HBMAT or Matrix-Market exports of K, M, C (implementer decision).

---

## 16. Persistent configuration and database files (summary) [IO]

| File | Location | Contents | Written when |
|---|---|---|---|
| `SASSIini.xml` | OS default config dir (e.g. platformdirs `user_config_dir("ACS_SASSI")`) | module executable paths, module input/output extensions, environment/display settings | Location `Ok`, Extension `Ok`, Exit |
| `SASSIdb.xml` | OS default data dir | group to model tree (model name + location only) | Load Model dialog changes |
| `<model>.<binary>` (name TBD) | model directory | full binary model (SAVE/RESUME) | SAVE |
| `<model>.pre` | model directory (default) | text commands | WRITE / Model > Output |
| `<model>.inp` | model directory | ANSYS APDL | ANSYS / Export to ANSYS |
| module inputs `<model>.equ/.soi/.sit/.poi/.hou/.frc/.anl/.mot/.str/.rdi` | model directory | module input decks | AFWRITE |
| `FILE1, FILE2, FILE3, FILE4 (<model>.N4), FILE8, FILE9, FILE13, FILE14, FILE15, FILE77, FILE81, FILE82` | model directory | inter-module binary/formatted files | modules |
| `modelname_Excv.map` | model directory | excavation old to new node map (Option AA) | MERGESOIL |

Check Options (Options > Check) are **not** persisted and reset at every start (manual 6.5.3; see
the Options spec).

---

## 17. Suggested verification tests for this section

1. **Macros**: the basic, graphing (check generated file names) and nested examples. The nested
   driver must produce the same files as the graphing driver.
2. **Loops**: Loop.pre gives 125 nodes with the coordinate law from Section 6.3; the counter
   operator test.
3. **Interpreter robustness**: an unknown command prints `<X> Command not found` and later lines
   still run; `<blank>` arguments take defaults; `*` comments are echoed during INP.
4. **MERGE**: two one-element models with known ids. Check the id offsets, the translation of
   Mdl2 nodes, and that Mdl1 data is unchanged.
5. **MERGESOIL**: structure and excavation sharing an interface of n coincident nodes.
   - Mode 1 gives n merged nodes (higher id unused) and a mapping file with n pairs.
   - Mode 2 gives n springs with k = 1e7.
   - Mode 3 splits the springs by SepLevel (Stiff below, 10 above).
6. **CSECT/CALCPAR**:
   - A solid block cut by a plane: area = expected.
   - An oblique shell: width = t/sin(alpha).
   - Uniform normal stress s: axial force = s*A and moments = 0.
   - Linearly varying stress: M = s_max*I/c.
7. **ANSYS converter** (needs synthetic `.cdb` fixtures, one per supported element):
   - nodes, connectivity, groups (renumbered from 1), materials (weight = DENS*g), thickness,
     spring component, masses with MUNITS 0;
   - a warning for each unsupported element or option;
   - a beam axis-mapping test with B not equal to H;
   - round trip SASSI -> Export to ANSYS -> (CDWRITE in ANSYS, if available) -> convert ->
     compare.
8. **Module runner**: refuse without MDL; refuse when the input file is missing (AFWRITE not
   run); outputs land in the model directory; three-line stdin protocol.

---

## 18. Open questions / ambiguities (implementer must decide)

1. **Model > New numbering**: is the "unused" model number the lowest free one or max+1? (Spec
   uses lowest.)
2. **Model > Open target**: which model number receives the opened model (active, new, or 0)?
   What does "Add Model" ask for (name, path, title?) and does it create the directory?
3. **Binary SAVE format**: not specified (proprietary). Choose a versioned format.
4. **Behavior with no MDL**: the manual says SAVE/RESUME/AFWRITE "may run" to an unknown place.
   Spec recommends a hard error. Export to ANSYS documents writing `.inp` (blank name) in a
   platform working directory; keep that, but warn?
5. **Macro**: missing `$k$` arguments, recursion depth, and path resolution for relative macro
   files (current directory set by `CD`/`MDL`? active model path? the calling file's folder?).
   The basic example uses `.\Node-Macro.pre`. MKDIR's text says relative paths are based on the
   working directory, which defaults to the UI install directory and is changed by `CD` or `MDL`.
6. **Variables**: whether `@expr` may be embedded inside a longer token (e.g. `Node@X[#].rs`);
   exact grammar of `#` with postfix operators (`#+1`?) and whether it mutates anything; which
   loop `#` refers to outside `@V[#]`; behavior of `@X[i]` out of range; whether the counter
   operators accept a variable on the right (`@X+@Y`).
7. **FOREACH**: whether the body is re-parsed each iteration (spec: yes, deferred substitution),
   and whether a body that is a FOREACH over a different variable sees the outer `#`.
8. **MERGE**: exact offset rule (max id vs. count); handling of soil layers, frequencies, model
   and analysis options, interaction nodes and forces from Mdl2.
9. **MERGESOIL**:
   - The manual example puts the mapping file in the 4th position, which conflicts with the
     formal syntax (`[Stiff]` is 4th).
   - Tolerance for "coincident" nodes.
   - Whether only nodes at or below ground elevation are interface nodes.
   - Spring rotational components and damping.
   - Mode 0 semantics (no coupling at all?).
   - How materials become soil layers: layer thickness, and Vp/Vs vs E/nu storage.
   - Which node is kept when merging (manual: the higher number becomes unused).
10. **ROTATE**: order in which the three rotations apply and sign convention (assume right-hand
    rule, applied in argument order rxy, ryz, rzx).
11. **CUTVOL**: membership test (all nodes inside, centroid inside, or any node inside).
    **Cuts**: global to the session or owned by a model?
12. **CUT2SUB**: what is copied besides elements and nodes (materials, properties, masses, BCs,
    interaction flags, options)? Are element numbers preserved or compressed?
13. **CSECT/CALCPAR**:
    - The definition of the "right" vector when it is not in the section plane.
    - Handling of beams and springs cut by the plane.
    - The `.ess` content for shells (stresses vs. resultants).
    - The order of the six quantities in the non-verbose output.
    - The example `calcpar,1,0,0,0,1,0` omits the required `<sysno>` (default?).
    - Example 2 cuts a wall with `cutvol` X in [52.5, 52.8] but then sets the CSECT plane at
      x = 55.6 with normal (1,0,0), which looks inconsistent (perhaps the normal was meant to be
      (0,0,1) at z = 4).
14. **ANSYS converter**:
    - which MPDATA labels are read (EX/PRXY/NUXY/DENS/DMPR/DAMP);
    - how damping is set;
    - what "RBLOCK" means exactly (assume RLBLOCK);
    - mapping of COMBIN14 KEYOPT(3)=1/2 and of oblique springs;
    - BEAM44 taper and offsets;
    - SHELL63 variable thickness;
    - PIPE288 section type (manual says ASEC or RECT only, yet pipes are "converted to equivalent
      straight beams");
    - SOLID185 "KEYOPT4" (the manual's naming may not match ANSYS's KEYOPT numbering; check with
      the ANSYS docs);
    - whether D constraints and nodal coordinate rotations are converted;
    - group/element numbering policy;
    - whether gravity is also written as the model GRAVITY constant;
    - whether masses are converted to weight units.
15. **ANSYS beam axis convention**: confirm against the ANSYS docs for the user's version that K
    defines the element x-z plane for BEAM4/BEAM44/BEAM188. The SASSI convention (K in the 1-2
    plane) is from manual Fig. 9.4. Also confirm the BEAM44 KEYOPT(7)/(8) digit order.
16. **Option AA flag**: how a converted model is "flagged" (a model attribute? ANSYSMODELTYPE?),
    and which AFWRITE outputs change.
17. **Export to ANSYS**: the manual gives no element list. The table in 12.3 is a recommendation.
    Also: whether implicit-type (ETYPE 0) elements above ground count as "structural" for
    export.
18. **Module extensions**: COMBIN's input extension, and every output extension except EQUAKE's
    (spec derives `_<module>.out`). Is NONLINEAR's input `.eql`?
19. **Location window**: "Currently, there are no modules for use with the UI" (manual sentence,
    unclear). The SSI2ANSYS executable is labeled "SASSIANSYS Module".
20. **LOADGEN**: the exact LOADGEN input-file format is in the separate Integration manual (not
    available). The Static dialog's 2nd "Trans. Acceleration Results" example shows a `.pre`
    file, which looks like a placeholder. The dynamic-load text says to copy the files to the
    "SASSI Model and Results Input" path while also calling it the ANSYS folder (contradictory).
21. **SSI2ANSYS**: the "SE Utility File Name" field named in the text is missing from the dialog.
    The `.sub` binary format and the COO* file formats are not documented here.
22. **Export Image**: the menu writes `.bmp`. Does CAPTUREPLOT choose the format from the
    extension? (The macro example uses `.png`.)
23. **Export Table**: the CSV layout is unspecified.
